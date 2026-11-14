"""Neural posterior estimation (simulation-based Bayesian inference) of the Svartsengi
source from NISAR LOS rates (159A, 082D) with and without GNSS.

Simulator: Mogi source (east, north, depth, volume rate) projected on each cell's line of
sight, plus per-track offset and plane (nuisance, marginalised), plus spatially correlated
noise whose covariance is estimated from the residuals of a least-squares fit to the real
data, plus GNSS model error. A network maps the observation vector to a mixture of
Gaussians over the 4 parameters; after training it gives the posterior for any data in
milliseconds. Calibration is checked on held-out simulations (do 68 % / 95 % intervals
contain the truth 68 % / 95 % of the time?).

Usage: python scripts/reykjanes/ml_npe_source.py [n_sims=80000]
"""
import sys
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import torch
from scipy.optimize import least_squares

import reykjanes_mogi as M
from reykjanes_gnss import REF, STATIONS
from reykjanes_joint import SIGMA_GNSS
from reykjanes_sources import Data

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "output" / "reykjanes"
ML = OUT / "ml"
BIN_KM = 1.2
PRIOR = dict(e=(-6, 6), n=(-6, 6), depth=(1.0, 10.0), logv=(np.log10(2e-3), np.log10(0.08)))  # dV km^3/yr
K_MDN = 8
torch.manual_seed(0)
RNG = np.random.default_rng(0)


def mogi_enu(e, n, p):
    e0, n0, d, dv = p[..., 0:1], p[..., 1:2], p[..., 2:3], p[..., 3:4]
    de, dn = e[None] - e0, n[None] - n0
    R3 = (de ** 2 + dn ** 2 + d ** 2) ** 1.5
    k = 3 / (4 * np.pi) * 1e5 * dv
    return k * de / R3, k * dn / R3, k * d / R3


class Obs:
    """Binned real observations and the fixed geometry used by the simulator."""
    def __init__(self, D):
        self.tr = {}
        for t in M.TRACKS:
            d = D.ins[t]
            key = np.c_[np.round(d["e"] / BIN_KM), np.round(d["n"] / BIN_KM)].astype(int)
            uk, inv = np.unique(key, axis=0, return_inverse=True)
            inv = inv.ravel()
            cnt = np.bincount(inv)
            keep = cnt >= 6
            avg = lambda v: np.bincount(inv, weights=v)[keep] / cnt[keep]
            self.tr[t] = dict(e=avg(d["e"]), n=avg(d["n"]), y=avg(d["y"]),
                              vec=np.stack([avg(d["vec"][:, c]) for c in range(3)], -1))
        self.gn = np.array([D.gn[s]["en"] for s in STATIONS])
        self.gref = np.array(D.gn[REF]["en"])
        self.gv = np.concatenate([D.gn[s]["v"] for s in STATIONS])
        self.n_ins = sum(self.tr[t]["y"].size for t in M.TRACKS)
        print("observation:", {t: self.tr[t]["y"].size for t in M.TRACKS}, "InSAR bins +", self.gv.size, "GNSS values")

    def forward(self, p, nuis=None):
        """p: (N, 4) [e, n, depth, dV]; nuis: (N, 6) offset, plane E, plane N per track. Returns (N, n_ins + 12)."""
        out = []
        for k, t in enumerate(M.TRACKS):
            d = self.tr[t]
            ue, un, uz = mogi_enu(d["e"], d["n"], p)
            los = ue * d["vec"][:, 0] + un * d["vec"][:, 1] + uz * d["vec"][:, 2]
            if nuis is not None:
                los = los + nuis[:, 3 * k:3 * k + 1] + nuis[:, 3 * k + 1:3 * k + 2] * d["e"] / 10 + nuis[:, 3 * k + 2:3 * k + 3] * d["n"] / 10
            out.append(los)
        ge, gnn, gz = mogi_enu(self.gn[:, 0], self.gn[:, 1], p)
        re, rn, rz = mogi_enu(self.gref[:1], self.gref[1:], p)
        g = np.stack([ge - re, gnn - rn, gz - rz], -1).reshape(p.shape[0], -1)
        out.append(g)
        return np.concatenate(out, 1)


def lsq_fit(O):
    y = np.concatenate([O.tr[t]["y"] for t in M.TRACKS])
    def res(q):
        p = np.array([[q[0], q[1], q[2], 10 ** q[3]]])
        pred = O.forward(p, q[4:][None])[0][:O.n_ins]
        return pred - y
    q0 = np.r_[-1.5, -0.3, 4.0, np.log10(0.015), np.zeros(6)]
    r = least_squares(res, q0)
    return r.x, res(r.x)


def exp_cov_fit(O, resid):
    """Empirical covariance of the residuals vs distance, fitted with sigma^2 exp(-d / L), per track."""
    out, i0 = {}, 0
    for t in M.TRACKS:
        d = O.tr[t]
        n = d["y"].size
        r = resid[i0:i0 + n]
        i0 += n
        D = np.hypot(d["e"][:, None] - d["e"][None], d["n"][:, None] - d["n"][None])
        bins = np.arange(0, 20, 1.2)
        prod = np.outer(r, r)
        cv = [prod[(D >= a) & (D < a + 1.2)].mean() for a in bins]
        best = None
        for L in np.linspace(1, 20, 39):
            s2 = max(cv[0], 1e-3)
            err = np.nansum((np.array(cv) - s2 * np.exp(-(bins + 0.6) / L)) ** 2)
            if best is None or err < best[0]:
                best = (err, s2, L)
        _, s2, L = best
        C = s2 * np.exp(-D / L) + 0.05 * s2 * np.eye(n)
        out[t] = dict(L=L, sigma=np.sqrt(s2), chol=np.linalg.cholesky(C))
        print(f"noise model {t}: sigma {np.sqrt(s2):.2f} cm/yr, correlation length {L:.1f} km")
    return out


def sample_prior(N):
    u = RNG.random((N, 4))
    lo = np.array([PRIOR[k][0] for k in ("e", "n", "depth", "logv")])
    hi = np.array([PRIOR[k][1] for k in ("e", "n", "depth", "logv")])
    q = lo + u * (hi - lo)
    p = q.copy()
    p[:, 3] = 10 ** q[:, 3]
    nuis = np.c_[RNG.uniform(-10, 10, N), RNG.normal(0, 2, N), RNG.normal(0, 2, N),
                 RNG.uniform(-10, 10, N), RNG.normal(0, 2, N), RNG.normal(0, 2, N)]
    return q, p, nuis


def simulate(O, noise, N):
    q, p, nuis = sample_prior(N)
    x = O.forward(p, nuis)
    i0 = 0
    for t in M.TRACKS:
        n = O.tr[t]["y"].size
        x[:, i0:i0 + n] += RNG.standard_normal((N, n)) @ noise[t]["chol"].T
        i0 += n
    x[:, i0:] += RNG.normal(0, SIGMA_GNSS, (N, x.shape[1] - i0))
    return q.astype(np.float32), x.astype(np.float32)


class MDN(torch.nn.Module):
    def __init__(self, nin, nout=4, K=K_MDN, h=256):
        super().__init__()
        self.K, self.d = K, nout
        self.net = torch.nn.Sequential(torch.nn.Linear(nin, h), torch.nn.GELU(), torch.nn.Linear(h, h), torch.nn.GELU(),
                                       torch.nn.Linear(h, h), torch.nn.GELU())
        self.logit = torch.nn.Linear(h, K)
        self.mu = torch.nn.Linear(h, K * nout)
        self.ldiag = torch.nn.Linear(h, K * nout)
        self.off = torch.nn.Linear(h, K * nout * (nout - 1) // 2)
        self.tril = torch.tril_indices(nout, nout, -1)

    def dist(self, x):
        z = self.net(x)
        B = x.shape[0]
        Lm = torch.zeros(B, self.K, self.d, self.d)
        Lm[:, :, range(self.d), range(self.d)] = torch.exp(self.ldiag(z).view(B, self.K, self.d).clamp(-7, 3))
        Lm[:, :, self.tril[0], self.tril[1]] = self.off(z).view(B, self.K, -1)
        comp = torch.distributions.MultivariateNormal(self.mu(z).view(B, self.K, self.d), scale_tril=Lm)
        return torch.distributions.MixtureSameFamily(torch.distributions.Categorical(logits=self.logit(z)), comp)


def train(x, q, cols, label, epochs=25):
    xs = x[:, cols]
    xm, xsd = xs.mean(0), xs.std(0) + 1e-6
    qm, qsd = q.mean(0), q.std(0)
    X = torch.tensor((xs - xm) / xsd)
    Q = torch.tensor((q - qm) / qsd)
    nval = len(X) // 10
    net = MDN(X.shape[1])
    opt = torch.optim.Adam(net.parameters(), lr=1e-3)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, epochs)
    best, best_state = 1e9, None
    t0 = time.time()
    for ep in range(epochs):
        net.train()
        perm = torch.randperm(len(X) - nval) + nval
        for i in range(0, len(perm), 512):
            b = perm[i:i + 512]
            loss = -net.dist(X[b]).log_prob(Q[b]).mean()
            opt.zero_grad(); loss.backward(); opt.step()
        sched.step()
        net.eval()
        with torch.no_grad():
            v = -net.dist(X[:nval]).log_prob(Q[:nval]).mean().item()
        if v < best:
            best, best_state = v, {k: w.clone() for k, w in net.state_dict().items()}
        if ep % 10 == 9:
            print(f"  [{label}] epoch {ep + 1}: validation -log p {v:.3f}  ({time.time() - t0:.0f} s)")
    net.load_state_dict(best_state)
    return dict(net=net, xm=xm, xsd=xsd, qm=qm, qsd=qsd, cols=cols)


def posterior(model, xobs, n=20000):
    with torch.no_grad():
        X = torch.tensor(((xobs[:, model["cols"]] - model["xm"]) / model["xsd"]).astype(np.float32))
        s = model["net"].dist(X).sample((n,)).numpy()  # (n, B, 4)
    return s * model["qsd"] + model["qm"]


def calibration(model, x, q, n_test=600):
    s = posterior(model, x[:n_test], 2000)                  # (2000, n_test, 4)
    cov = {}
    for lev in (0.68, 0.95):
        lo, hi = np.quantile(s, [(1 - lev) / 2, (1 + lev) / 2], axis=0)
        cov[lev] = ((q[:n_test] >= lo) & (q[:n_test] <= hi)).mean(0)
    return cov


def summarize(s, label):
    names = ("east km", "north km", "depth km", "dV 1e6 m3/yr")
    vals = s.copy()
    vals[..., 3] = 10 ** vals[..., 3] * 1e3
    out = []
    for i, nm in enumerate(names):
        q = np.percentile(vals[:, 0, i], [2.5, 16, 50, 84, 97.5])
        out.append(q)
        print(f"  [{label}] {nm:13s} median {q[2]:7.2f}   68 % {q[1]:7.2f} .. {q[3]:7.2f}   95 % {q[0]:7.2f} .. {q[4]:7.2f}")
    return vals[:, 0], out


def main(n_sims=80000):
    global N_SIMS
    N_SIMS = n_sims
    ML.mkdir(parents=True, exist_ok=True)
    D = Data()
    O = Obs(D)
    qhat, resid = lsq_fit(O)
    print(f"least squares (InSAR only, for the noise model): east {qhat[0]:+.2f} north {qhat[1]:+.2f} km, depth {qhat[2]:.2f} km, "
          f"dV {10 ** qhat[3] * 1e3:.1f}e6 m3/yr, residual rms {np.sqrt(np.mean(resid ** 2)):.2f} cm/yr")
    noise = exp_cov_fit(O, resid)
    t0 = time.time()
    q, x = simulate(O, noise, n_sims)
    print(f"simulated {n_sims} data sets in {time.time() - t0:.0f} s")
    xobs = np.concatenate([np.concatenate([O.tr[t]["y"] for t in M.TRACKS]), O.gv])[None].astype(np.float32)
    nin = O.n_ins
    models = {"NISAR + GNSS": train(x, q, np.arange(x.shape[1]), "NISAR + GNSS"),
              "NISAR only": train(x, q, np.arange(nin), "NISAR only")}
    res = {}
    qt, xt = simulate(O, noise, 600)
    for label, mdl in models.items():
        cov = calibration(mdl, xt, qt)
        print(f"calibration [{label}]: 68 % intervals contain the truth "
              f"{', '.join(f'{c:.0%}' for c in cov[0.68])}; 95 %: {', '.join(f'{c:.0%}' for c in cov[0.95])} "
              f"(east, north, depth, volume)")
        t1 = time.time()
        s = posterior(mdl, xobs)
        print(f"posterior for the real data in {time.time() - t1:.2f} s:")
        vals, qs = summarize(s, label)
        res[label] = (vals, qs, cov)
    np.savez(ML / "npe_posterior.npz", **{k.replace(" ", "_").replace("+", "and"): v[0] for k, v in res.items()})
    plot(res, D)


def plot(res, D):
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))
    cols = {"NISAR + GNSS": "tab:red", "NISAR only": "tab:blue"}
    for label, (v, qs, cov) in res.items():
        c = cols[label]
        ax[0].scatter(v[:4000, 2], v[:4000, 3], s=2, alpha=.15, color=c)
        ax[0].scatter([], [], color=c, label=label)
        lon, lat = D.lonlat(v[:4000, 0], v[:4000, 1])
        ax[1].scatter(lon, lat, s=2, alpha=.15, color=c)
        ax[2].hist(v[:, 2], bins=60, range=(1, 10), histtype="step", color=c, density=True, label=label)
    ax[0].set_xlabel("depth, km"), ax[0].set_ylabel("volume rate, 10^6 m^3/yr"), ax[0].legend(fontsize=8, frameon=False)
    ax[0].set_title("posterior: depth vs volume (the trade-off)", fontsize=10)
    ax[1].set_title("posterior: source position", fontsize=10), ax[1].set_xlabel("lon"), ax[1].set_ylabel("lat")
    ax[2].axvspan(3.5, 4.5, color="0.85", label="IMO ~4 km")
    ax[2].set_title("posterior of depth", fontsize=10), ax[2].set_xlabel("km"), ax[2].legend(fontsize=8, frameon=False)
    fig.suptitle(f"Neural posterior estimation of the Svartsengi source (Mogi), trained on {N_SIMS} simulations", fontsize=11)
    fig.tight_layout()
    fig.savefig(ML / "npe_posterior.png", dpi=110, facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 80000)


CACHE = ML / "obs_cache.npz"


def save_obs_cache(O, noise):
    """Store the binned observation vector, geometry and noise model, so the inference scripts
    run without the raw GUNW granules (python -c "import ml_npe_source as m; m.build_cache()")."""
    out = {}
    for t in M.TRACKS:
        for k in ("e", "n", "y", "vec"):
            out[f"{t}_{k}"] = O.tr[t][k]
        out[f"{t}_chol"] = noise[t]["chol"]
        out[f"{t}_sigma"], out[f"{t}_L"] = noise[t]["sigma"], noise[t]["L"]
    out.update(gn=O.gn, gref=O.gref, gv=O.gv)
    np.savez(CACHE, **out)
    print("wrote", CACHE)


def build_cache():
    D = Data()
    O = Obs(D)
    _, resid = lsq_fit(O)
    save_obs_cache(O, exp_cov_fit(O, resid))


def cached(cls):
    """An Obs (or subclass) instance filled from obs_cache.npz instead of the raw data, plus the
    Gaussian noise model of the NPE: returns (O, noise)."""
    z = np.load(CACHE)
    O = cls.__new__(cls)
    O.tr = {t: {k: z[f"{t}_{k}"] for k in ("e", "n", "y", "vec")} for t in M.TRACKS}
    O.gn, O.gref, O.gv = z["gn"], z["gref"], z["gv"]
    O.n_ins = sum(O.tr[t]["y"].size for t in M.TRACKS)
    if hasattr(cls, "post_init"):
        O.post_init()
    noise = {t: dict(chol=z[f"{t}_chol"], sigma=float(z[f"{t}_sigma"]), L=float(z[f"{t}_L"])) for t in M.TRACKS}
    return O, noise
