"""Neural posterior estimation of the Svartsengi source as an elongated horizontal body
(compound dislocation model: a horizontal and a vertical opening crack along one axis),
from NISAR LOS rates (159A, 082D) with and without GNSS.

Parameters (priors uniform unless noted):
  east, north   centre, km from Svartsengi            [-6, 6]
  depth         centre depth, km                      [1.5, 8]
  half_len      half-length a of the body, km         [0.3, 6]   (a small = point-like)
  strike        direction of the long axis, deg       [0, 180)
  log10 V       volume rate, km^3/yr                  [2e-3, 8e-2]  (2-80 million m^3/yr)
  w             share of V in the horizontal crack    [0, 1]  (1 = sill-like, 0 = dike-like)
The cross-section is fixed (semi-axes 0.25 km): the data cannot resolve it.
Nuisance per track (offset, plane) and correlated noise as in ml_npe_source.py.
The network models logit-transformed parameters, so posterior samples respect the prior bounds.

Usage: python scripts/reykjanes/ml_npe_ellipsoid.py [n_sims=40000]
"""
import sys
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import torch

import reykjanes_mogi as M
from okada import tensile
from ml_npe_source import MDN, Obs, exp_cov_fit, lsq_fit
from reykjanes_joint import SIGMA_GNSS
from reykjanes_sources import Data

ROOT = Path(__file__).resolve().parents[2]
ML = ROOT / "output" / "reykjanes" / "ml"
B = C = 0.25  # km, fixed semi-axes of the cross-section
NAMES = ("east km", "north km", "depth km", "half-length km", "strike deg", "log10 V", "sill share w")
LO = np.array([-6, -6, 1.5, 0.3, 0.0, np.log10(2e-3), 0.0])
HI = np.array([6, 6, 8.0, 6.0, 180.0, np.log10(8e-2), 1.0])
RNG = np.random.default_rng(3)
EPS = 1e-4


def to_z(q):
    """Bounded parameters -> unbounded: logit of the position inside the prior box. The network
    models z, so no posterior sample can fall outside the prior bounds after mapping back."""
    u = np.clip((q - LO) / (HI - LO), EPS, 1 - EPS)
    return np.log(u / (1 - u)).astype(np.float32)


def from_z(z):
    return LO + (HI - LO) / (1 + np.exp(-z))
torch.manual_seed(3)


class EllObs(Obs):
    def __init__(self, D):
        super().__init__(D)
        self.post_init()

    def post_init(self):
        self.pe = np.concatenate([self.tr[t]["e"] for t in M.TRACKS] + [self.gn[:, 0], self.gref[:1]])
        self.pn = np.concatenate([self.tr[t]["n"] for t in M.TRACKS] + [self.gn[:, 1], self.gref[1:]])
        self.sizes = [self.tr[t]["y"].size for t in M.TRACKS]

    def forward_one(self, q, nuis):
        e0, n0, d, a, strike, lv, w = q
        V = 10 ** lv * 1e9                                   # m^3/yr
        Uz = w * V / (4 * a * B * 1e6)                        # m/yr, horizontal crack 2a x 2b
        Uy = (1 - w) * V / (4 * a * C * 1e6)                  # vertical crack 2a x 2c along the axis
        u = tensile(self.pe, self.pn, e0, n0, d, 2 * a, 2 * B, strike, 0.0, Uz) \
            + tensile(self.pe, self.pn, e0, n0, d, 2 * a, 2 * C, strike, 90.0, Uy)
        u = u * 100.0                                         # cm/yr, (points, 3)
        out, i0 = [], 0
        for k, t in enumerate(M.TRACKS):
            n = self.sizes[k]
            dd = self.tr[t]
            los = (u[i0:i0 + n] * dd["vec"]).sum(1)
            los += nuis[3 * k] + nuis[3 * k + 1] * dd["e"] / 10 + nuis[3 * k + 2] * dd["n"] / 10
            out.append(los)
            i0 += n
        g = u[i0:i0 + len(self.gn)] - u[-1]                   # GNSS relative to NYLA
        out.append(g.ravel())
        return np.concatenate(out)


def simulate(O, noise, N):
    q = LO + RNG.random((N, 7)) * (HI - LO)
    nuis = np.c_[RNG.uniform(-10, 10, N), RNG.normal(0, 2, N), RNG.normal(0, 2, N),
                 RNG.uniform(-10, 10, N), RNG.normal(0, 2, N), RNG.normal(0, 2, N)]
    x = np.empty((N, sum(O.sizes) + 3 * len(O.gn)), np.float32)
    t0 = time.time()
    for i in range(N):
        x[i] = O.forward_one(q[i], nuis[i])
        if i == 199:
            print(f"  {(time.time() - t0) / 200 * 1e3:.1f} ms per simulation -> ~{(time.time() - t0) / 200 * N / 60:.1f} min for {N}")
    i0 = 0
    for t in M.TRACKS:
        n = O.tr[t]["y"].size
        x[:, i0:i0 + n] += (RNG.standard_normal((N, n)) @ noise[t]["chol"].T).astype(np.float32)
        i0 += n
    x[:, i0:] += RNG.normal(0, SIGMA_GNSS, (N, x.shape[1] - i0)).astype(np.float32)
    return q.astype(np.float32), x


def train(x, q, cols, label, epochs=25):
    xs = x[:, cols]
    xm, xsd = xs.mean(0), xs.std(0) + 1e-6
    qm, qsd = q.mean(0), q.std(0)
    X, Q = torch.tensor((xs - xm) / xsd), torch.tensor((q - qm) / qsd)
    nval = len(X) // 10
    net = MDN(X.shape[1], nout=q.shape[1])
    opt = torch.optim.Adam(net.parameters(), lr=1e-3)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, epochs)
    best, state, t0 = 1e9, None, time.time()
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
            best, state = v, {k: w.clone() for k, w in net.state_dict().items()}
        if ep % 5 == 4:
            print(f"  [{label}] epoch {ep + 1}: validation -log p {v:.3f} ({time.time() - t0:.0f} s)")
    net.load_state_dict(state)
    return dict(net=net, xm=xm, xsd=xsd, qm=qm, qsd=qsd, cols=cols)


def post(m, x, n=20000):
    """Posterior samples in the physical (bounded) parameters."""
    with torch.no_grad():
        X = torch.tensor(((x[:, m["cols"]] - m["xm"]) / m["xsd"]).astype(np.float32))
        s = m["net"].dist(X).sample((n,)).numpy()
    return from_z(s * m["qsd"] + m["qm"])


def main(n_sims=40000):
    ML.mkdir(parents=True, exist_ok=True)
    D = Data()
    O = EllObs(D)
    _, resid = lsq_fit(O)
    noise = exp_cov_fit(O, resid)
    print(f"simulating {n_sims} elongated sources ...")
    q, x = simulate(O, noise, n_sims)
    qt, xt = simulate(O, noise, 600)
    xobs = np.concatenate([np.concatenate([O.tr[t]["y"] for t in M.TRACKS]), O.gv])[None].astype(np.float32)
    nin = sum(O.sizes)
    results = {}
    for label, cols in (("NISAR + GNSS", np.arange(x.shape[1])), ("NISAR only", np.arange(nin))):
        mdl = train(x, to_z(q), cols, label)
        s = post(mdl, xt, 2000)
        cov = {lev: ((qt >= np.quantile(s, (1 - lev) / 2, 0)) & (qt <= np.quantile(s, (1 + lev) / 2, 0))).mean(0) for lev in (.68, .95)}
        print(f"calibration [{label}] 68 %: " + ", ".join(f"{c:.0%}" for c in cov[.68]) + " | 95 %: " + ", ".join(f"{c:.0%}" for c in cov[.95]))
        sr = post(mdl, xobs)[:, 0]
        inside = np.all((sr >= LO) & (sr <= HI), axis=1).mean()
        print(f"   samples inside the prior bounds: {inside:.1%}")
        vals = sr.copy()
        vals[:, 5] = 10 ** vals[:, 5] * 1e3                  # -> million m^3/yr
        print(f"posterior [{label}] for the real data:")
        for i, nm in enumerate(NAMES):
            nm2 = "V 1e6 m3/yr" if i == 5 else nm
            qq = np.percentile(vals[:, i], [2.5, 16, 50, 84, 97.5])
            print(f"   {nm2:15s} median {qq[2]:7.2f}   68 % {qq[1]:7.2f} .. {qq[3]:7.2f}   95 % {qq[0]:7.2f} .. {qq[4]:7.2f}")
        # how sure are we that the body is elongated, and in which direction?
        el = vals[:, 3] > 1.5
        print(f"   P(half-length > 1.5 km, i.e. longer than 3 km) = {el.mean():.0%}; "
              f"strike among those: median {np.median(vals[el, 4]):.0f} deg, 68 % "
              f"{np.percentile(vals[el, 4], 16):.0f}..{np.percentile(vals[el, 4], 84):.0f}")
        results[label] = (vals, cov)
    np.savez(ML / "npe_ellipsoid.npz", **{k.replace(" ", "_").replace("+", "and"): v[0] for k, v in results.items()})
    plot(results, D, n_sims)


def plot(results, D, n_sims):
    fig, ax = plt.subplots(1, 4, figsize=(19, 4.6))
    col = {"NISAR + GNSS": "tab:red", "NISAR only": "tab:blue"}
    for label, (v, cov) in results.items():
        c = col[label]
        ax[0].hist(v[:, 3] * 2, bins=60, range=(0, 12), histtype="step", density=True, color=c, label=label)
        el = v[:, 3] > 1.5
        ax[1].hist(v[el, 4], bins=36, range=(0, 180), histtype="step", density=True, color=c, label=label)
        ax[2].scatter(v[:4000, 2], v[:4000, 5], s=2, alpha=.15, color=c)
        ax[3].hist(v[:, 6], bins=40, range=(0, 1), histtype="step", density=True, color=c, label=label)
    ax[0].set_title("length of the body, km (2a)", fontsize=10); ax[0].set_xlabel("km"); ax[0].legend(fontsize=8, frameon=False)
    ax[1].set_title("strike of the long axis (when > 3 km long)", fontsize=10); ax[1].set_xlabel("deg from north")
    ax[2].set_title("depth vs volume rate", fontsize=10); ax[2].set_xlabel("depth, km"); ax[2].set_ylabel("10^6 m^3/yr")
    ax[3].set_title("share of the volume in horizontal opening (1 = sill)", fontsize=10)
    fig.suptitle(f"Neural posterior estimation, elongated source (CDM), {n_sims} simulations, logit-bounded parameters", fontsize=11)
    fig.tight_layout()
    fig.savefig(ML / "npe_ellipsoid.png", dpi=110, facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 40000)
