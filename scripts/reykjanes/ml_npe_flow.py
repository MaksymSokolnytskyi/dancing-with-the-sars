"""A neural posterior estimator built to match MCMC, for the same problem and likelihood as
ml_mcmc_source.py --noise screens (exponential + white covariance fitted to the real screen
library, nuisance offset / plane marginalised as Gaussian).

What changes against the first NPE (ml_npe_source / ml_npe_ellipsoid):
  1. physics-aware input: every track is whitened with the same matrix the MCMC likelihood uses
     (noise covariance + nuisance offset/plane), GNSS divided by its sigma, and the whitened
     vector projected on the K leading principal components of the whitened noise-free signals.
     In these coordinates the noise is exactly N(0, I): the network sees a K-dim, well-scaled
     summary instead of 713 raw numbers with correlated noise.
  2. fresh noise every epoch on the same noise-free simulations (infinite augmentation).
  3. a conditional neural spline flow (zuko NSF) in logit-bounded parameters instead of an
     8-component Gaussian mixture.
  4. a second, truncated round (TSNPE): simulations only inside the box that the first-round
     posterior of the real data occupies (99.9 % range, widened by 25 %), prior otherwise unchanged.
  5. an ensemble of 3 flows in round 2 (samples pooled).

Compared with the MCMC posterior (mcmc_<kind>_screens.npz): shift of every median in MCMC
standard deviations, width ratio, and calibration on held-out simulations from the truncated box.

Usage: python scripts/reykjanes/ml_npe_flow.py {mogi|ellipsoid} [--data gnss|nisar] [--n1 60000] [--n2 120000]
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
import torch
import zuko

import reykjanes_mogi as M
import ml_npe_source as NS
import ml_npe_ellipsoid as NE
from ml_mcmc_source import SD_OFFSET, SD_PLANE, Likelihood
from reykjanes_joint import SIGMA_GNSS

ML = Path(__file__).resolve().parents[2] / "output" / "reykjanes" / "ml"
K_PCA = 48
torch.manual_seed(7)
RNG = np.random.default_rng(7)
torch.set_num_threads(int(__import__("os").environ.get("NPE_THREADS", "8")))
BATCH = int(__import__("os").environ.get("NPE_BATCH", "2048"))   # big batches: far fewer, cheaper steps per epoch on a CPU


def arg(name, default):
    return type(default)(sys.argv[sys.argv.index(name) + 1]) if name in sys.argv else default


class Problem:
    def __init__(self, kind, use_gnss):
        self.kind, self.use_gnss = kind, use_gnss
        if kind == "mogi":
            self.O, _ = NS.cached(NS.Obs)
            self.lo = np.array([NS.PRIOR[k][0] for k in ("e", "n", "depth", "logv")])
            self.hi = np.array([NS.PRIOR[k][1] for k in ("e", "n", "depth", "logv")])
            self.names = ("east km", "north km", "depth km", "dV 1e6 m3/yr")
            self.vcol = 3
        else:
            self.O, _ = NS.cached(NE.EllObs)
            self.lo, self.hi = NE.LO.copy(), NE.HI.copy()
            self.names = NE.NAMES[:5] + ("V 1e6 m3/yr", "sill share w")
            self.vcol = 5
        z = np.load(ML / "screens.npz")
        cov = {t: z[f"cov_par_{t}"] for t in M.TRACKS}
        self.like = Likelihood(self.O, cov, use_gnss)              # whitening matrices, as in MCMC
        self.sizes = [self.O.tr[t]["y"].size for t in M.TRACKS]

    def clean(self, q):
        """Noise-free simulator output (no nuisance: it is inside the whitening covariance)."""
        if self.kind == "mogi":
            p = q.copy(); p[:, 3] = 10 ** q[:, 3]
            return self.O.forward(p)
        z6 = np.zeros(6)
        from progress import Progress
        with Progress(f"симуляції Окади ({self.kind})", total=len(q), unit="сим") as p:
            out = []
            for i, qi in enumerate(q):
                out.append(self.O.forward_one(qi, z6))
                if i % 200 == 199:
                    p.update(200)
        return np.array(out)

    def whiten(self, x):
        out, i0 = [], 0
        for n, W in zip(self.sizes, self.like.W):
            out.append(x[:, i0:i0 + n] @ W.T)
            i0 += n
        if self.use_gnss:
            out.append(x[:, i0:] / SIGMA_GNSS)
        return np.concatenate(out, 1)

    def observed(self):
        x = np.concatenate([np.concatenate([self.O.tr[t]["y"] for t in M.TRACKS]), self.O.gv])[None]
        return self.whiten(x)


def to_z(q, lo, hi):
    u = np.clip((q - lo) / (hi - lo), 1e-5, 1 - 1e-5)
    return np.log(u / (1 - u))


def from_z(z, lo, hi):
    return lo + (hi - lo) / (1 + np.exp(-z))


class Net(torch.nn.Module):
    def __init__(self, nin, npar):
        super().__init__()
        self.emb = torch.nn.Sequential(torch.nn.Linear(nin, 256), torch.nn.GELU(), torch.nn.Linear(256, 256), torch.nn.GELU(),
                                       torch.nn.Linear(256, 128))
        self.flow = zuko.flows.NSF(features=npar, context=128, transforms=5, hidden_features=(192, 192))

    def log_prob(self, z, x):
        return self.flow(self.emb(x)).log_prob(z)

    def sample(self, x, n):
        with torch.no_grad():
            return self.flow(self.emb(x)).sample((n,))


def train(S, Q, box, epochs, label, seed, scale=None):
    """S: whitened-projected noise-free signals (N, K); Q: logit params (N, d). Fresh N(0, I) noise each epoch.
    scale: optional per-component divisor of the network input (large signals, e.g. dikes)."""
    torch.manual_seed(seed)
    N, K = S.shape
    nval = N // 10
    St, Qt = torch.tensor(S, dtype=torch.float32), torch.tensor(Q, dtype=torch.float32)
    gval = torch.Generator().manual_seed(123)
    sc = torch.ones(K) if scale is None else torch.tensor(scale, dtype=torch.float32)
    Xv = (St[:nval] + torch.randn(nval, K, generator=gval)) / sc
    net = Net(K, Q.shape[1])
    opt = torch.optim.AdamW(net.parameters(), lr=2e-3, weight_decay=1e-5)
    sched = torch.optim.lr_scheduler.OneCycleLR(opt, max_lr=2e-3, total_steps=epochs * ((N - nval) // BATCH + 1))
    best, state, t0 = 1e9, None, time.time()
    from progress import Progress
    prog = Progress(f"навчання: {label}", total=epochs, unit="епоха")
    for ep in range(epochs):
        net.train()
        perm = torch.randperm(N - nval) + nval
        for i in range(0, len(perm), BATCH):
            b = perm[i:i + BATCH]
            x = (St[b] + torch.randn(len(b), K)) / sc            # fresh noise every time
            loss = -net.log_prob(Qt[b], x).mean()
            opt.zero_grad(); loss.backward()
            torch.nn.utils.clip_grad_norm_(net.parameters(), 5.0)
            opt.step(); sched.step()
        net.eval()
        with torch.no_grad():
            v = -net.log_prob(Qt[:nval], Xv).mean().item()
        if v < best:
            best, state = v, {k: w.clone() for k, w in net.state_dict().items()}
        prog.update(1, note=f"перевірочна -log p {v:.3f} (найкраща {best:.3f})")
        if ep % 10 == 9 or ep == epochs - 1:
            print(f"  [{label}] epoch {ep + 1}/{epochs}: validation -log p {v:.3f} ({time.time() - t0:.0f} s)", flush=True)
    net.load_state_dict(state)
    prog.close()
    return net, time.time() - t0


def posterior(nets, xk, lo, hi, n=20000, scale=None):
    X = torch.tensor(xk if scale is None else xk / scale, dtype=torch.float32)
    per = n // len(nets)
    s = np.concatenate([net.sample(X, per).numpy()[:, 0] for net in nets])
    return from_z(s.astype(float), lo, hi)


def physical(P, v):
    v = v.copy()
    v[:, P.vcol] = 10 ** v[:, P.vcol] * 1e3
    return v


def main():
    kind = sys.argv[1]
    data = arg("--data", "gnss")
    use_gnss = data == "gnss"
    label = "NISAR + GNSS" if use_gnss else "NISAR only"
    n1, n2 = arg("--n1", 60000), arg("--n2", 120000)
    ep1, ep2 = arg("--ep1", 40), arg("--ep2", 60)
    P = Problem(kind, use_gnss)
    lo, hi = P.lo, P.hi
    K_PCA = arg("--k", globals()["K_PCA"])     # number of kept components (signal must be << noise beyond it)
    T = {}
    # ---- round 1: prior
    t0 = time.time()
    q1 = lo + RNG.random((n1, len(lo))) * (hi - lo)
    X1 = P.whiten(P.clean(q1))
    T["sim1"] = time.time() - t0
    # PCA basis of the whitened noise-free signals (noise stays N(0, I) after an orthonormal projection)
    Xc = X1 - X1.mean(0)
    _, sv, Vt = np.linalg.svd(Xc[:20000], full_matrices=False)
    U = Vt[:K_PCA].T
    kept = (sv[:K_PCA] ** 2).sum() / (sv ** 2).sum()
    print(f"[{kind} {label}] {n1} prior simulations in {T['sim1']:.0f} s; {K_PCA} components keep {kept:.4%} of the signal variance; "
          f"component {K_PCA} still has {sv[K_PCA - 1] / np.sqrt(len(Xc[:20000])):.2f} x the noise")
    xo = P.observed() @ U
    net1, T["train1"] = train(X1 @ U, to_z(q1, lo, hi), None, ep1, f"{kind} round 1", 1)
    s1 = posterior([net1], xo, lo, hi)
    # ---- round 2: truncated box
    a, b = np.quantile(s1, [0.0005, 0.9995], 0)
    w = b - a
    lo2, hi2 = np.maximum(lo, a - 0.25 * w), np.minimum(hi, b + 0.25 * w)
    if kind == "ellipsoid":                      # strike is periodic: keep the full range if the box touches an end
        if lo2[4] <= lo[4] + 1 or hi2[4] >= hi[4] - 1:
            lo2[4], hi2[4] = lo[4], hi[4]
    print("   truncated box: " + ", ".join(f"{nm.split()[0]} {l:.2f}..{h:.2f}" for nm, l, h in zip(P.names, lo2, hi2)))
    t0 = time.time()
    q2 = lo2 + RNG.random((n2, len(lo))) * (hi2 - lo2)
    X2 = P.whiten(P.clean(q2)) @ U
    T["sim2"] = time.time() - t0
    nets, T["train2"] = [], 0.0
    for k in range(arg("--ens", 1)):
        net, tt = train(X2, to_z(q2, lo2, hi2), None, ep2, f"{kind} round 2 #{k + 1}", 10 + k)
        nets.append(net)
        T["train2"] += tt
    t0 = time.time()
    s2 = posterior(nets, xo, lo2, hi2)
    T["infer"] = time.time() - t0
    # ---- calibration on held-out simulations from the truncated box
    qt = lo2 + RNG.random((400, len(lo))) * (hi2 - lo2)
    xt = (P.whiten(P.clean(qt)) @ U) + RNG.standard_normal((400, K_PCA))
    # sampled in chunks of 25 test sets: all 400 at once needs > 600 MB inside the flow
    st = np.concatenate([np.concatenate([from_z(net.sample(torch.tensor(xt[j:j + 25], dtype=torch.float32), 2100 // len(nets)).numpy().astype(float), lo2, hi2)
                                         for net in nets]) for j in range(0, 400, 25)], axis=1)   # (~2100, 400, d)
    cal = {lev: ((qt >= np.quantile(st, (1 - lev) / 2, 0)) & (qt <= np.quantile(st, (1 + lev) / 2, 0))).mean(0).round(3).tolist()
           for lev in (.68, .95)}
    # ---- compare with MCMC
    mc_file = ML / f"mcmc_{kind}_screens.npz"
    key = label.replace(" ", "_").replace("+", "and")
    mc = np.load(mc_file)[key] if mc_file.exists() else None
    v1, v2 = physical(P, s1), physical(P, s2)
    old_file = ML / ("npe_screens_mogi.npz" if kind == "mogi" else "npe_screens_ellipsoid.npz")
    old = np.load(old_file)[f"raw__{key}"] if old_file.exists() else None
    rows = {}
    print(f"[{kind} {label}] calibration 68 %: {cal[.68]}; 95 %: {cal[.95]}")
    for i, nm in enumerate(P.names):
        r = {}
        for tag, v in (("flow_round1", v1), ("flow_round2", v2), ("npe_old", old), ("mcmc", mc)):
            if v is not None:
                r[tag] = np.round(np.percentile(v[:, i], [2.5, 16, 50, 84, 97.5]), 3).tolist()
        if mc is not None:
            m = r["mcmc"]
            sd = (m[3] - m[1]) / 2
            for tag in ("flow_round1", "flow_round2", "npe_old"):
                if tag in r:
                    r[f"{tag}_shift_sd"] = round((r[tag][2] - m[2]) / sd, 2)
                    r[f"{tag}_width_ratio"] = round((r[tag][3] - r[tag][1]) / (m[3] - m[1]), 2)
        rows[nm] = r
        print(f"   {nm:15s} MCMC {r.get('mcmc', [0]*5)[2]:7.2f} [{r.get('mcmc', [0]*5)[1]:.2f}..{r.get('mcmc', [0]*5)[3]:.2f}] | "
              f"flow r2 {r['flow_round2'][2]:7.2f} [{r['flow_round2'][1]:.2f}..{r['flow_round2'][3]:.2f}] "
              f"shift {r.get('flow_round2_shift_sd', 0):+.2f} sd, width x{r.get('flow_round2_width_ratio', 0):.2f} | "
              f"r1 shift {r.get('flow_round1_shift_sd', 0):+.2f}, x{r.get('flow_round1_width_ratio', 0):.2f} | "
              f"old NPE shift {r.get('npe_old_shift_sd', 0):+.2f}, x{r.get('npe_old_width_ratio', 0):.2f}")
    extra = {}
    if kind == "ellipsoid":
        extra = {tag: round(float((v[:, 3] > 1.5).mean()), 3) for tag, v in (("flow_round2", v2), ("mcmc", mc), ("npe_old", old)) if v is not None}
        print(f"   P(longer than 3 km): {extra}")
    T = {k: round(v, 2) for k, v in T.items()}
    print(f"   time: simulations {T['sim1'] + T['sim2']:.0f} s, training {T['train1'] + T['train2']:.0f} s, "
          f"posterior for new data {T['infer']:.2f} s")
    out = dict(kind=kind, data=label, calibration=cal, params=rows, p_long=extra, timing=T, n1=n1, n2=n2, k_pca=K_PCA,
               box=[lo2.tolist(), hi2.tolist()])
    tag = f"npe_flow_{kind}_{data}"
    (ML / f"{tag}.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    np.savez(ML / f"{tag}.npz", **{key: v2})


if __name__ == "__main__":
    main()
