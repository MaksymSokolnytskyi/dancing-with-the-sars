"""MCMC cross-check of the neural posterior (NPE) for the Svartsengi source.

Same physics (Mogi point source, or the elongated compound-dislocation body), same binned
NISAR LOS rates and GNSS velocities, same priors and the same noise model as
ml_npe_source.py / ml_npe_ellipsoid.py, but the posterior is sampled directly with an
affine-invariant ensemble sampler (Goodman & Weare 2010 stretch move) instead of being
learned by a network.

Likelihood: per track, Gaussian with covariance C_t + A_t P A_t^T, where C_t is the
exponential noise covariance fitted to least-squares residuals (as in the NPE simulator)
and A_t P A_t^T marginalises the nuisance offset and plane analytically (offset sd 5.77 cm/yr
= the sd of the NPE's uniform(-10, 10) prior, plane sd 2 cm/yr per 10 km as in the NPE).
GNSS: independent, sd SIGMA_GNSS.

Options:
  --noise gauss      exponential covariance from residuals (default; identical to the NPE)
  --noise screens    exponential + white covariance fitted to the real atmospheric screen library
                     (output/reykjanes/ml/screens.npz from ml_screens.py)

Usage: python scripts/reykjanes/ml_mcmc_source.py {mogi|ellipsoid|both} [--steps N] [--walkers W] [--noise gauss|screens] [--data gnss|nisar]
"""
import json
import sys
import time
from pathlib import Path

import numpy as np

import reykjanes_mogi as M
from ml_npe_source import PRIOR, Obs, exp_cov_fit, lsq_fit
from ml_npe_ellipsoid import HI as E_HI, LO as E_LO, EllObs
from reykjanes_joint import SIGMA_GNSS
from reykjanes_sources import Data

ML = Path(__file__).resolve().parents[2] / "output" / "reykjanes" / "ml"
SD_OFFSET, SD_PLANE = 20 / np.sqrt(12), 2.0
M_LO = np.array([PRIOR[k][0] for k in ("e", "n", "depth", "logv")])
M_HI = np.array([PRIOR[k][1] for k in ("e", "n", "depth", "logv")])


def arg(name, default):
    return type(default)(sys.argv[sys.argv.index(name) + 1]) if name in sys.argv else default


class Likelihood:
    def __init__(self, O, cov, use_gnss):
        self.O, self.use_gnss = O, use_gnss
        self.y = [O.tr[t]["y"] for t in M.TRACKS]
        self.W = []
        for t in M.TRACKS:
            d = O.tr[t]
            A = np.c_[np.ones_like(d["e"]), d["e"] / 10, d["n"] / 10]
            S = cov[t] + A @ np.diag([SD_OFFSET ** 2, SD_PLANE ** 2, SD_PLANE ** 2]) @ A.T
            self.W.append(np.linalg.inv(np.linalg.cholesky(S)))   # whitening matrix
        self.n_ins = sum(v.size for v in self.y)

    def __call__(self, pred):
        """pred: (N, n_ins + 12) noise-free simulator output without nuisances -> log-likelihood (N,)."""
        ll, i0 = np.zeros(len(pred)), 0
        for y, W in zip(self.y, self.W):
            n = y.size
            r = (pred[:, i0:i0 + n] - y) @ W.T
            ll -= 0.5 * (r ** 2).sum(1)
            i0 += n
        if self.use_gnss:
            ll -= 0.5 * (((pred[:, i0:] - self.O.gv) / SIGMA_GNSS) ** 2).sum(1)
        return ll


class MogiModel:
    names = ("east km", "north km", "depth km", "dV 1e6 m3/yr")
    lo, hi = M_LO, M_HI

    def __init__(self, O):
        self.O = O

    def predict(self, q):
        p = q.copy()
        p[:, 3] = 10 ** q[:, 3]
        return self.O.forward(p)

    @staticmethod
    def physical(q):
        v = q.copy()
        v[:, 3] = 10 ** v[:, 3] * 1e3
        return v

    def wrap(self, q):
        return q


class EllipsoidModel:
    names = ("east km", "north km", "depth km", "half-length km", "strike deg", "V 1e6 m3/yr", "sill share w")
    lo, hi = E_LO, E_HI

    def __init__(self, O):
        self.O = O

    def predict(self, q):
        z = np.zeros(6)
        return np.array([self.O.forward_one(qi, z) for qi in q])

    @staticmethod
    def physical(q):
        v = q.copy()
        v[:, 5] = 10 ** v[:, 5] * 1e3
        return v

    def wrap(self, q):
        q = q.copy()
        q[:, 4] = np.mod(q[:, 4], 180.0)   # a horizontal or vertical crack rotated by 180 deg is the same body
        return q


def log_post(model, like, q):
    inside = np.all((q >= model.lo) & (q <= model.hi), axis=1)
    out = np.full(len(q), -np.inf)
    if inside.any():
        out[inside] = like(model.predict(q[inside]))
    return out


def stretch_sampler(model, like, start, steps, walkers, rng, a=2.0, label=""):
    """Goodman & Weare affine-invariant ensemble sampler, two half-ensembles updated in turn."""
    ndim = start.shape[1]
    x = start.copy()
    lp = log_post(model, like, x)
    chain = np.empty((steps, walkers, ndim))
    lps = np.empty((steps, walkers))
    acc = 0
    half = walkers // 2
    t0 = time.time()
    from progress import Progress
    prog = Progress(f"MCMC {label}", total=steps, unit="крок", note=f"{walkers} ланцюжків")
    for s in range(steps):
        for first in (True, False):
            S = np.arange(half) if first else np.arange(half, walkers)
            Cc = np.arange(half, walkers) if first else np.arange(half)
            z = ((a - 1) * rng.random(len(S)) + 1) ** 2 / a
            partner = x[rng.choice(Cc, len(S))]
            prop = model.wrap(partner + z[:, None] * (x[S] - partner))
            lpp = log_post(model, like, prop)
            lnr = (ndim - 1) * np.log(z) + lpp - lp[S]
            ok = np.log(rng.random(len(S))) < lnr
            x[S[ok]], lp[S[ok]] = prop[ok], lpp[ok]
            acc += ok.sum()
        chain[s], lps[s] = x, lp
        prog.update(1, note=f"{walkers} ланцюжків, прийнято {acc / ((s + 1) * walkers):.0%}, найкраща log p {lp.max():.1f}")
        if (s + 1) % max(steps // 5, 1) == 0:
            print(f"  [{label}] step {s + 1}/{steps}, acceptance {acc / ((s + 1) * walkers):.2f} ({time.time() - t0:.0f} s)", flush=True)
    prog.close()
    return chain, lps, acc / (steps * walkers)


def autocorr_time(chain):
    """Integrated autocorrelation time per parameter, averaged over walkers (Sokal's window, c=5)."""
    n, w, d = chain.shape
    out = []
    for k in range(d):
        f = np.zeros(n)
        for j in range(w):
            x = chain[:, j, k] - chain[:, j, k].mean()
            ft = np.fft.rfft(x, 2 * n)
            ac = np.fft.irfft(ft * np.conj(ft))[:n]
            f += ac / ac[0] if ac[0] > 0 else 0
        f /= w
        taus = 2 * np.cumsum(f) - 1
        m = np.arange(n) < 5 * taus
        win = np.argmin(m) if not m.all() else n - 1
        out.append(float(taus[win]))
    return np.array(out)


def rhat(chain):
    """Gelman-Rubin on the two halves of every walker's chain (split R-hat over walkers)."""
    n = chain.shape[0] // 2
    c = np.concatenate([chain[:n], chain[n:2 * n]], 1)           # (n, 2w, d)
    m = c.mean(0)
    B = n * m.var(0, ddof=1)
    Wv = c.var(0, ddof=1).mean(0)
    return np.sqrt(((n - 1) / n * Wv + B / n) / Wv)


ONLY = None


def run(kind, noise_kind, steps, walkers):
    D = Data()
    O = (Obs if kind == "mogi" else EllObs)(D)
    _, resid = lsq_fit(O)
    noise = exp_cov_fit(O, resid)
    if noise_kind == "gauss":
        cov = {t: noise[t]["chol"] @ noise[t]["chol"].T for t in M.TRACKS}
    else:
        z = np.load(ML / "screens.npz")
        cov = {t: z[f"cov_par_{t}"] for t in M.TRACKS}
        print("noise: exponential + white covariance fitted to the real screen library (raw)")
    model = (MogiModel if kind == "mogi" else EllipsoidModel)(O)
    npe_file = ML / ("npe_posterior.npz" if kind == "mogi" else "npe_ellipsoid.npz")
    npe = np.load(npe_file) if npe_file.exists() else None
    rng = np.random.default_rng(11)
    results = {}
    configs = [c for c in (("NISAR + GNSS", True), ("NISAR only", False)) if ONLY is None or c[0] == ONLY]
    prev = ML / f"mcmc_{kind}{'' if noise_kind == 'gauss' else '_screens'}.json"
    keep_prev = json.loads(prev.read_text(encoding="utf-8")) if ONLY and prev.exists() else {}
    for label, use_gnss in configs:
        like = Likelihood(O, cov, use_gnss)
        # start the walkers from a prior draw filtered to the best of many, so nothing is assumed
        cand = model.lo + rng.random((4000 if kind == "mogi" else 1500, len(model.lo))) * (model.hi - model.lo)
        lp = log_post(model, like, cand)
        best = cand[np.argsort(lp)[-walkers:]]
        t0 = time.time()
        chain, lps, acc = stretch_sampler(model, like, best, steps, walkers, rng, label=f"{kind} {label}")
        burn = steps // 3
        post = chain[burn:]
        # walkers stranded in a far, much less probable mode (a known weakness of stretch moves):
        # dropped when their mean log-posterior is > 10 below the ensemble median
        mlp = lps[burn:].mean(0)
        good = mlp > np.median(mlp) - 10
        dropped = int((~good).sum())
        post = post[:, good]
        if dropped:
            print(f"   dropped {dropped} of {walkers} walkers stuck in a mode {np.median(mlp) - mlp[~good].max():.0f}+ log units less probable")
        tau = autocorr_time(post)
        rh = rhat(post)
        flat = post.reshape(-1, post.shape[-1])
        vals = model.physical(flat)
        ess = flat.shape[0] / np.maximum(tau, 1)
        n_eval = steps * walkers
        print(f"[{kind} {label}] {n_eval} forward runs in {time.time() - t0:.0f} s; acceptance {acc:.2f}; "
              f"tau {np.round(tau, 1).tolist()}; ESS {np.round(ess).astype(int).tolist()}; R-hat {np.round(rh, 3).tolist()}")
        summ = {}
        for i, nm in enumerate(model.names):
            q = np.percentile(vals[:, i], [2.5, 16, 50, 84, 97.5])
            summ[nm] = np.round(q, 3).tolist()
            line = f"   {nm:15s} MCMC median {q[2]:7.2f}  68 % {q[1]:7.2f}..{q[3]:7.2f}  95 % {q[0]:7.2f}..{q[4]:7.2f}"
            key = label.replace(" ", "_").replace("+", "and")
            if npe is not None and key in npe:
                s = npe[key][:, i]
                qn = np.percentile(s, [2.5, 16, 50, 84, 97.5])
                width = (q[3] - q[1]) / 2
                line += (f" | NPE median {qn[2]:7.2f}  68 % {qn[1]:7.2f}..{qn[3]:7.2f}"
                         f" | shift {(qn[2] - q[2]) / max(width, 1e-9):+.2f} sd, width ratio NPE/MCMC {(qn[3] - qn[1]) / max(q[3] - q[1], 1e-9):.2f}")
            print(line)
        extra = {}
        if kind == "ellipsoid":
            el = vals[:, 3] > 1.5
            extra = dict(p_long=round(float(el.mean()), 3),
                         strike_long=np.round(np.percentile(vals[el, 4], [16, 50, 84]), 1).tolist() if el.any() else None)
            print(f"   P(longer than 3 km) MCMC {extra['p_long']:.0%}"
                  + (f", NPE {np.mean(npe[key][:, 3] > 1.5):.0%}" if npe is not None and key in npe else ""))
        results[label] = dict(vals=vals, summary=summ, tau=tau.tolist(), rhat=rh.tolist(), ess=ess.tolist(),
                              acceptance=acc, n_forward=n_eval, dropped_walkers=dropped, seconds=time.time() - t0, best_lp=float(lps.max()), **extra)
    tag = f"mcmc_{kind}" + ("" if noise_kind == "gauss" else "_screens")
    rng2 = np.random.default_rng(0)
    old = dict(np.load(ML / f"{tag}.npz")) if keep_prev and (ML / f"{tag}.npz").exists() else {}
    old.update({k.replace(" ", "_").replace("+", "and"):
                v["vals"][rng2.choice(len(v["vals"]), min(20000, len(v["vals"])), replace=False)] for k, v in results.items()})
    np.savez(ML / f"{tag}.npz", **old)
    meta = {**keep_prev, **{k: {kk: vv for kk, vv in v.items() if kk != "vals"} for k, v in results.items()}}
    (ML / f"{tag}.json").write_text(json.dumps(meta, indent=1), encoding="utf-8")
    return results


if __name__ == "__main__":
    kinds = ["mogi", "ellipsoid"] if sys.argv[1] == "both" else [sys.argv[1]]
    noise_kind = arg("--noise", "gauss")
    ONLY = {"gnss": "NISAR + GNSS", "nisar": "NISAR only"}.get(arg("--data", "both"))
    for k in kinds:
        steps = arg("--steps", 4000 if k == "mogi" else 2500)
        walkers = arg("--walkers", 64 if k == "mogi" else 48)
        run(k, noise_kind, steps, walkers)
