"""Library of real atmospheric screens for the source inference, per track (159A, 082D),
on the same 1.2 km bins as the NPE / MCMC observation vector.

  raw   per date: cumulative LOS along the pair chain minus a per-bin quadratic trend in
        time (offset + rate + deceleration, so the slowing uplift does not leak into the
        "atmosphere"); what is left is troposphere + residual ionosphere + unwrapping /
        decorrelation noise + any non-quadratic motion. Variance corrected for the 3
        fitted terms by sqrt(N / (N - 3)).
  gp    the smooth part of each raw screen: spatial Gaussian-process posterior mean
        (exponential kernel, hyper-parameters learned on bins > 8 km from Svartsengi so
        the dome cannot be taken for atmosphere, length >= 10 km), plus the learned
        white nugget, which the simulator adds as independent noise.

The rate in the observation is (cumulative at the last date - at the first date) / span,
so one noise realisation of the rate map is (s_j - s_i) / span for two dates i != j, with
a random sign. Pairs that involve the first or last date are kept out of the training
library (that difference is the real data's own noise).

Writes output/reykjanes/ml/screens.npz (+ screens.png).
Usage: python scripts/reykjanes/ml_screens.py
"""
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize

import reykjanes_mogi as M
from ml_npe_source import BIN_KM, Obs, exp_cov_fit, lsq_fit
from reykjanes_analysis import MIN_PAIR_COH, chain, los_cm, pair_coherence, pairs_by_track
from reykjanes_sources import Data

ML = Path(__file__).resolve().parents[2] / "output" / "reykjanes" / "ml"
DOME_KM = 8.0


def bin_index(D, t):
    """Same binning as Obs: returns (inverse index of each pixel, keep mask over bins, counts)."""
    d = D.ins[t]
    key = np.c_[np.round(d["e"] / BIN_KM), np.round(d["n"] / BIN_KM)].astype(int)
    _, inv = np.unique(key, axis=0, return_inverse=True)
    inv = inv.ravel()
    cnt = np.bincount(inv)
    return inv, cnt >= 6


def bin_nanmean(v, inv, keep):
    ok = np.isfinite(v)
    s = np.bincount(inv[ok], weights=v[ok], minlength=keep.size)
    c = np.bincount(inv[ok], minlength=keep.size)
    with np.errstate(invalid="ignore"):
        return np.where(c[keep] > 0, s[keep] / np.maximum(c[keep], 1), np.nan)


def cumulative(D, t):
    """Per-date cumulative LOS (cm) on the observation pixels. A pair missing for a pixel is
    bridged with that pixel's mean rate, so a gap adds no spurious screen change."""
    geo, P = D.geo, D.P
    mask = D.obs[t]["mask"]
    sl = (slice(None, None, M.STEP), slice(None, None, M.STEP))
    pairs = []
    for t0, t1, p in chain(pairs_by_track()[t]):
        if pair_coherence(p, geo) < MIN_PAIR_COH:
            continue
        d, _ = los_cm(p, geo)
        if d is not None:
            pairs.append((t0, t1, d[sl][mask]))
    rate = D.obs[t]["y"] / 365.25                   # cm/day, the observation itself
    dates = [pairs[0][0]]
    c = [np.zeros(mask.sum())]
    for t0, t1, d in pairs:
        if t0 != dates[-1]:                          # chain broken by a dropped pair: bridge with the rate
            c.append(c[-1] + rate * (t0 - dates[-1]).days)
            dates.append(t0)
        c.append(c[-1] + np.where(np.isfinite(d), d, rate * (t1 - t0).days))
        dates.append(t1)
    return dates, np.array(c), len(pairs)


def raw_screens(tdays, Cb):
    """Cb: (T, bins). Residual of a per-bin quadratic in time, variance-corrected."""
    x = (tdays - tdays.mean()) / 365.25
    A = np.c_[np.ones_like(x), x, x ** 2]
    coef, *_ = np.linalg.lstsq(A, Cb, rcond=None)
    r = Cb - A @ coef
    return r * np.sqrt(len(x) / (len(x) - 3))


def gp_fit(e, n, S, far):
    """Exponential-kernel GP shared by all dates, learned on far bins. Returns (sigma, L, nugget)."""
    Dm = np.hypot(e[far, None] - e[None, far], n[far, None] - n[None, far])
    Y = S[:, far]
    Y = Y - Y.mean(1, keepdims=True)

    def nll(th):
        s, L, w = np.exp(th)
        K = s ** 2 * np.exp(-Dm / L) + w ** 2 * np.eye(far.sum())
        try:
            Lc = np.linalg.cholesky(K)
        except np.linalg.LinAlgError:
            return 1e12
        a = np.linalg.solve(Lc, Y.T)
        return 0.5 * (a ** 2).sum() + Y.shape[0] * np.log(np.diag(Lc)).sum()
    r = minimize(nll, np.log([1.0, 15.0, 0.3]), method="L-BFGS-B",
                 bounds=[(np.log(0.05), np.log(20)), (np.log(10.0), np.log(80)), (np.log(0.02), np.log(10))])
    return np.exp(r.x)


def gp_smooth(e, n, S, far, s, L, w):
    """Posterior mean of the smooth screen at every bin, conditioned on the far bins only."""
    Dff = np.hypot(e[far, None] - e[None, far], n[far, None] - n[None, far])
    Daf = np.hypot(e[:, None] - e[None, far], n[:, None] - n[None, far])
    K = s ** 2 * np.exp(-Dff / L) + w ** 2 * np.eye(far.sum())
    Ks = s ** 2 * np.exp(-Daf / L)
    mu = S[:, far].mean(1, keepdims=True)
    return mu + (Ks @ np.linalg.solve(K, (S[:, far] - mu).T)).T


def ledoit_wolf(X):
    """Shrunk second-moment matrix of realisations X (n, p), target mu*I (Ledoit & Wolf 2004):
    the library has fewer realisations than bins, so the raw matrix is singular."""
    n, p = X.shape
    S = X.T @ X / n
    mu = np.trace(S) / p
    d2 = ((S - mu * np.eye(p)) ** 2).sum() / p
    b2 = min(sum(((np.outer(x, x) - S) ** 2).sum() for x in X) / n ** 2 / p, d2)
    a = b2 / d2
    return a * mu * np.eye(p) + (1 - a) * S, a


def param_cov(X, e, n):
    """Exponential covariance + white nugget fitted to the library's second moments vs distance:
    a full-rank, smooth stand-in for the empirical matrix (which has fewer realisations than bins)."""
    S = X.T @ X / len(X)
    Dm = np.hypot(e[:, None] - e[None], n[:, None] - n[None])
    off = ~np.eye(e.size, dtype=bool)
    edges = np.arange(0, 26, 1.2)
    mid = edges[:-1] + 0.6
    emp = np.array([S[off & (Dm >= a) & (Dm < a + 1.2)].mean() for a in edges[:-1]])
    var0 = np.diag(S).mean()
    best = None
    for L in np.linspace(2, 150, 149):
        k = np.exp(-mid / L)
        s2 = max(float(k @ emp / (k @ k)), 1e-4)
        err = ((emp - s2 * k) ** 2).sum()
        if best is None or err < best[0]:
            best = (err, s2, L)
    _, s2, L = best
    nug = max(var0 - s2, 0.05 * var0)
    return s2 * np.exp(-Dm / L) + nug * np.eye(e.size), (np.sqrt(s2), L, np.sqrt(nug))


def pair_library(T, exclude_ends=True):
    ends = {0, T - 1}
    return np.array([(i, j) for i in range(T) for j in range(i + 1, T)
                     if not (exclude_ends and (i in ends or j in ends))])


def main():
    ML.mkdir(parents=True, exist_ok=True)
    D = Data()
    O = Obs(D)
    _, resid = lsq_fit(O)
    noise = exp_cov_fit(O, resid)
    out, fig_rows = {}, []
    for t in M.TRACKS:
        t0 = time.time()
        dates, C, npairs = cumulative(D, t)
        inv, keep = bin_index(D, t)
        Cb = np.array([bin_nanmean(c, inv, keep) for c in C])          # (T, bins)
        assert Cb.shape[1] == O.tr[t]["y"].size
        tdays = np.array([(d - dates[0]).days for d in dates], float)
        span = tdays[-1]
        e, n = O.tr[t]["e"], O.tr[t]["n"]
        S = raw_screens(tdays, Cb)
        far = np.hypot(e, n) > DOME_KM
        s, L, w = gp_fit(e, n, S, far)
        G = gp_smooth(e, n, S, far, s, L, w)
        lib = pair_library(len(dates))
        diff_raw = (S[lib[:, 1]] - S[lib[:, 0]]) * 365.25 / span        # cm/yr, rate-noise realisations
        diff_gp = (G[lib[:, 1]] - G[lib[:, 0]]) * 365.25 / span
        # with sign flips the library mean is zero; covariance of +/- realisations = second moment
        cov_raw, a_raw = ledoit_wolf(diff_raw)
        cov_gp, a_gp = ledoit_wolf(diff_gp)
        cov_gp += (2 * w ** 2 * (365.25 / span) ** 2) * np.eye(e.size)
        cov_par, (ps, pl, pn) = param_cov(diff_raw, e, n)
        real = (S[-1] - S[0]) * 365.25 / span
        exp_sig = noise[t]["sigma"]
        print(f"{t}: {npairs} pairs, {len(dates)} dates over {span:.0f} days, {e.size} bins ({time.time() - t0:.0f} s)")
        print(f"   raw screen std per date: median {np.median(S.std(1)):.2f} cm, worst {S.std(1).max():.2f} cm on {dates[int(S.std(1).argmax())]:%Y-%m-%d}")
        print(f"   GP: sigma {s:.2f} cm, length {L:.1f} km, nugget {w:.2f} cm; smooth part carries "
              f"{np.mean(G.var(1) / np.maximum(S.var(1), 1e-9)):.0%} of the screen variance")
        print(f"   rate-noise std (cm/yr): raw library {np.sqrt(np.diag(cov_raw)).mean():.2f}, GP library "
              f"{np.sqrt(np.diag(cov_gp)).mean():.2f}, Gaussian exp model {exp_sig:.2f} (L {noise[t]['L']:.1f} km); "
              f"the real data's own end-to-end screen {real.std():.2f}")
        kurt = np.mean(((diff_raw - diff_raw.mean(0)) / diff_raw.std(0)) ** 4)
        print(f"   Ledoit-Wolf shrinkage: raw {a_raw:.2f}, GP {a_gp:.2f}; parametric fit to the library: "
              f"sigma {ps:.2f} cm/yr, length {pl:.0f} km, white {pn:.2f} cm/yr")
        print(f"   kurtosis of the raw rate noise {kurt:.1f} (Gaussian = 3); {len(lib)} training pairs (ends excluded)")
        out.update({f"raw_{t}": S.astype(np.float32), f"gp_{t}": G.astype(np.float32), f"days_{t}": tdays,
                    f"span_{t}": span, f"nugget_{t}": w, f"gp_hyper_{t}": np.array([s, L, w]),
                    f"cov_raw_{t}": cov_raw, f"cov_par_{t}": cov_par, f"cov_gp_{t}": cov_gp, f"dates_{t}": np.array([str(d.date()) for d in dates]),
                    f"e_{t}": e, f"n_{t}": n, f"real_{t}": real})
        fig_rows.append((t, e, n, S, G, dates, cov_raw, noise[t], diff_raw))
    np.savez(ML / "screens.npz", **out)
    plot(fig_rows)


def plot(rows):
    fig, ax = plt.subplots(len(rows), 4, figsize=(18, 4.3 * len(rows)))
    for r, (t, e, n, S, G, dates, cov_raw, nz, diff) in enumerate(rows):
        k = int(S.std(1).argmax())
        for a, v, title in ((ax[r, 0], S[k], f"{t} raw screen {dates[k]:%Y-%m-%d}"), (ax[r, 1], G[k], f"{t} GP screen {dates[k]:%Y-%m-%d}")):
            sc = a.scatter(e, n, c=v - np.nanmean(v), s=22, cmap="PuOr", vmin=-4, vmax=4, marker="s")
            a.set_title(title, fontsize=9); a.set_aspect("equal"); fig.colorbar(sc, ax=a, label="cm")
        ax[r, 2].bar(range(len(dates)), S.std(1), color="tab:purple", label="raw")
        ax[r, 2].bar(range(len(dates)), G.std(1), color="tab:orange", alpha=.7, label="GP smooth")
        ax[r, 2].set_title(f"{t} screen strength per date (std, cm)", fontsize=9); ax[r, 2].legend(fontsize=8, frameon=False)
        Dm = np.hypot(e[:, None] - e[None], n[:, None] - n[None])
        bins = np.arange(0, 24, 1.2)
        emp = [cov_raw[(Dm >= b) & (Dm < b + 1.2)].mean() for b in bins]
        ax[r, 3].plot(bins + .6, emp, "o-", label="real screens (rate noise)")
        ax[r, 3].plot(bins + .6, nz["sigma"] ** 2 * np.exp(-(bins + .6) / nz["L"]), "--", label="Gaussian exp model (NPE)")
        ax[r, 3].set_title(f"{t} noise covariance vs distance, (cm/yr)^2", fontsize=9); ax[r, 3].set_xlabel("km")
        ax[r, 3].legend(fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(ML / "screens.png", dpi=105, facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    main()
