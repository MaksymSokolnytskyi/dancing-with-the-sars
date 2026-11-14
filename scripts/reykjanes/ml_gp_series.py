"""Gaussian processes in time and space for the Reykjanes NISAR vertical series.

Model for one cell:   y(t) = c + b t + f(t) + a_t(x) + e
  * c, b       offset and linear rate (Gaussian priors, wide);
  * f(t)       smooth non-linear motion, GP with squared-exponential kernel (sigma_f, l_t);
  * a_t(x)     the atmospheric / ionospheric screen of date t: uncorrelated between dates,
               smooth in space -> estimated per date with a spatial GP (sigma_s, l_s) fitted to
               the residuals of the temporal model, then removed;
  * e          white noise (sigma_n).
Hyperparameters are shared by all cells (temporal) or all dates (spatial) and found by
maximising the marginal likelihood. Two rounds: temporal -> spatial -> temporal again.

Output: output/reykjanes/ml/gp_series.npz and figures; prints the independent GNSS check
(rates) before and after.

Usage: python scripts/reykjanes/ml_gp_series.py
"""
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "output" / "reykjanes"
ML = OUT / "ml"
RNG = np.random.default_rng(1)


# ------------------------------------------------------------------ temporal GP
def k_time(t, sf, lt):
    d = t[:, None] - t[None]
    return sf ** 2 * np.exp(-0.5 * (d / lt) ** 2)


def temporal_cov(t, sf, lt, sn, sb=50.0, sc=10.0):
    """Covariance of y over the dates: offset + linear trend + smooth GP + white noise."""
    return sc ** 2 + sb ** 2 * np.outer(t, t) + k_time(t, sf, lt) + sn ** 2 * np.eye(t.size)


def temporal_nll(logp, t, Y):
    sf, lt, sn = np.exp(logp)
    C = temporal_cov(t, sf, lt, sn)
    L = np.linalg.cholesky(C)
    a = np.linalg.solve(L, Y.T)
    return 0.5 * np.sum(a ** 2) / Y.shape[0] + np.sum(np.log(np.diag(L))) + 0.5 * t.size * np.log(2 * np.pi)


def temporal_fit(t, Y):
    sub = Y[RNG.choice(Y.shape[0], min(3000, Y.shape[0]), replace=False)]
    best = None
    for x0 in ([1.0, 0.2, 1.5], [3.0, 0.1, 2.0], [0.5, 0.4, 1.0]):
        r = minimize(temporal_nll, np.log(x0), args=(t, sub), method="L-BFGS-B",
                     bounds=[(np.log(0.05), np.log(30)), (np.log(0.03), np.log(2.0)), (np.log(0.1), np.log(10))])
        if best is None or r.fun < best.fun:
            best = r
    return np.exp(best.x)


def temporal_posterior(t, Y, sf, lt, sn, sb=50.0, sc=10.0):
    """Posterior of the rate b, of the smooth curve (c + b t + f) and the residual, per row of Y.
    Rows may contain NaN: grouped by missing-data pattern."""
    n = t.size
    rate = np.full(Y.shape[0], np.nan)
    rsd = np.full(Y.shape[0], np.nan)
    smooth = np.full(Y.shape, np.nan)
    pats = {}
    for i, row in enumerate(np.isfinite(Y)):
        if row.sum() >= 0.7 * n:
            pats.setdefault(row.tobytes(), []).append(i)
    Cs = sc ** 2 + sb ** 2 * np.outer(t, t) + k_time(t, sf, lt)   # signal covariance
    for key, rows in pats.items():
        m = np.frombuffer(key, dtype=bool)
        C = Cs[np.ix_(m, m)] + sn ** 2 * np.eye(m.sum())
        Ci = np.linalg.inv(C)
        y = Y[np.ix_(rows, np.where(m)[0])]
        # rate: Cov(b, y) = sb^2 t
        kb = sb ** 2 * t[m]
        rate[rows] = y @ (Ci @ kb)
        rsd[rows] = np.sqrt(sb ** 2 - kb @ Ci @ kb)
        smooth[np.ix_(rows, np.arange(n))] = y @ (Ci @ Cs[m, :])
    return rate, rsd, smooth


# ------------------------------------------------------------------ spatial GP
def k_space(D, ss, ls):
    return ss ** 2 * np.exp(-D / ls)  # exponential (Matern 1/2): atmosphere is rough


def spatial_nll(logp, D, R):
    ss, ls, se = np.exp(logp)
    C = k_space(D, ss, ls) + se ** 2 * np.eye(D.shape[0])
    L = np.linalg.cholesky(C)
    a = np.linalg.solve(L, R.T)
    return 0.5 * np.sum(a ** 2) / R.shape[0] + np.sum(np.log(np.diag(L)))


def gp_clean(Y0, t, Xk, Yk, ny, nx, verbose=True, protect=None):
    """Temporal GP -> spatial GP screens per date -> temporal GP again. Y0: (cells, T) cm.
    protect: boolean (ny, nx) mask of cells NOT used to learn the screens (the dome), so a real
    local signal there cannot be mistaken for atmosphere; the screen is interpolated into it."""
    T = t.size
    # round 1: temporal GP
    full = np.isfinite(Y0).all(1)
    sf, lt, sn = temporal_fit(t, Y0[full])
    if verbose: print(f"temporal GP round 1: sigma_f {sf:.2f} cm, length {lt * 365.25:.0f} days, white noise {sn:.2f} cm")
    r0, rs0, sm0 = temporal_posterior(t, Y0, sf, lt, sn)
    resid = Y0 - sm0                                            # date-specific part

    # spatial GP per date on a ~1.2 km subsample, shared hyperparameters
    sub = np.zeros((ny, nx), bool)
    sub[::5, ::5] = True
    keep = sub.ravel() & full if protect is None else sub.ravel() & full & ~protect.ravel()
    sidx = np.where(keep)[0]
    P = np.c_[Xk.ravel()[sidx], Yk.ravel()[sidx]]
    D = np.hypot(P[:, None, 0] - P[None, :, 0], P[:, None, 1] - P[None, :, 1])
    Rs = resid[sidx].T[1:]                                      # dates x points (date 0 is the zero)
    best = None
    for x0 in ([1.5, 5.0, 1.0], [2.5, 10.0, 1.0], [1.0, 2.0, 1.5]):
        r = minimize(spatial_nll, np.log(x0), args=(D, Rs), method="L-BFGS-B",
                     bounds=[(np.log(0.1), np.log(20)), (np.log(10.0), np.log(80)), (np.log(0.05), np.log(10))])  # screens are long-wavelength
        if best is None or r.fun < best.fun:
            best = r
    ss, ls, se = np.exp(best.x)
    if verbose: print(f"spatial GP: screen sigma {ss:.2f} cm, length {ls:.1f} km, local noise {se:.2f} cm "
          f"({len(sidx)} points per date)")
    C = k_space(D, ss, ls) + se ** 2 * np.eye(len(sidx))
    alpha = np.linalg.solve(C, Rs.T)                            # points x dates
    allpts = np.c_[Xk.ravel(), Yk.ravel()]
    screens = np.zeros((T, ny * nx))
    for i0 in range(0, ny * nx, 4000):                           # predict in blocks (memory)
        blk = allpts[i0:i0 + 4000]
        Ks = k_space(np.hypot(blk[:, None, 0] - P[None, :, 0], blk[:, None, 1] - P[None, :, 1]), ss, ls)
        screens[1:, i0:i0 + 4000] = (Ks @ alpha).T
    clean = Y0 - screens.T
    # round 2: temporal GP on the cleaned series
    full2 = np.isfinite(clean).all(1)
    sf2, lt2, sn2 = temporal_fit(t, clean[full2])
    if verbose: print(f"temporal GP round 2: sigma_f {sf2:.2f} cm, length {lt2 * 365.25:.0f} days, white noise {sn2:.2f} cm")
    rate, rsd, smooth = temporal_posterior(t, clean, sf2, lt2, sn2)

    return clean, smooth, rate, rsd, screens, (sf, lt, sn, ss, ls, se, sf2, lt2, sn2), r0, rs0


def main():
    ML.mkdir(parents=True, exist_ok=True)
    d = np.load(OUT / "breathing_data.npz")
    up = d["up"].astype(float)                       # (T, ny, nx) cm, NYLA-referenced
    dates = [datetime.fromisoformat(str(x)) for x in d["dates"]]
    t = np.array([(x - dates[0]).days / 365.25 for x in dates])
    T, ny, nx = up.shape
    lon, lat = d["lon"].astype(float), d["lat"].astype(float)
    kx, ky = 111.32 * np.cos(np.deg2rad(63.88)), 111.2
    Xk, Yk = (lon - lon.mean()) * kx, (lat - lat.mean()) * ky   # km
    Y0 = up.reshape(T, -1).T                                    # (cells, T)

    model = d["model_up"].astype(float)
    pi, pj = np.unravel_index(np.nanargmax(model), model.shape)
    yy, xx = np.mgrid[:ny, :nx]
    protect = np.hypot(yy - pi, xx - pj) * 0.24 < 8.0           # dome: screens learned outside it
    clean, smooth, rate, rsd, screens, (sf, lt, sn, ss, ls, se, sf2, lt2, sn2), r0, rs0 = gp_clean(Y0, t, Xk, Yk, ny, nx, protect=protect)

    # independent GNSS check before/after (rates; GNSS never used in processing)
    import json
    vec = json.loads((OUT / "breathing_vectors.json").read_text(encoding="utf-8"))
    from reykjanes_gnss import STATIONS
    print("\nindependent GNSS check, vertical rate cm/yr (NISAR 3x3 cells around the station):")
    rows = []
    for s in STATIONS:
        g = vec["gnss"][s]
        i = int(np.argmin((lat - g["lat"]) ** 2 + ((lon - g["lon"]) * 0.44) ** 2))
        ii, jj = np.unravel_index(i, (ny, nx))
        cells = [(a * nx + b) for a in range(ii - 1, ii + 2) for b in range(jj - 1, jj + 2)]
        def stat(series):
            y = np.nanmedian(series[cells], 0)
            r_, s_, _ = temporal_posterior(t, y[None], sf2, lt2, sn2)
            return r_[0], s_[0], y
        ra, sa, ya = stat(Y0)
        rb, sb_, yb = stat(clean)
        gv = g["val"]
        za = (ra - gv["gnss"]) / np.hypot(sa, gv["gnss_sigma"])
        zb = (rb - gv["gnss"]) / np.hypot(sb_, gv["gnss_sigma"])
        # per-date rms against GNSS
        gt = [datetime.fromisoformat(x) for x in g["t"]]
        gd = np.array([np.median([u for q, u in zip(gt, g["up"]) if abs((q - x).days) <= 3] or [np.nan]) for x in dates])
        gd -= gd[0]
        rms_a = np.sqrt(np.nanmean((ya - gd) ** 2))
        rms_b = np.sqrt(np.nanmean((yb - gd) ** 2))
        rows.append((s, ra, sa, rb, sb_, gv["gnss"], za, zb, rms_a, rms_b, ya, yb, gd))
        print(f"  {s}: raw {ra:+.1f} ± {sa:.1f} ({za:+.1f}σ) -> GP-cleaned {rb:+.1f} ± {sb_:.1f} ({zb:+.1f}σ) "
              f"vs GNSS {gv['gnss']:+.1f}; per-date rms {rms_a:.1f} -> {rms_b:.1f} cm")

    np.savez_compressed(ML / "gp_series.npz", clean=clean.T.reshape(T, ny, nx).astype(np.float32),
                        smooth=smooth.T.reshape(T, ny, nx).astype(np.float32),
                        screens=screens.reshape(T, ny, nx).astype(np.float32),
                        rate=rate.reshape(ny, nx).astype(np.float32), rate_sigma=rsd.reshape(ny, nx).astype(np.float32),
                        rate_raw=r0.reshape(ny, nx).astype(np.float32), rate_sigma_raw=rs0.reshape(ny, nx).astype(np.float32),
                        t=t, dates=d["dates"], hyper=np.array([sf, lt, sn, ss, ls, se, sf2, lt2, sn2]))
    plot(rows, dates, rate.reshape(ny, nx), rsd.reshape(ny, nx), screens.reshape(T, ny, nx), d)
    print("saved", ML / "gp_series.npz")


def plot(rows, dates, rate, rsd, screens, d):
    fig, axes = plt.subplots(2, 4, figsize=(17, 7.5))
    for ax, (s, ra, sa, rb, sb_, gr, za, zb, rmsa, rmsb, ya, yb, gd) in zip(axes[0], rows):
        ax.plot(dates, gd, "k-", lw=1.2, label="GPS (not used)")
        ax.plot(dates, ya, "o-", ms=3, color="0.6", lw=1, label=f"NISAR raw ({rmsa:.1f} cm rms)")
        ax.plot(dates, yb, "o-", ms=3, color="tab:red", lw=1.4, label=f"NISAR GP-cleaned ({rmsb:.1f})")
        ax.set_title(f"{s}: rate {ra:+.1f}→{rb:+.1f} vs GPS {gr:+.1f} cm/yr", fontsize=9)
        ax.legend(fontsize=7, frameon=False)
        ax.grid(alpha=.3)
        ax.tick_params(axis="x", labelsize=7, rotation=30)
    a = axes[1]
    im = a[0].imshow(rate, cmap="RdBu_r", vmin=-25, vmax=25); a[0].set_title("GP rate, cm/yr", fontsize=9); fig.colorbar(im, ax=a[0])
    im = a[1].imshow(rsd, cmap="magma_r", vmin=0, vmax=4); a[1].set_title("GP rate sigma, cm/yr", fontsize=9); fig.colorbar(im, ax=a[1])
    k = int(np.nanargmax([np.nanstd(sc) for sc in screens]))
    im = a[2].imshow(screens[k], cmap="PuOr", vmin=-6, vmax=6); a[2].set_title(f"spatial-GP screen, worst date {d['dates'][k]}", fontsize=9); fig.colorbar(im, ax=a[2])
    a[3].bar(range(len(dates)), [np.nanstd(sc) for sc in screens], color="tab:purple")
    a[3].set_title("screen strength per date (std, cm)", fontsize=9)
    a[3].set_xticks(range(0, len(dates), 4)); a[3].set_xticklabels([str(x)[:10] for x in d["dates"][::4]], rotation=40, fontsize=7)
    for x in a[:3]:
        x.set_xticks([]), x.set_yticks([])
    fig.suptitle("Gaussian processes: temporal trend + smooth motion per cell, spatial atmosphere screen per date", fontsize=11)
    fig.tight_layout()
    fig.savefig(ML / "gp_series.png", dpi=110, facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    main()
