"""Step 0 of the synthetic-data plan: how big is the atmosphere in ONE pair, at which scales?

For every pair of the chain: residual = pair - (rate from the OTHER pairs) x pair span.
The residual is atmosphere + noise (+ any non-linear motion), measured without peeking.
Outputs per-pair stats, structure functions (noise vs distance), DEM-correlated part,
and a figure. These numbers parameterise the GPU noise generator.

    uv run python scripts/atmo_stats.py                       # Mexico, descending chain
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from hachure import dem as demmod, gridding, nisar

ROOT = Path(__file__).resolve().parents[1]
BBOX = (-99.32, 19.17, -98.82, 19.58)
REF_BBOX = (-99.14, 19.53, -99.09, 19.57)
LAGS_KM = np.array([0.25, 0.5, 1, 2, 3, 5, 7, 10, 15, 20])


def structure_fn(r, dx_km, rng, n=200_000):
    """sqrt(0.5 * mean((r(p) - r(p+h))^2)) for each lag, random directions, valid pixels only."""
    h, w = r.shape
    out = []
    for lag in LAGS_KM:
        L = lag / dx_km
        a = rng.uniform(0, 2 * np.pi, n)
        y0 = rng.integers(0, h, n)
        x0 = rng.integers(0, w, n)
        y1 = np.round(y0 + L * np.sin(a)).astype(int)
        x1 = np.round(x0 + L * np.cos(a)).astype(int)
        ok = (y1 >= 0) & (y1 < h) & (x1 >= 0) & (x1 < w)
        d = r[y0[ok], x0[ok]] - r[y1[ok], x1[ok]]
        d = d[np.isfinite(d)]
        out.append(np.sqrt(0.5 * np.mean(d ** 2)) if d.size > 500 else np.nan)
    return np.array(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--direction", default="DESCENDING")
    ap.add_argument("--name", default="mexico-city")
    ap.add_argument("--height", type=float, default=2240.0)
    args = ap.parse_args()

    found = nisar.search((BBOX[0] + BBOX[2]) / 2, (BBOX[1] + BBOX[3]) / 2, direction=args.direction)
    chain = nisar.sequential_chain(nisar.dominant_track(found))
    cache = ROOT / "data" / "raw" / "nisar_cache" / ("bbox_" + "_".join(f"{v:.3f}" for v in BBOX)
                                                     + f"_h{int(args.height)}")
    pairs, names = [], []
    for g in chain:
        p = cache / f"{g.name}.npz"
        if p.exists():
            pairs.append(nisar.process(nisar.RawWindow.load(p), ref_bbox=REF_BBOX, coh_min=0.3))
            names.append(g)
    print(f"{len(pairs)} cached pairs of {len(chain)} in the chain")
    disps, cohs, x, y = nisar.align(pairs)
    days = np.array([p.days for p in pairs])
    dx_km = abs(x[1] - x[0]) / 1000

    # DEM on the pair grid, for the height-correlated (stratified) part of the atmosphere
    raw_dem = ROOT / "data" / "raw" / "copdem_mexico-city_50m.tif"
    zd, Td, crsd = demmod.load(raw_dem)
    t_pair = gridding.transform_from_centres(x, y)
    z = gridding.to_grid(zd, Td, crsd, t_pair, f"EPSG:{pairs[0].epsg}", disps.shape[1:])

    ok = np.isfinite(disps)
    tot = np.where(ok, disps, 0).sum(0)
    span = np.where(ok, days[:, None, None], 0).sum(0)
    rng = np.random.default_rng(0)
    rows, sfs, sfs_flat, examples = [], [], [], []
    for i, g in enumerate(names):
        okk = ok[i] & (span - days[i] > 200)       # LOO rate needs most of the year
        with np.errstate(all="ignore"):
            rate_loo = (tot - np.where(ok[i], disps[i], 0)) / (span - days[i])   # m/day
        sig = rate_loo * days[i]
        res = np.where(okk, disps[i] - sig, np.nan)
        res -= np.nanmedian(res)
        # height-correlated part: res ~ a + k * z  (robust enough with 10^5 pixels)
        m = np.isfinite(res) & np.isfinite(z)
        A = np.c_[np.ones(m.sum()), z[m] - 2240]
        coef, *_ = np.linalg.lstsq(A, res[m], rcond=None)
        strat = np.full(res.shape, np.nan)
        strat[m] = A @ coef
        var_strat = 1 - np.nanvar(res[m] - strat[m]) / np.nanvar(res[m])
        res_flat = res - np.where(np.isfinite(strat), strat, 0)
        # ... and a plane (orbit / long-wave), after removing stratification
        yy, xx = np.nonzero(np.isfinite(res_flat))
        P = np.c_[np.ones(len(yy)), xx, yy]
        pc, *_ = np.linalg.lstsq(P, res_flat[yy, xx], rcond=None)
        res_flat[yy, xx] -= P @ pc
        bowl = np.isfinite(sig) & (rate_loo * 365.25 < -0.15) & okk     # sinking > 15 cm/yr
        rows.append(dict(
            pair=f"{g.ref_date:%Y-%m-%d}->{g.sec_date:%Y-%m-%d}", days=float(days[i]),
            valid_pct=float(np.isfinite(res).mean() * 100),
            signal_bowl_cm=float(np.nanmedian(sig[bowl]) * 100) if bowl.any() else None,
            noise_std_cm=float(np.nanstd(res) * 100),
            strat_cm_per_100m=float(coef[1] * 100 * 100), strat_var_pct=float(var_strat * 100),
            noise_after_strat_plane_cm=float(np.nanstd(res_flat) * 100)))
        sfs.append(structure_fn(res, dx_km, rng))
        sfs_flat.append(structure_fn(res_flat, dx_km, rng))
        examples.append((res, sig))
        r = rows[-1]
        print(f"  {r['pair']} ({r['days']:.0f}d) signal in bowl {r['signal_bowl_cm']:+.2f} cm | "
              f"noise std {r['noise_std_cm']:.2f} cm, height-linked {r['strat_var_pct']:.0f}% "
              f"({r['strat_cm_per_100m']:+.2f} cm/100 m), after strat+plane {r['noise_after_strat_plane_cm']:.2f} cm")

    sfs, sfs_flat = np.array(sfs) * 100, np.array(sfs_flat) * 100
    # power law fit of the median structure function between 1 and 10 km: SF ~ h^(beta/2)
    med = np.nanmedian(sfs_flat, 0)
    sel = (LAGS_KM >= 1) & (LAGS_KM <= 10)
    slope = np.polyfit(np.log(LAGS_KM[sel]), np.log(med[sel]), 1)[0]
    summary = dict(
        n_pairs=len(rows), pixel_km=dx_km,
        median_signal_bowl_cm=float(np.nanmedian([r["signal_bowl_cm"] for r in rows])),
        median_noise_std_cm=float(np.median([r["noise_std_cm"] for r in rows])),
        median_noise_after_strat_plane_cm=float(np.median([r["noise_after_strat_plane_cm"] for r in rows])),
        median_strat_var_pct=float(np.median([r["strat_var_pct"] for r in rows])),
        lags_km=LAGS_KM.tolist(), sf_median_cm=np.nanmedian(sfs, 0).tolist(),
        sf_flat_median_cm=med.tolist(), sf_flat_p90_cm=np.nanpercentile(sfs_flat, 90, 0).tolist(),
        sf_loglog_slope_1_10km=float(slope))
    out = ROOT / "data" / "processed" / args.name
    json.dump(dict(summary=summary, pairs=rows), open(out / "atmo_stats.json", "w"), indent=1)
    print(json.dumps(summary, indent=1))
    figure(out / "atmo_stats.png", rows, sfs, sfs_flat, examples, summary)
    print(out / "atmo_stats.png")


def figure(path, rows, sfs, sfs_flat, examples, s):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(15, 9))
    ax = fig.add_subplot(2, 3, 1)
    for sf in sfs_flat:
        ax.plot(LAGS_KM, sf, color="0.7", lw=0.8)
    ax.plot(LAGS_KM, np.nanmedian(sfs, 0), "k--", lw=1.5, label="median, raw residual")
    ax.plot(LAGS_KM, np.nanmedian(sfs_flat, 0), "k", lw=2, label="median, minus height part + plane")
    ax.axhline(abs(s["median_signal_bowl_cm"]), color="#a8322a", lw=2,
               label=f"12-day sinking in the bowl ({abs(s['median_signal_bowl_cm']):.1f} cm)")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("distance between two points, km"); ax.set_ylabel("typical difference, cm")
    ax.set_title("Atmosphere noise vs distance (one pair)"); ax.legend(fontsize=8)
    ax = fig.add_subplot(2, 3, 2)
    d = [r["pair"][:10] for r in rows]
    ax.bar(range(len(rows)), [r["noise_std_cm"] for r in rows], color="0.6", label="noise std")
    ax.bar(range(len(rows)), [r["noise_after_strat_plane_cm"] for r in rows], color="0.3",
           label="after height part + plane")
    ax.plot(range(len(rows)), [abs(r["signal_bowl_cm"] or 0) for r in rows], "o", color="#a8322a",
            label="signal in bowl")
    ax.set_xticks(range(0, len(rows), 3), d[::3], rotation=60, fontsize=7)
    ax.set_ylabel("cm"); ax.set_title("Per pair: noise vs signal"); ax.legend(fontsize=8)
    ax = fig.add_subplot(2, 3, 3)
    ax.bar(range(len(rows)), [r["strat_var_pct"] for r in rows], color="#5f594d")
    ax.set_xticks(range(0, len(rows), 3), d[::3], rotation=60, fontsize=7)
    ax.set_ylabel("% of noise variance"); ax.set_title("Share explained by terrain height")
    order = np.argsort([r["noise_std_cm"] for r in rows])
    for k, i in enumerate([order[0], order[len(order) // 2], order[-1]]):
        res, sig = examples[i]
        ax = fig.add_subplot(2, 3, 4 + k)
        v = np.nanpercentile(np.abs(res), 98) * 100
        ax.imshow(res * 100, cmap="RdBu", vmin=-v, vmax=v, interpolation="nearest")
        ax.contour(np.nan_to_num(sig * 100), levels=[-0.8], colors="k", linewidths=0.8)
        ax.set_title(f"{rows[i]['pair']}: residual ±{v:.1f} cm\n(black line = bowl, 0.8 cm in this pair)",
                     fontsize=9)
        ax.axis("off")
    fig.tight_layout()
    fig.savefig(path, dpi=90)
    plt.close(fig)


if __name__ == "__main__":
    main()
