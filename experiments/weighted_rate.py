"""Experiment A: noise-weighted rate + per-pixel error, checked against GNSS.

Pair noise = a(date k+1) - a(date k), so var(pair k) = s_k^2 + s_{k+1}^2. Solve the per-DATE
noise s_k (non-negative least squares), then fit cum(t) = c + v t per pixel with weights 1/s_k^2.
Error map: the rate error from the fit, scaled by how noise grows with distance from the
reference rock (measured structure function).

    uv run python scripts/weighted_rate.py
"""
from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path

import numpy as np
import rasterio
from rasterio.warp import transform as warp_transform
from scipy.optimize import nnls

from hachure import gridding, nisar

ROOT = Path(__file__).resolve().parents[1]
BBOX = (-99.32, 19.17, -98.82, 19.58)
REF_BBOX = (-99.14, 19.53, -99.09, 19.57)
OUT = ROOT / "data" / "processed" / "mexico-city"


def load_chain():
    found = nisar.search((BBOX[0] + BBOX[2]) / 2, (BBOX[1] + BBOX[3]) / 2, direction="DESCENDING")
    chain = nisar.sequential_chain(nisar.dominant_track(found))
    cache = ROOT / "data" / "raw" / "nisar_cache" / ("bbox_" + "_".join(f"{v:.3f}" for v in BBOX) + "_h2240")
    pairs = [nisar.process(nisar.RawWindow.load(cache / f"{g.name}.npz"), ref_bbox=REF_BBOX, coh_min=0.3)
             for g in chain]
    return chain, pairs


def main():
    chain, pairs = load_chain()
    stats = json.load(open(OUT / "atmo_stats.json"))
    pv = np.array([r["noise_std_cm"] for r in stats["pairs"]]) ** 2          # cm^2 per pair
    assert len(pv) == len(pairs)
    n = len(pairs)
    # var(pair k) = s_k^2 + s_{k+1}^2  ->  NNLS for the n+1 date variances
    M = np.zeros((n, n + 1))
    for k in range(n):
        M[k, k] = M[k, k + 1] = 1
    s2_nnls, _ = nnls(M, pv)
    # The exact solve is ill-posed on a chain (it zeroes every few dates and then trusts those
    # dates completely). Stable estimator instead: a date's noise = half the mean variance of
    # the two pairs that touch it.
    pad = np.r_[pv[0], pv, pv[-1]]
    s2 = (pad[:-1] + pad[1:]) / 4
    print(f"  (exact NNLS solve zeroed {np.sum(s2_nnls < 0.05)} of {n + 1} dates - not used)")
    dates = [chain[0].ref_date] + [g.sec_date for g in chain]
    for d, v in zip(dates, s2):
        print(f"  {d:%Y-%m-%d}  date noise {np.sqrt(v):.2f} cm")

    rate_raw, coh, frac, x, y, disps = nisar.stack_rate(pairs)
    days = [p.days for p in pairs]
    cum = nisar.cumulative(disps, days, rate_raw)               # (n+1, H, W) m
    t = np.array([(d - dates[0]).days for d in dates]) / 365.25

    def fit(w):
        w = w[:, None, None] * np.isfinite(cum)
        c = np.nan_to_num(cum)
        sw, st, sy = w.sum(0), (w * t[:, None, None]).sum(0), (w * c).sum(0)
        stt, sty = (w * t[:, None, None] ** 2).sum(0), (w * t[:, None, None] * c).sum(0)
        den = sw * stt - st ** 2
        with np.errstate(all="ignore"):
            v = (sw * sty - st * sy) / den
            a = (sy - v * st) / sw
            sig = np.sqrt(sw / den)                             # m/yr if weights are 1/m^2
            # honest error: inflate by the actual scatter around the line (reduced chi^2),
            # which also absorbs seasonal errors the per-date model treats as independent
            r = c - a - v * t[:, None, None]
            chi2 = (w * r ** 2).sum(0) / np.maximum((w > 0).sum(0) - 2, 1)
            sig = sig * np.sqrt(np.maximum(chi2, 1.0))
        v[~np.isfinite(rate_raw)] = np.nan
        return v, sig

    v_ols, _ = fit(np.ones(n + 1))
    v_wls, sig_unit = fit(1 / (s2 / 1e4))                       # cm^2 -> m^2
    chain_sum = rate_raw

    # per-pixel error from each pixel's own scatter around its line (chi^2-inflated), so it
    # grows with distance from the reference by itself, as the atmosphere does
    rx, ry = warp_transform("EPSG:4326", f"EPSG:{pairs[0].epsg}",
                            [(REF_BBOX[0] + REF_BBOX[2]) / 2], [(REF_BBOX[1] + REF_BBOX[3]) / 2])
    xx, yy = np.meshgrid(x, y)
    r_km = np.hypot(xx - rx[0], yy - ry[0]) / 1000
    sigma = sig_unit
    ok_s = np.isfinite(sigma) & np.isfinite(v_wls)
    for lo, hi in [(0, 5), (5, 10), (10, 20), (20, 40)]:
        m = ok_s & (r_km >= lo) & (r_km < hi)
        if m.any():
            print(f"  error at {lo}-{hi} km from reference: median {np.median(sigma[m]) * 100:.1f} cm/yr")

    # GNSS at MMX1
    lat, lon = 19.4317, -99.0684
    sx, sy = warp_transform("EPSG:4326", f"EPSG:{pairs[0].epsg}", [lon], [lat])
    near = np.hypot(xx - sx[0], yy - sy[0]) <= 250
    g_rate = -26.7
    res = {}
    for name, v in [("chain sum (map so far)", chain_sum), ("plain fit, 28 dates", v_ols),
                    ("noise-weighted fit", v_wls)]:
        val = np.nanmedian(v[near]) * 100
        res[name] = val
        print(f"MMX1 {name:24s} {val:+6.1f} cm/yr   diff vs GNSS {val - g_rate:+5.1f}")
    sg = np.nanmedian(sigma[near]) * 100
    print(f"MMX1 1-sigma error (weighted) {sg:.1f} cm/yr, distance to reference {np.median(r_km[near]):.1f} km;"
          f" GNSS within {abs(res['noise-weighted fit'] - g_rate) / sg:.1f} sigma")
    ok = np.isfinite(v_wls)
    print(f"map: median error {np.nanmedian(sigma[ok]) * 100:.1f} cm/yr, 90th pct {np.nanpercentile(sigma[ok], 90) * 100:.1f}")
    print(f"map change weighted - chain sum: median {np.nanmedian((v_wls - chain_sum)[ok]) * 100:+.1f}, "
          f"1st/99th pct {np.nanpercentile((v_wls - chain_sum)[ok], 1) * 100:+.1f} / "
          f"{np.nanpercentile((v_wls - chain_sum)[ok], 99) * 100:+.1f} cm/yr")

    tf = gridding.transform_from_centres(x, y)
    with rasterio.open(OUT / "nisar_rate_weighted.tif", "w", driver="GTiff", height=v_wls.shape[0],
                       width=v_wls.shape[1], count=3, dtype="float32", crs=f"EPSG:{pairs[0].epsg}",
                       transform=tf, nodata=np.nan, compress="deflate") as ds:
        ds.write(v_wls.astype("float32"), 1)
        ds.write(coh.astype("float32"), 2)
        ds.write(sigma.astype("float32"), 3)
        for i, d in enumerate(["vertical_rate_m_per_yr", "mean_coherence", "rate_1sigma_m_per_yr"], 1):
            ds.set_band_description(i, d)
    json.dump(dict(mmx1=res, gnss=g_rate, mmx1_sigma=sg,
                   date_noise_cm={f"{d:%Y-%m-%d}": float(np.sqrt(v)) for d, v in zip(dates, s2)}),
              open(OUT / "weighted_rate.json", "w"), indent=1)
    figure(v_wls, chain_sum, sigma, dates, s2, res, g_rate, sg)


def figure(v, v0, sigma, dates, s2, res, g, sg):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 4, figsize=(19, 4.8))
    ax[0].bar(range(len(s2)), np.sqrt(s2), color=["#a8322a" if d.month in (5, 6, 7, 8, 9, 10) else "0.5" for d in dates])
    ax[0].set_xticks(range(0, len(dates), 3), [f"{d:%y-%m-%d}" for d in dates[::3]], rotation=60, fontsize=7)
    ax[0].set_ylabel("cm"); ax[0].set_title("Noise per date (red = May-Oct)")
    for a, im, title, cmap, lim in [(ax[1], v0 * 100, "Rate, chain sum (old)", "RdBu", (-35, 35)),
                                    (ax[2], v * 100, "Rate, noise-weighted (new)", "RdBu", (-35, 35)),
                                    (ax[3], sigma * 100, "1-sigma error, cm/yr", "magma_r", (0, 5))]:
        h = a.imshow(im, cmap=cmap, vmin=lim[0], vmax=lim[1], interpolation="nearest")
        fig.colorbar(h, ax=a, shrink=0.8)
        a.set_title(title); a.axis("off")
    fig.suptitle("MMX1 (airport): GNSS %.1f | chain sum %+.1f | plain fit %+.1f | weighted %+.1f ± %.1f cm/yr"
                 % (g, res["chain sum (map so far)"], res["plain fit, 28 dates"], res["noise-weighted fit"], sg))
    fig.tight_layout()
    fig.savefig(OUT / "weighted_rate.png", dpi=85)
    plt.close(fig)
    print(OUT / "weighted_rate.png")


if __name__ == "__main__":
    main()
