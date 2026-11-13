"""P-trust: per-pixel error of the yearly rate from the pixel's own residual series, plus calibration
against the independent morning-evening difference (pre-registered).

    uv run python scripts/trust_map.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import rasterio
from rasterio.warp import transform as warp_transform

import sys
sys.path[:0] = [str(Path(__file__).parent), str(Path(__file__).parents[1] / "scripts")]
from h6_horizontal import smooth  # noqa: E402

from hachure import gridding, nisar

ROOT = Path(__file__).resolve().parents[1]
BBOX = (-99.32, 19.17, -98.82, 19.58)
REF_BBOX = (-99.14, 19.53, -99.09, 19.57)
CACHE = ROOT / "data" / "raw" / "nisar_cache" / ("bbox_" + "_".join(f"{v:.3f}" for v in BBOX) + "_h2240")


def sigma_map(direction):
    found = nisar.search((BBOX[0] + BBOX[2]) / 2, (BBOX[1] + BBOX[3]) / 2, direction=direction)
    chain = nisar.sequential_chain(nisar.dominant_track(found))
    pairs = [nisar.process(nisar.RawWindow.load(CACHE / f"{g.name}.npz"), ref_bbox=REF_BBOX, coh_min=0.3)
             for g in chain]
    v, coh, frac, x, y, _ = nisar.stack_rate(pairs)
    disps, _, x2, y2 = nisar.align(pairs)
    assert disps.shape[1:] == v.shape
    days = np.array([p.days for p in pairs])
    T = days.sum() / 365.25
    ok = np.isfinite(disps)
    tot = np.where(ok, disps, 0).sum(0)
    span = np.where(ok, days[:, None, None], 0).sum(0)
    r = []
    for i in range(len(pairs)):
        with np.errstate(all="ignore"):
            rl = (tot - np.where(ok[i], disps[i], 0)) / (span - days[i])
        r.append(np.where(ok[i], disps[i] - rl * days[i], np.nan) * 100)
    r = np.array(r)
    q = np.nansum(r ** 2, 0) + 2 * np.nansum(r[:-1] * r[1:], 0)
    q[~np.isfinite(v)] = np.nan
    qs = smooth(q, 1000 / abs(x[1] - x[0]) / 2.355)
    sig = np.sqrt(np.clip(qs, 0.01, None)) / T
    sig[~np.isfinite(v)] = np.nan
    print(f"{direction}: sigma_rate median {np.nanmedian(sig):.2f} cm/yr, 90th pct {np.nanpercentile(sig, 90):.2f}, "
          f"Q<=0 before clip {np.mean(qs[np.isfinite(qs)] <= 0) * 100:.1f}% of pixels")
    return v * 100, sig, x, y, pairs[0].epsg


def main():
    va, sa, xa, ya, ea = sigma_map("ASCENDING")
    vd, sd, xd, yd, ed = sigma_map("DESCENDING")
    src, dst = gridding.transform_from_centres(xa, ya), gridding.transform_from_centres(xd, yd)
    va_d = gridding.to_grid(va, src, f"EPSG:{ea}", dst, f"EPSG:{ed}", vd.shape)
    sa_d = gridding.to_grid(sa, src, f"EPSG:{ea}", dst, f"EPSG:{ed}", vd.shape)
    z = (va_d - vd) / np.sqrt(sa_d ** 2 + sd ** 2)
    m = np.isfinite(z)
    p1, p2 = np.mean(np.abs(z[m]) < 1), np.mean(np.abs(z[m]) < 2)
    rx, ry = warp_transform("EPSG:4326", f"EPSG:{ed}", [(REF_BBOX[0] + REF_BBOX[2]) / 2], [(REF_BBOX[1] + REF_BBOX[3]) / 2])
    XX, YY = np.meshgrid(xd, yd)
    rk = np.hypot(XX - rx[0], YY - ry[0]) / 1000
    near, far = np.nanmedian(np.abs(z[m & (rk < 5)])), np.nanmedian(np.abs(z[m & (rk >= 20)]))
    print(f"calibration on {m.sum()} common pixels: |z|<1 {p1 * 100:.0f}% (55-80), |z|<2 {p2 * 100:.0f}% (88-98); "
          f"median |z| 0-5 km {near:.2f}, 20+ km {far:.2f} (ratio {max(near, far) / min(near, far):.2f}, need < 1.5)")
    print(f"median offset of (asc - desc) {np.nanmedian((va_d - vd)[m]):+.2f} cm/yr; std of z {np.nanstd(z[m]):.2f}")
    for lo, hi in [(0, 5), (5, 10), (10, 20), (20, 40)]:
        mm = m & (rk >= lo) & (rk < hi)
        print(f"  {lo:2d}-{hi:2d} km from reference: sigma asc {np.nanmedian(sa_d[mm]):.2f}, desc {np.nanmedian(sd[mm]):.2f}, "
              f"|z|<1 {np.mean(np.abs(z[mm]) < 1) * 100:.0f}%")
    out = ROOT / "data/processed/mexico-city-asc/nisar_rate_sigma.tif"
    with rasterio.open(out, "w", driver="GTiff", height=va.shape[0], width=va.shape[1], count=2, dtype="float32",
                       crs=f"EPSG:{ea}", transform=src, nodata=np.nan, compress="deflate") as ds:
        ds.write((va / 100).astype("float32"), 1)
        ds.write((sa / 100).astype("float32"), 2)
        ds.set_band_description(1, "vertical_rate_m_per_yr_chain_sum")
        ds.set_band_description(2, "rate_1sigma_m_per_yr")
    json.dump(dict(p1=float(p1), p2=float(p2), near=float(near), far=float(far),
                   sigma_asc_median=float(np.nanmedian(sa)), sigma_desc_median=float(np.nanmedian(sd))),
              open(ROOT / "data/processed/mexico-city/trust_calibration.json", "w"), indent=1)
    np.savez_compressed(ROOT / "data/processed/mexico-city/trust_maps.npz", va=va_d, vd=vd, sa=sa_d, sd=sd, z=z)
    print(out)


if __name__ == "__main__":
    main()
