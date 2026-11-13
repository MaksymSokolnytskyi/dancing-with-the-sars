"""H9 (pre-registered): is the small-scale part of the CV morning iono screens tied to cropland and rain?

    uv run python scripts/h9_iono_landcover.py
"""
from __future__ import annotations

import json
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import rasterio

from hachure import gridding, nisar

ROOT = Path(__file__).resolve().parents[1]
BBOX = (-120.20, 35.80, -119.00, 36.50)
CACHE = ROOT / "data/raw/nisar_cache" / ("bbox_" + "_".join(f"{v:.3f}" for v in BBOX) + "_h80")
B = 25   # 2 km blocks of 80 m pixels


def blocks(a):
    H, W = (a.shape[0] // B) * B, (a.shape[1] // B) * B
    with np.errstate(all="ignore"):
        return np.nanmean(a[:H, :W].reshape(H // B, B, W // B, B), (1, 3))


def main():
    found = nisar.search((BBOX[0] + BBOX[2]) / 2, (BBOX[1] + BBOX[3]) / 2, direction="ASCENDING")
    chain = nisar.sequential_chain(nisar.dominant_track(found))
    with rasterio.open(ROOT / "data/processed/central-valley-asc/nisar_rate.tif") as ds:
        tr, crs, shape = ds.transform, ds.crs, ds.shape
    lc = np.load(ROOT / "data/raw/context/cv_worldcover_frac.npz")
    crop, built = blocks(lc["cropland"]), blocks(lc["builtup"])
    rain = json.load(open(ROOT / "data/raw/context/cv_rain_openmeteo.json"))["daily"]
    rmap = {d: (p or 0.0) for d, p in zip(rain["time"], rain["precipitation_sum"])}

    def wet(d):
        return sum(rmap.get((d - timedelta(days=i)).strftime("%Y-%m-%d"), 0) for i in range(0, 4)) >= 2

    win, summ, stds = None, None, []
    for g in chain:
        raw = nisar.RawWindow.load(CACHE / f"{g.name}.npz")
        s = raw.iono.copy()
        ok = np.isfinite(s)
        yy, xx = np.mgrid[:s.shape[0], :s.shape[1]]
        c = np.linalg.lstsq(np.column_stack([np.ones(ok.sum()), xx[ok], yy[ok]]), s[ok], rcond=None)[0]
        s = s - (c[0] + c[1] * xx + c[2] * yy)
        b = -nisar.LOS_SIGN * raw.wavelength / (4 * np.pi) * (-s) / np.cos(np.radians(raw.incidence)) * 100   # cm
        b = gridding.to_grid(b, gridding.transform_from_centres(raw.x, raw.y), f"EPSG:{raw.epsg}", tr, crs, shape)
        bb = blocks(b)
        r0, r1 = nisar._parse_dates(g.name)
        m = r1.month
        if m in (11, 12, 1, 2, 3):
            win = bb if win is None else win + bb
        if m in (5, 6, 7, 8, 9):
            summ = bb if summ is None else summ + bb
        stds.append((np.nanstd(bb), wet(r0) or wet(r1), r0.date(), r1.date()))
    ok = np.isfinite(win) & np.isfinite(summ) & np.isfinite(crop)
    rw = float(np.corrcoef(win[ok], crop[ok])[0, 1])
    rs = float(np.corrcoef(summ[ok], crop[ok])[0, 1])
    cm, bm = ok & (crop > 0.6), ok & (built > 0.3)
    d = float(np.mean(win[cm]) - np.mean(win[bm])) if bm.sum() >= 3 else float("nan")
    sw = np.median([s for s, w, *_ in stds if w]) if any(w for _, w, *_ in stds) else np.nan
    sd = np.median([s for s, w, *_ in stds if not w])
    print(f"blocks {ok.sum()} (cropland>0.6: {cm.sum()}, built-up>0.3: {bm.sum()})")
    print(f"(a) corr(B_win, cropland) {rw:+.2f} (>= +0.30), corr(B_sum, cropland) {rs:+.2f} (|.| <= {abs(rw) / 2:.2f}) -> "
          f"{'PASS' if rw >= 0.3 and abs(rs) <= abs(rw) / 2 else 'fail'}")
    print(f"(b) mean B_win cropland - built-up {d:+.2f} cm (>= +3) -> {'PASS' if d >= 3 else 'fail'}")
    print(f"(c) block std of b_k: wet pairs {sw:.2f} cm (n={sum(w for _, w, *_ in stds)}), dry {sd:.2f} cm "
          f"(n={sum(not w for _, w, *_ in stds)}), ratio {sw / sd:.2f} (>= 1.5) -> {'PASS' if sw / sd >= 1.5 else 'fail'}")
    print(f"context: B_win block mean {np.nanmean(win):+.2f} cm, std {np.nanstd(win):.2f}; B_sum mean {np.nanmean(summ):+.2f}, std {np.nanstd(summ):.2f}")
    for s, w, a, e in stds:
        print(f"   {a} -> {e}  {'WET' if w else '   '}  std {s:.2f} cm")
    np.savez_compressed(ROOT / "data/processed/central-valley-asc/h9_blocks.npz", win=win, summ=summ, crop=crop, built=built)
    json.dump(dict(r_win=rw, r_sum=rs, crop_minus_built=d, std_wet=float(sw), std_dry=float(sd)),
              open(ROOT / "data/processed/central-valley-asc/h9.json", "w"), indent=1)


if __name__ == "__main__":
    main()
