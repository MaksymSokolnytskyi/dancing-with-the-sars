"""L-retro (pre-registered in notes/experiments.md): did NISAR asc see the Chhochen Khola source slope before 26 Aug 2026?

    uv run python scripts/l_retro.py
Reads the cached GUNW windows (nisar_gunw_rate.py --name langtang-asc) and the GLO-30 DEM.
Outputs data/processed/landslide/langtang_lretro.json and a quick-look PNG.
"""
from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path

import numpy as np
import rasterio
from rasterio.warp import Resampling, reproject, transform as warp_transform
from scipy import ndimage

from hachure import gridding, nisar

ROOT = Path(__file__).resolve().parents[1]
CACHE = ROOT / "data/raw/nisar_cache/bbox_85.370_28.170_85.670_28.470_h5000"
DEM = ROOT / "data/raw/dem/langtang_cop30_utm.tif"
OUT = ROOT / "data/processed/landslide"
SUMMIT = (85.5170, 28.2563)          # Langtang Lirung, snapped to the DEM maximum within 3 km
EVENT = datetime(2026, 8, 26)
BASE = (datetime(2025, 12, 10), datetime(2026, 6, 8))


def los(raw):
    """LOS displacement (cm, + toward satellite) with the plane-only iono screen; NaN where invalid (coh < 0.3)."""
    phase = raw.phase.copy()
    if raw.iono is not None:
        ok = np.isfinite(raw.iono)
        yy, xx = np.mgrid[:raw.iono.shape[0], :raw.iono.shape[1]]
        A = np.column_stack([np.ones(ok.sum()), xx[ok], yy[ok]])
        c = np.linalg.lstsq(A, raw.iono[ok], rcond=None)[0]
        phase -= c[0] + c[1] * xx + c[2] * yy
    d = nisar.LOS_SIGN * raw.wavelength / (4 * np.pi) * phase * 100
    d[~(np.isfinite(phase) & (raw.comp > 0) & (raw.coherence >= 0.3))] = np.nan
    return d


def main():
    pairs = []
    for f in sorted(CACHE.glob("*.npz")):
        r, s = nisar._parse_dates(f.name)
        pairs.append((r, s, f))
    pairs.sort(key=lambda p: p[0])
    raws = [nisar.RawWindow.load(f) for _, _, f in pairs]
    shape = raws[0].phase.shape
    keep = [i for i, rw in enumerate(raws) if rw.phase.shape == shape]
    if len(keep) < len(raws):
        print(f"dropping {len(raws) - len(keep)} pairs on a different grid")
    pairs = [pairs[i] for i in keep]
    raws = [raws[i] for i in keep]
    x, y, epsg = raws[0].x, raws[0].y, raws[0].epsg
    print(f"{len(pairs)} pairs on grid {shape}, EPSG {epsg}")

    # DEM on the radar grid
    z = np.full(shape, np.nan, "float32")
    tr = gridding.transform_from_centres(x, y)
    with rasterio.open(DEM) as ds:
        reproject(rasterio.band(ds, 1), z, dst_transform=tr, dst_crs=f"EPSG:{epsg}", resampling=Resampling.average)
    gy, gx = np.gradient(z, y, x)                 # dz/dnorth, dz/deast (coordinate-aware)
    north_facing = -gy > 0                        # downslope has a north component
    slope = np.degrees(np.arctan(np.hypot(gx, gy)))
    XX, YY = np.meshgrid(x, y)
    sx, sy = warp_transform("EPSG:4326", f"EPSG:{epsg}", [SUMMIT[0]], [SUMMIT[1]])
    near = np.hypot(XX - sx[0], YY - sy[0]) <= 3000
    k = int(np.nanargmax(np.where(near & np.isfinite(z), z, -np.inf)))
    sx0, sy0 = XX.flat[k], YY.flat[k]
    print(f"summit snapped to DEM max {z.flat[k]:.0f} m")
    dist = np.hypot(XX - sx0, YY - sy0)

    coh = np.array([rw.coherence for rw in raws], float)
    base = [i for i, (r, s, _) in enumerate(pairs) if r >= BASE[0] and s <= BASE[1]]
    co = [i for i, (r, s, _) in enumerate(pairs) if r < EVENT < s]
    pre = [i for i, (r, s, _) in enumerate(pairs) if s < EVENT]
    assert len(co) == 1, co
    co = co[0]
    gbase = np.nanmedian(coh[base], 0)
    decor = coh[co] < 0.5 * gbase
    cand = decor & (dist <= 5000) & north_facing & (z >= 4500) & (z <= 6000)
    lab, n = ndimage.label(cand)
    OUT.mkdir(parents=True, exist_ok=True)
    res = dict(n_pairs=len(pairs), base_pairs=len(base), summit_dem_m=float(z.flat[k]))
    if n == 0 or np.bincount(lab.ravel())[1:].max() < 20:
        res["verdict"] = "could not locate source zone"
        print(res)
        json.dump(res, open(OUT / "langtang_lretro.json", "w"), indent=1)
        return
    big = np.bincount(lab.ravel())[1:].argmax() + 1
    src = lab == big
    hmed = float(np.nanmedian(z[src]))
    px = abs(x[1] - x[0])
    dsrc = ndimage.distance_transform_edt(~src) * px
    ctrl = (np.abs(z - hmed) <= 300) & north_facing & (dsrc >= 3000) & (dsrc <= 15000) & ~decor
    cy, cx = np.argwhere(src).mean(0)
    clon, clat = warp_transform(f"EPSG:{epsg}", "EPSG:4326", [x[int(cx)]], [y[int(cy)]])
    res.update(source_px=int(src.sum()), source_km2=float(src.sum() * px ** 2 / 1e6),
               source_centroid=(round(clon[0], 4), round(clat[0], 4)), source_h_median=hmed,
               source_slope_median=float(np.nanmedian(slope[src])), control_px=int(ctrl.sum()))
    print(res)
    np.savez_compressed(OUT / "langtang_zones.npz", src=src, ctrl=ctrl, x=x, y=y, z=z, north=north_facing,
                        dist=dist, decor=decor)

    rows = []
    for i, ((r, s, _), rw) in enumerate(zip(pairs, raws)):
        d = los(rw)
        comp = rw.comp.astype(int)
        cv = ctrl & np.isfinite(d)
        main = np.bincount(comp[cv].ravel()).argmax() if cv.any() else -1
        sm = src & np.isfinite(d) & (comp == main)
        cm = cv & (comp == main)
        dd = float(np.median(d[sm]) - np.median(d[cm])) if sm.sum() >= 5 and cm.sum() >= 20 else float("nan")
        role = "co" if i == co else "base" if i in base else "pre" if i in pre else "post"
        rows.append(dict(ref=f"{r:%Y-%m-%d}", sec=f"{s:%Y-%m-%d}", days=float(rw.days), d_cm=dd, n_src=int(sm.sum()),
                         n_ctrl=int(cm.sum()), coh_src=float(np.nanmedian(coh[i][src])),
                         coh_ctrl=float(np.nanmedian(coh[i][ctrl])), role=role))
        rr = rows[-1]
        print(f"  {rr['ref']} -> {rr['sec']} ({rr['days']:.0f} d) {role:4s} d {dd:+7.2f} cm  "
              f"[{rr['n_src']} src px, {rr['n_ctrl']} ctrl px]  coh src {rr['coh_src']:.2f} ctrl {rr['coh_ctrl']:.2f}")

    db = np.array([rows[i]["d_cm"] for i in base])
    db = db[np.isfinite(db)]
    sig = float(1.4826 * np.median(np.abs(db - np.median(db)))) if len(db) >= 5 else float("nan")
    last = rows[pre[-1]]
    L1 = bool(np.isfinite(last["d_cm"]) and last["d_cm"] <= -3 * sig)

    def per12(ix):
        return float(np.nanmean([abs(rows[i]["d_cm"]) * 12 / rows[i]["days"] for i in ix]))

    l2_last, l2_base = per12(pre[-3:]), per12(base)
    L2 = bool(l2_last >= 2 * l2_base)

    def ratio(i):
        return rows[i]["coh_src"] / rows[i]["coh_ctrl"]

    bmin = min(ratio(i) for i in base)
    L3 = bool(ratio(pre[-1]) < bmin)
    res.update(sigma_d=sig, n_base_valid=int(len(db)), L1=dict(d_last=last["d_cm"], thresh=-3 * sig, passed=L1),
               L2=dict(last3_per12=l2_last, base_per12=l2_base, passed=L2),
               L3=dict(ratio_last=ratio(pre[-1]), base_min=bmin,
                       base_median=float(np.median([ratio(i) for i in base])), passed=L3), rows=rows)
    print(f"sigma_d {sig:.2f} cm from {len(db)} baseline pairs")
    print(f"L1: d(last pre) {last['d_cm']:+.2f} vs threshold {-3 * sig:+.2f} -> {'PASS' if L1 else 'FAIL'}")
    print(f"L2: last-3 mean |d|/12d {l2_last:.2f} vs 2 x baseline {2 * l2_base:.2f} -> {'PASS' if L2 else 'FAIL'}")
    print(f"L3: coh ratio last {ratio(pre[-1]):.2f} vs baseline min {bmin:.2f} -> {'PASS' if L3 else 'FAIL'}")
    json.dump(res, open(OUT / "langtang_lretro.json", "w"), indent=1, default=float)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(15, 5))
    ax[0].imshow(np.ma.masked_invalid(z), cmap="gray")
    ax[0].contour(src, [0.5], colors="r")
    ax[0].contour(ctrl, [0.5], colors="c", linewidths=0.5)
    ax[0].set_title("DEM, source zone (red), control (cyan)")
    ax[1].imshow(coh[co] / np.maximum(gbase, 1e-3), vmin=0, vmax=1.5, cmap="magma")
    ax[1].set_title("co-event / baseline coherence")
    t = [datetime.fromisoformat(r["sec"]) for r in rows]
    ax[2].plot(t, [r["d_cm"] for r in rows], "o-", label="d = source - control, cm")
    ax[2].plot(t, [r["coh_src"] / r["coh_ctrl"] * 10 for r in rows], "s--", label="coh ratio x10")
    ax[2].axvline(EVENT, color="r")
    ax[2].axhline(-3 * sig, color="k", ls=":")
    ax[2].legend()
    ax[2].tick_params(axis="x", rotation=45)
    fig.tight_layout()
    fig.savefig(OUT / "langtang_lretro.png", dpi=110)


if __name__ == "__main__":
    main()
