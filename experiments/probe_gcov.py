"""GC: first look at NISAR GCOV over Central Valley (criteria in notes/experiments.md, section GC).

uv run python experiments/probe_gcov.py   -> data/processed/central-valley-asc/gcov_gc.json + gcov_gc.png
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
import rasterio
from rasterio.warp import Resampling, reproject

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src")]
from hachure import engine, gcov  # noqa: E402

cfg = engine.load_site(ROOT / "sites/central-valley-asc.json")
bbox = tuple(cfg["bbox"])
lon, lat = (bbox[0] + bbox[2]) / 2, (bbox[1] + bbox[3]) / 2
scenes = [s for s in gcov.search(lon, lat, direction="ASCENDING", start="2026-01-01", end="2026-10-31")]
top = max({(s.path, s.frame) for s in scenes}, key=lambda k: sum((s.path, s.frame) == k for s in scenes))
monthly = {}
for s in scenes:
    if (s.path, s.frame) == top:
        monthly.setdefault(s.date.month, s)
picked = [monthly[m] for m in sorted(monthly)]
print(f"track {top}: {len(picked)} dates", [f"{s.date:%m-%d}" for s in picked])
t0 = time.time()
paths = gcov.fetch(picked, bbox)
print(f"fetched in {time.time() - t0:.0f}s")

frac = np.load(ROOT / cfg["fields"]["landcover"]["npz"])
with rasterio.open(ROOT / cfg["fields"]["landcover"]["grid_tif"]) as ds:
    dst_tf, dst_crs, shape = ds.transform, ds.crs, (ds.height, ds.width)
crop = frac["cropland"] >= 0.6
built = frac["builtup"] >= 0.3
print("classes px:", int(crop.sum()), int(built.sum()))


def on_grid(win, term):
    v = (win.hv if term == "hv" else win.hh).astype("float64")
    bad = ~np.isfinite(v) | (v <= 0) | (win.mask == 0) | (win.mask == 255)
    v[bad] = np.nan
    dx, dy = win.x[1] - win.x[0], win.y[1] - win.y[0]
    src_tf = rasterio.Affine(dx, 0, win.x[0] - dx / 2, 0, dy, win.y[0] - dy / 2)
    out = np.full(shape, np.nan)
    reproject(v, out, src_transform=src_tf, src_crs=f"EPSG:{win.epsg}", dst_transform=dst_tf, dst_crs=dst_crs,
              resampling=Resampling.average, src_nodata=np.nan, dst_nodata=np.nan)
    return 10 * np.log10(out)


rows = []
for p in paths:
    w = gcov.Window.load(p)
    if not rows:
        u, c = np.unique(w.mask, return_counts=True)
        print("mask values:", dict(zip(u.tolist(), c.tolist())), "window", w.hh.shape, "dx", w.x[1] - w.x[0])
    hh, hv = on_grid(w, "hh"), on_grid(w, "hv")
    r = {"date": w.date.date().isoformat()}
    for cls, m in (("crop", crop), ("built", built)):
        for term, a in (("hh", hh), ("hv", hv)):
            vals = a[m & np.isfinite(a)]
            r[f"{cls}_{term}"] = round(float(np.median(vals)), 2) if vals.size > 100 else None
            r[f"{cls}_n"] = int(vals.size)
    rows.append(r)
    print(r)


def rng(key):
    v = [r[key] for r in rows if r[key] is not None]
    return round(max(v) - min(v), 2) if v else None


res = {"GC1_crop_hv_range_db": rng("crop_hv"), "GC2_built_hv_range_db": rng("built_hv"),
       "GC3_crop_hh_range_db": rng("crop_hh"), "rows": rows, "track": top}
res["GC1"] = "pass" if res["GC1_crop_hv_range_db"] is not None and res["GC1_crop_hv_range_db"] >= 3 else "fail"
res["GC2"] = "pass" if res["GC2_built_hv_range_db"] is not None and res["GC2_built_hv_range_db"] <= 1.5 else "fail"
res["GC3"] = "pass" if (res["GC1_crop_hv_range_db"] or 0) > (res["GC3_crop_hh_range_db"] or 99) else "fail"
out = ROOT / "data/processed/central-valley-asc"
(out / "gcov_gc.json").write_text(json.dumps(res, indent=1), encoding="utf-8")
print({k: v for k, v in res.items() if k != "rows"})

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
fig, ax = plt.subplots(1, 3, figsize=(15, 4.5))
d = [r["date"][5:] for r in rows]
for cls, c in (("crop", "tab:green"), ("built", "tab:gray")):
    ax[0].plot(d, [r[f"{cls}_hv"] for r in rows], "o-", c=c, label=f"{cls} HV")
    ax[0].plot(d, [r[f"{cls}_hh"] for r in rows], "s--", c=c, label=f"{cls} HH", alpha=.6)
ax[0].legend(); ax[0].set_ylabel("median gamma0, dB"); ax[0].set_title("Central Valley, GCOV monthly")
w0, w1 = gcov.Window.load(paths[0]), gcov.Window.load(paths[len(paths) // 2])
ax[1].imshow(w0.db("hv"), vmin=-28, vmax=-10, cmap="gray"); ax[1].set_title(f"HV {w0.date:%Y-%m-%d}")
ax[2].imshow(w1.db("hv"), vmin=-28, vmax=-10, cmap="gray"); ax[2].set_title(f"HV {w1.date:%Y-%m-%d}")
fig.tight_layout(); fig.savefig(out / "gcov_gc.png", dpi=90)
print("wrote", out / "gcov_gc.png")
