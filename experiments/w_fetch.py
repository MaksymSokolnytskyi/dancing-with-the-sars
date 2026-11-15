"""W: fetch everything the Wisconsin quad-pol tests need (criteria in notes/experiments.md, section W).

- NISAR GCOV band B (80 m, HH + HV) windows for every QPDH scene of the four tracks over the full window;
- band A (20 m) full covariance (10 terms) over the sub-window, one scene per track and calendar month;
- class layers: ESA WorldCover 2021 fractions and Sentinel-2 L2A mosaics (water, NDVI, NDWI) for the leaf-on and the
  leaf-off date, on one 80 m UTM grid for band B and on each band A grid;
- Open-Meteo hourly 2 m air temperature at every overpass.

uv run python experiments/w_fetch.py      -> data/raw/gcov_cache/..., data/processed/wisconsin/{grid.npz, gridA_*.npz, temps.json}
Resumable: cached GCOV windows are not read again.
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
import requests

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src"), str(ROOT / "experiments")]
import k_fetch  # noqa: E402  (also points GDAL at an ASCII CA bundle)
import rasterio  # noqa: E402
from rasterio.transform import from_origin  # noqa: E402
from rasterio.warp import transform_bounds  # noqa: E402

from hachure import gcov  # noqa: E402

BBOX = (-91.00, 45.02, -90.56, 45.28)        # full window, band B
SUB = (-90.77, 45.17, -90.67, 45.24)         # quad-pol sub-window, band A (chosen from WorldCover only)
HEIGHT = 350.0                               # terrain above the ellipsoid, m
TRACKS = {76, 148, 41, 113}
S2_DATES = ("2026-07-14", "2025-11-03")      # leaf-on, leaf-off (snow-free)
EPSG, RES = 32615, 80.0
OUT = ROOT / "data/processed/wisconsin"
WC_URL = k_fetch.WC_URL.replace("N45E033", "N45W093")


def target_grid():
    w, s, e, n = transform_bounds("EPSG:4326", f"EPSG:{EPSG}", *BBOX)
    w, n = np.floor(w / RES) * RES, np.ceil(n / RES) * RES
    W, H = int(np.ceil((e - w) / RES)), int(np.ceil((n - s) / RES))
    return from_origin(w, n, RES, RES), (H, W)


def s2_mosaic(date, dst_tf, shape, epsg=EPSG):
    """Sentinel-2 layers of one date; scenes with less no-data first, later scenes only fill the gaps."""
    r = requests.post(k_fetch.STAC, json={"collections": ["sentinel-2-l2a"], "bbox": list(BBOX),
                                          "datetime": f"{date}T00:00:00Z/{date}T23:59:59Z", "limit": 20},
                      headers=k_fetch.UA, timeout=60)
    items = sorted(r.json()["features"], key=lambda f: (f["properties"].get("s2:nodata_pixel_percentage", 100),
                                                        f["properties"]["eo:cloud_cover"]))
    out = None
    for it in items:
        a = it["assets"]
        valid = lambda v: np.where(v == 0, np.nan, v.astype("float32"))  # noqa: E731
        lay = {"water": k_fetch.read_cog(a["scl"]["href"], dst_tf, shape,
                                         lambda v: np.where(v == 0, np.nan, (v == 6).astype("float32")), epsg=epsg)}
        for b in ("red", "nir", "green"):
            lay[b] = k_fetch.read_cog(a[b]["href"], dst_tf, shape, valid, epsg=epsg)
        lay["ndvi"] = (lay["nir"] - lay["red"]) / (lay["nir"] + lay["red"])
        lay["ndwi"] = (lay["green"] - lay["nir"]) / (lay["green"] + lay["nir"])
        lay = {k: lay[k] for k in ("water", "ndvi", "ndwi")}
        print(f"  S2 {date}: {it['id']} filled {np.isfinite(lay['ndvi']).mean():.0%}", flush=True)
        if out is None:
            out = lay
        else:
            for k in out:
                out[k] = np.where(np.isfinite(out[k]), out[k], lay[k])
        if np.isfinite(out["ndvi"]).all():
            break
    return out


def layers(dst_tf, shape, epsg=EPSG):
    wc = lambda c: (lambda v: np.where(v == 0, np.nan, (v == c).astype("float32")))  # noqa: E731
    out = {f"wc_{n}": k_fetch.read_cog(WC_URL, dst_tf, shape, wc(c), epsg=epsg)
           for n, c in (("water", 80), ("crop", 40), ("tree", 10))}
    for d in S2_DATES:
        for k, v in s2_mosaic(d, dst_tf, shape, epsg).items():
            out[f"s2_{d}_{k}"] = v
    return out


def monthly(scenes):
    first = {}
    for s in sorted(scenes, key=lambda s: s.date):
        first.setdefault((s.path, s.date.year, s.date.month), s)
    return sorted(first.values(), key=lambda s: (s.path, s.date))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    lon, lat = (BBOX[0] + BBOX[2]) / 2, (BBOX[1] + BBOX[3]) / 2
    scenes = [s for s in gcov.search(lon, lat) if s.path in TRACKS and "_QPDH_" in s.name]
    print(f"{len(scenes)} QPDH scenes:", {p: sum(s.path == p for s in scenes) for p in sorted(TRACKS)}, flush=True)
    t0 = time.time()
    paths = gcov.fetch(scenes, BBOX, freq="B", height_m=HEIGHT, log=lambda m: print(m, flush=True))
    print(f"band B: {len(paths)} windows in {time.time() - t0:.0f}s", flush=True)
    sub = [s for s in gcov.search((SUB[0] + SUB[2]) / 2, (SUB[1] + SUB[3]) / 2)
           if s.path in TRACKS and "_QPDH_" in s.name]
    a_scenes = monthly(sub)
    t0 = time.time()
    a_paths = gcov.fetch(a_scenes, SUB, freq="A", height_m=HEIGHT, offdiag=True, log=lambda m: print(m, flush=True))
    print(f"band A: {len(a_paths)} windows in {time.time() - t0:.0f}s", flush=True)

    (OUT / "temps.json").write_text(json.dumps(k_fetch.temps(scenes + a_scenes, BBOX), indent=1, default=str),
                                    encoding="utf-8")
    print("temperatures saved", flush=True)
    tf, shape = target_grid()
    np.savez_compressed(OUT / "grid.npz", transform=np.array(tf)[:6], epsg=EPSG, **layers(tf, shape))
    done = {}                                     # band A windows on the same grid share their class layers
    for p in a_paths:
        w = gcov.Window.load(p)                   # on the product's own UTM zone (32616 here), not EPSG above
        dx, dy = w.x[1] - w.x[0], w.y[1] - w.y[0]
        atf = rasterio.Affine(dx, 0, w.x[0] - dx / 2, 0, dy, w.y[0] - dy / 2)
        key = (w.epsg, float(w.x[0]), float(w.y[0]), w.hh.shape)
        if key not in done:
            done[key] = layers(atf, w.hh.shape, w.epsg)
        np.savez_compressed(OUT / f"gridA_{p.stem}.npz", **done[key])
    print(f"class layers saved in {OUT}", flush=True)


if __name__ == "__main__":
    main()
