"""B: fetch everything the Bakota tests need (criteria in notes/experiments.md, section B; run b_geometry.py first).

- NISAR GCOV band A (20 m, HH + HV) windows for every DHDH scene of the four tracks;
- on each band A grid: WorldCover 2021 fractions, the forest-slope masks from b_geometry.py, and Sentinel-2 water
  (SCL = 6) for the two class dates and for every clear S2 date within ±5 days of a radar date (test B3).

uv run python experiments/b_fetch.py      -> data/raw/gcov_cache/..., data/processed/bakota/gridA_*.npz, s2_dates.json
Resumable: cached GCOV windows are not read again.
"""
import json
import sys
import time
from datetime import datetime
from pathlib import Path

import numpy as np
import requests

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src"), str(ROOT / "experiments")]
import k_fetch  # noqa: E402  (also points GDAL at an ASCII CA bundle)
import rasterio  # noqa: E402
from b_geometry import BBOX, GEOID, OUT, TRACKS  # noqa: E402
from rasterio.warp import Resampling  # noqa: E402

from hachure import gcov  # noqa: E402

CLASS_DATES = ("2026-08-15", "2026-03-10")
WC_URL = k_fetch.WC_URL.replace("N45E033", "N48E024")
MAX_GAP_DAYS = 5


def s2_items():
    r = requests.post(k_fetch.STAC, json={"collections": ["sentinel-2-l2a"], "bbox": list(BBOX),
                                          "datetime": "2025-10-15T00:00:00Z/2026-11-10T23:59:59Z", "limit": 300,
                                          "query": {"eo:cloud_cover": {"lt": 5}, "s2:nodata_pixel_percentage": {"lt": 5}}},
                      headers=k_fetch.UA, timeout=60)
    best = {}
    for f in r.json()["features"]:
        d = f["properties"]["datetime"][:10]
        key = (f["properties"].get("s2:nodata_pixel_percentage", 100), f["properties"]["eo:cloud_cover"])
        if d not in best or key < best[d][0]:
            best[d] = (key, f)
    return {d: v[1] for d, v in sorted(best.items())}


def s2_water(item, tf, shape, epsg):
    return k_fetch.read_cog(item["assets"]["scl"]["href"], tf, shape,
                            lambda v: np.where(v == 0, np.nan, (v == 6).astype("float32")), epsg=epsg)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    lon, lat = (BBOX[0] + BBOX[2]) / 2, (BBOX[1] + BBOX[3]) / 2
    scenes = [s for s in gcov.search(lon, lat) if s.path in TRACKS and s.pols == "DHDH"]
    print(f"{len(scenes)} GCOV scenes:", {p: sum(s.path == p for s in scenes) for p in sorted(TRACKS)}, flush=True)
    G = dict(np.load(OUT / "geometry.npz"))
    zh = float(np.nanmedian(G["z"])) + GEOID
    t0 = time.time()
    paths = gcov.fetch(scenes, BBOX, freq="A", height_m=zh, log=lambda m: print(m, flush=True))
    print(f"band A: {len(paths)} windows in {time.time() - t0:.0f}s", flush=True)

    items = s2_items()
    s2dates = [datetime.fromisoformat(d) for d in items]
    need = set(CLASS_DATES)
    for s in scenes:
        near = min(s2dates, key=lambda d: abs(d - s.date))
        if abs((near - s.date).total_seconds()) <= MAX_GAP_DAYS * 86400:
            need.add(near.date().isoformat())
    (OUT / "s2_dates.json").write_text(json.dumps({"clear": list(items), "used": sorted(need)}, indent=1), encoding="utf-8")
    print(f"S2: {len(items)} clear dates, {len(need)} used", flush=True)

    dem_tf, dem_crs = rasterio.Affine(*G["transform"]), "EPSG:32635"
    done = {}
    for p in paths:
        w = gcov.Window.load(p)
        dx, dy = w.x[1] - w.x[0], w.y[1] - w.y[0]
        atf = rasterio.Affine(dx, 0, w.x[0] - dx / 2, 0, dy, w.y[0] - dy / 2)
        key = (w.epsg, float(w.x[0]), float(w.y[0]), w.hh.shape)
        if key not in done:
            wc = lambda c: (lambda v: np.where(v == 0, np.nan, (v == c).astype("float32")))  # noqa: E731
            lay = {f"wc_{n}": k_fetch.read_cog(WC_URL, atf, w.hh.shape, wc(c), epsg=w.epsg)
                   for n, c in (("water", 80), ("crop", 40), ("tree", 10))}
            for k in [k for k in G if k.startswith(("slope_toward_", "shadow_", "layover_"))] + ["z"]:
                lay[f"geo_{k}"] = k_fetch.warp(G[k].astype("float64"), dem_tf, dem_crs, atf, w.hh.shape,
                                               how=Resampling.nearest, epsg=w.epsg)
            for d in sorted(need):
                lay[f"s2_{d}_water"] = s2_water(items[d], atf, w.hh.shape, w.epsg)
            done[key] = lay
            print(f"  class layers for grid {key[1:3]}", flush=True)
        np.savez_compressed(OUT / f"gridA_{p.stem}.npz", **done[key])
    print(f"class layers saved in {OUT}", flush=True)


if __name__ == "__main__":
    main()
