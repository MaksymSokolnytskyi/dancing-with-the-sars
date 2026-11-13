"""H9 context data: ESA WorldCover (cropland / built-up fraction on the CV rate grid) + Open-Meteo daily rain.

    uv run python scripts/fetch_h9_context.py
"""
from __future__ import annotations

import json
import urllib.request
from pathlib import Path

import numpy as np
import rasterio
from rasterio.enums import Resampling
from rasterio.warp import reproject

ROOT = Path(__file__).resolve().parents[1]
UA = {"User-Agent": "space-apps-research"}
TILES = ["N33W120", "N36W120", "N33W123", "N36W123"]
URL = "https://esa-worldcover.s3.eu-central-1.amazonaws.com/v200/2021/map/ESA_WorldCover_10m_2021_v200_{}_Map.tif"


def main():
    out = ROOT / "data" / "raw" / "context"
    out.mkdir(parents=True, exist_ok=True)
    with rasterio.open(ROOT / "data/processed/central-valley-asc/nisar_rate.tif") as ds:
        tr, crs, shape = ds.transform, ds.crs, ds.shape
    frac = {}
    for cls, name in ((40, "cropland"), (50, "builtup"), (30, "grass"), (10, "tree")):
        acc = np.zeros(shape, "float32")
        for t in TILES:
            with rasterio.Env(GDAL_HTTP_USERAGENT="space-apps-research", GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR"):
                with rasterio.open("/vsicurl/" + URL.format(t)) as src:
                    # read at ~80 m via an overview level, binary mask of the class, averaged onto the rate grid
                    ov = src.overviews(1)
                    f = next((o for o in ov if o >= 8), ov[-1] if ov else 1)
                    h, w = src.height // f, src.width // f
                    a = src.read(1, out_shape=(h, w), resampling=Resampling.nearest)
                    t2 = src.transform * src.transform.scale(src.width / w, src.height / h)
                    dst = np.full(shape, np.nan, "float32")
                    reproject((a == cls).astype("float32"), dst, src_transform=t2, src_crs=src.crs, dst_transform=tr,
                              dst_crs=crs, resampling=Resampling.average, src_nodata=None, dst_nodata=np.nan)
                    acc = np.where(np.isfinite(dst), np.maximum(acc, dst), acc) if False else np.where(np.isfinite(dst), dst, acc)
        frac[name] = acc
        print(name, f"mean fraction {np.nanmean(acc):.2f}", flush=True)
    np.savez_compressed(out / "cv_worldcover_frac.npz", **frac)
    q = ("https://archive-api.open-meteo.com/v1/archive?latitude=36.15&longitude=-119.6&start_date=2025-10-15"
         "&end_date=2026-10-25&daily=precipitation_sum&timezone=UTC")
    with urllib.request.urlopen(urllib.request.Request(q, headers=UA), timeout=120) as r:
        rain = json.loads(r.read())
    json.dump(rain, open(out / "cv_rain_openmeteo.json", "w"))
    p = np.array(rain["daily"]["precipitation_sum"], float)
    print(f"rain: {len(p)} days, total {np.nansum(p):.0f} mm, wet days {np.sum(p >= 1)}")


if __name__ == "__main__":
    main()
