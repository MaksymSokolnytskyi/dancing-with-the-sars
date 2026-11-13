"""Clip ESA WorldCover 2021 (10 m, CC-BY 4.0) to a site's bbox: reads only the needed windows of the public COG tiles.

    uv run python scripts/fetch_worldcover.py sites/mexico-city-asc.json

Writes data/raw/context/<site id>_worldcover_10m.tif (EPSG:4326, uint8 classes: 10 tree, 20 shrub, 30 grass, 40 crop,
50 built-up, 60 bare, 80 water, 90 wetland, ...). Attribution: © ESA WorldCover project 2021 / Contains modified Copernicus
Sentinel data (2021) processed by ESA WorldCover consortium.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np
import rasterio
from rasterio.transform import from_origin
from rasterio.windows import from_bounds

ROOT = Path(__file__).resolve().parents[1]
URL = "https://esa-worldcover.s3.eu-central-1.amazonaws.com/v200/2021/map/ESA_WorldCover_10m_2021_v200_{}_Map.tif"
RES = 1 / 12000                                   # WorldCover pixel, degrees


def tile_name(lat0, lon0):
    return f"{'N' if lat0 >= 0 else 'S'}{abs(lat0):02d}{'E' if lon0 >= 0 else 'W'}{abs(lon0):03d}"


def main():
    cfg = json.loads(Path(sys.argv[1]).read_text(encoding="utf-8"))
    w, s, e, n = cfg["bbox"]
    W = int(round((e - w) / RES)); H = int(round((n - s) / RES))
    out = np.zeros((H, W), "uint8")
    T = from_origin(w, n, RES, RES)
    tiles = {(math.floor(la / 3) * 3, math.floor(lo / 3) * 3) for la in (s, n - 1e-9) for lo in (w, e - 1e-9)}
    for la, lo in sorted(tiles):
        name = tile_name(la, lo)
        with rasterio.Env(GDAL_HTTP_USERAGENT="space-apps-research", GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR"):
            with rasterio.open("/vsicurl/" + URL.format(name)) as src:
                bw, bs = max(w, lo), max(s, la)
                be, bn = min(e, lo + 3), min(n, la + 3)
                win = from_bounds(bw, bs, be, bn, src.transform).round_offsets().round_lengths()
                a = src.read(1, window=win)
                c0 = int(round((bw - w) / RES)); r0 = int(round((n - bn) / RES))
                h, ww = min(a.shape[0], H - r0), min(a.shape[1], W - c0)
                out[r0:r0 + h, c0:c0 + ww] = np.where(a[:h, :ww] > 0, a[:h, :ww], out[r0:r0 + h, c0:c0 + ww])
                print(f"{name}: {a.shape[1]} x {a.shape[0]} px", flush=True)
    dst = ROOT / "data/raw/context" / f"{cfg['id']}_worldcover_10m.tif"
    dst.parent.mkdir(parents=True, exist_ok=True)
    with rasterio.open(dst, "w", driver="GTiff", height=H, width=W, count=1, dtype="uint8", crs="EPSG:4326", transform=T,
                       compress="deflate") as ds:
        ds.write(out, 1)
    vals, cnt = np.unique(out, return_counts=True)
    print(dst, f"{dst.stat().st_size / 1e6:.1f} MB", {int(v): round(c / out.size, 3) for v, c in zip(vals, cnt)})


if __name__ == "__main__":
    main()
