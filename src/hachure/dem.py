"""Fetch Copernicus GLO-30 DEM windows (public AWS COGs) and project them to UTM."""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import rasterio
from rasterio.crs import CRS
from rasterio.merge import merge
from rasterio.warp import Resampling, calculate_default_transform, reproject

COPDEM_URL = "https://copernicus-dem-30m.s3.amazonaws.com/{name}/{name}.tif"


def _tile_name(lat: int, lon: int) -> str:
    ns = "N" if lat >= 0 else "S"
    ew = "E" if lon >= 0 else "W"
    return f"Copernicus_DSM_COG_10_{ns}{abs(lat):02d}_00_{ew}{abs(lon):03d}_00_DEM"


def tiles_for_bbox(w: float, s: float, e: float, n: float) -> list[str]:
    return [
        _tile_name(lat, lon)
        for lat in range(math.floor(s), math.ceil(n))
        for lon in range(math.floor(w), math.ceil(e))
    ]


def utm_crs(lon: float, lat: float) -> CRS:
    zone = int((lon + 180) // 6) + 1
    return CRS.from_epsg((32600 if lat >= 0 else 32700) + zone)


def _trim_nan_edges(z: np.ndarray, transform):
    """Drop edge rows/cols containing NaN (left over from the rotation into UTM)."""
    r0, r1, c0, c1 = 0, z.shape[0], 0, z.shape[1]
    while True:
        sub = np.isnan(z[r0:r1, c0:c1])
        edges = {
            "top": sub[0].sum(), "bottom": sub[-1].sum(),
            "left": sub[:, 0].sum(), "right": sub[:, -1].sum(),
        }
        worst = max(edges, key=edges.get)
        if edges[worst] == 0:
            break
        if worst == "top":
            r0 += 1
        elif worst == "bottom":
            r1 -= 1
        elif worst == "left":
            c0 += 1
        else:
            c1 -= 1
    new_t = transform * rasterio.Affine.translation(c0, r0)
    return z[r0:r1, c0:c1], new_t


def fetch_copdem(bbox: tuple[float, float, float, float], out_path: str | Path,
                 res: float = 30.0) -> Path:
    """Read bbox (lon/lat W,S,E,N) from GLO-30 COGs over HTTP, reproject to UTM at `res` m.

    Only the needed window is streamed, not whole tiles. Result is cached at out_path.
    """
    out_path = Path(out_path)
    if out_path.exists():
        return out_path
    w, s, e, n = bbox
    with rasterio.Env(GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR",
                      CPL_VSIL_CURL_ALLOWED_EXTENSIONS=".tif"):
        srcs = []
        for name in tiles_for_bbox(w, s, e, n):
            try:
                srcs.append(rasterio.open("/vsicurl/" + COPDEM_URL.format(name=name)))
            except rasterio.errors.RasterioIOError:
                print(f"  (no tile {name} - ocean?)")
        if not srcs:
            raise RuntimeError("No Copernicus DEM tiles for this bbox")
        src_crs = srcs[0].crs
        mosaic, mtrans = merge(srcs, bounds=(w, s, e, n), nodata=np.nan, dtype="float32")
        for src in srcs:
            src.close()

    h, wd = mosaic.shape[1:]
    left, top = mtrans.c, mtrans.f
    right, bottom = left + wd * mtrans.a, top + h * mtrans.e
    dst_crs = utm_crs((w + e) / 2, (s + n) / 2)
    transform, width, height = calculate_default_transform(
        src_crs, dst_crs, wd, h, left, bottom, right, top, resolution=res)
    dst = np.full((height, width), np.nan, dtype="float32")
    reproject(mosaic[0], dst, src_transform=mtrans, src_crs=src_crs,
              dst_transform=transform, dst_crs=dst_crs,
              resampling=Resampling.bilinear, src_nodata=np.nan, dst_nodata=np.nan)
    dst, transform = _trim_nan_edges(dst, transform)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with rasterio.open(out_path, "w", driver="GTiff", height=dst.shape[0], width=dst.shape[1],
                       count=1, dtype="float32", crs=dst_crs, transform=transform,
                       nodata=np.nan, compress="deflate") as dst_ds:
        dst_ds.write(dst, 1)
    return out_path


def load(path: str | Path):
    with rasterio.open(path) as ds:
        return ds.read(1).astype("float64"), ds.transform, ds.crs
