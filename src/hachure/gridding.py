"""Put rasters from other grids (NISAR UTM 80 m) onto the DEM working grid."""
from __future__ import annotations

import numpy as np
from rasterio.transform import from_origin
from rasterio.warp import Resampling, reproject


def transform_from_centres(x: np.ndarray, y: np.ndarray):
    """Affine for a north-up grid given pixel-centre coordinate vectors."""
    dx, dy = x[1] - x[0], y[0] - y[1]
    return from_origin(x[0] - dx / 2, y[0] + dy / 2, dx, dy)


def to_grid(arr: np.ndarray, src_transform, src_crs, dst_transform, dst_crs, shape,
            resampling=Resampling.bilinear) -> np.ndarray:
    out = np.full(shape, np.nan, dtype="float32")
    reproject(arr.astype("float32"), out, src_transform=src_transform, src_crs=src_crs,
              dst_transform=dst_transform, dst_crs=dst_crs, resampling=resampling,
              src_nodata=np.nan, dst_nodata=np.nan)
    return out.astype("float64")
