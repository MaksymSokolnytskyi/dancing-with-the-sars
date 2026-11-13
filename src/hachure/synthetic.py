"""SYNTHETIC test fields (stand-ins until we stream real NISAR GUNW).

Mexico City: subsidence of the lacustrine clay basin, hotspots placed at the
areas known for the fastest sinking (Nezahualcoyotl/airport, Iztapalapa-Tlahuac,
Chalco, Xochimilco). Magnitudes are order-of-magnitude realistic (tens of cm/yr),
NOT measurements. Never present these numbers as data.
"""
from __future__ import annotations

import numpy as np
from rasterio.warp import transform as warp_transform
from scipy.ndimage import gaussian_filter

# (lon, lat, peak rate m/yr, sigma km)
MEXICO_HOTSPOTS = [
    (-99.03, 19.42, -0.35, 4.5),   # Nezahualcoyotl / airport
    (-99.05, 19.33, -0.30, 4.0),   # Iztapalapa - Tlahuac
    (-98.92, 19.27, -0.33, 4.0),   # Chalco
    (-99.11, 19.27, -0.18, 3.0),   # Xochimilco
    (-99.13, 19.43, -0.10, 3.5),   # historic centre
]


def _xy_grid(transform, shape):
    rows, cols = np.mgrid[0:shape[0], 0:shape[1]]
    x = transform.c + (cols + 0.5) * transform.a
    y = transform.f + (rows + 0.5) * transform.e
    return x, y


def subsidence_rate(transform, crs, shape, elevation: np.ndarray,
                    basin_floor_m: float = 2280.0, seed: int = 7):
    """Vertical rate (m/yr, negative = down) + synthetic coherence (0..1)."""
    x, y = _xy_grid(transform, shape)
    rate = np.zeros(shape)
    for lon, lat, peak, sig_km in MEXICO_HOTSPOTS:
        (hx,), (hy,) = warp_transform("EPSG:4326", crs, [lon], [lat])
        r2 = (x - hx) ** 2 + (y - hy) ** 2
        rate += peak * np.exp(-r2 / (2 * (sig_km * 1000) ** 2))
    # only compressible lake clays sink: fade out above the basin floor
    on_clay = 1 / (1 + np.exp((elevation - basin_floor_m) / 12.0))
    rate *= on_clay

    rng = np.random.default_rng(seed)
    res = abs(transform.a)
    noise = gaussian_filter(rng.standard_normal(shape), 1500 / res)   # atmosphere-like
    rate += 0.004 * noise / noise.std()

    # coherence: high on city (stable scatterers), lower on vegetated steep slopes + patches
    slope = np.hypot(*np.gradient(elevation, res))
    patches = gaussian_filter(rng.standard_normal(shape), 800 / res)
    coh = 0.85 - 1.2 * np.clip(slope, 0, 0.5) + 0.25 * patches / patches.std()
    return rate, np.clip(coh, 0, 1)
