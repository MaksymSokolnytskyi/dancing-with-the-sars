"""Derived fields: generalized surface, downhill direction, Lehmann slope tone, Dufour shading.

Everything lives in pixel space of a square metric grid: x = column (east),
y = row (south). `res` is metres per pixel.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.ndimage import distance_transform_edt, gaussian_filter


@dataclass
class Field:
    z: np.ndarray        # generalized surface used for contours/stopping
    ux: np.ndarray       # unit downhill direction, column component
    uy: np.ndarray       # unit downhill direction, row component
    slope: np.ndarray    # slope angle, degrees (from the tone surface)
    dark: np.ndarray     # 0..1 ink darkness (Lehmann x Dufour light)
    res: float


def prepare(z: np.ndarray, res: float, *, sigma_dir_m: float = 90.0,
            sigma_tone_m: float = 45.0, light_azimuth: float = 315.0,
            light_altitude: float = 45.0, lehmann_max_deg: float = 45.0,
            lit_floor: float = 0.06, gamma: float = 0.7,
            z_factor: float = 1.0, valid: np.ndarray | None = None,
            fill: float | None = None, extra_shadow: np.ndarray | None = None) -> Field:
    """Build the fields hachures need.

    Two generalization scales (Dufour's trick): stroke *direction* follows the
    smoother macro-surface, stroke *tone* keeps more micro-relief.
    z_factor scales the field (use it to exaggerate small SAR displacements).
    valid: boolean mask; outside it no stroke is seeded or continued
           (e.g. InSAR coherence < threshold -> white paper = "radar unsure").
    """
    nan = np.isnan(z)
    if valid is not None:
        nan |= ~valid
    if nan.any():
        if fill is not None:
            z = np.where(nan, fill, z)
        else:
            # nearest valid value: no artificial cliffs along mask edges
            idx = distance_transform_edt(nan, return_distances=False, return_indices=True)
            z = z[tuple(idx)]
    z = z * z_factor
    z_dir = gaussian_filter(z, sigma_dir_m / res)
    z_tone = gaussian_filter(z, sigma_tone_m / res)

    gy, gx = np.gradient(z_dir, res)            # d/drow, d/dcol (m/m)
    mag = np.hypot(gx, gy)
    safe = np.where(mag > 1e-9, mag, 1.0)
    ux, uy = -gx / safe, -gy / safe              # downhill unit vector
    ux[mag <= 1e-9] = 0.0
    uy[mag <= 1e-9] = 0.0

    ty, tx = np.gradient(z_tone, res)
    slope = np.degrees(np.arctan(np.hypot(tx, ty)))

    # Lambert shading. Map axes: east = +col, north = -row.
    dzdx, dzdy = tx, -ty
    norm = np.sqrt(1 + dzdx**2 + dzdy**2)
    az, alt = np.radians(light_azimuth), np.radians(light_altitude)
    lx, ly, lz = np.sin(az) * np.cos(alt), np.cos(az) * np.cos(alt), np.sin(alt)
    illum = np.clip((-dzdx * lx - dzdy * ly + lz) / norm, 0, 1)

    lehmann = np.clip(slope / lehmann_max_deg, 0, 1)
    shadow = np.clip(1 - illum, 0, 1)
    # normalize shadow so a flat plain sits at 0 and a slope facing away from light at 1
    flat_shadow = 1 - np.sin(alt)
    shadow = np.clip((shadow - flat_shadow * 0.5) / (1 - flat_shadow * 0.5), 0, 1)
    if extra_shadow is not None:                 # e.g. cast_shadow(): a rim shading the floor
        shadow = np.maximum(shadow, extra_shadow)
    dark = (lehmann * (lit_floor + (1 - lit_floor) * shadow)) ** gamma
    if valid is not None:
        slope = np.where(valid, slope, 0.0)
        dark = np.where(valid, dark, 0.0)
    return Field(z=z_dir, ux=ux, uy=uy, slope=slope, dark=np.clip(dark, 0, 1), res=res)


def cast_shadow(z: np.ndarray, res: float, *, azimuth: float = 315.0, altitude: float = 30.0,
                max_dist_m: float | None = None, soft_m: float = 0.0) -> np.ndarray:
    """Cast shadows by ray marching toward the sun: 0 = lit, 1 = fully shaded.

    A pixel is shaded when something along the ray toward the sun rises above the ray.
    `soft_m` > 0 gives a soft edge: shade grows from 0 to 1 as the blocker pokes
    0..soft_m metres above the ray (a cheap penumbra that also hides small noise bumps).
    z must be in metres (apply any exaggeration first); north = -row, east = +col.
    """
    from scipy.ndimage import map_coordinates

    h, w = z.shape
    rel = float(np.nanmax(z) - np.nanmin(z))
    tan = np.tan(np.radians(altitude))
    if max_dist_m is None:
        max_dist_m = rel / max(tan, 1e-6)
    n = int(np.ceil(max_dist_m / res))
    az = np.radians(azimuth)
    dc, dr = np.sin(az), -np.cos(az)             # one pixel step toward the sun
    rr, cc = np.mgrid[0:h, 0:w].astype("float64")
    best = np.full(z.shape, -np.inf)
    for k in range(1, n + 1):
        zk = map_coordinates(z, [rr + k * dr, cc + k * dc], order=1, mode="nearest")
        np.maximum(best, zk - z - k * res * tan, out=best)
    if soft_m <= 0:
        return (best > 0).astype("float64")
    return np.clip(best / soft_m, 0, 1)


def prepare_scalar(d: np.ndarray, res: float, *, band: float, target_deg: float = 32.0,
                   valid: np.ndarray | None = None, sigma_dir_m: float = 300.0,
                   sigma_tone_m: float = 150.0, lit_floor: float = 0.45,
                   tone_pct: float = 85.0, **kw) -> tuple[Field, float]:
    """Hachures for a non-elevation field (displacement, velocity, backscatter change).

    The field is rescaled so a high-percentile gradient (`tone_pct`) reads as a `target_deg` slope,
    which lets the DEM machinery (Lehmann tone, shading, stop thresholds) work unchanged.
    `band` is the contour band in the field's own units (e.g. 0.02 m of subsidence);
    returns (Field, interval_in_scaled_units) to pass to trace.hachures(interval_m=...).
    Light is kept but with a high lit_floor: shading gives the bowl volume, yet lit
    sides must not lose their information.
    """
    ok = ~np.isnan(d) if valid is None else (valid & ~np.isnan(d))
    fill = float(np.nanmedian(d[ok]))
    tmp = gaussian_filter(np.where(ok, d, fill), sigma_tone_m / res)
    gy, gx = np.gradient(tmp, res)
    g = np.hypot(gx, gy)[ok]
    # tone_pct < 100 lets the strongest gradients saturate (full ink) while
    # keeping contrast in the interior of the signal
    g = g[g > np.percentile(g, 50)]  # ignore the flat/noisy half when picking the scale
    scale = np.tan(np.radians(target_deg)) / max(float(np.percentile(g, tone_pct)), 1e-12)
    f = prepare(d, res, sigma_dir_m=sigma_dir_m, sigma_tone_m=sigma_tone_m,
                lit_floor=lit_floor, z_factor=scale, valid=valid, **kw)
    return f, band * scale
