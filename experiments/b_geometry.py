"""B (prep, radar-independent): where should NISAR see radar shadow and layover in the Dniester canyon at Bakota?

Uses only the Copernicus DEM and the viewing geometry (line-of-sight unit vector and incidence angle cubes from the GCOV
metadata, one scene per track); no backscatter value is read. Writes the predicted masks for the pre-registered tests.

uv run python experiments/b_geometry.py   -> data/processed/bakota/geometry.npz, geometry.json
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.interpolate import RegularGridInterpolator

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src"), str(ROOT / "experiments")]
import k_fetch  # noqa: E402,F401  (ASCII CA bundle for GDAL; must come before rasterio)
from hachure import dem, gcov, nisar  # noqa: E402

BBOX = (26.90, 48.55, 27.06, 48.64)          # Bakota bay and the canyon meanders around it
TRACKS = {100, 172, 94, 166}
OUT = ROOT / "data/processed/bakota"
GEOID = 30.0                                 # approximate EGM2008 geoid height here, m (DEM heights -> ellipsoid)


def cube_layer(h5, name, height, x, y):
    """A radarGrid cube at the layer nearest `height` (m above the ellipsoid), bilinear onto the DEM grid (x, y)."""
    g = "science/LSAR/GCOV/metadata/radarGrid"
    hz = h5[f"{g}/heightAboveEllipsoid"][()]
    k = int(np.argmin(np.abs(hz - height)))
    cx, cy, lay = h5[f"{g}/xCoordinates"][()], h5[f"{g}/yCoordinates"][()], h5[f"{g}/{name}"][k].astype("float64")
    if cy[0] > cy[-1]:
        cy, lay = cy[::-1], lay[::-1]
    f = RegularGridInterpolator((cy, cx), lay, bounds_error=False, fill_value=None)
    yy, xx = np.meshgrid(y, x, indexing="ij")
    return f(np.column_stack([yy.ravel(), xx.ravel()])).reshape(yy.shape)


def shadow_layover(z, res, ux, uy, inc_deg):
    """Per-pixel radar shadow (terrain toward the sensor rises above the ray) and layover (slope toward the sensor
    steeper than the incidence angle). (ux, uy) = horizontal unit vector from the target toward the sensor."""
    H, W = z.shape
    tan_elev = np.tan(np.radians(90.0 - inc_deg))          # ray elevation above the horizon
    gy, gx = np.gradient(z, res)                            # dz/dy (rows run south), dz/dx
    slope_toward = gx * ux + (-gy) * uy                     # rise per metre toward the sensor
    layover = slope_toward > np.tan(np.radians(inc_deg))
    shadow = np.zeros_like(z, bool)
    step = np.array([ux, -uy]) / np.hypot(ux, uy)          # (dx, drow-sign) in pixels per unit
    rr, cc = np.mgrid[0:H, 0:W].astype("float64")
    for d in range(1, 120):                                 # up to 3.6 km toward the sensor
        r2 = np.rint(rr + step[1] * d).astype(int)
        c2 = np.rint(cc + step[0] * d).astype(int)
        ok = (r2 >= 0) & (r2 < H) & (c2 >= 0) & (c2 < W)
        rise = np.full_like(z, -np.inf)
        rise[ok] = (z[r2[ok], c2[ok]] - z[ok]) / (d * res)
        shadow |= rise > tan_elev
    return shadow, layover, slope_toward


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    tif = dem.fetch_copdem(BBOX, ROOT / "data/raw/dem/bakota_30m.tif", res=30.0)
    z, tf, crs = dem.load(tif)
    res = tf.a
    x = tf.c + res * (np.arange(z.shape[1]) + 0.5)
    y = tf.f + tf.e * (np.arange(z.shape[0]) + 0.5)
    print("DEM", z.shape, crs, f"height {np.nanmin(z):.0f}–{np.nanmax(z):.0f} m")
    lon, lat = (BBOX[0] + BBOX[2]) / 2, (BBOX[1] + BBOX[3]) / 2
    scenes = [s for s in gcov.search(lon, lat) if s.path in TRACKS and s.pols == "DHDH"]
    pick = {}
    for s in sorted(scenes, key=lambda s: s.date):
        pick[s.path] = s                                     # latest scene per track (geometry repeats)
    out, info = {"z": z, "transform": np.array(tf)[:6]}, {"bbox": BBOX, "dem_shape": z.shape, "tracks": {}}
    for t, s in sorted(pick.items()):
        with nisar.open_h5(s.url) as h5:
            epsg = int(h5["science/LSAR/GCOV/grids/frequencyB/projection"][()])
            assert epsg == crs.to_epsg(), (epsg, crs)
            h = float(np.nanmedian(z)) + GEOID
            lx, ly = cube_layer(h5, "losUnitVectorX", h, x, y), cube_layer(h5, "losUnitVectorY", h, x, y)
            inc = cube_layer(h5, "incidenceAngle", h, x, y)
        n = np.hypot(lx, ly)
        ux, uy = float(np.median(lx / n)), float(np.median(ly / n))
        theta = float(np.median(inc))
        sh, lo, st = shadow_layover(z, res, ux, uy, theta)
        out[f"shadow_{t}"], out[f"layover_{t}"], out[f"slope_toward_{t}"] = sh, lo, st.astype("float32")
        info["tracks"][t] = {"scene": s.name, "dir": s.direction, "toward_sensor_EN": [round(ux, 3), round(uy, 3)],
                             "incidence_deg": round(theta, 2), "shadow_px": int(sh.sum()), "layover_px": int(lo.sum()),
                             "shadow_frac": round(float(sh.mean()), 5), "layover_frac": round(float(lo.mean()), 5)}
        print(t, info["tracks"][t])
    gy, gx = np.gradient(z, res)
    slope = np.degrees(np.arctan(np.hypot(gx, gy)))
    info["slope_deg_p50_p99_max"] = [round(float(np.nanpercentile(slope, q)), 1) for q in (50, 99, 100)]
    print("slope p50/p99/max", info["slope_deg_p50_p99_max"])
    np.savez_compressed(OUT / "geometry.npz", **out)
    (OUT / "geometry.json").write_text(json.dumps(info, indent=1), encoding="utf-8")


if __name__ == "__main__":
    main()
