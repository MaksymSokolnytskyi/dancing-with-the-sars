"""GOES-19 ABI band 13 (10.3 um) brightness temperature over a box, at given UTC times.

Public NOAA bucket on AWS (no account). Reads only the window from each 4 MB CONUS file
(h5py over HTTP), falls back to the Full Disk sector if CONUS does not cover the box.

    uv run python scripts/fetch_goes.py                 # all NISAR overpass times, Mexico City
    uv run python scripts/fetch_goes.py 10 15           # other bands (7.3 um water vapour, 12.3 um split window)
"""
from __future__ import annotations

import re
import sys
from datetime import datetime, timedelta
from pathlib import Path

import fsspec
import h5py
import numpy as np
import requests

ROOT = Path(__file__).resolve().parents[1]
BUCKET = "https://noaa-goes19.s3.amazonaws.com"
BOX = (-99.32, 19.17, -98.82, 19.58)
OUT = ROOT / "data" / "raw" / "goes19"


def latlon(x, y, p):
    """GOES fixed grid (scan angles, rad) -> lat/lon, GOES-R PUG formulas."""
    req, rpol, H, lon0 = p["semi_major_axis"], p["semi_minor_axis"], \
        p["perspective_point_height"] + p["semi_major_axis"], np.radians(p["longitude_of_projection_origin"])
    X, Y = np.meshgrid(x, y)
    a = np.sin(X) ** 2 + np.cos(X) ** 2 * (np.cos(Y) ** 2 + (req / rpol) ** 2 * np.sin(Y) ** 2)
    b = -2 * H * np.cos(X) * np.cos(Y)
    c = H ** 2 - req ** 2
    with np.errstate(invalid="ignore"):
        rs = (-b - np.sqrt(b ** 2 - 4 * a * c)) / (2 * a)
    sx, sy, sz = rs * np.cos(X) * np.cos(Y), -rs * np.sin(X), rs * np.cos(X) * np.sin(Y)
    lat = np.degrees(np.arctan((req / rpol) ** 2 * sz / np.sqrt((H - sx) ** 2 + sy ** 2)))
    lon = np.degrees(lon0 - np.arctan(sy / (H - sx)))
    return lat, lon


def find(product, t, band=13):
    doy = t.timetuple().tm_yday
    keys = []
    for hh in {t.hour, (t - timedelta(minutes=10)).hour}:
        r = requests.get(BUCKET, params={"list-type": 2,
                                         "prefix": f"{product}/{t.year}/{doy:03d}/{hh:02d}/OR_{product.replace('ABI-', 'ABI-')}-M6C{band:02d}"},
                         timeout=60)
        keys += re.findall(r"<Key>([^<]+)</Key>", r.text)
    best, bd = None, 1e9
    for k in keys:
        s = re.search(r"_s(\d{4})(\d{3})(\d{2})(\d{2})(\d{2})", k)
        st = datetime(int(s[1]), 1, 1) + timedelta(days=int(s[2]) - 1, hours=int(s[3]), minutes=int(s[4]), seconds=int(s[5]))
        d = abs((st - t).total_seconds())
        if d < bd:
            best, bd = k, d
    return best, bd


def read_box(key):
    with fsspec.open(f"{BUCKET}/{key}", "rb", block_size=2 ** 20) as f, h5py.File(f, "r") as h:
        p = {k: float(np.atleast_1d(v)[0]) for k, v in h["goes_imager_projection"].attrs.items()
             if np.issubdtype(np.asarray(v).dtype, np.number)}
        xs = h["x"][()] * h["x"].attrs["scale_factor"] + h["x"].attrs["add_offset"]
        ys = h["y"][()] * h["y"].attrs["scale_factor"] + h["y"].attrs["add_offset"]
        # coarse lat/lon on every 20th pixel to find the window, then exact on the window
        lat, lon = latlon(xs[::20], ys[::20], p)
        m = (lat > BOX[1] - 0.2) & (lat < BOX[3] + 0.2) & (lon > BOX[0] - 0.2) & (lon < BOX[2] + 0.2)
        if not m.any():
            return None
        r, c = np.nonzero(m)
        r0, r1 = max(r.min() * 20 - 20, 0), min(r.max() * 20 + 40, len(ys))
        c0, c1 = max(c.min() * 20 - 20, 0), min(c.max() * 20 + 40, len(xs))
        rad = h["Rad"][r0:r1, c0:c1].astype(float)
        sf, ao = h["Rad"].attrs["scale_factor"][0], h["Rad"].attrs["add_offset"][0]
        rad = rad * sf + ao
        fk1, fk2, bc1, bc2 = (float(h[k][()]) for k in ("planck_fk1", "planck_fk2", "planck_bc1", "planck_bc2"))
        with np.errstate(all="ignore"):
            bt = (fk2 / np.log(fk1 / rad + 1) - bc1) / bc2
        lat, lon = latlon(xs[c0:c1], ys[r0:r1], p)
    return bt, lat, lon


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    for band in [int(b) for b in sys.argv[1:]] or [13]:
        fetch_band(band)


def fetch_band(band):
    times = []
    for name in ("mexico-city", "mexico-city-asc"):
        z = np.load(ROOT / "data" / "processed" / name / "nisar_series.npz")
        times += [datetime.fromisoformat(str(d)) for d in z["dates"]]
    for t in sorted(set(times)):
        out = OUT / f"bt{band:02d}_{t:%Y%m%dT%H%M}.npz"
        if out.exists():
            continue
        got = None
        for product in ("ABI-L1b-RadC", "ABI-L1b-RadF"):
            key, dt = find(product, t, band)
            if key is None or dt > 900:
                continue
            got = read_box(key)
            if got is not None:
                break
        if got is None:
            print(f"{t:%Y-%m-%d %H:%M} no GOES scan covering the box")
            continue
        bt, lat, lon = got
        inbox = (lat >= BOX[1]) & (lat <= BOX[3]) & (lon >= BOX[0]) & (lon <= BOX[2])
        np.savez_compressed(out, bt=bt, lat=lat, lon=lon, key=key, dt=dt)
        print(f"B{band} {t:%Y-%m-%d %H:%M} {product[-5:]} scan {dt:.0f}s off | box BT median {np.nanmedian(bt[inbox]):.0f} K, "
              f"cold<235K {np.mean(bt[inbox] < 235) * 100:.0f}%", flush=True)


if __name__ == "__main__":
    main()
