"""Read NISAR L2 GUNW layers over Mexico City (descending track 113, frame 079)
straight from the ASF archive with HTTP range requests.

By default only the window over the city is read (~5 MB of each 2.3 GB granule).

Usage:
  python scripts/reykjanes/download_gunw.py            # city window, every pair of the track/frame
  python scripts/reykjanes/download_gunw.py --full 1   # full 80 m frame for the first N pairs
"""
import json
import sys
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import h5py
import numpy as np
from pyproj import Transformer

from inspect_gunw import ROOT, open_remote

COLLECTION = "C2854335566-ASF"  # NISAR_L2_GUNW_PROVISIONAL_V1
POINT = "-99.07,19.43"          # Mexico City
TRACK_FRAME = "_113_D_079_"
WINDOW = (-99.37, 19.09, -98.74, 19.65)  # lon0, lat0, lon1, lat1 (figure box + margin)
UNW = "science/LSAR/GUNW/grids/frequencyA/unwrappedInterferogram"
RG = "science/LSAR/GUNW/metadata/radarGrid"
LAYERS = ["unwrappedPhase", "coherenceMagnitude", "connectedComponents", "ionospherePhaseScreen"]


def find_granules():
    url = ("https://cmr.earthdata.nasa.gov/search/granules.umm_json"
           f"?collection_concept_id={COLLECTION}&point={POINT}&page_size=500&sort_key=start_date")
    items = json.load(urllib.request.urlopen(url))["items"]
    out = []
    for it in items:
        name = it["umm"]["GranuleUR"]
        if TRACK_FRAME not in name:
            continue
        href = next(r["URL"] for r in it["umm"]["RelatedUrls"]
                    if r["URL"].startswith("https://") and r["URL"].endswith(".h5"))
        out.append((name, href))
    return out


def fetch(item, out_dir, full):
    name, href = item
    dst = out_dir / f"{name}.h5"
    if dst.exists():
        return f"skip {name[17:60]}"
    tmp = dst.with_suffix(".part")
    with open_remote(href) as src, h5py.File(tmp, "w") as out:
        g = src[UNW]
        x, y = g["HH/xCoordinates"][()], g["HH/yCoordinates"][()]
        epsg = int(g["HH/projection"][()])
        if full:
            sy, sx = slice(None), slice(None)
        else:
            tr = Transformer.from_crs(4326, epsg, always_xy=True)
            cx, cy = tr.transform([WINDOW[0], WINDOW[2]] * 2, [WINDOW[1]] * 2 + [WINDOW[3]] * 2)
            ix = np.where((x >= min(cx)) & (x <= max(cx)))[0]
            iy = np.where((y >= min(cy)) & (y <= max(cy)))[0]
            sy, sx = slice(iy[0], iy[-1] + 1), slice(ix[0], ix[-1] + 1)
        for k in LAYERS:
            out.create_dataset(k, data=g["HH"][k][sy, sx], compression="gzip", compression_opts=4)
        out["mask"] = g["mask"][sy, sx]
        out["x"], out["y"] = x[sx], y[sy]
        out.attrs["epsg"] = epsg
        out.attrs["centerFrequency"] = float(src["science/LSAR/GUNW/grids/frequencyA/centerFrequency"][()])
        out.attrs["granule"] = name
        # incidence angle: coarse radar-grid cube, layer closest to the basin floor (~2240 m)
        r = src[RG]
        h = r["heightAboveEllipsoid"][()]
        i = int(np.argmin(np.abs(h - 2240.0)))
        out["incidenceAngle"] = r["incidenceAngle"][i]
        out["inc_x"], out["inc_y"] = r["xCoordinates"][()], r["yCoordinates"][()]
        out.attrs["inc_height"] = float(h[i])
    tmp.rename(dst)
    return f"done {name[17:60]} {dst.stat().st_size / 1e6:.1f} MB"


if __name__ == "__main__":
    full = "--full" in sys.argv
    granules = find_granules()
    if full:
        granules = granules[:int(sys.argv[sys.argv.index("--full") + 1])]
    out_dir = ROOT / "data" / ("gunw" if full else "window")
    out_dir.mkdir(parents=True, exist_ok=True)
    print(len(granules), "granules ->", out_dir)
    with ThreadPoolExecutor(6) as ex:
        for msg in ex.map(lambda it: fetch(it, out_dir, full), granules):
            print(msg, flush=True)
