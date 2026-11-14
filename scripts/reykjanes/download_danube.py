"""Read NISAR L2 GCOV backscatter (HH, HV) over the Danube Delta straight from
the ASF archive with HTTP range requests: only the delta window is transferred,
not the 2 GB granule. 20 m posting, linear gamma0 power.

Usage: python scripts/reykjanes/download_danube.py [--list]
"""
import json
import re
import sys
import urllib.request
from concurrent.futures import ThreadPoolExecutor

import h5py
import numpy as np
from pyproj import Transformer

from inspect_gunw import ROOT, open_remote

COLLECTIONS = ["C2854338529-ASF",  # NISAR_L2_GCOV_PROVISIONAL_V1 (Jun-Oct 2026)
               "C3622241360-ASF",  # NISAR_L2_GCOV_V1
               "C3622214170-ASF"]  # NISAR_L2_GCOV_BETA_V1 (Nov 2025 - Jan 2026: the 2026 winter scenes)
POINT = "29.6,45.2"                # Sulina
WINDOW = (28.60, 44.30, 30.20, 45.60)  # lon0, lat0, lon1, lat1: whole delta + Razim-Sinoie lagoons + margin
YEAR = ""  # all NISAR acquisitions over the delta (Nov 2025 onwards)
OUT = ROOT / "data" / "danube_nisar"
GRID = "science/LSAR/GCOV/grids/frequencyA"
NAME = re.compile(r"GCOV_(\d+)_(\d+)_([AD])_(\d+)_\w+?_(\w+?)_A_(\d{8})T(\d{6})")


def find_granules():
    found = {}
    for coll in COLLECTIONS:
        url = ("https://cmr.earthdata.nasa.gov/search/granules.umm_json"
               f"?collection_concept_id={coll}&point={POINT}&page_size=500&sort_key=start_date")
        for it in json.load(urllib.request.urlopen(url))["items"]:
            name = it["umm"]["GranuleUR"]
            m = NAME.search(name)
            if not m or not m.group(6).startswith(YEAR):
                continue
            href = next(r["URL"] for r in it["umm"]["RelatedUrls"]
                        if r["URL"].startswith("https://") and r["URL"].endswith(".h5"))
            cycle, track, direction, frame, _, date, _ = m.groups()
            found[(cycle, track, frame)] = dict(name=name, href=href, cycle=int(cycle), track=int(track),
                                                direction=direction, frame=int(frame), date=date)
    return sorted(found.values(), key=lambda d: (d["track"], d["frame"], d["date"]))


def fetch(item):
    dst = OUT / f"{item['name']}.h5"
    tag = f"{item['date']} track {item['track']}{item['direction']} frame {item['frame']}"
    if dst.exists():
        return f"skip {tag}"
    tmp = dst.with_suffix(".part")
    with open_remote(item["href"]) as src:
        g = src[GRID]
        epsg = int(g["projection"][()])
        x, y = g["xCoordinates"][()], g["yCoordinates"][()]
        tr = Transformer.from_crs(4326, epsg, always_xy=True)
        cx, cy = tr.transform([WINDOW[0], WINDOW[2]] * 2, [WINDOW[1]] * 2 + [WINDOW[3]] * 2)
        ix = np.where((x >= min(cx)) & (x <= max(cx)))[0]
        iy = np.where((y >= min(cy)) & (y <= max(cy)))[0]
        if ix.size == 0 or iy.size == 0:
            return f"no overlap {tag}"
        sy, sx = slice(iy[0], iy[-1] + 1), slice(ix[0], ix[-1] + 1)
        with h5py.File(tmp, "w") as out:
            for k in ("HHHH", "HVHV"):
                out.create_dataset(k, data=g[k][sy, sx], compression="gzip", compression_opts=4, shuffle=True,
                                   chunks=(256, 256))
            out.create_dataset("mask", data=g["mask"][sy, sx], compression="gzip")
            out["x"], out["y"] = x[sx], y[sy]
            out.attrs.update(epsg=epsg, granule=item["name"], date=item["date"], track=item["track"],
                             frame=item["frame"], direction=item["direction"], cycle=item["cycle"],
                             units="linear gamma0 power (RTC)", window_lonlat=WINDOW)
            hh = out["HHHH"][()]
        full = ((max(cx) - min(cx)) / 20) * ((max(cy) - min(cy)) / 20)  # pixels if the granule covered the whole box
        cover = (np.isfinite(hh) & (hh > 0)).sum() / full
    tmp.rename(dst)
    return f"done {tag} {hh.shape} {dst.stat().st_size / 1e6:.0f} MB, valid share of box {min(cover, 1):.0%}"


if __name__ == "__main__":
    granules = find_granules()
    print(len(granules), "granules in", YEAR)
    for d in granules:
        print(" ", d["date"], f"track {d['track']}{d['direction']} frame {d['frame']} cycle {d['cycle']}")
    if "--list" in sys.argv:
        sys.exit()
    OUT.mkdir(parents=True, exist_ok=True)
    with ThreadPoolExecutor(5) as ex:
        for msg in ex.map(fetch, granules):
            print(msg, flush=True)
