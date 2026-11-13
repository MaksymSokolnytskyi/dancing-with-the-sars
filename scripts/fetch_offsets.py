"""Read GUNW pixelOffsets (slant range, along track, correlation peak, mask) for a bbox, per pair, into a small cache.

    uv run python scripts/fetch_offsets.py --name langtang --bbox 85.37 28.17 85.67 28.47 --direction ASCENDING \
        --start 2025-10-01 --end 2026-09-30
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from rasterio.warp import transform as warp_transform

from hachure import nisar

ROOT = Path(__file__).resolve().parents[1]
G = "science/LSAR/GUNW/grids/frequencyA/pixelOffsets"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", required=True)
    ap.add_argument("--bbox", nargs=4, type=float, required=True)
    ap.add_argument("--direction", default="ASCENDING")
    ap.add_argument("--start")
    ap.add_argument("--end")
    a = ap.parse_args()
    lon, lat = (a.bbox[0] + a.bbox[2]) / 2, (a.bbox[1] + a.bbox[3]) / 2
    chain = nisar.sequential_chain(nisar.dominant_track(nisar.search(lon, lat, direction=a.direction, start=a.start, end=a.end)))
    out = ROOT / "data/raw/nisar_offsets" / a.name
    out.mkdir(parents=True, exist_ok=True)
    for g in chain:
        f = out / f"{g.name}.npz"
        if f.exists():
            continue
        with nisar.open_h5(g.url) as h5:
            x, y = h5[f"{G}/xCoordinates"][()], h5[f"{G}/yCoordinates"][()]
            epsg = nisar._epsg(h5, G)
            xs, ys = warp_transform("EPSG:4326", f"EPSG:{epsg}", [a.bbox[0], a.bbox[2]], [a.bbox[1], a.bbox[3]])
            cs, rs = nisar._slice(x, min(xs), max(xs)), nisar._slice(y, min(ys), max(ys))
            d = {k: h5[f"{G}/HH/{k}"][rs, cs] for k in ("slantRangeOffset", "alongTrackOffset", "correlationSurfacePeak")}
            d["mask"] = h5[f"{G}/mask"][rs, cs]
            units = {k: h5[f"{G}/HH/{k}"].attrs.get("units", b"?") for k in ("slantRangeOffset", "alongTrackOffset")}
            np.savez_compressed(f, x=x[cs], y=y[rs], epsg=epsg, days=g.days, **d)
        print(f"{g.ref_date:%Y-%m-%d} -> {g.sec_date:%Y-%m-%d}: {d['slantRangeOffset'].shape}, units {units}", flush=True)


if __name__ == "__main__":
    main()
