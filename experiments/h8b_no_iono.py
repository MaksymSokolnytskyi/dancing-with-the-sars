"""H8b: rate maps with vs without the per-pair ionosphere correction (pre-registered).

    uv run python scripts/h8b_no_iono.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from rasterio.warp import transform as warp_transform

sys.path[:0] = [str(Path(__file__).parent), str(Path(__file__).parents[1] / "scripts")]
from h6_horizontal import smooth  # noqa: E402

from hachure import gridding, nisar

ROOT = Path(__file__).resolve().parents[1]
BBOX = (-99.32, 19.17, -98.82, 19.58)
REF_BBOX = (-99.14, 19.53, -99.09, 19.57)
CACHE = ROOT / "data" / "raw" / "nisar_cache" / ("bbox_" + "_".join(f"{v:.3f}" for v in BBOX) + "_h2240")
MMX1 = (-99.0684, 19.4317)


def rate(direction, iono):
    found = nisar.search((BBOX[0] + BBOX[2]) / 2, (BBOX[1] + BBOX[3]) / 2, direction=direction)
    chain = nisar.sequential_chain(nisar.dominant_track(found))
    pairs = [nisar.process(nisar.RawWindow.load(CACHE / f"{g.name}.npz"), ref_bbox=REF_BBOX, coh_min=0.3, iono=iono)
             for g in chain]
    v, coh, frac, x, y, _ = nisar.stack_rate(pairs)
    return v * 100, x, y, pairs[0].epsg


def at_mmx1(v, x, y, epsg):
    sx, sy = warp_transform("EPSG:4326", f"EPSG:{epsg}", [MMX1[0]], [MMX1[1]])
    xx, yy = np.meshgrid(x, y)
    return float(np.nanmedian(v[np.hypot(xx - sx[0], yy - sy[0]) <= 250]))


def main():
    out = {}
    for iono in (True, False):
        vd, xd, yd, ed = rate("DESCENDING", iono)
        va, xa, ya, ea = rate("ASCENDING", iono)
        dst = gridding.transform_from_centres(xd, yd)
        va_d = gridding.to_grid(va, gridding.transform_from_centres(xa, ya), f"EPSG:{ea}", dst, f"EPSG:{ed}", vd.shape)
        dx = abs(xd[1] - xd[0])
        row = dict(mmx1_desc=at_mmx1(vd, xd, yd, ed), mmx1_asc=at_mmx1(va, xa, ya, ea))
        for km in (1, 10):
            D = smooth(va_d - vd, km * 1000 / dx / 2.355)
            row[f"D_std_{km}km"] = float(np.nanstd(D - np.nanmedian(D)))
        key = "with_iono" if iono else "no_iono"
        out[key] = row
        print(f"{key:10s} MMX1 desc {row['mmx1_desc']:+.1f}  asc {row['mmx1_asc']:+.1f}  (GNSS -26.7) | "
              f"std(asc-desc) 1 km {row['D_std_1km']:.2f}, 10 km {row['D_std_10km']:.2f} cm/yr")
    ch = out["no_iono"]["D_std_10km"] / out["with_iono"]["D_std_10km"] - 1
    print(f"P1: 10 km D std change without iono {ch * 100:+.0f}% (need <= -15%)")
    json.dump(out, open(ROOT / "data/processed/mexico-city/h8b_no_iono.json", "w"), indent=1)


if __name__ == "__main__":
    main()
