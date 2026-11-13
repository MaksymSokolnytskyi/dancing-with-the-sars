"""H9-P: per-pair iono screen -> its plane only. Mexico part (H8b code path) + any P-port site.

    uv run python scripts/h9_iono_plane.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path[:0] = [str(Path(__file__).parent), str(Path(__file__).parents[1] / "scripts")]
from h6_horizontal import smooth  # noqa: E402
from h8b_no_iono import ROOT, at_mmx1, rate  # noqa: E402

from hachure import gridding


def main():
    out = {}
    for iono in (True, "plane"):
        vd, xd, yd, ed = rate("DESCENDING", iono)
        va, xa, ya, ea = rate("ASCENDING", iono)
        dst = gridding.transform_from_centres(xd, yd)
        va_d = gridding.to_grid(va, gridding.transform_from_centres(xa, ya), f"EPSG:{ea}", dst, f"EPSG:{ed}", vd.shape)
        dx = abs(xd[1] - xd[0])
        row = dict(mmx1_desc=at_mmx1(vd, xd, yd, ed), mmx1_asc=at_mmx1(va, xa, ya, ea))
        for km in (1, 10):
            D = smooth(va_d - vd, km * 1000 / dx / 2.355)
            row[f"D_std_{km}km"] = float(np.nanstd(D - np.nanmedian(D)))
        key = "full_screen" if iono is True else "plane_only"
        out[key] = row
        print(f"{key:11s} MMX1 desc {row['mmx1_desc']:+.1f}  asc {row['mmx1_asc']:+.1f}  (GNSS -26.7) | "
              f"std(asc-desc) 1 km {row['D_std_1km']:.2f}, 10 km {row['D_std_10km']:.2f} cm/yr", flush=True)
    f, p = out["full_screen"], out["plane_only"]
    gap_f, gap_p = abs(f["mmx1_asc"] + 26.7), abs(p["mmx1_asc"] + 26.7)
    print(f"Mexico: D10 {p['D_std_10km']:.2f} (<= 1.16) -> {'PASS' if p['D_std_10km'] <= 1.16 else 'fail'}; "
          f"MMX1 asc gap {gap_f:.1f} -> {gap_p:.1f} (worse by <= 0.5) -> {'PASS' if gap_p - gap_f <= 0.5 else 'fail'}")
    json.dump(out, open(ROOT / "data/processed/mexico-city/h9_iono_plane.json", "w"), indent=1)


if __name__ == "__main__":
    main()
