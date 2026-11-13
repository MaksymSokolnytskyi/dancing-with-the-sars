"""P-tie (pre-registered): GNSS-anchored rates, leave-one-station-out. Uses P-port rates (plane removed, full screen).

    uv run python scripts/p_tie.py --site central-valley
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path[:0] = [str(Path(__file__).parent), str(Path(__file__).parents[1] / "scripts")]
import p_port  # noqa: E402


def design(lon, lat, quad):
    x, y = lon - lon.mean(), lat - lat.mean()
    cols = [np.ones_like(x), x, y] + ([x * x, x * y, y * y] if quad else [])
    return np.column_stack(cols)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--site", required=True)
    ap.add_argument("--variant", default="")
    args = ap.parse_args()
    p_port.VAR = args.variant
    import csv
    pos = {r["station"]: (float(r["lon"]), float(r["lat"])) for r in csv.DictReader(open(p_port.P / "gnss" / f"{args.site}_rates.csv"))}
    for pas in ("asc", "desc"):
        rate, tr, crs = p_port.load_rate(args.site, pas)
        st = [s for s in p_port.at_stations(args.site, pas, rate, tr, crs) if np.isfinite(s[3])]
        g = np.array([s[1] for s in st]); n = np.array([s[3] for s in st])
        lon = np.array([pos[s[0]][0] for s in st]); lat = np.array([pos[s[0]][1] for s in st])
        d = n - g
        p1 = np.sqrt(np.mean((d - np.median(d)) ** 2))
        out = {}
        for quad in (False, True):
            if quad and len(st) - 1 < 15:
                continue
            res = []
            for i in range(len(st)):
                m = np.arange(len(st)) != i
                c, *_ = np.linalg.lstsq(design(lon[m], lat[m], quad), d[m], rcond=None)
                # predict at the left-out station with the same centring as the training set
                x0, y0 = lon[m].mean(), lat[m].mean()
                xi, yi = lon[i] - x0, lat[i] - y0
                row = np.array([1, xi, yi] + ([xi * xi, xi * yi, yi * yi] if quad else []))
                res.append(d[i] - row @ c)
            res = np.array(res)
            out["quadratic" if quad else "plane"] = (float(np.sqrt(np.mean(res ** 2))), float(np.median(np.abs(res))))
        print(f"{args.site}{args.variant} {pas}: {len(st)} stations | P1 (offset only) RMS {p1:.2f} | " +
              " | ".join(f"LOO {k}: RMS {v[0]:.2f}, median|res| {v[1]:.2f}" for k, v in out.items()))


if __name__ == "__main__":
    main()
