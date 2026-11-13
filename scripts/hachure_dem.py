"""Render a Dufour-style hachure map from Copernicus GLO-30 DEM.

Example (default = Dents du Midi / Martigny, sheet XVII of the Dufour map):
    uv run python scripts/hachure_dem.py
    uv run python scripts/hachure_dem.py --name everest --bbox 86.80 27.92 87.02 28.06
"""
from __future__ import annotations

import argparse
import time
from pathlib import Path

from hachure import dem, export, field, trace

ROOT = Path(__file__).resolve().parents[1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", default="dents-du-midi")
    ap.add_argument("--bbox", nargs=4, type=float, default=[6.85, 46.08, 7.12, 46.25],
                    metavar=("W", "S", "E", "N"))
    ap.add_argument("--res", type=float, default=30.0, help="grid metres per pixel")
    ap.add_argument("--interval", type=float, default=50.0, help="contour band, metres")
    ap.add_argument("--spacing", type=float, default=1.6, help="stroke spacing, grid px")
    ap.add_argument("--sigma-dir", type=float, default=90.0)
    ap.add_argument("--sigma-tone", type=float, default=45.0)
    ap.add_argument("--scale", type=float, default=2.0, help="output px per grid px")
    ap.add_argument("--wmax", type=float, default=1.6)
    args = ap.parse_args()

    raw = ROOT / "data" / "raw" / f"copdem_{args.name}_{int(args.res)}m.tif"
    out = ROOT / "data" / "processed" / args.name
    out.mkdir(parents=True, exist_ok=True)

    t = time.time()
    print(f"DEM -> {dem.fetch_copdem(tuple(args.bbox), raw, res=args.res)}  ({time.time() - t:.1f}s)")
    z, transform, crs = dem.load(raw)
    print(f"grid {z.shape}, z {z.min():.0f}..{z.max():.0f} m")

    t = time.time()
    f = field.prepare(z, args.res, sigma_dir_m=args.sigma_dir, sigma_tone_m=args.sigma_tone)
    strokes = trace.hachures(f, interval_m=args.interval, spacing_px=args.spacing)
    print(f"{len(strokes)} strokes  ({time.time() - t:.1f}s)")

    print(export.to_png(strokes, z.shape, out / "hachures.png", scale=args.scale, w_max=args.wmax))
    print(export.to_svg(strokes, z.shape, out / "hachures.svg", scale=args.scale, w_max=args.wmax))
    print(export.to_geojson(strokes, transform, crs, out / "hachures.geojson"))


if __name__ == "__main__":
    main()
