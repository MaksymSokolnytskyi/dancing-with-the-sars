"""Test O3b (notes/experiments.md, pre-registered 2026-11-14): do the NISAR bowls show up in EGMS 2019–2023?

    uv run python scripts/ostrava/egms_o3b.py

Inputs: output/ostrava_karvina/zones_vertical_osm.geojson (the study's bowls) and the EGMS L3 Ortho vertical tiles from
scripts/ostrava/egms_fetch.py (data/ostrava/egms/*.tiff, mean velocity in mm/yr, 100 m, EPSG:3035).
For each bowl: the mean EGMS vertical velocity over the 100 m cells inside the polygon (all_touched=False), the number of
cells with data and their share of the polygon.

Pre-registered prediction: >= 8 of the 10 PL bowls of O1 (M01–M09, M18) have EGMS mean <= -2 cm/yr, and <= 1 of the 4 CZ
bowls (M13, M14, M16, M20) does. Pass: both halves. Partial: one. Fail: neither.
A bowl with no EGMS cell inside counts as "not <= -2 cm/yr" (no confirmation); this reading is written here, not in O3b,
so it is reported separately.
POST-HOC (added after seeing that most bowls have no EGMS cell inside): the same statistics in a ring 0–500 m outside
each bowl (columns ring_*), and the share of all EGMS cells in the two tiles at <= -2 cm/yr, for comparison.
Writes output/ostrava_karvina/egms_o3b.csv.
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np
import rasterio
from rasterio.features import geometry_mask
from rasterio.warp import transform_geom
from scipy.ndimage import binary_dilation

ROOT = Path(__file__).resolve().parents[2]
ZONES = ROOT / "output" / "ostrava_karvina" / "zones_vertical_osm.geojson"
EGMS = sorted((ROOT / "data" / "ostrava" / "egms").glob("EGMS_L3_*_U_*.tiff"))
OUT = ROOT / "output" / "ostrava_karvina" / "egms_o3b.csv"
PL10 = ["M01", "M02", "M03", "M04", "M05", "M06", "M07", "M08", "M09", "M18"]
CZ4 = ["M13", "M14", "M16", "M20"]
THR_MM = -20.0                       # -2 cm/yr


def main() -> int:
    if not EGMS:
        print("no EGMS tiles: run scripts/ostrava/egms_fetch.py first")
        return 1
    arrays = []
    for path in EGMS:
        t = rasterio.open(path)
        a = t.read(1).astype(float)
        a[a == t.nodata] = np.nan
        arrays.append((t, a))
    rows = []
    for f in json.loads(ZONES.read_text(encoding="utf-8"))["features"]:
        p = f["properties"]
        g = transform_geom("EPSG:4326", "EPSG:3035", f["geometry"])
        vals, cells, rvals, rcells = [], 0, [], 0
        for t, a in arrays:
            inside = ~geometry_mask([g], out_shape=a.shape, transform=t.transform, invert=False, all_touched=False)
            cells += int(inside.sum())
            vals.extend(a[inside & np.isfinite(a)].tolist())
            if inside.any():
                ring = binary_dilation(inside, iterations=5) & ~inside          # 5 cells of 100 m
                rcells += int(ring.sum())
                rvals.extend(a[ring & np.isfinite(a)].tolist())
        vals, rvals = np.array(vals), np.array(rvals)
        mean = float(vals.mean()) if vals.size else None
        rows.append({"id": p["id"], "country": p["country"], "nisar_up_cm_yr": p["up_mean_cm_yr"],
                     "n_tracks": p["n_tracks_confirm"], "egms_cells_with_data": int(vals.size), "polygon_cells": cells,
                     "egms_share": round(vals.size / cells, 2) if cells else 0.0,
                     "egms_mean_cm_yr": round(mean / 10, 2) if mean is not None else None,
                     "egms_min_cm_yr": round(float(vals.min()) / 10, 2) if vals.size else None,
                     "sinks_le_2cm": bool(mean is not None and mean <= THR_MM),
                     "ring_cells_with_data": int(rvals.size), "ring_cells": rcells,
                     "ring_median_cm_yr": round(float(np.median(rvals)) / 10, 2) if rvals.size else None})
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    by = {r["id"]: r for r in rows}
    for r in rows:
        print(f"{r['id']} {r['country']}  NISAR {r['nisar_up_cm_yr']:7.1f}  EGMS mean "
              f"{'   n/a' if r['egms_mean_cm_yr'] is None else format(r['egms_mean_cm_yr'], '6.2f')} cm/yr  "
              f"cells {r['egms_cells_with_data']}/{r['polygon_cells']}  {'<= -2' if r['sinks_le_2cm'] else '     '}  "
              f"| post-hoc ring 0-500 m: median {'n/a' if r['ring_median_cm_yr'] is None else r['ring_median_cm_yr']} cm/yr "
              f"({r['ring_cells_with_data']}/{r['ring_cells']} cells)")
    pl = sum(by[i]["sinks_le_2cm"] for i in PL10)
    cz = sum(by[i]["sinks_le_2cm"] for i in CZ4)
    pl_nodata = [i for i in PL10 if by[i]["egms_mean_cm_yr"] is None]
    cz_nodata = [i for i in CZ4 if by[i]["egms_mean_cm_yr"] is None]
    pl_ok, cz_ok = pl >= 8, cz <= 1
    verdict = "PASS" if pl_ok and cz_ok else "PARTIAL" if pl_ok or cz_ok else "FAIL"
    print(f"\nPL half: {pl}/10 bowls with EGMS mean <= -2 cm/yr (need >= 8) -> {'pass' if pl_ok else 'fail'}; no EGMS data: {pl_nodata or 'none'}")
    print(f"CZ half: {cz}/4 bowls with EGMS mean <= -2 cm/yr (need <= 1) -> {'pass' if cz_ok else 'fail'}; no EGMS data: {cz_nodata or 'none'}")
    allv = np.concatenate([a[np.isfinite(a)] for _, a in arrays])
    print(f"O3b: {verdict}   (written to {OUT.relative_to(ROOT)})")
    print(f"post-hoc: all EGMS cells in the two tiles: median {np.median(allv) / 10:+.2f} cm/yr, "
          f"{(allv <= THR_MM).mean():.1%} at <= -2 cm/yr")
    return 0


if __name__ == "__main__":
    sys.exit(main())
