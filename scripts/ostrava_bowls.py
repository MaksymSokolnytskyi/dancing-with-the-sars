"""Oleh's subsidence bowls (PR #5, block O) -> two outline layers for the Ostrava-Karvina plate.

    uv run python scripts/ostrava_bowls.py [path/to/zones_vertical_osm.geojson]

Input: zones_vertical_osm.geojson from Oleh's three-track run (A057, A158, D152), kept outside git in
data/raw/from-team/ostrava/ostrava_karvina/. Output, into web/data/ostrava-asc/:
  bowls.geojson            reliable bowls by the fixed rule O3a (dashed ink outline + "M01" label)
  bowls_candidate.geojson  the other bowls (dotted, fainter)
Rule O3a (notes/experiments.md, fixed 2026-11-14): confirmed by 3 tracks AND no-data share < 0.6 AND water share < 0.10.
The plate uses this rule, not the hand-picked "reliable" column of the first run (they differ in 8 bowls).
The rings become LineStrings (the old-map sign for a mining boundary is a broken line), plus one labelled point per bowl.
Field values in Oleh's file are partly Ukrainian ("так"/"ні", land-use text); only numbers and English labels are kept.
The layers are listed under "extra_layers" in sites/ostrava-asc.json, so `sars build ostrava-asc` keeps them.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "data/raw/from-team/ostrava/ostrava_karvina/zones_vertical_osm.geojson"
OUT = ROOT / "web/data/ostrava-asc"


def rings(geom):
    polys = [geom["coordinates"]] if geom["type"] == "Polygon" else geom["coordinates"]
    return [[[round(x, 5), round(y, 5)] for x, y in ring] for p in polys for ring in p]


def main(src: Path = SRC):
    fc = json.loads(src.read_text(encoding="utf-8"))
    out = {"confirmed": [], "candidate": []}
    for f in fc["features"]:
        p = f["properties"]
        sure = int(p["n_tracks_confirm"]) >= 3 and p["no_data_frac"] < 0.6 and p["water_frac"] < 0.10      # rule O3a
        props = {"id": p["id"], "country": p.get("country"), "tracks": int(p["n_tracks_confirm"]),
                 "up_mean_cm_yr": p["up_mean_cm_yr"], "up_p10_cm_yr": p["up_p10_cm_yr"], "area_km2": p["area_km2"],
                 "no_data_frac": p["no_data_frac"], "buildings": p.get("buildings"), "status": "confirmed" if sure else "candidate"}
        key = "confirmed" if sure else "candidate"
        for r in rings(f["geometry"]):
            out[key].append({"type": "Feature", "properties": props, "geometry": {"type": "LineString", "coordinates": r}})
        out[key].append({"type": "Feature", "properties": {**props, "name": p["id"], "rank": 1},
                         "geometry": {"type": "Point", "coordinates": [round(p["lon"], 5), round(p["lat"], 5)]}})
    for key, name in (("confirmed", "bowls.geojson"), ("candidate", "bowls_candidate.geojson")):
        (OUT / name).write_text(json.dumps({"type": "FeatureCollection", "features": out[key]}, separators=(",", ":")),
                                encoding="utf-8")
        pts = [f["properties"] for f in out[key] if f["geometry"]["type"] == "Point"]
        for c in ("PL", "CZ"):
            q = [x for x in pts if x["country"] == c]
            a = sum(x["area_km2"] for x in q)
            wm = sum(x["area_km2"] * x["up_mean_cm_yr"] for x in q) / a if a else float("nan")
            print(f"{name} {c}: {len(q)} bowls {' '.join(x['id'] for x in q)}; {a:.2f} km2, area-weighted {wm:.1f} cm/yr, "
                  f"buildings {sum(x['buildings'] or 0 for x in q)}")


if __name__ == "__main__":
    main(Path(sys.argv[1]) if len(sys.argv) > 1 else SRC)
