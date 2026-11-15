"""Wash the two countries of a border plate in two tones, as old political maps and strategy-game maps do.

    uv run python scripts/add_countries.py ostrava-karvina-asc

Splits the plate's frame along the drawn border (web/data/<plate>/ctx_border.geojson, from OpenStreetMap) and writes
ctx_countries.geojson: a pale wash over each country and a stronger ribbon just inside the border (tone 0/1 = wash,
2/3 = ribbon; the theme's polygons.country.fills give the colours). The country names go into the places layer.
Which side is which is told by one town on each side (SIDES below). Adds "countries" to the manifest and its base list.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

from shapely.geometry import Point, box, mapping, shape
from shapely.ops import linemerge, split, unary_union

ROOT = Path(__file__).resolve().parents[1]
# per plate: [name, a town on that side, where the name is written]
SIDES = {"ostrava-karvina-asc": [[{"en": "CZECHIA", "uk": "ЧЕХІЯ"}, (18.28, 49.83), (18.45, 49.775)],
                                 [{"en": "POLAND", "uk": "ПОЛЬЩА"}, (18.60, 49.95), (18.73, 49.905)]]}
RIBBON_DEG = 0.006                                                # ribbon width, ~450 m


def main():
    for sid in sys.argv[1:]:
        out = ROOT / "web" / "data" / sid
        man = json.loads((out / "manifest.json").read_text(encoding="utf-8"))
        border = json.loads((out / "ctx_border.geojson").read_text(encoding="utf-8"))
        line = linemerge(unary_union([shape(f["geometry"]) for f in border["features"]]))
        frame = box(*man["bbox"])
        parts = list(split(frame, line).geoms)
        ribbon = line.buffer(RIBBON_DEG)
        feats, places = [], json.loads((out / "ctx_places.geojson").read_text(encoding="utf-8"))
        places["features"] = [f for f in places["features"] if not f["properties"].get("country")]
        for k, (name, town, at) in enumerate(SIDES[sid]):
            side = unary_union([g for g in parts if g.contains(Point(town)) or g.distance(Point(town)) < 1e-9])
            if side.is_empty:                                     # the border does not cut the frame cleanly: all parts not of the other town
                raise SystemExit(f"{sid}: could not find the side of {name['en']}")
            feats.append({"type": "Feature", "properties": {"tone": k, "name": name["en"]}, "geometry": mapping(side.simplify(1e-4))})
            feats.append({"type": "Feature", "properties": {"tone": k + 2}, "geometry": mapping(side.intersection(ribbon).simplify(1e-4))})
            places["features"].append({"type": "Feature", "properties": {"name": name["en"], "rank": 5, "country": True},
                                       "geometry": {"type": "Point", "coordinates": list(at)}})
        (out / "ctx_countries.geojson").write_text(json.dumps({"type": "FeatureCollection", "features": feats}, separators=(",", ":")), encoding="utf-8")
        (out / "ctx_places.geojson").write_text(json.dumps(places, separators=(",", ":"), ensure_ascii=False), encoding="utf-8")
        man["layers"] = [l for l in man["layers"] if l["id"] != "countries"]
        man["layers"].insert(0, {"id": "countries", "kind": "polygons", "file": "ctx_countries.geojson", "style": "country", "z": 0.5,
                                 "label": {"en": "Countries", "uk": "Країни"}, "attribution": "© OpenStreetMap contributors"})
        man["base"] = ["countries"] + [b for b in man.get("base", []) if b != "countries"]
        (out / "manifest.json").write_text(json.dumps(man, indent=1, ensure_ascii=False), encoding="utf-8")
        print(f"{sid}: {len(parts)} parts of the frame, countries + ribbons written")


if __name__ == "__main__":
    main()
