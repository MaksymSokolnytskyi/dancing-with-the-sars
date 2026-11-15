"""Ostrava–Karviná study -> context files for the atlas plate `sites/ostrava-karvina-asc.json`.

    uv run python scripts/ostrava/plate_layers.py

Reads the study outputs (`output/ostrava_karvina/zones_vertical_osm.geojson`, `osm_heaps.geojson`, written by
`ostrava_karvina.ipynb` and `ostrava_osm.mjs`) and writes `data/ostrava/plate/` (ignored by git):

- `bowls_reliable.geojson`, `bowls_candidate.geojson`: bowl outlines as lines + one label point per bowl ("M01 −68");
- `heaps.geojson`: spoil-heap outlines (OSM) as lines.

"Reliable" is the rule-based set of notes/experiments.md O3a, fixed before it was applied to the plate: confirmed by
3 tracks AND no-data share < 0.6 AND water share < 0.10. Everything else is a candidate. Prints the per-country numbers
the plate's pages quote.
"""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "output" / "ostrava_karvina"
OUT = ROOT / "data" / "ostrava" / "plate"


def reliable(p: dict) -> bool:
    """Rule O3a (notes/experiments.md)."""
    return p["n_tracks_confirm"] >= 3 and p["no_data_frac"] < 0.6 and p["water_frac"] < 0.10


def rings(geom: dict) -> list:
    if geom["type"] == "Polygon":
        return [geom["coordinates"][0]]
    if geom["type"] == "MultiPolygon":
        return [poly[0] for poly in geom["coordinates"]]
    return []


def line(coords, props):
    pts = [[round(x, 5), round(y, 5)] for x, y, *_ in coords]
    return {"type": "Feature", "properties": props, "geometry": {"type": "LineString", "coordinates": pts}}


def main():
    zones = json.loads((SRC / "zones_vertical_osm.geojson").read_text(encoding="utf-8"))["features"]
    OUT.mkdir(parents=True, exist_ok=True)
    sets = {"reliable": [], "candidate": []}
    stats = {}
    for f in zones:
        p = f["properties"]
        cls = "reliable" if reliable(p) else "candidate"
        rate = round(p["up_mean_cm_yr"])
        for r in rings(f["geometry"]):
            sets[cls].append(line(r, {"id": p["id"]}))
        sets[cls].append({"type": "Feature",
                          "properties": {"name": f"{p['id']} {'−' if rate < 0 else ''}{abs(rate)}", "rank": 2 if cls == "reliable" else 1},
                          "geometry": {"type": "Point", "coordinates": [round(p["lon"], 5), round(p["lat"], 5)]}})
        s = stats.setdefault((cls, p["country"]), {"bowls": [], "area_km2": 0.0, "buildings": 0, "roads_km": 0.0,
                                                   "rail_km": 0.0, "w": 0.0})
        s["bowls"].append(p["id"])
        s["area_km2"] += p["area_km2"]
        s["buildings"] += p["buildings"]
        s["roads_km"] += p["roads_km"]
        s["rail_km"] += p["rail_km"]
        s["w"] += p["area_km2"] * p["up_mean_cm_yr"]
    for cls, feats in sets.items():
        (OUT / f"bowls_{cls}.geojson").write_text(json.dumps({"type": "FeatureCollection", "features": feats},
                                                             ensure_ascii=False), encoding="utf-8")
    heaps = json.loads((SRC / "osm_heaps.geojson").read_text(encoding="utf-8"))["features"]
    hf = [line(r, {}) for f in heaps if f.get("geometry") for r in rings(f["geometry"])]
    (OUT / "heaps.geojson").write_text(json.dumps({"type": "FeatureCollection", "features": hf}), encoding="utf-8")
    print(f"wrote {OUT.relative_to(ROOT)}: {sum(1 for f in sets['reliable'] if f['geometry']['type'] == 'Point')} reliable, "
          f"{sum(1 for f in sets['candidate'] if f['geometry']['type'] == 'Point')} candidate bowls, {len(hf)} heap outlines")
    for (cls, c), s in sorted(stats.items()):
        print(f"  {cls:9s} {c}: {len(s['bowls'])} bowls ({' '.join(s['bowls'])}), {s['area_km2']:.2f} km2, "
              f"{s['buildings']} buildings, {s['roads_km']:.1f} km roads, {s['rail_km']:.1f} km rail, "
              f"area-weighted vertical {s['w'] / s['area_km2']:.1f} cm/yr")


if __name__ == "__main__":
    main()
