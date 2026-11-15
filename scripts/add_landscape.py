"""Add the static landscape (relief hachures, woods, built-up land, roads, places) to an already built plate.

    uv run python scripts/add_landscape.py wisconsin [central-valley-asc ...]

The block plates (Wisconsin, Central Valley) are built by their own scripts and carried only their radar blocks and a
few lines; this draws the same land as the Mexico City plate around them, from the site config (`terrain`, `landcover`,
`context`), and adds the layer ids to the manifest's `base` list, so every story page keeps them on. Radar layers and
the story are left as they are. Needs the Copernicus DEM, the WorldCover clip and the OSM files (`sars fetch`,
scripts/fetch_worldcover.py, scripts/fetch_osm.py).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from hachure import engine  # noqa: E402

BASE = ["terrain", "trees", "built", "roads", "rail", "water", "rivers", "canals", "places", "border", "peaks"]


def main():
    for sid in sys.argv[1:]:
        cfg = engine.load_site(ROOT / "sites" / f"{sid}.json")
        out = ROOT / "web" / "data" / sid
        mp = out / "manifest.json"
        man = json.loads(mp.read_text(encoding="utf-8"))
        new = []
        if cfg.get("terrain") is not None:
            z, T, crs = engine._dem(cfg)
            fc = engine.terrain_layer(cfg, z, T, crs)
            engine._write(out, "terrain.geojson", fc)
            new.append({"id": "terrain", "kind": "hachures-static", "file": "terrain.geojson",
                        "label": {"en": "Relief", "uk": "Рельєф"}, "props": ["dark"]})
            print(f"{sid}: terrain {len(fc['features'])} strokes")
        missing = [c["id"] for c in cfg.get("context", []) if not (ROOT / c["src"]).exists()]
        if missing:
            print(f"{sid}: no source yet for {missing}, left out")
        cfg["context"] = [c for c in cfg.get("context", []) if c["id"] not in missing]
        new += engine.context_layers(cfg, out)
        if cfg.get("landcover") and (ROOT / cfg["landcover"]["src"]).exists():
            new += engine.landcover_layers(cfg, out)
        ids = {l["id"] for l in new}
        man["layers"] = [l for l in man["layers"] if l["id"] not in ids] + new
        have = {l["id"] for l in man["layers"]}
        man["base"] = [b for b in BASE if b in have]
        mp.write_text(json.dumps(man, indent=1, ensure_ascii=False), encoding="utf-8")
        print(f"{sid}: base {man['base']}")


if __name__ == "__main__":
    main()
