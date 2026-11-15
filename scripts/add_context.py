"""Add (or refresh) a plate's OSM context layers without rebuilding its radar layers.

    uv run sars fetch <id> --only osm           # Overpass -> data/raw/osm/<id>/
    uv run python scripts/add_context.py <id>   # -> web/data/<id>/ctx_*.geojson + manifest entries

For plates whose radar layers come from other scripts (panama-quake, wisconsin) or take long to rebuild (central-valley-asc).
Also copies the site's `frame_sub` (cartouche subtitle) into the manifest.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from hachure import engine  # noqa: E402


def main(sid: str):
    cfg = engine.load_site(ROOT / "sites" / f"{sid}.json")
    cfg.setdefault("dem_res", 30)
    out = ROOT / "web/data" / sid
    man_p = out / "manifest.json"
    man = json.loads(man_p.read_text(encoding="utf-8"))
    new = engine.context_layers(cfg, out)
    ids = {e["id"] for e in new}
    man["layers"] = [l for l in man["layers"] if l["id"] not in ids] + new
    if cfg.get("frame_sub"):
        man["frame_sub"] = cfg["frame_sub"]
    man_p.write_text(json.dumps(man, ensure_ascii=False, indent=1), encoding="utf-8")
    print(f"{sid}: {len(new)} context layers -> {man_p.relative_to(ROOT)}")


if __name__ == "__main__":
    main(sys.argv[1])
