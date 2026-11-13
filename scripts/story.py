"""Publish a plate's pages without rebuilding its layers (seconds instead of minutes): edit the text, run this, reload.

    uv run python scripts/story.py sites/mexico-city-asc.json

Copies the file named by the site's "story" key (e.g. sites/mexico-city-asc.story.json) to web/data/<id>/story.json and
points the existing manifest at it. A full scripts/build_site.py run does the same.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from hachure import engine

ROOT = Path(__file__).resolve().parents[1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("site")
    ap.add_argument("--out", default=str(ROOT / "web" / "data"))
    a = ap.parse_args()
    cfg = engine.load_site(a.site)
    out = Path(a.out) / cfg["id"]
    mf = out / "manifest.json"
    assert mf.exists(), f"{mf} missing: run scripts/build_site.py first"
    name = engine.write_story(cfg, out)
    assert name, f'{a.site} has no "story" key'
    man = json.loads(mf.read_text(encoding="utf-8"))
    man["story"] = name
    mf.write_text(json.dumps(man, indent=1, ensure_ascii=False), encoding="utf-8")
    pages = json.loads((out / name).read_text(encoding="utf-8"))["pages"]
    ids = {l["id"] for l in man["layers"]} | {l.get("group") for l in man["layers"]}
    for i, p in enumerate(pages, 1):
        unknown = [x for x in p.get("layers", []) if x not in ids]
        if unknown:
            print(f"  page {i}: unknown layers {unknown}")
    print(f"{out / name}: {len(pages)} pages")


if __name__ == "__main__":
    main()
