"""Re-pack web/reykjanes-scroll.html on the current engine without its source data.

    uv run python scripts/repack_reykjanes.py

scripts/build_reykjanes_scroll.py needs the Reykjanes analysis products (--src), which live on one machine only. When
just the engine changes (web/app.js, style.css, index.html, themes/scroll.json) or the math window
(docs/reykjanes-math.html), this keeps the packed data blocks of the existing page and rebuilds everything around
them with the same patches the full build applies. The math window is kept only if the page already had it;
--no-math makes the jury build out of a test build (no math window, no page IX, no test badge):

    uv run python scripts/repack_reykjanes.py --no-math
"""
from __future__ import annotations

import argparse
import base64
import gzip
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import build_reykjanes_scroll as B  # noqa: E402

# story text written before the math window was translated (the full build writes the new wording itself)
STORY_FIX = [("(in Ukrainian)", "(in English)"), ("(українською)", "(англійською)"),
             ("Open the math window (Ukrainian)", "Open the math window")]


def unpack(s):
    return json.loads(gzip.decompress(base64.b64decode(s)).decode("utf-8"))


def towns(packs):
    """Town streets (OSM) and the built-up hatch (WorldCover) for Grindavík, Njarðvík, Vogar: the same town sign as on
    the other plates. Drawn from data/raw (scripts/fetch_osm.py, scripts/fetch_worldcover.py); skipped if missing."""
    import tempfile
    sys.path.insert(0, str(B.ROOT / "src"))
    from hachure import engine
    key = f"data/{B.SID}/manifest.json"
    man = unpack(packs[key])
    osm, wc = "data/raw/osm/reykjanes/streets.json", "data/raw/context/reykjanes_worldcover_10m.tif"
    cfg = {"bbox": man["bbox"], "context": [], "landcover": None}
    if (B.ROOT / osm).exists():
        cfg["context"].append({"id": "streets", "kind": "lines", "src": osm, "style": "street", "visible": True, "z": 1.15,
                               "label": {"en": "Town streets", "uk": "Вулиці міст"}, "attribution": "© OpenStreetMap contributors"})
    if (B.ROOT / wc).exists():
        cfg["landcover"] = {"src": wc, "attribution": "© ESA WorldCover 2021 (CC-BY 4.0)",
                            "built": {"class": 50, "block_m": 40, "min_share": 0.4, "spacing_m": 120, "min_run_m": 80, "visible": True}}
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp)
        new = engine.context_layers(cfg, out) + (engine.landcover_layers(cfg, out) if cfg["landcover"] else [])
        for e in new:
            packs[f"data/{B.SID}/{e['file']}"] = B.pack(json.loads((out / e["file"]).read_text(encoding="utf-8")))
    ids = {e["id"] for e in new}
    man["layers"] = [l for l in man["layers"] if l["id"] not in ids] + new
    man["base"] = sorted(ids)                                     # kept on through every page, like the landscape elsewhere
    packs[key] = B.pack(man)
    return sorted(ids)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-math", action="store_true", help="drop the math window and its page (jury build)")
    a = ap.parse_args()
    html = B.DST.read_text(encoding="utf-8")
    packs = dict(re.findall(r'<script type="text/plain" data-pack="([^"]+)">([^<]*)</script>', html))
    assert packs, "no packed data in " + str(B.DST)
    with_math = 'id="math-open"' in html and not a.no_math
    head, app, tail, theme = B.engine()
    packs["themes/scroll.json"] = B.pack(B.theme_for(theme))
    # the atlas index as it is now, so the last page leads on to the next plate (that plate opens in index.html)
    packs["data/scenes.json"] = B.pack(json.loads((B.ROOT / "web" / "data" / "scenes.json").read_text(encoding="utf-8")))
    # the land of the locator globe in the corner
    packs["data/world-land.json"] = B.pack(json.loads((B.ROOT / "web" / "data" / "world-land.json").read_text(encoding="utf-8")))
    key = f"data/{B.SID}/story.json"
    if key in packs:
        txt = json.dumps(unpack(packs[key]), ensure_ascii=False)
        for a, b in STORY_FIX:
            txt = txt.replace(a, b)
        st = B.fix_text(json.loads(txt))
        if not with_math:
            st["pages"] = [p for p in st["pages"] if p.get("section") != "math"]
        packs[key] = B.pack(st)
    added = towns(packs)
    blocks ="\n".join(f'  <script type="text/plain" data-pack="{k}">{v}</script>' for k, v in packs.items())
    math = B.MATH.read_text(encoding="utf-8") if with_math else ""
    out = B.patch_head(head, with_math) + blocks + '\n  <script type="module">' + B.patch_app(app, with_math) + "</script>\n" + math + tail
    B.DST.write_text(out, encoding="utf-8", newline="\n")
    print(f"re-packed {B.DST.relative_to(B.ROOT)}: {len(packs)} data blocks kept, math window {'kept' if with_math else 'none'}, towns {added}, "
          f"{B.DST.stat().st_size / 1e6:.1f} MB")


if __name__ == "__main__":
    main()
