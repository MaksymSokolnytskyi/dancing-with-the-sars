"""Pack one scene in one look into a single self-contained HTML file you can send to someone.

    uv run python scripts/pack_scene.py mexico-city-asc --theme scroll --drop built --show city
    -> share/mexico-city-asc-scroll.html

Everything the app fetches (scenes list, theme, manifest, layers) is gzip+base64'd into <script type="text/plain">
blocks and served to the app by a small getJSON shim (DecompressionStream). The scene and style pickers are hidden.
Still needs internet for the map/deck.gl libraries (unpkg) and the Google font; no server is needed: open the file.

--drop  layer ids to leave out (smaller file); --show  layer ids or group names to switch on by default.
"""
from __future__ import annotations

import argparse
import base64
import gzip
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
WEB = ROOT / "web"


def blob(path: str, obj) -> str:
    raw = json.dumps(obj, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    b64 = base64.b64encode(gzip.compress(raw, 9, mtime=0)).decode("ascii")
    return f'<script type="text/plain" data-pack="{path}">{b64}</script>\n'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("scene")
    ap.add_argument("--theme", default="scroll")
    ap.add_argument("--drop", nargs="*", default=[])
    ap.add_argument("--show", nargs="*", default=[])
    ap.add_argument("--out")
    a = ap.parse_args()

    idx = json.loads((WEB / "data/scenes.json").read_text(encoding="utf-8"))
    idx["scenes"] = [s for s in idx["scenes"] if s["id"] == a.scene]
    assert idx["scenes"], f"scene {a.scene} not in web/data/scenes.json"
    man = json.loads((WEB / f"data/{a.scene}/manifest.json").read_text(encoding="utf-8"))
    man["layers"] = [l for l in man["layers"] if l["id"] not in a.drop]
    for l in man["layers"]:
        if l["id"] in a.show or l.get("group") in a.show:
            l["visible"] = True
    for l in man["layers"]:                          # rasters travel inside the manifest as data URLs
        if l["kind"] == "image":
            f = WEB / f"data/{a.scene}/{l['file']}"
            mime = "image/png" if f.suffix == ".png" else "image/jpeg"
            l["url"] = f"data:{mime};base64," + base64.b64encode(f.read_bytes()).decode("ascii")
    theme = json.loads((WEB / f"themes/{a.theme}.json").read_text(encoding="utf-8"))

    packs = [blob("data/scenes.json", idx), blob(f"themes/{a.theme}.json", theme),
             blob(f"data/{a.scene}/manifest.json", man)]
    if man.get("story"):                             # the plate's pages travel too
        packs.append(blob(f"data/{a.scene}/{man['story']}",
                          json.loads((WEB / f"data/{a.scene}/{man['story']}").read_text(encoding="utf-8"))))
    sizes = {}
    for l in man["layers"]:
        if l["kind"] == "image":
            sizes[l["id"]] = len(l["url"]); continue
        f = WEB / f"data/{a.scene}/{l['file']}"
        b = blob(f"data/{a.scene}/{l['file']}", json.loads(f.read_text(encoding="utf-8")))
        packs.append(b); sizes[l["id"]] = len(b)

    js = (WEB / "app.js").read_text(encoding="utf-8")
    old = re.search(r"^const getJSON = .*$", js, re.M)
    assert old, "getJSON line not found in app.js"
    js = js.replace(old.group(0), '''const getJSON = async (u) => {           // packed build: data lives in <script data-pack> blocks
  const el = document.querySelector(`script[data-pack="${u.split("?")[0]}"]`);
  if (!el) throw new Error("not packed: " + u);
  const bytes = Uint8Array.from(atob(el.textContent.trim()), (c) => c.charCodeAt(0));
  return JSON.parse(await new Response(new Blob([bytes]).stream().pipeThrough(new DecompressionStream("gzip"))).text());
};''', 1)
    js = re.sub(r"^const THEMES = \[.*?\];", f'const THEMES = ["{a.theme}"];', js, count=1, flags=re.M)
    assert f'const THEMES = ["{a.theme}"];' in js
    css = (WEB / "style.css").read_text(encoding="utf-8")
    css += "\n/* packed build: one place, one look */\n#panel section:has(#scene) { display: none; }\n"

    html = (WEB / "index.html").read_text(encoding="utf-8")
    html, n1 = re.subn(r'<link rel="stylesheet" href="style\.css[^"]*">', lambda m: f"<style>\n{css}</style>", html)
    html, n2 = re.subn(r'<script type="module" src="app\.js[^"]*"></script>',
                       lambda m: "".join(packs) + f'<script type="module">\n{js}</script>', html)
    assert n1 == 1 and n2 == 1, (n1, n2)
    title = man["title"]["en"] if isinstance(man["title"], dict) else man["title"]
    html = re.sub(r"<title>.*?</title>", f"<title>{title} · Dancing with the SARs</title>", html, count=1)

    out = Path(a.out) if a.out else ROOT / "share" / f"{a.scene}-{a.theme}.html"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(html, encoding="utf-8")
    top = sorted(sizes.items(), key=lambda kv: -kv[1])[:5]
    print(out, f"{out.stat().st_size / 1e6:.1f} MB;", "biggest:", ", ".join(f"{k} {v / 1e6:.1f}" for k, v in top))


if __name__ == "__main__":
    main()
