"""Build a site for the web map: sites/<id>.json -> data/processed/<id>/engine_velocity.npz -> web/data/<id>/.

    uv run python scripts/build_site.py sites/mexico-city-asc.json            # reuse cached velocity if present
    uv run python scripts/build_site.py sites/mexico-city-asc.json --rebuild  # recompute rate + trust from the GUNW cache
    uv run python scripts/build_site.py sites/mexico-city-asc.json --preview  # also a quick-look PNG (default style)

New place = new JSON in sites/ + cached GUNW windows (scripts/nisar_gunw_rate.py fills the cache). No code changes.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

from hachure import engine

ROOT = Path(__file__).resolve().parents[1]


def preview(man_path: Path, out_png: Path):
    """Quick-look in the default engraving style: grey relief, red strokes, faded where trust is weak."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    d = man_path.parent
    man = json.loads(man_path.read_text(encoding="utf-8"))
    fig, ax = plt.subplots(figsize=(10, 9), dpi=110)
    fig.patch.set_facecolor("#f3eedf")
    ax.set_facecolor("#f3eedf")
    for lay in man["layers"]:
        if lay["kind"] in ("cells", "polygons", "card"):
            continue
        fc = json.loads((d / lay["file"]).read_text(encoding="utf-8"))
        if lay["kind"] == "hachures-static":
            segs = [f["geometry"]["coordinates"] for f in fc["features"]]
            ax.add_collection(LineCollection(segs, colors="#5f594d", linewidths=[0.1 + 0.9 * f["properties"]["dark"] for f in fc["features"]]))
        else:
            for cls, alpha in (("strong", 1.0), ("weak", 0.35)):
                fs = [f for f in fc["features"] if f["properties"]["trust"] == cls]
                ax.add_collection(LineCollection([f["geometry"]["coordinates"] for f in fs], colors="#a8322a", alpha=alpha,
                                                 linewidths=[0.3 + 1.6 * f["properties"]["dark"] for f in fs]))
    w, s, e, n = man["bbox"]
    ax.set_xlim(w, e)
    ax.set_ylim(s, n)
    ax.set_aspect(1 / np.cos(np.radians((s + n) / 2)))
    ax.axis("off")
    st = man["stats"].get("velocity", {})
    gn = man["anchor"]["gnss"]
    lab = "radar only" if man["anchor"]["type"] == "radar-only" else "tied to GNSS"
    chk = "; ".join(f"{g['station']}: NISAR {g['nisar_cm_yr']:+.1f} vs GNSS {g['rate_cm_yr']:+.1f} cm/yr"
                    for g in gn if g.get("nisar_cm_yr") is not None)
    ax.set_title(f"{man['title']} · {man['dates'][0]} → {man['dates'][-1]} · {man['n_pairs']} pairs\n"
                 f"{lab}; {chk}\nstrokes {st.get('strokes')}, trust {st.get('trust')}, faded = 1 ≤ SNR < 3, "
                 f"hidden = SNR < 1", fontsize=9, family="serif")
    fig.savefig(out_png, facecolor=fig.get_facecolor(), bbox_inches="tight")
    plt.close(fig)
    return out_png


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("site")
    ap.add_argument("--rebuild", action="store_true")
    ap.add_argument("--preview", action="store_true")
    ap.add_argument("--keep-terrain", action="store_true", help="reuse web/data/<id>/terrain.geojson if present")
    ap.add_argument("--out", default=str(ROOT / "web" / "data"))
    a = ap.parse_args()
    cfg = engine.load_site(a.site)
    cfg["_keep_terrain"] = a.keep_terrain
    proc = ROOT / "data/processed" / cfg["id"]
    proc.mkdir(parents=True, exist_ok=True)
    vf = proc / "engine_velocity.npz"
    t = time.time()
    v = None
    if {"velocity", "grid"} & set(cfg["layers"]):
        if vf.exists() and not a.rebuild:
            v = engine.load_velocity(vf)
            print(f"{cfg['id']}: cached velocity ({v.n_pairs} pairs)")
        else:
            print(f"{cfg['id']}: computing velocity + trust from the GUNW cache")
            v = engine.velocity(cfg)
            engine.save_velocity(v, vf)
    if v is not None:
        v = engine.tropo_correct(cfg, v)        # no-op unless the site asks for it
        v = engine.highpass_correct(cfg, v)     # no-op unless the site asks for it
    man = engine.export_site(cfg, v, Path(a.out) / cfg["id"])
    # scene index for the front-end
    idx_p = Path(a.out) / "scenes.json"
    idx = json.loads(idx_p.read_text(encoding="utf-8")) if idx_p.exists() else {"scenes": []}
    entry = {"id": cfg["id"], "title": cfg.get("title", cfg["id"]), "subtitle": cfg.get("subtitle", ""), "kind": "+".join(cfg["layers"])}
    ids = [s["id"] for s in idx["scenes"]]
    if cfg["id"] in ids:                       # rebuilding keeps the scene's place in the menu
        old = idx["scenes"][ids.index(cfg["id"])]   # keep keys the build does not write (blurb, href)
        idx["scenes"][ids.index(cfg["id"])] = {**old, **entry}
    else:
        idx["scenes"].append(entry)
    idx_p.write_text(json.dumps(idx, indent=1, ensure_ascii=False), encoding="utf-8")
    print(f"{man}  ({time.time() - t:.0f}s)")
    if a.preview:
        print(preview(man, proc / "engine_preview.png"))


if __name__ == "__main__":
    main()
