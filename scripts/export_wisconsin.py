"""Wisconsin plate (sections W, WS, TW of notes/experiments.md) -> web/data/wisconsin/.

    uv run python scripts/export_wisconsin.py [dir with story_cells.json, story_harvest.json]

The input files are written by `experiments/ws_tests.py` into data/processed/wisconsin/ws/ (not in git; see
notes/wisconsin-story.md for how to regenerate them). A copy sent by the analyst sits in data/raw/from-team/wisconsin/,
which is the default here. No radar is read: this only reshapes the test outputs for the map.

Writes:
  blocks.json    cells layer: 1 km blocks, HV and HH per date relative to the block's own Jun-Sep 2026 median on the
                 same track (time modes), and the block's land class (categorical)
  harvest.json   cells layer: crop blocks, "standing" until the harvest date, "cut" after; radar date = median over the
                 tracks that saw one, optical date = Sentinel-2 (WS4)
  story.json     copied from sites/wisconsin.story.json
  manifest.json
"""
from __future__ import annotations

import datetime as dt
import json
import statistics
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "data/raw/from-team/wisconsin"
OUT = ROOT / "web/data/wisconsin"
BBOX = [-91.00, 45.02, -90.56, 45.28]
CLASSES = [("mixed", "Mixed (no class over half the block)", "Змішаний (жоден клас не займає й половини)"),
           ("crop", "Fields", "Поля"), ("corn", "Corn", "Кукурудза"), ("soy", "Soybean", "Соя"),
           ("upland", "Upland forest", "Ліс на сухому"), ("swamp", "Swamp forest", "Болотний ліс"), ("lake", "Lake", "Озеро")]


def _day(s):
    return dt.date.fromisoformat(s)


def main(src: Path = SRC):
    cells = json.loads((src / "story_cells.json").read_text(encoding="utf-8"))
    harv = json.loads((src / "story_harvest.json").read_text(encoding="utf-8"))
    dates = cells["dates"]
    OUT.mkdir(parents=True, exist_ok=True)
    keys = [c[0] for c in CLASSES]
    rows = [[c[0], c[1], keys.index(c[2]), c[4], c[3]] for c in cells["cells"]]          # hv first: it is the default series
    blocks = {
        "cell_km": cells["cell_km"], "columns": ["lon", "lat", "class", "hv_anom_db", "hh_anom_db"],
        "dates": dates, "tracks": cells["tracks"], "states": cells["states"], "note": cells["note"],
        "series": "hv_anom_db", "series_unit": "dB vs own summer", "default_mode": "hv_anom_db",
        "modes": [
            {"key": "hv_anom_db", "series": "hv_anom_db", "ramp": "anom_db", "time": True,
             "label": {"en": "HV echo vs its own summer (dB)", "uk": "Відлуння HV проти власного літа (дБ)"}},
            {"key": "hh_anom_db", "series": "hh_anom_db", "ramp": "anom_db", "time": True,
             "label": {"en": "HH echo vs its own summer (dB)", "uk": "Відлуння HH проти власного літа (дБ)"}},
            {"key": "class", "categorical": "wis_class", "label": {"en": "Land class", "uk": "Клас поверхні"},
             "names": [{"en": c[1], "uk": c[2]} for c in CLASSES]},
        ],
        "cells": rows,
    }
    (OUT / "blocks.json").write_text(json.dumps(blocks, ensure_ascii=False, separators=(",", ":")), encoding="utf-8")

    hrows, n_radar = [], 0
    for lon, lat, by_track, optical in harv["cells"]:
        rd = [_day(v) for v in by_track.values() if v]
        radar = None
        if rd:
            n_radar += 1
            radar = dt.date.fromordinal(int(round(statistics.median(d.toordinal() for d in rd))))
        cut = lambda when: [None if when is None else int(_day(s) >= when) for s in dates]    # noqa: E731
        hrows.append([lon, lat, cut(radar), cut(_day(optical) if optical else None)])
    harvest = {
        "cell_km": cells["cell_km"], "columns": ["lon", "lat", "radar_cut", "optical_cut"], "dates": dates,
        "default_mode": "radar_cut",
        "modes": [
            {"key": "radar_cut", "series": "radar_cut", "categorical": "harvest", "time": True,
             "label": {"en": "harvest seen by the radar", "uk": "жнива, які бачить радар"},
             "names": [{"en": "standing", "uk": "стоїть"}, {"en": "cut", "uk": "скошено"}]},
            {"key": "optical_cut", "series": "optical_cut", "categorical": "harvest", "time": True,
             "label": {"en": "harvest seen by Sentinel-2", "uk": "жнива, які бачить Sentinel-2"},
             "names": [{"en": "green", "uk": "зелене"}, {"en": "turned", "uk": "пожовкло"}]},
        ],
        "cells": hrows,
    }
    (OUT / "harvest.json").write_text(json.dumps(harvest, ensure_ascii=False, separators=(",", ":")), encoding="utf-8")

    story = json.loads((ROOT / "sites/wisconsin.story.json").read_text(encoding="utf-8"))
    (OUT / "story.json").write_text(json.dumps(story, ensure_ascii=False, indent=1), encoding="utf-8")
    w, s, e, n = BBOX
    man = {
        "engine_version": 1, "id": "wisconsin", "title": story["title"],
        "subtitle": "NISAR L-band GCOV, four tracks: a year of forest, fields, swamps and lakes",
        "bbox": BBOX, "center": [(w + e) / 2, (s + n) / 2], "dates": sorted(set(dates)), "n_pairs": 0,
        "direction": "ASCENDING+DESCENDING", "reference": {"bbox": None, "label": ""},
        "anchor": {"type": "radar-only", "gnss": []}, "trust": {"snr_strong": 3.0, "snr_weak": 1.0, "segments": True},
        "stats": {"blocks": {"cells": len(rows)}, "harvest": {"cells": len(hrows), "with_radar_date": n_radar}},
        "layers": [
            {"id": "blocks", "kind": "cells", "file": "blocks.json", "label": {"en": "1 km blocks", "uk": "Блоки 1 км"},
             "props": blocks["columns"]},
            {"id": "harvest", "kind": "cells", "file": "harvest.json", "label": {"en": "Crop blocks: harvest", "uk": "Блоки полів: жнива"},
             "props": harvest["columns"]},
        ],
        "story": "story.json",
        "notes": ["Built by scripts/export_wisconsin.py from the WS test outputs (notes/experiments.md, sections W, WS, TW)."],
        "generated": dt.datetime.now().isoformat(timespec="seconds"),
    }
    (OUT / "manifest.json").write_text(json.dumps(man, ensure_ascii=False, indent=1), encoding="utf-8")
    sc_path = ROOT / "web/data/scenes.json"
    sc = json.loads(sc_path.read_text(encoding="utf-8"))
    entry = {"id": "wisconsin", "title": story["title"], "subtitle": man["subtitle"], "kind": "cells"}
    ids = [x["id"] for x in sc["scenes"]]
    if "wisconsin" in ids:                                 # keep its place and hand-written keys (title, blurb)
        sc["scenes"][ids.index("wisconsin")] = {**entry, **sc["scenes"][ids.index("wisconsin")]}
    else:
        sc["scenes"].append(entry)
    sc_path.write_text(json.dumps(sc, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
    import add_context                                       # OSM context + cartouche subtitle from sites/wisconsin.json
    add_context.main("wisconsin")
    print(f"wisconsin: {len(rows)} blocks, {len(hrows)} crop blocks ({n_radar} with a radar harvest date), {len(dates)} dates")


if __name__ == "__main__":
    main(Path(sys.argv[1]) if len(sys.argv) > 1 else SRC)
