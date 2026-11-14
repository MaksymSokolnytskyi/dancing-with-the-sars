"""K diagnostic (post-hoc): a control lake for the K levels (notes/experiments.md, section K).

Open water of Lake Balaton (Hungary), same product type (GCOV DHDH band B) and season as the Kakhovka summer dates.
uv run python experiments/k_diag.py      -> data/processed/kakhovka/k_diag.json
"""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src")]
from hachure import gcov  # noqa: E402

BOX = (17.72, 46.84, 17.78, 46.88)           # open water only, away from the shores
OUT = ROOT / "data/processed/kakhovka"


def main():
    lon, lat = (BOX[0] + BOX[2]) / 2, (BOX[1] + BOX[3]) / 2
    scenes = [s for s in gcov.search(lon, lat, start="2026-07-01", end="2026-09-30") if s.pols == "DHDH"][:3]
    rows = []
    for p in gcov.fetch(scenes, BOX, freq="B"):
        w = gcov.Window.load(p)
        ok = (w.mask > 0) & (w.mask < 255)
        hh, hv, neb = (10 * np.log10(np.median(v[ok])) for v in (w.hh, w.hv, w.extra["neb_hv"]))
        rows.append({"date": w.date.date().isoformat(), "track": int(w.name.split("_")[5]),
                     "water_hh_db": round(float(hh), 1), "water_hv_db": round(float(hv), 1), "neb_hv_db": round(float(neb), 1),
                     "hv_minus_neb_db": round(float(hv - neb), 1)})
        print(rows[-1])
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "k_diag.json").write_text(json.dumps({"box": BOX, "rows": rows}, indent=1), encoding="utf-8")


if __name__ == "__main__":
    main()
