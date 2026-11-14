"""NASA FIRMS thermal hotspots (VIIRS S-NPP / NOAA-20 375 m, MODIS 1 km) over the Svartsengi area, Nov 2023 - Aug 2025.

Standard-processing archive via the FIRMS area API (key in the environment variable FIRMS_MAP_KEY, sent only to firms.modaps.eosdis.nasa.gov),
5-day requests (API limit). Checks against the known eruptions on the Sundhnukur crater row: detected or not, delay of the first
hotspot after the eruption start, distance of the hotspots to the Sundhnukur crater row line and to the dike trace from NPE.

Writes data/firms/firms_reykjanes.csv, output/reykjanes/ml/firms_hotspots.json and firms_hotspots.png.
Usage: python scripts/reykjanes/firms_hotspots.py
"""
from _secrets import read_env
import csv
import io
import json
import subprocess
import time
from datetime import date, datetime, timedelta
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from pyproj import Transformer

from progress import Progress
from reykjanes_analysis import PLACES

ROOT = Path(__file__).resolve().parents[2]
ML = ROOT / "output" / "reykjanes" / "ml"
OUT = ROOT / "data" / "firms"
BBOX = (-22.75, 63.76, -22.10, 63.98)                      # W, S, E, N
SOURCES = ("VIIRS_SNPP_SP", "VIIRS_NOAA20_SP", "MODIS_SP")
START, END = date(2023, 11, 1), date(2025, 8, 31)
ERUPTIONS = [("2023-12-18", "2023-12-21"), ("2024-01-14", "2024-01-16"), ("2024-02-08", "2024-02-09"), ("2024-03-16", "2024-05-09"),
             ("2024-05-29", "2024-06-22"), ("2024-08-22", "2024-09-05"), ("2024-11-20", "2024-12-09"), ("2025-04-01", "2025-04-01"),
             ("2025-07-16", "2025-08-05")]
# Sundhnukur crater row, approximate line (lon, lat): from north of Stora-Skogfell to Hagafell (geological maps; 2023-2024 fissures)
CRATER_ROW = ((-22.355, 63.905), (-22.420, 63.858))


def fetch():
    env = read_env()
    key = env["FIRMS_MAP_KEY"].strip().strip('"')
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    days = [(START + timedelta(days=d)) for d in range(0, (END - START).days + 1, 5)]
    with Progress("NASA FIRMS: теплові точки 2023–2025", total=len(days) * len(SOURCES), unit="запит") as pr:
        for src in SOURCES:
            for d in days:
                url = f"https://firms.modaps.eosdis.nasa.gov/api/area/csv/{key}/{src}/{','.join(map(str, BBOX))}/5/{d.isoformat()}"
                r = subprocess.run(["curl", "-s", "--max-time", "60", url], capture_output=True, text=True)
                txt = r.stdout
                if txt.startswith("latitude") or txt.startswith("country_id"):
                    for row in csv.DictReader(io.StringIO(txt)):
                        row["source"] = src
                        rows.append(row)
                elif txt.strip():
                    print(f"  {src} {d}: {txt[:120]!r}", flush=True)
                pr.update(1, note=f"{len(rows)} точок")
                time.sleep(0.3)
    del key
    if rows:
        keys = sorted({k for r in rows for k in r})
        with open(OUT / "firms_reykjanes.csv", "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=keys)
            w.writeheader()
            w.writerows(rows)
    return rows


def main():
    import sys
    rows = list(csv.DictReader(open(OUT / 'firms_reykjanes.csv', encoding='utf-8'))) if '--cached' in sys.argv else fetch()
    print(len(rows), "hotspots", flush=True)
    to_utm = Transformer.from_crs(4326, 32627, always_xy=True)
    sx, sy = to_utm.transform(*PLACES["Svartsengi"])
    lon = np.array([float(r["longitude"]) for r in rows]); lat = np.array([float(r["latitude"]) for r in rows])
    ux, uy = to_utm.transform(lon, lat)
    E, N = (ux - sx) / 1000, (uy - sy) / 1000
    t = np.array([datetime.strptime(r["acq_date"], "%Y-%m-%d").date() for r in rows])
    frp = np.array([float(r.get("frp") or 0) for r in rows])
    conf = np.array([str(r.get("confidence", "")) for r in rows])
    src = np.array([r["source"] for r in rows])
    # line distances
    (a0, b0), (a1, b1) = (to_utm.transform(*CRATER_ROW[0]), to_utm.transform(*CRATER_ROW[1]))
    x0, y0, x1, y1 = (a0 - sx) / 1000, (b0 - sy) / 1000, (a1 - sx) / 1000, (b1 - sy) / 1000

    def dist(px, py, x0, y0, x1, y1):
        vx, vy = x1 - x0, y1 - y0
        tt = np.clip(((px - x0) * vx + (py - y0) * vy) / (vx * vx + vy * vy), 0, 1)
        return np.hypot(px - (x0 + tt * vx), py - (y0 + tt * vy))
    d_row = dist(E, N, x0, y0, x1, y1)
    trace = None
    if (ML / "s1_npe_v3.npz").exists():
        pm = np.median(np.load(ML / "s1_npe_v3.npz")["ps"], 0)
        u = np.array([np.sin(np.deg2rad(pm[2])), np.cos(np.deg2rad(pm[2]))])
        trace = (pm[0] - u[0] * pm[3] / 2, pm[1] - u[1] * pm[3] / 2, pm[0] + u[0] * pm[3] / 2, pm[1] + u[1] * pm[3] / 2)
    d_dike = dist(E, N, *trace) if trace else np.full(len(E), np.nan)
    ev = []
    for s, e in ERUPTIONS:
        s, e = date.fromisoformat(s), date.fromisoformat(e)
        m = (t >= s) & (t <= e + timedelta(days=1))
        first = min(t[m]) if m.any() else None
        ev.append(dict(start=s.isoformat(), end=e.isoformat(), n=int(m.sum()), detected=bool(m.any()),
                       first_delay_days=(first - s).days if first else None, frp_max=round(float(frp[m].max()), 1) if m.any() else None,
                       median_dist_crater_row_km=round(float(np.median(d_row[m])), 2) if m.any() else None,
                       median_dist_dike_km=round(float(np.nanmedian(d_dike[m])), 2) if m.any() and trace else None))
    in_erupt = np.zeros(len(t), bool)
    for s, e in ERUPTIONS:
        in_erupt |= (t >= date.fromisoformat(s)) & (t <= date.fromisoformat(e) + timedelta(days=1))
    out = dict(n=len(rows), by_source={s_: int((src == s_).sum()) for s_ in SOURCES}, eruptions=ev,
               outside_eruptions=int((~in_erupt).sum()), outside_far_from_row=int((~in_erupt & (d_row > 3)).sum()),
               note="hotspots outside the eruption windows: cooling lava fields, Svartsengi power plant (industrial), false alarms")
    (ML / "firms_hotspots.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    for x in ev:
        print(x, flush=True)
    print("outside eruption windows:", out["outside_eruptions"], "of which > 3 km from the crater row:", out["outside_far_from_row"])
    # figure
    fig, ax = plt.subplots(1, 2, figsize=(16, 5.5), gridspec_kw=dict(width_ratios=[1, 1.5]))
    a = ax[0]
    mv = np.char.startswith(src.astype(str), "VIIRS")
    sc = a.scatter(E[mv], N[mv], c=[(x - START).days for x in t[mv]], s=4 + np.sqrt(frp[mv]), cmap="viridis", alpha=.5)
    a.plot([x0, x1], [y0, y1], "r-", lw=2, label="кратерний ряд Сундгнукур (геол. карти)")
    if trace:
        a.plot([trace[0], trace[2]], [trace[1], trace[3]], "k--", lw=1.5, label="слід дайки 10.11.2023 (NPE v3)")
    for k in ("Svartsengi", "Grindavik", "Sundhnukur"):
        px, py = to_utm.transform(*PLACES[k])
        a.plot((px - sx) / 1000, (py - sy) / 1000, "k^" if k == "Svartsengi" else "ko", ms=5)
        a.text((px - sx) / 1000 + .3, (py - sy) / 1000 - .5, {"Svartsengi": "Свартсенгі", "Grindavik": "Ґріндавік", "Sundhnukur": "Сундгнукур"}[k], fontsize=8)
    cb = plt.colorbar(sc, ax=a, fraction=.04); cb.set_label("днів від 01.11.2023", fontsize=8)
    a.set_aspect("equal"); a.set_xlim(-10, 12); a.set_ylim(-9, 11); a.grid(alpha=.3)
    a.legend(fontsize=7, loc="lower right"); a.set_xlabel("схід від Свартсенгі, км"); a.set_ylabel("північ, км")
    a.set_title(f"NASA FIRMS: {int(mv.sum())} точок VIIRS 375 м (на графіку праворуч — ще MODIS)", fontsize=10)
    a = ax[1]
    for s, e in ERUPTIONS:
        a.axvspan(date.fromisoformat(s), date.fromisoformat(e) + timedelta(days=1), color="orange", alpha=.25, lw=0)
    for s_, c in zip(SOURCES, ("tab:red", "tab:blue", "0.4")):
        m = src == s_
        a.scatter(t[m], frp[m], s=8, color=c, alpha=.6, label=s_.replace("_SP", ""))
    a.set_yscale("symlog", linthresh=10)
    a.set_ylabel("потужність випромінювання FRP, МВт"); a.legend(fontsize=8)
    a.set_title("теплові точки в часі; помаранчеве — відомі виверження", fontsize=10)
    fig.tight_layout()
    fig.savefig(ML / "firms_hotspots.png", dpi=100, facecolor="white")
    print("wrote", ML / "firms_hotspots.png")


if __name__ == "__main__":
    main()
