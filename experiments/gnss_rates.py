"""GNSS vertical rates for every station of a site, over each NISAR pass window (P-port).

    uv run python scripts/gnss_rates.py --site central-valley --window asc 2025-11-24 2026-10-14 --window desc 2025-11-30 2026-10-20
Writes data/processed/gnss/<site>_rates.csv. Rules fixed in the P-port pre-registration:
>= 70% daily coverage of the window, sigma = max(0.2 cm/yr, 3 x formal), equipment changes flagged.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def read(path):
    t, up, ant, lat, lon = [], [], [], None, None
    for line in open(path):
        f = line.split()
        if f[0] == "site" or len(f) < 20:
            continue
        t.append(float(f[2]))
        up.append(float(f[11]) + float(f[12]))
        ant.append(float(f[13]))
        lat, lon = float(f[20]), (float(f[21]) + 180) % 360 - 180
    return np.array(t), np.array(up), np.array(ant), lon, lat


def decyear(s):
    d = datetime.fromisoformat(s)
    y0 = datetime(d.year, 1, 1)
    return d.year + (d - y0).total_seconds() / (datetime(d.year + 1, 1, 1) - y0).total_seconds()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--site", required=True)
    ap.add_argument("--window", nargs=3, action="append", required=True, metavar=("NAME", "START", "END"))
    args = ap.parse_args()
    rows = []
    for p in sorted((ROOT / "data" / "raw" / "gnss" / args.site).glob("*.tenv3")):
        t, up, ant, lon, lat = read(p)
        for name, a, b in args.window:
            t0, t1 = decyear(a), decyear(b)
            m = (t >= t0) & (t <= t1)
            cover = m.sum() / ((t1 - t0) * 365.25)
            row = dict(station=p.stem, lon=lon, lat=lat, window=name, n=int(m.sum()), coverage=round(cover, 2),
                       rate_cm_yr=np.nan, sigma_cm_yr=np.nan, equipment_change=False, used=False)
            if m.sum() >= 30:
                A = np.c_[np.ones(m.sum()), t[m] - t0]
                coef, res, *_ = np.linalg.lstsq(A, up[m], rcond=None)
                r = up[m] - A @ coef
                cov = np.linalg.inv(A.T @ A) * r.var(ddof=2)
                row["rate_cm_yr"] = round(coef[1] * 100, 2)
                row["sigma_cm_yr"] = round(max(0.2, 3 * np.sqrt(cov[1, 1]) * 100), 2)
                row["equipment_change"] = bool(np.ptp(ant[m]) > 0)
                row["used"] = bool(cover >= 0.7)
            rows.append(row)
    out = ROOT / "data" / "processed" / "gnss"
    out.mkdir(parents=True, exist_ok=True)
    with open(out / f"{args.site}_rates.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    for name, *_ in args.window:
        u = [r for r in rows if r["window"] == name and r["used"]]
        print(f"{args.site} {name}: {len(u)} usable stations; rates (cm/yr): " +
              ", ".join(f"{r['station']} {r['rate_cm_yr']:+.1f}{'*' if r['equipment_change'] else ''}"
                        for r in sorted(u, key=lambda r: r["rate_cm_yr"])))
    print(out / f"{args.site}_rates.csv", "(* = equipment change inside the window)")


if __name__ == "__main__":
    main()
