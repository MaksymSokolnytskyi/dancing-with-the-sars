"""H14 (pre-registered in notes/experiments.md): per-pair plane tie to daily GNSS, leave-one-out scoring.

    uv run python scripts/h14_pair_tie.py houston desc
"""
from __future__ import annotations

import csv
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
from rasterio.warp import transform as wt

sys.path[:0] = [str(Path(__file__).parent), str(Path(__file__).parents[1] / "scripts")]
from gnss_rates import read  # noqa: E402
from trust_audit import SITES  # noqa: E402

from hachure import nisar

ROOT = Path(__file__).resolve().parents[1]


def decyear(d):
    y0 = datetime(d.year, 1, 1)
    return d.year + (d - y0).total_seconds() / (datetime(d.year + 1, 1, 1) - y0).total_seconds()


def gnss_at(t, up, when, half=2 / 365.25):
    m = np.abs(t - decyear(when)) <= half
    return float(np.mean(up[m])) if m.any() else np.nan


def plane_fit(E, N, r):
    ok = np.isfinite(r)
    if ok.sum() >= 5:
        A = np.column_stack([np.ones(ok.sum()), E[ok], N[ok]])
        c = np.linalg.lstsq(A, r[ok], rcond=None)[0]
        return lambda e, n: c[0] + c[1] * e + c[2] * n
    if ok.sum() >= 2:
        m = float(np.mean(r[ok]))
        return lambda e, n: m + 0 * e
    return lambda e, n: 0 * e


def main():
    site, pas = sys.argv[1], sys.argv[2]
    s = SITES[site]
    bbox = s["bbox"]
    cache = ROOT / "data/raw/nisar_cache" / ("bbox_" + "_".join(f"{v:.3f}" for v in bbox) + f"_h{s['h']}")
    found = nisar.search((bbox[0] + bbox[2]) / 2, (bbox[1] + bbox[3]) / 2,
                         direction="ASCENDING" if pas == "asc" else "DESCENDING")
    chain = [g for g in nisar.sequential_chain(nisar.dominant_track(found)) if (cache / f"{g.name}.npz").exists()]
    pairs, kept = [], []
    for g in chain:
        try:
            pairs.append(nisar.process(nisar.RawWindow.load(cache / f"{g.name}.npz"), ref_bbox=s["ref"], coh_min=0.3,
                                       iono="plane"))
            kept.append(g)
        except ValueError:
            pass
    rate, coh, frac, x, y, disps = nisar.stack_rate(pairs)
    XX, YY = np.meshgrid(x, y)
    st = [r for r in csv.DictReader(open(ROOT / f"data/processed/gnss/{site}_rates.csv"))
          if r["window"] == pas and r["used"] == "True"]
    xs, ys = wt("EPSG:4326", f"EPSG:{pairs[0].epsg}", [float(r["lon"]) for r in st], [float(r["lat"]) for r in st])
    E, N = np.array(xs) / 1e3, np.array(ys) / 1e3
    E, N = E - E.mean(), N - N.mean()
    nis = np.full((len(pairs), len(st)), np.nan)
    gn = np.full((len(pairs), len(st)), np.nan)
    series = {r["station"]: read(ROOT / "data/raw/gnss" / site / f"{r['station']}.tenv3") for r in st}
    for j, (r, sx, sy) in enumerate(zip(st, xs, ys)):
        m = np.hypot(XX - sx, YY - sy) <= 250
        t, up = series[r["station"]][0], series[r["station"]][1]
        for k, (g, dd) in enumerate(zip(kept, disps)):
            if np.isfinite(dd[m]).sum() >= 3:
                nis[k, j] = np.nanmedian(dd[m])
            gn[k, j] = gnss_at(t, up, g.sec_date) - gnss_at(t, up, g.ref_date)
    days = np.array([p.days for p in pairs], float)
    grate = np.array([float(r["rate_cm_yr"]) for r in st])

    def chain_rate(d):
        ok = np.isfinite(d)
        return np.where(ok, d, 0).sum(0) / np.where(ok, days[:, None], 0).sum(0) * 365.25 * 100

    base = chain_rate(nis)
    loo = np.full(len(st), np.nan)
    for i in range(len(st)):
        corr = nis.copy()
        for k in range(len(pairs)):
            rr = nis[k] - gn[k]
            rr[i] = np.nan
            f = plane_fit(E, N, rr)
            corr[k] = nis[k] - f(E, N)
        loo[i] = chain_rate(corr)[i]
    for label, v in (("chain sum (no tie)", base), ("per-pair tie, LOO", loo)):
        d = v - grate
        ok = np.isfinite(d)
        d = d[ok] - np.median(d[ok])
        print(f"{site} {pas} {label:20s}: {ok.sum()} stations, RMS {np.sqrt(np.mean(d ** 2)):.2f} cm/yr, "
              f"median|d| {np.median(np.abs(d)):.2f}, max|d| {np.abs(d).max():.1f}")
    cover = np.isfinite(gn).sum(1)
    print(f"  stations with GNSS per pair: min {cover.min()}, median {int(np.median(cover))}, max {cover.max()}")
    for j, r in enumerate(st):
        print(f"   {r['station']}: GNSS {grate[j]:+.1f}  chain {base[j] - np.nanmedian(base - grate):+.1f}  "
              f"LOO tie {loo[j] - np.nanmedian(loo - grate):+.1f}")


if __name__ == "__main__":
    main()
