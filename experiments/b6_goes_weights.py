"""B6: GOES cloudiness as side information for the per-date noise (pre-registered).

sigma^2(date) = a + b * g(date), g = GOES band-13 BT std in the box at overpass time.
var(pair k) = 2a + b (g_k + g_{k+1})  -> 2 unknowns from 27 pairs (well-posed, unlike A).
Weighted per-pixel line fit; compare plain / A / B6 at MMX1 (GNSS) and against the morning map.

    uv run python scripts/b6_goes_weights.py
"""
from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path

import numpy as np
from rasterio.warp import transform as warp_transform

from hachure import gridding

ROOT = Path(__file__).resolve().parents[1]
P = ROOT / "data" / "processed"
BBOX = (-99.32, 19.17, -98.82, 19.58)
GNSS = -26.7
MMX1 = (-99.0684, 19.4317)


def goes_g(t):
    fs = sorted((ROOT / "data/raw/goes19").glob("bt13_*.npz"))
    ts = [datetime.strptime(f.stem[5:], "%Y%m%dT%H%M") for f in fs]
    k = int(np.argmin([abs((s - t).total_seconds()) for s in ts]))
    if abs((ts[k] - t).total_seconds()) > 3600:
        return np.nan
    g = np.load(fs[k])
    inbox = (g["lat"] >= BBOX[1]) & (g["lat"] <= BBOX[3]) & (g["lon"] >= BBOX[0]) & (g["lon"] <= BBOX[2])
    return float(np.nanstd(g["bt"][inbox]))


def fit(cum, t, w):
    w = w[:, None, None] * np.isfinite(cum)
    c = np.nan_to_num(cum)
    T = t[:, None, None]
    sw, st, sy = w.sum(0), (w * T).sum(0), (w * c).sum(0)
    stt, sty = (w * T * T).sum(0), (w * T * c).sum(0)
    with np.errstate(all="ignore"):
        v = (sw * sty - st * sy) / (sw * stt - st ** 2)
    v[sw == 0] = np.nan
    return v


def main():
    z = np.load(P / "mexico-city" / "nisar_series.npz")
    dates = [datetime.fromisoformat(str(d)) for d in z["dates"]]
    cum = z["cum"].astype("float64")
    t = np.array([(d - dates[0]).days for d in dates]) / 365.25
    pv = np.array([r["noise_std_cm"] for r in json.load(open(P / "mexico-city" / "atmo_stats.json"))["pairs"]]) ** 2
    g = np.array([goes_g(d) for d in dates])
    print(f"GOES BT std per date: {np.sum(np.isfinite(g))} of {len(g)} found; missing -> median")
    g = np.where(np.isfinite(g), g, np.nanmedian(g))
    X = np.c_[2 * np.ones(len(pv)), g[:-1] + g[1:]]
    (a, b), *_ = np.linalg.lstsq(X, pv, rcond=None)
    rng = np.random.default_rng(3)
    boot = []
    for _ in range(4000):
        i = rng.integers(0, len(pv), len(pv))
        boot.append(np.linalg.lstsq(X[i], pv[i], rcond=None)[0][1])
    lo, hi = np.percentile(boot, [2.5, 97.5])
    print(f"sigma^2 = a + b*g:  a = {a:.2f} cm^2, b = {b:.3f} cm^2/K  (95% CI of b {lo:.3f} .. {hi:.3f})")
    s2 = np.maximum(a + b * g, 0.05)
    for d, gg, s in zip(dates, g, np.sqrt(s2)):
        print(f"  {d:%Y-%m-%d}  BT std {gg:5.1f} K -> date noise {s:.2f} cm")

    wa = json.load(open(P / "mexico-city" / "weighted_rate.json"))["date_noise_cm"]
    s2_A = np.array([wa[f"{d:%Y-%m-%d}"] for d in dates]) ** 2
    rates = {"plain": fit(cum, t, np.ones(len(t))), "A (adjacent pairs)": fit(cum, t, 1 / s2_A),
             "B6 (GOES)": fit(cum, t, 1 / s2)}

    # morning reference on the same grid
    za = np.load(P / "mexico-city-asc" / "nisar_series.npz")
    ad = [datetime.fromisoformat(str(d)) for d in za["dates"]]
    ta = np.array([(d - ad[0]).days for d in ad]) / 365.25
    dst = (gridding.transform_from_centres(z["x"], z["y"]), f"EPSG:{int(z['epsg'])}", cum.shape[1:])
    acum = np.stack([gridding.to_grid(c, gridding.transform_from_centres(za["x"], za["y"]),
                                      f"EPSG:{int(za['epsg'])}", *dst) for c in za["cum"].astype("float64")])
    v_asc = fit(acum, ta, np.ones(len(ta)))

    xx, yy = np.meshgrid(z["x"], z["y"])
    sx, sy = warp_transform("EPSG:4326", f"EPSG:{int(z['epsg'])}", [MMX1[0]], [MMX1[1]])
    near = np.hypot(xx - sx[0], yy - sy[0]) <= 250
    common = np.isfinite(v_asc) & np.all([np.isfinite(v) for v in rates.values()], 0)
    out = dict(a=float(a), b=float(b), b_ci=[float(lo), float(hi)], asc_mmx1=float(np.nanmedian(v_asc[near]) * 100))
    print(f"\nmorning (asc) plain fit at MMX1: {out['asc_mmx1']:+.1f} cm/yr; GNSS {GNSS}")
    base = None
    for name, v in rates.items():
        m = np.nanmedian(v[near]) * 100
        d = (v_asc - v)[common] * 100
        sd = float(np.std(d - np.median(d)))
        base = base or sd
        out[name] = dict(mmx1=float(m), diff_gnss=float(m - GNSS), asc_minus_desc_std=sd,
                         change_vs_plain_pct=float((sd / base - 1) * 100))
        print(f"  {name:20s} MMX1 {m:+6.1f} (vs GNSS {m - GNSS:+5.1f})  | std(asc - desc) {sd:.2f} cm/yr "
              f"({(sd / base - 1) * 100:+.1f}% vs plain)")
    json.dump(out, open(P / "mexico-city" / "b6_goes_weights.json", "w"), indent=1)


if __name__ == "__main__":
    main()
