"""B2: morning (ascending) vs evening (descending) pass over the whole city, dry vs wet season.

    uv run python scripts/asc_desc_check.py
Measure 1: per-date scatter around each pixel's yearly line, per pass and season.
Measure 2: dry-window and wet-window rates per pass; std of (asc - desc) difference maps.
"""
from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path

import numpy as np

from hachure import gridding

ROOT = Path(__file__).resolve().parents[1]
P = ROOT / "data" / "processed"


def load(name, dst=None):
    z = np.load(P / name / "nisar_series.npz")
    dates = [datetime.fromisoformat(str(d)) for d in z["dates"]]
    cum = z["cum"].astype("float64")
    t = gridding.transform_from_centres(z["x"], z["y"])
    crs = f"EPSG:{int(z['epsg'])}"
    if dst is not None:
        cum = np.stack([gridding.to_grid(c, t, crs, dst[0], dst[1], dst[2]) for c in cum])
    return dates, cum, (t, crs, cum.shape[1:])


def season(d):
    if d.month in (11, 12, 1, 2, 3, 4):
        return "dry"
    if (d.month == 5 and d.day >= 15) or d.month in (6, 7, 8, 9, 10):
        return "wet"
    return None


def linfit(t, y):
    """Per-pixel least-squares line through the dates in y (n, H, W); returns slope, residuals."""
    ok = np.isfinite(y)
    w = ok.astype(float)
    yc = np.nan_to_num(y)
    T = t[:, None, None]
    sw, st, sy = w.sum(0), (w * T).sum(0), (w * yc).sum(0)
    stt, sty = (w * T * T).sum(0), (w * T * yc).sum(0)
    with np.errstate(all="ignore"):
        v = (sw * sty - st * sy) / (sw * stt - st ** 2)
        a = (sy - v * st) / sw
    v[sw < 4] = np.nan
    return v, y - a - v * T


def main():
    d_dates, d_cum, grid = load("mexico-city")
    a_dates, a_cum, _ = load("mexico-city-asc", dst=grid)
    common = np.isfinite(d_cum[-1]) & np.isfinite(a_cum[-1])
    print(f"common valid pixels: {common.sum()} ({common.mean() * 100:.0f}% of grid)")

    out = {}
    per_date = {}
    for name, dates, cum in [("asc 06:05", a_dates, a_cum), ("desc 19:13", d_dates, d_cum)]:
        t = np.array([(d - dates[0]).days for d in dates]) / 365.25
        _, res = linfit(t, np.where(common, cum, np.nan))
        sd = []
        for k, d in enumerate(dates):
            r = res[k][common]
            sd.append(float(np.nanstd(r - np.nanmedian(r)) * 100))
        per_date[name] = (dates, sd)
        for s in ("dry", "wet"):
            v = [x for d, x in zip(dates, sd) if season(d) == s]
            out[f"{name} {s} scatter cm"] = float(np.median(v))
            print(f"  {name:10s} {s}: median per-date scatter {np.median(v):.2f} cm over {len(v)} dates")
    for s in ("dry", "wet"):
        ratio = out[f"desc 19:13 {s} scatter cm"] / out[f"asc 06:05 {s} scatter cm"]
        out[f"desc/asc {s}"] = ratio
        print(f"  desc/asc scatter ratio, {s}: {ratio:.2f}")

    maps = {}
    for s in ("dry", "wet"):
        rates = []
        for dates, cum in [(a_dates, a_cum), (d_dates, d_cum)]:
            idx = [k for k, d in enumerate(dates) if season(d) == s]
            t = np.array([(dates[k] - dates[idx[0]]).days for k in idx]) / 365.25
            v, _ = linfit(t, cum[idx])
            rates.append(v)
        diff = np.where(common, rates[0] - rates[1], np.nan) * 100
        diff -= np.nanmedian(diff)
        maps[s] = diff
        out[f"asc-desc {s} map std cm/yr"] = float(np.nanstd(diff))
        print(f"  {s} window: asc - desc rate map std {np.nanstd(diff):.1f} cm/yr "
              f"(rates from {len(idx)} desc dates)")
    out["wet/dry map std ratio"] = out["asc-desc wet map std cm/yr"] / out["asc-desc dry map std cm/yr"]
    print(f"  wet/dry ratio of difference-map std: {out['wet/dry map std ratio']:.2f}")
    json.dump(out, open(P / "mexico-city" / "asc_desc_check.json", "w"), indent=1)
    figure(per_date, maps, out)


def figure(per_date, maps, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(16, 5))
    ax = fig.add_subplot(1, 3, 1)
    for name, c in [("asc 06:05", "#3b6ea8"), ("desc 19:13", "#a8322a")]:
        dates, sd = per_date[name]
        ax.plot(dates, sd, "o-", color=c, ms=4, label=f"{name} local time")
    ax.axvspan(datetime(2026, 5, 15), datetime(2026, 10, 31), color="0.9", zorder=0)
    ax.text(datetime(2026, 7, 1), ax.get_ylim()[1] * 0.95, "wet season", ha="center", fontsize=9)
    ax.set_ylabel("scatter around the yearly line, cm")
    ax.set_title("Per-date atmosphere: morning vs evening pass")
    ax.legend()
    ax.tick_params(axis="x", rotation=45)
    lim = np.nanpercentile(np.abs(maps["wet"]), 98)
    for k, s in enumerate(("dry", "wet")):
        ax = fig.add_subplot(1, 3, 2 + k)
        h = ax.imshow(maps[s], cmap="RdBu", vmin=-lim, vmax=lim, interpolation="nearest")
        fig.colorbar(h, ax=ax, shrink=0.8, label="cm/yr")
        ax.set_title(f"{s} season: morning rate − evening rate\nstd {out[f'asc-desc {s} map std cm/yr']:.1f} cm/yr")
        ax.axis("off")
    fig.tight_layout()
    fig.savefig(P / "mexico-city" / "asc_desc_check.png", dpi=85)
    print(P / "mexico-city" / "asc_desc_check.png")


if __name__ == "__main__":
    main()
