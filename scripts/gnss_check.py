"""Validate the NISAR time series against a GNSS station (NGL tenv3, IGS20 frame).

    uv run python scripts/gnss_check.py --station MMX1
NISAR is relative to the Guadalupe reference rock, GNSS is absolute; both are shifted to
0 at the first NISAR date, so only the motion over the year is compared.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
from rasterio.warp import transform as warp_transform

ROOT = Path(__file__).resolve().parents[1]


def read_tenv3(path):
    t, up, lon, lat = [], [], None, None
    for line in open(path):
        f = line.split()
        if f[0] == "site" or len(f) < 13:
            continue
        yr = float(f[2])
        y0 = int(yr)
        d0 = datetime(y0, 1, 1)
        t.append(d0 + (datetime(y0 + 1, 1, 1) - d0) * (yr - y0))
        up.append(float(f[11]) + float(f[12]))
        lat, lon = float(f[-3]), float(f[-2])
    return np.array(t), np.array(up), lon, lat


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--station", default="MMX1")
    ap.add_argument("--name", default="mexico-city")
    ap.add_argument("--radius-m", type=float, default=250, help="NISAR pixels averaged around the station")
    args = ap.parse_args()

    t, up, lon, lat = read_tenv3(ROOT / "data" / "raw" / "gnss" / f"{args.station}.tenv3")
    z = np.load(ROOT / "data" / "processed" / args.name / "nisar_series.npz")
    dates = np.array([datetime.fromisoformat(str(d)) for d in z["dates"]])
    sx, sy = warp_transform("EPSG:4326", f"EPSG:{int(z['epsg'])}", [lon], [lat])
    xx, yy = np.meshgrid(z["x"], z["y"])
    near = np.hypot(xx - sx[0], yy - sy[0]) <= args.radius_m
    cum = z["cum"][:, near]
    nis = np.nanmedian(cum, axis=1)
    print(f"{args.station} at {lat:.4f},{lon:.4f}: {near.sum()} NISAR pixels within {args.radius_m:.0f} m, "
          f"{np.isfinite(cum[-1]).sum()} valid")

    # GNSS: 7-day running median at each NISAR date, zeroed at the first date
    g_at = []
    for d in dates:
        m = np.abs(np.array([(x - d).total_seconds() for x in t])) <= 3.5 * 86400
        g_at.append(np.median(up[m]) if m.sum() >= 3 else np.nan)
    g_at = np.array(g_at) - (g_at[0] if np.isfinite(g_at[0]) else np.nan)
    days = np.array([(d - dates[0]).days for d in dates]) / 365.25
    ok = np.isfinite(g_at) & np.isfinite(nis)
    rg = np.polyfit(days[ok], g_at[ok], 1)[0] * 100
    rn = np.polyfit(days[ok], nis[ok], 1)[0] * 100
    diff = (nis - g_at)[ok] * 100
    print(f"rate over the year: GNSS {rg:+.1f} cm/yr, NISAR {rn:+.1f} cm/yr")
    print(f"NISAR - GNSS at the {ok.sum()} dates: mean {diff.mean():+.2f} cm, std {diff.std():.2f} cm, "
          f"max |diff| {np.abs(diff).max():.2f} cm")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    w = (t >= dates[0] - timedelta(days=30)) & (t <= dates[-1] + timedelta(days=30))
    base = np.interp(dates[0].timestamp(), [x.timestamp() for x in t[w]], up[w])
    fig, ax = plt.subplots(figsize=(9, 4.5))
    ax.plot(t[w], (up[w] - base) * 100, ".", color="0.6", ms=3, label=f"GNSS {args.station}, daily")
    ax.plot(dates, g_at * 100, "-", color="0.2", lw=1.5, label="GNSS, 7-day median at NISAR dates")
    ax.plot(dates, nis * 100, "o-", color="#a8322a", ms=4, label="NISAR, 27-pair chain")
    ax.set_ylabel("vertical motion since 25 Oct 2025, cm")
    ax.set_title(f"{args.station}: GNSS {rg:+.1f} cm/yr vs NISAR {rn:+.1f} cm/yr  "
                 f"(diff std {diff.std():.1f} cm)")
    ax.grid(alpha=0.3)
    ax.legend()
    fig.tight_layout()
    out = ROOT / "data" / "processed" / args.name / f"gnss_{args.station}.png"
    fig.savefig(out, dpi=100)
    print(out)


if __name__ == "__main__":
    main()
