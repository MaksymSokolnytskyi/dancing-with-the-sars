"""B4: is the evening (01 UTC) troposphere over Mexico City more variable than the morning (12 UTC)?

Uses GNSS zenith delay at MMX1 (NGL, 5-min) and reanalysis water vapour (Open-Meteo, hourly,
5x5 grid), then correlates with the radar residual noise per pair.

    uv run python scripts/tropo_external.py
"""
from __future__ import annotations

import gzip
import json
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
WET_MM_PER_KGM2 = 6.5          # zenith wet delay per mm of precipitable water (rule of thumb)


def load_ztd():
    t, ztd, ge, gn = [], [], [], []
    for f in sorted((ROOT / "data/raw/gnss/trop/x").glob("MMX1.*.trop.gz")):
        on = False
        for line in gzip.open(f, "rt"):
            if line.startswith("+TROP/SOLUTION"):
                on = True
                continue
            if line.startswith("-TROP/SOLUTION"):
                break
            if not on or line.startswith("*"):
                continue
            p = line.split()
            yy, doy, sec = p[1].split(":")
            t.append(datetime(2000 + int(yy), 1, 1) + timedelta(days=int(doy) - 1, seconds=int(sec)))
            ztd.append(float(p[2]))
            # NGL README: the two gradient columns are swapped in the files
            gn.append(float(p[5]))
            ge.append(float(p[7]))
    return np.array(t), np.array(ztd), np.hypot(np.array(ge), np.array(gn))


def at_hour(t, v, hour, minute=0, win_min=15):
    """Daily value nearest to hh:mm UTC (median within +-win), keyed by date."""
    out = {}
    secs = np.array([(x - x.replace(hour=0, minute=0, second=0)).total_seconds() for x in t])
    target = hour * 3600 + minute * 60
    m = np.abs(secs - target) <= win_min * 60
    days = np.array([x.date() for x in t])
    for d in np.unique(days[m]):
        sel = m & (days == d)
        out[d] = np.median(v[sel])
    return out


def main():
    t, ztd, grad = load_ztd()
    print(f"GNSS MMX1: {len(t)} 5-min epochs, {t[0]:%Y-%m-%d} -> {t[-1]:%Y-%m-%d}")
    res = {}
    print("\n1. GNSS at MMX1, all days of the year:")
    print("   hour UTC | day-to-day ZTD change std (mm) | 12-day ZTD change std (mm) | median gradient (mm) | 1-hour ZTD std (mm)")
    curve = []
    for h in range(24):
        a = at_hour(t, ztd, h)
        g = at_hour(t, grad, h)
        ds = sorted(a)
        d1 = [a[d] - a[d - timedelta(days=1)] for d in ds if d - timedelta(days=1) in a]
        d12 = [a[d] - a[d - timedelta(days=12)] for d in ds if d - timedelta(days=12) in a]
        # within-hour variability (turbulence-ish): std of 5-min ZTD within +-30 min, median over days
        secs = np.array([(x - x.replace(hour=0, minute=0, second=0)).total_seconds() for x in t])
        m = np.abs(secs - h * 3600) <= 1800
        days = np.array([x.date() for x in t])[m]
        z = ztd[m]
        hs = [np.std(z[days == d]) for d in np.unique(days) if np.sum(days == d) >= 8]
        row = (h, np.std(d1), np.std(d12), np.median(list(g.values())), np.median(hs))
        curve.append(row)
        mark = "  <- morning pass" if h == 12 else ("  <- evening pass" if h == 1 else "")
        print(f"   {h:02d}       | {row[1]:6.1f} | {row[2]:6.1f} | {row[3]:5.2f} | {row[4]:4.2f}{mark}")
    c = {r[0]: r for r in curve}
    res["gnss_12day_change_ratio_evening_morning"] = c[1][2] / c[12][2]
    res["gnss_gradient_ratio"] = c[1][3] / c[12][3]
    res["gnss_1h_std_ratio"] = c[1][4] / c[12][4]
    print(f"   evening/morning: 12-day change {res['gnss_12day_change_ratio_evening_morning']:.2f}, "
          f"gradient {res['gnss_gradient_ratio']:.2f}, 1-hour std {res['gnss_1h_std_ratio']:.2f}")

    # reanalysis grid
    js = json.load(open(ROOT / "data/raw/weather/openmeteo_grid.json"))
    times = np.array([datetime.fromisoformat(x) for x in js[0]["hourly"]["time"]])
    tcwv = np.array([p["hourly"]["total_column_integrated_water_vapour"] for p in js], float)  # (25, T)
    blh = np.array([p["hourly"]["boundary_layer_height"] for p in js], float)
    pr = np.array([p["hourly"]["precipitation"] for p in js], float)
    elev = np.array([p["elevation"] for p in js])
    low = elev < 2600                                    # valley points, like the radar's coherent area
    hours = np.array([x.hour for x in times])
    print(f"\n2. Reanalysis, {len(js)} grid points ({low.sum()} in the valley), hourly:")
    for h in (1, 12):
        sel = hours == h
        spread = np.nanstd(tcwv[low][:, sel], axis=0) * WET_MM_PER_KGM2
        dd = tcwv[low][:, sel][:, 12:] - tcwv[low][:, sel][:, :-12]
        dspread = np.nanstd(dd - np.nanmean(dd, 0), axis=0) * WET_MM_PER_KGM2
        print(f"   {h:02d} UTC: valley spread of wet delay {np.nanmedian(spread):.1f} mm; spread of 12-day change "
              f"{np.nanmedian(dspread):.1f} mm; boundary layer {np.nanmedian(blh[low][:, sel]):.0f} m; "
              f"rain in hour (any point) {np.mean(np.nanmax(pr[:, sel], 0) > 0.2) * 100:.0f}% of days")
        res[f"era_spread12_{h:02d}"] = float(np.nanmedian(dspread))
        res[f"blh_{h:02d}"] = float(np.nanmedian(blh[low][:, sel]))
    res["era_spread_ratio"] = res["era_spread12_01"] / res["era_spread12_12"]
    print(f"   evening/morning ratio of 12-day spread: {res['era_spread_ratio']:.2f}")

    # 3. per pair: radar residual noise vs external variability
    iv = json.load(open(ROOT / "data/processed/mexico-city/iono_vs_tropo.json"))
    print("\n3. Per pair, radar residual noise vs external change between the pair's two dates:")
    for name, hh, mm in [("asc 06:05", 12, 5), ("desc 19:13", 1, 13)]:
        rows = iv[name]["rows"]
        zd = at_hour(t, ztd, hh, mm)
        x1, x2, y = [], [], []
        for r in rows:
            d0 = datetime.fromisoformat(r["pair"]).date() + (timedelta(days=1) if hh == 1 else timedelta(0))
            d1 = d0 + timedelta(days=12)
            k0 = np.argmin(np.abs(times - datetime.combine(d0, datetime.min.time()).replace(hour=hh)))
            k1 = np.argmin(np.abs(times - datetime.combine(d1, datetime.min.time()).replace(hour=hh)))
            dfield = (tcwv[low][:, k1] - tcwv[low][:, k0]) * WET_MM_PER_KGM2
            if d0 in zd and d1 in zd:
                x1.append(abs(zd[d1] - zd[d0]))
                x2.append(np.nanstd(dfield))
                y.append(r["res_cm"])
        r1, p1 = spearmanr(x1, y)
        r2, p2 = spearmanr(x2, y)
        res[f"spearman_gnss_{name}"] = float(r1)
        res[f"spearman_era_{name}"] = float(r2)
        print(f"   {name}: n={len(y)}  vs |GNSS ZTD change| rho={r1:+.2f} (p={p1:.2f});  "
              f"vs reanalysis spread of change rho={r2:+.2f} (p={p2:.2f})")
    json.dump(res, open(ROOT / "data/processed/mexico-city/tropo_external.json", "w"), indent=1)


if __name__ == "__main__":
    main()
