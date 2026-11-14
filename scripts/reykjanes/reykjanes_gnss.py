"""Check NISAR at Reykjanes against GNSS (Nevada Geodetic Laboratory daily IGS20 series).

1. NISAR LOS time series at the station pixels vs GNSS east/north/up projected on the
   same line of sight; both relative to the far field (GNSS minus station NYLA, NISAR
   referenced to pixels > 14 km from Svartsengi).
2. Mogi model (output/reykjanes/mogi_fit_offset.npz) predicted E/N/U velocities at the
   stations vs GNSS velocities over the same period.

Usage: python scripts/reykjanes/reykjanes_gnss.py
"""
from datetime import datetime, timedelta
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np

from reykjanes_analysis import (MIN_PAIR_COH, OUT, Geo, chain, los_cm, pair_coherence, pairs_by_track)

ROOT = Path(__file__).resolve().parents[2]
GNSS = ROOT / "data" / "gnss"
STATIONS = ["SENG", "AUSV", "GEVK", "VOGC"]
REF = "NYLA"
TRACKS = ("159A", "082D")
HALF = 2  # 5 x 5 pixels (400 m) around a station


def read_tenv3(name):
    t, e, n, u, ll = [], [], [], [], None
    for line in open(GNSS / f"{name}.tenv3"):
        p = line.split()
        if p[0] == "site":
            continue
        t.append(datetime(2000, 1, 1) + timedelta(days=float(p[3]) - 51544))  # MJD 51544 = 2000-01-01
        e.append(float(p[7]) + float(p[8]))
        n.append(float(p[9]) + float(p[10]))
        u.append(float(p[11]) + float(p[12]))
        ll = (float(p[21]) - (360 if float(p[21]) > 180 else 0), float(p[20]))
    return np.array(t), np.c_[e, n, u] * 100, ll  # cm


def relative(name, ref):
    t, x, ll = read_tenv3(name)
    tr, xr, _ = read_tenv3(ref)
    common, ia, ib = np.intersect1d(t, tr, return_indices=True)
    return common, x[ia] - xr[ib], ll


def at_dates(t, x, dates, half_days=3):
    out = []
    for d in dates:
        m = np.abs(np.array([(s - d).days for s in t])) <= half_days
        out.append(np.median(x[m], 0) if m.any() else np.full(3, np.nan))
    out = np.array(out)
    return out - out[0]


def nisar_station_series(geo, track, stations):
    prs = chain(pairs_by_track()[track])
    dates = [prs[0][0]] + [t1 for _, t1, _ in prs]
    steps = {s: [] for s in stations}
    good, dts, los = [], [], None
    for t0, t1, p in prs:
        ok = pair_coherence(p, geo) >= MIN_PAIR_COH
        d, g = los_cm(p, geo)
        ok = ok and d is not None
        los = g if g is not None else los
        for s, (i, j) in stations.items():
            v = np.nanmedian(d[i - HALF:i + HALF + 1, j - HALF:j + HALF + 1]) if d is not None else np.nan
            steps[s].append(v)
        good.append(ok)
        dts.append((t1 - t0).days)
    good, dts = np.array(good), np.array(dts)
    series, vec = {}, {}
    for s, (i, j) in stations.items():
        v = np.array(steps[s])
        g = good & np.isfinite(v)
        rate = np.sum(v[g]) / np.sum(dts[g])
        series[s] = np.r_[0, np.cumsum(np.where(g, v, rate * dts))]
        e = np.nanmedian(los["losUnitVectorX"][i - HALF:i + HALF + 1, j - HALF:j + HALF + 1])
        n = np.nanmedian(los["losUnitVectorY"][i - HALF:i + HALF + 1, j - HALF:j + HALF + 1])
        vec[s] = np.array([e, n, np.sqrt(1 - e ** 2 - n ** 2)])
    return dates, series, vec, good


def mogi_enu(fit, lon, lat):
    """Predicted E, N, U velocity (cm/yr) of the fitted Mogi source at a lon/lat."""
    k = 3 / (4 * np.pi) * fit["lin"][0] * 1e5
    de = (lon - float(fit["lon"])) * 111.32 * np.cos(np.deg2rad(lat))
    dn = (lat - float(fit["lat"])) * 111.2
    d = float(fit["theta"][2])
    R3 = (de ** 2 + dn ** 2 + d ** 2) ** 1.5
    return np.array([k * de / R3, k * dn / R3, k * d / R3])


def rate(t, x, t0, t1):
    m = np.array([(s >= t0) & (s <= t1) for s in t])
    yr = np.array([(s - t0).days / 365.25 for s in t[m]])
    return np.array([np.polyfit(yr, x[m, c], 1)[0] for c in range(3)])


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    geo = Geo(pairs_by_track()["159A"][0][2])
    gn = {s: relative(s, REF) for s in STATIONS}
    pix = {}
    for s, (_, _, ll) in gn.items():
        i, j = geo.ij(*ll)
        if 0 <= i < geo.y.size and 0 <= j < geo.x.size:
            pix[s] = (i, j)
    print("stations on the NISAR grid:", {s: f"{gn[s][2][1]:.4f}N {gn[s][2][0]:.4f}E, {geo.dist_km[pix[s]]:.1f} km from Svartsengi"
                                          for s in pix})
    fig, axes = plt.subplots(len(pix), len(TRACKS), figsize=(12, 2.6 * len(pix)), sharex="col", squeeze=False)
    rows = []
    for c, tr in enumerate(TRACKS):
        dates, ns, vec, good = nisar_station_series(geo, tr, pix)
        for r, s in enumerate(pix):
            t, x, _ = gn[s]
            g = at_dates(t, x, dates)
            g_los = g @ vec[s]
            diff = ns[s] - g_los
            ok = np.isfinite(diff)
            span = (dates[-1] - dates[0]).days / 365.25
            rows.append((tr, s, ns[s][-1], g_los[ok][-1] if ok.any() else np.nan,
                         np.sqrt(np.nanmean((diff - np.nanmean(diff)) ** 2)), span))
            ax = axes[r, c]
            win = np.array([(d >= dates[0] - timedelta(days=10)) & (d <= dates[-1] + timedelta(days=10)) for d in t])
            gl = x[win] @ vec[s]
            ax.plot(t[win], gl - np.interp(mdates.date2num(dates[0]), mdates.date2num(t[win]), gl),
                    color="0.75", lw=0.6, label="GNSS daily, projected on LOS")
            ax.plot(dates, g_los, "k.-", ms=3, lw=0.8, label="GNSS at SAR dates (±3 d)")
            ax.plot(dates, ns[s], "o-", ms=3, lw=1.2, color="tab:red", label="NISAR")
            bad = [i + 1 for i, q in enumerate(good) if not q]
            ax.plot([dates[i] for i in bad], ns[s][bad], "x", color="tab:red", ms=7)
            ax.set_xlim(dates[0] - timedelta(days=10), dates[-1] + timedelta(days=10))
            ax.set_title(f"{s} · track {tr}: NISAR {ns[s][-1]:+.1f}, GNSS {rows[-1][3]:+.1f} cm; rms diff {rows[-1][4]:.1f} cm",
                         fontsize=8)
            ax.grid(alpha=0.3)
            ax.set_ylabel("LOS cm", fontsize=8)
            if r == 0 and c == 0:
                ax.legend(fontsize=7, frameon=False)
            ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %y"))
    fig.suptitle(f"NISAR vs GNSS (relative to {REF}) at station pixels; x = low-coherence pair bridged", fontsize=10)
    fig.tight_layout()
    fig.savefig(OUT / "gnss_vs_nisar.png", dpi=120, facecolor="white")
    plt.close(fig)
    print("\ntrack station  NISAR_total  GNSS_total  rms_diff  (cm LOS)")
    for tr, s, a, b, e, sp in rows:
        print(f"{tr:5s} {s:6s} {a:+10.1f} {b:+11.1f} {e:9.1f}   over {sp:.2f} yr")

    # Mogi predictions vs GNSS velocities, Nov 2025 - Oct 2026
    fit = np.load(OUT / "mogi_fit_offset.npz")
    t0, t1 = datetime(2025, 11, 3), datetime(2026, 10, 17)
    print(f"\nMogi (depth {float(fit['theta'][2]):.1f} km, dV {fit['lin'][0] * 1e3:.1f}e6 m3/yr) vs GNSS rate "
          f"{t0:%Y-%m-%d} - {t1:%Y-%m-%d}, cm/yr relative to {REF}")
    print("station   GNSS  E     N     U   |  Mogi  E     N     U")
    out = []
    for s in STATIONS:
        t, x, ll = gn[s]
        g = rate(t, x, t0, t1)
        m = mogi_enu(fit, *ll) - mogi_enu(fit, *read_tenv3(REF)[2])
        out.append((s, g, m))
        print(f"{s:6s} {g[0]:+6.1f}{g[1]:+6.1f}{g[2]:+6.1f}   | {m[0]:+6.1f}{m[1]:+6.1f}{m[2]:+6.1f}")
    plot_mogi_vs_gnss(out, fit)


def plot_mogi_vs_gnss(out, fit):
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.6))
    for c, (ax, lab) in enumerate(zip(axes, ("east", "north", "up"))):
        g = [o[1][c] for o in out]
        m = [o[2][c] for o in out]
        ax.scatter(g, m, color="tab:red")
        for o in out:
            ax.annotate(o[0], (o[1][c], o[2][c]), xytext=(4, 3), textcoords="offset points", fontsize=8)
        lim = max(map(abs, g + m)) * 1.2 + 0.5
        ax.plot([-lim, lim], [-lim, lim], color="0.6", lw=0.8)
        ax.set_xlim(-lim, lim), ax.set_ylim(-lim, lim)
        ax.set_xlabel(f"GNSS {lab}, cm/yr"), ax.set_ylabel(f"Mogi {lab}, cm/yr")
        ax.grid(alpha=0.3)
    fig.suptitle(f"Mogi source from NISAR (depth {float(fit['theta'][2]):.1f} km) vs GNSS velocities, "
                 f"Nov 2025 - Oct 2026, relative to {REF}", fontsize=10)
    fig.tight_layout()
    fig.savefig(OUT / "mogi_vs_gnss.png", dpi=120, facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    main()
