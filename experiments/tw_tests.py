"""TW1-TW3: forest HV at dawn and dusk against the afternoon vapour pressure deficit (criteria in notes/experiments.md,
section TW). Uses the cached GCOV windows of sections W (Wisconsin, band B) and B (Bakota, band A); run w_fetch.py,
b_geometry.py and b_fetch.py first. Fetches hourly VPD and rain from the Open-Meteo archive (cached as JSON).

uv run python experiments/tw_tests.py     -> data/processed/tree_water/{tw_results.json, tw_tests.png, weather_*.json}
"""
import json
import sys
from datetime import datetime, timedelta
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import requests  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src"), str(ROOT / "experiments")]
import k_fetch  # noqa: E402  (ASCII CA bundle for GDAL; must come before rasterio)
import rasterio  # noqa: E402
from rasterio.warp import Resampling  # noqa: E402

import b_geometry  # noqa: E402
import w_fetch  # noqa: E402
from hachure import gcov  # noqa: E402
from k_tests import MIN_PX, db, on_grid  # noqa: E402
from w_tests import classes  # noqa: E402

OUT = ROOT / "data/processed/tree_water"
SEASON = ("2026-06-15", "2026-09-30")       # local solar dates
RAIN_MM = 1.0
N_PERM, N_BOOT = 10_000, 5_000
SITES = {"wisconsin": w_fetch.BBOX, "bakota": b_geometry.BBOX}


def center(bbox):
    return (bbox[0] + bbox[2]) / 2, (bbox[1] + bbox[3]) / 2


def weather(site):
    f = OUT / f"weather_{site}.json"
    if not f.exists():
        lon, lat = center(SITES[site])
        r = requests.get("https://archive-api.open-meteo.com/v1/archive",
                         params={"latitude": lat, "longitude": lon, "start_date": "2026-06-10", "end_date": "2026-10-03",
                                 "hourly": "vapour_pressure_deficit,precipitation", "timezone": "UTC"},
                         headers=k_fetch.UA, timeout=120)
        r.raise_for_status()
        f.write_text(json.dumps(r.json()["hourly"]), encoding="utf-8")
    h = json.loads(f.read_text(encoding="utf-8"))
    vpd = {datetime.fromisoformat(t): v for t, v in zip(h["time"], h["vapour_pressure_deficit"])}
    rain = {datetime.fromisoformat(t): v for t, v in zip(h["time"], h["precipitation"])}
    return vpd, rain


def utc_time(name):
    return datetime.strptime(name.split("_")[11], "%Y%m%dT%H%M%S")


def afternoon_vpd(vpd, local_day, lon):
    """Mean hourly VPD over 12:00-18:00 local solar time of `local_day` (hours stamped in UTC)."""
    off = timedelta(hours=lon / 15)
    t0 = datetime.combine(local_day, datetime.min.time()) + timedelta(hours=12) - off
    v = [vpd.get(t0.replace(minute=0, second=0, microsecond=0) + timedelta(hours=i)) for i in range(1, 7)]
    v = [x for x in v if x is not None]
    return float(np.mean(v)) if len(v) >= 4 else None


def pass_weather(t_utc, lon, vpd, rain):
    loc = t_utc + timedelta(hours=lon / 15)
    day = loc.date() if loc.hour >= 18 else loc.date() - timedelta(days=1)
    h = t_utc.replace(minute=0, second=0)
    r24 = [rain.get(h - timedelta(hours=i)) for i in range(24)]
    r48 = [rain.get(h - timedelta(hours=i)) for i in range(48)]
    return {"local": loc.isoformat(timespec="minutes"), "local_date": loc.date().isoformat(),
            "vpd": afternoon_vpd(vpd, day, lon), "vpd_after": afternoon_vpd(vpd, day + timedelta(days=1), lon),
            # hours absent from the Open-Meteo file (or NaN) add nothing to the sums below, so a gap could make a
            # wet pass look dry; count them so the output shows whether that happened
            "rain24_missing": sum(x is None or np.isnan(x) for x in r24),
            "rain48_missing": sum(x is None or np.isnan(x) for x in r48),
            "rain24": float(np.nansum([x for x in r24 if x is not None])),
            "rain48": float(np.nansum([x for x in r48 if x is not None]))}


def rows_wisconsin(vpd, rain):
    G = dict(np.load(w_fetch.OUT / "grid.npz"))
    tf, shape = rasterio.Affine(*G["transform"]), G["wc_water"].shape
    C = classes(G)
    lon = center(w_fetch.BBOX)[0]
    rows = []
    for p in sorted(gcov.cache_dir(w_fetch.BBOX, "B").glob("NISAR_*.npz")):
        t = utc_time(p.stem)
        r = {"site": "wisconsin", "name": p.stem, "track": int(p.stem.split("_")[5]), "dir": p.stem.split("_")[6],
             "release": p.stem.split("_")[13], **pass_weather(t, lon, vpd, rain)}
        if not SEASON[0] <= r["local_date"] <= SEASON[1]:
            continue
        L = on_grid(gcov.Window.load(p), tf, shape, epsg=w_fetch.EPSG)
        for c, key in (("decid", "forest"), ("conif", "conif")):
            for pol in ("hv", "hh"):
                v = L[pol][C[c]]
                v = v[np.isfinite(v)]
                r[f"{key}_{pol}"] = float(np.median(db(v))) if v.size >= MIN_PX else None
                r[f"n_{key}"] = int(v.size)
        rows.append(r)
    return rows


def rows_bakota(vpd, rain):
    G = dict(np.load(b_geometry.OUT / "geometry.npz"))
    dem_tf = rasterio.Affine(*G["transform"])
    gy, gx = np.gradient(G["z"], dem_tf.a)
    slope = np.degrees(np.arctan(np.hypot(gx, gy)))
    lon = center(b_geometry.BBOX)[0]
    rows = []
    for p in sorted(gcov.cache_dir(b_geometry.BBOX, "A").glob("NISAR_*.npz")):
        t = utc_time(p.stem)
        tr = int(p.stem.split("_")[5])
        r = {"site": "bakota", "name": p.stem, "track": tr, "dir": p.stem.split("_")[6],
             "release": p.stem.split("_")[13], **pass_weather(t, lon, vpd, rain)}
        if not SEASON[0] <= r["local_date"] <= SEASON[1]:
            continue
        w = gcov.Window.load(p)
        L = dict(np.load(b_geometry.OUT / f"gridA_{p.stem}.npz"))
        dx, dy = w.x[1] - w.x[0], w.y[1] - w.y[0]
        atf = rasterio.Affine(dx, 0, w.x[0] - dx / 2, 0, dy, w.y[0] - dy / 2)
        sl = k_fetch.warp(slope, dem_tf, "EPSG:32635", atf, w.hh.shape, how=Resampling.nearest, epsg=w.epsg)
        valid = (w.mask > 0) & (w.mask < 255) & (w.hv > 0) & (w.hh > 0)
        forest = (valid & (L["wc_tree"] >= 0.9) & (sl < 10)
                  & ~(L[f"geo_shadow_{tr}"] > 0.5) & ~(L[f"geo_layover_{tr}"] > 0.5))
        r["n_forest"] = int(forest.sum())
        for pol, a in (("hv", w.hv), ("hh", w.hh)):
            r[f"forest_{pol}"] = float(np.median(db(a[forest]))) if forest.sum() >= MIN_PX else None
        rows.append(r)
    return rows


def anomalies(rows, key, xkey="vpd", centre_only=False):
    """Per (site, track): detrend y = rows[key] and x = rows[xkey] with a line in time (or just centre them)."""
    out = []
    groups = {}
    for r in rows:
        if r.get(key) is not None and r.get(xkey) is not None:
            groups.setdefault((r["site"], r["track"]), []).append(r)
    for g, rr in sorted(groups.items()):
        if len(rr) < 3:
            continue
        d = np.array([(datetime.fromisoformat(r["local"]) - datetime(2026, 6, 15)).total_seconds() / 86400 for r in rr])
        res = []
        for k in (key, xkey):
            v = np.array([r[k] for r in rr])
            res.append(v - v.mean() if centre_only else v - np.polyval(np.polyfit(d, v, 1), d))
        out += [{"site": g[0], "track": g[1], "date": r["local_date"], "y": float(a), "x": float(b)}
                for r, a, b in zip(rr, *res)]
    return out


def rho_perm(an, rng):
    if len(an) < 3:
        return None, None
    y = np.array([a["y"] for a in an])
    x = np.array([a["x"] for a in an])
    rho = float(spearmanr(x, y)[0])
    idx = {}
    for i, a in enumerate(an):
        idx.setdefault((a["site"], a["track"]), []).append(i)
    xp = x.copy()
    hits = 0
    for _ in range(N_PERM):
        for ii in idx.values():
            xp[ii] = rng.permutation(x[ii])
        hits += spearmanr(xp, y)[0] <= rho
    return rho, (hits + 1) / (N_PERM + 1)


def slope_ci(an, rng):
    if len(an) < 5:
        return None
    x = np.array([a["x"] for a in an])
    y = np.array([a["y"] for a in an])
    b = [np.polyfit(x[i], y[i], 1)[0] for i in (rng.integers(0, len(x), len(x)) for _ in range(N_BOOT))]
    return [float(np.polyfit(x, y, 1)[0]), float(np.percentile(b, 5)), float(np.percentile(b, 95))]


def analyse(rows, rng, key="forest_hv", rain_key="rain24", centre_only=False, xkey="vpd"):
    dry = [r for r in rows if r[rain_key] < RAIN_MM]
    res = {}
    for when, d in (("evening", "D"), ("morning", "A")):
        sel = [r for r in dry if r["dir"] == d]
        an = anomalies(sel, key, xkey, centre_only)
        rho, p = rho_perm(an, rng)
        e = {"n": len(an), "rho": rho, "p": p, "slope_db_per_kpa": slope_ci(an, rng), "anom": an}
        for site in SITES:
            a = [x for x in an if x["site"] == site]
            e[f"n_{site}"] = len(a)
            e[f"rho_{site}"] = float(spearmanr([x["x"] for x in a], [x["y"] for x in a])[0]) if len(a) >= 3 else None
        res[when] = e
    return res


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    for site, f in (("wisconsin", rows_wisconsin), ("bakota", rows_bakota)):
        vpd, rain = weather(site)
        rows += f(vpd, rain)
    rows = [r for r in rows if r["vpd"] is not None]
    print(len(rows), "passes in season;", sum(r["rain24"] < RAIN_MM for r in rows), "dry (24 h)")
    print("passes with missing rain hours: 24 h window", sum(r["rain24_missing"] > 0 for r in rows),
          "(of them counted dry:", sum(r["rain24_missing"] > 0 and r["rain24"] < RAIN_MM for r in rows), ");",
          "48 h window", sum(r["rain48_missing"] > 0 for r in rows))
    rng = np.random.default_rng(2026)
    main_res = analyse(rows, rng)
    ev, mo = main_res["evening"], main_res["morning"]
    v1 = ("inconclusive" if ev["n"] < 15 else
          "pass" if ev["rho"] <= -0.3 and ev["p"] <= 0.05 else "fail")
    small = [s for s in SITES if ev[f"n_{s}"] < 6]
    v2 = ("inconclusive" if small else
          "pass" if all(ev[f"rho_{s}"] < 0 for s in SITES) else "fail")
    v3 = ("not judged (TW1 did not pass)" if v1 != "pass" else
          "inconclusive" if mo["n"] < 15 else "pass" if mo["rho"] >= ev["rho"] + 0.2 else "fail")
    res = {"rows": rows, "main": main_res,
           "TW1": {"n": ev["n"], "rho": ev["rho"], "p": ev["p"], "verdict": v1},
           "TW2": {s: {"n": ev[f"n_{s}"], "rho": ev[f"rho_{s}"]} for s in SITES} | {"verdict": v2},
           "TW3": {"n": mo["n"], "rho_m": mo["rho"], "p_m": mo["p"], "rho_e": ev["rho"], "verdict": v3},
           "reported": {
               "slope_evening_db_per_kpa_[est,p5,p95]": ev["slope_db_per_kpa"],
               "centre_only": analyse(rows, rng, centre_only=True),
               "rain48": analyse(rows, rng, rain_key="rain48"),
               "hh": analyse(rows, rng, key="forest_hh"),
               "conifer_hv_wisconsin": analyse([r for r in rows if r["site"] == "wisconsin"], rng, key="conif_hv"),
               "placebo_vpd_after": analyse(rows, rng, xkey="vpd_after")}}
    for k, v in res["reported"].items():
        if isinstance(v, dict):
            for when in ("evening", "morning"):
                v[when].pop("anom", None)
    (OUT / "tw_results.json").write_text(json.dumps(res, indent=1, default=str), encoding="utf-8")
    for k in ("TW1", "TW2", "TW3"):
        print(k, json.dumps(res[k]))
    print("slope evening dB/kPa [est, p5, p95]:", ev["slope_db_per_kpa"])
    for k, v in res["reported"].items():
        if isinstance(v, dict):
            print(f"  {k}: " + "; ".join(
                f"{w} n={v[w]['n']} rho={v[w]['rho'] if v[w]['rho'] is None else round(v[w]['rho'], 2)} "
                f"p={v[w]['p'] if v[w]['p'] is None else round(v[w]['p'], 3)} "
                f"WI={v[w]['rho_wisconsin'] if v[w]['rho_wisconsin'] is None else round(v[w]['rho_wisconsin'], 2)} "
                f"BK={v[w]['rho_bakota'] if v[w]['rho_bakota'] is None else round(v[w]['rho_bakota'], 2)}"
                for w in ("evening", "morning")))
    plot(res)


def plot(res):
    fig, ax = plt.subplots(1, 2, figsize=(13, 5.2), sharey=True)
    for a, when in zip(ax, ("evening", "morning")):
        e = res["main"][when]
        for site, m, c in (("wisconsin", "o", "C0"), ("bakota", "s", "C3")):
            an = [x for x in e["anom"] if x["site"] == site]
            a.scatter([x["x"] for x in an], [x["y"] for x in an], marker=m, color=c, s=28,
                      label=f"{site.capitalize()} (ρ = {e[f'rho_{site}']:.2f}, n = {len(an)})"
                      if e[f"rho_{site}"] is not None else site)
        if e["slope_db_per_kpa"]:
            xs = np.linspace(min(x["x"] for x in e["anom"]), max(x["x"] for x in e["anom"]), 2)
            a.plot(xs, e["slope_db_per_kpa"][0] * xs, "k-", lw=0.8)
        a.axhline(0, color="0.8", lw=0.6)
        a.axvline(0, color="0.8", lw=0.6)
        t = "dusk passes (~19:05 local)" if when == "evening" else "dawn passes (~05:00 local)"
        a.set(title=f"{t}: pooled ρ = {e['rho']:.2f}, p = {e['p']:.3f}, n = {e['n']}",
              xlabel="afternoon VPD anomaly, kPa")
        a.legend(fontsize=8)
    ax[0].set_ylabel("forest HV anomaly, dB (per track, detrended)")
    fig.suptitle("TW: does the forest HV echo drop after a dry afternoon? NISAR GCOV, dry passes, 15 Jun – 30 Sep 2026 "
                 "(pre-registered, pushed before the data)", fontsize=10)
    fig.tight_layout()
    fig.savefig(OUT / "tw_tests.png", dpi=120)
    print("figure:", OUT / "tw_tests.png")


if __name__ == "__main__":
    main()
