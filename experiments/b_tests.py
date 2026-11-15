"""B1-B3: Bakota tests (criteria in notes/experiments.md, section B; run b_geometry.py and b_fetch.py first).

uv run python experiments/b_tests.py      -> data/processed/bakota/{b_results.json, b_tests.png}
"""
import json
import sys
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import rasterio  # noqa: E402
from rasterio.warp import Resampling  # noqa: E402
from scipy.ndimage import binary_erosion  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src"), str(ROOT / "experiments")]
import k_fetch  # noqa: E402
from b_fetch import CLASS_DATES, MAX_GAP_DAYS  # noqa: E402
from b_geometry import BBOX, OUT  # noqa: E402
from hachure import gcov  # noqa: E402

MIN_PX = 50
TAN15 = np.tan(np.radians(15))


def db(v):
    with np.errstate(divide="ignore", invalid="ignore"):
        return 10 * np.log10(v)


def med_db(v):
    v = v[np.isfinite(v) & (v > 0)]
    return float(np.median(db(v))) if v.size >= MIN_PX else None


def main():
    G = dict(np.load(OUT / "geometry.npz"))
    dem_tf = rasterio.Affine(*G["transform"])
    gy, gx = np.gradient(G["z"], dem_tf.a)
    slope = np.degrees(np.arctan(np.hypot(gx, gy)))
    s2 = json.loads((OUT / "s2_dates.json").read_text(encoding="utf-8"))
    used = [datetime.fromisoformat(d) for d in s2["used"]]
    rows = []
    for p in sorted(gcov.cache_dir(BBOX, "A").glob("NISAR_*.npz")):
        w = gcov.Window.load(p)
        L = dict(np.load(OUT / f"gridA_{p.stem}.npz"))
        t = int(w.name.split("_")[5])
        dx, dy = w.x[1] - w.x[0], w.y[1] - w.y[0]
        atf = rasterio.Affine(dx, 0, w.x[0] - dx / 2, 0, dy, w.y[0] - dy / 2)
        sl = k_fetch.warp(slope, dem_tf, "EPSG:32635", atf, w.hh.shape, how=Resampling.nearest, epsg=w.epsg)
        st = L[f"geo_slope_toward_{t}"]
        valid = (w.mask > 0) & (w.mask < 255) & (w.hh > 0) & (w.hv > 0)
        forest = valid & (L["wc_tree"] >= 0.9) & (sl >= 15)
        facing, away = forest & (st > TAN15), forest & (st < -TAN15)
        water = binary_erosion((L[f"s2_{CLASS_DATES[0]}_water"] == 1) & (L[f"s2_{CLASS_DATES[1]}_water"] == 1)) & valid
        crop = valid & (L["wc_crop"] >= 0.9)
        g2s = w.extra.get("g2s")
        r = {"name": w.name, "date": w.date.date().isoformat(), "track": t, "dir": w.name.split("_")[6],
             "n_facing": int(facing.sum()), "n_away": int(away.sum()), "n_water": int(water.sum())}
        for pol, a in (("hh", w.hh), ("hv", w.hv)):
            r[f"facing_{pol}"], r[f"away_{pol}"] = med_db(a[facing]), med_db(a[away])
            if g2s is not None:
                r[f"facing_s_{pol}"], r[f"away_s_{pol}"] = med_db((a * g2s)[facing]), med_db((a * g2s)[away])
            r[f"water_{pol}"], r[f"crop_{pol}"] = med_db(a[water]), med_db(a[crop])
        neb = w.extra.get("neb_hv")
        nb = None if neb is None else neb[water]
        r["neb_ok"] = bool(nb is not None and nb.size >= MIN_PX and np.isfinite(nb).all() and (nb > 0).all())
        r["water_hv_over_neb"] = float(np.median(db(w.hv[water] / nb))) if r["neb_ok"] else None
        r["water_neb_hv"] = float(np.median(db(nb))) if r["neb_ok"] else None
        if r["water_hh"] is not None and r["crop_hh"] is not None:          # B3: radar and optical water area
            thr = 0.5 * (r["water_hh"] + r["crop_hh"])
            px = abs(dx * dy) / 1e6
            r["radar_water_km2"] = float((valid & (db(w.hh) < thr)).sum() * px)
            near = min(used, key=lambda d: abs(d - w.date))
            if abs((near - w.date).total_seconds()) <= MAX_GAP_DAYS * 86400:
                ow = L[f"s2_{near.date().isoformat()}_water"]
                r["s2_date"] = near.date().isoformat()
                r["optical_water_km2"] = float(np.nansum(np.where(valid, ow, 0)) * px)
        rows.append(r)
    print(f"{len(rows)} band A windows")
    res = {"rows": rows}
    tracks = sorted({r["track"] for r in rows})

    b1 = {}
    for t in tracks:
        rr = [r for r in rows if r["track"] == t and "2026-06-01" <= r["date"] <= "2026-09-30"]
        e = {"n_dates": len(rr)}
        for pol in ("hh", "hv"):
            for tag, f, a in (("gamma", "facing_", "away_"), ("sigma", "facing_s_", "away_s_")):
                d = [r[f + pol] - r[a + pol] for r in rr if r.get(f + pol) is not None and r.get(a + pol) is not None]
                e[f"d_{tag}_{pol}"] = float(np.median(d)) if d else None
        b1[t] = e
    ok = [e for e in b1.values() if e["d_gamma_hh"] is not None and e["d_gamma_hv"] is not None]
    res["B1"] = {"per_track": b1, "verdict": "inconclusive" if len(ok) < 3 else
                 ("pass" if sum(abs(e["d_gamma_hh"]) <= 1 and abs(e["d_gamma_hv"]) <= 1 for e in ok) >= 3 else "fail")}

    k = [r for r in rows if r["neb_ok"]]
    res["B2a"] = {"n_neb_ok": len(k), "n_neb_broken": sum(not r["neb_ok"] for r in rows if r["n_water"] >= MIN_PX),
                  "share": float(np.mean([abs(r["water_hv_over_neb"]) <= 2 for r in k])) if k else None,
                  "values": [round(r["water_hv_over_neb"], 2) for r in k]}
    res["B2a"]["verdict"] = "inconclusive" if len(k) < 5 else ("pass" if res["B2a"]["share"] >= 0.8 else "fail")
    k = [r for r in rows if r["water_hh"] is not None and r["water_hv"] is not None]
    s = float(np.mean([r["water_hh"] >= r["water_hv"] - 0.5 for r in k])) if k else None
    res["B2b"] = {"n": len(k), "share": s, "verdict": "inconclusive" if s is None else ("pass" if s >= 0.8 else "fail")}

    m = [r for r in rows if r.get("optical_water_km2")]
    ratio = [r["radar_water_km2"] / r["optical_water_km2"] for r in m]
    res["B3a"] = {"n_pairs": len(m), "share": float(np.mean([0.85 <= q <= 1.15 for q in ratio])) if m else None,
                  "ratios": [round(q, 3) for q in ratio]}
    res["B3a"]["verdict"] = "inconclusive" if len(m) < 6 else ("pass" if res["B3a"]["share"] >= 0.7 else "fail")
    rho = float(spearmanr([r["radar_water_km2"] for r in m], [r["optical_water_km2"] for r in m])[0]) if len(m) >= 6 else None
    res["B3b"] = {"n_pairs": len(m), "spearman": rho,
                  "verdict": "inconclusive" if rho is None else ("pass" if rho >= 0.5 else "fail")}

    (OUT / "b_results.json").write_text(json.dumps(res, indent=1, default=str), encoding="utf-8")
    for k in ("B1", "B2a", "B2b", "B3a", "B3b"):
        print(k, json.dumps(res[k], default=str)[:600])
    plot(res, tracks)


def plot(res, tracks):
    rows = res["rows"]
    fig, ax = plt.subplots(1, 3, figsize=(18, 5.2))
    a = ax[0]
    xs = np.arange(len(tracks))
    for i, (tag, pol, c) in enumerate((("gamma", "hh", "C0"), ("gamma", "hv", "C1"), ("sigma", "hh", "0.6"),
                                       ("sigma", "hv", "0.8"))):
        v = [res["B1"]["per_track"][t][f"d_{tag}_{pol}"] or 0 for t in tracks]
        a.bar(xs + (i - 1.5) * 0.2, v, 0.2, color=c, label=f"{'γ⁰ (RTC)' if tag == 'gamma' else 'σ⁰'} {pol.upper()}")
    a.axhspan(-1, 1, color="0.92", zorder=0)
    a.set_xticks(xs, [f"{t}" for t in tracks])
    a.set(xlabel="track", ylabel="facing − away slopes, dB", title="B1: forested slopes facing vs away from the radar")
    a.legend(fontsize=7)
    a = ax[1]
    k = [r for r in rows if r["neb_ok"]]
    a.scatter([r["water_neb_hv"] for r in k], [r["water_hv"] for r in k], s=18, label="water HV (NEB finite)")
    lo, hi = -40, -5
    a.plot([lo, hi], [lo, hi], "k-", lw=0.8)
    a.fill_between([lo, hi], [lo - 2, hi - 2], [lo + 2, hi + 2], color="0.88", zorder=0)
    a.set(xlabel="product NEB HV, dB", ylabel="canyon water γ⁰ HV, dB", title="B2: canyon water vs the stated noise floor")
    a = ax[2]
    m = sorted((r for r in rows if r.get("optical_water_km2")), key=lambda r: r["date"])
    d = [datetime.fromisoformat(r["date"]) for r in m]
    a.plot(d, [r["radar_water_km2"] for r in m], "o-", label="radar (HH threshold)")
    a.plot(d, [r["optical_water_km2"] for r in m], "s--", label="Sentinel-2 (≤ 5 days)")
    a.set(ylabel="water area in the window, km²", title="B3: water area, radar vs optical")
    a.legend(fontsize=7)
    a.tick_params(axis="x", rotation=30)
    fig.suptitle("Bakota, Dniester canyon, NISAR GCOV band A: tests B1–B3 (pre-registered, pushed before the data)")
    fig.tight_layout()
    fig.savefig(OUT / "b_tests.png", dpi=120)
    print("figure:", OUT / "b_tests.png")


if __name__ == "__main__":
    main()
