"""K0-K4: Kakhovka bed radio-physics tests (criteria in notes/experiments.md, section K; run experiments/k_fetch.py first).

uv run python experiments/k_tests.py      -> data/processed/kakhovka/{k_results.json, k_tests.png}
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

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src"), str(ROOT / "experiments")]
from hachure import gcov  # noqa: E402
from k_fetch import BBOX, EPSG, OUT, S2_DATES, warp  # noqa: E402

AUG, FEB = S2_DATES
MIN_PX = 50                                   # fewer class pixels in a window: the window is skipped for that class


def db(v):
    with np.errstate(divide="ignore", invalid="ignore"):
        return 10 * np.log10(v)


def classes(L, prefix=""):
    """Class masks from the radar-independent layers (definitions: notes/experiments.md, section K)."""
    g = lambda k: L[f"{prefix}s2_{k}"]  # noqa: E731
    dry = (g(f"{AUG}_water") == 0) & (g(f"{FEB}_water") == 0)
    water = (g(f"{AUG}_water") == 1) & (g(f"{FEB}_water") == 1) & (g(f"{AUG}_ndwi") >= 0.3)
    if "wc_water" not in L:
        return {"water": water}
    bed = L["wc_water"] >= 0.9
    return {"water": water,
            "veg": bed & dry & (g(f"{AUG}_ndvi") >= 0.6),
            "bare": bed & dry & (g(f"{AUG}_ndvi") < 0.25),
            "crop": (L["wc_crop"] >= 0.9) & (L["wc_water"] == 0)}


def on_grid(w, tf, shape):
    """Window layers (linear) averaged onto the 80 m analysis grid; invalid samples dropped."""
    dx, dy = w.x[1] - w.x[0], w.y[1] - w.y[0]
    src_tf = rasterio.Affine(dx, 0, w.x[0] - dx / 2, 0, dy, w.y[0] - dy / 2)
    bad = (w.mask == 0) | (w.mask == 255)
    out = {}
    for k, v in [("hh", w.hh), ("hv", w.hv)] + list(w.extra.items()):
        v = v.astype("float64").copy()
        v[bad | ~np.isfinite(v)] = np.nan
        if k in ("hh", "hv", "vv", "vh"):
            v[v <= 0] = np.nan
        out[k] = warp(v, src_tf, f"EPSG:{w.epsg}", tf, shape)
    return out


def frozen_state(t):
    if t["t_overpass"] is None or t["t_24h_mean"] is None:
        return None
    if t["t_overpass"] < -2 and t["t_24h_mean"] < 0:
        return "frozen"
    return "unfrozen" if t["t_overpass"] > 2 else None


def main():
    G = dict(np.load(OUT / "grid.npz"))
    tf, shape = rasterio.Affine(*G["transform"]), G["wc_water"].shape
    C = classes(G)
    print("class pixels:", {k: int(v.sum()) for k, v in C.items()})
    T = json.loads((OUT / "temps.json").read_text(encoding="utf-8"))
    cache = gcov.cache_dir(BBOX, "B")
    rows, audit = [], []
    for p in sorted(cache.glob("NISAR_*.npz")):
        w = gcov.Window.load(p)
        L = on_grid(w, tf, shape)
        track = int(w.name.split("_")[5])
        r = {"name": w.name, "date": w.date.date().isoformat(), "track": track,
             "dir": w.name.split("_")[6], "state": frozen_state(T[w.name]), "t": T[w.name]["t_overpass"]}
        for c, m in C.items():
            for k in ("hh", "hv", "inc", "neb_hh", "neb_hv"):
                v = L[k][m & np.isfinite(L["hh"])]
                v = v[np.isfinite(v)]
                r[f"{c}_{k}"] = (float(np.median(db(v) if k in ("hh", "hv", "neb_hh", "neb_hv") else v))
                                 if v.size >= MIN_PX else None)
            r[f"{c}_n"] = int((m & np.isfinite(L["hh"])).sum())
        m = C["water"] & np.isfinite(L["hv"]) & np.isfinite(L["neb_hv"])
        r["water_hv_over_neb"] = float(np.median(db(L["hv"][m] / L["neb_hv"][m]))) if m.sum() >= MIN_PX else None
        m = C["water"] & np.isfinite(L["hh"]) & np.isfinite(L["neb_hh"])
        r["water_hh_over_neb"] = float(np.median(db(L["hh"][m] / L["neb_hh"][m]))) if m.sum() >= MIN_PX else None
        rows.append(r)
        audit.append({"name": w.name, **w.meta})
    print(f"{len(rows)} band B windows")
    res = {"bbox": BBOX, "class_px": {k: int(v.sum()) for k, v in C.items()}, "rows": rows, "audit": audit}

    # K0: product audit
    def vals(key):
        return [a[key] for a in audit if key in a]
    res["K0"] = {k: sorted({str(v) for v in vals(k)}) for k in
                 ("noiseCorrectionApplied", "faradayRotationApplied", "rfiMitigationApplied", "isFullCovariance",
                  "crosstalk_all_zero", "antpat_A_all_one", "antpat_B_all_one", "faraday_A", "faraday_B",
                  "productVersion", "compositeReleaseId")}
    for k in ("rfi_A_HH", "rfi_A_HV", "rfi_B_HH", "rfi_B_HV"):
        v = np.array(vals(k), float)
        res["K0"][k] = {"min": float(np.nanmin(v)), "median": float(np.nanmedian(v)), "max": float(np.nanmax(v)),
                        "outside_0_1": int(np.sum((v < 0) | (v > 1))), "n": int(v.size)}

    # K1a / K1b
    k1 = [r for r in rows if r["water_hv_over_neb"] is not None]
    ok_a = [abs(r["water_hv_over_neb"]) <= 2 for r in k1]
    ok_b = [r["water_hh"] >= r["water_hv"] - 0.5 for r in k1 if r["water_hh"] is not None and r["water_hv"] is not None]
    res["K1a"] = {"n": len(k1), "share_within_2dB": float(np.mean(ok_a)) if k1 else None,
                  "median_ratio_db": float(np.median([r["water_hv_over_neb"] for r in k1])) if k1 else None,
                  "verdict": ("pass" if np.mean(ok_a) >= 0.8 else "fail") if k1 else "inconclusive"}
    res["K1b"] = {"n": len(ok_b), "share": float(np.mean(ok_b)) if ok_b else None,
                  "verdict": ("pass" if np.mean(ok_b) >= 0.8 else "fail") if ok_b else "inconclusive"}

    # K1c: band A, one scene per track, S2 water on the band A grid
    k1c = []
    for p in sorted(gcov.cache_dir(BBOX, "A").glob("NISAR_*.npz")):
        w = gcov.Window.load(p)
        A = dict(np.load(OUT / f"gridA_{p.stem}.npz"))
        m = classes(A)["water"] & (w.mask > 0) & (w.mask < 255) & (w.hv > 0) & np.isfinite(w.extra["neb_hv"])
        e = {"name": w.name, "date": w.date.date().isoformat(), "n": int(m.sum())}
        if m.sum() >= MIN_PX:
            e["hv_over_neb_db"] = float(np.median(db(w.hv[m] / w.extra["neb_hv"][m])))
            e["hh_over_neb_db"] = float(np.median(db(w.hh[m] / w.extra["neb_hh"][m])))
            e["water_hv_db"], e["neb_hv_db"] = float(np.median(db(w.hv[m]))), float(np.median(db(w.extra["neb_hv"][m])))
            e["within_2dB"] = abs(e["hv_over_neb_db"]) <= 2
        k1c.append(e)
    res["K1c"] = k1c

    tracks = sorted({r["track"] for r in rows})

    # K2: seasonal amplitude on unfrozen dates
    k2 = {}
    for t in tracks:
        rr = [r for r in rows if r["track"] == t and r["state"] == "unfrozen"]
        e = {"n_dates": len(rr)}
        for c in ("veg", "crop"):
            for k in ("hv", "hh"):
                v = [r[f"{c}_{k}"] for r in rr if r[f"{c}_{k}"] is not None]
                e[f"{c}_{k}_amp"] = float(max(v) - min(v)) if len(v) >= 2 else None
        k2[t] = e
    elig = [e for e in k2.values() if e["veg_hv_amp"] is not None]
    res["K2"] = {"per_track": k2, "eligible_tracks": len(elig),
                 "K2a": "inconclusive" if len(elig) < 3 else
                 ("pass" if sum(e["veg_hv_amp"] <= 2.0 for e in elig) >= 3 else "fail"),
                 "K2b": "inconclusive" if len(elig) < 3 else
                 ("pass" if sum(e["crop_hv_amp"] is not None and e["veg_hv_amp"] < e["crop_hv_amp"]
                                for e in elig) >= 3 else "fail")}

    # K3: angular signature, Jun-Sep 2026, unfrozen
    k3 = {}
    for t in tracks:
        rr = [r for r in rows if r["track"] == t and r["state"] == "unfrozen" and "2026-06-01" <= r["date"] <= "2026-09-30"]
        e = {"n_dates": len(rr)}
        for c in ("veg", "bare", "water"):
            for k in ("hh", "hv", "inc"):
                v = [r[f"{c}_{k}"] for r in rr if r[f"{c}_{k}"] is not None]
                e[f"{c}_{k}"] = float(np.median(v)) if v else None
        k3[t] = e
    th = [e[f"{c}_inc"] for e in k3.values() for c in ("veg", "bare") if e[f"{c}_inc"] is not None]
    spread = float(max(th) - min(th)) if th else 0.0
    out3 = {"per_track": k3, "theta_spread_deg": spread}
    if spread < 5:
        out3["verdict"] = "inconclusive"
    else:
        def rng(c, k):
            v = [e[f"{c}_{k}"] for e in k3.values() if e[f"{c}_{k}"] is not None]
            return float(max(v) - min(v)) if len(v) >= 2 else None

        def slope(c):
            pts = [(e[f"{c}_inc"], e[f"{c}_hh"]) for e in k3.values() if e[f"{c}_inc"] is not None and e[f"{c}_hh"] is not None]
            if len(pts) < 2:
                return None
            x, y = np.array(pts).T
            return float(np.polyfit(x, y, 1)[0])
        out3.update(veg_hv_range=rng("veg", "hv"), veg_hh_range=rng("veg", "hh"), bare_hh_slope=slope("bare"),
                    veg_hh_slope=slope("veg"), water_hh_slope=slope("water"))
        out3["K3a"] = ("pass" if out3["veg_hv_range"] is not None and out3["veg_hv_range"] <= 1.5
                       and out3["veg_hh_range"] <= 1.5 else "fail")
        out3["K3b"] = ("pass" if out3["bare_hh_slope"] is not None and out3["veg_hh_slope"] is not None
                       and out3["bare_hh_slope"] <= -0.10 and out3["bare_hh_slope"] < out3["veg_hh_slope"] else "fail")
    res["K3"] = out3

    # K4: frozen vs unfrozen, Nov-Mar
    k4 = []
    for t in tracks:
        rr = [r for r in rows if r["track"] == t and (r["date"][5:7] in ("11", "12", "01", "02", "03"))]
        fr = [r for r in rr if r["state"] == "frozen"]
        un = [r for r in rr if r["state"] == "unfrozen"]
        for c in ("bare", "crop"):
            f = [r[f"{c}_hh"] for r in fr if r[f"{c}_hh"] is not None]
            u = [r[f"{c}_hh"] for r in un if r[f"{c}_hh"] is not None]
            if len(f) >= 1 and len(u) >= 2:
                k4.append({"track": t, "class": c, "n_frozen": len(f), "n_unfrozen": len(u),
                           "diff_db": float(np.median(f) - np.median(u))})
    res["K4"] = {"cases": k4, "verdict": "inconclusive" if len(k4) < 2 else
                 ("pass" if sum(e["diff_db"] <= -1.5 for e in k4) >= len(k4) / 2 else "fail")}

    (OUT / "k_results.json").write_text(json.dumps(res, indent=1, default=str), encoding="utf-8")
    for k in ("K0", "K1a", "K1b", "K1c", "K2", "K3", "K4"):
        v = res[k]
        print(k, json.dumps({kk: vv for kk, vv in v.items() if kk not in ("per_track", "cases")} if isinstance(v, dict)
                            else v, default=str)[:900])
    plot(res, tracks)


def plot(res, tracks):
    rows = res["rows"]
    col = {t: c for t, c in zip(tracks, plt.cm.tab10.colors)}
    fig, ax = plt.subplots(2, 3, figsize=(16, 9.5))
    a = ax[0, 0]                                         # K1: water vs stated noise floor
    for t in tracks:
        rr = [r for r in rows if r["track"] == t and r["water_hv"] is not None]
        a.scatter([r["water_neb_hv"] for r in rr], [r["water_hv"] for r in rr], s=18, color=col[t], label=f"B, track {t}")
    for e in res["K1c"]:
        if "water_hv_db" in e:
            a.scatter(e["neb_hv_db"], e["water_hv_db"], marker="*", s=140, color="k")
    lim = a.get_xlim()[0] - 3, a.get_xlim()[1] + 3
    lo, hi = min(lim[0], a.get_ylim()[0]), max(lim[1], a.get_ylim()[1])
    a.plot([lo, hi], [lo, hi], "k-", lw=0.8)
    a.fill_between([lo, hi], [lo - 2, hi - 2], [lo + 2, hi + 2], color="0.85", zorder=0)
    a.set(xlabel="product noise-equivalent backscatter HV, dB", ylabel="calm water γ⁰ HV, dB",
          title="K1: water vs the stated noise floor (★ = band A)")
    a.legend(fontsize=7)
    for j, (c, k, ttl) in enumerate([("veg", "hv", "K2: vegetated bed"), ("crop", "hv", "K2: cropland control")]):
        a = ax[0, 1 + j]
        for t in tracks:
            rr = [r for r in rows if r["track"] == t and r[f"{c}_{k}"] is not None]
            d = [datetime.fromisoformat(r["date"]) for r in rr]
            a.plot(d, [r[f"{c}_hv"] for r in rr], "-o", ms=3, color=col[t], label=f"HV {t}")
            a.plot(d, [r[f"{c}_hh"] for r in rr], "--", lw=0.8, color=col[t])
            fz = [(dd, r[f"{c}_hh"]) for dd, r in zip(d, rr) if r["state"] == "frozen"]
            if fz:
                a.scatter(*zip(*fz), marker="x", color="b", s=30)
        a.set(title=f"{ttl}: HV (solid), HH (dashed), × frozen", ylabel="γ⁰, dB")
        a.tick_params(axis="x", rotation=30)
    a = ax[1, 0]                                         # K3: gamma0 vs incidence
    for c, m in (("veg", "o"), ("bare", "s"), ("water", "^")):
        for k, fill in (("hh", True), ("hv", False)):
            pts = [(e[f"{c}_inc"], e[f"{c}_{k}"]) for e in res["K3"]["per_track"].values()
                   if e[f"{c}_inc"] is not None and e[f"{c}_{k}"] is not None]
            if pts:
                x, y = zip(*sorted(pts))
                a.plot(x, y, marker=m, mfc=None if fill else "none", label=f"{c} {k.upper()}")
    a.set(xlabel="incidence angle θ, deg", ylabel="γ⁰, dB (Jun–Sep 2026 median)", title="K3: angular signature")
    a.legend(fontsize=7)
    a = ax[1, 1]                                         # K4
    cs = res["K4"]["cases"]
    a.bar(range(len(cs)), [e["diff_db"] for e in cs], color=["0.4" if e["class"] == "bare" else "tan" for e in cs])
    a.axhline(-1.5, color="b", ls="--", lw=0.8)
    a.set_xticks(range(len(cs)), [f"{e['track']}\n{e['class']}" for e in cs], fontsize=7)
    a.set(ylabel="HH frozen − unfrozen, dB", title="K4: frozen ground (dashed = pre-registered −1.5 dB)")
    a = ax[1, 2]                                         # K0
    a.axis("off")
    txt = "K0 product audit\n" + "\n".join(f"{k}: {v}" for k, v in res["K0"].items())
    a.text(0, 1, txt, va="top", family="monospace", fontsize=6.5, wrap=True)
    fig.suptitle("Kakhovka bed, NISAR GCOV L-band: radio-physics tests K0–K4 (pre-registered)")
    fig.tight_layout()
    fig.savefig(OUT / "k_tests.png", dpi=130)
    print("figure:", OUT / "k_tests.png")


if __name__ == "__main__":
    main()
