"""W1-W8: Wisconsin quad-pol radio-physics tests (criteria in notes/experiments.md, section W; run experiments/w_fetch.py first).

uv run python experiments/w_tests.py      -> data/processed/wisconsin/{w_results.json, w_tests.png}
"""
import json
import math
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
from k_tests import MIN_PX, db, frozen_state, on_grid  # noqa: E402
from w_fetch import BBOX, EPSG, OUT, S2_DATES, SUB  # noqa: E402

JUL, NOV = S2_DATES
CLS = ("water", "decid", "conif", "crop")


def classes(L):
    """Class masks from WorldCover + Sentinel-2 (definitions: notes/experiments.md, section W)."""
    g = lambda k: L[f"s2_{k}"]  # noqa: E731
    leafy = (L["wc_tree"] >= 0.9) & (g(f"{JUL}_ndvi") >= 0.75)
    return {"water": (g(f"{JUL}_water") == 1) & (g(f"{NOV}_water") == 1) & (g(f"{JUL}_ndwi") >= 0.3),
            "decid": leafy & (g(f"{NOV}_ndvi") <= 0.45),
            "conif": leafy & (g(f"{NOV}_ndvi") >= 0.60),
            "crop": (L["wc_crop"] >= 0.9) & (L["wc_water"] == 0)}


def at_least(k, n):
    return n > 0 and k >= math.ceil(2 * n / 3)


def band_b(T):
    G = dict(np.load(OUT / "grid.npz"))
    tf, shape = rasterio.Affine(*G["transform"]), G["wc_water"].shape
    C = classes(G)
    rows = []
    for p in sorted(gcov.cache_dir(BBOX, "B").glob("NISAR_*.npz")):
        w = gcov.Window.load(p)
        L = on_grid(w, tf, shape, epsg=EPSG)
        r = {"name": w.name, "date": w.date.date().isoformat(), "track": int(w.name.split("_")[5]),
             "state": frozen_state(T[w.name]), "t": T[w.name]["t_overpass"]}
        for c, m in C.items():
            for k in ("hh", "hv", "inc", "neb_hv"):
                v = L[k][m & np.isfinite(L["hh"])]
                v = v[np.isfinite(v)]
                r[f"{c}_{k}"] = float(np.median(v if k == "inc" else db(v))) if v.size >= MIN_PX else None
        m = C["water"] & np.isfinite(L["hv"])
        neb = L["neb_hv"][m]
        r["water_neb_broken"] = bool(m.sum() >= MIN_PX and (~np.isfinite(neb) | (neb <= 0)).any())
        r["water_hv_over_neb"] = (float(np.median(db(L["hv"][m] / neb))) if m.sum() >= MIN_PX and not r["water_neb_broken"]
                                  else None)
        rows.append(r)
    return rows, {k: int(v.sum()) for k, v in C.items()}


def coherency(c):
    """3x3 Pauli coherency matrix from the class-mean 4x4 covariance terms (term XY = <X Y*>), cross-pol symmetrised."""
    cj = np.conj
    t11 = c["hhhh"] + c["vvvv"] + 2 * c["hhvv"].real
    t22 = c["hhhh"] + c["vvvv"] - 2 * c["hhvv"].real
    t33 = c["hvhv"] + c["vhvh"] + 2 * c["hvvh"].real
    t12 = c["hhhh"] - c["vvvv"] - c["hhvv"] + cj(c["hhvv"])
    t13 = c["hhhv"] + c["hhvh"] + cj(c["hvvv"]) + cj(c["vhvv"])
    t23 = c["hhhv"] + c["hhvh"] - cj(c["hvvv"]) - cj(c["vhvv"])
    T = 0.5 * np.array([[t11, t12, t13], [cj(t12), t22, t23], [cj(t13), cj(t23), t33]], complex)
    return T


def h_alpha(T):
    lam, vec = np.linalg.eigh(T)
    lam = np.clip(lam.real, 1e-12, None)
    p = lam / lam.sum()
    H = float(-(p * np.log(p) / np.log(3)).sum())
    alpha = float(np.degrees((p * np.arccos(np.clip(np.abs(vec[0, :]), 0, 1))).sum()))
    return H, alpha


def faraday_deg(c):
    """Bickel-Bates estimate from the covariance: a = HV - VH, b = HH + VV.

    Sign convention: a reciprocal target S seen through a one-way Faraday rotation Omega is M = R(Omega) S R(Omega) with
    R(Omega) = [[cos, sin], [-sin, cos]] in the (H, V) basis; this function returns +Omega for such an M
    (experiments/w_selfcheck.py)."""
    aa = c["hvhv"] + c["vhvh"] - 2 * c["hvvh"].real
    bb = c["hhhh"] + c["vvvv"] + 2 * c["hhvv"].real
    ab = np.conj(c["hhhv"]) + c["hvvv"] - np.conj(c["hhvh"]) - c["vhvv"]
    return float(-np.degrees(0.25 * np.angle(bb - aa - 2j * ab.real)))


def band_a(T):
    rows = []
    for p in sorted(gcov.cache_dir(SUB, "A", offdiag=True).glob("NISAR_*.npz")):
        w = gcov.Window.load(p)
        f = OUT / f"gridA_{p.stem}.npz"
        if not f.exists():                         # e.g. the test read outside the monthly pick
            continue
        C = classes(dict(np.load(f)))
        terms = {"hhhh": w.hh, "hvhv": w.hv, "vhvh": w.extra["vh"], "vvvv": w.extra["vv"],
                 **{k: w.extra[k] for k in ("hhhv", "hhvh", "hhvv", "hvvh", "hvvv", "vhvv")}}
        ok = (w.mask > 0) & (w.mask < 255) & (w.hh > 0) & (w.hv > 0) & (terms["vhvh"] > 0) & (terms["vvvv"] > 0)
        r = {"name": w.name, "date": w.date.date().isoformat(), "track": int(w.name.split("_")[5]),
             "state": frozen_state(T[w.name]) if w.name in T else None}
        C["forest"] = C["decid"] | C["conif"]
        for c, m in C.items():
            m = m & ok
            r[f"{c}_n"] = int(m.sum())
            if m.sum() < MIN_PX:
                continue
            mean = {k: v[m].astype(complex).mean() for k, v in terms.items()}
            for k, key in (("hhhh", "hh"), ("hvhv", "hv"), ("vhvh", "vh"), ("vvvv", "vv")):
                r[f"{c}_{key}"] = float(np.median(db(terms[k][m])))
                nb = w.extra.get(f"neb_{key}")
                r[f"{c}_neb_{key}"] = float(np.median(db(nb[m]))) if nb is not None else None
                ok_nb = None if nb is None else (np.isfinite(nb[m]) & (nb[m] > 0))
                r[f"{c}_{key}_over_neb"] = (float(np.median(db(terms[k][m] / nb[m]))) if ok_nb is not None
                                            and ok_nb.all() else None)       # median of 10*log10(term / NEB), as in K1c
            r[f"{c}_hv_vh_db"] = float(np.median(db(w.hv[m] / terms["vhvh"][m])))
            r[f"{c}_xcoh"] = float(abs(mean["hvvh"]) / math.sqrt(mean["hvhv"].real * mean["vhvh"].real))
            r[f"{c}_H"], r[f"{c}_alpha"] = h_alpha(coherency(mean))
            r[f"{c}_faraday_deg"] = faraday_deg(mean)
        rows.append(r)
    return rows


def share(xs):
    return float(np.mean(xs)) if xs else None


def main():
    T = json.loads((OUT / "temps.json").read_text(encoding="utf-8"))
    rows, cls_px = band_b(T)
    arows = band_a(T)
    print("band B class px:", cls_px, "| windows B", len(rows), "A", len(arows))
    res = {"class_px_B": cls_px, "rows_B": rows, "rows_A": arows}
    tracks = sorted({r["track"] for r in rows})

    k = [r for r in rows if r["water_hv"] is not None]          # scenes with a water class
    fin = [r["water_hv_over_neb"] for r in k if r["water_hv_over_neb"] is not None]
    res["W1a"] = {"n": len(k), "share": share([v is not None and abs(v) <= 2 for v in (r["water_hv_over_neb"] for r in k)]),
                  "n_neb_broken": sum(r["water_neb_broken"] for r in k), "n_finite": len(fin),
                  "n_finite_within_2dB": sum(abs(v) <= 2 for v in fin),
                  "finite_range_db": [round(min(fin), 2), round(max(fin), 2)] if fin else None}
    k = [r for r in rows if r["water_hh"] is not None and r["water_hv"] is not None]
    res["W1b"] = {"n": len(k), "share": share([r["water_hh"] >= r["water_hv"] - 0.5 for r in k])}
    k = [r for r in arows if r.get("water_hv_over_neb") is not None and r.get("water_vh_over_neb") is not None]
    res["W1c"] = {"n": len(k), "share": share([abs(r["water_hv_over_neb"]) <= 2 and abs(r["water_vh_over_neb"]) <= 2
                                               for r in k]),
                  "hv_over_neb": [round(r["water_hv_over_neb"], 2) for r in k]}
    for t in ("W1a", "W1b", "W1c"):
        s = res[t]["share"]
        res[t]["verdict"] = "inconclusive" if s is None else ("pass" if s >= 0.8 else "fail")

    w2 = {}
    for t in tracks:
        rr = [r for r in rows if r["track"] == t and r["state"] == "unfrozen"]
        months = [r["date"] for r in rr]
        elig = (len(rr) >= 4 and any(d[5:7] in ("06", "07", "08") for d in months)
                and any("10-15" <= d[5:] <= "11-30" for d in months))
        e = {"n": len(rr), "eligible": elig}
        for c in ("decid", "conif", "crop"):
            v = [r[f"{c}_hv"] for r in rr if r[f"{c}_hv"] is not None]
            e[f"{c}_hv_amp"] = float(max(v) - min(v)) if len(v) >= 2 else None
        w2[t] = e
    el = [e for e in w2.values() if e["eligible"] and None not in (e["decid_hv_amp"], e["conif_hv_amp"], e["crop_hv_amp"])]
    n = len(el)
    res["W2"] = {"per_track": w2, "eligible": n}
    for name, ok in (("W2a", [e["decid_hv_amp"] <= 2.0 for e in el]),
                     ("W2b", [e["decid_hv_amp"] < e["crop_hv_amp"] for e in el]),
                     ("W2c", [e["decid_hv_amp"] - e["conif_hv_amp"] <= 1.0 for e in el])):
        res["W2"][name] = "inconclusive" if n < 2 else ("pass" if at_least(sum(ok), n) else "fail")

    w3 = {}
    for t in tracks:
        rr = [r for r in rows if r["track"] == t and r["state"] == "unfrozen" and "2026-06-01" <= r["date"] <= "2026-09-30"]
        if not rr:
            continue
        e = {"n": len(rr)}
        for c in ("decid", "conif"):
            for q in ("hh", "hv", "inc"):
                v = [r[f"{c}_{q}"] for r in rr if r[f"{c}_{q}"] is not None]
                e[f"{c}_{q}"] = float(np.median(v)) if v else None
        w3[t] = e
    th = [e[f"{c}_inc"] for e in w3.values() for c in ("decid", "conif") if e[f"{c}_inc"] is not None]
    spread = float(max(th) - min(th)) if th else 0.0
    res["W3"] = {"per_track": w3, "theta_spread_deg": spread}
    if spread < 5:
        res["W3"]["verdict"] = "inconclusive"
    else:
        rg = {f"{c}_{q}": float(np.ptp([e[f"{c}_{q}"] for e in w3.values() if e[f"{c}_{q}"] is not None]))
              for c in ("decid", "conif") for q in ("hh", "hv")}
        res["W3"].update(ranges_db=rg, verdict="pass" if all(v <= 1.5 for v in rg.values()) else "fail")

    w4 = []
    for t in tracks:
        rr = [r for r in rows if r["track"] == t and "2025-10-01" <= r["date"] <= "2026-01-31"]
        for c in ("decid", "conif", "crop"):
            f = [r[f"{c}_hh"] for r in rr if r["state"] == "frozen" and r[f"{c}_hh"] is not None]
            u = [r[f"{c}_hh"] for r in rr if r["state"] == "unfrozen" and r[f"{c}_hh"] is not None]
            if len(f) >= 1 and len(u) >= 2:
                w4.append({"track": t, "class": c, "n_frozen": len(f), "n_unfrozen": len(u),
                           "diff_db": float(np.median(f) - np.median(u))})
    res["W4"] = {"cases": w4, "verdict": "inconclusive" if len(w4) < 2 else
                 ("pass" if sum(e["diff_db"] <= -1.5 for e in w4) >= len(w4) / 2 else "fail")}

    for c in ("forest", "crop"):
        k = [r for r in arows if r.get(f"{c}_hv_vh_db") is not None]
        res[f"W5a_{c}"] = {"n": len(k), "share": share([abs(r[f"{c}_hv_vh_db"]) <= 0.5 for r in k]),
                           "values": [round(r[f"{c}_hv_vh_db"], 2) for r in k]}
        res[f"W5b_{c}"] = {"n": len(k), "share": share([r[f"{c}_xcoh"] >= 0.9 for r in k]),
                           "values": [round(r[f"{c}_xcoh"], 3) for r in k]}
    k = [r for r in arows if r.get("forest_vv") is not None]
    res["W6a"] = {"n": len(k), "share": share([abs(r["forest_vv"] - r["forest_hh"]) <= 1.0 for r in k]),
                  "values": [round(r["forest_vv"] - r["forest_hh"], 2) for r in k]}
    k = [r for r in arows if r.get("water_vv") is not None and r["water_neb_hh"] is not None and r["water_neb_vv"] is not None
         and r["water_hh"] >= r["water_neb_hh"] + 3 and r["water_vv"] >= r["water_neb_vv"] + 3]
    res["W6b"] = {"n_gated": len(k), "share": share([r["water_vv"] - r["water_hh"] >= 1.0 for r in k]),
                  "values": [round(r["water_vv"] - r["water_hh"], 2) for r in k]}
    k = [r for r in arows if r["state"] == "unfrozen" and r.get("forest_H") is not None]
    res["W7"] = {"n": len(k), "share": share([r["forest_H"] >= 0.7 and 40 <= r["forest_alpha"] <= 55 for r in k]),
                 "forest": [(round(r["forest_H"], 3), round(r["forest_alpha"], 1)) for r in k]}
    for name in ("W5a_forest", "W5a_crop", "W5b_forest", "W5b_crop", "W6a", "W7"):
        s = res[name]["share"]
        res[name]["verdict"] = "inconclusive" if s is None else ("pass" if s >= 0.8 else "fail")
    s = res["W6b"]["share"]
    res["W6b"]["verdict"] = "inconclusive" if res["W6b"]["n_gated"] < 5 else ("pass" if s >= 0.8 else "fail")
    res["W8"] = {"forest_faraday_deg": [(r["date"], r["track"], round(r["forest_faraday_deg"], 2))
                                        for r in arows if r.get("forest_faraday_deg") is not None]}

    (OUT / "w_results.json").write_text(json.dumps(res, indent=1, default=str), encoding="utf-8")
    for k in sorted(x for x in res if x.startswith("W")):
        v = res[k]
        print(k, json.dumps({kk: vv for kk, vv in v.items() if kk not in ("per_track", "cases")}, default=str)[:700])
    plot(res, tracks)


def plot(res, tracks):
    rows, arows = res["rows_B"], res["rows_A"]
    col = {t: c for t, c in zip(tracks, plt.cm.tab10.colors)}
    fig, ax = plt.subplots(2, 4, figsize=(21, 9.5))
    a = ax[0, 0]
    for t in tracks:
        rr = [r for r in rows if r["track"] == t and r["water_hv"] is not None]
        a.scatter([r["water_neb_hv"] for r in rr], [r["water_hv"] for r in rr], s=16, color=col[t], label=f"B {t}")
    for r in arows:
        if r.get("water_hv") is not None:
            a.scatter(r["water_neb_hv"], r["water_hv"], marker="*", s=90, color="k")
    lo, hi = -35, -5
    a.plot([lo, hi], [lo, hi], "k-", lw=0.8)
    a.fill_between([lo, hi], [lo - 2, hi - 2], [lo + 2, hi + 2], color="0.85", zorder=0)
    a.set(xlabel="product NEB HV, dB", ylabel="calm water γ⁰ HV, dB", title="W1: water vs the stated noise floor (★ band A)")
    a.legend(fontsize=7)
    for j, c in enumerate(("decid", "conif", "crop")):
        if j == 2:
            a = ax[0, 3]
        else:
            a = ax[0, 1 + j]
        for t in tracks:
            rr = [r for r in rows if r["track"] == t and r[f"{c}_hv"] is not None]
            d = [datetime.fromisoformat(r["date"]) for r in rr]
            a.plot(d, [r[f"{c}_hv"] for r in rr], "-o", ms=3, color=col[t], label=f"HV {t}")
            a.plot(d, [r[f"{c}_hh"] for r in rr], "--", lw=0.8, color=col[t])
            fz = [(dd, r[f"{c}_hh"]) for dd, r in zip(d, rr) if r["state"] == "frozen"]
            if fz:
                a.scatter(*zip(*fz), marker="x", color="b", s=30)
        a.set(title=f"W2/W4 {c}: HV solid, HH dashed, × frozen", ylabel="γ⁰, dB")
        a.tick_params(axis="x", rotation=30)
    ax[0, 3].legend(fontsize=7)
    a = ax[1, 0]
    for c, m in (("decid", "o"), ("conif", "s")):
        for q, fill in (("hh", True), ("hv", False)):
            pts = sorted((e[f"{c}_inc"], e[f"{c}_{q}"]) for e in res["W3"]["per_track"].values()
                         if e[f"{c}_inc"] is not None and e[f"{c}_{q}"] is not None)
            if pts:
                a.plot(*zip(*pts), marker=m, mfc=None if fill else "none", label=f"{c} {q.upper()}")
    a.set(xlabel="θ, deg", ylabel="γ⁰, dB (Jun–Sep median)", title="W3: angular signature")
    a.legend(fontsize=7)
    a = ax[1, 1]
    dd = lambda rr: [datetime.fromisoformat(r["date"]) for r in rr]  # noqa: E731
    for c, mk in (("forest", "o"), ("crop", "s")):
        rr = [r for r in arows if r.get(f"{c}_hv_vh_db") is not None]
        a.plot(dd(rr), [r[f"{c}_hv_vh_db"] for r in rr], mk, label=f"{c} HV−VH, dB")
        a.plot(dd(rr), [r[f"{c}_vv"] - r[f"{c}_hh"] for r in rr], mk, mfc="none", label=f"{c} VV−HH, dB")
    rr = [r for r in arows if r.get("water_vv") is not None]
    a.plot(dd(rr), [r["water_vv"] - r["water_hh"] for r in rr], "^", label="water VV−HH, dB")
    a.axhspan(-0.5, 0.5, color="0.9", zorder=0)
    a.set(title="W5/W6: reciprocity and co-pol ratio (band A)")
    a.legend(fontsize=7)
    a.tick_params(axis="x", rotation=30)
    a = ax[1, 2]
    hh = np.linspace(0, 1, 50)
    for c, mk in (("forest", "o"), ("crop", "s"), ("water", "^")):
        rr = [r for r in arows if r.get(f"{c}_H") is not None]
        a.scatter([r[f"{c}_H"] for r in rr], [r[f"{c}_alpha"] for r in rr], marker=mk, label=c)
    a.add_patch(plt.Rectangle((0.7, 40), 0.3, 15, fill=False, ls="--", color="g"))
    a.set(xlim=(0, 1), ylim=(0, 90), xlabel="entropy H", ylabel="ᾱ, deg",
          title="W7: Cloude–Pottier H/α (dashed = forest prediction)")
    a.legend(fontsize=7)
    a = ax[1, 3]
    rr = [r for r in arows if r.get("forest_faraday_deg") is not None]
    for t in tracks:
        q = [r for r in rr if r["track"] == t]
        a.plot(dd(q), [r["forest_faraday_deg"] for r in q], "o", color=col[t], label=f"track {t}")
    a.axhspan(-3, 3, color="0.9", zorder=0)
    a.set(ylabel="Ω, deg", title="W8: Faraday rotation estimate (forest, band A)")
    a.legend(fontsize=7)
    a.tick_params(axis="x", rotation=30)
    fig.suptitle("Northern Wisconsin, NISAR GCOV quad-pol L-band: tests W1–W8 (pre-registered)")
    fig.tight_layout()
    fig.savefig(OUT / "w_tests.png", dpi=120)
    print("figure:", OUT / "w_tests.png")


if __name__ == "__main__":
    main()
