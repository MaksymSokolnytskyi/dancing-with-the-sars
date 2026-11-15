"""WS1-WS10: Wisconsin story tests (criteria in notes/experiments.md, section WS; run ws_refs.py and ws_fetch.py first).

uv run python experiments/ws_tests.py     -> data/processed/wisconsin/ws/{ws_results.json, ws_*.png, story_*.json}
"""
import gzip  # noqa: F401
import json
import sys
from datetime import datetime, timedelta
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy.ndimage import binary_erosion, uniform_filter  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src"), str(ROOT / "experiments")]
import k_fetch  # noqa: E402  (ASCII CA bundle for GDAL; must come before rasterio)
import rasterio  # noqa: E402
from rasterio.warp import Resampling, transform as warp_transform  # noqa: E402

from hachure import gcov, nisar  # noqa: E402
from k_tests import MIN_PX, db, frozen_state, on_grid  # noqa: E402
from w_fetch import BBOX, EPSG, HEIGHT, OUT, SUB, TRACKS  # noqa: E402
from w_tests import classes, coherency, faraday_deg, h_alpha  # noqa: E402
from ws_fetch import RAW, SUB2  # noqa: E402

RES = OUT / "ws"
N_PERM = 10_000
RNG = np.random.default_rng(2026)
TERMS = ("hhhv", "hhvh", "hhvv", "hvvh", "hvvv", "vhvv")


# ---------------------------------------------------------------- common
def weather():
    h = json.loads((RAW / "weather.json").read_text(encoding="utf-8"))
    t = {datetime.fromisoformat(k): v for k, v in zip(h["time"], h["temperature_2m"])}
    p = {datetime.fromisoformat(k): v for k, v in zip(h["time"], h["precipitation"])}
    return t, p


def state_at(t_utc, temp):
    h = t_utc.replace(minute=0, second=0, microsecond=0)
    last = [temp.get(h - timedelta(hours=i)) for i in range(24)]
    last = [v for v in last if v is not None]
    return frozen_state({"t_overpass": temp.get(h), "t_24h_mean": float(np.mean(last)) if last else None})


def rain_before(t_utc, prec, hours):
    h = t_utc.replace(minute=0, second=0, microsecond=0)
    return float(np.nansum([prec.get(h - timedelta(hours=i)) or 0 for i in range(hours)]))


def utc(name):
    return datetime.strptime(name.split("_")[11], "%Y%m%dT%H%M%S")


def masks80():
    G = dict(np.load(OUT / "grid.npz"))
    R = dict(np.load(OUT / "ws_refs.npz"))
    C = classes(G)
    wet = np.nan_to_num(R["nwi_pfo"]) + np.nan_to_num(R["nwi_pss"]) + np.nan_to_num(R["nwi_pem"])
    M = {"corn": R["cdl_corn"] >= 0.9, "soy": R["cdl_soy"] >= 0.9, "swamp": R["nwi_pfo"] >= 0.9,
         "upland": (C["decid"] | C["conif"]) & (wet == 0), "forest": C["decid"] | C["conif"],
         "lake": binary_erosion(np.nan_to_num(G["wc_water"]) >= 0.9), "crop": C["crop"]}
    return M, rasterio.Affine(*G["transform"]), G["wc_water"].shape


def to_grid(mask, tf, epsg_to, dst_tf, shape):
    v = k_fetch.warp(mask.astype("float64"), tf, f"EPSG:{EPSG}", dst_tf, shape, how=Resampling.nearest, epsg=epsg_to)
    return np.nan_to_num(v) > 0.5


def perm_rho(x, y, groups, rng=RNG, two_sided=False):
    x, y = np.asarray(x, float), np.asarray(y, float)
    rho = float(spearmanr(x, y)[0])
    idx = {}
    for i, g in enumerate(groups):
        idx.setdefault(g, []).append(i)
    xp, hits = x.copy(), 0
    for _ in range(N_PERM):
        for ii in idx.values():
            xp[ii] = rng.permutation(x[ii])
        r = spearmanr(xp, y)[0]
        hits += (abs(r) >= abs(rho)) if two_sided else (r >= rho)
    return rho, (hits + 1) / (N_PERM + 1)


def frac(ok):
    return float(np.mean(ok)) if ok else None


# ---------------------------------------------------------------- band B rows
def band_b(M, tf, shape, temp, prec):
    rows, blocks = [], {}
    for p in sorted(gcov.cache_dir(BBOX, "B").glob("NISAR_*.npz")):
        w = gcov.Window.load(p)
        L = on_grid(w, tf, shape, epsg=EPSG)
        t = utc(p.stem)
        r = {"name": p.stem, "date": t.date().isoformat(), "t": t.isoformat(), "track": int(p.stem.split("_")[5]),
             "dir": p.stem.split("_")[6], "release": p.stem.split("_")[13], "state": state_at(t, temp),
             "rain72": rain_before(t, prec, 72)}
        for c, m in M.items():
            for pol in ("hh", "hv"):
                v = L[pol][m]
                v = v[np.isfinite(v)]
                r[f"{c}_{pol}"] = float(np.median(db(v))) if v.size >= MIN_PX else None
        rows.append(r)
        blocks[p.stem] = {k: L[k] for k in ("hh", "hv")}
    return rows, blocks


# ---------------------------------------------------------------- band A rows (SUB, SUB2)
def band_a(M, tf, temp):
    rows = {}
    for tag, box in (("SUB", SUB), ("SUB2", SUB2)):
        for p in sorted(gcov.cache_dir(box, "A", offdiag=True).glob("NISAR_*.npz")):
            w = gcov.Window.load(p)
            if "hhvv" not in w.extra or "vv" not in w.extra:
                continue
            dx, dy = w.x[1] - w.x[0], w.y[1] - w.y[0]
            atf = rasterio.Affine(dx, 0, w.x[0] - dx / 2, 0, dy, w.y[0] - dy / 2)
            terms = {"hhhh": w.hh, "hvhv": w.hv, "vhvh": w.extra["vh"], "vvvv": w.extra["vv"],
                     **{k: w.extra[k] for k in TERMS}}
            ok = ((w.mask > 0) & (w.mask < 255) & (w.hh > 0) & (w.hv > 0) & (terms["vhvh"] > 0) & (terms["vvvv"] > 0))
            t = utc(p.stem)
            key = p.stem
            r = rows.setdefault(key, {"name": key, "date": t.date().isoformat(), "t": t.isoformat(),
                                      "track": int(key.split("_")[5]), "dir": key.split("_")[6],
                                      "release": key.split("_")[13], "state": state_at(t, temp), "sums": {}})
            for c in ("corn", "soy", "swamp", "upland", "forest"):
                m = to_grid(M[c], tf, w.epsg, atf, w.hh.shape) & ok
                if not m.any():
                    continue
                if tag == "SUB2" and c in ("corn", "soy"):
                    r[f"{c}_hh_vv_db"] = float(np.median(db(w.hh[m] / terms["vvvv"][m]))) if m.sum() >= MIN_PX else None
                    r[f"{c}_hv"] = float(np.median(db(w.hv[m]))) if m.sum() >= MIN_PX else None
                    r[f"{c}_n"] = int(m.sum())
                if tag == "SUB" and c == "forest" and m.sum() >= MIN_PX:
                    r["omega_sub"] = faraday_deg({k: v[m].astype(complex).mean() for k, v in terms.items()})
                s = r["sums"].setdefault(c, {"n": 0, **{k: 0j for k in terms}})
                s["n"] += int(m.sum())
                for k, v in terms.items():
                    s[k] += v[m].astype(complex).sum()
            r.setdefault("epsg", {})[tag] = w.epsg
    out = []
    for r in rows.values():
        for c, s in r.pop("sums").items():
            if s["n"] >= MIN_PX:
                mean = {k: v / s["n"] for k, v in s.items() if k != "n"}
                r[f"{c}_phi_deg"] = float(np.degrees(np.angle(mean["hhvv"])))
                r[f"{c}_H"], r[f"{c}_alpha"] = h_alpha(coherency(mean))
                r[f"{c}_nA"] = s["n"]
        out.append(r)
    return sorted(out, key=lambda r: r["t"])


# ---------------------------------------------------------------- tests on band A / band B
def ws1(arows):
    sel = [r for r in arows if "2025-10-22" <= r["date"] <= "2025-11-30" and r["state"] == "unfrozen"
           and r.get("corn_hh_vv_db") is not None and r.get("soy_hh_vv_db") is not None]
    d1 = [r["corn_hh_vv_db"] - r["soy_hh_vv_db"] for r in sel]
    d2 = [r["corn_hv"] - r["soy_hv"] for r in sel]
    later = [{"date": r["date"], "track": r["track"], "d_hh_vv": r["corn_hh_vv_db"] - r["soy_hh_vv_db"],
              "d_hv": r["corn_hv"] - r["soy_hv"]} for r in arows if r["date"] >= "2026-06-01"
             and r.get("corn_hh_vv_db") is not None and r.get("soy_hh_vv_db") is not None]
    v = lambda d: "inconclusive" if len(d) < 3 else ("pass" if np.mean([x >= 1 for x in d]) >= 2 / 3 else "fail")  # noqa: E731
    return {"n": len(sel), "dates": [r["date"] + f" {r['track']}" for r in sel],
            "corn_minus_soy_hh_vv_db": [round(x, 2) for x in d1], "corn_minus_soy_hv_db": [round(x, 2) for x in d2],
            "WS1a": v(d1), "WS1b": v(d2), "reported_2026": later}


def ws2(rows, arows):
    a = [r["swamp_hh"] - r["upland_hh"] for r in rows if r["state"] == "unfrozen"
         and r.get("swamp_hh") is not None and r.get("upland_hh") is not None]
    b = [(r["date"], abs(r["swamp_phi_deg"]) - abs(r["upland_phi_deg"])) for r in arows if r["state"] == "unfrozen"
         and r.get("swamp_phi_deg") is not None and r.get("upland_phi_deg") is not None]
    v = lambda d, thr: "inconclusive" if len(d) < 5 else ("pass" if np.mean([x >= thr for x in d]) >= 2 / 3 else "fail")  # noqa: E731
    season = lambda lo, hi, xs: [x for d, x in xs if lo <= d[5:] <= hi]  # noqa: E731
    bb = [x for _, x in b]
    phis = [(r["date"], r.get("swamp_phi_deg"), r.get("upland_phi_deg"), r.get("swamp_alpha"), r.get("upland_alpha"))
            for r in arows if r["state"] == "unfrozen"]
    return {"WS2a": {"n": len(a), "median_db": float(np.median(a)) if a else None,
                     "share_ge_1db": frac([x >= 1 for x in a]), "verdict": v(a, 1.0)},
            "WS2b": {"n": len(bb), "median_deg": float(np.median(bb)) if bb else None,
                     "share_ge_10deg": frac([x >= 10 for x in bb]), "verdict": v(bb, 10.0),
                     "leaf_off_Oct_Nov": season("10-01", "11-30", b), "leaf_on_Jun_Sep": season("06-01", "09-30", b),
                     "phi_swamp_upland_alpha": phis}}


def ws3(rows):
    days = {}
    for f in sorted((RAW / "s2").glob("lake_*.npz")):
        days[f.stem[5:]] = np.load(f)["scl"]
    M, _, _ = masks80()
    lake = M["lake"]
    s2 = []
    for d, scl in days.items():
        v = scl[lake]
        clear = np.isin(v, (6, 11))
        if clear.sum() < max(MIN_PX, 0.2 * lake.sum()):
            continue
        s2.append({"date": d, "clear": int(clear.sum()), "ice": float(np.mean(v[clear] == 11)),
                   "water": float(np.mean(v[clear] == 6))})
    ice = next((s["date"] for s in s2 if s["date"] >= "2025-11-01" and s["ice"] >= 0.5), None)
    open_ = None
    if ice:
        before = [s["date"] for s in s2 if s["date"] < ice and s["water"] >= 0.5]
        open_ = before[-1] if before else None
    res = {"lake_px": int(lake.sum()), "s2": s2, "ice_date": ice, "open_date": open_}
    if not ice or not open_:
        res["verdict"] = "inconclusive"
        return res
    per = {}
    for r in rows:
        if r.get("lake_hh") is None or not r["date"] >= "2025-10-01" or r["date"] > "2026-02-28":
            continue
        tag = "ice" if r["date"] > ice else "open" if r["date"] < open_ else None
        if tag:
            per.setdefault(r["track"], {"ice": [], "open": []})[tag].append(r["lake_hh"])
    d = {t: float(np.median(v["ice"]) - np.median(v["open"])) for t, v in per.items() if v["ice"] and v["open"]}
    res["per_track_db"] = d
    res["series"] = [(r["date"], r["track"], r.get("lake_hh"), r.get("lake_hv")) for r in rows
                     if "2025-10-01" <= r["date"] <= "2026-02-28"]
    res["verdict"] = ("inconclusive" if len(d) < 2 else
                      "pass" if np.mean([x >= 3 for x in d.values()]) >= 0.5 else "fail")
    return res


def ws4(rows, blocks, M, tf, shape):
    B = 12
    H, W = (shape[0] // B) * B, (shape[1] // B) * B
    cropm = (M["corn"] | M["soy"])[:H, :W]

    def blk(a, m):
        a = np.where(m, a[:H, :W], np.nan).reshape(H // B, B, W // B, B).transpose(0, 2, 1, 3).reshape(H // B, W // B, -1)
        with np.errstate(all="ignore"):
            import warnings
            warnings.simplefilter("ignore", RuntimeWarning)
            return np.nanmedian(a, -1), np.isfinite(a).sum(-1)
    _, ncrop = blk(np.ones(shape), cropm)
    ok = ncrop >= 50
    # optical series
    opt = []
    for f in sorted((RAW / "s2").glob("ndvi_*.npz")):
        z = np.load(f)
        clear = np.isin(z["scl"], (4, 5))
        with np.errstate(invalid="ignore", divide="ignore"):
            nd = np.where(clear, (z["nir"] - z["red"]) / (z["nir"] + z["red"]), np.nan)
        med, n = blk(nd, cropm)
        opt.append((f.stem[5:], np.where(n >= 0.5 * ncrop, med, np.nan)))
    odates = [d for d, _ in opt]
    O = np.array([a for _, a in opt])                                    # (dates, by, bx)
    o_date = np.full(ok.shape, None, object)
    for i, j in zip(*np.nonzero(ok)):
        s = O[:, i, j]
        if not np.isfinite(s).any() or np.nanmax(s) < 0.7:
            continue
        k = int(np.nanargmax(s))
        after = [q for q in range(k + 1, len(s)) if np.isfinite(s[q]) and s[q] < 0.4]
        if after:
            o_date[i, j] = odates[after[0]]
    # radar series per track (2026)
    pairs, rmap = [], {}
    for t in sorted({r["track"] for r in rows}):
        rr = sorted([r for r in rows if r["track"] == t and r["date"] >= "2026-06-01"], key=lambda r: r["t"])
        if len(rr) < 4:
            continue
        S = np.array([blk(db(blocks[r["name"]]["hv"]), cropm)[0] for r in rr])
        dates = [r["date"] for r in rr]
        rd = np.full(ok.shape, None, object)
        for i, j in zip(*np.nonzero(ok)):
            s = S[:, i, j]
            if not np.isfinite(s).sum() >= 4:
                continue
            k = int(np.nanargmax(s))
            mx = s[k]
            for q in range(k + 1, len(s) - 1):
                if s[q] <= mx - 3 and s[q + 1] <= mx - 3:
                    rd[i, j] = dates[q]
                    break
            if rd[i, j] and o_date[i, j]:
                lag = (datetime.fromisoformat(rd[i, j]) - datetime.fromisoformat(o_date[i, j])).days
                pairs.append({"track": t, "i": int(i), "j": int(j), "radar": rd[i, j], "optical": o_date[i, j],
                              "lag_days": lag})
        rmap[t] = rd
    lags = [p["lag_days"] for p in pairs]
    res = {"blocks": int(ok.sum()), "optical_dates": odates, "n_optical": int(sum(x is not None for x in o_date.ravel())),
           "n_pairs": len(pairs), "share_within_24d": frac([abs(x) <= 24 for x in lags]),
           "median_lag_days": float(np.median(lags)) if lags else None, "pairs": pairs,
           "radar_dates_per_track": {t: int(sum(x is not None for x in v.ravel())) for t, v in rmap.items()}}
    res["verdict"] = ("inconclusive" if len(pairs) < 30 else "pass" if res["share_within_24d"] >= 0.6 else "fail")
    # block centres for the map
    xs = tf.c + tf.a * (np.arange(W // B) * B + B / 2)
    ys = tf.f + tf.e * (np.arange(H // B) * B + B / 2)
    XX, YY = np.meshgrid(xs, ys)
    lon, lat = warp_transform(f"EPSG:{EPSG}", "EPSG:4326", XX.ravel().tolist(), YY.ravel().tolist())
    res["_map"] = {"lon": np.array(lon).reshape(XX.shape), "lat": np.array(lat).reshape(XX.shape), "ok": ok,
                   "radar": rmap, "optical": o_date}
    return res


def ws5(rows):
    sel = [r for r in rows if r["state"] == "unfrozen" and r.get("crop_hh") is not None and
           ("2025-10-15" <= r["date"] <= "2025-11-30" or "2026-10-01" <= r["date"] <= "2026-11-10")]
    g = {}
    for r in sel:
        g.setdefault((r["track"], r["date"][:4]), []).append(r)
    xs, ys, gs = [], [], []
    for k, rr in g.items():
        if len(rr) < 2:
            continue
        m = np.mean([r["crop_hh"] for r in rr])
        for r in rr:
            xs.append(r["rain72"])
            ys.append(r["crop_hh"] - m)
            gs.append(k)
    if len(xs) < 12:
        return {"n": len(xs), "verdict": "inconclusive"}
    rho, p = perm_rho(xs, ys, gs)
    return {"n": len(xs), "rho": rho, "p": p, "verdict": "pass" if rho >= 0.4 and p <= 0.05 else "fail",
            "points": [[round(a, 1), round(b, 2)] for a, b in zip(xs, ys)]}


# ---------------------------------------------------------------- GUNW
def gunw_rows(M, tf, temp):
    rows = []
    for f in sorted((RAW / "gunw").glob("*.npz")):
        raw = nisar.RawWindow.load(f)
        ref, sec = nisar._parse_dates(f.stem)
        dx, dy = raw.x[1] - raw.x[0], raw.y[1] - raw.y[0]
        gtf = rasterio.Affine(dx, 0, raw.x[0] - dx / 2, 0, dy, raw.y[0] - dy / 2)
        r = {"name": f.stem, "track": int(f.stem.split("_")[5]), "dir": f.stem.split("_")[6],
             "ref": ref.isoformat(), "sec": sec.isoformat(), "days": round(raw.days),
             "state": (state_at(ref, temp), state_at(sec, temp))}
        valid = np.isfinite(raw.phase) & (raw.comp > 0)
        main = np.bincount(raw.comp[valid].ravel()).argmax() if valid.any() else 0
        los = -raw.wavelength / (4 * np.pi) * (raw.phase - (raw.iono if raw.iono is not None else 0)) * 1000  # mm, + = towards
        for c in ("upland", "crop", "forest"):
            m = to_grid(M[c], tf, raw.epsg, gtf, raw.phase.shape)
            coh = raw.coherence[m & np.isfinite(raw.coherence)]
            r[f"{c}_coh"] = float(np.median(coh)) if coh.size >= MIN_PX else None
            mm = m & valid & (raw.comp == main) & np.isfinite(los)
            r[f"{c}_los_mm"] = float(np.median(los[mm])) if mm.sum() >= MIN_PX else None
        # iono screen median over SUB
        bx, by = warp_transform("EPSG:4326", f"EPSG:{raw.epsg}", [SUB[0], SUB[2]], [SUB[1], SUB[3]])
        yy, xx = np.meshgrid(raw.y, raw.x, indexing="ij")
        sub = (xx >= min(bx)) & (xx <= max(bx)) & (yy >= min(by)) & (yy <= max(by))
        if raw.iono is not None:
            v = raw.iono[sub & np.isfinite(raw.iono)]
            r["iono_sub_rad"] = float(np.median(v)) if v.size else None
            w = raw.iono[np.isfinite(raw.iono)]
            r["iono_window_rad"] = float(np.median(w)) if w.size else None
        r["wavelength"] = float(raw.wavelength)
        rows.append(r)
    return rows


def ws6(g):
    res, cases = {}, []
    for t in sorted({r["track"] for r in g}):
        for c in ("upland", "crop"):
            fz = [r[f"{c}_coh"] for r in g if r["track"] == t and r["days"] == 12 and r["state"] == ("frozen", "frozen")
                  and r.get(f"{c}_coh") is not None]
            uf = [r[f"{c}_coh"] for r in g if r["track"] == t and r["days"] == 12
                  and r["state"] == ("unfrozen", "unfrozen") and r.get(f"{c}_coh") is not None]
            if fz and len(uf) >= 2:
                d = float(np.median(fz) - np.median(uf))
                res[f"{t}_{c}"] = {"frozen": [round(x, 3) for x in fz], "unfrozen_median": round(float(np.median(uf)), 3),
                                   "diff": round(d, 3)}
                cases.append(d)
    res["verdict"] = ("inconclusive" if len(cases) < 2 else
                      "pass" if np.mean([x >= 0.10 for x in cases]) >= 2 / 3 else "fail")
    return res


def ws7(g):
    res, sums = {}, []
    for t in sorted({r["track"] for r in g}):
        rr = [r for r in g if r["track"] == t and r["days"] == 12 and "2025-11-01" <= r["ref"][:10]
              and r["sec"][:10] <= "2026-01-31" and r.get("crop_los_mm") is not None and r.get("upland_los_mm") is not None]
        if not rr:
            continue
        d = [r["crop_los_mm"] - r["upland_los_mm"] for r in rr]
        res[t] = {"n_pairs": len(rr), "per_pair_mm": [round(x, 1) for x in d], "cum_mm": round(float(np.sum(d)), 1)}
        if len(rr) >= 4:
            sums.append(float(np.sum(d)))
    ok = len(sums) >= 2
    res["verdict"] = ("inconclusive" if not ok else
                      "pass" if sum(x >= 5 for x in sums) >= 2 and not any(x <= -5 for x in sums) else "fail")
    return res


def ws8(g, arows):
    om = {(r["track"], r["date"]): r["omega_sub"] for r in arows if r.get("omega_sub") is not None}
    pts = []
    for r in g:
        a, b = (r["track"], r["ref"][:10]), (r["track"], r["sec"][:10])
        if a in om and b in om and r.get("iono_sub_rad") is not None:
            pts.append({"pair": r["name"][:0] + f"{r['track']} {a[1]}→{b[1]}", "d_omega": om[b] - om[a],
                        "iono": r["iono_sub_rad"], "iono_window": r.get("iono_window_rad")})
    res = {"n": len(pts), "points": pts}
    if all((p["iono_window"] or 0) == 0 for p in pts) and pts:
        res["verdict"] = "inconclusive (screen median 0)"
        return res
    if len(pts) < 8:
        res["verdict"] = "inconclusive"
        return res
    rho, p = perm_rho([q["iono"] for q in pts], [q["d_omega"] for q in pts], [0] * len(pts), two_sided=True)
    res.update({"rho": rho, "p": p, "verdict": "pass" if abs(rho) >= 0.5 and p <= 0.05 else "fail",
                "slope_deg_per_rad": float(np.polyfit([q["iono"] for q in pts], [q["d_omega"] for q in pts], 1)[0])})
    return res


# ---------------------------------------------------------------- WS9: Faraday against CODE GIM
def read_ionex(path):
    lines = Path(path).read_text(errors="ignore").splitlines()
    exp, maps, cur, lat, row = -1, {}, None, None, []
    i = 0
    while i < len(lines):
        L = lines[i]
        lab = L[60:].strip()
        if lab == "EXPONENT":
            exp = int(L[:6])
        elif lab == "START OF TEC MAP":
            cur = []
        elif lab == "EPOCH OF CURRENT MAP" and cur is not None:
            y, mo, d, h, mi, s = (int(v) for v in L[:36].split())
            ep = datetime(y, mo, d) + timedelta(hours=h, minutes=mi, seconds=s)
        elif lab == "LAT/LON1/LON2/DLON/H" and cur is not None:
            lat = float(L[2:8])
            row, i = [], i + 1
            while len(row) < 73:
                row += [int(L2) for L2 in lines[i].split()]
                i += 1
            cur.append((lat, row))
            continue
        elif lab == "END OF TEC MAP" and cur is not None:
            maps[ep] = (np.array([a for a, _ in cur]), np.array([b for _, b in cur], float) * 10.0 ** exp)
            cur = None
        elif lab == "START OF RMS MAP":
            break
        i += 1
    return maps                                                   # epoch -> (lats, TECU[lat, lon -180..180 step 5])


def vtec(maps, t, lat, lon):
    eps = sorted(maps)
    k = max(0, min(len(eps) - 2, int(np.searchsorted(eps, t)) - 1))
    a, b = eps[k], eps[k + 1]
    wgt = (t - a).total_seconds() / (b - a).total_seconds()
    out = []
    for e in (a, b):
        lats, M = maps[e]
        lons = np.arange(-180, 180.1, 5.0)
        li = np.interp(lat, lats[::-1], np.arange(len(lats))[::-1])
        lj = np.interp(lon, lons, np.arange(len(lons)))
        i0, j0 = int(li), int(lj)
        fi, fj = li - i0, lj - j0
        v = (M[i0, j0] * (1 - fi) * (1 - fj) + M[min(i0 + 1, len(lats) - 1), j0] * fi * (1 - fj)
             + M[i0, j0 + 1] * (1 - fi) * fj + M[min(i0 + 1, len(lats) - 1), j0 + 1] * fi * fj)
        out.append(v)
    return out[0] * (1 - wgt) + out[1] * wgt


def track_geometry():
    """Per track: LOS unit vector (target -> sensor, E/N/U) and band A centre frequency at the SUB centre."""
    f = RES / "geometry.json"
    if f.exists():
        return {int(k): v for k, v in json.loads(f.read_text()).items()}
    lon, lat = (SUB[0] + SUB[2]) / 2, (SUB[1] + SUB[3]) / 2
    pick = {}
    for s in gcov.search(lon, lat):
        if s.path in TRACKS and s.pols == "QPDH":
            pick[s.path] = s
    out = {}
    g = "science/LSAR/GCOV/metadata/radarGrid"
    for t, s in pick.items():
        with nisar.open_h5(s.url) as h5:
            hz = h5[f"{g}/heightAboveEllipsoid"][()]
            k = int(np.argmin(np.abs(hz - HEIGHT)))
            cx, cy = h5[f"{g}/xCoordinates"][()], h5[f"{g}/yCoordinates"][()]
            epsg = int(h5["science/LSAR/GCOV/grids/frequencyA/projection"][()])
            (x,), (y,) = warp_transform("EPSG:4326", f"EPSG:{epsg}", [lon], [lat])
            i, j = int(np.argmin(np.abs(cy - y))), int(np.argmin(np.abs(cx - x)))
            ux, uy = float(h5[f"{g}/losUnitVectorX"][k, i, j]), float(h5[f"{g}/losUnitVectorY"][k, i, j])
            fp = nisar._find(h5, "centerFrequency", "frequencyA")
            fr = float(np.atleast_1d(h5[fp][()])[0])
        out[t] = {"los_en": [ux, uy], "f_hz": fr}
    RES.mkdir(parents=True, exist_ok=True)
    f.write_text(json.dumps(out, indent=1))
    return out


def omega_pred(t_utc, geo, maps_by_day):
    R, Hs = 6371.0, 450.0
    lon0, lat0 = np.radians((SUB[0] + SUB[2]) / 2), np.radians((SUB[1] + SUB[3]) / 2)
    ux, uy = geo["los_en"]
    uz = np.sqrt(max(0.0, 1 - ux * ux - uy * uy))
    e = np.arcsin(uz)
    A = np.arctan2(ux, uy)
    psi = np.pi / 2 - e - np.arcsin(R / (R + Hs) * np.cos(e))
    lat = np.arcsin(np.sin(lat0) * np.cos(psi) + np.cos(lat0) * np.sin(psi) * np.cos(A))
    lon = lon0 + np.arcsin(np.sin(psi) * np.sin(A) / np.cos(lat))
    F = 1 / np.sqrt(1 - (R * np.cos(e) / (R + Hs)) ** 2)
    maps = maps_by_day.get(t_utc.date())
    if maps is None:
        return None, None
    tecu = vtec(maps, t_utc, np.degrees(lat), np.degrees(lon))
    stec = tecu * F * 1e16
    Bt, I, D = 45400e-9, np.radians(69.5), np.radians(-1.5)
    B = np.array([Bt * np.cos(I) * np.sin(D), Bt * np.cos(I) * np.cos(D), -Bt * np.sin(I)])
    bpar = float(B @ np.array([ux, uy, uz]))
    om = 2.365e4 * bpar * stec / geo["f_hz"] ** 2
    return float(np.degrees(om)), float(tecu)


def ws9(arows):
    geo = track_geometry()
    maps = {}
    for f in sorted((RAW / "gim").glob("*.INX")):
        y, doy = int(f.name[11:15]), int(f.name[15:18])
        maps[(datetime(y, 1, 1) + timedelta(days=doy - 1)).date()] = read_ionex(f)
    pts = []
    for r in arows:
        if r.get("omega_sub") is None or r["track"] not in geo:
            continue
        p, tecu = omega_pred(datetime.fromisoformat(r["t"]), geo[r["track"]], maps)
        if p is not None:
            pts.append({"date": r["date"], "track": r["track"], "release": r["release"], "omega": r["omega_sub"],
                        "omega_pred": p, "vtec_tecu": tecu})
    res = {"n": len(pts), "points": pts, "geometry": geo}
    if len(pts) < 15:
        res["WS9a"] = "inconclusive"
    else:
        rho, p = perm_rho([abs(q["omega_pred"]) for q in pts], [abs(q["omega"]) for q in pts], [0] * len(pts))
        res.update({"rho_abs": rho, "p": p, "WS9a": "pass" if rho >= 0.5 and p <= 0.05 else "fail"})
    big = [q for q in pts if abs(q["omega_pred"]) >= 1]
    res["sign_share"] = frac([np.sign(q["omega"]) == np.sign(q["omega_pred"]) for q in big])
    res["n_big"] = len(big)
    res["WS9b"] = "inconclusive" if not big else ("pass" if res["sign_share"] >= 0.8 else "fail")
    if pts:
        res["slope"] = float(np.polyfit([q["omega_pred"] for q in pts], [q["omega"] for q in pts], 1)[0])
    return res


# ---------------------------------------------------------------- WS10 maps + story data
def ws10_maps(M, tf):
    out = {}
    for tag, box in (("SUB", SUB), ("SUB2", SUB2)):
        first = {}
        for p in sorted(gcov.cache_dir(box, "A", offdiag=True).glob("NISAR_*.npz")):
            ym = p.stem.split("_")[11][:6]
            first.setdefault(ym, p)
        panels = []
        for ym, p in sorted(first.items()):
            w = gcov.Window.load(p)
            if "hhvv" not in w.extra:
                continue
            c = {"hhhh": w.hh, "hvhv": w.hv, "vhvh": w.extra["vh"], "vvvv": w.extra["vv"], **{k: w.extra[k] for k in TERMS}}
            sm = {k: (uniform_filter(np.nan_to_num(v.real), 5) + 1j * uniform_filter(np.nan_to_num(v.imag), 5)
                      if np.iscomplexobj(v) else uniform_filter(np.nan_to_num(v.astype(float)), 5)) for k, v in c.items()}
            T = coherency(sm)                                                  # (3, 3, ny, nx)
            Tm = np.moveaxis(T, (0, 1), (-2, -1))
            lam, vec = np.linalg.eigh(Tm)
            lam = np.clip(lam.real, 1e-12, None)
            pr = lam / lam.sum(-1, keepdims=True)
            Hh = -(pr * np.log(pr) / np.log(3)).sum(-1)
            al = np.degrees((pr * np.arccos(np.clip(np.abs(vec[..., 0, :]), 0, 1))).sum(-1))
            pauli = np.stack([sm["hhhh"] + sm["vvvv"] - 2 * sm["hhvv"].real, 2 * sm["hvhv"],
                              sm["hhhh"] + sm["vvvv"] + 2 * sm["hhvv"].real], -1) / 2
            bad = ~((w.mask > 0) & (w.mask < 255))
            Hh[bad], al[bad] = np.nan, np.nan
            panels.append((ym, Hh, al, pauli))
        if not panels:
            continue
        n = len(panels)
        fig, ax = plt.subplots(3, n, figsize=(2.1 * n, 6.6))
        ax = np.atleast_2d(ax).reshape(3, n)
        for k, (ym, Hh, al, pa) in enumerate(panels):
            ax[0, k].imshow(Hh, vmin=0.3, vmax=1, cmap="magma")
            ax[1, k].imshow(al, vmin=10, vmax=60, cmap="viridis")
            rgb = np.clip((10 * np.log10(np.clip(pa, 1e-6, None)) + 25) / 25, 0, 1)
            ax[2, k].imshow(rgb)
            ax[0, k].set_title(f"{ym[:4]}-{ym[4:]}", fontsize=8)
            for a in ax[:, k]:
                a.set_xticks([]), a.set_yticks([])
        ax[0, 0].set_ylabel("entropy H"), ax[1, 0].set_ylabel("ᾱ, deg"), ax[2, 0].set_ylabel("Pauli RGB")
        fig.suptitle(f"WS10 {tag}: per-pixel polarimetry, first scene of each month (5×5 boxcar)", fontsize=9)
        fig.tight_layout()
        fig.savefig(RES / f"ws10_{tag}.png", dpi=110)
        plt.close(fig)
        out[tag] = [ym for ym, *_ in panels]
    return out


def story_cells(rows, blocks, M, tf, shape):
    """1 km blocks: dominant class and HH / HV per date relative to the block's own Jun-Sep 2026 median on that track."""
    B = 12
    H, W = (shape[0] // B) * B, (shape[1] // B) * B
    rr = sorted(rows, key=lambda r: r["t"])

    def blk(a):
        import warnings
        warnings.simplefilter("ignore", RuntimeWarning)
        return np.nanmedian(a[:H, :W].reshape(H // B, B, W // B, B).transpose(0, 2, 1, 3).reshape(H // B, W // B, -1), -1)
    S = {pol: np.array([blk(db(blocks[r["name"]][pol])) for r in rr]) for pol in ("hh", "hv")}
    anom = {}
    for pol, A in S.items():
        out = np.full_like(A, np.nan)
        for t in {r["track"] for r in rr}:
            idx = [i for i, r in enumerate(rr) if r["track"] == t]
            summer = [i for i in idx if "2026-06-01" <= rr[i]["date"] <= "2026-09-30"]
            if len(summer) >= 3:
                out[idx] = A[idx] - np.nanmedian(A[summer], 0)
        anom[pol] = out
    names = ["lake", "swamp", "upland", "corn", "soy", "crop"]
    frac = np.array([blk(M[n].astype(float)) for n in names])
    dom = np.argmax(frac, 0)
    xs = tf.c + tf.a * (np.arange(W // B) * B + B / 2)
    ys = tf.f + tf.e * (np.arange(H // B) * B + B / 2)
    XX, YY = np.meshgrid(xs, ys)
    lon, lat = warp_transform(f"EPSG:{EPSG}", "EPSG:4326", XX.ravel().tolist(), YY.ravel().tolist())
    cells = []
    for k in range(XX.size):
        i, j = divmod(k, XX.shape[1])
        hv = anom["hv"][:, i, j]
        if not np.isfinite(hv).any():
            continue
        cells.append([round(lon[k], 5), round(lat[k], 5), names[dom[i, j]] if frac[dom[i, j], i, j] >= 0.3 else "mixed",
                      [None if not np.isfinite(v) else round(float(v), 2) for v in anom["hh"][:, i, j]],
                      [None if not np.isfinite(v) else round(float(v), 2) for v in hv]])
    return {"cell_km": B * abs(tf.a) / 1000, "columns": ["lon", "lat", "class", "hh_anom_db", "hv_anom_db"],
            "dates": [r["date"] for r in rr], "tracks": [r["track"] for r in rr], "states": [r["state"] for r in rr],
            "note": "per-block HH / HV (dB) relative to the block's own Jun-Sep 2026 median on the same track", "cells": cells}


# ---------------------------------------------------------------- main
def main():
    RES.mkdir(parents=True, exist_ok=True)
    temp, prec = weather()
    M, tf, shape = masks80()
    print("class px (80 m):", {k: int(v.sum()) for k, v in M.items()}, flush=True)
    rows, blocks = band_b(M, tf, shape, temp, prec)
    arows = band_a(M, tf, temp)
    g = gunw_rows(M, tf, temp)
    print(f"band B {len(rows)}, band A {len(arows)}, GUNW {len(g)}", flush=True)
    res = {"rows_B": rows, "rows_A": arows, "rows_gunw": g}
    res["WS1"] = ws1(arows)
    res["WS2"] = ws2(rows, arows)
    res["WS3"] = ws3(rows)
    w4 = ws4(rows, blocks, M, tf, shape)
    mp = w4.pop("_map")
    res["WS4"] = w4
    res["WS5"] = ws5(rows)
    res["WS6"] = ws6(g)
    res["WS7"] = ws7(g)
    res["WS8"] = ws8(g, arows)
    res["WS9"] = ws9(arows)
    res["WS10"] = ws10_maps(M, tf)
    (RES / "ws_results.json").write_text(json.dumps(res, indent=1, default=str), encoding="utf-8")
    (RES / "story_cells.json").write_text(json.dumps(story_cells(rows, blocks, M, tf, shape), default=str), encoding="utf-8")
    harvest = [[round(float(mp["lon"][i, j]), 5), round(float(mp["lat"][i, j]), 5),
                {str(t): v[i, j] for t, v in mp["radar"].items()}, mp["optical"][i, j]]
               for i, j in zip(*np.nonzero(mp["ok"]))]
    (RES / "story_harvest.json").write_text(json.dumps({"columns": ["lon", "lat", "radar_by_track", "optical"],
                                                         "cells": harvest}, default=str), encoding="utf-8")
    for k in ("WS1", "WS2", "WS3", "WS4", "WS5", "WS6", "WS7", "WS8", "WS9"):
        v = {a: b for a, b in res[k].items() if a not in ("points", "pairs", "series", "s2", "phi_swamp_upland_alpha",
                                                          "reported_2026", "geometry")} if isinstance(res[k], dict) else res[k]
        print(k, json.dumps(v, default=str)[:700], flush=True)


if __name__ == "__main__":
    main()
