"""Pre-registered tests L1, L2, C1, (C2), H: old lakebed and ground cracks vs the NISAR rate (Mexico City).

    uv run python scripts/test_lake_cracks.py sites/mexico-city-asc.json

Pre-registration: notes/experiments.md, section "L / C". Thresholds are fixed there; this script only reports.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from rasterio.features import rasterize
from rasterio.warp import transform_geom
from scipy.ndimage import map_coordinates, median_filter
from scipy.stats import spearmanr

from hachure import engine, gridding, nisar

ROOT = Path(__file__).resolve().parents[1]
ATLAS = ROOT / "data/raw/cdmx_atlas"
ZONE = {"Zona I": 1, "Zona II": 2, "Zona III": 3}          # Lomas, Transición, Lago (mapping fixed in the pre-registration)
NAMES = {1: "Lomas", 2: "Transición", 3: "Lago"}
NO_TROPO = "--no-tropo" in sys.argv


def load(name):
    return json.loads((ATLAS / name).read_bytes().decode("utf-8", errors="replace"))


def fields(cfg):
    v = engine.load_velocity(ROOT / "data/processed" / cfg["id"] / "engine_velocity.npz")
    if not NO_TROPO:
        v = engine.tropo_correct(cfg, v)
    z, T, crs = engine._dem(cfg)
    res, shape, p = cfg["dem_res"], z.shape, cfg["velocity"]
    src_t, src_crs = gridding.transform_from_centres(v.x, v.y), f"EPSG:{v.epsg}"

    def grid(a):
        return gridding.to_grid(np.asarray(a, float), src_t, src_crs, T, crs, shape)

    years = (v.dates[-1] - v.dates[0]).days / 365.25
    coh = np.nan_to_num(grid(v.coherence)) >= p.get("coh_min_draw", 0.45)
    r_raw = grid(v.rate * 100)
    r_raw[~(coh & np.isfinite(r_raw))] = np.nan
    # R_s exactly as velocity_layer builds the final frame (smoothed 400 m, deramped), as a rate
    sig = 400 / res
    g = grid(v.cum[-1] - v.cum[0])
    ok = coh & np.isfinite(g)
    g = median_filter(np.where(ok, g, np.nanmedian(g[ok])), size=5)
    s, w = engine._smooth_nan(g, sig, ok)
    if p.get("deramp", True):
        s = s - nisar.fit_ramp(np.where(ok & (w > 0.3), s, np.nan), cut=-max(0.005, 0.03 * years))
    support = engine._smooth_nan(np.ones(shape), sig, coh)[1] > 0.3
    valid_s = coh & support
    r_s = np.where(valid_s, s / years * 100, np.nan)
    gy, gx = np.gradient(r_s)
    G = np.hypot(gx, gy) / (res / 1000)                       # cm/yr per km
    return dict(T=T, crs=crs, shape=shape, res=res, r_raw=r_raw, r_s=r_s, G=G, valid_s=valid_s)


def burn(fc, src_crs, F, value):
    shapes = []
    for f in fc["features"]:
        if f["geometry"] is None:
            continue
        val = value(f["properties"])
        if val:
            shapes.append((transform_geom(src_crs, F["crs"], f["geometry"]), val))
    return rasterize(shapes, out_shape=F["shape"], transform=F["T"], fill=0, dtype="float32")


def crack_points(F, step_m=100.0):
    fc = load("Fracturas.geojson")
    pts = []
    for f in fc["features"]:
        geom = transform_geom("EPSG:3857", F["crs"], f["geometry"])
        lines = geom["coordinates"] if geom["type"] == "MultiLineString" else [geom["coordinates"]]
        for ln in lines:
            a = np.asarray(ln, float)
            if len(a) < 2:
                continue
            L = np.concatenate([[0], np.cumsum(np.hypot(*np.diff(a, axis=0).T))])
            t = np.arange(0, L[-1] + 1e-9, step_m) if L[-1] > step_m else np.array([L[-1] / 2])
            pts.append(np.column_stack([np.interp(t, L, a[:, 0]), np.interp(t, L, a[:, 1])]))
    xy = np.vstack(pts)
    T = F["T"]
    col = (xy[:, 0] - T.c) / T.a - 0.5
    row = (xy[:, 1] - T.f) / T.e - 0.5
    return np.column_stack([row, col]), len(fc["features"])


def sample(a, rc):
    H, W = a.shape
    r, c = np.round(rc[:, 0]).astype(int), np.round(rc[:, 1]).astype(int)
    inside = (r >= 0) & (r < H) & (c >= 0) & (c < W)
    out = np.full(len(rc), np.nan)
    out[inside] = a[r[inside], c[inside]]
    return out


def main():
    cfg = engine.load_site(sys.argv[1])
    F = fields(cfg)
    zones = burn(load("SPC_Zonificacion_Geotecnica.geojson"), "EPSG:6365", F,
                 lambda pr: ZONE.get(str(pr.get("INTENSIDAD")).strip(), 0))
    r_raw, r_s, G = F["r_raw"], F["r_s"], F["G"]
    out = {}

    # ---------------- L1, L2
    vr = np.isfinite(r_raw)
    fast = vr & (r_raw <= -10)
    share_lago = float((zones[fast] == 3).mean())
    out["L1"] = dict(fast_px=int(fast.sum()), share_in_lago=round(share_lago, 3),
                     share_no_zone=round(float((zones[fast] == 0).mean()), 3), passed=share_lago >= 0.80)
    med = {NAMES[k]: float(np.median(r_raw[vr & (zones == k)])) if (vr & (zones == k)).any() else None for k in (1, 2, 3)}
    n = {NAMES[k]: int((vr & (zones == k)).sum()) for k in (1, 2, 3)}
    order_ok = med["Lago"] < med["Transición"] < med["Lomas"]
    out["L2"] = dict(median_cm_yr={k: round(x, 2) for k, x in med.items()}, n_px=n,
                     lago_minus_lomas=round(med["Lago"] - med["Lomas"], 2),
                     passed=bool(order_ok and med["Lago"] - med["Lomas"] <= -5))

    # ---------------- C1 (gradient), C2 (rate, exploratory)
    rc, n_lines = crack_points(F)
    Gp, Rp = sample(G, rc), sample(r_raw, rc)
    okp = np.isfinite(Gp) & sample(F["valid_s"].astype(float), rc).astype(bool)
    obs_G = float(np.median(Gp[okp])) if okp.sum() >= 50 else None
    okr = np.isfinite(Rp)
    obs_R = float(np.median(Rp[okr])) if okr.sum() >= 50 else None
    rng = np.random.default_rng(0)
    nullG, nullR = [], []
    allowed = np.isin(zones, (2, 3))
    while len(nullG) < 500:
        dx, dy = rng.uniform(-3000, 3000, 2)
        if np.hypot(dx, dy) < 1000:
            continue
        sh = rc + np.array([-dy / F["res"], dx / F["res"]])        # row grows southwards (T.e < 0)
        keep = sample(allowed.astype(float), sh) == 1
        g = sample(G, sh)
        kg = keep & np.isfinite(g) & (sample(F["valid_s"].astype(float), sh) == 1)
        nullG.append(np.median(g[kg]) if kg.sum() >= 50 else np.nan)
        r = sample(r_raw, sh)
        kr = keep & np.isfinite(r)
        nullR.append(np.median(r[kr]) if kr.sum() >= 50 else np.nan)
    nullG, nullR = np.array(nullG), np.array(nullR)
    p95 = float(np.nanpercentile(nullG, 95))
    out["C1"] = dict(lines=n_lines, points=len(rc), points_valid=int(okp.sum()),
                     zone_of_points={NAMES.get(k, "none"): int(((sample(zones, rc) == k) & okp).sum()) for k in (0, 1, 2, 3)},
                     median_G_cm_yr_per_km=None if obs_G is None else round(obs_G, 3),
                     null_median=round(float(np.nanmedian(nullG)), 3), null_p95=round(p95, 3),
                     share_null_ge_obs=None if obs_G is None else round(float(np.nanmean(nullG >= obs_G)), 4),
                     null_shifts_used=int(np.isfinite(nullG).sum()),
                     passed=None if obs_G is None else bool(obs_G > p95))
    out["C2_exploratory"] = dict(median_rate_cm_yr=None if obs_R is None else round(obs_R, 2),
                                 null_median=round(float(np.nanmedian(nullR)), 2),
                                 share_null_le_obs=None if obs_R is None else round(float(np.nanmean(nullR <= obs_R)), 4))

    # ---------------- H: their subsidence classes vs our smoothed rate
    mid = {"02-10": 6.0, "11-20": 15.5, "21-30": 25.5}
    cls = burn(load("SPC_Hundimientos_ITRF.geojson"), "EPSG:6365", F,
               lambda pr: mid.get(str(pr.get("MAGNI_NUM")).strip(), 0))
    m = (cls > 0) & np.isfinite(r_s)
    rho = spearmanr(cls[m], -r_s[m]).statistic if m.sum() > 10 else np.nan
    out["H"] = dict(px=int(m.sum()), spearman=round(float(rho), 3),
                    median_r_s_by_class={k: round(float(np.median(r_s[m & (cls == v)])), 2) if (m & (cls == v)).any() else None
                                         for k, v in mid.items()},
                    passed=bool(rho >= 0.6))
    print(json.dumps(out, indent=1, ensure_ascii=False))
    dst = ROOT / "data/processed" / cfg["id"] / ("test_lake_cracks.json" if NO_TROPO or not cfg.get("tropo") else "test_lake_cracks_tropo.json")
    dst.write_text(json.dumps(out, indent=1, ensure_ascii=False), encoding="utf-8")


if __name__ == "__main__":
    main()
