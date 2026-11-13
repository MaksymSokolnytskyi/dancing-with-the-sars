"""H4 (leftover ionosphere) vs H5 (evening tropospheric turbulence): three discriminators.

For each pass, residual of every pair = pair - (rate from the other pairs) x span (as atmo_stats).
 1. Pattern test: spatial correlation between a pair's residual and its own ionosphere screen.
    A leftover of an imperfect correction looks like a scaled copy of the screen -> |corr| high.
 2. Scale test: noise vs distance for both passes. Ionosphere lives at long wavelengths
    (tens of km), turbulence at short ones -> where is the evening excess?
 3. Height test: share of the residual explained by terrain height (stratified troposphere).

    uv run python scripts/iono_vs_tropo.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from hachure import dem as demmod, gridding, nisar

ROOT = Path(__file__).resolve().parents[1]
BBOX = (-99.32, 19.17, -98.82, 19.58)
REF_BBOX = (-99.14, 19.53, -99.09, 19.57)
LAGS = np.array([0.5, 1, 2, 5, 10, 20, 30])


def sf(r, rng, n=150_000, dx=0.08):
    h, w = r.shape
    out = []
    for lag in LAGS:
        L = lag / dx
        a = rng.uniform(0, 2 * np.pi, n)
        y0, x0 = rng.integers(0, h, n), rng.integers(0, w, n)
        y1, x1 = np.round(y0 + L * np.sin(a)).astype(int), np.round(x0 + L * np.cos(a)).astype(int)
        ok = (y1 >= 0) & (y1 < h) & (x1 >= 0) & (x1 < w)
        d = r[y0[ok], x0[ok]] - r[y1[ok], x1[ok]]
        d = d[np.isfinite(d)]
        out.append(np.sqrt(0.5 * np.mean(d ** 2)) if d.size > 500 else np.nan)
    return np.array(out)


def plane_free(a):
    yy, xx = np.nonzero(np.isfinite(a))
    P = np.c_[np.ones(len(yy)), xx, yy]
    c, *_ = np.linalg.lstsq(P, a[yy, xx], rcond=None)
    out = a.copy()
    out[yy, xx] -= P @ c
    return out


def run(direction):
    found = nisar.search((BBOX[0] + BBOX[2]) / 2, (BBOX[1] + BBOX[3]) / 2, direction=direction)
    chain = nisar.sequential_chain(nisar.dominant_track(found))
    cache = ROOT / "data" / "raw" / "nisar_cache" / ("bbox_" + "_".join(f"{v:.3f}" for v in BBOX) + "_h2240")
    raws = [nisar.RawWindow.load(cache / f"{g.name}.npz") for g in chain]
    pairs = [nisar.process(r, ref_bbox=REF_BBOX, coh_min=0.3) for r in raws]
    disps, _, x, y = nisar.align(pairs)
    # ionosphere screens on the same aligned grid, converted like the displacement (m, vertical)
    ionos = []
    for r, p in zip(raws, pairs):
        k = nisar.LOS_SIGN * r.wavelength / (4 * np.pi) / np.cos(np.radians(r.incidence))
        ionos.append(nisar.PairWindow(r.iono * k, r.coherence, r.x, r.y, r.epsg, r.days))
    iono_al = nisar.align(ionos)[0]
    days = np.array([p.days for p in pairs])
    zd, Td, crsd = demmod.load(ROOT / "data" / "raw" / "copdem_mexico-city_50m.tif")
    z = gridding.to_grid(zd, Td, crsd, gridding.transform_from_centres(x, y), f"EPSG:{pairs[0].epsg}",
                         disps.shape[1:])
    ok = np.isfinite(disps)
    tot = np.where(ok, disps, 0).sum(0)
    span = np.where(ok, days[:, None, None], 0).sum(0)
    rng = np.random.default_rng(1)
    rows, sfs = [], []
    for i, g in enumerate(chain):
        with np.errstate(all="ignore"):
            rl = (tot - np.where(ok[i], disps[i], 0)) / (span - days[i])
        res = np.where(ok[i] & (span - days[i] > 200), disps[i] - rl * days[i], np.nan)
        res -= np.nanmedian(res)
        io = np.where(np.isfinite(res), iono_al[i], np.nan)
        m = np.isfinite(res) & np.isfinite(io) & np.isfinite(z)
        rp, ip = plane_free(np.where(m, res, np.nan)), plane_free(np.where(m, io, np.nan))
        corr = np.corrcoef(rp[m], ip[m])[0, 1]
        corr_raw = np.corrcoef(res[m], io[m])[0, 1]
        A = np.c_[np.ones(m.sum()), z[m]]
        cz, *_ = np.linalg.lstsq(A, res[m], rcond=None)
        h_share = 1 - np.var(res[m] - A @ cz) / np.var(res[m])
        rows.append(dict(pair=f"{g.ref_date:%Y-%m-%d}", res_cm=float(np.nanstd(res) * 100),
                         iono_cm=float(np.nanstd(io) * 100), corr=float(corr_raw),
                         corr_noplane=float(corr), height_share=float(h_share)))
        sfs.append(sf(res, rng) * 100)
    return rows, np.array(sfs)


def main():
    out = {}
    for d, name in [("ASCENDING", "asc 06:05"), ("DESCENDING", "desc 19:13")]:
        rows, sfs = run(d)
        out[name] = dict(rows=rows, sf=np.nanmedian(sfs, 0).tolist())
        c = np.array([r["corr"] for r in rows])
        cn = np.array([r["corr_noplane"] for r in rows])
        h = np.array([r["height_share"] for r in rows])
        print(f"{name}: {len(rows)} pairs")
        print(f"  1. residual vs its iono screen: median corr {np.median(c):+.2f} (|corr| median "
              f"{np.median(np.abs(c)):.2f}); without planes {np.median(cn):+.2f} (|.| {np.median(np.abs(cn)):.2f})")
        print(f"     pairs with |corr| > 0.3: {np.sum(np.abs(c) > 0.3)} of {len(c)}")
        print(f"  3. residual share explained by height: median {np.median(h) * 100:.0f}%")
    a, dd = np.array(out["asc 06:05"]["sf"]), np.array(out["desc 19:13"]["sf"])
    print("  2. noise vs distance (median, cm):")
    for L, x, y in zip(LAGS, a, dd):
        print(f"     {L:5.1f} km  asc {x:.2f}  desc {y:.2f}  ratio {y / x:.2f}")
    json.dump(out, open(ROOT / "data" / "processed" / "mexico-city" / "iono_vs_tropo.json", "w"), indent=1)


if __name__ == "__main__":
    main()
