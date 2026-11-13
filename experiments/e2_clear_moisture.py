"""E2: in clear sky, does GOES moisture (split window 10.3-12.3 um, 7.3 um water vapour) predict
the radar atmosphere?  (pre-registered in notes/experiments.md)

Zero shift (low clear-air vapour has little parallax). Null = shuffled maps among clear pairs.
Applied test: per pair, fit coefficients on the west half, apply to the east half (and back).

    uv run python scripts/e2_clear_moisture.py
"""
from __future__ import annotations

import json
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
from rasterio.warp import transform as warp_transform
from scipy.interpolate import LinearNDInterpolator

sys.path[:0] = [str(Path(__file__).parent), str(Path(__file__).parents[1] / "scripts")]
from e1_clouds_pattern import BBOX, blocks, plane_free, residuals  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
GOES = ROOT / "data" / "raw" / "goes19"
OUT = ROOT / "data" / "processed" / "mexico-city"


def band_map(band, t, epsg, BX, BY):
    fs = sorted(GOES.glob(f"bt{band:02d}_*.npz"))
    ts = [datetime.strptime(f.stem[5:], "%Y%m%dT%H%M") for f in fs]
    k = int(np.argmin([abs((s - t).total_seconds()) for s in ts]))
    if abs((ts[k] - t).total_seconds()) > 3600:
        return None, None
    g = np.load(fs[k])
    bt, lat, lon = g["bt"], g["lat"], g["lon"]
    m = np.isfinite(bt) & np.isfinite(lat)
    gx, gy = warp_transform("EPSG:4326", f"EPSG:{epsg}", lon[m].tolist(), lat[m].tolist())
    f = LinearNDInterpolator(np.c_[gx, gy], bt[m])
    inbox = (lat >= BBOX[1]) & (lat <= BBOX[3]) & (lon >= BBOX[0]) & (lon <= BBOX[2])
    return f(np.c_[BX.ravel(), BY.ravel()]).reshape(BX.shape), float(np.mean(bt[inbox] < 235))


def r_of(a, b, mask):
    x, y = a[mask], b[mask]
    return float(np.corrcoef(x, y)[0, 1])


def run(direction, series, rng):
    res, x, y, epsg = residuals(direction, series)
    _, bx, by = blocks(res[0][2], x, y)
    BX, BY = np.meshgrid(bx, by)
    full = np.ones(BX.shape, bool)
    maps = {}
    rows = []
    for d0, d1, r in res:
        for d in (d0, d1):
            if d not in maps:
                maps[d] = {b: band_map(b, d, epsg, BX, BY) for b in (10, 13, 15)}
        if any(maps[d][b][0] is None for d in (d0, d1) for b in (10, 13, 15)):
            continue
        cold = max(maps[d0][13][1], maps[d1][13][1])
        rb, _, _ = blocks(r, x, y)
        mask = np.isfinite(rb)
        if mask.sum() < 40:
            continue
        dsplit = (maps[d1][13][0] - maps[d1][15][0]) - (maps[d0][13][0] - maps[d0][15][0])
        rows.append(dict(d0=d0, cold=cold, rp=plane_free(rb, mask), mask=mask,
                         split=plane_free(dsplit, full),
                         b10=plane_free(maps[d1][10][0] - maps[d0][10][0], full),
                         b13=plane_free(maps[d1][13][0] - maps[d0][13][0], full)))
    clear = [r for r in rows if r["cold"] <= 0.05]
    n = len(clear)
    out = dict(n_clear=n)
    print(f"\n{direction}: {n} clear pairs (of {len(rows)})")
    for key, sign in [("split", -1), ("b10", +1), ("b13", +1)]:
        R = np.array([[r_of(ri["rp"], rj[key], ri["mask"]) for rj in clear] for ri in clear])
        own = np.diag(R)
        med = float(np.median(own))
        null = []
        for _ in range(2000):
            p = rng.permutation(n)
            while np.any(p == np.arange(n)):
                p = rng.permutation(n)
            null.append(abs(np.median(R[np.arange(n), p])))
        null = np.array(null)
        out[key] = dict(median_r=med, expected_sign=sign, null95=float(np.percentile(null, 95)),
                        p=float(np.mean(null >= abs(med))), frac_expected_sign=float(np.mean(np.sign(own) == sign)))
        print(f"  {key:5s}: median r {med:+.2f} (expected sign {sign:+d}) | null 95th {out[key]['null95']:.2f} | "
              f"p {out[key]['p']:.3f} | pairs with expected sign {np.mean(np.sign(own) == sign) * 100:.0f}%")
    # applied: west half -> east half and back, three predictors
    red = []
    for r in clear:
        H, W = r["rp"].shape
        cols = np.arange(W)[None, :].repeat(H, 0)
        Xall = np.stack([r["split"], r["b10"], r["b13"]], -1)
        vals = []
        for train, test in [(cols < W // 2, cols >= W // 2), (cols >= W // 2, cols < W // 2)]:
            tr, te = train & r["mask"], test & r["mask"]
            if tr.sum() < 20 or te.sum() < 20:
                continue
            A = np.c_[np.ones(tr.sum()), Xall[tr]]
            c, *_ = np.linalg.lstsq(A, r["rp"][tr], rcond=None)
            pred = np.c_[np.ones(te.sum()), Xall[te]] @ c
            vals.append(1 - np.var(r["rp"][te] - pred) / np.var(r["rp"][te]))
        if vals:
            red.append(np.mean(vals))
    out["holdout_var_reduction_median"] = float(np.median(red))
    out["holdout_var_reduction_all"] = [float(v) for v in red]
    print(f"  hold-out (west<->east) variance reduction: median {np.median(red) * 100:+.0f}% "
          f"(pairs improved {np.mean(np.array(red) > 0) * 100:.0f}%)")
    # same correlations in cloudy pairs, for curiosity (not pre-registered)
    cl = [r for r in rows if r["cold"] > 0.05]
    if cl:
        print(f"  [post-hoc] cloudy pairs ({len(cl)}): median r split {np.median([r_of(r['rp'], r['split'], r['mask']) for r in cl]):+.2f},"
              f" b10 {np.median([r_of(r['rp'], r['b10'], r['mask']) for r in cl]):+.2f}")
    return out


def main():
    rng = np.random.default_rng(11)
    out = {d: run(d, s, rng) for d, s in [("DESCENDING", "mexico-city"), ("ASCENDING", "mexico-city-asc")]}
    json.dump(out, open(OUT / "e2_clear_moisture.json", "w"), indent=1)


if __name__ == "__main__":
    main()
