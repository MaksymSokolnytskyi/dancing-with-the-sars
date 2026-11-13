"""Where is Mexico City sinking fastest, and which number is honest? (diagnostic, not a test)

    uv run python scripts/diag_depth.py sites/mexico-city-asc.json

Compares three ways of naming "the deepest": raw stack rate, smoothed + deramped cumulative (what the strokes draw),
and its 1st percentile (what the manifest called "deepest"). For the top spots prints sigma, chain pieces and the series.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from rasterio.warp import transform as warp_transform
from scipy import ndimage

from hachure import engine

ROOT = Path(__file__).resolve().parents[1]


def main():
    cfg = engine.load_site(sys.argv[1])
    v = engine.load_velocity(ROOT / "data/processed" / cfg["id"] / "engine_velocity.npz")
    years = (v.dates[-1] - v.dates[0]).days / 365.25
    ok = np.isfinite(v.rate) & (v.coherence >= cfg["velocity"].get("coh_min_draw", 0.45))
    r = np.where(ok, v.rate * 100, np.nan)
    print(f"{v.n_pairs} pairs, {years:.2f} yr, valid {ok.mean():.2f}")
    print(f"raw stack rate cm/yr: p0.1 {np.nanpercentile(r, 0.1):.1f}  p1 {np.nanpercentile(r, 1):.1f}  "
          f"p5 {np.nanpercentile(r, 5):.1f}  median {np.nanmedian(r):.1f}")
    # 500 m median blocks of the raw rate: robust 'neighbourhood' rate
    px = abs(v.x[1] - v.x[0])
    B = max(1, int(round(500 / px)))
    H, W = (r.shape[0] // B) * B, (r.shape[1] // B) * B
    with np.errstate(all="ignore"):
        blk = r[:H, :W].reshape(H // B, B, W // B, B)
        med = np.nanmedian(blk, axis=(1, 3))
        nval = np.isfinite(blk).sum((1, 3))
        sg = np.nanmedian((v.sigma * 100)[:H, :W].reshape(H // B, B, W // B, B), axis=(1, 3))
        seg = np.nanmedian(np.where(v.segments > 0, v.segments, np.nan)[:H, :W].reshape(H // B, B, W // B, B), axis=(1, 3))
        cum = np.nanmedian(((v.cum - v.cum[:1]) * 100)[:, :H, :W].reshape(len(v.dates), H // B, B, W // B, B), axis=(2, 4))
    med[nval < 0.5 * B * B] = np.nan
    print(f"500 m blocks (>=50% valid): p1 {np.nanpercentile(med, 1):.1f}, min {np.nanmin(med):.1f} cm/yr")
    order = np.argsort(np.nan_to_num(med, nan=1e9).ravel())
    xs, ys = v.x[:W].reshape(-1, B).mean(1), v.y[:H].reshape(-1, B).mean(1)
    seen = []
    print("deepest 500 m blocks (distinct spots >= 2 km apart):")
    for k in order:
        i, j = np.unravel_index(k, med.shape)
        if not np.isfinite(med[i, j]):
            break
        if any(np.hypot(xs[j] - a, ys[i] - b) < 2000 for a, b in seen):
            continue
        seen.append((xs[j], ys[i]))
        lon, lat = warp_transform(f"EPSG:{v.epsg}", "EPSG:4326", [xs[j]], [ys[i]])
        c = cum[:, i, j]
        jumps = np.abs(np.diff(c[np.isfinite(c)]))
        print(f"  {lon[0]:.4f}, {lat[0]:.4f}: rate {med[i, j]:6.1f} ± {sg[i, j]:.1f} cm/yr, pieces {seg[i, j]:.0f}, "
              f"cum end {c[-1]:6.1f} cm, biggest step {jumps.max() if len(jumps) else np.nan:.1f} cm, "
              f"series {np.round(c[::4], 0).tolist()}")
        if len(seen) >= 8:
            break


if __name__ == "__main__":
    main()
