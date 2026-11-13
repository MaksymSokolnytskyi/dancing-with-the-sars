"""G1': per-date atmospheric noise from the radar alone, via the covariance of adjacent pairs.

Pair k = phi_{k+1} - phi_k (date potentials on a path graph). Edge covariance C = B Sigma B^T:
off-diagonal of adjacent edges = -sigma^2 of the shared date. Pre-registered in notes/experiments.md.

    uv run python scripts/g1_date_noise.py
"""
from __future__ import annotations

import json
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

sys.path[:0] = [str(Path(__file__).parent), str(Path(__file__).parents[1] / "scripts")]
from b6_goes_weights import goes_g  # noqa: E402
from e1_clouds_pattern import BBOX, REF_BBOX, plane_free  # noqa: E402

from hachure import nisar

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "processed" / "mexico-city"


def run(direction, series):
    found = nisar.search((BBOX[0] + BBOX[2]) / 2, (BBOX[1] + BBOX[3]) / 2, direction=direction)
    chain = nisar.sequential_chain(nisar.dominant_track(found))
    cache = ROOT / "data" / "raw" / "nisar_cache" / ("bbox_" + "_".join(f"{v:.3f}" for v in BBOX) + "_h2240")
    pairs = [nisar.process(nisar.RawWindow.load(cache / f"{g.name}.npz"), ref_bbox=REF_BBOX, coh_min=0.3)
             for g in chain]
    disps, _, x, y = nisar.align(pairs)
    days = np.array([p.days for p in pairs])
    z = np.load(ROOT / "data" / "processed" / series / "nisar_series.npz")
    dates = [datetime.fromisoformat(str(d)) for d in z["dates"]]
    ok = np.isfinite(disps)
    tot = np.where(ok, disps, 0).sum(0)
    span = np.where(ok, days[:, None, None], 0).sum(0)
    n = len(pairs)
    out = {"median": [], "plane": []}
    var_k = {"median": [], "plane": []}
    for k in range(n - 1):
        both = ok[k] & ok[k + 1]
        with np.errstate(all="ignore"):
            v = (tot - np.where(ok[k], disps[k], 0) - np.where(ok[k + 1], disps[k + 1], 0)) / \
                (span - np.where(ok[k], days[k], 0) - np.where(ok[k + 1], days[k + 1], 0))
        m = both & np.isfinite(v)
        r1 = np.where(m, disps[k] - v * days[k], np.nan)
        r2 = np.where(m, disps[k + 1] - v * days[k + 1], np.nan)
        for mode in ("median", "plane"):
            if mode == "plane":
                a, b = plane_free(r1, m), plane_free(r2, m)
            else:
                a, b = r1 - np.nanmedian(r1), r2 - np.nanmedian(r2)
            av, bv = a[m] * 100, b[m] * 100
            out[mode].append(-np.mean((av - av.mean()) * (bv - bv.mean())))
            var_k[mode].append((np.var(av), np.var(bv)))
    s2 = np.array(out["median"])                     # interior dates 1..n-1
    s2p = np.array(out["plane"])
    inner = dates[1:n]
    g = np.array([goes_g(d) for d in inner])
    print(f"\n{direction}: {n} pairs, {len(inner)} interior dates")
    for d, a, b, gg in zip(inner, s2, s2p, g):
        print(f"  {d:%Y-%m-%d}  sigma^2 {a:+6.2f} cm^2 (plane-free {b:+6.2f})  GOES BT std {gg:5.1f} K")
    pos = float(np.mean(s2 > 0)), float(np.mean(s2p > 0))
    sig = np.sqrt(np.clip(s2, 0, None))
    gm = np.isfinite(g)
    rho, p = spearmanr(sig[gm], g[gm])
    X = np.c_[np.ones(gm.sum()), g[gm]]
    (a, b), *_ = np.linalg.lstsq(X, s2[gm], rcond=None)
    rng = np.random.default_rng(5)
    boot = np.array([np.linalg.lstsq(X[i], s2[gm][i], rcond=None)[0]
                     for i in (rng.integers(0, gm.sum(), gm.sum()) for _ in range(4000))])
    ca, cb = np.percentile(boot[:, 0], [2.5, 97.5]), np.percentile(boot[:, 1], [2.5, 97.5])
    # consistency: var(r_k) vs sigma_k^2 + sigma_{k+1}^2 for pairs with both ends interior
    cons = [var_k["median"][k][1] / (s2[k] + s2[k + 1]) for k in range(len(s2) - 1) if s2[k] + s2[k + 1] > 0]
    print(f"  P1 sigma^2 > 0: {pos[0] * 100:.0f}% (plane-free {pos[1] * 100:.0f}%)")
    print(f"  P2 Spearman(sigma, GOES BT std) = {rho:+.2f} (p {p:.3f}, n {gm.sum()})")
    print(f"  P3 sigma^2 = a + b g: a {a:+.2f} [{ca[0]:+.2f}, {ca[1]:+.2f}] cm^2, b {b:+.3f} [{cb[0]:+.3f}, {cb[1]:+.3f}]")
    print(f"  consistency var(pair) / (sigma_k^2 + sigma_k+1^2): median {np.median(cons):.2f}")
    print(f"  median sigma {np.median(sig):.2f} cm")
    return dict(dates=[f"{d:%Y-%m-%d}" for d in inner], s2=s2.tolist(), s2_plane=s2p.tolist(), goes_std=g.tolist(),
                frac_pos=pos, spearman=[float(rho), float(p)], a=float(a), a_ci=ca.tolist(), b=float(b),
                b_ci=cb.tolist(), consistency=float(np.median(cons)), median_sigma=float(np.median(sig)))


def main():
    res = {d: run(d, s) for d, s in [("DESCENDING", "mexico-city"), ("ASCENDING", "mexico-city-asc")]}
    ratio = res["DESCENDING"]["median_sigma"] / res["ASCENDING"]["median_sigma"]
    res["P4_ratio"] = ratio
    print(f"\nP4 evening/morning median sigma ratio: {ratio:.2f}")
    json.dump(res, open(OUT / "g1_date_noise.json", "w"), indent=1)


if __name__ == "__main__":
    main()
