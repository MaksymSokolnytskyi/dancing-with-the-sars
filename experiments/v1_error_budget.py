"""V1: does the node+edge (tridiagonal) noise model predict the observed asc-desc rate difference?

    uv run python scripts/v1_error_budget.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path[:0] = [str(Path(__file__).parent), str(Path(__file__).parents[1] / "scripts")]
from h6_horizontal import smooth  # noqa: E402

from hachure import nisar

ROOT = Path(__file__).resolve().parents[1]
BBOX = (-99.32, 19.17, -98.82, 19.58)
REF_BBOX = (-99.14, 19.53, -99.09, 19.57)
CACHE = ROOT / "data" / "raw" / "nisar_cache" / ("bbox_" + "_".join(f"{v:.3f}" for v in BBOX) + "_h2240")


def predicted(direction, scales=(1, 10)):
    found = nisar.search((BBOX[0] + BBOX[2]) / 2, (BBOX[1] + BBOX[3]) / 2, direction=direction)
    chain = nisar.sequential_chain(nisar.dominant_track(found))
    pairs = [nisar.process(nisar.RawWindow.load(CACHE / f"{g.name}.npz"), ref_bbox=REF_BBOX, coh_min=0.3)
             for g in chain]
    disps, _, x, y = nisar.align(pairs)
    days = np.array([p.days for p in pairs])
    T = days.sum() / 365.25
    ok = np.isfinite(disps)
    tot = np.where(ok, disps, 0).sum(0)
    span = np.where(ok, days[:, None, None], 0).sum(0)
    dx = abs(x[1] - x[0])
    out = {}
    for L in scales:
        rs = []
        for i in range(len(pairs)):
            with np.errstate(all="ignore"):
                rl = (tot - np.where(ok[i], disps[i], 0)) / (span - days[i])
            r = np.where(ok[i], disps[i] - rl * days[i], np.nan) * 100
            r = smooth(r, L * 1000 / dx / 2.355)
            rs.append(r - np.nanmedian(r))
        n = len(rs)
        C = np.full((n, n), np.nan)
        for j in range(n):
            for k in range(j, min(n, j + 6)):
                m = np.isfinite(rs[j]) & np.isfinite(rs[k])
                C[j, k] = C[k, j] = np.mean(rs[j][m] * rs[k][m]) - rs[j][m].mean() * rs[k][m].mean()
        diag = np.trace(C)
        lag1 = np.nansum(np.diag(C, 1))
        far = sum(np.nansum(np.diag(C, l)) for l in range(2, 6))
        out[L] = dict(var_pred=float((diag + 2 * lag1) / T ** 2), diag=float(diag), lag1=float(lag1),
                      lag2to5=float(far), T=float(T))
        print(f"  {direction[:4]} L={L:2d} km: sum var {diag:7.2f}, 2*sum lag1 {2 * lag1:+7.2f}, "
              f"2*sum lags 2-5 {2 * far:+7.2f} cm^2 -> predicted rate std {np.sqrt(out[L]['var_pred']):.2f} cm/yr")
    return out


def main():
    obs = json.load(open(ROOT / "data/processed/mexico-city/h8b_no_iono.json"))["with_iono"]
    pa, pd = predicted("ASCENDING"), predicted("DESCENDING")
    res = {}
    for L in (1, 10):
        o = obs[f"D_std_{L}km"] ** 2
        p = pa[L]["var_pred"] + pd[L]["var_pred"]
        res[L] = dict(observed_std=float(np.sqrt(o)), predicted_std=float(np.sqrt(p)), R=float(o / p),
                      asc=pa[L], desc=pd[L])
        print(f"L={L:2d} km: observed std(D) {np.sqrt(o):.2f}, predicted {np.sqrt(p):.2f} cm/yr -> R = {o / p:.2f}")
    json.dump(res, open(ROOT / "data/processed/mexico-city/v1_error_budget.json", "w"), indent=1)


if __name__ == "__main__":
    main()
