"""T1: turbulent-layer depth from the break of the residual structure function (pre-registered).

    uv run python scripts/t1_layer_depth.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

from hachure import nisar

ROOT = Path(__file__).resolve().parents[1]
BBOX = (-99.32, 19.17, -98.82, 19.58)
REF_BBOX = (-99.14, 19.53, -99.09, 19.57)
CACHE = ROOT / "data" / "raw" / "nisar_cache" / ("bbox_" + "_".join(f"{v:.3f}" for v in BBOX) + "_h2240")
LAGS = np.array([0.16, 0.24, 0.32, 0.48, 0.64, 0.96, 1.28, 2.0, 3.0, 5.0, 8.0, 12.0, 20.0])


def sf(r, rng, n=200_000, dx=0.08):
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


def model(th, r):
    nug, A, rb = np.exp(th)
    x = r / rb
    # smooth broken power law for the std: x^(5/6) below rb, x^(1/3) above
    turb = A * x ** (5 / 6) / (1 + x) ** (5 / 6 - 1 / 3)
    return np.sqrt(nug ** 2 + turb ** 2)


def run(direction):
    found = nisar.search((BBOX[0] + BBOX[2]) / 2, (BBOX[1] + BBOX[3]) / 2, direction=direction)
    chain = nisar.sequential_chain(nisar.dominant_track(found))
    pairs = [nisar.process(nisar.RawWindow.load(CACHE / f"{g.name}.npz"), ref_bbox=REF_BBOX, coh_min=0.3)
             for g in chain]
    disps, _, x, y = nisar.align(pairs)
    days = np.array([p.days for p in pairs])
    ok = np.isfinite(disps)
    tot = np.where(ok, disps, 0).sum(0)
    span = np.where(ok, days[:, None, None], 0).sum(0)
    rng = np.random.default_rng(2)
    sfs = []
    for i in range(len(pairs)):
        with np.errstate(all="ignore"):
            rl = (tot - np.where(ok[i], disps[i], 0)) / (span - days[i])
        r = np.where(ok[i] & (span - days[i] > 200), disps[i] - rl * days[i], np.nan) * 100
        sfs.append(sf(r - np.nanmedian(r), rng))
    sfs = np.array(sfs)
    med = np.nanmedian(sfs, 0)
    m = np.isfinite(med)
    fit = least_squares(lambda th: np.log(model(th, LAGS[m])) - np.log(med[m]), x0=np.log([0.3, 1.0, 1.0]),
                        bounds=(np.log([1e-3, 1e-3, 0.05]), np.log([5, 20, 50])))
    nug, A, rb = np.exp(fit.x)
    # bootstrap over pairs for the break scale
    rbs = []
    for _ in range(300):
        s = sfs[rng.integers(0, len(sfs), len(sfs))]
        mm = np.nanmedian(s, 0)
        f = least_squares(lambda th: np.log(model(th, LAGS[m])) - np.log(mm[m]), x0=fit.x,
                          bounds=(np.log([1e-3, 1e-3, 0.05]), np.log([5, 20, 50])))
        rbs.append(np.exp(f.x[2]))
    sl = lambda a, b: np.polyfit(np.log(LAGS[(LAGS >= a) & (LAGS <= b)]),            # noqa: E731
                                 np.log(med[(LAGS >= a) & (LAGS <= b)]), 1)[0]
    s_small, s_large = sl(0.5, 2.0), sl(5.0, 20.0)
    print(f"{direction}: std vs lag (cm): " + " ".join(f"{L:g}:{v:.2f}" for L, v in zip(LAGS, med)))
    print(f"  slopes: 0.5-2 km {s_small:.2f}, 5-20 km {s_large:.2f} | fit: nugget {nug:.2f} cm, A {A:.2f} cm, "
          f"break r_b {rb:.2f} km (bootstrap 16-84%: {np.percentile(rbs, 16):.2f}-{np.percentile(rbs, 84):.2f})")
    return dict(sf=med.tolist(), slope_small=s_small, slope_large=s_large, nugget=nug, A=A, rb=rb,
                rb_ci=[float(np.percentile(rbs, 16)), float(np.percentile(rbs, 84))],
                rb_boot=[float(v) for v in rbs])


def main():
    out = {d: run(d) for d in ("ASCENDING", "DESCENDING")}
    a, d = out["ASCENDING"], out["DESCENDING"]
    ratio = d["rb"] / a["rb"]
    pboot = np.mean(np.array(d["rb_boot"])[:300] > 1.5 * np.array(a["rb_boot"])[:300])
    print(f"P1: asc {a['slope_small']:.2f} > {a['slope_large']:.2f}: {a['slope_small'] > a['slope_large']}, "
          f"desc {d['slope_small']:.2f} > {d['slope_large']:.2f}: {d['slope_small'] > d['slope_large']}")
    print(f"P2: evening/morning break ratio {ratio:.2f} (need >= 1.5); bootstrap share > 1.5: {pboot:.2f}")
    print(f"P3: evening r_b {d['rb']:.2f} km in 1-4: {1 <= d['rb'] <= 4}")
    json.dump(out, open(ROOT / "data/processed/mexico-city/t1_layer_depth.json", "w"), indent=1)


if __name__ == "__main__":
    main()
