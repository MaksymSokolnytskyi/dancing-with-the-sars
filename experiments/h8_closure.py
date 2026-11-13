"""H8: triangle closure of independently processed GUNW pairs (12 + 12 vs 24 days). Pre-registered.

    uv run python scripts/h8_closure.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.ndimage import gaussian_filter

from hachure import nisar

ROOT = Path(__file__).resolve().parents[1]
BBOX = (-99.32, 19.17, -98.82, 19.58)
REF_BBOX = (-99.14, 19.53, -99.09, 19.57)
CACHE = ROOT / "data" / "raw" / "nisar_cache" / ("bbox_" + "_".join(f"{v:.3f}" for v in BBOX) + "_h2240")
OUT = ROOT / "data" / "processed" / "mexico-city"


def raw(g):
    p = CACHE / f"{g.name}.npz"
    if not p.exists():
        print(f"  downloading window of {g.name[:60]}...", flush=True)
        with nisar.open_h5(g.url) as h5:
            nisar.read_raw(h5, BBOX, g.days, height_m=2240).save(p)
    return nisar.RawWindow.load(p)


def smooth_std(a, km=2.0, dx=80.0):
    ok = np.isfinite(a)
    s = km * 1000 / dx / 2.355
    num, den = gaussian_filter(np.where(ok, a, 0), s), gaussian_filter(ok.astype(float), s)
    with np.errstate(all="ignore"):
        sm = np.where(ok & (den > 0.3), num / den, np.nan)
    return float(np.nanstd(sm))


def main():
    out = {}
    for d in ("DESCENDING", "ASCENDING"):
        found = nisar.dominant_track(nisar.search((BBOX[0] + BBOX[2]) / 2, (BBOX[1] + BBOX[3]) / 2, direction=d))
        twelve = {g.ref_date.date(): g for g in found if round(g.days) == 12}
        for g24 in [g for g in found if round(g.days) == 24]:
            a = twelve.get(g24.ref_date.date())
            if a is None or a.sec_date.date() not in twelve:
                continue
            b = twelve[a.sec_date.date()]
            if b.sec_date.date() != g24.sec_date.date():
                continue
            print(f"{d}: triangle {g24.ref_date:%Y-%m-%d} -> {a.sec_date:%Y-%m-%d} -> {g24.sec_date:%Y-%m-%d}")
            raws = [raw(x) for x in (a, b, g24)]
            res = {}
            for iono in (True, False):
                pw = [nisar.process(r, ref_bbox=REF_BBOX, coh_min=0.3, iono=iono) for r in raws]
                disps, _, x, y = nisar.align(pw)
                c = (disps[2] - disps[0] - disps[1]) * 100
                c -= np.nanmedian(c)
                key = "iono" if iono else "no_iono"
                res[key] = dict(std=float(np.nanstd(c)), std_2km=smooth_std(c, dx=abs(x[1] - x[0])),
                                p99=float(np.nanpercentile(np.abs(c), 99)), valid=float(np.isfinite(c).mean()))
                if iono:
                    io = [r.iono for r in raws]
                    res["c_map"] = c
            ci, cn = res["iono"], res["no_iono"]
            print(f"  closure with iono    std {ci['std']:.2f} cm, after 2 km smoothing {ci['std_2km']:.2f} "
                  f"(ratio {ci['std_2km'] / ci['std']:.2f}), |c| 99th pct {ci['p99']:.1f}")
            print(f"  closure without iono std {cn['std']:.2f} cm, after 2 km smoothing {cn['std_2km']:.2f} "
                  f"(ratio {cn['std_2km'] / cn['std']:.2f})")
            print(f"  P1 {ci['std']:.2f} >= 0.5: {ci['std'] >= 0.5} | P2 ratio >= 0.7: {ci['std_2km'] / ci['std'] >= 0.7}"
                  f" | P3 no-iono/iono = {cn['std'] / ci['std']:.2f} < 0.5: {cn['std'] / ci['std'] < 0.5}")
            np.save(OUT / f"h8_closure_{d[:4].lower()}.npy", res.pop("c_map"))
            out[d] = res
    json.dump(out, open(OUT / "h8_closure.json", "w"), indent=1)


if __name__ == "__main__":
    main()
