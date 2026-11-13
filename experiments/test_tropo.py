"""Criteria T1-T3 for the height-correlated delay fix (pre-registered in notes/experiments.md, section "T").

    uv run python scripts/test_tropo.py sites/mexico-city-asc.json
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from rasterio.warp import transform as warp_transform

from hachure import engine

ROOT = Path(__file__).resolve().parents[1]
SW = (-99.31, 19.20, -99.26, 19.30)


def sw_median(v):
    lon, lat = SW[0::2], SW[1::2]
    xs, ys = warp_transform("EPSG:4326", f"EPSG:{v.epsg}", list(lon), list(lat))
    m = ((v.x[None, :] >= min(xs)) & (v.x[None, :] <= max(xs)) & (v.y[:, None] >= min(ys)) & (v.y[:, None] <= max(ys)))
    ok = m & np.isfinite(v.rate) & (v.coherence >= 0.45)
    return float(np.median(v.rate[ok]) * 100), int(ok.sum())


def summary(cfg, v):
    sw, n = sw_median(v)
    g = engine.gnss_check(cfg, v)
    h = engine.headline(cfg, v)
    return dict(sw_median_cm_yr=round(sw, 2), sw_px=n,
                mmx1={x["station"]: x.get("nisar_cm_yr") for x in g},
                headline=dict(rate=h["rate_cm_yr"], sigma=h["sigma_cm_yr"], lon=h["lon"], lat=h["lat"]))


def main():
    cfg = engine.load_site(sys.argv[1])
    v0 = engine.load_velocity(ROOT / "data/processed" / cfg["id"] / "engine_velocity.npz")
    v1 = engine.tropo_correct(cfg, v0)
    b, a = summary(cfg, v0), summary(cfg, v1)
    gn = {x["station"]: x["rate_cm_yr"] for x in cfg.get("gnss", [])}
    st = next(iter(gn))
    res = {"before": b, "after": a, "tropo": getattr(v1, "tropo", None),
           "T1": {"value": a["sw_median_cm_yr"], "passed": abs(a["sw_median_cm_yr"]) <= 4},
           "T2": {"value": a["mmx1"][st], "gnss": gn[st], "passed": abs(a["mmx1"][st] - gn[st]) <= 3},
           "T3": {"shift": round(a["headline"]["rate"] - b["headline"]["rate"], 2),
                  "passed": abs(a["headline"]["rate"] - b["headline"]["rate"]) < 3}}
    print(json.dumps(res, indent=1, ensure_ascii=False))
    (ROOT / "data/processed" / cfg["id"] / "test_tropo.json").write_text(json.dumps(res, indent=1), encoding="utf-8")


if __name__ == "__main__":
    main()
