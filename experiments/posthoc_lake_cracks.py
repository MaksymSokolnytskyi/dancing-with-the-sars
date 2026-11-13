"""POST-HOC follow-ups to test_lake_cracks.py (not pre-registered; labelled so in notes/experiments.md).

    uv run python scripts/posthoc_lake_cracks.py sites/mexico-city-asc.json

a) L1 detail: in which zone are the fast pixels that are not in Lago?
b) C1 symmetric: observed crack points restricted to Lago/Transición, like the null.
c) C1 at a finer scale: gradient of the raw rate smoothed only 150 m (cracks are ~100 m features).
"""
from __future__ import annotations

import json
import sys

import numpy as np

sys.path[:0] = [str(__import__("pathlib").Path(__file__).parent), str(__import__("pathlib").Path(__file__).parents[1] / "scripts")]
import test_lake_cracks as t  # noqa: E402
from hachure import engine  # noqa: E402


def shift_test(rc, G, valid, allowed, res, obs_keep, n=500, seed=0):
    Gp = t.sample(G, rc)
    ok = np.isfinite(Gp) & (t.sample(valid.astype(float), rc) == 1) & obs_keep
    obs = float(np.median(Gp[ok]))
    rng, null = np.random.default_rng(seed), []
    while len(null) < n:
        dx, dy = rng.uniform(-3000, 3000, 2)
        if np.hypot(dx, dy) < 1000:
            continue
        sh = rc + np.array([-dy / res, dx / res])
        g = t.sample(G, sh)
        k = (t.sample(allowed.astype(float), sh) == 1) & np.isfinite(g) & (t.sample(valid.astype(float), sh) == 1)
        null.append(np.median(g[k]) if k.sum() >= 50 else np.nan)
    null = np.array(null)
    return dict(points=int(ok.sum()), obs=round(obs, 3), null_median=round(float(np.nanmedian(null)), 3),
                null_p95=round(float(np.nanpercentile(null, 95)), 3), share_null_ge_obs=round(float(np.nanmean(null >= obs)), 4))


def main():
    cfg = engine.load_site(sys.argv[1])
    F = t.fields(cfg)
    zones = t.burn(t.load("SPC_Zonificacion_Geotecnica.geojson"), "EPSG:6365", F,
                   lambda pr: t.ZONE.get(str(pr.get("INTENSIDAD")).strip(), 0))
    out = {}
    r = F["r_raw"]
    fast = np.isfinite(r) & (r <= -10)
    out["a_fast_px_by_zone"] = {t.NAMES.get(k, "none"): round(float((zones[fast] == k).mean()), 3) for k in (0, 1, 2, 3)}
    out["a_fast_px_by_zone_le_-15"] = {t.NAMES.get(k, "none"): round(float((zones[np.isfinite(r) & (r <= -15)] == k).mean()), 3)
                                       for k in (0, 1, 2, 3)}
    rc, _ = t.crack_points(F)
    allowed = np.isin(zones, (2, 3))
    zp = t.sample(zones, rc)
    out["b_symmetric_smoothed400"] = shift_test(rc, F["G"], F["valid_s"], allowed, F["res"], np.isin(zp, (2, 3)))
    out["b_lago_only_smoothed400"] = shift_test(rc, F["G"], F["valid_s"], zones == 3, F["res"], zp == 3)
    # c) finer gradient: raw rate, 150 m smoothing
    valid = np.isfinite(r)
    s, _ = engine._smooth_nan(np.nan_to_num(r), 150 / F["res"], valid)
    s[~valid] = np.nan
    gy, gx = np.gradient(s)
    G150 = np.hypot(gx, gy) / (F["res"] / 1000)
    out["c_symmetric_raw150"] = shift_test(rc, G150, valid, allowed, F["res"], np.isin(zp, (2, 3)))
    out["c_lago_only_raw150"] = shift_test(rc, G150, valid, zones == 3, F["res"], zp == 3)
    print(json.dumps(out, indent=1, ensure_ascii=False))


if __name__ == "__main__":
    main()
