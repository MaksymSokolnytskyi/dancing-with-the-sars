"""D, unblinding step: compare the frozen SBI / least-squares results with the USGS solutions.

Run only after data/processed/events/sbi_<event>.json is final (results recorded in notes first).

    uv run python scripts/sbi_compare.py
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
IDS = {"nagqu": "us6000u43g", "kumamoto": "us6000tgb9", "panama": "us6000u18k", "panama-wide": "us6000u18k"}


def usgs(eid):
    d = json.load(open(ROOT / "data" / "raw" / "usgs" / f"{eid}.json"))
    mt = d["properties"]["products"]["moment-tensor"][0]["properties"]
    planes = [tuple(float(mt[f"nodal-plane-{i}-{k}"]) for k in ("strike", "dip", "rake")) for i in (1, 2)]
    return dict(planes=planes, mw=float(mt.get("derived-magnitude", d["properties"]["mag"])),
                lon=d["geometry"]["coordinates"][0], lat=d["geometry"]["coordinates"][1],
                centroid_depth=float(mt.get("derived-depth", "nan")))


def angdiff(a, b):
    return abs((a - b + 180) % 360 - 180)


def plane_err(st, dp, plane):
    s2, d2, _ = plane
    ds = angdiff(st, s2)
    if dp > 75 and d2 > 75:                      # near-vertical: strike and strike+180 are the same plane
        ds = min(ds, angdiff(st + 180, s2))
    return ds, abs(dp - d2)


def surface_dist_to_epicentre(p):
    """Distance (km) from the epicentre (origin of our E/N) to the fault's surface projection."""
    st, dp = math.radians(p["strike"]), math.radians(p["dip"])
    s = np.array([math.sin(st), math.cos(st)])
    dd = np.array([math.cos(st), -math.sin(st)]) * math.cos(dp)
    c = np.array([p["e0"], p["n0"]])
    # coordinates of the origin in the fault frame
    rel = -c
    ua = max(abs(rel @ s) - p["L"] / 2, 0)
    hw = p["W"] / 2 * math.cos(dp)                 # half-width of the surface projection
    perp = abs(rel @ np.array([math.cos(st), -math.sin(st)]))
    ub = max(perp - hw, 0)
    return math.hypot(ua, ub)


def judge(name, p, ref):
    errs = [plane_err(p["strike"], p["dip"], pl) for pl in ref["planes"]]
    i = int(np.argmin([e[0] + e[1] for e in errs]))
    ds, dd = errs[i]
    dmw = abs(p["Mw"] - ref["mw"])
    dist = surface_dist_to_epicentre(p)
    ok = dict(strike=ds <= 15, dip=dd <= 15, mw=dmw <= 0.2, location=dist <= 10)
    print(f"  {name:4s}: strike {p['strike']:.0f} vs {ref['planes'][i][0]:.0f} (diff {ds:.0f}) | dip {p['dip']:.0f} vs "
          f"{ref['planes'][i][1]:.0f} (diff {dd:.0f}) | rake {p['rake']:.0f} vs {ref['planes'][i][2]:.0f} | Mw "
          f"{p['Mw']:.2f} vs {ref['mw']:.2f} | fault-to-epicentre {dist:.1f} km  ->  "
          + ", ".join(f"{k} {'PASS' if v else 'fail'}" for k, v in ok.items()))
    return dict(plane=i + 1, d_strike=ds, d_dip=dd, d_mw=dmw, dist_km=dist, passes=ok)


def main():
    out = {}
    for ev, eid in IDS.items():
        f = ROOT / "data" / "processed" / "events" / f"sbi_{ev}.json"
        if not f.exists():
            continue
        r = json.load(open(f))
        ref = usgs(eid)
        print(f"{ev}: USGS planes {ref['planes']}, Mw {ref['mw']}, centroid depth {ref['centroid_depth']} km")
        sbi = {k: v[0] for k, v in r["sbi"].items()}
        out[ev] = dict(usgs=ref, sbi=judge("SBI", sbi, ref), lsq=judge("LSQ", r["lsq"], ref),
                       inference_ms=r["inference_ms"])
        for arm, h in r.get("hybrid", {}).items():
            hp = h["params"]
            out[ev][f"hybrid_{arm}"] = judge(arm, hp, ref)
        print(f"  SBI inference {r['inference_ms']:.0f} ms, training {r['train_s'] / 60:.1f} min, LSQ {r['lsq_s']:.0f} s")
    json.dump(out, open(ROOT / "data" / "processed" / "events" / "sbi_vs_usgs.json", "w"), indent=1)


if __name__ == "__main__":
    main()
