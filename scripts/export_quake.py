"""Quake card for the web map: two-look fault estimate (Q-2look v1) + confidence flag + map layers.

    SBI_BLOCK=24 SBI_LMAX=200e3 SBI_XYMAX=60e3 HYB_XYMAX=80e3 HYB_LMAX=300e3 \
        uv run --group ml python scripts/export_quake.py sites/panama-quake.json

Writes web/data/<id>/: card.json (kind "card"), fault.geojson (polygons, style "fault"),
los_asc.json / los_desc.json (cells: observed, model, residual LOS in cm), then the engine adds
relief + manifest from the site config. Confidence = share of the top-16 proposals that agree with the
best one within 15 deg strike and dip (D-conf, AUC 0.93 on synthetic events); threshold from the config.
The reference (USGS) is optional, shown only as a labelled comparison, never used by the estimate.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
import torch
from rasterio.warp import transform as warp_transform

sys.path.insert(0, str(Path(__file__).parent))
from sbi_conf import angdiff, q2look  # noqa: E402
from sbi_fault import BLOCK, EVENTS, OUT, Net, Sim, dev, load_event  # noqa: E402
from sbi_hybrid2 import set_heading  # noqa: E402

from hachure import engine, nisar

ROOT = Path(__file__).resolve().parents[1]


def fault_outline(q, epi, epsg):
    """Surface projection of the fault rectangle (+ top edge) in lon/lat. Okada convention: (e0, n0, depth) = centroid."""
    ex, ey = warp_transform("EPSG:4326", f"EPSG:{epsg}", [epi[0]], [epi[1]])
    st, dp = math.radians(q["strike"]), math.radians(q["dip"])
    s = np.array([math.sin(st), math.cos(st)])                     # along strike (E, N)
    d = np.array([math.cos(st), -math.sin(st)])                    # down-dip, horizontal
    c = np.array([ex[0] + q["e0"] * 1e3, ey[0] + q["n0"] * 1e3])
    hl, hw = q["L"] * 1e3 / 2, q["W"] * 1e3 / 2 * math.cos(dp)
    ring = [c - hl * s - hw * d, c + hl * s - hw * d, c + hl * s + hw * d, c - hl * s + hw * d]
    ring.append(ring[0])
    xs, ys = zip(*ring)
    lon, lat = warp_transform(f"EPSG:{epsg}", "EPSG:4326", list(xs), list(ys))
    pts = np.round(np.column_stack([lon, lat]), 5).tolist()
    return pts, pts[:2]                                             # rectangle, top edge (shallow side)


def los_cells(name, sim, model, epi):
    """Block centres in lon/lat with observed / model / residual LOS (cm), all de-planed like the fit."""
    ev = EVENTS[name]
    co = sorted((ROOT / "data/raw/nisar_events" / ev["dir"]).glob(ev["co"]))[0]
    r = nisar.RawWindow.load(co)
    ex, ey = warp_transform("EPSG:4326", f"EPSG:{r.epsg}", [epi[0]], [epi[1]])
    e, n = sim.e.cpu().numpy() + ex[0], sim.n.cpu().numpy() + ey[0]
    lon, lat = warp_transform(f"EPSG:{r.epsg}", "EPSG:4326", e.tolist(), n.tolist())
    obs = sim.deplane(sim.obs[None])[0].cpu().numpy() * 100
    mod = model * 100
    res = obs - mod
    cells = [[round(lon[i], 5), round(lat[i], 5), round(float(obs[i]), 1), round(float(mod[i]), 1), round(float(res[i]), 1)]
             for i in range(len(obs))]
    px = abs(r.x[1] - r.x[0]) * BLOCK
    return {"cell_km": px / 1000, "columns": ["lon", "lat", "los_cm", "model_cm", "residual_cm"], "value": "los_cm",
            "default_mode": "los_cm", "epsg": r.epsg, "granule": co.name,
            "modes": [{"key": "los_cm", "label": {"en": "radar: observed", "uk": "радар: виміряно"}, "ramp": "los"},
                      {"key": "model_cm", "label": {"en": "fault model", "uk": "модель розлому"}, "ramp": "los"},
                      {"key": "residual_cm", "label": {"en": "residual", "uk": "залишок"}, "ramp": "los"}],
            "rms_obs_cm": round(float(np.sqrt((obs ** 2).mean())), 2), "rms_res_cm": round(float(np.sqrt((res ** 2).mean())), 2),
            "cells": cells}, r.epsg


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("site")
    ap.add_argument("--out", default=str(ROOT / "web" / "data"))
    a = ap.parse_args()
    cfg = engine.load_site(a.site)
    qc = cfg["quake"]
    out = Path(a.out) / cfg["id"]
    out.mkdir(parents=True, exist_ok=True)
    sA, sD = Sim(load_event(qc["asc"])), Sim(load_event(qc["desc"]))
    set_heading(sD, EVENTS[qc["desc"]]["heading"])
    net = Net().to(dev)
    net.load_state_dict(torch.load(ROOT / qc["net"]))
    net.eval()
    sims = [sA, sD]
    # several seeds: one random draw must not decide the card. Shown = lowest two-look residual; spread reported.
    runs = []
    for seed in range(qc.get("seeds", 12)):
        torch.manual_seed(seed)
        runs.append(q2look(net, sims, [s.deplane(s.obs[None]) for s in sims], n_net=qc.get("n_net", 2000), top=qc.get("top", 16)))
    q, conf_run, rel = min(runs, key=lambda r: r[2])
    agree = [r for r in runs if angdiff(r[0]["strike"], q["strike"]) <= 15 and abs(r[0]["dip"] - q["dip"]) <= 15]
    conf = float(np.median([r[1] for r in runs]))
    spread = {"runs": len(runs), "runs_agreeing": len(agree),
              "strike": [round(min(r[0]["strike"] for r in runs)), round(max(r[0]["strike"] for r in runs))],
              "dip": [round(min(r[0]["dip"] for r in runs)), round(max(r[0]["dip"] for r in runs))],
              "Mw": [round(min(r[0]["Mw"] for r in runs), 2), round(max(r[0]["Mw"] for r in runs), 2)],
              "confidence": [round(min(r[1] for r in runs), 2), round(max(r[1] for r in runs), 2)]}
    for r in runs:
        print(f"    strike {r[0]['strike']:5.0f} dip {r[0]['dip']:4.0f} Mw {r[0]['Mw']:.2f} conf {r[1]:.2f} rel {r[2]:.2f}")
    p = {k: torch.tensor([q[k] * (1e3 if k in ("e0", "n0", "depth", "L", "W") else 1)], device=dev)
         for k in ("e0", "n0", "depth", "strike", "dip", "L", "W", "rake", "slip")}
    epi = qc["epicentre"]
    layers = []
    for tag, nm, s in (("asc", qc["asc"], sA), ("desc", qc["desc"], sD)):
        with torch.no_grad():
            model = s.deplane(s.los(p))[0].cpu().numpy()
        g, epsg = los_cells(nm, s, model, epi)
        (out / f"los_{tag}.json").write_text(json.dumps(g, separators=(",", ":"), ensure_ascii=False), encoding="utf-8")
        layers.append({"id": f"los_{tag}", "kind": "cells", "file": f"los_{tag}.json",
                       "label": {"en": f"Radar LOS, {'morning' if tag == 'asc' else 'evening'} pass",
                                 "uk": f"Радар, {'ранковий' if tag == 'asc' else 'вечірній'} прохід"},
                       "props": g["columns"], "visible": tag == "asc"})
        if tag == "asc":
            epsg_asc = epsg
        print(f"  {tag}: {len(g['cells'])} cells, rms obs {g['rms_obs_cm']} cm, residual {g['rms_res_cm']} cm")
    rect, top = fault_outline(q, epi, epsg_asc)
    fc = {"type": "FeatureCollection", "features": [
        {"type": "Feature", "properties": {"part": "plane", "strike": round(q["strike"]), "dip": round(q["dip"])},
         "geometry": {"type": "Polygon", "coordinates": [rect]}},
        {"type": "Feature", "properties": {"part": "top_edge"}, "geometry": {"type": "LineString", "coordinates": top}},
        {"type": "Feature", "properties": {"part": "epicentre"}, "geometry": {"type": "Point", "coordinates": epi}}]}
    (out / "fault.geojson").write_text(json.dumps(fc, separators=(",", ":")), encoding="utf-8")
    thr = qc.get("conf_threshold", 0.75)
    card = {"kind": "card", "event": qc.get("event", cfg.get("title")), "date": qc.get("date"),
            "estimate": {"strike": round(q["strike"]), "dip": round(q["dip"]), "rake": round(q["rake"]),
                         "Mw": round(q["Mw"], 2), "length_km": round(q["L"], 1), "width_km": round(q["W"], 1),
                         "slip_m": round(q["slip"], 2), "centroid_depth_km": round(q["depth"], 1)},
            "confidence": round(conf, 2), "threshold": thr, "confident": bool(conf >= thr), "spread": spread,
            "relative_residual": round(rel, 2),
            "method": "Q-2look v1: neural proposals on the morning pass + exact slip solve on both passes; "
                      "confidence = top-16 agreement (15 deg strike & dip), median over seeds; shown run = lowest residual of the seeds",
            "provenance": qc.get("provenance", "")}
    ref = qc.get("reference")
    if ref:
        best = min(ref["planes"], key=lambda pl: angdiff(pl[0], q["strike"]) + abs(pl[1] - q["dip"]))
        card["reference"] = {"source": ref.get("source", "USGS"), "Mw": ref.get("mw"), "plane": best,
                             "d_strike": round(angdiff(best[0], q["strike"])), "d_dip": round(abs(best[1] - q["dip"])),
                             "note": "comparison only; not used by the estimate"}
    (out / "card.json").write_text(json.dumps(card, indent=1, ensure_ascii=False), encoding="utf-8")
    print(f"  card: {card['estimate']} conf {conf:.2f} -> {'confident' if card['confident'] else 'NOT confident'}")
    cfg["extra_layers"] = layers + [
        {"id": "fault", "kind": "polygons", "file": "fault.geojson", "style": "fault",
         "label": {"en": "Fault estimate", "uk": "Оцінка розлому"}, "props": ["part", "strike", "dip"]},
        {"id": "card", "kind": "card", "file": "card.json", "label": {"en": "Quake card", "uk": "Картка землетрусу"}}]
    man = engine.export_site(cfg, None, out)
    idx_p = Path(a.out) / "scenes.json"
    idx = json.loads(idx_p.read_text(encoding="utf-8")) if idx_p.exists() else {"scenes": []}
    entry = {"id": cfg["id"], "title": cfg.get("title", cfg["id"]), "subtitle": cfg.get("subtitle", ""), "kind": "quake"}
    ids = [s["id"] for s in idx["scenes"]]
    if cfg["id"] in ids:                       # rebuilding keeps the scene's place in the menu
        old = idx["scenes"][ids.index(cfg["id"])]   # keep keys the build does not write (blurb, href)
        idx["scenes"][ids.index(cfg["id"])] = {**old, **entry}
    else:
        idx["scenes"].append(entry)
    idx_p.write_text(json.dumps(idx, indent=1, ensure_ascii=False), encoding="utf-8")
    print(man)


if __name__ == "__main__":
    main()
