"""D-hyb: network as an instant proposal, then a short least-squares polish from its best samples.

Arms (pre-registered in notes/experiments.md, settings fixed before running on real data):
  H1 (primary): 16 starts = the 16 best of 2000 network geometry samples (slip solved exactly, v3)
  H2:           the same 16 + 16 random prior starts (guards against a too-narrow network posterior)
Both: Adam, lr 0.05, 60 iterations (chosen on noisy synthetics only: same accuracy as 300 at lr 0.02, 5x faster), same forward model, parameterisation and bounds as the LSQ baseline.

    uv run --group ml python scripts/sbi_hybrid.py --event kumamoto
    uv run --group ml python scripts/sbi_hybrid.py --event kumamoto --timing   # synthetic scene only
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import torch

sys.path.insert(0, str(Path(__file__).parent))
from sbi_fault import OUT, Net, Sim, dev, load_event, physical, sample, untarget  # noqa: E402
from sbi_post3 import greens, solve  # noqa: E402

from hachure import okada

import os
N_NET, TOP = 2000, 16
ITERS, LR = int(os.environ.get("HYB_ITERS", 60)), float(os.environ.get("HYB_LR", 0.05))
_XY = float(os.environ.get("HYB_XYMAX", 40e3))     # D-3w diagnostic: 80 km
_LM = float(os.environ.get("HYB_LMAX", 150e3))     # D-3w diagnostic: 300 km
LO = dict(e0=-_XY, n0=-_XY, ztop=50, strike=-720, dip=5, rake=-720, L=2e3, W=1e3, slip=0.02)
HI = dict(e0=_XY, n0=_XY, ztop=30e3, strike=720, dip=89.5, rake=720, L=_LM, W=60e3, slip=30)


def proposals(net, sim, yo):
    g = untarget(sample(net, sim.image(yo), N_NET))
    top = g["depth"] - g["W"] / 2 * torch.sin(torch.deg2rad(g["dip"]))
    g["depth"] = torch.where(top < 0.1, g["depth"] - top + 0.1, g["depth"])
    p = {k: v * (1e3 if k in ("e0", "n0", "depth", "L", "W") else 1) for k, v in g.items()}
    a, b = greens(sim, p)
    s1, s2 = solve(a, b, yo.expand(a.shape[0], -1))
    rms = (yo - (s1[:, None] * a + s2[:, None] * b)).pow(2).mean(1).sqrt().nan_to_num(1e9)
    p["rake"] = torch.rad2deg(torch.atan2(s2, s1))
    p["slip"] = torch.hypot(s1, s2).clamp(0.03, 29)
    idx = torch.argsort(rms)[:TOP]
    return {k: v[idx] for k, v in p.items()}


def polish(sim, yo, p0):
    p0 = dict(p0)
    p0["ztop"] = (p0.pop("depth") - p0["W"] / 2 * torch.sin(torch.deg2rad(p0["dip"]))).clamp(min=100)
    p0["strike"] = p0["strike"] % 360
    z = {k: torch.logit(((v - LO[k]) / (HI[k] - LO[k])).clamp(1e-4, 1 - 1e-4)).clone().requires_grad_(True)
         for k, v in p0.items()}
    opt = torch.optim.Adam(z.values(), lr=LR)

    def unpack():
        p = {k: LO[k] + (HI[k] - LO[k]) * torch.sigmoid(v) for k, v in z.items()}
        p["depth"] = p.pop("ztop") + p["W"] / 2 * torch.sin(torch.deg2rad(p["dip"]))
        return p

    for _ in range(ITERS):
        r = sim.deplane(sim.los(unpack())) - yo
        loss = torch.nan_to_num((r ** 2).mean(1), nan=0.0)
        opt.zero_grad()
        loss.sum().backward()
        for v in z.values():
            v.grad.nan_to_num_(0.0)
        opt.step()
    with torch.no_grad():
        p = unpack()
        loss = torch.nan_to_num(((sim.deplane(sim.los(p)) - yo) ** 2).mean(1), nan=1e9)
        b = int(torch.argmin(loss))
    return {k: float(v[b]) for k, v in p.items()}, float(loss[b].sqrt())


def run(net, sim, yo, arm):
    torch.cuda.synchronize()
    t0 = time.time()
    with torch.no_grad():
        p0 = proposals(net, sim, yo)
    if arm == "H2":
        pr = sim.prior(TOP)
        p0 = {k: torch.cat([p0[k], pr[k]]) for k in p0}
    best, rms = polish(sim, yo, p0)
    torch.cuda.synchronize()
    q = physical(best)
    q["Mw"] = okada.moment_magnitude(best["L"], best["W"], best["slip"])
    return q, rms, time.time() - t0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--event", required=True)
    ap.add_argument("--timing", action="store_true", help="time both arms on one synthetic scene, no real data")
    args = ap.parse_args()
    torch.manual_seed(2)
    sim = Sim(load_event(args.event))
    net = Net().to(dev)
    net.load_state_dict(torch.load(OUT / f"sbi_{args.event}_net.pt"))
    net.eval()
    if args.timing:
        _, _, p = sim.batch(1)
        yo = sim.deplane(sim.los(p))
        for arm in ("H1", "H2", "H1", "H2"):            # second pass = warm GPU
            q, rms, dt = run(net, sim, yo, arm)
            print(f"synthetic {arm}: {dt:.2f} s, rms {rms * 100:.2f} cm, truth strike {float(p['strike'][0]):.0f} "
                  f"dip {float(p['dip'][0]):.0f} -> {q['strike']:.0f} / {q['dip']:.0f}")
        return
    yo = sim.deplane(sim.obs[None])
    run(net, sim, yo, "H1")                               # warm-up, discarded
    res = json.load(open(OUT / f"sbi_{args.event}.json"))
    res["hybrid"] = {}
    for arm in ("H1", "H2"):
        q, rms, dt = run(net, sim, yo, arm)
        print(f"{args.event} {arm}: {dt:.2f} s, rms {rms * 100:.1f} cm (data std {float(yo.std()) * 100:.1f}) | "
              + ", ".join(f"{k} {v:.2f}" for k, v in q.items()))
        res["hybrid"][arm] = dict(params=q, rms_cm=rms * 100, seconds=dt)
    json.dump(res, open(OUT / f"sbi_{args.event}.json", "w"), indent=1)


if __name__ == "__main__":
    main()
