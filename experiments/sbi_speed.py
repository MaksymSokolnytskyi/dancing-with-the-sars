"""D-2look-speed (pre-registered): choose the two-look hybrid speed-up on synthetic events only.

    SBI_BLOCK=24 SBI_LMAX=200e3 SBI_XYMAX=60e3 HYB_XYMAX=80e3 HYB_LMAX=300e3 \
        uv run --group ml python scripts/sbi_speed.py --asc panama-wide --desc panama-wide-desc
"""
from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path

import numpy as np
import torch

sys.path[:0] = [str(Path(__file__).parent), str(Path(__file__).parents[1] / "scripts")]
from sbi_fault import EVENTS, OUT, Net, Sim, dev, load_event, physical, sample, untarget  # noqa: E402
from sbi_hybrid import HI, LO, LR  # noqa: E402
from sbi_hybrid2 import set_heading  # noqa: E402
from sbi_post3 import greens, solve  # noqa: E402

from hachure import okada

SETTINGS = [(2000, 16, 60), (2000, 16, 30), (1000, 8, 30), (1000, 8, 15), (500, 8, 15), (2000, 16, 0)]


def hybrid(net, sims, ys, n_net, top, iters):
    w = [1 / math.sqrt(s.obs.numel()) for s in sims]
    torch.cuda.synchronize()
    t0 = time.time()
    with torch.no_grad():
        g = untarget(sample(net, sims[0].image(ys[0]), n_net))
        tp = g["depth"] - g["W"] / 2 * torch.sin(torch.deg2rad(g["dip"]))
        g["depth"] = torch.where(tp < 0.1, g["depth"] - tp + 0.1, g["depth"])
        p = {k: v * (1e3 if k in ("e0", "n0", "depth", "L", "W") else 1) for k, v in g.items()}
        ab = [greens(s, p) for s in sims]
        A = torch.cat([ab[i][0] * w[i] for i in range(2)], 1)
        B = torch.cat([ab[i][1] * w[i] for i in range(2)], 1)
        Y = torch.cat([ys[i] * w[i] for i in range(2)], 1).expand(A.shape[0], -1)
        s1, s2 = solve(A, B, Y)
        rms = (Y - (s1[:, None] * A + s2[:, None] * B)).pow(2).sum(1).nan_to_num(1e9)
        p["rake"] = torch.rad2deg(torch.atan2(s2, s1))
        p["slip"] = torch.hypot(s1, s2).clamp(0.03, 29)
        idx = torch.argsort(rms)[:top]
        p0 = {k: v[idx] for k, v in p.items()}
    p0["ztop"] = (p0.pop("depth") - p0["W"] / 2 * torch.sin(torch.deg2rad(p0["dip"]))).clamp(min=100)
    p0["strike"] = p0["strike"] % 360
    z = {k: torch.logit(((v - LO[k]) / (HI[k] - LO[k])).clamp(1e-4, 1 - 1e-4)).clone().requires_grad_(True)
         for k, v in p0.items()}
    opt = torch.optim.Adam(z.values(), lr=LR)

    def unpack():
        q = {k: LO[k] + (HI[k] - LO[k]) * torch.sigmoid(v) for k, v in z.items()}
        q["depth"] = q.pop("ztop") + q["W"] / 2 * torch.sin(torch.deg2rad(q["dip"]))
        return q

    def loss_of(q):
        return sum(torch.nan_to_num(((s.deplane(s.los(q)) - y) ** 2).mean(1), nan=0.0) for s, y in zip(sims, ys))

    for _ in range(iters):
        loss = loss_of(unpack())
        opt.zero_grad()
        loss.sum().backward()
        for v in z.values():
            v.grad.nan_to_num_(0.0)
        opt.step()
    with torch.no_grad():
        q = unpack()
        b = int(torch.argmin(loss_of(q)))
    torch.cuda.synchronize()
    best = {k: float(v[b]) for k, v in q.items()}
    out = physical(best)
    out["Mw"] = okada.moment_magnitude(best["L"], best["W"], best["slip"])
    return out, time.time() - t0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--asc", required=True)
    ap.add_argument("--desc", required=True)
    ap.add_argument("--n", type=int, default=24)
    a = ap.parse_args()
    torch.manual_seed(7)
    sA, sD = Sim(load_event(a.asc)), Sim(load_event(a.desc))
    set_heading(sD, EVENTS[a.desc]["heading"])
    sims = [sA, sD]
    net = Net().to(dev)
    net.load_state_dict(torch.load(OUT / f"sbi_{a.asc}_net.pt"))
    net.eval()
    truths = []
    while len(truths) < a.n:
        p = sA.prior(4096)
        mw = (2 / 3) * (torch.log10(3e10 * p["L"] * p["W"] * p["slip"]) - 9.1)
        keep = (mw >= 7.0) & (mw <= 7.9)
        for i in torch.nonzero(keep).ravel()[: a.n - len(truths)].tolist():
            truths.append({k: v[i:i + 1] for k, v in p.items()})
    hybrid(net, sims, [s.deplane(s.obs[None]) for s in sims], 500, 8, 2)        # warm-up only, not judged
    res = {str(s): [] for s in SETTINGS}
    for j, tp in enumerate(truths):
        ys = []
        for s in sims:
            k = torch.randint(0, len(s.noise), (1,), device=dev)
            sign = 1.0 if torch.rand(1).item() < 0.5 else -1.0
            ys.append(s.deplane(s.los(tp) + sign * s.noise[k] + s.turbulence(1)))
        tq = physical({k: float(v) for k, v in tp.items()})
        tmw = okada.moment_magnitude(float(tp["L"]), float(tp["W"]), float(tp["slip"]))
        marks = []
        for st in SETTINGS:
            q, dt = hybrid(net, sims, ys, *st)
            dst = abs((q["strike"] - tq["strike"] + 180) % 360 - 180)
            ddip = abs(q["dip"] - tq["dip"])
            dmw = abs(q["Mw"] - tmw)
            ok = bool(dst <= 15 and ddip <= 15 and dmw <= 0.2)
            res[str(st)].append(dict(dstrike=dst, ddip=ddip, dmw=dmw, t=dt, ok=ok))
            marks.append(f"{st[2]}it/{st[0]}:" + ("ok" if ok else "--"))
        print(f"event {j + 1}/{len(truths)}: Mw {tmw:.2f} strike {tq['strike']:.0f} dip {tq['dip']:.0f} | " + " ".join(marks),
              flush=True)
    summ = {}
    for st in SETTINGS:
        r = res[str(st)]
        summ[str(st)] = dict(success=float(np.mean([x["ok"] for x in r])),
                             med_dstrike=float(np.median([x["dstrike"] for x in r])),
                             med_ddip=float(np.median([x["ddip"] for x in r])),
                             med_dmw=float(np.median([x["dmw"] for x in r])),
                             med_t=float(np.median([x["t"] for x in r])))
        s = summ[str(st)]
        print(f"N_NET {st[0]:5d} TOP {st[1]:2d} ITERS {st[2]:2d}: success {s['success']:.2f}, "
              f"med |dstrike| {s['med_dstrike']:.1f}, |ddip| {s['med_ddip']:.1f}, |dMw| {s['med_dmw']:.3f}, "
              f"time {s['med_t']:.2f} s")
    full = summ[str(SETTINGS[0])]
    ok = [st for st in SETTINGS if summ[str(st)]["success"] >= full["success"] - 0.10
          and summ[str(st)]["med_dmw"] <= full["med_dmw"] + 0.05]
    pick = min(ok, key=lambda st: summ[str(st)]["med_t"])
    nop = summ[str(SETTINGS[-1])]
    print(f"selected {pick}: {summ[str(pick)]['med_t']:.2f} s -> S1 {'PASS' if summ[str(pick)]['med_t'] < 10 else 'FAIL'}")
    print(f"S2: no-polish success {nop['success']:.2f} vs full {full['success']:.2f} -> "
          f"{'PASS' if full['success'] - nop['success'] >= 0.2 else 'FAIL'}")
    json.dump(dict(summary=summ, picked=list(pick), raw=res), open(OUT / "sbi_speed.json", "w"), indent=1)


if __name__ == "__main__":
    main()
