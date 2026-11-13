"""D-2look-lite (pre-registered): asc network proposals, slip solve + polish on asc AND desc jointly.

    SBI_BLOCK=24 SBI_LMAX=200e3 SBI_XYMAX=60e3 HYB_XYMAX=80e3 HYB_LMAX=300e3 \
        uv run --group ml python scripts/sbi_hybrid2.py --asc panama-wide --desc panama-wide-desc
"""
from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path

import torch

sys.path.insert(0, str(Path(__file__).parent))
from sbi_compare import IDS, judge, usgs  # noqa: E402
from sbi_fault import EVENTS, OUT, Net, Sim, dev, load_event, physical, sample, untarget  # noqa: E402
from sbi_hybrid import HI, ITERS, LO, LR, N_NET, TOP  # noqa: E402
from sbi_post3 import greens, solve  # noqa: E402

from hachure import okada


def set_heading(sim, heading):
    th = torch.acos(sim.lv[2])
    az = math.radians(heading + 90)
    sim.lv = torch.stack([torch.sin(th) * math.sin(az), torch.sin(th) * math.cos(az), torch.cos(th)])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--asc", required=True)
    ap.add_argument("--desc", required=True)
    a = ap.parse_args()
    torch.manual_seed(2)
    sA, sD = Sim(load_event(a.asc)), Sim(load_event(a.desc))
    set_heading(sD, EVENTS[a.desc]["heading"])
    sims = [sA, sD]
    w = [1 / math.sqrt(s.obs.numel()) for s in sims]
    ys = [s.deplane(s.obs[None]) for s in sims]
    net = Net().to(dev)
    net.load_state_dict(torch.load(OUT / f"sbi_{a.asc}_net.pt"))
    net.eval()

    def fwd(p):
        return [s.deplane(s.los(p)) for s in sims]

    def run(arm):
        torch.cuda.synchronize()
        t0 = time.time()
        with torch.no_grad():
            g = untarget(sample(net, sA.image(ys[0]), N_NET))
            top = g["depth"] - g["W"] / 2 * torch.sin(torch.deg2rad(g["dip"]))
            g["depth"] = torch.where(top < 0.1, g["depth"] - top + 0.1, g["depth"])
            p = {k: v * (1e3 if k in ("e0", "n0", "depth", "L", "W") else 1) for k, v in g.items()}
            ab = [greens(s, p) for s in sims]
            A = torch.cat([ab[i][0] * w[i] for i in range(2)], 1)
            B = torch.cat([ab[i][1] * w[i] for i in range(2)], 1)
            Y = torch.cat([ys[i] * w[i] for i in range(2)], 1).expand(A.shape[0], -1)
            s1, s2 = solve(A, B, Y)
            rms = (Y - (s1[:, None] * A + s2[:, None] * B)).pow(2).sum(1).nan_to_num(1e9)
            p["rake"] = torch.rad2deg(torch.atan2(s2, s1))
            p["slip"] = torch.hypot(s1, s2).clamp(0.03, 29)
            idx = torch.argsort(rms)[:TOP]
            p0 = {k: v[idx] for k, v in p.items()}
            if arm == "H2":
                pr = sA.prior(TOP)
                p0 = {k: torch.cat([p0[k], pr[k]]) for k in p0}
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
            return sum(torch.nan_to_num(((f - y) ** 2).mean(1), nan=0.0) for f, y in zip(fwd(q), ys))

        for _ in range(ITERS):
            loss = loss_of(unpack())
            opt.zero_grad()
            loss.sum().backward()
            for v in z.values():
                v.grad.nan_to_num_(0.0)
            opt.step()
        with torch.no_grad():
            q = unpack()
            L = loss_of(q)
            b = int(torch.argmin(L))
        torch.cuda.synchronize()
        best = {k: float(v[b]) for k, v in q.items()}
        out = physical(best)
        out["Mw"] = okada.moment_magnitude(best["L"], best["W"], best["slip"])
        return out, float(L[b].sqrt()), time.time() - t0

    run("H1")                                         # warm-up
    ref = usgs(IDS[a.asc])
    res = json.load(open(OUT / f"sbi_{a.asc}.json"))
    res["hybrid2look"] = {}
    for arm in ("H1", "H2"):
        q, rms, dt = run(arm)
        print(f"{arm}: {dt:.2f} s, rms {rms * 100:.2f} cm | " + ", ".join(f"{k} {v:.2f}" for k, v in q.items()))
        res["hybrid2look"][arm] = dict(params=q, rms_cm=rms * 100, seconds=dt, judge=judge(arm, q, ref))
    json.dump(res, open(OUT / f"sbi_{a.asc}.json", "w"), indent=1)


if __name__ == "__main__":
    main()
