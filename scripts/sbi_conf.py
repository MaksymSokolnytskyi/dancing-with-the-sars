"""D-conf (pre-registered): does agreement among the top proposals flag a wrong quake card? Plus a rake check.

    SBI_BLOCK=24 SBI_LMAX=200e3 SBI_XYMAX=60e3 HYB_XYMAX=80e3 HYB_LMAX=300e3 \
        uv run --group ml python scripts/sbi_conf.py --asc panama-wide --desc panama-wide-desc
Q-2look v1 (N_NET 2000, TOP 16, no polish), synthetic events only.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
import torch

sys.path.insert(0, str(Path(__file__).parent))
from sbi_fault import EVENTS, OUT, Net, Sim, dev, load_event, physical, sample, untarget  # noqa: E402
from sbi_hybrid2 import set_heading  # noqa: E402
from sbi_post3 import greens, solve  # noqa: E402

from hachure import okada


def angdiff(a, b, period=360.0):
    return abs((a - b + period / 2) % period - period / 2)


def auc(score, label):
    """P(score of a success > score of a failure), ties count 1/2."""
    s, f = score[label], score[~label]
    if len(s) == 0 or len(f) == 0:
        return float("nan")
    return float(((s[:, None] > f[None]).sum() + 0.5 * (s[:, None] == f[None]).sum()) / (len(s) * len(f)))


def q2look(net, sims, ys, n_net=2000, top=16):
    w = [1 / math.sqrt(s.obs.numel()) for s in sims]
    with torch.no_grad():
        g = untarget(sample(net, sims[0].image(ys[0]), n_net))
        tp = g["depth"] - g["W"] / 2 * torch.sin(torch.deg2rad(g["dip"]))
        g["depth"] = torch.where(tp < 0.1, g["depth"] - tp + 0.1, g["depth"])
        p = {k: v * (1e3 if k in ("e0", "n0", "depth", "L", "W") else 1) for k, v in g.items()}
        ab = [greens(s, p) for s in sims]
        A = torch.cat([ab[i][0] * w[i] for i in range(2)], 1)
        B = torch.cat([ab[i][1] * w[i] for i in range(2)], 1)
        Yd = torch.cat([ys[i] * w[i] for i in range(2)], 1)
        Y = Yd.expand(A.shape[0], -1)
        s1, s2 = solve(A, B, Y)
        rms = (Y - (s1[:, None] * A + s2[:, None] * B)).pow(2).sum(1).nan_to_num(1e9)
        p["rake"] = torch.rad2deg(torch.atan2(s2, s1))
        p["slip"] = torch.hypot(s1, s2).clamp(0.03, 29)
        idx = torch.argsort(rms)[:top]
        cand = [physical({k: float(v[i]) for k, v in p.items()}) for i in idx.tolist()]
        for c, i in zip(cand, idx.tolist()):
            c["Mw"] = okada.moment_magnitude(float(p["L"][i]), float(p["W"][i]), float(p["slip"][i]))
        rel = float(rms[idx[0]].sqrt() / Yd.pow(2).sum().sqrt())
    best = cand[0]
    conf = float(np.mean([angdiff(c["strike"], best["strike"]) <= 15 and abs(c["dip"] - best["dip"]) <= 15 for c in cand]))
    return best, conf, rel


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--asc", required=True)
    ap.add_argument("--desc", required=True)
    ap.add_argument("--n", type=int, default=48)
    a = ap.parse_args()
    torch.manual_seed(11)
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
        for i in torch.nonzero((mw >= 7.0) & (mw <= 7.9)).ravel()[: a.n - len(truths)].tolist():
            truths.append({k: v[i:i + 1] for k, v in p.items()})
    rows = []
    for j, tp in enumerate(truths):
        ys = []
        for s in sims:
            k = torch.randint(0, len(s.noise), (1,), device=dev)
            sign = 1.0 if torch.rand(1).item() < 0.5 else -1.0
            ys.append(s.deplane(s.los(tp) + sign * s.noise[k] + s.turbulence(1)))
        tq = physical({k: float(v) for k, v in tp.items()})
        tmw = okada.moment_magnitude(float(tp["L"]), float(tp["W"]), float(tp["slip"]))
        q, conf, rel = q2look(net, sims, ys)
        dst, ddip, dmw = angdiff(q["strike"], tq["strike"]), abs(q["dip"] - tq["dip"]), abs(q["Mw"] - tmw)
        ok = bool(dst <= 15 and ddip <= 15 and dmw <= 0.2)
        rows.append(dict(conf=conf, rel=rel, ok=ok, dstrike=dst, ddip=ddip, dmw=dmw, drake=angdiff(q["rake"], tq["rake"]),
                         t_strike=tq["strike"], t_dip=tq["dip"], t_rake=tq["rake"], t_mw=tmw))
        print(f"event {j + 1}/{len(truths)}: {'ok' if ok else '--'} conf {conf:.2f} rel {rel:.2f} | dstrike {dst:.0f} "
              f"ddip {ddip:.0f} dMw {dmw:.2f} drake {rows[-1]['drake']:.0f}", flush=True)
    c = np.array([r["conf"] for r in rows])
    rl = np.array([r["rel"] for r in rows])
    ok = np.array([r["ok"] for r in rows])
    f1 = auc(c, ok)
    hi = c >= 0.5
    f2 = float(ok[hi].mean()) if hi.sum() else float("nan")
    f3 = auc(-rl, ok)
    dr = np.array([r["drake"] for r in rows])[ok]
    r1 = float(np.median(dr)) if len(dr) else float("nan")
    print(f"success {ok.mean():.2f} ({ok.sum()}/{len(ok)})")
    print(f"F1: AUC(conf) {f1:.2f} (>= 0.75) -> {'PASS' if f1 >= 0.75 else 'FAIL'}")
    print(f"F2: success when conf >= 0.5: {f2:.2f} on {hi.sum()} events (>= 0.70, n >= 10) -> "
          f"{'PASS' if hi.sum() >= 10 and f2 >= 0.7 else 'FAIL'}")
    print(f"F3: AUC(-relative residual) {f3:.2f} (>= 0.70) -> {'PASS' if f3 >= 0.7 else 'FAIL'}")
    print(f"R1: median |drake| among successes {r1:.1f} deg (<= 20) -> {'PASS' if r1 <= 20 else 'FAIL'}")
    json.dump(dict(rows=rows, F1=f1, F2=f2, n_hi=int(hi.sum()), F3=f3, R1=r1), open(OUT / "sbi_conf.json", "w"), indent=1)


if __name__ == "__main__":
    main()
