"""D v3 post-processing (still blind to USGS): geometry from the network, slip vector solved exactly.

Why: the v1 and v2 networks both ignore the sign of the image (posterior mean of the slip components
correlates 0.02-0.07 with the truth on synthetics; flipping the image leaves every output unchanged), so
their rake is the prior. But Okada surface displacement is LINEAR in the slip vector:
    los = s_ss * G_ss(geometry) + s_ds * G_ds(geometry),  s_ss = slip cos(rake), s_ds = slip sin(rake).
So for every geometry sample from the network we solve the 2-unknown least-squares problem in closed form
(variable projection, Golub-Pereyra 1973), and get the slip uncertainty by pushing the same simulator noise
(real pre-event maps x flips x random scale + power-law turbulence) through the same 2x2 solve.

    uv run --group ml python scripts/sbi_post3.py --event nagqu
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

sys.path.insert(0, str(Path(__file__).parent))
from sbi_fault import OUT, Net, Sim, dev, load_event, sample, untarget  # noqa: E402

from hachure import okada

GEOM = ("e0", "n0", "depth", "strike", "dip", "L", "W")


def circ(v):
    c = math.degrees(math.atan2(np.mean(np.sin(np.radians(v))), np.mean(np.cos(np.radians(v)))))
    d = (v - c + 180) % 360 - 180
    return [c, float(c + np.percentile(d, 16)), float(c + np.percentile(d, 84))]


def greens(sim, g):
    """Deplaned LOS for unit strike-slip and unit dip-slip on each geometry: (B, P) x 2."""
    p = {k: g[k] for k in GEOM}
    one = torch.ones_like(g["dip"])
    a = sim.deplane(sim.los({**p, "rake": 0 * one, "slip": one}))
    b = sim.deplane(sim.los({**p, "rake": 90 * one, "slip": one}))
    return a, b


def solve(a, b, y):
    """Least-squares (s_ss, s_ds) for y ~ s_ss a + s_ds b, batched over rows."""
    aa, bb, ab = (a * a).sum(1), (b * b).sum(1), (a * b).sum(1)
    ay, by = (a * y).sum(1), (b * y).sum(1)
    det = (aa * bb - ab * ab).clamp(min=1e-12)
    return (bb * ay - ab * by) / det, (aa * by - ab * ay) / det


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--event", required=True)
    ap.add_argument("--n", type=int, default=20000)
    args = ap.parse_args()
    torch.manual_seed(1)
    ev = load_event(args.event)
    sim = Sim(ev)
    net = Net().to(dev)
    net.load_state_dict(torch.load(OUT / f"sbi_{args.event}_net.pt"))
    net.eval()
    yo = sim.deplane(sim.obs[None])
    torch.cuda.synchronize()
    t0 = time.time()
    with torch.no_grad():
        g = untarget(sample(net, sim.image(yo), args.n))
        top = g["depth"] - g["W"] / 2 * torch.sin(torch.deg2rad(g["dip"]))
        g["depth"] = torch.where(top < 0.1, g["depth"] - top + 0.1, g["depth"])     # no faults above ground
        gk = {k: v * (1e3 if k in ("e0", "n0", "depth", "L", "W") else 1) for k, v in g.items()}
        ss, ds, rms = [], [], []
        for i in range(0, args.n, 1000):
            sl = {k: v[i:i + 1000] for k, v in gk.items()}
            a, b = greens(sim, sl)
            B = a.shape[0]
            # exact solve on the data + one draw of simulator noise pushed through the same solve
            k = torch.randint(0, len(sim.noise), (B,), device=dev)
            scale = (0.7 + 0.8 * torch.rand(B, 1, device=dev)) * torch.where(torch.rand(B, 1, device=dev) < 0.5, -1.0, 1.0)
            noise = sim.deplane(scale * sim.noise[k] + sim.turbulence(B))
            s1, s2 = solve(a, b, yo.expand(B, -1))
            n1, n2 = solve(a, b, noise)
            r = yo - (s1[:, None] * a + s2[:, None] * b)
            ss.append(s1 + n1)
            ds.append(s2 + n2)
            rms.append(r.pow(2).mean(1).sqrt())
        ss, ds, rms = torch.cat(ss), torch.cat(ds), torch.cat(rms)
    torch.cuda.synchronize()
    inf_s = time.time() - t0
    s = {k: v.cpu().numpy() for k, v in g.items()}
    s["rake"] = np.degrees(np.arctan2(ds.cpu().numpy(), ss.cpu().numpy()))
    s["slip"] = np.hypot(ss.cpu().numpy(), ds.cpu().numpy())
    rms = rms.cpu().numpy()
    rms = np.where(np.isfinite(rms), rms, np.inf)
    mw = np.array([okada.moment_magnitude(a * 1e3, b * 1e3, c) for a, b, c in zip(s["L"], s["W"], s["slip"])])
    post = {}
    for k, v in s.items():
        post[k] = circ(v) if k in ("strike", "rake") else [float(np.median(v)), float(np.percentile(v, 16)),
                                                           float(np.percentile(v, 84))]
    post["strike"] = [v % 360 for v in post["strike"]]
    post["Mw"] = [float(np.median(mw)), float(np.percentile(mw, 16)), float(np.percentile(mw, 84))]
    rk = s["rake"]
    modes = {"left-lat (|r|<45)": np.mean(np.abs(rk) < 45), "thrust (45..135)": np.mean((rk >= 45) & (rk <= 135)),
             "right-lat (|r|>135)": np.mean(np.abs(rk) > 135), "normal (-135..-45)": np.mean((rk <= -45) & (rk >= -135))}
    b = int(np.argmin(rms))
    best = {k: float(s[k][b]) for k in ("e0", "n0", "depth", "strike", "dip", "rake", "L", "W", "slip")}
    best["Mw"] = float(mw[b])
    print(f"{args.event}: {args.n} geometry samples + exact slip solve in {inf_s * 1e3:.0f} ms")
    print("rake mass: " + ", ".join(f"{k} {v * 100:.0f}%" for k, v in modes.items()))
    print("posterior (median [16, 84]):")
    for k, (m, lo, hi) in post.items():
        print(f"  {k:7s} {m:8.2f}   [{lo:8.2f}, {hi:8.2f}]")
    print(f"per-sample rms (data, best slip): median {np.median(rms) * 100:.1f} cm, best {rms[b] * 100:.1f} cm, "
          f"data std {float(yo.std()) * 100:.1f} cm")
    print("best sample: " + ", ".join(f"{k} {v:.2f}" for k, v in best.items()))
    res = json.load(open(OUT / f"sbi_{args.event}.json"))
    if "sbi_method" not in res:                     # keep the v2 network-only summary on re-runs
        res["sbi_v2_net_only"] = res.get("sbi")
    res["sbi"] = post
    res["sbi_method"] = "v3: network geometry + exact linear slip solve with simulator-noise push-through"
    res["inference_ms"] = inf_s * 1e3
    res["sbi_best_sample"] = best
    res["sbi_best_sample_rms_cm"] = float(rms[b] * 100)
    res["sbi_rms_median_cm"] = float(np.median(rms) * 100)
    res["rake_modes"] = {k: float(v) for k, v in modes.items()}
    json.dump(res, open(OUT / f"sbi_{args.event}.json", "w"), indent=1)


if __name__ == "__main__":
    main()
