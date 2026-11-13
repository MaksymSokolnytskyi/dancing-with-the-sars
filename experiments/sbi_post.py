"""D, post-processing of the frozen networks (still blind to USGS): circular summaries for strike AND
rake, misfit of the summary model, and the best-fitting of 4000 posterior samples.

    uv run --group ml python scripts/sbi_post.py --event nagqu
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
import torch

sys.path[:0] = [str(Path(__file__).parent), str(Path(__file__).parents[1] / "scripts")]
from sbi_fault import OUT, Net, Sim, dev, load_event, sample, untarget  # noqa: E402

from hachure import okada

KM = ("e0", "n0", "depth", "L", "W")


def circ(v):
    c = math.degrees(math.atan2(np.mean(np.sin(np.radians(v))), np.mean(np.cos(np.radians(v)))))
    d = (v - c + 180) % 360 - 180
    return c, float(c + np.percentile(d, 16)), float(c + np.percentile(d, 84))


def to_torch(p):
    return {k: torch.tensor(np.atleast_1d(p[k]) * (1e3 if k in KM else 1), dtype=torch.float32, device=dev)
            for k in ("e0", "n0", "depth", "strike", "dip", "rake", "L", "W", "slip")}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--event", required=True)
    args = ap.parse_args()
    torch.manual_seed(1)
    ev = load_event(args.event)
    sim = Sim(ev)
    net = Net().to(dev)
    net.load_state_dict(torch.load(OUT / f"sbi_{args.event}_net.pt"))
    net.eval()
    yo = sim.deplane(sim.obs[None])
    s = {k: v.cpu().numpy() for k, v in untarget(sample(net, sim.image(yo), 20000)).items()}
    post = {}
    for k, v in s.items():
        if k in ("strike", "rake"):
            post[k] = list(circ(v))
        else:
            post[k] = [float(np.median(v)), float(np.percentile(v, 16)), float(np.percentile(v, 84))]
    post["strike"][0] %= 360
    mw = np.array([okada.moment_magnitude(a * 1e3, b * 1e3, c) for a, b, c in zip(s["L"], s["W"], s["slip"])])
    post["Mw"] = [float(np.median(mw)), float(np.percentile(mw, 16)), float(np.percentile(mw, 84))]
    rk = s["rake"]
    modes = {"left-lat (|r|<45)": np.mean(np.abs(rk) < 45), "thrust (45..135)": np.mean((rk >= 45) & (rk <= 135)),
             "right-lat (|r|>135)": np.mean(np.abs(rk) > 135), "normal (-135..-45)": np.mean((rk <= -45) & (rk >= -135))}
    print("rake posterior mass: " + ", ".join(f"{k} {v * 100:.0f}%" for k, v in modes.items()))
    with torch.no_grad():
        summ = {k: post[k][0] for k in ("e0", "n0", "depth", "strike", "dip", "rake", "L", "W", "slip")}
        rms_med = float(((sim.deplane(sim.los(to_torch(summ)))[0] - yo[0]) ** 2).mean().sqrt())
        idx = np.random.default_rng(0).choice(len(rk), 4000, replace=False)
        sub = {k: s[k][idx].copy() for k in summ}
        top = sub['depth'] - sub['W'] / 2 * np.sin(np.radians(sub['dip']))
        sub['depth'] = np.where(top < 0.1, sub['depth'] - top + 0.1, sub['depth'])
        for k in summ:
            s[k][idx] = sub[k]
        P = to_torch(sub)
        r = torch.cat([((sim.deplane(sim.los({k: v[i:i + 500] for k, v in P.items()})) - yo) ** 2).mean(1).sqrt()
                       for i in range(0, 4000, 500)]).cpu().numpy()
        r = np.where(np.isfinite(r), r, np.inf)
        b = int(np.argmin(r))
        rk_i = s['rake'][idx]
        for name, m in [('left-lat', np.abs(rk_i) < 45), ('thrust', (rk_i >= 45) & (rk_i <= 135)),
                        ('right-lat', np.abs(rk_i) > 135), ('normal', (rk_i <= -45) & (rk_i >= -135))]:
            if m.any():
                print(f'  mode {name:9s}: {m.sum():4d} samples, rms median {np.median(r[m]) * 100:5.1f} cm, best {r[m].min() * 100:5.1f} cm')
        best = {k: float(s[k][idx][b]) for k in summ}
        best["Mw"] = okada.moment_magnitude(best["L"] * 1e3, best["W"] * 1e3, best["slip"])
        model = sim.deplane(sim.los(to_torch(best)))[0]
    print("circular posterior summary (median [16, 84]):")
    for k, (m, lo, hi) in post.items():
        print(f"  {k:7s} {m:8.2f}   [{lo:8.2f}, {hi:8.2f}]")
    print(f"summary-model rms {rms_med * 100:.1f} cm | best of 4000 posterior samples rms {r[b] * 100:.1f} cm "
          f"(sample rms median {np.median(r) * 100:.1f}) | data std {float(yo.std()) * 100:.1f} cm")
    print("best posterior sample: " + ", ".join(f"{k} {v:.2f}" for k, v in best.items()))
    res = json.load(open(OUT / f"sbi_{args.event}.json"))
    res["sbi_v1_linear_rake_summary"] = res.pop("sbi")
    res["sbi"] = post
    res["sbi_rms_cm"] = rms_med * 100
    res["sbi_best_sample"] = best
    res["sbi_best_sample_rms_cm"] = float(r[b] * 100)
    res["rake_modes"] = {k: float(v) for k, v in modes.items()}
    json.dump(res, open(OUT / f"sbi_{args.event}.json", "w"), indent=1)
    np.savez_compressed(OUT / f"sbi_{args.event}_maps.npz", obs=yo[0].cpu().numpy(), sbi=model.cpu().numpy(),
                        mask=sim.mask, shape=np.array(sim.shape), E=ev["E"], N=ev["N"])


if __name__ == "__main__":
    main()
