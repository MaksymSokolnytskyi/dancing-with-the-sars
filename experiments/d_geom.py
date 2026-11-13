"""D-geom (a), pre-registered in notes/experiments.md: Fisher spectrum of the Okada -> asc LOS map vs centroid depth.

    uv run --group ml python scripts/d_geom.py
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import torch

sys.path[:0] = [str(Path(__file__).parent), str(Path(__file__).parents[1] / "scripts")]
from d_ident import KEYS, Look, fisher  # noqa: E402
from sbi_fault import Sim, dev, load_event  # noqa: E402

DEPTHS = [10, 15, 20, 30, 40, 60, 80]


def main():
    sim = Sim(load_event("nagqu"))
    asc = Look(sim, -12.0)
    slip = 3.33e19 / (3e10 * 30e3 * 15e3)
    spec, vecs = [], []
    for d in DEPTHS:
        v0 = torch.tensor([0, 0, d, 192, 50, -119, math.log(30), math.log(15), math.log(slip)], device=dev)
        w, V = fisher([asc], v0)
        spec.append(w)
        vecs.append(V)
        print(f"depth {d:3d} km: eigenvalues " + " ".join(f"{x:.2e}" for x in w), flush=True)
    spec = np.array(spec)
    lo = [i for i, d in enumerate(DEPTHS) if d >= 30]
    ld = np.log(np.array(DEPTHS)[lo])
    slopes = np.array([np.polyfit(ld, np.log(np.maximum(spec[lo, k], 1e-300)), 1)[0] for k in range(9)])
    print("log-log slope vs depth (30-80 km), stiffest -> sloppiest: " + " ".join(f"{s:+.2f}" for s in slopes))
    g1 = np.median(slopes[:7]) - slopes[7:].max()
    print(f"G1: median slope of top 7 minus the shallower of the 2 sloppiest = {g1:.2f} (need >= 3) -> "
          f"{'PASS' if g1 >= 3 else 'FAIL'}")
    V = vecs[-1]
    ok = True
    for k in (7, 8):
        v = V[:, k]
        shape = v[6:9]
        w_share = float(np.sum(shape ** 2) / np.sum(v ** 2))
        msum = abs(float(shape.sum())) / float(np.linalg.norm(shape))
        ok &= (w_share >= 0.7) and (msum < 0.3)
        print(f"  eigvec {k} at 80 km: weight on logL/logW/logslip {w_share:.2f}, |sum|/norm {msum:.2f}; " +
              " ".join(f"{v[j]:+.2f}*{KEYS[j]}" for j in np.argsort(-np.abs(v))[:4]))
    print(f"G2 -> {'PASS' if ok else 'FAIL'}")


if __name__ == "__main__":
    main()
