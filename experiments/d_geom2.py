"""D-geom (b), pre-registered: Fisher bands vs ramp-removal degree k = 0, 1, 2 (forward-mode Jacobian).

    uv run --group ml python scripts/d_geom2.py
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import torch

sys.path[:0] = [str(Path(__file__).parent), str(Path(__file__).parents[1] / "scripts")]
from d_ident import KEYS, phys  # noqa: E402
from sbi_fault import Sim, dev, load_event  # noqa: E402

from hachure import okada

DEPTHS = [30, 40, 60, 80]
SCALE = torch.tensor([1, 1, 1, 57.2958, 57.2958, 57.2958, 1, 1, 1])   # angles per radian


def projector(sim, k):
    e, n = sim.e / 3e4, sim.n / 3e4
    cols = [torch.ones_like(e)]
    if k >= 1:
        cols += [e, n]
    if k >= 2:
        cols += [e * e, e * n, n * n]
    A = torch.stack(cols, 1)
    return A, torch.linalg.pinv(A)


def main():
    sim = Sim(load_event("nagqu"))
    th = torch.acos(sim.lv[2])
    az = math.radians(-12.0 + 90)
    lv = torch.stack([torch.sin(th) * math.sin(az), torch.sin(th) * math.cos(az), torch.cos(th)])
    slip = 3.33e19 / (3e10 * 30e3 * 15e3)
    out = {}
    for k in (0, 1, 2):
        A, Ap = projector(sim, k)

        def f(v):
            p = phys(v[None])
            e = sim.e[None, :] - p["e0"][:, None]
            n = sim.n[None, :] - p["n0"][:, None]
            uE, uN, uZ = okada.displacement(e, n, p["depth"], p["strike"], p["dip"], p["L"], p["W"], p["rake"], p["slip"])
            y = (uE * lv[0] + uN * lv[1] + uZ * lv[2])[0]
            return y - A @ (Ap @ y)

        spec, vecs = [], []
        for d in DEPTHS:
            v0 = torch.tensor([0, 0, d, 192, 50, -119, math.log(30), math.log(15), math.log(slip)], device=dev,
                              dtype=torch.float32)
            J = torch.autograd.functional.jacobian(f, v0, vectorize=True, strategy="forward-mode")
            J = (J * SCALE.to(dev)).double()
            w, V = torch.linalg.eigh(J.T @ J / 0.01 ** 2)
            spec.append(w.cpu().numpy()[::-1])
            vecs.append(V.cpu().numpy()[:, ::-1])
        spec = np.array(spec)
        ld = np.log(DEPTHS)
        slopes = np.array([np.polyfit(ld, np.log(np.maximum(spec[:, j], 1e-300)), 1)[0] for j in range(9)])
        steep = int((slopes < -7).sum())
        V = vecs[-1]
        pos = float(np.mean([np.sum(V[:2, j] ** 2) / np.sum(V[:, j] ** 2) for j in (7, 8)]))
        out[k] = (steep, pos)
        print(f"k={k}: slopes " + " ".join(f"{s:+.2f}" for s in slopes) + f" | steep band {steep} | "
              f"(e0,n0) weight in 2 sloppiest at 80 km {pos:.2f}", flush=True)
        for j in (6, 7, 8):
            vj = V[:, j]
            print(f"     eig {j}: " + " ".join(f"{vj[i]:+.2f}*{KEYS[i]}" for i in np.argsort(-np.abs(vj))[:4]))
    q1 = out[0][0] < out[1][0] and out[0][1] < 0.3
    q2 = out[2][0] > out[1][0]
    print(f"Q1 {'PASS' if q1 else 'FAIL'}  Q2 {'PASS' if q2 else 'FAIL'}")


if __name__ == "__main__":
    main()
