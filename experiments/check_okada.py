"""Validate src/hachure/okada.py against cutde (Nikkhoo & Walter triangular dislocations, half-space)
and against basic physics (left-lateral sense, hanging-wall uplift).

    uv run python scripts/check_okada.py
"""
from __future__ import annotations

import itertools

import cutde.halfspace as hs
import numpy as np
import torch

from hachure import okada


def okada_np(e, n, **p):
    t = {k: torch.tensor([v], dtype=torch.float64) for k, v in p.items()}
    u = okada.displacement(torch.tensor(e), torch.tensor(n), **t)
    return np.stack([c[0].numpy() for c in u], 1)


def cutde_disp(e, n, depth, strike, dip, L, W, rake, slip, order):
    st, dp = np.radians(strike), np.radians(dip)
    s = np.array([np.sin(st), np.cos(st), 0.0])
    dd = np.array([np.cos(dp) * np.cos(st), -np.cos(dp) * np.sin(st), -np.sin(dp)])   # down-dip, right of strike
    c = np.array([0, 0, -depth])
    A, B = c - L / 2 * s - W / 2 * dd, c + L / 2 * s - W / 2 * dd       # top edge
    C, D = c + L / 2 * s + W / 2 * dd, c - L / 2 * s + W / 2 * dd       # bottom edge
    tris = np.array([[A, B, C], [A, C, D]]) if order == 0 else np.array([[A, C, B], [A, D, C]])
    obs = np.c_[e, n, np.zeros_like(e)]
    slip_v = np.array([[slip * np.cos(np.radians(rake)), slip * np.sin(np.radians(rake)), 0]] * 2)
    return hs.disp_free(obs, tris, slip_v, 0.25)


def main():
    rng = np.random.default_rng(0)
    e, n = rng.uniform(-40e3, 40e3, 3000), rng.uniform(-40e3, 40e3, 3000)
    worst = 0.0
    for strike, dip, rake in itertools.product((10, 135, 250), (20, 50, 80), (0, 90, -60, 180)):
        p = dict(depth=8e3, strike=float(strike), dip=float(dip), L=30e3, W=12e3, rake=float(rake), slip=2.0)
        uo = okada_np(e, n, **p)
        best = None
        for order in (0, 1):
            for sgn in (1, -1):
                uc = sgn * cutde_disp(e, n, order=order, **p)
                err = np.abs(uo - uc).max() / np.abs(uc).max()
                if best is None or err < best[0]:
                    best = (err, order, sgn)
        worst = max(worst, best[0])
        print(f"strike {strike:3d} dip {dip:2d} rake {rake:4d}: max rel err {best[0] * 100:.3f}% "
              f"(cutde vertex order {best[1]}, sign {best[2]:+d})")
    print(f"worst {worst * 100:.3f}%")
    # physics: vertical left-lateral N-striking fault -> point east of the fault moves north
    u = okada_np(np.array([2e3]), np.array([0.0]), depth=5e3, strike=0.0, dip=89.0, L=20e3, W=10e3, rake=0.0, slip=1.0)
    print(f"left-lateral, point 2 km east: uN = {u[0, 1]:+.3f} m (expect > 0)")
    # thrust dipping east (strike 0, dip to the right = east): uplift above the hanging wall (east side)
    u = okada_np(np.array([5e3, -5e3]), np.array([0.0, 0.0]), depth=6e3, strike=0.0, dip=30.0, L=20e3, W=10e3,
                 rake=90.0, slip=1.0)
    print(f"thrust dipping east: uZ east {u[0, 2]:+.3f}, west {u[1, 2]:+.3f} (expect east uplift > west)")


if __name__ == "__main__":
    main()
