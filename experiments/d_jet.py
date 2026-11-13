"""D-jet (pre-registered): Fisher bands of a deplaned fault map vs graded ranks of the jet map.

    uv run python scripts/d_jet.py
CPU, float64, noise-free. See notes/experiments.md (D-jet) and notes/math.md (M9).
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import torch

from hachure import okada

torch.set_default_dtype(torch.float64)
ROOT = Path(__file__).resolve().parents[1]
D = 60e3
KEYS = ["e0", "n0", "depth", "strike", "dip", "rake", "logL", "logW", "logslip"]
th, az = math.radians(35.0), math.radians(-12.0 + 90)
LV = (math.sin(th) * math.sin(az), math.sin(th) * math.cos(az), math.cos(th))
SLIP = 3.33e19 / (3e10 * 30e3 * 15e3)
V0 = torch.tensor([0.0, 0.0, D / 1e3, 192.0, 50.0, -119.0, math.log(30), math.log(15), math.log(SLIP)])


def field(v, x, y):
    """LOS displacement at points (x, y) [m] for natural-unit parameters v (9,)."""
    e, n = x - v[0] * 1e3, y - v[1] * 1e3
    uE, uN, uZ = okada.displacement(e[None], n[None], v[2:3] * 1e3, v[3:4], v[4:5], torch.exp(v[6:7]) * 1e3,
                                    torch.exp(v[7:8]) * 1e3, v[5:6], torch.exp(v[8:9]))
    return (uE * LV[0] + uN * LV[1] + uZ * LV[2])[0]


def mono(x, y, deg):
    return torch.stack([x ** (m - b) * y ** b for m in range(deg + 1) for b in range(m + 1)], 1)


def jac(x, y):
    return torch.autograd.functional.jacobian(lambda v: field(v, x, y), V0)         # (P, 9)


def grid(l, n):
    t = torch.linspace(-l, l, n)
    X, Y = torch.meshgrid(t, t, indexing="xy")
    return X.ravel(), Y.ravel()


def fisher_orders(k):
    svals = []
    for r in (0.05, 0.1, 0.2):
        l = r * D
        x, y = grid(l, 41)
        J = jac(x, y)
        A = mono(x / l, y / l, k)
        P = torch.eye(len(x)) - A @ torch.linalg.pinv(A)
        svals.append(torch.linalg.svdvals(P @ J).numpy())
    s = np.log(np.array(svals))
    slopes = np.polyfit(np.log([0.05, 0.1, 0.2]), s, 1)[0]
    return slopes


def jet_ranks(k, top):
    l = 0.05 * D
    x, y = grid(l, 61)
    J = jac(x, y)
    deg = k + 10
    A = mono(x / l, y / l, deg)                       # well-conditioned unit-window basis
    C = torch.linalg.lstsq(A, J, driver="gelsd").solution   # (n_mono, 9)
    sc = torch.cat([torch.full((m + 1,), (D / l) ** m) for m in range(deg + 1)])
    C = C * sc[:, None]                               # coefficients in units of d^m
    rows, ranks, start = [], [], 0
    blocks = {}
    for m in range(deg + 1):
        blocks[m] = C[start:start + m + 1]
        start += m + 1
    prev = 0
    out = {}
    for m in range(k + 1, top + 1):
        rows.append(blocks[m])
        S = torch.linalg.svdvals(torch.cat(rows))
        rk = int((S > 1e-7 * S[0]).sum())
        out[m] = rk - prev
        prev = rk
    return out


def main():
    res = {}
    for k in (0, 1, 2):
        sl = fisher_orders(k)
        ordA = sorted(int(round(v)) for v in sl)
        rB = jet_ranks(k, k + 6)
        ordB = sorted(m for m, r in rB.items() for _ in range(r))[:9]
        tot = sum(rB.values())
        reach = next((m for m in sorted(rB) if sum(rB[j] for j in rB if j <= m) >= 9), None)
        j1 = ordA == ordB
        j2 = reach is not None and reach <= k + 5
        print(f"k={k}: Fisher slopes {np.round(np.sort(sl), 2).tolist()}")
        print(f"      A orders {ordA}")
        print(f"      B incremental jet ranks {rB} (total {tot}) -> orders {ordB}")
        print(f"      J1 {'PASS' if j1 else 'FAIL'}; J2 rank 9 reached at order {reach} (<= {k + 5}) -> {'PASS' if j2 else 'FAIL'}")
        res[k] = dict(slopes=sl.tolist(), ordersA=ordA, jet_ranks=rB, ordersB=ordB, J1=j1, J2=j2, reach=reach)
    json.dump(res, open(ROOT / "data/processed/events/d_jet.json", "w"), indent=1, default=int)


if __name__ == "__main__":
    main()
