"""D-ident (exploratory, not a pre-registered test): what can one radar map say about a fault?

A. Sloppiness: eigen-decomposition of the Fisher information J^T J of the Okada -> deplaned-LOS map
   (white 1 cm noise) in natural units, for the ascending look, a descending look, and both.
B. Mechanism confusion: for a Nagqu-like normal fault (USGS finite-fault geometry), the best fit with the
   rake FORCED into the strike-slip family, for asc only / desc only / both looks. Ratio of misfits says
   whether one look can tell the two mechanisms apart.

    uv run --group ml python scripts/d_ident.py
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import torch

sys.path[:0] = [str(Path(__file__).parent), str(Path(__file__).parents[1] / "scripts")]
from sbi_fault import Sim, dev, load_event  # noqa: E402

from hachure import okada

KEYS = ["e0", "n0", "depth", "strike", "dip", "rake", "logL", "logW", "logslip"]
UNIT = dict(e0="km", n0="km", depth="km", strike="rad", dip="rad", rake="rad", logL="ln", logW="ln", logslip="ln")


class Look:
    def __init__(self, sim, heading):
        self.sim = sim
        th = torch.acos(sim.lv[2])                              # per-pixel incidence from the asc vector
        az = math.radians(heading + 90)
        self.lv = torch.stack([torch.sin(th) * math.sin(az), torch.sin(th) * math.cos(az), torch.cos(th)])

    def __call__(self, p):
        s = self.sim
        e = s.e[None, :] - p["e0"][:, None]
        n = s.n[None, :] - p["n0"][:, None]
        uE, uN, uZ = okada.displacement(e, n, p["depth"], p["strike"], p["dip"], p["L"], p["W"], p["rake"], p["slip"])
        return s.deplane(uE * self.lv[0] + uN * self.lv[1] + uZ * self.lv[2])


def phys(v):
    """natural-unit vector (B, 9) -> Okada dict in metres / degrees."""
    return dict(e0=v[:, 0] * 1e3, n0=v[:, 1] * 1e3, depth=v[:, 2] * 1e3, strike=v[:, 3], dip=v[:, 4], rake=v[:, 5],
                L=torch.exp(v[:, 6]) * 1e3, W=torch.exp(v[:, 7]) * 1e3, slip=torch.exp(v[:, 8]))


def forward(looks, v):
    return torch.cat([lk(phys(v)) for lk in looks], 1)


def fisher(looks, v0):
    J = torch.autograd.functional.jacobian(lambda x: forward(looks, x[None])[0], v0)   # (P, 9), metres per unit
    J = J * torch.tensor([1, 1, 1, 57.2958, 57.2958, 57.2958, 1, 1, 1], device=J.device)  # angles per radian
    F = (J.T @ J) / 0.01 ** 2                                                            # 1 cm white noise
    w, V = torch.linalg.eigh(F.double())
    return w.cpu().numpy()[::-1], V.cpu().numpy()[:, ::-1]


def fit(looks, y, v0, rake_lo, rake_hi, iters=400):
    """Multi-start Adam with rake squeezed into [rake_lo, rake_hi]; other params free (depth top >= 0.1 km)."""
    B = v0.shape[0]
    z = v0.clone()
    r0 = ((z[:, 5] - rake_lo) / (rake_hi - rake_lo)).clamp(0.02, 0.98)
    z[:, 5] = torch.logit(r0)
    z = z.requires_grad_(True)
    opt = torch.optim.Adam([z], lr=0.03)

    def unpack(z):
        rake = rake_lo + (rake_hi - rake_lo) * torch.sigmoid(z[:, 5])
        dip = 5 + 84.5 * torch.sigmoid(z[:, 4] / 20)                          # dip in (5, 89.5)
        half = torch.exp(z[:, 7]) / 2 * torch.sin(torch.deg2rad(dip))
        depth = torch.maximum(z[:, 2], half + 0.1)                            # no fault above ground
        return torch.stack([z[:, 0], z[:, 1], depth, z[:, 3], dip, rake, z[:, 6], z[:, 7], z[:, 8]], 1)

    for _ in range(iters):
        r = forward(looks, unpack(z)) - y
        loss = torch.nan_to_num((r ** 2).mean(1), nan=1.0)
        opt.zero_grad()
        loss.sum().backward()
        z.grad.nan_to_num_(0.0)
        opt.step()
    with torch.no_grad():
        v = unpack(z)
        rms = torch.nan_to_num(((forward(looks, v) - y) ** 2).mean(1), nan=1e9).sqrt()
        b = int(torch.argmin(rms))
    return v[b], float(rms[b])


def main():
    torch.manual_seed(0)
    sim = Sim(load_event("nagqu"))
    asc, desc = Look(sim, -12.0), Look(sim, -168.0)
    combos = {"asc": [asc], "desc": [desc], "asc+desc": [asc, desc]}
    # Nagqu-like truth: USGS finite fault strike 192, dip 50, rake -119, M0 3.3e19; centroid ~9 km, 30 x 15 km
    slip = 3.33e19 / (3e10 * 30e3 * 15e3)
    truth = torch.tensor([0, 0, 9, 192, 50, -119, math.log(30), math.log(15), math.log(slip)], device=dev)
    print(f"truth: normal fault 192/50/-119, 30 x 15 km, slip {slip:.2f} m, Mw "
          f"{okada.moment_magnitude(30e3, 15e3, slip):.2f}")

    print("\nA. Fisher information (1 cm white noise): eigenvalues and the stiffest / sloppiest directions")
    for name, looks in combos.items():
        w, V = fisher(looks, truth.clone())
        sd = 1 / np.sqrt(np.maximum(w, 1e-30))
        print(f"  {name:9s} 1-sigma along eigen-directions: " + " ".join(f"{s:.2g}" for s in sd) +
              f"   (condition number {w[0] / max(w[-1], 1e-30):.1e})")
        for i in (0, 8):
            top = np.argsort(-np.abs(V[:, i]))[:4]
            print(f"     {'stiffest' if i == 0 else 'sloppiest'}: " +
                  " ".join(f"{V[j, i]:+.2f}*{KEYS[j]}" for j in top))

    print("\nB. Best strike-slip imitation of the normal fault (noise-free), misfit vs signal")
    starts = 32
    for name, looks in combos.items():
        y = forward(looks, truth[None])
        sig = float(y.std())
        for fam, lo, hi in (("left-lateral", -30.0, 30.0), ("right-lateral", 150.0, 210.0), ("normal (control)", -170.0, -10.0)):
            v0 = torch.stack([torch.empty(starts, device=dev).uniform_(-8, 8), torch.empty(starts, device=dev).uniform_(-8, 8),
                              torch.empty(starts, device=dev).uniform_(4, 14), torch.empty(starts, device=dev).uniform_(0, 360),
                              torch.empty(starts, device=dev).uniform_(-30, 30),  # logit-scaled dip
                              torch.empty(starts, device=dev).uniform_(lo, hi),
                              torch.empty(starts, device=dev).uniform_(2.5, 4), torch.empty(starts, device=dev).uniform_(1.8, 3.2),
                              torch.empty(starts, device=dev).uniform_(-0.5, 1.5)], 1)
            v, rms = fit(looks, y, v0, lo, hi)
            print(f"  {name:9s} {fam:16s}: rms {rms * 100:5.2f} cm = {rms / sig * 100:4.1f}% of signal ({sig * 100:.1f} cm) | "
                  f"strike {float(v[3]) % 360:.0f} dip {float(v[4]):.0f} rake {float(v[5]):.0f} "
                  f"Mw {okada.moment_magnitude(math.exp(float(v[6])) * 1e3, math.exp(float(v[7])) * 1e3, math.exp(float(v[8]))):.2f}",
                  flush=True)


if __name__ == "__main__":
    main()
