"""Okada (1985) surface displacement of a uniform-slip rectangular fault in an elastic half-space.

Batched PyTorch (works on GPU): parameters have shape (B,), observation points (P,) -> (B, P).
Convention (as Beauducel's okada85): east/north relative to the fault CENTROID (m), depth of the
centroid (m, positive down), strike clockwise from north, dip to the right of strike, rake
Aki-Richards (0 = left-lateral, 90 = thrust), slip in m. Validated against cutde in
scripts/check_okada.py.
"""
from __future__ import annotations

import math

import torch

EPS = 1e-9


def _R(xi, eta, q):
    return torch.sqrt(xi * xi + eta * eta + q * q)


def _atan_term(xi, eta, q, R):
    out = torch.atan(xi * eta / (q * R + EPS * torch.sign(q * R + EPS)))
    return torch.where(q.abs() > EPS, out, torch.zeros_like(out))


def _I(xi, eta, q, sd, cd, R, k):
    """Okada's I1..I5 for cos(dip) != 0 (k = mu/(lambda+mu) = 1-2nu)."""
    db = eta * sd - q * cd
    yb = eta * cd + q * sd
    X = torch.sqrt(xi * xi + q * q)
    i5 = k * 2 / cd * torch.atan((eta * (X + q * cd) + X * (R + X) * sd) / (xi * (R + X) * cd + EPS))
    i5 = torch.where(xi.abs() > EPS, i5, torch.zeros_like(i5))
    i4 = k / cd * (torch.log(R + db) - sd * torch.log(R + eta))
    i3 = k * (yb / (cd * (R + db)) - torch.log(R + eta)) + sd / cd * i4
    i2 = k * (-torch.log(R + eta)) - i3
    i1 = k * (-xi / (cd * (R + db))) - sd / cd * i5
    return i1, i2, i3, i4, i5, db, yb


def _f(xi, eta, q, sd, cd, k):
    R = _R(xi, eta, q)
    i1, i2, i3, i4, i5, db, yb = _I(xi, eta, q, sd, cd, R, k)
    at = _atan_term(xi, eta, q, R)
    ux_ss = xi * q / (R * (R + eta)) + at + i1 * sd
    uy_ss = yb * q / (R * (R + eta)) + q * cd / (R + eta) + i2 * sd
    uz_ss = db * q / (R * (R + eta)) + q * sd / (R + eta) + i4 * sd
    ux_ds = q / R - i3 * sd * cd
    uy_ds = yb * q / (R * (R + xi)) + cd * at - i1 * sd * cd
    uz_ds = db * q / (R * (R + xi)) + sd * at - i5 * sd * cd
    return ux_ss, uy_ss, uz_ss, ux_ds, uy_ds, uz_ds


def displacement(e, n, depth, strike, dip, L, W, rake, slip, nu=0.25):
    """e, n: (P,) or (B, P) metres. Others: (B,) tensors (degrees for angles). Returns uE, uN, uZ (B, P)."""
    if e.dim() == 1:
        e, n = e[None, :], n[None, :]
    col = lambda t: t[:, None]                                    # noqa: E731
    st, dp, rk = (col(torch.deg2rad(a)) for a in (strike, dip, rake))
    depth, L, W, slip = col(depth), col(L), col(W), col(slip)
    # avoid the cos(dip) = 0 singularity of the I-terms
    dp = torch.clamp(dp, max=math.radians(89.9))
    sd, cd = torch.sin(dp), torch.cos(dp)
    cs, ss = torch.cos(st), torch.sin(st)
    U1, U2 = torch.cos(rk) * slip, torch.sin(rk) * slip
    d = depth + sd * W / 2
    ec = e + cs * cd * W / 2
    nc = n - ss * cd * W / 2
    x = cs * nc + ss * ec + L / 2
    y = ss * nc - cs * ec + cd * W
    p = y * cd + d * sd
    q = y * sd - d * cd
    k = 1 - 2 * nu
    acc = [0.0] * 6
    for xi, eta, s in ((x, p, 1), (x, p - W, -1), (x - L, p, -1), (x - L, p - W, 1)):
        for i, v in enumerate(_f(xi, eta, q, sd, cd, k)):
            acc[i] = acc[i] + s * v
    ux = -(U1 * acc[0] + U2 * acc[3]) / (2 * math.pi)
    uy = -(U1 * acc[1] + U2 * acc[4]) / (2 * math.pi)
    uz = -(U1 * acc[2] + U2 * acc[5]) / (2 * math.pi)
    uE = ss * ux - cs * uy
    uN = cs * ux + ss * uy
    return uE, uN, uz


def moment_magnitude(L, W, slip, mu=3e10):
    m0 = mu * L * W * abs(slip)
    return (2 / 3) * (math.log10(m0) - 9.1)
