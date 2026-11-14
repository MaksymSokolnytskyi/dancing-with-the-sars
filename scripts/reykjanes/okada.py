"""Surface displacement of a rectangular tensile (opening) dislocation in an elastic
half-space, after Okada (1985, BSSA 75, 1135), with a convenient centre-based
geometry. Units: km for geometry, any unit for opening (displacement comes out in
the same unit).

    enu = tensile(e, n, ec, nc, depth_c, length, width, strike_deg, dip_deg, opening, nu=0.25)

e, n: observation points (east, north, km); ec, nc, depth_c: centre of the plane
(depth positive down, km); strike clockwise from north; dip 0 = horizontal sill,
90 = vertical dike. Returns an array (..., 3) of east, north, up.

Run this file to execute the self-checks.
"""
import numpy as np

EPS = 1e-12


def _chinnery(f, x, p, L, W, q, sd, cd, a):
    return f(x, p, q, sd, cd, a) - f(x, p - W, q, sd, cd, a) - f(x - L, p, q, sd, cd, a) + f(x - L, p - W, q, sd, cd, a)


def _terms(xi, eta, q, sd, cd, a):
    """a = mu / (lambda + mu) = 1 - 2 nu ... here a = mu/(lambda+mu)."""
    R = np.sqrt(xi ** 2 + eta ** 2 + q ** 2)
    yt = eta * cd + q * sd
    dt = eta * sd - q * cd
    Re = R + eta
    Rx = R + xi
    if abs(cd) > 1e-8:
        X = np.sqrt(xi ** 2 + q ** 2)
        I5 = a * 2 / cd * np.arctan((eta * (X + q * cd) + X * (R + X) * sd) / (xi * (R + X) * cd + EPS))
        I4 = a / cd * (np.log(R + dt) - sd * np.log(Re))
        I3 = a * (yt / (cd * (R + dt)) - np.log(Re)) + sd / cd * I4
        I1 = a * (-xi / (cd * (R + dt))) - sd / cd * I5
    else:
        I5 = -a * xi * sd / (R + dt)
        I3 = a / 2 * (eta / (R + dt) + yt * q / (R + dt) ** 2 - np.log(Re))
        I1 = -a / 2 * xi * q / (R + dt) ** 2
    at = np.arctan(xi * eta / (q * R + EPS))
    ux = q ** 2 / (R * Re) - I3 * sd ** 2
    uy = -dt * q / (R * Rx) - sd * (xi * q / (R * Re) - at) - I1 * sd ** 2
    uz = yt * q / (R * Rx) + cd * (xi * q / (R * Re) - at) - I5 * sd ** 2
    return np.stack([ux, uy, uz], -1)


def _terms_ss_ds(xi, eta, q, sd, cd, a):
    """Okada (1985) surface terms for strike-slip and dip-slip (per unit slip, before the 1/2pi and the Chinnery sum)."""
    R = np.sqrt(xi ** 2 + eta ** 2 + q ** 2)
    yt = eta * cd + q * sd
    dt = eta * sd - q * cd
    Re = R + eta
    Rx = R + xi
    if abs(cd) > 1e-8:
        X = np.sqrt(xi ** 2 + q ** 2)
        I5 = a * 2 / cd * np.arctan((eta * (X + q * cd) + X * (R + X) * sd) / (xi * (R + X) * cd + EPS))
        I4 = a / cd * (np.log(R + dt) - sd * np.log(Re))
        I3 = a * (yt / (cd * (R + dt)) - np.log(Re)) + sd / cd * I4
        I1 = a * (-xi / (cd * (R + dt))) - sd / cd * I5
    else:
        I5 = -a * xi * sd / (R + dt)
        I4 = -a * q / (R + dt)
        I3 = a / 2 * (eta / (R + dt) + yt * q / (R + dt) ** 2 - np.log(Re))
        I1 = -a / 2 * xi * q / (R + dt) ** 2
    I2 = a * (-np.log(Re)) - I3
    at = np.arctan(xi * eta / (q * R + EPS))
    ss = np.stack([xi * q / (R * Re) + at + I1 * sd,
                   yt * q / (R * Re) + q * cd / Re + I2 * sd,
                   dt * q / (R * Re) + q * sd / Re + I4 * sd], -1)
    ds = np.stack([q / R - I3 * sd * cd,
                   yt * q / (R * Rx) + cd * at - I1 * sd * cd,
                   dt * q / (R * Rx) + sd * at - I5 * sd * cd], -1)
    return ss, ds


def _terms_all(xi, eta, q, sd, cd, a):
    """All three Okada surface solutions in one pass (shared R, I-terms): tensile, strike-slip, dip-slip -> (..., 9)."""
    R = np.sqrt(xi ** 2 + eta ** 2 + q ** 2)
    yt = eta * cd + q * sd
    dt = eta * sd - q * cd
    Re = R + eta
    Rx = R + xi
    lRe = np.log(Re)
    if abs(cd) > 1e-8:
        X = np.sqrt(xi ** 2 + q ** 2)
        I5 = a * 2 / cd * np.arctan((eta * (X + q * cd) + X * (R + X) * sd) / (xi * (R + X) * cd + EPS))
        I4 = a / cd * (np.log(R + dt) - sd * lRe)
        I3 = a * (yt / (cd * (R + dt)) - lRe) + sd / cd * I4
        I1 = a * (-xi / (cd * (R + dt))) - sd / cd * I5
    else:
        Rd = R + dt
        I5 = -a * xi * sd / Rd
        I4 = -a * q / Rd
        I3 = a / 2 * (eta / Rd + yt * q / Rd ** 2 - lRe)
        I1 = -a / 2 * xi * q / Rd ** 2
    I2 = -a * lRe - I3
    at = np.arctan(xi * eta / (q * R + EPS))
    qRRe, qRRx = q / (R * Re), q / (R * Rx)
    xqRRe = xi * qRRe
    return np.stack([q * qRRe - I3 * sd ** 2, -dt * qRRx - sd * (xqRRe - at) - I1 * sd ** 2, yt * qRRx + cd * (xqRRe - at) - I5 * sd ** 2,
                     xqRRe + at + I1 * sd, yt * qRRe + q * cd / Re + I2 * sd, dt * qRRe + q * sd / Re + I4 * sd,
                     q / R - I3 * sd * cd, yt * qRRx + cd * at - I1 * sd * cd, dt * qRRx + sd * at - I5 * sd * cd], -1)


def dislocation_fast(e, n, ec, nc, depth_c, length, width, strike_deg, dip_deg, strike_slip=0.0, dip_slip=0.0, opening=0.0, nu=0.25):
    """Same as dislocation(), one pass for all slip types (about 2x faster when opening and slip are both present)."""
    e, n = np.asarray(e, float), np.asarray(n, float)
    phi, dip = np.deg2rad(strike_deg), np.deg2rad(dip_deg)
    sd, cd = np.sin(dip), np.cos(dip)
    if abs(cd) < 1e-8:
        cd = 0.0
    s = np.array([np.sin(phi), np.cos(phi)])
    nrm = np.array([-np.cos(phi), np.sin(phi)])
    de, dn = e - ec, n - nc
    x = de * s[0] + dn * s[1] + length / 2
    y = de * nrm[0] + dn * nrm[1] + width * cd / 2
    d = depth_c + width * sd / 2
    p = y * cd + d * sd
    q = y * sd - d * cd
    T = _chinnery(_terms_all, x, p, length, width, q, sd, cd, 1 - 2 * nu)
    ex = lambda v: np.asarray(v)[..., None] if np.ndim(v) else v
    u = (T[..., 0:3] * ex(opening) - T[..., 3:6] * ex(strike_slip) - T[..., 6:9] * ex(dip_slip)) / (2 * np.pi)
    ue = u[..., 0] * s[0] + u[..., 1] * nrm[0]
    un = u[..., 0] * s[1] + u[..., 1] * nrm[1]
    return np.stack([ue, un, u[..., 2]], -1)


def dislocation(e, n, ec, nc, depth_c, length, width, strike_deg, dip_deg, strike_slip=0.0, dip_slip=0.0, opening=0.0, nu=0.25):
    """Okada (1985) surface displacement (east, north, up) of a rectangular dislocation with strike-slip (positive =
    left-lateral), dip-slip (positive = reverse / thrust; negative = normal) and opening, centre-based geometry as in
    tensile(). ec, nc, depth_c, length, width and the slips may be arrays that broadcast against e, n (many patches at once)."""
    e, n = np.asarray(e, float), np.asarray(n, float)
    phi, dip = np.deg2rad(strike_deg), np.deg2rad(dip_deg)
    sd, cd = np.sin(dip), np.cos(dip)
    if abs(cd) < 1e-8:
        cd = 0.0
    s = np.array([np.sin(phi), np.cos(phi)])
    nrm = np.array([-np.cos(phi), np.sin(phi)])
    de, dn = e - ec, n - nc
    x = de * s[0] + dn * s[1] + length / 2
    y = de * nrm[0] + dn * nrm[1] + width * cd / 2
    d = depth_c + width * sd / 2
    if np.any(d - width * sd < -1e-9):
        raise ValueError("plane reaches above the surface")
    p = y * cd + d * sd
    q = y * sd - d * cd
    a = 1 - 2 * nu
    u = np.zeros(np.broadcast(x, p, q).shape + (3,))
    if np.any(opening):
        u = u + _chinnery(_terms, x, p, length, width, q, sd, cd, a) * (np.asarray(opening)[..., None] if np.ndim(opening) else opening)
    if np.any(strike_slip) or np.any(dip_slip):
        f = lambda xi, eta, q_, sd_, cd_, a_: np.concatenate(_terms_ss_ds(xi, eta, q_, sd_, cd_, a_), -1)
        sdd = _chinnery(f, x, p, length, width, q, sd, cd, a)
        U1 = np.asarray(strike_slip)[..., None] if np.ndim(strike_slip) else strike_slip
        U2 = np.asarray(dip_slip)[..., None] if np.ndim(dip_slip) else dip_slip
        u = u - sdd[..., :3] * U1 - sdd[..., 3:] * U2
    u = u / (2 * np.pi)
    ue = u[..., 0] * s[0] + u[..., 1] * nrm[0]
    un = u[..., 0] * s[1] + u[..., 1] * nrm[1]
    return np.stack([ue, un, u[..., 2]], -1)


def tensile(e, n, ec, nc, depth_c, length, width, strike_deg, dip_deg, opening, nu=0.25):
    e, n = np.asarray(e, float), np.asarray(n, float)
    phi, dip = np.deg2rad(strike_deg), np.deg2rad(dip_deg)
    sd, cd = np.sin(dip), np.cos(dip)
    if abs(cd) < 1e-8:
        cd = 0.0
    s = np.array([np.sin(phi), np.cos(phi)])     # strike direction (east, north)
    nrm = np.array([-np.cos(phi), np.sin(phi)])  # Okada y: 90 deg anticlockwise from strike
    de, dn = e - ec, n - nc
    x = de * s[0] + dn * s[1] + length / 2
    y = de * nrm[0] + dn * nrm[1] + width * cd / 2
    d = depth_c + width * sd / 2                 # Okada's d: depth of the lower edge
    if d - width * sd < 0:
        raise ValueError("plane reaches above the surface")
    p = y * cd + d * sd
    q = y * sd - d * cd
    lam_over = 1 - 2 * nu                        # mu/(lambda+mu) = 1 - 2 nu
    u = _chinnery(_terms, x, p, length, width, q, sd, cd, lam_over) * opening / (2 * np.pi)
    ue = u[..., 0] * s[0] + u[..., 1] * nrm[0]
    un = u[..., 0] * s[1] + u[..., 1] * nrm[1]
    return np.stack([ue, un, u[..., 2]], -1)


def _checks():
    # 1) small horizontal crack -> point tensile source: uz = 3 dV d^3 / (2 pi R^5), ur = 3 dV r d^2 / (2 pi R^5)
    d, a, U = 3.0, 0.05, 1.0
    r = np.array([0.0, 1.0, 2.0, 4.0, 8.0])
    out = tensile(r, 0 * r, 0, 0, d, a, a, 30, 0, U)
    dV = U * a * a
    R = np.sqrt(r ** 2 + d ** 2)
    uz_pt, ur_pt = 3 * dV * d ** 3 / (2 * np.pi * R ** 5), 3 * dV * r * d ** 2 / (2 * np.pi * R ** 5)
    print("sill vs point crack, uz ratio:", np.round(out[:, 2] / uz_pt, 3), " ur ratio:", np.round(out[1:, 0] / ur_pt[1:], 3))
    # 2) vertical dike striking north: horizontal motion away from the dike, symmetric uplift lobes, no net rotation
    e = np.array([-3.0, -1.0, 1.0, 3.0])
    out = tensile(e, 0 * e, 0, 0, 2.5, 6, 4, 0, 90, 1.0)
    print("dike (strike N), east motion at e = -3,-1,1,3:", np.round(out[:, 0], 4), " up:", np.round(out[:, 2], 4))
    # 3) a large shallow sill uplifts its centre by about its opening
    out = tensile(np.array([0.0]), np.array([0.0]), 0, 0, 0.5, 40, 40, 0, 0, 1.0)
    print("40 x 40 km sill at 0.5 km: centre uplift / opening =", np.round(out[0, 2], 3))
    # 4) a sill rotated by its strike gives the same uplift field rotated
    a1 = tensile(np.array([2.0]), np.array([0.0]), 0, 0, 3, 6, 2, 90, 0, 1.0)[0]
    a2 = tensile(np.array([0.0]), np.array([2.0]), 0, 0, 3, 6, 2, 0, 0, 1.0)[0]
    print("rotation check (uz equal, horizontal swapped):", np.round(a1, 4), np.round(a2, 4))
    # 5) dislocation() reproduces tensile()
    pts = np.array([-3.0, 0.5, 2.0])
    t1 = tensile(pts, pts * 0.7, 0.3, -0.2, 2.5, 6, 3, 35, 80, 1.3)
    t2 = dislocation(pts, pts * 0.7, 0.3, -0.2, 2.5, 6, 3, 35, 80, opening=1.3)
    print("dislocation == tensile:", np.allclose(t1, t2))
    # 6) long faults reaching the surface: the jump across the fault equals the slip
    L, W, y = 400.0, 100.0, 0.01                       # strike north; Okada y (= -east here) positive on the west side
    for nm, kw, comp in (("strike-slip", dict(strike_slip=1.0), 1), ("dip-slip", dict(dip_slip=1.0), 2), ("opening", dict(opening=1.0), 0)):
        up = dislocation(np.array([-y, y]), np.zeros(2), 0, 0, W / 2, L, W, 0, 90, **kw)
        print(f"{nm:11s}: jump across a vertical surface-breaking fault ({'ENU'[comp]}) = {up[0, comp] - up[1, comp]:+.3f} (|jump| should be 1)")
    # 7) vectorized patches == sum of single patches
    ec = np.array([[-1.0, 1.0]]); dep = np.array([[2.0, 3.0]])
    v = dislocation(pts[:, None], pts[:, None] * 0, ec, ec * 0, dep, 2.0, 1.0, 20, 90, strike_slip=np.array([[0.5, -0.2]]), opening=np.array([[1.0, 2.0]])).sum(1)
    s1 = sum(dislocation(pts, pts * 0, ec[0, i], 0, dep[0, i], 2.0, 1.0, 20, 90, strike_slip=[0.5, -0.2][i], opening=[1.0, 2.0][i]) for i in range(2))
    print("vectorized patches == loop:", np.allclose(v, s1))


if __name__ == "__main__":
    _checks()
