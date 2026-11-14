"""Axisymmetric finite-element model of an inflating magma body in a LAYERED elastic crust,
with a first-order topography correction - the "expensive physics" counterpart of Mogi.

Crust (assumption, labelled as such in the outputs):
  Vp near the surface and to ~2.5 km from the station corrections of the VOG model in
  Ruzek et al. 2023 (EGUsphere-2023-2465, Table 1: Vp ~2.0-3.3 km/s at the surface,
  ~3.7-5.6 km/s at 0.5-2.6 km); below 2.5 km a generic Icelandic upper/middle crust
  (6.2 km/s at 4 km, 6.5 at 6 km, 6.8 at 10 km, 7.0 below 15 km).
  Vs = Vp / 1.763 (as in that paper), density from Gardner's rule, moduli = dynamic
  (rho Vs^2), so only their depth CONTRAST matters for the shape of the surface field.
Source: a small sphere (radius 0.3 km) with an isotropic eigenstrain (Eshelby inclusion);
  its volume change is the flux of u through the sphere, i.e. the dV of an equivalent
  Mogi cavity - so dV here and in the Mogi fits mean the same physical quantity.
Topography: the "varying depth" approximation (Williams & Wadge 1998): at every point the
  source depth is counted from that point's elevation (ASTER GDEM v3, 30 m; NASADEM/SRTM do not reach 63.9 N).

FE: P2 triangles on a graded tensor mesh of the (r, z) half plane, 150 x 150 km, rollers on the
axis, the outer side and the bottom.

  green(depth_km) -> r_km, ur, uz (m per m^3 of dV) at the surface
  validate()      -> homogeneous FE vs analytic Mogi
Usage: python scripts/reykjanes/fem_layered.py   (validation + timing)
"""
import time
from functools import lru_cache
from pathlib import Path

import numpy as np
from skfem import Basis, BilinearForm, ElementTriP2, ElementVector, LinearForm, MeshTri, condense, solve
from skfem.helpers import grad

ROOT = Path(__file__).resolve().parents[2]
DEM_COPERNICUS = ROOT / "data" / "dem" / "Copernicus_DSM_COG_10_N63_00_W023_00_DEM.tif"
DEM_ASTER = ROOT / "data" / "dem" / "aster" / "ASTGTMV003_N63W023_dem.tif"     # NASA/METI ASTER GDEM v3 (LP DAAC)
DEM = DEM_ASTER if DEM_ASTER.exists() else DEM_COPERNICUS                       # NASADEM stops at 60 N
R_INC = 0.3      # km
SIZE = 150.0     # km (far boundaries: < 1 % effect within 20 km)
VPVS = 1.763


def vp_profile(z_km):
    """Vp (km/s) at depth z_km (positive down). Assumed profile, see module docstring."""
    zk = np.array([0.0, 2.5, 4.0, 6.0, 10.0, 15.0, 200.0])
    vk = np.array([2.6, 5.5, 6.2, 6.5, 6.8, 7.0, 7.0])
    return np.interp(z_km, zk, vk)


def lame(z_km, layered=True):
    vp = vp_profile(z_km) if layered else np.full_like(np.asarray(z_km, float), 6.2)
    vs = vp / VPVS
    rho = 1.741 * vp ** 0.25 * 1000                  # Gardner, kg/m^3
    mu = rho * (vs * 1000) ** 2
    lam = rho * (vp * 1000) ** 2 - 2 * mu
    return lam / 1e9, mu / 1e9                       # GPa (units cancel in the displacement per dV)


def graded(a, b, h0, ratio=1.08, hmax=4.0):
    """Coordinates from a to b starting with step h0, growing geometrically."""
    out, x, h = [a], a, h0
    while x + h < b:
        x += h
        out.append(x)
        h = min(h * ratio, hmax)
    out.append(b)
    return np.array(out)


def mesh_for(depth):
    """Graded tensor mesh: 0.12 km cells at the axis, the source and the surface, growing 13-20 %
    per cell outwards (checked: < 0.2 % misfit to analytic Mogi in a homogeneous half-space)."""
    r = graded(0.0, SIZE, 0.12, 1.13)
    up = depth - graded(0.0, depth, 0.12, 1.15)[::-1]               # from the surface down to the source
    down = depth + graded(0.0, SIZE - depth, 0.12, 1.15)
    zs = np.unique(np.r_[up, down])
    zs = np.unique(np.r_[zs, graded(0.0, min(depth, 3.0), 0.15, 1.2)])   # refine near the surface too
    return MeshTri.init_tensor(r, -zs[::-1])                       # z <= 0 (up is +z)


def green(depth, layered=True, mesh=None):
    """Surface response to an inflating sphere at `depth` km. Returns r (km), ur, uz in m per m^3."""
    m = mesh if mesh is not None else mesh_for(depth)
    e = ElementVector(ElementTriP2())
    ib = Basis(m, e, intorder=4)

    @BilinearForm
    def a(u, v, w):
        r, z = w.x[0], w.x[1]
        lam, mu = lame(-z, layered)
        gu, gv = grad(u), grad(v)
        eu = (gu[0][0], gu[1][1], u[0] / r, 0.5 * (gu[0][1] + gu[1][0]))
        ev = (gv[0][0], gv[1][1], v[0] / r, 0.5 * (gv[0][1] + gv[1][0]))
        tu, tv = eu[0] + eu[1] + eu[2], ev[0] + ev[1] + ev[2]
        return (lam * tu * tv + 2 * mu * (eu[0] * ev[0] + eu[1] * ev[1] + eu[2] * ev[2] + 2 * eu[3] * ev[3])) * r

    @LinearForm
    def f(v, w):
        r, z = w.x[0], w.x[1]
        lam, mu = lame(-z, layered)
        inside = (r ** 2 + (z + depth) ** 2) < R_INC ** 2
        gv = grad(v)
        tv = gv[0][0] + gv[1][1] + v[0] / r
        return inside * (3 * lam + 2 * mu) / 3 * tv * r            # unit volumetric eigenstrain

    A, b = a.assemble(ib), f.assemble(ib)
    dofs = ib.get_dofs(lambda x: np.isclose(x[0], 0.0)).nodal["u^1"]
    dofs = np.r_[dofs, ib.get_dofs(lambda x: np.isclose(x[0], SIZE)).nodal["u^1"],
                 ib.get_dofs(lambda x: np.isclose(x[1], -SIZE)).nodal["u^2"]]
    u = solve(*condense(A, b, D=np.unique(dofs)))
    # surface displacement
    top = ib.get_dofs(lambda x: np.isclose(x[1], 0.0)).nodal
    xr = ib.doflocs[0][top["u^1"]]
    order = np.argsort(xr)
    rs, ur, uz = xr[order], u[top["u^1"]][order], u[top["u^2"]][order]
    # volume change of the inclusion: flux of u through a sphere of radius R_INC (numerical, 2D quadrature)
    th = np.linspace(0, np.pi, 181)[1:-1]
    RF = 2.5 * R_INC                                                # outside the inclusion the field is divergence-free
    pr, pz = RF * np.sin(th), -depth + RF * np.cos(th)
    P = np.vstack([pr, pz])
    (ur_u, ur_b), (uz_u, uz_b) = ib.split(u)
    urv = ur_b.probes(P) @ ur_u
    uzv = uz_b.probes(P) @ uz_u
    nr, nz = np.sin(th), np.cos(th)
    dS = 2 * np.pi * pr * RF * (th[1] - th[0])
    dV = float(np.sum((urv * nr + uzv * nz) * dS))                  # km^3 per unit eigenstrain (u in km since x in km)
    # u has length units of the mesh (km); per m^3 of dV: divide by dV (km^3 -> m^3) and convert u to m
    scale = 1000.0 / (dV * 1e9)
    return rs, ur * scale, uz * scale


def mogi_surface(r, depth, nu=0.25):
    """Analytic Mogi (homogeneous half-space), m per m^3 of dV, r and depth in km."""
    R3 = ((r ** 2 + depth ** 2) * 1e6) ** 1.5
    k = (1 - nu) / np.pi
    return k * r * 1000 / R3, k * depth * 1000 / R3


def validate():
    for d in (3.0, 5.0):
        t0 = time.time()
        r, ur, uz = green(d, layered=False)
        t1 = time.time() - t0
        vp, vs = 6.2, 6.2 / VPVS
        nu = (vp ** 2 - 2 * vs ** 2) / (2 * (vp ** 2 - vs ** 2))
        mr, mz = mogi_surface(r, d, nu)
        sel = r < 15
        err_z = np.max(np.abs(uz[sel] - mz[sel])) / mz.max()
        err_r = np.max(np.abs(ur[sel] - mr[sel])) / mz.max()
        print(f"homogeneous, depth {d} km: FE solve {t1:.1f} s, {len(r)} surface nodes; max misfit vs Mogi: uz {err_z:.1%}, ur {err_r:.1%} of peak uz "
              f"(nu {nu:.3f})")
        t0 = time.time()
        r2, ur2, uz2 = green(d, layered=True)
        print(f"layered,     depth {d} km: FE solve {time.time() - t0:.1f} s; peak uz {uz2.max() / mz.max():.2f} x homogeneous, "
              f"half-width {np.interp(uz2.max() / 2, uz2[::-1], r2[::-1]):.2f} km vs {np.interp(mz.max() / 2, mz[::-1], r[::-1]):.2f} km")


if __name__ == "__main__":
    validate()
