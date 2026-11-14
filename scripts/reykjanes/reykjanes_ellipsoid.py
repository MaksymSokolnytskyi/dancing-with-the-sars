"""Elongated ellipsoidal reservoir under Svartsengi as a finite compound dislocation
model (CDM, Nikkhoo et al. 2017, GJI 208): three mutually perpendicular rectangular
opening cracks through the centre of an ellipsoid with semi-axes a (horizontal, along
strike), b (horizontal, across strike), c (vertical):

    crack Z  horizontal,          2a x 2b, opening Uz
    crack Y  vertical, along strike, 2a x 2c, opening Uy
    crack X  vertical, across strike, 2b x 2c, opening Ux

Openings are solved linearly with Ux, Uy, Uz >= 0 (an inflating cavity), the
geometry (centre, top depth, a, b, c, strike) non-linearly. Volume change =
4 (Ux b c + Uy a c + Uz a b). Data and weights as in reykjanes_sources.py.

Usage: python scripts/reykjanes/reykjanes_ellipsoid.py
"""
import numpy as np
from scipy.optimize import least_squares, lsq_linear

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import reykjanes_mogi as M
from okada import tensile
from reykjanes_analysis import OUT
from reykjanes_gnss import REF, STATIONS
from reykjanes_joint import SIGMA_GNSS
from reykjanes_sources import Data, mogi

MIN_TOP = 0.3  # km below the surface
RAMP = "offset"  # "offset" or "plane": free InSAR offset, or offset + plane, per track
INSAR_ONLY = False  # True: no GNSS rows, so GNSS stays an independent check
PROFILE = True


def out_name():
    return "ellipsoid_fit" + ("" if RAMP == "offset" else f"_{RAMP}") + ("_insar" if INSAR_ONLY else "") + ".npz"


def cdm(e, n, p):
    """ENU (cm/yr) for 1 m/yr opening of each crack: array (..., 3 components, 3 cracks)."""
    e0, n0, top, a, b, c, strike = p
    d = top + c
    z = tensile(e, n, e0, n0, d, 2 * a, 2 * b, strike, 0.0, 1.0)
    y = tensile(e, n, e0, n0, d, 2 * a, 2 * c, strike, 90.0, 1.0)
    x = tensile(e, n, e0, n0, d, 2 * b, 2 * c, strike + 90.0, 90.0, 1.0)
    return np.stack([x, y, z], -1) * 100.0


def system(D, p, sphere_only=False):
    rows, ys = [], []
    per = 3 if RAMP == "plane" else 1
    nlin = 3 + per * len(M.TRACKS)
    for k, t in enumerate(M.TRACKS):
        d = D.ins[t]
        G = cdm(d["e"], d["n"], p)                       # (npix, 3 comp, 3 cracks)
        A = np.zeros((d["y"].size, nlin))
        A[:, :3] = np.einsum("pcj,pc->pj", G, d["vec"])
        A[:, 3 + per * k] = 1
        if RAMP == "plane":
            A[:, 4 + per * k] = d["e"] / 10
            A[:, 5 + per * k] = d["n"] / 10
        rows.append(A * d["w"])
        ys.append(d["y"] * d["w"])
    ref = D.gn[REF]["en"]
    Gr = cdm(np.array([ref[0]]), np.array([ref[1]]), p)[0]
    for s in ([] if INSAR_ONLY else STATIONS):
        en = D.gn[s]["en"]
        A = np.zeros((3, nlin))
        A[:, :3] = cdm(np.array([en[0]]), np.array([en[1]]), p)[0] - Gr
        rows.append(A / SIGMA_GNSS)
        ys.append(D.gn[s]["v"] / SIGMA_GNSS)
    A, y = np.vstack(rows), np.concatenate(ys)
    if sphere_only:  # one shared opening for the three cracks
        A = np.c_[A[:, :3].sum(1), A[:, 3:]]
    return A, y


def solve(D, p, sphere_only=False):
    A, y = system(D, p, sphere_only)
    k = 1 if sphere_only else 3
    lb = np.r_[np.zeros(k), -np.inf * np.ones(A.shape[1] - k)]
    sol = lsq_linear(A, y, bounds=(lb, np.inf * np.ones(A.shape[1])), method="bvls")
    return sol.x, A @ sol.x - y


def fit(D, starts, lo, hi, transform, sphere_only=False):
    best = None
    for s in starts:
        try:
            r = least_squares(lambda q: solve(D, transform(q), sphere_only)[1], s, bounds=(lo, hi), diff_step=1e-3)
        except ValueError:
            continue
        if best is None or r.cost < best.cost:
            best = r
    return transform(best.x)


def describe(D, label, p, sphere_only=False):
    m, r = solve(D, p, sphere_only)
    e0, n0, top, a, b, c, strike = p
    U = np.repeat(m[0], 3) if sphere_only else m[:3]
    vol = 4 * (U[0] * b * c + U[1] * a * c + U[2] * a * b) * 1e6  # m/yr * km^2 -> m^3/yr
    chi2 = float(np.sum(r ** 2))
    k = 7 + (1 if sphere_only else 3) + (3 if RAMP == "plane" else 1) * len(M.TRACKS)
    if sphere_only:
        k -= 2  # a = b = c
    lon, lat = D.lonlat(e0, n0)
    print(f"\n== {label}: chi2 {chi2:.1f}, {k} parameters, AIC {chi2 + 2 * k:.1f}")
    print(f"   centre {lat:.4f}N {lon:.4f}E, depth {top + c:.2f} km (top {top:.2f}, bottom {top + 2 * c:.2f} km)")
    print(f"   semi-axes a {a:.2f} (along {strike % 180:.0f} deg), b {b:.2f}, c {c:.2f} km; a/b {a / b:.1f}, a/c {a / c:.1f}")
    print(f"   openings Ux {U[0] * 100:.1f}, Uy {U[1] * 100:.1f}, Uz {U[2] * 100:.1f} cm/yr; "
          f"volume {vol / 1e6:.1f}e6 m3/yr ({vol / 3.156e7:.2f} m3/s)")
    ref = D.gn[REF]["en"]
    Gr = (cdm(np.array([ref[0]]), np.array([ref[1]]), p)[0] * U).sum(-1)
    gres = []
    for s in STATIONS:
        en = D.gn[s]["en"]
        pr = (cdm(np.array([en[0]]), np.array([en[1]]), p)[0] * U).sum(-1) - Gr
        gres.append(pr - D.gn[s]["v"])
    gres = np.array(gres)
    print(f"   GNSS rms {np.sqrt(np.mean(gres ** 2)):.2f} cm/yr; residual E/N/U:",
          ", ".join(f"{s} {g[0]:+.1f}/{g[1]:+.1f}/{g[2]:+.1f}" for s, g in zip(STATIONS, gres)))
    if not sphere_only:
        nl = 3 + (3 if RAMP == "plane" else 1) * len(M.TRACKS)
        for k, t in enumerate(M.TRACKS):
            d = D.ins[t]
            G = cdm(d["e"], d["n"], p)
            pred = np.einsum("pcj,pc->pj", G, d["vec"]) @ U
            per = 3 if RAMP == "plane" else 1
            pred = pred + m[3 + per * k] + (m[4 + per * k] * d["e"] / 10 + m[5 + per * k] * d["n"] / 10 if per == 3 else 0)
            ramp = f", plane {m[4 + per * k]:+.1f} E / {m[5 + per * k]:+.1f} N cm/yr per 10 km" if per == 3 else ""
            print(f"   {t}: offset {m[3 + per * k]:+.1f} cm/yr{ramp}, InSAR rms {np.sqrt(np.mean((pred - d['y']) ** 2)):.2f}")
    return dict(m=m, chi2=chi2, k=k, vol=vol, p=p, U=U)


def validate():
    """Small cube with equal openings = isotropic point source -> same pattern as Mogi."""
    r = np.array([0.5, 1.0, 2.0, 4.0, 8.0])
    g = cdm(r, 0 * r, (0, 0, 2.9, 0.05, 0.05, 0.05, 30)).sum(-1)
    mg = mogi(r, 0 * r, (0, 0, 2.95))
    print("validation, small equal-opening CDM / Mogi (should be constant):",
          "up", np.round(g[:, 2] / mg[:, 2] / (g[0, 2] / mg[0, 2]), 3),
          "east", np.round(g[:, 0] / mg[:, 0] / (g[0, 0] / mg[0, 0]), 3))


def main():
    validate()
    D = Data()
    # q = (e0, n0, top, a, b, c, strike); bounds keep the body below MIN_TOP
    lo = [-8, -8, MIN_TOP, 0.2, 0.2, 0.2, -90]
    hi = [8, 8, 8, 8, 8, 4, 270]
    ident = lambda q: tuple(q)
    starts = [np.array([-1.5, -0.2, 2.5, a, b, 0.8, s]) for a, b in ((2.5, 1.0), (1.5, 1.5), (4.0, 0.8))
              for s in (30, 60, 120)]
    res = {}
    p = fit(D, starts, lo, hi, ident)
    res["ellipsoid"] = describe(D, "ellipsoid (CDM, free axes)", p)
    # sphere-like reference: a = b = c, one opening
    sph = lambda q: (q[0], q[1], q[2], q[3], q[3], q[3], 0.0)
    ps = fit(D, [np.array([-1.5, -0.2, 2.5, 0.8])], [-8, -8, MIN_TOP, 0.2], [8, 8, 8, 4], sph, sphere_only=True)
    res["sphere"] = describe(D, "sphere-like CDM (a = b = c)", ps, sphere_only=True)
    if not PROFILE:
        np.savez(OUT / out_name(), p=np.array(p), m=res["ellipsoid"]["m"])
        print("saved", out_name())
        return
    # robustness of the elongation: refit with the strike fixed at a set of values
    print("\nstrike profile (all else free): strike -> chi2")
    prof = []
    for s in range(0, 180, 15):
        fx = lambda q, s=s: (q[0], q[1], q[2], q[3], q[4], q[5], float(s))
        pp = fit(D, [np.r_[np.array(p[:6])]], lo[:6], hi[:6], fx)
        m, r = solve(D, pp)
        prof.append((s, float(np.sum(r ** 2)), pp[3] / pp[4]))
        print(f"   {s:3d} deg: chi2 {prof[-1][1]:7.1f}, a/b {prof[-1][2]:.1f}")
    np.savez(OUT / out_name(), p=np.array(p), m=res["ellipsoid"]["m"], profile=np.array(prof))
    plot(D, res, prof)


def plot(D, res, prof):
    r = res["ellipsoid"]
    fig, axes = plt.subplots(1, 3, figsize=(16, 5.2))
    ee, nn = np.meshgrid(np.linspace(-14, 14, 140), np.linspace(-12, 12, 120))
    G = cdm(ee, nn, r["p"])
    up = (G[..., 2, :] * r["U"]).sum(-1)
    ax = axes[0]
    im = ax.contourf(ee, nn, up, levels=np.linspace(-6, 24, 16), cmap="RdBu_r", extend="both")
    e0, n0, top, a, b, c, s = r["p"]
    t = np.linspace(0, 2 * np.pi, 100)
    sv, nv = np.array([np.sin(np.deg2rad(s)), np.cos(np.deg2rad(s))]), np.array([-np.cos(np.deg2rad(s)), np.sin(np.deg2rad(s))])
    ring = np.array([e0, n0])[:, None] + np.outer(sv, a * np.cos(t)) + np.outer(nv, b * np.sin(t))
    ax.plot(*ring, "k-", lw=2)
    ref = D.gn[REF]["en"]
    Gr = (cdm(np.array([ref[0]]), np.array([ref[1]]), r["p"])[0] * r["U"]).sum(-1)
    for st in STATIONS:
        en = D.gn[st]["en"]
        pr = (cdm(np.array([en[0]]), np.array([en[1]]), r["p"])[0] * r["U"]).sum(-1) - Gr
        v = D.gn[st]["v"]
        ax.plot(*en, "^", color="gold", mec="k", ms=7)
        ax.quiver(*en, v[0], v[1], color="k", angles="xy", scale_units="xy", scale=2.0, width=0.007)
        ax.quiver(*en, pr[0], pr[1], color="tab:red", angles="xy", scale_units="xy", scale=2.0, width=0.007)
        ax.annotate(f"{st} U {v[2]:+.0f}/{pr[2]:+.0f}", en, xytext=(5, -12), textcoords="offset points", fontsize=7)
    ax.set_aspect("equal"), ax.set_xlim(-14, 14), ax.set_ylim(-12, 12)
    ax.set_xlabel("km east of Svartsengi"), ax.set_ylabel("km north")
    ax.set_title("ellipsoid (map outline) and model uplift, cm/yr\nGNSS black, model red (1 km = 2 cm/yr)", fontsize=9)
    fig.colorbar(im, ax=ax, fraction=0.04)
    ax = axes[1]  # vertical section along the long axis
    ax.add_patch(matplotlib.patches.Ellipse((0, -(top + c)), 2 * a, 2 * c, fc="tab:red", alpha=0.35, ec="k"))
    ax.add_patch(matplotlib.patches.Ellipse((0, -(top + c)), 2 * b, 2 * c, fc="none", ec="tab:blue", ls="--"))
    ax.axhline(0, color="k")
    ax.set_xlim(-8, 8), ax.set_ylim(-8, 1), ax.set_aspect("equal")
    ax.set_xlabel("km along (red) / across (blue dashed) the long axis"), ax.set_ylabel("depth, km")
    ax.set_title(f"section: a {a:.1f}, b {b:.1f}, c {c:.1f} km; openings "
                 f"Ux {r['U'][0] * 100:.0f}, Uy {r['U'][1] * 100:.0f}, Uz {r['U'][2] * 100:.0f} cm/yr", fontsize=9)
    ax = axes[2]
    pr = np.array(prof)
    ax.plot(pr[:, 0], pr[:, 1] - pr[:, 1].min(), "o-")
    ax.axhline(4, color="0.6", ls="--", lw=0.8)
    ax.set_xlabel("strike of the long axis, deg"), ax.set_ylabel("delta chi2")
    ax2 = ax.twinx()
    ax2.plot(pr[:, 0], pr[:, 2], "s--", color="tab:red", ms=4)
    ax2.set_ylabel("a / b", color="tab:red")
    ax.set_title("how well is the elongation direction determined?", fontsize=9)
    fig.suptitle(f"Svartsengi ellipsoidal source (CDM), NISAR + GNSS: depth {top + c:.1f} km, "
                 f"{2 * a:.1f} x {2 * b:.1f} x {2 * c:.1f} km, volume {r['vol'] / 1e6:.1f}e6 m^3/yr; "
                 f"AIC {r['chi2'] + 2 * r['k']:.0f} (sphere-like {res['sphere']['chi2'] + 2 * res['sphere']['k']:.0f})",
                 fontsize=10)
    fig.tight_layout()
    fig.savefig(OUT / ("ellipsoid_fit.png" if RAMP == "offset" else f"ellipsoid_fit_{RAMP}.png"), dpi=115, facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    import sys
    RAMP = sys.argv[1] if len(sys.argv) > 1 else "offset"
    INSAR_ONLY = "--insar-only" in sys.argv
    PROFILE = "--no-profile" not in sys.argv
    print("InSAR ramp per track:", RAMP, "| GNSS used in the fit:", not INSAR_ONLY)
    main()
