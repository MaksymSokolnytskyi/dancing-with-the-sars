"""Svartsengi source models fitted jointly to NISAR LOS rates and GNSS velocities:

  mogi       point pressure source                         (x, y, depth)          -> dV
  sill       horizontal rectangular opening crack (Okada)  (x, y, depth, L, W, strike) -> opening
  sill+dike  the sill plus a vertical opening dike along the Sundhnukur crater row
             (x, y, strike, length, bottom depth; top fixed at 0.5 km)            -> two openings

Same data, weights and per-track InSAR offsets as reykjanes_joint.py, so the misfits
are directly comparable. Geometry in local east/north km around the Svartsengi marker.

Usage: python scripts/reykjanes/reykjanes_sill.py
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import least_squares

import reykjanes_mogi as M
from okada import tensile
from reykjanes_analysis import OUT, PLACES, Geo, pairs_by_track
from reykjanes_gnss import REF, STATIONS
from reykjanes_joint import N_EFF, PERIOD, SIGMA_GNSS, SIGMA_INSAR, load_gnss

DIKE_TOP = 0.5  # km


class Data:
    def __init__(self):
        M.RAMP = "offset"
        self.geo = Geo(pairs_by_track()["159A"][0][2])
        data = M.load(self.geo)
        self.P = M.Problem(self.geo, data)
        self.obs = self.P.observations({t: M.rate_map(data[t]["per"]) for t in M.TRACKS})
        sx, sy = self.geo.tr.transform(*PLACES["Svartsengi"])
        self.o = np.array([sx, sy]) / 1000
        E, N = self.enu_xy(self.P.X, self.P.Y)
        self.pix = {}
        for t in M.TRACKS:
            m = self.obs[t]["mask"]
            e, n, u = (c[m] for c in self.P.vec[t])
            self.pix[t] = dict(E=E[m], N=N[m], los=np.c_[e, n, u], y=self.obs[t]["y"],
                               w=np.sqrt(N_EFF / m.sum()) / SIGMA_INSAR)
        g = load_gnss()
        self.gnss = []
        for s in STATIONS + [REF]:
            x, y = self.geo.tr.transform(*g[s]["ll"])
            self.gnss.append((s, *self.enu_xy(x / 1000, y / 1000), g[s]["v"]))

    def enu_xy(self, X, Y):
        dx, dy = X - self.o[0], Y - self.o[1]
        e = dx * self.P.east[0] + dy * self.P.east[1]
        n = dx * self.P.north[0] + dy * self.P.north[1]
        return e, n

    def to_lonlat(self, e, n):
        x = self.o[0] + e * self.P.east[0] + n * self.P.north[0]
        y = self.o[1] + e * self.P.east[1] + n * self.P.north[1]
        from pyproj import Transformer
        return Transformer.from_crs(self.geo.epsg, 4326, always_xy=True).transform(x * 1000, y * 1000)


def mogi_unit(E, N, e0, n0, d):
    de, dn = E - e0, N - n0
    R3 = (de ** 2 + dn ** 2 + d ** 2) ** 1.5
    k = 3 / (4 * np.pi) * 1e5  # 1 km^3/yr -> cm/yr
    return np.stack([k * de / R3, k * dn / R3, k * d / R3], -1)


MODELS = {
    # name: (parameter names, start, lower, upper, unit-response builder -> list of ENU arrays per linear term)
    "mogi": (["e0", "n0", "depth"], [-1.5, -0.3, 3.8], [-8, -8, 0.5], [8, 8, 12],
             lambda p, E, N: [mogi_unit(E, N, *p)]),
    "sill": (["e0", "n0", "depth", "length", "width", "strike"], [-1.5, -0.3, 3.8, 4, 3, 30],
             [-8, -8, 1.0, 0.5, 0.5, -90], [8, 8, 10, 20, 20, 90],
             lambda p, E, N: [tensile(E, N, p[0], p[1], p[2], p[3], p[4], p[5], 0, 1e5)]),  # 1 km opening/yr -> cm/yr
    "sill+dike": (["e0", "n0", "depth", "length", "width", "strike", "de0", "dn0", "dstrike", "dlength", "dbottom"],
                  [-1.5, -0.3, 3.8, 4, 3, 30, 2.0, -1.0, 30, 6, 3.0],
                  [-8, -8, 1.0, 0.5, 0.5, -90, -2, -6, -10, 2, 1.5], [8, 8, 10, 20, 20, 90, 5, 4, 70, 14, 6],
                  lambda p, E, N: [tensile(E, N, p[0], p[1], p[2], p[3], p[4], p[5], 0, 1e5),
                                   tensile(E, N, p[6], p[7], (DIKE_TOP + p[10]) / 2, p[9], p[10] - DIKE_TOP, p[8], 90, 1e5)]),
}


def system(D, model, p):
    names, _, _, _, build = MODELS[model]
    nlin = None
    rows, ys = [], []
    for k, t in enumerate(M.TRACKS):
        px = D.pix[t]
        terms = build(p, px["E"], px["N"])
        nlin = len(terms) + len(M.TRACKS)
        A = np.zeros((px["y"].size, nlin))
        for j, u in enumerate(terms):
            A[:, j] = np.sum(u * px["los"], 1)
        A[:, len(terms) + k] = 1
        rows.append(A * px["w"])
        ys.append(px["y"] * px["w"])
    E = np.array([g[1] for g in D.gnss])
    N = np.array([g[2] for g in D.gnss])
    terms = build(p, E, N)
    for i, (s, _, _, v) in enumerate(D.gnss[:-1]):  # last entry is the reference station
        A = np.zeros((3, nlin))
        for j, u in enumerate(terms):
            A[:, j] = u[i] - u[-1]
        rows.append(A / SIGMA_GNSS)
        ys.append(v / SIGMA_GNSS)
    return np.vstack(rows), np.concatenate(ys)


def solve(D, model, starts):
    names, start, lo, hi, _ = MODELS[model]
    best = None

    def res(p):
        A, y = system(D, model, p)
        m, *_ = np.linalg.lstsq(A, y, rcond=None)
        return A @ m - y
    for s0 in starts:
        s0 = np.clip(s0, np.array(lo) + 1e-3, np.array(hi) - 1e-3)
        sol = least_squares(res, s0, bounds=(lo, hi), x_scale="jac", max_nfev=400)
        if best is None or sol.cost < best.cost:
            best = sol
    A, y = system(D, model, best.x)
    m, *_ = np.linalg.lstsq(A, y, rcond=None)
    return best.x, m, float(np.sum((A @ m - y) ** 2))


def predict_gnss(D, model, p, m):
    _, _, _, _, build = MODELS[model]
    E = np.array([g[1] for g in D.gnss])
    N = np.array([g[2] for g in D.gnss])
    terms = build(p, E, N)
    u = sum(m[j] * terms[j] for j in range(len(terms)))
    return {g[0]: u[i] - u[-1] for i, g in enumerate(D.gnss[:-1])}


def insar_rms(D, model, p, m):
    _, _, _, _, build = MODELS[model]
    out = {}
    for k, t in enumerate(M.TRACKS):
        px = D.pix[t]
        terms = build(p, px["E"], px["N"])
        pred = sum(m[j] * np.sum(terms[j] * px["los"], 1) for j in range(len(terms))) + m[len(terms) + k]
        out[t] = float(np.sqrt(np.mean((pred - px["y"]) ** 2)))
    return out


def main():
    D = Data()
    rng = np.random.default_rng(1)
    fits = {}
    for model in MODELS:
        names, start, lo, hi, _ = MODELS[model]
        starts = [np.array(start)] + [np.array(start) + rng.normal(0, 1, len(start)) * (np.array(hi) - np.array(lo)) * 0.08
                                      for _ in range(7)]
        if model == "sill+dike":  # seed with the best sill
            ps = fits["sill"][0]
            starts = [np.r_[ps, s[6:]] for s in starts]
        p, m, chi2 = solve(D, model, starts)
        fits[model] = (p, m, chi2)
        g = predict_gnss(D, model, p, m)
        ins = insar_rms(D, model, p, m)
        rms_g = np.sqrt(np.mean([(g[s] - v) ** 2 for s, _, _, v in D.gnss[:-1]]))
        k = len(names) + len(m)
        lon, lat = D.to_lonlat(p[0], p[1])
        print(f"\n== {model}: chi2 {chi2:.1f}, parameters {k}, AIC {chi2 + 2 * k:.1f}; "
              f"rms InSAR {ins['159A']:.2f}/{ins['082D']:.2f}, GNSS {rms_g:.2f} cm/yr")
        print("   " + ", ".join(f"{n} {v:.2f}" for n, v in zip(names, p)) + f"  (source at {lat:.4f}N {lon:.4f}E)")
        if model == "mogi":
            print(f"   dV {m[0] * 1e9 / 1e6:.1f}e6 m3/yr")
        else:
            op = m[0] * 1e5  # km/yr -> cm/yr
            dv = m[0] * p[3] * p[4] * 1e9
            print(f"   sill opening {op:.1f} cm/yr over {p[3]:.1f} x {p[4]:.1f} km, strike {p[5]:.0f} deg, "
                  f"dV {dv / 1e6:.1f}e6 m3/yr ({dv / 3.156e7:.2f} m3/s)")
            if model == "sill+dike":
                ddv = m[1] * p[9] * (p[10] - DIKE_TOP) * 1e9
                dl1, dl2 = D.to_lonlat(p[6], p[7])
                print(f"   dike opening {m[1] * 1e5:.1f} cm/yr, centre {dl2:.4f}N {dl1:.4f}E, strike {p[8]:.0f}, "
                      f"length {p[9]:.1f} km, {DIKE_TOP}-{p[10]:.1f} km deep, dV {ddv / 1e6:.1f}e6 m3/yr")
        print("   station  GNSS E/N/U         model E/N/U")
        for s, _, _, v in D.gnss[:-1]:
            print(f"   {s:6s} {v[0]:+5.1f} {v[1]:+5.1f} {v[2]:+5.1f}    {g[s][0]:+5.1f} {g[s][1]:+5.1f} {g[s][2]:+5.1f}")
    np.savez(OUT / "source_models.npz", **{f"{k}_p": v[0] for k, v in fits.items()},
             **{f"{k}_lin": v[1] for k, v in fits.items()})
    plot(D, fits)


def plot(D, fits):
    E, N = D.enu_xy(D.P.X, D.P.Y)
    fig, axes = plt.subplots(len(fits), 3, figsize=(14, 4.2 * len(fits)), squeeze=False)
    for r, (model, (p, m, chi2)) in enumerate(fits.items()):
        build = MODELS[model][4]
        terms = build(p, E.ravel(), N.ravel())
        u = sum(m[j] * terms[j] for j in range(len(terms))).reshape(E.shape + (3,))
        g = predict_gnss(D, model, p, m)
        for c, t in enumerate(M.TRACKS):
            ax = axes[r, c]
            e, n, uu = D.P.vec[t]
            pred = u[..., 0] * e + u[..., 1] * n + u[..., 2] * uu
            data = np.full(E.shape, np.nan)
            data[D.obs[t]["mask"]] = D.obs[t]["y"] - m[len(terms) + c]
            res = data - pred
            im = ax.imshow(res, cmap="RdBu_r", vmin=-12, vmax=12,
                           extent=[E.min(), E.max(), N.min(), N.max()], interpolation="nearest")
            draw_sources(ax, model, p)
            ax.set_title(f"{model}: {t} residual (rms {np.sqrt(np.nanmean(res ** 2)):.1f} cm/yr)", fontsize=9)
            ax.set_xlim(-12, 12), ax.set_ylim(-10, 10)
            fig.colorbar(im, ax=ax, fraction=0.04)
        ax = axes[r, 2]
        for s, e0, n0, v in D.gnss[:-1]:
            ax.quiver(e0, n0, v[0], v[1], color="k", angles="xy", scale_units="xy", scale=2.5, width=0.008)
            ax.quiver(e0, n0, g[s][0], g[s][1], color="tab:red", angles="xy", scale_units="xy", scale=2.5, width=0.008)
            ax.annotate(f"{s} U {v[2]:+.0f}/{g[s][2]:+.0f}", (e0, n0), xytext=(5, -10), textcoords="offset points",
                        fontsize=7)
        draw_sources(ax, model, p)
        ax.set_xlim(-8, 10), ax.set_ylim(-8, 13), ax.set_aspect("equal"), ax.grid(alpha=0.3)
        ax.set_title(f"{model}: GNSS (black) vs model (red); 1 km arrow = 2.5 cm/yr; chi2 {chi2:.0f}", fontsize=9)
    for ax in axes.ravel():
        ax.set_xlabel("east, km from Svartsengi", fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "source_models.png", dpi=110, facecolor="white")
    plt.close(fig)


def outline(e0, n0, L, W, strike):
    s = np.array([np.sin(np.deg2rad(strike)), np.cos(np.deg2rad(strike))])
    q = np.array([-s[1], s[0]])
    c = [np.array([e0, n0]) + a * L / 2 * s + b * W / 2 * q for a, b in ((-1, -1), (1, -1), (1, 1), (-1, 1), (-1, -1))]
    return np.array(c)


def draw_sources(ax, model, p):
    if model == "mogi":
        ax.plot(p[0], p[1], "k*", ms=10)
        return
    o = outline(p[0], p[1], p[3], p[4], p[5])
    ax.plot(o[:, 0], o[:, 1], "k-", lw=1.2)
    if model == "sill+dike":
        s = np.array([np.sin(np.deg2rad(p[8])), np.cos(np.deg2rad(p[8]))])
        a, b = np.array([p[6], p[7]]) - s * p[9] / 2, np.array([p[6], p[7]]) + s * p[9] / 2
        ax.plot([a[0], b[0]], [a[1], b[1]], "m-", lw=2.5)


if __name__ == "__main__":
    main()
