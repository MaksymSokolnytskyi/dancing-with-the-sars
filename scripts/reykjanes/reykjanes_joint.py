"""Joint Mogi inversion of NISAR LOS rates (159A, 082D) and GNSS velocities at Svartsengi.

InSAR rows: rate maps from coherent pairs (see reykjanes_mogi.py), with a free offset
(and optionally a plane) per track, because InSAR is only relative. GNSS rows: east,
north and up velocities of SENG, AUSV, GEVK, VOGC relative to NYLA over the same
period, modelled as Mogi(station) - Mogi(NYLA), with no free offset: GNSS fixes the
absolute amplitude that InSAR alone cannot.

Weights: an InSAR map carries far fewer independent values than pixels (atmospheric
noise is spatially correlated), so the InSAR block counts as N_EFF values with sigma =
its rms; GNSS values have sigma = SIGMA_GNSS (model error, not the mm formal error).

Usage: python scripts/reykjanes/reykjanes_joint.py [offset|plane]
"""
import sys
from datetime import datetime

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from pyproj import Transformer
from scipy.optimize import least_squares

import reykjanes_mogi as M
from reykjanes_analysis import OUT, PLACES, Geo, pairs_by_track
from reykjanes_gnss import REF, STATIONS, rate, read_tenv3, relative

N_EFF = 200        # effective independent InSAR values per track
SIGMA_INSAR = 4.6  # cm/yr, rms of the InSAR-only fit
SIGMA_GNSS = 1.0   # cm/yr
PERIOD = (datetime(2025, 11, 3), datetime(2026, 10, 17))


class Joint(M.Problem):
    def __init__(self, geo, data, gnss, ramp):
        super().__init__(geo, data)
        self.ramp = ramp
        tr = geo.tr
        self.gnss = []
        rx, ry = tr.transform(*gnss[REF]["ll"])
        self.ref_xy = (rx / 1000, ry / 1000)
        for s in STATIONS:
            px, py = tr.transform(*gnss[s]["ll"])
            self.gnss.append((s, px / 1000, py / 1000, gnss[s]["v"]))

    def enu(self, x0, y0, depth, px, py):
        """Mogi E, N, U for unit volume rate (1 km^3/yr) in cm/yr at map-km points."""
        dx, dy = np.asarray(px) - x0, np.asarray(py) - y0
        R3 = (dx ** 2 + dy ** 2 + depth ** 2) ** 1.5
        k = 3 / (4 * np.pi) * 1e5
        ux, uy, uz = k * dx / R3, k * dy / R3, k * depth / R3
        e = ux * self.east[0] + uy * self.east[1]
        n = ux * self.north[0] + uy * self.north[1]
        return np.stack([e, n, uz], -1)

    def system(self, theta, obs, w_insar, w_gnss):
        x0, y0, depth = theta
        nlin = 1 + 3 * len(M.TRACKS)
        rows, ys = [], []
        for k, t in enumerate(M.TRACKS):
            m = obs[t]["mask"]
            A = np.zeros((m.sum(), nlin))
            A[:, 0] = self.green(x0, y0, depth, t, m)
            A[:, 1 + 3 * k] = 1
            if self.ramp == "plane":
                A[:, 2 + 3 * k] = self.X[m] - self.X.mean()
                A[:, 3 + 3 * k] = self.Y[m] - self.Y.mean()
            rows.append(A * w_insar[t])
            ys.append(obs[t]["y"] * w_insar[t])
        ref = self.enu(x0, y0, depth, *self.ref_xy)
        for s, px, py, v in self.gnss:
            g = self.enu(x0, y0, depth, px, py) - ref
            A = np.zeros((3, nlin))
            A[:, 0] = g
            rows.append(A * w_gnss)
            ys.append(v * w_gnss)
        return np.vstack(rows), np.concatenate(ys)

    def fit(self, obs, start, w_insar, w_gnss):
        def res(theta):
            A, y = self.system(theta, obs, w_insar, w_gnss)
            m, *_ = np.linalg.lstsq(A, y, rcond=None)
            return A @ m - y
        sol = least_squares(res, start, bounds=([start[0] - 8, start[1] - 8, 0.5], [start[0] + 8, start[1] + 8, 15]))
        A, y = self.system(sol.x, obs, w_insar, w_gnss)
        m, *_ = np.linalg.lstsq(A, y, rcond=None)
        return sol.x, m

    def misfits(self, theta, m, obs):
        out = {}
        for k, t in enumerate(M.TRACKS):
            msk = obs[t]["mask"]
            pred = self.green(*theta, t, msk) * m[0] + m[1 + 3 * k]
            if self.ramp == "plane":
                pred += m[2 + 3 * k] * (self.X[msk] - self.X.mean()) + m[3 + 3 * k] * (self.Y[msk] - self.Y.mean())
            out[t] = float(np.sqrt(np.mean((pred - obs[t]["y"]) ** 2)))
        ref = self.enu(*theta, *self.ref_xy)
        pred = {s: (self.enu(*theta, px, py) - ref) * m[0] for s, px, py, _ in self.gnss}
        out["gnss"] = float(np.sqrt(np.mean([(pred[s] - v) ** 2 for s, _, _, v in self.gnss])))
        return out, pred


def load_gnss():
    g = {}
    for s in STATIONS + [REF]:
        t, x, ll = (relative(s, REF) if s != REF else (*read_tenv3(REF)[:2], read_tenv3(REF)[2]))
        g[s] = dict(ll=ll, v=rate(t, x, *PERIOD) if s != REF else np.zeros(3))
    return g


def main(ramp="offset"):
    M.RAMP = ramp
    geo = Geo(pairs_by_track()["159A"][0][2])
    data = M.load(geo)
    gnss = load_gnss()
    P = Joint(geo, data, gnss, ramp)
    obs = P.observations({t: M.rate_map(data[t]["per"]) for t in M.TRACKS})
    w_insar = {t: np.sqrt(N_EFF / obs[t]["y"].size) / SIGMA_INSAR for t in M.TRACKS}
    sx, sy = geo.tr.transform(*PLACES["Svartsengi"])
    start = np.array([sx / 1000, sy / 1000, 4.0])
    to_ll = Transformer.from_crs(geo.epsg, 4326, always_xy=True)

    def describe(label, theta, m):
        mis, pred = P.misfits(theta, m, obs)
        lon, lat = to_ll.transform(theta[0] * 1000, theta[1] * 1000)
        dv = m[0] * 1e9  # m^3/yr
        print(f"{label:28s} {lat:.4f}N {lon:.4f}E depth {theta[2]:5.2f} km  dV {dv / 1e6:5.1f}e6 m3/yr "
              f"({dv / 3.156e7:.2f} m3/s)  rms InSAR {mis['159A']:.2f}/{mis['082D']:.2f}, GNSS {mis['gnss']:.2f} cm/yr")
        return mis, pred, (lon, lat)

    print(f"ramp per track: {ramp}; InSAR counted as {N_EFF} values per track, sigma {SIGMA_INSAR}; GNSS sigma {SIGMA_GNSS}")
    th0, m0 = P.fit(obs, start, w_insar, 0.0)
    describe("InSAR only", th0, m0)
    results = {}
    for f in (0.3, 1.0, 3.0):
        th, m = P.fit(obs, start, w_insar, f / SIGMA_GNSS)
        results[f] = (th, m, *describe(f"joint, GNSS weight x{f}", th, m))
    thg, mg = P.fit(obs, start, {t: 1e-6 for t in M.TRACKS}, 1 / SIGMA_GNSS)
    describe("GNSS only (12 values)", thg, mg)

    th, m, mis, pred, ll = results[1.0]
    print("\nstation  GNSS E/N/U          joint model E/N/U    (cm/yr, relative to NYLA)")
    for s, _, _, v in P.gnss:
        print(f"{s:6s} {v[0]:+5.1f} {v[1]:+5.1f} {v[2]:+5.1f}     {pred[s][0]:+5.1f} {pred[s][1]:+5.1f} {pred[s][2]:+5.1f}")

    print("\ndepth profile (joint, GNSS weight x1): depth -> total misfit, dV")
    prof = []
    for d in (2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8):
        def res(xy):
            A, y = P.system([xy[0], xy[1], d], obs, w_insar, 1 / SIGMA_GNSS)
            mm, *_ = np.linalg.lstsq(A, y, rcond=None)
            return A @ mm - y
        s = least_squares(res, th[:2])
        A, y = P.system([*s.x, d], obs, w_insar, 1 / SIGMA_GNSS)
        mm, *_ = np.linalg.lstsq(A, y, rcond=None)
        chi2 = float(np.sum((A @ mm - y) ** 2))
        prof.append((d, chi2, mm[0] * 1e3))
    cmin = min(p[1] for p in prof)
    for d, c, dv in prof:
        print(f"  {d:4.1f} km: chi2 {c:7.1f} (+{c - cmin:5.1f}), dV {dv:5.1f}e6 m3/yr")
    ok = [d for d, c, _ in prof if c - cmin <= 4.0]  # ~2 sigma for 1-3 effective parameters
    print(f"depths within delta-chi2 <= 4 of the best: {min(ok)} - {max(ok)} km")
    np.savez(OUT / f"joint_fit_{ramp}.npz", theta=th, lin=m, lon=ll[0], lat=ll[1], profile=np.array(prof))
    plot(P, geo, obs, th, m, pred, ll, prof, ramp)


def plot(P, geo, obs, th, m, pred, ll, prof, ramp):
    fig = plt.figure(figsize=(15, 8.5))
    gs = fig.add_gridspec(2, 4)
    ext = [P.X.min(), P.X.max(), P.Y.min(), P.Y.max()]
    full = np.ones(P.X.shape, bool)
    for r, t in enumerate(M.TRACKS):
        k = M.TRACKS.index(t)
        mod = (P.green(*th, t, full) * m[0]).reshape(P.X.shape) + m[1 + 3 * k]
        if ramp == "plane":
            mod += m[2 + 3 * k] * (P.X - P.X.mean()) + m[3 + 3 * k] * (P.Y - P.Y.mean())
        data = np.full(P.X.shape, np.nan)
        data[obs[t]["mask"]] = obs[t]["y"]
        off = m[1 + 3 * k]  # the free InSAR offset is removed from data and model for display
        for c, (img, title) in enumerate(((data - off, f"{t} data (offset removed)"), (mod - off, f"{t} joint model"),
                                          (data - mod, f"{t} residual"))):
            ax = fig.add_subplot(gs[r, c])
            im = ax.imshow(img, cmap="RdBu_r", vmin=-16, vmax=16, extent=ext, interpolation="nearest")
            ax.plot(th[0], th[1], "k*", ms=10)
            for s, px, py, v in P.gnss:
                ax.plot(px, py, "^", color="gold", mec="k", ms=6)
            ax.set_title(title + ", LOS cm/yr", fontsize=9)
            ax.set_xticks([]), ax.set_yticks([])
            fig.colorbar(im, ax=ax, fraction=0.045)
    ax = fig.add_subplot(gs[0, 3])
    xs = [px for _, px, _, _ in P.gnss] + [th[0]]
    ys = [py for _, _, py, _ in P.gnss] + [th[1]]
    for s, px, py, v in P.gnss:
        # arrows in map axes: rotate east/north back onto the EPSG:3413 x/y directions
        gx, gy = v[0] * P.east + v[1] * P.north
        mx, my = pred[s][0] * P.east + pred[s][1] * P.north
        ax.quiver(px, py, gx, gy, color="k", angles="xy", scale_units="xy", scale=2.5, width=0.008)
        ax.quiver(px, py, mx, my, color="tab:red", angles="xy", scale_units="xy", scale=2.5, width=0.008)
        ax.annotate(f"{s}  U {v[2]:+.0f}/{pred[s][2]:+.0f}", (px, py), xytext=(6, -10), textcoords="offset points", fontsize=7)
    ax.plot(th[0], th[1], "k*", ms=10)
    ax.set_xlim(min(xs) - 6, max(xs) + 6), ax.set_ylim(min(ys) - 4, max(ys) + 4)
    ax.set_aspect("equal")
    ax.set_title("GNSS (black) vs model (red) horizontal;\n1 km of arrow = 2.5 cm/yr; U = up GNSS/model, cm/yr", fontsize=9)
    ax.set_xticks([]), ax.set_yticks([])
    ax = fig.add_subplot(gs[1, 3])
    pr = np.array(prof)
    ax.plot(pr[:, 0], pr[:, 1] - pr[:, 1].min(), "o-")
    ax.axhline(4, color="0.6", ls="--", lw=0.8)
    ax.set_xlabel("source depth, km"), ax.set_ylabel("delta chi2")
    ax2 = ax.twinx()
    ax2.plot(pr[:, 0], pr[:, 2], "s--", color="tab:red", ms=4)
    ax2.set_ylabel("dV, 10^6 m^3/yr", color="tab:red")
    ax.set_title("depth profile", fontsize=9)
    fig.suptitle(f"Joint NISAR + GNSS Mogi source: {ll[1]:.3f}N {abs(ll[0]):.3f}W, depth {th[2]:.1f} km, "
                 f"dV {m[0] * 1e3:.1f}x10^6 m^3/yr ({m[0] * 1e9 / 3.156e7:.2f} m^3/s); "
                 f"{PERIOD[0]:%b %Y}-{PERIOD[1]:%b %Y}", fontsize=11)
    fig.tight_layout()
    fig.savefig(OUT / f"joint_fit_{ramp}.png", dpi=115, facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "offset")
