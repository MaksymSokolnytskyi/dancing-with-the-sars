"""Mogi point-source model of the Svartsengi inflation from NISAR LOS rates.

Data: for the long ascending (159A) and descending (082D) tracks, a LOS rate map
= sum of LOS change over coherent pairs / sum of their time spans (cm/yr), so
winter pairs lost to snow and ionospheric fringes simply drop out.

Model (elastic half-space, Poisson ratio 0.25):
    u_r = 3 dV / (4 pi) * r / R^3,   u_z = 3 dV / (4 pi) * d / R^3,   R^2 = r^2 + d^2
projected on each pixel's target-to-sensor unit vector, plus one offset and a plane
per track (residual long-wavelength atmosphere / referencing). Source x, y, depth are
fitted non-linearly; dV and the ramps are solved linearly inside. Uncertainty from a
bootstrap over pairs (resampling whole pairs keeps the atmospheric noise structure).

Usage: python scripts/reykjanes/reykjanes_mogi.py [n_bootstrap=40] [offset|plane]
"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from pyproj import Transformer
from scipy.optimize import least_squares

from reykjanes_analysis import (MIN_PAIR_COH, OUT, PLACES, Geo, chain, los_cm, pair_coherence,
                                pairs_by_track)

TRACKS = ("159A", "082D")
RADIUS_KM = 15.0  # pixels used around Svartsengi
STEP = 3          # 240 m sampling of the 80 m grid
RNG = np.random.default_rng(0)
RAMP = "offset"  # "offset" or "plane" per track (a plane trades off against a deep source)


def load(geo):
    """Per track: list of (LOS cm map, days) for coherent pairs, and the LOS unit vectors."""
    data = {}
    for t in TRACKS:
        prs = chain(pairs_by_track()[t])
        per, los = [], None
        for t0, t1, p in prs:
            if pair_coherence(p, geo) < MIN_PAIR_COH:
                continue
            d, g = los_cm(p, geo)
            if d is not None:
                per.append((d, (t1 - t0).days))
                los = g
        data[t] = dict(per=per, los=los, span=(prs[0][0], prs[-1][1]))
        print(f"{t}: {len(per)} coherent pairs of {len(prs)}, {prs[0][0]:%Y-%m-%d} -> {prs[-1][1]:%Y-%m-%d}")
    return data


def rate_map(per, idx=None):
    idx = range(len(per)) if idx is None else idx
    num = np.zeros(per[0][0].shape)
    den = np.zeros(per[0][0].shape)
    for k in idx:
        d, dt = per[k]
        ok = np.isfinite(d)
        num += np.where(ok, d, 0)
        den += ok * dt
    frac = den / sum(per[k][1] for k in idx)
    return np.where(frac >= 0.7, num / np.maximum(den, 1) * 365.25, np.nan)


class Problem:
    def __init__(self, geo, data):
        self.geo = geo
        sl = (slice(None, None, STEP), slice(None, None, STEP))
        self.X, self.Y = geo.X[sl] / 1000, geo.Y[sl] / 1000  # km
        self.sel = geo.dist_km[sl] < RADIUS_KM
        self.vec = {}
        for t in TRACKS:
            l = data[t]["los"]
            e, n = l["losUnitVectorX"][sl], l["losUnitVectorY"][sl]
            self.vec[t] = (e, n, np.sqrt(np.clip(1 - e ** 2 - n ** 2, 0, 1)))
        # EPSG:3413 axes are rotated from east/north; angle of the local east direction
        lon0 = -45.0
        lam = np.deg2rad(PLACES["Svartsengi"][0] - lon0)
        self.east = np.array([np.cos(lam), np.sin(lam)])   # map-axis components of the unit east vector
        self.north = np.array([-np.sin(lam), np.cos(lam)])

    def green(self, x0, y0, depth, t, mask):
        """LOS response of a unit (1 km^3/yr) source, in cm/yr, at the masked pixels of track t."""
        dx, dy = self.X[mask] - x0, self.Y[mask] - y0
        r2 = dx ** 2 + dy ** 2
        R3 = (r2 + depth ** 2) ** 1.5
        k = 3 / (4 * np.pi) * 1e5  # km^3 / km^2 -> km, * 1e5 = cm
        ue_m, un_m = k * dx / R3, k * dy / R3  # horizontal motion along the map axes
        ue = ue_m * self.east[0] + un_m * self.east[1]
        un = ue_m * self.north[0] + un_m * self.north[1]
        uz = k * depth / R3
        e, n, u = (c[mask] for c in self.vec[t])
        return ue * e + un * n + uz * u

    def design(self, theta, obs):
        """Columns: source, then offset + plane for each track; rows stacked over tracks."""
        x0, y0, depth = theta
        rows, cols_lin, yv = [], [], []
        for k, t in enumerate(TRACKS):
            m = obs[t]["mask"]
            g = self.green(x0, y0, depth, t, m)
            ramp = np.zeros((m.sum(), 3 * len(TRACKS)))
            ramp[:, 3 * k] = 1
            if RAMP == "plane":
                ramp[:, 3 * k + 1] = self.X[m] - self.X.mean()
                ramp[:, 3 * k + 2] = self.Y[m] - self.Y.mean()
            rows.append(np.c_[g, ramp])
            yv.append(obs[t]["y"])
        return np.vstack(rows), np.concatenate(yv)

    def fit(self, obs, start):
        def resid(theta):
            A, y = self.design(theta, obs)
            m, *_ = np.linalg.lstsq(A, y, rcond=None)
            return A @ m - y
        sol = least_squares(resid, start, bounds=([start[0] - 8, start[1] - 8, 0.5], [start[0] + 8, start[1] + 8, 15]),
                            x_scale=[1, 1, 1])
        A, y = self.design(sol.x, obs)
        m, *_ = np.linalg.lstsq(A, y, rcond=None)
        rms = np.sqrt(np.mean((A @ m - y) ** 2))
        return sol.x, m, rms

    def observations(self, maps):
        obs = {}
        for t in TRACKS:
            v = maps[t][::STEP, ::STEP]
            mask = self.sel & np.isfinite(v)
            obs[t] = dict(mask=mask, y=v[mask])
        return obs


def main(n_boot=40, ramp="offset"):
    global RAMP
    RAMP = ramp
    print("ramp per track:", RAMP)
    OUT.mkdir(parents=True, exist_ok=True)
    tracks = pairs_by_track()
    geo = Geo(tracks["159A"][0][2])
    data = load(geo)
    P = Problem(geo, data)
    maps = {t: rate_map(data[t]["per"]) for t in TRACKS}
    sx, sy = geo.tr.transform(*PLACES["Svartsengi"])
    start = np.array([sx / 1000, sy / 1000, 4.0])
    theta, m, rms = P.fit(P.observations(maps), start)
    to_ll = Transformer.from_crs(geo.epsg, 4326, always_xy=True)
    lon, lat = to_ll.transform(theta[0] * 1000, theta[1] * 1000)
    print(f"best fit: lon {lon:.4f} lat {lat:.4f}, depth {theta[2]:.2f} km, dV {m[0] * 1e3:.2f} x10^6 m^3/yr, rms {rms:.2f} cm/yr")
    print(f"  source is {np.hypot(theta[0] * 1000 - sx, theta[1] * 1000 - sy) / 1000:.1f} km from the Svartsengi marker")

    # bootstrap over pairs
    boot = []
    for b in range(n_boot):
        bm = {t: rate_map(data[t]["per"], RNG.integers(0, len(data[t]["per"]), len(data[t]["per"]))) for t in TRACKS}
        th, mb, _ = P.fit(P.observations(bm), theta)
        boot.append([*th, mb[0] * 1e3])
    boot = np.array(boot)
    lo, hi = np.percentile(boot, [5, 95], axis=0)
    med = np.median(boot, axis=0)
    print(f"bootstrap median: depth {med[2]:.2f} km, dV {med[3]:.2f} x10^6 m^3/yr; "
          f"{(boot[:, 2] > 14.5).sum()} of {len(boot)} runs ran to the 15 km depth bound")
    blon, blat = to_ll.transform(boot[:, 0] * 1000, boot[:, 1] * 1000)
    print(f"bootstrap ({n_boot}, 5-95 %): depth {lo[2]:.2f}-{hi[2]:.2f} km, dV {lo[3]:.2f}-{hi[3]:.2f} x10^6 m^3/yr, "
          f"horizontal scatter {np.std(boot[:, 0]) * 1000:.0f} m E-W axis, {np.std(boot[:, 1]) * 1000:.0f} m N-S axis")

    # maximum uplift rate implied at the surface above the source
    uz0 = 3 / (4 * np.pi) * m[0] / theta[2] ** 2 * 1e5
    print(f"implied uplift rate above the source: {uz0:.1f} cm/yr")
    np.savez(OUT / f"mogi_fit_{RAMP}.npz", theta=theta, lin=m, rms=rms, boot=boot, lon=lon, lat=lat)
    plot(geo, P, maps, theta, m, rms, (lon, lat), boot, (blon, blat), data)


def plot(geo, P, maps, theta, m, rms, ll, boot, bll, data):
    fig, axes = plt.subplots(2, 3, figsize=(14, 8.6))
    ext = [P.X.min(), P.X.max(), P.Y.min(), P.Y.max()]
    for row, t in enumerate(TRACKS):
        full = np.ones(P.X.shape, bool)
        k = TRACKS.index(t)
        model = (P.green(*theta, t, full) * m[0]).reshape(P.X.shape)
        ramp = m[1 + 3 * k] + m[2 + 3 * k] * (P.X - P.X.mean()) + m[3 + 3 * k] * (P.Y - P.Y.mean())
        obs = maps[t][::STEP, ::STEP]
        show = np.where(P.sel, 1, np.nan)
        sp = data[t]["span"]
        for ax, img, title in zip(axes[row], (obs - ramp, model, obs - ramp - model),
                                  (f"{t} data (ramp removed)\n{sp[0]:%d.%m.%y}-{sp[1]:%d.%m.%y}", f"{t} Mogi model",
                                   f"{t} residual")):
            im = ax.imshow(img * show, cmap="RdBu_r", vmin=-16, vmax=16, extent=ext, origin="upper",
                           interpolation="nearest")
            ax.plot(theta[0], theta[1], "k*", ms=10)
            for n in ("Svartsengi", "Grindavik", "Fagradalsfjall"):
                px, py = geo.tr.transform(*PLACES[n])
                ax.plot(px / 1000, py / 1000, "k+", ms=7)
                ax.annotate(n, (px / 1000, py / 1000), xytext=(4, 3), textcoords="offset points", fontsize=7)
            ax.set_title(title + ", LOS cm/yr", fontsize=9)
            ax.set_xticks([]), ax.set_yticks([])
            fig.colorbar(im, ax=ax, fraction=0.045)
    fig.suptitle(f"Svartsengi Mogi source (best fit): {ll[1]:.3f}N {abs(ll[0]):.3f}W, depth {theta[2]:.1f} km, "
                 f"dV {m[0] * 1e3:.1f} x10^6 m^3/yr; rms {rms:.1f} cm/yr; ramp per track: {RAMP}   (* = source)",
                 fontsize=10)
    fig.tight_layout()
    fig.savefig(OUT / f"mogi_fit_{RAMP}.png", dpi=120, facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 40, sys.argv[2] if len(sys.argv) > 2 else "offset")
