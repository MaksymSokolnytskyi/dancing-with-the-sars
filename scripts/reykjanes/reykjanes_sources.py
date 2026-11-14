"""Compare source models for the Svartsengi inflation, fitted jointly to NISAR LOS
rates (159A, 082D) and GNSS velocities (SENG, AUSV, GEVK, VOGC relative to NYLA),
Nov 2025 - Oct 2026, with the same data weights as reykjanes_joint.py.

  mogi        point source
  sill        horizontal rectangular opening crack (Okada 1985): centre, depth,
              length, width, strike free -> also answers "is the reservoir elongated?"
  sill+dike   the sill plus a vertical opening dike along the Sundhnukur crater row
              (centre and strike free, 7 km long, 1-4 km depth): slow dike opening?
  sill+bg     the sill plus the pre-unrest background velocity of SENG relative to NYLA
              (plate spreading / geothermal subsidence) removed from all GNSS stations

Usage: python scripts/reykjanes/reykjanes_sources.py
"""
from datetime import datetime

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from pyproj import Transformer
from scipy.optimize import least_squares

import reykjanes_mogi as M
from okada import tensile
from reykjanes_analysis import OUT, PLACES, Geo, pairs_by_track
from reykjanes_gnss import REF, STATIONS, rate, relative
from reykjanes_joint import N_EFF, PERIOD, SIGMA_GNSS, SIGMA_INSAR, load_gnss

BACKGROUND = (datetime(2015, 7, 1), datetime(2019, 12, 31))  # before the 2020 unrest at Thorbjorn
SUNDHNUKUR = (-22.395, 63.885)  # middle of the crater row
TO_KM3 = 1e-9


class Data:
    def __init__(self):
        self.geo = geo = Geo(pairs_by_track()["159A"][0][2])
        data = M.load(geo)
        P = M.Problem(geo, data)
        obs = P.observations({t: M.rate_map(data[t]["per"]) for t in M.TRACKS})
        self.P, self.obs = P, obs
        # local east/north km around Svartsengi for every pixel and station
        sx, sy = geo.tr.transform(*PLACES["Svartsengi"])
        self.o = (sx / 1000, sy / 1000)
        self.ins = {}
        for t in M.TRACKS:
            m = obs[t]["mask"]
            e, n = self.local(P.X[m], P.Y[m])
            vec = np.stack([c[m] for c in P.vec[t]], -1)
            self.ins[t] = dict(e=e, n=n, vec=vec, y=obs[t]["y"], w=np.sqrt(N_EFF / m.sum()) / SIGMA_INSAR)
        g = load_gnss()
        self.gn = {}
        for s in STATIONS + [REF]:
            px, py = geo.tr.transform(*g[s]["ll"])
            self.gn[s] = dict(en=self.local(px / 1000, py / 1000), v=g[s]["v"])
        t, x, _ = relative("SENG", REF)
        self.bg = rate(t, x, *BACKGROUND)
        print(f"background SENG - {REF} {BACKGROUND[0]:%Y}-{BACKGROUND[1]:%Y}: E {self.bg[0]:+.2f} N {self.bg[1]:+.2f} "
              f"U {self.bg[2]:+.2f} cm/yr")
        self.to_ll = Transformer.from_crs(geo.epsg, 4326, always_xy=True)

    def local(self, X, Y):
        dx, dy = np.asarray(X) - self.o[0], np.asarray(Y) - self.o[1]
        P = self.P
        return dx * P.east[0] + dy * P.east[1], dx * P.north[0] + dy * P.north[1]

    def lonlat(self, e, n):
        P = self.P
        x = self.o[0] + e * P.east[0] + n * P.north[0]
        y = self.o[1] + e * P.east[1] + n * P.north[1]
        return self.to_ll.transform(x * 1000, y * 1000)


# --- sources: each returns ENU (cm/yr) per unit amplitude at points (e, n) --------------------
def mogi(e, n, p):
    e0, n0, d = p
    de, dn = e - e0, n - n0
    R3 = (de ** 2 + dn ** 2 + d ** 2) ** 1.5
    k = 3 / (4 * np.pi) * 1e5  # 1 km^3/yr -> cm/yr
    return np.stack([k * de / R3, k * dn / R3, k * d / R3], -1)


def sill(e, n, p):
    e0, n0, d, L, W, strike = p
    return tensile(e, n, e0, n0, d, L, W, strike, 0.0, 1.0)  # 1 m/yr opening -> m/yr ... scaled below


def dike(e, n, p):
    e0, n0, strike = p
    return tensile(e, n, e0, n0, 2.5, 7.0, 3.0, strike, 90.0, 1.0)


MODELS = {
    # name: list of (source fn, n params, unit -> cm/yr scale, label of amplitude)
    "mogi": [(mogi, 3, 1.0, "dV km^3/yr")],
    "sill": [(sill, 6, 100.0, "opening m/yr")],
    "sill+dike": [(sill, 6, 100.0, "opening m/yr"), (dike, 3, 100.0, "dike opening m/yr")],
    "sill+bg": [(sill, 6, 100.0, "opening m/yr")],
    "mogi+sundhnukur": [(mogi, 3, 1.0, "dV km^3/yr"), (dike, 3, 100.0, "dike opening m/yr")],
    "sill+sundhnukur": [(sill, 6, 100.0, "opening m/yr"), (dike, 3, 100.0, "dike opening m/yr")],
}


def system(D, name, theta):
    srcs = MODELS[name]
    split = np.cumsum([0] + [k for _, k, _, _ in srcs])
    ps = [theta[split[i]:split[i + 1]] for i in range(len(srcs))]
    nlin = len(srcs) + len(M.TRACKS)
    rows, ys = [], []
    for k, t in enumerate(M.TRACKS):
        d = D.ins[t]
        A = np.zeros((d["y"].size, nlin))
        for j, (fn, _, sc, _) in enumerate(srcs):
            A[:, j] = np.sum(fn(d["e"], d["n"], ps[j]) * sc * d["vec"], -1)
        A[:, len(srcs) + k] = 1
        rows.append(A * d["w"])
        ys.append(d["y"] * d["w"])
    ref = D.gn[REF]["en"]
    bg = D.bg if name == "sill+bg" else np.zeros(3)
    for s in STATIONS:
        A = np.zeros((3, nlin))
        en = D.gn[s]["en"]
        for j, (fn, _, sc, _) in enumerate(srcs):
            A[:, j] = (fn(np.array([en[0]]), np.array([en[1]]), ps[j])[0]
                       - fn(np.array([ref[0]]), np.array([ref[1]]), ps[j])[0]) * sc
        rows.append(A / SIGMA_GNSS)
        ys.append((D.gn[s]["v"] - bg) / SIGMA_GNSS)
    return np.vstack(rows), np.concatenate(ys), ps


def solve(D, name, theta):
    A, y, ps = system(D, name, theta)
    m, *_ = np.linalg.lstsq(A, y, rcond=None)
    return m, A @ m - y, ps


def fit(D, name, starts, lo, hi):
    best = None
    for s0 in starts:
        try:
            sol = least_squares(lambda th: solve(D, name, th)[1], s0, bounds=(lo, hi))
        except ValueError:
            continue
        if best is None or sol.cost < best.cost:
            best = sol
    return best.x


def report(D, name, theta):
    m, r, ps = solve(D, name, theta)
    srcs = MODELS[name]
    chi2 = float(np.sum(r ** 2))
    n_ins = sum(1 for _ in M.TRACKS)
    rms = {}
    for k, t in enumerate(M.TRACKS):
        d = D.ins[t]
        pred = sum(np.sum(fn(d["e"], d["n"], ps[j]) * sc * d["vec"], -1) * m[j] for j, (fn, _, sc, _) in enumerate(srcs))
        rms[t] = float(np.sqrt(np.mean((pred + m[len(srcs) + k] - d["y"]) ** 2)))
    ref = D.gn[REF]["en"]
    bg = D.bg if name == "sill+bg" else np.zeros(3)
    gpred = {}
    for s in STATIONS:
        en = D.gn[s]["en"]
        gpred[s] = sum((fn(np.array([en[0]]), np.array([en[1]]), ps[j])[0]
                        - fn(np.array([ref[0]]), np.array([ref[1]]), ps[j])[0]) * sc * m[j]
                       for j, (fn, _, sc, _) in enumerate(srcs))
    gres = np.array([gpred[s] - (D.gn[s]["v"] - bg) for s in STATIONS])
    k_par = len(theta) + len(m)
    print(f"\n== {name}: chi2 {chi2:.1f} with {k_par} parameters (AIC {chi2 + 2 * k_par:.1f}); "
          f"InSAR rms {rms['159A']:.2f}/{rms['082D']:.2f}, GNSS rms {np.sqrt(np.mean(gres ** 2)):.2f} cm/yr")
    for j, (fn, _, sc, lab) in enumerate(srcs):
        p = ps[j]
        lon, lat = D.lonlat(p[0], p[1])
        if fn is mogi:
            print(f"   Mogi {lat:.4f}N {lon:.4f}E depth {p[2]:.2f} km, dV {m[j] * 1e3:.1f}e6 m3/yr ({m[j] * 1e9 / 3.156e7:.2f} m3/s)")
        elif fn is sill:
            vol = m[j] * p[3] * p[4] * 1e6  # m/yr * km^2 -> m^3/yr
            print(f"   sill {lat:.4f}N {lon:.4f}E depth {p[2]:.2f} km, {p[3]:.1f} x {p[4]:.1f} km, strike {p[5] % 180:.0f} deg, "
                  f"opening {m[j] * 100:.1f} cm/yr, volume {vol / 1e6:.1f}e6 m3/yr ({vol / 3.156e7:.2f} m3/s)")
        else:
            vol = m[j] * 7.0 * 3.0 * 1e6
            print(f"   dike {lat:.4f}N {lon:.4f}E strike {p[2] % 180:.0f} deg, 7 x 3 km (1-4 km depth), "
                  f"opening {m[j] * 100:.1f} cm/yr, volume {vol / 1e6:.1f}e6 m3/yr")
    print("   GNSS residual (model - data) E/N/U cm/yr:",
          ", ".join(f"{s} {g[0]:+.1f}/{g[1]:+.1f}/{g[2]:+.1f}" for s, g in zip(STATIONS, gres)))
    return dict(chi2=chi2, k=k_par, rms=rms, gres=gres, m=m, ps=ps, gpred=gpred)


def main():
    D = Data()
    res = {}
    # mogi
    th = fit(D, "mogi", [np.array([-1.5, -0.2, 4.0])], [-8, -8, 0.5], [8, 8, 15])
    res["mogi"] = (th, report(D, "mogi", th))
    # sill: several starts in size and strike
    lo, hi = [-8, -8, 0.8, 0.3, 0.3, -90], [8, 8, 10, 20, 20, 270]
    starts = [np.array([-1.5, -0.2, 4.0, L, W, s]) for L, W in ((3, 3), (8, 3), (12, 5)) for s in (0, 45, 90, 135)]
    th_s = fit(D, "sill", starts, lo, hi)
    res["sill"] = (th_s, report(D, "sill", th_s))
    # sill + dike near the Sundhnukur crater row
    de, dn = D.local(*[c / 1000 for c in D.geo.tr.transform(*SUNDHNUKUR)])
    starts = [np.r_[th_s, de, dn, s] for s in (0, 20, 40)]
    th_d = fit(D, "sill+dike", starts, lo + [de - 6, dn - 5, -30], hi + [de + 3, dn + 5, 90])  # dike free to move
    res["sill+dike"] = (th_d, report(D, "sill+dike", th_d))
    # sill with the pre-unrest background removed from GNSS
    th_b = fit(D, "sill+bg", [th_s], lo, hi)
    res["sill+bg"] = (th_b, report(D, "sill+bg", th_b))
    # the dike pinned to the Sundhnukur crater row (+-1 km, strike N10W - N40E): is there slow opening there?
    dl, dh = [de - 1, dn - 1.5, -10], [de + 1, dn + 1.5, 40]
    th_ms = fit(D, "mogi+sundhnukur", [np.r_[res["mogi"][0], de, dn, s] for s in (0, 20)],
                [-8, -8, 0.5] + dl, [8, 8, 15] + dh)
    res["mogi+sundhnukur"] = (th_ms, report(D, "mogi+sundhnukur", th_ms))
    th_ss = fit(D, "sill+sundhnukur", [np.r_[th_s, de, dn, s] for s in (0, 20)], lo + dl, hi + dh)
    res["sill+sundhnukur"] = (th_ss, report(D, "sill+sundhnukur", th_ss))
    np.savez(OUT / "source_models.npz", **{f"{k}_theta": v[0] for k, v in res.items()},
             **{f"{k}_lin": v[1]["m"] for k, v in res.items()})
    plot(D, res)


def plot(D, res):
    names = list(res)
    fig, axes = plt.subplots(2, len(names), figsize=(4.4 * len(names), 8.4))
    ext_e = (-14, 14)
    for c, name in enumerate(names):
        th, r = res[name]
        srcs = MODELS[name]
        ee, nn = np.meshgrid(np.linspace(*ext_e, 140), np.linspace(-12, 12, 120))
        up = sum(fn(ee, nn, r["ps"][j])[..., 2] * sc * r["m"][j] for j, (fn, _, sc, _) in enumerate(srcs))
        ax = axes[0, c]
        im = ax.contourf(ee, nn, up, levels=np.linspace(-6, 24, 16), cmap="RdBu_r", extend="both")
        for j, (fn, _, _, _) in enumerate(srcs):
            p = r["ps"][j]
            if fn is mogi:
                ax.plot(p[0], p[1], "k*", ms=12)
            else:
                L, W, s = (p[3], p[4], p[5]) if fn is sill else (7.0, 0.0, p[2])
                a = np.deg2rad(s)
                sv, nv = np.array([np.sin(a), np.cos(a)]), np.array([-np.cos(a), np.sin(a)])
                cs = [np.array(p[:2]) + sv * L / 2 * i + nv * W / 2 * k for i, k in ((-1, -1), (1, -1), (1, 1), (-1, 1), (-1, -1))]
                ax.plot(*np.array(cs).T, "k-" if fn is sill else "m-", lw=2)
        for s in STATIONS + [REF]:
            ax.plot(*D.gn[s]["en"], "^", color="gold", mec="k", ms=6)
            g = D.gn[s]["v"] - (D.bg if name == "sill+bg" else 0)
            if s in r["gpred"]:
                ax.quiver(*D.gn[s]["en"], g[0], g[1], color="k", angles="xy", scale_units="xy", scale=2.0, width=0.007)
                ax.quiver(*D.gn[s]["en"], r["gpred"][s][0], r["gpred"][s][1], color="tab:red", angles="xy",
                          scale_units="xy", scale=2.0, width=0.007)
        ax.set_aspect("equal"), ax.set_xlim(*ext_e), ax.set_ylim(-12, 12)
        ax.set_title(f"{name}: model uplift, cm/yr\nGNSS black, model red (1 km = 2 cm/yr)", fontsize=9)
        ax.set_xlabel("km east of Svartsengi"), ax.set_ylabel("km north")
        fig.colorbar(im, ax=ax, fraction=0.04)
        ax = axes[1, c]
        x = np.arange(len(STATIONS))
        for k, (lab, col) in enumerate((("E", "tab:blue"), ("N", "tab:orange"), ("U", "tab:green"))):
            ax.bar(x + (k - 1) * 0.25, r["gres"][:, k], 0.25, color=col, label=lab)
        ax.axhline(0, color="k", lw=0.8)
        ax.set_xticks(x), ax.set_xticklabels(STATIONS)
        ax.set_ylim(-6, 6)
        ax.set_title(f"GNSS residual (model - data), cm/yr\nInSAR rms {r['rms']['159A']:.2f}/{r['rms']['082D']:.2f}; "
                     f"AIC {r['chi2'] + 2 * r['k']:.0f}", fontsize=9)
        ax.legend(fontsize=8, frameon=False)
    fig.suptitle("Svartsengi source models, NISAR + GNSS, Nov 2025 - Oct 2026", fontsize=11)
    fig.tight_layout()
    fig.savefig(OUT / "source_models.png", dpi=115, facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    main()
