"""Reykjanes (Svartsengi - Sundhnukur) deformation from NISAR GUNW pairs.

Per pair: subtract the ionosphere screen, keep coherent pixels of the unwrapping
component that holds the Svartsengi area, convert to line-of-sight centimetres
(positive = towards the satellite), and reference to the far field (median of
coherent pixels more than FAR_KM from Svartsengi, after removing a plane fitted
to them). Pairs of one track are chained into a time series.

Usage: python scripts/reykjanes/reykjanes_analysis.py
"""
import re
from collections import defaultdict
from datetime import datetime
from pathlib import Path

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
from pyproj import Transformer
from scipy.interpolate import RegularGridInterpolator

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data" / "reykjanes" / "gunw"
OUT = ROOT / "output" / "reykjanes"
C = 299792458.0
MIN_COH = 0.3
FAR_KM = 14.0
PLACES = {  # lon, lat
    "Svartsengi": (-22.433, 63.878),
    "Grindavik": (-22.432, 63.842),
    "Sundhnukur": (-22.390, 63.880),
    "Blue Lagoon": (-22.449, 63.880),
    "Fagradalsfjall": (-22.270, 63.900),
    "Keilir": (-22.172, 63.942),
}
NAME = re.compile(r"GUNW_\d{3}_(\d{3})_([AD])_\d{3}_\d{3}_.*?_(\d{8})T\d{6}_\d{8}T\d{6}_(\d{8})T")


def pairs_by_track():
    out = defaultdict(list)
    for p in DATA.glob("*.h5"):
        m = NAME.search(p.name)
        out[m.group(1) + m.group(2)].append((datetime.strptime(m.group(3), "%Y%m%d"),
                                              datetime.strptime(m.group(4), "%Y%m%d"), p))
    for v in out.values():
        v.sort()
    return dict(sorted(out.items()))


def chain(pairs):
    """Path of consecutive pairs (end of one = start of the next) from the earliest start
    to the latest reachable end, using as many (short, more coherent) pairs as possible."""
    pairs = [p for p in pairs if not empty(p[2])]
    if not pairs:
        return []
    start = min(p[0] for p in pairs)
    best = {start: []}  # date -> longest pair path reaching it
    for t0, t1, p in sorted(pairs, key=lambda q: (q[1], q[0])):
        if t0 in best and (t1 not in best or len(best[t0]) + 1 > len(best[t1])):
            best[t1] = best[t0] + [(t0, t1, p)]
    return best[max(best)]


def empty(path):
    """Some granules carry an all-zero unwrapped phase (unwrapping not delivered)."""
    with h5py.File(path, "r") as f:
        return float(np.nanstd(f["unw/unwrappedPhase"][::8, ::8])) == 0.0


class Geo:
    def __init__(self, path):
        with h5py.File(path, "r") as f:
            self.x, self.y = f["unw/x"][()], f["unw/y"][()]
            self.epsg = int(f.attrs["epsg"])
        self.tr = Transformer.from_crs(4326, self.epsg, always_xy=True)
        self.X, self.Y = np.meshgrid(self.x, self.y)
        sx, sy = self.tr.transform(*PLACES["Svartsengi"])
        self.dist_km = np.hypot(self.X - sx, self.Y - sy) / 1000
        self.extent = [self.x[0] - 40, self.x[-1] + 40, self.y[-1] - 40, self.y[0] + 40]

    def ij(self, lon, lat):
        px, py = self.tr.transform(lon, lat)
        return int(np.argmin(np.abs(self.y - py))), int(np.argmin(np.abs(self.x - px)))


def los_cm(path, geo, iono=True):
    with h5py.File(path, "r") as f:
        ph = f["unw/unwrappedPhase"][()].astype(float)
        coh = f["unw/coherenceMagnitude"][()]
        cc = f["unw/connectedComponents"][()]
        if iono:
            ph -= f["unw/ionospherePhaseScreen"][()]
        wl = C / float(f.attrs["centerFrequency"])
        cube = {k: f[f"cube/{k}"][()] for k in ("incidenceAngle", "losUnitVectorX", "losUnitVectorY")}
        cx, cy = f["cube/x"][()], f["cube/y"][()]
    ok = np.isfinite(ph) & (cc > 0) & (coh >= MIN_COH)
    near = ok & (geo.dist_km < 4)
    if near.sum() < 50:
        return None, None
    comp = np.bincount(cc[near]).argmax()
    ok &= cc == comp
    d = np.where(ok, -ph * wl / (4 * np.pi) * 100, np.nan)  # cm, + = towards the satellite
    far = ok & (geo.dist_km > FAR_KM)
    if far.sum() < 200:
        return None, None
    A = np.c_[(geo.X[far] - geo.x.mean()) / 1e4, (geo.Y[far] - geo.y.mean()) / 1e4, np.ones(far.sum())]
    coef, *_ = np.linalg.lstsq(A, d[far], rcond=None)
    d -= coef[0] * (geo.X - geo.x.mean()) / 1e4 + coef[1] * (geo.Y - geo.y.mean()) / 1e4 + coef[2]
    pts = np.stack([geo.Y.ravel(), geo.X.ravel()], -1)
    los = {k: RegularGridInterpolator((cy[::-1], cx), v[::-1], bounds_error=False, fill_value=None)(pts).reshape(d.shape)
           for k, v in cube.items()}
    return d, los


def local(geo, img, name="Svartsengi", r_in=2.0, ring=(8.0, 14.0)):
    """Mean within r_in km of the place minus the median of a ring around it: cancels
    the long-wavelength atmosphere that a single plane leaves behind."""
    px, py = geo.tr.transform(*PLACES[name])
    dk = np.hypot(geo.X - px, geo.Y - py) / 1000
    inner = img[(dk < r_in) & np.isfinite(img)]
    outer = img[(dk > ring[0]) & (dk < ring[1]) & np.isfinite(img)]
    if inner.size < 20 or outer.size < 100:
        return np.nan
    return float(np.median(inner) - np.median(outer))


MIN_PAIR_COH = 0.35  # winter pairs below this are dominated by snow decorrelation and ionospheric fringes


def pair_coherence(path, geo):
    with h5py.File(path, "r") as f:
        c = f["unw/coherenceMagnitude"][()]
    return float(np.nanmedian(c[geo.dist_km < 9]))


def at(geo, img, name, r=4):
    i, j = geo.ij(*PLACES[name])
    return np.nanmedian(img[max(i - r, 0):i + r + 1, max(j - r, 0):j + r + 1])


COMMON = (datetime(2026, 6, 20), datetime(2026, 10, 20))  # all usable tracks observe this summer


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    tracks = pairs_by_track()
    geo = Geo(next(iter(tracks.values()))[0][2])
    res = {}
    for t, prs in tracks.items():
        ch = chain(prs)
        if not ch:
            print(f"{t}: {len(prs)} pairs, none usable (unwrapped phase is empty in every granule)")
            continue
        per, los = [], None
        for t0, t1, p in ch:
            d, g = los_cm(p, geo)
            per.append((t0, t1, d if d is not None else np.full(geo.X.shape, np.nan)))
            los = g if g is not None else los
        dates = [ch[0][0]] + [t1 for _, t1, _ in ch]
        # cumulative maps; a pixel missing in one pair is bridged with 0 there and flagged by count
        steps = np.array([d for *_, d in per])
        cum = np.concatenate([np.zeros((1,) + geo.X.shape), np.cumsum(np.nan_to_num(steps), axis=0)])
        valid = np.concatenate([np.ones((1,) + geo.X.shape, bool), np.cumsum(~np.isfinite(steps), axis=0) == 0])
        coh = np.array([pair_coherence(p, geo) for *_, p in ch])
        step = np.array([local(geo, d) for *_, d in per])
        dt = np.array([(t1 - t0).days for t0, t1, _ in per])
        good = (coh >= MIN_PAIR_COH) & np.isfinite(step)
        rate = np.median(step[good] / dt[good]) if good.any() else np.nan  # cm/day from good pairs
        bridged = np.where(good, step, rate * dt)  # bad pair -> the track's typical rate
        loc = np.r_[0, np.cumsum(bridged)]
        res[t] = dict(pairs=per, dates=dates, cum=cum, valid=valid, local=loc, los=los, n_all=len(prs),
                      good=good, coh=coh, rate_cm_yr=rate * 365.25)
        print(f"{t}: {len(prs)} pairs, chain {len(ch)} {dates[0]:%Y-%m-%d} -> {dates[-1]:%Y-%m-%d}; "
              f"Svartsengi vs ring: total {loc[-1]:+.1f} cm LOS, {(~good).sum()} pairs bridged (coherence < {MIN_PAIR_COH}), "
              f"median rate of good pairs {rate * 365.25:+.1f} cm/yr")
        print("   per pair:", " ".join(f"{t0:%m%d}:{v:+.1f}{'' if g else '*'}" for (t0, _, _), v, g in zip(per, step, good)))
    plot_series(res)
    for t in res:
        plot_pairs(geo, t, res[t])
    vel = velocities(geo, res)
    decompose(geo, res, vel)
    return geo, res


def velocities(geo, res):
    """Per-pixel LOS velocity (cm/yr) over COMMON by a least-squares line through the dates."""
    vel = {}
    for t, r in res.items():
        d = np.array([(x - r["dates"][0]).days / 365.25 for x in r["dates"]])
        sel = np.array([COMMON[0] <= x <= COMMON[1] for x in r["dates"]])
        if sel.sum() < 6:
            continue
        cum = np.where(r["valid"][sel], r["cum"][sel], np.nan)
        tt = d[sel] - d[sel].mean()
        y = cum - np.nanmean(cum, axis=0)
        n = np.isfinite(y).sum(0)
        v = np.nansum(y * tt[:, None, None], 0) / np.nansum(np.where(np.isfinite(y), tt[:, None, None] ** 2, 0), 0)
        vel[t] = np.where(n >= 0.8 * sel.sum(), v, np.nan)
        print(f"velocity {t}: {sel.sum()} dates in common window, Svartsengi-ring {local(geo, vel[t]):+.1f} cm/yr")
    return vel


def decompose(geo, res, vel):
    """Ascending and descending mean velocities -> east and up (north motion assumed small:
    both geometries see north with the same sign, so it cannot be separated)."""
    asc = [t for t in vel if t.endswith("A")]
    dsc = [t for t in vel if t.endswith("D")]
    va, vd = np.nanmean([vel[t] for t in asc], 0), np.nanmean([vel[t] for t in dsc], 0)
    ea = np.nanmean([res[t]["los"]["losUnitVectorX"] for t in asc], 0)
    ed = np.nanmean([res[t]["los"]["losUnitVectorX"] for t in dsc], 0)
    ua = np.cos(np.deg2rad(np.nanmean([res[t]["los"]["incidenceAngle"] for t in asc], 0)))
    ud = np.cos(np.deg2rad(np.nanmean([res[t]["los"]["incidenceAngle"] for t in dsc], 0)))
    det = ea * ud - ed * ua
    east = (va * ud - vd * ua) / det
    up = (ea * vd - ed * va) / det
    np.savez_compressed(OUT / "reykjanes_velocity.npz", x=geo.x, y=geo.y, epsg=geo.epsg, east=east, up=up,
                        **{f"los_{t}": v for t, v in vel.items()})
    print(f"decomposition from asc {asc} and desc {dsc}: Svartsengi up {local(geo, up):+.1f}, "
          f"east {local(geo, east):+.1f} cm/yr (vs ring)")
    fig, axes = plt.subplots(1, 4, figsize=(18, 4.6))
    for ax, img, title, v in zip(axes, (va, vd, up, east),
                                 (f"ascending LOS {'+'.join(asc)}", f"descending LOS {'+'.join(dsc)}",
                                  "vertical (up +)", "east-west (east +)"), (30, 30, 30, 15)):
        im = ax.imshow(img, cmap="RdBu_r", vmin=-v, vmax=v, extent=geo.extent, interpolation="nearest")
        for n in ("Svartsengi", "Grindavik", "Fagradalsfjall", "Keilir"):
            px, py = geo.tr.transform(*PLACES[n])
            ax.plot(px, py, "k+", ms=7)
            ax.annotate(n, (px, py), xytext=(4, 3), textcoords="offset points", fontsize=7)
        ax.set_title(title + ", cm/yr", fontsize=9)
        ax.set_axis_off()
        fig.colorbar(im, ax=ax, fraction=0.04)
    fig.suptitle(f"Reykjanes, NISAR velocity {COMMON[0]:%d %b} - {COMMON[1]:%d %b %Y} (far field > {FAR_KM:.0f} km = 0)",
                 fontsize=10)
    fig.tight_layout()
    fig.savefig(OUT / "velocity_decomposition.png", dpi=120, facecolor="white")
    plt.close(fig)


def plot_series(res):
    fig, ax = plt.subplots(figsize=(10, 4.6))
    for t, r in res.items():
        line, = ax.plot(r["dates"], r["local"], marker="o", ms=3, lw=1.3,
                        label=f"{t} ({r['good'].sum()}/{len(r['pairs'])} pairs used)")
        bad = [i + 1 for i, g in enumerate(r["good"]) if not g]
        ax.plot([r["dates"][i] for i in bad], r["local"][bad], "x", color=line.get_color(), ms=7)
    ax.axhline(0, color="0.6", lw=0.8)
    ax.set_ylabel("Svartsengi minus 8-14 km ring\nLOS cm (+ towards satellite)")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, ncol=3, frameon=False)
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %Y"))
    ax.set_title("NISAR: Svartsengi relative to its surroundings, cumulative by track "
                 f"(x = pair with coherence < {MIN_PAIR_COH}, bridged with the track's median rate)", fontsize=9)
    fig.tight_layout()
    fig.savefig(OUT / "timeseries_local.png", dpi=130, facecolor="white")
    plt.close(fig)


def plot_pairs(geo, t, r):
    per = r["pairs"]
    n = len(per)
    cols = min(7, n)
    rows = int(np.ceil(n / cols))
    fig, axes = plt.subplots(rows, cols, figsize=(2.6 * cols, 2.4 * rows + 0.6), squeeze=False)
    for ax in axes.ravel():
        ax.axis("off")
    sx, sy = geo.tr.transform(*PLACES["Svartsengi"])
    for ax, (t0, t1, d) in zip(axes.ravel(), per):
        ax.imshow(d, cmap="RdBu_r", vmin=-6, vmax=6, extent=geo.extent, interpolation="nearest")
        ax.plot(sx, sy, "k+", ms=6)
        ax.set_title(f"{t0:%d.%m}-{t1:%d.%m.%y}  {local(geo, d):+.1f}", fontsize=7)
    fig.suptitle(f"Track {t}: LOS change per pair, cm (red = towards satellite, ±6); Svartsengi minus ring", fontsize=9)
    fig.tight_layout()
    fig.savefig(OUT / f"pairs_{t}.png", dpi=110, facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    main()
