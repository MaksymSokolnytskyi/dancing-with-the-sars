"""NISAR GUNW pairs over Mexico City -> vertical subsidence rate, time series,
and a Dufour-style hachure map.

Per pair: subtract the ionosphere screen, drop low-coherence pixels and foreign
unwrapping components, zero at the Sierra de Guadalupe reference, convert phase
to centimetres and from line of sight to vertical. Then stack the chain.

Usage: python scripts/reykjanes/make_figures.py
"""
import re
from datetime import datetime
from pathlib import Path

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from pyproj import Transformer
from scipy.interpolate import RegularGridInterpolator
from scipy.ndimage import gaussian_filter

ROOT = Path(__file__).resolve().parents[2]
WINDOW_DIR = ROOT / "data" / "window"
FRAME_DIR = ROOT / "data" / "gunw"
OUT = ROOT / "output"
C = 299792458.0

BBOX = (-99.32, 19.14, -98.79, 19.60)  # lon0, lat0, lon1, lat1
PLACES = {  # name: (lon, lat)
    "REF Guadalupe": (-99.105, 19.568),
    "Centro": (-99.133, 19.433),
    "Neza/airport": (-99.030, 19.425),
    "Iztapalapa": (-99.050, 19.325),
    "Xochimilco": (-99.105, 19.257),
    "Chalco": (-98.900, 19.257),
}
FRAME_PLACES = {
    "Mexico City": (-99.133, 19.433),
    "Cuernavaca": (-99.234, 18.922),
    "Popocatepetl": (-98.622, 19.023),
}
MIN_COH = 0.3
REF_HALF = 5        # reference window half-size, pixels (80 m posting)
MIN_PAIR_FRAC = 0.8  # a pixel must be valid in this share of the pairs


def dates(name):
    d = re.findall(r"(\d{8})T\d{6}", name)
    return datetime.strptime(d[0], "%Y%m%d"), datetime.strptime(d[2], "%Y%m%d")


def load(path):
    with h5py.File(path, "r") as f:
        return dict(
            name=path.stem, t0=dates(path.name)[0], t1=dates(path.name)[1],
            x=f["x"][()], y=f["y"][()], epsg=int(f.attrs["epsg"]),
            wl=C / float(f.attrs["centerFrequency"]),
            phase=f["unwrappedPhase"][()], coh=f["coherenceMagnitude"][()],
            cc=f["connectedComponents"][()], iono=f["ionospherePhaseScreen"][()],
            inc=f["incidenceAngle"][()], inc_x=f["inc_x"][()], inc_y=f["inc_y"][()],
        )


def load_chain():
    """Pairs sorted by date, keeping only those that continue the chain (end = next start)."""
    pairs = [load(p) for p in sorted(WINDOW_DIR.glob("*.h5"), key=lambda p: dates(p.name))]
    chain = [pairs[0]]
    for p in pairs[1:]:
        if p["t0"] == chain[-1]["t1"]:
            chain.append(p)
    return chain


class Grid:
    def __init__(self, epsg):
        self.tr = Transformer.from_crs(4326, epsg, always_xy=True)
        cx, cy = self.tr.transform([BBOX[0], BBOX[2]] * 2, [BBOX[1]] * 2 + [BBOX[3]] * 2)
        self.x = np.arange(np.ceil(min(cx) / 80) * 80, max(cx), 80.0)
        self.y = np.arange(np.floor(max(cy) / 80) * 80, min(cy), -80.0)
        self.shape = (self.y.size, self.x.size)
        self.extent = [self.x[0] - 40, self.x[-1] + 40, self.y[-1] - 40, self.y[0] + 40]

    def index(self, lon, lat):
        px, py = self.tr.transform(lon, lat)
        return int(np.argmin(np.abs(self.y - py))), int(np.argmin(np.abs(self.x - px)))


def pair_vertical_cm(p, g, iono=True):
    """Vertical displacement of one pair on grid g, cm, NaN where not trusted."""
    j = np.round((g.x - p["x"][0]) / (p["x"][1] - p["x"][0])).astype(int)
    i = np.round((g.y - p["y"][0]) / (p["y"][1] - p["y"][0])).astype(int)
    assert i.min() >= 0 and j.min() >= 0 and i.max() < p["y"].size and j.max() < p["x"].size
    sub = np.ix_(i, j)
    ph, coh, cc = p["phase"][sub].astype(float), p["coh"][sub], p["cc"][sub]
    if iono:
        ph = ph - p["iono"][sub]
    valid = np.isfinite(ph) & (cc > 0) & (coh >= MIN_COH)
    ri, rj = g.index(*PLACES["REF Guadalupe"])
    w = (slice(ri - REF_HALF, ri + REF_HALF + 1), slice(rj - REF_HALF, rj + REF_HALF + 1))
    if not valid[w].any():
        return None
    ref_cc = np.bincount(cc[w][valid[w]]).argmax()  # keep the component holding the reference
    valid &= cc == ref_cc
    ph = np.where(valid, ph, np.nan)
    ph -= np.nanmedian(ph[w])
    f = RegularGridInterpolator((p["inc_y"][::-1], p["inc_x"]), p["inc"][::-1], bounds_error=False, fill_value=None)
    inc = np.deg2rad(f(np.stack(np.meshgrid(g.y, g.x, indexing="ij"), -1)))
    # positive unwrapped phase = range increase (motion away from the radar)
    return -ph * p["wl"] / (4 * np.pi) / np.cos(inc) * 100.0


def stack(pairs, g, iono=True):
    """Returns rate (cm/yr) and cumulative displacement (cm) at every date of the chain."""
    total, tsum = np.zeros(g.shape), np.zeros(g.shape)
    disp, dts = [], []
    for p in pairs:
        d = pair_vertical_cm(p, g, iono)
        dt = (p["t1"] - p["t0"]).days / 365.25
        if d is None:
            print(f"  {p['t0']:%Y-%m-%d} -> {p['t1']:%Y-%m-%d}  reference not valid, pair skipped")
            d = np.full(g.shape, np.nan)
        ok = np.isfinite(d)
        total += np.where(ok, d, 0.0)
        tsum += ok * dt
        disp.append(d)
        dts.append(dt)
    count = sum(np.isfinite(d) for d in disp)
    rate = np.where(count >= MIN_PAIR_FRAC * len(pairs), total / np.where(tsum > 0, tsum, np.nan), np.nan)
    # time series: gaps of a pixel are bridged with its own mean rate
    cum = np.cumsum([np.where(np.isfinite(d), d, rate * dt) for d, dt in zip(disp, dts)], axis=0)
    cum = np.concatenate([np.where(np.isfinite(rate), 0.0, np.nan)[None], cum])
    return rate, cum


def deramp(g, img):
    X, Y = np.meshgrid((g.x - g.x.mean()) / 1e4, (g.y - g.y.mean()) / 1e4)
    m = np.isfinite(img)
    coef, *_ = np.linalg.lstsq(np.c_[X[m], Y[m], np.ones(m.sum())], img[m], rcond=None)
    out = img - (coef[0] * X + coef[1] * Y + coef[2])
    ri, rj = g.index(*PLACES["REF Guadalupe"])  # keep the reference at zero after the fit
    ref = np.nanmedian(out[ri - 6:ri + 7, rj - 6:rj + 7])
    return out - (ref if np.isfinite(ref) else 0.0)


def label(ax, g, places=PLACES):
    for name, (lon, lat) in places.items():
        px, py = g.tr.transform(lon, lat)
        ax.plot(px, py, "+", color="k", ms=9, mew=1.5)
        ax.annotate(name, (px, py), xytext=(5, 0), textcoords="offset points", fontsize=8, va="center",
                    bbox=dict(fc="white", ec="none", alpha=0.6, pad=1.5))


def point_values(g, img):
    out = {}
    for name, (lon, lat) in PLACES.items():
        i, j = g.index(lon, lat)
        out[name] = np.nanmedian(img[max(i - 6, 0):i + 7, max(j - 6, 0):j + 7])
    return out


def fig_rate(pairs, g, fname, iono=True):
    rate, _ = stack(pairs, g, iono)
    flat = deramp(g, rate)
    fig, axes = plt.subplots(1, 2, figsize=(11, 5.2))
    kind = "iono-corrected" if iono else "raw"
    for ax, img, title in zip(axes, (rate, flat), (
            f"NISAR vertical rate, cm/yr ({len(pairs)} desc pairs, {kind})", "planar ramp removed")):
        v = np.nanpercentile(np.abs(img), 99)
        im = ax.imshow(img, cmap="RdBu", vmin=-v, vmax=v, extent=g.extent, interpolation="nearest")
        label(ax, g)
        ax.set_title(title, fontsize=10)
        ax.set_axis_off()
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    fig.tight_layout()
    fig.savefig(OUT / fname, dpi=150, facecolor="white")
    plt.close(fig)
    a, b = point_values(g, rate), point_values(g, flat)
    print(f"{fname}: {pairs[0]['t0']:%Y-%m-%d} -> {pairs[-1]['t1']:%Y-%m-%d}, valid {np.isfinite(rate).mean():.0%}")
    for k in PLACES:
        print(f"  {k:14s} {a[k]:7.1f}  deramped {b[k]:7.1f} cm/yr")
    return rate, flat


def fig_timeseries(pairs, g, cum):
    """Cumulative vertical motion at the labelled places (planar ramp removed per date)."""
    days = [pairs[0]["t0"]] + [p["t1"] for p in pairs]
    flat = np.array([deramp(g, c) if np.isfinite(c).any() else c for c in cum])
    ri, rj = g.index(*PLACES["REF Guadalupe"])
    fig, ax = plt.subplots(figsize=(8, 4.2))
    for name, (lon, lat) in PLACES.items():
        i, j = g.index(lon, lat)
        s = np.nanmedian(flat[:, max(i - 6, 0):i + 7, max(j - 6, 0):j + 7], axis=(1, 2))
        ref = np.nanmedian(flat[:, ri - 6:ri + 7, rj - 6:rj + 7], axis=(1, 2))
        ax.plot(days, s - ref, marker="o", ms=3, lw=1.4, label=name)
    ax.axhline(0, color="0.6", lw=0.8)
    ax.set_ylabel("cumulative vertical motion, cm")
    ax.set_title(f"NISAR time series relative to Sierra de Guadalupe ({len(pairs)} pairs)", fontsize=10)
    ax.legend(fontsize=8, frameon=False)
    ax.grid(alpha=0.3)
    fig.autofmt_xdate()
    fig.tight_layout()
    fig.savefig(OUT / "mexico_city_timeseries.png", dpi=150, facecolor="white")
    plt.close(fig)


def smooth_nan(img, sigma):
    m = np.isfinite(img)
    num = gaussian_filter(np.where(m, img, 0.0), sigma)
    den = gaussian_filter(m.astype(float), sigma)
    return np.where(den > 0.3, num / np.maximum(den, 1e-6), np.nan)


def hachures(g, field, step_cm, spacing_m=240.0, max_len_m=1100.0, min_contour_m=4000.0):
    """Dufour-style hachures: strokes run down the slope of `field` from one
    contour to the next, so steep parts of the bowl get short, dense, heavy strokes."""
    from contourpy import contour_generator
    z = smooth_nan(field, 6)
    gy, gx = np.gradient(z, g.y, g.x)  # d/dy (northing), d/dx (easting), cm per m
    fz = RegularGridInterpolator((g.y[::-1], g.x), np.dstack([gx, gy])[::-1], bounds_error=False, fill_value=np.nan)
    cg = contour_generator(x=g.x, y=g.y, z=np.ma.masked_invalid(z))
    lo, hi = np.nanmin(z), np.nanmax(z)
    segs, widths = [], []
    for k, level in enumerate(np.arange(np.ceil(lo / step_cm) * step_cm, min(hi, -2.5 * step_cm), step_cm)):
        for line in cg.lines(level):
            if len(line) < 2:
                continue
            d = np.r_[0, np.cumsum(np.hypot(*np.diff(line, axis=0).T))]
            if d[-1] < min_contour_m:  # skip specks of noise
                continue
            s = np.arange((k % 2) * spacing_m / 2, d[-1], spacing_m)  # stagger alternate rows
            pts = np.c_[np.interp(s, d, line[:, 0]), np.interp(s, d, line[:, 1])]
            grad = fz(pts[:, ::-1])
            slope = np.hypot(grad[:, 0], grad[:, 1])
            # flat ground (next contour further than max_len_m away) gets no stroke
            ok = np.isfinite(slope) & (slope > step_cm / max_len_m)
            # downslope = towards more subsidence; stroke reaches the next contour
            length = step_cm / slope[ok] * 0.85
            direction = -grad[ok] / slope[ok, None]
            a = pts[ok]
            segs.extend(np.stack([a, a + direction * length[:, None]], axis=1))
            widths.extend(slope[ok])
    return np.array(segs), np.array(widths), z


def fig_hachures(pairs, g, flat_rate, years):
    total = flat_rate * years  # accumulated over the chain, cm
    step = 3.0
    segs, slope, z = hachures(g, total, step)
    paper, ink = "#f4efe1", "#8c1c13"
    fig, ax = plt.subplots(figsize=(9.5, 8.6), facecolor=paper)
    ax.set_facecolor(paper)
    lw = 0.25 + 1.9 * np.clip(slope / np.percentile(slope, 97), 0, 1)
    ax.add_collection(LineCollection(segs, colors=ink, linewidths=lw, capstyle="butt"))
    for name, (lon, lat) in PLACES.items():
        px, py = g.tr.transform(lon, lat)
        ax.plot(px, py, "+", color="#222", ms=8, mew=1.2)
        i, j = g.index(lon, lat)
        val = np.nanmedian(total[max(i - 6, 0):i + 7, max(j - 6, 0):j + 7])
        ax.annotate(f"{name}\n{val:+.0f} cm", (px, py), xytext=(6, 0), textcoords="offset points",
                    fontsize=8.5, va="center", family="serif", style="italic", color="#222",
                    path_effects=[pe.withStroke(linewidth=2.5, foreground=paper)])
    ax.set_xlim(g.extent[0], g.extent[1])
    ax.set_ylim(g.extent[2], g.extent[3])
    ax.set_aspect("equal")
    # 10 km scale bar
    x0, y0 = g.extent[0] + 2500, g.extent[2] + 2500
    ax.plot([x0, x0 + 10000], [y0, y0], color="#222", lw=2, solid_capstyle="butt")
    ax.text(x0 + 5000, y0 + 700, "10 km", ha="center", fontsize=8, family="serif")
    ax.set_xticks([]), ax.set_yticks([])
    for s in ax.spines.values():
        s.set_color("#222"), s.set_linewidth(1.2)
    ax.set_title("Mexico City subsidence bowl, NISAR L-band\n"
                 f"{pairs[0]['t0']:%d %b %Y} to {pairs[-1]['t1']:%d %b %Y}; hachures run downslope, one row per {step:.0f} cm from -{3 * step:.0f} cm,"
                 " heavier where sinking changes fastest", fontsize=10, family="serif")
    fig.tight_layout()
    fig.savefig(OUT / "mexico_city_hachures.png", dpi=170, facecolor=paper)
    plt.close(fig)
    print("hachures:", len(segs), "strokes; deepest", np.nanmin(z).round(1), "cm")


def fig_frame():
    """Full-frame interferogram of one pair, phase re-wrapped to [-pi, pi)."""
    files = sorted(FRAME_DIR.glob("*.h5"), key=lambda p: dates(p.name))
    if not files:
        return
    pair = load(files[0])
    ph = np.where(np.isfinite(pair["phase"]) & (pair["cc"] > 0), pair["phase"], np.nan)
    x, y = pair["x"], pair["y"]
    tr = Transformer.from_crs(4326, pair["epsg"], always_xy=True)
    fig, ax = plt.subplots(figsize=(9, 9))
    ax.imshow(np.angle(np.exp(1j * ph)), cmap="twilight_shifted", vmin=-np.pi, vmax=np.pi,
              extent=[x[0], x[-1], y[-1], y[0]], interpolation="nearest")
    for name, (lon, lat) in FRAME_PLACES.items():
        px, py = tr.transform(lon, lat)
        ax.plot(px, py, "o", mfc="none", mec="yellow", ms=10, mew=1.2,
                path_effects=[pe.withStroke(linewidth=2.5, foreground="k")])
        ax.annotate(name, (px, py), xytext=(9, 0), textcoords="offset points", va="center", fontsize=11,
                    color="white", path_effects=[pe.withStroke(linewidth=2.5, foreground="k")])
    ax.set_title(f"NISAR L-band interferogram, track 113 D frame 079, {pair['t0']:%Y-%m-%d} to {pair['t1']:%Y-%m-%d}",
                 fontsize=10)
    ax.set_axis_off()
    fig.tight_layout()
    fig.savefig(OUT / "frame_wrapped_interferogram.png", dpi=150, facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True)
    chain = load_chain()
    g = Grid(chain[0]["epsg"])
    print(len(chain), "chained pairs", chain[0]["t0"].date(), "->", chain[-1]["t1"].date())

    # reproduction of the reference figure: first 12 pairs, no ionosphere correction
    fig_rate(chain[:12], g, "mexico_city_vertical_rate_12pairs_raw.png", iono=False)
    fig_rate(chain[:12], g, "mexico_city_vertical_rate_12pairs_iono.png", iono=True)

    # full chain
    rate, flat = fig_rate(chain, g, "mexico_city_vertical_rate_full.png", iono=True)
    _, cum = stack(chain, g, iono=True)
    years = (chain[-1]["t1"] - chain[0]["t0"]).days / 365.25
    np.savez_compressed(OUT / "mexico_city_subsidence.npz", x=g.x, y=g.y, epsg=chain[0]["epsg"],
                        rate_cm_yr=rate, rate_deramped_cm_yr=flat, cumulative_cm=cum.astype(np.float32),
                        dates=[f"{d:%Y-%m-%d}" for d in [chain[0]["t0"]] + [p["t1"] for p in chain]])
    fig_timeseries(chain, g, cum)
    fig_hachures(chain, g, flat, years)
    fig_frame()
