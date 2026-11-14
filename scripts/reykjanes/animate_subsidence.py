"""Animation: Dufour-style hachures of the Mexico City subsidence bowl growing
date by date. Needs output/mexico_city_subsidence.npz from make_figures.py.

Usage: python scripts/reykjanes/animate_subsidence.py
"""
from datetime import datetime

import matplotlib
matplotlib.use("Agg")
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.collections import LineCollection
from scipy.ndimage import gaussian_filter1d

from make_figures import OUT, PLACES, Grid, deramp, hachures

STEP_CM = 3.0
PAPER, INK = "#f4efe1", "#8c1c13"


def load():
    d = np.load(OUT / "mexico_city_subsidence.npz")
    g = Grid(int(d["epsg"]))
    days = [datetime.strptime(s, "%Y-%m-%d") for s in d["dates"]]
    cum = np.array([deramp(g, c.astype(float)) for c in d["cumulative_cm"]])
    # Single dates carry centimetres of atmosphere. Keep each pixel's linear trend and
    # smooth only the departure from it over neighbouring dates, so the ends stay unbiased.
    t = np.array([(x - days[0]).days / 365.25 for x in days])
    rate = cum[-1] / t[-1]
    trend = rate[None] * t[:, None, None]
    cum = trend + gaussian_filter1d(cum - trend, 1.5, axis=0, mode="nearest")
    return g, days, cum


def main():
    g, days, cum = load()
    frames = [hachures(g, c, STEP_CM) if i else (np.empty((0, 2, 2)), np.empty(0), c) for i, c in enumerate(cum)]
    wmax = np.percentile(frames[-1][1], 97)  # one stroke-weight scale for every frame

    fig, ax = plt.subplots(figsize=(7.6, 7.1), facecolor=PAPER)
    ax.set_facecolor(PAPER)
    lc = LineCollection([], colors=INK, capstyle="butt")
    ax.add_collection(lc)
    labels = {}
    for name, (lon, lat) in PLACES.items():
        px, py = g.tr.transform(lon, lat)
        ax.plot(px, py, "+", color="#222", ms=7, mew=1.1)
        labels[name] = ax.annotate("", (px, py), xytext=(6, 0), textcoords="offset points", fontsize=8, va="center",
                                   family="serif", style="italic", color="#222",
                                   path_effects=[pe.withStroke(linewidth=2.5, foreground=PAPER)])
    ax.set_xlim(g.extent[0], g.extent[1])
    ax.set_ylim(g.extent[2], g.extent[3])
    ax.set_aspect("equal")
    x0, y0 = g.extent[0] + 2500, g.extent[2] + 2500
    ax.plot([x0, x0 + 10000], [y0, y0], color="#222", lw=2, solid_capstyle="butt")
    ax.text(x0 + 5000, y0 + 700, "10 km", ha="center", fontsize=8, family="serif")
    ax.set_xticks([]), ax.set_yticks([])
    title = ax.set_title("x\nx", fontsize=10, family="serif")  # two lines so tight_layout leaves room
    fig.tight_layout()

    def draw(k):
        k = min(k, len(days) - 1)  # extra frames hold the last date
        segs, slope, _ = frames[k]
        lc.set_segments(segs)
        lc.set_linewidths(0.25 + 1.9 * np.clip(slope / wmax, 0, 1))
        for name, (lon, lat) in PLACES.items():
            i, j = g.index(lon, lat)
            v = np.nanmedian(cum[k, max(i - 6, 0):i + 7, max(j - 6, 0):j + 7])
            labels[name].set_text(f"{name}\n{round(v) + 0:+d} cm")
        title.set_text("Mexico City subsidence, NISAR L-band\n"
                       f"{days[0]:%d %b %Y} to {days[k]:%d %b %Y}   (one hachure row per {STEP_CM:.0f} cm from -{3 * STEP_CM:.0f} cm)")
        return [lc, title, *labels.values()]

    anim = FuncAnimation(fig, draw, frames=len(days) + 8, blit=False)
    out = OUT / "mexico_city_hachures.gif"
    anim.save(out, writer=PillowWriter(fps=4), dpi=110, savefig_kwargs={"facecolor": PAPER})
    print(out, f"{out.stat().st_size / 1e6:.1f} MB,", len(days), "dates")


if __name__ == "__main__":
    main()
