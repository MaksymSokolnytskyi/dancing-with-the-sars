"""Figure for the claim: after every 2026 fountain episode, NISAR L-band backscatter (GCOV, HH,
100 m) darkens downwind (SW) of Halema'uma'u - consistent with fine tephra fall smoothing the
surface - and brightens inside the caldera where rough new lava covers the floor.

All 10 before/after pairs (2 tracks x 5 episodes) are resampled to one 100 m grid around the
vent; each change map is referenced to the median change > 7 km from the vent (rain, soil
moisture, calibration shifts the whole scene).

Panels: context (HH before), mean HH change, how many episodes darken / brighten each cell,
mean change by direction from the vent in a 1.5-6 km ring (trade winds blow from the NE, so
tephra should fall to the SW, azimuth ~225 deg), HH and HV change along the SW and NE axes.

Writes output/kilauea/tephra_evidence.png and tephra_evidence.json.
Usage: python scripts/reykjanes/kilauea_tephra_figure.py
"""
import json
import re
from pathlib import Path

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from pyproj import Transformer
from scipy.interpolate import RegularGridInterpolator

from kilauea_episodes import EVENTS, STAMP, VENT, track

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data" / "kilauea" / "gcov"
OUT = ROOT / "output" / "kilauea"
LOOK = 10            # 10 m -> 100 m
HALF = 8000.0        # m around the vent
STEP = 100.0


def load(path, band):
    with h5py.File(path, "r") as f:
        a = f[f"{band}_db100"][()].astype(np.float32)
        x, y, epsg = f["x"][()], f["y"][()], int(f.attrs["epsg"])
    a[a == -32768] = np.nan
    ny, nx = a.shape[0] // LOOK * LOOK, a.shape[1] // LOOK * LOOK
    p = 10 ** (a[:ny, :nx] / 1000)
    p = np.nanmean(p.reshape(ny // LOOK, LOOK, nx // LOOK, LOOK), axis=(1, 3))
    return 10 * np.log10(p), x[:nx].reshape(-1, LOOK).mean(1), y[:ny].reshape(-1, LOOK).mean(1), epsg


def on_grid(img, x, y, epsg, gx, gy):
    """Sample a product grid (its own UTM zone: 072D is in 4N, 151A in 5N) on the common 4N grid."""
    yy, xx = (y, img) if y[0] < y[-1] else (y[::-1], img[::-1])
    f = RegularGridInterpolator((yy, x), xx, bounds_error=False, fill_value=np.nan)
    GY, GX = np.meshgrid(gy, gx, indexing="ij")
    if epsg != 32604:
        GX, GY = Transformer.from_crs(32604, epsg, always_xy=True).transform(GX, GY)
    return f(np.stack([GY, GX], -1))


def main():
    vx, vy = Transformer.from_crs(4326, 32604, always_xy=True).transform(*VENT)
    gx = vx + np.arange(-HALF, HALF + 1, STEP)
    gy = vy + np.arange(HALF, -HALF - 1, -STEP)
    E, N = np.meshgrid(gx - vx, gy - vy)
    R = np.hypot(E, N) / 1000
    AZ = np.degrees(np.arctan2(E, N)) % 360
    files = sorted(DATA.glob("*.h5"))
    ch = {"HHHH": [], "HVHV": []}
    labels, before0 = [], None
    for ev, lab in EVENTS.items():
        for t in sorted({track(f.name) for f in files}):
            acq = sorted((STAMP.findall(f.name)[0], f) for f in files if track(f.name) == t)
            b = [f for s, f in acq if s < ev][-1:]
            a = [f for s, f in acq if s >= ev][:1]
            if not (b and a):
                continue
            for band in ch:
                A0 = on_grid(*load(b[0], band), gx, gy)
                A1 = on_grid(*load(a[0], band), gx, gy)
                d = A1 - A0
                far = R > 7
                ch[band].append(d - np.nanmedian(d[far]))
                if band == "HHHH" and before0 is None:
                    before0 = A0
            labels.append(f"{lab} {t}")
    H, V = np.array(ch["HHHH"]), np.array(ch["HVHV"])
    n = len(H)
    meanH = np.nanmean(H, 0)
    dark = np.sum(H < -1.5, 0)
    bright = np.sum(H > 1.5, 0)
    # direction from the vent, 1.5-6 km ring
    sect = np.arange(0, 360, 15)
    ring = (R > 1.5) & (R < 6)
    az_mean = np.array([[np.nanmean(h[ring & (AZ >= s) & (AZ < s + 15)]) for s in sect] for h in H])
    # profiles along SW (225) and NE (45) axes
    dist = np.arange(-6, 6.01, 0.2)
    def profile(img, az):
        u = np.deg2rad(az)
        px, py = vx + 1000 * dist * np.sin(u), vy + 1000 * dist * np.cos(u)
        f = RegularGridInterpolator((gy[::-1], gx), img[::-1], bounds_error=False, fill_value=np.nan)
        return f(np.c_[py, px])
    prof_H = np.array([profile(h, 225) for h in H])
    prof_V = np.array([profile(v, 225) for v in V])
    inner = R < 1.0
    sw = ring & (AZ >= 195) & (AZ < 255)
    ne = ring & (AZ >= 15) & (AZ < 75)
    stats = dict(
        pairs=labels,
        caldera_mean_dB=[round(float(np.nanmean(h[inner])), 2) for h in H],
        sw_mean_dB=[round(float(np.nanmean(h[sw])), 2) for h in H],
        ne_mean_dB=[round(float(np.nanmean(h[ne])), 2) for h in H],
        sw_hv_mean_dB=[round(float(np.nanmean(v[sw])), 2) for v in V],
        darkest_sector_deg=[int(sect[np.nanargmin(a)] + 7.5) for a in az_mean])
    (OUT / "tephra_evidence.json").write_text(json.dumps(stats, indent=1, ensure_ascii=False), encoding="utf-8")
    for k in ("caldera_mean_dB", "sw_mean_dB", "ne_mean_dB", "sw_hv_mean_dB", "darkest_sector_deg"):
        print(f"{k:20s} {stats[k]}")

    ext = (-HALF / 1000, HALF / 1000, -HALF / 1000, HALF / 1000)
    fig = plt.figure(figsize=(17, 10))
    gs = fig.add_gridspec(2, 4, hspace=.32, wspace=.38)
    ax = [fig.add_subplot(gs[0, i]) for i in range(4)] + [fig.add_subplot(gs[1, 0:2], projection=None),
                                                          fig.add_subplot(gs[1, 2:4])]
    ax[0].imshow(before0, extent=ext, cmap="gray", vmin=-22, vmax=-2)
    ax[0].set_title("контекст: HH до епізоду 41 (дБ)\nкальдера і кратер видно як темну пляму", fontsize=9)
    im = ax[1].imshow(meanH, extent=ext, cmap="PuOr", vmin=-3, vmax=3)
    ax[1].set_title(f"середня зміна HH\n({n} пар «після − до»)", fontsize=9)
    fig.colorbar(im, ax=ax[1], shrink=.75, label="дБ")
    im = ax[2].imshow(dark, extent=ext, cmap="Oranges", vmin=0, vmax=n)
    ax[2].set_title(f"пар, де темніше на >1.5 дБ\n(з {n})", fontsize=9)
    fig.colorbar(im, ax=ax[2], shrink=.75)
    im = ax[3].imshow(bright, extent=ext, cmap="Purples", vmin=0, vmax=n)
    ax[3].set_title(f"пар, де світліше на >1.5 дБ\n(з {n})", fontsize=9)
    fig.colorbar(im, ax=ax[3], shrink=.75)
    for a in ax[:4]:
        a.plot(0, 0, "k+", ms=10, mew=1.5)
        a.add_patch(plt.Circle((0, 0), 1.5, fill=False, ls=":", lw=.8, color="k"))
        a.add_patch(plt.Circle((0, 0), 6, fill=False, ls=":", lw=.8, color="k"))
        a.set_xlabel("схід, км від Халемаумау", fontsize=8); a.tick_params(labelsize=7)
    ax[0].set_ylabel("північ, км", fontsize=8)
    for a in ax[1:4]:
        a.annotate("", xy=(-4.2, -4.2), xytext=(-1.2, -1.2) if False else (4.5, 4.5),
                   arrowprops=dict(arrowstyle="-|>", color="tab:blue", lw=1.6, alpha=.6))
    ax[1].text(4.6, 4.9, "пасат з ПнСх", color="tab:blue", fontsize=8, ha="right")
    # direction plot
    a = ax[4]
    for row in az_mean:
        a.plot(sect + 7.5, row, color="0.6", lw=.8)
    a.plot(sect + 7.5, np.nanmean(az_mean, 0), color="tab:orange", lw=2.4, label="середнє за всіма парами")
    a.axvspan(195, 255, color="tab:orange", alpha=.12, label="за вітром (ПдЗх)")
    a.axvspan(15, 75, color="tab:blue", alpha=.10, label="проти вітру (ПнСх)")
    a.axhline(0, color="k", lw=.6)
    a.set_xlim(0, 360); a.set_xticks(range(0, 361, 45)); a.set_xticklabels(["Пн", "ПнСх", "Сх", "ПдСх", "Пд", "ПдЗх", "Зх", "ПнЗх", "Пн"])
    a.set_ylabel("зміна HH, дБ (кільце 1.5–6 км)"); a.legend(fontsize=8, frameon=False, loc="lower left")
    a.set_title("зміна яскравості залежно від напрямку від кратера: сірі лінії — окремі пари", fontsize=9)
    # SW-NE profile
    a = ax[5]
    a.plot(dist, np.nanmean(prof_H, 0), color="tab:red", lw=2.2, label="HH (поверхневе відбиття)")
    a.fill_between(dist, *np.nanpercentile(prof_H, [25, 75], 0), color="tab:red", alpha=.15)
    a.plot(dist, np.nanmean(prof_V, 0), color="tab:green", lw=2.2, label="HV (об’ємне розсіяння)")
    a.fill_between(dist, *np.nanpercentile(prof_V, [25, 75], 0), color="tab:green", alpha=.15)
    a.axvspan(-1.0, 1.0, color="tab:purple", alpha=.10, label="кальдера / нова лава (~1 км)")
    a.axhline(0, color="k", lw=.6)
    a.set_xlabel("км від кратера: ← ПнСх (проти вітру)   |   ПдЗх (за вітром) →"); a.set_ylabel("зміна, дБ")
    a.set_xlim(-6, 6)
    a.legend(fontsize=8, frameon=False)
    a.set_title("профіль ПнСх → ПдЗх через кратер (середнє і 25–75 % по 10 парах)", fontsize=9)
    fig.suptitle("Кілауеа, NISAR L-band GCOV (100 м): після епізодів фонтанування 2026 поверхня темніє за вітром і світлішає в кальдері",
                 fontsize=12)
    fig.savefig(OUT / "tephra_evidence.png", dpi=110, facecolor="white", bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
