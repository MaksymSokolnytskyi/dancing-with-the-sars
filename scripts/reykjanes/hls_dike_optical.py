"""Optical check of the 10 Nov 2023 dike: where did the dike reach the surface?

NASA HLS (Harmonized Landsat Sentinel-2, HLSS30 v2.0, 30 m, LP DAAC), tile T27VVL. Each post-eruption scene
is compared with a pre-dike scene of the SAME season (moss reflectance and sun angle change a lot over the year):
  2023-03-08 -> 2024-03-02  after the eruptions of 18 Dec 2023, 14 Jan and 8 Feb 2024 on the Sundhnukur crater row
  2023-04-14 -> 2024-04-26  during the 16 Mar - 9 May 2024 eruption (active vents glow in SWIR)

  - SWIR false colour (B12, B11, B4): fresh basalt is black, moss-covered old lava is green-grey,
    active lava saturates B12 (thermal emission at 2.2 um) -> bright orange/white
  - new lava: NIR (B8A) reflectance dropped by > 0.06 and is now < 0.10 (fresh basalt is dark
    in all bands, moss is bright in NIR); clouds / shadow / snow removed with Fmask
  - active vents (2024-04-26): B12 > 0.25 and B12 > 1.2 x B11 (thermal excess), clustered
  - comparison with the dike from Sentinel-1 + NPE: distance from the vents and from the eruptive
    fissures (new-lava pixels at the upstream end) to each dike trace.
The eruptions of 2024 came from later dikes in the same corridor, so this checks the corridor and the
strike rather than the exact November 2023 dike.

Writes output/reykjanes/ml/hls_dike_optical.png and .json.
Usage: python scripts/reykjanes/hls_dike_optical.py
"""
import glob
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import rasterio
from pyproj import Transformer
from rasterio.windows import from_bounds
from scipy import ndimage

from reykjanes_analysis import PLACES

ROOT = Path(__file__).resolve().parents[2]
HLS = ROOT / "data" / "hls"
ML = ROOT / "output" / "reykjanes" / "ml"
SCENES = {"2023-03-08": "HLS.S30.T27VVL.2023067T131251", "2023-04-14": "HLS.S30.T27VVL.2023104T130301",
          "2024-03-02": "HLS.S30.T27VVL.2024062T131301", "2024-04-26": "HLS.S30.T27VVL.2024117T131259"}
PAIRS = {"2024-03-02": "2023-03-08", "2024-04-26": "2023-04-14"}       # post -> same-season pre-dike reference
BOX = (-10, 10, -9, 9)        # km from Svartsengi (east min, max, north min, max)


def load(scene, band, sx, sy):
    f = HLS / scene / f"{scene}.v2.0.{band}.tif"
    with rasterio.open(f) as r:
        w = from_bounds(sx + BOX[0] * 1000, sy + BOX[2] * 1000, sx + BOX[1] * 1000, sy + BOX[3] * 1000, r.transform).round_offsets().round_lengths()
        v = r.read(1, window=w).astype(float)
        if band != "Fmask":
            v[v == r.nodata] = np.nan
            v *= 1e-4
        return v


def bad(fm):
    fm = fm.astype(int)
    return ((fm >> 1) & 1) | ((fm >> 2) & 1) | ((fm >> 3) & 1) | ((fm >> 4) & 1)     # cloud, adjacent, shadow, snow


def rgb(s):
    out = np.dstack([s["B12"], s["B11"], s["B04"]])
    lo, hi = np.nanpercentile(out, 2, (0, 1)), np.nanpercentile(out, 98, (0, 1))
    return np.clip((out - lo) / (hi - lo), 0, 1)


def dist_to_segment(px, py, x0, y0, x1, y1):
    vx, vy = x1 - x0, y1 - y0
    t = np.clip(((px - x0) * vx + (py - y0) * vy) / (vx * vx + vy * vy), 0, 1)
    return np.hypot(px - (x0 + t * vx), py - (y0 + t * vy))


def dike_traces():
    tr = {}
    def line(p):
        u = np.array([np.sin(np.deg2rad(p[2])), np.cos(np.deg2rad(p[2]))])
        return (p[0] - u[0] * p[3] / 2, p[1] - u[1] * p[3] / 2, p[0] + u[0] * p[3] / 2, p[1] + u[1] * p[3] / 2)
    if (ML / "s1_grindavik_dike.npz").exists():
        z = np.load(ML / "s1_grindavik_dike.npz")
        tr["NPE, однакове розкриття"] = line(np.median(z["ps"], 0))
        tr["MCMC, однакове розкриття"] = line(np.median(z["pmc"], 0))
    if (ML / "s1_dike_variable.npz").exists():
        z = np.load(ML / "s1_dike_variable.npz")
        tr["NPE, змінне розкриття"] = line(np.median(z["ps"], 0))
        tr["MCMC, змінне розкриття"] = line(np.median(z["pmc"], 0))
    return tr


def main():
    to_utm = Transformer.from_crs(4326, 32627, always_xy=True)
    sx, sy = to_utm.transform(*PLACES["Svartsengi"])
    D = {}
    for d, sc in SCENES.items():
        D[d] = {b: load(sc, b, sx, sy) for b in ("B04", "B8A", "B11", "B12", "Fmask")}
        D[d]["bad"] = bad(D[d]["Fmask"]).astype(bool)
        print(d, f"masked (cloud/shadow/snow) {D[d]['bad'].mean():.1%}", flush=True)
    ny, nx = D["2023-03-08"]["B04"].shape
    E, N = np.meshgrid(BOX[0] + (np.arange(nx) + .5) * 0.03, BOX[3] - (np.arange(ny) + .5) * 0.03)
    new_lava = {}
    for d, d0 in PAIRS.items():
        s, pre = D[d], D[d0]
        m = (s["B8A"] - pre["B8A"] < -0.06) & (s["B8A"] < 0.10) & ~s["bad"] & ~pre["bad"]
        m = ndimage.binary_opening(m, iterations=1)
        lab, nl = ndimage.label(m)
        sizes = ndimage.sum(m, lab, range(1, nl + 1))
        m = np.isin(lab, 1 + np.flatnonzero(sizes >= 30))                    # >= 30 pixels = 2.7 ha
        new_lava[d] = m
        print(d, f"new lava {m.sum() * 0.0009:.1f} km2", flush=True)
    s = D["2024-04-26"]
    hot = (s["B12"] > 0.25) & (s["B12"] > 1.2 * s["B11"]) & ~s["bad"]
    lab, nh = ndimage.label(ndimage.binary_dilation(hot, iterations=2))
    vents = []
    for i in range(1, nh + 1):
        k = (lab == i) & hot
        if k.sum() >= 2:
            w = s["B12"][k]
            vents.append((float((E[k] * w).sum() / w.sum()), float((N[k] * w).sum() / w.sum()), int(k.sum())))
    vents.sort(key=lambda v: -v[2])
    print("hot vents (east, north km, pixels):", [(round(a, 2), round(b, 2), c) for a, b, c in vents[:8]], flush=True)
    traces = dike_traces()
    # upstream ends of the new lava: lava pixels within 400 m of a hot vent or, for the March scene,
    # the highest-NIR-drop pixels nearest to the crater row are not separable -> use all new-lava pixels
    lava_pts = np.c_[E[new_lava["2024-03-02"]], N[new_lava["2024-03-02"]]]
    stats = {}
    for nm, (x0, y0, x1, y1) in traces.items():
        dv = [float(dist_to_segment(v[0], v[1], x0, y0, x1, y1)) for v in vents[:5]]
        dl = dist_to_segment(lava_pts[:, 0], lava_pts[:, 1], x0, y0, x1, y1) if len(lava_pts) else np.array([np.nan])
        stats[nm] = dict(vent_distance_km=[round(v, 2) for v in dv], lava_median_distance_km=round(float(np.median(dl)), 2),
                         lava_within_1km=round(float(np.mean(dl < 1.0)), 3))
        print(f"{nm:28s} vents {stats[nm]['vent_distance_km']} km; new lava median distance {stats[nm]['lava_median_distance_km']} km, "
              f"{stats[nm]['lava_within_1km']:.0%} within 1 km", flush=True)
    # orientation of the new-lava field (principal axis) vs dike strike
    if len(lava_pts) > 50:
        c = np.cov((lava_pts - lava_pts.mean(0)).T)
        w, v = np.linalg.eigh(c)
        ax = v[:, -1]
        strike_lava = float(np.degrees(np.arctan2(ax[0], ax[1])) % 180)
    else:
        strike_lava = None
    out = dict(scenes=SCENES, new_lava_km2={d: round(float(m.sum() * 0.0009), 2) for d, m in new_lava.items()},
               vents=vents[:8], lava_axis_deg=strike_lava, dikes=stats,
               note="HLS S30 v2.0 (NASA LP DAAC); new lava = NIR drop > 0.06 and NIR < 0.10; vents = B12 > 0.25 and B12 > 1.2 B11")
    (ML / "hls_dike_optical.json").write_text(json.dumps(out, indent=1, ensure_ascii=False), encoding="utf-8")
    print("principal axis of the new lava field:", strike_lava)

    # ---- figure
    fig, axs = plt.subplots(1, 4, figsize=(19, 5.6))
    show = ("2023-03-08", "2024-03-02", "2024-04-26")
    ext = [BOX[0], BOX[1], BOX[2], BOX[3]]
    pl = {k: ((to_utm.transform(*PLACES[k])[0] - sx) / 1000, (to_utm.transform(*PLACES[k])[1] - sy) / 1000)
          for k in ("Svartsengi", "Grindavik", "Sundhnukur")}
    ua = {"Svartsengi": "Свартсенгі", "Grindavik": "Ґріндавік", "Sundhnukur": "Сундгнукур"}
    ttl = {"2023-03-08": "08.03.2023 — до дайки", "2024-03-02": "02.03.2024 — після 3 вивержень",
           "2024-04-26": "26.04.2024 — під час виверження"}
    for a, d in zip(axs[:3], show):
        a.imshow(rgb(D[d]), extent=ext, origin="upper")
        a.set_title(f"{ttl[d]}\nSWIR: B12, B11, B4", fontsize=9)
    a = axs[3]
    a.imshow(rgb(D["2024-03-02"]) * 0.35 + 0.65, extent=ext, origin="upper")
    a.contourf(E, N, new_lava["2024-03-02"].astype(float), levels=[0.5, 1.5], colors=["#222222"], alpha=.85)
    a.contourf(E, N, (new_lava["2024-04-26"] & ~new_lava["2024-03-02"]).astype(float), levels=[0.5, 1.5], colors=["#7a5195"], alpha=.6)
    sty = {"NPE, однакове розкриття": dict(color="tab:blue", ls="-"), "MCMC, однакове розкриття": dict(color="tab:blue", ls="--"),
           "NPE, змінне розкриття": dict(color="tab:red", ls="-"), "MCMC, змінне розкриття": dict(color="tab:red", ls="--")}
    for nm, (x0, y0, x1, y1) in traces.items():
        a.plot([x0, x1], [y0, y1], lw=2, label=nm, **sty.get(nm, {}))
    if vents:
        v = np.array(vents[:5])
        a.plot(v[:, 0], v[:, 1], "*", color="orange", mec="k", ms=13, label="активні жерла 26.04.2024")
    a.set_title("нова лава (чорна: до 02.03; фіолетова: до 26.04)\nі дайка з Sentinel-1", fontsize=9)
    a.legend(fontsize=7, loc="lower left", framealpha=.9)
    for a in axs:
        for k, (x0, y0) in pl.items():
            a.plot(x0, y0, "w^" if k == "Svartsengi" else "wo", mec="k", ms=6)
            a.text(x0 + (-0.35 if k == "Svartsengi" else 0.35), y0 + (0.6 if k == "Sundhnukur" else -0.9 if k == "Svartsengi" else 0.35), ua[k],
                   ha="right" if k == "Svartsengi" else "left", fontsize=8, color="k", bbox=dict(fc="white", alpha=.6, lw=0, pad=1))
        a.set_xlim(BOX[0], BOX[1]); a.set_ylim(BOX[2], BOX[3])
        a.set_xlabel("схід від Свартсенгі, км", fontsize=8)
        a.tick_params(labelsize=7)
    axs[0].set_ylabel("північ, км", fontsize=8)
    fig.suptitle("Оптична перевірка дайки: NASA HLS (Sentinel-2, 30 м) — свіжа лава на кратерному ряді Сундгнукур", fontsize=12)
    fig.savefig(ML / "hls_dike_optical.png", dpi=105, facecolor="white", bbox_inches="tight")
    print("wrote", ML / "hls_dike_optical.png")


if __name__ == "__main__":
    main()
