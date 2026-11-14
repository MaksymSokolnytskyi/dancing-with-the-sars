"""Fissure / graben maps from the loss of InSAR coherence (HyP3 Sentinel-1 _corr.tif), verified against the Sundhnukur crater row.

For each event pair the coherence is compared with the quiet pair of the same track (7/5 -> 19/17 Oct 2023) on a common 100 m
grid: drop = corr(event) - corr(quiet). Cells with drop < -0.25 and corr(event) < 0.35 are "lost". Connected lost zones are
classified by shape (principal axes of the cell cloud): ELONGATED (length > 1 km, axis ratio > 5) -> fissure / graben candidate,
compact and large -> lava field / other. The two tracks are combined (a cell counts if lost on either track).

Verification: distance of every fissure candidate to the Sundhnukur crater row line (geological maps), to the 10 Nov 2023 dike
trace (NPE v3) and to the HLS hot vents; orientation compared with the NE-SW trend of the fissure swarm (~35-40 deg).
Writes output/reykjanes/ml/coherence_fissures.json and coherence_fissures.png.
"""
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from pyproj import Transformer
from scipy import ndimage

import s1_grindavik_dike as S
from firms_hotspots import CRATER_ROW
from reykjanes_analysis import PLACES
from s1_common_v2 import product_dir

ML = S.ML
CELL = 0.1                       # km
GRID = (-12, 12, -12, 10)        # km from Svartsengi
EVENTS = {"real_20231110": None, "event_20240208": None, "event_20240302": None, "event_20240316": None, "event_20240529": None}
LAB = {"real_20231110": "дайка 10.11.2023", "event_20240208": "08.02.2024", "event_20240302": "02.03.2024 (без виверження)",
       "event_20240316": "16.03.2024", "event_20240529": "29.05.2024"}


def to_grid(d):
    v, x, y, epsg = S.read(d, "corr")
    to_utm = Transformer.from_crs(4326, 32627, always_xy=True)
    sx, sy = to_utm.transform(*PLACES["Svartsengi"])
    ux, uy = Transformer.from_crs(epsg, 32627, always_xy=True).transform(x, y)
    e, n = (ux - sx) / 1000, (uy - sy) / 1000
    ok = np.isfinite(v) & (v > 0)
    nx, ny = int((GRID[1] - GRID[0]) / CELL), int((GRID[3] - GRID[2]) / CELL)
    ix = ((e[ok] - GRID[0]) / CELL).astype(int); iy = ((GRID[3] - n[ok]) / CELL).astype(int)
    m = (ix >= 0) & (ix < nx) & (iy >= 0) & (iy < ny)
    s = np.zeros((ny, nx)); c = np.zeros((ny, nx))
    np.add.at(s, (iy[m], ix[m]), v[ok][m]); np.add.at(c, (iy[m], ix[m]), 1)
    with np.errstate(invalid="ignore"):
        return np.where(c >= 4, s / c, np.nan)


def seg_dist(px, py, x0, y0, x1, y1):
    vx, vy = x1 - x0, y1 - y0
    t = np.clip(((px - x0) * vx + (py - y0) * vy) / (vx * vx + vy * vy), 0, 1)
    return np.hypot(px - (x0 + t * vx), py - (y0 + t * vy))


def main():
    idx = json.loads((S.S1 / "pairs_index.json").read_text())
    pairs = {"real_20231110": [dict(track="A016", ref=S.PAIRS["A016"].split("_")[0], sec=S.PAIRS["A016"].split("_")[1]),
                               dict(track="D155", ref=S.PAIRS["D155"].split("_")[0], sec=S.PAIRS["D155"].split("_")[1])]}
    pairs.update({k: v for k, v in idx.items() if k in EVENTS})
    quiet = {p["track"]: to_grid(product_dir(p["track"], p["ref"], p["sec"])) for p in idx["quiet_2023"]}
    to_utm = Transformer.from_crs(4326, 32627, always_xy=True)
    sx, sy = to_utm.transform(*PLACES["Svartsengi"])
    (a0, b0), (a1, b1) = to_utm.transform(*CRATER_ROW[0]), to_utm.transform(*CRATER_ROW[1])
    row = ((a0 - sx) / 1000, (b0 - sy) / 1000, (a1 - sx) / 1000, (b1 - sy) / 1000)
    pm = np.median(np.load(ML / "s1_npe_v3.npz")["ps"], 0)
    u = np.array([np.sin(np.deg2rad(pm[2])), np.cos(np.deg2rad(pm[2]))])
    dike = (pm[0] - u[0] * pm[3] / 2, pm[1] - u[1] * pm[3] / 2, pm[0] + u[0] * pm[3] / 2, pm[1] + u[1] * pm[3] / 2)
    ny, nx = quiet["A016"].shape
    E, N = np.meshgrid(GRID[0] + (np.arange(nx) + .5) * CELL, GRID[3] - (np.arange(ny) + .5) * CELL)
    res, maps = {}, {}
    for name, pl in pairs.items():
        lost = np.zeros((ny, nx), bool)
        for p in pl:
            c = to_grid(product_dir(p["track"], p["ref"], p["sec"]))
            q = quiet[p["track"]]
            lost |= np.nan_to_num((c - q) < -0.35) & np.nan_to_num(c < 0.25) & np.isfinite(q) & (q > 0.5)
        lost = ndimage.binary_opening(lost, iterations=1)
        lab, nl = ndimage.label(ndimage.binary_closing(lost, iterations=1))
        fis, blobs = [], []
        for i in range(1, nl + 1):
            k = lab == i
            if k.sum() < 50:
                continue
            pts = np.c_[E[k], N[k]]
            c0 = pts.mean(0)
            w, v = np.linalg.eigh(np.cov((pts - c0).T))
            ratio = np.sqrt(w[1] / max(w[0], 1e-6))
            ax = v[:, 1]
            proj = (pts - c0) @ ax
            length = float(proj.max() - proj.min())
            strike = float(np.degrees(np.arctan2(ax[0], ax[1])) % 180)
            item = dict(e=round(float(c0[0]), 2), n=round(float(c0[1]), 2), length_km=round(length, 2), axis_ratio=round(float(ratio), 1),
                        strike_deg=round(strike, 0), area_km2=round(float(k.sum() * CELL ** 2), 2),
                        dist_crater_row_km=round(float(np.median(seg_dist(pts[:, 0], pts[:, 1], *row))), 2),
                        dist_dike_km=round(float(np.median(seg_dist(pts[:, 0], pts[:, 1], *dike))), 2))
            (fis if (length > 3.0 and ratio > 2.5) else blobs).append(item)
        fis.sort(key=lambda d: -d["length_km"])
        dr = seg_dist(E, N, *row)
        buf = dr < 1.5
        valid = np.isfinite(quiet["A016"]) & (quiet["A016"] > 0.5)
        frac_lost_in_buf = float(lost[buf & valid].sum() / max(lost[valid].sum(), 1))
        frac_area_buf = float((buf & valid).sum() / valid.sum())
        res[name] = dict(enrichment_crater_row_1p5km=round(frac_lost_in_buf / frac_area_buf, 1), lost_share_within_1p5km=round(frac_lost_in_buf, 3),
                         area_share_within_1p5km=round(frac_area_buf, 3), n_fissures=len(fis), fissures=fis[:12], total_fissure_length_km=round(sum(f["length_km"] for f in fis), 1),
                         n_other=len(blobs), lost_area_km2=round(float(lost.sum() * CELL ** 2), 1),
                         fissure_median_dist_crater_row_km=round(float(np.median([f["dist_crater_row_km"] for f in fis])), 2) if fis else None,
                         fissure_median_dist_dike_km=round(float(np.median([f["dist_dike_km"] for f in fis])), 2) if fis else None,
                         fissure_strike_median=round(float(np.median([f["strike_deg"] for f in fis])), 0) if fis else None)
        maps[name] = lab > 0
        print(f"{name}: enrichment within 1.5 km of the crater row x{res[name]['enrichment_crater_row_1p5km']} "
              f"({frac_lost_in_buf:.0%} of the lost area in {frac_area_buf:.0%} of the map); {len(fis)} elongated zones ({res[name]['total_fissure_length_km']} km), median distance to the crater row "
              f"{res[name]['fissure_median_dist_crater_row_km']} km, to the dike {res[name]['fissure_median_dist_dike_km']} km, "
              f"strike {res[name]['fissure_strike_median']}; lost area {res[name]['lost_area_km2']} km2", flush=True)
    # random-line control: how close would randomly placed lines be to the crater row?
    rng = np.random.default_rng(0)
    rnd = seg_dist(rng.uniform(GRID[0], GRID[1], 5000), rng.uniform(GRID[2], GRID[3], 5000), *row)
    out = dict(cell_km=CELL, events=res, control_random_points_median_dist_crater_row_km=round(float(np.median(rnd)), 2),
               crater_row=CRATER_ROW, dike_trace_npe_v3=[round(float(v), 2) for v in dike])
    (ML / "coherence_fissures.json").write_text(json.dumps(out, indent=1, ensure_ascii=False), encoding="utf-8")
    print("control: random points median distance to the crater row", out["control_random_points_median_dist_crater_row_km"], "km")
    fig, axs = plt.subplots(1, len(maps), figsize=(4.2 * len(maps), 4.6))
    for a, (name, m) in zip(axs, maps.items()):
        a.imshow(m, extent=GRID, cmap="Greys", vmin=0, vmax=1.4, origin="upper")
        for f in res[name]["fissures"]:
            th = np.deg2rad(f["strike_deg"])
            a.plot([f["e"] - np.sin(th) * f["length_km"] / 2, f["e"] + np.sin(th) * f["length_km"] / 2],
                   [f["n"] - np.cos(th) * f["length_km"] / 2, f["n"] + np.cos(th) * f["length_km"] / 2], color="tab:blue", lw=2)
        a.plot([row[0], row[2]], [row[1], row[3]], "r-", lw=1.5, label="кратерний ряд Сундгнукур")
        a.plot([dike[0], dike[2]], [dike[1], dike[3]], "k--", lw=1, label="слід дайки (NPE v3)")
        for k in ("Svartsengi", "Grindavik"):
            px, py = to_utm.transform(*PLACES[k])
            a.plot((px - sx) / 1000, (py - sy) / 1000, "k^" if k == "Svartsengi" else "ko", ms=4)
        a.set_title(f"{LAB[name]}\nу смузі ±1.5 км від кратерного ряду: {res[name]['lost_share_within_1p5km']:.0%} втраченої площі (×{res[name]['enrichment_crater_row_1p5km']})",
                    fontsize=8)
        a.set_xlim(GRID[0], GRID[1]); a.set_ylim(GRID[2], GRID[3]); a.tick_params(labelsize=7)
    axs[0].legend(fontsize=7, loc="lower left")
    fig.suptitle("Тріщини з втрати когерентності Sentinel-1 (сірі — зони втрати; сині — осі видовжених зон); ×N — у скільки разів частіше, ніж випадково", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(ML / "coherence_fissures.png", dpi=100, facecolor="white")
    print("wrote", ML / "coherence_fissures.png")


if __name__ == "__main__":
    main()
