"""Near-field data for the 10 Nov 2023 dike: Sentinel-1 pixel offsets (ASF HyP3 AUTORIFT, 120 m), binned on 1 km cells.

AUTORIFT gives the velocity in radar range (vr) and azimuth (va) directions (m/yr) over the 12-day pair; displacement =
v x dt / 365.25. Offsets work where interferometry loses coherence (the graben, Grindavik), and azimuth offsets add the
along-track horizontal direction that InSAR cannot see. Errors (product metadata): range ~8.9 m/yr (~0.29 m per 120 m
pixel), azimuth ~48 m/yr (~1.6 m per pixel); 1-km medians of ~70 pixels are ~10x better.

Steps: (1) bin both components on 1 km cells around Svartsengi (same frame as the InSAR bins); (2) sign / scale check:
range offsets vs the unwrapped InSAR LOS outside the graben (expected slope +-1); (3) unit vectors: range = ground-to-
satellite look vector (sign from step 2), azimuth = horizontal flight direction (right-looking SAR: look direction rotated
90 deg counter-clockwise). Writes data/sentinel1/s1_offsets_bins.npz and output/reykjanes/ml/s1_offsets.json / .png.
"""
import glob
import json

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from pyproj import Transformer

import s1_grindavik_dike as S
from reykjanes_analysis import PLACES

ML = S.ML
FILES = {"A016": "*20231031T185905*_X_*20231112T185905*.nc", "D155": "*20231110T075906*_X_*20231122T075906*.nc"}


def load(track):
    f = sorted(glob.glob(str(S.S1 / FILES[track])))
    if not f:
        return None
    h = h5py.File(f[0], "r")
    x, y = h["x"][:], h["y"][:]
    dt = float(h["img_pair_info"].attrs["date_dt"])
    out = {}
    for k in ("vr", "va"):
        v = h[k][0].astype(float)
        v[v == h[k].attrs["_FillValue"]] = np.nan
        out[k] = v * dt / 365.25 * 100                         # cm over the pair
    X, Y = np.meshgrid(x, y)
    to_utm = Transformer.from_crs(3413, 32627, always_xy=True)
    ux, uy = to_utm.transform(X, Y)
    sx, sy = Transformer.from_crs(4326, 32627, always_xy=True).transform(*PLACES["Svartsengi"])
    return dict(e=(ux - sx) / 1000, n=(uy - sy) / 1000, dt=dt, err=dict(vr=float(h["vr"].attrs["error"]) * dt / 365.25 * 100,
                                                                         va=float(h["va"].attrs["error"]) * dt / 365.25 * 100), **out)


def main():
    z = np.load(S.S1 / "s1_bins.npz")
    res, bins = {}, {}
    fig, axs = plt.subplots(2, 3, figsize=(17, 10))
    for k, t in enumerate(S.TRACKS):
        d = load(t)
        if d is None:
            print("no offsets for", t); continue
        box = (d["e"] > S.BOX[0]) & (d["e"] < S.BOX[1]) & (d["n"] > S.BOX[2]) & (d["n"] < S.BOX[3])
        key = np.c_[np.floor(d["e"][box] / S.BIN), np.floor(d["n"][box] / S.BIN)].astype(int)
        cells = {}
        for comp in ("vr", "va"):
            v = d[comp][box]
            ok = np.isfinite(v)
            kk, vv = key[ok], v[ok]
            order = np.lexsort((kk[:, 1], kk[:, 0]))
            kk, vv = kk[order], vv[order]
            brk = np.flatnonzero(np.any(np.diff(kk, axis=0) != 0, axis=1)) + 1
            for a, b in zip(np.split(kk, brk), np.split(vv, brk)):
                if len(b) >= 20:
                    cells.setdefault((a[0, 0], a[0, 1]), {})[comp] = (float(np.median(b)), 1.2533 * d["err"][comp] / np.sqrt(len(b)), len(b))
        keys = sorted(c for c in cells if "vr" in cells[c] and "va" in cells[c])
        E = np.array([(a + 0.5) * S.BIN for a, _ in keys]); N = np.array([(b + 0.5) * S.BIN for _, b in keys])
        vr = np.array([cells[c]["vr"][0] for c in keys]); va = np.array([cells[c]["va"][0] for c in keys])
        sr = np.array([cells[c]["vr"][1] for c in keys]); sa = np.array([cells[c]["va"][1] for c in keys])
        # InSAR at the same cells (LOS toward the satellite, cm)
        ins = z[t]
        ikey = {(int(np.floor(e / S.BIN)), int(np.floor(n / S.BIN))): (los, vec) for e, n, los, *vec in ins[:, :6]}
        both = np.array([c in ikey for c in keys])
        L = np.array([ikey[c][0] for c, b in zip(keys, both) if b])
        vr_b = vr[both]
        A = np.c_[L - np.median(L), np.ones(len(L))]
        slope = float(np.linalg.lstsq(A, vr_b - np.median(vr_b), rcond=None)[0][0])
        corr = float(np.corrcoef(L, vr_b)[0, 1])
        vec_mean = ins[:, 3:6].mean(0)
        h = vec_mean[:2] / np.linalg.norm(vec_mean[:2])                    # horizontal ground-to-satellite
        look = -h
        flight = np.array([-look[1], look[0]])                              # rotate look 90 deg counter-clockwise
        res[t] = dict(n_cells=len(keys), n_with_insar=int(both.sum()), n_offsets_only=int((~both).sum()),
                      slope_range_vs_los=round(slope, 3), corr_range_vs_los=round(corr, 3),
                      range_sigma_cell_cm_median=round(float(np.median(sr)), 1), azimuth_sigma_cell_cm_median=round(float(np.median(sa)), 1),
                      range_p1_p99_cm=np.percentile(vr, [1, 99]).round(1).tolist(), azimuth_p1_p99_cm=np.percentile(va, [1, 99]).round(1).tolist(),
                      flight_unit_en=flight.round(3).tolist())
        sign = -1.0 if slope < 0 else 1.0                                   # range offset = sign * LOS(toward)
        bins[t] = np.c_[E, N, vr, sr, va, sa]
        bins[t + "_meta"] = np.r_[sign, flight, vec_mean]
        print(f"{t}: {len(keys)} offset cells ({(~both).sum()} where InSAR has no data); range offset vs InSAR LOS: slope {slope:+.2f}, "
              f"corr {corr:+.2f}; cell sigma range {np.median(sr):.1f} cm, azimuth {np.median(sa):.1f} cm; range p1/p99 {res[t]['range_p1_p99_cm']} cm", flush=True)
        lim = np.percentile(np.abs(vr), 98)
        for j, (v, nm) in enumerate(((vr, "зсув пікселів по дальності"), (va, "зсув пікселів по азимуту"))):
            a = axs[k, j]
            sc = a.scatter(E, N, c=v, cmap="RdBu_r", vmin=-lim, vmax=lim, s=12, marker="s", linewidths=0)
            a.scatter(E[~both], N[~both], facecolors="none", edgecolors="k", s=12, marker="s", linewidths=.4)
            a.set_aspect("equal"); a.set_xlim(S.BOX[0], S.BOX[1]); a.set_ylim(S.BOX[2], S.BOX[3]); a.tick_params(labelsize=7)
            a.set_title(f"{t}: {nm}, см (обведені — де InSAR без даних)", fontsize=9)
            plt.colorbar(sc, ax=a, fraction=.035).set_label("см", fontsize=8)
        a = axs[k, 2]
        a.scatter(L, vr_b, s=4, alpha=.5)
        xx = np.linspace(L.min(), L.max(), 10)
        a.plot(xx, sign * (xx - np.median(L)) + np.median(vr_b), "k--", lw=.8, label="нахил ±1")
        a.set_xlabel("InSAR LOS до супутника, см", fontsize=8); a.set_ylabel("зсув по дальності, см", fontsize=8)
        a.set_title(f"{t}: нахил {slope:+.2f}, кореляція {corr:+.2f}", fontsize=9); a.legend(fontsize=7); a.tick_params(labelsize=7)
    np.savez(S.S1 / "s1_offsets_bins.npz", **bins)
    (ML / "s1_offsets.json").write_text(json.dumps(res, indent=1), encoding="utf-8")
    fig.suptitle("Зсуви пікселів Sentinel-1 (ASF HyP3 AUTORIFT) для дайки 10.11.2023, клітинки 1 км", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(ML / "s1_offsets.png", dpi=95, facecolor="white")
    print("wrote", ML / "s1_offsets.png")


if __name__ == "__main__":
    main()
