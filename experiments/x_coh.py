"""X-coh (pre-registered in notes/experiments.md): is cropland coherence a crop calendar? Central Valley asc.

    uv run python scripts/x_coh.py
"""
from __future__ import annotations

import json
from datetime import timedelta
from pathlib import Path

import numpy as np
import rasterio

from hachure import gridding, nisar

ROOT = Path(__file__).resolve().parents[1]
BBOX = (-120.20, 35.80, -119.00, 36.50)
CACHE = ROOT / "data/raw/nisar_cache" / ("bbox_" + "_".join(f"{v:.3f}" for v in BBOX) + "_h80")
OUT = ROOT / "data/processed/central-valley-asc"
B = 25


def blocks(a):
    H, W = (a.shape[0] // B) * B, (a.shape[1] // B) * B
    with np.errstate(all="ignore"):
        return np.nanmedian(a[:H, :W].reshape(H // B, B, W // B, B).transpose(0, 2, 1, 3).reshape(H // B, W // B, -1), -1)


def main():
    found = nisar.search((BBOX[0] + BBOX[2]) / 2, (BBOX[1] + BBOX[3]) / 2, direction="ASCENDING")
    chain = nisar.sequential_chain(nisar.dominant_track(found))
    with rasterio.open(OUT / "nisar_rate.tif") as ds:
        tr, crs, shape = ds.transform, ds.crs, ds.shape
    lc = np.load(ROOT / "data/raw/context/cv_worldcover_frac.npz")
    crop, built = blocks(lc["cropland"]), blocks(lc["builtup"])
    rain = json.load(open(ROOT / "data/raw/context/cv_rain_openmeteo.json"))["daily"]
    rmap = {d: (p or 0.0) for d, p in zip(rain["time"], rain["precipitation_sum"])}

    def wet(d):
        return sum(rmap.get((d - timedelta(days=i)).strftime("%Y-%m-%d"), 0) for i in range(0, 4)) >= 2

    cm, bm = crop > 0.6, built > 0.3
    rows, curves = [], []
    for g in chain:
        f = CACHE / f"{g.name}.npz"
        if not f.exists():
            continue
        raw = nisar.RawWindow.load(f)
        c = gridding.to_grid(raw.coherence.astype(float), gridding.transform_from_centres(raw.x, raw.y),
                             f"EPSG:{raw.epsg}", tr, crs, shape)
        cb = blocks(c)
        r0, r1 = nisar._parse_dates(g.name)
        mid = r0 + (r1 - r0) / 2
        R = float(np.nanmedian(cb[cm]) / np.nanmedian(cb[bm]))
        rows.append(dict(ref=f"{r0:%Y-%m-%d}", sec=f"{r1:%Y-%m-%d}", mid_month=mid.month, R=R,
                         crop=float(np.nanmedian(cb[cm])), built=float(np.nanmedian(cb[bm])), wet=bool(wet(r0) or wet(r1))))
        curves.append(cb[cm])
        r = rows[-1]
        print(f"  {r['ref']} -> {r['sec']} {'WET' if r['wet'] else '   '} crop {r['crop']:.3f} built {r['built']:.3f} R {R:.3f}")
    R = np.array([r["R"] for r in rows])
    mm = np.array([r["mid_month"] for r in rows])
    summer, winter = R[np.isin(mm, [5, 6, 7, 8, 9])], R[np.isin(mm, [11, 12, 1, 2, 3])]
    x1 = bool(np.median(summer) < np.median(winter))
    ac = float(np.corrcoef(R[:-1], R[1:])[0, 1])
    x2 = ac >= 0.3
    w = np.array([r["wet"] for r in rows])
    dwet = float(R[~w].mean() - R[w].mean()) if w.any() and (~w).any() else float("nan")
    x3 = bool(dwet >= 0.05)
    print(f"{len(rows)} pairs; cropland blocks {cm.sum()}, built-up blocks {bm.sum()}")
    print(f"X1: median R summer {np.median(summer):.3f} (n {len(summer)}) vs winter {np.median(winter):.3f} (n {len(winter)}) -> {'PASS' if x1 else 'FAIL'}")
    print(f"X2: lag-1 autocorrelation {ac:+.2f} (>= 0.3) -> {'PASS' if x2 else 'FAIL'}")
    print(f"X3: mean R dry - wet {dwet:+.3f} (>= 0.05; wet n {w.sum()}) -> {'PASS' if x3 else 'FAIL'}")
    # post-hoc exploratory: k-means on per-block cropland curves
    Y = np.array(curves).T                       # blocks x pairs
    ok = np.isfinite(Y).all(1)
    Y = Y[ok]
    Z = (Y - Y.mean(1, keepdims=True))
    rng = np.random.default_rng(0)
    C = Z[rng.choice(len(Z), 4, replace=False)]
    for _ in range(50):
        lab = np.argmin(((Z[:, None] - C[None]) ** 2).sum(-1), 1)
        C = np.array([Z[lab == k].mean(0) if (lab == k).any() else C[k] for k in range(4)])
    print("post-hoc k-means (k=4) on mean-removed cropland block curves:")
    for k in range(4):
        cur = Y[lab == k].mean(0)
        print(f"  cluster {k}: {int((lab == k).sum())} blocks, mean coh {cur.mean():.2f}, "
              f"min at {rows[int(np.argmin(cur))]['sec']}, max at {rows[int(np.argmax(cur))]['sec']}")
    json.dump(dict(rows=rows, X1=x1, X2=dict(ac=ac, passed=x2), X3=dict(dry_minus_wet=dwet, passed=x3),
                   clusters=[dict(n=int((lab == k).sum()), curve=Y[lab == k].mean(0).tolist()) for k in range(4)]),
              open(OUT / "x_coh.json", "w"), indent=1)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from datetime import datetime
    t = [datetime.fromisoformat(r["sec"]) for r in rows]
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.5))
    ax[0].plot(t, R, "o-", label="R = cropland / built-up coherence")
    for i, r in enumerate(rows):
        if r["wet"]:
            ax[0].axvline(t[i], color="b", alpha=0.15)
    ax[0].legend(); ax[0].tick_params(axis="x", rotation=45)
    for k in range(4):
        ax[1].plot(t, Y[lab == k].mean(0), label=f"cluster {k} ({(lab == k).sum()})")
    ax[1].legend(); ax[1].tick_params(axis="x", rotation=45); ax[1].set_title("post-hoc: cropland coherence 'dances'")
    fig.tight_layout(); fig.savefig(OUT / "x_coh.png", dpi=110)


if __name__ == "__main__":
    main()
