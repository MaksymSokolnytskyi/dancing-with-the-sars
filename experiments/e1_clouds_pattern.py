"""E1: does a pair's radar residual map look like the GOES cloud-change map? (pre-registered)

Residual (leave-one-out, as atmo_stats) -> 2 km blocks; GOES band-13 dBT = BT(date2) - BT(date1)
interpolated to the same blocks with a shift (dx, dy) applied to the GOES map; plane removed
from both; Pearson r. One common shift for all cloudy pairs; null = shuffled dBT maps.

    uv run python scripts/e1_clouds_pattern.py
"""
from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path

import numpy as np
from rasterio.warp import transform as warp_transform
from scipy.interpolate import LinearNDInterpolator

from hachure import nisar

ROOT = Path(__file__).resolve().parents[1]
BBOX = (-99.32, 19.17, -98.82, 19.58)
REF_BBOX = (-99.14, 19.53, -99.09, 19.57)
GOES = ROOT / "data" / "raw" / "goes19"
OUT = ROOT / "data" / "processed" / "mexico-city"
BLOCK = 25                                    # 25 x 80 m = 2 km
SHIFTS = np.arange(-15, 16) * 1000.0          # m
EXPECT_AZ = 126.9                             # toward the GOES-19 sub-satellite point
ORIGIN = "upper"


def residuals(direction, series):
    found = nisar.search((BBOX[0] + BBOX[2]) / 2, (BBOX[1] + BBOX[3]) / 2, direction=direction)
    chain = nisar.sequential_chain(nisar.dominant_track(found))
    cache = ROOT / "data" / "raw" / "nisar_cache" / ("bbox_" + "_".join(f"{v:.3f}" for v in BBOX) + "_h2240")
    pairs = [nisar.process(nisar.RawWindow.load(cache / f"{g.name}.npz"), ref_bbox=REF_BBOX, coh_min=0.3)
             for g in chain]
    disps, _, x, y = nisar.align(pairs)
    days = np.array([p.days for p in pairs])
    ok = np.isfinite(disps)
    tot = np.where(ok, disps, 0).sum(0)
    span = np.where(ok, days[:, None, None], 0).sum(0)
    z = np.load(ROOT / "data" / "processed" / series / "nisar_series.npz")
    dates = [datetime.fromisoformat(str(d)) for d in z["dates"]]
    assert len(dates) == len(pairs) + 1, (len(dates), len(pairs))
    out = []
    for i in range(len(pairs)):
        with np.errstate(all="ignore"):
            rl = (tot - np.where(ok[i], disps[i], 0)) / (span - days[i])
        r = np.where(ok[i] & (span - days[i] > 200), disps[i] - rl * days[i], np.nan)
        out.append((dates[i], dates[i + 1], r - np.nanmedian(r)))
    return out, x, y, pairs[0].epsg


def blocks(a, x, y):
    H, W = (a.shape[0] // BLOCK) * BLOCK, (a.shape[1] // BLOCK) * BLOCK
    v = a[:H, :W].reshape(H // BLOCK, BLOCK, W // BLOCK, BLOCK)
    n = np.isfinite(v).sum((1, 3))
    with np.errstate(all="ignore"):
        m = np.nanmean(v, (1, 3))
    m[n < 0.4 * BLOCK * BLOCK] = np.nan
    bx = x[:W].reshape(-1, BLOCK).mean(1)
    by = y[:H].reshape(-1, BLOCK).mean(1)
    return m, bx, by


def goes_interp(t, epsg):
    fs = sorted(GOES.glob("bt13_*.npz"))
    ts = [datetime.strptime(f.stem[5:], "%Y%m%dT%H%M") for f in fs]
    k = int(np.argmin([abs((s - t).total_seconds()) for s in ts]))
    if abs((ts[k] - t).total_seconds()) > 3600:
        return None
    g = np.load(fs[k])
    bt, lat, lon = g["bt"], g["lat"], g["lon"]
    m = np.isfinite(bt) & np.isfinite(lat)
    gx, gy = warp_transform("EPSG:4326", f"EPSG:{epsg}", lon[m].tolist(), lat[m].tolist())
    return LinearNDInterpolator(np.c_[gx, gy], bt[m]), bt, lat, lon


def plane_free(a, mask):
    yy, xx = np.nonzero(mask)
    P = np.c_[np.ones(len(yy)), xx, yy]
    c, *_ = np.linalg.lstsq(P, a[yy, xx], rcond=None)
    out = np.full(a.shape, np.nan)
    out[yy, xx] = a[yy, xx] - P @ c
    return out


def cold_frac(bt, lat, lon):
    inbox = (lat >= BBOX[1]) & (lat <= BBOX[3]) & (lon >= BBOX[0]) & (lon <= BBOX[2])
    return float(np.mean(bt[inbox] < 235))


def run(direction, series):
    res, x, y, epsg = residuals(direction, series)
    _, bx, by = blocks(res[0][2], x, y)
    global ORIGIN
    ORIGIN = "upper" if by[0] > by[-1] else "lower"
    BX, BY = np.meshgrid(bx, by)
    full = np.ones(BX.shape, bool)
    cache = {}
    rows = []
    for d0, d1, r in res:
        for d in (d0, d1):
            if d not in cache:
                cache[d] = goes_interp(d, epsg)
        if cache[d0] is None or cache[d1] is None:
            continue
        rb, _, _ = blocks(r, x, y)
        mask = np.isfinite(rb)
        if mask.sum() < 40:
            continue
        rp = plane_free(rb, mask)
        # dBT for every shift: GOES map moved by (sx, sy) => sample GOES at (X - sx, Y - sy)
        dbt = np.empty((len(SHIFTS), len(SHIFTS)) + BX.shape)
        for a, sy in enumerate(SHIFTS):
            for b, sx in enumerate(SHIFTS):
                q = np.c_[(BX - sx).ravel(), (BY - sy).ravel()]
                v = cache[d1][0](q) - cache[d0][0](q)
                dbt[a, b] = plane_free(v.reshape(BX.shape), full)
        cf = max(cold_frac(*cache[d0][1:]), cold_frac(*cache[d1][1:]))
        rows.append(dict(d0=d0, rp=rp, mask=mask, dbt=dbt, cold=cf))
    return rows


def corr_cube(rp, mask, dbt):
    """Pearson r of rp (masked) with every shifted dbt map -> (S, S)."""
    a = rp[mask]
    a = (a - a.mean()) / a.std()
    b = dbt[:, :, mask]
    b = (b - b.mean(-1, keepdims=True)) / b.std(-1, keepdims=True)
    return (b * a).mean(-1)


def analyse(rows, name, rng):
    cloudy = [r for r in rows if r["cold"] > 0.05]
    clear = [r for r in rows if r["cold"] <= 0.05]
    n = len(cloudy)
    R = np.array([[corr_cube(ri["rp"], ri["mask"], rj["dbt"]) for rj in cloudy] for ri in cloudy])  # (n,n,S,S)
    own = np.median(R[np.arange(n), np.arange(n)], 0)                       # (S,S)
    k = np.unravel_index(np.argmax(np.abs(own)), own.shape)
    best_r = own[k]
    sy, sx = SHIFTS[k[0]], SHIFTS[k[1]]
    dist = np.hypot(sx, sy) / 1000
    az = np.degrees(np.arctan2(sx, sy)) % 360
    zero = own[len(SHIFTS) // 2, len(SHIFTS) // 2]
    null = []
    for _ in range(2000):
        p = rng.permutation(n)
        while n > 1 and np.any(p == np.arange(n)):
            p = rng.permutation(n)
        null.append(np.max(np.abs(np.median(R[np.arange(n), p], 0))))
    null = np.array(null)
    per_pair = R[np.arange(n), np.arange(n)][:, k[0], k[1]]
    clear_r = [corr_cube(r["rp"], r["mask"], r["dbt"])[len(SHIFTS) // 2, len(SHIFTS) // 2] for r in clear]
    # each cloudy pair's own best shift (exploratory)
    own_best = []
    for i in range(n):
        c = R[i, i]
        kk = np.unravel_index(np.argmax(c), c.shape)
        own_best.append((SHIFTS[kk[1]] / 1000, SHIFTS[kk[0]] / 1000, float(c[kk])))
    out = dict(n_cloudy=n, n_clear=len(clear), best_median_r=float(best_r), shift_east_km=sx / 1000,
               shift_north_km=sy / 1000, shift_km=float(dist), shift_azimuth=float(az),
               r_at_zero_shift=float(zero), null95=float(np.percentile(null, 95)),
               p_value=float(np.mean(null >= abs(best_r))), frac_pairs_positive=float(np.mean(per_pair > 0)),
               clear_median_r=float(np.median(clear_r)) if clear_r else None,
               per_pair=[dict(pair=f"{r['d0']:%Y-%m-%d}", cold=r["cold"], r=float(v))
                         for r, v in zip(cloudy, per_pair)],
               own_best_shift_km=own_best)
    print(f"\n{name}: {n} cloudy pairs, {len(clear)} clear")
    print(f"  common best shift: {sx / 1000:+.0f} km E, {sy / 1000:+.0f} km N  (= {dist:.1f} km toward {az:.0f} deg;"
          f" expected 3-12 km toward {EXPECT_AZ:.0f} +-45)")
    print(f"  median r at best shift {best_r:+.2f} | at zero shift {zero:+.2f} | null 95th pct {out['null95']:.2f}"
          f" | p = {out['p_value']:.3f}")
    print(f"  pairs with r > 0 at best shift: {np.mean(per_pair > 0) * 100:.0f}%  | clear pairs median r (no shift) "
          f"{out['clear_median_r']}")
    for r, v in sorted(zip(cloudy, per_pair), key=lambda t: -t[0]["cold"]):
        print(f"    {r['d0']:%Y-%m-%d}  cold {r['cold'] * 100:4.0f}%  r {v:+.2f}")
    return out, own, cloudy, k


def figure(own, cloudy, k, name, fn):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    top = sorted(range(len(cloudy)), key=lambda i: -cloudy[i]["cold"])[:3]
    fig, ax = plt.subplots(2, 4, figsize=(17, 8))
    s = SHIFTS / 1000
    h = ax[0, 0].imshow(own, origin="lower", extent=(s[0], s[-1], s[0], s[-1]), cmap="RdBu_r", vmin=-0.4, vmax=0.4)
    ax[0, 0].plot([0, 9 * np.sin(np.radians(EXPECT_AZ))], [0, 9 * np.cos(np.radians(EXPECT_AZ))], "k--", lw=1)
    ax[0, 0].plot(s[k[1]], s[k[0]], "k*", ms=12)
    ax[0, 0].set_xlabel("GOES map moved east, km"); ax[0, 0].set_ylabel("moved north, km")
    ax[0, 0].set_title("median r over cloudy pairs\n(dashed = parallax direction)")
    fig.colorbar(h, ax=ax[0, 0], shrink=0.8)
    ax[1, 0].axis("off")
    for j, i in enumerate(top):
        r = cloudy[i]
        a = ax[0, j + 1].imshow(r["rp"] * 100, origin=ORIGIN, cmap="RdBu", vmin=-4, vmax=4)
        ax[0, j + 1].set_title(f"radar residual {r['d0']:%Y-%m-%d}, cm"); fig.colorbar(a, ax=ax[0, j + 1], shrink=0.8)
        b = ax[1, j + 1].imshow(r["dbt"][k[0], k[1]], origin=ORIGIN, cmap="RdBu", vmin=-30, vmax=30)
        ax[1, j + 1].set_title(f"GOES dBT (shifted), K, cold {r['cold'] * 100:.0f}%"); fig.colorbar(b, ax=ax[1, j + 1], shrink=0.8)
    for a in ax.ravel()[1:]:
        a.set_xticks([]); a.set_yticks([])
    fig.suptitle(f"E1 {name}: radar residual vs cloud change (north up)")
    fig.tight_layout()
    fig.savefig(fn, dpi=80)
    plt.close(fig)
    print(fn)


def main():
    rng = np.random.default_rng(7)
    out = {}
    for d, series, name in [("DESCENDING", "mexico-city", "desc 19:13"), ("ASCENDING", "mexico-city-asc", "asc 06:05")]:
        rows = run(d, series)
        o, own, cloudy, k = analyse(rows, name, rng)
        out[name] = o
        figure(own, cloudy, k, name, OUT / f"e1_{name.split()[0]}.png")
    json.dump(out, open(OUT / "e1_clouds_pattern.json", "w"), indent=1)


if __name__ == "__main__":
    main()
