"""H6: is the morning-minus-evening 'vertical' rate difference east-west motion toward the bowls?
H3b: is the MMX1 bias local (rate changes fast around the station)?   (both pre-registered)

    uv run python scripts/h6_horizontal.py
"""
from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path

import numpy as np
from rasterio.warp import transform as warp_transform
from scipy.ndimage import gaussian_filter

from hachure import gridding, nisar

ROOT = Path(__file__).resolve().parents[1]
P = ROOT / "data" / "processed"
MMX1 = (-99.0684, 19.4317)
GNSS_UP = -26.7


def rate(name, dst=None):
    z = np.load(P / name / "nisar_series.npz")
    d = [datetime.fromisoformat(str(x)) for x in z["dates"]]
    t = np.array([(x - d[0]).days for x in d]) / 365.25
    cum = z["cum"].astype("float64")
    if dst is not None:
        src = gridding.transform_from_centres(z["x"], z["y"])
        cum = np.stack([gridding.to_grid(c, src, f"EPSG:{int(z['epsg'])}", *dst) for c in cum])
    w = np.isfinite(cum).astype(float)
    c = np.nan_to_num(cum)
    T = t[:, None, None]
    sw, st, sy, stt, sty = w.sum(0), (w * T).sum(0), (w * c).sum(0), (w * T * T).sum(0), (w * T * c).sum(0)
    with np.errstate(all="ignore"):
        v = (sw * sty - st * sy) / (sw * stt - st ** 2)
    v[sw < 10] = np.nan
    return v * 100, z


def incidence(direction):
    BBOX = (-99.32, 19.17, -98.82, 19.58)
    found = nisar.search((BBOX[0] + BBOX[2]) / 2, (BBOX[1] + BBOX[3]) / 2, direction=direction)
    g = nisar.sequential_chain(nisar.dominant_track(found))[0]
    cache = ROOT / "data" / "raw" / "nisar_cache" / ("bbox_" + "_".join(f"{v:.3f}" for v in BBOX) + "_h2240")
    r = nisar.RawWindow.load(cache / f"{g.name}.npz")
    return float(np.nanmedian(r.incidence))


def gnss_en():
    t, e, n = [], [], []
    for line in open(ROOT / "data/raw/gnss/MMX1.tenv3"):
        p = line.split()
        if not p or p[0] != "MMX1":
            continue
        y = float(p[2])
        if 2025.81 <= y <= 2026.81:
            t.append(y); e.append(float(p[7]) + float(p[8])); n.append(float(p[9]) + float(p[10]))
    t = np.array(t)
    return np.polyfit(t, e, 1)[0] * 100, np.polyfit(t, n, 1)[0] * 100, len(t)


def smooth(a, sigma):
    ok = np.isfinite(a)
    num = gaussian_filter(np.where(ok, a, 0), sigma)
    den = gaussian_filter(ok.astype(float), sigma)
    with np.errstate(all="ignore"):
        out = num / den
    out[den < 0.3] = np.nan
    return out


def main():
    vd, z = rate("mexico-city")
    dst = (gridding.transform_from_centres(z["x"], z["y"]), f"EPSG:{int(z['epsg'])}", vd.shape)
    va, _ = rate("mexico-city-asc", dst)
    th_a, th_d = incidence("ASCENDING"), incidence("DESCENDING")
    th = np.radians((th_a + th_d) / 2)
    print(f"incidence: asc {th_a:.1f}, desc {th_d:.1f} deg")
    dx = abs(z["x"][1] - z["x"][0])
    sig = 1000 / dx / 2.355                                   # 1 km FWHM
    D = smooth(va - vd, sig)
    D -= np.nanmedian(D)
    vu = smooth((va + vd) / 2, sig)
    xsign = np.sign(z["x"][1] - z["x"][0])
    ysign = np.sign(z["y"][1] - z["y"][0])
    gy, gx = np.gradient(vu, dx)                               # per m along rows/cols
    dvdx, dvdy = gx * xsign * 1000, gy * ysign * 1000          # cm/yr per km, east / north
    m = np.isfinite(D) & np.isfinite(dvdx) & np.isfinite(dvdy)
    cx = np.corrcoef(D[m], dvdx[m])[0, 1]
    cy = np.corrcoef(D[m], dvdy[m])[0, 1]
    # sub-sample for an honest effective-n note: 1 km blocks
    print(f"H6: corr(D, dv/dx) = {cx:+.2f}   corr(D, dv/dy) = {cy:+.2f}   (pixels {m.sum()})")
    E = D / (2 * np.tan(th))
    k = np.polyfit(dvdx[m], E[m], 1)
    print(f"    E = {k[0]:+.2f} * dv/dx {k[1]:+.2f}   (E cm/yr, dv/dx cm/yr per km)")

    xx, yy = np.meshgrid(z["x"], z["y"])
    sx, sy = warp_transform("EPSG:4326", f"EPSG:{int(z['epsg'])}", [MMX1[0]], [MMX1[1]])
    r = np.hypot(xx - sx[0], yy - sy[0])
    ge, gn, ng = gnss_en()
    near = r <= 500
    print(f"    MMX1 GNSS: east {ge:+.1f}, north {gn:+.1f} cm/yr ({ng} days);  radar-implied east "
          f"{np.nanmedian(E[near]):+.1f} cm/yr (D at MMX1 {np.nanmedian(D[near]):+.1f})")
    print(f"    note: D is median-referenced, so radar-implied E is relative to the city median")

    print("H3b: rate around MMX1 (cm/yr), GNSS up", GNSS_UP)
    out = dict(inc_asc=th_a, inc_desc=th_d, corr_dx=float(cx), corr_dy=float(cy), slope=k.tolist(),
               gnss_e=ge, gnss_n=gn, radar_E_mmx1=float(np.nanmedian(E[near])), rings={})
    for lo, hi in [(0, 60), (0, 250), (0, 500), (500, 1000), (1000, 2000)]:
        mm = (r >= lo) & (r < hi)
        a, d = np.nanmedian(va[mm]), np.nanmedian(vd[mm])
        out["rings"][f"{lo}-{hi} m"] = [float(a), float(d)]
        print(f"    {lo:5d}-{hi:5d} m  asc {a:+6.1f}  desc {d:+6.1f}  (pixels {np.sum(mm & np.isfinite(va))})")
    near1 = r <= 1000
    grad = np.hypot(dvdx, dvdy)
    print(f"    local gradient within 1 km: median {np.nanmedian(grad[near1]):.1f} cm/yr per km")
    out["grad_mmx1"] = float(np.nanmedian(grad[near1]))
    json.dump(out, open(P / "mexico-city" / "h6_horizontal.json", "w"), indent=1)
    figure(D, dvdx, E, ge, sx[0], sy[0], z, cx)


def figure(D, dvdx, E, ge, sx, sy, z, cx):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ext = (z["x"][0], z["x"][-1], z["y"][-1], z["y"][0])
    fig, ax = plt.subplots(1, 2, figsize=(14, 6))
    lim = np.nanpercentile(np.abs(D), 98)
    h = ax[0].imshow(D, cmap="RdBu", vmin=-lim, vmax=lim, extent=ext)
    ax[0].set_title("morning − evening 'vertical' rate, cm/yr")
    fig.colorbar(h, ax=ax[0], shrink=0.8)
    lim2 = np.nanpercentile(np.abs(dvdx), 98)
    h = ax[1].imshow(-dvdx, cmap="RdBu", vmin=-lim2, vmax=lim2, extent=ext)
    ax[1].set_title(f"−∂(vertical rate)/∂x, cm/yr per km  (corr with left: {-cx:+.2f})")
    fig.colorbar(h, ax=ax[1], shrink=0.8)
    for a in ax:
        a.plot(sx, sy, "k^", ms=8)
        a.set_xticks([]); a.set_yticks([])
    fig.tight_layout()
    fig.savefig(P / "mexico-city" / "h6_horizontal.png", dpi=85)
    plt.close(fig)


if __name__ == "__main__":
    main()
