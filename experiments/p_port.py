"""P-port evaluation (rules fixed in notes/experiments.md before any map or GNSS series was looked at).

    uv run python scripts/p_port.py --site central-valley            # whatever passes are ready
P1/P2: NISAR pipeline rate (nisar_rate.tif, chain sum) at GNSS stations, one constant offset removed.
P3:    evening/morning noise = median over pairs of the spatial std of (pair displacement - rate x dt) on 2 km blocks.
P4:    std of (asc - desc) on common pixels after 1 km smoothing.
"""
from __future__ import annotations

import argparse
import csv
import json
from datetime import datetime
from pathlib import Path

import numpy as np
import rasterio
from rasterio.warp import transform as warp_transform
from scipy.ndimage import gaussian_filter

from hachure import gridding

ROOT = Path(__file__).resolve().parents[1]
P = ROOT / "data" / "processed"
THRESH = {"central-valley": dict(P1=2.5, P2=0.9, P4=2.0), "houston": dict(P1=1.0, P2=0.7, P4=1.5)}


VAR = ""
PLANE = True      # pre-registered: "plane removed", as in the Mexico maps (robust one-sided fit on a 400 m-smoothed rate)


def load_rate(site, pas):
    with rasterio.open(P / f"{site}-{pas}{VAR}" / "nisar_rate.tif") as ds:
        v, tr, crs = ds.read(1).astype(float), ds.transform, ds.crs
    if PLANE:
        from hachure.nisar import fit_ramp
        ok = np.isfinite(v)
        wsum = gaussian_filter(ok.astype(float), 400 / 80)
        sm = gaussian_filter(np.where(ok, v, 0.0), 400 / 80) / np.maximum(wsum, 1e-6)
        v = v - fit_ramp(np.where(ok & (wsum > 0.3), sm, np.nan))
    return v * 100, tr, crs


def at_stations(site, pas, rate, transform, crs):
    rows = [r for r in csv.DictReader(open(P / "gnss" / f"{site}_rates.csv")) if r["window"] == pas and r["used"] == "True"]
    xs, ys = warp_transform("EPSG:4326", crs, [float(r["lon"]) for r in rows], [float(r["lat"]) for r in rows])
    H, W = rate.shape
    cols = (np.arange(W) + 0.5) * transform.a + transform.c
    rws = (np.arange(H) + 0.5) * transform.e + transform.f
    out = []
    for r, x, y in zip(rows, xs, ys):
        ci = np.abs(cols - x) <= 250
        ri = np.abs(rws - y) <= 250
        sub = rate[np.ix_(ri, ci)]
        cc, rr = np.meshgrid(cols[ci], rws[ri])
        sub = sub[np.hypot(cc - x, rr - y) <= 250]
        v = np.nanmedian(sub) if np.isfinite(sub).sum() >= 3 else np.nan
        out.append((r["station"], float(r["rate_cm_yr"]), float(r["sigma_cm_yr"]), v, int(np.isfinite(sub).sum())))
    return out


def p12(site, pas):
    rate, tr, crs = load_rate(site, pas)
    st = [s for s in at_stations(site, pas, rate, tr, crs) if np.isfinite(s[3])]
    g = np.array([s[1] for s in st])
    n = np.array([s[3] for s in st])
    off = np.median(n - g)
    res = n - off - g
    rms = float(np.sqrt(np.mean(res ** 2)))
    r = float(np.corrcoef(n, g)[0, 1])
    th = THRESH[site]
    print(f"{site} {pas}: {len(st)} stations with NISAR pixels; offset {off:+.2f} cm/yr")
    for (name, gg, sg, nn, k), e in sorted(zip(st, res), key=lambda a: a[0][1]):
        print(f"   {name}: GNSS {gg:+6.1f} ± {sg:.1f}   NISAR-offset {nn - off:+6.1f}   diff {e:+5.1f}   ({k} px)")
    print(f"   P1 RMS {rms:.2f} cm/yr (<= {th['P1']}) -> {'PASS' if rms <= th['P1'] else 'fail'}; "
          f"P2 corr {r:.3f} (>= {th['P2']}) -> {'PASS' if r >= th['P2'] else 'fail'}; GNSS rate std {g.std():.2f}")
    return dict(n=len(st), offset=float(off), rms=rms, corr=r, gnss_std=float(g.std()),
                stations={s[0]: dict(gnss=s[1], sigma=s[2], nisar=s[3] - off) for s in st})


def p3_noise(site, pas):
    z = np.load(P / f"{site}-{pas}{VAR}" / "nisar_series.npz")
    d = [datetime.fromisoformat(str(x)) for x in z["dates"]]
    rate, *_ = load_rate(site, pas)
    b = 25                                                    # 25 x 80 m = 2 km blocks
    stds = []
    for k in range(len(d) - 1):
        dt = (d[k + 1] - d[k]).days / 365.25
        r = (z["cum"][k + 1] - z["cum"][k]) * 100 - rate * dt
        H, W = (r.shape[0] // b) * b, (r.shape[1] // b) * b
        blk = r[:H, :W].reshape(H // b, b, W // b, b)
        with np.errstate(all="ignore"):
            m = np.nanmean(blk, (1, 3))
            cnt = np.isfinite(blk).sum((1, 3))
        m[cnt < b * b * 0.3] = np.nan
        stds.append(np.nanstd(m))
    return float(np.median(stds)), stds


def smooth(a, sigma):
    ok = np.isfinite(a)
    num = gaussian_filter(np.where(ok, a, 0), sigma)
    den = gaussian_filter(ok.astype(float), sigma)
    with np.errstate(all="ignore"):
        out = num / den
    out[den < 0.3] = np.nan
    return out


def p4(site):
    a, tra, crs = load_rate(site, "asc")
    d, trd, crsd = load_rate(site, "desc")
    d = gridding.to_grid(d, trd, crsd, tra, crs, a.shape)
    sa, sd = smooth(a, 1000 / 80 / 2), smooth(d, 1000 / 80 / 2)   # gaussian sigma 0.5 km = ~1 km smoothing
    diff = sa - sd
    ok = np.isfinite(diff)
    return float(np.nanstd(diff)), int(ok.sum()), float(np.nanmedian(diff))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--site", required=True, choices=list(THRESH))
    ap.add_argument("--variant", default="", help="suffix of the processed dirs, e.g. -ionoplane")
    ap.add_argument("--no-plane", action="store_true", help="transparency only: skip the pre-registered plane removal")
    args = ap.parse_args()
    global PLANE, VAR
    VAR = args.variant
    PLANE = not args.no_plane
    res = {}
    ready = [p for p in ("asc", "desc") if (P / f"{args.site}-{p}{VAR}" / "nisar_rate.tif").exists()]
    for pas in ready:
        res[pas] = p12(args.site, pas)
        res[pas]["noise_cm"], _ = p3_noise(args.site, pas)
        print(f"   pair noise (2 km blocks, median over pairs): {res[pas]['noise_cm']:.2f} cm")
    if len(ready) == 2:
        ratio = res["desc"]["noise_cm"] / res["asc"]["noise_cm"]
        want = ">= 1.5" if args.site == "houston" else "< 1.5"
        ok = ratio >= 1.5 if args.site == "houston" else ratio < 1.5
        print(f"P3 evening/morning noise ratio {ratio:.2f} ({want}) -> {'PASS' if ok else 'fail'}")
        s, n, med = p4(args.site)
        print(f"P4 std(asc - desc), 1 km smoothing: {s:.2f} cm/yr on {n} pixels (median diff {med:+.2f}) "
              f"(<= {THRESH[args.site]['P4']}) -> {'PASS' if s <= THRESH[args.site]['P4'] else 'fail'}")
        res["P3_ratio"], res["P4_std"] = ratio, s
    json.dump(res, open(P / f"p_port_{args.site}{args.variant}{'_noplane' if args.no_plane else ''}.json", "w"), indent=1)


if __name__ == "__main__":
    main()
