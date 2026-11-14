"""Per-date NISAR displacement over Reykjanes for the "Earth breathes" web map,
with per-cell uncertainty and an independent GNSS check.

Processing (GNSS stations SENG, AUSV, GEVK, VOGC are NOT used anywhere here):
  * ascending 159A and descending 082D: cumulative LOS at every date (coherent pairs
    summed; a low-coherence winter pair is bridged with the pixel's mean rate);
  * zero = GNSS station NYLA (median of pixels within 1.5 km of it at every date), the
    same reference as the GNSS series;
  * minus the per-track plane and the north-south motion of an ellipsoid source fitted
    to NISAR ONLY (ellipsoid_fit_plane_insar.npz): both geometries see north with the
    same sign, so it would leak into the vertical;
  * common 12-day dates, vertical + east-west decomposition, 240 m cells;
  * per cell: linear rate and its 1-sigma error from the scatter about the line,
    inflated for lag-1 autocorrelation of the residuals; SNR = |rate| / sigma;
  * hachures that grow with time (one stroke row per 2 cm of total motion), each with
    rate, sigma, SNR and a trust class (strong >= 3 sigma, weak 1-3, below < 1);
  * headline: the ~0.7 km neighbourhood (3 x 3 cells) with the largest rate - 2 sigma.

Output: output/reykjanes/breathing_data.npz and breathing_vectors.json

Usage: python scripts/reykjanes/build_breathing_data.py [--model FILE.npz]
"""
import json
import sys
from datetime import datetime, timedelta

import numpy as np
from scipy.ndimage import gaussian_filter

from reykjanes_analysis import MIN_PAIR_COH, OUT, PLACES, Geo, chain, los_cm, pair_coherence, pairs_by_track
from reykjanes_gnss import REF, STATIONS, read_tenv3, relative

TRACKS = ("159A", "082D")
BLOCK = 3            # 80 m -> 240 m
REF_RADIUS_KM = 1.5
SNR_STRONG, SNR_WEAK = 3.0, 1.0
STEP_CM = 2.0        # one hachure row per 2 cm of motion over the period


def block_median(a, b=BLOCK):
    h, w = (a.shape[-2] // b) * b, (a.shape[-1] // b) * b
    s = a[..., :h, :w].reshape(a.shape[:-2] + (h // b, b, w // b, b))
    with np.errstate(all="ignore"):
        return np.nanmedian(np.nanmedian(s, axis=-1), axis=-2)


def local_km(geo, X, Y):
    """Local east/north km around Svartsengi from EPSG:3413 metres (same convention as reykjanes_mogi)."""
    sx, sy = geo.tr.transform(*PLACES["Svartsengi"])
    lam = np.deg2rad(PLACES["Svartsengi"][0] + 45.0)
    east, north = np.array([np.cos(lam), np.sin(lam)]), np.array([-np.sin(lam), np.cos(lam)])
    dx, dy = (X - sx) / 1000, (Y - sy) / 1000
    return dx * east[0] + dy * east[1], dx * north[0] + dy * north[1]


def track_series(geo, t, nyla_mask):
    prs = chain(pairs_by_track()[t])
    steps, good, dts, los = [], [], [], None
    for t0, t1, p in prs:
        d, g = los_cm(p, geo)
        ok = d is not None and pair_coherence(p, geo) >= MIN_PAIR_COH
        steps.append(d if d is not None else np.full(geo.X.shape, np.nan))
        good.append(ok)
        dts.append((t1 - t0).days)
        los = g if g is not None else los
    steps, good, dts = np.array(steps), np.array(good), np.array(dts, float)
    g = good[:, None, None] & np.isfinite(steps)
    rate = np.where(g, steps, 0).sum(0) / np.maximum((g * dts[:, None, None]).sum(0), 1)  # cm/day
    filled = np.where(g, steps, rate[None] * dts[:, None, None])
    cum = np.concatenate([np.zeros((1,) + geo.X.shape), np.cumsum(filled, 0)])
    cum[:, g.sum(0) / len(prs) < 0.6] = np.nan
    cum -= np.array([np.nanmedian(c[nyla_mask]) for c in cum])[:, None, None]
    dates = [prs[0][0]] + [t1 for _, t1, _ in prs]
    return dates, cum, los, good


def rate_sigma(series, yrs):
    """Least-squares rate per cell and its 1-sigma error (AR(1)-inflated). series: (T, ...)."""
    T = yrs - yrs.mean()
    y = series - np.nanmean(series, 0)
    ok = np.isfinite(y)
    n = ok.sum(0)
    Tm = np.where(ok, T[:, None, None], 0)
    rate = np.nansum(y * Tm, 0) / np.maximum(np.sum(Tm ** 2, 0), 1e-9)
    res = np.where(ok, y - rate[None] * T[:, None, None], np.nan)
    s2 = np.nansum(res ** 2, 0) / np.maximum(n - 2, 1)
    r0, r1 = res[:-1], res[1:]
    rho = np.nansum(r0 * r1, 0) / np.maximum(np.nansum(res ** 2, 0), 1e-9)
    rho = np.clip(rho, 0, 0.95)
    sigma = np.sqrt(s2 / np.maximum(np.sum(Tm ** 2, 0), 1e-9)) * np.sqrt((1 + rho) / (1 - rho))
    bad = n < 0.7 * len(yrs)
    return np.where(bad, np.nan, rate), np.where(bad, np.nan, sigma), rho


def smooth_nan(img, s):
    m = np.isfinite(img)
    num = gaussian_filter(np.where(m, img, 0.0), s)
    den = gaussian_filter(m.astype(float), s)
    return np.where(den > 0.4, num / np.maximum(den, 1e-6), np.nan)


def hachures(x, y, total, cum, rate, sigma, step=STEP_CM, spacing_m=260.0, max_len_m=1300.0, min_contour_m=2500.0):
    """Strokes on contours of total motion, running downslope to the next contour. x, y: 1-D
    metres (ascending); arrays row-major in (y, x). Returns stroke dicts in map metres."""
    from contourpy import contour_generator
    from scipy.interpolate import RegularGridInterpolator
    gy, gx = np.gradient(total, y, x)
    fg = RegularGridInterpolator((y, x), np.dstack([gx, gy]), bounds_error=False, fill_value=np.nan)
    fr = RegularGridInterpolator((y, x), rate, bounds_error=False, fill_value=np.nan)
    fs = RegularGridInterpolator((y, x), sigma, bounds_error=False, fill_value=np.nan)
    fc = [RegularGridInterpolator((y, x), c, bounds_error=False, fill_value=np.nan) for c in cum]
    cg = contour_generator(x=x, y=y, z=np.ma.masked_invalid(total))
    lo, hi = np.nanmin(total), np.nanmax(total)
    levels = sorted([l for l in np.arange(step, hi, step)] + [-l for l in np.arange(step, -lo, step)])
    out = []
    for k, level in enumerate(levels):
        for line in cg.lines(level):
            if len(line) < 3:
                continue
            d = np.r_[0, np.cumsum(np.hypot(*np.diff(line, axis=0).T))]
            if d[-1] < min_contour_m:
                continue
            s = np.arange((k % 2) * spacing_m / 2, d[-1], spacing_m)
            pts = np.c_[np.interp(s, d, line[:, 0]), np.interp(s, d, line[:, 1])]
            g = fg(pts[:, ::-1])
            slope = np.hypot(g[:, 0], g[:, 1])
            ok = np.isfinite(slope) & (slope > step / max_len_m)
            for p0, gg, sl in zip(pts[ok], g[ok], slope[ok]):
                length = step / sl * 0.85
                p1 = p0 - gg / sl * length                      # downslope = towards less motion
                q = p0[::-1][None]
                r, sg = float(fr(q)[0]), float(fs(q)[0])
                snr = abs(r) / sg if np.isfinite(sg) and sg > 0 else 0.0
                out.append(dict(p0=p0, p1=p1, level=float(level), rate=r, sigma=sg, snr=snr,
                                slope=float(sl * 1000), t=[float(f(q)[0]) for f in fc]))
    return out


def main(model_file="ellipsoid_fit_plane_insar.npz"):
    geo = Geo(pairs_by_track()["159A"][0][2])
    nyla_ll = read_tenv3(REF)[2]
    nx_, ny_ = geo.tr.transform(*nyla_ll)
    nyla_mask = np.hypot(geo.X - nx_, geo.Y - ny_) < REF_RADIUS_KM * 1000
    E, N = local_km(geo, geo.X, geo.Y)
    ne, nn = local_km(geo, nx_, ny_)

    # corrections from a source model fitted to NISAR only
    from reykjanes_ellipsoid import cdm
    fit = np.load(OUT / model_file)
    m = fit["m"]
    U = m[:3]
    G = cdm(E.ravel(), N.ravel(), tuple(fit["p"]))
    model = (G * U).sum(-1).reshape(E.shape + (3,))
    G0 = (cdm(np.array([ne]), np.array([nn]), tuple(fit["p"]))[0] * U).sum(-1)
    model_up, model_north = model[..., 2] - G0[2], model[..., 1] - G0[1]
    has_plane = m.size >= 3 + 3 * len(TRACKS)
    print(f"corrections from {model_file} (GNSS used in that fit: {'insar' not in model_file})")

    raw = {}
    for k, t in enumerate(TRACKS):
        dates, cum, los, good = track_series(geo, t, nyla_mask)
        yrs = np.array([(d - dates[0]).days / 365.25 for d in dates])
        if has_plane:
            pe, pn = m[4 + 3 * k], m[5 + 3 * k]
            cum = cum - yrs[:, None, None] * ((pe * (E - ne) + pn * (N - nn)) / 10)[None]
        cum = cum - yrs[:, None, None] * (los["losUnitVectorY"] * model_north)[None]
        raw[t] = dict(dates=dates, cum=cum, los=los, good=good)
        print(f"{t}: {len(dates) - 1} pairs, {good.sum()} coherent")

    t0 = max(raw[t]["dates"][0] for t in TRACKS)
    t1 = min(raw[t]["dates"][-1] for t in TRACKS)
    grid, d = [], t0
    while d <= t1:
        grid.append(d)
        d += timedelta(days=12)
    if grid[-1] != t1:
        grid.append(t1)

    def interp(t):
        ds = np.array([(x - raw[t]["dates"][0]).days for x in raw[t]["dates"]], float)
        q = np.array([(x - raw[t]["dates"][0]).days for x in grid], float)
        c = raw[t]["cum"]
        i = np.clip(np.searchsorted(ds, q) - 1, 0, len(ds) - 2)
        w = ((q - ds[i]) / (ds[i + 1] - ds[i]))[:, None, None]
        out = c[i] * (1 - w) + c[i + 1] * w
        return out - out[0]

    L = {t: interp(t) for t in TRACKS}
    e = {t: raw[t]["los"]["losUnitVectorX"] for t in TRACKS}
    u = {t: np.cos(np.deg2rad(raw[t]["los"]["incidenceAngle"])) for t in TRACKS}
    det = e["159A"] * u["082D"] - e["082D"] * u["159A"]
    up = block_median((e["159A"] * L["082D"] - e["082D"] * L["159A"]) / det)
    east = block_median((L["159A"] * u["082D"] - L["082D"] * u["159A"]) / det)
    Xd, Yd = block_median(geo.X), block_median(geo.Y)
    model_up_d = block_median(model_up)
    yrs = np.array([(x - grid[0]).days / 365.25 for x in grid])

    # per-cell rate and sigma
    rate, sigma, rho = rate_sigma(up, yrs)
    snr = np.abs(rate) / sigma
    trusted = snr >= SNR_STRONG
    print(f"cells: {np.isfinite(rate).sum()} with data; strong {np.nansum(snr >= SNR_STRONG)}, "
          f"weak {np.nansum((snr >= SNR_WEAK) & (snr < SNR_STRONG))}, below {np.nansum(snr < SNR_WEAK)}; "
          f"median sigma {np.nanmedian(sigma):.1f} cm/yr, median lag-1 rho {np.nanmedian(rho):.2f}")

    # independent GNSS check, by rates
    val = {}
    for s in STATIONS:
        tt, x, ll = relative(s, REF)
        sx, sy = geo.tr.transform(*ll)
        i = int(np.argmin(np.abs(Yd[:, 0] - sy)))
        j = int(np.argmin(np.abs(Xd[0] - sx)))
        w = (slice(max(i - 1, 0), i + 2), slice(max(j - 1, 0), j + 2))
        nis = np.array([np.nanmedian(a[w]) for a in up])
        nr, ns, _ = rate_sigma(nis[:, None, None], yrs)
        msk = np.array([(q >= grid[0]) & (q <= grid[-1]) for q in tt])
        gy_ = np.array([(q - grid[0]).days / 365.25 for q in tt[msk]])
        gu = x[msk, 2]
        gr = np.polyfit(gy_, gu, 1)[0]
        gres = gu - np.polyval(np.polyfit(gy_, gu, 1), gy_)
        gs = np.std(gres) / np.sqrt(np.sum((gy_ - gy_.mean()) ** 2)) * np.sqrt(30)  # daily noise is correlated: ~30-day blocks
        nr, ns = float(nr[0, 0]), float(ns[0, 0])
        z = (nr - gr) / np.hypot(ns, gs)
        val[s] = dict(nisar=round(nr, 1), nisar_sigma=round(ns, 1), gnss=round(float(gr), 1), gnss_sigma=round(float(gs), 2),
                      diff=round(nr - float(gr), 1), z=round(float(z), 1))
        print(f"independent check {s}: NISAR {nr:+.1f} +- {ns:.1f} vs GNSS {gr:+.1f} +- {gs:.2f} cm/yr; "
              f"diff {nr - gr:+.1f} ({z:+.1f} sigma)")

    # headline: 3 x 3 neighbourhood with the largest rate - 2 sigma (sigma not reduced: errors are correlated)
    from scipy.ndimage import uniform_filter
    okc = np.isfinite(rate) & np.isfinite(sigma)
    cnt = uniform_filter(okc.astype(float), 3, mode="constant") * 9
    rmean = uniform_filter(np.where(okc, rate, 0), 3, mode="constant") * 9 / np.maximum(cnt, 1)
    smean = uniform_filter(np.where(okc, sigma, 0), 3, mode="constant") * 9 / np.maximum(cnt, 1)
    score = np.where(cnt >= 8, rmean - 2 * smean, -np.inf)
    hi, hj = np.unravel_index(np.argmax(score), score.shape)
    hlon, hlat = geo.tr.transform(Xd[hi, hj], Yd[hi, hj], direction="INVERSE")
    p99 = float(np.nanpercentile(rate[trusted], 99))
    headline = dict(rate=round(float(rmean[hi, hj]), 1), sigma=round(float(smean[hi, hj]), 1),
                    at_least=round(float(rmean[hi, hj] - 2 * smean[hi, hj]), 1), lon=float(hlon), lat=float(hlat),
                    cum=round(float(np.nanmedian(up[-1][max(hi - 1, 0):hi + 2, max(hj - 1, 0):hj + 2])), 1),
                    p99=round(p99, 1), area_km2=round(float(np.sum(trusted & (rate > 5)) * 0.24 ** 2), 1))
    print("headline:", headline)

    # time-growing hachures with trust
    xs, ys = Xd[0], Yd[:, 0]
    flip = ys[1] < ys[0]
    sl = slice(None, None, -1) if flip else slice(None)
    total = smooth_nan(up[-1], 1.5)[sl]
    cum_s = np.array([smooth_nan(c, 1.5) for c in up])[:, sl]
    strokes = hachures(xs, ys[sl], total, cum_s, smooth_nan(rate, 1.5)[sl], smooth_nan(sigma, 1.5)[sl])
    vec = []
    for st in strokes:
        a = geo.tr.transform(*st["p0"], direction="INVERSE")
        b = geo.tr.transform(*st["p1"], direction="INVERSE")
        trust = "strong" if st["snr"] >= SNR_STRONG else "weak" if st["snr"] >= SNR_WEAK else "below"
        vec.append(dict(a=[round(a[0], 6), round(a[1], 6)], b=[round(b[0], 6), round(b[1], 6)], level=st["level"],
                        rate=round(st["rate"], 1), sigma=round(st["sigma"], 1), snr=round(st["snr"], 1), trust=trust,
                        slope=round(st["slope"], 1), t=[round(v, 1) if np.isfinite(v) else None for v in st["t"]]))
    tc = {k: sum(v["trust"] == k for v in vec) for k in ("strong", "weak", "below")}
    print(f"hachures: {len(vec)} strokes, trust {tc}")

    lon, lat = geo.tr.transform(Xd, Yd, direction="INVERSE")
    out = dict(lon=lon.astype(np.float32), lat=lat.astype(np.float32),
               dates=np.array([x.strftime("%Y-%m-%d") for x in grid]),
               up=up.astype(np.float32), east=east.astype(np.float32), model_up=model_up_d.astype(np.float32),
               rate=rate.astype(np.float32), sigma=sigma.astype(np.float32))
    gnss = {}
    for s in STATIONS:
        tt, x, ll = relative(s, REF)
        msk = np.array([(q >= grid[0] - timedelta(days=3)) & (q <= grid[-1] + timedelta(days=3)) for q in tt])
        x0 = np.median(x[msk][:5], 0)
        gnss[s] = dict(t=[q.strftime("%Y-%m-%d") for q in tt[msk]][::2], up=np.round(x[msk][::2, 2] - x0[2], 2).tolist(),
                       lon=float(ll[0]), lat=float(ll[1]), val=val[s])
    np.savez_compressed(OUT / "breathing_data.npz", **out)
    (OUT / "breathing_vectors.json").write_text(json.dumps(dict(
        hachures=vec, trust_counts=tc, headline=headline, gnss=gnss, model_file=model_file,
        snr=dict(strong=SNR_STRONG, weak=SNR_WEAK), step_cm=STEP_CM,
        sections=[dict(name="Північ–південь: Гріндавік → Свартсенгі → Вогар", a=[-22.445, 63.825], b=[-22.395, 63.965]),
                  dict(name="Уздовж осі резервуару (NE–SW)", a=[-22.60, 63.835], b=[-22.30, 63.915])])), encoding="utf-8")
    print("saved breathing_data.npz", up.shape, "and breathing_vectors.json")


if __name__ == "__main__":
    mf = sys.argv[sys.argv.index("--model") + 1] if "--model" in sys.argv else "ellipsoid_fit_plane_insar.npz"
    main(mf)
