"""Site engine: config -> style-free data layers for any front-end.

A *site* is a JSON file in `sites/` (bbox, pass, reference area, processing and hachure parameters).
Everything here produces DATA, never colours: hachure geometry with physical properties
(rate, 1-sigma, SNR, chain breaks, per-vertex time), a coarse grid with time series, a manifest.
Styling lives in the front-end themes (`web/themes/*.json`), so a new place or a new look
never needs engine changes.

    from hachure import engine
    cfg = engine.load_site("sites/mexico-city-asc.json")
    v = engine.velocity(cfg)                 # rate, sigma, segments, cumulative series on the radar grid
    engine.export_site(cfg, v, out_dir)      # GeoJSON / JSON layers + manifest.json
"""
from __future__ import annotations

import json
from collections import Counter
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

import numpy as np
from rasterio.warp import transform as warp_transform
from scipy.ndimage import gaussian_filter, map_coordinates, median_filter

from . import dem, field, gridding, nisar, trace, trust

ROOT = Path(__file__).resolve().parents[2]
ENGINE_VERSION = 1


# ------------------------------------------------------------------ config

def load_site(path) -> dict:
    cfg = json.loads(Path(path).read_text(encoding="utf-8"))
    cfg.setdefault("layers", ["terrain", "velocity", "grid"])
    cfg.setdefault("trust", {"snr_strong": 3.0, "snr_weak": 1.0, "segments": True})
    cfg.setdefault("grid_km", 0.5)
    cfg.setdefault("gnss", [])
    cfg.setdefault("anchor", "radar-only")
    return cfg


def cache_dir(cfg) -> Path:
    b = cfg["bbox"]
    return ROOT / "data/raw/nisar_cache" / ("bbox_" + "_".join(f"{v:.3f}" for v in b) + f"_h{int(cfg['cache_height'])}")


def cached_chain(cfg) -> list[nisar.Granule]:
    """Sequential chain of the dominant track among CACHED windows (works offline)."""
    letter = cfg["direction"][0].upper()
    grans = []
    for f in sorted(cache_dir(cfg).glob("*.npz")):
        parts = f.stem.split("_")
        if len(parts) < 8 or parts[6] != letter:
            continue
        ref, sec = nisar._parse_dates(f.stem)
        grans.append(nisar.Granule(f.stem, "", cfg["direction"].upper(), int(parts[5]), int(parts[7]), ref, sec, None))
    if not grans:
        raise FileNotFoundError(f"no cached {cfg['direction']} windows in {cache_dir(cfg)}")
    top = Counter((g.path, g.frame) for g in grans).most_common(1)[0][0]
    return nisar.sequential_chain([g for g in grans if (g.path, g.frame) == top])


# ------------------------------------------------------------------ radar products

@dataclass
class Velocity:
    rate: np.ndarray        # m/yr, vertical, chain sum
    sigma: np.ndarray       # m/yr, 1-sigma (trust audit, optional chain-break term)
    segments: np.ndarray    # runs of valid pairs per pixel (1 = unbroken, 0 = no data)
    coherence: np.ndarray   # mean coherence
    cum: np.ndarray         # (n_dates, H, W) cumulative vertical displacement, m
    dates: list[datetime]
    x: np.ndarray
    y: np.ndarray
    epsg: int
    n_pairs: int


def velocity(cfg, verbose: bool = True) -> Velocity:
    chain = cached_chain(cfg)
    iono = {"full": True, "none": False}.get(cfg.get("iono", "plane"), cfg.get("iono", "plane"))
    pairs, kept = [], []
    for g in chain:
        try:
            pairs.append(nisar.process(nisar.RawWindow.load(cache_dir(cfg) / f"{g.name}.npz"),
                                       ref_bbox=tuple(cfg["ref_bbox"]), coh_min=cfg.get("coh_min", 0.3), iono=iono))
            kept.append(g)
        except ValueError as e:             # unusable product: becomes a gap in the chain
            if verbose:
                print(f"  skip {g.name[:48]}...: {e}")
    rate, coh, frac, x, y, disps = nisar.stack_rate(pairs)
    days = [p.days for p in pairs]
    dates = [kept[0].ref_date] + [g.sec_date for g in kept]
    a = trust.audit(pairs, [d.isoformat() for d in dates], segments=cfg["trust"].get("segments", True))
    ok = np.isfinite(disps)
    seg = ok[0].astype(int) + (ok[1:] & ~ok[:-1]).sum(0)
    seg = np.where(ok.any(0), seg, 0)
    cum = nisar.cumulative(disps, days, rate)
    if verbose:
        print(f"  {len(pairs)} pairs, {dates[0]:%Y-%m-%d} -> {dates[-1]:%Y-%m-%d}, "
              f"median sigma {np.nanmedian(a.sigma):.2f} cm/yr, unbroken pixels {np.mean(seg[seg > 0] == 1) * 100:.0f}%")
    return Velocity(rate=rate, sigma=a.sigma / 100, segments=seg, coherence=coh, cum=cum, dates=dates,
                    x=x, y=y, epsg=pairs[0].epsg, n_pairs=len(pairs))


def save_velocity(v: Velocity, path) -> None:
    np.savez_compressed(path, rate=v.rate, sigma=v.sigma, segments=v.segments, coherence=v.coherence,
                        cum=v.cum.astype("float32"), dates=np.array([d.isoformat() for d in v.dates]),
                        x=v.x, y=v.y, epsg=v.epsg, n_pairs=v.n_pairs)


def load_velocity(path) -> Velocity:
    z = np.load(path)
    return Velocity(rate=z["rate"], sigma=z["sigma"], segments=z["segments"], coherence=z["coherence"],
                    cum=z["cum"].astype("float64"), dates=[datetime.fromisoformat(str(d)) for d in z["dates"]],
                    x=z["x"], y=z["y"], epsg=int(z["epsg"]), n_pairs=int(z["n_pairs"]))


def _height_term(d, xx, yy, h, use, clips=2):
    """Fit d = a + b x + c y + e h on `use`, with 3-sigma clips; return e (per metre)."""
    m = use & np.isfinite(d)
    e = 0.0
    for _ in range(clips + 1):
        if m.sum() < 200:
            return 0.0
        A = np.column_stack([np.ones(m.sum()), xx[m], yy[m], h[m]])
        c = np.linalg.lstsq(A, d[m], rcond=None)[0]
        e = c[3]
        res = d - (c[0] + c[1] * xx + c[2] * yy + c[3] * h)
        s = np.nanstd(res[m])
        m = m & (np.abs(res) < 3 * s)
    return float(e)


def tropo_correct(cfg, v: Velocity, verbose: bool = True) -> Velocity:
    """Remove a height-correlated delay (stratified troposphere): cfg['tropo'] = {'mode': 'height', 'stable_above_m': ...}.

    Fitted per date on coherent pixels above `stable_above_m` (hills that should not move), jointly with a plane;
    only the height term is removed. h0 = median height in ref_bbox, so the reference area stays at zero.
    """
    p = cfg.get("tropo") or {}
    if p.get("mode") != "height":
        return v
    z, T, crs = _dem(cfg)
    h = gridding.to_grid(z, T, crs, gridding.transform_from_centres(v.x, v.y), f"EPSG:{v.epsg}", v.rate.shape)
    lon0, lat0, lon1, lat1 = cfg["ref_bbox"]
    rx, ry = warp_transform("EPSG:4326", f"EPSG:{v.epsg}", [lon0, lon1], [lat0, lat1])
    inref = ((v.x[None, :] >= min(rx)) & (v.x[None, :] <= max(rx)) & (v.y[:, None] >= min(ry)) & (v.y[:, None] <= max(ry)))
    h0 = float(np.nanmedian(h[inref])) if np.isfinite(h[inref]).any() else float(np.nanmedian(h))
    dh = h - h0
    yy, xx = np.mgrid[:h.shape[0], :h.shape[1]].astype(float)
    use = (h >= p.get("stable_above_m", 2400)) & (v.coherence >= cfg["velocity"].get("coh_min_draw", 0.45))
    cum = v.cum.copy()
    es = [0.0]
    for k in range(1, len(cum)):
        e = _height_term(cum[k] - cum[0], xx, yy, dh, use)
        cum[k] = cum[k] - e * np.nan_to_num(dh)
        es.append(e)
    e_rate = _height_term(v.rate, xx, yy, dh, use)
    rate = v.rate - e_rate * np.nan_to_num(dh)
    if verbose:
        print(f"  tropo: h0 {h0:.0f} m, stable px {int(use.sum())}, rate term {e_rate * 1e4:+.2f} cm/yr per 100 m, "
              f"last-date term {es[-1] * 1e4:+.2f} cm per 100 m")
    out = Velocity(rate=rate, sigma=v.sigma, segments=v.segments, coherence=v.coherence, cum=cum, dates=v.dates,
                   x=v.x, y=v.y, epsg=v.epsg, n_pairs=v.n_pairs)
    out.tropo = {"h0_m": round(h0, 1), "rate_cm_yr_per_100m": round(e_rate * 1e4, 3),
                 "stable_above_m": p.get("stable_above_m", 2400)}
    return out


def highpass_correct(cfg, v: Velocity, verbose: bool = True) -> Velocity:
    """Remove long-wavelength residuals (atmosphere, ionosphere left after the plane): cfg['highpass'] = {sigma_km, mask_cm_yr, grow_m, detect_m}.

    A Gaussian low-pass (sigma_km) of the rate and of each date's displacement is subtracted. Small fast features must not
    leak into the low-pass, so pixels whose rate departs from a first, unmasked low-pass by more than mask_cm_yr (grown by
    grow_m) are left out of it; broad patches are not masked and are removed. Off unless the site asks for it.
    """
    p = cfg.get("highpass") or {}
    if not p:
        return v
    from scipy.ndimage import binary_dilation
    px = abs(float(v.x[1] - v.x[0]))
    sig = p.get("sigma_km", 5) * 1000 / px

    def lowpass(a, skip, s=sig):
        w = np.isfinite(a) & ~skip
        num = gaussian_filter(np.where(w, a, 0.0), s, mode="nearest")
        den = gaussian_filter(w.astype(float), s, mode="nearest")
        return np.where(den > 0.05, num / np.maximum(den, 1e-6), 0.0)

    none = np.zeros(v.rate.shape, bool)
    # the per-pixel rate is noisy, so fast features are found on a lightly smoothed rate (detect_m)
    resid = (lowpass(v.rate, none, p.get("detect_m", 400) / px) - lowpass(v.rate, none)) * 100
    moving = binary_dilation(np.abs(np.nan_to_num(resid)) > p.get("mask_cm_yr", 8),
                             iterations=max(1, int(round(p.get("grow_m", 500) / px))))
    lp_rate = lowpass(v.rate, moving)
    cum = v.cum.copy()
    for k in range(1, len(cum)):
        cum[k] = cum[k] - lowpass(cum[k] - cum[0], moving)
    if verbose:
        r = lp_rate[np.isfinite(v.rate)] * 100
        print(f"  highpass: sigma {p.get('sigma_km', 5)} km, masked {moving.mean():.1%} of pixels, "
              f"removed rate p5/p50/p95 {np.percentile(r, 5):+.1f} / {np.percentile(r, 50):+.1f} / {np.percentile(r, 95):+.1f} cm/yr")
    out = Velocity(rate=v.rate - lp_rate, sigma=v.sigma, segments=v.segments, coherence=v.coherence, cum=cum, dates=v.dates,
                   x=v.x, y=v.y, epsg=v.epsg, n_pairs=v.n_pairs)
    if getattr(v, "tropo", None):
        out.tropo = v.tropo
    out.highpass = {"sigma_km": p.get("sigma_km", 5), "mask_cm_yr": p.get("mask_cm_yr", 8), "grow_m": p.get("grow_m", 500),
                    "detect_m": p.get("detect_m", 400)}
    return out


# ------------------------------------------------------------------ helpers

def _smooth_nan(a, sig, valid):
    w = gaussian_filter(valid.astype(float), sig)
    return gaussian_filter(np.where(valid, a, 0.0), sig) / np.maximum(w, 1e-6), w


def _at(a, pts):
    return map_coordinates(a, [pts[:, 1], pts[:, 0]], order=1, mode="nearest")


def _lonlat(pts, T, crs):
    xs = T.c + (pts[:, 0] + 0.5) * T.a
    ys = T.f + (pts[:, 1] + 0.5) * T.e
    lon, lat = warp_transform(crs, "EPSG:4326", xs.tolist(), ys.tolist())
    return np.round(np.column_stack([lon, lat]), 5)


def _fc(feats):
    return {"type": "FeatureCollection", "features": feats}


# ------------------------------------------------------------------ layers

def terrain_layer(cfg, z, T, crs, simplify=2):
    p = cfg.get("terrain", {})
    f = field.prepare(z, cfg["dem_res"], sigma_dir_m=p.get("sigma_dir_m", 150), sigma_tone_m=p.get("sigma_tone_m", 80))
    strokes = trace.hachures(f, interval_m=p.get("interval_m", 60), spacing_px=p.get("spacing_px", 1.5))
    feats = []
    for s in strokes:
        q = s.pts[::simplify] if len(s.pts) > 2 * simplify else s.pts
        feats.append({"type": "Feature", "properties": {"dark": round(s.dark, 3)},
                      "geometry": {"type": "LineString", "coordinates": _lonlat(q, T, crs).tolist()}})
    return _fc(feats)


def velocity_layer(cfg, v: Velocity, z, T, crs, simplify=2):
    """Hachures of the final cumulative displacement; every stroke carries its physics and its timing.

    properties: rate_cm_yr, sigma_cm_yr, snr, segments, trust ('strong' | 'weak' | 'below'), dark (0..1 tone),
                t (days since the first date at which each vertex is reached, for growth animation;
                   null = never within the series).
    """
    p = cfg["velocity"]
    res = cfg["dem_res"]
    shape = z.shape
    src_t = gridding.transform_from_centres(v.x, v.y)
    src_crs = f"EPSG:{v.epsg}"

    def grid(a):
        return gridding.to_grid(np.asarray(a, float), src_t, src_crs, T, crs, shape)

    years = (v.dates[-1] - v.dates[0]).days / 365.25
    coherent = np.nan_to_num(grid(v.coherence)) >= p.get("coh_min_draw", 0.45)
    sig = 400 / res
    frames_s, final = [], None
    for k in range(len(v.dates)):
        g = grid(v.cum[k] - v.cum[0])
        ok = coherent & np.isfinite(g)
        g = median_filter(np.where(ok, g, np.nanmedian(g[ok]) if ok.any() else 0), size=5)
        s, w = _smooth_nan(g, sig, ok)
        if k > 0 and p.get("deramp", True):
            yrs = (v.dates[k] - v.dates[0]).days / 365.25
            plane = nisar.fit_ramp(np.where(ok & (w > 0.3), s, np.nan), cut=-max(0.005, 0.03 * yrs))
            g, s = g - plane, s - plane
        frames_s.append(s)
        final = np.where(ok, g, np.nan)
    final_s = frames_s[-1]
    support = _smooth_nan(np.ones(shape), sig, coherent)[1] > 0.3
    floor = p.get("noise_floor_m_yr", 0.05) * years
    drawn = (final_s <= -floor) if p.get("sink_only", True) else (np.abs(final_s) >= floor)
    valid = coherent & drawn & support
    fsub, interval = field.prepare_scalar(np.where(np.isfinite(final), final, 0.0), res, band=p.get("band_m", 0.02),
                                          valid=valid, sigma_dir_m=p.get("sigma_dir_m", 500),
                                          sigma_tone_m=p.get("sigma_tone_m", 250))
    strokes = trace.hachures(fsub, interval_m=interval, spacing_px=1.5, max_len_px=p.get("max_len_px", 30),
                             min_slope_deg=1.0, min_dark=0.12)

    rate_g = final_s / years * 100                                     # cm/yr, smoothed, deramped
    sigma_g = grid(v.sigma * 100)
    sigma_g = _smooth_nan(np.nan_to_num(sigma_g), 250 / res, np.isfinite(sigma_g))[0]
    seg_g = grid(np.where(v.segments > 0, v.segments, np.nan))
    days = np.array([(d - v.dates[0]).days for d in v.dates], float)
    mid = np.array([s.pts[len(s.pts) // 2] for s in strokes]) if strokes else np.zeros((0, 2))
    mids_t = np.array([_at(np.nan_to_num(f), mid) for f in frames_s])          # (n_dates, n_strokes)
    final_mid = np.minimum(mids_t[-1], -1e-4) if p.get("sink_only", True) else np.where(
        np.abs(mids_t[-1]) < 1e-4, 1e-4, mids_t[-1])
    frac = np.maximum.accumulate(np.clip(mids_t / final_mid, 0, 1), axis=0)
    frac[-1] = 1.0
    ts = cfg["trust"]
    feats = []
    for i, s in enumerate(strokes):
        q = s.pts[::simplify] if len(s.pts) > 2 * simplify else s.pts
        # orient every stroke from the high end (less sinking) to the low end, so a front-end can taper it
        # (thick at the top, like Swiss wedge hachures): a bowl then shows converging wedges, a dome diverging ones
        ends = _at(np.nan_to_num(final_s), q[[0, -1]])
        if ends[0] < ends[-1]:
            q = q[::-1]
        sense = "down" if _at(np.nan_to_num(final_s), mid[i:i + 1])[0] < 0 else "up"
        L = np.concatenate([[0], np.cumsum(np.hypot(*np.diff(q, axis=0).T))])
        share = L / max(L[-1], 1e-9)
        tv = np.interp(share, frac[:, i], days, left=0.0, right=days[-1])
        r, sg = float(_at(rate_g, mid[i:i + 1])[0]), float(_at(sigma_g, mid[i:i + 1])[0])
        sn = abs(r) / sg if np.isfinite(sg) and sg > 0 else 0.0
        segv = _at(np.nan_to_num(seg_g, nan=0), mid[i:i + 1])[0]
        feats.append({"type": "Feature", "properties": {
            "rate_cm_yr": round(r, 2), "sigma_cm_yr": round(sg, 2) if np.isfinite(sg) else None,
            "snr": round(sn, 2), "segments": int(round(segv)),
            "trust": "strong" if sn >= ts["snr_strong"] else "weak" if sn >= ts["snr_weak"] else "below",
            "dark": round(s.dark, 3), "sense": sense, "t": np.round(tv, 1).tolist()},
            "geometry": {"type": "LineString", "coordinates": _lonlat(q, T, crs).tolist()}})
    return _fc(feats), dict(interval=float(interval), deepest_cm=float(np.nanpercentile(final_s[valid], 1) * 100)
                            if valid.any() else 0.0, frames_s=frames_s, valid=valid)


def contours_layer(cfg, v: Velocity, frames_s, valid, T, crs, simplify=2):
    """Labelled iso-lines of displacement (like depth lines on a sea chart), one set per date.

    cfg['contours'] = {"step_cm": 5, "major_cm": 10}. Features: {level_cm, major, d (date index)}.
    Lines are drawn only where the strokes are allowed (coherent, supported, above the noise floor),
    so a contour never claims more than the hachures do.
    """
    import contourpy
    p = cfg.get("contours", {})
    step, major = p.get("step_cm", 5), p.get("major_cm", 10)
    sink_only = cfg["velocity"].get("sink_only", True)
    feats = []
    for k, fr in enumerate(frames_s):
        if k == 0:
            continue
        z = np.where(valid, fr * 100, np.nan)
        if not np.isfinite(z).any():
            continue
        lo, hi = np.nanmin(z), np.nanmax(z)
        levels = [l for l in np.arange(np.floor(lo / step) * step, hi + step, step) if l != 0 and (l < 0 or not sink_only)]
        gen = contourpy.contour_generator(z=z, line_type=contourpy.LineType.Separate)
        for lvl in levels:
            for ln in gen.lines(float(lvl)):
                if len(ln) < 6:
                    continue
                q = ln[::simplify] if len(ln) > 2 * simplify else ln
                feats.append({"type": "Feature", "properties": {"level_cm": int(lvl), "major": bool(int(lvl) % major == 0),
                                                                "d": k},
                              "geometry": {"type": "LineString", "coordinates": _lonlat(q, T, crs).tolist()}})
    return _fc(feats)


def _gnss_series(row, dates, half_days=3.0):
    """GNSS vertical at each radar date (median of +-3 days of daily solutions), cm, relative to the first date."""
    path = ROOT / row["series"]
    t, up = [], []
    for line in open(path):
        f = line.split()
        if f[0] == "site" or len(f) < 13:
            continue
        yr = float(f[2])
        y0 = int(yr)
        d0 = datetime(y0, 1, 1)
        t.append(d0 + (datetime(y0 + 1, 1, 1) - d0) * (yr - y0))
        up.append(float(f[11]) + float(f[12]))
    t, up = np.array(t), np.array(up)
    vals = []
    for d in dates:
        m = np.array([abs((ti - d).total_seconds()) <= half_days * 86400 for ti in t])
        vals.append(float(np.median(up[m])) if m.any() else np.nan)
    vals = (np.array(vals) - vals[0]) * 100
    return [round(x, 1) if np.isfinite(x) else None for x in vals]


def surface_layer(cfg, v: Velocity, frames_s, mask, T, crs):
    """The smoothed, deramped displacement surface on a light lon/lat grid, one frame per date (mm, int).

    It is exactly the field the hachures and depth lines are drawn from, so a cross-section in the browser
    agrees with the map. Also carries a 1-sigma rate grid (for an error band) and GNSS series as ground truth.
    cfg['surface'] = {"step_m": 500, "section": [[lon, lat], [lon, lat]]  (default cut for the story)}.
    """
    p = cfg.get("surface", {})
    step = p.get("step_m", 500)
    w, s_, e, n = cfg["bbox"]
    lat_c = (s_ + n) / 2
    dlat = step / 111320.0
    dlon = step / (111320.0 * np.cos(np.radians(lat_c)))
    lons = np.arange(w, e + dlon / 2, dlon)
    lats = np.arange(s_, n + dlat / 2, dlat)
    LO, LA = np.meshgrid(lons, lats)
    X, Y = warp_transform("EPSG:4326", crs, LO.ravel().tolist(), LA.ravel().tolist())
    col = (np.array(X) - T.c) / T.a - 0.5
    row = (np.array(Y) - T.f) / T.e - 0.5
    pts = np.column_stack([col, row])
    inside = (col >= 0) & (row >= 0) & (col <= mask.shape[1] - 1) & (row <= mask.shape[0] - 1)
    ok = inside.copy()
    ok[inside] = mask[np.round(row[inside]).astype(int), np.round(col[inside]).astype(int)]
    frames = []
    for fr in frames_s:
        val = _at(np.nan_to_num(fr), pts) * 1000
        frames.append([int(round(x)) if k else None for x, k in zip(val, ok)])
    sg = gridding.to_grid(np.asarray(v.sigma * 1000, float), gridding.transform_from_centres(v.x, v.y),
                          f"EPSG:{v.epsg}", T, crs, mask.shape)
    sg = _smooth_nan(np.nan_to_num(sg), 250 / cfg["dem_res"], np.isfinite(sg))[0]
    sval = _at(sg, pts)
    gn = [{"station": g["station"], "lon": g["lon"], "lat": g["lat"], "cm": _gnss_series(g, v.dates)}
          for g in cfg.get("gnss", []) if g.get("series")]
    return {"kind": "surface", "lon0": round(lons[0], 6), "lat0": round(lats[0], 6), "dlon": dlon, "dlat": dlat,
            "nx": len(lons), "ny": len(lats), "unit": "mm", "dates": [d.date().isoformat() for d in v.dates],
            "frames": frames, "sigma_mm_yr": [int(round(x)) if k else None for x, k in zip(sval, ok)],
            "years": (v.dates[-1] - v.dates[0]).days / 365.25, "section": p.get("section"), "gnss": gn}


def headline(cfg, v: Velocity, block_m=500.0):
    """The number we dare to say out loud: the 500 m neighbourhood where we are most sure the ground sinks fast.

    Blocks = median of the raw stack rate (>= 80% valid pixels, coherent), then the median of 3x3 blocks (~1.5 km,
    >= 6 good blocks). Headline = the most negative 95% upper bound (rate + 2 sigma): 'here it sinks at least this fast'.
    (First version used single blocks with >= 50% valid and picked a sparse low-coherence pocket at -53 cm/yr
    surrounded by -30..-35; see notes/experiments.md D-depth.) Also the 1st percentile
    over all pixels ('1% of the area sinks faster than ...'). Diagnostic rule fixed 11.11 after scripts/diag_depth.py.
    """
    import warnings
    ok = np.isfinite(v.rate) & (v.coherence >= cfg["velocity"].get("coh_min_draw", 0.45))
    r = np.where(ok, v.rate * 100, np.nan)
    px = abs(v.x[1] - v.x[0])
    B = max(1, int(round(block_m / px)))
    H, W = (r.shape[0] // B) * B, (r.shape[1] // B) * B
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        blk = r[:H, :W].reshape(H // B, B, W // B, B)
        med = np.nanmedian(blk, axis=(1, 3))
        sg = np.nanmedian((v.sigma * 100)[:H, :W].reshape(H // B, B, W // B, B), axis=(1, 3))
        cum = np.nanmedian(((v.cum - v.cum[:1]) * 100)[:, :H, :W].reshape(len(v.dates), H // B, B, W // B, B), axis=(2, 4))
    med[np.isfinite(blk).sum((1, 3)) < 0.8 * B * B] = np.nan          # well-covered blocks only
    # neighbourhood = 3x3 blocks (~1.5 km), needs >= 6 of 9 good blocks: one sparse pocket cannot set the headline
    from numpy.lib.stride_tricks import sliding_window_view
    pad = np.pad(med, 1, constant_values=np.nan)
    win = sliding_window_view(pad, (3, 3)).reshape(*med.shape, 9)
    neigh = np.nanmedian(win, -1)
    neigh[np.isfinite(win).sum(-1) < 6] = np.nan
    med = neigh
    bound = med + 2 * sg
    if not np.isfinite(bound).any():
        return {}
    i, j = np.unravel_index(np.nanargmin(bound), bound.shape)
    xs, ys = v.x[:W].reshape(-1, B).mean(1), v.y[:H].reshape(-1, B).mean(1)
    lon, lat = warp_transform(f"EPSG:{v.epsg}", "EPSG:4326", [xs[j]], [ys[i]])
    return {"rate_cm_yr": round(float(med[i, j]), 1), "sigma_cm_yr": round(float(sg[i, j]), 1),
            "at_least_cm_yr": round(float(bound[i, j]), 1), "lon": round(lon[0], 4), "lat": round(lat[0], 4),
            "cum_cm": round(float(cum[-1, i, j]), 1), "block_m": block_m,
            "p1_cm_yr": round(float(np.nanpercentile(r, 1)), 1),
            "rule": "~1.5 km neighbourhood (3x3 well-covered 500 m blocks) with the most negative rate + 2 sigma; single pass, LOS taken as vertical"}


def grid_layer(cfg, v: Velocity):
    """Coarse cells (grid_km) for click-a-point and for raster-style views: rate, sigma, segments, series."""
    px = abs(v.x[1] - v.x[0])
    B = max(1, int(round(cfg["grid_km"] * 1000 / px)))
    H, W = (v.rate.shape[0] // B) * B, (v.rate.shape[1] // B) * B

    def blk(a, fn=np.nanmedian):
        with np.errstate(all="ignore"):
            return fn(a[..., :H, :W].reshape(*a.shape[:-2], H // B, B, W // B, B), axis=(-3, -1))
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        rate = blk(v.rate) * 100
        sigma = blk(v.sigma) * 100
        seg = blk(np.where(v.segments > 0, v.segments, np.nan))
        cum = blk(v.cum - v.cum[:1]) * 100
    xs = v.x[:W].reshape(-1, B).mean(1)
    ys = v.y[:H].reshape(-1, B).mean(1)
    XX, YY = np.meshgrid(xs, ys)
    ok = np.isfinite(rate)
    lon, lat = warp_transform(f"EPSG:{v.epsg}", "EPSG:4326", XX[ok].tolist(), YY[ok].tolist())
    cells = []
    for k, (i, j) in enumerate(zip(*np.nonzero(ok))):
        cells.append([round(lon[k], 5), round(lat[k], 5), round(float(rate[i, j]), 2),
                      round(float(sigma[i, j]), 2) if np.isfinite(sigma[i, j]) else None,
                      int(seg[i, j]) if np.isfinite(seg[i, j]) else None,
                      [round(float(c), 1) if np.isfinite(c) else None for c in cum[:, i, j]]])
    return {"cell_km": cfg["grid_km"], "columns": ["lon", "lat", "rate_cm_yr", "sigma_cm_yr", "segments", "cum_cm"],
            "dates": [d.date().isoformat() for d in v.dates], "cells": cells}


def rhythm_layer(cfg, g):
    """Rows of the "partitura": one time series per named place, taken from the grid cells (cum_cm) around it.

    cfg['rhythm'] = {r_km, min_cells, rows: [{name, lon, lat, ref?}]}. Per row: median curve of the cells within r_km,
    12-day steps, least-squares rate, sigma = median cell sigma (cells this close share their atmosphere, so no 1/sqrt(n)),
    rates of the first and second half, and `still` when |rate| < 2 sigma (drawn and played as silence).
    Rows come out sorted from the fastest sinking down. Style-free: the front-end decides how a row looks and sounds."""
    p = cfg["rhythm"]
    r_km, min_cells = p.get("r_km", 1.0), p.get("min_cells", 4)
    cols = g["columns"]
    iR, iS, iC = cols.index("rate_cm_yr"), cols.index("sigma_cm_yr"), cols.index("cum_cm")
    C = g["cells"]
    lon = np.array([c[0] for c in C]); lat = np.array([c[1] for c in C])
    days = np.array([(datetime.fromisoformat(d) - datetime.fromisoformat(g["dates"][0])).days for d in g["dates"]], float)
    yrs = days / 365.25

    def fit(t, y):
        ok = np.isfinite(y)
        return float(np.polyfit(t[ok], y[ok], 1)[0]) if ok.sum() >= 3 else None

    rows = []
    for spec in p["rows"]:
        d = np.hypot((lon - spec["lon"]) * np.cos(np.radians(spec["lat"])) * 111.32, (lat - spec["lat"]) * 110.57)
        idx = np.nonzero(d <= r_km)[0]
        if len(idx) < min_cells:
            continue
        S = np.array([[np.nan if x is None else x for x in C[i][iC]] for i in idx], float)
        with np.errstate(all="ignore"):
            ser = np.nanmedian(S, 0)
        rate = fit(yrs, ser)
        sig = [C[i][iS] for i in idx if C[i][iS] is not None]
        sigma = float(np.median(sig)) if sig else None
        h = len(yrs) // 2
        rows.append({"name": spec["name"], "lon": spec["lon"], "lat": spec["lat"], "r_km": r_km, "cells": int(len(idx)),
                     "ref": bool(spec.get("ref", False)),
                     "series_cm": [round(float(x), 2) if np.isfinite(x) else None for x in ser],
                     "rate_cm_yr": round(rate, 2) if rate is not None else None,
                     "sigma_cm_yr": round(sigma, 2) if sigma is not None else None,
                     "rate_first_half": round(fit(yrs[:h + 1], ser[:h + 1]) or 0, 2),
                     "rate_second_half": round(fit(yrs[h:], ser[h:]) or 0, 2),
                     "still": bool(rate is None or sigma is None or abs(rate) < 2 * sigma)})
    rows.sort(key=lambda r: r["rate_cm_yr"] if r["rate_cm_yr"] is not None else 0)
    return {"dates": g["dates"], "unit": "cm", "r_km": r_km, "rows": rows,
            "rule": "row = median of the 500 m cells within r_km; still = |rate| < 2 sigma (set before looking); "
                    "a ref row is the zero point of the map, so it is flat by construction"}


def gnss_check(cfg, v: Velocity, radius_m=250.0):
    rows = []
    if not cfg["gnss"]:
        return rows
    xs, ys = warp_transform("EPSG:4326", f"EPSG:{v.epsg}", [g["lon"] for g in cfg["gnss"]], [g["lat"] for g in cfg["gnss"]])
    XX, YY = np.meshgrid(v.x, v.y)
    for g, sx, sy in zip(cfg["gnss"], xs, ys):
        m = np.hypot(XX - sx, YY - sy) <= radius_m
        n = float(np.nanmedian(v.rate[m]) * 100) if np.isfinite(v.rate[m]).sum() >= 3 else None
        rows.append({**g, "nisar_cm_yr": round(n, 2) if n is not None else None,
                     "diff_cm_yr": round(n - g["rate_cm_yr"], 2) if n is not None else None})
    return rows


# ------------------------------------------------------------------ site export

def _contour_mask(cfg, v, st):
    """Contours may cross the small holes between strokes, but only inside the region the strokes cover."""
    from scipy.ndimage import binary_closing, binary_fill_holes
    return binary_fill_holes(binary_closing(st["valid"], iterations=3))


def _dem(cfg):
    if cfg.get("dem_path"):
        return dem.load(ROOT / cfg["dem_path"])
    demf = ROOT / "data/raw" / f"copdem_{cfg.get('dem_name', cfg['id'])}_{int(cfg['dem_res'])}m.tif"
    dem.fetch_copdem(tuple(cfg["bbox"]), demf, res=cfg["dem_res"])
    return dem.load(demf)


# ------------------------------------------------------------------ context layers (lake, cracks, city): style-free lines + labels

def _match(props, where):
    return all(str(props.get(k, "")).strip() in (v if isinstance(v, list) else [v]) for k, v in (where or {}).items())


def _in_bbox(coords, bbox, pad=0.02):
    w, s, e, n = bbox
    return any(w - pad <= x <= e + pad and s - pad <= y <= n + pad for x, y in coords)


def _rings_as_lines(geom):
    t, c = geom["type"], geom["coordinates"]
    if t == "LineString":
        return [c]
    if t == "MultiLineString" or t == "Polygon":
        return list(c)
    if t == "MultiPolygon":
        return [r for poly in c for r in poly]
    return []


def _read_features(spec):
    """GeoJSON (any CRS) or raw Overpass JSON -> list of (geometry in lon/lat, properties)."""
    path = ROOT / spec["src"]
    raw = json.loads(path.read_bytes().decode("utf-8", errors="replace"))
    out = []
    if "elements" in raw:                                  # Overpass with `out geom`
        for el in raw["elements"]:
            tags = el.get("tags", {})
            if el["type"] == "node":
                out.append(({"type": "Point", "coordinates": [el["lon"], el["lat"]]}, tags))
            elif el["type"] == "way" and el.get("geometry"):
                out.append(({"type": "LineString", "coordinates": [[p["lon"], p["lat"]] for p in el["geometry"]]}, tags))
            elif el["type"] == "relation":
                ls = [[[p["lon"], p["lat"]] for p in m["geometry"]] for m in el.get("members", [])
                      if m.get("type") == "way" and m.get("geometry") and m.get("role", "") in ("", "outer", "forward", "backward")]
                if ls:
                    out.append(({"type": "MultiLineString", "coordinates": ls}, tags))
        return out
    src_crs = spec.get("src_crs") or raw.get("crs", {}).get("properties", {}).get("name", "EPSG:4326")
    src_crs = src_crs.replace("urn:ogc:def:crs:", "").replace("::", ":")
    from rasterio.warp import transform_geom
    for f in raw["features"]:
        if f.get("geometry") is None:
            continue
        g = f["geometry"] if src_crs in ("EPSG:4326", "OGC:1.3:CRS84") else transform_geom(src_crs, "EPSG:4326", f["geometry"])
        out.append((g, f.get("properties") or {}))
    return out


def _thin(line, step=1e-4):
    """Drop vertices closer than `step` degrees (~10 m) to the previous kept one; keep the ends."""
    a = np.asarray(line, float)
    if len(a) <= 2:
        return np.round(a, 5).tolist()
    keep, last = [0], a[0]
    for i in range(1, len(a) - 1):
        if abs(a[i, 0] - last[0]) + abs(a[i, 1] - last[1]) >= step:
            keep.append(i)
            last = a[i]
    keep.append(len(a) - 1)
    return np.round(a[keep], 5).tolist()


def _pick_name(pr, keys, latin_only=False):
    """First non-empty name among `keys` (a key or a list, e.g. ["name:en", "name"]); latin_only skips other scripts."""
    for k in keys if isinstance(keys, list) else [keys]:
        v = pr.get(k)
        if v and (not latin_only or all(ord(ch) < 0x0370 or 0x1E00 <= ord(ch) <= 0x1EFF or ch in "’Ⅱ" for ch in v)):
            return v
    return None


def _num(v):
    try:
        return float(str(v).replace(",", "."))
    except (TypeError, ValueError):
        return None


def lines_layer(cfg, spec):
    """Context lines + label points, filtered by `where`, clipped to the bbox (whole features kept).

    Points: `name_key` (key or list, first found wins), `latin_only`, `min_tag` ({tag: lowest number}, e.g. {"ele": 6000}),
    `point_format` (e.g. "{name} {ele}"). Lines: `label_from` may also be a list.
    """
    feats, labels, seen = [], {}, set()
    for g, pr in _read_features(spec):
        if not _match(pr, spec.get("where")):
            continue
        if spec.get("dedupe_key"):                       # e.g. a metro line mapped once per direction
            key = pr.get(spec["dedupe_key"])
            if key in seen:
                continue
            seen.add(key)
        if g["type"] == "Point":
            if any((_num(pr.get(k)) or -1e9) < v for k, v in (spec.get("min_tag") or {}).items()):
                continue
            if spec.get("points_as") == "dots":            # e.g. stations: a mark, no name (names would bury the map)
                if _in_bbox([g["coordinates"]], cfg["bbox"], 0):
                    feats.append({"type": "Feature", "properties": {"dot": 1},
                                  "geometry": {"type": "Point", "coordinates": np.round(g["coordinates"], 5).tolist()}})
                continue
            nm = _pick_name(pr, spec.get("name_key", "name"), spec.get("latin_only"))
            if spec.get("points_as_labels", True) and nm and _in_bbox([g["coordinates"]], cfg["bbox"], 0):
                if spec.get("point_format"):
                    ele = _num(pr.get("ele"))
                    nm = spec["point_format"].format(name=nm, ele=f"{int(round(ele)):,}".replace(",", " ") if ele else "")
                feats.append({"type": "Feature", "properties": {"name": nm.strip(),
                                                               "rank": spec.get("rank", {}).get(str(pr.get(spec.get("rank_key", "place"))), 1)},
                              "geometry": {"type": "Point", "coordinates": np.round(g["coordinates"], 5).tolist()}})
            continue
        lines = [ln for ln in _rings_as_lines(g) if len(ln) >= 2 and _in_bbox(ln, cfg["bbox"])]
        if not lines:
            continue
        if spec.get("min_len_m"):
            L = sum(np.hypot(np.diff(np.asarray(ln)[:, 0]) * 104e3, np.diff(np.asarray(ln)[:, 1]) * 111e3).sum() for ln in lines)
            if L < spec["min_len_m"]:
                continue
        feats.append({"type": "Feature", "properties": {},
                      "geometry": {"type": "MultiLineString", "coordinates": [_thin(ln, spec.get("thin_deg", 1e-4)) for ln in lines]}})
        lk = spec.get("label_from")                        # one label per name, at the middle of its longest piece
        lv = _pick_name(pr, lk, spec.get("latin_only")) if lk else None
        if lv and (not spec.get("label_only") or lv in spec["label_only"]):
            name = spec.get("label_format", "{}").format(lv)
            longest = max(lines, key=len)
            if name not in labels or len(longest) > labels[name][1]:
                labels[name] = (longest[len(longest) // 2], len(longest))
    for name, (p, _) in labels.items():
        feats.append({"type": "Feature", "properties": {"name": name, "rank": 2},
                      "geometry": {"type": "Point", "coordinates": np.round(p, 5).tolist()}})
    return _fc(feats)


def _join_rings(lines, tol=1e-7):
    """Join line pieces that share end points into closed rings (OSM multipolygon members); open leftovers are dropped."""
    todo = [list(map(tuple, ln)) for ln in lines if len(ln) >= 2]
    rings = []
    while todo:
        cur = todo.pop(0)
        grown = True
        while grown and (abs(cur[0][0] - cur[-1][0]) > tol or abs(cur[0][1] - cur[-1][1]) > tol):
            grown = False
            for i, ln in enumerate(todo):
                if ln[0] == cur[-1]:
                    cur += ln[1:]
                elif ln[-1] == cur[-1]:
                    cur += ln[::-1][1:]
                elif ln[-1] == cur[0]:
                    cur = ln[:-1] + cur
                elif ln[0] == cur[0]:
                    cur = ln[::-1][:-1] + cur
                else:
                    continue
                todo.pop(i)
                grown = True
                break
        if len(cur) >= 4 and abs(cur[0][0] - cur[-1][0]) <= tol and abs(cur[0][1] - cur[-1][1]) <= tol:
            rings.append([list(q) for q in cur])
    return rings


def waterlines_layer(cfg, spec):
    """Engraver's water lines inside polygons: the shore plus offset lines whose spacing grows away from it.

    `levels_m`: distances from the shore; property k = 0 at the shore .. 1 at the last line (front-end fades with k).
    """
    from contourpy import contour_generator
    from rasterio.features import rasterize
    from scipy.ndimage import distance_transform_edt
    z, T, crs = _dem(cfg)
    from rasterio.warp import transform_geom
    polys = []
    for g, pr in _read_features(spec):
        if not _match(pr, spec.get("where")):
            continue
        if g["type"] in ("Polygon", "MultiPolygon"):
            polys.append(g)
        elif g["type"] in ("LineString", "MultiLineString"):         # Overpass: closed ways, relation member ways
            rings = _join_rings(_rings_as_lines(g))
            if rings:
                polys.append({"type": "MultiPolygon", "coordinates": [[r] for r in rings]})
    shapes = [(transform_geom("EPSG:4326", crs, g), 1) for g in polys]
    mask = rasterize(shapes, out_shape=z.shape, transform=T, fill=0, dtype="uint8").astype(bool)
    dist = distance_transform_edt(mask) * cfg["dem_res"]
    levels = spec.get("levels_m", [0, 80, 180, 300, 450, 650, 900])
    gen = contour_generator(z=np.where(mask, dist, -1.0))
    feats = []
    for i, lv in enumerate(levels):
        for ln in gen.lines(max(lv, cfg["dem_res"] * 0.5)):
            if len(ln) < spec.get("min_pts", 6):
                continue
            q = _lonlat(np.asarray(ln)[::2] if len(ln) > 12 else np.asarray(ln), T, crs)
            feats.append({"type": "Feature", "properties": {"k": round(i / max(1, len(levels) - 1), 3)},
                          "geometry": {"type": "LineString", "coordinates": q.tolist()}})
    return _fc(feats)


def relief_info(cfg, out, verbose=True):
    """Height above sea level for a map frame / cursor readout, plus spot heights (summits) as labelled points.

    cfg['relief_info'] = {step_m, peaks, window_m, prominence_m, spacing_m}. Writes elev.json and ctx_peaks.geojson.
    """
    from scipy.ndimage import maximum_filter, minimum_filter
    p = cfg.get("relief_info")
    if not p:
        return []
    z, T, crs = _dem(cfg)
    res = cfg["dem_res"]
    w, s, e, n = cfg["bbox"]
    step = p.get("step_m", 300)
    dlat = step / 111_320
    dlon = step / (111_320 * np.cos(np.radians((s + n) / 2)))
    lons, lats = np.arange(w, e + 1e-9, dlon), np.arange(s, n + 1e-9, dlat)
    LO, LA = np.meshgrid(lons, lats)
    xs, ys = warp_transform("EPSG:4326", crs, LO.ravel().tolist(), LA.ravel().tolist())
    cols = (np.array(xs) - T.c) / T.a - 0.5
    rows = (np.array(ys) - T.f) / T.e - 0.5
    zz = map_coordinates(np.where(np.isfinite(z), z, -9999), [rows, cols], order=1, mode="constant", cval=-9999)
    grid = [None if v < -1000 else int(round(v)) for v in zz]
    _write(out, "elev.json", {"kind": "elevation", "unit": "m", "lon0": w, "lat0": s, "dlon": dlon, "dlat": dlat,
                              "nx": len(lons), "ny": len(lats), "z": grid})
    # summits: local maxima over a window, prominent over their surroundings, spaced apart
    win = max(3, int(p.get("window_m", 3000) / res) | 1)
    zf = np.where(np.isfinite(z), z, -9999)
    peak = (zf == maximum_filter(zf, size=win)) & (zf - minimum_filter(np.where(np.isfinite(z), z, 9999), size=win * 2 + 1)
                                                   >= p.get("prominence_m", 300))
    edge = win // 2 + 2                                    # a cut-off slope at the DEM edge is not a summit
    peak[:edge], peak[-edge:], peak[:, :edge], peak[:, -edge:] = False, False, False, False
    r, c = np.nonzero(peak)
    order = np.argsort(-zf[r, c])
    chosen = []
    for i in order:
        if len(chosen) >= p.get("peaks", 12):
            break
        if all(np.hypot(r[i] - rr, c[i] - cc) * res >= p.get("spacing_m", 4000) for rr, cc in chosen):
            chosen.append((r[i], c[i]))
    pts = np.array([[cc, rr] for rr, cc in chosen], float) if chosen else np.zeros((0, 2))
    ll = _lonlat(pts, T, crs) if len(pts) else []
    feats = [{"type": "Feature", "properties": {"name": f"▲ {int(round(zf[rr, cc])):,}".replace(",", " "), "rank": 1},
              "geometry": {"type": "Point", "coordinates": list(map(float, q))}} for (rr, cc), q in zip(chosen, ll)]
    if verbose:
        print(f"  relief: elev {len(lons)} x {len(lats)} at {step} m, {len(feats)} summits")
    elev = [{"id": "elevation", "kind": "elevation", "file": "elev.json", "label": {"en": "Height", "uk": "Висота"}}]
    if not feats:                                          # peaks: 0 -> named summits come from elsewhere (e.g. OSM)
        return elev
    _write(out, "ctx_peaks.geojson", _fc(feats))
    return elev + [
            {"id": "peaks", "kind": "lines", "file": "ctx_peaks.geojson", "style": "peak",
             "label": {"en": "Summits, m above sea level", "uk": "Вершини, м над рівнем моря"}, "visible": True, "z": 3.6}]


def landcover_layers(cfg, out, verbose=True):
    """Map symbols from a land-cover raster (ESA WorldCover classes), style-free:

    trees: jittered points where the tree share of a cell is high (front-end draws a little tree);
    built: 45-degree hatch lines clipped to built-up land (the old-map "tone" of a town).
    cfg['landcover'] = {src, attribution, trees: {class, cell_m, min_share}, built: {class, block_m, min_share, spacing_m}}
    """
    import rasterio
    from scipy.ndimage import binary_closing, binary_opening
    p = cfg.get("landcover")
    if not p:
        return []
    with rasterio.open(ROOT / p["src"]) as ds:
        a = ds.read(1)
        T = ds.transform
    w, s, e, n = cfg["bbox"]
    lat_m = 111_320.0
    lon_m = lat_m * np.cos(np.radians((s + n) / 2))
    px_m = abs(T.a) * lat_m                                       # ~9 m for WorldCover at 1/12000 degree
    entries = []
    rng = np.random.default_rng(0)

    def blocks(cls, block_m):
        k = max(1, int(round(block_m / px_m)))
        H, W = (a.shape[0] // k) * k, (a.shape[1] // k) * k
        return (a[:H, :W] == cls).reshape(H // k, k, W // k, k).mean((1, 3)), k

    q = p.get("trees")
    if q:
        frac, k = blocks(q.get("class", 10), q.get("cell_m", 300))
        r, c = np.nonzero(frac >= q.get("min_share", 0.6))
        jit = rng.uniform(-0.3, 0.3, (len(r), 2))
        lon = T.c + (c + 0.5 + jit[:, 0]) * k * T.a
        lat = T.f + (r + 0.5 + jit[:, 1]) * k * T.e
        feats = [{"type": "Feature", "properties": {"tree": round(float(f), 2)},
                  "geometry": {"type": "Point", "coordinates": [round(float(x), 5), round(float(y), 5)]}}
                 for x, y, f in zip(lon, lat, frac[r, c])]
        _write(out, "ctx_trees.geojson", _fc(feats))
        entries.append({"id": "trees", "kind": "lines", "file": "ctx_trees.geojson", "style": "tree",
                        "label": q.get("label", {"en": "Woods and parks", "uk": "Ліси й парки"}), "visible": q.get("visible", False), "by_theme": True,
                        "z": q.get("z", 1.25), "attribution": p.get("attribution")})
        if verbose:
            print(f"  landcover trees: {len(feats)} symbols")
    q = p.get("built")
    if q:
        frac, k = blocks(q.get("class", 50), q.get("block_m", 60))
        m = frac >= q.get("min_share", 0.5)
        m = binary_closing(binary_opening(m, iterations=1), iterations=2)
        bm = k * px_m                                             # block size, m
        sp = q.get("spacing_m", 220) / bm                         # hatch spacing in blocks
        H, W = m.shape
        feats = []
        # lines x - y = const (45 degrees on the ground, rows grow southwards); walk each diagonal, keep runs inside the mask
        for d in np.arange(-H, W, sp):
            i0 = max(0, int(np.ceil(-d)))
            i = np.arange(i0, H)
            j = np.round(i + d).astype(int)
            ok = (j >= 0) & (j < W)
            i, j = i[ok], j[ok]
            if len(i) < 2:
                continue
            inside = m[i, j]
            edges = np.flatnonzero(np.diff(np.r_[0, inside.astype(int), 0]))
            for a0, a1 in zip(edges[::2], edges[1::2]):
                if a1 - a0 < max(2, int(q.get("min_run_m", 120) / bm)):
                    continue
                x0, y0 = T.c + (j[a0] + 0.5) * k * T.a, T.f + (i[a0] + 0.5) * k * T.e
                x1, y1 = T.c + (j[a1 - 1] + 0.5) * k * T.a, T.f + (i[a1 - 1] + 0.5) * k * T.e
                feats.append([[round(x0, 5), round(y0, 5)], [round(x1, 5), round(y1, 5)]])
        fc = _fc([{"type": "Feature", "properties": {}, "geometry": {"type": "MultiLineString", "coordinates": feats}}])
        _write(out, "ctx_built.geojson", fc)
        entries.append({"id": "built", "kind": "lines", "file": "ctx_built.geojson", "style": "built",
                        "label": q.get("label", {"en": "Built-up land", "uk": "Забудова"}), "visible": q.get("visible", False), "by_theme": True,
                        "z": q.get("z", 1.1), "attribution": p.get("attribution")})
        if verbose:
            print(f"  landcover built: {len(feats)} hatch segments, {(out / 'ctx_built.geojson').stat().st_size / 1e6:.2f} MB")
    q = p.get("blocks")
    if q:
        entries.append(_town_blocks(a, T, q, px_m, cfg, out, p.get("attribution"), verbose))
    return entries


def _town_blocks(a, T, q, px_m, cfg, out, attribution, verbose=True):
    """Old town-plan sign as an image: built-up land (WorldCover) with the streets (OSM) cut out as paper-coloured gaps.

    Every street is used, so the blocks come out at their real size; the raster cell (cell_m, ~18 m) gives the
    stair-stepped "little squares" of a Siegfried sheet up close and a Dufour-like town tone zoomed out.
    Writes ctx_blocks.png: white + alpha (alpha = how built-up, 0 on streets); the front-end tints it with the theme ink.
    q = {class, cell_m, min_share, roads, streets, gap_m {highway class: width m, "street": width m}}
    """
    import rasterio.features
    from PIL import Image
    from rasterio.transform import Affine
    from scipy.ndimage import binary_closing, binary_dilation, uniform_filter
    k = max(1, int(round(q.get("cell_m", 18) / px_m)))
    H, W = (a.shape[0] // k) * k, (a.shape[1] // k) * k
    frac = (a[:H, :W] == q.get("class", 50)).reshape(H // k, k, W // k, k).mean((1, 3))
    T2 = T * Affine.scale(k)
    cell = k * px_m
    built = binary_closing(frac >= q.get("min_share", 0.5), iterations=1)
    gaps = dict(q.get("gap_m", {}))
    street_w = gaps.pop("street", 0)

    def burn(geoms, width_m):
        m = rasterio.features.rasterize(geoms, out_shape=built.shape, transform=T2, all_touched=width_m > cell).astype(bool)
        it = max(0, int(round((width_m / cell - 1) / 2)))
        return binary_dilation(m, iterations=it) if it else m

    streets = np.zeros(built.shape, bool)
    roads = _read_features({"src": q["roads"]})
    for cls, gap in gaps.items():
        geoms = [g for g, t in roads if t.get("highway") == cls and g["type"] == "LineString"]
        if geoms:
            streets |= burn(geoms, gap)
    if q.get("streets") and street_w:
        streets |= burn([g for g, _ in _read_features({"src": q["streets"]}) if g["type"] == "LineString"], street_w)
    blocks = built & ~streets
    tone = uniform_filter(frac, 5)                                   # neighbourhood density: old centre darker than suburbs
    alpha = np.where(blocks, np.clip(110 + 145 * tone, 0, 255), 0).astype("uint8")
    img = Image.fromarray(np.dstack([np.full_like(alpha, 255), alpha]), "LA")
    img.save(out / "ctx_blocks.png", optimize=True)
    w, n = T2.c, T2.f
    e, s_ = w + T2.a * alpha.shape[1], n + T2.e * alpha.shape[0]
    if verbose:
        print(f"  landcover blocks: {alpha.shape[1]} x {alpha.shape[0]} px at {cell:.0f} m, "
              f"{blocks.mean():.1%} block, {streets.mean():.1%} street, {(out / 'ctx_blocks.png').stat().st_size / 1e6:.2f} MB")
    return {"id": "blocks", "kind": "image", "file": "ctx_blocks.png", "bounds": [w, s_, e, n], "style": "blocks",
            "label": q.get("label", {"en": "Blocks and streets", "uk": "Квартали й вулиці"}),
            "visible": q.get("visible", False), "by_theme": True, "z": q.get("z", 1.1),
            "attribution": (attribution or "") + " · © OpenStreetMap contributors",
            "note": q.get("note", {"en": "Built-up land (ESA WorldCover, 2021) with the streets (OpenStreetMap) left white, "
                                         "drawn like an old town plan; darker = denser.",
                                   "uk": "Забудова (ESA WorldCover, 2021), а вулиці (OpenStreetMap) лишені білими, "
                                         "як на старому плані міста; темніше = щільніше."})}


def context_layers(cfg, out, verbose=True):
    """cfg['context'] = [{id, kind: 'lines'|'waterlines', src, ...}] -> files + manifest entries (kind 'lines')."""
    entries = []
    for spec in cfg.get("context", []):
        fc = waterlines_layer(cfg, spec) if spec.get("kind") == "waterlines" else lines_layer(cfg, spec)
        name = f"ctx_{spec['id']}.geojson"
        _write(out, name, fc)
        e = {"id": spec["id"], "kind": "lines", "file": name, "style": spec.get("style", spec["id"]),
             "label": spec.get("label", spec["id"]), "visible": spec.get("visible", True), "z": spec.get("z", 1.2)}
        for k in ("group", "group_label", "note", "attribution"):
            if spec.get(k):
                e[k] = spec[k]
        entries.append(e)
        if verbose:
            print(f"  context {spec['id']}: {len(fc['features'])} features, {(out / name).stat().st_size / 1e6:.2f} MB")
    return entries


def _write(out, name, obj):
    (out / name).write_text(json.dumps(obj, separators=(",", ":"), ensure_ascii=False), encoding="utf-8")


def export_site(cfg, v: Velocity | None, out_dir, verbose: bool = True) -> Path:
    """Write every layer listed in cfg['layers'] + manifest.json. `v` is needed only for velocity/grid/GNSS."""
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    L = cfg["layers"]
    if {"terrain", "velocity"} & set(L):
        z, T, crs = _dem(cfg)
    layers, stats, dates, n_pairs = [], {}, [], 0
    if v is not None:
        dates, n_pairs = [d.date().isoformat() for d in v.dates], v.n_pairs
    if "terrain" in L:
        tf = out / "terrain.geojson"
        if cfg.get("_keep_terrain") and tf.exists():          # relief depends only on the DEM + terrain{} settings
            fc = json.loads(tf.read_text(encoding="utf-8"))
        else:
            fc = terrain_layer(cfg, z, T, crs)
            _write(out, "terrain.geojson", fc)
        layers.append({"id": "terrain", "kind": "hachures-static", "file": "terrain.geojson",
                       "label": {"en": "Relief", "uk": "Рельєф"}, "props": ["dark"]})
        if verbose:
            print(f"  terrain: {len(fc['features'])} strokes")
    if "velocity" in L:
        fc, st = velocity_layer(cfg, v, z, T, crs)
        _write(out, "velocity.geojson", fc)
        n = Counter(f["properties"]["trust"] for f in fc["features"])
        stats["velocity"] = {"strokes": len(fc["features"]), "trust": dict(n), "deepest_cm": round(st["deepest_cm"], 1),
                             "band_cm": cfg["velocity"].get("band_m", 0.02) * 100}
        # "headline": false in velocity{} = no "sure it sinks at least" number (e.g. features smaller than its 1.5 km window)
        stats["velocity"]["headline"] = headline(cfg, v) if cfg["velocity"].get("headline", True) else {}
        layers.append({"id": "velocity", "kind": "hachures-time", "file": "velocity.geojson",
                       "label": {"en": "Ground motion", "uk": "Рух ґрунту"},
                       "props": ["rate_cm_yr", "sigma_cm_yr", "snr", "segments", "trust", "dark", "sense", "t"]})
        if "contours" in L:
            cf = contours_layer(cfg, v, st["frames_s"], st["valid"] | _contour_mask(cfg, v, st), T, crs)
            _write(out, "contours.geojson", cf)
            layers.append({"id": "contours", "kind": "contours", "file": "contours.geojson",
                           "label": {"en": "Depth lines", "uk": "Лінії глибини"}, "props": ["level_cm", "major", "d"],
                           "unit": "cm", "step_cm": cfg.get("contours", {}).get("step_cm", 5)})
            if verbose:
                print(f"  contours: {len(cf['features'])} lines over {len(v.dates) - 1} dates")
        if "surface" in L:
            sf = surface_layer(cfg, v, st["frames_s"], st["valid"] | _contour_mask(cfg, v, st), T, crs)
            _write(out, "surface.json", sf)
            layers.append({"id": "surface", "kind": "surface", "file": "surface.json",
                           "label": {"en": "Cross-section", "uk": "Розріз"}, "visible": False})
            if verbose:
                print(f"  surface: {sf['nx']} x {sf['ny']} x {len(sf['frames'])} dates, gnss {[g['station'] for g in sf['gnss']]}")
        if verbose:
            print(f"  velocity: {len(fc['features'])} strokes, trust {dict(n)}, deepest {st['deepest_cm']:.0f} cm")
    if "grid" in L:
        g = grid_layer(cfg, v)
        g.update({"value": "rate_cm_yr", "series": "cum_cm", "series_unit": "cm",
                  "modes": [{"key": "sigma_cm_yr", "label": {"en": "error σ", "uk": "похибка σ"}, "ramp": "sigma"},
                            {"key": "segments", "label": {"en": "chain breaks", "uk": "розриви ланцюжка"}, "ramp": "segments"}]})
        _write(out, "grid.json", g)
        layers.append({"id": "grid", "kind": "cells", "file": "grid.json",
                       "label": {"en": "Trust cells", "uk": "Клітинки довіри"}, "props": g["columns"]})
        if verbose:
            print(f"  grid: {len(g['cells'])} cells of {cfg['grid_km']} km")
        if cfg.get("rhythm"):
            rh = rhythm_layer(cfg, g)
            _write(out, "rhythm.json", rh)
            layers.append({"id": "rhythm", "kind": "rhythm", "file": "rhythm.json",
                           "label": {"en": "Score of the districts", "uk": "Партитура районів"},
                           "visible": cfg["rhythm"].get("visible", True)})
            if verbose:
                print(f"  rhythm: {[((r["name"].get("en") if isinstance(r["name"], dict) else r["name"]), r['rate_cm_yr'], r['still']) for r in rh['rows']]}")
    if "fields" in L:
        g = fields_layer(cfg, k=cfg.get("fields", {}).get("k", 4))
        g.update({"series": "rel_coh", "series_unit": "rel. coherence", "default_mode": "rel_coh",
                  "modes": [{"key": "rel_coh", "label": {"en": "coherence on the date", "uk": "когерентність на дату"},
                             "ramp": "coherence", "time": True},
                            {"key": "cluster", "label": {"en": "field rhythm (4 types)", "uk": "ритм полів (4 типи)"},
                             "categorical": "dance"},
                            {"key": "mean_rel_coh", "label": {"en": "mean coherence", "uk": "середня когерентність"},
                             "ramp": "coherence"}]})
        _write(out, "fields.json", g)
        layers.append({"id": "fields", "kind": "cells", "file": "fields.json",
                       "label": {"en": "Rhythm of the fields", "uk": "Ритм полів"}, "props": g["columns"]})
        stats["fields"] = {"cells": len(g["cells"]), "clusters": g["cluster_sizes"]}
        dates = dates or g["dates"]
        n_pairs = n_pairs or len(g["dates"])
        if verbose:
            print(f"  fields: {len(g['cells'])} blocks, clusters {g['cluster_sizes']}")
    if "footprint" in L:
        fc, st = footprint_layer(cfg)
        _write(out, "footprint.geojson", fc)
        layers.append({"id": "footprint", "kind": "polygons", "file": "footprint.geojson", "style": "change",
                       "label": {"en": "Change footprint", "uk": "Слід змін"}, "props": ["area_km2", "pair"]})
        stats["footprint"] = {"polygons": len(fc["features"]),
                              "largest_km2": fc["features"][0]["properties"]["area_km2"] if fc["features"] else 0,
                              "event_date": st["event_date"], "false_alarm_check": st["false_alarm_check"]}
        dates = dates or st["dates"]
        n_pairs = n_pairs or len(st["dates"]) - 1
        if verbose:
            print(f"  footprint: {stats['footprint']}")
        if cfg["footprint"].get("series"):
            fc2, strip = change_series(cfg)
            _write(out, "changes.geojson", fc2)
            _write(out, "change_strip.json", strip)
            layers.append({"id": "changes", "kind": "polygons", "file": "changes.geojson", "style": "change_pair", "z": 2.1,
                           "label": {"en": "Lost radar fingerprint, pair by pair", "uk": "Втрачений відбиток, пара за парою"},
                           "props": ["d", "area_km2", "pair", "event"], "time": "d"})
            layers.append({"id": "change_strip", "kind": "rhythm", "file": "change_strip.json",
                           "label": {"en": "Change through the year", "uk": "Зміни протягом року"}, "visible": False})
            if verbose:
                print(f"  changes: {len(fc2['features'])} patches; largest per pair {strip['rows'][0]['bars']}")
    layers += context_layers(cfg, out, verbose)
    layers += relief_info(cfg, out, verbose)
    layers += landcover_layers(cfg, out, verbose)
    for extra in cfg.get("extra_layers", []):        # layers produced by other scripts (e.g. a quake card)
        layers.append(extra)
    if v is not None and getattr(v, "tropo", None):
        stats["tropo"] = v.tropo
    if v is not None and getattr(v, "highpass", None):
        stats["highpass"] = v.highpass
    story = write_story(cfg, out)
    gn = gnss_check(cfg, v) if v is not None else []
    w, s, e, n_ = cfg["bbox"]
    man = {
        "engine_version": ENGINE_VERSION, "id": cfg["id"], "title": cfg.get("title", cfg["id"]),
        "subtitle": cfg.get("subtitle", ""), "frame_sub": cfg.get("frame_sub"), "bbox": cfg["bbox"], "center": [(w + e) / 2, (s + n_) / 2],
        "dates": dates, "n_pairs": n_pairs, "direction": cfg["direction"],
        "reference": {"bbox": cfg.get("ref_bbox"), "label": cfg.get("ref_label", "")},
        "anchor": {"type": cfg["anchor"], "gnss": gn},
        "trust": cfg["trust"], "stats": stats, "layers": layers, "notes": cfg.get("notes", []),
        "generated": datetime.now().isoformat(timespec="seconds"),
    }
    if story:
        man["story"] = story
    (out / "manifest.json").write_text(json.dumps(man, indent=1, ensure_ascii=False), encoding="utf-8")
    return out / "manifest.json"


def write_story(cfg, out) -> str | None:
    """Copy the plate's pages (cfg['story'] = path to sites/<id>.story.json) next to the layers; checks layer ids loosely.

    Pages are written by people, so this only validates the shape and returns the file name for the manifest."""
    src = cfg.get("story")
    if not src:
        return None
    p = Path(src)
    if not p.is_absolute():
        p = Path(__file__).resolve().parents[2] / p
    story = json.loads(p.read_text(encoding="utf-8"))
    assert isinstance(story.get("pages"), list) and story["pages"], f"{p}: no pages"
    for i, pg in enumerate(story["pages"]):
        assert "title" in pg and "camera" in pg, f"{p}: page {i + 1} needs title and camera"
    _write(out, "story.json", story)
    return "story.json"


# ------------------------------------------------------------------ change layers from coherence (no unwrapping)

def chain_coherence(cfg):
    """(pairs as RawWindow list, granules) of the cached chain: coherence needs no phase processing."""
    chain = cached_chain(cfg)
    raws = [nisar.RawWindow.load(cache_dir(cfg) / f"{g.name}.npz") for g in chain]
    shape = raws[0].coherence.shape
    keep = [i for i, r in enumerate(raws) if r.coherence.shape == shape]
    return [raws[i] for i in keep], [chain[i] for i in keep]


def fields_layer(cfg, k=4, seed=0):
    """'Rhythm of the fields' (X-coh): per-block coherence through the year, relative to a stable class, + k-means 'dances'.

    cfg['fields'] = {"block_km": 2, "landcover": {"npz": ..., "grid_tif": ..., "crop": "cropland", "ref": "builtup",
                     "crop_min": 0.6, "ref_min": 0.3}}   (landcover optional: then all blocks, no normalisation)
    """
    import rasterio
    p = cfg.get("fields", {})
    raws, chain = chain_coherence(cfg)
    x, y, epsg = raws[0].x, raws[0].y, raws[0].epsg
    px = abs(x[1] - x[0])
    B = max(1, int(round(p.get("block_km", 2) * 1000 / px)))
    H, W = (len(y) // B) * B, (len(x) // B) * B

    def blk(a):
        with np.errstate(all="ignore"):
            return np.nanmedian(a[:H, :W].reshape(H // B, B, W // B, B).transpose(0, 2, 1, 3).reshape(H // B, W // B, -1), -1)
    import warnings
    warnings.simplefilter("ignore", RuntimeWarning)
    C = np.array([blk(r.coherence.astype(float)) for r in raws])           # (n_pairs, h, w)
    crop = np.ones(C.shape[1:], bool)
    norm = np.ones(len(raws))
    lc = p.get("landcover")
    if lc:
        z = np.load(ROOT / lc["npz"])
        with rasterio.open(ROOT / lc["grid_tif"]) as ds:
            tr, crs = ds.transform, ds.crs
        dst = gridding.transform_from_centres(x, y)
        frac = {key: blk(gridding.to_grid(z[lc[key]].astype(float), tr, crs, dst, f"EPSG:{epsg}", (len(y), len(x))))
                for key in ("crop", "ref")}
        crop = frac["crop"] > lc.get("crop_min", 0.6)
        ref = frac["ref"] > lc.get("ref_min", 0.3)
        norm = np.array([np.nanmedian(c[ref]) for c in C])
    R = C / norm[:, None, None]
    ok = crop & np.isfinite(R).all(0)
    Y = R[:, ok].T                                                         # blocks x pairs
    Z = Y - Y.mean(1, keepdims=True)
    rng = np.random.default_rng(seed)
    cen = Z[rng.choice(len(Z), k, replace=False)]
    for _ in range(60):
        lab = np.argmin(((Z[:, None] - cen[None]) ** 2).sum(-1), 1)
        cen = np.array([Z[lab == j].mean(0) if (lab == j).any() else cen[j] for j in range(k)])
    order = np.argsort([-(lab == j).sum() for j in range(k)])               # cluster 0 = most common dance
    relab = np.empty(k, int)
    relab[order] = np.arange(k)
    lab = relab[lab]
    xs, ys = x[:W].reshape(-1, B).mean(1), y[:H].reshape(-1, B).mean(1)
    XX, YY = np.meshgrid(xs, ys)
    lon, lat = warp_transform(f"EPSG:{epsg}", "EPSG:4326", XX[ok].tolist(), YY[ok].tolist())
    mids = [(g.ref_date + (g.sec_date - g.ref_date) / 2).date().isoformat() for g in chain]
    cells = [[round(lon[i], 5), round(lat[i], 5), int(lab[i]), round(float(Y[i].mean()), 3),
              [round(float(v), 3) for v in Y[i]]] for i in range(len(Y))]
    curves = [[round(float(v), 3) for v in Y[lab == j].mean(0)] for j in range(k)]
    return {"cell_km": B * px / 1000, "columns": ["lon", "lat", "cluster", "mean_rel_coh", "rel_coh"],
            "dates": mids, "pair_days": [round(g.days) for g in chain], "normalised_by": lc["ref"] if lc else None,
            "cluster_curves": curves, "cluster_sizes": [int((lab == j).sum()) for j in range(k)], "cells": cells}


def change_series(cfg):
    """The same footprint rule on every pair (vs the pair before it): where the ground lost its radar 'fingerprint', and when.

    Returns (GeoJSON patches with {d, area_km2, pair, event}, a rhythm-style strip {dates, unit, rows[{name, bars}]}).
    d is the index (in the scene dates) of the pair's second date, so the map shows a pair's patches while the slider sits
    in it. The pair that contains cfg['footprint']['event_date'] is marked event=True: the scar the plate is about.
    Everything here inherits the post-hoc caveat of the rule itself (found on this event, not yet tested blind)."""
    from rasterio import features
    from scipy import ndimage
    p = cfg["footprint"]
    ev = datetime.fromisoformat(p["event_date"])
    raws, chain = chain_coherence(cfg)
    x, y, epsg = raws[0].x, raws[0].y, raws[0].epsg
    T = gridding.transform_from_centres(x, y)
    px2 = abs(x[1] - x[0]) * abs(y[1] - y[0]) / 1e6
    feats, largest, total = [], [None], [None]
    for i in range(1, len(raws)):
        a, b = raws[i - 1].coherence.astype(float), raws[i].coherence.astype(float)
        m = ndimage.binary_opening((b / np.maximum(a, 1e-3) < p.get("ratio", 0.35)) & (a > p.get("pre_coh", 0.4)), iterations=1)
        lab, n = ndimage.label(m)
        sizes = np.bincount(lab.ravel()); sizes[0] = 0
        big = sizes * px2 >= p.get("min_km2", 0.05)
        largest.append(round(float(sizes.max() * px2), 3))
        total.append(round(float(sizes[big].sum() * px2), 3))
        is_ev = chain[i].ref_date < ev < chain[i].sec_date
        keep = np.isin(lab, np.nonzero(big)[0])
        for geom, val in features.shapes(lab.astype("int32"), mask=keep, transform=T):
            for k, ring in enumerate(geom["coordinates"]):
                r = np.array(ring)
                lon, lat = warp_transform(f"EPSG:{epsg}", "EPSG:4326", r[:, 0].tolist(), r[:, 1].tolist())
                geom["coordinates"][k] = np.round(np.column_stack([lon, lat]), 5).tolist()
            feats.append({"type": "Feature", "geometry": geom,
                          "properties": {"d": i + 1, "area_km2": round(float(sizes[int(val)] * px2), 3), "event": bool(is_ev),
                                         "pair": f"{chain[i].ref_date:%Y-%m-%d} → {chain[i].sec_date:%Y-%m-%d}"}})
    dates = [chain[0].ref_date.date().isoformat()] + [g.sec_date.date().isoformat() for g in chain]
    strip = {"dates": dates, "unit": "km²", "kind": "bars",
             "note": {"en": "area where the radar lost its grip in each 12-day pair, compared with the pair before",
                      "uk": "площа, де радар втратив «відбиток» землі в кожній 12-денній парі, порівняно з попередньою"},
             "rows": [{"name": {"en": "largest patch", "uk": "найбільша пляма"}, "bars": largest},
                      {"name": {"en": "all patches", "uk": "усі плями разом"}, "bars": total}],
             "event_date": p["event_date"],
             "rule": "co / previous-pair coherence < ratio where the previous pair > pre_coh; patches >= min_km2 (post-hoc rule)"}
    return {"type": "FeatureCollection", "features": feats}, strip


def footprint_layer(cfg):
    """Change footprint (L-dpm): co-event coherence / adjacent pre-event coherence < ratio where the pre pair was coherent.

    cfg['footprint'] = {"event_date": "2026-08-26", "ratio": 0.35, "pre_coh": 0.4, "min_km2": 0.05}
    Returns GeoJSON polygons (largest first) with area and pair dates, plus the same rule on two pre-event
    pairs as a false-alarm check.
    """
    from rasterio import features
    from scipy import ndimage
    p = cfg["footprint"]
    ev = datetime.fromisoformat(p["event_date"])
    raws, chain = chain_coherence(cfg)
    co = next(i for i, g in enumerate(chain) if g.ref_date < ev < g.sec_date)
    x, y, epsg = raws[0].x, raws[0].y, raws[0].epsg
    T = gridding.transform_from_centres(x, y)
    px2 = abs(x[1] - x[0]) * abs(y[1] - y[0]) / 1e6

    def mask(i):
        a, b = raws[i - 1].coherence.astype(float), raws[i].coherence.astype(float)
        m = (b / np.maximum(a, 1e-3) < p.get("ratio", 0.35)) & (a > p.get("pre_coh", 0.4))
        return ndimage.binary_opening(m, iterations=1)
    m = mask(co)
    lab, n = ndimage.label(m)
    sizes = np.bincount(lab.ravel())
    sizes[0] = 0
    keep = np.isin(lab, np.nonzero(sizes * px2 >= p.get("min_km2", 0.05))[0])
    feats = []
    for geom, val in features.shapes(lab.astype("int32"), mask=keep, transform=T):
        ring = np.array(geom["coordinates"][0])
        lon, lat = warp_transform(f"EPSG:{epsg}", "EPSG:4326", ring[:, 0].tolist(), ring[:, 1].tolist())
        feats.append({"type": "Feature", "properties": {"area_km2": round(float(sizes[int(val)] * px2), 2),
                                                         "pair": f"{chain[co].ref_date:%Y-%m-%d} → {chain[co].sec_date:%Y-%m-%d}"},
                      "geometry": {"type": "Polygon", "coordinates": [np.round(np.column_stack([lon, lat]), 5).tolist()]}})
    feats.sort(key=lambda f: -f["properties"]["area_km2"])
    checks = []
    for i in (co - 1, co - 2):
        if i >= 1:
            s = np.bincount(ndimage.label(mask(i))[0].ravel())
            s[0] = 0
            checks.append({"pair": f"{chain[i].ref_date:%Y-%m-%d} → {chain[i].sec_date:%Y-%m-%d}",
                           "largest_km2": round(float(s.max() * px2), 3)})
    return {"type": "FeatureCollection", "features": feats}, {"event_date": p["event_date"], "false_alarm_check": checks,
                                                               "dates": [chain[0].ref_date.date().isoformat()] +
                                                                        [g.sec_date.date().isoformat() for g in chain]}
