"""Animated hachure map: one year of NISAR subsidence growing as red ink (MP4 + GIF).

Input: data/processed/<name>/nisar_series.npz from scripts/nisar_gunw_rate.py.

Idea: stroke GEOMETRY is fixed (hachures of the final cumulative displacement, so frames
don't flicker with noise); TIME drives the ink. Each stroke grows from its top end, and its
drawn length is the share of the year's total sinking already done at the stroke's middle
(cumulative now / cumulative at the last date). A pause or speed-up in sinking shows as
strokes stalling or jumping.

    uv run python scripts/animate_subsidence.py                  # Mexico City
    uv run python scripts/animate_subsidence.py --name fake-test --label FAKE
"""
from __future__ import annotations

import argparse
import time
from datetime import datetime
from pathlib import Path

import numpy as np
from scipy.ndimage import (distance_transform_edt, gaussian_filter, gaussian_filter1d,
                           map_coordinates, median_filter, zoom)

from hachure import dem, field, gridding, trace
from hachure.nisar import fit_ramp

ROOT = Path(__file__).resolve().parents[1]
PAPER, GREY, RED = "#f3eedf", "#5f594d", "#a8322a"


def smooth_nan(a, sig, valid):
    w = gaussian_filter(valid.astype(float), sig)
    return gaussian_filter(np.where(valid, a, 0.0), sig) / np.maximum(w, 1e-6), w


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", default="mexico-city")
    ap.add_argument("--bbox", nargs=4, type=float, default=[-99.32, 19.17, -98.82, 19.58])
    ap.add_argument("--res", type=float, default=50.0)
    ap.add_argument("--band", type=float, default=0.02, help="m of subsidence per ink band")
    ap.add_argument("--coh-min", type=float, default=0.45)
    ap.add_argument("--noise-floor", type=float, default=0.05, help="m/yr; slower = not drawn")
    ap.add_argument("--sigma-dir", type=float, default=500.0)
    ap.add_argument("--sigma-tone", type=float, default=250.0)
    ap.add_argument("--sub", type=int, default=4, help="interpolated frames between dates")
    ap.add_argument("--fps", type=int, default=12)
    ap.add_argument("--hold", type=int, default=24, help="extra frames on the last date")
    ap.add_argument("--scale", type=float, default=1.0)
    ap.add_argument("--label", default="NISAR L-band, descending track")
    ap.add_argument("--light", choices=["static", "dynamic"], default="static",
                    help="dynamic = stroke weight re-lit every frame from the current depth")
    ap.add_argument("--shadow-km", type=float, default=0.0,
                    help="cast shadows: length (km) of the shadow thrown by the deepest point at "
                         "the last date; sets the vertical exaggeration. 0 = off")
    ap.add_argument("--sun-alt", type=float, default=30.0, help="sun altitude for cast shadows")
    ap.add_argument("--shadow-soft-cm", type=float, default=2.0,
                    help="blockers rising less than this (real cm) above the ray give partial shade")
    ap.add_argument("--terrain-shadow", action="store_true",
                    help="mountains cast real (unexaggerated) shadows on the valley, same sun")
    ap.add_argument("--terrain-soft-m", type=float, default=40.0,
                    help="soft shadow edge for the terrain, metres of blocker height")
    ap.add_argument("--hatch-spacing", type=int, default=2, help="px between shadow lines")
    ap.add_argument("--stills", nargs="+", type=float,
                    help="only render PNG stills at these fractions of the year (0..1), no video")
    ap.add_argument("--tag", default="", help="suffix for output file names")
    args = ap.parse_args()

    out = ROOT / "data" / "processed" / args.name
    z_raw = np.load(out / "nisar_series.npz")
    cum, coh_src = z_raw["cum"].astype("float64"), z_raw["coherence"]
    dates = [datetime.fromisoformat(str(d)) for d in z_raw["dates"]]
    src_t = gridding.transform_from_centres(z_raw["x"], z_raw["y"])
    src_crs = f"EPSG:{int(z_raw['epsg'])}"
    years = (dates[-1] - dates[0]).days / 365.25
    print(f"series: {len(dates)} dates, {dates[0]:%Y-%m-%d} -> {dates[-1]:%Y-%m-%d}")

    # ---- DEM grid + static grey terrain
    raw = ROOT / "data" / "raw" / f"copdem_mexico-city_{int(args.res)}m.tif"
    dem.fetch_copdem(tuple(args.bbox), raw, res=args.res)
    zdem, T, crs = dem.load(raw)
    shape = zdem.shape
    t0 = time.time()
    terrain_cast = None
    if args.terrain_shadow:
        terrain_cast = field.cast_shadow(zdem, args.res, azimuth=315.0, altitude=args.sun_alt,
                                         soft_m=args.terrain_soft_m)
        print(f"terrain shadows: {terrain_cast.mean() * 100:.0f}% of the map shaded "
              f"({time.time() - t0:.0f}s)")
    terrain = trace.hachures(field.prepare(zdem, args.res, sigma_dir_m=150, sigma_tone_m=80,
                                           extra_shadow=terrain_cast),
                             interval_m=60, spacing_px=1.5)
    print(f"terrain {len(terrain)} strokes ({time.time() - t0:.0f}s)")

    # ---- series onto the DEM grid, temporal + spatial denoising, per-frame deramp
    cum = gaussian_filter1d(cum, 1.0, axis=0, mode="nearest")       # atmosphere flicker
    cum -= cum[0]
    coh = gridding.to_grid(coh_src, src_t, src_crs, T, crs, shape)
    coherent = np.nan_to_num(coh) >= args.coh_min
    sig = 400 / args.res
    frames_d, frames_s = [], []
    for k in range(len(dates)):
        g = gridding.to_grid(cum[k], src_t, src_crs, T, crs, shape)
        ok = coherent & np.isfinite(g)
        g = median_filter(np.where(ok, g, np.nanmedian(g[ok]) if ok.any() else 0), size=5)
        s, w = smooth_nan(g, sig, ok)
        if k > 0:
            yrs = (dates[k] - dates[0]).days / 365.25
            plane = fit_ramp(np.where(ok & (w > 0.3), s, np.nan), cut=-max(0.005, 0.03 * yrs))
            g, s = g - plane, s - plane
        frames_d.append(np.where(ok, g, np.nan))
        frames_s.append(s)
    final, final_s = frames_d[-1], frames_s[-1]
    support = smooth_nan(np.ones(shape), sig, coherent)[1] > 0.3
    sinking = (final_s <= -args.noise_floor * years) & support
    valid = coherent & sinking
    print(f"final: 1st pct {np.nanpercentile(final[valid], 1) * 100:.1f} cm over {years:.2f} yr")

    # ---- strokes on the FINAL displacement
    f_sub, interval = field.prepare_scalar(np.where(np.isfinite(final), final, 0.0), args.res,
                                           band=args.band, valid=valid,
                                           sigma_dir_m=args.sigma_dir,
                                           sigma_tone_m=args.sigma_tone)
    strokes = trace.hachures(f_sub, interval_m=interval, spacing_px=1.5, max_len_px=30,
                             min_slope_deg=1.0, min_dark=0.12)
    print(f"subsidence {len(strokes)} strokes")
    fs = np.where(np.isfinite(final_s), final_s, 0.0)

    def at(a, p):
        return map_coordinates(a, [p[:, 1], p[:, 0]], order=1, mode="nearest")

    mids = np.array([at(fs, s.pts[len(s.pts) // 2: len(s.pts) // 2 + 1])[0] for s in strokes])
    mid_pts = np.array([s.pts[len(s.pts) // 2] for s in strokes])
    cum_len = [np.concatenate([[0], np.cumsum(np.hypot(*np.diff(s.pts, axis=0).T))]) for s in strokes]
    all_pts = np.vstack([s.pts for s in strokes])
    owner = np.repeat(np.arange(len(strokes)), [len(s.pts) for s in strokes])
    n_pts = np.bincount(owner, minlength=len(strokes))

    def per_stroke_mean(a):
        return np.bincount(owner, weights=at(a, all_pts), minlength=len(strokes)) / n_pts

    # ---- light: same field scaling as the strokes, so a re-lit frame matches the static map
    scale = interval / args.band
    res = args.res
    depth_final = abs(float(np.nanpercentile(final_s[valid], 1)))
    tan_sun = np.tan(np.radians(args.sun_alt))
    exag = args.shadow_km * 1000 * tan_sun / max(depth_final, 1e-3)
    fill_idx = distance_transform_edt(~support, return_distances=False, return_indices=True)
    if args.shadow_km > 0:
        print(f"cast shadows: sun {args.sun_alt:.0f} deg, vertical exaggeration x{exag:.0f}")

    def lighting(k, a):
        """-> (stroke widths, cast-shadow raster or None) for time step (k, a)."""
        s_t = (1 - a) * frames_s[k] + a * frames_s[k + 1]
        cast = None
        if args.shadow_km > 0:
            # outside the radar's support take the nearest valid depth (no fake cliffs at the
            # mask edge), and ray-march on a 2x coarser grid (shadows are smooth anyway)
            surf = np.minimum(s_t, 0.0)[tuple(fill_idx)] * exag
            c = field.cast_shadow(surf[::2, ::2], res * 2, azimuth=315.0, altitude=args.sun_alt,
                                  max_dist_m=args.shadow_km * 1300,
                                  soft_m=args.shadow_soft_cm / 100 * exag)
            cast = zoom(c, 2, order=1)[: shape[0], : shape[1]]
            cast = np.pad(cast, [(0, shape[0] - cast.shape[0]), (0, shape[1] - cast.shape[1])],
                          mode="edge") * support
        if terrain_cast is not None:
            cast = terrain_cast if cast is None else np.maximum(cast, terrain_cast)
        if args.light == "static":
            return widths, cast
        d_t = (1 - a) * frames_d[k] + a * frames_d[k + 1]
        f_t = field.prepare(d_t, res, sigma_dir_m=args.sigma_dir, sigma_tone_m=args.sigma_tone,
                            lit_floor=0.45, z_factor=scale, valid=valid, extra_shadow=cast)
        return 0.25 + 2.2 * per_stroke_mean(f_t.dark) * 72 / dpi, cast

    def hatch(cast, chunk=4, thr=0.08):
        """Engraver's shadow: thin horizontal lines whose weight follows the shade."""
        rows = np.arange(1, cast.shape[0], args.hatch_spacing)
        nc = cast.shape[1] // chunk
        v = cast[rows, : nc * chunk].reshape(len(rows), nc, chunk).mean(-1)
        ri, ci = np.nonzero(v > thr)
        y, x0 = rows[ri].astype(float), ci * chunk
        segs = np.stack([np.stack([x0, y], 1), np.stack([x0 + chunk, y], 1)], 1)
        return list(segs), list((0.15 + 0.75 * v[ri, ci]) * 72 / dpi)

    # ---- render
    import imageio.v2 as imageio
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection

    h, w = shape
    dpi = 100
    fig = plt.figure(figsize=(w * args.scale / dpi, (h * args.scale + 70) / dpi), dpi=dpi)
    ax = fig.add_axes([0, 0, 1, h * args.scale / (h * args.scale + 70)])
    fig.patch.set_facecolor(PAPER)
    ax.set_facecolor(PAPER)
    ax.set_xlim(0, w)
    ax.set_ylim(h, 0)
    ax.axis("off")
    tw = [(0.12 + 1.9 * s.dark) * 72 / dpi for s in terrain]
    ax.add_collection(LineCollection([s.pts for s in terrain], linewidths=tw, colors=GREY))
    shade = LineCollection([], colors=GREY, alpha=0.85, capstyle="butt")
    ax.add_collection(shade)
    red = LineCollection([], colors=RED, capstyle="round")
    ax.add_collection(red)
    widths = np.array([(0.25 + 2.2 * s.dark) * 72 / dpi for s in strokes])
    title = fig.text(0.02, 1 - 30 / (h * args.scale + 70), "", fontsize=17, family="serif",
                     color="#1b1a17", va="center")
    sub = fig.text(0.02, 1 - 56 / (h * args.scale + 70), "", fontsize=11, family="serif",
                   color="#5f594d", va="center")
    bar_ax = fig.add_axes([0.62, 1 - 50 / (h * args.scale + 70), 0.35, 12 / (h * args.scale + 70)])
    bar_ax.axis("off")

    # frame schedule: sub-steps between dates, then hold the last
    sched = []
    for k in range(len(dates) - 1):
        for j in range(args.sub):
            sched.append((k, j / args.sub))
    sched += [(len(dates) - 2, 1.0)] * (args.hold + 1)

    mid_frames = [at(np.where(np.isfinite(s), s, 0.0), mid_pts) for s in frames_s]
    # stroke length = share of the year's sinking already done at the stroke's middle
    # (an earlier "draw when the ground reaches the stroke's ring depth" rule left the map
    # blank until the last weeks, because every ring is only reached at the very end)
    final_mid = np.minimum(mids, -1e-4)

    def draw(k, a):
        c_mid = (1 - a) * mid_frames[k] + a * mid_frames[k + 1]
        frac = np.clip(c_mid / final_mid, 0, 1)
        w_now, cast = lighting(k, a)
        if cast is not None:
            hs, hw = hatch(cast)
            shade.set_segments(hs)
            shade.set_linewidths(hw if hw else [0])
        segs, lw = [], []
        for i in np.flatnonzero(frac > 0.02):
            L = cum_len[i]
            cut = frac[i] * L[-1]
            m = np.searchsorted(L, cut)
            p = strokes[i].pts[: max(m, 1) + 1]
            if m < len(L) and m > 0:
                tt = (cut - L[m - 1]) / max(L[m] - L[m - 1], 1e-9)
                p = np.vstack([strokes[i].pts[:m], strokes[i].pts[m - 1] + tt * (strokes[i].pts[m] - strokes[i].pts[m - 1])])
            segs.append(p)
            lw.append(w_now[i])
        red.set_segments(segs)
        red.set_linewidths(lw)
        d = dates[k] + (dates[k + 1] - dates[k]) * a
        depth = (1 - a) * frames_s[k] + a * frames_s[k + 1]
        worst = np.nanpercentile(depth[valid], 1) * 100 if valid.any() else 0
        title.set_text(f"Mexico City is sinking  ·  {d:%d %b %Y}")
        sub.set_text(f"deepest so far: {abs(worst):.0f} cm   ·   each red ring = {args.band * 100:.0f} cm"
                     f"   ·   {args.label}")
        bar_ax.clear()
        bar_ax.axis("off")
        prog = (d - dates[0]).total_seconds() / (dates[-1] - dates[0]).total_seconds()
        bar_ax.barh([0], [1], color="#d9d2c0", height=1)
        bar_ax.barh([0], [prog], color=RED, height=1)
        bar_ax.set_xlim(0, 1)
        fig.canvas.draw()
        return np.asarray(fig.canvas.buffer_rgba())[..., :3].copy()

    tag = f"_{args.tag}" if args.tag else ""
    if args.stills:
        span = (dates[-1] - dates[0]).total_seconds()
        secs = [(d - dates[0]).total_seconds() for d in dates]
        for q in args.stills:
            t = q * span
            k = min(int(np.searchsorted(secs, t, side="right")) - 1, len(dates) - 2)
            a = float(np.clip((t - secs[k]) / (secs[k + 1] - secs[k]), 0, 1))
            t0 = time.time()
            p = out / f"subsidence_still{tag}_{int(round(q * 100)):03d}.png"
            imageio.imwrite(p, draw(k, a))
            print(f"{p}  ({time.time() - t0:.1f}s)")
        plt.close(fig)
        return

    mp4 = out / f"subsidence_animation{tag}.mp4"
    writer = imageio.get_writer(mp4, fps=args.fps, codec="libx264", quality=8,
                                macro_block_size=8)
    gif_frames = []
    t0 = time.time()
    for n, (k, a) in enumerate(sched):
        img = draw(k, a)
        writer.append_data(img)
        if n % 2 == 0 or n == len(sched) - 1:
            gif_frames.append(img[::2, ::2])
        if n == len(sched) - 1:
            imageio.imwrite(out / f"subsidence_animation{tag}_last.png", img)
    writer.close()
    gif = out / f"subsidence_animation{tag}.gif"
    imageio.mimsave(gif, gif_frames, duration=2 / args.fps, loop=0)
    plt.close(fig)
    print(f"{len(sched)} frames ({time.time() - t0:.0f}s)")
    print(mp4)
    print(gif)


if __name__ == "__main__":
    main()
