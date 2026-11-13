"""Two-ink map: grey terrain hachures + red subsidence hachures (Mexico City).

The subsidence field is SYNTHETIC (src/hachure/synthetic.py) until we stream
real NISAR GUNW unwrapped phase. Strokes in the red layer run "downhill" on the
displacement surface, i.e. they point to where the ground sinks fastest; each red
band = `--band` metres of subsidence (default 2 cm/yr), like a meaningful fringe.

    uv run python scripts/hachure_subsidence.py
"""
from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np

from hachure import dem, export, field, synthetic, trace

ROOT = Path(__file__).resolve().parents[1]


def _rate_on_grid(path, transform, crs, shape):
    """Reproject (rate, coherence) bands of a nisar_rate.tif onto the DEM grid."""
    import rasterio
    from rasterio.warp import Resampling, reproject
    out = np.full((2, *shape), np.nan, dtype="float32")
    with rasterio.open(path) as ds:
        for b in (1, 2):
            reproject(rasterio.band(ds, b), out[b - 1], dst_transform=transform, dst_crs=crs,
                      resampling=Resampling.bilinear, dst_nodata=np.nan)
    return out[0].astype("float64"), np.nan_to_num(out[1]).astype("float64")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bbox", nargs=4, type=float, default=[-99.32, 19.17, -98.82, 19.58])
    ap.add_argument("--res", type=float, default=50.0)
    ap.add_argument("--band", type=float, default=0.02, help="m/yr per red band")
    ap.add_argument("--coh-min", type=float, default=0.35)
    ap.add_argument("--noise-floor", type=float, default=0.01, help="m/yr; weaker = not drawn")
    ap.add_argument("--scale", type=float, default=1.3)
    ap.add_argument("--rate", help="nisar_rate.tif from scripts/nisar_gunw_rate.py "
                                   "(default: SYNTHETIC field)")
    ap.add_argument("--label", default="NISAR", help="source label written into outputs")
    ap.add_argument("--sigma-dir", type=float, default=300.0,
                    help="generalization of stroke direction, m (raise for noisy/few pairs)")
    ap.add_argument("--sigma-tone", type=float, default=150.0)
    ap.add_argument("--no-deramp", action="store_true",
                    help="keep the planar trend of a real rate map (default: removed robustly)")
    ap.add_argument("--median", type=int, default=0,
                    help="median filter size (px) against speckle outliers; 5 is good for real data")
    ap.add_argument("--sink-only", action="store_true",
                    help="draw only where smoothed rate <= -noise_floor (subsidence map)")
    args = ap.parse_args()

    name = "mexico-city"
    raw = ROOT / "data" / "raw" / f"copdem_{name}_{int(args.res)}m.tif"
    out = ROOT / "data" / "processed" / name
    out.mkdir(parents=True, exist_ok=True)

    dem.fetch_copdem(tuple(args.bbox), raw, res=args.res)
    z, transform, crs = dem.load(raw)
    print(f"grid {z.shape}, z {z.min():.0f}..{z.max():.0f} m")

    t = time.time()
    f_terrain = field.prepare(z, args.res, sigma_dir_m=150, sigma_tone_m=80)
    terrain = trace.hachures(f_terrain, interval_m=60, spacing_px=1.5)
    print(f"terrain: {len(terrain)} strokes ({time.time() - t:.1f}s)")

    if args.rate:
        rate, coh = _rate_on_grid(args.rate, transform, crs, z.shape)
        label = args.label
        suffix = "_" + label.lower().replace(" ", "-")
    else:
        rate, coh = synthetic.subsidence_rate(transform, crs, z.shape, z)
        label, suffix = "SYNTHETIC", ""
    coherent = (coh >= args.coh_min) & np.isfinite(rate)
    from scipy.ndimage import gaussian_filter, median_filter
    if args.median:
        # NaN-aware enough: fill invalid with the median of valid pixels before filtering
        fillv = np.nanmedian(rate[coherent]) if coherent.any() else 0.0
        rate = median_filter(np.where(coherent, rate, fillv), size=args.median)
    # only draw signal above the noise floor (smoothed so the mask edge is clean);
    # normalized convolution so masked pixels don't pull the average toward 0
    sig = 400 / args.res
    wsum = gaussian_filter(coherent.astype(float), sig)
    smooth = gaussian_filter(np.where(coherent, rate, 0.0), sig) / np.maximum(wsum, 1e-6)
    if args.rate and not args.no_deramp:
        from hachure.nisar import fit_ramp
        plane = fit_ramp(np.where(coherent & (wsum > 0.3), smooth, np.nan))
        rate, smooth = rate - plane, smooth - plane
        print(f"ramp removed: {np.ptp(plane) * 100:.1f} cm/yr across the window")
    signal =(smooth <= -args.noise_floor) if args.sink_only else (np.abs(smooth) >= args.noise_floor)
    signal &= wsum > 0.3
    valid = coherent & signal
    rate = np.where(np.isfinite(rate), rate, 0.0)
    print(f"{label} rate {np.nanmin(rate) * 100:.1f}..{np.nanmax(rate) * 100:.1f} cm/yr, "
          f"{(~coherent).mean() * 100:.0f}% low coherence, {(~signal).mean() * 100:.0f}% below noise floor")
    t = time.time()
    f_sub, interval = field.prepare_scalar(rate, args.res, band=args.band, valid=valid,
                                           sigma_dir_m=args.sigma_dir,
                                           sigma_tone_m=args.sigma_tone)
    sinking = trace.hachures(f_sub, interval_m=interval, spacing_px=1.5, max_len_px=30,
                             min_slope_deg=1.0, min_dark=0.12)
    print(f"subsidence: {len(sinking)} strokes ({time.time() - t:.1f}s)")

    layers = [
        export.Layer(terrain, ink="#5f594d", w_min=0.12, w_max=2.0, name="terrain"),
        export.Layer(sinking, ink="#a8322a", w_min=0.25, w_max=2.4, name="subsidence"),
    ]
    print(export.to_png(layers, z.shape, out / f"subsidence{suffix}.png", scale=args.scale))
    print(export.to_svg(layers, z.shape, out / f"subsidence{suffix}.svg", scale=args.scale))
    print(export.to_geojson(terrain, transform, crs, out / "terrain.geojson"))
    print(export.to_geojson(sinking, transform, crs, out / f"subsidence{suffix}.geojson",
                            extra={"source": label}))


if __name__ == "__main__":
    main()
