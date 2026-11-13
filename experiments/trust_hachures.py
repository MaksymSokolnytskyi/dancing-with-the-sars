"""Prototype: hachures that show how much to trust them (morning pass, one still).

Same field as hachure_subsidence.py; each red stroke gets the signal-to-noise of the rate under it
(|rate| / 1-sigma from nisar_rate_sigma.tif): SNR >= 2 -> normal red; 1..2 -> pale and thin;
< 1 -> not drawn. Writes the plain and the trust version side by side for comparison.

    uv run python scripts/trust_hachures.py
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
from scipy.ndimage import gaussian_filter, map_coordinates, median_filter

import sys
sys.path[:0] = [str(Path(__file__).parent), str(Path(__file__).parents[1] / "scripts")]
from hachure_subsidence import _rate_on_grid  # noqa: E402

from hachure import dem, export, field, trace
from hachure.nisar import fit_ramp

ROOT = Path(__file__).resolve().parents[1]
BBOX = (-99.32, 19.17, -98.82, 19.58)
RES = 50.0
ASC = ROOT / "data" / "processed" / "mexico-city-asc"
OUT = ROOT / "data" / "processed" / "mexico-city"


def main():
    z, transform, crs = dem.load(ROOT / "data" / "raw" / f"copdem_mexico-city_{int(RES)}m.tif")
    terrain = trace.hachures(field.prepare(z, RES, sigma_dir_m=150, sigma_tone_m=80), interval_m=60, spacing_px=1.5)
    rate, coh = _rate_on_grid(ASC / "nisar_rate.tif", transform, crs, z.shape)
    import rasterio
    from rasterio.warp import Resampling, reproject
    sig = np.full(z.shape, np.nan, dtype="float32")
    with rasterio.open(ASC / "nisar_rate_sigma.tif") as ds:
        reproject(rasterio.band(ds, 2), sig, dst_transform=transform, dst_crs=crs, resampling=Resampling.bilinear,
                  dst_nodata=np.nan)
    sig = sig.astype(float)
    coherent = (coh >= 0.35) & np.isfinite(rate)
    rate = median_filter(np.where(coherent, rate, np.nanmedian(rate[coherent])), size=5)
    s = 400 / RES
    wsum = gaussian_filter(coherent.astype(float), s)
    smooth = gaussian_filter(np.where(coherent, rate, 0.0), s) / np.maximum(wsum, 1e-6)
    plane = fit_ramp(np.where(coherent & (wsum > 0.3), smooth, np.nan))
    rate, smooth = rate - plane, smooth - plane
    valid = coherent & (smooth <= -0.01) & (wsum > 0.3)
    f_sub, interval = field.prepare_scalar(np.where(np.isfinite(rate), rate, 0.0), RES, band=0.02, valid=valid,
                                           sigma_dir_m=300, sigma_tone_m=150)
    strokes = trace.hachures(f_sub, interval_m=interval, spacing_px=1.5, max_len_px=30, min_slope_deg=1.0,
                             min_dark=0.12)
    # SNR per stroke: smoothed |rate| over smoothed sigma at the stroke midpoint
    sig_s = gaussian_filter(np.nan_to_num(sig, nan=np.nanmedian(sig)), s)
    snr_map = np.abs(smooth) / np.maximum(sig_s, 1e-4)
    sure, unsure, dropped = [], [], 0
    for st in strokes:
        c, r = st.pts[len(st.pts) // 2]
        snr = float(map_coordinates(snr_map, [[r], [c]], order=1)[0])
        if snr >= 2:
            sure.append(st)
        elif snr >= 1:
            unsure.append(trace.Stroke(st.pts, st.dark * 0.6))
        else:
            dropped += 1
    print(f"{len(strokes)} red strokes: SNR>=2 {len(sure)}, 1-2 {len(unsure)}, <1 dropped {dropped}; "
          f"median sigma {np.nanmedian(sig[valid]) * 100:.1f} cm/yr where drawn")
    plain = [export.Layer(terrain, ink="#5f594d", w_min=0.12, w_max=2.0, name="terrain"),
             export.Layer(strokes, ink="#a8322a", w_min=0.25, w_max=2.4, name="subsidence")]
    trust = [export.Layer(terrain, ink="#5f594d", w_min=0.12, w_max=2.0, name="terrain"),
             export.Layer(unsure, ink="#d9a49c", w_min=0.15, w_max=1.0, name="subsidence_unsure"),
             export.Layer(sure, ink="#a8322a", w_min=0.25, w_max=2.4, name="subsidence_sure")]
    a = export.to_png(plain, z.shape, OUT / "trust_plain_asc.png", scale=1.0)
    b = export.to_png(trust, z.shape, OUT / "trust_faded_asc.png", scale=1.0)
    from PIL import Image
    ia, ib = Image.open(a), Image.open(b)
    both = Image.new("RGB", (ia.width + ib.width + 20, ia.height), "white")
    both.paste(ia, (0, 0))
    both.paste(ib, (ia.width + 20, 0))
    both.save(OUT / "trust_compare_asc.png")
    print(OUT / "trust_compare_asc.png")


if __name__ == "__main__":
    main()
