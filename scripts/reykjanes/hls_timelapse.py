"""SWIR time-lapse of the Sundhnukur crater row, Nov 2023 - Aug 2025, from NASA HLS (Sentinel-2 + Landsat).

All HLSS30 / HLSL30 v2.0 granules of tile T27VVL with scene cloud cover <= 50 % are checked; for each,
only the study window is read (cloud-optimized GeoTIFF range requests via a short-lived signed URL from
LP DAAC, so the Earthdata token is sent only to data.lpdaac.earthdatacloud.nasa.gov). A date is kept if
>= 85 % of the window is clear (Fmask: no cloud, adjacent cloud, shadow or snow) and the sun is high
enough (the window is not dark). SWIR false colour: B12, B11, B4 (Landsat: B07, B06, B04).

Writes output/reykjanes/ml/hls_timelapse.png (grid of the clear dates) and hls_timelapse.json.
Usage: python scripts/reykjanes/hls_timelapse.py
"""
from _secrets import read_env
import json
import subprocess
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import rasterio
from pyproj import Transformer
from rasterio.windows import from_bounds

from progress import Progress
from reykjanes_analysis import PLACES

ROOT = Path(__file__).resolve().parents[2]
ML = ROOT / "output" / "reykjanes" / "ml"
CACHE = ROOT / "data" / "hls" / "timelapse"
BOX = (-8, 8, -7, 7)
BANDS = {"S30": ("B12", "B11", "B04", "Fmask"), "L30": ("B07", "B06", "B04", "Fmask")}
ERUPTIONS = [("2023-12-18", "2023-12-21"), ("2024-01-14", "2024-01-16"), ("2024-02-08", "2024-02-09"), ("2024-03-16", "2024-05-09"),
             ("2024-05-29", "2024-06-22"), ("2024-08-22", "2024-09-05"), ("2024-11-20", "2024-12-09"), ("2025-04-01", "2025-04-01"),
             ("2025-07-16", "2025-08-05")]


def cmr():
    out = []
    for c in ("HLSS30", "HLSL30"):
        page = 1
        while True:
            r = subprocess.run(["curl", "-s", f"https://cmr.earthdata.nasa.gov/search/granules.json?short_name={c}&version=2.0"
                                f"&point=-22.40,63.87&temporal=2023-11-01T00:00:00Z,2025-08-31T23:59:59Z&cloud_cover=0,50"
                                f"&page_size=200&page_num={page}"], capture_output=True, text=True)
            e = json.loads(r.stdout)["feed"]["entry"]
            out += [x["producer_granule_id"] for x in e]
            if len(e) < 200:
                break
            page += 1
    return sorted(out, key=lambda g: g.split(".")[3])


def signed(tok, url):
    r = subprocess.run(["curl", "-s", "-o", "/dev/null", "-w", "%{redirect_url}", "-H", f"Authorization: Bearer {tok}", url],
                       capture_output=True, text=True)
    return r.stdout.strip()


def window(url, sx, sy):
    with rasterio.open(url) as r:
        w = from_bounds(sx + BOX[0] * 1000, sy + BOX[2] * 1000, sx + BOX[1] * 1000, sy + BOX[3] * 1000, r.transform).round_offsets().round_lengths()
        return r.read(1, window=w)


def main():
    env = read_env()
    tok = env["EARTHDATA_API_KEY"].strip().strip('"')
    sx, sy = Transformer.from_crs(4326, 32627, always_xy=True).transform(*PLACES["Svartsengi"])
    CACHE.mkdir(parents=True, exist_ok=True)
    grans = cmr()
    print(len(grans), "granules with scene cloud <= 50 %", flush=True)
    kept = []
    with Progress("HLS: перевірка дат (Fmask у вікні)", total=len(grans), unit="знімок") as pr, \
            rasterio.Env(GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR", CPL_VSIL_CURL_ALLOWED_EXTENSIONS=".tif"):
        for g in grans:
            kind = g.split(".")[1]
            cf = CACHE / f"{g}.npz"
            base = f"https://data.lpdaac.earthdatacloud.nasa.gov/lp-prod-protected/HLS{kind}.020/{g}.v2.0/{g}.v2.0"
            try:
                if cf.exists():
                    z = dict(np.load(cf))
                else:
                    fm = window(signed(tok, f"{base}.Fmask.tif"), sx, sy).astype(int)
                    bad = ((fm >> 1) & 1) | ((fm >> 2) & 1) | ((fm >> 3) & 1) | ((fm >> 4) & 1) | (fm == 255)
                    z = dict(clear=np.array(1 - bad.mean()))
                    if z["clear"] >= 0.85:
                        for b, nm in zip(BANDS[kind][:3], ("swir2", "swir1", "red")):
                            z[nm] = window(signed(tok, f"{base}.{b}.tif"), sx, sy).astype(np.int16)
                    np.savez_compressed(cf, **z)
                if float(z["clear"]) >= 0.85 and "red" in z:
                    red = z["red"].astype(float) * 1e-4
                    if np.nanmedian(np.where(z["red"] == -9999, np.nan, red)) > 0.02:       # skip polar-night darkness
                        kept.append((g, z))
                pr.update(1, note=f"{len(kept)} чистих дат")
            except Exception as ex:                                       # a broken granule: skip it, keep going
                print("skip", g, type(ex).__name__, str(ex)[:80], flush=True)
                pr.update(1)
    del tok
    # one image per calendar day (S30 and L30 may both exist)
    by_day = {}
    for g, z in kept:
        day = g.split(".")[3][:7]
        if day not in by_day or float(z["clear"]) > float(by_day[day][1]["clear"]):
            by_day[day] = (g, z)
    days = sorted(by_day)
    from datetime import datetime, timedelta
    label = lambda doy: (datetime(int(doy[:4]), 1, 1) + timedelta(days=int(doy[4:7]) - 1)).date()
    print("clear dates:", [str(label(d)) for d in days], flush=True)
    n = len(days)
    cols = 6
    rows = int(np.ceil(n / cols))
    fig, axs = plt.subplots(rows, cols, figsize=(3.1 * cols, 3.0 * rows + 0.6))
    axs = np.atleast_2d(axs)
    to_utm = Transformer.from_crs(4326, 32627, always_xy=True)
    pl = {k: ((to_utm.transform(*PLACES[k])[0] - sx) / 1000, (to_utm.transform(*PLACES[k])[1] - sy) / 1000) for k in ("Svartsengi", "Grindavik")}
    for a in axs.ravel():
        a.axis("off")
    for a, d in zip(axs.ravel(), days):
        g, z = by_day[d]
        im = np.dstack([z["swir2"], z["swir1"], z["red"]]).astype(float) * 1e-4
        im[im < -0.5] = np.nan
        im = np.clip((im - 0.0) / 0.35, 0, 1)                            # fixed stretch: dates comparable
        a.imshow(np.nan_to_num(im), extent=BOX, origin="upper")
        for k, (x0, y0) in pl.items():
            a.plot(x0, y0, "w^" if k == "Svartsengi" else "wo", mec="k", ms=4)
        dt = label(d)
        active = any(s <= str(dt) <= e for s, e in ERUPTIONS)
        a.set_title(f"{dt} {'(виверження)' if active else ''} · {g.split('.')[1]}", fontsize=8, color="firebrick" if active else "k")
        a.axis("on"); a.set_xticks([]); a.set_yticks([])
    fig.suptitle("Кратерний ряд Сундгнукур, листопад 2023 – серпень 2025: SWIR (B12, B11, B4), NASA HLS 30 м; лише дати, "
                 "де ≥ 85 % вікна без хмар і снігу; ▲ Свартсенгі, ● Ґріндавік; червоні заголовки — під час виверження", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(ML / "hls_timelapse.png", dpi=90, facecolor="white")
    (ML / "hls_timelapse.json").write_text(json.dumps(dict(n_checked=len(grans), clear_dates=[str(label(d)) for d in days],
                                                           granules=[by_day[d][0] for d in days]), indent=1), encoding="utf-8")
    print("wrote", ML / "hls_timelapse.png")


if __name__ == "__main__":
    main()
