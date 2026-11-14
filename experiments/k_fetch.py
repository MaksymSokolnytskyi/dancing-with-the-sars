"""K: fetch everything the Kakhovka radio-physics tests need (criteria in notes/experiments.md, section K).

- NISAR GCOV band B (80 m) windows for every scene of the four tracks over the upper bed, with incidence angle,
  noise-equivalent backscatter and the calibration metadata (src/hachure/gcov.py);
- band A (20 m) windows for one scene per track, the one closest to the Sentinel-2 summer date (test K1c);
- class layers on one 80 m UTM grid: ESA WorldCover 2021 (bed = reservoir water before 2023, cropland) and
  Sentinel-2 L2A (water, NDVI) on 2026-08-23 and 2026-03-11, plus the same S2 layers on each band A grid;
- Open-Meteo hourly 2 m air temperature at every overpass (frozen / unfrozen dates).

uv run python experiments/k_fetch.py      -> data/raw/gcov_cache/..., data/processed/kakhovka/{grid.npz, temps.json}
Resumable: cached GCOV windows are not read again.
"""
import json
import os
import shutil
import sys
import tempfile
import time
from datetime import datetime, timedelta
from pathlib import Path

import certifi
import numpy as np


def _ascii_ca_bundle():
    """GDAL's curl cannot open a CA file under a non-ASCII path (e.g. a Cyrillic Desktop): use an ASCII copy.
    Must run before rasterio is imported, which points GDAL_CURL_CA_BUNDLE at certifi."""
    src = Path(certifi.where())
    if str(src).isascii():
        return
    dst = Path(tempfile.gettempdir()) / "space-apps-cacert.pem"
    if not dst.exists() or dst.stat().st_size != src.stat().st_size:
        shutil.copyfile(src, dst)
    os.environ["GDAL_CURL_CA_BUNDLE"] = os.environ["PROJ_CURL_CA_BUNDLE"] = str(dst)


_ascii_ca_bundle()
import rasterio  # noqa: E402
import requests  # noqa: E402
from rasterio.transform import from_origin  # noqa: E402
from rasterio.warp import Resampling, reproject, transform_bounds  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src")]
from hachure import gcov  # noqa: E402

BBOX = (34.75, 47.35, 35.20, 47.70)          # upper bed (Velykyi Luh); natural surfaces only (see section K)
HEIGHT = 45.0                                 # bed height above the ellipsoid, m (~16 m a.s.l. + geoid ~29 m)
TRACKS = {42, 143, 36, 108}
S2_DATES = ("2026-08-23", "2026-03-11")      # 02-27 was pre-registered; see the K amendment
EPSG, RES = 32636, 80.0
OUT = ROOT / "data/processed/kakhovka"
UA = {"User-Agent": "space-apps-research"}
WC_URL = "https://esa-worldcover.s3.eu-central-1.amazonaws.com/v200/2021/map/ESA_WorldCover_10m_2021_v200_N45E033_Map.tif"
STAC = "https://earth-search.aws.element84.com/v1/search"


def target_grid():
    w, s, e, n = transform_bounds("EPSG:4326", f"EPSG:{EPSG}", *BBOX)
    w, n = np.floor(w / RES) * RES, np.ceil(n / RES) * RES
    W, H = int(np.ceil((e - w) / RES)), int(np.ceil((n - s) / RES))
    return from_origin(w, n, RES, RES), (H, W)


def warp(src, src_tf, src_crs, dst_tf, shape, how=Resampling.average):
    out = np.full(shape, np.nan, "float32")
    reproject(src.astype("float32"), out, src_transform=src_tf, src_crs=src_crs, dst_transform=dst_tf,
              dst_crs=f"EPSG:{EPSG}", resampling=how, src_nodata=np.nan, dst_nodata=np.nan)
    return out


def read_cog(url, dst_tf, shape, fn):
    """Window of a remote COG covering the grid, mapped by fn (per-pixel), then averaged onto the grid."""
    with rasterio.Env(GDAL_HTTP_USERAGENT="space-apps-research", GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR"):
        with rasterio.open("/vsicurl/" + url) as src:
            gb = rasterio.transform.array_bounds(*shape, dst_tf)
            b = transform_bounds(f"EPSG:{EPSG}", src.crs, gb[0], gb[1], gb[2], gb[3])
            win = rasterio.windows.from_bounds(*b, src.transform).round_offsets().round_lengths()
            a = src.read(1, window=win, boundless=True, fill_value=0)
            return warp(fn(a), src.window_transform(win), src.crs, dst_tf, shape)


def s2_item(date):
    r = requests.post(STAC, json={"collections": ["sentinel-2-l2a"], "bbox": list(BBOX),
                                  "datetime": f"{date}T00:00:00Z/{date}T23:59:59Z", "limit": 10}, headers=UA, timeout=60)
    items = sorted(r.json()["features"], key=lambda f: (f["properties"].get("s2:nodata_pixel_percentage", 100),
                                                          f["properties"]["eo:cloud_cover"]))
    return items[0]


def s2_layers(date, dst_tf, shape):
    it = s2_item(date)
    a = it["assets"]
    print(f"  S2 {date}: {it['id']} cloud {it['properties']['eo:cloud_cover']:.1f}%", flush=True)
    water = read_cog(a["scl"]["href"], dst_tf, shape, lambda v: np.where(v == 0, np.nan, (v == 6).astype("float32")))
    red = read_cog(a["red"]["href"], dst_tf, shape, lambda v: np.where(v == 0, np.nan, v.astype("float32")))
    nir = read_cog(a["nir"]["href"], dst_tf, shape, lambda v: np.where(v == 0, np.nan, v.astype("float32")))
    green = read_cog(a["green"]["href"], dst_tf, shape, lambda v: np.where(v == 0, np.nan, v.astype("float32")))
    return {"water": water, "ndvi": (nir - red) / (nir + red), "ndwi": (green - nir) / (green + nir)}


def grid_layers(dst_tf, shape):
    out = {}
    wc = lambda c: (lambda v: np.where(v == 0, np.nan, (v == c).astype("float32")))  # noqa: E731
    out["wc_water"] = read_cog(WC_URL, dst_tf, shape, wc(80))
    out["wc_crop"] = read_cog(WC_URL, dst_tf, shape, wc(40))
    for d in S2_DATES:
        for k, v in s2_layers(d, dst_tf, shape).items():
            out[f"s2_{d}_{k}"] = v
    return out


def temps(scenes):
    days = sorted({s.date.date() for s in scenes})
    lon, lat = (BBOX[0] + BBOX[2]) / 2, (BBOX[1] + BBOX[3]) / 2
    r = requests.get("https://archive-api.open-meteo.com/v1/archive",
                     params={"latitude": lat, "longitude": lon, "start_date": str(days[0] - timedelta(days=1)),
                             "end_date": str(days[-1]), "hourly": "temperature_2m", "timezone": "UTC"},
                     headers=UA, timeout=120).json()["hourly"]
    t = {datetime.fromisoformat(k): v for k, v in zip(r["time"], r["temperature_2m"])}
    out = {}
    for s in scenes:
        h = s.date.replace(minute=0, second=0)
        last = [t.get(h - timedelta(hours=i)) for i in range(24)]
        last = [v for v in last if v is not None]
        out[s.name] = {"t_overpass": t.get(h), "t_24h_mean": float(np.mean(last)) if last else None}
    return out


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    lon, lat = (BBOX[0] + BBOX[2]) / 2, (BBOX[1] + BBOX[3]) / 2
    scenes = [s for s in gcov.search(lon, lat) if s.path in TRACKS and s.pols == "DHDH"]
    print(f"{len(scenes)} GCOV scenes:", {p: sum(s.path == p for s in scenes) for p in sorted(TRACKS)}, flush=True)

    t0 = time.time()
    paths = gcov.fetch(scenes, BBOX, freq="B", height_m=HEIGHT)
    print(f"band B: {len(paths)} windows in {time.time() - t0:.0f}s", flush=True)
    s2day = datetime.fromisoformat(S2_DATES[0])
    a_scenes = [min((s for s in scenes if s.path == p), key=lambda s: abs(s.date - s2day)) for p in sorted(TRACKS)]
    a_paths = gcov.fetch(a_scenes, BBOX, freq="A", height_m=HEIGHT)
    print(f"band A: {len(a_paths)} windows", [s.date.date().isoformat() for s in a_scenes], flush=True)

    (OUT / "temps.json").write_text(json.dumps(temps(scenes), indent=1, default=str), encoding="utf-8")
    print("temperatures saved", flush=True)

    tf, shape = target_grid()
    lay = grid_layers(tf, shape)
    np.savez_compressed(OUT / "grid.npz", transform=np.array(tf)[:6], epsg=EPSG, **lay)
    for p in a_paths:                                   # S2 layers on each band A grid (K1c)
        w = gcov.Window.load(p)
        dx, dy = w.x[1] - w.x[0], w.y[1] - w.y[0]
        atf = rasterio.Affine(dx, 0, w.x[0] - dx / 2, 0, dy, w.y[0] - dy / 2)
        assert w.epsg == EPSG, w.epsg
        a = {}
        for d in S2_DATES:
            for k, v in s2_layers(d, atf, w.hh.shape).items():
                a[f"s2_{d}_{k}"] = v
        np.savez_compressed(OUT / f"gridA_{p.stem}.npz", **a)
    print(f"class layers saved in {OUT}", flush=True)


if __name__ == "__main__":
    main()
