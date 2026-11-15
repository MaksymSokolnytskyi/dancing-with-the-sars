"""WS: fetch what the Wisconsin story tests need (criteria in notes/experiments.md, section WS).

uv run python experiments/ws_fetch.py gcov [SUB|SUB2] [track]   band A full covariance, every QPDH scene
uv run python experiments/ws_fetch.py gunw      GUNW HH windows (phase, iono screen, coherence, components), every pair
uv run python experiments/ws_fetch.py s2        Sentinel-2 on the band B grid: NDVI 2026 (WS4), SCL of the lakes (WS3)
uv run python experiments/ws_fetch.py gim       CODE final IONEX maps for every band A scene day (CDDIS, Earthdata token)
uv run python experiments/ws_fetch.py weather   Open-Meteo hourly temperature and precipitation
Every step is resumable (cached files are not fetched again).
"""
import gzip
import json
import sys
import time
from pathlib import Path

import numpy as np
import requests

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src"), str(ROOT / "experiments")]
import k_fetch  # noqa: E402  (ASCII CA bundle for GDAL; must come before rasterio)
import rasterio  # noqa: E402
from rasterio.warp import Resampling, reproject, transform_bounds  # noqa: E402

from hachure import gcov, nisar  # noqa: E402
from w_fetch import BBOX, EPSG, HEIGHT, OUT, SUB, TRACKS  # noqa: E402

SUB2 = (-90.863, 45.140, -90.759, 45.214)
RAW = ROOT / "data/raw/ws"
CENTER = ((BBOX[0] + BBOX[2]) / 2, (BBOX[1] + BBOX[3]) / 2)


def scenes():
    return [s for s in gcov.search(*CENTER) if s.path in TRACKS and s.pols == "QPDH"]


def do_gcov():
    """Optional args: SUB | SUB2 and a track number, to run several readers in parallel."""
    sc = scenes()
    if len(sys.argv) > 3:
        sc = [s for s in sc if s.path == int(sys.argv[3])]
    for name, box in (("SUB", SUB), ("SUB2", SUB2)):
        if len(sys.argv) > 2 and sys.argv[2] != name:
            continue
        t0 = time.time()
        p = gcov.fetch(sc, box, freq="A", height_m=HEIGHT, offdiag=True, log=lambda m: print(m, flush=True))
        print(f"{name}: {len(p)} band A windows in {time.time() - t0:.0f}s", flush=True)


def do_gunw():
    d = RAW / "gunw"
    d.mkdir(parents=True, exist_ok=True)
    gs = [g for g in nisar.search(*CENTER) if g.path in TRACKS]
    print(len(gs), "GUNW pairs", flush=True)
    order = list(enumerate(gs))
    if len(sys.argv) > 2 and sys.argv[2] == "rev":          # a second reader from the other end
        order = order[::-1]
    for i, g in order:
        f = d / f"{g.name}.npz"
        if f.exists():
            continue
        t0 = time.time()
        with nisar.open_h5(g.url) as h5:
            raw = nisar.read_raw(h5, BBOX, g.days, pol="HH", height_m=HEIGHT)
        raw.save(f)
        print(f"  [{i + 1}/{len(gs)}] {g.path}{g.direction[0]} {g.ref_date:%Y-%m-%d} +{g.days:.0f}d "
              f"{f.stat().st_size / 1e6:.1f} MB {time.time() - t0:.0f}s", flush=True)


def read_onto(url, tf, shape, how, res=80.0):
    """A remote COG resampled onto the band B grid, read from its overviews."""
    dst = np.full(shape, np.nan, "float32")
    with rasterio.Env(GDAL_HTTP_USERAGENT="space-apps-research", GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR"):
        with rasterio.open("/vsicurl/" + url) as src:
            b = transform_bounds(f"EPSG:{EPSG}", src.crs, *rasterio.transform.array_bounds(*shape, tf))
            win = rasterio.windows.from_bounds(*b, src.transform).intersection(
                rasterio.windows.Window(0, 0, src.width, src.height)).round_offsets().round_lengths()
            f = max(1, int(res / abs(src.res[0]) / 2))                 # keep ~2 source samples per grid cell
            out = (max(1, int(win.height // f)), max(1, int(win.width // f)))
            a = src.read(1, window=win, out_shape=out, resampling=Resampling.nearest).astype("float32")
            wt = src.window_transform(win) * rasterio.Affine.scale(win.width / out[1], win.height / out[0])
            a[a == 0] = np.nan
            reproject(a, dst, src_transform=wt, src_crs=src.crs, dst_transform=tf, dst_crs=f"EPSG:{EPSG}",
                      resampling=how, src_nodata=np.nan, dst_nodata=np.nan)
    return dst


def s2_items(start, end, cloud):
    r = requests.post(k_fetch.STAC, json={"collections": ["sentinel-2-l2a"], "bbox": list(BBOX),
                                          "datetime": f"{start}T00:00:00Z/{end}T23:59:59Z", "limit": 300,
                                          "query": {"eo:cloud_cover": {"lt": cloud}}}, headers=k_fetch.UA, timeout=60)
    by = {}
    for f in r.json()["features"]:
        by.setdefault(f["properties"]["datetime"][:10], []).append(f)
    return dict(sorted(by.items()))


def do_s2():
    G = dict(np.load(OUT / "grid.npz"))
    tf, shape = rasterio.Affine(*G["transform"]), G["wc_water"].shape
    d = RAW / "s2"
    d.mkdir(parents=True, exist_ok=True)
    for tag, (start, end, cloud, bands) in {"ndvi": ("2026-06-01", "2026-11-10", 60, ("scl", "red", "nir")),
                                             "lake": ("2025-11-01", "2026-01-31", 90, ("scl",))}.items():
        for day, its in s2_items(start, end, cloud).items():
            f = d / f"{tag}_{day}.npz"
            if f.exists():
                continue
            lay = {}
            for it in sorted(its, key=lambda x: x["properties"].get("s2:nodata_pixel_percentage", 100)):
                a = it["assets"]
                one = {b: read_onto(a[b]["href"], tf, shape, Resampling.nearest if b == "scl" else Resampling.average)
                       for b in bands}
                for b, v in one.items():
                    lay[b] = v if b not in lay else np.where(np.isfinite(lay[b]), lay[b], v)
            np.savez_compressed(f, **lay)
            print(f"  S2 {tag} {day}: {len(its)} item(s), filled {np.isfinite(lay['scl']).mean():.0%}", flush=True)


def do_gim():
    import os
    d = RAW / "gim"
    d.mkdir(parents=True, exist_ok=True)
    tok = (Path(os.path.expanduser("~/.edl_token")).read_text().strip() if not os.environ.get("EARTHDATA_TOKEN")
           else os.environ["EARTHDATA_TOKEN"])
    days = sorted({s.date.date() for s in scenes()})
    for day in days:
        doy = day.timetuple().tm_yday
        name = f"COD0OPSFIN_{day.year}{doy:03d}0000_01D_01H_GIM.INX.gz"
        f = d / name.replace(".gz", "")
        if f.exists():
            continue
        url = f"https://cddis.nasa.gov/archive/gnss/products/ionex/{day.year}/{doy:03d}/{name}"
        r = requests.get(url, headers=k_fetch.UA | {"Authorization": f"Bearer {tok}"}, timeout=120)
        if r.status_code != 200 or r.content[:2] != b"\x1f\x8b":
            print(f"  GIM {day}: not available ({r.status_code})", flush=True)
            continue
        f.write_bytes(gzip.decompress(r.content))
        print(f"  GIM {day}: {f.stat().st_size / 1e3:.0f} kB", flush=True)


def do_weather():
    r = requests.get("https://archive-api.open-meteo.com/v1/archive",
                     params={"latitude": CENTER[1], "longitude": CENTER[0], "start_date": "2025-10-01",
                             "end_date": "2026-11-12", "hourly": "temperature_2m,precipitation,snow_depth",
                             "timezone": "UTC"}, headers=k_fetch.UA, timeout=120)
    r.raise_for_status()
    RAW.mkdir(parents=True, exist_ok=True)
    (RAW / "weather.json").write_text(json.dumps(r.json()["hourly"]), encoding="utf-8")
    print("weather hours:", len(r.json()["hourly"]["time"]))


if __name__ == "__main__":
    {"gcov": do_gcov, "gunw": do_gunw, "s2": do_s2, "gim": do_gim, "weather": do_weather}[sys.argv[1]]()
