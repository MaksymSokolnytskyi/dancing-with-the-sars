"""Danube Delta: multi-year context + before/after data, read as windows from the
Microsoft Planetary Computer open archives (COGs, anonymous access, only the
delta window is transferred).

  s1       Sentinel-1 RTC (VV, VH), relative orbit 36 descending, 40 m, Oct 2014 - now.
           Monthly sample 2015-2023, every acquisition 2024-2026 (flood Sep/Oct 2024,
           the 2026 low-water season). Stored as uint8: dB = (v-1)/4 - 40, v = 0 is no data.
  landsat  Landsat 4-9 C2 L2, 30 m, 1984 - now: per year one best scene in the high-water
           season (15 Apr - 30 Jun) and one in the low-water season (15 Aug - 31 Oct).
           Stored as MNDWI and NDVI, int8 = round(index * 100), -128 = cloud / no data.
  static   JRC Global Surface Water 1984-2020 (occurrence, change, seasonality, recurrence,
           transitions, extent), Copernicus DEM 30 m, ESA WorldCover 2021.

Usage: python scripts/reykjanes/download_danube_history.py {s1|landsat|static|all} [--limit N]
"""
import csv
import sys
import threading
import time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from datetime import date

import h5py
import numpy as np
import rasterio
import requests
from pyproj import Transformer
from rasterio.enums import Resampling
from rasterio.transform import Affine
from rasterio.windows import from_bounds
from shapely.geometry import box, shape
from shapely.ops import unary_union

from inspect_gunw import ROOT

STAC = "https://planetarycomputer.microsoft.com/api/stac/v1"
WINDOW = (28.60, 44.30, 30.20, 45.60)  # lon0, lat0, lon1, lat1
OUT = ROOT / "data" / "danube"
UTM = "EPSG:32635"
WIN_POLY = box(*WINDOW)
CORE_POLY = box(28.70, 44.40, 29.95, 45.50)  # the delta proper: scenes are chosen by how much of it they cover


# ---------------------------------------------------------------- STAC + signing
def stac_search(**body):
    body.setdefault("limit", 500)
    out, url = [], f"{STAC}/search"
    while url:
        r = requests.post(url, json=body, timeout=120)
        r.raise_for_status()
        j = r.json()
        out += j["features"]
        nxt = [lk for lk in j.get("links", []) if lk["rel"] == "next"]
        url = nxt[0]["href"] if nxt else None
        body = nxt[0].get("body", body) if nxt else body
    return out


_tok, _lock = {}, threading.Lock()


def signed(href, coll):
    with _lock:
        t = _tok.get(coll)
        if not t or time.time() - t[1] > 1200:
            r = requests.get(f"https://planetarycomputer.microsoft.com/api/sas/v1/token/{coll}", timeout=60)
            r.raise_for_status()
            _tok[coll] = t = (r.json()["token"], time.time())
    return f"{href}?{t[0]}"


# ---------------------------------------------------------------- target grids
class Target:
    def __init__(self, crs, res, bounds, origin_offset=0.0):
        self.crs, self.res = crs, res
        minx, miny, maxx, maxy = bounds
        self.x0 = np.floor((minx - origin_offset) / res) * res + origin_offset
        self.y0 = np.ceil((maxy - origin_offset) / res) * res + origin_offset
        self.w = int(np.ceil((maxx - self.x0) / res))
        self.h = int(np.ceil((self.y0 - miny) / res))
        self.transform = Affine(res, 0, self.x0, 0, -res, self.y0)

    def attrs(self):
        return dict(crs=self.crs, transform=list(self.transform)[:6], res=self.res)


def utm_target(res, offset=0.0):
    b = Transformer.from_crs(4326, UTM, always_xy=True).transform_bounds(*WINDOW, densify_pts=21)
    return Target(UTM, res, b, offset)


def read_item(href, coll, tgt, resampling, fill):
    """One band of one item on the target grid (None if the item does not overlap or has another CRS)."""
    with rasterio.open(signed(href, coll)) as ds:
        if str(ds.crs).upper() != tgt.crs.upper():
            return None
        l, b, r, t = ds.bounds
        c0 = max(int(np.floor((l - tgt.x0) / tgt.res)), 0)
        c1 = min(int(np.ceil((r - tgt.x0) / tgt.res)), tgt.w)
        r0 = max(int(np.floor((tgt.y0 - t) / tgt.res)), 0)
        r1 = min(int(np.ceil((tgt.y0 - b) / tgt.res)), tgt.h)
        if c1 <= c0 or r1 <= r0:
            return None
        win = from_bounds(tgt.x0 + c0 * tgt.res, tgt.y0 - r1 * tgt.res, tgt.x0 + c1 * tgt.res,
                          tgt.y0 - r0 * tgt.res, transform=ds.transform)
        a = ds.read(1, window=win, out_shape=(r1 - r0, c1 - c0), resampling=resampling,
                    boundless=True, fill_value=fill)
        return a, (r0, r1, c0, c1)


def paste(dest, got, valid_fn=lambda a: np.ones(a.shape, bool)):
    if got is None:
        return
    a, (r0, r1, c0, c1) = got
    d = dest[r0:r1, c0:c1]
    m = np.isnan(d) & valid_fn(a) if d.dtype.kind == "f" else (d == 0) & valid_fn(a)
    d[m] = a[m]


def write_common(out, tgt, **attrs):
    out.attrs.update(tgt.attrs())
    out.attrs.update(attrs)


def save(path, arrays, tgt, **attrs):
    tmp = path.with_suffix(".part")
    with h5py.File(tmp, "w") as out:
        for k, v in arrays.items():
            out.create_dataset(k, data=v, compression="gzip", compression_opts=4, shuffle=True,
                               chunks=(min(256, v.shape[0]), min(256, v.shape[1])))
        write_common(out, tgt, **attrs)
    tmp.rename(path)


# ---------------------------------------------------------------- Sentinel-1
S1_ORBIT, S1_RES = 36, 40.0


def s1_dates():
    items = stac_search(collections=["sentinel-1-rtc"], bbox=list(WINDOW), datetime="2014-10-01/2026-12-31")
    by = defaultdict(list)
    for it in items:
        p = it["properties"]
        if p["sat:relative_orbit"] == S1_ORBIT and p["proj:epsg"] == 32635:
            by[p["datetime"][:10]].append(it)
    days = sorted(by)
    pick, seen = [], set()
    for d in days:
        if d >= "2024-01-01" or d[:7] not in seen:  # monthly before 2024, everything after
            if d < "2024-01-01":
                seen.add(d[:7])
            pick.append(d)
    return pick, by


def s1_one(d, its, tgt):
    dst = OUT / "s1" / f"s1_{d.replace('-', '')}_o{S1_ORBIT}D.h5"
    if dst.exists():
        return [d, "skip", "", ""]
    pols = {}
    for pol in ("vv", "vh"):
        arr = np.full((tgt.h, tgt.w), np.nan, np.float32)
        for it in its:
            paste(arr, read_item(it["assets"][pol]["href"], "sentinel-1-rtc", tgt, Resampling.average, -32768.0),
                  lambda a: a > 0)
        db = 10 * np.log10(np.where(arr > 0, arr, np.nan))
        v = np.clip(np.round((db + 40) * 4) + 1, 1, 255)
        pols[pol] = np.where(np.isfinite(v), v, 0).astype(np.uint8)
    share = float((pols["vv"] > 0).mean())
    save(dst, pols, tgt, date=d, items=",".join(i["id"] for i in its), valid_share=share,
         platform=its[0]["properties"]["platform"], encoding="dB = (v-1)/4 - 40; 0 = no data")
    return [d, "done", f"{share:.2f}", f"{dst.stat().st_size / 1e6:.0f} MB"]


def run_s1(limit):
    (OUT / "s1").mkdir(parents=True, exist_ok=True)
    pick, by = s1_dates()
    pick = pick[:limit] if limit else pick
    print("s1:", len(pick), "dates", pick[0], "->", pick[-1], flush=True)
    tgt = utm_target(S1_RES)
    rows = []
    with ThreadPoolExecutor(6) as ex:
        for row in ex.map(lambda d: s1_one(d, by[d], tgt), pick):
            rows.append(row)
            print("s1", *row, flush=True)
    write_manifest("s1", rows)


# ---------------------------------------------------------------- Landsat
SEASONS = {"high": ("04-15", "06-30"), "low": ("08-15", "10-31")}


def landsat_picks():
    items = stac_search(collections=["landsat-c2-l2"], bbox=list(WINDOW), datetime="1984-01-01/2026-12-31",
                        query={"eo:cloud_cover": {"lt": 60}})
    by = defaultdict(list)
    for it in items:
        if it["properties"]["proj:epsg"] == 32635:
            by[it["properties"]["datetime"][:10]].append(it)
    scored = []
    for d, its in by.items():
        cover = unary_union([shape(i["geometry"]) for i in its]).intersection(CORE_POLY).area / CORE_POLY.area
        cloud = float(np.mean([i["properties"]["eo:cloud_cover"] for i in its]))
        scored.append((d, cover, cloud, its))
    picks = []
    for year in range(1984, 2027):
        for season, (a, b) in SEASONS.items():
            c = [s for s in scored if f"{year}-{a}" <= s[0] <= f"{year}-{b}" and s[1] >= 0.85]
            if c:
                best = min(c, key=lambda s: s[2])
                picks.append((year, season, *best))
    return picks


def landsat_one(pick, tgt):
    year, season, d, cover, cloud, its = pick
    dst = OUT / "landsat" / f"landsat_{d.replace('-', '')}_{season}.h5"
    if dst.exists():
        return [d, season, "skip", "", ""]
    res = {}
    for band in ("red", "green", "nir08", "swir16", "qa_pixel"):
        arr = np.zeros((tgt.h, tgt.w), np.uint16)
        for it in its:
            paste(arr, read_item(it["assets"][band]["href"], "landsat-c2-l2", tgt, Resampling.nearest, 0))
        res[band] = arr
    qa = res["qa_pixel"]
    bad = (qa & 0b11110) != 0  # dilated cloud, cirrus, cloud, cloud shadow
    refl = {k: res[k].astype(np.float32) * 2.75e-5 - 0.2 for k in ("red", "green", "nir08", "swir16")}
    ok = ~bad & (res["green"] > 0) & (res["swir16"] > 0) & (res["nir08"] > 0) & (res["red"] > 0)
    with np.errstate(invalid="ignore", divide="ignore"):
        mndwi = (refl["green"] - refl["swir16"]) / (refl["green"] + refl["swir16"])
        ndvi = (refl["nir08"] - refl["red"]) / (refl["nir08"] + refl["red"])
    q = lambda x: np.where(ok & np.isfinite(x), np.clip(np.round(x * 100), -100, 100), -128).astype(np.int8)
    save(dst, {"mndwi": q(mndwi), "ndvi": q(ndvi)}, tgt, date=d, season=season, platform=its[0]["properties"]["platform"],
         items=",".join(i["id"] for i in its), scene_cloud_cover=cloud, core_cover=cover,
         valid_share=float(ok.mean()), encoding="int8 = round(index*100); -128 = cloud / no data")
    return [d, season, "done", f"{ok.mean():.2f}", f"{dst.stat().st_size / 1e6:.0f} MB"]


def run_landsat(limit):
    (OUT / "landsat").mkdir(parents=True, exist_ok=True)
    picks = landsat_picks()
    picks = picks[-limit:] if limit else picks
    print("landsat:", len(picks), "scenes", picks[0][2], "->", picks[-1][2], flush=True)
    tgt = utm_target(30.0, 15.0)
    rows = []
    with ThreadPoolExecutor(5) as ex:
        for row in ex.map(lambda p: landsat_one(p, tgt), picks):
            rows.append(row)
            print("landsat", *row, flush=True)
    write_manifest("landsat", rows)


# ---------------------------------------------------------------- static layers
def mosaic_latlon(coll, asset, res, resampling, dtype, fill, items):
    tgt = Target("EPSG:4326", res, WINDOW)
    tgt.transform = Affine(res, 0, tgt.x0, 0, -res, tgt.y0)
    dest = np.full((tgt.h, tgt.w), np.nan if np.dtype(dtype).kind == "f" else 0, dtype)
    for it in items:
        paste(dest, read_item(it["assets"][asset]["href"], coll, tgt, resampling, fill),
              (lambda a: a != fill) if fill is not None else (lambda a: np.ones(a.shape, bool)))
    return dest, tgt


def run_static():
    (OUT / "static").mkdir(parents=True, exist_ok=True)
    jobs = [("jrc-gsw", a, 0.00025, Resampling.nearest, np.uint8, 255)
            for a in ("occurrence", "change", "seasonality", "recurrence", "transitions", "extent")]
    for coll, asset, res, rs, dt, fill in jobs:
        items = stac_search(collections=[coll], bbox=list(WINDOW))
        dst = OUT / "static" / f"gsw_{asset}.h5"
        if dst.exists():
            continue
        arr, tgt = mosaic_latlon(coll, asset, res, rs, dt, fill, items)
        save(dst, {asset: arr}, tgt, source="JRC Global Surface Water v1.4, 1984-2020", asset=asset)
        print("static", dst.name, arr.shape, f"{(arr > 0).mean():.2f} nonzero", flush=True)
    items = stac_search(collections=["cop-dem-glo-30"], bbox=list(WINDOW))
    dst = OUT / "static" / "cop_dem_glo30.h5"
    if not dst.exists():
        arr, tgt = mosaic_latlon("cop-dem-glo-30", "data", 1 / 3600, Resampling.bilinear, np.float32, None, items)
        save(dst, {"elevation_m": arr}, tgt, source="Copernicus DEM GLO-30 (EGM2008 heights)")
        print("static", dst.name, arr.shape, np.nanmin(arr), np.nanmax(arr), flush=True)
    items = stac_search(collections=["esa-worldcover"], bbox=list(WINDOW), datetime="2021-01-01/2021-12-31")
    dst = OUT / "static" / "worldcover_2021.h5"
    if not dst.exists():
        arr, tgt = mosaic_latlon("esa-worldcover", "map", 0.000333333, Resampling.nearest, np.uint8, 0, items)
        save(dst, {"class": arr}, tgt, source="ESA WorldCover 2021 v200 (10 m, subsampled to ~37 m)")
        print("static", dst.name, arr.shape, np.unique(arr), flush=True)


def write_manifest(name, rows):
    with open(OUT / f"manifest_{name}.csv", "w", newline="") as f:
        csv.writer(f).writerows(rows)


if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "all"
    limit = int(sys.argv[sys.argv.index("--limit") + 1]) if "--limit" in sys.argv else 0
    OUT.mkdir(parents=True, exist_ok=True)
    if what in ("static", "all"):
        run_static()
    if what in ("landsat", "all"):
        run_landsat(limit)
    if what in ("s1", "all"):
        run_s1(limit)
