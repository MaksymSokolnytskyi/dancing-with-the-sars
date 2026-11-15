"""Why does EGMS have no data inside most Ostrava bowls? A post-hoc diagnostic (notes/experiments.md, block O).

    uv run python scripts/ostrava/egms_nodata.py

Puts on the EGMS 100 m grid (EPSG:3035), over the study box: whether the EGMS L3 vertical cell has data, the dominant
ESA WorldCover class, the share of built-up land, the NISAR 2026 vertical rate of the study, and the distance to the nearest
bowl. Then prints the share of EGMS cells with data:
  1. by land cover, far from the bowls (> 3 km): what EGMS sees on ordinary ground of each kind;
  2. inside the bowls, by land cover, against the same classes far away;
  3. by NISAR 2026 vertical rate, for cropland/grassland cells only (land cover held fixed);
  4. by distance to the nearest bowl, cropland/grassland only.
Also the C-band (Sentinel-1) limit: the largest LOS change between two acquisitions that phase can follow, lambda/4.
Inputs: data/ostrava/egms/*.tiff (egms_fetch.py), data/raw/context/ostrava-karvina-asc_worldcover_10m.tif (sars fetch),
output/ostrava_karvina/{vertical_cm_yr.tif, zones_vertical_osm.geojson} (the study).
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import rasterio
from rasterio.enums import Resampling
from rasterio.features import rasterize
from rasterio.transform import from_origin
from rasterio.warp import reproject, transform_bounds, transform_geom
from scipy.ndimage import distance_transform_edt

ROOT = Path(__file__).resolve().parents[2]
BBOX = (18.13, 49.68, 18.83, 50.07)
WC = ROOT / "data/raw/context/ostrava-karvina-asc_worldcover_10m.tif"
NISAR = ROOT / "output/ostrava_karvina/vertical_cm_yr.tif"
ZONES = ROOT / "output/ostrava_karvina/zones_vertical_osm.geojson"
CLASSES = {10: "tree cover", 20: "shrubland", 30: "grassland", 40: "cropland", 50: "built-up", 60: "bare", 80: "water", 90: "wetland"}


def main():
    l, b, r, t = transform_bounds("EPSG:4326", "EPSG:3035", *BBOX)
    l, b, r, t = np.floor(l / 100) * 100, np.floor(b / 100) * 100, np.ceil(r / 100) * 100, np.ceil(t / 100) * 100
    W, H = int((r - l) / 100), int((t - b) / 100)
    T = from_origin(l, t, 100, 100)
    egms = np.full((H, W), np.nan)
    for p in sorted((ROOT / "data/ostrava/egms").glob("EGMS_L3_*_U_*.tiff")):
        with rasterio.open(p) as s:
            a = np.full((H, W), np.nan)
            reproject(s.read(1), a, src_transform=s.transform, src_crs="EPSG:3035", src_nodata=s.nodata,
                      dst_transform=T, dst_crs="EPSG:3035", dst_nodata=np.nan, resampling=Resampling.nearest)
            egms = np.where(np.isfinite(a), a, egms)
    with rasterio.open(WC) as s:
        wc = s.read(1)
        dom = np.zeros((H, W), "uint8")
        reproject(wc, dom, src_transform=s.transform, src_crs=s.crs, dst_transform=T, dst_crs="EPSG:3035",
                  resampling=Resampling.mode)
        built = np.zeros((H, W), "float32")
        reproject((wc == 50).astype("float32"), built, src_transform=s.transform, src_crs=s.crs, dst_transform=T,
                  dst_crs="EPSG:3035", resampling=Resampling.average)
    with rasterio.open(NISAR) as s:
        up = np.full((H, W), np.nan, "float32")
        reproject(s.read(1), up, src_transform=s.transform, src_crs=s.crs, src_nodata=np.nan, dst_transform=T,
                  dst_crs="EPSG:3035", dst_nodata=np.nan, resampling=Resampling.average)
    feats = json.loads(ZONES.read_text(encoding="utf-8"))["features"]
    bowls = rasterize([(transform_geom("EPSG:4326", "EPSG:3035", f["geometry"]), 1) for f in feats], out_shape=(H, W),
                      transform=T, fill=0).astype(bool)
    dist = distance_transform_edt(~bowls) * 100            # m to the nearest bowl cell (0 inside)
    has = np.isfinite(egms)
    land = dom != 0
    print(f"grid {W} x {H} cells of 100 m; EGMS data in {has[land].mean():.1%} of land cells; bowls {bowls.sum()} cells\n")

    print("1. Far from the bowls (> 3 km): share of EGMS cells with data, by WorldCover class")
    far = land & (dist > 3000)
    for c, nm in CLASSES.items():
        m = far & (dom == c)
        if m.sum() > 50:
            print(f"   {nm:11s} {has[m].mean():6.1%}  ({m.sum()} cells)")

    print("\n2. Inside the bowls vs the same class far away")
    for c, nm in CLASSES.items():
        m = bowls & (dom == c)
        if m.sum():
            print(f"   {nm:11s} inside {has[m].mean():6.1%} ({m.sum()} cells)   far {has[far & (dom == c)].mean():6.1%}")
    m = bowls & (built >= 0.3)
    mf = far & (built >= 0.3)
    print(f"   built-up >= 30 % of the cell: inside {has[m].mean() if m.sum() else float('nan'):6.1%} ({m.sum()} cells)   "
          f"far {has[mf].mean():6.1%}")

    print("\n3. Cropland / grassland only: share with EGMS data by NISAR 2026 vertical rate (study, two tracks)")
    ag = land & np.isin(dom, [30, 40]) & np.isfinite(up)
    for lo, hi in [(-5, 99), (-10, -5), (-15, -10), (-25, -15), (-40, -25), (-999, -40)]:
        m = ag & (up > lo) & (up <= hi)
        if m.sum():
            print(f"   {lo:5d} < up <= {hi:3d} cm/yr: {has[m].mean():6.1%}  ({m.sum()} cells)")

    print("\n4. Cropland / grassland only: share with EGMS data by distance to the nearest bowl")
    ag2 = land & np.isin(dom, [30, 40])
    for lo, hi in [(-1, 0), (0, 300), (300, 600), (600, 1000), (1000, 2000), (2000, 3000), (3000, 1e9)]:
        m = ag2 & (dist > lo) & (dist <= hi)
        lab = "inside" if hi == 0 else f"{lo:.0f}-{hi:.0f} m" if hi < 1e8 else f"> {lo:.0f} m"
        print(f"   {lab:12s} {has[m].mean():6.1%}  ({m.sum()} cells)")

    time_series(l, t, W, H, T, feats)

    lam = 5.547                                              # cm, Sentinel-1 C-band
    print(f"\nC-band: lambda/4 = {lam / 4:.2f} cm of LOS change between two acquisitions; "
          f"as a rate: {lam / 4 / 6 * 365.25:.0f} cm/yr with a 6-day revisit (2019-2021), "
          f"{lam / 4 / 12 * 365.25:.0f} cm/yr with 12 days (from Dec 2021)")


def time_series(l, t, W, H, T, feats):
    """5. EGMS cells that survive near the PL bowls: rate by distance to the bowl edge and by year (needs pandas)."""
    try:
        import pandas as pd
    except ImportError:
        print("\n5. skipped: run with `uv run --with pandas==3.1.1 python scripts/ostrava/egms_nodata.py` for the time series")
        return
    import zipfile
    pl = [f for f in feats if f["properties"]["country"] == "PL"]
    m = rasterize([(transform_geom("EPSG:4326", "EPSG:3035", f["geometry"]), 1) for f in pl], out_shape=(H, W),
                  transform=T, fill=0).astype(bool)
    dist = distance_transform_edt(~m) * 100
    parts = []
    for z in sorted((ROOT / "data/ostrava/egms").glob("EGMS_L3_*_U_*.zip")):
        with zipfile.ZipFile(z) as zf:
            name = next(n for n in zf.namelist() if n.endswith(".csv"))
            for ch in pd.read_csv(zf.open(name), chunksize=50000):
                j = ((ch.easting.values - l) // 100).astype(int)
                i = ((t - ch.northing.values) // 100).astype(int)
                ok = (i >= 0) & (i < H) & (j >= 0) & (j < W)
                ch = ch[ok].copy()
                ch["dist"] = dist[i[ok], j[ok]]
                parts.append(ch[ch["dist"] <= 600])
    df = pd.concat(parts)
    dates = [c for c in df.columns if c.isdigit()]
    dd = pd.to_datetime(dates)
    yrs = (dd - dd[0]).days.values / 365.25
    X = df[dates].values
    for y in range(2019, 2024):
        k = dd.year == y
        df[f"rate_{y}"] = [np.polyfit(yrs[k][np.isfinite(r[k])], r[k][np.isfinite(r[k])], 1)[0]
                           if np.isfinite(r[k]).sum() > 10 else np.nan for r in X]
    df["band"] = pd.cut(df["dist"], [-1, 0, 200, 400, 600], labels=["inside", "0-200 m", "200-400 m", "400-600 m"])
    cols = ["mean_velocity", "rmse"] + [f"rate_{y}" for y in range(2019, 2024)]
    tab = df.groupby("band", observed=True)[cols].median() / 10
    tab["rmse"] *= 10
    tab["cells"] = df.groupby("band", observed=True).size()
    print("\n5. EGMS cells near the PL bowls that have data: median vertical rate (cm/yr; rmse in mm) by distance to the bowl")
    print(tab.round(2).to_string())


if __name__ == "__main__":
    main()
