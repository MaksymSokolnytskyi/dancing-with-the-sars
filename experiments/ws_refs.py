"""WS (prep, radar-independent): reference layers for the Wisconsin story tests (notes/experiments.md, section WS).

- USDA NASS Cropland Data Layer 2025 (Wisconsin state file from the CropScape cache), 30 m -> corn / soybean fractions;
- USFWS National Wetlands Inventory polygons (map service query) -> forested (PFO), scrub-shrub (PSS), emergent (PEM)
  wetland fractions;
both on the band B 80 m grid of section W (data/processed/wisconsin/grid.npz), plus the choice of a second band A
sub-window (SUB2) from these layers only: the 8 x 8 km box with the most corn + soybean that keeps >= 5 % forested wetland.

uv run python experiments/ws_refs.py      -> data/processed/wisconsin/ws_refs.npz, ws_refs.json, data/raw/ws/*
"""
import json
import re
import sys
from pathlib import Path

import numpy as np
import requests

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src"), str(ROOT / "experiments")]
import k_fetch  # noqa: E402  (ASCII CA bundle for GDAL; must come before rasterio)
import rasterio  # noqa: E402
from rasterio import features  # noqa: E402
from rasterio.warp import Resampling, transform_bounds  # noqa: E402

from w_fetch import BBOX, EPSG, OUT  # noqa: E402

RAW = ROOT / "data/raw/ws"
CDL_ZIP = "https://nassgeodata.gmu.edu/nass_data_cache/byfips/CDL_{year}_55.zip"
NWI = "https://fwspublicservices.wim.usgs.gov/wetlandsmapservice/rest/services/Wetlands/MapServer/0/query"
CORN, SOY = 1, 5
BOX_KM = 8.0


def cdl(year=2025):
    """Wisconsin state CDL (static CropScape cache, ~41 MB); the GetCDLFile service returned errors on 15.11."""
    f = RAW / f"CDL_{year}_55.zip"
    if not f.exists():
        f.write_bytes(requests.get(CDL_ZIP.format(year=year), headers=k_fetch.UA, timeout=900).content)
    return f"/vsizip/{f.as_posix()}/CDL_{year}_55.tif"


def nwi():
    f = RAW / "nwi.geojson"
    if not f.exists():
        feats, off = [], 0
        while True:
            r = requests.get(NWI, params={"geometry": ",".join(map(str, BBOX)), "geometryType": "esriGeometryEnvelope",
                                          "inSR": 4326, "outSR": 4326, "spatialRel": "esriSpatialRelIntersects",
                                          "outFields": "ATTRIBUTE,WETLAND_TYPE", "returnGeometry": "true",
                                          "resultOffset": off, "resultRecordCount": 500, "f": "geojson"},
                             headers=k_fetch.UA, timeout=180).json()
            feats += r["features"]
            print(f"  NWI {len(feats)} polygons", flush=True)
            if not r["features"] or not r.get("properties", {}).get("exceededTransferLimit", len(r["features"]) == 500):
                break
            off += len(r["features"])
        f.write_text(json.dumps({"type": "FeatureCollection", "features": feats}), encoding="utf-8")
    return json.loads(f.read_text(encoding="utf-8"))


def main():
    RAW.mkdir(parents=True, exist_ok=True)
    G = dict(np.load(OUT / "grid.npz"))
    tf, shape = rasterio.Affine(*G["transform"]), G["wc_water"].shape
    out = {}
    with rasterio.open(cdl()) as src:
        b = transform_bounds(f"EPSG:{EPSG}", src.crs, *rasterio.transform.array_bounds(*shape, tf))
        win = rasterio.windows.from_bounds(*b, src.transform).round_offsets().round_lengths()
        a, wtf = src.read(1, window=win), src.window_transform(win)
        for name, code in (("corn", CORN), ("soy", SOY)):
            out[f"cdl_{name}"] = k_fetch.warp((a == code).astype("float64"), wtf, src.crs, tf, shape, epsg=EPSG)
    fc = nwi()
    fine = 10.0                                                       # rasterise at 10 m, then average to 80 m
    k = int(round(abs(tf.a) / fine))
    ftf = rasterio.Affine(tf.a / k, 0, tf.c, 0, tf.e / k, tf.f)
    from rasterio.warp import transform_geom
    for name, pre in (("pfo", "PFO"), ("pss", "PSS"), ("pem", "PEM")):
        geoms = [transform_geom("EPSG:4326", f"EPSG:{EPSG}", f["geometry"]) for f in fc["features"]
                 if (f["properties"].get("ATTRIBUTE") or "").startswith(pre)]
        m = features.rasterize(((g, 1) for g in geoms), out_shape=(shape[0] * k, shape[1] * k), transform=ftf,
                               fill=0, dtype="uint8") if geoms else np.zeros((shape[0] * k, shape[1] * k), "uint8")
        out[f"nwi_{name}"] = m.reshape(shape[0], k, shape[1], k).mean((1, 3))
    n = int(round(BOX_KM * 1000 / abs(tf.a)))
    best = None
    for i in range(0, shape[0] - n + 1, 4):
        for j in range(0, shape[1] - n + 1, 4):
            s = (slice(i, i + n), slice(j, j + n))
            crop = float(np.nanmean(out["cdl_corn"][s] + out["cdl_soy"][s]))
            pfo = float(np.nanmean(out["nwi_pfo"][s]))
            if pfo >= 0.05 and (best is None or crop > best[0]):
                best = (crop, pfo, i, j)
    crop, pfo, i, j = best
    x0, y0 = tf * (j, i)
    x1, y1 = tf * (j + n, i + n)
    sub2 = transform_bounds(f"EPSG:{EPSG}", "EPSG:4326", min(x0, x1), min(y0, y1), max(x0, x1), max(y0, y1))
    info = {"sub2": [round(v, 3) for v in sub2], "sub2_corn_soy_share": round(crop, 3), "sub2_pfo_share": round(pfo, 3),
            "window_shares": {k2: round(float(np.nanmean(v)), 4) for k2, v in out.items()},
            "nwi_polygons": len(fc["features"])}
    np.savez_compressed(OUT / "ws_refs.npz", **out)
    (OUT / "ws_refs.json").write_text(json.dumps(info, indent=1), encoding="utf-8")
    print(json.dumps(info, indent=1))


if __name__ == "__main__":
    main()
