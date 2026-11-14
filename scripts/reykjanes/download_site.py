"""Read NISAR L2 windows over a study site straight from the ASF archive with
HTTP range requests (only the site window is transferred, not the 2 GB granule).

  GUNW: unwrapped phase, coherence, connected components, ionosphere screen, water
        mask (80 m); optionally pixel offsets (80 m) and the 20 m wrapped
        interferogram; incidence angle, LOS unit vectors and troposphere screens
        from the radar-grid cube at the site height.
  GCOV: HH and HV gamma0 (RTC) at 20 m, stored as int16 = round(dB * 100),
        -32768 = no data.

Every pair / acquisition over the site is taken, from all tracks and frames.
When a granule exists in several processing versions, the newest is kept.

Usage: python scripts/reykjanes/download_site.py {reykjanes|mississippi|kilauea} [--list] [--products GUNW,GCOV]
       [--events YYYYMMDD,...]  only granules around these dates
"""
import csv
import os
import json
import re
import sys
import urllib.request
from concurrent.futures import ProcessPoolExecutor
from functools import partial

import h5py
import numpy as np
from pyproj import Transformer

from inspect_gunw import ROOT, open_remote

COLLECTIONS = {
    "GUNW": ["C3622247503-ASF", "C2854335566-ASF", "C2850261892-ASF"],  # V1, PROVISIONAL_V1, BETA_V1
    "GCOV": ["C3622241360-ASF", "C2854338529-ASF", "C3622214170-ASF"],
}
SITES = {
    # Svartsengi - Sundhnukur - Grindavik - Fagradalsfjall dike / uplift system
    "reykjanes": dict(point=(-22.43, 63.88), window=(-22.80, 63.76, -22.05, 64.02), height=60.0,
                      products=["GUNW"], offsets=True, wrapped=True),
    # Atchafalaya + Wax Lake deltas, Terrebonne / Barataria marshes
    "mississippi": dict(point=(-91.20, 29.50), window=(-91.60, 29.20, -90.50, 29.85), height=0.0,
                        products=["GCOV", "GUNW"], offsets=False, wrapped=False),
    # Kilauea summit (Halema'uma'u) + upper East and Southwest rift zones
    "kilauea": dict(point=(-155.287, 19.407), window=(-155.40, 19.30, -155.15, 19.48), height=1100.0,
                    products=["GUNW", "GCOV"], offsets=True, wrapped=False),
}
VERSION = re.compile(r"_([PXB])(\d{5})_")


def find(site, product):
    lon, lat = SITES[site]["point"]
    best = {}
    for coll in COLLECTIONS[product]:
        url = ("https://cmr.earthdata.nasa.gov/search/granules.umm_json"
               f"?collection_concept_id={coll}&point={lon},{lat}&page_size=2000")
        for it in json.load(urllib.request.urlopen(url))["items"]:
            name = it["umm"]["GranuleUR"]
            m = VERSION.search(name)
            key = re.sub(r"^NISAR_L2_\w+?_(GUNW|GCOV)_", r"\1_", name[:m.start()]) if m else name
            href = next(r["URL"] for r in it["umm"]["RelatedUrls"]
                        if r["URL"].startswith("https://") and r["URL"].endswith(".h5"))
            ver = int(m.group(2)) if m else 0
            if key not in best or ver > best[key][0]:
                best[key] = (ver, name, href)
    return sorted((v[1], v[2]) for v in best.values())


STAMP = re.compile(r"_(\d{8})T\d{6}")


def by_events(items, product, events):
    """Keep only granules around event dates (YYYYMMDD): GUNW pairs whose interval contains the
    event (reference date < event <= secondary date); GCOV: per track, the last acquisition
    before and the first on/after each event."""
    track = lambda n: re.search(r"_(\d{3}_[AD])_\d{3}_", n).group(1)
    keep = set()
    for ev in events:
        if product == "GUNW":
            for n, h in items:
                d = STAMP.findall(n)
                if d[0] < ev <= d[2]:
                    keep.add((n, h))
        else:
            for tr in {track(n) for n, _ in items}:
                acq = sorted((STAMP.findall(n)[0], n, h) for n, h in items if track(n) == tr)
                before = [a for a in acq if a[0] < ev]
                after = [a for a in acq if a[0] >= ev]
                for a in (before[-1:] + after[:1]):
                    keep.add(a[1:])
    return sorted(keep)


def window_slices(x, y, epsg, win):
    tr = Transformer.from_crs(4326, epsg, always_xy=True)
    cx, cy = tr.transform_bounds(*win, densify_pts=21)[::2], tr.transform_bounds(*win, densify_pts=21)[1::2]
    ix = np.where((x >= min(cx)) & (x <= max(cx)))[0]
    iy = np.where((y >= min(cy)) & (y <= max(cy)))[0]
    if ix.size == 0 or iy.size == 0:
        return None
    return slice(iy[0], iy[-1] + 1), slice(ix[0], ix[-1] + 1)


def gz(out, k, a):
    out.create_dataset(k, data=a, compression="gzip", compression_opts=4, shuffle=True,
                       chunks=(min(256, a.shape[0]), min(256, a.shape[1])) if a.ndim == 2 else None)


def radar_cube(src_root, out, height):
    r = src_root["metadata/radarGrid"]
    h = r["heightAboveEllipsoid"][()]
    i = int(np.argmin(np.abs(h - height)))
    for k in ("incidenceAngle", "losUnitVectorX", "losUnitVectorY", "wetTroposphericPhaseScreen",
              "hydrostaticTroposphericPhaseScreen", "slantRangeSolidEarthTidesPhase"):
        if k in r:
            gz(out, f"cube/{k}", r[k][i])
    out["cube/x"], out["cube/y"] = r["xCoordinates"][()], r["yCoordinates"][()]
    out["cube"].attrs["height"] = float(h[i])


def fetch_gunw(src, out, cfg):
    root = src["science/LSAR/GUNW"]
    g = root["grids/frequencyA/unwrappedInterferogram"]
    x, y, epsg = g["HH/xCoordinates"][()], g["HH/yCoordinates"][()], int(g["HH/projection"][()])
    s = window_slices(x, y, epsg, cfg["window"])
    if s is None:
        return None
    for k in ("unwrappedPhase", "coherenceMagnitude", "connectedComponents", "ionospherePhaseScreen"):
        gz(out, f"unw/{k}", g["HH"][k][s])
    gz(out, "unw/mask", g["mask"][s])
    out["unw/x"], out["unw/y"] = x[s[1]], y[s[0]]
    valid = float(np.isfinite(out["unw/unwrappedPhase"][()]).mean())
    if cfg["offsets"]:
        o = root["grids/frequencyA/pixelOffsets"]
        so = window_slices(o["HH/xCoordinates"][()], o["HH/yCoordinates"][()], epsg, cfg["window"])
        for k in ("alongTrackOffset", "slantRangeOffset", "correlationSurfacePeak"):
            gz(out, f"off/{k}", o["HH"][k][so])
        out["off/x"], out["off/y"] = o["HH/xCoordinates"][so[1]], o["HH/yCoordinates"][so[0]]
    if cfg["wrapped"]:
        w = root["grids/frequencyA/wrappedInterferogram"]
        sw = window_slices(w["HH/xCoordinates"][()], w["HH/yCoordinates"][()], epsg, cfg["window"])
        gz(out, "wrp/wrappedInterferogram", w["HH/wrappedInterferogram"][sw])
        gz(out, "wrp/coherenceMagnitude", w["HH/coherenceMagnitude"][sw])
        out["wrp/x"], out["wrp/y"] = w["HH/xCoordinates"][sw[1]], w["HH/yCoordinates"][sw[0]]
    radar_cube(root, out, cfg["height"])
    out.attrs["centerFrequency"] = float(root["grids/frequencyA/centerFrequency"][()])
    out.attrs["epsg"] = epsg
    return valid


def fetch_gcov(src, out, cfg):
    g = src["science/LSAR/GCOV/grids/frequencyA"]
    x, y, epsg = g["xCoordinates"][()], g["yCoordinates"][()], int(g["projection"][()])
    s = window_slices(x, y, epsg, cfg["window"])
    if s is None:
        return None
    terms = [t.decode() for t in g["listOfCovarianceTerms"][()]]
    shape = (s[0].stop - s[0].start, s[1].stop - s[1].start)
    step = 512  # rows per block: keeps a worker to ~tens of MB on a 6 GB machine
    for k in ("HHHH", "HVHV", "VVVV", "VHVH"):
        if k in terms:
            d = out.create_dataset(k + "_db100", shape, np.int16, compression="gzip", compression_opts=4,
                                   shuffle=True, chunks=(min(256, shape[0]), min(256, shape[1])))
            for r0 in range(0, shape[0], step):
                r1 = min(r0 + step, shape[0])
                a = g[k][s[0].start + r0:s[0].start + r1, s[1]]
                with np.errstate(divide="ignore", invalid="ignore"):
                    db = 10 * np.log10(np.where(a > 0, a, np.nan))
                d[r0:r1] = np.where(np.isfinite(db), np.round(db * 100), -32768).astype(np.int16)
                del a, db
    gz(out, "mask", g["mask"][s])
    out["x"], out["y"] = x[s[1]], y[s[0]]
    out.attrs.update(epsg=epsg, encoding="dB = value / 100; -32768 = no data; gamma0 (RTC)",
                     covariance_terms=",".join(terms))
    if "HHHH_db100" not in out:
        return 0.0
    h = out["HHHH_db100"]
    return float(sum(int((h[r:r + 512] != -32768).sum()) for r in range(0, shape[0], 512)) / (shape[0] * shape[1]))


def fetch(site, product, item, attempts=4):
    """One granule; a network failure is retried with a fresh signed URL and never stops the run."""
    import time
    for k in range(attempts):
        try:
            return fetch_once(site, product, item)
        except Exception as e:  # timeouts, dropped connections, expired signatures
            tmp = ROOT / "data" / site / product.lower() / f"{item[0]}.part"
            try:
                tmp.unlink(missing_ok=True)
            except OSError:  # still held open by the failed h5py handle; fetch_once rewrites it
                pass
            err = f"{type(e).__name__}: {e}"[:200]
            time.sleep(10 * (k + 1))
    return [item[0], "error", err]


def fetch_once(site, product, item):
    name, href = item
    cfg = dict(SITES[site])
    if os.environ.get("NISAR_LITE"):  # unwrapped phase + cube only: no 20 m wrapped, no offsets
        cfg.update(offsets=False, wrapped=False)
    dst = ROOT / "data" / site / product.lower() / f"{name}.h5"
    moved = dst.parent / "on_phone.txt"  # written by sync_phone.ps1 after verifying the copy
    if dst.exists() or (moved.exists() and name in moved.read_text(encoding="utf-8-sig").split()):
        return [name, "skip", ""]
    tmp = dst.with_suffix(".part")
    with open_remote(href) as src, h5py.File(tmp, "w") as out:
        out.attrs.update(granule=name, site=site, window_lonlat=cfg["window"])
        valid = (fetch_gunw if product == "GUNW" else fetch_gcov)(src, out, cfg)
    if valid is None:
        tmp.unlink()
        return [name, "no overlap", ""]
    tmp.rename(dst)
    return [name, "done", f"{valid:.2f}"]


if __name__ == "__main__":
    site = sys.argv[1]
    products = (sys.argv[sys.argv.index("--products") + 1].split(",") if "--products" in sys.argv
                else SITES[site]["products"])
    for product in products:
        items = find(site, product)
        if "--tracks" in sys.argv:  # e.g. 159_A,082_D
            want = sys.argv[sys.argv.index("--tracks") + 1].split(",")
            items = [it for it in items if any(f"_{w}_" in it[0] for w in want)]
        if "--events" in sys.argv:
            items = by_events(items, product, sys.argv[sys.argv.index("--events") + 1].split(","))
        tracks = sorted({re.search(r"_(\d{3}_[AD]_\d{3})_", n).group(1) for n, _ in items})
        print(f"{site} {product}: {len(items)} granules, track_dir_frame {tracks}", flush=True)
        if "--list" in sys.argv:
            continue
        (ROOT / "data" / site / product.lower()).mkdir(parents=True, exist_ok=True)
        rows = []
        # processes, not threads: h5py serialises every HDF5 call behind one lock, so
        # threads in one process read one request at a time
        with ProcessPoolExecutor(int(os.environ.get("NISAR_WORKERS", "8"))) as ex:
            for row in ex.map(partial(fetch, site, product), items):
                rows.append(row)
                print(product, row[1], row[2], row[0][17:80], flush=True)
        with open(ROOT / "data" / site / f"manifest_{product.lower()}.csv", "w", newline="") as f:
            csv.writer(f).writerows([["granule", "status", "valid_share"]] + rows)
