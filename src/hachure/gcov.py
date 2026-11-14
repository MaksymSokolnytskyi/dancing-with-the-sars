"""NISAR L2 GCOV (radiometrically terrain-corrected backscatter, gamma0) - window reads.

Same idea as GUNW: search ASF, read only the plate's window of each ~2-7 GB HDF5 over HTTP, cache it as a small .npz.
Grids: frequencyA (fine, 10-20 m) and frequencyB (coarse, ~80 m, about 8x smaller). Diagonal terms are power (linear);
dB = 10*log10. The `mask` layer: 0 = invalid, 255 = fill, 1..N = valid (subswath number; checked against the file's
attribute text). Noise is NOT subtracted (`noiseCorrectionApplied = False`), so calm water sits at the noise floor.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path

import numpy as np
import requests
from rasterio.warp import transform as warp_transform

from hachure.nisar import ASF_SEARCH, _epsg, _find, _incidence, _slice

ROOT = Path(__file__).resolve().parents[2]
CACHE = ROOT / "data/raw/gcov_cache"


@dataclass
class Scene:
    name: str
    url: str
    direction: str
    path: int
    frame: int
    date: datetime
    pols: str            # e.g. DHDH (dual-pol H transmit)


def search(lon: float, lat: float, *, direction: str | None = None, start: str | None = None,
           end: str | None = None, max_results: int = 500) -> list[Scene]:
    """GCOV scenes covering a point, one per (track, date), preferring PR over UR and the latest version."""
    params = {"dataset": "NISAR", "processingLevel": "GCOV", "intersectsWith": f"POINT({lon} {lat})",
              "maxResults": max_results, "output": "jsonlite"}
    if direction:
        params["flightDirection"] = direction.upper()
    if start:
        params["start"] = start
    if end:
        params["end"] = end
    r = requests.get(ASF_SEARCH, params=params, headers={"User-Agent": "space-apps-research"}, timeout=60)
    r.raise_for_status()
    best: dict[tuple, Scene] = {}
    for x in r.json()["results"]:
        if not x["downloadUrl"].endswith(".h5"):
            continue
        name = x["granuleName"]
        parts = name.split("_")
        date = datetime.strptime(x["startTime"][:19], "%Y-%m-%dT%H:%M:%S")
        s = Scene(name, x["downloadUrl"], x["flightDirection"], int(x["path"]), int(x["frame"]), date,
                  parts[9] if len(parts) > 9 else "?")
        key = (s.direction, s.path, s.frame, date.date())
        rank = ("_PR_" in name, "_F_F_" not in name, name[-3:])     # PR, not a "F" (diagnostic?) flag, latest version
        old = best.get(key)
        if old is None or rank > ("_PR_" in old.name, "_F_F_" not in old.name, old.name[-3:]):
            best[key] = s
    return sorted(best.values(), key=lambda s: (s.direction, s.path, s.frame, s.date))


@dataclass
class Window:
    hh: np.ndarray       # gamma0, linear power
    hv: np.ndarray | None
    mask: np.ndarray     # uint8 from the file
    looks: np.ndarray
    x: np.ndarray
    y: np.ndarray
    epsg: int
    date: datetime
    name: str
    # on the window grid: vv, vh (quad-pol only), inc (incidence, deg), neb_<pol> (noise-equivalent backscatter,
    # linear), g2s (gamma0 -> sigma0 factor); absent keys = not in the file or an old cache
    extra: dict = field(default_factory=dict)
    meta: dict = field(default_factory=dict)   # calibration scalars and processing flags (see _meta)

    def save(self, path) -> None:
        np.savez_compressed(path, hh=self.hh, hv=self.hv if self.hv is not None else np.zeros(0), mask=self.mask,
                            looks=self.looks, x=self.x, y=self.y, epsg=self.epsg, date=self.date.isoformat(),
                            name=self.name, meta=json.dumps(self.meta),
                            **{f"x_{k}": v for k, v in self.extra.items()})

    @classmethod
    def load(cls, path) -> "Window":
        d = np.load(path)
        hv = d["hv"] if d["hv"].size else None
        extra = {k[2:]: d[k] for k in d.files if k.startswith("x_")}
        meta = json.loads(str(d["meta"])) if "meta" in d.files else {}
        return cls(d["hh"], hv, d["mask"], d["looks"], d["x"], d["y"], int(d["epsg"]),
                   datetime.fromisoformat(str(d["date"])), str(d["name"]), extra, meta)

    def db(self, term: str = "hv") -> np.ndarray:
        v = {"hh": self.hh, "hv": self.hv}.get(term)
        v = (self.extra[term] if v is None else v).astype("float64")
        ok = np.isfinite(v) & (v > 0) & (self.mask > 0) & (self.mask < 255)
        out = np.full(v.shape, np.nan)
        out[ok] = 10 * np.log10(v[ok])
        return out


def read_window(h5, bbox, *, freq: str = "B", name: str = "", date: datetime | None = None,
                height_m: float | None = None) -> Window:
    g = f"science/LSAR/GCOV/grids/frequency{freq}"
    xs, ys = h5[f"{g}/xCoordinates"][()], h5[f"{g}/yCoordinates"][()]
    epsg = _epsg(h5, g)
    w, s, e, n = bbox
    bx, by = warp_transform("EPSG:4326", f"EPSG:{epsg}", [w, e, w, e], [s, s, n, n])
    cs, rs = _slice(xs, min(bx), max(bx)), _slice(ys, min(by), max(by))
    hh = h5[f"{g}/HHHH"][rs, cs] if f"{g}/HHHH" in h5 else h5[f"{g}/VVVV"][rs, cs]
    hv = h5[f"{g}/HVHV"][rs, cs] if f"{g}/HVHV" in h5 else None
    x, y = xs[cs], ys[rs]
    extra = {}
    for key, term in (("vv", "VVVV"), ("vh", "VHVH")):
        if f"{g}/{term}" in h5 and f"{g}/HHHH" in h5:
            extra[key] = h5[f"{g}/{term}"][rs, cs].astype("float32")
    if f"{g}/rtcGammaToSigmaFactor" in h5:
        extra["g2s"] = h5[f"{g}/rtcGammaToSigmaFactor"][rs, cs].astype("float32")
    if _find(h5, "incidenceAngle", "radarGrid"):
        extra["inc"] = _incidence(h5, x, y, height_m).astype("float32")
    neb = f"science/LSAR/GCOV/metadata/calibrationInformation/frequency{freq}/noiseEquivalentBackscatter"
    for pol in ("HH", "HV", "VH", "VV"):
        if f"{neb}/{pol}" in h5:
            extra[f"neb_{pol.lower()}"] = _on_grid(h5, f"{neb}/{pol}", x, y)
    return Window(hh.astype("float32"), None if hv is None else hv.astype("float32"),
                  h5[f"{g}/mask"][rs, cs], h5[f"{g}/numberOfLooks"][rs, cs].astype("float32"),
                  x, y, epsg, date or datetime.min, name, extra, _meta(h5))


def _on_grid(h5, path: str, x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """A 2-D metadata layer on its coarse grid (sibling x/yCoordinates), bilinear onto the window pixels."""
    from scipy.interpolate import RegularGridInterpolator
    grp = path.rsplit("/", 1)[0]
    cx, cy, lay = h5[f"{grp}/xCoordinates"][()], h5[f"{grp}/yCoordinates"][()], h5[path][()].astype("float64")
    if cy[0] > cy[-1]:
        cy, lay = cy[::-1], lay[::-1]
    f = RegularGridInterpolator((cy, cx), lay, bounds_error=False, fill_value=None)
    yy, xx = np.meshgrid(y, x, indexing="ij")
    return f(np.column_stack([yy.ravel(), xx.ravel()])).reshape(yy.shape).astype("float32")


def _meta(h5) -> dict:
    """Calibration scalars and processing flags, plus checks on the calibration layers (for the product audit)."""
    m = "science/LSAR/GCOV/metadata/"
    out = {}
    for k in ("noiseCorrectionApplied", "faradayRotationApplied", "rfiMitigationApplied", "isFullCovariance",
              "polarimetricSymmetrizationApplied"):
        p = m + "processingInformation/parameters/" + k
        if p in h5:
            out[k] = h5[p][()].decode()
    for p in ("productVersion", "compositeReleaseId", "isDithered"):
        if f"science/LSAR/identification/{p}" in h5:
            out[p] = h5[f"science/LSAR/identification/{p}"][()].decode()
    c = m + "calibrationInformation/"
    for f in ("A", "B"):
        if f"{c}frequency{f}" not in h5:
            continue
        out[f"faraday_{f}"] = float(h5[f"{c}frequency{f}/faradayRotation"][()])
        for pol in ("HH", "HV", "VH", "VV"):
            p = f"{c}frequency{f}/{pol}/rfiLikelihood"
            if p in h5:
                out[f"rfi_{f}_{pol}"] = float(h5[p][()])
        ap = f"{c}frequency{f}/elevationAntennaPattern/HH"
        if ap in h5:
            a = np.abs(h5[ap][()])
            a = a[np.isfinite(a)]
            out[f"antpat_{f}_all_one"] = bool(a.size and np.allclose(a, 1.0))
    xt = c + "crosstalk/txHorizontalCrosspol"
    if xt in h5:
        out["crosstalk_all_zero"] = bool(np.all(np.nan_to_num(np.abs(h5[xt][()])) == 0))
    return out


def cache_dir(bbox, freq: str = "B") -> Path:
    return CACHE / ("bbox_" + "_".join(f"{v:.3f}" for v in bbox) + f"_f{freq}")


def fetch(scenes: list[Scene], bbox, *, freq: str = "B", height_m: float | None = None, log=print) -> list[Path]:
    """Read and cache each scene's window (skips cached ones). Returns the cache paths that exist."""
    import time
    from hachure import nisar
    d = cache_dir(bbox, freq)
    d.mkdir(parents=True, exist_ok=True)
    out = []
    for i, s in enumerate(scenes, 1):
        p = d / f"{s.name}.npz"
        if not p.exists():
            t = time.time()
            try:
                with nisar.open_h5(s.url) as h5:
                    win = read_window(h5, bbox, freq=freq, name=s.name, date=s.date, height_m=height_m)
                tmp = d / f"_part_{i}.npz"
                win.save(tmp)
                tmp.replace(p)
                log(f"  [{i}/{len(scenes)}] {s.date:%Y-%m-%d} {p.stat().st_size / 2**20:.1f} MB {time.time() - t:.0f}s")
            except Exception as e:                                   # noqa: BLE001
                log(f"  [{i}/{len(scenes)}] FAILED {s.name[:60]}: {e}")
                continue
        out.append(p)
    return out
