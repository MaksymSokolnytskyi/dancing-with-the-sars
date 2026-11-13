"""NISAR L2 GCOV (radiometrically terrain-corrected backscatter, gamma0) - window reads.

Same idea as GUNW: search ASF, read only the plate's window of each ~2-7 GB HDF5 over HTTP, cache it as a small .npz.
Grids: frequencyA (fine, 10-20 m) and frequencyB (coarse, ~80 m, about 8x smaller). Diagonal terms are power (linear);
dB = 10*log10. The `mask` layer: 0 = invalid, 255 = fill, 1..N = valid (subswath count) - VERIFY against the product spec.
"""
from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

import numpy as np
import requests
from rasterio.warp import transform as warp_transform

from hachure.nisar import ASF_SEARCH, _epsg, _slice

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

    def save(self, path) -> None:
        np.savez_compressed(path, hh=self.hh, hv=self.hv if self.hv is not None else np.zeros(0), mask=self.mask,
                            looks=self.looks, x=self.x, y=self.y, epsg=self.epsg, date=self.date.isoformat(),
                            name=self.name)

    @classmethod
    def load(cls, path) -> "Window":
        d = np.load(path)
        hv = d["hv"] if d["hv"].size else None
        return cls(d["hh"], hv, d["mask"], d["looks"], d["x"], d["y"], int(d["epsg"]),
                   datetime.fromisoformat(str(d["date"])), str(d["name"]))

    def db(self, term: str = "hv") -> np.ndarray:
        v = (self.hv if term == "hv" else self.hh).astype("float64")
        ok = np.isfinite(v) & (v > 0) & (self.mask > 0) & (self.mask < 255)
        out = np.full(v.shape, np.nan)
        out[ok] = 10 * np.log10(v[ok])
        return out


def read_window(h5, bbox, *, freq: str = "B", name: str = "", date: datetime | None = None) -> Window:
    g = f"science/LSAR/GCOV/grids/frequency{freq}"
    xs, ys = h5[f"{g}/xCoordinates"][()], h5[f"{g}/yCoordinates"][()]
    epsg = _epsg(h5, g)
    w, s, e, n = bbox
    bx, by = warp_transform("EPSG:4326", f"EPSG:{epsg}", [w, e, w, e], [s, s, n, n])
    cs, rs = _slice(xs, min(bx), max(bx)), _slice(ys, min(by), max(by))
    hh = h5[f"{g}/HHHH"][rs, cs] if f"{g}/HHHH" in h5 else h5[f"{g}/VVVV"][rs, cs]
    hv = h5[f"{g}/HVHV"][rs, cs] if f"{g}/HVHV" in h5 else None
    return Window(hh.astype("float32"), None if hv is None else hv.astype("float32"),
                  h5[f"{g}/mask"][rs, cs], h5[f"{g}/numberOfLooks"][rs, cs].astype("float32"),
                  xs[cs], ys[rs], epsg, date or datetime.min, name)


def cache_dir(bbox, freq: str = "B") -> Path:
    return CACHE / ("bbox_" + "_".join(f"{v:.3f}" for v in bbox) + f"_f{freq}")


def fetch(scenes: list[Scene], bbox, *, freq: str = "B", log=print) -> list[Path]:
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
                    win = read_window(h5, bbox, freq=freq, name=s.name, date=s.date)
                tmp = d / f"_part_{i}.npz"
                win.save(tmp)
                tmp.replace(p)
                log(f"  [{i}/{len(scenes)}] {s.date:%Y-%m-%d} {p.stat().st_size / 2**20:.1f} MB {time.time() - t:.0f}s")
            except Exception as e:                                   # noqa: BLE001
                log(f"  [{i}/{len(scenes)}] FAILED {s.name[:60]}: {e}")
                continue
        out.append(p)
    return out
