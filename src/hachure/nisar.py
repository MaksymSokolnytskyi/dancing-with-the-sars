"""NISAR L2 GUNW (geocoded unwrapped interferogram) loader.

search()      ASF search API (no auth) -> list[Granule]
open_h5()     local path, or HTTPS streamed with an Earthdata token (env EARTHDATA_TOKEN);
              only the HDF5 chunks we slice are fetched, not the ~2.3 GB file.
read_pair()   window over a lon/lat bbox -> vertical displacement (m), coherence, mask
stack_rate()  many pairs -> mean vertical rate (m/yr) + per-pair displacements for animation

HDF5 paths are DISCOVERED by dataset name (unwrappedPhase, coherenceMagnitude, ...),
not hard-coded, because the product spec is still provisional. Run
`scripts/nisar_gunw_rate.py --inspect` on the first real file and check the tree.

ASSUMPTIONS TO VERIFY ON REAL DATA (marked VERIFY below):
- LOS sign: d_los = LOS_SIGN * lambda/(4 pi) * phase, LOS_SIGN=-1 -> positive = toward
  satellite (uplift). Mexico City must come out NEGATIVE; if not, use flip_sign.
- incidence angle cube at metadata/radarGrid/incidenceAngle (height, y, x), degrees.
"""
from __future__ import annotations

import os
import re
from dataclasses import dataclass
from datetime import datetime

import h5py
import numpy as np
import requests
from rasterio.warp import transform as warp_transform

ASF_SEARCH = "https://api.daac.asf.alaska.edu/services/search/param"
C = 299_792_458.0
DEFAULT_WAVELENGTH = C / 1.2575e9      # L-band, used only if the file has no centerFrequency
DEFAULT_INCIDENCE_DEG = 37.0
LOS_SIGN = -1.0                         # VERIFY (see module docstring)

_DATES = re.compile(r"_(\d{8}T\d{6})_(\d{8}T\d{6})_(\d{8}T\d{6})_(\d{8}T\d{6})_")


@dataclass
class Granule:
    name: str
    url: str
    direction: str          # ASCENDING / DESCENDING
    path: int
    frame: int
    ref_date: datetime
    sec_date: datetime
    size_bytes: int | None

    @property
    def days(self) -> float:
        return (self.sec_date - self.ref_date).total_seconds() / 86400


def _parse_dates(name: str):
    m = _DATES.search(name)
    if not m:
        raise ValueError(f"cannot parse dates from {name}")
    ref, _, sec, _ = (datetime.strptime(x, "%Y%m%dT%H%M%S") for x in m.groups())
    return ref, sec


def search(lon: float, lat: float, *, direction: str | None = None,
           start: str | None = None, end: str | None = None,
           max_results: int = 500) -> list[Granule]:
    """Find GUNW pairs covering a point. Dedupes PR/UR variants and versions per date pair."""
    params = {"dataset": "NISAR", "processingLevel": "GUNW",
              "intersectsWith": f"POINT({lon} {lat})", "maxResults": max_results,
              "output": "jsonlite"}
    if direction:
        params["flightDirection"] = direction.upper()
    if start:
        params["start"] = start
    if end:
        params["end"] = end
    r = requests.get(ASF_SEARCH, params=params, timeout=60)
    r.raise_for_status()
    best: dict[tuple, Granule] = {}
    for x in r.json()["results"]:
        if x.get("productType") != "GUNW" or not x["downloadUrl"].endswith(".h5"):
            continue
        ref, sec = _parse_dates(x["granuleName"])
        size = None
        sizes = (x.get("nisar") or {}).get("sizeMB") or {}
        if x["fileName"] in sizes:
            size = sizes[x["fileName"]].get("bytes")
        g = Granule(x["granuleName"], x["downloadUrl"], x["flightDirection"],
                    int(x["path"]), int(x["frame"]), ref, sec, size)
        key = (g.direction, g.path, g.frame, ref.date(), sec.date())
        # prefer standard "PR" over "UR", then the latest version string
        rank = ("_PR_" in g.name, g.name[-3:])
        if key not in best or rank > ("_PR_" in best[key].name, best[key].name[-3:]):
            best[key] = g
    return sorted(best.values(), key=lambda g: (g.direction, g.path, g.frame, g.ref_date))


def sequential_chain(granules: list[Granule]) -> list[Granule]:
    """Longest run of consecutive pairs (each secondary date = next reference date).

    In a chain the atmosphere of every intermediate date enters once with + and once
    with -, so it cancels in the cumulative sum; only the end dates' noise remains.
    """
    by_ref: dict = {}
    for g in sorted(granules, key=lambda g: g.days):          # prefer the shortest pair
        by_ref.setdefault(g.ref_date.date(), g)
    best: list[Granule] = []
    for start in sorted(by_ref):
        chain, d = [], start
        while d in by_ref:
            chain.append(by_ref[d])
            d = by_ref[d].sec_date.date()
        if len(chain) > len(best):
            best = chain
    return best


def dominant_track(granules: list[Granule]) -> list[Granule]:
    """Keep the (direction, path, frame) with the most pairs - LOS geometry must be consistent."""
    from collections import Counter
    if not granules:
        return []
    top = Counter((g.direction, g.path, g.frame) for g in granules).most_common(1)[0][0]
    return [g for g in granules if (g.direction, g.path, g.frame) == top]


# ---------------------------------------------------------------- opening files

def _signed_url(url: str, token: str) -> str:
    """ASF TEA: bearer token -> redirect to a short-lived presigned S3 URL (range-readable)."""
    r = requests.get(url, headers={"Authorization": f"Bearer {token}"},
                     allow_redirects=False, timeout=60)
    if r.status_code in (301, 302, 303, 307, 308):
        loc = r.headers["Location"]
        if "urs.earthdata.nasa.gov" in loc:
            raise PermissionError("Earthdata rejected the token (expired? NISAR access not granted?)")
        return loc
    if r.status_code == 200:
        return url
    raise PermissionError(f"ASF returned {r.status_code} for {url}")


TOKEN_FILE = os.path.join(os.path.expanduser("~"), ".edl_token")


def _token() -> str | None:
    """Env EARTHDATA_TOKEN, else ~/.edl_token (outside the repo, so never committed)."""
    tok = os.environ.get("EARTHDATA_TOKEN")
    if not tok and os.path.exists(TOKEN_FILE):
        with open(TOKEN_FILE, encoding="utf-8-sig") as fh:
            tok = fh.read().strip()
    return tok or None


def open_h5(src: str, token: str | None = None, block_mb: int = 8) -> h5py.File:
    """Local path or https URL. Token: argument, env EARTHDATA_TOKEN, or ~/.edl_token."""
    if not src.startswith("http"):
        return h5py.File(src, "r")
    import fsspec
    token = token or _token()
    if not token:
        raise PermissionError(f"No Earthdata token: put it in {TOKEN_FILE} or env EARTHDATA_TOKEN "
                              "(urs.earthdata.nasa.gov -> Generate Token)")
    f = fsspec.open(_signed_url(src, token), mode="rb", block_size=block_mb * 2**20,
                    cache_type="blockcache").open()
    return h5py.File(f, "r")


def tree(h5: h5py.File, max_lines: int = 400) -> str:
    """Printable tree of datasets (name, shape, dtype) for checking the real layout."""
    out: list[str] = []

    def visit(name, obj):
        if isinstance(obj, h5py.Dataset) and len(out) < max_lines:
            out.append(f"{name}  {obj.shape} {obj.dtype}")
    h5.visititems(visit)
    return "\n".join(out)


def _find(h5: h5py.File, leaf: str, must_contain: str = "") -> str | None:
    hits: list[str] = []

    def visit(name, obj):
        if isinstance(obj, h5py.Dataset) and name.split("/")[-1] == leaf and must_contain in name:
            hits.append(name)
    h5.visititems(visit)
    return sorted(hits, key=len)[0] if hits else None


# ---------------------------------------------------------------- reading a window

@dataclass
class PairWindow:
    disp_v: np.ndarray        # vertical displacement, m (NaN = invalid)
    coherence: np.ndarray
    x: np.ndarray             # projected x of columns (pixel centres)
    y: np.ndarray             # projected y of rows
    epsg: int
    days: float


def _epsg(h5, grid_group: str) -> int:
    for cand in (f"{grid_group}/projection", _find(h5, "projection")):
        if cand and cand in h5:
            d = h5[cand]
            for key in ("epsg_code", "epsg"):
                if key in d.attrs:
                    return int(np.atleast_1d(d.attrs[key])[0])
            try:
                return int(d[()])
            except (TypeError, ValueError):
                pass
    raise KeyError("no projection/EPSG found - inspect the tree")


def _slice(coords: np.ndarray, lo: float, hi: float) -> slice:
    idx = np.flatnonzero((coords >= lo) & (coords <= hi))
    if len(idx) == 0:
        raise ValueError("bbox does not overlap this granule")
    return slice(int(idx[0]), int(idx[-1]) + 1)


def _incidence(h5, x: np.ndarray, y: np.ndarray, height_m: float | None = None) -> np.ndarray:
    path = _find(h5, "incidenceAngle", "radarGrid")
    if not path:
        print("  (no incidenceAngle cube - using 37 deg)")
        return np.full((len(y), len(x)), DEFAULT_INCIDENCE_DEG)
    cube = h5[path]
    grp = path.rsplit("/", 1)[0]
    cx, cy = h5[f"{grp}/xCoordinates"][()], h5[f"{grp}/yCoordinates"][()]
    k = cube.shape[0] // 2
    if cube.ndim == 3 and height_m is not None and f"{grp}/heightAboveEllipsoid" in h5:
        k = int(np.argmin(np.abs(h5[f"{grp}/heightAboveEllipsoid"][()] - height_m)))
    layer = cube[k] if cube.ndim == 3 else cube[()]
    from scipy.interpolate import RegularGridInterpolator
    ys = cy if cy[0] < cy[-1] else cy[::-1]
    lay = layer if cy[0] < cy[-1] else layer[::-1]
    f = RegularGridInterpolator((ys, cx), lay, bounds_error=False, fill_value=None)
    yy, xx = np.meshgrid(y, x, indexing="ij")
    return f(np.column_stack([yy.ravel(), xx.ravel()])).reshape(yy.shape)


@dataclass
class RawWindow:
    """Everything read from one GUNW for a bbox, before any processing (cacheable)."""
    phase: np.ndarray         # unwrapped phase, rad (NaN = fill)
    iono: np.ndarray | None   # ionosphere phase screen, rad
    coherence: np.ndarray
    comp: np.ndarray          # connected components, 0 = invalid
    incidence: np.ndarray     # deg, at height_m
    x: np.ndarray
    y: np.ndarray
    epsg: int
    wavelength: float
    days: float

    def save(self, path) -> None:
        np.savez_compressed(path, **{k: (v if v is not None else np.array([]))
                                     for k, v in self.__dict__.items()})

    @classmethod
    def load(cls, path) -> "RawWindow":
        z = np.load(path)
        d = {k: z[k] for k in z.files}
        d["iono"] = d["iono"] if d["iono"].size else None
        d["epsg"], d["wavelength"], d["days"] = int(d["epsg"]), float(d["wavelength"]), float(d["days"])
        return cls(**d)


def read_raw(h5: h5py.File, bbox: tuple[float, float, float, float], days: float, *,
             pol: str = "HH", height_m: float | None = None) -> RawWindow:
    unw_path = _find(h5, "unwrappedPhase", f"/{pol}/") or _find(h5, "unwrappedPhase")
    if not unw_path:
        raise KeyError("no unwrappedPhase dataset - inspect the tree")
    pol_grp = unw_path.rsplit("/", 1)[0]
    grid_grp = pol_grp.rsplit("/", 1)[0]
    xs, ys = h5[f"{grid_grp}/xCoordinates"][()], h5[f"{grid_grp}/yCoordinates"][()]
    epsg = _epsg(h5, grid_grp)

    w, s, e, n = bbox
    bx, by = warp_transform("EPSG:4326", f"EPSG:{epsg}", [w, e, w, e], [s, s, n, n])
    cs = _slice(xs, min(bx), max(bx))
    rs = _slice(ys, min(by), max(by))
    x, y = xs[cs], ys[rs]

    phase = h5[unw_path][rs, cs].astype("float64")
    fill = h5[unw_path].attrs.get("_FillValue")
    if fill is not None and np.isfinite(np.atleast_1d(fill)[0]):
        phase[phase == np.atleast_1d(fill)[0]] = np.nan
    coh = h5[f"{pol_grp}/coherenceMagnitude"][rs, cs].astype("float64")
    cc_path = f"{pol_grp}/connectedComponents"
    if cc_path in h5:
        comp = h5[cc_path][rs, cs].astype(np.int32)
        cfill = h5[cc_path].attrs.get("_FillValue")
        if cfill is not None:
            comp[comp == int(np.atleast_1d(cfill)[0])] = 0
    else:
        comp = np.ones(phase.shape, np.int32)
    iono_path = f"{pol_grp}/ionospherePhaseScreen"
    iono = h5[iono_path][rs, cs].astype("float64") if iono_path in h5 else None

    freq_path = _find(h5, "centerFrequency", "frequencyA")
    lam = C / float(np.atleast_1d(h5[freq_path][()])[0]) if freq_path else DEFAULT_WAVELENGTH
    inc = _incidence(h5, x, y, height_m)
    return RawWindow(phase, iono, coh, comp, inc, x, y, epsg, lam, days)


def robust_screen(screen: np.ndarray, dx_m: float, block_m: float = 2000.0, win_blocks: int = 5) -> np.ndarray:
    """Median of block_m blocks, then a NaN-aware win_blocks x win_blocks median, nearest fill, bilinear upsample."""
    from scipy.interpolate import RegularGridInterpolator
    from scipy.ndimage import distance_transform_edt, generic_filter
    b = max(1, int(round(block_m / dx_m)))
    H, W = screen.shape
    nb, mb = -(-H // b), -(-W // b)
    pad = np.full((nb * b, mb * b), np.nan)
    pad[:H, :W] = screen
    blk = pad.reshape(nb, b, mb, b).transpose(0, 2, 1, 3).reshape(nb, mb, b * b)
    frac = np.isfinite(blk).mean(-1)
    with np.errstate(all="ignore"):
        med = np.nanmedian(blk, -1)
    med[frac < 0.25] = np.nan
    sm = generic_filter(med, lambda v: np.nanmedian(v) if np.isfinite(v).any() else np.nan, size=win_blocks,
                        mode="nearest")
    bad = ~np.isfinite(sm)
    if bad.all():
        return np.zeros_like(screen)
    if bad.any():
        idx = distance_transform_edt(bad, return_distances=False, return_indices=True)
        sm = sm[tuple(idx)]
    cy = (np.arange(nb) + 0.5) * b - 0.5
    cx = (np.arange(mb) + 0.5) * b - 0.5
    f = RegularGridInterpolator((cy, cx), sm, bounds_error=False, fill_value=None)
    yy, xx = np.mgrid[:H, :W]
    return f(np.column_stack([yy.ravel(), xx.ravel()])).reshape(H, W)


def process(raw: RawWindow, *, ref_bbox: tuple[float, float, float, float] | None = None,
            coh_min: float = 0.3, iono: bool | str = True, flip_sign: bool = False) -> "PairWindow":
    """Raw window -> vertical displacement (m).

    Checked on a real GUNW (2026-10-08/20, Mexico frame): subtracting ionospherePhaseScreen
    halves the phase scatter (4.76 -> 2.32 rad, corr 0.89), so iono=True is right.
    The tropospheric cubes at a single height did NOT help -> not applied (needs DEM heights).
    """
    phase = raw.phase.copy()
    if iono == "plane" and raw.iono is not None:
        # P-iono / H9: keep only the planar part of the per-pair iono screen (the ramps), drop its small-scale part
        ok = np.isfinite(raw.iono)
        yy, xx = np.mgrid[:raw.iono.shape[0], :raw.iono.shape[1]]
        A = np.column_stack([np.ones(ok.sum()), xx[ok], yy[ok]])
        c = np.linalg.lstsq(A, raw.iono[ok], rcond=None)[0]
        phase -= c[0] + c[1] * xx + c[2] * yy
    elif iono == "rplane" and raw.iono is not None:
        # H11: robust plane of the screen (islands of many cycles must not tilt the plane)
        ok = np.isfinite(raw.iono)
        yy, xx = np.mgrid[:raw.iono.shape[0], :raw.iono.shape[1]]
        keep = ok.copy()
        for _ in range(5):
            A = np.column_stack([np.ones(keep.sum()), xx[keep], yy[keep]])
            c = np.linalg.lstsq(A, raw.iono[keep], rcond=None)[0]
            r = raw.iono - (c[0] + c[1] * xx + c[2] * yy)
            mad = 1.4826 * np.nanmedian(np.abs(r[keep] - np.nanmedian(r[keep])))
            keep = ok & (np.abs(r - np.nanmedian(r[keep])) <= max(3 * mad, 2 * np.pi))
        phase -= c[0] + c[1] * xx + c[2] * yy
    elif iono == "robust" and raw.iono is not None:
        # H10: the screen has its own unwrapping islands; keep only its robust smooth part (10 km median of 2 km blocks)
        phase -= robust_screen(raw.iono, abs(float(raw.x[1] - raw.x[0])))
    elif iono and raw.iono is not None:
        phase -= raw.iono
    elif iono:
        print("  (no ionospherePhaseScreen - not corrected)")
    comp = raw.comp
    valid = np.isfinite(phase) & (comp > 0) & (raw.coherence >= coh_min)

    # reference region: pick its dominant component, zero its median phase
    yy, xx = np.meshgrid(raw.y, raw.x, indexing="ij")
    if ref_bbox:
        rx, ry = warp_transform("EPSG:4326", f"EPSG:{raw.epsg}",
                                [ref_bbox[0], ref_bbox[2]], [ref_bbox[1], ref_bbox[3]])
        ref = valid & (xx >= min(rx)) & (xx <= max(rx)) & (yy >= min(ry)) & (yy <= max(ry))
    else:
        ref = valid
    if ref.sum() < 20:
        raise ValueError("reference region has <20 valid pixels")
    main = np.bincount(comp[ref].ravel()).argmax()
    valid &= comp == main
    ref &= comp == main
    phase -= np.median(phase[ref])

    sign = -LOS_SIGN if flip_sign else LOS_SIGN
    d_los = sign * raw.wavelength / (4 * np.pi) * phase
    d_v = d_los / np.cos(np.radians(raw.incidence))      # assumes purely vertical motion
    d_v[~valid] = np.nan
    return PairWindow(d_v, raw.coherence, raw.x, raw.y, raw.epsg, raw.days)


def read_pair(h5: h5py.File, bbox, days: float, *, pol: str = "HH",
              height_m: float | None = None, **kw) -> "PairWindow":
    return process(read_raw(h5, bbox, days, pol=pol, height_m=height_m), **kw)


def fit_ramp(field: np.ndarray, iters: int = 4, cut: float = -0.03) -> np.ndarray:
    """Planar trend (residual long-wavelength atmosphere/iono/orbit) of a SMOOTHED rate field.

    Robust for one-sided signals like subsidence: each iteration refits the plane only on
    pixels whose residual is above `cut` (m/yr), so sinking bowls don't drag the plane down.
    Pass a smoothed field: on raw pixels the noise (~10 cm/yr for 12 pairs) exceeds `cut`
    and the plane creeps upward (learned the hard way on real NISAR data).
    Returns the plane; subtract it from both raw and smoothed fields.
    """
    m0 = np.isfinite(field)
    yy, xx = np.mgrid[:field.shape[0], :field.shape[1]]
    use = m0.copy()
    plane = np.zeros_like(field)
    for _ in range(iters):
        a = np.column_stack([np.ones(use.sum()), xx[use], yy[use]])
        c = np.linalg.lstsq(a, field[use], rcond=None)[0]
        plane = c[0] + c[1] * xx + c[2] * yy
        use = m0 & (field - plane > cut)
    return plane


def _common(axes: list[np.ndarray]) -> np.ndarray:
    common = np.round(axes[0], 3)
    for a in axes[1:]:
        common = np.intersect1d(common, np.round(a, 3))
    return common


def align(pairs: list[PairWindow]):
    """Put all pairs on the coordinates they share. Returns disps (N,H,W) m, cohs, x, y."""
    cx = _common([p.x for p in pairs])
    cy = _common([p.y for p in pairs])
    if len(cx) < 2 or len(cy) < 2:
        raise ValueError("pairs share no common grid - different frames? use dominant_track()")
    cy = cy[::-1] if pairs[0].y[0] > pairs[0].y[-1] else cy   # keep north-up order
    disps, cohs = [], []
    for p in pairs:
        ci = np.searchsorted(np.round(p.x, 3), cx)
        ry = np.round(p.y, 3)
        ri = (len(ry) - 1 - np.searchsorted(ry[::-1], cy)) if ry[0] > ry[-1] else np.searchsorted(ry, cy)
        disps.append(p.disp_v[np.ix_(ri, ci)])
        cohs.append(p.coherence[np.ix_(ri, ci)])
    return np.stack(disps), np.stack(cohs), cx, cy


def stack_rate(pairs: list[PairWindow], min_frac: float = 0.6):
    """Rate = sum of displacements / sum of their time spans, per pixel over valid pairs.

    For a sequential chain the intermediate dates' atmosphere cancels in the sum.
    Returns rate (m/yr), mean coherence, valid fraction, x, y, aligned disps (N,H,W).
    """
    disps, cohs, cx, cy = align(pairs)
    days = np.array([p.days for p in pairs])[:, None, None]
    ok = np.isfinite(disps)
    frac = ok.mean(axis=0)
    span = np.where(ok, days, 0).sum(axis=0)
    with np.errstate(all="ignore"):
        rate = np.where(ok, disps, 0).sum(axis=0) / span * 365.25
    rate[frac < min_frac] = np.nan
    return rate, cohs.mean(axis=0), frac, cx, cy, disps


def cumulative(disps: np.ndarray, days: list[float], rate: np.ndarray) -> np.ndarray:
    """Cumulative vertical displacement at each acquisition date of a chain (N+1, H, W), m.

    Gaps (invalid pixel in a pair) are filled with rate * dt so curves stay continuous;
    pixels with no rate stay NaN. First date = 0.
    """
    fill = rate[None] * (np.array(days)[:, None, None] / 365.25)
    d = np.where(np.isfinite(disps), disps, fill)
    cum = np.concatenate([np.zeros((1, *d.shape[1:])), np.cumsum(d, axis=0)])
    cum[:, ~np.isfinite(rate)] = np.nan
    return cum
