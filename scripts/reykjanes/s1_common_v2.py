"""Shared geometry for the second round of the 10 Nov 2023 dike analysis (graben mask + GNSS).

GeomV2 behaves like s1_grindavik_dike.Geometry (same .e/.n/.vec/.trk/.h/.y, same source_enu and clean), but
  - InSAR bins within MASK_KM of the dike trace are dropped (graben: low coherence, unwrapping errors,
    unmodelled normal faults); the trace is the uniform-dike NPE median (s1_grindavik_dike.npz)
  - GNSS (data/gnss/s1_gnss_2023.json) is appended as pseudo-bins: every valid component (E, N, U) of every
    station in every window is one point with a unit "look vector", and trk = window (0 = A016, 1 = D155),
    so the Svartsengi source automatically takes the volume change of the right window
  - whitening: InSAR with the exponential covariance per track (+ offset/plane marginalized), GNSS / sigma
  - noise: from the residual fit (cached) or, when available, from the quiet 2023 pairs (s1_noise_quiet.json)
"""
import json
from pathlib import Path

import numpy as np
from pyproj import Transformer

import s1_grindavik_dike as S
from reykjanes_analysis import PLACES

ML = S.ML
ROOT = Path(__file__).resolve().parents[2]
MASK_KM = 1.0
WIN = {"A016": 0, "D155": 1}


def trace_dist(e, n, e0, n0, strike, length):
    u = np.array([np.sin(np.deg2rad(strike)), np.cos(np.deg2rad(strike))])
    de, dn = e - e0, n - n0
    along = de * u[0] + dn * u[1]
    across = -de * u[1] + dn * u[0]
    return np.where(np.abs(along) <= length / 2, np.abs(across), np.hypot(np.abs(along) - length / 2, across))


def noise_info(source="auto"):
    """Exponential covariance parameters per track: quiet pairs if available (and asked), else residual fit (cached)."""
    tot = ML / "s1_noise_total.json"
    if source in ("auto", "total") and tot.exists():
        return json.loads(tot.read_text())["tracks"], "atmosphere (quiet pair) + model error"
    q = ML / "s1_noise_quiet.json"
    if source == "quiet" and q.exists():
        return json.loads(q.read_text())["tracks"], "quiet"
    c = ML / "s1_noise_residual.json"
    if not c.exists():
        G = S.Geometry()
        _, info, _ = S.fit_noise(G)
        c.write_text(json.dumps({k: {kk: float(vv) for kk, vv in v.items()} for k, v in info.items()}, indent=1))
    return json.loads(c.read_text()), "residual"


class GeomV2(S.Geometry):
    def __init__(self, gnss=True, mask_km=MASK_KM, noise="auto", trace=None):
        super().__init__()
        if not hasattr(self, "h"):                       # Mogi source: no elevations needed
            self.h = np.zeros(len(self.e))
        if trace is None:
            pm = np.median(np.load(ML / "s1_grindavik_dike.npz")["ps"], 0)
            trace = (pm[0], pm[1], pm[2], pm[3])
        self.trace = trace
        keep = trace_dist(self.e, self.n, *trace) > mask_km if mask_km > 0 else np.ones(len(self.e), bool)
        self.n_masked = int((~keep).sum())
        for k in ("e", "n", "vec", "y", "trk", "h"):
            setattr(self, k, getattr(self, k)[keep])
        self.sizes = [int((self.trk == k).sum()) for k in range(len(S.TRACKS))]
        self.n_ins = len(self.e)
        # GNSS pseudo-bins
        self.gnss = gnss
        gs, ge, gn, gv, gy, gt, gsig, glab, gsb = [], [], [], [], [], [], [], [], []
        if gnss:
            J = json.loads((ROOT / "data" / "gnss" / "s1_gnss_2023.json").read_text())
            for s, st in J.items():
                for w, d in st["windows"].items():
                    for c in range(3):
                        if not d["valid"][c]:
                            continue
                        ge.append(st["e_km"]); gn.append(st["n_km"]); gv.append(np.eye(3)[c]); gy.append(d["enu_cm"][c])
                        gt.append(WIN[w]); gsig.append(d["sigma_cm"][c]); glab.append(f"{s} {w} {'ENU'[c]}")
                        gsb.append(d.get("sigma_base_cm", d["sigma_cm"])[c])
        self.n_gnss = len(gy)
        self.gnss_label, self.gnss_sigma, self.gnss_sigma_base = glab, np.array(gsig), np.array(gsb)
        if self.n_gnss:
            import rasterio
            from fem_layered import DEM
            to_utm = Transformer.from_crs(4326, 32627, always_xy=True)
            sx, sy = to_utm.transform(*PLACES["Svartsengi"])
            lon, lat = Transformer.from_crs(32627, 4326, always_xy=True).transform(np.array(ge) * 1000 + sx, np.array(gn) * 1000 + sy)
            with rasterio.open(DEM) as src:
                h = np.array([v[0] for v in src.sample(zip(lon, lat))], float)
            h = np.clip(np.where(np.isfinite(h) & (h > -100), h, 0), 0, None) / 1000
            self.e = np.r_[self.e, ge]; self.n = np.r_[self.n, gn]; self.vec = np.r_[self.vec, np.array(gv)]
            self.y = np.r_[self.y, gy]; self.trk = np.r_[self.trk, gt]; self.h = np.r_[self.h, h]
        if noise == "ref":
            info, self.noise_source = {t: dict(REF_NOISE) for t in S.TRACKS}, "reference (whitening only; noise level unknown)"
            if self.n_gnss:
                self.gnss_sigma = self.gnss_sigma_base.copy()          # whiten GNSS with the base sigma; the 5 % term is simulated
        else:
            info, self.noise_source = noise_info(noise)
        self.noise_info = info
        covs = []
        for k, t in enumerate(S.TRACKS):
            m = np.arange(self.n_ins)[self.trk[:self.n_ins] == k]
            p = info[t]
            C = S.exp_cov(self.e[m], self.n[m], p["sill_cm2"], p["range_km"], p["white_cm2"] + 0.05)
            if "sill2_cm2" in p:                       # second exponential: model error
                D = np.hypot(self.e[m][:, None] - self.e[m][None], self.n[m][:, None] - self.n[m][None])
                C = C + p["sill2_cm2"] * np.exp(-D / p["range2_km"]) + p.get("white2_cm2", 0.0) * np.eye(len(m))
            covs.append(C)
        self.set_noise(covs)
        print(f"GeomV2: {self.n_ins} InSAR bins ({self.n_masked} masked within {mask_km} km of the trace), {self.n_gnss} GNSS values, "
              f"noise from {self.noise_source}", flush=True)

    def design(self, k):
        m = np.arange(self.n_ins)[self.trk[:self.n_ins] == k]
        return np.c_[np.ones(len(m)), self.e[m] / 10, self.n[m] / 10]

    def whiten(self, x):
        out, i0 = [], 0
        for n, W in zip(self.sizes, self.W):
            out.append(x[:, i0:i0 + n] @ W.T)
            i0 += n
        if self.n_gnss:
            out.append(x[:, self.n_ins:] / self.gnss_sigma)
        return np.concatenate(out, 1)

    def sample_noise(self, N, rng=None):
        rng = rng or S.RNG
        parts = []
        for k in range(len(S.TRACKS)):
            Lc = np.linalg.cholesky(self.cov[k] + 1e-9 * np.eye(self.sizes[k]))
            A = self.design(k)
            nu = rng.normal(0, 1, (N, 3)) * [S.SD_OFFSET / 3, S.SD_PLANE / 3, S.SD_PLANE / 3]
            parts.append(rng.normal(0, 1, (N, self.sizes[k])) @ Lc.T + nu @ A.T)
        if self.n_gnss:
            parts.append(rng.normal(0, 1, (N, self.n_gnss)) * self.gnss_sigma)
        return np.concatenate(parts, 1)

    def insar_order(self):
        """Indices that sort the InSAR part by track (the whitening blocks assume track order)."""
        return np.argsort(self.trk[:self.n_ins], kind="stable")


# ---------------------------------------------------------------- other Sentinel-1 pairs on the SAME cells
def product_dir(track, ref, sec):
    """HyP3 product directory for a pair given the two SLC granule names (or yyyymmdd strings)."""
    import glob
    a = ref[17:32] if ref.startswith("S1") else ref
    b = sec[17:32] if sec.startswith("S1") else sec
    d = [p for p in glob.glob(str(S.S1 / f"S1*_{a}*_{b}*")) if Path(p).is_dir()]
    if not d:
        raise FileNotFoundError(f"no HyP3 product {track} {a} {b}")
    return Path(d[0])


def pair_vector(G, dirs, min_px=20):
    """LOS (cm) of a pair binned onto G's InSAR cells; empty cells filled by inverse-distance weighting of the
    6 nearest filled cells. Returns (y_insar, fraction_filled)."""
    from scipy.spatial import cKDTree
    to_utm = Transformer.from_crs(4326, 32627, always_xy=True)
    sx, sy = to_utm.transform(*PLACES["Svartsengi"])
    out, filled = np.empty(G.n_ins), 0
    for k, t in enumerate(S.TRACKS):
        m = np.flatnonzero(G.trk[:G.n_ins] == k)
        cells = np.c_[np.floor(G.e[m] / S.BIN), np.floor(G.n[m] / S.BIN)].astype(int)
        los, x, y, epsg = S.read(dirs[t], "los_disp")
        corr = S.read(dirs[t], "corr")[0]
        ux, uy = Transformer.from_crs(epsg, 32627, always_xy=True).transform(x, y)
        e, n = (ux - sx) / 1000, (uy - sy) / 1000
        ok = np.isfinite(los) & (los != 0) & (corr > 0.35)
        key = np.c_[np.floor(e[ok] / S.BIN), np.floor(n[ok] / S.BIN)].astype(int)
        v = los[ok] * 100
        d = {}
        order = np.lexsort((key[:, 1], key[:, 0]))
        key, v = key[order], v[order]
        brk = np.flatnonzero(np.any(np.diff(key, axis=0) != 0, axis=1)) + 1
        for kk, vv in zip(np.split(key, brk), np.split(v, brk)):
            if len(vv) >= min_px:
                d[(kk[0, 0], kk[0, 1])] = np.median(vv)
        vals = np.array([d.get((a, b), np.nan) for a, b in cells])
        miss = ~np.isfinite(vals)
        if miss.any() and (~miss).sum() >= 6:
            tree = cKDTree(np.c_[G.e[m][~miss], G.n[m][~miss]])
            dist, idx = tree.query(np.c_[G.e[m][miss], G.n[m][miss]], 6)
            w = 1 / np.maximum(dist, 0.3) ** 2
            vals[miss] = (vals[~miss][idx] * w).sum(1) / w.sum(1)
        filled += miss.sum()
        out[m] = vals
    return out, filled / G.n_ins


def gnss_vector(G, ref_dates):
    """GNSS values in G's order for windows given as {track: (ref 'YYYY-MM-DD', sec 'YYYY-MM-DD')}.
    Evening path 16 (18:59): before = ref-1..ref, after = sec..sec+1; morning path 155 (07:59): before = ref-2..ref-1,
    after = sec-1..sec. Same detrending as s1_gnss_2023 (trend 1 Jan - 20 Oct 2023). Missing -> NaN."""
    from datetime import datetime, timedelta
    from reykjanes_gnss import read_tenv3
    import s1_gnss_2023 as GN
    cache = {}
    vals = np.full(G.n_gnss, np.nan)
    for i, lab in enumerate(G.gnss_label):
        s, w, c = lab.split()
        if s not in cache:
            t, x, _ = read_tenv3(s)
            dec = np.array([(a - datetime(2023, 1, 1)).days for a in t], float)
            fit = (t >= datetime(2023, 1, 1)) & (t <= datetime(2023, 10, 20))
            coef = np.linalg.lstsq(np.c_[np.ones(fit.sum()), dec[fit]], x[fit], rcond=None)[0]
            cache[s] = (t, x - np.c_[np.ones_like(dec), dec] @ coef)
        t, xd = cache[s]
        r, q = (datetime.fromisoformat(z) for z in ref_dates[w])
        if w == "A016":
            b, a = (r - timedelta(days=1), r), (q, q + timedelta(days=1))
        else:
            b, a = (r - timedelta(days=2), r - timedelta(days=1)), (q - timedelta(days=1), q)
        mb = np.array([b[0] <= z <= b[1] for z in t]); ma = np.array([a[0] <= z <= a[1] for z in t])
        if mb.any() and ma.any():
            vals[i] = (xd[ma].mean(0) - xd[mb].mean(0))["ENU".index(c)]
    return vals


# ---------------------------------------------------------------- noise of UNKNOWN level (hierarchical), for NPE v3
REF_NOISE = dict(sill_cm2=0.2, range_km=8.0, white_cm2=0.05)     # reference covariance used only for whitening
HIER_SILL = (0.1, 40.0)                                           # cm2, log-uniform: quiet pairs ~0.2-0.5, Nov 2023 ~30
HIER_RANGES = (3.0, 5.0, 8.0, 12.0, 16.0)                         # km


class HierNoise:
    """Per-simulation atmospheric noise with its own level and correlation length + offset/plane per track;
    GNSS: base sigma + 5 % of the (clean) signal."""
    def __init__(self, G):
        self.G = G
        self.L = {}
        for k in range(len(S.TRACKS)):
            m = np.flatnonzero(G.trk[:G.n_ins] == k)
            D = np.hypot(G.e[m][:, None] - G.e[m][None], G.n[m][:, None] - G.n[m][None])
            self.L[k] = [np.linalg.cholesky(np.exp(-D / r) + 1e-6 * np.eye(len(m))).astype(np.float32) for r in HIER_RANGES]

    def __call__(self, clean, rng):
        G, N = self.G, len(clean)
        out = np.zeros_like(clean, dtype=np.float32)
        sill = np.exp(rng.uniform(np.log(HIER_SILL[0]), np.log(HIER_SILL[1]), N)).astype(np.float32)
        white = sill * rng.uniform(0, 0.2, N).astype(np.float32)
        ri = rng.integers(0, len(HIER_RANGES), N)
        i0 = 0
        for k in range(len(S.TRACKS)):
            n = G.sizes[k]
            A = G.design(k)
            z = rng.standard_normal((N, n)).astype(np.float32)
            for j, Lc in enumerate(self.L[k]):
                s = ri == j
                if s.any():
                    out[s, i0:i0 + n] = (z[s] @ Lc.T) * np.sqrt(sill[s])[:, None]
            out[:, i0:i0 + n] += rng.standard_normal((N, n)).astype(np.float32) * np.sqrt(white)[:, None]
            nu = rng.normal(0, 1, (N, 3)) * [S.SD_OFFSET / 3, S.SD_PLANE / 3, S.SD_PLANE / 3]
            out[:, i0:i0 + n] += (nu @ A.T).astype(np.float32)
            i0 += n
        if G.n_gnss:
            sb = G.gnss_sigma_base
            sd = np.sqrt(sb[None] ** 2 + (0.05 * clean[:, G.n_ins:]) ** 2)
            out[:, G.n_ins:] = rng.standard_normal((N, G.n_gnss)) * sd
        return out, sill
