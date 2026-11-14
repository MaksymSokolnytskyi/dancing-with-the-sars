"""A real dike for the NPE: the 10 November 2023 Grindavik intrusion seen by Sentinel-1.

NISAR data start in November 2025, so no NISAR pair contains a dike yet. To show the real-time
dike inversion on a REAL event we use two Sentinel-1 interferograms processed on demand by
ASF HyP3 (INSAR_GAMMA, 10x2 looks = 40 m pixels), both spanning the intrusion of 10 Nov 2023:

  path 16 ascending   2023-10-31 -> 2023-11-12   (C-band, wavelength 5.55 cm)
  path 155 descending 2023-11-10 07:59 UTC (hours BEFORE the dike) -> 2023-11-22

Model (12 parameters): one vertical rectangular tensile dislocation (Okada 1985, homogeneous
half-space) with uniform opening + a point source near Svartsengi (the reservoir drained into the
dike). The two pairs cover different time windows, so its volume change is separate for each track.
Svartsengi source (--source):
  fem   (default) LAYERED-CRUST finite-element Green's functions with topography, the same physics
        that agreed best with MCMC for the 2025-2026 NISAR inflation (scripts/fem_layered.py:
        table output/reykjanes/ml/fem_green_table.npz, 1-10 km b.s.l. every 0.1 km, linear in depth;
        ASTER GDEM v3 elevations, variable-depth approximation); depth is below sea level
  mogi  homogeneous half-space Mogi (written to s1_grindavik_dike_mogi.*)
The parameter keys keep the name "Mogi ..." for both variants.

Steps
  download  fetch the finished HyP3 products (public URLs; job ids in data/sentinel1/jobs_submitted.json)
  prep      LOS displacement -> 1 km bins in the local frame (km from Svartsengi), LOS unit vectors from
            the HyP3 look-vector rasters; noise covariance (exponential + white) fitted to the residuals
            of a quick least-squares fit, offset + plane per track marginalized as in the NISAR code
  npe       100 000 simulations, whitening + PCA, neural spline flow (scripts/ml_npe_flow.train);
            calibration on 400 synthetic events; posterior for the REAL data, timed; MCMC on the same
            likelihood for comparison.

Writes output/reykjanes/ml/s1_grindavik_dike.json and .png.
Usage: python scripts/reykjanes/s1_grindavik_dike.py [download] [prep] [npe]   (default: all)
"""
from _secrets import read_env
import glob
import json
import subprocess
import sys
import time
import zipfile
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import torch
from pyproj import Transformer
from scipy.optimize import least_squares

from ml_mcmc_source import log_post, stretch_sampler
from ml_npe_flow import from_z, posterior, to_z, train
from okada import tensile
from progress import Progress
from reykjanes_analysis import PLACES

ROOT = Path(__file__).resolve().parents[2]
S1 = ROOT / "data" / "sentinel1"
ML = ROOT / "output" / "reykjanes" / "ml"
TRACKS = ("A016", "D155")
PAIRS = {"A016": "20231031T185905_20231112T185905", "D155": "20231110T075906_20231122T075906"}   # frame 204 (south) for A016
NAMES = ("dike east km", "dike north km", "strike deg", "length km", "top km", "height km", "log10 opening m",
         "Mogi east km", "Mogi north km", "Mogi depth km", "Mogi dV A016 1e6 m3", "Mogi dV D155 1e6 m3")
LO = np.array([-8, -8, 0, 2, 0, 1, -1.0, -4, -4, 2, -150, -150])
HI = np.array([8, 8, 180, 22, 3, 8, 1.0, 4, 4, 8, 30, 30])
BIN = 1.0           # km
BOX = (-20, 20, -14, 14)
SD_OFFSET, SD_PLANE = 30.0, 10.0      # cm: wide priors on the marginalized offset / plane per track
K = 256
SOURCE = sys.argv[sys.argv.index("--source") + 1] if "--source" in sys.argv else "fem"
TAG = "" if SOURCE == "fem" else "_" + SOURCE          # output file suffix
RNG = np.random.default_rng(23)
PUBLISHED = dict(length_km=15, volume_1e6m3=(130, 139), max_opening_m=8, svartsengi_dV_1e6m3=(-82, -76),
                 source="Sigmundsson et al. 2024, Science 383:1228 (doi 10.1126/science.adn2838) and EGU24-14829")


def arg(name, default):
    return type(default)(sys.argv[sys.argv.index(name) + 1]) if name in sys.argv else default


# ------------------------------------------------------------------ download
def download():
    jobs = json.loads((S1 / "jobs_submitted.json").read_text())["jobs"]
    env = read_env()
    tok = env["EARTHDATA_API_KEY"].strip().strip('"')
    for j in jobs:
        r = subprocess.run(["curl", "-s", "-H", f"Authorization: Bearer {tok}",
                            f"https://hyp3-api.asf.alaska.edu/jobs/{j['job_id']}"], capture_output=True, text=True)
        d = json.loads(r.stdout)
        if d["status_code"] != "SUCCEEDED":
            print(j["job_id"][:8], d["status_code"]); continue
        f = d["files"][0]
        z = S1 / f["filename"]
        if not z.exists():
            print("downloading", f["filename"], f"{f['size'] / 1e6:.0f} MB", flush=True)
            subprocess.run(["curl", "-s", "-L", "-o", str(z), f["url"]], check=True)    # public product URL
        with zipfile.ZipFile(z) as zf:
            zf.extractall(S1)
    del tok


# ------------------------------------------------------------------ prep
def product(track):
    a, b = PAIRS[track].split("_")
    dirs = [Path(p) for p in glob.glob(str(S1 / f"S1*_{a}_{b}_*")) if Path(p).is_dir()]
    if not dirs:
        raise FileNotFoundError(f"no HyP3 product for {track} {PAIRS[track]} in {S1}")
    return dirs[0]


def read(d, suffix):
    """Raster values and pixel-centre coordinates, cropped to the study box around Svartsengi."""
    import rasterio
    from rasterio.windows import from_bounds
    f = next(d.glob(f"*_{suffix}.tif"))
    with rasterio.open(f) as r:
        tr = Transformer.from_crs(4326, r.crs.to_epsg(), always_xy=True)
        cx, cy = tr.transform(*PLACES["Svartsengi"])
        m = 2000
        w = from_bounds(cx + BOX[0] * 1000 - m, cy + BOX[2] * 1000 - m, cx + BOX[1] * 1000 + m, cy + BOX[3] * 1000 + m, r.transform)
        w = w.round_offsets().round_lengths()
        v = r.read(1, window=w, boundless=True, fill_value=np.nan).astype(float)
        if r.nodata is not None:
            v[v == r.nodata] = np.nan
        T = r.window_transform(w)
        rows, cols = np.mgrid[0:v.shape[0], 0:v.shape[1]]
        x, y = rasterio.transform.xy(T, rows.ravel(), cols.ravel())
        return v.ravel(), np.asarray(x), np.asarray(y), r.crs.to_epsg()


def prep():
    to_utm = Transformer.from_crs(4326, 32627, always_xy=True)
    sx, sy = to_utm.transform(*PLACES["Svartsengi"])
    out = {}
    for t in TRACKS:
        d = product(t)
        los, x, y, epsg = read(d, "los_disp")
        corr = read(d, "corr")[0]
        th, ph = read(d, "lv_theta")[0], read(d, "lv_phi")[0]
        tr = Transformer.from_crs(epsg, 32627, always_xy=True)
        ux, uy = tr.transform(x, y)
        e, n = (ux - sx) / 1000, (uy - sy) / 1000
        ok = (np.isfinite(los) & (los != 0) & (corr > 0.35) & (e > BOX[0]) & (e < BOX[1]) & (n > BOX[2]) & (n < BOX[3]))
        e, n, los, th, ph = e[ok], n[ok], los[ok] * 100, th[ok], ph[ok]
        # HyP3 look vector (ground -> satellite): theta = elevation angle above the horizontal, phi = angle from east, anticlockwise
        vec = np.c_[np.cos(th) * np.cos(ph), np.cos(th) * np.sin(ph), np.sin(th)]
        key = np.c_[np.floor(e / BIN), np.floor(n / BIN)].astype(int)
        uk, inv, cnt = np.unique(key, axis=0, return_inverse=True, return_counts=True)
        inv = inv.ravel()
        good = np.flatnonzero(cnt >= 40)
        order = np.argsort(inv)
        starts = np.r_[0, np.cumsum(cnt)]
        rows = []
        for g in good:
            idx = order[starts[g]:starts[g + 1]]
            rows.append([np.median(e[idx]), np.median(n[idx]), np.median(los[idx]), *np.median(vec[idx], 0),
                         1.4826 * np.median(np.abs(los[idx] - np.median(los[idx])))])
        rows = np.array(rows)
        out[t] = rows
        print(f"{t}: {ok.sum()} coherent pixels -> {len(rows)} bins; LOS {np.percentile(rows[:, 2], [1, 50, 99]).round(1)} cm; "
              f"mean look vector {rows[:, 3:6].mean(0).round(2)}", flush=True)
    np.savez(S1 / "s1_bins.npz", **out)


class Geometry:
    def __init__(self):
        z = np.load(S1 / "s1_bins.npz")
        self.tr = {t: z[t] for t in TRACKS}
        self.sizes = [len(self.tr[t]) for t in TRACKS]
        self.e = np.concatenate([self.tr[t][:, 0] for t in TRACKS])
        self.n = np.concatenate([self.tr[t][:, 1] for t in TRACKS])
        self.vec = np.concatenate([self.tr[t][:, 3:6] for t in TRACKS])
        self.y = np.concatenate([self.tr[t][:, 2] for t in TRACKS])
        self.trk = np.repeat(np.arange(len(TRACKS)), self.sizes)
        if SOURCE == "fem":
            import rasterio
            from fem_layered import DEM
            t = np.load(ML / "fem_green_table.npz")
            self.fd, self.fr, self.fur, self.fuz = t["depths"], t["r"], t["ur"], t["uz"]     # m per m^3 on r (km)
            to_utm = Transformer.from_crs(4326, 32627, always_xy=True)
            sx, sy = to_utm.transform(*PLACES["Svartsengi"])
            lon, lat = Transformer.from_crs(32627, 4326, always_xy=True).transform(self.e * 1000 + sx, self.n * 1000 + sy)
            with rasterio.open(DEM) as src:
                h = np.array([v[0] for v in src.sample(zip(lon, lat))], float)
            h[(h < -100) | ~np.isfinite(h)] = 0.0
            self.h = np.clip(h, 0, None) / 1000                                               # km
            print(f"FE source: layered crust + topography ({Path(DEM).name}), elevations {np.percentile(self.h * 1000, [5, 50, 95]).round(0)} m", flush=True)

    def source_enu(self, p):
        """ENU displacement (cm) of the Svartsengi point source: layered FE Green's function or Mogi."""
        de, dn = self.e - p[7], self.n - p[8]
        dv = np.where(self.trk == 0, p[10], p[11])                      # 1e6 m3
        if SOURCE == "fem":
            D = p[9]
            f = np.clip((D - self.fd[0]) / (self.fd[1] - self.fd[0]), 0, len(self.fd) - 1.000001)
            i, w = int(f), f - int(f)
            Gr = (1 - w) * self.fur[i] + w * self.fur[i + 1]
            Gz = (1 - w) * self.fuz[i] + w * self.fuz[i + 1]
            sc = D / (D + self.h)
            r = np.hypot(de, dn)
            ur = sc ** 2 * np.interp(r * sc, self.fr, Gr) * dv * 1e6 * 100
            uz = sc ** 2 * np.interp(r * sc, self.fr, Gz) * dv * 1e6 * 100
            rr = np.maximum(r, 1e-6)
            return np.c_[ur * de / rr, ur * dn / rr, uz]
        R3 = (de ** 2 + dn ** 2 + p[9] ** 2) ** 1.5
        k = 3 / (4 * np.pi) * 1e5 * dv * 1e-3                           # 1e6 m3 -> km3 -> cm
        return np.c_[k * de / R3, k * dn / R3, k * p[9] / R3]

    def clean(self, q):
        out = np.empty((len(q), len(self.e)))
        for i, p in enumerate(q):
            u = tensile(self.e, self.n, p[0], p[1], p[4] + p[5] / 2, p[3], p[5], p[2], 90.0, 10 ** p[6]) * 100
            u = u + self.source_enu(p)
            out[i] = (u * self.vec).sum(1)
        return out

    def design(self, k):
        m = self.trk == k
        return np.c_[np.ones(m.sum()), self.e[m] / 10, self.n[m] / 10]

    def set_noise(self, cov):
        self.cov, self.W = cov, []
        for k in range(len(TRACKS)):
            A = self.design(k)
            S = cov[k] + A @ np.diag([SD_OFFSET ** 2, SD_PLANE ** 2, SD_PLANE ** 2]) @ A.T
            self.W.append(np.linalg.inv(np.linalg.cholesky(S)))

    def whiten(self, x):
        out, i0 = [], 0
        for n, W in zip(self.sizes, self.W):
            out.append(x[:, i0:i0 + n] @ W.T)
            i0 += n
        return np.concatenate(out, 1)

    def sample_noise(self, N):
        parts = []
        for k in range(len(TRACKS)):
            Lc = np.linalg.cholesky(self.cov[k] + 1e-9 * np.eye(self.sizes[k]))
            A = self.design(k)
            nu = RNG.normal(0, 1, (N, 3)) * [SD_OFFSET / 3, SD_PLANE / 3, SD_PLANE / 3]
            parts.append(RNG.normal(0, 1, (N, self.sizes[k])) @ Lc.T + nu @ A.T)
        return np.concatenate(parts, 1)


def exp_cov(e, n, sill, rng_km, white):
    D = np.hypot(e[:, None] - e[None], n[:, None] - n[None])
    return sill * np.exp(-D / rng_km) + white * np.eye(len(e))


def fit_noise(G):
    """Quick least-squares fit (white noise, offset + plane per track), then an exponential covariance fitted
    to the binned semivariogram of the residuals of each track."""
    def resid(p):
        pred = G.clean(p[None])[0]
        r, i0 = [], 0
        for k, n in enumerate(G.sizes):
            A = G.design(k)
            rr = G.y[i0:i0 + n] - pred[i0:i0 + n]
            c, *_ = np.linalg.lstsq(A, rr, rcond=None)
            r.append(rr - A @ c)
            i0 += n
        return np.concatenate(r)
    cand = LO + RNG.random((3000, 12)) * (HI - LO)
    cost = [np.sum(resid(c) ** 2) for c in cand]
    best = None
    with Progress("стартова підгонка МНК (для моделі шуму)", total=8, unit="старт") as pr:
        for c in cand[np.argsort(cost)[:8]]:
            f = least_squares(resid, c, bounds=(LO, HI), x_scale=(HI - LO) / 10, max_nfev=400)
            if best is None or f.cost < best.cost:
                best = f
            pr.update(1, note=f"найкраща сума квадратів {2 * best.cost:.0f} см²")
    r = resid(best.x)
    print("least-squares start:", dict(zip(NAMES, best.x.round(2))), f"rms {np.sqrt(np.mean(r ** 2)):.2f} cm", flush=True)
    covs, info, i0 = [], {}, 0
    for k, n in enumerate(G.sizes):
        rr, e, nn = r[i0:i0 + n], G.e[i0:i0 + n], G.n[i0:i0 + n]
        D = np.hypot(e[:, None] - e[None], nn[:, None] - nn[None])
        V = 0.5 * (rr[:, None] - rr[None]) ** 2
        edges = np.arange(0, 16.5, 1.0)
        h = 0.5 * (edges[1:] + edges[:-1])
        g = np.array([V[(D > a) & (D <= b)].mean() for a, b in zip(edges[:-1], edges[1:])])
        def f_res(p):
            s, L, w = np.exp(p)
            return (w + s * (1 - np.exp(-h / L))) - g
        p = np.exp(least_squares(f_res, np.log([np.var(rr), 3.0, 0.3 * np.var(rr)])).x)
        covs.append(exp_cov(e, nn, *p))
        info[TRACKS[k]] = dict(sill_cm2=round(p[0], 2), range_km=round(p[1], 2), white_cm2=round(p[2], 2))
        i0 += n
    print("noise covariance:", info, flush=True)
    return covs, info, best.x


class Model:
    lo, hi = LO, HI

    def __init__(self, G):
        self.G = G

    def predict(self, q):
        return self.G.clean(q)

    def wrap(self, q):
        q = q.copy()
        q[:, 2] = np.mod(q[:, 2], 180)
        return q


def volume(s):
    return s[..., 3] * s[..., 5] * 10 ** s[..., 6]          # km * km * m = 1e6 m3


def npe():
    G = Geometry()
    covs, noise_info, lsq = fit_noise(G)
    G.set_noise(covs)
    T = {}
    n_sims = arg("--sims", 100000)
    yobs = G.y[None]

    def one_round(lo, hi, label, seed):
        t0 = time.time()
        q = lo + RNG.random((n_sims, 12)) * (hi - lo)
        X = np.empty((n_sims, len(G.e)))
        with Progress(f"симуляції Окади + Моґі, {label} ({n_sims})", total=n_sims, unit="сим") as pr:
            for j in range(0, n_sims, 5000):
                X[j:j + 5000] = G.clean(q[j:j + 5000])
                pr.update(5000)
        S = G.whiten(X)
        del X
        ts = time.time() - t0
        _, sv, Vt = np.linalg.svd(S[:20000] - S[:20000].mean(0), full_matrices=False)
        U = Vt[:K].T
        SP = S @ U
        scale = SP.std(0) + 1.0
        print(f"[{label}] {n_sims} simulations in {ts:.0f} s; {K} components keep {(sv[:K] ** 2).sum() / (sv ** 2).sum():.4%} "
              f"of the signal variance; component {K} = {sv[K - 1] / np.sqrt(20000):.2f} x noise", flush=True)
        del S
        net, tt = train(SP, to_z(q, LO, HI), None, arg("--epochs", 30), f"дайка Sentinel-1, {label}", seed, scale=scale)
        t0 = time.time()
        ps = posterior([net], G.whiten(yobs) @ U, LO, HI, 6000, scale=scale)
        ti = time.time() - t0
        return dict(net=net, U=U, scale=scale, ps=ps, sim=ts, train=tt, infer=ti, lo=lo, hi=hi)

    R1 = one_round(LO, HI, "раунд 1", 7)
    w = HI - LO
    a, b = np.quantile(R1["ps"], [0.0005, 0.9995], 0)
    lo2, hi2 = np.maximum(LO, a - 0.25 * (b - a)), np.minimum(HI, b + 0.25 * (b - a))
    if lo2[2] - LO[2] < 0.02 * w[2] or HI[2] - hi2[2] < 0.02 * w[2]:      # strike is periodic
        lo2[2], hi2[2] = LO[2], HI[2]
    print("truncated box: " + ", ".join(f"{nm} {l:.2f}..{h:.2f}" for nm, l, h in zip(NAMES, lo2, hi2)), flush=True)
    R = one_round(lo2, hi2, "раунд 2", 8)
    net, U, scale = R["net"], R["U"], R["scale"]
    T.update(simulate=R1["sim"] + R["sim"], train=R1["train"] + R["train"], infer_real=R["infer"])
    # calibration on held-out simulations from the truncated box (Gaussian noise of the fitted covariance + offset/plane)
    nt = 400
    qt = lo2 + RNG.random((nt, 12)) * (hi2 - lo2)
    xt = G.clean(qt) + G.sample_noise(nt)
    Xt = torch.tensor(G.whiten(xt) @ U / scale, dtype=torch.float32)
    with torch.no_grad():
        st = np.concatenate([from_z(net.flow(net.emb(Xt[j:j + 50])).sample((1000,)).numpy().astype(float), LO, HI)
                             for j in range(0, nt, 50)], axis=1)
    cal = {}
    for lev in (.68, .95):
        a, b = np.quantile(st, [(1 - lev) / 2, (1 + lev) / 2], 0)
        cal[lev] = ((qt >= a) & (qt <= b)).mean(0).round(3).tolist()
    print("calibration 68 %:", cal[.68], "\n            95 %:", cal[.95], flush=True)
    ps = R["ps"]
    pm = np.median(ps, 0)
    print(f"REAL Sentinel-1 posterior in {T['infer_real']:.2f} s", flush=True)
    for i, nm in enumerate(NAMES):
        print(f"   {nm:22s} round 1 {np.percentile(R1['ps'][:, i], [16, 50, 84]).round(2)}   round 2 {np.percentile(ps[:, i], [16, 50, 84]).round(2)}")
    V = volume(ps)
    print(f"   dike volume 1e6 m3     {np.percentile(V, [16, 50, 84]).round(0)}  (published {PUBLISHED['volume_1e6m3']})")
    # ------------- MCMC on the same likelihood, started from the best of 4000 prior draws
    model = Model(G)
    yw = G.whiten(G.y[None])[0]
    like = lambda pred: -0.5 * ((G.whiten(pred) - yw) ** 2).sum(1)
    cand = LO + RNG.random((4000, 12)) * (HI - LO)
    lp = log_post(model, like, cand)
    steps = arg("--mcmc-steps", 3000)
    t0 = time.time()
    chain, lps, acc = stretch_sampler(model, like, cand[np.argsort(lp)[-48:]], steps, 48, RNG, label="дайка Sentinel-1")
    T["mcmc"] = time.time() - t0
    post = chain[steps // 2:]
    mlp = lps[steps // 2:].mean(0)
    pmc = post[:, mlp > mlp.max() - 15].reshape(-1, 12)
    lp_npe, lp_mc = like(G.clean(pm[None]))[0], like(G.clean(np.median(pmc, 0)[None]))[0]
    print(f"MCMC {T['mcmc']:.0f} s, acceptance {acc:.2f}; log-likelihood at medians: NPE {lp_npe:.1f}, MCMC {lp_mc:.1f}", flush=True)
    for i, nm in enumerate(NAMES):
        print(f"   {nm:22s} NPE {np.median(ps[:, i]):8.2f}   MCMC {np.median(pmc[:, i]):8.2f}")
    T = {k: round(v, 3) for k, v in T.items()}
    pct = lambda s: {nm: np.percentile(s[:, i], [16, 50, 84]).round(3).tolist() for i, nm in enumerate(NAMES)}
    out = dict(source=SOURCE, pairs=PAIRS, bins={t: n for t, n in zip(TRACKS, G.sizes)}, noise=noise_info, n_sims=n_sims, k_pca=K,
               timing=T, calibration=cal, npe=pct(ps), npe_round1=pct(R1['ps']), box=[lo2.tolist(), hi2.tolist()], mcmc=pct(pmc),
               npe_volume=np.percentile(V, [16, 50, 84]).round(1).tolist(),
               mcmc_volume=np.percentile(volume(pmc), [16, 50, 84]).round(1).tolist(),
               opening_m=np.percentile(10 ** ps[:, 6], [16, 50, 84]).round(2).tolist(),
               loglike=dict(npe=round(float(lp_npe), 1), mcmc=round(float(lp_mc), 1)),
               published=PUBLISHED)
    (ML / f"s1_grindavik_dike{TAG}.json").write_text(json.dumps(out, indent=1, ensure_ascii=False), encoding="utf-8")
    np.savez_compressed(ML / f"s1_grindavik_dike{TAG}.npz", ps=ps[:4000], pmc=pmc[RNG.choice(len(pmc), min(4000, len(pmc)), replace=False)])
    plot(G, out)


def dike_line(s):
    u = np.array([np.sin(np.deg2rad(s[2])), np.cos(np.deg2rad(s[2]))])
    return [s[0] - u[0] * s[3] / 2, s[0] + u[0] * s[3] / 2], [s[1] - u[1] * s[3] / 2, s[1] + u[1] * s[3] / 2]


def plot(G=None, out=None):
    G = G or Geometry()
    out = out or json.loads((ML / f"s1_grindavik_dike{TAG}.json").read_text(encoding="utf-8"))
    Z = np.load(ML / f"s1_grindavik_dike{TAG}.npz")
    ps = Z["ps"]
    pm = np.median(ps, 0)
    pred = G.clean(pm[None])[0]
    to_utm = Transformer.from_crs(4326, 32627, always_xy=True)
    sx, sy = to_utm.transform(*PLACES["Svartsengi"])
    pl = {k: ((to_utm.transform(*PLACES[k])[0] - sx) / 1000, (to_utm.transform(*PLACES[k])[1] - sy) / 1000)
          for k in ("Svartsengi", "Grindavik", "Sundhnukur")}
    fig = plt.figure(figsize=(18, 10.5))
    gs = fig.add_gridspec(3, 4, height_ratios=[1, 1, .62], hspace=.38, wspace=.3)
    pmm = np.median(Z["pmc"], 0)
    i0 = 0
    ttl_tr = {"A016": "висхідний трек 16, 31.10→12.11.2023", "D155": "низхідний трек 155, 10.11 (ранок)→22.11.2023"}
    for k, t in enumerate(TRACKS):
        n = G.sizes[k]
        e, nn, y, p = G.e[i0:i0 + n], G.n[i0:i0 + n], G.y[i0:i0 + n], pred[i0:i0 + n]
        A = G.design(k)
        c, *_ = np.linalg.lstsq(A, y - p, rcond=None)           # offset + plane (marginalized in the inference)
        p = p + A @ c
        lim = np.percentile(np.abs(y - np.median(y)), 99)
        for j, (v, nm) in enumerate(((y, "дані Sentinel-1"), (p, "модель за медіаною NPE"), (y - p, "залишок"))):
            a = fig.add_subplot(gs[k, j])
            sc = a.scatter(e, nn, c=v, cmap="RdBu_r", vmin=-lim + (np.median(y) if j < 2 else 0), vmax=lim + (np.median(y) if j < 2 else 0),
                           s=14, marker="s", linewidths=0)
            a.plot(*dike_line(pm), color="k", lw=1.6)
            a.plot(*dike_line(pmm), color="0.45", lw=1.4, ls="--")
            for nm2, (x0, y0) in pl.items():
                a.plot(x0, y0, "k^" if nm2 == "Svartsengi" else "ko", ms=4)
            a.set_aspect("equal"); a.set_xlim(BOX[0], BOX[1]); a.set_ylim(BOX[2], BOX[3])
            a.set_title(f"{ttl_tr[t]}\n{nm}", fontsize=9)
            a.tick_params(labelsize=7)
            plt.colorbar(sc, ax=a, fraction=0.035, pad=0.02).set_label("см", fontsize=8)
        i0 += n
    a = fig.add_subplot(gs[0:2, 3])
    for s in ps[:300]:
        a.plot(*dike_line(s), color="tab:blue", alpha=.05)
    for nm2, (x0, y0) in pl.items():
        a.plot(x0, y0, "k^" if nm2 == "Svartsengi" else "ko", ms=5)
        a.text(x0 + 0.5, y0 + (0.9 if nm2 == "Sundhnukur" else -0.9 if nm2 == "Svartsengi" else 0), {"Svartsengi": "Свартсенгі", "Grindavik": "Ґріндавік", "Sundhnukur": "Сундгнукур"}[nm2], fontsize=8)
    a.plot(*dike_line(pm), color="k", lw=1.8, label="NPE (медіана)")
    a.plot(*dike_line(pmm), color="0.45", lw=1.6, ls="--", label="MCMC (медіана)")
    a.plot(np.median(ps[:, 7]), np.median(ps[:, 8]), "r*", ms=10, label="джерело Свартсенгі (NPE)")
    a.set_aspect("equal"); a.set_xlim(-12, 12); a.set_ylim(-12, 12)
    a.legend(fontsize=7, frameon=False, loc="lower right")
    a.set_title(f"NPE: 300 можливих дайок за {out['timing']['infer_real']:.1f} с", fontsize=9)
    a.set_xlabel("схід від Свартсенгі, км", fontsize=8)
    a = fig.add_subplot(gs[2, :])
    a.axis("off")
    f = lambda v: f"{v[1]:.1f} [{v[0]:.1f}–{v[2]:.1f}]"
    N, Mc = out["npe"], out["mcmc"]
    rows = [["довжина, км", f(N["length km"]), f(Mc["length km"]), "~15"],
            ["розкриття, м", f(out["opening_m"]), f([10 ** v for v in Mc["log10 opening m"]]), "до ~8 (макс.)"],
            ["об’єм дайки, млн м³", f(out["npe_volume"]), f(out["mcmc_volume"]), "130–139"],
            ["простягання, °", f(N["strike deg"]), f(Mc["strike deg"]), "—"],
            ["Свартсенгі dV (A016), млн м³", f(N["Mogi dV A016 1e6 m3"]), f(Mc["Mogi dV A016 1e6 m3"]), "−76…−82"],
            ["Свартсенгі dV (D155), млн м³", f(N["Mogi dV D155 1e6 m3"]), f(Mc["Mogi dV D155 1e6 m3"]), "—"],
            ["лог-правдоподібність у медіані", f"{out['loglike']['npe']:.0f}", f"{out['loglike']['mcmc']:.0f}", "—"],
            ["час відповіді", f"{out['timing']['infer_real']:.1f} с", f"{out['timing']['mcmc'] / 60:.0f} хв", "дні–тижні"]]
    tb = a.table(cellText=rows, colLabels=["", "NPE", "MCMC", "опубліковано*"], loc="center", cellLoc="center", colWidths=[.25, .2, .2, .2])
    tb.auto_set_font_size(False); tb.set_fontsize(9); tb.scale(1, 1.45)
    a.text(0.5, -0.32, "медіана [16–84 %];  * Sigmundsson et al. 2024 (Science 383:1228), EGU24-14829; dV Свартсенгі опубліковано до 12.11", fontsize=9, ha="center", va="top", transform=a.transAxes)
    fig.suptitle("Справжня дайка 10.11.2023 під Ґріндавіком: Sentinel-1 (ASF HyP3) → NPE за секунди; Свартсенгі — "
                 + ("шарувата кора + рельєф (МСЕ)" if out.get("source", "mogi") == "fem" else "Моґі, однорідний півпростір"), fontsize=13)
    fig.savefig(ML / f"s1_grindavik_dike{TAG}.png", dpi=105, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    print("wrote", ML / f"s1_grindavik_dike{TAG}.png")


if __name__ == "__main__":
    steps = [a for a in sys.argv[1:] if a in ("download", "prep", "npe", "plot")] or ["download", "prep", "npe"]
    for s in steps:
        globals()[s]()
