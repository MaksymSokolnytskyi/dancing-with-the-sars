"""Source inference with the expensive physics of fem_layered.py (layered crust + topography),
compared with the homogeneous Mogi results.

  table   Green's-function table of the FE model on a depth grid (1-10 km, 0.1 km): 91 FE
          solves, run once. With it the forward model costs milliseconds again (an emulator),
          which is the only way MCMC is affordable here: MCMC straight on the FE model would
          need ~256 000 solves x ~5 s = ~15 days.
  mcmc    MCMC (same sampler, priors and likelihood as ml_mcmc_source.py --noise screens) on
          the emulator, for: layered + topography, layered without topography; NISAR + GNSS.
  npe     NPE trained on DIRECT FE simulations (one solve per simulation, no emulator), with the
          whitened/PCA flow of ml_npe_flow.py and a truncated second round; compared with the
          emulator MCMC. This is the expensive-physics case NPE is built for.

Topography: varying-depth approximation; a point at elevation h sees the source at depth
D + h (D below sea level). For the small shifts here (h < 0.4 km) the depth-D field is
mapped self-similarly: G(r, D + h) ~ (D / (D + h))^2 G(r D / (D + h), D).

Usage: python scripts/reykjanes/ml_fem_inference.py table|mcmc|npe [--sims1 1500 --sims2 3000 --workers 6]
"""
import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import rasterio

import reykjanes_mogi as M
import ml_npe_source as NS
from build_breathing_app import local_to_lonlat
from fem_layered import DEM, green
from ml_mcmc_source import Likelihood, autocorr_time, rhat, stretch_sampler
from reykjanes_joint import SIGMA_GNSS

ML = Path(__file__).resolve().parents[2] / "output" / "reykjanes" / "ml"
TABLE = ML / "fem_green_table.npz"
DEPTHS = np.round(np.arange(1.0, 10.01, 0.1), 2)
R_GRID = np.r_[np.linspace(0, 10, 101), np.linspace(10.25, 40, 120)]
LO = np.array([-6, -6, 1.2, np.log10(2e-3)])
HI = np.array([6, 6, 9.5, np.log10(0.08)])
NAMES = ("east km", "north km", "depth km b.s.l.", "dV 1e6 m3/yr")


def arg(name, default):
    return type(default)(sys.argv[sys.argv.index(name) + 1]) if name in sys.argv else default


def solve_depth(d):
    r, ur, uz = green(float(d), layered=True)
    return np.interp(R_GRID, r, ur), np.interp(R_GRID, r, uz)


def build_table(workers):
    t0 = time.time()
    from progress import Progress
    res = []
    with Progress("таблиця Гріна МСЕ (91 глибина)", total=len(DEPTHS), unit="розв'язок", note=f"{workers} процесів") as p,             ProcessPoolExecutor(workers) as ex:
        for r in ex.map(solve_depth, DEPTHS):
            res.append(r)
            p.update(1)
    np.savez(TABLE, depths=DEPTHS, r=R_GRID, ur=np.array([a for a, _ in res]), uz=np.array([b for _, b in res]))
    print(f"{len(DEPTHS)} FE solves in {time.time() - t0:.0f} s with {workers} workers -> {TABLE.name}")


class Geometry:
    """Observation points (InSAR bins + GNSS + reference), their LOS vectors and elevations."""
    def __init__(self, use_topo=True):
        O, _ = NS.cached(NS.Obs)
        self.O = O
        self.e = np.concatenate([O.tr[t]["e"] for t in M.TRACKS] + [O.gn[:, 0], O.gref[:1]])
        self.n = np.concatenate([O.tr[t]["n"] for t in M.TRACKS] + [O.gn[:, 1], O.gref[1:]])
        self.sizes = [O.tr[t]["y"].size for t in M.TRACKS]
        self.vec = np.concatenate([O.tr[t]["vec"] for t in M.TRACKS])
        lon, lat = local_to_lonlat(self.e, self.n)
        with rasterio.open(DEM) as src:
            h = np.array([v[0] for v in src.sample(zip(lon, lat))], float)
        h[(h < -100) | ~np.isfinite(h)] = 0.0
        self.h = np.clip(h, 0, None) / 1000 if use_topo else np.zeros_like(h)
        self.n_ins = sum(self.sizes)

    def enu(self, q, G):
        """q = (e0, n0, D, log10 dV km^3); G = (ur(r), uz(r)) per m^3 for source depth D (no topo shift)."""
        e0, n0, D, lv = q
        dV = 10 ** lv * 1e9
        s = D / (D + self.h)
        de, dn = self.e - e0, self.n - n0
        r = np.hypot(de, dn)
        ur = s ** 2 * np.interp(r * s, R_GRID, G[0]) * dV * 100      # cm/yr
        uz = s ** 2 * np.interp(r * s, R_GRID, G[1]) * dV * 100
        rr = np.maximum(r, 1e-6)
        return np.stack([ur * de / rr, ur * dn / rr, uz], -1)

    def predict(self, q, G):
        u = self.enu(q, G)
        los = (u[:self.n_ins] * self.vec).sum(1)
        g = (u[self.n_ins:-1] - u[-1]).ravel()
        return np.r_[los, g]


class Emulator:
    def __init__(self):
        z = np.load(TABLE)
        self.d, self.ur, self.uz = z["depths"], z["ur"], z["uz"]

    def G(self, D):
        i = np.clip(np.searchsorted(self.d, D) - 1, 0, len(self.d) - 2)
        w = (D - self.d[i]) / (self.d[i + 1] - self.d[i])
        return (1 - w) * self.ur[i] + w * self.ur[i + 1], (1 - w) * self.uz[i] + w * self.uz[i + 1]


class EmuModel:
    """Model interface for ml_mcmc_source.stretch_sampler."""
    lo, hi = LO, HI

    def __init__(self, geo, emu):
        self.geo, self.emu = geo, emu

    def predict(self, q):
        return np.array([self.geo.predict(qi, self.emu.G(qi[2])) for qi in q])

    def wrap(self, q):
        return q


def likelihood(geo, use_gnss=True):
    z = np.load(ML / "screens.npz")
    return Likelihood(geo.O, {t: z[f"cov_par_{t}"] for t in M.TRACKS}, use_gnss)


def run_mcmc(steps=3000, walkers=64):
    from ml_mcmc_source import log_post
    emu = Emulator()
    rng = np.random.default_rng(5)
    out = {}
    for label, topo in (("layered + topography", True), ("layered, flat", False)):
        geo = Geometry(topo)
        model = EmuModel(geo, emu)
        like = likelihood(geo)
        cand = LO + rng.random((3000, 4)) * (HI - LO)
        lp = log_post(model, like, cand)
        t0 = time.time()
        chain, lps, acc = stretch_sampler(model, like, cand[np.argsort(lp)[-walkers:]], steps, walkers, rng, label=label)
        post = chain[steps // 3:]
        mlp = lps[steps // 3:].mean(0)
        post = post[:, mlp > np.median(mlp) - 10]
        flat = post.reshape(-1, 4).copy()
        flat[:, 3] = 10 ** flat[:, 3] * 1e3
        summ = {nm: np.round(np.percentile(flat[:, i], [2.5, 16, 50, 84, 97.5]), 3).tolist() for i, nm in enumerate(NAMES)}
        out[label] = dict(summary=summ, rhat=rhat(post).tolist(), seconds=round(time.time() - t0), n_forward=steps * walkers,
                          samples=flat[np.random.default_rng(0).choice(len(flat), 6000, replace=False)])
        print(f"[MCMC emulator, {label}] {time.time() - t0:.0f} s, R-hat {np.round(rhat(post), 3).tolist()}")
        for nm, v in summ.items():
            print(f"   {nm:17s} {v[2]:7.2f} [{v[1]:.2f}..{v[3]:.2f}]")
    np.savez(ML / "fem_mcmc.npz", **{k.replace(" ", "_").replace(",", "").replace("+", "and"): v["samples"] for k, v in out.items()})
    (ML / "fem_mcmc.json").write_text(json.dumps({k: {kk: vv for kk, vv in v.items() if kk != "samples"} for k, v in out.items()}, indent=1),
                                      encoding="utf-8")


def _solve(d):
    """One DIRECT FE solve (worker process) at source depth d: surface Green's functions on R_GRID."""
    r, ur, uz = green(float(d), layered=True)
    return np.interp(R_GRID, r, ur), np.interp(R_GRID, r, uz)


def simulate_depths(depths, workers):
    """A full FE solve for each sampled depth (no interpolation between depths)."""
    t0 = time.time()
    from progress import Progress
    out = []
    with Progress(f"прямі розв'язки МСЕ ({len(depths)})", total=len(depths), unit="розв'язок", note=f"{workers} процесів") as p,             ProcessPoolExecutor(workers) as ex:
        for g in ex.map(_solve, depths, chunksize=4):
            out.append(g)
            p.update(1)
    return np.array(out), time.time() - t0          # (n, 2, len(R_GRID))


def expand(geo, depths, G, lo, hi, per, rng):
    """Each FE solve fixes the field shape for its depth; position and volume only translate and scale it,
    so one solve gives `per` exact training simulations (east, north, dV drawn from the box)."""
    q, X = [], []
    for d, g in zip(depths, G):
        qq = lo + rng.random((per, 4)) * (hi - lo)
        qq[:, 2] = d
        q.append(qq)
        X.append(np.array([geo.predict(qi, g) for qi in qq]))
    return np.concatenate(q), np.concatenate(X)


def run_npe(n1, n2, workers):
    import torch
    from ml_npe_flow import from_z, posterior, to_z, train
    K = arg("--k", 64)
    per = arg("--per", 60)
    epochs = arg("--epochs", 80)
    rng = np.random.default_rng(9)
    geo = Geometry(True)
    like = likelihood(geo)
    cache = ML / "fem_direct_solves.npz"
    old = dict(np.load(cache)) if cache.exists() else {"d": np.zeros(0), "G": np.zeros((0, 2, len(R_GRID)))}

    def whiten(x):
        out, i0 = [], 0
        for n, W in zip(geo.sizes, like.W):
            out.append(x[:, i0:i0 + n] @ W.T)
            i0 += n
        out.append(x[:, i0:] / SIGMA_GNSS)
        return np.concatenate(out, 1)

    def solves(depths, tag):
        G, t = simulate_depths(depths, workers)
        old["d"], old["G"] = np.r_[old["d"], depths], np.concatenate([old["G"], G])
        np.savez(cache, **old)
        T[tag] = t
        return G

    T = {}
    d1 = LO[2] + rng.random(n1) * (HI[2] - LO[2])
    G1 = solves(d1, "sim1")
    print(f"round 1: {n1} direct FE solves in {T['sim1']:.0f} s ({T['sim1'] / n1 * workers:.1f} s per solve per worker), "
          f"x{per} positions/volumes each", flush=True)
    q1, X1 = expand(geo, d1, G1, LO, HI, per, rng)
    W1 = whiten(X1)
    _, sv, Vt = np.linalg.svd(W1[:20000] - W1[:20000].mean(0), full_matrices=False)
    U = Vt[:K].T
    print(f"   {K} components; component {K} still has {sv[K - 1] / np.sqrt(min(len(W1), 20000)):.2f} x the noise", flush=True)
    xobs = whiten(np.concatenate([np.concatenate([geo.O.tr[t]["y"] for t in M.TRACKS]), geo.O.gv])[None]) @ U
    net1, T["train1"] = train(W1 @ U, to_z(q1, LO, HI), None, epochs, "FE round 1", 1)
    s1 = posterior([net1], xobs, LO, HI)
    a, b = np.quantile(s1, [0.0005, 0.9995], 0)
    lo2, hi2 = np.maximum(LO, a - 0.25 * (b - a)), np.minimum(HI, b + 0.25 * (b - a))
    print("   truncated box: " + ", ".join(f"{nm.split()[0]} {l:.2f}..{h:.2f}" for nm, l, h in zip(NAMES, lo2, hi2)), flush=True)
    d2 = lo2[2] + rng.random(n2) * (hi2[2] - lo2[2])
    G2 = solves(d2, "sim2")
    print(f"round 2: {n2} direct FE solves in {T['sim2']:.0f} s", flush=True)
    q2, X2 = expand(geo, d2, G2, lo2, hi2, per, rng)
    nets, T["train2"] = [], 0.0
    for k in range(arg("--ens", 1)):
        net, tt = train(whiten(X2) @ U, to_z(q2, lo2, hi2), None, epochs, f"FE round 2 #{k + 1}", 20 + k)
        nets.append(net)
        T["train2"] += tt
    t0 = time.time()
    s2 = posterior(nets, xobs, lo2, hi2)
    T["infer"] = time.time() - t0
    # calibration on held-out direct FE solves from the truncated box (one position/volume each)
    ncal = arg("--ncal", 150)
    dt = lo2[2] + rng.random(ncal) * (hi2[2] - lo2[2])
    Gt, _ = simulate_depths(dt, workers)
    qt, Xt = expand(geo, dt, Gt, lo2, hi2, 1, rng)
    xt = whiten(Xt) @ U + rng.standard_normal((ncal, K))
    st = np.concatenate([np.concatenate([from_z(net.sample(torch.tensor(xt[j:j + 25], dtype=torch.float32), 2000 // len(nets)).numpy().astype(float), lo2, hi2)
                                         for net in nets]) for j in range(0, ncal, 25)], axis=1)
    cal = {lev: ((qt >= np.quantile(st, (1 - lev) / 2, 0)) & (qt <= np.quantile(st, (1 + lev) / 2, 0))).mean(0).round(3).tolist()
           for lev in (.68, .95)}
    print(f"   calibration on held-out FE solves, 68 %: {cal[.68]}; 95 %: {cal[.95]}")
    s2[:, 3] = 10 ** s2[:, 3] * 1e3
    summ = {nm: np.round(np.percentile(s2[:, i], [2.5, 16, 50, 84, 97.5]), 3).tolist() for i, nm in enumerate(NAMES)}
    mc = json.loads((ML / "fem_mcmc.json").read_text(encoding="utf-8")).get("layered + topography") if (ML / "fem_mcmc.json").exists() else None
    print("[NPE on direct FE simulations] vs [MCMC on the emulator]:")
    cmp = {}
    for nm, v in summ.items():
        line = f"   {nm:17s} NPE {v[2]:7.2f} [{v[1]:.2f}..{v[3]:.2f}]"
        if mc:
            m = mc["summary"][nm]
            sd = (m[3] - m[1]) / 2
            cmp[nm] = dict(shift_sd=round((v[2] - m[2]) / sd, 2), width_ratio=round((v[3] - v[1]) / (m[3] - m[1]), 2))
            line += f" | MCMC {m[2]:7.2f} [{m[1]:.2f}..{m[3]:.2f}] | shift {cmp[nm]['shift_sd']:+.2f} sd, width x{cmp[nm]['width_ratio']:.2f}"
        print(line)
    T = {k: round(v, 1) for k, v in T.items()}
    print(f"   time: {n1 + n2} FE solves {T['sim1'] + T['sim2']:.0f} s on {workers} workers ({(n1 + n2) * per} training simulations), "
          f"training {T['train1'] + T['train2']:.0f} s, posterior for new data {T['infer']:.2f} s")
    (ML / "fem_npe.json").write_text(json.dumps(dict(summary=summ, compare=cmp, calibration=cal, timing=T, n_solves=n1 + n2, per_solve=per,
                                                     k_pca=K, workers=workers, box=[lo2.tolist(), hi2.tolist()]), indent=1), encoding="utf-8")
    np.savez(ML / "fem_npe.npz", samples=s2[:6000])


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "table":
        build_table(arg("--workers", 6))
    elif cmd == "mcmc":
        run_mcmc()
    else:
        run_npe(arg("--sims1", 1500), arg("--sims2", 3000), arg("--workers", 6))
