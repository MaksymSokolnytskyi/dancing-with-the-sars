"""NPE v4: richer dike kinematics + a BLIND set-up (knowledge frozen before 10 Nov 2023).

Physics as in Sigmundsson et al. 2024 (homogeneous elastic half-space): Okada dislocations for the dike and the graben,
Mogi point source for Svartsengi. Kinematics (20 parameters):
  plane     centre east/north (km from Svartsengi), strike, length, top depth, height, dip (70-110 deg),
            two straight segments with a bend at the centre (+-30 deg)
  opening   2-D Gaussian patch on the plane, 6 x 4 patches per segment: peak (log10 m, up to 20 m), centre and width
            along strike and along depth (fractions of length / height)
  shear     uniform strike-slip on the dike plane (+-2 m, positive = left-lateral)
  graben    two normal faults from the surface down to the top of the dike, dipping 65 deg towards the axis, slip 0-2 m
  source    Svartsengi (Mogi): east, north, depth, volume change per window (A016, D155)
BLIND: the prior covers the whole study area (centre +-15 / +-12 km, any strike, length 2-30 km); no mask built from the
result (all 1 km cells used: the graben is MODELLED instead); the noise level is unknown (hierarchical noise as v3);
25 % of simulations without a dike. Training uses simulations only; nothing measured after 9 Nov 2023 enters the design.

Tests: SBC + posterior contraction per parameter (which kinematic parameters the data constrain); recovery of a
synthetic "paper-like" dike (8 m at 2-4 km between Grindavik and Sundhnukur) under Nov-2023-level noise; the real pair;
null pairs and 2024 events.
Writes output/reykjanes/ml/s1_npe_v4.json / .npz / .pt and s1_npe_v4_*.png.
Usage: python scripts/reykjanes/s1_dike_npe_v4.py [--sims 100000] [--epochs 25] [--workers 3] | --post
       --v5: robust version - 160 000 simulations x 3 noise realizations, single-threaded workers, and NUISANCES not in the
       physics added to half of the examples (secondary small source anywhere, stratified troposphere, unwrapping-error
       patches, local lava compaction, GNSS steps); SBC on 500 cases split by nuisance yes / no. Outputs s1_npe_v5.*
"""
import os
import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime

if "--source" not in sys.argv:
    sys.argv += ["--source", "mogi"]
import numpy as np
import torch

import s1_grindavik_dike as S
from ml_npe_flow import Net, from_z, to_z, train
from okada import dislocation, dislocation_fast
from progress import Progress
from s1_common_v2 import GeomV2, HierNoise, gnss_vector, pair_vector, product_dir

ML = S.ML
NAMES = ("dike east km", "dike north km", "strike deg", "length km", "top km", "height km", "dip deg",
         "log10 peak opening m", "opening centre along", "opening width along", "opening centre depth", "opening width depth",
         "strike-slip m", "graben slip m", "bend deg",
         "src east km", "src north km", "src depth km", "src dV A016 1e6 m3", "src dV D155 1e6 m3")
LO = np.array([-15, -12, 0, 2, 0, 1, 70, -1.5, -0.5, 0.1, 0.0, 0.1, -2, 0, -30, -6, -6, 2, -200, -200], float)
HI = np.array([15, 12, 180, 30, 3, 8, 110, 1.3, 0.5, 0.8, 1.0, 0.8, 2, 2, 30, 6, 6, 8, 30, 30], float)
NA, ND = (5, 3) if "--v5" in sys.argv else (6, 4)      # patches per segment along strike / along depth
NODIKE = 0.25
K = 256
V5 = "--v5" in sys.argv
TAG = "_v5" if V5 else "_v4"
RNG = np.random.default_rng(80 if V5 else 70)
TRACE0 = (0.0, 0.0, 0.0, 1.0)  # dummy trace: no mask is built from any result
_G = None


def geom():
    return GeomV2(gnss=True, mask_km=0, noise="ref", trace=TRACE0)


def patches(p):
    """Centres (east, north), depths, opening, strike and dip for both segments (NA x ND each)."""
    e0, n0, st, L, top, H, dip = p[:7]
    peak, xc, sx, zc, sz = 10 ** p[7], p[8], p[9], p[10], p[11]
    out = []
    for side, s_strike in ((+1, st + p[14] / 2), (-1, st - p[14] / 2)):
        u = np.array([np.sin(np.deg2rad(s_strike)), np.cos(np.deg2rad(s_strike))])
        r = np.array([np.cos(np.deg2rad(st)), -np.sin(np.deg2rad(st))])            # right of the mean strike
        a = (np.arange(NA) + 0.5) * (L / 2) / NA                                  # along-strike distance from the centre
        z = top + (np.arange(ND) + 0.5) * H / ND                                  # patch-centre depths
        A, Z = np.meshgrid(a, z)
        off = (Z - top) / np.tan(np.deg2rad(dip))                                 # dip < 90: plane goes to the right
        ce = e0 + side * A * u[0] + off * r[0]
        cn = n0 + side * A * u[1] + off * r[1]
        xn, zn = side * A / L, (Z - top) / H
        op = peak * np.exp(-((xn - xc) ** 2 / (2 * sx ** 2) + (zn - zc) ** 2 / (2 * sz ** 2)))
        out.append((ce.ravel(), cn.ravel(), Z.ravel(), op.ravel(), s_strike))
    return out


def forward(G, p):
    e0, n0, st, L, top, H, dip = p[:7]
    sd_eff, sgn = (dip, 1.0) if dip <= 90 else (180 - dip, -1.0)
    u = np.zeros((len(G.e), 3))
    E, N = G.e[:, None], G.n[:, None]
    for ce, cn, z, op, s_strike in patches(p):
        strike_eff = s_strike if dip <= 90 else s_strike + 180
        w = (H / ND) / np.sin(np.deg2rad(sd_eff))
        u += dislocation_fast(E, N, ce[None], cn[None], z[None], (L / 2) / NA, w, strike_eff, sd_eff,
                         strike_slip=sgn * p[12] * np.ones((1, len(op))), opening=op[None]).sum(1)
    if p[13] > 0 and top > 0.2:                                                   # graben: two normal faults
        r = np.array([np.cos(np.deg2rad(st)), -np.sin(np.deg2rad(st))])
        wsurf = top / np.tan(np.deg2rad(65))
        Wf = top / np.sin(np.deg2rad(65))
        for side, fstrike in ((-1, st), (+1, st + 180)):                          # left fault dips right, right fault dips left
            c = np.array([e0, n0]) + side * r * wsurf / 2
            u += dislocation(G.e, G.n, c[0], c[1], top / 2, L, Wf, fstrike, 65, dip_slip=-p[13])
    u = u * 100 + G.source_enu(np.r_[np.zeros(7), p[15:20]])
    return (u * G.vec).sum(1)


def volume(s):
    s = np.atleast_2d(s)
    out = []
    for p in s:
        dip = p[6] if p[6] <= 90 else 180 - p[6]
        area = ((p[3] / 2) / NA) * ((p[5] / ND) / np.sin(np.deg2rad(dip)))
        out.append(sum(op.sum() for *_, op, _ in patches(p)) * area)
    return np.array(out)


def peak_depth(p):
    return p[4] + p[10] * p[5]


def draw(n, lo=LO, hi=HI):
    q = lo + RNG.random((n, len(lo))) * (hi - lo)
    nd = RNG.random(n) < NODIKE
    q[nd, 7] = RNG.uniform(-1.5, -1.2, nd.sum())
    q[nd, 12] = RNG.uniform(-0.05, 0.05, nd.sum())
    q[nd, 13] = 0.0
    return q, nd


def _init():
    global _G
    torch.set_num_threads(1)
    _G = geom()


LAMBDA_HALF = 2.77         # cm: Sentinel-1 C-band half wavelength = one unwrapping cycle


def nuisance(G, clean, rng, frac=0.5):
    """Deformation / artefacts NOT in the physics, added to a fraction of the examples (data space, cm)."""
    N, n = clean.shape
    out = np.zeros_like(clean)
    use = rng.random(N) < frac
    idx = np.flatnonzero(use)
    if not len(idx):
        return out, use
    E, Nn, vec = G.e[None], G.n[None], G.vec
    # 1) secondary small Mogi source anywhere (geothermal / other magma): +-5e6 m3, depth 1-6 km
    k = idx[rng.random(len(idx)) < 0.5]
    if len(k):
        e0, n0 = rng.uniform(-15, 15, (len(k), 1)), rng.uniform(-12, 12, (len(k), 1))
        d = rng.uniform(1, 6, (len(k), 1)); dv = rng.uniform(-5, 5, (len(k), 1)) * 1e-3
        de, dn = E - e0, Nn - n0
        R3 = (de ** 2 + dn ** 2 + d ** 2) ** 1.5
        c = 3 / (4 * np.pi) * 1e5 * dv
        out[k] += (c * de / R3) * vec[:, 0] + (c * dn / R3) * vec[:, 1] + (c * d / R3) * vec[:, 2]
    ni = G.n_ins
    for t in range(len(S.TRACKS)):
        m = np.flatnonzero(G.trk[:ni] == t)
        # 2) stratified troposphere: delay proportional to elevation, +-1.5 cm/km (1 sigma)
        k = idx[rng.random(len(idx)) < 0.5]
        out[np.ix_(k, m)] += rng.normal(0, 1.5, (len(k), 1)) * G.h[m][None]
        # 3) unwrapping errors: 1-2 discs of +-1 or 2 cycles
        k = idx[rng.random(len(idx)) < 0.3]
        for i in k:
            for _ in range(rng.integers(1, 3)):
                c0 = G.e[m][rng.integers(len(m))], G.n[m][rng.integers(len(m))]
                r = rng.uniform(1, 3)
                disc = np.hypot(G.e[m] - c0[0], G.n[m] - c0[1]) < r
                out[i, m[disc]] += rng.choice([-2, -1, 1, 2]) * LAMBDA_HALF
    # 4) local subsidence of fresh lava / compaction: 1-5 cm, sigma 0.5-2 km, anywhere within 10 km
    k = idx[rng.random(len(idx)) < 0.3]
    if len(k):
        e0, n0 = rng.uniform(-10, 10, (len(k), 1)), rng.uniform(-10, 10, (len(k), 1))
        sg, a = rng.uniform(0.5, 2, (len(k), 1)), rng.uniform(1, 5, (len(k), 1))
        g = -a * np.exp(-((E[:, :ni] - e0) ** 2 + (Nn[:, :ni] - n0) ** 2) / (2 * sg ** 2))
        out[np.ix_(k, np.arange(ni))] += g * vec[None, :ni, 2]
    # 5) GNSS step in one component (like VOGS on 20 Nov 2023)
    if G.n_gnss:
        k = idx[rng.random(len(idx)) < 0.1]
        for i in k:
            out[i, ni + rng.integers(G.n_gnss)] += rng.choice([-1, 1]) * rng.uniform(3, 30)
    return out, use


def _chunk(q):
    return np.array([forward(_G, p) for p in q], dtype=np.float32)


def simulate(q, workers, label):
    for v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
        os.environ[v] = "1"                                          # one thread per worker: no oversubscription
    out = []
    chunks = np.array_split(q, max(1, len(q) // 250))
    with Progress(f"NPE v4: симуляції, {label} ({len(q)})", total=len(q), unit="сим", note=f"{workers} процесів") as pr, \
            ProcessPoolExecutor(workers, initializer=_init) as ex:
        for r in ex.map(_chunk, chunks):
            out.append(r)
            pr.update(len(r))
    return np.concatenate(out)


class Feat:
    def __init__(self, Xw, ni):
        _, sv, Vt = np.linalg.svd(Xw[:, :ni] - Xw[:, :ni].mean(0), full_matrices=False)
        self.U, self.ni = Vt[:K].T.astype(np.float32), ni
        self.info = f"{K} components keep {(sv[:K] ** 2).sum() / (sv ** 2).sum():.4%}"
        self.scale = None

    def __call__(self, Xw):
        return np.c_[Xw[:, :self.ni] @ self.U, Xw[:, self.ni:]].astype(np.float32)


def sample(net, F, Xw, n):
    X = torch.tensor(F(Xw) / F.scale, dtype=torch.float32)
    with torch.no_grad():
        return np.concatenate([from_z(net.flow(net.emb(X[j:j + 10])).sample((n,)).numpy().astype(np.float32), LO, HI)
                               for j in range(0, len(X), 10)], axis=1)


def summarize(ps):
    V = volume(ps[:3000])
    return dict(pct={nm: np.percentile(ps[:, i], [16, 50, 84]).round(3).tolist() for i, nm in enumerate(NAMES)},
                volume=np.percentile(V, [2.5, 16, 50, 84, 97.5]).round(1).tolist(), p_nodike=round(float(np.mean(V < 5)), 3),
                peak_opening=np.percentile(10 ** ps[:, 7], [16, 50, 84]).round(2).tolist(),
                peak_depth=np.percentile([peak_depth(p) for p in ps[:3000]], [16, 50, 84]).round(2).tolist())


def main():
    workers = int(S.arg("--workers", 4 if V5 else 3))
    n_sims, epochs, reps = int(S.arg("--sims", 120000 if V5 else 100000)), int(S.arg("--epochs", 20 if V5 else 25)), 3 if V5 else 2
    G = geom()
    HN = HierNoise(G)
    T, res = {}, {}
    if "--post" in sys.argv:
        ck = torch.load(ML / f"s1_npe{TAG}.pt", weights_only=False)
        F = Feat.__new__(Feat); F.U, F.ni, F.scale = ck["U"], ck["ni"], ck["scale"]
        net = Net(F.U.shape[1] + G.n_gnss, len(LO)); net.load_state_dict(ck["state"]); net.eval()
        T = json.loads((ML / f"s1_npe{TAG}_timing.json").read_text())
    else:
        t0 = time.time()
        q, nd = draw(n_sims)
        clean = simulate(q, workers, "навчання")
        T["simulate"] = round(time.time() - t0)
        X, F = [], None
        for rep in range(reps):
            for j in range(0, len(clean), 5000):
                c = clean[j:j + 5000]
                nu = nuisance(G, c, RNG)[0] if V5 else 0.0
                xw = G.whiten(c + nu + HN(c, RNG)[0]).astype(np.float32)
                if F is None:
                    F = Feat(xw, G.n_ins)
                X.append(F(xw))
        del clean
        X = np.concatenate(X)
        F.scale = X.std(0) + 1.0
        print(F.info, flush=True)
        Q = np.concatenate([to_z(q, LO, HI)] * reps)
        net, tt = train(X, Q, None, epochs, f"NPE {TAG[1:]} (кінематика, сліпий)", 81 if V5 else 71, scale=F.scale)
        T["train"] = round(tt)
        del X, Q
        torch.save(dict(state=net.state_dict(), U=F.U, scale=F.scale, ni=F.ni), ML / f"s1_npe{TAG}.pt")
        (ML / f"s1_npe{TAG}_timing.json").write_text(json.dumps(T))
    # ---- real 10 Nov 2023 pair
    t0 = time.time()
    ps = sample(net, F, G.whiten(G.y[None]), 8000)[:, 0].astype(float)
    T["infer_real"] = round(time.time() - t0, 1)
    res["real_20231110"] = summarize(ps)
    np.savez_compressed(ML / f"s1_npe{TAG}.npz", ps=ps)
    r = res["real_20231110"]
    print(f"REAL {TAG[1:]}: V {r['volume']}; peak {r['peak_opening']} m at {r['peak_depth']} km; strike {r['pct']['strike deg']}; dip {r['pct']['dip deg']}; "
          f"strike-slip {r['pct']['strike-slip m']}; graben {r['pct']['graben slip m']}; src depth {r['pct']['src depth km']}, dV_A {r['pct']['src dV A016 1e6 m3']}", flush=True)
    # ---- SBC + contraction (hierarchical noise)
    nt = 500 if V5 else 300
    qt, ndt = draw(nt)
    ct = simulate(qt, workers, "SBC")
    nu, used = nuisance(G, ct, RNG) if V5 else (0.0, np.zeros(nt, bool))
    st = sample(net, F, G.whiten(ct + nu + HN(ct, RNG)[0]), 800)
    ranks = (st < qt[None]).sum(0)
    from scipy.stats import chisquare
    unif = lambda x: float(chisquare(np.histogram(x, bins=16, range=(0, 800))[0]).pvalue)
    dk = ~ndt
    contraction = {nm: round(float(1 - np.median(st[:, dk, i].var(0)) / ((HI[i] - LO[i]) ** 2 / 12)), 3) for i, nm in enumerate(NAMES)}
    res["sbc"] = dict(n=nt, p_uniform={nm: round(unif(ranks[:, i]), 3) for i, nm in enumerate(NAMES)}, contraction=contraction, ranks=ranks.tolist())
    if V5:
        cov68 = lambda m: {nm: round(float(np.mean((qt[m, i] >= np.quantile(st[:, m, i], .16, 0)) & (qt[m, i] <= np.quantile(st[:, m, i], .84, 0)))), 3) for i, nm in enumerate(NAMES)}
        Vt, Vs = volume(qt), np.array([volume(st[:200, j].astype(float)) for j in range(nt)]).T
        vcov = lambda m: round(float(np.mean((Vt[m] >= np.quantile(Vs[:, m], .16, 0)) & (Vt[m] <= np.quantile(Vs[:, m], .84, 0)))), 3)
        pnd = np.mean(Vs < 5, 0)
        res["sbc"]["coverage68_with_nuisance"] = cov68(used & dk); res["sbc"]["coverage68_clean"] = cov68(~used & dk)
        res["sbc"]["volume_coverage68"] = dict(with_nuisance=vcov(used & dk), clean=vcov(~used & dk))
        res["sbc"]["p_nodike_median"] = dict(nodike_with_nuisance=round(float(np.median(pnd[ndt & used])), 3), nodike_clean=round(float(np.median(pnd[ndt & ~used])), 3),
                                             dikes_gt20=round(float(np.median(pnd[dk & (Vt > 20)])), 3))
        print("v5 coverage68 (volume): with nuisance", res["sbc"]["volume_coverage68"], "; P(no dike):", res["sbc"]["p_nodike_median"], flush=True)
    print("SBC v4 p:", res["sbc"]["p_uniform"], flush=True)
    print("contraction (1 - posterior var / prior var), dikes only:", contraction, flush=True)
    # ---- synthetic paper-like dike, Nov-2023-level noise, 20 realizations
    truth = np.array([1.5, -1.0, 35, 15, 1.0, 4.5, 90, np.log10(8.0), 0.1, 0.2, 0.55, 0.25, 0.5, 0.8, 0.0, -2.5, -1.0, 4.1, -80, -80])
    c = np.repeat(forward(G, truth)[None].astype(np.float32), 20, 0)
    Lc = HN.L
    noise = np.zeros_like(c)
    i0 = 0
    for k in range(len(S.TRACKS)):
        n = G.sizes[k]
        noise[:, i0:i0 + n] = RNG.standard_normal((20, n)).astype(np.float32) @ Lc[k][1].T * np.sqrt(30.0)     # range 5 km, sill 30 cm2
        i0 += n
    noise[:, G.n_ins:] = RNG.standard_normal((20, G.n_gnss)) * np.sqrt(G.gnss_sigma_base ** 2 + (0.05 * c[:, G.n_ins:]) ** 2)
    if V5:
        noise = noise + nuisance(G, c, RNG, frac=1.0)[0]
    sy = sample(net, F, G.whiten(c + noise), 1500)
    Vtrue = float(volume(truth)[0])
    rec = [dict(V=np.percentile(volume(sy[:1000, j]), [16, 50, 84]).round(1).tolist(), peak=np.percentile(10 ** sy[:, j, 7], [16, 50, 84]).round(2).tolist(),
                peak_depth=np.percentile([peak_depth(p) for p in sy[:1000, j]], [16, 50, 84]).round(2).tolist(),
                strike=np.percentile(sy[:, j, 2], [16, 50, 84]).round(1).tolist()) for j in range(20)]
    cover = lambda key, tv: float(np.mean([a[key][0] <= tv <= a[key][2] for a in rec]))
    res["synthetic_paper_like"] = dict(truth=dict(zip(NAMES, truth.round(3).tolist())), V_true=round(Vtrue, 1), peak_true=8.0, peak_depth_true=round(peak_depth(truth), 2),
                                       noise="exponential sill 30 cm2, range 5 km (Nov 2023 level) + offset/plane absent",
                                       median_V=round(float(np.median([a["V"][1] for a in rec])), 1),
                                       median_peak=round(float(np.median([a["peak"][1] for a in rec])), 2),
                                       median_peak_depth=round(float(np.median([a["peak_depth"][1] for a in rec])), 2),
                                       cover68=dict(V=cover("V", Vtrue), peak=cover("peak", 8.0), peak_depth=cover("peak_depth", peak_depth(truth)), strike=cover("strike", 35.0)),
                                       realizations=rec)
    print("synthetic paper-like:", {k: v for k, v in res["synthetic_paper_like"].items() if k != "realizations"}, flush=True)
    # ---- null pairs and 2024 events (no retraining)
    idx = json.loads((S.S1 / "pairs_index.json").read_text())
    post = {"real_20231110": ps[:2000]}
    for name, pairs in idx.items():
        dirs = {p["track"]: product_dir(p["track"], p["ref"], p["sec"]) for p in pairs}
        yi, _ = pair_vector(G, dirs)
        dates = {p["track"]: (datetime.strptime(p["ref"][17:25], "%Y%m%d").date().isoformat(), datetime.strptime(p["sec"][17:25], "%Y%m%d").date().isoformat()) for p in pairs}
        y = np.r_[yi, gnss_vector(G, dates)]
        xw = G.whiten(y[None]); xw[~np.isfinite(xw)] = 0
        pe = sample(net, F, xw, 6000)[:, 0].astype(float)
        res[name] = summarize(pe)
        post[name] = pe[:2000]
        print(f"{name:16s} P(no dike) {res[name]['p_nodike']:.2f}; V {res[name]['volume'][1:4]}; centre ({res[name]['pct']['dike east km'][1]:.1f}, "
              f"{res[name]['pct']['dike north km'][1]:.1f}); strike {res[name]['pct']['strike deg'][1]:.0f}; src dV {res[name]['pct']['src dV A016 1e6 m3'][1]:.0f}", flush=True)
    np.savez_compressed(ML / f"s1_npe{TAG}_events.npz", **post)
    res["timing"], res["prior"], res["names"] = T, [LO.tolist(), HI.tolist()], list(NAMES)
    (ML / f"s1_npe{TAG}.json").write_text(json.dumps(res, indent=1, ensure_ascii=False), encoding="utf-8")


if __name__ == "__main__":
    main()
