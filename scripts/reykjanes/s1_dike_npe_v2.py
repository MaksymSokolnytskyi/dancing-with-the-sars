"""NPE v2 for the 10 Nov 2023 Grindavik dike: InSAR (graben masked) + GNSS, a "no dike" class, SBC.

Model as s1_dike_variable.py: 8 vertical segments along strike, log10 opening at 5 nodes, Svartsengi source in
the layered crust (FE), separate source volume change for each window. Data: s1_common_v2.GeomV2.

  prior     corridor around the dike trace already found (uniform-dike NPE, confirmed by the HLS vents):
            centre +-4 km, strike 0-60 deg, length 2-24 km; --wide: centre +-15 km, strike 0-180, length 2-30
  no dike   NODIKE of the simulations have all opening nodes in [-1.5, -1.2] (<= 6 cm: no signal), so the
            network can answer "there is no dike"; reported as P(V_dike < 5e6 m3)
  network   one round, N simulations, whitened InSAR -> 256 PCA components, whitened GNSS appended,
            neural spline flow (ml_npe_flow.train), fresh N(0, I) noise every epoch
  checks    SBC on 500 prior draws (ranks + coverage); real 10 Nov 2023 posterior; MCMC (start: NPE samples)
            with and without GNSS -> H1; posterior predictive residuals -> PPC

Writes output/reykjanes/ml/s1_npe_v2{_wide}.json / .npz / .pt.
Usage: python scripts/reykjanes/s1_dike_npe_v2.py [--wide] [--sims 60000] [--epochs 20] [--workers 4] [--mcmc-steps 1500]
"""
import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import torch

sys.argv += [] if "--source" in sys.argv else ["--source", "fem"]
import s1_grindavik_dike as S
import s1_dike_variable as V
from ml_mcmc_source import stretch_sampler
from ml_npe_flow import from_z, to_z, train
from progress import Progress
from s1_common_v2 import GeomV2

ML = S.ML
WIDE = "--wide" in sys.argv
TAG = "_wide" if WIDE else ""
NODIKE = 0.25
K = 256
RNG = np.random.default_rng(41 if WIDE else 40)
NAMES = V.NAMES
_G = None


def arg(name, default):
    return type(default)(sys.argv[sys.argv.index(name) + 1]) if name in sys.argv else default


def prior_box(G):
    e0, n0, st = G.trace[0], G.trace[1], G.trace[2]
    if WIDE:
        lo = [-15, -12, 0, 2, 0, 1]; hi = [15, 12, 180, 30, 3, 8]
    else:
        lo = [e0 - 4, n0 - 4, max(0, st - 30), 2, 0, 1]; hi = [e0 + 4, n0 + 4, st + 30, 24, 3, 8]
    lo += [-1.5] * V.NNODE + [-4, -4, 2, -200, -200]
    hi += [1.0] * V.NNODE + [4, 4, 8, 30, 30]
    return np.array(lo, float), np.array(hi, float)


def draw(lo, hi, n):
    q = lo + RNG.random((n, len(lo))) * (hi - lo)
    nd = RNG.random(n) < NODIKE
    q[nd, 6:6 + V.NNODE] = RNG.uniform(-1.5, -1.2, (nd.sum(), V.NNODE))
    return q, nd


def _init():
    global _G
    _G = GeomV2(gnss=True)


def _chunk(q):
    return np.array([V.forward(_G, p) for p in q])


def simulate(q, workers, label):
    out = []
    chunks = np.array_split(q, max(1, len(q) // 1000))
    with Progress(f"NPE v2{TAG}: симуляції, {label} ({len(q)})", total=len(q), unit="сим", note=f"{workers} процесів") as pr, \
            ProcessPoolExecutor(workers, initializer=_init) as ex:
        for r in ex.map(_chunk, chunks):
            out.append(r)
            pr.update(len(r))
    return np.concatenate(out)


class Features:
    """Whitened InSAR -> PCA (U), whitened GNSS appended; per-feature scale for the network."""
    def __init__(self, G, Xw):
        ni = G.n_ins
        _, sv, Vt = np.linalg.svd(Xw[:20000, :ni] - Xw[:20000, :ni].mean(0), full_matrices=False)
        self.U, self.ni = Vt[:K].T, ni
        self.info = f"{K} components keep {(sv[:K] ** 2).sum() / (sv ** 2).sum():.4%}; component {K} = {sv[K - 1] / np.sqrt(20000):.2f} x noise"
        F = self(Xw)
        self.scale = F.std(0) + 1.0

    def __call__(self, Xw):
        return np.c_[Xw[:, :self.ni] @ self.U, Xw[:, self.ni:]]


def sample(net, F, Xw, lo, hi, n):
    X = torch.tensor(F(Xw) / F.scale, dtype=torch.float32)
    with torch.no_grad():
        return np.concatenate([from_z(net.flow(net.emb(X[j:j + 10])).sample((n,)).numpy().astype(np.float32), lo, hi)
                               for j in range(0, len(X), 10)], axis=1)                      # (n, N, d); small chunks: low memory


class SavedFeatures:
    def __init__(self, ck):
        self.U, self.ni, self.scale = ck["U"], ck["ni"], ck["scale"]
        self.info = f"loaded {K} components"

    def __call__(self, Xw):
        return np.c_[Xw[:, :self.ni] @ self.U, Xw[:, self.ni:]]


def main():
    workers, n_sims = arg("--workers", 4), arg("--sims", 60000 if not WIDE else 30000)
    G = GeomV2(gnss=True)
    lo, hi = prior_box(G)
    T = {}
    if "--post" in sys.argv:                       # reuse the trained network: SBC, MCMC, PPC only
        from ml_npe_flow import Net
        ck = torch.load(ML / f"s1_npe_v2{TAG}.pt", weights_only=False)
        lo, hi = ck["lo"], ck["hi"]
        F = SavedFeatures(ck)
        net = Net(F.U.shape[1] + G.n_gnss, len(lo))
        net.load_state_dict(ck["state"]); net.eval()
        old = json.loads((ML / "s1_npe_v2_log_timing.json").read_text()) if (ML / "s1_npe_v2_log_timing.json").exists() else {}
        T.update(old)
    else:
        t0 = time.time()
        q, nd = draw(lo, hi, n_sims)
        Xw = G.whiten(simulate(q, workers, "навчання"))
        T["simulate"] = time.time() - t0
        F = Features(G, Xw)
        print(F.info, flush=True)
        net, T["train"] = train(F(Xw), to_z(q, lo, hi), None, arg("--epochs", 20 if not WIDE else 12), f"NPE v2{TAG}", 51, scale=F.scale)
        del Xw, q
        torch.save(dict(state=net.state_dict(), U=F.U, scale=F.scale, lo=lo, hi=hi, ni=F.ni), ML / f"s1_npe_v2{TAG}.pt")
        if not WIDE:
            (ML / "s1_npe_v2_log_timing.json").write_text(json.dumps(T))
    # ---- real data
    t0 = time.time()
    ps = sample(net, F, G.whiten(G.y[None]), lo, hi, 8000)[:, 0]
    T["infer_real"] = time.time() - t0
    Vd = V.volume(ps)
    res = dict(prior=[lo.tolist(), hi.tolist()], wide=WIDE, n_sims=n_sims, nodike_fraction=NODIKE, features=F.info, noise=G.noise_source,
               n_insar=G.n_ins, n_gnss=G.n_gnss, timing=T)
    res["real"] = dict(pct={nm: np.percentile(ps[:, i], [16, 50, 84]).round(3).tolist() for i, nm in enumerate(NAMES)},
                       volume=np.percentile(Vd, [2.5, 16, 50, 84, 97.5]).round(1).tolist(), p_nodike=round(float(np.mean(Vd < 5)), 4),
                       max_opening=np.percentile([V.seg_open(p).max() for p in ps[:3000]], [16, 50, 84]).round(2).tolist())
    print(f"REAL: posterior in {T['infer_real']:.1f} s; V {res['real']['volume']} 1e6 m3; P(no dike) {res['real']['p_nodike']}; "
          f"strike {res['real']['pct']['strike deg']}; src depth {res['real']['pct']['src depth km']}; dV_A {res['real']['pct']['src dV A016 1e6 m3']}", flush=True)
    np.savez_compressed(ML / f"s1_npe_v2{TAG}.npz", ps=ps)
    if WIDE:
        (ML / f"s1_npe_v2{TAG}.json").write_text(json.dumps(res, indent=1, ensure_ascii=False), encoding="utf-8")
        return
    # ---- SBC on 500 prior draws (incl. no-dike), Gaussian noise of the same covariance
    nt = 500
    qt, ndt = draw(lo, hi, nt)
    Xt_w = G.whiten(simulate(qt, workers, "SBC") + G.sample_noise(nt, RNG))
    st = sample(net, F, Xt_w, lo, hi, 1000)                                            # (1000, nt, d), float32
    del Xt_w
    ranks = (st < qt[None]).sum(0)                                                     # (nt, d)
    Vt_true = V.volume(qt)
    Vs = np.array([V.volume(st[:, j]) for j in range(nt)]).T                           # (1000, nt)
    ranks_V = (Vs < Vt_true[None]).sum(0)
    from scipy.stats import chisquare
    def unif_p(r, bins=20):
        h = np.histogram(r, bins=bins, range=(0, 1000))[0]
        return float(chisquare(h).pvalue)
    levels = np.linspace(0.05, 0.95, 19)
    cover = []
    for lev in levels:
        a, b = np.quantile(st, [(1 - lev) / 2, (1 + lev) / 2], 0)
        cover.append(((qt >= a) & (qt <= b)).mean(0).tolist())
    Va, Vb = np.quantile(Vs, [0.16, 0.84], 0)
    res["sbc"] = dict(n=nt, ranks=ranks.tolist(), ranks_volume=ranks_V.tolist(), p_uniform={nm: unif_p(ranks[:, i]) for i, nm in enumerate(NAMES)},
                      p_uniform_volume=unif_p(ranks_V), levels=levels.round(2).tolist(), coverage=cover,
                      coverage_volume_68=round(float(np.mean((Vt_true >= Va) & (Vt_true <= Vb))), 3),
                      nodike=dict(n=int(ndt.sum()), p_nodike_median=round(float(np.median(np.mean(Vs[:, ndt] < 5, 0))), 3),
                                  p_nodike_on_dikes_median=round(float(np.median(np.mean(Vs[:, ~ndt & (Vt_true > 20)] < 5, 0))), 3)))
    print("SBC uniformity p:", {k: round(v, 3) for k, v in res["sbc"]["p_uniform"].items()}, "volume", round(res["sbc"]["p_uniform_volume"], 3), flush=True)
    print(f"SBC no-dike cases: median P(no dike) {res['sbc']['nodike']['p_nodike_median']}; on dikes > 20e6 m3: {res['sbc']['nodike']['p_nodike_on_dikes_median']}", flush=True)
    # ---- MCMC with and without GNSS (start: NPE samples), same model
    steps = arg("--mcmc-steps", 1500)
    res["mcmc"] = {}
    for use_g in (True, False):
        Gm = G if use_g else GeomV2(gnss=False)
        yw = Gm.whiten(Gm.y[None])[0]
        like = lambda pred, Gm=Gm, yw=yw: -0.5 * ((Gm.whiten(pred) - yw) ** 2).sum(1)
        model = V.Model(Gm)
        model.lo, model.hi = lo, hi
        start = ps[RNG.choice(len(ps), 48, replace=False)]
        t0 = time.time()
        ch, lp, acc = stretch_sampler(model, like, start, steps, 48, RNG, label=f"NPE v2: MCMC {'з GNSS' if use_g else 'без GNSS'}")
        post = ch[steps // 2:]
        mlp = lp[steps // 2:].mean(0)
        pmc = post[:, mlp > mlp.max() - 15].reshape(-1, len(lo))
        sub = pmc[RNG.choice(len(pmc), min(3000, len(pmc)), replace=False)]
        Vm = V.volume(sub)
        key = "gnss" if use_g else "insar_only"
        res["mcmc"][key] = dict(seconds=round(time.time() - t0), acceptance=round(float(acc), 3),
                                pct={nm: np.percentile(pmc[:, i], [16, 50, 84]).round(3).tolist() for i, nm in enumerate(NAMES)},
                                volume=np.percentile(Vm, [2.5, 16, 50, 84, 97.5]).round(1).tolist(),
                                corr_V_depth=round(float(np.corrcoef(Vm, sub[:, 13])[0, 1]), 3),
                                corr_V_dVA=round(float(np.corrcoef(Vm, sub[:, 14])[0, 1]), 3),
                                loglike_median=round(float(like(model.predict(np.median(pmc, 0)[None]))[0]), 1))
        np.save(ML / f"s1_npe_v2_mcmc_{key}.npy", sub)
        print(f"MCMC {key}: {res['mcmc'][key]['seconds']} s; V {res['mcmc'][key]['volume']}; corr(V, depth) {res['mcmc'][key]['corr_V_depth']}, "
              f"corr(V, dV_A) {res['mcmc'][key]['corr_V_dVA']}; logL {res['mcmc'][key]['loglike_median']}", flush=True)
    Vn = V.volume(ps[:3000])
    res["real"]["corr_V_depth"] = round(float(np.corrcoef(Vn, ps[:3000, 13])[0, 1]), 3)
    res["real"]["corr_V_dVA"] = round(float(np.corrcoef(Vn, ps[:3000, 14])[0, 1]), 3)
    yw = G.whiten(G.y[None])[0]
    res["real"]["loglike_median"] = round(float(-0.5 * ((G.whiten(V.Model(G).predict(np.median(ps, 0)[None])) - yw) ** 2).sum()), 1)
    # ---- PPC: 200 posterior samples -> residuals; projection on the dike Jacobian vs pure noise
    pp = ps[RNG.choice(len(ps), 200, replace=False)]
    pred = np.array([V.forward(G, p) for p in pp])
    rw = G.whiten(G.y[None] - pred)                                                    # whitened residuals (200, n)
    pm = np.median(ps, 0)
    J = []
    for i in range(6 + V.NNODE):                                                       # dike parameters only
        d = np.zeros_like(pm); d[i] = 1e-3 * (hi[i] - lo[i])
        J.append((G.whiten(V.forward(G, pm + d)[None]) - G.whiten(V.forward(G, pm - d)[None]))[0] / (2 * d[i]))
    J = np.array(J).T
    Qj, _ = np.linalg.qr(J)
    proj = np.linalg.norm(rw @ Qj, axis=1)                                             # residual energy along dike directions
    proj_noise = np.linalg.norm(G.whiten(G.sample_noise(500, RNG)) @ Qj, axis=1)
    res["ppc"] = dict(chi2_median=round(float(np.median((rw ** 2).sum(1))), 1), n_data=int(rw.shape[1]),
                      dike_projection_median=round(float(np.median(proj)), 2),
                      dike_projection_noise_95=round(float(np.percentile(proj_noise, 95)), 2),
                      median_residual=np.median(G.y[None] - pred, 0).round(2).tolist())
    print(f"PPC: chi2 {res['ppc']['chi2_median']} / {res['ppc']['n_data']}; residual along dike directions {res['ppc']['dike_projection_median']} "
          f"(pure noise 95 %: {res['ppc']['dike_projection_noise_95']})", flush=True)
    (ML / "s1_npe_v2.json").write_text(json.dumps(res, indent=1, ensure_ascii=False), encoding="utf-8")


if __name__ == "__main__":
    main()
