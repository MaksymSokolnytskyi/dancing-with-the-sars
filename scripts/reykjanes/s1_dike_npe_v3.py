"""NPE v3: the v2 dike model, but the atmospheric noise level is UNKNOWN and learned from the data.

Why: the noise of a 12-day Sentinel-1 pair here varies ~10x between pairs (quiet Oct 2023: 0.3-0.4 cm rms;
Nov 2023: ~4 cm, which does NOT scale with the deformation signal). NPE v2 used one fixed covariance fitted to the
Nov 2023 residuals, so it was blind to the small 2024 dikes. v3 draws, for every training simulation, its own noise
level (sill log-uniform 0.1-40 cm2, white 0-20 % of it) and correlation length (3-16 km) + offset/plane; GNSS:
base sigma + 5 % of the signal. Inputs are whitened with a fixed reference covariance (s1_common_v2.REF_NOISE);
the network infers the noise level from the data itself. Two noise realizations per simulation.

Outputs: s1_npe_v3.pt / .json / .npz (real 10 Nov 2023 posterior, SBC with random noise levels, PPC with the total
covariance), then the null tests and 2024 events (s1_npe_v3_apply.json, s1_null_tests_v3.png, s1_events_2024_v3.png).
Usage: python scripts/reykjanes/s1_dike_npe_v3.py [--sims 60000] [--epochs 20] [--workers 3]
"""
import json
import sys
import time

import numpy as np
import torch

sys.argv += [] if "--source" in sys.argv else ["--source", "fem"]
import s1_dike_npe_v2 as P2
import s1_dike_variable as V
import s1_grindavik_dike as S
from ml_npe_flow import Net, to_z, train
from s1_common_v2 import GeomV2, HierNoise

ML = S.ML
K = 256
RNG = np.random.default_rng(60)


class Feat:
    def __init__(self, Xw_fit, ni):
        _, sv, Vt = np.linalg.svd(Xw_fit[:, :ni] - Xw_fit[:, :ni].mean(0), full_matrices=False)
        self.U, self.ni = Vt[:K].T.astype(np.float32), ni
        self.info = f"{K} components keep {(sv[:K] ** 2).sum() / (sv ** 2).sum():.4%}"
        self.scale = None

    def __call__(self, Xw):
        return np.c_[Xw[:, :self.ni] @ self.U, Xw[:, self.ni:]].astype(np.float32)


def noisy_features(G, HN, clean, F=None, reps=2):
    out, chunks = [], []
    for r in range(reps):
        for j in range(0, len(clean), 5000):
            c = clean[j:j + 5000]
            n, _ = HN(c, RNG)
            xw = G.whiten(c + n).astype(np.float32)
            if F is None:
                F = Feat(xw, G.n_ins)
            out.append(F(xw))
    return np.concatenate(out), F


def main():
    workers, n_sims, epochs = P2.arg("--workers", 3), P2.arg("--sims", 60000), P2.arg("--epochs", 20)
    G = GeomV2(gnss=True, noise="ref")
    HN = HierNoise(G)
    lo, hi = P2.prior_box(G)
    T = {}
    t0 = time.time()
    q, nd = P2.draw(lo, hi, n_sims)
    clean = P2.simulate(q, workers, "v3, навчання").astype(np.float32)
    T["simulate"] = time.time() - t0
    X, F = noisy_features(G, HN, clean)
    del clean
    F.scale = X.std(0) + 1.0
    print(F.info, flush=True)
    Q = np.r_[to_z(q, lo, hi), to_z(q, lo, hi)]
    net, T["train"] = train(X, Q, None, epochs, "NPE v3 (рівень шуму невідомий)", 61, scale=F.scale)
    del X, Q
    torch.save(dict(state=net.state_dict(), U=F.U, scale=F.scale, lo=lo, hi=hi, ni=F.ni), ML / "s1_npe_v3.pt")
    # ---- real data
    t0 = time.time()
    ps = P2.sample(net, F, G.whiten(G.y[None]), lo, hi, 8000)[:, 0].astype(float)
    T["infer_real"] = time.time() - t0
    Vd = V.volume(ps)
    res = dict(prior=[lo.tolist(), hi.tolist()], n_sims=n_sims, reps=2, nodike_fraction=P2.NODIKE, features=F.info,
               noise="hierarchical: sill log-U(0.1, 40) cm2, range 3-16 km; whitening with the reference covariance", timing=T)
    res["real"] = dict(pct={nm: np.percentile(ps[:, i], [16, 50, 84]).round(3).tolist() for i, nm in enumerate(V.NAMES)},
                       volume=np.percentile(Vd, [2.5, 16, 50, 84, 97.5]).round(1).tolist(), p_nodike=round(float(np.mean(Vd < 5)), 4))
    print(f"REAL v3: {T['infer_real']:.1f} s; V {res['real']['volume']}; P(no dike) {res['real']['p_nodike']}; strike {res['real']['pct']['strike deg']}; "
          f"src depth {res['real']['pct']['src depth km']}; dV_A {res['real']['pct']['src dV A016 1e6 m3']}", flush=True)
    np.savez_compressed(ML / "s1_npe_v3.npz", ps=ps)
    # ---- SBC with random noise levels
    nt = 500
    qt, ndt = P2.draw(lo, hi, nt)
    ct = P2.simulate(qt, workers, "v3, SBC")
    nz, sill = HN(ct, RNG)
    st = P2.sample(net, F, G.whiten(ct + nz), lo, hi, 1000)
    ranks = (st < qt[None]).sum(0)
    Vt_true = V.volume(qt)
    Vs = np.array([V.volume(st[:, j].astype(float)) for j in range(nt)]).T
    from scipy.stats import chisquare
    unif = lambda r: float(chisquare(np.histogram(r, bins=20, range=(0, 1000))[0]).pvalue)
    pnd = np.mean(Vs < 5, 0)
    quiet = sill < 1.0
    res["sbc"] = dict(n=nt, ranks=ranks.tolist(), ranks_volume=(Vs < Vt_true[None]).sum(0).tolist(),
                      p_uniform={nm: unif(ranks[:, i]) for i, nm in enumerate(V.NAMES)}, p_uniform_volume=unif((Vs < Vt_true[None]).sum(0)),
                      levels=np.linspace(0.05, 0.95, 19).round(2).tolist(),
                      coverage=[((qt >= a) & (qt <= b)).mean(0).tolist() for a, b in
                                (np.quantile(st, [(1 - l) / 2, (1 + l) / 2], 0) for l in np.linspace(0.05, 0.95, 19))],
                      nodike=dict(n=int(ndt.sum()), p_nodike_median=round(float(np.median(pnd[ndt])), 3),
                                  p_nodike_median_quiet_noise=round(float(np.median(pnd[ndt & quiet])), 3) if (ndt & quiet).any() else None,
                                  p_nodike_on_dikes_median=round(float(np.median(pnd[~ndt & (Vt_true > 20)])), 3)),
                      small_dikes=dict(n=int(np.sum(~ndt & (Vt_true > 1) & (Vt_true < 30))),
                                       detected_frac_quiet=round(float(np.mean(pnd[~ndt & (Vt_true > 1) & (Vt_true < 30) & quiet] < 0.5)), 3)
                                       if (~ndt & (Vt_true > 1) & (Vt_true < 30) & quiet).any() else None))
    print("SBC v3 p:", {k: round(v, 3) for k, v in res["sbc"]["p_uniform"].items()}, "volume", round(res["sbc"]["p_uniform_volume"], 3), flush=True)
    print("SBC v3 no-dike / small dikes:", res["sbc"]["nodike"], res["sbc"]["small_dikes"], flush=True)
    # ---- PPC with the total (atmosphere + model error) covariance, compared with the v2 MCMC on the same covariance
    Gt = GeomV2(gnss=True, noise="auto")
    pp = ps[RNG.choice(len(ps), 200, replace=False)]
    pred = np.array([V.forward(Gt, p) for p in pp])
    rw = Gt.whiten(Gt.y[None] - pred)
    res["ppc"] = dict(chi2_median=round(float(np.median((rw ** 2).sum(1))), 1), n_data=int(rw.shape[1]),
                      chi2_at_median=round(float((Gt.whiten(Gt.y[None] - V.forward(Gt, np.median(ps, 0))[None]) ** 2).sum()), 1))
    v2 = json.loads((ML / "s1_npe_v2.json").read_text(encoding="utf-8"))
    res["mcmc_from_v2"] = v2["mcmc"]
    print("PPC v3:", res["ppc"], "| MCMC (v2, same covariance) chi2", -2 * v2["mcmc"]["gnss"]["loglike_median"], flush=True)
    (ML / "s1_npe_v3.json").write_text(json.dumps(res, indent=1, ensure_ascii=False), encoding="utf-8")


if __name__ == "__main__":
    main()
