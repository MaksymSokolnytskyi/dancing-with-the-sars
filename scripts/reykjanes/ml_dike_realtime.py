"""'Thousands of identical problems': dike inversion seconds after a new NISAR pair arrives.

Setting: a 12-day NISAR pair on both tracks (159A ascending, 082D descending) over the
Svartsengi - Sundhnukur - Grindavik area, plus 12-day GNSS displacements. A dike intrudes
somewhere in the window between the two acquisitions. Question: where is it, how long,
how deep, how much did it open, what volume?

Model: vertical rectangular tensile dislocation (Okada 1985, scripts/okada.py)
  east, north (centre, km from Svartsengi)  [-8, 8]
  strike (deg)                              [0, 180)
  length (km)                               [1, 16]
  top depth (km)                            [0, 3]
  height (km, bottom - top)                 [1, 7]
  log10 opening (m)                         [-1, 0.8]   (0.1 - 6.3 m)
Noise: real NISAR atmosphere. Training uses the Gaussian approximation (exponential + white
covariance fitted to the real per-date screens, rescaled from rate to a single pair);
TESTING uses the real screens themselves (two random dates, randomly rotated), so the test
is on noise the network never saw. GNSS: 0.5 cm per component.

Steps
  1. train once: 150 000 Okada simulations, whitened + PCA input, neural spline flow ensemble;
  2. 'live' test: 1000 new synthetic events with real atmosphere -> posterior for each,
     timed; accuracy and calibration;
  3. MCMC on two of these events (same likelihood), timed, for comparison;
  4. a showcase event of the scale of the 10 Nov 2023 Grindavik dike (synthetic!).

Writes output/reykjanes/ml/dike_realtime.json and dike_realtime.png.
Usage: python scripts/reykjanes/ml_dike_realtime.py [--sims 150000] [--events 1000]
"""
import json
import sys
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import torch
from scipy.spatial import cKDTree

import reykjanes_mogi as M
import ml_npe_source as NS
from ml_mcmc_source import Likelihood, log_post, rhat, stretch_sampler
from ml_npe_flow import K_PCA, from_z, posterior, to_z, train
from okada import tensile

ML = Path(__file__).resolve().parents[2] / "output" / "reykjanes" / "ml"
NAMES = ("east km", "north km", "strike deg", "length km", "top km", "height km", "opening m")
LO = np.array([-8, -8, 0, 1, 0, 1, -1.0])
HI = np.array([8, 8, 180, 16, 3, 7, 0.8])
SIGMA_GNSS_PAIR = 0.5
RNG = np.random.default_rng(17)
K_DIKE = 128     # dikes have sharp near-field patterns: more components than for a smooth inflation source


def arg(name, default):
    return type(default)(sys.argv[sys.argv.index(name) + 1]) if name in sys.argv else default


class DikeGeometry:
    def __init__(self):
        O, _ = NS.cached(NS.Obs)
        self.O = O
        self.sizes = [O.tr[t]["y"].size for t in M.TRACKS]
        self.e = np.concatenate([O.tr[t]["e"] for t in M.TRACKS] + [O.gn[:, 0], O.gref[:1]])
        self.n = np.concatenate([O.tr[t]["n"] for t in M.TRACKS] + [O.gn[:, 1], O.gref[1:]])
        self.vec = np.concatenate([O.tr[t]["vec"] for t in M.TRACKS])
        self.n_ins = sum(self.sizes)
        z = np.load(ML / "screens.npz")
        # covariance of ONE pair (cm^2) from the rate-noise fit: rate = pair * 365.25 / span
        self.cov = {t: z[f"cov_par_{t}"] * (float(z[f"span_{t}"]) / 365.25) ** 2 for t in M.TRACKS}
        self.screens = {t: z[f"raw_{t}"] for t in M.TRACKS}
        self.like = Likelihood(O, self.cov, True)
        self.like_gn = SIGMA_GNSS_PAIR
        self.maps = {}
        for t in M.TRACKS:
            e, n = O.tr[t]["e"], O.tr[t]["n"]
            tree = cKDTree(np.c_[e, n])
            self.maps[t] = np.array([tree.query(np.c_[np.cos(a) * f * e - np.sin(a) * n, np.sin(a) * f * e + np.cos(a) * n])[1]
                                     for f in (1, -1) for a in np.deg2rad(np.arange(0, 360, 5))])

    def clean(self, q):
        out = np.empty((len(q), self.n_ins + 3 * len(self.O.gn)))
        prog = None
        if len(q) >= 1000:
            from progress import Progress
            prog = Progress(f"симуляції дайок ({len(q)})", total=len(q), unit="сим")
        for i, (e0, n0, strike, L, top, H, lo) in enumerate(q):
            if prog and i % 500 == 499:
                prog.update(500)
            u = tensile(self.e, self.n, e0, n0, top + H / 2, L, H, strike, 90.0, 10 ** lo) * 100   # cm
            out[i, :self.n_ins] = (u[:self.n_ins] * self.vec).sum(1)
            out[i, self.n_ins:] = (u[self.n_ins:-1] - u[-1]).ravel()
        if prog:
            prog.close()
        return out

    def whiten(self, x):
        out, i0 = [], 0
        for n, W in zip(self.sizes, self.like.W):
            out.append(x[:, i0:i0 + n] @ W.T)
            i0 += n
        out.append(x[:, i0:] / SIGMA_GNSS_PAIR)
        return np.concatenate(out, 1)

    def real_noise(self, N):
        """Real atmosphere of one pair: two random dates of each track, rotated / mirrored, random sign."""
        parts = []
        for t in M.TRACKS:
            S = self.screens[t]
            i = RNG.integers(0, len(S), (N, 2))
            i[:, 1] = np.where(i[:, 1] == i[:, 0], (i[:, 0] + 1) % len(S), i[:, 1])
            d = S[i[:, 1]] - S[i[:, 0]]
            mp = self.maps[t][RNG.integers(0, len(self.maps[t]), N)]
            parts.append(np.take_along_axis(d, mp, 1) * RNG.choice([-1, 1], (N, 1)))
        parts.append(RNG.normal(0, SIGMA_GNSS_PAIR, (N, 3 * len(self.O.gn))))
        return np.concatenate(parts, 1)


class DikeModel:
    lo, hi = LO, HI

    def __init__(self, G):
        self.G = G

    def predict(self, q):
        return self.G.clean(q)

    def wrap(self, q):
        q = q.copy()
        q[:, 2] = np.mod(q[:, 2], 180)
        return q


class GaussLike:
    """Likelihood for MCMC on one event: InSAR whitening as in training, GNSS sigma 0.5 cm."""
    def __init__(self, G, x):
        self.G, self.xw = G, G.whiten(x[None])[0]

    def __call__(self, pred):
        r = self.G.whiten(pred) - self.xw
        return -0.5 * (r ** 2).sum(1)


def main():
    n_sims, n_events = arg("--sims", 150000), arg("--events", 1000)
    G = DikeGeometry()
    T = {}
    # ---------------- 1. train once
    t0 = time.time()
    q = LO + RNG.random((n_sims, 7)) * (HI - LO)
    S = G.whiten(G.clean(q))
    T["simulate"] = time.time() - t0
    print(f"{n_sims} Okada dike simulations in {T['simulate']:.0f} s", flush=True)
    _, sv, Vt = np.linalg.svd(S[:20000] - S[:20000].mean(0), full_matrices=False)
    U = Vt[:K_DIKE].T
    SP = S @ U
    scale = SP.std(0) + 1.0
    print(f"   {K_DIKE} components keep {(sv[:K_DIKE] ** 2).sum() / (sv ** 2).sum():.4%} of the signal variance", flush=True)
    nets, T["train"] = [], 0.0
    for k in range(arg("--ens", 1)):
        net, tt = train(SP, to_z(q, LO, HI), None, arg("--epochs", 40), f"дайка, потік #{k + 1}", 30 + k, scale=scale)
        nets.append(net)
        T["train"] += tt
    del S
    # ---------------- 2. live stream of new events with REAL atmosphere
    qe = LO + RNG.random((n_events, 7)) * (HI - LO)
    qe[:, 3] = RNG.uniform(3, 16, n_events)                          # events big enough to matter
    qe[:, 6] = RNG.uniform(-0.3, 0.8, n_events)
    xe = G.clean(qe) + G.real_noise(n_events)
    xk = G.whiten(xe) @ U
    t0 = time.time()
    X = torch.tensor(xk / scale, dtype=torch.float32)
    with torch.no_grad():                                            # in chunks of 50 events (memory)
        smp = np.concatenate([np.concatenate([from_z(net.flow(net.emb(X[j:j + 50])).sample((1200 // len(nets),)).numpy().astype(float), LO, HI)
                                              for net in nets]) for j in range(0, len(X), 50)], axis=1)   # (1200, n, 7)
    T["infer_all"] = time.time() - t0
    T["infer_one"] = T["infer_all"] / n_events
    print(f"posteriors for {n_events} new events in {T['infer_all']:.1f} s ({T['infer_one'] * 1000:.0f} ms per event)", flush=True)
    # strike is periodic: errors modulo 180
    med = np.median(smp, 0)
    err = med - qe
    err[:, 2] = (err[:, 2] + 90) % 180 - 90
    cover = {}
    for lev in (.68, .95):
        lo_q, hi_q = np.quantile(smp, [(1 - lev) / 2, (1 + lev) / 2], 0)
        cover[lev] = ((qe >= lo_q) & (qe <= hi_q)).mean(0).round(3).tolist()
    vol_true = qe[:, 3] * qe[:, 5] * 10 ** qe[:, 6] * 1e-3            # 1e6 m^3 (km * km * m)
    vol_post = smp[..., 3] * smp[..., 5] * 10 ** smp[..., 6] * 1e-3
    vol_err = np.median(vol_post, 0) / vol_true
    acc = {nm: round(float(np.median(np.abs(err[:, i]))), 3) for i, nm in enumerate(NAMES)}
    print("median absolute error: " + ", ".join(f"{k} {v}" for k, v in acc.items()))
    print(f"volume: median ratio estimate/true {np.median(vol_err):.2f}, 68 % of events within "
          f"{np.percentile(np.abs(np.log(vol_err)), 68):.2f} in log")
    print(f"calibration on REAL atmosphere, 68 %: {cover[.68]}; 95 %: {cover[.95]}")
    # ---------------- 3. MCMC on two events
    mc = []
    for j in range(2):
        like = GaussLike(G, xe[j])
        model = DikeModel(G)
        cand = LO + RNG.random((2000, 7)) * (HI - LO)
        lp = log_post(model, like, cand)
        t0 = time.time()
        steps = arg("--mcmc-steps", 2500)
        chain, lps, a = stretch_sampler(model, like, cand[np.argsort(lp)[-48:]], steps, 48, RNG, label=f"дайка, подія {j}")
        tt = time.time() - t0
        post = chain[steps // 3:]
        mlp = lps[steps // 3:].mean(0)
        post = post[:, mlp > np.median(mlp) - 10].reshape(-1, 7)
        mc.append(dict(seconds=round(tt), n_forward=steps * 48, true=qe[j].round(3).tolist(),
                       mcmc_median=np.median(post, 0).round(3).tolist(), npe_median=med[j].round(3).tolist(),
                       mcmc_68=np.percentile(post, [16, 84], 0).round(3).tolist(),
                       npe_68=np.percentile(smp[:, j], [16, 84], 0).round(3).tolist()))
        print(f"MCMC event {j}: {tt:.0f} s; true {np.round(qe[j], 2).tolist()}\n   MCMC {np.round(np.median(post, 0), 2).tolist()}\n"
              f"   NPE  {np.round(med[j], 2).tolist()}", flush=True)
    T["mcmc_one"] = float(np.mean([m["seconds"] for m in mc]))
    # ---------------- 4. showcase: a dike of the scale of 10 Nov 2023 (synthetic)
    show = np.array([[-0.5, -2.0, 35.0, 14.0, 0.8, 4.0, np.log10(2.0)]])
    xs = G.clean(show) + G.real_noise(1)
    t0 = time.time()
    ps = posterior(nets, G.whiten(xs) @ U, LO, HI, 6000, scale=scale)
    T["infer_show"] = time.time() - t0
    T = {k: round(v, 3) for k, v in T.items()}
    print(f"showcase: posterior in {T['infer_show']:.2f} s; length {np.percentile(ps[:, 3], [16, 50, 84]).round(1)} km, "
          f"opening {np.percentile(10 ** ps[:, 6], [16, 50, 84]).round(2)} m, strike {np.percentile(ps[:, 2], [16, 50, 84]).round(0)}")
    out = dict(timing=T, n_sims=n_sims, n_events=n_events, median_abs_error=acc, calibration=cover,
               volume_ratio_median=round(float(np.median(vol_err)), 3), mcmc=mc,
               mcmc_for_all_events_hours=round(T["mcmc_one"] * n_events / 3600, 1),
               showcase=dict(true=show[0].tolist(), p16_50_84=np.percentile(ps, [16, 50, 84], 0).round(3).tolist()))
    (ML / "dike_realtime.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    lo_q, hi_q = np.quantile(smp, [.16, .84], 0)
    np.savez_compressed(ML / "dike_realtime_plot.npz", show=show, xs=xs, ps=ps, qe=qe, med=med, q16=lo_q, q84=hi_q)
    plot(G)


def dike_line(s):
    u = np.array([np.sin(np.deg2rad(s[2])), np.cos(np.deg2rad(s[2]))])
    return [s[0] - u[0] * s[3] / 2, s[0] + u[0] * s[3] / 2], [s[1] - u[1] * s[3] / 2, s[1] + u[1] * s[3] / 2]


def plot(G=None):
    """Figure from the saved results (dike_realtime_plot.npz + .json): no retraining needed."""
    G = G or DikeGeometry()
    D = np.load(ML / "dike_realtime_plot.npz")
    J = json.loads((ML / "dike_realtime.json").read_text(encoding="utf-8"))
    show, xs, ps, qe, med, q16, q84 = (D[k] for k in ("show", "xs", "ps", "qe", "med", "q16", "q84"))
    T, n_events = J["timing"], J["n_events"]
    pm = np.median(ps, 0)
    model = G.clean(pm[None])[0]
    fig = plt.figure(figsize=(18, 10.5))
    gs = fig.add_gridspec(2, 4, hspace=.42, wspace=.45)
    i0 = 0
    lim = np.percentile(np.abs(xs[0, :G.n_ins]), 98)
    for k, t in enumerate(M.TRACKS):
        n = G.sizes[k]
        e, nn = G.O.tr[t]["e"], G.O.tr[t]["n"]
        for col, (vals, ttl) in enumerate(((xs[0, i0:i0 + n], "дані (дайка + атмосфера)"), (model[i0:i0 + n], "модель за медіаною NPE"))):
            if k == 1 and col == 1:
                continue
            a = fig.add_subplot(gs[0, k if col == 0 else 2])
            sc = a.scatter(e, nn, c=vals, cmap="RdBu_r", vmin=-lim, vmax=lim, s=24, marker="s", linewidths=0)
            a.plot(*dike_line(show[0]), "k-", lw=2.2, label="справжня дайка")
            a.plot(*dike_line(pm), color="gold", ls="--", lw=1.8, label="NPE (медіана)")
            a.plot(0, 0, "k^", ms=6)
            a.set_aspect("equal"); a.set_xlim(-12, 12); a.set_ylim(-12, 12)
            a.set_title(f"{t}: {ttl}", fontsize=9)
            a.set_xlabel("схід від Свартсенгі, км", fontsize=8)
            if k == 0 and col == 0:
                a.set_ylabel("північ, км", fontsize=8)
                a.legend(fontsize=7, loc="lower left", framealpha=.85)
            plt.colorbar(sc, ax=a, fraction=0.046, pad=0.03).set_label("см", fontsize=8)
        i0 += n
    a = fig.add_subplot(gs[0, 3])
    for s in ps[:250]:
        a.plot(*dike_line(s), color="tab:blue", alpha=.06)
    a.plot(*dike_line(show[0]), "r-", lw=2.5, label="справжня дайка")
    a.plot(0, 0, "k^", label="Свартсенгі")
    a.set_aspect("equal"); a.set_xlim(-12, 12); a.set_ylim(-12, 12); a.legend(fontsize=7, frameon=False, loc="lower left")
    a.set_title(f"250 можливих дайок за {T['infer_show']:.1f} с\n"
                f"довжина {np.percentile(ps[:, 3], 50):.1f} км [{np.percentile(ps[:, 3], 16):.1f}–{np.percentile(ps[:, 3], 84):.1f}]\n"
                f"розкриття {np.percentile(10 ** ps[:, 6], 50):.2f} м (правда: 14 км, 2 м)", fontsize=9)
    a.set_xlabel("схід, км", fontsize=8)
    for j, (i, nm, f) in enumerate(((3, "довжина, км", lambda v: v), (6, "розкриття, м", lambda v: 10 ** v), (4, "верх дайки, км", lambda v: v))):
        a = fig.add_subplot(gs[1, j])
        m, lo_, hi_ = f(med[:300, i]), f(q16[:300, i]), f(q84[:300, i])
        a.errorbar(f(qe[:300, i]), m, yerr=[m - lo_, hi_ - m], fmt=".", ms=3, alpha=.4, elinewidth=.6)
        lim2 = [f(qe[:, i]).min(), f(qe[:, i]).max()]
        a.plot(lim2, lim2, "k--", lw=.8)
        a.set_xlabel(f"правда: {nm}", fontsize=8); a.set_ylabel("оцінка NPE (медіана, 68 %)", fontsize=8)
        a.set_title(f"{nm}: 300 з {n_events} тестових подій", fontsize=9)
    a = fig.add_subplot(gs[1, 3])
    bars = [("NPE: 1 подія", T["infer_one"]), (f"NPE: {n_events} подій", T["infer_all"]),
            ("MCMC: 1 подія", T["mcmc_one"]), (f"MCMC: {n_events} подій", T["mcmc_one"] * n_events),
            ("навчання NPE (раз)", T["simulate"] + T["train"])]
    vals = [b[1] for b in bars]
    a.barh([b[0] for b in bars], vals, color=["tab:blue", "tab:blue", "tab:red", "tab:red", "0.6"])
    a.set_xscale("log"); a.set_xlim(min(vals) / 3, max(vals) * 60)
    a.set_xlabel("секунди (лог. шкала)", fontsize=8)
    a.tick_params(labelsize=8)
    for k, v in enumerate(vals):
        a.text(v * 1.4, k, f"{v:.2f} с" if v < 10 else f"{v:.0f} с" if v < 60 else f"{v / 60:.0f} хв" if v < 7200 else f"{v / 3600:.0f} год", va="center", fontsize=8)
    a.set_title("скільки часу коштує відповідь", fontsize=9)
    fig.suptitle("Інверсія дайки за секунди (синтетична подія масштабу Ґріндавіка 10.11.2023, справжня атмосфера NISAR): "
                 "колір = незгорнуте зміщення, см", fontsize=12)
    fig.savefig(ML / "dike_realtime.png", dpi=105, facecolor="white", bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    plot() if "--plot" in sys.argv else main()
