"""The 10 Nov 2023 Grindavik dike with VARIABLE opening along strike (Sentinel-1, ASF HyP3).

The uniform-opening rectangle (s1_grindavik_dike.py) is too simple for the real dike: NPE and MCMC
disagreed. Here the dike is split into NSEG vertical segments along strike; the opening is given at
NNODE equally spaced nodes (log10 m, so a segment can almost close) and linearly interpolated to the
segment centres. The Svartsengi source uses the layered-crust FE Green's functions with topography
(as in s1_grindavik_dike.py --source fem). Same data, bins, noise covariance and offset/plane
marginalization as the uniform model.

Parameters (16): dike east, north, strike, length, top, height; log10 opening at 5 nodes (SW -> NE);
source east, north, depth b.s.l., dV for the A016 pair, dV for the D155 pair.

Steps: two NPE rounds (100 000 simulations each, in parallel processes, 256 PCA components),
calibration on held-out simulations, real-data posterior (timed), MCMC on the same likelihood (timed).
Writes output/reykjanes/ml/s1_dike_variable.json / .npz / .png.
Usage: python scripts/reykjanes/s1_dike_variable.py [--sims 100000] [--workers 4] [--mcmc-steps 2500] | --plot
"""
import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import torch
from pyproj import Transformer

sys.argv += [] if "--source" in sys.argv else ["--source", "fem"]
import s1_grindavik_dike as S
from ml_mcmc_source import log_post, stretch_sampler
from ml_npe_flow import from_z, posterior, to_z, train
from okada import tensile
from progress import Progress
from reykjanes_analysis import PLACES

ML = S.ML
NSEG, NNODE = 8, 5
NAMES = ("dike east km", "dike north km", "strike deg", "length km", "top km", "height km") + \
        tuple(f"log10 opening node {i + 1}" for i in range(NNODE)) + \
        ("src east km", "src north km", "src depth km", "src dV A016 1e6 m3", "src dV D155 1e6 m3")
LO = np.r_[[-8, -8, 0, 2, 0, 1], [-1.5] * NNODE, [-4, -4, 2, -150, -150]].astype(float)
HI = np.r_[[8, 8, 180, 24, 3, 8], [1.0] * NNODE, [4, 4, 8, 30, 30]].astype(float)
K = 256
RNG = np.random.default_rng(31)
_G = None


def arg(name, default):
    return type(default)(sys.argv[sys.argv.index(name) + 1]) if name in sys.argv else default


def seg_open(p):
    """Opening (m) of the NSEG segments, SW -> NE, from the NNODE log10 node values."""
    xs = (np.arange(NSEG) + 0.5) / NSEG
    return 10 ** np.interp(xs, np.linspace(0, 1, NNODE), p[6:6 + NNODE])


def forward(G, p):
    e0, n0, strike, L, top, H = p[:6]
    u_s = np.array([np.sin(np.deg2rad(strike)), np.cos(np.deg2rad(strike))])
    op = seg_open(p)
    dl = L / NSEG
    u = np.zeros((len(G.e), 3))
    for k in range(NSEG):
        c = (k + 0.5) * dl - L / 2                                    # along strike, SW (-) -> NE (+)
        u += tensile(G.e, G.n, e0 + c * u_s[0], n0 + c * u_s[1], top + H / 2, dl, H, strike, 90.0, op[k])
    u = u * 100 + G.source_enu(np.r_[np.zeros(7), p[11:16]])           # source params at indices 7..11 of the 12-vector
    return (u * G.vec).sum(1)


def _init():
    global _G
    _G = S.Geometry()


def _chunk(q):
    return np.array([forward(_G, p) for p in q])


def simulate(q, workers, label):
    out = []
    chunks = np.array_split(q, max(1, len(q) // 1000))
    with Progress(f"симуляції дайки зі змінним розкриттям, {label} ({len(q)})", total=len(q), unit="сим",
                  note=f"{workers} процесів") as pr, ProcessPoolExecutor(workers, initializer=_init) as ex:
        for r in ex.map(_chunk, chunks):
            out.append(r)
            pr.update(len(r))
    return np.concatenate(out)


class Model:
    lo, hi = LO, HI

    def __init__(self, G):
        self.G = G

    def predict(self, q):
        return np.array([forward(self.G, p) for p in q])

    def wrap(self, q):
        q = q.copy()
        q[:, 2] = np.mod(q[:, 2], 180)
        return q


def volume(s):
    s = np.atleast_2d(s)
    return np.array([(p[3] / NSEG) * p[5] * seg_open(p).sum() for p in s])     # km * km * m = 1e6 m3


def main():
    workers, n_sims = arg("--workers", 4), arg("--sims", 100000)
    G = S.Geometry()
    covs, noise_info, _ = S.fit_noise(G)
    G.set_noise(covs)
    yobs = G.y[None]
    T = {}

    def one_round(lo, hi, label, seed):
        t0 = time.time()
        q = lo + RNG.random((n_sims, len(LO))) * (hi - lo)
        Sx = G.whiten(simulate(q, workers, label))
        ts = time.time() - t0
        _, sv, Vt = np.linalg.svd(Sx[:20000] - Sx[:20000].mean(0), full_matrices=False)
        U = Vt[:K].T
        SP = Sx @ U
        scale = SP.std(0) + 1.0
        print(f"[{label}] {n_sims} simulations in {ts:.0f} s; {K} components keep {(sv[:K] ** 2).sum() / (sv ** 2).sum():.4%}; "
              f"component {K} = {sv[K - 1] / np.sqrt(20000):.2f} x noise", flush=True)
        del Sx
        net, tt = train(SP, to_z(q, LO, HI), None, arg("--epochs", 30), f"дайка зі змінним розкриттям, {label}", seed, scale=scale)
        t0 = time.time()
        ps = posterior([net], G.whiten(yobs) @ U, LO, HI, 6000, scale=scale)
        return dict(net=net, U=U, scale=scale, ps=ps, sim=ts, train=tt, infer=time.time() - t0)

    R1 = one_round(LO, HI, "раунд 1", 11)
    w = HI - LO
    a, b = np.quantile(R1["ps"], [0.0005, 0.9995], 0)
    lo2, hi2 = np.maximum(LO, a - 0.25 * (b - a)), np.minimum(HI, b + 0.25 * (b - a))
    if lo2[2] - LO[2] < 0.02 * w[2] or HI[2] - hi2[2] < 0.02 * w[2]:
        lo2[2], hi2[2] = LO[2], HI[2]
    print("truncated box: " + ", ".join(f"{nm} {l:.2f}..{h:.2f}" for nm, l, h in zip(NAMES, lo2, hi2)), flush=True)
    R = one_round(lo2, hi2, "раунд 2", 12)
    net, U, scale = R["net"], R["U"], R["scale"]
    T.update(simulate=R1["sim"] + R["sim"], train=R1["train"] + R["train"], infer_real=R["infer"])
    nt = 400
    qt = lo2 + RNG.random((nt, len(LO))) * (hi2 - lo2)
    xt = simulate(qt, workers, "калібрування") + G.sample_noise(nt)
    Xt = torch.tensor(G.whiten(xt) @ U / scale, dtype=torch.float32)
    with torch.no_grad():
        st = np.concatenate([from_z(net.flow(net.emb(Xt[j:j + 50])).sample((1000,)).numpy().astype(float), LO, HI)
                             for j in range(0, nt, 50)], axis=1)
    cal = {}
    for lev in (.68, .95):
        lo_, hi_ = np.quantile(st, [(1 - lev) / 2, (1 + lev) / 2], 0)
        cal[lev] = ((qt >= lo_) & (qt <= hi_)).mean(0).round(3).tolist()
    print("calibration 68 %:", cal[.68], flush=True)
    print("calibration 95 %:", cal[.95], flush=True)
    ps = R["ps"]
    model = Model(G)
    yw = G.whiten(yobs)[0]
    like = lambda pred: -0.5 * ((G.whiten(pred) - yw) ** 2).sum(1)
    pm = np.median(ps, 0)
    print(f"REAL posterior in {T['infer_real']:.2f} s; log-likelihood at the NPE median {like(model.predict(pm[None]))[0]:.1f}", flush=True)
    # MCMC, started from the best of 3000 prior draws (as for the uniform dike)
    cand = LO + RNG.random((3000, len(LO))) * (HI - LO)
    lp = log_post(model, like, cand)
    steps = arg("--mcmc-steps", 2500)
    t0 = time.time()
    chain, lps, acc = stretch_sampler(model, like, cand[np.argsort(lp)[-48:]], steps, 48, RNG, label="дайка зі змінним розкриттям")
    T["mcmc"] = time.time() - t0
    post = chain[steps // 2:]
    mlp = lps[steps // 2:].mean(0)
    pmc = post[:, mlp > mlp.max() - 15].reshape(-1, len(LO))
    pmm = np.median(pmc, 0)
    ll = dict(npe=round(float(like(model.predict(pm[None]))[0]), 1), mcmc=round(float(like(model.predict(pmm[None]))[0]), 1),
              npe_best_of_200=round(float(like(model.predict(ps[:200])).max()), 1))
    print(f"MCMC {T['mcmc']:.0f} s, acceptance {acc:.2f}; log-likelihood: {ll}", flush=True)
    for i, nm in enumerate(NAMES):
        print(f"   {nm:24s} NPE {np.percentile(ps[:, i], [16, 50, 84]).round(2)}   MCMC {np.percentile(pmc[:, i], [16, 50, 84]).round(2)}")
    Vn, Vm = volume(ps[:3000]), volume(pmc[RNG.choice(len(pmc), min(3000, len(pmc)), replace=False)])
    print(f"   dike volume 1e6 m3: NPE {np.percentile(Vn, [16, 50, 84]).round(0)}, MCMC {np.percentile(Vm, [16, 50, 84]).round(0)}")
    T = {k: round(v, 3) for k, v in T.items()}
    pct = lambda s: {nm: np.percentile(s[:, i], [16, 50, 84]).round(3).tolist() for i, nm in enumerate(NAMES)}
    omax = lambda s: np.array([seg_open(p).max() for p in s[:3000]])
    out = dict(model=f"{NSEG} segments, opening at {NNODE} nodes (log10 m), layered FE source", bins=dict(zip(S.TRACKS, G.sizes)),
               noise=noise_info, n_sims=n_sims, k_pca=K, timing=T, calibration=cal, npe=pct(ps), mcmc=pct(pmc),
               box=[lo2.tolist(), hi2.tolist()], loglike=ll,
               npe_volume=np.percentile(Vn, [16, 50, 84]).round(1).tolist(), mcmc_volume=np.percentile(Vm, [16, 50, 84]).round(1).tolist(),
               npe_max_opening=np.percentile(omax(ps), [16, 50, 84]).round(2).tolist(),
               mcmc_max_opening=np.percentile(omax(pmc), [16, 50, 84]).round(2).tolist(),
               published=S.PUBLISHED)
    (ML / "s1_dike_variable.json").write_text(json.dumps(out, indent=1, ensure_ascii=False), encoding="utf-8")
    np.savez_compressed(ML / "s1_dike_variable.npz", ps=ps[:4000], pmc=pmc[RNG.choice(len(pmc), min(4000, len(pmc)), replace=False)])
    plot(G, out)


def segs(p):
    u = np.array([np.sin(np.deg2rad(p[2])), np.cos(np.deg2rad(p[2]))])
    c = (np.arange(NSEG + 1) / NSEG - 0.5) * p[3]
    return p[0] + c * u[0], p[1] + c * u[1]


def plot(G=None, out=None):
    G = G or S.Geometry()
    out = out or json.loads((ML / "s1_dike_variable.json").read_text(encoding="utf-8"))
    Z = np.load(ML / "s1_dike_variable.npz")
    ps, pmc = Z["ps"], Z["pmc"]
    pm, pmm = np.median(ps, 0), np.median(pmc, 0)
    pred = forward(G, pm)
    to_utm = Transformer.from_crs(4326, 32627, always_xy=True)
    sx, sy = to_utm.transform(*PLACES["Svartsengi"])
    pl = {k: ((to_utm.transform(*PLACES[k])[0] - sx) / 1000, (to_utm.transform(*PLACES[k])[1] - sy) / 1000)
          for k in ("Svartsengi", "Grindavik", "Sundhnukur")}
    fig = plt.figure(figsize=(18, 11))
    gs = fig.add_gridspec(3, 4, height_ratios=[1, 1, .62], hspace=.4, wspace=.3)
    ttl_tr = {"A016": "висхідний трек 16, 31.10→12.11.2023", "D155": "низхідний трек 155, 10.11 (ранок)→22.11.2023"}
    i0 = 0
    for k, t in enumerate(S.TRACKS):
        n = G.sizes[k]
        e, nn, y, p = G.e[i0:i0 + n], G.n[i0:i0 + n], G.y[i0:i0 + n], pred[i0:i0 + n]
        A = G.design(k)
        c, *_ = np.linalg.lstsq(A, y - p, rcond=None)
        p = p + A @ c
        lim = np.percentile(np.abs(y - np.median(y)), 99)
        for j, (v, nm) in enumerate(((y, "дані Sentinel-1"), (p, "модель за медіаною NPE"), (y - p, "залишок"))):
            a = fig.add_subplot(gs[k, j])
            off = np.median(y) if j < 2 else 0
            sc = a.scatter(e, nn, c=v, cmap="RdBu_r", vmin=off - lim, vmax=off + lim, s=14, marker="s", linewidths=0)
            xs, ys = segs(pm)
            a.plot(xs, ys, "k-", lw=1.6)
            xs, ys = segs(pmm)
            a.plot(xs, ys, color="0.45", lw=1.3, ls="--")
            for nm2, (x0, y0) in pl.items():
                a.plot(x0, y0, "k^" if nm2 == "Svartsengi" else "ko", ms=4)
            a.set_aspect("equal"); a.set_xlim(S.BOX[0], S.BOX[1]); a.set_ylim(S.BOX[2], S.BOX[3])
            a.set_title(f"{ttl_tr[t]}\n{nm}", fontsize=9)
            a.tick_params(labelsize=7)
            plt.colorbar(sc, ax=a, fraction=0.035, pad=0.02).set_label("см", fontsize=8)
        i0 += n
    # map of the opening along the dike (NPE median), with the MCMC dike
    a = fig.add_subplot(gs[0, 3])
    for s in ps[:200]:
        xs, ys = segs(s)
        a.plot([xs[0], xs[-1]], [ys[0], ys[-1]], color="tab:blue", alpha=.04)
    xs, ys = segs(pm)
    op = seg_open(pm)
    for k in range(NSEG):
        a.plot(xs[k:k + 2], ys[k:k + 2], color=plt.cm.inferno_r(min(op[k] / 6, 1)), lw=2 + 2.2 * min(op[k], 6))
    xs, ys = segs(pmm)
    a.plot(xs, ys, color="0.45", lw=1.4, ls="--", label="MCMC (медіана)")
    for nm2, (x0, y0) in pl.items():
        a.plot(x0, y0, "k^" if nm2 == "Svartsengi" else "ko", ms=5)
        a.text(x0 + 0.6, y0 + (0.9 if nm2 == "Sundhnukur" else -0.9 if nm2 == "Svartsengi" else 0),
               {"Svartsengi": "Свартсенгі", "Grindavik": "Ґріндавік", "Sundhnukur": "Сундгнукур"}[nm2], fontsize=8)
    a.set_aspect("equal"); a.set_xlim(-12, 12); a.set_ylim(-12, 12)
    a.legend(fontsize=7, frameon=False, loc="lower right")
    a.set_title(f"NPE: дайка за {out['timing']['infer_real']:.1f} с\n(товщина й колір = розкриття)", fontsize=9)
    a.set_xlabel("схід від Свартсенгі, км", fontsize=8)
    # opening profile along strike
    a = fig.add_subplot(gs[1, 3])
    x = (np.arange(NSEG) + 0.5) / NSEG
    O = np.array([seg_open(p) for p in ps[:3000]])
    Om = np.array([seg_open(p) for p in pmc[:3000]])
    Ln, Lm = np.median(ps[:, 3]), np.median(pmc[:, 3])
    a.fill_between((x - .5) * Ln, *np.percentile(O, [16, 84], 0), color="tab:blue", alpha=.25, step="mid")
    a.step((x - .5) * Ln, np.median(O, 0), where="mid", color="tab:blue", label="NPE (медіана, 68 %)")
    a.fill_between((x - .5) * Lm, *np.percentile(Om, [16, 84], 0), color="0.5", alpha=.25, step="mid")
    a.step((x - .5) * Lm, np.median(Om, 0), where="mid", color="0.3", ls="--", label="MCMC")
    a.set_xlabel("вздовж дайки від центру, км (ПдЗх → ПнСх)", fontsize=8)
    a.set_ylabel("розкриття, м", fontsize=8)
    a.legend(fontsize=7, frameon=False)
    a.set_title("розкриття вздовж простягання", fontsize=9)
    a.tick_params(labelsize=7)
    a = fig.add_subplot(gs[2, :])
    a.axis("off")
    f = lambda v: f"{v[1]:.1f} [{v[0]:.1f}–{v[2]:.1f}]"
    N, Mc = out["npe"], out["mcmc"]
    rows = [["довжина, км", f(N["length km"]), f(Mc["length km"]), "~15"],
            ["об’єм дайки, млн м³", f(out["npe_volume"]), f(out["mcmc_volume"]), "130–139"],
            ["найбільше розкриття, м", f(out["npe_max_opening"]), f(out["mcmc_max_opening"]), "до ~8"],
            ["простягання, °", f(N["strike deg"]), f(Mc["strike deg"]), "—"],
            ["Свартсенгі dV (A016), млн м³", f(N["src dV A016 1e6 m3"]), f(Mc["src dV A016 1e6 m3"]), "−76…−82"],
            ["лог-правдоподібність у медіані", f"{out['loglike']['npe']:.0f}", f"{out['loglike']['mcmc']:.0f}", "—"],
            ["час відповіді", f"{out['timing']['infer_real']:.1f} с", f"{out['timing']['mcmc'] / 60:.0f} хв", "—"]]
    tb = a.table(cellText=rows, colLabels=["", "NPE", "MCMC", "опубліковано*"], loc="center", cellLoc="center", colWidths=[.25, .2, .2, .2])
    tb.auto_set_font_size(False); tb.set_fontsize(9); tb.scale(1, 1.45)
    a.text(0.5, -0.3, "медіана [16–84 %];  * Sigmundsson et al. 2024 (Science 383:1228), EGU24-14829; dV Свартсенгі опубліковано до 12.11",
           fontsize=9, ha="center", va="top", transform=a.transAxes)
    fig.suptitle(f"Дайка 10.11.2023 зі змінним розкриттям ({NSEG} сегментів): Sentinel-1 → NPE; Свартсенгі — шарувата кора + рельєф (МСЕ)", fontsize=13)
    fig.savefig(ML / "s1_dike_variable.png", dpi=105, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    print("wrote", ML / "s1_dike_variable.png")


if __name__ == "__main__":
    plot() if "--plot" in sys.argv else main()
