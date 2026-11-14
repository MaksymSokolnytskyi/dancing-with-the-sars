"""Change points in the GP-cleaned NISAR vertical series over Reykjanes, as Bayes factors.

Null model per cell (the Gaussian process of ml_gp_series.py, round 2):
    y = c + b t + f(t) + e,   f ~ GP(sigma_f, l_t),  e white (sigma_n)          -> covariance C0
Alternatives add one more term with a Gaussian prior on its size:
    step   s * H(t - tau),        s ~ N(0, 5 cm)          sudden offset (an intrusion deflates
                                                         Svartsengi within hours)
    hinge  d * (t - tau)+,        d ~ N(0, 20 cm/yr)     change of rate
The log Bayes factor has a closed form (matrix determinant lemma / Sherman-Morrison), so all
cells and all candidate dates are tested at once. We take the best tau and subtract log(#tau)
for having searched; a change is reported when the corrected log BF > log(20) (strong).

Injection test: a synthetic intrusion (Svartsengi deflates by A cm with the shape of the source
model) is added to the raw series on one date; the full GP cleaning is rerun and we record
whether the dome centre is flagged at the right date, and the false-alarm rate far from the
dome (> 12 km) and on the real data.

Usage: python scripts/reykjanes/ml_changepoint.py
"""
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from ml_gp_series import gp_clean, temporal_cov

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "output" / "reykjanes"
ML = OUT / "ml"
EDGE = 2
PRIOR = {"step": 5.0, "hinge": 20.0}
THRESH = np.log(20.0)


def bayes_factors(Y, t, hyper):
    """Y: (cells, T) complete rows. Returns dict kind -> (logBF_corrected, tau index, posterior size)."""
    sf, lt, sn = hyper
    C0 = temporal_cov(t, sf, lt, sn)
    Ci = np.linalg.inv(C0)
    taus = list(range(EDGE, t.size - EDGE))
    out = {}
    for kind, s0 in PRIOR.items():
        best = np.full(Y.shape[0], -np.inf)
        btau = np.zeros(Y.shape[0], int)
        bsize = np.zeros(Y.shape[0])
        for tau in taus:
            h = (t >= t[tau]).astype(float) if kind == "step" else np.maximum(t - t[tau], 0)
            a = h @ Ci @ h
            b = Y @ (Ci @ h)
            lbf = -0.5 * np.log1p(s0 ** 2 * a) + 0.5 * s0 ** 2 * b ** 2 / (1 + s0 ** 2 * a)
            size = s0 ** 2 * b / (1 + s0 ** 2 * a)              # posterior mean of s or d
            better = lbf > best
            best[better], btau[better], bsize[better] = lbf[better], tau, size[better]
        out[kind] = (best - np.log(len(taus)), btau, bsize)
    return out


def analyse(up, t, d):
    T, ny, nx = up.shape
    model = d["model_up"].astype(float)
    pi, pj = np.unravel_index(np.nanargmax(model), model.shape)
    yy, xx = np.mgrid[:ny, :nx]
    protect = np.hypot(yy - pi, xx - pj) * 0.24 < 8.0
    lon, lat = d["lon"].astype(float), d["lat"].astype(float)
    kx, ky = 111.32 * np.cos(np.deg2rad(63.88)), 111.2
    Xk, Yk = (lon - lon.mean()) * kx, (lat - lat.mean()) * ky
    clean, smooth, rate, rsd, screens, hyper, *_ = gp_clean(up.reshape(T, -1).T, t, Xk, Yk, ny, nx, verbose=False, protect=protect)
    ok = np.isfinite(clean).all(1)
    Y = clean[ok] - clean[ok][:, :1]
    bf = bayes_factors(Y, t, hyper[6:9])
    res = {}
    for kind, (lbf, tau, size) in bf.items():
        for name, v in (("lbf", lbf), ("tau", tau), ("size", size)):
            a = np.full(ny * nx, np.nan)
            a[ok] = v
            res[f"{kind}_{name}"] = a.reshape(ny, nx)
    return res, hyper


def report(res, d, label, dome, far):
    lines = []
    for kind in ("step", "hinge"):
        lbf = res[f"{kind}_lbf"]
        m = lbf > THRESH
        fin = np.isfinite(lbf)
        txt = f"{kind}: flagged {np.sum(m)} of {np.sum(fin)} cells ({np.mean(m[fin]):.1%}); far field (> 12 km) {np.mean(m[far & fin]):.1%}"
        if m.any():
            top = np.bincount(res[f"{kind}_tau"][m].astype(int)).argmax()
            txt += f"; most common date {d['dates'][top]}, median size {np.nanmedian(res[f'{kind}_size'][m]):+.1f} {'cm' if kind == 'step' else 'cm/yr'}"
        dl = np.nanmedian(lbf[dome])
        txt += f"; dome centre log BF {dl:+.1f} at {d['dates'][int(np.nanmedian(res[f'{kind}_tau'][dome]))]}" \
               f" (size {np.nanmedian(res[f'{kind}_size'][dome]):+.1f})"
        lines.append(txt)
        print(f"  [{label}] {txt}")
    return lines


def main():
    ML.mkdir(parents=True, exist_ok=True)
    d = np.load(OUT / "breathing_data.npz")
    up = d["up"].astype(float)
    dates = [datetime.fromisoformat(str(x)) for x in d["dates"]]
    t = np.array([(x - dates[0]).days / 365.25 for x in dates])
    model = d["model_up"].astype(float)
    shape = np.where(model > 0.05 * np.nanmax(model), model / np.nanmax(model), 0.0)
    pi, pj = np.unravel_index(np.nanargmax(model), model.shape)
    ny, nx = model.shape
    yy, xx = np.mgrid[:ny, :nx]
    dist = np.hypot(yy - pi, xx - pj) * 0.24
    dome, far = dist < 1.0, dist > 12

    print(f"log Bayes factor threshold {THRESH:.1f} (BF 20) after the look-elsewhere correction")
    res, hyper = analyse(up, t, d)
    print(f"GP hyperparameters used: sigma_f {hyper[6]:.2f} cm, length {hyper[7] * 365.25:.0f} d, white {hyper[8]:.2f} cm")
    report(res, d, "real data", dome, far)

    inj = int(np.argmin(np.abs(t - t[-1] * 0.6)))
    print(f"\ninjection tests: synthetic intrusion (Svartsengi deflates) between {d['dates'][inj - 1]} and {d['dates'][inj]}")
    tests = []
    for amp in (1.0, 2.0, 4.0, 8.0):
        up_i = up.copy()
        up_i[inj:] -= amp * shape[None]
        r_i, _ = analyse(up_i, t, d)
        lbf = np.nanmedian(r_i["step_lbf"][dome])
        tau = int(np.nanmedian(r_i["step_tau"][dome]))
        size = np.nanmedian(r_i["step_size"][dome])
        found = lbf > THRESH and abs(tau - inj) <= 1
        fa = np.nanmean((r_i["step_lbf"] > THRESH)[far])
        tests.append((amp, lbf, tau, size, found, fa))
        print(f"  {amp:.0f} cm: dome log BF {lbf:+.1f}, date {d['dates'][tau]}, estimated step {size:+.1f} cm -> "
              f"{'FOUND' if found else 'not found'}; far-field false alarms {fa:.1%}")
    np.savez_compressed(ML / "changepoint.npz", **{k: v.astype(np.float32) for k, v in res.items()},
                        tests=np.array([x[:4] + (float(x[4]), x[5]) for x in tests]), inj=inj)
    plot(res, tests, d, inj)


def plot(res, tests, d, inj):
    fig, ax = plt.subplots(1, 4, figsize=(18, 4.4))
    for a, kind, unit in ((ax[0], "step", "cm"), (ax[1], "hinge", "cm/yr")):
        v = np.where(res[f"{kind}_lbf"] > THRESH, res[f"{kind}_size"], np.nan)
        lim = 6 if kind == "step" else 30
        im = a.imshow(v, cmap="RdBu_r", vmin=-lim, vmax=lim)
        a.set_title(f"significant {kind} (BF > 20): size, {unit}", fontsize=9)
        fig.colorbar(im, ax=a); a.set_xticks([]); a.set_yticks([])
    m = res["hinge_lbf"] > THRESH
    if m.any():
        ax[2].hist(res["hinge_tau"][m].ravel(), bins=np.arange(-.5, len(d["dates"]) + .5), color="tab:purple")
    ax[2].set_title("dates of significant rate changes", fontsize=9)
    ax[2].set_xticks(range(0, len(d["dates"]), 4)); ax[2].set_xticklabels([str(x)[:10] for x in d["dates"][::4]], rotation=40, fontsize=7)
    amps = [x[0] for x in tests]
    ax[3].plot(amps, [x[1] for x in tests], "o-", label="dome log Bayes factor")
    ax[3].plot(amps, [-x[3] for x in tests], "s-", label="estimated deflation, cm")
    ax[3].plot(amps, amps, "k:", lw=.8, label="injected deflation")
    ax[3].axhline(THRESH, ls="--", color="tab:blue", lw=.7)
    ax[3].set_xscale("log", base=2); ax[3].set_xlabel("synthetic deflation at the dome centre, cm")
    ax[3].set_title(f"injection test at {d['dates'][inj]}", fontsize=9); ax[3].legend(fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(ML / "changepoint.png", dpi=110, facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    main()
