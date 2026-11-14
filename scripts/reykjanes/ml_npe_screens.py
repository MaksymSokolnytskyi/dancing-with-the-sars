"""Neural posterior estimation with REAL atmospheric screens as the simulator noise, compared
with the Gaussian noise model used so far.

Noise variants added to the same noise-free simulations (source + nuisance offset/plane):
  gauss   exponential covariance fitted to least-squares residuals (as ml_npe_source.py)
  raw     rate noise built from two real per-date screens of the track, (s_j - s_i) / span,
          rotated about Svartsengi by a random angle, possibly mirrored, random sign
          (ml_screens.py). Without the rotation the ~20 date screens span a ~20-dimensional
          subspace that the network learns to project out, and its posteriors are far too narrow
  gp      the same with the GP-smoothed screens + the GP white nugget
Training uses screens from 75 % of the dates; the held-out dates (always including the
first and last date, whose difference is the real data's own noise) build the test noise,
so calibration is measured on atmospheres the network has never seen.

Cross test: every network is also calibrated on held-out simulations with RAW screens,
which says how well a network trained on idealised noise survives the real atmosphere.

Usage: python scripts/reykjanes/ml_npe_screens.py {mogi|ellipsoid} [n_sims] [--epochs 25]
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
import torch
from scipy.spatial import cKDTree

import reykjanes_mogi as M
import ml_npe_source as NS
import ml_npe_ellipsoid as NE
from reykjanes_joint import SIGMA_GNSS
from reykjanes_sources import Data

ML = Path(__file__).resolve().parents[2] / "output" / "reykjanes" / "ml"
RNG = np.random.default_rng(21)
N_ROT = 72   # rotation angles (5 deg); x2 with mirror images
torch.manual_seed(21)


def split_dates(T, rng):
    inner = np.arange(1, T - 1)
    test = set(rng.choice(inner, max(2, int(round(0.25 * T)) - 2), replace=False).tolist()) | {0, T - 1}
    train = [(i, j) for i in range(T) for j in range(i + 1, T) if i not in test and j not in test]
    held = [(i, j) for i in range(T) for j in range(i + 1, T) if (i in test or j in test) and (i, j) != (0, T - 1)]
    return np.array(train), np.array(held)


class Noise:
    def __init__(self, O, gauss_noise):
        z = np.load(ML / "screens.npz")
        self.O, self.g = O, gauss_noise
        self.lib = {}
        for t in M.TRACKS:
            T = z[f"raw_{t}"].shape[0]
            tr, te = split_dates(T, RNG)
            k = 365.25 / float(z[f"span_{t}"])
            e, n = O.tr[t]["e"], O.tr[t]["n"]
            tree = cKDTree(np.c_[e, n])
            maps = []
            for f in (1.0, -1.0):
                for a in np.deg2rad(np.arange(N_ROT) * 360 / N_ROT):
                    er, nr = np.cos(a) * f * e - np.sin(a) * n, np.sin(a) * f * e + np.cos(a) * n
                    maps.append(tree.query(np.c_[er, nr])[1])
            self.lib[t] = dict(raw=z[f"raw_{t}"] * k, gp=z[f"gp_{t}"] * k, nug=float(z[f"nugget_{t}"]) * k,
                               train=tr, test=te, maps=np.array(maps))
            print(f"{t}: {T} dates, {len(tr)} training / {len(te)} held-out screen pairs")

    def add(self, x, kind, part):
        x = x.copy()
        N, i0 = len(x), 0
        for t in M.TRACKS:
            n = self.O.tr[t]["y"].size
            if kind == "gauss":
                x[:, i0:i0 + n] += (RNG.standard_normal((N, n)) @ self.g[t]["chol"].T).astype(np.float32)
            else:
                L = self.lib[t]
                S = L[kind]
                # one real screen difference, rotated about Svartsengi by a random angle and
                # possibly mirrored: ~20 date screens span only a ~20-dimensional subspace,
                # which a network learns to project out (over-confident posteriors)
                pr = L[part][RNG.integers(0, len(L[part]), N)]
                mp = L["maps"][RNG.integers(0, len(L["maps"]), N)]
                sign = RNG.choice([-1.0, 1.0], N)[:, None]
                d = S[pr[:, 1]] - S[pr[:, 0]]
                nz = sign * np.take_along_axis(d, mp, 1)
                if kind == "gp":
                    nz = nz + RNG.normal(0, np.sqrt(2) * L["nug"], nz.shape)
                x[:, i0:i0 + n] += nz.astype(np.float32)
            i0 += n
        x[:, i0:] += RNG.normal(0, SIGMA_GNSS, (N, x.shape[1] - i0)).astype(np.float32)
        return x


def clean_sims(kind, O, N):
    if kind == "mogi":
        q, p, nuis = NS.sample_prior(N)
        return q.astype(np.float32), O.forward(p, nuis).astype(np.float32)
    q = NE.LO + RNG.random((N, 7)) * (NE.HI - NE.LO)
    nuis = np.c_[RNG.uniform(-10, 10, N), RNG.normal(0, 2, N), RNG.normal(0, 2, N),
                 RNG.uniform(-10, 10, N), RNG.normal(0, 2, N), RNG.normal(0, 2, N)]
    x = np.empty((N, sum(O.sizes) + 3 * len(O.gn)), np.float32)
    t0 = time.time()
    for i in range(N):
        x[i] = O.forward_one(q[i], nuis[i])
    print(f"  {N} elongated-body simulations in {time.time() - t0:.0f} s")
    return q.astype(np.float32), x


def fit(kind, x, q, cols, label, epochs):
    if kind == "mogi":
        return NS.train(x, q, cols, label, epochs)
    return NE.train(x, NE.to_z(q), cols, label, epochs)


def sample(kind, m, x, n):
    if kind == "mogi":
        return NS.posterior(m, x, n)                      # (n, B, 4), log10 V in km^3
    return NE.post(m, x, n)                               # (n, B, 7)


def coverage(kind, m, xt, qt):
    s = sample(kind, m, xt, 2000)
    return {lev: ((qt >= np.quantile(s, (1 - lev) / 2, 0)) & (qt <= np.quantile(s, (1 + lev) / 2, 0))).mean(0).round(3).tolist()
            for lev in (.68, .95)}


def physical(kind, v):
    v = v.copy()
    j = 3 if kind == "mogi" else 5
    v[:, j] = 10 ** v[:, j] * 1e3
    return v


def main(kind, n_sims, epochs):
    D = Data()
    O = NS.Obs(D) if kind == "mogi" else NE.EllObs(D)
    _, resid = NS.lsq_fit(O)
    gnoise = NS.exp_cov_fit(O, resid)
    NZ = Noise(O, gnoise)
    names = ("east km", "north km", "depth km", "dV 1e6 m3/yr") if kind == "mogi" else \
        ("east km", "north km", "depth km", "half-length km", "strike deg", "V 1e6 m3/yr", "sill share w")
    q, x0 = clean_sims(kind, O, n_sims)
    qt, xt0 = clean_sims(kind, O, 800)
    test_raw = NZ.add(xt0, "raw", "test")
    xobs = np.concatenate([np.concatenate([O.tr[t]["y"] for t in M.TRACKS]), O.gv])[None].astype(np.float32)
    nin = sum(O.tr[t]["y"].size for t in M.TRACKS)
    res = {}
    for nk in ("gauss", "raw", "gp"):
        x = NZ.add(x0, nk, "train")
        test_own = NZ.add(xt0, nk, "test")
        for label, cols in (("NISAR + GNSS", np.arange(x.shape[1])), ("NISAR only", np.arange(nin))):
            tag = f"{nk} | {label}"
            m = fit(kind, x, q, cols, tag, epochs)
            cov_own = coverage(kind, m, test_own, qt)
            cov_raw = coverage(kind, m, test_raw, qt)
            s = physical(kind, sample(kind, m, xobs, 20000)[:, 0])
            summ = {nm: np.round(np.percentile(s[:, i], [2.5, 16, 50, 84, 97.5]), 3).tolist() for i, nm in enumerate(names)}
            row = dict(summary=summ, cal_own=cov_own, cal_on_raw=cov_raw, samples=s)
            if kind == "ellipsoid":
                row["p_long"] = round(float((s[:, 3] > 1.5).mean()), 3)
            res[tag] = row
            print(f"[{kind} {tag}] calibration on its own noise 68 %: {cov_own[.68]}; on held-out REAL screens 68 %: "
                  f"{cov_raw[.68]}, 95 %: {cov_raw[.95]}")
            print("   " + "; ".join(f"{nm} {v[2]:.2f} [{v[1]:.2f}..{v[3]:.2f}]" for nm, v in summ.items())
                  + (f"; P(>3 km) {row['p_long']:.0%}" if kind == "ellipsoid" else ""), flush=True)
    np.savez(ML / f"npe_screens_{kind}.npz", **{k.replace(" | ", "__").replace(" ", "_").replace("+", "and"): v["samples"][:6000]
                                                for k, v in res.items()})
    (ML / f"npe_screens_{kind}.json").write_text(json.dumps({k: {kk: vv for kk, vv in v.items() if kk != "samples"}
                                                             for k, v in res.items()}, indent=1), encoding="utf-8")


if __name__ == "__main__":
    kind = sys.argv[1]
    n = int(sys.argv[2]) if len(sys.argv) > 2 and not sys.argv[2].startswith("--") else (80000 if kind == "mogi" else 40000)
    ep = int(sys.argv[sys.argv.index("--epochs") + 1]) if "--epochs" in sys.argv else 25
    main(kind, n, ep)
