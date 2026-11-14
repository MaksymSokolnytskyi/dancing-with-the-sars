"""Noise covariance of 12-day Sentinel-1 pairs from a QUIET pair (7/5 Oct -> 19/17 Oct 2023, before the unrest of 25 Oct).

The covariance used so far was fitted to the residuals of a simple dike model, so it also contains model error
(hypothesis H3: it is too large). Here the quiet pair is binned onto the same 1 km cells (no mask), an offset +
plane is removed per track, and an exponential + white covariance is fitted to its semivariogram.
The residual quiet maps are also saved as one realization of REAL noise for the tests.

Writes output/reykjanes/ml/s1_noise_quiet.json (read by s1_common_v2.noise_info) and s1_noise_quiet.npz.
"""
import json

import numpy as np
from scipy.optimize import least_squares

import s1_grindavik_dike as S
from s1_common_v2 import GeomV2, pair_vector, product_dir

ML = S.ML


def main():
    idx = json.loads((S.S1 / "pairs_index.json").read_text())["quiet_2023"]
    G = GeomV2(gnss=False, mask_km=0, noise="residual")
    dirs = {p["track"]: product_dir(p["track"], p["ref"], p["sec"]) for p in idx}
    y, filled = pair_vector(G, dirs)
    out, maps = {}, {}
    for k, t in enumerate(S.TRACKS):
        m = np.flatnonzero(G.trk[:G.n_ins] == k)
        e, n, v = G.e[m], G.n[m], y[m]
        A = np.c_[np.ones(len(m)), e / 10, n / 10]
        r = v - A @ np.linalg.lstsq(A, v, rcond=None)[0]
        D = np.hypot(e[:, None] - e[None], n[:, None] - n[None])
        Vg = 0.5 * (r[:, None] - r[None]) ** 2
        edges = np.arange(0, 16.5, 1.0)
        h = 0.5 * (edges[1:] + edges[:-1])
        g = np.array([Vg[(D > a) & (D <= b)].mean() for a, b in zip(edges[:-1], edges[1:])])
        f = lambda p: (np.exp(p[2]) + np.exp(p[0]) * (1 - np.exp(-h / np.exp(p[1])))) - g
        p = np.exp(least_squares(f, np.log([np.var(r), 3.0, 0.2 * np.var(r)])).x)
        out[t] = dict(sill_cm2=round(float(p[0]), 3), range_km=round(float(p[1]), 2), white_cm2=round(float(p[2]), 3),
                      rms_cm=round(float(r.std()), 2), semivariogram=dict(h_km=h.tolist(), gamma_cm2=g.round(3).tolist()))
        maps[t] = r
        print(f"{t}: quiet-pair rms {r.std():.2f} cm; covariance sill {p[0]:.2f} cm2, range {p[1]:.1f} km, white {p[2]:.2f} cm2 "
              f"(residual-fit covariance: sill {G.noise_info[t]['sill_cm2']:.1f}, range {G.noise_info[t]['range_km']:.1f}, white {G.noise_info[t]['white_cm2']:.1f})", flush=True)
    (ML / "s1_noise_quiet.json").write_text(json.dumps(dict(pairs=idx, filled_fraction=round(float(filled), 3), tracks=out,
                                                            residual_fit=G.noise_info), indent=1), encoding="utf-8")
    np.savez(ML / "s1_noise_quiet.npz", e=G.e[:G.n_ins], n=G.n[:G.n_ins], trk=G.trk[:G.n_ins], **maps)
    print(f"cells filled by interpolation: {filled:.1%}")


if __name__ == "__main__":
    main()
