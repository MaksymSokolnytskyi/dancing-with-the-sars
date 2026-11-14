"""Total noise covariance = atmosphere (quiet pair, s1_noise_quiet.json) + MODEL ERROR.

H3 showed that the atmosphere of a 12-day Sentinel-1 pair here is only ~0.3-0.4 cm rms, while the residuals of the
dike models are ~5 cm: the dominant error is the model (simple dike, no graben faults, no strike-slip, 1 km patches).
Training NPE on the atmosphere alone would make it overconfident on real data. The model-error part is fitted to the
semivariogram of the residuals of the best model (distributed opening, InSAR only, s1_dike_distributed_nognss_resid.npz)
minus the semivariogram of the quiet pair: second exponential (sill2, range2) + white2.

Writes output/reykjanes/ml/s1_noise_total.json (read first by s1_common_v2.noise_info).
"""
import json

import numpy as np
from scipy.optimize import least_squares

import s1_grindavik_dike as S

ML = S.ML


def semivariogram(e, n, r, edges):
    D = np.hypot(e[:, None] - e[None], n[:, None] - n[None])
    V = 0.5 * (r[:, None] - r[None]) ** 2
    return np.array([V[(D > a) & (D <= b)].mean() for a, b in zip(edges[:-1], edges[1:])])


def main():
    Q = json.loads((ML / "s1_noise_quiet.json").read_text())["tracks"]
    R = np.load(ML / "s1_dike_distributed_nognss_resid.npz")
    edges = np.arange(0, 16.5, 1.0)
    h = 0.5 * (edges[1:] + edges[:-1])
    out = {}
    for k, t in enumerate(S.TRACKS):
        m = R["trk"] == k
        e, n, r = R["e"][m], R["n"][m], R["resid"][m]
        A = np.c_[np.ones(m.sum()), e / 10, n / 10]
        r = r - A @ np.linalg.lstsq(A, r, rcond=None)[0]
        g = semivariogram(e, n, r, edges)
        q = Q[t]
        gq = q["white_cm2"] + q["sill_cm2"] * (1 - np.exp(-h / q["range_km"]))
        f = lambda p: (np.exp(p[2]) + np.exp(p[0]) * (1 - np.exp(-h / np.exp(p[1])))) - (g - gq)
        p = np.exp(least_squares(f, np.log([max(np.var(r) - q["sill_cm2"], 0.1), 3.0, 0.3])).x)
        out[t] = dict(sill_cm2=q["sill_cm2"], range_km=q["range_km"], white_cm2=q["white_cm2"],
                      sill2_cm2=round(float(p[0]), 3), range2_km=round(float(p[1]), 2), white2_cm2=round(float(p[2]), 3),
                      residual_rms_cm=round(float(r.std()), 2))
        print(f"{t}: residual rms {r.std():.2f} cm; atmosphere sill {q['sill_cm2']:.2f} cm2 / {q['range_km']:.1f} km; "
              f"model error sill {p[0]:.2f} cm2 / {p[1]:.1f} km + white {p[2]:.2f} cm2", flush=True)
    (ML / "s1_noise_total.json").write_text(json.dumps(dict(tracks=out, note=__doc__.split("\n")[0]), indent=1), encoding="utf-8")


if __name__ == "__main__":
    main()
