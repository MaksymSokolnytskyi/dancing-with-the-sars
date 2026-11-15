"""Self-check of the W polarimetry helpers on synthetic scattering vectors (no data, no network).

uv run python experiments/w_selfcheck.py      -> prints each check; exits 1 if one fails

- coherency(): the Pauli coherency matrix from the 4x4 covariance terms equals the one computed directly from k;
- h_alpha(): a pure surface (HH = VV, no cross-pol) gives H ~ 0, alpha ~ 0 deg; a random volume gives H ~ 1, alpha ~ 55-60 deg;
- faraday_deg(): a reflection-symmetric target rotated as M = R(Omega) S R(Omega) gives back +Omega (sign convention).
"""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src"), str(ROOT / "experiments")]
from w_tests import coherency, faraday_deg, h_alpha  # noqa: E402

N = 200_000
POLS = ("hh", "hv", "vh", "vv")


def cov(S):
    return {a + b: (S[a] * np.conj(S[b])).mean() for a in POLS for b in POLS}


def cn(rng, scale=1.0):
    return scale * (rng.normal(size=N) + 1j * rng.normal(size=N)) / np.sqrt(2)


def rotate(S, omega_deg):
    o = np.radians(omega_deg)
    R = np.array([[np.cos(o), np.sin(o)], [-np.sin(o), np.cos(o)]])
    M = np.einsum("ij,jkn,kl->iln", R, np.array([[S["hh"], S["hv"]], [S["vh"], S["vv"]]]), R)
    return {"hh": M[0, 0], "hv": M[0, 1], "vh": M[1, 0], "vv": M[1, 1]}


def main():
    rng = np.random.default_rng(0)
    checks = []
    hv = cn(rng, 0.5)
    vol = {"hh": cn(rng), "hv": hv, "vh": hv, "vv": cn(rng)}
    k = np.array([vol["hh"] + vol["vv"], vol["hh"] - vol["vv"], vol["hv"] + vol["vh"]]) / np.sqrt(2)
    err = float(np.abs(coherency(cov(vol)) - (k @ k.conj().T) / N).max())
    checks.append(("coherency equals direct Pauli T", err < 1e-9, f"max |diff| {err:.1e}"))

    s = cn(rng)
    zero = np.zeros(N, complex)
    H, a = h_alpha(coherency(cov({"hh": s, "hv": zero, "vh": zero, "vv": s})))
    checks.append(("surface: H < 0.05, alpha < 5 deg", H < 0.05 and a < 5, f"H {H:.3f}, alpha {a:.1f}"))
    H, a = h_alpha(coherency(cov(vol)))
    checks.append(("random volume: H > 0.9, 50 < alpha < 62 deg", H > 0.9 and 50 < a < 62, f"H {H:.3f}, alpha {a:.1f}"))

    for om in (0.0, 3.0, 5.0, -4.0):
        est = faraday_deg(cov(rotate(vol, om)))
        checks.append((f"Faraday {om:+.0f} deg recovered with its sign", abs(est - om) < 0.2, f"estimate {est:+.2f}"))

    bad = 0
    for name, ok, info in checks:
        print(f"  {'ok' if ok else '!!'}  {name}: {info}")
        bad += not ok
    print("all good" if not bad else f"{bad} failed")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
