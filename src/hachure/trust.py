"""Trust audit for a NISAR 12-day chain: how much to believe a rate map, from the radar itself.

Pieces (each validated on Mexico City, see notes/experiments.md: G1', V1, P-trust):
  date_noise   per-date (node) noise from the covariance of adjacent pairs; the pair network is a path graph,
               adjacent edges share one date, so cov(r_k, r_k+1) = -sigma^2(date k+1)          (G1')
  edge_share   how much of a pair's variance is NOT explained by its two dates: 1 - (s_k^2 + s_k+1^2) / var(r_k)
               edge noise (processing, per-pair iono estimate) random-walks into rates; node noise telescopes
  sigma_map    per-pixel 1-sigma of the chain-sum rate: Q = sum r^2 + 2 sum r_k r_k+1, smoothed 1 km   (P-trust)
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.ndimage import gaussian_filter


def smooth(a, sigma):
    ok = np.isfinite(a)
    num = gaussian_filter(np.where(ok, a, 0), sigma)
    den = gaussian_filter(ok.astype(float), sigma)
    with np.errstate(all="ignore"):
        out = num / den
    out[den < 0.3] = np.nan
    return out


def residuals(disps: np.ndarray, days: np.ndarray, skip: tuple[int, ...]) -> tuple[np.ndarray, np.ndarray]:
    """Leave-out rate (excluding pairs in `skip`) per pixel and the valid mask of all pairs."""
    ok = np.isfinite(disps)
    tot = np.where(ok, disps, 0).sum(0)
    span = np.where(ok, days[:, None, None], 0).sum(0)
    for k in skip:
        tot = tot - np.where(ok[k], disps[k], 0)
        span = span - np.where(ok[k], days[k], 0)
    with np.errstate(all="ignore"):
        return tot / span, ok


@dataclass
class Audit:
    dates: list[str]
    node_var: np.ndarray            # cm^2, interior dates
    pair_var: np.ndarray            # cm^2, all pairs (median-removed residual variance)
    edge_share: float               # median over pairs with both ends interior
    sigma: np.ndarray = field(repr=False)    # per-pixel 1-sigma of the rate, cm/yr
    rate: np.ndarray = field(repr=False)     # chain-sum rate, cm/yr


def date_noise(disps: np.ndarray, days: np.ndarray):
    """sigma^2 of each interior date (cm^2) and the residual variance of each pair (cm^2)."""
    n = len(disps)
    node, pv = [], []
    for k in range(n - 1):
        v, ok = residuals(disps, days, (k, k + 1))
        m = ok[k] & ok[k + 1] & np.isfinite(v)
        a = (disps[k] - v * days[k])[m] * 100
        b = (disps[k + 1] - v * days[k + 1])[m] * 100
        a, b = a - np.median(a), b - np.median(b)
        node.append(-np.mean((a - a.mean()) * (b - b.mean())))
    for k in range(n):
        v, ok = residuals(disps, days, (k,))
        m = ok[k] & np.isfinite(v)
        a = (disps[k] - v * days[k])[m] * 100
        pv.append(np.var(a - np.median(a)))
    return np.array(node), np.array(pv)


def edge_share(node: np.ndarray, pair_var: np.ndarray) -> float:
    """Pairs k = 1..n-2 have both dates interior: share of their variance not explained by node noise."""
    s = [1 - (max(node[k - 1], 0) + max(node[k], 0)) / pair_var[k] for k in range(1, len(node)) if pair_var[k] > 0]
    return float(np.median(s)) if s else float("nan")


def sigma_map(disps: np.ndarray, days: np.ndarray, rate_valid: np.ndarray, dx_m: float) -> np.ndarray:
    T = days.sum() / 365.25
    r = []
    for i in range(len(disps)):
        v, ok = residuals(disps, days, (i,))
        r.append(np.where(ok[i], disps[i] - v * days[i], np.nan) * 100)
    r = np.array(r)
    q = np.nansum(r ** 2, 0) + 2 * np.nansum(r[:-1] * r[1:], 0)
    q[~rate_valid] = np.nan
    qs = smooth(q, 1000 / dx_m / 2.355)
    sig = np.sqrt(np.clip(qs, 0.01, None)) / T
    sig[~rate_valid] = np.nan
    return sig


def segment_var(disps: np.ndarray, node: np.ndarray) -> np.ndarray:
    """M2b: extra node variance (cm^2) from chain gaps. node[j-1] = variance of interior date j (j = 1..n-1)."""
    ok = np.isfinite(disps)
    nv = np.clip(np.nan_to_num(node), 0, None)
    extra = np.zeros(disps.shape[1:])
    n = len(disps)
    for k in range(1, n):                      # segment starts at pair k (date k) after a missing pair k-1
        extra += np.where(ok[k] & ~ok[k - 1], nv[k - 1], 0)
    for k in range(n - 1):                     # segment ends at pair k (date k+1) before a missing pair k+1
        extra += np.where(ok[k] & ~ok[k + 1], nv[k], 0)
    return extra


def audit(pairs, dates: list[str], segments: bool = False) -> Audit:
    """pairs: list of nisar.PairWindow of one sequential chain (same track)."""
    from hachure import nisar
    rate, coh, frac, x, y, disps = nisar.stack_rate(pairs)
    days = np.array([p.days for p in pairs], float)
    node, pv = date_noise(disps, days)
    sig = sigma_map(disps, days, np.isfinite(rate), abs(x[1] - x[0]))
    if segments:                               # P5b: add the gap term of M2b
        T = days.sum() / 365.25
        sig = np.sqrt(sig ** 2 + segment_var(disps, node) / T ** 2)
    return Audit(dates=dates[1:len(pairs)], node_var=node, pair_var=pv, edge_share=edge_share(node, pv),
                 sigma=sig, rate=rate * 100)


def report(name: str, a: Audit, gnss_rows: list[dict] | None = None) -> str:
    s = np.sqrt(np.clip(a.node_var, 0, None))
    worst = np.argsort(-a.node_var)[:3]
    lines = [f"# Trust audit: {name}", "",
             f"- Dates: {len(a.dates) + 2} ({len(a.dates)} interior). Median per-date noise **{np.median(s):.2f} cm**; "
             f"negative estimates (noise below detection) on {np.mean(a.node_var <= 0) * 100:.0f}% of dates.",
             f"- Noisiest dates: " + ", ".join(f"{a.dates[i]} ({s[i]:.1f} cm)" for i in worst) + ".",
             f"- **Edge share** (pair noise not explained by its two dates; it accumulates along the chain): "
             f"**{a.edge_share * 100:.0f}%**.",
             f"- Per-pixel rate 1σ: median **{np.nanmedian(a.sigma):.2f} cm/yr**, 90th percentile "
             f"{np.nanpercentile(a.sigma, 90):.2f}.", ""]
    if gnss_rows:
        lines += ["| Station | GNSS cm/yr | NISAR − offset | σ NISAR | z |", "|---|---|---|---|---|"]
        for r in gnss_rows:
            lines.append(f"| {r['station']} | {r['gnss']:+.1f} ± {r['sigma_gnss']:.1f} | {r['nisar']:+.1f} | "
                         f"{r['sigma_nisar']:.1f} | {r['z']:+.2f} |")
        z = np.array([r["z"] for r in gnss_rows])
        lines += ["", f"- |z| < 1 for **{np.mean(np.abs(z) < 1) * 100:.0f}%** of stations (ideal 68%), "
                      f"|z| < 2 for {np.mean(np.abs(z) < 2) * 100:.0f}% (ideal 95%)."]
    return "\n".join(lines) + "\n"
