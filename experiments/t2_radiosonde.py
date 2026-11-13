"""T2: radiosonde wet-delay scale height at NISAR overpass times (pre-registered).

    uv run python scripts/t2_radiosonde.py
"""
from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
F = ROOT / "data/raw/igra/por/MXM00076679-data.txt"


def read(start=datetime(2025, 10, 1)):
    out, cur, rows = {}, None, []
    with open(F) as f:
        for line in f:
            if line.startswith("#"):
                if cur is not None and rows:
                    out[cur] = np.array(rows)
                y, mo, d, h = int(line[13:17]), int(line[18:20]), int(line[21:23]), int(line[24:26])
                cur = datetime(y, mo, d, h) if y >= start.year and h != 99 and datetime(y, mo, d) >= start else None
                rows = []
            elif cur is not None:
                gph, t, dpd = int(line[16:21]), int(line[22:27]), int(line[34:39])
                if gph > -8000 and t > -8000 and dpd > -8000:
                    rows.append((gph, t / 10, dpd / 10))
    if cur is not None and rows:
        out[cur] = np.array(rows)
    return out


def h50(prof):
    z, t, dpd = prof[:, 0].astype(float), prof[:, 1], prof[:, 2]
    o = np.argsort(z)
    z, t, dpd = z[o], t[o], dpd[o]
    td = t - dpd
    e = 6.112 * np.exp(17.67 * td / (td + 243.5))
    nw = 3.73e5 * e / (t + 273.15) ** 2
    z0 = z[0]
    m = z <= z0 + 12000
    z, nw = z[m], nw[m]
    if len(z) < 6 or z[-1] - z0 < 6000:
        return np.nan, np.nan
    cum = np.concatenate([[0], np.cumsum(0.5 * (nw[1:] + nw[:-1]) * np.diff(z))]) * 1e-6     # m
    return float(np.interp(0.5 * cum[-1], cum, z) - z0) / 1000, float(cum[-1] * 1000)       # km, mm


def main():
    snd = read()
    print(f"{len(snd)} soundings since 2025-10")
    res = {}
    for name, series, hour in [("evening", "mexico-city", 0), ("morning", "mexico-city-asc", 12)]:
        z = np.load(ROOT / "data/processed" / series / "nisar_series.npz")
        dates = [datetime.fromisoformat(str(d)) for d in z["dates"]]
        rows = []
        for d in dates:
            k = datetime(d.year, d.month, d.day, hour)
            if k in snd:
                h, zwd = h50(snd[k])
                rows.append((d, h, zwd))
        hs = np.array([r[1] for r in rows])
        res[name] = dict(n=int(np.isfinite(hs).sum()), h50_median=float(np.nanmedian(hs)),
                         zwd_median=float(np.nanmedian([r[2] for r in rows])),
                         rows=[(f"{r[0]:%Y-%m-%d}", r[1], r[2]) for r in rows])
        print(f"{name}: {res[name]['n']} soundings matched of {len(dates)} dates; median h50 {np.nanmedian(hs):.2f} km "
              f"(16-84%: {np.nanpercentile(hs, 16):.2f}-{np.nanpercentile(hs, 84):.2f}), median ZWD {res[name]['zwd_median']:.0f} mm")
    ratio = res["evening"]["h50_median"] / res["morning"]["h50_median"]
    print(f"P1: evening/morning h50 ratio {ratio:.2f} (need >= 1.5)")
    print(f"P2: evening median h50 {res['evening']['h50_median']:.2f} km in 1-4: {1 <= res['evening']['h50_median'] <= 4}")
    g = json.load(open(ROOT / "data/processed/mexico-city/g1_date_noise.json"))["DESCENDING"]
    sig = dict(zip(g["dates"], np.sqrt(np.clip(g["s2"], 0, None))))
    pairs = [(sig[d], h, zwd) for d, h, zwd in res["evening"]["rows"] if d in sig and np.isfinite(h)]
    rho, p = spearmanr([a for a, _, _ in pairs], [b for _, b, _ in pairs])
    rz, pz = spearmanr([a for a, _, _ in pairs], [c for _, _, c in pairs])
    print(f"P3: evening sigma_date vs h50: rho {rho:+.2f} (p {p:.2f}, n {len(pairs)})   [also vs ZWD: rho {rz:+.2f}, p {pz:.2f}]")
    res["ratio"], res["p3"] = ratio, [float(rho), float(p), len(pairs)]
    res["p3_zwd"] = [float(rz), float(pz)]
    json.dump(res, open(ROOT / "data/processed/mexico-city/t2_radiosonde.json", "w"), indent=1, default=str)


if __name__ == "__main__":
    main()
