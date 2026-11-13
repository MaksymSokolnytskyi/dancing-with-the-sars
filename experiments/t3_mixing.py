"""T3: radiosonde mixed-layer depth (parcel method) vs radar noise (pre-registered).

    uv run python scripts/t3_mixing.py
"""
from __future__ import annotations

import json
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

sys.path[:0] = [str(Path(__file__).parent), str(Path(__file__).parents[1] / "scripts")]
from t2_radiosonde import read  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
F = ROOT / "data/raw/igra/por/MXM00076679-data.txt"


def read_p():
    """Same as t2.read but keeps pressure (Pa) for potential temperature."""
    out, cur, rows = {}, None, []
    with open(F) as f:
        for line in f:
            if line.startswith("#"):
                if cur is not None and rows:
                    out[cur] = np.array(rows)
                y, mo, d, h = int(line[13:17]), int(line[18:20]), int(line[21:23]), int(line[24:26])
                cur = datetime(y, mo, d, h) if y >= 2025 and h != 99 and datetime(y, mo, d) >= datetime(2025, 10, 1) else None
                rows = []
            elif cur is not None:
                p, gph, t = int(line[9:15]), int(line[16:21]), int(line[22:27])
                if p > 0 and gph > -8000 and t > -8000:
                    rows.append((gph, p, t / 10))
    if cur is not None and rows:
        out[cur] = np.array(rows)
    return out


def zmix(prof, dtheta=1.0):
    z, p, t = prof[:, 0].astype(float), prof[:, 1].astype(float), prof[:, 2]
    o = np.argsort(z)
    z, p, t = z[o], p[o], t[o]
    th = (t + 273.15) * (1e5 / p) ** 0.2857
    z0 = z[0]
    above = np.nonzero(th > th[0] + dtheta)[0]
    if len(above) == 0 or z[-1] - z0 < 4000:
        return np.nan
    i = above[0]
    if i == 0:
        return 0.0
    # interpolate the crossing
    zc = np.interp(th[0] + dtheta, [th[i - 1], th[i]], [z[i - 1], z[i]])
    return float(zc - z0) / 1000


def main():
    snd = read_p()
    res = {}
    for name, series, hour in [("evening", "mexico-city", 0), ("morning", "mexico-city-asc", 12)]:
        z = np.load(ROOT / "data/processed" / series / "nisar_series.npz")
        dates = [datetime.fromisoformat(str(d)) for d in z["dates"]]
        rows = [(f"{d:%Y-%m-%d}", zmix(snd[k])) for d in dates if (k := datetime(d.year, d.month, d.day, hour)) in snd]
        v = np.array([r[1] for r in rows])
        res[name] = dict(rows=rows, median=float(np.nanmedian(v)))
        print(f"{name}: n {np.isfinite(v).sum()}, z_mix median {np.nanmedian(v):.2f} km "
              f"(16-84%: {np.nanpercentile(v, 16):.2f}-{np.nanpercentile(v, 84):.2f})")
    ratio = res["evening"]["median"] / max(res["morning"]["median"], 1e-3)
    print(f"P1: evening/morning ratio {ratio:.1f} (need >= 3)")
    print(f"P2: evening z_mix {res['evening']['median']:.2f} km in 1-4: {1 <= res['evening']['median'] <= 4}")
    g = json.load(open(ROOT / "data/processed/mexico-city/g1_date_noise.json"))["DESCENDING"]
    sig = dict(zip(g["dates"], np.sqrt(np.clip(g["s2"], 0, None))))
    pr = [(sig[d], h) for d, h in res["evening"]["rows"] if d in sig and np.isfinite(h)]
    rho, p = spearmanr(*zip(*pr))
    print(f"P3: evening sigma_date vs z_mix: rho {rho:+.2f} (p {p:.2f}, n {len(pr)})")
    res["ratio"], res["p3"] = ratio, [float(rho), float(p), len(pr)]
    json.dump(res, open(ROOT / "data/processed/mexico-city/t3_mixing.json", "w"), indent=1)


if __name__ == "__main__":
    main()
