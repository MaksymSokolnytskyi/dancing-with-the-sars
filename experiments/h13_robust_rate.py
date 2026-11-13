"""H13 (pre-registered in notes/experiments.md): Theil-Sen rate within chain segments vs the chain sum.

    uv run --group ml python scripts/h13_robust_rate.py houston desc
Same pairs, iono 'plane', same GNSS evaluation as H12.
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np
import torch
from rasterio.warp import transform as wt

sys.path[:0] = [str(Path(__file__).parent), str(Path(__file__).parents[1] / "scripts")]
from trust_audit import SITES, plane_removed  # noqa: E402

from hachure import nisar

ROOT = Path(__file__).resolve().parents[1]
MIN_DT = 60.0
dev = torch.device("cuda" if torch.cuda.is_available() else "cpu")


def theil_sen_segments(disps: np.ndarray, days: np.ndarray, chunk: int = 40000) -> np.ndarray:
    """disps (n, H, W) in m per sequential pair; returns rate m/yr (H, W): median within-segment slope, dt >= MIN_DT."""
    n, H, W = disps.shape
    t = np.concatenate([[0.0], np.cumsum(days)])
    I, J = np.triu_indices(n + 1, 1)
    keep = (t[J] - t[I]) >= MIN_DT
    I, J = I[keep], J[keep]
    dt = torch.tensor(t[J] - t[I], dtype=torch.float32, device=dev)[:, None]
    D = disps.reshape(n, -1)
    out = np.full(D.shape[1], np.nan, np.float32)
    for s in range(0, D.shape[1], chunk):
        d = torch.tensor(D[:, s:s + chunk], dtype=torch.float32, device=dev)
        ok = torch.isfinite(d)
        cum = torch.zeros(n + 1, d.shape[1], device=dev)
        seg = torch.zeros(n + 1, d.shape[1], dtype=torch.int32, device=dev)
        for k in range(n):
            cum[k + 1] = torch.where(ok[k], cum[k] + torch.nan_to_num(d[k]), torch.zeros_like(cum[k]))
            seg[k + 1] = seg[k] + (~ok[k]).int()
        same = seg[I] == seg[J]
        sl = torch.where(same, (cum[J] - cum[I]) / dt, torch.full_like(dt.expand(-1, d.shape[1]), float("nan")))
        cnt = same.sum(0)
        med = torch.nanmedian(sl, 0).values
        med[cnt < 3] = float("nan")
        out[s:s + chunk] = (med * 365.25).cpu().numpy()
    return out.reshape(H, W)


def main():
    site, pas = sys.argv[1], sys.argv[2]
    s = SITES[site]
    bbox = s["bbox"]
    cache = ROOT / "data/raw/nisar_cache" / ("bbox_" + "_".join(f"{v:.3f}" for v in bbox) + f"_h{s['h']}")
    found = nisar.search((bbox[0] + bbox[2]) / 2, (bbox[1] + bbox[3]) / 2,
                         direction="ASCENDING" if pas == "asc" else "DESCENDING")
    chain = [g for g in nisar.sequential_chain(nisar.dominant_track(found)) if (cache / f"{g.name}.npz").exists()]
    pairs = []
    for g in chain:
        try:
            pairs.append(nisar.process(nisar.RawWindow.load(cache / f"{g.name}.npz"), ref_bbox=s["ref"], coh_min=0.3,
                                       iono="plane"))
        except ValueError:
            pass
    rate, coh, frac, x, y, disps = nisar.stack_rate(pairs)
    days = np.array([p.days for p in pairs], float)
    ts = theil_sen_segments(disps, days)
    ts[~np.isfinite(rate)] = np.nan
    XX, YY = np.meshgrid(x, y)
    if site == "mexico-city":
        st = [dict(station="MMX1", lon=-99.0684, lat=19.4317, rate_cm_yr=-26.7)]
    else:
        st = [r for r in csv.DictReader(open(ROOT / f"data/processed/gnss/{site}_rates.csv"))
              if r["window"] == pas and r["used"] == "True"]
    xs, ys = wt("EPSG:4326", f"EPSG:{pairs[0].epsg}", [float(r["lon"]) for r in st], [float(r["lat"]) for r in st])
    for label, rr in (("chain sum", rate), ("Theil-Sen seg", ts)):
        v = plane_removed(rr * 100)
        d = []
        for r, sx, sy in zip(st, xs, ys):
            m = np.hypot(XX - sx, YY - sy) <= 250
            d.append(np.nanmedian(v[m]) - float(r["rate_cm_yr"]) if np.isfinite(v[m]).sum() >= 3 else np.nan)
        d = np.array(d)
        if site == "mexico-city":
            print(f"{site} {pas} {label:14s}: MMX1 NISAR {d[0] - 26.7:+.1f} vs GNSS -26.7 -> error {abs(d[0]):.2f}")
        else:
            d = d[np.isfinite(d)]
            d -= np.median(d)
            print(f"{site} {pas} {label:14s}: {len(d)} stations, P1 RMS {np.sqrt(np.mean(d ** 2)):.2f}, "
                  f"median|d| {np.median(np.abs(d)):.2f}, max|d| {np.abs(d).max():.1f}")
    np.save(ROOT / f"data/processed/trust/h13_{site}_{pas}_theilsen_rate.npy", ts)


if __name__ == "__main__":
    main()
