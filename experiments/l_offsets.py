"""L-offsets (pre-registered in notes/experiments.md): GUNW pixel offsets over the Chhochen Khola source zone.

    uv run python scripts/l_offsets.py
Needs: scripts/fetch_offsets.py --name langtang ..., and scripts/l_retro.py (writes the zone masks).
"""
from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path

import numpy as np

from hachure import nisar

ROOT = Path(__file__).resolve().parents[1]
OFF = ROOT / "data/raw/nisar_offsets/langtang"
Z = ROOT / "data/processed/landslide/langtang_zones.npz"
EVENT = datetime(2026, 8, 26)
BASE = (datetime(2025, 12, 10), datetime(2026, 6, 8))


def rstd(v):
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    return float(1.4826 * np.median(np.abs(v - np.median(v)))) if len(v) >= 5 else float("nan")


def on_grid(d, key, x0, y0):
    """Place an offsets layer onto the zone grid by matching pixel-centre coordinates (to 1 m)."""
    xi = np.round((d["x"] - x0[0]) / (x0[1] - x0[0])).astype(int)
    yi = np.round((d["y"] - y0[0]) / (y0[1] - y0[0])).astype(int)
    if np.abs(x0[np.clip(xi, 0, len(x0) - 1)] - d["x"]).max() > 1 or np.abs(y0[np.clip(yi, 0, len(y0) - 1)] - d["y"]).max() > 1:
        return None
    out = np.full((len(y0), len(x0)), np.nan)
    kx, ky = (xi >= 0) & (xi < len(x0)), (yi >= 0) & (yi < len(y0))
    out[np.ix_(yi[ky], xi[kx])] = d[key][np.ix_(ky, kx)].astype(float)
    return out


def main():
    zz = np.load(Z)
    src, ctrl, x0, y0, z, north = zz["src"], zz["ctrl"], zz["x"], zz["y"], zz["z"], zz["north"]
    XX, YY = np.meshgrid(x0, y0)
    sy, sx = np.argwhere(src).mean(0)
    cxm, cym = np.interp(sx, np.arange(len(x0)), x0), np.interp(sy, np.arange(len(y0)), y0)
    disc = (np.hypot(XX - cxm, YY - cym) <= 1000) & north & (z >= 4500) & (z <= 6000)
    rows = []
    for f in sorted(OFF.glob("*.npz")):
        r, s = nisar._parse_dates(f.name)
        d = np.load(f)
        rg, az, pk, mk = (on_grid(d, k, x0, y0) for k in ("slantRangeOffset", "alongTrackOffset", "correlationSurfacePeak", "mask"))
        if rg is None:
            print(f"grid mismatch, skip {f.name[:50]}")
            continue
        mi = np.nan_to_num(mk).astype(int)
        ok = np.isfinite(rg) & np.isfinite(az) & (mi % 10 > 0) & ((mi // 10) % 10 > 0) & (mi // 100 == 0)
        thr = np.nanmedian(pk[ctrl & ok]) if (ctrl & ok).any() else np.inf
        good = ok & (pk >= thr)
        row = dict(ref=f"{r:%Y-%m-%d}", sec=f"{s:%Y-%m-%d}", days=float(d["days"]),
                   role="co" if r < EVENT < s else "base" if (r >= BASE[0] and s <= BASE[1]) else "pre" if s < EVENT else "post")
        for name, zone in (("src", src), ("disc", disc)):
            a, c = zone & good, ctrl & good
            if a.sum() >= 3 and c.sum() >= 20:
                row[f"{name}_rg"] = float(np.median(rg[a]) - np.median(rg[c]))
                row[f"{name}_az"] = float(np.median(az[a]) - np.median(az[c]))
            else:
                row[f"{name}_rg"] = row[f"{name}_az"] = float("nan")
            row[f"{name}_n"] = int(a.sum())
        row["ctrl_rg_std"] = rstd(rg[ctrl & good])
        rows.append(row)
        print(f"  {row['ref']} -> {row['sec']} {row['role']:4s} src rg {row['src_rg']:+6.3f} az {row['src_az']:+6.3f} "
              f"(n {row['src_n']}) | disc rg {row['disc_rg']:+6.3f} az {row['disc_az']:+6.3f} (n {row['disc_n']}) | "
              f"ctrl pixel scatter rg {row['ctrl_rg_std']:.3f}")
    base = [r for r in rows if r["role"] == "base"]
    pre = [r for r in rows if r["role"] in ("pre", "base")]
    last = sorted(pre, key=lambda r: r["sec"])[-1]
    res = dict(rows=rows, last_pre=last["ref"] + "->" + last["sec"])
    for name in ("src", "disc"):
        s_rg, s_az = rstd([r[f"{name}_rg"] for r in base]), rstd([r[f"{name}_az"] for r in base])
        dmag = [np.hypot(r[f"{name}_rg"], r[f"{name}_az"]) for r in base]
        dl = float(np.hypot(last[f"{name}_rg"], last[f"{name}_az"]))
        o1 = bool(abs(last[f"{name}_rg"]) >= 3 * s_rg or abs(last[f"{name}_az"]) >= 3 * s_az)
        o2 = bool(dl >= 3 * np.nanmedian(dmag))
        res[name] = dict(sigma_rg=s_rg, sigma_az=s_az, last_rg=last[f"{name}_rg"], last_az=last[f"{name}_az"],
                         last_mag=dl, base_median_mag=float(np.nanmedian(dmag)), O1=o1, O2=o2)
        tag = "JUDGED" if name == "src" else "secondary"
        print(f"{name} ({tag}): sigma rg {s_rg:.3f} az {s_az:.3f} m; last pre rg {last[f'{name}_rg']:+.3f} az "
              f"{last[f'{name}_az']:+.3f}; O1 {'PASS' if o1 else 'FAIL'}; |d| {dl:.3f} vs 3 x {np.nanmedian(dmag):.3f} "
              f"-> O2 {'PASS' if o2 else 'FAIL'}")
    json.dump(res, open(ROOT / "data/processed/landslide/langtang_offsets.json", "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
