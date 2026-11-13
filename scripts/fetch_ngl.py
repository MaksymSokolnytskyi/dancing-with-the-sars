"""Download NGL (Nevada Geodetic Lab) GNSS daily series (tenv3, IGS20) for all stations in a bbox.

    uv run python scripts/fetch_ngl.py --name central-valley --bbox -120.20 35.80 -119.00 36.50
"""
from __future__ import annotations

import argparse
import time
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
UA = {"User-Agent": "space-apps-research"}
HOLD = "https://geodesy.unr.edu/NGLStationPages/DataHoldings.txt"
SERIES = "https://geodesy.unr.edu/gps_timeseries/IGS20/tenv3/IGS20/{}.tenv3"


def get(url):
    with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=120) as r:
        return r.read()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", required=True)
    ap.add_argument("--bbox", nargs=4, type=float, required=True)
    ap.add_argument("--min-end", default="2026-06-01", help="station must have data after this date")
    args = ap.parse_args()
    out = ROOT / "data" / "raw" / "gnss" / args.name
    out.mkdir(parents=True, exist_ok=True)
    hold = ROOT / "data" / "raw" / "gnss" / "DataHoldings.txt"
    if not hold.exists():
        hold.write_bytes(get(HOLD))
    lon0, lat0, lon1, lat1 = args.bbox
    sel = []
    for line in hold.read_text(errors="replace").splitlines()[1:]:
        f = line.split()
        if len(f) < 9:
            continue
        try:
            lat, lon = float(f[1]), float(f[2])
        except ValueError:
            continue
        lon = (lon + 180) % 360 - 180
        if lon0 <= lon <= lon1 and lat0 <= lat <= lat1 and f[8] >= args.min_end:
            sel.append((f[0], lon, lat, f[7], f[8]))
    print(f"{args.name}: {len(sel)} stations with data after {args.min_end}")
    for sta, lon, lat, d0, d1 in sel:
        p = out / f"{sta}.tenv3"
        if p.exists():
            continue
        try:
            p.write_bytes(get(SERIES.format(sta)))
            print(f"  {sta} {lon:.3f} {lat:.3f} {d0}..{d1}", flush=True)
        except Exception as e:  # noqa: BLE001
            print(f"  {sta} failed: {e}")
        time.sleep(0.5)
    (out / "stations.txt").write_text("\n".join(f"{s} {lo:.5f} {la:.5f} {a} {b}" for s, lo, la, a, b in sel))


if __name__ == "__main__":
    main()
