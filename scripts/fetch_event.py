"""Grab NISAR GUNW windows around an earthquake (co-seismic pairs) + a quick-look PNG.

Finds the GUNW pairs whose dates span the event, reads only a box around the epicentre
(no full 2.3 GB download) and caches it as RawWindow .npz, same format as the Mexico cache.
Also caches up to --pre pre-event pairs of the same track (no quake inside = pure noise
samples for the synthetic-data generator).

    uv run python scripts/fetch_event.py --name nagqu-2026 --lon 93.7 --lat 31.9 --date 2026-10-22
"""
from __future__ import annotations

import argparse
import time
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np

from hachure import nisar

ROOT = Path(__file__).resolve().parents[1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", required=True)
    ap.add_argument("--lon", type=float, required=True)
    ap.add_argument("--lat", type=float, required=True)
    ap.add_argument("--date", required=True, help="event date YYYY-MM-DD (UTC)")
    ap.add_argument("--half", type=float, default=0.35, help="half box size, degrees")
    ap.add_argument("--max-days", type=float, default=36)
    ap.add_argument("--pre", type=int, default=2, help="pre-event pairs per track to cache")
    args = ap.parse_args()

    t = datetime.fromisoformat(args.date)
    bbox = (args.lon - args.half, args.lat - args.half, args.lon + args.half, args.lat + args.half)
    found = nisar.search(args.lon, args.lat, start=(t - timedelta(days=120)).strftime("%Y-%m-%d"),
                         end=(t + timedelta(days=60)).strftime("%Y-%m-%d"))
    co = [g for g in found if g.ref_date < t + timedelta(days=1) and g.sec_date > t
          and g.days <= args.max_days]
    pick = list(co)
    for track in {(g.direction, g.path, g.frame) for g in co}:
        pre = [g for g in found if (g.direction, g.path, g.frame) == track and g.sec_date < t
               and g.days <= 24]
        pick += sorted(pre, key=lambda g: g.sec_date)[-args.pre:]
    print(f"{args.name}: {len(co)} co-seismic pairs, {len(pick) - len(co)} pre-event pairs")

    out = ROOT / "data" / "raw" / "nisar_events" / args.name
    out.mkdir(parents=True, exist_ok=True)
    for g in pick:
        kind = "co" if g in co else "pre"
        p = out / f"{kind}_{g.name}.npz"
        t0 = time.time()
        if not p.exists():
            with nisar.open_h5(g.url) as h5:
                raw = nisar.read_raw(h5, bbox, g.days)
            raw.save(p)
        raw = nisar.RawWindow.load(p)
        ph = raw.phase - (raw.iono if raw.iono is not None else 0)
        ok = np.isfinite(ph) & (raw.comp > 0)
        print(f"  {kind:3s} {g.direction[:4]} p{g.path} f{g.frame} {g.ref_date:%Y-%m-%d}->"
              f"{g.sec_date:%Y-%m-%d} ({g.days:.0f} d) {ph.shape} valid {ok.mean() * 100:.0f}% "
              f"phase span {np.nanpercentile(ph[ok], 99) - np.nanpercentile(ph[ok], 1):.0f} rad "
              f"({time.time() - t0:.0f}s)", flush=True)
        quicklook(raw, out / f"{kind}_{g.direction[:4]}_p{g.path}_{g.ref_date:%m%d}-{g.sec_date:%m%d}.png",
                  f"{args.name}  {kind}  {g.ref_date:%Y-%m-%d} -> {g.sec_date:%Y-%m-%d}")


def quicklook(raw, path, title):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ph = raw.phase - (raw.iono if raw.iono is not None else 0)
    ph = np.where(raw.comp > 0, ph, np.nan)
    fig, ax = plt.subplots(1, 2, figsize=(12, 5.6))
    ax[0].imshow(np.angle(np.exp(1j * ph)), cmap="twilight", interpolation="nearest")
    ax[0].set_title("rewrapped phase (1 cycle = 12 cm LOS)")
    lo, hi = np.nanpercentile(ph, [2, 98])
    im = ax[1].imshow(ph * raw.wavelength / (4 * np.pi) * 100, cmap="RdBu",
                      vmin=lo * raw.wavelength / (4 * np.pi) * 100,
                      vmax=hi * raw.wavelength / (4 * np.pi) * 100, interpolation="nearest")
    ax[1].set_title("unwrapped, cm along line of sight")
    fig.colorbar(im, ax=ax[1], shrink=0.8)
    for a in ax:
        a.axis("off")
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(path, dpi=80)
    plt.close(fig)


if __name__ == "__main__":
    main()
