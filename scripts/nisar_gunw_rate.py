"""NISAR GUNW pairs -> vertical subsidence rate GeoTIFF (input for hachure_subsidence.py --rate).

    # 0) once: create a token at https://urs.earthdata.nasa.gov  (Generate Token), then
    #    set EARTHDATA_TOKEN in your own shell. Never paste it into chats/docs.
    uv run python scripts/nisar_gunw_rate.py --list                 # what ASF has (no auth)
    uv run python scripts/nisar_gunw_rate.py --inspect              # HDF5 tree of 1st pair
    uv run python scripts/nisar_gunw_rate.py --max-pairs 6          # rate from 6 pairs
    uv run python scripts/nisar_gunw_rate.py --local f1.h5 f2.h5    # already-downloaded files
"""
from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np
import rasterio
from rasterio.transform import from_origin

from hachure import nisar

ROOT = Path(__file__).resolve().parents[1]

# Mexico City defaults. Reference = Sierra de Guadalupe volcanic rock (assumed stable).
BBOX = (-99.32, 19.17, -98.82, 19.58)
REF_BBOX = (-99.14, 19.53, -99.09, 19.57)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", default="mexico-city")
    ap.add_argument("--bbox", nargs=4, type=float, default=BBOX)
    ap.add_argument("--ref-bbox", nargs=4, type=float, default=REF_BBOX)
    ap.add_argument("--direction", choices=["ASCENDING", "DESCENDING"], default=None)
    ap.add_argument("--start")
    ap.add_argument("--end")
    ap.add_argument("--max-pairs", type=int, default=6)
    ap.add_argument("--coh-min", type=float, default=0.3)
    ap.add_argument("--no-iono", action="store_true")
    ap.add_argument("--iono-plane", action="store_true", help="subtract only the planar part of each iono screen")
    ap.add_argument("--iono-rplane", action="store_true", help="H11: subtract a robust (island-proof) plane of the iono screen")
    ap.add_argument("--iono-robust", action="store_true", help="H10: subtract a 10 km median-smoothed iono screen")
    ap.add_argument("--flip-sign", action="store_true")
    ap.add_argument("--height", type=float, default=2240.0,
                    help="terrain height (m) for the incidence-angle layer; Mexico basin ~2240")
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--inspect", action="store_true")
    ap.add_argument("--local", nargs="+", help="local GUNW .h5 files instead of ASF")
    args = ap.parse_args()

    if args.local:
        sources = [(p, nisar._parse_dates(Path(p).name)) for p in args.local]
        sources = [(p, (sec - ref).total_seconds() / 86400) for p, (ref, sec) in sources]
    else:
        lon = (args.bbox[0] + args.bbox[2]) / 2
        lat = (args.bbox[1] + args.bbox[3]) / 2
        found = nisar.search(lon, lat, direction=args.direction, start=args.start, end=args.end)
        track = nisar.dominant_track(found)
        print(f"ASF: {len(found)} GUNW pairs; dominant track {len(track)} pairs")
        if args.list:
            for g in found:
                gb = f"{g.size_bytes / 1e9:.1f} GB" if g.size_bytes else "?"
                print(f"  {g.direction[:4]} p{g.path} f{g.frame}  {g.ref_date:%Y-%m-%d} -> "
                      f"{g.sec_date:%Y-%m-%d} ({g.days:.0f} d)  {gb}")
            return
        chain = nisar.sequential_chain(track)
        print(f"sequential chain: {len(chain)} pairs, {chain[0].ref_date:%Y-%m-%d} -> "
              f"{chain[-1].sec_date:%Y-%m-%d}" if chain else "no chain")
        chain = chain[-args.max_pairs:]
        sources = [(g.url, g.days) for g in chain]
        print(f"using {len(chain)} pairs: {chain[0].ref_date:%Y-%m-%d} -> {chain[-1].sec_date:%Y-%m-%d}")

    if args.inspect:
        with nisar.open_h5(sources[0][0]) as h5:
            print(nisar.tree(h5))
        return

    bbox_key = "bbox_" + "_".join(f"{v:.3f}" for v in args.bbox) + f"_h{int(args.height)}"
    cache = ROOT / "data" / "raw" / "nisar_cache" / bbox_key
    cache.mkdir(parents=True, exist_ok=True)
    pairs, dates = [], []
    for src, days in sources:
        t = time.time()
        name = Path(src).name.replace(".h5", "")
        cpath = cache / f"{name}.npz"
        if cpath.exists():
            raw = nisar.RawWindow.load(cpath)
            how = "cache"
        else:
            with nisar.open_h5(src) as h5:
                raw = nisar.read_raw(h5, tuple(args.bbox), days, height_m=args.height)
            raw.save(cpath)
            how = "ASF" if src.startswith("http") else "local"
        try:
            p = nisar.process(raw, ref_bbox=tuple(args.ref_bbox), coh_min=args.coh_min,
                              iono=("rplane" if args.iono_rplane else "robust" if args.iono_robust else "plane" if args.iono_plane else not args.no_iono), flip_sign=args.flip_sign)
        except ValueError as e:            # unusable product (e.g. 0% valid): skip the pair, chain gets a gap
            print(f"  SKIP {name[:60]}...: {e}", flush=True)
            continue
        ok = np.isfinite(p.disp_v)
        print(f"  {name[:60]}...  {ok.mean() * 100:.0f}% valid, median "
              f"{np.nanmedian(p.disp_v) * 100:+.2f} cm / {days:.0f} d  ({how}, {time.time() - t:.0f}s)",
              flush=True)
        pairs.append(p)
        ref_d, sec_d = nisar._parse_dates(name)
        dates.append((ref_d, sec_d))

    rate, coh, frac, x, y, disps = nisar.stack_rate(pairs)
    out = ROOT / "data" / "processed" / args.name
    out.mkdir(parents=True, exist_ok=True)
    dx, dy = x[1] - x[0], y[0] - y[1]
    transform = from_origin(x[0] - dx / 2, y[0] + dy / 2, dx, dy)
    tif = out / "nisar_rate.tif"
    with rasterio.open(tif, "w", driver="GTiff", height=rate.shape[0], width=rate.shape[1],
                       count=2, dtype="float32", crs=f"EPSG:{pairs[0].epsg}",
                       transform=transform, nodata=np.nan, compress="deflate") as ds:
        ds.write(rate.astype("float32"), 1)
        ds.write(coh.astype("float32"), 2)
        ds.set_band_description(1, "vertical_rate_m_per_yr")
        ds.set_band_description(2, "mean_coherence")
    print(tif)

    # time series for animation: cumulative vertical displacement at every chain date
    chain_ok = all(dates[i][1].date() == dates[i + 1][0].date() for i in range(len(dates) - 1))
    if not chain_ok:
        print("  (pairs are not a sequential chain - series is approximate)")
    cum = nisar.cumulative(disps, [p.days for p in pairs], rate)
    series_dates = [dates[0][0]] + [d[1] for d in dates]
    np.savez_compressed(out / "nisar_series.npz", cum=cum.astype("float32"),
                        dates=np.array([d.isoformat() for d in series_dates]),
                        x=x, y=y, epsg=pairs[0].epsg, coherence=coh.astype("float32"))
    print(out / "nisar_series.npz", f"({len(series_dates)} dates)")


if __name__ == "__main__":
    main()
