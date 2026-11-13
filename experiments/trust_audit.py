"""Trust audit for any site / pass from the cached GUNW windows (agreed plan point 2).

    uv run python scripts/trust_audit.py --site houston --pass asc [--iono plane]
Writes data/processed/trust/<site>-<pass>[-ionoplane]/{report.md, audit.json, sigma.tif}.
With a GNSS rate table (data/processed/gnss/<site>_rates.csv) it adds the station z-table (P-port P5).
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
import rasterio
from rasterio.warp import transform as warp_transform
from scipy.ndimage import gaussian_filter

from hachure import gridding, nisar, trust

ROOT = Path(__file__).resolve().parents[1]
SITES = {
    "mexico-city": dict(bbox=(-99.32, 19.17, -98.82, 19.58), ref=(-99.14, 19.53, -99.09, 19.57), h=2240),
    "central-valley": dict(bbox=(-120.20, 35.80, -119.00, 36.50), ref=(-119.10, 36.33, -119.02, 36.43), h=80),
    "houston": dict(bbox=(-95.90, 29.55, -95.05, 30.30), ref=(-95.20, 29.60, -95.10, 29.68), h=30),
}


def plane_removed(v_cm):
    """Same robust one-sided plane removal as P-port / the Mexico maps (on a 400 m-smoothed field)."""
    v = v_cm / 100
    ok = np.isfinite(v)
    wsum = gaussian_filter(ok.astype(float), 400 / 80)
    sm = gaussian_filter(np.where(ok, v, 0.0), 400 / 80) / np.maximum(wsum, 1e-6)
    return (v - nisar.fit_ramp(np.where(ok & (wsum > 0.3), sm, np.nan))) * 100


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--site", required=True, choices=list(SITES))
    ap.add_argument("--pass", dest="pas", required=True, choices=["asc", "desc"])
    ap.add_argument("--iono", default="full", choices=["full", "plane"])
    ap.add_argument("--segments", action="store_true", help="P5b: add the chain-gap node term (math M2b) to sigma")
    args = ap.parse_args()
    s = SITES[args.site]
    bbox = s["bbox"]
    cache = ROOT / "data/raw/nisar_cache" / ("bbox_" + "_".join(f"{v:.3f}" for v in bbox) + f"_h{s['h']}")
    found = nisar.search((bbox[0] + bbox[2]) / 2, (bbox[1] + bbox[3]) / 2,
                         direction="ASCENDING" if args.pas == "asc" else "DESCENDING")
    chain = [g for g in nisar.sequential_chain(nisar.dominant_track(found)) if (cache / f"{g.name}.npz").exists()]
    iono = "plane" if args.iono == "plane" else True
    pairs, kept = [], []
    for g in chain:
        try:
            pairs.append(nisar.process(nisar.RawWindow.load(cache / f"{g.name}.npz"), ref_bbox=s["ref"], coh_min=0.3,
                                       iono=iono))
            kept.append(g)
        except ValueError as e:                     # unusable product: skip (the chain gets a gap)
            print(f"SKIP {g.name[:60]}: {e}")
    chain = kept
    dates = [f"{chain[0].ref_date:%Y-%m-%d}"] + [f"{g.sec_date:%Y-%m-%d}" for g in chain]
    a = trust.audit(pairs, dates, segments=args.segments)
    _, _, _, x, y, disps = nisar.stack_rate(pairs)
    okd = np.isfinite(disps)
    nseg = okd[0].astype(int) + (okd[1:] & ~okd[:-1]).sum(0)      # runs of valid pairs per pixel (1 = unbroken)
    nseg = np.where(okd.any(0), nseg, 0)
    epsg = pairs[0].epsg
    name = f"{args.site}-{args.pas}" + ("-ionoplane" if args.iono == "plane" else "") + ("-seg" if args.segments else "")
    out = ROOT / "data/processed/trust" / name
    out.mkdir(parents=True, exist_ok=True)

    rows = None
    gcsv = ROOT / "data/processed/gnss" / f"{args.site}_rates.csv"
    if gcsv.exists():
        g = [r for r in csv.DictReader(open(gcsv)) if r["window"] == args.pas and r["used"] == "True"]
        v = plane_removed(a.rate)
        xs, ys = warp_transform("EPSG:4326", f"EPSG:{epsg}", [float(r["lon"]) for r in g], [float(r["lat"]) for r in g])
        XX, YY = np.meshgrid(x, y)
        tmp = []
        for r, sx, sy in zip(g, xs, ys):
            near = np.hypot(XX - sx, YY - sy) <= 250
            if np.isfinite(v[near]).sum() >= 3:
                tmp.append(dict(station=r["station"], gnss=float(r["rate_cm_yr"]), sigma_gnss=float(r["sigma_cm_yr"]),
                                nis=float(np.nanmedian(v[near])), sigma_nisar=float(np.nanmedian(a.sigma[near]))))
        off = np.median([t["nis"] - t["gnss"] for t in tmp])
        rows = []
        for t in tmp:
            t["nisar"] = t.pop("nis") - off
            t["z"] = (t["nisar"] - t["gnss"]) / np.hypot(t["sigma_nisar"], t["sigma_gnss"])
            rows.append(t)
    rep = trust.report(name, a, rows)
    (out / "report.md").write_text(rep, encoding="utf-8")
    print(rep)
    json.dump(dict(dates=a.dates, node_var=a.node_var.tolist(), pair_var=a.pair_var.tolist(), edge_share=a.edge_share,
                   sigma_median=float(np.nanmedian(a.sigma)), gnss=rows), open(out / "audit.json", "w"), indent=1)
    tr = gridding.transform_from_centres(x, y)
    with rasterio.open(out / "sigma.tif", "w", driver="GTiff", height=a.rate.shape[0], width=a.rate.shape[1], count=3,
                       dtype="float32", crs=f"EPSG:{epsg}", transform=tr, nodata=np.nan, compress="deflate") as ds:
        ds.write((a.rate / 100).astype("float32"), 1)
        ds.write((a.sigma / 100).astype("float32"), 2)
        ds.set_band_description(1, "vertical_rate_m_per_yr_chain_sum")
        ds.set_band_description(2, "rate_1sigma_m_per_yr")
        ds.write(nseg.astype("float32"), 3)
        ds.set_band_description(3, "chain_segments (1 = unbroken, 0 = no data)")
    share = {int(v): float((nseg == v).sum() / max((nseg > 0).sum(), 1)) for v in range(1, 6)}
    print("chain segments per valid pixel: " + ", ".join(f"{v}: {100 * f:.0f}%" for v, f in share.items()) +
          f", >5: {100 * (nseg > 5).sum() / max((nseg > 0).sum(), 1):.0f}%")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7, 6))
    im = ax.imshow(np.ma.masked_equal(nseg, 0), cmap="magma_r", vmin=1, vmax=8, interpolation="nearest")
    plt.colorbar(im, ax=ax, label="times the chain is broken + 1")
    ax.set_title(f"{name}: chain segments per pixel (1 = complete year)")
    fig.tight_layout()
    fig.savefig(out / "segments.png", dpi=110)


if __name__ == "__main__":
    main()
