"""Write small FAKE GUNW-like HDF5 files to test the loader end-to-end without Earthdata.

Layout mirrors what nisar.py expects (it discovers paths by dataset name, so this is
a plausible guess, not the official spec). Phase = synthetic Mexico City rate * dt,
plus a per-pair constant offset (unwrapping ambiguity), ionosphere ramp (stored in
ionospherePhaseScreen so correction should remove it), noise, and a bad component.

    uv run python tests/make_fake_gunw.py   -> data/raw/fake_gunw/*.h5
"""
from __future__ import annotations

from datetime import datetime, timedelta
from pathlib import Path

import h5py
import numpy as np
from rasterio.transform import from_origin

from hachure import dem, nisar, synthetic

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "raw" / "fake_gunw"


def main():
    dem_path = dem.fetch_copdem((-99.40, 19.10, -98.75, 19.65),
                                ROOT / "data" / "raw" / "copdem_fake-gunw_80m.tif", res=80.0)
    z, t, crs = dem.load(dem_path)
    rate, coh = synthetic.subsidence_rate(t, crs, z.shape, z)
    epsg = crs.to_epsg()
    h, w = z.shape
    xs = t.c + (np.arange(w) + 0.5) * t.a
    ys = t.f + (np.arange(h) + 0.5) * t.e
    lam = nisar.C / 1.2575e9
    inc = 37.0
    rng = np.random.default_rng(1)
    OUT.mkdir(parents=True, exist_ok=True)
    start = datetime(2026, 9, 2, 1, 13, 29)
    for k in range(4):
        ref, sec = start + timedelta(days=12 * k), start + timedelta(days=12 * (k + 1))
        dv = rate * 12 / 365.25
        d_los = dv * np.cos(np.radians(inc))
        phase = d_los / (nisar.LOS_SIGN * lam / (4 * np.pi))
        iono = np.broadcast_to(3.0 * (np.linspace(-1, 1, w)[None, :] + 0.5 * k), z.shape)
        phase = phase + iono + 7.3 * (k - 1.5) + 0.15 * rng.standard_normal(z.shape)
        comp = np.ones(z.shape, np.uint16)
        comp[: h // 8, : w // 8] = 2                                       # isolated island
        phase[comp == 2] += 2 * np.pi * 3                                  # wrong ambiguity
        f = lambda d: d.strftime("%Y%m%dT%H%M%S")
        name = (f"NISAR_L2_PR_GUNW_FAKE_113_D_079_0{k}_4000_SH_{f(ref)}_{f(ref)}_"
                f"{f(sec)}_{f(sec)}_P00000_N_F_J_001.h5")
        with h5py.File(OUT / name, "w") as hf:
            g = hf.create_group("science/LSAR/GUNW/grids/frequencyA/unwrappedInterferogram")
            g["xCoordinates"], g["yCoordinates"] = xs, ys
            pr = g.create_dataset("projection", data=np.uint32(epsg))
            pr.attrs["epsg_code"] = epsg
            hf["science/LSAR/GUNW/grids/frequencyA/centerFrequency"] = 1.2575e9
            p = g.create_group("HH")
            p.create_dataset("unwrappedPhase", data=phase.astype("float32"), chunks=(128, 128))
            p.create_dataset("coherenceMagnitude", data=coh.astype("float32"), chunks=(128, 128))
            p.create_dataset("connectedComponents", data=comp, chunks=(128, 128))
            p.create_dataset("ionospherePhaseScreen", data=iono.astype("float32"), chunks=(128, 128))
            rg = hf.create_group("science/LSAR/GUNW/metadata/radarGrid")
            rg["xCoordinates"] = xs[::20]
            rg["yCoordinates"] = ys[::20]
            rg["incidenceAngle"] = np.full((3, len(ys[::20]), len(xs[::20])), inc)
        print(OUT / name)
    from rasterio import open as ropen
    with ropen(OUT / "truth_rate.tif", "w", driver="GTiff", height=h, width=w, count=1,
               dtype="float32", crs=crs, transform=from_origin(t.c, t.f, t.a, -t.e)) as ds:
        ds.write(rate.astype("float32"), 1)


if __name__ == "__main__":
    main()
