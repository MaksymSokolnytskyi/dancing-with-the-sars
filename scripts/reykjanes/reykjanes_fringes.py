"""20 m wrapped NISAR interferograms around Svartsengi - Sundhnukur for one track.

Usage: python scripts/reykjanes/reykjanes_fringes.py TRACK [from YYYYMMDD] [to YYYYMMDD]
"""
import sys
from datetime import datetime

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from pyproj import Transformer

from reykjanes_analysis import OUT, PLACES, pairs_by_track

HALF_KM = 9.0


def main(track, d0="19000101", d1="21000101"):
    d0, d1 = datetime.strptime(d0, "%Y%m%d"), datetime.strptime(d1, "%Y%m%d")
    prs = [p for p in pairs_by_track()[track] if d0 <= p[0] <= d1]
    cols = min(6, len(prs))
    rows = int(np.ceil(len(prs) / cols))
    fig, axes = plt.subplots(rows, cols, figsize=(3.0 * cols, 3.0 * rows + 0.5), squeeze=False)
    for ax in axes.ravel():
        ax.axis("off")
    for ax, (t0, t1, p) in zip(axes.ravel(), prs):
        with h5py.File(p, "r") as f:
            x, y = f["wrp/x"][()], f["wrp/y"][()]
            tr = Transformer.from_crs(4326, int(f.attrs["epsg"]), always_xy=True)
            cx, cy = tr.transform(*PLACES["Svartsengi"])
            ix = np.where(np.abs(x - cx) < HALF_KM * 1000)[0]
            iy = np.where(np.abs(y - cy) < HALF_KM * 1000)[0]
            w = f["wrp/wrappedInterferogram"][iy[0]:iy[-1] + 1, ix[0]:ix[-1] + 1]
            c = f["wrp/coherenceMagnitude"][iy[0]:iy[-1] + 1, ix[0]:ix[-1] + 1]
        ph = np.where(c > 0.2, np.angle(w), np.nan)
        ext = [x[ix[0]], x[ix[-1]], y[iy[-1]], y[iy[0]]]
        ax.imshow(ph, cmap="twilight_shifted", vmin=-np.pi, vmax=np.pi, extent=ext, interpolation="nearest")
        for n in ("Svartsengi", "Grindavik", "Sundhnukur"):
            px, py = tr.transform(*PLACES[n])
            ax.plot(px, py, "k+", ms=6)
        ax.set_title(f"{t0:%d.%m.%y}-{t1:%d.%m.%y}  coh {np.nanmedian(c):.2f}", fontsize=8)
    fig.suptitle(f"Track {track}: 20 m wrapped interferograms, ±{HALF_KM:.0f} km around Svartsengi "
                 "(one colour cycle = 12 cm LOS; pixels with coherence < 0.2 blank)", fontsize=9)
    fig.tight_layout()
    out = OUT / f"fringes_{track}_{d0:%Y%m%d}_{d1:%Y%m%d}.png"
    fig.savefig(out, dpi=110, facecolor="white")
    print(out)


if __name__ == "__main__":
    main(*sys.argv[1:])
