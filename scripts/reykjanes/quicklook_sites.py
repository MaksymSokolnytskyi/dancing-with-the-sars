"""Quick looks of the downloaded site windows.

Usage: python scripts/reykjanes/quicklook_sites.py {reykjanes|mississippi}
"""
import re
import sys
from pathlib import Path

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "output"


def tag(name):
    m = re.search(r"_(\d{3})_([AD])_(\d{3})_.*?_(\d{8})T\d{6}_\d{8}T\d{6}(?:_(\d{8})T)?", name)
    return f"{m.group(1)}{m.group(2)} f{m.group(3)} {m.group(4)}" + (f"-{m.group(5)}" if m.group(5) else "")


def reykjanes():
    files = sorted((ROOT / "data/reykjanes/gunw").glob("*.h5"), key=lambda p: re.findall(r"\d{8}T", p.name)[0])
    track = "_082_D_"  # one descending track, all its pairs
    sel = [p for p in files if track in p.name][:12]
    n = len(sel)
    fig, axes = plt.subplots(3, (n + 2) // 3, figsize=(4 * ((n + 2) // 3), 9))
    for ax in axes.ravel():
        ax.axis("off")
    for ax, p in zip(axes.ravel(), sel):
        with h5py.File(p, "r") as f:
            w = f["wrp/wrappedInterferogram"][::2, ::2]
            c = f["wrp/coherenceMagnitude"][::2, ::2]
        ph = np.where((c > 0.25) & np.isfinite(c), np.angle(w), np.nan)
        ax.imshow(ph, cmap="twilight_shifted", vmin=-np.pi, vmax=np.pi, interpolation="nearest")
        ax.set_title(tag(p.name), fontsize=8)
    fig.suptitle("Reykjanes: NISAR 20 m wrapped interferograms, track 082 D (one fringe = 12 cm LOS)", fontsize=10)
    fig.tight_layout()
    fig.savefig(OUT / "reykjanes_wrapped_082D.png", dpi=110, facecolor="white")
    print("pairs per track:", {t: sum(t in p.name for p in files) for t in sorted({re.search(r"_(\d{3}_[AD])_", p.name).group(1) for p in files})})


def mississippi():
    files = sorted((ROOT / "data/mississippi/gcov").glob("*.h5"), key=lambda p: re.findall(r"\d{8}T", p.name)[0])
    track = "_105_A_016_"
    sel = [p for p in files if track in p.name]
    pick = [sel[0], sel[len(sel) // 2], sel[-1]]
    fig, axes = plt.subplots(1, 3, figsize=(16, 5.5))
    for ax, p in zip(axes, pick):
        with h5py.File(p, "r") as f:
            hh = f["HHHH_db100"][::2, ::2].astype(float)
            hv = f["HVHV_db100"][::2, ::2].astype(float)
        hh[hh == -32768], hv[hv == -32768] = np.nan, np.nan
        rgb = np.dstack([np.clip((hh / 100 + 25) / 25, 0, 1), np.clip((hv / 100 + 32) / 22, 0, 1),
                         np.clip((hh - hv) / 100 / 15, 0, 1)])
        ax.imshow(np.where(np.isfinite(rgb), rgb, 1))
        ax.set_title(tag(p.name) + "  R=HH G=HV B=HH/HV", fontsize=9)
        ax.axis("off")
    fig.suptitle("Mississippi Delta (Atchafalaya - Wax Lake - Terrebonne): NISAR GCOV", fontsize=10)
    fig.tight_layout()
    fig.savefig(OUT / "mississippi_gcov_rgb.png", dpi=110, facecolor="white")
    print("scenes per track:", {t: sum(t in p.name for p in files) for t in sorted({re.search(r"_(\d{3}_[AD]_\d{3})_", p.name).group(1) for p in files})})


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True)
    {"reykjanes": reykjanes, "mississippi": mississippi}[sys.argv[1]]()
