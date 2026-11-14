"""Kilauea summit, five lava-fountain episodes of 2026 seen by NISAR (tracks 072D, 151A).

For every GUNW pair whose interval contains an episode: LOS displacement (ionosphere removed,
coherence >= 0.3, the connected component of the summit, referenced to the median of pixels
> 7 km from Halema'uma'u), the coherence, and the HH backscatter change (GCOV, the last
acquisition before vs the first after the episode, same track).

Writes output/kilauea/episodes.png and episodes.json.
Usage: python scripts/reykjanes/kilauea_episodes.py
"""
import json
import re
from datetime import datetime
from pathlib import Path

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from pyproj import Transformer

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data" / "kilauea"
OUT = ROOT / "output" / "kilauea"
C = 299792458.0
VENT = (-155.287, 19.405)          # Halema'uma'u
EVENTS = {"20260124": "еп. 41 (24.01)", "20260215": "еп. 42 (15.02)", "20260310": "еп. 43 (10.03)",
          "20260627": "еп. 50 (27.06)", "20260716": "еп. 51 (16.07)"}
LOOK = 10
STAMP = re.compile(r"_(\d{8})T\d{6}")
TRACK = re.compile(r"_(\d{3})_([AD])_\d{3}_")


def track(n):
    m = TRACK.search(n)
    return m.group(1) + m.group(2)


def grid(x, y, epsg):
    X, Y = np.meshgrid(x, y)
    vx, vy = Transformer.from_crs(4326, epsg, always_xy=True).transform(*VENT)
    return X, Y, np.hypot(X - vx, Y - vy) / 1000, (vx, vy)


def los_pair(path):
    with h5py.File(path, "r") as f:
        ph = f["unw/unwrappedPhase"][()].astype(float) - f["unw/ionospherePhaseScreen"][()]
        coh = f["unw/coherenceMagnitude"][()]
        cc = f["unw/connectedComponents"][()]
        x, y, epsg = f["unw/x"][()], f["unw/y"][()], int(f.attrs["epsg"])
        wl = C / float(f.attrs["centerFrequency"])
    X, Y, dist, v = grid(x, y, epsg)
    ok = np.isfinite(ph) & (cc > 0) & (coh >= 0.3)
    near = ok & (dist < 3)
    comp = np.bincount(cc[near]).argmax() if near.any() else np.bincount(cc[ok]).argmax()
    ok &= cc == comp
    d = np.where(ok, -ph * wl / (4 * np.pi) * 100, np.nan)        # cm, + = towards the satellite
    far = ok & (dist > 7)
    d -= np.nanmedian(d[far])
    return dict(d=d, coh=coh, dist=dist, ext=(x[0], x[-1], y[-1], y[0]), vent=v, epsg=epsg)


def gcov(path):
    with h5py.File(path, "r") as f:
        a = f["HHHH_db100"][()].astype(np.float32)
        x, y, epsg = f["x"][()], f["y"][()], int(f.attrs["epsg"])
    a[a == -32768] = np.nan
    # 10 x 10 multilook in linear power (10 m -> 100 m): single-look speckle alone swings +/-3 dB
    k = LOOK
    ny, nx = a.shape[0] // k * k, a.shape[1] // k * k
    p = 10 ** (a[:ny, :nx] / 1000)
    p = np.nanmean(p.reshape(ny // k, k, nx // k, k), axis=(1, 3))
    return 10 * np.log10(p), x[:nx].reshape(-1, k).mean(1), y[:ny].reshape(-1, k).mean(1), epsg


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    gunw = sorted(DATA.glob("gunw/*.h5"))
    gc = sorted(DATA.glob("gcov/*.h5"))
    rows, panels = [], {}
    for ev, lab in EVENTS.items():
        for p in gunw:
            d = STAMP.findall(p.name)
            if not d[0] < ev <= d[2]:
                continue
            t = track(p.name)
            r = los_pair(p)
            dist = r["dist"]
            inner = np.isfinite(r["d"]) & (dist < 1.5)
            ring = np.isfinite(r["d"]) & (dist > 2.5) & (dist < 4.5)
            row = dict(event=ev, label=lab, track=t, ref=d[0], sec=d[2],
                       days=(datetime.strptime(d[2], "%Y%m%d") - datetime.strptime(d[0], "%Y%m%d")).days,
                       summit_los_cm=round(float(np.nanmedian(r["d"][ring])), 2) if ring.sum() > 50 else None,
                       caldera_valid=round(float(inner.mean() / max((dist < 1.5).mean(), 1e-9)), 2),
                       coh_caldera=round(float(np.nanmedian(r["coh"][dist < 1.5])), 2),
                       coh_far=round(float(np.nanmedian(r["coh"][dist > 7])), 2),
                       los_min_cm=round(float(np.nanpercentile(r["d"][dist < 6], 1)), 1),
                       los_max_cm=round(float(np.nanpercentile(r["d"][dist < 6], 99)), 1))
            # backscatter change, same track: last acquisition before vs first on/after the episode
            acq = sorted((STAMP.findall(g.name)[0], g) for g in gc if track(g.name) == t)
            b = [g for s, g in acq if s < ev][-1:]
            a = [g for s, g in acq if s >= ev][:1]
            if b and a:
                A0, x, y, epsg = gcov(b[0])
                A1, *_ = gcov(a[0])
                dd = A1 - A0
                _, _, gd, _ = grid(x, y, epsg)
                dd -= np.nanmedian(dd[gd > 7])          # whole-scene shifts (rain, soil moisture, calibration)
                big = np.abs(dd) > 2
                row.update(gcov_before=STAMP.findall(b[0].name)[0], gcov_after=STAMP.findall(a[0].name)[0],
                           gcov_changed_caldera=round(float(np.nanmean(big[gd < 2])), 3),
                           gcov_changed_far=round(float(np.nanmean(big[gd > 7])), 3))
                panels[(ev, t, "gcov")] = (dd, (x[0], x[-1], y[-1], y[0]))
            panels[(ev, t, "los")] = (r["d"], r["ext"], r["vent"])
            panels[(ev, t, "coh")] = (r["coh"], r["ext"], r["vent"])
            rows.append(row)
            print(f"{lab} {t} {d[0]}->{d[2]} ({row['days']} d): summit ring LOS {row['summit_los_cm']} cm, "
                  f"LOS 1-99 % within 6 km {row['los_min_cm']}..{row['los_max_cm']} cm, coherence caldera "
                  f"{row['coh_caldera']} / far {row['coh_far']}, unwrapped in caldera {row['caldera_valid']:.0%}"
                  + (f"; HH change >2 dB: caldera {row['gcov_changed_caldera']:.0%}, far {row['gcov_changed_far']:.1%}"
                     if "gcov_changed_caldera" in row else ""))
    (OUT / "episodes.json").write_text(json.dumps(rows, indent=1, ensure_ascii=False), encoding="utf-8")
    plot(rows, panels)


def plot(rows, panels):
    keys = sorted({(r["event"], r["track"]) for r in rows})
    fig, ax = plt.subplots(3, len(keys), figsize=(2.6 * len(keys), 8.4))
    for j, (ev, t) in enumerate(keys):
        r = next(x for x in rows if x["event"] == ev and x["track"] == t)
        d, ext, v = panels[(ev, t, "los")]
        im0 = ax[0, j].imshow(d, extent=ext, cmap="RdBu_r", vmin=-8, vmax=8)
        ax[0, j].set_title(f"{EVENTS[ev]}\n{t} {r['ref'][4:]}→{r['sec'][4:]}", fontsize=8)
        c, ext, v = panels[(ev, t, "coh")]
        im1 = ax[1, j].imshow(c, extent=ext, cmap="gray", vmin=0, vmax=1)
        if (ev, t, "gcov") in panels:
            g, ext2 = panels[(ev, t, "gcov")]
            im2 = ax[2, j].imshow(g, extent=ext2, cmap="PuOr", vmin=-4, vmax=4)
        for i in range(3):
            ax[i, j].plot(*v, "k+", ms=8)
            ax[i, j].set_xticks([]); ax[i, j].set_yticks([])
    ax[0, 0].set_ylabel("LOS, см (+ до супутника)", fontsize=8)
    ax[1, 0].set_ylabel("когерентність", fontsize=8)
    ax[2, 0].set_ylabel("зміна HH, дБ (після − до, 100 м)", fontsize=8)
    fig.colorbar(im0, ax=ax[0], shrink=.8, label="см"); fig.colorbar(im1, ax=ax[1], shrink=.8)
    fig.colorbar(im2, ax=ax[2], shrink=.8, label="дБ")
    fig.suptitle("Кілауеа: епізоди фонтанування 2026 у даних NISAR (+ = кратер Халемаумау)", fontsize=11)
    fig.savefig(OUT / "episodes.png", dpi=110, facecolor="white", bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
