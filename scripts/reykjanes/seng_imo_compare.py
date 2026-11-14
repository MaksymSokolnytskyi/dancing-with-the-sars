"""SENG (Svartsengi) GNSS in the layout of the Icelandic Met Office plot, with NISAR on top.

The Icelandic Met Office (Vedurstofa Islands) publishes daily north / east / up plots of station SENG
in the North-American-plate frame; the user supplied the plot of Nov 2025 - 9 Nov 2026 (up ~230 mm,
north ~62 mm, east ~40 mm). Here:
  - GNSS: the same receiver, Nevada Geodetic Laboratory daily series (IGS20), shown (a) as is, zeroed at
    the first NISAR date, and (b) relative to NYLA (the far-field station that NISAR is referenced to)
  - NISAR: vertical and east components from the two-track decomposition (159A + 082D, 240 m cells,
    median of 3 x 3 cells at SENG; output/reykjanes/breathing_data.npz), referenced to the NYLA area
  - IMO: end values read off the supplied plot (markers), for scale only.
NISAR cannot give north (both tracks look east-west). NISAR is referenced to the NYLA area, so the
"SENG minus NYLA" curves are the fair comparison for NISAR; frames differ for east / north
(NGL: IGS20, IMO: North-American plate), not for up.

Writes output/reykjanes/seng_imo_compare.png and .json.
"""
import json
from datetime import datetime

import matplotlib
matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np

from reykjanes_gnss import OUT, read_tenv3, relative

IMO_END = dict(date=datetime(2026, 11, 9), north=62, east=40, up=230)     # mm, read off the IMO plot (user-supplied)


def main():
    d = np.load(OUT / "breathing_data.npz", allow_pickle=True)
    dates = [datetime.fromisoformat(str(x)[:10]) for x in d["dates"]]
    t, x, ll = read_tenv3("SENG")
    tr, xr, _ = relative("SENG", "NYLA")
    i = np.unravel_index(np.argmin((d["lon"] - ll[0]) ** 2 + ((d["lat"] - ll[1]) * 2.2) ** 2), d["lon"].shape)
    sl = (slice(i[0] - 1, i[0] + 2), slice(i[1] - 1, i[1] + 2))
    nis = {c: np.array([np.nanmedian(d[c][k][sl]) for k in range(len(dates))]) * 10 for c in ("up", "east")}   # mm
    nis = {c: v - v[0] for c, v in nis.items()}
    t0 = dates[0]

    def zero(tt, xx):
        k = np.argmin(np.abs(np.array([(a - t0).days for a in tt])))
        m = (tt >= t0 - (tt[k] - t0) if False else tt >= tt[k])
        return tt[m], (xx[m] - xx[k]) * 10                                                    # mm

    ta, xa = zero(t, x)
    tb, xb = zero(tr, xr)
    fig, axs = plt.subplots(3, 1, figsize=(8.5, 11), sharex=True)
    stats = {}
    for a, (k, nm) in zip(axs, ((1, "north"), (0, "east"), (2, "up"))):
        a.plot(ta, xa[:, k], ".", ms=3, color="tab:red", label="GPS SENG (NGL, IGS20)")
        a.plot(tb, xb[:, k], ".", ms=2.5, color="0.55", label="GPS SENG мінус NYLA")
        if nm in nis:
            a.plot(dates, nis[nm], "o-", color="tab:blue", ms=4, lw=1.2, label="NISAR (159A + 082D), відносно району NYLA")
            g = np.array([np.nanmedian(xb[np.abs(np.array([(q - dd).days for q in tb])) <= 3, k]) for dd in dates])
            ok = np.isfinite(g)
            stats[nm] = dict(nisar_end_mm=round(float(nis[nm][-1]), 1), gps_rel_nyla_end_mm=round(float(g[ok][-1]), 1),
                             rms_diff_mm=round(float(np.sqrt(np.mean((nis[nm][ok] - g[ok]) ** 2))), 1),
                             nisar_last_date=dates[-1].date().isoformat())
        a.plot(IMO_END["date"], IMO_END[nm], "D", color="limegreen", mec="k", ms=8, label="IMO, кінцеве значення з графіка")
        a.axhline(0, color="k", lw=.6)
        a.set_ylabel({"north": "північ, мм", "east": "схід, мм", "up": "вгору, мм"}[nm])
        a.grid(alpha=.3)
        stats.setdefault(nm, {})["gps_end_mm"] = round(float(xa[-1, k]), 1)
    axs[0].legend(fontsize=8, loc="upper left")
    axs[2].legend(fontsize=8, loc="upper left")
    axs[2].xaxis.set_major_formatter(mdates.DateFormatter("%b\n%Y"))
    s = stats["up"]
    fig.suptitle(f"Свартсенгі (SENG): GPS і NISAR, як на графіку Ісландського метеобюро\n"
                 f"вгору: GPS {s['gps_end_mm']:.0f} мм (до {ta[-1].date()}),\nNISAR {s['nisar_end_mm']:.0f} мм "
                 f"(до {s['nisar_last_date']}), IMO ~{IMO_END['up']} мм (до 09.11.2026); розбіжність NISAR − GPS {s['rms_diff_mm']:.0f} мм (rms)",
                 fontsize=10)
    fig.text(0.5, 0.005, "Схід/північ GPS — у системі ITRF (NGL), IMO — відносно Північноамериканської плити; NISAR не бачить руху на північ.\n"
             "Сірі точки (SENG мінус NYLA) — те, що має бачити NISAR, бо він прив’язаний до району NYLA.", ha="center", fontsize=8)
    fig.tight_layout(rect=(0, 0.03, 1, 0.95))
    fig.savefig(OUT / "seng_imo_compare.png", dpi=105, facecolor="white")
    (OUT / "seng_imo_compare.json").write_text(json.dumps(dict(stats=stats, imo_end=dict(IMO_END, date="2026-11-09")), indent=1), encoding="utf-8")
    print(json.dumps(stats, indent=1))


if __name__ == "__main__":
    main()
