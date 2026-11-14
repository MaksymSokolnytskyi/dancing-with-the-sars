"""Which NISAR polarization shows the 2026 Kilauea fountain-episode changes best?

NISAR GCOV over Kilauea is dual-pol DH: HH (co-pol: surface + double-bounce scattering) and HV
(cross-pol: volume / multiple scattering). The products carry only the diagonal terms HHHH and
HVHV (no HH-HV correlation), so a full H/alpha decomposition is not possible; what we can compare:

  HH         co-pol backscatter change
  HV         cross-pol backscatter change
  HV/HH      cross-pol ratio (dB difference): independent of the absolute calibration of a date and
             of soil-moisture changes that scale both channels, sensitive to the scattering MECHANISM
             (smooth tephra blanket -> less volume scattering; rough a'a lava -> more)
  span       HH + HV total power

For every before/after pair (2 tracks x 5 episodes, 100 m multilook) and every channel:
  - contrast-to-noise CNR = median change in a zone / robust std of the change > 7 km from the vent
    (zones: caldera floor < 1 km = new lava; SW tephra sector 195-255 deg, 1.5-6 km; NE sector as a
    control - upwind, no tephra expected)
  - detection: share of zone cells beyond the threshold that flags 5 % of background cells.

Writes output/kilauea/polarization.png and polarization.json.
Usage: python scripts/reykjanes/kilauea_polarization.py
"""
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from pyproj import Transformer

from kilauea_episodes import EVENTS, STAMP, VENT, track
from kilauea_tephra_figure import DATA, HALF, OUT, STEP, load, on_grid

CH = ("HH", "HV", "HV/HH", "HH+HV")
CH_UA = {"HH": "HH (ко-поляризація)", "HV": "HV (крос-поляризація)", "HV/HH": "відношення HV/HH", "HH+HV": "сумарна потужність HH+HV"}


def main():
    vx, vy = Transformer.from_crs(4326, 32604, always_xy=True).transform(*VENT)
    gx = vx + np.arange(-HALF, HALF + 1, STEP)
    gy = vy + np.arange(HALF, -HALF - 1, -STEP)
    E, N = np.meshgrid(gx - vx, gy - vy)
    R = np.hypot(E, N) / 1000
    AZ = np.degrees(np.arctan2(E, N)) % 360
    ring = (R > 1.5) & (R < 6)
    zones = {"caldera": R < 1.0, "SW": ring & (AZ >= 195) & (AZ < 255), "NE": ring & (AZ >= 15) & (AZ < 75)}
    far = R > 7
    files = sorted(DATA.glob("*.h5"))
    ch = {c: [] for c in CH}
    labels, spans = [], []
    for ev, lab in EVENTS.items():
        for t in sorted({track(f.name) for f in files}):
            acq = sorted((STAMP.findall(f.name)[0], f) for f in files if track(f.name) == t)
            b = [(s, f) for s, f in acq if s < ev][-1:]
            a = [(s, f) for s, f in acq if s >= ev][:1]
            if not (b and a):
                continue
            img = {}
            for when, (s, f) in (("b", b[0]), ("a", a[0])):
                hh = on_grid(*load(f, "HHHH"), gx, gy)
                hv = on_grid(*load(f, "HVHV"), gx, gy)
                img[when] = {"HH": hh, "HV": hv, "HV/HH": hv - hh,
                             "HH+HV": 10 * np.log10(10 ** (hh / 10) + 10 ** (hv / 10))}
            for c in CH:
                d = img["a"][c] - img["b"][c]
                ch[c].append(d - np.nanmedian(d[far]))
            days = (np.datetime64(f"{a[0][0][:4]}-{a[0][0][4:6]}-{a[0][0][6:]}") - np.datetime64(f"{b[0][0][:4]}-{b[0][0][4:6]}-{b[0][0][6:]}")).astype(int)
            labels.append(f"{lab} {t} ({days} дн.)")
            spans.append(int(days))
    D = {c: np.array(v) for c, v in ch.items()}
    n = len(labels)

    def robust_sd(x):
        x = x[np.isfinite(x)]
        return 1.4826 * np.median(np.abs(x - np.median(x)))

    stats = {}
    for c in CH:
        cnr = {z: [] for z in zones}
        det = {z: [] for z in zones}
        for d in D[c]:
            sd = robust_sd(d[far])
            lo, hi = np.nanquantile(d[far], [0.025, 0.975])          # 5 % of background cells flagged
            for z, m in zones.items():
                v = d[m]
                v = v[np.isfinite(v)]
                cnr[z].append(float(np.median(v) / sd))
                det[z].append(float(np.mean((v < lo) | (v > hi))))
        stats[c] = dict(cnr=cnr, detect=det, noise_sd_dB=[round(robust_sd(d[far]), 2) for d in D[c]],
                        mean_change_dB={z: [round(float(np.nanmedian(d[m])), 2) for d in D[c]] for z, m in zones.items()})
    summary = {c: {z: dict(cnr_median=round(float(np.median(stats[c]["cnr"][z])), 2),
                           cnr_abs_median=round(float(np.median(np.abs(stats[c]["cnr"][z]))), 2),
                           detect_median=round(float(np.median(stats[c]["detect"][z])) * 100, 1),
                           sign_consistent=int(max(np.sum(np.array(stats[c]["cnr"][z]) > 0), np.sum(np.array(stats[c]["cnr"][z]) < 0))))
                   for z in zones} for c in CH}
    for c in CH:
        print(f"{c:6s} noise {np.median(stats[c]['noise_sd_dB']):.2f} dB | " + " | ".join(
            f"{z}: CNR {summary[c][z]['cnr_median']:+.2f}, detected {summary[c][z]['detect_median']:.0f} %, same sign {summary[c][z]['sign_consistent']}/{n}"
            for z in zones))
    out = dict(pairs=labels, span_days=spans, channels=list(CH), summary=summary,
               per_pair={c: {k: v for k, v in s.items()} for c, s in stats.items()},
               note="CNR = median change in zone / robust std of change > 7 km from the vent; detection threshold flags 5 % of background")
    (OUT / "polarization.json").write_text(json.dumps(out, indent=1, ensure_ascii=False), encoding="utf-8")

    # ---- figure
    long_pairs = [i for i, s in enumerate(spans) if s >= 20] or list(range(n))
    fig = plt.figure(figsize=(15, 9.2))
    gs = fig.add_gridspec(2, 4, height_ratios=[1, 0.9], hspace=0.38, wspace=0.28)
    ext = [-HALF / 1000, HALF / 1000, -HALF / 1000, HALF / 1000]
    th = np.linspace(0, 2 * np.pi, 200)
    for k, c in enumerate(CH):
        a = fig.add_subplot(gs[0, k])
        m = np.nanmean(D[c][long_pairs], 0)
        lim = 2.0 if c != "HV/HH" else 1.5
        im = a.imshow(m, extent=ext, cmap="RdBu_r", vmin=-lim, vmax=lim, origin="upper")
        a.plot(np.cos(th), np.sin(th), "k-", lw=0.7)
        for r0 in (1.5, 6):
            a.plot(r0 * np.cos(th), r0 * np.sin(th), "k:", lw=0.5)
        for az in (195, 255):
            u = np.deg2rad(az)
            a.plot([1.5 * np.sin(u), 6 * np.sin(u)], [1.5 * np.cos(u), 6 * np.cos(u)], "k:", lw=0.5)
        a.set_title(f"{CH_UA[c]}\nсередня зміна, {len(long_pairs)} довгих пар (≥ 20 дн.)", fontsize=9)
        a.set_xlabel("схід від жерла, км", fontsize=8)
        if k == 0:
            a.set_ylabel("північ, км", fontsize=8)
            a.text(-4.6, -4.9, "SW:\nтефра", fontsize=8)
            a.text(0.0, 1.1, "кальдера", fontsize=8, ha="center")
        plt.colorbar(im, ax=a, fraction=0.046, pad=0.03, label="дБ")
    colors = {"HH": "tab:blue", "HV": "tab:green", "HV/HH": "tab:purple", "HH+HV": "0.45"}
    for j, (z, ttl) in enumerate((("caldera", "дно кальдери (< 1 км): нова лава"), ("SW", "сектор SW 1.5–6 км: тефра"),
                                   ("NE", "сектор NE 1.5–6 км: контроль (проти вітру)"))):
        a = fig.add_subplot(gs[1, j])
        x = np.arange(n)
        w = 0.2
        for k, c in enumerate(CH):
            a.bar(x + (k - 1.5) * w, stats[c]["cnr"][z], w, color=colors[c], label=c)
        a.axhline(0, color="k", lw=0.6)
        for s in (-2, 2):
            a.axhline(s, color="k", lw=0.4, ls=":")
        a.set_xticks(x)
        a.set_xticklabels([f"{l.split()[1]} {l.split()[3]}, {d} дн." for l, d in zip(labels, spans)], fontsize=6, rotation=90)
        a.set_title(ttl, fontsize=9)
        a.set_ylabel("зміна / шум фону (CNR)", fontsize=8)
        if j == 0:
            a.legend(fontsize=7, ncol=2)
    a = fig.add_subplot(gs[1, 3])
    a.axis("off")
    rows = [[c] + [f"{summary[c][z]['cnr_abs_median']:.1f} · {summary[c][z]['detect_median']:.0f} %" for z in zones] for c in CH]
    tb = a.table(cellText=rows, colLabels=["канал", "кальдера", "SW", "NE"], loc="center", cellLoc="center")
    tb.auto_set_font_size(False)
    tb.set_fontsize(8)
    tb.scale(1.25, 1.8)
    a.set_title("медіана |CNR| · частка виявлених клітинок\n(поріг: 5 % хибних тривог на фоні)", fontsize=9)
    fig.suptitle("Кілауеа 2026: яка поляризація NISAR краще показує тефру і нову лаву (GCOV, 100 м, 10 пар до/після)", fontsize=12)
    fig.savefig(OUT / "polarization.png", dpi=105, facecolor="white", bbox_inches="tight")
    print("wrote", OUT / "polarization.png")


if __name__ == "__main__":
    main()
