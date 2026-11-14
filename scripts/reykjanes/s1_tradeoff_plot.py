"""H1 / H8 figure: does GNSS remove the trade-off between the dike volume and the Svartsengi source?

Reads the depth profiles of s1_dike_distributed.py (InSAR + GNSS with the layered FE source; InSAR only; InSAR + GNSS
with Mogi) and, if present, the NPE v2 / MCMC posteriors (s1_npe_v2.npz, s1_npe_v2_mcmc_*.npy).
Writes output/reykjanes/ml/s1_tradeoff.png.
"""
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import s1_dike_variable as V
from s1_grindavik_dike import ML

RUNS = {"": ("InSAR + GNSS, шарувата кора (МСЕ)", "tab:red", "-"), "_nognss": ("лише InSAR, МСЕ", "tab:blue", "-"),
        "_mogi": ("InSAR + GNSS, Моґі", "tab:red", "--")}


def main():
    fig, axs = plt.subplots(1, 4, figsize=(19, 4.8))
    for tag, (lab, col, ls) in RUNS.items():
        f = ML / f"s1_dike_distributed{tag}.json"
        if not f.exists():
            continue
        d = json.loads(f.read_text(encoding="utf-8"))
        P = d["depth_profile"]
        z = np.array([p["depth"] for p in P])
        c2 = np.array([p["chi2"] for p in P])
        axs[0].plot(z, c2 - c2.min(), ls, color=col, marker="o", ms=3, label=lab)
        axs[1].plot(z, [p["V"] for p in P], ls, color=col, marker="o", ms=3, label=lab)
        axs[2].plot(z, [p["dV_A016"] for p in P], ls, color=col, marker="o", ms=3, label=lab)
    axs[0].axhline(1, color="0.5", lw=.7, ls=":"); axs[0].axhline(4, color="0.5", lw=.7, ls=":")
    axs[0].set_ylim(0, 60)
    axs[0].set_ylabel("Δχ² від найкращої глибини", fontsize=9)
    axs[0].set_title("наскільки дані визначають глибину Свартсенгі\n(Δχ² = 1 і 4 — межі 1σ і 2σ)", fontsize=9)
    axs[1].axhspan(130, 139, color="0.85", label="опубліковано 130–139")
    axs[1].set_ylabel("об’єм дайки, млн м³", fontsize=9)
    axs[1].set_title("об’єм дайки за кожної глибини джерела", fontsize=9)
    axs[2].axhspan(-82, -76, color="0.85", label="опубліковано −76…−82")
    axs[2].set_ylabel("dV Свартсенгі (вікно A016), млн м³", fontsize=9)
    axs[2].set_title("об’єм джерела за кожної глибини", fontsize=9)
    for a in axs[:3]:
        a.axvspan(4.0, 4.2, color="gold", alpha=.35, lw=0)
        a.set_xlabel("глибина джерела Свартсенгі, км (жовте: стаття, Моґі)", fontsize=8)
        a.grid(alpha=.3); a.tick_params(labelsize=8)
    axs[1].legend(fontsize=7)
    # NPE v2 / MCMC joint posterior V_dike vs source depth
    a = axs[3]
    have = False
    for nm, f, col in (("NPE v2 (InSAR + GNSS)", ML / "s1_npe_v2.npz", "tab:purple"),
                       ("MCMC з GNSS", ML / "s1_npe_v2_mcmc_gnss.npy", "tab:red"), ("MCMC без GNSS", ML / "s1_npe_v2_mcmc_insar_only.npy", "tab:blue")):
        if not f.exists():
            continue
        s = np.load(f)["ps"] if f.suffix == ".npz" else np.load(f)
        s = s[np.random.default_rng(0).choice(len(s), min(1500, len(s)), replace=False)]
        a.scatter(s[:, 13], V.volume(s), s=3, alpha=.25, color=col, label=nm)
        have = True
    a.axhspan(130, 139, color="0.85"); a.axvspan(4.0, 4.2, color="gold", alpha=.35, lw=0)
    a.set_xlabel("глибина джерела Свартсенгі, км", fontsize=8); a.set_ylabel("об’єм дайки, млн м³", fontsize=8)
    a.set_title("апостеріорні зразки: об’єм дайки ↔ глибина джерела" if have else "NPE v2 / MCMC: ще не пораховано", fontsize=9)
    if have:
        a.legend(fontsize=7, markerscale=4)
    a.grid(alpha=.3); a.tick_params(labelsize=8)
    fig.suptitle("Гіпотеза H1: компроміс «об’єм дайки ↔ джерело Свартсенгі» з GNSS і без; H8: МСЕ проти Моґі", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(ML / "s1_tradeoff.png", dpi=100, facecolor="white")
    print("wrote", ML / "s1_tradeoff.png")


if __name__ == "__main__":
    main()
