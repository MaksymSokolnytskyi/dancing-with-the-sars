"""SBC, PPC and prior-sensitivity figures for NPE v2 (reads s1_npe_v2.json / .npz, s1_npe_v2_wide.npz).
Writes output/reykjanes/ml/s1_sbc.png, s1_ppc.png, s1_prior_sensitivity.png.
"""
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import s1_dike_variable as V
from s1_grindavik_dike import ML

UA = {"dike east km": "схід дайки", "dike north km": "північ дайки", "strike deg": "простягання", "length km": "довжина",
      "top km": "верх", "height km": "висота", "src east km": "схід джерела", "src north km": "північ джерела",
      "src depth km": "глибина джерела", "src dV A016 1e6 m3": "dV джерела A016", "src dV D155 1e6 m3": "dV джерела D155"}


def ua(nm):
    return UA.get(nm, nm.replace("log10 opening node", "розкриття, вузол"))


def main():
    import sys
    v3 = "--v3" in sys.argv
    R = json.loads((ML / ("s1_npe_v3.json" if v3 else "s1_npe_v2.json")).read_text(encoding="utf-8"))
    sb = R["sbc"]
    ranks = np.array(sb["ranks"])
    names = list(sb["p_uniform"])
    fig, axs = plt.subplots(3, 6, figsize=(19, 8.5))
    for a, i in zip(axs.ravel(), list(range(len(names))) + [-1]):
        r, nm, p = (ranks[:, i], ua(names[i]), sb["p_uniform"][names[i]]) if i >= 0 else (np.array(sb["ranks_volume"]), "об’єм дайки", sb["p_uniform_volume"])
        a.hist(r, bins=20, range=(0, 1000), color="tab:blue" if p > 0.01 else "tab:red", alpha=.8)
        a.axhline(len(r) / 20, color="k", lw=.8, ls="--")
        a.set_title(f"{nm}\np = {p:.3f}", fontsize=8, color="k" if p > 0.01 else "firebrick")
        a.set_xticks([]); a.tick_params(labelsize=7)
    a = axs.ravel()[-1]
    a.cla()
    lev = np.array(sb["levels"])
    cov = np.array(sb["coverage"])
    for j in range(cov.shape[1]):
        a.plot(lev, cov[:, j], color="0.6", lw=.7)
    a.plot([0, 1], [0, 1], "k--", lw=.8)
    a.set_xlabel("номінальний рівень", fontsize=7); a.set_ylabel("частка правди в інтервалі", fontsize=7)
    a.set_title("покриття (усі параметри)", fontsize=8); a.tick_params(labelsize=7)
    fig.suptitle(f"SBC: ранги правди серед 1000 зразків NPE {'v3 (випадковий рівень шуму)' if v3 else 'v2'} на {sb['n']} симуляціях (25 % без дайки); рівномірна гістограма = калібрування. "
                 f"Червоне — p < 0.01. Без дайки: медіана P(без дайки) = {sb['nodike']['p_nodike_median']:.2f}", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(ML / ("s1_sbc_v3.png" if v3 else "s1_sbc.png"), dpi=95, facecolor="white")
    if v3:
        print("wrote SBC v3"); return
    # PPC
    pp = R["ppc"]
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.2))
    a = axs[0]
    labels = ["NPE v2 (зразки)", "MCMC з GNSS (медіана)", "MCMC без GNSS (медіана)", "очікування для шуму"]
    vals = [pp["chi2_median"], -2 * R["mcmc"]["gnss"]["loglike_median"], -2 * R["mcmc"]["insar_only"]["loglike_median"], pp["n_data"]]
    a.barh(labels, vals, color=["tab:purple", "tab:red", "tab:blue", "0.6"])
    for k, v in enumerate(vals):
        a.text(v * 1.02, k, f"{v:.0f}", va="center", fontsize=8)
    a.set_xlabel("χ² (вибілені залишки)", fontsize=8)
    a.set_title("наскільки модель пояснює дані", fontsize=9)
    a.tick_params(labelsize=8)
    a = axs[1]
    a.barh(["залишок NPE уздовж\nнапрямків дайки", "чистий шум, 95 %"], [pp["dike_projection_median"], pp["dike_projection_noise_95"]], color=["tab:purple", "0.6"])
    a.set_xlabel("норма проєкції вибіленого залишку на якобіан дайки", fontsize=8)
    a.set_title("системні «крила» в залишках (PPC)", fontsize=9)
    a.tick_params(labelsize=8)
    fig.suptitle("Перевірка за апостеріорним передбаченням (PPC) для 10.11.2023", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(ML / "s1_ppc.png", dpi=100, facecolor="white")
    # prior sensitivity
    if (ML / "s1_npe_v2_wide.npz").exists():
        s, w = np.load(ML / "s1_npe_v2.npz")["ps"], np.load(ML / "s1_npe_v2_wide.npz")["ps"]
        W = json.loads((ML / "s1_npe_v2_wide.json").read_text(encoding="utf-8"))
        fig, axs = plt.subplots(1, 4, figsize=(17, 4))
        for a, (i, nm) in zip(axs[:3], ((2, "простягання, °"), (0, "схід центру, км"), (1, "північ центру, км"))):
            a.hist(s[:, i], bins=50, density=True, histtype="step", lw=1.8, color="tab:purple", label="стандартне апріорне (±4 км, ±30°)")
            a.hist(w[:, i], bins=50, density=True, histtype="step", lw=1.8, color="tab:orange", label="широке (±15 км, 0–180°)")
            a.set_xlabel(nm, fontsize=8); a.tick_params(labelsize=7)
        axs[0].legend(fontsize=7)
        a = axs[3]
        a.hist(np.clip(V.volume(s[:3000]), 0.5, 500), bins=np.logspace(-0.3, 2.7, 50), histtype="step", lw=1.8, color="tab:purple")
        a.hist(np.clip(V.volume(w[:3000]), 0.5, 500), bins=np.logspace(-0.3, 2.7, 50), histtype="step", lw=1.8, color="tab:orange")
        a.set_xscale("log"); a.set_xlabel("об’єм дайки, млн м³", fontsize=8); a.tick_params(labelsize=7)
        fig.suptitle(f"Чутливість до апріорного розподілу (H6): широке — {W['n_sims'] // 1000} тис. симуляцій, менша мережа", fontsize=11)
        fig.tight_layout(rect=(0, 0, 1, 0.92))
        fig.savefig(ML / "s1_prior_sensitivity.png", dpi=100, facecolor="white")
    print("wrote SBC / PPC / prior-sensitivity figures")


if __name__ == "__main__":
    main()
