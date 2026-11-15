"""WS: one summary figure of the Wisconsin story tests (reads data/processed/wisconsin/ws/ws_results.json).

uv run python experiments/ws_plot.py      -> data/processed/wisconsin/ws/ws_tests.png
"""
import json
import sys
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src"), str(ROOT / "experiments")]
RES = ROOT / "data/processed/wisconsin/ws"
D = lambda s: datetime.fromisoformat(s[:10])  # noqa: E731


def main():
    r = json.loads((RES / "ws_results.json").read_text(encoding="utf-8"))
    fig, ax = plt.subplots(2, 3, figsize=(18, 10))
    # WS9
    a = ax[0, 0]
    p = r["WS9"]["points"]
    for d, m, c in (("A", "o", "C0"), ("D", "s", "C3")):
        q = [x for x in p if x["track"] in (76, 148)] if d == "A" else [x for x in p if x["track"] in (41, 113)]
        a.scatter([x["omega_pred"] for x in q], [x["omega"] for x in q], marker=m, color=c, s=22,
                  label="ascending (dawn)" if d == "A" else "descending (dusk)")
    lim = max([abs(x["omega_pred"]) for x in p] + [abs(x["omega"]) for x in p] + [1]) * 1.1
    a.plot([-lim, lim], [-lim, lim], "k:", lw=0.8)
    a.axhline(0, color="0.85", lw=0.6), a.axvline(0, color="0.85", lw=0.6)
    a.set(xlabel="Ω predicted from CODE TEC × B∥, deg", ylabel="Ω measured (Bickel–Bates), deg",
          title=f"WS9 Faraday: ρ(|Ω|,|Ωp|) = {r['WS9'].get('rho_abs', float('nan')):.2f}, sign agree "
                f"{(r['WS9'].get('sign_share') or 0):.0%}")
    a.legend(fontsize=8)
    # WS1
    a = ax[0, 1]
    rows = [x for x in r["rows_A"] if x.get("corn_hh_vv_db") is not None and x.get("soy_hh_vv_db") is not None]
    a.plot([D(x["date"]) for x in rows], [x["corn_hh_vv_db"] - x["soy_hh_vv_db"] for x in rows], "o", label="HH/VV")
    a.plot([D(x["date"]) for x in rows], [x["corn_hv"] - x["soy_hv"] for x in rows], "s", label="HV")
    a.axhline(1, color="0.6", ls="--", lw=0.8), a.axhline(0, color="0.85", lw=0.6)
    a.axvspan(D("2025-10-22"), D("2025-11-30"), color="0.93", zorder=0)
    a.set(ylabel="corn − soybean, dB", title=f"WS1 corn stalks (shaded = judged): a {r['WS1']['WS1a']}, b {r['WS1']['WS1b']}")
    a.legend(fontsize=8), a.tick_params(axis="x", rotation=30)
    # WS3
    a = ax[0, 2]
    s = [x for x in r["WS3"].get("series", []) if x[2] is not None]
    for t in sorted({x[1] for x in s}):
        q = [x for x in s if x[1] == t]
        a.plot([D(x[0]) for x in q], [x[2] for x in q], "o-", ms=4, label=f"track {t}")
    for k, c in (("open_date", "C2"), ("ice_date", "C0")):
        if r["WS3"].get(k):
            a.axvline(D(r["WS3"][k]), color=c, ls="--", lw=1, label=k.replace("_", " ") + " (S2)")
    a.set(ylabel="lake HH, dB", title=f"WS3 lake ice: {r['WS3']['verdict']}")
    a.legend(fontsize=7), a.tick_params(axis="x", rotation=30)
    # WS4
    a = ax[1, 0]
    lags = [x["lag_days"] for x in r["WS4"]["pairs"]]
    if lags:
        a.hist(lags, bins=np.arange(-90, 91, 12), color="C1")
    a.axvspan(-24, 24, color="0.92", zorder=0)
    a.set(xlabel="radar harvest date − optical harvest date, days", ylabel="(block, track) pairs",
          title=f"WS4 harvest: {r['WS4']['n_pairs']} pairs, within ±24 d "
                f"{(r['WS4'].get('share_within_24d') or 0):.0%} → {r['WS4']['verdict']}")
    # WS6
    a = ax[1, 1]
    g = [x for x in r["rows_gunw"] if x["days"] == 12 and x.get("upland_coh") is not None]
    col = {"frozen": "C0", "unfrozen": "C3"}
    for x in g:
        st = x["state"]
        c = col.get(st[0]) if st[0] == st[1] else "0.6"
        a.plot(D(x["ref"]), x["upland_coh"], "o", color=c or "0.6", ms=4)
        a.plot(D(x["ref"]), x["crop_coh"] if x.get("crop_coh") is not None else np.nan, "^", color=c or "0.6", ms=4)
    a.set(ylabel="12-day coherence (o forest, ^ cropland)", title=f"WS6 frozen (blue) vs unfrozen (red): {r['WS6']['verdict']}")
    a.tick_params(axis="x", rotation=30)
    # WS5 + WS8
    a = ax[1, 2]
    pts = r["WS5"].get("points", [])
    a.scatter([x[0] for x in pts], [x[1] for x in pts], color="C4", s=22)
    a.set(xlabel="rain in the 72 h before the pass, mm", ylabel="cropland HH anomaly, dB",
          title=f"WS5 rain on bare fields: ρ = {r['WS5'].get('rho', float('nan')):.2f} → {r['WS5']['verdict']}")
    fig.suptitle("Northern Wisconsin, NISAR L-band: story tests WS1–WS9 (pre-registered, pushed before the data)")
    fig.tight_layout()
    fig.savefig(RES / "ws_tests.png", dpi=110)
    print("figure:", RES / "ws_tests.png")


if __name__ == "__main__":
    main()
