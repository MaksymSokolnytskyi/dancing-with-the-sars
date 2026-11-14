"""Apply the trained NPE v2 (no retraining) to other Sentinel-1 pairs: null tests (H4) and 2024 events (H7).

  quiet_2023      7/5 Oct -> 19/17 Oct 2023: no dike, no inflation  -> expect P(no dike) high
  inflation_2024  21/19 Jul -> 2 Aug / 31 Jul 2024: Svartsengi inflation only -> expect no dike, source dV > 0
  event_20240208, event_20240302 (dike without eruption), event_20240316, event_20240529: dikes on the Sundhnukur row

Each pair is binned onto the SAME 1 km cells (s1_common_v2.pair_vector; empty cells interpolated, share reported)
with the same graben mask (the 2023 trace; the 2024 dikes run along the same crater row); GNSS for the same
stations / windows (s1_common_v2.gnss_vector; a missing value is set to 0 after whitening = "no information").
Hot vents from the HLS time-lapse cache (B12 > 0.25 and B12 > 1.2 x B11) for the eruptions of 16 Mar
(scene 2024-03-30) and 29 May 2024 (scene 2024-06-10) -> distance to the NPE dike trace.

Writes output/reykjanes/ml/s1_npe_v2_apply.json and s1_events_2024.png, s1_null_tests.png.
"""
import json
from datetime import datetime

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import torch
from pyproj import Transformer

import s1_dike_variable as V
import s1_grindavik_dike as S
from ml_npe_flow import Net, from_z
from reykjanes_analysis import PLACES
from s1_common_v2 import GeomV2, gnss_vector, pair_vector, product_dir

ML = S.ML
import sys
V3 = "--v3" in sys.argv
TAG = "_v3" if V3 else ""
HLS_VENT_SCENES = {"event_20240316": "HLS.L30.T27VVL.2024090", "event_20240529": "HLS.S30.T27VVL.2024162"}


def load_net(G):
    ck = torch.load(ML / f"s1_npe{TAG or '_v2'}.pt", weights_only=False)
    net = Net(ck["U"].shape[1] + G.n_gnss, len(ck["lo"]))
    net.load_state_dict(ck["state"])
    net.eval()
    return net, ck


def posterior(net, ck, G, y, n=6000):
    xw = G.whiten(y[None])
    xw[~np.isfinite(xw)] = 0.0
    f = np.c_[xw[:, :ck["ni"]] @ ck["U"], xw[:, ck["ni"]:]] / ck["scale"]
    with torch.no_grad():
        z = net.flow(net.emb(torch.tensor(f, dtype=torch.float32))).sample((n,)).numpy()[:, 0].astype(float)
    return from_z(z, ck["lo"], ck["hi"])


def vents(prefix):
    import glob
    f = sorted(glob.glob(str(S.ROOT / "data" / "hls" / "timelapse" / f"{prefix}*.npz")))
    if not f:
        return []
    z = np.load(f[0])
    if "swir2" not in z:
        return []
    b12, b11 = z["swir2"] * 1e-4, z["swir1"] * 1e-4
    hot = (b12 > 0.25) & (b12 > 1.2 * b11)
    from scipy import ndimage
    lab, nl = ndimage.label(ndimage.binary_dilation(hot, iterations=2))
    ny, nx = b12.shape
    E, N = np.meshgrid(-8 + (np.arange(nx) + .5) * 16 / nx, 7 - (np.arange(ny) + .5) * 14 / ny)
    out = []
    for i in range(1, nl + 1):
        k = (lab == i) & hot
        if k.sum() >= 2:
            out.append((float(E[k].mean()), float(N[k].mean()), int(k.sum())))
    return sorted(out, key=lambda v: -v[2])[:3]


def seg_dist(px, py, p):
    u = np.array([np.sin(np.deg2rad(p[2])), np.cos(np.deg2rad(p[2]))])
    x0, y0, x1, y1 = p[0] - u[0] * p[3] / 2, p[1] - u[1] * p[3] / 2, p[0] + u[0] * p[3] / 2, p[1] + u[1] * p[3] / 2
    vx, vy = x1 - x0, y1 - y0
    t = np.clip(((px - x0) * vx + (py - y0) * vy) / (vx * vx + vy * vy), 0, 1)
    return float(np.hypot(px - (x0 + t * vx), py - (y0 + t * vy)))


def main():
    G = GeomV2(gnss=True, noise="ref" if V3 else "auto")
    net, ck = load_net(G)
    idx = json.loads((S.S1 / "pairs_index.json").read_text())
    groups = {"real_20231110": None, **{k: v for k, v in idx.items()}}
    out, post = {}, {}
    for name, pairs in groups.items():
        if pairs is None:
            y, filled = G.y.copy(), 0.0
        else:
            try:
                dirs = {p["track"]: product_dir(p["track"], p["ref"], p["sec"]) for p in pairs}
            except FileNotFoundError as ex:
                print("skip", name, ex); continue
            yi, filled = pair_vector(G, dirs)
            dates = {p["track"]: (datetime.strptime(p["ref"][17:25], "%Y%m%d").date().isoformat(),
                                  datetime.strptime(p["sec"][17:25], "%Y%m%d").date().isoformat()) for p in pairs}
            y = np.r_[yi, gnss_vector(G, dates)]
        ps = posterior(net, ck, G, y)
        Vd = V.volume(ps[:3000])
        pm = np.median(ps, 0)
        pred = V.forward(G, pm)
        r = G.whiten((y - pred)[None])
        r[~np.isfinite(r)] = 0
        o = dict(filled_fraction=round(float(filled), 3), p_nodike=round(float(np.mean(Vd < 5)), 3),
                 volume=np.percentile(Vd, [16, 50, 84]).round(1).tolist(),
                 pct={nm: np.percentile(ps[:, i], [16, 50, 84]).round(2).tolist() for i, nm in enumerate(V.NAMES)},
                 chi2=round(float((r ** 2).sum()), 1), n_data=int(r.shape[1]),
                 gnss_missing=int(np.sum(~np.isfinite(y[G.n_ins:]))))
        if name in HLS_VENT_SCENES:
            vv = vents(HLS_VENT_SCENES[name])
            o["hls_vents"] = vv
            o["vent_distance_km"] = [round(seg_dist(a, b, pm), 2) for a, b, _ in vv]
        out[name], post[name] = o, ps[:2000]
        print(f"{name:16s} filled {filled:.0%}; P(no dike) {o['p_nodike']:.2f}; V {o['volume']} 1e6 m3; centre "
              f"({pm[0]:.1f}, {pm[1]:.1f}) km, strike {pm[2]:.0f}, length {pm[3]:.1f}; src dV {o['pct']['src dV A016 1e6 m3'][1]:.0f}/"
              f"{o['pct']['src dV D155 1e6 m3'][1]:.0f}; chi2 {o['chi2']:.0f}/{o['n_data']}" +
              (f"; vents {o.get('vent_distance_km')}" if "vent_distance_km" in o else ""), flush=True)
    (ML / f"s1_npe_v2_apply{TAG}.json").write_text(json.dumps(out, indent=1, ensure_ascii=False), encoding="utf-8")
    np.savez_compressed(ML / f"s1_npe_v2_apply{TAG}.npz", **post)
    plot(out, post)


def plot(out, post):
    to_utm = Transformer.from_crs(4326, 32627, always_xy=True)
    sx, sy = to_utm.transform(*PLACES["Svartsengi"])
    pl = {k: ((to_utm.transform(*PLACES[k])[0] - sx) / 1000, (to_utm.transform(*PLACES[k])[1] - sy) / 1000)
          for k in ("Svartsengi", "Grindavik", "Sundhnukur")}
    lab = {"real_20231110": "10.11.2023 (навчальна подія)", "quiet_2023": "тиша, жовт. 2023", "inflation_2024": "лише інфляція, лип. 2024",
           "event_20240208": "08.02.2024", "event_20240302": "02.03.2024 (без виверження)", "event_20240316": "16.03.2024", "event_20240529": "29.05.2024"}
    # null tests: volume distributions
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.4))
    for k, col in (("real_20231110", "tab:red"), ("quiet_2023", "tab:green"), ("inflation_2024", "tab:blue")):
        if k in post:
            v = V.volume(post[k])
            ax[0].hist(np.clip(v, 0.1, 400), bins=np.logspace(-1, np.log10(400), 50), histtype="step", lw=1.8, color=col,
                       label=f"{lab[k]}: P(без дайки) {out[k]['p_nodike']:.2f}")
            ax[1].hist(post[k][:, 14], bins=40, histtype="step", lw=1.8, color=col, label=lab[k])
    ax[0].set_xscale("log"); ax[0].axvline(5, color="k", ls=":", lw=.8)
    ax[0].set_xlabel("об’єм дайки, млн м³ (лог.)"); ax[0].set_title("нульовий тест: чи «бачить» NPE дайку там, де її немає", fontsize=10)
    ax[0].legend(fontsize=8)
    ax[1].axvline(0, color="k", lw=.6)
    ax[1].set_xlabel("dV Свартсенгі, трек A016, млн м³"); ax[1].set_title("джерело: інфляція має дати dV > 0", fontsize=10)
    ax[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(ML / f"s1_null_tests{TAG}.png", dpi=100, facecolor="white")
    # events 2024: map of traces
    fig, ax = plt.subplots(1, 2, figsize=(14, 6), gridspec_kw=dict(width_ratios=[1, 1.15]))
    a = ax[0]
    cols = {"real_20231110": "k", "event_20240208": "tab:blue", "event_20240302": "tab:green", "event_20240316": "tab:orange", "event_20240529": "tab:red"}
    for k, c in cols.items():
        if k not in post:
            continue
        for p in post[k][:120]:
            u = np.array([np.sin(np.deg2rad(p[2])), np.cos(np.deg2rad(p[2]))])
            a.plot([p[0] - u[0] * p[3] / 2, p[0] + u[0] * p[3] / 2], [p[1] - u[1] * p[3] / 2, p[1] + u[1] * p[3] / 2], color=c, alpha=.04)
        pm = np.median(post[k], 0)
        u = np.array([np.sin(np.deg2rad(pm[2])), np.cos(np.deg2rad(pm[2]))])
        a.plot([pm[0] - u[0] * pm[3] / 2, pm[0] + u[0] * pm[3] / 2], [pm[1] - u[1] * pm[3] / 2, pm[1] + u[1] * pm[3] / 2], color=c, lw=2.2, label=lab[k])
        for v in out[k].get("hls_vents", []):
            a.plot(v[0], v[1], "*", color=c, mec="k", ms=13)
    for k, (x0, y0) in pl.items():
        a.plot(x0, y0, "k^" if k == "Svartsengi" else "ko", ms=5)
        a.text(x0 + .4, y0 + (.6 if k == "Sundhnukur" else -.8 if k == "Svartsengi" else 0), {"Svartsengi": "Свартсенгі", "Grindavik": "Ґріндавік", "Sundhnukur": "Сундгнукур"}[k], fontsize=8)
    a.set_aspect("equal"); a.set_xlim(-10, 12); a.set_ylim(-11, 11); a.grid(alpha=.3)
    a.legend(fontsize=8, loc="lower right"); a.set_title(f"сліди дайок за NPE {'v3' if V3 else 'v2'} без перенавчання (★ — гарячі жерла HLS)", fontsize=10)
    a.set_xlabel("схід від Свартсенгі, км"); a.set_ylabel("північ, км")
    a = ax[1]
    a.axis("off")
    rows = []
    for k in ("real_20231110", "event_20240208", "event_20240302", "event_20240316", "event_20240529", "quiet_2023", "inflation_2024"):
        if k not in out:
            continue
        o = out[k]
        f = lambda v: f"{v[1]:.0f} [{v[0]:.0f}–{v[2]:.0f}]"
        rows.append([lab[k], f(o["volume"]), f"{o['p_nodike']:.2f}", f"{o['pct']['length km'][1]:.0f}", f"{o['pct']['strike deg'][1]:.0f}°",
                     f"{o['pct']['dike north km'][1]:+.1f}", f"{o['pct']['src dV A016 1e6 m3'][1]:.0f}", ", ".join(f"{d:.1f}" for d in o.get("vent_distance_km", [])) or "—"])
    tb = a.table(cellText=rows, colLabels=["подія", "об’єм,\nмлн м³", "P(без\nдайки)", "довж.,\nкм", "прост.", "центр\nN, км", "dV Св.\nA016", "жерла,\nкм"],
                 loc="center", cellLoc="center", colWidths=[.27, .14, .09, .08, .08, .08, .09, .14])
    tb.auto_set_font_size(False); tb.set_fontsize(7.5); tb.scale(1, 1.9)
    for (r, c), cell in tb.get_celld().items():
        if r == 0:
            cell.set_height(cell.get_height() * 1.5)
    a.set_title(f"NPE {'v3 (рівень шуму невідомий)' if V3 else 'v2'} на подіях, яких мережа не бачила", fontsize=10)
    fig.tight_layout()
    fig.savefig(ML / f"s1_events_2024{TAG}.png", dpi=100, facecolor="white")
    print("wrote null-test and events figures", TAG)


if __name__ == "__main__":
    main()
