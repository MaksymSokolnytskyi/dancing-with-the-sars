"""Build output/reykjanes/earth_breathes.html: "Земля дихає" web map with a time slider,
layer toggles, per-cell trust, time-growing hachures, a cross-section tool, a headline
number with its rule and an independent GNSS check.

Inputs: output/reykjanes/breathing_data.npz, breathing_vectors.json (build_breathing_data.py)
and map_layers.json (build_reykjanes_report.py).

Usage: python scripts/reykjanes/build_breathing_app.py
"""
import base64
import json
from pathlib import Path

import numpy as np
from pyproj import Transformer

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "output" / "reykjanes"
BOX = (-22.80, 63.76, -22.05, 64.02)
NLON, NLAT = 420, 290


def b64(a):
    return base64.b64encode(np.ascontiguousarray(a).tobytes()).decode()


def mm16(a):
    return np.where(np.isfinite(a), np.clip(np.round(a * 10), -32000, 32000), -32768).astype(np.int16)  # cm -> mm


ML = OUT / "ml"
SVARTSENGI = (-22.433, 63.878)


def local_to_lonlat(e, n):
    """Local east/north km around Svartsengi (the frame of the source models) -> lon, lat."""
    tr = Transformer.from_crs(4326, 3413, always_xy=True)
    sx, sy = tr.transform(*SVARTSENGI)
    lam = np.deg2rad(SVARTSENGI[0] + 45.0)
    east, north = np.array([np.cos(lam), np.sin(lam)]), np.array([-np.sin(lam), np.cos(lam)])
    X = sx + 1000 * (e * east[0] + n * north[0])
    Y = sy + 1000 * (e * east[1] + n * north[1])
    return tr.transform(X, Y, direction="INVERSE")


def calibration_from_log(path):
    """{label: {"68": [...], "95": [...]}} parsed from an NPE log."""
    import re
    out = {}
    if not path.exists():
        return out
    for line in path.read_text(encoding="utf-8").splitlines():
        m = re.match(r"calibration \[(.+?)\]", line)
        if not m:
            continue
        parts = re.split(r"95 %", line)
        p68 = [int(v) for v in re.findall(r"(\d+)%", parts[0].split("68 %")[1])]
        p95 = [int(v) for v in re.findall(r"(\d+)%", parts[1])]
        out[m.group(1)] = {"68": p68, "95": p95}
    return out


def npe_block(path, names, logpath, n_keep=3000, prefix="", cal_json=None):
    """Posterior samples + summary per data set. prefix selects one variant inside a multi-variant
    file (e.g. "raw__"); cal_json gives calibration from a json written by the screen experiment."""
    if not path.exists():
        return None
    z = np.load(path)
    cal = calibration_from_log(logpath) if logpath else {}
    rng = np.random.default_rng(0)
    out = {}
    for key, label in (("NISAR_and_GNSS", "NISAR + GNSS"), ("NISAR_only", "NISAR only")):
        key = prefix + key
        if key not in z:
            continue
        s = z[key]
        s = s[rng.choice(len(s), min(n_keep, len(s)), replace=False)]
        lon, lat = local_to_lonlat(s[:, 0], s[:, 1])
        summ = {nm: np.round(np.percentile(s[:, i], [2.5, 16, 50, 84, 97.5]), 2).tolist() for i, nm in enumerate(names)}
        block = dict(samples=np.round(s, 3).tolist(), lon=np.round(lon, 5).tolist(), lat=np.round(lat, 5).tolist(),
                     summary=summ, calibration=cal.get(label))
        if cal_json is not None and label in cal_json:
            cj = cal_json[label]
            if "cal_on_raw" in cj:
                block["calibration_real"] = {"68": [round(v * 100) for v in cj["cal_on_raw"]["0.68"]],
                                             "95": [round(v * 100) for v in cj["cal_on_raw"]["0.95"]]}
            if "rhat" in cj:
                block["mcmc"] = dict(rhat=max(cj["rhat"]), ess=int(min(cj["ess"])), n_forward=cj["n_forward"],
                                     seconds=round(cj["seconds"]))
        if len(names) == 7:  # elongated body: end points of the long axis for a subset of samples
            k = s[:300]
            a, strike = k[:, 3], np.deg2rad(k[:, 4])
            e1, n1 = k[:, 0] + a * np.sin(strike), k[:, 1] + a * np.cos(strike)
            e2, n2 = k[:, 0] - a * np.sin(strike), k[:, 1] - a * np.cos(strike)
            l1, b1 = local_to_lonlat(e1, n1)
            l2, b2 = local_to_lonlat(e2, n2)
            block["axes"] = np.round(np.c_[b1, l1, b2, l2], 5).tolist()
            block["p_long"] = round(float(np.mean(s[:, 3] > 1.5)), 3)
        out[label] = block
    return out


def ml_payload(d, ny, nx):
    """GP, change-point and NPE results packed for the page (everything optional)."""
    meta = {}
    res = dict(gp_clean=None, gp_screen=None, gp_smooth=None, cp_step=None, cp_hinge=None)
    if (ML / "gp_series.npz").exists():
        g = np.load(ML / "gp_series.npz")
        res["gp_clean"] = b64(mm16(g["clean"] - g["clean"][:1]))
        res["gp_smooth"] = b64(mm16(g["smooth"] - g["smooth"][:1]))
        res["gp_screen"] = b64(mm16(g["screens"]))
        h = g["hyper"]
        meta["gp"] = dict(sigma_f=round(float(h[6]), 2), length_days=round(float(h[7]) * 365.25), white=round(float(h[8]), 2),
                          screen_sigma=round(float(h[3]), 2), screen_len_km=round(float(h[4]), 1),
                          screen_std=[round(float(np.nanstd(s)), 2) for s in g["screens"]])
        # GP-cleaned station rates (vertical), for the validation table
        lon, lat = d["lon"], d["lat"]
        t = g["t"]
        vec = json.loads((OUT / "breathing_vectors.json").read_text(encoding="utf-8"))
        gr = {}
        for s, gg in vec["gnss"].items():
            i = int(np.argmin((lat - gg["lat"]) ** 2 + ((lon - gg["lon"]) * 0.44) ** 2))
            ii, jj = np.unravel_index(i, (ny, nx))
            y = np.nanmedian(g["clean"][:, max(ii - 1, 0):ii + 2, max(jj - 1, 0):jj + 2], axis=(1, 2))
            ok = np.isfinite(y)
            gr[s] = round(float(np.polyfit(t[ok], y[ok], 1)[0]), 1)
        meta["gp"]["station_rate"] = gr
    if (ML / "changepoint.npz").exists():
        c = np.load(ML / "changepoint.npz")
        thr = np.log(20.0)
        res["cp_step"] = b64(mm16(np.where(c["step_lbf"] > thr, c["step_size"], np.nan)))
        res["cp_hinge"] = b64(mm16(np.where(c["hinge_lbf"] > thr, c["hinge_size"], np.nan)))
        tests = [dict(amp=float(x[0]), lbf=round(float(x[1]), 1), date=str(d["dates"][int(x[2])]), size=round(float(x[3]), 1),
                      found=bool(x[4]), false_alarm=round(float(x[5]) * 100, 1)) for x in c["tests"]]
        fin = np.isfinite(c["step_lbf"])
        meta["cp"] = dict(tests=tests, inj_date=str(d["dates"][int(c["inj"])]),
                          step_share=round(float(np.mean((c["step_lbf"] > thr)[fin]) * 100), 1),
                          hinge_share=round(float(np.mean((c["hinge_lbf"] > thr)[fin]) * 100), 1))
    meta["npe"] = dict(
        mogi=npe_block(ML / "npe_posterior.npz", ("east", "north", "depth", "volume"), ML / "npe_log.txt"),
        ellipsoid=npe_block(ML / "npe_ellipsoid.npz", ("east", "north", "depth", "half_length", "strike", "volume", "sill_share"),
                            ML / "npe_ellipsoid_log.txt"))
    meta["methods"] = methods_payload()
    tj = ML / "timings.json"
    meta["timings"] = flow_timings(json.loads(tj.read_text(encoding="utf-8")) if tj.exists() else None)
    meta["expensive"] = expensive_payload()
    res["meta"] = meta
    return res


def methods_payload():
    """Same physics, different inference / noise: NPE (Gaussian noise), MCMC (same likelihood),
    MCMC with the real-screen covariance, NPE trained on raw and on GP real screens."""
    names = dict(mogi=("east", "north", "depth", "volume"),
                 ellipsoid=("east", "north", "depth", "half_length", "strike", "volume", "sill_share"))
    out = {}
    for kind in ("mogi", "ellipsoid"):
        m = {}
        for meth, fname in (("mcmc", f"mcmc_{kind}"), ("mcmc_screens", f"mcmc_{kind}_screens")):
            j = ML / f"{fname}.json"
            m[meth] = npe_block(ML / f"{fname}.npz", names[kind], None, 2000,
                                cal_json=json.loads(j.read_text(encoding="utf-8")) if j.exists() else None)
        j = ML / f"npe_screens_{kind}.json"
        if j.exists():
            sj = json.loads(j.read_text(encoding="utf-8"))
            for meth, pre in (("npe_gauss2", "gauss"), ("npe_raw", "raw"), ("npe_gp", "gp")):
                cj = {lab: sj[f"{pre} | {lab}"] for lab in ("NISAR + GNSS", "NISAR only") if f"{pre} | {lab}" in sj}
                m[meth] = npe_block(ML / f"npe_screens_{kind}.npz", names[kind], None, 2000, prefix=pre + "__", cal_json=cj)
        j = ML / f"npe_flow_{kind}_gnss.json"
        if j.exists():      # competitive NPE: whitened + PCA flow, truncated second round
            fj = json.loads(j.read_text(encoding="utf-8"))
            blk = npe_block(ML / f"npe_flow_{kind}_gnss.npz", names[kind], None, 2000)
            if blk:
                for b in blk.values():
                    b["calibration_flow"] = {"68": [round(v * 100) for v in fj["calibration"]["0.68"]],
                                             "95": [round(v * 100) for v in fj["calibration"]["0.95"]]}
                    b["vs_mcmc"] = {nm: dict(shift=r.get("flow_round2_shift_sd"), width=r.get("flow_round2_width_ratio"))
                                    for nm, r in fj["params"].items()}
                m["npe_flow"] = blk
        out[kind] = {k: v for k, v in m.items() if v}
    return out


def flow_timings(tim):
    """Add the flow NPE to the timing table (training = both rounds incl. simulations)."""
    if tim is None:
        return None
    for kind in ("mogi", "ellipsoid"):
        j = ML / f"npe_flow_{kind}_gnss.json"
        if j.exists():
            fj = json.loads(j.read_text(encoding="utf-8"))
            T = fj["timing"]
            tim.setdefault(kind, {}).setdefault("NISAR + GNSS", {})["npe_flow"] = dict(
                train_s=round(T["sim1"] + T["sim2"] + T["train1"] + T["train2"]), n_sims=fj["n1"] + fj["n2"], infer_s=T["infer"])
    return tim


def png64(path):
    return "data:image/png;base64," + base64.b64encode(path.read_bytes()).decode() if path.exists() else None


def expensive_payload():
    """Layered-crust FE (MCMC on the emulator vs NPE on direct FE solves), real-time dike, Kilauea tephra figure."""
    out = {}
    if (ML / "fem_mcmc.json").exists():
        out["fem_mcmc"] = json.loads((ML / "fem_mcmc.json").read_text(encoding="utf-8"))
    if (ML / "fem_npe.json").exists():
        out["fem_npe"] = json.loads((ML / "fem_npe.json").read_text(encoding="utf-8"))
    mj = ML / "mcmc_mogi_screens.json"
    if mj.exists():
        z = np.load(ML / "mcmc_mogi_screens.npz")["NISAR_and_GNSS"]
        out["mogi_ref"] = {nm: np.round(np.percentile(z[:, i], [2.5, 16, 50, 84, 97.5]), 3).tolist()
                           for i, nm in enumerate(("east km", "north km", "depth km b.s.l.", "dV 1e6 m3/yr"))}
    if (ML / "dike_realtime.json").exists():
        d = json.loads((ML / "dike_realtime.json").read_text(encoding="utf-8"))
        d.pop("mcmc", None)
        out["dike"] = d
        out["dike_png"] = png64(ML / "dike_realtime.png")
    if (ML / "s1_grindavik_dike.json").exists():             # the REAL 10 Nov 2023 dike, Sentinel-1 via ASF HyP3
        s1 = json.loads((ML / "s1_grindavik_dike.json").read_text(encoding="utf-8"))
        mo = ML / "s1_grindavik_dike_mogi_r1.json"
        if mo.exists():
            m = json.loads(mo.read_text(encoding="utf-8"))
            s1["mogi_r1"] = dict(npe_volume=m["npe_volume"], length=m["npe"]["length km"], loglike=m["loglike"])
        out["s1"] = s1
        out["s1_png"] = png64(ML / "s1_grindavik_dike.png")
    out["tephra_png"] = png64(OUT.parent / "kilauea" / "tephra_evidence.png")
    out["pol_png"] = png64(OUT.parent / "kilauea" / "polarization.png")
    pj = OUT.parent / "kilauea" / "polarization.json"
    if pj.exists():
        pz = json.loads(pj.read_text(encoding="utf-8"))
        out["pol"] = dict(summary=pz["summary"], pairs=pz["pairs"], span_days=pz["span_days"],
                          sw_change=pz["per_pair"] and {c: pz["per_pair"][c]["mean_change_dB"]["SW"] for c in pz["channels"]},
                          noise={c: float(np.median(pz["per_pair"][c]["noise_sd_dB"])) for c in pz["channels"]})
    return out


def main():
    d = np.load(OUT / "breathing_data.npz")
    vec = json.loads((OUT / "breathing_vectors.json").read_text(encoding="utf-8"))
    lon, lat = d["lon"], d["lat"]
    ny, nx = lon.shape
    tr = Transformer.from_crs(4326, 3413, always_xy=True)
    X, Y = tr.transform(lon.astype(float), lat.astype(float))
    x0, dx = X[0, 0], X[0, 1] - X[0, 0]
    y0, dy = Y[0, 0], Y[1, 0] - Y[0, 0]
    GL, GA = np.meshgrid(np.linspace(BOX[0], BOX[2], NLON), np.linspace(BOX[3], BOX[1], NLAT))
    gx, gy = tr.transform(GL, GA)
    j = np.round((gx - x0) / dx).astype(int)
    i = np.round((gy - y0) / dy).astype(int)
    ok = (i >= 0) & (i < ny) & (j >= 0) & (j < nx)
    idx = np.where(ok, i * nx + j, -1).astype(np.int32)
    layers = json.loads((OUT / "map_layers.json").read_text(encoding="utf-8"))
    ml = ml_payload(d, ny, nx)
    payload = dict(
        ml=ml["meta"], gp_clean=ml["gp_clean"], gp_screen=ml["gp_screen"],
        cp_step=ml["cp_step"], cp_hinge=ml["cp_hinge"],
        box=BOX, nlon=NLON, nlat=NLAT, nx=nx, ny=ny, dates=[str(x) for x in d["dates"]],
        grid=dict(x0=float(x0), y0=float(y0), dx=float(dx), dy=float(dy)),
        idx=b64(idx), up=b64(mm16(d["up"])), east=b64(mm16(d["east"])), model=b64(mm16(d["model_up"])),
        rate=b64(mm16(d["rate"])), sigma=b64(mm16(d["sigma"])),
        lon=b64(lon.astype(np.float32)), lat=b64(lat.astype(np.float32)),
        hach=vec["hachures"], trust_counts=vec["trust_counts"], headline=vec["headline"], gnss=vec["gnss"],
        sections=vec["sections"], snr=vec["snr"], step=vec["step_cm"], model_file=vec["model_file"], layers=layers)
    tpl = (Path(__file__).with_name("earth_breathes_template.html")).read_text(encoding="utf-8")
    html = tpl.replace("/*__DATA__*/", "const DATA = " + json.dumps(payload) + ";")
    (OUT / "earth_breathes.html").write_text(html, encoding="utf-8")
    print("wrote", OUT / "earth_breathes.html", f"{len(html) / 1e6:.1f} MB; grid {ny}x{nx}, "
          f"{len(payload['dates'])} dates, {len(payload['hach'])} strokes")


if __name__ == "__main__":
    main()
