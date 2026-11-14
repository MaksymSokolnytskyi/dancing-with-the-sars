"""Distributed opening of the 10 Nov 2023 Grindavik dike (Sentinel-1 + GNSS), graben masked.

For a FIXED dike plane the surface displacement is linear in the opening of each patch, so the opening
distribution is a regularized linear least-squares problem that solves in a fraction of a second
(the classical approach, e.g. Sigmundsson et al. 2024).

  data      InSAR bins outside the graben mask + GNSS pseudo-bins (s1_common_v2.GeomV2); --no-gnss drops GNSS
  plane     vertical, centred on the uniform-dike NPE corridor, LENGTH x DEPTH km, PATCH x PATCH km patches
  unknowns  patch openings (>= 0), Svartsengi source volume change per window (layered FE or Mogi, --source),
            offset + plane per InSAR track
  smoothing Laplacian, weight from 5-fold cross-validation over 3-km blocks of InSAR bins (GNSS always kept)
  geometry  grid over strike / shift of the plane; the Svartsengi depth is PROFILED (H1): for every depth
            3-8 km the best solution, its chi2, dike volume and source dV, with and without GNSS
  errors    100 synthetic noise draws; analytic covariance of the active set -> corr(V_dike, dV_source)

Writes output/reykjanes/ml/s1_dike_distributed{_nognss}{_mogi}.json / .png and s1_tradeoff_profile*.json.
Usage: python scripts/reykjanes/s1_dike_distributed.py [--no-gnss] [--source fem|mogi]
"""
import itertools
import json
import sys
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from pyproj import Transformer
from scipy.optimize import lsq_linear

sys.argv += [] if "--source" in sys.argv else ["--source", "fem"]
import s1_grindavik_dike as S
from okada import tensile
from progress import Progress
from reykjanes_analysis import PLACES
from s1_common_v2 import GeomV2

ML = S.ML
LENGTH, DEPTH, PATCH = 20.0, 8.0, 1.0
GNSS = "--no-gnss" not in sys.argv
TAG = ("" if GNSS else "_nognss") + ("" if S.SOURCE == "fem" else "_" + S.SOURCE)
RNG = np.random.default_rng(5)


def plane(e0, n0, strike):
    u = np.array([np.sin(np.deg2rad(strike)), np.cos(np.deg2rad(strike))])
    ns, nd = int(LENGTH / PATCH), int(DEPTH / PATCH)
    return u, ns, nd, (np.arange(ns) + 0.5) * PATCH - LENGTH / 2, (np.arange(nd) + 0.5) * PATCH


def dike_cols(G, e0, n0, strike):
    u, ns, nd, cs, cd = plane(e0, n0, strike)
    cols = []
    for j in range(nd):
        for i in range(ns):
            d = tensile(G.e, G.n, e0 + cs[i] * u[0], n0 + cs[i] * u[1], cd[j], PATCH, PATCH, strike, 90.0, 1.0) * 100
            cols.append((d * G.vec).sum(1))
    return np.array(cols).T, (ns, nd)


def other_cols(G, src):
    cols = []
    for k in range(2):
        p = np.r_[np.zeros(7), src[0], src[1], src[2], 1.0 if k == 0 else 0.0, 1.0 if k == 1 else 0.0]
        cols.append((G.source_enu(p) * G.vec).sum(1))
    n = len(G.e)
    for k in range(2):
        m = np.zeros(n)
        m[:G.n_ins] = (G.trk[:G.n_ins] == k)
        cols += [m, m * G.e / 10, m * G.n / 10]
    return np.array(cols).T


def laplacian(ns, nd):
    L = []
    for j in range(nd):
        for i in range(ns):
            r = np.zeros(ns * nd)
            r[j * ns + i] = -4
            for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                ii, jj = i + di, j + dj
                if 0 <= ii < ns and 0 <= jj < nd:
                    r[jj * ns + ii] = 1
                elif jj < 0:
                    r[j * ns + i] += 1                      # free surface: no penalty for reaching the top
            L.append(r)
    return np.array(L)


def wmatrix(G):
    n = len(G.e)
    Wh = np.zeros((n, n))
    i0 = 0
    for k, W in enumerate(G.W):
        Wh[i0:i0 + len(W), i0:i0 + len(W)] = W                  # note: W includes the offset/plane marginalization;
        i0 += len(W)                                             # nuisance columns are also fitted (harmless, slightly conservative)
    if G.n_gnss:
        Wh[G.n_ins:, G.n_ins:] = np.diag(1 / G.gnss_sigma)
    return Wh


def solve(A, y, Wh, L, lam, npatch, rows=None):
    if rows is not None:
        A, y, Wh = A[rows], y[rows], Wh[np.ix_(rows, rows)]
    na = A.shape[1]
    Aw = np.vstack([Wh @ A, np.c_[lam * L, np.zeros((L.shape[0], na - npatch))]])
    yw = np.r_[Wh @ y, np.zeros(L.shape[0])]
    lb = np.r_[np.zeros(npatch), -np.inf * np.ones(na - npatch)]
    return lsq_linear(Aw, yw, bounds=(lb, np.inf * np.ones(na)), method="bvls").x


def chi2(A, x, y, Wh):
    r = Wh @ (y - A @ x)
    return float(r @ r)


def main():
    t_all = time.time()
    G = GeomV2(gnss=GNSS)
    Wh = wmatrix(G)
    y = G.y
    e_c, n_c, st_c = G.trace[0], G.trace[1], G.trace[2]
    z = np.load(ML / "s1_grindavik_dike.npz")
    src0 = np.median(z["ps"], 0)[7:10]
    # CV folds on InSAR blocks; GNSS always in the training set
    blk = np.floor(G.e[:G.n_ins] / 3).astype(int) * 1000 + np.floor(G.n[:G.n_ins] / 3).astype(int)
    ub = np.unique(blk)
    fold_of = dict(zip(ub, RNG.permutation(len(ub)) % 5))
    folds = np.r_[[fold_of[b] for b in blk], -np.ones(G.n_gnss, int)]
    D, (ns, nd) = dike_cols(G, e_c, n_c, st_c)
    A = np.c_[D, other_cols(G, src0)]
    L = laplacian(ns, nd)
    npatch = ns * nd
    lams = np.array([0.3, 0.6, 1, 2, 4, 8])
    cvs = []
    with Progress(f"розподілене розкриття{TAG}: вибір згладжування (CV)", total=len(lams), unit="λ") as pr:
        for lam in lams:
            err = 0.0
            for f in range(5):
                tr, te = folds != f, folds == f
                x = solve(A, y, Wh, L, lam, npatch, rows=tr)
                r = Wh[np.ix_(te, te)] @ (y[te] - A[te] @ x)
                err += float(r @ r)
            cvs.append(err)
            pr.update(1, note=f"λ={lam}: CV {err:.0f}")
    lam = float(lams[int(np.argmin(cvs))])
    print("CV:", dict(zip(lams.tolist(), np.round(cvs, 0).tolist())), "-> lambda", lam, flush=True)
    # geometry grid: strike and shift of the plane (source at the NPE position, depth profiled below)
    grid = list(itertools.product((-10, -5, 0, 5, 10), (-1.0, 0.0, 1.0), (-1.0, 0.0, 1.0)))
    res = []
    with Progress(f"розподілене розкриття{TAG}: сітка геометрії", total=len(grid), unit="варіант") as pr:
        for dst, de, dn in grid:
            Dg, _ = dike_cols(G, e_c + de, n_c + dn, st_c + dst)
            Ag = np.c_[Dg, other_cols(G, src0)]
            x = solve(Ag, y, Wh, L, lam, npatch)
            res.append((chi2(Ag, x, y, Wh), dst, de, dn))
            pr.update(1, note=f"найкраща χ² {min(q[0] for q in res):.0f}")
    res.sort()
    _, dst, de, dn = res[0]
    geo = dict(east=e_c + de, north=n_c + dn, strike=st_c + dst)
    D, _ = dike_cols(G, geo["east"], geo["north"], geo["strike"])
    # H1: profile over the Svartsengi source depth
    depths = np.arange(3.0, 8.01, 0.5)
    prof = []
    for zz in depths:
        Az = np.c_[D, other_cols(G, np.r_[src0[:2], zz])]
        x = solve(Az, y, Wh, L, lam, npatch)
        prof.append(dict(depth=float(zz), chi2=chi2(Az, x, y, Wh), V=float(x[:npatch].sum() * PATCH ** 2),
                         dV_A016=float(x[npatch]), dV_D155=float(x[npatch + 1]),
                         chi2_gnss=float(((y[G.n_ins:] - (Az @ x)[G.n_ins:]) / G.gnss_sigma) @ ((y[G.n_ins:] - (Az @ x)[G.n_ins:]) / G.gnss_sigma)) if G.n_gnss else None))
    best = min(prof, key=lambda p: p["chi2"])
    zbest = best["depth"]
    print("depth profile: " + "; ".join(f"z {p['depth']:.1f}: chi2 {p['chi2']:.0f}, V {p['V']:.0f}, dV {p['dV_A016']:.0f}/{p['dV_D155']:.0f}" for p in prof), flush=True)
    geo.update(src_east=float(src0[0]), src_north=float(src0[1]), src_depth=zbest)
    A = np.c_[D, other_cols(G, np.r_[src0[:2], zbest])]
    t0 = time.time()
    x = solve(A, y, Wh, L, lam, npatch)
    t_solve = time.time() - t0
    op = x[:npatch].reshape(nd, ns)
    vol = float(op.sum() * PATCH ** 2)
    # analytic covariance on the active set (open patches + free parameters)
    act = np.r_[np.flatnonzero(x[:npatch] > 0.02), np.arange(npatch, A.shape[1])]
    J = np.vstack([Wh @ A[:, act], np.c_[lam * L, np.zeros((L.shape[0], A.shape[1] - npatch))][:, act]])
    C = np.linalg.pinv(J.T @ J)
    g = np.zeros(len(act))
    g[:np.sum(act < npatch)] = PATCH ** 2                         # dV_dike / d(opening)
    iA = int(np.flatnonzero(act == npatch)[0])
    var_V, cov_VA = g @ C @ g, (C @ g)[iA]
    corr_V_dVA = float(cov_VA / np.sqrt(var_V * C[iA, iA]))
    # noise draws
    noise = G.sample_noise(100, RNG)
    vols, omax, dvs = [], [], []
    with Progress(f"розподілене розкриття{TAG}: 100 реалізацій шуму", total=100, unit="реаліз.") as pr:
        for k in range(100):
            xb = solve(A, A @ x + noise[k], Wh, L, lam, npatch)
            vols.append(xb[:npatch].sum() * PATCH ** 2); omax.append(xb[:npatch].max()); dvs.append(xb[npatch:npatch + 2])
            pr.update(1)
    vols, omax, dvs = np.array(vols), np.array(omax), np.array(dvs)
    pred = A @ x
    out = dict(gnss=GNSS, source=S.SOURCE, n_insar=G.n_ins, n_masked=G.n_masked, n_gnss=G.n_gnss, noise=G.noise_source,
               lambda_cv=lam, cv=dict(zip([str(v) for v in lams], [round(c, 1) for c in cvs])), geometry={k: round(float(v), 3) for k, v in geo.items()},
               volume_1e6m3=round(vol, 1), volume_noise_16_84=np.percentile(vols, [16, 84]).round(1).tolist(),
               volume_sd_analytic=round(float(np.sqrt(var_V)), 1), corr_Vdike_dVsrcA016=round(corr_V_dVA, 3),
               corr_Vdike_dVsrcA016_noise=round(float(np.corrcoef(vols, dvs[:, 0])[0, 1]), 3),
               max_opening_m=round(float(op.max()), 2), max_opening_noise_16_84=np.percentile(omax, [16, 84]).round(2).tolist(),
               length_open_gt_0p5m_km=float((op.max(0) > 0.5).sum() * PATCH),
               depth_of_max_opening_km=float((np.argmax(op.max(1)) + 0.5) * PATCH),
               src_dV_1e6m3=dict(A016=round(float(x[npatch]), 1), D155=round(float(x[npatch + 1]), 1)),
               src_dV_noise_16_84=dict(A016=np.percentile(dvs[:, 0], [16, 84]).round(1).tolist(), D155=np.percentile(dvs[:, 1], [16, 84]).round(1).tolist()),
               depth_profile=prof, chi2=round(chi2(A, x, y, Wh), 1), n_data=len(y),
               chi2_insar=round(float(((G.whiten(np.atleast_2d(y - pred)))[0, :G.n_ins] ** 2).sum()), 1),
               gnss_fit=[dict(label=l, obs=round(float(o), 2), pred=round(float(p), 2), sigma=round(float(s), 2))
                         for l, o, p, s in zip(G.gnss_label, y[G.n_ins:], pred[G.n_ins:], G.gnss_sigma)],
               seconds_one_solve=round(t_solve, 2), seconds_total=round(time.time() - t_all, 1), opening=op.round(3).tolist(), published=S.PUBLISHED)
    np.savez(ML / f"s1_dike_distributed{TAG}_resid.npz", e=G.e[:G.n_ins], n=G.n[:G.n_ins], trk=G.trk[:G.n_ins], resid=(y - pred)[:G.n_ins])
    (ML / f"s1_dike_distributed{TAG}.json").write_text(json.dumps(out, indent=1, ensure_ascii=False), encoding="utf-8")
    print(f"[{TAG or 'gnss+fem'}] V {vol:.0f} [{out['volume_noise_16_84']}] 1e6 m3 (analytic sd {np.sqrt(var_V):.0f}); max opening {op.max():.2f} m at "
          f"{out['depth_of_max_opening_km']} km; Svartsengi z {zbest} km, dV {x[npatch]:.0f} / {x[npatch + 1]:.0f}; corr(V, dV_A) {corr_V_dVA:+.2f}; "
          f"chi2 {out['chi2']:.0f} / {len(y)}; one solve {t_solve:.2f} s", flush=True)
    plot(G, A, x, op, geo, out)


def plot(G, A, x, op, geo, out):
    to_utm = Transformer.from_crs(4326, 32627, always_xy=True)
    sx, sy = to_utm.transform(*PLACES["Svartsengi"])
    pl = {k: ((to_utm.transform(*PLACES[k])[0] - sx) / 1000, (to_utm.transform(*PLACES[k])[1] - sy) / 1000)
          for k in ("Svartsengi", "Grindavik", "Sundhnukur")}
    ua = {"Svartsengi": "Свартсенгі", "Grindavik": "Ґріндавік", "Sundhnukur": "Сундгнукур"}
    pred = A @ x
    u, ns, nd, cs, cd = plane(geo["east"], geo["north"], geo["strike"])
    tr = (geo["east"] - u[0] * LENGTH / 2, geo["north"] - u[1] * LENGTH / 2, geo["east"] + u[0] * LENGTH / 2, geo["north"] + u[1] * LENGTH / 2)
    fig = plt.figure(figsize=(18, 11))
    gs = fig.add_gridspec(3, 4, height_ratios=[1, 1, .9], hspace=.42, wspace=.3)
    ttl_tr = {"A016": "висхідний 16, 31.10→12.11", "D155": "низхідний 155, 10.11→22.11"}
    i0 = 0
    for k, t in enumerate(S.TRACKS):
        n = G.sizes[k]
        sl = slice(i0, i0 + n)
        yy = G.y[sl]
        lim = np.percentile(np.abs(yy - np.median(yy)), 99)
        for j, (v, nm) in enumerate(((yy, "дані (поза маскою грабену)"), (pred[sl], "модель"), (yy - pred[sl], "залишок"))):
            a = fig.add_subplot(gs[k, j])
            off = np.median(yy) if j < 2 else 0
            sc = a.scatter(G.e[sl], G.n[sl], c=v, cmap="RdBu_r", vmin=off - lim, vmax=off + lim, s=12, marker="s", linewidths=0)
            a.plot([tr[0], tr[2]], [tr[1], tr[3]], "k-", lw=1)
            for nm2, (x0, y0) in pl.items():
                a.plot(x0, y0, "k^" if nm2 == "Svartsengi" else "ko", ms=4)
            a.set_aspect("equal"); a.set_xlim(S.BOX[0], S.BOX[1]); a.set_ylim(S.BOX[2], S.BOX[3])
            a.set_title(f"{ttl_tr[t]}: {nm}", fontsize=9)
            a.tick_params(labelsize=7)
            plt.colorbar(sc, ax=a, fraction=0.035, pad=0.02).set_label("см", fontsize=8)
        i0 += n
    # GNSS arrows (horizontal), both windows
    a = fig.add_subplot(gs[0:2, 3])
    a.plot([tr[0], tr[2]], [tr[1], tr[3]], "k-", lw=1)
    if G.n_gnss:
        gl = np.array(G.gnss_label)
        for w, col in (("A016", "tab:blue"), ("D155", "tab:orange")):
            for s in sorted({l.split()[0] for l in gl}):
                ie = np.flatnonzero(gl == f"{s} {w} E")
                inn = np.flatnonzero(gl == f"{s} {w} N")
                if len(ie) and len(inn):
                    i_e, i_n = G.n_ins + ie[0], G.n_ins + inn[0]
                    a.quiver(G.e[i_e], G.n[i_e], G.y[i_e], G.y[i_n], color=col, angles="xy", scale_units="xy", scale=8, width=.006)
                    a.quiver(G.e[i_e], G.n[i_e], pred[i_e], pred[i_n], color="k", angles="xy", scale_units="xy", scale=8, width=.003)
                    a.text(G.e[i_e] + .5, G.n[i_e] + .5, s, fontsize=7)
        a.quiver([], [], [], [], color="tab:blue", label="GNSS A016 (31.10→12.11)")
        a.plot([], [], color="tab:orange", label="GNSS D155 (10.11→22.11)")
        a.plot([], [], color="k", label="модель")
        a.legend(fontsize=7, loc="lower left")
    for nm2, (x0, y0) in pl.items():
        a.plot(x0, y0, "k^" if nm2 == "Svartsengi" else "ko", ms=4)
    a.set_aspect("equal"); a.set_xlim(-25, 30); a.set_ylim(-15, 32)
    a.set_title("GNSS: горизонтальні зміщення\n(1 км на мапі = 8 см)" if G.n_gnss else "без GNSS", fontsize=9)
    a.tick_params(labelsize=7)
    a = fig.add_subplot(gs[2, :3])
    im = a.imshow(op, extent=[-LENGTH / 2, LENGTH / 2, DEPTH, 0], cmap="magma_r", vmin=0, aspect="auto")
    plt.colorbar(im, ax=a, fraction=0.02, pad=0.01).set_label("розкриття, м", fontsize=8)
    for nm2, (x0, y0) in pl.items():
        along = (x0 - geo["east"]) * u[0] + (y0 - geo["north"]) * u[1]
        a.axvline(along, color="0.5", lw=.6, ls=":")
        a.text(along, -0.15, ua[nm2], fontsize=8, ha="center", va="bottom")
    a.set_xlabel("вздовж дайки від центру, км (ПдЗх → ПнСх)", fontsize=8); a.set_ylabel("глибина, км", fontsize=8)
    a.set_title(f"розподіл розкриття (патчі {PATCH:.0f} × {PATCH:.0f} км, простягання {geo['strike']:.0f}°)", fontsize=9, pad=14)
    a.tick_params(labelsize=7)
    a = fig.add_subplot(gs[2, 3])
    a.axis("off")
    rows = [["об’єм дайки, млн м³", f"{out['volume_1e6m3']:.0f} [{out['volume_noise_16_84'][0]:.0f}–{out['volume_noise_16_84'][1]:.0f}]", "130–139"],
            ["найбільше розкриття, м", f"{out['max_opening_m']:.1f} (на {out['depth_of_max_opening_km']:.1f} км)", "~8 (2–4 км)"],
            ["довжина з розкриттям > 0.5 м, км", f"{out['length_open_gt_0p5m_km']:.0f}", "~15"],
            ["Свартсенгі: глибина, км", f"{geo['src_depth']:.1f}", "4.0–4.2*"],
            ["Свартсенгі dV A016, млн м³", f"{out['src_dV_1e6m3']['A016']:.0f}", "−76…−82"],
            ["corr(V дайки, dV Св.)", f"{out['corr_Vdike_dVsrcA016']:+.2f}", "—"],
            ["χ² / дані", f"{out['chi2']:.0f} / {out['n_data']}", "—"]]
    tb = a.table(cellText=rows, colLabels=["", "ця модель", "опубл."], loc="center", cellLoc="center", colWidths=[.5, .32, .22])
    tb.auto_set_font_size(False); tb.set_fontsize(8); tb.scale(1, 1.45)
    a.text(0.5, -0.04, "[16–84 %]: 100 реалізацій шуму; опубл.: Sigmundsson et al. 2024;\n* однорідний півпростір (ми — шарувата кора)" if S.SOURCE == "fem"
           else "[16–84 %]: 100 реалізацій шуму; опубл.: Sigmundsson et al. 2024", fontsize=7.5, ha="center", va="top", transform=a.transAxes)
    fig.suptitle(f"Дайка 10.11.2023: розподілене розкриття, маска грабену ±1 км, " + ("InSAR + GNSS" if GNSS else "лише InSAR") + ", Свартсенгі — "
                 + ("шарувата кора (МСЕ)" if S.SOURCE == "fem" else "Моґі"), fontsize=12)
    fig.savefig(ML / f"s1_dike_distributed{TAG}.png", dpi=100, facecolor="white", bbox_inches="tight")
    print("wrote", ML / f"s1_dike_distributed{TAG}.png")


if __name__ == "__main__":
    main()
