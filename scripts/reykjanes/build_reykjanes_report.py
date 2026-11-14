"""Build output/reykjanes/report.html: a self-contained report of the Reykjanes NISAR
analysis (figures embedded, data overlays and OpenStreetMap layers on a Leaflet map).

Needs the outputs of reykjanes_analysis / _gnss / _joint / _sources / _ellipsoid and
data/osm/reykjanes_layers.json (osm_extract.py).

Usage: python scripts/reykjanes/build_reykjanes_report.py
"""
import base64
import io
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colormaps
from scipy.interpolate import RegularGridInterpolator

import reykjanes_mogi as M
from reykjanes_analysis import PLACES, Geo, pairs_by_track
from reykjanes_ellipsoid import cdm, solve as cdm_solve
from reykjanes_gnss import REF, STATIONS
from reykjanes_sources import Data

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "output" / "reykjanes"
BOX = (-22.80, 63.76, -22.05, 64.02)
NLON, NLAT = 760, 520  # overlay raster on a lon/lat grid (~60 m)


def smooth_nan(img, sigma):
    from scipy.ndimage import gaussian_filter
    m = np.isfinite(img)
    num = gaussian_filter(np.where(m, img, 0.0), sigma)
    den = gaussian_filter(m.astype(float), sigma)
    return np.where(den > 0.4, num / np.maximum(den, 1e-6), np.nan)


def hachure_layers(field, lons, lats, step, start, spacing_m=220.0, max_len_m=1200.0, min_contour_m=2500.0):
    """Dufour-style hachures of a rate field (cm/yr) on a lon/lat grid, as GeoJSON.

    Strokes start on contours every `step` from +-`start` outwards and run downslope to the
    next contour; uplift levels are drawn red, subsidence levels blue, and strokes are grouped
    into five width classes by slope (cm/yr per km). Also returns the contours themselves."""
    from contourpy import contour_generator
    lat0, lon0 = float(np.mean(lats)), float(np.mean(lons))
    kx, ky = 111320.0 * np.cos(np.deg2rad(lat0)), 111200.0
    x = (lons - lon0) * kx                       # metres, increasing
    y = (lats[::-1] - lat0) * ky                 # metres, increasing (rows flipped to south-up order)
    z = field[::-1]
    gy, gx = np.gradient(z, y, x)                # cm/yr per metre
    fg = RegularGridInterpolator((y, x), np.dstack([gx, gy]), bounds_error=False, fill_value=np.nan)
    cg = contour_generator(x=x, y=y, z=np.ma.masked_invalid(z))
    zmin, zmax = np.nanmin(z), np.nanmax(z)
    levels = [l for l in np.arange(start, zmax, step)] + [-l for l in np.arange(start, -zmin, step)]
    strokes, contours = [], []
    for k, level in enumerate(sorted(levels)):
        for line in cg.lines(level):
            if len(line) < 3:
                continue
            d = np.r_[0, np.cumsum(np.hypot(*np.diff(line, axis=0).T))]
            if d[-1] < min_contour_m:
                continue
            contours.append((level, line))
            s = np.arange((k % 2) * spacing_m / 2, d[-1], spacing_m)
            pts = np.c_[np.interp(s, d, line[:, 0]), np.interp(s, d, line[:, 1])]
            g = fg(pts[:, ::-1])
            slope = np.hypot(g[:, 0], g[:, 1])
            ok = np.isfinite(slope) & (slope > step / max_len_m)
            length = step / slope[ok] * 0.85
            a = pts[ok]
            b = a - g[ok] / slope[ok, None] * length[:, None]   # downslope
            for p0, p1, sl in zip(a, b, slope[ok] * 1000):     # slope in cm/yr per km
                strokes.append((level, p0, p1, sl))
    to_ll = lambda p: [round(lon0 + p[0] / kx, 6), round(lat0 + p[1] / ky, 6)]
    if strokes:
        smax = np.percentile([s[3] for s in strokes], 97)
    feats = []
    for sign, color in ((1, "#b2182b"), (-1, "#2166ac")):
        for c in range(5):
            lo, hi = c / 5, (c + 1) / 5
            lines = [[to_ll(p0), to_ll(p1)] for lev, p0, p1, sl in strokes
                     if np.sign(lev) == sign and lo <= min(sl / smax, 0.999) < hi]
            if lines:
                feats.append(dict(type="Feature", geometry=dict(type="MultiLineString", coordinates=lines),
                                  properties=dict(color=color, w=round(0.6 + 2.6 * (c + 0.5) / 5, 2),
                                                  sign=sign, slope_lo=round(lo * smax, 1), slope_hi=round(hi * smax, 1))))
    cfeats = [dict(type="Feature", geometry=dict(type="LineString", coordinates=[to_ll(p) for p in line[::2]]),
                   properties=dict(level=round(float(level), 1))) for level, line in contours]
    return (dict(type="FeatureCollection", features=feats), dict(type="FeatureCollection", features=cfeats),
            dict(n=len(strokes), slope_max=float(smax) if strokes else 0.0))


def png_b64(fig_or_path):
    if isinstance(fig_or_path, (str, Path)):
        return base64.b64encode(Path(fig_or_path).read_bytes()).decode()
    buf = io.BytesIO()
    fig_or_path.savefig(buf, format="png", dpi=110, facecolor="white", bbox_inches="tight")
    plt.close(fig_or_path)
    return base64.b64encode(buf.getvalue()).decode()


def rgba_png(img, cmap, vmin, vmax, alpha=0.75):
    c = colormaps[cmap]((np.clip(img, vmin, vmax) - vmin) / (vmax - vmin))
    c[..., 3] = np.where(np.isfinite(img), alpha, 0)
    buf = io.BytesIO()
    plt.imsave(buf, c, format="png")
    return base64.b64encode(buf.getvalue()).decode()


def main():
    geo = Geo(pairs_by_track()["159A"][0][2])
    data = M.load(geo)
    maps = {t: M.rate_map(data[t]["per"]) for t in M.TRACKS}
    D = Data()
    fit = np.load(OUT / "ellipsoid_fit.npz")
    p_ell = tuple(fit["p"])
    m_ell, _ = cdm_solve(D, p_ell)
    U = m_ell[:3]
    off = {t: m_ell[3 + k] for k, t in enumerate(M.TRACKS)}

    # lon/lat grid and its position on the NISAR grid / local km frame
    lons = np.linspace(BOX[0], BOX[2], NLON)
    lats = np.linspace(BOX[3], BOX[1], NLAT)
    LO, LA = np.meshgrid(lons, lats)
    X, Y = geo.tr.transform(LO, LA)
    E, N = D.local(X / 1000, Y / 1000)

    def sample(img):
        f = RegularGridInterpolator((geo.y[::-1], geo.x), img[::-1], method="nearest", bounds_error=False,
                                    fill_value=np.nan)
        return f(np.stack([Y, X], -1))

    # Tie each track to GNSS: shift it by the mean of (NISAR - GNSS projected on the LOS) at SENG, AUSV,
    # GEVK, VOGC and NYLA (= 0, the GNSS reference). 159A keeps a broad residual ramp (+6..+9 cm/yr at
    # 8-14 km, the VOGC drift) that one offset cannot remove.
    from reykjanes_gnss import read_tenv3
    ties = {}
    for t in M.TRACKS:
        rows, res = [], []
        for s in STATIONS + [REF]:
            lon_s, lat_s = read_tenv3(s)[2]
            i, j = geo.ij(lon_s, lat_s)
            w = (slice(max(i - 2, 0), i + 3), slice(max(j - 2, 0), j + 3))
            nis = np.nanmedian(maps[t][w])
            l = data[t]["los"]
            vec_s = np.array([np.nanmedian(l["losUnitVectorX"][w]), np.nanmedian(l["losUnitVectorY"][w])])
            vec_s = np.r_[vec_s, np.sqrt(1 - np.sum(vec_s ** 2))]
            g_los = float(D.gn[s]["v"] @ vec_s)
            if np.isfinite(nis):
                rows.append([1.0, geo.x[j] / 1e4, geo.y[i] / 1e4])
                res.append(nis - g_los)
        # an offset only: the stations cluster around Svartsengi, so a fitted plane would be
        # extrapolated wildly across the rest of the frame
        shift = float(np.mean(res))
        after = np.array(res) - shift
        ties[t] = dict(before=np.round(res, 1).tolist(), after=np.round(after, 1).tolist(), shift=round(shift, 2))
        maps[t] = maps[t] - shift
        print(f"GNSS tie {t}: NISAR-GNSS at {STATIONS + [REF]} before {np.round(res, 1)}, after {np.round(after, 1)}")
    # LOS rate maps (tied to GNSS), vertical + east from the two geometries
    los = {t: sample(maps[t]) for t in M.TRACKS}
    vec = {t: {k: sample(data[t]["los"][k]) for k in ("losUnitVectorX", "incidenceAngle")} for t in M.TRACKS}
    ea, ed = vec["159A"]["losUnitVectorX"], vec["082D"]["losUnitVectorX"]
    ua, ud = np.cos(np.deg2rad(vec["159A"]["incidenceAngle"])), np.cos(np.deg2rad(vec["082D"]["incidenceAngle"]))
    det = ea * ud - ed * ua
    east = (los["159A"] * ud - los["082D"] * ua) / det
    up = (ea * los["082D"] - ed * los["159A"]) / det
    # ellipsoid model on the same grid (+ tilt)
    G = cdm(E.ravel(), N.ravel(), p_ell)
    mod = (G * U).sum(-1).reshape(E.shape + (3,))
    mup = mod[..., 2]
    gy, gx = np.gradient(mup)
    # grid spacing in km along rows / columns (lon/lat grid -> local km)
    dx_km = np.hypot(np.gradient(E, axis=1), np.gradient(N, axis=1))
    dy_km = np.hypot(np.gradient(E, axis=0), np.gradient(N, axis=0))
    tilt = np.hypot(gx / dx_km, gy / dy_km) * 10.0  # cm/yr per km -> microradian/yr
    resid = {t: los[t] - np.sum(np.stack([mod[..., 0], mod[..., 1], mod[..., 2]], -1) *
                                np.stack([vec[t]["losUnitVectorX"], sample(data[t]["los"]["losUnitVectorY"]),
                                          np.cos(np.deg2rad(vec[t]["incidenceAngle"]))], -1), -1)
             for t in M.TRACKS}

    # hachures: NISAR vertical smoothed (~0.5 km) and referenced to its own far field; the model as is
    up_s = smooth_nan(up, 8)  # tied to GNSS; hachures only where the signal stands above the noise
    up_s = np.where(np.hypot(E, N) < 12, up_s, np.nan)
    hach_obs, cont_obs, h_obs = hachure_layers(up_s, lons, lats, step=2.0, start=4.0, min_contour_m=3000.0)
    hach_mod, cont_mod, h_mod = hachure_layers(np.where(mup > 0.2, mup, np.nan), lons, lats, step=2.0, start=2.0)
    print("hachures: NISAR", h_obs, "model", h_mod)

    overlays = {
        "NISAR vertical (asc + desc), cm/yr": rgba_png(up, "RdBu_r", -20, 20),
        "NISAR east-west (east +), cm/yr": rgba_png(east, "PuOr_r", -10, 10),
        "NISAR LOS 159A ascending, cm/yr": rgba_png(los["159A"], "RdBu_r", -20, 20),
        "NISAR LOS 082D descending, cm/yr": rgba_png(los["082D"], "RdBu_r", -20, 20),
        "Model uplift (ellipsoid), cm/yr": rgba_png(np.where(mup > 0.5, mup, np.nan), "RdBu_r", -20, 20, 0.6),
        "Model tilt, microradian/yr": rgba_png(np.where(tilt > 2, tilt, np.nan), "magma_r", 0, 60, 0.65),
        "Residual 159A (data - model), cm/yr": rgba_png(resid["159A"], "RdBu_r", -10, 10),
    }

    # ---- infrastructure exposure from OSM
    osm = json.loads((ROOT / "data" / "osm" / "reykjanes_layers.json").read_text())
    f_up = RegularGridInterpolator((lats[::-1], lons), mup[::-1], bounds_error=False, fill_value=np.nan)
    f_tilt = RegularGridInterpolator((lats[::-1], lons), tilt[::-1], bounds_error=False, fill_value=np.nan)
    f_obs = RegularGridInterpolator((lats[::-1], lons), up[::-1], method="nearest", bounds_error=False, fill_value=np.nan)

    def line_stats(feats):
        out = []
        for ft in feats:
            c = np.array(ft["coords"])
            if len(c) < 2:
                continue
            pts = c[:, ::-1]
            u, tl = f_up(pts), f_tilt(pts)
            seg = np.hypot(np.diff(c[:, 0]) * 111.32 * np.cos(np.deg2rad(63.88)), np.diff(c[:, 1]) * 111.2)
            ft = dict(ft)
            ft.update(up_max=float(np.nanmax(u)) if np.isfinite(u).any() else None,
                      tilt_max=float(np.nanmax(tl)) if np.isfinite(tl).any() else None, km=float(seg.sum()))
            out.append(ft)
        return out

    infra = {k: line_stats(osm[k]) for k in ("pipelines", "power", "barriers")}
    roads_main = [r for r in osm["roads"] if r.get("highway") in ("motorway", "trunk", "primary", "secondary",
                                                                   "tertiary", "unclassified")]
    infra["roads"] = line_stats(roads_main)
    b = np.array([[x["lat"], x["lon"]] for x in osm["buildings"]])
    bu, bt, bo = f_up(b), f_tilt(b), f_obs(b)
    bins = [(-1e9, 2), (2, 5), (5, 10), (10, 20), (20, 1e9)]
    bstats = [(lo, hi, int(np.sum((bu >= lo) & (bu < hi)))) for lo, hi in bins]
    grind = (b[:, 0] > 63.83) & (b[:, 0] < 63.86) & (b[:, 1] > -22.46) & (b[:, 1] < -22.40)

    def summarize(name, feats, key=None):
        fs = [f for f in feats if f["up_max"] is not None]
        groups = {}
        for f in fs:
            g = f.get(key) or "other" if key else name
            groups.setdefault(g, []).append(f)
        rows = []
        for g, v in sorted(groups.items(), key=lambda kv: -sum(f["km"] for f in kv[1])):
            km = sum(f["km"] for f in v)
            km_hot = sum(f["km"] for f in v if f["up_max"] > 10)
            rows.append(dict(group=g, n=len(v), km=km, up=max(f["up_max"] for f in v),
                             tilt=max(f["tilt_max"] for f in v), km10=km_hot))
        return rows

    exposure = {
        "pipelines": summarize("pipelines", infra["pipelines"], "substance"),
        "power": summarize("power lines", infra["power"]),
        "barriers": summarize("barriers", infra["barriers"]),
        "roads": summarize("roads", infra["roads"], "highway"),
    }
    named_roads = {}
    for r in infra["roads"]:
        if r.get("name") and r["up_max"] is not None:
            nr = named_roads.setdefault(r["name"], dict(km=0, up=0, tilt=0))
            nr["km"] += r["km"]
            nr["up"] = max(nr["up"], r["up_max"])
            nr["tilt"] = max(nr["tilt"], r["tilt_max"])
    named_roads = sorted(named_roads.items(), key=lambda kv: -kv[1]["up"])[:8]

    # ---- GeoJSON for the map
    def gj_lines(feats, kind):
        return [dict(type="Feature", geometry=dict(type="LineString", coordinates=f["coords"]),
                     properties=dict(kind=kind, name=f.get("name", ""), sub=f.get("substance", f.get("highway", "")),
                                     op=f.get("operator", ""), up=round(f["up_max"], 1) if f["up_max"] is not None else None,
                                     tilt=round(f["tilt_max"], 0) if f["tilt_max"] is not None else None))
                for f in feats]

    layers = {k: dict(type="FeatureCollection", features=gj_lines(infra[k], k)) for k in infra}
    bpts = [[round(x["lat"], 5), round(x["lon"], 5), round(float(u), 1) if np.isfinite(u) else None]
            for x, u in zip(osm["buildings"], bu)]
    places = [dict(name=p["name"], lat=p["lat"], lon=p["lon"], place=p["place"]) for p in osm["places"] if p["name"]]

    # volcanic features and key sites, with the model value at each
    def at_point(lat, lon):
        u, tl = f_up([[lat, lon]])[0], f_tilt([[lat, lon]])[0]
        return (round(float(u), 1) if np.isfinite(u) else None, round(float(tl), 0) if np.isfinite(tl) else None)

    volc = []
    for v in osm["volcanoes"]:
        u, tl = at_point(v["lat"], v["lon"])
        volc.append(dict(name=v["name"] or "вулканічна форма", en=v["name_en"], lat=v["lat"], lon=v["lon"],
                         status=v["status"], ele=v["ele"], up=u, tilt=tl,
                         recent=v["status"] == "active" or any(y in v["name"] for y in ("2021", "2022", "2023", "2024", "2025"))))
    vents = [[v["lat"], v["lon"]] for v in osm["vents"]]
    # Sundhnukur crater row = the OSM vents between Stora-Skogfell and Hagafell, joined south -> north
    row = sorted([v for v in vents if 63.875 < v[0] < 63.90 and -22.40 < v[1] < -22.36], key=lambda p: p[0])
    peaks = [dict(name=p["name"], lat=p["lat"], lon=p["lon"], ele=p["ele"]) for p in osm["peaks"]]
    key = []
    for s in osm["sites"]:
        if s["kind"] == "hotel":
            continue
        name = s["name"] or ("Keflavík airport (KEF)" if s["kind"] == "airport" else "")
        if s["kind"] == "substation" and not s["name"] and abs(s["lat"] - 63.879) < 0.01 and abs(s["lon"] + 22.433) < 0.02:
            if any(k["kind"] == "power plant" and "Свартсенгі" in k["name"] for k in key):
                continue
            name = "Підстанція біля електростанції Свартсенгі (HS Orka)"  # the plant itself is not tagged in OSM
            s = dict(s, kind="power plant")
        if not name:
            continue
        u, tl = at_point(s["lat"], s["lon"])
        key.append(dict(name=name, kind=s["kind"], lat=s["lat"], lon=s["lon"], up=u, tilt=tl))
    for p in osm["places"]:
        if p["name"] in ("Grindavík", "Vogar", "Reykjanesbær", "Hafnir"):
            u, tl = at_point(p["lat"], p["lon"])
            key.append(dict(name=p["name"], kind="town", lat=p["lat"], lon=p["lon"], up=u, tilt=tl))
    main_roads = [r for r in infra["roads"] if r.get("name") in ("Grindavíkurvegur", "Nesvegur", "Bláalónsvegur",
                                                                  "Norðurljósavegur", "Reykjanesbraut")]
    print("key sites:", [(k["name"], k["kind"], k["up"], k["tilt"]) for k in key])
    print("crater row vents:", len(row), "volcanoes:", len(volc), "recent:", [v["name"] for v in volc if v["recent"]])

    # sources
    e0, n0, top, a, bb, c, strike = p_ell
    t = np.linspace(0, 2 * np.pi, 72)
    s_, n_ = np.array([np.sin(np.deg2rad(strike)), np.cos(np.deg2rad(strike))]), \
        np.array([-np.cos(np.deg2rad(strike)), np.sin(np.deg2rad(strike))])
    ring = np.array([e0, n0])[:, None] + np.outer(s_, a * np.cos(t)) + np.outer(n_, max(bb, 0.4) * np.sin(t))
    ring_ll = [list(D.lonlat(x, y))[::-1] for x, y in ring.T]
    joint = np.load(OUT / "joint_fit_offset.npz")
    sources = dict(ellipsoid=ring_ll, mogi=[float(joint["lat"]), float(joint["lon"])],
                   ell_centre=list(D.lonlat(e0, n0))[::-1])

    # GNSS with model predictions
    ref = D.gn[REF]["en"]
    Gr = (cdm(np.array([ref[0]]), np.array([ref[1]]), p_ell)[0] * U).sum(-1)
    gnss = []
    from reykjanes_gnss import read_tenv3
    for s in STATIONS + [REF]:
        ll = read_tenv3(s)[2]
        en = D.gn[s]["en"]
        pred = (cdm(np.array([en[0]]), np.array([en[1]]), p_ell)[0] * U).sum(-1) - Gr
        v = D.gn[s]["v"]
        gnss.append(dict(name=s, lat=ll[1], lon=ll[0], v=[round(float(x), 1) for x in v],
                         m=[round(float(x), 1) for x in pred], ref=s == REF))

    figs = {k: png_b64(OUT / f) for k, f in {
        "timeseries": "timeseries_local.png", "gnss": "gnss_vs_nisar.png", "fringes": "fringes_159A_20251201_20260320.png",
        "pairs": "pairs_159A.png", "mogi": "mogi_fit_offset.png", "joint": "joint_fit_offset.png",
        "sources": "source_models.png", "ellipsoid": "ellipsoid_fit.png", "mogi_gnss": "mogi_vs_gnss.png"}.items()
        if (OUT / f).exists()}

    stats = dict(n_buildings=int(len(b)), bstats=bstats, grind_n=int(grind.sum()),
                 grind_up=[float(np.nanmin(bu[grind])), float(np.nanmax(bu[grind]))],
                 grind_obs=[float(np.nanpercentile(bo[grind], 10)), float(np.nanpercentile(bo[grind], 90))],
                 grind_tilt=float(np.nanmax(bt[grind])),
                 exposure=exposure, named_roads=named_roads,
                 ell=dict(depth=top + c, length=2 * a, strike=strike % 180, vol=float(4 * (U[0] * bb * c + U[1] * a * c + U[2] * a * bb) * 1e6)),
                 max_tilt=float(np.nanmax(tilt)), max_up=float(np.nanmax(mup)))
    layers["keyroads"] = dict(type="FeatureCollection", features=gj_lines(main_roads, "keyroads"))
    layers["hach_obs"], layers["cont_obs"] = hach_obs, cont_obs
    layers["hach_mod"], layers["cont_mod"] = hach_mod, cont_mod
    stats["hach"] = dict(obs=h_obs, mod=h_mod)
    stats["ties"] = ties
    volcanic = dict(volcanoes=volc, vents=vents, row=row, peaks=peaks, key=key)
    # vector layers for other pages (the "Earth breathes" web map)
    (OUT / "map_layers.json").write_text(json.dumps(dict(box=BOX, layers=layers, buildings=bpts, places=places,
                                                         sources=sources, gnss=gnss, volcanic=volcanic)), encoding="utf-8")
    html = render(overlays, layers, bpts, places, sources, gnss, figs, stats, volcanic)
    (OUT / "report.html").write_text(html, encoding="utf-8")
    print("wrote", OUT / "report.html", f"{len(html) / 1e6:.1f} MB")
    print(json.dumps({k: v for k, v in stats.items() if k not in ("exposure",)}, default=str)[:1500])
    for k, rows in exposure.items():
        print(k, [(r["group"], round(r["km"], 1), round(r["up"], 1), round(r["tilt"]), round(r["km10"], 1)) for r in rows])


def render(overlays, layers, bpts, places, sources, gnss, figs, stats, volcanic):
    tpl = (Path(__file__).with_name("reykjanes_report_template.html")).read_text(encoding="utf-8")
    payload = dict(box=BOX, overlays=overlays, layers=layers, buildings=bpts, places=places, sources=sources,
                   gnss=gnss, stats=stats, volcanic=volcanic)
    html = tpl.replace("/*__DATA__*/", "const DATA = " + json.dumps(payload) + ";")
    for k, v in figs.items():
        html = html.replace(f"__FIG_{k}__", f"data:image/png;base64,{v}")
    return html


if __name__ == "__main__":
    main()
