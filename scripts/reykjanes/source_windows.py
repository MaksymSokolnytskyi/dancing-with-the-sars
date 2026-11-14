"""Does the Svartsengi source move in time? Joint NISAR + GNSS Mogi fits (per-track
offset + plane) in four ~3-month windows, with a bootstrap over the pairs of each window.

Usage: python scripts/reykjanes/source_windows.py [n_boot=25]
"""
import json
import sys
from datetime import datetime

import numpy as np
from pyproj import Transformer

import reykjanes_mogi as M
from reykjanes_analysis import MIN_PAIR_COH, OUT, PLACES, Geo, chain, los_cm, pair_coherence, pairs_by_track
from reykjanes_gnss import REF, STATIONS, rate, read_tenv3, relative
from reykjanes_joint import N_EFF, SIGMA_GNSS, SIGMA_INSAR, Joint

WINDOWS = [(datetime(2025, 11, 3), datetime(2026, 1, 31)), (datetime(2026, 2, 1), datetime(2026, 4, 30)),
           (datetime(2026, 5, 1), datetime(2026, 7, 31)), (datetime(2026, 8, 1), datetime(2026, 10, 17))]
RNG = np.random.default_rng(5)


def load_pairs(geo):
    out = {}
    for t in M.TRACKS:
        per, los = [], None
        for t0, t1, p in chain(pairs_by_track()[t]):
            if pair_coherence(p, geo) < MIN_PAIR_COH:
                continue
            d, g = los_cm(p, geo)
            if d is not None:
                per.append((t0 + (t1 - t0) / 2, d, (t1 - t0).days))
                los = g
        out[t] = dict(per=per, los=los)
    return out


def main(n_boot=25):
    M.RAMP = "plane"
    geo = Geo(pairs_by_track()["159A"][0][2])
    pairs = load_pairs(geo)
    to_ll = Transformer.from_crs(geo.epsg, 4326, always_xy=True)
    sx, sy = geo.tr.transform(*PLACES["Svartsengi"])
    rows = []
    for w0, w1 in WINDOWS:
        sel = {t: [(d, dt) for mid, d, dt in pairs[t]["per"] if w0 <= mid <= w1] for t in M.TRACKS}
        n = {t: len(v) for t, v in sel.items()}
        gn = {}
        for s in STATIONS + [REF]:
            if s == REF:
                gn[s] = dict(ll=read_tenv3(REF)[2], v=np.zeros(3))
            else:
                tt, x, ll = relative(s, REF)
                gn[s] = dict(ll=ll, v=rate(tt, x, w0, w1))
        data = {t: dict(per=sel[t], los=pairs[t]["los"], span=(w0, w1)) for t in M.TRACKS}
        P = Joint(geo, data, gn, "plane")
        start = np.array([sx / 1000, sy / 1000, 4.0])

        def fit(maps):
            obs = P.observations(maps)
            w_ins = {t: np.sqrt(N_EFF / obs[t]["y"].size) / SIGMA_INSAR for t in M.TRACKS}
            th, m = P.fit(obs, start, w_ins, 1 / SIGMA_GNSS)
            return th, m

        th, m = fit({t: M.rate_map(sel[t]) for t in M.TRACKS})
        boot = []
        for b in range(n_boot):
            maps = {t: M.rate_map(sel[t], RNG.integers(0, len(sel[t]), len(sel[t]))) for t in M.TRACKS}
            tb, mb = fit(maps)
            boot.append([tb[0] * 1000 - sx, tb[1] * 1000 - sy, tb[2], mb[0] * 1e3])
        boot = np.array(boot)
        lon, lat = to_ll.transform(th[0] * 1000, th[1] * 1000)
        q = lambda i: np.percentile(boot[:, i], [16, 84])
        row = dict(window=f"{w0:%Y-%m-%d}..{w1:%Y-%m-%d}", pairs=n, lon=round(lon, 4), lat=round(lat, 4),
                   dx_m=round(th[0] * 1000 - sx), dy_m=round(th[1] * 1000 - sy), depth=round(th[2], 2), dv=round(m[0] * 1e3, 1),
                   dx68=np.round(q(0)).tolist(), dy68=np.round(q(1)).tolist(), depth68=np.round(q(2), 2).tolist(),
                   dv68=np.round(q(3), 1).tolist(), gnss_seng_up=round(float(gn["SENG"]["v"][2]), 1))
        rows.append(row)
        print(f"{row['window']}: pairs {n}; centre {row['dx_m']:+5d} m E-axis {row['dy_m']:+5d} m N-axis of the Svartsengi marker "
              f"(68 %: {row['dx68']}, {row['dy68']}), depth {row['depth']:.2f} km {row['depth68']}, "
              f"dV {row['dv']:.1f} {row['dv68']} e6 m3/yr; GNSS SENG up {row['gnss_seng_up']:+.1f} cm/yr", flush=True)
    (OUT / "source_windows.json").write_text(json.dumps(rows, indent=1), encoding="utf-8")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 25)
