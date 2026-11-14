"""GNSS displacements for the two Sentinel-1 windows of the 10 Nov 2023 dike (NGL daily series, IGS20).

Stations with data in autumn 2023 within ~40 km of Svartsengi: SENG (0.5 km), NYLA (18 km), VOGS (36 km),
REYK (37 km). (AUSV, GEVK, VOGC started later; the ~20 IMO stations used by Sigmundsson et al. 2024 are
not in the NGL archive.)

  windows   A016: 31 Oct 2023 18:59 -> 12 Nov 18:59   before = mean of 30-31 Oct, after = mean of 12-13 Nov
            D155: 10 Nov 2023 07:59 -> 22 Nov 07:59   before = mean of 08-09 Nov (the daily position of
                  10 Nov mixes pre- and post-dike), after = mean of 21-22 Nov
  secular   a linear trend per station and component fitted to 1 Jan - 20 Oct 2023 is removed
            (plate motion; the unrest started on 25 Oct)
  sigma     scatter of the detrended daily positions (x sqrt(1/n_before + 1/n_after)) combined with a model
            error of 0.3 cm horizontal / 0.8 cm vertical + 5 % of the displacement

Writes data/gnss/s1_gnss_2023.json (local frame: km from Svartsengi, UTM 27N; displacements in cm).
"""
import json
from datetime import datetime
from pathlib import Path

import numpy as np
from pyproj import Transformer

from reykjanes_analysis import PLACES
from reykjanes_gnss import read_tenv3

ROOT = Path(__file__).resolve().parents[2]
STATIONS = ("SENG", "NYLA", "VOGS", "REYK")
WINDOWS = {"A016": (("2023-10-30", "2023-10-31"), ("2023-11-12", "2023-11-13")),
           "D155": (("2023-11-08", "2023-11-09"), ("2023-11-21", "2023-11-22"))}
MODEL_ERR = np.array([0.3, 0.3, 0.8])          # cm (E, N, U)
REL_ERR = 0.05                                 # + 5 % of the signal: model error near the dike (patches, no graben faults)


def day(s):
    return datetime.fromisoformat(s)


def main():
    to_utm = Transformer.from_crs(4326, 32627, always_xy=True)
    sx, sy = to_utm.transform(*PLACES["Svartsengi"])
    out = {}
    for s in STATIONS:
        t, x, ll = read_tenv3(s)                                    # cm, ll = (lon, lat)
        dec = np.array([(a - datetime(2023, 1, 1)).days for a in t], float)
        fit = (t >= datetime(2023, 1, 1)) & (t <= datetime(2023, 10, 20))
        A = np.c_[np.ones(fit.sum()), dec[fit]]
        coef = np.linalg.lstsq(A, x[fit], rcond=None)[0]           # (2, 3)
        res = x[fit] - A @ coef
        scatter = 1.4826 * np.median(np.abs(res - np.median(res, 0)), 0)
        xd = x - np.c_[np.ones_like(dec), dec] @ coef               # detrended
        ux, uy = to_utm.transform(*ll)
        st = dict(e_km=round((ux - sx) / 1000, 3), n_km=round((uy - sy) / 1000, 3), lon=ll[0], lat=ll[1],
                  n_fit_days=int(fit.sum()), scatter_cm=scatter.round(3).tolist(), windows={})
        for w, (b, a) in WINDOWS.items():
            mb = np.array([day(b[0]) <= q <= day(b[1]) for q in t])
            ma = np.array([day(a[0]) <= q <= day(a[1]) for q in t])
            if mb.sum() == 0 or ma.sum() == 0:
                continue
            d = xd[ma].mean(0) - xd[mb].mean(0)
            sig = np.sqrt(scatter ** 2 * (1 / mb.sum() + 1 / ma.sum()) + MODEL_ERR ** 2 + (REL_ERR * d) ** 2)
            # equipment steps: a day-to-day jump > 10 cm at a station > 10 km from Svartsengi cannot be this event
            span = np.array([day(b[0]) <= q <= day(a[1]) for q in t])
            jumps = np.abs(np.diff(xd[span], axis=0)).max(0) if span.sum() > 1 else np.zeros(3)
            valid = ~((jumps > 10) & (np.hypot(st["e_km"], st["n_km"]) > 10))
            if not valid.all():
                print(f"   {s} {w}: equipment step in {[c for c, v in zip('ENU', valid) if not v]} (jump {jumps.round(1)} cm) -> dropped", flush=True)
            sig_base = np.sqrt(scatter ** 2 * (1 / mb.sum() + 1 / ma.sum()) + MODEL_ERR ** 2)
            st["windows"][w] = dict(enu_cm=d.round(2).tolist(), sigma_cm=sig.round(2).tolist(), sigma_base_cm=sig_base.round(3).tolist(), valid=valid.tolist(),
                                    n_before=int(mb.sum()), n_after=int(ma.sum()))
        out[s] = st
        print(f"{s}: ({st['e_km']:+.1f}, {st['n_km']:+.1f}) km, daily scatter {st['scatter_cm']} cm; " +
              "; ".join(f"{w} ENU {v['enu_cm']} ± {v['sigma_cm']}" for w, v in st["windows"].items()))
    (ROOT / "data" / "gnss" / "s1_gnss_2023.json").write_text(json.dumps(out, indent=1), encoding="utf-8")


if __name__ == "__main__":
    main()
