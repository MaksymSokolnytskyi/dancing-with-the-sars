"""C (pre-registered): injection-recovery detection threshold on the Mexico City asc rate map.

    uv run python scripts/c_threshold.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.ndimage import gaussian_filter

sys.path[:0] = [str(Path(__file__).parent), str(Path(__file__).parents[1] / "scripts")]
from trust_audit import SITES, plane_removed  # noqa: E402

from hachure import nisar

ROOT = Path(__file__).resolve().parents[1]
VS = [0.5, 1, 2, 4, 8]
SS = [0.5, 1, 2, 4]
MATCHED = "--matched" in sys.argv          # C-prime: least-squares bowl amplitude instead of the smoothed value


def smooth(a, sig_px):
    ok = np.isfinite(a)
    num = gaussian_filter(np.where(ok, a, 0), sig_px)
    den = gaussian_filter(ok.astype(float), sig_px)
    with np.errstate(all="ignore"):
        return num / den


def main():
    s = SITES["mexico-city"]
    bbox = s["bbox"]
    cache = ROOT / "data/raw/nisar_cache" / ("bbox_" + "_".join(f"{v:.3f}" for v in bbox) + f"_h{s['h']}")
    found = nisar.search((bbox[0] + bbox[2]) / 2, (bbox[1] + bbox[3]) / 2, direction="ASCENDING")
    chain = [g for g in nisar.sequential_chain(nisar.dominant_track(found)) if (cache / f"{g.name}.npz").exists()]
    pairs = []
    for g in chain:
        try:
            pairs.append(nisar.process(nisar.RawWindow.load(cache / f"{g.name}.npz"), ref_bbox=s["ref"], coh_min=0.3,
                                       iono="plane"))
        except ValueError as e:
            print("SKIP", e)
    rate, _, _, x, y, _ = nisar.stack_rate(pairs)
    rate *= 100                                                   # cm/yr, valid where finite
    dx = abs(x[1] - x[0])
    print(f"{len(pairs)} pairs, grid {rate.shape}, {dx:.0f} m")
    H, W = rate.shape
    YY, XX = np.mgrid[:H, :W] * dx
    base = plane_removed(rate)
    quiet = np.isfinite(base) & (np.abs(smooth(base, 2000 / dx / 2.355)) < 2)
    rng = np.random.default_rng(0)
    qidx = np.argwhere(quiet)
    res = {}

    def stat(field, cy, cx, sk):
        sm = smooth(field, sk * 1000 / dx / 2.355)
        out = []
        for i, (yy_, xx_) in enumerate(zip(cy, cx)):
            rr = np.hypot(YY[max(0, yy_ - R): yy_ + R, max(0, xx_ - R): xx_ + R] - yy_ * dx,
                          XX[max(0, yy_ - R): yy_ + R, max(0, xx_ - R): xx_ + R] - xx_ * dx)
            sub = field[max(0, yy_ - R): yy_ + R, max(0, xx_ - R): xx_ + R]
            ann = sub[(rr >= 2 * sk * 1000) & (rr <= 4 * sk * 1000)]
            if np.isfinite(ann).sum() < 10:
                out.append(np.nan)
            elif MATCHED:
                g = np.exp(-rr ** 2 / (2 * (sk * 1000 / 2.355) ** 2))
                m = (rr <= 2 * sk * 1000) & np.isfinite(sub)
                v = sub[m] - np.nanmedian(ann)
                out.append(float(np.dot(g[m], v) / np.dot(g[m], g[m])) if m.sum() >= 5 else np.nan)
            else:
                out.append(sm[yy_, xx_] - np.nanmedian(ann))
        return np.array(out)

    for sk in SS:
        R = int(4 * sk * 1000 / dx) + 2
        # noise points and injection centres, separated
        pick = qidx[rng.permutation(len(qidx))]
        centres = []
        for p in pick:
            if all(np.hypot(*(p - c)) * dx >= 3 * sk * 1000 for c in centres):
                centres.append(p)
            if len(centres) == 40:
                break
        centres = np.array(centres)
        far = qidx[rng.choice(len(qidx), 200, replace=False)]      # noise is measured on the un-injected map
        noise = stat(base, far[:, 0], far[:, 1], sk)
        sd = 1.4826 * np.nanmedian(np.abs(noise - np.nanmedian(noise)))
        fa = float(np.nanmean(noise < -3 * sd))
        sig = sk * 1000 / 2.355
        bowl = np.zeros_like(rate)
        for c in centres:
            bowl += -np.exp(-((YY - c[0] * dx) ** 2 + (XX - c[1] * dx) ** 2) / (2 * sig ** 2))
        row = {}
        for v in VS:
            inj = plane_removed(rate + v * bowl)
            dstat = stat(inj, centres[:, 0], centres[:, 1], sk)
            row[v] = float(np.nanmean(dstat < -3 * sd))
        res[sk] = dict(sigma_D=float(sd), false_alarm=fa, detect=row, n=len(centres))
        print(f"s={sk} km ({len(centres)} bowls): sigma_D {sd:.2f} cm/yr, false alarm {fa * 100:.1f}%, detection " +
              " ".join(f"{v}:{row[v] * 100:.0f}%" for v in VS), flush=True)

    def thr50(sk):
        d = res[sk]["detect"]
        lv, p = np.log(VS), np.array([d[v] for v in VS])
        for i in range(len(VS) - 1):
            if p[i] < 0.5 <= p[i + 1]:
                return float(np.exp(lv[i] + (0.5 - p[i]) / (p[i + 1] - p[i]) * (lv[i + 1] - lv[i])))
        return float(VS[0]) if p[0] >= 0.5 else float("inf")

    t = {sk: thr50(sk) for sk in SS}
    print("50% thresholds (cm/yr): " + " ".join(f"s={sk}: {t[sk]:.2f}" for sk in SS))
    if MATCHED:
        c1 = t[1] <= 3.5
        c2 = all(t[a] <= t[b] for a, b in zip(SS[:-1], SS[1:]))
        c3 = all(res[sk]["false_alarm"] <= 0.05 for sk in SS)
        print(f"C'1 {'PASS' if c1 else 'FAIL'}  C'2 {'PASS' if c2 else 'FAIL'}  C'3 {'PASS' if c3 else 'FAIL'}")
        name = "c_threshold_mexico_matched.json"
    else:
        c1 = all(res[sk]["detect"][v] >= 0.9 for sk in SS if sk >= 1 for v in VS if v >= 4)
        c2 = t[0.5] > t[2]
        c3 = all(res[sk]["false_alarm"] <= 0.05 for sk in SS)
        print(f"C1 {'PASS' if c1 else 'FAIL'}  C2 {'PASS' if c2 else 'FAIL'}  C3 {'PASS' if c3 else 'FAIL'}")
        name = "c_threshold_mexico.json"
    json.dump(dict(res=res, thr50=t, C1=c1, C2=c2, C3=c3), open(ROOT / "data/processed" / name, "w"), indent=1)


if __name__ == "__main__":
    main()
