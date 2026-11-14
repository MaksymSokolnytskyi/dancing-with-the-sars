"""Collect the measured run times of the source-inference methods from their logs / json files
into output/reykjanes/ml/timings.json, for the MCMC vs NPE table of the app.

  MCMC  no training; every new data set costs a full sampler run (seconds, forward runs).
  NPE   training once (simulations + network training); a new data set costs one pass
        through the network (posterior samples in a fraction of a second).

Training times are read from the epoch lines of the logs (the Mogi logs print every 10
epochs, so the 25-epoch time is extrapolated from epoch 20). Inference was timed in the
first NPE runs (20 000 posterior samples); the screen NPEs use the same network.
Usage: python scripts/reykjanes/ml_timings.py
"""
import json
import re
from pathlib import Path

ML = Path(__file__).resolve().parents[2] / "output" / "reykjanes" / "ml"
EP = re.compile(r"\[(?:(\w+) \| )?(NISAR \+ GNSS|NISAR only)\] epoch (\d+): validation -log p [-\d.]+ +\((\d+) s\)")


def train_seconds(log, epochs=25):
    out = {}
    if not log.exists():
        return out
    for noise, data, ep, sec in EP.findall(log.read_text(encoding="utf-8")):
        out[(noise or "gauss", data)] = int(sec) * epochs / int(ep)       # last line wins
    return out


def first(pattern, text, cast=float):
    m = re.search(pattern, text)
    return cast(m.group(1)) if m else None


def main():
    res = {"mogi": {}, "ellipsoid": {}}
    inf = []
    for f in ("npe_log.txt", "npe_ellipsoid_log.txt"):
        p = ML / f
        if p.exists():
            inf += [float(v) for v in re.findall(r"posterior for the real data in ([\d.]+) s", p.read_text(encoding="utf-8"))]
    infer = max(inf) if inf else 0.1
    sims = {"mogi": 80000 * 9 / 40000, "ellipsoid": None}                 # vectorised Mogi: 9 s per 40 000
    el = (ML / "npe_screens_ellipsoid_log.txt")
    if el.exists():
        sims["ellipsoid"] = first(r"40000 elongated-body simulations in (\d+) s", el.read_text(encoding="utf-8"))
    n_sims = {"mogi": 80000, "ellipsoid": 40000}
    for kind in ("mogi", "ellipsoid"):
        first_log = ML / ("npe_log.txt" if kind == "mogi" else "npe_ellipsoid_log.txt")
        tr_first = train_seconds(first_log)
        tr_scr = train_seconds(ML / f"npe_screens_{kind}_log.txt")
        txt = first_log.read_text(encoding="utf-8") if first_log.exists() else ""
        m = re.search(r"simulated (\d+) data sets in (\d+) s", txt)
        first_sims = (int(m.group(1)), float(m.group(2))) if m else (40000, sims[kind] or 0)
        for data in ("NISAR + GNSS", "NISAR only"):
            row = {}
            for meth, src, ns, st in (("npe", tr_first.get(("gauss", data)), *first_sims),
                                      ("npe_raw", tr_scr.get(("raw", data)), n_sims[kind], sims[kind]),
                                      ("npe_gp", tr_scr.get(("gp", data)), n_sims[kind], sims[kind])):
                if src:
                    row[meth] = dict(train_s=round(src + (st or 0)), n_sims=ns, infer_s=infer)
            for meth, f in (("mcmc", f"mcmc_{kind}.json"), ("mcmc_screens", f"mcmc_{kind}_screens.json")):
                p = ML / f
                if p.exists():
                    j = json.loads(p.read_text(encoding="utf-8")).get(data)
                    if j:
                        row[meth] = dict(train_s=0, infer_s=round(j["seconds"]), n_forward=j["n_forward"])
            res[kind][data] = row
    res["note"] = ("one laptop CPU (8 threads); the NISAR-only ellipsoid MCMC runs shared the CPU with a second run, "
                   "so their times are about 2x the stand-alone cost")
    (ML / "timings.json").write_text(json.dumps(res, indent=1, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(res, indent=1, ensure_ascii=False))


if __name__ == "__main__":
    main()
