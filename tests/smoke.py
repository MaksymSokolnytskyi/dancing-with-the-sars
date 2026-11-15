"""Offline smoke test: no network, no token, no downloads. Run before every push.

    uv run python tests/smoke.py

Checks:
  1. every engine module imports;
  2. every plate config in sites/ is well formed, and its lock file (if any) parses;
  3. every built plate in web/data has a manifest whose layer files exist, and scenes.json lists only built plates;
  4. the NISAR reader works end to end on a tiny synthetic GUNW-like HDF5 (read_raw -> process);
  5. the hachure tracer draws strokes on a synthetic hill.
"""
from __future__ import annotations

import importlib
import importlib.util
import json
import sys
import tempfile
import traceback
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

failures: list[str] = []


def check(name):
    def wrap(fn):
        try:
            fn()
            print(f"  ok  {name}")
        except Exception:                                        # noqa: BLE001
            failures.append(name)
            print(f"  !!  {name}\n" + "".join("        " + l for l in traceback.format_exc().splitlines(True)[-6:]))
        return fn
    return wrap


@check("engine modules import")
def _imports():
    for m in ("cli", "dem", "engine", "export", "field", "gcov", "gridding", "nisar", "synthetic", "trace", "trust"):
        importlib.import_module(f"hachure.{m}")
    # okada is PyTorch only (used by scripts/sbi_*); torch comes with `uv sync --group ml`, so check it only when present
    if importlib.util.find_spec("torch") is not None:
        importlib.import_module("hachure.okada")
    else:
        print("  (torch not installed - hachure.okada not checked; it needs `uv sync --group ml`)")


@check("plate configs and lock files")
def _sites():
    from hachure import engine, nisar
    for p in sorted((ROOT / "sites").glob("*.json")):
        if p.name.endswith((".story.json", ".lock.json")):
            continue
        cfg = engine.load_site(p)
        assert cfg["id"] == p.stem, f"{p.name}: id {cfg['id']!r} differs from the file name"
        w, s, e, n = cfg["bbox"]
        assert -180 <= w < e <= 180 and -90 <= s < n <= 90, f"{p.name}: bad bbox"
        story = p.with_name(f"{cfg['id']}.story.json")
        if story.exists():
            json.loads(story.read_text(encoding="utf-8"))
        lock = p.with_name(f"{cfg['id']}.lock.json")
        if lock.exists():
            lk = json.loads(lock.read_text(encoding="utf-8"))
            assert lk["granules"], f"{lock.name}: empty"
            for g in lk["granules"]:
                assert g["url"].startswith("https://") and g["url"].endswith(".h5"), f"{lock.name}: bad url"
                nisar._parse_dates(g["name"])


@check("built plates in web/data")
def _web():
    idx = json.loads((ROOT / "web/data/scenes.json").read_text(encoding="utf-8"))
    scenes, ids = idx["scenes"], {s["id"] for s in idx["scenes"]}
    for ch in idx.get("chapters", []):                            # the atlas contents name only listed plates
        assert set(ch["scenes"]) <= ids, f"chapter {ch['id']}: unknown {set(ch['scenes']) - ids}"
    gl = json.loads((ROOT / "web/data/glossary.json").read_text(encoding="utf-8"))["terms"]
    assert len({g["id"] for g in gl}) == len(gl) and all(g["title"].get("uk") and g["text"].get("uk") for g in gl), "glossary"
    for sc in scenes:
        if sc.get("href"):                                       # a packed page of its own (Reykjanes)
            assert (ROOT / "web" / sc["href"]).exists(), f"{sc['id']}: missing {sc['href']}"
            continue
        d = ROOT / "web/data" / sc["id"]
        man = json.loads((d / "manifest.json").read_text(encoding="utf-8"))
        for layer in man["layers"]:
            if layer.get("file"):
                assert (d / layer["file"]).exists(), f"{sc['id']}: missing {layer['file']}"


@check("NISAR reader on a synthetic GUNW-like file")
def _reader():
    import h5py
    from hachure import nisar
    epsg, h, w = 32614, 60, 80                                   # UTM 14N, Mexico City
    xs = 480_000 + 80.0 * np.arange(w)
    ys = 2_150_000 - 80.0 * np.arange(h)
    rate = -0.3 * np.exp(-(((np.arange(w)[None, :] - 40) / 15) ** 2 + ((np.arange(h)[:, None] - 30) / 12) ** 2))
    lam = nisar.C / 1.2575e9
    phase = (rate * 12 / 365.25 * np.cos(np.radians(37))) / (nisar.LOS_SIGN * lam / (4 * np.pi))
    with tempfile.TemporaryDirectory() as tmp:
        f = Path(tmp) / "NISAR_L2_PR_GUNW_001_001_A_001_002_4000_SH_20260101T000000_20260101T000000_20260113T000000_20260113T000000_P00000_N_F_J_001.h5"
        with h5py.File(f, "w") as hf:
            g = hf.create_group("science/LSAR/GUNW/grids/frequencyA/unwrappedInterferogram")
            g["xCoordinates"], g["yCoordinates"] = xs, ys
            g.create_dataset("projection", data=np.uint32(epsg)).attrs["epsg_code"] = epsg
            hf["science/LSAR/GUNW/grids/frequencyA/centerFrequency"] = 1.2575e9
            p = g.create_group("HH")
            p["unwrappedPhase"] = phase.astype("float32")
            p["coherenceMagnitude"] = np.full((h, w), 0.8, "float32")
            p["connectedComponents"] = np.ones((h, w), np.uint16)
        from rasterio.warp import transform as warp
        lo, la = warp(f"EPSG:{epsg}", "EPSG:4326", [xs[2], xs[-3]], [ys[-3], ys[2]])
        bbox = (lo[0], la[0], lo[1], la[1])                      # just inside the synthetic grid
        with nisar.open_h5(str(f)) as hf:
            raw = nisar.read_raw(hf, bbox, 12.0)
        assert raw.phase.size > 0 and np.isfinite(raw.phase).any()
        pair = nisar.process(raw, iono=False)
        d = pair.disp_v
        bowl = float(np.nanmin(d) - np.nanmedian(d[:, :3]))      # centre minus the far edge
        want = -0.3 * 12 / 365.25                                # 30 cm/yr over 12 days, vertical
        assert abs(bowl - want) < 0.2 * abs(want), f"bowl {bowl * 100:.2f} cm, expected {want * 100:.2f} cm"


@check("hachures on a synthetic hill")
def _hachures():
    from hachure import field, trace
    yy, xx = np.mgrid[0:120, 0:120]
    z = 400 * np.exp(-((xx - 60) ** 2 + (yy - 60) ** 2) / (2 * 25 ** 2))
    strokes = trace.hachures(field.prepare(z, 30.0), interval_m=40.0)
    assert len(strokes) > 50, f"only {len(strokes)} strokes"


print("all good" if not failures else f"{len(failures)} failed: {', '.join(failures)}")
sys.exit(1 if failures else 0)
