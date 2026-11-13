"""One command for the whole pipeline, so every teammate fetches and builds the same way.

    uv run sars doctor                              # python, packages, Earthdata token (never printed), network
    uv run sars sites                               # the plates in sites/ and what is cached for each
    uv run sars search sites/mexico-city-asc.json   # NISAR GUNW pairs over the plate (no login needed)
    uv run sars fetch  sites/mexico-city-asc.json   # everything the plate needs; resumable; window reads only
    uv run sars fetch  sites/mexico-city-asc.json --dry-run      # what would be read, nothing downloaded
    uv run sars build  sites/mexico-city-asc.json   # layers for the web map (web/data/<id>/)
    uv run sars story  sites/mexico-city-asc.json   # republish the plate's pages only (seconds)
    uv run sars pack   mexico-city-asc              # one-file share pack in share/
    uv run sars serve                               # http://localhost:8765/?scene=<id>&theme=scroll
    uv run sars new my-place --bbox W S E N --direction ASCENDING   # a new plate from the template

NISAR pairs are pinned in sites/<id>.lock.json (granule names + URLs). `fetch` reads the lock if it exists, so the whole
team works on the same pairs; `fetch --update` searches ASF again and rewrites the lock. Only the plate's bbox window of
each GUNW is read over HTTP (a few MB per pair instead of ~2 GB per file), and cached in data/raw/nisar_cache/.

The Earthdata token stays outside the repo: env EARTHDATA_TOKEN or a one-line file ~/.edl_token.
"""
from __future__ import annotations

import argparse
import json
import os
import runpy
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def _run_script(name: str, *args: str):
    """Run scripts/<name>.py as if from the command line (keeps those scripts usable on their own)."""
    path = ROOT / "scripts" / f"{name}.py"
    old = sys.argv
    sys.argv = [str(path), *map(str, args)]
    sys.path.insert(0, str(path.parent))
    try:
        runpy.run_path(str(path), run_name="__main__")
    finally:
        sys.argv = old
        sys.path.remove(str(path.parent))


def _site(path: str) -> dict:
    from hachure import engine
    p = Path(path)
    if not p.exists() and (ROOT / "sites" / f"{path}.json").exists():
        p = ROOT / "sites" / f"{path}.json"                     # `sars fetch mexico-city-asc` works too
    cfg = engine.load_site(p)
    cfg["_path"] = p
    return cfg


def _lock_path(cfg) -> Path:
    return Path(cfg["_path"]).with_name(f"{cfg['id']}.lock.json")


def _mb(n: float) -> str:
    return f"{n / 1e6:,.1f} MB".replace(",", " ")


# ---------------------------------------------------------------- doctor
# Windows keeps paths under 260 characters unless long paths are switched on. NISAR file names are ~125 characters,
# so a repo cloned deep (e.g. inside OneDrive or a temp folder) breaks the cache. Check before fetching.
GUNW_NAME_LEN = 130


def _long_paths_on() -> bool:
    if os.name != "nt":
        return True
    try:
        import winreg
        k = winreg.OpenKey(winreg.HKEY_LOCAL_MACHINE, r"SYSTEM\CurrentControlSet\Control\FileSystem")
        return bool(winreg.QueryValueEx(k, "LongPathsEnabled")[0])
    except OSError:
        return False


def _path_problem(folder: Path) -> str | None:
    n = len(str(folder.resolve())) + 1 + GUNW_NAME_LEN
    if n <= 259 or _long_paths_on():
        return None
    over = n - 259
    return (f"cache paths would be ~{n} characters, Windows allows 259. Clone the repo into a shorter folder "
            f"(at least {over} characters shorter, e.g. C:/sars), or ask an admin to enable Windows long paths "
            f"(LongPathsEnabled).")


def cmd_doctor(a):
    ok = True
    print(f"python {sys.version.split()[0]} at {sys.executable}")
    for mod in ("numpy", "scipy", "rasterio", "h5py", "fsspec", "contourpy", "requests"):
        try:
            __import__(mod)
            print(f"  ok  {mod}")
        except Exception as e:                                   # noqa: BLE001
            ok = False
            print(f"  !!  {mod}: {e}")
    from hachure import nisar
    tok = nisar._token()
    src = "env EARTHDATA_TOKEN" if os.environ.get("EARTHDATA_TOKEN") else nisar.TOKEN_FILE
    if not tok:
        ok = False
        print(f"  !!  no Earthdata token. Make one at https://urs.earthdata.nasa.gov (Generate Token) and save it as one "
              f"line in {nisar.TOKEN_FILE}, or set EARTHDATA_TOKEN. Never paste it into chats or docs.")
    else:
        print(f"  ok  Earthdata token found ({src}, {len(tok)} characters; not shown)")
    try:
        t = time.time()
        found = nisar.search(-99.13, 19.43, start="2026-01-01", end="2026-01-31")
        print(f"  ok  ASF search answers ({len(found)} GUNW over Mexico City in Jan 2026, {time.time() - t:.1f}s)")
        if tok and found and a.deep:
            t = time.time()
            nisar._signed_url(found[0].url, tok)
            print(f"  ok  token accepted by ASF ({time.time() - t:.1f}s)")
    except Exception as e:                                       # noqa: BLE001
        ok = False
        print(f"  !!  ASF: {e}")
    worst = ROOT / "data/raw/nisar_cache" / ("bbox_" + "_".join(["-120.000"] * 4) + "_h5000")
    if (msg := _path_problem(worst)):
        ok = False
        print(f"  !!  {msg}")
    else:
        print("  ok  path length (NISAR cache file names fit)")
    print("all good" if ok else "fix the !! lines above")
    return 0 if ok else 1


# ---------------------------------------------------------------- sites
def cmd_sites(a):
    from hachure import engine
    for p in sorted((ROOT / "sites").glob("*.json")):
        if p.name.endswith((".story.json", ".lock.json")):
            continue
        cfg = engine.load_site(p)
        cache = engine.cache_dir(cfg) if "cache_height" in cfg else None
        n = len(list(cache.glob("*.npz"))) if cache and cache.exists() else 0
        lock = p.with_name(f"{cfg['id']}.lock.json")
        nl = len(json.loads(lock.read_text(encoding="utf-8"))["granules"]) if lock.exists() else None
        built = (ROOT / "web/data" / cfg["id"] / "manifest.json").exists()
        title = cfg.get("title")
        title = title.get("en") if isinstance(title, dict) else title
        print(f"{cfg['id']:22s} {str(title or '')[:34]:34s} cached {n:3d} pairs"
              f"{f', lock {nl}' if nl is not None else ', no lock'}{', built' if built else ''}")


# ---------------------------------------------------------------- search / fetch
def _search(cfg, a):
    from hachure import nisar
    w, s, e, n = cfg["bbox"]
    found = nisar.search((w + e) / 2, (s + n) / 2, direction=cfg.get("direction"),
                         start=a.start or cfg.get("start"), end=a.end or cfg.get("end"))
    track = nisar.dominant_track(found)
    chain = nisar.sequential_chain(track) if not a.all_pairs else sorted(track, key=lambda g: g.ref_date)
    return found, track, chain


def cmd_search(a):
    cfg = _site(a.site)
    found, track, chain = _search(cfg, a)
    tracks = {}
    for g in found:
        tracks.setdefault((g.direction, g.path, g.frame), []).append(g)
    print(f"{cfg['id']}: {len(found)} GUNW pairs over the bbox centre")
    for (d, p, f), gs in sorted(tracks.items(), key=lambda kv: -len(kv[1])):
        print(f"  {d[:4]} path {p} frame {f}: {len(gs)} pairs, {gs[0].ref_date:%Y-%m-%d} -> {gs[-1].sec_date:%Y-%m-%d}")
    if chain:
        print(f"chain used by fetch: {len(chain)} pairs, {chain[0].ref_date:%Y-%m-%d} -> {chain[-1].sec_date:%Y-%m-%d}")
        if a.verbose:
            for g in chain:
                print(f"    {g.ref_date:%Y-%m-%d} -> {g.sec_date:%Y-%m-%d} ({g.days:.0f} d)  {g.name}")


def _write_lock(cfg, chain, path: Path):
    d = chain[0] if chain else None
    lock = {"site": cfg["id"], "bbox": cfg["bbox"], "direction": cfg.get("direction"),
            "track": {"path": d.path, "frame": d.frame} if d else None,
            "searched": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "granules": [{"name": g.name, "url": g.url, "ref": g.ref_date.isoformat(), "sec": g.sec_date.isoformat(),
                          "days": round(g.days, 3)} for g in chain]}
    path.write_text(json.dumps(lock, indent=1), encoding="utf-8")


def cmd_lock(a):
    """Pin the pairs this machine already has in its cache (what the built plate used) as the plate's lock."""
    from hachure import engine
    cfg = _site(a.site)
    used = [g.name for g in engine.cached_chain(cfg)]             # exactly the pairs `build` reads
    a.start = a.end = None
    a.all_pairs = True
    found, _, _ = _search(cfg, a)
    url = {g.name: g for g in found}
    chain = [url[n] for n in used if n in url]
    missing = len(used) - len(chain)
    if not chain:
        print("nothing in the cache matches an ASF granule of this plate's direction")
        return 1
    chain.sort(key=lambda g: g.ref_date)
    _write_lock(cfg, chain, _lock_path(cfg))
    print(f"{_lock_path(cfg).relative_to(ROOT)}: {len(chain)} pairs, {chain[0].ref_date:%Y-%m-%d} -> "
          f"{chain[-1].sec_date:%Y-%m-%d}" + (f"; {missing} cached names not found on ASF (renamed or withdrawn)" if missing > 0 else ""))
    return 0


def _locked(path: Path):
    from hachure import nisar
    lock = json.loads(path.read_text(encoding="utf-8"))
    out = []
    for x in lock["granules"]:
        ref, sec = datetime.fromisoformat(x["ref"]), datetime.fromisoformat(x["sec"])
        out.append(nisar.Granule(x["name"], x["url"], lock.get("direction") or "", lock["track"]["path"],
                                 lock["track"]["frame"], ref, sec, None))
    return out


def fetch_nisar(cfg, a) -> int:
    from hachure import engine, nisar
    if "cache_height" not in cfg:
        print("  nisar: this plate reads no GUNW windows (no cache_height), skipped")
        return 0
    lock = _lock_path(cfg)
    if lock.exists() and not a.update:
        chain = _locked(lock)
        how = f"lock {lock.name}"
    else:
        _, _, chain = _search(cfg, a)
        if not chain:
            print("  nisar: ASF has no pairs here (yet)")
            return 1
        if not a.dry_run:
            _write_lock(cfg, chain, lock)
        how = "ASF search" + (" (lock written)" if not a.dry_run else "")
    if a.max_pairs:
        chain = chain[-a.max_pairs:]
    cache = engine.cache_dir(cfg)
    todo = [g for g in chain if not (cache / f"{g.name}.npz").exists()]
    have = [p.stat().st_size for p in cache.glob("*.npz")] if cache.exists() else []
    per = sum(have) / len(have) if have else None
    print(f"  nisar: {len(chain)} pairs from {how}; {len(chain) - len(todo)} cached, {len(todo)} to read"
          + (f" (~{_mb(per * len(todo))} on disk at this plate's window size)" if per and todo else ""))
    if a.dry_run or not todo:
        return 0
    if not nisar._token():
        print("  !! no Earthdata token: see `uv run sars doctor`")
        return 1
    if (msg := _path_problem(cache)):
        print(f"  !! {msg}")
        return 1
    cache.mkdir(parents=True, exist_ok=True)
    bad = 0
    for i, g in enumerate(todo, 1):
        t = time.time()
        try:
            with nisar.open_h5(g.url) as h5:
                raw = nisar.read_raw(h5, tuple(cfg["bbox"]), g.days, height_m=cfg["cache_height"])
            tmp = cache / f"_part_{i}.npz"                         # short: the final name is long
            raw.save(tmp)
            tmp.replace(cache / f"{g.name}.npz")                   # a broken run never leaves a half file behind
            print(f"    [{i}/{len(todo)}] {g.ref_date:%Y-%m-%d} -> {g.sec_date:%Y-%m-%d}  "
                  f"{_mb((cache / f'{g.name}.npz').stat().st_size)}  {time.time() - t:.0f}s", flush=True)
        except Exception as e:                                   # noqa: BLE001  (one bad product must not stop the rest)
            bad += 1
            print(f"    [{i}/{len(todo)}] FAILED {g.name[:60]}...: {e}", flush=True)
    return 1 if bad else 0


def fetch_dem(cfg, a) -> int:
    from hachure import dem
    if cfg.get("dem_path"):
        p = ROOT / cfg["dem_path"]
    elif "dem_res" not in cfg:
        print("  dem: not used by this plate")
        return 0
    else:
        p = ROOT / "data/raw" / f"copdem_{cfg.get('dem_name', cfg['id'])}_{int(cfg['dem_res'])}m.tif"
    if p.exists():
        print(f"  dem: {p.relative_to(ROOT)} present")
        return 0
    print(f"  dem: Copernicus GLO-30 window -> {p.relative_to(ROOT)}" + (" (dry run)" if a.dry_run else ""))
    if not a.dry_run:
        p.parent.mkdir(parents=True, exist_ok=True)
        dem.fetch_copdem(tuple(cfg["bbox"]), p, res=cfg["dem_res"])
    return 0


def fetch_osm(cfg, a) -> int:
    srcs = [c.get("src", "") for c in cfg.get("context", [])]
    if not cfg.get("osm") and not any("data/raw/osm/" in s for s in srcs):
        return 0
    missing = [s for s in srcs if "data/raw/osm/" in s and not (ROOT / s).exists()]
    print(f"  osm: Overpass layers -> data/raw/osm/{cfg['id']}/ ({len(missing)} missing)" + (" (dry run)" if a.dry_run else ""))
    if not a.dry_run and (missing or not (ROOT / "data/raw/osm" / cfg["id"]).exists()):
        _run_script("fetch_osm", cfg["_path"])                      # skips files that already exist
    return 0


def fetch_landcover(cfg, a) -> int:
    lc = cfg.get("landcover")
    if not lc:
        return 0
    p = ROOT / lc["src"] if lc.get("src") else None
    if p and p.exists():
        print(f"  worldcover: {p.relative_to(ROOT)} present")
        return 0
    print("  worldcover: ESA WorldCover 10 m window" + (" (dry run)" if a.dry_run else ""))
    if not a.dry_run:
        _run_script("fetch_worldcover", cfg["_path"])
    return 0


def _download(url: str, dst: Path) -> None:
    import requests
    dst.parent.mkdir(parents=True, exist_ok=True)
    tmp = dst.with_name(dst.name + ".part")
    with requests.get(url, headers={"User-Agent": "space-apps-research"}, stream=True, timeout=120) as r:
        r.raise_for_status()
        with open(tmp, "wb") as f:
            for chunk in r.iter_content(1 << 20):
                f.write(chunk)
    tmp.replace(dst)


def fetch_context(cfg, a) -> int:
    """Context files that carry a "url" in the site config (e.g. the CDMX risk atlas layers)."""
    todo = [c for c in cfg.get("context", []) if c.get("url") and c.get("src") and not (ROOT / c["src"]).exists()]
    bad = 0
    for c in todo:
        print(f"  context: {c['url']} -> {c['src']}" + (" (dry run)" if a.dry_run else ""))
        if not a.dry_run:
            try:
                _download(c["url"], ROOT / c["src"])
            except Exception as e:                                   # noqa: BLE001
                bad += 1
                print(f"    FAILED: {e}")
    return 1 if bad else 0


NGL_SERIES = "https://geodesy.unr.edu/gps_timeseries/IGS20/tenv3/IGS20/{}.tenv3"


def fetch_gnss(cfg, a) -> int:
    """GNSS daily series (Nevada Geodetic Laboratory) for the stations a plate checks against."""
    bad = 0
    for g in cfg.get("gnss", []):
        if not g.get("series") or (ROOT / g["series"]).exists():
            continue
        url = NGL_SERIES.format(g["station"])
        print(f"  gnss: {g['station']} -> {g['series']}" + (" (dry run)" if a.dry_run else ""))
        if not a.dry_run:
            try:
                _download(url, ROOT / g["series"])
            except Exception as e:                                   # noqa: BLE001
                bad += 1
                print(f"    FAILED: {e}")
    return bad and 1


def cmd_fetch(a):
    cfg = _site(a.site)
    print(f"{cfg['id']}:")
    parts = {"dem": fetch_dem, "nisar": fetch_nisar, "osm": fetch_osm, "worldcover": fetch_landcover,
             "context": fetch_context, "gnss": fetch_gnss}
    only = a.only.split(",") if a.only else list(parts)
    rc = 0
    for k in only:
        rc |= parts[k](cfg, a)
    other = [c["src"] for c in cfg.get("context", []) if c.get("src") and "data/raw/osm/" not in c["src"]
             and not (ROOT / c["src"]).exists()]
    for s in other:
        print(f"  !! context file not fetched by sars (get it by hand, see docs/data-guide.md): {s}")
    lc = (cfg.get("fields") or {}).get("landcover") or {}
    for key in ("npz", "grid_tif"):
        if lc.get(key) and not (ROOT / lc[key]).exists():
            print(f"  !! fields.landcover.{key} missing: {lc[key]} (see docs/data-guide.md, Central Valley)")
    return rc


# ---------------------------------------------------------------- build / story / pack / serve / new
def cmd_build(a):
    args = [str(_site(a.site)["_path"])] + (["--rebuild"] if a.rebuild else []) + (["--keep-terrain"] if a.keep_terrain else [])
    _run_script("build_site", *args)


def cmd_story(a):
    _run_script("story", str(_site(a.site)["_path"]))


def cmd_pack(a):
    _run_script("pack_scene", a.scene, "--theme", a.theme, *a.extra)


def cmd_serve(a):
    import functools
    import http.server
    h = functools.partial(http.server.SimpleHTTPRequestHandler, directory=str(ROOT / "web"))
    print(f"http://localhost:{a.port}/?scene=mexico-city-asc&theme=scroll   (Ctrl+C to stop)")
    http.server.ThreadingHTTPServer(("127.0.0.1", a.port), h).serve_forever()


TEMPLATE = {
    "id": None, "title": {"en": None, "uk": None}, "subtitle": "", "bbox": None, "direction": "ASCENDING",
    "cache_height": 0, "anchor": "radar-only", "dem_res": 30,
    "terrain": {"interval_m": 50, "spacing_px": 3.0, "sigma_dir_m": 250, "sigma_tone_m": 120},
    "layers": ["terrain", "velocity", "grid"],
    "notes": ["New plate: check the reference area (ref_bbox) before trusting any rate."],
}


def cmd_new(a):
    p = ROOT / "sites" / f"{a.id}.json"
    if p.exists():
        print(f"{p} exists")
        return 1
    cfg = json.loads(json.dumps(TEMPLATE))
    cfg.update(id=a.id, bbox=a.bbox, direction=a.direction, cache_height=a.height)
    cfg["title"] = {"en": a.title or a.id, "uk": a.title or a.id}
    if a.ref_bbox:
        cfg["ref_bbox"] = a.ref_bbox
    p.write_text(json.dumps(cfg, indent=1, ensure_ascii=False), encoding="utf-8")
    print(f"wrote {p.relative_to(ROOT)}. Next: uv run sars search {a.id}; uv run sars fetch {a.id}; uv run sars build {a.id}")
    return 0


def main(argv=None):
    for s in (sys.stdout, sys.stderr):                           # Windows consoles: names and arrows print, not crash
        try:
            s.reconfigure(encoding="utf-8")
        except Exception:                                        # noqa: BLE001
            pass
    ap = argparse.ArgumentParser(prog="sars", description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    d = sub.add_parser("doctor", help="check the environment and the Earthdata token")
    d.add_argument("--deep", action="store_true", help="also ask ASF to accept the token")
    sub.add_parser("sites", help="list plates and what is cached")
    for name in ("search", "fetch"):
        s = sub.add_parser(name, help="NISAR pairs over a plate" if name == "search" else "download what a plate needs")
        s.add_argument("site")
        s.add_argument("--start")
        s.add_argument("--end")
        s.add_argument("--all-pairs", action="store_true", help="every pair of the track, not only one 12-day chain")
        if name == "search":
            s.add_argument("-v", "--verbose", action="store_true")
        else:
            s.add_argument("--dry-run", action="store_true")
            s.add_argument("--update", action="store_true", help="search ASF again and rewrite the lock")
            s.add_argument("--max-pairs", type=int)
            s.add_argument("--only", help="comma list of: dem,nisar,osm,worldcover,context,gnss")
    lk = sub.add_parser("lock", help="pin the cached pairs of a plate (run where the plate was built)")
    lk.add_argument("site")
    b = sub.add_parser("build", help="build the web layers of a plate")
    b.add_argument("site")
    b.add_argument("--rebuild", action="store_true")
    b.add_argument("--keep-terrain", action="store_true")
    s = sub.add_parser("story", help="republish a plate's pages")
    s.add_argument("site")
    p = sub.add_parser("pack", help="one-file share pack")
    p.add_argument("scene")
    p.add_argument("--theme", default="scroll")
    p.add_argument("extra", nargs=argparse.REMAINDER)
    v = sub.add_parser("serve", help="serve web/ locally")
    v.add_argument("--port", type=int, default=8765)
    n = sub.add_parser("new", help="a new plate from the template")
    n.add_argument("id")
    n.add_argument("--bbox", nargs=4, type=float, required=True, metavar=("W", "S", "E", "N"))
    n.add_argument("--direction", choices=["ASCENDING", "DESCENDING"], default="ASCENDING")
    n.add_argument("--height", type=float, default=0, help="typical terrain height, m (incidence angle)")
    n.add_argument("--ref-bbox", nargs=4, type=float)
    n.add_argument("--title")
    a = ap.parse_args(argv)
    fn = {"doctor": cmd_doctor, "sites": cmd_sites, "search": cmd_search, "fetch": cmd_fetch, "lock": cmd_lock, "build": cmd_build,
          "story": cmd_story, "pack": cmd_pack, "serve": cmd_serve, "new": cmd_new}[a.cmd]
    return fn(a) or 0


if __name__ == "__main__":
    sys.exit(main())
