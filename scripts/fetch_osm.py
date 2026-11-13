"""Fetch a few OpenStreetMap context layers for a site's bbox via Overpass (stdlib only, polite retries).

    uv run python scripts/fetch_osm.py sites/mexico-city-asc.json

Writes data/raw/osm/<site id>/<name>.json (raw Overpass JSON with geometry). Which layers: the site's
"osm" block {name: overpass-ql body}, or the defaults below. Attribution: © OpenStreetMap contributors (ODbL).
"""
from __future__ import annotations

import json
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SERVERS = ["https://overpass-api.de/api/interpreter", "https://maps.mail.ru/osm/tools/overpass/api/interpreter",
           "https://overpass.kumi.systems/api/interpreter"]
DEFAULT = {
    "metro": 'relation["route"="subway"]({b});out geom;',
    "rail": 'way["railway"~"^(subway|light_rail|rail)$"]({b});out geom;',
    "runways": 'way["aeroway"="runway"]({b});out geom;',
    "canals": 'way["waterway"="canal"]({b});out geom;',
    "water": '(way["natural"="water"]({b});relation["natural"="water"]({b}););out geom;',
    "places": 'node["place"~"^(city|town|borough|suburb)$"]({b});out;',
    "stations": 'node["railway"="station"]["station"~"^(subway|light_rail)$"]({b});out;',
    # street network that cuts the built-up mask into blocks (Dufour/Siegfried town sign); residential streets left out
    "roads": 'way["highway"~"^(motorway|trunk|primary|secondary|tertiary)$"]({b});out tags geom;',
    "streets": 'way["highway"~"^(residential|living_street|unclassified|pedestrian)$"]({b});out skel geom;',   # geometry only
}


def query(body, timeout=180):
    data = urllib.parse.urlencode({"data": f"[out:json][timeout:{timeout}];{body}"}).encode()
    last = None
    for attempt in range(6):
        for url in SERVERS:
            try:
                req = urllib.request.Request(url, data=data, headers={"User-Agent": "space-apps-research"})
                with urllib.request.urlopen(req, timeout=timeout + 30) as r:
                    return json.loads(r.read().decode("utf-8"))
            except Exception as e:  # busy server, timeout, gateway error: try the next one, then back off
                last = e
        wait = 20 * 2 ** attempt
        print(f"  all servers busy ({last}); retry in {wait}s", flush=True)
        time.sleep(wait)
    raise RuntimeError(f"Overpass unavailable: {last}")


def main():
    cfg = json.loads(Path(sys.argv[1]).read_text(encoding="utf-8"))
    w, s, e, n = cfg["bbox"]
    b = f"{s},{w},{n},{e}"
    out = ROOT / "data/raw/osm" / cfg["id"]
    out.mkdir(parents=True, exist_ok=True)
    for name, body in (cfg.get("osm") or DEFAULT).items():
        dst = out / f"{name}.json"
        if dst.exists():
            print(f"{name}: cached")
            continue
        res = query(body.format(b=b))
        dst.write_text(json.dumps(res, ensure_ascii=False), encoding="utf-8")
        print(f"{name}: {len(res.get('elements', []))} elements, {dst.stat().st_size / 1e6:.1f} MB", flush=True)


if __name__ == "__main__":
    main()
