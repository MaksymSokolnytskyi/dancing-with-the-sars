"""Extract infrastructure around Svartsengi from the Geofabrik Iceland extract
(data/osm/iceland-latest.osm.pbf) with pyosmium -> data/osm/reykjanes_layers.json.

Layers: roads, power lines, pipelines (steam / hot water / other), barriers and dykes
(the lava protection walls), buildings (centroids), places.

Usage: python scripts/reykjanes/osm_extract.py
"""
import json
from pathlib import Path

import osmium

ROOT = Path(__file__).resolve().parents[2]
PBF = ROOT / "data" / "osm" / "iceland-latest.osm.pbf"
OUT = ROOT / "data" / "osm" / "reykjanes_layers.json"
BOX = (-22.80, 63.76, -22.05, 64.02)  # lon0, lat0, lon1, lat1 (the NISAR window)
ROADS = {"motorway", "trunk", "primary", "secondary", "tertiary", "unclassified", "residential", "service", "track"}


def inside(lon, lat):
    return BOX[0] <= lon <= BOX[2] and BOX[1] <= lat <= BOX[3]


class Handler(osmium.SimpleHandler):
    def __init__(self):
        super().__init__()
        self.layers = {k: [] for k in ("roads", "power", "pipelines", "barriers", "buildings", "places",
                                       "volcanoes", "vents", "peaks", "sites")}

    def point(self, t, lon, lat):
        """Volcanic features and key sites, from nodes or way centroids."""
        name = t.get("name", "")
        base = dict(name=name, name_en=t.get("name:en", ""), lon=lon, lat=lat, ele=t.get("ele", ""))
        if t.get("natural") == "volcano":
            self.layers["volcanoes"].append(dict(base, status=t.get("volcano:status", ""), type=t.get("volcano:type", "")))
        elif t.get("geological") == "volcanic_vent":
            self.layers["vents"].append(base)
        elif t.get("natural") == "peak" and name:
            self.layers["peaks"].append(base)
        kind = None
        near_svartsengi = abs(lat - 63.879) < 0.01 and abs(lon + 22.433) < 0.02  # the plant itself is untagged
        if t.get("power") == "plant" or (t.get("power") == "substation" and (name or near_svartsengi)):
            kind = "power plant" if t.get("power") == "plant" else "substation"
        elif t.get("aeroway") == "aerodrome":
            kind = "airport"
        elif t.get("leisure") == "resort" or (t.get("tourism") == "attraction" and "lón" in name):
            kind = "resort"
        elif t.get("tourism") == "hotel" and name:
            kind = "hotel"
        if kind:
            self.layers["sites"].append(dict(base, kind=kind, operator=t.get("operator", "")))

    def node(self, n):
        if not n.location.valid() or not inside(n.location.lon, n.location.lat):
            return
        if "place" in n.tags:
            self.layers["places"].append(dict(name=n.tags.get("name", ""), place=n.tags["place"],
                                              lon=n.location.lon, lat=n.location.lat))
        self.point(n.tags, n.location.lon, n.location.lat)

    def way(self, w):
        t = w.tags
        if any(k in t for k in ("power", "aeroway", "leisure", "tourism", "natural", "geological")) and \
                (t.get("power") in ("plant", "substation") or t.get("aeroway") == "aerodrome" or
                 t.get("leisure") == "resort" or t.get("natural") == "volcano" or "geological" in t):
            try:
                pts = [(nd.lon, nd.lat) for nd in w.nodes if nd.location.valid()]
            except osmium.InvalidLocationError:
                pts = []
            if pts:
                lon, lat = sum(p[0] for p in pts) / len(pts), sum(p[1] for p in pts) / len(pts)
                if inside(lon, lat):
                    self.point(t, lon, lat)
        kind = None
        if t.get("highway") in ROADS:
            kind = "roads"
        elif t.get("power") in ("line", "minor_line", "cable"):
            kind = "power"
        elif t.get("man_made") == "pipeline":
            kind = "pipelines"
        elif t.get("man_made") in ("dyke", "embankment") or t.get("barrier") in ("wall", "retaining_wall", "fence") \
                and ("lava" in (t.get("name", "") + t.get("description", "") + t.get("note", "")).lower() or
                     t.get("barrier") != "fence"):
            kind = "barriers"
        elif "building" in t:
            kind = "buildings"
        if kind is None:
            return
        try:
            pts = [(nd.lon, nd.lat) for nd in w.nodes if nd.location.valid()]
        except osmium.InvalidLocationError:
            return
        if not pts or not any(inside(*p) for p in pts):
            return
        props = {k: t.get(k) for k in ("name", "highway", "substance", "operator", "voltage", "man_made", "barrier",
                                       "building") if t.get(k)}
        if kind == "buildings":
            lon = sum(p[0] for p in pts) / len(pts)
            lat = sum(p[1] for p in pts) / len(pts)
            self.layers[kind].append(dict(lon=lon, lat=lat, **props))
        else:
            self.layers[kind].append(dict(coords=[[round(a, 6), round(b, 6)] for a, b in pts], **props))


if __name__ == "__main__":
    h = Handler()
    h.apply_file(str(PBF), locations=True)
    OUT.write_text(json.dumps(h.layers))
    for k, v in h.layers.items():
        print(k, len(v))
    print("barrier names:", sorted({b.get("name", "") for b in h.layers["barriers"]})[:20])
    print("sites:", [(s["kind"], s["name"]) for s in h.layers["sites"] if s["kind"] != "hotel"])
    print("vents:", [(round(v["lon"], 3), round(v["lat"], 3), v["name"]) for v in h.layers["vents"]])
    print("volcanoes:", [(v["name"], v["status"]) for v in h.layers["volcanoes"] if v["name"]])
    print("pipeline substances:", sorted({p.get("substance", "?") for p in h.layers["pipelines"]}))
