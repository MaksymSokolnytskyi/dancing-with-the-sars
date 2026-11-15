// OSM context for the Ostrava–Karviná / Upper Silesian subsidence bowls (NISAR vertical clusters).
// Plain Node (>= 18, built-in fetch), no npm packages.
//
// In : output/ostrava_karvina/zones_vertical_clustered.geojson
// Out: .../zones_vertical_osm.geojson, .../zones_vertical_osm.csv, .../osm_summary.json
// Raw Overpass answers are cached in data/ostrava/osm_cache/ (delete to refresh).
//
// Per bowl: buildings (count, footprint m², by centroid-in-polygon), roads km by class and rail km
// (segments densified to 10 m, midpoint-in-polygon), landuse composition (25 m sample grid),
// nearest mine object, nearest named settlement, country (OSM admin_level=2 via is_in).
// Usage: node scripts/ostrava/ostrava_osm.mjs
import fs from "node:fs";
import path from "node:path";
import { fileURLToPath } from "node:url";

const ROOT = path.resolve(path.dirname(fileURLToPath(import.meta.url)), "..", "..");   // the repo root
const DIR = path.join(ROOT, "output", "ostrava_karvina");
// "Reliable" = a HAND-PICKED list (post-hoc) from the first notebook run, 14.11.2026: bowls confirmed by 3 tracks and judged
// reliable by eye. M11, M12, M16, M17, M21 and M22 are also confirmed by 3 tracks but were kept as candidates. The rest are
// candidates. See notes/experiments.md, block O.
const ROBUST = new Set(["M01", "M02", "M03", "M04", "M05", "M06", "M07", "M08", "M09", "M18"]);
const CACHE = path.join(ROOT, "data", "ostrava", "osm_cache");
const BBOX = [49.70, 18.15, 50.05, 18.80]; // S W N E
const OVERPASS = "https://overpass-api.de/api/interpreter";
fs.mkdirSync(CACHE, { recursive: true });

// ---------- Overpass ----------
async function overpass(name, q) {
  const f = path.join(CACHE, name + ".json");
  if (fs.existsSync(f)) return JSON.parse(fs.readFileSync(f, "utf8"));
  for (let attempt = 0; attempt < 8; attempt++) {
    const r = await fetch(OVERPASS, {
      method: "POST",
      headers: { "Content-Type": "application/x-www-form-urlencoded", "User-Agent": "space-apps-research", Accept: "application/json" },
      body: "data=" + encodeURIComponent(q),
    });
    const t = await r.text();
    if (r.ok && t.startsWith("{")) { fs.writeFileSync(f, t); return JSON.parse(t); }
    console.log(`  overpass ${name}: HTTP ${r.status}, retry in ${10 * (attempt + 1)} s`);
    await new Promise((res) => setTimeout(res, 10000 * (attempt + 1)));
  }
  throw new Error("overpass failed: " + name);
}

// ---------- geometry in local metres ----------
const LAT0 = 49.9, KX = 111320 * Math.cos((LAT0 * Math.PI) / 180), KY = 110574;
const P = (lon, lat) => [lon * KX, lat * KY];
const ringArea = (r) => { let a = 0; for (let i = 0, j = r.length - 1; i < r.length; j = i++) a += (r[j][0] + r[i][0]) * (r[j][1] - r[i][1]); return Math.abs(a) / 2; };
function inRing(p, r) {
  let c = false;
  for (let i = 0, j = r.length - 1; i < r.length; j = i++) {
    if ((r[i][1] > p[1]) !== (r[j][1] > p[1]) && p[0] < ((r[j][0] - r[i][0]) * (p[1] - r[i][1])) / (r[j][1] - r[i][1]) + r[i][0]) c = !c;
  }
  return c;
}
// poly = array of polygons, each [outer, ...holes] in metres
const inPoly = (p, poly) => poly.some((rings) => inRing(p, rings[0]) && !rings.slice(1).some((h) => inRing(p, h)));
function segDist(p, a, b) {
  const dx = b[0] - a[0], dy = b[1] - a[1], L = dx * dx + dy * dy;
  const t = L ? Math.max(0, Math.min(1, ((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / L)) : 0;
  return Math.hypot(p[0] - a[0] - t * dx, p[1] - a[1] - t * dy);
}
function distToPoly(p, poly) {
  if (inPoly(p, poly)) return 0;
  let d = Infinity;
  for (const rings of poly) for (const r of rings) for (let i = 1; i < r.length; i++) d = Math.min(d, segDist(p, r[i - 1], r[i]));
  return d;
}
function lengthInside(pts, poly) { // pts in metres
  let L = 0;
  for (let i = 1; i < pts.length; i++) {
    const a = pts[i - 1], b = pts[i], s = Math.hypot(b[0] - a[0], b[1] - a[1]), n = Math.max(1, Math.ceil(s / 10));
    for (let k = 0; k < n; k++) {
      const t = (k + 0.5) / n;
      if (inPoly([a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])], poly)) L += s / n;
    }
  }
  return L;
}
const geomToPoly = (g) => (g.type === "Polygon" ? [g.coordinates] : g.coordinates).map((rings) => rings.map((r) => r.map(([x, y]) => P(x, y))));
const wayPts = (el) => (el.geometry || []).map((n) => P(n.lon, n.lat));
function elRings(el) { // closed outer rings of a way or multipolygon relation (coarse: closed outer members only)
  if (el.type === "way") { const r = wayPts(el); return r.length > 3 ? [r] : []; }
  if (el.type === "relation") return (el.members || []).filter((m) => m.role === "outer" && m.geometry && m.geometry.length > 3)
    .map((m) => m.geometry.map((n) => P(n.lon, n.lat))).filter((r) => Math.hypot(r[0][0] - r.at(-1)[0], r[0][1] - r.at(-1)[1]) < 1);
  return [];
}
const centroid = (r) => { let x = 0, y = 0; for (const p of r) { x += p[0]; y += p[1]; } return [x / r.length, y / r.length]; };
const elCenter = (el) => (el.center ? P(el.center.lon, el.center.lat) : el.lat != null ? P(el.lon, el.lat) : null);

// ---------- classes ----------
const ROAD_CLASS = { motorway: "major", motorway_link: "major", trunk: "major", trunk_link: "major", primary: "major", primary_link: "major",
  secondary: "major", secondary_link: "major", tertiary: "major", tertiary_link: "major",
  residential: "local", unclassified: "local", living_street: "local", service: "local", road: "local" };
const RAIL = new Set(["rail", "light_rail", "tram", "narrow_gauge", "subway"]);
function landClass(t) {
  if (t.man_made === "spoil_heap" || t.landuse === "landfill") return "spoil heap / landfill";
  if (t.landuse === "construction") return "construction site";
  if (t.natural === "water" || t.landuse === "reservoir" || t.landuse === "basin" || t.water) return "water";
  if (t.natural === "wetland") return "wetland";
  if (t.landuse === "industrial" || t.industrial || t.landuse === "railway") return "industrial";
  if (t.landuse === "quarry") return "quarry";
  if (t.landuse === "residential" || t.landuse === "retail" || t.landuse === "commercial") return "residential / commercial";
  if (t.landuse === "farmland" || t.landuse === "meadow" || t.landuse === "grass" || t.landuse === "farmyard" || t.landuse === "orchard" || t.natural === "grassland") return "farmland / meadow";
  if (t.landuse === "forest" || t.natural === "wood" || t.natural === "scrub") return "forest / scrub";
  if (t.leisure || t.landuse === "recreation_ground" || t.landuse === "allotments" || t.landuse === "cemetery") return "park / allotments / cemetery";
  if (t.landuse === "brownfield" || t.landuse === "greenfield") return "brownfield";
  return null;
}
const LAND_PRIORITY = ["spoil heap / landfill", "construction site", "water", "wetland", "quarry", "industrial", "residential / commercial", "park / allotments / cemetery", "brownfield", "forest / scrub", "farmland / meadow"];

// ---------- main ----------
const zones = JSON.parse(fs.readFileSync(path.join(DIR, "zones_vertical_clustered.geojson"), "utf8"));
const bb = BBOX.join(",");
const mines = await overpass("mines", `[out:json][timeout:120];(
  nwr["man_made"="mineshaft"](${bb}); nwr["industrial"="mine"](${bb}); nwr["landuse"="industrial"]["name"~"KWK|[Kk]opalni|Důl|ZG |Ruch |[Mm]ine|Doly"](${bb}); );out center tags;`);
const places = await overpass("places", `[out:json][timeout:120];nwr["place"~"^(city|town|village|suburb|quarter|hamlet)$"]["name"](${bb});out center tags;`);
const mineList = mines.elements.map((e) => ({ p: elCenter(e), name: e.tags.name || "", kind: e.tags.man_made === "mineshaft" ? "mine shaft" : e.tags.industrial === "mine" ? "mine" : "mine industrial area", id: `${e.type}/${e.id}`, disused: !!(e.tags.disused || e.tags["disused:man_made"] || e.tags.abandoned) })).filter((m) => m.p);
const placeList = places.elements.map((e) => ({ p: elCenter(e), name: e.tags.name, place: e.tags.place })).filter((m) => m.p);
console.log(`mines/shafts: ${mineList.length}, named settlements: ${placeList.length}`);
// country for all bowls in one Overpass request (one is_in block per centroid, outputs come in order)
const ctryKey = "countries_" + zones.features.map((f) => `${f.properties.lon.toFixed(4)}_${f.properties.lat.toFixed(4)}`).join("").replace(/\./g, "").slice(0, 40) + zones.features.length;
const ctryQ = "[out:json][timeout:300];" + zones.features.map((f) => `is_in(${f.properties.lat},${f.properties.lon})->.a;area.a["admin_level"="2"];out tags;make sep;out;`).join("");
const ctry = []; let cur = null;
for (const el of (await overpass(ctryKey, ctryQ)).elements) {
  if (el.type === "sep") { ctry.push(cur || "?"); cur = null; } else cur = cur || el.tags?.["ISO3166-1"] || null;
}
if (ctry.length !== zones.features.length) throw new Error(`country answers ${ctry.length} != ${zones.features.length}`);
const buildingPts = []; // centroids of buildings inside bowls, for the map

for (const [zi, ft] of zones.features.entries()) {
  const pr = ft.properties, poly = geomToPoly(ft.geometry);
  const cs = ft.geometry.type === "Polygon" ? [ft.geometry.coordinates] : ft.geometry.coordinates;
  let [w, s, e, n] = [180, 90, -180, -90];
  for (const rings of cs) for (const [x, y] of rings[0]) { w = Math.min(w, x); e = Math.max(e, x); s = Math.min(s, y); n = Math.max(n, y); }
  const pad = 0.0007, b = [s - pad, w - pad, n + pad, e + pad].map((v) => v.toFixed(5)).join(",");
  const d = await overpass("bbox_" + b.replace(/[.,]/g, "_"), `[out:json][timeout:180];(
    nwr["building"](${b}); way["highway"](${b}); way["railway"](${b});
    nwr["landuse"](${b}); nwr["natural"~"water|wetland|wood|scrub|grassland"](${b}); nwr["man_made"="spoil_heap"](${b}); nwr["leisure"](${b});
    way["waterway"](${b}); nwr["man_made"~"mineshaft|spoil_heap"](${b}); );out geom;`);

  let nb = 0, barea = 0; const roads = { major: 0, local: 0 }; let rail = 0, railDisused = 0, waterway = 0;
  const landEls = [], railNames = new Set(), roadNames = new Set(), inside = [];
  for (const el of d.elements) {
    const t = el.tags || {};
    if (t.building) {
      const rings = elRings(el); if (!rings.length) continue;
      const c0 = centroid(rings[0]);
      if (inPoly(c0, poly)) { nb++; barea += rings.reduce((a, r) => a + ringArea(r), 0); buildingPts.push([+(c0[0] / KX).toFixed(6), +(c0[1] / KY).toFixed(6), pr.id]); }
      continue;
    }
    if (el.type === "way" && t.highway && ROAD_CLASS[t.highway]) {
      const L = lengthInside(wayPts(el), poly); roads[ROAD_CLASS[t.highway]] += L;
      if (L > 50 && ROAD_CLASS[t.highway] === "major") roadNames.add(t.ref || t.name || t.highway);
    }
    if (el.type === "way" && t.railway) {
      const L = lengthInside(wayPts(el), poly);
      if (RAIL.has(t.railway)) { rail += L; if (L > 50) railNames.add(t.name || t.usage || t.service || t.railway); }
      else if (t.railway === "disused" || t.railway === "abandoned") railDisused += L;
    }
    if (el.type === "way" && t.waterway) waterway += lengthInside(wayPts(el), poly);
    const lc = landClass(t);
    if (lc) for (const r of elRings(el)) landEls.push({ r, lc, bbox: r.reduce((a, p) => [Math.min(a[0], p[0]), Math.min(a[1], p[1]), Math.max(a[2], p[0]), Math.max(a[3], p[1])], [Infinity, Infinity, -Infinity, -Infinity]) });
    if (t.man_made === "mineshaft" || t.man_made === "spoil_heap" || t.landuse === "construction" || t.landuse === "landfill") {
      const c = elCenter(el) || (elRings(el)[0] && centroid(elRings(el)[0]));
      if (c && distToPoly(c, poly) < 1) inside.push(`${t.man_made || t.landuse}${t.name ? " «" + t.name + "»" : ""}`);
    }
  }
  // landuse composition on a 25 m grid
  const [X0, Y0] = P(w, s), [X1, Y1] = P(e, n), comp = {}; let ns = 0;
  for (let x = X0; x <= X1; x += 25) for (let y = Y0; y <= Y1; y += 25) {
    if (!inPoly([x, y], poly)) continue; ns++;
    let best = null;
    for (const L of landEls) if (x >= L.bbox[0] && x <= L.bbox[2] && y >= L.bbox[1] && y <= L.bbox[3] && inRing([x, y], L.r))
      if (!best || LAND_PRIORITY.indexOf(L.lc) < LAND_PRIORITY.indexOf(best)) best = L.lc;
    comp[best || "no landuse"] = (comp[best || "no landuse"] || 0) + 1;
  }
  const landuse = Object.entries(comp).sort((a, b) => b[1] - a[1]).map(([k, v]) => `${k} ${Math.round((100 * v) / ns)}%`).join(", ");

  const near = (list, f = () => true) => list.filter(f).map((m) => ({ ...m, d: distToPoly(m.p, poly) })).sort((a, b) => a.d - b.d)[0];
  const mine = near(mineList), mineNamed = near(mineList, (m) => m.name);
  const town = near(placeList, (p) => ["city", "town", "village"].includes(p.place)), local = near(placeList);
  const ctr = P(pr.lon, pr.lat);
  const country = ctry[zi];

  Object.assign(pr, {
    country, reliable: ROBUST.has(pr.id) ? "yes" : "candidate",
    buildings: nb, built_m2: Math.round(barea), built_frac: +(barea / (pr.area_km2 * 1e6)).toFixed(3),
    roads_km: +((roads.major + roads.local) / 1000).toFixed(2), roads_major_km: +(roads.major / 1000).toFixed(2), roads_local_km: +(roads.local / 1000).toFixed(2),
    major_roads: [...roadNames].slice(0, 6).join("; "),
    rail_km: +(rail / 1000).toFixed(2), rail_disused_km: +(railDisused / 1000).toFixed(2), rail_lines: [...railNames].slice(0, 4).join("; "),
    waterway_km: +(waterway / 1000).toFixed(2),
    landuse,
    osm_inside: [...new Set(inside)].slice(0, 6).join("; "),
    mine: mineNamed ? mineNamed.name : "", mine_kind: mineNamed ? mineNamed.kind : "", mine_dist_km: mineNamed ? +(mineNamed.d / 1000).toFixed(2) : null, mine_osm: mineNamed ? mineNamed.id : "",
    nearest_mine_obj_km: mine ? +(mine.d / 1000).toFixed(2) : null, nearest_mine_obj: mine ? `${mine.kind}${mine.name ? " " + mine.name : ""}${mine.disused ? " (disused)" : ""}` : "",
    settlement: town ? town.name : "", settlement_km: town ? +(Math.hypot(town.p[0] - ctr[0], town.p[1] - ctr[1]) / 1000).toFixed(1) : null,
    locality: local ? local.name : "", locality_place: local ? local.place : "",
  });
  console.log(`${pr.id} ${country} ${pr.area_km2} km² | bldg ${nb} | roads ${pr.roads_km} km | rail ${pr.rail_km} km | mine: ${pr.mine || "—"} ${pr.mine_dist_km ?? ""} km | ${pr.locality} | ${landuse}`);
}

// ---------- outputs ----------
fs.writeFileSync(path.join(DIR, "zones_vertical_osm.geojson"), JSON.stringify(zones));
const cols = Object.keys(zones.features[0].properties);
const csvEsc = (v) => (v == null ? "" : /[",\n]/.test(String(v)) ? `"${String(v).replace(/"/g, '""')}"` : String(v));
fs.writeFileSync(path.join(DIR, "zones_vertical_osm.csv"), "﻿" + [cols.join(","), ...zones.features.map((f) => cols.map((c) => csvEsc(f.properties[c])).join(","))].join("\n"));

fs.writeFileSync(path.join(DIR, "osm_buildings_in_bowls.json"), JSON.stringify(buildingPts));
const summary = {};
for (const [label, filt] of [["all", () => true], ["robust", (p) => ROBUST.has(p.id)], ["candidates", (p) => !ROBUST.has(p.id)]]) {
  for (const c of ["PL", "CZ"]) {
    const ps = zones.features.map((f) => f.properties).filter((p) => p.country === c && filt(p));
    summary[`${label}_${c}`] = {
      zones: ps.length, ids: ps.map((p) => p.id).join(" "),
      area_km2: +ps.reduce((a, p) => a + p.area_km2, 0).toFixed(2),
      buildings: ps.reduce((a, p) => a + p.buildings, 0), built_m2: ps.reduce((a, p) => a + p.built_m2, 0),
      roads_km: +ps.reduce((a, p) => a + p.roads_km, 0).toFixed(1), rail_km: +ps.reduce((a, p) => a + p.rail_km, 0).toFixed(1),
      up_mean_area_weighted_cm_yr: ps.length ? +(ps.reduce((a, p) => a + p.up_mean_cm_yr * p.area_km2, 0) / ps.reduce((a, p) => a + p.area_km2, 0)).toFixed(1) : null,
    };
  }
}
// historical bowls: parts of one former cluster share `group`
const groups = {};
for (const p of zones.features.map((f) => f.properties)) {
  const g = (groups[p.group] ??= { group: p.group, ids: [], country: p.country, area_km2: 0, buildings: 0, roads_km: 0, rail_km: 0, _w: 0, reliable_parts: 0 });
  g.ids.push(p.id); g.area_km2 += p.area_km2; g.buildings += p.buildings; g.roads_km += p.roads_km; g.rail_km += p.rail_km;
  g._w += p.up_mean_cm_yr * p.area_km2; if (ROBUST.has(p.id)) g.reliable_parts++;
  if (g.country !== p.country) g.country += "/" + p.country;
}
const groupList = Object.values(groups).map((g) => ({ ...g, ids: g.ids.join(" "), area_km2: +g.area_km2.toFixed(2), roads_km: +g.roads_km.toFixed(2), rail_km: +g.rail_km.toFixed(2),
  up_mean_area_weighted_cm_yr: +(g._w / g.area_km2).toFixed(1), _w: undefined })).sort((a, b) => a.up_mean_area_weighted_cm_yr - b.up_mean_area_weighted_cm_yr);
summary.groups = groupList.filter((g) => g.ids.includes(" "));
fs.writeFileSync(path.join(DIR, "osm_summary.json"), JSON.stringify(summary, null, 1));
const gcols = ["group", "ids", "country", "reliable_parts", "area_km2", "buildings", "roads_km", "rail_km", "up_mean_area_weighted_cm_yr"];
fs.writeFileSync(path.join(DIR, "zones_vertical_osm_groups.csv"), "﻿" + [gcols.join(","), ...groupList.map((g) => gcols.map((c) => csvEsc(g[c])).join(","))].join("\n"));
const { groups: _g, ...flat } = summary;
console.table(flat); console.table(summary.groups);
