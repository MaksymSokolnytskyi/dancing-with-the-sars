// Results page for the Ostrava–Karviná study: OSM exposure of the NISAR bowls (watershed split, water mask),
// PL vs CZ, historical groups, spoil heaps, map. Plain Node (>= 18), no npm packages. A team report, in English.
// Out: output/ostrava_karvina/results_ostrava.html
// Usage: node scripts/ostrava/build_ostrava_results.mjs   (after ostrava_osm.mjs and ostrava_map.mjs)
import fs from "node:fs";
import path from "node:path";
import { fileURLToPath } from "node:url";

const ROOT = path.resolve(path.dirname(fileURLToPath(import.meta.url)), "..", "..");   // the repo root
const DIR = path.join(ROOT, "output", "ostrava_karvina");
const rd = (f) => JSON.parse(fs.readFileSync(path.join(DIR, f), "utf8"));
const Z = rd("zones_vertical_osm.geojson").features.map((f) => f.properties);
const S = rd("osm_summary.json"), MAP = rd("map_vertical_osm_labels.json");
const esc = (s) => String(s ?? "").replace(/[&<>"]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" })[c]);
// thousands with a comma, independent of the Node build's locale data
const n0 = (v) => (v == null || v === "" ? "—" : String(Number(v)).replace(/\B(?=(\d{3})+(?!\d))/g, ","));
const pct = (v) => (v == null || v === "" ? "—" : Math.round(v * 100) + "%");
const s = (k) => S[k] || { zones: 0, area_km2: 0, buildings: 0, built_m2: 0, roads_km: 0, rail_km: 0, up_mean_area_weighted_cm_yr: null, ids: "" };
// the period the numbers use: the pairs of the two tracks behind the vertical (A158_F027 + D152_F063)
const PERIOD = "24 Jun – 11 Oct 2026";
const PROV = `<span class="prov">provisional, not validated</span>`;

const sumRows = [["robust", "Reliable (hand-picked, post-hoc)"], ["candidates", "Candidates"], ["all", "All"]].flatMap(([k, lab]) => ["PL", "CZ"].map((c) => {
  const x = s(`${k}_${c}`);
  return `<tr><td>${lab}</td><td><b>${c}</b></td><td>${x.zones}</td><td>${x.area_km2}</td><td>${n0(x.buildings)}</td><td>${n0(Math.round(x.built_m2 / 1000))}</td><td>${x.roads_km}</td><td>${x.rail_km}</td><td>${x.up_mean_area_weighted_cm_yr ?? "—"}</td><td class="small mut">${esc(x.ids)}</td></tr>`;
})).join("");

// p10 near the unwrapping limit is shown as a lower bound ("≤"), faded
const p10 = (v) => (v != null && v !== "" && v < -100 ? `<span class="faint" title="near the L-band unwrapping limit: a lower bound, not a measurement">≤ ${v}</span>` : v);
const zoneRows = Z.map((p) => `<tr class="${p.reliable === "yes" ? "" : "cand"}"><td><b>${p.id}</b></td><td>${p.group}</td><td>${p.country}</td><td>${esc(p.locality)}</td>
<td>${p.area_km2}</td><td><b>${p.up_mean_cm_yr}</b></td><td>${p10(p.up_p10_cm_yr)}</td><td>${p.n_tracks_confirm}</td><td>${pct(p.no_data_frac)}</td><td>${pct(p.heap_frac)}</td><td>${pct(p.water_frac)}</td>
<td>${n0(p.buildings)}</td><td>${p.roads_km}</td><td>${p.rail_km}</td><td>${esc(p.mine) || "—"}${p.mine ? ` <span class="mut">${p.mine_dist_km} km</span>` : ""}</td></tr>`).join("");

const groupRows = (S.groups || []).map((g) => `<tr><td>${g.group}</td><td>${esc(g.ids)}</td><td>${g.country}</td><td>${g.area_km2}</td><td><b>${g.up_mean_area_weighted_cm_yr}</b></td><td>${n0(g.buildings)}</td><td>${g.roads_km}</td><td>${g.rail_km}</td></tr>`).join("");

const heapRows = Z.filter((p) => p.up_mean_on_heap !== "" && p.up_mean_on_heap != null && p.heap_frac > 0.01).map((p) =>
  `<tr><td>${p.id}</td><td>${pct(p.heap_frac)}</td><td>${p.up_mean_on_heap}</td><td>${p.up_mean_off_heap ?? "—"}</td></tr>`).join("");

const mapFig = (m, alt) => `<div class="map"><img src="data:image/png;base64,${fs.readFileSync(path.join(DIR, m.file)).toString("base64")}" alt="${esc(alt)}" width="${m.w}" height="${m.h}"><div class="lbl">`
  + (m.places || []).map((p) => `<span class="pl ${p.place}" style="left:${p.x}%;top:${p.y}%">${esc(p.name)}</span>`).join("")
  + (m.mines || []).map((p) => `<span class="mn" style="left:${p.x}%;top:${p.y}%" title="${esc(p.full || p.name)} (OSM label)">${esc(p.name)}</span>`).join("")
  + m.labels.map((l) => `<span class="${l.rel ? "lr" : "lc"}" style="left:${l.x}%;top:${l.y}%" title="${l.id}: ${l.v} cm/yr (${l.country})">${l.id}</span>`).join("") + `</div></div>`;

// ---------- OSM details ----------
const osmRows = Z.map((p) => `<tr class="${p.reliable === "yes" ? "" : "cand"}"><td><b>${p.id}</b></td><td>${esc(p.settlement)} <span class="mut">${p.settlement_km ?? ""} km</span><br><span class="small mut">${esc(p.locality)} (${esc(p.locality_place)})</span></td>
<td class="small">${esc(p.landuse)}</td><td class="small">${esc(p.major_roads) || "—"}</td><td class="small">${esc(p.rail_lines) || "—"}${p.rail_disused_km ? `<br><span class="mut">disused: ${p.rail_disused_km} km</span>` : ""}</td>
<td>${p.built_frac != null ? pct(p.built_frac) : "—"}</td><td>${p.waterway_km}</td><td class="small">${esc(p.nearest_mine_obj)} <span class="mut">${p.nearest_mine_obj_km ?? ""} km</span></td><td class="small">${esc(p.osm_inside) || "—"}</td></tr>`).join("");
// named mine objects from OSM and the nearest reliable bowl (distance between centroids)
const MINES = JSON.parse(fs.readFileSync(path.join(ROOT, "data", "ostrava", "osm_cache", "mines.json"), "utf8")).elements
  .filter((e) => e.tags.name && e.tags.man_made !== "mineshaft" && e.tags.industrial !== "salt_pond")
  .map((e) => ({ name: e.tags.name, lon: e.center ? e.center.lon : e.lon, lat: e.center ? e.center.lat : e.lat, kind: e.tags.industrial === "mine" ? "mine" : e.tags.landuse === "industrial" ? "mine industrial area" : e.tags.industrial || "" }));
const km = (a, b) => Math.hypot((a.lon - b.lon) * 111.32 * Math.cos(49.9 * Math.PI / 180), (a.lat - b.lat) * 110.57);
const REL = Z.filter((p) => p.reliable === "yes");
const mineRows = MINES.map((m) => { const n = REL.map((p) => ({ p, d: km(m, p) })).sort((a, b) => a.d - b.d)[0]; return { ...m, near: n.p.id, d: n.d, v: n.p.up_mean_cm_yr }; })
  .sort((a, b) => a.d - b.d).map((m) => `<tr><td>${esc(m.name)}</td><td>${m.kind}</td><td>${m.lat.toFixed(3)}, ${m.lon.toFixed(3)}</td><td>${m.near}</td><td>${m.d.toFixed(1)}</td><td>${m.v}</td></tr>`).join("");

// ---------- notebook figures (kept as separate files next to the page: 2.9–3 MB each) ----------
const FIGS = [
  ["map_vertical_east.png", "Vertical (A158 + D152) and east–west motion with bowl outlines. Vertical: negative = subsidence; horizontal motion in the bowls points to their centres."],
  ["timeseries_A158_F027.png", "LOS time series, track A158_F027 (ascending). Bowl 1 ≈ M01, bowl 2 ≈ M02, bowl 3 ≈ M07. This is LOS, not vertical."],
  ["timeseries_D152_F063.png", "LOS time series, track D152_F063 (descending)."],
  ["timeseries_A057_F027.png", "LOS time series, track A057_F027 (ascending, confirmation only)."],
  ["map_A158_F027.png", "Track A158_F027: LOS velocity, Isolation Forest score and anomalous zones."],
  ["map_D152_F063.png", "Track D152_F063: LOS velocity, Isolation Forest score and anomalous zones."],
  ["map_A057_F027.png", "Track A057_F027: LOS velocity, Isolation Forest score and anomalous zones."],
].filter(([f]) => fs.existsSync(path.join(DIR, f)));
const figHtml = FIGS.map(([f, cap]) => `<figure><a href="${f}" target="_blank"><img src="${f}" alt="${esc(cap)}" loading="lazy"></a><figcaption>${esc(cap)}</figcaption></figure>`).join("");
const legendGrad = MAP.legend.stops.map(([v, c]) => `${c} ${((v - MAP.legend.min) / (MAP.legend.max - MAP.legend.min)) * 100}%`).join(",");

const R = s("robust_PL"), RC = s("robust_CZ");
const top = Z.filter((p) => p.reliable === "yes").sort((a, b) => a.up_mean_cm_yr - b.up_mean_cm_yr)[0];
// CZ bowls confirmed by 3 tracks but not in the hand-picked set (the "0 reliable in CZ" caveat)
const cz3 = Z.filter((p) => p.country === "CZ" && p.n_tracks_confirm >= 3);

const html = `<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Ostrava–Karviná bowls</title>
<style>
:root{--bg:#f6f5f2;--card:#fff;--ink:#1d1d1b;--mut:#6b6a66;--line:#e2e0da;--acc:#b23a2b}
@media (prefers-color-scheme:dark){:root:not([data-theme=light]){--bg:#151514;--card:#1f1f1d;--ink:#ecebe7;--mut:#a09f99;--line:#34332f;--acc:#e0705f}}
:root[data-theme=dark]{--bg:#151514;--card:#1f1f1d;--ink:#ecebe7;--mut:#a09f99;--line:#34332f;--acc:#e0705f}
body{margin:0;background:var(--bg);color:var(--ink);font:15px/1.55 system-ui,Segoe UI,sans-serif}
main{max-width:1150px;margin:0 auto;padding:24px 16px}
h1{font-size:30px;margin:0 0 4px}h2{font-size:20px;margin:0 0 10px}h3{font-size:16px;margin:16px 0 6px}
section{background:var(--card);border:1px solid var(--line);border-radius:12px;padding:18px 20px;margin:16px 0}
.mut{color:var(--mut)}.small{font-size:12.5px}.wrap{overflow-x:auto}.faint{opacity:.55}
.prov{display:inline-block;font-size:11px;font-weight:600;color:var(--acc);border:1px dashed var(--acc);border-radius:4px;padding:0 5px;margin:2px 0}
table{border-collapse:collapse;width:100%;font-size:13.5px}th,td{border-bottom:1px solid var(--line);padding:6px 8px;text-align:left;vertical-align:top}
th{font-weight:600;color:var(--mut);white-space:nowrap}tr.cand td{opacity:.7}
.kpi{display:flex;gap:14px;flex-wrap:wrap}.kpi div{flex:1 1 220px;border:1px solid var(--line);border-radius:10px;padding:10px 14px}
.kpi b{display:block;font-size:26px;color:var(--acc)}
.map{position:relative;background:#fff;border:1px solid var(--line);border-radius:8px;overflow:hidden}
.map img{display:block;width:100%;height:auto}.map svg{display:none}
.lbl span{position:absolute;transform:translate(-50%,-135%);font:600 11px/1 system-ui;padding:2px 4px;border-radius:4px;white-space:nowrap;cursor:default}
.lbl .lr{background:#1d1d1b;color:#fff}.lbl .lc{background:rgba(255,255,255,.85);color:#5b5a56;border:1px solid #c9c7c0;font-weight:500}
.legend{display:flex;align-items:center;gap:10px;flex-wrap:wrap;font-size:13px;margin:10px 0 0}
.bar{width:260px;height:12px;border-radius:3px;background:linear-gradient(90deg,${legendGrad})}
.sw{display:inline-block;width:14px;height:10px;vertical-align:middle;margin-right:4px;border-radius:2px}
.lbl .pl{background:none;color:#1d1d1b;font-weight:500;text-shadow:0 0 3px #fff,0 0 3px #fff,0 0 2px #fff;transform:translate(-50%,-50%);pointer-events:none}
.lbl .pl.city,.lbl .pl.town{font-weight:700;font-size:12px}.lbl .pl.village{font-size:10px;color:#3a3936}
.lbl .mn{background:none;color:#111;font:italic 600 9.5px/1 system-ui;text-shadow:0 0 3px #fff,0 0 3px #fff;transform:translate(6px,-50%)}
ol.steps{padding-left:20px;margin:6px 0}ol.steps li{margin:4px 0}
.two{display:grid;grid-template-columns:1fr 1fr;gap:18px}@media (max-width:760px){.two{grid-template-columns:1fr}}
figure{margin:0 0 18px}figure img{width:100%;height:auto;border:1px solid var(--line);border-radius:8px;background:#fff}figcaption{font-size:13px;color:var(--mut);margin-top:4px}
</style></head><body><main>
<h1>Ostrava–Karviná: a mining border seen from space</h1>
<p class="mut">NISAR L2 GUNW (PROVISIONAL), pairs ${PERIOD}. Vertical from A158_F027 + D152_F063, confirmation by A057_F027. Bowls: threshold −15 cm/yr, OSM water mask, watershed split. OSM as of 14 Nov 2026. All rates are <b>vertical</b>, bowl means, cm/yr (negative = subsidence), unless stated otherwise. ${PROV}</p>

<section><h2>Problem</h2>
<p>The Ostrava–Karviná (CZ) and Jastrzębie–Rybnik (PL) districts of the Upper Silesian Coal Basin are a dense urban area with a long history of underground hard-coal mining. Over the mined-out longwalls the void slowly closes and the ground above it sinks, forming <b>mining subsidence bowls</b>: local subsidence from tens of centimetres to metres per year. The damage: cracked houses, bent roads and railways, broken utilities, flooded hollows and induced seismicity.</p>
<p>The idea for the atlas: <b>a border seen from space</b>. On the Czech side mining is winding down, so the ground should be calmer. On the Polish side the mines still work, so the ground should drop visibly. The questions: is this difference visible in the first NISAR data, how many buildings, roads and railways stand on sinking ground, and what is there: a mine, a spoil heap or water.</p></section>

<section><h2>Idea and method</h2>
<div class="two"><div>
<ol class="steps">
<li><b>Data.</b> All 12-day NISAR L2 GUNW interferograms (L-band, collection <code>NISAR_L2_GUNW_PROVISIONAL_V1</code>) over the box 18.15–18.80°E, 49.70–50.05°N: 57 granules, of which 3 frames are used (A057_F027, D152_F063, A158_F027).</li>
<li><b>Crops without downloading frames.</b> HTTPS range reads of the HDF5 (earthaccess + h5py), 80 m pixels.</li>
<li><b>Phase cleaning.</b> Coherence mask &lt; 0.3, swath edge, ionosphere screen, deramp, long-wavelength noise filter (σ = 5 km), reference to the scene median.</li>
<li><b>LOS velocity</b>: mean over pairs (displacement / dt), cm/yr.</li>
<li><b>Isolation Forest</b> on pixel features (velocity, pair scatter, gradient, coherence, share of valid pairs): an anomaly score.</li>
</ol></div><div>
<ol class="steps" start="6">
<li><b>Vertical and east–west</b> from an ascending + descending track pair (A158 + D152). LOS vectors from the product metadata (not flipped); NISAR looks <b>left</b>.</li>
<li><b>Bowls.</b> Vertical &lt; −15 cm/yr minus OSM water → add white centres (coherence loss) → merge fragments closer than 400 m → watershed split of bowls with several centres → keep bowls confirmed by at least 2 LOS tracks.</li>
<li><b>OSM context.</b> Country (admin_level=2), buildings, roads, railways, spoil heaps, mines, settlements, land use, for each bowl.</li>
<li><b>Consistency check.</b> Horizontal motion in a bowl should point to its centre (82 % of bowls). This is a consistency check, not an independent test: the method and its settings were developed on these data. Ponds, spoil heaps and linear features were checked separately.</li>
</ol></div></div>
<p class="small mut">Limits: PROVISIONAL data and only ~4 months, so slow millimetre trends are not visible yet. The fastest centres (more than ~12 cm LOS in 12 days) are near the limit of unambiguous L-band phase unwrapping, so their values are lower bounds. Mine names and mining status are OSM labels, not checked at the source. No EGMS validation yet. All settings were chosen on these same data (post-hoc, see notes/experiments.md, block O).</p></section>

<section><h2>Headline numbers ${PROV}</h2>
<div class="kpi">
<div><b>${n0(R.buildings)}</b>buildings stand in ${R.zones} reliable Polish bowls (${R.area_km2} km²), which sink by ${Math.abs(R.up_mean_area_weighted_cm_yr)} cm/yr on average (area-weighted). ${PROV}</div>
<div><b>${RC.zones}</b>reliable bowls on the Czech side${RC.zones ? ` (${n0(RC.buildings)} buildings)` : " after the water mask"}. The reliable set is hand-picked; ${cz3.length} Czech bowl${cz3.length === 1 ? " is" : "s are"} also confirmed by 3 tracks (${cz3.map((p) => p.id).join(", ") || "none"}).</div>
<div><b>${R.roads_km} km</b>of roads and <b style="display:inline;font-size:inherit">${R.rail_km} km</b> of railways in the reliable Polish bowls. ${PROV}</div>
</div>
<p class="small mut">Most active reliable bowl: ${top.id} (${esc(top.locality)}, ${top.country}): ${top.up_mean_cm_yr} cm/yr on average; its p10 (${top.up_p10_cm_yr} cm/yr) is a lower bound near the unwrapping limit. Peak values (~−2 m/yr LOS, ~−3 m/yr in the smoothed vertical) are not measurements, only provisional lower bounds until validation.</p>
</section>

<section><h2>PL vs CZ</h2><div class="wrap"><table>
<tr><th>Set</th><th>Country</th><th>Bowls</th><th>Area, km²</th><th>Buildings</th><th>Built-up, 1000 m²</th><th>Roads, km</th><th>Railways, km</th><th>Vertical (weighted)</th><th>id</th></tr>${sumRows}</table></div>
<p class="small mut">Country from the OSM admin_level=2 border at the bowl centroid. A building counts if its centroid is inside the polygon. Roads: motorway…tertiary + residential/unclassified/service. Lengths measured in 10 m steps. "Reliable" is a hand-picked set (post-hoc); a rule-based set is pre-registered in notes/experiments.md (O3a).</p></section>

<section><h2>Map</h2>
${mapFig(MAP.zoom, "Vertical rate, bowls, spoil heaps and buildings: the Jastrzębie–Rybnik cluster")}
<div class="legend"><span>−40</span><span class="bar"></span><span>+40 cm/yr</span>
<span><span class="sw" style="background:#8a6a3b88;border:1px solid #6b4f27"></span>spoil heap (OSM)</span>
<span><span class="sw" style="background:#1d1d1b;border-radius:50%;width:8px;height:8px"></span>building in a bowl</span>
<span><span class="sw" style="background:#fff;border:2px solid #1d1d1b"></span>reliable bowl (hand-picked)</span>
<span><span class="sw" style="background:#fff;border:1px dashed #5b5a56"></span>candidate</span>
<span><span class="sw" style="background:#fff;border:1px solid #ccc"></span>no data (in bowl centres: coherence lost to fast motion)</span>
<span><span class="sw" style="background:#55534e;height:3px"></span>roads (OSM)</span>
<span><span class="sw" style="background:repeating-linear-gradient(90deg,#2b2a27 0 4px,#fff 4px 7px);height:3px"></span>railways</span>
<span><span class="sw" style="background:#6d3fc0;height:3px"></span>CZ/PL border</span>
<span><span class="sw" style="background:#111;width:9px;height:9px;outline:1px solid #fff"></span>mine / mine industrial area (OSM label; small square = shaft)</span></div>
<h3>Whole area</h3>
${mapFig(MAP.full, "Vertical rate over the whole area 18.15–18.80°E, 49.70–50.05°N")}
</section>

<section><h2>Bowls ${PROV}</h2><div class="wrap"><table>
<tr><th>id</th><th>group</th><th>Country</th><th>Near</th><th>km²</th><th>Vert. mean</th><th>p10</th><th>Tracks</th><th>No data</th><th>Spoil heap</th><th>Water</th><th>Buildings</th><th>Roads, km</th><th>Rail, km</th><th>Nearest named mine (OSM label)</th></tr>${zoneRows}</table></div>
<p class="small mut">Faded rows are candidates: do not publish without validation (EGMS). p10 values below −100 cm/yr are shown as lower bounds (near the unwrapping limit). Mine names are OSM labels, not checked at the source.</p></section>

<section><h2>OSM: what is inside and around the bowls</h2><div class="wrap"><table>
<tr><th>id</th><th>Settlement (centre, km) / nearest place</th><th>Land use (share of area)</th><th>Major roads</th><th>Railways</th><th>Built-up</th><th>Waterways, km</th><th>Nearest mine object</th><th>OSM objects inside</th></tr>${osmRows}</table></div>
<p class="small mut">Land use from a 25 m grid inside the polygon (landuse / natural / man_made). "Built-up" is the share of the bowl under buildings. The nearest mine object includes shafts (man_made=mineshaft).</p>
<h3>Mines in OSM and the nearest reliable bowl</h3>
<div class="wrap"><table><tr><th>Name (OSM label)</th><th>Type</th><th>lat, lon</th><th>Nearest reliable bowl</th><th>Distance to its centroid, km</th><th>Its vertical, cm/yr</th></tr>${mineRows}</table></div>
<p class="small mut">Centroid distances are rough, without mining-lease boundaries. Statuses ("w likwidacji", "nieczynna", "bývalý") come from OSM names, not checked. Lease boundaries: PL — MIDAS PIG-PIB, CZ — ČGS geoportal (not done yet).</p></section>

<section><h2>Notebook figures</h2>
<p class="small mut">Files in <code>output/ostrava_karvina/</code>, made by <code>scripts/ostrava/ostrava_karvina.ipynb</code>. Click for full size. Time series are LOS, not vertical.</p>
${figHtml}</section>

<section><h2>Historical bowls (group)</h2>
<p>The watershed split some earlier clusters into several centres. Group totals are for comparison with earlier figures.</p>
<div class="wrap"><table><tr><th>group</th><th>Parts</th><th>Country</th><th>Area, km²</th><th>Vertical (weighted)</th><th>Buildings</th><th>Roads, km</th><th>Rail, km</th></tr>${groupRows}</table></div></section>

<section><h2>Spoil heaps: on the heap or next to it?</h2>
<div class="wrap"><table><tr><th>id</th><th>Heap share</th><th>Vertical on the heap</th><th>Vertical off the heap</th></tr>${heapRows}</table></div>
<p>Where a heap covers part of a bowl, the heap sinks more slowly than the ground next to it. So the strongest subsidence is not heap compaction but the ground over the workings.</p></section>

<script>
// hide map labels that collide with a higher-priority one (bowls > towns/cities > villages > mines); rerun on resize
function declutter(){document.querySelectorAll('.map .lbl').forEach(function(box){
  var rank=function(e){return e.classList.contains('lr')?0:e.classList.contains('lc')?1:e.classList.contains('city')||e.classList.contains('town')?2:e.classList.contains('village')?3:4;};
  var els=[].slice.call(box.children).sort(function(a,b){return rank(a)-rank(b);}),kept=[];
  els.forEach(function(e){e.style.visibility='';var r=e.getBoundingClientRect();
    var hit=kept.some(function(k){return r.left<k.right+2&&k.left<r.right+2&&r.top<k.bottom+1&&k.top<r.bottom+1;});
    if(hit)e.style.visibility='hidden';else kept.push(r);});});}
window.addEventListener('load',declutter);window.addEventListener('resize',declutter);
</script>
<p class="small mut">Scripts: scripts/ostrava/ (notebook, ostrava_osm.mjs, ostrava_map.mjs, build_ostrava_results.mjs). OSM data © OpenStreetMap contributors, ODbL. NISAR: NASA/ISRO, ASF DAAC, PROVISIONAL product.</p>
</main></body></html>`;
fs.writeFileSync(path.join(DIR, "results_ostrava.html"), html);
console.log("written results_ostrava.html");
