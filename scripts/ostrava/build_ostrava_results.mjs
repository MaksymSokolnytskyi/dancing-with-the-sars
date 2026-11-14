// Results page for the Ostrava–Karviná atlas place: OSM exposure of the NISAR bowls (watershed split, water mask),
// PL vs CZ, historical groups, spoil heaps, map. Plain Node (>= 18), no npm packages.
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
const n0 = (v) => (v == null || v === "" ? "—" : Number(v).toLocaleString("uk-UA"));
const pct = (v) => (v == null || v === "" ? "—" : Math.round(v * 100) + "%");
const s = (k) => S[k] || { zones: 0, area_km2: 0, buildings: 0, built_m2: 0, roads_km: 0, rail_km: 0, up_mean_area_weighted_cm_yr: null, ids: "" };

const sumRows = [["robust", "Надійні (відібрані вручну, post-hoc)"], ["candidates", "Кандидати"], ["all", "Усі"]].flatMap(([k, lab]) => ["PL", "CZ"].map((c) => {
  const x = s(`${k}_${c}`);
  return `<tr><td>${lab}</td><td><b>${c}</b></td><td>${x.zones}</td><td>${x.area_km2}</td><td>${n0(x.buildings)}</td><td>${n0(Math.round(x.built_m2 / 1000))}</td><td>${x.roads_km}</td><td>${x.rail_km}</td><td>${x.up_mean_area_weighted_cm_yr ?? "—"}</td><td class="small mut">${esc(x.ids)}</td></tr>`;
})).join("");

const zoneRows = Z.map((p) => `<tr class="${p.reliable === "так" ? "" : "cand"}"><td><b>${p.id}</b></td><td>${p.group}</td><td>${p.country}</td><td>${esc(p.locality)}</td>
<td>${p.area_km2}</td><td><b>${p.up_mean_cm_yr}</b></td><td>${p.up_p10_cm_yr}</td><td>${p.n_tracks_confirm}</td><td>${pct(p.no_data_frac)}</td><td>${pct(p.heap_frac)}</td><td>${pct(p.water_frac)}</td>
<td>${n0(p.buildings)}</td><td>${p.roads_km}</td><td>${p.rail_km}</td><td>${esc(p.mine) || "—"}${p.mine ? ` <span class="mut">${p.mine_dist_km} км</span>` : ""}</td></tr>`).join("");

const groupRows = (S.groups || []).map((g) => `<tr><td>${g.group}</td><td>${esc(g.ids)}</td><td>${g.country}</td><td>${g.area_km2}</td><td><b>${g.up_mean_area_weighted_cm_yr}</b></td><td>${n0(g.buildings)}</td><td>${g.roads_km}</td><td>${g.rail_km}</td></tr>`).join("");

const heapRows = Z.filter((p) => p.up_mean_on_heap !== "" && p.up_mean_on_heap != null && p.heap_frac > 0.01).map((p) =>
  `<tr><td>${p.id}</td><td>${pct(p.heap_frac)}</td><td>${p.up_mean_on_heap}</td><td>${p.up_mean_off_heap ?? "—"}</td></tr>`).join("");

const mapFig = (m, alt) => `<div class="map"><img src="data:image/png;base64,${fs.readFileSync(path.join(DIR, m.file)).toString("base64")}" alt="${esc(alt)}" width="${m.w}" height="${m.h}"><div class="lbl">`
  + (m.places || []).map((p) => `<span class="pl ${p.place}" style="left:${p.x}%;top:${p.y}%">${esc(p.name)}</span>`).join("")
  + (m.mines || []).map((p) => `<span class="mn" style="left:${p.x}%;top:${p.y}%" title="${esc(p.full || p.name)} (мітка OSM)">${esc(p.name)}</span>`).join("")
  + m.labels.map((l) => `<span class="${l.rel ? "lr" : "lc"}" style="left:${l.x}%;top:${l.y}%" title="${l.id}: ${l.v} см/рік (${l.country})">${l.id}</span>`).join("") + `</div></div>`;

// ---------- OSM details ----------
const osmRows = Z.map((p) => `<tr class="${p.reliable === "так" ? "" : "cand"}"><td><b>${p.id}</b></td><td>${esc(p.settlement)} <span class="mut">${p.settlement_km ?? ""} км</span><br><span class="small mut">${esc(p.locality)} (${esc(p.locality_place)})</span></td>
<td class="small">${esc(p.landuse)}</td><td class="small">${esc(p.major_roads) || "—"}</td><td class="small">${esc(p.rail_lines) || "—"}${p.rail_disused_km ? `<br><span class="mut">закинуті: ${p.rail_disused_km} км</span>` : ""}</td>
<td>${p.built_frac != null ? pct(p.built_frac) : "—"}</td><td>${p.waterway_km}</td><td class="small">${esc(p.nearest_mine_obj)} <span class="mut">${p.nearest_mine_obj_km ?? ""} км</span></td><td class="small">${esc(p.osm_inside) || "—"}</td></tr>`).join("");
// named mine objects from OSM and the nearest reliable bowl (distance between centroids)
const MINES = JSON.parse(fs.readFileSync(path.join(ROOT, "data", "ostrava", "osm_cache", "mines.json"), "utf8")).elements
  .filter((e) => e.tags.name && e.tags.man_made !== "mineshaft" && e.tags.industrial !== "salt_pond")
  .map((e) => ({ name: e.tags.name, lon: e.center ? e.center.lon : e.lon, lat: e.center ? e.center.lat : e.lat, kind: e.tags.industrial === "mine" ? "шахта" : e.tags.landuse === "industrial" ? "промзона шахти" : e.tags.industrial || "" }));
const km = (a, b) => Math.hypot((a.lon - b.lon) * 111.32 * Math.cos(49.9 * Math.PI / 180), (a.lat - b.lat) * 110.57);
const REL = Z.filter((p) => p.reliable === "так");
const mineRows = MINES.map((m) => { const n = REL.map((p) => ({ p, d: km(m, p) })).sort((a, b) => a.d - b.d)[0]; return { ...m, near: n.p.id, d: n.d, v: n.p.up_mean_cm_yr }; })
  .sort((a, b) => a.d - b.d).map((m) => `<tr><td>${esc(m.name)}</td><td>${m.kind}</td><td>${m.lat.toFixed(3)}, ${m.lon.toFixed(3)}</td><td>${m.near}</td><td>${m.d.toFixed(1)}</td><td>${m.v}</td></tr>`).join("");

// ---------- notebook figures (kept as separate files next to the page: 2.9–3 MB each) ----------
const FIGS = [
  ["map_vertical_east.png", "Вертикаль (A158 + D152) і рух схід–захід з контурами мульд. Вертикаль: від'ємне = осідання; горизонталь у центрах мульд спрямована до центру."],
  ["timeseries_A158_F027.png", "Часові ряди LOS, трек A158_F027 (висхідний). Мульда 1 ≈ M01, мульда 2 ≈ M02, мульда 3 ≈ M07. Це LOS, не вертикаль."],
  ["timeseries_D152_F063.png", "Часові ряди LOS, трек D152_F063 (низхідний)."],
  ["timeseries_A057_F027.png", "Часові ряди LOS, трек A057_F027 (висхідний, лише для підтвердження)."],
  ["map_A158_F027.png", "Трек A158_F027: швидкість LOS, бал Isolation Forest та аномальні зони."],
  ["map_D152_F063.png", "Трек D152_F063: швидкість LOS, бал Isolation Forest та аномальні зони."],
  ["map_A057_F027.png", "Трек A057_F027: швидкість LOS, бал Isolation Forest та аномальні зони."],
].filter(([f]) => fs.existsSync(path.join(DIR, f)));
const figHtml = FIGS.map(([f, cap]) => `<figure><a href="${f}" target="_blank"><img src="${f}" alt="${esc(cap)}" loading="lazy"></a><figcaption>${esc(cap)}</figcaption></figure>`).join("");
const legendGrad = MAP.legend.stops.map(([v, c]) => `${c} ${((v - MAP.legend.min) / (MAP.legend.max - MAP.legend.min)) * 100}%`).join(",");

const R = s("robust_PL"), RC = s("robust_CZ");
const top = Z.filter((p) => p.reliable === "так").sort((a, b) => a.up_mean_cm_yr - b.up_mean_cm_yr)[0];

const html = `<!doctype html><html lang="uk"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Острава–Карвіна: мульди і OSM</title>
<style>
:root{--bg:#f6f5f2;--card:#fff;--ink:#1d1d1b;--mut:#6b6a66;--line:#e2e0da;--acc:#b23a2b}
@media (prefers-color-scheme:dark){:root:not([data-theme=light]){--bg:#151514;--card:#1f1f1d;--ink:#ecebe7;--mut:#a09f99;--line:#34332f;--acc:#e0705f}}
:root[data-theme=dark]{--bg:#151514;--card:#1f1f1d;--ink:#ecebe7;--mut:#a09f99;--line:#34332f;--acc:#e0705f}
body{margin:0;background:var(--bg);color:var(--ink);font:15px/1.55 system-ui,Segoe UI,sans-serif}
main{max-width:1150px;margin:0 auto;padding:24px 16px}
h1{font-size:30px;margin:0 0 4px}h2{font-size:20px;margin:0 0 10px}h3{font-size:16px;margin:16px 0 6px}
section{background:var(--card);border:1px solid var(--line);border-radius:12px;padding:18px 20px;margin:16px 0}
.mut{color:var(--mut)}.small{font-size:12.5px}.wrap{overflow-x:auto}
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
<h1>Острава–Карвіна: кордон видобутку, видимий з космосу</h1>
<p class="mut">NISAR L2 GUNW (PROVISIONAL), 17.06–17.10.2026. Вертикаль з A158_F027 + D152_F063, підтвердження A057_F027. Мульди: поріг −15 см/рік, маска води OSM, поділ watershed. OSM станом на 14.11.2026. Усі швидкості — <b>вертикаль</b>, середні по мульді, см/рік (від'ємне = осідання), якщо не сказано інше.</p>

<section><h2>Проблема</h2>
<p>Остравсько-Карвінський (CZ) та Ястшембсько-Рибницький (PL) райони Верхньосілезького вугільного басейну — густонаселена агломерація з довгою історією підземного видобутку кам'яного вугілля. Над вибраними лавами порожнина поступово закривається, і земля над нею опускається, утворюючи <b>шахтні мульди</b>: локальне осідання від десятків сантиметрів до метрів на рік. Наслідки — тріщини в будинках, деформація доріг і залізниць, розриви комунікацій, підтоплення понижень і техногенна сейсмічність.</p>
<p>Ідея для атласу: <b>кордон, видимий з космосу</b>. На чеському боці видобуток згортається, і земля має бути спокійнішою. На польському шахти працюють, і земля має помітно провалюватися. Питання дослідження: чи видно цю різницю в перших даних NISAR, скільки будинків, доріг і залізниць стоїть на землі, що опускається, і що там: шахта, відвал чи вода.</p></section>

<section><h2>Ідея і метод</h2>
<div class="two"><div>
<ol class="steps">
<li><b>Дані.</b> Усі 12-денні інтерферограми NISAR L2 GUNW (L-діапазон, колекція <code>NISAR_L2_GUNW_PROVISIONAL_V1</code>) над bbox 18.15–18.80°E, 49.70–50.05°N: 57 гранул, з них використано 3 кадри (A057_F027, D152_F063, A158_F027).</li>
<li><b>Вирізка без завантаження кадрів.</b> HTTPS range-читання HDF5 (earthaccess + h5py), піксель 80 м.</li>
<li><b>Очистка фази.</b> Маска когерентності &lt; 0.3, край смуги, іоносферний екран, deramp, фільтр довгохвильового шуму (σ = 5 км), прив'язка до медіани сцени.</li>
<li><b>Швидкість LOS</b> — середнє по парах (зміщення / dt), см/рік.</li>
<li><b>Isolation Forest</b> на ознаках пікселя (швидкість, розкид пар, градієнт, когерентність, частка валідних пар) — бал аномальності.</li>
</ol></div><div>
<ol class="steps" start="6">
<li><b>Вертикаль і схід–захід</b> з пари висхідний + низхідний трек (A158 + D152). Вектори LOS з метаданих; NISAR дивиться <b>вліво</b>.</li>
<li><b>Мульди.</b> Вертикаль &lt; −15 см/рік мінус вода OSM → додаються білі центри (втрата когерентності) → злиття фрагментів ближче 400 м → поділ watershed для мульд з кількома центрами → лишаються підтверджені щонайменше 2 LOS-треками.</li>
<li><b>OSM-контекст.</b> Країна (admin_level=2), будівлі, дороги, залізниці, відвали, шахти, населені пункти, землекористування — для кожної мульди.</li>
<li><b>Перевірка.</b> Горизонталь у мульді має бути спрямована до центру (82% мульд). Окремо перевірено ставки, відвали та лінійні об'єкти.</li>
</ol></div></div>
<p class="small mut">Обмеження: дані PROVISIONAL і лише ~4 місяці, тож повільні міліметрові тренди ще не видно. Найшвидші центри (понад ~12 см LOS за 12 днів) близькі до межі однозначного розгортання фази L-діапазону. Назви шахт і статуси видобутку — мітки OSM, джерелом не перевірено. Валідації EGMS ще немає.</p></section>

<section><h2>Головна цифра</h2>
<div class="kpi">
<div><b>${n0(R.buildings)}</b>будівель стоять у ${R.zones} надійних польських мульдах (${R.area_km2} км²), які опускаються в середньому на ${Math.abs(R.up_mean_area_weighted_cm_yr)} см/рік (зважено на площу)</div>
<div><b>${RC.zones}</b>надійних мульд на чеському боці${RC.zones ? ` (${n0(RC.buildings)} будівель)` : " після маскування води"}</div>
<div><b>${R.roads_km} км</b>доріг і <b style="display:inline;font-size:inherit">${R.rail_km} км</b> залізниць у надійних польських мульдах</div>
</div>
<p class="small mut">Найактивніша надійна мульда: ${top.id} (${esc(top.locality)}, ${top.country}): ${top.up_mean_cm_yr} см/рік у середньому, p10 ${top.up_p10_cm_yr} см/рік. Пікові значення (~−2 м/рік LOS, ~−3 м/рік у згладженій вертикалі) — лише попередньо, до валідації.</p>
</section>

<section><h2>PL проти CZ</h2><div class="wrap"><table>
<tr><th>Набір</th><th>Країна</th><th>Мульд</th><th>Площа, км²</th><th>Будівель</th><th>Забудова, тис. м²</th><th>Дороги, км</th><th>Залізниці, км</th><th>Вертикаль (зваж.)</th><th>id</th></tr>${sumRows}</table></div>
<p class="small mut">Країна визначена за межею OSM admin_level=2 у центроїді мульди. Будівля рахується, якщо її центроїд у полігоні. Дороги: motorway…tertiary + residential/unclassified/service. Довжини рахувалися з кроком 10 м.</p></section>

<section><h2>Мапа</h2>
${mapFig(MAP.zoom, "Вертикальна швидкість, мульди, відвали і будівлі: кластер Ястшембе–Рибник")}
<div class="legend"><span>−40</span><span class="bar"></span><span>+40 см/рік</span>
<span><span class="sw" style="background:#8a6a3b88;border:1px solid #6b4f27"></span>відвал (OSM)</span>
<span><span class="sw" style="background:#1d1d1b;border-radius:50%;width:8px;height:8px"></span>будівля в мульді</span>
<span><span class="sw" style="background:#fff;border:2px solid #1d1d1b"></span>надійна мульда</span>
<span><span class="sw" style="background:#fff;border:1px dashed #5b5a56"></span>кандидат</span>
<span><span class="sw" style="background:#fff;border:1px solid #ccc"></span>немає даних (у центрах мульд — втрата когерентності через швидкий рух)</span>
<span><span class="sw" style="background:#55534e;height:3px"></span>дороги (OSM)</span>
<span><span class="sw" style="background:repeating-linear-gradient(90deg,#2b2a27 0 4px,#fff 4px 7px);height:3px"></span>залізниці</span>
<span><span class="sw" style="background:#6d3fc0;height:3px"></span>кордон CZ/PL</span>
<span><span class="sw" style="background:#111;width:9px;height:9px;outline:1px solid #fff"></span>шахта / промзона шахти (мітка OSM; малий квадрат — ствол)</span></div>
<h3>Увесь район</h3>
${mapFig(MAP.full, "Вертикальна швидкість по всьому району 18.15–18.80°E, 49.70–50.05°N")}
</section>

<section><h2>Мульди</h2><div class="wrap"><table>
<tr><th>id</th><th>group</th><th>Країна</th><th>Поблизу</th><th>км²</th><th>Верт. сер.</th><th>p10</th><th>Треків</th><th>Без даних</th><th>Відвал</th><th>Вода</th><th>Будівель</th><th>Дороги, км</th><th>Залізн., км</th><th>Найближча названа шахта (мітка OSM)</th></tr>${zoneRows}</table></div>
<p class="small mut">Напівпрозорі рядки — кандидати, не публікувати без валідації (EGMS). Назви шахт — мітки OSM, джерелом не перевірено.</p></section>

<section><h2>OSM: що всередині і навколо мульд</h2><div class="wrap"><table>
<tr><th>id</th><th>Населений пункт (центр, км) / найближче місце</th><th>Землекористування (частка площі)</th><th>Великі дороги</th><th>Залізниці</th><th>Забудова</th><th>Водотоки, км</th><th>Найближчий шахтний об'єкт</th><th>Об'єкти OSM всередині</th></tr>${osmRows}</table></div>
<p class="small mut">Землекористування — за сіткою 25 м усередині полігона (landuse / natural / man_made). «Забудова» — частка площі мульди під будівлями. Найближчий шахтний об'єкт включає стволи (man_made=mineshaft).</p>
<h3>Шахти з OSM і найближча надійна мульда</h3>
<div class="wrap"><table><tr><th>Назва (мітка OSM)</th><th>Тип</th><th>lat, lon</th><th>Найближча надійна мульда</th><th>Відстань до її центроїда, км</th><th>Її вертикаль, см/рік</th></tr>${mineRows}</table></div>
<p class="small mut">Відстань між центроїдами — грубо, без меж гірничих відводів. Статуси («w likwidacji», «nieczynna», «bývalý») — з назв в OSM, не перевірено. Межі відводів: PL — MIDAS PIG-PIB, CZ — геопортал ČGS (ще не зроблено).</p></section>

<section><h2>Графіки з ноутбука</h2>
<p class="small mut">Файли з <code>output/ostrava_karvina/</code>, згенеровані ноутбуком <code>scripts/ostrava/ostrava_karvina.ipynb</code>. Клік відкриває повний розмір. Часові ряди — у LOS, не у вертикалі.</p>
${figHtml}</section>

<section><h2>Історичні мульди (group)</h2>
<p>Watershed розділив колишні кластери на кілька центрів. Сума по групі — для порівняння з попередніми цифрами.</p>
<div class="wrap"><table><tr><th>group</th><th>Частини</th><th>Країна</th><th>Площа, км²</th><th>Вертикаль (зваж.)</th><th>Будівель</th><th>Дороги, км</th><th>Залізн., км</th></tr>${groupRows}</table></div></section>

<section><h2>Відвали: на відвалі чи поруч?</h2>
<div class="wrap"><table><tr><th>id</th><th>Частка відвалу</th><th>Вертикаль на відвалі</th><th>Вертикаль поза відвалом</th></tr>${heapRows}</table></div>
<p>Там, де відвал займає частину мульди, на ньому земля осідає повільніше, ніж поруч. Отже, найсильніше осідання — це не ущільнення відвалу, а земля над виробками.</p></section>

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
<p class="small mut">Скрипти: scripts/ostrava_osm.mjs, scripts/ostrava_map.mjs, scripts/build_ostrava_results.mjs. Дані OSM © учасники OpenStreetMap, ODbL. NISAR: NASA/ISRO, ASF DAAC, продукт PROVISIONAL.</p>
</main></body></html>`;
fs.writeFileSync(path.join(DIR, "results_ostrava.html"), html);
console.log("written results_ostrava.html");
