// Map for the Ostrava atlas place in plain Node (>= 18, no npm packages, no GDAL):
// vertical_cm_yr.tif (diverging red=subsidence / blue=uplift, clipped ±40 cm/yr) + spoil heaps (brown fill)
// + bowl outlines (reliable: thick dark, candidates: thin grey) + OSM buildings inside bowls (dots).
// Writes PNGs plus a JSON with label positions (pixels) so the HTML page can put SVG labels on top.
//   map_vertical_osm_full.png  - whole bbox, 2 px per 80 m cell
//   map_vertical_osm_zoom.png  - Jastrzębie–Rybnik cluster, 6 px per cell
// Usage: node scripts/ostrava/ostrava_map.mjs   (after ostrava_osm.mjs)
import fs from "node:fs";
import path from "node:path";
import { fileURLToPath } from "node:url";
import zlib from "node:zlib";

const ROOT = path.resolve(path.dirname(fileURLToPath(import.meta.url)), "..", "..");   // the repo root
const DIR = path.join(ROOT, "output", "ostrava_karvina");

// ---------- GeoTIFF (uncompressed float32) ----------
function readTif(f) {
  const b = fs.readFileSync(f), le = b.toString("ascii", 0, 2) === "II";
  const u16 = (o) => (le ? b.readUInt16LE(o) : b.readUInt16BE(o)), u32 = (o) => (le ? b.readUInt32LE(o) : b.readUInt32BE(o));
  const f64 = (o) => (le ? b.readDoubleLE(o) : b.readDoubleBE(o)), f32 = (o) => (le ? b.readFloatLE(o) : b.readFloatBE(o));
  const ifd = u32(4), n = u16(ifd), T = {};
  for (let i = 0; i < n; i++) {
    const o = ifd + 2 + i * 12, tag = u16(o), type = u16(o + 2), cnt = u32(o + 4);
    const size = { 3: 2, 4: 4, 12: 8, 2: 1 }[type] * cnt, at = size <= 4 ? o + 8 : u32(o + 8);
    const rd = (k) => (type === 3 ? u16(at + 2 * k) : type === 4 ? u32(at + 4 * k) : type === 12 ? f64(at + 8 * k) : 0);
    T[tag] = type === 2 ? b.toString("ascii", at, at + cnt - 1) : Array.from({ length: cnt }, (_, k) => rd(k));
  }
  const W = T[256][0], H = T[257][0];
  if ((T[259] || [1])[0] !== 1 || T[258][0] !== 32) throw new Error("only uncompressed float32 supported");
  const data = new Float32Array(W * H), rps = (T[278] || [H])[0];
  T[273].forEach((off, s) => { for (let k = 0; k < rps * W && s * rps * W + k < W * H; k++) data[s * rps * W + k] = f32(off + 4 * k); });
  const nodata = T[42113] ? parseFloat(T[42113]) : NaN;
  for (let i = 0; i < data.length; i++) if (data[i] === nodata) data[i] = NaN;
  const [sx, sy] = T[33550], tp = T[33922];
  return { W, H, data, x0: tp[3] - tp[0] * sx, y0: tp[4] + tp[1] * sy, sx, sy };
}
// ---------- UTM 33N ----------
const A = 6378137, F = 1 / 298.257223563, E2 = F * (2 - F), K0 = 0.9996;
function ll2utm(lon, lat, zone = 33) {
  const l0 = ((zone - 1) * 6 - 180 + 3) * Math.PI / 180, p = lat * Math.PI / 180, l = lon * Math.PI / 180;
  const ep2 = E2 / (1 - E2), N = A / Math.sqrt(1 - E2 * Math.sin(p) ** 2), T = Math.tan(p) ** 2, C = ep2 * Math.cos(p) ** 2, Aa = Math.cos(p) * (l - l0);
  const M = A * ((1 - E2 / 4 - 3 * E2 ** 2 / 64 - 5 * E2 ** 3 / 256) * p - (3 * E2 / 8 + 3 * E2 ** 2 / 32 + 45 * E2 ** 3 / 1024) * Math.sin(2 * p)
    + (15 * E2 ** 2 / 256 + 45 * E2 ** 3 / 1024) * Math.sin(4 * p) - (35 * E2 ** 3 / 3072) * Math.sin(6 * p));
  return [K0 * N * (Aa + (1 - T + C) * Aa ** 3 / 6 + (5 - 18 * T + T * T + 72 * C - 58 * ep2) * Aa ** 5 / 120) + 500000,
    K0 * (M + N * Math.tan(p) * (Aa * Aa / 2 + (5 - T + 9 * C + 4 * C * C) * Aa ** 4 / 24 + (61 - 58 * T + T * T + 600 * C - 330 * ep2) * Aa ** 6 / 720))];
}
// ---------- PNG writer ----------
const CRC = new Int32Array(256).map((_, n) => { let c = n; for (let k = 0; k < 8; k++) c = c & 1 ? 0xedb88320 ^ (c >>> 1) : c >>> 1; return c; });
const crc32 = (buf) => { let c = -1; for (const x of buf) c = CRC[(c ^ x) & 255] ^ (c >>> 8); return (c ^ -1) >>> 0; };
function chunk(type, data) { const len = Buffer.alloc(4); len.writeUInt32BE(data.length); const td = Buffer.concat([Buffer.from(type), data]); const c = Buffer.alloc(4); c.writeUInt32BE(crc32(td)); return Buffer.concat([len, td, c]); }
function writePng(f, w, h, rgba) {
  const raw = Buffer.alloc((w * 4 + 1) * h);
  for (let y = 0; y < h; y++) { raw[y * (w * 4 + 1)] = 0; rgba.copy(raw, y * (w * 4 + 1) + 1, y * w * 4, (y + 1) * w * 4); }
  const ihdr = Buffer.alloc(13); ihdr.writeUInt32BE(w, 0); ihdr.writeUInt32BE(h, 4); ihdr[8] = 8; ihdr[9] = 6;
  fs.writeFileSync(f, Buffer.concat([Buffer.from([137, 80, 78, 71, 13, 10, 26, 10]), chunk("IHDR", ihdr), chunk("IDAT", zlib.deflateSync(raw, { level: 9 })), chunk("IEND", Buffer.alloc(0))]));
}
// ---------- colours: red (subsidence) – neutral grey – blue (uplift) ----------
const hex = (h) => [1, 3, 5].map((i) => parseInt(h.slice(i, i + 2), 16));
const STOPS = [[-40, hex("#7a1414")], [-20, hex("#e34948")], [0, hex("#f0efec")], [20, hex("#2a78d6")], [40, hex("#123f7a")]];
function colour(v) {
  v = Math.max(-40, Math.min(40, v));
  for (let i = 1; i < STOPS.length; i++) if (v <= STOPS[i][0]) { const [a, ca] = STOPS[i - 1], [b, cb] = STOPS[i], t = (v - a) / (b - a); return ca.map((c, k) => Math.round(c + t * (cb[k] - c))); }
  return STOPS.at(-1)[1];
}

const R = readTif(path.join(DIR, "vertical_cm_yr.tif"));
const zones = JSON.parse(fs.readFileSync(path.join(DIR, "zones_vertical_osm.geojson"), "utf8"));
const heaps = JSON.parse(fs.readFileSync(path.join(DIR, "osm_heaps.geojson"), "utf8"));
const bld = JSON.parse(fs.readFileSync(path.join(DIR, "osm_buildings_in_bowls.json"), "utf8"));
const ringsOf = (g) => (g.type === "Polygon" ? [g.coordinates] : g.type === "MultiPolygon" ? g.coordinates : []);

// ---------- OSM context layers (cached Overpass answers) ----------
const CACHE = path.join(ROOT, "data", "ostrava", "osm_cache");
async function overpass(name, q) {
  const f = path.join(CACHE, name + ".json");
  if (fs.existsSync(f)) return JSON.parse(fs.readFileSync(f, "utf8"));
  for (let a = 0; a < 8; a++) {
    const r = await fetch("https://overpass-api.de/api/interpreter", { method: "POST", headers: { "Content-Type": "application/x-www-form-urlencoded", "User-Agent": "space-apps-research", Accept: "application/json" }, body: "data=" + encodeURIComponent(q) });
    const t = await r.text();
    if (r.ok && t.startsWith("{")) { fs.writeFileSync(f, t); return JSON.parse(t); }
    console.log(`  overpass ${name}: HTTP ${r.status}, retry`); await new Promise((s) => setTimeout(s, 15000 * (a + 1)));
  }
  throw new Error("overpass failed " + name);
}
const BB = "49.70,18.15,50.05,18.80", ZOOM = [18.50, 49.945, 18.73, 50.04], ZB = `${ZOOM[1]},${ZOOM[0]},${ZOOM[3]},${ZOOM[2]}`;
const lines = (j) => j.elements.filter((e) => e.geometry).map((e) => ({ t: e.tags || {}, pts: e.geometry.map((n) => [n.lon, n.lat]) }));
const L = {
  mainroads: lines(await overpass("map_mainroads", `[out:json][timeout:180];way["highway"~"^(motorway|trunk|primary)$"](${BB});out tags geom;`)),
  rail: lines(await overpass("map_rail", `[out:json][timeout:180];way["railway"~"^(rail|narrow_gauge)$"](${BB});out tags geom;`)),
  border: lines(await overpass("map_border", `[out:json][timeout:180];way["boundary"="administrative"]["admin_level"="2"](${BB});out tags geom;`)),
  zoomroads: lines(await overpass("map_zoomroads", `[out:json][timeout:180];way["highway"~"^(motorway|trunk|primary|secondary|tertiary)(_link)?$"](${ZB});out tags geom;`)),
};
const ctr = (e) => (e.center ? [e.center.lon, e.center.lat] : e.lat != null ? [e.lon, e.lat] : null);
const MINES = JSON.parse(fs.readFileSync(path.join(CACHE, "mines.json"), "utf8")).elements.map((e) => ({ p: ctr(e), t: e.tags })).filter((m) => m.p && m.t.man_made !== "salt_pond" && m.t.industrial !== "salt_pond");
const PLACES = JSON.parse(fs.readFileSync(path.join(CACHE, "places.json"), "utf8")).elements.map((e) => ({ p: ctr(e), name: e.tags.name, place: e.tags.place })).filter((m) => m.p);
console.log(`OSM layers: main roads ${L.mainroads.length}, rail ${L.rail.length}, border ${L.border.length}, zoom roads ${L.zoomroads.length}, mines ${MINES.length}`);

function render(name, win, scale) {
  const [lon0, lat0, lon1, lat1] = win || [];
  // window in raster cells
  let i0 = 0, j0 = 0, i1 = R.W - 1, j1 = R.H - 1;
  if (lon0 != null) {
    const pts = [[lon0, lat0], [lon1, lat0], [lon0, lat1], [lon1, lat1]].map(([a, b]) => ll2utm(a, b));
    i0 = Math.max(0, Math.floor((Math.min(...pts.map((p) => p[0])) - R.x0) / R.sx)); i1 = Math.min(R.W - 1, Math.ceil((Math.max(...pts.map((p) => p[0])) - R.x0) / R.sx));
    j0 = Math.max(0, Math.floor((R.y0 - Math.max(...pts.map((p) => p[1]))) / R.sy)); j1 = Math.min(R.H - 1, Math.ceil((R.y0 - Math.min(...pts.map((p) => p[1]))) / R.sy));
  }
  const w = (i1 - i0 + 1) * scale, h = (j1 - j0 + 1) * scale, img = Buffer.alloc(w * h * 4);
  const toPx = (lon, lat) => { const [x, y] = ll2utm(lon, lat); return [((x - R.x0) / R.sx - i0) * scale, ((R.y0 - y) / R.sy - j0) * scale]; };
  const put = (x, y, c, a = 255) => { x = Math.round(x); y = Math.round(y); if (x < 0 || y < 0 || x >= w || y >= h) return; const k = (y * w + x) * 4, t = a / 255;
    img[k] = img[k] * (1 - t) + c[0] * t; img[k + 1] = img[k + 1] * (1 - t) + c[1] * t; img[k + 2] = img[k + 2] * (1 - t) + c[2] * t; img[k + 3] = Math.max(img[k + 3], a); };
  // raster
  for (let y = 0; y < h; y++) for (let x = 0; x < w; x++) {
    const v = R.data[(j0 + Math.floor(y / scale)) * R.W + i0 + Math.floor(x / scale)], k = (y * w + x) * 4;
    if (isNaN(v)) { img[k] = img[k + 1] = img[k + 2] = 255; img[k + 3] = 255; } else { const c = colour(v); img[k] = c[0]; img[k + 1] = c[1]; img[k + 2] = c[2]; img[k + 3] = 255; }
  }
  // polygon fill (scanline, even-odd) and outlines
  function fillPoly(rings, c, a) {
    const pr = rings.map((r) => r.map(([lo, la]) => toPx(lo, la)));
    for (let y = 0; y < h; y++) {
      const xs = []; for (const r of pr) for (let i = 0, j = r.length - 1; i < r.length; j = i++) { const [xa, ya] = r[i], [xb, yb] = r[j]; if ((ya > y + 0.5) !== (yb > y + 0.5)) xs.push(xa + ((y + 0.5 - ya) / (yb - ya)) * (xb - xa)); }
      xs.sort((p, q) => p - q); for (let k = 0; k + 1 < xs.length; k += 2) for (let x = Math.ceil(xs[k]); x <= xs[k + 1]; x++) put(x, y, c, a);
    }
  }
  function line(rings, c, width, dash = 0) {
    for (const r of rings) for (let i = 1; i < r.length; i++) {
      const [xa, ya] = toPx(...r[i - 1]), [xb, yb] = toPx(...r[i]), n = Math.ceil(Math.hypot(xb - xa, yb - ya) * 2) || 1;
      for (let s = 0; s <= n; s++) { const x = xa + (s / n) * (xb - xa), y = ya + (s / n) * (yb - ya);
        if (dash && Math.floor((i * 7 + s / 2) / dash) % 2) continue;
        for (let dy = -width / 2; dy <= width / 2; dy += 0.5) for (let dx = -width / 2; dx <= width / 2; dx += 0.5) put(x + dx, y + dy, c); }
    }
  }
  for (const f of heaps.features) for (const rings of ringsOf(f.geometry)) { fillPoly(rings, hex("#8a6a3b"), 120); line(rings, hex("#6b4f27"), Math.max(1, scale / 3)); }
  // roads (grey), railways (dark, dashed), state border (violet)
  const zoomed = name === "zoom";
  for (const r of zoomed ? L.zoomroads : L.mainroads) { const major = /motorway|trunk|primary/.test(r.t.highway); line([r.pts], hex(major ? "#55534e" : "#8d8b85"), zoomed ? (major ? 2.5 : 1.5) : 1); }
  for (const r of L.rail) { line([r.pts], hex("#ffffff"), zoomed ? 3 : 2); line([r.pts], hex("#2b2a27"), zoomed ? 1.5 : 1, 4); }
  for (const r of L.border) { line([r.pts], hex("#ffffff"), zoomed ? 5 : 3.5); line([r.pts], hex("#6d3fc0"), zoomed ? 3 : 2); }
  const dot = Math.max(1, scale / 2.5);
  for (const [lo, la] of bld) { const [x, y] = toPx(lo, la); for (let dy = -dot; dy <= dot; dy += 0.5) for (let dx = -dot; dx <= dot; dx += 0.5) if (dx * dx + dy * dy <= dot * dot) put(x + dx, y + dy, hex("#1d1d1b"), 230); }
  const labels = [];
  for (const f of zones.features) {
    const p = f.properties, rel = p.reliable === "так";
    for (const rings of ringsOf(f.geometry)) {
      if (rel) { line(rings, hex("#ffffff"), Math.max(2, scale / 1.2)); line(rings, hex("#1d1d1b"), Math.max(1, scale / 2.5)); }
      else line(rings, hex("#5b5a56"), Math.max(1, scale / 3), 3);
    }
    const [x, y] = toPx(p.lon, p.lat);
    if (x >= 0 && y >= 0 && x < w && y < h) labels.push({ id: p.id, x: +(x / w * 100).toFixed(2), y: +(y / h * 100).toFixed(2), rel, v: p.up_mean_cm_yr, country: p.country });
  }
  // mines: shafts = small squares, named mine / mine plant = bigger square with white ring
  const inside = (x, y) => x >= 0 && y >= 0 && x < w && y < h;
  const pct2 = (x, y) => ({ x: +(x / w * 100).toFixed(2), y: +(y / h * 100).toFixed(2) });
  const mines = [];
  for (const m of MINES) {
    const [x, y] = toPx(...m.p); if (!inside(x, y)) continue;
    const big = m.t.man_made !== "mineshaft", s = (big ? (zoomed ? 6 : 4) : (zoomed ? 3 : 2));
    for (let dy = -s - 1; dy <= s + 1; dy++) for (let dx = -s - 1; dx <= s + 1; dx++) put(x + dx, y + dy, Math.abs(dx) > s || Math.abs(dy) > s ? hex("#ffffff") : hex(big ? "#111111" : "#3a3936"));
    if (big && m.t.name) {
      const full = m.t.name.replace(/Kopalnia Węgla Kamiennego/g, "KWK").replace(/[„”"]/g, "");
      const short = full.replace(/\s*\(.*?\)\s*/g, " ").replace(/\s*-\s*Ruch\s+/i, " / ").replace(/Borynia-Zofiówka-Bzie/, "BZB").trim();
      if (!mines.some((q) => q.name === short)) mines.push({ name: short, full, ...pct2(x, y) });
    }
  }
  const want = zoomed ? new Set(["town", "village"]) : new Set(["city", "town"]);
  const seenP = new Set();
  const places = PLACES.filter((p) => want.has(p.place) && !seenP.has(p.name) && seenP.add(p.name)).map((p) => { const [x, y] = toPx(...p.p); return inside(x, y) ? { name: p.name, place: p.place, ...pct2(x, y) } : null; }).filter(Boolean);
  const f = `map_vertical_osm_${name}.png`;
  writePng(path.join(DIR, f), w, h, img);
  console.log(`${f}: ${w}x${h}`);
  return { file: f, w, h, labels, mines, places };
}
const out = { full: render("full", null, 2), zoom: render("zoom", ZOOM, 6), legend: { min: -40, max: 40, stops: STOPS.map(([v, c]) => [v, "#" + c.map((x) => x.toString(16).padStart(2, "0")).join("")]) } };
fs.writeFileSync(path.join(DIR, "map_vertical_osm_labels.json"), JSON.stringify(out, null, 1));
