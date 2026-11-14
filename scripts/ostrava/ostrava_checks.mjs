// Checks for odd objects in the Ostrava NISAR vertical map (task 3):
//   1) M02: one bowl or two?  (connected components of the vertical field inside the polygon at several levels,
//      local minima of a 3x3-smoothed field)
//   2) the linear feature near E 771 km / N 5545 km (UTM33): which OSM lines run there
//   3) M01: profile of vertical rate + OSM objects (landuse composition is in zones_vertical_osm.csv)
// Reads vertical_cm_yr.tif (uncompressed float32 GeoTIFF, EPSG:32633) without GDAL.
// Usage: node scripts/ostrava/ostrava_checks.mjs
import fs from "node:fs";
import path from "node:path";
import { fileURLToPath } from "node:url";

const ROOT = path.resolve(path.dirname(fileURLToPath(import.meta.url)), "..", "..");   // the repo root
const DIR = path.join(ROOT, "output", "ostrava_karvina");
const CACHE = path.join(ROOT, "data", "ostrava", "osm_cache");

// ---------- minimal GeoTIFF reader ----------
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
  const data = new Float32Array(W * H), offs = T[273], rps = (T[278] || [H])[0];
  offs.forEach((off, s) => { for (let k = 0; k < rps * W && s * rps * W + k < W * H; k++) data[s * rps * W + k] = f32(off + 4 * k); });
  const nodata = T[42113] ? parseFloat(T[42113]) : NaN;
  for (let i = 0; i < data.length; i++) if (data[i] === nodata) data[i] = NaN;
  const [sx, sy] = T[33550], tp = T[33922];
  return { W, H, data, x0: tp[3] - tp[0] * sx, y0: tp[4] + tp[1] * sy, sx, sy };
}

// ---------- UTM (WGS84) ----------
const A = 6378137, F = 1 / 298.257223563, E2 = F * (2 - F), K0 = 0.9996;
function ll2utm(lon, lat, zone = 33) {
  const l0 = ((zone - 1) * 6 - 180 + 3) * Math.PI / 180, p = lat * Math.PI / 180, l = lon * Math.PI / 180;
  const ep2 = E2 / (1 - E2), N = A / Math.sqrt(1 - E2 * Math.sin(p) ** 2), T = Math.tan(p) ** 2, C = ep2 * Math.cos(p) ** 2, Aa = Math.cos(p) * (l - l0);
  const M = A * ((1 - E2 / 4 - 3 * E2 ** 2 / 64 - 5 * E2 ** 3 / 256) * p - (3 * E2 / 8 + 3 * E2 ** 2 / 32 + 45 * E2 ** 3 / 1024) * Math.sin(2 * p)
    + (15 * E2 ** 2 / 256 + 45 * E2 ** 3 / 1024) * Math.sin(4 * p) - (35 * E2 ** 3 / 3072) * Math.sin(6 * p));
  const x = K0 * N * (Aa + (1 - T + C) * Aa ** 3 / 6 + (5 - 18 * T + T * T + 72 * C - 58 * ep2) * Aa ** 5 / 120) + 500000;
  const y = K0 * (M + N * Math.tan(p) * (Aa * Aa / 2 + (5 - T + 9 * C + 4 * C * C) * Aa ** 4 / 24 + (61 - 58 * T + T * T + 600 * C - 330 * ep2) * Aa ** 6 / 720));
  return [x, y];
}
function utm2ll(x, y, zone = 33) { // Newton on ll2utm (accurate to mm, fine here)
  let lon = (zone - 1) * 6 - 180 + 3, lat = y / 111000;
  for (let i = 0; i < 20; i++) {
    const [x1, y1] = ll2utm(lon, lat, zone), [x2] = ll2utm(lon + 1e-5, lat, zone), [, y3] = ll2utm(lon, lat + 1e-5, zone);
    lon += (x - x1) / ((x2 - x1) / 1e-5); lat += (y - y1) / ((y3 - y1) / 1e-5);
  }
  return [lon, lat];
}

const R = readTif(path.join(DIR, "vertical_cm_yr.tif"));
const v = (i, j) => (i >= 0 && j >= 0 && i < R.W && j < R.H ? R.data[j * R.W + i] : NaN);
const pix = (x, y) => [Math.floor((x - R.x0) / R.sx), Math.floor((R.y0 - y) / R.sy)];
const xy = (i, j) => [R.x0 + (i + 0.5) * R.sx, R.y0 - (j + 0.5) * R.sy];
console.log(`vertical_cm_yr.tif ${R.W}x${R.H}, ${R.sx} m, UL ${R.x0.toFixed(0)} ${R.y0.toFixed(0)}`);

function inRing(p, r) { let c = false; for (let i = 0, j = r.length - 1; i < r.length; j = i++) if ((r[i][1] > p[1]) !== (r[j][1] > p[1]) && p[0] < ((r[j][0] - r[i][0]) * (p[1] - r[i][1])) / (r[j][1] - r[i][1]) + r[i][0]) c = !c; return c; }
function components(mask, w, h) { // 8-connected labels
  const lab = new Int32Array(w * h).fill(-1); let n = 0; const sizes = [];
  for (let s = 0; s < w * h; s++) {
    if (!mask[s] || lab[s] >= 0) continue; const st = [s]; lab[s] = n; let c = 0;
    while (st.length) { const q = st.pop(); c++; const qi = q % w, qj = (q - qi) / w;
      for (let dj = -1; dj <= 1; dj++) for (let di = -1; di <= 1; di++) { const a = qi + di, b2 = qj + dj, t = b2 * w + a;
        if (a >= 0 && b2 >= 0 && a < w && b2 < h && mask[t] && lab[t] < 0) { lab[t] = n; st.push(t); } } }
    sizes.push(c); n++;
  }
  return { lab, sizes };
}
const out = {};
const zones = JSON.parse(fs.readFileSync(path.join(DIR, "zones_vertical_clustered.geojson"), "utf8"));

// ---------- 1) M02 and the other big bowls ----------
for (const id of ["M02", "M04", "M03", "M06"]) {
  const ft = zones.features.find((f) => f.properties.id === id);
  const rings = (ft.geometry.type === "Polygon" ? [ft.geometry.coordinates] : ft.geometry.coordinates).map((r) => r[0].map(([lo, la]) => ll2utm(lo, la)));
  let [i0, j0, i1, j1] = [1e9, 1e9, -1e9, -1e9];
  for (const r of rings) for (const [x, y] of r) { const [i, j] = pix(x, y); i0 = Math.min(i0, i); i1 = Math.max(i1, i); j0 = Math.min(j0, j); j1 = Math.max(j1, j); }
  i0 -= 2; j0 -= 2; i1 += 2; j1 += 2; const w = i1 - i0 + 1, h = j1 - j0 + 1;
  const inside = new Uint8Array(w * h), sm = new Float32Array(w * h).fill(NaN);
  for (let j = 0; j < h; j++) for (let i = 0; i < w; i++) {
    const p = xy(i0 + i, j0 + j); inside[j * w + i] = rings.some((r) => inRing(p, r)) ? 1 : 0;
    let s = 0, c = 0; for (let dj = -1; dj <= 1; dj++) for (let di = -1; di <= 1; di++) { const z = v(i0 + i + di, j0 + j + dj); if (!isNaN(z)) { s += z; c++; } }
    if (c >= 3) sm[j * w + i] = s / c;
  }
  const polyParts = components(inside, w, h).sizes.length;
  const levels = {};
  for (const L of [-15, -30, -50, -80]) {
    const m = new Uint8Array(w * h); for (let k = 0; k < w * h; k++) m[k] = inside[k] && sm[k] < L ? 1 : 0;
    const { lab, sizes } = components(m, w, h);
    const big = sizes.map((s, k) => [s, k]).filter(([s]) => s >= 3).sort((a, b) => b[0] - a[0]);
    levels[L] = big.map(([s, k]) => { let sx = 0, sy = 0, mn = 0; for (let t = 0; t < w * h; t++) if (lab[t] === k) { const [x, y] = xy(i0 + (t % w), j0 + Math.floor(t / w)); sx += x; sy += y; mn = Math.min(mn, sm[t]); }
      const [lo, la] = utm2ll(sx / s, sy / s); return { px: s, km2: +(s * R.sx * R.sy / 1e6).toFixed(2), lon: +lo.toFixed(4), lat: +la.toFixed(4), min_sm: +mn.toFixed(0) }; });
  }
  // gap (no data) share and positions of minima of the smoothed field (window 5x5 = 400 m)
  const minima = [];
  for (let j = 2; j < h - 2; j++) for (let i = 2; i < w - 2; i++) {
    const k = j * w + i; if (!inside[k] || isNaN(sm[k]) || sm[k] > -40) continue; let isMin = true;
    for (let dj = -2; dj <= 2 && isMin; dj++) for (let di = -2; di <= 2; di++) { const t = (j + dj) * w + i + di; if (t !== k && sm[t] < sm[k]) { isMin = false; break; } }
    if (isMin) { const [lo, la] = utm2ll(...xy(i0 + i, j0 + j)); minima.push({ lon: +lo.toFixed(4), lat: +la.toFixed(4), v_cm_yr: +sm[k].toFixed(0) }); }
  }
  out[id] = { polygon_parts: polyParts, components_by_level: levels, minima_5x5: minima.sort((a, b) => a.v_cm_yr - b.v_cm_yr).slice(0, 6) };
  console.log(`\n${id}: polygon parts ${polyParts}`);
  for (const [L, cs] of Object.entries(levels)) console.log(`  < ${L} см/рік: ${cs.length} компонент(и) ` + cs.map((c) => `[${c.km2} км², ${c.lon},${c.lat}, min ${c.min_sm}]`).join(" "));
  console.log("  мінімуми:", minima.slice(0, 6).map((m) => `${m.lon},${m.lat} ${m.v_cm_yr}`).join(" | "));
}

// ---------- 2) linear feature near E 771 km, N 5545 km ----------
{
  const [cx, cy] = [771000, 5545000], [lo, la] = utm2ll(cx, cy);
  // vertical-rate profile along E-W and N-S through the point (±2 km) to see which way the line runs
  const prof = (dx, dy) => { const s = []; for (let d = -2000; d <= 2000; d += 80) { const [i, j] = pix(cx + d * dx, cy + d * dy); s.push(isNaN(v(i, j)) ? null : +v(i, j).toFixed(1)); } return s; };
  // local direction of the strongest gradient line: compare std along rows vs columns in a 3x3 km window
  console.log(`\nLinear feature: E771 N5545 = ${lo.toFixed(4)}E ${la.toFixed(4)}N`);
  console.log("  E-W profile (80 m step):", prof(1, 0).join(" "));
  console.log("  N-S profile (80 m step):", prof(0, 1).join(" "));
  const f = path.join(CACHE, "linear_771.json");
  let d;
  if (fs.existsSync(f)) d = JSON.parse(fs.readFileSync(f, "utf8"));
  else {
    const b = [la - 0.01, lo - 0.03, la + 0.01, lo + 0.03].map((x) => x.toFixed(5)).join(",");
    const q = `[out:json][timeout:90];(way["railway"](${b});way["highway"~"motorway|trunk|primary|secondary|tertiary"](${b});way["embankment"](${b});way["man_made"~"embankment|dyke"](${b});way["waterway"~"river|canal|stream"](${b}););out tags geom;`;
    for (let a = 0; a < 6 && !d; a++) {
      const r = await fetch("https://overpass-api.de/api/interpreter", { method: "POST", headers: { "Content-Type": "application/x-www-form-urlencoded", "User-Agent": "space-apps-research", Accept: "application/json" }, body: "data=" + encodeURIComponent(q) });
      const t = await r.text();
      if (r.ok && t.startsWith("{")) { d = JSON.parse(t); fs.writeFileSync(f, t); } else { console.log(`  overpass HTTP ${r.status}, retry`); await new Promise((s) => setTimeout(s, 15000 * (a + 1))); }
    }
    if (!d) d = { elements: [] };
  }
  const near = d.elements.map((e) => {
    let md = Infinity; for (const n of e.geometry || []) { const [x, y] = ll2utm(n.lon, n.lat); md = Math.min(md, Math.hypot(x - cx, y - cy)); }
    const t = e.tags; return { d: Math.round(md), what: t.railway ? `railway=${t.railway}${t.usage ? "/" + t.usage : ""}` : t.highway ? `highway=${t.highway}` : t.waterway ? `waterway=${t.waterway}` : t.power ? "power=line" : t.barrier ? `barrier=${t.barrier}` : "embankment", name: t.name || t.ref || "", embankment: t.embankment || "" };
  }).sort((a, b) => a.d - b.d);
  const seen = new Set(), list = near.filter((n) => { const k = n.what + n.name; if (seen.has(k)) return false; seen.add(k); return true; }).slice(0, 12);
  console.log("  OSM lines nearby:"); for (const n of list) console.log(`   ${n.d} m  ${n.what} ${n.name} ${n.embankment ? "(embankment)" : ""}`);
  out.linear_771 = { lon: +lo.toFixed(4), lat: +la.toFixed(4), osm_nearby: list, profile_EW: prof(1, 0), profile_NS: prof(0, 1) };
}

// ---------- 3) M01 profile ----------
{
  const [cx, cy] = ll2utm(18.35458, 49.78179); const s = [];
  for (let d = -1200; d <= 1200; d += 80) { const [i, j] = pix(cx + d, cy); s.push(isNaN(v(i, j)) ? null : +v(i, j).toFixed(0)); }
  console.log("\nM01 E-W vertical profile ±1.2 km:", s.join(" "));
  out.M01_profile_EW = s;
}
fs.writeFileSync(path.join(DIR, "checks_odd_objects.json"), JSON.stringify(out, null, 1));
