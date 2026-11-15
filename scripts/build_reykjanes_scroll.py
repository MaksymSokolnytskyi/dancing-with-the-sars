"""Reykjanes scroll plate: one self-contained HTML page (web/reykjanes-scroll.html) on the repo's map engine.

The engine (web/index.html, style.css, app.js, themes/scroll.json) is packed the same way scripts/pack_scene.py does it,
with a few patches so it speaks of uplift instead of sinking and opens an extra "math for the team" window.
The Reykjanes layers are built from products of the Reykjanes analysis pipeline (kept outside git, like all data),
found under --src DIR:

  breathing_data.npz      per-date vertical NISAR displacement (159A + 082D), rate, sigma (240 m, EPSG:3413)
  breathing_vectors.json  time-growing hachures, headline, GNSS series and the independent GNSS check
  map_layers.json         OSM roads / pipelines / lava barriers, source outlines
  ml/*.json               dike traces (NPE v3, blind v4), crater row, coherence-loss zones, FIRMS summary
  data/dem/Copernicus_DSM_COG_10_N63_00_W023_00_DEM.tif, data/firms/firms_reykjanes.csv, data/osm/reykjanes_layers.json

A separate "math behind the plate" window (docs/reykjanes-math.html, English, formulas typeset by KaTeX) is injected as an overlay.

Usage: uv run python scripts/build_reykjanes_scroll.py --src <folder with output/reykjanes and data/>
       -> web/reykjanes-scroll.html (open the file; needs internet only for MapLibre, deck.gl and the font)
"""
import argparse
import base64
import csv
import gzip
import json
import re
from datetime import date
from pathlib import Path

import numpy as np
import rasterio
from contourpy import contour_generator
from rasterio.warp import transform as warp
from rasterio.windows import from_bounds
from scipy.interpolate import RegularGridInterpolator
from scipy.ndimage import distance_transform_edt, gaussian_filter

ROOT = Path(__file__).resolve().parents[1]
WEB = ROOT / "web"
MATH = ROOT / "docs" / "reykjanes-math.html"   # the math window, English, KaTeX (Ukrainian original: docs/uk/reykjanes-math.html)
DST = WEB / "reykjanes-scroll.html"
SRC = ROOT                                       # set by --src: folder holding output/reykjanes/ and data/
OUT = SRC / "output" / "reykjanes"
ML = OUT / "ml"


def configure(src):
    global SRC, OUT, ML
    SRC = Path(src).resolve()
    OUT = SRC / "output" / "reykjanes"
    ML = OUT / "ml"
    need = [OUT / "breathing_data.npz", OUT / "breathing_vectors.json", OUT / "map_layers.json",
            ML / "coherence_fissures.json", ML / "s1_npe_v4_independent.json", ML / "firms_hotspots.json",
            SRC / "data/dem/Copernicus_DSM_COG_10_N63_00_W023_00_DEM.tif", SRC / "data/firms/firms_reykjanes.csv",
            SRC / "data/osm/reykjanes_layers.json"]
    miss = [str(f) for f in need if not f.exists()]
    if miss:
        raise SystemExit("missing inputs under --src:\n  " + "\n  ".join(miss))
SID = "reykjanes"
BBOX = [-22.78, 63.79, -22.03, 63.999]          # DEM tile ends at 64.0 N / 22.0 W
DLON, DLAT = 0.005, 0.0022                       # ~245 m at 63.9 N
SVARTSENGI = (-22.433, 63.878)
KX = 111.32 * np.cos(np.deg2rad(SVARTSENGI[1]))


def L(uk, en):
    return {"uk": uk, "en": en}


def local_to_ll(e, n):
    return [round(SVARTSENGI[0] + e / KX, 5), round(SVARTSENGI[1] + n / 111.32, 5)]


def smooth_nan(img, s):
    m = np.isfinite(img)
    num = gaussian_filter(np.where(m, img, 0.0), s)
    den = gaussian_filter(m.astype(float), s)
    return np.where(den > 0.4, num / np.maximum(den, 1e-6), np.nan)


def rate_sigma(y, yrs):
    """LSQ rate and AR(1)-inflated 1-sigma of a 1-D series (same rule as build_breathing_data.rate_sigma)."""
    ok = np.isfinite(y)
    if ok.sum() < 0.7 * len(y):
        return np.nan, np.nan
    T = yrs[ok] - yrs[ok].mean()
    v = y[ok] - y[ok].mean()
    r = (v * T).sum() / (T ** 2).sum()
    res = v - r * T
    s2 = (res ** 2).sum() / max(len(v) - 2, 1)
    rho = np.clip((res[:-1] * res[1:]).sum() / max((res ** 2).sum(), 1e-9), 0, 0.95)
    return r, np.sqrt(s2 / (T ** 2).sum()) * np.sqrt((1 + rho) / (1 - rho))


def fc(features):
    return {"type": "FeatureCollection", "features": features}


def line(coords, **props):
    return {"type": "Feature", "properties": props, "geometry": {"type": "LineString", "coordinates": coords}}


def point(lon, lat, **props):
    return {"type": "Feature", "properties": props, "geometry": {"type": "Point", "coordinates": [round(lon, 5), round(lat, 5)]}}


def clip_ok(lon, lat):
    return BBOX[0] <= lon <= BBOX[2] and BBOX[1] <= lat <= BBOX[3]


# ------------------------------------------------------------------ engine
PACKED_GETJSON = '''const getJSON = async (u) => {           // packed build: data lives in <script data-pack> blocks
  const el = document.querySelector(`script[data-pack="${u.split("?")[0]}"]`);
  if (!el) throw new Error("not packed: " + u);
  const bytes = Uint8Array.from(atob(el.textContent.trim()), (c) => c.charCodeAt(0));
  return JSON.parse(await new Response(new Blob([bytes]).stream().pipeThrough(new DecompressionStream("gzip"))).text());
};'''


def engine():
    """The repo's map app, packed like scripts/pack_scene.py: (html head up to the app, app js, html tail, scroll theme)."""
    js = (WEB / "app.js").read_text(encoding="utf-8")
    old = re.search(r"^const getJSON = .*$", js, re.M)
    assert old, "getJSON line not found in app.js"
    js = js.replace(old.group(0), PACKED_GETJSON, 1)
    js = re.sub(r"^const THEMES = \[.*?\];", 'const THEMES = ["scroll"];', js, count=1, flags=re.M)
    # map labels may be {en, uk}, like every other label (the page opens in English; ?lang=uk for Ukrainian)
    js = must(js, "s: letter(f.properties.name, f.properties.rank ?? 1)", "s: letter(label(f.properties.name), f.properties.rank ?? 1)")
    css = (WEB / "style.css").read_text(encoding="utf-8")
    css += "\n/* packed build: one place, one look */\n#panel section:has(#scene) { display: none; }\n"
    html = (WEB / "index.html").read_text(encoding="utf-8")
    html, n = re.subn(r'<link rel="stylesheet" href="style\.css[^"]*">', lambda m: f"<style>\n{css}</style>", html)
    assert n == 1, "style.css link not found in web/index.html"
    i = html.index('<script type="module" src="app.js')
    j = html.index("</script>", i) + len("</script>")
    theme = json.loads((WEB / "themes" / "scroll.json").read_text(encoding="utf-8"))
    return html[:i], js, html[j:], theme


def must(s, old, new, count=1):
    n = s.count(old)
    if n < 1:
        raise SystemExit(f"engine patch failed, not found: {old[:80]!r}")
    return s.replace(old, new) if count == 0 else s.replace(old, new, count)


# wording fixes from Oleh's review of 15 Nov 2026, applied to the built story (and by repack_reykjanes.py to a packed one)
TEXT_FIX = [
    ("rose by about 25 cm", "rose by about 23 cm"), ("піднявся приблизно на 25 см", "піднявся приблизно на 23 см"),
    ("As the year plays the wedges grow; each ring of strokes is 2 more cm.",
     "As the year plays the wedges grow. Strokes stand in rings, each ring is 2 more cm of uplift: the tighter the rings, the steeper "
     "the dome bends. Thin lines and bands of tone show how far the ground has risen since the first date: one every 10 cm."),
    ("Поки програється рік, клини ростуть; кожне кільце штрихів — ще 2 см.",
     "Поки програється рік, клини ростуть. Штрихи стоять кільцями, кожне кільце — ще 2 см підняття: що тісніші кільця, то крутіше "
     "вигинається купол. Тонкі лінії й смуги тону показують, наскільки земля піднялась від першої дати: одна на кожні 10 см."),
    ("a deeper chamber in our model needs a larger dike.",
     "a deeper chamber in our model needs a larger dike. For this event our chamber lies at about 6 km, deeper than the 4.7 km of "
     "page III: it is another event fitted to other data, and a layered crust puts the source about 1 km deeper."),
    ("глибша камера в нашій моделі потребує більшої дайки.",
     "глибша камера в нашій моделі потребує більшої дайки. Для цієї події наша камера лежить на глибині близько 6 км, глибше за 4,7 км "
     "зі сторінки III: це інша подія й інші дані, а шарувата кора опускає джерело приблизно на 1 км."),
]


def fix_text(st):
    txt = json.dumps(st, ensure_ascii=False)
    for a, b in TEXT_FIX:
        if b not in txt:
            txt = txt.replace(a, b)
    st = json.loads(txt)
    for p in st.get("pages", []):                      # the layer a page talks about blinks when the camera arrives
        if p.get("section") == "exhale":
            p.setdefault("focus", ["dike_blind"])
        if p.get("section") == "why" and "source" in p.get("layers", []):
            p.setdefault("focus", ["source"])
    return st


def patch_app(app, with_math):
    uk = {
        'stroke: "кожне кільце штрихів = {b} см опускання"': 'stroke: "кожне кільце штрихів = {b} см підняття"',
        'sinks: "опускається ≈{v} см/рік"': 'sinks: "піднімається ≈{v} см/рік"',
        'frame_sub: "опускання землі за радаром NISAR · {a}–{b}"': 'frame_sub: "підняття землі за радаром NISAR · {a}–{b}"',
        'lg_wedge: "клин: товстий кінець вище, тонкий вказує, куди земля опускається більше"':
            'lg_wedge: "клин: товстий кінець там, де земля піднялась більше; тонкий вказує, куди підняття згасає"',
        'lg_dense: "густі штрихи: земля прогинається найкрутіше"': 'lg_dense: "густі штрихи: купол вигинається найкрутіше"',
        'lg_contour: "лінії глибини кожні {s} см від першої дати"': 'lg_contour: "лінії підняття кожні {s} см від першої дати"',
        'headline: "Де ми певні, що тоне найшвидше: <b>{r} ± {s} см/рік</b> (щонайменше {b}). 1% площі опускається швидше за {p} см/рік"':
            'headline: "Де ми певні, що піднімається найшвидше: <b>{r} ± {s} см/рік</b> (щонайменше +{b}). 1% надійної площі піднімається швидше за {p} см/рік"',
        'lg_up: "синій: земля піднімається"': 'lg_up: "синій: земля опускається"',
        'score_note: "клин = опускання району за 12 днів · угорі найшвидші, внизу нерухомі"':
            'score_note: "клин = підняття місця за 12 днів · угорі найшвидші, внизу нерухомі"',
        'score_sound_t: "кожен дзвін — ще 2 см опускання: що швидше район сідає, то частіше дзвонить; нерухомі мовчать"':
            'score_sound_t: "кожен дзвін — ще 2 см підняття: що швидше місце здіймається, то частіше дзвонить; нерухомі мовчать"',
        'story_next: "що ще може NISAR"': 'story_next: "що ще може NISAR", story_exhale: "видих", story_math: "для команди"',
        'story_open: "Читати таблицю ▸",': 'story_open: "Читати таблицю ▸", math_open: "∑ Математика таблиці", '
                                           'test_build: "ТЕСТОВА ВЕРСІЯ",',
    }
    en = {
        'stroke: "each ring of strokes = {b} cm of sinking"': 'stroke: "each ring of strokes = {b} cm of uplift"',
        'sinks: "sinking ≈{v} cm/yr"': 'sinks: "rising ≈{v} cm/yr"',
        'frame_sub: "ground sinking seen by the NISAR radar · {a}–{b}"': 'frame_sub: "ground uplift seen by the NISAR radar · {a}–{b}"',
        'lg_wedge: "wedge: thick end is higher, thin end points to where the ground sinks more"':
            'lg_wedge: "wedge: thick end where the ground rose most, thin end points to where the uplift fades"',
        'headline: "Where we are sure it sinks fastest: <b>{r} ± {s} cm/yr</b> (at least {b}). 1% of the area sinks faster than {p} cm/yr"':
            'headline: "Where we are sure it rises fastest: <b>{r} ± {s} cm/yr</b> (at least +{b}). 1% of the reliable area rises faster than {p} cm/yr"',
        'lg_contour: "depth lines every {s} cm since the first date"': 'lg_contour: "uplift lines every {s} cm since the first date"',
        'score_note: "wedge = how far the district went down in 12 days · fastest on top, still at the bottom"':
            'score_note: "wedge = how far the place rose in 12 days · fastest on top, still at the bottom"',
        'score_sound_t: "each chime = 2 more cm of sinking: the faster a district sinks, the more often it chimes; still ground stays silent"':
            'score_sound_t: "each chime = 2 more cm of uplift: the faster a place rises, the more often it chimes; still ground stays silent"',
        'story_next: "what else NISAR can do"': 'story_next: "what else NISAR can do", story_exhale: "exhale", story_math: "for the team"',
        'story_open: "Read the plate ▸",': 'story_open: "Read the plate ▸", math_open: "∑ The math behind the plate", '
                                          'test_build: "TEST BUILD",',
    }
    for k, v in {**uk, **en}.items():
        app = must(app, k, v)
    # uplift is the main colour; sinking strokes take color_down
    app = must(app, "getColor: (p) => [...(p.up ? c.color_up || c.color : c.color),",
               "getColor: (p) => [...(p.up ? c.color : c.color_down || c.color_up || c.color),")
    app = must(app, 'if (v != null && v < 0) parts.push(T("sinks")', 'if (v != null && v > 0) parts.push(T("sinks")')
    # block diagram: the extreme by magnitude, not the most negative value
    app = must(app, "for (const v of fin) if (v != null && v < deepest) deepest = v;",
               "for (const v of fin) if (v != null && Math.abs(v) > Math.abs(deepest)) deepest = v;")
    app = must(app, "pxPerCm = 4 * gap / Math.max(1, -deepest / 10)", "pxPerCm = 4 * gap / Math.max(1, Math.abs(deepest) / 10)")
    # partitura: a rising site is "moving" (sense "up" in rhythm.json)
    app = must(app, "function scoreSteps(d) {", "const SG = (d) => (d.sense === \"up\" ? -1 : 1);\nfunction scoreSteps(d) {")
    app = must(app, "return v == null || a == null ? null : (v - a) * 12 / Math.max(1, days[j + 1] - days[j]);",
               "return v == null || a == null ? null : SG(d) * (v - a) * 12 / Math.max(1, days[j + 1] - days[j]);")
    app = must(app, "const lvl = Math.floor(-seriesAt(d.rows[i], days, state.day) / BELL_CM)",
               "const lvl = Math.floor(-SG(d) * seriesAt(d.rows[i], days, state.day) / BELL_CM)")
    # the math window: a link in every story card; arrows and Esc belong to the window while it is open
    if with_math:
        app = must(app, '`<div class="st-foot"><a href="#" id="st-explore">${T("story_explore")}</a></div>`',
                   '`<div class="st-foot"><a href="#" class="mathbtn">${T("math_open")}</a> · <a href="#" id="st-explore">${T("story_explore")}</a></div>`')
    # the packed page carries the whole atlas index (for the onward link) but opens on its own plate, not on the cover
    app = must(app, 'if (idx.atlas && (!q.get("scene") || q.get("cover") === "1")) showCover(true);',
               'if (idx.atlas && q.get("cover") === "1") showCover(true);')
    app = must(app, 'const scene = q.get("scene") || idx.scenes.find((s) => !s.href).id,',
               f'const scene = q.get("scene") || "{SID}",')
    app = must(app, "if (!state.storyOn || /INPUT|SELECT|TEXTAREA/.test(",
               "if (document.body.classList.contains(\"math-open\") || !state.storyOn || /INPUT|SELECT|TEXTAREA/.test(")
    return app


def patch_head(head, with_math):
    head = re.sub(r"<title>.*?</title>", "<title>Reykjanes · The ground breathes</title>", head, count=1)
    extra = ""                                           # the jury build (--no-math) carries no test badge
    if with_math:
        extra += '\n      <div class="testbuild" data-i18n="test_build">TEST BUILD</div>'
        extra +='\n      <button class="mathbtn" id="math-open" data-i18n="math_open">∑ The math behind the plate</button>'
    # the atlas link of index.html leads back to the atlas cover from this separate page
    head = must(head, '<a href="#" id="to-cover" class="small">☰ Edelweiss</a>', '<a href="index.html?cover=1" id="back-atlas" style="display:inline-block;margin-bottom:6px;color:var(--muted);text-decoration:none;letter-spacing:.04em;font-size:12px">☰ Edelweiss</a>')
    head = must(head, '<button id="story-open" class="hidden"></button>', '<button id="story-open" class="hidden"></button>' + extra)
    head = must(head, "</head>", "<style>.testbuild{display:inline-block;margin:8px 6px 0 0;padding:2px 7px;border:1px dashed #8c2e1e;"
                "color:#8c2e1e;font-size:11px;letter-spacing:.08em}</style>\n</head>")
    return head


def theme_for(th):
    th = json.loads(json.dumps(th))
    if "tint" in th:                                     # layer tints: uplift warm like the strokes, sinking slate
        th["tint"].update({"up": [[230, 198, 172], [132, 46, 24]], "down": [[192, 208, 216], [52, 84, 112]]})
    th["velocity"]["color"] = [150, 44, 26]              # uplift: magma red
    th["velocity"]["color_down"] = [52, 84, 110]         # sinking: the old blue
    th["ramps"]["sigma_r"] = [[0, [0, 0, 0, 0]], [1.0, [40, 70, 140, 15]], [2.0, [40, 70, 140, 80]], [3.5, [30, 50, 120, 160]]]
    th["lines"].update({
        "coast": {"color": [64, 96, 100], "alpha": 230, "width": 1.1, "legend": False},
        "sea": {"color": [64, 96, 100], "alpha": 170, "width": 0.8, "fade": True, "legend": False},
        "road": {"color": [58, 42, 28], "alpha": 200, "width": 1.6, "casing": 0.7},
        "pipe": {"color": [58, 42, 28], "alpha": 210, "width": 1.1, "dash": [5, 3]},
        "power": {"color": [58, 42, 28], "alpha": 150, "width": 0.8, "dash": [1, 3]},
        "wall": {"color": [120, 60, 20], "alpha": 240, "width": 2.4, "zebra": [3, 2], "casing": 0.5},
        "source": {"color": [150, 44, 26], "alpha": 235, "width": 1.6, "dash": [6, 3], "label_size": 13, "label_color": [150, 44, 26],
                   "italic": True, "dot": 3.2},
        "dike": {"color": [110, 30, 18], "alpha": 255, "width": 3.2, "label_size": 13, "label_color": [110, 30, 18], "italic": True},
        "dike_blind": {"color": [110, 30, 18], "alpha": 200, "width": 1.6, "dash": [3, 3]},
        "row": {"color": [58, 42, 28], "alpha": 230, "width": 1.4, "dash": [8, 3, 2, 3]},
        "fire": {"color": [190, 70, 20], "alpha": 255, "width": 0, "dot": 1.7},
        "lossband": {"color": [120, 60, 20], "alpha": 80, "width": 5.0},
    })
    th["lines"]["place"]["caps_rank"] = 3
    th["lines"]["peak"]["collide"] = True
    th["layers_on"] = []
    return th


# ------------------------------------------------------------------ the regular lon/lat grid
class Grid:
    def __init__(self):
        self.lon = np.arange(BBOX[0], BBOX[2] + 1e-9, DLON)
        self.lat = np.arange(BBOX[1], BBOX[3] + 1e-9, DLAT)
        self.LON, self.LAT = np.meshgrid(self.lon, self.lat)
        d = np.load(OUT / "breathing_data.npz")
        self.dates = [str(x) for x in d["dates"]]
        def tr(lon, lat):                        # EPSG:4326 -> EPSG:3413 (the CRS of breathing_data.npz)
            x, y = warp("EPSG:4326", "EPSG:3413", np.ravel(lon).tolist(), np.ravel(lat).tolist())
            return np.reshape(x, np.shape(lon)), np.reshape(y, np.shape(lat))
        X, Y = tr(d["lon"].astype(float), d["lat"].astype(float))
        xs, ys = X.mean(0), Y.mean(1)
        o = np.argsort(ys)
        self._xs, self._ys, self._o = xs, ys[o], o
        self.TX, self.TY = tr(self.LON, self.LAT)
        self.src = d

    def resample(self, a):
        f = RegularGridInterpolator((self._ys, self._xs), a[self._o], bounds_error=False, fill_value=np.nan)
        return f(np.dstack([self.TY, self.TX]))


def build_fields(G):
    d = G.src
    up = d["up"]
    up_s = np.array([smooth_nan(u, 1.5) for u in up])
    frames_s = np.array([G.resample(u) for u in up_s])                 # cm, smoothed (strokes, contours, cut)
    frames = np.array([G.resample(u) for u in up])                     # cm, raw (cells)
    rate = G.resample(d["rate"].astype(float))
    sigma = G.resample(d["sigma"].astype(float))
    sig_s = G.resample(smooth_nan(d["sigma"].astype(float), 1.5))
    return frames_s, frames, rate, sigma, sig_s


def years(dates):
    d0 = date.fromisoformat(dates[0])
    return np.array([(date.fromisoformat(x) - d0).days / 365.25 for x in dates])


# ------------------------------------------------------------------ layers
def velocity_layer(dates):
    V = json.loads((OUT / "breathing_vectors.json").read_text(encoding="utf-8"))
    days = np.array([(date.fromisoformat(x) - date.fromisoformat(dates[0])).days for x in dates], float)
    slopes = np.array([h["slope"] for h in V["hachures"]])
    s95 = np.percentile(slopes, 95)
    feats = []
    for h in V["hachures"]:
        a, b = np.array(h["a"]), np.array(h["b"])
        if not (clip_ok(*a) and clip_ok(*b)):
            continue
        lev = h["level"]
        sg = 1.0 if lev >= 0 else -1.0
        s = np.array([np.nan if v is None else v for v in h["t"]], float) * sg
        s = np.where(np.isfinite(s), s, -np.inf)
        s = np.maximum.accumulate(s)
        n = 7
        coords, ts = [], []
        for k in range(n):
            f = k / (n - 1)
            coords.append([round(float(a[0] + f * (b[0] - a[0])), 6), round(float(a[1] + f * (b[1] - a[1])), 6)])
            thr = f * abs(lev)
            j = int(np.argmax(s >= thr)) if np.any(s >= thr) else len(s) - 1
            if j == 0:
                t = 0.0
            elif s[j] >= thr and np.isfinite(s[j - 1]):
                w = (thr - s[j - 1]) / max(s[j] - s[j - 1], 1e-6)
                t = days[j - 1] + np.clip(w, 0, 1) * (days[j] - days[j - 1])
            else:
                t = days[j]
            ts.append(round(float(max(t, ts[-1] if ts else 0.0)), 1))
        feats.append({"type": "Feature", "properties": {
            "rate_cm_yr": h["rate"], "sigma_cm_yr": h["sigma"], "snr": h["snr"], "segments": 1, "trust": h["trust"],
            "dark": round(float(np.clip(h["slope"] / s95, 0.15, 1.0)), 3), "sense": "up" if lev >= 0 else "down", "t": ts},
            "geometry": {"type": "LineString", "coordinates": coords}})
    trust = {k: sum(f["properties"]["trust"] == k for f in feats) for k in ("strong", "weak", "below")}
    return fc(feats), V, trust


def contours_layer(G, frames_s, step=5):
    feats = []
    for di, fr in enumerate(frames_s):
        cg = contour_generator(x=G.lon, y=G.lat, z=np.ma.masked_invalid(fr))
        hi, lo = np.nanmax(fr), np.nanmin(fr)
        levels = [l for l in range(step, int(hi) + 1, step)] + [-l for l in range(step, int(-lo) + 1, step)]
        for lev in levels:
            for ln in cg.lines(lev):
                if len(ln) < 8:
                    continue
                ln = ln[::2] if len(ln) > 60 else ln
                feats.append(line([[round(float(x), 5), round(float(y), 5)] for x, y in ln],
                                  level_cm=lev, major=lev % 10 == 0, d=di))
    return fc(feats)


def surface_layer(G, frames_s, sig_s, V, dates):
    fr = [[None if not np.isfinite(v) else int(round(v * 10)) for v in f.ravel()] for f in frames_s]
    sg = [None if not np.isfinite(v) else round(float(v) * 10, 1) for v in sig_s.ravel()]
    gn = []
    for s, g in V["gnss"].items():
        t = np.array([(date.fromisoformat(x) - date.fromisoformat(dates[0])).days for x in g["t"]], float)
        q = np.array([(date.fromisoformat(x) - date.fromisoformat(dates[0])).days for x in dates], float)
        cm = np.interp(q, t, np.array(g["up"], float), left=np.nan, right=np.nan)
        cm = cm - (cm[np.isfinite(cm)][0] if np.isfinite(cm).any() else 0)
        gn.append({"station": s, "lon": g["lon"], "lat": g["lat"], "cm": [None if not np.isfinite(v) else round(float(v), 1) for v in cm]})
    return {"kind": "surface", "lon0": BBOX[0], "lat0": BBOX[1], "dlon": DLON, "dlat": DLAT, "nx": len(G.lon), "ny": len(G.lat),
            "unit": "mm", "dates": dates, "frames": fr, "sigma_mm_yr": sg, "years": float(years(dates)[-1]),
            "section": [[-22.445, 63.825], [-22.395, 63.965]], "gnss": gn}


def grid_layer(G, frames, rate, sigma, dates):
    cells = []
    ny, nx = rate.shape
    for i in range(0, ny - 1, 2):
        for j in range(0, nx - 1, 2):
            r = rate[i:i + 2, j:j + 2]
            if np.isfinite(r).sum() < 3:
                continue
            s = sigma[i:i + 2, j:j + 2]
            ser = np.nanmean(frames[:, i:i + 2, j:j + 2].reshape(len(dates), -1), 1)
            cells.append([round(float(G.lon[j] + DLON / 2), 5), round(float(G.lat[i] + DLAT / 2), 5),
                          round(float(np.nanmean(r)), 2), round(float(np.nanmean(s)), 2),
                          [None if not np.isfinite(v) else round(float(v), 1) for v in ser]])
    return {"cell_km": 0.49, "columns": ["lon", "lat", "rate_cm_yr", "sigma_cm_yr", "cum_cm"], "dates": dates,
            "value": "rate_cm_yr", "series": "cum_cm", "series_unit": "cm",
            "modes": [{"key": "sigma_cm_yr", "label": L("похибка σ", "error σ"), "ramp": "sigma_r"}], "cells": cells}


RHYTHM = [  # name, lon, lat. Hafnir (west edge, -13 cm/yr) left out: an edge/plane artefact, not a place we can vouch for
    (L("Свартсенгі (GPS SENG)", "Svartsengi (GPS SENG)"), -22.4338, 63.8823),
    (L("Блакитна лагуна", "Blue Lagoon"), -22.449, 63.880),
    (L("Ельдвьорп", "Eldvörp"), -22.535, 63.865),
    (L("Торбйорн", "Þorbjörn"), -22.426, 63.857),
    (L("Сундгнукур", "Sundhnúkur"), -22.390, 63.880),
    (L("Гріндавік", "Grindavík"), -22.432, 63.842),
    (L("Фаградальсфьядль", "Fagradalsfjall"), -22.270, 63.900),
    (L("Воґар", "Vogar"), -22.387, 63.983),
    (L("Ньярдвік (NYLA, нуль)", "Njarðvík (NYLA, zero)"), -22.7381, 63.9735),
]


def rhythm_layer(G, frames, dates, headline):
    yrs = years(dates)
    km_x, km_y = (G.LON - 0) * KX, G.LAT * 111.32
    rows = []
    for name, lon, lat in RHYTHM:
        if lon is None:
            lon, lat = headline["lon"], headline["lat"]
        r_km = 1.0
        m = np.hypot(km_x - lon * KX, km_y - lat * 111.32) <= r_km
        if m.sum() < 4:
            continue
        ser = np.array([np.nanmedian(f[m]) if np.isfinite(f[m]).sum() >= 3 else np.nan for f in frames])
        if np.isfinite(ser).sum() < 0.7 * len(ser):
            continue
        ser = ser - ser[0]
        ref = "NYLA" in name["uk"]
        rt, sg = rate_sigma(ser, yrs)
        if ref:
            ser = np.zeros_like(ser)
            rt = 0.0
        rows.append({"name": name, "lon": round(lon, 4), "lat": round(lat, 4), "r_km": r_km, "cells": int(m.sum()), "ref": ref,
                     "series_cm": [None if not np.isfinite(v) else round(float(v), 2) for v in ser],
                     "rate_cm_yr": round(float(rt), 2), "sigma_cm_yr": round(float(sg), 2),
                     "still": bool(ref or abs(rt) < 2 * sg)})
    rows.sort(key=lambda r: (r["ref"], -r["rate_cm_yr"]))
    return {"dates": dates, "unit": "cm", "r_km": 1.0, "sense": "up", "rows": rows,
            "rule": "row = median of the 245 m cells within r_km; still = |rate| < 2 sigma (AR(1)-inflated); the NYLA row is the zero "
                    "point of the map, so it is flat by construction"}


def dem_layers():
    with rasterio.open(SRC / "data" / "dem" / "Copernicus_DSM_COG_10_N63_00_W023_00_DEM.tif") as r:
        w = from_bounds(BBOX[0], BBOX[1], BBOX[2], BBOX[3], r.transform)
        z = r.read(1, window=w).astype(float)
        tf = r.window_transform(w)
    # ~90 m working grid
    k = 3
    h, wd = (z.shape[0] // k) * k, (z.shape[1] // k) * k
    zz = z[:h, :wd].reshape(h // k, k, wd // k, k).mean((1, 3))
    lon = tf.c + (np.arange(zz.shape[1]) + 0.5) * tf.a * k
    lat = tf.f + (np.arange(zz.shape[0]) + 0.5) * tf.e * k
    lat, zz = lat[::-1], zz[::-1]
    land = zz > 0.5
    dx = abs(tf.a) * k * KX * 1000
    dy = abs(tf.e) * k * 111320
    # coast and ghost lines out at sea
    cg = contour_generator(x=lon, y=lat, z=gaussian_filter(land.astype(float), 0.8))
    coast = [line([[round(float(a), 5), round(float(b), 5)] for a, b in ln[::2]]) for ln in cg.lines(0.5) if len(ln) > 20]
    dist = distance_transform_edt(~land, sampling=(dy, dx))
    cg2 = contour_generator(x=lon, y=lat, z=gaussian_filter(dist, 1.0))
    ghost = []
    for lev, kk in ((260, 0.15), (650, 0.45), (1200, 0.75)):
        ghost += [line([[round(float(a), 5), round(float(b), 5)] for a, b in ln[::3]], k=kk) for ln in cg2.lines(lev) if len(ln) > 30]
    # relief hachures: short strokes down the slope, thicker where steeper
    zs = gaussian_filter(np.where(land, zz, 0), 1.0)
    gy, gx = np.gradient(zs, dy, dx)
    rng = np.random.default_rng(7)
    sp = 230.0
    feats = []
    for yy in np.arange(0, lat.size * dy, sp):
        for xx in np.arange(0, lon.size * dx, sp):
            px, py = xx + rng.uniform(-0.35, 0.35) * sp, yy + rng.uniform(-0.35, 0.35) * sp
            j, i = int(px / dx), int(py / dy)
            if not (0 <= i < lat.size and 0 <= j < lon.size) or not land[i, j]:
                continue
            g = np.array([gx[i, j], gy[i, j]])
            s = np.hypot(*g)
            if s < 0.035:
                continue
            u = -g / s
            ln_m = sp * (0.45 + 0.4 * min(1.0, s / 0.25))
            p0 = np.array([px, py]) - u * ln_m / 2
            p1 = np.array([px, py]) + u * ln_m / 2
            c = [[round(float(lon[0] + p[0] / dx * (lon[1] - lon[0])), 5), round(float(lat[0] + p[1] / dy * (lat[1] - lat[0])), 5)]
                 for p in (p0, (p0 + p1) / 2, p1)]
            feats.append(line(c, dark=round(float(min(1.0, s / 0.3)), 3)))
    # elevation for the readout (~180 m)
    e = zz[::2, ::2]
    elev = {"kind": "elevation", "unit": "m", "lon0": float(lon[0]), "lat0": float(lat[0]),
            "dlon": float(lon[2] - lon[0]), "dlat": float(lat[2] - lat[0]), "nx": int(e.shape[1]), "ny": int(e.shape[0]),
            "z": [None if v <= 0.5 else int(round(v)) for v in e.ravel()]}
    return fc(feats), fc(coast), fc(ghost), elev


def osm_layers():
    M = json.loads((OUT / "map_layers.json").read_text(encoding="utf-8"))
    out = {}
    for key, kind in (("keyroads", "road"), ("pipelines", "pipe"), ("power", "power"), ("barriers", "wall")):
        out[kind] = fc([line([[round(x, 5), round(y, 5)] for x, y in f["geometry"]["coordinates"]])
                        for f in M["layers"][key]["features"]
                        if any(clip_ok(x, y) for x, y in f["geometry"]["coordinates"])])
    O = json.loads((SRC / "data" / "osm" / "reykjanes_layers.json").read_text(encoding="utf-8"))
    rank = {"town": 3, "village": 2, "suburb": 2, "hamlet": 1}
    pl = [point(p["lon"], p["lat"], name=p["name"], rank=rank[p["place"]]) for p in O["places"]
          if p["place"] in rank and clip_ok(p["lon"], p["lat"]) and p["name"]]
    for nm, lon, lat, rk in (("Svartsengi", -22.433, 63.878, 2), ("Bláa lónið", -22.449, 63.880, 1), ("Sundhnúkur", -22.390, 63.880, 2),
                             ("Fagradalsfjall", -22.270, 63.900, 2), ("Eldvörp", -22.535, 63.865, 1)):
        pl.append(point(lon, lat, name=nm, rank=rk))
    out["places"] = fc(pl)
    pk = []
    for p in O["peaks"]:
        try:
            ele = int(float(p["ele"]))
        except (TypeError, ValueError):
            continue
        if ele >= 120 and p["name"] and clip_ok(p["lon"], p["lat"]):
            pk.append((ele, point(p["lon"], p["lat"], name=f"▲ {p['name']} {ele}", rank=1)))
    pk.sort(key=lambda x: -x[0])
    out["peaks"] = fc([p for _, p in pk[:14]])
    s = M["sources"]
    ell = [[round(b, 5), round(a, 5)] for a, b in s["ellipsoid"]]
    ell.append(ell[0])
    out["source"] = fc([line(ell), point(s["mogi"][1], s["mogi"][0], dot=1),
                        point(s["mogi"][1] + 0.01, s["mogi"][0] - 0.0105, name=L("джерело: 4.7 км (шари) · 4.1 км (Моґі)", "source: 4.7 km (layers) · 4.1 km (Mogi)"), rank=1)])
    return out


def dike_layers():
    C = json.loads((ML / "coherence_fissures.json").read_text(encoding="utf-8"))
    I4 = json.loads((ML / "s1_npe_v4_independent.json").read_text(encoding="utf-8"))
    v3 = C["dike_trace_npe_v3"]
    v4 = I4["trace_v4"]
    a, b = local_to_ll(v3[0], v3[1]), local_to_ll(v3[2], v3[3])
    mid = [(a[0] + b[0]) / 2, (a[1] + b[1]) / 2]
    dike = fc([line([a, b]), point(mid[0] + 0.035, mid[1] - 0.004, name=L("дайка 10.11.2023 · NPE v3", "dike 10 Nov 2023 · NPE v3"), rank=1)])
    blind = fc([line([local_to_ll(v4[0], v4[1]), local_to_ll(v4[2], v4[3])])])
    row = fc([line([[round(x, 5), round(y, 5)] for x, y in C["crater_row"]])])
    bands = []
    for f in C["events"]["real_20231110"]["fissures"]:
        th = np.deg2rad(f["strike_deg"])
        h = f["length_km"] / 2
        e0, n0 = f["e"], f["n"]
        bands.append(line([local_to_ll(e0 - h * np.sin(th), n0 - h * np.cos(th)), local_to_ll(e0 + h * np.sin(th), n0 + h * np.cos(th))]))
    F = json.loads((ML / "firms_hotspots.json").read_text(encoding="utf-8"))
    wins = [(date.fromisoformat(e["start"]), date.fromisoformat(e["end"])) for e in F["eruptions"]]
    # one dot per 400 m cell that burned repeatedly (>= 8 VIIRS detections in the eruption windows): isolated single
    # detections (cooling lava elsewhere, the power plant, false alarms) would otherwise outnumber the lava itself
    cnt = {}
    with open(SRC / "data" / "firms" / "firms_reykjanes.csv", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            if r["instrument"] != "VIIRS" or r["confidence"] == "l":
                continue
            dd = date.fromisoformat(r["acq_date"])
            if not any(s <= dd <= e for s, e in wins):
                continue
            lon, lat = float(r["longitude"]), float(r["latitude"])
            if not clip_ok(lon, lat):
                continue
            key = (round(lon * KX / 0.4), round(lat * 111.32 / 0.4))
            c = cnt.setdefault(key, [0, 0.0, 0.0])
            c[0] += 1; c[1] += lon; c[2] += lat
    pts = [point(c[1] / c[0], c[2] / c[0], dot=1) for c in cnt.values() if c[0] >= 8]
    return dike, blind, row, fc(bands), fc(pts), F


# ------------------------------------------------------------------ story
STATIONS_ORDER = ("SENG", "AUSV", "GEVK", "VOGC")


def gnss_chart(V, lang):
    """Small inline SVG for a story card: NISAR vs GPS vertical rate (±1σ) at the four independent stations."""
    g = V["gnss"]
    s0 = g["SENG"]
    ink, red = "#3a2a1c", "#962c1a"
    W, top, rh, x0, x1 = 300, 18, 24, 92, 270
    lo, hi = -5.0, 25.0
    X = lambda v: x0 + (x1 - x0) * (v - lo) / (hi - lo)
    out = [f'<line x1="{X(0):.1f}" x2="{X(0):.1f}" y1="{top - 6}" y2="{top + 4 * rh - 6}" stroke="{ink}" stroke-width=".4" stroke-dasharray="2,2"/>']
    for v in (0, 10, 20):
        out.append(f'<text x="{X(v):.1f}" y="{top + 4 * rh + 6}" font-size="9" text-anchor="middle" fill="{ink}">{v}</text>')
    unit = "см/рік" if lang == "uk" else "cm/yr"
    out.append(f'<text x="{x1 + 26}" y="{top + 4 * rh + 6}" font-size="9" text-anchor="end" fill="{ink}">{unit}</text>')
    km = "км" if lang == "uk" else "km"
    for i, s in enumerate(STATIONS_ORDER):
        v = g[s]["val"]
        d = ((g[s]["lon"] - s0["lon"]) * 48.97) ** 2 + ((g[s]["lat"] - s0["lat"]) * 111.32) ** 2
        y = top + i * rh
        ok = abs(v["z"]) < 2
        out.append(f'<text x="0" y="{y + 3}" font-size="10.5" fill="{ink}"{" font-weight=\"bold\"" if ok else ""}>{s}</text>'
                   f'<text x="34" y="{y + 3}" font-size="9" fill="{ink}" opacity=".7">{d ** .5:.1f} {km}</text>')
        for val, sg, col, dy, fill in ((v["nisar"], v["nisar_sigma"], red, -4, red), (v["gnss"], v["gnss_sigma"], ink, 4, "#e8dcbf")):
            out.append(f'<line x1="{X(val - sg):.1f}" x2="{X(val + sg):.1f}" y1="{y + dy}" y2="{y + dy}" stroke="{col}" stroke-width="1.4"/>'
                       f'<circle cx="{X(val):.1f}" cy="{y + dy}" r="3" fill="{fill}" stroke="{col}" stroke-width="1.2"/>')
        out.append(f'<text x="{W}" y="{y + 3}" font-size="9.5" text-anchor="end" fill="{red if not ok else ink}">z {v["z"]:+.1f}</text>')
    leg = ("● радар NISAR  ○ GPS  (±1σ); жирним — збіг |z| &lt; 2" if lang == "uk"
           else "● NISAR radar  ○ GPS  (±1σ); bold — agreement |z| &lt; 2")
    H = top + 4 * rh + 22
    out.append(f'<text x="0" y="{H - 2}" font-size="9" fill="{ink}" opacity=".8">{leg}</text>')
    return f'<svg viewBox="0 0 {W} {H}" width="100%" style="display:block;margin:4px 0">{"".join(out)}</svg>'


def story(V, dates, with_math=True):
    H = V["headline"]
    cam_all = {"bbox": [-22.74, 63.80, -22.10, 63.995]}
    cam_dome = {"bbox": [-22.62, 63.815, -22.30, 63.93]}
    P = []
    P.append({"section": "intro", "title": L("Земля, що дихає", "The ground that breathes"),
              "text": [L("Під Свартсенгі на півострові Рейк'янес збирається магма, і земля над нею здіймається, наче груди на вдиху. "
                         "Кожні 12 днів радар NISAR вимірював цей рух з космосу, вдень і вночі, крізь полярну темряву й хмари.",
                         "Under Svartsengi on the Reykjanes peninsula magma gathers, and the ground above it rises like a chest breathing in. "
                         "Every 12 days the NISAR radar measured it from space, day and night, through polar dark and clouds."),
                       L("За рік центр купола піднявся приблизно на 25 см. Потім земля видихає: магма тріщиною-дайкою йде до поверхні, "
                         "і купол різко сідає.",
                         "In one year the centre of the dome rose by about 25 cm. Then the ground breathes out: magma escapes through a dike, "
                         "and the dome drops.")],
              "facts": [{"value": f"+{H['rate']}", "unit": L("см/рік", "cm/yr"), "label": L(f"найшвидше підняття, ±{H['sigma']}", f"fastest uplift, ±{H['sigma']}")},
                        {"value": f"{H['area_km2']:.0f}", "unit": L("км²", "km²"), "label": L("надійно піднімається (≥ 3σ, > 5 см/рік)", "rises reliably (≥ 3σ, > 5 cm/yr)")},
                        {"value": "4", "unit": L("мм", "mm"), "label": L("різниця з GPS SENG за рік (225 проти 229 мм); збігається 1 з 4 станцій", "difference from GPS SENG over the year (225 vs 229 mm); 1 of 4 stations agrees")}],
              "camera": cam_all, "layers": ["terrain", "velocity", "coast", "sea", "places", "peaks", "infra"], "day": "end"})
    P.append({"section": "what", "title": L("Дивіться, як здіймається купол", "Watch the dome rise"),
              "text": [L("Кожен червоний клин стоїть на схилі купола: товстий кінець там, де земля піднялась більше. Поки програється рік, "
                         "клини ростуть; кожне кільце штрихів — ще 2 см.",
                         "Each red wedge stands on the flank of the dome, thick where the ground rose more. As the year plays the wedges grow; "
                         "each ring of strokes is 2 more cm."),
                       L("Пунктир — там, де радар менш певний. Сині клини — земля, що трохи опускається (переважно краї й шум).",
                         "Dashes mark where the radar is less sure. Blue wedges are ground that sinks a little (mostly edges and noise).")],
              "camera": cam_dome, "layers": ["terrain", "velocity", "contours", "coast", "sea", "places", "peaks"], "day": "start",
              "play": True, "legend": True})
    P.append({"section": "why", "title": L("Чому тут: магма на глибині ~5 км", "Why here: magma ~5 km down"),
              "text": [L("Форма купола видає джерело. Найпростіша модель (точкове джерело Моґі) кладе його на 4.1 км; пружна кора з шарами, "
                         "яка м'якша згори, опускає його до 4.7 км. Видовжена форма — резервуар ~6 км, витягнутий на 58° від півночі, "
                         "уздовж рифтових тріщин півострова.",
                         "The shape of the dome gives away its source. The simplest model (a Mogi point source) puts it at 4.1 km; a layered "
                         "crust, softer on top, moves it to 4.7 km. The elongated fit is a ~6 km reservoir at 58° from north, along the rift."),
                       L("Щороку туди надходить 15–18 млн м³ магми — приблизно 20 олімпійських басейнів щодня.",
                         "15–18 million m³ of magma arrive every year — about 20 Olympic pools a day.")],
              "facts": [{"value": "4.7", "unit": L("км", "km"), "label": L("глибина джерела, шарувата кора (Моґі: 4.1)", "source depth, layered crust (Mogi: 4.1)")},
                        {"value": "15–18", "unit": L("млн м³/рік", "Mm³/yr"), "label": L("приплив магми", "magma inflow")},
                        {"value": "~6", "unit": L("км", "km"), "label": L("довжина резервуару, вісь 58°", "reservoir length, axis 58°")}],
              "camera": cam_dome, "layers": ["terrain", "velocity", "source", "coast", "sea", "places", "infra"], "day": "end",
              "marks": [{"lon": H["lon"], "lat": H["lat"], "text": L(f"центр купола: +{H['rate']} см/рік", f"dome centre: +{H['rate']} cm/yr"), "dx": -40, "dy": -30},
                        {"lon": -22.4338, "lat": 63.8823, "text": L("Свартсенгі, GPS SENG", "Svartsengi, GPS SENG"), "dx": 44, "dy": -26}]})
    P.append({"section": "why", "title": L("Розріз: купол над камерою", "The cut: a dome over the chamber"),
              "text": [L("Розріз іде з півдня на північ: Гріндавік → Свартсенгі → Воґар. Будиночки нахилені з величезним перебільшенням: "
                         "так видно, де купол вигинається найкрутіше і тягне труби геотермальної станції та дороги.",
                         "The cut runs south to north: Grindavík → Svartsengi → Vogar. The houses tilt with huge exaggeration, showing where "
                         "the dome bends most and pulls on the power-plant pipes and the roads."),
                       L("Кружечок на профілі — GPS-станція SENG на тій самій осі: вона незалежно підтверджує висоту купола. Захисні вали "
                         "(коричневі) зведено довкола Свартсенгі й Гріндавіка після вивержень 2023–2024 рр.",
                         "The circle on the profile is GPS station SENG on the same axis: it independently confirms the height of the dome. "
                         "The lava barriers (brown) were built around Svartsengi and Grindavík after the 2023–2024 eruptions.")],
              "camera": {"bbox": [-22.60, 63.82, -22.25, 63.975]}, "layers": ["terrain", "velocity", "surface", "walls", "infra", "coast", "places"],
              "day": "end", "cut": "cut", "line": [[-22.445, 63.825], [-22.395, 63.965]]})
    P.append({"section": "listen", "title": L("У кожного місця свій ритм вдиху", "Each place breathes at its own pace"),
              "text": [L("Кожен ряд унизу — одне місце, кожен клин — 12 днів. Натисніть ▶ з увімкненим звуком: місце дзвонить щоразу, "
                         "коли піднімається ще на 2 см.",
                         "Each row below is one place, each wedge 12 days. Press ▶ with sound on: a place chimes every time it rises 2 more cm."),
                       L("Блакитна лагуна й Свартсенгі над камерою дзвонять близько 10 разів за рік, Гріндавік — двічі, Фаградальсфьядль мовчить. "
                         "Ньярдвік — нуль карти, він мовчить за побудовою.",
                         "The Blue Lagoon and Svartsengi above the chamber chime about 10 times a year, Grindavík twice, Fagradalsfjall stays silent. "
                         "Njarðvík is the zero of the map and silent by construction.")],
              "camera": cam_dome, "layers": ["terrain", "velocity", "rhythm", "coast", "sea", "places"], "day": "start"})
    P.append({"section": "sure", "title": L("Наскільки ми певні?", "How sure are we?"),
              "text": [L("Чотири GPS-станції, які не брали участі в оцінці. <b>Збігається одна з чотирьох:</b> SENG у Свартсенгі "
                         "(радар 20.1 ± 1.7, GPS 21.0 ± 1.2 см/рік, z = −0.4; за рік 225 проти 229 мм). Три інші розходяться на 4–7 см/рік: "
                         "AUSV і GEVK стоять лише за 4–5 км від SENG, біля Гріндавіка, де GPS ще показує 7–8 см/рік, а радар — 0–4.",
                         "Four GPS stations that took no part in the estimate. <b>One of the four agrees:</b> SENG at Svartsengi "
                         "(radar 20.1 ± 1.7, GPS 21.0 ± 1.2 cm/yr, z = −0.4; over the year 225 vs 229 mm). The other three differ by 4–7 cm/yr: "
                         "AUSV and GEVK are only 4–5 km from SENG, by Grindavík, where GPS still shows 7–8 cm/yr and the radar 0–4."),
                       L(gnss_chart(V, "uk"), gnss_chart(V, "en")),
                       L("Отже, надійне число — центр купола. Похибка σ у клітинках (сині, 1–1.7 см/рік) замала поруч із розбіжностями 4–7 см/рік: "
                         "вона рахує лише розкид ряду в часі, а не систематику (атмосфера на великих відстанях, точка відліку). Причину ми ще не знайшли.",
                         "So the firm number is the dome centre. The per-cell error σ (blue, 1–1.7 cm/yr) is too small next to 4–7 cm/yr misfits: "
                         "it counts only the scatter of the series in time, not systematic errors (long-range atmosphere, the reference point). "
                         "We have not found the cause yet.")],
              "camera": cam_dome, "layers": ["terrain", "velocity", "grid", "coast", "sea", "places"], "day": "end",
              "modes": {"grid": "sigma_cm_yr"},
              "marks": [{"lon": -22.4338, "lat": 63.8823, "text": L("SENG: радар +20.1 · GPS +21.0 см/рік", "SENG: radar +20.1 · GPS +21.0 cm/yr"), "dx": 46, "dy": -34},
                        {"lon": -22.4169, "lat": 63.8436, "text": L("AUSV: радар +3.5 · GPS +7.6", "AUSV: radar +3.5 · GPS +7.6"), "dx": 40, "dy": 30},
                        {"lon": -22.4716, "lat": 63.8354, "text": L("GEVK: радар −0.3 · GPS +7.0", "GEVK: radar −0.3 · GPS +7.0"), "dx": -40, "dy": 34}]})
    P.append({"section": "exhale", "title": L("Видих: дайка 10 листопада 2023", "Breathing out: the dike of 10 November 2023"),
              "text": [L("Так виглядає видих. За кілька годин магма з-під Свартсенгі пішла тріщиною на ~15 км під Гріндавік, а купол "
                         "просів. NISAR тоді ще не літав, тож ми перевірили метод на двох радарних парах Sentinel-1.",
                         "This is what breathing out looks like. Within hours magma from under Svartsengi ran ~15 km along a crack under "
                         "Grindavík, and the dome dropped. NISAR was not flying yet, so we tested the method on two Sentinel-1 radar pairs."),
                       L("Нейромережа (NPE), навчена на симуляціях, за 7 секунд дала об'єм дайки <b>150 [131–170] млн м³</b> і просідання "
                         "Свартсенгі <b>−77 [−92…−63] млн м³</b> (68-відсоткові інтервали). Стаття в Science (~20 станцій GPS, чотири супутники): "
                         "130–139 і −76…−82. Інтервал просідання повністю покриває значення статті. Інтервал об'єму покриває майже весь її "
                         "діапазон (131–139 із 130–139), але наша медіана вища за нього на 8–15 %: глибша камера в нашій моделі потребує "
                         "більшої дайки.",
                         "A neural network (NPE) trained on simulations gave, in 7 seconds, a dike volume of <b>150 [131–170] million m³</b> and a "
                         "Svartsengi drop of <b>−77 [−92…−63] million m³</b> (68 % intervals). The Science paper (~20 GPS stations, four satellites): "
                         "130–139 and −76…−82. Our drop interval covers the paper's values entirely. Our volume interval covers almost all of the "
                         "paper's range (131–139 of 130–139), but our median is 8–15 % above it: a deeper chamber in our model needs a larger dike."),
                       L("Помаранчеві точки — місця, що горіли під час вивержень 2023–2025 рр. (NASA FIRMS): лава вийшла на кратерному ряду над знайденою дайкою і розтеклася.",
                         "Orange dots are places that burned in the 2023–2025 eruptions (NASA FIRMS): lava came out on the crater row above the dike we found and spread.")],
              "facts": [{"value": "150", "unit": L("млн м³", "Mm³"), "label": L("об'єм дайки [131–170]; стаття 130–139", "dike volume [131–170]; paper 130–139")},
                        {"value": "−77", "unit": L("млн м³", "Mm³"), "label": L("Свартсенгі [−92…−63]; стаття −76…−82", "Svartsengi [−92…−63]; paper −76…−82")},
                        {"value": "7", "unit": L("с", "s"), "label": L("відповідь NPE (MCMC: 7 хв)", "NPE answer (MCMC: 7 min)")}],
              "camera": {"bbox": [-22.58, 63.80, -22.28, 63.93]},
              "layers": ["terrain", "dike", "dike_blind", "crater_row", "lossband", "firms", "source", "walls", "coast", "sea", "places"],
              "day": "end",
              "marks": [{"lon": -22.432, "lat": 63.842, "text": L("Гріндавік: евакуйовано 10.11.2023", "Grindavík: evacuated 10.11.2023"), "dx": -44, "dy": 26}]})
    P.append({"section": "next", "title": L("Що ще NISAR може розповісти тут", "What else NISAR can tell here"),
              "text": [L("Нова пара NISAR приходить кожні 12 днів. Мережа, навчена один раз, оцінює нову дайку за 0.35 с замість 49 с для MCMC "
                         "(1000 подій: ~6 хв проти ~14 год). На тихих парах вона каже «дайки немає» з імовірністю 0.86. "
                         "<b>Поки це надійно лише для великих дайок:</b> на малих подіях 2024 р. мережа невпевнена або занижує "
                         "(29.05.2024: 7 [3–14] млн м³ проти 13–19 у літературі).",
                         "A new NISAR pair arrives every 12 days. A network trained once sizes a new dike in 0.35 s instead of 49 s for MCMC "
                         "(1000 events: ~6 min vs ~14 h). On quiet pairs it says 'no dike' with probability 0.86. "
                         "<b>So far this is reliable for large dikes only:</b> on the small 2024 events the network is unsure or too low "
                         "(29 May 2024: 7 [3–14] million m³ vs 13–19 in the literature)."),
                       L("Сліпа мережа, яка не знала, де шукати, знаходить дайку з точністю ~2 км (випадкова точка — 5.1 км). "
                         "А втрата когерентності в тих самих радарних парах показує, де земля рветься: для дайки 2023 р. біля кратерного ряду "
                         "вона трапляється в 4.9 раза частіше, ніж випадково.",
                         "A blind network that did not know where to look finds the dike within ~2 km (a random point: 5.1 km). "
                         "And loss of coherence in the same radar pairs shows where the ground tears: for the 2023 dike it is 4.9 times more "
                         "frequent near the crater row than by chance.")],
              "facts": [{"value": "0.35", "unit": L("с", "s"), "label": L("нова велика дайка після нової пари", "a new large dike after a new pair")},
                        {"value": "~2", "unit": L("км", "km"), "label": L("сліпий пошук дайки", "blind dike search")},
                        {"value": "0.86", "label": L("P(без дайки) на тихій парі", "P(no dike) on a quiet pair")}],
              "camera": cam_all, "layers": ["terrain", "velocity", "dike_blind", "firms", "coast", "sea", "places", "infra"], "day": "end"})
    if with_math:
        P.append({"section": "math", "title": L("Математика для команди", "Math for the team"),
                  "text": [L("Тестова версія: внутрішній розділ для команди (англійською). Як ми рахували, чому NPE і MCMC іноді розходяться, "
                             "і як наш NPE v3 порівнюється зі статтею Sigmundsson et al. 2024 (Science) по кожному числу.",
                             "Test build: an internal section for the team (in English). How we computed, why NPE and MCMC sometimes disagree, "
                             "and how our NPE v3 compares with Sigmundsson et al. 2024 (Science), number by number."),
                           L('<button class="mathbtn bigbtn">∑ Відкрити математичне вікно</button>', '<button class="mathbtn bigbtn">∑ Open the math window</button>'),
                           L("Клавіша M відкриває вікно з будь-якої сторінки; Esc закриває.", "Key M opens the window from any page; Esc closes it.")],
                  "camera": {"bbox": [-22.58, 63.80, -22.28, 63.93]}, "layers": ["terrain", "dike", "source", "coast", "sea", "places"], "day": "end"})
    return fix_text({"title": L("Рейк'янес", "Reykjanes"), "pages": P,
            "note": "Pages of the Reykjanes plate, written for the Mexico City scroll engine; numbers from output/reykjanes/results.html."})


# ------------------------------------------------------------------ assemble
def pack(obj):
    raw = obj if isinstance(obj, bytes) else json.dumps(obj, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    return base64.b64encode(gzip.compress(raw, 9, mtime=0)).decode()


def main():
    ap = argparse.ArgumentParser(description="Build web/reykjanes-scroll.html")
    ap.add_argument("--src", default=str(SRC), help="folder with output/reykjanes/ and data/{dem,firms,osm}/")
    ap.add_argument("--out", default=str(DST))
    ap.add_argument("--no-math", dest="math", action="store_false", help="leave out the internal math window (jury build)")
    a = ap.parse_args()
    configure(a.src)
    dst = Path(a.out)
    head, app, tail, theme = engine()
    G = Grid()
    dates = G.dates
    print("dates", dates[0], "->", dates[-1], len(dates))
    frames_s, frames, rate, sigma, sig_s = build_fields(G)
    vel, V, trust = velocity_layer(dates)
    H = V["headline"]
    print("strokes", len(vel["features"]), trust)
    cont = contours_layer(G, frames_s)
    surf = surface_layer(G, frames_s, sig_s, V, dates)
    grid = grid_layer(G, frames, rate, sigma, dates)
    rhy = rhythm_layer(G, frames, dates, H)
    print("rhythm", [(r["name"]["uk"], r["rate_cm_yr"], r["sigma_cm_yr"], r["still"]) for r in rhy["rows"]])
    terrain, coast, ghost, elev = dem_layers()
    print("terrain strokes", len(terrain["features"]), "coast", len(coast["features"]))
    osm = osm_layers()
    dike, blind, row, bands, firms, F = dike_layers()
    print("firms dots", len(firms["features"]))

    infra_g = {"group": "infra", "group_label": L("Інфраструктура (OpenStreetMap)", "Infrastructure (OpenStreetMap)"),
               "attribution": "© OpenStreetMap contributors"}
    files = {}

    def add(lid, kind, obj, fname, **kw):
        files[fname] = obj
        return {"id": lid, "kind": kind, "file": fname, **kw}

    layers = [
        add("terrain", "hachures-static", terrain, "terrain.geojson", label=L("Рельєф (Copernicus GLO-30)", "Relief (Copernicus GLO-30)"), props=["dark"]),
        add("velocity", "hachures-time", vel, "velocity.geojson", label=L("Підняття ґрунту", "Ground uplift"),
            props=["rate_cm_yr", "sigma_cm_yr", "snr", "segments", "trust", "dark", "sense", "t"]),
        add("contours", "contours", cont, "contours.geojson", label=L("Лінії підняття", "Uplift lines"), props=["level_cm", "major", "d"],
            unit="cm", step_cm=5, visible=False),
        add("surface", "surface", surf, "surface.json", label=L("Розріз", "Cross-section"), visible=False),
        add("grid", "cells", grid, "grid.json", label=L("Клітинки довіри", "Trust cells"),
            props=["lon", "lat", "rate_cm_yr", "sigma_cm_yr", "cum_cm"], visible=False),
        add("rhythm", "rhythm", rhy, "rhythm.json", label=L("Партитура місць", "Score of the places"), visible=False),
        add("sea", "lines", ghost, "ctx_sea.geojson", style="sea", label=L("Море", "Sea"), z=1.2),
        add("coast", "lines", coast, "ctx_coast.geojson", style="coast", label=L("Берег", "Coast"), z=1.3),
        add("roads", "lines", osm["road"], "ctx_roads.geojson", style="road", label=L("Дороги", "Roads"), z=3.4, **infra_g),
        add("pipes", "lines", osm["pipe"], "ctx_pipes.geojson", style="pipe", label=L("Труби геотермальної станції", "Geothermal pipes"), z=3.45, **infra_g),
        add("power", "lines", osm["power"], "ctx_power.geojson", style="power", label=L("ЛЕП", "Power lines"), z=3.45, **infra_g),
        add("walls", "lines", osm["wall"], "ctx_walls.geojson", style="wall", label=L("Захисні вали від лави", "Lava barriers"), z=3.5,
            visible=False, attribution="© OpenStreetMap contributors"),
        add("source", "lines", osm["source"], "ctx_source.geojson", style="source",
            label=L("Джерело: резервуар (еліпсоїд) і точка Моґі", "Source: reservoir (ellipsoid) and Mogi point"), z=3.7, visible=False),
        add("crater_row", "lines", row, "ctx_row.geojson", style="row", label=L("Кратерний ряд Сундгнукур (геол. карта)", "Sundhnúkur crater row (geol. map)"),
            z=3.55, visible=False),
        add("lossband", "lines", bands, "ctx_loss.geojson", style="lossband",
            label=L("Смуги втрати когерентності InSAR, 10.11.2023", "InSAR coherence-loss bands, 10.11.2023"), z=1.4, visible=False,
            note=L("Там, де земля рветься (грабен, тріщини), радар втрачає когерентність: ×4.9 частіше біля кратерного ряду, ніж випадково.",
                   "Where the ground tears (graben, cracks) the radar loses coherence: ×4.9 more often near the crater row than by chance.")),
        add("dike", "lines", dike, "ctx_dike.geojson", style="dike", label=L("Дайка 10.11.2023 · NPE v3 (Sentinel-1 + GNSS)", "Dike 10.11.2023 · NPE v3 (Sentinel-1 + GNSS)"),
            z=3.8, visible=False),
        add("dike_blind", "lines", blind, "ctx_dike_blind.geojson", style="dike_blind",
            label=L("Та сама дайка, сліпа мережа NPE v4", "The same dike, blind NPE v4"), z=3.75, visible=False,
            note=L("Простягання те саме (33° проти 34°). v4 довша на 3 км і зсунута на ~1.5 км на південь: її південний кінець іде в море, де радар не має даних, а апріорне дозволяло дайку до 30 км. Різниця в межах інтервалів обох мереж.",
                   "Same strike (33° vs 34°). v4 is 3 km longer and ~1.5 km further south: its southern end runs into the sea, where the radar has no data, and the prior allowed dikes up to 30 km. The difference is within both networks' intervals.")),
        add("firms", "lines", firms, "ctx_firms.geojson", style="fire",
            label=L("Вогні вивержень 2023–2025, NASA FIRMS VIIRS", "Eruption fires 2023–2025, NASA FIRMS VIIRS"), z=3.6, visible=False,
            attribution="NASA FIRMS (VIIRS 375 m)",
            note=L("9 з 9 вивержень виявлено, перша точка — через 0–1 добу.", "9 of 9 eruptions detected, first point within 0–1 day.")),
        add("places", "lines", osm["places"], "ctx_places.geojson", style="place", label=L("Назви", "Places"), z=3.9,
            attribution="© OpenStreetMap contributors"),
        add("peaks", "lines", osm["peaks"], "ctx_peaks.geojson", style="peak", label=L("Вершини, м над рівнем моря", "Summits, m a.s.l."), z=3.9),
        add("elevation", "elevation", elev, "elev.json", label=L("Висота", "Height")),
    ]
    gnss = [{"station": s, "lon": g["lon"], "lat": g["lat"], "rate_cm_yr": g["val"]["gnss"], "nisar_cm_yr": g["val"]["nisar"],
             "diff_cm_yr": g["val"]["diff"], "z": g["val"]["z"], "note": "independent: not used in the rate, reference or corrections"}
            for s, g in V["gnss"].items()]
    gnss.append({"station": "NYLA", "lon": -22.7381068685, "lat": 63.9735443819, "note": "zero point of the map"})
    man = {"engine_version": 1, "id": SID, "title": L("Рейк'янес", "Reykjanes"),
           "subtitle": "NISAR L-band GUNW, tracks 159A + 082D → vertical motion", "bbox": BBOX, "center": [-22.43, 63.88], "dates": dates,
           "n_pairs": len(dates) - 1, "direction": "ASCENDING+DESCENDING",
           "reference": {"label": "GNSS NYLA (Njarðvík)"},
           "anchor": {"type": "radar-only", "gnss": gnss[:1] + [{k: v for k, v in g.items() if k != "nisar_cm_yr"} for g in gnss[1:]]},
           "trust": {"snr_strong": 3.0, "snr_weak": 1.0, "segments": False},
           "stats": {"velocity": {"strokes": len(vel["features"]), "trust": trust, "band_cm": 2.0,
                                  "headline": {"rate_cm_yr": f"+{H['rate']}", "sigma_cm_yr": H["sigma"], "at_least_cm_yr": H["at_least"],
                                               "lon": H["lon"], "lat": H["lat"], "cum_cm": H["cum"], "p1_cm_yr": f"+{H['p99']}",
                                               "rule": "~0.7 km neighbourhood (3x3 cells of 240 m) with the largest rate - 2 sigma"}}},
           "layers": layers, "notes": [], "generated": f"data to {dates[-1]}", "story": "story.json"}
    pages = story(V, dates, a.math)
    th = theme_for(theme)

    packs = {"data/scenes.json": json.loads((ROOT / "web" / "data" / "scenes.json").read_text(encoding="utf-8")),   # whole index: onward link
             "themes/scroll.json": th, f"data/{SID}/manifest.json": man, f"data/{SID}/story.json": pages}
    for f, obj in files.items():
        packs[f"data/{SID}/{f}"] = obj
    blocks = "\n".join(f'  <script type="text/plain" data-pack="{k}">{pack(v)}</script>' for k, v in packs.items())

    math = MATH.read_text(encoding="utf-8") if a.math else ""
    html = patch_head(head, a.math) + blocks + "\n  <script type=\"module\">" + patch_app(app, a.math) + "</script>\n" + math + tail
    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_text(html, encoding="utf-8")
    print(f"saved {dst}: {dst.stat().st_size / 1e6:.1f} MB")
    for k, v in packs.items():
        print(f"  {k}: {len(pack(v)) / 1e3:.0f} kB")


if __name__ == "__main__":
    main()
