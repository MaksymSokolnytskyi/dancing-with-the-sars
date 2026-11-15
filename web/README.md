# Engine → web map: how to add a place, a look, or a layer

Two halves that never know about each other's details:

- **Engine** (`src/hachure/engine.py`, `scripts/build_site.py`, `scripts/export_quake.py`) turns NISAR products
  into **style-free data**: geometry plus physical properties (rate, σ, SNR, chain breaks, time). No colours.
- **Web map** (`web/`) draws whatever the manifest describes, in whatever theme is chosen. Nothing site-specific.

```
sites/<id>.json ──build──▶ web/data/<id>/manifest.json + layer files ──▶ web/app.js + web/themes/<theme>.json
```

Run locally: `uv run sars serve`, open http://localhost:8765/?scene=<id>&theme=<name> (add `&lang=uk` for the Ukrainian interface).

## Add a place (no code)

1. `uv run sars new <id> --bbox W S E N --direction ASCENDING --height <m>` writes `sites/<id>.json` from the template,
   then `uv run sars fetch <id>` fills the NISAR cache (`data/raw/nisar_cache/bbox_..._h<height>/`), the DEM and the context.
2. Edit the site file; the existing plates are worked examples:
   - `sites/mexico-city-asc.json`: slow motion (subsidence). Layers terrain + velocity + grid; GNSS check optional.
   - `sites/central-valley-asc.json`: fields rhythm (coherence through the year). Layer `fields`.
   - `sites/langtang-asc.json`: sudden change footprint (debris flow, landslide, flood edge). Layer `footprint`.
   - `sites/panama-quake.json`: earthquake card. Built by `export_quake.py`; needs a trained network for that area.
3. `uv run sars build <id>` (for a quake: `uv run --group ml python scripts/export_quake.py ...`).
   The scene appears in the menu (`web/data/scenes.json`).

Useful keys: `bbox`, `direction`, `cache_height`, `ref_bbox`, `iono`, `layers`, `terrain{}`, `velocity{}`, `trust{}`,
`grid_km`, `gnss[]`, `anchor` ("radar-only" | "gnss"), `fields{}`, `footprint{}`, `quake{}`, `notes[]`, `dem_path` | `dem_name`.
`landcover` ({src, attribution, trees{class, cell_m, min_share}, built{class, block_m, min_share, spacing_m}}) turns a land-cover
raster (`scripts/fetch_worldcover.py` clips ESA WorldCover to the bbox) into little-tree symbols and a hatch "tone" for built-up land.
These layers carry `by_theme`: they are on only in themes that list them in `layers_on` (the parchment style does).
`relief_info` ({step_m, peaks, window_m, prominence_m, spacing_m}) writes a coarse height grid (`elevation` kind, used by the
frame readout) and labelled summits. Context `points_as: "dots"` draws points as marks without names (stations).
`tropo` ({"mode": "height", "stable_above_m": ...}) removes a height-correlated delay fitted on hills (Mexico: Ajusco).
`context[]`: extra lines from any GeoJSON (any CRS) or raw Overpass JSON (`scripts/fetch_osm.py`), each
`{id, kind: "lines" | "waterlines", src, where?, style, label, visible?, z?, group?, group_label?, note?, attribution?,
label_from?, label_format?, dedupe_key?, min_len_m?, levels_m? (waterlines)}`. `name_key` / `label_from` can be a list
(first found wins, e.g. `["name:en", "name"]`), `latin_only` skips names in other scripts, `min_tag` ({"ele": 6000}) and
`point_format` ("{name} ▲ {ele}") make summit labels. `waterlines` also takes raw Overpass polygons (closed ways, relations).
`footprint.series: true` adds the change lens: patches for every pair (`changes`, polygons with `time: "d"`, theme style
`change_pair`: the engraving is scraped off inside and stippled like scree; the event's patches stay as a scar) and a bars strip
(`change_strip`). `frame_sub` ({en, uk}, with {a} {b} for the years) replaces the cartouche subtitle. `waterlines` turns polygons into engraved
water lines (shore + offset lines that spread out inwards).

## Add a look (no code)

Copy `web/themes/engraving.json`, change colours, fonts and widths, then add the name to `THEMES` in `web/app.js`.
A theme has: `paper/ink/muted/font`, optional raster `basemap`, `terrain`, `velocity` (with alpha per trust class),
named `ramps` (sigma, segments, coherence, los), `categorical` palettes (dance), `polygons` styles (change, fault), `panel`.
A layer asks for a ramp or style **by name**, so a theme only has to provide the names.
Themes now: `engraving` (default), `night`, `scroll` (parchment: `font_load` for a Google font, `texture {grain, stains, vignette}`
drawn in the browser, `frame {band_px}` = graticule border with degree-minute labels, a scale bar and a cursor readout of
position, height above sea level and sinking rate). `?frame=1` shows the frame in any theme.
Line styles (`lines{}`) can be plain, `dash`, `zebra` (old railway sign), `casing`, `fade` (water lines) or `dot` (stations).
Old-map lettering for place names: `italic`, `spaced` (letter-spacing), `caps_rank` + `caps_size` (capitals for the biggest
places), `collide` (names that would touch are dropped, redone on zoom). Trees: `airy` = minimum gap in icon sizes,
so zoomed out only every n-th tree shows. The scene title is set in a cartouche at the top of the frame.
Town blocks (`landcover.blocks` in the site config: WorldCover built-up minus every OSM street, ~18 m cells) come out
as a layer of kind `image` (`file` + `bounds`, white + alpha PNG) that the theme tints: `lines.blocks {color, alpha,
zoom_fade [z0, z1, min]}` — crisp Siegfried-style blocks up close, a lighter Dufour-style town tone zoomed out.
OSM `roads` and `streets` come from `scripts/fetch_osm.py`.

**Partitura (score of the districts).** `rhythm` in the site config ({r_km, min_cells, rows: [{name, lon, lat, ref?}]})
makes a `rhythm` layer: per named place, the median curve of the grid cells within r_km, its rate ± σ, and `still` when
|rate| < 2σ. The front-end draws a strip under the map that shares the time slider: one row per place, one wedge per
interval hanging from the staff (length = how far the ground went down, scaled to 12 days). Hover a row to ring the place on
the map, click to pin it and fly there. With the score on, ▶ plays the year in ~18 s; the ♪ switch (or `?sound=1`) makes
each moving row chime every 2 cm it goes down (`BELL_CM` in app.js), on its own note of a minor pentatonic (fastest = lowest),
and the row flashes at the cursor. Fast districts chime often, slow ones rarely, so together they make a polyrhythm
(Mexico: ~70 chimes a year, 14 for Nezahualcóyotl, 6 for the centre). A chime rings only at a new deepest level, so noise
that lifts a district and drops it back does not ring twice. Still rows stay silent.

## Atlas pages (the plate tells its story)

Page text: the first paragraph is the page's one idea and stands on the card; the rest opens in a detail sheet beside
the card ("Tell me more"; `keep: n` on a page keeps n paragraphs on the card). A "sure" page must say its main limit in
that first paragraph. `cells` layers are drawn as a smooth tint (theme `"cells_as": "tint"`): the block values are read as
one field, bilinear between block centres, and painted into a bitmap with the mode's ramp; categorical modes fall back to
the engraved hatch (theme `hatch`). A theme can also ask for `"hatch"`, `"dots"` or `"squares"`. On motion plates the
ground's total move since the first date (the `cum_cm` series of the `grid` cells) is laid under the strokes as layer
tints, one band per `tint.step` cm (10; a manifest may set `tint_step`, Ostrava-Karviná uses 5): umber for sinking, slate
for uplift (Reykjanes swaps them to match its strokes), darker = further. Contours are then drawn only at band edges, so
lines and tints agree. While motion or cells are on the map, the terrain, context lines and town blocks are drawn at
theme `quiet` strength (0.65 in scroll) and every label is drawn last, over the strokes. Build numbers under the
legend show only with `?debug=1`. The score
strip (layer kind `rhythm`, "the music") is left out of the atlas; `?score=1` brings it back.

A plate is played as pages. Write `sites/<id>.story.json`, point the site at it (`"story": "sites/<id>.story.json"`), then

    uv run sars story <id>      # seconds: copies the pages, no layer rebuild

and reload. Each page: `section` (intro | what | why | sure | next | listen, shown as the kicker), `title`, `text[]`
(paragraphs, `{en, uk}`), optional `facts[{value, unit, label}]` (big numbers), `camera` (`{bbox}` or `{center, zoom}`),
`layers` (the exact set of layer ids or groups that are on; theme-only layers still follow the theme), `modes`
(`{cells layer id: mode key}`), `cut` (cut | scan | block, for a `surface` layer), `day` (start | end | YYYY-MM-DD),
`play` (start the year when the camera arrives), `legend` (show the motion legend on the page), `pin` (a score row name),
`marks[{lon, lat, text, dx, dy}]` (engraved call-outs with a leader line; dx/dy = where the words sit, in pixels).

The plate opens on its pages (`?page=n` for a page, `?story=0` for the free map). ← → or PageUp/PageDown (a presentation
clicker) turn pages, Esc leaves; "explore the map yourself" opens the usual panel and "Read the plate" comes back.
Text is short on purpose: a page should read in about 20 seconds.

## The atlas: cover, chapters, glossary

`data/scenes.json` may carry an `atlas` block (cover words, `{en, uk}`) and `chapters` (`id`, `title`, `intro`, `scenes`
in reading order). With them, `index.html` without `?scene=` (or with `?cover=1`) opens the title page *Edelweiss*:
the engraved flower (`edelweiss.svg`, drawn by `scripts/draw_edelweiss.py`), the three chapters with their plates
(`blurb` per scene), the glossary and the credits. The last page of a plate leads on to the next plate in reading order,
and the last plate back to the contents; "☰ Contents" in the story card and in the panel reopen the cover.

`data/glossary.json` holds the hover help: per term the forms to find (`match.en`, `match.uk`), a title and a short
text. The first use of each term on a page gets a dotted underline; hovering, focusing or tapping it opens the card.
English forms take only a plural ending (`pass` matches "passes", not "passed"); Ukrainian forms are stems.

Publish: `uv run python scripts/publish_site.py` copies the committed `web/` (nothing else) to the public repository
`MaksymSokolnytskyi/edelweiss-atlas`, served at https://maksymsokolnytskyi.github.io/edelweiss-atlas/.

## Send one place as one file

    uv run sars pack mexico-city-asc --theme scroll   # extra pack_scene.py options pass through

writes `share/<scene>-<theme>.html` (~5 MB): all data inside, pickers hidden, opens from disk without a server.
Needs internet only for the map libraries and the font.

## Layer kinds (the contract)

| kind | file | what the front-end does |
|---|---|---|
| `hachures-static` | GeoJSON lines `{dark}` | relief strokes |
| `hachures-time` | GeoJSON lines `{rate_cm_yr, sigma_cm_yr, snr, segments, trust, dark, t[]}` | strokes grow with the time slider; fade by trust |
| `contours` | GeoJSON lines `{level_cm, major, d}` (d = date index) | depth lines of the current date, majors labelled |
| `surface` | `{lon0, lat0, dlon, dlat, nx, ny, frames[date][cell] mm, sigma_mm_yr, gnss[], section}` | cross-section A–B: profile per date, ±2σ band, tilted houses, GPS dots |
| `lines` | GeoJSON lines `{k?}` + points `{name, rank}`; `style`, `z`, optional `group` | context (old lakes, cracks, city); a `group` is one switch; points become labels |
| `cells` | `{columns, cells, cell_km, modes[], series?, dates?, cluster_curves?}` | coloured blocks; mode picker; click shows the cell's curve |
| `polygons` | GeoJSON (any geometry) + `style` | areas, edges and points in the theme style |
| `card` | `{estimate, confidence, threshold, confident, spread, reference?}` | panel card with a confidence flag |
| `rhythm` | `{dates, rows[{name, lon, lat, r_km, series_cm, rate_cm_yr, sigma_cm_yr, still, ref}]}` | partitura strip under the map + rings on the map + optional sound |
| `rhythm` with `kind: "bars"` | `{dates, unit, note, event_date?, rows[{name, bars[per interval]}]}` | the same strip as bars (each row scaled to its own max, the event interval in red), value of the shown pair on the right; no sound |

A cells `mode` is `{key, label, ramp | categorical, time?}`; `time: true` colours by the series value at the slider date.
Labels can be a string or `{"en": ..., "uk": ...}`.

A new layer kind means: one engine function that writes the file + one branch in `buildLayers()` (and a legend entry).

URL switches: `?section=1|scan|block` (cut, scanner, block diagram), `?city=1` (city group on), `?taper=top|bottom`
(wedge thick at the top, or a "drop" thick at the bottom; theme key `velocity.taper`).

Motion strokes come oriented from the high end to the low end; the front-end draws them as wedges (thick at the top)
and dashed when the signal is weak. `velocity.color_up` colours rising ground.

## Honesty rules the map keeps

- Strokes below 1σ are not drawn; 1–3σ are faded (theme decides how much).
- Badge says **"radar only"** unless the site is tied to GNSS; a GNSS check is shown when a station is inside.
- Quake card shows the confidence flag and the spread over seeds; the reference solution is labelled "comparison only".
- Post-hoc rules (footprint, fields clustering) carry a note in the manifest (`notes`), shown in the panel footer.
