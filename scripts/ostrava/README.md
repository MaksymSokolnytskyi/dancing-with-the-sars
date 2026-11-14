# Ostrava–Karviná: mining subsidence from NISAR

Research pipeline behind `results_ostrava.html` (a Ukrainian team report, not a jury page): NISAR GUNW over the Ostrava–Karviná
(CZ) and Jastrzębie–Rybnik (PL) coal district, vertical and east–west motion, subsidence bowls, and what OSM says stands on them.
Not part of the `sars` CLI yet. Tests and their status are logged in `notes/experiments.md`, block O.

## Steps (run from the repo root)

| step | what | writes (all under `output/ostrava_karvina/`, ignored by git) |
|---|---|---|
| 1 | `scripts/ostrava/ostrava_karvina.ipynb`: search NISAR GUNW (`NISAR_L2_GUNW_PROVISIONAL_V1`, 17 Jun – 14 Nov 2026), read only the box over HTTPS, stacks, LOS velocity, Isolation Forest, vertical + east from an ascending and a descending track, OSM water / spoil-heap masks, bowl polygons | `stack_*.nc`, `velocity_*.tif`, `anomaly_*.tif`, `zones_*.{csv,geojson}`, `vertical_cm_yr.tif`, `east_cm_yr.tif`, `zones_vertical_clustered.{csv,geojson}`, `osm_water.geojson`, `osm_heaps.geojson`, PNGs |
| 2 | `node scripts/ostrava/ostrava_osm.mjs`: per bowl, buildings, roads, rail, land use, nearest mine and settlement, country (Overpass API) | `zones_vertical_osm.{csv,geojson}`, `zones_vertical_osm_groups.csv`, `osm_summary.json`, `osm_buildings_in_bowls.json` |
| 3 | `node scripts/ostrava/ostrava_map.mjs`: the map (vertical rate, heaps, bowl outlines, buildings) without GDAL | `map_vertical_osm_full.png`, `map_vertical_osm_zoom.png`, `map_vertical_osm_labels.json` |
| 4 | `node scripts/ostrava/build_ostrava_results.mjs`: the results page | `results_ostrava.html` |
| – | `node scripts/ostrava/ostrava_checks.mjs` (optional): odd objects (is M02 one bowl or two, a linear feature, a profile across M01) | `checks_odd_objects.json` |

Step 1 needs Python with the extra packages listed in the notebook's first cell (run it with `uv run --with ...`, no change
to the project lock). Steps 2–4 need Node 18 or newer and nothing else (no npm packages). Overpass answers are cached in
`data/ostrava/osm_cache/` (ignored by git); delete it to refresh OSM. Requests send the User-Agent `space-apps-research`.

## Secrets

Only the Earthdata login, read from the environment: `EARTHDATA_TOKEN`, else `~/.edl_token`, else `~/.netrc`. Nothing secret
is stored in the notebook, the scripts or the outputs. Overpass needs no key.

## Checked

Steps 2–4 were re-run from the repo on the outputs of step 1: the CSV / JSON outputs are byte-identical to the first run and
the page differs only in two labels changed on purpose. Step 1 was not re-run from the repo (it downloads ~57 granule windows);
its code is the first run's code with paths, login and parallel start method changed, checked line by line.
