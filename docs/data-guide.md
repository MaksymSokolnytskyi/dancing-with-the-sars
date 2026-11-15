# Data guide

Every dataset the atlas uses: what it is, who publishes it, its licence, how to get it, where it lands on disk, its
format, and which plate needs it.

## Rules

- **No data in git.** `data/` is ignored. The one exception is `web/data/`, the built plates, which the web app and GitHub
  Pages need.
- **One way to fetch.** `uv run sars fetch <plate>` gets everything a plate needs. It is resumable: what is on disk is
  never read again. Start with `--dry-run` to see what it will do and roughly how big it is.
- **Windows, not files.** NISAR files are 1–25 GB each. We read only the plate's box over HTTP, and only the HDF5 chunks
  under it travel. A plate costs MBs per date, not GBs.
- **Pinned pairs.** `sites/<plate>.lock.json` lists the exact NISAR granules a plate is built from. Everyone fetching
  that plate gets the same granules, so everyone builds the same map. `sars fetch <plate> --update` searches ASF again
  and rewrites the lock: commit the new lock only if the plate should move to the new pairs.
- **Credentials stay outside the repo.** Only NISAR needs a login: an Earthdata token in `~/.edl_token` or
  `EARTHDATA_TOKEN`. Everything else is open and anonymous.

## Folder map

```
data/
  raw/
    nisar_cache/bbox_<W>_<S>_<E>_<N>_h<height>/<granule>.npz   GUNW windows (one per pair)        sars fetch
    gcov_cache/bbox_<W>_<S>_<E>_<N>_f<A|B>/<granule>.npz        GCOV windows (one per date)        src/hachure/gcov.py
    nisar_events/<event>/                                       GUNW windows around earthquakes    scripts/fetch_event.py
    nisar_offsets/<name>/                                       GUNW pixel offsets                 scripts/fetch_offsets.py
    copdem_<name>_<res>m.tif, dem/<name>.tif                    Copernicus DEM windows             sars fetch
    context/<plate>_worldcover_10m.tif                          ESA WorldCover clip                sars fetch
    osm/<plate>/<layer>.json                                    OpenStreetMap (Overpass JSON)      sars fetch
    cdmx_atlas/*.geojson                                        Mexico City risk atlas layers      sars fetch
    gnss/<STATION>.tenv3                                        GNSS daily positions               sars fetch
  processed/<plate>/                                            intermediate products (rate maps, series, test outputs)
web/data/<plate>/                                               the built plate: manifest + layers (in git)
```

## Datasets

### Used by the plates

| Dataset | Publisher | Licence | What we take | How to get it |
|---|---|---|---|---|
| **NISAR L2 GUNW** (geocoded unwrapped interferograms) | NASA/JPL and ISRO, via the ASF DAAC | NASA open data, no restrictions; cite the mission | unwrapped phase, coherence, connected components, ionosphere screen, incidence angle; 80 m grid, 12-day pairs | `sars fetch <plate>` (needs the Earthdata token) |
| **Copernicus DEM GLO-30** | ESA / Airbus | free licence with attribution | heights for the relief hachures, terrain shading and summits | `sars fetch` (public AWS bucket, anonymous) |
| **ESA WorldCover 2021, 10 m** | ESA | CC BY 4.0 | land-cover classes: built-up blocks, trees, cropland | `sars fetch` (public COG tiles, anonymous) |
| **OpenStreetMap** | OSM contributors | ODbL | rivers, lakes, glaciers, paths, villages, roads, metro | `sars fetch` (Overpass API) |
| **Atlas de Riesgos CDMX** | Government of Mexico City | open data portal; check the terms before redistributing | geotechnical zones (the old lakebed), mapped ground cracks | `sars fetch` (URL in the plate's `context` entry) |
| **GNSS daily series** | Nevada Geodetic Laboratory | free; cite Blewitt, Hammond & Kreemer (2018) | station MMX1 (Mexico City airport) as ground truth | `sars fetch` (station listed in the plate's `gnss`) |

### Used in research and experiments, not by the CLI

| Dataset | Licence | Used for | How to get it |
|---|---|---|---|
| **NISAR L2 GCOV** (backscatter gamma0, HH and HV; full quad-pol covariance where acquired) | NASA open data | test GC: crop seasons by polarization (`notes/research-resources-2026-11-13.md`); Kakhovka tests K (noise floor, scattering, product audit); Wisconsin quad-pol tests W (reciprocity, co-pol ratio, H/α, Faraday) and story tests WS; Bakota tests B (terrain flattening, noise floor, water area) | `experiments/probe_gcov.py`, `experiments/k_fetch.py`, `experiments/w_fetch.py`, `experiments/ws_fetch.py`, `experiments/b_fetch.py` via `src/hachure/gcov.py` |
| NISAR GUNW pixel offsets | NASA open data | landslide test L-offsets | `scripts/fetch_offsets.py` |
| NISAR GUNW around earthquakes | NASA open data | fault parameter tests (SBI), quake cards | `scripts/fetch_event.py --name ... --lon ... --lat ... --date ...` |
| USGS earthquake catalogue and finite faults | public domain | blind check of our fault estimates | `scripts/fetch_event.py`, `scripts/sbi_compare.py` |
| Open-Meteo daily rain (historical) | CC BY 4.0 | test H9 (Central Valley) | `experiments/fetch_h9_context.py` |
| Copernicus Sentinel-2 L2A (scene classification, red, green, NIR) | Copernicus, free and open; cite "contains modified Copernicus Sentinel data" | water and NDVI classes for the Kakhovka tests K, the Wisconsin tests W and the Bakota tests B; harvest dates and lake ice (WS3, WS4) (radar-independent) | `experiments/k_fetch.py`, `experiments/w_fetch.py`, `experiments/b_fetch.py`, `experiments/ws_fetch.py` (Element 84 Earth Search STAC, anonymous COG reads) |
| Open-Meteo hourly 2 m temperature (historical) | CC BY 4.0 | frozen / unfrozen overpass dates, tests K4 and W4 | `experiments/k_fetch.py`, `experiments/w_fetch.py` |
| Open-Meteo hourly vapour pressure deficit and rain (historical) | CC BY 4.0 | test TW (Wisconsin, Bakota) | `experiments/tw_tests.py` |
| USDA NASS Cropland Data Layer 2025 (Wisconsin state file, 30 m) | public domain (USDA NASS) | corn and soybean classes, Wisconsin story tests WS1, WS4 (radar-independent) | `experiments/ws_refs.py` (CropScape static cache, ~41 MB) |
| USFWS National Wetlands Inventory | public domain (U.S. Fish and Wildlife Service) | swamp forest (PFO) and other wetland classes, test WS2 | `experiments/ws_refs.py` (Wetlands map service query) |
| CODE final global ionosphere maps (IONEX, 1 h) | IGS / AIUB products, free; cite CODE (Dach et al.) | predicted Faraday rotation, test WS9 | `experiments/ws_fetch.py gim` (NASA CDDIS, Earthdata token) |
| NISAR GUNW HH windows, northern Wisconsin (phase, iono screen, coherence) | NASA open data | tests WS6–WS8 (frozen coherence, frost heave, two ionosphere meters) | `experiments/ws_fetch.py gunw` |
| GOES-19 water vapour, IGRA radiosondes | public | atmosphere tests (B6, E2, T2) | `experiments/fetch_goes.py`, `experiments/t2_radiosonde.py` |
| NISAR L2 GUNW PROVISIONAL (Ostrava–Karviná; pairs 24 Jun – 11 Oct 2026 behind the vertical, search window 17 Jun – 14 Nov) | NASA open data (provisional product) | mining subsidence bowls (block O) | `scripts/ostrava/ostrava_karvina.ipynb` (earthaccess search, HTTPS window reads; Earthdata token) |
| EGMS L3 Ortho vertical, release 2019–2023 (European Ground Motion Service, Copernicus Land Monitoring Service) | free with registration; © European Union, Copernicus Land Monitoring Service, EEA | test O3b: EGMS vertical inside and around the Ostrava bowls | `scripts/ostrava/egms_fetch.py` with the personal download key in `EGMS_TOKEN` (tiles E49N29, E49N30, ~150 MB zipped, in `data/ostrava/egms/`) |
| OpenStreetMap via the Overpass API and osmnx (Ostrava–Karviná) | ODbL | water and spoil-heap masks, buildings, roads, rail, land use, mines, borders per bowl | the notebook (osmnx) and `scripts/ostrava/ostrava_osm.mjs` (Overpass, cached in `data/ostrava/osm_cache/`) |
| Sentinel-1 interferograms, Nov 2023 – Jul 2024 (Reykjanes) | Copernicus, free and open; processed by ASF HyP3 (INSAR_GAMMA) | real-dike check of the NPE (10 Nov 2023), null pairs, 2024 events, coherence loss | `scripts/reykjanes/run_s1_dike.py`, `scripts/reykjanes/hyp3_wait.py` (Earthdata token) |
| GNSS daily series, Iceland (SENG, AUSV, GEVK, VOGC, NYLA, VOGS, REYK) | Nevada Geodetic Laboratory; free, cite Blewitt, Hammond & Kreemer (2018) | independent check of the uplift (R1), dike inversion (R3) | `.tenv3` files into `data/gnss/` (`scripts/fetch_ngl.py --name reykjanes --bbox -22.8 63.75 -22.0 64.05`, then copy) |
| NASA FIRMS active fires (VIIRS 375 m, MODIS 1 km) | NASA, open (cite FIRMS) | eruption detection, blind check of the dike trace | `scripts/reykjanes/firms_hotspots.py` (key `FIRMS_MAP_KEY` in the environment) |
| NASA HLS S30 (Harmonized Landsat Sentinel-2) | NASA LP DAAC, open | new lava and vents 2023–2024 | `scripts/reykjanes/hls_timelapse.py`, `hls_dike_optical.py` (Earthdata token) |
| ASTER GDEM v3 | NASA/METI, LP DAAC, open | topography for the layered-crust FE model | tile N63W023 into `data/dem/aster/` |
| Geofabrik Iceland OSM extract | OSM contributors, ODbL | roads, geothermal pipes, power lines, lava barriers, places | `data/osm/iceland-latest.osm.pbf`, then `scripts/reykjanes/osm_extract.py` |
| Sigmundsson et al. 2024, Science 383:1228 (doi 10.1126/science.adn2838) | journal article | reference values for the 2023 dike (R3) | cited; numbers typed in `docs/uk/reykjanes-math.html` |

## What each plate needs

| Plate | Kind | Needs | Status for a fresh clone |
|---|---|---|---|
| `mexico-city-asc` | slow motion (subsidence) | NISAR (28 pairs), DEM, WorldCover, OSM, CDMX atlas, GNSS MMX1 | `sars fetch` + `sars build` |
| `langtang-asc` | sudden change (debris flow) | NISAR (25 pairs), DEM, OSM | `sars fetch` + `sars build` (checked on a clean clone: same result) |
| `central-valley-asc` | fields rhythm (coherence) | NISAR (25 pairs), plus a land-cover grid | **two extra steps**, see below |
| `panama-quake` | earthquake card | quake windows, a trained network (GPU) | research pipeline: `scripts/fetch_event.py`, then `scripts/export_quake.py` with `uv sync --group ml` |
| `ostrava-karvina-asc` | slow motion (mining subsidence) + the study's bowls | NISAR (8 pairs, track 158 A), DEM, WorldCover, OSM, and the Ostrava study outputs for the bowls | run the study first (`scripts/ostrava/`), then `plate_layers.py`, `sars fetch` + `sars build` (below) |
| `reykjanes` (packed page `web/reykjanes-scroll.html`) | uplift + dike | NISAR GUNW 159A + 082D (also 010D, 087A, 111D), Copernicus DEM tile N63W023, GNSS (NGL), Sentinel-1 HyP3, FIRMS, HLS, OSM extract | research pipeline in `scripts/reykjanes/` (below), then `scripts/build_reykjanes_scroll.py` |

**Central Valley, extra steps.** The `fields` layer clusters cropland pixels, so it needs WorldCover fractions on the
plate's rate grid. Both files come from the research pipeline, not from `sars fetch` yet:

1. `uv run python scripts/nisar_gunw_rate.py --name central-valley-asc --bbox -120.20 35.80 -119.00 36.50 --ref-bbox -119.10 36.33 -119.02 36.43 --direction ASCENDING --height 80 --max-pairs 30`
   writes `data/processed/central-valley-asc/nisar_rate.tif`. It reads the windows that `sars fetch` already cached. The
   `fields` layer uses only this file's grid: 80 m, EPSG:32611, 1007 × 1376.
2. `uv run python experiments/fetch_h9_context.py` writes `data/raw/context/cv_worldcover_frac.npz` (and the rain series).

These two commands are reconstructed from the scripts' options and have not yet been run on a clean clone.
`sars fetch central-valley-asc` warns when either file is missing. Folding both steps into the CLI is open work.

**Reykjanes, research pipeline.** Not in the CLI yet. All scripts live in `scripts/reykjanes/` and write to `data/` and
`output/` at the repo root (both ignored by git). Secrets come from the environment only: `EARTHDATA_TOKEN` (or
`~/.edl_token`) and `FIRMS_MAP_KEY`. Order, with what each step writes (details in `scripts/reykjanes/README.md`):

1. `download_site.py reykjanes` → NISAR GUNW windows in `data/reykjanes/gunw/`.
2. GNSS `.tenv3` into `data/gnss/`; the Copernicus DEM tile into `data/dem/`; FIRMS via `firms_hotspots.py` → `data/firms/`;
   OSM via `osm_extract.py` → `data/osm/reykjanes_layers.json`.
3. `reykjanes_ellipsoid.py`, `reykjanes_gnss.py` → source fits; `build_breathing_data.py` → `output/reykjanes/breathing_data.npz`,
   `breathing_vectors.json`; `build_breathing_app.py` → `output/reykjanes/map_layers.json`.
4. Sentinel-1 and ML: `run_s1_dike.py`, `s1_dike_npe_v3.py`, `s1_dike_npe_v4.py`, `coherence_fissures.py`, `results_v4_section.py`
   → `output/reykjanes/ml/*.json`; `run_pipeline.py` runs the NISAR source, FE and real-time dike steps.
5. `uv run python scripts/build_reykjanes_scroll.py` (add `--no-math` for a jury build) → `web/reykjanes-scroll.html`.

Not yet checked on a clean clone. Extra packages used only here: `uv sync --group ml --group research` (torch, zuko, shapely, osmium, pyproj, scikit-fem).

**Ostrava–Karviná, research pipeline.** Not in the CLI yet. A notebook plus plain Node scripts in `scripts/ostrava/`
(order and requirements in `scripts/ostrava/README.md`); outputs go to `output/ostrava_karvina/` (ignored by git). The
Earthdata login comes from `EARTHDATA_TOKEN`, `~/.edl_token` or `~/.netrc`. EGMS (Copernicus, free after registration) is
used for test O3b (`egms_fetch.py`, `egms_o3b.py`; result in block O).

**Ostrava–Karviná plate** (`sites/ostrava-karvina-asc.json`). The motion is the engine's own (one ascending track, as for the
other plates); the bowl outlines come from the study. Steps: the study (steps 1–2 of `scripts/ostrava/README.md`), then
`uv run python scripts/ostrava/plate_layers.py` (writes `data/ostrava/plate/`), `uv run sars fetch ostrava-karvina-asc`,
`uv run sars build ostrava-karvina-asc`. The built plate in `web/data/ostrava-karvina-asc/` is in git, so the map opens
without the study. Checked once on this machine, not on a clean clone.

## File formats

**GUNW window** (`nisar_cache/.../<granule>.npz`, written by `nisar.RawWindow.save`):

| Key | Meaning |
|---|---|
| `phase` | unwrapped phase in radians (NaN = no data) |
| `iono` | ionosphere phase screen in radians (empty if the file has none) |
| `coherence` | interferometric coherence, 0 to 1 |
| `comp` | connected components (0 = invalid) |
| `incidence` | incidence angle in degrees, at the plate's typical height `cache_height` |
| `x`, `y` | projected pixel-centre coordinates (UTM, metres) |
| `epsg`, `wavelength` (m), `days` | grid projection, radar wavelength, pair length in days |

Read one with `nisar.RawWindow.load(path)`, then use `nisar.process(raw, ...)` to get vertical displacement in metres.

- The file name is the granule name.
- The folder name holds the plate box and `cache_height`, so two plates with the same box share windows.
- Size per pair ranges from about 2.6 MB (Langtang) to about 45 MB (Central Valley, a large box).

**GCOV window** (`gcov_cache/.../<granule>.npz`, `gcov.Window`):

| Key | Meaning |
|---|---|
| `hh`, `hv` | gamma0 as linear power |
| `mask` | validity mask (1 = valid in the windows we read) |
| `looks` | number of looks |
| `x`, `y`, `epsg` | pixel-centre coordinates and grid projection |
| `date`, `name` | acquisition date, granule name |

`Window.db("hv")` returns values in dB. Frequency B (80 m) is about 13 MB per date for the Central Valley box.

**Lock file** (`sites/<plate>.lock.json`):
- Contains `site`, `bbox`, `direction`, `track {path, frame}` and `searched` (UTC time of the ASF search).
- `granules[]` lists each pair as `{name, url, ref, sec, days}`.
- Written by `sars fetch --update` or `sars lock`.

**Built plate** (`web/data/<plate>/`):
- `manifest.json` lists the layers and their files. The files are GeoJSON or JSON with physical values only, no colours.
- The web app's contract for layer kinds is in [web/README.md](../web/README.md).

## Adding a dataset

1. Prefer sources that are open, anonymous and readable by window (COG, HDF5 or Zarr over HTTP, an API with a bbox).
2. Write the fetcher so it is idempotent:
   - skip what exists;
   - write to a temporary name, then rename;
   - send the User-Agent `space-apps-research`;
   - never send personal data.
3. Put the output under `data/raw/<dataset>/...` and reference it from the plate config. If a plain URL is enough, give
   the `context` entry a `url` and `sars fetch` downloads it.
4. Add a row to the tables above and to the credits in [README.md](../README.md): publisher, licence, what we use.
5. Pre-register any claim the dataset will support in `notes/experiments.md` before looking at the result.
