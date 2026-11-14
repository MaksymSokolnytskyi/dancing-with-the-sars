# Dancing with the SARs: an atlas of a moving Earth

NASA Space Apps Challenge 2026, challenge [Dancing with the SARs](https://www.spaceappschallenge.org/2026/challenges/dancing-with-the-sars/).

An interactive atlas of surface change seen by the NISAR radar. Each **plate** is one place and one kind of change
(a city sinking, a debris flow, fields through the year, an earthquake), drawn in the visual language of 19th-century
engraved maps (Dufour, Lehmann hachures) instead of rainbow fringes, and told as a few short pages. Every plate says how
sure it is: weak signal is drawn dashed, numbers carry their error, tests that failed are reported.

**New here?** Read this page, then [CONTRIBUTING.md](CONTRIBUTING.md) and [docs/data-guide.md](docs/data-guide.md).
Claude Code sessions start from [CLAUDE.md](CLAUDE.md). Ukrainian onboarding: [docs/uk/README.md](docs/uk/README.md).

## Quick start (the same for everyone)

1. Install [uv](https://docs.astral.sh/uv/) and clone the repo. On Windows clone into a short folder (e.g. `C:\sars`):
   NISAR file names are ~125 characters and Windows stops at 260 per path unless long paths are enabled.
2. `uv sync` (add `--group ml` only if you need the GPU/ML experiments; it pulls PyTorch).
3. Earthdata token, once: log in at https://urs.earthdata.nasa.gov, Generate Token, and save it as **one line** in
   `~/.edl_token` (Windows: `C:\Users\<you>\.edl_token`) or in the env var `EARTHDATA_TOKEN`.
   Use your own account. Never commit the token or paste it into chats, docs or issues.
4. Check: `uv run sars doctor` (add `--deep` to have ASF accept the token), then `uv run python tests/smoke.py`
   (offline, ~15 s).

Then, for any plate in `sites/`:

```bash
uv run sars sites                         # plates, cached pairs, lock, built?
uv run sars fetch mexico-city-asc --dry-run
uv run sars fetch mexico-city-asc         # DEM, NISAR windows, OSM, WorldCover, context, GNSS; resumable
uv run sars build mexico-city-asc         # -> web/data/mexico-city-asc/
uv run sars serve                         # http://localhost:8765/?scene=mexico-city-asc&theme=scroll
```

`fetch` reads only the plate's window of each NISAR GUNW file over HTTP (MBs, not the 2 GB files), caches it in
`data/raw/nisar_cache/`, and never re-reads what is cached. The pairs are pinned in `sites/<id>.lock.json`, so everyone
builds from the same granules. `uv run sars fetch <id> --update` searches ASF again and rewrites the lock (commit it if the
plate should move to the new pairs). Checked on a fresh clone: the NISAR windows and the DEM come out
byte-identical and the built Langtang plate matches. Only OpenStreetMap context can differ a little (it is fetched live,
e.g. 641 vs 649 glacier outlines two days apart).

Other commands: `sars search <id>` (what ASF has over the plate, no login), `sars story <id>` (republish a plate's pages in
seconds), `sars pack <id>` (one HTML file with all data inside, opens from disk), `sars new <id> --bbox W S E N`
(a new plate from the template), `sars lock <id>` (pin what your cache holds).

## Add a place

1. `uv run sars new my-place --bbox W S E N --direction ASCENDING --height <typical terrain height, m>`
2. `uv run sars search my-place`, then `uv run sars fetch my-place`, then `uv run sars build my-place`.
3. Open `http://localhost:8765/?scene=my-place&theme=scroll`, tune the site JSON (layers, reference area, context),
   write the pages in `sites/my-place.story.json`, `uv run sars story my-place`.

All site keys, layer kinds, themes and the page format are in [web/README.md](web/README.md).

## Layout

| folder | what |
|---|---|
| `src/hachure/` | the engine: NISAR access (`nisar.py` GUNW, `gcov.py` backscatter), hachures (`field.py`, `trace.py`), layers for the web (`engine.py`), CLI (`cli.py`) |
| `sites/` | one JSON per plate, its pages (`.story.json`) and pinned pairs (`.lock.json`) |
| `scripts/` | pipeline steps the CLI calls (build, story, pack, OSM, WorldCover, quake card); `scripts/reykjanes/` = the Reykjanes research pipeline (see `docs/data-guide.md`) |
| `web/` | the map app (MapLibre + deck.gl, no build step); `web/data/` holds built plates; `web/reykjanes-scroll.html` is the packed Reykjanes plate, a **test build** with the internal math window (open the file; rebuild with `scripts/build_reykjanes_scroll.py`, `--no-math` for the jury). It carries a frozen, patched copy of `app.js`: engine fixes reach it only after a rebuild |
| `experiments/` | one-off tests behind the numbers in `notes/experiments.md` (pre-registered thresholds, results) |
| `notes/` | research notes and the experiment log (index: `notes/README.md`) |
| `docs/` | data guide (datasets, licences, formats, what each plate needs); `docs/uk/` Ukrainian onboarding |
| `tests/` | `smoke.py`: offline check of the engine, plate configs and built plates |
| `data/` | downloads and caches (not in git) |

## Data and credits

- NISAR L2 GUNW (interferograms): NASA/JPL and ISRO, distributed by the Alaska Satellite Facility DAAC (NASA Earthdata).
- Copernicus DEM GLO-30: © DLR e.V. 2010–2014 and © Airbus Defence and Space GmbH 2014–2018, provided under COPERNICUS by
  the European Union and ESA.
- ESA WorldCover 10 m 2021 (CC BY 4.0): © ESA WorldCover project / contains modified Copernicus Sentinel data (2021).
- OpenStreetMap (ODbL): © OpenStreetMap contributors.
- Mexico City context (lakebed zone, mapped cracks): Atlas de Riesgos de la Ciudad de México.
- GNSS station series (Mexico, Reykjanes): Nevada Geodetic Laboratory.
- Kakhovka tests (`notes/experiments.md`, section K): contains modified Copernicus Sentinel-2 data (2026), via Element 84
  Earth Search; Open-Meteo historical weather (CC BY 4.0).
- Reykjanes plate (`web/reykjanes-scroll.html`): Copernicus Sentinel-1 interferograms of November 2023 processed by ASF HyP3
  (dike check only); NASA FIRMS VIIRS 375 m active fires; NASA HLS (LP DAAC) vents; ASTER GDEM v3 (NASA/METI) for the FE model; reference values from Sigmundsson et al.
  2024, Science 383:1228 (doi 10.1126/science.adn2838).

NASA does not endorse any non-U.S. Government entity. Each non-NASA source is used under its own terms.

## Licence

Code: MIT (see `LICENSE`). Data keep their own licences (above).
