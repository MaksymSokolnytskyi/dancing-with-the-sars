# Space Apps 2026 submission: Edelweiss

Text for the project form, ready to paste. Numbers come from the plates; the plates carry their errors and failed tests.

## Project Name

Edelweiss: an atlas of the moving ground

## Summary (≈190 words)

Edelweiss is an interactive atlas of surface change seen by NISAR, the new NASA–ISRO L-band radar that images the same
ground every 12 days, in daylight or darkness and through cloud. Each plate is one place and one kind of change:
Mexico City sinking into its old lakebed (in places about 30 cm in a year), shallow bowls over coal mines under the
Czech–Polish border, the volcanic dome of Reykjanes in Iceland rising as magma gathers, and the farm year of
California's Central Valley and northern Wisconsin heard in the radar echo. Instead of rainbow interferogram fringes we
draw the data the way 19th-century cartographers drew relief: engraved hachures, ink tints, outlined field parcels and
small signs, so anyone can read the map in seconds. Each plate tells its story in a few short pages. It plays the year
smoothly, makes the feature it talks about blink, and can be narrated aloud. Every plate also says how sure it is. Tests
were written down before we looked at the data, weak results are drawn faint, and failed tests are reported, not
hidden. We built it so that NISAR's measurements can be understood and used by city services, planners and the public,
not only by radar specialists.

## Project Demonstration

https://maksymsokolnytskyi.github.io/edelweiss-atlas/demo/edelweiss-demo.mp4

## Project

https://maksymsokolnytskyi.github.io/edelweiss-atlas/

## Project Details

**What it does.** Edelweiss turns NISAR Level-2 products into plates of an atlas. A plate opens on its story: three
to five short pages, each with one idea, a camera move, the layers it needs, and one or two big numbers. The year can
be played, and the map cross-fades from one 12-day date to the next. The layer a page talks about blinks a few times
so the eye finds it. Questions a reader may have ("How sure are we?", "What failed?") wait in an appendix after the
story. A narrator can read each page aloud and turn the pages by itself. The interface is bilingual (English and
Ukrainian).

**The plates.**
- **Mexico City.** Sinking of up to about 30 cm a year, and it follows the shore of the vanished lake almost exactly.
  The one GPS station that checks the radar measured −27.0 cm/yr; the radar gives −27.8.
- **Ostrava–Karviná.** 22 sinking patches over the coal field. 12 pass a reliability rule fixed before the plate was
  drawn, and about 1,140 buildings stand on the Polish ones. Early, provisional data; the check against Europe's
  ground-motion map came out only partial.
- **Reykjanes.** The dome over the Svartsengi magma body rises at up to 23.5 ± 1.0 cm/yr. A cross-section shows the
  tilt of the towns on its flank, and the 2023 dike is shown as the breath out. Only one of four independent GPS
  stations agrees, and the plate says so.
- **Central Valley.** The farm calendar in 12-day radar coherence: calm in winter, restless in spring, in four rhythm
  types. The rain test failed and is reported.
- **Wisconsin.** HV/HH backscatter of woods and fields through a year. Corn is 3–8 dB brighter than soybean, and the
  harvest is seen in the radar. Four tests failed (lake ice, frost heave, swamp forest, dry afternoons) and are listed.

**How it works.**
- **Reading the data.** A Python pipeline searches the ASF DAAC. It reads only small windows of the NISAR HDF5 files
  over the network (megabytes instead of 2 GB granules), using an Earthdata token.
- **Processing.** It turns GUNW interferogram stacks into vertical motion, with a per-pixel error budget behind every
  "how sure" statement, and computes GCOV backscatter statistics per 1 km block.
- **Drawing.** Our own hachure algorithm traces engraved strokes along the motion field. Context comes from the
  Copernicus DEM (relief), ESA WorldCover (woods, towns) and OpenStreetMap (roads, streets, borders).
- **The map app.** It is static HTML and JavaScript (MapLibre GL + deck.gl, no build step) hosted on GitHub Pages, so
  it opens anywhere without a login.

**Honesty rules.** Before any test that supports a claim on a plate, we wrote the prediction and its pass/fail
threshold into an experiment log. Thresholds were never moved after seeing a result. Anything found afterwards is
labelled as found after the tests. Weak signal is drawn faded or dashed, numbers carry their errors, and failed tests
stay on the plates.

**Benefits and impact.** Ground motion and land change are invisible from the street but matter to water networks,
metros, mines, towns on a volcano and farms. A 12-day NISAR plate drawn this way can be read by a city engineer or a
student without any radar training. The same pipeline runs for any new place from one small config file.

**Creativity.** We replaced rainbow fringes with the visual language of 19th-century survey maps (Dufour, Lehmann)
and political atlases: hachures, ink washes, outlined parcels, engraved signs, a locator globe. Together they read
like the map of a strategy game, and each plate tells its story like a field guide.

**Tools.**
- Python (uv, numpy, h5py, rasterio, scipy, shapely, earthaccess; PyTorch for one research test) and JavaScript
  (MapLibre GL, deck.gl).
- GitHub and GitHub Pages.
- Claude Code as a coding and writing assistant (see SI use).

**Factors we considered.**
- Small downloads and window reads, so the atlas could be built on laptops in 48 hours.
- Colour choices that also part for colour-blind readers.
- No personal data in the repository or in requests to outside services.
- Credits and licences for every dataset.

## Use of Super Intelligence (SI)

Suggested selections: code, data analysis, text, and audio (narration script). No SI-generated images or video. No
NASA branding was generated or modified.

Description to give:

- **Tool.** We used Claude (Anthropic) through Claude Code as an assistant for writing and reviewing Python and
  JavaScript, data processing scripts, documentation and the story text of the plates.
- **Who decided.** Team members set the questions, chose the places, reviewed the code and checked every number shown
  on a plate against the data. The pre-registered tests and their thresholds were decided by the team before the
  results were seen.
- **Narration.** The narrator's script was drafted with Claude and edited by the team. The voice is the viewer's own
  browser speech synthesis. The app says this on the narrator button.

## NASA Data

| Link Text | Link Url |
|---|---|
| NISAR L2 GUNW and GCOV (NASA/JPL and ISRO), ASF DAAC | https://asf.alaska.edu/datasets/daac/nisar/ |
| NASA Earthdata (search and access) | https://www.earthdata.nasa.gov/ |
| NASA FIRMS active fires (VIIRS 375 m) | https://firms.modaps.eosdis.nasa.gov/ |
| NASA HLS, Harmonized Landsat Sentinel-2 (LP DAAC) | https://www.earthdata.nasa.gov/data/projects/hls |
| ASTER GDEM v3 (NASA/METI, LP DAAC) | https://lpdaac.usgs.gov/products/astgtmv003/ |
| IGS ionosphere maps via NASA CDDIS | https://cddis.nasa.gov/ |
| ASF HyP3 processing (Sentinel-1 interferograms of Nov 2023) | https://hyp3-docs.asf.alaska.edu/ |

## Space Agency Partner & Other Data

| Link Text | Link Url |
|---|---|
| Copernicus DEM GLO-30 (ESA / EU, DLR, Airbus) | https://dataspace.copernicus.eu/explore-data/data-collections/copernicus-contributing-missions/collections-description/COP-DEM |
| ESA WorldCover 10 m 2021 (CC BY 4.0) | https://esa-worldcover.org/ |
| Copernicus Sentinel-2 L2A via Element 84 Earth Search | https://earth-search.aws.element84.com/v1 |
| European Ground Motion Service (EGMS), Copernicus Land Monitoring Service | https://egms.land.copernicus.eu/ |
| OpenStreetMap (ODbL), via the Overpass API | https://www.openstreetmap.org/copyright |
| Natural Earth (public domain) | https://www.naturalearthdata.com/ |
| GNSS time series, Nevada Geodetic Laboratory | https://geodesy.unr.edu/ |
| Open-Meteo historical weather (CC BY 4.0) | https://open-meteo.com/ |
| USDA NASS Cropland Data Layer (CropScape) | https://nassgeodata.gmu.edu/CropScape/ |
| U.S. Fish and Wildlife Service National Wetlands Inventory | https://www.fws.gov/program/national-wetlands-inventory |
| Atlas de Riesgos de la Ciudad de México | https://www.atlas.cdmx.gob.mx/ |
| Sigmundsson et al. 2024, Science 383:1228 (reference values) | https://doi.org/10.1126/science.adn2838 |
| MapLibre GL JS | https://maplibre.org/ |
| deck.gl | https://deck.gl/ |
| Claude Code (Anthropic), coding and writing assistant | https://www.anthropic.com/claude-code |
