# Research note: the official resources, what we take, where ML fits (13.11.2026)

Sources: the challenge page and its Resources tab (read 13.11), the NISAR Data User Guide (nisar-docs.asf.alaska.edu),
CMR and ASF search (queried 13.11), the NISAR StoryMaps and application white papers, UAVSAR postcards, ARSET, and the
partner pages (CSA, CONAE). No accounts were created, no forms sent. Small window reads only.

## 1. What the challenge text asks, and where we stand

| The page asks | Our answer today | Gap |
|---|---|---|
| Interactive app, NISAR, one or more kinds of surface change, places around the world | Atlas of plates: Mexico City subsidence, Langtang debris flow, Central Valley field rhythm, quake cards | none |
| Open source, user-friendly, "actionable insights" | MIT repo, one CLI, plain-language pages | public URL (GitHub Pages) |
| "Hypothesize the physical drivers" | Mexico: groundwater pumping, old lakebed clay | each plate needs one explicit "why" line |
| "Which polarizations highlight these changes on land and ice?" | GUNW is co-pol (HH) only | **GCOV HH vs HV** (section 4): the direct answer |
| Other sensors / frequencies | none | RADARSAT-1 (C-band, 1998–2008) as a historical baseline; Sentinel-1 is optional |
| Time series, ancillary data (population, vegetation type) | WorldCover, OSM, rain (Central Valley) | population (WorldPop or GHSL) for "people on sinking ground" |
| Human drivers (agriculture, industry) | irrigation calendar (coherence), pumping | say it on the pages |
| Personal or hometown context | none | a Ukraine plate (see 6.4) |
| Sliders, web maps | MapLibre and deck.gl, time slider | none |
| "Future applications of NISAR" | none | a closing page or slide |

## 2. NISAR products: what exists over our places (ASF, 13.11)

Mission facts worth quoting:
- L-band wavelength 23.8 cm, so one interferogram fringe ≈ 11.9 cm of line-of-sight motion.
- 12-day repeat.
- Launched 30 Jul 2025.
- Validated products released 6 Nov 2026 for data from 21 Oct 2026 on; earlier data stays "provisional" until reprocessing in 2027.
- Permanent data gap 27 Jul – 10 Aug 2026.

| Product | What it is | Size per file | Over our places | Use for us |
|---|---|---|---|---|
| GUNW (L2) | unwrapped interferogram, coherence, iono layer, offsets | ~2.3 GB | all four | **in use** (all plates) |
| GCOV (L2) | terrain-corrected backscatter gamma0, HH and HV; grid A 10–20 m, grid B ~80 m | 2–7 GB | Mexico 13, Langtang 22, Central Valley 10, Houston 28 (since June; more back to January) | **new lens**: crops, forest, water, burn scars; answers "which polarization" |
| GOFF (L2) | pixel offsets, 80 m | ~1 GB | Langtang 17 | glacier speed if a plate needs it (Langtang glaciers); Langtang L-offsets test already used GUNW offsets |
| GSLC (L2) | geocoded complex image | ~22 GB | Mexico, Langtang | own interferograms at 5–10 m; too heavy for the hackathon |
| RUNW, RIFG, ROFF (L1) | radar-geometry versions | 0.4–1.3 GB | Langtang | not needed (L2 is easier) |
| SME2 (L3) | soil moisture, 200 m | ~0.1 GB | all four | small and easy: "soil moisture" layer under the field plate; still provisional |
| Urgent Response L1/L2 | ~6 h products after big quakes and volcano alerts | kept 30 days | event-driven | quake cards: say how fast NISAR can answer |

Access recipe:
- **Search:** ASF search API (no login) or CMR.
  - Validated collections: `NISAR_L2_GUNW_V1` (C3622247503-ASF), `NISAR_L2_GCOV_V1` (C3622241360-ASF).
  - Provisional collections: `..._PROVISIONAL_V1` (GUNW C2854335566-ASF, GCOV C2854338529-ASF).
- **Read:** Earthdata bearer token, then a presigned S3 URL, then HTTP range reads of the HDF5. Only the chunks under the plate's window move.
- `src/hachure/nisar.py` (GUNW) and `src/hachure/gcov.py` (GCOV, new) do this; `uv run sars fetch` wraps it.

## 3. What others already show (so we do not repeat it)

- **NISAR StoryMaps** (five, by JPL): floods, wildfires, oil slicks, glaciers, subsidence.
  - All are built on pre-launch airborne UAVSAR or ALOS data.
  - All use the stock ArcGIS scroll template, rainbow fringes and grey backscatter swipes.
  - **The subsidence map already covers the Central Valley and names Mexico City.**
- How we differ, said on the project page:
  - real post-launch NISAR products, a year of 12-day pairs;
  - engraved hachures instead of rainbow fringes;
  - every number carries its error and failed tests are shown;
  - open source and runnable by anyone;
  - our Central Valley plate is about the farming rhythm, not subsidence.
- **UAVSAR postcards:** polarimetric false colour (red and pink = flooded vegetation, black = open water) and rainbow interferograms. This is the visual convention we answer.
- **ARSET:** "Harnessing NISAR" (Jul 2026) builds flood, landslide and volcano maps from GUNW in QGIS; a displacement training runs 12 and 19 Nov. Useful to cite as "the standard workflow is desktop GIS; ours is a web atlas".

## 4. GCOV first look: Central Valley (pre-registered as GC in notes/experiments.md)

What was read:
- 10 monthly GCOV scenes (Jan–Oct 2026), ascending track 34/20, frequency B grid (80 m), HH and HV.
- Only the plate window was read: ~12.8 MB cached per date, ~25 s per read.
- Code: `src/hachure/gcov.py` (search, window read, cache) and `experiments/probe_gcov.py`.
- Picture: `data/processed/central-valley-asc/gcov_gc.png`.

Pre-registered checks:

| check | result | verdict |
|---|---|---|
| GC1: valley-wide cropland HV swings ≥ 3 dB | 2.3 dB | fail |
| GC2: built-up HV stable within 1.5 dB | 1.0 dB | pass: dates are comparable, the calibration is consistent |
| GC3: HV swings more than HH on cropland | 2.3 vs 2.5 dB | fail |

Post-hoc, per pixel:
- A single cropland pixel swings about 5.3 dB in HH and 4.2 dB in HV over the year.
- The noise floor on built-up pixels is about 2.8 dB.
- So single fields clearly "breathe". The valley-wide median hides it, because fields on different calendars cancel out.

Plates:
- **For a plate:** a per-field change map or "field cards", not one valley curve.
- **Polarization story, honestly:** here **HH** carries most of the seasonal signal, not HV as textbooks suggest for canopy.
  - The likely reasons are irrigation, soil moisture and row structure (untested).
  - This is a real answer to the challenge's "which polarization" question, and it comes from data, not the brochure.

## 5. Partner data (CSA, CONAE): verdicts for 48 hours

| Source | Access | Over our places | Verdict |
|---|---|---|---|
| RADARSAT-1 (CSA), 1995–2013, C-band | open COGs on AWS S3 (`radarsat-r1-l1-cog`), Open Government Licence Canada, no account | Mexico City 3 scenes (Dec 2001); Houston before and after Hurricane Ike (12 and 15 Sep 2008); Langtang 1 (2001); Central Valley 1 (1998) | **usable** as a historical baseline ("the same city 25 years earlier"); single scenes, not InSAR pairs; window reads of the COG; georeferencing to check |
| RCM (CSA) | EODMS account; public users get ≥ 16 m over international land; no open demo scenes near us | none ready | not for this hackathon |
| CSA SAR Python tutorial | GitHub notebooks (Sentinel-1, generic) | — | little to reuse |
| SAOCOM (CONAE), L-band | account and licence; open only over Argentina | none of our places | not for this hackathon unless we add an Argentine plate and someone registers |

The "Space Agency Partner & Other Data" field gets: Copernicus DEM and WorldCover (ESA), RADARSAT-1 if used (CSA), OSM, CDMX risk atlas, Nevada Geodetic Laboratory GNSS.

## 6. Ideas that tie it together

1. **Two lenses per place.** GUNW says how the ground *moves*; GCOV says what the surface *is* and how it changes (crops grow, water spreads, forest burns).
   - The legend explains HH (double bounce: buildings, flooded trees) vs HV (volume: canopy, crops).
   - This is the challenge's polarization question answered on the map.
2. **Then and now.** RADARSAT-1 2001 over Mexico City next to NISAR 2026: the same lakebed, two radars, 25 years.
   - Honest framing: C-band amplitude only, no motion from it.
3. **People on moving ground.** Population raster (GHSL or WorldPop, open) × subsidence rate gives "how many people live on ground sinking faster than 10 cm/yr".
   - This is the "actionable insight" line for the jury.
4. **Hometown plate (Ukraine).** The page explicitly invites it. Candidate stories, all need a look at coverage first:
   - Kakhovka reservoir bed after the dam breach (GCOV: dry bed turning to forest; GUNW: settling ground);
   - Kryvyi Rih open-pit mines and tailings;
   - Dnipro city itself.
   - Coverage, checked 13.11 (ASF, no login), Nov 2025 – Nov 2026:

     | Place | GUNW | GCOV | Notes |
     |---|---|---|---|
     | Dnipro | 32 | 38 | |
     | Kakhovka bed | 44 | 54 | |
     | Kryvyi Rih | 40 | 50 | |

     - The best track everywhere is ascending 42/26, with 15 GUNW and 18 GCOV.
     - That is about half of Mexico's density, so expect chain breaks. Rates need the gap-aware rule (H12).
5. **Future applications page.** Urgent Response products, the planned 2027 reprocessing to "validated", and the S-band half of the mission (ISRO).

## 7. Where ML has room (honest)

**Already in the repo:**
- **SBI for quake fault parameters** (`scripts/sbi_*.py`): a CNN plus mixture-density posterior trained on GPU simulations, then a least-squares polish.
- Results are pre-registered and partly failed: the network alone missed by the letter, the hybrid works better.
- This is real ML with an honest scorecard, and it is our strongest "science" story.

**Cheap and real in 48 hours:**

| Idea | Method | Data | Why it fits |
|---|---|---|---|
| Field types from time series | unsupervised clustering (k-means or GMM) of per-pixel HV, HH and coherence curves over a year | GCOV + GUNW, Central Valley | turns the "rhythm of the fields" into a legend of crop calendars; X-coh did it with coherence only |
| Change points | per-pixel Bayesian or CUSUM change detection on GCOV log-ratio | GCOV, Houston and Langtang | "when did this pixel change": floods, landslide scars, new construction |
| Water and flood masks | Otsu or GMM threshold on HH and HV dB, checked against NISAR_WATERMASK | GCOV | the flood story with two polarizations |
| Gap filling with error bars | Gaussian process per pixel over dates | GUNW cumulative series | smooth time slider with honest uncertainty bands, covering the 27 Jul – 10 Aug gap |

**A stretch (only if someone owns it):**
- SAR foundation-model embeddings (e.g. a pretrained ViT on Sentinel-1) for similarity search: "find places that move like Mexico City". Heavy and C-band-trained.
- Learned denoising of interferograms. Hard to validate in 48 h.
- Landslide precursor classifier: we have one event (L-retro), too few labels.

**Rule for all of these:** pre-register the check before running, state the result on the plate.

## 8. Next steps (proposed)

- [ ] Decide whether the atlas gets a GCOV lens (Central Valley HV, or Houston water) for the demo.
- [x] Run an ASF coverage check over Ukraine: covered (section 6.4). Next: `sars new` for one Ukrainian place, then a dry-run fetch.
- [ ] Do a RADARSAT-1 window read over Mexico City 2001 (one COG, window only) and check its georeferencing.
- [ ] Add a population layer (GHSL 100 m) to the Mexico plate: people on ground sinking > 10 cm/yr.
- [ ] Add one "why" line per plate (physical driver) and a "future applications" page.
