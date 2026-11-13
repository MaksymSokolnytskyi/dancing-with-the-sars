# Research directions (living list)

Status: **done** / **next** / **open** / **idea**. Details and pre-registrations are in `experiments.md`; mathematics is in `math.md`.

## Product (goes into the app)

| ID | Direction | Status | Why it matters |
|---|---|---|---|
| P-asc | Rebuild the Mexico map and animation from the morning (asc) pass | open (render deferred) | 2× cleaner than the evening pass |
| V1 | Is the error budget complete? | done (R 0.88 at 1 km, 0.45 at 10 km) | Basis for the trust layer |
| P-trust | Per-pixel rate error map | partial (asc usable; desc estimator unstable; over-conservative far from reference) | Trust layer in the hachures |
| P-flag | GOES-based "bad pair" flag per date | idea | Cheap quality gate: σ² ∝ BT std |
| P-iono | Temporally regularised ionosphere correction | idea | Per-pair estimate error random-walks into rates (H8) |
| M6 | Persistent-homology significance per subsidence bowl | idea | Show only bowls above noise |

## Science threads

| ID | Direction | Status |
|---|---|---|
| B5 / E1 / E2 / G1′ | Atmosphere: per-date, GOES-predictable, 1–4 km layer in the evening | done |
| H8 / H8b | Per-pair ionosphere-estimate error is edge noise; correction still needed | done |
| H3b | MMX1 sits on a 2 cm/yr/km slope; part of the GNSS gap | done (weak) |
| H4b / H7 | Persistent asc − desc pattern, height-correlated on slopes | open |
| H7b | Does the slope–height relation flip with season? (troposphere vs geometry) | open |
| T1 / T2 / T3 | Layer depth: radar break 1.9 km evening vs < 0.2 km morning; noise = turbulence × gradient, not moisture amount | done (T1 P2–P3 pass; T2, T3 mostly fail, IGRA too coarse) |
| S1 | Radar + GOES stereo → cloud-top and moisture heights | idea |
| M1.2 | Covariance = Laplacian + low-rank seasonal + diagonal; full error model with no external data | idea |
| M2-SLAM | Pose-graph / loop-closure tools for InSAR networks | idea |
| M1-NT | Bring network-tomography estimators (Fourier-domain, Duffield) to InSAR chains: per-date AND per-pair variance maps | idea (lemma = corollary of Chen-Cao-Bu) |
| P-iono2 | Per-acquisition ionosphere (Fattahi 2017, Piromthong 2026) instead of per-pair GUNW screens | idea, cure for H8 |
| T2-hires | High-resolution soundings / lidar mixing height, Tatarskii C_n^2 per date | open (IGRA too coarse) |

## ML and simulation

| ID | Direction | Status |
|---|---|---|
| D | Okada forward model + simulation-based inference: fault parameters in seconds (Nagqu, Kumamoto) | done, unblinded: SBI fails by the letter (Kumamoto misses Mw by 0.01; Nagqu wrong mechanism); LSQ passes Kumamoto. Next: D-hyb |
| D-hyb | Network as an instant proposal + 60 LSQ steps from its top 16 samples | done: Kumamoto pass (by a hair), Nagqu strike fail vs W-phase plane but 2° from the USGS finite fault; 7 s, ~20× faster than LSQ-128 |
| D-ident | Identifiability map: Fisher spectrum of the Okada → LOS map; strike-slip vs normal confusion | done (math.md M9): condition number 2.9e3 per look, dip stiffest, along-strike position / depth sloppiest; strike-slip imitation is 2–3× worse, so the Nagqu mechanism error was a learning failure |
| D-geom | Sloppy-model follow-ups: (a) Fisher spectrum vs depth/size, power law = multipole boundaries (MBAM); (b) does the SBI net's posterior align with the Fisher sloppy directions (learned space → reduced model) | (a) done post-hoc: two bands 4 + 5, gap ~3.5 ≈ multipole 4; position falls into the sloppy band because of the deplaning gauge (math M9 addendum). (b) ramp degree: gauge story partly refuted. **D-jet (night 10–11.11): Fisher bands = graded ranks of the jet map, and the ranks fill the 2-D monomial staircase (k+2, k+3, …); exact for k = 0, 1, k = 2 only asymptotically. Candidate corollary of flat-limit kernel spectra (Wathen–Zhu, Barthelmé–Usevich), math M9 addendum 3.** Next: literature check; genericity proof sketch; an app rule "how many fault numbers this image can pin down" |
| D-eq | Z₂-equivariant head (even/odd features) as an alternative to variable projection | idea (math.md M8) |
| D-2look | Joint asc + desc SBI / hybrid with magnitude-scaled window (Panama has desc co-event pairs p141, p40) | lite version done (diagnostic): Panama strike 2°, dip 11°, ΔMw 0.13, 13.7 s. Next: freeze the recipe, speed to < 10 s, blind test on the next event |
| D-auto | Auto-trigger: USGS feed → NISAR pair → fault card in seconds, for every M6+ in the NISAR era | in progress: Panama M7.7 as third blind test |
| C | GPU synthetic events + U-Net, injection-recovery on real noise; detection-threshold map | open |
| Style | Dufour style transfer (pet project) | idea, notes in `ml-style-dufour.md` |

## Other data slices (2026-11-11): what else the same tools can read

We have used one product layer (unwrapped phase), one city, two events. The pipeline (GUNW window reader, chain rates, graph noise model, trust map, Okada + SBI, GOES / sounding links) transfers to these:

| ID | Slice | Reuses | Status |
|---|---|---|---|
| X-port | Portability: Central Valley (CA) and Houston, ~1 year of 12-day chains, both passes; rate + trust audit | whole rate + trust pipeline | **downloading** (radar + all NGL GNSS stations in both boxes) |
| X-gnss | Calibrate P-trust against **dozens** of GNSS stations instead of one (MMX1) | P-trust, gnss_check | done: P5 fail (σ collapses in CV); P-tie fail (errors are islands, not smooth) |
| X-season | Central Valley: elastic (seasonal, recovers) vs inelastic (permanent) compaction from the 1-year series; ratio → aquifer storage; compare with USGS well levels | chain series, graph noise model | idea, strong app story: "which part of the sinking is forever" |
| X-canal | Subsidence profile along the California Aqueduct / Delta-Mendota canal (loss of canal capacity) | rate map + trust | idea, concrete stakeholders |
| X-coh | Coherence time series as a change layer: crop cycles, irrigation, harvest, floods, burns, construction | GUNW coherence layer (already read) | idea, no new download |
| X-iono | Ionosphere screens as a space-weather product: per-date TEC gradients via the same graph inversion as G1′ (node/edge split), compared with GNSS TEC | GUNW iono layer, M1 | idea |
| X-hrrr | HRRR 3 km hourly weather model (AWS, CONUS only) at the US sites: T2/T3 with real vertical profiles at overpass time | T1–T3 scripts | idea, cures "IGRA too coarse" |
| X-cities | Global sinking-city scan: same pipeline on many NISAR frames, ranked by rate × trust | everything | later |
| X-gcov | GCOV (backscatter) product: floods and soil moisture, the challenge's other change types | new reader | later |

## Slope failures and debris flows (user idea, 2026-11-11)

| ID | Direction | Status |
|---|---|---|
| L-retro | Retrospective: did NISAR see the Chhochen Khola source slope (26 Aug 2026 ice-rock avalanche → debris flow, Rasuwa) move before failure? asc p98 chain Nov 2025 → Aug 2026 | next (night queue item 6) |
| L-offsets | Precursor via GUNW `pixelOffsets` (offset tracking) instead of phase: fast, icy slopes decorrelate | idea, after L-retro failed on phase |
| L-dpm | Debris-flow / landslide footprint from co-event coherence loss (damage proxy map), with a season guard | post-hoc works on Chhochen Khola (4.8 km² tongue, no pre-event false alarm, one autumn false alarm) |
| L-watch | App layer: slow-moving slopes as downslope-velocity hachures + acceleration flag (inverse-velocity time-to-failure, Fukuzono 1985) + rainfall (GPM IMERG) | idea |
| ML-5 | Slope state: per-slope downslope series + coherence + rain + slope geometry → acceleration probability / time-to-failure posterior; trained on a simulated rain-driven creep model + real noise (same SBI recipe as D) | idea |
| L-load | Channel loading: new scars and deposits from coherence loss / GCOV after storms and quakes; glacial-lake growth from GCOV water masks | idea |
| X-coh | "Rhythm of the fields": cropland vs built-up coherence through the year (Central Valley asc) | Night 10–11.11: X1 (summer lower) PASS, X2 (smooth calendar, lag-1 autocorrelation 0.85) PASS, X3 (rain) FAIL. Post-hoc: one counter-dancing field cluster. Next: check the clusters against USDA CDL (needs OK to download); turn into a demo layer |
| D-conf | Confidence flag for the quake card: agreement among the top-16 proposals | AUC 0.93 on 48 synthetic M7+ events; ≥ 0.9 agreement → 14/14 right. Panama is low-confidence (~0.15). Next: freeze threshold 0.75 for the next blind event |
| H14+ | Per-date error field tied to GNSS: plane per pair (H14) cut Houston desc 17.8 → 6.3 but failed the limits; the field is smooth but not planar | Next: a spline or kriged per-date surface with LOO-chosen stiffness; or learn the per-date field from radar (ML-2) |
