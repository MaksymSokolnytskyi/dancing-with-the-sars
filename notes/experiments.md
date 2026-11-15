# Experiments: expectations written BEFORE running (2026-11-10)

Rule: each experiment states up front what result counts as success, partial or failure.
The thresholds are not moved after seeing the data. Results go under each block with a date.

Baseline facts, measured by `scripts/atmo_stats.py` and `scripts/gnss_check.py`:
- 12-day sinking in the bowl: about 0.73 cm.
- Noise in one pair: about 2.0 cm std (1.5 cm after removing the height-linked part and a plane).
- Noise vs distance: 0.44 cm at 250 m, 0.73 cm at 2 km, 1.3–1.9 cm at 10–20 km. Log-log slope 0.34, close to Kolmogorov's 1/3.
- Noise is seasonal: 1–1.5 cm in the dry season, up to 3.6 cm in the wet season.
- MMX1 (airport), year Oct 2025 → Oct 2026:
  - GNSS rate −26.7 cm/yr.
  - NISAR rate −32.8 cm/yr (plain regression over the 28 dates). Bias −6.1 cm/yr.
  - NISAR follows GNSS until March, then runs low by 3–6 cm through the wet season.

---

## A. Noise-weighted rate and an error map
- **What:** fit a line through the 28 cumulative dates, weighting each date by 1/σ². σ per date is solved from the 27 per-pair noise variances (pair noise = noise(date k) − noise(date k+1)). The error per pixel grows with distance from the reference rock, following the measured noise-vs-distance curve.
- **Hypothesis:** the gap with GNSS is mainly noisy wet-season dates, so down-weighting them pulls NISAR toward GNSS.
- **Success:** |NISAR − GNSS| at MMX1 < 3 cm/yr, and GNSS lies within 2σ of the NISAR rate.
- **Partial:** 3–5 cm/yr.
- **Failure:** ≥ 5 cm/yr. Then the bias is systematic (seasonal or geometric), not random, and that is what experiment B is for.
- **Also expected:** σ_rate ≈ 1 cm/yr near the reference and 2–4 cm/yr 15–20 km away. The map pattern stays the same, with no new artefacts.
- **Result (2026-11-10): PARTIAL, close to failure.** `scripts/weighted_rate.py`, `data/processed/mexico-city/weighted_rate.png`.
  - MMX1, descending pass, against GNSS −26.7 cm/yr:

    | Estimate | Rate (cm/yr) | Gap vs GNSS |
    |---|---|---|
    | Chain sum | −32.4 | −5.7 |
    | Plain fit over 28 dates | −31.1 | −4.4 |
    | Noise-weighted fit | −31.2 | −4.5 |

  - Method fix, made openly: the exact per-date noise solve is ill-posed on a chain. It zeroed 6 of the 28 dates and then trusted them completely. That first run gave −28.2, but this is an artefact and must not be quoted. Replaced by a stable estimator (half the mean variance of the two adjacent pairs) plus a chi²-inflated error.
  - The error comes out at 0.9 cm/yr, yet GNSS sits 5σ away. So the gap is systematic, not random noise → experiment B.

## B quick check (2026-11-10): ascending vs descending at MMX1
- Ascending pass (06:05 local time), 29 dates: NISAR −29.8 vs GNSS −27.0 cm/yr, a gap of −2.8. The per-date difference has std 1.0 cm, max 2.9 cm.
- Descending pass (19:13 local time): −32.8 vs −26.7 cm/yr, a gap of −6.1. The per-date difference has std 2.5 cm, max 8.2 cm.
- Both passes go low in the wet season, but the evening pass much more. Wet-season afternoon and evening convection is the likely cause.
- This supports H1: water vapour, dependent on time of day. H2 (geometry) predicts opposite signs on the two passes; both are negative, so H2 is not the main cause.
- The remaining common −2.8 cm/yr gap is still unexplained. It could be the reference rock, the vertical assumption, or the monument vs the radar pixel.
- Caveat: this rests on ONE station. Next, a GNSS-free test: map the asc − desc difference in the wet season vs the dry season.

## B2. GNSS-free check over the whole city (prediction written 2026-11-10, before running)
- **Measure 1, scatter per date.** Fit a line per pixel per pass over the whole year. For each date, take the spatial std of the residual, then the median per season and pass.
- **Measure 2, window rates.** Fit rates separately in the dry window (Nov–Apr) and the wet window (mid-May to Oct) for each pass. Then take the std over the city of the asc − desc difference.
  - Horizontal motion and geometry are the same in both seasons. So only the CHANGE between seasons counts as evidence.
- **H1 (evening convection) predicts:**
  - Measure 1: dry-season scatter similar for both passes (desc/asc ratio 0.7–1.4), wet-season desc/asc > 1.5.
  - Measure 2: the wet-season difference map has std > 1.5× the dry-season one, in patches rather than a smooth tilt.
- **H1 is rejected if** the desc/asc ratio is about the same in both seasons, or the wet-season difference map is no noisier than the dry-season one.
- **Result (2026-11-10): H1 REJECTED.** `scripts/asc_desc_check.py`, `data/processed/mexico-city/asc_desc_check.png`.
  - Per-date scatter:

    | Pass | Dry season | Wet season |
    |---|---|---|
    | asc 06:05 | 0.75 cm | 0.86 cm |
    | desc 19:13 | 1.45 cm | 1.73 cm |
    | desc/asc ratio | 1.92 | 2.00 |

  - The difference-map std is 6.5 cm/yr in the dry window and 6.3 in the wet window (ratio 0.98). It has the same large-scale pattern in both seasons, so it is geometric or systematic, not weather.
- **What holds instead: the evening pass is about 2× noisier ALL YEAR.** This is robust whatever the mechanism. Practical rule: use the morning (ascending) pass for Mexico City.
- **New candidate H4: residual ionosphere after the split-spectrum correction.** Mexico City at 19°N is near the equatorial anomaly, and 19:13 local time is just after sunset.
  - Ionosphere screens in the GUNW files:

    | Pass | Screen std, median | Max | Gradient, 95th pct |
    |---|---|---|---|
    | asc | 1.1 cm | 4.5 cm | 0.27 cm/km |
    | desc | 2.9 cm | 8.8 cm | 0.47 cm/km |

  - The desc screens peak in May–Aug at 5–9 cm. That matches the window where descending runs away from GNSS at MMX1.
  - Against H4: inside the desc track, per-pair ionosphere size correlates only weakly with residual noise (Spearman 0.21, n=27). Pairs with screens > 4 cm have 2.49 cm residual noise vs 1.93 cm for screens < 2 cm.
- **Competing candidate H5:** a more turbulent evening troposphere after daytime heating.
- **Not yet decided between H4 and H5.** Next discriminators:
  - GNSS tropospheric delays at MMX1, morning vs evening variability.
  - Structure-function shape: the ionosphere has longer wavelengths than the troposphere.

## B3. H4 (ionosphere) vs H5 (troposphere): internal tests (2026-11-10)
Script: `scripts/iono_vs_tropo.py`, plus inline checks.

| Test | Result | Points to |
|---|---|---|
| Rescaling the GUNW ionosphere screen by 0–1.2 | Changes only the plane. After plane removal the noise is the same for any scale: asc 0.70–0.72 cm, desc 1.66–1.73 cm | The screen is too smooth to explain the evening excess. "Over-correction" ruled out |
| Noise vs distance | desc/asc ratio rises from 1.4 at 0.5 km to 2.3 at ≥ 10 km | Medium and long scales: either |
| Anisotropy (N–S vs E–W structure function) | desc median ratio 1.00 at 2, 5 and 10 km, isotropic | Against field-aligned plasma bubbles. Favours H5 |
| Height-linked share of the residual | asc 11%, desc 4% | Stratified troposphere is not the evening problem |

- **Current lean:** H5, a more turbulent evening troposphere (deep mixed layer after afternoon heating, possibly the urban heat island). Not proven.
- H4 is not dead either. Medium-scale ionospheric disturbances could look isotropic within a 50 km window.

## B4. External data test (prediction written before downloading)
- **Data:**
  - NGL 5-minute zenith tropospheric delay at MMX1.
  - Open-Meteo reanalysis hourly total-column water vapour, precipitation and CAPE on a 5×5 grid over the window.
- **H5 predicts all three:**
  1. At MMX1, the ZTD change between consecutive acquisition dates is larger at 01 UTC (evening pass) than at 12 UTC (morning), ratio > 1.5.
  2. The reanalysis water-vapour spread across the window is larger at 01 UTC than at 12 UTC.
  3. Per pair, the evening residual noise correlates with the water-vapour spread or change (Spearman > 0.4).
- **H4 is favoured if** ZTD and water vapour are NOT more variable in the evening (ratio < 1.2) while the radar noise still is.
- **Result (2026-11-10): MIXED.** `scripts/tropo_external.py`.
  - **Prediction 1 FAILED.** The 12-day ZTD change at MMX1 is about the same evening vs morning (38.3 vs 36.9 mm, ratio 1.04).
  - **Prediction 2 HELD.** Reanalysis valley spread of the 12-day water-vapour change, evening vs morning: 7.8 vs 5.3 mm of wet delay, ratio 1.46. Boundary layer 340 m vs 35 m (the morning is stable). Rain somewhere in the box in the overpass hour: 59% of evenings vs 18% of mornings.
  - **Prediction 3 FAILED.** No per-pair link. Radar noise vs |ΔZTD|: ρ = −0.20. Radar noise vs reanalysis spread: ρ = −0.12. n = 22, p > 0.3. The reanalysis (~25 km) is too coarse for a 1–30 km field.
  - **Post-hoc, not pre-registered.** The GNSS horizontal delay GRADIENT, which measures spatial unevenness and is exactly what the radar feels, is 2.16× larger at 01 UTC than at 12 UTC. Over the day, 1-hour ZTD std is 1.53× larger. The gradient is at its minimum at 10–16 UTC (04–10 local) and peaks at 22–01 UTC (16–19 local). That matches the radar ratio of about 2 remarkably well, but it was found after looking.
  - **Verdict:** H5 is supported in climatology (3 of 4 diurnal indicators). It is not shown event by event. H4 is still not excluded.

## B5. GOES-19 clouds at overpass time (prediction written before downloading)
- **Data:** GOES-19 ABI band 13 (10.3 µm), CONUS sector, 2 km, the scan nearest each NISAR overpass. Computed over the radar window: cold-cloud fraction (brightness temperature < 235 K = deep convection) and BT std.
- **H5 per-event predicts:** across descending pairs, the radar residual noise correlates with cloudiness on the pair's two dates (max of the two). Spearman ρ > 0.4.
- **If ρ ≈ 0:** the evening noise is not deep-convection-driven. What remains is clear-air boundary-layer turbulence or H4.
- **Result (2026-11-10): SUPPORTED, with a caveat.** 56 of 57 overpasses had a CONUS scan within 2.5 min.
  - Descending pass, radar noise vs cloudiness on the pair's two dates, n = 25:

    | Metric | Spearman ρ | p | Pre-registered? |
    |---|---|---|---|
    | BT std | +0.68 | < 0.001 | yes, passes |
    | Cold-cloud fraction (< 235 K) | +0.38 | 0.06 | yes, just misses 0.4 |
    | Cloud fraction (< 260 K) | +0.12 | – | no, added after |

  - Ascending pass, n = 28: cold-cloud fraction ρ = +0.43 (p = 0.02).
  - Evenings are cloudier: deep convection over > 5% of the box on 26% of evening overpasses vs 7% of morning ones. The median box cloud fraction is 15% in the evening vs 0% in the morning.
  - **Key nuance:** in clear-sky pairs (both dates < 5% cloud) the evening noise is still 1.80 cm (n = 5) vs 0.90 cm in the morning (n = 22). So the evening has a clear-air baseline excess of about 2×, and convection adds on top (2.27 cm in cloudy evening pairs).
- **Picture now:** evening noise = turbulent clear-air boundary layer (or ionosphere, not excluded for the baseline) + convective events. Per-pair cloudiness from GOES predicts radar quality.

## Next hypotheses (opened by B5)
### E1. Does the residual map look like the cloud map? (pre-registered 2026-11-10, before running)
- **What:** for each pair, take the leave-one-out residual map (descending primary, ascending reported).
  1. Average it into 2 km blocks (the GOES pixel size).
  2. Correlate it with ΔBT = BT(second date) − BT(first date) from GOES band 13 on the same grid.
  3. Remove a plane from both.
  4. Shift the GOES map over ±15 km in 1 km steps and look for one common shift for all cloudy pairs.
- **Cloudy pair:** cold cloud (< 235 K) over > 5% of the box on either date.
- **Sign expected:**
  - More cloud and water vapour on the second date → more delay → the radar sees "moving away", so the residual is negative.
  - Colder BT on the second date → ΔBT negative.
  - So the correlation is **positive**.
- **Parallax, a prediction from geometry alone:**
  - GOES-19 sees Mexico City at a 35.5° zenith angle from the south-east (azimuth 127° to the sub-satellite point).
  - A cloud top at height h shows up 0.71·h too far to the north-west.
  - So the best alignment should move the GOES map **toward the SE (127° ± 45°) by 3–12 km**, for tops at 4–16 km.
  - The radar sees low water vapour (below about 3 km), so its own parallax is ≤ 2 km and is ignored.
- **Null:** the same statistic, including the shift search, with ΔBT maps shuffled among the cloudy pairs (2000 permutations).
- **Success, all of:**
  - P1: median correlation at the common best shift is above the 95th percentile of the shuffled null.
  - P2: that correlation is positive.
  - P3: the best shift points SE (82–172°) and is 3–12 km long.
  - Control: clear pairs have |median r| < 0.1.
- **Partial:** P1 and P2 hold, but the shift is ~0 or in a wrong direction. That would mean the moisture is not under the cloud tops, for example spread out or outflow.
- **Failure:** P1 fails. Then GOES predicts *how much* noise a pair has (B5) but not *where* it is, and E as a product needs better data (GPM rain, GNSS network).
- **Gut expectation:** weak but real, median r ≈ 0.2–0.35 in cloudy evening pairs, with the shift roughly SE. Delay comes from water vapour in the lowest kilometres, and a cloud top is only a proxy for that.

### B6. GOES as side information for the per-date noise (pre-registered 2026-11-10, before running)
- **What:**
  - Experiment A had to solve "noise per date" from pair variances. On a chain that solve is ill-posed: n pairs, n + 1 unknowns.
  - Fix: parameterise σ²_date = a + b·g_date, with g = GOES BT std of that date in the box.
  - Fit (a, b) from var(pair k) = 2a + b·(g_k + g_{k+1}). That is 2 unknowns from 27 equations, so well-posed.
  - Weighted line fit per pixel with 1/σ²_date.
  - Compare with the plain fit and with A, at MMX1 and city-wide.
  - City-wide reference: the full-year morning (asc) rate, the cleaner pass.
- **Success, all of:**
  - P0: b > 0, with a bootstrap 95% interval not touching 0. This is B5 restated.
  - P1: MMX1 moves toward GNSS by ≥ 1 cm/yr vs the plain fit (−31.1).
  - P2: std(asc − desc rate map) drops by ≥ 5% vs the plain fit, and by at least as much as A.
- **Gut expectation:**
  - P0 yes.
  - P1 probably **no**. A showed the MMX1 bias is systematic (5σ), and random-date weighting cannot remove a systematic bias.
  - P2 small, about 3–8%.
  - If P1 fails while P2 holds, the bias at MMX1 is something other than random atmosphere: seasonal, horizontal motion, or the reference.

### H6. The morning−evening map difference is horizontal motion toward the bowls (pre-registered 2026-11-10, before running; opened by B6)
- **Why:**
  - Even over the full year, std(asc − desc) of the "vertical" rate maps is 6.6 cm/yr. Per-pixel errors are about 1–2 cm/yr.
  - Both maps assume purely vertical motion.
  - The radar sees east–west motion with **opposite signs** from the two passes: the morning pass looks west, the evening pass looks east.
  - So D = v_asc − v_desc ≈ 2·E·tan(θ), with E the eastward velocity and θ the incidence angle.
- **Physics:** an aquifer-compaction bowl pulls the ground horizontally toward its centre. E ∝ −∂(subsidence)/∂x, i.e. E has the sign of −∂v_up/∂x.
  - This is a field derived from the same gradient we draw as hachures: **hachures = predicted directions of horizontal flow**.
- **Prediction:**
  - P1: corr(D, ∂v_up/∂x) < −0.3 city-wide, both smoothed to 1 km.
  - P2: |corr(D, ∂v_up/∂y)| is smaller than |corr with ∂/∂x|, because the north component is nearly invisible to a polar orbit.
  - P3: GNSS MMX1 east velocity has the sign predicted by D at MMX1, and is within a factor of 2 of D / (2 tan θ).
- **Failure:** corr ≈ 0. Then the difference is atmosphere or processing, not motion.

### H3b. The MMX1 bias is local: GNSS monument vs the 250 m radar average (pre-registered, same time)
- **Why:** B6 and A cannot move the MMX1 bias. Morning and evening are both too negative: −29.9 and −31.1 vs GNSS −26.7. A random-noise problem would not do that.
- **Candidate explanations:**
  - (a) the radar average around the station mixes in faster-sinking soil, while a monument on a building or piles sinks slower;
  - (b) the reference rock is moving;
  - (c) a scale error.
- **Prediction for (a):**
  - The radar rate changes by > 1.5 cm/yr between the nearest pixel and the 250 m and 1 km averages.
  - The local gradient at MMX1 is > 2 cm/yr per km.
- **If the rate is flat around MMX1:** (a) is out, and (b) or (c) remain.

## B. Where the GNSS gap comes from
- **What:** compare the ascending pass (06:05 local time) with the descending pass (19:13) at MMX1.
- **Hypothesis H1:** wet-season water vapour, stratified between the reference rock (high) and the airport (low). Then the two passes have different errors in the wet season, and both agree with GNSS in the dry season.
- **Hypothesis H2:** geometry. Horizontal motion leaks in because we assume purely vertical motion. Then the two passes have opposite-sign errors all year round.
- **Hypothesis H3:** an unwrapping or reference problem. Then a jump appears in specific pairs, not as a seasonal drift.
- **Success:** the data pick one hypothesis clearly.

## C. GPU generator + injection-recovery detector
- **What:**
  1. Generate synthetic small events (Gaussian bowls, 0.3–3 km wide, 0.5–5 cm) on top of noise.
  2. Train a U-Net on synthetic data only.
  3. Test by injecting events into the REAL Mexico pairs, which carry real noise, and count how many are found.
- **Baseline:** Gaussian smoothing plus a threshold, at the same false-alarm rate.
- **Prediction:**
  - At 1 false alarm per scene, the threshold for a 1 km event is about 1.5–2 cm in the dry season and 3–4 cm in the wet season.
  - The network beats smoothing by 1.3–2×.
- **Success:** ≥ 1.3× lower detectable amplitude than smoothing, measured on real pairs.
- **Failure:** < 1.1×. Then the noise is simply too strong for a single pair, and the honest output is the threshold map itself.

## D. Fault parameters in seconds (simulation-based inference)
- **What:**
  1. Run an Okada forward model in PyTorch.
  2. Add noise from the pre-event pairs.
  3. Train a network that outputs fault parameters with uncertainty.
  4. Test on Nagqu M7.0 and Kumamoto M6.8 against USGS solutions.
- **Success:** strike and dip within 15° of USGS, and the centroid within 10 km.
- **Partial:** right on Nagqu (clean, 1 m) but wrong on Kumamoto (sea, weaker signal).
- **Failure:** wrong on both. Then the forward model or the geometry is mis-specified.
- **Not a claim:** earthquake prediction. This is fast characterisation after the event.

## E. Atmosphere as a signal (water vapour at 80 m)
- **What:** compare the residual maps from A with weather at overpass time: rain (GPM IMERG) and GNSS tropospheric delay, if we can get them.
- **Success:** wet-delay anomalies in the residuals sit where convection or rain is, clearly more often than chance. Then the residual maps are a water-vapour product.
- **Limit known in advance:** InSAR sees only differences between two dates, never absolute water vapour.

## Results 2026-11-10 (E1, B6, H6, H3b)

### E1: SUPPORTED on P1 and P2; P3 direction right, length just over the limit; control FAILED, in an interesting way
- **Data:** descending pass, 11 cloudy pairs.
- **Common best shift of the GOES map:** 12.6 km toward 108°. The pre-registered window was 3–12 km toward 82–172°, so the direction fits and the length overshoots by 0.6 km.
  - 12.6 km = 0.71·h gives h ≈ 18 km for pure GOES parallax. That is too high for cloud tops.
  - Post-hoc: the desc radar looks east, so its own parallax pushes moisture east by h_w·tan 44°. Adding 1–3 km of that gives a plausible split: tops about 14 km, moisture about 2 km.
- **P1:** median r = +0.36 at the best shift vs +0.10 at zero shift. Shuffled-null 95th percentile 0.33, p = 0.018. Passes.
- **P2:** positive in 10 of 11 pairs. Passes. The best single pair, 2026-09-02, has r = +0.65.
- **Control (clear pairs, |r| < 0.1):** FAILED. Clear evening pairs give median r = +0.27 even without clouds.
  - Post-hoc reading: in clear sky the 10.3 µm window still senses water vapour and surface cooling. GOES "sees" clear-air moisture too. This opens E2.
- **Ascending (morning):** only 4 cloudy pairs, not significant (p = 1.0). Best shift 6 km toward 99°.
- **Files:** `e1_desc.png`, `e1_asc.png`, `e1_clouds_pattern.json`.

### B6: P0 yes, P1 no, P2 no. As expected, the MMX1 bias is not random noise
- **P0:** σ²_date = −0.02 + 0.323·g, b 95% CI 0.18–0.46. The constant term a is about 0: with GOES as the explanatory variable, essentially all per-date noise variance scales with BT std.
- **P1, MMX1:**

  | Method | Rate (cm/yr) |
  |---|---|
  | plain | −31.1 |
  | A | −31.2 |
  | B6 | −31.0 |

  The change is 0.1 cm/yr, below the 1 cm/yr threshold.
- **P2, std(asc − desc):**

  | Method | std (cm/yr) | vs plain |
  |---|---|---|
  | plain | 6.64 | – |
  | A | 6.57 | −1.1% |
  | B6 | 6.48 | −2.5% |

  B6 is better than A but below the 5% threshold.

### H6: FAILED. The asc − desc difference is not east–west motion toward the bowls
- corr(D, ∂v/∂x) = −0.14. The sign is right but the pre-registered threshold was −0.3. corr(D, ∂v/∂y) = +0.05.
- MMX1 GNSS east velocity −1.0 cm/yr vs radar-implied +1.7 relative to the city median.
- Incidence angles: asc 38.5°, desc 44.0°.
- Open question: what is the persistent 6.6 cm/yr pattern? It was the same in dry and wet seasons (B2). → H4b.

### H3b: WEAKLY SUPPORTED. The station sits on the slope of the bowl
- Local gradient at MMX1 is 2.1 cm/yr per km (the threshold was 2).

  | Ring around MMX1 | asc (cm/yr) | desc (cm/yr) |
  |---|---|---|
  | 0–250 m | −29.9 | −31.1 |
  | 1–2 km | −28.0 | −28.0 |

  GNSS is −26.7.
- So a few hundred metres of placement explain about 2 cm/yr of the 3–4 cm/yr gap. Both passes are biased the same way, which points to a common cause (local, reference, or scale) rather than atmosphere.

## New hypotheses (from the results above)
- **E2 (pre-registered now, data downloading: GOES bands 10 = 7.3 µm water vapour, 15 = 12.3 µm):** clear-sky GOES moisture predicts the radar atmosphere. In clear evening pairs:
  - the split-window change Δ(BT13 − BT15) correlates **negatively** with the residual (more vapour → bigger split → more delay → negative residual);
  - ΔBT10 correlates **positively**;
  - |median r| ≥ 0.25 and above the shuffled null.
  - Applied test: fit one coefficient per pair on the west half of the map and apply it to the east half. Residual variance drops by ≥ 10% (median over pairs).
  - Meaning if true: a GACOS-like correction from geostationary data with 5-minute latency and 2 km pixels.
- **H4b (not yet run):** the persistent asc − desc pattern is a systematic ionosphere-correction bias (evening equatorial anomaly), not motion.
  - Predictions: ≥ 50% of the smoothed D variance is a plane or quadratic; recomputing both maps without the iono correction changes std(D) by ≥ 20%.
- **S1, radar–GOES stereo (post-hoc idea from E1):** the best-shift vector = GOES parallax (cloud-top height) + radar parallax (moisture height, opposite sign for asc and desc).
  - With enough cloudy pairs on both passes, the two unknowns separate: cloud-top height and the effective height of the moisture layer, from a 2-sensor "stereo" without any weather model.
- **H7 (pre-registered now, together with the first half of H4b):** the persistent D pattern is a DEM-error or height leak.
  - Prediction: |corr(D, DEM height)| > 0.3 within the valley.
  - H4b P1 is tested at the same time: ≥ 50% of the 1 km-smoothed D variance is a plane or quadratic surface.
  - If neither holds, D is small-scale and local: unwrapping, masks, or real local differences.
- **Result H7 / H4b-P1: both FAILED as pre-registered, with one post-hoc lead.**
  - Plane share of D variance 0.40, quadratic 0.44. The threshold was 0.5.
  - corr(D, height) in the valley −0.08.
  - Post-hoc: over ALL pixels, mountains included, corr(D, height) = −0.53. The two passes disagree mainly on the slopes, which is a height-dependent bias.
  - Candidates: stratified troposphere with a different seasonal phase morning vs evening that does not cancel in a 1-year fit, or DEM / layover effects on steep slopes that differ by look direction.
  - Scale of D:

    | Smoothed to | std (cm/yr) |
    |---|---|
    | none | 6.6 |
    | 1 km | 3.8 |
    | 10 km | 2.6 |

    Even the 10 km part is about 2× bigger than random atmosphere predicts (about 1 cm/yr), so something systematic remains.
  - Next: split D into valley and slopes. Fit D = k·height on the slopes per season half. Check whether k flips with season (troposphere) or stays fixed (geometry).

## Mathematical view (for spotting links to other fields)
- **The pair chain is a graph.**
  - Dates are nodes and pairs are edges. The atmosphere of a date is a node potential φ_k, and a pair sees φ_{k+1} − φ_k.
  - Per-date noise variances from pair variances means solving |B|ᵀσ² = var(pairs), with |B| the unsigned incidence matrix.
  - Theorem (signless Laplacian): |B| has full column rank iff every component of the graph contains an **odd cycle**. A pure 12-day chain is a path, hence bipartite, hence singular. The null vector is the alternating ±1 vector, exactly why NNLS zeroed every other date in A.
  - Adding 24-day pairs creates triangles (k, k+1, k+2), which are odd cycles, so the per-date noise becomes identifiable **without** GOES. This gives an independent check of B6 (b, a ≈ 0) using only radar. → **G1**.
- **Residual maps are a phase screen, just like in astronomy.**
  - The structure function slope of 1/3 (D ∝ r^{2/3}) is Kolmogorov. For a path-integrated screen, theory gives D ∝ r^{5/3} below the layer thickness and r^{2/3} above it.
  - The crossover scale is therefore a measure of the turbulent layer depth. Saturation at 10–20 km is the von Kármán outer scale L0.
  - This is exactly the telescope "seeing" problem: Fried parameter r0, outer scale L0, frozen flow (Taylor hypothesis), multi-guide-star tomography.
  - Adaptive-optics machinery should transfer directly: phase-screen priors, reconstruction from several "guide stars" (GNSS stations, GOES bands), and wind-driven prediction from the 5-min GOES sequence. → **T1:** crossover scale from the radar = boundary-layer depth, evening > morning (reanalysis BLH 340 m vs 35 m).
- **E1 / S1 are stereo vision.**
  - The best shift is a disparity. The two unknown heights (cloud top and moisture layer) are solved from two disparities with opposite radar look directions, as in epipolar geometry.
- **B6 is heteroscedastic regression with an observable volatility proxy.**
  - GOES BT std plays the role of "realized volatility" (a VIX for the radar), and a ≈ 0 says the variance is proportional to it.
  - This links to state-space / Kalman models with an exogenous, known R_k. That is the natural form for a near-real-time NISAR time-series filter.

### E2 result (2026-11-10): PARTIAL. GOES sees the evening moisture in the 7.3 µm water-vapour band, not in the split window
- **Evening (desc), 14 clear pairs:**

  | Predictor | Median r | Pairs with the predicted sign | Shuffled-null 95th pct | p | Verdict |
  |---|---|---|---|---|---|
  | ΔBT10, 7.3 µm water vapour | +0.30 | 13 of 14 | 0.16 | < 0.001 | passes, sign as predicted |
  | ΔBT13 | +0.27 | – | – | 0.002 | was the E1 control "failure" |
  | Split window Δ(BT13 − BT15) | −0.01 | – | – | 0.88 | fails |

  Hold-out (fit on the west half, apply on the east) reduces variance by only **5%**, median, with 71% of pairs improved. The pre-registered threshold was 10%, so this fails.
- **Morning (asc), 24 clear pairs:** nothing usable. BT10 r = +0.09 (n.s.). BT13 r = −0.14, the opposite sign, probably surface temperature. Hold-out makes things worse (−17%).
- **Reading:**
  - The evening delay noise lives where the 7.3 µm channel looks. Over a 2240 m city, at ~780 hPa, that weighting function peaks roughly 1–4 km above the ground.
  - It does not live in the near-surface vapour that the split window senses.
  - This points to a deep evening convective mixing layer rather than a thin surface layer. It fits T1 (layer depth from the structure-function break) and S1 (moisture-layer height): both should come out at a few km in the evening.
  - As a correction product, GOES is a weak predictor: it tells where the noise is, but explains only about 5% of it at 2 km.

### G1 → G1′. Per-date noise from the radar alone (pre-registered 2026-11-10, before running)
- **Data check:**
  - ASF has only nearest-neighbour GUNW pairs. Desc: 25 × 12 d, 2 × 24 d, 1 × 36 d. Asc: 27 × 12 d, 2 × 24 d.
  - Most longer pairs bridge missing dates, so there is about one triangle per pass.
  - One odd cycle formally makes |B| full rank, but the conditioning is hopeless far from the triangle.
- **Key realisation (graph view):**
  - The edge covariance matrix is C = B Σ Bᵀ, with B the signed incidence matrix and Σ = diag(σ²_date).
  - A only used diag(C) = |B| σ², which is bipartite-singular.
  - The **off-diagonal** of C is −σ² of the shared date for adjacent edges, which identifies each interior date directly.
  - A redundant 24-day pair from the same SLCs is exactly the sum of two 12-day ones, so it carries no information beyond this covariance. It needs no download.
- **Estimator:** σ̂²_{k+1} = −cov_pixels(r_k, r_{k+1}).
  - r = pair − v·span, with v the rate from all pairs except k and k+1.
  - Median removed; also computed with a plane removed.
- **Predictions:**
  - P1: covariances of adjacent pairs are negative (σ̂² > 0) for ≥ 90% of interior dates. If many are positive, adjacent pairs share a common signal (non-linear or seasonal motion), and the white-date model is wrong.
  - P2: Spearman ρ(σ̂_date, GOES BT std) > 0.4 on desc. This is independent of B6, which used pair variances only.
  - P3: σ̂² = a + b·g gives b > 0, and the 95% CI of a includes 0 or a < 0.3 cm². That would be consistent with B6.
  - P4: median evening σ̂ / median morning σ̂ lies in 1.5–2.5 (B2 found ≈ 1.9–2.0).
- **Consistency check, not a pass/fail:** var(r_k) ≈ σ̂²_k + σ̂²_{k+1}.

### G1′ result (2026-11-10): MOSTLY SUPPORTED. The radar alone recovers the per-date noise, and the GOES link is confirmed independently
- **Evening (desc), 26 interior dates:**

  | Prediction | Result | Verdict |
  |---|---|---|
  | P1: σ̂² > 0 for ≥ 90% of dates | 85% (median version), 92% (plane-free) | borderline; the negatives are all ≈ 0 |
  | P2: Spearman ρ(σ̂, GOES BT std) > 0.4 | **+0.56** (p 0.003) | ✓ |
  | P3: b > 0 and a ≈ 0 | b = +0.37 [0.13, 0.61], a = −0.89 [−2.5, +0.6] | ✓; matches B6 (b = 0.32) from a different part of the data |

  - The noisiest date is 2026-05-17, σ² = 10.9 cm² (σ ≈ 3.3 cm), which is also the cloudiest (BT std 17 K).
- **Morning (asc), 27 interior dates:** P1 85%. ρ = +0.14 (n.s.): morning atmosphere and GOES variability are both tiny.
- **P4:** evening/morning median σ̂ = 1.23 / 0.48 = **2.53**, just outside the 1.5–2.5 window. Borderline fail; same direction as B2.
- **Unexpected (post-hoc): the noise splits into a date part and a pair part.**
  - var(pair) / (σ̂²_k + σ̂²_{k+1}) is 1.21 in the evening but **1.76** in the morning.
  - So about 43% of morning pair variance is **pair-specific** (an edge term), not per-date atmosphere. In the evening it is about 17%.
  - Graph language: node noise (atmosphere, telescopes away in chain sums) vs edge noise (decorrelation, per-pair unwrapping, the per-pair split-spectrum ionosphere estimate). Edge noise does **not** telescope: it accumulates as a random walk along the chain.

### H8 (new, from G1′): edge noise random-walks along the 12-day chain and produces the persistent biases
- **Why:** A, B6 and E2 only deal with node noise. The MMX1 bias (both passes too negative) and the persistent asc − desc pattern D (2.6 cm/yr at 10 km, same in dry and wet seasons) are what a random walk of spatially smooth per-pair errors would leave.
- **Test:** the real independently processed 24-day GUNW pairs (desc 2025-11-18→12-12; asc candidate 2025-11-23 → 24 d) close triangles with two 12-day pairs. Closure map = d24 − (d12a + d12b).
- **Prediction, to be written in full before running:** closure std is ≈ the edge-noise level (about 0.5–0.7 cm in the morning) and has smooth large-scale structure, not pixel noise. If closure is ≈ 0, the edge term is pixel-level decorrelation and H8 is out.

### H8 test: triangle closure (pre-registered 2026-11-10, before downloading)
- **Triangles:**
  - desc 2025-11-18 → 11-30 → 12-12, closed by the real 24-day GUNW 11-18 → 12-12;
  - asc 2025-11-23 → 12-05 → 12-17, closed by the 24-day GUNW 11-23 → 12-17.
- **Closure** c = d24 − (d12a + d12b), vertical cm, referenced as usual and median removed.
  - Computed with the GUNW ionosphere correction (as used in all maps) and without it.
- **Expected size if G1′ edge noise is real:** with edge noise σ_e per pair (≈ 0.6 cm morning, ≈ 0.8 cm evening, from G1′), std(c) ≈ √3·σ_e, i.e. about 1.0 cm morning and 1.4 cm evening.
- **Predictions:**
  - P1: std(c) with the iono correction ≥ 0.5 cm on both passes. Below that, the edge term is too small to matter for H8.
  - P2, smoothness: std of c after 2 km Gaussian smoothing ≥ 0.7·std(c raw). The closure is a large-scale field, not pixel speckle. **This is what H8 needs** to bias rates over km scales.
  - P3, source: std(c without iono) < 0.5·std(c with iono). Raw unwrapped phase closes almost exactly; the per-pair split-spectrum ionosphere estimate breaks closure. If true, the ionosphere *estimate* is the main per-pair error that random-walks.
- **Gut:** P1 yes, P2 yes, P3 yes for the evening; less sure for the morning.
- **Caveat known in advance:** one triangle per pass, so this is a case study, not statistics.
- **H8 closure result: P1 split, P2 split, P3 FAILED as stated. Post-hoc reading at the km scale is clear.**

  | | Desc with iono | Desc without iono | Asc with iono | Asc without iono |
  |---|---|---|---|---|
  | std(c), cm | 0.78 | 0.59 | 0.48 | 0.32 |
  | std after 2 km smoothing, cm | 0.49 | 0.17 | 0.35 | 0.09 |
  | Ratio smoothed / raw | 0.64 | 0.29 | 0.72 | 0.27 |

  - P1 (std ≥ 0.5 cm): yes on desc, borderline no on asc.
  - P2 (ratio ≥ 0.7): asc yes, desc no.
  - P3 (raw no-iono/iono < 0.5): 0.76 and 0.66, so it fails.
  - **Post-hoc:** at the scale that biases rates (≥ 2 km), the raw phase closes to 0.1–0.2 cm, and the per-pair ionosphere estimate adds 3–4× more large-scale non-closure. The non-iono closure is pixel speckle (multilook / decorrelation).
  - Order of magnitude (post-hoc arithmetic): a smooth per-pair error of about 0.49/√3 ≈ 0.28 cm, random-walking over 27 pairs in a year, gives a rate error of order 1–1.5 cm/yr at km scales per pass. That is the same order as D at 10 km (2.6 cm/yr, both passes combined) and part of the MMX1 gap.
  - **Key physics:** the *true* ionosphere is node noise (it belongs to a date, so it telescopes away in a chain). The *error of its per-pair estimate* is edge noise (it accumulates). Correcting each pair can therefore make rates **worse** even while it makes each single pair look better.

### H8b. Rates without the per-pair ionosphere correction (pre-registered now, before running)
- **What:** rebuild both rate maps (chain sum, same method for every variant) with iono = False. Compare with iono = True.
- **Predictions:**
  - P1: std(asc − desc) at 10 km smoothing drops by ≥ 15% without the correction.
  - P2: MMX1 moves toward GNSS (−26.7) by ≥ 1 cm/yr on at least one pass.
- **If P1 fails while single-pair noise rises a lot:** the true evening ionosphere is not stationary node noise. A seasonal trend survives, and the correction is needed but should be regularised (plane-only / smoothed).
- **H8b result: P1 FAILED badly; P2 passes on the letter but cannot be interpreted. The correction is needed.**
  - Without the iono correction, std(asc − desc) at 10 km goes from 1.51 to **10.4 cm/yr** (+587%). At 1 km it goes from 3.0 to 12.9.
  - The true ionosphere does not telescope away. A chain sum keeps the end dates' screens, which are several cm in the evening, and the non-stationary part survives.
  - MMX1 asc moves from −29.7 to −26.0 (GNSS −26.7) and desc from −32.4 to −32.0. With a 10 cm/yr map error that asc "match" is luck, not a fix.
- **Picture of the error budget now:**
  - (1) per-date troposphere: node noise, telescopes in chains, predicted by GOES;
  - (2) true ionosphere: node noise but with large end-date and seasonal effects, so it must be corrected;
  - (3) the per-pair iono-estimate error: edge noise, smooth at ≥ 2 km, about 0.2–0.3 cm per pair, random-walks into rates (order 1 cm/yr);
  - (4) pixel speckle: harmless after averaging.
- **Next, highest impact:**
  - **V1: is this error budget complete?** Predict the std of D at 10 km from (1) + (3) alone and compare with the observed 1.5 cm/yr (chain-sum maps). If they agree, the per-pixel error map for the app is trustworthy.
  - Fix to try later: a temporally regularised ionosphere (smooth the per-pair estimates over neighbouring pairs before applying them).

### V1. Is the error budget complete? (pre-registered 2026-11-10, before running)
- **Model (math.md M1):** for each pass, the error of the chain-sum rate is Σ_k r_k / T, with r_k the leave-one-out residual of pair k.
  - If noise is only node (per-date) and edge (per-pair) noise, independent, then cov(r_j, r_k) = 0 for |j − k| ≥ 2.
  - So var(rate error) = [Σ_k var(r_k) + 2 Σ_k cov(r_k, r_{k+1})] / T².
  - Everything is computed on maps smoothed to L = 1 km and L = 10 km (spatial variance, median removed), exactly like D.
- **Prediction:** R(L) = observed var(D_L) / (pred_asc + pred_desc), where D = chain-sum asc − desc (H8b: std 3.00 at 1 km, 1.51 at 10 km).
  - **Complete model:** R in [0.5, 2] at both scales.
  - **Gut:** R ≈ 1–1.5 at 1 km, but R > 2 at 10 km. D was the same in the dry and wet seasons and is height-correlated on slopes, so I expect a persistent long-wavelength component the white model misses.
  - **Diagnostic, reported either way:** the sum of covariances at lags ≥ 2 (should be ≈ 0 under the model).
- **V1 result: PASSES at 1 km (R = 0.88). Just under the window at 10 km (R = 0.45, model OVER-predicts). The opposite of my gut.**

  | Scale | Observed std(D) | Predicted std(D) | R |
  |---|---|---|---|
  | 1 km | 3.00 | 3.20 | 0.88 |
  | 10 km | 1.51 | 2.27 | 0.45 |

  - Per pass, predicted std at 10 km: **asc 2.05**, desc 0.97 cm/yr. The "clean" morning pass has the larger long-wavelength rate error, because its noise is edge-heavy (43%, G1′) and edge noise random-walks (M2).
  - There is no sign of a missing persistent component. The asc − desc pattern D is explained by accumulated noise within a factor of about 2. H4b / H7 "systematics" are not needed to explain its size; the height correlation on slopes stays open.
  - **Caveat:** the prediction is a difference of large numbers (Σvar − |2Σlag1|). The leave-one-out residuals also carry a self-fit constraint that makes lags ≥ 2 negative (−2.8 to −11 cm²). So the 10 km number is fragile, factor-2 level.
  - Cleaner estimator for later: residualise each pass with the *other* pass's rate, which removes the self-fit constraint and adds a known rank-1 term.
  - **Product consequence:** a per-pixel error map can be built from adjacent-pair covariances. It must include the edge random walk, otherwise the morning map's long-wavelength error is underestimated by about 2×.

### D. Fault parameters in seconds: full pre-registration (2026-11-10, before building)
- **Blindness:** the USGS moment-tensor and finite-fault solutions for us6000u43g (Nagqu M7.0) and us6000tgb9 (Kumamoto M6.8) are saved in `data/raw/usgs/` and **not looked at** until the inference is final.
- **Forward model:**
  - Okada (1985) uniform-slip rectangle in an elastic half-space, ν = 0.25, surface displacement (E, N, U), written in PyTorch (batched, GPU).
  - Validated against an independent implementation (cutde, Nikkhoo & Walter triangular dislocations) to < 1% before use.
  - Projected to LOS with each pair's incidence and the look direction (left-looking).
- **Parameters (9):** centroid E/N offset, centroid depth, strike, dip, rake, length, width, slip. Nuisance: a random plane (orbit / iono ramp).
- **Noise:** real NISAR noise from the event's own pre-event pairs (random crops, flips, ±sign). This is the "atmosphere trick".
- **Inference:** CNN encoder → Gaussian posterior head (strike as sin/cos), trained on simulations only. Not a full normalising flow; stated honestly.
- **Baseline:** classic non-linear least squares with the same forward model (multi-start Adam), timed.
- **Success (per event), all of:**
  - strike and dip within 15° of the closer USGS nodal plane;
  - Mw (μ = 30 GPa) within 0.2 of USGS;
  - the surface projection of the fault within 10 km of the USGS epicentre;
  - inference < 1 s on GPU.
- **Calibration on held-out synthetics:** 68% intervals cover 60–76% of truths.
- **Partial:** Nagqu passes, Kumamoto fails (sea, weaker signal, 24-day pair).
- **Failure:** both fail. Then the forward model, geometry or sign conventions are mis-specified.
- **Gut:**
  - Nagqu: strike and Mw pass; dip is the risk (single look direction, dip–slip trade-off ±20°).
  - Kumamoto: 50/50.

### T1. Layer depth from the break in the structure function (pre-registered 2026-11-10, before running)
- **Theory (math.md M3):** for a path-integrated Kolmogorov screen of effective thickness h, the noise std grows as r^{5/6} for r ≪ h and as r^{1/3} for r ≫ h.
- **Measurement:**
  - Structure function of the leave-one-out residual maps at lags 0.16–20 km, median over pairs, per pass.
  - Model: SF²(r) = nugget² + A²·f(r / r_b), a smooth broken power law with fixed slopes 5/6 → 1/3. The nugget absorbs pixel speckle.
- **Predictions:**
  - P1: on both passes, the slope of std between 0.5 and 2 km is larger than between 5 and 20 km (curvature consistent with a break).
  - P2: evening r_b ≥ 1.5 × morning r_b.
  - P3: evening r_b in 1–4 km, matching E2's "moisture 1–4 km above the city".
- **Gut:** P1 yes, P2 likely, P3 uncertain. Note that r_b is an *effective* thickness and depends on the profile shape.
- **T1 result: P1 FAILED (a design flaw: raw slopes include the nugget). P2 and P3 PASS.**
  - Raw std slopes rise with lag on both passes: asc 0.20 → 0.27, desc 0.27 → 0.41. That is the opposite of P1.
    - Reason: the pixel nugget (0.28 / 0.42 cm) flattens the small lags. P1 should have been stated on nugget-corrected slopes. Kept as failed.
  - Fit with nugget and fixed theoretical slopes:

    | Pass | Break r_b | Bootstrap 16–84% |
    |---|---|---|
    | Evening | **1.9 km** | 1.35–3.15 km |
    | Morning | 0.06 km (at the 0.05 km lower bound) | – |

    The morning break is below the 160 m resolution, so it is really "< 0.2 km".
  - The ratio is ≫ 1.5 in 100% of bootstraps. P3: evening r_b inside 1–4 km.
  - **Three independent views agree:**
    - radar structure-function break, evening about 2 km;
    - GOES 7.3 µm channel sensing a 1–4 km layer (E2);
    - morning layer < 0.2 km vs reanalysis morning boundary layer of 35 m.
  - Reanalysis evening BLH (340 m) is lower. Reading: at 19:13 the daytime mixed layer (2–4 km over Mexico City in the afternoon) has just stopped being driven but is still there as a residual layer full of moisture structure. The radar sees the residual layer, not the instantaneous BLH.
  - **Radar as a boundary-layer meter:** one number per pass from the noise itself, at no extra cost.
  - New direction **T2:** per-pair r_b vs reanalysis BLH at 18 and 19 local (when the layer collapses), and vs radiosonde mixing height if available (Mexico City airport sondes at 12 UTC and 00 UTC, near both passes!).

### T2. Weather balloons at overpass time (pre-registered 2026-11-10, before running)
- **Data:** IGRA2 radiosondes from Mexico City airport (MXM00076679). The 00 UTC launch is 1 h before the evening pass; the 12 UTC launch is 5 min before the morning pass.
- **Per sounding:**
  - wet refractivity N_w = 3.73e5·e/T² → zenith wet delay ZWD;
  - h50 = the height above the station containing half of the ZWD. This is the "wet-delay scale height", the balloon's version of r_b.
- **Predictions:**
  - P1: median h50 is evening ≥ 1.5 × morning.
  - P2: evening median h50 within 1–4 km. Same layer as T1 (r_b 1.9 km) and E2 (1–4 km).
  - P3: across evening dates, Spearman ρ(σ_date from G1′, h50) > 0.3. A deeper moist layer means more turbulent delay.
- **Gut:** P1 yes. P2 probably about 1.5–2.5 km. P3 uncertain, because noise also depends on how *patchy* the moisture is, which one balloon cannot see.
- **T2 result: P1 FAILED, P2 PASSES (but is not specific), P3 FAILED. A useful refinement.**
  - Wet-delay scale height h50: evening 1.74 km, morning 1.59 km, ratio 1.10. The moisture sits at the same height in both.
  - ZWD: evening 88 mm, morning 73 mm.
  - Evening σ_date vs h50: ρ = −0.05. vs ZWD (the amount of moisture): ρ = **0.00**.
  - **Reading:** radar noise is not about *how much* or *where* the moisture is. It is about *turbulence that makes it patchy*.
    - The morning has the same moisture column, but stably stratified (a laminar sheet → almost no lateral structure).
    - The evening has the same moisture, stirred into patches up to about 2 km (T1 r_b = 1.9 km).
    - Consistent with B4: the GNSS ZTD 12-day change is the same morning and evening, ratio 1.04.
- **T3 (pre-registered now, before running): balloon temperature profile → mixed-layer depth.**
  - Parcel method: z_mix = the lowest height where potential temperature exceeds the surface θ + 1 K.
  - P1: evening median z_mix ≥ 3 × morning.
  - P2: evening median z_mix within 1–4 km. This is the T1 r_b, measured by a completely different instrument.
  - P3: across evening dates, Spearman ρ(σ_date, z_mix) > 0.3.
  - Gut: P1 yes (morning inversion), P2 likely, P3 maybe 0.2–0.4.
- **T3 result: P1 FAILED, P2 FAILED, P3 PASSES (weakly). The main limit is the balloon data's vertical resolution.**
  - z_mix: evening 0.30 km, morning 0.23 km.
  - At 00 UTC (18 local, sunset) the surface is already cooling, so the surface-parcel method hits the nascent surface inversion.
  - The IGRA file holds only mandatory and significant levels; the first level above the surface is about 250 m up. The 0.2–0.3 km values are at the resolution limit.
  - P3: ρ(σ_date, z_mix) = +0.38 (p 0.06, n 26). Passes the 0.3 threshold, not the conventional significance.
  - Post-hoc (residual-layer top, θ above 0.5 km + 1 K): 2.72 km for both passes with almost no spread. That is a fixed coarse level, so it carries no information.
  - **Conclusion:** IGRA levels are too coarse for boundary-layer structure. Needed: high-resolution soundings (BUFR, about 10 m steps), or ceilometer / lidar mixing heights in Mexico City (RUOA / UNAM network?). → data TODO.

### P-trust. Per-pixel rate error map + calibration (pre-registered 2026-11-10, before running)
- **Estimator (math.md M1, single realisation per pixel):** Q(x) = Σ_k r_k(x)² + 2 Σ_k r_k(x)·r_{k+1}(x), with r_k the leave-one-out residual of pair k.
  - Under the node + edge model, E[Q] = var(chain-sum error at x).
  - Smooth Q with a 1 km Gaussian (estimator noise only) and clip at > 0. σ_rate(x) = √Q / T.
  - Computed per pass, relative to the same reference rock, so it grows with distance from the reference by itself.
- **Calibration test:** z(x) = (v_asc − v_desc) / √(σ_asc² + σ_desc²) on common pixels (chain-sum rates, pixel level).
  - P1: |z| < 1 for 55–80% of pixels (68% ideal).
  - P2: |z| < 2 for 88–98% of pixels (95% ideal).
  - P3: no strong trend of |z| with distance from the reference: the median |z| in the 0–5 km and 20+ km rings differs by < 1.5×.
- **Gut:** slightly over-conservative (V1 over-predicted at 10 km), so |z| < 1 for about 75–80%. Real horizontal-motion differences add a few outliers.
- **P-trust result: PARTIAL. P1 and P2 PASS overall; P3 FAILS. The descending estimator is unstable.**
  - Calibration on 306k common pixels:

    | | Measured | Window | Verdict |
    |---|---|---|---|
    | \|z\| < 1 | 79% | 55–80 | pass, just inside |
    | \|z\| < 2 | 95% | 88–98 | pass |
    | std(z) | 1.21 | – | – |
    | Median \|z\|, 0–5 km vs 20+ km from the reference | 0.73 vs 0.43 (ratio 1.70) | ratio < 1.5 | **fail** |

    Near the reference the map is about right (|z| < 1 for 64%). Far away it is over-conservative (82%).
  - **σ maps:**
    - asc: median 3.9 cm/yr at pixel level (80 m, speckle included). Rings: 2.0 cm/yr at 0–10 km → 4.4 at 20–40 km. Q ≤ 0 on only 2.5% of pixels: usable.
    - desc: Q ≤ 0 on **22.8%** of pixels (clipped; for example the 5–10 km ring median is 0.10 cm/yr, which is nonsense).
    - Why: on desc the node noise is large and nearly cancels (V1: 81.7 − 80.7 cm²), so a single-realisation estimate of the small difference is swamped.
  - Product file: `data/processed/mexico-city-asc/nisar_rate_sigma.tif` (morning rate + 1σ).
  - **Improvement for later, needs a new pre-registration:**
    - estimate node and edge variances locally (G1′ in windows) and assemble σ from them instead of the raw Q;
    - or use a structure-function model σ(distance to reference) × local nugget.
- **Trust-hachure prototype (2026-11-10, morning pass, one still, not pre-registered: design, not a test).**
  - `scripts/trust_hachures.py` → `data/processed/mexico-city/trust_compare_asc.png` (left plain, right trust).
  - Per-stroke SNR = |smoothed rate| / smoothed σ. SNR ≥ 2 → normal red; 1–2 → pale, thin; < 1 → not drawn.
  - Of 23,789 red strokes: 58% sure, 16% pale, 26% dropped.
  - What disappears is exactly what looked suspicious in the plain map: the "flowers" in the mountains (SW, W, N edges), where coherence is low and the distance from the reference is large. The bowls remain intact.
  - To do: pale strokes are nearly invisible at this scale, so try dashes or jitter instead of colour. Also add a legend line ("drawn only where signal > 1σ").

### D result, part 1 (2026-11-11): a bug found before unblinding, and the fix (frozen before looking at USGS)
- **Bug: both networks ignore the sign of the image.**
  - v2 (signed slip components u1, u2 instead of rake) still gave posterior mass on the correct slip sense = **0.50** on synthetics. That is coin-flip.
  - Diagnosis on 512 synthetics:
    - posterior-mean vs truth correlation: e0 0.87, n0 0.86, depth 0.82, L 0.81, W 0.68, dip 0.64, strike (sin/cos) 0.6, **u1 0.07, u2 0.02**;
    - flipping the sign of the input image leaves every output unchanged (correlations identical to 2 decimals).
  - So the network learned only sign-blind features: where and how big, not which way.
  - Not a physical degeneracy: tested exactly that rake+180 negates the field (|a+b| = 0) and that strike+180 (with or without rake+180 / −rake) gives a different field (relative difference 0.4–3).
  - Most likely cause: a symmetric saddle. The data distribution is invariant under (image → −image, slip → −slip), and the network sits in the symmetric solution (see math.md M8).
- **Fix (v3, `scripts/sbi_post3.py`), made before unblinding: variable projection.**
  - Okada displacement is linear in the slip vector: los = s_ss·G_ss(geometry) + s_ds·G_ds(geometry).
  - The network gives only the geometry (7 nonlinear parameters). For each geometry sample, (s_ss, s_ds) is solved exactly by 2×2 least squares.
  - Slip uncertainty: the same simulator noise (real pre-event maps × flips × random scale + power-law turbulence) is pushed through the same solve.
  - Synthetic check (200 held-out, Nagqu simulator):

    | | v2 | v3 |
    |---|---|---|
    | rake-sense mass, all | 0.50 | **0.71** |
    | rake-sense mass, larger events (n = 146) | 0.50 | **0.85** |
    | 68% coverage, rake | – | 88% (too wide) |
    | 68% coverage, slip | – | 60% (slightly too narrow) |

    The remaining sense errors come from geometry uncertainty (strike) propagating into the slip solve.
- **Deviation from the pre-registration (stated):**
  - The posterior head is an 8-component Gaussian mixture, not one Gaussian.
  - v3 inference = 8 ms network + exact solve. With 20,000 samples it takes **8.8 s**, which breaks the "< 1 s" criterion as run. With 2,000 samples it would be < 1 s. The criterion is judged on the run as made: **fail on time** for v3 at 20k samples.
- **Nagqu, frozen (v3, before unblinding).** median [16, 84]:

  | Parameter | v3 SBI | LSQ (128 starts, 141 s) |
  |---|---|---|
  | strike | **212** [210, 215] | 194 |
  | dip | **49** [39, 59] | 46 |
  | rake | **+9** [−2, +21]: left-lateral 98% | −96 (normal) |
  | Mw | **6.97** [6.86, 7.04] | 6.99 |
  | centroid E, N | 13.3, 3.6 km | 7.9, 5.7 km |
  | depth | 5.2 km | 7.8 km |
  | L × W | 37 × 9.5 km | 15 × 17 km |
  | slip | 3.5 m | 4.9 m |

  - Fit: SBI per-sample rms median 18.8 cm, best 14.7 cm. LSQ 11.4 cm. Data std 21.2 cm.
  - The two methods **disagree on the mechanism** (strike-slip vs normal), and LSQ fits better. One uniform-slip rectangle fits Nagqu poorly either way.
  - Gut before unblinding: the LSQ normal mechanism is more likely right, and the SBI geometry posterior is too narrow (mode-collapsed).
- **Kumamoto, frozen (v3, before unblinding).** median [16, 84]:

  | Parameter | v3 SBI | LSQ (128 starts, 172 s) |
  |---|---|---|
  | strike | **213** [212, 214] | 206 |
  | dip | **67** [60, 75] | 74 |
  | rake | **179** [167, 191]: right-lateral 98% | 178 |
  | Mw | **6.59** [6.50, 6.67] | 6.61 |
  | centroid E, N | 1.8, −0.6 km | 4.9, 2.2 km |
  | L × W | 34 × 8 km | 21 × 11 km |
  | slip | 1.1 m | 1.5 m |

  - Fit: SBI per-sample rms median 3.5 cm, best 2.8 cm. LSQ 1.8 cm. Data std 5.0 cm.
  - SBI and LSQ agree on the mechanism: right-lateral, steep, NE–SW strike.
  - v2 network-only rake was −109 with a 16–84 range spanning all 360° (the prior). v3 pins it.
  - Synthetic check (Kumamoto simulator): sense mass 0.65 all / 0.81 larger events; 68% coverage rake 88%, slip 66%.
  - **Timing:** 2,000 samples 0.82 s (< 1 s, pass); 20,000 samples 3.0 s. The saved JSON is the 20,000 run. Verdict on time: pass at 2k, fail at 20k. The summaries agree to the second decimal.
  - Gut before unblinding: Kumamoto passes strike / dip / Mw. Mw might be a touch low (Mw 6.8 catalogue).

### D result, part 2 (2026-11-11): unblinded against USGS. **FAIL by the letter for SBI; the forward model is fine; the network is the weak link**
- USGS W-phase moment tensors:
  - Nagqu: planes 175 / 57 / −119 and 40 / 43 / −54 (normal), Mw 7.0.
  - Kumamoto: planes 212 / 75 / −164 and 118 / 74 / −16 (right-lateral on the NE–SW plane), Mw 6.8.

  | Event, method | Strike Δ | Dip Δ | Mechanism | Mw (USGS) | Fault–epicentre | Verdict (all of 4 + time) |
  |---|---|---|---|---|---|---|
  | Nagqu SBI v3 | 37° ✗ | 8° ✓ | strike-slip, **wrong** | 6.97 (7.0) ✓ | 6.2 km ✓ | **fail** |
  | Nagqu LSQ | 19° ✗ | 11° ✓ | normal ✓ (−96 vs −119) | 6.99 ✓ | 0.1 km ✓ | fail (strike, by 4°) |
  | Kumamoto SBI v3 | 0° ✓ | 8° ✓ | right-lateral ✓ | 6.59 (6.8) ✗ by 0.01 | 0.2 km ✓ | **fail** (Mw, by 0.01) |
  | Kumamoto LSQ | 6° ✓ | 1° ✓ | right-lateral ✓ | 6.61 ✓ | 1.9 km ✓ | **pass** (not held to the 1 s limit; 172 s) |

- **Calibration window (60–76%):**
  - pass: e0, n0, depth, dip, L;
  - fail: slip 83–85% and strike 89–93% (too wide);
  - v3 rake 88% (too wide).
- **Pre-registered branches:**
  - "Failure: both fail → forward model, geometry or sign conventions mis-specified." The SBI half of this is met by the letter.
  - The diagnosis does **not** follow: least squares with the *same* forward model, LOS geometry and signs passes Kumamoto on all four criteria and gets the Nagqu mechanism right.
  - So the conventions are right. The failure is in the network: the sign blindness (fixed by v3), plus a too-narrow geometry posterior on Nagqu (strike 212 [210, 215] while the truth is 175).
- **My gut was half right.** I expected Nagqu to be easier. It was the harder one: M7 normal faulting in Tibet, one look direction, and one rectangle fits it poorly (LSQ rms 11.4 of 21.2 cm).
- **Mw:** both methods are about 0.2 low on Kumamoto. That is common for geodetic vs seismic moments (shallow slip deficit, uniform slip, μ = 30 GPa). Not a bug.
- **What this means for the product.**
  - Use the network as a fast *proposal* (milliseconds), then polish with a few LSQ steps started from its top samples. Simple hybrid; expected to keep the speed and LSQ's accuracy.
  - Next test, pre-register first: LSQ started from the best 16 SBI samples, under 5 s total, judged on the same criteria.
  - Better simulator realism for Nagqu-like cases: distributed slip, two segments.

### D-hyb. Network as an instant proposal + short least-squares polish (pre-registered 2026-11-11, AFTER unblinding D, BEFORE running the hybrid on real data)
- **Caveat stated up front:** we have now seen the USGS answers. To keep this honest, every setting below was fixed on synthetics only, and the real data is run once. No reruns with other settings count.
- **Method (`scripts/sbi_hybrid.py`):**
  - 2,000 network geometry samples; slip solved exactly (v3); keep the 16 with the lowest misfit.
  - Polish all 16 together: Adam, lr 0.05, 60 iterations, same forward model, parameterisation and bounds as the 128-start LSQ baseline; report the best.
  - **H1 (primary):** the 16 network starts.
  - **H2:** the same 16 + 16 random prior starts (insurance against a too-narrow network posterior).
- **Settings chosen on 6 noisy synthetics (Kumamoto simulator):**
  - 60 iterations at lr 0.05 gave the same errors as 300 at lr 0.02, in 7 s instead of 32 s.
  - Typical errors: strike < 2°, dip < 5°, Mw within 0.03. One scene in 6 missed strike by 30° and Mw by 0.5 in both settings (a wrong basin).
  - Timing: ~0.1 s per iteration, dominated by GPU kernel-launch overhead in the Okada graph, not by arithmetic.
- **Time budget:** < 10 s total (network + solve + polish). The original "< 1 s" criterion is **not** claimed for the hybrid. It stays met by the network alone (8 ms) and by v3 at 2k samples (0.8 s).
- **Success per event:** the same 4 criteria as D (strike and dip within 15° of the closer USGS plane, Mw within 0.2, fault within 10 km of the epicentre) + time < 10 s.
- **Predictions:**
  - Kumamoto: H1 passes all (it lands where LSQ-128 landed: strike 206, dip 74, Mw 6.61).
  - Nagqu: **H1 fails.** The network posterior sits at strike 212, the polish stays in that basin and keeps the wrong (strike-slip) mechanism. H2 lands near LSQ-128 (strike ~194, normal) and therefore *also* fails strike by ~19°.
  - Misfit: H2 rms ≤ LSQ-128 rms + 10% on both events.
- **What would surprise me:** H1 on Nagqu switching to a normal mechanism. That would mean the polish crosses basins, and the network's strike error was cosmetic.
- **D-hyb result (2026-11-11, one run per event): Kumamoto PASS (H1, by a hair); Nagqu FAIL on strike as predicted. But the prediction about *why* was wrong.**

  | Event, arm | Strike Δ | Dip Δ | Mechanism | Mw Δ | Fault–epicentre | Time | Verdict |
  |---|---|---|---|---|---|---|---|
  | Kumamoto H1 | 6° | 0° | right-lateral ✓ | 0.196 | 2.0 km | 6.9 s | **pass** (Mw within 0.004 of the limit) |
  | Kumamoto H2 | 7° | 0° | ✓ | 0.208 | 2.1 km | 7.1 s | fail (Mw by 0.008) |
  | Nagqu H1 | 19° | 10° | **normal ✓** | 0.003 | 0.2 km | 6.9 s | fail (strike) |
  | Nagqu H2 | 19° | 11° | normal ✓ | 0.001 | 0.0 km | 6.9 s | fail (strike) |

  - Misfit equals LSQ-128 on both events (1.8 and 11.4 cm), in 7 s instead of 141–172 s, so ~20× faster. The H2 prediction passes.
  - **The surprise I named came true.** Nagqu H1 started from the network's strike-slip samples at strike 212, crossed basins in 60 steps and ended on the normal mechanism (rake −97). So the network's error was in the "nearby" sense: it put us in the right region, and the local polish fixed the mechanism.
  - Honest note: the Kumamoto pass is within 0.004 of the Mw limit, and H2 misses it by 0.008. Read this as "Mw is about 0.2 low, right at the threshold", not as a clear pass.
- **Post-hoc (not part of the verdict): USGS's own finite-fault model for Nagqu uses strike 192, dip 50, Mw 6.95.**
  - Our hybrid / LSQ: strike 194, dip 46, Mw 7.00. That is 2° and 4° from USGS's finite fault.
  - The pre-registered reference (the W-phase point-source nodal plane, strike 175) is itself 17° away from USGS's finite fault.
  - So the "19° strike failure" is mostly the spread between two USGS products. Kumamoto's finite fault (212 / 75, Mw 6.75) agrees with its moment tensor, and our Mw is 0.15 below it.
  - Lesson for future pre-registrations: name the reference product precisely and allow for its own uncertainty (a point-source nodal plane vs a finite fault can differ by ~15–20° for M7 events).

### P-port. Does the pipeline travel? Central Valley (CA) and Houston (pre-registered 2026-11-11, while downloading, BEFORE looking at any rate map or GNSS series)
- **What I have seen so far:**
  - the station lists (names, positions, date ranges);
  - four log lines from the first Central Valley asc pairs: 92–98% valid pixels, single-pair medians −0.4 to +2.1 cm.
  - Nothing else.
- **Data:** full 12-day chains, about one year, both passes, the same code and settings as Mexico (coh_min 0.3, iono correction on, plane removed). No retuning.

  | Site | Morning pass | Evening pass |
  |---|---|---|
  | Central Valley | asc p34 f20, 05:08 local | desc p114 f70, 18:48 local |
  | Houston | asc p62 f16, 05:43 local | desc p113 f73, 19:10 local |

  - Reference boxes fixed in advance: CV = Sierra foothills near Lemon Cove; Houston = SE Pasadena.
- **GNSS:**
  - NGL IGS20 daily vertical; linear rate over the same first-to-last chain dates.
  - Only stations with ≥ 70% daily coverage of that window.
  - σ_GNSS = max(0.2 cm/yr, 3 × formal error), where the ×3 allows for coloured noise.
- **NISAR at a station:** median of valid pixels within 250 m; vertical-only assumption, as in Mexico.
- **Offset:** our reference area may move, so one constant offset per map is removed (the median of NISAR − GNSS over stations). Nothing else is fitted to GNSS.
- **Predictions:**
  - **P1 (accuracy):** RMS over stations of (NISAR − GNSS − offset) is ≤ 2.5 cm/yr in CV and ≤ 1.0 cm/yr in Houston, for each pass.
  - **P2 (pattern):** correlation of NISAR vs GNSS rates over stations is ≥ 0.9 in CV (large range of rates) and ≥ 0.7 in Houston (small range).
  - **P3 (does the "morning is cleaner" rule travel?).** Noise = median over pairs of the spatial std of the single-pair residual (pair displacement − rate × Δt), on 2 km blocks.
    - **Houston:** evening / morning ≥ 1.5 (humid, convective evenings, like Mexico).
    - **CV:** evening / morning < 1.5 (dry, stable air most of the year).
    - So the rule is about climate, not about clocks. If both ratios come out ≥ 1.5, it is about clocks after all (boundary-layer cycle everywhere).
  - **P4 (asc vs desc):** std of (asc − desc) on common pixels after 1 km smoothing is ≤ 2 cm/yr in CV and ≤ 1.5 cm/yr in Houston.
  - **P5 (trust map), only once the trust module is packaged (plan point 2):** z = (NISAR − GNSS − offset) / √(σ_NISAR² + σ_GNSS²); |z| < 1 for 55–85% of station-pass pairs.
- **Gut:**
  - CV max subsidence −15 to −30 cm/yr around Corcoran / Tulare Lake bed.
  - Houston max −2 to −4 cm/yr in the NW (Katy, Cypress, Spring, Conroe), and about 0 in the SE (where pumping stopped decades ago).
  - P1 passes in CV and is at risk in Houston (small signal, humid air).
  - P3 is the one I am least sure of.
- **Failure means:**
  - P1 fails in CV: unwrapping or decorrelation over farmland breaks the pipeline.
  - P1 fails in Houston: the atmosphere swamps small signals, so the app needs the trust layer before it shows anything there.
- **P-port, first part (2026-11-11): Central Valley morning pass FAILS P1 and P2.**
  - With the pre-registered plane removal, 13 stations: RMS 3.78 cm/yr (limit 2.5), correlation 0.849 (needed 0.9). GNSS rates span −16.8 to +1.1 cm/yr.
  - Without the plane (transparency only): RMS 4.26, r 0.77, offset +7.0.
  - Station by station: the deep-bowl stations (CAFP, HANF, QCAP, CAHA, GNSS −11 to −14 cm/yr) come out 4–8 cm/yr too shallow in NISAR. The edges agree within about 2.
  - Single-pair noise on 2 km blocks: 2.2 cm (Mexico morning was about 0.9).
- **Post-hoc diagnosis (not pre-registered):**
  - Sum of the per-pair changes at each station, by season (NISAR vs GNSS, cm):

    | Station | Nov–Mar, with iono | Nov–Mar, no iono | May–Sep, with iono | May–Sep, no iono |
    |---|---|---|---|---|
    | HANF | **+11.3** vs −0.4 | −0.3 vs −0.4 | −3.1 vs −5.1 | −4.7 vs −5.1 |
    | QCAP | **+10.3** vs −0.4 | +2.2 | −3.7 vs −6.1 | −12.9 |
    | CAHA | **+10.8** vs −0.3 | +3.1 | −4.8 vs −5.7 | −13.1 |
    | P566 | **+9.4** vs −1.3 | +1.3 | +1.2 vs +0.2 | −1.3 |
    | TLBN | +4.2 vs −1.4 | **−35.0** | −7.0 vs −8.0 | −11.6 |
    | P056 | +1.6 vs −0.5 | **−19.6** | 0.0 vs −2.7 | −2.3 |

  - **Reading:** with the per-pair ionosphere correction, the valley floor gets about +10 cm of fake winter uplift. Without the correction that bias is mostly gone at the valley-floor stations, but large long-wavelength ramps appear (−20 to −35 cm at the low-coherence stations on the edges).
  - So the correction is needed for the ramps, and its own per-pair error accumulates along the chain. This is exactly Mexico's H8/H8b, now 10× bigger and seasonal.
  - Per-date atmosphere cannot do this: in a chain sum it telescopes (≤ a few cm end to end). This is edge noise / bias.
- **New hypothesis H9 (from this, to be pre-registered before testing):**
  - The GUNW ionosphere screen (split-spectrum, dispersive part) leaks a non-ionospheric, frequency-dependent surface phase over wet farmland: soil moisture and vegetation scattering differ slightly between the sub-bands.
  - Prediction: the winter bias is larger over irrigated or wet fields than over towns and orchards, and steps up after rain dates.
  - The product cure is P-iono: keep only the smooth part of each pair's iono screen (plane / low order), or estimate the ionosphere per date over the network (P-iono2).
  - Must be judged on sites it was NOT tuned on: Houston and Mexico GNSS.

### H9-P. Keep only the plane of each pair's iono screen (pre-registered 2026-11-11, BEFORE running the variant anywhere; Houston rates not yet seen in any variant)
- **Variant:** `--iono-plane` (`nisar.process(iono="plane")`): subtract the least-squares plane of the GUNW ionospherePhaseScreen instead of the full screen. Everything else is as in P-port.
- **Development set (already seen, so it cannot confirm anything):** Central Valley. Expect P1 RMS ≤ 2.5 cm/yr on both passes with the variant.
- **Confirmation sets (judged):**
  - **Houston (primary, unseen):**
    - on each pass, P1 RMS(plane) ≤ RMS(full screen);
    - the variant passes P1 (≤ 1.0 cm/yr) on at least the morning pass.
  - **Mexico:**
    - D = std(asc − desc) at 10 km smoothing falls from 2.6 to ≤ 2.0 cm/yr;
    - the MMX1 morning gap |NISAR − GNSS| gets no worse than +0.5 cm/yr.
- **Gut:** big help in CV, some in Mexico, unclear in Houston (small signal; the iono screen may hold real small-scale ionosphere at 30° N).
- **Failure meaning:** if the plane variant is worse in Houston, the small-scale part of the screen is real ionosphere there. Then the cure is not "drop it" but "estimate it per date over the network" (P-iono2).
- **Correction to the H9-P pre-registration, written before running the variant.** The Mexico baseline "D at 10 km = 2.6" came from the H7 computation, not from the H8b code I will reuse. The H8b code gives 1.51 cm/yr for the full screen. I keep the *relative* criterion I meant (2.0 / 2.6 = a 23% drop), so the Mexico pass threshold is **D at 10 km ≤ 1.16 cm/yr** with the H8b code. The MMX1 criterion is unchanged (asc gap no worse than +0.5).
- **H9-P results so far (2026-11-11):**
  - **Mexico (confirmation set):**

    | | Full screen | Plane only | Criterion |
    |---|---|---|---|
    | D at 10 km | 1.51 | 1.44 | ≤ 1.16 → **fail** |
    | D at 1 km | 3.00 | 2.91 | – |
    | MMX1 morning (GNSS −26.7) | −29.7 | **−26.9** | gap 3.0 → 0.2 → **pass** |
    | MMX1 evening | −32.4 | **−28.9** | – |

    The long-standing MMX1 gap (A/B, B6, H3b) **almost vanishes** on the morning pass and shrinks by 3.5 on the evening pass. This is one station, so it is not proof. But it is the first change that moved that number at all (B6 weights: no; H3b slope: explained ~2 of 3–4).
  - **Central Valley morning (development set):** RMS 3.78 → **2.78** (limit 2.5, still fail), correlation 0.849 → **0.911** (pass). HANF is still +5.4 too shallow; P540 on the west edge is now −6.2 (leftover curvature of the ionosphere, which a plane cannot catch).
  - **Houston: pending** (primary judgement).

### H9. Is the winter iono-screen bias in Central Valley tied to wet farmland and rain? (pre-registered 2026-11-11 01:20, before downloading land cover or rain data)
- **Quantity.** For each Central Valley morning pair k, take the small-scale iono residual s_k = ionospherePhaseScreen − its least-squares plane.
  - Convert it to the vertical displacement it *adds* to the corrected map: b_k = −LOS_SIGN·λ/(4π)·s_k / cos(incidence), the same conversion as `nisar.process`, applied to −s_k.
  - Season sums: B_win = Σ b_k over pairs ending Nov–Mar; B_sum = Σ over pairs ending May–Sep.
  - Everything on 2 km blocks over the whole box.
- **Land cover:** ESA WorldCover 2021 (10 m), cropland fraction per 2 km block (class 40). Built-up = class 50.
- **Rain:** Open-Meteo archive daily precipitation at the box centre (36.15 N, 119.6 W). A pair is "wet" if either acquisition date had ≥ 2 mm of rain in the 3 days before it.
- **Predictions:**
  - **(a)** corr(B_win, cropland fraction) over blocks ≥ +0.3, and |corr(B_sum, cropland)| ≤ half of that.
  - **(b)** cropland-dominated blocks (> 60% cropland) have mean B_win at least 3 cm higher than built-up-dominated blocks (> 30% built-up).
  - **(c)** median over pairs of the block std of s_k (in cm, as b_k) is ≥ 1.5× larger in wet pairs than in dry pairs.
- **Gut:** (a) yes, (b) yes, (c) unsure: the iono screen is heavily smoothed, so a rain step may be blurred.
- **Alternative it tests against:** a real winter ionosphere structure. That would not follow land cover, and it would show the same pattern over the foothills and towns.
- **H9 result (2026-11-11 ~01:40): FAILS on (a) and (b). (c) passes but is confounded. The "moisture over farmland leaks into the iono screen" idea is rejected.**
  - (a) corr(B_win, cropland) = **−0.01** (needed ≥ +0.30); corr(B_sum, cropland) = −0.13. 2,200 blocks of 2 km.
  - (b) cropland − built-up winter mean = **−0.21 cm** (needed ≥ +3). 645 cropland blocks, 82 built-up.
  - (c) block std of the small-scale iono residual: wet pairs 1.10 cm (n = 7), dry 0.54 cm (n = 18), ratio 2.05 → pass.
    - Caveat: wet pairs are mostly winter. The largest values (1.3–2.8 cm) are Nov–Mar, with some dry winter pairs just as large. So (c) cannot separate rain from season.
  - **What this says:** the small-scale part of the iono screens sums over the winter to only ±1.6 cm (block std), with mean ≈ 0. The +10 cm winter bias at the valley stations is therefore **long-wavelength (tens of km)**: the difference between the reference box on the east edge and the valley centre 30–60 km away.
  - That fits a split-spectrum estimation error (noisy, then heavily low-pass filtered, so the error is smooth and big-scale). It does not fit surface leakage. It also explains why plane-only removal helped only partly: a plane cannot catch curvature at this scale.
  - **Next idea (P-iono3, not yet registered):** quadratic surface instead of plane, or a reference network (several stable areas, or GNSS ties) instead of one corner box. The latter is standard practice and would remove long-wavelength error regardless of its cause.

### Trust-audit module (agreed plan point 2), 2026-11-11 ~01:55
- `src/hachure/trust.py` + `scripts/trust_audit.py --site … --pass … [--iono plane]` → `data/processed/trust/<site>-<pass>/{report.md, audit.json, sigma.tif}`.
- **Validation on Mexico asc:** it reproduces the earlier separate scripts: per-date noise 0.48 cm (G1′: 0.484), edge share 44% (G1′: ~43%), per-pixel σ median 3.94 cm/yr (P-trust: 3.9).
- **P-port P5 (pre-registered, judged now that the module exists), Central Valley morning, full iono screen: FAIL.** |z| < 1 for 31% of stations (window 55–85%). Plane-only variant: 46% (still fail).
- **Why, and it matters for the product:** the per-pixel Q estimator collapses in CV.
  - Its median σ is 0.11 cm/yr, which is exactly the clip floor, so Q ≤ 0 on most pixels.
  - CV per-date (node) noise is large (median 1.5 cm, 3× Mexico morning) and the edge share is small (11%). Q is then a difference of nearly equal numbers, as on Mexico desc (V1, M1.3), but worse.
  - Deeper point: the dominant CV error is a **systematic** per-pair iono error with one sign through the winter. Any self-consistency estimator (Q, G1′) partly absorbs it into the rate, because a steady bias looks exactly like deformation. **The radar alone cannot certify itself against a common-mode drift.** That needs an outside anchor: GNSS ties, or an independent iono estimate.
- **Consequence for the app:** the trust layer must combine (1) radar-only σ, where it is valid (Q > 0 and the node / edge model OK), and (2) a GNSS-tie check where stations exist, and say which one it is using.

### D-3. Third blind earthquake: Panama M7.7, 2026-10-09 (us6000u18k) (pre-registered 2026-11-11 ~02:05, before fetching any NISAR data or the USGS solution)
- **Why this event.** Of the M ≥ 6.5 events in the NISAR era, it is the shallow (10 km) onshore one with a clean ascending co-event pair: p76, 10-05 → 10-17.
  - Rejected: Palu M6.7 (no NISAR pairs); Venezuela M7.2 (an M7.5 doublet 32 s later, 75 km away); Flores M7.8 (offshore).
- **What I know going in:** epicentre (−80.757, 7.607), catalogue M7.7, depth 10 km. The magnitude is not blind, as it wasn't for Nagqu / Kumamoto. The mechanism and planes are not looked at.
- **Known risks, stated in advance:**
  - (1) An M7.7 rupture can be ~100–150 km long, so the 1° box (±0.5°) may clip it. The prior allows L ≤ 80 km.
  - (2) An M6.6 aftershock (us6000u19w) ~2.5 h later, 76 km west, just outside the box, may leak deformation into the west edge.
  - (3) Tropical vegetation: lower coherence.
  - (4) Heading −12° is the 32° N value. At 7.6° N it is about −11°; the 1° difference is negligible.
- **Method (frozen, unchanged from D / D-hyb):**
  - `fetch_event.py --half 0.5`, then train `sbi_fault.py` (8,000 steps, v2 net);
  - `sbi_post3.py` (v3 network + exact slip) and `sbi_hybrid.py` (H1 primary, H2), run once each.
- **Reference, named precisely (lesson from D-hyb):**
  - **primary:** the USGS W-phase moment tensor (Mww) nodal planes and Mww;
  - **secondary, reported but not judged:** the USGS finite-fault model, if one exists.
- **Success (hybrid H1):** strike and dip within 15° of the closer Mww plane, |ΔMw| ≤ 0.2, fault within 10 km of the epicentre, hybrid time < 10 s.
- **Gut:**
  - Location and Mw pass, unless the rupture is clipped by the box (then Mw comes out low).
  - Strike / dip 50 / 50.
  - v3 network alone: right mechanism with probability ~0.6.

### P-tie. GNSS-anchored rate maps: remove a smooth surface fitted to NISAR − GNSS, judged leave-one-station-out (pre-registered 2026-11-11 ~02:15; Houston not yet seen)
- **Why:** H9 and P5 show that the dominant error in CV is a smooth, one-sign drift (per-pair iono estimate) that the radar cannot see in itself. The standard cure is to anchor the map to GNSS.
- **Method:** for each site and pass, take the P-port NISAR rate (full iono screen, plane removed) at the stations.
  - Leaving each station out in turn: fit a surface to (NISAR − GNSS) at the other stations, then predict the left-out station.
  - **Surface:** plane (3 parameters) for Houston and CV. Quadratic (6 parameters) only if ≥ 15 stations remain, so Houston only.
  - Report the LOO RMS of (corrected NISAR − GNSS).
- **Predictions:**
  - CV: LOO RMS with a plane ≤ 2.0 cm/yr, on both passes (full-screen P1 was 3.78 for asc).
  - Houston: LOO RMS ≤ 0.6 cm/yr with the better of plane / quadratic, and ≤ the P1 RMS (constant offset only) of the same pass.
- **Gut:** CV improves but HANF stays an outlier: its error looks local (+5 to +8 in every variant). Houston passes.
- **Product meaning:** if it passes, the app shows "GNSS-anchored" maps where stations exist, and radar-only maps with a wider trust band elsewhere.
- **P-port deviation (02:20):** the Central Valley desc pair 2026-05-05 → 05-17 is an unusable product (0% valid pixels in the whole window), and it crashed the chain at the reference check. The pipeline now **skips** pairs it cannot reference and prints SKIP; the chain then has a gap there. No thresholds or settings changed. This applies to all sites from now on.
- **P-port Houston (2026-11-11 ~02:40): FAIL on P1, P2 and P4, both passes; P3 "passes" but is meaningless here.**

  | Pass | P1 RMS (limit 1.0) | P2 corr | Offset |
  |---|---|---|---|
  | asc | 11.9 cm/yr | −0.26 | −9 cm/yr |
  | desc | 95 cm/yr | −0.02 | −118 cm/yr |

  - Single-pair scene medians on desc run from −138 to +108 cm per 12 days.
  - **Diagnosis (post-hoc):** phase relative to station ADKS, in cycles, desc pairs:
    - most stations sit at 0 ± 0.3 cycles in every pair, so the bulk of the city is fine;
    - TXET jumps by −5 to +11 whole cycles from pair to pair; THSU / WEPD by −11 / −15 in the first pair.
    - So these are **unwrapping islands**: integer-cycle jumps in patches, inside what the product labels as one connected component. If the reference box (SE Pasadena, by the ship channel and bay) sits in or near such a patch, every pair gets a random multi-cycle offset. That gives the absurd scene medians and offsets.
  - **Implication:** the pipeline needs (1) a reference chosen for connectivity and stability (or the GNSS network), and (2) per-pixel detection of integer-cycle jumps along the chain. This is the clearest ML target found so far (see `notes/app-goal.md`, ML-2).
- **P-port, Central Valley both passes (full screen, pre-registered rules), 2026-11-11 ~02:00:**

  | Pass | P1 RMS | P2 corr |
  |---|---|---|
  | asc | 3.78 | 0.849 |
  | desc | 3.62 | 0.858 |

  - P1 and P2 fail on both passes. The desc chain skipped 2 unusable pairs (05-05, 05-17).
  - **P3 PASS: evening/morning noise ratio 0.81 (< 1.5).** In dry Central Valley the evening pass is *not* noisier. The Mexico "morning is cleaner" rule is about climate, not clocks, as predicted. Houston's 5.3 is contaminated by iono-screen jumps (below), so the humid-side half cannot be judged yet.
  - P4 fail: std(asc − desc) 5.44 cm/yr.
- **H9-P verdict (plane-only iono): FAIL as pre-registered, but a big gain in Houston.**

  | Site, pass | RMS full screen | RMS plane only | Corr full → plane |
  |---|---|---|---|
  | Houston asc | 11.86 | **5.65** | – |
  | Houston desc | 95.0 | **29.7** | – |
  | CV asc (development) | 3.78 | 2.78 | – |
  | CV desc (development) | 3.62 | 3.72 | 0.858 → 0.854 |

  - Houston: "plane ≤ full on each pass" holds; "P1 ≤ 1.0 on the morning pass" fails. With Mexico D10 failing earlier, the verdict is FAIL.
  - CV desc: no gain, unlike CV asc.
- **KEY DIAGNOSIS (post-hoc, Houston desc): the multi-cycle jumps are in the GUNW `ionospherePhaseScreen`, not in the interferometric phase.**
  - Relative to station ADKS, per pair:
    - phase alone: all stations within ±0.8 cycle (one TXET pair at 2.7);
    - iono screen alone: TXET −10.5 / −7.2 / +3.1 cycles, THSU +11.4, WEPD +15.0.
  - So the split-spectrum ionosphere estimate has its own unwrapping errors (sub-band phases are noisier, and the iono estimate amplifies them), and subtracting it *injects* islands of whole-cycle-scale error into a clean phase.
  - This unifies Mexico H8 (iono adds edge noise), CV (+10 cm winter bias) and Houston (islands): **in all three cities the dominant error source we found is the per-pair GUNW ionosphere screen.**
  - The product cure has to be a *robust* ionosphere: detect and mask the jump islands in the screen, then fit a smooth surface (or estimate the ionosphere per date over the network). That is an ML-2 target with clean labels: the screen minus a robust smooth fit, quantised.
- **P-tie result (2026-11-11 ~02:10): FAIL on every site and pass.** Run with `scripts/p_tie.py`.

  | Site, pass | Stations | P1 RMS (offset only) | LOO plane RMS | LOO quadratic RMS | Median abs residual (plane) |
  |---|---|---|---|---|---|
  | CV asc | 13 | 3.78 | 2.95 (limit 2.0) | – | 2.48 |
  | CV desc | 12 | 3.62 | 3.75 | – | 2.71 |
  | Houston asc | 26 | 11.86 | 10.13 | 13.18 | 2.71 |
  | Houston desc | 33 | 95.05 | 103.1 | 119.0 | 30.2 |

  - Gut was wrong for Houston: a smooth tie cannot help when the error is islands.
- **Post-hoc (labelled, not a test): P-tie on the plane-only iono variant.**

  | Site, pass | P1 RMS | LOO plane RMS | Median abs residual |
  |---|---|---|---|
  | CV asc | 2.78 | 2.61 | **1.07** |
  | CV desc | 3.72 | 4.72 | 3.14 |
  | Houston asc | 5.65 | 6.01 | 1.93 |
  | Houston desc | 29.7 | 24.0 | 12.1 |

- **Reading.**
  - The errors that remain are **not smooth**. A plane through the stations removes almost nothing, and the gap between RMS and median shows a few stations carry most of the error.
  - The GNSS tie is still the right *label* ("anchored" vs "radar-only"), but it is not a cure.
  - The cure has to come before the rate: robust per-pair iono screen (mask the jump islands), then the per-pixel cycle-jump check (ML-2).
  - CV desc is consistently worse than CV asc in every variant, so something pass-specific is going on there. Not diagnosed yet.

### L-retro. Did NISAR see the Chhochen Khola source slope move before the 26 Aug 2026 failure? (pre-registered 2026-11-11 ~02:30, before any download of the area)
- **Event.** Ice-rock avalanche, north face of Langtang Lirung (7,227 m), Rasuwa, Nepal. 26 Aug 2026, 08:37 local. ~67 Mm³; cascaded into a debris flow down the Lhende Khola.
  - The source is at 5,000–5,500 m (ICIMOD). No published coordinates found (the Springer paper is paywalled).
- **Not fully blind.** ICIMOD's FAQ cites a scientific coordination group: "satellite radar shows the slope had moved by up to 30 cm in the three to four weeks before it failed". It gives no product or method. So we know a precursor existed in *some* radar data; the question is whether *our* pipeline sees it in NISAR GUNW at 80 m.
- **Data.** NISAR asc p98 f16, bbox 85.37–85.67 E, 28.17–28.47 N.
  - 25 pairs, 2025-11-28 → 2026-10-06. All 12-day except 2026-07-26 → 08-19 (24 days: the 08-07 acquisition is missing).
  - Pre-event: 21 pairs, the last one being 07-26 → 08-19, which covers ~4 to ~1 weeks before failure. Co-event: 08-19 → 08-31.
  - DEM: Copernicus GLO-30.
- **Source zone (defined from post-event data, which is allowed because it uses location only).**
  - Pixels within 5 km of the summit, on the north side (DEM aspect within 90° of north), DEM height 4,500–6,000 m.
  - Must lose coherence in the co-event pair: γ(08-19 → 08-31) < 0.5 × median γ of the same pixel over pre-event pairs Dec–Jun.
  - Must be the largest connected such patch (≥ 20 pixels).
  - If no patch qualifies, the source zone is undefined, and the test is reported as "could not locate".
- **Control.**
  - Pixels at the same DEM height band (source median ± 300 m), north-facing, 3–15 km from the source zone, outside the co-event decorrelation mask.
  - Statistic per pair: d = median LOS(source) − median LOS(control), in cm.
  - The same-height control removes most of the stratified troposphere.
  - Ionosphere: plane-only screen (the H9-P lesson). LOS sign is the pipeline's convention, positive = toward the satellite.
- **Noise level σ_d.** Robust std (1.4826 × MAD) of d over the pairs with ref date ≥ 2025-12-10 and sec date ≤ 2026-06-08.
  - That is a winter–spring baseline. It excludes the first pair and the monsoon.
- **Predictions.**
  - **L1 (displacement precursor):** in the 07-26 → 08-19 pair, |d| ≥ 3 σ_d, with the sign of downslope motion. Downslope on a north face is mostly north + down; for the asc east-looking geometry "down" projects to *away* from the satellite, so the sign is negative.
  - **L2 (acceleration):** mean |d| per 12 days over the last 3 pre-event pairs (07-02 → 08-19, scaled by days) ≥ 2× the baseline mean |d| per 12 days.
  - **L3 (coherence precursor):** median coherence in the source zone divided by median in the control, in the 07-26 → 08-19 pair, is below the minimum of that ratio over the baseline pairs.
- **Gut.**
  - L1 fails: 30 cm in 3–4 weeks is > 2 L-band cycles per 24 days on a small, steep, icy patch, so it most likely decorrelates or unwraps wrong instead of showing as clean displacement.
  - L3 passes (decorrelation is the visible symptom), but it is confounded by monsoon snow and melt, so even a pass is weak evidence.
  - L2 is a coin toss.
- **Product meaning.**
  - If L1 or L3 passes, a "slope-watch" layer is plausible: flag steep high slopes whose coherence or LOS departs from same-height neighbours.
  - If all fail, NISAR GUNW at 80 m is too coarse for this failure type, and we say so in the app (only full-resolution RSLC or another product could do it).
- **Implementation details fixed before any data was seen** (`scripts/l_retro.py`, 02:45).
  - LOS pixels need coherence ≥ 0.3 and a nonzero component.
  - Source and control pixels must be in the same connected component (the one dominant among valid control pixels), so no 2π offsets between them.
  - d needs ≥ 5 valid source and ≥ 20 control pixels, else NaN. A NaN in the last pre-event pair means L1 FAILS ("not measurable").
  - Coherence ratio uses all pixels (no threshold).
  - Summit is snapped to the DEM maximum within 3 km of 28.2563 N, 85.5170 E.

### H10. Robust ionosphere screen: median-smooth the GUNW iono screen before subtracting (pre-registered 2026-11-11 ~02:55, before the variant has run anywhere)
- **Why.** The Houston diagnosis shows the per-pair `ionospherePhaseScreen` itself has multi-cycle-scale jump islands (its own unwrapping errors). The true ionosphere is smooth on scales ≫ 10 km.
  - Plane-only (H9-P) throws away the real curved part of the ionosphere.
  - A robust smooth keeps the curve but rejects islands smaller than ~5 km.
- **Variant `--iono-robust`** (`nisar.process(iono="robust")`):
  1. 2 km block medians of the screen (blocks need ≥ 25% valid pixels);
  2. 5 × 5-block median filter (10 km window, NaN-aware);
  3. fill gaps by nearest block;
  4. bilinear upsample back to the pixel grid;
  5. subtract.
  - Everything else is as in P-port (plane removed, same reference, same chains, skips).
- **Predictions (P-port metrics from `p_port.py`, variant `-ionorobust`).**
  - **R1:** Houston P1 RMS robust < plane-only on both passes (plane: asc 5.65, desc 29.7).
  - **R2:** Houston asc P1 RMS ≤ 1.0 cm/yr (the original G1 target).
  - **R3 (CV, development site, seen data):** P1 RMS asc ≤ 2.78 and desc ≤ 3.62 (the best earlier variant on each pass).
  - **R4 (sanity on Mexico):** at MMX1, asc rate within 2 cm/yr of GNSS (−26.7). Reported only if time allows.
- **Gut.**
  - R1 pass.
  - R2 fail: TXET also has a 2.7-cycle jump in the interferometric phase itself, and some stations sit on islands > 5 km.
  - R3: asc pass, desc coin toss.

### D-geom (a). Are the sloppy directions of a fault the multipole truncation? (pre-registered 2026-11-11 ~03:05, exploratory, synthetic only)
- **Setup.**
  - Same Okada → deplaned asc LOS map and Fisher matrix as D-ident (Nagqu grid, 1 cm white noise).
  - Normal fault 192/50/−119, 30 × 15 km.
  - Hold the moment fixed (Mw 6.95) and move the centroid depth over 10, 15, 20, 30, 40, 60, 80 km.
- **Theory.**
  - Far from the source (depth ≫ size) the surface field is a point moment tensor plus corrections of relative size (size/depth)².
  - At fixed moment the two "shape" directions (log L vs log W vs log slip, moment-preserving) enter only through those corrections, so their Fisher eigenvalues should carry an *extra* factor ~ (size/depth)⁴ relative to the rest.
- **Predictions.**
  - **G1:** over depth 30–80 km, the log–log slope of the two smallest eigenvalues vs depth is steeper than the median slope of the other seven by ≥ 3 (theory says ~4).
  - **G2:** at 80 km, the two sloppiest eigenvectors are moment-preserving shape combinations: |coefficient sum over logL, logW, logslip| < 0.3 × their norm, and ≥ 70% of the vector's weight on those three.
- **Known confound.** The window is fixed (Nagqu grid). Deep sources make broad signals that the window and the deplaning truncate, which can steepen *all* slopes. That is why G1 compares slopes between eigenvalues, not absolute slopes.
- **Gut.** G2 pass. G1 pass but with an exponent nearer 2–3 than 4, because of the window.
- **L-retro result (2026-11-11 ~03:25).**
  - **Source zone located, barely.** 20 pixels (the minimum), 0.13 km², centroid 85.532 E, 28.282 N, median height 5,270 m, ~3 km NNE of the summit (DEM max 7,207 m). That fits ICIMOD's 5,000–5,500 m.
  - **L1 FAIL (not measurable).** The source zone has coherence ≥ 0.3 in only 1 of 25 pairs (Feb). Every pre-event pair from mid-April on sits at ~0.09, the estimator floor, so there is no LOS displacement to measure. σ_d is undefined (1 baseline value).
  - **L2 FAIL (not measurable).**
  - **L3 PASS by the letter, but uninformative.**
    - Coherence ratio source/control: 0.17 in 07-26 → 08-19 vs the baseline minimum 0.47.
    - But the source sits at the coherence floor from mid-April, four months before failure. The ratio falls in June–August because the *control* recovers (0.29 → 0.67) while the source stays at the floor.
    - Likely confound: the source is glacier ice, and melting ice decorrelates every summer. A glacier-matched control (WorldCover snow/ice) would be needed to say more. Not done.
  - Gut was right on L1 and L2. L3 passed as guessed, and for the confounded reason.
- **Post-hoc (labelled): the debris-flow path is clearly mapped by coherence loss.**
  - Co-event pair 08-19 → 08-31 vs the adjacent pre-event pair: largest patch with coherence ratio < 0.35 (pre ≥ 0.4) is 4.8 km².
  - Straight-line extent 5.8 km, from ~3,990 m down to ~2,650 m along the valley, heading NW (85.50 E 28.30 N → 85.46 E 28.34 N). The quick-look shows a longer, curved tongue.
  - **Same rule on pre-event pairs** (07-14 → 07-26, 07-26 → 08-19): largest patch 0 and 18 px, so there is no false alarm in the weeks before the event.
  - **But 09-12 → 09-24 gives 644 px** (widespread autumn decorrelation, control coherence 0.25). The rule needs a season or scene-wide guard.
- **Reading for the app.**
  1. NISAR GUNW at 80 m cannot watch an icy 5,000 m face for precursors: no coherence.
  2. It *can* map where a debris flow went within one 12-day revisit. This is the known "damage proxy map" idea (coherence change, ARIA DPM). It is not novel, but it is a solid, honest demo layer.
  3. ICIMOD's "30 cm in 3–4 weeks" was probably measured by offset tracking or another sensor, not interferometric phase. NISAR GUNW carries a `pixelOffsets` layer; that is the next thing to try (direction L-offsets).
- **H10 result (2026-11-11 ~03:20): FAIL on every prediction. Robust (10 km median) screen is *worse* than plane-only.**

  | Site, pass | Full screen | Plane only | Robust (H10) | Prediction |
  |---|---|---|---|---|
  | Houston asc | 11.86 | 5.65 | 7.96 | R1 fail, R2 fail |
  | Houston desc | 95 | 29.7 | 87.2 | R1 fail |
  | CV asc | 3.78 | 2.78 | 2.90 | R3 fail (P2 corr 0.906 passes the P2 limit) |
  | CV desc | 3.62 | 3.72 | 3.82 | R3 fail |

  - **Reading.** The bad parts of the screen are not small islands that a 10 km median can reject. In Houston desc the station errors after the robust screen span −316 to +154 cm/yr, so the screen errors cover areas comparable to the window.
  - A per-pair spatial fix cannot clean them. The cure has to use **time**: the ionosphere estimated per *date* over the network (P-iono2, Fattahi 2017), or closure / jump detection along the chain (ML-2).
  - Plane-only stays the best simple option.
- **D-3 frozen results (written 2026-11-11 ~03:30, BEFORE opening any USGS mechanism or finite fault).**

  | Method | Strike | Dip | Rake | Mw | L × W, km | Depth, km | e0, n0, km | rms, cm (data std 16.0) | Time |
  |---|---|---|---|---|---|---|---|---|---|
  | Hybrid H1 (primary) | 77.6 | 29.6 | 110.8 | 7.28 | 107 × 16 | 5.0 | −8.1, 0.9 | 9.9 | 6.9 s |
  | Hybrid H2 | 78.0 | 23.7 | 105.5 | 7.27 | 120 × 17 | 4.7 | −3.5, 1.4 | 9.9 | 7.0 s |
  | SBI v3 (net + exact slip) | 76.9 [74, 80] | 35 [21, 50] | 135 [114, 157] | 7.15 [7.03, 7.27] | 77 × 18 | 7.0 | – | median 12.5 | 1.4 s |
  | LSQ, 128 starts | 241.6 | 63.6 | −149.3 | 7.66 | 61 × 27 | 24.5 | −8.2, −9.7 | 9.5 | 157 s |

  - SBI v3 rake mass: thrust 55%, right-lateral 42%.
  - All network-based answers: ENE-striking, shallow-dipping thrust, Mw 7.15–7.28, with L at or past the prior edge (80 km for v3; the hybrid went to 107–120 km).
  - LSQ instead found a steep oblique right-lateral fault (the conjugate-like geometry) with Mw 7.66.
  - Data coverage is 47% valid (tropical vegetation).
- **D-3 unblinded (2026-11-11 ~03:40): FAIL.**
  - USGS reference:
    - **Mww 7.7**, centroid depth 13.5 km;
    - W-phase planes 25/75/−155 and **288/66/−17**: a strike-slip earthquake on a WNW-striking plane;
    - finite fault 295/60/−16, three segments striking 275–295°, sub-fault lengths 40 + 30 + 30 km, Mw 7.66.

  | Method | Strike diff | Dip diff | Mw diff | Location | Verdict |
  |---|---|---|---|---|---|
  | **H1 (primary)** | 52° | 45° | 0.42 | pass | **FAIL**: thrust instead of strike-slip |
  | H2 | 53° | 51° | 0.43 | pass | FAIL |
  | SBI v3 | 52° | 40° | 0.55 | pass | FAIL |
  | LSQ | 47° | 3° | 0.04 | pass | fails strike only; right Mw, steep plane, wrong rake |

  - **Gut check:**
    - location passed as predicted;
    - Mw failed, and the gut said it would fail exactly in this case ("unless the rupture is clipped by the box");
    - strike/dip, called 50/50, failed;
    - v3 mechanism (60% predicted right) was wrong.
  - **G3 status: 1 of 3 blind events pass (Kumamoto). The G3 target (2 of 3) is NOT met.**
- **Post-hoc diagnosis (labelled):**
  1. **The rupture is bigger than both the box and the prior, which was risk (1) stated in advance.**
     - The finite fault is ~100 km long (three segments); Wells & Coppersmith gives ~140 km surface length for M7.7 strike-slip.
     - The window is ±55 km and the prior allows L ≤ 80 km.
     - The network then explains a clipped, edge-touching signal with a long, shallow thrust (the hybrid pushed L to 107–120 km and Mw stayed low). LSQ, with a free moment, got Mw right.
  2. Only 47% of the window is valid (vegetation), and the M6.6 aftershock 76 km W may leak in at the west edge.
  - **Product fix (not a test):** D-auto must size the box and the length prior from the catalogue magnitude (box half-size ≥ 0.75 × the scaling-law length; L prior up to 1.5× it). The training set must contain clipped ruptures, or the window must be large enough.
  - Any rerun on Panama is now unblinded and can only be a diagnostic.

### L-offsets. Can GUNW pixel offsets see the Chhochen Khola source move before failure? (pre-registered after L-retro, before any offset data is read)
- **Why.** Interferometric phase fails on the icy face (no coherence). Offsets (amplitude cross-correlation) do not need phase coherence.
  - The GUNW has `pixelOffsets/HH/{slantRangeOffset, alongTrackOffset, correlationSurfacePeak}` on the same 80 m grid.
- **Data.** Same 25 asc pairs and bbox as L-retro; offsets window read into a separate cache (`data/raw/nisar_offsets/langtang/`).
- **Zones.** Exactly as L-retro: the same source zone (20 px) and the same control rule.
  - Secondary zone, reported but not judged: a 1 km disc around the source centroid, north-facing, 4,500–6,000 m.
- **Statistic per pair.** d_rg, d_az = median offset (source) − median offset (control), in metres, using pixels with mask OK and correlation peak above the per-pair median of the control. Combined: d = hypot(d_rg, d_az).
- **Noise.** Robust std (1.4826 × MAD) of each component over the same baseline pairs as L-retro (ref ≥ 2025-12-10, sec ≤ 2026-06-08).
- **Predictions.**
  - **O1:** in 07-26 → 08-19, |d_rg| ≥ 3 σ_rg or |d_az| ≥ 3 σ_az.
  - **O2:** d in 07-26 → 08-19 ≥ 3 × the baseline median of d.
- **Gut.** Both fail.
  - Offset noise at 80 m posting is likely tens of cm.
  - 30 cm is a *peak* value spread over the 3–4 weeks.
  - The source zone is tiny (20 px).
  - Chance of a pass about 0.2.
- **L-offsets implementation note (before reading any offset values).**
  - "Mask OK" means both sub-swath digits are nonzero and the water digit is 0. The GUNW mask is 3 digits: water / reference sub-swath / secondary sub-swath; this window has values 11 and 111.
  - The offsets grid is the same 80 m grid shifted by 4 rows, so layers are placed onto the zone grid by coordinates.

### D-3w. Panama diagnostic: was the window the problem? (pre-registered after unblinding, so a DIAGNOSTIC, not a blind test)
- **Hypothesis.** D-3 failed because the ~100 km rupture did not fit the ±55 km window or the L ≤ 80 km prior.
- **Change only the size settings; everything else stays frozen.**
  - Window ±1.0° (fetch_event `--half 1.0`), 2 km blocks (`SBI_BLOCK=24`), so the grid stays ~110 × 110.
  - Length prior up to 200 km (`SBI_LMAX`); centroid prior ±60 km (`SBI_XYMAX`).
  - Hybrid bounds ±80 km and L ≤ 300 km.
  - Same network, steps, hybrid H1 / H2.
- **Prediction (judged against the same Mww planes and rules as D-3).** H1 gets strike and dip within 15° of the closer plane AND |ΔMw| ≤ 0.2.
  - If it does, the window explanation holds and D-auto must size windows by magnitude.
  - If the mechanism is still thrust, the problem is elsewhere: vegetation gaps (47% valid), aftershock leakage, or the prior.
- **Gut.** Mw improves (≥ 7.45). The mechanism is right with probability ~0.5.
- **L-offsets result: O1 and O2 FAIL (not measurable under the pre-registered rule).**
  - The source zone's correlation peak is 0.22–0.36 against 0.77 in the control, so the "peak ≥ control median" filter removed every source pixel in every pair.
  - The along-track offset layer is NaN throughout this window; only slant range is populated.
- **Post-hoc (labelled; no quality filter, mask = 11 only, slant range only).** Source minus control, metres, + means range increased (moving away from the satellite).
  - Winter (Dec–Feb): +0.12 to +0.26 per 12 days, every pair.
  - Spring: −0.08 to +0.41.
  - 07-14 → 07-26: +0.28.
  - **07-26 → 08-19 (24 days): +0.37**, i.e. +0.19 per 12 days.
  - Co-event: +0.18.
  - Post-event: +0.07, +0.64, +0.52.
  - Baseline robust std 0.12 m per pair. The 1 km disc gives the same picture.
- **Reading.**
  - The source ice was moving **steadily** away from the satellite at roughly 0.15 m per 12 days (~4–5 m/yr in range) all winter: a hanging glacier flowing.
  - The last pre-event pair (0.37 m in 24 days) matches ICIMOD's "up to 30 cm in 3–4 weeks" in size, **but it is not faster than the winter rate.** In our data the "precursor" is indistinguishable from normal glacier flow.
  - The control pixel scatter is only 2–7 cm, so offsets are precise on rock. On the icy source the correlation is low, so per-pair noise is ~12 cm.
  - So we cannot claim acceleration. Being able to measure motion on ice where phase fails is still a real capability.
- **Product meaning.** A "slope-watch" layer from GUNW offsets is feasible as a *velocity* map of fast ice and rock (decimetres per 12 days).
  - An *acceleration alarm* would need a longer, cleaner series or a smaller noise floor (both passes, more pixels).
  - Next steps: L-watch over all of High Mountain Asia's steep ice faces, and comparison with the asc+desc velocity of known fast glaciers.
- **D-geom (a) result: G1 FAIL (2.48 < 3), G2 FAIL.**

  Fisher eigenvalues, 1 cm noise:

  | Centroid depth | Largest eigenvalue | Smallest eigenvalue |
  |---|---|---|
  | 10 km | 4.8e6 | 1.7e3 |
  | 80 km | 3.1e3 | 7.2e−4 |

  Log–log slopes vs depth (30–80 km), stiffest → sloppiest: −4.6, −5.6, −5.3, −5.5 | −7.9, −9.3, −8.5, −8.1, −9.6.

  The 80 km sloppiest eigenvectors:
  - #8 is mostly north position plus depth (shape weight 0.09);
  - #9 is east position plus slip/width (shape weight 0.50).

  Gut said G2 pass and G1 pass with a weaker exponent; both were wrong as stated.
- **Post-hoc (labelled), but this is the interesting part.**
  - **The spectrum splits into two bands, 4 + 5, not 7 + 2.** Four directions decay like depth^−5.2 and five like depth^−8.7.
  - The gap of ~3.5 is close to the predicted extra (size/depth)⁴. So the multipole hierarchy is there, but the second band holds 5 directions, and it includes **horizontal position**, not only the moment-preserving shape.
  - **Why position (hypothesis):** we deplane every map (remove a best-fit plane). For a source much deeper than the window is wide, moving it sideways changes the map by ≈ its gradient, which over the window is ≈ a plane. Deplaning projects that onto the ramp space and kills it.
  - So deplaning is a *quotient by a gauge group* (planar ramps), and position becomes a near-gauge direction at large depth.
  - The pre-registered count was wrong because it ignored the gauge, not the multipoles.

### D-geom (b). Does the ramp-removal degree move directions between bands? (pre-registered after D-geom (a), synthetic)
- **Setup.** Same Fisher setup as D-geom (a), asc look, depths 30, 40, 60, 80 km. The deplaning projector removes polynomials of degree k:
  - k = 0: constant only;
  - k = 1: plane (= D-geom (a));
  - k = 2: quadratic.
  - "Steep band" = eigen-directions whose log–log slope (30–80 km) is steeper than −7.
- **Predictions.**
  - **Q1:** k = 0 has fewer steep-band directions than k = 1 (5), and at 80 km the two sloppiest eigenvectors carry < 30% of their weight on (e0, n0).
  - **Q2:** k = 2 has more steep-band directions than k = 1.
- **Gut.** Q1 pass, Q2 pass.
- **D-geom (b) result: Q1 FAIL, Q2 PASS.**

  | Ramp removal | Steep-band count | (e0, n0) weight in the 2 sloppiest at 80 km |
  |---|---|---|
  | k = 0 (constant) | 5 | 0.33 (limit < 0.30) |
  | k = 1 (plane) | 5 | 0.61 |
  | k = 2 (quadratic) | 6 | 0.29 |

  - The whole spectrum steepens by ~1.2 in exponent per polynomial degree, rather than single directions hopping between bands.
  - Even with only a constant removed, horizontal position stays in the sloppy band (eigenvectors 7 and 9 at 80 km are n0- and e0-dominated).
- **Reading.** The deplaning-gauge story from D-geom (a) is at most partial.
  - The better description: the window is small compared with the depth, so what the map holds is a *Taylor expansion* of the far field across the window (u0 + x·∇u + ½ x·∇²u·x + …).
  - Each polynomial removed deletes one Taylor order. That shifts every eigenvalue's exponent, more or less uniformly.
  - Position is sloppy because moving a broad field sideways mostly changes its *low* Taylor orders, the ones we remove or that are weak.
  - The two-band structure (4 + 5, gap ~3.5) survives in all three cases.

### C. Detection threshold: what subsidence bowl can a NISAR rate map of Mexico City see? (pre-registered 2026-11-11, before running)
- **Data.** Mexico City asc chain (cached), plane-only iono screen (the best simple variant tonight), chain-sum rate. The same robust plane is removed *after* injection, as in the product.
- **Injection.** Gaussian subsidence bowls with peak rate v ∈ {0.5, 1, 2, 4, 8} cm/yr and FWHM s ∈ {0.5, 1, 2, 4} km.
  - Because the chain-sum rate is linear in the pair displacements, injecting v·Δt into every pair equals adding the bowl to the rate on valid pixels. That is what is done.
  - Centres: 40 per (v, s), drawn from "quiet" pixels (|rate smoothed at 2 km| < 2 cm/yr, valid), ≥ 3 s apart, seed 0.
- **Detector.** D = (rate smoothed with a Gaussian of FWHM s) at the centre, minus the median rate in the annulus 2s–4s.
  - Detected if D < −3 σ_D. σ_D = robust std of D at 200 other quiet points with nothing injected, same s.
- **Predictions.**
  - **C1:** detection ≥ 90% for v ≥ 4 cm/yr at every s ≥ 1 km.
  - **C2:** the 50% detection rate (interpolated in log v) is higher for s = 0.5 km than for s = 2 km.
  - **C3:** false-alarm rate (D < −3 σ_D at the 200 noise points) ≤ 5% for every s. A check of the robust σ.
- **Gut.**
  - 50% thresholds ≈ 1–2 cm/yr.
  - C1 pass. C2 pass. C3 pass.
  - 4 km bowls may be *harder* than 2 km: atmosphere and long-wavelength iono error grow with scale, and the plane removal eats part of a wide bowl.
- **C deviation (technical, before any s ≥ 2 km result existed).** The first run required noise points > 4s from every injected centre. At s = 2 and 4 km no such points exist in the 45 × 52 km map, so σ_D was NaN.
  - The noise statistic is computed on the *un-injected* map, so that requirement was unnecessary. Noise points are now 200 random quiet pixels.
  - First-run values, kept for the record (s = 0.5 / 1 km): σ_D 0.95 / 0.91 cm/yr, false alarms 1.5 / 2.5%, 50% thresholds 5.4 / 5.2 cm/yr.
- **C result (2026-11-11): C1 FAIL, C2 FAIL, C3 PASS.**

  | Bowl FWHM | Bowls | σ_D, cm/yr | False alarms | Detection at 2 / 4 / 8 cm/yr | 50% threshold, cm/yr |
  |---|---|---|---|---|---|
  | 0.5 km | 40 | 0.90 | 2.0% | 2 / 20 / 95% | 5.3 |
  | 1 km | 40 | 0.96 | 1.5% | 5 / 15 / 75% | 6.0 |
  | 2 km | 40 | 1.08 | 0% | 2 / 5 / 38% | > 8 |
  | 4 km | 13 | 1.23 | 0% | 0 / 0 / 15% | > 8 |

  - C2 came out the *other* way: small bowls are easier. The gut's side note ("4 km may be harder") was right.
- **Reading (post-hoc, labelled).**
  1. **The noise floor is real and modest:** σ_D ≈ 0.9 cm/yr at 0.5–1 km, rising to 1.2 at 4 km. Long-wavelength atmosphere and iono residuals grow with scale, as in the Kolmogorov picture (M3).
  2. **The registered detector wastes half the signal.** Smoothing a Gaussian bowl with a same-width Gaussian halves its peak.
     - With 3 s spacing, neighbouring bowls also sit in the 2s–4s annulus and lift the background. Plane removal eats part of the wide bowls.
     - So "≈ 6 σ_D" is the detector's threshold, not the data's.
     - A proper matched-filter amplitude estimate (fit the bowl shape, compare with its own noise) should reach ≈ 3 σ_D ≈ 3 cm/yr at 1 km. Not tested; it would need a new pre-registration.
  3. **Product meaning.**
     - The app can honestly say "bowls under ~3–5 cm/yr and wider than ~2 km are below what this map can confirm on its own".
     - Wide, slow regional subsidence needs GNSS or a longer series. Sharp local bowls (wells, construction, sinkholes) are where NISAR alone shines.
     - That is the ML-3 target, stated in numbers.

### C′. Same injection test with a matched-filter detector (pre-registered after C, before running)
- **Change: only the detector.** For each centre:
  - subtract the median of the 2s–4s annulus;
  - fit the amplitude a of the bowl template g (peak 1, FWHM s) by least squares on the disc r ≤ 2s.
- **Rule.** Detected if a < −3 σ_a. σ_a = robust std of a at the same 200 un-injected quiet points. Same injections, spacing, seed and plane removal as C.
- **Predictions.**
  - **C′1:** 50% threshold at s = 1 km ≤ 3.5 cm/yr.
  - **C′2:** thresholds do not decrease with size (0.5 ≤ 1 ≤ 2 ≤ 4 km, ties allowed).
  - **C′3:** false alarms ≤ 5% at every s.
- **Gut.** All three pass.
- **C′ result: C′1 FAIL, C′2 PASS, C′3 PASS.**
  - Amplitude noise σ_a: 1.77 / 1.92 / 2.12 / 2.45 cm/yr. Detection rates and 50% thresholds (5.3 / 6.0 / > 8 / > 8) are **identical** to C.
- **Correction of my C reading (point 2 is retracted).**
  - The smoothed-centre value and the least-squares amplitude are both Gaussian-weighted linear filters of the same data, so they have the same signal-to-noise. "Half the signal wasted" was wrong: the matched estimate doubles the amplitude *and* its noise.
  - **The ~5–6 cm/yr threshold for 0.5–1 km bowls is set by the data (one year, 28 pairs, morning pass), not by the detector.**
  - What still depresses detection for 2–4 km bowls is the test design: crowded neighbours in the annulus, and plane removal on wide bowls. That part of point 2 stands.
- **Product numbers to show:** one year of NISAR asc over Mexico City confirms, on its own, a 1 km bowl only if it sinks faster than ~6 cm/yr at 50% (~8 cm/yr at ~80%).
  - That is fine for Mexico City's 20–40 cm/yr, but marginal for most cities (1–3 cm/yr). There, the honest app must say "radar-only: not detectable; needs longer series or GNSS".
  - With N years the threshold should fall like 1/√N for node noise and like 1/N^(1/2) or slower for edge noise (M2). To be checked when data allow.
- **Development check (CV only, not a test, ~03:10): "skip the per-pair iono screen, remove a smooth surface from the rate" (from M2: per-date iono telescopes in the chain sum).**

  | Variant | Plane P1 RMS | Quadratic P1 RMS | Cubic P1 RMS |
  |---|---|---|---|
  | No iono | 6.26 | 7.39 | 6.92 |
  | Full screen | 3.74 | 6.71 | 6.87 |
  | Plane-only screen | 3.14 | 6.17 | 5.68 |

  (Robust whole-map fits.)
  - Higher-degree surfaces absorb the real ~100 km subsidence bowl, so the idea fails where the signal is as wide as the window.
  - The no-iono error is not a smooth end-date screen either (plane 6.26 vs 3.1–3.7 with a screen).
  - **Not taken to Houston.**
  - Side note: the full-screen plane fit here (3.74) differs slightly from P-port (3.78) because the fit is different. P-port numbers stay the reference.
- **P-iono2 feasibility note.** ASF lists only sequential 12-day GUNW pairs (plus a few 24-day pairs where an acquisition is missing). Without redundant pairs, per-date ionosphere cannot be separated from per-pair estimate errors.
  - So the proper cure needs either redundant pairs (24/36-day GUNW, if NISAR ever releases them) or our own processing from RSLC.
  - Likely out of scope for the hackathon. Worth a question to the mentors on Friday: are longer-baseline GUNW pairs planned?
- **D-3w result (diagnostic, unblinded): FAIL by the D-3 rule, but the mechanism moved toward the truth.**

  | Method | Strike / dip / rake | Mw | Verdict vs Mww plane 25/75/−155 |
  |---|---|---|---|
  | H1 | 37 / 5 / 180 | 7.46 | strike PASS (12°), dip FAIL (70°), Mw FAIL (0.24), location PASS; 7.5 s |
  | H2 | 38 / 12 / −172 | 7.43 | same pattern |
  | SBI v3 alone | 86 / 27 / −57 | 7.23 | fail |
  | LSQ | 238 / 50 / −148 | 7.88 | Mw pass, strike off 50° |

  - Gut: "Mw improves (≥ 7.45)" right, just (7.46); "mechanism right p ≈ 0.5": half (strike-slip and strike right, dip absurd).
- **Reading.**
  - Widening the window turned thrust into strike-slip, put the strike within 12° of a Mww plane, and raised Mw from 7.28 to 7.46. So **the window size was a real part of the D-3 failure.**
  - The dip is not constrained: 5° is unphysical for this event. Only 33% of the wider window is valid (vegetation), and there is a single asc look.
  - D-ident says two looks separate mechanisms 9–11× better, and Panama has desc co-event pairs (p141, p40).
  - **Next for G3:** a two-look (asc + desc) SBI and hybrid, plus a magnitude-scaled window. Recorded as D-2look.

### D-2look-lite. Panama with asc + desc in the hybrid step (diagnostic, unblinded; pre-registered before running)
- **Method.**
  - Proposals: 2,000 samples from the panama-wide asc network (unchanged). Slip is solved exactly on the *concatenated* asc + desc deplaned LOS, each look weighted 1/√(its pixel count), so both looks count equally.
  - Selection and polishing: top 16, then 60 Adam steps on the equally weighted sum of both looks' mean squared residuals.
  - Same bounds as D-3w. H1 primary; H2 adds 16 prior starts.
- **Desc data.** p40 co-event pair 2026-10-02 → 10-14 (4000 mode, same window ±1°, 2 km blocks), heading −169°, left-looking.
- **Judged with the D-3 rule** against the Mww planes.
- **Predictions.**
  - **K1:** H1 dip within 15° of the plane it matches in strike (the D-3w failure point).
  - **K2:** H1 passes all four D-3 criteria.
- **Gut.** K1 pass (p ≈ 0.6). K2 fail (p ≈ 0.6), because Mw tends to stay low.
- **D-2look-lite result (diagnostic, unblinded): K1 PASS, K2 PASS as worded.**
  - Desc p40 co-event map: 22% valid, incidence 36°.

  | Arm | Strike / dip / rake | Mw | vs Mww plane 288/66/−17 | Time |
  |---|---|---|---|---|
  | H1 | 290 / 55 / 36 | 7.57 | strike 2°, dip 11°, ΔMw 0.13, location 0 km: all four pass | 13.7 s |
  | H2 | 288 / 49 / 44 | 7.38 | dip 17° (fail), ΔMw 0.32 (fail) | 13.7 s |

  - Against the USGS **finite fault (295/60/−16)**: H1 is 5° off in strike and 5° in dip.
  - The rms values printed are on the equal-weight scale, not cm.
- **Honest caveats.**
  1. D-3's success list also had "hybrid time < 10 s". Two looks doubled the forward cost, so **13.7 s fails that criterion.** K2 named only the four judge criteria. Engineering fix: fewer polish steps or a single fused forward call.
  2. **Rake is off by ~50°** (36 vs −17: left-lateral with thrust instead of a slight normal component). The mechanism family (left-lateral strike-slip on the WNW plane) is right; the dip-slip share is not. Rake is not in the judged criteria.
  3. **Unblinded.** The two-look idea came from D-ident before this run, but we had already seen USGS for Panama. This is strong evidence for the method, not a blind pass.
- **Meaning for G3 / the app.**
  - **Window scaled by magnitude + both passes + network proposals + 60-step polish** recovers a M7.7 strike-slip fault to within a few degrees, in seconds.
  - That is the pipeline to freeze and test blind on the next M6.5+ NISAR-era event (D-auto).
  - The D-ident prediction (two looks separate mechanisms 9–11× better) is borne out on real data.

### D-2look-K. No-harm check on Kumamoto (diagnostic, unblinded; pre-registered before running)
- **Setup.** Same `sbi_hybrid2.py` with the frozen Kumamoto asc network (default 1 km blocks, default priors). The desc co-event pair is p3 2026-07-08 → 08-13 (36 days), heading −168°.
- **Prediction K3.** H1 still passes all four judge criteria (strike, dip, Mw, location).
- **Gut.** Pass. The 36-day desc pair may add post-seismic signal and noise.
- **D-2look-K: could not run.** The Kumamoto desc p3 window has 0% valid 1 km blocks (the frame misses the epicentre window, or everything is masked). K3 is untested. Nagqu has no desc co-event pair. So D-2look evidence rests on Panama alone.
- **Quantised-jump idea for the iono screen (checked ~03:10, negative).**
  - The GUNW metadata gives f0 = 1239 MHz and a 40 MHz range bandwidth. The sub-band bandwidths are NaN (not recorded).
  - Assuming the usual split into thirds (f_L,H = f0 ∓ B/3), a 2π unwrapping error gives an iono-screen jump of:
    - 23.5 main-phase cycles (low band);
    - 23.0 cycles (high band);
    - 0.5 cycles (both bands).
  - The Houston jumps (−10.5, −7.2, +3.1, +11.4, +15 cycles vs ADKS) are not near these quanta. The screen is also low-pass filtered, which smears any quantum.
  - **No simple integer correction is possible.** The fix has to be detect-and-mask (ML-2) or redundant pairs.
- **Look at the Houston desc iono screens (post-hoc picture, `data/processed/houston_desc_ionoscreens.png`).**
  - The screen minus its plane is flat (std 0.0–0.1 cycles) in 12 of the first 18 pairs.
  - The other 6 (Oct–Jan) carry compact islands of 10–20 km, up to ±15 cycles, near the SE edge of the frame (bay and water side), plus range-direction streaks.
  - **This explains why plane-only still left 30 cm/yr:** a least-squares plane is dragged by a 15-cycle island, so in the island pairs we subtracted a bogus tilt.

### H11. Robust plane of the iono screen (pre-registered 2026-11-11 ~03:15, before running)
- **Variant `--iono-rplane`.**
  - Fit a plane to the screen, drop pixels whose residual exceeds max(3 × 1.4826 × MAD, 1 cycle); refit four times.
  - Subtract only that plane.
  - Everything else as P-port.
- **Predictions.**
  - **R1′:** Houston P1 RMS below plane-only on both passes (asc < 5.65, desc < 29.7).
  - **R2′:** Houston asc P1 ≤ 1.0 cm/yr.
  - **R3′:** CV asc P1 ≤ 2.78 (no worse than plane-only; CV screens have no big islands).
- **Gut.**
  - R1′ pass, desc especially.
  - R2′ fail but close (1–3).
  - R3′ pass (tie within 0.1).
- **H11 result: R1′ PASS (barely), R2′ FAIL, R3′ PASS (tie).**

  | Site, pass | Robust plane (H11) | Plane only | Note |
  |---|---|---|---|
  | Houston asc | 5.63 | 5.65 | – |
  | Houston desc | 27.8 | 29.7 | offset −18.6 |
  | CV asc | 2.78 | 2.78 | P2 0.911 passes |
  | CV desc | 3.73 | 3.72 | – |

  - CV P3 is now 1.09: still a pass.
  - Gut: R1′ pass but not "desc especially"; R2′ fail but not "close" (5.6 ≫ 1–3).
  - **The island-tilted plane explains only a small part.** Houston's remaining error lives elsewhere (next: per-pair station diagnosis).
- **Houston desc per-pair diagnosis (post-hoc, robust-plane iono).**
  - Station spread per pair is 1–3 cm in Oct–Jan and in a few spring pairs. In most pairs from 01-17 on, **one station, TXET, swings by ±35 to ±146 cm**, with alternating signs.
  - TXET's chain sum is only +12 cm, so these are **per-date (node) errors that telescope**. Other stations' chain sums are within ±10 cm.
  - Entries over 0.6 cycle are not integer-like: median fractional part 0.29.
- **KEY FINDING (post-hoc, labelled): gaps in the chain break the telescoping.**
  - `stack_rate` sums each pixel's displacements over the pairs where it is valid. If a pixel misses even one pair, its chain splits into segments. Each segment contributes the node errors of *its* end dates, which no longer cancel.
  - 71% of Houston desc rate pixels have gaps (median 2 segments). With ±1 m local date errors, a single gap ruins the rate.
  - Restricting to pixels valid in every pair:

  | Site, pass, iono | All pixels (P-port) | No-gap pixels | ≤ 1 segment |
  |---|---|---|---|
  | **Houston desc, rplane** | 27.84 (33 stations) | **2.80 (29 stations, median 2.30)** | 3.22 |
  | Houston asc, rplane | 5.63 (median 0.98) | 7.77 (median 1.08) | – |
  | CV asc, plane | 2.78 | 3.03 (11 stations) | – |

  - In Houston asc and CV the problem is a few outlier stations, not gaps.
- **Product consequences.**
  1. The rate estimator must be **gap-aware**: either require a complete chain, or carry a per-pixel segment count into the trust layer. For a path graph with s segments, node noise contributes 2s·σ²_node / T² to the rate variance (M2 extension, see math.md M2b).
  2. The trust audit's Q estimator uses nansum and is blind to this. It must add the segment term.
  3. The Houston map is usable on complete-chain pixels (~29% of pixels): 2.8 cm/yr against GNSS, still above the 1.0 target. Houston desc is the evening pass in a humid city, so node noise is large.

### H12. Gap-aware rate: keep only pixels valid in every pair of the chain (pre-registered 2026-11-11 ~03:20, on data not yet seen with this rule)
- **Rule.** Rate from `stack_rate` as before, but only on pixels whose chain has no gaps (valid fraction = 1). Plane removed as in P-port.
- **Test data, not yet seen with this rule.**
  - (a) CV desc, plane-only iono: P1 RMS over the stations that keep pixels; all-pixel reference 3.72.
  - (b) Mexico City desc (evening pass, the humid-side case), plane-only iono: MMX1 error |NISAR − GNSS|, GNSS −26.7. The all-pixel plane-only value was −28.9, error 2.2. If MMX1 loses its pixels, (b) is "not measurable".
- **Predictions.**
  - **G-a:** CV desc RMS with no-gap pixels < 3.72.
  - **G-b:** Mexico desc MMX1 error with no-gap pixels < 2.2 cm/yr.
- **Gut.**
  - G-a fail: CV asc did not improve, and CV's problem is outlier stations.
  - G-b pass: the evening humid pass is node-noise dominated, like Houston desc.
- **Houston desc is excluded:** it is where the rule was found.
- **H12 result.**

  | Test | All pixels | No-gap pixels | Verdict |
  |---|---|---|---|
  | G-a: CV desc, P1 RMS | 3.72 (12 stations) | 3.04 (10 stations; 59% of pixels have gaps) | **PASS**, but 2 stations dropped out |
  | G-b: Mexico desc, MMX1 | −27.3 (error 0.59) | −25.0 (error 1.68) | **PASS by the letter (< 2.2)** |

  - **G-b is a registration mistake on my side.** The reference "2.2" came from the H9-P code path; the same-code all-pixel value is 0.59. Like for like, the no-gap rule made MMX1 *worse*: 56% of pixels have gaps, and the nearest complete-chain pixels sit elsewhere on the bowl's steep slope. So in substance G-b fails.
  - Gut: G-a fail was wrong; G-b pass was right by the letter, wrong in substance.
- **Overall reading.**
  - The gap rule rescues Houston desc (post-hoc, 27.8 → 2.8) and helps CV desc (3.7 → 3.0), but it costs coverage (56–71% of pixels are gappy).
  - At a single station on a steep bowl it can hurt, because it moves the sample to other pixels.
  - So the product should *not* simply drop gappy pixels. It should keep them and **carry the segment count into σ** (math M2b), plus prefer complete-chain pixels when sampling near a station.
  - Next engineering step: segment-aware σ in `trust.py`.

### P5b. Segment-aware σ: does adding the gap term (M2b) calibrate the trust layer? (pre-registered ~03:25, before running)
- **Change in `trust.py`.** σ² gets an extra term: the sum of the node variances (G1′ estimates, negatives clipped to 0) of every interior date where a pixel's chain is broken (segment starts and ends, excluding the chain's first and last dates), divided by T².
  - Opt-in flag `segments=True`; the default stays as before.
- **Runs.** `trust_audit.py --iono plane` for CV asc, CV desc, Houston asc, Houston desc. Same GNSS z-table as P5.
- **Prediction P5b.** |z| < 1 for 55–85% of stations in at least 3 of the 4 site/pass combinations.
- **Gut.** Fail.
  - CV's σ collapses for a separate reason: the Q estimator ≤ 0 clips to the floor.
  - Houston asc is outlier-dominated.
  - Houston desc may move into range.
- **Result P5b (~03:45): FAIL.** Only 1 of 4 combinations falls in the 55–85% band (3 were needed).

  | site/pass | |z|<1 without gap term | with gap term | median σ without → with (cm/yr) |
  |---|---|---|---|
  | CV asc | 46% | 54% | 0.11 → 3.71 |
  | CV desc | 83% | 92% | 13.6 → 14.1 |
  | Houston asc | 81% | 81% | 6.7 → 7.1 |
  | Houston desc | 9% | 21% | 0.10 → 17.1 |

  - Without the gap term, 2 of 4 were in range. With it, CV desc overshoots to 92%, so only 1 of 4 is.
- **Reading (post-hoc).**
  - The gap term does what M2b says. Where the old σ sat on the floor (CV asc, Houston desc, median about 0.1 cm/yr), it lifts σ to a realistic level and moves |z|<1 the right way: +8 and +12 points.
  - It cannot rescue Houston desc. The big errors there (tens to 150 cm/yr) are per-date jumps at stations that are *inside* complete segments too. The gap term only counts the breaks.
  - CV desc was already over-conservative. Its σ is dominated by the per-pair term, about 14 cm/yr against real errors of about 3.7.
  - So the trust layer has two separate faults:
    - (a) the Q estimator clips to the floor where the chain has few complete pixels;
    - (b) the per-pair variance is inflated in CV desc.
  - The gap term fixes neither directly. It is a correct but small add-on.
- **Decision.**
  - Keep `segments=True` available. It is physically right and moves things the right way.
  - Do not call the trust layer calibrated.
  - For the app: show the segment count as its own layer ("this pixel's chain is broken N times") rather than folding it into one σ. The user can then see why a pixel is uncertain.

### X-coh. Is field coherence a crop calendar? (pre-registered ~03:35, before looking at any coherence time series)
- **Data**
  - Central Valley asc sequential chain: the same pairs as the rate run.
  - Each pair's coherence is regridded to the `central-valley-asc/nisar_rate.tif` grid.
  - WorldCover fractions in 2 km blocks, the same file as H9.
  - Cropland blocks: cropland fraction > 0.6. Built-up blocks: built-up fraction > 0.3.
- **Quantity.** R_k = median block coherence over cropland blocks ÷ median over built-up blocks, for pair k. Built-up is the stable reference, so R removes whole-scene effects.
- **X1 (season).** Median R over pairs whose midpoint falls in May–Sep is lower than median R over pairs with midpoint in Nov–Mar.
  - Expected mechanism: summer irrigated canopy and harvest change the fields.
- **X2 (a calendar, not noise).** Lag-1 autocorrelation of R_k along the chain ≥ 0.3.
- **X3 (rain, secondary).** Mean R in wet pairs (≥ 2 mm within the 4 days ending at either date, the same rule as H9) is at least 0.05 lower than in dry pairs.
- **Gut.** X2 pass, X3 pass, X1 coin-flip.
  - Winter rain and winter wheat may make winter the worse season.
- **Post-hoc only.** k-means (k = 4) on per-block cropland coherence curves, to see whether fields split into "dances", for app layer ML-1. Labelled exploratory.
- **Result X-coh (~03:45): X1 PASS, X2 PASS, X3 FAIL.**
  - 25 pairs; 820 cropland blocks and 54 built-up blocks.
  - **X1.** Median R is 0.68 in summer against 0.87 in winter.
  - **X2.** Lag-1 autocorrelation is +0.85, so the curve is smooth: a calendar, not noise.
  - **X3.** Dry minus wet is only +0.015. Rain does not visibly cut field coherence at L-band.
- **What the curve looks like.**
  - R ≈ 0.85–0.95 from December to February.
  - It falls in steps through March and early April, reaching 0.61–0.62 from mid-April to May.
  - It recovers partly to about 0.73 in July–September.
  - The lowest points (0.56 and 0.64) are the two **24-day** pairs. That is longer time, not season, so read the August and October dips with that in mind.
  - Built-up coherence stays at 0.65–0.78 all year.
- **Post-hoc (exploratory) k-means, k = 4.**
  - Three clusters (722 blocks) follow the main calendar: best in January, worst in spring or late summer.
  - Cluster 2 (98 blocks) dances the opposite way: worst on 12 March, best on 20 September.
  - Possible reading: winter-grown fields (grain or forage, active Feb–Apr, bare after harvest) versus irrigated summer and orchard land. Not checked against a crop map (USDA CDL would be the check).
- **For the app.** This is a cheap, robust "rhythm of the fields" layer straight from the GUNW coherence: no unwrapping, no ionosphere. It answers the challenge's "agricultural lands transform" line directly.
  - Next check (to pre-register): do the clusters match USDA Cropland Data Layer classes?
- Files: `data/processed/central-valley-asc/x_coh.json`, `x_coh.png`.

### D-2look-speed. Freeze the two-look quake recipe and choose a speed-up on synthetics only (pre-registered ~03:33, before running)
- **Frozen recipe (Q-2look v1).**
  - The window and priors scale with the catalogue magnitude M, which is known at alert time. RLD = 10^(−2.44 + 0.59·M) km (Wells & Coppersmith subsurface rupture length).
  - Settings:
    - L_MAX = clamp(2·RLD, 80, 300) km;
    - XY_MAX (net) = max(25, 0.5·RLD) km;
    - HYB_XYMAX = 1.3·XY_MAX;
    - HYB_LMAX = min(1.5·L_MAX, 300) km;
    - BLOCK = 12 if M < 7.3, else 24.
  - Steps:
    - the network is trained on the asc look;
    - top-TOP of N_NET proposals by exact two-look slip solve;
    - Adam polish of ITERS steps on asc + desc jointly, each look weighted 1/√P;
    - pick the lowest misfit.
  - For Panama (M 7.7) this gives L ≈ 250 km and XY ≈ 63 km. That is close to the D-2look-lite settings (200 km, 60 km), so the existing Panama-wide net is reused for the synthetic study.
- **Speed study, on synthetics only, Panama-wide geometry.**
  - 24 synthetic events: prior draws filtered to Mw 7.0–7.9.
  - Each look: Okada LOS + a random real pre-event noise map of that look (random sign and flips) + a power-law screen; then deplaned.
  - Settings (N_NET, TOP, ITERS): full (2000, 16, 60); (2000, 16, 30); (1000, 8, 30); (1000, 8, 15); (500, 8, 15); (2000, 16, 0) = no polish.
  - Metrics against the synthetic truth:
    - success = |Δstrike| ≤ 15°, |Δdip| ≤ 15° and |ΔMw| ≤ 0.2 (strike error taken mod 360);
    - also the median wall time.
- **Selection rule.** The cheapest setting (by median time) with success ≥ full − 0.10 and median |ΔMw| ≤ full + 0.05 becomes Q-2look v1.
- **Prediction S1.** The selected setting runs in < 10 s median.
- **Prediction S2.** No-polish loses ≥ 0.2 success against full (the polish matters).
- **Gut.** S1 pass (≈ 5–7 s), S2 pass.
- **Not done here.** No new real-event run counts. The next blind event uses Q-2look v1 as frozen above.

### D-jet. Are the Fisher bands the graded ranks of the jet map? (pre-registered ~03:38, before running; math.md M9 addendum 2)
- **Setup.**
  - Nagqu-like Okada source (192/50/−119, 30 × 15 km) at centroid depth d = 60 km, directly below the window centre.
  - Asc LOS: incidence 35°, heading −12°. Noise-free, float64 on the CPU.
  - Square window [−ℓ, ℓ]², 41 × 41 points (fixed point count).
  - Polynomial of degree k removed (k = 0, 1, 2). Parameters as in D-ident (9).
- **A (Fisher, by window size).**
  - Singular values s_i of the Jacobian of the deplaned map, for ℓ/d = 0.05, 0.1, 0.2.
  - Order of direction i = slope of log s_i vs log ℓ, rounded to an integer. With a fixed point count, s ∝ ℓ^m for jet order m.
- **B (jet, pure algebra).**
  - Taylor coefficients at the window centre of every ∂u/∂θ_j, from a degree-(k+10) polynomial fit on a small window (ℓ/d = 0.05, 61 × 61).
  - Stack orders k+1, k+2, …; record the incremental ranks r_m (relative SVD tolerance 1e-7 after scaling order m by d^m).
- **Prediction J1.** For each k, the multiset of A orders equals the multiset {m repeated r_m times} from B, for all 9 directions.
  - This is the "band multiplicities = graded ranks of the jet map" statement.
- **Prediction J2.** Rank 9 is reached by order k + 5, i.e. any window identifies the full fault given enough noise-free orders.
- **Gut.** J1 pass (it is the singular-value theory of graded matrices if the blocks are generic); J2 pass.
  - Risk: symmetric geometry (source right below the centre) makes some jet blocks non-generic. That would still be informative.
- **Result D-jet (~03:42): J1 FAIL (2 of 3 k pass), J2 PASS.**
  - k = 0: A orders [1,1,2,2,2,3,3,3,3] = B. Incremental jet ranks 2, 3, 4. PASS.
  - k = 1: A orders [2,2,2,3,3,3,3,4,4] = B. Ranks 3, 4, 2. PASS.
  - k = 2: A gives [3,3,3,3,4,4,4,4,**5**] against B [3,3,3,3,4,4,4,4,4]. The last slope is 4.67, which rounds to 5. FAIL by the letter.
  - J2: rank 9 is reached at order 3, 4 and 4 for k = 0, 1, 2. PASS.
- **Post-hoc (labelled), smaller windows ℓ/d = 0.0125–0.05.**
  - k = 0 and k = 1 slopes converge to integers: 1.00, 1.01, 1.98, … and 1.99, 2.00, …, 4.02.
  - k = 2 now rounds to the B multiset, but four slopes (3.32, 3.67, 4.25, 4.42) are still mid-way. The order-3 and order-4 bands mix until ℓ/d is smaller still.
  - So the k = 2 miss looks like "not yet asymptotic", not a wrong count. It stays a FAIL as registered.
- **The clean pattern (post-hoc, the main finding).** The incremental ranks are **maximal at every order**: r_m = min(m + 1, 9 − Σ_{j<m} r_j). Each jet order m has m + 1 coefficients (2-D monomials), and the 9 fault parameters fill them like water fills a staircase. Nothing about Okada beyond "generic" is needed.
  - Example for k = 1 (ramp removed): 3 parameter combinations are seen at order 2 (curvature of the map), 4 more at order 3, and the last 2 only at order 4.
  - Each order costs a factor (ℓ/d)² in Fisher information.
- **Why it matters for the app.**
  - Given window size, source depth and noise, the count tells how many fault combinations an image can pin down.
  - Example: if only order 2 clears the noise, the image constrains 3 numbers, not 9, so the card must show the rest as "not determined".
- Details and the cross-field link: math.md M9 addendum 3.
- **Result D-2look-speed (~03:48): S1 PASS, S2 FAIL.** 24 synthetic Mw 7.0–7.9 events, Panama-wide asc + desc geometry, real pre-event noise.

  | N_NET / TOP / ITERS | success | median |Δstrike| | median |Δdip| | median |ΔMw| | median time |
  |---|---|---|---|---|---|
  | 2000 / 16 / 60 (full) | 0.33 | 9.2° | 9.4° | 0.099 | 12.96 s |
  | 2000 / 16 / 30 | 0.33 | 17.8° | 7.2° | 0.081 | 6.85 s |
  | 1000 / 8 / 30 | 0.38 | 18.2° | 11.1° | 0.103 | 6.65 s |
  | 1000 / 8 / 15 | 0.42 | 12.9° | 9.8° | 0.102 | 3.62 s |
  | 500 / 8 / 15 | 0.46 | 2.9° | 6.7° | 0.112 | 3.33 s |
  | 2000 / 16 / 0 (no polish) | 0.38 | 2.9° | 8.9° | 0.089 | **0.54 s** |

  - By the selection rule, **Q-2look v1 = (2000, 16, 0): network proposals + exact two-look slip solve, no gradient polish, 0.54 s.**
  - S2 fails: the polish does not help on synthetics.
  - With 24 events the success rates carry about ±0.1 sampling error, so all six settings are statistically the same. The rule then simply picks the cheapest, as intended.
- **Reading (post-hoc).**
  - Errors are bimodal: either within a few degrees or wildly off. Strike errors are 0–12° for 13 of 24 events, then 27–180°.
  - Four of the misses are 177–180° strike flips. Five truths have dip 83–88°, where strike+180 is almost the same plane, so some "misses" are near-equivalent geometry. Not re-scored, because the rule did not allow equivalence.
  - Overall success of only ~0.4 on synthetic M7+ is the honest headline. The quake card must show a confidence flag, e.g. agreement between the top proposals, not a single answer.
- **Diagnostic on real Panama with Q-2look v1 (unblinded, not a test).**
  - H1: 286/56/57, Mw 7.51, **0.62 s**. Against USGS 288/66/−17, Mw 7.7: strike 2°, dip 11°, ΔMw 0.19 (just inside 0.2).
  - That is practically the same as the 13.7 s polished run (290/55/36, Mw 7.57).
  - **The rake is wrong in both:** +57 (oblique reverse) against −17 (strike-slip). Strike and dip are recovered, but the mechanism label would be wrong.
  - Next diagnostic to pre-register: is the rake error a two-look limitation for a long, partly-covered rupture (desc valid only 22%), or a bias in the slip solve?
  - The JSON before this run is kept as `sbi_panama-wide.before_nopolish.json`.

### D-conf. Can the quake card know when it is wrong? Plus a rake check (pre-registered ~03:49, before running)
- **Setup.**
  - Q-2look v1 as frozen: N_NET 2000, TOP 16, no polish.
  - Same Panama-wide asc + desc geometry and noise model as D-2look-speed.
  - **48 new synthetic events** (seed 11, Mw 7.0–7.9).
- **Confidence score c.** The share of the top-16 proposals (ranked by two-look misfit after the exact slip solve) whose strike is within 15° (mod 360) and dip within 15° of the best proposal.
- **Success.** As in D-2look-speed: |Δstrike| ≤ 15°, |Δdip| ≤ 15°, |ΔMw| ≤ 0.2.
- **F1.** AUC of c for predicting success ≥ 0.75.
- **F2.** Among events with c ≥ 0.5, success ≥ 0.70, with at least 10 such events.
- **F3 (secondary).** AUC of the relative residual (best two-look misfit ÷ data std; lower = more confident) ≥ 0.70.
- **R1 (rake).** Among successful events, median |Δrake| ≤ 20°.
  - If it passes, the Panama rake error is specific to that event: a multi-segment rupture, or partial coverage (desc 22% valid). It is not a general blind spot of the recipe.
- **Gut.** F1 pass (~0.8), F2 pass, F3 coin-flip, R1 pass.
- **Result D-conf (~03:50): F1 PASS, F2 FAIL (just), F3 PASS, R1 PASS.**
  - 48 new synthetic events. Q-2look v1 success 0.42 (20/48), the same as in D-2look-speed.
  - **F1.** AUC of the agreement score c = **0.93**.
  - **F2.** Success when c ≥ 0.5 is 0.67 on 27 events (needed 0.70): FAIL.
  - **F3.** AUC of the relative residual = 0.91.
  - **R1.** When strike, dip and Mw are right, the rake is right too: median |Δrake| = 3.5°.
- **Post-hoc thresholds (labelled).**

  | condition | events | success |
  |---|---|---|
  | c ≥ 0.75 | 20 | 0.80 (below: 0.14) |
  | c ≥ 0.9 | 14 | 1.00 (below: 0.18) |
  | c ≥ 0.75 and residual ≤ 0.8 | 15 | 0.87 |

  - So the card *can* know when it is wrong. Agreement among the top proposals is a strong flag, and ≥ 0.75 is the natural threshold to pre-register for the next blind event.
- **Panama under this lens (unblinded diagnostic).** Twelve random seeds of Q-2look v1 on the real Panama pair:
  - 9 of 12 give strike 280–299° (USGS 288°);
  - 3 give 112–116° with a shallow dip (a flipped mode);
  - the dip wanders 29–70°, Mw 7.32–7.62, and the rake is +25 to +54 (USGS −17) in the main mode;
  - **confidence c is 0.06–0.50, median ≈ 0.15.**
- **Reading.** Panama is a *low-confidence* event for this recipe. The earlier "4 of 4 pass" (D-2look-lite, and the 0.62 s no-polish run) was partly a lucky draw for dip. Strike is robust; dip, Mw and rake are not.
  - The confidence flag would have said "not sure", which is the honest output.
  - R1 passing on synthetics means the rake error is specific to Panama (a long multi-segment rupture, desc 22% valid), not a general blind spot.
- **Correction to an earlier claim.** The night report called D-2look-lite "the strongest result of the night for quake mode". It should read: strike is robustly recovered; the rest was luck, and the confidence flag catches this.

### H13. Robust per-pixel rate: Theil–Sen within chain segments (pre-registered ~03:51, before running)
- **Idea.**
  - The chain-sum rate uses only the end dates of each segment, so per-date (node) errors at those ends go straight into the rate (M2b).
  - A node error is an *isolated outlier* in the cumulative time series: it shifts date j alone, not the dates after it.
  - A robust slope that uses all dates (Theil–Sen) should ignore isolated outliers, including TXET-style ±0.5–1.5 m date jumps, and also use the interior dates.
- **Method.**
  - Iono "plane", same pairs and pixel mask as the chain sum (valid fraction ≥ 0.6).
  - Per pixel, cumulative displacement is restarted at each gap. Slopes are taken over all date pairs (i, j) in the same segment with t_j − t_i ≥ 60 days. Rate = median slope (at least 3 slopes, else NaN).
- **Evaluation.** Identical to H12: `plane_removed`, 250 m station disc, one median offset removed; Mexico = MMX1 error.
- **Predictions.**
  - **H13a.** Houston desc all-pixel P1 RMS ≤ 10 cm/yr (chain sum: 27.84).
  - **H13b.** CV desc P1 RMS ≤ 3.4 cm/yr (chain sum 3.72).
  - **H13c (no-harm).** Mexico desc MMX1 error ≤ 1.5 cm/yr (chain sum 0.59).
- **Also reported (not judged).** Houston asc and CV asc.
- **Gut.** H13a pass, H13c pass, H13b coin-flip.
  - Risk: CV's big seasonal swing biases the median slope of a one-year series.
- **Result H13 (~03:53): H13a FAIL, H13b FAIL, H13c PASS.** The robust rate is worse almost everywhere.

  | site/pass | chain sum RMS | Theil–Sen RMS |
  |---|---|---|
  | Houston desc | 29.71 | 46.10 |
  | CV desc | 3.72 | 4.54 |
  | Mexico desc MMX1 error | 0.59 | 0.58 |
  | Houston asc | 5.65 | 7.03 |
  | CV asc | 2.78 | 3.18 |

  - (The Houston desc chain sum here is 29.7. The 27.8 in H11 was the rplane variant.)
- **Why (post-hoc look at three stations, iono plane).**
  - TXET is now clean: ±1–2 cm per pair.
  - UHRI and CFJV show **the same** per-pair series to within a few cm: −14, +4, …, +59, −88, +46, +57, −102, +34, +7, +79, −149, +112, …
  - So the errors are not isolated outliers. Almost *every* date carries a ±0.5–1.5 m error, and the error is **common to a large region** (the stations are tens of km apart) relative to the reference area.
  - A robust slope needs a clean majority of dates, and there is none. Hence the failure.
  - It is a spatially smooth error field *per date*. That is exactly what a per-pair tie to GNSS can remove, unlike the rate-level tie in P-tie, which failed.
- **Lesson.** "Node errors telescope" (M1) and "gaps break it" (M2b) are both true. But in Houston desc the node errors are so large and so common that even partial telescoping leaves metres.

### H14. Per-pair GNSS tie: remove a plane fitted to (NISAR pair − GNSS displacement) at the stations, pair by pair (pre-registered ~03:54, before running)
- **Method.**
  - Iono "plane", the same pairs and stations as H12/H13.
  - For each pair, every used station gives r = NISAR pair displacement (median in the 250 m disc, vertical) − GNSS vertical change between the two dates. The GNSS value is the mean of daily solutions within ±2 days of each date; a station without data is dropped from that pair.
  - Fit a plane a + b·E + c·N by least squares (at least 5 stations, else a constant; at least 2, else no correction).
  - Subtract the plane at the station locations, then form the chain-sum rate at each station over the valid pairs.
- **Honest scoring.** Leave-one-out: the held-out station is excluded from every pair's fit. P1 = RMS of the LOO errors after removing one median offset, the same as P-tie.
- **Predictions.**
  - **H14a.** Houston desc LOO RMS ≤ 3 cm/yr (P-tie 103; chain sum 29.7).
  - **H14b.** CV desc LOO RMS ≤ 3.0 (chain sum 3.72; P-tie 3.75).
  - **H14c.** Houston asc LOO RMS ≤ 3.0 (chain sum 5.65).
- **Gut.** H14a pass (the error is common-mode and smooth), H14c pass, H14b coin-flip (CV has only 12–13 stations).
- **Result H14 (~03:58): H14a FAIL, H14b FAIL, H14c FAIL.**

  | site/pass | stations | no tie, RMS (cm/yr) | per-pair tie LOO, RMS (cm/yr) |
  |---|---|---|---|
  | Houston desc | 33 | 17.75 | **6.33** |
  | CV desc | 13 | 7.73 | 5.80 |
  | Houston asc | 27 | 7.85 | 8.74 |
  | CV asc | 13 | 2.81 | 2.32 (not judged) |

  - The "no tie" baselines here are station-level chain sums without the `plane_removed` step, so they differ from the P-port numbers (e.g. CV desc 3.72 there). The station sets differ slightly too (≥ 3 valid pixels per pair).
  - The pre-registered limits were the absolute values, and all three judged cases miss them.
- **Reading (post-hoc).**
  - In Houston desc the tie removes the monster errors: ZHU1 +92.9 → +7.9, TXET +25.8 → +14.8, TXNH and TXTG +14.3 → −1.5 / −1.3. The worst station drops from 94 to 15 cm/yr.
  - But it adds 5–12 cm/yr at stations that were fine (WEPD −10.7, LKHU −12.1, GSEC +10.6, AULT +10.9).
  - So the per-date error field is smooth on large scales but **not planar**. A plane through 30 stations removes its average tilt and leaves ±10 cm/yr bumps between them.
  - In Houston asc, where errors were small to begin with, the tie only adds noise.
- **What would work better (not tested).**
  - A smoother but richer per-date field: a low-order 2-D spline or kriging of the GNSS residuals per date, with LOO choosing its stiffness.
  - Or a date-level field estimated from the radar itself (ML-2). Both are on the list for the hackathon, not tonight.

### D-depth. Which "deepest" number may the Mexico map say out loud? (diagnostic 11.11 afternoon, not a pre-registered test)
- **Problem.** The manifest said "deepest 28 cm" (1st percentile of the smoothed, deramped field); single strokes carried up to −46 cm/yr; one grid cell −72 cm/yr.
- **Diagnostic** (`scripts/diag_depth.py`): raw stack rate p1 −30.7 cm/yr; deepest 500 m block (≥ 50% valid) −59 ± 12; the "most sure" block by rate + 2σ was −53.5 ± 3.8 at −99.031, 19.439.
- **But** that block sits in a low-coherence pocket (valid share 0.1–0.5 in a 2 × 2 km patch) surrounded by −30…−35 on every side. Smoothed field there: −20…−25. At MMX1 raw, smoothed and GNSS all agree (−26.9 / −26.4 / −27.0).
- **Rule now (labelled post-hoc):** headline = ~1.5 km neighbourhood (median of 3 × 3 blocks, each ≥ 80% valid, ≥ 6 good blocks) with the most negative rate + 2σ. Mexico asc: **−31.2 ± 2.5 cm/yr** (at least −26.1) at −99.044, 19.457; p1 −30.7.
- **Reading.** Honest public number: "the east of the city sinks about 30 cm a year". Faster local spots (−38…−43 in the SE, Chalco side) exist in the data, but with weaker coverage. Not claimed until a second pass or GNSS confirms.

### L / C. Old lake and cracks — pre-registration (written 11.11 evening, BEFORE downloading the layers)
- **Data (official, open):** Atlas de Riesgos CDMX, `www.atlas.cdmx.gob.mx/datosAbiertos/`:
  `SPC_Zonificacion_Geotecnica.geojson` (zones Lomas / Transición / Lago), `Fracturas.geojson` (mapped ground fractures),
  `SPC_Hundimientos_ITRF.geojson` (their own subsidence layer), `SPC_CG_Fracturamiento.geojson` (Iztapalapa 2013).
  Not yet looked at. Zone names are mapped to Lomas / Transición / Lago from the attribute table **before** any NISAR number is read; the mapping is written down here.
- **Radar fields (Mexico asc, as on the map):** `R_raw` = stack rate on the DEM grid (the one the headline uses);
  `R_s` = smoothed (400 m), deramped rate behind the strokes. Gradient `G = |∇R_s|` in cm/yr per km.
- **L1.** Of all valid pixels with `R_raw ≤ −10 cm/yr`, at least **80%** lie in the Lago zone.
- **L2.** Median `R_raw`: Lago < Transición < Lomas, and Lago − Lomas ≤ **−5 cm/yr**.
  (L1–L2 are low-risk: the literature says the lakebed clay sinks. They check that our map tells the known story.)
- **C1 (the real test).** Cracks sit where the ground sinks *unevenly*, not just fast.
  Points every 100 m along `Fracturas` lines, inside the bbox, on valid pixels (fewer than 50 points → "not judged").
  Statistic: median `G` at the points. Null: 500 rigid random shifts of the whole crack set, dx, dy uniform in [−3, 3] km, |shift| ≥ 1 km,
  keeping only shifted points on valid pixels in Lago or Transición. **Pass:** observed median > 95th percentile of the null.
- **C2 (exploratory, no claim):** the same test with `R_raw` instead of `G` (are cracks also where it is fastest?). Reported, not judged.
- **H (exploratory with a threshold):** if `SPC_Hundimientos_ITRF` carries numeric rates, Spearman ρ between them and `R_s` at the same places ≥ **0.6**
  (ranks only: their epoch and reference differ from ours).
- Thresholds are not moved after the run. Anything else found is labelled post-hoc.
- **Mapping fixed after reading the attribute tables, before any NISAR number** (11.11 evening):
  `INTENSIDAD` "Zona I" (Lomas) → Lomas, "Zona II" (Transición) → Transición, "Zona III" (Lacustre) → Lago. 103 polygons, EPSG:6365.
  `Fracturas`: 6964 lines (EPSG:3857); all used (Fractura 6036, Grietas 922, Fractura Inferida 6).
  `SPC_Hundimientos_ITRF`: 5 polygons with classes "02-10", "11-20", "21-30" cm/yr → H is run per valid pixel with class midpoints 6 / 15.5 / 25.5
  against −`R_s` (Spearman, threshold 0.6 unchanged).
- **Results (11.11 evening; `scripts/test_lake_cracks.py`, output `data/processed/mexico-city-asc/test_lake_cracks.json`):**
  | test | number | threshold | verdict |
  |---|---|---|---|
  | L1 fast (≤ −10 cm/yr) pixels in Lago | **75.0%** (Lomas 13.0, Transición 10.3, no zone 1.7) | ≥ 80% | **fail** |
  | L2 medians Lomas / Transición / Lago | −1.9 / −3.9 / **−15.2** cm/yr; Lago − Lomas −13.3 | order + ≤ −5 | **pass** |
  | C1 median gradient at cracks (10 880 pts) | 2.82 cm/yr/km; null median 2.50, p95 3.31; 27% of shifts ≥ obs | > p95 | **fail** |
  | C2 (exploratory) median rate at cracks | −7.4 vs null −8.4 cm/yr | — | cracks are *not* in the faster places |
  | H their classes 2–10 / 11–20 / 21–30 vs our R_s | ρ = **0.70**; our medians −4.0 / −15.1 / −20.9 | ≥ 0.6 | **pass** |
- **Post-hoc (`scripts/posthoc_lake_cracks.py`, not pre-registered):**
  - Faster than −15 cm/yr: 86.6% in Lago. The 13% of ≤ −10 pixels on "Lomas" are unexplained yet (hill-slope artefacts? zoning edges outside CDMX?). Worth a look before the map says "only the lakebed sinks".
  - C1 made symmetric (observed points only in Lago/Transición, like the null): 24% of shifts ≥ obs. Lago only: 39%. Finer gradient (raw rate, 150 m): 22% / 42%. **No variant shows cracks tied to steep differential subsidence.**
  - Likely reasons (not tested): the crack map follows where field teams looked (81 survey sources, dense in the western alcaldías on firm ground), and cracks are 10–100 m features while our field is 50 m pixels smoothed to 150–400 m.
- **Reading.** The map tells the known story (old lakebed sinks ~15 cm/yr, hills ~2; independent CDMX classes agree, ρ 0.70). The crack link is **not shown** by these data; the map must not imply it. Cracks may still be drawn as context, labelled "mapped fractures (CDMX atlas)", without a causal caption.
  - **Where the "sinking hills" are (post-hoc):** almost all of them sit in the SW corner, −99.27…−99.31, 19.21–19.30
    (Ajusco / Sierra de las Cruces, ~2 800–3 500 m), a uniform −13 cm/yr over 2 km blocks. Mountains of rock do not sink that fast:
    this looks like a height-correlated atmosphere delay (or the scene-edge ramp), not motion. **To do:** fit rate vs DEM height
    on Lomas only and remove it (a standard tropospheric fix), then check that the strokes vanish there and L1 rises; label it post-hoc.

### T. Height-correlated delay fix for the "sinking mountains" — criteria written BEFORE the run (11.11 evening)
- **Origin is post-hoc** (found while reading L1), so the fix is judged only by criteria fixed here, before it is run.
- **Method** (config `tropo: {"mode": "height", "stable_above_m": 2400}`): per date, on coherent pixels higher than 2400 m
  (the basin floor is ~2240 m), fit `d = a + b·x + c·y + e·(h − h0)` (least squares, two 3σ clips); subtract only `e·(h − h0)`
  from that date's cumulative displacement; same for the stack rate. `h0` = median DEM height inside `ref_bbox`, so the reference stays at zero.
- **T1.** SW mountain blocks (lon −99.31…−99.26, lat 19.20…19.30): median `R_raw` within ±4 cm/yr after the fix (before: about −13).
- **T2.** MMX1 radar rate within 3 cm/yr of GNSS −27.0 after the fix (before: −26.9).
- **T3.** Headline moves by less than 3 cm/yr (before: −31.2 ± 2.5).
- **T4 (reported, post-hoc by nature).** L1 re-run.
- If T2 or T3 fails, the fix is **not** used for the map.
- **Results T (11.11 evening, `scripts/test_tropo.py`):** fit on 102 938 px above 2400 m, h0 = 2401 m; height term **−0.51 cm/yr per 100 m**.
  | criterion | before | after | verdict |
  |---|---|---|---|
  | T1 SW mountains median | −9.9 | **−5.1** cm/yr | **fail** (limit ±4): only half of the mountain signal is linear in height |
  | T2 MMX1 vs GNSS −27.0 | −26.9 | −27.8 | pass |
  | T3 headline shift | −31.2 | −32.1 (Δ −0.9) | pass |
  | T4 L1 re-run (post-hoc) | 75.0% | **80.7%** | above the 80% line; L2 −1.6 / −4.7 / −16.1; H ρ 0.70; C1 unchanged (28% of shifts ≥ obs) |
- **Decision (by the rule above):** T2 and T3 pass → the fix is used for the map. T1 failing means the mountains still carry about −5 cm/yr
  that is not explained by height: the map should keep strokes there faded or say "uncertain in the mountains". Not claimed as motion.

### GC. GCOV first look: do the fields "breathe" in HV? (13.11, criteria written BEFORE any GCOV value was read)
- **Why.** The challenge asks which polarizations show which change. GCOV (RTC gamma0, HH and HV) is the product for crops, forest
  and water; we have used only GUNW so far. Central Valley already has a coherence plate, so GCOV would add the second lens.
- **Data.** NISAR L2 GCOV over the `central-valley-asc` bbox, dominant ascending track, about one date per month (Jan–Oct 2026,
  at most 10 dates), frequency B grid (coarse, small reads), terms HHHH and HVHV, values in dB = 10·log10(gamma0).
  Pixels where `mask` flags invalid or no data are dropped.
- **Classes.** WorldCover fractions already on the plate grid (`cv_worldcover_frac.npz`): cropland ≥ 0.6, built-up ≥ 0.3 (as in X).
- **GC1.** Cropland: max − min over dates of the per-date median HV ≥ 3 dB.
- **GC2.** Built-up: same range ≤ 1.5 dB (a stable target; if it fails, dates are not comparable and no plate uses absolute GCOV levels).
- **GC3.** Cropland: HV range > HH range (HV is the volume / canopy channel).
- Reported, not judged: per-date medians, timing and MB read per date.
- **Results GC (13.11, `experiments/probe_gcov.py`):** track asc 34/20, 10 dates (one per month, Jan–Oct 2026),
  frequency B grid 80 m, window 1007 × 1376 px, ~12.8 MB cached per date (~25 s each). The `mask` layer was 1 on every
  window pixel (no invalid pixels inside this bbox).
  | criterion | value | limit | verdict |
  |---|---|---|---|
  | GC1 cropland median HV, max − min | 2.34 dB (−20.7 Oct … −18.4 Jul) | ≥ 3 | **fail** |
  | GC2 built-up median HV, max − min | 0.97 dB | ≤ 1.5 | **pass**: dates are comparable |
  | GC3 cropland HV range > HH range | 2.34 vs 2.45 dB | HV > HH | **fail** |
- **Post-hoc (not pre-registered), per pixel instead of the class median:** median per-pixel range over the 10 dates.
  | class | HH | HV |
  |---|---|---|
  | cropland | 5.3 dB | 4.2 dB |
  | built-up | 2.8 dB | 2.9 dB |
  - The built-up figure is the noise floor of a single 80 m pixel.
  - 86% of cropland pixels move ≥ 3 dB in HH, 74% in HV; built-up 43%.
  - HV beats HH in only 30% of cropland pixels.
- **Reading.**
  - Single fields do breathe, well above the noise floor: about 5 dB against 3 dB. The class median hides it, because fields
    on different calendars average out (the pictures show fields switching on and off between January and June).
  - The "HV is the crop channel" expectation failed here: HH moves more, likely because of irrigation and soil moisture and
    row structure (an untested explanation).
  - For a plate: show per-field HH and HV curves or a change map, not one valley-wide curve. Say plainly that HH, not HV,
    carries most of the seasonal signal in this valley.

## R. Reykjanes (Svartsengi uplift, 10 Nov 2023 dike): log of tests, recorded 2026-11-14

**Status of this block: post-hoc record.** The Reykjanes tests below were run in a separate working folder before this
log existed, so their predictions and thresholds were **not** written here in advance. What was fixed before each test is
stated where the scripts (now in `scripts/reykjanes/`) show it; everything else is labelled **post-hoc**.
Numbers come from `output/reykjanes/results.html` and `output/reykjanes/ml/*.json`. The plate is `web/reykjanes-scroll.html`.

### R1. Independent GNSS check of the NISAR vertical rate (2025-11 → 2026-10)
- Fixed before: stations SENG, AUSV, GEVK, VOGC are not used in the rate, the reference (NYLA) or any correction
  (`build_breathing_data.py` docstring); the NISAR value is the median of 3 × 3 cells at the station.
- Post-hoc: the agreement threshold |z| < 2 is written down only now.
- Result: **1 of 4 agrees.** SENG 20.1 ± 1.7 vs 21.0 ± 1.2 cm/yr (z −0.4); AUSV 3.5 vs 7.6 (z −2.7); GEVK −0.3 vs 7.0
  (z −4.0); VOGC 4.3 vs −1.8 (z +2.8). AUSV and GEVK are 4.4 and 5.5 km from SENG, by Grindavík. The per-cell σ
  (1–1.7 cm/yr) is too small for these misfits: it does not hold systematic errors. Cause not found.

### R2. Source of the uplift: MCMC vs NPE, homogeneous vs layered crust
- Fixed before: MCMC with the NISAR screen covariance is the reference; NPE is judged by 68 % coverage on unseen
  simulations (ideal 68 %) and by its shift from MCMC in σ.
- Result: Mogi depth 4.1 [4.0–4.3] km, 18 [17–19] Mm³/yr; layered crust (FE) 4.7 [4.5–4.8] km, 15 [14–16] Mm³/yr. The NPE
  flow matches MCMC for Mogi (shifts ≤ 0.5σ); for the elongated reservoir its depth is 1.7σ shallower (open issue).
  NPE with GP screens is badly calibrated (27–52 %) and is not used.

### R3. Dike of 10 Nov 2023 from two Sentinel-1 pairs + 4 GNSS, against Sigmundsson et al. 2024 (H1–H8)
- Reference: Sigmundsson et al. 2024, Science 383:1228: dike 130–139 Mm³, Svartsengi −76…−82 Mm³, Mogi depth 4.0–4.2 km.
- Fixed before running NPE v3: the model (variable opening, graben mask ±1 km, FE source), the 25 % no-dike class.
  **Post-hoc:** the hierarchical noise range (0.1–40 cm²) was chosen after seeing ~30 cm² in Nov 2023; the 5 % GNSS error
  and the 256 PCA components after the first inversions; v2 → v3 → v4 were each made after the previous one failed on
  this event.
- Verdicts (no thresholds were set in advance, so these are descriptive): H1 GNSS removes the dike ↔ chamber trade-off:
  partly; H2 the NPE–MCMC gap comes from the model, not NPE: **no** (χ² of NPE samples 2787 vs 865); H3 the noise covariance
  was too high: yes; H4 NPE does not invent dikes: yes (quiet pair P(no dike) 0.86); H5 calibrated (SBC): mostly (15 of 16);
  H6 the signal is not from the prior: partly; H7 transfers to the 2024 events: partly (29 May 2024: 7 [3–14] vs 13–19 Mm³);
  H8 layered crust vs the paper: deeper by ~1 km, as expected.
- Result NPE v3: dike 150 [131–170] Mm³ (median above the paper's range; the interval covers 131–139 of it), Svartsengi
  −77 [−92…−63] Mm³ (covers the paper's values), depth ~6 km (paper 4.0–4.2, Mogi).
- v4 blind: trace 1.9–2.9 km from FIRMS hotspots (random points 5.1 km). "Almost blind": see the leakage table in the
  math window (`docs/uk/reykjanes-math.html`, section 9).

### R4. Plate choices made after looking (post-hoc)
- The Hafnir row was left out of the score: −13.5 cm/yr at the west edge, judged an edge/plane artefact after seeing it.
- FIRMS dots: one dot per 400 m cell with ≥ 8 VIIRS detections in the eruption windows; the threshold (first 1, then 4,
  then 8) was raised after seeing single detections cover the map.
- The σ colour ramp of the trust cells was softened after seeing it (display only, no number changed).
- Coherence-loss enrichment: the plate quotes the measured ×4.9 for the 2023 dike only.

### R5. Next tests, pre-registered 2026-11-14 (before running)
- **R5a. v3 with Mogi instead of the FE source.** Prediction: Svartsengi depth median 5.0–6.0 km (H8 gave 5.5 with the
  distributed model). Pass: median in 5.0–6.0. Fail: outside. The dike volume is reported, not judged.
- **R5b. v4 kinematics + v3 corridor and mask.** Prediction: the v3/v4 gap in the dike centre (north) shrinks below 0.5 km
  and the volume median falls to 130–160 Mm³. Pass: both. Partial: one. Fail: neither. This says how much of the gap is
  the prior.
- **R5c. Truly blind test on future NISAR pairs.** Model and hyper-parameters frozen at the commit that adds this line.
  Pass on a quiet pair: P(no dike) ≥ 0.8. Pass on a dike pair: trace within 2 km of the vents later mapped by IMO.

## K. Kakhovka reservoir bed: L-band radio physics with GCOV (`edoTenseii` + Claude)

Goal: use the drained Kakhovka reservoir bed (dry since June 2023) as a natural test range for L-band radio physics:
mirror (water), volume (young willow and poplar forest that grew on the bed), rough bare ground (sand and silt), all side
by side, seen by four NISAR tracks at different incidence angles. Afterwards the same code runs on a quad-pol site
(northern Wisconsin, track D113/065) for validation.

**Pre-registered 2026-11-14 ~19:30 (Dnipro time), before any GCOV backscatter value over Kakhovka was read.**
What had been seen before writing this, and nothing else:
- ASF search counts: ~50 GCOV scenes (`2005_DHDH`, HH + HV, band A 20 MHz on a 20 m grid, band B 5 MHz on 80 m) from
  4 tracks (A42/26, A143/25, D36/64, D108/65), Nov 2025 – Nov 2026; Sentinel-2 L2A cloud-free scenes exist (tile 36TXT).
- The metadata of one scene (A42, 2026-10-27): `noiseCorrectionApplied = False`, `faradayRotationApplied = False`,
  `rfiMitigationApplied = False`, `isFullCovariance = False`; product noise-equivalent backscatter (`noiseEquivalentBackscatter`)
  medians −15 dB HH / −14 dB HV in band A and −21 / −19 dB in band B; `rfiLikelihood` 2.4e6 (A) and 9.2 (B) although the
  attribute says the value lies in [0, 1]; crosstalk all zeros; elevation antenna pattern exactly 1.

**Window.** W 34.75, S 47.35, E 35.20, N 47.70 (the upper bed, Velykyi Luh). Chosen to leave out settlements and
industrial sites; only natural surface classes are analysed. (This paragraph was reworded on 2026-11-14 before the repo
goes public; the window and every rule below are unchanged.)

**Classes** (independent of the radar):
- *bed*: ESA WorldCover 2021 class "permanent water" (the reservoir before it was drained).
- *water*: Sentinel-2 L2A scene classification SCL = 6 on both 2026-08-23 and 2026-02-27, NDWI ≥ 0.3 in August, and
  every S2 pixel under a GCOV pixel water (no mixed edge pixels).
- *vegetated bed*: bed, NDVI ≥ 0.6 on 2026-08-23, not water on either S2 date.
- *bare bed*: bed, NDVI < 0.25 on 2026-08-23, not water on either date.
- *cropland control*: WorldCover 2021 cropland outside the bed.
- Frozen / unfrozen dates: Open-Meteo hourly 2 m air temperature at the overpass hour, plus the 24 h mean before it.
  Frozen: overpass T < −2 °C and 24 h mean < 0 °C. Unfrozen: overpass T > +2 °C. Others are left out of K2 and K4.

**K0. Product audit (reported, not judged).** For every scene: `rfiLikelihood` per band and polarization, Faraday rotation,
the processing flags, whether crosstalk is all zero and the antenna pattern all one.

**K1. Is the water at the noise floor the product states?** dB values = 10·log10 of gamma0 (linear power).
NEB interpolated from its coarse grid to each GCOV pixel.
- K1a. Band B, every scene: median over water of 10·log10(γ⁰_HV / NEB_HV) within ±2 dB in ≥ 80 % of scenes.
  Pass means: noise is not subtracted, calm water returns nothing in HV, and the stated NEB is the real floor in γ⁰ units.
- K1b. Band B: median water HH ≥ median water HV − 0.5 dB in ≥ 80 % of scenes (HH carries some Bragg return from
  wind ripples, HV only noise).
- K1c. Band A, one scene per track (4 scenes): the same ratio as K1a. Reported per scene; judged with the same ±2 dB.

**K2. Does the 24 cm wave see leaves?** Band B, per track, unfrozen dates only (rule above), per-date class median in dB,
seasonal amplitude = max − min of these medians.
- K2a. Vegetated bed: HV amplitude ≤ 2.0 dB on ≥ 3 of 4 tracks (leaves ≈ 5 cm ≪ λ are weak Rayleigh scatterers;
  trunks and branches dominate and do not change with season).
- K2b. Vegetated-bed HV amplitude < cropland HV amplitude on ≥ 3 of 4 tracks.

**K3. Angular signature γ⁰(θ).** Band B, dates 2026-06-01 … 2026-09-30, unfrozen. Per track: class median γ⁰ (dB) and
class median incidence angle θ (from the `radarGrid/incidenceAngle` cube at the bed height).
- If the four tracks span < 5° of θ over the classes, K3 is **inconclusive** and not judged.
- K3a. Vegetated bed: max − min of the four track medians ≤ 1.5 dB in HV and in HH (volume scattering, γ⁰ almost flat in θ).
- K3b. Bare bed: least-squares slope of HH γ⁰ against θ ≤ −0.10 dB/deg, and steeper (more negative) than the
  vegetated-bed HH slope. Ascending and descending tracks also look from opposite sides; this confound is stated, not removed.

**K4. Frozen ground lowers backscatter (permittivity of ice ≈ 3.2 against ≈ 80 for liquid water).** Band B, bare bed and
cropland control, per track with ≥ 1 frozen and ≥ 2 unfrozen dates in Nov 2025 – Mar 2026:
- K4. Median HH on frozen dates is lower than the median of the unfrozen dates of the same months by ≥ 1.5 dB on at
  least half of the eligible track × class cases. Fewer than 2 eligible cases: inconclusive.

**Amendment (2026-11-14 ~20:50, before any class-based radar value was computed; only the K0 audit had been printed).**
The class layers came out empty: the Sentinel-2 scene of 2026-02-27 (S2B_36TXT_20260227_0) covers only ~21 % of the
window (79 % no-data at the swath edge), and the code had picked the edge-of-swath scene of 2026-08-23 too. Changes:
the second S2 date is **2026-03-11** (S2C_36TXT_20260311_0: 0 % no-data, 0.1 % cloud, no snow; early spring, leaf-off),
and each date uses the scene with the least no-data. Every class rule and threshold above stays as written.
Implementation details fixed now, also before any class-based radar value: *bed* = WorldCover water fraction ≥ 0.9 of
the 80 m cell; *cropland control* = cropland fraction ≥ 0.9 and no water; a window is used for a class only with ≥ 50
class pixels.

**Timing note.** The pre-registration above and the results below went into git in one commit (efbd7ee, 21:31), so
git cannot show that the predictions were written at ~19:30; that time is our own record. The planned Wisconsin run will
commit its pre-registration on its own, before any data is read.

**Results K (2026-11-14 ~21:30, `experiments/k_fetch.py`, `experiments/k_tests.py`).** 48 band B windows (A42: 17 dates
Nov 2025 – Oct 2026; A143, D36, D108: Jun – Nov 2026 only), 4 band A windows (Aug 2026). Class pixels on the 80 m grid:
water 2 631, vegetated bed 41 187, bare bed 2 758, cropland 73 973.

| test | value | limit | verdict |
|---|---|---|---|
| K0 `rfiLikelihood` outside [0, 1] | A: 41 of 48 scenes (max 3.3e12); B: 39–40 of 48 (max 46) | reported | product field broken in later processing (X05009 Nov 2025 scenes are inside [0, 1]) |
| K0 calibration layers | crosstalk all zero, elevation antenna pattern all one, Faraday 0 and not applied, noise not subtracted, in all 48 | reported | placeholders in this provisional product |
| K1a water HV / NEB_HV within ±2 dB | 10 % of scenes; median ratio +6.4 dB | ≥ 80 % | **fail** |
| K1b water HH ≥ HV − 0.5 dB | 40 % of scenes | ≥ 80 % | **fail** (HV above HH over water in most scenes) |
| K1c band A water HV / NEB_HV | −5.9, +7.2, +6.7, +7.1 dB | ±2 dB | **fail** on 4 of 4 |
| K2a vegetated-bed HV amplitude ≤ 2 dB | 4.5, 13.5, 5.1, 6.2 dB (D36, A42, D108, A143) | ≥ 3 of 4 | **fail** |
| K2b vegetated-bed HV amp < cropland HV amp | 4 of 4 | ≥ 3 of 4 | pass |
| K3 θ spread | 11.0° (34.6–45.4°) | ≥ 5° | judged |
| K3a vegetated-bed range across tracks | HV 3.0 dB, HH 0.8 dB | ≤ 1.5 both | **fail** |
| K3b bare-bed HH slope | +0.12 dB/deg (vegetated +0.07) | ≤ −0.10 and steeper | **fail** |
| K4 frozen − unfrozen HH (A42 only, 2 cases) | bare +1.4 dB, cropland +1.0 dB | ≤ −1.5 on half | **fail** |

**Diagnostic (post-hoc, not pre-registered).**
- Over this window the backscatter levels are not physical for these surfaces: HV is close to or above HH in every
  class, water included (median HV − HH +0.5 dB on cropland and +1.5 dB on water; water HV −19 … 0 dB), and all classes
  move together from date to date by up to ~10 dB.
- Control lake (`experiments/k_diag.py`): open water of Lake Balaton (Hungary), box 17.72–17.78 E, 46.84–46.88 N,
  three band B scenes Jul 2026, same product type: water HH −15.4 / −19.2 / −16.0 dB, HV −22.3 / −23.7 / −23.2 dB,
  HV − NEB_HV +2.2 … +3.2 dB. There the physics holds (HV well below HH, close to the stated floor). Central Valley
  cropland (test GC) also sat at HV ≈ −20 dB.
- The cause was not investigated. Confound noted: the A42 winter dates were processed as X05009/X05010, the others as
  P05023/P05036.
- **Consequence:** K1–K4 here say nothing about the surfaces, fails or pass (the K2b pass included). The surface tests
  move to the Wisconsin site (planned, not run yet) with the same code and thresholds. The K0 audit stands: it describes the public
  product, not the site.
- Code note (2026-11-14, review of PR #4): the K4 date filter now matches the text (2025-11-01 … 2026-03-31); the K4
  result is unchanged.

## W. Northern Wisconsin: L-band radio physics with full quad-pol GCOV (`edoTenseii` + Claude)

Goal: run the surface tests that the Kakhovka bed could not answer (section K: the backscatter levels there were not
physical) on a clean site, and add the tests that need all four polarizations. NISAR GCOV `2005_QPDH` here carries the **full
covariance matrix** in band A (10 terms, `isFullCovariance = True`, complex off-diagonals) and HH/HV in band B.
(The goal sentence and the heading were reworded on 2026-11-14 to match section K; no test or threshold changed.)

**Pre-registered 2026-11-14, written ~22:00 (Dnipro time), before any GCOV backscatter value over this site was read.**
What git shows (review of PR #6): the text was committed locally on its own at 22:02 (b5fe3df; it said "~23:50" by
mistake). That commit was rewritten by two rebases during the data run (22:50 and ~23:11) to reword the goal and the
heading; the branch was first pushed at 23:45, after the results commit (23:31). So, as in section K, the timing is our
own record. The original commit is kept unchanged at `kiberatlet/w-prereg-original`; `git diff b5fe3df <this commit>
-- notes/experiments.md` shows that only the heading, the goal sentence and this timing line changed, no test or
threshold.
Seen before writing, and nothing else: ASF scene lists (52 QPDH scenes, tracks A76/24: 18, A148/24: 7, D41/65: 10,
D113/66: 16; no QPDH scenes Feb – mid-June 2026), the term lists and flags of one file (D113, 2026-11-01:
`rfiLikelihood` = inf and 1e35), ESA WorldCover shares (used to choose the windows), Sentinel-2 scene lists, and the
array shapes (not values) of one test read (A76, 2026-10-29, band A).

**Windows.** Full window W −91.00, S 45.02, E −90.56, N 45.28 (WorldCover: tree 67 %, cropland 20 %, grass 11 %,
water 1.2 %), band B 80 m, all 52 scenes. Sub-window W −90.77, S 45.17, E −90.67, N 45.24, chosen from WorldCover only as the
~8 × 8 km box with most open water that keeps ≥ 30 % tree and ≥ 10 % cropland (water 12 %, tree 67 %, cropland 12 %):
band A 20 m, full covariance, one scene per track and calendar month (the earliest of the month).
Terrain height for the incidence cube 350 m above the ellipsoid.

**Classes** (radar-independent; Sentinel-2 L2A 2026-07-14 = leaf-on, 2025-11-03 = leaf-off and snow-free, tiles 15TXL + 15TXK):
- *water*: SCL = 6 on both dates for every S2 pixel under the radar pixel, and July NDWI ≥ 0.3.
- *deciduous forest*: WorldCover tree fraction ≥ 0.9, July NDVI ≥ 0.75, November NDVI ≤ 0.45.
- *conifer forest*: WorldCover tree fraction ≥ 0.9, July NDVI ≥ 0.75, November NDVI ≥ 0.60.
- *cropland*: WorldCover cropland fraction ≥ 0.9, no water.
- Frozen / unfrozen dates: the section K rule (Open-Meteo hourly 2 m temperature). A window is used for a class only with
  ≥ 50 class pixels. dB = 10·log10 of gamma0.
- Corn against soybean (USDA CDL) is **not tested**: the CDL service returned errors on 14.11.

**W1. Noise floor** (repeat of K1). W1a band B: water HV / NEB_HV within ±2 dB in ≥ 80 % of scenes. W1b band B: water
HH ≥ water HV − 0.5 dB in ≥ 80 %. W1c band A: water HV / NEB_HV and VH / NEB_VH both within ±2 dB in ≥ 80 % of scenes.

**W2. Does the 24 cm wave see leaves?** Band B, unfrozen dates. A track is eligible with ≥ 4 unfrozen dates including
≥ 1 in Jun–Aug and ≥ 1 in 15 Oct – 30 Nov. Fewer than 2 eligible tracks: inconclusive.
- W2a. Deciduous-forest HV amplitude (max − min of per-date medians) ≤ 2.0 dB on ≥ 2/3 of eligible tracks.
- W2b. Deciduous-forest HV amplitude < cropland HV amplitude on ≥ 2/3 of eligible tracks.
- W2c. Deciduous HV amplitude − conifer HV amplitude ≤ 1.0 dB on ≥ 2/3 of eligible tracks (losing all leaves adds little).

**W3. Angular signature.** Band B, Jun–Sep 2026, unfrozen. Per track, class medians of γ⁰ and θ. Spread of θ over
tracks < 5°: inconclusive. W3: for deciduous and for conifer forest, the range over tracks ≤ 1.5 dB in HV and in HH.

**W4. Frozen ground and frozen trees.** Band B, Oct 2025 – Jan 2026, per track and class (deciduous, conifer, cropland),
eligible with ≥ 1 frozen and ≥ 2 unfrozen dates. W4: median HH frozen − median HH unfrozen ≤ −1.5 dB on ≥ half of the
eligible cases. Fewer than 2 eligible cases: inconclusive.

**W5. Reciprocity HV = VH** (band A). For forest (deciduous ∪ conifer) and cropland, per scene:
- W5a. |median 10·log10(HV / VH)| ≤ 0.5 dB in ≥ 80 % of scenes.
- W5b. Cross-pol coherence |⟨HV·VH*⟩| / √(⟨|HV|²⟩⟨|VH|²⟩) of the class means ≥ 0.9 in ≥ 80 % of scenes.

**W6. Co-pol ratio against scattering theory** (band A).
- W6a. Forest: |median VV − median HH| ≤ 1.0 dB in ≥ 80 % of scenes (volume scattering, no preferred polarization).
- W6b. Water: VV − HH ≥ +1.0 dB (Bragg / small-perturbation model) in ≥ 80 % of the scenes where water HH and VV are both
  ≥ 3 dB above their NEB. Fewer than 5 such scenes: inconclusive.

**W7. Cloude–Pottier H/α** (band A, unfrozen scenes, from the class-mean 3 × 3 coherency matrix). W7: forest entropy
H ≥ 0.7 and 40° ≤ ᾱ ≤ 55° in ≥ 80 % of scenes. Cropland and water: reported, not judged.

**W8. Faraday rotation** (band A, forest class mean, Bickel–Bates estimator from the full covariance:
Ω = ¼·arg(⟨|b|²⟩ − ⟨|a|²⟩ − 2j·Re⟨a b*⟩) with a = HV − VH, b = HH + VV). Reported, not judged; expectation |Ω| ≤ 3°.
Crosstalk is not corrected in the product (all-zero layers, section K), so a non-zero Ω may be crosstalk, not ionosphere.

**Results W (2026-11-14 ~23:40, `experiments/w_fetch.py`, `experiments/w_tests.py`).** 52 band B windows (A76: 18,
D113: 16, D41: 11, A148: 7), 27 band A full-covariance windows (one per track and month). Class pixels on the 80 m grid:
water 76, deciduous 15 684, conifer 14 990, cropland 21 538. Fix before the run, no rule changed: the products are on UTM
zone 16 (EPSG 32616), so the band A class layers are built on each window's own zone.

| test | value | limit | verdict |
|---|---|---|---|
| W1a water HV / NEB_HV within ±2 dB (band B) | 0 of 42: NEB layer zero or NaN under the water pixels in 30; finite in 12, all outside (+3.3 … +5.1 dB) | ≥ 80 % | **fail** |
| W1b water HH ≥ HV − 0.5 dB | 42 of 42 | ≥ 80 % | pass |
| W1c band A water vs NEB | no water pixels on the band A grid (see below) | | inconclusive |
| W2a deciduous HV amplitude ≤ 2 dB | 0.8, 2.8, 1.8 dB (D41, A76, D113; A148 not eligible) | ≥ 2/3 | pass |
| W2b deciduous HV amp < cropland HV amp | 3 of 3 (cropland 2.6, 4.4, 2.7 dB) | ≥ 2/3 | pass |
| W2c deciduous − conifer HV amp ≤ 1 dB | −0.2, −1.6, +0.1 dB | ≥ 2/3 | pass |
| W3 forest range over tracks, θ spread 12.7° | deciduous HH 0.5 / HV 0.7 dB, conifer HH 0.9 / HV 0.4 dB | ≤ 1.5 dB | pass |
| W4 frozen − unfrozen HH | deciduous −2.0 / −2.2, conifer −2.6 / −2.3, cropland −1.2 / −1.8 dB (A76 / D113) | ≤ −1.5 on half | pass (4 of 6) |
| W5a \|HV − VH\| ≤ 0.5 dB | forest 16 of 27, cropland 17 of 27 | ≥ 80 % | **fail** |
| W5b cross-pol coherence ≥ 0.9 | forest 10 of 27 (0.76–0.95), cropland 0 of 27 (0.25–0.77) | ≥ 80 % | **fail** |
| W6a forest \|VV − HH\| ≤ 1 dB | 5 of 27 (VV − HH −0.5 … −4.2 dB) | ≥ 80 % | **fail** |
| W6b water VV − HH | no scene passes the NEB gate (no band A water class) | | inconclusive |
| W7 forest H ≥ 0.7 and 40° ≤ ᾱ ≤ 55° | 18 of 22 unfrozen scenes (H 0.65–0.97, ᾱ 48.8–53.1°) | ≥ 80 % | pass |
| W8 Faraday Ω (forest, reported; sign convention in `faraday_deg`) | Oct 2025 – Jan 2026: −1.2 … +1.8°; Jun – Nov 2026: −5.2 … +0.1° | expectation \|Ω\| ≤ 3° | reported |

**Post-hoc (not pre-registered), what the fails and gaps show.**
- *Processing version.* W5a and W6a split by `compositeReleaseId`: X05009/X05010 (Oct 2025 – Jan 2026, 10 scenes) give
  forest HV − VH −0.9 … −1.4 dB (median) and VV − HH −2.8 dB; P05023/P05036 (Jun – Nov 2026, 17 scenes) give HV − VH
  −0.1 … −0.3 dB and VV − HH −0.9 … −1.1 dB. Reciprocity (HV = VH) does not depend on season, so the early releases carry a
  ≈ 1 dB cross-pol channel imbalance that the 2026 releases mostly do not. Season and release are confounded for VV − HH.
- *Cross-pol coherence and noise.* Noise is not subtracted (section K), and receiver noise in HV and VH is independent, so
  the coherence falls where the cross-pol signal is weak: cropland (low HV) 0.25–0.77, forest (high HV) 0.76–0.95. W5b as
  written asked more than an un-denoised product can give.
- *The NEB layer itself.* `noiseEquivalentBackscatter` HV has zeros inside the window in 28 of 52 band B scenes, and on the
  band A grid it is all NaN in 15 of 27 scenes; it is complete only in some P05023/P05036 scenes. Where it is finite, calm
  water HV (about −22 dB, like the Lake Balaton control in section K) sits 3–6 dB above it. Either the stated floor is
  too low or the water returns some HV; we cannot tell which from these data. Add to the K0 audit.
- *Empty water class on the band A grid.* The July NDWI rule (≥ 0.3) removes every lake pixel of the sub-window
  (NDWI there is below 0.3 in July; ≥ 0.3 in November). So W1c and W6b test nothing; the rule was too strict for a
  summer lake, not a radar result.
- *Faraday.* The estimate is near zero in winter and 2–5° in magnitude in summer (negative in the convention of
  `faraday_deg`: M = R(Ω)·S·R(Ω), checked in `experiments/w_selfcheck.py`), largest on the descending (evening) tracks.
  That fits more ionospheric electron content in summer evenings, but crosstalk is not corrected, and no TEC data were
  compared. Not a claim.
- **Reading.** On this clean site the surface physics behaves as textbook L-band: the wave does not see leaves (W2), forest
  γ⁰ is flat in angle (W3), freezing drops backscatter by ~2 dB (W4), forest is a high-entropy volume scatterer at
  ᾱ ≈ 50° (W7). What fails is the product's polarimetric calibration in the early releases (W5a, W6a) and its noise layer
  (W1a), plus our own too-strict water rule (W1c, W6b).
- *Review of PR #6 (2026-11-15), code only, no verdict changed:* W1c now uses the median of 10·log10(HV / NEB) as the
  text says (it was median dB minus median dB); W1a reports broken-layer and finite counts separately; the Faraday sign
  convention is fixed and stated (the W8 values above flipped sign, magnitudes unchanged); `experiments/w_selfcheck.py`
  reruns the synthetic checks of the polarimetry helpers.

## B. Bakota, Dniester canyon: water and terrain geometry with GCOV (`edoTenseii` + Claude)

Goal: a Ukrainian site far from any front, with a drowned village (Bakota, flooded in 1981 by the Dniester reservoir)
and a ~150 m deep canyon. Three questions: does the product's radiometric terrain correction (RTC) flatten slopes that
face the radar; does calm canyon water sit at the stated noise floor; does the radar see the reservoir's water area
change the way optical images do.

**Pre-registered before any GCOV backscatter value over this site was read; this commit is pushed to GitHub before the
data run (the push time is the check).** Seen before writing, and nothing else:
- ASF lists: 48 GCOV scenes (`2005_DHDH`, HH + HV) from 4 tracks (A100/26: 11, A172/26: 11, D94/64: 11, D166/64: 15),
  Oct 2025 – Nov 2026; 39 GUNW pairs. One file's layer list: no layover/shadow layer in GCOV.
- `experiments/b_geometry.py` (DEM + the line-of-sight and incidence cubes only, no backscatter): Copernicus DEM 30 m,
  heights 117–307 m, slope p50 2.5°, p99 30.5°, max 48.4°. Predicted radar shadow and layover cover 0–32 of 130 869
  pixels per track, so a shadow/layover test is **not** run here: the canyon walls are not steep enough for NISAR's
  36–47° incidence. Forested slopes ≥ 15° facing the sensor / facing away: 1 351 / 2 170 (D94), 2 437 / 1 992 (A100),
  1 365 / 2 188 (D166), 2 443 / 2 004 (A172) pixels.
- ESA WorldCover 2021 shares, and the list of clear, full-coverage Sentinel-2 L2A dates (31, Oct 2025 – Nov 2026).

**Window.** W 26.90, S 48.55, E 27.06, N 48.64 (Bakota bay and the canyon around it; the Dniester HPP dam, ~30 km
downstream, is outside and is not studied). GCOV band A (20 m), all 48 scenes; terrain height 260 m above the ellipsoid
for the incidence cube. dB = 10·log10 of gamma0. A window is used for a class only with ≥ 50 class pixels.

**Classes** (radar-independent):
- *forest slope, facing / away*: on the 30 m DEM grid, WorldCover tree fraction ≥ 0.9, slope ≥ 15°, and the slope
  component toward the sensor (per track, from its line-of-sight vector) > tan 15° (*facing*) or < −tan 15° (*away*).
- *water*: Sentinel-2 SCL = 6 for every S2 pixel under the radar pixel on both 2026-08-15 and 2026-03-10, eroded by one
  band A pixel (no NDWI rule: section W showed it can empty a summer water class).

**B1. Terrain flattening.** Dates 2026-06-01 … 2026-09-30. Per track, Δ = median over dates of
[median γ⁰(facing) − median γ⁰(away)] in dB, for HH and for HV.
- B1: |Δ| ≤ 1.0 dB in HH and in HV on ≥ 3 of 4 tracks (RTC removes the brightening of slopes that face the radar).
- Reported, not judged: the same difference for σ⁰ = γ⁰ × `rtcGammaToSigmaFactor`, to show how much RTC changes.

**B2. Noise floor on canyon water** (repeat of K1 / W1 on a Ukrainian site, band A).
- B2a. Water HV / NEB_HV (median of 10·log10 of the per-pixel ratio) within ±2 dB in ≥ 80 % of the scenes where the
  NEB layer is finite and > 0 under all water pixels; fewer than 5 such scenes: inconclusive. The number of scenes with a
  broken NEB layer is reported.
- B2b. Water HH ≥ water HV − 0.5 dB in ≥ 80 % of scenes.

**B3. Does the radar see the reservoir's water area change?** Radar water per scene: band A HH (dB) below a scene
threshold set halfway between the median HH of the *water* class and the median HH of WorldCover cropland
(fraction ≥ 0.9) in the same scene. Optical water: S2 SCL = 6 on the clear S2 date nearest the radar date, used only
if within ±5 days. Both as area (km²) over the whole window on the 20 m grid.
- B3a. Radar area / optical area within 0.85–1.15 in ≥ 70 % of the matched pairs (≥ 6 pairs, else inconclusive).
- B3b. Spearman ρ between the radar and the optical areas over the matched pairs ≥ 0.5 (the radar follows the same
  rise and fall of the water area).

**Results B (2026-11-15, `experiments/b_fetch.py`, `experiments/b_tests.py`; pre-registration pushed at 13:54 before
the data run).** 47 band A windows (one SHSH scene left out: no HV), all on one 20 m grid: A100 10, A172 11, D94 11,
D166 15; releases P05023/P05036, plus 5 X05009 scenes on D166. 16 Sentinel-2 dates matched a radar date within ±5 days.

| test | value | limit | verdict |
|---|---|---|---|
| B1 facing − away γ⁰, HH | 0.79, 1.25, 0.42, 1.89 dB (D94, A100, D166, A172) | \|Δ\| ≤ 1 | 2 of 4 |
| B1 facing − away γ⁰, HV | 3.07, 2.29, 2.47, 2.79 dB | \|Δ\| ≤ 1 | 0 of 4 |
| B1 (both pols, ≥ 3 of 4 tracks) | | | **fail** |
| B2a water HV / NEB_HV within ±2 dB | 43 of 47 (median +1.1 dB; outside: +2.2, +3.8, +4.7, +10.0); NEB layer complete under the water in all 47 | ≥ 80 % | **pass** |
| B2b water HH ≥ HV − 0.5 dB | 34 of 47 | ≥ 80 % | **fail** |
| B3a radar / optical water area within 0.85–1.15 | 0 of 33 (ratios 1.2–1.7 on 17 pairs, 6.6–9.0 on 13, 3.4 once, ~60 twice) | ≥ 70 % | **fail** |
| B3b Spearman radar vs optical area | −0.20 (33 pairs) | ≥ 0.5 | **fail** |

Reported (B1): σ⁰ = γ⁰ × `rtcGammaToSigmaFactor`, facing − away: HH −1.5 … −2.4 dB, HV −0.8 … +0.1 dB.

**Post-hoc (not pre-registered).**
- *B1 sign check (after the results, no change made).* The sides are as intended: NISAR looks **left** (`lookDirection`)
  and the product defines its line-of-sight vector from the target to the sensor, which is how `b_geometry.py` used it.
  So, after RTC, forested slopes that face the radar stay brighter than slopes that face away, by ~2.3–3.1 dB in HV and
  0.4–1.9 dB in HH. The terrain correction under-corrects these slopes, more in the volume channel. The σ⁰ numbers are
  reported only: we did not check what `rtcGammaToSigmaFactor` normalises to, so we draw nothing from their sign.
- *B2: the first site where the stated noise floor holds.* Across the 47 scenes the product NEB (HV) ranges from
  −22.9 to −10.9 dB (it changes by ~12 dB between tracks and dates), and calm canyon water follows it within ±2 dB in 43.
  Compare Wisconsin (W1a: layer often broken, water 3–5 dB above where finite) and the Balaton control (+2.2 … +3.2 dB).
  B2b fails because the water HH also sits near the floor: there is little wind-driven Bragg return in the canyon.
- *B3: the optical reference is not stable.* Sentinel-2 SCL water in the window jumps between ~19 km² (close to
  WorldCover's ~17 km²) and ~3.6 km² (once 0.5 km²) from one date to the next, which the reservoir cannot do; the SCL
  class drops most of the water on those dates (likely turbid or algae-rich summer water, not checked). The radar area is
  steady at 24–33 km² (one outlier, 64 km², on a date where water HH and cropland HH were equal). On the 17 pairs with a
  plausible optical area, the radar threshold gives 1.2–1.7 × more water: dark, smooth ground and mixed shore pixels
  pass the HH threshold. B3 as designed cannot show a water-level "breathing"; it would need a stable optical water mask
  (e.g. NDWI with a per-date threshold) or gauge data.

## TW. Do trees dry out by evening? Forest HV against the air's drying power (`edoTenseii` + Claude)

**Idea.** Two natural processes and a link between them. (1) On a hot, dry afternoon the air pulls water out of leaves
faster than roots refill the stems: the vapour pressure deficit (VPD) drives transpiration, and tree water content falls
through the day and recovers at night. (2) Water sets the permittivity of wood (ε of fresh wood ~ 20 – 30 at 1.26 GHz,
of dry wood ~ 2 – 4), and the L-band HV return of a forest is volume scattering from branches and stems. So drier wood
should give a weaker HV echo. NISAR passes this latitude at dawn (ascending, ~05:00 local solar time) and at dusk
(descending, ~19:05), so it sees the same forest after the night's recovery and after the day's drying. Tower radars
(e.g. BorealScat in Sweden) have reported diurnal L-band forest backscatter changes linked to tree water; we test whether
a satellite with a 12-day repeat can see it.

**Pre-registered 2026-11-15 ~15:50 (Dnipro time), committed and pushed on its own before any value used below is computed.**
Seen before writing:
- Wisconsin: the per-date median HV of the deciduous class (section W, test W2) and the 2 m temperature at each
  overpass (`temps.json`, frozen / unfrozen rule), but never set against VPD, humidity or rain.
- Bakota: per-date HV medians of forest on slopes ≥ 15° (test B1); forest on flat ground was never computed.
- Scene metadata only (track, direction, UTC time, release): every scene in the window below is release P05023.
- Pixel counts of the Bakota forest class (from WorldCover and the DEM, no radar value).
- No VPD or rain value for any site has been read.

**Data.**
- *Wisconsin*: band B 80 m, full window of section W, class *deciduous forest* exactly as in W.
  Evening tracks D41, D113; morning track A76 (A148 has no scene in the window).
- *Bakota*: band A 20 m, window of section B. Class *forest*: WorldCover tree fraction ≥ 0.9, DEM slope < 10° (keeps the
  terrain-correction error of B1 out), not in the radar shadow or layover of that track. Evening tracks D94, D166;
  morning tracks A100, A172.
- Valid pixels: `mask` in 1 … 254 and HV > 0. Class value per scene = median of 10·log10(γ⁰ HV); fewer than 50 pixels
  drops the scene.
- *Season*: local solar dates 15 Jun – 30 Sep 2026 (leaf-on, unfrozen). Local solar time = UTC + longitude / 15°.
  Wisconsin evening passes are at ~01:06 – 01:14 UTC, so their local date is the UTC date minus one day.
- *Weather*: Open-Meteo historical archive, hourly at the window centre: `vapour_pressure_deficit` (kPa) and
  `precipitation` (mm).

**Variables.**
- *Afternoon VPD* of a pass: mean hourly VPD over 12:00 – 18:00 local solar time of the latest afternoon that ended
  before the pass (same day for an evening pass, previous day for a morning pass).
- *Dry pass*: precipitation over the 24 h before the pass < 1.0 mm (a wet canopy adds water that is not in the wood).
  Only dry passes are used in the judged tests.
- *Anomalies*: within each (site, track), the dry-pass HV values and their afternoon VPD values are each detrended with a
  straight line in time (removes slow phenology, growth and the seasonal VPD cycle). Each track is compared only with
  itself, so incidence angle and look side cancel.

**Tests.**
- **TW1 (main).** Evening passes, both sites pooled: Spearman ρ between the HV anomaly and the afternoon-VPD anomaly.
  Pass: ρ ≤ −0.3 and one-sided permutation p ≤ 0.05 (10 000 permutations of VPD within each track). Otherwise fail.
  Fewer than 15 dry evening passes: inconclusive.
- **TW2 (replication on two continents).** Evening passes, the same ρ at each site on its own. Pass: ρ < 0 at both
  sites. Fail: ρ ≥ 0 at either. A site with fewer than 6 dry evening passes: inconclusive.
- **TW3 (dawn against dusk).** Morning passes, pooled ρ_m by the same recipe. Pass: ρ_m ≥ ρ_e + 0.2 (the night refills
  the trees, so the morning link is weaker). Judged only if TW1 passes; otherwise reported.
- **Reported, not judged:** slope of HV anomaly on VPD anomaly (dB per kPa, evening, pooled, with a bootstrap 90 %
  interval); the same tests with centring only (no linear trend), with a 48 h rain rule, for HH, and for the Wisconsin
  conifer class; a placebo: the evening HV anomaly against the VPD of the afternoon *after* the pass (expected weaker
  than TW1; it is not zero because weather persists for days).

**What would go wrong.** The expected change is small: a few % of wood moisture moves ε by ~10 %, which moves a
branch-scattering HV by roughly 0.5 – 1 dB, while the per-date class median scatters by ~0.3 dB plus calibration drift.
In a dense canopy the echo depends on the scattering-to-extinction ratio, which is less sensitive to ε than either
term alone. A fail is a usable result: it bounds what NISAR can see of tree water from orbit.

**Timing note.** The pre-registration text says "~15:50"; its commit (ebab20e) is stamped 15:43 Dnipro time and was pushed
right after, before `experiments/tw_tests.py` existed. The first run of the script returned no passes at all (an hour key
kept fractional microseconds, so no VPD value matched); the fix touched only that key, no rule or threshold.

**Results TW (2026-11-15 ~16:10, `experiments/tw_tests.py`).** 53 passes in season, 38 dry (24 h rule): evening 21
(Wisconsin D41 4, D113 5; Bakota D94 6, D166 6), morning 17 (Wisconsin A76 5; Bakota A100 5, A172 7). All release P05023.
Class pixels per scene: Wisconsin deciduous 15 684 (80 m), Bakota flat forest ≥ 7 291 (20 m). Afternoon VPD ranged
0.38 – 2.62 kPa; its detrended anomalies ±0.6 kPa.

| test | value | limit | verdict |
|---|---|---|---|
| TW1 evening, pooled ρ (HV anomaly, VPD anomaly) | −0.11, permutation p = 0.34, n = 21 | ρ ≤ −0.3 and p ≤ 0.05 | **fail** |
| TW2 evening ρ per site | Wisconsin +0.10 (n = 9), Bakota −0.08 (n = 12) | both < 0 | **fail** |
| TW3 morning against evening | ρ_m = −0.31 (p = 0.13, n = 17) | judged only if TW1 passes | not judged |

Reported:
- Slope, evening pooled: +0.03 dB/kPa, bootstrap 90 % interval −0.25 … +0.43 dB/kPa.
- Centring only (no trend): evening ρ +0.09, morning +0.15.
- 48 h rain rule: evening n = 5 (too few), morning ρ −0.44 (p = 0.07).
- HH: evening +0.12, morning −0.09.
- Wisconsin conifer HV: evening +0.17 (n = 9).
- Placebo (VPD of the afternoon after the pass): evening −0.04, morning −0.07.

**Post-hoc (not pre-registered).**
- *What the fail bounds.* Wisconsin evening alone: the deciduous HV anomaly scatters by only 0.11 dB (s.d.) while the
  afternoon VPD anomaly spans ±0.5 kPa. Slope +0.03 dB/kPa, bootstrap 90 % −0.15 … +0.33. So a drier-than-usual
  afternoon lowers the dusk HV of this forest by less than ~0.15 dB per kPa (one-sided 95 %). That is several times below
  the 0.5 – 1 dB we expected from wood permittivity. Two readings, not tested:
  1. The day-to-day change in wood water at dusk is smaller than a few %, because stems buffer it.
  2. Most of the daily water swing is in the leaves and thin twigs, which the 24 cm wave does not see (W2).
- *This is not the dawn-to-dusk swing itself.* Morning and evening tracks look from opposite sides at different angles,
  so the test compares dusk with dusk only. A tower or a same-geometry morning/evening pair would be needed for the
  full diurnal cycle.
- *Bakota evening tracks sit close to the stated noise floor.* The flat-forest HV on D94 is about −9 dB with NEB_HV
  −12 … −13 dB, so the floor is only 3 – 4 dB below the signal. On D166 the gap is ~6 dB, on the morning tracks ~11 dB
  (NEB from section B). Noise not subtracted in GCOV can shift those medians. The HV anomaly does not follow the NEB
  across dates (Spearman +0.15 evening, +0.01 morning).
- *Two Bakota outliers.* D166 on 17 and 29 Sep: +1.18 dB, then −1.14 dB (forest HV −13.4 dB on 29 Sep, 2 dB below the
  other dates). Both passes were dry by the 24 h rule. Cause not checked; it may be the start of leaf fall or a
  processing change. They do not decide TW2: Wisconsin alone already has ρ > 0.
- **Reading.** From orbit, at a 12-day repeat, NISAR does not see an afternoon's dryness in the dusk echo of a temperate
  forest above ~0.15 dB/kPa. The stability is itself useful: summer dusk HV of a deciduous forest repeats within ~0.1 dB
  in Wisconsin, a usable reference target for calibration checks.
