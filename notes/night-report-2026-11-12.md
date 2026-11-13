# Night report, night of 10→11.11.2026 (queue completed, report closed ~04:05; the file name with the date 12.11 is wrong)

Briefly, in simple words: what we checked, what came out against the pre-registered predictions, what surprised us, and what you need to decide.

## The main things from the night

**Velocity maps (Mexico City, valley, Houston).** The main enemy is errors tied to *dates* (atmosphere and the ionosphere correction for each pair). In a full chain of pairs they cancel out, but:
- **gaps in the chain "break" the cancellation** (71% of Houston desc pixels have gaps; without them the error goes 28 → 2.8 cm/year);
- in Houston desc these errors are also **huge and shared across a large area** (±0.5–1.5 m on almost every date, the same at stations tens of km apart).
- Three cures were checked, and none passed the threshold:
  - a robust median of slopes (H13) is worse, because bad dates are the majority;
  - tying each pair to GNSS with a plane (H14) cut down the monsters (94 → 15 cm/year), but a plane is too rigid;
  - smoothing the screen (H10, H11) gives almost nothing.
- Conclusion: into the app goes an **honest trust layer** (number of chain breaks + a "tied to GNSS" label), and the real cures (a flexible per-date error field from GNSS or from the radar itself, ML-2) go to the hackathon.

**Earthquakes.**
- The recipe is frozen (Q-2look v1): network + exact slip solution from both passes, **0.5 s** instead of 13 s; "polishing" adds nothing.
- On synthetics it is right in only ~40% of cases, **but it knows itself when it is right**: agreement among the 16 best proposals predicts success with AUC 0.93; at agreement ≥ 0.9, 14 of 14 are correct.
- Panama is a low-confidence case (strike stable, the rest not), and the flag would have shown it.

**Fields.** Field coherence in the valley has a clear yearly rhythm (winter 0.87, spring-summer 0.61–0.73); even a group of fields that "dances in reverse" is visible. A cheap, reliable layer that fits the challenge name.

**Debris flow in Nepal.** The phase does not see a precursor on the ice wall; the debris-flow path is clearly visible: a ready demo layer.

**Mathematics.** How many fault parameters a small radar window sees is a count of monomials (steps of m + 1). It looks like a known theorem about kernel matrices in the "flat limit". The best candidate for "space → statement" (math.md M9, addendum 3).

## What was checked (verdicts against the pre-registrations)

| Test | Question | Verdict | What it means for the app |
|---|---|---|---|
| D-hyb | Earthquake: the network proposes, 60 fit steps refine | Kumamoto passed, Tibet (Nagqu) did not: strike 19° from W-phase, but only 2° from the USGS finite fault | "Earthquake" mode is real, ~7 seconds per event |
| P-port, valley | Does the Mexico City pipeline transfer to another place | Failed on accuracy (RMS 3.6–3.8 cm/year with a limit of 2.5). **Passed:** in the dry valley, evening is not worse than morning (0.81) | "Morning is cleaner" is about climate, not about the clock. Automatic pass choice must look at climate |
| P-port, Houston | The same | Failed (RMS 12 and 95 cm/year) | Cause found (below), not a pipeline bug |
| H9 | Does the ionosphere "leak" through wet fields | Failed (a, b). (c) passed, but mixed with the season | The ionosphere error is large-scale (tens of km), not patchy from fields |
| H9-P | Take only the slope of the ionosphere correction | Failed by the letter. But Houston 11.9 → 5.7 and 95 → 30, valley asc 3.8 → 2.8 | A useful option, not a full cure |
| Trust module | Is σ honest at GNSS stations | Checked on Mexico City (reproduces). In the valley, failed: the σ estimator "drops" to zero | The trust layer works where the noise comes from dates. Where the noise comes from pairs, a different formula is needed |
| P-tie | Tie the map to GNSS with a plane, check with a held-out station | Failed everywhere (valley 2.95 / 3.75 with a limit of 2.0, Houston 10 / 103) | The "tied to GNSS" label stays, but it is not a cure: the errors come in islands |
| H10 | Ionosphere screen smoothed with a 10 km median | Failed everywhere, even worse than "slope only" (Houston 7.96 / 87, valley 2.90 / 3.82) | A spatial filter does not cure it. We need the ionosphere per date (P-iono2) or ML-2 |
| D-3 Panama | Third blind earthquake, M7.7 | Failed: the network gave reverse 78/30/111, Mw 7.28; USGS gives strike-slip 288/66/−17, Mw 7.7. Least squares got Mw right (7.66), but not the strike | "Earthquake" mode must choose the window and the prior by magnitude. G3: 1 of 3 |
| L-retro | Did NISAR see slope motion before the collapse in Nepal | L1, L2 not measurable (coherence at rock bottom since April). L3 "passed", but because of ice, not because of a precursor | Monitoring icy slopes by phase is impossible. Next: pixel offsets (pixelOffsets) |
| C | What subsidence bowl is visible on the Mexico City map (we inject artificial bowls) | Failed by the letter: the 50% threshold ≈ 5–6 cm/year for 0.5–1 km bowls. Wider ones (2–4 km) are seen worse, not better. False alarms ≤ 2% | A matched filter (C′) gave exactly the same thresholds, so the limit is in the data, not in the detector (my post-hoc claim about "the lost half of the signal" was wrong). For the app: "confirm slow, wide subsidence with GNSS; NISAR sees sharp local bowls on its own" |
| D-geom | How a fault "disappears" with depth (mathematics) | Failed by the letter, but two tiers of directions were found (4 and 5) with a gap ≈ multipolar | A scientific bonus, not for the demo (math.md M9) |
| H11 | Robust plane of the ionosphere screen (islands do not skew the slope) | Formally passed, but the gain is tiny (Houston desc 29.7 → 27.8) | Islands are not the main problem; the main one is gaps (see above) |
| H12 | Only pixels without gaps in the chain | Valley desc 3.72 → 3.04 (passed). Mexico City: passed by the letter, worse in substance (my mistake in the reference) | Do not discard them, account for them in trust |
| P5b | Will the trust layer become honest if gaps are accounted for (M2b) | Failed: within bounds only 1 of 4 (3 needed). The direction is right: valley asc 46 → 54%, Houston desc 9 → 21% | Gaps are a correct but small addition. Show them as a separate layer "how many times the chain breaks" |
| X-coh | Do the valley fields "breathe" with the calendar (field coherence vs the city) | X1 and X2 passed (summer 0.68 vs winter 0.87; the curve is smooth, autocorrelation 0.85). X3 (rain) failed | A cheap, reliable "field rhythm" layer straight from coherence, without unwrapping or ionosphere. Post-hoc: there is a group of fields that "dances in reverse" |
| D-jet | Does the number of Fisher "tiers" equal the ranks of the jet of derivatives (mathematics) | k = 0 and 1 matched exactly, k = 2 only on even smaller windows (failed by the letter). Rank 9 is reached (passed) | How many fault parameters a window sees is a count of monomials. A theorem candidate (below) |
| D-2look-speed | Freeze the earthquake recipe and speed it up, on synthetics only | S1 passed: the chosen variant 0.54 s. S2 failed: "polishing" adds nothing | Recipe Q-2look v1 = network + exact slip solution, no fitting. On Panama 0.62 s, strike 2°, dip 11°. But on synthetics, success only ~40%, and the Panama rake is wrong |
| D-conf | Does the earthquake card know when it is wrong | F1 passed (AUC 0.93), F3 passed (0.91), rake on success is accurate (3.5°). F2 failed by a hair (0.67 vs a threshold of 0.70 for agreement ≥ 0.5) | The confidence flag works: at agreement ≥ 0.75 success is 80%, at ≥ 0.9 it is 100% (post-hoc thresholds; freeze ≥ 0.75 for the next blind test) |
| H13 | Robust velocity (median of slopes over chain pieces) | Failed (Houston desc 29.7 → 46, valley desc 3.7 → 4.5); Mexico City unchanged | Houston errors are not isolated but shared across the region on almost every date |
| H14 | Tie each pair to GNSS with a plane (check with a held-out station) | Failed 0 of 3. But Houston desc 17.8 → 6.3, worst station 94 → 15 | The per-date error field is smooth, but not a plane. Next: a more flexible surface (spline / kriging) or ML-2 |
| D-ident | What a single radar image can say about a fault | Fault dip is seen best; position along strike / depth worst. Normal vs strike-slip differ by a factor of 2–3 with one pass and 9–11 with two | The mechanism error at Nagqu was a network error, not a limit of the data |

## The most important discovery of the second half of the night: gaps break the telescoping

- In Houston desc, almost all the error comes from one station (TXET). There, date after date, there are jumps of ±0.5–1.5 m, but they are "nodal" (per date), so in a full chain they cancel out.
- The problem is that **71% of pixels have gaps** (at least one invalid pair). Each gap breaks the chain into pieces, and the noise of the end dates of each piece leaks into the velocity.
- If we take only pixels without gaps, the Houston desc error drops from **27.8 to 2.8 cm/year** (post-hoc). For valley desc, checked in advance: 3.7 → 3.0. In Mexico City it got worse at one station, because the nearest full pixels lie on the steep slope of the bowl.
- The math is simple (math.md M2b): node noise gives a variance of 2s·σ²/T², where s is the number of pieces.
- **For the app:** do not discard such pixels; carry the number of pieces into the trust layer.

## Second half of the night: three more things worth knowing

**1. The earthquake recipe became 20 times faster, but we must be honest about accuracy.**
- On 24 synthetic M7–7.9 earthquakes (with real Panama noise), 60 "polishing" steps gave nothing. The network plus the exact slip solution give the same in half a second.
- The recipe is frozen as Q-2look v1; the next event will be a blind test of exactly this recipe.
- Bad news: on synthetics, full success is only ~40% of cases. The errors are "either accurate or completely off"; some are 180° flips for nearly vertical faults, which is almost the same plane.
- On Panama, strike and dip are correct, but **the rake is wrong** (+57 instead of −17), so the "fault type" label would have been wrong.
- For the app: the earthquake card must show a confidence flag, not a single answer.

**2. The valley fields have a rhythm, and it is visible without any complex corrections.**
- Field coherence relative to the city: ~0.9 in winter, drops in steps in March–April to ~0.61, partly recovers in summer.
- Rain has almost no effect on it (L-band).
- One group of fields (98 blocks of 2 km) "dances in reverse": worst in March, best in September. Possibly winter crops vs irrigated summer crops and orchards; not checked against the USDA crop map.
- A ready "field choreography" layer, literally fitting the challenge name. Plot: `data/processed/central-valley-asc/x_coh.png`.

**3. Mathematics: a theorem candidate in the spirit of "space → statement".**
- How many combinations of fault parameters a small radar window can determine is a **count of monomials**.
- Each order of derivatives of the field has m + 1 coefficients, and the 9 parameters fill them like water fills steps: with the plane removed, 3 combinations are visible at the level of map curvature, 4 more at the next order, the last 2 one order higher.
- Each order costs a factor of (ℓ/d)² in information.
- Exactly this grouping of eigenvalues is known for kernel matrices in the "flat limit" (radial basis function interpolation: Fornberg, Wathen–Zhu, Barthelmé–Usevich; names from memory, need checking).
- So this is probably a consequence of known theory, read in a new place. Details in math.md (M9, addendum 3).
- Practical output: a rule "how many numbers about the fault this image actually determines" for the earthquake card.

**The "how many times the chain breaks" layer** is now written as the third channel in `data/processed/trust/*-seg/sigma.tif` together with `segments.png`. Full pixels (chain without breaks): Houston desc 29%, Houston asc 54%, valley asc 44%, valley desc 38%.

## Surprises

1. **The islands in Houston sit in the ionosphere screen, not in the phase.** Relative to station ADKS: phase ±0.8 cycle, ionosphere screen from −10 to +15 "cycles". This changes the goal of ML-2: look for jumps in the screen too.
2. **Fake winter uplift of +10 cm in the valley** from the same correction. The radar cannot see this by itself, because the error has a single sign.
3. **Dry climate = evening is not worse.** The "morning is cleaner" hypothesis was confirmed as a climate effect, not a time-of-day effect.

## New directions (written in notes/research-directions.md)

- **P-iono2 (now the main cure candidate):** the ionosphere per date over the pair network, not per pair.
- **L-dpm:** a map of the debris-flow or landslide scar from coherence loss. Post-hoc on Chhochen Khola: a 4.8 km² "tongue", about 6 km down the valley to the northwest, from ~4000 to ~2650 m. On two pairs before the event there are no false alarms; on one autumn pair (September) there is a false alarm, so a seasonal safeguard is needed.
- **L-offsets:** a precursor via pixel offsets (offset tracking), not phase; the phase does not catch fast motion on ice.
- **D-geom (mathematics, M9), now with a result:**
  - The Fisher spectrum does split into two tiers (4 and 5 directions), the gap ≈ multipolar.
  - But the fault position falls into the sloppy tier. The first explanation (plane removal as "calibration") was confirmed only partly.
  - A better explanation: a small window sees the Taylor series of the far field. Hence the hypothesis that the number of directions in each tier is the rank of the map "parameters → field derivatives" at each order.
  - A neighboring known field is Euler deconvolution in gravity and magnetic prospecting; the literature needs checking.
  - This is a candidate for an exact statement in the spirit of your "space → theorem".
- **Earlier version of this item (before the result):** the fault as a "sloppy model" (Sethna / Transtrum).
  - Hypothesis: the boundaries of the Okada model manifold are multipole approximations (point source, line source). Then the "sloppy" eigenvalues should decay as a power of the size/depth ratio. Can be checked in an hour.
  - Link to your idea "learned space → theorem": do the network's uncertainty directions match the Fisher directions? If so, a simplified model can be read from the network's latent space.
- **L-retro (debris flow / collapse, Nepal):** ongoing, see below.

## Debris flow and collapse at Chhochen Khola (Nepal, 26.08.2026)

- What is known: a collapse of ice and rock from the north wall of Langtang Lirung (5000–5500 m), then a debris flow down the Lhende Khola.
- ICIMOD writes that, according to radar data, the slope moved up to 30 cm in the 3–4 weeks before the collapse. So our test is not fully blind: we know the precursor was present in *some* radar data.
- Our question: does NISAR GUNW at 80 m see it, and could it be a "slope monitoring" layer in the app.
- Pre-registration done before the download: L1 (displacement ≥ 3σ), L2 (acceleration), L3 (coherence drop).
- Data: 25 pairs of the ascending pass (November 2025 to October 2026). The last pair before the event is 26.07 → 19.08 (24 days).
- Result:
  - The source zone was found at an altitude of ~5270 m, ~3 km north of the summit, but it is tiny (0.13 km²).
  - Coherence there has been at rock bottom (~0.09) since mid-April, 4 months before the collapse, so motion cannot be measured. Most likely this is just summer ice melt.
  - Post-hoc: the debris-flow path in the pair 19.08 → 31.08 is very clearly visible (see `data/processed/landslide/langtang_lretro.png`).
  - **Pixel offsets (L-offsets):**
  - By the registered rule the test is not measurable: correlation quality at the source is low, and the filter rejected everything.
  - Post-hoc without the filter: the ice at the source moved steadily, ~15 cm per 12 days, all winter. This is simply a hanging glacier flowing.
  - In the 24 days before the collapse it moved 37 cm. This matches ICIMOD's "up to 30 cm", **but it is not faster than in winter**. In our data the "precursor" does not differ from ordinary ice motion.
  - Important: on rock the offset accuracy is 2–7 cm, so a velocity map of fast slopes and ice from GUNW is real. There is no alarm signal for acceleration yet.
- Honest conclusion for the app: "where the debris flow went": yes, within 12 days. "Warn in advance" on an ice wall: no, at least not from GUNW phase.

## Earthquake in Panama (third blind test)

- The frozen results were recorded before unblinding.
- Hybrid H1: reverse, strike 78°, dip 30°, Mw 7.28, ~7 s.
- USGS: strike-slip, plane 288/66/−17 (finite fault 295/60/−16, three segments, ~100 km), Mw 7.7.
- Failed on strike, dip and magnitude. The location was right.
- Cause (post-hoc, but named as a risk in advance): the rupture is longer than the ±55 km window and the prior L ≤ 80 km.
- For the product: the window size and the length prior must depend on the catalog magnitude.
- **Diagnostics with a window twice as wide (D-3w, no longer blind).** The hybrid found the correct type (strike-slip), strike 12° from USGS, Mw 7.46. But the dip is absurd (5° instead of 75°), so by the rule it is still a failure.
- Conclusion: the window was part of the problem. Dip is not determined from one pass with 33% valid pixels. Next step for "earthquake" mode: ascending + descending pass together (D-2look); for Panama, descending pairs exist.
- **D-2look-lite (done during the night, also not blind): both passes together.**
  - The two-pass hybrid gave strike 290° vs 288° from USGS, dip 55° vs 66°, Mw 7.57 vs 7.7. Four of four criteria passed. Within 5° of the USGS finite fault (295/60).
  - Caveats: 13.7 s instead of < 10 s; rake off by ~50°; the test is not blind.
  - ~~This is the strongest result of the night for "earthquake" mode.~~ **Correction (D-conf):** over 12 random runs the strike is stable (9 of 12: 280–299°), but dip wanders 29–70°, Mw 7.3–7.6, rake wrong; agreement among proposals is low (~0.15). So "4 of 4" was partly luck. The recipe is frozen as Q-2look v1 (no fitting, 0.5 s) together with the confidence flag.

## What you need to decide

1. **Product decision for velocity maps.**
   - I propose that the app shows velocity together with two trust layers:
     - "how many times the chain breaks" (done, third channel of `sigma.tif`);
     - the "tied to GNSS / radar only" label.
   - Honestly mark Houston desc as "unreliable". Agree?
2. **Earthquakes: frozen recipe + flag.**
   - The next blind test is Q-2look v1 with an agreement threshold ≥ 0.75 (if lower, the card says "not confident").
   - We need to pick an event: M ≥ 6.5 with both NISAR passes. Search now, or wait for the official resources on 13.11?
3. **"Field rhythm" layer (X-coh).** Include it in the demo? Checking the field groups needs the USDA CDL crop map. It is not on the night's list of allowed sources, so your OK for the download is needed (public, no login).
4. **Debris flow: "where the debris flow went" layer** (coherence loss): include it in the demo? There is no precursor from phase; we should say so honestly.
5. **Questions for the mentors on Friday.**
   - Will there be GUNW pairs with a longer baseline (24/36 days) to compute the ionosphere per date (P-iono2)?
   - Why does Houston desc have meter-level errors per date?
6. **Mathematics (D-jet).** It is worth spending an hour checking the literature (flat-limit kernel spectra: Wathen–Zhu, Barthelmé–Usevich; Euler deconvolution) before calling this new. Do it now or after the hackathon?
