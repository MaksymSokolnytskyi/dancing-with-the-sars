# End goal: what the user sees, what ML gives, and which experiments serve it

Status as of 11.11.2026, night. A living document: every new experiment must point to one of the goals G1–G4 below, otherwise we do not do it.

## The app in one sentence

**"Honest map of ground motion":** NISAR shows where the ground sinks or shifts, drawn with engraving strokes (hachures), and every stroke shows how far it can be trusted.

## What the user sees (demo)

1. **Choosing a place:** Mexico City, the California valley, Houston (or any city with NISAR coverage).
2. **Hachure map:** red strokes run down into subsidence bowls; thickness shows the rate.
   - Time slider: the strokes grow date by date (the animation already exists for Mexico City).
3. **Trust layer**, on by default:
   - a stroke where the signal is smaller than the noise is not drawn;
   - a borderline stroke is drawn pale or dashed;
   - a quality label next to the map: "tied to N GNSS stations" or "radar only: error ±X cm/year".
4. **"Why so" panel:**
   - the best pass for this place (morning or evening, depending on climate);
   - bad dates (thunderstorm per GOES, a broken pair);
   - warnings about known traps (ionosphere in winter, unwrapping islands).
5. **"Earthquake" mode:** list of recent M6+ → click → a fault card in ~7 seconds.
   - Fault rectangle on the map, direction, dip, type (strike-slip, normal, reverse), magnitude, with uncertainty.
   - Next to it, hachures of the coseismic displacement.

## What from the night's work actually changes the demo

| Finding | How it affects the app | Novelty (honestly) |
|---|---|---|
| Hachures instead of colored fringes | The main visual hook | High: we have not seen InSAR in the Dufour style |
| Trust layer in the strokes themselves (they disappear or fade by SNR) | The user does not see "flowers" in the mountains that are only noise | High for visualization. The σ method itself is known in parts |
| Per-date noise from the pair network itself (graph, lemma M1.1) + noise share of a pair | A "bad dates" list without any external data | Medium: this follows from known network tomography, but we did not find it in this form for InSAR chains |
| **The radar cannot check itself against a same-sign drift** (valley: +10 cm of fake winter uplift from the ionosphere correction) | The "tied to GNSS / radar only" label is mandatory. Without it the map lies confidently | A practical finding for NISAR GUNW users: the product is new, there are few such reports yet |
| The ionosphere correction, computed separately for each pair, accumulates error along the chain (Mexico City H8, valley) | A "correction slope only" option. Next step: estimate the ionosphere per date | The idea itself is known (Fattahi 2017), but measuring this effect on NISAR is new |
| The islands in Houston sit **in the GUNW ionosphere screen**, not in the phase (found out during the night) | Without detecting them the Houston map is unusable. Smoothing the screen (H10) and tying to GNSS (P-tie) do not help. We need the ionosphere per date, or ML-2 | High for NISAR: the product is new, the flaw is in the product itself |
| The debris-flow path is visible as coherence loss over a single 12-day pass (Nepal, post-hoc) | A "where the debris flow / landslide went" layer | Low (known damage proxy map), but clear and honest |
| Earthquake in 7 seconds: the network gives the geometry, algebra gives the slip, a short fit gives the precision. Blind tests: 1 of 3 (Panama failed because the rupture is longer than the window) | "Earthquake" mode in the demo, but the window and the prior must depend on magnitude | Medium: SBI for faults exists. The novelty is in combining it with real NISAR noise and an exact solution for the slip (M8) |
| "Morning is cleaner" depends on climate (yes for Mexico City, still checking elsewhere) | Automatic pass choice for each city | Small, but useful and clear |

## ML: inputs → outputs

| ID | Why (goal) | Input | Output | Training | Status |
|---|---|---|---|---|---|
| **ML-1 Fault in seconds** | G3 | Coseismic LOS map on a ~1 km grid (2 km for large events) + mask + coordinates + look angle | Posterior of the fault geometry (7 parameters) → exact solution for the slip → 60 fit steps | Okada synthetics + real NISAR noise from pre-event pairs + turbulence | Works. Kumamoto passed, Tibet borderline, Panama in progress |
| **ML-2 Unwrapping island detector** | G1 | Per pixel: unwrapped phase of the pair, coherence, connected components, ionosphere, phase of neighboring pairs (chain context), closure triangles where available | Mask "pixel shifted by k cycles" (classification k ∈ {−2..2}) or a correction | Inject artificial whole-cycle jumps into real clean pairs (Mexico City morning, valley) + labels from closures and GNSS | New. Needed, because Houston is unusable without it. **Update after the night:** the jumps sit in the ionosphere screen, so the input must contain the screen separately from the phase, and the output must be a mask "the screen is broken here" |
| **ML-3 Detection threshold map** (C) | G2 | Real pairs without an event + a synthetic bowl of given size and rate | Probability of detecting a bowl of X cm/year and size Y km at this place | Injection-recovery | Planned. Gives "what this map can see" |
| ML-4 Cleaning a pair of atmosphere (optional) | G1 | One pair + GOES water vapor + DEM | Pair without turbulence | Kolmogorov synthetics + real noise | Idea. Least important for the demo |

Everything else (trust layer, GNSS tie, pass choice) is classical, without ML, and should stay that way: there, being right matters more than being clever.

## Goals and experiments

- **G1. The map is correct** (accuracy against GNSS in ≥ 2 cities).
  - Done: H8, H8b, H9-P, P-port.
  - During the night: P-tie failed, H10 failed, CV P3 passed. The best simple option is ionosphere "slope only" (valley asc 2.78).
  - Next: ionosphere per date (P-iono2), ML-2 for Houston.
  - "Done" criterion: RMS against GNSS ≤ 2.5 cm/year in the valley and ≤ 1 in Houston, with held-out stations.
- **G2. Trust is honest** (|z| < 1 for ~68% of stations).
  - Done: V1, P-trust, trust module.
  - Next: σ = radar + GNSS tie; ML-3.
  - "Done" criterion: calibration 55–85% in two cities.
- **G3. Earthquake in seconds.**
  - Done: D, D-hyb.
  - During the night: Panama failed (1 of 3). Next: window and prior by magnitude, a new blind test; D-auto (card from a USGS event ID).
  - "Done" criterion: 2 of 3 events pass the blind test, time < 10 s.
- **G4. Beautiful and clear.**
  - Hachures, trust hachures, animation, "why so" panel.
  - "Done" criterion: an untrained person says within 30 seconds where it is sinking and whether to trust it.
- **Science "for later"** (outside the demo, but in the document): D-ident / sloppy models, T1 radar as meteorology, mathematics of the noise graph.

## Priorities before the hackathon starts (Friday 13.11)

1. G1 + G2 on Mexico City and the valley: GNSS tie, "radar only / tied" label.
2. G3: Panama + the earthquake card as a separate module.
3. Start ML-2 on Mexico City and the valley (the truth is known there), apply it to Houston.
4. Write the demo scenario for the team.
