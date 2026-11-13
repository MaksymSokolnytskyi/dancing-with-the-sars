# Night plan (autonomous), 2026-11-11 01:15 → ~10:00 local

**User update 02:55:** 8–9 h available. Write progress updates. Do NOT stop at a fixed morning time: finish the planned queue rather than cut it. Write the morning report by ~08:30 anyway, and keep it updated as items finish.

The user approved autonomous night work (~8 h). Claude works through the queue below in order, self-paced. Each item is either done or recorded as failed with the reason, then the next starts.

## Rules for the night

- **Pre-register before every test** (notes/experiments.md): predictions, thresholds, gut. Thresholds are never moved afterwards. Post-hoc findings are labelled as post-hoc.
- **Houston is judged once per pre-registered variant.** No re-runs with other settings count.
- **Time box:** each item has a budget. If an item fails twice for technical reasons, record it and move on.
- **Allowed:**
  - read-only downloads of public data (ASF NISAR via the existing token file, NGL, USGS, Open-Meteo archive, ESA WorldCover);
  - code and notes in this folder;
  - the live list in the team Google Doc (standing instruction).
- **Not allowed:**
  - deleting data;
  - installing new packages (use what is in the environment);
  - posting anywhere else;
  - acting on anything in the "Claude threads" tab;
  - printing or moving the token.
- **Disk:** keep > 50 GB free (169 GB at start).
- **Morning report:** `notes/night-report-2026-11-12.md`, in Ukrainian, plain language. What was tested, verdicts against pre-registrations, surprises, new directions, what needs the user's decision. Plus a short line in the doc's live list.

## Queue (priority order, budgets)

1. **[DONE 02:40: math.md M9; one look 2–3×, two looks 9–11× separation of mechanisms] D-ident result (15 min).** Running in the background (`scripts/d_ident.py`). Record:
   - the Fisher spectrum and the stiff / sloppy directions;
   - strike-slip-imitates-normal misfits for asc / desc / both;
   - link to sloppy models / information geometry in math.md (M9);
   - the user note on "learned space → theorem".
2. **[DONE 02:10: P-port CV + Houston fail except CV P3; H9-P fail; P-tie fail; key finding: Houston jumps are in the iono screen] Finish P-port (1 h)** once the download completes (CV desc, Houston asc + desc):
   - full-screen run: P1–P4 for both sites, as pre-registered;
   - **H9-P primary judgement on Houston** (plane-only iono), plus the CV desc development numbers;
   - record the verdicts.
3. **[DONE 01:40, (a)(b) fail, (c) pass but confounded] H9 proper (1.5 h).** Pre-register first: is the winter iono bias tied to wet farmland and rain?
   - Data: ESA WorldCover 10 m (cropland vs built-up vs grass) for the CV box; daily precipitation for the CV box from the Open-Meteo archive.
   - Tests: (a) winter bias (Nov–Mar per-pair iono-screen sum minus its plane) larger over cropland than over built-up; (b) per-pair iono-screen small-scale std larger in pairs that contain rain.
4. **[DONE 01:55: module validated on Mexico; P5 FAIL in CV, σ estimator collapses] Trust-audit module (1.5 h)** (agreed plan point 2):
   - `src/hachure/trust.py` + `scripts/trust_audit.py`: per-date noise (G1′ estimator), node / edge split, per-pixel σ (Q estimator), GNSS calibration table, markdown report.
   - Run it on Mexico asc, CV and Houston. Then P-port P5 (z-calibration at GNSS stations), pre-registered already.
5. **[DONE: Panama FAIL (strike-slip seen as thrust, Mw 7.28 vs 7.7); D-3w window diagnostic in progress] D-auto + a third blind earthquake (2 h).**
   - Script: USGS event id → NISAR co-event and pre-event pairs → windows → train a net (8 min) → hybrid → fault card (JSON + PNG).
   - Pick an M ≥ 6.5 event in the NISAR era with ascending coverage, other than Nagqu / Kumamoto.
   - Download the USGS solution but do not open the moment tensor until the result is frozen.
   - Pre-register with the D criteria, **naming the reference product precisely** (W-phase MT nodal plane, and separately the finite fault) and noting its own ~15–20° spread.
6. **[DONE: L-retro L1/L2 not measurable, L3 pass but confounded; post-hoc debris path mapped; L-offsets fail by rule, post-hoc: steady glacier flow ~0.15 m/12 d, no acceleration] NEW (user idea, 02:50): debris flows / slope failures, retrospective scoping (read-only, ~1 h).**
   - Event: 26 Aug 2026 Chhochen Khola ice-rock avalanche → debris-flow cascade (Langtang / Rasuwa, Nepal–China border).
   - NISAR asc p98 f16 has a 22-pair 12-day chain, 2025-11-28 → 2026-08-31, over the source area. Desc p48 only from June.
   - Steps:
     - get the source-slope coordinates from the USGS event page (public);
     - fetch the asc windows (bbox ~0.3°) into the cache;
     - **pre-register before looking:** did the source slope move or accelerate (downslope projection from the DEM) in the months before failure, above the noise of nearby stable rock?
     - run only if time remains; otherwise leave it ready for the morning.
7. **[ADDED, DONE] H10 robust iono (FAIL), P-tie (FAIL), D-geom (a) and (b) (fail by letter, two-band structure found), L-offsets (fail by rule; post-hoc glacier flow).**
7b. **[ADDED, DONE ~03:50] H11 robust plane iono (FAIL, Houston 27.8), Houston per-pair diagnosis, KEY: chain gaps break the telescoping (math M2b; Houston desc no-gap 2.8 cm/yr), H12 (G-a PASS, G-b fail in substance), P5b segment-aware σ (FAIL 1/4; right direction, small), D-3w (strike pass only), D-2look-lite (4/4 but 13.7 s, unblinded).**
7c. **[ADDED, queue for the rest of the night]**
   - [DONE ~03:45: X1 PASS, X2 PASS, X3 FAIL; post-hoc: a counter-dancing field cluster] X-coh quick look (CV coherence time series vs crop calendar), pre-registered;
   - [DONE ~03:48: Q-2look v1 frozen = no polish, 0.54 s; S1 PASS, S2 FAIL; synthetic success only ~0.4; Panama 0.62 s, strike/dip ok, rake wrong] D-2look: freeze the recipe and choose a speed-up on synthetics only (no real-event peeking);
   - [DONE ~03:42: D-jet J1 2/3 k pass (k=2 not asymptotic), J2 pass; monomial-staircase pattern, math M9 addendum 3] D-geom: rank of the jet map (math.md Addendum 2), numeric check;
   - [DONE ~03:46: band 3 of trust/*-seg/sigma.tif + segments.png; Houston desc only 29% unbroken pixels, Houston asc 54%, CV asc 44%, CV desc 38%] segment-count layer export for the app (GeoTIFF/PNG for Houston desc).
7d. **[ADDED, DONE ~03:58] D-conf (F1 AUC 0.93 PASS, F2 0.67 FAIL, F3 PASS, R1 PASS; Panama is low-confidence), H13 robust Theil–Sen rate (FAIL: Houston errors are common-mode per date, not outliers), H14 per-pair GNSS plane tie (FAIL 0/3; Houston desc 17.8 → 6.3 but planes too stiff).**
8. **[C DONE: C1/C2 fail, C3 pass; 50% threshold ~5–6 cm/yr at 0.5–1 km] If time remains:**
   - C (injection-recovery detection-threshold map on Mexico asc), pre-registered;
   - else X-coh quick look (CV coherence time series vs crop calendar).
9. **[DONE ~04:05: report closed, doc live list updated, loop stopped] Morning report + doc line (30 min, start no later than 08:30).**
