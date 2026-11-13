# Mathematical structures behind the findings

What the SAR experiments look like as mathematics, and where else the same equations live.
Each entry says how sure we are: **proved here**, **known (cite)**, **conjecture**, or **analogy**.
Literature checks are marked TODO; nothing here is claimed as new until checked.

---

## M1. The pair network is a graph, and the noise is a cochain decomposition

**Setup.**
- Dates are vertices V, interferometric pairs are oriented edges E, and B is the |E|×|V| signed incidence matrix.
- A pair measures d = B φ + e:
  - φ_v is the per-date signal (atmosphere, true ionosphere, and deformation evaluated at that date): a 0-cochain.
  - e is per-pair noise (decorrelation, unwrapping, the per-pair ionosphere *estimate*): a 1-cochain.
  - In the language of discrete exterior calculus, d = δφ + e, with δ the coboundary.
- **Closure** of a triangle (u, v, w) is (δ d)(u, v, w) = (δ e)(u, v, w), because δδ = 0. Closure sees **only** edge noise. This is how H8 isolated the ionosphere-estimate error.
- Hodge decomposition of edge signals (on the clique complex): **gradient** (date-borne) ⊕ **curl** (closure-borne) ⊕ **harmonic** (cycles not filled by triangles).
  - A pure 12-day chain is a tree, so its gradient space is the whole edge space. First moments cannot separate node and edge noise there.

**Second moments do separate them.** Let node noise be independent with variances σ_v² and edge noise independent with variances τ_e². Then the edge covariance is

  C = B Σ_V Bᵀ + Σ_E,

so C is a vertex-weighted **lower Hodge 1-Laplacian** plus a diagonal.
- Off-diagonal entries: C_ef = ±σ_v² if edges e and f share vertex v (the sign depends on orientation), and 0 otherwise.
- Diagonal entries: C_ee = σ_u² + σ_w² + τ_e² for e = (u, w).

**Lemma M1.1** (proved here, elementary; TODO literature: LS-VCE, Teunissen & Amiri-Simkooei 2008; InSAR variance-component estimation).
- (a) σ_v² is identifiable from C iff deg(v) ≥ 2.
- (b) τ_e² is identifiable iff both endpoints have degree ≥ 2.
- (c) For a leaf v with edge e, only σ_v² + τ_e² is identifiable.
- (d) Nevertheless, the variance of any **chain sum** 1ᵀd along a path is identifiable for any graph, because it equals 1ᵀC1.

*Proof sketch.*
- (a) If deg(v) ≥ 2, pick two edges at v; their covariance is ±σ_v². If deg(v) = 1, σ_v² appears only in C_ee, together with τ_e².
- (b) and (c) follow by subtracting known vertex terms from the diagonal.
- (d) holds because 1ᵀC1 needs only sums of entries. ∎

*Used in:*
- G1′ (per-date noise from adjacent covariances, which A's diagonal-only NNLS could not identify);
- the morning "edge share" (≈ 43% of pair variance);
- V1.

**Contrast, the diagonal-only problem.** Using only diag(C) means solving |B| σ² = diag(C) with the **unsigned** incidence |B|. The signless Laplacian |B|ᵀ|B| is singular iff the graph is bipartite (known; Cvetković–Rowlinson–Simić). A path is bipartite, and the null vector is the alternating ±1. That is exactly the failure seen in experiment A, where NNLS zeroed every other date.

**Conjecture / direction M1.2.**
- When node noise is *temporally* correlated (seasonal water vapour), C gains non-adjacent entries.
- The structure of C then becomes "graph Laplacian + low-rank (seasonal) + diagonal". Separating these is a structured covariance-decomposition problem, in the spirit of robust PCA / factor analysis on graphs.
- For a product, this decomposition would give honest per-pixel error bars with no external data.

## M2. Edge noise accumulates as a random walk; node noise telescopes

- Chain sum S = Σ_k d_k = (φ_N − φ_0) + Σ_k e_k.
- Node noise contributes only its end dates. Edge noise contributes a random walk with variance N·τ².
- So a correction that turns node noise into edge noise (a per-pair estimate of a per-date quantity, like the split-spectrum ionosphere) can make **rates** worse while making every **pair** better (H8 / H8b).
- General principle: *never estimate a node quantity on edges if you will later integrate along a path.*
  - Known analogues: dead reckoning / odometry drift (per-step errors vs per-landmark errors); loop closure in SLAM is exactly the fix (closure = cycles in the pose graph).
  - So SLAM pose-graph optimisation is a ready-made toolbox for InSAR networks. **Analogy, worth trying.**

## M3. Residual maps are a turbulent phase screen (Kolmogorov / adaptive optics)

- Structure function std ∝ r^{1/3}, i.e. D(r) ∝ r^{2/3}: Kolmogorov in the regime r ≫ layer thickness.
  - For r ≪ thickness, the path-integrated screen gives D ∝ r^{5/3}.
  - The crossover scale ≈ effective turbulent layer depth (T1).
  - Saturation at 10–20 km ≈ von Kármán outer scale L0.
- Same mathematics as astronomical seeing (Fried parameter r0, frozen flow / Taylor hypothesis, tomography with several guide stars).
  - Directly reusable: phase-screen priors (von Kármán covariance) for kriging / GP correction, and multi-source tomography (GNSS + GOES bands as "guide stars").
  - **Known in AO; application to NISAR L-band urban noise: TODO check.**

## M4. Radar + geostationary imager = two-camera stereo

- E1's best shift is a disparity: the GOES parallax 0.71·h_cloud toward 307° plus the radar's own parallax h_wet·tan θ, toward the radar.
- The two passes look in opposite directions, so two disparities give two heights (cloud top and moisture layer). This is epipolar geometry with two very different sensors.
- **Conjecture S1**, not yet tested.

## M5. GOES BT variability is an observable volatility proxy

- σ²_date ≈ b·std(BT), with a ≈ 0 (B6, G1′). This is heteroscedastic regression with a known variance driver, like realised volatility / VIX.
- Natural home: a Kalman filter / state-space model with exogenous, time-varying measurement noise R_k. That is the form of a near-real-time NISAR time-series filter.

## M6. Hachures as Morse theory, plus topological persistence for "is this bowl real?"

**Analogy → direction.**
- Hachure strokes follow the gradient of the subsidence surface. The Morse–Smale complex of that surface partitions the city into basins of attraction, one per compaction centre.
- Persistent homology (0-dim, sublevel sets) ranks bowls by depth. With the noise model from M1 and M3, a persistence threshold of k·σ(scale) gives a **significance test for each bowl**.
- Product: show only bowls that persist above the noise.
- TODO: check TDA + InSAR literature.

## M7. Phase unwrapping is a topology problem (known)

- Residues are winding numbers, i.e. H1 obstructions. Branch cuts form a minimum-cost matching. Coherence masks behave like percolation clusters.
- Known (Goldstein 1988; Costantini network flow). Listed for completeness.
- Possible link: percolation threshold of the coherence mask as a predictor of unwrapping risk.

## M1.3. Result from V1: the chain-sum error is a difference of two large numbers

- var(rate) = [Σ var(r_k) + 2 Σ cov(r_k, r_{k+1})] / T². For per-date-dominated noise the two terms nearly cancel. On desc at 10 km: 81.7 − 80.7 cm².
- What survives is the edge part plus the end dates. This is the M2 telescoping seen numerically.
- Estimator problem (open):
  - Residualising with a rate fitted from the same pairs imposes Σ r ≈ 0, which adds −var/(n−1)-type covariance at every lag.
  - A cleaner design is cross-pass residualisation. An unbiased estimator of 1ᵀ C 1 under a self-fit constraint is a small, well-posed statistics problem. **TODO: derive.**

## M3b. Noise = turbulence × gradient, not amount (T2 / T3)

- Radar noise does not follow the moisture amount (ρ = 0.00 with balloon ZWD) or the moisture height (h50 the same morning and evening).
- The standard theory says the same thing. **Tatarskii (known):** C_n² = a²·L0^{4/3}·M², where M is the vertical gradient of the potential refractive index (in practice the humidity gradient) and L0 is the turbulent outer scale.
  - The delay variance is the path integral of C_n².
  - So noise needs **both** turbulence (L0, mixing) **and** a moisture gradient for the turbulence to act on. A full but laminar column is quiet.
- **Direction:** estimate ∫ L0^{4/3} M² dz per date from high-resolution soundings (M from humidity gradients; L0 from the mixed-layer depth). This should predict σ_date better than any single variable.
  - IGRA levels are too coarse; this needs BUFR high-resolution soundings.
- **Same equation elsewhere:** optical scintillation / astronomical seeing (C_n² profiles from SCIDAR), radio-telescope phase calibration (ALMA water-vapour radiometers), and acoustic propagation in the ocean (sound-speed structure).
  - ALMA's water-vapour-radiometer phase correction is the closest technical analogue of a GOES-based InSAR correction. TODO: read how ALMA scales WVR → phase.

## M1.4. Literature check of Lemma M1.1 (2026-11-10, subagent; abstracts plus one full text)

**Verdict.** Parts (a)–(c) are a special case of a **known** identifiability criterion from **network tomography**. The degree form itself was not found stated anywhere. Treat it as a clean corollary or folklore, **not a new theorem**.

**Network tomography (closest prior art).**
- Chen, Cao & Bu, *Network Tomography: Identifiability and Fourier Domain Estimation*, arXiv:0712.3618 / IEEE TSP 2010. Full text read.
  - Theorem 2: for Y = AX with independent X, the variances are identifiable iff the matrix of A's rows plus the elementwise products of row pairs has full column rank.
  - This condition is necessary for Gaussian X: Cao, Davis, Vander Wiel & Yu 2000, JASA 95:1063.
- Our model is d = [B I]·[φ; e]. A leaf's column of B equals ±(its edge's column of I), which gives (a)–(c) at once.
- Duffield & Lo Presti 2004, IEEE/ACM ToN 12(6):978: link variances from path covariances. This is the same idea as G1′, but in computer networks.

**Diagonal-only / bipartite contrast is classical.**
- Grubbs 1948, JASA 43:243: per-instrument variances from pairwise differences. Three instruments (an odd cycle) are solvable; two are not.
- Signless Laplacian singular iff the graph has a bipartite component: Cvetković, Rowlinson & Simić 2007, LAA 423:155.

**InSAR side.**
- C = A·Q_SLC·Aᵀ is standard: Hanssen 2001; Kampes 2006; Agram & Simons 2015, JGR 120:2752.
- Variance-component estimation (LS-VCE): Teunissen & Amiri-Simkooei 2008, J. Geod. 82:65; Cao et al. 2018, J. Geod. 92:185.
- Closure as a cycle-space residual: Zheng et al. 2022, IEEE TGRS 60; De Zan 2015; Ansari 2021.
- Per-date screens from pairs that share a date: Tymofyeyeva & Fialko 2015, JGR, common-scene stacking.
- **Not found:**
  - a Hodge-decomposition framing;
  - an explicit statement that one 12-day chain suffices (via adjacent covariances) for per-date variances;
  - "per-pair estimation of a per-date quantity turns node noise into random-walking edge noise".
- Fattahi, Simons & Agram 2017 and Piromthong, Hooper & Lazecky 2026 (IEEE TGRS 64:5203015) moved to *per-acquisition* ionosphere estimation because it is less noisy. That is the practical cure for M2, already in use.

**What is still plausibly ours:** two observations, not theorems.
1. Explicitly carrying the network-tomography criterion over to InSAR chains: the degree form, per-date atmospheric variance from a pure sequential chain, and the morning/evening "edge share".
2. The M2 framing that explains *why* per-acquisition ionosphere estimation is better for rates.

**Caveats:** VCE runs on residuals of an unknown deformation model (see V1's self-fit constraint); closure bias makes τ correlated across pairs.

**Lesson:** "find the field where this is a corollary" took about 4 minutes of search and saved us from claiming a theorem. Keep doing this for every candidate.

## M8. A Z₂ symmetry the network failed to break; linear parameters solved exactly (from D, 2026-11-11)
- **Setting.**
  - The simulator has an exact symmetry g: (slip vector s, image y) → (−s, −y), because y = G(θ)·s + noise, linear in s, and the noise distribution is symmetric.
  - The joint distribution p(θ, s, y) is therefore Z₂-invariant. The true posterior is Z₂-equivariant: p(θ, −s | −y) = p(θ, s | y).
- **What happened.** The trained network became Z₂-*invariant* instead of equivariant: q(· | y) = q(· | −y).
  - The only way to be invariant and calibrated is to make the posterior over s symmetric, so the sign information is thrown away.
  - This is a symmetric critical point: the "average over the group orbit" solution.
  - Related known phenomena:
    - symmetric saddles in phase retrieval (x vs −x);
    - sign ambiguity in ICA / PCA;
    - "simplicity bias": magnitude features are learned first, sign features need a product (sign × location relative to the fault), which is a harder, higher-order feature.
- **Two cures.**
  1. **Equivariance by construction (representation theory of Z₂).** Split the head: geometry from even features, s from odd features.
     - e.g. h_even = f(y) + f(−y), h_odd = f(y) − f(−y);
     - s-mean = W·h_odd ⊙ (gate from h_even).
     - This is the smallest case of group-equivariant networks (Cohen & Welling 2016): the irreps of Z₂ are "even" and "odd".
  2. **Variable projection (Golub & Pereyra 1973; separable nonlinear least squares).**
     - Never learn the linear parameters. Learn θ, and for each θ solve s exactly: ŝ(θ) = G⁺(θ) y.
     - The conditional p(s | θ, y) is Gaussian for Gaussian noise. With our non-Gaussian real noise, it is the push-forward of noise samples through G⁺.
     - In Monte Carlo language this is **Rao–Blackwellisation**: integrate out analytically what you can. Lower variance, and equivariance is automatic.
     - **Implemented (v3).** Sense mass 0.50 → 0.71 (0.85 for larger events).
- **General lesson, maybe worth a short note.** In SBI for linear-in-amplitude forward models (any Green's-function problem: geodesy, seismology, EEG/MEG sources, gravity), *do not* put amplitudes in the network.
  - The network learns geometry; amplitudes come from linear algebra.
  - Separable structure + learned proposal = exact conditional. Cheap, and immune to the sign trap.
  - Literature check needed: "simulation-based inference variable projection", "marginalised linear parameters neural posterior estimation". Probably done in gravitational waves, where distance and phase are marginalised analytically (Dax et al. DINGO).

## M9. A fault seen by radar is a (moderately) sloppy model (from D-ident, 2026-11-11)
**Setting.** Okada rectangle → deplaned LOS map, Nagqu-like normal fault (192/50/−119, 30 × 15 km).
- 9 parameters: e0, n0, depth (km), strike, dip, rake (rad), log L, log W, log slip.
- Fisher information F = JᵀJ / σ² with 1 cm white noise. Script: `scripts/d_ident.py`.

**Spectrum (1σ along each eigen-direction).**

| Looks | Stiffest σ | Sloppiest σ | Condition number |
|---|---|---|---|
| asc | 0.0004 | 0.021 | 2.9e3 |
| desc | 0.0004 | 0.023 | 2.9e3 |
| asc + desc | 0.0003 | 0.011 | 1.4e3 |

- The 9 eigenvalues spread over 3.5 decades, roughly evenly in log: the signature of a sloppy model (Sethna, Transtrum et al.). This one is mild; systems biology models span 6–12 decades.
- **Stiffest direction:** always dip (−0.9 to −0.96), with some strike. The radar measures fault dip best.
- **Sloppiest direction:** where the fault sits along its own length and how deep it is.
  - asc: depth + north position (+0.75 depth, +0.61 n0);
  - desc: east position;
  - asc + desc: north position (the strike is ~N–S, so this is "slide the fault along its own strike").
  - Each look has its own blind direction, and two looks cut the sloppiest σ in half: stereo for parameters.

**Mechanism identifiability (part B, noise-free).** We forced the rake into strike-slip families and fitted the normal-fault map.

| Looks | Left-lateral residual (share of signal) | Right-lateral residual (share of signal) | Normal-rake control residual (share of signal) |
|---|---|---|---|
| asc | 5.07 cm (38%) | 4.47 cm (34%) | 2.23 cm (17%) |
| desc | 4.56 cm (31%) | 4.10 cm (28%) | 1.37 cm (9%) |
| asc + desc | 7.37 cm (53%) | 5.68 cm (41%) | 0.64 cm (5%), recovers 193/48/−119 |

- The control does not reach zero because the optimiser is imperfect. The ratios still say a strike-slip imitation is 2–3× worse with one look, and **9–11× worse with two looks** (the control then recovers the truth to 1°).
- With real Nagqu noise of ~1–2 cm, one look **can** tell the mechanisms apart. So SBI v3 calling Nagqu strike-slip (rake +9) was a *learning* failure, not an identifiability limit, and the hybrid fixed it (normal fault).

**Cross-field links.**
1. **Sloppy models and information geometry** (Transtrum, Machta, Sethna 2010–2015). The model manifold has a hierarchy of widths, and its boundaries are simpler models. They argue sloppiness arises generically from multipole-like expansions.
   - **Candidate statement for us:** the boundary faces of the Okada manifold are the multipole truncations of the source:
     - depth → ∞ or L, W → 0: a point moment tensor;
     - W → 0: a line source.
   - The sloppy eigenvalues should then decay like (fault size / centroid depth)^(2k), with k the multipole order.
   - Testable in an hour: Fisher spectrum vs depth/size, look for a power law. This is very likely an *instance* of Transtrum's general argument, not a new theorem, but a clean, checkable instance in geodesy. The MBAM (manifold boundary approximation method) would then *derive* the point-source and line-source approximations geodesists use by hand.
2. **The user's "learned space → theorem" idea.**
   - The NPE network's learned summary s(y) induces a metric on parameter space: the pull-back of the Fisher metric.
   - Known neighbour: IMNN (Charnock, Lavaux, Wandelt 2018; Alsing & Wandelt) trains summaries to maximise Fisher information.
   - **Concrete probe:** compare the net's posterior covariance eigen-directions with the Fisher sloppy directions above.
     - If they align, the net has "learned the information geometry" and its latent space can be read as a coordinate chart where sloppy combinations are explicit.
     - Then a reduced model (fewer effective parameters) can be read off the learned space, exactly the "learned space → simpler theory" move.
   - Not tested yet. Recorded as research direction D-geom.

**M9 addendum: D-geom (a) result (post-hoc reading).**
- The Fisher spectrum vs centroid depth splits into two bands: 4 directions ∝ d^−5.2 and 5 directions ∝ d^−8.7. The gap is ≈ 3.5, close to the multipole prediction of 4.
- The second band contains horizontal position, not only shape. Mechanism: the data are deplaned, i.e. we observe the image of the field under the projector P onto the orthogonal complement of the ramp space R = span{1, x, y}.
  - A translation changes the field by −a·∇u. When the source depth d ≫ window size ℓ, ∇u over the window is a constant plus O(ℓ/d) corrections, so P(a·∇u) = O((ℓ/d) · |∇u|).
  - Translation is pushed down one order: a near-gauge direction.
- **Candidate statement (checkable):** for a source of depth d observed through a window of size ℓ with polynomial-ramp removal of degree k, the Fisher eigenvalues group into bands by (multipole order) + (order of the polynomial suppression), with band exponents shifted by 2 per order of (ℓ/d).
  - That would make "which fault parameters a deplaned radar map can see" a counting problem: multipole order × the polynomial gauge.
  - Test: rerun with degree-0 removal (constant only) and degree-2 removal. Position should move back to the stiff band with degree 0, and more directions should drop with degree 2.
  - This is the kind of "learned space → statement" the user asked about; here it came from the Fisher geometry, not from a learned latent space.

**M9 addendum 2: D-geom (b) partly refutes the gauge story.**
- Varying the removed polynomial degree k = 0, 1, 2 steepens *all* exponents by ~1.2 per degree.
- Position stays sloppy even with only a constant removed. Band count: 5, 5, 6.
- Revised picture: the observation is a window of size ℓ ≪ d onto a far field, so the data are effectively the Taylor jet of the field at the window centre, up to some order.
  - Removing degree-k polynomials deletes the first k+1 orders of the jet. Every parameter derivative is a jet too, so each eigenvalue scales like (ℓ/d)^(2m) with m set by the lowest surviving order its direction reaches.
  - The 4 + 5 split is then "directions that reach the first surviving order" vs "directions that need one order more". The second group includes position, because translation acts on the jet as differentiation, which mixes orders.
- **Candidate statement (still unproven):** for a point-like source seen through a small window with degree-k removal, the Fisher eigenvalue exponents in d are 2·(order) + const. The band multiplicities equal the ranks of the parameter-to-jet maps at successive orders.
  - That is a finite linear-algebra computation (rank of the jet map), so the band counts can be predicted exactly and checked against the numbers above.
  - This is the kind of "learned space → statement" the user likes, but here it comes from Fisher geometry. Literature check needed: far-field identifiability of point sources, "jet" or moment methods in potential-field geophysics (Euler deconvolution uses exactly such local derivatives).

**M9 addendum 3: D-jet, the bands are monomial counting (2026-11-12 night).**
- **Setting.**
  - θ ∈ ℝ⁹ are the fault parameters.
  - u(x; θ) is the LOS field on a window B_ℓ (side 2ℓ) around a point at distance ~d from the source, with ℓ ≪ d.
  - P_k removes polynomials of degree ≤ k. The Fisher matrix is F(ℓ) = J(ℓ)ᵀJ(ℓ), where J's columns are P_k ∂_θj u on B_ℓ.
- **Mechanism.** Write each ∂_θj u as its Taylor jet at the window centre, c_j = Σ_m c_j^(m) · (x/d)^m, where c^(m) has m + 1 entries. On B_ℓ the order-m part has size (ℓ/d)^m.
  - So J(ℓ) = M · diag((ℓ/d)^m) · C up to an O(1) well-conditioned basis change. M is the monomial basis on the unit window, C stacks the jet blocks C_m (rows: coefficients of order m; columns: parameters).
  - For a graded product like this, the singular values split into groups s ~ (ℓ/d)^m. The group sizes are the incremental ranks r_m = rank[C_{k+1}; …; C_m] − rank[C_{k+1}; …; C_{m−1}].
  - This is the standard flat-limit / graded-matrix argument.
- **Numerical result (D-jet).**
  - For a Nagqu-like Okada source, A (Fisher slopes in ℓ) and B (algebraic jet ranks) agree exactly for k = 0 and 1. They agree for k = 2 once ℓ/d ≤ 0.05 (post-hoc).
  - The incremental ranks are always maximal: r_m = min(m + 1, n − Σ_{j<m} r_j).
- **Statement (conjecture, checked numerically for one source and three k).**
  - For a generic n-parameter source observed through a small window with degree-k polynomial removal, the Fisher eigenvalues fall into bands λ ~ (ℓ/d)^(2m), m = k+1, k+2, ….
  - The band sizes are those of the 2-D monomial staircase: k + 2, k + 3, …, truncated so they sum to n.
  - Equivalently: **how many source parameters a far-field radar window can resolve at a given SNR is a counting problem on monomials.** n_resolved(SNR) = Σ_{m ≤ m*} (m + 1) − (k+1)(k+2)/2, capped at n, where m* is the largest order whose signal (ℓ/d)^m · |u| exceeds the noise.
- **Cross-field link (literature check needed; names from memory).**
  - Exactly this grouping, eigenvalues of size ε^(2m) with multiplicities equal to the number of degree-m monomials, is known for **kernel matrices in the flat limit** of radial-basis-function interpolation:
    - Fornberg, Driscoll, and Schaback studied the flat limit;
    - Wathen & Zhu gave the spectral distribution;
    - Barthelmé & Usevich (SIAM J. Matrix Anal. Appl., ~2021) treat the general smooth-kernel case with graded ranks of Wronskian-like matrices.
  - Our F(ℓ) is a Gram matrix of smooth functions on a shrinking domain, which is the same object. So the fault-identifiability band structure would be a corollary of the flat-limit theory, read in a new place: the "points" are fault parameters, and the "flatness" is the window shrinking relative to the source distance.
  - It also gives a clean bridge to **Euler deconvolution** (potential-field geophysics), which estimates source position from the first few local derivatives: the staircase says exactly how many source parameters each extra derivative order buys.
- **What is still missing for a theorem.**
  - (1) A proof that the Okada jet map is generic (full graded rank) away from a measure-zero set of geometries; symmetric geometries could drop rank.
  - (2) Uniform asymptotics (the k = 2 bands mix slowly).
  - (3) The depth-scaling version (D-geom (a) bands were in d at a finite window, a different limit).
- **For the user:** this is the strongest "learned / measured space → statement" candidate so far.
  - The statement is short, it is checkable by algebra (ranks), and it connects to a known theorem elsewhere (flat-limit kernel spectra).
  - Worth a proper literature check before calling it new.

## M2b. Gaps break the telescoping (from the Houston diagnosis, 2026-11-11)
- **Setup.** Path graph of dates t_0 … t_N. Pair k observes d_k = Δdef_k + n_{k+1} − n_k + e_k (node noise n, edge noise e).
  - The chain-sum rate over a complete path is (def_N − def_0 + n_N − n_0 + Σ e_k) / T. Node noise enters only through the two end dates.
- **With missing edges** (pixel invalid in some pairs) the path splits into s connected segments [a_j, b_j]. The pooled estimator Σ_valid d_k / Σ_valid Δt_k becomes:
  - rate + Σ_j (n_{b_j} − n_{a_j}) / T_valid + Σ e / T_valid.
  - **Var_node = 2s · σ²_node / T²** (for independent, equal node variances).
  - Each gap adds two unpaired node terms. In graph language, the rank of the boundary operator restricted to the valid edges drops, and the node noise's projection onto the cokernel grows with the number of components.
- **Consequence.** In a humid city, where a local node (one date) can be off by ~1 m (TXET), one gap is enough to ruin the rate. The trust layer must use s per pixel.
- **Better estimators:**
  - (a) complete-chain pixels only;
  - (b) weight segments by span² (least-squares slope per segment), which down-weights short segments;
  - (c) bridge a gap from spatial neighbours *only if* the node errors are spatially smooth. At TXET they are not.
- This is a direct corollary of Lemma M1.1 (cochain decomposition): the node part lives in the image of the coboundary of the *valid subgraph*, and telescoping needs that subgraph connected.
