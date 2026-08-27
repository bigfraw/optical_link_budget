# Backlog — the unimplemented and unwired work

Date: 2026-08-26. Two repository sweeps made this list: one sweep of the
documentation (the READMEs, docs/, CLAUDE.md) and one sweep of the code (the
stubs, the warning flags, the debt comments). The rule for the next sessions:
**fix the items on this list before you add a feature.** Each new feature
extends this list faster than the fixes close it.

Groups: fidelity 0 (analytic), fidelity 1 (statistical / Monte Carlo),
fidelity 2 (wave optics), infrastructure and code debt, documentation debt,
external dependencies. Inside each group: wiring steps (built, but no budget
consumes it), physics gaps, and numerical or validation issues. Line numbers
are from 2026-08-26 and can drift.

---

## Top of the stack (the recommended order)

1. **Gap 2 — the pre-compensated uplink models NO scintillation.** The one
   item both sweeps call MAJOR. See 0-W1.
2. **HIGH (owner-flagged, 2026-08-27) — stop the reliance on the
   `DEFAULT_HS` 20-layer array.** HV5/7 is a continuous profile; the planner
   and the physics must take a callable, not a hand-discretised grid. See
   2-I2.
3. **DONE — the turbulent screen-count floor `min_screens`.** Work package 7
   resolved it. See 2-N1.
4. **Gap 3 — thread the beam curvature f0 into the Fried call site.** Small,
   already prepared at the physics layer. See 0-W2.
5. **The stale docs that contradict the code.** Cheap, and they mislead every
   later session. See the documentation-debt group.
6. **The kernel repo commit.** One `git add` in `my_analysis_modules`, plus
   the KR-24 constants. See the external group.
7. **The owner decisions.** `downlink_budget` default, the fidelity-2 Term,
   the FAST-versus-field reference model. Decide once; then wire. See 0-W5,
   2-W1.

---

## Fidelity 0 — analytic

### Wiring steps (built, no budget consumes it)

- **0-W1. Gap 2: add the residual-scintillation Term to the pre-compensated
  uplink (MAJOR).** The beacon + AO budget is phase-only; its Terms flag
  `NO SCINTILLATION` (olb/links/uplink.py:285, :398). The missing number is
  built: `andrews.paths.uplink_scintillation_index(tracked=True)` (Ch. 12,
  Eqs. (57)–(60); see olb/turbulence/andrews/paths.py:625). No budget reads
  it. Do not trust the corrected uplink fade until this Term exists. Docs:
  CLAUDE.md:207, docs/api-budget.md:241, docs/physics.md:760.
- **0-W2. Gap 3: thread the curvature f0 into the Fried parameter.**
  `andrews.beam.beam_params` takes any f0; `gaussian_fried_parameter` stays
  collimated, and the call site in olb/models/coupling/terrestrial.py passes
  no f0. A deliberately diverged beam gets a wrong r0. Docs:
  docs/physics.md:594, crosscheck GF-01.
- **0-W3. Gap 1: the aperture angle-of-arrival tilt feeds no Term.**
  `andrews.structure.angle_of_arrival_variance` is built and delegated to;
  no coupling Term adds contribution C, so the received tip-tilt is a lower
  bound. Note conflict C-04: olb holds two tilt conventions (gradient 0.174
  vs the Noll route in ao.py); a caller that adds them must say which.
- **0-W4. Gap 6 and Gap 7: `l0`/`L0` and the temporal faces have no
  consumer.** No Term passes an inner or outer scale. No Term reads the
  Greenwood frequency, tau0, the fade rate, or the fade duration
  (andrews/temporal.py). The roadmap wants a tracking-bandwidth / servo-lag
  Term (README node NT7; also the TODO at
  olb/models/coupling/terrestrial.py:108).
- **0-W5. `downlink_budget` still defaults to `model="lognormal"`.** The
  `model="auto"` selector exists and is opt-in. The switch moves the
  strong-regime total by several dB. OWNER DECISION.
- **0-W6. The Andrews aperture-averaging chain is unused.** The downlink
  Term keeps the numerical Airy-filter integral; the terrestrial Term keeps
  the Churnside fit (optimistic 5–13 % against the book's exact Eq. (60)).
  Conflict C-06: Ch. 12, Eq. (39) is the closed form that should supersede
  the integral when a work package moves the Term over.
- **0-W7. The K distribution and the lognormal-Rician PDF are unused**
  (andrews/distributions.py). The three new spectrum labels in
  olb/assumptions.py (Tatarskii, exponential, modified) have no Term.
- **0-W8. The Andrews wander route stays on the shelf BY DECISION** (C-01
  closed via Belmonte; the budgets keep the Dios/Belmonte 2.07 kernel; the
  Andrews 7.25 route sits in andrews/wander.py for measurement only). Listed
  so nobody re-opens it.

### Physics gaps

- **0-P1. Laser guide star: focal (cone) anisoplanatism.** `LaserGuideStar`
  is a placeholder; `uplink_budget` raises (olb/links/uplink.py:495,
  olb/scenario.py:77). Its cone anisoplanatism differs from the point-ahead
  angular form.
- **0-P2. The short terrestrial retro link has no module.**
  `retro_space_budget` assumes a long slant range and independent legs. The
  book also gives a backscatter amplification and a coupled double passage
  that olb does not model (crosscheck RT-01..03).
- **0-P3. The strong / moderate aperture-averaged downlink index.** The
  gamma-gamma Term is a POINT receiver (olb/links/downlink.py:209); its fade
  is deeper than the true aperture fade. The book gives no slant
  aperture-averaged strong index; a second source or a derivation is needed.
  (Owner earlier decided to keep point-receiver; see the memory
  `downlink-strong-aperture-decision`.)
- **0-P4. Gap 8: the annular (obscured) receive aperture.** No book source
  exists; one central-obscuration gap touches many Terms
  (olb/links/downlink.py:109, uplink.py:122, terrestrial.py:183, the
  coupling modules). Needs a new reference.
- **0-P5. The terrestrial scintillation Term carries no pointing jitter.**
  It is on-axis (r=0) only; the jitter-into-beta fold sits in the Monte
  Carlo path only (TODO at olb/links/terrestrial.py:132). Fold the jitter
  into the off-axis radius, and add the mean-power loss Term.
- **0-P6. No tracking-loop bandwidth model.** The tip-tilt correction is
  all-or-nothing (olb/models/coupling/terrestrial.py:108). A finite servo
  leaves a residual tilt. Pairs with 0-W4.
- **0-P7. The AO cutoff `f_c = sqrt(n_modes)/(2D)` is a heuristic with no
  source** (crosscheck AO-06, unmatched).
- **0-P8. The pointing Term uses a small-aperture on-axis approximation**
  with no aperture averaging of the fade (olb/models/pointing.py:111). The
  exponential-in-dB fade law is an olb extension with no book counterpart
  (crosscheck PT-02..04); the book route is Ch. 12, Eq. (53) (row G-140).
  Conflict C-09: the `beta2 += 2 (sigma_theta L)^2` jitter fold is correct
  arithmetic but uncited — label it and compare.
- **0-P9. The extinction model is one-parameter** (one zenith tau or one
  dB/km). Wavelength-resolved / MODTRAN stays planned (README NA1). The
  airmass model breaks below 5 deg elevation (olb/models/extinction.py:103).
- **0-P10. Only one Gaussian transmit beam; `Transmitter.m2` is dead.** No
  model reads m2. Flat-top beams and incoherent aperture diversity stay
  planned (README TB3, TB4).
- **0-P11. `eta_max(a)` assumes a uniform, flat-wavefront aperture.** A
  near-field terrestrial link breaks both; the error runs safe. The docs
  call this the open Gap-3 upgrade (docs/physics.md:969).
- **0-P12. The MMF NA gate is a flat factor.** No turbulence re-broadening
  of the focal spot, no mode-count saturation; `optimal_focus` is geometric
  (docs/physics.md:1134).
- **0-P13. The uncorrected coupling curve extrapolates past its D/r0
  limit** and only warns (olb/models/coupling/downlink.py:145,
  terrestrial.py:323).
- **0-P14. The strong-turbulence effective beam parameters are deferred.**
  `effective_beam_params` (Theta_e, Lambda_e) is coded and has no caller
  (gaussian_fried.py:330; crosscheck GF-17). Deliberate, to match the Dios
  weak regime.

### Numerical and validation issues

- **0-N1. TL-05: the terrestrial weak gate tests ONE criterion.** Ch. 5,
  Eq. (16), printed p. 140 needs `sigma_R^2 < 1` AND
  `sigma_R^2 Lambda^{5/6} < 1`. A focused or strongly diffracted beam can
  pass a gate it must fail. Fix with one shared helper (Table 2 row G-20).
- **0-N2. The strong-fluctuation parameter q carries 1.22 one level too
  high** at olb/turbulence/plane_wave_scintillation.py:284 (crosscheck
  ledger line 742). Verify and fix.
- **0-N3. The budget adds dB losses / per-term p-quantiles.** Not a book
  form; it is a conservative upper bound (crosscheck RS-02, RS-03). Record
  the bound; consider a joint-sample check.
- **0-N4. The near-field truncation flag is not conservative.** The true
  value can sit above or below the far-field form
  (olb/models/gaussian_efficiency.py:189). The fidelity-2 vacuum layer is
  the verifier; consider an auto-route.
- **0-N5. The quasi-frequency has no upper limit of its own** (`b_2` grows
  as `f_max^{1/3}`). Any caller must set the band from the detector or an
  inner scale. Pairs with the inner-scale memory.
- **0-N6. The fade rate and the fade time have no external check** — the
  book gives no worked example; the faces are checked against internal
  identities only.

### Blocked by the source (documented refusals — need a second source)

Each of these raises `NotImplementedError` with a citation; the book prints
no form, and olb guesses no coefficient. Full list: docs/physics.md:852–897.
The path forward for each is a second reference or a derivation.

- The Gaussian two-scale STRONG branch (the unresolved eta_X of Ch. 9,
  Eq. (109)) and its aperture chain (andrews/scintillation.py:424,
  andrews/aperture.py:351).
- The weak Gaussian-beam aperture flux variance (Ch. 10, Eq. (78) is a
  numerical double integral; andrews/aperture.py:368).
- The convergent-beam Rytov variance (Theta0 < 1;
  andrews/scintillation.py:340).
- Inner/outer scale on any Ch. 12 slant form (andrews/paths.py:401); the
  full Gaussian-beam downlink index Ch. 12, Eqs. (36)–(37) (row G-133).
- The temporal spectrum with a finite inner or outer scale; the strong
  spherical/Gaussian temporal spectra; the weak-only aperture-averaged
  temporal spectrum (andrews/temporal.py:424, :436, :446).
- The Gaussian row of the modified spectrum (App. III, Table III;
  andrews/structure.py:271).
- The K distribution below sigma_I^2 = 1.

---

## Fidelity 1 — statistical / Monte Carlo

- **1-1. The terrestrial statistical coupling fade does not exist.** The
  SMF Term is mean-only, so it LOCKS the budget out of a fade margin
  (olb/models/coupling/terrestrial.py:302, olb/results.py:299). FAST is
  far-field only; a near-field Gaussian beam needs the split-step model —
  so the real fix is the fidelity-2 wiring (2-W1).
- **1-2. FAST limits NT1–NT4.** Point-ahead is off (`DTHETA=0`;
  olb/models/coupling/fast.py:239) — compute it from the orbit and pass it.
  Scalar elevation only (fast.py:132). No obscuration in the coupled flux
  or the mean-only fibre model. No tip-tilt wander removal.
- **1-3. Strong-fluctuation routing.** The uplink and terrestrial links
  only WARN when the weak-fluctuation limit is exceeded; the downlink
  already routes to gamma-gamma. Route the other links to a strong-regime
  model (a Monte Carlo or the fidelity-2 layer; see the memory
  `strong-fluctuation-numerical`).
- **1-4. The Dios duplicate (ponytail DEBT).** The analytic
  beam_wave_scintillation path and the coupled-flux MC duplicate the same
  equations; the jitter correction sits in the MC path only
  (olb/turbulence/beam_wave_scintillation.py:46). Converge on ONE
  implementation; do not make a third copy.

---

## Fidelity 2 — wave optics

- **2-W1. No fidelity-2 Term reaches a budget (OWNER DECISION).** The layer
  is built and self-checked (vacuum: 0.011 dB agreement; turbulent uplink
  reciprocity: 0.18–0.54 dB against the coupled-flux MC). The blocker: the
  FAST model and the field disagree on the fibre coupling (the field reads
  0.7–2.9 dB LESS loss than FAST, and 2.5 dB less than the terrestrial
  analytic Term). Pick the reference model; then wire the Term. The
  downlink `"montecarlo"` slot in olb/links/downlink.py:247 is reserved for
  exactly this.
- **2-N1. `min_screens` and `_merge_layers` — DONE (work package 7).**
  `_merge_layers` now clamps a weak path UP to EXACTLY `min_screens`
  contiguous Cn2-weighted groups, through the new `_equal_weight_groups`.
  The old bail-out, which returned one screen per Cn2 layer, is gone, so a
  200-layer profile no longer gives 200 screens. A profile that has fewer
  layers than `min_screens` warns and keeps its layers, because the planner
  does not split a layer. The Rytov cap still raises the count above the
  floor on a strong path. The integers 15 / 9 / 5 are CONFIRMED by an olb
  convergence sweep, and the docstring says that the source is olb and not
  Schmidt, which gives no floor. The absolute lower bound is 4, the moment
  count of Eq. (9.65), printed p. 164. The three turbulent examples were
  re-run at the new screen counts. The sweep table and the re-run numbers
  are in the WP7 note of docs/schmidt-crosscheck.md.
  DEFERRED, and still open: the planner does not SOLVE Eq. (9.65), and
  `SamplingReport` does not carry the moment error. The module self-check
  measures it instead, against `olb/waveoptics/schmidt/turbulence.py`, and
  the Cn2-weighted centroid grouping holds every moment inside 1 percent.
  olb/waveoptics/turbulence/sampling.py:271, :311.
- **2-I2. Continuous Cn2 profiles — drop the `DEFAULT_HS` crutch (HIGH,
  owner-flagged 2026-08-27).** HV5/7 and the other Cn2 models are continuous
  functions; the 20-layer `DEFAULT_HS` array is a hand-made discretisation,
  and it leaks into the physics wherever a decision reads the grid instead of
  the profile. Work package 7 removed the worst leak (the screen count), but
  the screen PLACEMENT still comes from the array. The post-WP7 matched-seed
  measurement (the WP7 note in docs/schmidt-crosscheck.md) sharpened the
  question: the bottom screen HEIGHT is a null, and the live variable is
  whether the near-ground `Cn2` is spread over many thin screens or lumped
  into one — a resolution question that only the continuous ground layer can
  answer. The same measurement shows the placement moves the deep SMF fade
  tail (about 2 dB at p5, direction consistent, not yet resolved above the
  Monte-Carlo noise at 200 trials). The change, in two separate steps: (1) `turbulent_grid` and
  `_plan_space` in olb/waveoptics/turbulence/sampling.py accept a callable
  `cn2(h)` and compute the group integrals, centroids, Rytov shares, and the
  Eq. (9.65) moments by quadrature on the callable; `DEFAULT_HS` stays only
  as the fallback for an array caller. This also makes the Gauss-quadrature
  screen placement (tracker candidate, S-22) implementable. (2) LATER, and
  separately: the fidelity-0/1 modules that integrate over `hs` arrays
  (slant extinction and scintillation, uplink flux, FAST) move to callables;
  that step is wide, mechanical, and must move no numbers. The owner decided
  on 2026-08-27 to flag this here and NOT build it yet.
- **2-S1. The Schmidt cross-check gaps S-01 to S-28.** The Schmidt
  foundation layer (`olb/waveoptics/schmidt/`) is validation only, and its
  tracker holds 28 numbered gaps between the book and the production
  wave-optics code. Do NOT copy the rows here; read Table 2 of
  docs/schmidt-crosscheck.md. The HIGH-priority rows are: S-13 (the
  co-moving route has no book equation; compare it against the scaled flat
  grid of Eq. (6.65) on the 600 km uplink — the WP5 example now measures
  1.7e-3 soft and 2.3e-2 hard at m = 247), S-14 (the split step holds ONE
  flat pitch; Eq. (8.18) gives each step its own, so the grid cannot grow
  with a diverging beam), S-16 (constraints 1 and 2 of Ch. 7 are
  implemented nowhere; `guard = 4.0` and `pixels_per_feature = 16` have no
  source), S-21 (the three turbulent geometry constraints of Eqs. (9.86) to
  (9.88) are never evaluated), and S-22 (the layer moment rule; the same
  item as 2-N1). The MEDIUM rows S-15 (the absorber shape), S-17 (the
  Fresnel minimum distance), S-20 (the phase pitch rule) and S-27 (the
  aotools subharmonic screen) are recorded there too. Each one is an owner
  decision, because each one moves a production number.
- **2-P1. The temporal (frozen-flow) axis is a stub.** `TemporalScreens`
  raises (olb/waveoptics/turbulence/temporal.py:54); the layer gives
  snapshots only — no fade rate, no fade duration. The design note lives in
  the class docstring.
- **2-P2. The folded / retro double pass is a stub.** `folded_terrestrial`
  and the `"retro"` direction raise (run.py:111, :252, :393). The two
  passes share screens, so they are correlated; that needs its own design.
- **2-P3. No co-moving (spherical) screen.** `split_step` takes a flat grid
  only; a long slant path pays the pixel cost.
- **2-P4. The reciprocity route carries no point-ahead anisoplanatism**
  (the uplink and downlink read the same screens;
  docs/api-waveoptics.md:644).
- **2-I1. `TurbWaveResult` is a minimal scalar record — DO NOT extend it
  piece by piece.** The rich record (the E-field inside the receive
  aperture) gets its own design session (run.py:83; memory
  `waveoptics-results-deferred`).
- **2-N2. Known numerical readings to keep in view:** the Fourier screen
  structure function reads up to 15 % low over r/r0 0.3–1.6 (ratios only);
  the aperture-averaged analytic factor fails when the aperture holds the
  beam (the 100 mm bucket case); the grid sizer warns past `forvard_max_z`
  and under the `n_max` clamp (olb/waveoptics/grid.py:181, :200); the
  runner warns when the receive aperture reaches the absorbing band
  (turbulence/run.py:282) — an automatic grid-from-aperture size would
  close that one.
- **2-N3. Speed: tune the grid size and resolution along the path**,
  validated against the well-sampled reference runs (memory
  `waveoptics-speed-exploration`).

---

## Infrastructure and code debt

- **I-1. Scalar-elevation limits.** The gamma-gamma Term
  (olb/links/downlink.py:185) and the FAST Term
  (olb/models/coupling/fast.py:132) refuse an elevation array. Vectorise or
  loop internally.
- **I-2. Duplicate physics copies (Gap 10).** The Rytov standard deviation,
  the plane-wave coherence radius, and the Fried parameter each exist in
  three places; the lognormal dB faces exist twice (crosscheck TL-01..04,
  GF-05, GF-10, KR-23, KR-25). Converge.
- **I-3. The diverged (Theta, Lambda) coupled-flux feed has no cross-check
  test** against the closed-form on-axis index (TODO at
  olb/turbulence/uplink_flux.py:94). The measurement exists in the
  crosscheck (+3.06 % / −1.57 % / −0.02 %); turn it into a test.
- **I-4. Dead or decorative parameters.** `geometry` is unused in every
  terrestrial Term (kept for the signature); `precompensation` is silently
  ignored on a downlink or retro link (olb/scenario.py:165) — refuse it or
  document it as uplink-only.

---

## Documentation debt

- **DD-1. Three docs still call the aperture angle-of-arrival function a
  raising stub** — the code now delegates (docs/architecture.md:44,
  docs/physics.md:1152, docs/andrews-crosscheck.md:53). Fix the three.
- **DD-2. README node NT5 still shows the diverged-feed check as planned**
  — it is measured and closed (README.md:211).
- **DD-3. DONE.** docs/api-waveoptics.md now carries the `min_screens`
  caveat, the `rmax` factor-4 note and the `fresnel_weight_min` note in the
  `QualityPreset` table.
- **DD-4. Crosscheck Table 3 is partly stale** — the "not found in olb"
  spectra rows predate WP3's andrews/spectra.py.
- **DD-5. Citation faults:** ao.py:151 credits "Andrews Ch. 3" for a Noll
  1976 result (AO-07); the aperture-averaged integral cites Ch. 12 where
  Ch. 9 Eq. (25) / Ch. 10 Eq. (59) are closer (PW-05); the C-02 reference-
  plane docstring note is still owed.
- **DD-6. The point-ahead decorrelation framing is fixed in code, not in
  the docs** (memory `pointahead-decorrelation-framing`).
- **DD-7. The wired-versus-available status lives in three places**
  (CLAUDE.md, examples/andrews/README.md, the crosscheck). Each change
  needs three edits; consider one home.
- **DD-8. Crosscheck batch 2 waits on owner input** (the "(owner to
  specify)" columns).

---

## External dependencies

- **X-1. The kernel repo holds uncommitted fixes.** `coupled_flux.py` in
  D:\repos\my_analysis_modules is untracked; the Dios-verified fixes sit in
  the working tree only. Commit them.
- **X-2. KR-24: the kernel keeps three wrong constants**
  (general_atmospherics.py:23 uses 0.54 / 1.22 / 0.509; the book uses
  0.49 / 1.11 / 0.51). The olb half is fixed; the kernel half is open.
- **X-3. Gap 8 and the strong aperture-averaged index need a second
  source** (a literature task, not a code task).

---

## Big reference buckets (do not work these as one item)

- **Crosscheck Table 2: 166 book-capability gap rows, 84 at priority P1**
  (docs/andrews-crosscheck.md:365–543). Named still-open rows: G-20, G-41,
  G-42, G-71, G-75, G-97, G-98, G-125, G-133, G-140, G-151. Mine per work
  package.
- **Crosscheck Table 3: about 40 book constants absent from olb** (partly
  stale, see DD-4).
- **Inherent limits (recorded, not fixable from this book):** the
  plane-parallel slant atmosphere (zenith-angle limit 60 deg); the
  two-scale large-scale branch does not reduce to Kolmogorov.
