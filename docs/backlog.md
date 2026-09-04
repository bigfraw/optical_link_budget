# Backlog — the unimplemented and unwired work

Date: 2026-08-26. Refreshed: 2026-09-04 (a full validation sweep against
commit 9316ba6; each item was checked against the code, marked DONE or PARTLY
DONE where the code closed it, and its line references were corrected; no item
was removed). Two repository sweeps made this list: one sweep of the
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

1. **Gap 2 is DECIDED and WIRED (2026-08-27): the pre-compensated uplink
   fade comes from FAST.** No trustworthy closed form exists; the model of
   record is the fidelity-1 FAST route, and the wiring is DONE:
   `uplink_fast_term` plus `uplink_budget(fidelity=1)` (the default for a
   pre-compensated scenario). See 0-W1 for the decision record and 1-2 for the
   remaining FAST limits.
2. **HIGH (owner-flagged, 2026-08-27) — stop the reliance on the
   `DEFAULT_HS` 20-layer array. STEP 1 DONE (2026-09-02).** HV5/7 is a
   continuous profile; the turbulent planner now takes a callable and
   integrates it, and the default fidelity-2 SPACE budgets use the continuous
   plan. Step 2 (the fidelity-0/1 `hs`-array modules move to callables) is
   still open. See 2-I2.
2a. **VERY IMPORTANT — DO SOON (owner-flagged 2026-09-02): the tail-convergence
   study.** The continuous planner from step 1 now lets the near-ground
   resolution be dialled independently. Run the study WP7 could not: does the
   deep SMF fade tail (p5/p1 — the link availability margin) converge as the
   near-ground `Cn2` is resolved? See 2-I2T.
3. **DONE — the turbulent screen-count floor `min_screens`.** Work package 7
   resolved it. See 2-N1.
4. **DONE — Gap 3, thread the beam curvature f0 into the Fried call site.**
   Closed 2026-08-27. See 0-W2.
5. **PARTLY DONE — the stale docs that contradict the code.** Cheap, and they
   mislead every later session. DD-3 and DD-6 are closed and DD-2 is half
   closed (2026-09-04 sweep). DD-1, DD-4, DD-5, DD-7, DD-8 and DD-9 stay open.
   See the documentation-debt group.
6. **DONE — olb is self-contained.** The `my_analysis_modules` kernels are
   vendored into olb and `_deps.py` is deleted (2026-08-28). See X-1. The
   kernel-repo commit and the KR-24 constants remain a concern for that repo
   only, not for olb.
7. **The owner decisions.** `downlink_budget` default (0-W5), and the
   FAST-versus-field reference model — whether wave optics ever becomes a
   DEFAULT. The turbulent fidelity-2 Term is now WIRED as opt-in (2-W1, done
   2026-08-28); the reference-model choice stays open. Next owner-requested
   step: an automatic fidelity selector. See 0-W5, 2-W1.
8. **Closed since the last refresh (2026-09-04 sweep):** 0-N1 (the shared
   `rytov_weak` gate), 0-N2 (the Rytov constants), 1-6 (the aperture-averaged
   lognormal certification), I-3, I-4, I-5, DD-3 and DD-6. Partly closed: 1-3
   (fidelity 2 is the manual strong-regime route; the AUTOMATIC selection
   remains), 2-N3 (the grid ideas are measured and buried; the co-moving
   long-path schedule remains), I-2 and DD-2.

---

## Fidelity 0 — analytic

### Wiring steps (built, no budget consumes it)

- **0-W1. Gap 2 — DECIDED 2026-08-27: the pre-compensated uplink gets NO
  analytic scintillation Term.** The owner rejected the earlier plan to wire
  `andrews.paths.uplink_scintillation_index(tracked=True)`. Three reasons:
  the tracked form removes the wander fully, which is a perfect tilt
  correction, and the beacon tilt decorrelates from the uplink path over the
  point-ahead angle; the same decorrelation applies to each higher corrected
  order; and a decorrelated correction reshapes the beam, so the Ch. 12
  normalisation by the vacuum-diffraction beam radius breaks. The tracked
  form is OPTIMISTIC, not a bound. No trustworthy closed form exists; the
  literature computes this case numerically. Resolution: the beacon + AO
  budget stays phase-only and mean-only, returns with LOUD flags
  (`NO SCINTILLATION, NO FADE`, plus the new extended-Marechal limit flag at
  sigma2 > 1 rad^2, T. S. Ross, DOI 10.1364/AO.48.001812), and stays useful
  for the geometric-only path (turbulence=False). The model of record is the
  fidelity-1 FAST Monte Carlo with the point-ahead offset, and that wiring is
  DONE (2026-08-27, same day): `uplink_fast_term` in
  olb/models/fast.py computes DTHETA from `geometry.point_ahead_rad`
  and returns the pure turbulence penalty with a real fade;
  `uplink_budget(fidelity=1)` (the default for a pre-compensated scenario)
  consumes it, and `fidelity=0` keeps the analytic phase-only pair as the
  no-dependency fallback. The FAST 0.1.7 Monte Carlo is direction-agnostic
  (the commented `PROP_DIR="up"` branches touch only the analytic budget olb
  does not read), so the reciprocity mapping needs no static-floor
  correction. Docs updated 2026-08-27: paths.py docstring, CLAUDE.md,
  docs/physics.md, docs/api-budget.md.
- **0-W2. Gap 3: thread the curvature f0 into the Fried parameter — DONE
  (2026-08-27).** The terrestrial SMF coupling call site in
  olb/models/coupling/terrestrial.py now reads the launch curvature from the
  transmitter divergence through the new `olb.beam.launch_curvature` and passes
  it as `f0` to `gaussian_fried_parameter_profile`. So a deliberately diverged
  beam gets its own r0. `launch_curvature` is ONE shared implementation: the
  Dios scintillation feed in olb/turbulence/uplink_flux.py now calls it too
  (it held a private copy of the same f0 algebra). The self-check asserts a
  diverged beam gives a larger r0 and a smaller coupling loss than a collimated
  one. STILL OPEN (a tidy-up, not a gap): the single-path
  `gaussian_fried_parameter` keeps the collimated signature; the budgets use
  the profile form, which is fixed. Docs: docs/physics.md GF-01 note updated;
  crosscheck GF-01.
- **0-W3. Gap 1: the aperture angle-of-arrival tilt feeds no Term.**
  `andrews.structure.angle_of_arrival_variance` is built, and
  `olb.turbulence.angle_of_arrival.aperture_angle_of_arrival_variance`
  delegates to it. No coupling Term adds contribution C: the terrestrial
  walk-off Term reads `wander_arrival_angle_variance` only
  (olb/models/coupling/terrestrial.py:41). So the received tip-tilt is a lower
  bound. Note conflict C-04: olb holds two tilt conventions (gradient 0.174
  vs the Noll route in ao.py); a caller that adds them must say which.
- **0-W4. Gap 6 and Gap 7: `l0`/`L0` and the temporal faces have no
  consumer.** No Term passes an inner or outer scale. No Term reads the
  Greenwood frequency, tau0, the fade rate, or the fade duration
  (andrews/temporal.py). The roadmap wants a tracking-bandwidth / servo-lag
  Term (README node NT7; also the deferred TODO at
  olb/models/coupling/terrestrial.py:254).
- **0-W5. `downlink_budget` still defaults to `model="lognormal"`.** The
  `model="auto"` selector exists and is opt-in (olb/links/downlink.py:470 and
  the budget keyword `scint_model="lognormal"` at :521). The switch moves the
  strong-regime total by several dB. OWNER DECISION.
- **0-W6. The Andrews aperture-averaging chain is unused.** The downlink
  Term keeps the numerical Airy-filter integral
  (`plane_wave_scintillation.aperture_averaged_scintillation_index`); the
  terrestrial Term keeps the Churnside fit
  (`aperture_averaging_factor_weak`, olb/links/terrestrial.py:148; optimistic
  5–13 % against the book's exact Eq. (60)).
  Conflict C-06: Ch. 12, Eq. (39) is the closed form that should supersede
  the integral when a work package moves the Term over.
- **0-W7. The K distribution and the lognormal-Rician PDF are unused**
  (andrews/distributions.py). Of the three spectrum labels in
  olb/assumptions.py, Tatarskii now labels
  `plane_wave_scintillation.aperture_averaging_factor_weak_inner`
  (olb/turbulence/plane_wave_scintillation.py:470), but no Term calls that
  function; the exponential and the modified labels stay inside
  andrews/spectra.py. So no Term reads any of the three.
- **0-W8. The Andrews wander route stays on the shelf BY DECISION** (C-01
  closed via Belmonte; the budgets keep the Dios/Belmonte 2.07 kernel; the
  Andrews 7.25 route sits in andrews/wander.py for measurement only). Listed
  so nobody re-opens it.

### Physics gaps

- **0-P1. Laser guide star: focal (cone) anisoplanatism.** `LaserGuideStar`
  is a placeholder; `uplink_budget` raises (olb/links/uplink.py:648,
  olb/scenario.py:82). Its cone anisoplanatism differs from the point-ahead
  angular form.
- **0-P2. The short terrestrial retro link has no module.**
  `retro_space_budget` assumes a long slant range and independent legs. The
  book also gives a backscatter amplification and a coupled double passage
  that olb does not model (crosscheck RT-01..03).
- **0-P3. The strong / moderate aperture-averaged downlink index.** The
  gamma-gamma Term is a POINT receiver (the flag sits in `_gamma_gamma_term`,
  olb/links/downlink.py:235–247); its fade
  is deeper than the true aperture fade. The book gives no slant
  aperture-averaged strong index; a second source or a derivation is needed.
  (Owner earlier decided to keep point-receiver; see the memory
  `downlink-strong-aperture-decision`.)
- **0-P4. Gap 8: the annular (obscured) receive aperture.** No book source
  exists; one central-obscuration gap touches many Terms
  (olb/links/downlink.py:126, uplink.py:157,
  olb/models/coupling/downlink.py:142,
  olb/models/coupling/terrestrial.py:577). Needs a new reference.
- **0-P5. The terrestrial scintillation Term carries no pointing jitter.**
  It is on-axis (r=0) only; the jitter-into-beta fold sits in the Monte
  Carlo path only (TODO at olb/links/terrestrial.py:153). Fold the jitter
  into the off-axis radius, and add the mean-power loss Term.
- **0-P6. No tracking-loop bandwidth model.** The tip-tilt correction is
  all-or-nothing (olb/models/coupling/terrestrial.py:254). A finite servo
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
  airmass model breaks below 5 deg elevation (olb/models/extinction.py:102).
- **0-P10. Only one Gaussian transmit beam; `Transmitter.m2` is dead.** No
  model reads m2 (it is declared at olb/terminal.py:75 and read nowhere
  else). Flat-top beams and incoherent aperture diversity stay
  planned (README TB3, TB4).
- **0-P11. `eta_max(a)` assumes a uniform, flat-wavefront aperture — the
  CURVATURE half is CLOSED (2026-08-31), the ILLUMINATION half stays open.**
  The terrestrial coupling Terms now charge the received-beam curvature always:
  the true focus of a diverging received beam sits at `dz_curv = f^2/(R_rx - f)`
  BEYOND the focal plane (S. A. Self, DOI 10.1364/AO.22.000658), the detector is
  `dz_eff = defocus_m - dz_curv` from it, and the SMF mean penalty uses the
  defocus-aberrated closed form `smf_eta_defocused(a, c)` (Ruilier and Cassaing,
  DOI 10.1364/JOSAA.18.000143). What STAYS open is the ILLUMINATION half: a
  near-field received Gaussian tapers across the aperture, and `eta_max` assumes
  a uniform one. That error runs SAFE (the constant is then conservative). See
  docs/physics.md section 6a and validation/defocus/.
- **0-P15. The terrestrial SMF walk-off DISPLACEMENT response stays
  geometric.** The MEAN modal defocus penalty is modelled (0-P11), but
  `terrestrial_smf_walkoff_term` still answers a displacement with a
  two-Gaussian overlap of the defocused spot and the fibre mode. It does not
  model how the defocus phase reshapes the modal overlap against a
  displacement, so the walk-off fade is OPTIMISTIC off focus. The Term flags
  itself loudly when `defocus_m` is not zero. The full treatment is Ruilier and
  Cassaing, DOI 10.1364/JOSAA.18.000143; an MMF or fidelity 2 sidesteps it.
- **0-P16. The bidirectional wrapper models the DIVERGING launch only.**
  `olb/links/bidirectional.py` maps the collimator defocus to a transmit
  divergence through `|dz|`, and a `Transmitter` cannot hold a converging beam.
  So `dz > 0` (a converging launch) is outside the fidelity-0 model, and the
  wrapper gives it the mirror-image diverging divergence. One `dz` also drives
  BOTH sides of a monostatic terminal, so a deliberately diverged terminal pays
  `|dz| + dz_curv` of receive defocus. Both limits are in the module docstrings;
  a converging launch needs fidelity 2.
- **0-P12. The MMF NA gate is a flat factor.** No turbulence re-broadening
  of the focal spot, no mode-count saturation; `optimal_focus` is geometric
  (docs/physics.md:1354–1358).
- **0-P13. The uncorrected coupling curve extrapolates past its D/r0
  limit** and only warns (olb/models/coupling/downlink.py:150,
  olb/models/coupling/terrestrial.py:588).
- **0-P14. The strong-turbulence effective beam parameters are deferred.**
  `effective_beam_params` (Theta_e, Lambda_e) is coded and has no caller in
  `links/` or `models/`. The physics now lives in
  olb/turbulence/andrews/beam.py:163, and
  olb/turbulence/gaussian_fried.py:181 keeps the old signature and delegates
  to it (crosscheck GF-17). Deliberate, to match the Dios weak regime.

### Numerical and validation issues

- **0-N1. DONE (2026-08-29). TL-05: the terrestrial weak gate tested ONE
  criterion.** Ch. 5, Eq. (16), printed p. 140 needs `sigma_R^2 < 1` AND
  `sigma_R^2 Lambda^{5/6} < 1`. Commit 3470190 added the one shared
  beam-aware helper `rytov_weak(sigma2_R, Lambda)` in
  olb/turbulence/andrews/scintillation.py, which applies both conditions
  through the binding strength `sigma2_R * max(1, Lambda**(5/6))`. The
  terrestrial Term (olb/links/terrestrial.py) and the uplink coupled-flux
  path (olb/turbulence/uplink_flux.py, olb/links/uplink.py) call it, so a
  focused beam now trips the gate. The lognormal-PDF house rule is the
  distinct `LOGNORMAL_PDF_LIMIT` (Table 2 row G-20; crosscheck TL-05).
- **0-N2. DONE (2026-08-25). The strong-fluctuation parameter q carried 1.22
  one level too high.** Commit b0a25a1 corrected the constants. The point
  plane-wave index now lives in
  olb/turbulence/andrews/scintillation.py with the book values 0.49, 1.11 and
  0.51, and `plane_wave_scintillation.plane_wave_scintillation_index_closed`
  delegates to it. No source file holds 0.54, 1.22 or 0.509 any more
  (crosscheck ledger PW-09, KR-24).
- **0-N3. The budget adds dB losses / per-term p-quantiles.** Not a book
  form; it is a conservative upper bound (crosscheck RS-02, RS-03). Record
  the bound; consider a joint-sample check.
- **0-N4. The near-field truncation flag is not conservative.** The true
  value can sit above or below the far-field form
  (olb/models/gaussian_efficiency.py:186–193). The fidelity-2 vacuum layer is
  the verifier; consider an auto-route.
- **0-N5. The quasi-frequency has no upper limit of its own** (`b_2` grows
  as `f_max^{1/3}`; olb/turbulence/andrews/temporal.py:608). Any caller must
  set the band from the detector or an inner scale. Pairs with the
  inner-scale memory.
- **0-N6. The fade rate and the fade time have no external check** — the
  book gives no worked example; the faces are checked against internal
  identities only.

### Blocked by the source (documented refusals — need a second source)

Each of these raises `NotImplementedError` with a citation; the book prints
no form, and olb guesses no coefficient. Full list: docs/physics.md:917–962.
The path forward for each is a second reference or a derivation.

- The Gaussian two-scale STRONG branch (the unresolved eta_X of Ch. 9,
  Eq. (109)) and its aperture chain (andrews/scintillation.py:609,
  andrews/aperture.py:422).
- The weak Gaussian-beam aperture flux variance (Ch. 10, Eq. (78) is a
  numerical double integral; andrews/aperture.py:439).
- The convergent-beam Rytov variance (Theta0 < 1;
  andrews/scintillation.py:463).
- Inner/outer scale on any Ch. 12 slant form (andrews/paths.py:485); the
  full Gaussian-beam downlink index Ch. 12, Eqs. (36)–(37) (row G-133). The
  aperture-averaged STRONG downlink index is refused in the same module
  (andrews/paths.py:586; see 0-P3).
- The temporal spectrum with a finite inner or outer scale; the strong
  spherical/Gaussian temporal spectra; the weak-only aperture-averaged
  temporal spectrum (andrews/temporal.py:538, :550, :560).
- The Gaussian row of the modified spectrum (App. III, Table III;
  andrews/structure.py:345).
- The K distribution below sigma_I^2 = 1 (andrews/distributions.py:528; it
  raises `ValueError`, not `NotImplementedError`).

---

## Fidelity 1 — statistical / Monte Carlo

- **1-1. The terrestrial statistical coupling fade — RESOLVED at fidelity 2
  (2026-08-28).** The default (fidelity 0) SMF Term is still mean-only and LOCKS
  the budget out of a fade margin (olb/models/coupling/terrestrial.py). But
  `terrestrial_budget(fidelity=2, wave=...)` replaces it with the two wave-optics
  Terms (vacuum-optics + turbulence), which carry a real fade, so the budget then
  gives a fade margin. See 2-W1. (FAST is far-field only; the near-field Gaussian
  beam needs the split-step model, which is what fidelity 2 uses. Fidelity 1 is
  UNAVAILABLE for a terrestrial link and raises.) STILL default fidelity 0 by
  owner decision (the field reads less coupling loss than the incumbent; the
  reference-model gap of 2-W1 stays open, although the terrestrial MMF part of
  that gap fell to about 1.2 dB once the received curvature was charged, see
  0-P11 and 2-W1). A PROPOSED way to fill the empty fidelity-1 rung for a
  terrestrial link, without FAST, is the calibrated lognormal draw of 1-8.
- **1-2. FAST limits NT1–NT4 — CARRIES GAP 2 (see 0-W1); the uplink entry
  point is DONE (2026-08-27).** `uplink_fast_term` in
  olb/models/fast.py is the pre-compensated uplink model of record:
  it computes `DTHETA` from `geometry.point_ahead_rad` (NT1 closed for the
  uplink), sets the numeric launch waist, and returns the pure turbulence
  penalty by reciprocity (Shapiro, DOI 10.1364/JOSA.61.000492; Farley,
  DOI 10.1364/OE.458659) with the weak-regime amplitude gate.
  `uplink_budget(fidelity=1)` (the default for a pre-compensated scenario)
  consumes it. Still
  open from the first cut: the DOWNLINK Term keeps `DTHETA=[0,0]`, which is
  correct for a receive-side AO that senses the same downlink beam; scalar
  elevation only (both Terms); no obscuration in the coupled flux or the
  mean-only fibre model; no tip-tilt wander removal; the FAST servo/WFS
  defaults (DSUBAP=0.02 m, TLOOP=TEXP=1 ms, ALIAS on) pass through
  unreviewed — override with `fast_params`.
- **1-3. Strong-fluctuation routing.** The uplink and terrestrial links
  only WARN when the weak-fluctuation limit is exceeded. The downlink can route
  to gamma-gamma, but that route is OPT-IN: `downlink_budget` still defaults to
  `model="lognormal"`, and only `model="auto"` selects the gamma-gamma Term
  above `LOGNORMAL_PDF_LIMIT`. Fidelity 2 gives a MANUAL strong-regime route for
  a terrestrial link and for an uncorrected uplink; a pre-compensated uplink has
  no fidelity-2 route (it raises). What remains is the AUTOMATIC selection: route
  the other links to a strong-regime model (a Monte Carlo or the fidelity-2
  layer; see the memory `strong-fluctuation-numerical`).
- **1-5. Validate the FAST point-ahead residual against the analytic Stone
  decorrelation — DONE (2026-09-02). VERDICT: MATCH, both routes validated.**
  The study is `validation/fast_stone_pointahead/` and the entry of record is
  physics.md Section 9j. With the servo off (`TLOOP=0`, `TEXP=0`, zero wind)
  the PAOLA filter reduces exactly to `2 - 2cos(delta_r . kappa)` (verified to
  9e-16), and at MATCHED MODE SETS the two routes agree to about 5 % across
  the full sweep (point-ahead 0.25x to 2x, ZMAX 1 to 66 plus the exact-zero
  uncorrected anchor, elevation 30 to 90 deg); the fitting sides agree to
  0.6 %. THE MODE SETS MATTER: the FAST modal
  mask keeps the piston and the tilts, so its analytic partner is Stone
  `remove='none'`; the production pairing (against `piston_tilt`) reads 3.5x
  and the whole factor is that convention. The backlog first reading is
  reproduced (3.04 vs 3.79 dB) and the Term-level gap decomposes into the mode
  set, the FAST auto-grid truncation, and the Marechal-vs-Monte-Carlo mapping.
  Two OPEN FAST cautions, measured: `sim.aniso_servo_error` leaks
  `mask(1-mask)` of the uncorrected band (0.061 rad^2 at theta = 0, truth 0),
  and the FAST grid misses 29 to 48 % of the whole-plane Kolmogorov residual
  (the shipped Term's auto grid, `df` = 3.11 rad/m, misses more); the missing
  scales sit far above the aperture, so the coupled-flux effect is damped, and
  that damping is not quantified. The original recipe follows. The two routes compute the
  same quantity: the residual phase of a finite-aperture pre-compensation
  that measures off-axis by the point-ahead angle. FAST integrates the
  PAOLA aniso-servo filter over the corrected spatial-frequency mask and
  exposes the number (`sim.aniso_servo_error`, per layer;
  `ao_power_spectra.G_AO_PAOLA`). `uplink_point_ahead_term` sums the Stone
  modal decorrelation residual 2 sigma_n^2 (1 - rho_n) over the corrected
  Zernike orders (DOI 10.1364/JOSAA.11.000347). Compare them at matched
  conditions: the same Cn2 profile, aperture, and corrected order; the
  servo and sensor effects off (`TLOOP=0`, `TEXP=0` reduces the PAOLA
  filter to the pure two-path form 2 - 2cos(delta_r . kappa); `ALIAS`
  False, `NOISE=0`). A match validates both routes. A mismatch measures the
  real difference between the two finite-aperture treatments: PAOLA masks
  spatial frequencies, Stone projects Zernike modes. Also compare the
  fitting side (`sim.fitting_error` against the Noll residual of
  `uplink_fitting_term`), and then the full Terms (the FAST mean against
  the Marechal sum — first reading at AO(60), 60 deg, 1.5 m: 3.1 dB against
  3.79 dB). Sweep the point-ahead angle and the corrected order before the
  FAST Term is trusted at other operating points. (No validation script for this
  comparison exists in `validation/` yet.)
- **1-4. The Dios duplicate (ponytail DEBT).** The analytic
  beam_wave_scintillation path and the coupled-flux MC duplicate the same
  equations; the jitter correction sits in the MC path only (the `ponytail: DEBT`
  comment in olb/turbulence/beam_wave_scintillation.py). Converge on ONE
  implementation; do not make a third copy.
- **1-6. DONE (2026-09-01). Certify the aperture-averaged lognormal power draw
  against fidelity 2 (owner-flagged 2026-08-29).** The `--full` certification run
  closed this item. The script is
  `validation/lognormal_certification/lognormal_certification.py` (commit
  927c952), and the certification of record is that folder README plus
  docs/physics.md Section 9e. The verdict is PASS to the 1 percent fade for a
  weak horizontal path. The record below stays for the history.
  The question: in WEAK turbulence,
  does the received-power distribution under APERTURE AVERAGING still follow the
  lognormal that the cheap analytic route assumes? The route takes the
  aperture-averaged index sigma2_P = A sigma2_I and draws from a lognormal (see
  the `sigma2_P = A * sigma2_I` step in olb/links/terrestrial.py and the
  lognormal build in olb/links/downlink.py). But an
  aperture integrates a correlated lognormal field, and a sum of lognormals is
  NOT lognormal: as D grows the power drifts toward Gaussian (thinner tails), and
  for a finite Gaussian beam any beam wander adds a pointing tail the single
  index cannot describe. The test: run the fidelity-2 split-step Monte Carlo
  (olb/waveoptics/turbulence), histogram the aperture-collected power across
  trials, and compare the empirical PDF (and, crucially, the deep-fade tail
  quantiles) against the lognormal built from the analytic mean + sigma2_P at the
  matched sigma2_I. Sweep D/rho_0 (point aperture to strong averaging) and the
  beam geometry (collimated, diverged, focused). A pass certifies that the
  'easy' analytic weak-turbulence calc can generate a trustworthy power
  distribution; a fail bounds where it may be used and points to the composite
  (lognormal x pointing) or a direct empirical sampler.
  NOTE ON LADDER LABELLING: producing a power DISTRIBUTION (draws) from an
  analytic sigma2_I is STATISTICAL, so by the olb ladder it is FIDELITY 1, even
  though the sigma2_I itself is an easy fidelity-0 analytic quantity. The owner
  confirms this labelling is intended and fine. The fidelity-2 sim here is the
  reference that certifies the fidelity-1 draw; it is not a new budget path.
  See the memory `aperture-averaged-lognormal-certification` and the discussion
  of aperture averaging and beam wander in the C-05 / TL-05 thread.
  UPDATE (2026-09-01): the certification script EXISTS at
  `validation/lognormal_certification/`. It sweeps `D/rho_0` from 0.20 to 7.89 on
  one firmly weak 2 km horizontal path (`sigma_R^2 = 0.21`), for a collimated and
  a diverged launch, and it reports the index, the fade quantiles and the skew of
  `ln P` apart, so an INDEX error and a SHAPE error do not mix. QUICK-MODE first
  reading (150 trials, `rapid` preset): the lognormal FAMILY HOLDS -- with the
  index refit to the measured value every case agrees inside 0.12 dB at the 5 %
  fade, and the skew of `ln P` stays in [-0.38, +0.19] with no trend against D,
  so no drift to a Gaussian power and no pointing tail is visible in this band.
  The fault is the INDEX: near `D/rho_0 = 3` the analytic `sigma2_P` reads 2.1x
  (collimated) to 2.7x (diverged) LOW, so the analytic fade is optimistic by 0.24
  to 0.26 dB at the 5 % fade. The analytic on-axis index also takes the waist
  only, so it gives the collimated and the diverged launch the SAME number while
  the field does not.
  UPDATE (2026-09-01, second pass): the script now SPLITS that index error in
  two, and the answer is the FILTER. One propagation for each trial
  (`propagate_turbulent_field`) now serves the whole aperture sweep, so the point
  estimator and every diameter read the SAME atmosphere; a matched-seed check
  against `propagate_turbulent_scenario` on the shared grid agrees BIT FOR BIT,
  and the quick run fell from about 9 minutes to about 2.5 minutes. The POINT
  index (the mean irradiance in an 8 mm on-axis disc, 0.14 of the Fresnel scale)
  reads `sigma2_I` = 0.0615 collimated and 0.0676 diverged against the analytic
  0.0744, so the Dios on-axis form is 10 to 20 % HIGH -- a modest error. The
  Churnside filter is the fault: `A_eff = sigma2_P_sim / sigma2_I_sim` runs 1.4x
  the analytic A at `D/rho_0` = 1, 2.6x to 2.9x at `D/rho_0` = 3, and 2.5x at
  `D/rho_0` = 7.9. So A OVER-AVERAGES across the whole band, and the point index
  partly HIDES it (the two errors pull in opposite directions). Two more
  findings. (a) The D = 40 cm COLLIMATED column of the first reading was a
  BEAM-FILLING artifact, not physics: `w(L)` = 19.7 cm there, so the aperture
  catches `eta_fill` = 0.87 of the beam and measures near-total power, which
  fluctuates little (A ratio 0.53). The diverged launch at the same diameter has
  `w(L)` = 40.4 cm, `eta_fill` = 0.39, and it behaves like every other unfilled
  case (A ratio 2.54). The script now computes `eta_fill` for every case and
  FLAGS a case past 0.5 as BEAM-FILLING-LIMITED, in the log, the JSON and the
  figure. That is backlog 2-N2 measured. (b) The ABSOLUTE impact is small: the
  fade spread falls from 1.07 dB (D = 1 cm) to 0.06 to 0.13 dB (D = 40 cm), so
  the WORST relative index error (2.9x) moves the 5 % fade by 0.26 dB only.
  DONE for this path (2026-09-01, the `--full` run: 1500 trials for each
  launch, `standard` preset, about 2.8 hours). VERDICT PASS to the 1 % fade:
  the refit lognormal agrees inside 0.128 dB at the 5 % fade and 0.210 dB at
  the 1 % fade in every case; the whole analytic route inside 0.289 dB and
  0.413 dB (both worst cases the diverged 15 cm receiver); the skew of `ln P`
  sits near -0.2 with no trend. The full run also RETIRES one quick-mode
  finding: the 10 to 20 % point-index bias was a `rapid`-preset artifact -- at
  `standard` the analytic `sigma2_I` reads only 3 to 4 % high. The filter
  fault stands, milder: `A_eff/A` about 1.2 at `D/rho_0 = 1`, 1.8 to 2.5 at 3,
  2.4 at 7.9 (unfilled). The certification of record is in the folder README
  and physics.md Section 9e. WHAT REMAINS is the 1-8 gate (b) sweep, not this
  item: a focused launch, a stronger Cn2, and a longer path.
- **1-8. Terrestrial fidelity 1 = the calibrated lognormal draw (PROPOSED
  2026-09-01).** A terrestrial link has NO fidelity-1 rung: FAST is far-field
  only, so `terrestrial_budget(fidelity=1)` raises (1-1). The proposal fills that
  rung with a CALIBRATED DRAW instead of a new analytic model:
  1. Run a SHORT fidelity-2 batch for the scenario (approximately 100 to 200
     trials, minutes on the fast `ScreenFactory`) and measure the received-power
     mean and the aperture-averaged index `sigma2_P` empirically.
  2. Draw the fade from the lognormal REFIT to those two measured moments.
  3. Cache the calibration (a `Campaign`,
     `olb/waveoptics/turbulence/campaign.py`, keeps a run on disk), so a sweep or an optimiser
     pays the simulation one time and the draw after that.
  The exact static Terms (extinction, geometric, pointing) STAY analytic, and the
  analytic scintillation Term stays as the free sanity anchor and the regime gate.
  RATIONALE: the quick-mode run of `validation/lognormal_certification/`
  (2026-09-01, see 1-6) certified the lognormal FAMILY in the weak regime. With
  the index refit to the measured value, every case agrees inside 0.12 dB at the
  5 % fade and 0.30 dB at the 1 % fade out to `D/rho_0 = 7.9`, with no skew trend,
  even with the beam wander fully uncorrected. So the weak link is the analytic
  FEED, not the distribution: the Churnside plane-wave `A` and the waist-only
  Dios point index misread the index by 2.1x to 2.7x near `D/rho_0 = 3`, and the
  analytic chain cannot tell a collimated launch from a diverged one. A measured
  mean and index remove that feed.
  BEAM-FILLING CAVEAT: the `D = 40 cm` column of the quick-mode sweep is
  BEAM-FILLING-limited (the collimated launch fills approximately 87 % of the
  aperture, the diverged launch approximately 39 %), so the aperture holds the
  beam. That is the known failure of the analytic averaging factor (2-N2), and
  that column therefore does NOT test the averaging filter.
  GATES before this becomes a default: (a) the `--full` certification run of 1-6
  (DONE 2026-09-01, PASS to the 1 % fade; see 1-6 and physics.md 9e); (b) at
  least two more scenarios (a stronger `Cn2`, a longer path); (c) the 2-W1 owner
  reference-model decision. NOTE on gate (b): the lognormal SHAPE is certified
  in the WEAK band only. In strong fluctuation even the REFIT lognormal can
  fail (the fade PDF moves toward gamma-gamma), so the stronger-Cn2 case of
  gate (b) must certify the SHAPE again, not only the index — a draw from a
  wrong family stays wrong at any calibration. This is a PROPOSED design,
  approved for the backlog only. It is NOT built, and it must not be started
  yet. See the memory `terrestrial-calibrated-draw-plan`.
- **1-7. REFERENCE for the residual scintillation of a pre-compensated
  uplink.** Gap 2 (0-W1) decided that NO trustworthy analytic scintillation
  form exists for a beacon + AO pre-compensated ground-to-space beam, so the
  budget stays mean-only there and the fade comes from the fidelity-1 FAST
  Monte Carlo. Look at 'Phase estimation at the point-ahead angle for AO
  pre-compensated ground to GEO satellite telecoms' for a treatment of the
  RESIDUAL scintillation fluctuations that survive the pre-compensation. It
  can inform a future residual-scintillation model or a cross-check of the
  FAST point-ahead residual (see 1-5). A reading task, not a build task.

---

## Fidelity 2 — wave optics

- **2-AO. VERY IMPORTANT — the fidelity-2 sims model NO adaptive optics. This
  includes tip-tilt correction (2026-08-31).** The split-step layer
  (`olb/waveoptics/turbulence/`) propagates the field through the raw phase
  screens and reads the receive plane with NO wavefront correction applied. So a
  fidelity-2 budget shows the UNCORRECTED atmosphere on EVERY link:
  - No tip-tilt (beam-wander / angle-of-arrival) removal. A tracked terminal
    removes the received tilt in the real system; the sim does not, so the
    fidelity-2 fade and the coupling loss are PESSIMISTIC for a tracked link,
    and the walk-off is fully counted.
  - No higher-order AO (the deformable-mirror correction of the residual phase).
  - No uplink pre-compensation. A beacon + AO uplink (the fidelity-1 model of
    record, `uplink_fast_term`) has NO fidelity-2 equivalent: `uplink_budget`
    RAISES at `fidelity=2` for a pre-compensated scenario, and the reciprocity
    route carries no point-ahead correction either (see 2-P4).
  This means a fidelity-2 result is directly comparable to the fidelity-0/1
  UNCORRECTED case only. It is NOT a like-for-like reference for the
  AO-corrected fidelity-1 terms (`smf_fast_term`, `uplink_fast_term`) or for any
  tracked terrestrial coupling Term until AO is modelled in the field solve.
  The fix is a correction stage in the split-step receive path: at minimum a
  tip-tilt removal (subtract the measured wavefront tilt, or the centroid shift),
  then a modal / zonal higher-order correction, and for the uplink a
  pre-compensation phase applied at the launch plane crossed with the point-ahead
  decorrelation. Each stage changes budget numbers, so each is an owner-gated
  step. NOTE (2026-09-02): `olb/waveoptics/camera.py` `spot_metrics` MEASURES the
  spot centroid, but it is diagnostic only. No runner and no budget uses that
  measurement to remove the tilt. Pairs with 2-P4 (point-ahead anisoplanatism in
  the reciprocity route) and the reference-model gap of 2-W1 (the field reads
  less coupling loss than FAST — part of that gap is this missing correction, so
  the two are NOT yet comparable).
- **2-W1. Fidelity-2 is WIRED whole-path via `fidelity=0|1|2` (2026-08-28,
  branch `fidelity2-budget-wiring`). BOTH the turbulent split step AND the
  vacuum core are now consumed.** A fidelity-2 budget shows TWO Terms: a
  DETERMINISTIC `waveoptics_vacuum_term` (the full no-turbulence loss from
  `propagate_scenario`) and a STOCHASTIC `waveoptics_turbulence_term` (the fade);
  only the analytic extinction and pointing Terms stay. The caller precomputes
  both records once with `run_fidelity2` -> `Fidelity2Bundle` and passes `wave=`;
  the budget never runs the sim. The old per-component knobs
  (`smf_fidelity="waveoptics"`, `scint_model="montecarlo"`, budget-level
  `smf_fidelity`/`precomp_fidelity`) are REMOVED at the budget level and folded
  into `fidelity`. It closes 1-1 (terrestrial fade lock), and addresses 0-P11 /
  Gap-3 and 0-N4 (the field solve sidesteps the analytic `eta_max` and the
  near-field truncation flag). All in olb/models/waveoptics.py (a FIDELITY-named
  module at the `models/` level; the coupling package re-exports the two coupling
  faces); see `examples/waveoptics/budget_wiring.py`. STILL OPEN and
  OWNER-GATED: whether wave optics ever becomes a DEFAULT. That is the
  reference-model gap — the field reads LESS fibre coupling loss than the
  incumbents (0.7–2.9 dB less than FAST, ~2.5 dB less than the terrestrial
  analytic Term). FOLLOW-UP (owner-requested 2026-08-28): an AUTOMATIC fidelity
  selector, the way `model="auto"` picks a distribution.

  UPDATE (2026-08-31): a SPACE link now takes the ANALYTIC geometric Term by
  default (`run_fidelity2(vacuum="analytic")`, so `wave.vacuum` is None), and
  only a TERRESTRIAL link keeps the wave vacuum Term. `run_fidelity2(vacuum=
  "wave")` opts a space link back in.

  QUANTIFIED for the TERRESTRIAL MMF leg (2026-08-31, validation/defocus/). Most
  of the old terrestrial disagreement was the missing received-curvature defocus,
  not a wave-optics gap: with the curvature charged (0-P11) the fidelity-0
  against fidelity-2 MMF coupling gap of the report scenario falls from about
  7 dB to about 1.2 dB (8.54 dB analytic against 7.08 dB field, at 1550 nm,
  L = 5 km, w0 = 0.02 m, D = 0.2 m, 25 um core, f = 4.524 m). The residual is the
  Airy-versus-Gaussian SPOT SHAPE: the truncated pupil makes an Airy pattern
  whose slow rings a Gaussian spot model omits. The gap is NOT closed: the SPACE
  half (downlink SMF against FAST, the 0.7–2.9 dB rows above) is untested against
  this correction, and the residual ~1 dB spot-shape term is not chased.
- **2-W2. The fidelity-2 SMF path ignores `defocus_m` — DONE (2026-09-04).**
  `olb.waveoptics.smf.coupling_efficiency` now takes `defocus_m` and
  `focal_length_m`. It multiplies the received field with the SAME quadratic
  pupil phase the multimode leg uses, `exp(-i*pi*defocus_m*rho^2/(lam*f^2))`,
  which moved into the shared helper `olb.waveoptics.mmf.defocus_phase`. So the
  two fidelity-2 legs read ONE defocus convention: a DIVERGING received beam
  couples best at a POSITIVE `defocus_m`. `olb.waveoptics.run._smf_eta` resolves
  the focal length (an explicit `SMF.focal_length_m` wins, else
  `SMF.optimal_focus` gives `f = pi*(D/2)*w_m/(lambda*1.12)`) and it feeds
  `SMF.defocus_m` to every call site: the vacuum run, the turbulent trials, and
  the multi-arm `detector_etas` path. A defocus with no focal length raises. The
  default `defocus_m=0.0` keeps every old number bit-identical. The `smf.py`
  self-check now measures the field overlap against the closed form
  `smf_eta_defocused(a=1.12, c)` (D = 0.1 m, f = 0.5 m): 0.8145/0.8145 at
  c = 0, 0.7537/0.7537 at c = 1, 0.5936/0.5936 at c = 2 and 0.2076/0.2076 at
  c = 4. So the aberrated single-mode closed form HAS a field reference now.
- **2-W3. The power-to-pixel-brightness conversion for the Camera detector is
  NOT built (owner-deferred 2026-09-02).** The pieces exist: the `Camera`
  dataclass (olb/terminal.py: `pixel_pitch_m`, `n_pixels`, `focal_length_m`,
  `defocus_m`), `olb.waveoptics.camera.camera_image` (the focal spot binned onto
  the pixel grid, normalised to the collected power) and `spot_metrics`, and the
  vacuum route `run_fidelity2(turbulence=False)` -> budget. The missing step is
  the scale from the budget received power to watts (or photons, or counts) per
  pixel. A demonstration script did this with one multiplication and was REMOVED
  on purpose: the owner wants a HOLISTIC camera model, not a bare scale. That
  model adds the camera-specific parameters (for example the quantum efficiency,
  the integration time, the read noise, the dark current, the full-well and
  saturation limit, the bit depth) so the output is a real detector signal, not
  an ideal power map. Design it in one pass with the owner before any wiring;
  each parameter changes what "pixel brightness" means.
- **2-W4. Retire `cache.py` in favour of `Campaign` — DONE (2026-09-04).**
  `olb/waveoptics/turbulence/cache.py` (the P4 opt-in disk cache of scalar
  runs, block sub-seeds, no field) is DELETED, with its self-check and
  `validation/waveoptics_speed/cache_check.py`. `Campaign` replaces it: its
  blocks are bit-identical slices of ONE seeded native run (the runner's
  `start_index`), and it stores the receive-field patch. The one piece the
  campaign imported, the `cache_key` content fingerprint, moved unchanged to
  `olb/waveoptics/turbulence/fingerprint.py`, so an existing campaign manifest
  still matches. No budget ever called the cache. The measured numbers of the
  cache stay as a record in docs/waveoptics-efficiency-plan.md Section 8.
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
  `_equal_weight_groups` and `_merge_layers` are in
  olb/waveoptics/turbulence/sampling.py.
- **2-I2. Continuous Cn2 profiles — drop the `DEFAULT_HS` crutch (HIGH,
  owner-flagged 2026-08-27). STEP 1 DONE (2026-09-02).** HV5/7 and the other
  Cn2 models are continuous functions; the 20-layer `DEFAULT_HS` array is a
  hand-made discretisation, and it leaks into the physics wherever a decision
  reads the grid instead of the profile. Work package 7 removed the worst leak
  (the screen count), but the screen PLACEMENT still came from the array. The
  post-WP7 matched-seed measurement (the WP7 note in docs/schmidt-crosscheck.md)
  sharpened the question: the bottom screen HEIGHT is a null, and the live
  variable is whether the near-ground `Cn2` is spread over many thin screens or
  lumped into one — a resolution question that only the continuous ground layer
  can answer. The same measurement showed a hint that the placement moves the
  deep SMF fade tail (about 2 dB at p5, direction consistent, NOT resolved above
  the Monte-Carlo noise at 200 trials; the mean is a null, 0.23 +/- 0.48 dB).
  The change, in two separate steps:
  - **(1) DONE (2026-09-02).** `turbulent_grid` and `_plan_space` in
    olb/waveoptics/turbulence/sampling.py now take a callable `cn2(h)` and
    INTEGRATE it. The DEFAULT (no hs/cn2_profile) builds the site HV5/7 callable
    and integrates it; `DEFAULT_HS` is now the fallback for an explicit array
    caller ONLY (`_plan_space_array`, frozen behaviour). The continuous planner
    (`_plan_space_continuous`) places screens by EQUAL RYTOV WEIGHT
    (N = max(min_screens, ceil(sigma2_R_total / cap)) equal-weight slabs) with a
    Cn2-weighted centroid per slab. It is grid-free by construction: a finer
    internal integration grid does not move the plan. Validated in the module
    self-check — the continuous r0 matches a fine-grid analytic to <1% (the
    coarse 20-layer trapezoid is biased ~2% low in r0, the crutch bias), every
    profile moment holds inside 1% (Schmidt Eq. (9.65)), and the sampling stays
    good at 10/30/90 deg. `cn2` and `h_top_m` are threaded through
    `propagate_turbulent_scenario`, `propagate_turbulent_field`, `run_waveoptics`,
    `run_fidelity2`, and the opt-in cache (the cache fingerprints the callable by
    sampling it). The default fidelity-2 SPACE budgets now use the continuous
    plan; only the screen PLACEMENT moves (the total turbulence is conserved),
    and the mean is validated flat. The Gauss-quadrature screen placement
    (tracker candidate, S-22) is now implementable on this base. FOLLOW-UP: the
    TAIL-CONVERGENCE STUDY — see 2-I2T, flagged VERY IMPORTANT / DO SOON by the
    owner on 2026-09-02.
  - **(2) OPEN, LATER, and separately:** the fidelity-0/1 modules that integrate
    over `hs` arrays (slant extinction and scintillation, uplink flux, FAST) move
    to callables; that step is wide, mechanical, and must move no numbers. Still
    NOT built.
- **2-I2T. The tail-convergence study (VERY IMPORTANT — DO SOON, owner-flagged
  2026-09-02).** Now that the continuous planner (2-I2 step 1) lets the
  near-ground resolution be dialled INDEPENDENTLY of the physics, run the study
  that WP7 could not. THE QUESTION: does the deep SMF fade tail (p5, p1) CONVERGE
  as the near-ground `Cn2` is resolved with more, thinner screens — or is the
  ~2 dB p5 sensitivity seen in the post-WP7 matched-seed re-test (docs/schmidt-
  crosscheck.md:1273) a real, convergent effect that the default screen count
  under-resolves? WHY IT MATTERS: the fade tail sets the LINK AVAILABILITY
  margin, so an under-resolved tail biases the availability the budget reports;
  the mean is already validated flat, so the tail is the open risk in the
  fidelity-2 space budgets. THE METHOD: hold the grid and the seed set fixed,
  sweep the effective near-ground screen resolution (e.g. raise `min_screens`,
  or add a near-ground refinement to the equal-weight cut), and measure p50 /
  p10 / p5 / p1 of the SMF (point-receiver) fade against screen count at 30 deg
  and a low elevation. Resolve the p5 gap ABOVE the Monte-Carlo noise — the
  re-test showed that needs about 4x the 200 trials (so about 800), which the
  fast `ScreenFactory` now makes cheap. Record the convergence curve and, if the
  tail moves, RE-TIER the default screen count for the tail (this feeds 2-I3,
  the per-channel preset revision, and the receiver-kind floor question). Land it
  as a `validation/` study with a written note in docs/schmidt-crosscheck.md.
- **2-I3. Revise the `QualityPreset` approach (owner-flagged 2026-08-27;
  scope widened 2026-08-29).**
  One preset table serves two channel families that measure differently, and
  the convergence data says they deserve different numbers. The evidence, all
  in the WP7 and post-WP7 notes of docs/schmidt-crosscheck.md: the
  terrestrial case converges by 4 to 5 screens while the slant case needs 7,
  so `min_screens` 15 / 9 / 5 over-serves a horizontal path and the
  terrestrial floor could sit lower; the `rapid` floor of 5 reads 10 percent
  low on an APERTURE index but 22 percent low on a POINT index, so the right
  floor also depends on the receiver (a fibre pays the point figure); and the
  space planner reads the cap per LAYER GROUP while the terrestrial planner
  reads it per EQUAL SLAB, so one `sigma2_r_screen_max` value binds the two
  families differently. The revision: split the preset table per channel
  family (terrestrial / space), and decide whether the floor keys on the
  receiver kind; source every number from the existing sweep data or a new
  sweep, and record it in the tracker. Owner decision on the shape; flagged,
  not built.

  THE GOAL, RESTATED (2026-08-29). The presets are chiefly an INTERNAL
  VALIDATION tool. An end user wants ONE accurate simulation, not a quality
  knob to turn for accuracy. So the deliverable is not only a re-tiered table.
  It is three settings, each sourced from data: (1) a validated MINIMUM per
  channel family, the floor that still hits the reference; (2) a probably-safe
  INTERIM default to ship until the sweep is complete; and (3) the well-sampled
  REFERENCE itself. The knob stays for validation, but the shipped default must
  not ask the user to trade accuracy.

  THE TEST CATALOGUE (2026-08-29). The current floors rest on a NARROW sweep: a
  30 deg slab and a 2 km horizontal path (docs/schmidt-crosscheck.md WP7). That
  is too thin to certify a minimum. Build a bigger catalogue of conditions and
  verify the minimum against ALL of them, especially the EXTREME links: a low
  elevation / high airmass slant, a strong Cn2, a long or turbulent terrestrial
  path, a small and a large aperture, and the three receiver kinds (point,
  aperture, fibre) across uplink and downlink. The fast `ScreenFactory` (about
  10x per screen, validated in validation/waveoptics_speed/) makes the broad
  sweep cheap, so the reason for the narrow one is gone. Source every floor from
  this catalogue and record it in the tracker.
- **2-S1. The Schmidt cross-check gaps S-01 to S-28.** The Schmidt
  foundation layer (`olb/waveoptics/schmidt/`) is validation only, and its
  tracker holds 28 numbered gaps between the book and the production
  wave-optics code. Do NOT copy the rows here; read Table 2 of
  docs/schmidt-crosscheck.md. The HIGH-priority rows are: S-13 (the
  co-moving route has no book equation; compare it against the scaled flat
  grid of Eq. (6.65) — the WP5 example
  `examples/schmidt/propagator_kernels.py` now measures 1.7e-3 soft and
  2.3e-2 hard at m = 247, on a 500 m hop of a 1 mm waist; the tracker still
  asks for the same comparison on the 600 km uplink), S-14 (the split step
  holds ONE flat pitch; Eq. (8.18) gives each step its own, so the grid
  cannot grow with a diverging beam), S-16 (constraints 1 and 2 of Ch. 7 are
  implemented nowhere; `guard = 4.0` and `pixels_per_feature = 16` have no
  source), and S-21 (the three turbulent geometry constraints of Eqs. (9.86)
  to (9.88) are never evaluated). S-22 (the layer moment rule; the same item
  as 2-N1) is now LOW in the tracker, because WP7 measured it. The MEDIUM
  rows S-15 (the absorber shape), S-17 (the
  Fresnel minimum distance), S-20 (the phase pitch rule) and S-27 (the
  aotools subharmonic screen) are recorded there too. Each one is an owner
  decision, because each one moves a production number.
- **2-P1. The temporal (frozen-flow) axis is a stub.** `TemporalScreens`
  raises (olb/waveoptics/turbulence/temporal.py:54); the layer gives
  snapshots only — no fade rate, no fade duration. The design note lives in
  the class docstring.
- **2-P2. The folded / retro double pass is a stub.** `folded_terrestrial`
  and the `"retro"` direction raise (run.py:231, :443, :608). The two
  passes share screens, so they are correlated; that needs its own design.
- **2-P3. No co-moving (spherical) screen.** `split_step` takes a flat grid
  only; a long slant path pays the pixel cost.
- **2-P4. The reciprocity route carries no point-ahead anisoplanatism**
  (the uplink and downlink read the same screens;
  docs/api-waveoptics.md:824).
- **2-I1. `TurbWaveResult` — the rich record is DONE (2026-09-04).** The rule
  was: a minimal scalar record, do NOT extend it piece by piece; the E-field
  inside the receive aperture gets its own design session
  (memory `waveoptics-results-deferred`). That session ran. THE DECISION: the
  per-trial SCALARS stay exactly as they are, and the record gains ONE optional
  field pair — `TurbWaveResult.fields` (the masked receive-plane field of each
  trial, complex64, BEFORE the receive clip) and `TurbWaveResult.patch` (the
  `FieldPatch` that says which grid pixels those values are). The runner stores
  them only when the caller gives `patch_radius_m`, so the old record is
  bit-identical and a budget never reads a field. The store pays for
  `recouple`/`recollect`: a smaller receive aperture, an obscuration, another
  detector, another focal length and another defocus are then a POST-HOC crop,
  with no new propagation. See olb/waveoptics/turbulence/run.py and
  `Campaign` (campaign.py). EARLIER (2026-09-02): `TurbTrial.detector_etas`
  holds the per-arm coupling efficiencies of a multi-detector run, and it stays
  None for a single detector.
- **2-N2. Known numerical readings to keep in view:** the Fourier screen
  structure function reads up to 15 % low over r/r0 0.3–1.6 (ratios only);
  the aperture-averaged analytic factor fails when the aperture holds the
  beam (the 100 mm bucket case); the grid sizer warns past `forvard_max_z`
  and under the `n_max` clamp (olb/waveoptics/grid.py:208, :227); the
  runner warns when the receive aperture reaches the absorbing band
  (turbulence/run.py:476) — an automatic grid-from-aperture size would
  close that one.
- **2-N3. PARTLY DONE (2026-08-29). Speed: tune the grid size and resolution
  along the path**, validated against the well-sampled reference runs (memory
  `waveoptics-speed-exploration`). Speed-plan P2 measured the two obvious ideas
  and BURIED both for the wired scenarios
  (`validation/waveoptics_speed/coarse_screen_experiment.py` and
  `beam_grid_experiment.py`): a coarse-then-interpolate screen loses the
  Fresnel-scale phase and fails the sigma2_I kill line, and a beam-sized screen
  or a per-plane-pitch chain (gap S-14) buys no pixel operations on a space slab
  or a terrestrial path. WHAT REMAINS: a resolution schedule for a long
  CO-MOVING space path, which the runner avoids today because it simulates the
  ~20 km slab alone with a plane-wave input.
- **2-N5. Investigate the auto grid sizer — does it discriminate elevation?
  (owner-flagged 2026-08-29).** In `presentation/gen_data.py` the two hero
  elevations, 90 deg (zenith) and 30 deg, produce the IDENTICAL grid under the
  standard preset: N = 512 pixels and 9 screens for BOTH, although 30 deg is the
  harder path (slant 40 vs 20 km, sigma2_R 0.231 vs 0.065, r0 12.4 vs 18.8 cm).
  The 30 deg side is wider (3.54 vs 2.78 m) and its pixel coarser (6.92 vs
  5.42 mm), but the power-of-two rounding lands both on N = 512, and the screen
  count is pinned at `min_screens = 9` for both because the per-screen Rytov cap
  never binds at these geometries. So the sizer is FLOOR-limited and
  ROUNDING-limited here, not physics-limited, and it gives the harder path no
  finer sampling. QUESTION to settle: is that correct (the zenith case is simply
  over-sampled, so 30 deg needs nothing more), or is the min_screens floor plus
  the next-power-of-two step MASKING a real sampling difference that the 30 deg
  case should pay? Check the achieved `pixels_per_r0` and `fresnel_pixels_min`
  for both against a converged reference. NOTE: the timing puzzle that raised
  this (30 deg ran FASTER than 90 deg in gen_data) is a WARM-UP / run-order
  artefact, NOT sampling — a matched-seed re-run gives round 0 (cold) 90 deg
  1452 ms/trial and 30 deg 1443 ms/trial (equal, as the equal grids predict),
  and only the warm round 1 diverges (90 deg 689, 30 deg 1041 ms) as FFT plans
  and screen caches settle. Per-trial cost tracks the grid, and the grids are
  equal; the presentation numbers reflect which elevation ran first. Pairs with
  2-I3 (the preset revision) and 2-I2 (continuous profiles, which drive the
  screen placement). See `turbulent_grid` in
  olb/waveoptics/turbulence/sampling.py.
- **2-N4a. Run WHOLE fidelity-2 sims in parallel, the NON-TEMPORAL case — DONE
  (2026-09-04).** The snapshot trials are independent, so there is no limit.
  `olb.waveoptics.turbulence.campaign.Campaign` is the answer: it keeps the
  trials on disk in fixed BLOCKS, and `Campaign.run(n, workers=W)` opens ONE
  warm `ProcessPoolExecutor` for the whole call and runs each block SERIALLY
  inside its process. The parallelism lives at ONE level only, because threads
  inside processes over-subscribe the cores. The blocks are bit-identical SLICES
  of one seeded native run, through the runner's new `start_index`, and a
  manifest rebuilds the grid and the plan, so a resumed campaign never re-sizes.
  The fair P3 rerun (2026-09-04) says WHY a pool must stay warm: threads and
  processes TIE on the wall time of ONE run (0.99x space, 1.04x terrestrial),
  and processes win 1.14x to 1.74x in steady state only, because the Windows
  spawn costs 2.5 to 4.4 s. So there is NO automatic selector between the two
  routes, and threads stay the default of one run. See
  docs/waveoptics-efficiency-plan.md Section 8.6 and
  `validation/waveoptics_speed/fair_scaling_rerun.py`.
- **2-N4b. Run WHOLE fidelity-2 sims in parallel, the TEMPORAL case — STILL
  OPEN.** A frozen-flow time axis (2-P1, still a stub) needs the screen arrays
  in a fixed order, so a naive whole-sim parallel split breaks the time
  correlation. Two ways out: (a) build the screen arrays BEFORE the parallel
  fan-out, then hand each worker its ready arrays; or (b) the leaning choice —
  simulate about 1 second of link time for each worker (this holds MANY
  coherence times), and still multiprocess across the 1-second blocks. Decide
  the block length from the coherence time and the wind, and record the choice.
  Pairs with 2-N3 and the P3 scaling data (`validation/waveoptics_speed/`);
  needs the temporal axis (2-P1) first.
- **2-N6. A large-campaign validation is the NEXT STEP (2026-09-04).** The
  campaign store is built and its self-check and demo run at tens of trials
  only. Nobody has run thousands of trials through it yet. Measure: the wall
  time and the disk size of a real campaign, the resume after a kill, the
  behaviour of the fade quantiles as the trial count grows, and the memory of
  the streamed `recouple`/`recollect`. The `EmpiricalSampler` tail rule (ten
  samples past the availability) sets the count: 1,000 trials for 99 percent
  and 10,000 for 99.9 percent. NO numbers yet. NOTE the owner REFUSAL: a
  parametric tail fitted to the simulated bulk, to extrapolate a deeper fade, is
  rejected. The owner does not want extrapolation, and 99.99 percent
  availability is out of scope.

---

## Infrastructure and code debt

- **I-1. Scalar-elevation limits — RESOLVED by a sweep helper (2026-09-02).**
  The gamma-gamma Term (olb/links/downlink.py) and the FAST Term
  (olb/models/fast.py) each model ONE line of sight, so an internal
  vectorisation is not possible: FAST runs one Monte Carlo per geometry, and the
  gamma-gamma Term carries one (alpha, beta) pair. This was CONFIRMED against the
  fast-aosim source, whose own multi-elevation driver
  (`complete_orbit_simulation.py`) builds one `fast.Fast(...)` per zenith angle in
  a loop. So the correct fix is a LOOP, not vectorisation. The new top-level
  helper `olb.budgets_vs_elevation(scenario, elevations, **kwargs)`
  (`olb/sweep.py`) builds a scalar-elevation `CircularOrbit` for each angle from
  `scenario.channel.altitude_m`, calls the family/direction budget function
  (reusing `multidetector._budget_function`), and returns
  `list[(elevation_deg, Budget)]`. It mirrors `multi_detector_budgets` and sits
  ABOVE `olb.links`. A `TerrestrialScenario` has no elevation axis and raises.
  The self-check asserts the I-1 regression: a gamma-gamma array-elevation call
  raises, but the sweep runs it one angle at a time. Docs: docs/api-budget.md.
- **I-2. Duplicate physics copies (Gap 10) — DONE for the production paths
  (lognormal faces 2026-09-02, Rytov std 2026-09-04).** The four crosscheck forms
  now each have ONE canonical home that every budget-facing copy reads:
  1. The lognormal dB FACES (TL-01..04, DL-01..04, crosscheck G-24). The four
     expressions (mean_db, quantile, sampler, the -sigma_l2/2 offset) were
     written inline in BOTH `downlink._lognormal_term` and
     `terrestrial_scintillation_term`. Both now route through the ONE shared
     adapter `olb/models/fade.py` `irradiance_fade_term` with the andrews
     lognormal helpers (`lognormal_params`/`_mean_log`/`_quantile`/`_rvs`), the
     SAME path the gamma-gamma Term already used. The mean and the sampler are
     BYTE identical to the retired inline code; the quantile matches to machine
     precision (~1e-16 dB: the adapter takes -10 log10(exp(x)), the inline code
     took -10 x / ln10 directly). A raw-formula parity guard in each module's
     self-check (and in `fade.py`) enforces this mechanically. The dead
     `scipy.stats.norm` and `_LN10` imports were removed from both link modules.
  2. The plane-wave RYTOV standard deviation (GF-05, KR-23). Two hard-coded
     production copies were left. `plane_wave_scintillation.sigma1_rytov` (missed
     when its neighbour `coherence_radius` was delegated) and the inline
     `olb/links/terrestrial.py` regime gate now both read
     `andrews.scintillation.rytov_variance(wave="plane")`. This MOVES NO NUMBERS:
     the Rytov constant 1.23 is book-identical in the copies and in andrews (the
     0.423-vs-0.4240 constant question is a FRIED matter, not a Rytov one). A
     full-precision before/after harness over 113 values (mean_db, three
     quantiles, seeded 1000-draw samples, every meta field and the assumptions
     frame of the downlink and terrestrial scintillation Terms, plus
     `sigma1_rytov`) shows mean_db, ALL samples, meta and assumptions
     BIT-IDENTICAL; the only movement is the quantile face, worst 3.6e-15 dB (the
     shared-adapter round-trip of point 1). `sigma1_rytov` keeps its `@assumes`
     weak-regime Constraint, no production Term traces it, and the terrestrial
     call sits outside the factory trace, so no assumptions frame moves.
  3+4. The plane-wave COHERENCE RADIUS (GF-10) and single-path FRIED parameter
     (GF-11, KR-25) were ALREADY converged by the Andrews-foundation refactor and
     the de-vendoring: `gaussian_fried.plane_wave_coherence_radius`,
     `gaussian_fried.plane_wave_fried_parameter`, `gaussian_fried.rytov_std`,
     `plane_wave_scintillation.coherence_radius` and `ao.py` all DELEGATE to
     `andrews.structure` (the 0.4240 book chain), and the old
     `my_analysis_modules` copy is gone. So there is NO open production Fried or
     coherence duplicate.
  DELIBERATELY LEFT (not a production duplicate): `andrews/wander.py` keeps its
  own `spherical_fried_parameter` (0.16) and `plane_fried_parameter_slant` (0.42)
  forms. That module is the parked, measurement-only Andrews wander route (0-W8,
  Conflict C-01): no budget consumes it, its self-check validates it
  independently, and its constant convention differs, so it is not a
  trivially-safe converge. Leave it until the wander route is ever wired.
- **I-3. The diverged (Theta, Lambda) coupled-flux feed cross-check test —
  DONE (2026-09-01).** The `uplink_flux.py` module self-check now drives the
  `_scintillation_beam` diverged feed through the coupled-flux on-axis kernel
  (Path A) and compares it against the closed-form Dios beam-wave index
  `beam_wave_scintillation.on_axis_scintillation_index` (Path B), on a weak,
  homogeneous 2 km horizontal-equivalent path. Both read the same launch
  curvature `f0`, so their Theta and Lambda match. The two agree to ~0.002 %
  for the collimated AND the diverged (5x) beam, far inside the 15 % gate. The
  old `ponytail: TODO` in `_scintillation_beam` is retired. (Path A trims the
  single z = L integrand node, a removable 0/0 the coupled-flux kernel makes
  there; production grids never reach it.)
- **I-4. Dead or decorative parameters — DONE (2026-09-01).** The
  `precompensation` half is REFUSED: `SpaceScenario.__post_init__` raises a
  `ValueError` when the field is set on a downlink or a retro link, so a user
  error no longer passes silently; the self-check asserts both refusals. The
  `geometry` half is DOCUMENTED as deliberate: every terrestrial Term keeps the
  unused parameter for the uniform `f(scenario, geometry) -> Term` signature,
  and each docstring says so.
- **I-5. A `TerrestrialScenario` is one direction only — DONE (2026-09-01).**
  `TerrestrialScenario` now holds a `direction` field
  (`TerrestrialDirection`, "forward" | "reverse", default "forward"), and
  `tx_terminal` / `rx_terminal` read it: forward gives tx=near, rx=far (the old
  behaviour, so no caller changes), reverse swaps the two. The channel is
  symmetric, so only the role mapping changes, and the model interface stays
  the same. The terrestrial direction is a DIFFERENT Literal from the space
  `Direction`. The wave-optics layer told the two families apart with
  `hasattr(scenario, "direction")`; that test now reads `"ground"`, because a
  terrestrial scenario also has a direction. Docs updated: CLAUDE.md, docs/architecture.md,
  docs/api-terminal-scenario.md, docs/getting-started.md.

---

## Documentation debt

- **DD-1. DONE (2026-08-28).** The docs no longer call the aperture
  angle-of-arrival function a raising stub. `aperture_arrival_angle_variance`
  delegates to `andrews.structure.angle_of_arrival_variance`; the docs now say
  so and point to 0-W3 (no Term consumes it). Fixed in docs/architecture.md,
  docs/physics.md, and docs/andrews-crosscheck.md (the batch-2 note, the Table 2
  legend, the G-34 row, and the 2.91 constants-ledger row).
- **DD-2. PARTLY DONE (verified 2026-09-04).** The README "Next / planned" graph
  dropped the two closed nodes NT5 (validate the diverged coupled-flux feed;
  measured and closed) and NT8 (thread f0 into the terrestrial Fried call; 0-W2),
  and corrected NT1, which flatly said "DTHETA = 0 today" although the uplink now
  computes the point-ahead offset (1-2). The fidelity-ladder diagram node F2 is
  ALSO corrected: it now reads "SMF and MMF Terms", not "NO Term yet". STILL to
  sweep (owner presentation choice): the README node NP2 frames pre-comp uplink
  scintillation as a "MAJOR GAP" to fill, although Gap 2 is DECIDED (no analytic
  Term; FAST is the model of record).
- **DD-3. DONE.** docs/api-waveoptics.md now carries the `min_screens`
  caveat, the `rmax` factor-4 note and the `fresnel_weight_min` note in the
  `QualityPreset` table.
- **DD-4. DONE (2026-09-02).** The Ch. 3 spectrum-model rows of Crosscheck
  Table 3 that read "not found in olb" are reconciled against
  `olb/turbulence/andrews/spectra.py` (WP3). Four rows now point at the code and
  read IMPLEMENTED, with a "Reconciled 2026-09-02" tag: the Tatarskii inner-scale
  wavenumber (`TATARSKII_KM = 5.92`), the modified-atmospheric inner-scale
  wavenumber (`MODIFIED_KL = 3.3`), the high-wavenumber bump terms (`1.802`,
  `0.254` in `modified_atmospheric`), and the outer-scale conventions
  (`VON_KARMAN_C0 = 2 pi`, `EXPONENTIAL_C0 = 4 pi`, `MODIFIED_EQ23_C0 = 4 pi`).
  A dated header note records the reconciliation. Rows that are still genuinely
  absent from the code are unchanged. NOTE (a follow-up, not DD-4): the
  structure-function / coherence-radius constants (2.914/1.093, 1.64/1.87,
  0.55/0.62) now live in `andrews/structure.py`, not `spectra.py`, so they were
  left for a separate structure.py reconciliation.
- **DD-5. Citation faults — AO-07 addressed (2026-08-28); two left, owner-gated.**
  AO-07: the "Andrews Ch. 3 for a Noll 1976 result" fault is GONE from the code
  (a refactor since 2026-08-26 left the one remaining `ao.py` "Ch. 3" citation on
  the genuine Kolmogorov phase PSD, which is a correct attribution). The Noll
  residual-coefficient citations (`ao.py` module docstring and the constants
  block) were missing the required DOI; added `10.1364/JOSA.66.000207`
  (already used elsewhere). The C-02 reference-plane note is now PAID by the
  assumptions refactor: `gaussian_fried.TX_REFERRED_WEIGHT` and
  `coupled_flux.PATH_WEIGHT` are named Constraints that state the
  transmitter-referred path weight and cite Dios Eq. (3),
  DOI 10.1364/AO.43.003866. STILL OPEN, needs owner physics judgment: the
  aperture-averaged integral in `plane_wave_scintillation.py` still cites
  Ch. 12, Eq. (38) where Ch. 9 Eq. (25) / Ch. 10 Eq. (59) may be closer (PW-05).
- **DD-6. DONE (verified 2026-08-28).** The point-ahead decorrelation framing
  is fixed in BOTH the code and the docs (committed ee23223, 2026-08-24). This
  backlog entry was stale: physics.md 5g (the "NOT a penalty for correcting"
  paragraph), api-budget.md, and the `uplink.py` docstrings all describe the
  per-order decorrelation residual 2 sigma_n^2 (1 - rho_n). The old "correcting
  more modes injects error" wording is present nowhere. Memory
  `pointahead-decorrelation-framing` already reads RESOLVED.
- **DD-7. The wired-versus-available status lives in three places**
  (CLAUDE.md, examples/andrews/README.md, the crosscheck). Each change
  needs three edits; consider one home.
- **DD-8. Crosscheck batch 2 waits on owner input** (the "(owner to
  specify)" columns).
- **DD-9. The measured-validity home EXISTS (2026-09-01).** docs/physics.md
  Section 9, "Measured validity: what the validation scripts certify", collects
  where the fidelity-0 and fidelity-1 models HOLD. It has one entry for each
  physics question (9a to 9i), and each entry gives the model, the reference,
  the measured numbers with a date, the verdict, and the script.
  STANDING RULE: a validation script that reaches a verdict adds an entry there,
  or it updates the entry that it supersedes. A result that lives only in a run
  log, a memory note or a backlog aside is not documented. Entry 9e holds the
  FULL-run numbers of 1-6 (updated 2026-09-01, same day).

---

## External dependencies

- **X-1. RESOLVED for olb (2026-08-28): the coupled-flux kernels are vendored.**
  olb copied them into `olb/turbulence/coupled_flux.py`, cross-validated
  bit-for-bit against the `my_analysis_modules` working tree (which held the
  Dios-verified fixes). olb no longer depends on `my_analysis_modules`. The
  kernel-repo owner may still want to commit its own working tree, but that is
  no longer an olb blocker. (One stale trace is left: the `pyproject.toml`
  Pylance block still names the deleted `olb/_deps.py` and keeps
  `extraPaths = ["../my_analysis_modules"]`. It has no runtime effect.)
- **X-2. KR-24: the kernel keeps three wrong constants**
  (general_atmospherics.py:23 uses 0.54 / 1.22 / 0.509; the book uses
  0.49 / 1.11 / 0.51). The olb half is fixed; the kernel half is open.
- **X-3. Gap 8 and the strong aperture-averaged index need a second
  source** (a literature task, not a code task).

---

## Big reference buckets (do not work these as one item)

- **Crosscheck Table 2: 166 book-capability gap rows, 84 at priority P1**
  (docs/andrews-crosscheck.md, the rows between the Table 2 heading and Table 3).
  Named still-open rows: G-41, G-42, G-71, G-75, G-97, G-98, G-125, G-133,
  G-140, G-151. G-20 (the Gaussian-beam weak-fluctuation gate) is CLOSED in the
  code: `andrews.scintillation.rytov_weak(sigma2_R, Lambda)` applies both
  Ch. 5, Eq. (16) conditions, and the terrestrial and uplink Terms call it. The
  Table 2 row for G-20 is not yet marked. Mine the rest per work package.
- **Crosscheck Table 3: 40 book constants absent from olb** (partly
  stale, see DD-4).
- **Inherent limits (recorded, not fixable from this book):** the
  plane-parallel slant atmosphere (zenith-angle limit 60 deg); the
  two-scale large-scale branch does not reduce to Kolmogorov.
