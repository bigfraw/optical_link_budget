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

1. **Gap 2 is DECIDED and WIRED (2026-08-27): the pre-compensated uplink
   fade comes from FAST.** No trustworthy closed form exists; the model of
   record is the fidelity-1 FAST route, and the wiring is DONE:
   `uplink_fast_term` plus `uplink_budget(fidelity=1)` (the default for a
   pre-compensated scenario). See 0-W1 for the decision record and 1-2 for the
   remaining FAST limits.
2. **HIGH (owner-flagged, 2026-08-27) — stop the reliance on the
   `DEFAULT_HS` 20-layer array.** HV5/7 is a continuous profile; the planner
   and the physics must take a callable, not a hand-discretised grid. See
   2-I2.
3. **DONE — the turbulent screen-count floor `min_screens`.** Work package 7
   resolved it. See 2-N1.
4. **DONE — Gap 3, thread the beam curvature f0 into the Fried call site.**
   Closed 2026-08-27. See 0-W2.
5. **The stale docs that contradict the code.** Cheap, and they mislead every
   later session. See the documentation-debt group.
6. **DONE — olb is self-contained.** The `my_analysis_modules` kernels are
   vendored into olb and `_deps.py` is deleted (2026-08-28). See X-1. The
   kernel-repo commit and the KR-24 constants remain a concern for that repo
   only, not for olb.
7. **The owner decisions.** `downlink_budget` default (0-W5), and the
   FAST-versus-field reference model — whether wave optics ever becomes a
   DEFAULT. The turbulent fidelity-2 Term is now WIRED as opt-in (2-W1, done
   2026-08-28); the reference-model choice stays open. Next owner-requested
   step: an automatic fidelity selector. See 0-W5, 2-W1.

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
  only WARN when the weak-fluctuation limit is exceeded; the downlink
  already routes to gamma-gamma. Route the other links to a strong-regime
  model (a Monte Carlo or the fidelity-2 layer; see the memory
  `strong-fluctuation-numerical`).
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
  FAST Term is trusted at other operating points.
- **1-4. The Dios duplicate (ponytail DEBT).** The analytic
  beam_wave_scintillation path and the coupled-flux MC duplicate the same
  equations; the jitter correction sits in the MC path only
  (olb/turbulence/beam_wave_scintillation.py:46). Converge on ONE
  implementation; do not make a third copy.
- **1-6. Certify the aperture-averaged lognormal power draw against
  fidelity 2 (owner-flagged 2026-08-29).** The question: in WEAK turbulence,
  does the received-power distribution under APERTURE AVERAGING still follow the
  lognormal that the cheap analytic route assumes? The route takes the
  aperture-averaged index sigma2_P = A sigma2_I and draws from a lognormal (see
  olb/links/terrestrial.py:162, olb/links/downlink.py:88). But an
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
  the field does not. The item stays OPEN: the `--full` deep-tail run (1500
  trials, `standard` preset) is still to run, and the aperture-averaging factor A
  is the quantity to look at next.
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
  the WORST relative index error (2.9x) moves the 5 % fade by 0.26 dB only. The
  item stays OPEN on the `--full` deep-tail run; the named target is now the weak
  aperture-averaging factor A over `1 <= D/rho_0 <= 8`, BELOW the beam-filling
  limit.
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
  3. Cache the calibration (the P4 disk cache,
     `olb/waveoptics/turbulence/cache.py`, exists), so a sweep or an optimiser
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
  step. Pairs with 2-P4 (point-ahead anisoplanatism in the reciprocity route) and
  the reference-model gap of 2-W1 (the field reads less coupling loss than FAST —
  part of that gap is this missing correction, so the two are NOT yet comparable).
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
  near-field truncation flag). All in olb/models/coupling/waveoptics.py; see
  `examples/waveoptics/budget_wiring.py`. STILL OPEN and OWNER-GATED: whether
  wave optics ever becomes a DEFAULT. That is the reference-model gap — the field
  reads LESS fibre coupling loss than the incumbents (0.7–2.9 dB less than FAST,
  ~2.5 dB less than the terrestrial analytic Term). FOLLOW-UP
  (owner-requested 2026-08-28): an AUTOMATIC fidelity selector, the way
  `model="auto"` picks a distribution.

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
- **2-W2. The fidelity-2 SMF path ignores `defocus_m` (2026-08-31).** The
  fidelity-2 MMF leg now reads `MMF.defocus_m` (the plane `z = f + defocus_m`, a
  quadratic pupil phase), so a non-focal-plane light bucket is
  simulated. The single-mode leg (`olb/waveoptics/smf.py`, the pupil-mode
  overlap) takes no defocus, so a fidelity-2 SMF budget always reads the
  focal-plane coupling. A fidelity-2 cross-check of the analytic
  `smf_eta_defocused(a, c)` therefore has no field reference yet. The fix is a
  defocus phase on the back-projected fibre mode, the same quadratic factor the
  MMF path uses.
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
  screen placement). olb/waveoptics/turbulence/sampling.py:494.
- **2-N4. Run WHOLE fidelity-2 sims in parallel (owner-flagged 2026-08-29).**
  The current campaign parallelises the trials WITHIN one run (the `Threader`
  thread pool, P3 measured that processes beat threads and threads saturate at 8
  to 16 workers). But the observed CPU AND memory load stays far below the
  machine capacity, so a whole sim (a full `propagate_turbulent_scenario` call)
  can run in parallel with other whole sims to fill the machine. Two cases:
  1. NON-TEMPORAL (the snapshot layer of today): the trials are independent, so
     there is NO limit. Run many whole sims side by side (different scenarios,
     seeds, or blocks), across processes.
  2. TEMPORAL (2-P1, still a stub): a frozen-flow time axis needs the screen
     arrays in a fixed order, so a naive whole-sim parallel split breaks the time
     correlation. Two ways out: (a) build the screen arrays BEFORE the parallel
     fan-out, then hand each worker its ready arrays; or (b) the leaning choice —
     simulate about 1 second of link time per worker (this holds MANY coherence
     times), and still multiprocess across the 1-second blocks. Decide the block
     length from the coherence time and the wind, and record the choice.
  Pairs with 2-N3 and the P3 scaling data (`validation/waveoptics_speed/`);
  needs the temporal axis (2-P1) before case 2 is real.

---

## Infrastructure and code debt

- **I-1. Scalar-elevation limits.** The gamma-gamma Term
  (olb/links/downlink.py:185) and the FAST Term
  (olb/models/fast.py:132) refuse an elevation array. Vectorise or
  loop internally.
- **I-2. Duplicate physics copies (Gap 10).** The Rytov standard deviation,
  the plane-wave coherence radius, and the Fried parameter each exist in
  three places; the lognormal dB faces exist twice (crosscheck TL-01..04,
  GF-05, GF-10, KR-23, KR-25). Converge.
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
- **DD-2. DONE (2026-08-28).** The README "Next / planned" graph dropped the
  two closed nodes NT5 (validate the diverged coupled-flux feed; measured and
  closed) and NT8 (thread f0 into the terrestrial Fried call; 0-W2), and
  corrected NT1, which flatly said "DTHETA = 0 today" although the uplink now
  computes the point-ahead offset (1-2). STILL to sweep (found while fixing this,
  NOT yet done, owner presentation choice): the fidelity-ladder diagram node F2
  reads "NO Term yet", but 2-W1 wired `waveoptics_turbulence_term`; and NP2
  frames pre-comp uplink scintillation as a "MAJOR GAP" to fill, although Gap 2
  is DECIDED (no analytic Term; FAST is the model of record).
- **DD-3. DONE.** docs/api-waveoptics.md now carries the `min_screens`
  caveat, the `rmax` factor-4 note and the `fresnel_weight_min` note in the
  `QualityPreset` table.
- **DD-4. Crosscheck Table 3 is partly stale** — the "not found in olb"
  spectra rows predate WP3's andrews/spectra.py.
- **DD-5. Citation faults — AO-07 addressed (2026-08-28); two left, owner-gated.**
  AO-07: the "Andrews Ch. 3 for a Noll 1976 result" fault is GONE from the code
  (a refactor since 2026-08-26 left the one remaining `ao.py` "Ch. 3" citation on
  the genuine Kolmogorov phase PSD, which is a correct attribution). The Noll
  residual-coefficient citations (`ao.py` module docstring and the constants
  block) were missing the required DOI; added `10.1364/JOSA.66.000207`
  (already used elsewhere). STILL OPEN, need owner physics judgment: the
  aperture-averaged integral cites Ch. 12 where Ch. 9 Eq. (25) / Ch. 10 Eq. (59)
  may be closer (PW-05); and the C-02 reference-plane docstring note is still
  owed.
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
  no longer an olb blocker.
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
