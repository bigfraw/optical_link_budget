# CLAUDE.md — optical_link_budget (olb)

Guidance for Claude Code that works in this repository.

## Purpose

The package builds optical (laser) ground-to-space link budgets with
atmospheric propagation, fade statistics, and Monte Carlo. It models uplink,
downlink, and retroreflected links to a LEO satellite.

## Architecture (one-way dependency: turbulence <- models and links)

- `olb/terminal.py` — pure data. ALL terminal hardware lives here. A `Terminal`
  holds `aperture_m`, `obscuration_ratio`, `wavelength_m`, `pointing_jitter_rad`,
  an optional `Transmitter` (`waist_m`, `power_dbm`, `m2`, `divergence_rad`), an
  optional `Detector` (`Aperture`, `SMF`, or `MMF`, each with `sensitivity_dbm`;
  `SMF` also carries `focal_length_m`, `mode_field_radius_m`, and `optimal_focus`;
  `MMF` is a light bucket with `core_radius_m`, `focal_length_m`, an optional
  `numerical_aperture` (the angular acceptance gate; None keeps the old
  spatial-only coupling), and `optimal_focus`), and a
  `compensation` stack (`TipTilt`, `AO`). A terminal parameter can only be set
  through a Terminal.
- `olb/scenario.py` — pure data. Two scenario families, one interface. A
  `SpaceScenario` holds two terminals (`ground`, `space`), a `Channel`, the
  `direction` ("uplink" | "downlink" | "retro"), and `availability_target`. A
  `Channel` is the space propagation channel: `site` plus the orbit
  `altitude_m`. A `TerrestrialScenario` is the horizontal (ground-to-ground)
  family. It holds two terminals (`near`, `far`) and a `TerrestrialChannel`
  (`site`, `path_length_m`, `attenuation_db_per_km`, `cn2`). It has NO
  `direction`, because "terrestrial" is a channel family, not a tx/rx geometry.
  A channel holds NO hardware. Both families expose the SAME interface that the
  models read: `scenario.tx_terminal`, `scenario.rx_terminal`,
  `scenario.channel`. So no model changes between the families. A SpaceScenario
  sets the roles from the direction: uplink -> tx=ground, rx=space; downlink ->
  tx=space, rx=ground; retro -> tx=rx=ground. A TerrestrialScenario is one-way:
  tx=near, rx=far. There is NO `Scenario` alias and NO `Link` dataclass. `Site`
  stays. A SpaceScenario also holds an optional uplink `precompensation` source
  (`DownlinkBeacon`, `LaserGuideStar` (a placeholder), or None); the uplink
  budget reads it to select the turbulence physics.
- `olb/turbulence/` — pure physics. It imports only numpy, scipy, and `_deps`.
  It does not import a scenario or Term. Files: `profiles.py` (Cn2 profiles,
  `default_cn2_profile`, `DEFAULT_HS`), `plane_wave_scintillation.py`
  (plane-wave scintillation indices, aperture-averaging integral; the
  space-to-ground downlink model), `beam_wave_scintillation.py` (Dios
  Gaussian-beam scintillation, on and off axis; the uplink model),
  `anisoplanatism.py` (Stone 1994 angular
  anisoplanatic phase variance, with the finite adaptive-optics band and
  `max_radial_order`), `uplink_flux.py` (the LEO-uplink coupled-flux Monte
  Carlo wrapper),
  `angle_of_arrival.py` (the received tip-tilt of a Gaussian beam: the
  beam-wander arrival tilt is the working model; the aperture angle-of-arrival
  tilt now delegates to `andrews/structure.py`),
  and `andrews/` — the Andrews and Phillips foundation layer, nine modules of
  pure book physics (`aperture.py`, `beam.py`, `distributions.py`, `paths.py`,
  `scintillation.py`, `spectra.py`, `structure.py`, `temporal.py`,
  `wander.py`). Each function cites its chapter, equation number and printed
  page from DOI 10.1117/3.626196. The files above KEEP their names and their
  signatures and call it. Put new book physics there, not in a link module.
  `olb/models/fade.py` turns one irradiance model into the three Term faces.
- `olb/models/` — Term factories `f(scenario, geometry) -> Term`. Each factory
  is named for the physics it computes. Some use a link-specific simplification,
  and the name says so: `geometric.py`, `extinction.py` (`slant_extinction_term`
  for the slant airmass path AND `terrestrial_extinction_term` for the
  horizontal Beer-Lambert path), `pointing.py`, and the `coupling/` package
  (`_common.py` holds the shared SMF physics; `downlink.py` holds
  `downlink_coupling_term`; `terrestrial.py` holds `terrestrial_smf_coupling_term`,
  `terrestrial_smf_walkoff_term`, and `terrestrial_mmf_coupling_term`; `fast.py`
  holds the FAST fibre coupling. `from olb.models.coupling import <name>` still
  works).
- `olb/links/` — per-link Terms and budget assembly: `uplink.py`
  (`uplink_turbulence_term`, `uplink_point_ahead_term`, `uplink_fitting_term`,
  `uplink_budget`; the budget dispatches on the scenario `precompensation`
  source, so a DownlinkBeacon + AO replaces the coupled-flux Term with the AO
  error budget = fitting error (Noll) + point-ahead anisoplanatism (Stone). BUT
  that corrected budget is PHASE-ONLY: it drops the scintillation that the
  coupled-flux Term carried, so it understates the deep fade. This is a MAJOR
  known gap. The Terms flag it (`NO SCINTILLATION`). Fix it before the corrected
  uplink fade is trusted), `downlink.py`
  (`downlink_scintillation_term`, `downlink_budget`), `retro_space.py`
  (`retro_space_budget`; retroreflection as a retransmission, SPACE only).
  `retro_budget` is a backward-compatible alias of `retro_space_budget`, kept in
  `olb/links/__init__.py` (there is no `retro.py` file). A short terrestrial
  retro link needs its own module.
  `terrestrial.py` (`terrestrial_budget`; horizontal ground-to-ground link;
  the geometric, horizontal-extinction, and pointing Terms are exact.
  `terrestrial_scintillation_term` gives a real lognormal fade with three faces.
  It uses the Dios on-axis Gaussian-beam scintillation index and the weak
  aperture-averaging factor. `terrestrial_budget` turns it on by default for an
  aperture or no-detector receiver. An `SMF` detector takes the mean-only
  fibre-coupling Term, plus the receive tip-tilt walk-off fade Term
  (`terrestrial_smf_walkoff_term`) when the coupling optics are set. An `MMF`
  (light bucket) takes the spot-in-core coupling Term plus the same walk-off fade
  (`terrestrial_mmf_coupling_term`). The walk-off reads the received tip-tilt from
  `olb.turbulence.angle_of_arrival` (beam wander) plus the receive jitter; the
  coupling Term keeps the higher-order residual only, so the tip-tilt is not
  counted two times. `terrestrial_budget` also takes a master `turbulence` switch
  that drops every turbulence quantity but keeps the static and jitter parts).
- `olb/waveoptics/` — the fidelity-2 field propagation layer. The CORE carries no
  turbulence. A trimmed port of LightPipes (BSD-3-Clause, see `LIGHTPIPES_LICENSE.txt` in the
  package) that keeps the LightPipes names and call order: `field.py` (Field,
  Begin, Normal, Power, Intensity, Phase, SubIntensity), `sources.py` (GaussBeam,
  PlaneWave, CircAperture, CircScreen), `propagators.py` (Forvard, Fresnel,
  GForvard; the three take a FLAT grid only, and each one raises on a spherical
  field), and `lenses.py` (Lens, LensForvard, LensFresnel, Convert; the thin lens
  and the spherical (co-moving) coordinate route, which moves the grid with the
  beam so a long space link stays sampled on a small pixel count). Three
  olb-native modules sit on that core: `smf.py` (the fibre mode
  and the overlap coupling efficiency), `grid.py` (`GridSpec.for_scenario`, the
  automatic grid sizer with a manual override, `beam_magnification`, and
  `forvard_max_z`), and `run.py`
  (`propagate_scenario` -> `WaveResult`, one end-to-end propagation). The sizer
  selects the ROUTE and the runner obeys it: `for_scenario` tries a flat grid
  first and falls back to the scaled (co-moving) grid, `GridSpec.scaled` records
  the choice, and `propagate_scenario` runs GForvard (an almost untouched
  Gaussian), the flat Fresnel convolution, or the three-call lens recipe
  (Lens -> LensFresnel -> Convert) on a scaled grid. The core is
  pure numpy and scipy. It imports nothing from the rest of olb, so the turbulent
  split-step layer uses the same propagators. That layer is the sub-package
  `olb/waveoptics/turbulence/`: `screens.py` (aotools-backed phase screens, a
  lazy LGPL import, `screen_r0`, `phase_screen`, `Screen`), `splitstep.py`
  (`super_gaussian_boundary`, `split_step`), `sampling.py` (the turbulent grid
  sizer and screen planner: `QualityPreset`, `PRESETS`
  reference/standard/rapid, `ScreenPlan`, `SamplingReport`, `turbulent_grid`),
  `run.py` (`TurbTrial`, `TurbWaveResult`, `propagate_turbulent_scenario`,
  `propagate_turbulent_field` (one snapshot as a complex receive-plane Field,
  for a plot; it does NOT extend the scalar record), the
  `folded_terrestrial` stub), and `temporal.py` (the `TemporalScreens`
  NotImplementedError stub). It gives SNAPSHOTS: one atmosphere per seed, no time
  axis. The trials are independent, so `propagate_turbulent_scenario` takes an
  optional `threader` (`olb.waveoptics.Threader`, a general thread pool in
  `threader.py`) that runs them across threads; the FFT releases the GIL, so the
  threads give a real speed-up. A space scenario always propagates the DOWNLINK
  slab, and an uplink reads the same field through the Shapiro reciprocity
  overlap (DOI 10.1364/JOSA.61.000492). Neither part builds a Term, and neither
  changes a budget.
- `olb/results.py` — `Term` (three faces: mean_db, quantile, sampler) and
  `Budget`. Monte Carlo is not a separate path. The Budget asks each Term for
  samples, not means.
- `olb/assumptions.py` — each Term declares its beam type, turbulence regime,
  and spectrum. `Budget.check()` flags a scenario that breaks an assumption.
- `olb/_deps.py` — the ONLY module that imports the shared physics kernels from
  `my_analysis_modules`. Set the `MY_ANALYSIS_MODULES` environment variable, or
  place that repo at `D:\repos\my_analysis_modules`. The `fast` package is
  optional (HV57 Cn2); without it, use `default_cn2_profile`.

## Conventions

- All documentation uses ASD-STE100 Simplified Technical English. See
  `CONVENTIONS.md`. This applies to docstrings, comments, and commit messages.
- Loss is positive dB. Gain is negative dB.
- EVERY equation needs a DOI. Each formula in the code or documentation must
  cite the source paper or book by DOI (in the docstring or a comment next to
  the equation). No uncited physics.
- Run a module with `python -m olb.<...>` from the repository root. The package
  uses package-relative imports.
- Each module has an `if __name__ == '__main__':` self-check.

## Working preferences

- Delegate substantial code writing to Opus 5 subagents, guided by the ponytail
  skill (the laziest solution that works; borrow the shared kernels, do not
  duplicate them; no speculative abstraction). Keep the thin interface
  consistent across the models.
- Every subagent prompt that writes code or documentation must include the
  ASD-STE100 rule.

## Current state

The Andrews foundation layer EXISTS. `olb/turbulence/andrews/` holds nine
modules of pure book physics: `aperture.py`, `beam.py`, `distributions.py`,
`paths.py`, `scintillation.py`, `spectra.py`, `structure.py`, `temporal.py` and
`wander.py`. Every equation cites its chapter, equation number and printed page
from Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196. The older
turbulence modules keep their names and their signatures, and their bodies call
the new layer. `olb/models/fade.py` turns one irradiance model into the three
Term faces. See `docs/physics.md` Section 5h and `docs/andrews-crosscheck.md`.

The downlink budget now selects its distribution:
`downlink_scintillation_term(..., model="auto")` gives the lognormal Term below
sigma2_I = 0.25 and the gamma-gamma Term at or above it. The gamma-gamma Term is
valid at every fluctuation strength, but it models a POINT receiver, because the
book gives no aperture-averaged downlink index in that regime.

The Schmidt foundation layer EXISTS on the branch `schmidt`.
`olb/waveoptics/schmidt/` holds four modules of pure book physics:
`fourier.py` (Chs. 2 and 3), `fresnel.py` (Chs. 6 and 8), `sampling.py`
(Chs. 7 and 8) and `turbulence.py` (Ch. 9). Every equation cites its chapter,
equation number and printed page from Schmidt (2010), DOI 10.1117/3.866274.
Andrews owns the ANALYTIC value of a quantity; Schmidt owns the SIMULATION
rule. The layer is VALIDATION ONLY: no budget, no Term, no sizer and no runner
reads it, the sub-package exports nothing, and the LightPipes production code
keeps its bodies. Three example scripts in `examples/schmidt/` measure the
production layer against the book. See `docs/physics.md` Section 8,
`docs/api-waveoptics.md` Section 10, `examples/schmidt/README.md`, and the
tracker `docs/schmidt-crosscheck.md` (the chapter index, the glossary, the
42-row forward map, the 28 gaps S-01 to S-28, and the constants ledger). The
production modules now carry the book equation numbers in their docstrings, and
`olb/waveoptics/grid.py` `forvard_max_z` is CORRECTLY cited: it is constraint 4,
Ch. 7, Eq. (7.59), printed p. 127, at m = 1, not "Ch. 6".

Open items:

- **The turbulent screen-count floor `min_screens` is UNJUSTIFIED, and the
  merge fallback couples the count to the profile sampling. TO REVISE.** In
  `olb/waveoptics/turbulence/sampling.py` a `QualityPreset` carries two screen
  rules. `sigma2_r_screen_max` (0.05 / 0.10 / 0.25) is the Rytov thin-screen
  UPPER bound; it cites Schmidt Ch. 9 and it sets the count on strong paths.
  `min_screens` (15 / 9 / 5) is a LOWER floor for weak paths, where the Rytov
  cap passes with one screen but one screen gives phase only, no scintillation.
  BUT `min_screens` has NO derivation and NO DOI (the one field in the preset
  table without one), and it was never exercised: `_merge_layers` does not clamp
  UP to `min_screens` when the natural merge undershoots it; it BAILS OUT and
  returns one screen per Cn2 layer. So every weak space case gets 20 screens —
  the layer count of `DEFAULT_HS`, not physics. A finely sampled real profile
  would explode the count (a 200-layer profile gives 200 screens for a slab that
  the Rytov cap says needs one). THE Ch. 9 EVIDENCE IS NOW IN, and it settles
  two of the three parts. (a) Schmidt gives NO screen-count floor. Eq. (9.90),
  printed p. 174, is a sampling floor of the FFT method; the "5 to 10" of
  Sec. 9.2.5, printed p. 165, counts the unknowns of Eq. (9.75); and the 11
  planes of Sec. 9.5.2, printed p. 177, come with no formula. So 15 / 9 / 5
  cannot be cited to Schmidt, and neither can any other integer. (b) The
  principled replacement is the layer MOMENT rule, Ch. 9, Eq. (9.65), printed
  p. 164: the layered Cn2 must match the continuous profile for the moments
  0 <= m <= 7. It fixes the positions and the strengths together, it decouples
  the count from the profile sampling, and it gives a real floor of 4 (8
  equations, 2 free numbers per screen). (c) The `rmax = 0.1` cap of Listing
  9.5, printed p. 175, is on the LOG-AMPLITUDE variance `sigma_chi^2`, and
  `sigma2_r_screen_max` is on the plane-wave Rytov variance, which is 4 times
  it (measured 3.9994). So the book cap is 0.4 on the olb number, and the
  presets are 8x / 4x / 1.6x stricter. THE REVISION ITSELF IS STILL OPEN (work
  package 7): make `_merge_layers` place and weight the screens to minimise
  `moment_error`, take the count from the larger of `min_planes` and the
  `rmax` cap with a floor of 4, and report the achieved moment error in
  `SamplingReport`. Do NOT keep 15 / 9 / 5 with a Schmidt citation. The three
  turbulent examples measured their agreement numbers at 20 screens, so re-run
  them after the change. The evidence is the WP7 gate verdict in
  `docs/schmidt-crosscheck.md`.
- **Gap 2, the pre-compensated uplink, is STILL open.**
  `andrews.paths.uplink_scintillation_index(tracked=True)` gives the floor of the
  residual scintillation (Ch. 12, Eqs. (57) to (60)), but NO budget reads it yet.
  The beacon-plus-adaptive-optics uplink budget stays phase-only, and its Terms
  still say `NO SCINTILLATION`. Wiring that Term is the next physics task.
- **Gap 3 is closed at the physics layer only.** `andrews.beam.beam_params`
  takes any input curvature f0, and `andrews.structure.coherence_radius` takes
  the beam. But the single-path `gaussian_fried.gaussian_fried_parameter` keeps
  its collimated signature, and the terrestrial fibre-coupling call site in
  `olb/models/coupling/terrestrial.py` still passes no curvature. Thread f0 into
  that call to make a deliberately diverged beam drive the Fried parameter.
- **Gap 8, the annular (obscured) receive aperture, needs another source.** A
  full-text search of the book finds no obscured-aperture filter.
- **TL-05**: the terrestrial weak gate tests one plane-wave threshold on a
  Gaussian beam. Ch. 5, Eq. (16), printed p. 140, needs two.
- **`downlink_budget` still defaults to `model="lognormal"`.** The selector
  `model="auto"` exists but is opt-in. The switch is an owner decision,
  because the gamma-gamma Term is point-receiver (see above) and the change
  moves the strong-regime total by several dB.
- **Built but NOT consumed by any budget yet** (each is a deliberate,
  owner-gated wiring step, because each changes budget numbers):
  `andrews/temporal.py` (Greenwood, tau0, fade rate and duration — no Term
  reads them); the inner/outer-scale branches (no Term passes `l0`/`L0`);
  the Andrews Ch. 6 wander route in `andrews/wander.py` (the uplink budget
  keeps the Dios/Belmonte kernel route, per Conflict C-01); the K
  distribution; the fidelity-2 `olb/waveoptics/` layer, the vacuum core AND the
  turbulent `olb/waveoptics/turbulence/` split step (see below).
- **The fidelity-2 wave-optics layer is BUILT and SELF-CHECKED, but NO budget
  consumes it.** `olb/waveoptics/` propagates a real complex field, and it agrees
  with the fidelity-0 analytic Terms in the far field with a light truncation, to
  0.02 dB. It is the no-turbulence validator that answers the near-field flag in
  `olb/models/gaussian_efficiency.py`: a hard-truncated beam inside the Rayleigh
  range breaks the far-field truncation efficiency. The layer now covers the
  SPACE link too, through the co-moving grid of `lenses.py`: a 600 km uplink with
  a hard truncation and a central obscuration runs on a 2048-pixel grid, and its
  total agrees with the fidelity-0 total to 0.011 dB. Compare the TOTAL only. The
  per-Term numbers do NOT compare, because `tx_efficiency_loss_db` is an on-axis
  far-field gain ratio and `tx_truncation_db` is a power ratio. There is no
  fidelity-2 Term,
  because such a Term moves the totals of an existing budget. That wiring is an
  owner decision. The TURBULENT split-step part of fidelity 2 is now BUILT and
  SELF-CHECKED too, at `olb/waveoptics/turbulence/`. It is SNAPSHOT-only (one
  atmosphere per seed; `temporal.py` is a NotImplementedError stub). It imports
  `aotools` for the screens (LGPL-3.0, the optional `screens` extra; olb imports
  it, olb does not copy it). A space scenario ALWAYS propagates the downlink slab
  and an uplink reads it through the Shapiro reciprocity overlap. Three
  comparison examples run it: the uplink reciprocity mean agrees with the
  coupled-flux MC to 0.18 dB at the zenith and 0.54 dB at 30 deg; the downlink
  aperture scintillation agrees with the analytic index at 30/60/90 deg (ratios
  1.10 to 1.27); and the fibre-coupling comparisons read LESS loss than the
  analytic and FAST models (the terrestrial SMF Term is about 2.5 dB higher than
  the field, and FAST is 0.7 to 3 dB higher across elevation). STILL: no budget
  consumes it, no Term exists, and the wiring is an owner decision. Deliberately
  deferred: the results record is minimal scalars (a results-processing design
  chat is planned — do NOT extend `TurbWaveResult` piece by piece), the temporal
  frozen-flow axis, a co-moving (spherical) screen, and the folded/retro double
  pass (the two passes share screens, so they are correlated).
  `examples/waveoptics/` demonstrates the layer with six scripts, three vacuum
  and three turbulent.
- **The kernel repo has uncommitted fixes.** `coupled_flux.py` in
  `D:\repos\my_analysis_modules` is untracked there; the Dios-verified
  fixes sit in its working tree only. The owner must commit them.
- **`examples/andrews/`** demonstrates the layer script by script; its
  README repeats this wired-versus-available status.
