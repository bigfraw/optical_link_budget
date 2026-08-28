# CLAUDE.md — optical_link_budget (olb)

Guidance for Claude Code that works in this repository.

## Purpose

The package builds optical (laser) ground-to-space link budgets with
atmospheric propagation, fade statistics, and Monte Carlo. It models uplink,
downlink, and retroreflected links to a LEO satellite.

The organising idea is a **fidelity ladder**, chosen with one `fidelity=0|1|2`
argument per budget. **Fidelity 0** is analytic (closed-form, the most
assumptions). **Fidelity 1** is statistical (a distribution or Monte Carlo — FAST
coupling, coupled-flux uplink — a real fade, some assumptions). **Fidelity 2** is
wave optics (a split-step field solve, assumption-free, the most expensive); it
appears as TWO Terms, a deterministic vacuum-optics Term (the full no-turbulence
loss launch to detector) and a stochastic turbulence Term, with only the analytic
extinction and pointing Terms alongside. Fidelity 1 does not exist for a
terrestrial link (FAST is far-field; a near-field Gaussian beam needs fidelity 2).
Fidelity 2 needs a precomputed `wave` bundle from
`olb.models.waveoptics.run_fidelity2`; the budget never runs the sim itself. See the
README fidelity ladder.

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
- `olb/turbulence/` — pure physics. It imports only numpy, scipy, and the olb
  leaf modules (`units`, `beam`). It does not import a scenario or Term. Files:
  `profiles.py` (Cn2 profiles, the Hufnagel-Valley `get_c2n`, the Bufton wind
  `v_wind`, `default_cn2_profile`, `DEFAULT_HS`), `coupled_flux.py` (the vendored
  Dios coupled-flux kernels for the uplink MC), `plane_wave_scintillation.py`
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
  horizontal Beer-Lambert path), `pointing.py`, and the two FIDELITY-named
  modules that sit at the `models/` level (not inside a category package) because
  each spans several Term categories and is named for its fidelity, not one
  physics: `fast.py` (fidelity 1: `smf_fast_term`, the FAST downlink fibre
  coupling, AND `uplink_fast_term`, the pre-compensated uplink turbulence Term;
  the two share the FAST loader, the Cn2 layering, and the AO mapping) and
  `waveoptics.py` (fidelity 2: `run_fidelity2`/`run_waveoptics` (the runners),
  `Fidelity2Bundle`, `waveoptics_vacuum_term` (the deterministic geometric loss),
  `waveoptics_turbulence_term` (the fade), and `waveoptics_smf_coupling_term` (the
  turbulent fibre-coupling face)). The `coupling/` package holds the
  category-native coupling Terms (`_common.py` holds the shared SMF physics;
  `downlink.py` holds `downlink_coupling_term`; `terrestrial.py` holds
  `terrestrial_smf_coupling_term`, `terrestrial_smf_walkoff_term`, and
  `terrestrial_mmf_coupling_term`), and it RE-EXPORTS the coupling-category Terms
  that a fidelity module owns (`smf_fast_term` from `fast.py`,
  `waveoptics_smf_coupling_term` from `waveoptics.py`), so a coupling Term is
  discoverable in the coupling namespace whatever its fidelity. `from
  olb.models.coupling import <name>` still works for every coupling Term.
- `olb/links/` — per-link Terms and budget assembly. Every budget takes one
  whole-path `fidelity=0|1|2` argument (the fidelity ladder; see `## Purpose`).
  `uplink.py` (`uplink_turbulence_term`, `uplink_point_ahead_term`,
  `uplink_fitting_term`, `uplink_budget`; the budget dispatches on the scenario
  `precompensation` source crossed with `fidelity`. A DownlinkBeacon + AO
  pre-compensates: `fidelity=1` (the default, the model of record) is ONE FAST
  Monte-Carlo Term (`uplink_fast_term`) with the point-ahead decorrelation and a
  real fade; `fidelity=0` is the analytic pair = fitting error (Noll) +
  point-ahead anisoplanatism (Stone), PHASE-ONLY and MEAN-ONLY: no scintillation
  and no fade; `fidelity=2` raises (reciprocity has no AO/point-ahead). An
  UNCORRECTED uplink: `fidelity=1` (default) = coupled flux; `fidelity=0` raises
  (no analytic mean-only); `fidelity=2` = the two wave-optics Terms. That
  mean-only limit is a DECISION (2026-08-27): no trustworthy analytic form exists
  for the scintillation of a pre-compensated beam. The analytic Terms carry loud
  flags (`NO SCINTILLATION, NO FADE`, plus the extended-Marechal limit flag). See
  backlog 0-W1 and 1-2), `downlink.py`
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
- olb is SELF-CONTAINED: it no longer depends on `my_analysis_modules`. The
  physics kernels once borrowed through `olb/_deps.py` are now VENDORED into olb:
  the unit conversions in `olb/units.py`, the Gaussian-beam `gaussz`/`zR` in
  `olb/beam.py`, the `Satellite`/`SatellitePass` geometry in `olb/geometry.py`,
  the Hufnagel-Valley Cn2 and Bufton wind in `olb/turbulence/profiles.py`, and
  the Dios coupled-flux kernels in `olb/turbulence/coupled_flux.py`. `_deps.py`
  is deleted. The `fast` package (FAST fibre coupling / HV57 Cn2) stays an
  optional external dependency, imported lazily; without it, use
  `default_cn2_profile`.

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

- **The turbulent screen-count floor `min_screens` is RESOLVED (work package
  7).** In `olb/waveoptics/turbulence/sampling.py`, `_merge_layers` now clamps a
  weak path UP to EXACTLY `min_screens` contiguous Cn2-weighted groups, through
  the new `_equal_weight_groups`. The old bail-out, which returned one screen
  per Cn2 layer, is gone. So the screen count follows the PRESET and not the
  layer count: a 20-layer `DEFAULT_HS` profile and a 200-layer profile of the
  same atmosphere both give `min_screens` screens on a weak slab. The Rytov cap
  `sigma2_r_screen_max` still RAISES the count above the floor on a strong path,
  unchanged. A profile that has fewer layers than `min_screens` warns, because
  the planner does not split a layer, and it keeps its layers. The integers
  15 / 9 / 5 are CONFIRMED, and their source is an olb convergence sweep, not
  Schmidt: the book gives no screen-count floor. The sweep holds the grid fixed
  and it moves the count only; the aperture scintillation index of a 30 deg
  downlink slab is 19 percent low at 3 screens, 10 percent low at 5, and flat
  from 7 up, and the mean collected power holds inside 0.11 dB everywhere. So 9
  and 15 sit on the plateau, and 5 is the stated rapid compromise. No preset may
  go under 4, the moment floor of Ch. 9, Eq. (9.65), printed p. 164 (8 equations
  against 2 free numbers per screen). The grouping does not SOLVE Eq. (9.65),
  but the Cn2-weighted centroid matches all 8 moments of the default profile to
  better than 1 percent; the module self-check measures that against the
  `olb.waveoptics.schmidt` reference layer. See WP7 in
  `docs/schmidt-crosscheck.md` for the full sweep table.
- **Gap 2 is DECIDED (2026-08-27): the pre-compensated uplink gets NO analytic
  scintillation Term.** `andrews.paths.uplink_scintillation_index(tracked=True)`
  is OPTIMISTIC there, not a floor: it models a perfect tilt removal, the
  correction decorrelates over the point-ahead angle mode by mode, and a
  decorrelated correction reshapes the beam, so the Ch. 12 normalisation
  breaks. The beacon-plus-adaptive-optics budget stays phase-only and
  mean-only, and its Terms carry loud flags (`NO SCINTILLATION, NO FADE`, plus
  the extended-Marechal flag past sigma2 = 1 rad^2, T. S. Ross,
  DOI 10.1364/AO.48.001812). The model of record is the fidelity-1 FAST Monte
  Carlo with the point-ahead DTHETA, and its uplink entry point EXISTS
  (2026-08-27): `uplink_fast_term` in `olb/models/fast.py`, consumed
  by `uplink_budget(fidelity=1)` (the default for a pre-compensated scenario).
  The remaining FAST limits are backlog 1-2.
- **Gap 3 is WIRED (2026-08-27).** The terrestrial fibre-coupling call site in
  `olb/models/coupling/terrestrial.py` now reads the launch curvature f0 from
  the transmitter divergence through `olb.beam.launch_curvature` and passes it
  to `gaussian_fried_parameter_profile`, so a deliberately diverged beam drives
  its own r0. `launch_curvature` is one shared implementation (the Dios feed in
  `olb/turbulence/uplink_flux.py` calls it too). The single-path
  `gaussian_fried.gaussian_fried_parameter` keeps its collimated signature (a
  tidy-up; the budgets use the profile form, which is general in f0).
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
  distribution.
- **The fidelity-2 wave-optics layer is WIRED into the budgets as
  `fidelity=2` (2026-08-28, branch `fidelity2-budget-wiring`).** A fidelity-2
  budget shows TWO Terms, both from `olb/waveoptics/`: a DETERMINISTIC
  vacuum-optics Term (`waveoptics_vacuum_term`, the full no-turbulence loss
  launch to detector from `propagate_scenario` — truncation + geometric spread +
  aperture capture + vacuum fibre coupling) and a STOCHASTIC turbulence Term
  (`waveoptics_turbulence_term`, from the split-step Monte Carlo). Only the
  analytic extinction (absorption) and pointing (mechanical jitter) Terms stay;
  the analytic geometric, launch-truncation, scintillation, and coupling Terms
  DROP. The caller precomputes both records ONCE with
  `olb.models.waveoptics.run_fidelity2` -> `Fidelity2Bundle`, and passes `wave=`;
  the budget never runs the sim. TERRESTRIAL is simulated end to end on one flat
  grid (the vacuum run shares that grid, so the turbulence penalty = turbulent /
  vacuum is exact); SPACE cannot be (the turbulent runner does only the ~20 km
  slab with a plane-wave input), so the vacuum run uses its own co-moving grid
  over the full slant range and the two Terms add (the slab outputs are
  vacuum-limit-1.0 penalties). The vacuum-optics Term matches the analytic
  geometric closely (uplink 33.67 vs 33.63 dB). All default budgets are UNCHANGED
  (terrestrial fidelity=0, downlink/uplink fidelity=1). Fidelity 1 is
  UNAVAILABLE for terrestrial (raises, backlog 1-1); fidelity 0 is unavailable
  for an uncorrected uplink (raises); fidelity 2 is unavailable for a
  pre-compensated uplink and for retro (raises — the folded double pass shares
  screens). The turbulence Term carries a SNAPSHOT-ONLY flag (fade depth, not
  rate/duration) and an under-sampled-tail quantile warning
  (`olb.results.EmpiricalSampler`). `examples/waveoptics/budget_wiring.py`
  demonstrates all three. STILL owner-gated: whether wave optics ever becomes a
  DEFAULT (the 2-W1 fibre-coupling reference gap stays open — the field reads 1
  to 3 dB LESS coupling loss than FAST/analytic). OWNER FOLLOW-UP (2026-08-28):
  an AUTOMATIC fidelity selector, the way `model="auto"` picks a distribution.
  The turbulent layer is SNAPSHOT-only (`temporal.py` is a NotImplementedError
  stub). It imports `aotools` for the screens (LGPL-3.0, the optional `screens`
  extra). Deliberately deferred: the results record is minimal scalars (do NOT
  extend `TurbWaveResult` piece by piece), the temporal frozen-flow axis, a
  co-moving (spherical) screen, and the folded/retro double pass (correlated
  screens). `examples/waveoptics/` demonstrates the layer with seven scripts
  (three vacuum, three turbulent, and the budget-wiring demo).
- **The coupled-flux kernels are VENDORED (2026-08-28).** olb copied them into
  `olb/turbulence/coupled_flux.py`, cross-validated bit-for-bit against the
  `my_analysis_modules` working tree (which held the Dios-verified fixes). So
  olb now carries the fixed version regardless of whether the kernel repo ever
  commits them, and olb no longer depends on `my_analysis_modules` at all.
- **`examples/andrews/`** demonstrates the layer script by script; its
  README repeats this wired-versus-available status.
