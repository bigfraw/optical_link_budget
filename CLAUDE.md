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
appears as a stochastic turbulence Term plus a deterministic geometric loss, with
only the analytic extinction and pointing Terms alongside. The geometric loss is
the wave-optics vacuum Term for a TERRESTRIAL link, but the ANALYTIC geometric
Term for a SPACE link (the default; a ground-space link is far field, so the
full-path wave vacuum run is skipped — it is slow and grid-noise-limited;
`run_fidelity2(vacuum="wave")` opts back in). Fidelity 1 does not exist for a
terrestrial link (FAST is far-field; a near-field Gaussian beam needs fidelity 2).
Fidelity 2 needs a precomputed `wave` bundle from
`olb.models.waveoptics.run_fidelity2`; the budget never runs the sim itself. See the
README fidelity ladder.

## Architecture (one-way dependency: turbulence <- models and links)

- `olb/terminal.py` — pure data. ALL terminal hardware lives here. A `Terminal`
  holds `aperture_m`, `obscuration_ratio`, `wavelength_m`, `pointing_jitter_rad`,
  an optional `Transmitter` (`waist_m`, `power_dbm`, `m2`, `divergence_rad`), an
  optional `Detector` (`Aperture`, `SMF`, or `MMF`, each with `sensitivity_dbm`;
  `SMF` also carries `focal_length_m`, `mode_field_radius_m`, `optimal_focus`, and
  `defocus_m`;
  `MMF` is a light bucket with `core_radius_m`, `focal_length_m`, an optional
  `numerical_aperture` (the angular acceptance gate; None keeps the old
  spatial-only coupling), `optimal_focus`, and `defocus_m`; `Camera` is a
  focal-plane array with `pixel_pitch_m`, `n_pixels`, `focal_length_m`, and
  `defocus_m`), and a
  `compensation` stack (`TipTilt`, `AO`). `defocus_m` puts the detector at
  z = f + defocus_m; 0.0 is the nominal focal plane. `Detector = Union[Aperture,
  SMF, MMF, Camera]`. A `Camera` is DIAGNOSTIC: no budget builds a coupling Term
  for it, `terrestrial_budget` and `downlink_budget(fidelity=2)` treat it like an
  `Aperture`, and `downlink_budget` at fidelity 0 or 1 raises. All four detector
  dataclasses carry `frac: Optional[float] = None`, the beamsplitter power
  fraction of that arm (see `olb/models/splitter.py`); None on every detector
  keeps the single-detector behaviour. A terminal parameter can only
  be set through a Terminal.
- `olb/scenario.py` — pure data. Two scenario families, one interface. A
  `SpaceScenario` holds two terminals (`ground`, `space`), a `Channel`, the
  `direction` ("uplink" | "downlink" | "retro"), and `availability_target`. A
  `Channel` is the space propagation channel: `site` plus the orbit
  `altitude_m`. A `TerrestrialScenario` is the horizontal (ground-to-ground)
  family. It holds two terminals (`near`, `far`) and a `TerrestrialChannel`
  (`site`, `path_length_m`, `attenuation_db_per_km`, `cn2`). It holds its OWN
  `direction` (`TerrestrialDirection`, "forward" | "reverse"), a DIFFERENT type
  from the space `Direction`, because "terrestrial" is a channel family, not a
  tx/rx geometry. A channel holds NO hardware. Both families expose the SAME interface that the
  models read: `scenario.tx_terminal`, `scenario.rx_terminal`,
  `scenario.channel`. So no model changes between the families. A SpaceScenario
  sets the roles from the direction: uplink -> tx=ground, rx=space; downlink ->
  tx=space, rx=ground; retro -> tx=rx=ground. A TerrestrialScenario sets the
  roles the same way: forward (the default) -> tx=near, rx=far; reverse ->
  tx=far, rx=near. The channel is symmetric, so only the roles change. There is NO `Scenario` alias and NO `Link` dataclass. `Site`
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
  horizontal Beer-Lambert path), `pointing.py`, `splitter.py` (the receive
  beamsplitter: `resolve_fracs` is the frac autosolve — at most ONE detector may
  leave `frac=None` and it takes the remainder, a lone None takes 1.0, the given
  fractions must not add up to more than 1.0 (a sum below 1.0 is the splitter
  excess loss), and a violation raises; `splitter_term(frac)` is the fixed
  `-10*log10(frac)` dB Term of category "system"; `arm_scenario(scenario,
  detector)` is the `dataclasses.replace` of one arm's receive terminal, aware of
  the scenario family and the direction), and the two FIDELITY-named
  modules that sit at the `models/` level (not inside a category package) because
  each spans several Term categories and is named for its fidelity, not one
  physics: `fast.py` (fidelity 1: `smf_fast_term`, the FAST downlink fibre
  coupling, AND `uplink_fast_term`, the pre-compensated uplink turbulence Term;
  the two share the FAST loader, the Cn2 layering, and the AO mapping) and
  `waveoptics.py` (fidelity 2: `run_fidelity2`/`run_waveoptics` (the runners),
  `Fidelity2Bundle`, `waveoptics_vacuum_term` (the deterministic geometric loss),
  `waveoptics_turbulence_term` (the fade), `waveoptics_smf_coupling_term` (the
  turbulent single-mode fibre-coupling face), and `waveoptics_mmf_coupling_term`
  (the turbulent multimode light-bucket coupling face), and
  `waveoptics_vacuum_mmf_term` (the DETERMINISTIC multimode core-capture face of
  a vacuum-only bundle)). The `coupling/` package holds the
  category-native coupling Terms (`_common.py` holds the shared SMF physics: the
  flat-wavefront `smf_eta_max_from_a(a)` AND the defocus-aberrated closed form
  `smf_eta_defocused(a, c) = 2 a^2 |(1-e^-(a^2-ic))/(a^2-ic)|^2` (Shaklan and
  Roddier DOI 10.1364/AO.27.002334; Ruilier and Cassaing
  DOI 10.1364/JOSAA.18.000143);
  `downlink.py` holds `downlink_coupling_term`; `terrestrial.py` holds
  `terrestrial_smf_coupling_term`, `terrestrial_smf_walkoff_term`,
  `terrestrial_mmf_coupling_term`, and the public
  `curvature_focus_shift(scenario)`). The terrestrial Terms ALWAYS charge the
  RECEIVED CURVATURE: a horizontal received beam is a diverging Gaussian of
  phase-front radius R_rx (`olb.beam.phase_front_radius`), so its true focus sits
  at dz_curv = f^2/(R_rx - f) BEYOND the focal plane (S. A. Self, Appl. Opt. 22,
  658 (1983), DOI 10.1364/AO.22.000658), and the Terms evaluate the detector at
  dz_eff = defocus_m - dz_curv. `optimal_focus` is a focal-LENGTH rule and NEVER
  moves the detector; set `detector.defocus_m = curvature_focus_shift(scenario)`
  for a tracked (aligned) coupler. The package RE-EXPORTS the coupling-category Terms
  that a fidelity module owns (`smf_fast_term` from `fast.py`,
  `waveoptics_smf_coupling_term` and `waveoptics_mmf_coupling_term` from
  `waveoptics.py`), so a coupling Term is
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
  that drops every turbulence quantity but keeps the static and jitter parts.
  Every terrestrial coupling Term charges the received-curvature defocus, in BOTH
  branches, because that curvature is static optics, not turbulence. At
  `fidelity=2` an MMF receiver gets one more Term, the wave-optics light-bucket
  core coupling, which already holds the detector defocus).
  `bidirectional.py` (`defocused_terminal`, `bidirectional_terrestrial`,
  `BidirectionalBudget`; a monostatic collimator has ONE defocus dz that drives
  BOTH the transmit divergence and the receive coupling, so the wrapper returns
  the forward and the reverse budget of one horizontal path. Two fidelity-0
  limits: dz > 0 (a converging launch) is outside the divergence model, and a
  diverged monostatic terminal pays |dz| + dz_curv of receive defocus).
- `olb/waveoptics/` — the fidelity-2 field propagation layer. The CORE carries no
  turbulence. A trimmed port of LightPipes (BSD-3-Clause, see `LIGHTPIPES_LICENSE.txt` in the
  package) that keeps the LightPipes names and call order: `field.py` (Field,
  Begin, Normal, Power, Intensity, Phase, SubIntensity), `sources.py` (GaussBeam,
  PlaneWave, CircAperture, CircScreen), `propagators.py` (Forvard, Fresnel,
  GForvard; the three take a FLAT grid only, and each one raises on a spherical
  field), and `lenses.py` (Lens, LensForvard, LensFresnel, Convert; the thin lens
  and the spherical (co-moving) coordinate route, which moves the grid with the
  beam so a long space link stays sampled on a small pixel count). Five
  olb-native modules sit on that core: `smf.py` (the fibre mode
  and the overlap coupling efficiency; it takes NO defocus, see backlog 2-W2),
  `mmf.py` (the multimode light-bucket
  coupling: `focal_intensity` and `mmf_coupling_efficiency`, both of which take a
  `defocus_m` (the plane z = f + defocus_m, a quadratic pupil phase of SIGN
  `exp(-i*pi*defocus_m*rho^2/(lam*f^2))`, so a DIVERGING received beam couples
  best at a POSITIVE defocus_m)), `camera.py` (the focal-plane array:
  `camera_image` bins the focused spot onto the square camera pixels, and
  `spot_metrics` -> `SpotMetrics` gives the centroid, the second-moment radius
  and the on-sensor power fraction; it reuses `mmf.focal_intensity`, so it is a
  DIAGNOSTIC layer that builds no Term), `grid.py`
  (`GridSpec.for_scenario`, the
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
  `olb/waveoptics/turbulence/`: `screens.py` (the phase screens: `screen_r0`,
  `phase_screen`, `Screen`, and TWO generators — the DEFAULT `ScreenFactory`, a
  fast self-contained generator (cached sqrt-PSD filter, separable
  outer-product subharmonics, two screens per FFT; numpy and scipy only), and
  the opt-in `aotools` wrapper as the reference path, a lazy LGPL import),
  `splitstep.py` (`super_gaussian_boundary`, `split_step`), `sampling.py` (the
  turbulent grid sizer and screen planner: `QualityPreset`, `PRESETS`
  reference/standard/rapid, `ScreenPlan`, `SamplingReport`, `turbulent_grid`.
  The DEFAULT space planner is CONTINUOUS (item 2-I2 step 1): it INTEGRATES a
  Cn2 callable `cn2(h)` (the site HV5/7 when none is given) and cuts the slab
  into equal-Rytov-weight screens at Cn2-weighted centroids. An explicit
  `hs`/`cn2_profile` array takes the LEGACY discrete planner; `DEFAULT_HS` is
  the fallback for that array caller ONLY, no longer the physics grid of the
  default budget. `cn2`/`h_top_m` thread through `run_waveoptics`,
  `run_fidelity2`, the runners, and the cache),
  `run.py` (`TurbTrial`, `TurbWaveResult`, `propagate_turbulent_scenario`,
  `propagate_turbulent_field` (one snapshot as a complex receive-plane Field,
  for a plot; it does NOT extend the scalar record); both take
  `screen_generator="olb"` (the default) | "aotools"; the two draw DIFFERENT
  atmospheres for the same seed, and the statistics agree; the
  `folded_terrestrial` stub), `cache.py` (`cached_propagate_turbulent_scenario`,
  an opt-in, off-by-default disk cache of whole runs, extendable by block; no
  budget calls it), and `temporal.py` (the `TemporalScreens`
  NotImplementedError stub). It gives SNAPSHOTS: one atmosphere per seed, no time
  axis. The trials are independent, so `propagate_turbulent_scenario` takes an
  optional `threader` (`olb.waveoptics.Threader`, a general thread pool in
  `threader.py`, default `min(16, cores)` workers) that runs them across threads;
  the FFT releases the GIL, so the threads give a real speed-up. A space scenario
  always propagates the DOWNLINK slab, and an uplink reads the same field through
  the Shapiro reciprocity overlap (DOI 10.1364/JOSA.61.000492). Neither part
  builds a Term, and neither changes a budget.
  `propagate_turbulent_scenario` also takes an optional `detectors` sequence (the
  beamsplitter arms): each trial then computes the coupling efficiency of EVERY
  arm from the SAME clipped receive field and reports them in the ONE new tuple
  field `TurbTrial.detector_etas`, so N arms cost ONE Monte Carlo. The default
  `detectors=None` leaves that field None, so the single-detector record is
  bit-identical. The `frac` of a detector is IGNORED there.
- `olb/multidetector.py` — the top-level per-arm budget helper,
  `multi_detector_budgets(scenario, geometry, detectors, wave=None, **kwargs)`.
  It gives one `(detector, Budget)` pair for each arm: it copies the scenario
  with `arm_scenario`, it calls the budget function of the scenario family and
  direction, and it adds the fixed `splitter_term` (an arm of fraction 1.0 gets
  no splitter row). The input scenario does not change, and a per-arm error is
  NOT caught. The split is CROSS-CUTTING, not the physics of one link, so this
  module sits ABOVE `links/` in the one-way dependency order: it imports
  `olb.links`, `olb.models` and `olb.terminal`, and nothing in `links/` or
  `models/` imports it back. `olb/__init__.py` exports it.
- `olb/sweep.py` — the top-level elevation-sweep helper,
  `budgets_vs_elevation(scenario, elevations, geometry_factory=None, **kwargs)`.
  It gives one `(elevation_deg, Budget)` pair for each angle: it builds a
  scalar-elevation `CircularOrbit` from `scenario.channel.altitude_m` (or a
  supplied `geometry_factory`), and it calls the family/direction budget function
  (it reuses `multidetector._budget_function`). It exists because some Terms model
  ONE line of sight and refuse an elevation ARRAY (FAST runs one Monte Carlo per
  geometry; the gamma-gamma Term carries one `(alpha, beta)` pair), so the correct
  fix is a LOOP, not vectorisation (backlog I-1). A `TerrestrialScenario` has no
  elevation axis and raises. Like `multidetector.py` it is CROSS-CUTTING and sits
  ABOVE `links/`. `olb/__init__.py` exports it.
- `olb/results.py` — `Term` (three faces: mean_db, quantile, sampler) and
  `Budget`. Monte Carlo is not a separate path. The Budget asks each Term for
  samples, not means.
- `olb/assumptions.py` — each Term declares its beam type, turbulence regime,
  and spectrum. `Budget.check()` flags a scenario that breaks an assumption. A
  physics function OWNS its assumptions through an `@assumes(...)` decorator that
  attaches a `FuncAssumptions` record and optional `Constraint` runtime checks; a
  Term factory opens `with trace_assumptions() as trace:` around its physics
  calls, and every decorated function that runs registers its source and any
  violation, so the Term inherits the union through `trace.merge(...)`.
  `merge_assumptions(*records)` recomposes finished Terms (retro).
  `Assumptions` gained `constraints` (`(source, Constraint)` pairs) and
  `provenance` (traced source names); `flag(reason, source=)` tags a
  scenario-level fact. `results.py` gained `provenance` and `n_constraints`
  columns on `assumptions_frame()`, a new `constraints_frame()`, and a
  `Budget.check()` untraced-Term guard (a turbulence or coupling Term with empty
  provenance is reported; a legitimately untraced Term self-declares
  `provenance=["untraced: ..."]`). Scope: `olb/turbulence/**` plus the link and
  model factories; `olb/waveoptics/` is deferred. See `docs/architecture.md`
  Section 5 and `docs/api-budget.md`.
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

The function-owned assumptions refactor is MERGED (2026-08-29). Every public
physics function in `olb/turbulence/**` carries an `@assumes(...)` decorator (88
decorated functions across 18 modules; the `olb/assumptions.py` self-check asserts
a floor of 85, so a dropped decorator fails mechanically). A prose-only limit is
now a runtime `Constraint` check, and a Term factory opens `trace_assumptions()`
so the Term inherits the union. Newly ENFORCED checks flip `ok` to not-ok in cases
that read ok before, and this is the INTENDED effect of the refactor, NOT a
regression. Two flip a CURRENT budget: the Gaussian second weak condition (a
focused beam) and the extended-Marechal limit (a strong AO residual). The first
zenith enforcement is real but LATENT: the production factories trace the parallel
`plane_wave_scintillation` and `uplink_flux` feeders, not `andrews.paths`, so the
zenith check fires only when a factory wires the `andrews.paths` slant integrators.
Honest status items:
- The 0.25 house rule keeps ONE canonical definition, `LOGNORMAL_PDF_LIMIT = 0.25`
  in `andrews/scintillation.py`. The old `WEAK_FLUCTUATION_LIMIT` name is fully
  retired: no source file references it.
- The terrestrial SMF walk-off weak-limit gap is CLOSED by a FACTORY regime flag,
  not by a function-owned check, because the vendored Dios wander kernel
  `coupled_flux.beam_wander_variance` has no runtime check to inherit. A
  function-owned weak-regime check in that kernel is an OPEN follow-up.
- `MARECHAL_SIGMA2_MAX = 1.0` is DUPLICATED in `olb/turbulence/ao.py` and
  `olb/links/uplink.py` (the one-way `turbulence <- links` dependency forbids
  importing up). A centralise-down is an OPEN follow-up.
- The terrestrial scintillation hard-flag MIGRATED from the `sigma_R^2 >= 1.0`
  axis to the beam-wave index axis (`sigma_I^2 >= 2.4`, owned by
  `beam_wave_scintillation.on_axis_scintillation_index`). The two coincide on the
  tested triggers; a narrow band can now read ok where it flagged, and the tighter
  lognormal-PDF flag (0.25) backstops the common cases.
- `olb/waveoptics/` is DEFERRED (a different, numerical-sampling assumption
  family). Its Terms self-declare `provenance=["untraced: ..."]` to pass the
  `Budget.check()` guard.

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

- **Several detectors, the master turbulence switch, and the Camera are BUILT
  (2026-09-02).** A `Terminal` still holds ONE detector: about twenty detector
  dispatch sites read that one field, so a receive path that feeds more than one
  detector makes one budget for each ARM. Each detector carries `frac`, the
  beamsplitter power fraction, and `olb.models.splitter.resolve_fracs` autosolves
  it (at most ONE `frac=None` takes the remainder; a lone None takes 1.0; a sum
  above 1.0 raises). `olb.multidetector.multi_detector_budgets` gives one
  `(detector, Budget)` pair for each arm, and it adds the fixed `splitter_term`.
  At fidelity 2 the arms share ONE Monte Carlo:
  `run_fidelity2(..., detectors=[...])` returns one `Fidelity2Bundle` for each
  arm from one run, through the new `TurbTrial.detector_etas`. This is EXACT: a
  beamsplitter scales the field of an arm by a constant, and every coupling
  efficiency is power-normalised, so `frac` never touches eta and it enters one
  time as the fixed dB Term. The fidelity-2 MASTER TURBULENCE SWITCH is also
  wired: `run_fidelity2(turbulence=False)` makes no screens and no trials and it
  gives a vacuum-only bundle (`turbulent=None`; the EMPTY bundle of a space link
  is valid), and all three fidelity-2 budgets honour `turbulence=False` and keep
  the deterministic Terms (a terrestrial or downlink MMF receiver keeps one
  deterministic core-capture Term, `waveoptics_vacuum_mmf_term`). This mirrors
  the fidelity-0 master switch, so the ladder reads the same at every rung. The
  new `Camera` detector and `olb/waveoptics/camera.py` (`camera_image`,
  `spot_metrics`, `SpotMetrics`) are DIAGNOSTIC only: they measure the spot shape
  and the spot position for a tracking loop, and they build NO Term. The
  power-to-pixel-brightness (a holistic camera model) is DEFERRED to backlog
  2-W3.
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
- **TL-05 / C-05 code half DONE (2026-08-29).** The weak-fluctuation gate is
  now one shared, beam-aware helper: `olb.turbulence.andrews.scintillation`
  `rytov_weak(sigma2_R, Lambda=None)` returns `"weak"|"soft"|"hard"` on the
  Rytov-variance axis with tiers `RYTOV_CONFIDENT_WEAK=0.3` (soft) and
  `WEAK_REGIME_LIMIT=1.0` (hard). For a Gaussian beam it reads the receiver-plane
  `Lambda` and applies BOTH Ch. 5, Eq. (16) conditions (the binding strength is
  `sigma2_R * max(1, Lambda**(5/6))`), so a FOCUSED beam trips a gate a
  plane-wave test would pass (TL-05). The terrestrial scintillation Term
  (`olb/links/terrestrial.py`) and the uplink coupled-flux path
  (`olb/turbulence/uplink_flux.py`, `olb/links/uplink.py`) both call it; the
  uplink uses the Dios reliability edge `UPLINK_SIGMA2X_LIMIT=0.6` on the
  log-amplitude variance (hard_limit = 4*0.6 on the sigma2_R axis), MORE generous
  than the book because the two-scale coupled-flux index saturates gracefully.
  The lognormal-PDF house rule is now a DISTINCT named limit
  `LOGNORMAL_PDF_LIMIT=0.25` on sigma2_I (fade-PDF shape), no longer conflated
  with the regime gate. The old `WEAK_FLUCTUATION_LIMIT` name is RETIRED from the
  code: the four constants live only in `andrews/scintillation.py`. The downlink
  (`links/downlink.py`) now uses `LOGNORMAL_PDF_LIMIT` for the lognormal Term and
  the model="auto" switch (a PDF decision) and `WEAK_REGIME_LIMIT` for the
  gamma-gamma Term's regime flag (was the factor-of-4 error: it tested the true
  sigma2_R against 0.25). `fast.py` (both the SMF and the pre-compensated uplink
  Terms) gives the amplitude log-normal two flags: a REGIME hard-flag at 1.0 and
  a lognormal-PDF caution at 0.25 (meta `amplitude_rytov_regime`). `waveoptics.py`
  uses `WEAK_REGIME_LIMIT` for its weak/strong regime LABEL (the split-step solver
  is valid at all strengths, so 0.25 there was also too tight).
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
  slab with a plane-wave input), so the slab outputs are vacuum-limit-1.0
  penalties and the geometric loss is a SEPARATE additive Term. **The SPACE
  geometric loss is ANALYTIC by default (2026-08-31, `run_fidelity2` default
  `vacuum="analytic"`).** A ground-space link is far field, so the analytic
  geometric Term (`geometric_loss_term` + the opt-in `tx_gaussian_efficiency_term`
  truncation) is exact AND cheap, and the budget uses it (`wave.vacuum` is None).
  The wave-optics vacuum run over the full slant range is SKIPPED: it costs ~14 s
  and is grid-noise-limited (it cannot resolve the mm-scale aperture edges over a
  ~2000 km path, so the loss scatters +/- 1 to 4 dB and does not converge at a
  practical grid size — this ALSO removes the coarse-grid `GridSpec` warning that
  full-path run emitted). `run_fidelity2(vacuum="wave")` opts a space link back
  into the wave-optics vacuum Term (research / cross-check). A TERRESTRIAL link
  keeps the wave vacuum (its penalty is turbulent / vacuum on the SAME grid, an
  exact baseline; `vacuum="analytic"` raises for terrestrial). The analytic
  geometric Term matches the fidelity-0/1 geometric exactly (uplink 33.63 dB),
  and `validation/vacuum_loss/vacuum_loss_validation.py` cross-checks it against a
  well-resolved wave solve (terrestrial far field agrees to ~0.15 dB) and shows
  the space full-path scatter. All default budgets are UNCHANGED (terrestrial
  fidelity=0, downlink/uplink fidelity=1). Fidelity 1 is
  UNAVAILABLE for terrestrial (raises, backlog 1-1); fidelity 0 is unavailable
  for an uncorrected uplink (raises); fidelity 2 is unavailable for a
  pre-compensated uplink and for retro (raises — the folded double pass shares
  screens). The turbulence Term carries a SNAPSHOT-ONLY flag (fade depth, not
  rate/duration) and an under-sampled-tail quantile warning
  (`olb.results.EmpiricalSampler`). `examples/waveoptics/budget_wiring.py`
  demonstrates all three. STILL owner-gated: whether wave optics ever becomes a
  DEFAULT (the 2-W1 fibre-coupling reference gap stays open — the field reads 1
  to 3 dB LESS coupling loss than FAST/analytic. The TERRESTRIAL MMF half of that
  gap is now QUANTIFIED and mostly explained: with the received curvature charged
  (2026-08-31) it falls from about 7 dB to about 1.2 dB, and the residual is the
  Airy-versus-Gaussian spot shape. The SPACE half is untested against that
  correction, so the gap is NOT closed). OWNER FOLLOW-UP (2026-08-28):
  an AUTOMATIC fidelity selector, the way `model="auto"` picks a distribution.
  The turbulent layer is SNAPSHOT-only (`temporal.py` is a NotImplementedError
  stub). Its DEFAULT screen generator is self-contained (numpy and scipy only);
  `aotools` is now the opt-in reference generator only (LGPL-3.0, the optional
  `screens` extra). Deliberately deferred: the results record is minimal scalars
  (do NOT extend `TurbWaveResult` piece by piece), the temporal frozen-flow axis,
  a co-moving (spherical) screen, and the folded/retro double pass (correlated
  screens). `examples/waveoptics/` demonstrates the layer with seven scripts
  (three vacuum, three turbulent, and the budget-wiring demo).
- **The coupled-flux kernels are VENDORED (2026-08-28).** olb copied them into
  `olb/turbulence/coupled_flux.py`, cross-validated bit-for-bit against the
  `my_analysis_modules` working tree (which held the Dios-verified fixes). So
  olb now carries the fixed version regardless of whether the kernel repo ever
  commits them, and olb no longer depends on `my_analysis_modules` at all.
- **The fidelity-2 speed campaign is DONE (2026-08-29; P0 to P4, see
  `docs/waveoptics-efficiency-plan.md` Section 8 and `validation/waveoptics_speed/`).**
  P0 found screen generation was ~80% of a trial. P1 added the fast
  `ScreenFactory`. TWO owner-decided DEFAULTS then changed (commit e8c7f77):
  `screen_generator` defaults to `"olb"` (the fast generator; `"aotools"` is the
  opt-in reference), and `Threader` caps at `min(16, cores)`. The olb default
  changes the random draws of a seeded fidelity-2 run, so it changes the
  fidelity-2 budget numbers that read `run_fidelity2`; pass
  `screen_generator="aotools"` to reproduce an old aotools run bit-identically.
  The broad validity pass (`generator_validation.py`) shows the olb generator
  agrees with aotools and the analytic index across geometries, presets, the
  outer scale, and the FADE TAIL. P2 measured and BURIED two grid ideas (coarse
  screens, beam-following grid). P3 measured the parallel scaling (processes beat
  threads; threads saturate at 8 to 16 workers). P4 added an opt-in, off-by-
  default disk cache (`olb/waveoptics/turbulence/cache.py`), extendable by block.
  OWNER FOLLOW-UP: a true single-seed tail extension needs a start-index argument
  in the runner (the cache uses block sub-seeds, so a cached run is not the
  bit-identical trials of a native run).
- **The non-focal-plane (defocus) sensing and the received curvature are WIRED
  (2026-08-31, see `validation/defocus/`).** `SMF`/`MMF` carry `defocus_m`; the
  terrestrial coupling Terms grow the spot over `dz_eff = defocus_m - dz_curv`
  and displace it with the ray-optics chief-ray tilt lever `(f+dz)*theta`, which
  keeps the PHYSICAL dz.
  The fidelity-2 `defocus_m` SIGN was inverted and is FIXED. The SMF MEAN defocus
  penalty is now MODELLED (`smf_eta_defocused`), so only the SMF walk-off
  DISPLACEMENT response stays geometric (a loud flag, backlog 0-P15). OPEN: the
  fidelity-2 SMF leg reads no defocus (backlog 2-W2); a converging monostatic
  launch is outside the bidirectional model (backlog 0-P16); the deterministic
  (non-jitter) pointing offset is still not modelled.
- **`examples/andrews/`** demonstrates the layer script by script; its
  README repeats this wired-versus-available status.
