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
  tilt is a DEFERRED stub that raises `NotImplementedError`, owner to specify).
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

## Next task (ASAP)

Generalise the Gaussian-beam Fried parameter past the collimated case. The
single-path form `gaussian_fried_parameter` in
`olb/turbulence/gaussian_fried.py` fixes the input curvature `Theta0 = 1`
through the constant `COLLIMATED_THETA0`. So it holds only for a collimated
beam. The profile form `gaussian_fried_parameter_profile` already accepts a
phase-front radius of curvature `f0` and computes `Theta0 = 1 - L/f0`. But its
default `f0 = inf` keeps the beam collimated, and the one call site (the
terrestrial single-mode-fibre coupling in `olb/models/coupling/terrestrial.py`) does not
pass `f0`.

Plan: add a curvature argument to `gaussian_fried_parameter`. Thread it through
`output_beam_params`, the way the profile form already does. A diverged beam has
`f0 < 0`, so `Theta0 > 1`. A focused beam has `0 < f0`, so `Theta0 < 1`. Then
feed the diverged-beam curvature into the coupling call. So a deliberately
diverged uplink beam also drives the Fried parameter. The package already
recasts that beam through a virtual waist in `olb/beam.py` for the geometric and
the scintillation Terms. See the flag in `docs/physics.md` section 5e.
