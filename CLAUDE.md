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

Open items:

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
- **Conflict C-01** needs Belmonte 2000 to close the 3.50 beam-wander gap.
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
- **The kernel repo has uncommitted fixes.** `coupled_flux.py` in
  `D:\repos\my_analysis_modules` is untracked there; the Dios-verified
  fixes sit in its working tree only. The owner must commit them.
- **`examples/andrews/`** demonstrates the layer script by script; its
  README repeats this wired-versus-available status.
