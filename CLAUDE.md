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
  optional `Detector` (`Aperture` or `SMF`, each with `sensitivity_dbm`), and a
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
  stays.
- `olb/turbulence/` — pure physics. It imports only numpy, scipy, and `_deps`.
  It does not import a scenario or Term. Files: `profiles.py` (Cn2 profiles,
  `default_cn2_profile`, `DEFAULT_HS`), `scintillation.py` (scintillation
  indices, aperture-averaging integral), `coupled_flux.py` (coupled-flux Monte
  Carlo wrapper).
- `olb/models/` — direction-agnostic Term factories `f(scenario, geometry) ->
  Term`: `geometric.py`, `transmittance.py` (slant airmass AND horizontal
  Beer-Lambert), `pointing.py`.
- `olb/links/` — per-link Terms and budget assembly: `uplink.py`
  (`uplink_turbulence_term`, `uplink_budget`), `downlink.py`
  (`downlink_scintillation_term`, `downlink_budget`), `retro_space.py`
  (`retro_space_budget`; retroreflection as a retransmission, SPACE only).
  `retro.py` is a backward-compatible alias that re-exports `retro_budget =
  retro_space_budget`. A short terrestrial retro link needs its own module.
  `terrestrial.py` (`terrestrial_budget`; horizontal ground-to-ground link;
  the geometric, horizontal-extinction, and pointing Terms are exact, but
  `terrestrial_scintillation_term` is a RESERVED SLOT that raises
  NotImplementedError until the Andrews Gaussian-beam forms are added).
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

Add deliberate uplink beam divergence, the way `tn2_kepler` did it. Recast a
diverged beam of aperture radius `w0` and far-field half-angle divergence
`theta` as a Gaussian from a virtual waist behind the aperture:

    w_v = lambda / (pi * theta)
    d   = zR(w_v) * sqrt((w0/w_v)^2 - 1)      # virtual distance behind aperture
    free-space radius at range z = gaussz(w_v, d + z)

`theta` must be at least the diffraction limit `lambda/(pi*w0)`. Reference:
`D:/misc code/tn2_kepler/fso_spot_size.py` (`virtual_waist`,
`free_space_radius`, `spot_sizes`).

Plan: add `divergence_rad` (half-angle) to `Link`; port `virtual_waist` and
`free_space_radius` into `olb/turbulence/beam.py`; wire the diverged beam size
into the geometric, pointing, and turbulence uplink terms. The shared
`coupled_flux_montecarlo` has no divergence argument, but the lower-level shared
functions accept a `w_free` override. Recommended: reimplement the short uplink
Monte Carlo loop in olb by composing those lower-level functions with
`free_space_radius`, without editing `my_analysis_modules`.
