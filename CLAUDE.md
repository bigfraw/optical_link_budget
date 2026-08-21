# CLAUDE.md — optical_link_budget (olb)

Guidance for Claude Code that works in this repository.

## Purpose

The package builds optical (laser) ground-to-space link budgets with
atmospheric propagation, fade statistics, and Monte Carlo. It models uplink,
downlink, and retroreflected links to a LEO satellite.

## Architecture (one-way dependency: turbulence <- models and links)

- `olb/turbulence/` — pure physics. It imports only numpy, scipy, and `_deps`.
  It does not import Scenario or Term. Files: `profiles.py` (Cn2 profiles,
  `default_cn2_profile`, `DEFAULT_HS`), `scintillation.py` (scintillation
  indices, aperture-averaging integral), `coupled_flux.py` (coupled-flux Monte
  Carlo wrapper).
- `olb/models/` — direction-agnostic Term factories `f(scenario, geometry) ->
  Term`: `geometric.py`, `transmittance.py`, `pointing.py`.
- `olb/links/` — per-direction Terms and budget assembly: `uplink.py`
  (`uplink_turbulence_term`, `uplink_budget`), `downlink.py`
  (`downlink_scintillation_term`, `downlink_budget`), `retro.py`
  (`retro_budget`).
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
