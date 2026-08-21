# optical_link_budget (olb)

Optical (laser) ground-to-space link budgets with atmospheric propagation, fade
statistics, and Monte Carlo. The package models uplink, downlink, and
retroreflected links at optical wavelengths (for example 1550 nm) to a LEO
satellite.

## Dependency

This package reuses proven physics kernels from a sibling repository,
`my_analysis_modules`. `olb/_deps.py` is the only module that imports them.

- Set the environment variable `MY_ANALYSIS_MODULES` to the location of
  `my_analysis_modules`, or place it at `D:\repos\my_analysis_modules` (the
  default path).
- `fast-aosim` is optional (`pip install fast-aosim`). It supplies the
  Hufnagel-Valley HV57 Cn2 profile, and the fidelity-1 single-mode-fibre modal
  coupling (`smf_fidelity="fast"`). Without it, pass an explicit `cn2_profile`,
  or use the built-in `default_cn2_profile` (which uses `get_c2n`), and use the
  default `smf_fidelity="reciprocity"`.

## Install

```
pip install -e .
```

## Quickstart

```python
import numpy as np
from olb import Scenario, Link, Site, CircularOrbit, uplink_budget, retro_budget

# Uplink budget over an analytic circular orbit.
scenario = Scenario(
    link=Link(tx_power_dbm=40.0, rx_sensitivity_dbm=-40.0, pointing_jitter_rad=2e-6),
    altitude_m=600e3,
)
budget = uplink_budget(scenario, CircularOrbit(600e3, 60.0))
print(budget.to_frame())                 # itemised terms
print(budget.assumptions_frame())        # model constraints per term
print(budget.check())                    # flags a scenario that breaks an assumption

mc = budget.monte_carlo(20000, rng=np.random.default_rng(0), availabilities=(0.99,))
print(mc["margin_db"][0.99])             # 99 % link margin [dB]
```

## Structure

The package uses one-way dependencies: `turbulence/` <- `models/` and `links/`.

- `olb/turbulence/` — pure physics (no Scenario, no Term): Cn2 profiles,
  scintillation indices, aperture averaging, coupled-flux Monte Carlo.
- `olb/models/` — direction-agnostic Term factories: geometric spreading,
  atmospheric extinction, pointing jitter.
- `olb/links/` — per-direction Terms and budget assembly: `uplink_budget`,
  `downlink_budget`, `retro_budget`.
- `olb/results.py` — `Term` (mean / analytic quantile / sampler) and `Budget`.
  Monte Carlo is not a separate path. The Budget asks each Term for samples.
- `olb/assumptions.py` — the model constraints (beam type, turbulence regime,
  spectrum) that each Term declares, and the check that flags a broken
  assumption.

## Documentation convention

All documentation uses ASD-STE100 Simplified Technical English. See
[CONVENTIONS.md](CONVENTIONS.md).
