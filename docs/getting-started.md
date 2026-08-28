# Getting started with olb

`olb` (optical_link_budget) builds optical laser link budgets. It models
ground-to-space uplink, downlink, and retroreflected links, and horizontal
ground-to-ground links. Each budget gives a mean loss, an analytic fade margin,
or a Monte Carlo of the joint distribution.

## 1. Install and the dependency

Install the package from the repository root:

```
pip install -e .
```

`olb` reuses proven physics kernels from a sibling repository,
`my_analysis_modules`. The module `olb/_deps.py` is the only module that imports
them. Give `olb` the location of that repository in one of two ways:

- Set the environment variable `MY_ANALYSIS_MODULES` to its path.
- Or place it at the default path `D:\repos\my_analysis_modules`.

`fast-aosim` is an optional extra (`pip install fast-aosim`). It adds two things:

- The Hufnagel-Valley HV57 Cn2 profile (turbulence strength against height).
- The fidelity-1 single-mode-fibre coupling. Set `downlink_budget(fidelity=1)`.
  This is the only statistical SMF model. It gives a mean, a quantile, and a
  fade.

Without `fast-aosim`, pass an explicit `cn2_profile`, or use the built-in
`default_cn2_profile`. For fibre coupling, use `fidelity=0`. The analytic model
gives the mean coupling loss only. It gives no fade.

`aotools` is a second optional extra (`pip install aotools`, or
`pip install -e .[screens]`). It draws the random phase screens of the fidelity-2
turbulent split-step layer, `olb.waveoptics.turbulence`. `olb` imports it. `olb`
does not copy it, because `aotools` is LGPL-3.0. No budget needs it: the layer
builds no Term.

## 2. The mental model

All hardware lives on a `Terminal`. A `Terminal` owns a telescope aperture, an
operating wavelength, a tracking jitter, and up to three optional parts:

- a `Transmitter` (it launches a beam; it carries the waist and the launch
  power),
- a `Detector` (an `Aperture` bucket, or an `SMF` single-mode fibre; it carries
  the sensitivity),
- a `compensation` stack (receive-side wavefront correction: `TipTilt`, `AO`).

A scenario pairs two terminals with a channel. There are two scenario families:

- A `SpaceScenario` holds a `ground` terminal, a `space` terminal, a `Channel`
  (a `Site` plus an orbit `altitude_m`), and a `direction`. The direction is
  `"uplink"`, `"downlink"`, or `"retro"`. It resolves which terminal transmits
  and which receives.
- A `TerrestrialScenario` holds a `near` terminal, a `far` terminal, and a
  `TerrestrialChannel` (a `Site`, a `path_length_m`, an extinction, and a Cn2).
  The link is one-way: `near` transmits, `far` receives.

Both families expose the same interface: `scenario.tx_terminal`,
`scenario.rx_terminal`, and `scenario.channel`. So the models read one interface
and do not change between the families. A channel holds no hardware.

The geometry is separate from the scenario. Use `CircularOrbit(altitude_m,
elevation_deg)` for a space link, or `HorizontalPath(path_length_m)` for a
terrestrial link. A budget function takes a scenario and a geometry, and returns
a `Budget`.

Each budget is built at one of three levels of rigour. Select the level with a
single whole-path `fidelity` argument on the budget (`fidelity=0|1|2`):

- Fidelity 0 is analytic. A Term gives a closed-form loss. It carries no fade.
- Fidelity 1 is statistical. A Term gives samples, so the budget gives a real
  Monte Carlo fade (the FAST modal coupling, the coupled-flux uplink).
- Fidelity 2 is wave optics. It appears as two Terms: a deterministic
  vacuum-optics Term (the full no-turbulence loss) and a stochastic turbulence
  Term (the fade). It needs a precomputed `wave` bundle from
  `olb.models.coupling.run_fidelity2`. The budget never runs the split-step
  simulation itself.

Fidelity 1 does not exist for a terrestrial link.

## 3. A minimal uplink example

This example launches a 1 W beam from a ground beam director to a satellite
bucket receiver.

```python
import numpy as np
from olb import (SpaceScenario, Channel, Site, CircularOrbit,
                 Terminal, Transmitter, Aperture, uplink_budget)

# The launch power sits on the Transmitter. The sensitivity sits on the Detector.
ground_tx = Terminal(
    aperture_m=0.15, wavelength_m=1550e-9, pointing_jitter_rad=1e-6,
    transmitter=Transmitter(waist_m=0.06, power_dbm=30.0),   # 1 W launch
)
space_rx = Terminal(
    aperture_m=0.05, wavelength_m=1550e-9,
    detector=Aperture(sensitivity_dbm=-40.0),                # power-in-bucket
)
scenario = SpaceScenario(
    ground=ground_tx, space=space_rx, direction="uplink",
    channel=Channel(site=Site(cn2_ground=1.7e-14), altitude_m=600e3),
)

budget = uplink_budget(scenario, CircularOrbit(600e3, elevation_deg=60.0))
```

Read the budget four ways:

```python
print(budget.to_frame())                 # the itemised terms (one row per Term)
print(budget.assumptions_frame())        # the model constraints per Term
print(budget.check())                    # the list of broken assumptions

mc = budget.monte_carlo(20000, rng=np.random.default_rng(0),
                        availabilities=(0.99,))
print(mc["mean_loss_db"])                # the mean total loss [dB]
print(mc["fade_db"][0.99])               # the 99% fade level [dB]
print(mc["margin_db"][0.99])             # the 99% link margin [dB]
```

`to_frame()` gives a table with the mean loss of each Term. `assumptions_frame()`
gives the beam type, the turbulence regime, and the validity limit of each Term.
`check()` returns the pairs that break a model assumption; it warns for each.
`monte_carlo()` draws joint samples and gives the mean, the fade, and the margin.

Loss is positive dB. Gain is negative dB.

The uplink budget folds the tracking jitter into the coupled-flux turbulence
Term. So it does not add a separate pointing Term when turbulence is on. For the
full runnable script, see [../validation/uplink_divergence.py](../validation/uplink_divergence.py).

## 4. The four link families

### Uplink (ground to space)

Build a ground transmit terminal and a space receive terminal. Set the direction
to `"uplink"`. Call `uplink_budget(scenario, geometry)`. Pass a diverged beam with
`Transmitter(divergence_rad=...)`; a wider beam scintillates less and points more
easily, but spreads more. See
[../validation/uplink_divergence.py](../validation/uplink_divergence.py).

### Downlink (space to ground)

Build a space transmit terminal and a ground receive terminal. Set the direction
to `"downlink"`. Call `downlink_budget(scenario, geometry)`. The receive detector
selects the receive-side model:

- An `Aperture` (bucket) detector gives the plane-wave scintillation fade.
- An `SMF` (single-mode fibre) detector gives the fibre-coupling loss and fade.
  Add a `compensation` stack (`TipTilt`, `AO`) to clean the wavefront and buy
  back the coupling. Select the fidelity with `fidelity=1` (the statistical
  default) or `fidelity=0` (mean-only, no fade).

See [../examples/downlink_terminal.py](../examples/downlink_terminal.py).

### Retro (space, long slant path)

Build one ground terminal that both transmits and receives. Build a space
terminal as a passive retroreflector (an aperture only). Set the direction to
`"retro"`. Call `retro_space_budget(scenario, geometry)`. The retro re-emits the
power it catches, so the budget stacks an up-leg and a down-leg. For different
transmit and receive apertures on the one ground terminal, set `aperture_m` on
the `Transmitter`. See [../examples/retro_link.py](../examples/retro_link.py).

### Terrestrial (horizontal ground-to-ground path)

Build a near terminal and a far terminal. Use a `TerrestrialChannel` and a
`HorizontalPath` geometry. Call `terrestrial_budget(scenario, geometry)`. An
`Aperture` receiver gives the horizontal Gaussian-beam scintillation fade. An
`SMF` receiver gives the fidelity-0 mean-only coupling loss instead. See
[../examples/terrestrial_link.py](../examples/terrestrial_link.py).

For a bistatic station with different transmit and receive apertures on the
one-way links, see [../examples/build_a_link.py](../examples/build_a_link.py).

## 5. What a fade margin needs

A fade margin needs a stochastic Term. A Term has three faces: a mean
(`mean_db`), an analytic quantile (`quantile_db`), and a sampler (`sample_db`). A
deterministic Term, such as geometric spreading, sets only the mean. A
statistical Term, such as scintillation, gives all three faces. So it carries a
real fade.

A budget refuses a fade when any Term is mean-only (fidelity 0). A mean-only Term
gives the expected loss of a quantity that really fluctuates, but it models no
fade. The terrestrial SMF coupling Term is one example. It locks the whole budget
to fidelity 0. Then `fade_margin_db()` raises a `ValueError`, and `monte_carlo()`
suppresses the fade and reports the mean only. Read `budget.provides_fade` to
test this. Read `budget.total_loss_db()` for the mean loss. Use a statistical
(fidelity-1) coupling model to get the coupling fade.
