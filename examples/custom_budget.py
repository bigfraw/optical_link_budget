'''
How to assemble a CUSTOM Budget by hand, Term by Term.

The package gives four pre-built budgets in `olb/links`: `uplink_budget`,
`downlink_budget`, `retro_space_budget`, and `terrestrial_budget`. Use one of
them first. Each one selects the correct Terms for its link family.

Build a Budget by hand only when a pre-built budget does not fit:

- The link is not one of the four families. An example is a new geometry, or a
  chain of legs that no factory covers.
- You must SWAP one Term. An example is a different turbulence model, or a Term
  from your own code.
- You want a what-if. Drop a Term, or add a Term, and see the effect on the
  total.

A Budget is only a list of Terms plus the scenario. Each Term factory has the
same shape: `f(scenario, geometry) -> Term`. So you can put any set of Terms in
one Budget. The Budget then gives the same faces as a pre-built budget:
`to_frame()`, `check()`, and `monte_carlo()`.

This script assembles the uplink from its four Terms: the geometric spreading,
the slant extinction, the pointing jitter, and the coupled-flux turbulence.

CAUTION: a hand-built Budget does not protect you from a double count. The
pre-built `uplink_budget` folds the tracking jitter INTO the coupled-flux
turbulence Term when the turbulence is on, and it then drops the separate
pointing Term. This script keeps the two Terms apart, so it sets the ground
jitter to zero.

Run from the repo root:
    python -m examples.custom_budget

It builds a FAST-free Cn2 profile from `get_c2n`. It runs whether or not the
optional `fast` package is installed.
'''

import numpy as np

from olb import (SpaceScenario, Site, Channel, CircularOrbit, Budget,
                 Terminal, Transmitter, Aperture)
from olb.turbulence.profiles import DEFAULT_HS, get_c2n
from olb.models.geometric import geometric_loss_term
from olb.models.extinction import slant_extinction_term
from olb.models.pointing import pointing_loss_term
from olb.links.uplink import uplink_turbulence_term


def main():
    # A FAST-free Hufnagel-Valley Cn2 profile, so this script runs anywhere.
    cn2 = get_c2n(DEFAULT_HS, 21, 1.7e-14)

    scenario = SpaceScenario(
        ground=Terminal(aperture_m=0.5, wavelength_m=1550e-9,
                        pointing_jitter_rad=0,
                        transmitter=Transmitter(waist_m=0.1, power_dbm=40.0)),
        space=Terminal(aperture_m=0.08, wavelength_m=1550e-9,
                       detector=Aperture(sensitivity_dbm=-40.0)),
        direction="uplink",
        channel=Channel(site=Site(cn2_ground=1.7e-14), altitude_m=1500e3),
    )
    geom = CircularOrbit(1500e3, 60.0)   # single elevation

    # The Budget is the list of Terms. Add a Term, or remove a Term, to make
    # your own budget. Each factory reads the same scenario and geometry.
    budget = Budget([
        geometric_loss_term(scenario, geom),
        slant_extinction_term(scenario, geom),
        pointing_loss_term(scenario, geom),
        uplink_turbulence_term(scenario, geom, n_samples=4000, cn2_profile=cn2),
    ], scenario=scenario)

    print("Itemised custom budget at 60 deg elevation:")
    print(budget.to_frame().to_string(index=False))

    # check() flags a scenario that breaks a Term assumption. A hand-built
    # Budget needs this more than a pre-built one, because no factory selected
    # the Terms for you.
    flags = budget.check(warn=False)
    print("\nBroken assumptions:", flags if flags else "none")
    print()

    # The turbulence Term is a Monte Carlo model, so read the fade from
    # monte_carlo(), not from the analytic quantile.
    mc = budget.monte_carlo(4000, rng=np.random.default_rng(0),
                            availabilities=(0.5, 0.99))
    print(f"mean total loss dB : {float(mc['mean_loss_db']):.2f}")
    print(f"fade 99% dB        : {float(mc['fade_db'][0.99]):.2f}")
    print(f"margin 99% dB      : {float(mc['margin_db'][0.99]):.2f}")

    # check: the itemised table has all four terms and the margin is finite
    assert len(budget.terms) == 4
    assert np.isfinite(mc["margin_db"][0.99])
    print("\ncustom budget OK")


if __name__ == "__main__":
    main()
