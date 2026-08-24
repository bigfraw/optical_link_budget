'''
End-to-end example: the code composes all four models into one Budget.

This is the exact check that produced the integration table. It assembles the
geometric, atmospheric, pointing, and coupled-flux turbulence terms into a
single Budget. It then evaluates the Budget by Monte Carlo for mean loss, fade,
and margin.

Run from the repo root:
    python -m examples.end_to_end

It builds a FAST-free Cn2 profile from general_atmospherics.get_c2n. It runs
whether or not the `fast` package is installed.
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
    # FAST-free Cn2 profile (Hufnagel-Valley via get_c2n) so this runs anywhere
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

    budget = Budget([
        geometric_loss_term(scenario, geom),
        slant_extinction_term(scenario, geom),
        pointing_loss_term(scenario, geom),
        uplink_turbulence_term(scenario, geom, n_samples=4000, cn2_profile=cn2),
    ], scenario=scenario)

    print(budget.to_frame().to_string(index=False))
    print()

    mc = budget.monte_carlo(4000, rng=np.random.default_rng(0),
                            availabilities=(0.5, 0.99))
    print(f"mean total loss dB : {float(mc['mean_loss_db']):.2f}")
    print(f"fade 99% dB        : {float(mc['fade_db'][0.99]):.2f}")
    print(f"margin 99% dB      : {float(mc['margin_db'][0.99]):.2f}")

    # check: the itemised table has all four terms and the margin is finite
    assert len(budget.terms) == 4
    assert np.isfinite(mc["margin_db"][0.99])
    print("\nend-to-end OK")


if __name__ == "__main__":
    main()
