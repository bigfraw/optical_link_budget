'''
Uplink simulation example.

This script sets up a ground-to-satellite optical uplink and evaluates it. It
shows the current API:
- build a Scenario and an analytic orbit geometry,
- assemble the uplink budget from olb.links,
- read the itemised terms and the model assumptions,
- check whether the scenario breaks any assumption,
- run a Monte Carlo for the fade and the link margin,
- sweep elevation for the 99 % margin.

The uplink turbulence term is a Monte Carlo model with no closed-form quantile,
so the uplink budget is evaluated with monte_carlo(), not the analytic fade.

Run from the repo root:
    python -m examples.uplink_sim
'''

import warnings

import numpy as np

from olb import Scenario, Link, Site, CircularOrbit, uplink_budget

# The weak-fluctuation guard warns at low elevation. This script reports the
# same information through budget.check(), so silence the duplicate warning.
warnings.simplefilter("ignore")


def main():
    # The ground station transmits. The satellite receives.
    scenario = Scenario(
        link=Link(
            direction="uplink",
            wavelength_m=1550e-9,
            tx_waist_m=0.06,             # ground transmit beam waist w0 [m]
            tx_power_dbm=42,          # 10 W launch power
            rx_diameter_m=0.05,         # satellite receive aperture [m]
            pointing_jitter_rad=2e-6,   # 2 urad tracking jitter
            rx_sensitivity_dbm=-40.0,   # required received power
        ),
        site=Site(cn2_ground=5.7e-14),
        altitude_m=1500e3,
    )
    rng = np.random.default_rng(0)

    # --- one elevation: itemised budget, assumptions, Monte Carlo ----------
    budget = uplink_budget(scenario, CircularOrbit(1500e3, 60.0), n_samples=4000)

    print("Itemised uplink budget at 60 deg elevation:")
    print(budget.to_frame().to_string(index=False))

    print("\nModel assumptions:")
    print(budget.assumptions_frame().to_string(index=False))

    flags = budget.check(warn=False)
    print("\nBroken assumptions:", flags if flags else "none")

    mc = budget.monte_carlo(8000, rng=rng, availabilities=(0.99, 0.999))
    print(f"\nMean loss:            {float(mc['mean_loss_db']):6.2f} dB")
    for a in (0.99, 0.999):
        print(f"Fade / margin at {a*100:5.1f} %:  "
              f"{float(mc['fade_db'][a]):6.2f} dB  /  {float(mc['margin_db'][a]):6.2f} dB")

    # --- elevation sweep: the 99 % margin ----------------------------------
    print("\nElevation sweep (99 % availability):")
    print(" elev   mean-loss   99%-fade   99%-margin   assumptions")
    for elevation_deg in [20, 30, 45, 60, 90]:
        b = uplink_budget(scenario, CircularOrbit(1500e3, float(elevation_deg)),
                          n_samples=3000)
        m = b.monte_carlo(5000, rng=rng, availabilities=(0.99,))
        state = "flagged" if b.check(warn=False) else "ok"
        print(f" {elevation_deg:4d}   {float(m['mean_loss_db']):8.2f}   "
              f"{float(m['fade_db'][0.99]):8.2f}   {float(m['margin_db'][0.99]):9.2f}"
              f"   {state}")


if __name__ == "__main__":
    main()
