'''
How to build a link with the terminal model.

This example shows the shape of the package. You build two Terminals and a
Channel, put them in a Scenario, and set a direction. Hardware lives ONLY on the
Terminals: a Terminal owns its telescope aperture, its operating wavelength, its
tracking jitter, and up to three optional parts:

- a Transmitter (it can launch a beam),
- a Detector (it can receive: an Aperture bucket or a single-mode fibre),
- a Compensation stack (receive-side wavefront correction: tip-tilt, AO).

The Channel is the propagation path (ground site + satellite orbit). The
direction ("uplink" or "downlink") chooses which terminal transmits and which
receives. So the SAME two terminals describe both links. This script builds one
ground station and one satellite, then runs the link BOTH ways by flipping only
the direction -- no new hardware.

Run from the repo root:
    python -m examples.build_a_link
'''

import warnings
from dataclasses import replace

import numpy as np

from olb import (Scenario, Channel, Site, CircularOrbit, Terminal, Transmitter,
                 Aperture, SMF, TipTilt, AO, uplink_budget, downlink_budget)

warnings.simplefilter("ignore")   # budget.check() reports the same weak-fluctuation flags


def main():
    wavelength_m = 1550e-9

    # --- 1. The ground station (optical ground station, OGS) ---------------
    # A big telescope. It transmits the uplink beam, and it receives the
    # downlink into a single-mode fibre behind an adaptive-optics stage.
    ground = Terminal(
        aperture_m=0.7,
        wavelength_m=wavelength_m,
        pointing_jitter_rad=2e-6,
        transmitter=Transmitter(waist_m=0.06, power_dbm=42.0,   # ~16 W launch
                                divergence_rad=15e-6),          # deliberate divergence
        detector=SMF(sensitivity_dbm=-45.0),                    # coherent / fibre front end
        compensation=[TipTilt(), AO(n_modes=60)],               # clean the wavefront for the fibre
    )

    # --- 2. The satellite --------------------------------------------------
    # A small terminal. It transmits the downlink, and it receives the uplink
    # with a plain aperture (power-in-bucket) detector. No adaptive optics.
    space = Terminal(
        aperture_m=0.1,
        wavelength_m=wavelength_m,
        pointing_jitter_rad=1e-6,
        transmitter=Transmitter(waist_m=0.05, power_dbm=30.0),  # 1 W downlink
        detector=Aperture(sensitivity_dbm=-40.0),
    )

    # --- 3. The channel and the scenario -----------------------------------
    channel = Channel(site=Site(cn2_ground=1.7e-14), altitude_m=600e3)
    uplink = Scenario(ground=ground, space=space, direction="uplink",
                      channel=channel)
    geom = CircularOrbit(channel.altitude_m, elevation_deg=45.0)

    # The direction resolves the roles. Read them back to make the point:
    assert uplink.tx_terminal is ground and uplink.rx_terminal is space

    # --- 4. Run the link BOTH ways -----------------------------------------
    # Uplink: the ground transmits, the satellite receives. The beam fights the
    # ground turbulence right at launch (wander + scintillation).
    up = uplink_budget(uplink, geom, n_samples=4000)

    # The SAME terminals, the SAME channel -- only the direction flips.
    # Downlink: the satellite transmits, the ground receives with AO + fibre.
    downlink = replace(uplink, direction="downlink")
    assert downlink.tx_terminal is space and downlink.rx_terminal is ground
    down = downlink_budget(downlink, geom)

    rng = np.random.default_rng(0)
    for name, budget in (("UPLINK  (ground -> satellite)", up),
                         ("DOWNLINK (satellite -> ground)", down)):
        print("=" * 62)
        print(name)
        print(budget.to_frame()[["name", "mean_db"]].to_string(index=False))
        mc = budget.monte_carlo(8000, rng=rng, availabilities=(0.99,))
        print(f"mean loss {float(mc['mean_loss_db']):6.2f} dB   "
              f"99% margin {float(mc['margin_db'][0.99]):6.2f} dB\n")

    print("=" * 62)
    print("One scenario, one flag flipped. The uplink and the downlink use the "
          "same two terminals but see very different budgets: the uplink pays "
          "for launching through the ground turbulence, while the downlink pays "
          "for coupling the collected light into the fibre (which the AO buys "
          "back).")


if __name__ == "__main__":
    main()
