'''
How to build a retroreflected ground-to-space link.

A retro link is different from the one-way uplink and downlink. ONE ground
station both transmits the up-leg and receives the return, and the satellite is a
passive retroreflector (an aperture only, with no transmitter or detector). The
retro direction makes tx_terminal and rx_terminal the SAME object, so the ground
is ONE Terminal that carries both a Transmitter and a Detector.

That single Terminal can still be BISTATIC: a small beam director transmits the
up-leg while a large telescope receives the return. You cannot use two Terminals
here (retro forces one ground object), so the Transmitter carries its own
aperture_m. The Terminal aperture_m is the receive telescope; the Transmitter
aperture_m is the beam director. See olb.terminal.Transmitter. (For the two-
Terminal bistatic pattern on the one-way links, see examples/build_a_link.py.)

Retroreflection is a retransmission. The ground launches a beam up. The retro
captures the power over its aperture and re-emits it back down as a new beam. So
the budget is the up-leg plus the down-leg -- roughly the two one-way budgets
stacked. This is the SPACE model (long slant range, fully diverged return,
independent turbulence on the two legs). It does not hold for a short terrestrial
retro link.

Run from the repo root:
    python -m examples.retro_link
'''
import numpy as np

from olb import (Scenario, Channel, Site, CircularOrbit, Terminal, Transmitter,
                 SMF, TipTilt, AO, retro_space_budget)


def main():
    wavelength_m = 1550e-9

    # --- 1. The ground station (one Terminal, transmit AND receive) -------
    # A BISTATIC OGS in one Terminal: a large telescope (aperture_m,
    # obscuration_ratio) receives the return into a single-mode fibre, and a
    # separate small beam director launches the up-leg. The Transmitter carries
    # the beam-director aperture, so the launch truncation reads 0.15 m, not the
    # 0.7 m receive telescope. Retro forces one ground object, so this is the
    # only way to give the station different transmit and receive apertures.
    # The return couples into a single-mode fibre. The budget below runs with
    # smf_fidelity="fast", so the coupling is the FAST fidelity-1 true LP01 modal
    # overlap (needs fast-aosim). A compensation stack maps to the FAST correction:
    # TipTilt -> tip-tilt mode, AO(n_modes) -> modal AO with ZMAX=n_modes, and an
    # empty stack -> no correction. For the bucket-detector return, swap the
    # detector for Aperture(sensitivity_dbm=-50); for the analytic mean-only
    # coupling loss (no fade), set smf_fidelity="mean".
    ground = Terminal(
        aperture_m=0.7,                  # receive telescope
        obscuration_ratio=0.3,           # 30% central obscuration (receive)
        wavelength_m=wavelength_m,
        pointing_jitter_rad=0,
        transmitter=Transmitter(waist_m=0.06, power_dbm=40.0, divergence_rad=None,
                                aperture_m=0.15,          # small beam director
                                obscuration_ratio=0.3),   # obscured director
        detector=SMF(sensitivity_dbm=-110.0),             # coherent / fibre return
        compensation=[TipTilt()]#, AO(n_modes=2)],         # AO for fibre coupling
    )

    # --- 2. The satellite (a passive retroreflector) ----------------------
    # Only the aperture and wavelength are used. Its transmitter/detector, if
    # any, are ignored: the retro just re-emits the power it catches.
    space = Terminal(
        aperture_m=0.05,
        wavelength_m=wavelength_m,
    )

    # --- 3. The channel and the scenario ----------------------------------
    channel = Channel(site=Site(cn2_ground=1.7e-14), altitude_m=1500e3)
    retro = Scenario(ground=ground, space=space, direction="retro",
                     channel=channel)

    for elevation_deg in (30.0, 45.0, 90.0):
        geom = CircularOrbit(channel.altitude_m, elevation_deg=elevation_deg)

        # The retro direction makes the ground both transmit and receive.
        assert retro.tx_terminal is ground and retro.rx_terminal is ground

        # --- 4. Run the retro link --------------------------------------------
        # The budget carries both legs: the up-leg (ground -> retro) and the
        # down-leg (retro -> ground), because the retro re-transmits the power it
        # catches.
        retro_b = retro_space_budget(retro, geom, n_samples=4000, smf_fidelity="fast")

        rng = np.random.default_rng(0)
        print("=" * 62)
        print(f"RETRO   (ground -> retro -> ground), {elevation_deg:5.1f} deg")
        print(retro_b.to_frame()[["name", "mean_db"]].to_string(index=False))
        mc = retro_b.monte_carlo(8000, rng=rng, availabilities=(0.99,))
        print(f"mean loss {float(mc['mean_loss_db']):6.2f} dB   "
            f"fade 99% {float(mc['fade_db'][0.99]):6.2f} dB\n")

        print("=" * 62)

        af = retro_b.assumptions_frame()
        broken = af[~af["ok"]]
        print("assumptions:", "all hold" if broken.empty else "BROKEN below")
        if not broken.empty:
            print(broken[["name", "violations"]].to_string(index=False))


if __name__ == "__main__":
    main()
