'''
How to build a link when the transmit and receive apertures are DIFFERENT.

This example shows the shape of the package for a BISTATIC station: the terminal
that transmits and the terminal that receives are separate hardware with
different apertures. Hardware lives ONLY on a Terminal: a Terminal owns its
telescope aperture, its operating wavelength, its tracking jitter, and up to
three optional parts:

- a Transmitter (it can launch a beam),
- a Detector (it can receive: an Aperture bucket or a single-mode fibre),
- a Compensation stack (receive-side wavefront correction: tip-tilt, AO).

A monostatic station shares ONE aperture for both directions, so ONE Terminal
serves both links and you flip only the direction. A BISTATIC station does not:
the transmit beam director and the receive telescope are different apertures. An
aperture is a Terminal parameter, so a different aperture needs a different
Terminal. So each role gets its own Terminal, and each direction gets its own
SpaceScenario.

This script builds four terminals -- a transmit and a receive terminal at the
ground, and a transmit and a receive terminal in space -- then wires the correct
pair into each direction:

    uplink   : ground_tx  transmits  ->  space_rx  receives
    downlink : space_tx   transmits  ->  ground_rx receives

For the retroreflected link (one ground telescope both transmits and receives),
see examples/retro_link.py.

Run from the repo root:
    python -m examples.build_a_link
'''

import numpy as np

from olb import (SpaceScenario, Channel, Site, CircularOrbit, Terminal, Transmitter,
                 Aperture, SMF, TipTilt, AO, uplink_budget, downlink_budget)
from olb.models.waveoptics import run_fidelity2
from olb.waveoptics.threader import Threader


def main():
    wavelength_m = 1550e-9

    # --- 1. The ground station: two separate telescopes -------------------
    # A bistatic optical ground station (OGS). A small beam director launches
    # the uplink. A separate, much larger telescope receives the downlink into
    # a single-mode fibre behind an adaptive-optics stage. The two apertures
    # differ, so they are two Terminals.
    ground_tx = Terminal(
        aperture_m=0.15,                 # small beam director
        obscuration_ratio=0.3,
        wavelength_m=wavelength_m,
        pointing_jitter_rad=1e-6,
        transmitter=Transmitter(waist_m=0.06, power_dbm=30.0,
                                divergence_rad=None),   # deliberate divergence
    )
    ground_rx = Terminal(
        aperture_m=0.7,                  # large receive telescope
        obscuration_ratio=0.3,           # 30% central obscuration
        wavelength_m=wavelength_m,
        pointing_jitter_rad=1e-6,
        detector=SMF(sensitivity_dbm=-110.0),            # coherent / fibre front end
        compensation=[TipTilt()],        # clean the wavefront for the fibre
    )

    # --- 2. The satellite: two separate terminals -------------------------
    # A small downlink transmitter and a small uplink receiver (a plain
    # power-in-bucket detector, no adaptive optics). Again, two apertures, so
    # two Terminals.
    space_tx = Terminal(
        aperture_m=0.08,
        wavelength_m=wavelength_m,
        pointing_jitter_rad=1e-6,
        transmitter=Transmitter(waist_m=0.03, power_dbm=30.0),  # 1 W downlink
    )
    space_rx = Terminal(
        aperture_m=0.05,
        wavelength_m=wavelength_m,
        pointing_jitter_rad=1e-6,
        detector=Aperture(sensitivity_dbm=-40.0),
    )

    # The apertures differ. That is why the roles cannot share a Terminal.
    assert ground_tx.aperture_m != ground_rx.aperture_m
    assert space_tx.aperture_m != space_rx.aperture_m

    # --- 3. The channel and one SpaceScenario per direction ---------------
    # The channel (site + orbit) is the same for both directions. The terminals
    # are not: each direction wires in its own transmit and receive Terminal.
    channel = Channel(site=Site(cn2_ground=1.7e-14), altitude_m=1500e3)

    geom = CircularOrbit(channel.altitude_m, elevation_deg=45.0)

    uplink = SpaceScenario(ground=ground_tx, space=space_rx, direction="uplink",
                    channel=channel)
    downlink = SpaceScenario(ground=ground_rx, space=space_tx, direction="downlink",
                        channel=channel)

    # The direction resolves the roles onto the terminals you supplied.
    assert uplink.tx_terminal is ground_tx and uplink.rx_terminal is space_rx
    assert downlink.tx_terminal is space_tx and downlink.rx_terminal is ground_rx

    # A small helper: print each direction's Term table and a Monte Carlo
    # summary. Both budget families read the SAME interface, so one reporter
    # serves fidelity 0/1 and fidelity 2 with no change.
    rng = np.random.default_rng(0)

    def report(title, up, down):
        print("#" * 62)
        print(f"# {title}")
        for name, budget in (("UPLINK  (ground beam director -> satellite)", up),
                             ("DOWNLINK (satellite -> ground telescope)", down)):
            print("=" * 62)
            print(f'{geom.elevation_deg:5.0f} deg elevation | {name}')
            print(budget.to_frame()[["name", "mean_db"]].to_string(index=False))
            print(budget.assumptions_frame())
            mc = budget.monte_carlo(8000, rng=rng, availabilities=(0.99,))
            # A mean-only (fidelity-0) coupling Term, for example the analytic
            # single-mode fibre coupling, models the expected loss but not the
            # fade. The budget then suppresses the fade rather than understate the
            # tail. monte_carlo reports this with the fade_available flag
            # (Budget.provides_fade), so read the flag, not fade_db is None.
            if not mc["fade_available"]:
                print(f"mean loss {float(mc['mean_loss_db']):6.2f} dB   "
                      f"fade 99% N/A -- mean-only fidelity-0 coupling, no fade\n")
            else:
                print(f"mean loss {float(mc['mean_loss_db']):6.2f} dB   "
                      f"fade 99% {float(mc['fade_db'][0.99]):6.2f} dB\n")

    # --- 4. The DEFAULT budgets (fidelity 0/1) ----------------------------
    # The everyday path: no wave-optics simulation. The uplink is the
    # coupled-flux Monte Carlo (fidelity 1, a real fade for the satellite
    # bucket, no FAST needed). The downlink is the FAST modal-overlap fibre
    # coupling (fidelity 1); if fast-aosim is not installed it falls back to the
    # analytic mean-only coupling (fidelity 0).
    up1 = uplink_budget(uplink, geom, fidelity=1, n_samples=4000)
    try:
        down1 = downlink_budget(downlink, geom, fidelity=1)
    except ImportError:
        print("fast-aosim not installed; the downlink uses fidelity 0 "
              "(analytic mean-only coupling).")
        down1 = downlink_budget(downlink, geom, fidelity=0)
    report("fidelity 0/1  (analytic / statistical, no wave optics)", up1, down1)

    # --- 5. The WAVE-OPTICS budgets (fidelity 2) --------------------------
    # fidelity=2 does NOT run the split-step sim implicitly. Precompute the wave
    # bundle ONCE per direction, then pass it as wave=. For a SPACE link the
    # geometric loss is ANALYTIC by default (a ground-space link is far field),
    # so run_fidelity2 SKIPS the vacuum propagation and only the turbulent Monte
    # Carlo runs. Pass vacuum="wave" to opt back into the wave-optics vacuum Term.
    # print('running fidelity=2 waveoptics for uplink...')

    th = Threader(max_workers=8)
    up2_wave = run_fidelity2(uplink, geom, n_trials=200, seed=0, threader=th)
    up2 = uplink_budget(uplink, geom, fidelity=2, wave=up2_wave)

    print('running fidelity=2 waveoptics for downlink...')
    down2_wave = run_fidelity2(downlink, geom, n_trials=200, seed=0, threader=th)
    down2 = downlink_budget(downlink, geom, fidelity=2, wave=down2_wave)
    report("fidelity 2  (wave optics)", up2, down2)

    print("=" * 62)
    print("A bistatic station transmits and receives through DIFFERENT "
          "apertures. An aperture is a Terminal parameter, so each role is its "
          "own Terminal: a small ground beam director launches the uplink, and "
          "a large ground telescope receives the downlink. Because the roles do "
          "not share hardware, each direction is its own SpaceScenario -- you wire in "
          "the transmit and receive Terminals for that link, not flip one flag.")


if __name__ == "__main__":
    main()
