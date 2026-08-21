'''
Downlink receive-terminal example.

The satellite transmits. The ground station detects with a telescope, an optional
adaptive-optics stage, and a detector front end. This script shows how the
receive Terminal changes the downlink budget:

- Aperture (power-in-bucket): parity with the plain downlink. The detector is
  phase-insensitive, so wavefront compensation does not change it.
- SMF (single-mode fibre) with NO correction: the fibre couples only the field
  that matches the fibre mode. Uncorrected turbulence distorts the wavefront, so
  the coupling loss is large (the Dikmelik-Davidson D/r0 regime).
- SMF with tip-tilt, then SMF with adaptive optics: each stage removes more
  wavefront error, so the residual variance falls and the coupling recovers
  (the extended-Marechal regime, eta = eta_max * exp(-residual variance)).

The receive Terminal is opt-in. When scenario.rx_terminal has a detector, the
receive-coupling Term owns the receive-side turbulence and REPLACES the plain
scintillation Term. When rx_terminal is None, the downlink budget is unchanged.

Run from the repo root:
    python -m examples.downlink_terminal
'''

import warnings

import numpy as np

from olb import (Scenario, Site, Channel, CircularOrbit, downlink_budget,
                 Terminal, Transmitter, Aperture, SMF, TipTilt, AO)

warnings.simplefilter("ignore")

APERTURE_M = 0.7   # ground telescope diameter [m]
WAVELENGTH_M = 1550e-9
SENSITIVITY_DBM = -45.0


def _ground(detector, compensation=()):
    '''Build the ground receive terminal: shared telescope, varied detector.'''
    return Terminal(APERTURE_M, wavelength_m=WAVELENGTH_M,
                    detector=detector, compensation=list(compensation))


def main():
    # The satellite transmits; the ground station receives. All hardware lives
    # on the two terminals.
    space = Terminal(
        aperture_m=0.05,            # satellite telescope [m]
        wavelength_m=WAVELENGTH_M,
        pointing_jitter_rad=1e-6,
        transmitter=Transmitter(waist_m=0.035, power_dbm=30.0),  # 1 W downlink
    )
    site = Site(cn2_ground=1.7e-14)
    altitude_m = 600e3
    geom = CircularOrbit(altitude_m, elevation_deg=45.0)

    # Each receive terminal shares the same telescope; only the detector and the
    # compensation stack change.
    terminals = [
        ("aperture (bucket)",   _ground(Aperture(sensitivity_dbm=SENSITIVITY_DBM))),
        ("SMF, no correction",  _ground(SMF(sensitivity_dbm=SENSITIVITY_DBM))),
        ("SMF + tip-tilt",      _ground(SMF(sensitivity_dbm=SENSITIVITY_DBM),
                                        [TipTilt()])),
        ("SMF + AO(60)",        _ground(SMF(sensitivity_dbm=SENSITIVITY_DBM),
                                        [TipTilt(), AO(n_modes=60)])),
        ("SMF + AO(200)",       _ground(SMF(sensitivity_dbm=SENSITIVITY_DBM),
                                        [TipTilt(), AO(n_modes=200)])),
    ]

    rows = []
    for label, ground in terminals:
        scn = Scenario(ground=ground, space=space, direction="downlink",
                       channel=Channel(site=site, altitude_m=altitude_m))
        budget = downlink_budget(scn, geom)

        # The receive-coupling Term carries the receive-side turbulence.
        coupling = next(t for t in budget.terms if t.category == "coupling")

        print("=" * 62)
        print(label)
        print(budget.to_frame()[["name", "mean_db"]].to_string(index=False))

        mc = budget.monte_carlo(8000, rng=np.random.default_rng(0),
                                availabilities=(0.99,))
        fade99 = float(mc["fade_db"][0.99])
        print(f"coupling loss {coupling.mean_db:6.2f} dB   "
              f"99% fade {fade99:6.2f} dB\n")
        rows.append((label, coupling.mean_db, float(mc["mean_loss_db"]), fade99))

    print("=" * 62)
    print("Summary (45 deg elevation, 0.7 m telescope)")
    print(f"{'receive terminal':<22}{'coupling':>10}{'total':>9}{'99% fade':>12}")
    for label, coupling_db, total_db, fade99 in rows:
        print(f"{label:<22}{coupling_db:>10.2f}{total_db:>9.2f}{fade99:>12.2f}")
    
    print("\nA single-mode fibre needs the wavefront cleaned up: with no "
          "correction it costs many dB, but adaptive optics recovers most of it.")
    


if __name__ == "__main__":
    main()
