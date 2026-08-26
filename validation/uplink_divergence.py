'''
Uplink beam-divergence example.

A ground station can widen (diverge) its transmit beam on purpose. This script
shows the trade that the divergence makes on a satellite uplink:

- geometric spreading loss RISES (a wider beam spills more power past the small
  receive aperture),
- pointing-jitter loss FALLS (a wider beam is less sensitive to tracking error),
- turbulence loss FALLS (a wider, more spherical-wave-like beam broadens and
  scintillates less; the divergence feeds the beam broadening AND the Dios
  scintillation index).

The deep-fade tail (the 99 % fade) tightens most, so a moderate divergence can
IMPROVE the 99 % link margin even though the mean loss goes up.

The transmit divergence is the far-field 1/e^2 HALF-angle. It cannot be smaller
than the diffraction limit lambda / (pi * w0). None means collimated.

Run from the repo root:
    python -m validation.uplink_divergence
'''

import warnings
from dataclasses import replace

import numpy as np

from olb import (SpaceScenario, Site, Channel, CircularOrbit, uplink_budget,
                 Terminal, Transmitter, Aperture)

# The weak-fluctuation guard warns at low elevation. budget.check() reports the
# same information, so silence the duplicate warning.
warnings.simplefilter("ignore")


def main():
    wavelength_m = 1550e-9
    waist_m = 0.06                       # ground transmit beam waist w0 [m]
    # The ground station transmits; the satellite receives. Only the transmit
    # divergence changes between cases.
    ground = Terminal(
        aperture_m=0.4,                  # ground telescope; wide, no launch truncation
        wavelength_m=wavelength_m,
        pointing_jitter_rad=2e-6,        # 2 urad tracking jitter
        transmitter=Transmitter(waist_m=waist_m, power_dbm=42.0),  # ~16 W launch
    )
    space = Terminal(
        aperture_m=0.05,                 # satellite receive aperture [m]
        wavelength_m=wavelength_m,
        detector=Aperture(sensitivity_dbm=-40.0),   # required received power
    )
    site = Site(cn2_ground=5.7e-14)
    altitude_m = 1500e3
    geom = CircularOrbit(altitude_m, elevation_deg=30.0)

    # The diffraction-limited half-angle sets the floor for the divergence.
    theta_min = wavelength_m / (np.pi * waist_m)
    print(f"Aperture w0 = {waist_m} m  ->  diffraction-limited "
          f"half-angle {theta_min * 1e6:.2f} urad\n")

    # Collimated, then two deliberate divergences.
    cases = [("collimated", None), ("15 urad", 15e-6), ("30 urad", 30e-6), ("60 urad", 60e-6)]

    rows = []
    for label, divergence in cases:
        g = replace(ground,
                    transmitter=replace(ground.transmitter, divergence_rad=divergence))
        scn = SpaceScenario(ground=g, space=space, direction="uplink",
                       channel=Channel(site=site, altitude_m=altitude_m))
        budget = uplink_budget(scn, geom, n_samples=4000)

        print("=" * 62)
        header = label if divergence is None else \
            f"{label}  ({divergence / theta_min:.1f}x the diffraction limit)"
        print(header)
        frame = budget.to_frame()[["name", "mean_db"]]
        print(frame.to_string(index=False))

        flags = budget.check(warn=False)
        print("broken assumptions:", "none" if not flags else flags)

        mc = budget.monte_carlo(8000, rng=np.random.default_rng(0),
                                availabilities=(0.99,))
        mean_loss = float(mc["mean_loss_db"])
        fade99 = float(mc["fade_db"][0.99])
        margin99 = float(mc["margin_db"][0.99])
        print(f"mean loss {mean_loss:6.2f} dB   99% fade {fade99:6.2f} dB   "
              f"99% margin {margin99:6.2f} dB\n")
        rows.append((header, mean_loss, fade99, margin99))

    # Compact side-by-side of the trade.
    print("=" * 62)
    print("Summary (60 deg elevation)")
    print(f"{'case':<34}{'mean':>8}{'99% fade':>11}{'99% margin':>12}")
    for header, mean_loss, fade99, margin99 in rows:
        print(f"{header:<34}{mean_loss:>8.2f}{fade99:>11.2f}{margin99:>12.2f}")

    best = max(rows, key=lambda r: r[3])
    print("\nDiverging trades a higher mean loss for a tighter fade tail, so the "
          "best 99% margin sits at an OPTIMUM divergence, not the widest beam.")
    print(f"Here '{best[0]}' wins the best 99% margin ({best[3]:.2f} dB); "
          "widening further over-spreads the beam and the geometric loss takes over.")


if __name__ == "__main__":
    main()
