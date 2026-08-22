'''
A terrestrial (horizontal-path) link: ground-to-ground over a fixed range.

This example builds the fidelity-zero terrestrial budget. The link runs along a
horizontal path between two ground terminals. There is no satellite, no orbit,
and no elevation angle. The range is the path length.

A terrestrial link uses a TerrestrialScenario, not a SpaceScenario. Both ends
are on the ground, so the two terminals are named for the path ends:

    near  -- the local end   (the transmitter)
    far   -- the remote end  (the receiver)

The channel is a TerrestrialChannel. It carries the path length, the horizontal
extinction coefficient (dB/km), and the constant Cn2 along the path. The
geometry is a HorizontalPath. It gives the range only.

The budget holds three deterministic Terms: the geometric spreading, the
horizontal Beer-Lambert extinction, and the pointing jitter. These Terms are
exact and direction-agnostic, so olb reuses them from the space links.

SCINTILLATION IS NOT INCLUDED. A horizontal-path scintillation model is steered
by the Gaussian-beam parameters at the receiver, so it is not the plane-wave
slant-path model that the downlink uses. That analytic Term is a reserved slot
(olb.links.terrestrial.terrestrial_scintillation_term). It raises
NotImplementedError and names the Andrews equations that it needs. So this
example keeps scintillation off (the default).

Run from the repo root:
    python -m examples.terrestrial_link
'''

import numpy as np

from olb import (TerrestrialScenario, TerrestrialChannel, Site, HorizontalPath,
                 Terminal, Transmitter, Aperture, terrestrial_budget)


def main():
    wavelength_m = 1550e-9
    path_length_m = 3e3          # 3 km horizontal link

    # --- 1. The two ground terminals --------------------------------------
    # The near (local) terminal launches the beam through a small director. A
    # wide aperture (0.2 m for a 0.02 m waist) does not truncate the beam, so
    # the launch-truncation Term does not fire.
    near = Terminal(
        aperture_m=0.2,
        wavelength_m=wavelength_m,
        pointing_jitter_rad=5e-6,                     # 5 urad tracking jitter
        transmitter=Transmitter(waist_m=0.02, power_dbm=30.0),   # 1 W launch
    )
    # The far (remote) terminal receives the beam into a power-in-bucket
    # detector.
    far = Terminal(
        aperture_m=0.2,
        wavelength_m=wavelength_m,
        detector=Aperture(sensitivity_dbm=-40.0),
    )

    # --- 2. The horizontal channel and the geometry -----------------------
    # The extinction coefficient is weather- and visibility-dependent. Set it
    # per site. 0.5 dB/km is a clear-air value. The constant Cn2 feeds the
    # (pending) scintillation Term only, so it does not change this budget.
    channel = TerrestrialChannel(
        site=Site(),
        path_length_m=path_length_m,
        attenuation_db_per_km=0.5,
        cn2=1e-14,
    )
    scenario = TerrestrialScenario(near=near, far=far, channel=channel)

    # The direction resolves the roles: near transmits, far receives.
    assert scenario.tx_terminal is near and scenario.rx_terminal is far

    geometry = HorizontalPath(path_length_m)

    # --- 3. Run the budget ------------------------------------------------
    budget = terrestrial_budget(scenario, geometry)   # scintillation off (default)

    print("=" * 62)
    print(f"Terrestrial link: {path_length_m / 1e3:.1f} km horizontal path")
    print("=" * 62)
    print(budget.to_frame()[["name", "category", "mean_db"]].to_string(index=False))

    total = float(budget.total_loss_db())
    fade_99 = float(budget.fade_margin_db(0.99))
    print(f"\ntotal mean loss   {total:6.2f} dB")
    print(f"99% analytic fade {fade_99:6.2f} dB")

    # The near terminal carries the launch power; the far terminal carries the
    # sensitivity. So the budget reports the link margin.
    mc = budget.monte_carlo(20000, rng=np.random.default_rng(0),
                            availabilities=(0.99,))
    print(f"99% link margin   {float(mc['margin_db'][0.99]):6.2f} dB")

    # --- 4. A path-length sweep -------------------------------------------
    # The geometry and the channel both take an array, so the budget broadcasts
    # over the range. A longer path costs more spread AND more extinction.
    lengths_km = np.array([1.0, 3.0, 5.0, 10.0])
    sweep_scn = TerrestrialScenario(
        near=near, far=far,
        channel=TerrestrialChannel(path_length_m=lengths_km * 1e3,
                                   attenuation_db_per_km=0.5, cn2=1e-14))
    sweep = terrestrial_budget(sweep_scn, HorizontalPath(lengths_km * 1e3))
    print("\n" + "=" * 62)
    print("Path-length sweep (mean loss):")
    for L, loss in zip(lengths_km, np.atleast_1d(sweep.total_loss_db())):
        print(f"  {L:5.1f} km  ->  {loss:6.2f} dB")

    print("\n" + "=" * 62)
    print("Scintillation is OFF. A horizontal-path scintillation model needs the "
          "Andrews Gaussian-beam forms (see terrestrial_scintillation_term). Add "
          "them, then set scintillation=True to include the turbulence fade.")


if __name__ == "__main__":
    main()
