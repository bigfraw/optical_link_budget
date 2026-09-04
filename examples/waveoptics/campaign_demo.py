"""A campaign of turbulent trials, stored on disk, read as a budget.

THE FLOW. One scenario, one geometry, one campaign, then the budgets:

    campaign = Campaign(scenario, orbit, ROOT, seed=2024).run(1000, workers=4)
    downlink_budget(scenario, orbit, fidelity=2, wave=campaign)
    multi_detector_budgets(scenario, orbit, arms, fidelity=2, wave=campaign)

A Campaign IS a fidelity-2 wave record, so it goes straight into the `wave`
slot of a budget. The budget calls olb.models.waveoptics.resolve_wave, which
turns the campaign into the bundle (or, with the arms, into the per-arm list).
There is no hand-built bundle and no recouple call in the budget path.

The script builds 1000 downlink snapshots in ten blocks of 100, on a warm pool
of four processes. Run it two times: the second run computes NOTHING, because
the blocks are already on disk. It then reads the store two ways, and it ends
with one DIAGNOSTIC section (a smaller receive aperture, outside the budget
flow).

Run it with:
    python examples/waveoptics/campaign_demo.py
"""

import os
import time

from olb.geometry import CircularOrbit
from olb.links import downlink_budget
from olb.multidetector import multi_detector_budgets
from olb.scenario import Channel, SpaceScenario
from olb.terminal import MMF, SMF, Terminal, Transmitter
from olb.waveoptics.turbulence import Campaign

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    "_campaign_demo2")


def main():
    """Build the campaign, then make the budgets from it."""
    lam = 1550e-9
    design_m = 0.40                       # the design receive aperture.
    scenario = SpaceScenario(
        ground=Terminal(aperture_m=design_m, wavelength_m=lam, detector=SMF(),
                        transmitter=Transmitter(waist_m=0.06)),
        space=Terminal(aperture_m=0.30, wavelength_m=lam,
                       transmitter=Transmitter(waist_m=0.05, power_dbm=30.0)),
        direction="downlink", channel=Channel())
    orbit = CircularOrbit(altitude_m=600e3, elevation_deg=30.0)

    # The patch holds the DESIGN aperture, so every smaller aperture is a crop.
    campaign = Campaign(scenario, orbit, ROOT, seed=2024, preset="rapid",
                        block_size=100, sizing_aperture_m=design_m)
    t0 = time.time()
    n = campaign.run(1000, workers=4, progress=True)
    print(f"{n} trials on disk after {time.time() - t0:.1f} s "
          f"({campaign.grid.n} px, {campaign.plan.z_m.size} screens)")

    # 1. ONE budget of the scenario receive path. The campaign is the record.
    budget = downlink_budget(scenario, orbit, fidelity=2, wave=campaign)
    print("\nfidelity-2 downlink budget from the campaign:")
    print(budget.to_frame().to_string(index=False))
    print(f"  total {float(budget.total_loss_db()):.2f} dB")

    # 2. ONE budget for each beamsplitter arm, from the SAME campaign. The
    # `frac` enters one time, as the fixed splitter Term.
    arms = [SMF(frac=0.7),
            MMF(core_radius_m=25e-6, focal_length_m=0.5)]   # frac None: 0.3
    print("\nper-arm fidelity-2 budgets from the same campaign:")
    for detector, arm_budget in multi_detector_budgets(scenario, orbit, arms,
                                                       fidelity=2,
                                                       wave=campaign):
        split = [t.mean_db for t in arm_budget.terms if t.name == "beamsplitter"]
        print(f"  {type(detector).__name__:<4s} splitter "
              f"{(split[0] if split else 0.0):5.2f} dB   total "
              f"{float(arm_budget.total_loss_db()):7.2f} dB")

    # 3. DIAGNOSTIC, outside the budget flow. The stored field does not know the
    # receive aperture, so a smaller aperture is a post-hoc crop of that field.
    eta = campaign.recouple(SMF(), aperture_m=0.20)
    print(f"\ndiagnostic (not a budget): D = 20.0 cm SMF coupling, mean eta "
          f"{eta.mean():.5f}, min {eta.min():.5f}")


if __name__ == '__main__':
    main()
