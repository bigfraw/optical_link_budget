"""The point-ahead anisoplanatism of a pre-compensated uplink (fidelity 2).

THE PHYSICS. A ground station senses the DOWNLINK beacon, it pre-distorts the
uplink beam with that measurement, and it launches the beam AHEAD of the
beacon, to where the satellite will be. The two directions cross the atmosphere
along two different paths, so the correction does not fit the uplink path. That
difference is the point-ahead anisoplanatism (Stone, Hu, Mills and Ma,
DOI 10.1364/JOSAA.11.000347). The runner models it with a laterally SHIFTED
WINDOW of the SAME screens, and it reads the uplink through the Shapiro
reciprocity overlap (DOI 10.1364/JOSA.61.000492).

THE FLOW. One scenario, one geometry, one campaign, then the budget:

    campaign = Campaign(scn, orbit, ROOT, seed=..., compensation="terminal",
                        point_ahead_rad=(0.0, theta, 10 arcsec)).run(40)
    uplink_budget(scn, orbit, fidelity=2, wave=campaign)

The campaign IS the fidelity-2 wave record, so it goes straight into the `wave`
slot. The budget reads the record angle that matches the geometry.

THE GRID IS SMALL and the trial count is low, so the script runs on a laptop in
about a minute. The numbers are a DEMONSTRATION, not a result: the study of
record is `validation/waveoptics_pointahead/`. Run the script two times: the
second run computes NOTHING, because the blocks are already on disk.

Run it with:
    python examples/waveoptics/uplink_point_ahead.py
"""

import os
import time

import numpy as np

from olb.geometry import CircularOrbit
from olb.links import uplink_budget
from olb.scenario import Channel, DownlinkBeacon, Site, SpaceScenario
from olb.terminal import AO, SMF, Terminal, TipTilt, Transmitter
from olb.waveoptics.turbulence import Campaign

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    "_campaigns", "uplink_point_ahead")

ARCSEC = np.pi / (180.0 * 3600.0)        # one arcsecond, in rad.
N_TRIALS = 40


def build_scenario():
    """Give the pre-compensated hero uplink and its orbit.

    The ground terminal launches from the FULL aperture, it carries an AO(10)
    stack, and the scenario names a DownlinkBeacon, so the ground stack senses
    the beacon direction.

    Returns:
        The pair (SpaceScenario, CircularOrbit).
    """
    lam = 1550e-9
    site = Site(cn2_ground=1.7e-14, wind_rms_m_s=21.0, outer_scale_m=25.0)
    ground = Terminal(aperture_m=0.40, wavelength_m=lam,
                      pointing_jitter_rad=2e-6,
                      detector=SMF(sensitivity_dbm=-45.0),
                      compensation=[AO(10)],
                      transmitter=Transmitter(waist_m=0.20, power_dbm=30.0))
    space = Terminal(aperture_m=0.10, wavelength_m=lam,
                     pointing_jitter_rad=1e-6)
    scn = SpaceScenario(ground=ground, space=space, direction="uplink",
                        channel=Channel(site=site, altitude_m=500e3),
                        precompensation=DownlinkBeacon())
    geom = CircularOrbit(altitude_m=500e3, elevation_deg=30.0)
    return scn, geom


def main():
    """Run the campaign, print the losses, then make the budget."""
    scn, geom = build_scenario()

    # The angle of the geometry. The budget reads THIS angle in the record.
    theta = float(np.asarray(geom.point_ahead_rad, dtype=float).ravel()[0])

    # THE ANGLE LIST. 0.0 repeats the beacon path, so it is the control column;
    # the geometry angle is the one the budget reads; 10 arcsec is a wide
    # reference. The runner sizes the oversize screen from the widest angle.
    angles = (0.0, theta, 10.0 * ARCSEC)

    # `store_screen_phase=True` keeps the summed screen phase of each trial, so
    # the POST-HOC section below can sense the beacon with another stack.
    campaign = Campaign(scn, geom, ROOT, seed=20260911, preset="rapid",
                        block_size=20, compensation="terminal",
                        point_ahead_rad=angles, store_screen_phase=True)
    t0 = time.time()
    n = campaign.run(N_TRIALS, progress=True)
    print(f"\n{n} trials on disk after {time.time() - t0:.1f} s "
          f"({campaign.grid.n} px grid, {campaign.screen_n} px screens, "
          f"{campaign.plan.z_m.size} screens, "
          f"{campaign.n_modes_corrected} Noll modes corrected)")

    # 1. THE LOSS OF EACH ANGLE, from the stored overlaps. The loss is
    # -10*log10(eta), and p5 is the loss the link exceeds 5 percent of the time.
    result = campaign.load()
    eta = np.array([t.eta_turb_pa for t in result.trials], dtype=float)
    loss = -10.0 * np.log10(eta)
    print("\nuplink turbulence loss by point-ahead angle "
          f"({N_TRIALS} trials, AO(10) pre-compensation):")
    print("  angle [arcsec]      mean [dB]      p5 [dB]")
    for i, a in enumerate(campaign.point_ahead_rad):
        print(f"  {a / ARCSEC:>10.2f}    {loss[:, i].mean():>11.2f}  "
              f"{np.percentile(loss[:, i], 95.0):>11.2f}")
    print("  (the 0 arcsec row is the BEACON direction, the control case)")

    # 2. THE BUDGET. It reads the angle of its own geometry.
    budget = uplink_budget(scn, geom, fidelity=2, wave=campaign)
    turb = next(t for t in budget.terms if t.meta.get("model") == "waveoptics")
    print("\nfidelity-2 uplink budget from the campaign:")
    print(budget.to_frame().to_string(index=False))
    print(f"  total {float(budget.total_loss_db()):.2f} dB")
    print(f"\nthe turbulence Term reads the point-ahead angle "
          f"{turb.meta['point_ahead_arcsec']:.2f} arcsec "
          f"(record column {turb.meta['point_ahead_index']}), "
          f"screen margin {turb.meta['screen_margin_m']:.2f} m, "
          f"screen {turb.meta['screen_n']} px")
    for flag in turb.assumptions.violations:
        print(f"  FLAG: {flag.split(':')[0]}")
    modelled = any("POINT-AHEAD ANISOPLANATISM MODELLED" in v
                   for v in turb.assumptions.violations)
    absent = not any("NO ANISOPLANATISM" in v
                     for v in turb.assumptions.violations)
    print(f"  POINT-AHEAD ANISOPLANATISM MODELLED present: {modelled}")
    print(f"  NO ANISOPLANATISM absent:                    {absent}")

    # 3. POST HOC, with no new propagation. The stored planes answer ANY
    # compensation stack, so a tip-tilt terminal costs no trial.
    eta_tt = campaign.recouple_point_ahead([TipTilt()])
    loss_tt = -10.0 * np.log10(eta_tt)
    print("\npost hoc (no propagation): the SAME planes with a TipTilt stack:")
    for i, a in enumerate(campaign.point_ahead_rad):
        print(f"  {a / ARCSEC:>10.2f} arcsec   mean "
              f"{loss_tt[:, i].mean():>6.2f} dB")


if __name__ == '__main__':
    main()
