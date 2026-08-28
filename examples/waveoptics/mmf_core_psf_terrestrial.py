'''
The focused spot on a multimode-fibre core, on a SHORT terrestrial link.

This is the terrestrial sibling of examples/waveoptics/mmf_core_psf.py. Both
scripts draw the focused spot on the light-bucket core, turbulent against no
turbulence. The two links sit in OPPOSITE corners of the turbulence.

THE POINT OF THIS SCRIPT is one PICTURE, not a statistics run. It focuses ONE
turbulent snapshot and ONE still (no-turbulence) snapshot on the SAME grid, and
it draws each focal spot on the core. It runs a single trial of each. It does
NOT run a Monte Carlo, so it gives NO fade distribution.

THE REGIME. A 5 km horizontal link at Cn2 = 5e-15 m^(-2/3), 1550 nm, into a
small 25 mm receiver.

- The plane-wave Fried parameter r0 is about 4.5 cm. The 25 mm aperter is
  SMALLER than one coherence cell, so D/r0 is about 0.55. The wavefront across
  the pupil carries only a mild phase aberration.
- The plane-wave Rytov variance sigma_R^2 is about 1.9. The scintillation is
  STRONG.

So the two effects point two ways. Because D/r0 < 1, the focused spot stays
fairly compact and it mostly sits inside the big 50 um-core bucket. The core is
forgiving. The dominant turbulence loss on this link is therefore the
SCINTILLATION of the collected power and the beam wander, NOT the spot breaking
up and spilling past the core. This is the OPPOSITE corner from the big-aperture
space downlink of mmf_core_psf.py, where D/r0 is about 4 and the broadened spot
spills the core.

HONESTY. The script prints the two measured core-capture fractions (turbulent
against no turbulence). If they are close, that is the physics result of this
corner, and the script says so. It does NOT shrink the core to manufacture a
contrast. If the strong scintillation puts amplitude holes or a displaced,
haloed spot in the turbulent snapshot, the picture shows it faithfully.

The single-trial Term is a light touch: it proves the fidelity-2 MMF coupling
Term is buildable for a terrestrial MMF. It prints the single-snapshot coupling
loss and the still-atmosphere static floor only. A fade distribution needs many
trials, and that is not the purpose here.

The figure goes to examples/waveoptics/figures/mmf_core_psf_terrestrial.png.

The script builds NO budget change. It uses the shared focal-plane helper
olb.waveoptics.mmf.focal_intensity, so it duplicates no FFT.

Run from the repo root:
    python -m examples.waveoptics.mmf_core_psf_terrestrial
'''

import dataclasses
import time
import warnings

import matplotlib.pyplot as plt
import numpy as np

from olb import Terminal, Transmitter
from olb.geometry import HorizontalPath
from olb.models.waveoptics import waveoptics_mmf_coupling_term
from olb.scenario import TerrestrialChannel, TerrestrialScenario
from olb.terminal import MMF
from olb.turbulence.gaussian_fried import plane_wave_fried_parameter
from olb.waveoptics.mmf import focal_intensity
from olb.waveoptics.run import _clip

WAVELENGTH_M = 1550e-9
PATH_M = 5000.0
CN2 = 5e-15                     # r0 about 4.5 cm, sigma_R^2 about 1.9

RX_APERTURE_M = 0.025          # a small 25 mm receiver, smaller than r0
RX_OBSCURATION = 0.0           # a clean pupil, so the still spot is a clean Airy

# A clean launch that OVERFILLS the small receiver. A 25 mm waist gives about a
# 100 mm beam at 5 km, so the beam overfills the 25 mm receiver about 4 times.
# The receiver then sees a near-flat wavefront across the pupil, so the picture
# carries the TURBULENCE, not the launch shape.
WAIST_M = 0.025
TX_APERTURE_M = 0.10           # 4 waists: a negligible launch clip

# The multimode fibre, a light bucket. A standard 50 um-core MMF (25 um radius).
# The focal length matches a 5.2 um diffraction spot radius to the core: the
# still spot 1/e^2 radius is w_s = lambda*f/(pi*(D/2)), so f gives w_s = 5.2 um
# well inside the 25 um core. Source: Goodman, Introduction to Fourier Optics,
# ISBN 978-0974707723.
CORE_RADIUS_M = 25e-6
CORE_RADIUS_M_SMF = 5.2e-6
FOCAL_LENGTH_M = (np.pi * (RX_APERTURE_M / 2.0) * CORE_RADIUS_M_SMF
                  / WAVELENGTH_M)
FIBRE_NA = 0.2                 # a common step-index MMF; it does not gate here

PRESET = "rapid"               # a spot picture, not a statistics run
N_TRIALS = 1
SEED = 20260828
QUIET_R0_M = 1e6               # a huge Fried parameter makes the screens flat

# The window of the focal-plane pictures, a few core radii wide.
WINDOW_M = 3.0 * CORE_RADIUS_M

PNG = "examples/waveoptics/figures/mmf_core_psf_terrestrial.png"


def build_scenario():
    '''Build the horizontal scenario: a clean launch into a light-bucket MMF.

    The transmitter is on the NEAR terminal, and the MMF is on the FAR terminal.
    The pointing jitter is zero on both terminals, so the picture carries the
    TURBULENCE only. The receive mechanical jitter is a separate analytic Term,
    and it is not in this demonstration. The path attenuation is zero, because
    the extinction loss does not change the focal spot.
    '''
    near = Terminal(aperture_m=TX_APERTURE_M, wavelength_m=WAVELENGTH_M,
                    pointing_jitter_rad=0.0,
                    transmitter=Transmitter(waist_m=WAIST_M, power_dbm=30.0))
    far = Terminal(aperture_m=RX_APERTURE_M, obscuration_ratio=RX_OBSCURATION,
                   wavelength_m=WAVELENGTH_M, pointing_jitter_rad=0.0,
                   detector=MMF(core_radius_m=CORE_RADIUS_M,
                                focal_length_m=FOCAL_LENGTH_M,
                                numerical_aperture=FIBRE_NA,
                                sensitivity_dbm=-110.0))
    return TerrestrialScenario(
        near=near, far=far,
        channel=TerrestrialChannel(path_length_m=PATH_M,
                                   attenuation_db_per_km=0.0, cn2=CN2))


def static_floor_db(scenario, geometry, grid, plan):
    '''Give the still-atmosphere MMF coupling loss, in dB.

    One trial on FLAT screens gives the STATIC encircled-energy floor: the
    fraction of the still diffraction spot that lands inside the core. The
    trial uses the SAME grid and the SAME screen positions, with the Fried
    parameter set very large, so the screens are flat.
    '''
    from olb.waveoptics.turbulence import propagate_turbulent_scenario
    quiet_plan = dataclasses.replace(
        plan, r0_m=np.full_like(plan.r0_m, QUIET_R0_M))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        result = propagate_turbulent_scenario(
            scenario, geometry, n_trials=1, seed=1, preset=PRESET,
            grid=grid, plan=quiet_plan)
    return -10.0 * np.log10(result.trials[0].mmf_eta)


def core_eta(field, path_desc):
    '''Focus one clipped field and give (If, dx_focal, eta) for the core disk.

    The field is already clipped to the receive aperture. focal_intensity
    focuses it and applies the numerical-aperture gate (the shared helper of
    olb.waveoptics.mmf). eta is the power fraction inside the on-axis core disk,
    the same encircled-energy definition the Term uses.
    '''
    If, dx_focal = focal_intensity(field, FOCAL_LENGTH_M,
                                   numerical_aperture=FIBRE_NA)
    n = If.shape[0]
    c = n // 2
    idx = np.arange(n) - c
    xx, yy = np.meshgrid(idx, idx)
    r2 = (xx ** 2 + yy ** 2) * dx_focal ** 2
    eta = float(If[r2 <= CORE_RADIUS_M ** 2].sum() / If.sum())
    n_core = CORE_RADIUS_M / dx_focal
    print(f"    {path_desc}: focal pixel {dx_focal * 1e6:.2f} um, "
          f"{n_core:.1f} pixels per core radius, core capture {eta:.3f}")
    return If, dx_focal, eta


def draw(If_turb, If_still, dx_focal, eta_turb, eta_still):
    '''Draw the focal-plane intensity on the core, turbulent and still.

    The two panels share one colour scale. In this corner (D/r0 < 1) the two
    spots look ALIKE, and both mostly sit inside the big core. The dashed circle
    is the core edge. The panel title gives the power fraction inside the core.
    See Goodman, Introduction to Fourier Optics, ISBN 978-0974707723, for the
    focal-plane intensity of a pupil field.
    '''
    n = If_turb.shape[0]
    c = n // 2
    # The focal field of view is limited. A small receive aperture in the
    # beam-sized grid gives a SHORT focal length, so the focal FOV half-width
    # (c*dx_focal) is small here, about +/-52 um. The window must stay INSIDE
    # the FOV, clear of the periodic FFT replica that sits at the FOV edge. So
    # the code clamps the half window to 0.72 of the half FOV. A window past the
    # FOV would wrap and show the replica, not the on-axis spot.
    half_px = min(int(round(WINDOW_M / dx_focal)), int(0.72 * c))
    lo, hi = c - half_px, c + half_px + 1
    win_um = half_px * dx_focal * 1e6
    extent = [-win_um, win_um, -win_um, win_um]
    core_um = CORE_RADIUS_M * 1e6

    turb = If_turb[lo:hi, lo:hi]
    still = If_still[lo:hi, lo:hi]
    vmax = float(max(turb.max(), still.max()))

    fig, (a0, a1) = plt.subplots(1, 2, figsize=(11.2, 5.3),
                                 constrained_layout=True)
    for ax, img, eta, title in (
            (a0, still, eta_still, "no turbulence (the static floor)"),
            (a1, turb, eta_turb, "one turbulent snapshot")):
        im = ax.imshow(img, extent=extent, origin="lower", cmap="inferno",
                       vmin=0.0, vmax=vmax)
        fig.colorbar(im, ax=ax, shrink=0.85, label="focal intensity, arb.")
        ax.add_patch(plt.Circle((0, 0), core_um, fill=False, edgecolor="cyan",
                                 linewidth=1.8, linestyle="--",
                                 label="core edge"))
        ax.set_xlabel("x, um")
        ax.set_ylabel("y, um")
        ax.set_title(f"{title}\ncore capture eta = {eta:.3f} "
                     f"({-10 * np.log10(eta):.2f} dB)", fontsize=10)
        ax.legend(loc="upper right", fontsize=8, framealpha=0.7)

    fig.suptitle("Fidelity-2 focused spot on a multimode-fibre core, "
                 f"{PATH_M * 1e-3:.0f} km horizontal link\n"
                 f"receive aperture {RX_APERTURE_M * 1e3:.0f} mm, core radius "
                 f"{CORE_RADIUS_M * 1e6:.0f} um, Cn2 = {CN2:g} m^(-2/3), "
                 f"{WAVELENGTH_M * 1e9:.0f} nm", fontsize=12)
    fig.savefig(PNG, dpi=150)
    plt.close(fig)


def main():
    t_start = time.time()
    scenario = build_scenario()
    geometry = HorizontalPath(PATH_M)

    w_s = WAVELENGTH_M * FOCAL_LENGTH_M / (np.pi * RX_APERTURE_M / 2.0)

    # The regime, computed before the run. The plane-wave Fried parameter r0
    # (Fried, DOI 10.1364/JOSA.56.001372), the seeing ratio D/r0, and the
    # plane-wave Rytov variance sigma_R^2 = 1.23 Cn2 k^(7/6) L^(11/6), with
    # k = 2*pi/lambda (Andrews and Phillips, 2nd ed. 2005,
    # DOI 10.1117/3.626196, Ch. 8).
    r0 = float(plane_wave_fried_parameter(PATH_M, CN2, WAVELENGTH_M))
    d_over_r0 = RX_APERTURE_M / r0
    k = 2.0 * np.pi / WAVELENGTH_M
    sigma2_R = 1.23 * CN2 * k ** (7.0 / 6.0) * PATH_M ** (11.0 / 6.0)

    print("=" * 78)
    print(f"{PATH_M * 1e-3:.0f} km horizontal link, fidelity-2 multimode-fibre "
          f"coupling and the focal spot")
    print("=" * 78)
    print(f"  wavelength              {WAVELENGTH_M * 1e9:11.1f} nm")
    print(f"  path length             {PATH_M * 1e-3:11.3f} km")
    print(f"  Cn2                     {CN2:11.3g} m^(-2/3)")
    print(f"  transmit waist          {WAIST_M * 1e3:11.1f} mm")
    print(f"  launch aperture         {TX_APERTURE_M * 1e3:11.1f} mm "
          f"(alpha = {(TX_APERTURE_M / 2) / WAIST_M:.2f}, a clean launch)")
    print(f"  receive aperture        {RX_APERTURE_M * 1e3:11.1f} mm, "
          f"obscuration {RX_OBSCURATION:g}, MMF light bucket")
    print(f"  core radius             {CORE_RADIUS_M * 1e6:11.1f} um")
    print(f"  focal length            {FOCAL_LENGTH_M:11.4f} m "
          f"(spot radius {w_s * 1e6:.1f} um)")
    print(f"  fibre NA                {FIBRE_NA:11.2f}")
    print(f"  preset                  {PRESET:>11}")
    print(f"  snapshots               {N_TRIALS:11d}")
    print("")
    print("the regime:")
    print(f"  r0, plane wave          {r0 * 1e2:11.2f} cm")
    print(f"  D / r0                  {d_over_r0:11.2f}   "
          f"(smaller than 1: a mild pupil phase)")
    print(f"  sigma_R^2, plane wave   {sigma2_R:11.2f}   "
          f"(larger than 1: strong scintillation)")
    print("  So the loss is FADE-dominated, not a spot that spills the core. The")
    print("  small aperture is inside one coherence cell, so the focused spot")
    print("  stays compact and the big light bucket holds it.")

    # The whole run needs aotools for the phase screens. Skip gracefully when it
    # is absent (the optional "screens" extra).
    try:
        from olb.waveoptics.turbulence import (propagate_turbulent_field,
                                               propagate_turbulent_scenario,
                                               turbulent_grid)
        with warnings.catch_warnings(record=True):
            warnings.simplefilter("always")
            grid, plan, report = turbulent_grid(scenario, geometry,
                                                preset=PRESET)
        print("")
        print(f"  grid {grid.n} px, {grid.size_m:.3f} m; {plan.z_m.size} "
              f"screens; r0(plan) = {plan.r0_total_m * 1e2:.2f} cm")

        floor_db = static_floor_db(scenario, geometry, grid, plan)
        print(f"  static encircled-energy floor {floor_db:.3f} dB "
              f"(one still-atmosphere trial)")

        # ---- the single-trial Term (a light touch) ----
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            result = propagate_turbulent_scenario(
                scenario, geometry, n_trials=N_TRIALS, seed=SEED, preset=PRESET,
                grid=grid, plan=plan)
            term = waveoptics_mmf_coupling_term(result, sigma2_I=sigma2_R)

        print("")
        print("the fidelity-2 MMF coupling Term (one snapshot):")
        print(f"  category                {term.category:>11}")
        print(f"  snapshot coupling loss  {term.mean_db:11.3f} dB")
        print(f"  static floor            {floor_db:11.3f} dB")
        print(f"  turbulence part         {term.mean_db - floor_db:11.3f} dB "
              f"(snapshot minus floor)")
        print("  This is ONE snapshot, so it is a single sample, not a mean. A")
        print("  fade distribution needs many trials, and that is not the")
        print("  purpose here. No 99% fade is printed.")

        # ---- the focal-plane pictures (the point of the script) ----
        # Trial 0 of the run, the turbulent snapshot. Then the SAME grid and the
        # SAME screen positions with flat screens (a huge Fried parameter), so
        # the two spots are directly comparable. See static_floor_db.
        print("")
        print("the focused spot on the core:")
        quiet_plan = dataclasses.replace(
            plan, r0_m=np.full_like(plan.r0_m, QUIET_R0_M))
        F_turb, _, _ = propagate_turbulent_field(
            scenario, geometry, seed=SEED, trial=0, grid=grid, plan=plan)
        F_still, _, _ = propagate_turbulent_field(
            scenario, geometry, seed=SEED, trial=0, grid=grid, plan=quiet_plan)

        clip_turb = _clip(F_turb, RX_APERTURE_M, RX_OBSCURATION)
        clip_still = _clip(F_still, RX_APERTURE_M, RX_OBSCURATION)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            If_turb, dx_focal, eta_turb = core_eta(clip_turb, "turbulent   ")
            If_still, _, eta_still = core_eta(clip_still, "no turbulence")

        print(f"  implied core capture: no turbulence {eta_still:.3f} "
              f"({-10 * np.log10(eta_still):.2f} dB), turbulent {eta_turb:.3f} "
              f"({-10 * np.log10(eta_turb):.2f} dB)")
        if abs(eta_turb - eta_still) < 0.05:
            print("  The two capture fractions are CLOSE. That is the physics of")
            print("  this corner: D/r0 < 1 keeps the spot compact, so the big")
            print("  light bucket holds the turbulent spot as well as the still")
            print("  spot. The turbulence loss lives in the collected-power fade")
            print("  and the wander, not in the spot spilling the core.")
        else:
            print("  The turbulent spot loses more than the still spot: the")
            print("  strong scintillation puts amplitude holes or a displaced,")
            print("  haloed spot in the snapshot, so some power spills the core.")

        draw(If_turb, If_still, dx_focal, eta_turb, eta_still)
        print("")
        print(f"figure saved: {PNG}")
        print(f"(elapsed {time.time() - t_start:.1f} s)")

    except ImportError as exc:
        print("")
        print(f"NOTE: aotools is absent ({exc.__class__.__name__}). The turbulent "
              "run needs it (the")
        print("optional 'screens' extra). Install aotools to run this example. "
              "Skipping.")


if __name__ == "__main__":
    main()
