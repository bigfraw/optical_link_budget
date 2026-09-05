'''
The fidelity-2 multimode-fibre coupling, and the focused spot on the core.

A multimode fibre is a LIGHT BUCKET. The core is a hard disk in the fibre plane,
fixed on the axis. The fibre collects the power of the focused spot that lands
inside that disk. So the coupling is the ENCIRCLED ENERGY of the focal spot
inside the core. The turbulent field carries the instantaneous tilt, so the
focused spot walks off the on-axis core on its own. The fade is therefore
intrinsic to the field, and it needs no added tilt term.

The script does two things on a 600 km downlink into a ground light-bucket
receiver.

1. THE FIDELITY-2 TERM. It runs a turbulent Monte Carlo at 30 degrees (strong
   turbulence) as a CAMPAIGN, and it builds the fidelity-2 MMF coupling Term
   (olb.models.waveoptics.waveoptics_mmf_coupling_term) from the stored trials.
   It prints the three Term faces: the mean coupling loss, the 99 percent fade,
   and the static encircled-energy floor (one still-atmosphere trial). This is
   the only statistical MMF coupling model in olb; there is no analytic or FAST
   sibling.

2. THE PSF IMAGES. It focuses one turbulent snapshot and one no-turbulence
   snapshot on the SAME grid, and it draws the focal-plane intensity on the core
   for each. A drawn circle marks the core edge. The no-turbulence spot sits
   mostly inside the core (high capture). The turbulent spot broadens and it
   spills OUTSIDE the core (a real loss). The panel titles give the power
   fraction inside the core, so the picture and the numbers agree.

WHAT DRIVES THE CONTRAST. The turbulence broadens the focused spot by about the
seeing ratio D/r0. This case has a large D/r0 (a 700 mm aperture at 30 degrees),
so the turbulent spot is several diffraction spots wide. The focal length sets
the still diffraction spot to CORE_RADIUS_M_SMF (5.2 um), well inside the 50 um
core, so the still spot fits and the broadened turbulent spot spills the core. A
weaker path (a higher elevation or a smaller aperture) shrinks D/r0 and it hides
the contrast.

THE CAMPAIGN. Both Monte Carlo runs are an olb.waveoptics.turbulence.Campaign,
the on-disk store of trials. The turbulent campaign holds the snapshots of the
Term. A second campaign of ONE trial on FLAT screens holds the static floor: it
is a one-trial store, but it keeps the same entry point and the second run of
the script reads it from disk. The stores go under
examples/waveoptics/_campaigns/mmf_core_psf/. So a second run of the script
computes no trial. The two PICTURES stay on propagate_turbulent_field, the
single-snapshot diagnostic entry point, because a picture needs the wide field
and the campaign keeps the receive disc only.

The figure goes to examples/waveoptics/figures/mmf_core_psf.png.

The script builds NO budget change. It demonstrates the Term and the shared
focal-plane helper (olb.waveoptics.mmf.focal_intensity).

Run from the repo root:
    python -m examples.waveoptics.mmf_core_psf
'''

import dataclasses
import os
import time
import warnings

import matplotlib.pyplot as plt
import numpy as np

from olb import CircularOrbit, Terminal, Transmitter
from olb.terminal import MMF
from olb.models.waveoptics import waveoptics_mmf_coupling_term
from olb.scenario import Channel, SpaceScenario
from olb.turbulence.plane_wave_scintillation import \
    aperture_averaged_scintillation_index
from olb.turbulence.profiles import DEFAULT_HS, default_cn2_profile
from olb.waveoptics.mmf import focal_intensity
from olb.waveoptics.turbulence import Campaign
from olb.waveoptics.run import _clip

WAVELENGTH_M = 1550e-9
ALTITUDE_M = 600e3
GROUND_APERTURE_M = 0.70
GROUND_OBSCURATION = 0.0        # a clean pupil, so the still spot is a clean Airy
SPACE_APERTURE_M = 0.10
SPACE_WAIST_M = 0.05
ELEVATION_DEG = 30            # a low elevation, for strong turbulence
PRESET = "rapid"
N_TRIALS = 60
BLOCK_SIZE = 15                 # 60 trials in four blocks, one for each worker
WORKERS = 4
SEED = 20260828
FIBRE_NA = 0.2                  # a common step-index MMF; it does not gate here

# The core radius, and the focal length. The still diffraction spot 1/e^2 radius
# is w_s = lambda*f/(pi*(D/2)). The focal length is set so the spot radius is
# CORE_RADIUS_M_SMF (a tight, single-mode-sized 5.2 um spot), well inside the
# 50 um core. So the still spot sits deep in the core, and the focal plane stays
# well sampled. Source: Goodman, Introduction to Fourier Optics,
# ISBN 978-0974707723.
CORE_RADIUS_M = 50e-6
CORE_RADIUS_M_SMF = 5.2e-6
FOCAL_LENGTH_M = (np.pi * (GROUND_APERTURE_M / 2.0) * CORE_RADIUS_M_SMF
                  / (WAVELENGTH_M))

# The Cn2 of a still atmosphere, for the static encircled-energy floor. 1e-24 is
# four decades below any real site, so one trial gives the still-atmosphere spot.
QUIET_CN2 = 1e-24
QUIET_R0_M = 1e6               # a huge Fried parameter makes the screens flat

# The window of the focal-plane pictures, a few core radii wide.
WINDOW_M = 3.0 * CORE_RADIUS_M

# The campaign store. Each case keeps its own directory under this root.
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "_campaigns",
                    "mmf_core_psf")

PNG = "examples/waveoptics/figures/mmf_core_psf.png"


def build_scenario():
    '''Build the downlink scenario: a light-bucket receiver on the ground.

    The pointing jitter is zero on both terminals, so the picture carries the
    TURBULENCE only. The receive mechanical jitter is a separate analytic Term,
    and it is not in this demonstration.
    '''
    ground = Terminal(aperture_m=GROUND_APERTURE_M,
                      obscuration_ratio=GROUND_OBSCURATION,
                      wavelength_m=WAVELENGTH_M, pointing_jitter_rad=0.0,
                      detector=MMF(core_radius_m=CORE_RADIUS_M,
                                   focal_length_m=FOCAL_LENGTH_M,
                                   numerical_aperture=FIBRE_NA,
                                   sensitivity_dbm=-110.0))
    space = Terminal(aperture_m=SPACE_APERTURE_M, wavelength_m=WAVELENGTH_M,
                     pointing_jitter_rad=0.0,
                     transmitter=Transmitter(waist_m=SPACE_WAIST_M,
                                             power_dbm=30.0))
    return SpaceScenario(ground=ground, space=space, direction="downlink",
                         channel=Channel(altitude_m=ALTITUDE_M))


def quiet_plan_of(plan):
    '''Give the same screen plan with a very large Fried parameter.

    The screens are then flat, so the trial carries no turbulence.
    '''
    return dataclasses.replace(plan, r0_m=np.full_like(plan.r0_m, QUIET_R0_M))


def static_floor_db(scenario, orbit, grid, plan, hs, cn2_profile):
    '''Give the still-atmosphere MMF coupling loss, in dB.

    One trial on FLAT screens gives the STATIC encircled-energy floor: the
    fraction of the still diffraction spot that lands inside the core. The
    turbulence part of any other number is that number minus this floor. The
    trial uses the SAME grid and the SAME screen positions, with the Fried
    parameter set very large.

    The trial is a one-trial Campaign, so it is stored on disk like the
    turbulent trials and a second run of the script reads it back.
    '''
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        still = Campaign(scenario, orbit, os.path.join(ROOT, "still"), seed=1,
                         preset=PRESET, block_size=1, grid=grid,
                         plan=quiet_plan_of(plan), hs=hs,
                         cn2_profile=cn2_profile)
        still.run(1)
        result = still.load(fields=False)
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

    The two panels share one colour scale, so the reader sees the still spot
    concentrated inside the core and the turbulent spot spread past it. The
    dashed circle is the core edge. The panel title gives the power fraction
    inside the core. See Goodman, Introduction to Fourier Optics,
    ISBN 978-0974707723, for the focal-plane intensity of a pupil field.
    '''
    n = If_turb.shape[0]
    c = n // 2
    half_px = int(round(WINDOW_M / dx_focal))
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
                 f"{ALTITUDE_M * 1e-3:.0f} km downlink, {ELEVATION_DEG:.0f} deg\n"
                 f"ground aperture {GROUND_APERTURE_M * 1e3:.0f} mm, core radius "
                 f"{CORE_RADIUS_M * 1e6:.0f} um, {WAVELENGTH_M * 1e9:.0f} nm",
                 fontsize=12)
    fig.savefig(PNG, dpi=150)
    plt.close(fig)


def main():
    t_start = time.time()
    scenario = build_scenario()
    hs = DEFAULT_HS
    cn2_profile = default_cn2_profile(scenario.channel.site, hs)
    orbit = CircularOrbit(ALTITUDE_M, elevation_deg=[ELEVATION_DEG])

    w_s = WAVELENGTH_M * FOCAL_LENGTH_M / (np.pi * GROUND_APERTURE_M / 2.0)

    print("=" * 78)
    print(f"{ALTITUDE_M * 1e-3:.0f} km downlink, fidelity-2 multimode-fibre "
          f"coupling and the focal spot")
    print("=" * 78)
    print(f"  wavelength              {WAVELENGTH_M * 1e9:11.1f} nm")
    print(f"  ground aperture         {GROUND_APERTURE_M * 1e3:11.1f} mm, "
          f"obscuration {GROUND_OBSCURATION:g}, MMF light bucket")
    print(f"  core radius             {CORE_RADIUS_M * 1e6:11.1f} um")
    print(f"  focal length            {FOCAL_LENGTH_M:11.2f} m "
          f"(spot radius {w_s * 1e6:.1f} um")
    print(f"  fibre NA                {FIBRE_NA:11.2f}")
    print(f"  elevation               {ELEVATION_DEG:11.1f} deg")
    print(f"  preset                  {PRESET:>11}")
    print(f"  snapshots               {N_TRIALS:11d}")
    print(f"  pool workers            {WORKERS:11d}")

    # The whole run needs the phase screens. Skip gracefully when an import
    # fails (the "olb" generator is self-contained, so this normally runs).
    try:
        from olb.waveoptics.turbulence import propagate_turbulent_field

        # THE CAMPAIGN. It sizes the grid and it plans the screens one time, it
        # keeps the manifest, and it stores every trial. A second run of the
        # script finds the blocks on disk and it computes nothing.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            campaign = Campaign(scenario, orbit,
                                os.path.join(ROOT, "turbulent"), seed=SEED,
                                preset=PRESET, block_size=BLOCK_SIZE, hs=hs,
                                cn2_profile=cn2_profile)
            n_done = campaign.run(N_TRIALS, workers=WORKERS, progress=True)
        grid, plan = campaign.grid, campaign.plan
        print("")
        print(f"  campaign {campaign.root_dir}")
        print(f"  {n_done} trials on disk")
        print(f"  grid {grid.n} px, {grid.size_m:.3f} m; {plan.z_m.size} "
              f"screens; r0 = {plan.r0_total_m * 1e2:.2f} cm; "
              f"D/r0 = {GROUND_APERTURE_M / plan.r0_total_m:.2f}")

        floor_db = static_floor_db(scenario, orbit, grid, plan, hs,
                                   cn2_profile)
        print(f"  static encircled-energy floor {floor_db:.3f} dB "
              f"(one still-atmosphere trial)")

        # The Term reads the stored record. The fields stay on disk, because
        # the Term needs the mmf_eta column only.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            result = campaign.load(N_TRIALS, fields=False)

            # The plane-wave scintillation index, for the Term regime flag.
            sigma2 = float(aperture_averaged_scintillation_index(
                GROUND_APERTURE_M, ELEVATION_DEG, WAVELENGTH_M, hs, cn2_profile))
            term = waveoptics_mmf_coupling_term(result, sigma2_I=sigma2)

        print("")
        print("the fidelity-2 MMF coupling Term:")
        print(f"  category                {term.category:>11}")
        print(f"  mean coupling loss      {term.mean_db:11.3f} dB")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")     # the deep tail may be under-sampled
            q99 = term.quantile_db(0.99)
        print(f"  99% fade                {q99:11.3f} dB")
        print(f"  static floor            {floor_db:11.3f} dB")
        print(f"  turbulence part (mean)  {term.mean_db - floor_db:11.3f} dB "
              f"(mean minus floor)")

        # ---- the focal-plane pictures ----
        # Trial 0 of the campaign, the turbulent snapshot. Then the SAME grid
        # and the SAME screen positions with flat screens (a huge Fried
        # parameter), so the two spots are directly comparable. A picture needs
        # the wide field, so the two snapshots come from the single-snapshot
        # entry point and not from the stored patch.
        print("")
        print("the focused spot on the core:")
        quiet_plan = quiet_plan_of(plan)
        F_turb, _, _ = propagate_turbulent_field(
            scenario, orbit, seed=SEED, trial=0, grid=grid, plan=plan)
        F_still, _, _ = propagate_turbulent_field(
            scenario, orbit, seed=SEED, trial=0, grid=grid, plan=quiet_plan)

        clip_turb = _clip(F_turb, GROUND_APERTURE_M, GROUND_OBSCURATION)
        clip_still = _clip(F_still, GROUND_APERTURE_M, GROUND_OBSCURATION)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            If_turb, dx_focal, eta_turb = core_eta(clip_turb, "turbulent  ")
            If_still, _, eta_still = core_eta(clip_still, "no turbulence")

        print(f"  implied core capture: no turbulence {eta_still:.3f} "
              f"({-10 * np.log10(eta_still):.2f} dB), turbulent {eta_turb:.3f} "
              f"({-10 * np.log10(eta_turb):.2f} dB)")
        print("  The no-turbulence eta is the static encircled energy. The "
              "turbulent")
        print("  spot broadens by about D/r0 and it spills past the core, so it "
              "reads")
        print("  the larger loss. The two panels share one colour scale.")

        draw(If_turb, If_still, dx_focal, eta_turb, eta_still)
        print("")
        print(f"figure saved: {PNG}")
        print(f"(elapsed {time.time() - t_start:.1f} s)")

    except ImportError as exc:
        print("")
        print(f"NOTE: a screen-generator import failed "
              f"({exc.__class__.__name__}: {exc}). The turbulent run needs the "
              "phase screens.")
        print("Skipping.")


if __name__ == "__main__":
    main()
