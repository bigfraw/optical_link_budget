'''
The fidelity-2 focal spot on a tracking camera: pixels, centroid and jitter.

A tracking camera does not see the continuous focal intensity. It sees the POWER
IN EACH PIXEL. The tracking loop then measures the CENTROID of that pixelated
image, and it drives the fine-steering mirror with it. So the camera image is the
sensor signal, and its frame-to-frame centroid scatter is the tracking jitter the
loop must follow.

The script takes a 600 km downlink at 30 degrees into a 0.7 m ground telescope
with a Camera detector. It runs a handful of turbulent snapshots as a CAMPAIGN
(olb.waveoptics.turbulence.Campaign), it rebuilds the stored receive field of
each trial, it clips that field at the ground aperture, and it bins the focal
spot onto the camera pixels (olb.waveoptics.camera.camera_image). For each
snapshot it prints:
  * the measured centroid, in pixels and in microradians on the sky;
  * the second-moment spot radius, in pixels;
  * the fraction of the collected power that lands on the sensor.
It prints one STILL-atmosphere row too, from the same grid and the same screen
positions with a very large Fried parameter. The still row is the instrument
floor: the centroid is on the axis, and the spot is the diffraction spot. That
row is a second, one-trial campaign.

THE CAMPAIGN IS THE RECORD. Each campaign stores the receive-plane field of
every trial on a disc of the receive-aperture radius. The script rebuilds the
full grid of one stored trial with the same helper the post-hoc coupling uses
(olb.waveoptics.turbulence.run), so the camera image is the image of the run.
The stores go under examples/waveoptics/_campaigns/camera_tracking/, and a
second run of the script computes nothing.

THE PLATE SCALE. A camera measures a position, and the focal length turns that
position into an ANGLE: theta = x/f. So the centroid scatter in microradians is
the arrival-angle jitter of the beam.

The script builds NO budget change and NO Term. It demonstrates the camera
discretisation alone.

The figure goes to examples/waveoptics/figures/camera_tracking.png.

Run from the repo root:
    python -m examples.waveoptics.camera_tracking
'''

import dataclasses
import os
import time
import warnings

import matplotlib.pyplot as plt
import numpy as np

from olb import CircularOrbit, Terminal, Transmitter
from olb.terminal import Camera
from olb.scenario import Channel, SpaceScenario
from olb.turbulence.profiles import DEFAULT_HS, default_cn2_profile
from olb.waveoptics.camera import camera_image, spot_metrics
from olb.waveoptics.run import _clip
from olb.waveoptics.turbulence import Campaign
from olb.waveoptics.turbulence.run import _patch_field, _rebuilt_fields

WAVELENGTH_M = 1550e-9
ALTITUDE_M = 600e3
GROUND_APERTURE_M = 0.70
GROUND_OBSCURATION = 0.0        # a clean pupil, so the still spot is a clean Airy
SPACE_APERTURE_M = 0.10
SPACE_WAIST_M = 0.05
ELEVATION_DEG = 30              # a low elevation, for strong turbulence
PRESET = "rapid"
N_SNAPSHOTS = 5
BLOCK_SIZE = 1                  # one trial for each block, so the pool fills
WORKERS = 4
SEED = 20260902

# The camera. A 15 um pitch on a 256 x 256 array is a common short-wave infrared
# tracking sensor.
PIXEL_PITCH_M = 15e-6
N_PIXELS = 256

# THE SIZING RULE. The still diffraction spot 1/e^2 radius is
#     w_s = lambda*f/(pi*(D/2))
# (Goodman, Introduction to Fourier Optics, ISBN 978-0974707723). The focal
# length is set so that w_s is TWO pixels, so the still spot spans about four
# pixels across. A tracking camera needs several pixels across the spot, because
# the centroid of a one-pixel spot carries no sub-pixel information.

SPOT_RADIUS_PIXELS = 2.0
FOCAL_LENGTH_M = (np.pi * (GROUND_APERTURE_M / 2.0)
                  * SPOT_RADIUS_PIXELS * PIXEL_PITCH_M / WAVELENGTH_M)

# THE GRID RULE. The fine focal sample is dx_focal = lambda*f/size, so the pupil
# grid SIDE sets the focal sampling. camera_image needs about three fine samples
# per camera pixel, and it needs the focal window (N//2)*dx_focal to cover the
# sensor half-side. The automatic turbulent grid is sized for the PUPIL, not for
# the focal plane, so this script enlarges it by GRID_ZOOM at a CONSTANT pixel
# pitch. A larger window with the same pixel keeps the turbulence sampling
# exactly as the sizer planned it, and it refines the focal plane by GRID_ZOOM.
GRID_ZOOM = 3

QUIET_R0_M = 1e6                # a huge Fried parameter makes the screens flat

# The panels show the middle of the sensor only. The spot and its motion are a
# few tens of pixels wide, and the whole 256 x 256 sensor hides them.
WINDOW_PIXELS = 48

# The campaign store. Each case keeps its own directory under this root.
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "_campaigns",
                    "camera_tracking")

PNG = "examples/waveoptics/figures/camera_tracking.png"


def build_scenario():
    '''Build the downlink scenario: a tracking camera on the ground.

    The pointing jitter is zero on both terminals, so the measured centroid
    scatter is the TURBULENCE alone. The receive mechanical jitter is a separate
    analytic Term, and it is not in this demonstration.
    '''
    ground = Terminal(aperture_m=GROUND_APERTURE_M,
                      obscuration_ratio=GROUND_OBSCURATION,
                      wavelength_m=WAVELENGTH_M, pointing_jitter_rad=0.0,
                      detector=Camera(pixel_pitch_m=PIXEL_PITCH_M,
                                      n_pixels=N_PIXELS,
                                      focal_length_m=FOCAL_LENGTH_M))
    space = Terminal(aperture_m=SPACE_APERTURE_M, wavelength_m=WAVELENGTH_M,
                     pointing_jitter_rad=0.0,
                     transmitter=Transmitter(waist_m=SPACE_WAIST_M,
                                             power_dbm=30.0))
    return SpaceScenario(ground=ground, space=space, direction="downlink",
                         channel=Channel(altitude_m=ALTITUDE_M))


def measure(result, row):
    '''Rebuild one STORED snapshot, clip it, and measure the camera image.

    The campaign keeps the receive-plane field of a trial on a disc of the
    receive-aperture radius. The rebuild scatters that disc back on the FULL
    grid, because the focal-plane pixel scale reads the whole grid. It uses the
    helpers of olb.waveoptics.turbulence.run, the same rebuild the post-hoc
    coupling (recouple) uses, so the image is the image of the run.

    Args:
        result: the TurbWaveResult of a campaign, loaded with fields=True.
        row:    the trial index.

    Returns:
        The pair (image, metrics). The clip applies the ground aperture, so the
        image holds the light the telescope collects.
    '''
    _, array = next(_rebuilt_fields(result, GROUND_APERTURE_M, [row]))
    field = _patch_field(result.patch, array, WAVELENGTH_M)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        clipped = _clip(field, GROUND_APERTURE_M, GROUND_OBSCURATION)
        image, _ = camera_image(clipped, FOCAL_LENGTH_M, PIXEL_PITCH_M,
                                N_PIXELS)
    return image, spot_metrics(image, PIXEL_PITCH_M)


def draw(img_still, img_turb, m_still, m_turb):
    '''Draw the still and the turbulent pixelated image, on one colour scale.

    Each panel marks the MEASURED centroid. The still spot sits on the axis and
    it is a few pixels wide. The turbulent spot is broader, and its centroid has
    moved: that motion is the tracking signal.
    '''
    # The spot is small against the whole sensor, so the panels show the middle
    # WINDOW_PIXELS x WINDOW_PIXELS pixels only. The metrics use the whole image.
    lo = N_PIXELS // 2 - WINDOW_PIXELS // 2
    hi = lo + WINDOW_PIXELS
    img_still = img_still[lo:hi, lo:hi]
    img_turb = img_turb[lo:hi, lo:hi]
    # The window edge, in um, on the zero-centred pixel grid of camera_image.
    half = (lo - 0.5 * N_PIXELS) * PIXEL_PITCH_M * 1e6
    extent = [half, -half, half, -half]
    vmax = float(max(img_still.max(), img_turb.max()))

    fig, (a0, a1) = plt.subplots(1, 2, figsize=(11.2, 5.3),
                                 constrained_layout=True)
    for ax, img, m, title in (
            (a0, img_still, m_still, "no turbulence (the instrument floor)"),
            (a1, img_turb, m_turb, "one turbulent snapshot")):
        im = ax.imshow(img, extent=extent, origin="lower", cmap="inferno",
                       vmin=0.0, vmax=vmax)
        fig.colorbar(im, ax=ax, shrink=0.85, label="power fraction per pixel")
        ax.plot(m.centroid_x_m * 1e6, m.centroid_y_m * 1e6, "+", color="cyan",
                markersize=14, markeredgewidth=2.0, label="measured centroid")
        ax.set_xlabel("x, um")
        ax.set_ylabel("y, um")
        ax.set_title(f"{title}\ncentroid ({m.centroid_x_m / PIXEL_PITCH_M:+.2f}, "
                     f"{m.centroid_y_m / PIXEL_PITCH_M:+.2f}) px, rms radius "
                     f"{m.rms_radius_m / PIXEL_PITCH_M:.2f} px", fontsize=10)
        ax.legend(loc="upper right", fontsize=8, framealpha=0.7)

    fig.suptitle("Fidelity-2 focal spot on a tracking camera, "
                 f"{ALTITUDE_M * 1e-3:.0f} km downlink, {ELEVATION_DEG:.0f} deg\n"
                 f"ground aperture {GROUND_APERTURE_M * 1e3:.0f} mm, "
                 f"{N_PIXELS} x {N_PIXELS} pixels of "
                 f"{PIXEL_PITCH_M * 1e6:.0f} um, {WAVELENGTH_M * 1e9:.0f} nm",
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
    urad_per_pixel = PIXEL_PITCH_M / FOCAL_LENGTH_M * 1e6

    print("=" * 78)
    print(f"{ALTITUDE_M * 1e-3:.0f} km downlink, the fidelity-2 focal spot on a "
          f"tracking camera")
    print("=" * 78)
    print(f"  wavelength              {WAVELENGTH_M * 1e9:11.1f} nm")
    print(f"  ground aperture         {GROUND_APERTURE_M * 1e3:11.1f} mm, "
          f"obscuration {GROUND_OBSCURATION:g}")
    print(f"  camera                  {N_PIXELS:6d} x {N_PIXELS} pixels of "
          f"{PIXEL_PITCH_M * 1e6:.1f} um")
    print(f"  focal length            {FOCAL_LENGTH_M:11.2f} m")
    print(f"  still spot radius w_s   {w_s * 1e6:11.2f} um "
          f"({w_s / PIXEL_PITCH_M:.2f} pixels)")
    print(f"  plate scale             {urad_per_pixel:11.3f} urad per pixel")
    print(f"  elevation               {ELEVATION_DEG:11.1f} deg")
    print(f"  preset                  {PRESET:>11}")
    print(f"  snapshots               {N_SNAPSHOTS:11d}")
    print(f"  pool workers            {WORKERS:11d}")

    # The turbulent run needs the phase screens. The default "olb" generator is
    # self-contained, so this normally runs. Skip gracefully if an import fails.
    try:
        from olb.waveoptics.turbulence import turbulent_grid
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            grid, plan, report = turbulent_grid(
                scenario, orbit, preset=PRESET, hs=hs, cn2_profile=cn2_profile)

        # Enlarge the grid at a CONSTANT pixel pitch, so the focal plane is well
        # sampled for the camera. See THE GRID RULE at the top of the module.
        grid = dataclasses.replace(grid, n=grid.n * GRID_ZOOM,
                                   size_m=grid.size_m * GRID_ZOOM)
        dx_focal = WAVELENGTH_M * FOCAL_LENGTH_M / grid.size_m
        print("")
        print(f"  grid {grid.n} px, {grid.size_m:.3f} m ({GRID_ZOOM}x the sizer "
              f"grid, same pixel); {plan.z_m.size} screens")
        print(f"  r0 = {plan.r0_total_m * 1e2:.2f} cm; "
              f"D/r0 = {GROUND_APERTURE_M / plan.r0_total_m:.2f}")
        print(f"  fine focal sample {dx_focal * 1e6:.2f} um "
              f"({PIXEL_PITCH_M / dx_focal:.2f} per camera pixel; 3 is the floor)")

        quiet_plan = dataclasses.replace(
            plan, r0_m=np.full_like(plan.r0_m, QUIET_R0_M))

        # THE TURBULENT CAMPAIGN. It stores the receive field of every trial,
        # so a second run of the script computes nothing.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            campaign = Campaign(scenario, orbit,
                                os.path.join(ROOT, "turbulent"), seed=SEED,
                                preset=PRESET, block_size=BLOCK_SIZE,
                                grid=grid, plan=plan, hs=hs,
                                cn2_profile=cn2_profile)
            campaign.run(N_SNAPSHOTS, workers=WORKERS)
            turb = campaign.load(N_SNAPSHOTS)

            # The still atmosphere: the same grid and the same screen
            # positions, with a very large Fried parameter, so the screens are
            # flat. One trial is enough, so this campaign holds one trial.
            still = Campaign(scenario, orbit, os.path.join(ROOT, "still"),
                             seed=1, preset=PRESET, block_size=1, grid=grid,
                             plan=quiet_plan, hs=hs, cn2_profile=cn2_profile)
            still.run(1)
        print(f"  campaign {campaign.root_dir}")

        img_still, m_still = measure(still.load(1), 0)

        images, metrics = [], []
        for trial in range(N_SNAPSHOTS):
            img, m = measure(turb, trial)
            images.append(img)
            metrics.append(m)

        print("")
        print("snapshot   centroid x   centroid y   centroid x   centroid y   "
              "rms radius   on sensor")
        print("             pixels       pixels        urad         urad       "
              "  pixels               ")
        rows = [("still   ", m_still)] + [(f"turb {i}  ", m)
                                          for i, m in enumerate(metrics)]
        for name, m in rows:
            # The plate scale: theta = x/f. See the module docstring.
            tx = m.centroid_x_m / FOCAL_LENGTH_M * 1e6
            ty = m.centroid_y_m / FOCAL_LENGTH_M * 1e6
            print(f"{name} {m.centroid_x_m / PIXEL_PITCH_M:11.3f}  "
                  f"{m.centroid_y_m / PIXEL_PITCH_M:11.3f}  {tx:11.3f}  "
                  f"{ty:11.3f}  {m.rms_radius_m / PIXEL_PITCH_M:11.3f}  "
                  f"{m.on_sensor_fraction:10.4f}")

        cx = np.array([m.centroid_x_m for m in metrics]) / FOCAL_LENGTH_M * 1e6
        cy = np.array([m.centroid_y_m for m in metrics]) / FOCAL_LENGTH_M * 1e6
        print("")
        print(f"  centroid scatter (std over {N_SNAPSHOTS} snapshots): "
              f"x {cx.std(ddof=1):.3f} urad, y {cy.std(ddof=1):.3f} urad")
        print("  That scatter is the arrival-angle jitter the camera measures, "
              "and the")
        print("  tracking loop must follow it. The still row is the instrument "
              "floor.")

        draw(img_still, images[0], m_still, metrics[0])
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
