'''
The pictures of ONE temporal record: the power against time, and the field.

Gate (e) (`hero_temporal.py`) gives the fade RATE and the fade DURATION as
numbers. This script gives the same record as PICTURES. It computes NO frame:
it only READS the record that is on disk, so it never starts a run.

THE RECORD is the hero downlink at 30 deg, record 0: 4000 frames at 0.5 ms
(2.0 s), the `standard` preset, 512 px at 6.86 mm, with a 0.35 m radius
receive-plane field patch stored for each frame.

WHAT IT DRAWS.
  1. figures/record_timeseries_30.png. Three panels over the full 2 s: the
     bucket power, the fibre-coupled power with its 5 percent level, and a
     100 ms zoom of both around the deepest fibre fade.
  2. figures/record_field_30.gif. 200 consecutive frames (0.1 s) that start 40
     frames before the deepest fibre fade. Each frame holds four panels: the
     wrapped phase on the 0.7 m aperture, the intensity on the aperture, the
     intensity at the fibre tip, and the fibre-coupled power with a marker at
     the current frame.

THE TWO POWERS come from the stored scalars: the bucket power is
`collected_power`, and the fibre-coupled power is `collected_power * smf_eta`.

THE POST-HOC READ RULE (see `docs/api-waveoptics.md`): a PUPIL quantity reads
the compact crop, and a FOCAL-PLANE quantity reads the padded full grid,
because the focal pixel size lambda*f/siz reads the grid EXTENT. This script
reads the FULL grid one time for each frame and it cuts the pupil panels out of
it, so one frame costs one block read, not two.

THE FIBRE OPTIC. The hero ground detector is a bare `SMF()`: it sets no focal
length and `optimal_focus` is False, so there is no focal length to read. The
fibre-tip panel needs one, so the script takes the OPTIMAL-FOCUS rule of the
package, f = pi*(D/2)*w_m/(lambda*1.12) with the SMF-28 mode field radius
(5.2 um). That is the design focal length of the same detector, and it is the
focal length that the `eta_max = 0.8145` of the record assumes. Source:
Shaklan and Roddier, Appl. Opt. 27, 2334 (1988), DOI 10.1364/AO.27.002334.

THE INTENSITY SCALE. The slab of a space link starts from a PLANE WAVE, so the
vacuum limit is a flat intensity I0 over the aperture. Turbulence moves the
power around but it keeps the total, so the MEAN intensity over the aperture
and over the animation window IS I0. The intensity panel divides by that mean,
so 1.0 is the vacuum limit.

Sources:
- Taylor, Proc. R. Soc. Lond. A 164, pp. 476 to 490 (1938),
  DOI 10.1098/rspa.1938.0032. The frozen-flow hypothesis of the time axis.
- Goodman, Introduction to Fourier Optics, ISBN 978-0974707723. The focal-plane
  field is the Fourier transform of the pupil field. See
  olb.waveoptics.mmf.focal_intensity.
- Shaklan and Roddier, DOI 10.1364/AO.27.002334. The single-mode-fibre coupling
  parameter a and the optimal-focus rule.

Run from the repository root:

    python -m validation.temporal_screens.record_plots
'''

import argparse
import dataclasses
import os
import time

import numpy as np

from olb.waveoptics.mmf import focal_intensity
from olb.waveoptics.run import SMF28_MODE_FIELD_RADIUS_M, _smf_focal_length
from validation.temporal_screens.hero_temporal import FIGS, campaign_of

ELEVATION_DEG = 30.0
RECORD = 0
N_ANIM = 200          # the frames of the animation (0.1 s at 0.5 ms).
LEAD_FRAMES = 40      # the frames before the deepest fade that it starts at.
ZOOM_MS = 100.0       # the width of the zoom panel, in ms.
FADE_FRACTION = 0.05  # the 5 percent fade level of gate (e).
FOV_MODE_RADII = 10.0  # the half-width of the fibre-tip panel, in mode radii.
LOG_DECADES = 4.0     # the decades of the fibre-tip log colour scale.
FPS = 20


class _Args:
    '''The campaign settings of the record, as `campaign_of` wants them.

    `hero_temporal.campaign_of` reads a parsed command line. This record was
    written with the DEFAULTS of that script, so the values are fixed here.
    They must match the manifest, or the campaign refuses the root.
    '''

    dt = 5e-4
    frames = 4000
    preset = 'standard'
    block_size = 500
    fft_backend = 'cupy'
    smoke = False
    store_screen_phase = False


def db_rel(power, reference):
    '''Give a power as dB over a reference power.

    Args:
        power:     the power, an array or a float.
        reference: the reference power.

    Returns:
        10*log10(power / reference).
    '''
    return 10.0 * np.log10(np.asarray(power, dtype=float) / float(reference))


def window_start(index, n_total, n_window, lead):
    '''Give the first frame of a window that holds one frame of interest.

    The window starts `lead` frames before the index, and it stays inside the
    record at both ends.

    Args:
        index:    the frame of interest.
        n_total:  the frames on disk.
        n_window: the frames of the window.
        lead:     the wanted frames before the index.

    Returns:
        The first frame of the window, as an int.
    '''
    return int(min(max(0, int(index) - int(lead)), max(0, int(n_total) - int(n_window))))


def fibre_focal_length(scenario):
    '''Give the design focal length of the hero fibre coupler, in m.

    The detector of the record sets no focal length, so the function asks for
    the optimal-focus rule of the package (a = 1.12, the SMF-28 mode field
    radius). Source: Shaklan and Roddier, DOI 10.1364/AO.27.002334.

    Args:
        scenario: the SpaceScenario of the record.

    Returns:
        The focal length in m.
    '''
    detector = scenario.ground.detector
    design = dataclasses.replace(detector, optimal_focus=True)
    return float(_smf_focal_length(design, scenario.ground.aperture_m,
                                   scenario.ground.wavelength_m))


def pupil_cut(field, aperture_m):
    '''Give the square cut of the receive aperture out of a full grid.

    The cut is the same for every frame of a record, so a caller builds it ONE
    time and reuses it. `record_ao_plots.py` uses the same cut, so the two
    scripts draw the same pupil panels.

    Args:
        field:      one frame, as a Field on the full grid.
        aperture_m: the receive aperture DIAMETER, in m.

    Returns:
        The tuple (slice, mask, half_width_m): the slice of both axes, the bool
        aperture mask of the cut, and the half-width of the cut, in m.
    '''
    n = field.field.shape[0]
    half = int(np.ceil(0.5 * aperture_m / field.dx)) + 1
    sl = slice(n // 2 - half, n // 2 + half + 1)
    axis = (np.arange(sl.start, sl.stop) - n // 2) * field.dx
    xx, yy = np.meshgrid(axis, axis)
    return sl, (xx ** 2 + yy ** 2) <= (0.5 * aperture_m) ** 2, half * field.dx


def read_panels(camp, rows, aperture_m, focal_length_m, half_px):
    '''Read the panel arrays of every frame of the animation window.

    The function reads the FULL grid of one frame (a focal-plane quantity needs
    the padded grid), it focuses it to the fibre tip, and it cuts the pupil
    panels out of the same array. So one frame costs ONE block read.

    ponytail: it holds the whole window in memory (about 20 MB at 200 frames).
    Stream it frame by frame if a window ever gets long.

    Args:
        camp:           the Campaign.
        rows:           the frame indices, in order.
        aperture_m:     the receive aperture DIAMETER, in m.
        focal_length_m: the focal length of the fibre coupler, in m.
        half_px:        the half-width of the fibre-tip panel, in focal pixels.

    Returns:
        A dict of the stacked panels and the two pixel sizes.
    '''
    phases, pupils, focals = [], [], []
    mask = None
    for row in rows:
        field = camp.field(int(row), compact=False)
        E = field.field
        if mask is None:
            sl, mask, extent_m = pupil_cut(field, aperture_m)
        cut = E[sl, sl]
        phases.append(np.where(mask, np.angle(cut), np.nan))
        pupils.append(np.where(mask, np.abs(cut) ** 2, np.nan))
        image, dx_focal = focal_intensity(field, focal_length_m)
        mid = image.shape[0] // 2
        focals.append(image[mid - half_px:mid + half_px,
                            mid - half_px:mid + half_px])
    return {'phase': np.array(phases), 'pupil': np.array(pupils),
            'focal': np.array(focals), 'dx_pupil': field.dx,
            'dx_focal': dx_focal, 'extent_m': extent_m}


def draw_timeseries(t_s, bucket, fibre, level, deep, path):
    '''Draw the three-panel power plot of the record. Give the path back.

    Args:
        t_s:    the time of every frame, in s.
        bucket: the bucket power of every frame.
        fibre:  the fibre-coupled power of every frame.
        level:  the 5 percent fibre power level.
        deep:   the frame of the deepest fibre fade.
        path:   the output file.

    Returns:
        The path.
    '''
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    b_db = db_rel(bucket, np.median(bucket))
    f_db = db_rel(fibre, np.median(fibre))
    dt_s = float(t_s[1] - t_s[0])
    half = int(round(0.5 * ZOOM_MS * 1e-3 / dt_s))
    lo = window_start(deep, t_s.size, 2 * half, half)
    sl = slice(lo, lo + 2 * half)

    fig, ax = plt.subplots(3, 1, figsize=(9.0, 8.0))
    ax[0].plot(t_s, b_db, lw=0.4, color='C0')
    ax[0].set_ylabel('bucket / median [dB]')
    ax[0].set_title(f'the hero downlink at {ELEVATION_DEG:.0f} deg, record '
                    f'{RECORD}: {t_s.size} frames at {dt_s * 1e3:.2f} ms')
    ax[1].plot(t_s, f_db, lw=0.4, color='0.3')
    ax[1].axhline(db_rel(level, np.median(fibre)), color='C3', ls='--', lw=1.0,
                  label=f'{FADE_FRACTION * 100:.0f} percent level')
    ax[1].set_ylabel('fibre / median [dB]')
    ax[1].legend(fontsize=8, loc='lower right')
    for a in ax[:2]:
        a.set_xlabel('time [s]')
        a.set_xlim(t_s[0], t_s[-1])
        a.grid(alpha=0.3)
    ax[2].plot(t_s[sl] * 1e3, b_db[sl], lw=0.8, color='C0', label='bucket')
    ax[2].plot(t_s[sl] * 1e3, f_db[sl], lw=0.8, color='0.3', label='fibre')
    ax[2].axvline(t_s[deep] * 1e3, color='C3', lw=0.8,
                  label=f'deepest fade, {f_db[deep]:.1f} dB')
    ax[2].set_xlabel('time [ms]')
    ax[2].set_ylabel('power / median [dB]')
    ax[2].set_title(f'the {ZOOM_MS:.0f} ms around the deepest fibre fade')
    ax[2].grid(alpha=0.3)
    ax[2].legend(fontsize=8, loc='lower right')
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)
    return path


def draw_animation(panels, fibre_db, t_ms, mode_radius_m, path):
    '''Draw the four-panel animation of the window. Give the path back.

    ponytail: the writer is Pillow at a small figure size (about 700 px wide),
    which keeps the file well under 15 MB with no panel downsampling. Move to
    ffmpeg (an mp4) if a longer window is ever wanted.

    Args:
        panels:        the dict of read_panels.
        fibre_db:      the fibre-coupled power of the window, in dB.
        t_ms:          the time of every frame of the window, in ms.
        mode_radius_m: the fibre mode field radius, in m.
        path:          the output file.

    Returns:
        The path.
    '''
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation, PillowWriter
    from matplotlib.colors import LogNorm

    pupil = panels['pupil'] / np.nanmean(panels['pupil'])
    focal = panels['focal']
    ext = panels['extent_m']
    pupil_extent = [-ext, ext, -ext, ext]
    fext = 0.5 * focal.shape[1] * panels['dx_focal'] * 1e6
    focal_extent = [-fext, fext, -fext, fext]
    vmax = float(focal.max())

    fig, ax = plt.subplots(1, 4, figsize=(10.0, 2.8))
    im0 = ax[0].imshow(panels['phase'][0], cmap='twilight', vmin=-np.pi,
                       vmax=np.pi, extent=pupil_extent, origin='lower')
    ax[0].set_title('phase [rad]', fontsize=8)
    im1 = ax[1].imshow(pupil[0], cmap='inferno', vmin=0.0,
                       vmax=float(np.nanpercentile(pupil, 99.9)),
                       extent=pupil_extent, origin='lower')
    ax[1].set_title('intensity / vacuum', fontsize=8)
    im2 = ax[2].imshow(focal[0], cmap='viridis', extent=focal_extent,
                       origin='lower',
                       norm=LogNorm(vmin=vmax / 10.0 ** LOG_DECADES,
                                    vmax=vmax))
    ax[2].add_patch(plt.Circle((0.0, 0.0), mode_radius_m * 1e6, fill=False,
                               color='w', lw=0.8))
    ax[2].set_title('fibre tip [um]', fontsize=8)
    for a in ax[:2]:
        a.set_xlabel('x [m]', fontsize=7)
    ax[2].set_xlabel('x [um]', fontsize=7)
    for a, im in zip(ax[:3], (im0, im1, im2)):
        a.tick_params(labelsize=6)
        fig.colorbar(im, ax=a, fraction=0.046).ax.tick_params(labelsize=6)
    ax[3].plot(t_ms, fibre_db, lw=0.8, color='0.3')
    marker, = ax[3].plot([t_ms[0]], [fibre_db[0]], 'o', color='C3', ms=4)
    ax[3].set_xlabel('time [ms]', fontsize=7)
    ax[3].set_ylabel('fibre / median [dB]', fontsize=7)
    ax[3].tick_params(labelsize=6)
    ax[3].grid(alpha=0.3)
    title = fig.suptitle('', fontsize=9)
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.93))

    def update(k):
        '''Put frame k on every panel.'''
        im0.set_data(panels['phase'][k])
        im1.set_data(pupil[k])
        im2.set_data(focal[k])
        marker.set_data([t_ms[k]], [fibre_db[k]])
        title.set_text(f't = {t_ms[k]:.1f} ms, fibre {fibre_db[k]:.1f} dB')
        return im0, im1, im2, marker, title

    anim = FuncAnimation(fig, update, frames=focal.shape[0], blit=False)
    anim.save(path, writer=PillowWriter(fps=FPS), dpi=70)
    plt.close(fig)
    return path


def _self_check():
    '''Check the two helpers against hand cases. It runs on every call.'''
    assert np.isclose(db_rel(2.0, 1.0), 3.0102999566), db_rel(2.0, 1.0)
    assert np.isclose(db_rel(0.5, 1.0), -3.0102999566)
    # The window keeps the lead when it fits, and it stays inside the record.
    assert window_start(100, 1000, 200, 40) == 60
    assert window_start(10, 1000, 200, 40) == 0
    assert window_start(990, 1000, 200, 40) == 800


def main():
    '''Read the record, draw the two figures, and print the numbers.'''
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument('--frames', type=int, default=N_ANIM,
                        help='the frames of the animation.')
    args = parser.parse_args()

    _self_check()
    t_start = time.time()
    camp = campaign_of(ELEVATION_DEG, RECORD, _Args())
    n = int(camp.n_stored)
    print(f'record {camp.root_dir}')
    print(f'  {n} frames at {camp.dt_s * 1e3:.2f} ms, grid {camp.grid.n} px at '
          f'{camp.grid.pixel_m * 1e3:.2f} mm, patch {camp.patch.radius_m:.2f} m')

    result = camp.load(n, fields=False)
    bucket = np.array([t.collected_power for t in result.trials], dtype=float)
    fibre = bucket * np.array([t.smf_eta for t in result.trials], dtype=float)
    t_s = camp.t_s()
    level = float(np.quantile(fibre, FADE_FRACTION))
    deep = int(np.argmin(fibre))
    depth_db = float(db_rel(fibre[deep], np.median(fibre)))
    print(f'  the 5 percent fibre level is {db_rel(level, np.median(fibre)):.2f}'
          f' dB under the median')
    print(f'  the deepest fibre fade is {depth_db:.2f} dB under the median, at '
          f'frame {deep} (t = {t_s[deep] * 1e3:.1f} ms)')

    os.makedirs(FIGS, exist_ok=True)
    path = draw_timeseries(t_s, bucket, fibre, level, deep,
                           os.path.join(FIGS,
                                        f'record_timeseries_{ELEVATION_DEG:.0f}.png'))
    print(f'  figure saved: {path} ({os.path.getsize(path) / 2 ** 20:.2f} MB)')

    focal_length_m = fibre_focal_length(camp.scenario)
    dx_focal = camp.scenario.ground.wavelength_m * focal_length_m / camp.grid.size_m
    half_px = int(np.ceil(FOV_MODE_RADII * SMF28_MODE_FIELD_RADIUS_M / dx_focal))
    print(f'  the fibre coupler is f = {focal_length_m:.3f} m, so the fibre-tip '
          f'pixel is {dx_focal * 1e6:.2f} um and the panel is '
          f'{2 * half_px} px wide')

    lo = window_start(deep, n, int(args.frames), LEAD_FRAMES)
    rows = np.arange(lo, lo + int(args.frames))
    panels = read_panels(camp, rows, camp.scenario.ground.aperture_m,
                         focal_length_m, half_px)
    path = draw_animation(panels, db_rel(fibre[rows], np.median(fibre)),
                          t_s[rows] * 1e3, SMF28_MODE_FIELD_RADIUS_M,
                          os.path.join(FIGS,
                                       f'record_field_{ELEVATION_DEG:.0f}.gif'))
    print(f'  animation saved: {path} '
          f'({os.path.getsize(path) / 2 ** 20:.2f} MB, {len(rows)} frames at '
          f'{FPS} fps, t = {t_s[rows[0]] * 1e3:.1f} to '
          f'{t_s[rows[-1]] * 1e3:.1f} ms)')
    print(f'(elapsed {time.time() - t_start:.1f} s)')


if __name__ == '__main__':
    main()
