'''
The pictures of ONE temporal record under PERFECT AO: four correction stacks.

`record_plots.py` draws the UNCORRECTED record. This script draws the SAME
record four times: with no correction, with a tip-tilt corrector, with a
10-mode corrector and with a 50-mode corrector. It computes NO frame. It only
READS the record that is on disk and it corrects each stored field post hoc.

THE RECORD is the hero downlink at 30 deg, record 0 (`hero_temporal.py`): 4000
frames at 0.5 ms (2.0 s), the `standard` preset, 512 px at 6.86 mm, with a
0.35 m radius receive-plane field patch stored for each frame.

THE CORRECTION is the post-hoc perfect-AO read of the package,
`Campaign.recouple_compensated`. It removes the first N Noll modes over the
0.7 m receive aperture of each frame, and it then couples the corrected field
into the single-mode fibre. The mode count of a stack is the count rule of
`olb.waveoptics.compensation.modal.modes_from_stack`: a TipTilt stage removes
the first 3 Noll modes, and an AO(n) stage removes the first n. Source:
R. J. Noll, J. Opt. Soc. Am. 66(3), 207 to 211 (1976),
DOI 10.1364/JOSA.66.000207, Table I and Eqs. (2) to (4).

PERFECT AO. The fit is an ideal modal fit of one snapshot. There is no
wavefront-sensor noise, no servo lag and no anisoplanatism. So every corrected
curve is the UPPER BOUND of the benefit of a corrector of that mode count.

THE SENSING SOURCE. A SPACE link senses the SUMMED SCREEN PHASE, because the
slab starts from a plane wave. The record must then hold that phase
(`hero_temporal.py --store-screen-phase`). A record that holds none falls back
to the WRAPPED SLOPES of the receive field, which is the TERRESTRIAL source.
The slope route reads an AO(21) correction 4 to 6 percent LOW on a space link
(see CLAUDE.md and `validation/waveoptics_ao/`), so it is a FALLBACK, not the
route of record. `--source auto` picks the screens when the record holds them,
and it prints a loud warning when it falls back.

WHAT IT DRAWS.
  1. figures/record_ao_timeseries_30.png. Five panels. The first four hold the
     fibre-coupled power of one stack over the full 2 s, all against the median
     of the UNCORRECTED record and on one y axis. The fifth holds a 100 ms zoom
     of all four around the deepest UNCORRECTED fade.
  2. figures/record_ao_field_30.gif. The SAME 200-frame window as
     `record_plots.py`, one ROW for each stack. The columns are the wrapped
     phase on the aperture, the intensity on the aperture and the intensity at
     the fibre tip. One shared panel holds the four power series with a marker
     at the current frame. The colour scale of a column is shared by every row,
     so the rows compare directly.

THE POWERS. The bucket power is the stored `collected_power`. The fibre-coupled
power is `collected_power * eta`, with eta the coupling efficiency of the
corrected field. The UNCORRECTED eta is the stored `smf_eta`.

THE CROP RULE (see `docs/api-waveoptics.md`): a PUPIL quantity reads the
compact crop, and a FOCAL-PLANE quantity reads the padded full grid, because
the focal pixel size lambda*f/siz reads the grid EXTENT. So this script
corrects on the CROP, which is what the package does, and it pads the corrected
crop back for the fibre-tip panel. The correction gives the same pixels either
way, because the Noll modes are zero outside the aperture mask.

Sources:
- Noll, DOI 10.1364/JOSA.66.000207. The Zernike mode order and the residual.
- Taylor, Proc. R. Soc. Lond. A 164, pp. 476 to 490 (1938),
  DOI 10.1098/rspa.1938.0032. The frozen-flow hypothesis of the time axis.
- Goodman, Introduction to Fourier Optics, ISBN 978-0974707723. The focal-plane
  field is the Fourier transform of the pupil field. See
  olb.waveoptics.mmf.focal_intensity.
- Shaklan and Roddier, Appl. Opt. 27, 2334 (1988), DOI 10.1364/AO.27.002334.
  The single-mode-fibre coupling parameter a and the optimal-focus rule.

Run from the repository root:

    python -m validation.temporal_screens.record_ao_plots
    python -m validation.temporal_screens.record_ao_plots --source screens
'''

import argparse
import json
import os
import time

import numpy as np

from olb.terminal import AO, TipTilt
from olb.waveoptics.compensation.modal import modes_from_stack
from olb.waveoptics.mmf import focal_intensity
from olb.waveoptics.run import SMF28_MODE_FIELD_RADIUS_M
from olb.waveoptics.turbulence.run import (_crop_array, _patch_field,
                                           _PostCorrector)
from validation.temporal_screens.hero_temporal import (campaign_of,
                                                       campaigns_root,
                                                       fade_durations)
from validation.temporal_screens.record_plots import (FADE_FRACTION, FIGS,
                                                      FPS, LEAD_FRAMES,
                                                      LOG_DECADES,
                                                      FOV_MODE_RADII, N_ANIM,
                                                      ZOOM_MS, _Args, db_rel,
                                                      fibre_focal_length,
                                                      pupil_cut, window_start)

ELEVATION_DEG = 30.0
RECORD = 0
DPI_GIF = 60         # the writer dpi. It holds the file under about 15 MB.

# The four correction stacks, in the plot order. The first one is the record as
# it was run. `modes_from_stack` gives the Noll mode count of each.
STACKS = (('none', ()),
          ('TipTilt', (TipTilt(),)),
          ('AO(10)', (AO(10),)),
          ('AO(50)', (AO(50),)))
COLOURS = {'none': '0.3', 'TipTilt': 'C0', 'AO(10)': 'C2', 'AO(50)': 'C1'}


def record_settings():
    '''Give the campaign settings of the record on disk.

    The settings are the defaults of `hero_temporal.py`, plus the
    `store_screen_phase` flag of the manifest. That flag enters the campaign
    fingerprint, so a wrong value makes the campaign refuse the root.

    Returns:
        An `_Args` instance.
    '''
    settings = _Args()
    path = os.path.join(campaigns_root(settings),
                        f'el{ELEVATION_DEG:02.0f}', f'r{RECORD}',
                        'manifest.json')
    if os.path.isfile(path):
        with open(path, encoding='utf-8') as fh:
            settings.store_screen_phase = bool(
                json.load(fh).get('store_screen_phase', False))
    return settings


def resolve_source(camp, wanted):
    '''Give the sensing source of the correction, and say when it is a fallback.

    A SPACE link senses the summed screen phase. A record that holds none can
    only use the wrapped slopes of the field, which reads a high-order
    correction low. So the fallback prints a loud warning.

    Args:
        camp:   the Campaign.
        wanted: "screens", "slopes" or "auto".

    Returns:
        The source name, a string.
    '''
    if wanted != 'auto':
        return wanted
    if camp.store_screen_phase:
        return 'screens'
    print('  WARNING  this record holds NO summed screen phase, so the read '
          'falls back to the')
    print('           wrapped SLOPES of the field. That is the TERRESTRIAL '
          'source. On a SPACE link it')
    print('           reads an AO(21) correction 4 to 6 percent LOW. Run '
          'hero_temporal.py with')
    print('           --store-screen-phase for the route of record.')
    return 'slopes'


def eta_of_stacks(camp, source, n_frames):
    '''Give the fibre coupling efficiency of every frame, for every stack.

    The uncorrected value is the STORED `smf_eta`, which the runner computed.
    Every other value comes from the post-hoc perfect-AO read of the package.

    Args:
        camp:     the Campaign.
        source:   "screens" or "slopes".
        n_frames: the number of frames to read.

    Returns:
        The pair (the bucket power of each frame, a dict of the coupling
        efficiency of each frame for each stack).
    '''
    result = camp.load(int(n_frames), fields=False)
    bucket = np.array([t.collected_power for t in result.trials], dtype=float)
    etas = {'none': np.array([t.smf_eta for t in result.trials], dtype=float)}
    print(f"  none:      0 Noll modes, mean eta {etas['none'].mean():.4f} "
          f"(the stored value)")
    for name, stack in STACKS[1:]:
        t0 = time.perf_counter()
        etas[name] = camp.recouple_compensated(
            list(stack), camp.scenario.ground.detector, n_trials=int(n_frames),
            source=source)
        print(f'  {name}: {modes_from_stack(stack):3d} Noll modes, mean eta '
              f'{etas[name].mean():.4f} ({time.perf_counter() - t0:.1f} s)',
              flush=True)
    return bucket, etas


def stack_row(name, fibre, reference, level_none, dt_s):
    '''Give the table row of one stack.

    Args:
        name:       the stack name.
        fibre:      the fibre-coupled power of every frame.
        reference:  the median fibre power of the UNCORRECTED record.
        level_none: the 5 percent fibre power level of the UNCORRECTED record.
        dt_s:       the time step, in s.

    Returns:
        A dict of the numbers of the row.
    '''
    own = float(np.quantile(fibre, FADE_FRACTION))
    shared = fade_durations(fibre, level_none, dt_s)
    private = fade_durations(fibre, own, dt_s)
    return {'stack': name,
            'median_db': float(db_rel(np.median(fibre), reference)),
            'p5_db': float(db_rel(own, reference)),
            'n_shared': int(shared.size),
            'mean_shared_ms': float(shared.mean() * 1e3) if shared.size else 0.0,
            'n_own': int(private.size),
            'mean_own_ms': float(private.mean() * 1e3) if private.size else 0.0}


def print_table(rows):
    '''Print the fade table of every stack.'''
    print('')
    print('  the level columns are dB over the UNCORRECTED median. A fade '
          'event is a maximal run of')
    print('  frames under a level. The "shared" columns use the UNCORRECTED 5 '
          'percent level, so they')
    print('  count the SAME fades; the "own" columns use the 5 percent level '
          'of that stack.')
    print(f"  {'stack':>8}{'median [dB]':>14}{'5 pc [dB]':>12}"
          f"{'shared N':>10}{'shared [ms]':>13}{'own N':>8}{'own [ms]':>11}")
    for r in rows:
        print(f"  {r['stack']:>8}{r['median_db']:>14.2f}{r['p5_db']:>12.2f}"
              f"{r['n_shared']:>10d}{r['mean_shared_ms']:>13.3f}"
              f"{r['n_own']:>8d}{r['mean_own_ms']:>11.3f}")
    print('')


def draw_timeseries(t_s, fibre, deep, path):
    '''Draw the five-panel power plot of the four stacks. Give the path back.

    Args:
        t_s:   the time of every frame, in s.
        fibre: a dict of the fibre-coupled power of each stack.
        deep:  the frame of the deepest UNCORRECTED fade.
        path:  the output file.

    Returns:
        The path.
    '''
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    reference = float(np.median(fibre['none']))
    db = {name: db_rel(fibre[name], reference) for name, _ in STACKS}
    dt_s = float(t_s[1] - t_s[0])
    half = int(round(0.5 * ZOOM_MS * 1e-3 / dt_s))
    lo = window_start(deep, t_s.size, 2 * half, half)
    sl = slice(lo, lo + 2 * half)

    fig, ax = plt.subplots(5, 1, figsize=(9.0, 12.0))
    for a, (name, stack) in zip(ax[:4], STACKS):
        a.plot(t_s, db[name], lw=0.4, color=COLOURS[name])
        a.set_ylabel(f'{name} [dB]')
        a.text(0.01, 0.06, f'{modes_from_stack(stack)} Noll modes, median '
                           f'{np.median(db[name]):+.2f} dB, 5 percent '
                           f'{np.quantile(db[name], FADE_FRACTION):+.2f} dB',
               transform=a.transAxes, fontsize=8)
        a.set_xlim(t_s[0], t_s[-1])
        a.set_ylim(min(d.min() for d in db.values()) - 1.0,
                   max(d.max() for d in db.values()) + 1.0)
        a.grid(alpha=0.3)
    ax[0].set_title(f'the hero downlink at {ELEVATION_DEG:.0f} deg, record '
                    f'{RECORD}: four perfect-AO stacks,\nthe fibre power over '
                    f'the UNCORRECTED median')
    ax[3].set_xlabel('time [s]')
    for name, _ in STACKS:
        ax[4].plot(t_s[sl] * 1e3, db[name][sl], lw=0.9, color=COLOURS[name],
                   label=name)
    ax[4].axvline(t_s[deep] * 1e3, color='C3', lw=0.8, ls='--',
                  label='the deepest uncorrected fade')
    ax[4].set_xlabel('time [ms]')
    ax[4].set_ylabel('fibre power [dB]')
    ax[4].set_title(f'the {ZOOM_MS:.0f} ms around the deepest uncorrected fade')
    ax[4].grid(alpha=0.3)
    ax[4].legend(fontsize=8, loc='lower right', ncol=2)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)
    return path


def read_window(camp, rows, correctors, aperture_m, focal_length_m, half_px):
    '''Read the panel arrays of every frame of the window, for every stack.

    The function reads ONE BLOCK FILE at a time, and it corrects the same
    stored field with each corrector. So one frame costs one crop rebuild, not
    four block reads. It corrects on the CROP and it pads the result back for
    the fibre-tip panel (the crop rule, see the module docstring).

    ponytail: it holds the whole window in memory (about 80 MB at 200 frames
    and four stacks). Stream it frame by frame if a window ever gets long.

    Args:
        camp:           the Campaign.
        rows:           the frame indices, in order.
        correctors:     a dict of the `_PostCorrector` of each stack. The value
                        None leaves the field as it was stored.
        aperture_m:     the receive aperture DIAMETER, in m.
        focal_length_m: the focal length of the fibre coupler, in m.
        half_px:        the half-width of the fibre-tip panel, in focal pixels.

    Returns:
        A dict of the stacked panels of each stack, and the two pixel sizes.
    '''
    lam = camp.scenario.ground.wavelength_m
    patch, side = camp.patch, camp.patch.crop().side
    offset = camp.patch.crop().offset
    panels = {name: {'phase': [], 'pupil': [], 'focal': []}
              for name in correctors}
    cut = None
    rows = [int(r) for r in rows]
    for b in sorted({r // camp.block_size for r in rows}):
        # ponytail: the block reader is private. The public `Campaign.field`
        # gives no screen phase, and it re-reads the whole block file for each
        # frame. Make it public if a second study needs the pair.
        cols = camp._read_block(b, fields=True)
        for row in [r for r in rows if r // camp.block_size == b]:
            k = row - b * camp.block_size
            crop = _crop_array(patch, cols['fields'][k], compact=True)
            phase = (None if cols['screen_phase'] is None
                     else cols['screen_phase'][k])
            for name, corrector in correctors.items():
                E = (crop if corrector is None
                     else corrector.correct(crop, phase))
                full = np.zeros((int(patch.n), int(patch.n)), dtype=E.dtype)
                full[offset:offset + side, offset:offset + side] = E
                field = _patch_field(patch, full, lam)
                if cut is None:
                    cut = pupil_cut(field, aperture_m)
                sl, mask, extent_m = cut
                box = full[sl, sl]
                panels[name]['phase'].append(
                    np.where(mask, np.angle(box), np.nan))
                panels[name]['pupil'].append(
                    np.where(mask, np.abs(box) ** 2, np.nan))
                image, dx_focal = focal_intensity(field, focal_length_m)
                mid = image.shape[0] // 2
                panels[name]['focal'].append(
                    image[mid - half_px:mid + half_px,
                          mid - half_px:mid + half_px])
    out = {name: {k: np.array(v) for k, v in d.items()}
           for name, d in panels.items()}
    return {'panels': out, 'dx_focal': dx_focal, 'extent_m': extent_m}


def draw_animation(window, fibre_db, t_ms, mode_radius_m, path, dpi=DPI_GIF):
    '''Draw the four-row animation of the window. Give the path back.

    Each ROW is one correction stack, and the colour scale of a COLUMN is
    shared by every row, so the rows compare directly.

    ponytail: the writer is Pillow at a small figure size, which holds the file
    under about 15 MB with no panel downsampling. Move to ffmpeg (an mp4) if a
    longer window is ever wanted.

    Args:
        window:        the dict of `read_window`.
        fibre_db:      a dict of the fibre power of each stack over the window,
                       in dB.
        t_ms:          the time of every frame of the window, in ms.
        mode_radius_m: the fibre mode field radius, in m.
        path:          the output file.
        dpi:           the writer dpi.

    Returns:
        The path.
    '''
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation, PillowWriter
    from matplotlib.colors import LogNorm

    panels = window['panels']
    names = [name for name, _ in STACKS]
    # The vacuum limit of the aperture intensity is its mean over the window
    # (the slab starts from a plane wave), so 1.0 is that limit.
    pupil = {n: panels[n]['pupil'] / np.nanmean(panels[n]['pupil'])
             for n in names}
    ext = window['extent_m']
    pupil_extent = [-ext, ext, -ext, ext]
    fext = 0.5 * panels[names[0]]['focal'].shape[1] * window['dx_focal'] * 1e6
    focal_extent = [-fext, fext, -fext, fext]
    pupil_max = float(max(np.nanpercentile(pupil[n], 99.9) for n in names))
    focal_max = float(max(panels[n]['focal'].max() for n in names))
    n_frames = panels[names[0]]['focal'].shape[0]

    fig = plt.figure(figsize=(9.2, 7.4))
    gs = fig.add_gridspec(4, 4, width_ratios=(1.0, 1.0, 1.0, 1.35))
    images, axes = {}, {}
    for r, name in enumerate(names):
        axes[name] = [fig.add_subplot(gs[r, c]) for c in range(3)]
        a0, a1, a2 = axes[name]
        im0 = a0.imshow(panels[name]['phase'][0], cmap='twilight',
                        vmin=-np.pi, vmax=np.pi, extent=pupil_extent,
                        origin='lower')
        im1 = a1.imshow(pupil[name][0], cmap='inferno', vmin=0.0,
                        vmax=pupil_max, extent=pupil_extent, origin='lower')
        im2 = a2.imshow(panels[name]['focal'][0], cmap='viridis',
                        extent=focal_extent, origin='lower',
                        norm=LogNorm(vmin=focal_max / 10.0 ** LOG_DECADES,
                                     vmax=focal_max))
        a2.add_patch(plt.Circle((0.0, 0.0), mode_radius_m * 1e6, fill=False,
                                color='w', lw=0.8))
        images[name] = (im0, im1, im2)
        a0.set_ylabel(name, fontsize=8)
        for a in (a0, a1, a2):
            a.tick_params(labelsize=5)
        if r == 0:
            a0.set_title('phase [rad]', fontsize=8)
            a1.set_title('intensity / vacuum', fontsize=8)
            a2.set_title('fibre tip [um]', fontsize=8)
        if r == len(names) - 1:
            a0.set_xlabel('x [m]', fontsize=7)
            a1.set_xlabel('x [m]', fontsize=7)
            a2.set_xlabel('x [um]', fontsize=7)
    strip = fig.add_subplot(gs[:, 3])
    markers = {}
    for name in names:
        strip.plot(t_ms, fibre_db[name], lw=0.7, color=COLOURS[name],
                   label=name)
        markers[name], = strip.plot([t_ms[0]], [fibre_db[name][0]], 'o',
                                    color=COLOURS[name], ms=4)
    strip.set_xlabel('time [ms]', fontsize=7)
    strip.set_ylabel('fibre power over the uncorrected median [dB]',
                     fontsize=7)
    strip.tick_params(labelsize=6)
    strip.grid(alpha=0.3)
    strip.legend(fontsize=7, loc='lower right')
    title = fig.suptitle('', fontsize=9)
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.95))

    def update(k):
        '''Put frame k on every panel.'''
        out = [title]
        for name in names:
            im0, im1, im2 = images[name]
            im0.set_data(panels[name]['phase'][k])
            im1.set_data(pupil[name][k])
            im2.set_data(panels[name]['focal'][k])
            markers[name].set_data([t_ms[k]], [fibre_db[name][k]])
            out += [im0, im1, im2, markers[name]]
        title.set_text('t = {:.1f} ms, fibre '.format(t_ms[k]) + ', '.join(
            f'{n} {fibre_db[n][k]:.1f} dB' for n in names))
        return out

    anim = FuncAnimation(fig, update, frames=n_frames, blit=False)
    anim.save(path, writer=PillowWriter(fps=FPS), dpi=int(dpi))
    plt.close(fig)
    return path


def _self_check():
    '''Check the stacks and the table row. It runs on every call.'''
    assert [modes_from_stack(s) for _, s in STACKS] == [0, 3, 10, 50]
    # A flat series of ones, with two frames pushed under both levels. The
    # median is 1.0, so the median column is 0.0 dB.
    dt = 1e-3
    x = np.ones(100)
    x[10] = x[50] = 0.5
    row = stack_row('t', x, 1.0, 0.75, dt)
    assert np.isclose(row['median_db'], 0.0), row
    assert row['n_shared'] == 2 and np.isclose(row['mean_shared_ms'], 1.0), row
    # The own 5 percent level of that series is 1.0, so every frame of the
    # 98 frames at 1.0 is NOT under it and only the two dips count.
    assert row['n_own'] >= 2, row


def main():
    '''Read the record, correct it four ways, and draw the two figures.'''
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument('--source', default='auto',
                        choices=['auto', 'screens', 'slopes'],
                        help='the sensing source of the correction. auto takes '
                             'the screens when the record holds them.')
    parser.add_argument('--frames', type=int, default=N_ANIM,
                        help='the frames of the animation.')
    parser.add_argument('--trials', type=int, default=None,
                        help='the frames of the time series. The default is '
                             'every stored frame.')
    parser.add_argument('--dpi', type=int, default=DPI_GIF,
                        help='the writer dpi of the animation.')
    args = parser.parse_args()

    _self_check()
    t_start = time.time()
    camp = campaign_of(ELEVATION_DEG, RECORD, record_settings())
    n = int(camp.n_stored) if args.trials is None else int(args.trials)
    print(f'record {camp.root_dir}')
    print(f'  {n} frames at {camp.dt_s * 1e3:.2f} ms, grid {camp.grid.n} px at '
          f'{camp.grid.pixel_m * 1e3:.2f} mm, patch {camp.patch.radius_m:.2f} m')
    aperture_m = camp.scenario.ground.aperture_m
    obscuration = camp.scenario.ground.obscuration_ratio
    print(f'  the aperture is {aperture_m:.2f} m, which is '
          f'{aperture_m / camp.grid.pixel_m:.1f} px, and the crop is '
          f'{camp.patch.crop().side} px')
    source = resolve_source(camp, args.source)
    print(f'  the sensing source is {source!r}')

    bucket, etas = eta_of_stacks(camp, source, n)
    fibre = {name: bucket * etas[name] for name, _ in STACKS}
    reference = float(np.median(fibre['none']))
    level_none = float(np.quantile(fibre['none'], FADE_FRACTION))
    print_table([stack_row(name, fibre[name], reference, level_none,
                           camp.dt_s) for name, _ in STACKS])

    t_s = camp.t_s(np.arange(n))
    deep = int(np.argmin(fibre['none']))
    print(f'  the deepest uncorrected fade is '
          f'{db_rel(fibre["none"][deep], reference):.2f} dB under the median, '
          f'at frame {deep} (t = {t_s[deep] * 1e3:.1f} ms)')
    os.makedirs(FIGS, exist_ok=True)
    path = draw_timeseries(t_s, fibre, deep,
                           os.path.join(FIGS, f'record_ao_timeseries_'
                                              f'{ELEVATION_DEG:.0f}.png'))
    print(f'  figure saved: {path} ({os.path.getsize(path) / 2 ** 20:.2f} MB)')

    focal_length_m = fibre_focal_length(camp.scenario)
    dx_focal = (camp.scenario.ground.wavelength_m * focal_length_m
                / camp.grid.size_m)
    half_px = int(np.ceil(FOV_MODE_RADII * SMF28_MODE_FIELD_RADIUS_M
                          / dx_focal))
    correctors = {'none': None}
    for name, stack in STACKS[1:]:
        correctors[name] = _PostCorrector(camp.patch, list(stack), aperture_m,
                                          obscuration, source, compact=True)
    lo = window_start(deep, n, int(args.frames), LEAD_FRAMES)
    rows = np.arange(lo, lo + int(args.frames))
    window = read_window(camp, rows, correctors, aperture_m, focal_length_m,
                         half_px)
    path = draw_animation(
        window, {name: db_rel(fibre[name][rows], reference)
                 for name, _ in STACKS},
        t_s[rows] * 1e3, SMF28_MODE_FIELD_RADIUS_M,
        os.path.join(FIGS, f'record_ao_field_{ELEVATION_DEG:.0f}.gif'),
        dpi=args.dpi)
    print(f'  animation saved: {path} '
          f'({os.path.getsize(path) / 2 ** 20:.2f} MB, {len(rows)} frames at '
          f'{FPS} fps, t = {t_s[rows[0]] * 1e3:.1f} to '
          f'{t_s[rows[-1]] * 1e3:.1f} ms)')
    print(f'(elapsed {time.time() - t_start:.1f} s)')


if __name__ == '__main__':
    main()
