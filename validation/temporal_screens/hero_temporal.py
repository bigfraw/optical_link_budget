'''
Gate (e) of the temporal strip screens: the fade RATE and the fade DURATION.

The four earlier gates show that a strip crop carries the right atmosphere
(gate c) and that it changes at the right rate (gates b and d). This script is
the PRODUCT of that work: it measures how OFTEN the hero downlink fades below
its 5 percent level, and how LONG each fade lasts. A snapshot campaign cannot
give either number, because its trials have no order in time.

THE CASE. The hero downlink of `validation/waveoptics_ao/`: 1550 nm, 500 km, a
700 mm ground telescope with a single-mode fibre, a 100 mm space terminal. It
runs at 30 deg and at 20 deg, at the site outer scale (25 m), in single
precision, with the `standard` preset.

THE RECORD. One record is 4000 frames at 0.5 ms, so it is 2.0 s long. A 2 s
record holds only about six independent 5 percent fades, so the study runs
EIGHT records of each elevation. Each record draws its own strips, so the eight
are independent, and they pool into about 50 events (a Poisson bar of about
14 percent).

WHAT IT REPORTS, for each elevation and for each of the two receivers:
  - the 5 percent fade LEVEL, the 5th percentile of the pooled power,
  - the number of fade EVENTS, where an event is a maximal run of consecutive
    frames under that level,
  - the fade RATE in events per second, with its Poisson bar sqrt(N) / T,
  - the MEAN and the MEDIAN fade duration in ms, each with a 2 SE bar from a
    bootstrap over the events.
THE TWO RECEIVERS are the single-mode fibre, p(t) = collected_power * smf_eta,
and the bucket, p(t) = collected_power. The fibre fade is phase-dominated and
the bucket fade is the aperture-averaged scintillation, so the two carry
different time scales.

THE CONTEXT NUMBERS. The script also prints the Greenwood frequency and the
coherence time tau0 of the SAME per-layer velocities that the strips use
(`StripPlan`), and the layer speed table. They are CONTEXT, not a test: the
Greenwood frequency is a phase quantity that does not know the pupil, so it
sits well above the tilt corner of a 0.7 m telescope (see gate (d)).

THE GATE has NO numeric band on the durations, because no reference model of
record gives them. Its requirement is that the event count and its Poisson bar
are REPORTED. The one PASS line is `n_events >= 20` for each elevation, so a
reader knows when the record set is too short to say anything.

THE CEILINGS of this script are marked `ponytail:` where they sit.

Sources:
- Taylor, Proc. R. Soc. Lond. A 164, pp. 476 to 490 (1938),
  DOI 10.1098/rspa.1938.0032. The frozen-flow hypothesis that gives the time
  axis.
- Greenwood, J. Opt. Soc. Am. 67(3), pp. 390 to 393 (1977),
  DOI 10.1364/JOSA.67.000390. The Greenwood frequency.
- Andrews and Phillips, Laser Beam Propagation through Random Media, 2nd ed.,
  DOI 10.1117/3.626196, Ch. 14, Eqs. (38) and (39), printed p. 622. The
  Greenwood frequency of a slant path and the coherence time tau0 = 0.314 / f_G.
- Schmidt, DOI 10.1117/3.866274, Ch. 9. The split-step method of every frame.

Outputs, next to this script:
    data/hero_temporal_<el>.csv       one row for each record and receiver
    figures/hero_temporal_<el>.png    p(t) of record 0 and the duration
                                      histogram

Run from the repository root. The smoke run takes the host backend:

    python -m validation.temporal_screens.hero_temporal --smoke
    python -m validation.temporal_screens.hero_temporal
    python -m validation.temporal_screens.hero_temporal --analyse
'''

import argparse
import csv
import os
import shutil
import time
import warnings

import numpy as np

from olb.turbulence.andrews.temporal import coherence_time, greenwood_frequency
from olb.turbulence.profiles import default_cn2_profile
from olb.waveoptics.turbulence import Campaign
from olb.waveoptics.turbulence.temporal import TemporalSpec, strip_plan
from validation.temporal_screens.tilt_spectrum import (layer_heights,
                                                       layer_speeds)
from validation.waveoptics_ao.waveoptics_ao import hero_scenario

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, 'data')
FIGS = os.path.join(HERE, 'figures')
ENV_ROOT = 'OLB_TEMPORAL_ROOT'

# ---- the case ----
ELEVATIONS_DEG = (30.0, 20.0)
N_RECORDS = 8
DT_S = 5e-4                  # the sample step, in s. It resolves tau0 (7 ms).
N_FRAMES = 4000              # 2.0 s of record.
PRESET = 'standard'
PRECISION = 'single'
SEED = 20260913
PATCH_RADIUS_M = 0.35        # half the hero ground aperture, as gate (d) uses.
FADE_FRACTION = 0.05         # the 5 percent fade level.
MIN_EVENTS = 20              # the PASS level of the pooled event count.
N_BOOT = 400                 # the bootstrap resamples of a duration bar.
BOOT_SEED = 12345

# The two receivers. The value is the pair of scalar columns that multiply into
# the received power of that receiver.
RECEIVERS = {'smf': ('collected_power', 'smf_eta'),
             'bucket': ('collected_power',)}


# ---------------------------------------------------------------------------
# The campaigns
# ---------------------------------------------------------------------------

def campaigns_root(args):
    '''Give the parent directory of every campaign of this study.

    Args:
        args: the parsed command line.

    Returns:
        A path. The smoke run keeps its own sub-directory, because it holds a
        different preset and a different frame count, and a campaign refuses a
        root whose manifest does not match.
    '''
    root = os.environ.get(ENV_ROOT) or os.path.join(HERE, 'campaigns')
    return os.path.join(root, 'smoke') if args.smoke else root


def campaign_of(elevation_deg, record, args):
    '''Build (or reopen) the campaign of one elevation and one record.

    The frames of the record ARE the trials of the campaign: frame k is trial k
    and it sits at the time k * dt_s. The strips go under the campaign root, so
    a deleted campaign takes its cache with it.

    Args:
        elevation_deg: the elevation of the line of sight, in deg.
        record:        the record index, from 0.
        args:          the parsed command line.

    Returns:
        A Campaign.
    '''
    scn, geom = hero_scenario(elevation_deg)
    root = os.path.join(campaigns_root(args), f'el{elevation_deg:02.0f}',
                        f'r{int(record)}')
    spec = TemporalSpec(dt_s=float(args.dt), n_frames=int(args.frames),
                        strip_dir='unused', record=int(record))
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        # L0_m stays out: None reads the site outer scale, which is the owner
        # value of 25 m (backlog 2-P5).
        return Campaign(scn, geom, root, seed=SEED, preset=args.preset,
                        block_size=int(args.block_size),
                        patch_radius_m=PATCH_RADIUS_M, precision=PRECISION,
                        fft_backend=args.fft_backend, temporal=spec,
                        store_screen_phase=bool(args.store_screen_phase))


def ensure_trials(camp, n_frames, args):
    '''Run the missing frames of one record, and report the wall time.

    Args:
        camp:     the Campaign.
        n_frames: the wanted frame count.
        args:     the parsed command line.

    Returns:
        The wall time of the call, in s.
    '''
    missing = max(0, int(n_frames) - int(camp.n_stored))
    if missing == 0:
        print(f'    {camp.n_stored} frames already on disk.', flush=True)
        return 0.0
    t0 = time.perf_counter()
    n_done = camp.run(int(n_frames), workers=args.workers, progress=False)
    wall = time.perf_counter() - t0
    print(f'    {n_done} frames on disk, {wall:.1f} s wall for {missing} new '
          f'frames ({wall / missing:.3f} s/frame)', flush=True)
    return float(wall)


def drop_strips(camp):
    '''Delete the strips of a complete record.

    A strip is a deletable cache: the seed of the campaign rebuilds it, bit for
    bit. One record of the hero downlink holds about 0.5 GB of strips, so eight
    records of two elevations would hold about 8 GB for nothing.

    Args:
        camp: the Campaign.

    Returns:
        The number of bytes that the call freed.
    '''
    path = os.path.join(camp.root_dir, 'strips')
    if not os.path.isdir(path):
        return 0
    freed = sum(os.path.getsize(os.path.join(path, f))
                for f in os.listdir(path)
                if os.path.isfile(os.path.join(path, f)))
    shutil.rmtree(path, ignore_errors=True)
    return int(freed)


# ---------------------------------------------------------------------------
# The fade statistics
# ---------------------------------------------------------------------------

def series_of(camp, columns, n_frames):
    '''Give the received power of every frame of one record.

    Args:
        camp:     the Campaign.
        columns:  the scalar column names that multiply into the power.
        n_frames: the number of frames to read.

    Returns:
        A float array of the power of each frame, in the frame order.
    '''
    result = camp.load(int(n_frames), fields=False)
    out = np.ones(len(result.trials), dtype=float)
    for name in columns:
        out *= np.array([getattr(t, name) for t in result.trials], dtype=float)
    return out


def fade_durations(power, level, dt_s):
    '''Give the duration of every fade event of one record, in s.

    A FADE EVENT is a maximal run of consecutive frames whose power sits under
    the level, and its duration is the frame count times the step. An event that
    touches the start or the end of the record is KEPT and it is therefore cut
    short.

    ponytail: the truncated edge events bias the mean duration LOW by about
    one mean duration for each record (about 2 percent here, at 8 records of
    2 s against a mean of a few ms). Drop them if a record ever gets short
    against its events.

    Args:
        power: the power of each frame.
        level: the fade level. A frame under it is in a fade.
        dt_s:  the time step, in s.

    Returns:
        A float array of the durations, in s.
    '''
    below = np.asarray(power, dtype=float) < float(level)
    if below.size == 0 or not below.any():
        return np.zeros(0, dtype=float)
    edge = np.diff(below.astype(np.int8))
    starts = np.flatnonzero(edge == 1) + 1
    ends = np.flatnonzero(edge == -1) + 1
    if below[0]:
        starts = np.concatenate(([0], starts))
    if below[-1]:
        ends = np.concatenate((ends, [below.size]))
    return (ends - starts).astype(float) * float(dt_s)


def bootstrap_half(values, fn, n_boot=N_BOOT, seed=BOOT_SEED):
    '''Give the 2 SE half-width of a statistic, from a bootstrap.

    The function resamples the events with replacement, it takes the statistic
    of each resample, and it gives two times the standard deviation of those
    draws. The seed is FIXED, so a rerun gives the same bar.

    Args:
        values: the sample.
        fn:     the statistic.
        n_boot: the number of resamples.
        seed:   the seed of the resampler.

    Returns:
        The half-width, as a float. It is NaN for an empty sample.
    '''
    x = np.asarray(values, dtype=float)
    if x.size < 2:
        return float('nan')
    rng = np.random.default_rng(seed)
    draws = np.array([fn(x[rng.integers(0, x.size, x.size)])
                      for _ in range(int(n_boot))], dtype=float)
    return float(2.0 * draws.std())


def fade_row(durations, total_s):
    '''Give the rate and the duration statistics of one pooled event set.

    The rate is N / T and its bar is the Poisson counting error sqrt(N) / T.

    Args:
        durations: the duration of every event, in s.
        total_s:   the total record time, in s.

    Returns:
        A dict of the statistics. The durations are in ms.
    '''
    d = np.asarray(durations, dtype=float)
    n = int(d.size)
    return {
        'n_events': n,
        'total_s': float(total_s),
        'rate_hz': n / float(total_s),
        'rate_half_hz': np.sqrt(n) / float(total_s),
        'mean_ms': float(d.mean() * 1e3) if n else float('nan'),
        'mean_half_ms': bootstrap_half(d * 1e3, np.mean),
        'median_ms': float(np.median(d) * 1e3) if n else float('nan'),
        'median_half_ms': bootstrap_half(d * 1e3, np.median),
        'max_ms': float(d.max() * 1e3) if n else float('nan'),
    }


# ---------------------------------------------------------------------------
# The context numbers
# ---------------------------------------------------------------------------

def timescales(camp, spec, elevation_deg):
    '''Give the Greenwood frequency, tau0 and the layer table of one record.

    The Greenwood integral reads the site Cn2 on a fine height grid and the
    SAME per-layer velocity that the strips use, interpolated onto that grid.
    The default Bufton wind of `greenwood_frequency` is NOT used, because it
    carries its own slew term and this record already holds the slew in the
    layer velocity (see gate (d)).

    Sources: Greenwood, DOI 10.1364/JOSA.67.000390; Andrews and Phillips,
    DOI 10.1117/3.626196, Ch. 14, Eqs. (38) and (39), printed p. 622.

    Args:
        camp:          the Campaign.
        spec:          the TemporalSpec of the record.
        elevation_deg: the elevation, in deg.

    Returns:
        A dict of f_G, tau0 and the per-layer height, speed and r0.
    '''
    sp = strip_plan(camp.plan, camp.grid, spec, camp.geometry, camp.L0_m)
    speeds = layer_speeds(sp)
    heights = layer_heights(camp.plan, elevation_deg)
    lam = camp.scenario.ground.wavelength_m
    hs = np.concatenate(([0.0], np.geomspace(1.0, float(heights.max()), 512)))
    cn2 = default_cn2_profile(camp.scenario.channel.site, hs)
    wind = np.interp(hs, heights[::-1], speeds[::-1])
    f_g = float(greenwood_frequency(hs, cn2, lam,
                                    elevation_deg=float(elevation_deg),
                                    wind_profile=wind))
    return {'f_greenwood_hz': f_g,
            'tau0_s': float(coherence_time(f_g)),
            'heights_m': heights, 'speeds_m_s': speeds,
            'r0_m': np.asarray(camp.plan.r0_m, dtype=float)}


# ---------------------------------------------------------------------------
# The outputs
# ---------------------------------------------------------------------------

def write_csv(elevation_deg, rows):
    '''Write one row for each record and receiver. Give the path back.'''
    os.makedirs(DATA, exist_ok=True)
    path = os.path.join(DATA, f'hero_temporal_{elevation_deg:02.0f}.csv')
    fields = ['record', 'receiver', 'n_events', 'total_s', 'rate_hz',
              'rate_half_hz', 'mean_ms', 'mean_half_ms', 'median_ms',
              'median_half_ms', 'max_ms', 'level']
    with open(path, 'w', newline='', encoding='utf-8') as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, '') for k in fields})
    return path


def draw_figure(elevation_deg, times, power, level, durations_ms):
    '''Draw p(t) of record 0 and the duration histogram. Give the path back.'''
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    os.makedirs(FIGS, exist_ok=True)
    fig, (top, low) = plt.subplots(2, 1, figsize=(9.0, 7.0))
    top.plot(times, 10.0 * np.log10(power / np.median(power)), lw=0.6,
             color='0.3', label='record 0')
    top.axhline(10.0 * np.log10(level / np.median(power)), color='C3', ls='--',
                lw=1.2, label=f'{FADE_FRACTION * 100:.0f} percent level')
    top.set_xlabel('time [s]')
    top.set_ylabel('fibre-coupled power / median [dB]')
    top.set_title(f'gate (e): the hero downlink at {elevation_deg:.0f} deg, '
                  f'one 2 s record')
    top.grid(alpha=0.3)
    top.legend(fontsize=8)
    if durations_ms.size:
        low.hist(durations_ms, bins=24, color='C0')
    low.set_xlabel('fade duration [ms]')
    low.set_ylabel('events')
    low.set_title(f'the pooled fade durations, {durations_ms.size} events')
    low.grid(alpha=0.3)
    fig.tight_layout()
    path = os.path.join(FIGS, f'hero_temporal_{elevation_deg:02.0f}.png')
    fig.savefig(path, dpi=140)
    plt.close(fig)
    return path


# ---------------------------------------------------------------------------
# The driver
# ---------------------------------------------------------------------------

def analyse(elevation_deg, camps, args):
    '''Pool the records of one elevation and report the fade statistics.

    Args:
        elevation_deg: the elevation, in deg.
        camps:         the Campaign of each record.
        args:          the parsed command line.

    Returns:
        True when the pooled single-mode-fibre event count holds the PASS band.
    '''
    n_frames = min([int(c.n_stored) for c in camps] + [int(args.frames)])
    if n_frames < int(args.frames):
        print(f'  INFO  only {n_frames} of {args.frames} frames are on every '
              f'record. The analysis reads {n_frames}.')
    if n_frames == 0:
        print('  INFO  no frame is stored. There is nothing to analyse.')
        return False
    dt_s = float(args.dt)
    total_s = len(camps) * n_frames * dt_s

    scales = timescales(camps[0], camps[0].temporal, elevation_deg)
    print(f'  grid                     {camps[0].grid.n:9d} px, '
          f'{camps[0].grid.pixel_m * 1e3:.2f} mm')
    print(f'  screens                  {camps[0].plan.z_m.size:9d}')
    print(f'  outer scale              {camps[0].L0_m:9.1f} m')
    print(f'  records                  {len(camps):9d} of {n_frames} frames '
          f'at {dt_s * 1e3:.2f} ms ({total_s:.2f} s in total)')
    print(f'  Greenwood f_G            '
          f'{scales["f_greenwood_hz"]:9.1f} Hz (context only)')
    print(f'  coherence time tau0      {scales["tau0_s"] * 1e3:9.3f} ms '
          f'(context only)')
    print('')
    print('  the layer velocity of the strips (slew + Bufton, ws = 0)')
    print(f"  {'layer':>6}{'h [m]':>12}{'v [m/s]':>10}{'r0 [cm]':>10}")
    for j, (h, v, r0) in enumerate(zip(scales['heights_m'],
                                       scales['speeds_m_s'],
                                       scales['r0_m'])):
        print(f'  {j:>6}{h:>12.0f}{v:>10.2f}{r0 * 1e2:>10.2f}')
    print('')

    rows, pooled, passed = [], {}, False
    for name, columns in RECEIVERS.items():
        series = [series_of(c, columns, n_frames) for c in camps]
        level = float(np.quantile(np.concatenate(series), FADE_FRACTION))
        per_record = [fade_durations(s, level, dt_s) for s in series]
        for record, durations in enumerate(per_record):
            row = fade_row(durations, n_frames * dt_s)
            row.update({'record': record, 'receiver': name, 'level': level})
            rows.append(row)
        total = fade_row(np.concatenate(per_record), total_s)
        total.update({'record': 'pooled', 'receiver': name, 'level': level})
        rows.append(total)
        pooled[name] = {'row': total, 'series': series,
                        'durations_ms': np.concatenate(per_record) * 1e3,
                        'level': level}

        ok = total['n_events'] >= MIN_EVENTS
        print(f'  {name}: the {FADE_FRACTION * 100:.0f} percent fade level is '
              f'{level:.4g} (a loss of '
              f'{-10.0 * np.log10(level / np.median(np.concatenate(series))):.2f}'
              f' dB under the median)')
        print(f'    events                 {total["n_events"]:9d}')
        print(f'    rate                   {total["rate_hz"]:9.2f} +/- '
              f'{total["rate_half_hz"]:.2f} /s (the bar is Poisson)')
        print(f'    mean duration          {total["mean_ms"]:9.3f} +/- '
              f'{total["mean_half_ms"]:.3f} ms (the bar is 2 SE)')
        print(f'    median duration        {total["median_ms"]:9.3f} +/- '
              f'{total["median_half_ms"]:.3f} ms (the bar is 2 SE)')
        print(f'    longest event          {total["max_ms"]:9.3f} ms')
        print(f'  {"PASS" if ok else "INFO"}  {name} events '
              f'{total["n_events"]} against the {MIN_EVENTS} needed for a '
              f'usable bar')
        if name == 'smf':
            passed = ok
        print('')

    print(f'  file saved: {write_csv(elevation_deg, rows)}')
    if args.figures:
        smf = pooled['smf']
        times = np.arange(n_frames) * dt_s
        print(f'  figure saved: '
              f'{draw_figure(elevation_deg, times, smf["series"][0], smf["level"], smf["durations_ms"])}')
        print('    Caption: the fibre-coupled power of one 2 s record with the '
              '5 percent fade level, and')
        print('    the histogram of the pooled fade durations.')
    return passed


def _self_check():
    '''Check the event finder against hand cases. It runs on every call.'''
    dt = 0.5
    # Two events, one of them at the END of the record, so it is truncated.
    d = fade_durations([1.0, 0.1, 0.1, 1.0, 1.0, 0.1], 0.5, dt)
    assert np.allclose(d, [1.0, 0.5]), d
    # An event at the START of the record.
    assert np.allclose(fade_durations([0.1, 1.0], 0.5, dt), [0.5])
    # No frame under the level, and every frame under it.
    assert fade_durations([1.0, 1.0], 0.5, dt).size == 0
    assert np.allclose(fade_durations([0.1, 0.1], 0.5, dt), [1.0])
    row = fade_row(np.array([0.001, 0.003]), 2.0)
    assert row['n_events'] == 2 and np.isclose(row['rate_hz'], 1.0), row
    assert np.isclose(row['mean_ms'], 2.0), row


def main():
    '''Run the campaigns, then the analysis, and print the gate lines.'''
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument('--elevations', nargs='+', type=float,
                        default=list(ELEVATIONS_DEG),
                        help='the elevations of the records [deg].')
    parser.add_argument('--records', type=int, default=N_RECORDS,
                        help='the records of each elevation.')
    parser.add_argument('--frames', type=int, default=N_FRAMES,
                        help='the frames of one record.')
    parser.add_argument('--dt', type=float, default=DT_S,
                        help='the time step of one frame [s].')
    parser.add_argument('--preset', default=PRESET,
                        help='the sampling preset of the campaigns.')
    parser.add_argument('--block-size', type=int, default=500,
                        help='the frames in one block file.')
    parser.add_argument('--fft-backend', default='cupy',
                        choices=['numpy', 'scipy', 'cupy'],
                        help='the FFT backend of the campaigns.')
    parser.add_argument('--workers', default=None,
                        help='the campaign pool size. Leave it out for the '
                             'CUDA backend: one device runs one stream.')
    parser.add_argument('--analyse', action='store_true',
                        help='read what is stored and compute no frame.')
    parser.add_argument('--store-screen-phase', action='store_true',
                        help='keep the summed screen phase of every frame. '
                             'The perfect-AO read of a SPACE link senses that '
                             'phase (record_ao_plots.py). It adds one float32 '
                             'array of the patch pixels to each frame, and it '
                             'enters the campaign fingerprint, so a record '
                             'that holds it is a NEW campaign.')
    parser.add_argument('--keep-strips', action='store_true',
                        help='keep the strips of a complete record.')
    parser.add_argument('--no-figures', dest='figures', action='store_false',
                        help='skip the figures.')
    parser.add_argument('--smoke', action='store_true',
                        help='one short record at 30 deg on the host, to '
                             'check the script end to end.')
    args = parser.parse_args()
    if args.smoke:
        args.elevations = [30.0]
        args.records = 1
        args.frames = 200
        args.preset = 'rapid'
        args.fft_backend = 'numpy'
    if args.workers not in (None, 'auto'):
        args.workers = int(args.workers)
    # A campaign runs WHOLE blocks, so a block larger than the record would run
    # frames past the travel that the strips hold. The clamp keeps every frame
    # inside its strip.
    # ponytail: the block must also DIVIDE the frame count, or the last block
    # runs past the end. 500 divides 4000, so the default pair is safe.
    args.block_size = min(int(args.block_size), int(args.frames))

    _self_check()
    t_start = time.time()
    print(f'gate (e), the hero downlink, {args.preset} preset, {PRECISION} '
          f'precision, seed {SEED}')
    print(f'  backend {args.fft_backend}, {args.records} records of '
          f'{args.frames} frames at {args.dt * 1e3:.2f} ms')
    print(f'  campaigns under {campaigns_root(args)}')
    print('')

    passes = []
    for elevation in args.elevations:
        print(f'ELEVATION {elevation:.0f} deg')
        camps = []
        for record in range(int(args.records)):
            camp = campaign_of(elevation, record, args)
            if not args.analyse:
                print(f'  record {record}: {camp.root_dir}')
                ensure_trials(camp, int(args.frames), args)
                if (not args.keep_strips
                        and camp.n_stored >= int(args.frames)):
                    freed = drop_strips(camp)
                    if freed:
                        print(f'    strips deleted, {freed / 2 ** 20:.0f} MB '
                              f'freed (the seed rebuilds them).')
            camps.append(camp)
        print('')
        passes.append(analyse(elevation, camps, args))

    print(f'  {sum(passes)} of {len(passes)} elevations hold the '
          f'{MIN_EVENTS} event band')
    print(f'(elapsed {time.time() - t_start:.1f} s)')


if __name__ == '__main__':
    main()
