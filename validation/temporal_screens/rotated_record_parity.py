'''
The ROTATED record against independent snapshots (backlog 2-P1b item 10).

THE TEST. A frozen-flow record must hold the SAME one-point statistics as a set
of INDEPENDENT snapshots of the same atmosphere. The strip only correlates the
frames in TIME (Taylor, DOI 10.1098/rspa.1938.0032); it must not change what
one frame is. `frame0_parity.py` (gate (c)) showed that for FRAME 0 of the
ALONG-TRACK route. This script shows it for the WHOLE ROTATED record, where
every frame passes through the three-shear rotation of Unser, Thevenaz and
Yaroslavsky, DOI 10.1109/83.469963.

THE TWO ARMS, both on the hero downlink at 30 deg (a 0.7 m ground telescope
with a single-mode fibre; `validation/waveoptics_ao/`):

  R. ONE rotated record of `--frames` frames at `--dt`, with a CROSSWIND
     (`wind_dir_deg = 90`), `rotated=True` and the default taper
     roll-off `rot_margin = 0.07`.
  S. The SAME number of independent snapshot trials, on the SAME grid, the
     same screen plan, the same preset and the same seed stream.

WHAT IT COMPARES, for the two receivers of gate (e) (the single-mode fibre and
the bucket) and for a near-POINT aperture:

  - the scintillation index sigma2_I = var(P) / mean(P)^2,
  - the mean power in dB,
  - the 5 percent and the 1 percent quantiles, in dB under the median of that
    arm. That pair IS the fade depth.

THE ERROR BAR IS A BOOTSTRAP of 400 resamples. The frames of a record are
CORRELATED in time, so arm R takes a moving BLOCK bootstrap of about 20 ms
blocks (Kunsch, Ann. Statist. 17(3), pp. 1217 to 1241 (1989),
DOI 10.1214/aos/1176347265). Arm S takes the ordinary bootstrap.

THE BANDS. 5 percent on a scintillation index, and 0.3 dB on a mean or a fade
quantile. The script prints every number whatever the verdict.

THE COST is reported: the wall time of each arm, the strip bytes on disk and
the peak working set of the process.

Outputs, next to this script:
    data/rotated_record_parity.csv   one row for each receiver and statistic

Run from the repository root. The smoke run takes the host backend:

    python -m validation.temporal_screens.rotated_record_parity --smoke
    python -m validation.temporal_screens.rotated_record_parity
    python -m validation.temporal_screens.rotated_record_parity --analyse

Sources:
- Taylor, Proc. R. Soc. Lond. A 164, pp. 476 to 490 (1938),
  DOI 10.1098/rspa.1938.0032. The frozen-flow hypothesis of the record.
- Unser, Thevenaz and Yaroslavsky, IEEE Trans. Image Process. 4, 1371 (1995),
  DOI 10.1109/83.469963. The three-shear rotation of every frame.
- Andrews and Phillips, Laser Beam Propagation through Random Media, 2nd ed.,
  DOI 10.1117/3.626196, Ch. 8. The scintillation index of an aperture.
- Kunsch, The jackknife and the bootstrap for general stationary observations,
  Ann. Statist. 17(3), pp. 1217 to 1241 (1989), DOI 10.1214/aos/1176347265.
  The moving block bootstrap of a correlated series.
- Schmidt, DOI 10.1117/3.866274, Ch. 9. The split-step method of every frame.
'''

import argparse
import csv
import os
import shutil
import time
import warnings

import numpy as np

from olb.waveoptics.turbulence import Campaign
from olb.waveoptics.turbulence.temporal import TemporalSpec, strip_plan
from validation.temporal_screens.table_precision import peak_working_set_bytes
from validation.waveoptics_ao.waveoptics_ao import hero_scenario

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, 'data')
ENV_ROOT = 'OLB_TEMPORAL_ROOT'

# ---- the case ----
ELEVATION_DEG = 30.0
DT_S = 5e-4                  # the sample step, in s (gate (e)).
N_FRAMES = 4000              # 2.0 s of record.
PRESET = 'standard'
PRECISION = 'single'
SEED = 20260913
PATCH_RADIUS_M = 0.35        # half the hero ground aperture.
WIND_DIR_DEG = 90.0          # the CROSSWIND that makes the rotation matter.
ROT_MARGIN = 0.07            # the taper roll-off past the footprint.
POINT_APERTURE_M = 0.03      # the near-POINT read, about 4 px at 6.86 mm.
BLOCK_MS = 20.0              # the block of the moving block bootstrap.
N_BOOT = 400
BOOT_SEED = 12345
BAND_INDEX = 0.05            # the band on a scintillation index (a ratio).
BAND_DB = 0.3                # the band on a mean or a fade quantile, in dB.

BANDS = []


def band(name, value, low, high):
    '''Record one PASS or FAIL line, and print it.'''
    ok = bool(low <= value <= high)
    line = (f"  {'PASS' if ok else 'FAIL'}  {name:48s} {value:9.4f}  "
            f"band [{low:.4f}, {high:.4f}]")
    BANDS.append(line)
    print(line)
    return ok


# ---------------------------------------------------------------------------
# The two campaigns
# ---------------------------------------------------------------------------

def campaigns_root(args):
    '''Give the parent directory of the two campaigns of this study.'''
    root = os.environ.get(ENV_ROOT) or os.path.join(HERE, 'campaigns')
    name = 'rotated_parity_smoke' if args.smoke else 'rotated_parity'
    return os.path.join(root, name)


def campaign_of(arm, args):
    '''Build (or reopen) the campaign of one arm.

    Args:
        arm:  'record' for the rotated frozen-flow record, or 'snapshot' for
              the set of independent trials.
        args: the parsed command line.

    Returns:
        A Campaign.
    '''
    scn, geom = hero_scenario(ELEVATION_DEG)
    spec = None
    if arm == 'record':
        spec = TemporalSpec(dt_s=float(args.dt), n_frames=int(args.frames),
                            strip_dir='unused', record=0,
                            wind_dir_deg=WIND_DIR_DEG, rotated=True,
                            rot_margin=ROT_MARGIN)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        # L0_m stays out: None reads the site outer scale of 25 m (2-P5).
        return Campaign(scn, geom, os.path.join(campaigns_root(args), arm),
                        seed=SEED, preset=args.preset,
                        block_size=int(args.block_size),
                        patch_radius_m=PATCH_RADIUS_M, precision=PRECISION,
                        fft_backend=args.fft_backend, temporal=spec)


def strip_bytes(camp):
    '''Give the bytes that the strips of a campaign hold on disk.'''
    path = os.path.join(camp.root_dir, 'strips')
    if not os.path.isdir(path):
        return 0
    return int(sum(os.path.getsize(os.path.join(path, f))
                   for f in os.listdir(path)
                   if os.path.isfile(os.path.join(path, f))))


def ensure_trials(camp, n_trials, args):
    '''Run the missing trials of one arm, and report the wall time.'''
    missing = max(0, int(n_trials) - int(camp.n_stored))
    if missing == 0:
        print(f'    {camp.n_stored} trials already on disk.', flush=True)
        return 0.0
    t0 = time.perf_counter()
    n_done = camp.run(int(n_trials), workers=args.workers, progress=False)
    wall = time.perf_counter() - t0
    print(f'    {n_done} trials on disk, {wall:.1f} s wall for {missing} new '
          f'({wall / missing:.4f} s/trial)', flush=True)
    return float(wall)


# ---------------------------------------------------------------------------
# The statistics
# ---------------------------------------------------------------------------

def series_of(camp, n_trials):
    '''Give the three power series of one arm.

    Args:
        camp:     the Campaign.
        n_trials: the number of trials to read.

    Returns:
        A dict of the receiver name against its power series. `bucket` is the
        collected power of the full aperture, `smf` is that power times the
        single-mode coupling efficiency, and `point` is the collected power of
        a near-point aperture, which the post-hoc `recollect` gives.
    '''
    result = camp.load(int(n_trials), fields=False)
    bucket = np.array([t.collected_power for t in result.trials], dtype=float)
    eta = np.array([t.smf_eta for t in result.trials], dtype=float)
    point = np.asarray(camp.recollect(aperture_m=POINT_APERTURE_M,
                                      obscuration_ratio=0.0,
                                      n_trials=int(n_trials)), dtype=float)
    return {'point': point, 'bucket': bucket, 'smf': bucket * eta}


def statistics(power):
    '''Give the four statistics of one power series.

    Args:
        power: the received power of each trial, in grid units.

    Returns:
        A dict. `sigma2_I` is the scintillation index var/mean^2 (Andrews and
        Phillips, DOI 10.1117/3.626196, Ch. 8). `mean_db` is 10 log10 of the
        mean. `p5_db` and `p1_db` are the fade depth of the 5 percent and the
        1 percent quantile, in dB UNDER the median of the same series, so they
        do not carry the absolute level.
    '''
    p = np.asarray(power, dtype=float)
    median = float(np.median(p))
    return {'sigma2_I': float(p.var() / p.mean() ** 2),
            'mean_db': float(10.0 * np.log10(p.mean())),
            'p5_db': float(10.0 * np.log10(np.quantile(p, 0.05) / median)),
            'p1_db': float(10.0 * np.log10(np.quantile(p, 0.01) / median))}


def bootstrap_se(power, block, rng, n_boot=N_BOOT):
    '''Give the bootstrap standard error of every statistic of one series.

    A BLOCK of one keeps the ordinary bootstrap, which is right for the
    independent snapshot arm. A block above one takes the MOVING BLOCK
    bootstrap of Kunsch, DOI 10.1214/aos/1176347265: it draws whole runs of
    consecutive frames, so the resample carries the time correlation of the
    record and the bar does not read too small.

    Args:
        power:  the power of each trial.
        block:  the block length, in trials.
        rng:    the resampler.
        n_boot: the number of resamples.

    Returns:
        A dict of the standard error of each statistic.
    '''
    p = np.asarray(power, dtype=float)
    n, b = p.size, max(1, int(block))
    draws = []
    for _ in range(int(n_boot)):
        if b == 1:
            sample = p[rng.integers(0, n, n)]
        else:
            starts = rng.integers(0, n - b + 1, int(np.ceil(n / b)))
            sample = np.concatenate([p[s:s + b] for s in starts])[:n]
        draws.append(statistics(sample))
    return {k: float(np.std([d[k] for d in draws], ddof=1)) for k in draws[0]}


# ---------------------------------------------------------------------------
# The report
# ---------------------------------------------------------------------------

def write_csv(rows):
    '''Write one row for each receiver and statistic. Give the path back.'''
    os.makedirs(DATA, exist_ok=True)
    path = os.path.join(DATA, 'rotated_record_parity.csv')
    fields = ['receiver', 'statistic', 'record', 'record_se', 'snapshot',
              'snapshot_se', 'comparison', 'value', 'two_se']
    with open(path, 'w', newline='', encoding='utf-8') as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    return path


def compare(series_r, series_s, n_frames, dt_s):
    '''Compare the two arms and print every table. Give the CSV rows back.'''
    rng = np.random.default_rng(BOOT_SEED)
    block = max(1, int(round(BLOCK_MS * 1e-3 / float(dt_s))))
    print(f'  the moving block bootstrap holds {block} frames '
          f'({block * dt_s * 1e3:.1f} ms), {N_BOOT} resamples')
    rows = []
    for name in ('point', 'bucket', 'smf'):
        stat_r = statistics(series_r[name])
        stat_s = statistics(series_s[name])
        se_r = bootstrap_se(series_r[name], block, rng)
        se_s = bootstrap_se(series_s[name], 1, rng)
        print('')
        print(f'  receiver: {name}')
        print(f'    {"statistic":12s}{"record":>12s}{"+/- 2 SE":>11s}'
              f'{"snapshot":>12s}{"+/- 2 SE":>11s}{"comparison":>20s}'
              f'{"value":>10s}{"+/- 2 SE":>11s}')
        for key in ('sigma2_I', 'mean_db', 'p5_db', 'p1_db'):
            two_se = 2.0 * float(np.hypot(se_r[key], se_s[key]))
            if key == 'sigma2_I':
                # A RATIO, because an index is a positive scale quantity.
                value = stat_r[key] / stat_s[key]
                two_se = value * 2.0 * float(np.hypot(
                    se_r[key] / stat_r[key], se_s[key] / stat_s[key]))
                what, low, high = ('record / snapshot',
                                   1.0 - BAND_INDEX, 1.0 + BAND_INDEX)
            else:
                # A DIFFERENCE, because a dB is already a logarithm.
                value = stat_r[key] - stat_s[key]
                what, low, high = ('record - snapshot [dB]',
                                   -BAND_DB, BAND_DB)
            print(f'    {key:12s}{stat_r[key]:12.5f}{2 * se_r[key]:11.5f}'
                  f'{stat_s[key]:12.5f}{2 * se_s[key]:11.5f}'
                  f'{what.split("[")[0].strip():>20s}{value:10.4f}'
                  f'{two_se:11.4f}')
            band(f'{name} {key} ({what})', value, low, high)
            rows.append({'receiver': name, 'statistic': key,
                         'record': stat_r[key], 'record_se': se_r[key],
                         'snapshot': stat_s[key], 'snapshot_se': se_s[key],
                         'comparison': what, 'value': value,
                         'two_se': two_se})
    return rows


def _self_check():
    '''Check the two estimators against hand cases. It runs on every call.'''
    # A constant series has no index, no fade and a known mean.
    stat = statistics(np.full(64, 4.0))
    assert abs(stat['sigma2_I']) < 1e-24, stat
    assert abs(stat['mean_db'] - 10.0 * np.log10(4.0)) < 1e-12, stat
    assert abs(stat['p5_db']) < 1e-12 and abs(stat['p1_db']) < 1e-12, stat
    # A lognormal series of a known sigma: sigma2_I = exp(sigma^2) - 1.
    rng = np.random.default_rng(3)
    x = np.exp(0.5 * rng.standard_normal(200000) - 0.125)
    want = np.exp(0.25) - 1.0
    assert abs(statistics(x)['sigma2_I'] / want - 1.0) < 0.02, statistics(x)
    # The BLOCK bootstrap of a correlated series must give a WIDER bar than
    # the ordinary one. A repeated series is the extreme correlated case.
    corr = np.repeat(rng.standard_normal(200) ** 2, 20)
    wide = bootstrap_se(corr, 20, np.random.default_rng(1), n_boot=60)
    tight = bootstrap_se(corr, 1, np.random.default_rng(1), n_boot=60)
    assert wide['mean_db'] > tight['mean_db'], (wide, tight)


def main():
    '''Run the two arms, compare them, and print the bands and the cost.'''
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument('--frames', type=int, default=N_FRAMES,
                        help='the frames of the record and the trials of the '
                             'snapshot arm.')
    parser.add_argument('--dt', type=float, default=DT_S,
                        help='the time step of one frame [s].')
    parser.add_argument('--preset', default=PRESET,
                        help='the sampling preset of both campaigns.')
    parser.add_argument('--block-size', type=int, default=500,
                        help='the trials in one block file.')
    parser.add_argument('--fft-backend', default='cupy',
                        choices=['numpy', 'scipy', 'cupy'],
                        help='the FFT backend of both campaigns.')
    parser.add_argument('--workers', default=None,
                        help='the campaign pool size. Leave it out for the '
                             'CUDA backend: one device runs one stream.')
    parser.add_argument('--analyse', action='store_true',
                        help='read what is stored and compute no trial.')
    parser.add_argument('--keep-strips', action='store_true',
                        help='keep the strips of the finished record.')
    parser.add_argument('--smoke', action='store_true',
                        help='one short pair on the host, to check the script '
                             'end to end.')
    args = parser.parse_args()
    if args.smoke:
        args.frames = 120
        args.preset = 'rapid'
        args.fft_backend = 'numpy'
    if args.workers not in (None, 'auto'):
        args.workers = int(args.workers)
    # A campaign runs WHOLE blocks, so a block larger than the record would run
    # frames past the travel that the strips hold, and a block that does not
    # DIVIDE the frame count runs past the end of the last one.
    args.block_size = min(int(args.block_size), int(args.frames))

    _self_check()
    t_start = time.time()
    print(f'the ROTATED record against independent snapshots, the hero '
          f'downlink at {ELEVATION_DEG:.0f} deg')
    print(f'  {args.preset} preset, {PRECISION} precision, seed {SEED}, '
          f'backend {args.fft_backend}')
    print(f'  {args.frames} frames at {args.dt * 1e3:.2f} ms, crosswind '
          f'{WIND_DIR_DEG:.0f} deg, taper roll-off {ROT_MARGIN:.2f} n')
    print(f'  campaigns under {campaigns_root(args)}')
    print('')

    walls, camps = {}, {}
    for arm in ('record', 'snapshot'):
        camps[arm] = campaign_of(arm, args)
        print(f'  arm {arm}: {camps[arm].root_dir}')
        walls[arm] = (0.0 if args.analyse
                      else ensure_trials(camps[arm], int(args.frames), args))
    n_read = min(int(c.n_stored) for c in camps.values())
    n_read = min(n_read, int(args.frames))
    if n_read == 0:
        print('  INFO  nothing is stored. There is nothing to analyse.')
        return
    if n_read < int(args.frames):
        print(f'  INFO  only {n_read} of {args.frames} trials are on both '
              f'arms. The analysis reads {n_read}.')

    grid = camps['record'].grid
    print('')
    print(f'  grid                     {grid.n:9d} px, '
          f'{grid.pixel_m * 1e3:.2f} mm')
    print(f'  screens                  {camps["record"].plan.z_m.size:9d}')
    print(f'  outer scale              {camps["record"].L0_m:9.1f} m')
    sp_rec = strip_plan(camps['record'].plan, grid,
                        TemporalSpec(dt_s=float(args.dt),
                                     n_frames=int(args.frames),
                                     strip_dir='unused',
                                     wind_dir_deg=WIND_DIR_DEG, rotated=True,
                                     rot_margin=ROT_MARGIN),
                        hero_scenario(ELEVATION_DEG)[1],
                        camps['record'].L0_m)
    print('  layer angles [deg]       '
          + ' '.join(f'{np.rad2deg(t):.0f}' for t in sp_rec.theta))
    print('  crop side of each layer  '
          + ' '.join(str(m) for m in sp_rec.m_crop) + ' px')
    strips = strip_bytes(camps['record'])
    print(f'  strips on disk           {strips / 2 ** 20:9.0f} MB')
    print(f'  record wall              {walls["record"]:9.1f} s')
    print(f'  snapshot wall            {walls["snapshot"]:9.1f} s')
    print(f'  peak working set         '
          f'{peak_working_set_bytes() / 2 ** 30:9.2f} GiB')

    series_r = series_of(camps['record'], n_read)
    series_s = series_of(camps['snapshot'], n_read)
    rows = compare(series_r, series_s, n_read, float(args.dt))

    print('')
    print(f'  file saved: {write_csv(rows)}')
    if not args.keep_strips and camps['record'].n_stored >= int(args.frames):
        path = os.path.join(camps['record'].root_dir, 'strips')
        shutil.rmtree(path, ignore_errors=True)
        print(f'  strips deleted, {strips / 2 ** 20:.0f} MB freed (the seed '
              f'rebuilds them).')

    failed = [line for line in BANDS if 'FAIL' in line]
    print('')
    print(f'  {len(BANDS) - len(failed)} of {len(BANDS)} pass bands hold')
    for line in failed:
        print(line)
    print('')
    print(f'(elapsed {time.time() - t_start:.1f} s)')


if __name__ == '__main__':
    main()
