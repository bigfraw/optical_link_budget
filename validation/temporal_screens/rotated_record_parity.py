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

  R. `--records N` rotated records of `--frames` frames at `--dt`, POOLED.
     Each record has its OWN index, so it has its OWN strips and its OWN
     atmosphere. The wind is a CROSSWIND (`wind_dir_deg = 90`), so the spec
     resolves to the ROTATED route on its own (the default from 2026-09-13),
     with the default taper roll-off `rot_margin = 0.07`.
  S. `N * --frames` independent snapshot trials, on the SAME grid, the same
     screen plan, the same preset and the same seed stream.

WHY SEVERAL RECORDS. One record of 2 s holds only about 75 independent
realisations of the TILT, because the tilt decorrelates in tens of
milliseconds. The fade quantiles of the fibre read that tilt, so their bar
falls only as the number of INDEPENDENT tilt times, not as the number of
frames. N records multiply that count by N at the same strip cost per record.
The script measures the count: see `effective_n`.

WHAT IT COMPARES, for the two receivers of gate (e) (the single-mode fibre and
the bucket) and for a near-POINT aperture:

  - the scintillation index sigma2_I = var(P) / mean(P)^2,
  - the mean power in dB,
  - the 5 percent and the 1 percent quantiles, in dB under the median of that
    arm. That pair IS the fade depth.

THE ERROR BAR IS A BOOTSTRAP of 400 resamples. The frames of a record are
CORRELATED in time, so arm R takes a moving BLOCK bootstrap of about 20 ms
blocks (Kunsch, Ann. Statist. 17(3), pp. 1217 to 1241 (1989),
DOI 10.1214/aos/1176347265), and it draws the blocks WITHIN each record, never
across two records. Arm S takes the ordinary bootstrap.

THE BANDS. 5 percent on a scintillation index, and 0.3 dB on a mean or a fade
quantile. The script prints every number whatever the verdict.

THE COST is reported: the wall time of each arm, the strip bytes on disk and
the peak working set of the process.

Outputs, next to this script:
    data/rotated_record_parity.csv   one row for each receiver and statistic

Run from the repository root. The smoke run takes the host backend:

    python -m validation.temporal_screens.rotated_record_parity --smoke
    python -m validation.temporal_screens.rotated_record_parity
    python -m validation.temporal_screens.rotated_record_parity \
        --dt 2e-3 --frames 1000 --records 5
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
N_RECORDS = 1                # the records that the analysis POOLS.
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


def record_spec(args, record):
    '''Give the TemporalSpec of one record.

    `rotated` stays out: the crosswind resolves the ROTATED route on its own
    (`TemporalSpec.resolve_rotated`, the default from 2026-09-13).

    Args:
        args:   the parsed command line.
        record: the record index, from 0. It gives each record its OWN strips.

    Returns:
        A TemporalSpec.
    '''
    return TemporalSpec(dt_s=float(args.dt), n_frames=int(args.frames),
                        strip_dir='unused', record=int(record),
                        wind_dir_deg=WIND_DIR_DEG, rot_margin=ROT_MARGIN)


def campaign_of(arm, args, record=0):
    '''Build (or reopen) the campaign of one arm.

    The record arms hold the time step in their directory name, because two
    steps are two different atmospheres of the same scenario and a campaign
    directory holds ONE manifest. The snapshot arm has no time axis, so its
    directory is shared by every step and a longer run only adds trials.

    Args:
        arm:    'record' for a rotated frozen-flow record, or 'snapshot' for
                the set of independent trials.
        args:   the parsed command line.
        record: the record index, for the record arm.

    Returns:
        A Campaign.
    '''
    scn, geom = hero_scenario(ELEVATION_DEG)
    spec = None
    name = arm
    if arm == 'record':
        spec = record_spec(args, record)
        name = f'record{int(record):02d}_dt{int(round(args.dt * 1e6))}us'
    n_trials = int(args.frames) * (1 if arm == 'record' else int(args.records))
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        # L0_m stays out: None reads the site outer scale of 25 m (2-P5).
        return Campaign(scn, geom, os.path.join(campaigns_root(args), name),
                        seed=SEED, preset=args.preset,
                        block_size=min(int(args.block_size), n_trials),
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


def bootstrap_se(segments, block, rng, n_boot=N_BOOT):
    '''Give the bootstrap standard error of every statistic of one arm.

    A BLOCK of one keeps the ordinary bootstrap, which is right for the
    independent snapshot arm. A block above one takes the MOVING BLOCK
    bootstrap of Kunsch, DOI 10.1214/aos/1176347265: it draws whole runs of
    consecutive frames, so the resample carries the time correlation of the
    record and the bar does not read too small.

    THE BLOCKS STAY INSIDE ONE RECORD. Two records are two atmospheres, so a
    block that crossed the join would join two unrelated frames. Each segment
    is resampled to its own length and the pieces are then joined, so every
    resample has the length of the pooled series.

    Args:
        segments: a list of the power series, ONE for each record. The
                  snapshot arm gives one segment.
        block:    the block length, in trials.
        rng:      the resampler.
        n_boot:   the number of resamples.

    Returns:
        A dict of the standard error of each statistic.
    '''
    segs = [np.asarray(s, dtype=float) for s in segments]
    b = max(1, int(block))
    draws = []
    for _ in range(int(n_boot)):
        parts = []
        for p in segs:
            n = p.size
            if b == 1 or n < 2 * b:
                parts.append(p[rng.integers(0, n, n)])
            else:
                starts = rng.integers(0, n - b + 1, int(np.ceil(n / b)))
                parts.append(
                    np.concatenate([p[s:s + b] for s in starts])[:n])
        draws.append(statistics(np.concatenate(parts)))
    return {k: float(np.std([d[k] for d in draws], ddof=1)) for k in draws[0]}


def effective_n(x):
    '''Give the number of INDEPENDENT samples that one time series holds.

    N_eff = N / (1 + 2 sum_k rho_k), with rho the autocorrelation and the sum
    cut at the first lag where rho goes below zero (the initial positive
    sequence of Geyer, Statist. Sci. 7(4), pp. 473 to 483 (1992),
    DOI 10.1214/ss/1177011137). For white noise rho is 0 and N_eff is N.

    Args:
        x: the series of one record.

    Returns:
        A float.
    '''
    p = np.asarray(x, dtype=float)
    d = p - p.mean()
    n = d.size
    var = float(d @ d)
    if var <= 0.0:
        return float(n)
    # The autocorrelation through the FFT: it is the same sum, and it costs
    # n log n instead of n^2 (Schmidt, DOI 10.1117/3.866274, Ch. 2).
    f = np.fft.rfft(d, 2 * n)
    rho = np.fft.irfft(f * np.conj(f), 2 * n)[:n] / var
    cut = np.flatnonzero(rho[1:] <= 0.0)
    k = int(cut[0]) if cut.size else n - 1
    return float(n / (1.0 + 2.0 * float(rho[1:1 + k].sum())))


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


def compare(segments_r, series_s, dt_s):
    '''Compare the two arms and print every table. Give the CSV rows back.

    Args:
        segments_r: a list of the series dict of each record of arm R.
        series_s:   the series dict of arm S.
        dt_s:       the time step of a frame, in s.

    Returns:
        The CSV rows.
    '''
    rng = np.random.default_rng(BOOT_SEED)
    block = max(1, int(round(BLOCK_MS * 1e-3 / float(dt_s))))
    print(f'  the moving block bootstrap holds {block} frames '
          f'({block * dt_s * 1e3:.1f} ms), {N_BOOT} resamples, drawn WITHIN '
          f'each of the {len(segments_r)} records')
    rows = []
    for name in ('point', 'bucket', 'smf'):
        segs_r = [s[name] for s in segments_r]
        stat_r = statistics(np.concatenate(segs_r))
        stat_s = statistics(series_s[name])
        se_r = bootstrap_se(segs_r, block, rng)
        se_s = bootstrap_se([series_s[name]], 1, rng)
        n_eff = sum(effective_n(s) for s in segs_r)
        print('')
        print(f'  independent tilt times in the pooled record: '
              f'{n_eff:.0f} of {sum(s.size for s in segs_r)} frames '
              f'({n_eff / len(segs_r):.0f} for each record)')
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
    wide = bootstrap_se([corr], 20, np.random.default_rng(1), n_boot=60)
    tight = bootstrap_se([corr], 1, np.random.default_rng(1), n_boot=60)
    assert wide['mean_db'] > tight['mean_db'], (wide, tight)
    # Two segments give a resample of the POOLED length.
    two = bootstrap_se([corr, corr], 20, np.random.default_rng(1), n_boot=8)
    assert set(two) == set(wide), (two, wide)
    # `effective_n` reads N for white noise and about N/20 for the series
    # above, which repeats every value 20 times.
    white = rng.standard_normal(20000)
    assert 0.8 < effective_n(white) / white.size < 1.2, effective_n(white)
    assert 0.5 < effective_n(corr) / (corr.size / 20.0) < 2.0, \
        effective_n(corr)


def main():
    '''Run the two arms, compare them, and print the bands and the cost.'''
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument('--frames', type=int, default=N_FRAMES,
                        help='the frames of the record and the trials of the '
                             'snapshot arm.')
    parser.add_argument('--dt', type=float, default=DT_S,
                        help='the time step of one frame [s].')
    parser.add_argument('--records', type=int, default=N_RECORDS,
                        help='the number of records to POOL. Each one holds '
                             'its own strips and its own atmosphere.')
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
        args.records = 2
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
    print(f'  {args.records} record(s) of {args.frames} frames at '
          f'{args.dt * 1e3:.2f} ms = '
          f'{args.records * args.frames * args.dt:.1f} s of atmosphere')
    print(f'  crosswind {WIND_DIR_DEG:.0f} deg, taper roll-off '
          f'{ROT_MARGIN:.2f} n, route '
          f'{"rotated" if record_spec(args, 0).resolve_rotated() else "box"} '
          f'(auto)')
    print(f'  campaigns under {campaigns_root(args)}')
    print('')

    n_snap = int(args.frames) * int(args.records)
    camps_r, walls_r, builds_r, strips_r = [], [], [], []
    for r in range(int(args.records)):
        camp = campaign_of('record', args, record=r)
        camps_r.append(camp)
        print(f'  arm record {r}: {camp.root_dir}')
        if args.analyse:
            walls_r.append(0.0)
            builds_r.append(0.0)
        else:
            # The strip build is idempotent, so this call TIMES it and the
            # runner then skips it.
            t_b = time.perf_counter()
            camp._build_strips()
            builds_r.append(time.perf_counter() - t_b)
            print(f'    strips built in {builds_r[-1]:.1f} s', flush=True)
            walls_r.append(ensure_trials(camp, int(args.frames), args))
        strips_r.append(strip_bytes(camp))
    camp_s = campaign_of('snapshot', args)
    print(f'  arm snapshot: {camp_s.root_dir}')
    wall_s = (0.0 if args.analyse else ensure_trials(camp_s, n_snap, args))

    n_read = min(int(c.n_stored) for c in camps_r)
    n_read = min(n_read, int(args.frames),
                 int(camp_s.n_stored) // int(args.records))
    if n_read == 0:
        print('  INFO  nothing is stored. There is nothing to analyse.')
        return
    if n_read < int(args.frames):
        print(f'  INFO  only {n_read} of {args.frames} frames are on every '
              f'arm. The analysis reads {n_read} for each record.')

    grid = camps_r[0].grid
    print('')
    print(f'  grid                     {grid.n:9d} px, '
          f'{grid.pixel_m * 1e3:.2f} mm')
    print(f'  screens                  {camps_r[0].plan.z_m.size:9d}')
    print(f'  outer scale              {camps_r[0].L0_m:9.1f} m')
    sp_rec = strip_plan(camps_r[0].plan, grid, record_spec(args, 0),
                        hero_scenario(ELEVATION_DEG)[1], camps_r[0].L0_m)
    print('  layer angles [deg]       '
          + ' '.join(f'{np.rad2deg(t):.0f}' for t in sp_rec.theta))
    print('  crop side of each layer  '
          + ' '.join(str(m) for m in sp_rec.m_crop) + ' px')
    print(f'  strips on disk           {sum(strips_r) / 2 ** 20:9.0f} MB '
          f'({sum(strips_r) / 2 ** 20 / len(strips_r):.0f} MB for each '
          f'record)')
    print(f'  strip build              {np.mean(builds_r):9.1f} s for each '
          f'record')
    print(f'  record wall              {sum(walls_r):9.1f} s '
          f'({sum(walls_r) / max(1, n_read * len(walls_r)):.4f} s/frame)')
    print(f'  snapshot wall            {wall_s:9.1f} s '
          f'({wall_s / max(1, n_snap):.4f} s/trial)')
    print(f'  peak working set         '
          f'{peak_working_set_bytes() / 2 ** 30:9.2f} GiB')

    segments_r = [series_of(c, n_read) for c in camps_r]
    series_s = series_of(camp_s, n_read * int(args.records))
    rows = compare(segments_r, series_s, float(args.dt))

    print('')
    print(f'  file saved: {write_csv(rows)}')
    if not args.keep_strips:
        for camp, size in zip(camps_r, strips_r):
            if int(camp.n_stored) >= int(args.frames):
                shutil.rmtree(os.path.join(camp.root_dir, 'strips'),
                              ignore_errors=True)
        print(f'  strips deleted, {sum(strips_r) / 2 ** 20:.0f} MB freed (the '
              f'seed rebuilds them).')

    failed = [line for line in BANDS if 'FAIL' in line]
    print('')
    print(f'  {len(BANDS) - len(failed)} of {len(BANDS)} pass bands hold')
    for line in failed:
        print(line)
    print('')
    print(f'(elapsed {time.time() - t_start:.1f} s)')


if __name__ == '__main__':
    main()
