'''
Gate (c) of the temporal strip screens: frame 0 IS a drawn screen stack.

THIS IS THE GATE THAT PROVES THE ROUTE. A frozen-flow record crops its frames
out of ONE long STRIP screen for each layer (see
`olb/waveoptics/turbulence/temporal.py`). A strip is a Fourier screen of a
RECTANGULAR grid, so its low-frequency content is built on a different
frequency grid from the SQUARE screen that a snapshot trial draws. If that
difference moved the turbulence, every number of a temporal record would be
wrong from the first frame.

THE MEASUREMENT. The script runs the hero downlink at 30 deg two ways, with the
SAME base seed:

  1. N independent SNAPSHOT trials, the run of record.
  2. N temporal RECORDS of 2 frames each, and it keeps FRAME 0 of each one.
     Record r draws its strips from the seed of trial r, so the two arms sit on
     the same seed stream.

Two statistics then have to agree: the APERTURE scintillation index of the
collected power, and the MEAN single-mode fibre coupling efficiency. The two
arms are NOT bit-identical and they cannot be: a rectangular draw and a square
draw of one seed give different random fields. So the band is statistical.

THERE ARE TWO BANDS FOR EACH STATISTIC. The first is the fixed 5 percent on the
ratio. The second is 2 standard errors of the difference, from a bootstrap.
THE SECOND IS THE MEANINGFUL ONE at a few hundred trials: the fibre coupling
and the aperture index are both heavy-tailed, so 256 trials give a 2 SE bar of
about 20 to 25 percent and the fixed 5 percent band cannot be resolved. The
script prints both, so the reader sees the size of the error bar next to the
fixed band.

WHAT THIS GATE FOUND (2026-09-13, 256 trials of each arm, `rapid`, 30 deg): the
index ratio is 0.939 with a 2 SE bar of 0.251, and the coupling ratio is 1.039
with a 2 SE bar of 0.215. Both sit well inside the error bar, so a strip crop
is the same atmosphere as a drawn screen stack. A run of 128 trials of the same
seed stream read 0.64 on the index; that was the tail of the index estimator,
not a deficit, and the 2 SE bar said so at the time.

THE PRESET IS `rapid`. The gate is a STATISTICAL comparison of two arms that
share the preset, so the preset only sets the cost. `rapid` keeps the run near
five minutes on eight cores. The physics of the preset (5 screens, 256 px) is
the physics of both arms.

THE STRIPS ARE DELETED after each record. One record of this case holds 25 MB
of strips, and 256 records would hold 6.4 GB. The strips are a cache that the
seed rebuilds, so nothing is lost.

Sources:
- Taylor, Proc. R. Soc. Lond. A 164, pp. 476 to 490 (1938),
  DOI 10.1098/rspa.1938.0032. The frozen-flow hypothesis.
- Schmidt, Numerical Simulation of Optical Wave Propagation with Examples in
  MATLAB, DOI 10.1117/3.866274, Ch. 9, Eqs. (9.78) to (9.81), printed pp. 166
  to 169. The Fourier screen and its subharmonics.
- Andrews and Phillips, Laser Beam Propagation through Random Media, 2nd ed.,
  DOI 10.1117/3.626196, Ch. 8. The scintillation index of an aperture.

Outputs, next to this script:
    data/frame0_parity.csv      the two scalars of every trial of both arms
    figures/frame0_parity.png   the two distributions, as CDFs

Run from the repo root:
    python -m validation.temporal_screens.frame0_parity
'''

import argparse
import csv
import os
import shutil
import time
import warnings

import numpy as np

from olb.waveoptics.threader import Threader
from olb.waveoptics.turbulence.run import propagate_turbulent_scenario
from olb.waveoptics.turbulence.temporal import TemporalSpec
from validation.waveoptics_ao.waveoptics_ao import hero_scenario

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, 'data')
FIGS = os.path.join(HERE, 'figures')

# ---- the case ----
ELEVATION_DEG = 30.0
PRESET = 'rapid'
PRECISION = 'single'
GENERATOR = 'olb'
SEED = 20260913
DT_S = 5e-4                  # the step of the two-frame probe record.
BAND = 0.05                  # the pass band on each ratio.
N_BOOT = 2000                # the bootstrap resamples of the error bar.

BANDS = []


def band(name, value, low, high):
    '''Record one PASS or FAIL line, and print it.'''
    ok = low <= value <= high
    line = (f"  {'PASS' if ok else 'FAIL'}  {name:46s} {value:9.4f}  "
            f"band [{low:.4f}, {high:.4f}]")
    BANDS.append(line)
    print(line)
    return ok


def scint_index(power):
    '''Give the scintillation index sigma2_I = var(P) / mean(P)^2.

    Source: Andrews and Phillips, DOI 10.1117/3.626196, Ch. 8. The aperture
    index of a collected power.
    '''
    p = np.asarray(power, dtype=float)
    return float(p.var() / p.mean() ** 2)


def boot_se(values, statistic, rng, n_boot=N_BOOT):
    '''Give the bootstrap standard error of one statistic.'''
    v = np.asarray(values, dtype=float)
    idx = rng.integers(0, v.size, size=(int(n_boot), v.size))
    return float(np.std([statistic(v[row]) for row in idx], ddof=1))


def run_snapshots(scn, geom, n_trials):
    '''Run the independent snapshot arm.

    Returns:
        The pair (collected_power, smf_eta), two float arrays.
    '''
    res = propagate_turbulent_scenario(
        scn, geom, n_trials=n_trials, seed=SEED, preset=PRESET,
        precision=PRECISION, screen_generator=GENERATOR, threader=Threader())
    return (np.array([t.collected_power for t in res.trials]),
            np.array([t.smf_eta for t in res.trials]))


def run_frame0(scn, geom, n_records, strip_root):
    '''Run the frame-0 arm: one two-frame record for each seed index.

    The strips of a record are removed as soon as that record is done, so the
    disk holds one record at a time.

    Returns:
        The pair (collected_power, smf_eta), two float arrays.
    '''
    power, eta = [], []
    for r in range(int(n_records)):
        spec = TemporalSpec(dt_s=DT_S, n_frames=2, strip_dir=strip_root,
                            record=r)
        res = propagate_turbulent_scenario(
            scn, geom, n_trials=1, seed=SEED, preset=PRESET,
            precision=PRECISION, screen_generator=GENERATOR, temporal=spec)
        power.append(res.trials[0].collected_power)
        eta.append(res.trials[0].smf_eta)
        shutil.rmtree(strip_root, ignore_errors=True)
    return np.asarray(power, dtype=float), np.asarray(eta, dtype=float)


def write_csv(snap, frame):
    '''Write the per-trial scalars of both arms. Give the path back.'''
    os.makedirs(DATA, exist_ok=True)
    path = os.path.join(DATA, 'frame0_parity.csv')
    with open(path, 'w', newline='', encoding='utf-8') as fh:
        w = csv.writer(fh)
        w.writerow(['trial', 'snapshot_power', 'snapshot_smf_eta',
                    'frame0_power', 'frame0_smf_eta'])
        for k in range(snap[0].size):
            w.writerow([k, snap[0][k], snap[1][k], frame[0][k], frame[1][k]])
    return path


def draw_cdf(snap, frame):
    '''Draw the two distributions of both arms. Give the path back.'''
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    os.makedirs(FIGS, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 4.0))
    titles = ('collected power / vacuum', 'SMF coupling efficiency')
    for ax, i, title in zip(axes, (0, 1), titles):
        for values, label, style in ((snap[i], 'snapshot', '-'),
                                     (frame[i], 'frame 0 of a record', '--')):
            v = np.sort(np.asarray(values, dtype=float))
            ax.step(v, np.arange(1, v.size + 1) / v.size, style, label=label)
        ax.set_xlabel(title)
        ax.set_ylabel('cumulative fraction')
        ax.grid(True, alpha=0.3)
        ax.legend(loc='lower right', fontsize=8)
    fig.suptitle(f'gate (c): frame 0 against a drawn snapshot, hero downlink '
                 f'{ELEVATION_DEG:.0f} deg, {PRESET} preset')
    fig.tight_layout()
    path = os.path.join(FIGS, 'frame0_parity.png')
    fig.savefig(path, dpi=140)
    plt.close(fig)
    return path


def main():
    '''Run the gate, print the bands, and write the outputs.'''
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--trials', type=int, default=256,
                        help='the trials of each arm (default 256).')
    args = parser.parse_args()

    t_start = time.time()
    scn, geom = hero_scenario(ELEVATION_DEG)
    strip_root = os.path.join(HERE, 'strips_frame0')
    shutil.rmtree(strip_root, ignore_errors=True)

    print(f'gate (c), the hero downlink at {ELEVATION_DEG:.0f} deg, '
          f'{PRESET} preset, {PRECISION} precision, seed {SEED}')
    print(f'  trials of each arm      {args.trials:11d}')
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        t0 = time.time()
        snap = run_snapshots(scn, geom, args.trials)
        snap_s = time.time() - t0
        t0 = time.time()
        try:
            frame = run_frame0(scn, geom, args.trials, strip_root)
        finally:
            shutil.rmtree(strip_root, ignore_errors=True)
        frame_s = time.time() - t0

    print(f'  snapshot arm            {snap_s:11.1f} s')
    print(f'  frame-0 arm             {frame_s:11.1f} s')
    print('')

    rng = np.random.default_rng(7)
    print('1. the aperture scintillation index of the collected power')
    s_snap, s_frame = scint_index(snap[0]), scint_index(frame[0])
    se = np.hypot(boot_se(snap[0], scint_index, rng),
                  boot_se(frame[0], scint_index, rng))
    print(f'  snapshot     sigma2_I = {s_snap:.4f}')
    print(f'  frame 0      sigma2_I = {s_frame:.4f}')
    print(f'  2 SE of the difference  {2 * se / s_snap:.4f} '
          f'(relative to the snapshot arm)')
    r_s = 2 * se / s_snap
    band('sigma2_I, frame 0 / snapshot', s_frame / s_snap,
         1.0 - BAND, 1.0 + BAND)
    band('sigma2_I, inside 2 SE of 1', s_frame / s_snap, 1.0 - r_s, 1.0 + r_s)

    print('')
    print('2. the mean single-mode fibre coupling efficiency')
    m_snap, m_frame = float(snap[1].mean()), float(frame[1].mean())
    se_m = np.hypot(snap[1].std(ddof=1), frame[1].std(ddof=1)) / \
        np.sqrt(args.trials)
    print(f'  snapshot     <eta>    = {m_snap:.6f}')
    print(f'  frame 0      <eta>    = {m_frame:.6f}')
    print(f'  2 SE of the difference  {2 * se_m / m_snap:.4f} '
          f'(relative to the snapshot arm)')
    r_m = 2 * se_m / m_snap
    band('<smf eta>, frame 0 / snapshot', m_frame / m_snap,
         1.0 - BAND, 1.0 + BAND)
    band('<smf eta>, inside 2 SE of 1', m_frame / m_snap, 1.0 - r_m, 1.0 + r_m)

    print('')
    print(f'  file saved: {write_csv(snap, frame)}')
    print(f'  figure saved: {draw_cdf(snap, frame)}')
    print('    Caption: the cumulative distribution of the collected power and '
          'of the')
    print('    fibre coupling, for an independent snapshot and for frame 0 of '
          'a record.')

    failed = [line for line in BANDS if 'FAIL' in line]
    print('')
    print(f'  {len(BANDS) - len(failed)} of {len(BANDS)} pass bands hold')
    for line in failed:
        print(line)
    print('')
    print(f'(elapsed {time.time() - t_start:.1f} s)')


if __name__ == '__main__':
    main()
