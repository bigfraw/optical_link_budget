'''
Arm 6: does the recursion KEEP a subharmonic start on a narrow frame?

THE QUESTION. Q7 shows that a subharmonic initial frame removes the optical
spin-up at row 0. Q6 shows that the recursion has its OWN stationary level on
a frame of 0.2 outer scales: the extrusion axis over-correlates at lags of
0.1 to 1 L0 at every `n_columns`. So a correct start may drift toward that
level as rows are added. This arm reads the drift, at the production cell of
Q6, from BOTH starts and with NO spin-up, so the two are comparable row for
row.

THE METHOD. N = 512, dx = 1 cm, r0 = 10 cm, L0 = 25 m (the frame is 0.2 L0),
`n_columns` in {2, 8}, two starts (the stock plain Fourier frame, and the olb
subharmonic frame through `subharmonic_start.SubharmonicStart`), 8 seeds.
Each seed extrudes 10000 rows (4 L0, the same record as the Q6 cell). The
reads are:

  1. rho(k) along the extrusion axis over the whole record, through the Q6
     estimator (`extrusion_stationarity.row_lag_sums` in column blocks), so
     the number compares with the Q6 excess of +0.228 (2 columns) and +0.177
     (8 columns) at 0.5 L0, which came from a plain start with a 5000-row
     spin-up.
  2. at rows 0, 512, 2048, 5000 and 10000: the piston-removed frame variance
     against B(0) minus the piston theory, and the Z-tilt over a 1.0 m pupil
     against the Noll filter integral (Noll, DOI 10.1364/JOSA.66.000207).

THE READING. If the rho excess is the same from both starts, the Q6 defect is
the frame width and the spin-up is decoupled from it. If the checkpoints of
the subharmonic start fall toward the plain level, the number of rows over
which they fall is the horizon a parallel worker gets from its start.

Sources: Assemat and Wilson, DOI 10.1364/OE.14.000988 (the recursion and its
covariance); Lane, Glindemann and Dainty, DOI 10.1088/0959-7174/2/3/003 (the
subharmonics).

Outputs, next to this script:
    data/kept_start_rho.csv      start, n_columns, k, r_m, rho_meas, rho_theory
    data/kept_start_frames.csv   start, n_columns, rows, metric, ratio
    figures/kept_start.png       rho(k) per start and the two checkpoint ratios

Run from the repository root (about 8 minutes):
    python -m validation.screens.kept_start
'''

import csv
import io
import contextlib
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '..', '..'))

from olb.waveoptics.schmidt.turbulence import (             # noqa: E402
    von_karman_phase_psd)
from validation.screens import helpers                     # noqa: E402
from validation.screens import extrusion_stationarity as es  # noqa: E402
from validation.screens.infinite_screen_stats import (     # noqa: E402
    N, DX_M, R0_M, L0_WIDE_M, WAVELENGTH_M, K_RAD_M, PUPIL_D_M,
    grid_piston_variance)
from validation.screens.subharmonic_start import SubharmonicStart  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
RHO_CSV = os.path.join(HERE, 'data', 'kept_start_rho.csv')
FRAMES_CSV = os.path.join(HERE, 'data', 'kept_start_frames.csv')
PNG_PATH = os.path.join(HERE, 'figures', 'kept_start.png')

N_SEEDS = 32                    # 8 in the first run (2026-09-08); 32 gives the error bars
RECORD_ROWS = 10000
CHECKPOINTS = (0, 512, 2048, 5000, 10000)
STARTS = ('plain', 'olb')
N_COLUMNS = (2, 8)
MASTER_SEED = 20260909
K_MAX = 5000                    # 2 L0 of lag
COLUMN_BLOCK = 64               # the FFT block of the row-lag estimator
Q6_EXCESS = {2: 0.228, 8: 0.177}  # the Q6 plain-start excess at 0.5 L0


def make_screen(start, n_columns, seed):
    '''Build one extruded screen with the named start and stencil depth.'''
    klass = type(f'Start_{start}', (SubharmonicStart,), {'start': start})
    with contextlib.redirect_stdout(io.StringIO()):
        return klass(N, DX_M, R0_M, L0_WIDE_M, random_seed=int(seed),
                     n_columns=n_columns)


def row_lag(rows, k_max):
    '''The Q6 route: raise the module K_MAX, sum over column blocks.'''
    keep = es.K_MAX
    es.K_MAX = int(k_max)
    try:
        total = np.zeros(k_max + 1)
        n_blocks = 0
        for start in range(0, rows.shape[1], COLUMN_BLOCK):
            total += es.row_lag_sums(rows[:, start:start + COLUMN_BLOCK])
            n_blocks += 1
    finally:
        es.K_MAX = keep
    m = rows.shape[0]
    pairs = m - np.arange(k_max + 1)
    return total / n_blocks / pairs


def main():
    t0 = time.time()
    mask = helpers.pupil_mask(N, DX_M, PUPIL_D_M)
    d_eff = helpers.mask_diameter(mask, DX_M)
    psd = lambda f: von_karman_phase_psd(f, R0_M, L0_WIDE_M)   # noqa: E731
    z_theory = helpers.tilt_filter_variance(psd, d_eff) / (
        K_RAD_M * d_eff / 4.0) ** 2
    b0 = float(helpers.vk_covariance_closed(np.array([0.0]), R0_M,
                                            L0_WIDE_M)[0])
    var_np_theory = b0 - grid_piston_variance(L0_WIDE_M)
    lags = np.arange(K_MAX + 1)
    rho_theory = helpers.vk_covariance_closed(lags * DX_M, R0_M,
                                              L0_WIDE_M) / b0
    seeds = helpers.spawn_seeds(MASTER_SEED, N_SEEDS)
    k_half = int(round(0.5 * L0_WIDE_M / DX_M))

    print('ARM 6. IS A SUBHARMONIC START KEPT ON A 0.2 L0 FRAME?')
    print(f'   N = {N}, L0 = {L0_WIDE_M:.0f} m, r0 = {R0_M * 1e2:.0f} cm, '
          f'{N_SEEDS} seeds, {RECORD_ROWS} rows each, no spin-up')
    print('')

    rho, rho_se = {}, {}
    frames, frames_se = {}, {}
    rho_rows, frame_rows = [], []
    for n_col in N_COLUMNS:
        for start in STARTS:
            t1 = time.time()
            rho_seeds = []
            var_np = {k: [] for k in CHECKPOINTS}
            angles = {k: [] for k in CHECKPOINTS}
            for seed in seeds:
                screen = make_screen(start, n_col, seed)
                record = np.empty((RECORD_ROWS, N))
                done = 0
                for k in CHECKPOINTS:
                    while done < k:
                        screen.add_row()
                        record[done] = screen._scrn[0]
                        done += 1
                    frame = screen.scrn
                    var_np[k].append(float(np.mean((frame - frame.mean())
                                                   ** 2)))
                    zx, zy, _, _ = helpers.zernike_tilt(frame, mask, DX_M,
                                                        WAVELENGTH_M)
                    angles[k].extend([zx, zy])
                c = row_lag(record, K_MAX)
                rho_seeds.append(c / c[0])
            rho_seeds = np.array(rho_seeds)
            # The per-seed rho is one sample of the record statistic, so the
            # standard error of the mean over the seeds is the error bar.
            rho[(start, n_col)] = rho_seeds.mean(axis=0)
            rho_se[(start, n_col)] = rho_seeds.std(axis=0, ddof=1) / np.sqrt(
                N_SEEDS)
            for k in CHECKPOINTS:
                v = np.array(var_np[k]) / var_np_theory
                a = np.array(angles[k]) ** 2 / z_theory
                frames[(start, n_col, k, 'nopiston_var')] = float(v.mean())
                frames_se[(start, n_col, k, 'nopiston_var')] = float(
                    v.std(ddof=1) / np.sqrt(v.size))
                frames[(start, n_col, k, 'ztilt_var')] = float(a.mean())
                frames_se[(start, n_col, k, 'ztilt_var')] = float(
                    a.std(ddof=1) / np.sqrt(a.size))
                for metric in ('nopiston_var', 'ztilt_var'):
                    frame_rows.append((start, n_col, k, metric,
                                       frames[(start, n_col, k, metric)],
                                       frames_se[(start, n_col, k, metric)]))
            for k in lags[::25]:
                rho_rows.append((start, n_col, int(k), k * DX_M,
                                 rho[(start, n_col)][k],
                                 rho_se[(start, n_col)][k], rho_theory[k]))
            print(f'   n_columns {n_col}, start {start:<6} done in '
                  f'{time.time() - t1:6.1f} s')

    print('')
    print('1. THE EXTRUSION-AXIS rho EXCESS AT 0.5 L0 (theory '
          f'{rho_theory[k_half]:.3f})')
    print(f'   {"n_col":>6}{"start":>8}{"rho":>9}{"excess":>9}{"2 SE":>8}'
          f'{"Q6 plain, spun up":>20}')
    for n_col in N_COLUMNS:
        for start in STARTS:
            r = rho[(start, n_col)][k_half]
            print(f'   {n_col:>6}{start:>8}{r:>9.3f}'
                  f'{r - rho_theory[k_half]:>+9.3f}'
                  f'{2 * rho_se[(start, n_col)][k_half]:>8.3f}'
                  f'{Q6_EXCESS[n_col]:>+20.3f}')
    print('')
    print('   THE TWO STARTS AGAINST EACH OTHER at 0.1, 0.5 and 1 L0: the '
          'difference and its 2 SE')
    for n_col in N_COLUMNS:
        for frac in (0.1, 0.5, 1.0):
            k = int(round(frac * L0_WIDE_M / DX_M))
            d = rho[('olb', n_col)][k] - rho[('plain', n_col)][k]
            se = np.hypot(rho_se[('olb', n_col)][k],
                          rho_se[('plain', n_col)][k])
            print(f'   n_columns {n_col}, lag {frac:.1f} L0: olb - plain = '
                  f'{d:+.3f} +- {2 * se:.3f}; olb - theory = '
                  f'{rho[("olb", n_col)][k] - rho_theory[k]:+.3f} +- '
                  f'{2 * rho_se[("olb", n_col)][k]:.3f}')
    print('')
    print('2. THE CHECKPOINTS, ratio to theory, with 2 SE')
    print(f'   {"n_col":>6}{"start":>8}{"rows":>7}{"no-piston var":>15}'
          f'{"2 SE":>8}{"Z-tilt":>9}{"2 SE":>8}')
    for n_col in N_COLUMNS:
        for start in STARTS:
            for k in CHECKPOINTS:
                print(f'   {n_col:>6}{start:>8}{k:>7}'
                      f'{frames[(start, n_col, k, "nopiston_var")]:>15.3f}'
                      f'{2 * frames_se[(start, n_col, k, "nopiston_var")]:>8.3f}'
                      f'{frames[(start, n_col, k, "ztilt_var")]:>9.3f}'
                      f'{2 * frames_se[(start, n_col, k, "ztilt_var")]:>8.3f}')

    os.makedirs(os.path.dirname(RHO_CSV), exist_ok=True)
    with open(RHO_CSV, 'w', newline='') as handle:
        w = csv.writer(handle)
        w.writerow(['start', 'n_columns', 'k', 'r_m', 'rho_meas', 'rho_se',
                    'rho_theory'])
        w.writerows(rho_rows)
    with open(FRAMES_CSV, 'w', newline='') as handle:
        w = csv.writer(handle)
        w.writerow(['start', 'n_columns', 'rows', 'metric', 'ratio',
                    'ratio_se'])
        w.writerows(frame_rows)

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.5))
    ax = axes[0]
    ax.plot(lags[1:] * DX_M / L0_WIDE_M, rho_theory[1:], 'k', lw=2.5,
            label='von Karman theory')
    for n_col in N_COLUMNS:
        for start in STARTS:
            x = lags[1:] * DX_M / L0_WIDE_M
            line, = ax.plot(x, rho[(start, n_col)][1:],
                            '--' if start == 'plain' else '-',
                            label=f'{start}, n_columns {n_col}')
            ax.fill_between(x,
                            rho[(start, n_col)][1:]
                            - 2 * rho_se[(start, n_col)][1:],
                            rho[(start, n_col)][1:]
                            + 2 * rho_se[(start, n_col)][1:],
                            color=line.get_color(), alpha=0.15)
    ax.set_xscale('log')
    ax.set_xlabel('row lag, outer scales')
    ax.set_ylabel('rho(k)')
    ax.set_title('extrusion axis, no spin-up, 10000 rows')
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    for ax, metric, title in zip(axes[1:], ('nopiston_var', 'ztilt_var'),
                                 ('piston-removed variance / theory',
                                  'Z-tilt over 1 m / Noll')):
        for n_col in N_COLUMNS:
            for start in STARTS:
                ax.errorbar([max(k, 1) for k in CHECKPOINTS],
                            [frames[(start, n_col, k, metric)]
                             for k in CHECKPOINTS],
                            yerr=[2 * frames_se[(start, n_col, k, metric)]
                                  for k in CHECKPOINTS],
                            fmt='o--' if start == 'plain' else 'o-',
                            capsize=3,
                            label=f'{start}, n_columns {n_col}')
        ax.axhline(1.0, color='k', lw=1)
        ax.set_xscale('log')
        ax.set_xlabel('add_row calls (row 0 drawn at 1)')
        ax.set_title(title, fontsize=10)
        ax.grid(alpha=0.3)
    axes[1].legend(fontsize=8)
    fig.suptitle(f'IS THE START KEPT? N {N}, L0 {L0_WIDE_M:.0f} m (frame '
                 f'{N * DX_M / L0_WIDE_M:.2f} L0), {N_SEEDS} seeds')
    fig.tight_layout()
    os.makedirs(os.path.dirname(PNG_PATH), exist_ok=True)
    fig.savefig(PNG_PATH, dpi=150)
    print(f'   figure: {PNG_PATH}')
    print(f'   total {time.time() - t0:.0f} s')


if __name__ == '__main__':
    main()
