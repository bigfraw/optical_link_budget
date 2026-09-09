'''
Arm 5: does a SUBHARMONIC initial frame remove the spin-up of the extrusion?

THE QUESTION. `PhaseScreenVonKarman.make_initial_screen` starts the
extrusion from a plain Fourier screen (`ft_phase_screen`, no subharmonics).
That frame holds no power below the grid fundamental, so its variance, its
piston and its tilt are low, and the recursion needs about 2 L0 / dx rows to
forget it (FINDINGS Q4 and Q5 point 2; 5000 rows at L0 = 25 m and dx = 1 cm).
The recursion is Markov in the last `n_columns` rows (Assemat and Wilson,
DOI 10.1364/OE.14.000988, Eq. (5) gives its covariance), so a start frame
that already holds the outer-scale band should need no spin-up at all.

THE METHOD. Three starts of the SAME class, at N = 512, dx = 1 cm,
r0 = 10 cm, L0 = 25 m, `n_columns = 2`, 16 seeds each:

  plain     the stock start, `ft_phase_screen`.
  olb       the olb production generator, `olb.waveoptics.turbulence.screens.
            phase_screen`, three subharmonic levels (Lane, Glindemann and
            Dainty, DOI 10.1088/0959-7174/2/3/003).
  aotools   `aotools.turbulence.phasescreen.ft_sh_phase_screen`, the same
            three-level recipe (it carries the shared-seed quirk of FINDINGS).

Each start reads three metrics against the theory at row 0 (the initial
frame), and again after 512 and 2048 `add_row` calls, so the test also shows
whether the recursion KEEPS a correct start or decays it toward the plain
level (the frame-width effect of FINDINGS Q6 would show as a decay):

  1. the frame variance <phi^2> against B(0) (Assemat and Wilson Eq. (5)).
     REPORT ONLY: see the piston note.
  2. the frame piston variance across the seeds against the double
     integral of B(r) over the frame (`infinite_screen_stats.
     grid_piston_variance`). REPORT ONLY.
  3. the PISTON-REMOVED frame variance against B(0) minus the piston
     theory. This is the variance the optics sees.
  4. the Z-tilt angle variance over a 1.0 m pupil against the Noll filter
     integral (Noll, DOI 10.1364/JOSA.66.000207), pooled over x and y.

THE PISTON NOTE. Both subharmonic generators SUBTRACT the mean of the
low-frequency part (aotools `ft_sh_phase_screen`, and `ScreenFactory` in
olb), so a subharmonic start holds a zero frame piston by construction. A
frame-wide constant phase does nothing to a propagation, so metrics 1 and 2
are reports and the verdict reads metrics 3 and 4 only. The first run of this
script (2026-09-08, 16 seeds) judged metrics 1 and 2 and read FAIL for that
reason alone; its Z-tilt read 1.11 (olb) and 0.89 (aotools) at row 0 against
0.61 for the plain start.

A ratio near 1.0 at row 0 means the spin-up is gone. A ratio that holds from
row 0 to row 2048 means the recursion keeps it.

Outputs, next to this script:
    data/subharmonic_start.csv     start, rows, metric, measured, theory, ratio
    figures/subharmonic_start.png  the three ratios against the row count

Run from the repository root:
    python -m validation.screens.subharmonic_start
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

from aotools.turbulence.infinitephasescreen import (        # noqa: E402
    PhaseScreenVonKarman)
from aotools.turbulence.phasescreen import ft_sh_phase_screen   # noqa: E402

from olb.waveoptics.turbulence.screens import phase_screen  # noqa: E402
from olb.waveoptics.schmidt.turbulence import (             # noqa: E402
    von_karman_phase_psd)
from validation.screens import helpers                     # noqa: E402
from validation.screens.infinite_screen_stats import (     # noqa: E402
    N, DX_M, R0_M, L0_WIDE_M, WAVELENGTH_M, K_RAD_M, PUPIL_D_M,
    grid_piston_variance)

HERE = os.path.dirname(os.path.abspath(__file__))
CSV_PATH = os.path.join(HERE, 'data', 'subharmonic_start.csv')
PNG_PATH = os.path.join(HERE, 'figures', 'subharmonic_start.png')

N_SEEDS = 32
ROW_CHECKPOINTS = (0, 512, 2048)
MASTER_SEED = 20260908
STARTS = ('plain', 'olb', 'aotools')
BAND = 0.20                     # the ratio band of the spin-up verdict


class SubharmonicStart(PhaseScreenVonKarman):
    '''The stock class with a selectable initial frame.

    The class attribute `start` names the generator of the initial frame.
    The A and B matrices, the stencil and `add_row` are the stock ones.
    '''

    start = 'plain'

    def make_initial_screen(self):
        self._R = np.random.default_rng(self.random_seed)
        if self.start == 'plain':
            return super().make_initial_screen()
        if self.start == 'olb':
            child = int(self._R.integers(0, 2 ** 31 - 1))
            scrn = phase_screen(self.r0, self.nx_size, self.pixel_scale,
                                L0_m=self.L0, seed=child, subharmonics=True)
        elif self.start == 'aotools':
            scrn = ft_sh_phase_screen(self.r0, self.nx_size, self.pixel_scale,
                                      self.L0, 1e-10, seed=self._R)
        else:
            raise ValueError(self.start)
        self._scrn = np.asarray(scrn, dtype=float)[:self.stencil_length,
                                                   :self.nx_size]


def make_screen(start, seed):
    '''Build one extruded screen with the named start. Silence the stock print.'''
    klass = type(f'Start_{start}', (SubharmonicStart,), {'start': start})
    with contextlib.redirect_stdout(io.StringIO()):
        return klass(N, DX_M, R0_M, L0_WIDE_M, random_seed=int(seed),
                     n_columns=2)


def main():
    t0 = time.time()
    mask = helpers.pupil_mask(N, DX_M, PUPIL_D_M)
    d_eff = helpers.mask_diameter(mask, DX_M)
    psd = lambda f: von_karman_phase_psd(f, R0_M, L0_WIDE_M)   # noqa: E731
    a2_var = helpers.tilt_filter_variance(psd, d_eff)
    z_theory = a2_var / (K_RAD_M * d_eff / 4.0) ** 2
    b0_theory = float(helpers.vk_covariance_closed(np.array([0.0]), R0_M,
                                                   L0_WIDE_M)[0])
    piston_theory = grid_piston_variance(L0_WIDE_M)
    seeds = helpers.spawn_seeds(MASTER_SEED, N_SEEDS)

    print('ARM 5. A SUBHARMONIC START AGAINST THE PLAIN START')
    print(f'   N = {N}, dx = {DX_M * 1e2:.0f} cm, r0 = {R0_M * 1e2:.0f} cm, '
          f'L0 = {L0_WIDE_M:.0f} m, n_columns = 2, {N_SEEDS} seeds')
    print(f'   theory: B(0) = {b0_theory:.1f} rad^2, frame piston = '
          f'{piston_theory:.1f} rad^2, Z-tilt var = {z_theory:.3e} rad^2 '
          f'over {d_eff:.3f} m')
    print('')

    # Extrude every seed one time. Read the frame at each checkpoint.
    rows_out = []
    results = {}
    for start in STARTS:
        frames = {k: [] for k in ROW_CHECKPOINTS}
        t1 = time.time()
        for seed in seeds:
            screen = make_screen(start, seed)
            done = 0
            for k in ROW_CHECKPOINTS:
                while done < k:
                    screen.add_row()
                    done += 1
                frames[k].append(screen.scrn.copy())
        for k in ROW_CHECKPOINTS:
            fr = np.array(frames[k])
            var = float(np.mean(fr * fr))
            means = fr.mean(axis=(1, 2))
            piston = float(np.mean(means ** 2))
            var_np = float(np.mean((fr - means[:, None, None]) ** 2))
            angles = []
            for f in fr:
                zx, zy, _, _ = helpers.zernike_tilt(f, mask, DX_M,
                                                    WAVELENGTH_M)
                angles.extend([zx, zy])
            angles = np.array(angles)
            ztilt = float(np.mean(angles ** 2))
            for metric, meas, want in (('frame_var', var, b0_theory),
                                       ('piston_var', piston, piston_theory),
                                       ('nopiston_var', var_np,
                                        b0_theory - piston_theory),
                                       ('ztilt_var', ztilt, z_theory)):
                results[(start, k, metric)] = meas / want
                rows_out.append((start, k, metric, meas, want, meas / want))
        print(f'   {start:<8} extruded in {time.time() - t1:6.1f} s')

    print('')
    print(f'   {"start":<9}{"rows":>6}{"var/B0":>10}{"piston":>10}'
          f'{"no-piston":>11}{"Z-tilt":>10}')
    for start in STARTS:
        for k in ROW_CHECKPOINTS:
            print(f'   {start:<9}{k:>6}'
                  f'{results[(start, k, "frame_var")]:>10.3f}'
                  f'{results[(start, k, "piston_var")]:>10.3f}'
                  f'{results[(start, k, "nopiston_var")]:>11.3f}'
                  f'{results[(start, k, "ztilt_var")]:>10.3f}')
    print(f'   relative standard error: the piston about '
          f'{100 * np.sqrt(2.0 / N_SEEDS):.0f} percent, the Z-tilt about '
          f'{100 * np.sqrt(2.0 / (2 * N_SEEDS)):.0f} percent, the '
          f'no-piston variance smaller (many pixels per frame)')
    print('   var/B0 and piston are REPORTS: the subharmonic generators '
          'zero the frame mean.')
    print('')

    # The verdict bands read the optical metrics only.
    judged = ('nopiston_var', 'ztilt_var')
    for start in ('olb', 'aotools'):
        ok0 = all(abs(results[(start, 0, m)] - 1.0) < BAND for m in judged)
        okk = all(abs(results[(start, 2048, m)] - 1.0) < BAND
                  for m in judged)
        print(f'   [{"PASS" if ok0 else "FAIL"}] {start}: row 0 inside '
              f'{BAND:.2f} of theory on the no-piston variance and the '
              f'Z-tilt (no spin-up)')
        print(f'   [{"PASS" if okk else "FAIL"}] {start}: row 2048 inside '
              f'{BAND:.2f} of theory on the same two (kept)')

    os.makedirs(os.path.dirname(CSV_PATH), exist_ok=True)
    with open(CSV_PATH, 'w', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(['start', 'rows', 'metric', 'measured', 'theory',
                         'ratio'])
        writer.writerows(rows_out)

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 4, figsize=(17, 4))
    for ax, metric, title in zip(axes,
                                 ('frame_var', 'piston_var', 'nopiston_var',
                                  'ztilt_var'),
                                 ('frame variance / B(0) (report)',
                                  'frame piston variance / theory (report)',
                                  'piston-removed variance / theory',
                                  'Z-tilt variance over 1 m / Noll')):
        for start in STARTS:
            ax.plot([max(k, 1) for k in ROW_CHECKPOINTS],
                    [results[(start, k, metric)] for k in ROW_CHECKPOINTS],
                    'o-', label=start)
        ax.axhline(1.0, color='k', lw=1)
        ax.axhspan(1 - BAND, 1 + BAND, color='0.85', zorder=0)
        ax.set_xscale('log')
        ax.set_xlabel('add_row calls (row 0 drawn at 1)')
        ax.set_title(title, fontsize=10)
        ax.grid(alpha=0.3)
    axes[0].legend()
    fig.suptitle(f'THE START OF THE EXTRUSION. N {N}, L0 {L0_WIDE_M:.0f} m, '
                 f'r0 {R0_M * 1e2:.0f} cm, n_columns 2, {N_SEEDS} seeds')
    fig.tight_layout()
    os.makedirs(os.path.dirname(PNG_PATH), exist_ok=True)
    fig.savefig(PNG_PATH, dpi=150)
    print(f'   csv: {CSV_PATH}')
    print(f'   figure: {PNG_PATH}')
    print(f'   total {time.time() - t0:.0f} s')


if __name__ == '__main__':
    main()
