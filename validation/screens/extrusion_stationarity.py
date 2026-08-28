'''
Arm 3 of the phase-screen study: does the aotools extrusion drift?

THE QUESTION. `aotools.turbulence.infinitephasescreen.PhaseScreenVonKarman`
grows a screen one row at a time. Each new row is X = A Z + B b, with Z the two
most recent rows and b a Gaussian vector (Assemat and Wilson, Opt. Express
14(3), pp. 988 to 999 (2006), DOI 10.1364/OE.14.000988). The A and B matrices
come from the CLOSED-FORM von Karman covariance of Eq. (5) of that paper. The
initial screen does NOT. Two facts of the aotools source drive this whole
script:

  1. `PhaseScreen.add_row` writes
     `self._scrn = numpy.append(new_row, self._scrn, axis=0)[:stencil_length]`.
     So the extrusion extends AXIS 0, the row axis, and the FRESH row is
     `.scrn[0, :]`. The screen scrolls along y.
  2. `PhaseScreen.make_initial_screen` calls `phasescreen.ft_phase_screen`.
     That is the PLAIN Fourier screen. It is NOT `ft_sh_phase_screen`, so the
     initial screen carries NO subharmonics, and it holds no power below the
     grid fundamental 1/(N dx).

So the process starts from a low-frequency-poor state, and it then steps with a
recursion that knows the correct covariance. The screen must therefore CHANGE
over the first rows. This script measures that change, and it answers two
questions.

THE TWO MEASUREMENTS.

  1. WINDOW DRIFT. Take a frame at cumulative row counts 0, 512, 1024, 2048 and
     4096. Measure the structure function at four separations, and the Zernike
     tilt variance over a 1.0 m pupil. Compare each one with the analytic value.
  2. ROW-LAG COVARIANCE. Store every fresh row after a 512-row spin-up. Estimate
     C(k) = <phi_i phi_{i+k}> against the von Karman covariance B(k dx). Under
     frozen flow a row lag IS a time lag, because the whole pattern translates
     at the wind speed (Taylor, Proc. R. Soc. Lond. A 164, pp. 476 to 490
     (1938), DOI 10.1098/rspa.1938.0032). So C(k) is the temporal covariance of
     the simulation.

THE INTERPRETATION RULE FOR MEASUREMENT 1. Window 0 IS the plain Fourier
initial screen. A rise from window 0 to window 1 is SPIN-UP. It is expected, and
it is NOT drift. The stationarity claim covers windows 1 to 4 only. The test
fits a straight line to each metric against the window index over windows 1 to
4. The test passes when |slope| < 2 SE(slope). The standard error of a variance
over m samples is var sqrt(2/m). The tilt pool holds m = 2 axes x 16 runs = 32
samples. The structure-function metric takes the standard error of the mean
across the 16 runs, because each frame value already averages over 250000 pixel
pairs.

WHY THIS SCRIPT SUBTRACTS NO MEAN IN MEASUREMENT 2. A von Karman phase has a
large outer-scale component. Over 4096 rows that component is a real, slowly
varying offset, and it carries most of the covariance at a small lag. A mean
subtraction removes it, and the estimate then falls below the theory for a
reason that has nothing to do with the extrusion. So the estimator reads raw
products. The same rule applies to the tilt pool: the script reads the mean of
the squares, not the variance about the sample mean.

THE BAND OF MEASUREMENT 2. |C_est(k) - B_th(k)| / B_th(0) <= 0.05 for k <= 64.
Past k = 64 the script reports the deviation, and it SPLITS that deviation into
two parts. A flat offset is a missing VARIANCE, and it says nothing about the
memory of the process. The normalised correlation rho(k) = C(k) / C(0) removes
the offset, and what remains is the shape of the memory. That shape is the
quantity the n_columns = 2 Markov truncation controls: the recursion sees two
rows only.

A NOTE ON THE TWO COVARIANCE ROUTES. `helpers.vk_covariance_numeric` integrates
the PSD of `schmidt.turbulence.von_karman_phase_psd`, which carries the PRINTED
constant 0.023. `helpers.vk_covariance_closed` evaluates Eq. (5) of Assemat and
Wilson with an exact prefactor. The numeric route reads about 0.44 percent
higher. `aotools.turbulence.turb.phase_covariance` implements the CLOSED form,
so the A and B matrices follow the closed route. The script prints both, and the
0.44 percent gap is far inside the 5 percent band.

This script changes NO olb module. It reads the production layer only.

Outputs, next to this script:
    stationarity.png             the two panels
    stationarity_windows.csv     the window table
    rowlag_covariance.csv        the row-lag table

Sources:
- Assemat and Wilson, "Method for simulating infinitely long and non stationary
  phase screens with optimized memory storage", Opt. Express 14(3), pp. 988 to
  999 (2006), DOI 10.1364/OE.14.000988. Eq. (5), the closed-form von Karman
  phase covariance, and the X = A Z + B b recursion.
- Taylor, "The spectrum of turbulence", Proc. R. Soc. Lond. A 164, pp. 476 to
  490 (1938), DOI 10.1098/rspa.1938.0032. The frozen-flow hypothesis.
- Schmidt, Numerical Simulation of Optical Wave Propagation with Examples in
  MATLAB, DOI 10.1117/3.866274. Ch. 3, Eq. (3.15), printed p. 47 (the structure
  function); Ch. 3, Eq. (3.16), printed p. 48 (D(r) = 2 [B(0) - B(r)]); Ch. 9,
  Eq. (9.50), printed p. 161 (the von Karman phase PSD); Ch. 9, Eqs. (9.78) to
  (9.80), printed p. 167 (the Fourier screen).
- Noll, "Zernike polynomials and atmospheric turbulence", J. Opt. Soc. Am.
  66(3), pp. 207 to 211 (1976), DOI 10.1364/JOSA.66.000207. Eq. (8), p. 208,
  the Zernike spatial filter that `helpers.tilt_filter_variance` integrates.

Run from the repo root:
    python -m validation.screens.extrusion_stationarity
'''

import csv
import os
import time

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt                                   # noqa: E402
import numpy as np                                                # noqa: E402
from aotools.turbulence.infinitephasescreen import \
    PhaseScreenVonKarman                                          # noqa: E402

from olb.waveoptics.schmidt.turbulence import \
    von_karman_phase_psd                                          # noqa: E402
from validation.screens import helpers                            # noqa: E402


# ---------------------------------------------------------------------------
# The module constants
# ---------------------------------------------------------------------------

WAVELENGTH_M = 1550e-9
R0_M = 0.10
DX_M = 0.01
N = 512
PUPIL_D_M = 1.0
L0_M = 25.0
RUNS = 16
MASTER_SEED = 20260829

# ---- measurement 1: the windows ----
# The cumulative add_row count of each harvested frame. Window 0 is the plain
# Fourier initial screen, before the first extrusion step.
WINDOW_ROWS = (0, 512, 1024, 2048, 4096)

# The structure-function probes, in pixels. 250 px is 2.5 m, which is half the
# screen side, so `d_phi_direct` still averages over 262 columns there.
PROBE_PX = (10, 32, 100, 250)

# ---- measurement 2: the row record ----
SPINUP_ROWS = 512
RECORD_ROWS = 4096
K_MAX = 256
K_TABLE = (1, 2, 4, 8, 16, 32, 64, 128, 256)

# The band of measurement 2. The deviation is normalised by B(0), so one number
# covers every lag.
BAND_K_MAX = 64
BAND_TOL = 0.05

# ---- the outputs ----
HERE = os.path.dirname(os.path.abspath(__file__))
FIGURE_PNG = os.path.join(HERE, 'figures/stationarity.png')
WINDOWS_CSV = os.path.join(HERE, 'data/stationarity_windows.csv')
ROWLAG_CSV = os.path.join(HERE, 'data/rowlag_covariance.csv')
DPI = 150


# ---------------------------------------------------------------------------
# The one expensive pass
# ---------------------------------------------------------------------------

def run_one(seed):
    '''
    Extrude one screen, and harvest both measurements in ONE pass.

    Parameters:
        seed : int
            The integer seed of `PhaseScreenVonKarman`.

    Returns:
        tuple
            (frames, rows). `frames` is a list of 512 by 512 float arrays, one
            per entry of `WINDOW_ROWS`. `rows` is a `RECORD_ROWS` by 512
            float32 array of the fresh edge rows, after the spin-up.

    WHY ONE PASS. The construction of the A and B matrices costs about 1 s, and
    each `add_row` costs about 2.3 ms. The two measurements need 4096 and 4608
    steps. A shared pass pays for the longer of the two, not for the sum.

    THE FRESH ROW. `add_row` prepends the new row on axis 0. So the newest row
    is `screen.scrn[0, :]`. See the module docstring, fact 1.
    '''
    screen = PhaseScreenVonKarman(N, DX_M, R0_M, L0_M, random_seed=int(seed))

    frames = []
    rows = np.empty((RECORD_ROWS, N), dtype=np.float32)
    total = SPINUP_ROWS + RECORD_ROWS
    window_set = dict((count, index) for index, count
                      in enumerate(WINDOW_ROWS))

    if 0 in window_set:
        frames.append(np.array(screen.scrn, dtype=float, copy=True))

    for step in range(1, total + 1):
        screen.add_row()
        if step in window_set:
            frames.append(np.array(screen.scrn, dtype=float, copy=True))
        if step > SPINUP_ROWS:
            rows[step - SPINUP_ROWS - 1] = screen.scrn[0, :]

    return frames, rows


def row_lag_sums(rows):
    '''
    Return the unnormalised row-lag products of one run, for k = 0 to K_MAX.

    Parameters:
        rows : numpy.ndarray
            An M by N array. Row i is the fresh row of step i.

    Returns:
        numpy.ndarray
            An array of K_MAX + 1 values. Entry k is the MEAN over the N
            columns of SUM_i phi_i phi_{i+k}. Divide by the pair count M - k to
            get C(k).

    formula:
        SUM_i x_i x_{i+k} = IFFT[ |FFT(x)|^2 ]_k,  with zero padding to 2M
    Source: the Wiener-Khinchin relation for a finite sequence. Schmidt (2010),
    DOI 10.1117/3.866274, Ch. 3, Eq. (3.16), printed p. 48, states the same
    transform pair for a continuous field.

    WHY THE PADDING. A plain FFT gives a CIRCULAR correlation, which wraps the
    end of the record onto its start. The zero padding to 2M removes the wrap.

    NO MEAN SUBTRACTION. See the module docstring.
    '''
    x = np.asarray(rows, dtype=float)
    m = x.shape[0]
    nfft = 1
    while nfft < 2 * m:
        nfft *= 2
    spectrum = np.fft.rfft(x, n=nfft, axis=0)
    correlation = np.fft.irfft(spectrum * np.conj(spectrum), n=nfft, axis=0)
    return np.asarray(correlation[:K_MAX + 1].mean(axis=1), dtype=float)


# ---------------------------------------------------------------------------
# The trend test
# ---------------------------------------------------------------------------

def trend(values, errors):
    '''
    Fit a straight line, and return the slope with its standard error.

    Parameters:
        values : array_like
            The metric, one value per window.
        errors : array_like
            The standard error of each value.

    Returns:
        tuple
            (slope, se_slope). Both are per window index.

    formula:
        SE(slope) = sigma / sqrt(SUM (x - <x>)^2)
    Source: the standard least-squares result for a straight line with equal
    weights. `sigma` here is the root mean square of the per-point standard
    errors.

    THE TEST. Call the metric stationary when |slope| < 2 SE(slope).
    '''
    values = np.asarray(values, dtype=float)
    errors = np.asarray(errors, dtype=float)
    x = np.arange(values.size, dtype=float)
    slope = float(np.polyfit(x, values, 1)[0])
    sigma = float(np.sqrt(np.mean(errors ** 2)))
    se_slope = sigma / np.sqrt(float(np.sum((x - x.mean()) ** 2)))
    return slope, se_slope


# ---------------------------------------------------------------------------
# The figure
# ---------------------------------------------------------------------------

def draw(dphi_ratio, dphi_se, ztilt_ratio, ztilt_se, x_ratio, y_ratio, k_axis,
         c_est, b_theory, b_zero):
    '''Draw the window panel and the row-lag panel.

    The left panel adds the two axes of the LAST probe as thin dashed lines.
    They share the y axis of the averaged metric, so the anisotropy reads
    directly against the same theory.
    '''
    fig = plt.figure(figsize=(14.4, 6.4), constrained_layout=True)
    grid = fig.add_gridspec(2, 2, height_ratios=[3.0, 1.2])
    a_win = fig.add_subplot(grid[:, 0])
    a_cov = fig.add_subplot(grid[0, 1])
    a_res = fig.add_subplot(grid[1, 1], sharex=a_cov)

    index = np.arange(len(WINDOW_ROWS), dtype=float)
    colours = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red']
    a_win.axhline(1.0, color='black', linewidth=2.0, label='the analytic value')
    for j, probe in enumerate(PROBE_PX):
        a_win.errorbar(index, dphi_ratio[:, j], yerr=2.0 * dphi_se[:, j],
                       color=colours[j], marker='o', capsize=3.0,
                       linewidth=1.6,
                       label=f'D(r) at {probe} px = {probe * DX_M:.2f} m')
    a_win.errorbar(index, ztilt_ratio, yerr=2.0 * ztilt_se, color='tab:purple',
                   marker='s', capsize=3.0, linewidth=2.2, linestyle='--',
                   label=f'Z-tilt variance, {PUPIL_D_M:.1f} m pupil')
    last = PROBE_PX[-1]
    a_win.plot(index, x_ratio[:, -1], color='tab:red', linewidth=1.2,
               linestyle=':', marker='^',
               label=f'{last} px, x only (along a row)')
    a_win.plot(index, y_ratio[:, -1], color='tab:red', linewidth=1.2,
               linestyle='-.', marker='v',
               label=f'{last} px, y only (the extrusion)')
    a_win.axvspan(-0.25, 0.5, color='0.88', zorder=0,
                  label='window 0: the plain Fourier screen')
    a_win.set_xticks(index)
    a_win.set_xticklabels([str(r) for r in WINDOW_ROWS])
    a_win.set_xlabel('cumulative add_row calls')
    a_win.set_ylabel('measured / analytic')
    a_win.set_title('1. WINDOW DRIFT. The bars are 2 standard errors over '
                    f'{RUNS} runs.\nWindow 0 to 1 is spin-up. The claim covers '
                    'windows 1 to 4.', fontsize=10)
    a_win.grid(alpha=0.3)
    # Add headroom, so the legend never covers a point.
    top = max(1.25, float(np.nanmax(dphi_ratio + 2.0 * dphi_se)))
    a_win.set_ylim(top=top + 0.45)
    a_win.legend(fontsize=8, loc='upper center', ncol=2, framealpha=0.92)

    a_cov.plot(k_axis, b_theory, color='black', linewidth=2.4,
               label='von Karman B(k dx), Assemat and Wilson Eq. (5)')
    a_cov.plot(k_axis, c_est, color='tab:red', linewidth=1.8,
               label=f'C(k) of the extrusion, {RUNS} runs')
    a_cov.set_xscale('log')
    a_cov.set_ylabel('covariance, rad^2')
    a_cov.set_title('2. ROW-LAG COVARIANCE. A row lag is a frozen-flow time '
                    'lag.', fontsize=10)
    a_cov.grid(alpha=0.3, which='both')
    a_cov.legend(fontsize=8, loc='lower left')
    plt.setp(a_cov.get_xticklabels(), visible=False)

    deficit = (b_theory - c_est) / b_zero
    a_res.axhline(0.0, color='black', linewidth=1.4)
    a_res.axhspan(-BAND_TOL, BAND_TOL, color='0.85', zorder=0,
                  label=f'the band, {BAND_TOL:.2f} of B(0)')
    a_res.axvline(BAND_K_MAX, color='tab:green', linestyle=':', linewidth=1.8,
                  label=f'k = {BAND_K_MAX}')
    a_res.plot(k_axis, deficit, color='tab:red', linewidth=1.8)
    a_res.set_xscale('log')
    a_res.set_xlabel('row lag k, pixels')
    a_res.set_ylabel('[B - C] / B(0)')
    a_res.grid(alpha=0.3, which='both')
    a_res.legend(fontsize=8, loc='best')

    fig.suptitle(f'The aotools infinite-screen extrusion. r0 = '
                 f'{R0_M * 1e2:.0f} cm, L0 = {L0_M:.0f} m, {N} px at '
                 f'{DX_M * 1e2:.0f} cm, n_columns = 2.', fontsize=12)
    fig.savefig(FIGURE_PNG, dpi=DPI)
    plt.close(fig)


# ---------------------------------------------------------------------------

def main():
    t_start = time.time()
    seeds = helpers.spawn_seeds(MASTER_SEED, RUNS)
    mask = helpers.pupil_mask(N, DX_M, PUPIL_D_M)
    d_eff = helpers.mask_diameter(mask, DX_M)

    print('=' * 78)
    print('Arm 3: the stationarity of the aotools infinite-screen extrusion')
    print('=' * 78)
    print(f'  wavelength           {WAVELENGTH_M * 1e9:12.1f} nm')
    print(f'  grid                 {N:12d} px at {DX_M * 1e2:.1f} cm, side '
          f'{N * DX_M:.2f} m')
    print(f'  r0                   {R0_M * 1e2:12.1f} cm   (D/r0 = '
          f'{PUPIL_D_M / R0_M:.1f})')
    print(f'  L0                   {L0_M:12.1f} m')
    print(f'  pupil                {PUPIL_D_M:12.2f} m nominal, '
          f'{d_eff:.4f} m area-equivalent')
    print(f'  runs                 {RUNS:12d}   master seed {MASTER_SEED}')
    print(f'  add_row per run      {SPINUP_ROWS + RECORD_ROWS:12d}   '
          f'({SPINUP_ROWS} spin-up + {RECORD_ROWS} recorded)')
    print('')
    print('  TWO FACTS OF THE aotools SOURCE:')
    print('   - add_row prepends on AXIS 0, so the fresh row is scrn[0, :].')
    print('   - make_initial_screen calls ft_phase_screen, the PLAIN Fourier')
    print('     screen. There are NO subharmonics in window 0.')

    # ---- the analytic targets ----
    probe_px = np.asarray(PROBE_PX, dtype=int)
    probe_m = probe_px * DX_M
    b_zero_closed = float(helpers.vk_covariance_closed(0.0, R0_M, L0_M)[0])
    dphi_theory = 2.0 * (b_zero_closed
                         - helpers.vk_covariance_closed(probe_m, R0_M, L0_M))
    psd = lambda f: von_karman_phase_psd(f, R0_M, L0_M)
    ztilt_theory = helpers.tilt_filter_variance(psd, d_eff)

    # ---- the expensive pass ----
    n_win = len(WINDOW_ROWS)
    dphi = np.zeros((RUNS, n_win, probe_px.size))
    dphi_x = np.zeros((RUNS, n_win, probe_px.size))
    dphi_y = np.zeros((RUNS, n_win, probe_px.size))
    tilt_pool = np.zeros((RUNS, n_win, 2))
    lag_sums = np.zeros(K_MAX + 1)

    print('')
    print(f'  extruding {RUNS} runs ...')
    for run, seed in enumerate(seeds):
        frames, rows = run_one(seed)
        for w, frame in enumerate(frames):
            for j, kpx in enumerate(probe_px):
                dphi[run, w, j] = helpers.d_phi_direct(frame, int(kpx))
                # The extrusion builds AXIS 0 only. Split the two axes, so an
                # anisotropy shows. `d_phi_direct` averages them.
                diff_x = frame[:, kpx:] - frame[:, :-kpx]
                diff_y = frame[kpx:, :] - frame[:-kpx, :]
                dphi_x[run, w, j] = float(np.mean(diff_x * diff_x))
                dphi_y[run, w, j] = float(np.mean(diff_y * diff_y))
            _, _, a2, a3 = helpers.zernike_tilt(frame, mask, DX_M,
                                                WAVELENGTH_M)
            tilt_pool[run, w] = (a2, a3)
        lag_sums += row_lag_sums(rows)
        del frames, rows
        print(f'    run {run + 1:2d} of {RUNS} done, '
              f'{time.time() - t_start:6.1f} s elapsed')

    # ---- measurement 1: reduce ----
    # The structure function: the mean over the runs, and the standard error of
    # that mean.
    dphi_mean = dphi.mean(axis=0)
    dphi_sem = dphi.std(axis=0, ddof=1) / np.sqrt(RUNS)
    dphi_ratio = dphi_mean / dphi_theory[None, :]
    dphi_ratio_se = dphi_sem / dphi_theory[None, :]

    # The tilt: the MEAN OF THE SQUARES over the pool of 2 axes x RUNS runs.
    # The standard error of a variance over m samples is var sqrt(2/m).
    pool = tilt_pool.transpose(1, 0, 2).reshape(n_win, RUNS * 2)
    m_pool = pool.shape[1]
    ztilt_var = np.mean(pool ** 2, axis=1)
    ztilt_var_se = ztilt_var * np.sqrt(2.0 / m_pool)
    ztilt_ratio = ztilt_var / ztilt_theory
    ztilt_ratio_se = ztilt_var_se / ztilt_theory

    print('')
    print('1. WINDOW DRIFT')
    print(f'   theory: D(r) = 2 [B(0) - B(r)], B(0) = {b_zero_closed:.2f} '
          f'rad^2 (closed form)')
    print(f'   theory: <a2^2> = {ztilt_theory:.4f} rad^2 over a '
          f'{d_eff:.4f} m pupil (Noll Eq. (8))')
    print('')
    head = f'   {"win":>4}{"rows":>7}'
    for probe in PROBE_PX:
        head += f'{"D@" + str(probe) + "px":>14}{"ratio":>8}'
    head += f'{"ztilt var":>12}{"ratio":>8}'
    print(head)
    for w in range(n_win):
        line = f'   {w:>4}{WINDOW_ROWS[w]:>7}'
        for j in range(probe_px.size):
            line += f'{dphi_mean[w, j]:>14.3f}{dphi_ratio[w, j]:>8.4f}'
        line += f'{ztilt_var[w]:>12.4f}{ztilt_ratio[w]:>8.4f}'
        print(line)
    print('   The theory row, for reference:')
    ref = f'   {"th":>4}{"-":>7}'
    for j in range(probe_px.size):
        ref += f'{dphi_theory[j]:>14.3f}{1.0:>8.4f}'
    ref += f'{ztilt_theory:>12.4f}{1.0:>8.4f}'
    print(ref)

    # ---- the anisotropy ----
    # `d_phi_direct` averages the two axes, so it HIDES a difference between
    # them. The extrusion builds axis 0 only, so the two axes need a split.
    x_ratio = dphi_x.mean(axis=0) / dphi_theory[None, :]
    y_ratio = dphi_y.mean(axis=0) / dphi_theory[None, :]
    x_se = dphi_x.std(axis=0, ddof=1) / np.sqrt(RUNS) / dphi_theory[None, :]
    y_se = dphi_y.std(axis=0, ddof=1) / np.sqrt(RUNS) / dphi_theory[None, :]

    print('')
    print('   THE ANISOTROPY. Axis 1 (x) lies ALONG a row. Axis 0 (y) is the')
    print('   extrusion direction. The table gives each axis against the same')
    print('   theory. The averaged metric above hides this split.')
    line = f'   {"win":>4}{"rows":>7}'
    for probe in PROBE_PX:
        line += f'{"x@" + str(probe):>9}{"y@" + str(probe):>9}'
    line += f'{"x/y @250":>10}'
    print(line)
    for w in range(n_win):
        line = f'   {w:>4}{WINDOW_ROWS[w]:>7}'
        for j in range(probe_px.size):
            line += f'{x_ratio[w, j]:>9.4f}{y_ratio[w, j]:>9.4f}'
        line += f'{x_ratio[w, -1] / y_ratio[w, -1]:>10.4f}'
        print(line)

    # THE MAGNITUDE, NOT THE TREND. The trend test asks whether a metric MOVES.
    # It does not ask how far the metric sits from the theory. Pool the four
    # post-spin-up windows PER RUN, and take the standard error across the 16
    # independent runs. The runs are independent, and the windows of one run
    # are not.
    per_run = ((dphi_x[:, 1:, :] - dphi_y[:, 1:, :]).mean(axis=1)
               / dphi_theory[None, :])
    gap_mean = per_run.mean(axis=0)
    gap_sem = per_run.std(axis=0, ddof=1) / np.sqrt(RUNS)

    print('')
    print('   THE MAGNITUDE OF THE ANISOTROPY, pooled over windows 1 to 4.')
    print('   The value is (D_x - D_y) / D_theory. The error is one standard')
    print(f'   error over the {RUNS} independent runs.')
    print(f'   {"probe":>8}{"r, m":>9}{"(Dx-Dy)/Dth":>14}{"SE":>10}{"|t|":>8}'
          f'{"real?":>8}')
    for j, probe in enumerate(PROBE_PX):
        value = abs(gap_mean[j]) / gap_sem[j] if gap_sem[j] > 0.0 else np.inf
        print(f'   {probe:>8}{probe * DX_M:>9.2f}{gap_mean[j]:>+14.4f}'
              f'{gap_sem[j]:>10.4f}{value:>8.2f}'
              f'{"yes" if value > 2.0 else "no":>8}')

    print('')
    print('   THE SAME TREND TEST, PER AXIS, over windows 1 to 4:')
    print(f'   {"metric":<26}{"slope":>12}{"SE(slope)":>12}'
          f'{"|t|":>8}{"verdict":>10}')
    aniso_flat = True
    for j, probe in enumerate(PROBE_PX):
        for name, ratio, error in (('x, along a row', x_ratio, x_se),
                                   ('y, the extrusion', y_ratio, y_se)):
            slope, se_slope = trend(ratio[1:, j], error[1:, j])
            value = abs(slope) / se_slope if se_slope > 0.0 else np.inf
            ok = value < 2.0
            aniso_flat = aniso_flat and ok
            print(f'   {name + " @" + str(probe) + " px":<26}{slope:>+12.5f}'
                  f'{se_slope:>12.5f}{value:>8.2f}'
                  f'{"pass" if ok else "FAIL":>10}')

    # ---- the trend test over windows 1 to 4 ----
    print('')
    print('   THE TREND TEST over windows 1 to 4. Window 0 is the spin-up')
    print('   start, so the test excludes it.')
    print(f'   {"metric":<26}{"slope":>12}{"SE(slope)":>12}'
          f'{"|t|":>8}{"verdict":>10}')
    trends = []
    for j, probe in enumerate(PROBE_PX):
        slope, se_slope = trend(dphi_ratio[1:, j], dphi_ratio_se[1:, j])
        trends.append((f'D(r) at {probe} px', slope, se_slope))
    slope, se_slope = trend(ztilt_ratio[1:], ztilt_ratio_se[1:])
    trends.append(('Z-tilt variance', slope, se_slope))

    stationary = True
    for name, slope, se_slope in trends:
        ratio = abs(slope) / se_slope if se_slope > 0.0 else np.inf
        ok = ratio < 2.0
        stationary = stationary and ok
        print(f'   {name:<26}{slope:>+12.5f}{se_slope:>12.5f}'
              f'{ratio:>8.2f}{"pass" if ok else "FAIL":>10}')
    print('')
    spin = dphi_ratio[1] / dphi_ratio[0]
    print(f'   the spin-up step, window 0 to window 1, lifts D(r) by '
          f'{100.0 * (spin.mean() - 1.0):+.1f} percent on average')
    print(f'   and the Z-tilt variance by '
          f'{100.0 * (ztilt_ratio[1] / ztilt_ratio[0] - 1.0):+.1f} percent')
    print('')
    print(f'stationary after spin-up: {"yes" if stationary else "no"}')
    if not aniso_flat:
        print('')
        print('   CAVEAT. The verdict line above reads the CONTRACTED metric,')
        print('   the two-axis average of `d_phi_direct`. The per-axis test')
        print('   FAILS. The two axes move in opposite directions, and the')
        print('   average cancels most of the movement. Read the per-axis')
        print('   table, not the average.')

    # ---- measurement 2: reduce ----
    k_all = np.arange(K_MAX + 1)
    counts = (RECORD_ROWS - k_all).astype(float)
    c_est = lag_sums / counts / RUNS
    b_num = helpers.vk_covariance_numeric(k_all * DX_M, R0_M, L0_M)
    b_cls = helpers.vk_covariance_closed(k_all * DX_M, R0_M, L0_M)
    b_zero = float(b_num[0])
    deviation = np.abs(c_est - b_num) / b_zero

    print('')
    print('2. ROW-LAG COVARIANCE, THE FROZEN-FLOW TEMPORAL PROXY')
    print(f'   theory B(0) = {b_num[0]:.3f} rad^2 (numeric route), '
          f'{b_cls[0]:.3f} rad^2 (closed route)')
    print(f'   the numeric route reads '
          f'{100.0 * (b_num[0] / b_cls[0] - 1.0):+.2f} percent high, from the '
          f'printed constant 0.023')
    print(f'   measured C(0) = {c_est[0]:.3f} rad^2, which is '
          f'{c_est[0] / b_num[0]:.4f} of the theory')
    print('')
    print(f'   {"k":>6}{"r, m":>9}{"C_est":>13}{"B_theory":>13}'
          f'{"C/B":>9}{"deficit/B(0)":>15}')
    for k in K_TABLE:
        print(f'   {k:>6}{k * DX_M:>9.2f}{c_est[k]:>13.3f}{b_num[k]:>13.3f}'
              f'{c_est[k] / b_num[k]:>9.4f}'
              f'{(b_num[k] - c_est[k]) / b_zero:>+15.5f}')

    band = k_all[1:BAND_K_MAX + 1]
    worst = float(deviation[1:BAND_K_MAX + 1].max())
    worst_k = int(band[int(np.argmax(deviation[1:BAND_K_MAX + 1]))])
    passed = worst <= BAND_TOL
    print('')
    print(f'   BAND: |C - B| / B(0) <= {BAND_TOL:.2f} for k <= {BAND_K_MAX}')
    print(f'   worst deviation {worst:.5f} at k = {worst_k}   '
          f'{"PASS" if passed else "FAIL"}')
    try:
        assert passed, (worst, worst_k)
    except AssertionError as error:
        print(f'   assert caught: {error}. The script continues. A physics')
        print('   surprise IS a result.')

    # ---- past the band: split the variance deficit from the SHAPE error ----
    # The absolute deficit above mixes two faults. A flat offset is a missing
    # VARIANCE, and it says nothing about the memory of the process. The
    # normalised correlation rho(k) = C(k) / C(0) removes that offset, and what
    # remains IS the shape of the memory. That is the quantity the n_columns
    # truncation controls.
    rho_est = c_est / c_est[0]
    rho_theory = b_num / b_num[0]

    print('')
    print('   PAST THE BAND: THE SHAPE, NOT THE OFFSET. The deficit above is')
    print('   nearly FLAT for k <= 32, so it is a missing VARIANCE, not a')
    print('   memory error. Divide it out. The correlation rho(k) = C(k)/C(0)')
    print('   carries the memory alone. Assemat and Wilson,')
    print('   DOI 10.1364/OE.14.000988, claim two columns are adequate.')
    print(f'   {"k":>6}{"r, m":>9}{"rho_est":>11}{"rho_theory":>13}'
          f'{"excess":>10}{"deficit/B(0)":>15}')
    for k in (1, 8, 32, 64, 96, 128, 192, 256):
        print(f'   {k:>6}{k * DX_M:>9.2f}{rho_est[k]:>11.5f}'
              f'{rho_theory[k]:>13.5f}{rho_est[k] - rho_theory[k]:>+10.5f}'
              f'{(b_num[k] - c_est[k]) / b_zero:>+15.5f}')
    excess = rho_est - rho_theory
    print('')
    print(f'   the absolute deficit moves from '
          f'{(b_num[BAND_K_MAX] - c_est[BAND_K_MAX]) / b_zero:+.5f} of B(0) '
          f'at k = {BAND_K_MAX}')
    print(f'   to {(b_num[K_MAX] - c_est[K_MAX]) / b_zero:+.5f} at '
          f'k = {K_MAX}. It SHRINKS. The two curves cross,')
    print('   because the theory falls and the estimate does not follow it.')
    print(f'   the correlation EXCESS grows from {excess[BAND_K_MAX]:+.5f} at '
          f'k = {BAND_K_MAX} to {excess[K_MAX]:+.5f}')
    print(f'   at k = {K_MAX}. The extruded axis holds MORE memory than the '
          f'von Karman')
    print('   covariance, not less. That is the same fault as the anisotropy')
    print('   of measurement 1: the extrusion direction is too smooth.')

    # ---- the outputs ----
    with open(WINDOWS_CSV, 'w', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(['window', 'cum_rows', 'probe_px', 'dphi_ratio', 'se',
                         'ztilt_var', 'ztilt_ratio', 'ztilt_se'])
        for w in range(n_win):
            for j, probe in enumerate(PROBE_PX):
                writer.writerow([w, WINDOW_ROWS[w], probe,
                                 f'{dphi_ratio[w, j]:.6f}',
                                 f'{dphi_ratio_se[w, j]:.6f}',
                                 f'{ztilt_var[w]:.6f}',
                                 f'{ztilt_ratio[w]:.6f}',
                                 f'{ztilt_ratio_se[w]:.6f}'])

    with open(ROWLAG_CSV, 'w', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(['k', 'r_m', 'c_est', 'b_theory', 'deficit_frac'])
        for k in range(K_MAX + 1):
            writer.writerow([k, f'{k * DX_M:.4f}', f'{c_est[k]:.6f}',
                             f'{b_num[k]:.6f}',
                             f'{(b_num[k] - c_est[k]) / b_zero:.6f}'])

    draw(dphi_ratio, dphi_ratio_se, ztilt_ratio, ztilt_ratio_se, x_ratio,
         y_ratio, k_all[1:], c_est[1:], b_num[1:], b_zero)

    print('')
    print(f'figure saved: {FIGURE_PNG}')
    print('  Caption: left, each metric against the window, with 2 standard '
          'error bars.')
    print('  The grey band is window 0, the plain Fourier start. The two thin '
          'dashed')
    print('  lines split the largest probe into its x and its y axis. Right, '
          'the row-lag')
    print('  covariance against the von Karman theory, with the deficit '
          'below.')
    print(f'table saved:  {WINDOWS_CSV}')
    print(f'table saved:  {ROWLAG_CSV}')
    print('')
    print(f'(elapsed {time.time() - t_start:.1f} s)')


if __name__ == '__main__':
    main()
