'''
Phases 1 and 2 of PLAN_n_columns.md: what does `n_columns` buy?

THE QUESTION. `aotools.turbulence.infinitephasescreen.PhaseScreenVonKarman`
extrudes a screen one row at a time with X = A Z + B b (Assemat and Wilson,
Opt. Express 14(3), pp. 988 to 999 (2006), DOI 10.1364/OE.14.000988). The
stencil Z holds `n_columns` rows only, so the recursion is a Markov process of
that order. The paper claims two rows are adequate. `FINDINGS.md` measured an
over-correlation along the extrusion axis and blamed that truncation, but it
measured it on aotools 1.0.7, whose covariance kernel ran in float32. The
installed build casts to float64 (upstream PR 111). So the question is open
again: on the FIXED kernel, does a larger `n_columns` buy anything?

PHASE 1, THE SWEEP. `n_columns` in {2, 4, 8, 16, 32}, L0 in {2.56 m, 25 m},
r0 = 0.10 m, dx = 0.01 m. N = 128 for every cell, plus N = 512 at `n_columns`
in {2, 8} as the production-size check. 8 seeds per cell. The spin-up is
2 L0 / dx rows. The metrics are:

  1. rho(k) along the EXTRUSION axis against the von Karman theory.
  2. The TRANSVERSE covariance at the same lags. That axis is the control,
     because the recursion does not touch it.
  3. The structure-function anisotropy D_ext(r) / D_trans(r).
  4. The row PISTON series m(t), the mean of each new row: its chunk-mean
     variance and its power spectrum against the transverse-averaged theory.
  5. The frame variance B(0) and the Z-tilt over a 1.0 m pupil.
  6. The cost: the setup time, the size of cov_zz, and the add_row time.
  7. The guards: no lstsq fallback, the minimum eigenvalue of BB^T, and the
     condition number of cov_zz.

PHASE 2, THE RECORD LENGTH. Two 2026-09-08 probes disagree: a 6144-row record
shows a flat rho floor of 0.13 to 0.16, and a 16384-row record reads the
theory. So the apparent over-correlation can be the WANDER of the row piston
across a short window, not the memory of the recursion. This phase extrudes one
record of 100 L0 / dx rows at `n_columns` 2 and 16, and it reads rho(1 L0) in
windows of 2, 5, 10, 20 and 100 L0. If the apparent rho falls to the theory as
the window grows, the record length is the cause and `n_columns` is not.

THE LAZY CHOICES, where the plan leaves the detail open:
  - The A and the B matrix do NOT depend on the seed. So each cell builds ONE
    screen object and it re-seeds the initial frame for each seed, through
    `make_initial_screen`. That pays the setup cost one time per cell.
  - The guards (the lstsq fallback, the BB^T eigenvalue, the cov_zz condition
    number) also do not depend on the seed. The script records them one time
    per cell. It SKIPS the condition number when cov_zz is larger than
    4096 x 4096, because the eigenvalue solve is too slow there.
  - The record at N = 512 with L0 = 25 m is CUT to 10000 rows (4 L0), not the
    50000 rows (20 L0) of the N = 128 cells. One add_row at N = 512 costs
    1.3 ms against 0.02 ms at N = 128, so the full record would cost more than
    the whole rest of the study. The plan asks for 20 L0 at N = 128 only.
  - The TRANSVERSE axis cannot reach a lag larger than half the grid side
    (0.64 m at N = 128, 2.56 m at N = 512). So the anisotropy band of the plan,
    r = 0.1 to 2 L0, is not reachable at L0 = 25 m. The script tests every
    reachable separation, and it flags which ones sit inside the plan window.
  - Phase 2 uses 4 seeds, the floor of the task, because its record is 5 times
    longer than a Phase 1 record.
  - The transverse covariance reads SQUARE N by N blocks of the row record
    through `infinite_screen_stats.axis_covariance`, so the estimator of the
    sibling script is reused with no change.
  - `extrusion_stationarity.row_lag_sums` caps its output at the module
    constant `K_MAX`. This script RAISES that constant, so it reuses the same
    Wiener-Khinchin estimator at the lags this study needs.

NO MEAN SUBTRACTION, anywhere. A von Karman phase carries real outer-scale
power in the sample mean. A mean subtraction removes it, and every estimate
then falls below the theory for a reason that has nothing to do with the
extrusion. This is the rule of the two sibling scripts.

A FAILED BAND IS A RESULT. No band raises. Every band prints PASS or FAIL.

Outputs, under `data/` and `figures/` next to this script:
    data/ncol_rowlag.csv      rho(k) per axis, per cell
    data/ncol_piston.csv      the chunk-mean variance and the piston spectrum
    data/ncol_cost.csv        the cost, the guards, the variance and the tilt
    data/ncol_window.csv      Phase 2, rho(1 L0) against the window length
    figures/ncol_rowlag.png   rho(k) per n_columns, the theory in black
    figures/ncol_anisotropy.png
    figures/ncol_piston.png
    figures/ncol_cost.png
    figures/ncol_window.png

This script changes NO olb module. It reads the production layer only.

Sources:
- Assemat and Wilson, "Method for simulating infinitely long and non stationary
  phase screens with optimized memory storage", Opt. Express 14(3), pp. 988 to
  999 (2006), DOI 10.1364/OE.14.000988. Eq. (5), the closed-form von Karman
  phase covariance, and the X = A Z + B b recursion.
- Noll, "Zernike polynomials and atmospheric turbulence", J. Opt. Soc. Am.
  66(3), pp. 207 to 211 (1976), DOI 10.1364/JOSA.66.000207. Eq. (8), p. 208,
  the Zernike tilt filter.
- Schmidt, Numerical Simulation of Optical Wave Propagation with Examples in
  MATLAB, DOI 10.1117/3.866274. Ch. 3, Eq. (3.14), printed p. 47 (the
  covariance); Ch. 3, Eq. (3.16), printed p. 48 (D(r) = 2 [B(0) - B(r)]);
  Ch. 9, Eq. (9.50), printed p. 161 (the von Karman phase PSD).
- Taylor, "The spectrum of turbulence", Proc. R. Soc. Lond. A 164, pp. 476 to
  490 (1938), DOI 10.1098/rspa.1938.0032. The frozen-flow hypothesis, which
  makes a row lag a time lag.

Run from the repo root:
    python -m validation.screens.n_columns_sweep
'''

import contextlib
import csv
import io
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
from validation.screens import extrusion_stationarity as es       # noqa: E402
from validation.screens import helpers                            # noqa: E402
from validation.screens.infinite_screen_stats import \
    axis_covariance                                               # noqa: E402


# ---------------------------------------------------------------------------
# The module constants
# ---------------------------------------------------------------------------

WAVELENGTH_M = 1550e-9
R0_M = 0.10
DX_M = 0.01
PUPIL_D_M = 1.0
MASTER_SEED = 20260908

# ---- the Phase 1 grid ----
N_COLUMNS = (2, 4, 8, 16, 32)
L0_TIGHT_M = 2.56
L0_WIDE_M = 25.0
N_SMALL = 128
N_BIG = 512
N_COLUMNS_BIG = (2, 8)          # the production-size check, N = 512
SEEDS_PHASE1 = 8

# The spin-up is 2 outer scales. The record is 20 outer scales, except at
# N = 512 with L0 = 25 m. See the docstring, THE LAZY CHOICES.
SPINUP_L0 = 2.0
RECORD_L0 = 20.0
RECORD_ROWS_BIG_WIDE = 10000    # the CUT record, 4 L0 at N = 512, L0 = 25 m

# ---- the metric settings ----
# The band of metric 1. The plan: |rho_meas - rho_theory| < 0.02 at every lag
# past 0.5 L0.
RHO_BAND = 0.02
RHO_BAND_FROM_L0 = 0.5
RHO_LAG_MAX_L0 = 4.0
RHO_TABLE_L0 = (0.1, 0.25, 0.5, 1.0, 2.0, 4.0)

# The band of metric 3, the anisotropy D_ext / D_trans.
ANISO_BAND = 0.05
ANISO_R_M = (0.05, 0.10, 0.20, 0.40, 0.60, 1.20, 2.50)
ANISO_WINDOW_L0 = (0.1, 2.0)    # the plan window, in outer scales

# The band of metric 4, the chunk-mean variance of the row piston. The chunk is
# ONE outer scale of rows. A ratio band is loose, because the statistic reads
# 20 chunks per seed only.
PISTON_CHUNK_L0 = 1.0
PISTON_BAND = 0.30

# The band of metric 5.
VARIANCE_BAND = 0.10
TILT_BAND = 0.20

# The size above which the condition number of cov_zz is SKIPPED. The
# eigenvalue solve is O(n^3), and 4096 already costs about 15 s.
COND_MAX_SIZE = 4096

# ---- the Phase 2 grid ----
N_COLUMNS_WINDOW = (2, 16)
SEEDS_PHASE2 = 4
WINDOW_RECORD_L0 = 100.0
WINDOW_LENGTHS_L0 = (2.0, 5.0, 10.0, 20.0, 100.0)
WINDOW_COLUMN_BLOCK = 16        # keep the padded FFT of a long record small

# ---- the outputs ----
HERE = os.path.dirname(os.path.abspath(__file__))
ROWLAG_CSV = os.path.join(HERE, 'data', 'ncol_rowlag.csv')
PISTON_CSV = os.path.join(HERE, 'data', 'ncol_piston.csv')
COST_CSV = os.path.join(HERE, 'data', 'ncol_cost.csv')
WINDOW_CSV = os.path.join(HERE, 'data', 'ncol_window.csv')
ROWLAG_PNG = os.path.join(HERE, 'figures', 'ncol_rowlag.png')
ANISO_PNG = os.path.join(HERE, 'figures', 'ncol_anisotropy.png')
PISTON_PNG = os.path.join(HERE, 'figures', 'ncol_piston.png')
COST_PNG = os.path.join(HERE, 'figures', 'ncol_cost.png')
WINDOW_PNG = os.path.join(HERE, 'figures', 'ncol_window.png')
DPI = 150

RESULTS = []


def band_check(label, ok, detail=''):
    '''
    Print one PASS or FAIL line, and record it for the summary.

    Parameters:
        label : str
            The name of the band.
        ok : bool
            True for a pass.
        detail : str
            An optional note after the label.

    Returns:
        bool
            The same verdict.

    THIS FUNCTION RAISES NOTHING. A failed band is a result of the study.
    '''
    RESULTS.append((label, bool(ok)))
    print(f'   [{"PASS" if ok else "FAIL"}] {label}'
          f'{"   " + detail if detail else ""}')
    return bool(ok)


# ---------------------------------------------------------------------------
# The theory
# ---------------------------------------------------------------------------

def piston_covariance(k_rows, n_side, l0_m):
    '''
    Return the covariance of the ROW-MEAN phase against the row lag.

    Parameters:
        k_rows : array_like
            The row lag, in pixels.
        n_side : int
            The number of columns of a row.
        l0_m : float
            The outer scale L0 [m].

    Returns:
        numpy.ndarray
            C_m(k) [rad^2], one value per lag.

    formula:
        C_m(k) = (1/N^2) SUM_{j,j'} B(sqrt((k dx)^2 + ((j - j') dx)^2))
               = (1/N^2) SUM_{d=-(N-1)}^{N-1} (N - |d|)
                 B(sqrt((k dx)^2 + (d dx)^2))
    Source: the variance of an average of a homogeneous random field. Schmidt
    (2010), DOI 10.1117/3.866274, Ch. 3, Eq. (3.14), printed p. 47, defines
    B(r). The covariance itself is Assemat and Wilson,
    DOI 10.1364/OE.14.000988, Eq. (5), through `helpers.vk_covariance_closed`.

    THE WEIGHT. A transverse separation of d pixels occurs (N - |d|) times in
    the double sum, so the double sum collapses to one sum over d.

    WHY THE STUDY NEEDS THIS. The row piston is the near-unit-root mode of the
    recursion. Its variance is NOT B(0): the transverse average removes every
    scale below the grid width.
    '''
    k = np.atleast_1d(np.asarray(k_rows, dtype=float))
    n_side = int(n_side)
    d = np.arange(-(n_side - 1), n_side, dtype=float)
    weight = (n_side - np.abs(d)) / float(n_side) ** 2
    r_m = np.hypot(k[:, None] * DX_M, d[None, :] * DX_M)
    cov = helpers.vk_covariance_closed(r_m.ravel(), R0_M,
                                       l0_m).reshape(r_m.shape)
    return cov.dot(weight)


def chunk_variance_theory(cov_m, length):
    '''
    Return the theoretical variance of the MEAN of one chunk of the piston.

    Parameters:
        cov_m : numpy.ndarray
            C_m(k) for k = 0 to `length`. Use `piston_covariance`.
        length : int
            The chunk length, in rows.

    Returns:
        float
            The variance of the chunk mean [rad^2].

    formula:
        var = (1/L^2) SUM_{k=-(L-1)}^{L-1} (L - |k|) C_m(k)
    Source: the variance of the average of L samples of a stationary sequence.
    It is the same triangular (Bartlett) weight as `piston_covariance`, one
    dimension lower.
    '''
    length = int(length)
    k = np.arange(1, length, dtype=float)
    return float((length * cov_m[0]
                  + 2.0 * np.sum((length - k) * cov_m[1:length]))
                 / length ** 2)


def piston_spectrum_theory(cov_m, length):
    '''
    Return the EXPECTED periodogram of one chunk of the row piston.

    Parameters:
        cov_m : numpy.ndarray
            C_m(k) for k = 0 to `length`. Use `piston_covariance`.
        length : int
            The chunk length, in rows.

    Returns:
        numpy.ndarray
            The expected |FFT|^2 / L at the `length // 2 + 1` real
            frequencies [rad^2 per cycle per row].

    formula:
        E[ |FFT(m)|^2 / L ] (f) = SUM_{k=-(L-1)}^{L-1} (1 - |k|/L) C_m(k)
                                  exp(-2 pi i f k)
    Source: the Wiener-Khinchin relation for a FINITE record. The triangular
    factor (1 - |k|/L) is the Bartlett weight of the biased lag estimator, and
    it is exact for the expectation of the raw periodogram. Schmidt (2010),
    DOI 10.1117/3.866274, Ch. 3, Eq. (3.16), printed p. 48, states the same
    transform pair for a continuous field.

    WHY THE EXPECTATION AND NOT THE PSD. A raw periodogram of L samples is
    NOT an unbiased estimate of the continuous spectrum. It is unbiased for the
    expression above. So the script compares like with like, and it needs no
    window and no Welch average.
    '''
    length = int(length)
    k = np.arange(length, dtype=float)
    weighted = (1.0 - k / length) * cov_m[:length]
    # The two-sided sum is symmetric. So it reads twice the positive lags, and
    # it then removes the double-counted zero lag.
    return np.fft.rfft(weighted).real * 2.0 - weighted[0]


# ---------------------------------------------------------------------------
# The extrusion
# ---------------------------------------------------------------------------

def build_screen(n_side, l0_m, n_columns, seed):
    '''
    Build one `PhaseScreenVonKarman`, and capture the aotools guards.

    Parameters:
        n_side : int
            The grid side, in pixels.
        l0_m : float
            The outer scale L0 [m].
        n_columns : int
            The number of stencil rows of the recursion.
        seed : int
            The integer seed of the initial frame.

    Returns:
        tuple
            (screen, guards, setup_seconds). `guards` is a dict with the
            lstsq-fallback flag, the minimum eigenvalue of BB^T, the condition
            number of cov_zz, and the size of cov_zz.

    THE STDOUT CAPTURE. `PhaseScreen.makeAMatrix` PRINTS "Cholesky solve
    failed..." and it falls back to a truncated pseudo-inverse. That fallback
    changes the statistics silently, so the study must know when it runs. The
    class emits a print, not a warning, so the only guard is a stdout capture.

    THE TWO EIGENVALUE GUARDS. `makeBMatrix` takes sqrt(W) of the singular
    values of BB^T = cov_xx - A cov_zx, so a negative eigenvalue becomes its
    own absolute value. The minimum eigenvalue tells if that happens. The
    condition number of cov_zz tells how much precision the inverse needs; it
    is the quantity that broke the float32 kernel of aotools 1.0.7.
    '''
    t0 = time.time()
    buffer = io.StringIO()
    with contextlib.redirect_stdout(buffer):
        screen = PhaseScreenVonKarman(int(n_side), DX_M, R0_M, float(l0_m),
                                      random_seed=int(seed),
                                      n_columns=int(n_columns))
    setup = time.time() - t0

    bbt = screen.cov_mat_xx - screen.A_mat.dot(screen.cov_mat_zx)
    min_eig = float(np.linalg.eigvalsh(bbt).min())

    # The eigenvalue solve can fail on a bad LAPACK build. A guard is a
    # diagnostic, so a failure must not stop the study.
    size = int(screen.cov_mat_zz.shape[0])
    cond = np.nan
    if size <= COND_MAX_SIZE:
        try:
            eigenvalues = np.linalg.eigvalsh(screen.cov_mat_zz)
            cond = float(abs(eigenvalues).max() / abs(eigenvalues).min())
        except np.linalg.LinAlgError:
            cond = np.nan

    guards = dict(fallback='Cholesky' in buffer.getvalue(), min_eig=min_eig,
                  cond=cond, size=size)
    return screen, guards, setup


def extrude(screen, seed, spinup, record):
    '''
    Re-seed one screen object, spin it up, and store the fresh rows.

    Parameters:
        screen : PhaseScreenVonKarman
            The object. This function CHANGES its state.
        seed : int
            The integer seed of the initial frame.
        spinup : int
            The number of add_row calls before the record starts.
        record : int
            The number of recorded rows.

    Returns:
        tuple
            (rows, frame, add_row_seconds). `rows` is a `record` by N float32
            array. `frame` is the last visible screen [rad].

    WHY THE RE-SEED. The A and the B matrix depend on the geometry only, not on
    the seed. `make_initial_screen` rebuilds the initial frame and the random
    generator from `screen.random_seed`, and it touches nothing else. So one
    setup serves every seed of a cell.

    THE FRESH ROW. `add_row` PREPENDS the new row on axis 0, so the fresh row
    is `screen._scrn[0, :]`. The public `.scrn` property crops to the requested
    size, and it holds the same row 0.
    '''
    screen.random_seed = int(seed)
    screen.make_initial_screen()

    n_side = int(screen.requested_nx_size)
    rows = np.empty((int(record), n_side), dtype=np.float32)
    t0 = time.time()
    for _ in range(int(spinup)):
        screen.add_row()
    for i in range(int(record)):
        screen.add_row()
        rows[i] = screen._scrn[0, :n_side]
    elapsed = time.time() - t0
    frame = np.array(screen.scrn, dtype=float)
    return rows, frame, elapsed / float(spinup + record)


def lag_sums_blocked(rows, k_max, block=WINDOW_COLUMN_BLOCK):
    '''
    Return the row-lag products of a long record, in column blocks.

    Parameters:
        rows : numpy.ndarray
            An M by N array of the fresh rows.
        k_max : int
            The largest lag, in rows.
        block : int
            The number of columns of one block.

    Returns:
        numpy.ndarray
            `k_max` + 1 values. Entry k is the mean over the N columns of
            SUM_i phi_i phi_{i+k}. Divide by (M - k) to get C(k).

    WHY THE BLOCKS. `extrusion_stationarity.row_lag_sums` pads the record to
    the next power of two and it transforms every column at one time. A 250000
    row record of 128 columns then needs a 524288 by 128 complex array, which
    is 1 GB. A column block keeps that array small. The columns are equal in
    number, so the block means average with equal weight.

    THE MODULE CONSTANT. `row_lag_sums` truncates at
    `extrusion_stationarity.K_MAX`. This function RAISES that constant, and it
    restores it. The estimator itself is unchanged.
    '''
    rows = np.asarray(rows)
    n_side = rows.shape[1]
    block = int(block)
    keep = es.K_MAX
    es.K_MAX = int(k_max)
    try:
        total = np.zeros(int(k_max) + 1)
        count = 0
        for start in range(0, n_side, block):
            total += es.row_lag_sums(rows[:, start:start + block])
            count += 1
        return total / count
    finally:
        es.K_MAX = keep


# ---------------------------------------------------------------------------
# One Phase 1 cell
# ---------------------------------------------------------------------------

def run_cell(n_side, l0_m, n_columns, seeds):
    '''
    Run one (N, L0, n_columns) cell over every seed, and reduce it.

    Parameters:
        n_side : int
            The grid side, in pixels.
        l0_m : float
            The outer scale L0 [m].
        n_columns : int
            The stencil depth of the recursion.
        seeds : list
            The integer seeds.

    Returns:
        dict
            Every measured quantity of the cell. See the keys below.

    THE ONE PASS. A cell extrudes each seed one time, and it harvests every
    metric from the SAME rows: the extrusion-axis correlation, the transverse
    covariance of square blocks of the record, the row piston, the frame
    variance and the frame tilt.
    '''
    spinup = int(round(SPINUP_L0 * l0_m / DX_M))
    if n_side == N_BIG and l0_m == L0_WIDE_M:
        record = RECORD_ROWS_BIG_WIDE
    else:
        record = int(round(RECORD_L0 * l0_m / DX_M))
    k_max = min(int(round(RHO_LAG_MAX_L0 * l0_m / DX_M)), record // 2)
    chunk = int(round(PISTON_CHUNK_L0 * l0_m / DX_M))
    k_trans = n_side // 2

    screen, guards, setup = build_screen(n_side, l0_m, n_columns, seeds[0])

    rho_seed = np.zeros((len(seeds), k_max + 1))
    b_trans = np.zeros(k_trans + 1)
    chunk_means = []
    periodogram = np.zeros(chunk // 2 + 1)
    n_periodogram = 0
    frame_var = np.zeros(len(seeds))
    tilt = np.zeros((len(seeds), 2))
    add_row_s = 0.0

    mask = helpers.pupil_mask(n_side, DX_M, PUPIL_D_M)

    for index, seed in enumerate(seeds):
        rows, frame, per_row = extrude(screen, seed, spinup, record)
        add_row_s += per_row

        # ---- metric 1: the extrusion axis ----
        sums = lag_sums_blocked(rows, k_max)
        c_est = sums / (record - np.arange(k_max + 1)).astype(float)
        rho_seed[index] = c_est / c_est[0]

        # ---- metric 2: the transverse axis, on SQUARE blocks ----
        n_blocks = min(20, record // n_side)
        blocks = [rows[b * n_side:(b + 1) * n_side].astype(float)
                  for b in range(n_blocks)]
        b_trans += axis_covariance(blocks, k_trans)[1] / len(seeds)
        del blocks

        # ---- metric 4: the row piston ----
        piston = rows.mean(axis=1).astype(float)
        n_chunks = piston.size // chunk
        pieces = piston[:n_chunks * chunk].reshape(n_chunks, chunk)
        chunk_means.append(pieces.mean(axis=1))
        periodogram += (np.abs(np.fft.rfft(pieces, axis=1)) ** 2).sum(axis=0) \
            / chunk
        n_periodogram += n_chunks

        # ---- metric 5: the frame ----
        frame_var[index] = float(np.mean(frame * frame))
        _, _, a2, a3 = helpers.zernike_tilt(frame, mask, DX_M, WAVELENGTH_M)
        tilt[index] = (a2, a3)
        del rows, frame

    del screen

    chunk_means = np.concatenate(chunk_means)
    periodogram /= n_periodogram

    # ---- the theory ----
    k_axis = np.arange(k_max + 1)
    b_theory = helpers.vk_covariance_closed(k_axis * DX_M, R0_M, l0_m)
    rho_theory = b_theory / b_theory[0]
    cov_m = piston_covariance(np.arange(chunk + 1), n_side, l0_m)
    d_eff = helpers.mask_diameter(mask, DX_M)
    psd = lambda f: von_karman_phase_psd(f, R0_M, l0_m)

    return dict(
        n_side=n_side, l0_m=l0_m, n_columns=n_columns, spinup=spinup,
        record=record, chunk=chunk, k_max=k_max, k_trans=k_trans,
        rho=rho_seed.mean(axis=0),
        rho_se=rho_seed.std(axis=0, ddof=1) / np.sqrt(len(seeds)),
        rho_theory=rho_theory,
        b_trans=b_trans,
        b_theory_trans=helpers.vk_covariance_closed(
            np.arange(k_trans + 1) * DX_M, R0_M, l0_m),
        chunk_var=float(np.mean(chunk_means ** 2)),
        chunk_var_se=float(np.mean(chunk_means ** 2)
                           * np.sqrt(2.0 / chunk_means.size)),
        chunk_var_theory=chunk_variance_theory(cov_m, chunk),
        n_chunks=int(chunk_means.size),
        periodogram=periodogram,
        periodogram_theory=piston_spectrum_theory(cov_m, chunk),
        frame_var=float(frame_var.mean()),
        frame_var_theory=float(b_theory[0]),
        tilt_var=float(np.mean(tilt ** 2)),
        tilt_var_theory=float(helpers.tilt_filter_variance(psd, d_eff)),
        guards=guards, setup_s=setup,
        add_row_ms=1e3 * add_row_s / len(seeds))


# ---------------------------------------------------------------------------
# One Phase 2 cell
# ---------------------------------------------------------------------------

def run_window_cell(n_side, l0_m, n_columns, seeds):
    '''
    Measure the APPARENT rho(1 L0) against the analysis window length.

    Parameters:
        n_side : int
            The grid side, in pixels.
        l0_m : float
            The outer scale L0 [m].
        n_columns : int
            The stencil depth of the recursion.
        seeds : list
            The integer seeds.

    Returns:
        dict
            The window lengths, the mean rho(1 L0) of each one, its standard
            error, and the theory.

    THE QUESTION. A short window cannot see a slow mode. The row piston of a
    von Karman screen wanders over many outer scales, so a short window reads
    that wander as a CONSTANT offset, and rho(k) then sits on a floor. This
    function cuts ONE long record into windows of several lengths, and it reads
    rho at a lag of one outer scale in each one.

    NO MEAN SUBTRACTION, in a window either. A mean subtraction inside a window
    is exactly the operation that would hide the effect this function looks
    for.
    '''
    spinup = int(round(SPINUP_L0 * l0_m / DX_M))
    record = int(round(WINDOW_RECORD_L0 * l0_m / DX_M))
    k_l0 = int(round(l0_m / DX_M))

    screen, guards, setup = build_screen(n_side, l0_m, n_columns, seeds[0])

    values = dict((w, []) for w in WINDOW_LENGTHS_L0)
    for seed in seeds:
        rows, _, _ = extrude(screen, seed, spinup, record)
        for w in WINDOW_LENGTHS_L0:
            length = int(round(w * l0_m / DX_M))
            for start in range(0, record - length + 1, length):
                sums = lag_sums_blocked(rows[start:start + length], k_l0)
                counts = (length - np.arange(k_l0 + 1)).astype(float)
                c_est = sums / counts
                values[w].append(float(c_est[k_l0] / c_est[0]))
        del rows
    del screen

    b_theory = helpers.vk_covariance_closed([0.0, l0_m], R0_M, l0_m)
    means, errors, counts = [], [], []
    for w in WINDOW_LENGTHS_L0:
        sample = np.asarray(values[w], dtype=float)
        means.append(float(sample.mean()))
        errors.append(float(sample.std(ddof=1) / np.sqrt(sample.size))
                      if sample.size > 1 else np.nan)
        counts.append(int(sample.size))
    return dict(n_side=n_side, l0_m=l0_m, n_columns=n_columns, record=record,
                windows=list(WINDOW_LENGTHS_L0), rho=np.array(means),
                rho_se=np.array(errors), counts=counts,
                rho_theory=float(b_theory[1] / b_theory[0]),
                guards=guards, setup_s=setup)


# ---------------------------------------------------------------------------
# The figures
# ---------------------------------------------------------------------------

def cell_label(cell):
    '''Return the short name of one cell, for a legend.'''
    return f'n_col {cell["n_columns"]:2d}, N {cell["n_side"]}'


def draw_rowlag(cells):
    '''Draw rho(k) per n_columns, one panel per (L0, N).'''
    panels = [(L0_TIGHT_M, N_SMALL), (L0_WIDE_M, N_SMALL),
              (L0_TIGHT_M, N_BIG), (L0_WIDE_M, N_BIG)]
    fig, axes = plt.subplots(2, 2, figsize=(13.6, 9.2),
                             constrained_layout=True)
    for axis, (l0_m, n_side) in zip(axes.ravel(), panels):
        chosen = [c for c in cells
                  if c['l0_m'] == l0_m and c['n_side'] == n_side]
        if not chosen:
            axis.set_visible(False)
            continue
        reference = chosen[0]
        lags = np.arange(reference['k_max'] + 1)[1:] * DX_M / l0_m
        axis.plot(lags, reference['rho_theory'][1:], color='black',
                  linewidth=2.6, label='von Karman theory, Eq. (5)')
        for cell in chosen:
            lag = np.arange(cell['k_max'] + 1)[1:] * DX_M / l0_m
            axis.plot(lag, cell['rho'][1:], linewidth=1.5,
                      label=f'n_columns = {cell["n_columns"]}')
        axis.axvline(RHO_BAND_FROM_L0, color='tab:green', linestyle=':',
                     linewidth=1.6, label=f'band from {RHO_BAND_FROM_L0} L0')
        axis.set_xscale('log')
        axis.set_xlabel('row lag, outer scales')
        axis.set_ylabel('rho(k) = C(k) / C(0)')
        axis.set_title(f'L0 = {l0_m:.2f} m, N = {n_side}, '
                       f'{reference["record"]} recorded rows', fontsize=10)
        axis.grid(alpha=0.3, which='both')
        axis.legend(fontsize=8)
    fig.suptitle('THE EXTRUSION-AXIS CORRELATION against n_columns. '
                 f'r0 = {R0_M * 1e2:.0f} cm, dx = {DX_M * 1e2:.0f} cm, '
                 f'{SEEDS_PHASE1} seeds.', fontsize=12)
    fig.savefig(ROWLAG_PNG, dpi=DPI)
    plt.close(fig)


def draw_anisotropy(cells, table):
    '''Draw D_ext(r) / D_trans(r) against the separation.'''
    fig, axes = plt.subplots(1, 2, figsize=(13.6, 5.4),
                             constrained_layout=True)
    for axis, l0_m in zip(axes, (L0_TIGHT_M, L0_WIDE_M)):
        for cell in cells:
            if cell['l0_m'] != l0_m:
                continue
            rows = [row for row in table if row['cell'] is cell]
            if not rows:
                continue
            r_m = np.array([row['r_m'] for row in rows])
            ratio = np.array([row['ratio'] for row in rows])
            style = '-' if cell['n_side'] == N_SMALL else '--'
            axis.plot(r_m, ratio, style, marker='o', markersize=4,
                      linewidth=1.5, label=cell_label(cell))
        axis.axhline(1.0, color='black', linewidth=2.2, label='isotropic')
        axis.axhspan(1.0 - ANISO_BAND, 1.0 + ANISO_BAND, color='0.85',
                     zorder=0, label=f'band, {ANISO_BAND:.2f}')
        axis.set_xscale('log')
        axis.set_xlabel('separation r, m')
        axis.set_ylabel('D_ext(r) / D_trans(r)')
        axis.set_title(f'L0 = {l0_m:.2f} m. A value above 1 means the '
                       'extrusion axis is ROUGHER.', fontsize=10)
        axis.grid(alpha=0.3, which='both')
        axis.legend(fontsize=8)
    fig.suptitle('THE STRUCTURE-FUNCTION ANISOTROPY. The transverse axis is '
                 'the control: the recursion does not touch it.', fontsize=12)
    fig.savefig(ANISO_PNG, dpi=DPI)
    plt.close(fig)


def draw_piston(cells):
    '''Draw the chunk-mean variance ratio and the piston spectrum.'''
    fig, axes = plt.subplots(1, 2, figsize=(13.6, 5.4),
                             constrained_layout=True)
    a0, a1 = axes

    for l0_m in (L0_TIGHT_M, L0_WIDE_M):
        for n_side in (N_SMALL, N_BIG):
            chosen = sorted([c for c in cells if c['l0_m'] == l0_m
                             and c['n_side'] == n_side],
                            key=lambda c: c['n_columns'])
            if not chosen:
                continue
            x = [c['n_columns'] for c in chosen]
            y = [c['chunk_var'] / c['chunk_var_theory'] for c in chosen]
            error = [c['chunk_var_se'] / c['chunk_var_theory'] for c in chosen]
            a0.errorbar(x, y, yerr=error, marker='o', capsize=3.0,
                        linewidth=1.6,
                        label=f'L0 = {l0_m:.2f} m, N = {n_side}')
    a0.axhline(1.0, color='black', linewidth=2.2, label='theory')
    a0.axhspan(1.0 - PISTON_BAND, 1.0 + PISTON_BAND, color='0.85', zorder=0,
               label=f'band, {PISTON_BAND:.2f}')
    a0.set_xscale('log', base=2)
    a0.set_xticks(N_COLUMNS)
    a0.set_xticklabels([str(v) for v in N_COLUMNS])
    a0.set_xlabel('n_columns')
    a0.set_ylabel('measured / theory')
    a0.set_title(f'the CHUNK-MEAN variance of the row piston, chunk = '
                 f'{PISTON_CHUNK_L0:.0f} L0', fontsize=10)
    a0.grid(alpha=0.3)
    a0.legend(fontsize=8)

    chosen = sorted([c for c in cells if c['l0_m'] == L0_WIDE_M
                     and c['n_side'] == N_SMALL],
                    key=lambda c: c['n_columns'])
    if chosen:
        freq = np.fft.rfftfreq(chosen[0]['chunk'])
        a1.loglog(freq[1:], chosen[0]['periodogram_theory'][1:], color='black',
                  linewidth=2.6, label='theory, the expected periodogram')
        for cell in chosen:
            a1.loglog(freq[1:], cell['periodogram'][1:], linewidth=1.4,
                      label=f'n_columns = {cell["n_columns"]}')
    a1.set_xlabel('frequency, cycles per row')
    a1.set_ylabel('|FFT(m)|^2 / L, rad^2')
    a1.set_title(f'the piston SPECTRUM, L0 = {L0_WIDE_M:.0f} m, N = {N_SMALL}',
                 fontsize=10)
    a1.grid(alpha=0.3, which='both')
    a1.legend(fontsize=8)

    fig.suptitle('THE ROW PISTON m(t), the near-unit-root mode of the '
                 'recursion.', fontsize=12)
    fig.savefig(PISTON_PNG, dpi=DPI)
    plt.close(fig)


def draw_cost(cells):
    '''Draw the setup time, the add_row time and the size of cov_zz.'''
    fig, axes = plt.subplots(1, 3, figsize=(15.0, 4.8),
                             constrained_layout=True)
    for n_side in (N_SMALL, N_BIG):
        chosen = sorted([c for c in cells if c['n_side'] == n_side
                         and c['l0_m'] == L0_WIDE_M],
                        key=lambda c: c['n_columns'])
        if not chosen:
            continue
        x = [c['n_columns'] for c in chosen]
        axes[0].plot(x, [c['setup_s'] for c in chosen], marker='o',
                     linewidth=1.6, label=f'N = {n_side}')
        axes[1].plot(x, [c['add_row_ms'] for c in chosen], marker='o',
                     linewidth=1.6, label=f'N = {n_side}')
        axes[2].plot(x, [c['guards']['size'] ** 2 * 8.0 / 2 ** 20
                         for c in chosen], marker='o', linewidth=1.6,
                     label=f'N = {n_side}')
    for axis, name in zip(axes, ('setup wall time, s',
                                 'add_row time, ms',
                                 'cov_zz, MB of float64')):
        axis.set_xscale('log', base=2)
        axis.set_yscale('log')
        axis.set_xticks(N_COLUMNS)
        axis.set_xticklabels([str(v) for v in N_COLUMNS])
        axis.set_xlabel('n_columns')
        axis.set_ylabel(name)
        axis.grid(alpha=0.3, which='both')
        axis.legend(fontsize=8)
    fig.suptitle('THE COST OF n_columns. The stencil holds n_columns N points, '
                 'so cov_zz grows as (n_columns N)^2.', fontsize=12)
    fig.savefig(COST_PNG, dpi=DPI)
    plt.close(fig)


def draw_window(window_cells):
    '''Draw the apparent rho(1 L0) against the analysis window length.'''
    fig, axes = plt.subplots(1, 2, figsize=(13.6, 5.4),
                             constrained_layout=True)
    for axis, l0_m in zip(axes, (L0_TIGHT_M, L0_WIDE_M)):
        chosen = [c for c in window_cells if c['l0_m'] == l0_m]
        if not chosen:
            axis.set_visible(False)
            continue
        for cell in chosen:
            axis.errorbar(cell['windows'], cell['rho'], yerr=cell['rho_se'],
                          marker='o', capsize=3.0, linewidth=1.6,
                          label=f'n_columns = {cell["n_columns"]}')
        axis.axhline(chosen[0]['rho_theory'], color='black', linewidth=2.4,
                     label=f'theory, rho(1 L0) = '
                           f'{chosen[0]["rho_theory"]:.4f}')
        axis.set_xscale('log')
        axis.set_xlabel('analysis window, outer scales')
        axis.set_ylabel('apparent rho(1 L0)')
        axis.set_title(f'L0 = {l0_m:.2f} m, N = {N_SMALL}, '
                       f'{WINDOW_RECORD_L0:.0f} L0 of record', fontsize=10)
        axis.grid(alpha=0.3, which='both')
        axis.legend(fontsize=8)
    fig.suptitle('PHASE 2. Is the apparent over-correlation a RECORD LENGTH, '
                 'not a Markov memory?', fontsize=12)
    fig.savefig(WINDOW_PNG, dpi=DPI)
    plt.close(fig)


# ---------------------------------------------------------------------------
# The CSV writer
# ---------------------------------------------------------------------------

def write_csv(path, header, rows):
    '''Write one table. The csv module needs newline="".'''
    with open(path, 'w', newline='', encoding='utf-8') as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)


# ---------------------------------------------------------------------------

def main():
    t_start = time.time()
    seeds1 = helpers.spawn_seeds(MASTER_SEED, SEEDS_PHASE1)
    seeds2 = helpers.spawn_seeds(MASTER_SEED + 1, SEEDS_PHASE2)

    print('=' * 78)
    print('PLAN_n_columns.md Phases 1 and 2: the n_columns sweep')
    print('=' * 78)
    print(f'  wavelength           {WAVELENGTH_M * 1e9:12.1f} nm')
    print(f'  r0                   {R0_M * 1e2:12.1f} cm, dx = '
          f'{DX_M * 1e2:.1f} cm')
    print(f'  n_columns            {str(N_COLUMNS):>12}')
    print(f'  outer scales         {L0_TIGHT_M:12.2f} m and {L0_WIDE_M:.2f} m')
    print(f'  grids                {N_SMALL:12d} px, plus {N_BIG} px at '
          f'n_columns {N_COLUMNS_BIG}')
    print(f'  seeds                {SEEDS_PHASE1:12d} (Phase 1), '
          f'{SEEDS_PHASE2} (Phase 2)')
    print(f'  spin-up              {SPINUP_L0:12.1f} L0, record '
          f'{RECORD_L0:.0f} L0')
    print(f'  master seed          {MASTER_SEED:12d}')
    print('')
    print('  THE CUT. At N = 512 with L0 = 25 m the record is '
          f'{RECORD_ROWS_BIG_WIDE} rows')
    print(f'  ({RECORD_ROWS_BIG_WIDE * DX_M / L0_WIDE_M:.0f} L0), not '
          f'{RECORD_L0:.0f} L0. One add_row costs 1.3 ms there against 0.02 ms')
    print('  at N = 128. The seed count is NOT cut.')

    # ---- Phase 1 ----
    plan = [(N_SMALL, l0_m, n_columns)
            for l0_m in (L0_TIGHT_M, L0_WIDE_M)
            for n_columns in N_COLUMNS]
    plan += [(N_BIG, l0_m, n_columns)
             for l0_m in (L0_TIGHT_M, L0_WIDE_M)
             for n_columns in N_COLUMNS_BIG]

    cells = []
    print('')
    print(f'  running {len(plan)} Phase 1 cells ...')
    for n_side, l0_m, n_columns in plan:
        t0 = time.time()
        cell = run_cell(n_side, l0_m, n_columns, seeds1)
        cells.append(cell)
        print(f'    N {n_side:4d}  L0 {l0_m:6.2f} m  n_columns {n_columns:3d}'
              f'   {time.time() - t0:7.1f} s   '
              f'(total {time.time() - t_start:7.1f} s)')

    # ---- 1. the extrusion-axis correlation ----
    print('')
    print('1. THE EXTRUSION-AXIS CORRELATION rho(k) = C(k) / C(0)')
    print('   theory: B(k dx) / B(0), Assemat and Wilson Eq. (5), '
          'DOI 10.1364/OE.14.000988')
    print(f'   band: |rho_meas - rho_theory| < {RHO_BAND:.2f} at every lag '
          f'past {RHO_BAND_FROM_L0:.1f} L0')
    print('')
    head = f'   {"N":>5}{"L0, m":>8}{"n_col":>7}'
    for multiple in RHO_TABLE_L0:
        head += f'{"d@" + str(multiple) + "L0":>11}'
    head += f'{"worst":>9}{"at, L0":>9}{"2 SE":>9}{"verdict":>9}'
    print(head)
    rowlag_rows = []
    for cell in cells:
        lags = np.arange(cell['k_max'] + 1)
        deviation = cell['rho'] - cell['rho_theory']
        line = (f'   {cell["n_side"]:>5}{cell["l0_m"]:>8.2f}'
                f'{cell["n_columns"]:>7}')
        for multiple in RHO_TABLE_L0:
            k = int(round(multiple * cell['l0_m'] / DX_M))
            line += (f'{deviation[k]:>+11.4f}' if k <= cell['k_max']
                     else f'{"n/a":>11}')
        band = lags >= RHO_BAND_FROM_L0 * cell['l0_m'] / DX_M
        if np.any(band):
            index = int(lags[band][np.argmax(np.abs(deviation[band]))])
            worst = float(deviation[index])
            ok = float(np.abs(deviation[band]).max()) < RHO_BAND
            line += (f'{worst:>+9.4f}{index * DX_M / cell["l0_m"]:>9.2f}'
                     f'{2.0 * cell["rho_se"][index]:>9.4f}'
                     f'{"pass" if ok else "FAIL":>9}')
        else:
            ok = True
            line += f'{"n/a":>9}{"n/a":>9}{"n/a":>9}{"n/a":>9}'
        print(line)
        cell['rho_ok'] = ok

        # The CSV holds a log-spaced subset of the lags, not every lag.
        subset = np.unique(np.round(np.geomspace(
            1, max(cell['k_max'], 1), 50)).astype(int))
        for k in subset:
            rowlag_rows.append([cell['n_columns'], f'{cell["l0_m"]:.2f}',
                                cell['n_side'], 'extrusion', int(k),
                                f'{k * DX_M:.4f}', f'{cell["rho"][k]:.6f}',
                                f'{cell["rho_se"][k]:.6f}',
                                f'{cell["rho_theory"][k]:.6f}'])
        rho_trans = cell['b_trans'] / cell['b_trans'][0]
        theory_trans = cell['b_theory_trans'] / cell['b_theory_trans'][0]
        for k in range(1, cell['k_trans'] + 1, 4):
            rowlag_rows.append([cell['n_columns'], f'{cell["l0_m"]:.2f}',
                                cell['n_side'], 'transverse', int(k),
                                f'{k * DX_M:.4f}', f'{rho_trans[k]:.6f}', '',
                                f'{theory_trans[k]:.6f}'])
    for cell in cells:
        band_check(f'1. rho(k), N {cell["n_side"]}, L0 {cell["l0_m"]:.2f} m, '
                   f'n_columns {cell["n_columns"]}', cell['rho_ok'])

    # ---- 2. the transverse control ----
    print('')
    print('2. THE TRANSVERSE AXIS, THE CONTROL. The recursion does not touch')
    print('   it, so it must hold the theory. The band is the same 0.02 on')
    print(f'   rho, out to half the grid side.')
    print('')
    print(f'   {"N":>5}{"L0, m":>8}{"n_col":>7}{"reach, m":>10}'
          f'{"worst d":>10}{"at r, m":>10}{"verdict":>9}')
    for cell in cells:
        rho_trans = cell['b_trans'] / cell['b_trans'][0]
        theory = cell['b_theory_trans'] / cell['b_theory_trans'][0]
        deviation = rho_trans - theory
        index = int(np.argmax(np.abs(deviation[1:]))) + 1
        ok = float(np.abs(deviation[1:]).max()) < RHO_BAND
        print(f'   {cell["n_side"]:>5}{cell["l0_m"]:>8.2f}'
              f'{cell["n_columns"]:>7}{cell["k_trans"] * DX_M:>10.2f}'
              f'{deviation[index]:>+10.4f}{index * DX_M:>10.3f}'
              f'{"pass" if ok else "FAIL":>9}')
        cell['trans_ok'] = ok
    for cell in cells:
        band_check(f'2. transverse rho, N {cell["n_side"]}, '
                   f'L0 {cell["l0_m"]:.2f} m, n_columns {cell["n_columns"]}',
                   cell['trans_ok'])

    # ---- 3. the anisotropy ----
    # D(r) = 2 [B(0) - B(r)]. Schmidt (2010), DOI 10.1117/3.866274, Ch. 3,
    # Eq. (3.16), printed p. 48. The extrusion axis reads the row record, and
    # the transverse axis reads the square blocks of the same record.
    print('')
    print('3. THE STRUCTURE-FUNCTION ANISOTROPY D_ext(r) / D_trans(r)')
    print(f'   band: inside {100.0 * ANISO_BAND:.0f} percent of 1.0.')
    print('   The plan asks for r = 0.1 to 2 L0. The TRANSVERSE axis cannot')
    print('   reach past half the grid side, so the test reads every reachable')
    print('   separation and it flags the plan window.')
    print('')
    print(f'   {"N":>5}{"L0, m":>8}{"n_col":>7}{"r, m":>9}{"r / L0":>9}'
          f'{"D_ext":>11}{"D_trans":>11}{"ratio":>9}{"in plan":>9}')
    aniso_table = []
    for cell in cells:
        worst, worst_r = 0.0, np.nan
        for r_m in ANISO_R_M:
            k = int(round(r_m / DX_M))
            if k > cell['k_trans'] or k > cell['k_max']:
                continue
            d_ext = 2.0 * (cell['rho'][0] - cell['rho'][k])
            d_tra = 2.0 * (1.0 - cell['b_trans'][k] / cell['b_trans'][0])
            ratio = d_ext / d_tra
            multiple = r_m / cell['l0_m']
            inside = (ANISO_WINDOW_L0[0] <= multiple <= ANISO_WINDOW_L0[1])
            aniso_table.append(dict(cell=cell, r_m=r_m, ratio=ratio,
                                    d_ext=d_ext, d_tra=d_tra, inside=inside))
            print(f'   {cell["n_side"]:>5}{cell["l0_m"]:>8.2f}'
                  f'{cell["n_columns"]:>7}{r_m:>9.2f}{multiple:>9.3f}'
                  f'{d_ext:>11.4f}{d_tra:>11.4f}{ratio:>9.4f}'
                  f'{"yes" if inside else "no":>9}')
            if abs(ratio - 1.0) > worst:
                worst, worst_r = abs(ratio - 1.0), r_m
        cell['aniso_worst'] = worst
        cell['aniso_r'] = worst_r
    print('')
    for cell in cells:
        band_check(f'3. anisotropy, N {cell["n_side"]}, '
                   f'L0 {cell["l0_m"]:.2f} m, n_columns {cell["n_columns"]}',
                   cell['aniso_worst'] < ANISO_BAND,
                   f'worst {100.0 * cell["aniso_worst"]:.1f} percent at '
                   f'r = {cell["aniso_r"]:.2f} m')

    # ---- 4. the row piston ----
    print('')
    print('4. THE ROW PISTON m(t), THE MEAN OF EACH NEW ROW')
    print(f'   the chunk is {PISTON_CHUNK_L0:.0f} L0 of rows. The theory is '
          f'the transverse-averaged')
    print('   covariance of Assemat and Wilson Eq. (5), with the triangular')
    print(f'   weight of a finite average. band: inside '
          f'{100.0 * PISTON_BAND:.0f} percent.')
    print('')
    print(f'   {"N":>5}{"L0, m":>8}{"n_col":>7}{"chunk":>8}{"n":>7}'
          f'{"var meas":>12}{"var theory":>12}{"ratio":>9}{"SE":>9}')
    piston_rows = []
    for cell in cells:
        ratio = cell['chunk_var'] / cell['chunk_var_theory']
        print(f'   {cell["n_side"]:>5}{cell["l0_m"]:>8.2f}'
              f'{cell["n_columns"]:>7}{cell["chunk"]:>8}{cell["n_chunks"]:>7}'
              f'{cell["chunk_var"]:>12.4f}{cell["chunk_var_theory"]:>12.4f}'
              f'{ratio:>9.4f}{cell["chunk_var_se"]:>9.4f}')
        cell['piston_ok'] = abs(ratio - 1.0) < PISTON_BAND
        piston_rows.append([cell['n_columns'], f'{cell["l0_m"]:.2f}',
                            cell['n_side'], 'chunk_var', '',
                            f'{cell["chunk_var"]:.6e}',
                            f'{cell["chunk_var_theory"]:.6e}',
                            f'{cell["chunk_var_se"]:.6e}'])
        freq = np.fft.rfftfreq(cell['chunk'])
        subset = np.unique(np.round(np.geomspace(
            1, freq.size - 1, 40)).astype(int))
        for i in subset:
            piston_rows.append([cell['n_columns'], f'{cell["l0_m"]:.2f}',
                                cell['n_side'], 'psd', f'{freq[i]:.6e}',
                                f'{cell["periodogram"][i]:.6e}',
                                f'{cell["periodogram_theory"][i]:.6e}', ''])
    print('')
    for cell in cells:
        band_check(f'4. piston chunk variance, N {cell["n_side"]}, '
                   f'L0 {cell["l0_m"]:.2f} m, n_columns {cell["n_columns"]}',
                   cell['piston_ok'],
                   f'ratio {cell["chunk_var"] / cell["chunk_var_theory"]:.3f}')

    print('')
    print('   THE PISTON SPECTRUM, measured over theory, at four frequencies.')
    print('   A ratio above 1 at the LOW frequencies is the slow wander that')
    print('   the plan describes.')
    print(f'   {"N":>5}{"L0, m":>8}{"n_col":>7}{"f=1/L":>10}{"f=4/L":>10}'
          f'{"f=16/L":>10}{"f=L/4":>10}')
    for cell in cells:
        line = (f'   {cell["n_side"]:>5}{cell["l0_m"]:>8.2f}'
                f'{cell["n_columns"]:>7}')
        for index in (1, 4, 16, cell['chunk'] // 4):
            ratio = (cell['periodogram'][index]
                     / cell['periodogram_theory'][index])
            line += f'{ratio:>10.3f}'
        print(line)

    # ---- 5. the frame variance and the tilt ----
    print('')
    print('5. THE FRAME VARIANCE AND THE Z-TILT over a '
          f'{PUPIL_D_M:.1f} m pupil')
    print(f'   band: variance inside {100.0 * VARIANCE_BAND:.0f} percent, '
          f'tilt inside {100.0 * TILT_BAND:.0f} percent.')
    print('   The variance carries the SPIN-UP, not only the recursion: a '
          'frame')
    print('   of one grid side holds far less than one L0 at L0 = 25 m.')
    print('')
    print(f'   {"N":>5}{"L0, m":>8}{"n_col":>7}{"<phi^2>":>11}{"theory":>11}'
          f'{"ratio":>9}{"a2 var":>11}{"theory":>11}{"ratio":>9}')
    for cell in cells:
        print(f'   {cell["n_side"]:>5}{cell["l0_m"]:>8.2f}'
              f'{cell["n_columns"]:>7}{cell["frame_var"]:>11.3f}'
              f'{cell["frame_var_theory"]:>11.3f}'
              f'{cell["frame_var"] / cell["frame_var_theory"]:>9.4f}'
              f'{cell["tilt_var"]:>11.4f}{cell["tilt_var_theory"]:>11.4f}'
              f'{cell["tilt_var"] / cell["tilt_var_theory"]:>9.4f}')

    # ---- 6 and 7. the cost and the guards ----
    print('')
    print('6. THE COST, AND 7. THE GUARDS OF PLAN SECTION 2')
    print('   E3: no lstsq fallback. E4: the minimum eigenvalue of BB^T must')
    print('   be positive. The condition number of cov_zz is SKIPPED above')
    print(f'   {COND_MAX_SIZE} x {COND_MAX_SIZE}, where the eigenvalue solve '
          f'is too slow.')
    print('')
    print(f'   {"N":>5}{"L0, m":>8}{"n_col":>7}{"cov_zz":>9}{"MB":>9}'
          f'{"setup, s":>10}{"add_row, ms":>12}{"min eig BBt":>13}'
          f'{"cond cov_zz":>13}{"lstsq":>8}')
    cost_rows = []
    guard_ok = True
    for cell in cells:
        guards = cell['guards']
        megabytes = guards['size'] ** 2 * 8.0 / 2 ** 20
        cond = 'skipped' if np.isnan(guards['cond']) \
            else f'{guards["cond"]:.3e}'
        print(f'   {cell["n_side"]:>5}{cell["l0_m"]:>8.2f}'
              f'{cell["n_columns"]:>7}{guards["size"]:>9}{megabytes:>9.1f}'
              f'{cell["setup_s"]:>10.2f}{cell["add_row_ms"]:>12.4f}'
              f'{guards["min_eig"]:>13.3e}{cond:>13}'
              f'{"YES" if guards["fallback"] else "no":>8}')
        guard_ok = guard_ok and (not guards['fallback']) \
            and guards['min_eig'] > 0.0
        cost_rows.append([
            cell['n_columns'], f'{cell["l0_m"]:.2f}', cell['n_side'],
            guards['size'], f'{megabytes:.3f}', f'{cell["setup_s"]:.4f}',
            f'{cell["add_row_ms"]:.6f}', int(guards['fallback']),
            f'{guards["min_eig"]:.6e}',
            '' if np.isnan(guards['cond']) else f'{guards["cond"]:.6e}',
            f'{cell["frame_var"] / cell["frame_var_theory"]:.6f}',
            f'{cell["tilt_var"] / cell["tilt_var_theory"]:.6f}'])
    print('')
    band_check('7. the guards: no lstsq fallback, BB^T positive', guard_ok)

    # ---- Phase 2 ----
    print('')
    print('=' * 78)
    print('PHASE 2. THE RECORD LENGTH AND THE PISTON')
    print('=' * 78)
    print(f'   one record of {WINDOW_RECORD_L0:.0f} L0 per seed, '
          f'{SEEDS_PHASE2} seeds, N = {N_SMALL}')
    print(f'   n_columns {N_COLUMNS_WINDOW}, windows '
          f'{WINDOW_LENGTHS_L0} outer scales')
    print('   The question: does the APPARENT rho(1 L0) fall to the theory as')
    print('   the analysis window grows? If it does, the FINDINGS excess is a')
    print('   record length, not a Markov memory.')

    window_cells = []
    for l0_m in (L0_TIGHT_M, L0_WIDE_M):
        for n_columns in N_COLUMNS_WINDOW:
            t0 = time.time()
            window_cells.append(
                run_window_cell(N_SMALL, l0_m, n_columns, seeds2))
            print(f'    L0 {l0_m:6.2f} m  n_columns {n_columns:3d}   '
                  f'{time.time() - t0:7.1f} s   '
                  f'(total {time.time() - t_start:7.1f} s)')

    print('')
    print(f'   {"L0, m":>8}{"n_col":>7}{"window, L0":>12}{"windows":>9}'
          f'{"rho(1 L0)":>12}{"SE":>9}{"theory":>9}{"excess":>9}')
    window_rows = []
    for cell in window_cells:
        for index, w in enumerate(cell['windows']):
            excess = cell['rho'][index] - cell['rho_theory']
            print(f'   {cell["l0_m"]:>8.2f}{cell["n_columns"]:>7}{w:>12.0f}'
                  f'{cell["counts"][index]:>9}{cell["rho"][index]:>12.4f}'
                  f'{cell["rho_se"][index]:>9.4f}{cell["rho_theory"]:>9.4f}'
                  f'{excess:>+9.4f}')
            window_rows.append([cell['n_columns'], f'{cell["l0_m"]:.2f}',
                                cell['n_side'], f'{w:.0f}',
                                cell['counts'][index],
                                f'{cell["rho"][index]:.6f}',
                                f'{cell["rho_se"][index]:.6f}',
                                f'{cell["rho_theory"]:.6f}'])
    print('')
    for cell in window_cells:
        short, long = float(cell['rho'][0]), float(cell['rho'][-1])
        band_check(f'Phase 2, L0 {cell["l0_m"]:.2f} m, '
                   f'n_columns {cell["n_columns"]}: the longest window reads '
                   f'the theory',
                   abs(long - cell['rho_theory']) < RHO_BAND,
                   f'{cell["rho"][0]:.4f} at {WINDOW_LENGTHS_L0[0]:.0f} L0 -> '
                   f'{long:.4f} at {WINDOW_LENGTHS_L0[-1]:.0f} L0, theory '
                   f'{cell["rho_theory"]:.4f}')
        print(f'          the apparent rho falls by {short - long:+.4f} from '
              f'the shortest window to the longest.')

    # ---- the outputs ----
    write_csv(ROWLAG_CSV,
              ['n_columns', 'L0_m', 'N', 'axis', 'k', 'r_m', 'rho_meas',
               'rho_se', 'rho_theory'], rowlag_rows)
    write_csv(PISTON_CSV,
              ['n_columns', 'L0_m', 'N', 'quantity', 'freq_cycles_per_row',
               'measured', 'theory', 'se'], piston_rows)
    write_csv(COST_CSV,
              ['n_columns', 'L0_m', 'N', 'cov_zz_size', 'cov_zz_MB',
               'setup_s', 'add_row_ms', 'lstsq_fallback', 'bbt_min_eig',
               'cond_cov_zz', 'frame_var_ratio', 'ztilt_ratio'], cost_rows)
    write_csv(WINDOW_CSV,
              ['n_columns', 'L0_m', 'N', 'window_L0', 'n_windows',
               'rho_1L0_meas', 'se', 'rho_1L0_theory'], window_rows)

    draw_rowlag(cells)
    draw_anisotropy(cells, aniso_table)
    draw_piston(cells)
    draw_cost(cells)
    draw_window(window_cells)

    print('')
    print('OUTPUTS')
    for path in (ROWLAG_CSV, PISTON_CSV, COST_CSV, WINDOW_CSV, ROWLAG_PNG,
                 ANISO_PNG, PISTON_PNG, COST_PNG, WINDOW_PNG):
        print(f'   {os.path.basename(path):<26}'
              f'{os.path.getsize(path) / 1024:8.1f} kB')

    print('')
    print('SUMMARY')
    failed = [label for label, ok in RESULTS if not ok]
    for label, ok in RESULTS:
        print(f'   [{"PASS" if ok else "FAIL"}] {label}')
    print(f'   {len(RESULTS) - len(failed)} of {len(RESULTS)} bands pass')
    if failed:
        print('   A failed band is a RESULT. Read the tables above.')

    print('')
    print(f'(elapsed {time.time() - t_start:.1f} s)')


if __name__ == '__main__':
    main()
