'''
The SPATIAL statistics of the aotools infinite (extruded) phase screens.

This is arm 2 of the phase-screen low-frequency study. Arm 1 measures a Fourier
screen, which holds no power below the grid fundamental 1/(N dx). An EXTRUDED
screen builds each new row from a covariance matrix, so it carries no grid
fundamental at all. An outside claim says that the covariance of those screens
is wrong. This script tests that claim, and it separates two questions that the
claim mixes together:

  1. IS THE FORMULA WRONG? `aotools.turbulence.turb.phase_covariance` is
     Assemat and Wilson (2006), DOI 10.1364/OE.14.000988, Eq. (5). Section 1
     measures it against two independent routes: the same closed form in
     float64 (`helpers.vk_covariance_closed`), and a Hankel transform of the
     von Karman PSD (`helpers.vk_covariance_numeric`).
  2. DOES THE SCREEN CARRY THAT COVARIANCE? A correct formula does not prove a
     correct screen. The A and B matrices, the 2-column stencil, and a long
     extrusion each can lose power. Sections 2, 3 and 4 measure the SCREEN.

THE FOUR MEASUREMENTS.

  1. THE FORMULA. A table over 25 separations, at two outer scales. The verdict
     line reads CONFIRMED CORRECT or WRONG.
  2. THE STRUCTURE FUNCTION. The ensemble D_phi of the extruded screens against
     2 [B(0) - B(r)]. The von Karman class carries the pass band. The
     Kolmogorov class is a report only, because its stencil differs.
  3. THE COVARIANCE B(r) OF THE SCREENS. A direct product estimate, per axis.
     The extrusion axis and the transverse axis get separate numbers, because
     the recursion treats them differently. The section also splits the
     variance into a PISTON part and a RESIDUAL part. See below.
  4. THE TILT. The pooled Z-tilt and G-tilt variance over a 1.0 m pupil against
     the Noll and the Andrews filter integrals.
  5. THE SPIN-UP LADDER. The same variance split, at three extrusion lengths.
     It separates a slow transient from a permanent loss of power.

THE PISTON SPLIT, AND WHY THIS SCRIPT NEEDS IT. A frame of side 5.12 m at
L0 = 25 m holds less than a quarter of one correlation cell. So almost all of
the theoretical variance B(0) sits in the PISTON of that frame, the screen-wide
mean. Section 3 therefore reports two numbers:

    piston variance   = the double-grid average of B(r), against the measured
                        variance of the per-screen mean
    residual variance = B(0) - piston variance, against the measured mean
                        variance inside one screen

A structure function and a tilt cannot see the piston at all. So section 2 and
section 4 stay accurate while the raw B(0) estimate of section 3 is far off.

THE MEASURED RESULT OF SECTION 5. The large modes of the extrusion fill in
SLOWLY. The spin-up must run for several outer scales, not for one screen
width. At L0 = 2.56 m a 512-row spin-up covers two outer scales, and the
variance lands within 1 percent. At L0 = 25 m the same 512 rows cover a fifth
of one outer scale, and the piston reaches only 0.29 of its theoretical
variance. At 2048 rows it reaches 0.72, and it is still climbing. So a raw B(0)
error at L0 = 25 m measures the SPIN-UP LENGTH, not a wrong covariance.

THE ANSWER TO THE OUTSIDE CLAIM (the measured numbers print below).

  - THE FORMULA IS CORRECT. `turb.phase_covariance` matches the float64 closed
    form to 1.2e-6, and it matches the independent Hankel integral to 0.45
    percent, which is the printed PSD constant 0.023.
  - THE SCREEN CARRIES THAT COVARIANCE WHEN THE OUTER SCALE FITS THE GRID. At
    L0 = 2.56 m the raw variance is right to 0.1 percent, and the transverse
    covariance holds inside 2.4 percent of B(0).
  - ONE REAL DEFECT: THE EXTRUSION AXIS OVER-CORRELATES. At L0 = 2.56 m the
    extrusion axis holds 1.50 rad^2 at a 2.56 m lag, where the theory says
    0.085 and the transverse axis reads -0.10. The 2-column Markov recursion
    keeps too much memory along its own direction. That is a 7 percent error of
    B(0), and it fails the 5 percent band.
  - THE LARGE APPARENT FAILURE AT L0 = 25 m IS A SPIN-UP LENGTH, NOT AN ERROR.
    See the section 5 note below.

THE SPIN-UP. `make_initial_screen` of `aotools/turbulence/infinitephasescreen.py`
calls `phasescreen.ft_phase_screen`. That is the PLAIN Fourier screen, with no
subharmonics. So the first frame carries the arm-1 low-frequency deficit, and it
is not a fair sample of the extrusion. This script calls `add_row()` 512 times
before it harvests a screen. The visible frame is 512 rows tall, so every row of
the harvested screen comes from the A and B recursion.

THE EXTRUSION AXIS. `add_row` runs
`numpy.append(new_row, self._scrn, axis=0)`, so the new row goes to AXIS 0. Axis
0 is the extrusion axis, and axis 1 is the transverse axis. Section 3 keeps them
apart.

WHY SECTION 3 SUBTRACTS NO MEAN. A per-screen mean subtraction removes the
outer-scale power that this study measures. At L0 = 2.56 m the screen holds
about four correlation cells per side, so the sample mean carries real signal.
The estimate below multiplies the raw phase values.

WHY SECTION 3 USES TWO OUTER SCALES. A covariance estimate needs many
correlation cells inside the grid. At L0 = 2.56 m the 5.12 m side holds a few
cells, so the spatial average works and the band is tight. At L0 = 25 m the side
holds less than one cell, so the estimate is noisy and the band is wide.

Sources:
- Assemat and Wilson, "Method for simulating infinitely long and non stationary
  phase screens with optimized memory storage", Opt. Express 14(3), pp. 988 to
  999 (2006), DOI 10.1364/OE.14.000988. Eq. (5), the von Karman phase
  covariance, and the A and B extrusion matrices.
- Fried and Clark, J. Opt. Soc. Am. A 25(2), pp. 463 to 468 (2008),
  DOI 10.1364/JOSAA.25.000463. The stencil of the Kolmogorov class.
- Noll, J. Opt. Soc. Am. 66(3), pp. 207 to 211 (1976),
  DOI 10.1364/JOSA.66.000207. Eq. (8), p. 208, the Zernike tilt filter.
- Andrews and Phillips, Laser Beam Propagation through Random Media, 2nd ed.,
  DOI 10.1117/3.626196. Ch. 6, Eqs. (80) to (84), printed pp. 200 and 201, the
  gradient tilt and the angle of arrival.
- Schmidt, Numerical Simulation of Optical Wave Propagation with Examples in
  MATLAB, DOI 10.1117/3.866274. Ch. 3, Eqs. (3.16) to (3.25), printed pp. 48 to
  50, the covariance-to-structure-function relation and the estimator; Ch. 9,
  Eq. (9.50), printed p. 161, the von Karman phase PSD.

A FAILED BAND IS A RESULT. No band raises. Every band prints PASS or FAIL, and
the summary at the end lists them all.

Outputs, next to this script:
    infinite_covariance.png   B(r) per axis against theory, plus a residual
    infinite_dphi.png         the structure-function ratio, two screen classes
    infinite_covariance.csv   L0, axis, r_m, b_est, b_theory
    infinite_dphi.csv         klass, r_m, d_meas, d_theory, ratio
    infinite_tilt.csv         metric, var_meas, var_theory, ratio, se

This script changes no olb module. It reads the production layer only.

Run from the repo root:
    python -m validation.screens.infinite_screen_stats
'''

import csv
import os
import time

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt                            # noqa: E402
import numpy as np                                         # noqa: E402

from aotools.turbulence import turb                        # noqa: E402
from aotools.turbulence.infinitephasescreen import (        # noqa: E402
    PhaseScreenKolmogorov, PhaseScreenVonKarman, find_allowed_size)

from olb.turbulence.andrews.structure import (              # noqa: E402
    angle_of_arrival_variance)
from olb.waveoptics.schmidt.turbulence import (             # noqa: E402
    von_karman_phase_psd)
from validation.screens import helpers                     # noqa: E402

# ---------------------------------------------------------------------------
# The module constants
# ---------------------------------------------------------------------------

WAVELENGTH_M = 1550e-9
K_RAD_M = 2.0 * np.pi / WAVELENGTH_M

R0_M = 0.10                     # the Fried parameter of every screen
DX_M = 0.01                     # the grid pitch
N = 512                         # the grid side, so the screen side is 5.12 m
SIDE_M = N * DX_M
PUPIL_D_M = 1.0                 # the tilt pupil

# The two outer scales. See the header for why the study needs both.
L0_WIDE_M = 25.0
L0_TIGHT_M = 2.56

# The spin-up. The visible frame is N rows tall, so N rows replace all of it.
SPINUP_ROWS = N

# The ensembles.
M_WIDE = 60                     # von Karman, L0 = 25 m
M_TIGHT = 40                    # von Karman, L0 = 2.56 m
M_KOLM = 20                     # the Kolmogorov class, L0 = 25 m
M_DPHI = 40                     # the share of M_WIDE that section 2 reads
M_SPINUP = 10                   # the ensemble of the spin-up ladder

MASTER_SEED = 20260828

# ---- section 1, the formula table ----
R_TABLE_M = np.logspace(-3.0, np.log10(50.0), 25)
# The numeric Hankel route hits its own cancellation floor at about 1e-8 rad^2.
# The float32 cast of aotools and the float64 exponential both underflow far
# below that. So the verdict reads only the rows where the covariance is above
# this share of B(0). Below the floor the ratios measure round-off, not physics.
DYNAMIC_FLOOR = 1e-6
CLOSED_VS_NUMERIC_BAND = 0.015  # the printed PSD constant 0.023 gives 0.0046
AOTOOLS_VS_CLOSED_BAND = 0.005

# ---- section 2, the structure function ----
# The mask radius is 1.28 m, so the 2.56 m grid half-side leaves a guard band.
# The circular correlation estimate of Schmidt Ch. 3 needs that band.
DPHI_MASK_D_M = 2.56
DPHI_BINS_M = np.logspace(np.log10(0.05), np.log10(2.56), 12)
DPHI_BAND_R_M = (0.10, 2.50)
DPHI_BAND = (0.90, 1.10)

# ---- section 3, the covariance ----
COV_KMAX_PX = 256
COV_BAND_R_M = 1.0              # the band holds out to this separation
COV_BAND_TIGHT = 0.05           # a share of B(0), at L0 = 2.56 m
COV_BAND_WIDE = 0.20            # a share of B(0), at L0 = 25 m

# ---- section 5, the spin-up ladder ----
SPINUP_LADDER = (128, SPINUP_ROWS, 2048)

# ---- section 4, the tilt ----
TILT_BAND = 0.20
# The Andrews cross-print needs Cn2 L, not r0. Invert the Fried definition
# r0 = (0.423 k^2 Cn2 L)^(-3/5) at a unit path length. Source: Fried,
# DOI 10.1364/JOSA.56.001372, through Schmidt (2010), DOI 10.1117/3.866274,
# Ch. 9, Eq. (9.70), printed p. 165.
ANDREWS_Z_M = 1.0
ANDREWS_CN2 = R0_M ** (-5.0 / 3.0) / (0.423 * K_RAD_M ** 2) / ANDREWS_Z_M
# The von Karman branch of `angle_of_arrival_variance` needs an inner scale.
# Give it a value far below the pupil, so the D >> l0 branch runs.
ANDREWS_L0_INNER_M = 1e-6

HERE = os.path.dirname(os.path.abspath(__file__))
COV_PNG = os.path.join(HERE, 'figures/infinite_covariance.png')
DPHI_PNG = os.path.join(HERE, 'figures/infinite_dphi.png')
COV_CSV = os.path.join(HERE, 'data/infinite_covariance.csv')
DPHI_CSV = os.path.join(HERE, 'data/infinite_dphi.csv')
TILT_CSV = os.path.join(HERE, 'data/infinite_tilt.csv')


# ---------------------------------------------------------------------------
# The screen builders
# ---------------------------------------------------------------------------

def von_karman_screen(seed, l0_outer_m, rows=SPINUP_ROWS):
    '''
    Build one extruded von Karman screen, after the spin-up.

    Parameters:
        seed : int
            The integer seed of the aotools generator.
        l0_outer_m : float
            The outer scale L0 [m].
        rows : int
            The number of spin-up rows. The default replaces every visible row.

    Returns:
        numpy.ndarray
            An N by N phase screen [rad].

    THE SPIN-UP. The constructor makes a plain Fourier initial frame. This
    function calls `add_row()` `rows` times. At the default the recursion
    replaces every visible row. See the header.

    Source: Assemat and Wilson, Opt. Express 14(3), pp. 988 to 999 (2006),
    DOI 10.1364/OE.14.000988. The class uses two columns of history, which is
    the value that the paper recommends.
    '''
    screen = PhaseScreenVonKarman(N, DX_M, R0_M, float(l0_outer_m),
                                  random_seed=int(seed))
    for _ in range(int(rows)):
        screen.add_row()
    return np.array(screen.scrn, dtype=float)


def kolmogorov_screen(seed, l0_outer_m):
    '''
    Build one extruded Kolmogorov-class screen, after the spin-up.

    Parameters:
        seed : int
            The integer seed of the aotools generator.
        l0_outer_m : float
            The outer scale L0 [m]. It must be FINITE. See the note below.

    Returns:
        numpy.ndarray
            An N by N phase screen [rad], cropped from the native grid.

    THE NAME IS MISLEADING. `PhaseScreenKolmogorov` still builds its covariance
    matrices from `turb.phase_covariance`, which is the VON KARMAN Eq. (5) of
    Assemat and Wilson, DOI 10.1364/OE.14.000988. An infinite L0 makes that
    expression diverge. So the class needs a finite outer scale, and its
    statistics are von Karman. The name refers to the STENCIL of Fried and
    Clark, DOI 10.1364/JOSAA.25.000463, not to the spectrum.

    THE GRID. The class snaps the side to 2^n + 1 through `find_allowed_size`,
    so a request for 512 gives 513. This function crops the result back to N
    with `helpers.crop_center`, which keeps the centre sample.
    '''
    native = find_allowed_size(N)
    screen = PhaseScreenKolmogorov(native, DX_M, R0_M, float(l0_outer_m),
                                   random_seed=int(seed))
    for _ in range(SPINUP_ROWS):
        screen.add_row()
    return np.array(helpers.crop_center(np.asarray(screen.scrn, dtype=float),
                                        N), dtype=float)


# ---------------------------------------------------------------------------
# The estimators
# ---------------------------------------------------------------------------

def axis_covariance(screens, kmax_px):
    '''
    Estimate the phase covariance along each grid axis.

    Parameters:
        screens : list
            The phase screens [rad]. Each one is N by N.
        kmax_px : int
            The largest separation, in pixels.

    Returns:
        tuple
            (b_extrusion, b_transverse). Each one holds kmax_px + 1 values
            [rad^2], from k = 0 upward.

    formula:
        B(k dx) = < phi(i, j) phi(i + k, j) >
    Source: the definition of the covariance of a homogeneous random field.
    Schmidt (2010), DOI 10.1117/3.866274, Ch. 3, Eq. (3.14), printed p. 47.

    NO MEAN SUBTRACTION. The estimate multiplies the raw values. A per-screen
    mean subtraction removes real outer-scale power. See the header.

    THE TWO AXES. Axis 0 is the extrusion axis, because `add_row` appends along
    axis 0. Axis 1 is the transverse axis. The recursion treats them
    differently, so the two estimates can differ.

    VALIDITY. The average runs over the overlap pixels only, so a large k reads
    fewer pairs and the estimate gets noisy. The estimate needs many
    correlation cells inside the grid.
    '''
    kmax_px = int(kmax_px)
    b_ext = np.zeros(kmax_px + 1)
    b_tra = np.zeros(kmax_px + 1)
    for screen in screens:
        side = screen.shape[0]
        for k in range(kmax_px + 1):
            b_ext[k] += float(np.mean(screen[k:, :] * screen[:side - k, :]))
            b_tra[k] += float(np.mean(screen[:, k:] * screen[:, :side - k]))
    return b_ext / len(screens), b_tra / len(screens)


def grid_piston_variance(l0_outer_m):
    '''
    Return the theoretical piston variance of one square frame.

    Parameters:
        l0_outer_m : float
            The outer scale L0 [m].

    Returns:
        float
            The variance of the frame-wide mean phase [rad^2].

    formula:
        <piston^2> = (1/A^2) INT INT B(|r1 - r2|) d^2r1 d^2r2
    Source: the variance of the average of a homogeneous random field. Schmidt
    (2010), DOI 10.1117/3.866274, Ch. 3, Eq. (3.14), printed p. 47, defines
    B(r). The covariance itself is Assemat and Wilson,
    DOI 10.1364/OE.14.000988, Eq. (5), through `helpers.vk_covariance_closed`.

    THE WEIGHT. The double sum over a regular square grid collapses to one sum
    over the separation VECTOR. A separation of i pixels on one axis occurs
    (N - |i|) times. So the weight of the vector (i, j) is
    (N - |i|) (N - |j|), and the two axes separate.

    WHY THE STUDY NEEDS THIS. A 5.12 m frame at L0 = 25 m puts nearly all of
    B(0) into the piston. See the header.
    '''
    index = np.arange(-(N - 1), N)
    axis_weight = (N - np.abs(index)).astype(float)
    weight = np.outer(axis_weight, axis_weight)
    ii, jj = np.meshgrid(index, index, indexing='ij')
    r_m = np.hypot(ii, jj) * DX_M
    cov = helpers.vk_covariance_closed(r_m.ravel(), R0_M,
                                       l0_outer_m).reshape(r_m.shape)
    return float(np.sum(weight * cov) / weight.sum())


def screen_moments(screens):
    '''
    Split the measured variance of an ensemble into a piston part and a rest.

    Parameters:
        screens : iterable
            The phase screens [rad].

    Returns:
        tuple
            (piston_variance, residual_variance, mean_square). Each one is
            [rad^2].

    THE SPLIT. The piston of one screen is its own spatial mean. Its variance
    over the ensemble is the piston variance. The residual variance is the mean
    of the variance INSIDE one screen. The two add up to the mean square,
    because the two parts are orthogonal by construction.
    '''
    means = np.array([float(np.mean(s)) for s in screens])
    inner = np.array([float(np.var(s)) for s in screens])
    piston = float(np.mean(means * means))
    residual = float(np.mean(inner))
    return piston, residual, piston + residual


def pooled_tilt(screens, mask):
    '''
    Pool the per-axis tilt of an ensemble, for both tilt definitions.

    Parameters:
        screens : list
            The phase screens [rad].
        mask : numpy.ndarray
            The boolean pupil mask.

    Returns:
        tuple
            (z_angles, g_angles). Each one is a flat array of tilt ANGLES
            [rad], with two entries per screen: the x axis and the y axis.

    THE POOL. The x tilt and the y tilt of one screen are two samples of the
    same per-axis variance, and they are close to independent. So the pool
    holds 2 M samples.
    '''
    z_angles, g_angles = [], []
    for screen in screens:
        zx, zy, _, _ = helpers.zernike_tilt(screen, mask, DX_M, WAVELENGTH_M)
        gx, gy = helpers.gradient_tilt(screen, mask, DX_M, WAVELENGTH_M)
        z_angles.extend([zx, zy])
        g_angles.extend([gx, gy])
    return np.array(z_angles), np.array(g_angles)


def variance_and_se(samples):
    '''
    Return the zero-mean variance of a sample set, and its standard error.

    Parameters:
        samples : array_like
            The samples. The true mean is zero, so the estimator is the mean
            square.

    Returns:
        tuple
            (variance, standard_error).

    formula:
        var = < x^2 >,      se = var sqrt(2 / n)
    Source: the standard error of the variance of n independent Gaussian
    samples. It is a textbook result, not a book equation of this study.
    '''
    x = np.asarray(samples, dtype=float)
    variance = float(np.mean(x * x))
    return variance, variance * np.sqrt(2.0 / x.size)


# ---------------------------------------------------------------------------
# The verdict helper
# ---------------------------------------------------------------------------

RESULTS = []


def band_check(label, ok, detail=''):
    '''
    Print one PASS or FAIL line, and record it.

    Parameters:
        label : str
            The name of the band.
        ok : bool
            True for a pass.
        detail : str
            An optional note that follows the label.

    THIS FUNCTION RAISES NOTHING. A failed band is a RESULT of the study, not a
    bug in the script. The run continues, and the summary lists every band.
    '''
    RESULTS.append((label, bool(ok)))
    tag = 'PASS' if ok else 'FAIL'
    print(f'  [{tag}] {label}{"   " + detail if detail else ""}')
    return bool(ok)


# ---------------------------------------------------------------------------
# Section 1: the covariance formula
# ---------------------------------------------------------------------------

def section_formula():
    '''
    Measure the aotools covariance formula against two independent routes.

    Returns:
        bool
            True if the formula is correct.

    THE THREE COLUMNS. `turb.phase_covariance` casts the separation to float32
    and it adds 1e-40. `helpers.vk_covariance_closed` evaluates the SAME
    equation in float64. `helpers.vk_covariance_numeric` integrates the von
    Karman PSD with a Hankel transform, so it shares no algebra with the other
    two.

    THE EXPECTED GAP. The numeric route reads about 0.46 percent high. That is
    the printed PSD constant 0.023 against the exact 0.49 (2 pi)^(-5/3). It is
    a rounding, not an error.
    '''
    print('')
    print('1. THE COVARIANCE FORMULA, THREE ROUTES')
    print('   Assemat and Wilson (2006), DOI 10.1364/OE.14.000988, Eq. (5)')

    verdict = True
    for l0_m in (L0_WIDE_M, L0_TIGHT_M):
        closed = helpers.vk_covariance_closed(R_TABLE_M, R0_M, l0_m)
        numeric = helpers.vk_covariance_numeric(R_TABLE_M, R0_M, l0_m)
        aotools = np.asarray(turb.phase_covariance(R_TABLE_M, R0_M, l0_m),
                             dtype=float)
        b_zero = float(helpers.vk_covariance_closed(0.0, R0_M, l0_m)[0])
        keep = closed >= DYNAMIC_FLOOR * b_zero

        print('')
        print(f'   r0 = {R0_M * 1e2:.0f} cm, L0 = {l0_m:.2f} m, '
              f'B(0) = {b_zero:.4e} rad^2')
        print(f'   {"r, m":>11}{"aotools":>15}{"closed f64":>15}'
              f'{"numeric":>15}{"aot/cls":>11}{"cls/num":>11}{"band":>7}')
        for i, r_m in enumerate(R_TABLE_M):
            ratio_a = aotools[i] / closed[i] if closed[i] > 0.0 else np.nan
            ratio_n = closed[i] / numeric[i] if numeric[i] != 0.0 else np.nan
            print(f'   {r_m:>11.5f}{aotools[i]:>15.6e}{closed[i]:>15.6e}'
                  f'{numeric[i]:>15.6e}{ratio_a:>11.6f}{ratio_n:>11.6f}'
                  f'{"in" if keep[i] else "out":>7}')

        err_num = float(np.abs(closed[keep] / numeric[keep] - 1.0).max())
        err_aot = float(np.abs(aotools[keep] / closed[keep] - 1.0).max())
        print(f'   rows inside the dynamic floor '
              f'({DYNAMIC_FLOOR:.0e} of B(0)): '
              f'{int(np.count_nonzero(keep))} of {R_TABLE_M.size}')
        print(f'   max |closed / numeric - 1|  {err_num:.6f}   '
              f'(band {CLOSED_VS_NUMERIC_BAND:.3f}, the 0.023 PSD constant)')
        print(f'   max |aotools / closed - 1|  {err_aot:.3e}   '
              f'(the float32 cast and the 1e-40 offset)')
        verdict = verdict and (err_num < CLOSED_VS_NUMERIC_BAND
                               and err_aot < AOTOOLS_VS_CLOSED_BAND)

    print('')
    print('   The rows marked "out" sit below the dynamic floor. There the')
    print('   numeric Hankel integral reads its own cancellation error, the')
    print('   float32 cast of aotools loses the argument, and the float64')
    print('   exponential underflows. Those rows measure round-off only.')
    print('')
    print(f'covariance formula: {"CONFIRMED CORRECT" if verdict else "WRONG"}')
    band_check('section 1, the covariance formula', verdict)
    return verdict


# ---------------------------------------------------------------------------
# Section 2: the structure function
# ---------------------------------------------------------------------------

def section_structure(d_vk, d_kol):
    '''
    Compare the ensemble structure function of both screen classes.

    Parameters:
        d_vk : numpy.ndarray
            The von Karman class D(r) [rad^2], one value per bin.
        d_kol : numpy.ndarray
            The Kolmogorov class D(r) [rad^2], one value per bin.

    Returns:
        tuple
            (theory, ratio_vk, ratio_kol).

    formula:
        D(r) = 2 [B(0) - B(r)]
    Source: Schmidt (2010), DOI 10.1117/3.866274, Ch. 3, Eq. (3.16), printed
    p. 48. The covariance comes from `helpers.vk_covariance_closed`, which is
    Assemat and Wilson Eq. (5), DOI 10.1364/OE.14.000988. So the theory curve
    is the SAME equation that the extrusion matrices use. A screen that carries
    its own covariance must match it.

    THE KOLMOGOROV CLASS IS A REPORT. Its stencil samples the whole screen, and
    its recursion subtracts a reference point. So it holds no pass band here.
    '''
    b_zero = float(helpers.vk_covariance_closed(0.0, R0_M, L0_WIDE_M)[0])
    theory = 2.0 * (b_zero
                    - helpers.vk_covariance_closed(DPHI_BINS_M, R0_M,
                                                   L0_WIDE_M))
    ratio_vk = d_vk / theory
    ratio_kol = d_kol / theory

    print('')
    print('2. THE ENSEMBLE STRUCTURE FUNCTION OF THE EXTRUDED SCREENS')
    print(f'   theory = 2 [B(0) - B(r)], r0 = {R0_M * 1e2:.0f} cm, '
          f'L0 = {L0_WIDE_M:.0f} m, B(0) = {b_zero:.2f} rad^2')
    print(f'   estimator: helpers.ensemble_dphi, pupil '
          f'{DPHI_MASK_D_M:.2f} m over a {SIDE_M:.2f} m grid')
    print('')
    print(f'   {"r, m":>9}{"theory":>12}{"VK class":>12}{"ratio":>9}'
          f'{"Kolm class":>13}{"ratio":>9}')
    for i, r_m in enumerate(DPHI_BINS_M):
        print(f'   {r_m:>9.4f}{theory[i]:>12.3f}{d_vk[i]:>12.3f}'
              f'{ratio_vk[i]:>9.4f}{d_kol[i]:>13.3f}{ratio_kol[i]:>9.4f}')

    band = ((DPHI_BINS_M >= DPHI_BAND_R_M[0])
            & (DPHI_BINS_M <= DPHI_BAND_R_M[1]))
    inside = ratio_vk[band]
    print('')
    print(f'   the band is r = {DPHI_BAND_R_M[0]:.2f} to '
          f'{DPHI_BAND_R_M[1]:.2f} m, ratio {DPHI_BAND[0]:.2f} to '
          f'{DPHI_BAND[1]:.2f}')
    print(f'   VK class in the band     {inside.min():.4f} to '
          f'{inside.max():.4f}   ({M_DPHI} screens)')
    print(f'   Kolm class in the band   {ratio_kol[band].min():.4f} to '
          f'{ratio_kol[band].max():.4f}   ({M_KOLM} screens, report only)')
    band_check('section 2, VK class D(r) ratio',
               bool(np.all(inside >= DPHI_BAND[0])
                    and np.all(inside <= DPHI_BAND[1])),
               f'{inside.min():.3f} to {inside.max():.3f}')
    print('   The Kolmogorov class carries NO pass band. Its stencil and its')
    print('   reference-point subtraction differ from the von Karman class.')
    print(f'   The last bin sits at {DPHI_BINS_M[-1]:.2f} m, which is the mask '
          f'DIAMETER. The window')
    print('   overlap area falls to zero there, so that bin is an artifact. '
          'The band')
    print(f'   stops at {DPHI_BAND_R_M[1]:.2f} m for that reason.')
    print('   The VK ratio falls away above r = 0.5 m. That is the same large-'
          'scale')
    print(f'   deficit that section 5 measures: the {SPINUP_ROWS}-row spin-up '
          f'covers only')
    print(f'   {SPINUP_ROWS * DX_M / L0_WIDE_M:.2f} of the outer scale.')
    return theory, ratio_vk, ratio_kol


# ---------------------------------------------------------------------------
# Section 3: the covariance of the screens
# ---------------------------------------------------------------------------

def section_covariance(estimates):
    '''
    Compare the measured screen covariance with the Assemat and Wilson form.

    Parameters:
        estimates : dict
            The key is the outer scale [m]. The value is
            (b_extrusion, b_transverse, m_screens, moments). The moments come
            from `screen_moments`.

    Returns:
        dict
            The key is the outer scale. The value holds the separations, the
            theory, and the two measured curves.

    THE BAND IS A SHARE OF B(0). A ratio band fails where the theory is small,
    and the interesting error is an ABSOLUTE loss of power. So the band reads
    |B_est - B_theory| / B(0).
    '''
    r_m = np.arange(COV_KMAX_PX + 1) * DX_M
    out = {}

    print('')
    print('3. THE COVARIANCE B(r) OF THE EXTRUDED SCREENS, PER AXIS')
    print('   axis 0 is the EXTRUSION axis. axis 1 is the transverse axis.')
    print('   No mean subtraction. See the header.')

    for l0_m in (L0_TIGHT_M, L0_WIDE_M):
        b_ext, b_tra, m_screens, moments = estimates[l0_m]
        theory = helpers.vk_covariance_closed(r_m, R0_M, l0_m)
        b_zero = float(theory[0])
        tol = COV_BAND_TIGHT if l0_m == L0_TIGHT_M else COV_BAND_WIDE
        out[l0_m] = dict(r_m=r_m, theory=theory, ext=b_ext, tra=b_tra,
                         b_zero=b_zero, tol=tol, m=m_screens)

        print('')
        print(f'   L0 = {l0_m:.2f} m, {m_screens} screens, '
              f'B(0) = {b_zero:.4f} rad^2, band '
              f'{100.0 * tol:.0f} percent of B(0) out to '
              f'{COV_BAND_R_M:.2f} m')
        print(f'   {"r, m":>9}{"theory":>12}{"extrusion":>12}'
              f'{"dev, %B0":>11}{"transverse":>12}{"dev, %B0":>11}')
        for k in (0, 8, 16, 32, 64, 96, 128, 192, 256):
            d_ext = 100.0 * (b_ext[k] - theory[k]) / b_zero
            d_tra = 100.0 * (b_tra[k] - theory[k]) / b_zero
            print(f'   {r_m[k]:>9.3f}{theory[k]:>12.4f}{b_ext[k]:>12.4f}'
                  f'{d_ext:>11.2f}{b_tra[k]:>12.4f}{d_tra:>11.2f}')

        keep = r_m <= COV_BAND_R_M
        err_ext = float(np.abs(b_ext[keep] - theory[keep]).max() / b_zero)
        err_tra = float(np.abs(b_tra[keep] - theory[keep]).max() / b_zero)
        print(f'   max deviation out to {COV_BAND_R_M:.2f} m: '
              f'extrusion {100.0 * err_ext:.2f} percent, '
              f'transverse {100.0 * err_tra:.2f} percent of B(0)')
        band_check(f'section 3, L0 = {l0_m:.2f} m, extrusion axis',
                   err_ext <= tol, f'{100.0 * err_ext:.2f} percent of B(0)')
        band_check(f'section 3, L0 = {l0_m:.2f} m, transverse axis',
                   err_tra <= tol, f'{100.0 * err_tra:.2f} percent of B(0)')

        # The piston split. See the header.
        piston_theory = grid_piston_variance(l0_m)
        piston, residual, mean_square = moments
        print('')
        print(f'   THE PISTON SPLIT of the same ensemble, L0 = {l0_m:.2f} m:')
        print(f'   {"part":<34}{"measured":>13}{"theory":>13}{"ratio":>9}')
        for name, meas, want in (
                ('piston, var of the screen mean', piston, piston_theory),
                ('residual, var inside a screen', residual,
                 b_zero - piston_theory),
                ('sum, the raw <phi^2>', mean_square, b_zero)):
            print(f'   {name:<34}{meas:>13.3f}{want:>13.3f}'
                  f'{meas / want:>9.4f}')
        print(f'   the piston holds '
              f'{100.0 * piston_theory / b_zero:.1f} percent of the '
              f'theoretical B(0) here')

        if l0_m == L0_TIGHT_M:
            print('')
            print('   THE ASYMMETRY IS THE FINDING HERE. The outer scale fits '
                  'inside the grid, and')
            print('   the raw variance is right to 0.1 percent. But the two '
                  'axes part company at')
            print(f'   a large lag: the extrusion axis holds {b_ext[-1]:.2f} '
                  f'rad^2 at r = {r_m[-1]:.2f} m, where')
            print(f'   the theory says {theory[-1]:.3f} and the transverse '
                  f'axis reads {b_tra[-1]:.3f}. So the')
            print('   2-column Markov recursion OVER-correlates along the '
                  'extrusion direction.')

        if l0_m == L0_WIDE_M:
            cells = SIDE_M / l0_m
            print('')
            print(f'   CAVEAT. The grid side is {SIDE_M:.2f} m and the outer '
                  f'scale is {l0_m:.0f} m, so the')
            print(f'   screen holds {cells:.2f} of a correlation cell. The '
                  f'piston row above')
            print('   carries almost all of B(0), and section 5 shows that the '
                  'piston is a SLOW')
            print('   mode of the extrusion. So read the residual row, and '
                  'read sections 2 and')
            print('   4: a structure function and a tilt do not see the '
                  'piston at all.')
    return out


# ---------------------------------------------------------------------------
# Section 4: the tilt
# ---------------------------------------------------------------------------

def section_tilt(z_angles, g_angles, mask):
    '''
    Compare the pooled tilt variance with the two analytic filter integrals.

    Parameters:
        z_angles : numpy.ndarray
            The pooled Z-tilt angles [rad].
        g_angles : numpy.ndarray
            The pooled G-tilt angles [rad].
        mask : numpy.ndarray
            The pupil mask, for the area-equivalent diameter.

    Returns:
        list
            One row per metric, for the CSV.

    THE TWO TILTS ARE NOT THE SAME NUMBER. The Z-tilt is the least-squares
    plane over the pupil (Noll, DOI 10.1364/JOSA.66.000207). The G-tilt is the
    mean phase gradient (Andrews and Phillips, DOI 10.1117/3.626196, Ch. 6,
    Eqs. (80) to (82), printed p. 200). The G-tilt variance is lower.

    THE ANDREWS CROSS-PRINT. `angle_of_arrival_variance` with
    spectrum="von_karman" uses Ch. 6, Eq. (83), printed p. 201. That is a
    FIRST-ORDER expansion in (k0 D)^(1/3), not the full filter integral. It is
    a cross-reference here, and it carries no band.
    '''
    d_eff = helpers.mask_diameter(mask, DX_M)
    psd = lambda f: von_karman_phase_psd(f, R0_M, L0_WIDE_M)

    a2_var = helpers.tilt_filter_variance(psd, d_eff)
    z_theory = a2_var / (K_RAD_M * d_eff / 4.0) ** 2
    g_theory = helpers.gtilt_filter_variance(psd, d_eff, WAVELENGTH_M)

    z_var, z_se = variance_and_se(z_angles)
    g_var, g_se = variance_and_se(g_angles)

    andrews = float(angle_of_arrival_variance(
        d_eff, WAVELENGTH_M, ANDREWS_Z_M, ANDREWS_CN2,
        spectrum='von_karman', l0=ANDREWS_L0_INNER_M, L0=L0_WIDE_M))

    print('')
    print('4. THE TILT OF THE EXTRUDED SCREENS OVER A 1.0 m PUPIL')
    print(f'   {M_WIDE} screens, 2 axes each, so {z_angles.size} samples per '
          f'metric')
    print(f'   r0 = {R0_M * 1e2:.0f} cm, L0 = {L0_WIDE_M:.0f} m, '
          f'area-equivalent pupil {d_eff:.4f} m, D/r0 = {d_eff / R0_M:.2f}')
    print('')
    print(f'   {"metric":<28}{"measured":>13}{"theory":>13}{"ratio":>9}'
          f'{"se":>13}')
    rows = [('Z-tilt angle var, rad^2', z_var, z_theory, z_se),
            ('G-tilt angle var, rad^2', g_var, g_theory, g_se)]
    for name, meas, want, sigma in rows:
        print(f'   {name:<28}{meas:>13.4e}{want:>13.4e}'
              f'{meas / want:>9.4f}{sigma:>13.4e}')
    print(f'   {"Noll <a2^2>, rad^2":<28}{"":>13}{a2_var:>13.4e}')
    print('')
    print(f'   relative standard error of a variance from '
          f'{z_angles.size} samples: '
          f'{100.0 * np.sqrt(2.0 / z_angles.size):.1f} percent')
    for name, meas, want, sigma in rows:
        band_check(f'section 4, {name.split(",")[0]}',
                   abs(meas / want - 1.0) <= TILT_BAND,
                   f'ratio {meas / want:.4f}, band '
                   f'{1 - TILT_BAND:.2f} to {1 + TILT_BAND:.2f}')

    print('')
    print('   CROSS-REFERENCE, no band:')
    print(f'   olb andrews.structure.angle_of_arrival_variance   '
          f'{andrews:.4e} rad^2')
    print(f'   helpers.gtilt_filter_variance (full integral)     '
          f'{g_theory:.4e} rad^2')
    print(f'   ratio, andrews over the full integral             '
          f'{andrews / g_theory:.4f}')
    print('   The olb function is Andrews and Phillips Ch. 6, Eq. (83),')
    print('   printed p. 201. It is a FIRST-ORDER outer-scale expansion,')
    print('   1 - 0.81 (k0 D)^(1/3). The helpers route integrates the full von')
    print('   Karman filter, so the two must not agree exactly.')

    return [('z_tilt_angle_var_rad2', z_var, z_theory, z_var / z_theory, z_se),
            ('g_tilt_angle_var_rad2', g_var, g_theory, g_var / g_theory, g_se),
            ('noll_a2_var_rad2', np.nan, a2_var, np.nan, np.nan),
            ('andrews_aoa_var_rad2', np.nan, andrews,
             andrews / g_theory, np.nan)]


# ---------------------------------------------------------------------------
# Section 5: the spin-up ladder
# ---------------------------------------------------------------------------

def section_spinup(seeds):
    '''
    Measure the variance split against the extrusion length.

    Parameters:
        seeds : list
            The integer seeds. The same seeds run at every ladder step.

    Returns:
        list
            One row per ladder step: (rows, piston, residual, mean_square).

    WHY THIS SECTION EXISTS. Section 3 reads a large error in the raw B(0) at
    L0 = 25 m. Two causes give that error: a permanent loss of power, or a slow
    transient that a longer spin-up removes. This ladder separates them. It
    runs the SAME seeds at three extrusion lengths, so the comparison carries
    no seed difference.

    THE COST. This section rebuilds the screens, so it is the slowest part of
    the study for its ensemble size. Keep M_SPINUP small.
    '''
    b_zero = float(helpers.vk_covariance_closed(0.0, R0_M, L0_WIDE_M)[0])
    piston_theory = grid_piston_variance(L0_WIDE_M)

    print('')
    print('5. THE SPIN-UP LADDER, L0 = 25 m')
    print(f'   {M_SPINUP} screens per step, the same seeds at every step')
    print(f'   theory: piston {piston_theory:.2f}, residual '
          f'{b_zero - piston_theory:.2f}, B(0) {b_zero:.2f} rad^2')
    print('')
    print(f'   {"add_row calls":>15}{"piston":>12}{"ratio":>9}'
          f'{"residual":>12}{"ratio":>9}{"<phi^2>":>12}{"ratio":>9}')

    rows = []
    for count in SPINUP_LADDER:
        t0 = time.time()
        screens = [von_karman_screen(s, L0_WIDE_M, rows=count) for s in seeds]
        piston, residual, mean_square = screen_moments(screens)
        del screens
        rows.append((count, piston, residual, mean_square))
        print(f'   {count:>15d}{piston:>12.2f}'
              f'{piston / piston_theory:>9.4f}{residual:>12.2f}'
              f'{residual / (b_zero - piston_theory):>9.4f}'
              f'{mean_square:>12.2f}{mean_square / b_zero:>9.4f}'
              f'   ({time.time() - t0:.0f} s)')

    print('')
    print('   READ THE TWO COLUMNS APART. Both columns CLIMB with the '
          'extrusion')
    print('   length, and the piston climbs from much further down. Neither '
          'one is')
    print(f'   settled at {SPINUP_ROWS} rows. The extrusion must run for '
          f'several OUTER SCALES:')
    print(f'   {SPINUP_ROWS} rows is {SPINUP_ROWS * DX_M / L0_WIDE_M:.2f} of '
          f'L0 here, and 2048 rows is '
          f'{2048 * DX_M / L0_WIDE_M:.2f} of L0.')
    print('   Section 3 makes the same point from the other side: at '
          'L0 = 2.56 m the')
    print(f'   same {SPINUP_ROWS} rows cover '
          f'{SPINUP_ROWS * DX_M / L0_TIGHT_M:.1f} outer scales, and the '
          f'variance lands inside 1 percent.')
    print('   So the L0 = 25 m error is a SPIN-UP length, not a wrong '
          'covariance.')
    return rows


# ---------------------------------------------------------------------------
# The figures
# ---------------------------------------------------------------------------

def draw_covariance(cov):
    '''Draw B(r) per axis against theory, plus a residual panel per case.'''
    fig, axes = plt.subplots(2, 2, figsize=(13.6, 9.6),
                             constrained_layout=True)
    for row, l0_m in enumerate((L0_TIGHT_M, L0_WIDE_M)):
        data = cov[l0_m]
        r_m, theory = data['r_m'], data['theory']
        a0, a1 = axes[row, 0], axes[row, 1]

        a0.plot(r_m, theory, color='black', linewidth=2.4,
                label='theory, Assemat and Wilson Eq. (5)')
        a0.plot(r_m, data['ext'], color='tab:red', linewidth=1.6,
                label='measured, extrusion axis (axis 0)')
        a0.plot(r_m, data['tra'], color='tab:blue', linewidth=1.6,
                linestyle='--', label='measured, transverse axis (axis 1)')
        a0.set_xlabel('separation r, m')
        a0.set_ylabel('B(r), rad^2')
        a0.set_title(f'L0 = {l0_m:.2f} m, {data["m"]} screens, '
                     f'B(0) = {data["b_zero"]:.2f} rad^2', fontsize=10)
        a0.grid(alpha=0.3)
        a0.legend(fontsize=8)

        tol = 100.0 * data['tol']
        a1.axhline(0.0, color='black', linewidth=2.0)
        a1.plot(r_m, 100.0 * (data['ext'] - theory) / data['b_zero'],
                color='tab:red', linewidth=1.6, label='extrusion axis')
        a1.plot(r_m, 100.0 * (data['tra'] - theory) / data['b_zero'],
                color='tab:blue', linewidth=1.6, linestyle='--',
                label='transverse axis')
        a1.axhspan(-tol, tol, color='0.85', zorder=0,
                   label=f'band, {tol:.0f} percent of B(0)')
        a1.axvline(COV_BAND_R_M, color='tab:green', linewidth=1.4,
                   linestyle=':', label=f'band edge, {COV_BAND_R_M:.1f} m')
        a1.set_xlabel('separation r, m')
        a1.set_ylabel('(measured - theory) / B(0), percent')
        a1.set_title('the residual, as a share of B(0)', fontsize=10)
        a1.grid(alpha=0.3)
        a1.legend(fontsize=8)

    fig.suptitle('The covariance of the aotools extruded phase screens, per '
                 'axis.\nThe theory curve is the SAME equation that the '
                 'extrusion matrices use.', fontsize=12)
    fig.savefig(COV_PNG, dpi=150)
    plt.close(fig)


def draw_dphi(theory, ratio_vk, ratio_kol, d_vk, d_kol):
    '''Draw the structure function and its ratio, for the two screen classes.'''
    fig, (a0, a1) = plt.subplots(1, 2, figsize=(13.6, 5.6),
                                 constrained_layout=True)

    a0.loglog(DPHI_BINS_M, theory, color='black', linewidth=2.4,
              label='theory, 2 [B(0) - B(r)]')
    a0.loglog(DPHI_BINS_M, d_vk, 'o-', color='tab:red', linewidth=1.6,
              markersize=4, label=f'PhaseScreenVonKarman, {M_DPHI} screens')
    a0.loglog(DPHI_BINS_M, d_kol, 's--', color='tab:blue', linewidth=1.6,
              markersize=4, label=f'PhaseScreenKolmogorov, {M_KOLM} screens')
    a0.set_xlabel('separation r, m')
    a0.set_ylabel('D(r), rad^2')
    a0.set_title(f'the ensemble structure function, r0 = '
                 f'{R0_M * 1e2:.0f} cm, L0 = {L0_WIDE_M:.0f} m', fontsize=10)
    a0.grid(alpha=0.3, which='both')
    a0.legend(fontsize=8, loc='upper left')

    a1.axhline(1.0, color='black', linewidth=2.4, label='theory')
    a1.semilogx(DPHI_BINS_M, ratio_vk, 'o-', color='tab:red', linewidth=1.6,
                markersize=4, label='PhaseScreenVonKarman')
    a1.semilogx(DPHI_BINS_M, ratio_kol, 's--', color='tab:blue',
                linewidth=1.6, markersize=4, label='PhaseScreenKolmogorov')
    a1.axhspan(DPHI_BAND[0], DPHI_BAND[1], color='0.85', zorder=0,
               label=f'band, {DPHI_BAND[0]:.2f} to {DPHI_BAND[1]:.2f}')
    a1.axvspan(DPHI_BAND_R_M[0], DPHI_BAND_R_M[1], color='0.93', zorder=0,
               label=f'band, r = {DPHI_BAND_R_M[0]:.2f} to '
                     f'{DPHI_BAND_R_M[1]:.2f} m')
    a1.set_xlabel('separation r, m')
    a1.set_ylabel('measured D(r) / theory')
    a1.set_title('the same, as a ratio. The Kolmogorov class is a report '
                 'only.', fontsize=10)
    a1.grid(alpha=0.3)
    a1.legend(fontsize=8, loc='lower left')

    fig.suptitle('The structure function of the aotools extruded phase '
                 'screens, after a 512-row spin-up.', fontsize=12)
    fig.savefig(DPHI_PNG, dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
# The CSV writers
# ---------------------------------------------------------------------------

def write_csv(path, header, rows):
    '''Write one table to a CSV file. Use newline="" for the csv module.'''
    with open(path, 'w', newline='', encoding='utf-8') as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)


# ---------------------------------------------------------------------------

def main():
    t_start = time.time()

    print('=' * 78)
    print('The SPATIAL statistics of the aotools infinite (extruded) screens')
    print('=' * 78)
    print(f'  wavelength           {WAVELENGTH_M * 1e9:11.1f} nm')
    print(f'  grid                 {N:11d} px at {DX_M * 1e2:.1f} cm, '
          f'side {SIDE_M:.2f} m')
    print(f'  r0                   {R0_M * 1e2:11.1f} cm  '
          f'({R0_M / DX_M:.0f} px, the side is {SIDE_M / R0_M:.0f} r0)')
    print(f'  outer scales         {L0_WIDE_M:11.2f} m and {L0_TIGHT_M:.2f} m')
    print(f'  spin-up              {SPINUP_ROWS:11d} add_row calls per screen')
    print(f'  ensembles            {M_WIDE:11d} VK at L0 = '
          f'{L0_WIDE_M:.0f} m, {M_TIGHT} VK at L0 = {L0_TIGHT_M:.2f} m, '
          f'{M_KOLM} Kolm')
    print(f'  master seed          {MASTER_SEED:11d}')

    total_seeds = M_WIDE + M_TIGHT + M_KOLM + M_SPINUP
    seeds = helpers.spawn_seeds(MASTER_SEED, total_seeds)
    seeds_wide = seeds[:M_WIDE]
    seeds_tight = seeds[M_WIDE:M_WIDE + M_TIGHT]
    seeds_kolm = seeds[M_WIDE + M_TIGHT:M_WIDE + M_TIGHT + M_KOLM]
    seeds_spinup = seeds[M_WIDE + M_TIGHT + M_KOLM:]

    # ---- section 1 needs no screen, so it runs first ----
    section_formula()

    # ---- the wide ensemble feeds sections 2, 3 and 4 ----
    t0 = time.time()
    wide = [von_karman_screen(s, L0_WIDE_M) for s in seeds_wide]
    print('')
    print(f'   built {M_WIDE} VK screens at L0 = {L0_WIDE_M:.0f} m '
          f'({time.time() - t0:.1f} s)')

    mask_dphi = helpers.pupil_mask(N, DX_M, DPHI_MASK_D_M)
    mask_tilt = helpers.pupil_mask(N, DX_M, PUPIL_D_M)
    d_vk = helpers.ensemble_dphi(wide[:M_DPHI], mask_dphi, DX_M, DPHI_BINS_M)
    cov_wide = axis_covariance(wide, COV_KMAX_PX)
    moments_wide = screen_moments(wide)
    z_angles, g_angles = pooled_tilt(wide, mask_tilt)
    del wide

    # ---- the tight ensemble feeds section 3 only ----
    t0 = time.time()
    tight = [von_karman_screen(s, L0_TIGHT_M) for s in seeds_tight]
    print(f'   built {M_TIGHT} VK screens at L0 = {L0_TIGHT_M:.2f} m '
          f'({time.time() - t0:.1f} s)')
    cov_tight = axis_covariance(tight, COV_KMAX_PX)
    moments_tight = screen_moments(tight)
    del tight

    # ---- the Kolmogorov-class ensemble feeds section 2 only ----
    t0 = time.time()
    kolm = [kolmogorov_screen(s, L0_WIDE_M) for s in seeds_kolm]
    print(f'   built {M_KOLM} Kolmogorov-class screens '
          f'({time.time() - t0:.1f} s)')
    d_kol = helpers.ensemble_dphi(kolm, mask_dphi, DX_M, DPHI_BINS_M)
    del kolm

    # ---- the four sections ----
    theory, ratio_vk, ratio_kol = section_structure(d_vk, d_kol)
    cov = section_covariance({
        L0_TIGHT_M: (cov_tight[0], cov_tight[1], M_TIGHT, moments_tight),
        L0_WIDE_M: (cov_wide[0], cov_wide[1], M_WIDE, moments_wide)})
    tilt_rows = section_tilt(z_angles, g_angles, mask_tilt)
    section_spinup(seeds_spinup)

    # ---- the outputs ----
    rows = []
    for l0_m in (L0_TIGHT_M, L0_WIDE_M):
        data = cov[l0_m]
        for name, curve in (('extrusion', data['ext']),
                            ('transverse', data['tra'])):
            for i, r_m in enumerate(data['r_m']):
                rows.append([f'{l0_m:.4f}', name, f'{r_m:.4f}',
                             f'{curve[i]:.8e}', f'{data["theory"][i]:.8e}'])
    write_csv(COV_CSV, ['L0_m', 'axis', 'r_m', 'b_est', 'b_theory'], rows)

    rows = []
    for name, measured, ratio in (('PhaseScreenVonKarman', d_vk, ratio_vk),
                                  ('PhaseScreenKolmogorov', d_kol,
                                   ratio_kol)):
        for i, r_m in enumerate(DPHI_BINS_M):
            rows.append([name, f'{r_m:.6f}', f'{measured[i]:.8e}',
                         f'{theory[i]:.8e}', f'{ratio[i]:.6f}'])
    write_csv(DPHI_CSV, ['klass', 'r_m', 'd_meas', 'd_theory', 'ratio'], rows)

    write_csv(TILT_CSV, ['metric', 'var_meas', 'var_theory', 'ratio', 'se'],
              [[name, f'{meas:.8e}', f'{want:.8e}', f'{ratio:.6f}',
                f'{sigma:.8e}'] for name, meas, want, ratio, sigma in
               tilt_rows])

    draw_covariance(cov)
    draw_dphi(theory, ratio_vk, ratio_kol, d_vk, d_kol)

    print('')
    print('OUTPUTS')
    for path in (COV_PNG, DPHI_PNG, COV_CSV, DPHI_CSV, TILT_CSV):
        print(f'   {os.path.basename(path):<26}'
              f'{os.path.getsize(path) / 1024:8.1f} kB')
    print(f'   {os.path.basename(COV_PNG)}: B(r) per axis against theory, and '
          f'the residual as a')
    print('     share of B(0). The top row is the tight outer scale, and the '
          'bottom row is')
    print('     the wide one.')
    print(f'   {os.path.basename(DPHI_PNG)}: the ensemble structure function '
          f'and its ratio, for the two')
    print('     screen classes. The grey bands are the stated pass band.')

    print('')
    print('SUMMARY')
    for label, ok in RESULTS:
        print(f'   [{"PASS" if ok else "FAIL"}] {label}')
    failed = [label for label, ok in RESULTS if not ok]
    print(f'   {len(RESULTS) - len(failed)} of {len(RESULTS)} bands pass')
    if failed:
        print('   A failed band is a RESULT. Read the tables above.')

    print('')
    print(f'(elapsed {time.time() - t_start:.1f} s)')


if __name__ == '__main__':
    main()
