'''
The shared truth and the shared estimators of the phase-screen study.

A Fourier phase screen holds no power below the grid fundamental 1/(N dx). So
its low-order modes are too weak, and its TILT is the weakest of them all. This
module holds the analytic truth that a screen must match, and the estimators
that measure a screen against that truth. The study scripts import it. This
module builds no screen, and it writes no figure.

THE FOUR GROUPS.

  1. THE GRID HELPERS. `crop_center`, `pupil_mask` and `radial_average` put a
     screen, a pupil and a two-dimensional map on ONE centre convention: the
     grid centre sits at the index n//2. That is the convention of
     `schmidt.fourier.structure_function`, so every estimate below reads the
     same origin.
  2. THE SEEDS. `spawn_seeds` and `spawn_rngs` split one master seed into
     independent streams. Two study scripts that share a master seed then see
     the same atmosphere.
  3. THE TILT. `zernike_tilt` and `gradient_tilt` are the two standard
     definitions of the tilt of a wavefront over a circular pupil. Z-tilt is
     the least-squares plane, and G-tilt is the mean gradient. They are NOT the
     same number. `tilt_filter_variance` and `gtilt_filter_variance` give the
     analytic variance of each one. `captured_fraction` gives the share of the
     Z-tilt variance that a screen of a given side can hold.
  4. THE COVARIANCE. `vk_covariance_numeric` and `vk_covariance_closed` are two
     independent routes to the von Karman phase covariance B(r). The first
     integrates the PSD. The second evaluates a closed form. Their agreement
     bounds the error of the integration grid.

THE CONVENTION. Every spatial frequency is an ORDINARY frequency f [1/m], not
an angular frequency kappa [rad/m]. `schmidt.turbulence` uses the same
convention, so a PSD passes between the two modules with no factor.

UNITS. Every length is a metre. A frequency is 1/m. A phase is a radian. A tilt
ANGLE is a radian. A Noll coefficient is a radian. This module returns no
decibels.

THE INTEGRATION GRID. Every frequency integral runs on ONE log-spaced grid,
`F_GRID`, with the trapezium rule. A Gauss-Kronrod quadrature such as
`scipy.integrate.quad` fails here, because the Bessel filters oscillate over
five decades of frequency. See the `F_GRID` block for the measured truncation
error.

Sources:
- Schmidt, Numerical Simulation of Optical Wave Propagation with Examples in
  MATLAB, DOI 10.1117/3.866274. Ch. 3, Eqs. (3.19) to (3.25), printed pp. 49
  and 50 (the structure-function estimator); Ch. 9, Eq. (9.44), printed p. 160
  (D(r) = 6.88 (r/r0)^(5/3)); Ch. 9, Eqs. (9.50) to (9.52), printed p. 161 (the
  phase PSDs).
- Noll, "Zernike polynomials and atmospheric turbulence", J. Opt. Soc. Am.
  66(3), pp. 207 to 211 (1976), DOI 10.1364/JOSA.66.000207. The Zernike
  spatial filter and the residual table.
- Andrews and Phillips, Laser Beam Propagation through Random Media, 2nd ed.,
  DOI 10.1117/3.626196. Ch. 6, Eqs. (80) to (84), printed pp. 200 and 201. The
  gradient (G) tilt.
- Assemat and Wilson, "Method for simulating infinitely long and non stationary
  phase screens with optimized memory storage", Opt. Express 14(3), pp. 988 to
  999 (2006), DOI 10.1364/OE.14.000988. Eq. (5), the closed-form von Karman
  phase covariance.

Run from the repo root:
    python -m validation.screens.helpers
'''

import numpy as np
from scipy.special import gamma as _gamma
from scipy.special import j0 as _j0
from scipy.special import jv as _jv
from scipy.special import kv as _kv

from olb.waveoptics.schmidt.fourier import structure_function
from olb.waveoptics.schmidt.turbulence import (kolmogorov_phase_psd,
                                               kolmogorov_structure_function,
                                               von_karman_phase_psd)

# numpy renamed `trapz` to `trapezoid` in version 2.0. Keep both.
_trapz = getattr(np, 'trapezoid', None) or np.trapz


# ---------------------------------------------------------------------------
# The module constants
# ---------------------------------------------------------------------------

# ---- the frequency integration grid ----
# Every integral of this module runs on this grid. The band covers 13 decades,
# so the Kolmogorov PSD f^(-11/3) is inside its numerical range at both ends.
#
# THE TRUNCATION. The Z-tilt integrand falls as f^(-2/3) below 1/D, so the
# omitted band 0 to F_MIN carries about 0.3 percent of the Z-tilt variance. The
# omitted band above F_MAX carries less than 1e-7. The self-check measures the
# total against the Noll constant, so the truncation is inside the stated band.
F_MIN = 1e-8
F_MAX = 1e5
F_COUNT = 400000
F_GRID = np.logspace(np.log10(F_MIN), np.log10(F_MAX), F_COUNT)

# ---- the analytic tilt constants ----
# Noll (1976), DOI 10.1364/JOSA.66.000207, Table IV, p. 209. The residual after
# piston is Delta_1 = 1.0299 (D/r0)^(5/3), and after one tilt axis it is
# Delta_2 = 0.582. So ONE tilt axis carries 1.0299 - 0.582 = 0.4479.
NOLL_ZTILT_A2 = 0.4479

# The same Z-tilt as an ANGLE variance. Divide by (k D / 4)^2 and the constant
# becomes 0.4479 * 16 / (2 pi)^2 = 0.182.
NOLL_ZTILT_ANGLE = 0.182

# Andrews and Phillips, DOI 10.1117/3.626196, Ch. 6, Eq. (84), printed p. 201.
# The per-axis G-tilt ANGLE variance is 0.174 (D/r0)^(5/3) (lambda/D)^2.
AP_GTILT_ANGLE = 0.174


# ---------------------------------------------------------------------------
# 1. The grid helpers
# ---------------------------------------------------------------------------

def crop_center(a, n):
    '''
    Return the central n by n block of a square array.

    Parameters:
        a : numpy.ndarray
            A square two-dimensional array.
        n : int
            The side of the block. It must not be larger than the input side.

    Returns:
        numpy.ndarray
            A view of the central block.

    THE CENTRE. The input centre sits at the index N//2, and the output centre
    sits at the index n//2. This function keeps that sample. So a crop does not
    move the origin of `schmidt.fourier.structure_function`.
    '''
    a = np.asarray(a)
    n = int(n)
    if a.ndim != 2 or a.shape[0] != a.shape[1]:
        raise ValueError('crop_center needs a square two-dimensional array')
    if n > a.shape[0]:
        raise ValueError('crop_center cannot grow an array')
    start = a.shape[0] // 2 - n // 2
    return a[start:start + n, start:start + n]


def pupil_mask(n, dx_m, d_m):
    '''
    Return a centred circular pupil mask.

    Parameters:
        n : int
            Grid points per side.
        dx_m : float
            Grid pitch [m].
        d_m : float
            Pupil diameter [m].

    Returns:
        numpy.ndarray
            An n by n boolean array. It is True inside the pupil.

    THE CENTRE. The coordinate axis is (arange(n) - n//2) dx, so the origin
    sits at the index n//2. That matches `schmidt.fourier.structure_function`,
    which puts the zero separation at the grid centre.
    '''
    n = int(n)
    axis = (np.arange(n) - n // 2) * float(dx_m)
    xx, yy = np.meshgrid(axis, axis)
    return np.hypot(xx, yy) <= 0.5 * float(d_m)


def radial_average(map2d, dx_m, r_bins_m):
    '''
    Average a two-dimensional map over annuli about the grid centre.

    Parameters:
        map2d : numpy.ndarray
            A square two-dimensional map.
        dx_m : float
            Grid pitch [m].
        r_bins_m : array_like
            The bin CENTRES [m]. The function builds the edges from the
            midpoints between the centres.

    Returns:
        numpy.ndarray
            One mean per bin. A bin with no pixel gives numpy.nan.

    THE EDGES. The inner edge of bin i is the midpoint between centre i-1 and
    centre i. The outer edge of the last bin extends by the same half step. The
    first inner edge is clipped at zero.
    '''
    map2d = np.asarray(map2d, dtype=float)
    n = map2d.shape[0]
    axis = (np.arange(n) - n // 2) * float(dx_m)
    xx, yy = np.meshgrid(axis, axis)
    radius = np.hypot(xx, yy).ravel()

    centres = np.asarray(r_bins_m, dtype=float)
    mid = 0.5 * (centres[:-1] + centres[1:])
    first = max(0.0, centres[0] - 0.5 * (centres[1] - centres[0]))
    last = centres[-1] + 0.5 * (centres[-1] - centres[-2])
    edges = np.concatenate([[first], mid, [last]])

    values = map2d.ravel()
    out = np.full(centres.size, np.nan)
    index = np.digitize(radius, edges) - 1
    for i in range(centres.size):
        keep = index == i
        if np.any(keep):
            out[i] = float(np.mean(values[keep]))
    return out


# ---------------------------------------------------------------------------
# 2. The seeds
# ---------------------------------------------------------------------------

def spawn_seeds(master, count):
    '''
    Split one master seed into independent integer seeds.

    Parameters:
        master : int
            The master seed.
        count : int
            How many seeds to return.

    Returns:
        list
            `count` Python integers.

    Use this for a generator that takes an INTEGER seed, such as the aotools
    screen of `olb/waveoptics/turbulence/screens.py`.
    '''
    state = np.random.SeedSequence(int(master)).generate_state(int(count))
    return [int(s) for s in state]


def spawn_rngs(master, count):
    '''
    Split one master seed into independent numpy generators.

    Parameters:
        master : int
            The master seed.
        count : int
            How many generators to return.

    Returns:
        list
            `count` `numpy.random.Generator` objects.

    Use this for a generator that takes an `rng`, such as the book screens of
    `olb/waveoptics/schmidt/turbulence.py`. The streams are independent, so two
    generators never draw the same numbers.
    '''
    children = np.random.SeedSequence(int(master)).spawn(int(count))
    return [np.random.default_rng(c) for c in children]


# ---------------------------------------------------------------------------
# 3. The tilt: the two estimators
# ---------------------------------------------------------------------------

def _mask_coordinates(mask, dx_m):
    '''Return the x and the y of every True pixel of a centred mask.'''
    mask = np.asarray(mask, dtype=bool)
    n = mask.shape[0]
    axis = (np.arange(n) - n // 2) * float(dx_m)
    xx, yy = np.meshgrid(axis, axis)
    return xx[mask], yy[mask]


def mask_diameter(mask, dx_m):
    '''
    Return the area-equivalent diameter of a pixel mask.

    Parameters:
        mask : numpy.ndarray
            A boolean pupil mask.
        dx_m : float
            Grid pitch [m].

    Returns:
        float
            D = 2 sqrt(A / pi), with A the mask area [m].

    WHY THE AREA. `zernike_tilt` fits over the ACTUAL pixels of the mask, not
    over an ideal circle. So the Noll bridge must read the diameter of the
    same pixel set. A well-sampled circular mask gives the nominal diameter to
    better than 0.5 percent.
    '''
    area = float(np.count_nonzero(mask)) * float(dx_m) ** 2
    return 2.0 * np.sqrt(area / np.pi)


def zernike_tilt(screen_rad, mask, dx_m, wavelength_m):
    '''
    Return the Z-tilt of a phase screen over a pupil.

    Parameters:
        screen_rad : numpy.ndarray
            The phase screen [rad].
        mask : numpy.ndarray
            A boolean pupil mask. Use `pupil_mask`.
        dx_m : float
            Grid pitch [m].
        wavelength_m : float
            Wavelength [m].

    Returns:
        tuple
            (alpha_x, alpha_y, a2, a3). The first pair is the tilt ANGLE on
            each axis [rad]. The second pair is the Noll coefficient of Z2 and
            Z3 [rad].

    formula:
        fit  phi(x,y) = c0 + gx x + gy y   over the pupil, least squares
        alpha_x = gx / k,      alpha_y = gy / k,      k = 2 pi / lambda
        a2 = gx D / 4,         a3 = gy D / 4
    Source: Noll, J. Opt. Soc. Am. 66(3), pp. 207 to 211 (1976),
    DOI 10.1364/JOSA.66.000207. Table I, p. 208, gives Z2 = 2 r cos(theta) and
    Z3 = 2 r sin(theta), with r normalised to the pupil RADIUS D/2. So
    Z2 = 4 x / D, and phi = a2 Z2 gives gx = 4 a2 / D. The bridge to the angle
    is alpha_x = gx / k = 4 a2 / (k D).

    WHY A PLANE FIT IS THE PROJECTION. Z2 and Z3 are the only Zernike modes
    that are linear in x and y. They are orthogonal to every other mode over a
    circular pupil. So the least-squares plane over that pupil IS the Z2 and Z3
    projection. The fit needs no explicit Zernike basis.

    THE DISCRETE PUPIL. The fit reads the actual mask pixel coordinates, so it
    is exact for a phase that is exactly a plane. The diameter D comes from
    `mask_diameter`, the area-equivalent diameter of the same pixel set.

    VALIDITY. The result is ONE realisation. Compare its ENSEMBLE variance with
    `tilt_filter_variance`.
    '''
    mask = np.asarray(mask, dtype=bool)
    phase = np.asarray(screen_rad, dtype=float)
    x, y = _mask_coordinates(mask, dx_m)
    design = np.column_stack([np.ones_like(x), x, y])
    coeff = np.linalg.lstsq(design, phase[mask], rcond=None)[0]

    k = 2.0 * np.pi / float(wavelength_m)
    d_m = mask_diameter(mask, dx_m)
    gx, gy = float(coeff[1]), float(coeff[2])
    return gx / k, gy / k, gx * d_m / 4.0, gy * d_m / 4.0


def gradient_tilt(screen_rad, mask, dx_m, wavelength_m):
    '''
    Return the G-tilt of a phase screen over a pupil.

    Parameters:
        screen_rad : numpy.ndarray
            The phase screen [rad].
        mask : numpy.ndarray
            A boolean pupil mask. Use `pupil_mask`.
        dx_m : float
            Grid pitch [m].
        wavelength_m : float
            Wavelength [m].

    Returns:
        tuple
            (alpha_x, alpha_y), the tilt ANGLE on each axis [rad].

    formula:
        alpha_x = < d(phi)/dx >_pupil / k,   k = 2 pi / lambda
    Source: Andrews and Phillips, Laser Beam Propagation through Random Media,
    2nd ed., DOI 10.1117/3.626196, Ch. 6, Eqs. (80) to (82), printed pp. 200
    and 201. The angle of arrival is the pupil average of the phase gradient,
    divided by the wave number.

    THE GRADIENT AT THE EDGE. `numpy.gradient` uses a central difference in the
    interior of the GRID, and a one-sided difference at the grid edge. The
    pupil leaves a guard band, so every mask pixel takes a central difference.
    This function averages the full-grid gradient over the mask.

    Z-TILT AND G-TILT ARE NOT THE SAME. The Z-tilt is the least-squares plane,
    and the G-tilt is the mean gradient. They agree exactly for a phase that IS
    a plane. For a Kolmogorov phase the G-tilt variance is lower: 0.174 against
    0.182, per axis, in units of (D/r0)^(5/3) (lambda/D)^2.
    '''
    mask = np.asarray(mask, dtype=bool)
    phase = np.asarray(screen_rad, dtype=float)
    # numpy.gradient returns the derivative along axis 0 first. Axis 0 is y.
    grad_y, grad_x = np.gradient(phase, float(dx_m))
    k = 2.0 * np.pi / float(wavelength_m)
    return float(np.mean(grad_x[mask])) / k, float(np.mean(grad_y[mask])) / k


# ---------------------------------------------------------------------------
# 4. The tilt: the analytic variances
# ---------------------------------------------------------------------------

def _zernike_tilt_integrand(psd_func, d_m):
    '''Return the Z-tilt integrand on `F_GRID`. See `tilt_filter_variance`.'''
    x = np.pi * float(d_m) * F_GRID
    filt = 2.0 * _jv(2, x) / x                 # Noll (1976), Eq. (8), p. 208.
    return 2.0 * np.pi * F_GRID * psd_func(F_GRID) * 2.0 * filt ** 2


def tilt_filter_variance(psd_func, d_m):
    '''
    Return the per-axis Noll Z-tilt coefficient variance.

    Parameters:
        psd_func : callable
            It takes f [1/m] and it returns the PHASE PSD [rad^2 m^2]. Pass a
            lambda that closes over r0 and L0.
        d_m : float
            Pupil diameter [m].

    Returns:
        float
            <a2^2> [rad^2]. The Z3 variance is the same number.

    formula:
        <a2^2> = INTEGRAL 2 pi f Phi(f) 2 [2 J2(pi D f) / (pi D f)]^2 df
    Source: Noll, J. Opt. Soc. Am. 66(3), pp. 207 to 211 (1976),
    DOI 10.1364/JOSA.66.000207. Eq. (8), p. 208, gives the Zernike spatial
    filter |Q_j(f)|^2 = (n+1) [2 J_{n+1}(2 pi f R) / (2 pi f R)]^2, with R the
    pupil radius. A tilt has n = 1, so n + 1 = 2 and the Bessel order is 2.
    With R = D/2 the Bessel argument is pi D f. Eq. (13) and Table IV, p. 209,
    give the Kolmogorov answer 0.4479 (D/r0)^(5/3) per axis.

    THE ANGLE. Divide by (k D / 4)^2 to get the tilt ANGLE variance. For the
    Kolmogorov case that gives 0.182 (D/r0)^(5/3) (lambda/D)^2.

    THE QUADRATURE. The integral runs on `F_GRID` with the trapezium rule. Do
    NOT use `scipy.integrate.quad` here. The Bessel filter oscillates over five
    decades of frequency, and an adaptive rule misses the band.

    VALIDITY. The filter assumes a FULL circular pupil. An obscured pupil needs
    another filter. The result is the variance over an ensemble, not one draw.
    '''
    return float(_trapz(_zernike_tilt_integrand(psd_func, d_m), F_GRID))


def gtilt_filter_variance(psd_func, d_m, wavelength_m):
    '''
    Return the per-axis G-tilt ANGLE variance.

    Parameters:
        psd_func : callable
            It takes f [1/m] and it returns the PHASE PSD [rad^2 m^2].
        d_m : float
            Pupil diameter [m].
        wavelength_m : float
            Wavelength [m].

    Returns:
        float
            The angle-of-arrival variance on one axis [rad^2].

    formula:
        <alpha^2> = (1/k^2) INTEGRAL 2 pi f Phi(f) [(2 pi f)^2 / 2]
                    [2 J1(pi D f) / (pi D f)]^2 df
    Source: Andrews and Phillips, Laser Beam Propagation through Random Media,
    2nd ed., DOI 10.1117/3.626196, Ch. 6, Eqs. (82) to (84), printed pp. 200
    and 201. The pupil average of the gradient gives the aperture filter
    [2 J1(pi D f) / (pi D f)]^2, which is the transform of the circular
    pupil. The factor (2 pi f)^2 / 2 is the gradient on ONE axis: the full
    gradient gives (2 pi f)^2, and the azimuth average splits it in two. The
    closed form is 0.174 (D/r0)^(5/3) (lambda/D)^2.

    THE QUADRATURE. See `tilt_filter_variance`.

    THE MEASURED CONSTANT. The self-check reads 0.1694, which is 2.7 percent
    below the printed 0.174 once the 0.44 percent of the PSD constant is
    removed. Sasiela, Electromagnetic Wave Propagation in Turbulence, 2nd ed.,
    DOI 10.1007/978-3-642-59022-0, gives 0.1698 for the same quantity. So the
    integral agrees with Sasiela, and the book constant is a rounded value. The
    self-check holds a 5 percent band for that reason.

    VALIDITY. The filter assumes a FULL circular pupil. The result is an
    ensemble variance.
    '''
    k = 2.0 * np.pi / float(wavelength_m)
    x = np.pi * float(d_m) * F_GRID
    filt = 2.0 * _jv(1, x) / x
    grad = (2.0 * np.pi * F_GRID) ** 2 / 2.0
    integrand = 2.0 * np.pi * F_GRID * psd_func(F_GRID) * grad * filt ** 2
    return float(_trapz(integrand, F_GRID)) / (k * k)


def captured_fraction(psd_func, d_m, f_low):
    '''
    Return the share of the Z-tilt variance that a finite screen holds.

    Parameters:
        psd_func : callable
            It takes f [1/m] and it returns the PHASE PSD [rad^2 m^2].
        d_m : float
            Pupil diameter [m].
        f_low : float
            The lowest frequency that the screen carries [1/m]. For a square
            screen of side L, take f_low = 1 / L.

    Returns:
        float
            The captured share, from 0 to 1.

    formula:
        fraction = INTEGRAL_{f_low}^{inf} I(f) df / INTEGRAL_0^{inf} I(f) df
        I(f) = 2 pi f Phi(f) 2 [2 J2(pi D f) / (pi D f)]^2
    Source: the integrand is `tilt_filter_variance`, from Noll,
    DOI 10.1364/JOSA.66.000207, Eq. (8), p. 208. The lower limit models the
    grid fundamental of a Fourier screen: Schmidt (2010),
    DOI 10.1117/3.866274, Ch. 9, text below Listing 9.2, printed p. 167,
    states that the grid cannot sample a frequency below 1/(N dx).

    THE MODEL IS A SHARP CUTOFF. A real Fourier screen does not remove the low
    band and keep the rest. It replaces the band by one lumped sample at the
    grid fundamental, and the subharmonic method of Lane, Glindemann and Dainty
    (DOI 10.1088/0959-7174/2/3/003) puts part of the band back. So this
    fraction is a LOWER bound on what a subharmonic screen holds, and it is the
    right model for a plain Fourier screen.

    VALIDITY. The denominator carries the truncation of `F_GRID`. See the
    `F_GRID` block.
    '''
    integrand = _zernike_tilt_integrand(psd_func, d_m)
    total = _trapz(integrand, F_GRID)
    keep = F_GRID >= float(f_low)
    part = _trapz(integrand[keep], F_GRID[keep])
    return float(part / total)


# ---------------------------------------------------------------------------
# 5. The von Karman phase covariance, two ways
# ---------------------------------------------------------------------------

def vk_covariance_numeric(r_m, r0_m, L0_m):
    '''
    Return the von Karman phase covariance by numerical integration.

    Parameters:
        r_m : array_like
            Separation [m].
        r0_m : float
            Fried parameter [m].
        L0_m : float
            Outer scale [m].

    Returns:
        numpy.ndarray
            B(r) [rad^2].

    formula:
        B(r) = INTEGRAL 2 pi f Phi_vk(f) J0(2 pi f r) df
    Source: the Hankel transform pair of a two-dimensional isotropic spectrum.
    Schmidt (2010), DOI 10.1117/3.866274, Ch. 3, Eq. (3.16), printed p. 48,
    relates the structure function to the covariance for an isotropic field.
    The PSD is Ch. 9, Eq. (9.50), printed p. 161, through
    `schmidt.turbulence.von_karman_phase_psd`.

    THE LOOP. The function loops over r. An r array of this study is short, so
    the loop costs less memory than an outer product of 400000 columns.

    VALIDITY. The result carries the truncation of `F_GRID`. At r = 10 m the
    omitted tail above F_MAX is below 1e-6 of B(r), because the integrand falls
    as f^(-8/3). Compare with `vk_covariance_closed`.
    '''
    r = np.atleast_1d(np.asarray(r_m, dtype=float))
    weight = 2.0 * np.pi * F_GRID * von_karman_phase_psd(F_GRID, r0_m, L0_m)
    out = np.array([float(_trapz(weight * _j0(2.0 * np.pi * F_GRID * ri),
                                 F_GRID)) for ri in r])
    return out


def vk_covariance_closed(r_m, r0_m, L0_m):
    '''
    Return the von Karman phase covariance in closed form.

    Parameters:
        r_m : array_like
            Separation [m].
        r0_m : float
            Fried parameter [m].
        L0_m : float
            Outer scale [m].

    Returns:
        numpy.ndarray
            B(r) [rad^2].

    formula:
        B(r) = (L0/r0)^(5/3) Gamma(11/6) / (2^(5/6) pi^(8/3))
               [(24/5) Gamma(6/5)]^(5/6) x^(5/6) K_{5/6}(x),   x = 2 pi r / L0
    Source: Assemat and Wilson, Opt. Express 14(3), pp. 988 to 999 (2006),
    DOI 10.1364/OE.14.000988, Eq. (5).

    THE LIMIT AT ZERO. K_{5/6}(x) diverges at x = 0, and the product
    x^(5/6) K_{5/6}(x) does not. The small-argument form
    K_nu(x) -> Gamma(nu) (2/x)^nu / 2 gives the limit
    x^(5/6) K_{5/6}(x) -> 2^(-1/6) Gamma(5/6). This function uses that value at
    x = 0.

    THE PRECISION. `aotools.turbulence.turb.phase_covariance` implements the
    same equation. It casts r to float32 first, and it adds 1e-40 to remove the
    zero. This function stays in float64 and it handles the zero exactly.

    THE CONSTANT. The prefactor here is exact. The route through
    `schmidt.turbulence.von_karman_phase_psd` carries the PRINTED constant
    0.023, which is 0.44 percent above the exact 0.49 (2 pi)^(-5/3). So
    `vk_covariance_numeric` reads about 0.44 percent higher than this function.
    That gap is the rounding of a printed constant, not an error.

    VALIDITY. It is the von Karman spectrum, so l0 = 0. The outer scale must be
    finite: the expression diverges as L0 goes to infinity.
    '''
    r = np.atleast_1d(np.asarray(r_m, dtype=float))
    prefactor = ((float(L0_m) / float(r0_m)) ** (5.0 / 3.0)
                 * _gamma(11.0 / 6.0) / (2.0 ** (5.0 / 6.0) * np.pi ** (8 / 3))
                 * ((24.0 / 5.0) * _gamma(6.0 / 5.0)) ** (5.0 / 6.0))

    x = 2.0 * np.pi * r / float(L0_m)
    limit = 2.0 ** (-1.0 / 6.0) * _gamma(5.0 / 6.0)
    safe = np.where(x > 0.0, x, 1.0)
    bessel = np.where(x > 0.0,
                      safe ** (5.0 / 6.0) * _kv(5.0 / 6.0, safe), limit)
    return prefactor * bessel


# ---------------------------------------------------------------------------
# 6. The structure-function estimators
# ---------------------------------------------------------------------------

def ensemble_dphi(screens, mask, dx_m, r_bins_m):
    '''
    Return the ensemble mean structure function, binned against separation.

    Parameters:
        screens : iterable
            The phase screens [rad]. Each one is square and real.
        mask : numpy.ndarray
            The pupil window. Use `pupil_mask`.
        dx_m : float
            Grid pitch [m].
        r_bins_m : array_like
            The separation bin centres [m].

    Returns:
        numpy.ndarray
            D(r) [rad^2], one value per bin.

    Source: the estimator is `schmidt.fourier.structure_function`, from Schmidt
    (2010), DOI 10.1117/3.866274, Ch. 3, Eqs. (3.19) to (3.25), printed pp. 49
    and 50. The theory target is Ch. 9, Eq. (9.44), printed p. 160.

    VALIDITY. The estimator makes the correlation CIRCULAR, so the mask must
    leave a guard band at the grid edge. The estimate holds out to the mask
    diameter only, because the window overlap area falls to zero there.
    '''
    mask = np.asarray(mask, dtype=float)
    total = None
    count = 0
    for screen in screens:
        one = structure_function(screen, mask, dx_m)
        total = one if total is None else total + one
        count += 1
    if count == 0:
        raise ValueError('ensemble_dphi needs at least one screen')
    return radial_average(total / count, dx_m, r_bins_m)


def d_phi_direct(scr, kpx):
    '''
    Return the structure function at one separation, by shifted differences.

    Parameters:
        scr : numpy.ndarray
            The phase screen [rad].
        kpx : int
            The separation in pixels.

    Returns:
        float
            D(kpx dx) [rad^2], averaged over the two axes.

    formula:
        D = 0.5 [ < (phi(x+k) - phi(x))^2 > + < (phi(y+k) - phi(y))^2 > ]
    Source: the definition of the structure function, Schmidt (2010),
    DOI 10.1117/3.866274, Ch. 3, Eq. (3.15), printed p. 47. The theory target
    is Ch. 9, Eq. (9.44), printed p. 160.

    THIS IS THE SECOND ESTIMATOR. The self-check of
    `olb/waveoptics/turbulence/screens.py` uses this form. It reads the WHOLE
    grid, and it applies no window. `ensemble_dphi` applies a pupil window and
    it uses the Fourier estimator. The two give slightly different numbers for
    the same screen. Compare each one against the theory, not against the
    other.

    VALIDITY. It reads one realisation. Average it over an ensemble.
    '''
    scr = np.asarray(scr, dtype=float)
    kpx = int(kpx)
    dh = scr[:, kpx:] - scr[:, :-kpx]
    dv = scr[kpx:, :] - scr[:-kpx, :]
    return 0.5 * (float(np.mean(dh * dh)) + float(np.mean(dv * dv)))


# ---------------------------------------------------------------------------

if __name__ == '__main__':
    import time

    t_start = time.time()
    LAMBDA_M = 1550e-9
    K_RAD_M = 2.0 * np.pi / LAMBDA_M

    print('=' * 78)
    print('validation.screens.helpers self-check')
    print('=' * 78)
    print(f'  integration grid   {F_COUNT:11d} points, {F_MIN:.0e} to '
          f'{F_MAX:.0e} 1/m, log spaced')
    print(f'  wavelength         {LAMBDA_M * 1e9:11.1f} nm')

    # ---- 1. the two tilt estimators on an exact plane ----
    # A plane phase a + b x + c y has ONE tilt. The least-squares plane is
    # exact, and the mean gradient is exact. So the two must agree.
    N, DX = 256, 0.01                              # a 2.56 m side
    D_PUPIL = 1.0
    mask = pupil_mask(N, DX, D_PUPIL)
    axis = (np.arange(N) - N // 2) * DX
    XX, YY = np.meshgrid(axis, axis)

    B_RAD_M, C_RAD_M = 3.0, -1.75                  # the plane slopes [rad/m]
    plane = 0.7 + B_RAD_M * XX + C_RAD_M * YY
    want_x, want_y = B_RAD_M / K_RAD_M, C_RAD_M / K_RAD_M

    zx, zy, a2, a3 = zernike_tilt(plane, mask, DX, LAMBDA_M)
    gx, gy = gradient_tilt(plane, mask, DX, LAMBDA_M)
    d_eff = mask_diameter(mask, DX)

    print('')
    print('1. THE TWO TILT ESTIMATORS ON AN EXACT PLANE')
    print(f'   {"quantity":<26}{"exact":>15}{"measured":>15}{"rel err":>12}')
    rows = [('Z-tilt alpha_x, urad', want_x * 1e6, zx * 1e6),
            ('Z-tilt alpha_y, urad', want_y * 1e6, zy * 1e6),
            ('G-tilt alpha_x, urad', want_x * 1e6, gx * 1e6),
            ('G-tilt alpha_y, urad', want_y * 1e6, gy * 1e6),
            ('Noll a2, rad', B_RAD_M * d_eff / 4.0, a2),
            ('Noll a3, rad', C_RAD_M * d_eff / 4.0, a3)]
    for name, want, got in rows:
        print(f'   {name:<26}{want:>15.6f}{got:>15.6f}'
              f'{abs(got / want - 1.0):>12.2e}')
    print(f'   nominal pupil diameter    {D_PUPIL:>15.6f} m')
    print(f'   area-equivalent diameter  {d_eff:>15.6f} m'
          f'   (mask_diameter)')

    err_z = max(abs(zx / want_x - 1.0), abs(zy / want_y - 1.0))
    err_g = max(abs(gx / want_x - 1.0), abs(gy / want_y - 1.0))
    assert err_z < 1e-12, err_z
    assert err_g < 1e-9, err_g
    # The Noll bridge alpha_x = 4 a2 / (k D) must close on the same D.
    assert abs(4.0 * a2 / (K_RAD_M * d_eff) / zx - 1.0) < 1e-12
    assert abs(4.0 * a3 / (K_RAD_M * d_eff) / zy - 1.0) < 1e-12

    # ---- 2. the Z-tilt filter against the Noll constant ----
    R0_M = 0.10
    D_OVER_R0 = D_PUPIL / R0_M
    kol = lambda f: kolmogorov_phase_psd(f, R0_M)

    want_a2 = NOLL_ZTILT_A2 * D_OVER_R0 ** (5.0 / 3.0)
    got_a2 = tilt_filter_variance(kol, D_PUPIL)
    want_angle = NOLL_ZTILT_ANGLE * D_OVER_R0 ** (5.0 / 3.0) \
        * (LAMBDA_M / D_PUPIL) ** 2
    got_angle = got_a2 / (K_RAD_M * D_PUPIL / 4.0) ** 2

    print('')
    print('2. THE ZERNIKE TILT FILTER, Noll (1976) Eq. (8) and Table IV')
    print(f'   Kolmogorov, r0 = {R0_M * 1e2:.0f} cm, D = {D_PUPIL:.2f} m, '
          f'D/r0 = {D_OVER_R0:.1f}')
    print(f'   {"quantity":<26}{"Noll":>15}{"integral":>15}{"rel err":>12}')
    print(f'   {"<a2^2>, rad^2":<26}{want_a2:>15.4f}{got_a2:>15.4f}'
          f'{abs(got_a2 / want_a2 - 1.0):>12.2e}')
    print(f'   {"<alpha^2>, rad^2":<26}{want_angle:>15.4e}'
          f'{got_angle:>15.4e}{abs(got_angle / want_angle - 1.0):>12.2e}')
    assert abs(got_a2 / want_a2 - 1.0) < 0.01, (got_a2, want_a2)

    # ---- 3. the G-tilt filter against Andrews and Phillips ----
    want_g = AP_GTILT_ANGLE * D_OVER_R0 ** (5.0 / 3.0) \
        * (LAMBDA_M / D_PUPIL) ** 2
    got_g = gtilt_filter_variance(kol, D_PUPIL, LAMBDA_M)

    print('')
    print('3. THE GRADIENT TILT FILTER, Andrews and Phillips Ch. 6, Eq. (84)')
    print(f'   {"quantity":<26}{"book":>15}{"integral":>15}{"rel err":>12}')
    print(f'   {"G-tilt <alpha^2>, rad^2":<26}{want_g:>15.4e}{got_g:>15.4e}'
          f'{abs(got_g / want_g - 1.0):>12.2e}')
    print(f'   {"Z-tilt <alpha^2>, rad^2":<26}{want_angle:>15.4e}'
          f'{got_angle:>15.4e}{"":>12}')
    print(f'   G-tilt over Z-tilt, measured  {got_g / got_angle:>13.4f}'
          f'   (0.174 / 0.182 = '
          f'{AP_GTILT_ANGLE / NOLL_ZTILT_ANGLE:.4f})')
    assert abs(got_g / want_g - 1.0) < 0.05, (got_g, want_g)
    assert got_g < got_angle, (got_g, got_angle)

    # ---- 4. the two covariance routes ----
    L0_M = 25.0
    r_probe = np.logspace(np.log10(0.05), np.log10(10.0), 9)
    b_num = vk_covariance_numeric(r_probe, R0_M, L0_M)
    b_cls = vk_covariance_closed(r_probe, R0_M, L0_M)

    print('')
    print('4. THE VON KARMAN COVARIANCE, TWO ROUTES')
    print(f'   r0 = {R0_M * 1e2:.0f} cm, L0 = {L0_M:.0f} m')
    print(f'   {"r, m":>10}{"numeric":>14}{"closed":>14}{"ratio":>10}')
    for i, ri in enumerate(r_probe):
        print(f'   {ri:>10.4f}{b_num[i]:>14.4f}{b_cls[i]:>14.4f}'
              f'{b_num[i] / b_cls[i]:>10.5f}')
    ratio_b = b_num / b_cls
    print(f'   the numeric route reads {100.0 * (ratio_b.mean() - 1.0):+.2f} '
          f'percent high, from the printed constant 0.023')
    assert np.all(np.abs(ratio_b - 1.0) < 0.015), ratio_b

    # ---- 5. the covariance gives the Kolmogorov structure function ----
    # D(r) = 2 [B(0) - B(r)]. Schmidt (2010), DOI 10.1117/3.866274, Ch. 3,
    # Eq. (3.16), printed p. 48. A large L0 must give Ch. 9, Eq. (9.44).
    #
    # THE OUTER-SCALE CORRECTION. The small-argument series of the Bessel
    # function is x^nu K_nu(x) = C0 - Q x^(2nu) + C0 x^2 / [4 (1 - nu)] - ...
    # with nu = 5/6, C0 = 2^(-1/6) Gamma(5/6) and Q = pi / [2^(5/6)
    # Gamma(11/6)]. The x^(5/3) term IS the Kolmogorov law, and the x^2 term is
    # the first outer-scale correction. Their ratio is
    #
    #     D(r) / D_kol(r) = 1 - VK_DEFICIT (2 pi r / L0)^(1/3) + O((r/L0)^(4/3))
    #     VK_DEFICIT = 3 C0 / (2 Q) = 0.8049
    #
    # The correction falls only as the CUBE ROOT of the outer scale. So an
    # L0 of 1 km leaves a 15 percent deficit at r = 1 m, and an L0 of 1000 km
    # is needed for 1 percent. The check below uses L0 = 1e6 m, and it also
    # measures the deficit against the coefficient above.
    VK_DEFICIT = (3.0 * 2.0 ** (-1.0 / 6.0) * _gamma(5.0 / 6.0)
                  / (2.0 * np.pi / (2.0 ** (5.0 / 6.0) * _gamma(11.0 / 6.0))))
    # The leading term of the same series IS the Kolmogorov law, and its
    # constant collapses to 2 [(24/5) Gamma(6/5)]^(5/6). The book prints that
    # number rounded to 6.88 (Ch. 9, Eq. (9.44), printed p. 160).
    KOL_D_EXACT = 2.0 * ((24.0 / 5.0) * _gamma(6.0 / 5.0)) ** (5.0 / 6.0)
    L0_BIG_M = 1e6
    r_sf = np.array([0.1, 0.2, 0.4, 0.7, 1.0])
    b_zero = float(vk_covariance_closed(0.0, R0_M, L0_BIG_M)[0])
    d_from_b = 2.0 * (b_zero - vk_covariance_closed(r_sf, R0_M, L0_BIG_M))
    d_kol = kolmogorov_structure_function(r_sf, R0_M)
    d_exact = KOL_D_EXACT * (r_sf / R0_M) ** (5.0 / 3.0)
    ratio_d = d_from_b / d_kol
    coefficient = ((1.0 - d_from_b / d_exact)
                   / (2.0 * np.pi * r_sf / L0_BIG_M) ** (1.0 / 3.0))

    print('')
    print('5. 2 [B(0) - B(r)] AGAINST 6.88 (r/r0)^(5/3)')
    print(f'   r0 = {R0_M * 1e2:.0f} cm, L0 = {L0_BIG_M:.0e} m, '
          f'B(0) = {b_zero:.4e} rad^2')
    print(f'   {"r, m":>10}{"2[B(0)-B(r)]":>16}{"6.88 law":>14}'
          f'{"ratio":>10}{"deficit coeff":>15}')
    for i, ri in enumerate(r_sf):
        print(f'   {ri:>10.3f}{d_from_b[i]:>16.4f}{d_kol[i]:>14.4f}'
              f'{ratio_d[i]:>10.5f}{coefficient[i]:>15.4f}')
    print(f'   the exact structure constant, 2 [(24/5) Gamma(6/5)]^(5/6) = '
          f'{KOL_D_EXACT:.4f}')
    print(f'   the book prints it rounded, 6.88, which is '
          f'{100.0 * (6.88 / KOL_D_EXACT - 1.0):+.3f} percent')
    print('   The deficit is the OUTER SCALE, not an error. The analytic')
    print(f'   coefficient 3 C0 / (2 Q) is {VK_DEFICIT:.4f}, and the measured')
    print('   column above holds it. The correction falls as the CUBE ROOT of')
    print('   L0, so L0 = 1e3 m still leaves a 15 percent deficit at r = 1 m.')
    assert np.all(np.abs(ratio_d - 1.0) < 0.03), ratio_d
    assert np.all(np.abs(coefficient / VK_DEFICIT - 1.0) < 0.05), coefficient

    # ---- 6. the captured fraction of a finite screen ----
    SIDE_M = 5.12                                  # 512 px at 1 cm
    full = captured_fraction(kol, D_PUPIL, F_MIN)
    part = captured_fraction(kol, D_PUPIL, 1.0 / SIDE_M)

    print('')
    print('6. THE CAPTURED SHARE OF THE Z-TILT VARIANCE')
    print(f'   {"f_low, 1/m":>14}{"screen side, m":>18}{"captured":>12}')
    print(f'   {F_MIN:>14.1e}{"unbounded":>18}{full:>12.5f}')
    for side in (2.56, SIDE_M, 20.48, 81.92):
        print(f'   {1.0 / side:>14.4f}{side:>18.2f}'
              f'{captured_fraction(kol, D_PUPIL, 1.0 / side):>12.5f}')
    print(f'   A {SIDE_M:.2f} m screen holds {100.0 * part:.1f} percent of the '
          f'Z-tilt variance of a')
    print(f'   {D_PUPIL:.2f} m pupil, under the sharp-cutoff model. The rest '
          f'sits below 1/side.')
    assert abs(full - 1.0) < 1e-3, full
    assert 0.0 < part < 1.0, part

    # ---- 7. the grid helpers and the two structure-function estimators ----
    # A ramp phi = a x has an exact structure function: D(dr) = (a dr_x)^2 at
    # every separation. Bin the EXACT map through the same `radial_average`, so
    # the comparison carries no binning difference. A bin is an annulus of
    # finite width, so its mean is NOT (a r)^2 / 2, the value at the centre.
    A_RAMP = 3.0
    ramp = A_RAMP * XX
    r_bins = np.array([0.2, 0.3, 0.45, 0.6, 0.8])
    got_bins = ensemble_dphi([ramp, ramp], mask.astype(float), DX, r_bins)
    want_bins = radial_average((A_RAMP * XX) ** 2, DX, r_bins)
    naive = (A_RAMP * r_bins) ** 2 / 2.0

    print('')
    print('7. THE GRID HELPERS AND THE TWO ESTIMATORS')
    print(f'   {"r, m":>10}{"ensemble_dphi":>16}{"exact, binned":>16}'
          f'{"ratio":>10}{"(a r)^2/2":>13}')
    for i, ri in enumerate(r_bins):
        print(f'   {ri:>10.3f}{got_bins[i]:>16.4f}{want_bins[i]:>16.4f}'
              f'{got_bins[i] / want_bins[i]:>10.5f}{naive[i]:>13.4f}')
    print('   The last column is the value at the bin CENTRE. It is not the')
    print('   bin mean, because an annulus of finite width carries <r^2>.')
    assert np.all(np.abs(got_bins / want_bins - 1.0) < 1e-9), got_bins

    # `d_phi_direct` averages the two axes. The ramp has no y slope, so the
    # answer is one half of the x value.
    for kpx in (4, 10, 25):
        got_direct = d_phi_direct(ramp, kpx)
        want_direct = 0.5 * (A_RAMP * kpx * DX) ** 2
        assert abs(got_direct / want_direct - 1.0) < 1e-12, got_direct
    print(f'   d_phi_direct on the same ramp matches 0.5 (a k dx)^2 to 1e-12')

    # The grid helpers.
    block = crop_center(np.arange(64).reshape(8, 8), 4)
    assert block.shape == (4, 4), block.shape
    assert block[2, 2] == np.arange(64).reshape(8, 8)[4, 4]
    area_error = (np.count_nonzero(mask) * DX ** 2
                  / (np.pi * D_PUPIL ** 2 / 4.0) - 1.0)
    assert abs(area_error) < 0.01, area_error
    print(f'   crop_center keeps the centre sample; pupil_mask area is '
          f'{100.0 * area_error:+.2f} percent')

    # The seeds.
    seeds = spawn_seeds(2026, 5)
    assert seeds == spawn_seeds(2026, 5)
    assert seeds != spawn_seeds(2027, 5)
    assert len(set(seeds)) == 5, seeds
    first = spawn_rngs(2026, 3)[0].standard_normal(4)
    assert np.array_equal(first, spawn_rngs(2026, 3)[0].standard_normal(4))
    assert not np.allclose(first, spawn_rngs(2026, 3)[1].standard_normal(4))
    print(f'   spawn_seeds and spawn_rngs are repeatable and independent')

    print('')
    print(f'(elapsed {time.time() - t_start:.1f} s)')
    print('self-check passed')
