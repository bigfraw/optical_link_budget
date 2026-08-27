'''
Atmospheric turbulence for a split-step simulation, from Schmidt Chapter 9.

This module holds the REFERENCE physics of Chapter 9: the phase power spectral
densities, the Monte-Carlo (Fourier) phase screen with its subharmonic
low-frequency compensation, the per-screen Rytov bound, the layered-atmosphere
moment rule, the turbulent sampling bounds, and the Section 9.5 procedure as a
runnable checklist.

Source of every equation:
    J. D. Schmidt, "Numerical Simulation of Optical Wave Propagation with
    Examples in MATLAB", SPIE Press Monograph PM199 (2010).
    DOI: 10.1117/3.866274
Chapter 9, printed pp. 149 to 183. Each function names its chapter, its
equation number, and its printed page.

THE CONVENTION. The book writes the refractive-index spectra in ANGULAR spatial
frequency kappa [rad/m], and it writes the SCREEN spectra in ORDINARY spatial
frequency f [1/m] (Ch. 9, text above Eq. (9.52), printed p. 161). This module
uses ORDINARY frequency everywhere, because `schmidt.fourier` does.

THIS MODULE IS A VALIDATOR. It duplicates, on purpose, physics that
`olb/waveoptics/turbulence/` gets from `aotools`. The point is to have a
second, book-derived implementation to check the production layer against. It
imports nothing from olb outside this sub-package.

UNITS. `r0_m`, `dx_m`, `l0_m`, `L0_m` and every length are metres. A frequency
is 1/m. A phase is radians. This module returns no decibels.

This module holds physics only. It imports numpy, scipy and `schmidt.fourier`.
'''

import numpy as np
from scipy.optimize import lsq_linear

from .fourier import freq_pitch, ift2

# numpy renamed `trapz` to `trapezoid` in version 2.0. Keep both.
_trapz = getattr(np, 'trapezoid', None) or np.trapz

# The per-screen cap of Listing 9.5, printed p. 175. See `screen_rytov_share`
# for the exact quantity that it bounds.
RMAX = 0.1

# The weak-fluctuation threshold on the LOG-AMPLITUDE variance. Ch. 9, text
# below Eq. (9.64), printed p. 163.
WEAK_SIGMA2_CHI = 0.25


# ---------------------------------------------------------------------------
# 1. The phase power spectral densities (Secs. 9.2.3 and 9.3)
# ---------------------------------------------------------------------------

def phase_psd(f, r0_m, L0_m=np.inf, l0_m=0.0):
    '''
    Return the modified von Karman PHASE power spectral density.

    Parameters:
        f : array_like
            Ordinary spatial frequency [1/m]. It is the radial frequency
            sqrt(fx^2 + fy^2).
        r0_m : float
            Fried parameter of the screen [m].
        L0_m : float
            Outer scale [m]. Give numpy.inf for no outer scale.
        l0_m : float
            Inner scale [m]. Give 0.0 for no inner scale.

    Returns:
        numpy.ndarray
            Phi_phi(f) [rad^2 m^2].

    formula:
        Phi_phi(f) = 0.023 r0^(-5/3) exp(-(f/fm)^2) / (f^2 + f0^2)^(11/6)
        fm = 5.92 / (2 pi l0),   f0 = 1 / L0
    Source: Schmidt (2010), DOI 10.1117/3.866274, Ch. 9, Eq. (9.51), printed
    p. 161, gives the angular form 0.49 r0^(-5/3) exp(-kappa^2/kappa_m^2) /
    (kappa^2 + kappa_0^2)^(11/6). Eq. (9.52), printed p. 161, converts the
    Kolmogorov case to ordinary frequency and gives the constant 0.023.
    Listing 9.2, lines 11 to 15, printed p. 167, prints the same expression in
    code, with fm and f0 as above.

    THE CONSTANT. The book prints 0.49 in angular frequency and 0.023 in
    ordinary frequency. A two-dimensional spectrum transforms as
    Phi(f) = (2 pi)^2 Phi(kappa = 2 pi f), so the constant is
    0.49 (2 pi)^(-5/3) = 0.0229. The two printed constants agree.

    THE SCALE FREQUENCIES. Ch. 9, text below Eq. (9.18), printed p. 155, gives
    kappa_m = 5.92/l0 and kappa_0 = 2 pi/L0, both ANGULAR. In ordinary
    frequency these are fm = 5.92/(2 pi l0) and f0 = 1/L0.

    THE THREE SPECTRA. This one expression gives all three of Ch. 9:
    - Kolmogorov, Eq. (9.49): L0 = inf and l0 = 0.
    - von Karman, Eq. (9.50): l0 = 0.
    - Modified von Karman, Eq. (9.51): both scales finite.
    Ch. 9, text below Eq. (9.18), printed p. 155, states that Eq. (9.18)
    reduces to Eq. (9.16) for l0 = 0 and L0 = inf. The same reduction holds
    for the phase spectra.

    VALIDITY.
    - Eq. (9.48), printed p. 160, derives the phase PSD from the index PSD
      through Phi_phi = 2 pi^2 k^2 z Phi_n. That step assumes a PLANE wave in
      WEAK turbulence. The screen is thin.
    - The Kolmogorov branch is valid only in the inertial subrange
      1/L0 << kappa << 1/l0 (Ch. 9, Eq. (9.16), printed p. 155).
    - The value at f = 0 is infinite when L0 is infinite. The screen
      generators set that sample to zero, as Listing 9.2, line 16, printed
      p. 167, does. This function does NOT, because the divergence is real
      physics.
    '''
    f = np.asarray(f, dtype=float)
    fm = np.inf if l0_m <= 0.0 else 5.92 / (2.0 * np.pi * float(l0_m))
    f0 = 0.0 if not np.isfinite(L0_m) else 1.0 / float(L0_m)
    denom = (f * f + f0 * f0) ** (11.0 / 6.0)
    top = 0.023 * float(r0_m) ** (-5.0 / 3.0) * np.exp(-(f / fm) ** 2)
    return np.where(denom > 0.0, top / np.where(denom > 0.0, denom, 1.0),
                    np.inf)


def kolmogorov_phase_psd(f, r0_m):
    '''
    Return the Kolmogorov phase PSD: 0.023 r0^(-5/3) f^(-11/3).

    Source: Schmidt (2010), DOI 10.1117/3.866274, Ch. 9, Eq. (9.52), printed
    p. 161. It is `phase_psd` with L0 = inf and l0 = 0.
    '''
    return phase_psd(f, r0_m, np.inf, 0.0)


def von_karman_phase_psd(f, r0_m, L0_m):
    '''
    Return the von Karman phase PSD. It adds the outer scale.

    Source: Schmidt (2010), DOI 10.1117/3.866274, Ch. 9, Eq. (9.50), printed
    p. 161. It is `phase_psd` with l0 = 0.
    '''
    return phase_psd(f, r0_m, L0_m, 0.0)


def modified_von_karman_phase_psd(f, r0_m, L0_m, l0_m):
    '''
    Return the modified von Karman phase PSD. It adds both scales.

    Source: Schmidt (2010), DOI 10.1117/3.866274, Ch. 9, Eq. (9.51), printed
    p. 161. It is `phase_psd` itself.
    '''
    return phase_psd(f, r0_m, L0_m, l0_m)


def kolmogorov_structure_function(r_m, r0_m):
    '''
    Return the Kolmogorov phase structure function D(r) = 6.88 (r/r0)^(5/3).

    Parameters:
        r_m : array_like
            Separation [m].
        r0_m : float
            Fried parameter [m].

    Returns:
        numpy.ndarray
            D(r) [rad^2].

    formula:
        D(r) = 6.88 (r / r0)^(5/3)
    Source: Schmidt (2010), DOI 10.1117/3.866274, Ch. 9, Eq. (9.44), printed
    p. 160. The constant 6.88 is the r0 definition of Eq. (9.41), printed
    p. 159: D(r0) = 6.88 rad^2.

    VALIDITY. It is the PLANE-wave, Kolmogorov result. It assumes l0 = 0 and
    L0 = infinite (Ch. 9, text below Eq. (9.44), printed p. 160). It is the
    verification target for a phase screen (Sec. 9.5.5, printed p. 180).
    '''
    return 6.88 * (np.asarray(r_m, dtype=float) / float(r0_m)) ** (5.0 / 3.0)


# ---------------------------------------------------------------------------
# 2. The Monte-Carlo phase screens (Sec. 9.3)
# ---------------------------------------------------------------------------

def ft_phase_screen(r0_m, n, dx_m, L0_m=np.inf, l0_m=0.0, rng=None):
    '''
    Return one Fourier-transform (Monte-Carlo) phase screen.

    Parameters:
        r0_m : float
            Fried parameter of the screen [m].
        n : int
            Grid points per side. Use an even count.
        dx_m : float
            Grid pitch [m].
        L0_m : float
            Outer scale [m]. Default is infinite.
        l0_m : float
            Inner scale [m]. Default is zero.
        rng : numpy.random.Generator or None
            The random source. None takes `numpy.random.default_rng()`.

    Returns:
        numpy.ndarray
            An n by n real phase screen [rad].

    formula:
        phi(x,y)      = SUM_n SUM_m c_nm exp[i 2 pi (fxn x + fym y)]
        <|c_nm|^2>    = Phi_phi(fxn, fym) dfx dfy
        c_nm          = (g1 + i g2) sqrt(Phi_phi) df,   g1, g2 ~ N(0,1)
        phi           = Re{ IFT[c] }
    Source: Schmidt (2010), DOI 10.1117/3.866274, Ch. 9, Sec. 9.3. The Fourier
    series is Eq. (9.78), printed p. 167. The coefficient variance is
    Eq. (9.79), printed p. 167, and its FFT form with df = 1/L is Eq. (9.80),
    printed p. 167. Listing 9.2, printed p. 167, prints the same steps in code.
    The method is credited to McGlamery, J. Opt. Soc. Am. 57(3), pp. 293 to
    297 (1967), DOI 10.1364/JOSA.57.000293.

    THE FACTOR TWO. `g1 + i g2` has a variance of 2, not 1. The real part of
    the inverse transform then keeps one half of the power, so the screen
    variance is the one that Eq. (9.79) asks for. Ch. 9, text below Listing
    9.2, printed p. 167, states that the real and the imaginary part give two
    uncorrelated screens.

    THE ZERO FREQUENCY. The generator sets the f = 0 sample of the PSD to
    zero, as Listing 9.2, line 16, printed p. 167, does. That removes the
    piston, which no propagation reads.

    VALIDITY. This screen is NOT accurate on its own. Ch. 9, text below Listing
    9.2, printed p. 167, states that the grid cannot sample low enough
    frequencies, so the low-order modes, and tilt above all, are too weak.
    Figure 9.3, printed p. 169, shows the measured structure function BELOW the
    theory, and the gap grows with the separation. Use `ft_sh_phase_screen`.
    '''
    rng = np.random.default_rng() if rng is None else rng
    n = int(n)
    df = freq_pitch(n, dx_m)
    fx = (np.arange(n) - n // 2) * df
    FX, FY = np.meshgrid(fx, fx)
    psd = phase_psd(np.hypot(FX, FY), r0_m, L0_m, l0_m)
    psd[n // 2, n // 2] = 0.0                # Listing 9.2, line 16, p. 167.

    # Eq. (9.79), printed p. 167, through the FFT form of Eq. (9.80).
    cn = (rng.standard_normal((n, n))
          + 1j * rng.standard_normal((n, n))) * np.sqrt(psd) * df
    # Eq. (9.78), printed p. 167. `ift2` with df = 1 gives the bare series sum.
    return np.real(ift2(cn, 1.0))


def subharmonic_screen(r0_m, n, dx_m, L0_m=np.inf, l0_m=0.0, rng=None, n_p=3):
    '''
    Return the LOW-frequency part of a subharmonic phase screen.

    Parameters:
        r0_m : float
            Fried parameter of the screen [m].
        n : int
            Grid points per side.
        dx_m : float
            Grid pitch [m].
        L0_m : float
            Outer scale [m].
        l0_m : float
            Inner scale [m].
        rng : numpy.random.Generator or None
            The random source.
        n_p : int
            The number of subharmonic grids. The book uses 3.

    Returns:
        numpy.ndarray
            An n by n real phase screen [rad] with a zero mean.

    formula:
        phi_LF(x,y) = SUM_{p=1..Np} SUM_{n=-1..1} SUM_{m=-1..1}
                      c_nm exp[i 2 pi (fxn x + fym y)]
        df_p = 1 / (3^p L),   L = n dx
    Source: Schmidt (2010), DOI 10.1117/3.866274, Ch. 9, Eq. (9.81), printed
    p. 169. Listing 9.3, printed p. 170, prints the same steps in code, with a
    3 by 3 frequency set and Np = 3 grids.

    THE PRIMARY SOURCES. The book implements the method of Lane, Glindemann
    and Dainty, "Simulation of a Kolmogorov phase screen", Waves in Random
    Media 2, pp. 209 to 224 (1992), DOI 10.1088/0959-7174/2/3/003 (Ch. 9, text
    above Listing 9.3, printed p. 168). The book names three other
    subharmonic routes: Herman and Strugala, Proc. SPIE 1221, pp. 183 to 192
    (1990), DOI 10.1117/12.18326; Johansson and Gavel, "Simulation of stellar
    speckle imaging", Proc. SPIE 2200, pp. 372 to 383 (1994),
    DOI 10.1117/12.177254, whose screens the book calls the closest match to
    theory (Ch. 9, text above Sec. 9.4, printed p. 172); and Sedmak, Appl.
    Opt. 37, pp. 4605 to 4613 (1998), DOI 10.1364/AO.37.004605. Frehlich,
    Appl. Opt. 39(3), pp. 393 to 397 (2000), DOI 10.1364/AO.39.000393, showed
    that subharmonic screens give the correct irradiance variance, and that
    plain Fourier screens do not.

    THE MEAN. Listing 9.3, line 38, printed p. 170, removes the mean of the
    low-frequency screen. This function does the same.

    VALIDITY. Add this screen to `ft_phase_screen`. Alone it holds only the
    frequencies below the grid fundamental 1/(n dx).
    '''
    rng = np.random.default_rng() if rng is None else rng
    n = int(n)
    side = n * float(dx_m)
    x = (np.arange(n) - n // 2) * float(dx_m)
    X, Y = np.meshgrid(x, x)
    lo = np.zeros((n, n), dtype=complex)

    for p in range(1, int(n_p) + 1):
        df_p = 1.0 / (3.0 ** p * side)              # Eq. (9.81), p. 169.
        fx = np.array([-1.0, 0.0, 1.0]) * df_p
        FX, FY = np.meshgrid(fx, fx)
        psd = phase_psd(np.hypot(FX, FY), r0_m, L0_m, l0_m)
        psd[1, 1] = 0.0                     # Listing 9.3, line 26, p. 170.
        cn = (rng.standard_normal((3, 3))
              + 1j * rng.standard_normal((3, 3))) * np.sqrt(psd) * df_p
        for c, fxi, fyi in zip(cn.ravel(), FX.ravel(), FY.ravel()):
            lo += c * np.exp(1j * 2.0 * np.pi * (fxi * X + fyi * Y))

    out = np.real(lo)
    return out - out.mean()                 # Listing 9.3, line 38, p. 170.


def ft_sh_phase_screen(r0_m, n, dx_m, L0_m=np.inf, l0_m=0.0, rng=None, n_p=3):
    '''
    Return one phase screen with the subharmonic low-frequency compensation.

    It is `ft_phase_screen` plus `subharmonic_screen`. See both for the
    parameters and the equations.

    Source: Schmidt (2010), DOI 10.1117/3.866274, Ch. 9, Listing 9.3, printed
    p. 170, which sums the high-frequency and the low-frequency screens.
    Figure 9.3, printed p. 169, shows the result against the theory of
    Eq. (9.44).

    VALIDITY. Ch. 9, text above Sec. 9.4, printed p. 172, states that the
    match to the theoretical structure function is close, not exact. The book
    points to Johansson and Gavel for a closer method.
    '''
    rng = np.random.default_rng() if rng is None else rng
    return (ft_phase_screen(r0_m, n, dx_m, L0_m, l0_m, rng)
            + subharmonic_screen(r0_m, n, dx_m, L0_m, l0_m, rng, n_p))


# ---------------------------------------------------------------------------
# 3. The per-screen Rytov bound (Sec. 9.2.5 and Listing 9.5)
# ---------------------------------------------------------------------------

def screen_rytov_share(r0_i_m, alpha_i, z_total_m, wavelength_m,
                       wave='spherical'):
    '''
    Return one screen's share of the path LOG-AMPLITUDE variance.

    Parameters:
        r0_i_m : array_like
            Fried parameter of each screen [m].
        alpha_i : array_like
            Fractional distance of each screen from the SOURCE, z_i / z.
        z_total_m : float
            Total path length [m].
        wavelength_m : float
            Wavelength [m].
        wave : str
            "spherical" for a point source, "plane" for a plane wave.

    Returns:
        numpy.ndarray
            The share of each screen. It is dimensionless.

    formula:
        plane:     d(sigma_chi^2)_i = 1.33 k^(-5/6) z^(5/6) r0_i^(-5/3)
                                      (1 - alpha_i)^(5/6)
        spherical: d(sigma_chi^2)_i = 1.33 k^(-5/6) z^(5/6) r0_i^(-5/3)
                                      alpha_i^(5/6) (1 - alpha_i)^(5/6)
    Source: Schmidt (2010), DOI 10.1117/3.866274, Ch. 9, Eq. (9.73), printed
    p. 165, is the plane-wave sum. Eq. (9.74), printed p. 165, is the
    spherical-wave sum. This function returns ONE TERM of that sum.

    WHAT THE CAP BOUNDS. Listing 9.5, lines 37 and 38, printed p. 175, set
    `rmax = 0.1` and then bound each screen through
    `x2 = rmax/1.33*(k/Dz)^(5/6) ./ A(2,:)`, where row 2 of A holds
    alpha^(5/6) (1 - alpha)^(5/6) and the unknown x holds r0_i^(-5/3). Multiply
    both sides by A(2,i) and the bound reads exactly

        1.33 k^(-5/6) z^(5/6) r0_i^(-5/3) alpha^(5/6) (1-alpha)^(5/6) <= 0.1

    which is this function with `wave="spherical"`. So the capped quantity is
    the screen's term in the SPHERICAL-WAVE LOG-AMPLITUDE variance
    sigma_chi^2, Eq. (9.64), printed p. 163. The book's text calls it "the
    overall Rytov number" (Sec. 9.5.1, printed p. 176), but the algebra above
    shows that it is sigma_chi^2, not the Rytov variance sigma_R^2. The two
    differ by a factor of 4: sigma_R^2 = 4 sigma_chi^2. The book credits the
    guideline to Martin and Flatte, Appl. Opt. 27(11), pp. 2111 to 2126
    (1988), DOI 10.1364/AO.27.002111.

    VALIDITY. Eq. (9.70), printed p. 165, states that r0_i is the PLANE-wave
    r0 of the layer, so the layer must be thin. The whole Rytov framework
    holds only for weak fluctuations, sigma_chi^2 < 0.25 (Ch. 9, text below
    Eq. (9.64), printed p. 163).
    '''
    k = 2.0 * np.pi / float(wavelength_m)
    a = np.asarray(alpha_i, dtype=float)
    x = np.asarray(r0_i_m, dtype=float) ** (-5.0 / 3.0)
    base = 1.33 * k ** (-5.0 / 6.0) * float(z_total_m) ** (5.0 / 6.0) * x
    if wave == 'plane':
        return base * (1.0 - a) ** (5.0 / 6.0)
    if wave == 'spherical':
        return base * a ** (5.0 / 6.0) * (1.0 - a) ** (5.0 / 6.0)
    raise ValueError("screen_rytov_share: wave must be 'plane' or "
                     "'spherical'; see Schmidt (2010), DOI 10.1117/3.866274, "
                     "Ch. 9, Eqs. (9.73) and (9.74), printed p. 165")


def max_screen_strength(alpha_i, z_total_m, wavelength_m, rmax=RMAX,
                        wave='spherical'):
    '''
    Return the largest allowed r0_i^(-5/3) of each screen.

    Parameters:
        alpha_i : array_like
            Fractional distance of each screen from the source, z_i / z.
        z_total_m : float
            Total path length [m].
        wavelength_m : float
            Wavelength [m].
        rmax : float
            The cap on one screen's log-amplitude share. The book uses 0.1.
        wave : str
            "spherical" or "plane".

    Returns:
        numpy.ndarray
            The upper bound on r0_i^(-5/3) [m^(-5/3)]. It is infinite where
            the path weight is zero, that is at alpha = 0 and alpha = 1 for a
            spherical wave.

    formula:
        x2_i = rmax / [1.33 k^(-5/6) z^(5/6) w(alpha_i)]
    Source: Schmidt (2010), DOI 10.1117/3.866274, Ch. 9, Listing 9.5, lines 37
    to 39, printed p. 175. It is `screen_rytov_share` inverted for
    r0_i^(-5/3).

    VALIDITY. Listing 9.5, line 39, printed p. 175, replaces the infinite
    bound with the finite value 50^(-5/3), that is a screen r0 of 50 m. That
    is a numerical convenience of `fmincon`, not physics. This function
    returns the true infinity. The caller replaces it.
    '''
    unit = screen_rytov_share(np.ones_like(np.asarray(alpha_i, dtype=float)),
                              alpha_i, z_total_m, wavelength_m, wave)
    return np.where(unit > 0.0, float(rmax) / np.where(unit > 0.0, unit, 1.0),
                    np.inf)


def screen_strengths(alpha_i, r0_sw_m, sigma2_chi_sw, z_total_m, wavelength_m,
                     rmax=RMAX, r0_max_m=50.0):
    '''
    Solve for the screen Fried parameters at the given screen positions.

    Parameters:
        alpha_i : array_like
            Fractional distance of each screen from the source, z_i / z. The
            book spaces them uniformly, from 0 to 1.
        r0_sw_m : float
            The target spherical-wave Fried parameter of the whole path [m].
        sigma2_chi_sw : float
            The target spherical-wave log-amplitude variance of the path.
        z_total_m : float
            Total path length [m].
        wavelength_m : float
            Wavelength [m].
        rmax : float
            The per-screen log-amplitude cap. See `screen_rytov_share`.
        r0_max_m : float
            The r0 that stands for an infinitely weak screen [m].

    Returns:
        numpy.ndarray
            The Fried parameter of each screen [m].

    formula:
        A x = b,   x_i = r0_i^(-5/3)
        A[0,i] = alpha_i^(5/3)
        A[1,i] = alpha_i^(5/6) (1 - alpha_i)^(5/6)
        b = [ r0_sw^(-5/3),  sigma_chi,sw^2 / 1.33 * (k/z)^(5/6) ]
        0 <= x_i <= rmax / 1.33 * (k/z)^(5/6) / A[1,i]
    Source: Schmidt (2010), DOI 10.1117/3.866274, Ch. 9, Eq. (9.75), printed
    p. 165, gives the 5-screen matrix system, and the text below it, printed
    p. 166, names the two rows. Listing 9.5, lines 26 to 41, printed p. 175,
    prints the constrained least-squares solve.

    THE SOLVER. Listing 9.5 minimises `sum((A*X - b).^2)` with `fmincon` under
    a box constraint. This function calls `scipy.optimize.lsq_linear`, which
    solves the same bounded linear least-squares problem. The book's MATLAB
    listing is not ported.

    VALIDITY.
    - The system is UNDERDETERMINED: two equations for n unknowns (Ch. 9, text
      below Eq. (9.74), printed p. 165). It matches r0 and sigma_chi^2 only.
      It does NOT match the moments of Eq. (9.65). See `moment_error`.
    - A negative x_i is unphysical, so the lower bound is zero (Ch. 9, text
      below Eq. (9.75), printed p. 166). A zero x_i means an infinitely weak
      screen; this function returns `r0_max_m` for it.
    - The book fixes the screen positions first, because the free positions
      make the system far more underdetermined (Ch. 9, text below Eq. (9.74),
      printed p. 165).
    '''
    a = np.asarray(alpha_i, dtype=float)
    k = 2.0 * np.pi / float(wavelength_m)
    scale = (k / float(z_total_m)) ** (5.0 / 6.0)

    A = np.vstack([a ** (5.0 / 3.0),
                   a ** (5.0 / 6.0) * (1.0 - a) ** (5.0 / 6.0)])
    b = np.array([float(r0_sw_m) ** (-5.0 / 3.0),
                  float(sigma2_chi_sw) / 1.33 * scale])

    upper = max_screen_strength(a, z_total_m, wavelength_m, rmax, 'spherical')
    upper = np.where(np.isfinite(upper), upper,
                     float(r0_max_m) ** (-5.0 / 3.0))
    x = lsq_linear(A, b, bounds=(np.zeros_like(a), upper)).x
    floor = float(r0_max_m) ** (-5.0 / 3.0)
    return np.maximum(x, floor) ** (-3.0 / 5.0)


def composite_r0(r0_i_m, alpha_i=None, wave='plane'):
    '''
    Add the screen Fried parameters into the path Fried parameter.

    Parameters:
        r0_i_m : array_like
            Fried parameter of each screen [m].
        alpha_i : array_like or None
            Fractional distance of each screen from the source. It is needed
            for the spherical case only.
        wave : str
            "plane" or "spherical".

    Returns:
        float
            The path Fried parameter [m].

    formula:
        plane:     r0 = [ SUM r0_i^(-5/3) ]^(-3/5)
        spherical: r0 = [ SUM r0_i^(-5/3) alpha_i^(5/3) ]^(-3/5)
    Source: Schmidt (2010), DOI 10.1117/3.866274, Ch. 9, Eqs. (9.71) and
    (9.72), printed p. 165.

    VALIDITY. The screens must be independent. Each r0_i is the plane-wave
    value of Eq. (9.70), printed p. 165, so each layer must be thin.
    '''
    x = np.asarray(r0_i_m, dtype=float) ** (-5.0 / 3.0)
    if wave == 'spherical':
        x = x * np.asarray(alpha_i, dtype=float) ** (5.0 / 3.0)
    elif wave != 'plane':
        raise ValueError("composite_r0: wave must be 'plane' or 'spherical'; "
                         "see Schmidt (2010), DOI 10.1117/3.866274, Ch. 9, "
                         "Eqs. (9.71) and (9.72), printed p. 165")
    return float(np.sum(x) ** (-3.0 / 5.0))


def screen_r0(cn2_integral_m13, wavelength_m):
    '''
    Return the Fried parameter of one thin layer.

    Parameters:
        cn2_integral_m13 : array_like
            The integral of Cn2 over the layer [m^(1/3)].
        wavelength_m : float
            Wavelength [m].

    Returns:
        The Fried parameter [m]. The type follows the input.

    formula:
        r0_i = (0.423 k^2 Cn2_i dz_i)^(-3/5)
    Source: Schmidt (2010), DOI 10.1117/3.866274, Ch. 9, Eq. (9.70), printed
    p. 165. The book credits it to Roggemann, Welsh, Montera and Rhoadarmer,
    Appl. Opt. 34(20), pp. 4037 to 4051 (1995), DOI 10.1364/AO.34.004037.

    VALIDITY. Ch. 9, text below Eq. (9.70), printed p. 165, states that this
    is the PLANE-wave r0, so it is valid only when the layer is very thin. A
    phase screen is thin when its thickness is much less than the propagation
    distance after it (Ch. 9, Sec. 9.2.4, printed p. 164).
    '''
    k = 2.0 * np.pi / float(wavelength_m)
    return (0.423 * k * k
            * np.asarray(cn2_integral_m13, dtype=float)) ** (-3.0 / 5.0)


# ---------------------------------------------------------------------------
# 4. The layered-atmosphere moment rule (Sec. 9.2.5, Eq. (9.65))
# ---------------------------------------------------------------------------

def profile_moments(cn2, z_m, m_max=7):
    '''
    Return the low-order moments of a CONTINUOUS Cn2 profile.

    Parameters:
        cn2 : array_like
            Cn2 sampled on `z_m` [m^(-2/3)].
        z_m : array_like
            The distance grid [m]. It must go up.
        m_max : int
            The highest moment. The book uses 7.

    Returns:
        numpy.ndarray
            The moments for m = 0 to m_max. Moment m has the units
            m^(1/3 + m).

    formula:
        mu_m = INTEGRAL_0^z Cn2(z) z^m dz,   0 <= m <= 7
    Source: Schmidt (2010), DOI 10.1117/3.866274, Ch. 9, Eq. (9.65), printed
    p. 164, left side.

    THE QUADRATURE. This function uses the trapezium rule on the given grid.
    The book gives no quadrature. A coarse grid gives a coarse moment, and the
    high moments feel it most, because z^7 weights the far end of the path.

    VALIDITY. Ch. 9, text above Eq. (9.65), printed p. 164, states the purpose:
    match the moments and the layered model then reproduces r0, theta_0 and
    sigma_chi^2 of the bulk turbulence.
    '''
    cn2 = np.asarray(cn2, dtype=float)
    z = np.asarray(z_m, dtype=float)
    return np.array([float(_trapz(cn2 * z ** m, z))
                     for m in range(int(m_max) + 1)])


def layer_moments(cn2_integral_m13, z_m, m_max=7):
    '''
    Return the low-order moments of a LAYERED (screen) model.

    Parameters:
        cn2_integral_m13 : array_like
            The integral of Cn2 over each layer, Cn2_i dz_i [m^(1/3)].
        z_m : array_like
            The distance of each screen from the source [m].
        m_max : int
            The highest moment. The book uses 7.

    Returns:
        numpy.ndarray
            The moments for m = 0 to m_max.

    formula:
        mu_m = SUM_i Cn2_i z_i^m dz_i,   0 <= m <= 7
    Source: Schmidt (2010), DOI 10.1117/3.866274, Ch. 9, Eq. (9.65), printed
    p. 164, right side.

    VALIDITY. The rule treats each layer as a POINT at z_i that carries the
    whole slab integral. It says nothing about the thickness of the slab, so a
    thick slab must still obey the thin-screen rule of Eq. (9.70).
    '''
    w = np.asarray(cn2_integral_m13, dtype=float)
    z = np.asarray(z_m, dtype=float)
    return np.array([float(np.sum(w * z ** m)) for m in range(int(m_max) + 1)])


def moment_error(cn2, z_profile_m, cn2_integral_m13, z_screens_m, m_max=7):
    '''
    Return the relative moment error of a proposed layering.

    Parameters:
        cn2 : array_like
            The continuous Cn2 profile [m^(-2/3)].
        z_profile_m : array_like
            The distance grid of that profile [m].
        cn2_integral_m13 : array_like
            The integral of Cn2 over each proposed layer [m^(1/3)].
        z_screens_m : array_like
            The distance of each proposed screen from the source [m].
        m_max : int
            The highest moment. The book uses 7.

    Returns:
        numpy.ndarray
            (layered - continuous) / continuous, for m = 0 to m_max. The value
            is 0.0 where the continuous moment is zero.

    Source: Schmidt (2010), DOI 10.1117/3.866274, Ch. 9, Eq. (9.65), printed
    p. 164. This function measures how well a layering satisfies that
    equality. Sec. 9.5.5, printed p. 182, names the same equation as the fix
    for a simulation that does not verify.

    THE BOOK GIVES NO TOLERANCE. It states the equality and it gives no error
    budget. The caller sets the tolerance.

    VALIDITY. A layering with n screens has 2n free numbers (the position and
    the strength of each screen). Eq. (9.65) gives 8 equations for
    0 <= m <= 7, so 4 screens are the smallest set that CAN match all 8
    moments. The book does not say this; it is the count.
    '''
    want = profile_moments(cn2, z_profile_m, m_max)
    got = layer_moments(cn2_integral_m13, z_screens_m, m_max)
    safe = np.where(want != 0.0, want, 1.0)
    return np.where(want != 0.0, (got - want) / safe, 0.0)


# ---------------------------------------------------------------------------
# 5. The turbulent sampling bounds (Sec. 9.4)
# ---------------------------------------------------------------------------

def fresnel_pitch_max(wavelength_m, z_m):
    '''
    Return the grid pitch that samples the SCINTILLATION.

    Parameters:
        wavelength_m : float
            Wavelength [m].
        z_m : array_like
            The distance from the screen to the observation plane [m].

    Returns:
        The largest allowed pitch [m]. The type follows the input.

    formula:
        dx <= sqrt(lambda z) / 2
    Source: Schmidt (2010), DOI 10.1117/3.866274, Ch. 9, Sec. 9.4, printed
    p. 172. The book states: "The scale size of scintillation is given
    approximately by the Fresnel length (lambda z)^(1/2), so they set delta_i
    to be the smallest of [the phase pitch], (lambda z)^(1/2)/2, and the grid
    spacing that just barely avoids aliasing of the free-space point spread
    function." The book credits the rule to Johnston and Lane, Appl. Opt.
    39(26), pp. 4761 to 4769 (2000), DOI 10.1364/AO.39.004761.

    VALIDITY. The bound is a prose rule of thumb, not a derived inequality.
    The book gives no equation number for it. Two samples per Fresnel length
    is the bare Nyquist rate of the irradiance correlation scale.
    '''
    return np.sqrt(float(wavelength_m) * np.asarray(z_m, dtype=float)) / 2.0


def phase_pitch_max(r0_m, max_step_rad=np.pi, sigmas=3.0):
    '''
    Return the grid pitch that samples the turbulent PHASE.

    Parameters:
        r0_m : float
            The Fried parameter of the screen or the path [m].
        max_step_rad : float
            The largest phase step between two adjacent samples [rad].
        sigmas : float
            How many standard deviations the criterion holds to. 3.0 gives
            99.7%.

    Returns:
        float
            The largest allowed pitch [m].

    formula:
        sigmas * sqrt( 6.88 (dx/r0)^(5/3) ) <= max_step_rad
        dx <= r0 [ max_step_rad^2 / (sigmas^2 * 6.88) ]^(3/5)
    Source: Schmidt (2010), DOI 10.1117/3.866274, Ch. 9, Sec. 9.4, printed
    p. 172, states the criterion in prose: "they compute the grid spacing
    delta at which phase differences less than pi in adjacent grid points
    occur more than 99.7% of the time". The structure function is Ch. 9,
    Eq. (9.44), printed p. 160. The book credits the rule to Johnston and
    Lane, Appl. Opt. 39(26), pp. 4761 to 4769 (2000),
    DOI 10.1364/AO.39.004761.

    THE ALGEBRA IS OURS. The book gives NO equation for this rule. The phase
    difference between two samples is Gaussian with a variance D(dx), so
    99.7% of the draws stay inside 3 sqrt(D). Set 3 sqrt(D(dx)) = pi and solve
    for dx. The default arguments give dx <= 0.332 r0, that is 3.0 pixels per
    r0.

    VALIDITY. The Gaussian assumption is the central-limit argument of Ch. 9,
    text below Eq. (9.78), printed p. 167. The structure function is the
    Kolmogorov one, so l0 = 0 and L0 = infinite.
    '''
    ratio = float(max_step_rad) ** 2 / (float(sigmas) ** 2 * 6.88)
    return float(r0_m) * ratio ** (3.0 / 5.0)


def blurred_extent(d_m, wavelength_m, dz_m, r0_m, c=2.0):
    '''
    Return the turbulence-blurred aperture extent.

    Parameters:
        d_m : float
            The vacuum extent: the source size D1, or the observation region
            of interest D2 [m].
        wavelength_m : float
            Wavelength [m].
        dz_m : float
            The propagation distance [m].
        r0_m : float
            The Fried parameter for the direction of travel [m]. Use the
            REVERSE-path r0 for D1.
        c : float
            The blur sensitivity. The book allows 2 to 8.

    Returns:
        float
            The blurred extent [m].

    formula:
        D1' = D1 + c lambda dz / r0_rev
        D2' = D2 + c lambda dz / r0
    Source: Schmidt (2010), DOI 10.1117/3.866274, Ch. 9, Eqs. (9.84) and
    (9.85), printed p. 173. The book credits the model to Coy; see also
    Mansell and Praus (Ch. 9, text above Eq. (9.84), printed p. 173).

    THE CONSTANT c. Ch. 9, text below Eq. (9.85), printed p. 173: "Typical
    values of c range from 2 to 8. Choosing c = 2 typically captures 97% of
    the light, and choosing c = 4 typically captures 99% of the light."
    Listing 9.6, line 2, printed p. 177, uses c = 2.

    VALIDITY. The model treats the turbulent spread as a diffraction grating
    of period r0. It is a sampling aid, not a beam-spread model.
    '''
    blur = float(c) * float(wavelength_m) * float(dz_m) / float(r0_m)
    return float(d_m) + blur


def constraint1_pitch_max(pitch_in_m, d1_m, d2_m, wavelength_m, dz_m):
    '''
    Return the largest observation-plane pitch that CONSTRAINT 1 allows.

    Parameters:
        pitch_in_m : float
            The source-plane pitch Delta1 [m].
        d1_m : float
            The blurred source extent D1' [m].
        d2_m : float
            The blurred observation extent D2' [m].
        wavelength_m : float
            Wavelength [m].
        dz_m : float
            The propagation distance [m].

    Returns:
        float
            The largest allowed Delta_n [m].

    formula:
        Delta_n <= -(D2'/D1') Delta1 + lambda dz / D1'
    Source: Schmidt (2010), DOI 10.1117/3.866274, Ch. 9, Eq. (9.86), printed
    p. 173. It is the vacuum Constraint 1, Ch. 7, Eq. (7.14), printed p. 119,
    with the blurred extents of Eqs. (9.84) and (9.85).

    VALIDITY. Constraint 1 makes the source grid fine enough to hold every ray
    that lands inside the observation-plane region of interest (Ch. 9, text
    below Eq. (9.83), printed p. 173).
    '''
    return (-(float(d2_m) / float(d1_m)) * float(pitch_in_m)
            + float(wavelength_m) * float(dz_m) / float(d1_m))


def constraint2_n_min(pitch_in_m, pitch_out_m, d1_m, d2_m, wavelength_m,
                      dz_m):
    '''
    Return the smallest grid count that CONSTRAINT 2 allows.

    Parameters:
        pitch_in_m : float
            The source-plane pitch Delta1 [m].
        pitch_out_m : float
            The observation-plane pitch Delta_n [m].
        d1_m : float
            The blurred source extent D1' [m].
        d2_m : float
            The blurred observation extent D2' [m].
        wavelength_m : float
            Wavelength [m].
        dz_m : float
            The propagation distance [m].

    Returns:
        float
            The smallest allowed N.

    formula:
        N >= D1'/(2 Delta1) + D2'/(2 Delta_n) + lambda dz / (2 Delta1 Delta_n)
    Source: Schmidt (2010), DOI 10.1117/3.866274, Ch. 9, Eq. (9.87), printed
    p. 174. It is the vacuum Constraint 2, Ch. 7, Eq. (7.20), printed p. 120,
    with the blurred extents.

    VALIDITY. Constraint 2 makes the grid wide enough that the wrap-around of
    the periodic transform stays outside the region of interest.
    '''
    d1p, d2p = float(d1_m), float(d2_m)
    p1, pn = float(pitch_in_m), float(pitch_out_m)
    return (d1p / (2.0 * p1) + d2p / (2.0 * pn)
            + float(wavelength_m) * float(dz_m) / (2.0 * p1 * pn))


def constraint3_pitch_range(pitch_in_m, d1_m, wavelength_m, dz_m,
                            curvature_m=np.inf):
    '''
    Return the observation-plane pitch band that CONSTRAINT 3 allows.

    Parameters:
        pitch_in_m : float
            The source-plane pitch Delta1 [m].
        d1_m : float
            The blurred source extent D1' [m].
        wavelength_m : float
            Wavelength [m].
        dz_m : float
            The propagation distance [m].
        curvature_m : float
            The source wavefront radius R [m]. Infinite is collimated. A
            negative value diverges; see Ch. 7, Eq. (7.32), printed p. 122.

    Returns:
        tuple
            (lowest, highest) allowed Delta_n [m].

    formula:
        (1 + dz/R) Delta1 - lambda dz / D1'
            <= Delta_n <=
        (1 + dz/R) Delta1 + lambda dz / D1'
    Source: Schmidt (2010), DOI 10.1117/3.866274, Ch. 9, Eq. (9.88), printed
    p. 174. It is the vacuum Constraint 3, Ch. 7, Eq. (7.53), printed p. 126,
    with the blurred source extent.

    VALIDITY. Constraint 3 stops the aliasing of the quadratic phase factor of
    the angular-spectrum method (Ch. 7, Eqs. (7.46) and (7.52), printed
    pp. 125 and 126). Ch. 7, Eq. (7.60), printed p. 129, states that the
    constraint does not apply when 1 + dz/R < D2/D1, because the geometric
    beam then stays inside D2. This function does NOT apply that exemption.
    '''
    mag = 1.0 + float(dz_m) / float(curvature_m)
    band = float(wavelength_m) * float(dz_m) / float(d1_m)
    centre = mag * float(pitch_in_m)
    return centre - band, centre + band


def max_partial_step(pitch_in_m, pitch_out_m, n, wavelength_m):
    '''
    Return the longest partial propagation that the grid samples.

    Parameters:
        pitch_in_m : float
            The source-plane pitch Delta1 [m].
        pitch_out_m : float
            The observation-plane pitch Delta_n [m].
        n : int
            Grid points per side.
        wavelength_m : float
            Wavelength [m].

    Returns:
        float
            The longest step dz_max [m].

    formula:
        dz_max = min(Delta1, Delta_n)^2 N / lambda
    Source: Schmidt (2010), DOI 10.1117/3.866274, Ch. 9, Eq. (9.89), printed
    p. 174. It repeats Ch. 8, Eq. (8.24), printed p. 144.

    VALIDITY. It is Constraint 4 (Ch. 7, Eq. (7.59), printed p. 127) rearranged
    for one partial step. Ch. 9, text above Eq. (9.86), printed p. 173, states
    that turbulence does NOT change Constraint 4, because it is a rule of the
    numerical method, not of the geometry.
    '''
    return (min(float(pitch_in_m), float(pitch_out_m)) ** 2 * int(n)
            / float(wavelength_m))


def min_planes(dz_m, dz_max_m):
    '''
    Return the smallest number of partial-propagation planes.

    Parameters:
        dz_m : float
            The total propagation distance [m].
        dz_max_m : float
            The longest partial step [m]. Use `max_partial_step`.

    Returns:
        int
            The plane count, source and observation planes included.

    formula:
        n_min = ceil(dz / dz_max) + 1
    Source: Schmidt (2010), DOI 10.1117/3.866274, Ch. 9, Eq. (9.90), printed
    p. 174. It repeats the plane count of Ch. 8, text below Eq. (8.24),
    printed p. 144.

    VALIDITY. It is a SAMPLING floor only. Ch. 9, Sec. 9.5.2, printed p. 177,
    shows the gap: the 50 km example needs 2 planes by this rule, and the book
    then uses 11 planes "to represent the atmosphere properly". The book gives
    NO formula for that second, larger floor.
    '''
    return int(np.ceil(float(dz_m) / float(dz_max_m)) + 1)


# ---------------------------------------------------------------------------
# 6. The Section 9.5 procedure, as a checklist
# ---------------------------------------------------------------------------

def properly_sampled_checklist(*, wavelength_m, z_total_m, n, pixel_m,
                               z_m, r0_m, r0_total_m, d_tx_m, d_rx_m,
                               size_m=None, pixel_out_m=None,
                               curvature_m=np.inf, c=2.0, rmax=RMAX,
                               wave='spherical'):
    '''
    Check one simulation plan against the Section 9.5 procedure.

    Every argument is a PLAIN number or a plain array. This function imports
    no olb module. The names match the olb wave-optics layer one to one:
    `n`, `pixel_m` and `size_m` are `GridSpec` fields, and `z_m`, `r0_m`,
    `r0_total_m` and `z_total_m` are `ScreenPlan` fields.

    Parameters:
        wavelength_m : float
            Wavelength [m].
        z_total_m : float
            The length of the turbulent path [m].
        n : int
            Grid points per side.
        pixel_m : float
            The source-plane grid pitch Delta1 [m].
        z_m : array_like
            The distance of each screen from the SOURCE [m].
        r0_m : array_like
            The Fried parameter of each screen [m].
        r0_total_m : float
            The Fried parameter of the whole path [m].
        d_tx_m : float
            The source extent D1 [m].
        d_rx_m : float
            The observation region of interest D2 [m].
        size_m : float or None
            The grid side [m]. None takes n * pixel_m.
        pixel_out_m : float or None
            The observation-plane pitch Delta_n [m]. None takes pixel_m, which
            is the flat-grid case.
        curvature_m : float
            The source wavefront radius R [m].
        c : float
            The blur sensitivity of Eqs. (9.84) and (9.85).
        rmax : float
            The per-screen log-amplitude cap of Listing 9.5.
        wave : str
            "spherical" or "plane". It selects the path weight.

    Returns:
        list
            One tuple (rule, satisfied, bound, actual, citation) for each step.
            `satisfied` is True, False, or None. None marks an ADVISORY step
            that this function cannot test.

    Source: Schmidt (2010), DOI 10.1117/3.866274, Ch. 9, Sec. 9.5, printed
    pp. 174 to 182. The five subsections are 9.5.1 (geometry and turbulence),
    9.5.2 (the sampling constraints), 9.5.3 (a vacuum run), 9.5.4 (the
    turbulent runs), and 9.5.5 (verify the output).

    VALIDITY.
    - The function tests ONE pitch pair and ONE grid count. It does not draw
      the constraint plot of Fig. 9.6, printed p. 178.
    - Steps 9.5.3, 9.5.4 and 9.5.5 are procedures, not inequalities, so they
      come back as advisory rows with `satisfied = None`.
    - Constraint 3 is not exempted here. See `constraint3_pitch_range`.
    - The scintillation pitch rule skips a screen whose path weight is zero,
      because such a screen adds no scintillation. For a spherical wave those
      are the screens at alpha = 0 and alpha = 1. The book does not state this
      exemption; it follows from Eq. (9.74), printed p. 165.
    '''
    z = np.asarray(z_m, dtype=float)
    r0 = np.asarray(r0_m, dtype=float)
    size_m = float(n) * float(pixel_m) if size_m is None else float(size_m)
    pn = float(pixel_m) if pixel_out_m is None else float(pixel_out_m)
    alpha = z / float(z_total_m)
    out = []

    # ---- 9.5.1 the turbulence conditions ----
    share = screen_rytov_share(r0, alpha, z_total_m, wavelength_m, wave)
    sigma2_chi = float(np.sum(share))
    out.append((
        '9.5.1 weak fluctuation', sigma2_chi < WEAK_SIGMA2_CHI,
        WEAK_SIGMA2_CHI, sigma2_chi,
        'Schmidt (2010), DOI 10.1117/3.866274, Ch. 9, Eq. (9.64) and the text '
        'below it, printed p. 163'))
    out.append((
        '9.5.1 per-screen Rytov cap', float(np.max(share)) <= float(rmax),
        float(rmax), float(np.max(share)),
        'Schmidt (2010), DOI 10.1117/3.866274, Ch. 9, Listing 9.5, lines 37 '
        'and 38, printed p. 175'))

    # ---- 9.5.2 the sampling constraints ----
    d1p = blurred_extent(d_tx_m, wavelength_m, z_total_m, r0_total_m, c)
    d2p = blurred_extent(d_rx_m, wavelength_m, z_total_m, r0_total_m, c)
    c1 = constraint1_pitch_max(pixel_m, d1p, d2p, wavelength_m, z_total_m)
    out.append((
        '9.5.2 constraint 1', pn <= c1, c1, pn,
        'Schmidt (2010), DOI 10.1117/3.866274, Ch. 9, Eq. (9.86), printed '
        'p. 173'))

    n2 = constraint2_n_min(pixel_m, pn, d1p, d2p, wavelength_m, z_total_m)
    out.append((
        '9.5.2 constraint 2', float(n) >= n2, n2, float(n),
        'Schmidt (2010), DOI 10.1117/3.866274, Ch. 9, Eq. (9.87), printed '
        'p. 174'))

    lo, hi = constraint3_pitch_range(pixel_m, d1p, wavelength_m, z_total_m,
                                     curvature_m)
    out.append((
        '9.5.2 constraint 3', lo <= pn <= hi, (lo, hi), pn,
        'Schmidt (2010), DOI 10.1117/3.866274, Ch. 9, Eq. (9.88), printed '
        'p. 174'))

    dz_max = max_partial_step(pixel_m, pn, n, wavelength_m)
    planes_min = min_planes(z_total_m, dz_max)
    planes = int(z.size) + 2                 # the source and the observation
    out.append((
        '9.5.2 partial-propagation planes', planes >= planes_min, planes_min,
        planes,
        'Schmidt (2010), DOI 10.1117/3.866274, Ch. 9, Eqs. (9.89) and (9.90), '
        'printed p. 174'))

    # ---- 9.4 the two turbulent pitch rules ----
    # A screen with a zero path weight adds no scintillation, so it is exempt.
    # For a spherical wave those are the screens at alpha = 0 and alpha = 1.
    live = share > 0.0
    dx_fresnel = float(np.min(fresnel_pitch_max(
        wavelength_m, np.maximum(float(z_total_m) - z[live], 0.0)))
        if np.any(live) else np.inf)
    out.append((
        '9.4 scintillation pitch', float(pixel_m) <= dx_fresnel, dx_fresnel,
        float(pixel_m),
        'Schmidt (2010), DOI 10.1117/3.866274, Ch. 9, Sec. 9.4, printed '
        'p. 172 (Johnston and Lane, DOI 10.1364/AO.39.004761)'))

    dx_phase = phase_pitch_max(r0_total_m)
    out.append((
        '9.4 phase pitch', float(pixel_m) <= dx_phase, dx_phase,
        float(pixel_m),
        'Schmidt (2010), DOI 10.1117/3.866274, Ch. 9, Sec. 9.4, printed '
        'p. 172 (Johnston and Lane, DOI 10.1364/AO.39.004761)'))

    # ---- the advisory steps ----
    out.append((
        '9.5.3 vacuum run', None, None, None,
        'Schmidt (2010), DOI 10.1117/3.866274, Ch. 9, Sec. 9.5.3, printed '
        'p. 178: run the same grid with no screens and compare it to the '
        'analytic field.'))
    out.append((
        '9.5.4 independent realisations', None, None, None,
        'Schmidt (2010), DOI 10.1117/3.866274, Ch. 9, Sec. 9.5.4, printed '
        'p. 179: draw a new set of screens for each realisation. Move the '
        'screens for a temporal run.'))
    out.append((
        '9.5.5 verify the output', None, None, None,
        'Schmidt (2010), DOI 10.1117/3.866274, Ch. 9, Sec. 9.5.5, printed '
        'p. 180: check the screen structure function against Eq. (9.44), and '
        'the observation-plane coherence factor against Eqs. (9.32) and '
        '(9.44).'))
    return out


if __name__ == '__main__':
    import time

    from .fourier import structure_function

    t_start = time.time()

    # ---------------- 1. the PSD reductions ----------------
    f = np.logspace(-3, 3, 400)
    r0 = 0.10

    # The modified von Karman reduces to Kolmogorov for L0 = inf and l0 = 0.
    # Ch. 9, text below Eq. (9.18), printed p. 155.
    mvk = modified_von_karman_phase_psd(f, r0, np.inf, 0.0)
    kol = kolmogorov_phase_psd(f, r0)
    err_red = float(np.max(np.abs(mvk / kol - 1.0)))
    assert err_red < 1e-14, err_red

    # The bare Kolmogorov form, Eq. (9.52), printed p. 161.
    hand = 0.023 * r0 ** (-5.0 / 3.0) * f ** (-11.0 / 3.0)
    err_k = float(np.max(np.abs(kol / hand - 1.0)))
    assert err_k < 1e-14, err_k

    # The angular constant 0.49 of Eq. (9.49) and the ordinary constant 0.023
    # of Eq. (9.52) agree: 0.49 (2 pi)^(-5/3) = 0.023.
    conv = 0.49 * (2.0 * np.pi) ** (-5.0 / 3.0)
    assert abs(conv / 0.023 - 1.0) < 0.005, conv

    # The outer scale flattens the low frequencies, and the inner scale cuts
    # the high ones. Eqs. (9.50) and (9.51), printed p. 161.
    vk = von_karman_phase_psd(f, r0, 10.0)
    assert vk[0] < kol[0], (vk[0], kol[0])            # the low end is capped
    assert abs(vk[-1] / kol[-1] - 1.0) < 1e-6         # the high end is free
    mvk2 = modified_von_karman_phase_psd(f, r0, 10.0, 0.01)
    assert mvk2[-1] < 1e-6 * vk[-1], (mvk2[-1], vk[-1])
    print(f'REDUCTION modified von Karman (L0=inf, l0=0) vs Kolmogorov : '
          f'max rel err = {err_red:.3e}  (target 1e-14)')
    print(f'REDUCTION Kolmogorov PSD vs 0.023 r0^(-5/3) f^(-11/3) : '
          f'max rel err = {err_k:.3e}  (target 1e-14)')
    print(f'REDUCTION 0.49 (2 pi)^(-5/3) = {conv:.5f} vs the printed 0.023 : '
          f'rel err = {abs(conv / 0.023 - 1.0):.3e}  (target 5e-3)')
    print('')

    # ---------------- 2. the screen structure function ----------------
    # Ch. 9, Eq. (9.44), printed p. 160, and Fig. 9.3, printed p. 169.
    # The tolerance philosophy follows the self-check of
    # olb/waveoptics/turbulence/screens.py: assert a BAND, not a point.
    t0 = time.time()
    N, dx, r0s, M = 512, 0.01, 0.10, 24        # a 5.12 m side, that is 51 r0

    axis = (np.arange(N) - N // 2) * dx
    XX, YY = np.meshgrid(axis, axis)
    # A pupil of 1.2 m radius. It leaves a guard band, so the circular
    # correlation of `structure_function` stays valid out to 2.4 m = 24 r0.
    pupil = (np.hypot(XX, YY) <= 1.2).astype(float)

    acc_sh = np.zeros((N, N))
    acc_ft = np.zeros((N, N))
    for i in range(M):
        # The two screens share the seed, so the FT part is the SAME draw.
        # The only difference is the subharmonic addition.
        acc_ft += structure_function(
            ft_phase_screen(r0s, N, dx, rng=np.random.default_rng(5000 + i)),
            pupil, dx)
        acc_sh += structure_function(
            ft_sh_phase_screen(r0s, N, dx,
                               rng=np.random.default_rng(5000 + i)),
            pupil, dx)
    acc_sh /= M
    acc_ft /= M

    # Read the x axis of the mean structure function.
    row_sh = acc_sh[N // 2, :]
    row_ft = acc_ft[N // 2, :]
    probes = np.array([0.3, 0.5, 0.8, 1.2, 1.6, 3.2, 8.0])   # in units of r0
    idx = N // 2 + np.round(probes * r0s / dx).astype(int)
    theory = kolmogorov_structure_function((idx - N // 2) * dx, r0s)
    ratio_sh = row_sh[idx] / theory
    ratio_ft = row_ft[idx] / theory
    t_sf = time.time() - t0

    # THE STATED BAND is r/r0 = 0.3 to 1.6. It is the band that the self-check
    # of olb/waveoptics/turbulence/screens.py measures, and the tolerance
    # (0.85, 1.02) is that file's tolerance. The subharmonic screen of this
    # module lands inside it.
    band = probes <= 1.6
    assert np.all(ratio_sh[band] > 0.85), ratio_sh[band]
    assert np.all(ratio_sh[band] < 1.02), ratio_sh[band]
    # The plain Fourier screen is LOW everywhere, and the deficit grows with
    # the separation. Ch. 9, text below Listing 9.2, printed p. 167.
    assert np.all(ratio_ft < 0.85), ratio_ft
    assert np.all(ratio_ft < ratio_sh), (ratio_ft, ratio_sh)
    assert ratio_ft[-1] < 0.55, ratio_ft[-1]
    assert ratio_ft[-1] < ratio_ft[0], ratio_ft
    # The subharmonics do NOT close the gap at a large separation. Say so.
    assert ratio_sh[-1] < 0.85, ratio_sh[-1]

    print(f'structure function, N={N}, r0 = {r0s / dx:.0f} px, {M} screens, '
          f'pupil 1.2 m:')
    print(f"  {'r/r0':>6}{'D theory':>12}{'D subharm':>12}{'ratio':>8}"
          f"{'D FT only':>12}{'ratio':>8}")
    for p_, th, ds, rs, dfo, rf in zip(probes, theory, row_sh[idx], ratio_sh,
                                       row_ft[idx], ratio_ft):
        print(f'  {p_:>6.1f}{th:>12.2f}{ds:>12.2f}{rs:>8.3f}'
              f'{dfo:>12.2f}{rf:>8.3f}')
    print('  (the stated band is r/r0 = 0.3 to 1.6, tolerance 0.85 '
          'to 1.02)')
    print(f'  (elapsed {t_sf:.1f} s)')
    print('')

    # ---------------- 3. the moment rule, Eq. (9.65) ----------------
    # A uniform Cn2 slab. The continuous moments are C z^(m+1)/(m+1).
    Z, C = 50e3, 1e-16
    zc = np.linspace(0.0, Z, 20001)
    cn2c = np.full_like(zc, C)
    mu = profile_moments(cn2c, zc, 7)
    mu_hand = np.array([C * Z ** (m + 1) / (m + 1) for m in range(8)])
    err_mu = float(np.max(np.abs(mu / mu_hand - 1.0)))
    assert err_mu < 1e-6, err_mu

    # 4-point Gauss-Legendre is EXACT for a polynomial of degree 7, so 4
    # screens at the Gauss nodes match every moment 0 <= m <= 7 of a uniform
    # slab. That is the smallest exact case of Eq. (9.65), printed p. 164.
    nodes, wts = np.polynomial.legendre.leggauss(4)
    z_gl = 0.5 * Z * (nodes + 1.0)
    w_gl = 0.5 * Z * wts * C
    e_gl = moment_error(cn2c, zc, w_gl, z_gl, 7)
    assert np.max(np.abs(e_gl)) < 1e-6, e_gl

    # The uniform layering that the book's example uses (11 equal screens at
    # the propagation planes) does NOT match the high moments.
    n_u = 11
    z_u = np.linspace(0.0, Z, n_u)
    w_u = np.full(n_u, C * Z / n_u)
    e_u = moment_error(cn2c, zc, w_u, z_u, 7)
    assert np.max(np.abs(e_u)) > 0.05, e_u

    print('moment match of a uniform 50 km slab, Eq. (9.65), m = 0..7:')
    print(f"  {'m':>3}{'continuous':>14}{'GL-4 err':>12}"
          f"{'11 uniform err':>16}")
    for m in range(8):
        print(f'  {m:>3}{mu[m]:>14.4e}{e_gl[m]:>12.2e}{e_u[m]:>16.2e}')
    print(f'  (trapezium vs exact, max rel err {err_mu:.2e})')
    print('')

    # ---------------- 4. the per-screen cap and the solver ----------------
    # The 50 km example of Sec. 9.5.1, printed p. 176. The steps below repeat
    # Listing 9.5, lines 15 to 22, printed p. 175.
    lam, Dz = 1e-6, 50e3
    kk = 2.0 * np.pi / lam
    r0sw = (0.423 * kk ** 2 * C * 3.0 / 8.0 * Dz) ** (-3.0 / 5.0)
    pgrid = np.linspace(0.0, Dz, 1000)
    rytov = 0.563 * kk ** (7.0 / 6.0) * np.sum(
        C * (1.0 - pgrid / Dz) ** (5.0 / 6.0) * pgrid ** (5.0 / 6.0)
        * (pgrid[1] - pgrid[0]))
    # The log-amplitude variance agrees with the printed 0.436 to 3 digits.
    assert abs(rytov - 0.436) < 0.005, rytov
    # A BOOK DISCREPANCY. Sec. 9.5.1, printed p. 176, prints r0_sw = 17.7 cm.
    # Listing 9.5, lines 15 to 17, printed p. 175, gives 12.66 cm for the same
    # numbers, and Problem 2, printed p. 183, confirms the (3/8)^(-3/5)
    # factor. The listing and the printed sigma_chi^2 agree with each other,
    # so the printed r0_sw is the odd value. This check follows the listing.
    assert abs(r0sw - 0.1266) < 0.001, r0sw

    nscr = 11
    alpha = np.arange(nscr) / (nscr - 1.0)
    r0_scr = screen_strengths(alpha, r0sw, rytov, Dz, lam)
    share = screen_rytov_share(r0_scr, alpha, Dz, lam, 'spherical')
    assert np.max(share) <= RMAX * (1.0 + 1e-9), np.max(share)
    r0_back = composite_r0(r0_scr, alpha, 'spherical')
    assert abs(r0_back / r0sw - 1.0) < 0.02, (r0_back, r0sw)
    assert abs(float(np.sum(share)) / rytov - 1.0) < 0.02, (np.sum(share),
                                                            rytov)

    # The factor 4: the plane-wave share of this module is sigma_chi^2, and
    # the olb rule 2.25 k^(7/6) (INT Cn2 dz) (z - z_i)^(5/6) is sigma_R^2.
    # 2.25 / 0.563 = 3.997.
    i_probe = 5
    # Invert Eq. (9.70), printed p. 165, to get the layer integral Cn2_i dz_i.
    cn2_i = r0_scr[i_probe] ** (-5.0 / 3.0) / (0.423 * kk ** 2)
    chi_pw = float(screen_rytov_share(r0_scr[i_probe], alpha[i_probe], Dz,
                                      lam, 'plane'))
    rytov_pw = 2.25 * kk ** (7.0 / 6.0) * cn2_i * (Dz - alpha[i_probe]
                                                   * Dz) ** (5.0 / 6.0)
    assert abs(rytov_pw / chi_pw - 4.0) < 0.01, (rytov_pw / chi_pw)

    print('screen strengths of the Sec. 9.5.1 example, 11 screens:')
    print(f'  r0_sw target            {r0sw * 1e2:11.3f} cm '
          f'(the book prints 17.7)')
    print(f'  sigma_chi,sw^2 target   {rytov:11.4f}    '
          f'(the book prints 0.436)')
    print(f'  r0_sw from the screens  {r0_back * 1e2:11.3f} cm')
    print(f'  sigma_chi,sw^2 back     {float(np.sum(share)):11.4f}')
    print(f'  largest screen share    {float(np.max(share)):11.4f} '
          f'(cap {RMAX})')
    print(f'  sigma_R^2 / sigma_chi^2 {rytov_pw / chi_pw:11.4f} '
          f'(the factor 4)')
    print(f"  {'i':>3}{'alpha':>8}{'r0_i [cm]':>12}{'share':>10}")
    for i, (aa, rr, ss) in enumerate(zip(alpha, r0_scr, share)):
        print(f'  {i:>3}{aa:>8.3f}{rr * 1e2:>12.2f}{ss:>10.4f}')
    print('')

    # ---------------- 5. the sampling bounds ----------------
    # The 50 km example of Sec. 9.5.2, printed p. 177: c = 2, Delta1 =
    # Delta_n = 1 cm, N = 512, and the book reads n_min = 2 planes.
    D2 = 0.5
    D1 = lam * Dz / (4.0 * D2)
    d1p = blurred_extent(D1, lam, Dz, r0sw, 2.0)
    d2p = blurred_extent(D2, lam, Dz, r0sw, 2.0)
    d1, dn, Ngrid = 10e-3, 10e-3, 512
    assert dn <= constraint1_pitch_max(d1, d1p, d2p, lam, Dz)
    n_need = constraint2_n_min(d1, dn, d1p, d2p, lam, Dz)
    assert Ngrid >= n_need, (Ngrid, n_need)
    # The book: "the required number of grid points is more than 2^8".
    assert 256 < n_need <= 512, n_need
    lo3, hi3 = constraint3_pitch_range(d1, d1p, lam, Dz, Dz)
    assert lo3 <= dn <= hi3, (lo3, dn, hi3)
    zmax = max_partial_step(d1, dn, Ngrid, lam)
    assert min_planes(Dz, zmax) == 2, min_planes(Dz, zmax)

    # The two Sec. 9.4 pitch rules.
    assert abs(fresnel_pitch_max(lam, Dz) - np.sqrt(lam * Dz) / 2) < 1e-15
    # 3 sigma of D(dx) = pi gives dx = 0.332 r0, that is 3.0 pixels per r0.
    px_per_r0 = 1.0 / (phase_pitch_max(1.0) / 1.0)
    assert abs(px_per_r0 - 3.01) < 0.05, px_per_r0

    print('sampling bounds of the Sec. 9.5.2 example, c = 2:')
    print(f'  D1 (point-source lobe)  {D1 * 1e3:11.3f} mm')
    print(f"  D1'                     {d1p:11.4f} m")
    print(f"  D2'                     {d2p:11.4f} m")
    print(f'  constraint 2, N >=      {n_need:11.1f}   (the book picks 512)')
    print(f'  constraint 3 band       [{lo3 * 1e3:.3f}, {hi3 * 1e3:.3f}] mm')
    print(f'  dz_max                  {zmax * 1e-3:11.3f} km')
    print(f'  n_min planes            {min_planes(Dz, zmax):11d}   '
          f'(the book uses 11)')
    print(f'  Fresnel pitch cap       {fresnel_pitch_max(lam, Dz) * 1e3:11.3f}'
          f' mm')
    print(f'  phase pitch cap         {phase_pitch_max(r0sw) * 1e3:11.3f} mm '
          f'({px_per_r0:.2f} px per r0)')
    print('')

    # ---------------- 6. the checklist ----------------
    z_scr = alpha * Dz
    good = properly_sampled_checklist(
        wavelength_m=lam, z_total_m=Dz, n=Ngrid, pixel_m=d1, z_m=z_scr,
        r0_m=r0_scr, r0_total_m=r0sw, d_tx_m=D1, d_rx_m=D2, curvature_m=Dz)
    tested = [row for row in good if row[1] is not None]
    advisory = [row for row in good if row[1] is None]
    assert len(advisory) == 3, advisory
    # This plan fails the WEAK-fluctuation step on purpose: the book's own
    # example has sigma_chi,sw^2 = 0.436, past the 0.25 threshold.
    by_name = {row[0]: row for row in good}
    assert by_name['9.5.1 weak fluctuation'][1] is False
    assert by_name['9.5.1 per-screen Rytov cap'][1] is True
    assert by_name['9.5.2 constraint 1'][1] is True
    assert by_name['9.5.2 constraint 2'][1] is True
    assert by_name['9.5.2 constraint 3'][1] is True
    assert by_name['9.4 phase pitch'][1] is True
    assert by_name['9.4 scintillation pitch'][1] is True
    assert by_name['9.5.2 partial-propagation planes'][1] is True
    assert sum(1 for row in good if row[1] is False) == 1

    # A hand-built FAILING plan: a 10 times coarser pixel on a 4 times smaller
    # grid. It must break both geometry constraints and both pitch rules.
    bad = properly_sampled_checklist(
        wavelength_m=lam, z_total_m=Dz, n=128, pixel_m=10 * d1, z_m=z_scr,
        r0_m=r0_scr, r0_total_m=r0sw, d_tx_m=D1, d_rx_m=D2, curvature_m=Dz)
    bad_by_name = {row[0]: row for row in bad}
    assert bad_by_name['9.5.2 constraint 1'][1] is False
    assert bad_by_name['9.5.2 constraint 3'][1] is False
    assert bad_by_name['9.4 phase pitch'][1] is False
    assert bad_by_name['9.4 scintillation pitch'][1] is False
    assert sum(1 for row in bad if row[1] is False) == 5

    def show(title, rows):
        """Print one checklist."""
        print(title)
        for name, ok, bound, actual, _ in rows:
            if ok is None:
                print(f'  {name:<36} advisory')
                continue
            mark = 'pass' if ok else 'FAIL'
            b = (f'[{bound[0]:.4g}, {bound[1]:.4g}]'
                 if isinstance(bound, tuple) else f'{bound:.4g}')
            print(f'  {name:<36} {mark}   bound {b:>22}   '
                  f'actual {actual:.4g}')

    show('checklist, the book example (N=512, 10 mm pixel):', good)
    print('')
    show('checklist, a coarse plan (N=128, 100 mm pixel):', bad)
    print('')
    print(f'(elapsed {time.time() - t_start:.1f} s)')
    print('self-check passed')
