'''
The Fresnel propagation kernels of Schmidt, and the partial-propagation loop.

This module holds the three vacuum propagators of Chapter 6, the absorbing
boundary of Section 8.1, and the partial-propagation loop of Sections 8.2 and
8.3. Each kernel moves a complex field from a source plane to a parallel
observation plane.

Source of every equation:
    J. D. Schmidt, "Numerical Simulation of Optical Wave Propagation with
    Examples in MATLAB", SPIE Press Monograph PM199 (2010).
    DOI: 10.1117/3.866274
Chapter 6, printed pp. 87 to 113, and Chapter 8, printed pp. 133 to 147.
Each function names its chapter, its equation number, and its printed page.

THE FOUR KERNELS.

    one_step_fresnel     Ch. 6.3.1. One transform. The output pitch is FIXED.
    two_step_fresnel     Ch. 6.3.2. Two transforms. The output pitch is FREE.
    angular_spectrum     Ch. 6.4.   Two transforms. The output pitch is FREE.
    partial_propagations Ch. 8.3.   A loop of angular-spectrum steps.

WHICH KERNEL. The book gives the rule at Ch. 7, Sec. 7.3, printed p. 120: the
Fresnel-integral method of Sec. 6.3 is valid for a LONG propagation, and the
angular-spectrum method of Sec. 6.4 is valid for a SHORT propagation. The two
methods sample two different quadratic-phase factors, so they alias in opposite
limits.

THE PISTON PHASE. Every kernel here DROPS the factor exp(i k z) of Ch. 6,
Eq. (6.5), printed p. 88, and of Ch. 6, Eq. (6.32), printed p. 95. The book's
own listings drop it too (Listings 6.1, 6.3 and 6.5, printed pp. 91, 96 and
102). The factor is a constant across the plane, so it changes no irradiance
and no relative phase. The four kernels agree with each other because they all
drop it. Add exp(i k z) if you need the absolute phase.

THE SAMPLING RULES ARE NOT HERE. Each docstring names the Chapter 7 constraint
that governs the kernel, by equation number. This module does NOT test a
constraint and does NOT refuse a bad grid. The module `sampling.py` owns those
tests. A caller that ignores the constraints gets an aliased result and no
warning.

UNITS. `wavelength` is in metres. Every pitch and every distance is in metres.
A field carries the units of the square root of an irradiance.

THE GRID. Every function needs a SQUARE grid of N by N samples. The origin sits
at index N//2 on each axis, which matches the convention of `fourier.ft2`. The
book uses the index range -N/2 to N/2-1 (Listing 6.1, line 9, printed p. 91).
The two conventions are the same for an even N.

This module holds physics only. It imports numpy and `.fourier` only. It
returns no decibels.
'''

import numpy as np

from .fourier import freq_pitch, ft2, ift2


def _grid(n, dx):
    '''
    Return the squared radius of an n by n grid of the pitch dx.

    The axis runs from -N/2 dx to (N/2 - 1) dx, as Listing 6.1, line 9,
    printed p. 91, does.
    '''
    axis = (np.arange(n) - n // 2) * float(dx)
    return axis[:, None] ** 2 + axis[None, :] ** 2


def _check_square(Uin, name):
    '''Return the side count of a square two-dimensional array.'''
    Uin = np.asarray(Uin)
    if Uin.ndim != 2 or Uin.shape[0] != Uin.shape[1]:
        raise ValueError(f'{name} needs a square two-dimensional grid; see '
                         'Schmidt (2010), DOI 10.1117/3.866274, Ch. 2, '
                         'Sec. 2.6, printed p. 36')
    return Uin.shape[0]


def one_step_fresnel(Uin, wavelength, dx1, z):
    '''
    Propagate a field a distance z with ONE Fresnel transform.

    Parameters:
        Uin : numpy.ndarray
            The source-plane field. Square and complex.
        wavelength : float
            The optical wavelength [m].
        dx1 : float
            The source-plane grid pitch [m].
        z : float
            The propagation distance [m]. It must not be zero.

    Returns:
        (numpy.ndarray, float)
            The observation-plane field, and the observation-plane pitch dx2
            [m]. The caller does NOT choose dx2.

    formula:
        U(r2) = exp(i k r2^2 / (2 z)) / (i lambda z)
                * FT{ U(r1) exp(i k r1^2 / (2 z)) }
        dx2   = lambda z / (N dx1)
    Source: Schmidt (2010), DOI 10.1117/3.866274, Ch. 6, Eq. (6.5), printed
    p. 88, gives the transform form of the Fresnel integral. Ch. 6, Eq. (6.15),
    printed p. 90, gives the same thing as the operator chain Q V F Q. Ch. 6,
    Eq. (6.16), printed p. 90, gives the FIXED observation pitch. Ch. 7,
    Eq. (7.21), printed p. 120, repeats that pitch.

    VALIDITY.

    - THE FRESNEL APPROXIMATION. The kernel solves Ch. 1, Eq. (1.57), printed
      p. 10. That integral holds in the PARAXIAL approximation of Ch. 1,
      Eqs. (1.49) and (1.50), printed p. 8: the direction cosines to the axis
      are close to 1. The book states at Ch. 6, text, printed p. 87, that the
      planes must be PARALLEL and that the paraxial approximation is then good.
      The book gives NO numerical threshold for the approximation.
    - WHAT THE KERNEL FIXES. The observation pitch is dx2 = lambda z /(N dx1).
      One step gives NO control of dx2 (Ch. 6, text below Eq. (6.16), printed
      p. 92). To change dx2, change N, dx1, or z. Use `two_step_fresnel` or
      `angular_spectrum` to set dx2 directly.
    - WHY IT ALIASES. The quadratic phase exp(i k r^2 /(2 z)) is NOT band
      limited (Ch. 6, text, printed p. 87). Its local frequency grows with the
      radius (Ch. 7, Eq. (7.37), printed p. 122). The transform aliases the
      chirp when the phase step between two neighbour samples passes pi at the
      edge of the source field. A SHORT z aliases first, because the chirp then
      turns fastest. This is why the method suits a LONG propagation (Ch. 7,
      Sec. 7.3, printed p. 120).
    - THE GOVERNING CONSTRAINTS. Ch. 7, Eq. (7.14), printed p. 119, is
      CONSTRAINT 1. Ch. 7, Eq. (7.20), printed p. 120, is CONSTRAINT 2. With
      the fixed pitch of Eq. (7.21) the two collapse to the SAME minimum grid
      count, Ch. 7, Eqs. (7.25) and (7.31), printed p. 121. Ch. 7, Eqs. (7.41)
      and (7.42), printed p. 123, add the MINIMUM distance that keeps the
      source chirp sampled. Constraints 3 and 4 do NOT apply, because this
      method has no transfer function. `sampling.py` owns those tests.
    '''
    n = _check_square(Uin, 'one_step_fresnel')
    z = float(z)
    if z == 0.0:
        raise ValueError('one_step_fresnel: z must not be zero; see Schmidt '
                         '(2010), DOI 10.1117/3.866274, Ch. 6, Eq. (6.16), '
                         'printed p. 90')
    k = 2.0 * np.pi / float(wavelength)

    # Ch. 6, Eq. (6.16), printed p. 90: the fixed observation pitch.
    dx2 = float(wavelength) * z / (n * float(dx1))

    r1sq = _grid(n, dx1)
    r2sq = _grid(n, dx2)

    # Ch. 6, Eq. (6.5), printed p. 88, without the piston factor exp(i k z).
    Uout = (np.exp(1j * k / (2.0 * z) * r2sq) / (1j * float(wavelength) * z)
            * ft2(np.asarray(Uin) * np.exp(1j * k / (2.0 * z) * r1sq), dx1))
    return Uout, dx2


def two_step_fresnel(Uin, wavelength, dx1, dx2, z):
    '''
    Propagate a field a distance z with TWO Fresnel transforms.

    The method sends the field to an intermediate plane and then to the
    observation plane. The position of that plane sets the output pitch, so the
    caller chooses dx2.

    Parameters:
        Uin : numpy.ndarray
            The source-plane field. Square and complex.
        wavelength : float
            The optical wavelength [m].
        dx1 : float
            The source-plane grid pitch [m].
        dx2 : float
            The wanted observation-plane grid pitch [m]. It must differ from
            dx1.
        z : float
            The propagation distance [m]. It must not be zero.

    Returns:
        numpy.ndarray
            The observation-plane field, on a grid of the pitch dx2.

    formula:
        m   = dx2 / dx1
        z1  = z / (1 - m)                 the first step
        dxa = lambda |z1| / (N dx1)       the intermediate pitch
        z2  = z - z1                      the second step
        Each step is the one-step kernel above.
    Source: Schmidt (2010), DOI 10.1117/3.866274, Ch. 6, Sec. 6.3.2, printed
    pp. 92 to 95. The operator chain is Eq. (6.18), printed p. 93. The pitch
    chain is Eqs. (6.19) to (6.21), printed p. 94. The scaling parameter is
    Eq. (6.24), printed p. 94. The intermediate-plane position is Eq. (6.25),
    printed p. 94. Eq. (6.30), printed p. 95, proves that the two branches give
    the same |m|. The method is due to Coy, and to Rydberg and Bengtsson
    (Ch. 6, text, printed p. 92).

    THE m CONSTRAINTS.

    - m = 1 is FORBIDDEN. Eq. (6.25), printed p. 94, then puts the intermediate
      plane at infinity, or exactly halfway (Table 6.2, printed p. 95). This
      function refuses m = 1. Use `angular_spectrum`, which holds m = 1.
    - m must be POSITIVE. m is a ratio of two grid pitches, and a pitch is
      positive.
    - m < 1 puts the intermediate plane BEYOND the observation plane, so z1 > z
      and z2 < 0 (Fig. 6.2, printed p. 93).
    - m > 1 puts the intermediate plane BEFORE the source plane, so z1 < 0 and
      z2 > z (Fig. 6.3, printed p. 94).
      A negative step is correct. The Fresnel integral holds for a negative
      distance, because it only reverses the sign of the chirp.
    - This function takes the MINUS branch of Eq. (6.25), z1 = z / (1 - m), as
      Listing 6.3, line 13, printed p. 96, does. Eq. (6.30) shows that the plus
      branch gives the same |m|.

    VALIDITY.

    - THE FRESNEL APPROXIMATION. The same as `one_step_fresnel`. The kernel
      applies Ch. 1, Eq. (1.57), printed p. 10, two times, so the paraxial
      condition of Ch. 1, Eqs. (1.49) and (1.50), printed p. 8, must hold on
      BOTH steps.
    - WHAT THE KERNEL FIXES. It frees the output pitch. The cost is a second
      transform, and an intermediate grid whose pitch dxa the caller does not
      control.
    - WHY IT ALIASES. Each step aliases like one step. The chirp of a step is
      exp(i k r^2 /(2 z_i)), and z1 or z2 can be MUCH shorter than z when m is
      close to 1. Look at Eq. (6.25): z1 goes to infinity as m goes to 1, and
      the second step then carries the whole aliasing risk. The intermediate
      grid also has a side of N dxa, and the field must stay inside it.
    - THE GOVERNING CONSTRAINTS. The book analyses ONE step only (Ch. 7,
      Sec. 7.3.1.1, printed p. 120). It gives NO separate constraint set for
      two steps. Apply the one-step rules, Ch. 7, Eqs. (7.14), (7.20), (7.25),
      (7.41) and (7.42), printed pp. 119 to 123, to EACH of the two steps, with
      the pair (dx1, dxa, z1) and then the pair (dxa, dx2, z2). `sampling.py`
      owns those tests.
    '''
    n = _check_square(Uin, 'two_step_fresnel')
    z = float(z)
    if z == 0.0:
        raise ValueError('two_step_fresnel: z must not be zero')
    m = float(dx2) / float(dx1)
    if m <= 0.0:
        raise ValueError('two_step_fresnel: m = dx2/dx1 must be positive')
    if m == 1.0:
        raise ValueError('two_step_fresnel: m = dx2/dx1 must not be 1; the '
                         'intermediate plane then goes to infinity. See '
                         'Schmidt (2010), DOI 10.1117/3.866274, Ch. 6, '
                         'Eq. (6.25) and Table 6.2, printed pp. 94 and 95. '
                         'Use angular_spectrum instead')
    k = 2.0 * np.pi / float(wavelength)

    # Ch. 6, Eq. (6.25), printed p. 94: the minus branch.
    z1 = z / (1.0 - m)
    # Ch. 6, Eq. (6.19), printed p. 94: the intermediate pitch.
    dxa = float(wavelength) * abs(z1) / (n * float(dx1))
    # Ch. 6, Eq. (6.23), printed p. 94.
    z2 = z - z1

    r1sq = _grid(n, dx1)
    rasq = _grid(n, dxa)
    r2sq = _grid(n, dx2)

    # Step 1. Ch. 6, Eq. (6.18), printed p. 93, right half.
    Uitm = (np.exp(1j * k / (2.0 * z1) * rasq) / (1j * float(wavelength) * z1)
            * ft2(np.asarray(Uin) * np.exp(1j * k / (2.0 * z1) * r1sq), dx1))
    # Step 2. Ch. 6, Eq. (6.18), printed p. 93, left half.
    Uout = (np.exp(1j * k / (2.0 * z2) * r2sq) / (1j * float(wavelength) * z2)
            * ft2(Uitm * np.exp(1j * k / (2.0 * z2) * rasq), dxa))
    return Uout


def angular_spectrum(Uin, wavelength, dx, z, dx2=None):
    '''
    Propagate a field a distance z with the angular-spectrum method.

    The method transforms the field, multiplies the spectrum by the free-space
    transfer function, and transforms back. With `dx2` it also scales the
    output grid.

    Parameters:
        Uin : numpy.ndarray
            The source-plane field. Square and complex.
        wavelength : float
            The optical wavelength [m].
        dx : float
            The source-plane grid pitch [m].
        z : float
            The propagation distance [m]. It must not be zero.
        dx2 : float or None
            The wanted observation-plane grid pitch [m]. None keeps the source
            pitch, that is m = 1. This is the baseline form of Ch. 6,
            Eqs. (6.31) and (6.32).

    Returns:
        numpy.ndarray
            The observation-plane field, on a grid of the pitch dx2, or of the
            pitch dx when dx2 is None.

    formula, the baseline form with m = 1:
        U(r2) = IFT{ H(f) FT{ U(r1) } },
        H(f)  = exp(-i pi lambda z f^2)
    formula, the scaled form with m = dx2/dx1:
        Q1 = exp(i (k/2) (1 - m) r1^2 / z)
        Q2 = exp(-i 2 pi^2 z f^2 / (m k))
        Q3 = exp(i (k/2) (m - 1) r2^2 / (m z))
        U(r2) = Q3 IFT{ Q2 FT{ Q1 U(r1) / m } }
    Source: Schmidt (2010), DOI 10.1117/3.866274, Ch. 6, Sec. 6.4, printed
    pp. 95 to 102. Eq. (6.31), printed p. 95, is the convolution form.
    Eq. (6.32), printed p. 95, is the free-space transfer function; this
    function drops its exp(i k z) factor. Eqs. (6.34) to (6.49), printed
    pp. 98 and 99, introduce m. Eq. (6.50), printed p. 99, is the SCALED
    operator chain, and Eq. (6.65), printed p. 100, repeats it. The pitch chain
    is Eqs. (6.51) to (6.54), printed pp. 99 and 100. The scaling is due to
    Tyler and Fried, and to Roberts (Ch. 6, text, printed p. 98).

    Q2 above is Ch. 6, Eq. (6.12), printed p. 89, with the argument -z/m. With
    k = 2 pi / lambda it reduces to exp(-i pi lambda z f^2 / m), so m = 1 gives
    exactly H(f) of Eq. (6.32) without its piston factor. So the scaled form
    contains the baseline form.

    VALIDITY.

    - THE FRESNEL APPROXIMATION. The kernel solves the CONVOLUTION form of the
      Fresnel integral, Ch. 6, Eq. (6.6), printed p. 88. That integral is the
      same Ch. 1, Eq. (1.57), printed p. 10, so the paraxial condition of
      Ch. 1, Eqs. (1.49) and (1.50), printed p. 8, holds here too. The transfer
      function of Eq. (6.32) also assumes PARALLEL source and observation
      planes (Ch. 6, text, printed p. 87).
    - WHAT THE KERNEL FIXES. With dx2 = None the output pitch EQUALS the input
      pitch (Ch. 6, text, printed p. 96: "we would be stuck with dx1 = dx2").
      With dx2 the pitch chain of Eq. (6.54), printed p. 100, gives dx2 = m dx.
    - WHY IT ALIASES. The method must sample TWO quadratic phases (Ch. 7,
      Eqs. (7.46) and (7.47), printed p. 125). The first is the source chirp
      Q1, and the second is the transfer function Q2. Q2 turns faster as z
      grows, so a LONG z aliases the transfer function at the edge of the
      frequency grid. This is why the method suits a SHORT propagation (Ch. 7,
      Sec. 7.3, printed p. 120). Chapter 8 relaxes the limit by splitting the
      path.
    - THE GOVERNING CONSTRAINTS. Ch. 7, Eq. (7.14), printed p. 119, is
      CONSTRAINT 1, and Ch. 7, Eq. (7.20), printed p. 120, is CONSTRAINT 2.
      Both are geometric. Ch. 7, Eq. (7.53), printed p. 126, is CONSTRAINT 3,
      which samples Q1 and needs the source wavefront radius R. Ch. 7,
      Eq. (7.59), printed p. 127, is CONSTRAINT 4, N >= lambda z /(dx1 dx2),
      which samples Q2. Constraint 4 is the one that fails on a long path.
      `sampling.py` owns those tests.
    '''
    n = _check_square(Uin, 'angular_spectrum')
    z = float(z)
    if z == 0.0:
        raise ValueError('angular_spectrum: z must not be zero')
    dx = float(dx)
    m = 1.0 if dx2 is None else float(dx2) / dx
    if m <= 0.0:
        raise ValueError('angular_spectrum: m = dx2/dx must be positive')
    k = 2.0 * np.pi / float(wavelength)
    df = freq_pitch(n, dx)

    r1sq = _grid(n, dx)
    r2sq = _grid(n, m * dx)
    fsq = _grid(n, df)

    # Ch. 6, Eq. (6.65), printed p. 100. Listing 6.5, printed p. 102, writes
    # the same three factors.
    Q1 = np.exp(1j * k / 2.0 * (1.0 - m) / z * r1sq)
    Q2 = np.exp(-1j * np.pi ** 2 * 2.0 * z / m / k * fsq)
    Q3 = np.exp(1j * k / 2.0 * (m - 1.0) / (m * z) * r2sq)
    return Q3 * ift2(Q2 * ft2(Q1 * np.asarray(Uin) / m, dx), df)


def super_gaussian_absorber(n, sigma_frac=0.47, power=16):
    '''
    Return the super-Gaussian absorbing boundary of the book.

    The boundary removes the light that reaches the edge of the grid. The
    discrete propagators are periodic, so that light comes back at the opposite
    edge. The boundary is close to 1 in the middle of the grid, and close to 0
    at the edge.

    Parameters:
        n : int
            Grid points per side.
        sigma_frac : float
            The half-width sigma, as a fraction of the grid side N. The book
            uses 0.47.
        power : float
            The exponent of the super-Gaussian. The book uses 16. Eq. (8.1)
            needs a value above 2.

    Returns:
        numpy.ndarray
            An n by n array of real numbers in the interval (0, 1]. The value
            at the centre is exactly 1.

    formula:
        g_sg(x,y) = exp( -(r / sigma)^power ),   power > 2
    with r the radius in PIXELS and sigma = sigma_frac N pixels.
    Source: Schmidt (2010), DOI 10.1117/3.866274, Ch. 8, Eq. (8.1), printed
    p. 134. Listing 8.1, lines 10 to 12, printed p. 142, gives the reference
    values: `w = 0.47*N` and `sg = exp(-nsq.^8/w^16)`, with `nsq` the squared
    pixel radius. That is a super-Gaussian of the power 16 in r, and of the
    half-width 0.47 N pixels. Figure 8.1, printed p. 134, plots sigma = 0.45 L
    and n = 16. The listing value 0.47 N is the one that the book RUNS, so this
    function takes it as the default.

    THE BOOK'S OTHER BOUNDARY. Ch. 8, Eq. (8.2), printed p. 134, gives the
    Tukey (cosine-taper) window as an alternative. This module does not build
    it. Nothing needs it.

    A DIFFERENT SHAPE LIVES ELSEWHERE. The function
    `olb.waveoptics.turbulence.splitstep.super_gaussian_boundary` uses power 8
    and a taper band of 0.125 of the HALF-side. That is a different
    parameterisation and a different shape. Do not compare the two numbers.

    VALIDITY.

    - The boundary is a NUMERICAL device, not physics. It removes energy, so
      the field after it does not conserve power. Keep the region of interest
      well inside the flat part. Ch. 8, text, printed p. 134: "we must be
      careful not to alter light in the central region of the grid".
    - Eq. (8.1) needs power > 2. A power of 2 gives a plain Gaussian, which
      attenuates the middle of the grid.
    - The book applies the boundary at EVERY partial-propagation plane
      (Eq. (8.18), printed p. 139, and Listing 8.1, line 38, printed p. 142).
    - The book gives no rule that ties sigma to the region of interest. The
      values 0.47 and 16 are the book's own numbers, not a derivation.
    '''
    n = int(n)
    sigma_frac = float(sigma_frac)
    power = float(power)
    if sigma_frac <= 0.0:
        raise ValueError('super_gaussian_absorber: sigma_frac must be '
                         'positive')
    if power <= 2.0:
        raise ValueError('super_gaussian_absorber: power must be above 2; see '
                         'Schmidt (2010), DOI 10.1117/3.866274, Ch. 8, '
                         'Eq. (8.1), printed p. 134')
    # The radius in PIXELS, on the index range -N/2 to N/2-1 of Listing 8.1.
    idx = np.arange(n) - n // 2
    r_pix = np.hypot(idx[:, None], idx[None, :])
    sigma = sigma_frac * n
    # Ch. 8, Eq. (8.1), printed p. 134.
    return np.exp(-np.minimum((r_pix / sigma) ** power, 700.0))


def partial_propagations(Uin, wavelength, dx1, dxn, z_planes, absorber=None):
    '''
    Propagate a field through n planes with n-1 angular-spectrum steps.

    The loop splits one long path into short steps. Each step has its own grid
    pitch, so the grid grows or shrinks with the beam. The absorbing boundary
    goes on each plane after the first.

    Parameters:
        Uin : numpy.ndarray
            The source-plane field. Square and complex.
        wavelength : float
            The optical wavelength [m].
        dx1 : float
            The pitch of the FIRST plane [m].
        dxn : float
            The pitch of the LAST plane [m].
        z_planes : sequence of float
            The positions of the n planes along the axis [m]. Give the source
            plane first. The book's `z` argument starts at the SECOND plane and
            the code prepends a zero (Listing 8.1, line 14, printed p. 142).
            This function takes the full list, so the first value is normally
            0. The values must increase.
        absorber : numpy.ndarray or None
            The boundary mask, from `super_gaussian_absorber`. None applies no
            boundary, which gives the pure vacuum result.

    Returns:
        numpy.ndarray
            The field in the last plane, on a grid of the pitch dxn.

    formula:
        alpha_i = z_i / Delta_z
        Delta_i = (1 - alpha_i) dx1 + alpha_i dxn        the LINEAR pitch rule
        m_i     = Delta_i+1 / Delta_i
        U <- Q[(1 - m_1)/Delta_z_1, r1] U(r1)            once, at the start
        loop i: U <- A[r_i+1] IFT{ Q2[-Delta_z_i/m_i, f_i] FT{ U / m_i } }
        U <- Q[(m_n-1 - 1)/(m_n-1 Delta_z_n-1), rn] U    once, at the end
    Source: Schmidt (2010), DOI 10.1117/3.866274, Ch. 8, Sec. 8.3, printed
    pp. 138 and 139. Eq. (8.18), printed p. 139, is the general chain.
    Eq. (8.8), printed p. 136, is the linear pitch rule, and Eqs. (8.3) to
    (8.7), printed p. 136, prove it from similar triangles. Eq. (8.9), printed
    p. 137, is the two-step case with the boundary operator A. Eqs. (8.14) and
    (8.15), printed p. 138, show that the two middle-plane quadratic phases are
    inverse, so Eq. (8.16), printed p. 138, keeps ONE quadratic phase at each
    END only. Listing 8.1, printed p. 142, writes the loop.

    WHY THE MIDDLE PHASES DROP. Eq. (8.15), printed p. 138, proves that
    (m_i - 1)/(m_i Delta_z_i) is the same number for every step, when the
    pitches follow the linear rule of Eq. (8.8). So the exit phase of one step
    cancels the entry phase of the next. This is why the loop body carries only
    Q2, and why the pitch rule is not free.

    VALIDITY.

    - THE FRESNEL APPROXIMATION. Every step is an angular-spectrum step, so the
      paraxial condition of Ch. 1, Eqs. (1.49) and (1.50), printed p. 8, holds
      on every step.
    - THE PITCH RULE IS NOT FREE. Eq. (8.8), printed p. 136, makes the pitch
      LINEAR in the distance. A caller cannot give an arbitrary pitch per
      plane. The cancellation above needs the linear rule.
    - WHAT THE LOOP BUYS. It relaxes CONSTRAINT 4 only. Ch. 8, Eqs. (8.19) to
      (8.22), printed pp. 139 and 143, prove that CONSTRAINT 3 does NOT change
      with the plane count. Constraints 1 and 2 are geometric, so they do not
      change either (Ch. 8, text, printed p. 139).
    - THE GOVERNING CONSTRAINTS. Ch. 8, Eq. (8.23), printed p. 143, gives
      constraint 4 per step: N >= lambda Delta_z_i /(Delta_i Delta_i+1).
      Ch. 8, Eq. (8.24), printed p. 144, turns it into the step cap
      Delta_z_max = min(dx1, dxn)^2 N / lambda, and the plane count
      n >= ceil(Delta_z / Delta_z_max) + 1. More planes are always allowed
      (Ch. 8, text, printed p. 144). `sampling.py` owns those tests. This
      function does NOT check the plane count.
    - THE BOUNDARY REMOVES ENERGY. See `super_gaussian_absorber`. The result
      does not conserve power when a boundary is given.
    - This function is a REFERENCE. It is not a production runner. It holds no
      turbulence: Ch. 9, Eq. (9.3), printed p. 150, adds the screen operator T
      in place of A.
    '''
    n_pix = _check_square(Uin, 'partial_propagations')
    z = np.asarray(z_planes, dtype=float)
    if z.ndim != 1 or z.size < 2:
        raise ValueError('partial_propagations: give at least two planes; see '
                         'Schmidt (2010), DOI 10.1117/3.866274, Ch. 8, '
                         'Eq. (8.18), printed p. 139')
    if np.any(np.diff(z) <= 0.0):
        raise ValueError('partial_propagations: z_planes must increase')
    n_plane = z.size
    k = 2.0 * np.pi / float(wavelength)

    # Ch. 8, Eq. (8.8), printed p. 136: the linear pitch rule.
    alpha = (z - z[0]) / (z[-1] - z[0])
    delta = (1.0 - alpha) * float(dx1) + alpha * float(dxn)
    dz = np.diff(z)
    m = delta[1:] / delta[:-1]

    # Ch. 8, Eq. (8.18), printed p. 139: the entry quadratic phase.
    U = np.asarray(Uin) * np.exp(
        1j * k / 2.0 * (1.0 - m[0]) / dz[0] * _grid(n_pix, delta[0]))

    for i in range(n_plane - 1):
        df = freq_pitch(n_pix, delta[i])
        fsq = _grid(n_pix, df)
        # Ch. 6, Eq. (6.12), printed p. 89, with the argument -Delta_z_i/m_i.
        Q2 = np.exp(-1j * np.pi ** 2 * 2.0 * dz[i] / m[i] / k * fsq)
        U = ift2(Q2 * ft2(U / m[i], delta[i]), df)
        if absorber is not None:
            # The operator A of Ch. 8, Eq. (8.18), printed p. 139. The book
            # applies it on EVERY plane after the first, the last one included
            # (Listing 8.1, line 38, printed p. 142).
            U = U * absorber

    # Ch. 8, Eq. (8.18), printed p. 139: the exit quadratic phase.
    return U * np.exp(1j * k / 2.0 * (m[-1] - 1.0) / (m[-1] * dz[-1])
                      * _grid(n_pix, delta[-1]))


if __name__ == '__main__':
    # ---------------- physics self-checks ----------------
    #
    # THE ANALYTIC TARGET. The free-space Gaussian-beam wave is the only closed
    # form that all four kernels can hit. With the amplitude U(r,0) =
    # exp(-r^2/W0^2), the field after a distance z is
    #
    #     U(r,z) = (W0/W) exp(-r^2/W^2) exp(i k r^2/(2 R)) exp(-i arctan(z/zR))
    #     W  = W0 sqrt(1 + (z/zR)^2),  R = z (1 + (zR/z)^2),  zR = pi W0^2/lam
    #
    # Source: Schmidt (2010), DOI 10.1117/3.866274, Ch. 1, Eqs. (1.53) to
    # (1.56), printed p. 9. Eq. (1.53) is U_G = (A/q) exp(i k r^2 /(2 q)), and
    # Eq. (1.54) is 1/q = 1/R + i lambda /(pi W^2). Put Eq. (1.54) into
    # Eq. (1.53) to get the amplitude and the two phases above. Eqs. (1.55) and
    # (1.56) give W(z) and R(z). Andrews and Phillips, 2nd ed. (2005),
    # DOI 10.1117/3.626196, Ch. 4, Eqs. (37) and (38), printed p. 93, print the
    # same solution.
    # The piston factor exp(i k z) is REMOVED, because the kernels drop it.
    #
    # THE GEOMETRY. It is common to the three Chapter 6 kernels. The one-step
    # pitch of Eq. (6.16) is fixed, so the geometry is chosen to make that
    # pitch the target for the other two. Then m = dx2/dx1 is the same for all
    # three, and the three results sit on ONE grid.
    N = 1024
    lam = 1.0e-6                      # [m]
    dx1 = 8.0e-5                      # [m]
    z_prop = 10.0                     # [m]
    W0 = 3.0e-3                       # [m]

    dx2_fixed = lam * z_prop / (N * dx1)          # Ch. 6, Eq. (6.16)
    m_common = dx2_fixed / dx1

    r1sq_chk = _grid(N, dx1)
    U0 = np.exp(-r1sq_chk / W0 ** 2)

    zR = np.pi * W0 ** 2 / lam
    W_z = W0 * np.sqrt(1.0 + (z_prop / zR) ** 2)
    R_z = z_prop * (1.0 + (zR / z_prop) ** 2)
    k_chk = 2.0 * np.pi / lam
    r2sq_chk = _grid(N, dx2_fixed)
    U_exact = (W0 / W_z * np.exp(-r2sq_chk / W_z ** 2)
               * np.exp(1j * k_chk * r2sq_chk / (2.0 * R_z))
               * np.exp(-1j * np.arctan(z_prop / zR)))

    # Compare inside three beam radii. The field outside that is below 1e-4 of
    # the peak, so a relative reading there is noise.
    core = r2sq_chk <= (3.0 * W_z) ** 2
    assert core.sum() > 1000, 'the comparison region must hold real samples'
    peak = float(np.max(np.abs(U_exact)))

    def _err(U):
        '''Return the largest field error inside the core, over the peak.'''
        return float(np.max(np.abs(U[core] - U_exact[core]))) / peak

    print(f'GEOMETRY N = {N}, lambda = {lam * 1e6:.1f} um, dx1 = '
          f'{dx1 * 1e6:.1f} um, z = {z_prop:.1f} m, W0 = {W0 * 1e3:.1f} mm')
    print(f'         dx2 = lambda z /(N dx1) = {dx2_fixed * 1e6:.3f} um, '
          f'm = {m_common:.4f}, W(z) = {W_z * 1e3:.3f} mm')

    # 1. one_step_fresnel against the closed form.
    U_one, dx2_one = one_step_fresnel(U0, lam, dx1, z_prop)
    assert abs(dx2_one - dx2_fixed) < 1e-18
    err_one = _err(U_one)
    assert err_one < 1e-12, err_one
    print(f'REDUCTION one_step_fresnel vs the Gaussian closed form : '
          f'max rel err = {err_one:.3e}  (target 1e-12)')

    # 2. two_step_fresnel against the closed form, on the same output pitch.
    U_two = two_step_fresnel(U0, lam, dx1, dx2_fixed, z_prop)
    err_two = _err(U_two)
    assert err_two < 1e-12, err_two
    print(f'REDUCTION two_step_fresnel vs the Gaussian closed form : '
          f'max rel err = {err_two:.3e}  (target 1e-12)')

    # 3. angular_spectrum, the SCALED form, against the closed form.
    U_as = angular_spectrum(U0, lam, dx1, z_prop, dx2=dx2_fixed)
    err_as = _err(U_as)
    assert err_as < 1e-12, err_as
    print(f'REDUCTION angular_spectrum (m = {m_common:.4f}) vs the Gaussian '
          f'closed form : max rel err = {err_as:.3e}  (target 1e-12)')

    # 4. angular_spectrum, the BASELINE form with dx2 = dx1, against the closed
    # form on its own grid. Constraint 4 of Ch. 7, Eq. (7.59), printed p. 127,
    # reads N >= lambda z /(dx1 dx1) = 1563 here, so the transfer function
    # aliases at N = 1024. Shorten the step to satisfy it: z = 4 m gives
    # N >= 625.
    z_short = 4.0
    W_s = W0 * np.sqrt(1.0 + (z_short / zR) ** 2)
    R_s = z_short * (1.0 + (zR / z_short) ** 2)
    U_exact_s = (W0 / W_s * np.exp(-r1sq_chk / W_s ** 2)
                 * np.exp(1j * k_chk * r1sq_chk / (2.0 * R_s))
                 * np.exp(-1j * np.arctan(z_short / zR)))
    U_base = angular_spectrum(U0, lam, dx1, z_short)
    core_s = r1sq_chk <= (3.0 * W_s) ** 2
    err_base = (float(np.max(np.abs(U_base[core_s] - U_exact_s[core_s])))
                / float(np.max(np.abs(U_exact_s))))
    assert err_base < 1e-12, err_base
    n_need = lam * z_short / dx1 ** 2
    print(f'REDUCTION angular_spectrum (m = 1, z = {z_short:.1f} m) vs the '
          f'Gaussian closed form : max rel err = {err_base:.3e}  '
          f'(target 1e-12, constraint 4 needs N >= {n_need:.0f})')

    # 5. The three Chapter 6 kernels against EACH OTHER, on the common
    # geometry. They solve the same integral, so they must agree better than
    # any one of them agrees with the closed form.
    d_12 = float(np.max(np.abs(U_one[core] - U_two[core]))) / peak
    d_13 = float(np.max(np.abs(U_one[core] - U_as[core]))) / peak
    d_23 = float(np.max(np.abs(U_two[core] - U_as[core]))) / peak
    assert max(d_12, d_13, d_23) < 1e-12, (d_12, d_13, d_23)
    print(f'REDUCTION kernel against kernel : one-two = {d_12:.3e}, '
          f'one-AS = {d_13:.3e}, two-AS = {d_23:.3e}  (target 1e-12)')

    # 6. The absorber. Check the interior, the cited edge value, and the decay.
    sg = super_gaussian_absorber(N)
    assert sg[N // 2, N // 2] == 1.0
    idx_chk = np.arange(N) - N // 2
    r_pix_chk = np.hypot(idx_chk[:, None], idx_chk[None, :])
    inner = r_pix_chk <= 0.2 * N
    assert float(np.min(sg[inner])) > 1.0 - 1e-5, float(np.min(sg[inner]))
    # The middle of an edge sits at the pixel radius N/2. Listing 8.1, printed
    # p. 142, gives sigma = 0.47 N, so g = exp(-(0.5/0.47)^16).
    edge_book = float(np.exp(-(0.5 / 0.47) ** 16))
    edge_code = float(sg[N // 2, 0])
    assert abs(edge_code - edge_book) < 1e-12, (edge_code, edge_book)
    # The mask never grows along a radius from the centre.
    row = sg[N // 2, N // 2:]
    assert np.all(np.diff(row) <= 0.0)
    print(f'REDUCTION super_gaussian_absorber : centre = '
          f'{sg[N // 2, N // 2]:.6f}, min inside 0.2 N = '
          f'{float(np.min(sg[inner])):.6f}, edge = {edge_code:.6f} '
          f'(book exp(-(0.5/0.47)^16) = {edge_book:.6f}), monotone = True')

    # A power of 2 or below is refused. Eq. (8.1) needs a power above 2.
    try:
        super_gaussian_absorber(16, power=2.0)
    except ValueError:
        pass
    else:
        raise AssertionError('super_gaussian_absorber must refuse power = 2')

    # 7. Partial propagation reproduces the single long step, in vacuum.
    # Six planes, so five partial propagations over the same 10 m.
    z_planes_chk = np.linspace(0.0, z_prop, 6)
    U_multi = partial_propagations(U0, lam, dx1, dx2_fixed, z_planes_chk)
    d_multi = float(np.max(np.abs(U_multi[core] - U_as[core]))) / peak
    assert d_multi < 1e-12, d_multi
    err_multi = _err(U_multi)
    assert err_multi < 1e-12, err_multi
    print(f'REDUCTION partial_propagations ({z_planes_chk.size} planes) vs '
          f'one angular_spectrum step : max rel diff = {d_multi:.3e}  '
          f'(target 1e-12); vs the closed form = {err_multi:.3e} '
          f'(target 1e-12)')

    # The step cap of Ch. 8, Eq. (8.24), printed p. 144, on this geometry.
    dz_max = min(dx1, dx2_fixed) ** 2 * N / lam
    n_min = int(np.ceil(z_prop / dz_max)) + 1
    print(f'         Ch. 8, Eq. (8.24) gives Delta_z_max = {dz_max:.3f} m, '
          f'so n >= {n_min} planes. This run used {z_planes_chk.size}.')

    # 8. The absorber leaves the beam alone. It stands well inside the flat
    # part, so the result must not move.
    U_abs = partial_propagations(U0, lam, dx1, dx2_fixed, z_planes_chk,
                                 absorber=sg)
    d_abs = float(np.max(np.abs(U_abs[core] - U_multi[core]))) / peak
    assert d_abs < 1e-9, d_abs
    print(f'REDUCTION the absorber does not touch the beam : max rel diff = '
          f'{d_abs:.3e}  (target 1e-9)')

    # 9. The refusals.
    for call, why in (
            (lambda: one_step_fresnel(U0, lam, dx1, 0.0), 'z = 0'),
            (lambda: two_step_fresnel(U0, lam, dx1, dx1, z_prop), 'm = 1'),
            (lambda: angular_spectrum(U0[:-1], lam, dx1, z_prop), 'not square'),
            (lambda: partial_propagations(U0, lam, dx1, dx2_fixed, [0.0]),
             'one plane'),
            (lambda: partial_propagations(U0, lam, dx1, dx2_fixed,
                                          [0.0, 5.0, 4.0]), 'z not sorted')):
        try:
            call()
        except ValueError:
            continue
        raise AssertionError(f'the kernel must refuse: {why}')

    print('self-check passed')
