'''
Digital Fourier transforms and the structure function, from Schmidt.

This module holds the two-dimensional discrete Fourier transform pair with the
physical scaling, the critical-sampling relation that ties the two grid pitches
together, and the structure-function estimator that uses that pair.

Source of every equation:
    J. D. Schmidt, "Numerical Simulation of Optical Wave Propagation with
    Examples in MATLAB", SPIE Press Monograph PM199 (2010).
    DOI: 10.1117/3.866274
Chapter 2, printed pp. 15 to 38, and Chapter 3, Sec. 3.3, printed pp. 47 to 50.
Each function names its chapter, its equation number, and its printed page.

THE CONVENTION. The book defines the continuous transform pair with the
exponent -i2*pi*f*x and no constant in front (Ch. 2, Eqs. (2.1) and (2.2),
printed pp. 15 and 16). The discrete forms below carry the sample pitch, so the
result approximates the continuous integral, not the bare sum. A bare
`numpy.fft.fft2` does NOT carry that pitch.

UNITS. `dx` is the grid pitch [m]. `df` is the frequency pitch [1/m]. The
transform of a field with units U gives units U m^2. The structure function of a
phase in radians gives radians squared.

This module holds physics only. It imports numpy only. It returns no decibels.
'''

import numpy as np
from numpy.fft import fft2 as _fft2
from numpy.fft import fftshift as _fftshift
from numpy.fft import ifft2 as _ifft2
from numpy.fft import ifftshift as _ifftshift


def freq_pitch(n, dx):
    '''
    Return the frequency-domain grid pitch of an n by n grid.

    Parameters:
        n : int
            Grid points per side.
        dx : float
            Spatial grid pitch [m].

    Returns:
        float
            Frequency grid pitch df [1/m].

    formula:
        df = 1 / (N dx)
    Source: Schmidt (2010), DOI 10.1117/3.866274, Ch. 2, text below Eq. (2.3),
    printed p. 16. Ch. 6, Eq. (6.51), printed p. 99, restates it as Delta_f1.

    VALIDITY. The relation is exact for the discrete transform. It is a
    definition, not an approximation: the grid side length is L = N dx, and the
    lowest resolved frequency is 1/L. `ft2` and `ift2` are a true inverse pair
    only when the caller passes this df to `ift2`.
    '''
    return 1.0 / (float(n) * float(dx))


def ft2(g, dx):
    '''
    Return the two-dimensional discrete Fourier transform of g, with scaling.

    Parameters:
        g : numpy.ndarray
            A square two-dimensional array. The origin sits at the grid centre.
        dx : float
            Spatial grid pitch [m].

    Returns:
        numpy.ndarray
            The transform G, with the zero frequency at the grid centre.

    formula:
        G(fx,fy) = fftshift(fft2(ifftshift(g))) dx^2
    Source: Schmidt (2010), DOI 10.1117/3.866274, Ch. 2, Eq. (2.6), printed
    p. 16, gives the one-dimensional sum. Eq. (2.32), printed p. 36, gives the
    two-dimensional continuous pair, and Sec. 2.6, printed p. 36, gives the
    two-dimensional scaling dx^2.

    The shift pair moves the origin to the first sample before the transform,
    and moves it back after (Ch. 2, Eq. (2.5) and Fig. 2.1, printed pp. 16 and
    17). The pitch dx^2 turns the sum into the Riemann sum of Ch. 2, Eq. (2.3),
    printed p. 15.

    VALIDITY. The book restricts the discussion to an even grid count N (Ch. 2,
    Sec. 2.1.3, printed p. 18) and to the same count and pitch in x and y
    (Ch. 2, Sec. 2.6, printed p. 36). This function keeps the square-grid
    assumption. The `ifftshift` and `fftshift` pair above is correct for an odd
    count too.

    The result is a sampled, rippled, and aliased copy of the continuous
    transform (Ch. 2, Eqs. (2.39) and (2.40), printed p. 36). A smaller dx cuts
    the aliasing. A larger grid side L = N dx cuts the ripple (Ch. 2, Sec. 2.4,
    printed p. 26).
    '''
    g = np.asarray(g)
    return _fftshift(_fft2(_ifftshift(g))) * float(dx) ** 2


def ift2(G, df):
    '''
    Return the two-dimensional discrete inverse Fourier transform of G.

    Parameters:
        G : numpy.ndarray
            A square two-dimensional array. The zero frequency sits at the grid
            centre.
        df : float
            Frequency grid pitch [1/m]. Use `freq_pitch(N, dx)`.

    Returns:
        numpy.ndarray
            The spatial function g, with the origin at the grid centre.

    formula:
        g(x,y) = fftshift(ifft2(ifftshift(G))) (N df)^2
    Source: Schmidt (2010), DOI 10.1117/3.866274, Ch. 2, Eq. (2.9), printed
    p. 17, gives the one-dimensional sum. Eq. (2.33), printed p. 36, gives the
    two-dimensional continuous pair, and Sec. 2.6, printed p. 37, gives the
    two-dimensional scaling (N df)^2.

    The factor N^2 cancels the 1/N^2 that `ifft2` applies. The factor df^2 turns
    the sum into the Riemann sum of Ch. 2, Eq. (2.7), printed p. 17.

    VALIDITY. Same as `ft2`. The function reads the grid count from the first
    axis, so the grid must be square. With df = 1/(N dx) the pair `ft2` then
    `ift2` is an exact identity, because dx^2 (N df)^2 = 1.
    '''
    G = np.asarray(G)
    n = G.shape[0]
    return _fftshift(_ifft2(_ifftshift(G))) * (n * float(df)) ** 2


def structure_function(ph, mask, dx):
    '''
    Return the structure function of one realisation of a windowed field.

    The structure function measures the mean square difference between two
    points of a random field, as a function of their separation. Use it to
    verify a phase screen against the Kolmogorov law.

    Parameters:
        ph : numpy.ndarray
            The field, for example a phase screen [rad]. Square and real.
        mask : numpy.ndarray
            The window w(r). It is 1 inside the pupil and 0 outside. Give an
            array of ones to use the whole grid.
        dx : float
            Spatial grid pitch [m].

    Returns:
        numpy.ndarray
            D(dr), the structure function against the separation dr. It has the
            same shape as `ph`, and the zero separation sits at the grid centre.
            The result is zero where the window overlap area A(dr) is zero.

    formula:
        W(f) = FT{w(r)},   U(f) = FT{u'(r)},   S(f) = FT{u'(r)^2}
        A(dr)      = IFT{ |W(f)|^2 }
        D'(dr)     = 2 IFT{ Re[S(f) W*(f)] - |U(f)|^2 }
        D(dr)      = D'(dr) / A(dr)
    Source: Schmidt (2010), DOI 10.1117/3.866274, Ch. 3, Sec. 3.3. The
    definition is Eq. (3.15), printed p. 47. The window is Eq. (3.6), printed
    p. 45, and the windowed field is Eq. (3.7), printed p. 45. The overlap area
    A(dr) is Eq. (3.10), printed p. 45. The window normalisation is Eq. (3.17),
    printed p. 48, and the transform form is Eqs. (3.19) to (3.25), printed
    pp. 49 and 50.

    THE TRANSFORM DIRECTION. The book prints a FORWARD transform at Eq. (3.25),
    printed p. 50. This function uses the inverse transform, as Listing 3.7,
    printed p. 48, does. The two give the same result here, because the bracket
    is real and even in f. The inverse transform also puts the separation
    coordinate on the same grid as the input.

    VALIDITY.
    - The result is ONE realisation. The statistical structure function is the
      ensemble mean over many independent draws (Ch. 3, text below Eq. (3.15),
      printed p. 47).
    - The discrete transform makes the correlation CIRCULAR. So the window must
      leave a guard band at the grid edge. A window of ones over the whole grid
      gives a wrapped, incorrect result for a field that does not go to zero at
      the edge.
    - The result is valid only where A(dr) > 0, that is out to the window
      diameter. The caller selects the separation range.
    - Eq. (3.16), printed p. 48, relates D to the auto-correlation for a
      statistically isotropic field only.
    '''
    ph = np.asarray(ph, dtype=float)
    mask = np.asarray(mask, dtype=float)
    if ph.shape != mask.shape:
        raise ValueError('structure_function needs ph and mask of one shape; '
                         'see Schmidt (2010), DOI 10.1117/3.866274, Ch. 3, '
                         'Eq. (3.7), printed p. 45')
    windowed = ph * mask
    df = freq_pitch(ph.shape[0], dx)

    # Eqs. (3.19) to (3.21), printed p. 49.
    W = ft2(mask, dx)
    U = ft2(windowed, dx)
    S = ft2(windowed ** 2, dx)

    # Eq. (3.10), printed p. 45: the window overlap area A(dr).
    area = np.real(ift2(W * np.conj(W), df))
    # Eq. (3.25), printed p. 50: the windowed structure function D'(dr).
    numerator = 2.0 * np.real(
        ift2(np.real(S * np.conj(W)) - np.abs(U) ** 2, df))

    # Eq. (3.17), printed p. 48: divide out the overlap area.
    out = np.zeros_like(ph)
    good = area > 1e-12 * np.max(np.abs(area))
    out[good] = numerator[good] / area[good]
    return out


if __name__ == '__main__':
    # ---------------- physics self-checks ----------------
    # 1. The Gaussian case of Ch. 2, Sec. 2.5.2, printed pp. 31 and 32.
    # The pair, in two dimensions through Ch. 2, Eq. (2.32), printed p. 36:
    #     g(x,y) = exp(-a^2 (x^2 + y^2))
    #     G(fx,fy) = (pi/a^2) exp(-pi^2 (fx^2 + fy^2) / a^2)
    # The constants follow the transform convention of Ch. 2, Eq. (2.1),
    # printed p. 15. The book prints the one-dimensional pair at Eqs. (2.23)
    # and (2.24).
    a = 1.0
    N = 256
    dx = 12.0 / N                     # a grid side of 12 m, so the tails vanish
    df = freq_pitch(N, dx)
    axis = (np.arange(N) - N // 2) * dx
    X, Y = np.meshgrid(axis, axis)
    g = np.exp(-(a ** 2) * (X ** 2 + Y ** 2))

    faxis = (np.arange(N) - N // 2) * df
    FX, FY = np.meshgrid(faxis, faxis)
    G_exact = (np.pi / a ** 2) * np.exp(
        -(np.pi ** 2) * (FX ** 2 + FY ** 2) / a ** 2)

    G_dft = ft2(g, dx)
    assert np.max(np.abs(np.imag(G_dft))) < 1e-12, 'the transform must be real'

    # Compare over the band that carries the signal. Ch. 2, Eqs. (2.26) and
    # (2.27), printed p. 33, set the band edge at the frequency where the
    # spectrum falls to a fraction p of its peak. For the pair above that
    # frequency is f_p = (a/pi) sqrt(-ln p). Take p = 1e-3.
    f_p = (a / np.pi) * np.sqrt(-np.log(1e-3))
    band = np.hypot(FX, FY) <= f_p
    err_gauss = float(np.max(np.abs(
        np.real(G_dft[band]) / G_exact[band] - 1.0)))
    assert err_gauss < 1e-6, err_gauss
    print(f'REDUCTION ft2 of a Gaussian vs Ch. 2, Eq. (2.24) : '
          f'max rel err = {err_gauss:.3e}  (target 1e-6, band |f| <= '
          f'{f_p:.3f} 1/m)')

    # 2. The round trip is an identity, because dx^2 (N df)^2 = 1.
    back = ift2(ft2(g, dx), df)
    err_trip = float(np.max(np.abs(back - g)))
    assert err_trip < 1e-12, err_trip
    print(f'REDUCTION ift2(ft2(g)) - g : max abs err = {err_trip:.3e}  '
          f'(target 1e-12)')

    # The round trip holds for a complex field too.
    rng = np.random.default_rng(0)
    field = rng.standard_normal((N, N)) + 1j * rng.standard_normal((N, N))
    err_cplx = float(np.max(np.abs(ift2(ft2(field, dx), df) - field)))
    assert err_cplx < 1e-10, err_cplx
    print(f'REDUCTION ift2(ft2(complex noise)) : max abs err = '
          f'{err_cplx:.3e}  (target 1e-10)')

    # 3. The structure function of a linear phase ramp phi = a_ramp x.
    # The difference phi(r) - phi(r + dr) is -a_ramp dx_sep everywhere, so
    # D(dr) = a_ramp^2 dx_sep^2 exactly, at every separation.
    a_ramp = 3.0                                   # [rad/m]
    ramp = a_ramp * X
    radius = 1.5                                   # a pupil of 1.5 m radius
    pupil = (np.hypot(X, Y) <= radius).astype(float)
    D = structure_function(ramp, pupil, dx)

    # Read the x axis of D, out to the pupil diameter. Beyond that the overlap
    # area falls to zero and the estimate is not defined.
    row = D[N // 2, :]
    sep = axis
    inside = np.abs(sep) <= 1.8 * radius
    D_exact = (a_ramp * sep) ** 2
    err_sf = float(np.max(np.abs(row[inside] - D_exact[inside]))
                   / np.max(D_exact[inside]))
    assert err_sf < 1e-8, err_sf
    print(f'REDUCTION structure_function of phi = {a_ramp:.1f} x vs '
          f'D = a^2 r^2 : max rel err = {err_sf:.3e}  (target 1e-8, '
          f'|dr| <= {1.8 * radius:.2f} m)')

    # The structure function is zero at zero separation.
    assert abs(D[N // 2, N // 2]) < 1e-9, D[N // 2, N // 2]
    probe = N // 2 + int(round(1.0 / dx))
    print(f'D(0) = {D[N // 2, N // 2]:.3e} rad^2   '
          f'D(r = {sep[probe]:.4f} m) along x = {row[probe]:.4f} rad^2 '
          f'(exact {(a_ramp * sep[probe]) ** 2:.4f})')

    # 4. A constant phase gives a structure function of zero everywhere.
    D_flat = structure_function(np.full((N, N), 2.7), pupil, dx)
    err_flat = float(np.max(np.abs(D_flat)))
    assert err_flat < 1e-9, err_flat
    print(f'REDUCTION structure_function of a constant phase : '
          f'max abs = {err_flat:.3e} rad^2  (target 1e-9)')

    # 5. The critical-sampling relation ties the two pitches.
    assert abs(freq_pitch(N, dx) * (N * dx) - 1.0) < 1e-15

    # 6. A mask of the wrong shape is refused.
    try:
        structure_function(ramp, pupil[:-2, :-2], dx)
    except ValueError:
        pass
    else:
        raise AssertionError('structure_function must refuse a bad mask shape')

    print('self-check passed')
