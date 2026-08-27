"""The free-space propagators: Forvard, Fresnel and GForvard.

Ported and trimmed from LightPipes (https://github.com/opticspy/lightpipes),
BSD-3-Clause. See LIGHTPIPES_LICENSE.txt in this package.

The module is pure physics. It imports numpy and scipy only. It gives three
routes from one plane to another plane:

- Forvard:  the FFT spectral method. The grid keeps its side. The method is
            periodic, so energy that leaves one edge comes back at the
            opposite edge.
- Fresnel:  the convolution method on a doubled grid. The doubled grid
            absorbs the wrap of the spectral method.
- GForvard: the analytic ABCD route. It is exact, but it accepts a pure
            fundamental Gaussian beam only.

The port keeps the numpy FFT branch of LightPipes. It drops the pyFFTW
branch.

Sources:
- Schmidt, Numerical Simulation of Optical Wave Propagation with Examples
  in MATLAB, DOI 10.1117/3.866274. The angular-spectrum form is Ch. 6,
  Eqs. (6.31) and (6.32), printed p. 95. The Fresnel convolution form is
  Ch. 6, Eq. (6.6), printed p. 88. The sampling limits are Ch. 7,
  Eqs. (7.41), (7.42) and (7.59), printed pp. 123 and 127.
- Goodman, Introduction to Fourier Optics, ISBN 978-0974707723. The Fresnel
  diffraction integral.
- Siegman, Lasers, ISBN 978-0935702118. The ABCD law for the complex beam
  parameter q.
- LightPipes manual, https://opticspy.github.io/lightpipes/manual.html.
  The implementation lineage.
"""

import numpy as np
from numpy.fft import fft2 as _fft2
from numpy.fft import ifft2 as _ifft2
from scipy.special import fresnel as _fresnel

from .field import Field


def _reject_spherical(Fin, name):
    """Refuse a field that is in spherical (co-moving) coordinates.

    A flat-grid propagator cannot read a field that carries a curvature,
    because the grid side is a function of the coordinate system. The
    lens propagators in lenses.py set that curvature. Convert() removes it.
    See the LightPipes manual,
    https://opticspy.github.io/lightpipes/manual.html.

    LightPipes prints a message and gives the field back unchanged. This
    port raises, because a silent bad result is worse than a stop.
    """
    if Fin._curvature != 0.0:
        raise ValueError(f'{name}: the field is in spherical coordinates. '
                         'Use Convert() first.')


def Forvard(Fin, z):
    """Propagate the field with the FFT spectral method.

    The method multiplies the angular spectrum with the transfer function

        H(fx,fy) = exp(i*k*z) * exp(-i*pi*lam*z*(fx^2 + fy^2))

    That expression is exactly Schmidt (2010), DOI 10.1117/3.866274, Ch. 6,
    Eq. (6.32), printed p. 95, and the two-transform chain is Eq. (6.31) on
    the same page. See also Goodman, ISBN 978-0974707723. NOTE: this function
    KEEPS the piston factor exp(i*k*z). The book's Listings 6.1, 6.3 and 6.5
    (printed pp. 91, 96 and 102) drop it, and so does
    olb.waveoptics.schmidt.fresnel.angular_spectrum. The irradiance is the
    same; a phase comparison must add or remove the factor.

    The range limit of the method is constraint 4, Ch. 7, Eq. (7.59), printed
    p. 127. olb.waveoptics.grid.forvard_max_z gives the number.

    The grid keeps its side and its pitch. The method is periodic. A beam
    that becomes wider than the grid wraps around the edges. Give the grid
    a side of about 8 times the largest beam radius. THAT SIDE RULE IS A
    RULE OF THUMB WITH NO BOOK SOURCE. The book's bound is constraints 1 and
    2, Ch. 7, Eqs. (7.14) and (7.20), printed pp. 119 and 120, and both need
    the source extent D1, the observation extent D2 and the range z. See
    docs/schmidt-crosscheck.md, gap S-16.

    Args:
        Fin: the input field.
        z:   the propagation distance, in m. A negative z propagates back.

    Returns:
        A new Field.

    Raises:
        ValueError: the field is in spherical coordinates.
    """
    _reject_spherical(Fin, 'Forvard')
    if z == 0:
        return Field.copy(Fin)

    Fout = Field.shallowcopy(Fin)
    N = Fout.N
    size = Fout.siz
    lam = Fout.lam

    in_out = np.zeros((N, N), dtype=np.complex128)
    in_out[:, :] = Fin.field

    # The legacy value of 2*pi keeps the port equal to the C++ LightPipes.
    _2pi = 2. * 3.141592654
    zz = z
    z = abs(z)
    kz = _2pi / lam * z
    cokz = np.cos(kz)
    sikz = np.sin(kz)

    # The alternating sign pattern does the same as a double fftshift,
    # but it is faster. See the LightPipes manual.
    iiN = np.ones((N,), dtype=float)
    iiN[1::2] = -1
    iiij = np.outer(iiN, iiN)
    in_out *= iiij

    # Bus = lam*z/2 * (fx^2 + fy^2). The phase of the transfer function is
    # -2*pi*Bus. Schmidt, DOI 10.1117/3.866274, Ch. 6, Eq. (6.32), printed
    # p. 95.
    z1 = z * lam / 2
    No2 = int(N / 2)
    SW = np.arange(-No2, N - No2) / size
    SW *= SW
    SSW = SW.reshape((-1, 1)) + SW
    Bus = z1 * SSW
    Ir = Bus.astype(int)            # truncate, do not round
    Abus = _2pi * (Ir - Bus)        # the phase, wrapped into [-2pi, 0]
    CC = np.cos(Abus) + 1j * np.sin(Abus)

    if zz >= 0.0:
        in_out = _fft2(in_out)
        in_out *= CC
        in_out = _ifft2(in_out)
    else:
        in_out = _ifft2(in_out)
        in_out *= CC.conjugate()
        in_out = _fft2(in_out)

    in_out *= (cokz + 1j * sikz)
    in_out *= iiij                  # numpy normalises the ifft already
    Fout.field = in_out
    Fout._IsGauss = False
    return Fout


def Fresnel(Fin, z):
    """Propagate the field with the convolution method.

    The method convolves the field with the Fresnel kernel on a grid of
    twice the side. That is the CONVOLUTION form of the Fresnel integral,
    Schmidt (2010), DOI 10.1117/3.866274, Ch. 6, Eq. (6.6), printed p. 88.
    See also Goodman, ISBN 978-0974707723.

    THE PIXEL-INTEGRATED KERNEL IS NOT THE BOOK'S. The kernel integral over
    one pixel has a closed form in the Fresnel integrals C(x) and S(x), and
    this port uses it. The book multiplies by the analytic transfer function
    of Ch. 6, Eq. (6.49), printed p. 99, instead. The two solve the same
    Eq. (6.6). The pixel integral is a refinement of the book, not a
    departure from it.

    The doubled grid absorbs the periodic wrap of the spectral method.
    The method needs a field that is zero at the edges of the grid. THE BOOK
    HAS NO ZERO-PADDED CONVOLUTION: it controls the wrap with the absorbing
    boundary of Ch. 8, Eq. (8.1), printed p. 134, and with the grid-size rule
    of Ch. 7, Eq. (7.20), printed p. 120. Two cures for one problem.

    The method has a MINIMUM distance. The convolution does not give a
    valid result when z is comparable with, or less than, the size of the
    aperture that diffracts the field. Use Forvard for a short hop. See the
    LightPipes manual, https://opticspy.github.io/lightpipes/manual.html.
    Schmidt, DOI 10.1117/3.866274, Ch. 7, Eqs. (7.41) and (7.42), printed
    p. 123, give the number: z >= D1 dx1 R / (lambda R - D1 dx1) for a source
    of the wavefront radius R, and z >= D1 dx1 / lambda for a flat source.
    This function states the rule in words only and it does not check it.
    olb.waveoptics.schmidt.sampling.fresnel_min_distance gives the bound.

    Args:
        Fin: the input field.
        z:   the propagation distance, in m. It must not be negative.

    Returns:
        A new Field.

    Raises:
        ValueError: z is negative, or the field is in spherical coordinates.
    """
    _reject_spherical(Fin, 'Fresnel')
    if z < 0:
        raise ValueError('Fresnel does not support negative z')
    if z == 0:
        return Field.copy(Fin)
    Fout = Field.shallowcopy(Fin)
    Fout.field = _field_Fresnel(z, Fin.field, Fin.dx, Fin.lam)
    Fout._IsGauss = False
    return Fout


def _field_Fresnel(z, field, dx, lam):
    """Do the Fresnel convolution on the raw array.

    The port keeps the legacy pixel pitch dx = siz/(N-1) of the C++
    LightPipes, so the numbers match the reference package.
    """
    N = field.shape[0]

    kz = 2. * 3.141592654 / lam * z
    siz = N * dx
    dx = siz / (N - 1)              # the legacy pitch of the C++ code
    cokz = np.cos(kz)
    sikz = np.sin(kz)

    No2 = int(N / 2)

    in_outF = np.zeros((2 * N, 2 * N), dtype=np.complex128)
    in_outK = np.zeros((2 * N, 2 * N), dtype=np.complex128)

    # The alternating sign pattern replaces the double fftshift.
    ii2N = np.ones((2 * N), dtype=float)
    ii2N[1::2] = -1
    iiij2N = np.outer(ii2N, ii2N)
    iiij2No2 = iiij2N[:2 * No2, :2 * No2]
    iiijN = iiij2N[:N, :N]

    # The kernel pixel integral. C(x) and S(x) are the Fresnel integrals.
    # Goodman, ISBN 978-0974707723. The kernel itself is Schmidt,
    # DOI 10.1117/3.866274, Ch. 6, Eq. (6.6), printed p. 88. The book does
    # not integrate it over a pixel; it transforms it, Eq. (6.49), printed
    # p. 99.
    RR = np.sqrt(1 / (2 * lam * z)) * dx * 2
    io = np.arange(0, (2 * No2) + 1)    # one extra sample to stride
    R1 = RR * (io - No2)
    fs, fc = _fresnel(R1)
    fss = np.outer(fs, fs)
    fsc = np.outer(fs, fc)
    fcs = np.outer(fc, fs)
    fcc = np.outer(fc, fc)

    temp_re = (fsc[1:, 1:] + fcs[1:, 1:])
    temp_re -= fsc[:-1, 1:]
    temp_re -= fcs[:-1, 1:]
    temp_re -= fsc[1:, :-1]
    temp_re -= fcs[1:, :-1]
    temp_re += fsc[:-1, :-1]
    temp_re += fcs[:-1, :-1]

    temp_im = (-fcc[1:, 1:] + fss[1:, 1:])
    temp_im += fcc[:-1, 1:]
    temp_im -= fss[:-1, 1:]
    temp_im += fcc[1:, :-1]
    temp_im -= fss[1:, :-1]
    temp_im -= fcc[:-1, :-1]
    temp_im += fss[:-1, :-1]

    temp_K = 1j * temp_im
    temp_K += temp_re
    temp_K *= iiij2No2
    temp_K *= 0.5
    in_outK[(N - No2):(N + No2), (N - No2):(N + No2)] = temp_K

    in_outF[(N - No2):(N + No2), (N - No2):(N + No2)] \
        = field[(N - 2 * No2):N, (N - 2 * No2):N]   # cut the field if N is odd
    in_outF[(N - No2):(N + No2), (N - No2):(N + No2)] *= iiij2No2

    in_outK = _fft2(in_outK)
    in_outF = _fft2(in_outF)
    in_outF *= in_outK
    in_outF *= iiij2N
    in_outF = _ifft2(in_outF)

    Ftemp = (in_outF[No2:N + No2, No2:N + No2]
             - in_outF[No2 - 1:N + No2 - 1, No2:N + No2])
    Ftemp += in_outF[No2 - 1:N + No2 - 1, No2 - 1:N + No2 - 1]
    Ftemp -= in_outF[No2:N + No2, No2 - 1:N + No2 - 1]
    Ftemp *= 0.25 * complex(cokz, sikz)
    Ftemp *= iiijN
    return Ftemp


def GForvard(Fin, z):
    """Propagate a pure Gaussian beam with the ABCD matrix.

    The route is analytic, so it has no grid artefact. It accepts a field
    from GaussBeam() only. Each mask or each FFT propagator clears the
    Gaussian flag.

    The ABCD law gives the new complex beam parameter:

        q_out = (A*q_in + B) / (C*q_in + D),  with [A B; C D] = [1 z; 0 1]

    Then w^2 = -lam/pi * (Im(q) + Re(q)^2/Im(q)), and 1/R = Re(1/q).
    See Siegman, Lasers, ISBN 978-0935702118.

    Schmidt (2010), DOI 10.1117/3.866274, Ch. 6, Sec. 6.5, gives the ray
    matrices, Eq. (6.70), printed p. 103, the thin-lens matrix, Eq. (6.76),
    printed p. 104, and the generalized Huygens-Fresnel (ABCD) integral,
    Eq. (6.77), printed p. 104. It does NOT give the closed-form q transform
    above, so the source of that step is Siegman. Eq. (6.77) holds for an
    azimuthally symmetric field, which a pure Gaussian is.

    Args:
        Fin: the input field. It must come from GaussBeam().
        z:   the propagation distance, in m.

    Returns:
        A new Field.

    Raises:
        ValueError: the input field is not a pure Gaussian beam, or the
            field is in spherical coordinates.
    """
    _reject_spherical(Fin, 'GForvard')
    return _ABCD(Fin, [[1.0, z], [0.0, 1.0]])


def _ABCD(Fin, M):
    """Apply one ABCD matrix to a pure Gaussian beam."""
    Fout = Field.copy(Fin)
    A, B = M[0][0], M[0][1]
    C, D = M[1][0], M[1][1]
    if not Fin._IsGauss:
        raise ValueError('GForvard: the field is not a pure Gaussian beam')

    Fout._q = (A * Fin._q + B) / (C * Fin._q + D)
    Fout._z = Fin._z + B
    w2 = -Fin.lam / np.pi * (Fout._q.imag
                             + Fout._q.real * Fout._q.real / Fout._q.imag)
    w02 = Fin._w0 * Fin._w0
    w = np.sqrt(w2)
    inv_R = (1 / Fout._q).real

    # The Gouy phase of the fundamental mode is arctan(z/z_R).
    # Siegman, ISBN 978-0935702118.
    z0 = np.pi * w02 / Fin.lam
    k = 2 * np.pi / Fin.lam
    phase_z = k * Fout._z - np.arctan(Fout._z / z0)

    r2 = Fin.mgrid_Rsquared
    phase_trans = k / 2 * inv_R * r2
    w0w = Fin._w0 / w
    Fout.field = (Fin._A * w0w * np.exp(-r2 / w2)
                  * np.exp(1j * (phase_trans + phase_z)))
    Fout._IsGauss = True
    Fout._w0 = Fin._w0
    Fout._A = Fin._A
    return Fout


if __name__ == '__main__':
    from .field import Begin, Power
    from .sources import CircAperture, GaussBeam

    lam = 1550e-9
    w0 = 5e-3
    z = 200.0
    N = 512
    zR = np.pi * w0**2 / lam
    wz = w0 * np.sqrt(1 + (z / zR)**2)       # the analytic radius at z

    def bucket_power(F, R):
        """The power inside a circle of the radius R, normalised."""
        return Power(CircAperture(F, R)) / Power(F)

    # ---- GForvard: the analytic radius ----
    size = 8 * wz
    F0 = GaussBeam(Begin(size, lam, N), w0)
    FG = GForvard(F0, z)
    # Read w back from two amplitudes on the x axis: the Gaussian gives
    # w^2 = r^2 / ln(|E(0)|/|E(r)|).
    c = N // 2
    r = FG.xvalues[c + N // 8]
    a0 = abs(FG.field[c, c])
    ar = abs(FG.field[c, c + N // 8])
    w_read = np.sqrt(r * r / np.log(a0 / ar))
    assert abs(w_read - wz) / wz < 1e-9
    # The analytic route conserves power.
    assert abs(Power(FG) / Power(F0) - 1.0) < 1e-9

    # ---- Forvard conserves power ----
    FF = Forvard(F0, z)
    assert abs(Power(FF) / Power(F0) - 1.0) < 1e-12

    # ---- the three routes agree on a well sampled grid ----
    FR = Fresnel(F0, z)
    b_g = bucket_power(FG, wz)
    b_f = bucket_power(FF, wz)
    b_r = bucket_power(FR, wz)
    assert abs(b_f - b_g) / b_g < 1e-3, (b_f, b_g)
    assert abs(b_r - b_g) / b_g < 1e-3, (b_r, b_g)

    # ---- the documented failure mode of Forvard ----
    # A grid of only 2 times the final radius wraps the beam at the edges.
    small = 2 * wz
    F0s = GaussBeam(Begin(small, lam, N), w0)
    b_gs = bucket_power(GForvard(F0s, z), wz)
    b_fs = bucket_power(Forvard(F0s, z), wz)
    assert abs(b_fs - b_gs) / b_gs > 1e-2, (b_fs, b_gs)

    # ---- the spherical-coordinate guard ----
    Fsph = Field.copy(F0)
    Fsph._curvature = -1e-6
    for name, call in (("Forvard", lambda: Forvard(Fsph, z)),
                       ("Fresnel", lambda: Fresnel(Fsph, z)),
                       ("GForvard", lambda: GForvard(Fsph, z))):
        try:
            call()
            raise AssertionError(f"{name} must refuse a spherical field")
        except ValueError as exc:
            assert 'Convert' in str(exc), name

    print(f"wavelength              {lam * 1e9:9.1f} nm")
    print(f"waist radius w0         {w0 * 1e3:9.3f} mm")
    print(f"distance z              {z:9.1f} m")
    print(f"Rayleigh range zR       {zR:9.3f} m")
    print(f"analytic w(z)           {wz * 1e3:9.4f} mm")
    print(f"w(z) read from GForvard {w_read * 1e3:9.4f} mm")
    print("")
    print("bucket power in radius w(z), grid = 8 w(z):")
    print(f"  GForvard (analytic)   {b_g:9.6f}")
    print(f"  Forvard  (FFT)        {b_f:9.6f}")
    print(f"  Fresnel  (convol.)    {b_r:9.6f}")
    print("")
    print("bucket power, grid = 2 w(z), the periodic artefact:")
    print(f"  GForvard (analytic)   {b_gs:9.6f}")
    print(f"  Forvard  (FFT)        {b_fs:9.6f}")
    print(f"  relative difference   {abs(b_fs - b_gs) / b_gs:9.4f}")
    print("self-check passed")
