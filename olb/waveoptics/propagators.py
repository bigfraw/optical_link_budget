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
from numpy.fft import fft2 as _np_fft2
from numpy.fft import ifft2 as _np_ifft2
from scipy import fft as _scipy_fft
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


# THE FFT BACKEND (an OPT-IN, 2026-09-06). "numpy" (the default) is the
# backend of record: every stored campaign was made with it. "scipy" runs
# the same transforms through scipy.fft with overwrite_x=True: it writes the
# result into the work array, and on the tested machine it ran a 1024 px
# complex64 transform 2.7 times faster (15.8 against 42.6 ms). The two
# backends agree at the rounding level of the field precision (6e-7 relative
# on the collected power and the SMF eta of a single-precision trial), but a
# scipy run is NOT bit-identical to a numpy run of the same seed, so a
# campaign carries the backend in its fingerprint. The setting is process
# wide: a pool worker sets it for itself (Campaign does that through the
# initializer).
FFT_BACKENDS = ("numpy", "scipy")
_fft_backend = "numpy"


def set_fft_backend(name):
    """Select the FFT backend of Forvard and Fresnel for this process.

    Args:
        name: "numpy" (the default, the backend of record) or "scipy".

    Returns:
        The previous backend name, so a caller can restore it.

    Raises:
        ValueError: the name is unknown.
    """
    global _fft_backend
    if name not in FFT_BACKENDS:
        raise ValueError(f"set_fft_backend: name must be one of "
                         f"{FFT_BACKENDS}, not {name!r}.")
    previous = _fft_backend
    _fft_backend = name
    return previous


def get_fft_backend():
    """Give the FFT backend name of this process."""
    return _fft_backend


def _fft2(a):
    if _fft_backend == "scipy":
        return _scipy_fft.fft2(a, overwrite_x=True, workers=1)
    return _np_fft2(a)


def _ifft2(a):
    if _fft_backend == "scipy":
        return _scipy_fft.ifft2(a, overwrite_x=True, workers=1)
    return _np_ifft2(a)


# THE FORVARD CACHE. The sign pattern depends on (N, dtype) only, and the
# transfer function on (N, size, lam, |z|, dtype). A split-step Monte Carlo
# makes the SAME hops in every trial, so without a cache each hop rebuilt
# four N x N arrays (three of them in double precision) that never change.
# The cache keeps them, so a hop moves only the field itself. The entries
# are read-only, and a thread that misses at the same time as another one
# computes the same array two times, which is safe. The bound is in BYTES,
# so a caller can budget it for each process: a 1024 px transfer function is
# 8 MB in single and 16 MB in double precision, so the default holds 32 or
# 16 distinct hops. A sweep past the bound drops the oldest entry. Set
# FORVARD_CACHE_BYTES before a run to change the budget.
FORVARD_CACHE_BYTES = 256 * 2 ** 20
_forvard_cache = {}


def clear_forvard_cache():
    """Drop every cached Forvard factor. A test or a memory-tight caller
    may call it."""
    _forvard_cache.clear()


def forvard_cache_bytes():
    """Give the bytes the Forvard cache holds now."""
    return sum(v[1].nbytes for k, v in _forvard_cache.items() if len(k) == 5)


def _forvard_factors(N, size, lam, z, cdtype):
    """Give the sign pattern and the transfer function of one Forvard hop.

    The arrays are the ones the body of Forvard built before the cache.
    The values are bit-identical: the phase wrap runs in double precision
    exactly as before, and the cache stores the finished factor only.

    Args:
        N:      the pixel count of one side.
        size:   the grid side, in m.
        lam:    the wavelength, in m.
        z:      the hop length, in m, not negative.
        cdtype: the complex type of the field.

    Returns:
        The pair (iiij, CC): the N x N real sign pattern and the N x N
        complex transfer function, both read-only.
    """
    cdtype = np.dtype(cdtype)
    rdtype = np.float32 if cdtype == np.complex64 else np.float64
    key = (int(N), float(size), float(lam), float(z), cdtype.str)
    hit = _forvard_cache.get(key)
    if hit is not None:
        return hit

    sign_key = (int(N), rdtype)
    iiij = _forvard_cache.get(sign_key)
    if iiij is None:
        # The alternating sign pattern does the same as a double fftshift,
        # but it is faster. See the LightPipes manual.
        iiN = np.ones((N,), dtype=rdtype)
        iiN[1::2] = -1
        iiij = np.outer(iiN, iiN)
        iiij.flags.writeable = False
        _forvard_cache[sign_key] = iiij

    # Bus = lam*z/2 * (fx^2 + fy^2). The phase of the transfer function is
    # -2*pi*Bus. Schmidt, DOI 10.1117/3.866274, Ch. 6, Eq. (6.32), printed
    # p. 95.
    _2pi = 2. * 3.141592654
    z1 = z * lam / 2
    No2 = int(N / 2)
    SW = np.arange(-No2, N - No2) / size
    SW *= SW
    SSW = SW.reshape((-1, 1)) + SW
    Bus = z1 * SSW
    # KEEP Bus IN DOUBLE PRECISION. Bus reaches thousands on a long hop, and
    # the next line takes the fractional part. In single precision that
    # subtraction loses 4 digits of the phase. So the wrap runs in double
    # precision, and only the finished factor CC takes the field precision.
    Ir = Bus.astype(int)            # truncate, do not round
    Abus = _2pi * (Ir - Bus)        # the phase, wrapped into [-2pi, 0]
    CC = (np.cos(Abus) + 1j * np.sin(Abus)).astype(cdtype)
    CC.flags.writeable = False

    while (_forvard_cache
           and forvard_cache_bytes() + CC.nbytes > FORVARD_CACHE_BYTES):
        # Drop the oldest transfer function. The sign patterns stay.
        oldest = next((k for k in _forvard_cache if len(k) == 5), None)
        if oldest is None:
            break
        del _forvard_cache[oldest]
    if CC.nbytes <= FORVARD_CACHE_BYTES:
        _forvard_cache[key] = (iiij, CC)
    return iiij, CC


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

    # THE PRECISION FOLLOWS THE FIELD. A complex64 field keeps the whole hot
    # loop in single precision: the working array, the sign pattern and the
    # transfer function. That halves the bytes each FFT moves.
    cdtype = Fin.field.dtype

    in_out = np.array(Fin.field, dtype=cdtype)     # one copy, no zero fill.

    # The legacy value of 2*pi keeps the port equal to the C++ LightPipes.
    _2pi = 2. * 3.141592654
    zz = z
    z = abs(z)
    kz = _2pi / lam * z
    cokz = np.cos(kz)
    sikz = np.sin(kz)

    # The sign pattern and the transfer function come from the cache. A
    # split-step trial makes the same hops as every other trial of the
    # plan, so the two arrays are built ONE time for each process.
    iiij, CC = _forvard_factors(N, size, lam, z, cdtype)
    in_out *= iiij

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
    # The precision follows the field. See Forvard.
    cdtype = field.dtype
    rdtype = np.float32 if cdtype == np.complex64 else np.float64

    kz = 2. * 3.141592654 / lam * z
    siz = N * dx
    dx = siz / (N - 1)              # the legacy pitch of the C++ code
    cokz = np.cos(kz)
    sikz = np.sin(kz)

    No2 = int(N / 2)

    in_outF = np.zeros((2 * N, 2 * N), dtype=cdtype)
    in_outK = np.zeros((2 * N, 2 * N), dtype=cdtype)

    # The alternating sign pattern replaces the double fftshift.
    ii2N = np.ones((2 * N), dtype=rdtype)
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

    # THE ABCD ROUTE STAYS IN DOUBLE PRECISION. The transverse phase reaches
    # millions of radians on a space link, so a single-precision r2 loses the
    # phase. The `field` setter casts the finished array to the field
    # precision, and this route is analytic, so it runs one time only.
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

    # ---- the factor cache gives the same field warm as cold ----
    clear_forvard_cache()
    FF_cold = Forvard(F0, z)
    assert forvard_cache_bytes() > 0
    FF_warm = Forvard(F0, z)
    assert np.array_equal(FF_cold.field, FF_warm.field)
    assert np.array_equal(FF_cold.field, FF.field)

    # ---- the scipy backend agrees at the rounding level, and restores ----
    assert get_fft_backend() == "numpy"
    prev = set_fft_backend("scipy")
    try:
        FS = Forvard(F0, z)
    finally:
        set_fft_backend(prev)
    assert get_fft_backend() == "numpy"
    rel = np.sqrt(np.mean(np.abs(FS.field - FF.field) ** 2)
                  / np.mean(np.abs(FF.field) ** 2))
    assert rel < 1e-12, rel
    print(f"scipy vs numpy Forvard rel rms {rel:.1e}")
    try:
        set_fft_backend("fftw")
        raise AssertionError("an unknown backend must raise")
    except ValueError:
        pass

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

    # ---- the single-precision route ----
    # A complex64 field keeps complex64 through each propagator, and it agrees
    # with the complex128 route to about 1e-6 of the peak amplitude.
    def rel_rms(a, b):
        """The rms difference of two fields, over the peak of the reference."""
        return float(np.sqrt(np.mean(np.abs(a - b) ** 2)) / np.abs(b).max())

    F0_32 = GaussBeam(Begin(size, lam, N, dtype=np.complex64), w0)
    assert F0_32.field.dtype == np.complex64
    FF32 = Forvard(F0_32, z)
    FR32 = Fresnel(F0_32, z)
    FG32 = GForvard(F0_32, z)
    assert FF32.field.dtype == np.complex64
    assert FR32.field.dtype == np.complex64
    assert FG32.field.dtype == np.complex64
    e_forvard = rel_rms(FF32.field, FF.field)
    e_fresnel = rel_rms(FR32.field, FR.field)
    e_gauss = rel_rms(FG32.field, FG.field)
    assert e_forvard < 1e-5, e_forvard
    assert e_fresnel < 1e-5, e_fresnel
    assert e_gauss < 1e-5, e_gauss
    # The single-precision Forvard also conserves the power.
    assert abs(Power(FF32) / Power(F0_32) - 1.0) < 1e-5

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
    print("")
    print("complex64 against complex128, relative rms of the field:")
    print(f"  Forvard               {e_forvard:9.2e}")
    print(f"  Fresnel               {e_fresnel:9.2e}")
    print(f"  GForvard              {e_gauss:9.2e}")
    print("self-check passed")
