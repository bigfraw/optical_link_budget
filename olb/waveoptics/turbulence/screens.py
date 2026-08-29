"""The turbulent phase screens, and how to put one into a field.

The module makes one Kolmogorov (or von Karman) phase screen, and it applies
that screen to a Field. It is the first half of the split-step layer. The
second half, splitstep.py, moves the field between the screens.

A screen is a real array in RADIANS at the pixel pitch of the propagation
grid. The module makes the screen AT that pitch, on purpose. A coarse screen
that an FFT interpolates up to the grid is band-limited above the coarse
Nyquist frequency. Such a screen loses the structure near the Fresnel scale
sqrt(lambda*z), and that structure builds the scintillation. Then the result
becomes a function of the coarse grid, not of the atmosphere. See Schmidt,
DOI 10.1117/3.866274, Ch. 9.

The random draw comes from aotools. The package is a DEPENDENCY. This module
imports it. It does not copy it, because aotools is LGPL-3.0.

Sources:
- Fried, Optical resolution through a randomly inhomogeneous medium,
  DOI 10.1364/JOSA.56.001372. The coherence length r0 and the phase
  structure function D_phi(r) = 6.88 (r/r0)^(5/3).
- Andrews and Phillips, Laser Beam Propagation through Random Media, 2nd ed.,
  DOI 10.1117/3.626196, Ch. 12, Eq. (23). The same r0 for a path integral.
- Schmidt, Numerical Simulation of Optical Wave Propagation with Examples in
  MATLAB, DOI 10.1117/3.866274, Ch. 9. The modified von Karman phase PSD is
  Eq. (9.51), printed p. 161, and its ordinary-frequency form is Eq. (9.52) on
  the same page. The Fourier-series screen is Eqs. (9.78) to (9.80), printed
  pp. 166 and 167, with Listing 9.2, printed p. 168. The subharmonic screen is
  Eq. (9.81), printed p. 169, with Listing 9.3, printed p. 170. The per-screen
  Fried parameter is Eq. (9.70), printed p. 165. The pitch rules are Sec. 9.4,
  printed p. 172.
- Lane, Glindemann and Dainty, Simulation of a Kolmogorov phase screen,
  DOI 10.1088/0959-7174/2/3/003. The subharmonic method that Eq. (9.81) uses.
- Johansson and Gavel, Simulation of stellar speckle imaging,
  DOI 10.1117/12.177254. The subharmonic set that the book calls the closest
  match to theory (Ch. 9, text above Sec. 9.4, printed p. 172).
"""

import numpy as np

from ..field import Field


def _load_aotools():
    """Import aotools lazily, with a helpful error when it is absent.

    aotools is LGPL-3.0. This package IMPORTS it as a dependency. It does
    not copy the code into olb.

    Returns:
        The aotools module.

    Raises:
        ImportError: aotools is not installed.
    """
    try:
        import aotools
        return aotools
    except ImportError as e:
        raise ImportError(
            "the turbulent split-step layer needs the `aotools` package. "
            "Run `pip install aotools`, or `pip install olb[screens]`."
        ) from e


def screen_r0(cn2_integral_m13, wavelength_m):
    """Calculate the Fried parameter of one slab of the path.

        r0 = (0.423 * k^2 * INT Cn2 dz)^(-3/5),  k = 2*pi/lambda

    See Fried, DOI 10.1364/JOSA.56.001372, and Andrews and Phillips,
    DOI 10.1117/3.626196, Ch. 12, Eq. (23). The same expression, the constant
    included, is Schmidt (2010), DOI 10.1117/3.866274, Ch. 9, Eq. (9.70),
    printed p. 165. It is the PLANE-wave r0 of one layer, so the layer must be
    thin.

    The caller gives the INTEGRAL of Cn2 over the slab, in m^(1/3). That
    integral carries any slant factor already. For a slant path the caller
    multiplies the vertical integral by the airmass sec(zenith). This
    function adds no geometry.

    Args:
        cn2_integral_m13: the integrated Cn2 of the slab, in m^(1/3). It is
            a scalar or an ndarray.
        wavelength_m:     the wavelength, in m.

    Returns:
        The Fried parameter, in m. The type follows the input.
    """
    k = 2.0 * np.pi / wavelength_m
    return (0.423 * k * k * np.asarray(cn2_integral_m13, float)) ** (-3.0 / 5.0)


def phase_screen(r0_m, n, pixel_m, L0_m=np.inf, l0_m=1e-6, seed=None,
                 subharmonics=True):
    """Make one random phase screen at the pitch of the propagation grid.

    The screen holds the modified von Karman phase spectrum

        PHI(f) = 0.023 r0^(-5/3) exp(-(f/fm)^2) / (f^2 + f0^2)^(11/6)

    with fm = 5.92/(2*pi*l0) and f0 = 1/L0. That is Schmidt (2010),
    DOI 10.1117/3.866274, Ch. 9, Eq. (9.51), printed p. 161, written in
    ORDINARY frequency, as Eq. (9.52) on the same page does. The angular
    constant converts: 0.49 (2 pi)^(-5/3) = 0.02290, which is 0.42% from the
    printed 0.023. The default scales give the pure Kolmogorov spectrum: L0
    is infinite and l0 is 1 um.

    The screen has the SAME pitch and the SAME pixel count as the
    propagation grid. Do not make a coarse screen and interpolate it up. A
    coarse screen carries no power above its own Nyquist frequency. It
    misses the structure at the Fresnel scale sqrt(lambda*z), which is the
    scale that builds the scintillation. The result then follows the coarse
    grid, not the atmosphere. See Schmidt, DOI 10.1117/3.866274, Sec. 9.4,
    printed p. 172.

    The Fourier screen is the Fourier-series draw of Schmidt,
    DOI 10.1117/3.866274, Ch. 9, Eqs. (9.78) to (9.80), printed pp. 166 and
    167, with Listing 9.2, printed p. 168. It is band-limited in the two
    directions. It holds no power above the grid Nyquist frequency, and it
    holds too little power below 1/(n*pixel). The subharmonics of Eq. (9.81),
    printed p. 169, add three levels of low frequency back. That method comes
    from Lane, Glindemann and Dainty, DOI 10.1088/0959-7174/2/3/003. The
    subharmonics lift the structure function, but they do not close the gap.
    The self-check measures the residual deficit.

    THE aotools SUBHARMONIC SCREEN IS NOT THE BOOK'S GENERATOR. Against
    Eq. (9.44), printed p. 160, the two agree well in the band r/r0 = 0.3 to
    1.6: the book generator reaches 0.88 to 0.93 of theory there, and aotools
    reaches 1 to 3% above it. Both fall away at a larger separation. The book
    calls the subharmonic set of Johansson and Gavel, DOI 10.1117/12.177254,
    the closest match to theory (Ch. 9, text above Sec. 9.4, printed p. 172).
    See docs/schmidt-crosscheck.md, gap S-27, and
    examples/schmidt/screens_and_turbulence.py.

    aotools takes an INTEGER seed. It builds numpy.random.default_rng(seed)
    itself. A caller that runs many screens bridges from its own
    numpy.random.SeedSequence: draw the integers with
    SeedSequence(...).generate_state(count) and pass one integer per screen.

    Args:
        r0_m:         the Fried parameter of the slab, in m.
        n:            the number of pixels along one side.
        pixel_m:      the pixel pitch of the propagation grid, in m.
        L0_m:         the outer scale, in m. The default is infinite.
        l0_m:         the inner scale, in m. The default is 1 um.
        seed:         an integer seed, or None for a random screen.
        subharmonics: True adds the three subharmonic levels.

    Returns:
        An n x n array of the phase, in radians.
    """
    aotools = _load_aotools()
    ps = aotools.turbulence.phasescreen

    # aotools calculates f0 = 1/L0 and then (f^2 + f0^2)^(11/6). An infinite
    # L0 gives f0 = 0.0, so the f = 0 pixel divides by zero. aotools sets
    # that pixel to zero after the divide, so the screen is correct, but
    # numpy prints a RuntimeWarning. A very large finite L0 gives the same
    # screen and no warning.
    L0 = 1.0e6 if not np.isfinite(L0_m) else float(L0_m)

    gen = ps.ft_sh_phase_screen if subharmonics else ps.ft_phase_screen
    return np.asarray(gen(float(r0_m), int(n), float(pixel_m), L0,
                          float(l0_m), seed=seed), dtype=float)


def _phase_psd_unit(f, L0_m, l0_m):
    """Give the modified von Karman phase PSD at r0 = 1 m.

        PHI_1(f) = 0.023 exp(-(f/fm)^2) / (f^2 + f0^2)^(11/6)

    with fm = 5.92/(2*pi*l0) and f0 = 1/L0. It is the r0-free part of Schmidt
    (2010), DOI 10.1117/3.866274, Ch. 9, Eq. (9.51), printed p. 161, in the
    ORDINARY-frequency form of Eq. (9.52), printed p. 161. The full spectrum
    scales as r0^(-5/3), so a screen scales the sqrt of this by r0^(-5/6). The
    f = 0 value is infinite for an infinite L0; the caller sets that sample to
    zero, as Listing 9.2, line 16, printed p. 167, does.

    Args:
        f:    the radial ordinary frequency, in 1/m. A scalar or an ndarray.
        L0_m: the outer scale, in m. Infinite for no outer scale.
        l0_m: the inner scale, in m. Zero or less for no inner scale.

    Returns:
        The phase PSD at r0 = 1 m, in rad^2 m^2.
    """
    f = np.asarray(f, dtype=float)
    fm = np.inf if l0_m <= 0.0 else 5.92 / (2.0 * np.pi * float(l0_m))
    f0 = 0.0 if not np.isfinite(L0_m) else 1.0 / float(L0_m)
    denom = (f * f + f0 * f0) ** (11.0 / 6.0)
    top = 0.023 * np.exp(-(f / fm) ** 2)
    return np.where(denom > 0.0, top / np.where(denom > 0.0, denom, 1.0),
                    np.inf)


class ScreenFactory:
    """A fast, cached generator of von Karman phase screens for one grid.

    THE POINT. It gives the SAME physics as `phase_screen`, but it is fast. The
    baseline `phase_screen` calls aotools once for each screen. aotools makes
    the whole spectrum every time, and it builds the subharmonic sum from 27
    full-grid complex exponentials. The profiling of P0 (see
    `validation/waveoptics_speed/`) shows the subharmonic sum is about 70% of a
    trial. This factory removes that cost with three structural wins:

    1. THE FILTER IS CACHED. The square root of the phase PSD depends on
       (n, pixel, L0, l0) only. The Fried parameter r0 enters as the scalar
       factor r0^(-5/6). So the factory builds the filter ONE time and it
       scales the filter for each screen.
    2. THE SUBHARMONIC BASIS IS SEPARABLE. Each subharmonic mode
       exp(i 2 pi (fx x + fy y)) is the outer product of two length-n vectors.
       The factory builds the three per-level basis matrices ONE time, then it
       makes the whole 9-mode sum of each level with two small matrix products,
       not 9 full-grid exponentials. See Schmidt, DOI 10.1117/3.866274,
       Eq. (9.81), printed p. 169, and Lane, Glindemann and Dainty,
       DOI 10.1088/0959-7174/2/3/003.
    3. TWO SCREENS FOR ONE TRANSFORM. One complex Gaussian draw and one inverse
       transform give TWO independent real screens, the real part and the
       imaginary part. Ch. 9, text below Listing 9.2, printed p. 167, states
       that the two parts are uncorrelated. `make_stack` uses this pairing.
       `make` uses the real part only.

    THE DRAWS DIFFER FROM aotools. This generator and aotools do NOT give the
    same screen for the same seed. Each one is a correct, independent draw of
    the same random field. The wave-optics runner keeps aotools as the default,
    so an old run stays bit-identical.

    THE SPECTRUM. It is the modified von Karman phase PSD of Schmidt,
    DOI 10.1117/3.866274, Ch. 9, Eq. (9.51), printed p. 161, in the ordinary
    frequency of Eq. (9.52). The Fourier-series screen is Eqs. (9.78) to (9.80),
    printed pp. 166 and 167. The subharmonic screen is Eq. (9.81), printed
    p. 169. The default scales give the pure Kolmogorov spectrum.

    Attributes:
        n:            the number of pixels along one side.
        pixel_m:      the pixel pitch, in m.
        L0_m:         the outer scale, in m.
        l0_m:         the inner scale, in m.
        subharmonics: True if the factory adds the subharmonic levels.
    """

    _EXPONENT = -5.0 / 6.0                 # r0 enters the screen as r0^(-5/6).

    def __init__(self, n, pixel_m, L0_m=np.inf, l0_m=1e-6, subharmonics=True,
                 n_sub_levels=3, dtype=np.float64):
        """Build the cached filter and the subharmonic basis for one grid.

        Args:
            n:            the number of pixels along one side.
            pixel_m:      the pixel pitch of the propagation grid, in m.
            L0_m:         the outer scale, in m. The default is infinite.
            l0_m:         the inner scale, in m. The default is 1 um. It matches
                          the default of `phase_screen`.
            subharmonics: True adds the subharmonic levels. Keep it True: the
                          low-frequency content drives the beam wander.
            n_sub_levels: the number of subharmonic levels. The book uses 3.
            dtype:        the OUTPUT floating type, numpy.float64 (the default)
                          or numpy.float32. float32 halves the memory and it
                          measures a small error; see the module self-check.
        """
        self.n = int(n)
        self.pixel_m = float(pixel_m)
        self.L0_m = float(L0_m)
        self.l0_m = float(l0_m)
        self.subharmonics = bool(subharmonics)
        self.n_sub_levels = int(n_sub_levels)
        self._rdtype = np.float32 if dtype == np.float32 else np.float64
        self._cdtype = (np.complex64 if dtype == np.float32
                        else np.complex128)

        n = self.n
        dx = self.pixel_m

        # ---- win 1: the cached high-frequency filter, at r0 = 1 m ----
        # Schmidt, Ch. 9, Eqs. (9.78) to (9.80), printed pp. 166 and 167:
        # cn = (g1 + i g2) sqrt(PHI) df, then phi = Re{ift(cn)}. The df factor
        # and the r0 = 1 spectrum fold into one cached array.
        df = 1.0 / (n * dx)
        fx = (np.arange(n) - n // 2) * df
        FX, FY = np.meshgrid(fx, fx)
        psd = _phase_psd_unit(np.hypot(FX, FY), self.L0_m, self.l0_m)
        psd[n // 2, n // 2] = 0.0             # Listing 9.2, line 16, p. 167.
        self._filt = (np.sqrt(psd) * df).astype(self._rdtype)

        # ---- win 2: the separable subharmonic basis, at r0 = 1 m ----
        # Schmidt, Ch. 9, Eq. (9.81), printed p. 169. For each level p the
        # frequency set is a 3 by 3 grid at the pitch df_p = 1/(3^p L). The
        # basis exp(i 2 pi f_p m x) is the same on the two axes, so one n by 3
        # matrix E_p serves both. The level sum is E_p @ C @ E_p.T.
        self._sub_filt = []
        self._E = []
        if self.subharmonics:
            side = n * dx
            x = (np.arange(n) - n // 2) * dx
            for p in range(1, self.n_sub_levels + 1):
                df_p = 1.0 / (3.0 ** p * side)     # Eq. (9.81), p. 169.
                f3 = np.array([-1.0, 0.0, 1.0]) * df_p
                FX3, FY3 = np.meshgrid(f3, f3)
                psd3 = _phase_psd_unit(np.hypot(FX3, FY3), self.L0_m, self.l0_m)
                psd3[1, 1] = 0.0            # Listing 9.3, line 26, p. 170.
                self._sub_filt.append((np.sqrt(psd3) * df_p).astype(
                    self._rdtype))
                E = np.exp(1j * 2.0 * np.pi * np.outer(x, f3))   # n by 3.
                self._E.append(E.astype(self._cdtype))

    def _ift_series(self, cn):
        """Give the bare Fourier-series sum of the coefficient grid cn.

        It is `ift2(cn, 1.0)` of Schmidt, DOI 10.1117/3.866274, Ch. 2,
        Eq. (2.9), printed p. 17, with df = 1: the centred inverse transform
        times n^2. The result is Eq. (9.78), printed p. 167.
        """
        n = self.n
        shifted = np.fft.fftshift(np.fft.ifft2(np.fft.ifftshift(cn)))
        return (shifted * (n * n)).astype(self._cdtype)

    def _subharmonic(self, r0_m, rng):
        """Give the low-frequency subharmonic screen for one r0 and one rng.

        It is Schmidt, DOI 10.1117/3.866274, Ch. 9, Eq. (9.81), printed p. 169,
        with the mean removed as Listing 9.3, line 38, printed p. 170, does.
        The 9-mode sum of each level is the matrix product E @ C @ E.T of the
        separable basis.
        """
        if not self.subharmonics:
            return np.zeros((self.n, self.n), dtype=self._rdtype)
        scale = float(r0_m) ** self._EXPONENT
        lo = np.zeros((self.n, self.n), dtype=self._cdtype)
        for E, sub_filt in zip(self._E, self._sub_filt):
            g = (rng.standard_normal((3, 3)) + 1j * rng.standard_normal((3, 3)))
            c = (g.astype(self._cdtype) * sub_filt * scale)
            lo += E @ c @ E.T
        out = np.real(lo)
        return (out - out.mean()).astype(self._rdtype)

    def _base_pair(self, rng):
        """Draw one coefficient grid and return the paired real and imaginary
        screens at r0 = 1 m.

        The real part and the imaginary part of the inverse transform are two
        independent, correctly scaled screens (Ch. 9, text below Listing 9.2,
        printed p. 167). Scale each by r0^(-5/6) for the wanted r0.
        """
        n = self.n
        g = (rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n)))
        cn = g.astype(self._cdtype) * self._filt
        full = self._ift_series(cn)
        return np.real(full), np.imag(full)

    def make(self, r0_m, rng):
        """Make one phase screen, in radians.

        It uses win 1 (the cached filter) and win 2 (the separable
        subharmonics). It draws one coefficient grid and it keeps the real part.

        Args:
            r0_m: the Fried parameter of the slab, in m.
            rng:  a numpy.random.Generator. Build one per screen from the
                  runner seed, for example
                  numpy.random.default_rng(_screen_seed(entropy, k, j)).

        Returns:
            An n x n array of the phase, in radians.
        """
        base, _ = self._base_pair(rng)
        hi = (float(r0_m) ** self._EXPONENT) * base
        return (hi + self._subharmonic(r0_m, rng)).astype(self._rdtype)

    def make_stack(self, r0_m_array, rng):
        """Make a stack of phase screens, one for each r0, in radians.

        It uses win 3 (the pairing): one coefficient grid and one inverse
        transform give TWO independent base screens, so a stack of m screens
        needs only ceil(m/2) transforms. The two screens of a pair take
        different r0 values, because the r0^(-5/6) scale is applied after the
        transform.

        THE SEED. This call draws from ONE Generator for the whole stack. That
        is a different seed granularity from `make`, which takes one Generator
        per screen. The runner uses `make` to keep the per-screen seed contract;
        this call is the fast batch form for a study.

        Args:
            r0_m_array: the Fried parameter of each screen, in m.
            rng:        a single numpy.random.Generator for the whole stack.

        Returns:
            A list of n x n arrays of the phase, in radians.
        """
        r0 = np.asarray(r0_m_array, dtype=float).ravel()
        m = r0.size
        out = [None] * m
        k = 0
        while k < m:
            re, im = self._base_pair(rng)
            out[k] = ((r0[k] ** self._EXPONENT) * re
                      + self._subharmonic(r0[k], rng)).astype(self._rdtype)
            if k + 1 < m:
                out[k + 1] = ((r0[k + 1] ** self._EXPONENT) * im
                              + self._subharmonic(r0[k + 1], rng)).astype(
                                  self._rdtype)
            k += 2
        return out


def Screen(Fin, phase_rad):
    """Put a phase screen into the field.

        E_out(x,y) = E_in(x,y) * exp(i * phi(x,y))

    The screen is a thin, pure phase element. The power does not change.
    It is the refraction operator T = exp(-i psi) of Schmidt (2010),
    DOI 10.1117/3.866274, Ch. 9, Eq. (9.2), printed p. 150. The sign of the
    exponent follows the phase convention of the field.

    Args:
        Fin:       the input field. It must be on a flat grid.
        phase_rad: an N x N array of the phase, in radians.

    Returns:
        A new Field.

    Raises:
        ValueError: the field is in spherical coordinates, or the phase
            array has the wrong shape.
    """
    if Fin._curvature != 0.0:
        raise ValueError('Screen: the field is in spherical coordinates. '
                         'Use Convert() first. A co-moving screen is not '
                         'implemented.')
    phase_rad = np.asarray(phase_rad, dtype=float)
    if phase_rad.shape != (Fin.N, Fin.N):
        raise ValueError(f'Screen: the phase array is {phase_rad.shape}, '
                         f'but the field is ({Fin.N}, {Fin.N})')
    Fout = Field.copy(Fin)
    Fout.field = Fout.field * np.exp(1j * phase_rad)
    Fout._IsGauss = False
    return Fout


if __name__ == '__main__':
    import time

    from ..field import Begin, Power
    from ..sources import GaussBeam

    lam = 1550e-9

    # ---- 0. the checks that need no aotools ----
    # screen_r0 follows the (-3/5) power law, and it accepts an array.
    k = 2.0 * np.pi / lam
    cn2i = 1e-13
    r0_one = screen_r0(cn2i, lam)
    assert abs(r0_one - (0.423 * k * k * cn2i) ** (-0.6)) < 1e-15, r0_one
    # Two times the integral gives 2^(-3/5) times r0.
    assert abs(screen_r0(2 * cn2i, lam) / r0_one - 2 ** -0.6) < 1e-12
    r0_arr = screen_r0(np.array([cn2i, 2 * cn2i]), lam)
    assert r0_arr.shape == (2,), r0_arr.shape
    assert abs(r0_arr[0] - r0_one) < 1e-15, (r0_arr[0], r0_one)

    # Screen() on a hand-made phase array.
    n0 = 64
    F0 = GaussBeam(Begin(20e-3, lam, n0), 3e-3)
    tilt = 0.7 * np.linspace(-1, 1, n0)[None, :] * np.ones((n0, 1))
    FS = Screen(F0, tilt)
    p_in, p_out = Power(F0), Power(FS)
    assert abs(p_out / p_in - 1.0) < 1e-14, (p_in, p_out)
    assert FS._IsGauss is False
    assert np.allclose(np.angle(FS.field) - np.angle(F0.field), tilt, atol=1e-9)

    # The documented failure modes, asserted AS failures.
    try:
        Screen(F0, np.zeros((n0 + 1, n0 + 1)))
        raise AssertionError('Screen must refuse a wrong-shape phase')
    except ValueError as exc:
        assert 'phase array' in str(exc), str(exc)
    Fsph = Field.copy(F0)
    Fsph._curvature = -1e-6
    try:
        Screen(Fsph, np.zeros((n0, n0)))
        raise AssertionError('Screen must refuse a spherical field')
    except ValueError as exc:
        assert 'Convert' in str(exc), str(exc)

    def d_phi(scr, kpx):
        """Measure D_phi at the separation kpx pixels, on the two axes."""
        dh = scr[:, kpx:] - scr[:, :-kpx]
        dv = scr[kpx:, :] - scr[:-kpx, :]
        return 0.5 * (np.mean(dh * dh) + np.mean(dv * dv))

    # ---- F. the fast ScreenFactory (no aotools needed) ----
    # The factory gives the SAME physics as phase_screen, with a different
    # random draw. It must pass the structure-function bounds of cases 2, 3 and
    # 5 below. These checks run whether or not aotools is present.

    # Fa. the seed is deterministic, and a different seed gives a new screen.
    fac = ScreenFactory(128, 0.01)
    sa = fac.make(0.1, np.random.default_rng(7))
    sb = fac.make(0.1, np.random.default_rng(7))
    sc = fac.make(0.1, np.random.default_rng(8))
    assert np.array_equal(sa, sb), 'the same rng must give the same screen'
    assert not np.allclose(sa, sc), 'a different rng must give a new screen'
    assert sa.shape == (128, 128), sa.shape

    # Fb. the structure function (case 2 for the factory).
    # D_phi(r) = 6.88 (r/r0)^(5/3). Fried, DOI 10.1364/JOSA.56.001372.
    nF, dxF, r0F, MF = 512, 0.01, 0.10, 40
    ksF = np.array([3, 5, 8, 11, 16])          # r/r0 = 0.3 to 1.6.
    facF = ScreenFactory(nF, dxF)
    accF = np.zeros(ksF.size)
    for i in range(MF):
        s = facF.make(r0F, np.random.default_rng(5000 + i))
        for j, kpx in enumerate(ksF):
            accF[j] += d_phi(s, kpx)
    accF /= MF
    theoryF = 6.88 * ((ksF * dxF) / r0F) ** (5.0 / 3.0)
    ratio_fac = accF / theoryF
    assert np.all(ratio_fac < 1.02), ratio_fac
    assert np.all(ratio_fac > 0.85), ratio_fac

    # Fc. the subharmonics are necessary (case 3 for the factory).
    nF2, dxF2, r0F2, MF2 = 256, 0.01, 0.10, 30
    kbigF = 32                                 # r = 3.2 r0.
    fac_on = ScreenFactory(nF2, dxF2, subharmonics=True)
    fac_off = ScreenFactory(nF2, dxF2, subharmonics=False)
    d_on_f = d_off_f = 0.0
    for i in range(MF2):
        d_on_f += d_phi(fac_on.make(r0F2, np.random.default_rng(6000 + i)),
                        kbigF)
        d_off_f += d_phi(fac_off.make(r0F2, np.random.default_rng(6000 + i)),
                         kbigF)
    d_on_f /= MF2
    d_off_f /= MF2
    th_bigF = 6.88 * ((kbigF * dxF2) / r0F2) ** (5.0 / 3.0)
    assert d_off_f < 0.8 * d_on_f, (d_off_f, d_on_f)
    assert d_off_f < 0.6 * th_bigF, (d_off_f, th_bigF)

    # Fd. two screens add as r0^(-5/3) (case 5 for the factory).
    def fit_r0_f(scr_list, dx_fit, kpx):
        d = np.mean([d_phi(s, kpx) for s in scr_list])
        return (kpx * dx_fit) * (6.88 / d) ** 0.6

    r0Fa, r0Fb, kfitF, MF3 = 0.10, 0.16, 8, 30
    sFa, sFb, sFsum = [], [], []
    for i in range(MF3):
        a1 = fac_on.make(r0Fa, np.random.default_rng(7000 + i))
        b1 = fac_on.make(r0Fb, np.random.default_rng(8000 + i))
        sFa.append(a1)
        sFb.append(b1)
        sFsum.append(a1 + b1)
    faF = fit_r0_f(sFa, dxF2, kfitF)
    fbF = fit_r0_f(sFb, dxF2, kfitF)
    fabF = fit_r0_f(sFsum, dxF2, kfitF)
    assert abs(fabF / (faF ** (-5 / 3) + fbF ** (-5 / 3)) ** -0.6 - 1.0) \
        < 0.05, (fabF, faF, fbF)

    # Fe. make_stack (win 3, the pairing) gives the same statistics as make.
    # A stack of many screens at one r0 must give the same mean D_phi as make
    # at the same separation. The band-limit deficit is the same for both, so
    # the two agree even where each is biased against the theory.
    stack = fac_on.make_stack([r0F2] * 60, np.random.default_rng(9000))
    d_stack = np.mean([d_phi(s, kfitF) for s in stack])
    d_make = np.mean([d_phi(fac_on.make(r0F2, np.random.default_rng(9500 + i)),
                            kfitF) for i in range(60)])
    assert abs(d_stack / d_make - 1.0) < 0.10, (d_stack, d_make)

    # Ff. the float32 switch measures a small error against float64.
    fac64 = ScreenFactory(256, 0.01, dtype=np.float64)
    fac32 = ScreenFactory(256, 0.01, dtype=np.float32)
    s64 = fac64.make(0.1, np.random.default_rng(11))
    s32 = fac32.make(0.1, np.random.default_rng(11))
    assert s32.dtype == np.float32, s32.dtype
    rms32 = float(np.sqrt(np.mean((s32 - s64) ** 2)))
    rel32 = rms32 / float(np.std(s64))
    assert rel32 < 1e-3, rel32

    try:
        _load_aotools()
    except ImportError as exc:
        print('aotools is absent, so only the ScreenFactory is checked.')
        print(f'  {exc}')
        print("")
        print(f"ScreenFactory structure function, n={nF}, "
              f"r0 = {r0F / dxF:.0f} px, {MF} screens:")
        for kpx, dm, th, rt in zip(ksF, accF, theoryF, ratio_fac):
            print(f"  r/r0 {kpx * dxF / r0F:5.2f}  D meas {dm:8.3f}  "
                  f"D theory {th:8.3f}  ratio {rt:6.3f}")
        print(f"float32 vs float64 rel rms {rel32:.2e}")
        print('self-check passed')
        raise SystemExit(0)

    # ---- 1. the seed is deterministic ----
    a = phase_screen(0.1, 128, 0.01, seed=7)
    b = phase_screen(0.1, 128, 0.01, seed=7)
    c = phase_screen(0.1, 128, 0.01, seed=8)
    assert np.array_equal(a, b), 'the same seed must give the same screen'
    assert not np.allclose(a, c), 'a different seed must give a new screen'
    assert a.shape == (128, 128), a.shape

    # ---- 2. the phase structure function ----
    # D_phi(r) = 6.88 (r/r0)^(5/3). Fried, DOI 10.1364/JOSA.56.001372.
    t0 = time.time()
    n, dx, r0 = 512, 0.01, 0.10            # r0 = 10 pixels, side = 51 r0
    M = 40
    ks = np.array([3, 5, 8, 11, 16])       # r from 3 dx to side/32
    acc = np.zeros(ks.size)
    for i in range(M):
        s = phase_screen(r0, n, dx, seed=1000 + i)
        for j, kpx in enumerate(ks):
            acc[j] += d_phi(s, kpx)
    acc /= M
    theory = 6.88 * ((ks * dx) / r0) ** (5.0 / 3.0)
    ratio_sh = acc / theory
    t_sf = time.time() - t0
    # The screen holds no power above the grid Nyquist frequency and too
    # little power below 1/(n*dx). So D_phi stays BELOW the theory. The
    # three subharmonic levels keep the deficit inside 15% over this range.
    assert np.all(ratio_sh < 1.02), ratio_sh
    assert np.all(ratio_sh > 0.85), ratio_sh

    # ---- 3. the subharmonics are necessary (a failure mode) ----
    n2, dx2, r0_2, M2 = 256, 0.01, 0.10, 30
    kbig = 32                              # r = 3.2 r0
    d_on = d_off = 0.0
    for i in range(M2):
        d_on += d_phi(phase_screen(r0_2, n2, dx2, seed=2000 + i), kbig)
        d_off += d_phi(phase_screen(r0_2, n2, dx2, seed=2000 + i,
                                    subharmonics=False), kbig)
    d_on /= M2
    d_off /= M2
    th_big = 6.88 * ((kbig * dx2) / r0_2) ** (5.0 / 3.0)
    assert d_off < 0.8 * d_on, (d_off, d_on)
    assert d_off < 0.6 * th_big, (d_off, th_big)

    # ---- 4. screen_r0 agrees with aotools cn2_to_r0 ----
    import aotools as _aot
    r0_ref = _aot.cn2_to_r0(cn2i, lam)     # aotools defaults to 500 nm
    assert abs(screen_r0(cn2i, lam) / r0_ref - 1.0) < 1e-12, (r0_one, r0_ref)

    # ---- 5. two screens add as r0^(-5/3) ----
    def fit_r0(scr_list, dx_fit, kpx):
        """Read r0 back from D_phi at one separation."""
        d = np.mean([d_phi(s, kpx) for s in scr_list])
        return (kpx * dx_fit) * (6.88 / d) ** 0.6

    r0_a, r0_b = 0.10, 0.16
    r0_ab = (r0_a ** (-5.0 / 3.0) + r0_b ** (-5.0 / 3.0)) ** (-3.0 / 5.0)
    kfit, M3 = 8, 30
    s_a, s_b, s_sum = [], [], []
    for i in range(M3):
        sa = phase_screen(r0_a, n2, dx2, seed=3000 + i)
        sb = phase_screen(r0_b, n2, dx2, seed=4000 + i)
        s_a.append(sa)
        s_b.append(sb)
        s_sum.append(sa + sb)
    fa = fit_r0(s_a, dx2, kfit)
    fb = fit_r0(s_b, dx2, kfit)
    fab = fit_r0(s_sum, dx2, kfit)
    # The band-limit deficit lifts each read-back r0 by the same factor, so
    # it cancels between the parts and the sum. The tolerance stays loose.
    assert abs(fab / (fa ** (-5 / 3) + fb ** (-5 / 3)) ** -0.6 - 1.0) < 0.05, \
        (fab, fa, fb)
    # The sum carries the same deficit factor as one screen. It does NOT
    # match r0_ab, because the read-back is biased. See case 2.
    assert abs(fab / r0_ab - fa / r0_a) < 0.05, (fab / r0_ab, fa / r0_a)

    # ---- 6. Screen() on a real screen keeps the power ----
    nn = 256
    Fg = GaussBeam(Begin(0.2, lam, nn), 0.02)
    scr = phase_screen(0.05, nn, 0.2 / nn, seed=99)
    Fout = Screen(Fg, scr)
    assert abs(Power(Fout) / Power(Fg) - 1.0) < 1e-14

    print(f"wavelength                 {lam * 1e9:9.1f} nm")
    print(f"r0 for Cn2 int 1e-13       {r0_one * 1e2:9.3f} cm")
    print(f"aotools cn2_to_r0          {r0_ref * 1e2:9.3f} cm")
    print("")
    print(f"structure function, n={n}, r0 = {r0 / dx:.0f} px, {M} screens:")
    for kpx, dm, th, rt in zip(ks, acc, theory, ratio_sh):
        print(f"  r/r0 {kpx * dx / r0:5.2f}  D meas {dm:8.3f}  "
              f"D theory {th:8.3f}  ratio {rt:6.3f}")
    print(f"  (elapsed {t_sf:.1f} s)")
    print("")
    print(f"subharmonics at r/r0 {kbig * dx2 / r0_2:.1f}, n={n2}:")
    print(f"  on     {d_on:9.3f}   theory {th_big:9.3f}")
    print(f"  off    {d_off:9.3f}   deficit {1 - d_off / d_on:6.3f}")
    print("")
    print(f"two screens, r0 read back from D_phi at r/r0 "
          f"{kfit * dx2 / r0_a:.1f}, n={n2}:")
    print(f"  screen A       set {r0_a * 1e2:6.2f} cm   read {fa * 1e2:6.2f}"
          f" cm   bias {fa / r0_a:5.3f}")
    print(f"  screen B       set {r0_b * 1e2:6.2f} cm   read {fb * 1e2:6.2f}"
          f" cm   bias {fb / r0_b:5.3f}")
    print(f"  A + B          set {r0_ab * 1e2:6.2f} cm   read {fab * 1e2:6.2f}"
          f" cm   bias {fab / r0_ab:5.3f}")
    print("")
    print(f"ScreenFactory structure function, n={nF}, r0 = {r0F / dxF:.0f} px, "
          f"{MF} screens:")
    for kpx, dm, th, rt in zip(ksF, accF, theoryF, ratio_fac):
        print(f"  r/r0 {kpx * dxF / r0F:5.2f}  D meas {dm:8.3f}  "
              f"D theory {th:8.3f}  ratio {rt:6.3f}")
    print(f"ScreenFactory subharmonics at r/r0 {kbigF * dxF2 / r0F2:.1f}: "
          f"on {d_on_f:.3f}, off {d_off_f:.3f}, deficit "
          f"{1 - d_off_f / d_on_f:.3f}")
    print(f"ScreenFactory float32 vs float64 rel rms {rel32:.2e}")
    print("self-check passed")
