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

    try:
        _load_aotools()
    except ImportError as exc:
        print('aotools is absent, so the screen statistics are not checked.')
        print(f'  {exc}')
        print('self-check passed')
        raise SystemExit(0)

    def d_phi(scr, kpx):
        """Measure D_phi at the separation kpx pixels, on the two axes."""
        dh = scr[:, kpx:] - scr[:, :-kpx]
        dv = scr[kpx:, :] - scr[:-kpx, :]
        return 0.5 * (np.mean(dh * dh) + np.mean(dv * dv))

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
    print("self-check passed")
