"""The split-step engine: propagate, apply a screen, propagate again.

The module moves a field along a flat-grid path and puts the given phase
screens at the given distances. It owns NO random number generator. The
caller makes the screens (see screens.py) and gives them as radian arrays.
So one seed gives one repeatable path.

The module gives two functions:

- super_gaussian_boundary: an absorbing mask for the edge of the grid.
- split_step:              the propagate-screen-propagate loop.

Sources:
- Schmidt, Numerical Simulation of Optical Wave Propagation with Examples in
  MATLAB, DOI 10.1117/3.866274. The split-step chain is Ch. 9, Eqs. (9.1) to
  (9.3), printed p. 150. The general partial-propagation chain is Ch. 8,
  Eq. (8.18), printed p. 139. The absorbing boundary is Ch. 8, Eq. (8.1),
  printed p. 134. The step cap is Ch. 8, Eq. (8.24), printed p. 144, repeated
  as Ch. 9, Eq. (9.89), printed p. 174. NOTE: this module keeps ONE flat pitch
  for the whole path, and Eq. (8.18) gives each step its own pitch, from the
  linear rule of Eq. (8.8), printed p. 136. See docs/schmidt-crosscheck.md,
  gap S-14.
- Andrews and Phillips, Laser Beam Propagation through Random Media, 2nd ed.,
  DOI 10.1117/3.626196, Ch. 8. The plane-wave Rytov variance that the
  self-check reproduces.
"""

import numpy as np

from ..field import Field
from ..propagators import Forvard
from .screens import Screen


def super_gaussian_boundary(n, width_frac=0.125, power=8):
    """Make an absorbing mask for the edge of the grid.

    The mask is 1.0 inside the radius (1 - width_frac) of the half-side. It
    falls as a super-Gaussian in the band outside that radius:

        m(rho) = exp(-t^p),   t = (rho - r_flat) / width_frac,  t >= 0

    Here rho is the radius in units of the half-side, so rho is 1.0 at the
    middle of an edge and 1.41 at a corner. The mask is exp(-1) at the
    middle of an edge and it is zero at the corners.

    The mask removes the energy that reaches the edge of the grid. The FFT
    propagator is periodic, so that energy comes back at the opposite edge.
    The concept is Schmidt (2010), DOI 10.1117/3.866274, Ch. 8, Eq. (8.1),
    printed p. 134: g = exp(-(r/sigma)^n) with n > 2.

    THE NUMBERS ARE NOT THE BOOK'S, AND THE SHAPE IS NOT EITHER. Eq. (8.1)
    gives the family and no values. The book RUNS power 16 and one half-width
    sigma = 0.47 N pixels from the centre (Listing 8.1, line 11, printed
    p. 142; Listing 9.7, line 19, printed p. 179; Fig. 8.1, printed p. 134,
    prints 0.45 L). This function runs power 8 and a taper BAND of 0.125 of
    the half-side, with a hard flat region inside it. So the book mask is
    0.068 at the middle of an edge and this mask is exp(-1) = 0.368 there: the
    book absorbs about 5 times harder. The book itself records that Flatte and
    others used n = 8 (Ch. 8, text, printed p. 134), so the power has a source
    in the literature, not in the book's own runs. The conflict is RECORDED,
    not changed. See docs/schmidt-crosscheck.md, gap S-15.

    Args:
        n:          the number of pixels along one side.
        width_frac: the width of the absorbing band, as a fraction of the
                    half-side.
        power:      the exponent p of the super-Gaussian.

    Returns:
        An n x n array of real numbers in the interval [0, 1].

    Raises:
        ValueError: width_frac is not inside (0, 1], or power is not
            positive.
    """
    if not 0.0 < width_frac <= 1.0:
        raise ValueError('super_gaussian_boundary: width_frac must be '
                         'inside (0, 1]')
    if power <= 0:
        raise ValueError('super_gaussian_boundary: power must be positive')
    x = (np.arange(n) - n // 2) / (n / 2.0)         # -1.0 .. 1.0
    rho = np.hypot(x[:, None], x[None, :])
    r_flat = 1.0 - width_frac
    t = np.clip((rho - r_flat) / width_frac, 0.0, None)
    return np.exp(-(t ** power))


def _substeps(gap_m, max_step_m):
    """Cut one gap into equal sub-steps that obey the sampling limit.

    Args:
        gap_m:      the distance to cover, in m.
        max_step_m: the longest sub-step that the grid samples well, in m.

    Returns:
        A 1d array of equal sub-steps. It is empty when the gap is zero.
    """
    if gap_m <= 0.0:
        return np.zeros(0)
    count = int(np.ceil(gap_m / max_step_m))
    return np.full(count, gap_m / count)


def _apply_mask(Fin, mask):
    """Multiply the field by a real mask. Return a new Field."""
    Fout = Field.copy(Fin)
    Fout.field = Fout.field * mask
    Fout._IsGauss = False
    return Fout


def split_step(Fin, z_screens_m, screens, z_total_m, *, boundary=None,
               max_step_m=None):
    """Propagate a field along a path through a list of phase screens.

    The function makes one hop to each screen, applies that screen, and
    makes a last hop to z_total_m. Each hop uses Forvard, the FFT spectral
    propagator.

    A hop that is longer than max_step_m breaks into equal sub-steps. The
    default limit is

        max_step = N * dx^2 / lambda

    Past that range the quadratic phase of the transfer function turns
    faster than one sample, so the propagator aliases. The rule is Schmidt
    (2010), DOI 10.1117/3.866274, Ch. 8, Eq. (8.24), printed p. 144, and the
    turbulent statement of it is Ch. 9, Eq. (9.89), printed p. 174. Both come
    from constraint 4, Ch. 7, Eq. (7.59), printed p. 127. (The same formula is
    in olb.waveoptics.grid.forvard_max_z. This module does not import that
    module, because grid.py reads the rest of olb.)

    THE MASK IS NECESSARY. The sub-steps alone remove NO aliasing: the
    sampled transfer function of the full step is the product of the
    sampled transfer functions of the sub-steps, so a split hop gives the
    same array as one long hop. The book gives the same result: Ch. 8,
    Eqs. (8.19) to (8.22), printed pp. 139 and 143, show that the
    intermediate pitches cancel and that constraint 3 does not change with
    the number of partial propagations. The mask is the part that helps. It
    removes the energy at the edge of the grid between two sub-steps, before
    the periodic propagator brings that energy back at the opposite edge. The
    book puts the boundary operator A at each intermediate plane too, Ch. 8,
    Eq. (8.18), printed p. 139. Give a boundary from
    super_gaussian_boundary() for any path that spreads the beam.

    THE SCREEN PLACEMENT DIFFERS FROM THE BOOK. Ch. 9, Eq. (9.3), printed
    p. 150, puts one screen AT each partial-propagation plane. The olb planner
    puts each screen at the Cn2-weighted centre of a merged slab, so the two
    differ by half a slab. The book does not treat that placement.

    Args:
        Fin:         the input field. It must be on a flat grid.
        z_screens_m: the distance of each screen from the input plane, in m.
                     The distances go up, and they stay inside
                     [0, z_total_m].
        screens:     one N x N phase array, in radians, for each distance.
        z_total_m:   the distance from the input plane to the output plane,
                     in m.
        boundary:    an N x N mask, or None. The function applies it after
                     each sub-step and after each screen.
        max_step_m:  the longest single hop, in m. None takes the default.

    Returns:
        A new Field at z_total_m.

    Raises:
        ValueError: the field is in spherical coordinates, the distances
            are not sorted or not inside the path, the number of screens
            does not match, or an array has the wrong shape.
    """
    if Fin._curvature != 0.0:
        raise ValueError('split_step: the field is in spherical coordinates. '
                         'Use Convert() first. A co-moving split step is not '
                         'implemented.')
    z = np.asarray(z_screens_m, dtype=float).ravel()
    screens = list(screens)
    if len(screens) != z.size:
        raise ValueError(f'split_step: {len(screens)} screens for {z.size} '
                         'distances')
    if z_total_m < 0.0:
        raise ValueError('split_step: z_total_m must not be negative')
    if z.size:
        if np.any(np.diff(z) < 0.0):
            raise ValueError('split_step: z_screens_m must go up')
        if z[0] < 0.0 or z[-1] > z_total_m:
            raise ValueError('split_step: z_screens_m must stay inside '
                             f'[0, {z_total_m}]')
    shape = (Fin.N, Fin.N)
    for i, scr in enumerate(screens):
        if np.shape(scr) != shape:
            raise ValueError(f'split_step: screen {i} is {np.shape(scr)}, '
                             f'but the field is {shape}')
    if boundary is not None:
        boundary = np.asarray(boundary, dtype=float)
        if boundary.shape != shape:
            raise ValueError(f'split_step: the boundary is {boundary.shape}, '
                             f'but the field is {shape}')
    if max_step_m is None:
        max_step_m = Fin.N * Fin.dx ** 2 / Fin.lam
    if max_step_m <= 0.0:
        raise ValueError('split_step: max_step_m must be positive')

    def hop(F, gap_m):
        """Cover one gap, and apply the mask after each sub-step."""
        for dz in _substeps(gap_m, max_step_m):
            F = Forvard(F, dz)
            if boundary is not None:
                F = _apply_mask(F, boundary)
        return F

    Fout = Field.copy(Fin)
    here = 0.0
    for zi, scr in zip(z, screens):
        Fout = hop(Fout, zi - here)
        Fout = Screen(Fout, scr)
        if boundary is not None:
            Fout = _apply_mask(Fout, boundary)
        here = zi
    return hop(Fout, z_total_m - here)


if __name__ == '__main__':
    import time

    from ..field import Begin, Power
    from ..sources import GaussBeam

    lam = 1550e-9

    # ---- 1. no screens and no mask is a plain Forvard ----
    n, side, z = 256, 0.2, 100.0
    F0 = GaussBeam(Begin(side, lam, n), 0.01)
    max_z = n * F0.dx ** 2 / lam
    assert z < max_z, (z, max_z)            # one sub-step, so it is exact
    Fa = split_step(F0, [], [], z)
    Fb = Forvard(F0, z)
    assert np.allclose(Fa.field, Fb.field), 'the empty path must be Forvard'

    # ---- 2. the boundary mask ----
    mask = super_gaussian_boundary(n)
    assert mask.shape == (n, n)
    assert mask.max() <= 1.0 and mask.min() >= 0.0
    # The mask is exactly 1.0 over the central 70% of the half-side.
    xx = (np.arange(n) - n // 2) / (n / 2.0)
    rr = np.hypot(xx[:, None], xx[None, :])
    assert np.all(np.abs(mask[rr <= 0.70] - 1.0) < 1e-12), mask[rr <= 0.70].min()
    assert mask[0, 0] < 1e-12, mask[0, 0]   # the corner is absorbed

    # A well contained beam does not see the mask.
    F_no = split_step(F0, [], [], z)
    F_bd = split_step(F0, [], [], z, boundary=mask)
    cut = slice(n // 4, 3 * n // 4)
    assert np.allclose(F_no.field[cut, cut], F_bd.field[cut, cut],
                       atol=1e-9 * np.abs(F_no.field).max()), 'interior moved'

    # On a tight grid the beam reaches the edge. The mask removes it.
    w0_t, z_t = 2e-3, 400.0
    F_tight = GaussBeam(Begin(0.02, lam, n), w0_t)
    m_t = super_gaussian_boundary(n)
    edge = rr > 0.9
    e_no = np.abs(split_step(F_tight, [], [], z_t).field[edge]) ** 2
    e_bd = np.abs(split_step(F_tight, [], [], z_t,
                             boundary=m_t).field[edge]) ** 2
    assert e_bd.sum() < e_no.sum(), (e_bd.sum(), e_no.sum())

    # ---- 3. a long gap breaks into sub-steps ----
    assert _substeps(0.0, 10.0).size == 0
    assert _substeps(10.0, 10.0).size == 1
    s3 = _substeps(25.0, 10.0)
    assert s3.size == 3, s3
    assert abs(s3.sum() - 25.0) < 1e-12, s3
    assert np.all(s3 <= 10.0 + 1e-12), s3
    assert np.allclose(s3, s3[0]), s3        # the sub-steps are equal
    # A split path and one long hop agree, because Forvard composes.
    z_long = 3.5 * max_z
    Fs = split_step(F0, [], [], z_long)
    Fl = Forvard(F0, z_long)
    assert _substeps(z_long, max_z).size == 4
    assert np.allclose(Fs.field, Fl.field, atol=1e-9 * np.abs(Fl.field).max())

    # ---- 4. the guards (documented failure modes) ----
    Fsph = Field.copy(F0)
    Fsph._curvature = -1e-6
    bad = [
        ('spherical', lambda: split_step(Fsph, [], [], z)),
        ('count', lambda: split_step(F0, [10.0], [], z)),
        ('order', lambda: split_step(F0, [20.0, 10.0],
                                     [np.zeros((n, n))] * 2, z)),
        ('range', lambda: split_step(F0, [2 * z], [np.zeros((n, n))], z)),
        ('shape', lambda: split_step(F0, [10.0], [np.zeros((4, 4))], z)),
        ('mask', lambda: split_step(F0, [], [], z,
                                    boundary=np.ones((4, 4)))),
    ]
    for name, call in bad:
        try:
            call()
            raise AssertionError(f'split_step must refuse the {name} case')
        except ValueError:
            pass

    # ---- 5. the power bookkeeping ----
    scr5 = [0.3 * np.cos(np.linspace(0, 9, n))[None, :] * np.ones((n, 1))
            for _ in range(4)]
    z5 = [20.0, 40.0, 60.0, 80.0]
    p_in = Power(F0)
    p_no = Power(split_step(F0, z5, scr5, z))
    p_bd = Power(split_step(F0, z5, scr5, z, boundary=mask))
    assert abs(p_no / p_in - 1.0) < 1e-10, (p_no, p_in)
    assert p_bd <= p_no * (1.0 + 1e-12), (p_bd, p_no)

    # ---- 6. the plane-wave Rytov variance ----
    # sigma2_R = 1.23 Cn2 k^(7/6) L^(11/6). Andrews and Phillips,
    # DOI 10.1117/3.626196, Ch. 8 (the weak plane-wave index).
    sig2_meas = None
    try:
        from .screens import phase_screen, screen_r0
        phase_screen(0.1, 8, 0.01, seed=0)
    except ImportError as exc:
        print('aotools is absent, so the Rytov case is not checked.')
        print(f'  {exc}')
    else:
        t0 = time.time()
        L, cn2 = 1000.0, 5e-15
        kw = 2.0 * np.pi / lam
        sig2_R = 1.23 * cn2 * kw ** (7.0 / 6.0) * L ** (11.0 / 6.0)
        n6, side6, n_scr, trials = 256, 0.5, 8, 60
        dz = L / n_scr
        z_scr = (np.arange(n_scr) + 0.5) * dz       # the slab mid-points
        r0_scr = screen_r0(cn2 * dz, lam)
        F_pw = Begin(side6, lam, n6)                # a unit plane wave
        dx6 = F_pw.dx
        # The screens carry NO subharmonics here. Forvard is periodic, and
        # this case fills the whole grid, so it uses no absorbing mask. The
        # subharmonics are not periodic on the grid, so they break the wrap
        # and they raise the measured variance. The scintillation comes
        # from the Fresnel scale sqrt(lambda*L), which the high-frequency
        # screen samples well.
        cut6 = slice(n6 // 4, 3 * n6 // 4)
        tot = tot2 = 0.0
        count = 0
        for i in range(trials):
            scr = [phase_screen(r0_scr, n6, dx6, seed=10000 + 100 * i + j,
                                subharmonics=False) for j in range(n_scr)]
            Fout = split_step(F_pw, z_scr, scr, L)
            I = np.abs(Fout.field[cut6, cut6]) ** 2
            tot += I.sum()
            tot2 += (I * I).sum()
            count += I.size
        mean_I = tot / count
        sig2_meas = tot2 / count / mean_I ** 2 - 1.0
        t_ry = time.time() - t0
        assert abs(sig2_meas / sig2_R - 1.0) < 0.20, (sig2_meas, sig2_R)

    print(f"wavelength                {lam * 1e9:9.1f} nm")
    print(f"grid                      {n:9d} px, {side * 1e3:.1f} mm")
    print(f"max hop N dx^2 / lambda   {max_z:9.1f} m")
    print(f"sub-steps for {z_long:.0f} m      {_substeps(z_long, max_z).size:9d}")
    print("")
    print("boundary mask, tight grid, energy outside 0.9 of the half-side:")
    print(f"  no mask                 {e_no.sum():9.3e}")
    print(f"  with mask               {e_bd.sum():9.3e}")
    print("")
    print("power through 4 phase screens:")
    print(f"  input                   {p_in:9.6e} W")
    print(f"  no mask                 {p_no:9.6e} W")
    print(f"  with mask               {p_bd:9.6e} W")
    if sig2_meas is not None:
        print("")
        print(f"plane-wave scintillation, {n_scr} screens, {trials} trials:")
        print(f"  path length             {L:9.1f} m")
        print(f"  Cn2                     {cn2:9.2e} m^-2/3")
        print(f"  r0 per screen           {r0_scr * 1e2:9.2f} cm")
        print(f"  Fresnel scale           {np.sqrt(lam * L) * 1e3:9.2f} mm")
        print(f"  sigma2_I measured       {sig2_meas:9.4f}")
        print(f"  sigma2_R theory         {sig2_R:9.4f}")
        print(f"  ratio                   {sig2_meas / sig2_R:9.3f}")
        print(f"  (elapsed {t_ry:.1f} s)")
    print("self-check passed")
