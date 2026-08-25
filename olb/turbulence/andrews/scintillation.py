'''
Scintillation index of Andrews and Phillips, from weak to strong fluctuation.

This module gives the normalised irradiance variance (the scintillation index)
of a plane wave, a spherical wave, and a Gaussian beam. It covers the weak
regime (Ch. 8) and the strong regime (Ch. 9, the extended Rytov theory). It also
gives the two log-irradiance variances that feed the gamma-gamma distribution.

Source of every equation:
    L. C. Andrews and R. L. Phillips, "Laser Beam Propagation through Random
    Media", 2nd ed., SPIE Press (2005). DOI: 10.1117/3.626196
Each function names its section, its equation number, and its printed page.

The three quantities:

- `rytov_variance` gives sigma_R^2 (plane), beta_0^2 = 0.4 sigma_R^2
  (spherical), and sigma_B^2 (Gaussian beam). Each one is the WEAK-fluctuation
  scintillation index of that wave type. The book uses sigma_R^2 as the
  strength-of-turbulence measure everywhere.
- `large_scale_log_variance` and `small_scale_log_variance` give sigma_lnX^2 and
  sigma_lnY^2 of the extended Rytov theory. They feed
  `olb.turbulence.andrews.distributions.gamma_gamma_params(sigma2_lnX,
  sigma2_lnY)` with no change.
- `scintillation_index` gives the index itself. Weak = Ch. 8.2. Strong =
  exp(sigma_lnX^2 + sigma_lnY^2) - 1, Ch. 9, Eq. (28), printed p. 333.

REGIME BOUNDARY. The book calls the fluctuations weak when sigma_R^2 < 1 (Ch. 8,
text below Eq. (23), printed pp. 264-265; Ch. 12, Eq. (40), printed p. 497). The
"auto" regime uses that boundary. For a Gaussian beam the book adds a second
condition, sigma_R^2 Lambda^(5/6) < 1 (Ch. 5, Eq. (16), printed p. 140, quoted
again on printed p. 265). The caller must test that second condition. This
module does not gate on it.

PLANE OF REFERENCE. This module takes ONE path length L and ONE scalar Cn2, so
it makes no path integral and it picks no reference plane. The book path
variable is xi = 1 - z/L (Ch. 8, text at Eq. (4), printed p. 261), which is
measured from the RECEIVER. A caller that integrates a Cn2 profile must choose
the reference plane itself.

This module holds physics only. It returns no decibels.

SCOPE (deferred): the finite inner-scale and outer-scale branches (Ch. 9,
Secs. 9.4.2, 9.5.2 and 9.6.3) are NOT built. The `l0` and `L0` keywords exist,
and a value other than None raises NotImplementedError.
'''

import numpy as np

from .beam import BeamParams, beam_params, effective_beam_params, wavenumber

# Weak-fluctuation boundary on the plane-wave Rytov variance. Source: Andrews
# and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196, Ch. 8, text below Eq. (23),
# printed pp. 264-265. Ch. 12, Eq. (40), printed p. 497, repeats it.
WEAK_REGIME_LIMIT = 1.0

_WAVES = ('plane', 'spherical', 'gaussian')


def _reject_scales(l0, L0):
    '''Refuse a finite inner scale or outer scale. The branch is not built.'''
    if l0 is not None or L0 is not None:
        raise NotImplementedError(
            'l0/L0 branch - WP3. Andrews and Phillips, 2nd ed. (2005), '
            'DOI 10.1117/3.626196, Secs. 9.4.2 (printed p. 337), 9.5.2 '
            '(printed p. 343) and 9.6.3 (printed p. 354) give the two-scale '
            'forms. This module builds the zero-inner-scale, infinite-outer-'
            'scale branch only.')


def _check_wave(wave):
    '''Refuse a wave type that this module does not know.'''
    if wave not in _WAVES:
        raise ValueError(f'wave must be one of {_WAVES}, not {wave!r}')


def _need_beam(beam, wave):
    '''Return the beam parameters, or refuse if the caller gave none.'''
    if wave == 'gaussian' and beam is None:
        raise ValueError('wave="gaussian" needs beam=BeamParams(...)')
    return beam


def rytov_variance(wavelength, z, cn2, *, wave='plane', beam=None):
    '''
    Return the weak-fluctuation Rytov variance of the named wave type.

    Parameters:
        wavelength : float or numpy.ndarray
            Optical wavelength [m].
        z : float or numpy.ndarray
            Path length L [m].
        cn2 : float or numpy.ndarray
            Refractive-index structure constant [m^-2/3]. It is constant over
            the path.
        wave : str
            "plane", "spherical", or "gaussian".
        beam : BeamParams, optional
            The beam parameters at the receiver. Required for "gaussian".

    Returns:
        float or numpy.ndarray
            sigma_R^2 (plane), beta_0^2 (spherical), or sigma_B^2 (Gaussian).

    formula:
        plane      sigma_R^2 = 1.23 Cn2 k^(7/6) L^(11/6),   k = 2*pi/lambda
        spherical  beta_0^2  = 0.40 sigma_R^2
        gaussian   sigma_B^2 = 3.86 sigma_R^2
                     { 0.40 [(1+2 Theta)^2 + 4 Lambda^2]^(5/12)
                       cos[ (5/6) arctan( (1+2 Theta) / (2 Lambda) ) ]
                       - (11/16) Lambda^(5/6) }
    Source: Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196:
        plane and spherical   Ch. 8, Eq. (20), printed p. 264. Ch. 9, Eqs. (63)
                              and (64), printed p. 341, restate them.
        gaussian              Ch. 8, Eq. (23), printed p. 264, longitudinal
                              half. Ch. 9, Eq. (93), printed p. 350, restates
                              it. The Ch. 8 Summary prints it again as Eq. (130),
                              printed p. 303.

    RESTRICTION on "gaussian": the book states that Eq. (23) holds "in the case
    of a collimated or divergent beam" (Ch. 8, text above Eq. (23), printed
    p. 264). So Theta0 must be 1 or more. A convergent beam needs the exact
    hypergeometric form, Ch. 8, Eq. (19), printed p. 263, which this module does
    not build. A convergent beam raises NotImplementedError.
    '''
    _check_wave(wave)
    k = wavenumber(wavelength)
    sigma2_R = (1.23 * np.asarray(cn2, dtype=float) * k ** (7.0 / 6.0)
                * np.asarray(z, dtype=float) ** (11.0 / 6.0))
    if wave == 'plane':
        return sigma2_R
    if wave == 'spherical':
        return 0.40 * sigma2_R
    return beam_rytov_variance(sigma2_R, _need_beam(beam, wave))


def beam_rytov_variance(sigma2_R, beam):
    '''
    Return the Gaussian-beam Rytov variance sigma_B^2 from sigma_R^2.

    This is the longitudinal (on-axis) scintillation index of a Gaussian beam
    under weak fluctuations.

    Parameters:
        sigma2_R : float or numpy.ndarray
            The PLANE-wave Rytov variance over the same path.
        beam : BeamParams
            The beam parameters at the receiver.

    Returns:
        float or numpy.ndarray
            sigma_B^2.

    See `rytov_variance` for the formula, the citations, and the restriction to
    a collimated or a divergent beam.
    '''
    if np.any(np.asarray(beam.theta0, dtype=float) < 1.0 - 1e-12):
        raise NotImplementedError(
            'convergent beam (Theta0 < 1). Andrews and Phillips, 2nd ed. '
            '(2005), DOI 10.1117/3.626196, state at printed p. 264 that '
            'Ch. 8, Eq. (23) holds for a collimated or a divergent beam only. '
            'Use the exact Ch. 8, Eq. (19), printed p. 263.')
    sigma2_R = np.asarray(sigma2_R, dtype=float)
    a = 1.0 + 2.0 * beam.theta
    modulus = (a ** 2 + 4.0 * beam.lam ** 2) ** (5.0 / 12.0)
    phase = np.cos((5.0 / 6.0) * np.arctan2(a, 2.0 * beam.lam))
    return 3.86 * sigma2_R * (0.40 * modulus * phase
                              - (11.0 / 16.0) * beam.lam ** (5.0 / 6.0))


def large_scale_log_variance(sigma2_R, *, wave='plane', l0=None, L0=None,
                             beam=None, r=0.0):
    '''
    Return the large-scale log-irradiance variance sigma_lnX^2.

    This is the first of the two variances that the extended Rytov theory needs.
    Feed it, with `small_scale_log_variance`, straight into
    `gamma_gamma_params(sigma2_lnX, sigma2_lnY)`.

    Parameters:
        sigma2_R : float or numpy.ndarray
            The PLANE-wave Rytov variance. Every branch below scales from the
            plane-wave value, which is what the book does.
        wave : str
            "plane", "spherical", or "gaussian".
        l0, L0 : None
            Inner scale and outer scale [m]. Not built. A value raises
            NotImplementedError.
        beam : BeamParams, optional
            The beam parameters at the receiver. Required for "gaussian".
        r : float
            Off-axis radius [m]. It has NO effect. The book splits only the
            LONGITUDINAL component into a large-scale and a small-scale part.
            The radial component stays separate (Ch. 9, Eq. (103), printed
            p. 353). The keyword is here to match `scintillation_index`.

    Returns:
        float or numpy.ndarray
            sigma_lnX^2.

    formula:
        plane      0.49 s^2 / (1 + 1.11 s^(12/5))^(7/6),   s^2 = sigma_R^2
        spherical  0.20 s^2 / (1 + 0.19 s^(12/5))^(7/6)
        gaussian   0.49 b^2 / (1 + 0.56 (1 + Theta) b^(12/5))^(7/6),
                   b^2 = sigma_B^2
    Source: Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196:
        plane      Ch. 9, Eq. (41), printed p. 335
        spherical  Ch. 9, Eq. (69), printed p. 342
        gaussian   Ch. 9, Eq. (97), printed p. 352
    '''
    _check_wave(wave)
    _reject_scales(l0, L0)
    s2 = np.asarray(sigma2_R, dtype=float)
    if wave == 'plane':
        return 0.49 * s2 / (1.0 + 1.11 * s2 ** (6.0 / 5.0)) ** (7.0 / 6.0)
    if wave == 'spherical':
        return 0.20 * s2 / (1.0 + 0.19 * s2 ** (6.0 / 5.0)) ** (7.0 / 6.0)
    bm = _need_beam(beam, wave)
    b2 = beam_rytov_variance(s2, bm)
    denom = (1.0 + 0.56 * (1.0 + bm.theta) * b2 ** (6.0 / 5.0)) ** (7.0 / 6.0)
    return 0.49 * b2 / denom


def small_scale_log_variance(sigma2_R, *, wave='plane', l0=None, L0=None,
                             beam=None, r=0.0):
    '''
    Return the small-scale log-irradiance variance sigma_lnY^2.

    This is the second of the two variances that the extended Rytov theory
    needs. See `large_scale_log_variance` for the parameters. The keyword `r`
    has no effect for the same reason.

    formula:
        plane      0.51 s^2 / (1 + 0.69 s^(12/5))^(5/6),   s^2 = sigma_R^2
        spherical  0.20 s^2 / (1 + 0.23 s^(12/5))^(5/6)
        gaussian   0.51 b^2 / (1 + 0.69 b^(12/5))^(5/6),   b^2 = sigma_B^2
    Source: Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196:
        plane      Ch. 9, Eq. (46), printed p. 336
        spherical  Ch. 9, Eq. (72), printed p. 342
        gaussian   Ch. 9, Eq. (101), printed p. 352
    In the saturation regime each branch goes to ln 2, which gives the small-
    scale variance its limit sigma_Y^2 -> 1 (Ch. 9, text at Eq. (35), printed
    p. 334).
    '''
    _check_wave(wave)
    _reject_scales(l0, L0)
    s2 = np.asarray(sigma2_R, dtype=float)
    if wave == 'plane':
        return 0.51 * s2 / (1.0 + 0.69 * s2 ** (6.0 / 5.0)) ** (5.0 / 6.0)
    if wave == 'spherical':
        return 0.20 * s2 / (1.0 + 0.23 * s2 ** (6.0 / 5.0)) ** (5.0 / 6.0)
    b2 = beam_rytov_variance(s2, _need_beam(beam, wave))
    return 0.51 * b2 / (1.0 + 0.69 * b2 ** (6.0 / 5.0)) ** (5.0 / 6.0)


def _radial_component(sigma2_R, lam_eff, w_eff, r, tracked, pointing_error_m,
                      wander_rms_m):
    '''
    Return the radial (off-axis) component of the Gaussian-beam index.

    formula:
        tracked    4.42 s^2 Lambda^(5/6) (r - sqrt(<rc^2>))^2 / W^2,
                   for r > sqrt(<rc^2>), else 0
        untracked  4.42 s^2 Lambda^(5/6) [ (r + sigma_pe)^2 U(r - sigma_pe)
                                           + sigma_pe^2 ] / W^2
    The extra sigma_pe^2 of the untracked form is the wander-induced pointing
    error. The book puts it in the LONGITUDINAL component, but the total is the
    same. Source: Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196,
    Ch. 8, Eqs. (40), (41) and (44), printed pp. 274-276, and Ch. 9, Eqs. (88),
    (90), (91) and (103), printed pp. 350 and 353. In the weak regime use the
    free-space Lambda and W. In the strong regime use Lambda_e and W_LT
    (Ch. 9, Eq. (88), printed p. 350).
    '''
    r = np.asarray(r, dtype=float)
    coef = 4.42 * sigma2_R * lam_eff ** (5.0 / 6.0) / w_eff ** 2
    if tracked:
        rc = np.asarray(wander_rms_m, dtype=float)
        return coef * np.where(r > rc, (r - rc) ** 2, 0.0)
    pe = np.asarray(pointing_error_m, dtype=float)
    return coef * (np.where(r > pe, (r + pe) ** 2, 0.0) + pe ** 2)


def scintillation_index(wavelength, z, cn2, *, wave='plane', regime='auto',
                        l0=None, L0=None, beam=None, r=0.0, tracked=True,
                        pointing_error_m=0.0, wander_rms_m=0.0):
    '''
    Return the scintillation index sigma_I^2 for a single path.

    Parameters:
        wavelength : float or numpy.ndarray
            Optical wavelength [m].
        z : float or numpy.ndarray
            Path length L [m].
        cn2 : float or numpy.ndarray
            Refractive-index structure constant [m^-2/3], constant over the
            path.
        wave : str
            "plane", "spherical", or "gaussian".
        regime : str
            "weak", "strong", or "auto". "auto" uses the book boundary
            sigma_R^2 < 1 (see WEAK_REGIME_LIMIT).
        l0, L0 : None
            Inner scale and outer scale [m]. Not built. A value raises
            NotImplementedError.
        beam : BeamParams, optional
            The beam parameters at the receiver. Required for "gaussian".
        r : float or numpy.ndarray
            Off-axis radius in the receiver plane [m]. It acts on a Gaussian
            beam only. A plane wave and a spherical wave have no radial
            component (Ch. 8, text below Eq. (17), printed p. 263).
        tracked : bool
            True removes the beam wander, per the Ch. 8.3.2 tracked model. False
            keeps it, per the Ch. 8.3.1 untracked model.
        pointing_error_m : float
            The rms wander-induced pointing error sigma_pe [m]. It acts on the
            UNTRACKED model only. Andrews gives it in Ch. 8, Eq. (36), printed
            p. 273. That equation needs a beam-wander module, which this package
            does not hold yet, so the caller must supply the number. The default
            0 gives the plain first-order Rytov result.
        wander_rms_m : float
            The rms beam-wander displacement sqrt(<rc^2>) [m]. It acts on the
            TRACKED model only. Andrews gives it in Ch. 8, Eq. (33), printed
            p. 272. Same note as above. The default 0 gives the plain
            first-order Rytov result.

    Returns:
        float or numpy.ndarray
            sigma_I^2.

    formula:
        weak    plane      sigma_R^2                     Ch. 8, Eq. (20), p. 264
                spherical  0.40 sigma_R^2                Ch. 8, Eq. (20), p. 264
                gaussian   radial + sigma_B^2            Ch. 8, Eq. (23), p. 264
        strong  exp( sigma_lnX^2 + sigma_lnY^2 ) - 1     Ch. 9, Eq. (28), p. 333
                plus the radial component for a Gaussian beam, with Lambda_e
                and W_LT                                 Ch. 9, Eq. (103), p. 353
    The strong form covers 0 <= sigma_R^2 < infinity (Ch. 9, Eqs. (47), (73) and
    (102), printed pp. 336, 342 and 352). It reduces to the weak form as
    sigma_R -> 0. So "auto" only picks the SIMPLER expression below the
    boundary; it does not switch physics.
    '''
    _check_wave(wave)
    _reject_scales(l0, L0)
    bm = _need_beam(beam, wave)

    sigma2_R = rytov_variance(wavelength, z, cn2, wave='plane')

    def weak():
        if wave == 'plane':
            return sigma2_R
        if wave == 'spherical':
            return 0.40 * sigma2_R
        radial = _radial_component(sigma2_R, bm.lam, bm.w, r, tracked,
                                   pointing_error_m, wander_rms_m)
        return radial + beam_rytov_variance(sigma2_R, bm)

    def strong():
        x = large_scale_log_variance(sigma2_R, wave=wave, beam=bm)
        y = small_scale_log_variance(sigma2_R, wave=wave, beam=bm)
        out = np.exp(x + y) - 1.0
        if wave != 'gaussian':
            return out
        bm_e = effective_beam_params(bm, sigma2_R)
        return out + _radial_component(sigma2_R, bm_e.lam, bm_e.w, r, tracked,
                                       pointing_error_m, wander_rms_m)

    if regime == 'weak':
        return weak()
    if regime == 'strong':
        return strong()
    if regime != 'auto':
        raise ValueError(f'regime must be "weak", "strong" or "auto", '
                         f'not {regime!r}')
    if np.ndim(sigma2_R) == 0:
        return weak() if sigma2_R < WEAK_REGIME_LIMIT else strong()
    return np.where(sigma2_R < WEAK_REGIME_LIMIT, weak(), strong())


if __name__ == '__main__':
    # ---------------- physics self-checks ----------------
    lam_m = 1550e-9
    L = 2000.0
    cn2_weak = 3e-16

    s2_R = rytov_variance(lam_m, L, cn2_weak)
    assert s2_R < 0.1, s2_R
    assert np.isclose(rytov_variance(lam_m, L, cn2_weak, wave='spherical'),
                      0.4 * s2_R)

    # A very wide beam is a plane wave. A point source is a spherical wave.
    bp_plane = beam_params(50.0, lam_m, L)
    bp_sph = beam_params(1e-5, lam_m, L)
    b2_plane = rytov_variance(lam_m, L, cn2_weak, wave='gaussian',
                              beam=bp_plane)
    b2_sph = rytov_variance(lam_m, L, cn2_weak, wave='gaussian', beam=bp_sph)
    # The rounded book constants 3.86 and 0.40 give 0.998 and 0.3996, not the
    # exact 1 and 0.4. So the tolerance is 3 parts in 1000.
    assert abs(b2_plane / s2_R - 1.0) < 3e-3, b2_plane / s2_R
    assert abs(b2_sph / s2_R - 0.40) < 3e-3, b2_sph / s2_R

    # A convergent beam is refused, not guessed.
    try:
        rytov_variance(lam_m, L, cn2_weak, wave='gaussian',
                       beam=beam_params(0.05, lam_m, L, 4000.0))
    except NotImplementedError:
        pass
    else:
        raise AssertionError('a convergent beam must raise')

    # A finite inner scale or outer scale is refused, not guessed.
    for kwargs in ({'l0': 5e-3}, {'L0': 1.0}):
        try:
            large_scale_log_variance(s2_R, wave='plane', **kwargs)
        except NotImplementedError:
            pass
        else:
            raise AssertionError('l0/L0 must raise')

    # The strong index saturates near 1 and never runs away.
    s2_strong = scintillation_index(lam_m, L, 1e-12, wave='plane')
    assert 0.5 < s2_strong < 3.0, s2_strong
    # The index peaks in the focusing regime, and the spherical-wave peak comes
    # LATER than the plane-wave peak. Ch. 9, text below Eq. (73), printed p. 343,
    # puts the two peaks near sigma_R = 2 and sigma_R = 4 (Fig. 9.7).
    grid = np.logspace(-1.0, 1.5, 4000)
    peak_s = grid[np.argmax(np.exp(
        large_scale_log_variance(grid ** 2, wave='plane')
        + small_scale_log_variance(grid ** 2, wave='plane')) - 1.0)]
    peak_sp = grid[np.argmax(np.exp(
        large_scale_log_variance(grid ** 2, wave='spherical')
        + small_scale_log_variance(grid ** 2, wave='spherical')) - 1.0)]
    assert 2.0 < peak_s < 4.0, peak_s
    assert 3.0 < peak_sp < 6.0, peak_sp
    assert peak_sp > peak_s, (peak_sp, peak_s)

    # The radial component is zero on axis and grows off axis.
    bp = beam_params(0.05, lam_m, L)
    on_axis = scintillation_index(lam_m, L, cn2_weak, wave='gaussian', beam=bp)
    off_axis = scintillation_index(lam_m, L, cn2_weak, wave='gaussian', beam=bp,
                                   r=0.5 * bp.w)
    assert off_axis > on_axis, (off_axis, on_axis)
    # An untracked beam with a pointing error scintillates more than a tracked
    # one (Ch. 8, Fig. 8.8, printed p. 276).
    untracked = scintillation_index(lam_m, L, cn2_weak, wave='gaussian',
                                    beam=bp, tracked=False,
                                    pointing_error_m=0.2 * bp.w)
    assert untracked > on_axis, (untracked, on_axis)

    # "auto" picks weak below the boundary and strong above it.
    assert scintillation_index(lam_m, L, cn2_weak) == s2_R
    assert scintillation_index(lam_m, L, 1e-12) == s2_strong

    # ---------------- REDUCTION checks ----------------
    from .. import beam_wave_scintillation as bws
    from .. import plane_wave_scintillation as pws

    # 3. rytov_variance(plane) reproduces plane_wave_scintillation.sigma1_rytov.
    ref = pws.sigma1_rytov(cn2_weak, lam_m, L) ** 2
    err = abs(s2_R - ref) / ref
    assert err < 1e-12, err
    print(f'REDUCTION rytov_variance(plane) : rel err = {err:.3e}  '
          f'(target 1e-12)')

    # 4. The strong plane form reproduces the fixed closed form. The parent now
    # DELEGATES to this module, so this check confirms the wiring. The second
    # comparison is independent: the d -> 0 limit of the Andrews aperture-
    # averaged form (Ch. 10, Eq. (69), printed p. 413) carries its own copy of
    # the four constants 0.49, 1.11, 0.51 and 0.69.
    cn2_mid = 1e-14
    mine = scintillation_index(lam_m, L, cn2_mid, wave='plane', regime='strong')
    parent = pws.plane_wave_scintillation_index_closed(cn2_mid, lam_m, L)
    err_wire = abs(mine - parent)
    assert err_wire < 1e-9, err_wire
    d0 = pws.aperture_averaged_index_andrews(0.0, cn2_mid, lam_m, L)
    err_d0 = abs(mine - d0)
    assert err_d0 < 1e-9, err_d0
    print(f'REDUCTION strong plane closed form : parent err = {err_wire:.3e}  '
          f'independent d=0 err = {err_d0:.3e}  (target 1e-9)')

    # 6. The strong model reduces to sigma_R^2 in the weak limit.
    target = 0.01
    cn2_small = cn2_weak * target / s2_R
    weak_limit = scintillation_index(lam_m, L, cn2_small, wave='plane',
                                     regime='strong')
    pct = abs(weak_limit - target) / target * 100.0
    assert pct < 3.0, pct
    print(f'REDUCTION strong -> sigma_R^2 at sigma_R^2 = 0.01 : '
          f'{pct:.3f} % (target 3 %)')

    # 5. GAP 9. The weak Gaussian on-axis index of Ch. 8, Eq. (23) against the
    # Dios path integral of beam_wave_scintillation. Homogeneous Cn2, one
    # horizontal path, collimated, sigma_R^2 < 0.1.
    hs = np.linspace(0.0, L, 800)
    cn2_flat = np.full_like(hs, cn2_weak)
    w0 = 0.05

    mine_coll = scintillation_index(lam_m, L, cn2_weak, wave='gaussian',
                                    beam=beam_params(w0, lam_m, L),
                                    regime='weak')
    dios_coll = bws.on_axis_scintillation_index(hs, cn2_flat, w0, lam_m)
    gap9_coll = (mine_coll - dios_coll) / dios_coll * 100.0
    assert abs(gap9_coll) < 15.0, gap9_coll
    print(f'GAP 9 collimated w0={w0} m, sigma_R^2={s2_R:.4f} : '
          f'Andrews Eq. (23) = {mine_coll:.6f}  Dios = {dios_coll:.6f}  '
          f'diff = {gap9_coll:+.2f} %')

    # The same measurement for a divergent beam, which is the uplink_flux use.
    f0_div = -1000.0
    mine_div = scintillation_index(lam_m, L, cn2_weak, wave='gaussian',
                                   beam=beam_params(w0, lam_m, L, f0_div),
                                   regime='weak')
    dios_div = bws.on_axis_scintillation_index(hs, cn2_flat, w0, lam_m,
                                               f0=f0_div)
    gap9_div = (mine_div - dios_div) / dios_div * 100.0
    print(f'GAP 9 divergent f0={f0_div} m : '
          f'Andrews Eq. (23) = {mine_div:.6f}  Dios = {dios_div:.6f}  '
          f'diff = {gap9_div:+.2f} %  (no assert)')

    print('self-check passed')
