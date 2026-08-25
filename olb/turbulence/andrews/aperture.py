'''
Aperture averaging of irradiance fluctuations (Andrews).

A receive lens that is wider than the irradiance correlation width sees several
correlation patches at one time. The measured power then fluctuates less than the
irradiance at one point. This module gives that reduction.

- `averaged_index` gives the irradiance FLUX variance sigma_I^2(D_G) in the plane
  of the photodetector.
- `averaging_factor` gives A = sigma_I^2(D_G) / sigma_I^2(0).

Source of every equation:
    L. C. Andrews and R. L. Phillips, "Laser Beam Propagation through Random
    Media", 2nd ed., SPIE Press (2005). DOI: 10.1117/3.626196
Chapter 10, Sec. 10.3, printed pp. 409 to 421. Each function names its section,
its equation number, and its printed page.

SOFT APERTURE, NOT A HARD FILTER. The Andrews chain builds the flux variance
through the ABCD matrix of a thin lens with a GAUSSIAN limiting aperture of
radius W_G. The book then ties the "hard aperture" diameter D_G to that soft
radius by D_G^2 = 8 W_G^2 (Ch. 10, text below Eq. (57), printed p. 411). Every
integral in Sec. 10.3 carries the soft factor exp(-D_G^2 kappa^2/16), for example
Eq. (59), printed p. 412. So this module does NOT use the hard circular
modulation transfer function, and it does NOT use the Airy filter
[2 J1(x)/x]^2 that `olb/turbulence/plane_wave_scintillation.py` uses. The two
agree in the limits and differ in between. Andrews DOES print the Airy function
elsewhere, as the piston Zernike filter, Ch. 14, Eq. (86), printed p. 634, but
never for aperture averaging. See Conflict C-06 in docs/andrews-crosscheck.md.

NO ANNULAR APERTURE. A full-text search of all 809 pages finds no aperture
filter, no modulation transfer function and no flux variance for a centrally
obscured (annular) RECEIVE aperture. Secs. 10.3.1 to 10.3.6 use only the soft
Gaussian aperture, or the unobscured circular transfer function of Eq. (54),
printed p. 410. So olb gap 8 cannot be closed from this book. The Terms that
assume an unobscured aperture must keep saying so through `olb/assumptions.py`.
The negative result is recorded as row G-108 of docs/andrews-crosscheck.md.

GAP - the GAUSSIAN-BEAM TWO-SCALE chain is NOT built. Ch. 10, Eqs. (79) to (86),
printed pp. 419 and 420, need the parameter eta_X of Eq. (84), printed p. 420,
which is the same unresolved equation as Ch. 9, Eq. (109). See the docstring of
`olb.turbulence.andrews.scintillation`. The ZERO-scale Gaussian chain,
Eqs. (87) to (90), IS built.

This module holds physics only. It imports numpy, scipy and sibling andrews
modules. It returns no decibels.
'''

import numpy as np
from scipy.special import hyp2f1

from .beam import wavenumber
from .scintillation import (_eta_filter, beam_rytov_variance, rytov_variance,
                            two_scale_parameters, weak_two_scale_index)

_WAVES = ('plane', 'spherical', 'gaussian')
_REGIMES = ('weak', 'strong')


def d_param(D, wavelength, z):
    '''
    Return the nondimensional aperture parameter d = sqrt(k D_G^2 / (4 L)).

    Source: Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196, Ch. 10,
    Eq. (68), printed p. 413. The parameter is the aperture radius scaled by the
    Fresnel zone: d = (D_G/2) / sqrt(L/k).
    '''
    k = wavenumber(wavelength)
    return np.sqrt(k * np.asarray(D, dtype=float) ** 2
                   / (4.0 * np.asarray(z, dtype=float)))


def omega_g(D, wavelength, z):
    '''
    Return the lens spot parameter Omega_G = 2 L / (k W_G^2) = 16 L / (k D_G^2).

    Source: Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196, Ch. 10,
    text below Eq. (78), printed p. 419, with D_G^2 = 8 W_G^2 from the text below
    Eq. (57), printed p. 411. The book also writes Omega_G = 4/d^2 (Ch. 10, text
    at Fig. 10.13, printed p. 420).
    '''
    return 4.0 / d_param(D, wavelength, z) ** 2


def plane_weak_averaging_fit(D, wavelength, z):
    '''
    Return the book's own weak-fluctuation aperture-averaging fit for a plane
    wave.

    formula:
        A = [ 1 + 1.062 k D_G^2 / (4 L) ]^(-7/6) = (1 + 1.062 d^2)^(-7/6),
        sigma_R^2 < 1
    Source: Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196, Ch. 10,
    Eq. (61), printed p. 412. The book states the fit holds to about 7 % against
    the exact Eq. (60).

    This is NOT the Churnside fit that
    `olb.turbulence.plane_wave_scintillation.aperture_averaging_factor_weak`
    uses. That one is A = [1 + 1.07 (k D^2/(4 L))^(7/6)]^(-1), with the exponent
    INSIDE. Churnside, Applied Optics 30 (1991) 1982,
    DOI 10.1364/AO.30.001982. The two differ by more than 10 % for a large
    aperture. The module self-check prints the measured difference.
    '''
    return (1.0 + 1.062 * d_param(D, wavelength, z) ** 2) ** (-7.0 / 6.0)


def _weak_plane(d2, sigma2_R):
    '''
    Return the weak plane-wave flux variance, Kolmogorov spectrum.

    formula:
        u = k D_G^2/(16 L) = d^2/4
        sigma_I^2(D_G) = 3.86 sigma_R^2 { (1 + u^2)^(11/12)
                                          sin[(11/6) arctan(1/u)]
                                          - (11/6) u^(5/6) }
    Source: Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196, Ch. 10,
    Eq. (60), printed p. 412. At D_G = 0 it goes to 3.86 sin(11 pi/12)
    sigma_R^2 = 0.999 sigma_R^2.
    '''
    u = d2 / 4.0
    return 3.86 * sigma2_R * ((1.0 + u ** 2) ** (11.0 / 12.0)
                              * np.sin((11.0 / 6.0) * np.arctan2(1.0, u))
                              - (11.0 / 6.0) * u ** (5.0 / 6.0))


def _weak_spherical_factor(d2):
    '''
    Return the exact weak spherical-wave aperture-averaging factor.

    formula:
        u = k D_G^2/(16 L) = d^2/4
        A = 9.66 Re[ i^(5/6) 2F1(-5/6, 11/6; 17/6; 1 + i u) - (11/16) u^(5/6) ]
    Source: Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196, Ch. 10,
    Eq. (53), printed p. 409. The book scales Eq. (51) by
    sigma_I,sp^2 = 0.4 sigma_R^2, so A goes to 1 at D_G = 0. The rounded book
    constant 9.66 gives 1.0109 there, a 1.1 % offset.

    BRANCH CONVENTION. The book prints the argument as 1 - i k D_G^2/(16 L_1).
    With the scipy branch of 2F1 that argument, together with the printed factor
    i^(5/6), gives a factor that GROWS with the lens, which is not physical. The
    conjugate branch gives A(0) = 1 and a factor that falls to zero. Note that
    Re[i^(5/6) 2F1(...; 1 + i u)] equals Re[i^(-5/6) 2F1(...; 1 - i u)], because
    the three hypergeometric parameters are real. So this is a complex-branch
    convention, not a change of any coefficient. The signs of the coefficients
    9.66 and 11/16 are as printed.
    '''
    u = np.asarray(d2, dtype=float) / 4.0
    z_arg = 1.0 + 1j * u
    term = np.exp(1j * 5.0 * np.pi / 12.0) * hyp2f1(-5.0 / 6.0, 11.0 / 6.0,
                                                    17.0 / 6.0, z_arg)
    return 9.66 * (np.real(term) - (11.0 / 16.0) * u ** (5.0 / 6.0))


def _strong_plane(d2, s2):
    '''
    Return the all-regime plane-wave flux variance, zero inner scale.

    formula:
        sigma_I^2(D_G) = exp[ 0.49 s^2 / (1 + 0.65 d^2 + 1.11 s^(12/5))^(7/6)
                            + 0.51 s^2 (1 + 0.69 s^(12/5))^(-5/6)
                              / (1 + 0.90 d^2 + 0.62 d^2 s^(12/5)) ] - 1
    with s^2 = sigma_R^2. Source: Andrews and Phillips, 2nd ed. (2005),
    DOI 10.1117/3.626196, Ch. 10, Eq. (69), printed p. 413.
    '''
    s125 = s2 ** (6.0 / 5.0)
    large = 0.49 * s2 / (1.0 + 0.65 * d2 + 1.11 * s125) ** (7.0 / 6.0)
    small = (0.51 * s2 * (1.0 + 0.69 * s125) ** (-5.0 / 6.0)
             / (1.0 + 0.90 * d2 + 0.62 * d2 * s125))
    return np.exp(large + small) - 1.0


def _strong_spherical(d2, b0):
    '''
    Return the all-regime spherical-wave flux variance, zero inner scale.

    formula:
        sigma_I^2(D_G) = exp[ 0.49 b^2 / (1 + 0.18 d^2 + 0.56 b^(12/5))^(7/6)
                            + 0.51 b^2 (1 + 0.69 b^(12/5))^(-5/6)
                              / (1 + 0.90 d^2 + 0.62 d^2 b^(12/5)) ] - 1
    with b^2 = beta_0^2 = 0.4 sigma_R^2. Source: Andrews and Phillips, 2nd ed.
    (2005), DOI 10.1117/3.626196, Ch. 10, Eq. (77), printed p. 416.
    '''
    b125 = b0 ** (6.0 / 5.0)
    large = 0.49 * b0 / (1.0 + 0.18 * d2 + 0.56 * b125) ** (7.0 / 6.0)
    small = (0.51 * b0 * (1.0 + 0.69 * b125) ** (-5.0 / 6.0)
             / (1.0 + 0.90 * d2 + 0.62 * d2 * b125))
    return np.exp(large + small) - 1.0


def _strong_gaussian(D, wavelength, z, s2, bm):
    '''
    Return the all-regime Gaussian-beam flux variance, zero inner scale.

    formula:
        Omega_G = 16 L / (k D_G^2)
        sigma_lnX = 0.49 [(Omega_G - Lambda)/(Omega_G + Lambda)]^2 sigma_B^2
            / { 1 + 0.4 (2 - Theta) (sigma_B/sigma_R)^(12/7)
                    / [ (Omega_G + Lambda)
                        (1/3 - Theta/2 + Theta^2/5)^(6/7) ]
                  + 0.56 (1 + Theta) sigma_B^(12/5) }^(7/6)
        sigma_lnY = 0.51 sigma_B^2 (1 + 0.69 sigma_B^(12/5))^(-5/6)
            / { 1 + [ 1.20 (sigma_R/sigma_B)^(12/5) + 0.83 sigma_R^(12/5) ]
                    / (Omega_G - Lambda) }
        sigma_I^2(D_G) = exp(sigma_lnX + sigma_lnY) - 1,   Omega_G >= Lambda
    Source: Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196, Ch. 10,
    Eqs. (87), (88) and (89), printed p. 420. The Rytov variance sigma_B^2 is
    Eq. (90) on the same page, which this package holds as
    `scintillation.beam_rytov_variance`.

    THE CAPTURE EFFECT. Both faces carry (Omega_G - Lambda), so the flux
    variance goes to ZERO when the lens radius equals the incident beam radius.
    The book states that result below Eq. (78), printed p. 419. A plane-wave
    aperture-averaging model cannot reproduce it.

    RESTRICTION. The book prints Omega_G >= Lambda. A lens WIDER than the beam
    is outside the model, and this function refuses it.

    NOTE. Eq. (88) is an independent algebraic fit. It does NOT reduce exactly to
    the plane-wave Eq. (69) at Theta = 1, Lambda = 0, nor to the spherical-wave
    Eq. (77) at Theta = Lambda = 0. The module self-check measures both gaps.
    '''
    og = omega_g(D, wavelength, z)
    lm, th = bm.lam, bm.theta
    if np.any(og < lm):
        raise ValueError(
            'Omega_G < Lambda: the collecting lens is wider than the incident '
            'beam. Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196, '
            'print Ch. 10, Eqs. (87) to (89), printed p. 420, for '
            'Omega_G >= Lambda only.')
    sb2 = beam_rytov_variance(s2, bm)
    poly = (1.0 / 3.0 - th / 2.0 + th ** 2 / 5.0) ** (6.0 / 7.0)
    denom = (1.0 + 0.4 * (2.0 - th) * (sb2 / s2) ** (6.0 / 7.0)
             / ((og + lm) * poly)
             + 0.56 * (1.0 + th) * sb2 ** (6.0 / 5.0)) ** (7.0 / 6.0)
    large = 0.49 * ((og - lm) / (og + lm)) ** 2 * sb2 / denom
    small = (0.51 * sb2 * (1.0 + 0.69 * sb2 ** (6.0 / 5.0)) ** (-5.0 / 6.0)
             / (1.0 + (1.20 * (s2 / sb2) ** (6.0 / 5.0)
                       + 0.83 * s2 ** (6.0 / 5.0)) / (og - lm)))
    return np.exp(large + small) - 1.0


def _strong_two_scale(D, wavelength, z, cn2, wave, l0, L0, s2):
    '''
    Return the all-regime flux variance with a finite inner and outer scale.

    formula (plane, with Q_l and Q_0 from `two_scale_parameters`):
        eta_Xd  = 2.61 / (1 + 0.65 d^2 + 0.45 sigma_R^2 Q_l^(1/6))
        eta_Xd0 = eta_Xd Q_0 / (eta_Xd + Q_0)
        sigma_lnX = 0.16 sigma_R^2 [ F(eta_Xd) - F(eta_Xd0) ]
        sigma_lnY = 0.51 sigma_PL^2 (1 + 0.69 sigma_PL^(12/5))^(-5/6)
            / [ 1 + 0.90 d^2 (sigma_R/sigma_PL)^(12/5)
                  + 0.62 d^2 sigma_R^(12/5) ]
        sigma_I^2(D_G) = exp(sigma_lnX + sigma_lnY) - 1
    with F the shared filter group of `scintillation._eta_filter`. Source:
    Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196, Ch. 10,
    Eqs. (62) to (68), printed pp. 412 and 413.

    formula (spherical): the same shape with 0.04 beta_0^2, eta_X top 8.56,
    slope 0.20 beta_0^2 and the aperture term 0.18 d^2. Source: Ch. 10,
    Eqs. (71) to (76), printed p. 415.

    The weak two-scale index sigma_PL^2 (or sigma_SP^2) comes from
    `scintillation.weak_two_scale_index`.
    '''
    ql, q0 = two_scale_parameters(wavelength, z, l0, L0)
    d2 = d_param(D, wavelength, z) ** 2
    if wave == 'plane':
        coef, top, slope, ap = 0.16 * s2, 2.61, 0.45 * s2, 0.65
        point = s2
    else:
        b0 = 0.40 * s2
        coef, top, slope, ap = 0.04 * b0, 8.56, 0.20 * b0, 0.18
        point = b0
    eta_xd = top / (1.0 + ap * d2 + slope * ql ** (1.0 / 6.0))
    large = coef * _eta_filter(eta_xd, ql)
    if L0 is not None:
        eta_xd0 = eta_xd * q0 / (eta_xd + q0)
        large = large - coef * _eta_filter(eta_xd0, ql)
    weak = weak_two_scale_index(wavelength, z, cn2, wave=wave, l0=l0)
    small = (0.51 * weak * (1.0 + 0.69 * weak ** (6.0 / 5.0)) ** (-5.0 / 6.0)
             / (1.0 + 0.90 * d2 * (point / weak) ** (6.0 / 5.0)
                + 0.62 * d2 * point ** (6.0 / 5.0)))
    return np.exp(large + small) - 1.0


def averaged_index(D, wavelength, z, cn2, *, wave='plane', regime='weak',
                   spectrum='kolmogorov', l0=None, L0=None, beam=None):
    '''
    Return the irradiance flux variance sigma_I^2(D_G) behind a receive lens.

    Parameters:
        D : float or numpy.ndarray
            Hard-aperture diameter D_G [m]. The soft Gaussian radius is
            W_G = D_G / sqrt(8).
        wavelength : float or numpy.ndarray
            Optical wavelength [m].
        z : float or numpy.ndarray
            Path length L [m].
        cn2 : float or numpy.ndarray
            Refractive-index structure constant [m^-2/3], constant over the path.
        wave : str
            "plane", "spherical" or "gaussian".
        regime : str
            "weak" uses the book's weak-fluctuation closed forms, which hold for
            sigma_R^2 < 1. "strong" uses the extended Rytov chains, which the
            book states hold under ALL irradiance fluctuation conditions.
        spectrum : str
            "kolmogorov" or "modified". A finite l0 or L0 forces "modified",
            because the two-scale chains are built on that spectrum.
        l0, L0 : float, optional
            Inner scale and outer scale [m]. They need regime="strong". A finite
            L0 also needs a finite l0.
        beam : BeamParams, optional
            The beam parameters at the receive lens. Required for "gaussian".

    Returns:
        float or numpy.ndarray
            sigma_I^2(D_G).

    See the private helpers of this module for the formula and the citation of
    each branch:
        weak   plane      Ch. 10, Eq. (60), printed p. 412
               spherical  Ch. 10, Eq. (53), printed p. 409
        strong plane      Ch. 10, Eq. (69), printed p. 413 (zero scale)
                          Ch. 10, Eqs. (62) to (68), printed pp. 412 and 413
               spherical  Ch. 10, Eq. (77), printed p. 416 (zero scale)
                          Ch. 10, Eqs. (71) to (76), printed p. 415
               gaussian   Ch. 10, Eqs. (87) to (90), printed p. 420 (zero scale)
    '''
    if wave not in _WAVES:
        raise ValueError(f'wave must be one of {_WAVES}, not {wave!r}')
    if regime not in _REGIMES:
        raise ValueError(f'regime must be one of {_REGIMES}, not {regime!r}')
    if wave == 'gaussian' and beam is None:
        raise ValueError('wave="gaussian" needs beam=BeamParams(...)')

    s2 = rytov_variance(wavelength, z, cn2, wave='plane')
    d2 = d_param(D, wavelength, z) ** 2

    if l0 is not None or L0 is not None:
        if regime != 'strong':
            raise ValueError('a finite l0 or L0 needs regime="strong"')
        if spectrum not in ('modified', 'kolmogorov'):
            raise ValueError(f'spectrum must be "modified", not {spectrum!r}')
        if l0 is None:
            raise ValueError(
                'a finite L0 also needs a finite l0. Andrews and Phillips, '
                '2nd ed. (2005), DOI 10.1117/3.626196, write Ch. 10, Eqs. (62) '
                'to (68), printed pp. 412 and 413, in BOTH Q_l and Q_0, and the '
                'outer-scale term of Eq. (64) also carries Q_l.')
        if wave == 'gaussian':
            raise NotImplementedError(
                'the GAUSSIAN two-scale aperture chain is not built. Andrews '
                'and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196, Ch. 10, '
                'Eqs. (79) to (86), printed pp. 419 and 420, need eta_X of '
                'Eq. (84), printed p. 420. That equation is the same '
                'unresolved form as Ch. 9, Eq. (109), printed p. 355. The '
                'coefficient is not guessed.')
        return _strong_two_scale(D, wavelength, z, cn2, wave, l0, L0, s2)

    if spectrum != 'kolmogorov':
        raise ValueError('spectrum="modified" needs a finite l0')

    if regime == 'weak':
        if wave == 'plane':
            return _weak_plane(d2, s2)
        if wave == 'spherical':
            return _weak_spherical_factor(d2) * 0.40 * s2
        raise NotImplementedError(
            'the WEAK Gaussian-beam aperture chain is not built. Andrews and '
            'Phillips, 2nd ed. (2005), DOI 10.1117/3.626196, Ch. 10, Eq. (78), '
            'printed p. 419, is a numerical double integral, and the book '
            'prints no closed form for it. Use regime="strong", which the book '
            'states holds under all irradiance fluctuation conditions.')

    if wave == 'plane':
        return _strong_plane(d2, s2)
    if wave == 'spherical':
        return _strong_spherical(d2, 0.40 * s2)
    return _strong_gaussian(D, wavelength, z, s2, beam)


def averaging_factor(D, wavelength, z, cn2, *, wave='plane', regime='weak',
                     spectrum='kolmogorov', l0=None, L0=None, beam=None):
    '''
    Return the aperture-averaging factor A = sigma_I^2(D_G) / sigma_I^2(0).

    The parameters are the same as `averaged_index`. The point-aperture value in
    the denominator uses the SAME branch, so A goes to 1 as D_G goes to zero.
    Source: Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196, Ch. 10,
    Eq. (56), printed p. 410, which defines A as that ratio.

    For a Gaussian beam the point value uses a very small but non-zero aperture,
    because Omega_G = 16 L/(k D_G^2) is infinite at D_G = 0.
    '''
    kwargs = dict(wave=wave, regime=regime, spectrum=spectrum, l0=l0, L0=L0,
                  beam=beam)
    point_d = 0.0
    if wave == 'gaussian':
        # Omega_G must stay finite, so use one part in 10^4 of the Fresnel zone.
        point_d = 1e-4 * np.sqrt(np.asarray(z, dtype=float)
                                 / wavenumber(wavelength))
    return (averaged_index(D, wavelength, z, cn2, **kwargs)
            / averaged_index(point_d, wavelength, z, cn2, **kwargs))


if __name__ == '__main__':
    # ---------------- physics self-checks ----------------
    from .beam import beam_params

    lam_m = 1550e-9
    L = 2000.0
    cn2_weak = 3e-16
    cn2_mid = 1e-14
    s2_weak = rytov_variance(lam_m, L, cn2_weak)
    fresnel = np.sqrt(L / wavenumber(lam_m))

    # The averaging factor starts at 1 and falls as the lens grows.
    ds = np.array([1e-6, 0.01, 0.05, 0.2, 0.5])
    for wv, reg in (('plane', 'weak'), ('spherical', 'weak'),
                    ('plane', 'strong'), ('spherical', 'strong')):
        a = averaging_factor(ds, lam_m, L, cn2_weak, wave=wv, regime=reg)
        assert abs(a[0] - 1.0) < 1e-6, (wv, reg, a[0])
        assert np.all(np.diff(a) < 0.0), (wv, reg, a)
        assert np.all((a > 0.0) & (a <= 1.0)), (wv, reg, a)

    # A larger lens gives a smaller flux variance.
    idx_small = averaged_index(0.01, lam_m, L, cn2_mid, regime='strong')
    idx_big = averaged_index(0.5, lam_m, L, cn2_mid, regime='strong')
    assert idx_big < idx_small

    # The weak plane form goes to sigma_R^2 at D = 0. The rounded book constants
    # 3.86 and sin(11 pi/12) give 0.9990, not exactly 1.
    err_point = abs(_weak_plane(0.0, s2_weak) / s2_weak - 1.0)
    assert err_point < 2e-3, err_point

    # The exact spherical factor of Eq. (53) is 1 at D = 0.
    a_sp0 = _weak_spherical_factor(0.0)
    assert abs(a_sp0 - 1.0) < 2e-2, a_sp0

    # The large-aperture asymptote of the weak plane form is A ~ D^(-7/3).
    a1 = averaging_factor(2.0, lam_m, L, cn2_weak)
    a2 = averaging_factor(4.0, lam_m, L, cn2_weak)
    slope = np.log(a2 / a1) / np.log(2.0)
    assert abs(slope + 7.0 / 3.0) < 0.05, slope

    # The Gaussian chain: the flux variance goes to zero when the lens radius
    # reaches the incident beam radius, that is Omega_G -> Lambda.
    bp = beam_params(0.02, lam_m, L)
    d_match = np.sqrt(16.0 * L / (wavenumber(lam_m) * bp.lam))
    small_lens = averaged_index(0.2 * d_match, lam_m, L, cn2_weak,
                                wave='gaussian', regime='strong', beam=bp)
    near_match = averaged_index(0.999 * d_match, lam_m, L, cn2_weak,
                                wave='gaussian', regime='strong', beam=bp)
    assert near_match < 1e-3 * small_lens, (near_match, small_lens)
    # A lens wider than the beam is refused.
    try:
        averaged_index(1.5 * d_match, lam_m, L, cn2_weak, wave='gaussian',
                       regime='strong', beam=bp)
    except ValueError:
        pass
    else:
        raise AssertionError('Omega_G < Lambda must raise')

    # The two unbuilt chains are refused, not guessed.
    try:
        averaged_index(0.1, lam_m, L, cn2_weak, wave='gaussian', regime='weak',
                       beam=bp)
    except NotImplementedError:
        pass
    else:
        raise AssertionError('the weak Gaussian chain must raise')
    try:
        averaged_index(0.1, lam_m, L, cn2_weak, wave='gaussian',
                       regime='strong', spectrum='modified', l0=5e-3, beam=bp)
    except NotImplementedError:
        pass
    else:
        raise AssertionError('the Gaussian two-scale chain must raise')

    # The two-scale chain: a finite outer scale lowers the flux variance.
    kw = dict(regime='strong', spectrum='modified', l0=5e-3)
    two_no_L0 = averaged_index(0.1, lam_m, L, 1e-13, **kw)
    two_L0 = averaged_index(0.1, lam_m, L, 1e-13, L0=1.0, **kw)
    assert two_L0 < two_no_L0, (two_L0, two_no_L0)
    for wv in ('plane', 'spherical'):
        a = averaging_factor(np.array([1e-6, 0.05, 0.2]), lam_m, L, 1e-13,
                             wave=wv, **kw)
        assert abs(a[0] - 1.0) < 1e-6 and np.all(np.diff(a) < 0.0), (wv, a)

    # ---------------- REDUCTION checks ----------------
    from .. import plane_wave_scintillation as pws

    # 1. The strong plane form reproduces the parent closed form, which now
    # DELEGATES to this module.
    mine = averaged_index(0.7, lam_m, L, cn2_mid, regime='strong')
    parent = pws.aperture_averaged_index_andrews(0.7, cn2_mid, lam_m, L)
    err = abs(mine - parent)
    assert err < 1e-9, err
    print(f'REDUCTION strong plane Eq. (69) vs '
          f'aperture_averaged_index_andrews : abs err = {err:.3e}  '
          f'(target 1e-9)')

    # 2. At d = 0 the strong plane form gives the point index of Ch. 9.
    d0 = averaged_index(0.0, lam_m, L, cn2_mid, regime='strong')
    point = pws.plane_wave_scintillation_index_closed(cn2_mid, lam_m, L)
    err_d0 = abs(d0 - point)
    assert err_d0 < 1e-12, err_d0
    print(f'REDUCTION strong plane at d = 0 vs the point index : '
          f'abs err = {err_d0:.3e}  (target 1e-12)')

    # 3. The two-scale plane chain goes to the zero-scale chain as l0 -> 0 and
    # L0 -> infinity. It carries the same restriction as the point form: see the
    # note in `scintillation.large_scale_log_variance`. Measure, do not assert.
    print('MEASURED two-scale plane chain against the Eq. (69) chain at '
          'D = 0.1 m:')
    for l0_try in (1e-3, 3e-3, 5e-3, 1e-2):
        two = averaged_index(0.1, lam_m, L, 1e-13, regime='strong',
                             spectrum='modified', l0=l0_try)
        zero = averaged_index(0.1, lam_m, L, 1e-13, regime='strong')
        print(f'   l0 = {l0_try:7.1e} m : two-scale = {float(two):8.4f}  '
              f'zero-scale = {float(zero):8.4f}  '
              f'ratio = {float(two / zero):7.4f}')

    # 4. The book weak fit Eq. (61) and the exact Eq. (60) against the Churnside
    # fit that olb ships. They are DIFFERENT functions. Report, do not assert.
    print('MEASURED weak plane aperture averaging, three models:')
    print('       d     Eq. (60) exact   Eq. (61) fit   Churnside 1.07   '
          'exact/Churnside')
    for d_target in (0.5, 1.0, 2.0, 5.0):
        D_try = 2.0 * d_target * fresnel
        exact = averaging_factor(D_try, lam_m, L, cn2_weak)
        fit = plane_weak_averaging_fit(D_try, lam_m, L)
        churn = pws.aperture_averaging_factor_weak(D_try, lam_m, L)
        print(f'   {d_target:5.1f}   {float(exact):13.5f}   '
              f'{float(fit):12.5f}   {float(churn):14.5f}   '
              f'{float(exact / churn):15.4f}')

    # 5. The Gaussian chain against the plane and the spherical chains. Eq. (88)
    # is an independent fit, so the two limits do NOT agree exactly. Measure.
    bp_pl = beam_params(50.0, lam_m, L)
    bp_sp = beam_params(1e-4, lam_m, L)
    print('MEASURED Gaussian Eq. (88) chain against Eqs. (69) and (77):')
    for d_target in (0.5, 1.0, 2.0):
        D_try = 2.0 * d_target * fresnel
        g_pl = averaged_index(D_try, lam_m, L, cn2_mid, wave='gaussian',
                              regime='strong', beam=bp_pl)
        p_pl = averaged_index(D_try, lam_m, L, cn2_mid, wave='plane',
                              regime='strong')
        g_sp = averaged_index(D_try, lam_m, L, cn2_mid, wave='gaussian',
                              regime='strong', beam=bp_sp)
        p_sp = averaged_index(D_try, lam_m, L, cn2_mid, wave='spherical',
                              regime='strong')
        print(f'   d = {d_target:4.1f} : plane limit ratio = '
              f'{float(g_pl / p_pl):6.3f}   spherical limit ratio = '
              f'{float(g_sp / p_sp):6.3f}')

    print(f'A(D = 0.7 m, 2 km, weak) = '
          f'{float(averaging_factor(0.7, lam_m, L, cn2_weak)):.5f}')
    print('self-check passed')
