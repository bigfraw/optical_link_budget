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

INNER SCALE AND OUTER SCALE. The `l0` and `L0` keywords select the two-scale
branches of Ch. 9, Secs. 9.4.2 (plane, printed pp. 337 to 340) and 9.5.2
(spherical, printed pp. 343 to 345). Those branches need the wavelength and the
path length as well, because the book writes them in the nondimensional
parameters Q_l = 10.89 L/(k l0^2) and Q_0 = 64 pi^2 L/(k L0^2). So
`large_scale_log_variance` and `small_scale_log_variance` take `wavelength` and
`z` beside `sigma2_R` when a scale is set.

GAP - the GAUSSIAN two-scale STRONG branch is NOT built. Ch. 9, Eq. (109),
printed p. 355, gives the parameter eta_X of the Gaussian beam. That equation
could not be read unambiguously from the source PDF: no reading recovered gives
both the plane-wave value 2.61 (Ch. 9, Eq. (54), printed p. 339) and the
spherical-wave value 8.56 (Ch. 10, Eq. (74), printed p. 415) in the two limits.
The coefficient is NOT guessed. `wave="gaussian"` with a scale raises
NotImplementedError. The WEAK Gaussian two-scale index, Ch. 9, Eq. (104),
printed p. 354, IS built: `weak_two_scale_index` gives it, and its plane and
spherical limits are measured in the self-check.
'''

import numpy as np

from .beam import BeamParams, beam_params, effective_beam_params, wavenumber

# Weak-fluctuation boundary on the plane-wave Rytov variance. Source: Andrews
# and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196, Ch. 8, text below Eq. (23),
# printed pp. 264-265. Ch. 12, Eq. (40), printed p. 497, repeats it.
WEAK_REGIME_LIMIT = 1.0

# The three-tier regime gate reads two thresholds on the SAME axis, the Rytov
# variance sigma_R^2 (the turbulence-strength invariant). See `rytov_weak`.
#
#   sigma_R^2 <= RYTOV_CONFIDENT_WEAK   firmly weak; no warning.
#   ... < sigma_R^2 < WEAK_REGIME_LIMIT canonical weak; a SOFT warning.
#   sigma_R^2 >= WEAK_REGIME_LIMIT      leaving weak; a HARD warning.
#
# RYTOV_CONFIDENT_WEAK is a HOUSE value, not a book number: below it the Rytov
# line and the true index agree to a few percent (Andrews et al. 1999, Fig.
# behaviour, DOI 10.1364/JOSAA.16.001417; the book weak boundary is 1.0). The
# HARD limit WEAK_REGIME_LIMIT is the book boundary above.
RYTOV_CONFIDENT_WEAK = 0.3

# The uplink coupled-flux hard limit, expressed on the log-amplitude variance
# sigma_x^2. Dios et al. (DOI 10.1364/AO.43.003866) find the two-scale Gaussian
# beam-wave index stays reliable to about sigma_x^2 = 0.6 -- MORE generous than
# the book sigma_R^2 = 1 (which is sigma_x^2 = 0.25, via sigma_I^2 = 4 sigma_x^2,
# Ch. 8, Eq. (13)), because the index is a product of a large-scale and a
# small-scale factor that saturates gracefully. On the sigma_R^2 axis this is
# sigma_R^2 = 4 * 0.6 = 2.4, so the uplink passes hard_limit = 4 *
# UPLINK_SIGMA2X_LIMIT to `rytov_weak`.
UPLINK_SIGMA2X_LIMIT = 0.6

# A DISTINCT house rule, NOT a regime boundary: the largest scintillation INDEX
# sigma_I^2 for which the lognormal irradiance PDF is trusted for fade draws. It
# is 4x tighter than the Rytov regime boundary (sigma_I^2 = 1) because the
# lognormal tail goes optimistic against simulation well before the Rytov theory
# for the index itself fails (Andrews and Phillips, 2nd ed. (2005), Ch. 11, Sec.
# 11.3, printed p. 451). Keep this SEPARATE from the regime gate above: this
# gates the PDF SHAPE (lognormal vs gamma-gamma), the regime gate certifies the
# analytic INDEX. Do not conflate the two, and do not "fix" this to 1.0.
LOGNORMAL_PDF_LIMIT = 0.25

# Q_l = L kl^2/k = 10.89 L/(k l0^2), with kl = 3.3/l0 the inner-scale wavenumber
# of the modified atmospheric spectrum. Source: Andrews and Phillips, 2nd ed.
# (2005), DOI 10.1117/3.626196, Ch. 9, text below Eq. (48), printed p. 338, and
# Ch. 10, Eq. (68), printed p. 413.
QL_CONSTANT = 10.89

# Q_0 = L k0^2/k = 64 pi^2 L/(k L0^2), so this branch uses k0 = 8*pi/L0. Source:
# Ch. 9, text below Eq. (57), printed p. 339, and Ch. 10, Eq. (68), printed
# p. 413.
Q0_CONSTANT = 64.0 * np.pi ** 2

_WAVES = ('plane', 'spherical', 'gaussian')


def two_scale_parameters(wavelength, z, l0=None, L0=None):
    '''
    Return the nondimensional scale parameters (Q_l, Q_0).

    formula:
        Q_l = 10.89 L / (k l0^2),   Q_0 = 64 pi^2 L / (k L0^2),   k = 2*pi/lambda
    Source: Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196, Ch. 9,
    text below Eq. (48), printed p. 338, and text below Eq. (57), printed p. 339.
    Ch. 10, Eq. (68), printed p. 413, prints both together.

    A missing inner scale gives Q_l = infinity, which is the zero-inner-scale
    limit. A missing outer scale gives Q_0 = 0, which is the infinite-outer-scale
    limit.
    '''
    k = wavenumber(wavelength)
    z = np.asarray(z, dtype=float)
    ql = np.inf if l0 is None else QL_CONSTANT * z / (k * float(l0) ** 2)
    q0 = 0.0 if L0 is None else Q0_CONSTANT * z / (k * float(L0) ** 2)
    return ql, q0


def _need_path(wavelength, z, l0, L0):
    '''Refuse a two-scale call that gives no path.'''
    if wavelength is None or z is None:
        raise ValueError('a finite l0 or L0 also needs wavelength and z, '
                         'because the book writes the two-scale forms in '
                         'Q_l = 10.89 L/(k l0^2) and Q_0 = 64 pi^2 L/(k L0^2)')
    return two_scale_parameters(wavelength, z, l0, L0)


def _eta_filter(eta, ql):
    '''
    Return the shared large-scale filter group of the two-scale theory.

    formula:
        [ eta Q_l / (eta + Q_l) ]^(7/6)
        [ 1 - 1.75 (eta/(eta+Q_l))^(1/2) + 0.25 (eta/(eta+Q_l))^(7/12) ]
    Source: Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196, Ch. 9,
    Eq. (53), printed p. 338. Eqs. (55), (56), (80) and (81), printed pp. 339,
    340, 419 and 420, all reuse it. Appendix III, Tables VII(b) and VIII(b),
    printed pp. 769 and 770, print the same group.
    '''
    ratio = eta / (eta + ql)
    return ((eta * ql / (eta + ql)) ** (7.0 / 6.0)
            * (1.0 - 1.75 * ratio ** 0.5 + 0.25 * ratio ** (7.0 / 12.0)))


def weak_two_scale_index(wavelength, z, cn2, *, wave='plane', l0=None,
                         beam=None):
    '''
    Return the WEAK-fluctuation scintillation index on the modified spectrum.

    This is the input that the small-scale log variance needs when the inner
    scale is finite. With l0 = None it reduces to the Kolmogorov Rytov variance.

    Parameters:
        wavelength : float or numpy.ndarray
            Optical wavelength [m].
        z : float or numpy.ndarray
            Path length L [m].
        cn2 : float or numpy.ndarray
            Refractive-index structure constant [m^-2/3].
        wave : str
            "plane", "spherical" or "gaussian".
        l0 : float, optional
            Inner scale [m]. None gives the Kolmogorov limit.
        beam : BeamParams, optional
            The beam parameters at the receiver. Required for "gaussian".

    Returns:
        float or numpy.ndarray
            sigma_PL^2 (plane), sigma_SP^2 (spherical) or sigma_G^2 (Gaussian).

    formula (plane, Q_l = 10.89 L/(k l0^2)):
        sigma_PL^2 = 3.86 sigma_R^2 { (1 + 1/Q_l^2)^(11/12)
            [ sin((11/6) atan Q_l)
              + 1.51 (1+Q_l^2)^(-1/4) sin((4/3) atan Q_l)
              - 0.27 (1+Q_l^2)^(-7/24) sin((5/4) atan Q_l) ]
            - 3.50 Q_l^(-5/6) }
    Source: Ch. 9, Eq. (48), printed p. 338. Appendix III, Table VII(a), printed
    p. 769, prints the same row with 1.507 and 0.273.

    formula (spherical):
        sigma_SP^2 = 9.65 beta_0^2 { 0.40 (1 + 9/Q_l^2)^(11/12)
            [ sin((11/6) atan(Q_l/3))
              + 2.61 (9+Q_l^2)^(-1/4) sin((4/3) atan(Q_l/3))
              - 0.518 (9+Q_l^2)^(-7/24) sin((5/4) atan(Q_l/3)) ]
            - 3.50 Q_l^(-5/6) }
    Source: Ch. 9, Eq. (75), printed p. 343. Appendix III, Table VIII(a), printed
    p. 770, prints the same row.

    formula (Gaussian, with phi_1 and phi_2 of Ch. 9, Eq. (105)):
        sigma_G^2 = 3.86 sigma_R^2 {
            0.40 [(1+2 Theta)^2 + (2 Lambda + 3/Q_l)^2]^(11/12)
                 [(1+2 Theta)^2 + 4 Lambda^2]^(-1/2)
            [ sin((11/6) phi_2 + phi_1)
              + 2.61 [(1+2 Theta)^2 Q_l^2 + (3 + 2 Lambda Q_l)^2]^(-1/4)
                sin((4/3) phi_2 + phi_1)
              - 0.52 [(1+2 Theta)^2 Q_l^2 + (3 + 2 Lambda Q_l)^2]^(-7/24)
                sin((5/4) phi_2 + phi_1) ]
            - 13.40 Lambda / (Q_l^(11/6) [(1+2 Theta)^2 + 4 Lambda^2])
            - (11/6) Q_l^(-5/6) [ (1 + 0.31 Lambda Q_l)^(5/6)
                                  + 1.10 (1 + 0.27 Lambda Q_l)^(1/3)
                                  - 0.19 (1 + 0.24 Lambda Q_l)^(1/4) ] }
        phi_1 = atan[2 Lambda / (1 + 2 Theta)]
        phi_2 = atan[(1 + 2 Theta) Q_l / (3 + 2 Lambda Q_l)]
    Source: Ch. 9, Eqs. (104) and (105), printed pp. 354 and 355.

    The Gaussian row reduces to the plane row at Theta = 1, Lambda = 0, to the
    spherical row at Theta = Lambda = 0, and to Ch. 8, Eq. (23) as
    Q_l -> infinity. The module self-check measures all three reductions.
    '''
    _check_wave(wave)
    bm = _need_beam(beam, wave)
    sigma2_R = rytov_variance(wavelength, z, cn2, wave='plane')
    if l0 is None:
        if wave == 'plane':
            return sigma2_R
        if wave == 'spherical':
            return 0.40 * sigma2_R
        return beam_rytov_variance(sigma2_R, bm)

    ql, _ = two_scale_parameters(wavelength, z, l0=l0)
    if wave == 'plane':
        t = np.arctan(ql)
        inner = (np.sin((11.0 / 6.0) * t)
                 + 1.51 * (1.0 + ql ** 2) ** (-0.25) * np.sin((4.0 / 3.0) * t)
                 - 0.27 * (1.0 + ql ** 2) ** (-7.0 / 24.0)
                 * np.sin((5.0 / 4.0) * t))
        return 3.86 * sigma2_R * ((1.0 + 1.0 / ql ** 2) ** (11.0 / 12.0) * inner
                                  - 3.50 * ql ** (-5.0 / 6.0))
    if wave == 'spherical':
        t = np.arctan(ql / 3.0)
        inner = (np.sin((11.0 / 6.0) * t)
                 + 2.61 * (9.0 + ql ** 2) ** (-0.25) * np.sin((4.0 / 3.0) * t)
                 - 0.518 * (9.0 + ql ** 2) ** (-7.0 / 24.0)
                 * np.sin((5.0 / 4.0) * t))
        return 9.65 * (0.40 * sigma2_R) * (
            0.40 * (1.0 + 9.0 / ql ** 2) ** (11.0 / 12.0) * inner
            - 3.50 * ql ** (-5.0 / 6.0))

    th, lm = bm.theta, bm.lam
    p = 1.0 + 2.0 * th
    base = p ** 2 + 4.0 * lm ** 2
    mixed = p ** 2 * ql ** 2 + (3.0 + 2.0 * lm * ql) ** 2
    phi1 = np.arctan2(2.0 * lm, p)
    phi2 = np.arctan2(p * ql, 3.0 + 2.0 * lm * ql)
    inner = (np.sin((11.0 / 6.0) * phi2 + phi1)
             + 2.61 * mixed ** (-0.25) * np.sin((4.0 / 3.0) * phi2 + phi1)
             - 0.52 * mixed ** (-7.0 / 24.0)
             * np.sin((5.0 / 4.0) * phi2 + phi1))
    lq = lm * ql
    tail = (11.0 / 6.0) * ql ** (-5.0 / 6.0) * (
        (1.0 + 0.31 * lq) ** (5.0 / 6.0)
        + 1.10 * (1.0 + 0.27 * lq) ** (1.0 / 3.0)
        - 0.19 * (1.0 + 0.24 * lq) ** 0.25)
    return 3.86 * sigma2_R * (
        0.40 * (p ** 2 + (2.0 * lm + 3.0 / ql) ** 2) ** (11.0 / 12.0)
        * base ** (-0.5) * inner
        - 13.40 * lm / (ql ** (11.0 / 6.0) * base)
        - tail)


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

    SPECTRUM: Kolmogorov. The 7/6 and 11/6 exponents (and the 1.23 constant) are
    the signature of the -11/3 power law with a zero inner scale. This is the
    ROOT of the weak theory: nearly every weak-regime quantity in the layer is
    sigma_R^2 times a geometric factor, so it inherits Kolmogorov from here. For
    a finite inner or outer scale use `weak_two_scale_index`.

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


def rytov_weak(sigma2_R, Lambda=None, *, hard_limit=WEAK_REGIME_LIMIT,
               soft_limit=RYTOV_CONFIDENT_WEAK):
    '''
    Classify the fluctuation regime as "weak", "soft", or "hard".

    This is the one shared weak-fluctuation gate. It reads the Rytov variance
    (the turbulence-strength invariant) against two thresholds and returns a
    label the caller turns into no warning, a soft warning, or a hard warning.

        sigma_R^2 <= soft_limit    -> "weak"   firmly weak; trust the model.
        soft_limit < ... < hard    -> "soft"   canonical weak; a soft warning.
        sigma_R^2 >= hard_limit    -> "hard"   leaving weak; a hard warning.

    GAUSSIAN BEAM. The plane-wave threshold is not adequate for a beam wave.
    Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196, Ch. 5, Eq. (16),
    printed p. 140, need BOTH sigma_R^2 < 1 AND sigma_R^2 Lambda^(5/6) < 1, with
    Lambda the OUTPUT-plane (receiver) Gaussian-beam parameter. Pass `Lambda` and
    the gate reads the binding strength

        s = sigma_R^2 * max(1, Lambda^(5/6)),

    so a FOCUSED beam (Lambda > 1) trips the gate before sigma_R^2 alone would,
    and a collimated or divergent beam (Lambda <= 0.5) is unchanged (the second
    condition is looser and never binds first).

    The two limits are keyword-only so a caller on a different axis can move
    them. The uplink works on the log-amplitude variance sigma_x^2 (with
    sigma_R^2 = 4 sigma_x^2, Ch. 8, Eq. (13)) and the Dios two-scale index is
    reliable to sigma_x^2 = UPLINK_SIGMA2X_LIMIT; it calls this with the strength
    already on the sigma_R^2 axis and hard_limit = 4 * UPLINK_SIGMA2X_LIMIT.

    NOTE. This gate certifies the analytic INDEX (is Rytov theory valid?). It is
    NOT the lognormal-PDF house rule LOGNORMAL_PDF_LIMIT, which gates the fade
    PDF SHAPE on the scintillation index sigma_I^2. Keep the two separate.

    Parameters:
        sigma2_R : float
            The strength on the Rytov-variance axis. Scalar.
        Lambda : float, optional
            The output-plane Gaussian-beam parameter. None for a plane wave (no
            beam correction).
        hard_limit : float
            The hard-warning threshold (default WEAK_REGIME_LIMIT = 1.0, the
            book boundary).
        soft_limit : float
            The soft-warning threshold (default RYTOV_CONFIDENT_WEAK = 0.3).

    Returns:
        str
            "weak", "soft", or "hard".
    '''
    s = float(sigma2_R)
    if Lambda is not None:
        s *= max(1.0, float(Lambda) ** (5.0 / 6.0))
    if s >= hard_limit:
        return 'hard'
    if s > soft_limit:
        return 'soft'
    return 'weak'


def large_scale_log_variance(sigma2_R, *, wave='plane', l0=None, L0=None,
                             beam=None, r=0.0, wavelength=None, z=None):
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
        l0, L0 : float, optional
            Inner scale and outer scale [m]. A value selects the two-scale
            branch and then `wavelength` and `z` are required as well. The
            Gaussian two-scale branch is not built; see the module docstring.
        beam : BeamParams, optional
            The beam parameters at the receiver. Required for "gaussian".
        r : float
            Off-axis radius [m]. It has NO effect. The book splits only the
            LONGITUDINAL component into a large-scale and a small-scale part.
            The radial component stays separate (Ch. 9, Eq. (103), printed
            p. 353). The keyword is here to match `scintillation_index`.
        wavelength, z : float, optional
            Optical wavelength [m] and path length L [m]. Required only for the
            two-scale branch, which needs Q_l and Q_0.

    Returns:
        float or numpy.ndarray
            sigma_lnX^2.

    formula (zero inner scale, infinite outer scale):
        plane      0.49 s^2 / (1 + 1.11 s^(12/5))^(7/6),   s^2 = sigma_R^2
        spherical  0.20 s^2 / (1 + 0.19 s^(12/5))^(7/6)
        gaussian   0.49 b^2 / (1 + 0.56 (1 + Theta) b^(12/5))^(7/6),
                   b^2 = sigma_B^2
    Source: Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196:
        plane      Ch. 9, Eq. (41), printed p. 335
        spherical  Ch. 9, Eq. (69), printed p. 342
        gaussian   Ch. 9, Eq. (97), printed p. 352

    formula (two scales, the DIFFERENCE of an inner-scale and an outer-scale
    term, Ch. 9, Eq. (51), printed p. 338):
        sigma_lnX^2 = C F(eta_X) - C F(eta_X0)
        F(eta)   = [eta Q_l/(eta+Q_l)]^(7/6)
                   [1 - 1.75 (eta/(eta+Q_l))^(1/2)
                      + 0.25 (eta/(eta+Q_l))^(7/12)]
        eta_X0   = eta_X Q_0 / (eta_X + Q_0)
        plane      C = 0.16 sigma_R^2,  eta_X = 2.61/(1+0.45 sigma_R^2 Q_l^(1/6))
        spherical  C = 0.04 beta_0^2,   eta_X = 8.56/(1+0.20 beta_0^2 Q_l^(1/6))
    Source: plane, Ch. 9, Eqs. (51) to (57), printed pp. 338 and 339; spherical,
    Ch. 10, Eqs. (72) to (76), printed p. 415, taken at d = 0, which the book
    also prints in Appendix III, Table VIII(b), printed p. 770.

    NOTE on the zero-inner-scale limit. As Q_l goes to infinity the two-scale
    plane branch goes to 0.16 sigma_R^2 [2.61/(1+0.45 sigma_R^2 Q_l^(1/6))]^(7/6),
    NOT to the Kolmogorov branch. The two agree only where
    0.45 sigma_R^2 Q_l^(1/6) equals 1.11 sigma_R^(12/5). The reason is in the
    book: Ch. 9, Eq. (54), printed p. 339, states the substitution
    L/(k rho_0^2) = 1.02 sigma_R^2 Q_l^(1/6) for the case rho_0 << l0 ONLY. So
    the two-scale branch is a MODERATE-to-STRONG turbulence model with a real
    inner scale, not a superset of the Kolmogorov branch. The module self-check
    measures the gap.
    '''
    _check_wave(wave)
    s2 = np.asarray(sigma2_R, dtype=float)
    if l0 is not None or L0 is not None:
        if wave == 'gaussian':
            raise NotImplementedError(
                'the GAUSSIAN two-scale large-scale branch is not built. '
                'Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196, '
                'Ch. 9, Eq. (109), printed p. 355, gives the parameter eta_X, '
                'but that equation could not be read unambiguously from the '
                'source PDF: no recovered reading gives both the plane-wave '
                'value 2.61 (Ch. 9, Eq. (54), printed p. 339) and the '
                'spherical-wave value 8.56 (Ch. 10, Eq. (74), printed p. 415) '
                'in the two limits. The coefficient is not guessed.')
        ql, q0 = _need_path(wavelength, z, l0, L0)
        if wave == 'plane':
            coef, top, slope = 0.16 * s2, 2.61, 0.45 * s2
        else:
            b0 = 0.40 * s2
            coef, top, slope = 0.04 * b0, 8.56, 0.20 * b0
        eta_x = top / (1.0 + slope * ql ** (1.0 / 6.0))
        out = coef * _eta_filter(eta_x, ql)
        if L0 is not None:
            eta_x0 = eta_x * q0 / (eta_x + q0)
            out = out - coef * _eta_filter(eta_x0, ql)
        return out
    if wave == 'plane':
        return 0.49 * s2 / (1.0 + 1.11 * s2 ** (6.0 / 5.0)) ** (7.0 / 6.0)
    if wave == 'spherical':
        return 0.20 * s2 / (1.0 + 0.19 * s2 ** (6.0 / 5.0)) ** (7.0 / 6.0)
    bm = _need_beam(beam, wave)
    b2 = beam_rytov_variance(s2, bm)
    denom = (1.0 + 0.56 * (1.0 + bm.theta) * b2 ** (6.0 / 5.0)) ** (7.0 / 6.0)
    return 0.49 * b2 / denom


def small_scale_log_variance(sigma2_R, *, wave='plane', l0=None, L0=None,
                             beam=None, r=0.0, wavelength=None, z=None):
    '''
    Return the small-scale log-irradiance variance sigma_lnY^2.

    This is the second of the two variances that the extended Rytov theory
    needs. See `large_scale_log_variance` for the parameters. The keyword `r`
    has no effect for the same reason.

    formula (zero inner scale, infinite outer scale):
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

    formula (finite inner scale):
        plane      0.51 sigma_PL^2 / (1 + 0.69 sigma_PL^(12/5))^(5/6)
        spherical  0.51 sigma_SP^2 / (1 + 0.69 sigma_SP^(12/5))^(5/6)
        gaussian   0.51 sigma_G^2  / (1 + 0.69 sigma_G^(12/5))^(5/6)
    with the weak two-scale index from `weak_two_scale_index`. Source: Ch. 9,
    Eqs. (59) and (60), printed p. 340 (plane); Eqs. (82) and (83), printed
    p. 345 (spherical); Eqs. (111) and (112), printed p. 355 (Gaussian). The
    OUTER scale has a negligible effect here, which the book states below
    Eq. (60), printed p. 340. So a value of L0 changes nothing in this function.
    '''
    _check_wave(wave)
    s2 = np.asarray(sigma2_R, dtype=float)
    if l0 is not None:
        if wavelength is None or z is None:
            raise ValueError('a finite l0 also needs wavelength and z')
        # Rebuild Cn2 from the plane-wave Rytov variance, so that the caller can
        # keep the sigma2_R interface. sigma_R^2 = 1.23 Cn2 k^(7/6) L^(11/6).
        k = wavenumber(wavelength)
        cn2 = s2 / (1.23 * k ** (7.0 / 6.0)
                    * np.asarray(z, dtype=float) ** (11.0 / 6.0))
        idx = weak_two_scale_index(wavelength, z, cn2, wave=wave, l0=l0,
                                   beam=beam)
        return 0.51 * idx / (1.0 + 0.69 * idx ** (6.0 / 5.0)) ** (5.0 / 6.0)
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
        l0, L0 : float, optional
            Inner scale and outer scale [m]. A value selects the two-scale
            branch on the modified atmospheric spectrum. The Gaussian STRONG
            two-scale branch is not built; see the module docstring.
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
    bm = _need_beam(beam, wave)

    sigma2_R = rytov_variance(wavelength, z, cn2, wave='plane')

    def weak():
        point = weak_two_scale_index(wavelength, z, cn2, wave=wave, l0=l0,
                                     beam=bm)
        if wave != 'gaussian':
            return point
        radial = _radial_component(sigma2_R, bm.lam, bm.w, r, tracked,
                                   pointing_error_m, wander_rms_m)
        return radial + point

    def strong():
        x = large_scale_log_variance(sigma2_R, wave=wave, beam=bm, l0=l0,
                                     L0=L0, wavelength=wavelength, z=z)
        y = small_scale_log_variance(sigma2_R, wave=wave, beam=bm, l0=l0,
                                     L0=L0, wavelength=wavelength, z=z)
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

    # ---------------- the three-tier weak gate ----------------
    # Plane wave (no beam correction): the tiers sit at 0.3 and 1.0 on sigma_R^2.
    assert rytov_weak(0.2) == 'weak'
    assert rytov_weak(0.3) == 'weak'            # boundary is inclusive-weak
    assert rytov_weak(0.5) == 'soft'
    assert rytov_weak(1.0) == 'hard'
    assert rytov_weak(3.0) == 'hard'
    # A collimated / divergent beam has Lambda <= 0.5, so Lambda^(5/6) < 1 and
    # the beam condition never binds: the label is the plane-wave label.
    assert rytov_weak(0.5, Lambda=0.4) == 'soft'
    assert rytov_weak(0.5, Lambda=0.0) == 'soft'
    # A FOCUSED beam (Lambda > 1) trips the gate before sigma_R^2 alone does:
    # sigma_R^2 = 0.5 is "soft" as a plane wave but "hard" at Lambda = 4
    # (0.5 * 4^(5/6) = 1.66 >= 1). This is the TL-05 fix in one line.
    assert rytov_weak(0.5) == 'soft'
    assert rytov_weak(0.5, Lambda=4.0) == 'hard', 0.5 * 4.0 ** (5.0 / 6.0)
    # The uplink axis: pass the strength as sigma_R^2 = 4 sigma_x^2 and raise the
    # hard limit to 4 * UPLINK_SIGMA2X_LIMIT = 2.4. The Dios edge sigma_x^2 = 0.6
    # is then "hard", and the book sigma_x^2 = 0.25 (sigma_R^2 = 1) is "soft".
    up_hard = 4.0 * UPLINK_SIGMA2X_LIMIT
    assert rytov_weak(4.0 * 0.25, hard_limit=up_hard) == 'soft'
    assert rytov_weak(4.0 * 0.6, hard_limit=up_hard) == 'hard'
    assert rytov_weak(4.0 * 0.05, hard_limit=up_hard) == 'weak'
    # The two house rules are distinct numbers, not the same gate.
    assert LOGNORMAL_PDF_LIMIT == 0.25 and WEAK_REGIME_LIMIT == 1.0
    assert RYTOV_CONFIDENT_WEAK == 0.3 and UPLINK_SIGMA2X_LIMIT == 0.6
    print('[gate] rytov_weak tiers (plane, collimated, focused, uplink) ok')

    # The GAUSSIAN two-scale STRONG branch is refused, not guessed.
    try:
        large_scale_log_variance(s2_R, wave='gaussian',
                                 beam=beam_params(0.05, lam_m, L), l0=5e-3,
                                 wavelength=lam_m, z=L)
    except NotImplementedError:
        pass
    else:
        raise AssertionError('the Gaussian two-scale branch must raise')

    # A two-scale call with no path is refused.
    try:
        large_scale_log_variance(s2_R, wave='plane', l0=5e-3)
    except ValueError:
        pass
    else:
        raise AssertionError('l0 without wavelength and z must raise')

    # The weak two-scale index. A finite inner scale RAISES the weak plane index
    # above the Kolmogorov value, because the spectral bump adds power near
    # 1/l0 (Ch. 3, Sec. 3.3.3, printed p. 68).
    l0_ref = 5e-3
    idx_l0 = weak_two_scale_index(lam_m, L, cn2_weak, wave='plane', l0=l0_ref)
    assert idx_l0 > s2_R, (idx_l0, s2_R)
    # A finite outer scale lowers the strong plane index.
    strong_no_L0 = scintillation_index(lam_m, L, 1e-13, wave='plane',
                                       regime='strong', l0=l0_ref)
    strong_L0 = scintillation_index(lam_m, L, 1e-13, wave='plane',
                                    regime='strong', l0=l0_ref, L0=1.0)
    assert strong_L0 < strong_no_L0, (strong_L0, strong_no_L0)

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
    # 7. The two-scale branches against the zero-scale branches.
    #    (a) The weak two-scale index goes to the Kolmogorov index as l0 -> 0.
    for wv, ref in (('plane', s2_R), ('spherical', 0.40 * s2_R)):
        tiny = weak_two_scale_index(lam_m, L, cn2_weak, wave=wv, l0=1e-9)
        pct = abs(tiny - ref) / ref * 100.0
        assert pct < 1.0, (wv, pct)
        print(f'REDUCTION weak_two_scale_index({wv}, l0->0) : {pct:.4f} % '
              f'(target 1 %)')
    #    The Gaussian row also goes to the plane row, the spherical row and the
    #    Ch. 8, Eq. (23) row.
    g_pl = weak_two_scale_index(lam_m, L, cn2_weak, wave='gaussian',
                                l0=1e-9, beam=bp_plane)
    g_sp = weak_two_scale_index(lam_m, L, cn2_weak, wave='gaussian',
                                l0=1e-9, beam=bp_sph)
    pct_gpl = abs(g_pl - s2_R) / s2_R * 100.0
    pct_gsp = abs(g_sp - 0.40 * s2_R) / (0.40 * s2_R) * 100.0
    assert pct_gpl < 1.0 and pct_gsp < 1.0, (pct_gpl, pct_gsp)
    g_beam = weak_two_scale_index(lam_m, L, cn2_weak, wave='gaussian',
                                  l0=1e-9, beam=bp)
    pct_gb = abs(g_beam - beam_rytov_variance(s2_R, bp)) / abs(
        beam_rytov_variance(s2_R, bp)) * 100.0
    assert pct_gb < 1.0, pct_gb
    print(f'REDUCTION Ch. 9, Eq. (104) (l0->0) : plane {pct_gpl:.4f} %  '
          f'spherical {pct_gsp:.4f} %  collimated beam vs Ch. 8, Eq. (23) '
          f'{pct_gb:.4f} %  (target 1 %)')

    #    (b) The outer-scale term goes to zero as L0 -> infinity.
    s2_mid = rytov_variance(lam_m, L, 1e-13)
    with_L0 = large_scale_log_variance(s2_mid, wave='plane', l0=l0_ref,
                                       L0=1e6, wavelength=lam_m, z=L)
    no_L0 = large_scale_log_variance(s2_mid, wave='plane', l0=l0_ref,
                                     wavelength=lam_m, z=L)
    pct_L0 = abs(with_L0 - no_L0) / no_L0 * 100.0
    assert pct_L0 < 1.0, pct_L0
    print(f'REDUCTION outer-scale term (L0 -> inf) : {pct_L0:.4f} % '
          f'(target 1 %)')

    #    (c) The large-scale two-scale branch does NOT go to the Kolmogorov
    #    branch as l0 -> 0. Measure the gap and report it. Ch. 9, Eq. (54),
    #    printed p. 339, restricts its substitution to rho_0 << l0, so the two
    #    branches agree only where 0.45 sigma_R^2 Q_l^(1/6) = 1.11 sigma_R^(12/5).
    kol_x = large_scale_log_variance(s2_mid, wave='plane')
    print('MEASURED large-scale two-scale gap against the Kolmogorov branch '
          f'(sigma_R^2 = {float(s2_mid):.3f}):')
    for l0_try in (1e-3, 3e-3, 5e-3, 1e-2, 1e-9):
        two = large_scale_log_variance(s2_mid, wave='plane', l0=l0_try,
                                       wavelength=lam_m, z=L)
        ql, _ = two_scale_parameters(lam_m, L, l0=l0_try)
        print(f'   l0 = {l0_try:9.1e} m   Q_l = {float(ql):11.3e}   '
              f'two-scale/Kolmogorov = {float(two / kol_x):7.4f}')

    print(f'GAP 9 divergent f0={f0_div} m : '
          f'Andrews Eq. (23) = {mine_div:.6f}  Dios = {dios_div:.6f}  '
          f'diff = {gap9_div:+.2f} %  (no assert)')

    print('self-check passed')
