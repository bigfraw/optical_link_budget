'''
Slant-path and satellite-link forms of Andrews and Phillips, Chapter 12.

This module holds the ground-to-space and space-to-ground physics. It gives:
    the slant secant                            Ch. 12, Eq. (14)
    the Hufnagel-Valley Cn2(h) profile          Ch. 12, Eq. (1)
    the Bufton wind and its rms                 Ch. 12, Eqs. (2) and (3)
    the SCIDAR outer-scale profiles             Ch. 12, Eqs. (6) and (7)
    the path moments mu_0 to mu_3               Ch. 12, Eqs. (18)-(21), (25),
                                                (26), (37), (55)
    the downlink scintillation index            Ch. 12, Eqs. (38), (39), (40)
    the uplink scintillation index              Ch. 12, Eqs. (54), (56)-(61)
    the uplink spatial coherence radius         Ch. 12, Eqs. (24)-(27)
    the isoplanatic angle                       Ch. 12, Eqs. (29) and (30)
    the point-ahead angle                       Ch. 12, Sec. 12.3.3

Source of every equation: Andrews and Phillips, Laser Beam Propagation through
Random Media, 2nd ed. (SPIE Press, 2005), DOI 10.1117/3.626196.

This module holds physics only. It returns no decibels. It imports numpy, its
sibling andrews modules, and olb._deps. It imports no scenario, no terminal, no
Term and no link.

PLANE OF REFERENCE. The book uses one normalised path variable for both link
directions, Ch. 12, Eq. (14), printed p. 490:
    uplink     xi = 1 - (h - h0)/(H - h0)
    downlink   xi = (h - h0)/(H - h0)
where h0 is the ground height and H is the satellite altitude. In BOTH cases
xi = 1 at the TRANSMITTER and xi = 0 at the RECEIVER, so xi = 1 - z/L with z
measured from the transmitter. Read this before you compare a path weight in
this module with a path weight in another module. See Conflict C-02 in
docs/andrews-crosscheck.md.

GEOMETRY LIMIT. The book uses a plane-parallel atmosphere. It writes the
satellite altitude as H = h0 + L cos(zeta) (Ch. 12, text below Eq. (14), printed
p. 490) and puts sec(zeta) in front of each path integral. It gives NO
Earth-curvature correction. It also limits the weak-fluctuation results to
zenith angles that do not exceed 45 to 60 deg (Ch. 12, Sec. 12.1, printed
p. 478, and Ch. 12, Sec. 12.9, printed p. 521). ZENITH_LIMIT_DEG holds that
bound.

STRENGTH LIMIT. Each weak form below states its own regime. The strong forms
Ch. 12, Eq. (40) (downlink) and Ch. 12, Eqs. (59) to (61) (uplink) hold for all
values of the Rytov variance. The beam-wander part of the uplink forms stays a
weak-fluctuation result, which the book states below Eq. (59), printed p. 506.

SPECTRUM. Chapter 12 uses the Kolmogorov spectrum only (Ch. 12, Eq. (15),
printed p. 490). It gives no slant-path form with an inner scale or an outer
scale. So an l0 or an L0 argument below is REFUSED, not approximated. Use
olb.turbulence.andrews.scintillation.weak_two_scale_index for the two-scale
forms on a single homogeneous path.
'''

import numpy as np

from ..._deps import get_c2n, v_wind
from .beam import wavenumber
from .scintillation import (WEAK_REGIME_LIMIT, large_scale_log_variance,
                            small_scale_log_variance)
from .wander import (beam_wander_variance_slant, plane_fried_parameter_slant,
                     pointing_error_variance_slant)

# Zenith angle above which the book does not trust the weak-fluctuation slant
# results. Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196, Ch. 12,
# Sec. 12.1, printed p. 478: weak theory is "sufficient ... provided the zenith
# angle is sufficiently small (less than 60 deg in most cases but may be
# restricted to zenith angles less than 45 deg in cases where ground-level Cn2
# is large)". Ch. 12, Sec. 12.9, printed p. 521, repeats the 45 to 60 deg bound.
# The same pages carry the plane-parallel geometry: the book never corrects
# sec(zeta) for the curvature of the Earth.
ZENITH_LIMIT_DEG = 60.0

# Scintillation-index constant of the downlink and uplink path integrals, and
# the point plane-wave constant that it must reproduce. Andrews and Phillips,
# 2nd ed. (2005), DOI 10.1117/3.626196, Ch. 12, Eqs. (36), (38), (39), (54),
# (57) and (58), printed pp. 495, 496, 503 and 504. The book rounds 2.2517 to
# 2.25, because 8.70 cos(5 pi/12) = 2.2517.
_MU3_CONSTANT = 8.70
_POINT_CONSTANT = 2.25

# Constant of the uplink beam-wander scintillation term. Andrews and Phillips,
# 2nd ed. (2005), DOI 10.1117/3.626196, Ch. 12, Eqs. (54), (56), (57), (59) and
# (61), printed pp. 503, 504 and 506.
_WANDER_SCINT_CONSTANT = 5.95

# Constant of the wave structure function, and the coherence-radius constant
# that follows from D(rho_0) = 2. Andrews and Phillips, 2nd ed. (2005),
# DOI 10.1117/3.626196, Ch. 12, Eqs. (17), (20), (22), (24) and (27), printed
# pp. 491 and 492. The book prints 1.45 in Eq. (22) and 1.46 in Eq. (27); both
# are 2.91/2 = 1.455 rounded.
_WSF_CONSTANT = 2.91
_RHO0_CONSTANT = 1.46

# Weight of the Lambda^(11/6) term of the wave structure function. Andrews and
# Phillips, 2nd ed. (2005), DOI 10.1117/3.626196, Ch. 12, Eqs. (17) and (24),
# printed pp. 491 and 492.
_LAMBDA_TERM = 0.62

# Isoplanatic-angle constant. Andrews and Phillips, 2nd ed. (2005),
# DOI 10.1117/3.626196, Ch. 12, Eqs. (29) and (30), printed p. 493.
_ISOPLANATIC_CONSTANT = 2.91

# Speed of light [m/s], for the point-ahead angle.
_C = 299792458.0

_DIRECTIONS = ('uplink', 'downlink')


def sec_zeta(elevation_deg):
    '''
    Return the slant secant sec(zeta) of the plane-parallel atmosphere.

    Parameters:
        elevation_deg : float or numpy.ndarray
            Elevation angle above the horizon [deg]. 90 is the zenith.

    Returns:
        float or numpy.ndarray
            sec(zeta) = 1/sin(elevation).

    formula:
        z = (h - h0) sec(zeta),   dz = sec(zeta) dh,   H = h0 + L cos(zeta)
    Source: Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196, Ch. 12,
    text at Eq. (14), printed p. 490. Every Ch. 12 path integral carries this
    factor, for example Eqs. (11) to (13), printed pp. 489 and 490.

    VALIDITY. The book gives no Earth-curvature correction, and it limits the
    weak-fluctuation slant results to zenith angles that do not exceed 45 to
    60 deg (Ch. 12, Sec. 12.1, printed p. 478; Ch. 12, Sec. 12.9, printed
    p. 521). See ZENITH_LIMIT_DEG. This function does NOT refuse a low
    elevation, because the strong-fluctuation forms below stay valid there. It
    refuses only an elevation at or below the horizon, where the plane-parallel
    geometry has no meaning.
    '''
    elevation = np.asarray(elevation_deg, dtype=float)
    if np.any(elevation <= 0.0) or np.any(elevation > 90.0):
        raise ValueError('elevation_deg must be in the range (0, 90]')
    return 1.0 / np.sin(np.radians(elevation))


def hufnagel_valley(h, wind_rms=21.0, cn2_ground=1.7e-14):
    '''
    Return the Hufnagel-Valley Cn2(h) profile [m^-2/3].

    Parameters:
        h : float or numpy.ndarray
            Altitude above ground level [m].
        wind_rms : float
            The rms high-altitude wind speed (pseudowind) w [m/s]. See
            `rms_wind`.
        cn2_ground : float
            The nominal ground value A = Cn2(0) [m^-2/3].

    Returns:
        numpy.ndarray
            Cn2(h) [m^-2/3].

    formula:
        Cn2(h) = 0.00594 (w/27)^2 (1e-5 h)^10 exp(-h/1000)
                 + 2.7e-16 exp(-h/1500)
                 + A exp(-h/100)
    Source: Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196, Ch. 12,
    Eq. (1), printed p. 481. The default pair w = 21 m/s and A = 1.7e-14 is the
    H-V5/7 model that the book uses through the whole chapter (Ch. 12, text
    below Eq. (3), printed p. 481).

    This function DELEGATES to the shared kernel `get_c2n`, which reader R7
    verified is Eq. (1) exactly. It exists so that a caller of this package
    reads the profile from a cited place. It adds no physics.
    '''
    return get_c2n(np.asarray(h, dtype=float), wind_rms, cn2_ground)


def bufton_wind(h, slew_deg_s=0.0, ground_wind_m_s=10.0):
    '''
    Return the Bufton wind-speed profile V(h) [m/s].

    Parameters:
        h : float or numpy.ndarray
            Altitude above ground level [m].
        slew_deg_s : float
            Slew rate of the satellite as seen from the ground [deg/s].
        ground_wind_m_s : float
            Ground wind speed Vg [m/s].

    Returns:
        numpy.ndarray
            V(h) [m/s].

    formula:
        V(h) = v_s h + Vg + 30 exp(-((h - 9400)/4800)^2)
    Source: Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196, Ch. 12,
    Eq. (3), printed p. 481.

    This function DELEGATES to the shared kernel `v_wind`, which reader R7
    verified is Eq. (3) exactly. The kernel takes the slew rate in deg/s and
    converts it. It adds no physics.
    '''
    return v_wind(np.asarray(h, dtype=float), slew_deg_s, ground_wind_m_s)


def rms_wind(slew_deg_s=0.0, ground_wind_m_s=10.0, points=2001):
    '''
    Return the rms high-altitude wind speed w [m/s] of the Bufton profile.

    Parameters:
        slew_deg_s, ground_wind_m_s : as `bufton_wind`.
        points : int
            Number of grid points over the 5 km to 20 km band.

    Returns:
        float
            w [m/s]. Feed it to `hufnagel_valley`.

    formula:
        w = [ (1/15e3) INT_5e3^20e3 V^2(h) dh ]^(1/2)
    Source: Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196, Ch. 12,
    Eq. (2), printed p. 481.

    Note that w = 21 m/s of the H-V5/7 model is a CHOSEN value, not the value
    that Eq. (2) gives for the default Bufton parameters. The book selects
    w = 21 m/s and A = 1.7e-14 so that r0 = 5 cm and theta0 = 7 urad at
    lambda = 0.5 um (Ch. 12, text below Eq. (3), printed p. 481).
    '''
    hs = np.linspace(5.0e3, 20.0e3, int(points))
    v = bufton_wind(hs, slew_deg_s, ground_wind_m_s)
    return float(np.sqrt(np.trapz(v ** 2, hs) / 15.0e3))


# The two SCIDAR outer-scale models, as (peak L0 [m], peak altitude [m]). The
# denominator scale height is 2500 m in both. Source: Andrews and Phillips,
# 2nd ed. (2005), DOI 10.1117/3.626196, Ch. 12, Eq. (6) ("scidar4") and Eq. (7)
# ("scidar5"), printed p. 483.
OUTER_SCALE_MODELS = {
    'scidar4': (4.0, 8500.0),
    'scidar5': (5.0, 7500.0),
}


def outer_scale_profile(h, model='scidar4'):
    '''
    Return the outer scale L0(h) [m] of a SCIDAR altitude model.

    Parameters:
        h : float or numpy.ndarray
            Altitude above ground level [m].
        model : str
            A key of OUTER_SCALE_MODELS. "scidar4" is Eq. (6), which caps L0 at
            4 m. "scidar5" is Eq. (7), which caps L0 at 5 m.

    Returns:
        numpy.ndarray
            L0(h) [m].

    formula:
        L0(h) = A / [ 1 + ((h - h_c)/2500)^2 ]
        Eq. (6): A = 4 m, h_c = 8500 m
        Eq. (7): A = 5 m, h_c = 7500 m
    Source: Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196, Ch. 12,
    Sec. 12.2.2, Eqs. (6) and (7), printed p. 483.

    RESTRICTION stated by the book on the same page: Coulman et al. claim the
    two formulas hold above 2 km, but Tatarskii and Zavorotny say they hold in
    the surface layer only. The book leaves the question open. Near the ground
    the book quotes instead the rules L0 ~ 0.4 h or L0 ~ 0.5 h (Ch. 12, text
    above Eq. (6), printed p. 483), which this module does not build.

    Feed the result to `olb.turbulence.andrews.wander.beam_wander_variance_slant`
    as the array `L0`, which accepts a height-dependent outer scale.
    '''
    if model not in OUTER_SCALE_MODELS:
        raise ValueError(f'model must be one of {tuple(OUTER_SCALE_MODELS)}, '
                         f'not {model!r}')
    peak, h_c = OUTER_SCALE_MODELS[model]
    h = np.asarray(h, dtype=float)
    return peak / (1.0 + ((h - h_c) / 2500.0) ** 2)


def _xi(hs, altitude_m, direction, h0):
    '''
    Return the normalised path variable xi on the height grid.

    xi = 1 at the transmitter and xi = 0 at the receiver, for both directions.
    Source: Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196, Ch. 12,
    Eq. (14), printed p. 490.
    '''
    if direction not in _DIRECTIONS:
        raise ValueError(f'direction must be one of {_DIRECTIONS}, '
                         f'not {direction!r}')
    hs = np.asarray(hs, dtype=float)
    h0 = float(hs[0]) if h0 is None else float(h0)
    if altitude_m is None:
        raise ValueError('altitude_m (the satellite altitude H) is required')
    frac = (hs - h0) / (float(altitude_m) - h0)
    return np.clip(1.0 - frac if direction == 'uplink' else frac, 0.0, 1.0)


def mu(hs, cn2_profile, order, *, direction='uplink', beam=None,
       altitude_m=None, h0=None):
    '''
    Return one path moment mu_n of the Chapter 12 slant forms.

    Parameters:
        hs : numpy.ndarray
            Heights above the ground station [m], ascending.
        cn2_profile : numpy.ndarray
            ZENITH Cn2(h) profile on the hs grid [m^-2/3].
        order : int
            0, 1, 2 or 3. See the formulas below.
        direction : str
            "uplink" or "downlink". It selects xi through Ch. 12, Eq. (14).
            Orders 1 and 3 also read the beam.
        beam : BeamParams, optional
            The beam parameters in the RECEIVER plane, from
            `olb.turbulence.andrews.beam.beam_params`. Orders 1 and 3 need it.
        altitude_m : float, optional
            Satellite altitude H above the same datum as hs [m]. Orders 1, 2
            and 3 need it. Order 0 does not.
        h0 : float, optional
            Ground height h0 [m]. It defaults to hs[0].

    Returns:
        float
            mu_n [m^(1/3)].

    formula (xi from Ch. 12, Eq. (14); z/L = 1 - xi is the distance from the
    TRANSMITTER):
        mu_0 = INT Cn2(h) dh
        mu_1 = INT Cn2(h) |Theta + Theta_bar (1 - xi)|^(5/3) dh
        mu_2 = INT Cn2(h) xi^(5/3) dh
        mu_3 = Re INT Cn2(h) { xi^(5/6) [Lambda xi + i(1 - Theta_bar xi)]^(5/6)
                               - Lambda^(5/6) xi^(5/3) } dh
    Source: Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196:
        mu_0            Ch. 12, Eq. (21), printed p. 491; Eq. (85), printed
                        p. 522
        mu_1 downlink   Ch. 12, Eq. (18), printed p. 491; Eq. (86), printed
                        p. 522
        mu_1 uplink     Ch. 12, Eq. (25), printed p. 492; Eq. (94), printed
                        p. 523
        mu_2 downlink   Ch. 12, Eq. (19), printed p. 491; Eq. (87), printed
                        p. 522
        mu_2 uplink     Ch. 12, Eq. (26), printed p. 492; Eq. (95), printed
                        p. 523
        mu_3 downlink   Ch. 12, Eq. (37), printed p. 495; Eq. (88), printed
                        p. 522
        mu_3 uplink     Ch. 12, Eq. (55), printed p. 503; Eq. (96), printed
                        p. 523

    READING of mu_1. The book prints the mu_1 bracket with the plain height
    fraction (h - h0)/(H - h0) in BOTH Eq. (18) and Eq. (25). Read that way the
    two equations are identical, and the downlink coherence radius at the ground
    comes out near 900 m, which is absurd. This module uses |Theta + Theta_bar
    (1 - xi)|^(5/3) instead, that is the weight of the distance FROM THE
    TRANSMITTER. Three facts fix that reading:
        1. Ch. 6, Eq. (115), printed p. 209, gives the same moment on a general
           slant path as INT Cn2(z) |Theta + Theta_bar z/L|^(5/3) dz, with z
           measured from the transmitter, and Ch. 6, Eq. (116) confirms it in
           the spherical-wave limit (Theta = Lambda = 0).
        2. The book states below Eq. (19), printed p. 491, that mu_1d = mu_0 for
           a downlink from space. Only the (1 - xi) reading gives that.
        3. The book states below Eq. (27), printed p. 492, that the uplink
           coherence radius AT THE SATELLITE is many times larger than the
           satellite. Only the (1 - xi) reading gives that.
    The book's own Worked Example 2, printed p. 525, prints mu_1d = 1.98e-19,
    which is the literal-xi reading, and which contradicts facts 2 and 3. The
    self-check of this module measures both readings and prints both numbers.

    RESTRICTION. mu_3 uses the Kolmogorov spectrum only (Ch. 12, Eq. (15),
    printed p. 490).
    '''
    hs = np.asarray(hs, dtype=float)
    cn2 = np.asarray(cn2_profile, dtype=float)
    order = int(order)
    if order == 0:
        return float(np.trapz(cn2, hs))
    if order not in (1, 2, 3):
        raise ValueError(f'order must be 0, 1, 2 or 3, not {order}')

    xi = _xi(hs, altitude_m, direction, h0)
    if order == 2:
        return float(np.trapz(cn2 * xi ** (5.0 / 3.0), hs))

    if beam is None:
        raise ValueError(f'order {order} needs beam=BeamParams(...)')
    theta = float(beam.theta)
    theta_bar = float(beam.theta_bar)
    lam = float(beam.lam)

    if order == 1:
        weight = np.abs(theta + theta_bar * (1.0 - xi)) ** (5.0 / 3.0)
        return float(np.trapz(cn2 * weight, hs))

    inner = (lam * xi + 1j * (1.0 - theta_bar * xi)) ** (5.0 / 6.0)
    integrand = xi ** (5.0 / 6.0) * inner - lam ** (5.0 / 6.0) * xi ** (5.0 / 3.0)
    return float(np.trapz(cn2 * np.real(integrand), hs))


def _refuse_two_scale(l0, L0):
    '''Refuse an inner or outer scale, which Chapter 12 never gives.'''
    if l0 is not None or L0 is not None:
        raise NotImplementedError(
            'Chapter 12 of Andrews and Phillips, 2nd ed. (2005), '
            'DOI 10.1117/3.626196, uses the Kolmogorov spectrum only '
            '(Ch. 12, Eq. (15), printed p. 490). It gives no slant-path '
            'scintillation form with an inner scale or an outer scale. The '
            'coefficient is not guessed. Use '
            'olb.turbulence.andrews.scintillation.weak_two_scale_index for a '
            'single homogeneous path.')


def downlink_scintillation_index(hs, cn2_profile, wavelength, elevation_deg,
                                 *, D=None, regime='auto', l0=None, L0=None):
    '''
    Return the downlink scintillation index sigma_I^2 at the ground.

    Parameters:
        hs : numpy.ndarray
            Heights above the ground station [m], ascending. hs[0] is h0.
        cn2_profile : numpy.ndarray
            ZENITH Cn2(h) profile on the hs grid [m^-2/3].
        wavelength : float
            Optical wavelength [m].
        elevation_deg : float or numpy.ndarray
            Elevation angle above the horizon [deg].
        D : float, optional
            HARD receive-aperture diameter D_G [m]. None (or 0) gives the point
            receiver.
        regime : str
            "weak", "strong" or "auto". "auto" uses the book boundary
            sigma_R^2 < 1 (WEAK_REGIME_LIMIT).
        l0, L0 : float, optional
            REFUSED. See the module docstring.

    Returns:
        float or numpy.ndarray
            sigma_I^2.

    formula:
        point, weak
            sigma_R^2 = 2.25 k^(7/6) sec^(11/6)(zeta)
                        INT Cn2(h) (h - h0)^(5/6) dh
        hard aperture D_G, weak
            sigma_I^2(D_G) = 8.70 k^(7/6) sec^(11/6)(zeta)
                Re INT Cn2(h) { [a + i(h - h0)]^(5/6) - a^(5/6) } dh,
                a = k D_G^2 cos(zeta)/16
        point, weak to strong
            sigma_I^2 = exp[ 0.49 s^2/(1 + 1.11 s^(12/5))^(7/6)
                           + 0.51 s^2/(1 + 0.69 s^(12/5))^(5/6) ] - 1,
                s^2 = sigma_R^2
    Source: Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196:
        point weak      Ch. 12, Eq. (38), printed p. 495; repeated as Eq. (92),
                        printed p. 522
        hard aperture   Ch. 12, Eq. (39), printed p. 496
        weak to strong  Ch. 12, Eq. (40), printed p. 497; repeated as Eq. (93),
                        printed p. 522

    The aperture form is printed as
        8.70 k^(7/6) (H-h0)^(5/6) sec^(11/6)(zeta)
        Re INT Cn2(h) { [k D_G^2/(16 L) + i xi]^(5/6)
                        - (k D_G^2/(16 L))^(5/6) } dh,
    with xi = (h - h0)/(H - h0) and L = (H - h0) sec(zeta). The code pulls
    (H - h0)^(5/6) inside the bracket, which removes H and L from the
    expression, because (H - h0) k D_G^2/(16 L) = k D_G^2 cos(zeta)/16. So this
    function needs no satellite altitude.

    The book states below Eq. (39), printed p. 496, that Eq. (39) reduces to
    Eq. (38) at D_G = 0. It does so only to the rounding of the two constants:
    8.70 cos(5 pi/12) = 2.2517, and the book prints 2.25. The self-check
    measures that 0.08 % gap.

    RESTRICTION. The book uses the PLANE-wave reduction for a downlink, because
    a wave from space enters the atmosphere at 20 km already close to a plane
    wave (Ch. 12, text below Eq. (21), printed p. 491). The full Gaussian-beam
    downlink index, Ch. 12, Eqs. (36) and (37), printed p. 495, is not built
    here. Compose `mu(..., order=3, direction="downlink")` with the beam if you
    need it.

    REGIME. The strong branch has no aperture form in the book, so a strong
    regime with a finite D_G is refused.
    '''
    _refuse_two_scale(l0, L0)
    hs = np.asarray(hs, dtype=float)
    cn2 = np.asarray(cn2_profile, dtype=float)
    k = wavenumber(wavelength)
    sec = sec_zeta(elevation_deg)
    dh = hs - hs[0]

    point = (_POINT_CONSTANT * k ** (7.0 / 6.0) * sec ** (11.0 / 6.0)
             * np.trapz(cn2 * dh ** (5.0 / 6.0), hs))

    if regime not in ('weak', 'strong', 'auto'):
        raise ValueError(f'regime must be "weak", "strong" or "auto", '
                         f'not {regime!r}')
    strong = regime == 'strong' or (regime == 'auto'
                                    and np.all(point >= WEAK_REGIME_LIMIT))

    if strong:
        if D:
            raise NotImplementedError(
                'Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196, '
                'give no aperture-averaged downlink index in the strong '
                'regime. Ch. 12, Eq. (39), printed p. 496, is a weak-theory '
                'form, and Ch. 12, Eq. (40), printed p. 497, is a point form. '
                'Use olb.turbulence.andrews.aperture.averaged_index for a '
                'single homogeneous path.')
        return (np.exp(large_scale_log_variance(point)
                       + small_scale_log_variance(point)) - 1.0)

    if not D:
        return point

    # Eq. (39). The bracket carries the length a = k D_G^2 cos(zeta)/16 [m].
    a = k * float(D) ** 2 / (16.0 * sec)
    a_arr = np.atleast_1d(np.asarray(a, dtype=float))
    kernel = ((a_arr[None, :] + 1j * dh[:, None]) ** (5.0 / 6.0)
              - a_arr[None, :] ** (5.0 / 6.0))
    integral = np.trapz(cn2[:, None] * np.real(kernel), hs, axis=0)
    out = (_MU3_CONSTANT * k ** (7.0 / 6.0)
           * np.atleast_1d(sec) ** (11.0 / 6.0) * integral)
    return float(out[0]) if np.ndim(sec) == 0 else out.reshape(np.shape(sec))


def _uplink_longitudinal(sigma2_bu, theta, strong):
    '''
    Return the uplink longitudinal (on-axis Rytov) scintillation index.

    formula:
        weak    sigma_B_u^2
        strong  exp[ 0.49 b^2/(1 + 0.56 (1 + Theta) b^(12/5))^(7/6)
                   + 0.51 b^2/(1 + 0.69 b^(12/5))^(5/6) ] - 1,  b^2 = sigma_B_u^2
    Source: Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196, Ch. 12,
    Eq. (58), printed p. 504 (weak), and Ch. 12, Eq. (60), printed p. 506
    (weak to strong), repeated as Eq. (99), printed p. 524.

    The strong branch repeats the algebra of
    `olb.turbulence.andrews.scintillation.large_scale_log_variance` and
    `small_scale_log_variance` with wave="gaussian". It cannot call them,
    because they rebuild sigma_B^2 from the HOMOGENEOUS closed form Ch. 8,
    Eq. (23), while a slant path takes sigma_B_u^2 from the mu_3u integral.
    '''
    if not strong:
        return sigma2_bu
    b125 = sigma2_bu ** (6.0 / 5.0)
    x = 0.49 * sigma2_bu / (1.0 + 0.56 * (1.0 + theta) * b125) ** (7.0 / 6.0)
    y = 0.51 * sigma2_bu / (1.0 + 0.69 * b125) ** (5.0 / 6.0)
    return np.exp(x + y) - 1.0


def uplink_scintillation_index(hs, cn2_profile, wavelength, elevation_deg,
                               beam, *, altitude_m, r=0.0, tracked=True,
                               regime='auto', pointing_error_m=None,
                               wander_rms_m=None, r0=None, c_r=None,
                               l0=None, L0=None):
    '''
    Return the uplink scintillation index sigma_I^2 at the satellite.

    Parameters:
        hs : numpy.ndarray
            Heights above the ground station [m], ascending. hs[0] is h0.
        cn2_profile : numpy.ndarray
            ZENITH Cn2(h) profile on the hs grid [m^-2/3].
        wavelength : float
            Optical wavelength [m].
        elevation_deg : float
            Elevation angle above the horizon [deg].
        beam : BeamParams
            The beam parameters at the SATELLITE, from
            `olb.turbulence.andrews.beam.beam_params(W0, wavelength, L, f0)`
            with L the slant range.
        altitude_m : float
            Satellite altitude H above the same datum as hs [m].
        r : float or numpy.ndarray
            Off-axis radius in the satellite plane [m].
        tracked : bool
            True removes the beam wander (Ch. 12, Eqs. (57) and (59)). False
            keeps it (Ch. 12, Eqs. (54), (56) and (61)).
        regime : str
            "weak", "strong" or "auto". "auto" uses the boundary
            sigma_B_u^2 < 1 (WEAK_REGIME_LIMIT).
        pointing_error_m : float, optional
            The rms wander-induced pointing error sigma_pe [m]. It acts on the
            UNTRACKED model. None composes
            `olb.turbulence.andrews.wander.pointing_error_variance_slant`.
        wander_rms_m : float, optional
            The rms beam-wander displacement sqrt(<r_c^2>) [m]. It acts on the
            TRACKED model. None composes
            `olb.turbulence.andrews.wander.beam_wander_variance_slant`.
        r0 : float, optional
            Fried's parameter [m]. None composes
            `olb.turbulence.andrews.wander.plane_fried_parameter_slant`.
        c_r : float, optional
            Scaling constant of the pointing-error cut-off kappa_r = c_r/r0.
            None takes the wander module default. The book leaves it free
            (Ch. 12, text below Eq. (53), printed p. 503) and uses
            kappa_r = 3.86/r0 in Fig. 12.13 and pi/r0 in Figs. 12.14 to 12.17.
        l0, L0 : float, optional
            REFUSED. See the module docstring.

    Returns:
        float or numpy.ndarray
            sigma_I^2 at the radius r.

    formula:
        sigma_B_u^2 = 8.70 mu_3u k^(7/6) (H - h0)^(5/6) sec^(11/6)(zeta)
        C           = 5.95 (2 W0/r0)^(5/3) / W^2
        tracked     sigma_I^2 = C (r - sqrt(<r_c^2>))^2 U(r - sqrt(<r_c^2>))
                                + long(sigma_B_u^2)
        untracked   sigma_I^2 = C [ (r + sigma_pe)^2 U(r - sigma_pe)
                                    + sigma_pe^2 ] + long(sigma_B_u^2)
    with long(.) the identity in the weak regime and the Ch. 12, Eq. (60)
    transform in the strong regime. W is the free-space beam radius at the
    satellite and W0 is the transmitter beam radius, both from `beam`.
    Source: Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196:
        mu_3u             Ch. 12, Eq. (55), printed p. 503
        sigma_B_u^2       Ch. 12, Eq. (58), printed p. 504
        untracked, weak   Ch. 12, Eq. (54), printed p. 503 (on axis), and
                          Ch. 12, Eq. (56), printed p. 504 (off axis)
        tracked, weak     Ch. 12, Eq. (57), printed p. 504
        tracked, strong   Ch. 12, Eqs. (59) and (60), printed p. 506
        untracked, strong Ch. 12, Eq. (61), printed p. 506
    The book writes the radial factor with angles, (alpha_r/W)^2 with
    alpha_r = r/L, and puts (H - h0)^2 sec^2(zeta) = L^2 in front. The code
    cancels the two, which gives (r/W)^2.

    THIS DOES NOT MODEL A PRE-COMPENSATED UPLINK (OLB GAP 2, decision
    2026-08-27). `tracked=True` models a beam with the wander FULLY removed:
    a perfect tilt correction. A real beacon measurement decorrelates from
    the uplink path over the point-ahead angle, so a residual tilt survives,
    and this form charges nothing for it. The same argument applies to each
    higher corrected order. Also, a decorrelated higher-order correction
    reshapes the beam at the satellite, and these forms normalise by the
    vacuum-diffraction beam radius W. So `tracked=True` is OPTIMISTIC for a
    pre-compensated uplink, and it is not a bound in either direction. An
    earlier docstring called it the floor of the residual scintillation; that
    claim was wrong. Decision: no analytic Term models the pre-compensated
    scintillation. The model of record is the fidelity-1 FAST Monte Carlo
    with the point-ahead offset (olb/models/fast.py, DTHETA;
    backlog item 1-2). The tracked form stays valid for what it names: a
    tilt-tracked, otherwise uncorrected beam.

    RESTRICTION. The beam-wander term stays a weak-fluctuation result in every
    branch, which the book states below Eq. (59), printed p. 506. The book also
    limits the weak forms to a transmitter radius near 20 cm or less (Ch. 12,
    text below Eq. (56), printed p. 504).
    '''
    _refuse_two_scale(l0, L0)
    hs = np.asarray(hs, dtype=float)
    cn2 = np.asarray(cn2_profile, dtype=float)
    k = wavenumber(wavelength)
    sec = float(sec_zeta(elevation_deg))
    h0 = float(hs[0])
    height = float(altitude_m) - h0
    range_m = height * sec

    w = float(beam.w)
    w0 = w / np.sqrt(float(beam.theta0) ** 2 + float(beam.lambda0) ** 2)
    theta0 = float(beam.theta0)
    f0 = np.inf if theta0 == 1.0 else range_m / (1.0 - theta0)

    mu3u = mu(hs, cn2, 3, direction='uplink', beam=beam, altitude_m=altitude_m)
    sigma2_bu = (_MU3_CONSTANT * mu3u * k ** (7.0 / 6.0)
                 * height ** (5.0 / 6.0) * sec ** (11.0 / 6.0))

    if regime not in ('weak', 'strong', 'auto'):
        raise ValueError(f'regime must be "weak", "strong" or "auto", '
                         f'not {regime!r}')
    strong = regime == 'strong' or (regime == 'auto'
                                    and sigma2_bu >= WEAK_REGIME_LIMIT)
    longitudinal = _uplink_longitudinal(sigma2_bu, float(beam.theta), strong)

    if r0 is None:
        r0 = plane_fried_parameter_slant(wavelength, hs, cn2, elevation_deg)
    coef = (_WANDER_SCINT_CONSTANT * (2.0 * w0 / float(r0)) ** (5.0 / 3.0)
            / w ** 2)
    r = np.asarray(r, dtype=float)

    if tracked:
        if wander_rms_m is None:
            wander_rms_m = np.sqrt(beam_wander_variance_slant(
                w0, wavelength, hs, cn2, range_m, f0=f0,
                elevation_deg=elevation_deg))
        rc = float(wander_rms_m)
        radial = coef * np.where(r > rc, (r - rc) ** 2, 0.0)
        return radial + longitudinal

    if pointing_error_m is None:
        kwargs = {} if c_r is None else {'c_r': c_r}
        pointing_error_m = np.sqrt(pointing_error_variance_slant(
            w0, wavelength, hs, cn2, range_m, f0=f0,
            elevation_deg=elevation_deg, r0=r0, **kwargs))
    pe = float(pointing_error_m)
    radial = coef * (np.where(r > pe, (r + pe) ** 2, 0.0) + pe ** 2)
    return radial + longitudinal


def uplink_coherence_radius(hs, cn2_profile, wavelength, elevation_deg, beam,
                            *, altitude_m):
    '''
    Return the uplink spatial coherence radius rho_0 [m] AT THE SATELLITE.

    Parameters:
        hs, cn2_profile, wavelength, elevation_deg : as
            `uplink_scintillation_index`.
        beam : BeamParams
            The beam parameters at the satellite.
        altitude_m : float
            Satellite altitude H [m].

    Returns:
        float
            rho_0 [m].

    formula:
        D(rho, L) = 2.91 k^2 rho^(5/3) sec(zeta)
                    (mu_1u + 0.62 mu_2u Lambda^(11/6))
        rho_0     = [ 1.46 k^2 sec(zeta)
                      (mu_1u + 0.62 mu_2u Lambda^(11/6)) ]^(-3/5)
    Source: Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196, Ch. 12,
    Eqs. (24) to (27), printed p. 492. The constant 1.46 is 2.91/2 rounded,
    which follows from the definition D(rho_0) = 2 (Ch. 12, text above Eq. (22),
    printed p. 491). Eq. (27) is derived from Eq. (24), so the weight of the
    Lambda^(11/6) term is the 0.62 that Eq. (24) prints.

    THIS IS MATRIX ROW G-130. It is the book's OWN uplink form. It is NOT the
    mirror of the downlink pair Ch. 12, Eqs. (17) to (22). Read the plane of
    reference before you compare it with anything:
        This function returns the coherence radius of the uplink wave IN THE
        SATELLITE PLANE. Through mu_1u it weights the turbulence by the distance
        from the GROUND TRANSMITTER, so it is small only for a beam that is
        nearly a plane wave. For a satellite link it is hundreds of metres, and
        the book says so below Eq. (27), printed p. 492: "the spatial coherence
        radius at the satellite will be many times larger than the probable size
        of the satellite".
        The kernel `spherical_wave_coherence_diameter` in
        my_analysis_modules/coupled_flux.py returns a DIFFERENT quantity: the
        GROUND-referred Fried parameter, weighted by ((L - z)/L)^(5/3). On a
        satellite uplink that weight is 1 over the whole turbulent layer, so
        that kernel reduces to Andrews Ch. 12, Eq. (23), printed p. 492, which
        is the same r0 that the book itself feeds into the uplink beam-wander
        and pointing-error equations (50), (51) and (53). So the two are not in
        conflict. The self-check measures both ratios.
    See Conflict C-02 in docs/andrews-crosscheck.md.

    RESTRICTION. Kolmogorov spectrum only (Ch. 12, Eq. (15), printed p. 490).
    The book notes below Eq. (27) that non-Kolmogorov stratospheric turbulence
    changes this result (Gurvich and Belen'kii).
    '''
    k = wavenumber(wavelength)
    sec = sec_zeta(elevation_deg)
    mu1u = mu(hs, cn2_profile, 1, direction='uplink', beam=beam,
              altitude_m=altitude_m)
    mu2u = mu(hs, cn2_profile, 2, direction='uplink', altitude_m=altitude_m)
    weighted = mu1u + _LAMBDA_TERM * mu2u * float(beam.lam) ** (11.0 / 6.0)
    return float((_RHO0_CONSTANT * k ** 2 * sec * weighted) ** (-3.0 / 5.0))


def isoplanatic_angle(hs, cn2_profile, wavelength, elevation_deg=90.0, *,
                      beam=None, altitude_m=None):
    '''
    Return the isoplanatic angle theta_0 [rad] of an upward path.

    Parameters:
        hs, cn2_profile, wavelength, elevation_deg : as
            `uplink_scintillation_index`.
        beam : BeamParams, optional
            The beam parameters at the satellite. None takes the spherical-wave
            limit Theta = Lambda = 0.
        altitude_m : float, optional
            Satellite altitude H [m]. The Gaussian-beam branch needs it.

    Returns:
        float
            theta_0 [rad].

    formula:
        Gaussian beam  theta_0 = cos^(8/5)(zeta)
                                 / [ (H - h0)
                                     (2.91 k^2 (mu_1u
                                      + 0.62 mu_2u Lambda^(11/6)))^(3/5) ]
        spherical wave theta_0 = cos^(8/5)(zeta)
                                 [ 2.91 k^2 INT Cn2(h) (h - h0)^(5/3) dh ]^(-3/5)
    Source: Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196, Ch. 12,
    Eq. (29) (Gaussian beam) and Eq. (30) (spherical wave), printed p. 493.
    The book states below Eq. (29) that Eq. (30) is the special case
    Theta = Lambda = 0, and the algebra confirms it: mu_1u then becomes
    (H - h0)^(-5/3) INT Cn2 (h - h0)^(5/3) dh, and the (H - h0) prefactor
    cancels.

    olb already holds a second route, `olb.turbulence.anisoplanatism.
    isoplanatic_angle`, which uses the Stone et al. 1994 constant 2.914381
    (DOI 10.1364/JOSAA.11.000347) in place of the book's rounded 2.91. The two
    differ by 0.09 %, because (2.914381/2.91)^(-3/5) = 0.99910. The self-check
    measures the difference.

    RESTRICTION. The book notes in the footnote on printed p. 493 that one
    experimental study says the isoplanatic angle does not apply to a
    tilt-related quantity, so it may not describe a point-ahead tracking error.
    '''
    hs = np.asarray(hs, dtype=float)
    cn2 = np.asarray(cn2_profile, dtype=float)
    k = wavenumber(wavelength)
    cos_zeta = 1.0 / sec_zeta(elevation_deg)

    if beam is None:
        integral = np.trapz(cn2 * (hs - hs[0]) ** (5.0 / 3.0), hs)
        return float(cos_zeta ** (8.0 / 5.0)
                     * (_ISOPLANATIC_CONSTANT * k ** 2 * integral)
                     ** (-3.0 / 5.0))

    height = float(altitude_m) - float(hs[0])
    mu1u = mu(hs, cn2, 1, direction='uplink', beam=beam,
              altitude_m=altitude_m)
    mu2u = mu(hs, cn2, 2, direction='uplink', altitude_m=altitude_m)
    weighted = mu1u + _LAMBDA_TERM * mu2u * float(beam.lam) ** (11.0 / 6.0)
    return float(cos_zeta ** (8.0 / 5.0)
                 / (height * (_ISOPLANATIC_CONSTANT * k ** 2 * weighted)
                    ** (3.0 / 5.0)))


def point_ahead_angle(velocity_m_s):
    '''
    Return the point-ahead angle theta_p [rad] of a moving satellite.

    Parameters:
        velocity_m_s : float or numpy.ndarray
            Satellite speed perpendicular to the line of sight [m/s].

    Returns:
        float or numpy.ndarray
            theta_p [rad].

    formula:
        theta_p = 2 V / c
    Source: Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196, Ch. 12,
    Sec. 12.3.3, text at Fig. 12.4, printed p. 488, and Ch. 12, Sec. 12.4.3,
    printed p. 493.

    The book states on printed p. 493 that theta_p is USUALLY MUCH LARGER than
    the isoplanatic angle, so a wave-front measured along the tracking path
    does not correct the turbulence along the transmit path. A LEO satellite at
    7 km/s gives 47 urad.
    '''
    return 2.0 * np.asarray(velocity_m_s, dtype=float) / _C


if __name__ == '__main__':
    from .beam import beam_params

    # ================= part 1: physics self-checks =================
    # The book's own numbers come from Ch. 12, Worked Examples 1 and 2, printed
    # pp. 524 and 525: a GEO link at H = 38.5e6 m, lambda = 1.06 um, W0 = 2 cm,
    # collimated, zenith angle 30 deg, H-V5/7.
    H_GEO = 38.5e6
    LAM = 1.06e-6
    W0 = 0.02
    ELEV = 60.0                          # zenith angle 30 deg
    hs = np.linspace(0.0, 20.0e3, 200001)
    cn2 = hufnagel_valley(hs)
    L_geo = H_GEO * sec_zeta(ELEV)
    bp = beam_params(W0, LAM, L_geo)

    assert abs(sec_zeta(90.0) - 1.0) < 1e-15
    assert abs(sec_zeta(30.0) - 2.0) < 1e-12
    print(f'sec_zeta 30 deg elevation : {sec_zeta(30.0):.6f}  (exact 2)')

    m0 = mu(hs, cn2, 0)
    assert abs(m0 / 2.24e-12 - 1.0) < 0.01, m0
    print(f'BOOK mu_0 H-V5/7          : {m0:.4e} m^(1/3)  (book 2.24e-12)')

    m3u = mu(hs, cn2, 3, direction='uplink', beam=bp, altitude_m=H_GEO)
    assert abs(m3u / 3.70e-17 - 1.0) < 0.03, m3u
    print(f'BOOK mu_3u Example 1      : {m3u:.4e} m^(1/3)  (book 3.70e-17)')

    assert abs(bp.w - 750.0) < 1.0, bp.w
    print(f'BOOK W at GEO             : {bp.w:.1f} m  (book 750 m)')

    r0_geo = plane_fried_parameter_slant(LAM, hs, cn2, ELEV)
    assert abs(r0_geo / 0.1124 - 1.0) < 0.01, r0_geo
    print(f'BOOK r0 Example 1         : {r0_geo * 100:.2f} cm  (book 11.24 cm)')

    s2_track = uplink_scintillation_index(hs, cn2, LAM, ELEV, bp,
                                          altitude_m=H_GEO, tracked=True)
    assert abs(s2_track - 0.07) < 0.005, s2_track
    print(f'BOOK uplink tracked       : {s2_track:.4f}  (book 0.07)')

    s2_untrack = uplink_scintillation_index(hs, cn2, LAM, ELEV, bp,
                                            altitude_m=H_GEO, tracked=False)
    assert abs(s2_untrack - 0.095) < 0.006, s2_untrack
    print(f'BOOK uplink untracked     : {s2_untrack:.4f}  (book 0.095)')
    assert s2_untrack > s2_track

    s2_down = downlink_scintillation_index(hs, cn2, LAM, ELEV)
    assert abs(s2_down - 0.13) < 0.006, s2_down
    print(f'BOOK downlink on axis     : {s2_down:.4f}  (book 0.13)')

    th0 = isoplanatic_angle(hs, cn2, LAM, ELEV)
    assert abs(th0 * 1e6 / 13.5 - 1.0) < 0.03, th0
    print(f'BOOK isoplanatic angle    : {th0 * 1e6:.2f} urad  (book 13.5)')

    # Ch. 12, Eq. (29) must reduce to Eq. (30) in the spherical-wave limit.
    bp_sph = beam_params(1e-6, LAM, L_geo)
    th0_beam = isoplanatic_angle(hs, cn2, LAM, ELEV, beam=bp_sph,
                                 altitude_m=H_GEO)
    err_iso = abs(th0_beam / th0 - 1.0)
    assert err_iso < 1e-9, err_iso
    print(f'LIMIT Eq. (29) -> Eq. (30): |ratio - 1| = {err_iso:.3e} '
          f'(target 1e-9)')

    # Ch. 12, Eq. (39) must reduce to Eq. (38) at D = 0. The gap is the book
    # rounding of 8.70 cos(5 pi/12) = 2.2517 to 2.25.
    s2_ap0 = downlink_scintillation_index(hs, cn2, LAM, ELEV, D=1e-12)
    gap39 = s2_ap0 / s2_down - 1.0
    assert abs(gap39) < 1e-3, gap39
    print(f'LIMIT Eq. (39) D -> 0     : {gap39 * 100:+.4f} % against Eq. (38) '
          f'(book rounding 2.2517 -> 2.25)')

    # The uplink coherence radius at the satellite must be huge, which the book
    # states below Ch. 12, Eq. (27), printed p. 492.
    rho_up = uplink_coherence_radius(hs, cn2, LAM, ELEV, bp, altitude_m=H_GEO)
    assert rho_up > 100.0, rho_up
    print(f'BOOK uplink rho_0 at sat  : {rho_up:.1f} m  (book "many times '
          f'larger than the satellite")')

    # Outer-scale models. Ch. 12, Eqs. (6) and (7), printed p. 483.
    assert abs(outer_scale_profile(8500.0, 'scidar4') - 4.0) < 1e-12
    assert abs(outer_scale_profile(7500.0, 'scidar5') - 5.0) < 1e-12
    assert outer_scale_profile(0.0, 'scidar4') < 4.0
    print(f'Eq. (6) L0 peak / at 0 m  : '
          f'{outer_scale_profile(8500.0, "scidar4"):.3f} m / '
          f'{outer_scale_profile(0.0, "scidar4"):.3f} m')

    # Bufton wind and its rms. Ch. 12, Eqs. (2) and (3), printed p. 481.
    w_rms = rms_wind()
    assert 20.0 < w_rms < 40.0, w_rms
    print(f'Eq. (2) rms wind Vg=10    : {w_rms:.2f} m/s  (the H-V5/7 w = 21 '
          f'm/s is a CHOSEN value, not this one)')
    assert abs(bufton_wind(9400.0, 0.0, 10.0) - 40.0) < 1e-9

    # Point-ahead angle. Ch. 12, Sec. 12.3.3, printed p. 488.
    pa = point_ahead_angle(7.0e3)
    assert abs(pa * 1e6 - 46.7) < 0.5, pa
    print(f'LEO point-ahead 7 km/s    : {pa * 1e6:.1f} urad  (book "50 urad")')
    assert pa > th0, 'the book states the point-ahead angle exceeds theta_0'

    # The refusals. Chapter 12 has no two-scale slant form.
    for call in (lambda: downlink_scintillation_index(hs, cn2, LAM, ELEV,
                                                      l0=0.005),
                 lambda: uplink_scintillation_index(hs, cn2, LAM, ELEV, bp,
                                                    altitude_m=H_GEO, L0=10.0)):
        try:
            call()
        except NotImplementedError:
            pass
        else:
            raise AssertionError('a two-scale slant form must be refused')
    print('REFUSED inner/outer scale : Ch. 12 uses Kolmogorov only, Eq. (15)')

    # ================= part 2: REDUCTION checks =================
    from .. import anisoplanatism, beam_wave_scintillation
    from .. import plane_wave_scintillation as pws
    from ..._deps import get_c2n as _kernel_c2n
    from ..._deps import spherical_wave_coherence_diameter

    # 1. hufnagel_valley is the kernel, with no change at all.
    err_hv = float(np.max(np.abs(hufnagel_valley(hs, 21.0, 1.7e-14)
                                 - _kernel_c2n(hs, 21.0, 1.7e-14))))
    assert err_hv == 0.0, err_hv
    print(f'REDUCTION hufnagel_valley vs _deps.get_c2n : max |diff| = '
          f'{err_hv:.1e}  (target exact)')

    # 2. The weak point downlink index is the old plane-wave slant integral.
    old_point = pws.plane_wave_scintillation_index(ELEV, LAM, hs, cn2)
    err_point = abs(s2_down / old_point - 1.0)
    assert err_point < 1e-9, err_point
    print(f'REDUCTION downlink point vs plane_wave_scintillation_index : '
          f'|ratio - 1| = {err_point:.3e}  (target 1e-9)')

    # 3. The book hard-aperture form Ch. 12, Eq. (39) against the olb numerical
    #    integral with the hard Airy filter. See Conflict C-06: the two are not
    #    the same filter, so this is a MEASUREMENT, not a test.
    print('REDUCTION Eq. (39) vs aperture_averaged_scintillation_index :')
    print('    D [m]     Eq. (39)      olb Airy      ratio')
    for d_m in (0.05, 0.20, 1.00):
        book = downlink_scintillation_index(hs, cn2, LAM, ELEV, D=d_m)
        olb_num = pws.aperture_averaged_scintillation_index(d_m, ELEV, LAM, hs,
                                                            cn2)
        print(f'    {d_m:5.2f}   {book:.6e}  {olb_num:.6e}  '
              f'{book / olb_num:8.4f}')

    # 4. G-130. The book uplink coherence radius against the kernel.
    hs_up = np.linspace(0.0, 20.0e3, 4001)
    cn2_up = hufnagel_valley(hs_up)
    k_up = wavenumber(LAM)
    z_path = hs_up * sec_zeta(ELEV)
    kernel_r0s = spherical_wave_coherence_diameter(k_up, L_geo, cn2_up, z_path)
    rho_book = uplink_coherence_radius(hs_up, cn2_up, LAM, ELEV, bp,
                                       altitude_m=H_GEO)
    fried_23 = plane_fried_parameter_slant(LAM, hs_up, cn2_up, ELEV)
    rho_down = (1.45 * mu(hs_up, cn2_up, 0) * k_up ** 2
                * sec_zeta(ELEV)) ** (-3.0 / 5.0)
    print(f'REDUCTION G-130 uplink coherence radius :')
    print(f'    book Ch. 12, Eq. (27) rho_0 at the satellite : '
          f'{rho_book:.3f} m')
    print(f'    kernel spherical_wave_coherence_diameter     : '
          f'{kernel_r0s:.5f} m')
    print(f'    ratio kernel / Eq. (27)                      : '
          f'{kernel_r0s / rho_book:.3e}  (DIFFERENT quantities)')
    print(f'    book Ch. 12, Eq. (23) ground Fried r0        : '
          f'{fried_23:.5f} m')
    print(f'    ratio kernel / Eq. (23)                      : '
          f'{kernel_r0s / fried_23:.6f}  (target 1.000)')
    print(f'    ratio kernel / (2.1 x downlink Eq. (22))     : '
          f'{kernel_r0s / (2.1 * rho_down):.6f}')
    assert abs(kernel_r0s / fried_23 - 1.0) < 0.01, kernel_r0s / fried_23

    # 5. A narrow layer must give the homogeneous closed form of mu_2.
    #    mu_2 = INT Cn2 xi^(5/3) dh -> Cn2 dh xi(h_layer)^(5/3).
    h_lo, h_hi = 9.0e3, 9.2e3
    hs_layer = np.linspace(h_lo, h_hi, 4001)
    cn2_layer = np.full_like(hs_layer, 1e-16)
    mu2_layer = mu(hs_layer, cn2_layer, 2, direction='uplink',
                   altitude_m=H_GEO, h0=0.0)
    xi_mid = 1.0 - 0.5 * (h_lo + h_hi) / H_GEO
    closed = 1e-16 * (h_hi - h_lo) * xi_mid ** (5.0 / 3.0)
    err_mu = abs(mu2_layer / closed - 1.0)
    assert err_mu < 0.02, err_mu
    print(f'REDUCTION narrow-layer mu_2 vs homogeneous closed form : '
          f'|ratio - 1| = {err_mu:.3e}  (target 2e-2)')

    # 6. The slant uplink longitudinal index against the Dios path integral.
    #    This is the gap-9 twin on a real slant path. The TRACKED branch is the
    #    true twin, because the Dios route carries no beam wander.
    hs_d = np.linspace(0.0, 20.0e3, 4001)
    cn2_d = hufnagel_valley(hs_d)
    W0_D = 0.10
    L_leo = 600.0e3 * sec_zeta(ELEV)
    bp_d = beam_params(W0_D, LAM, L_leo)
    dios = beam_wave_scintillation.on_axis_scintillation_index(
        hs_d, cn2_d, W0_D, LAM, elevation_deg=ELEV, path_length_m=L_leo)
    andrews_tr = uplink_scintillation_index(hs_d, cn2_d, LAM, ELEV, bp_d,
                                            altitude_m=600.0e3, tracked=True,
                                            regime='weak')
    andrews_un = uplink_scintillation_index(hs_d, cn2_d, LAM, ELEV, bp_d,
                                            altitude_m=600.0e3, tracked=False,
                                            regime='weak')
    print(f'REDUCTION slant uplink index vs the Dios route '
          f'(600 km, W0 = 10 cm, 60 deg elevation) :')
    print(f'    Dios on_axis_scintillation_index      : {dios:.6e}')
    print(f'    Andrews Ch. 12, Eq. (58) tracked      : {andrews_tr:.6e}  '
          f'({(andrews_tr / dios - 1.0) * 100:+.2f} %)')
    print(f'    Andrews Ch. 12, Eq. (54) untracked    : {andrews_un:.6e}  '
          f'({(andrews_un / dios - 1.0) * 100:+.2f} %  the wander term)')

    # 7. The Andrews isoplanatic angle against the Stone route already in olb.
    th_stone = anisoplanatism.isoplanatic_angle(hs, cn2, LAM, ELEV)
    err_stone = abs(th0 / th_stone - 1.0)
    assert err_stone < 0.01, err_stone
    print(f'REDUCTION isoplanatic_angle vs anisoplanatism.isoplanatic_angle : '
          f'|ratio - 1| = {err_stone:.3e}  (target 1e-2; constants 2.91 '
          f'against 2.914381)')

    # 8. The two readings of mu_1. See the mu docstring.
    xi_d = _xi(hs, H_GEO, 'downlink', None)
    literal = float(np.trapz(cn2 * np.abs(bp.theta + bp.theta_bar * xi_d)
                             ** (5.0 / 3.0), hs))
    used = mu(hs, cn2, 1, direction='downlink', beam=bp, altitude_m=H_GEO)
    print(f'MEASURED mu_1d two readings : (1 - xi) reading = {used:.4e} '
          f'(= mu_0 {m0:.4e}, as the book text needs); literal-xi reading = '
          f'{literal:.4e} (= the book Worked Example 2 value 1.98e-19)')

    print('self-check passed')
