'''
Pure plane-wave scintillation physics for a LEO-to-ground downlink.

The satellite is far away. The downlink source is a plane wave at the top of the
atmosphere. This module gives the plane-wave scintillation index and the
aperture-averaging integral for the plane wave that propagates down through
turbulence to the ground aperture. The functions are pure. They take numeric
arrays and return numeric arrays. The Term factory lives in olb.links.downlink.

Physics (plane wave, weak fluctuation, isotropic turbulence):
    The plane-wave scintillation index (point receiver), integrated over the
    Cn2 slant path, is
        sigma2_I = 2.25 * k^(7/6) * sec(zeta)^(11/6) * integral[ Cn2(h) h^(5/6) dh ]
    Source: Andrews and Phillips, Laser Beam Propagation through Random Media,
    2nd ed. (2005), Ch. 12, Eq. (38), printed p. 495 (repeated as Ch. 12,
    Eq. (92), printed p. 522). DOI 10.1117/3.626196.
    Here k = 2*pi/lambda, zeta is the zenith angle,
    sec(zeta) = 1/sin(elevation_deg), and h is the height above the ground
    station. The constant 2.25 and the h^(5/6) weight are from that equation.

    A ground telescope of diameter D averages the scintillation. This module
    computes the aperture-averaged flux scintillation index sigma2_P from the
    distributed-path Rytov double integral over height h and spatial wavenumber
    kappa:
        sigma2_P = 8 * pi^2 * k^2 * sec(zeta)
            * INT_h Cn2(h) * [ INT_kappa kappa * 0.033 * kappa^(-11/3)
            * (1 - cos(kappa^2 * h * sec(zeta) / k))
            * (2*J1(kappa*D/2) / (kappa*D/2))^2 dkappa ] dh
    Source: Andrews and Phillips, 2nd ed. (2005), plane-wave scintillation index
    Ch. 12, Eq. (38) and aperture-averaging filter Ch. 10. The filter
    (2*J1(x)/x)^2 assumes a uniform circular aperture with no central
    obscuration. The aperture-averaging factor is A = sigma2_P / sigma2_I. It
    obeys 0 < A <= 1, A -> 1 as D -> 0, and the large-aperture asymptote
    A ~ D^(-7/3).
'''

import numpy as np
from scipy.special import j1

from .andrews.scintillation import (
    scintillation_index as _andrews_scintillation_index,
)
from .profiles import DEFAULT_HS, get_c2n

# Weak-fluctuation limit for the plane-wave scintillation index. The lognormal
# irradiance PDF is trusted for sigma2_I below 0.25.
#
# HOUSE RULE, NOT A BOOK NUMBER. Andrews and Phillips, Laser Beam Propagation
# through Random Media, 2nd ed. (2005), DOI 10.1117/3.626196, give the weak
# limit as sigma_R^2 < 1 (Ch. 5, Eq. (15) and the text after it, printed
# p. 140; also Ch. 10, Eq. (61), printed p. 412, and Ch. 12, Eq. (40), printed
# p. 497). The value 0.25 is 4 times stricter. It is kept deliberately, because
# Ch. 11, Sec. 11.3, printed p. 451, says the lognormal tail is optimistic
# against simulation, and this module reports fade depths from that tail. Do
# not change this value to 1.0 as a "book fix".
WEAK_FLUCTUATION_LIMIT = 0.25


def _sec_zeta(elevation_deg):
    '''Return sec(zeta) = 1/sin(elevation) for the slant path.'''
    return 1.0 / np.sin(np.radians(np.asarray(elevation_deg, dtype=float)))


def plane_wave_scintillation_index(elevation_deg, wavelength, hs, cn2_profile):
    '''
    Return the plane-wave point scintillation index sigma2_I.

    Integrate the Cn2 slant path for the downlink plane wave.

    Parameters:
        elevation_deg : float or numpy.ndarray
            Elevation angle above the horizon [deg].
        wavelength : float
            Optical wavelength [m].
        hs : numpy.ndarray
            Heights above the ground station [m].
        cn2_profile : numpy.ndarray
            Zenith Cn2(h) profile on the hs grid [m^-2/3].

    Returns:
        float or numpy.ndarray
            sigma2_I, broadcast over the elevation shape.

    formula:
        sigma2_I = 2.25 * k^(7/6) * sec(zeta)^(11/6)
                   * trapz( Cn2(h) * h^(5/6), hs )
        Andrews and Phillips, 2nd ed. (2005), Ch. 12, Eq. (38), printed
        p. 495. DOI 10.1117/3.626196.
    '''
    k = 2.0 * np.pi / wavelength
    integral = np.trapz(np.asarray(cn2_profile) * hs ** (5.0 / 6.0), hs)
    return 2.25 * k ** (7.0 / 6.0) * _sec_zeta(elevation_deg) ** (11.0 / 6.0) * integral


# Spatial-wavenumber grid for the aperture-averaging integral. The grid is
# log-spaced from 1e-2 to 1e4 rad/m with 2000 points. The integrand peaks near
# kappa ~ sqrt(k / (h sec)), that is tens of rad/m for the used heights. This
# range holds the peak and the tail. With the aperture filter set to 1 the
# integral reproduces the analytic point index within 5%. See the __main__
# convergence check.
_KAPPA = np.logspace(-2.0, 4.0, 2000)
_KAPPA2 = _KAPPA ** 2


def _aperture_filter(kappa, rx_diameter_m):
    '''
    Return the circular-aperture averaging filter for the wavenumber grid.

    The filter is [2*J1(x) / x]^2 with x = kappa * D / 2. Its limit as x -> 0
    is 1. The code sets the value at x = 0 to 1.

    The filter assumes a uniform circular aperture with no central obscuration.
    An annular (obscured) aperture is not modelled yet.
    Source: Andrews and Phillips, Laser Beam Propagation through Random Media,
    2nd ed. (2005), Ch. 10 (aperture-averaging filter function).
    '''
    x = kappa * rx_diameter_m / 2.0
    out = np.ones_like(x)
    nz = x > 0.0
    out[nz] = (2.0 * j1(x[nz]) / x[nz]) ** 2
    return out


def _scintillation_integral(rx_diameter_m, elevation_deg, wavelength, hs,
                            cn2_profile):
    '''
    Return the aperture-averaged plane-wave flux scintillation index sigma2_I(D).

    Integrate the Rytov double integral over height h and spatial wavenumber
    kappa for the downlink plane wave.

    formula:
        sigma2_I(D) = 8 * pi^2 * k^2 * sec(zeta)
            * INT_h Cn2(h) * [ INT_kappa kappa * 0.033 * kappa^(-11/3)
            * (1 - cos(kappa^2 * h * sec(zeta) / k))
            * (2*J1(kappa*D/2) / (kappa*D/2))^2 dkappa ] dh
        Source: Andrews and Phillips, 2nd ed. (2005), plane-wave scintillation
        index Ch. 12, Eq. (38), printed p. 495, and aperture-averaging filter
        Ch. 10. DOI 10.1117/3.626196. The term
        0.033 * kappa^(-11/3) is the Kolmogorov spectrum with Cn2 factored out.
        No inner scale and no outer scale are used. The distance from turbulence
        at height h to the ground aperture is z = h * sec(zeta). The Fresnel
        filter (1 - cos(kappa^2 z / k)) uses that distance.
    '''
    k = 2.0 * np.pi / wavelength
    sec = _sec_zeta(elevation_deg)
    sec_arr = np.atleast_1d(np.asarray(sec, dtype=float))
    cn2 = np.asarray(cn2_profile, dtype=float)

    # Kolmogorov shape times aperture filter. Both depend on kappa only.
    base = _KAPPA * 0.033 * _KAPPA ** (-11.0 / 3.0) \
        * _aperture_filter(_KAPPA, rx_diameter_m)

    # Loop over the ~20 heights. Vectorise the kappa integral over elevation.
    h_integrand = np.empty((len(hs), sec_arr.size))
    for i, h in enumerate(hs):
        z = h * sec_arr / k                                   # shape (E,)
        fresnel = 1.0 - np.cos(_KAPPA2[:, None] * z[None, :])  # shape (K, E)
        inner = np.trapz(base[:, None] * fresnel, _KAPPA, axis=0)
        h_integrand[i] = cn2[i] * inner
    h_int = np.trapz(h_integrand, hs, axis=0)                 # shape (E,)
    result = 8.0 * np.pi ** 2 * k ** 2 * sec_arr * h_int

    if np.ndim(sec) == 0:
        return float(result[0])
    return result.reshape(np.shape(sec))


def aperture_averaged_scintillation_index(rx_diameter_m, elevation_deg,
                                          wavelength, hs, cn2_profile):
    '''
    Return the aperture-averaged flux scintillation index sigma2_I(D).

    Compute the distributed-path Rytov double integral for a circular ground
    aperture of diameter D. See _scintillation_integral for the formula and the
    citations.

    Parameters:
        rx_diameter_m : float
            Ground receive aperture diameter D [m].
        elevation_deg : float or numpy.ndarray
            Elevation angle above the horizon [deg].
        wavelength : float
            Optical wavelength [m].
        hs : numpy.ndarray
            Heights above the ground station [m].
        cn2_profile : numpy.ndarray
            Zenith Cn2(h) profile on the hs grid [m^-2/3].

    Returns:
        float or numpy.ndarray
            sigma2_I(D), broadcast over the elevation shape.
    '''
    return _scintillation_integral(rx_diameter_m, elevation_deg, wavelength, hs,
                                   cn2_profile)


def aperture_averaging_factor(rx_diameter_m, elevation_deg, wavelength, hs,
                              cn2_profile):
    '''
    Return the plane-wave aperture-averaging factor A (0 < A <= 1).

    Compute the ratio of the aperture-averaged index to the point index. The
    numerator uses the distributed-path integral with the circular-aperture
    filter. The denominator uses the same integral with the aperture filter set
    to 1 (D = 0). That denominator equals the plane-wave point index. So A -> 1
    as D -> 0, and A stays in the range (0, 1].

    Parameters:
        rx_diameter_m : float
            Ground receive aperture diameter D [m].
        elevation_deg : float or numpy.ndarray
            Elevation angle above the horizon [deg].
        wavelength : float
            Optical wavelength [m].
        hs : numpy.ndarray
            Heights above the ground station [m].
        cn2_profile : numpy.ndarray
            Zenith Cn2(h) profile on the hs grid [m^-2/3].

    Returns:
        float or numpy.ndarray
            A, broadcast over the elevation shape.

    formula:
        A(D) = sigma2_I(D) / sigma2_I(0)
        Source: Andrews and Phillips, 2nd ed. (2005), Ch. 10. The large-aperture
        limit gives the known asymptote A ~ D^(-7/3).
    '''
    num = _scintillation_integral(rx_diameter_m, elevation_deg, wavelength, hs,
                                  cn2_profile)
    den = _scintillation_integral(0.0, elevation_deg, wavelength, hs,
                                  cn2_profile)
    return num / den


# ---------------------------------------------------------------------------
# Closed-form single-path aperture averaging.
#
# The functions above integrate the Cn2 profile over the slant path. They are
# the rigorous path. The functions below are the closed-form algebraic
# approximations. They take one path length L and one scalar Cn2. They give a
# fast answer without an integral. Use them for a horizontal path, or for a
# quick estimate. Source: Andrews and Phillips, Laser Beam Propagation through
# Random Media, 2nd ed. (2005), and Churnside, Applied Optics 30 (1991) 1982.
# ---------------------------------------------------------------------------


def _wavenumber(wavelength):
    '''Return the optical wavenumber k = 2*pi/lambda.'''
    return 2.0 * np.pi / wavelength


def sigma1_rytov(cn2, wavelength, path_length_m):
    '''
    Return the plane-wave Rytov standard deviation sigma_1 for a single path.

    formula:
        sigma_1 = ( 1.23 Cn2 k^(7/6) L^(11/6) )^0.5,   k = 2*pi/lambda
    The Rytov variance is sigma_1^2. Source: Andrews and Phillips, 2nd ed.
    (2005), Ch. 5.
    '''
    k = _wavenumber(wavelength)
    return (1.23 * np.asarray(cn2, dtype=float) * k ** (7.0 / 6.0)
            * np.asarray(path_length_m, dtype=float) ** (11.0 / 6.0)) ** 0.5


def coherence_radius(cn2, wavelength, path_length_m):
    '''
    Return the plane-wave coherence radius rho_c for a single path.

    formula:
        rho_c = ( 1.46 Cn2 k^2 L )^(-3/5),   k = 2*pi/lambda
    Source: Andrews and Phillips, 2nd ed. (2005), Ch. 6.
    '''
    k = _wavenumber(wavelength)
    return (1.46 * np.asarray(cn2, dtype=float) * k ** 2
            * np.asarray(path_length_m, dtype=float)) ** (-3.0 / 5.0)


def plane_wave_scintillation_index_closed(cn2, wavelength, path_length_m):
    '''
    Return the point plane-wave scintillation index sigma_I^2 for a single path.

    This is the Andrews closed form. It holds for any turbulence strength. It
    has no inner scale and no outer scale.

    formula:
        sigma_I^2 = exp[ 0.49 s^2 / (1 + 1.11 s^(12/5))^(7/6)
                       + 0.51 s^2 / (1 + 0.69 s^(12/5))^(5/6) ] - 1
    with s = sigma_1 (the Rytov standard deviation). Source: Andrews and
    Phillips, Laser Beam Propagation through Random Media, 2nd ed. (2005),
    Ch. 9, Eq. (47), printed p. 336. DOI 10.1117/3.626196. The same four
    constants are repeated in Ch. 12, Eqs. (40) and (93), and in App. III
    Table VII(b). The d = 0 limit of `aperture_averaged_index_andrews`
    (Ch. 10, Eq. (69)) gives the same four constants.

    NEW HOME: `olb.turbulence.andrews.scintillation.scintillation_index` with
    wave="plane" and regime="strong". That function builds the same result from
    the two log-irradiance variances of Ch. 9, Eqs. (41) and (46), printed
    pp. 335 and 336, which also feed the gamma-gamma distribution.
    '''
    return _andrews_scintillation_index(wavelength, path_length_m, cn2,
                                        wave='plane', regime='strong')


def _d_param(rx_diameter_m, wavelength, path_length_m):
    '''
    Return the aperture parameter d = ( k D^2 / (4 L) )^0.5.

    Source: Andrews and Phillips, 2nd ed. (2005), Ch. 10.
    '''
    k = _wavenumber(wavelength)
    return (k * np.asarray(rx_diameter_m, dtype=float) ** 2
            / (4.0 * np.asarray(path_length_m, dtype=float))) ** 0.5


def aperture_averaged_index_andrews(rx_diameter_m, cn2, wavelength,
                                    path_length_m):
    '''
    Return the aperture-averaged plane-wave flux scintillation index sigma_I^2(D).

    This is the Andrews closed form for a circular aperture of diameter D. It
    holds for any turbulence strength. It has no inner scale and no outer scale.

    formula:
        sigma_I^2(D) = exp[ 0.49 s^2 / (1 + 0.65 d^2 + 1.11 s^(12/5))^(7/6)
                          + 0.51 s^2 (1 + 0.69 s^(12/5))^(-5/6)
                            / (1 + 0.90 d^2 + 0.62 d^2 s^(12/5)) ] - 1
    with s = sigma_1 and d the aperture parameter. Source: Andrews and Phillips,
    2nd ed. (2005), Ch. 10.
    '''
    s = sigma1_rytov(cn2, wavelength, path_length_m)
    s2 = s ** 2
    s125 = s ** (12.0 / 5.0)
    d2 = _d_param(rx_diameter_m, wavelength, path_length_m) ** 2
    term1 = 0.49 * s2 / (1.0 + 0.65 * d2 + 1.11 * s125) ** (7.0 / 6.0)
    term2 = (0.51 * s2 * (1.0 + 0.69 * s125) ** (-5.0 / 6.0)
             / (1.0 + 0.90 * d2 + 0.62 * d2 * s125))
    return np.exp(term1 + term2) - 1.0


def aperture_averaging_factor_weak(rx_diameter_m, wavelength, path_length_m):
    '''
    Return the weak-turbulence aperture-averaging factor A for a Kolmogorov path.

    Use this factor for a small inner scale. It holds for weak turbulence.

    formula:
        A = ( 1 + 1.07 (k D^2 / (4 L))^(7/6) )^(-1)
    Source: Churnside, Applied Optics 30 (1991) 1982,
    DOI 10.1364/AO.30.001982. The constant 1.07 is not in Andrews and
    Phillips. The Andrews counterpart is a DIFFERENT function: 2nd ed. (2005),
    Ch. 10, Eq. (61), printed p. 412, DOI 10.1117/3.626196, which is
    A = [1 + 1.062 k D_G^2/(4 L)]^(-7/6), with the exponent 7/6 outside the
    bracket. The two fits differ by up to 12 %.
    '''
    d2 = _d_param(rx_diameter_m, wavelength, path_length_m) ** 2
    return (1.0 + 1.07 * d2 ** (7.0 / 6.0)) ** (-1.0)


def aperture_averaging_factor_weak_inner(rx_diameter_m, inner_scale_m):
    '''
    Return the weak-turbulence aperture-averaging factor A for a large inner scale.

    Use this factor when the inner scale is much larger than the Fresnel zone. It
    holds for weak turbulence.

    formula:
        A = ( 1 + 2.21 (D / l0)^(7/3) )^(-1)
    with l0 the inner scale. Source: Churnside, Applied Optics 30 (1991) 1982,
    DOI 10.1364/AO.30.001982. The constant 2.21 is not in Andrews and Phillips.
    The Andrews counterpart is a DIFFERENT function: 2nd ed. (2005), Ch. 10,
    Eqs. (62) to (68), printed pp. 412-413, DOI 10.1117/3.626196, which give a
    finite inner scale through the two-scale parameter Q_l.
    '''
    ratio = np.asarray(rx_diameter_m, dtype=float) / np.asarray(inner_scale_m,
                                                                dtype=float)
    return (1.0 + 2.21 * ratio ** (7.0 / 3.0)) ** (-1.0)


def aperture_averaging_factor_strong(rx_diameter_m, cn2, wavelength,
                                     path_length_m):
    '''
    Return the strong-turbulence aperture-averaging factor A for a small inner scale.

    Use this factor when the inner scale is much smaller than the coherence
    length. It holds for strong turbulence.

    formula:
        A = (sI2 + 1) / (2 sI2) * (1 + 0.908 (D / (2 rho_c))^2)^(-1)
          + (sI2 - 1) / (2 sI2) * (1 + 0.162 (k rho_c D / (2 L))^(7/3))^(-1)
    with sI2 the point plane-wave index and rho_c the coherence radius. Source:
    Churnside, Applied Optics 30 (1991) 1982, DOI 10.1364/AO.30.001982. The
    constants 0.908 and 0.162 are not in Andrews and Phillips. Andrews plots
    the Churnside curve (Figs. 10.11 and 10.12, printed p. 418) but does not
    print the formula. The Andrews counterpart is a DIFFERENT function:
    `aperture_averaged_index_andrews` (2nd ed. (2005), Ch. 10, Eq. (69),
    printed p. 413, DOI 10.1117/3.626196), which covers the same regime.
    '''
    k = _wavenumber(wavelength)
    L = np.asarray(path_length_m, dtype=float)
    D = np.asarray(rx_diameter_m, dtype=float)
    si2 = plane_wave_scintillation_index_closed(cn2, wavelength, path_length_m)
    rho_c = coherence_radius(cn2, wavelength, path_length_m)
    term1 = ((si2 + 1.0) / (2.0 * si2)
             * (1.0 + 0.908 * (D / (2.0 * rho_c)) ** 2) ** (-1.0))
    term2 = ((si2 - 1.0) / (2.0 * si2)
             * (1.0 + 0.162 * (k * rho_c * D / (2.0 * L)) ** (7.0 / 3.0)) ** (-1.0))
    return term1 + term2


if __name__ == '__main__':
    # Pure-physics self-check. Use plain numeric inputs; this module must not
    # import the scenario or the geometry. The Site defaults set the Cn2 profile
    # (Bufton wind rms 21 m/s, HV57 ground scale 1.7e-14).
    lam = 1550e-9
    hs = DEFAULT_HS
    cn2 = get_c2n(hs, 21.0, 1.7e-14)
    D = 0.7

    # Longer slant path -> larger scintillation index.
    s_30 = plane_wave_scintillation_index(30.0, lam, hs, cn2)
    s_90 = plane_wave_scintillation_index(90.0, lam, hs, cn2)
    assert s_30 > s_90, (s_30, s_90)

    # Aperture averaging reduces the scintillation: 0 < A < 1, sigma2_P < sigma2_I.
    A_30 = aperture_averaging_factor(D, 30.0, lam, hs, cn2)
    assert 0.0 < A_30 < 1.0, A_30
    assert A_30 * s_30 < s_30

    # Convergence: with D -> 0 the integral reproduces the analytic point index.
    # This proves the leading constant 8*pi^2*0.033 is right.
    idx0 = aperture_averaged_scintillation_index(1e-6, 30.0, lam, hs, cn2)
    conv_pct = abs(idx0 - s_30) / s_30 * 100.0
    assert conv_pct < 5.0, (idx0, s_30, conv_pct)

    # A(D) is in (0, 1], falls as D grows, and A -> 1 as D -> 0.
    A_small = aperture_averaging_factor(1e-6, 30.0, lam, hs, cn2)
    A_big = aperture_averaging_factor(2.0, 30.0, lam, hs, cn2)
    assert abs(A_small - 1.0) < 1e-3, A_small
    assert 0.0 < A_big < A_30 < A_small <= 1.0, (A_big, A_30, A_small)

    # Large-aperture asymptote: A ~ D^(-7/3). Test the local log-log slope.
    D1, D2 = 3.0, 6.0
    A_D1 = aperture_averaging_factor(D1, 30.0, lam, hs, cn2)
    A_D2 = aperture_averaging_factor(D2, 30.0, lam, hs, cn2)
    slope = np.log(A_D2 / A_D1) / np.log(D2 / D1)
    assert abs(slope - (-7.0 / 3.0)) < 0.3, slope

    # Closed-form single-path aperture averaging. Use one path and one Cn2.
    L = 2400.0
    cn2_flat = 1e-15
    A_w = aperture_averaging_factor_weak(0.7, lam, L)
    A_s = aperture_averaging_factor_strong(0.7, cn2_flat, lam, L)
    A_wi = aperture_averaging_factor_weak_inner(0.7, 5e-3)
    # Each factor stays in (0, 1] and a larger aperture averages more.
    for A in (A_w, A_s, A_wi):
        assert 0.0 < A <= 1.0, A
    assert aperture_averaging_factor_weak(1.4, lam, L) < A_w
    # The closed-form aperture-averaged index is below the point index.
    si2 = plane_wave_scintillation_index_closed(cn2_flat, lam, L)
    si2_D = aperture_averaged_index_andrews(0.7, cn2_flat, lam, L)
    assert 0.0 < si2_D < si2, (si2_D, si2)

    print(f"sigma2_I  30 deg = {s_30:.4f}   90 deg = {s_90:.4f}")
    print(f"index(D->0) 30 deg = {idx0:.4f}   point = {s_30:.4f}   "
          f"conv = {conv_pct:.2f}%")
    print(f"A large-aperture log-log slope = {slope:.3f} (target -2.333)")
    print(f"A(D=0.7m) 30 deg = {A_30:.4f}   sigma2_P = {A_30 * s_30:.4f}")
    print(f"closed-form A_weak={A_w:.4f} A_strong={A_s:.4f} "
          f"A_weak_inner={A_wi:.4f}")
    print("self-check passed")
