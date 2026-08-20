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
    2nd ed. (2005), Eq. (12.44). Here k = 2*pi/lambda, zeta is the zenith angle,
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
    Ch. 12 (Eq. 12.44) and aperture-averaging filter Ch. 10. The filter
    (2*J1(x)/x)^2 assumes a uniform circular aperture with no central
    obscuration. The aperture-averaging factor is A = sigma2_P / sigma2_I. It
    obeys 0 < A <= 1, A -> 1 as D -> 0, and the large-aperture asymptote
    A ~ D^(-7/3).
'''

import numpy as np
from scipy.special import j1

from .profiles import DEFAULT_HS, get_c2n

# Rytov weak-fluctuation limit for the plane-wave scintillation index. The
# lognormal irradiance PDF is trusted for sigma2_I below approximately 0.25.
# Above it focusing and saturation make the lognormal model depart from data.
# Source: Andrews and Phillips, Laser Beam Propagation through Random Media,
# 2nd ed. (2005), Ch. 5 (weak-fluctuation regime, sigma2_R < ~0.25).
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
        Andrews and Phillips, 2nd ed. (2005), Eq. (12.44).
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
        index Ch. 12 (Eq. 12.44) and aperture-averaging filter Ch. 10. The term
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

    print(f"sigma2_I  30 deg = {s_30:.4f}   90 deg = {s_90:.4f}")
    print(f"index(D->0) 30 deg = {idx0:.4f}   point = {s_30:.4f}   "
          f"conv = {conv_pct:.2f}%")
    print(f"A large-aperture log-log slope = {slope:.3f} (target -2.333)")
    print(f"A(D=0.7m) 30 deg = {A_30:.4f}   sigma2_P = {A_30 * s_30:.4f}")
    print("self-check passed")
