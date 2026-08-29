'''
Angular anisoplanatism over a finite aperture.

Two wavefronts come to one aperture from two slightly different angles. The
turbulence gives them different phase. This module gives the mean-square phase
difference between them. This is the angular anisoplanatic error. It sets the
limit of a system that measures the turbulence on one source and corrects a
second source at a small angle from it. A laser guide star, a beacon uplink, and
a point-ahead downlink all have this error.

The functions are pure. They take numeric values and numpy arrays, and they
return numbers. This module does not import the scenario, the models, or the
links.

The classical law is
    sigma^2 = (theta / theta0)^(5/3)                                    Eq. (26)
with theta0 the isoplanatic angle
    theta0 = ( C1 k0^2 INT Cn2(h) S^(5/3) dh )^(-3/5),   S = h * airmass  Eq. (27)
    C1 = 2 (2 pi)^(8/3) C_A |HJ1(8/3,0,1)| = 2.914381
That law is for an aperture of zero size. It keeps the PISTON of the phase
difference. A real aperture does not lose light to piston. A real aperture with
a tip-tilt mirror does not lose light to tilt either. So the classical law reads
too large an error. The paper shows the size of that error in Fig. 1. For a
small aperture and a large angle the classical law is up to about 1 order of
magnitude too large.

This module gives the finite-aperture result. It removes the piston, or the
piston and the two tilts, over the aperture:
    sigma^2 = 2 (2 pi)^(8/3) C_A k0^2 R^(5/3) airmass
              * INT dh Cn2(h) I( S theta / R )                       Eq. (29)
    I(b) = INT_0^inf du u^(-8/3) (1 - J0(b u)) * M(u)                 Eq. (36)
Here R = D/2 is the aperture radius, and M(u) is a modal weight.

Modal weight, decorrelation, and adaptive optics.
    The phase difference splits into Zernike radial orders n. The paper gives
    the weight of each order (Eq. A11):
        p_n(u) = 4 (n+1)^2 ( J_{n+1}(u) / u )^2
    The order n=0 is the piston. The order n=1 is the two tilts. The full weight
    of all orders is 1 (a closure relation), so
        M(u) = 1 - p_0(u) - p_1(u)
    is the piston-and-tilt-removed weight of Eq. (36).

    Each order's variance is the DECORRELATION residual between the two
    directions. Write it for equal per-path statistics:
        <[a_n(I) - a_n(II)]^2> = 2 sigma_n^2 ( 1 - rho_n )
    where rho_n is the correlation of order n across the angle. An adaptive-optics
    system senses order n on one source and applies it to the other. It removes
    the correlated part rho_n. The decorrelated part 2 sigma_n^2 (1 - rho_n)
    stays. So the error of a system that corrects the orders 2..max_order is the
    BAND weight
        M(u) = p_2(u) + p_3(u) + ... + p_{max_order}(u)
    The band grows with max_order because each added order brings its own
    decorrelation residual (small for a well-correlated low order, up to twice
    the mode variance for a fully decorrelated one; the paper notes this "twice"
    limit, p. 352). The band goes up to the infinite-order value above. Set
    max_order to an integer for a finite adaptive-optics system. Leave it None for
    the ideal, infinite-order limit. This is the paper's Fig. 2 and its
    frequency-restricted variance, Eq. (43).

The code integrates the Bessel kernel of Eq. (36) directly. It does not use the
hypergeometric 3F2 series of Eqs. (31)-(32) in the paper. The two give the same
number, but the direct integral is shorter and it has no series-convergence
limit.

Source of all equations in this module:
    J. Stone, P. H. Hu, S. P. Mills and S. Ma, "Anisoplanatic effects in
    finite-aperture optical systems," J. Opt. Soc. Am. A 11(1), 347-357 (1994).
    DOI: 10.1364/JOSAA.11.000347
'''

import warnings

import numpy as np
from scipy.integrate import quad
from scipy.special import jv

from ..assumptions import (BEAM_PLANE_WAVE, SPECTRUM_KOLMOGOROV, Constraint,
                           assumes)

# The three slant functions below share these assumptions. The phase-difference
# result is a wavefront quantity (a phase structure function), so it does not
# carry a scintillation regime. There is no numeric validity gate, so no
# constraint carries a check.
_ISOPLANATISM = Constraint(
    "isoplanatism",
    "The model describes angular anisoplanatism: the phase error between two "
    "directions at a small angle grows over the isoplanatic angle theta0.",
    "10.1364/JOSAA.11.000347", "Eqs. (26), (27)")

# The slant scaling uses the plane-parallel airmass 1/sin(elevation). It has no
# Earth curvature. Source: Andrews and Phillips, 2nd ed. (2005),
# DOI 10.1117/3.626196, Ch. 12, printed p. 481.
_PLANE_PARALLEL = Constraint(
    "geometry",
    "The slant path uses the plane-parallel airmass 1/sin(elevation). It has no "
    "Earth curvature, so it breaks near the horizon.",
    "10.1117/3.626196", "Ch. 12, printed p. 481")

# The finite-aperture result is the pure angular case: both sources are at
# infinity (plane waves), so the two apertures have equal radius R1 = R2 = R.
_PURE_ANGULAR = Constraint(
    "field-region",
    "Both sources are at infinity (plane waves), so the two apertures have equal "
    "radius R1 = R2 = R. It is the pure angular case, not a finite-range beacon.",
    "10.1364/JOSAA.11.000347", "Eqs. (29), (36), R1(S) = R2(S) = R")

# The classical law is a zero-size aperture. It keeps the piston and the tilt.
_ZERO_APERTURE = Constraint(
    "aperture-order",
    "The classical law is a zero-size aperture. It keeps the piston and the "
    "tilt, so it reads up to about 10 times too large an error over a real "
    "aperture. Use anisoplanatic_phase_variance for a finite aperture.",
    "10.1364/JOSAA.11.000347", "Eqs. (1), (26); Fig. 1")

# Turbulence constant of the phase structure function. Eq. (14) of the paper:
# C_A = (5/36) 2^(1/3) / ( pi^(5/3) Gamma(1/3) ).
# DOI: 10.1364/JOSAA.11.000347
C_A = 0.0096932

# |HJ1(8/3, 0, 1)|, the value of the u-integral of Eq. (36) with no mode
# removed. See Ref. note 15 of the paper. DOI: 10.1364/JOSAA.11.000347
HJ1_8_3 = 1.11833

# (2 pi)^(8/3). It comes from the spatial-frequency scale of Eq. (29).
_TWO_PI_83 = (2.0 * np.pi) ** (8.0 / 3.0)


def _order_power(n, u):
    '''
    Return the anisoplanatic weight of one Zernike radial order.

    formula:
        p_n(u) = 4 (n+1)^2 ( J_{n+1}(u) / u )^2
    Source: Eq. (A11) of Stone et al. (1994), for the pure angular case with
    equal aperture radii. DOI: 10.1364/JOSAA.11.000347

    The ratio J_{n+1}(u)/u has a limit at u = 0. For the piston (n = 0) the limit
    is 1/2. For every higher order the limit is 0.

    Parameters:
        n : int
            Zernike radial order. 0 is the piston. 1 is the two tilts.
        u : float or numpy.ndarray
            Spatial frequency times the aperture radius. It has no unit.

    Returns:
        float or numpy.ndarray
            p_n(u), same shape as u.
    '''
    u = np.asarray(u, dtype=float)
    safe_u = np.where(u == 0.0, 1.0, u)
    ratio = jv(n + 1, u) / safe_u
    if n == 0:
        ratio = np.where(u == 0.0, 0.5, ratio)   # J1(u)/u -> 1/2 as u -> 0
    return 4.0 * (n + 1) ** 2 * ratio ** 2


def _mode_factor(u, n_lo, max_order):
    '''
    Return the modal weight M(u) of the anisoplanatic integrand.

    n_lo is the lowest radial order that counts as error. Use 0 to keep the
    piston, 1 to keep the tilt, and 2 to remove the piston and the tilt.

    max_order is the highest radial order that the correction touches. Use None
    for an ideal, infinite-order correction; then the weight uses the closed form
        M(u) = 1 - p_0(u) - ... - p_{n_lo-1}(u)
    Use an integer for a finite adaptive-optics system; then the weight is the
    band sum
        M(u) = p_{n_lo}(u) + ... + p_{max_order}(u)
    Source: Eqs. (29), (36) and (A11) of Stone et al. (1994).
    DOI: 10.1364/JOSAA.11.000347

    Parameters:
        u : float or numpy.ndarray
            Spatial frequency times the aperture radius. It has no unit.
        n_lo : int
            Lowest radial order that counts as error (0, 1, or 2).
        max_order : int or None
            Highest corrected radial order, or None for all orders.

    Returns:
        float or numpy.ndarray
            M(u), same shape as u.
    '''
    u = np.asarray(u, dtype=float)
    if max_order is None:
        factor = np.ones_like(u)
        for n in range(n_lo):
            factor = factor - _order_power(n, u)
        return factor
    factor = np.zeros_like(u)
    for n in range(n_lo, max_order + 1):
        factor = factor + _order_power(n, u)
    return factor


def _inner_integral(beta, n_lo, max_order):
    '''
    Return the spatial-frequency integral of the anisoplanatic variance.

    formula:
        I(beta) = INT_0^inf du u^(-8/3) (1 - J0(beta u)) * M(u)
    Source: Eqs. (29) and (36) of Stone et al. (1994).
    DOI: 10.1364/JOSAA.11.000347

    The integral converges. Near u = 0 the integrand goes as u^(4/3) when the
    piston is removed, and as u^(-2/3) when it is not. At large u the integrand
    goes as u^(-8/3). The code splits the range into parts. The integrand
    oscillates, so a single quad call over the full range is not reliable. A high
    max_order puts the Bessel peak near u = max_order, so the code carries the
    upper limit past that peak.

    Parameters:
        beta : float
            The scaled angular offset S * theta / R. It has no unit.
        n_lo : int
            Lowest radial order that counts as error (0, 1, or 2).
        max_order : int or None
            Highest corrected radial order, or None for all orders.

    Returns:
        float
            I(beta). It has no unit.
    '''
    def f(u):
        return (u ** (-8.0 / 3.0) * (1.0 - jv(0, beta * u))
                * _mode_factor(u, n_lo, max_order))

    # The J_{n+1} peak sits near u = max_order, so carry the top past it.
    top = 500.0 if max_order is None else max(500.0, 4.0 * (max_order + 2))
    edges = [0.0, 1.0, 50.0, top]
    with warnings.catch_warnings():
        # quad reports roundoff on the oscillating tail. The result is good.
        warnings.simplefilter('ignore')
        total = 0.0
        for a, b in zip(edges[:-1], edges[1:]):
            v, _ = quad(f, a, b, limit=400)
            total += v
    return total


def max_radial_order(n_zernike_modes):
    '''
    Return the highest complete Zernike radial order in a set of the first
    n_zernike_modes Noll modes.

    The count of modes through radial order n is (n+1)(n+2)/2. So this returns the
    largest n with (n+1)(n+2)/2 <= n_zernike_modes. Use it to turn an
    adaptive-optics mode count into the max_order of anisoplanatic_phase_variance.
    Source: R. J. Noll, "Zernike polynomials and atmospheric turbulence,"
    J. Opt. Soc. Am. 66(3), 207-211 (1976). DOI: 10.1364/JOSA.66.000207

    Parameters:
        n_zernike_modes : int
            Number of Noll Zernike modes that the system corrects.

    Returns:
        int
            Highest complete radial order. 0 for the piston only.
    '''
    n = 0
    while (n + 2) * (n + 3) // 2 <= n_zernike_modes:
        n += 1
    return n


@assumes(_ISOPLANATISM, _PLANE_PARALLEL, spectrum=SPECTRUM_KOLMOGOROV)
def isoplanatic_angle(hs, cn2_profile, wavelength, elevation_deg=90.0):
    '''
    Return the classical isoplanatic angle theta0.

    Integrate the zenith Cn2 profile with the h^(5/3) weight. Scale the path to
    the slant path with the airmass.

    formula:
        theta0 = ( C1 k0^2 INT Cn2(h) S^(5/3) dh )^(-3/5)
        S      = h * airmass,   airmass = 1 / sin(elevation)
        k0     = 2 * pi / lambda
        C1     = 2 (2 pi)^(8/3) C_A |HJ1(8/3,0,1)| = 2.914381
    Source: Eq. (27) of Stone et al. (1994). The paper gives the value 2.914381
    for C1 in the text after that equation.
    DOI: 10.1364/JOSAA.11.000347

    Parameters:
        hs : numpy.ndarray
            Heights above the ground station [m], ascending.
        cn2_profile : numpy.ndarray
            Zenith Cn2(h) profile on the hs grid [m^-2/3].
        wavelength : float
            Optical wavelength [m].
        elevation_deg : float
            Elevation angle above the horizon [deg]. 90 is the zenith.

    Returns:
        float
            theta0 [rad].
    '''
    hs = np.asarray(hs, dtype=float)
    cn2 = np.asarray(cn2_profile, dtype=float)
    k0 = 2.0 * np.pi / wavelength
    airmass = 1.0 / np.sin(np.radians(elevation_deg))
    # S^(5/3) = (h * airmass)^(5/3). The path element adds one more airmass, so
    # the total slant factor is airmass^(8/3).
    integral = np.trapz(cn2 * hs ** (5.0 / 3.0), hs) * airmass ** (8.0 / 3.0)
    c1 = 2.0 * _TWO_PI_83 * C_A * HJ1_8_3
    return float((c1 * k0 ** 2 * integral) ** (-3.0 / 5.0))


# Lowest radial order that counts as error, for each value of `remove`. 0 keeps
# the piston, 1 keeps the tilt, 2 removes the piston and the two tilts.
# Source: Eq. (36) of Stone et al. (1994). DOI: 10.1364/JOSAA.11.000347
_REMOVE_NLO = {'none': 0, 'piston': 1, 'piston_tilt': 2}


@assumes(_ISOPLANATISM, _PLANE_PARALLEL, _PURE_ANGULAR,
         beam_type=BEAM_PLANE_WAVE, spectrum=SPECTRUM_KOLMOGOROV)
def anisoplanatic_phase_variance(D, theta, hs, cn2_profile, wavelength,
                                 remove='piston_tilt', max_order=None,
                                 elevation_deg=90.0):
    '''
    Return the angular anisoplanatic phase variance over a finite aperture.

    Two plane waves come to the same aperture from two directions. The angle
    between them is theta. This function gives the mean-square phase difference
    between the two wavefronts over the aperture, for the modes that `remove` and
    `max_order` select.

    formula:
        sigma^2 = 2 (2 pi)^(8/3) C_A k0^2 R^(5/3) airmass
                  * INT dh Cn2(h) I( S theta / R )
        R    = D / 2,   S = h * airmass,   airmass = 1 / sin(elevation)
        k0   = 2 * pi / lambda
        I(b) = the modal spatial-frequency integral, see _inner_integral
    Source: Eqs. (29) and (36) of Stone et al. (1994), with R1(S) = R2(S) = R.
    That is the pure angular case: both sources are at infinity.
    DOI: 10.1364/JOSAA.11.000347

    Parameters:
        D : float
            Aperture diameter [m].
        theta : float
            Angle between the two directions [rad].
        hs : numpy.ndarray
            Heights above the ground station [m], ascending. Keep hs[0] above 0,
            because the ground layer adds no anisoplanatic error.
        cn2_profile : numpy.ndarray
            Zenith Cn2(h) profile on the hs grid [m^-2/3].
        wavelength : float
            Optical wavelength [m].
        remove : str
            "none" keeps all the modes. It gives the classical result of
            Eq. (26). "piston" removes the piston. "piston_tilt" removes the
            piston and the two tilts. Use "piston_tilt" for a terminal that has
            a tip-tilt mirror or a fast steering mirror.
        max_order : int or None
            The highest Zernike radial order that the correction touches. None
            (the default) is the ideal, infinite-order correction; then the
            result is the residual with the `remove` modes taken out. An integer
            is a finite adaptive-optics system that corrects the radial orders up
            to max_order; then the result is the band sum of those orders (the
            paper's Fig. 2). The band grows with max_order because each added
            order brings its own decorrelation residual between the two
            directions. Use max_radial_order to turn an adaptive-optics mode count
            into this value.
        elevation_deg : float
            Elevation angle above the horizon [deg]. 90 is the zenith.

    Returns:
        float
            sigma^2 [rad^2].

    Raises:
        ValueError
            If `remove` is not "none", "piston", or "piston_tilt", or if
            max_order is a negative integer.
    '''
    if remove not in _REMOVE_NLO:
        raise ValueError(
            "remove must be 'none', 'piston', or 'piston_tilt', "
            f"not {remove!r}"
        )
    n_lo = _REMOVE_NLO[remove]
    if max_order is not None and max_order < 0:
        raise ValueError(f"max_order must be a non-negative int or None, "
                         f"not {max_order!r}")

    hs = np.asarray(hs, dtype=float)
    cn2 = np.asarray(cn2_profile, dtype=float)
    k0 = 2.0 * np.pi / wavelength
    R = D / 2.0
    airmass = 1.0 / np.sin(np.radians(elevation_deg))

    # beta = S * theta / R at each height of the grid.
    betas = hs * airmass * theta / R
    inner = np.array([_inner_integral(b, n_lo, max_order) for b in betas])

    prefactor = 2.0 * _TWO_PI_83 * C_A * k0 ** 2 * R ** (5.0 / 3.0) * airmass
    return float(prefactor * np.trapz(cn2 * inner, hs))


@assumes(_ISOPLANATISM, _PLANE_PARALLEL, _ZERO_APERTURE,
         beam_type=BEAM_PLANE_WAVE, spectrum=SPECTRUM_KOLMOGOROV)
def anisoplanatic_phase_variance_classic(theta, hs, cn2_profile, wavelength,
                                         elevation_deg=90.0):
    '''
    Return the classical anisoplanatic phase variance.

    This is the reference law. It keeps the piston and the tilt, so it reads a
    larger error than the finite aperture has. The result does not depend on the
    aperture diameter. Compare it with anisoplanatic_phase_variance.

    formula:
        sigma^2 = (theta / theta0)^(5/3)
    Source: Eqs. (1) and (26) of Stone et al. (1994), with theta0 from Eq. (27).
    DOI: 10.1364/JOSAA.11.000347

    Parameters:
        theta : float
            Angle between the two directions [rad].
        hs : numpy.ndarray
            Heights above the ground station [m], ascending.
        cn2_profile : numpy.ndarray
            Zenith Cn2(h) profile on the hs grid [m^-2/3].
        wavelength : float
            Optical wavelength [m].
        elevation_deg : float
            Elevation angle above the horizon [deg]. 90 is the zenith.

    Returns:
        float
            sigma^2 [rad^2].
    '''
    theta0 = isoplanatic_angle(hs, cn2_profile, wavelength, elevation_deg)
    return float((theta / theta0) ** (5.0 / 3.0))


if __name__ == '__main__':
    # Pure-physics self-check. Build the HV5/7 profile that the paper uses.
    from .profiles import DEFAULT_HS, get_c2n

    hs = DEFAULT_HS
    cn2 = get_c2n(hs, 21.0, 1.7e-14)          # HV5/7: 21 m/s wind, 1.7e-14 ground
    lam = 0.5e-6                              # the paper works at 0.5 um

    # The paper quotes theta0 of about 7 urad for HV5/7 at 0.5 um at the zenith.
    theta0 = isoplanatic_angle(hs, cn2, lam, elevation_deg=90.0)
    assert np.isclose(theta0, 7e-6, rtol=0.10), theta0

    # theta0 falls toward the horizon, because the slant path holds more
    # turbulence.
    theta0_30 = isoplanatic_angle(hs, cn2, lam, elevation_deg=30.0)
    assert theta0_30 < theta0, (theta0_30, theta0)

    # With no mode removed, Eq. (29) must give the classical Eq. (26).
    for D, th in ((1.0, 10e-6), (0.5, 20e-6)):
        full = anisoplanatic_phase_variance(D, th, hs, cn2, lam, remove='none')
        classic = anisoplanatic_phase_variance_classic(th, hs, cn2, lam)
        assert np.isclose(full, classic, rtol=0.01), (D, th, full, classic)

    # The paper's Fig. 1c: D = 0.5 m at 20 urad gives about 1 rad^2 after the
    # removal of the piston and the tilt.
    v_fig1c = anisoplanatic_phase_variance(0.5, 20e-6, hs, cn2, lam,
                                           remove='piston_tilt')
    assert np.isclose(v_fig1c, 1.0, rtol=0.50), v_fig1c

    # Each removed mode takes variance away, so the order is fixed.
    D, th = 1.0, 15e-6
    v_none = anisoplanatic_phase_variance(D, th, hs, cn2, lam, remove='none')
    v_pist = anisoplanatic_phase_variance(D, th, hs, cn2, lam, remove='piston')
    v_ptilt = anisoplanatic_phase_variance(D, th, hs, cn2, lam,
                                           remove='piston_tilt')
    assert v_none >= v_pist >= v_ptilt, (v_none, v_pist, v_ptilt)

    # The paper's Fig. 1: the classical law reads too large an error, and it
    # gets worse as the aperture gets smaller. A small aperture at a large angle
    # shows more than 1 order of magnitude of overprediction.
    classic_40 = anisoplanatic_phase_variance_classic(40e-6, hs, cn2, lam)
    ratios = [classic_40 / anisoplanatic_phase_variance(d, 40e-6, hs, cn2, lam)
              for d in (0.2, 0.5, 1.0, 2.0)]
    assert all(r > 1.0 for r in ratios), ratios
    assert ratios[0] > ratios[1] > ratios[2] > ratios[3], ratios
    assert ratios[0] > 10.0, ratios[0]

    # A bad value of remove is an error.
    try:
        anisoplanatic_phase_variance(D, th, hs, cn2, lam, remove='tilt')
    except ValueError:
        pass
    else:
        raise AssertionError('remove must reject an unknown mode set')

    # --- finite-order (adaptive-optics band) self-check ---------------------
    # Closure: the weights of all radial orders sum to 1 at every u. This is the
    # relation that lets 1 - p0 - p1 stand for the sum of all higher orders.
    for uu in (0.5, 2.0, 5.0, 10.0):
        s = sum(_order_power(n, uu) for n in range(0, 400))
        assert np.isclose(s, 1.0, atol=2e-2), (uu, s)

    # A finite adaptive-optics band 2..max_order grows with max_order and goes up
    # to the infinite-order piston+tilt result. Each added order brings its own
    # decorrelation residual between the two directions.
    D, th = 1.0, 10e-6
    inf_pt = anisoplanatic_phase_variance(D, th, hs, cn2, lam, remove='piston_tilt')
    orders = (2, 3, 5, 10, 30, 90)
    bands = [anisoplanatic_phase_variance(D, th, hs, cn2, lam,
             remove='piston_tilt', max_order=m) for m in orders]
    assert all(x < y for x, y in zip(bands, bands[1:])), bands   # grows with order
    assert all(b < inf_pt for b in bands), (bands, inf_pt)       # each below ideal
    assert np.isclose(bands[-1], inf_pt, rtol=0.10), (bands[-1], inf_pt)  # converges

    # Correcting nothing above the tilt leaves no anisoplanatic error.
    zero = anisoplanatic_phase_variance(D, th, hs, cn2, lam,
                                        remove='piston_tilt', max_order=1)
    assert np.isclose(zero, 0.0, atol=1e-9), zero

    # A negative max_order is an error.
    try:
        anisoplanatic_phase_variance(D, th, hs, cn2, lam, max_order=-1)
    except ValueError:
        pass
    else:
        raise AssertionError('max_order must reject a negative order')

    # Noll mode count -> highest complete radial order.
    assert max_radial_order(1) == 0    # piston only
    assert max_radial_order(3) == 1    # through tilt
    assert max_radial_order(6) == 2    # through defocus + astigmatism
    assert max_radial_order(10) == 3
    assert max_radial_order(20) == 4   # 21 modes fill order 5, so 20 stops at 4

    print(f"theta0 = {theta0 * 1e6:.2f} urad  (paper: ~7)")
    print(f"theta0 @30deg = {theta0_30 * 1e6:.2f} urad")
    print(f"{'theta[urad]':>11} {'none':>8} {'pist_rm':>8} {'pist+tilt':>10}"
          f"  D = {D:.1f} m")
    for t in (5e-6, 10e-6, 15e-6, 20e-6):
        a = anisoplanatic_phase_variance(D, t, hs, cn2, lam, remove='none')
        b = anisoplanatic_phase_variance(D, t, hs, cn2, lam, remove='piston')
        c = anisoplanatic_phase_variance(D, t, hs, cn2, lam,
                                         remove='piston_tilt')
        print(f"{t * 1e6:11.0f} {a:8.3f} {b:8.3f} {c:10.3f}")
    print(f"D=0.5 m at 20 urad, piston+tilt removed = {v_fig1c:.3f} rad^2  "
          "(paper Fig. 1c: ~1.0)")
    print("classical overprediction at 40 urad, D = 0.2/0.5/1/2 m: "
          + "  ".join(f"{r:.1f}x" for r in ratios))
    print(f"AO band 2..n at D={D:.1f} m, {th*1e6:.0f} urad (grows to "
          f"{inf_pt:.3f} at infinite order):")
    print("   n:     " + "  ".join(f"{m:>6d}" for m in orders) + "     inf")
    print("   var: " + "  ".join(f"{b:6.3f}" for b in bands)
          + f"   {inf_pt:6.3f}")

    # --- assumptions layer ---------------------------------------------------
    from ..assumptions import trace_assumptions

    # (1) Value parity: a decorated function returns the identical value with and
    #     without a collection context.
    outside = isoplanatic_angle(hs, cn2, lam, 90.0)
    with trace_assumptions():
        inside = isoplanatic_angle(hs, cn2, lam, 90.0)
    assert outside == inside, (outside, inside)

    # (2) Registration: inside a context the three slant functions register their
    #     sources and kinds, and the physics layer emits no warning. (There is no
    #     numeric validity gate, so no violation block.)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        with trace_assumptions() as trace:
            isoplanatic_angle(hs, cn2, lam, 90.0)
            anisoplanatic_phase_variance(1.0, 10e-6, hs, cn2, lam)
            anisoplanatic_phase_variance_classic(10e-6, hs, cn2, lam)
    for name in ("isoplanatic_angle", "anisoplanatic_phase_variance",
                 "anisoplanatic_phase_variance_classic"):
        assert f"{__name__}.{name}" in trace.records, trace.records
    kinds = {c.kind for rec in trace.records.values() for c in rec.constraints}
    assert {"isoplanatism", "geometry", "field-region", "aperture-order"} <= kinds, kinds
    assert len(caught) == 0, "the anisoplanatism physics must not warn"

    print('self-check passed')
