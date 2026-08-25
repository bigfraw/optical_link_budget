'''
Beam wander and the wander-induced pointing error of Andrews and Phillips.

Turbulent eddies that are LARGER than the beam refract the whole beam. They do
not break it up. So the instantaneous centre of the beam ("hot spot") moves
around the boresight in the receiver plane. This is beam wander. Eddies that are
SMALLER than the beam spread it and make it breathe. This module holds the
large-scale (refractive) part only.

The module gives four quantities:
    <r_c^2>    the beam-wander displacement variance      Ch. 6, Eqs. (93)-(99)
    sigma_pe^2 the wander-induced pointing-error variance Ch. 8, Eqs. (36)-(38)
    W_LT       the long-term beam radius                  Ch. 6, Eq. (86)
    W_ST       the short-term beam radius                 Ch. 6, Eq. (100)
It gives each of the first two on a HOMOGENEOUS (constant Cn2) path and on a
SLANT path with a Cn2(h) profile (Ch. 12, Eqs. (50) and (53)).

Source of every equation: Andrews and Phillips, Laser Beam Propagation through
Random Media, 2nd ed. (SPIE Press, 2005), DOI 10.1117/3.626196.

CONVENTION (read this before you use a number from here).
    <r_c^2> and sigma_pe^2 are TWO-DIMENSIONAL (RADIAL) displacement variances.
    They are the variance of the MAGNITUDE of the hot-spot displacement in the
    receiver plane, not the variance along one Cartesian axis. Andrews states
    the two choices in Ch. 6, Sec. 6.6, printed p. 201 ("the variance of the hot
    spot displacement along an axis or ... the variance of the magnitude of the
    hot spot displacement"), and then uses the magnitude form throughout. Two
    facts pin it down:
        1. Ch. 6, Eq. (100), printed p. 205, adds <r_c^2> to W_ST^2 with the
           FACTOR 1. A beam radius is a radial quantity, so <r_c^2> must be
           radial too.
        2. Worked Example 2, Ch. 6, printed p. 215, calls the number "the rms
           displacement of the beam hot spot" and gets 3.35 cm from
           2.42 Cn2 L^3 W0^(-1/3). The self-check reproduces it.
    For a per-axis (one Cartesian component) variance, divide by 2. This module
    NEVER does that division for you. See Conflict C-03 in
    docs/andrews-crosscheck.md.

PLANE OF REFERENCE. z (or L) is measured FROM THE TRANSMITTER. The normalised
variable is xi = 1 - z/L, so xi = 1 at the transmitter and xi = 0 at the
receiver. The weight xi^2 puts the wander on the turbulence NEAR THE
TRANSMITTER. This is correct for beam wander, because the tilt that moves the
beam is applied at the transmitter (Ch. 8, Sec. 8.3, printed p. 272, the
reciprocity argument).

APPROXIMATIONS. Every form here uses the geometrical-optics approximation of
Ch. 6, Eq. (92), printed p. 203, and drops the diffraction part of the large-
scale filter (Ch. 6, Eq. (91), printed p. 203). So the beam radius along the
path is the refractive radius W(z) = W0 |Theta0 + Theta0_bar xi|, and the
wavelength does NOT enter the Kolmogorov result. The forms hold under weak
irradiance fluctuations.

This module holds physics only. It returns no decibels. It imports numpy, scipy,
and the sibling module beam.py.
'''

import numpy as np
from scipy.integrate import quad
from scipy.special import hyp2f1

from .beam import beam_params, wavenumber

# The prefactor of Ch. 6, Eq. (93), printed p. 203. The analytic value is
# 2 * 4 pi^2 * 0.033 * (1/2) * Gamma(1/6) = 7.2520. The book prints 7.25.
#
# This is the Andrews beam-wave SPECTRAL-FILTER constant. It is 3.50 times the
# 2.07 IMAGE-MOTION constant of the Dios/Belmonte kernel that the olb uplink
# chain uses. Both are the same RADIAL quantity; the ratio is purely the leading
# constant (see the C-01 self-check below). Belmonte, Appl. Opt. 39, 5426 (2000),
# DOI 10.1364/AO.39.005426, Eq. (21), validates the 2.07 form against a split-
# step simulation, so a simulation-validated wander must use the Dios kernel
# route, NOT this 7.25 form. See Conflict C-01 in docs/andrews-crosscheck.md.
WANDER_CONSTANT = 7.25

# The collimated reduction of Ch. 6, Eq. (94), printed p. 204. The book prints
# 2.42 for WANDER_CONSTANT / 3 = 2.4167, a 0.07 % rounding.
WANDER_CONSTANT_COLLIMATED = WANDER_CONSTANT / 3.0

# The scaling constant of the jitter filter cut-off kappa_r = C_r / r0.
# Source: Ch. 8, text below Eq. (35), printed p. 273: "C_r is a scaling constant
# typically on the order C_r ~ 2 pi".
CR_DEFAULT = 2.0 * np.pi

# The scaling constant of the outer-scale cut-off kappa_0 = C_0 / L0.
# Source: Ch. 6, text below Eq. (90), printed p. 203: 1 <= C_0 <= 8 pi.
# Ch. 12, Eq. (50), printed p. 502, uses kappa_0(h) = 1 / L0(h), so C_0 = 1.
C0_DEFAULT = 1.0


def _outer_scale_bracket(a, kappa_w0):
    '''
    Return the outer-scale (or jitter) reduction factor of the wander integrand.

    formula:
        1 - [ (kappa W0 a)^2 / (1 + (kappa W0 a)^2) ]^(1/6)
    with a = |Theta0 + Theta0_bar xi| the refractive beam-radius ratio W(z)/W0.
    Source: Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196, Ch. 6,
    Eq. (93), printed p. 203. The a^2 in the numerator is not easy to read in
    the printed equation. Ch. 6, Eq. (98), printed p. 204, confirms it: the
    focused-beam (Theta0 = 0, a = xi) reduction of Eq. (93) gives the printed
    (8/3)(kappa_0 W0)^(1/3) INT xi^2 (1 + kappa_0^2 W0^2 xi^2)^(-1/6) dxi only
    when the numerator carries a^2.

    A zero kappa is the infinite outer scale. Then the factor is 1.
    '''
    if kappa_w0 <= 0.0:
        return 1.0
    x = (kappa_w0 * a) ** 2
    return 1.0 - (x / (1.0 + x)) ** (1.0 / 6.0)


def _wander_integral(theta0, kappa_w0):
    '''
    Return the dimensionless xi-integral of Ch. 6, Eq. (93), printed p. 203.

    formula:
        INT_0^1 xi^2 |Theta0 + Theta0_bar xi|^(-1/3)
                * { 1 - [ (kappa W0 a)^2 / (1 + (kappa W0 a)^2) ]^(1/6) } dxi
    The value is 1/3 for a collimated beam (Theta0 = 1) with an infinite outer
    scale, and 3/8 for a beam focused in the receiver plane (Theta0 = 0).

    A convergent beam that focuses INSIDE the path (Theta0 < 0) makes the
    refractive radius pass through zero at xi = -Theta0/Theta0_bar. The
    singularity is integrable, so the integral is split there.
    '''
    theta0 = float(theta0)
    theta_bar0 = 1.0 - theta0

    def f(xi):
        a = abs(theta0 + theta_bar0 * xi)
        if a <= 0.0:
            return 0.0
        return xi ** 2 * a ** (-1.0 / 3.0) * _outer_scale_bracket(a, kappa_w0)

    points = None
    if theta_bar0 != 0.0:
        xi_zero = -theta0 / theta_bar0
        if 0.0 < xi_zero < 1.0:
            points = [xi_zero]
    value, _ = quad(f, 0.0, 1.0, points=points, limit=200)
    return value


def _shape_factor(theta0, kappa_w0, spectrum):
    '''
    Return the dimensionless shape factor of the wander variance.

    The wander variance is WANDER_CONSTANT * Cn2 * L^3 * W0^(-1/3) times this
    factor. For the Kolmogorov branch the factor is the closed form of Ch. 6,
    Eq. (94), printed p. 204, divided by 3:
        (1/3) 2F1(1/3, 1; 4; 1 - |Theta0|)
    For the exponential branch it is the xi-integral of Eq. (93).

    The two routes AGREE to machine precision for Theta0 >= 0. They DIFFER for
    Theta0 < 0, because Eq. (94) takes the absolute value of Theta0 AFTER the
    integration, not inside it. Worked Example 4, Ch. 6, printed p. 216, settles
    which one the book means: a beam focused at 900 m over a 1 km path
    (Theta0 = -0.1111) gives 1.90 cm, and only the Eq. (94) route reproduces it.
    So the Kolmogorov branch always uses Eq. (94).
    '''
    if spectrum == 'kolmogorov':
        return np.real(hyp2f1(1.0 / 3.0, 1.0, 4.0, 1.0 - abs(theta0))) / 3.0
    return _wander_integral(theta0, kappa_w0)


def _spectrum_cutoff(spectrum, L0, c0):
    '''
    Return the spectral cut-off kappa [rad/m] of the wander filter.

    "kolmogorov" is the plain power law with an infinite outer scale, so
    kappa = 0 and L0 must stay None. "exponential" is the outer-scale spectrum
    of Ch. 3, Sec. 3.3.2, restated as Ch. 6, Eq. (90), printed p. 203, with
    kappa_0 = C_0 / L0.
    '''
    if spectrum == 'kolmogorov':
        if L0 is not None:
            raise ValueError(
                'spectrum="kolmogorov" has an infinite outer scale. Pass '
                'spectrum="exponential" to use L0.')
        return 0.0
    if spectrum == 'exponential':
        if L0 is None:
            raise ValueError('spectrum="exponential" needs an outer scale L0.')
        return c0 / np.asarray(L0, dtype=float)
    raise ValueError(
        f'unknown spectrum {spectrum!r}. Use "kolmogorov" or "exponential".')


def beam_wander_variance(w0, wavelength, z, cn2, *, f0=np.inf,
                         spectrum='kolmogorov', L0=None, c0=C0_DEFAULT):
    '''
    Return the RADIAL beam-wander displacement variance <r_c^2> [m^2].

    The path is homogeneous: Cn2 is constant over it.

    Parameters:
        w0 : float
            Beam RADIUS (1/e field) at the transmitter [m].
        wavelength : float
            Optical wavelength [m]. It sets the beam parameters only. The
            geometrical-optics form of Eq. (93) does not use it, so the
            Kolmogorov answer is independent of it.
        z : float
            Path length L [m].
        cn2 : float
            Refractive-index structure constant [m^-2/3].
        f0 : float
            Phase-front radius of curvature at the transmitter [m]. Use
            numpy.inf for a collimated beam, a negative value for a divergent
            beam, and a positive value for a convergent beam.
        spectrum : str
            "kolmogorov" (infinite outer scale) or "exponential" (finite outer
            scale, Ch. 6, Eq. (90)).
        L0 : float, optional
            Outer scale [m]. Only for spectrum="exponential".
        c0 : float
            Outer-scale scaling constant, kappa_0 = c0 / L0. The book allows
            1 <= C_0 <= 8 pi (Ch. 6, text below Eq. (90), printed p. 203).

    Returns:
        float
            <r_c^2> [m^2]. RADIAL (two-dimensional). Divide by 2 for a per-axis
            variance. See the module CONVENTION note.

    formula:
        <r_c^2> = 7.25 Cn2 L^3 W0^(-1/3)
                  INT_0^1 xi^2 |Theta0 + Theta0_bar xi|^(-1/3)
                       { 1 - [ (kappa_0 W(z)/W0... )^2 / (1 + ...) ]^(1/6) } dxi
        Kolmogorov limit:
        <r_c^2> = 2.42 Cn2 L^3 W0^(-1/3) 2F1(1/3, 1; 4; 1 - |Theta0|)
    Source: Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196, Ch. 6,
    Eq. (93), printed p. 203 (general), Eq. (94), printed p. 204 (infinite outer
    scale), Eq. (95), printed p. 204 (collimated), Eq. (96), printed p. 204
    (focused), Eq. (97), printed p. 204 (collimated with an outer scale).

    Limits: weak irradiance fluctuations, and an outer scale larger than the
    transmitter beam. Ch. 6, Fig. 6.8, printed p. 205, shows that a finite outer
    scale near the beam size almost removes the wander.
    '''
    kappa = _spectrum_cutoff(spectrum, L0, c0)
    bp = beam_params(w0, wavelength, z, f0)
    shape = _shape_factor(float(bp.theta0), float(kappa) * float(w0), spectrum)
    return (WANDER_CONSTANT * float(cn2) * float(z) ** 3
            * float(w0) ** (-1.0 / 3.0) * shape)


def spherical_fried_parameter(wavelength, z, cn2):
    '''
    Return Fried's parameter r0 [m] of a reciprocal point source over the path.

    formula:
        r0 = (0.16 Cn2 k^2 L)^(-3/5)
    Source: Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196, Ch. 8,
    text below Eq. (33), printed p. 272. The same constant appears as the
    worked value 0.16 in Ch. 9, printed p. 384. It is the SPHERICAL-wave
    (point-source) coherence width, because the reciprocity argument of
    Sec. 8.3 sends a point source from the receiver back to the transmitter.

    The pointing-error filter of Ch. 8, Eq. (35), printed p. 272, uses this r0.
    '''
    k = wavenumber(wavelength)
    return (0.16 * float(cn2) * k ** 2 * float(z)) ** (-3.0 / 5.0)


def pointing_error_variance(w0, wavelength, z, cn2, *, f0=np.inf,
                            c_r=CR_DEFAULT, r0=None):
    '''
    Return the RADIAL wander-induced pointing-error variance sigma_pe^2 [m^2].

    Beam wander has two parts. The hot spot dances, driven by every eddy up to
    the outer scale. And the WHOLE short-term beam moves around its unperturbed
    position, driven by eddies up to Fried's parameter r0 only. The second part
    is the beam jitter. It flattens the mean beam profile near the boresight,
    which acts as an effective pointing error sigma_pe. That pointing error
    raises the on-axis scintillation index of an UNTRACKED beam.

    The integral is the beam-wander integral of Ch. 6, Eq. (93), with the outer
    scale kappa_0 replaced by the jitter cut-off kappa_r = C_r / r0.

    Parameters:
        w0, wavelength, z, cn2, f0 : as beam_wander_variance.
        c_r : float
            Scaling constant of the jitter cut-off. The book says C_r ~ 2 pi
            (Ch. 8, printed p. 273). Ch. 12, Fig. 12.12, printed p. 505, plots
            kappa_r = 10 and kappa_r = 40 to show the spread.
        r0 : float, optional
            Fried's parameter [m]. Defaults to spherical_fried_parameter.

    Returns:
        float
            sigma_pe^2 [m^2]. RADIAL, the same convention as <r_c^2>.
            Its square root is the `pointing_error_m` argument of
            olb.turbulence.andrews.scintillation.scintillation_index.

    formula:
        sigma_pe^2 = 7.25 Cn2 L^3 W0^(-1/3)
                     INT_0^1 xi^2 |Theta0 + Theta0_bar xi|^(-1/3)
                          { 1 - [ (kappa_r W(z)/W0)^2/(1+...) ]^(1/6) } dxi
        collimated: sigma_pe^2 = 0.48 (lambda L / 2 W0)^2 (2 W0/r0)^(5/3)
                                 [ 1 - (C_r^2 W0^2/r0^2
                                        / (1 + C_r^2 W0^2/r0^2))^(1/6) ]
    Source: Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196, Ch. 8,
    Eqs. (34) and (35), printed p. 272 (the jitter filter), Eq. (36), printed
    p. 273 (the general form), Eq. (37), printed p. 273 (collimated), Eq. (38),
    printed p. 274 (focused).

    Note the asymptotic behaviour of Ch. 8, Eq. (39), printed p. 274:
    sigma_pe^2 falls to zero for BOTH a very small and a very large beam.
    '''
    if r0 is None:
        r0 = spherical_fried_parameter(wavelength, z, cn2)
    kappa_r = float(c_r) / float(r0)
    bp = beam_params(w0, wavelength, z, f0)
    shape = _wander_integral(float(bp.theta0), kappa_r * float(w0))
    return (WANDER_CONSTANT * float(cn2) * float(z) ** 3
            * float(w0) ** (-1.0 / 3.0) * shape)


def long_term_beam_radius(beam, sigma2_R):
    '''
    Return the long-term beam radius W_LT [m] under weak fluctuations.

    formula:
        W_LT = W sqrt(1 + 1.33 sigma_R^2 Lambda^(5/6))
    Source: Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196, Ch. 6,
    Eq. (86), printed p. 202 (restated from Eq. (46)).

    Parameters:
        beam : BeamParams
            The free-space beam parameters at the receiver, from
            olb.turbulence.andrews.beam.beam_params.
        sigma2_R : float or numpy.ndarray
            The plane-wave Rytov variance sigma_R^2 over the same path.

    This is the WEAK-fluctuation long-term radius. For the strong-fluctuation
    radius use effective_beam_params in beam.py (Ch. 7, Eq. (57)).
    '''
    sigma2_R = np.asarray(sigma2_R, dtype=float)
    return beam.w * np.sqrt(1.0 + 1.33 * sigma2_R * beam.lam ** (5.0 / 6.0))


def short_term_beam_radius(beam, sigma2_R, rc2):
    '''
    Return the short-term beam radius W_ST [m].

    formula:
        W_LT^2 = W_ST^2 + <r_c^2>,   so   W_ST = sqrt(W_LT^2 - <r_c^2>)
    Source: Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196, Ch. 6,
    Eq. (100), printed p. 205 (Fante's relation), restated as Ch. 8, Eq. (32),
    printed p. 271.

    THE FACTOR ON <r_c^2> IS 1, NOT 2. That is only consistent with a RADIAL
    (two-dimensional) <r_c^2>, which is the convention of this module. A code
    that adds 2 <beta^2> reads <beta^2> as a PER-AXIS variance. See Conflict
    C-03 in docs/andrews-crosscheck.md.

    Parameters:
        beam : BeamParams
            The free-space beam parameters at the receiver.
        sigma2_R : float or numpy.ndarray
            The plane-wave Rytov variance over the path.
        rc2 : float or numpy.ndarray
            The RADIAL beam-wander variance <r_c^2> [m^2].

    The result is clipped at zero. A negative argument means the wander variance
    is larger than the long-term spot, which breaks the weak-fluctuation limit
    of Eq. (86).
    '''
    w_lt = long_term_beam_radius(beam, sigma2_R)
    return np.sqrt(np.maximum(w_lt ** 2 - np.asarray(rc2, dtype=float), 0.0))


def _slant_prefactor_and_xi(hs, range_m, elevation_deg):
    '''
    Return (airmass, xi) for a slant path sampled on the altitude grid hs.

    The path coordinate is z = h * airmass, with airmass = sec(zenith), and
    xi = 1 - z/L. So xi = 1 at the ground transmitter and falls towards the
    receiver. Source: Andrews and Phillips, 2nd ed. (2005),
    DOI 10.1117/3.626196, Ch. 12, text below Eq. (55), printed p. 503:
    xi = 1 - (h - h0)/(H - h0).
    '''
    hs = np.asarray(hs, dtype=float)
    airmass = 1.0 / np.sin(np.radians(float(elevation_deg)))
    xi = 1.0 - hs * airmass / float(range_m)
    return airmass, np.clip(xi, 0.0, None)


def _slant_variance(w0, wavelength, hs, cn2_profile, range_m, f0,
                    elevation_deg, kappa):
    '''
    Return the slant-path wander variance for a given spectral cut-off kappa.

    formula:
        <r_c^2> = 7.25 (H-h0)^2 sec^3(zenith) W0^(-1/3)
                  INT_h0^H Cn2(h) xi^2 |Theta0 + Theta0_bar xi|^(-1/3)
                       { 1 - [ (kappa W0 a)^2 / (1 + (kappa W0 a)^2) ]^(1/6) }
                  dh
    Source: Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196, Ch. 12,
    Eq. (50), printed p. 502 (with kappa = kappa_0(h), the outer scale), and
    Eq. (53), printed p. 503 (with kappa = kappa_r, the jitter cut-off).

    The code writes the same integral with the slant range L and the airmass,
    because (H-h0)^2 sec^3(z) = L^2 sec(z) and dz = sec(z) dh. This matches the
    way the rest of olb passes a zenith Cn2(h) profile plus an elevation.

    kappa may be a scalar or an array on the hs grid, so a height-dependent
    outer scale L0(h) works.
    '''
    hs = np.asarray(hs, dtype=float)
    cn2 = np.asarray(cn2_profile, dtype=float)
    airmass, xi = _slant_prefactor_and_xi(hs, range_m, elevation_deg)
    bp = beam_params(w0, wavelength, range_m, f0)
    theta0 = float(bp.theta0)
    a = np.abs(theta0 + (1.0 - theta0) * xi)

    kappa = np.asarray(kappa, dtype=float)
    x = (kappa * float(w0) * a) ** 2
    bracket = np.where(kappa > 0.0, 1.0 - (x / (1.0 + x)) ** (1.0 / 6.0), 1.0)

    shape = np.where(a > 0.0, xi ** 2 * a ** (-1.0 / 3.0) * bracket, 0.0)
    integral = np.trapz(cn2 * shape, hs)
    return (WANDER_CONSTANT * float(range_m) ** 2 * airmass
            * float(w0) ** (-1.0 / 3.0) * integral)


def beam_wander_variance_slant(w0, wavelength, hs, cn2_profile, range_m, *,
                               f0=np.inf, elevation_deg=90.0,
                               spectrum='kolmogorov', L0=None, c0=C0_DEFAULT):
    '''
    Return the RADIAL beam-wander variance <r_c^2> [m^2] over a slant path.

    Parameters:
        w0 : float
            Beam radius at the ground transmitter [m].
        wavelength : float
            Optical wavelength [m].
        hs : numpy.ndarray
            Heights above the ground station [m], ascending.
        cn2_profile : numpy.ndarray
            ZENITH Cn2(h) profile on the hs grid [m^-2/3].
        range_m : float
            Slant range L from the transmitter to the receiver [m].
        f0 : float
            Phase-front radius of curvature at the transmitter [m].
        elevation_deg : float
            Elevation angle [deg]. 90 is zenith.
        spectrum, L0, c0 : as beam_wander_variance. L0 may be an array on hs,
            so a height-dependent outer scale works.

    Returns:
        float
            <r_c^2> [m^2], RADIAL.

    Source: Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196, Ch. 12,
    Eq. (50), printed p. 502. Ch. 12, Eq. (51), printed p. 502, gives the
    nominal-outer-scale algebraic form, which this module does not need.

    On an uplink the turbulence sits in the first 20 km of a path that is
    hundreds of km long. So xi is close to 1 over the whole turbulent layer, and
    the answer is NOT the homogeneous-path reduction 2.42 Cn2 L^3 W0^(-1/3).
    '''
    kappa = _spectrum_cutoff(spectrum, L0, c0)
    return _slant_variance(w0, wavelength, hs, cn2_profile, range_m, f0,
                           elevation_deg, kappa)


def plane_fried_parameter_slant(wavelength, hs, cn2_profile,
                                elevation_deg=90.0):
    '''
    Return Fried's parameter r0 [m] for a slant path through a Cn2(h) profile.

    formula:
        r0 = [ 0.42 sec(zenith) k^2 INT Cn2(h) dh ]^(-3/5)
    Source: Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196, Ch. 12,
    Eq. (23), printed p. 492. Ch. 12, Eq. (53), printed p. 503, names this r0
    for the uplink pointing-error cut-off.
    '''
    hs = np.asarray(hs, dtype=float)
    cn2 = np.asarray(cn2_profile, dtype=float)
    k = wavenumber(wavelength)
    airmass = 1.0 / np.sin(np.radians(float(elevation_deg)))
    mu0 = np.trapz(cn2, hs)
    return (0.42 * airmass * k ** 2 * mu0) ** (-3.0 / 5.0)


def pointing_error_variance_slant(w0, wavelength, hs, cn2_profile, range_m, *,
                                  f0=np.inf, elevation_deg=90.0,
                                  c_r=CR_DEFAULT, r0=None):
    '''
    Return the RADIAL pointing-error variance sigma_pe^2 [m^2] over a slant path.

    Parameters: as beam_wander_variance_slant, plus
        c_r : float
            Scaling constant of the jitter cut-off kappa_r = c_r / r0.
        r0 : float, optional
            Fried's parameter [m]. Defaults to plane_fried_parameter_slant.

    Returns:
        float
            sigma_pe^2 [m^2], RADIAL. Its square root is the
            `pointing_error_m` argument of
            olb.turbulence.andrews.scintillation.scintillation_index.

    Source: Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196, Ch. 12,
    Eq. (53), printed p. 503. It is Eq. (50) with kappa_0 replaced by kappa_r.

    Ch. 12, Eqs. (54) and (57), printed pp. 503 and 504, then put sigma_pe and
    sqrt(<r_c^2>) into the untracked and tracked uplink scintillation index.
    Those two equations belong to the scintillation module, not to this one.
    '''
    if r0 is None:
        r0 = plane_fried_parameter_slant(wavelength, hs, cn2_profile,
                                         elevation_deg)
    kappa_r = float(c_r) / float(r0)
    return _slant_variance(w0, wavelength, hs, cn2_profile, range_m, f0,
                           elevation_deg, kappa_r)


if __name__ == '__main__':
    # ================= part 1: physics self-checks =================
    lam_m = 1550e-9

    # Ch. 6, Eq. (95), printed p. 204: a collimated beam with an infinite outer
    # scale gives exactly 2.42 Cn2 L^3 W0^(-1/3).
    rc2 = beam_wander_variance(0.05, lam_m, 2000.0, 3e-16)
    closed = WANDER_CONSTANT_COLLIMATED * 3e-16 * 2000.0 ** 3 * 0.05 ** (-1 / 3)
    assert abs(rc2 / closed - 1.0) < 1e-12, (rc2, closed)
    print(f'Eq. (95) collimated       : <rc^2> = {rc2:.6e} m^2  '
          f'rms = {np.sqrt(rc2) * 1e3:.3f} mm')

    # Ch. 6, Eq. (96), printed p. 204: a beam focused in the RECEIVER plane
    # (f0 = L, so Theta0 = 0) gives 2.72 Cn2 L^3 W0^(-1/3).
    rc2_foc = beam_wander_variance(0.05, lam_m, 2000.0, 3e-16, f0=2000.0)
    const_foc = rc2_foc / (3e-16 * 2000.0 ** 3 * 0.05 ** (-1 / 3))
    assert abs(const_foc - 2.72) < 0.01, const_foc
    print(f'Eq. (96) focused constant : {const_foc:.4f}  (book 2.72)')

    # The wavelength does not enter the geometrical-optics Kolmogorov form.
    assert np.isclose(beam_wander_variance(0.05, 800e-9, 2000.0, 3e-16), rc2)

    # The exponential branch reduces to the Kolmogorov branch as L0 -> infinity,
    # for Theta0 >= 0. This is the collimated-limit check at 1e-9. The outer
    # scale must be very large, because the reduction factor goes as
    # (kappa_0 W0)^(1/3), which converges slowly.
    rc2_big_L0 = beam_wander_variance(0.05, lam_m, 2000.0, 3e-16,
                                      spectrum='exponential', L0=1e40)
    err_limit = abs(rc2_big_L0 / rc2 - 1.0)
    assert err_limit < 1e-9, err_limit
    print(f'LIMIT L0 -> inf           : |ratio - 1| = {err_limit:.3e} '
          f'(target 1e-9)')

    # Ch. 6, Eq. (97), printed p. 204, and Fig. 6.8, printed p. 205: at
    # kappa_0 W0 = 0.1 the wander is less than 70 % of the infinite-outer-scale
    # value.
    ratio_fig68 = beam_wander_variance(0.05, lam_m, 2000.0, 3e-16,
                                       spectrum='exponential',
                                       L0=0.05 / 0.1) / rc2
    assert ratio_fig68 < 0.70, ratio_fig68
    print(f'Fig. 6.8 kappa0 W0 = 0.1  : ratio = {ratio_fig68:.4f}  (book < 0.70)')

    # Ch. 8, Eq. (37), printed p. 273: the collimated closed form of sigma_pe^2.
    cn2_t, L_t, w0_t = 3e-16, 2000.0, 0.05
    pe2 = pointing_error_variance(w0_t, lam_m, L_t, cn2_t)
    r0_t = spherical_fried_parameter(lam_m, L_t, cn2_t)
    x_t = (CR_DEFAULT * w0_t / r0_t) ** 2
    pe2_book = (0.48 * (lam_m * L_t / (2 * w0_t)) ** 2
                * (2 * w0_t / r0_t) ** (5.0 / 3.0)
                * (1.0 - (x_t / (1.0 + x_t)) ** (1.0 / 6.0)))
    err_pe = abs(pe2 / pe2_book - 1.0)
    assert err_pe < 0.02, (pe2, pe2_book)
    print(f'Eq. (37) collimated pe    : integral {pe2:.4e} vs closed '
          f'{pe2_book:.4e} m^2, ratio {pe2 / pe2_book:.5f}')

    # The pointing error is a filtered part of the wander, so it is smaller.
    assert pe2 < rc2, (pe2, rc2)

    # ---------------- REDUCTION check: Worked Example 2 ----------------
    # Ch. 6, printed p. 215. Collimated, L = 1 km, W0 = 1 cm, lambda = 1.55 um,
    # Cn2 = 1e-13, l0 = 0, L0 = infinity. The book prints
    #   W = 5 cm, W_LT = 6.52 cm, sqrt(<rc^2>) = 3.35 cm, W_ST = 5.59 cm.
    from .beam import beam_params as _bp
    L_we2, w0_we2, cn2_we2 = 1000.0, 0.01, 1e-13
    bp_we2 = _bp(w0_we2, lam_m, L_we2)
    sigma2_R_we2 = 1.23 * cn2_we2 * wavenumber(lam_m) ** (7 / 6) * L_we2 ** (11 / 6)
    rc2_we2 = beam_wander_variance(w0_we2, lam_m, L_we2, cn2_we2)
    w_lt_we2 = long_term_beam_radius(bp_we2, sigma2_R_we2)
    w_st_we2 = short_term_beam_radius(bp_we2, sigma2_R_we2, rc2_we2)
    book = {'sigma_R^2': 1.99, 'W': 0.05, 'W_LT': 0.0652,
            'rms rc': 0.0335, 'W_ST': 0.0559}
    got = {'sigma_R^2': sigma2_R_we2, 'W': float(bp_we2.w),
           'W_LT': float(w_lt_we2), 'rms rc': float(np.sqrt(rc2_we2)),
           'W_ST': float(w_st_we2)}
    for name in book:
        rel = abs(got[name] / book[name] - 1.0)
        assert rel < 0.02, (name, got[name], book[name], rel)
        print(f'WORKED EXAMPLE 2 {name:>9} : {got[name]:.5g} vs book '
              f'{book[name]:.5g}  ({100 * rel:.2f} %, target 2 %)')

    # ---------------- REDUCTION check: Worked Example 4 ----------------
    # Ch. 6, printed p. 216. lambda = 1.55 um, 2W0 = 10 cm, L = 1 km,
    # Cn2 = 5e-14. Collimated -> 1.81 cm. Convergent, focused at 900 m -> 1.90 cm.
    rms_coll = np.sqrt(beam_wander_variance(0.05, lam_m, 1000.0, 5e-14))
    rms_conv = np.sqrt(beam_wander_variance(0.05, lam_m, 1000.0, 5e-14, f0=900.0))
    assert abs(rms_coll / 0.0181 - 1.0) < 0.02, rms_coll
    assert abs(rms_conv / 0.0190 - 1.0) < 0.02, rms_conv
    assert rms_conv > rms_coll
    print(f'WORKED EXAMPLE 4          : collimated {rms_coll * 100:.3f} cm '
          f'(book 1.81), convergent {rms_conv * 100:.3f} cm (book 1.90)')

    # ---------------- slant path ----------------
    from ..._deps import get_c2n

    hs_grid = np.logspace(0.0, np.log10(20e3), 20)
    cn2_hv = get_c2n(hs_grid, 21.0, 1.7e-14)
    L_up, w0_up = 600e3, 1.0
    rc2_up = beam_wander_variance_slant(w0_up, lam_m, hs_grid, cn2_hv, L_up)
    # A lower elevation puts more air in the path, so the wander grows.
    rc2_up30 = beam_wander_variance_slant(w0_up, lam_m, hs_grid, cn2_hv, L_up,
                                          elevation_deg=30.0)
    assert rc2_up30 > rc2_up, (rc2_up30, rc2_up)
    # A finite outer scale can only reduce the wander (Ch. 12, Fig. 12.11).
    rc2_up_L0 = beam_wander_variance_slant(w0_up, lam_m, hs_grid, cn2_hv, L_up,
                                           spectrum='exponential', L0=10.0)
    assert rc2_up_L0 < rc2_up, (rc2_up_L0, rc2_up)
    pe2_up = pointing_error_variance_slant(w0_up, lam_m, hs_grid, cn2_hv, L_up)
    assert pe2_up < rc2_up, (pe2_up, rc2_up)
    r0_up = plane_fried_parameter_slant(lam_m, hs_grid, cn2_hv)
    print(f'slant zenith 600 km       : rms wander {np.sqrt(rc2_up):.2f} m '
          f'({np.sqrt(rc2_up) / L_up * 1e6:.2f} urad), rms pe '
          f'{np.sqrt(pe2_up):.2f} m, r0 = {r0_up * 100:.1f} cm')

    # ================= part 2: C-01 / C-03 adjudication =================
    # Andrews against the shared Dios kernel. This module CHANGES NOTHING in the
    # kernel or in the olb Dios path. It only measures. See Conflicts C-01 and
    # C-03 in docs/andrews-crosscheck.md.
    from ..._deps import (beam_wander_variance as kernel_wander,
                          short_term_beam_waist, long_term_beam_waist,
                          spherical_wave_coherence_diameter, gaussz)

    print('')
    print('C-01/C-03 adjudication (Andrews wander.py against coupled_flux.py)')
    print('  case          quantity                              value')

    # --- terrestrial case: 1550 nm, 2 km, Cn2 3e-16, W0 5 cm, collimated ---
    z_t = np.linspace(0.0, L_t, 2001)
    cn2_arr_t = cn2_t * np.ones_like(z_t)
    ws_free_t = gaussz(w0_t, z_t, lam_m)            # free-space W(z), diffracting
    ws_gom_t = w0_t * np.ones_like(z_t)             # refractive W(z) of Eq. (93)
    k_free_t = kernel_wander(L_t, cn2_arr_t, ws_free_t, z_t)
    k_gom_t = kernel_wander(L_t, cn2_arr_t, ws_gom_t, z_t)
    red_t = WANDER_CONSTANT_COLLIMATED * cn2_t * L_t ** 3 * w0_t ** (-1 / 3)

    # --- uplink case: the defaults of uplink_flux.py's own self-check ---
    airmass_up = 1.0
    ws_free_up = gaussz(w0_up, hs_grid, lam_m)
    k_up = kernel_wander(L_up, cn2_hv * airmass_up, ws_free_up, hs_grid)
    mu0_up = np.trapz(cn2_hv, hs_grid)
    red_up = WANDER_CONSTANT_COLLIMATED * (mu0_up / L_up) * L_up ** 3 * w0_up ** (-1 / 3)

    rows = [
        ('terrestrial', 'Andrews Eq. (93)/(94)  <rc^2> [m^2]', rc2),
        ('terrestrial', 'kernel 2.07, free-space W(z) [m^2]', k_free_t),
        ('terrestrial', 'kernel 2.07, GOM W(z)=W0    [m^2]', k_gom_t),
        ('terrestrial', 'Eq. (94) reduced 2.42 form  [m^2]', red_t),
        ('terrestrial', 'ratio Andrews / kernel(free-space)', rc2 / k_free_t),
        ('terrestrial', 'ratio Andrews / kernel(GOM)', rc2 / k_gom_t),
        ('terrestrial', 'ratio Eq.(94) / Andrews', red_t / rc2),
        ('terrestrial', 'ratio Eq.(94) / kernel(GOM)', red_t / k_gom_t),
        ('uplink', 'Andrews Eq. (50) slant <rc^2> [m^2]', rc2_up),
        ('uplink', 'kernel 2.07, free-space W(z) [m^2]', k_up),
        ('uplink', 'ratio Andrews / kernel', rc2_up / k_up),
        ('uplink', 'homogeneous 2.42 surrogate  [m^2]', red_up),
        ('uplink', 'ratio surrogate / Andrews (xi=1 effect)', red_up / rc2_up),
    ]
    for case, name, value in rows:
        print(f'  {case:<13} {name:<38} {value:.6g}')

    # The pure constant ratio must be exactly 7.25 / 2.07 = 3.5024.
    assert abs(rc2 / k_gom_t - WANDER_CONSTANT / 2.07) < 1e-6, rc2 / k_gom_t

    # --- C-03: the long-term waist convention ---
    # Andrews Eq. (100): W_LT^2 = W_ST^2 + <rc^2>, factor 1, <rc^2> RADIAL.
    # Kernel long_term_beam_waist: W_LT^2 = W_ST^2 + 2 <beta^2>, factor 2, so
    # <beta^2> must be PER-AXIS. Feed BOTH rules the SAME short-term waist, so
    # the ratio carries the convention and the constant only.
    k_t = wavenumber(lam_m)
    r0s_t = spherical_wave_coherence_diameter(k_t, L_t, cn2_arr_t, z_t)
    w_st_kernel = float(short_term_beam_waist(w0_t, L_t, np.pi * w0_t ** 2 / lam_m,
                                              k_t, r0s_t))
    w_lt_andrews = float(np.sqrt(w_st_kernel ** 2 + rc2))       # factor 1, radial
    w_lt_kernel = float(long_term_beam_waist(w_st_kernel, k_free_t))  # factor 2
    print(f'  {"terrestrial":<13} {"W_ST (kernel, shared input) [m]":<38} '
          f'{w_st_kernel:.6g}')
    print(f'  {"terrestrial":<13} {"W_LT Andrews Eq.(100) f=1 radial [m]":<38} '
          f'{w_lt_andrews:.6g}')
    print(f'  {"terrestrial":<13} {"W_LT kernel f=2 per-axis     [m]":<38} '
          f'{w_lt_kernel:.6g}')
    print(f'  {"terrestrial":<13} {"ratio W_LT Andrews / kernel":<38} '
          f'{w_lt_andrews / w_lt_kernel:.6g}')
    # A per-axis reading of the kernel would need 2 * 2.07 = 4.14 against 7.25,
    # so it does NOT close the gap either.
    print(f'  {"both":<13} {"kernel factor-2 vs Andrews factor-1":<38} '
          f'{WANDER_CONSTANT / (2.0 * 2.07):.6g}')

    print('self-check passed')
