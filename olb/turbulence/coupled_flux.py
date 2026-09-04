'''
Dios coupled-flux kernels for the LEO uplink (beam wander + scintillation).

This module holds the lower-level coupled-flux physics that
olb.turbulence.uplink_flux composes into a short uplink Monte Carlo. It was
borrowed from the my_analysis_modules kernel (coupled_flux.py) and is vendored
here VERBATIM, so olb owns it and does not depend on that repository. The
formulas keep their citations.

The physics is Dios et al., Applied Optics 43(18), 3866 (2004),
DOI 10.1364/AO.43.003866, with the beam-wander constant 2.07 from Belmonte,
Applied Optics 39(27), 5426 (2000), DOI 10.1364/AO.39.005426. The beam wander
KEEPS the Dios/Belmonte 2.07, NOT the Andrews 7.25: the whole uplink chain is
internally consistent Dios, and two split-step simulations validate 2.07. See
Conflict C-01 in docs/andrews-crosscheck.md, and olb.turbulence.andrews.wander
for the independent Andrews measurement kept side by side.

The module needs numpy, scipy, and the olb assumptions decorator layer. It
imports nothing else from the rest of olb and nothing from my_analysis_modules.

Each public physics function declares its own validity through the `@assumes`
decorator (olb.assumptions). Three assumptions recur here: the coherence-diameter
path weight is transmitter-referred and must not flip (`PATH_WEIGHT`); the
beam-wander variance is a radial (two-axis) variance (`RADIAL_VARIANCE`); and the
beam-wander constant is the Dios/Belmonte 2.07, which conflicts with the Andrews
7.25 (`C01_WANDER`, Conflict C-01). Outside a collection context the decorator is
a no-op, so the numeric output does not change.
'''

import numpy as np
from scipy.special import gamma, hyp1f1

from ..assumptions import (assumes, Constraint, BEAM_GAUSSIAN,
                           BEAM_SPHERICAL_WAVE, REGIME_WEAK, SPECTRUM_KOLMOGOROV)
from .andrews.beam import beam_params
from .andrews.scintillation import (rytov_variance, rytov_weak,
                                    WEAK_REGIME_LIMIT)

# ----------------------------------------------------------------------------
# The module assumptions, as shared Constraint instances (see the docstring).
# ----------------------------------------------------------------------------

# The transmitter-referred path weight of the spherical-wave coherence diameter.
PATH_WEIGHT = Constraint(
    "path-weight",
    "The coherence-diameter path weight ((L-z)/L)^(5/3) is transmitter-referred, "
    "for the uplink. Do not flip it to the receiver-referred (z/L)^(5/3).",
    "10.1364/AO.43.003866", "Eq. (3), printed p. 3868")

# The beam-wander variance is a radial (two-axis) quantity, never per-axis.
RADIAL_VARIANCE = Constraint(
    "variance-convention",
    "The beam-wander variance <beta^2> is the radial (two-axis) displacement "
    "variance. A one-axis draw uses 0.5*<beta^2>.",
    "10.1364/AO.43.003866", "Eqs. (9) and (10), printed p. 3868")

# The Dios/Belmonte 2.07 versus the Andrews 7.25 (Conflict C-01).
C01_WANDER = Constraint(
    "conflict",
    "The beam-wander constant is the Dios/Belmonte 2.07, NOT the Andrews 7.25 "
    "(a factor 3.50 lower). Two split-step simulations validate 2.07. See "
    "Conflict C-01.",
    "10.1364/AO.39.005426",
    "Belmonte Eq. (21), printed p. 5435; against Andrews DOI 10.1117/3.626196, "
    "Ch. 6, Eq. (93), printed p. 203")

# The lognormal flux draw of the coupled-flux sample.
LOGNORMAL_DRAW = Constraint(
    "pdf-shape",
    "The flux fluctuation is drawn from a lognormal irradiance PDF (a Gaussian "
    "log-amplitude).",
    "10.1364/AO.43.003866", "Eqs. (25) to (27), printed p. 3870")

# The wander-removal correction of the short-term waist turns negative for
# r0s/W0 > (1/0.26)^3 ~ 57. That is a small aperture in very weak turbulence.
_CORRECTION_MAX_RATIO = (1.0 / 0.26) ** 3


def _correction_range_check(args, result):
    '''Return a reason when the short-term-waist correction leaves its range.

    The correction 1 - 0.26 (r0s/W0)^(1/3) turns negative for
    r0s/W0 > (1/0.26)^3. The check reads the bound arguments and never warns.
    '''
    r0s = float(np.max(np.asarray(args["r0s"], dtype=float)))
    w0 = float(np.min(np.asarray(args["W0"], dtype=float)))
    ratio = r0s / w0
    if ratio > _CORRECTION_MAX_RATIO:
        return (f"r0s/W0 = {ratio:.1f} > {_CORRECTION_MAX_RATIO:.0f}; the "
                "wander-removal correction 1-0.26(r0s/W0)^(1/3) is out of range "
                "(a small aperture in very weak turbulence).")
    return None


CORRECTION_RANGE = Constraint(
    "approximation",
    "The wander-removal correction 1-0.26(r0s/W0)^(1/3) is an aperture-size "
    "effect; it stays in range only for r0s/W0 < 57, and W0 must be the physical "
    "launch radius.",
    "10.1364/AO.43.003866", "Eqs. (4) to (6), printed p. 3868",
    check=_correction_range_check)


def _wander_weak_regime_check(args, result):
    '''Return a reason when the beam-wander path leaves the weak regime.

    The Dios beam-wander variance is a WEAK-fluctuation model. The gate needs
    the wavelength, which the kernel signature does not carry, so the check runs
    only when the caller gives the optional `wavelength`. Without it the check
    returns None and the call is unchanged.

    The gate is the shared, beam-aware `rytov_weak`. It reads BOTH weak
    conditions of a Gaussian beam, sigma_R^2 < 1 AND sigma_R^2 Lambda^(5/6) < 1
    (Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196, Ch. 5, Eq. (16),
    printed p. 140), so a focused beam trips a gate that a plane-wave test
    passes. Only the "hard" tier is a violation. The "soft" tier (past
    RYTOV_CONFIDENT_WEAK but inside the book limit) stays a factory warning,
    because a Constraint check must not warn.

    The plane-wave Rytov variance is
        sigma_R^2 = 1.23 Cn2 k^(7/6) L^(11/6)                (constant Cn2)
        sigma_R^2 = 2.25 k^(7/6) INT Cn2(z) (L-z)^(5/6) dz   (a profile)
    Source: Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196, Ch. 8,
    Eq. (10), printed p. 262; the path-weighted form is Ch. 12, Eq. (38),
    printed p. 495, with the weight measured from the RECEIVER at z = L. The two
    agree for a constant Cn2 (2.25 * 6/11 = 1.227). The check never warns and
    never raises.
    '''
    wavelength = args["wavelength"]
    if wavelength is None:
        return None
    L = float(args["L"])
    z_grid = np.asarray(args["z"], dtype=float)
    cn2_vals = np.asarray(args["cn2"], dtype=float)
    w0 = float(np.asarray(args["ws"], dtype=float).reshape(-1)[0])
    flat = cn2_vals.reshape(-1)
    if flat.size < 2 or np.all(flat == flat[0]):
        # The constant-Cn2 shortcut, the same call the terrestrial factory made.
        sigma2_R = float(rytov_variance(wavelength, L, float(flat[0])))
    else:
        k = 2.0 * np.pi / float(wavelength)
        weight = np.maximum(L - z_grid, 0.0) ** (5.0 / 6.0)
        sigma2_R = float(2.25 * k ** (7.0 / 6.0)
                         * np.trapezoid(cn2_vals * weight, z_grid))
    Lambda = float(beam_params(w0, wavelength, L).lam)
    if rytov_weak(sigma2_R, Lambda) == 'hard':
        return (f"sigma2_R={sigma2_R:.3f} (Lambda={Lambda:.3f}) meets or exceeds "
                f"the Gaussian-beam weak limit {WEAK_REGIME_LIMIT}; the "
                "beam-wander model is not trusted. Use the fidelity-2 Monte "
                "Carlo.")
    return None


WANDER_WEAK_REGIME = Constraint(
    "regime",
    "The beam-wander variance is a weak-fluctuation model. It holds while "
    "sigma_R^2 < 1 AND sigma_R^2 Lambda^(5/6) < 1. Give the optional "
    "`wavelength` and the check runs; without it the gate is the caller's duty.",
    "10.1117/3.626196",
    "Ch. 5, Eq. (16), printed p. 140; Ch. 8, Eq. (10), printed p. 262",
    check=_wander_weak_regime_check)


def _lambda_function(L, k0, wL):
    return 2 * L / (k0 * wL ** 2)


def _theta_function(L, Z0):
    return (1 + (L / Z0) ** 2) ** (-1)


def _A(z, L, k_0, wL):
    '''Equation (17) of Dios et al. 2004, DOI 10.1364/AO.43.003866.'''
    Lambda = _lambda_function(L, k_0, wL)
    return (Lambda * L / k_0) * ((L - z) / L) ** 2


def _B(z, L, k_0, Z0):
    '''Equation (18) of Dios et al. 2004, DOI 10.1364/AO.43.003866.'''
    Theta = _theta_function(L, Z0)
    return (L / k_0) * ((L - z) / L) * (Theta + (1 - Theta) * z / L)


@assumes(PATH_WEIGHT, beam_type=BEAM_SPHERICAL_WAVE, spectrum=SPECTRUM_KOLMOGOROV)
def spherical_wave_coherence_diameter(k, L, cn2_vals, z):
    '''
    Spherical-wave coherence diameter r0s.

    r0s = [ 0.42*k^2 * integral_0^L( Cn2(z) * ((L-z)/L)^(5/3) dz ) ]^(-3/5)

    This is equation (3) of Dios et al. 2004, DOI 10.1364/AO.43.003866 (printed
    page 3868). The paper states it "for the uplink", so the weight
    ((L-z)/L)^(5/3) is TRANSMITTER-referred. Andrews and Phillips 2nd ed.,
    chapter 6, equations (115) and (116), print the mirror weight (z/L)^(5/3)
    for the receiver-referred (downlink) radius. The two are a plane-of-
    reference difference, not a fault. Do not flip this weight.

    Parameters:
        k (float) : Optical wave number (2*pi/lambda).
        L (float) : Propagation distance [m].
        cn2_vals (numpy.ndarray) : Cn2 values sampled at `z`.
        z (numpy.ndarray) : Path coordinates for `cn2_vals` [m].

    Returns:
        float : Spherical-wave coherence diameter r0s.
    '''
    z = np.asarray(z)
    if L < np.max(z):
        z_grid = z[z <= L]
        cn2_vals = cn2_vals[z <= L]
    else:
        z_grid = z

    weight = ((L - z_grid) / L) ** (5.0 / 3.0)
    path_integral = np.trapezoid(cn2_vals * weight, z_grid)

    return (0.42 * (k ** 2) * path_integral) ** (-3 / 5)


@assumes(CORRECTION_RANGE, beam_type=BEAM_GAUSSIAN, turbulence_regime=REGIME_WEAK,
         spectrum=SPECTRUM_KOLMOGOROV)
def short_term_beam_waist(W0, L, Z0, k, r0s, return_squared=False, w_free=None):
    '''
    Short-term beam waist at propagation distance z=L.

    W_ST^2(L) = W0^2 * (1 + L^2/Z0^2) + 2*{ [4.2*L/(k*r0s)] * [1 - 0.26*(r0s/W0)^(1/3)] }^2

    The W0^2*(1 + L^2/Z0^2) term is the free-space width at L for a COLLIMATED
    beam. For a deliberately diverged (or focused) transmitter, pass the actual
    free-space width as `w_free` and it replaces that term; `Z0` is then unused.
    The turbulence term is unaffected -- it is a spreading angle set by r0s, not
    by the transmitter geometry -- but `W0` must stay the physical beam radius
    at the launch aperture, since the 0.26*(r0s/W0)^(1/3) correction (which
    removes the wander contribution) is an aperture-size effect.

    Note that correction turns negative for r0s/W0 > (1/0.26)^3 ~ 57, i.e. a
    small aperture in very weak turbulence. It is squared, so the result stays
    positive, and the 4.2*L/(k*r0s) prefactor shrinks faster than the correction
    grows (the term falls off as r0s^(-4/3)), so the waist still tends to the
    free-space value -- but the factor itself is outside its intended range there.

    Source: Dios et al. 2004, DOI 10.1364/AO.43.003866, equations (4) to (6).

    Parameters:
        W0 (float) : Beam radius at the transmit aperture [m].
        L (float) : Propagation distance [m].
        Z0 (float) : Rayleigh range [m]. Ignored when `w_free` is given.
        k (float) : Optical wave number (2*pi/lambda).
        r0s (float) : Spherical-wave coherence diameter [m].
        return_squared (bool, optional) : If True, return W_ST^2 instead of W_ST.
        w_free (float, optional) : Free-space (turbulence-free) beam radius at
            L [m], for a non-collimated transmitter.

    Returns:
        float : Short-term beam waist (or its square).
    '''
    if w_free is None:
        geometric_term = W0 ** 2 * (1.0 + (L ** 2) / (Z0 ** 2))
    else:
        geometric_term = w_free ** 2
    turbulence_factor = (4.2 * L / (k * r0s)) * (1.0 - 0.26 * (r0s / W0) ** (1 / 3))
    w_st_squared = geometric_term + 2.0 * turbulence_factor ** 2

    if return_squared:
        return w_st_squared
    return np.sqrt(np.maximum(w_st_squared, 0.0))


@assumes(RADIAL_VARIANCE, C01_WANDER, WANDER_WEAK_REGIME, beam_type=BEAM_GAUSSIAN,
         turbulence_regime=REGIME_WEAK, spectrum=SPECTRUM_KOLMOGOROV)
def beam_wander_variance(L, cn2, ws, z, *, wavelength=None):
    '''
    Beam-wander variance <beta^2>.

    <beta^2> = 2.07 * integral_0^L( Cn2(z) * (L-z)^2 * [1/Ws(z)]^(1/3) dz )

    This is equation (11) of Dios et al. 2004, DOI 10.1364/AO.43.003866
    (printed page 3868). The constant, the integrand and the path weight agree
    with the paper. Dios does not derive equation (11). He takes it from
    Belmonte, Applied Optics 39, 5426 (2000), DOI 10.1364/AO.39.005426.

    CONVENTION: <beta^2> is the RADIAL (two-axis) displacement variance. Dios
    equation (9) gives beta = sqrt(beta_x^2 + beta_y^2), and equation (10) gives
    <beta_x^2> = <beta_y^2> = 0.5*<beta^2>. So a caller that draws one Cartesian
    axis must use the variance 0.5*<beta^2>.

    KNOWN DIFFERENCE from Andrews and Phillips, 2nd ed., DOI 10.1117/3.626196.
    Chapter 6, equation (93) (printed page 203) has the SAME integrand with the
    constant 7.25, and its infinite-outer-scale form, equation (94) (printed
    page 204), gives <r_c^2> = 2.42 Cn2 L^3 W0^(-1/3) for a collimated beam.
    For a constant Cn2 and Ws = W0, the constant 2.07 above gives
    0.69 Cn2 L^3 W0^(-1/3), which is 3.50 times lower. The Andrews quantity is
    also a radial variance, so the axis convention does NOT explain the gap.

    KEEP 2.07. Two split-step wave-optics simulations validate this form. Dios
    figure 3 (printed page 3871) compares equation (11) with an FFT-BPM
    simulation of the same uplink, and the two agree closely. Belmonte 2000
    (DOI 10.1364/AO.39.005426, the source Dios takes equation (11) from) prints
    the same 2.07 form as his equation (21) (printed page 5435) and compares it
    with his own phase-screen simulation in figures 11 and 12; it matches in
    weak-to-moderate turbulence. A factor of 3.50 would be plain on either plot.
    See the C-01 closure in docs/andrews-crosscheck.md.

    Ws(z) is the beam radius at z. Dios prints the symbol W_s(z) but does not
    define it. This function uses the free-space (diffracted) radius. For an
    uplink that choice changes almost nothing, because all the turbulence is in
    the first 20 km, where Ws(z) stays near W0.

    Parameters:
        L (float) : Propagation distance [m].
        cn2 (numpy.ndarray) : Cn2 profile sampled at `z`.
        ws (numpy.ndarray) : Beam radius profile Ws(z) [m], sampled at `z`.
        z (numpy.ndarray) : Path coordinates [m].
        wavelength (float, optional) : Optical wavelength [m]. It does NOT
            change the result. It turns ON the weak-regime runtime check
            (WANDER_WEAK_REGIME), which needs the Rytov variance. None (the
            default) leaves the check off, so every old call is unchanged.

    Returns:
        float : Beam-wander variance <beta^2> [m^2].
    '''
    z_grid = np.asarray(z)
    cn2_vals = np.asarray(cn2, dtype=float)
    ws_vals = np.asarray(ws, dtype=float)
    integrand = cn2_vals * ((L - z_grid) ** 2) * ((1.0 / ws_vals) ** (1 / 3))
    return 2.07 * np.trapezoid(integrand, z_grid)


@assumes(RADIAL_VARIANCE, beam_type=BEAM_GAUSSIAN, turbulence_regime=REGIME_WEAK)
def long_term_beam_waist(w_st, beta2):
    '''
    Long-term beam waist, combining short-term spreading with beam wander.

    W_LT^2(L) = W_ST^2(L) + 2*<beta^2>

    This is equation (1) of Dios et al. 2004, DOI 10.1364/AO.43.003866 (printed
    page 3867). The paper repeats it in equation (29), printed page 3870. The
    factor 2 is the paper's own factor on a RADIAL <beta^2> (see
    `beam_wander_variance`). It is NOT a per-axis to radial conversion.

    Andrews and Phillips, 2nd ed., DOI 10.1117/3.626196, chapter 6, equation
    (100) (printed page 205), puts the factor 1 on a radial <r_c^2>. With the
    constant of each source, the wander part of W_LT^2 is
    1.38 Cn2 L^3 W0^(-1/3) by Dios and 2.42 Cn2 L^3 W0^(-1/3) by Andrews. So
    the two combination rules differ by 1.75, not by 3.50. The Dios factor 2 and
    the Dios constant 2.07 partially cancel the gap.

    Parameters:
        w_st (float) : Short-term beam waist [m].
        beta2 (float) : Beam-wander variance <beta^2> [m^2].

    Returns:
        float : Long-term beam waist [m].
    '''
    return np.sqrt(w_st ** 2 + 2 * beta2)


@assumes(beam_type=BEAM_GAUSSIAN, turbulence_regime=REGIME_WEAK,
         spectrum=SPECTRUM_KOLMOGOROV)
def off_axis_scintillation_index(L, k_0, wL, cn2s, z_points, r):
    '''
    Off-axis scintillation index sigma_r,L^2(r, L), equation (20) of Dios et al.
    2004, DOI 10.1364/AO.43.003866.

    Parameters:
        L (float) : Propagation distance [m].
        k_0 (float) : Optical wave number (2*pi/lambda).
        wL (float) : Beam waist at distance L [m].
        cn2s (numpy.ndarray) : Cn2 values sampled at `z_points`.
        z_points (numpy.ndarray) : Altitude points [m] for `cn2s`.
        r (float) : Off-axis distance at which to evaluate [m].

    Returns:
        float : Off-axis scintillation index at radius `r`.
    '''
    a_z = _A(z_points, L, k_0, wL)
    hyp_arg = 2 * r ** 2 / wL ** 2
    integrand = cn2s * a_z ** (5 / 6) * (hyp1f1(-5 / 6, 1, hyp_arg) - 1)
    result = np.trapezoid(integrand, z_points)

    coefficient = 4 * np.pi ** 2 * k_0 ** 2 * gamma(-5 / 6) * 0.033
    return coefficient * result


@assumes(beam_type=BEAM_GAUSSIAN, turbulence_regime=REGIME_WEAK,
         spectrum=SPECTRUM_KOLMOGOROV)
def on_axis_scintillation_index(L, k_0, wL, Z0, cn2s, z_points):
    '''
    On-axis scintillation index sigma_r^2(0, L), equation (16) of Dios et al.
    2004, DOI 10.1364/AO.43.003866.

    Parameters:
        L (float) : Propagation distance [m].
        k_0 (float) : Optical wave number (2*pi/lambda).
        wL (float) : Beam waist at distance L [m].
        Z0 (float) : Rayleigh range [m].
        cn2s (numpy.ndarray) : Cn2 values sampled at `z_points`.
        z_points (numpy.ndarray) : Altitude points [m] for `cn2s`.

    Returns:
        float : On-axis scintillation index.
    '''
    a_z = _A(z_points, L, k_0, wL)
    b_z = _B(z_points, L, k_0, Z0)
    ratio = b_z / a_z

    # The cosine multiplies ONLY the second term. Dios et al. 2004 equation (16)
    # (DOI 10.1364/AO.43.003866) and Andrews and Phillips 2nd ed. chapter 8,
    # equation (17) (printed page 263, DOI 10.1117/3.626196) both give
    # A^(5/6) - (A^2 + B^2)^(5/12) * cos[(5/6) arctan(B/A)]. Factor out A^(5/6)
    # to get the form below. Before 2026-08 a parenthesis closed too early, so
    # the cosine multiplied the full bracket.
    integrand = cn2s * a_z ** (5 / 6) * (1 - (1 + ratio ** 2) ** (5 / 12) * np.cos((5 / 6) * np.arctan(ratio)))
    result = np.trapezoid(integrand, z_points)

    coefficient = 4 * np.pi ** 2 * k_0 ** 2 * gamma(-5 / 6) * 0.033
    return coefficient * result


def mean_off_axis_irradiance(r, wlt_L):
    '''
    Mean off-axis irradiance profile (long-term beam spread only).

    Parameters:
        r (float or numpy.ndarray) : Off-axis distance [m].
        wlt_L (float) : Long-term beam waist at the receiver [m].

    Returns:
        float or numpy.ndarray : Normalized mean off-axis irradiance.
    '''
    return np.exp(-2 * r ** 2 / wlt_L ** 2)


@assumes(LOGNORMAL_DRAW, beam_type=BEAM_GAUSSIAN, turbulence_regime=REGIME_WEAK,
         spectrum=SPECTRUM_KOLMOGOROV)
def coupled_flux_sample(beta, cn2_profile, Z0, hs, L, k_0, wL, wL_lt):
    '''
    Draw a single realization of the turbulence-induced flux fluctuation at
    beam-wander displacement `beta`.

    Parameters:
        beta (float) : Instantaneous beam-wander displacement [m].
        cn2_profile (numpy.ndarray) : Cn2 profile sampled at `hs`.
        Z0 (float) : Rayleigh range [m].
        hs (numpy.ndarray) : Altitude points [m] for `cn2_profile`.
        L (float) : Propagation distance [m].
        k_0 (float) : Optical wave number (2*pi/lambda).
        wL (float) : Free-space (diffraction-limited) beam radius at L [m].
            This is the W(L) of Dios equation (15), which sets Lambda.
        wL_lt (float) : Long-term beam waist at distance L [m]. Equation (24)
            uses it for the mean-irradiance weight of equation (25).

    Returns:
        tuple : (xi, xi_on_axis, sigma2_x, sigma2_x_on_axis, sigma2_gauss,
        sigma2_gauss_on_axis)
    '''
    sigma2_off = off_axis_scintillation_index(L, k_0, wL, cn2_profile, hs, beta)
    sigma2_on = on_axis_scintillation_index(L, k_0, wL, Z0, cn2_profile, hs)

    # Dios et al. 2004 equation (25) (DOI 10.1364/AO.43.003866, printed page
    # 3870): sigma2_I,Gb = (sigma2_I + sigma2_I,r) * <I>^2, where <I> is the
    # mean irradiance of equation (24) at the wander position beta. Equations
    # (13), (16) and (20) normalize the index to the LOCAL mean irradiance;
    # equation (25) re-normalizes it to the mean irradiance at the BEAM CENTER,
    # which is the normalization that equation (26) needs. Section 5, step
    # (c)(ii) of the paper tells you to use equation (25) at this point.
    #
    # An earlier patch removed this weight, because Andrews and Phillips 2nd ed.
    # chapter 8, equations (9) and (15) (DOI 10.1117/3.626196) keep the local
    # normalization. But this module implements Dios, and the removal made it
    # disagree with the paper it cites. The weight is back (2026-08-25).
    I_off = mean_off_axis_irradiance(beta, wL_lt)
    sigma2_gauss = (sigma2_on + sigma2_off) * I_off ** 2
    sigma2_gauss_on_axis = sigma2_on * I_off ** 2

    sigma2_x = 0.25 * np.log(1 + sigma2_gauss)
    sigma2_x_on_axis = 0.25 * np.log(1 + sigma2_gauss_on_axis)

    xi = np.random.normal(-sigma2_x, np.sqrt(sigma2_x), 1)
    xi_on_axis = np.random.normal(-sigma2_x_on_axis, np.sqrt(sigma2_x_on_axis), 1)

    return xi, xi_on_axis, sigma2_x, sigma2_x_on_axis, sigma2_gauss, sigma2_gauss_on_axis


def on_axis_irradiance(beta, wst_L, xi_beta):
    '''
    Instantaneous on-axis irradiance, combining beam-wander displacement `beta`
    with a log-amplitude turbulence fluctuation `xi_beta`.

    Parameters:
        beta (float or numpy.ndarray) : Instantaneous beam-wander displacement [m].
        wst_L (float) : Short-term beam waist at the receiver [m].
        xi_beta (float or numpy.ndarray) : Log-amplitude fluctuation sample.

    Returns:
        float or numpy.ndarray : Normalized on-axis irradiance.
    '''
    return np.exp(2 * xi_beta) * np.exp(-2 * beta ** 2 / wst_L ** 2)


if __name__ == '__main__':
    # A light sanity check of the vendored kernels. The end-to-end numeric
    # cross-check against the my_analysis_modules original is in the commit that
    # added this module (a seeded _flux_result run matched bit-for-bit).
    hs = np.logspace(0, np.log10(20e3), 20)
    cn2 = 0.00594 * (21 / 27) ** 2 * (1e-5 * hs) ** 10 * np.exp(-hs / 1000) \
        + 2.7e-16 * np.exp(-hs / 1500) + 1.7e-14 * np.exp(-hs / 100)
    k = 2 * np.pi / 1550e-9
    L = 1075e3
    cn2_slant = cn2 / np.sin(np.radians(30.0))

    r0s = spherical_wave_coherence_diameter(k, L, cn2_slant, hs)
    assert 0.05 < r0s < 0.5, r0s
    w_st = short_term_beam_waist(0.05, L, np.pi * 0.05 ** 2 / 1550e-9, k, r0s,
                                 w_free=0.6)
    assert w_st > 0.6                                 # turbulence broadens it
    beta2 = beam_wander_variance(L, cn2_slant, hs * 0 + 0.6, hs)
    assert beta2 > 0
    w_lt = long_term_beam_waist(w_st, beta2)
    assert w_lt > w_st                                # wander adds spread

    np.random.seed(0)
    xi, _, s2x, _, s2g, _ = coupled_flux_sample(
        np.array([0.01]), cn2_slant, np.pi * 0.05 ** 2 / 1550e-9, hs, L, k,
        0.6, 0.7)
    assert np.isfinite(xi) and s2x > 0 and s2g > 0
    Is = on_axis_irradiance(np.array([0.0]), w_st, np.array([0.0]))
    assert abs(float(Is[0]) - 1.0) < 1e-12            # on axis, no fluctuation

    print(f"coupled_flux self-check: r0s={r0s * 100:.2f} cm, "
          f"w_st={w_st:.2f} m, beta2={beta2:.3f} m^2")

    # ---------------- assumption self-checks ----------------
    import warnings
    from ..assumptions import trace_assumptions

    # (1) Value parity: one representative call returns the identical float with
    #     and without a collection context.
    r0s_out = spherical_wave_coherence_diameter(k, L, cn2_slant, hs)
    with trace_assumptions():
        r0s_in = spherical_wave_coherence_diameter(k, L, cn2_slant, hs)
    assert r0s_out == r0s_in, (r0s_out, r0s_in)

    # (2) Registration: inside a context the expected sources and kinds register.
    with trace_assumptions() as tr:
        spherical_wave_coherence_diameter(k, L, cn2_slant, hs)
        beam_wander_variance(L, cn2_slant, hs * 0 + 0.6, hs)
        long_term_beam_waist(0.6, beta2)
    mod = __name__
    assert f"{mod}.spherical_wave_coherence_diameter" in tr.records
    assert f"{mod}.beam_wander_variance" in tr.records
    kinds = {c.kind for rec in tr.records.values() for c in rec.constraints}
    assert {"path-weight", "variance-convention", "conflict"} <= kinds, kinds
    # The transmitter-referred path weight and the C-01 conflict both register.
    wander_rec = tr.records[f"{mod}.beam_wander_variance"]
    assert any(c.kind == "conflict" for c in wander_rec.constraints)

    # (3) A deliberately out-of-range call yields a source-prefixed violation,
    #     and the decorator check itself emits NO warning (a separate matter from
    #     any pre-existing warnings.warn, of which this function has none).
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        with trace_assumptions() as tr_bad:
            # r0s/W0 = 1.0/0.01 = 100 > 57 -> the correction leaves its range.
            short_term_beam_waist(0.01, L, np.pi * 0.01 ** 2 / 1550e-9, k, 1.0)
    assert any(v.startswith(f"[{mod}.short_term_beam_waist]")
               for v in tr_bad.violations), tr_bad.violations
    assert len(caught) == 0, "a decorator check must not warn"

    # (4) The function-owned weak-regime gate of the beam wander (2026-09-04).
    #     The optional `wavelength` turns the check on. It does NOT change the
    #     value. A horizontal 3 km path with a 2 cm launch waist is the
    #     terrestrial case: weak at Cn2 = 1e-16, hard at Cn2 = 1e-13.
    lam_h = 1550e-9
    L_h = 3e3
    zs_h = np.linspace(0.0, L_h, 64)
    ws_h = np.full_like(zs_h, 0.02)

    plain = beam_wander_variance(L_h, np.full_like(zs_h, 1e-13), ws_h, zs_h)
    with_lam = beam_wander_variance(L_h, np.full_like(zs_h, 1e-13), ws_h, zs_h,
                                    wavelength=lam_h)
    assert plain == with_lam, (plain, with_lam)     # the value does not move

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        with trace_assumptions() as tr_weak:
            beam_wander_variance(L_h, np.full_like(zs_h, 1e-16), ws_h, zs_h,
                                 wavelength=lam_h)
        with trace_assumptions() as tr_hard:
            beam_wander_variance(L_h, np.full_like(zs_h, 1e-13), ws_h, zs_h,
                                 wavelength=lam_h)
        with trace_assumptions() as tr_off:
            beam_wander_variance(L_h, np.full_like(zs_h, 1e-13), ws_h, zs_h)
    assert not any("beam-wander model is not trusted" in v
                   for v in tr_weak.violations), tr_weak.violations
    assert any(v.startswith(f"[{mod}.beam_wander_variance]")
               and "beam-wander model is not trusted" in v
               for v in tr_hard.violations), tr_hard.violations
    assert not tr_off.violations, "no wavelength -> no regime check"
    assert len(caught) == 0, "the regime check must not warn"

    # A Cn2 PROFILE takes the path-weighted form. It agrees with the
    # constant-Cn2 shortcut to the book rounding (2.25*6/11 = 1.227 vs 1.23).
    with trace_assumptions() as tr_prof:
        beam_wander_variance(L_h, np.full(zs_h.shape, 1e-13) * (1.0 + 1e-12 * zs_h),
                             ws_h, zs_h, wavelength=lam_h)
    assert any("beam-wander model is not trusted" in v
               for v in tr_prof.violations), tr_prof.violations

    print("coupled_flux assumptions self-check passed")
    print("coupled_flux self-check passed")
