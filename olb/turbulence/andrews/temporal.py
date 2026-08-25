'''
Temporal statistics of the irradiance, from Andrews and Phillips.

This module puts a TIME axis on the turbulence physics. It gives the temporal
power spectrum of the irradiance, the quasi-frequency that the fade-rate
equations need, and the Greenwood frequency of an adaptive-optics servo.

Source of every equation:
    L. C. Andrews and R. L. Phillips, "Laser Beam Propagation through Random
    Media", 2nd ed., SPIE Press (2005). DOI: 10.1117/3.626196
Each function names its section, its equation number, and its printed page.

THE BRIDGE FROM SPACE TO TIME. The Taylor frozen-turbulence hypothesis moves an
eddy pattern across the path at the transverse wind speed V. It does not let the
pattern change shape. So a spatial separation rho becomes a time lag tau through
rho = V tau (Ch. 3, Eq. (27), printed p. 73). The book gives the limit of that
hypothesis on the same page: it fails when V is much less than the turbulent
wind fluctuations. That occurs when the mean wind blows along the line of
sight.

FREQUENCY CONVENTION. The book writes each spectrum against the ANGULAR
frequency omega [rad/s], and it normalises with Ch. 8, Eq. (55), printed p. 282:
    sigma_I^2 = (1/2pi) INTEGRAL(0..inf) S_I(omega) domega.
This module takes and returns the CYCLIC frequency f [Hz], with omega = 2 pi f.
The two forms hold the same numbers, because
    sigma_I^2 = INTEGRAL(0..inf) S_I(2 pi f) df.
So the spectrum of this module integrates DIRECTLY to the scintillation index,
and its unit is seconds (power per hertz on a unit-mean irradiance).

WHAT IS BUILT:

- `taylor_wavenumber` maps a temporal frequency to a spatial wavenumber.
- `fresnel_frequency` gives the transition frequency f_t of the spectrum.
- `irradiance_temporal_spectrum` gives S_I(f). The weak branch (Ch. 8.5) covers
  a plane wave, a spherical wave, and a Gaussian beam. The strong branch
  (Ch. 9.8 and Ch. 10.3.6) covers a plane wave, with an optional receiver
  aperture that averages the spectrum.
- `quasi_frequency` gives nu0 from the second moment of a spectrum.
- `greenwood_frequency` and `coherence_time` give the servo bandwidth of an
  adaptive-optics system on a slant path.

This module holds physics only. It returns no decibels.

WHAT IS NOT BUILT, AND WHY:

- A finite inner scale or outer scale in ANY temporal spectrum. Ch. 9.8, printed
  p. 364, states "We will also ignore the effects of a finite inner scale and
  outer scale". Ch. 10, printed p. 425, adds that a scale changes the peak
  values but does not shift the peak position. The book gives no closed temporal
  form with a scale, so the code refuses `l0` and `L0`.
- A strong-regime spherical wave or Gaussian beam. Ch. 9.8, printed p. 364,
  states "limit our analysis to the case of a plane wave". Ch. 9, printed
  p. 364, adds that no Gaussian-beam covariance has been computed.
- An off-axis (radial) weak spectrum. Ch. 8, Eqs. (66) and (67), printed
  pp. 286-287, give it. This module gives the longitudinal (on-axis) part only.
'''

import numpy as np
from scipy.integrate import quad
from scipy.special import gamma as _gamma_fn, hyp1f1, kve

from ..._deps import v_wind
from .aperture import d_param
from .beam import wavenumber
from .scintillation import rytov_variance

# Greenwood constant of Ch. 14, Eq. (38), printed p. 622, restated as Ch. 14,
# Eq. (98), printed p. 637. DOI: 10.1117/3.626196
GREENWOOD_CONSTANT = 2.91

# Greenwood time constant for a constant wind, tau0 = 0.32 r0 / V. Source:
# Ch. 14, Eq. (39), printed p. 622. DOI: 10.1117/3.626196
GREENWOOD_R0_CONSTANT = 0.32

# Amplitude of the longitudinal temporal spectrum. Source: Ch. 8, Eq. (65),
# printed p. 285. DOI: 10.1117/3.626196
_PSD_AMPLITUDE = 3.90

# Tail-matched amplitude of the second spectral group of Ch. 8, Eq. (65).
#
#     C = -Gamma(-1/3) Gamma(11/6) / [ Gamma(1/2) Gamma(7/3) ] = 1.810729
#
# WHY THE CODE DERIVES THIS CONSTANT INSTEAD OF PRINTING IT. Ch. 8, Eq. (65)
# writes the second group as 0.29 i^(4/3) a_j^(-4/3). The parameter a_j of
# Ch. 8, Eq. (64) is COMPLEX, so a_j^(-4/3) has a branch cut, and the printed
# form does not say which branch to take. The code removes that ambiguity: it
# writes every argument through q_j = 1/(4 i a_j), which is 1/2 for a plane wave
# and 2/9 for a spherical wave, exactly as Ch. 8, Eqs. (57) and (59) print them.
# The coefficient then becomes c_j = C q_j^(4/3), which gives 0.7186 for a plane
# wave and 0.2437 for a spherical wave. The book prints 0.72 and 0.24, so C
# reproduces both printed constants to the two figures the book gives.
#
# C follows from the book's own statement that the spectrum decays as
# omega^(-8/3) above the Fresnel frequency (Ch. 8, text below Eq. (57), printed
# p. 283). The two spectral groups each carry a 1/omega tail. The tails must
# cancel, and C is the value that cancels them. With the book's rounded 0.72 the
# cancellation is incomplete, and the residual 1/omega tail makes the spectral
# integral diverge slowly. See the self-check, which measures both.
# DOI: 10.1117/3.626196
_PSD_TAIL_CONSTANT = float(-_gamma_fn(-1.0 / 3.0) * _gamma_fn(11.0 / 6.0)
                           / (_gamma_fn(0.5) * _gamma_fn(7.0 / 3.0)))

# Switch point of the confluent hypergeometric function. Below it the code uses
# scipy. Above it scipy loses every significant figure on a large imaginary
# argument, so the code uses the asymptotic expansion. Source of the expansion:
# NIST DLMF 13.7.2 (https://dlmf.nist.gov/13.7.E2). A comparison with mpmath on
# the imaginary axis sets the switch point: scipy holds 4e-11 at |z| = 20 and
# loses to 4e-7 at |z| = 30, while the 20-term expansion holds 1e-10 at |z| = 20
# and improves above it. So |z| = 22 keeps the whole range inside about 1e-10.
_ASYMPTOTIC_Z = 22.0
_ASYMPTOTIC_TERMS = 20

_WAVES = ('plane', 'spherical', 'gaussian')
_REGIMES = ('weak', 'strong')


def _hyp1f1(a, b, z):
    '''
    Return the confluent hypergeometric function 1F1(a; b; z) of a complex z.

    scipy.special.hyp1f1 loses all precision when |z| passes about 40 on the
    imaginary axis. The temporal spectra need |z| up to 1e10. So this helper
    uses scipy below |z| = 22 and the large-argument expansion of NIST
    DLMF 13.7.2 above it:
        1F1(a;b;z) = Gamma(b) [ e^z z^(a-b)/Gamma(a) S1
                              + (-z)^(-a)/Gamma(b-a) S2 ]
        S1 = SUM_n (b-a)_n (1-a)_n / (n! z^n)
        S2 = SUM_n (a)_n (a-b+1)_n / (n! (-z)^n)
    A check against mpmath sets the switch point. Both paths hold about 1e-10
    at |z| = 22, which is the worst point of the whole range.
    '''
    z = np.asarray(z, dtype=complex)
    out = np.empty_like(z)
    small = np.abs(z) < _ASYMPTOTIC_Z
    if np.any(small):
        out[small] = hyp1f1(a, b, z[small])
    big = ~small
    if np.any(big):
        zb = z[big]
        s1 = np.zeros_like(zb)
        term = np.ones_like(zb)
        for n in range(_ASYMPTOTIC_TERMS):
            if n:
                term = term * ((b - a + n - 1) * (1.0 - a + n - 1)) / (n * zb)
            s1 = s1 + term
        s2 = np.zeros_like(zb)
        term = np.ones_like(zb)
        for n in range(_ASYMPTOTIC_TERMS):
            if n:
                term = term * ((a + n - 1) * (a - b + 1 + n - 1)) / (n * -zb)
            s2 = s2 + term
        out[big] = _gamma_fn(b) * (
            np.exp(zb) * zb ** (a - b) / _gamma_fn(a) * s1
            + (-zb) ** (-a) / _gamma_fn(b - a) * s2)
    return out


def taylor_wavenumber(freq, wind_speed):
    '''
    Map a temporal frequency to a spatial wavenumber by frozen flow.

    formula:
        kappa = 2 pi f / V
    Source: Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196, Ch. 3,
    Sec. 3.4, Eq. (27), printed p. 73. The hypothesis moves the eddy pattern
    across the path at the speed V and does not change its shape, so a spatial
    period 2 pi/kappa passes the receiver in a time 1/f.

    The book states the limit of the hypothesis on the same page: it fails when
    V is much less than the turbulent wind fluctuations. That occurs when the
    mean wind blows along the line of sight.

    Parameters:
        freq : float or numpy.ndarray
            Temporal frequency [Hz].
        wind_speed : float or numpy.ndarray
            Mean wind speed transverse to the path [m/s].

    Returns:
        float or numpy.ndarray
            Spatial wavenumber kappa [rad/m].
    '''
    return (2.0 * np.pi * np.asarray(freq, dtype=float)
            / np.asarray(wind_speed, dtype=float))


def fresnel_frequency(wind_speed, wavelength, z):
    '''
    Return the Fresnel (transition) frequency of the irradiance spectrum [Hz].

    formula:
        omega_t = V / sqrt(L/k),   f_t = omega_t / (2 pi),   k = 2 pi/lambda
    Source: Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196, Ch. 8,
    text below Eq. (57), printed p. 283. The book states there that the spectrum
    stays flat below omega_t and decays as omega^(-8/3) above it. Ch. 9, text
    below Eq. (128), printed p. 365, repeats the same reading.

    The frequency is the wind speed divided by the Fresnel scale sqrt(L/k). The
    Fresnel scale sets the correlation width of the irradiance in weak
    fluctuations.

    Parameters:
        wind_speed : float or numpy.ndarray
            Mean transverse wind speed V [m/s].
        wavelength : float
            Optical wavelength [m].
        z : float
            Path length L [m].

    Returns:
        float or numpy.ndarray
            The Fresnel frequency f_t [Hz].
    '''
    return _angular_fresnel(wind_speed, wavelength, z) / (2.0 * np.pi)


def _angular_fresnel(wind_speed, wavelength, z):
    '''Return omega_t = V / sqrt(L/k) [rad/s]. See `fresnel_frequency`.'''
    k = wavenumber(wavelength)
    return (np.asarray(wind_speed, dtype=float)
            / np.sqrt(np.asarray(z, dtype=float) / k))


def _weak_beam_group(wave, beam):
    '''
    Return the two spectral group parameters (d_t, q1, q2) of Ch. 8, Eq. (64).

    formula:
        d_t = 0.67 - 0.17 Theta
        a_1 = 1 / (4 i d_t [1 - (Theta_bar + i Lambda) d_t])
        a_2 = 1 / (4 Lambda d_t^2)
        q_j = 1 / (4 i a_j)
            q_1 = d_t [1 - (Theta_bar + i Lambda) d_t],   q_2 = -i Lambda d_t^2
    Source: Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196, Ch. 8,
    Eq. (64), printed p. 285.

    The q form removes the branch cut of a_j^(-4/3). It also makes the two
    limits of the book explicit:
        plane wave      Theta = 1, Theta_bar = 0, Lambda = 0  ->  d_t = 0.50,
                        q_1 = 0.50, q_2 = 0. Ch. 8, Eq. (57), printed p. 283,
                        prints the same argument -i omega^2/(2 omega_t^2).
        spherical wave  Theta = 0, Theta_bar = 1, Lambda = 0  ->  d_t = 0.67,
                        q_1 = 0.221, q_2 = 0. Ch. 8, Eq. (59), printed p. 284,
                        prints -2 i omega^2/(9 omega_t^2), so q_1 = 2/9 = 0.222.
    '''
    if wave == 'plane':
        theta, theta_bar, lam = 1.0, 0.0, 0.0
    elif wave == 'spherical':
        theta, theta_bar, lam = 0.0, 1.0, 0.0
    else:
        if beam is None:
            raise ValueError('wave="gaussian" needs beam=BeamParams(...)')
        theta = float(np.asarray(beam.theta))
        theta_bar = float(np.asarray(beam.theta_bar))
        lam = float(np.asarray(beam.lam))
    dt = 0.67 - 0.17 * theta
    q1 = dt * (1.0 - (theta_bar + 1j * lam) * dt)
    q2 = -1j * lam * dt ** 2
    return dt, q1, q2


def _weak_spectrum(u, sigma2_R, omega_t, dt, q1, q2):
    '''
    Return the longitudinal weak-fluctuation spectrum of Ch. 8, Eq. (65).

    formula (with u = omega/omega_t and z_j = -i q_j u^2):
        S(omega) = 3.90 sigma_R^2 / (omega_t d_t^(5/6))
                   Re{ u^(-8/3) [ 1F1(-5/6; -1/3; z_2) - 1F1(-5/6; -1/3; z_1) ]
                       + i^(4/3) [ c_2 1F1(1/2; 7/3; z_2)
                                   - c_1 1F1(1/2; 7/3; z_1) ] }
        c_j = C q_j^(4/3)
    Source: Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196, Ch. 8,
    Eq. (65), printed p. 285. See `_PSD_TAIL_CONSTANT` for C, and
    `_weak_beam_group` for q_j.

    THE CODE MOVES ONE FACTOR. Ch. 8, Eq. (65) prints the factor u^(-8/3) in
    front of BOTH groups. The plane-wave limit Ch. 8, Eq. (57), printed p. 283,
    prints it in front of the FIRST group only, and Ch. 8, Eq. (59) does the
    same for the spherical wave. The code follows Eqs. (57) and (59), because a
    u^(-8/3) on the second group makes that group vanish at high frequency and
    leaves the printed 0.72 and 0.24 constants with nothing to cancel. With the
    reading of this code, Eq. (65) gives Eq. (57) at Theta = 1, Lambda = 0 and
    Eq. (59) at Theta = 0, Lambda = 0. The self-check measures both.
    '''
    u = np.asarray(u, dtype=float)
    z1 = -1j * q1 * u ** 2
    z2 = -1j * q2 * u ** 2
    c1 = _PSD_TAIL_CONSTANT * q1 ** (4.0 / 3.0)
    c2 = _PSD_TAIL_CONSTANT * q2 ** (4.0 / 3.0)
    first = (_hyp1f1(-5.0 / 6.0, -1.0 / 3.0, z2)
             - _hyp1f1(-5.0 / 6.0, -1.0 / 3.0, z1))
    second = (c2 * _hyp1f1(0.5, 7.0 / 3.0, z2)
              - c1 * _hyp1f1(0.5, 7.0 / 3.0, z1))
    with np.errstate(divide='ignore', invalid='ignore'):
        total = u ** (-8.0 / 3.0) * first + (1j ** (4.0 / 3.0)) * second
    amp = _PSD_AMPLITUDE * sigma2_R / (omega_t * dt ** (5.0 / 6.0))
    return amp * np.real(total)


def _strong_eta(sigma2_R):
    '''
    Return the two filter cutoff parameters (eta_X, eta_Y) of a plane wave.

    formula:
        eta_X = 2.61 / (1 + 1.11 sigma_R^(12/5))
        eta_Y = 3 (1 + 0.69 sigma_R^(12/5))
    Source: Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196, Ch. 9,
    Eq. (40), printed p. 335, and Ch. 9, Eq. (45), printed p. 336. Ch. 9,
    Eqs. (119) and (120), printed p. 362, use the same pair in the covariance.
    '''
    s125 = np.asarray(sigma2_R, dtype=float) ** (6.0 / 5.0)
    return 2.61 / (1.0 + 1.11 * s125), 3.0 * (1.0 + 0.69 * s125)


def _strong_covariance(s, sigma2_R, d2):
    '''
    Return the all-regime plane-wave temporal covariance of the irradiance.

    formula (with s = omega_t tau and d^2 = k D_G^2/(4 L)):
        B_lnX(s) = 0.49 sigma_R^2 / (1 + 0.65 d^2 + 1.11 sigma_R^(12/5))^(7/6)
                   1F1(7/6; 1; -s^2 eta_X / (4 + d^2 eta_X))
        B_lnY(s) = 0.51 sigma_R^2 (1 + 0.69 sigma_R^(12/5))^(-5/6)
                   / (1 + 0.90 d^2 + 0.62 d^2 sigma_R^(12/5))
                   (s^2 eta_Y)^(5/12) K_(5/6)(s sqrt(eta_Y))
        B_I(s)   = exp[ B_lnX(s) + B_lnY(s) ] - 1
    Source: Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196, Ch. 10,
    Eqs. (93), (94) and (95), printed pp. 421-422. At d = 0 the pair becomes
    Ch. 9, Eq. (126), printed p. 365. Ch. 9, Eq. (126) prints 0.50 in the
    small-scale amplitude; Ch. 9, Eq. (46), printed p. 336, and Ch. 10, Eq. (95)
    both print 0.51, so the code uses 0.51.

    At s = 0 the covariance becomes the aperture-averaged scintillation index of
    Ch. 10, Eq. (69), printed p. 413. The small-scale limit carries a +0.57 %
    offset, because (x)^(5/12) K_(5/6)(sqrt(x)) goes to 0.5 Gamma(5/6) 2^(5/6) =
    1.0056 as x goes to zero, not to 1. The self-check measures that offset.
    '''
    s = np.asarray(s, dtype=float)
    s125 = sigma2_R ** (6.0 / 5.0)
    eta_x, eta_y = _strong_eta(sigma2_R)
    amp_x = 0.49 * sigma2_R / (1.0 + 0.65 * d2 + 1.11 * s125) ** (7.0 / 6.0)
    amp_y = (0.51 * sigma2_R * (1.0 + 0.69 * s125) ** (-5.0 / 6.0)
             / (1.0 + 0.90 * d2 + 0.62 * d2 * s125))
    b_x = amp_x * hyp1f1(7.0 / 6.0, 1.0,
                         -s ** 2 * eta_x / (4.0 + d2 * eta_x))
    x = s * np.sqrt(eta_y)
    with np.errstate(divide='ignore', invalid='ignore', under='ignore'):
        # K_nu(x) = kve(nu, x) exp(-x). The scaled form keeps a large x finite.
        bessel = np.where(x > 0.0, kve(5.0 / 6.0, x) * np.exp(-x), np.inf)
        b_y = amp_y * np.where(
            x > 0.0, x ** (5.0 / 6.0) * bessel,
            0.5 * _gamma_fn(5.0 / 6.0) * 2.0 ** (5.0 / 6.0))
    return np.exp(b_x + np.nan_to_num(b_y)) - 1.0


def irradiance_temporal_spectrum(freq, wind_speed, wavelength, z, cn2, *,
                                 wave='plane', regime='weak', l0=None,
                                 L0=None, beam=None, D=None):
    '''
    Return the temporal power spectrum of the irradiance S_I(f) [s].

    The spectrum is one-sided and it is normalised on a unit-mean irradiance:
        INTEGRAL(0..inf) S_I(f) df = sigma_I^2.
    Source of that normalisation: Andrews and Phillips, 2nd ed. (2005),
    DOI 10.1117/3.626196, Ch. 8, Eqs. (54) and (55), printed p. 282. See the
    module docstring for the omega-to-f conversion.

    Parameters:
        freq : float or numpy.ndarray
            Temporal frequency f [Hz]. It must be above zero.
        wind_speed : float
            Mean transverse wind speed V [m/s].
        wavelength : float
            Optical wavelength [m].
        z : float
            Path length L [m].
        cn2 : float
            Refractive-index structure parameter [m^-2/3], constant on the path.
        wave : str
            "plane", "spherical", or "gaussian". The strong regime takes
            "plane" only.
        regime : str
            "weak" selects Ch. 8.5. "strong" selects Ch. 9.8 and Ch. 10.3.6,
            which hold in every regime from weak to saturated.
        l0, L0 : float, optional
            Inner scale and outer scale [m]. NOT SUPPORTED. See Raises.
        beam : BeamParams, optional
            Beam parameters at the receiver. Required for wave="gaussian".
        D : float, optional
            Receiver aperture diameter [m]. It averages the spectrum over the
            aperture. The strong regime takes it. The weak regime does not.

    Returns:
        numpy.ndarray
            S_I(f) [s], the same shape as freq.

    Raises:
        ValueError
            If wave or regime is not a known name, or if a Gaussian beam comes
            with no beam parameters.
        NotImplementedError
            If l0 or L0 is set. Ch. 9.8, printed p. 364, states that the book
            ignores a finite inner scale and outer scale in every temporal
            spectrum. Ch. 10, printed p. 425, adds that a scale changes the peak
            values of the spectrum but does not shift the peak position. The
            book gives no closed form, so the code does not guess one.
            If regime="strong" comes with a wave other than "plane". Ch. 9.8,
            printed p. 364, limits the strong-regime analysis to a plane wave.
            If regime="weak" comes with D. Ch. 10.3.6, printed pp. 421-422,
            gives the aperture-averaged temporal covariance in the all-regime
            form only, so use regime="strong" with D.

    formula (weak):
        Ch. 8, Eq. (65), printed p. 285, with Ch. 8, Eqs. (57) and (59),
        printed pp. 283 and 284, as its plane and spherical limits.
    formula (strong):
        Ch. 10, Eqs. (93) to (97), printed pp. 421-422, which reduce at D = 0 to
        Ch. 9, Eqs. (126) to (128), printed p. 365.
    '''
    if wave not in _WAVES:
        raise ValueError(f'unknown wave {wave!r}. Use one of {_WAVES}.')
    if regime not in _REGIMES:
        raise ValueError(f'unknown regime {regime!r}. Use one of {_REGIMES}.')
    if l0 is not None or L0 is not None:
        raise NotImplementedError(
            'a temporal spectrum with a finite inner or outer scale is not in '
            'the book. Ch. 9.8, printed p. 364, states that it ignores both '
            'scales, and Ch. 10, printed p. 425, states only that a scale '
            'changes the peak values. DOI 10.1117/3.626196')

    freq = np.asarray(freq, dtype=float)
    omega_t = float(_angular_fresnel(wind_speed, wavelength, z))
    sigma2_R = float(rytov_variance(wavelength, z, cn2, wave='plane'))

    if regime == 'weak':
        if D is not None:
            raise NotImplementedError(
                'the book gives no weak-only aperture-averaged temporal '
                'spectrum. Ch. 10.3.6, printed pp. 421-422, gives the '
                'all-regime form, so call this with regime="strong" and D. '
                'DOI 10.1117/3.626196')
        dt, q1, q2 = _weak_beam_group(wave, beam)
        u = 2.0 * np.pi * freq / omega_t
        return _weak_spectrum(u, sigma2_R, omega_t, dt, q1, q2)

    if wave != 'plane':
        raise NotImplementedError(
            f'the strong-regime temporal spectrum holds for a plane wave only. '
            f'Ch. 9.8, printed p. 364, states "limit our analysis to the case '
            f'of a plane wave". Got wave={wave!r}. DOI 10.1117/3.626196')
    d2 = 0.0 if D is None else float(d_param(D, wavelength, z)) ** 2

    def cov(tau):
        return float(_strong_covariance(omega_t * tau, sigma2_R, d2))

    # Ch. 9, Eq. (127), printed p. 365, and Ch. 10, Eq. (96), printed p. 422:
    # S_I(omega) = 4 INTEGRAL(0..inf) B_I(tau) cos(omega tau) dtau.
    out = np.array([
        4.0 * quad(cov, 0.0, np.inf, weight='cos',
                   wvar=2.0 * np.pi * float(f), limit=200)[0]
        for f in np.atleast_1d(freq)])
    return out.reshape(freq.shape) if freq.shape else float(out[0])


def quasi_frequency(freq, spectrum):
    '''
    Return the quasi-frequency nu0 [Hz] of a temporal irradiance spectrum.

    formula:
        b_0 = INTEGRAL(0..inf) S_I(f) df
        b_2 = INTEGRAL(0..inf) (2 pi f)^2 S_I(f) df
        nu0 = (1/2 pi) sqrt(b_2 / b_0) = sqrt( INT f^2 S df / INT S df )
    Source: Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196:
        the moments b_0 and b_2   Ch. 11, Eq. (14), printed p. 448
        nu0 = sqrt(b_2/b_0)/(2 pi)  Ch. 11, Eq. (15), printed p. 448
        nu0 from the irradiance covariance
                                  Ch. 11, Eq. (35), printed p. 456, and
                                  Ch. 11, Eq. (38), printed p. 456
        nu0 as the spectral moment ratio
                                  Ch. 12, Eq. (73), printed p. 514
    Ch. 12, Eq. (73) writes exactly this ratio:
        nu0 = (1/2 pi) [ -B_I''(0) / B_I(0) ]^(1/2)
            = (1/2 pi) [ INT omega^2 S_I domega / INT S_I domega ]^(1/2).
    The 2 pi cancels when the moments use the cyclic frequency f, which is what
    this module returns.

    The book states at printed p. 456 that nu0 is roughly the standard deviation
    of the normalised irradiance spectrum read as a probability density, and
    that the spectrum is about 3 nu0 wide.

    nu0 DEPENDS ON THE BAND. The book states at printed p. 283 that the
    irradiance spectrum decays as omega^(-8/3) with a Kolmogorov spectrum and a
    zero inner scale. So the second moment integrand goes as f^(-2/3), and b_2
    grows as f_max^(1/3). The moment has no upper limit of its own. The caller
    must set the top of the grid, from the detector bandwidth or from an inner
    scale. This is why the book sets nu0 to a fixed 550 Hz for its figures
    (printed p. 457 and printed p. 514) instead of computing it. The self-check
    of this module measures the band dependence.

    Parameters:
        freq : numpy.ndarray
            Frequency grid [Hz], ascending. Use a DENSE log-spaced grid. The
            spectrum carries small oscillations above the Fresnel frequency,
            from the exponential part of the confluent hypergeometric function,
            and a coarse grid aliases them into the second moment.
        spectrum : numpy.ndarray
            S_I(f) on that grid, from `irradiance_temporal_spectrum`.

    Returns:
        float
            nu0 [Hz].
    '''
    f = np.asarray(freq, dtype=float)
    s = np.asarray(spectrum, dtype=float)
    b0 = np.trapz(s, f)
    b2 = np.trapz(f ** 2 * s, f)
    return float(np.sqrt(b2 / b0))


def greenwood_frequency(hs, cn2_profile, wavelength, elevation_deg=90.0,
                        wind_profile=None):
    '''
    Return the Greenwood frequency f_G [Hz] of a slant path.

    formula:
        tau0 = [ 2.91 k^2 airmass INTEGRAL Cn2(h) V^(5/3)(h) dh ]^(-3/5)
        f_G  = 1 / tau0,   k = 2 pi/lambda,   airmass = 1/sin(elevation)
    Source: Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196, Ch. 14,
    Eq. (38), printed p. 622, restated as Ch. 14, Eq. (98), printed p. 637. The
    book writes the integral along the path z. A slant path has dz = airmass dh,
    which gives the single airmass factor above. The book states at printed
    p. 623 that "the Greenwood frequency is simply the reciprocal of the time
    constant", and that the time constant is of the order of milliseconds.

    The Greenwood frequency is the servo bandwidth that an adaptive-optics
    system needs. Below it the correction lags the turbulence.

    Parameters:
        hs : numpy.ndarray
            Heights above the ground station [m], ascending.
        cn2_profile : numpy.ndarray
            Zenith Cn2(h) on the hs grid [m^-2/3].
        wavelength : float
            Optical wavelength [m].
        elevation_deg : float
            Elevation angle above the horizon [deg]. 90 is the zenith.
        wind_profile : numpy.ndarray, optional
            Transverse wind speed V(h) on the hs grid [m/s]. None takes the
            Bufton wind of `olb._deps.v_wind`, which Ch. 12, Eq. (3), printed
            p. 481, prints with the same three constants.

    Returns:
        float
            The Greenwood frequency [Hz].
    '''
    hs = np.asarray(hs, dtype=float)
    cn2 = np.asarray(cn2_profile, dtype=float)
    wind = v_wind(hs) if wind_profile is None else np.asarray(wind_profile,
                                                              dtype=float)
    k = wavenumber(wavelength)
    airmass = 1.0 / np.sin(np.radians(elevation_deg))
    integral = np.trapz(cn2 * wind ** (5.0 / 3.0), hs) * airmass
    tau0 = (GREENWOOD_CONSTANT * k ** 2 * integral) ** (-3.0 / 5.0)
    return float(1.0 / tau0)


def coherence_time(greenwood_hz):
    '''
    Return the Greenwood time constant tau0 [s] from the Greenwood frequency.

    formula:
        tau0 = 1 / f_G
    Source: Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196, Ch. 14,
    text below Eq. (39), printed p. 623: "the Greenwood frequency is simply the
    reciprocal of the time constant". Ch. 14, Eq. (38), printed p. 622, defines
    tau0 itself. The book states on printed p. 623 that tau0 is typically of the
    order of milliseconds.

    Parameters:
        greenwood_hz : float or numpy.ndarray
            The Greenwood frequency [Hz].

    Returns:
        float or numpy.ndarray
            tau0 [s].
    '''
    return 1.0 / np.asarray(greenwood_hz, dtype=float)


if __name__ == '__main__':
    from scipy.special import gamma as _g

    from .aperture import averaged_index
    from .beam import beam_params
    from .scintillation import scintillation_index
    from .wander import plane_fried_parameter_slant

    def _grid(top, n=200001, bottom=1e-4):
        '''A dense log-spaced frequency grid from bottom to top [Hz].'''
        return np.logspace(np.log10(bottom), np.log10(top), n)

    def _integrate(grid, values):
        '''
        Integrate a spectrum over (0, top] from its values on a log grid.

        The spectrum is flat below the Fresnel frequency (Ch. 8, text below
        Eq. (57), printed 283), so the head below the first grid point is the
        rectangle values[0] * grid[0].
        '''
        return float(np.trapz(values, grid) + values[0] * grid[0])

    # === physics ============================================================

    # Eq. (27), printed 73: the frozen-flow map. A 10 m/s wind carries a 1 m
    # eddy past the receiver in 0.1 s, so kappa = 2 pi f/V must give 2 pi.
    kap = taylor_wavenumber(10.0, 10.0)
    print(f"[physics] Eq.(27) Taylor kappa  err = "
          f"{abs(kap - 2.0 * np.pi) / (2.0 * np.pi):.3e}")
    assert np.isclose(kap, 2.0 * np.pi)

    # The case for every weak spectral check below.
    LAM, LEN, CN2, WIND = 1.55e-6, 2000.0, 1.0e-16, 10.0
    S2R = float(rytov_variance(LAM, LEN, CN2, wave='plane'))
    OMT = float(_angular_fresnel(WIND, LAM, LEN))
    print(f"[physics] case: lambda=1.55 um, L=2000 m, Cn2=1e-16, V=10 m/s "
          f"-> sigma_R^2={S2R:.4f}, f_t={OMT / (2 * np.pi):.1f} Hz")

    # Eq. (55), printed 282: the spectrum integrates to the scintillation index.
    # The weak plane-wave index is sigma_R^2 (Ch. 8, text at Eq. (57)).
    G9 = _grid(1e9)

    def s_plane(f):
        return irradiance_temporal_spectrum(f, WIND, LAM, LEN, CN2,
                                            wave='plane')

    got = _integrate(G9, s_plane(G9))
    print(f"[physics] Eq.(55) plane INT S df = {got / S2R:.6f} sigma_R^2   "
          f"err = {abs(got - S2R) / S2R:.4f}")
    assert abs(got - S2R) / S2R < 0.01

    # The weak spherical index is 0.4 sigma_R^2 (Ch. 8, text below Eq. (59),
    # printed 284).
    def s_sph(f):
        return irradiance_temporal_spectrum(f, WIND, LAM, LEN, CN2,
                                            wave='spherical')

    got = _integrate(G9, s_sph(G9))
    print(f"[physics] Eq.(55) sph   INT S df = {got / S2R:.6f} sigma_R^2   "
          f"err = {abs(got - 0.4 * S2R) / (0.4 * S2R):.4f}")
    assert abs(got - 0.4 * S2R) / (0.4 * S2R) < 0.01

    # Eq. (63), printed 285: the same identity for a Gaussian beam. The target
    # is the covariance of Eq. (63) at tau = 0.
    BEAM = beam_params(0.02, LAM, LEN)
    dt_b, q1_b, q2_b = _weak_beam_group('gaussian', BEAM)
    target_b = 3.87 * S2R * np.real(
        (1j ** (5.0 / 6.0))
        * (1.0 - (BEAM.theta_bar + 1j * BEAM.lam) * dt_b) ** (5.0 / 6.0)
        - (BEAM.lam * dt_b) ** (5.0 / 6.0))

    def s_beam(f):
        return irradiance_temporal_spectrum(f, WIND, LAM, LEN, CN2,
                                            wave='gaussian', beam=BEAM)

    got = _integrate(G9, s_beam(G9))
    print(f"[physics] Eq.(63) beam  INT S df = {got:.6e}   "
          f"B(0)={float(target_b):.6e}   err = "
          f"{abs(got - target_b) / target_b:.4f}")
    assert abs(got - target_b) / target_b < 0.01

    # Fig. 8.13, printed 283: the scaled plane spectrum reads 2.50 at low
    # frequency, and the spherical one reads about 1.67.
    lowf = np.array([1e-4 * OMT / (2 * np.pi)])
    flat_pl = float(s_plane(lowf)[0]) * OMT / S2R
    flat_sp = float(s_sph(lowf)[0]) * OMT / (0.4 * S2R)
    print(f"[physics] Fig.8.13 flat plane   = {flat_pl:.3f} (book 2.50)")
    print(f"[physics] Fig.8.13 flat sph     = {flat_sp:.3f} (book 1.67)")
    assert abs(flat_pl - 2.50) < 0.02 and abs(flat_sp - 1.67) < 0.04

    # The strong branch, Ch. 10, Eqs. (93)-(97), printed 421-422. Its integral
    # must give the aperture-averaged scintillation index of Ch. 10, Eq. (69),
    # printed 413.
    CN2_S, LEN_S = 5.0e-14, 2000.0
    S2R_S = float(rytov_variance(LAM, LEN_S, CN2_S, wave='plane'))
    G_S = _grid(1e7, n=331)
    for diam in (None, 0.10):
        vals = irradiance_temporal_spectrum(G_S, WIND, LAM, LEN_S, CN2_S,
                                            regime='strong', D=diam)
        want = float(averaged_index(0.0 if diam is None else diam, LAM, LEN_S,
                                    CN2_S, wave='plane', regime='strong'))
        got = _integrate(G_S, vals)
        tag = 'D=0   ' if diam is None else 'D=10cm'
        print(f"[physics] Eq.(97) strong {tag} INT S df = {got:.5f}   "
              f"Eq.(69)={want:.5f}   err = {abs(got - want) / want:.4f}")
        assert abs(got - want) / want < 0.02
    print(f"[physics] strong case: sigma_R^2 = {S2R_S:.2f} (saturation)")

    # Eq. (38), printed 622, against Eq. (39), printed 622: a constant wind ties
    # the Greenwood time constant to the Fried parameter, tau0 = 0.32 r0/V.
    HS = np.linspace(0.0, 20000.0, 4001)
    CN2P = (0.00594 * (21.0 / 27.0) ** 2 * (1e-5 * HS) ** 10
            * np.exp(-HS / 1000.0)
            + 2.7e-16 * np.exp(-HS / 1500.0)
            + 1.7e-14 * np.exp(-HS / 100.0))
    VCONST = 20.0
    fg_const = greenwood_frequency(HS, CN2P, LAM, 90.0,
                                   np.full_like(HS, VCONST))
    r0 = plane_fried_parameter_slant(LAM, HS, CN2P, 90.0)
    ratio = coherence_time(fg_const) * VCONST / r0
    print(f"[physics] Eq.(39) tau0 V/r0     = {ratio:.4f} (book 0.32)   "
          f"err = {abs(ratio - 0.32) / 0.32:.4f}")
    assert abs(ratio - 0.32) / 0.32 < 0.03

    # The Bufton wind of Ch. 12, Eq. (3), printed 481, on the same HV profile.
    # The book states at printed 623 that tau0 is of the order of milliseconds.
    fg = greenwood_frequency(HS, CN2P, LAM, 90.0)
    print(f"[physics] HV Greenwood f_G      = {fg:.1f} Hz, "
          f"tau0 = {coherence_time(fg) * 1e3:.3f} ms (book: order ms)")
    assert 0.1e-3 < coherence_time(fg) < 20e-3

    # === reduction ==========================================================

    # Eq. (65), printed 285, reduces to Eq. (57), printed 283, for a plane wave.
    # The book prints 6.95 and 0.72; this code derives 3.90/d_t^(5/6) and
    # C q^(4/3).
    dt_p, q1_p, _ = _weak_beam_group('plane', None)
    amp_p = _PSD_AMPLITUDE / dt_p ** (5.0 / 6.0)
    c_p = _PSD_TAIL_CONSTANT * q1_p ** (4.0 / 3.0)
    print(f"[reduce ] Eq.(57) amplitude     = {amp_p:.4f} (book 6.95), "
          f"q1 = {q1_p.real:.4f} (book 0.5000), c = {c_p.real:.4f} (book 0.72)")
    assert abs(amp_p - 6.95) < 0.05 and abs(c_p.real - 0.72) < 0.005

    dt_s, q1_s, _ = _weak_beam_group('spherical', None)
    amp_s = _PSD_AMPLITUDE / dt_s ** (5.0 / 6.0)
    c_s = _PSD_TAIL_CONSTANT * q1_s ** (4.0 / 3.0)
    print(f"[reduce ] Eq.(59) amplitude     = {amp_s:.4f} (book 5.47), "
          f"q1 = {q1_s.real:.4f} (book 0.2222), c = {c_s.real:.4f} (book 0.24)")
    assert abs(amp_s - 5.47) < 0.05 and abs(c_s.real - 0.24) < 0.005

    # The literal Eq. (57) of the book, against this module, over the band.
    def book_57(f):
        u = 2.0 * np.pi * np.asarray(f, dtype=float) / OMT
        zz = -1j * u ** 2 / 2.0
        val = (u ** (-8.0 / 3.0)
               * (1.0 - _hyp1f1(-5.0 / 6.0, -1.0 / 3.0, zz))
               - 0.72 * (1j ** (4.0 / 3.0)) * _hyp1f1(0.5, 7.0 / 3.0, zz))
        return 6.95 * S2R / OMT * np.real(val)

    core = np.logspace(-1, np.log10(3.0 * OMT / (2 * np.pi)), 40)
    rel = np.max(np.abs(irradiance_temporal_spectrum(core, WIND, LAM, LEN, CN2,
                                                     wave='plane')
                        - book_57(core)) / np.abs(book_57(core)))
    print(f"[reduce ] Eq.(57) literal book  max err = {rel:.4f} up to 3 f_t")
    assert rel < 0.005

    # Above about 100 Fresnel frequencies the printed 0.72 leaves a residual
    # 1/omega tail that takes the printed spectrum NEGATIVE. A power spectral
    # density cannot go below zero, so the printed two-figure constant is the
    # cause. The derived C keeps the spectrum positive everywhere.
    tail = np.array([1e3, 1e4, 1e5])
    print("[reduce ] Eq.(57) far tail      f[Hz]  this module    printed 0.72")
    for f, a, b in zip(tail,
                       irradiance_temporal_spectrum(tail, WIND, LAM, LEN, CN2,
                                                    wave='plane'),
                       book_57(tail)):
        print(f"[reduce ]                   {f:9.0f}  {a: .5e}  {b: .5e}")
    assert np.all(irradiance_temporal_spectrum(tail, WIND, LAM, LEN, CN2,
                                               wave='plane') > 0.0)
    assert book_57(tail)[-1] < 0.0

    fs = np.logspace(0, 4, 40)

    # A Gaussian beam with a large Lambda goes to the spherical wave, and a
    # small Lambda with Theta = 1 goes to the plane wave. Ch. 8, text below
    # Eq. (65), printed 285, states that Eq. (65) holds both limits.
    # The Fresnel ratio Lambda must be small enough that Lambda d_t^2 u^2 stays
    # below one over the whole test band, because that product is the argument
    # of the second hypergeometric group.
    thin = beam_params(2.0e-5, LAM, LEN)     # Lambda0 >> 1 -> spherical
    fat = beam_params(5.0, LAM, LEN)         # Lambda0 << 1 -> plane
    lim_fs = np.logspace(-1, np.log10(10.0 * OMT / (2 * np.pi)), 40)
    for tag, bm, ref in (('->sph  ', thin, s_sph), ('->plane', fat, s_plane)):
        a = irradiance_temporal_spectrum(lim_fs, WIND, LAM, LEN, CN2,
                                         wave='gaussian', beam=bm)
        b = ref(lim_fs)
        e = np.max(np.abs(a - b) / np.abs(b))
        print(f"[reduce ] Eq.(65) beam {tag}  max err = {e:.4f} up to 10 f_t")
        assert e < 0.05

    # Eq. (73), printed 514: the quasi-frequency against direct quadrature of
    # its defining moment ratio. The quadrature is scipy adaptive quadrature
    # over the SAME band, decade by decade, so it is an independent measure of
    # the trapezium sum that `quasi_frequency` makes.
    # The band stops at 1 kHz, about 14 Fresnel frequencies. Above that the
    # spectrum carries small oscillations from the exponential part of the
    # confluent hypergeometric function, and adaptive quadrature aliases them,
    # so the two methods stop agreeing at the 1e-6 level.
    LO, HI = 1e-3, 1e3
    grid = np.logspace(np.log10(LO), np.log10(HI), 200001)
    spec = irradiance_temporal_spectrum(grid, WIND, LAM, LEN, CN2, wave='plane')
    nu0_grid = quasi_frequency(grid, spec)

    def _moment(power):
        edges = np.logspace(np.log10(LO), np.log10(HI), 13)
        return sum(quad(lambda f: f ** power
                        * float(s_plane(np.array([f]))[0]),
                        lo, hi, limit=200, epsabs=0.0, epsrel=1e-12)[0]
                   for lo, hi in zip(edges[:-1], edges[1:]))

    nu0_quad = np.sqrt(_moment(2.0) / _moment(0.0))
    print(f"[reduce ] Eq.(73) nu0           = {nu0_grid:.6f} Hz vs quadrature "
          f"{nu0_quad:.6f} Hz   err = "
          f"{abs(nu0_grid - nu0_quad) / nu0_quad:.3e}")
    assert abs(nu0_grid - nu0_quad) / nu0_quad < 1e-6

    # The second moment of the plane spectrum has no upper limit of its own,
    # because the spectrum decays as f^(-8/3) and the moment integrand as
    # f^(-2/3). So b_2 grows as f_max^(1/3) and nu0 depends on the band. Print
    # that, so no caller reads nu0 as band-independent.
    prev = None
    for top in (1e4, 1e5, 1e6, 1e7):
        g = np.logspace(-3, np.log10(top), 400001)
        s = irradiance_temporal_spectrum(g, WIND, LAM, LEN, CN2, wave='plane')
        nu = quasi_frequency(g, s)
        step = '' if prev is None else f"  (x{nu / prev:.3f} per decade)"
        print(f"[reduce ] nu0 band 0-{top:.0e} Hz  = {nu:8.2f} Hz{step}")
        prev = nu

    # The strong-regime covariance at tau = 0 against Ch. 10, Eq. (69).
    b_zero = float(_strong_covariance(0.0, S2R_S, 0.0))
    want = float(averaged_index(0.0, LAM, LEN_S, CN2_S, wave='plane',
                                regime='strong'))
    print(f"[reduce ] Eq.(93) B(0)          = {b_zero:.5f} vs Eq.(69) "
          f"{want:.5f}   err = {abs(b_zero - want) / want:.4f}")
    assert abs(b_zero - want) / want < 0.02
    print(f"[reduce ] small-scale s->0 limit= "
          f"{0.5 * _g(5.0 / 6.0) * 2.0 ** (5.0 / 6.0):.5f} (exact 1 wanted)")

    # The strong branch against the weak branch at a small sigma_R^2. The two
    # carry the SAME power, but NOT the same shape. Ch. 8.5 comes from the
    # Rytov covariance; Ch. 9.8 comes from the two-scale extended Rytov
    # covariance. So only the integrals must agree. The loop measures both.
    weak_cn2 = 1.0e-17
    s2r_w = float(rytov_variance(LAM, LEN, weak_cn2))
    ga = _grid(1e8, n=40001)
    ia = _integrate(ga, irradiance_temporal_spectrum(ga, WIND, LAM, LEN,
                                                     weak_cn2, wave='plane'))
    gb = _grid(1e7, n=331)
    ib = _integrate(gb, irradiance_temporal_spectrum(gb, WIND, LAM, LEN,
                                                     weak_cn2,
                                                     regime='strong'))
    print(f"[reduce ] strong vs weak (s_R^2={s2r_w:.4f}) power "
          f"{ib:.6e} vs {ia:.6e}   err = {abs(ib - ia) / ia:.4f}")
    assert abs(ib - ia) / ia < 0.02
    shape = np.array([1.0, 100.0, 1000.0])
    ratio = (irradiance_temporal_spectrum(shape, WIND, LAM, LEN, weak_cn2,
                                          regime='strong')
             / irradiance_temporal_spectrum(shape, WIND, LAM, LEN, weak_cn2,
                                            wave='plane'))
    print(f"[reduce ] strong/weak shape at 1, 100, 1000 Hz = "
          f"{ratio[0]:.3f}, {ratio[1]:.3f}, {ratio[2]:.3f}")

    # scintillation_index agrees with the spectral integral in the weak regime.
    idx = float(scintillation_index(LAM, LEN, CN2, wave='plane',
                                    regime='weak'))
    print(f"[reduce ] index vs INT S df     = {idx:.6e} vs "
          f"{_integrate(G9, s_plane(G9)):.6e}")

    # Every refusal the book forces.
    for kwargs, why in ((dict(l0=0.005), 'inner scale'),
                        (dict(L0=10.0), 'outer scale'),
                        (dict(regime='strong', wave='spherical'),
                         'strong spherical'),
                        (dict(D=0.1), 'weak aperture')):
        try:
            irradiance_temporal_spectrum(100.0, WIND, LAM, LEN, CN2, **kwargs)
        except NotImplementedError:
            print(f"[reduce ] refuses {why:<18} ok")
        else:
            raise AssertionError(f'{why} must raise NotImplementedError')

    print("self-check passed")
