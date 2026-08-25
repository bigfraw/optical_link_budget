'''
Gaussian-beam parameters of Andrews and Phillips, for any input curvature.

This module gives the four nondimensional beam parameters that the book uses
through every chapter: the input-plane pair (Theta0, Lambda0) and the
output-plane pair (Theta, Lambda). It also gives the free-space beam radius W at
the output plane, and the strong-fluctuation "effective" pair (Theta_e,
Lambda_e).

The module is general in the input curvature. The phase-front radius f0 sets the
beam type:
    f0 = +infinity   collimated beam, Theta0 = 1
    f0 < 0           divergent beam,  Theta0 > 1
    0 < f0           convergent (focused) beam, Theta0 < 1
Source of the classification: Andrews and Phillips, Laser Beam Propagation
through Random Media, 2nd ed. (SPIE Press, 2005), DOI 10.1117/3.626196, Ch. 4,
Sec. 4.4.1, text below Eq. (38), printed p. 93.

This module holds physics only. It imports numpy only. It returns no decibels.

Plane of reference: all quantities are referred to the TRANSMITTER at z = 0 and
the OUTPUT plane at z. No path integral occurs in this module, so no path weight
and no reference-plane choice is made here.
'''

from typing import NamedTuple

import numpy as np


class BeamParams(NamedTuple):
    '''
    The Gaussian-beam parameters at one output plane.

    The field order is a frozen contract across the package. Do not change it.

    Fields:
        theta0 : input-plane curvature (refraction) parameter Theta0.
        lambda0 : input-plane Fresnel ratio Lambda0.
        theta : output-plane refraction parameter Theta.
        theta_bar : the complement 1 - Theta.
        lam : output-plane diffraction parameter Lambda.
        w : free-space beam radius at the output plane [m].

    Source: Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196, Ch. 4,
    Eqs. (33), (44), (45) and (47), printed pp. 92 and 95.
    '''

    theta0: np.ndarray
    lambda0: np.ndarray
    theta: np.ndarray
    theta_bar: np.ndarray
    lam: np.ndarray
    w: np.ndarray


def wavenumber(wavelength):
    '''
    Return the optical wavenumber k = 2*pi/lambda [rad/m].

    Source: Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196, Ch. 4,
    Sec. 4.5, text at Eq. (53), printed p. 96.
    '''
    return 2.0 * np.pi / np.asarray(wavelength, dtype=float)


def beam_params(w0, wavelength, z, f0=np.inf):
    '''
    Return the Gaussian-beam parameters at range z, for any input curvature.

    Parameters:
        w0 : float or numpy.ndarray
            Beam RADIUS (1/e field) at the transmitter [m].
        wavelength : float or numpy.ndarray
            Optical wavelength [m].
        z : float or numpy.ndarray
            Range from the transmitter to the output plane [m].
        f0 : float or numpy.ndarray
            Phase-front radius of curvature at the transmitter [m]. Use
            numpy.inf for a collimated beam, a negative value for a divergent
            beam, and a positive value for a convergent beam.

    Returns:
        BeamParams
            The six fields, each broadcast over the input shapes.

    formula:
        Theta0 = 1 - z/f0,   Lambda0 = 2 z / (k W0^2),   k = 2*pi/lambda
        Theta  = Theta0 / (Theta0^2 + Lambda0^2)
        Lambda = Lambda0 / (Theta0^2 + Lambda0^2)
        Theta_bar = 1 - Theta
        W = W0 sqrt(Theta0^2 + Lambda0^2)
    Source: Andrews and Phillips, Laser Beam Propagation through Random Media,
    2nd ed. (2005), DOI 10.1117/3.626196:
        Theta0, Lambda0  Ch. 4, Eq. (33), printed p. 92
        W                Ch. 4, Eq. (37), printed p. 93
        Theta, Lambda    Ch. 4, Eq. (44), printed p. 95
        Theta_bar        Ch. 4, Eq. (45), printed p. 95
        the output-plane reading Theta = 1 - z/F, Lambda = 2 z/(k W^2)
                         Ch. 4, Eq. (47), printed p. 95
    Ch. 8, Eqs. (5) and (6), printed p. 261, and Ch. 6, Eq. (6), printed p. 183,
    restate the same pair with z = L.

    The plane-wave limit is Theta = 1, Lambda = 0. The spherical-wave limit is
    Theta = Lambda = 0. Source: Ch. 8, text below Eq. (6), printed p. 261.
    '''
    w0 = np.asarray(w0, dtype=float)
    z = np.asarray(z, dtype=float)
    k = wavenumber(wavelength)

    theta0 = 1.0 - z / np.asarray(f0, dtype=float)
    lambda0 = 2.0 * z / (k * w0 ** 2)

    denom = theta0 ** 2 + lambda0 ** 2
    theta = theta0 / denom
    lam = lambda0 / denom
    w = w0 * np.sqrt(denom)
    return BeamParams(theta0, lambda0, theta, 1.0 - theta, lam, w)


def effective_beam_params(bp, sigma2_R):
    '''
    Return the strong-fluctuation effective beam parameters.

    Turbulence spreads the beam and flattens its mean phase front. The effective
    parameters carry that change. They let a weak-fluctuation formula reach into
    the strong-fluctuation regime.

    Parameters:
        bp : BeamParams
            The free-space parameters at the output plane, from beam_params.
        sigma2_R : float or numpy.ndarray
            The plane-wave Rytov variance sigma_R^2 over the same path.

    Returns:
        BeamParams
            The effective parameters. The fields theta and lam hold Theta_e and
            Lambda_e. The field theta_bar holds 1 - Theta_e. The field w holds
            the long-term beam radius W_LT. The input-plane fields theta0 and
            lambda0 PASS THROUGH unchanged, because turbulence does not change
            the transmitter. The book gives no effective input-plane pair.

    formula:
        q = 1.22 sigma_R^(12/5),   so 2q/3 = 0.81 and 4q/3 = 1.63
        W_LT     = W sqrt(1 + 1.63 sigma_R^(12/5) Lambda)
        Theta_e  = (Theta - 0.81 sigma_R^(12/5) Lambda)
                   / (1 + 1.63 sigma_R^(12/5) Lambda)
        Lambda_e = Lambda / (1 + 1.63 sigma_R^(12/5) Lambda)
    Source: Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196, Ch. 7,
    Eq. (57) (W_LT) and Eq. (58) (Theta_e, Lambda_e), printed p. 242. Ch. 9,
    Eqs. (85) and (86), printed p. 349, restate the same pair.

    Note the book writes Theta_e = (Theta_bar + 2 q Lambda/3)/(1 + 4 q Lambda/3)
    with Theta_bar = 1 - Theta. The form above is the same result written with
    Theta, because Theta_e as printed gives 1 - Theta_e. This module keeps the
    olb convention: the returned theta field is Theta_e itself.
    '''
    sigma2_R = np.asarray(sigma2_R, dtype=float)
    # sigma_R^(12/5) written from the VARIANCE: (sigma_R^2)^(6/5).
    s125 = sigma2_R ** (6.0 / 5.0)
    denom = 1.0 + 1.63 * s125 * bp.lam
    theta_e = (bp.theta - 0.81 * s125 * bp.lam) / denom
    lam_e = bp.lam / denom
    w_lt = bp.w * np.sqrt(denom)
    return BeamParams(bp.theta0, bp.lambda0, theta_e, 1.0 - theta_e, lam_e,
                      w_lt)


if __name__ == '__main__':
    # ---------------- physics self-checks ----------------
    lam_m = 1550e-9
    k = wavenumber(lam_m)

    # Plane-wave limit: a very wide beam gives Theta -> 1 and Lambda -> 0.
    bp_plane = beam_params(50.0, lam_m, 2000.0)
    assert abs(bp_plane.theta - 1.0) < 1e-6, bp_plane.theta
    assert bp_plane.lam < 1e-6, bp_plane.lam

    # Spherical-wave limit: a point source gives Theta -> 0 and Lambda -> 0.
    bp_sph = beam_params(1e-5, lam_m, 2000.0)
    assert bp_sph.theta < 1e-6, bp_sph.theta
    assert bp_sph.lam < 1e-4, bp_sph.lam

    # Beam type from the input curvature. Ch. 4, printed p. 93.
    assert beam_params(0.05, lam_m, 2000.0, np.inf).theta0 == 1.0
    assert beam_params(0.05, lam_m, 2000.0, -1000.0).theta0 > 1.0
    assert beam_params(0.05, lam_m, 2000.0, 4000.0).theta0 < 1.0

    # Theta_bar is the complement of Theta. Ch. 4, Eq. (45).
    bp = beam_params(0.05, lam_m, 2000.0)
    assert abs(bp.theta_bar - (1.0 - bp.theta)) < 1e-15

    # The output-plane reading of Eq. (47): Lambda = 2 z/(k W^2).
    lam_from_w = 2.0 * 2000.0 / (k * bp.w ** 2)
    assert abs(lam_from_w - bp.lam) / bp.lam < 1e-12, (lam_from_w, bp.lam)

    # A convergent beam focused at the output plane has the smallest radius.
    w_focus = beam_params(0.05, lam_m, 2000.0, 2000.0).w
    assert w_focus < bp.w, (w_focus, bp.w)

    # Effective parameters: weak turbulence leaves the beam alone.
    bp_e_weak = effective_beam_params(bp, 1e-8)
    assert abs(bp_e_weak.theta - bp.theta) < 1e-8
    assert abs(bp_e_weak.lam - bp.lam) < 1e-8
    # Strong turbulence spreads the beam, so W_LT > W and Lambda_e < Lambda.
    bp_e = effective_beam_params(bp, 4.0)
    assert bp_e.w > bp.w, (bp_e.w, bp.w)
    assert bp_e.lam < bp.lam, (bp_e.lam, bp.lam)
    # Ch. 7, Eq. (86): Lambda_e = 2 L/(k W_LT^2).
    lam_e_from_w = 2.0 * 2000.0 / (k * bp_e.w ** 2)
    assert abs(lam_e_from_w - bp_e.lam) / bp_e.lam < 1e-12

    # Arrays broadcast.
    zs = np.array([1000.0, 2000.0, 4000.0])
    bp_arr = beam_params(0.05, lam_m, zs)
    assert bp_arr.lam.shape == zs.shape
    assert abs(bp_arr.theta[1] - bp.theta) < 1e-15

    # ---------------- REDUCTION checks ----------------
    from .. import gaussian_fried as gf

    # 1. beam_params(f0 = inf) reproduces gaussian_fried.output_beam_params.
    th_old, lam_old = gf.output_beam_params(2000.0, 0.05, lam_m)
    err_theta = abs(bp.theta - th_old)
    err_lam = abs(bp.lam - lam_old)
    assert err_theta < 1e-12 and err_lam < 1e-12, (err_theta, err_lam)
    print(f'REDUCTION output_beam_params : dTheta = {err_theta:.3e}  '
          f'dLambda = {err_lam:.3e}  (target 1e-12)')

    # 2. effective_beam_params reproduces the old gaussian_fried values.
    cn2_ref = 1e-14
    sigma2_R_ref = gf.rytov_std(2000.0, cn2_ref, lam_m) ** 2
    th_e_old, lam_e_old = gf.effective_beam_params(2000.0, 0.05, cn2_ref, lam_m)
    bp_e_ref = effective_beam_params(bp, sigma2_R_ref)
    err_th_e = abs(bp_e_ref.theta - th_e_old)
    err_lam_e = abs(bp_e_ref.lam - lam_e_old)
    assert err_th_e < 1e-9 and err_lam_e < 1e-9, (err_th_e, err_lam_e)
    print(f'REDUCTION effective_beam_params : dTheta_e = {err_th_e:.3e}  '
          f'dLambda_e = {err_lam_e:.3e}  (target 1e-9)')

    print(f'collimated 2 km w0=5 cm : Theta = {bp.theta:.5f}  '
          f'Lambda = {bp.lam:.5f}  W = {bp.w * 100:.2f} cm')
    print('self-check passed')
