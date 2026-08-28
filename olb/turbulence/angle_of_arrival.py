'''
Received tip-tilt (angle of arrival) of a Gaussian beam over a turbulent path.

This module gives the received tip-tilt angle variance of a Gaussian beam. The
tip-tilt is the random slope of the arriving wavefront. A receive telescope
focuses the beam onto a fibre tip. A tip-tilt of angle theta moves the focal
spot by f*theta, with f the focal length. So the received tip-tilt drives the
fibre-coupling loss. The Term factories live in olb.models.coupling.

The module is pure physics. It imports only numpy and olb.turbulence.coupled_flux.
It does not import the scenario, the terminal, or the results.

Two contributions:

  A. Beam-wander arrival tilt (the DOMINANT term, and the one this module gives).
     The turbulence moves the beam centroid at the receiver by a random offset
     r_c. That offset is an apparent tilt r_c/L of the arriving beam, with L the
     path length. The received radial (2-axis) tilt variance is
         sigma2_theta = <r_c^2> / L^2.
     The kernel beam_wander_variance integrates the free-space beam WIDTH profile
     w(z) along the path, so the result is Gaussian-beam-correct.
     Source: Dios et al., Applied Optics 43 (2004) 3866. DOI 10.1364/AO.43.003866.

  B. Aperture angle-of-arrival "corrugation" tilt (a second, smaller term).
     The wavefront that arrives at the receive aperture is corrugated, so its
     mean slope across the pupil is not zero. The per-axis variance of that
     slope is
         sigma2_theta = 2.91 Cn2 L D^(-1/3) = 0.174 (D/r0)^(5/3)(lambda/D)^2.
     Source: Andrews and Phillips, Laser Beam Propagation through Random Media,
     2nd ed. (2005), DOI 10.1117/3.626196, Ch. 6, Eq. (84), printed p. 201.
     This is the GRADIENT tilt of a centroid tracker, not the Noll Zernike
     tilt. See aperture_arrival_angle_variance below.
'''

import numpy as np

from .coupled_flux import beam_wander_variance
from .andrews.structure import (
    angle_of_arrival_variance as _andrews_angle_of_arrival_variance,
)


def wander_arrival_angle_variance(L, cn2_slant, w_profile, hs):
    '''
    Return the received beam-wander tip-tilt variance (radial, 2-axis) [rad^2].

    Reuse the beam-wander kernel. It integrates the free-space beam width profile
    w(z) along the path, so it is Gaussian-beam-correct. The beam-wander offset
    variance <r_c^2> maps to an apparent arrival tilt through the path length L:
        sigma2_theta = <r_c^2> / L^2.
    The result is the radial (2-axis) variance. The per-axis variance is one half
    of it. Source: Dios et al., Applied Optics 43 (2004) 3866, DOI
    10.1364/AO.43.003866 (the beam-wander offset variance).

    Parameters:
        L : float
            Path length from the transmitter to the receiver [m].
        cn2_slant : numpy.ndarray
            Cn2 along the path on the hs grid [m^-2/3].
        w_profile : numpy.ndarray
            Free-space beam radius w(z) along the path on the hs grid [m].
        hs : numpy.ndarray
            Distance along the path from the transmitter [m].

    Returns:
        float
            Radial (2-axis) received tip-tilt variance [rad^2].
    '''
    r_c2 = beam_wander_variance(L, cn2_slant, w_profile, hs)
    return float(np.squeeze(r_c2)) / float(L) ** 2


def aperture_arrival_angle_variance(D, r0, wavelength):
    '''
    Return the aperture angle-of-arrival tip-tilt variance (per axis) [rad^2].

    This is the second, smaller aperture angle-of-arrival "corrugation" tilt. It
    is separate from the beam-wander arrival tilt above.

    TILT DEFINITION - THE OWNER MADE THIS CHOICE. This function returns the
    ANDREWS GRADIENT TILT (G-tilt), which is what a centroid tracker measures.
    Andrews defines the tilt as the total phase difference across the pupil
    divided by the pupil width. It is NOT the Noll Zernike tilt.

    formula:
        <beta_a^2> = 2.91 Cn2 L (2 W_G)^(-1/3)
                   = 0.174 (D/r0)^(5/3) (lambda/D)^2      per axis
    Source: Andrews and Phillips, Laser Beam Propagation through Random Media,
    2nd ed. (2005), DOI 10.1117/3.626196, Ch. 6, Eq. (84), printed p. 201, with
    the definition Ch. 6, Eq. (82), printed p. 200. The slant-path version is
    Ch. 12, Eq. (28), printed p. 492.

    THE RECAST. Put Cn2 L = r0^(-5/3)/(0.423 k^2) and k = 2*pi/lambda into
    Eq. (84). Then 2.91/(0.423 * 4 pi^2) = 0.1743.

    THE ALTERNATIVE. The Noll Zernike tilt gives 0.182 (D/r0)^(5/3)(lambda/D)^2
    (Noll, JOSA 66 (1976) 207, DOI 10.1364/JOSA.66.000207). A full-book search
    finds no 0.182 in Andrews and Phillips. Note that `olb/turbulence/ao.py`
    uses the NOLL convention (1.0299 and 0.134), so a caller that mixes the two
    must say which tilt it means. See Conflict C-04 in
    docs/andrews-crosscheck.md.

    Parameters:
        D : float or numpy.ndarray
            Receive aperture diameter [m]. The book writes it as 2 W_G.
        r0 : float or numpy.ndarray
            Fried atmospheric coherence width over the same path [m].
        wavelength : float
            Optical wavelength [m].

    Returns:
        float or numpy.ndarray
            PER-AXIS tilt variance [rad^2]. The radial (2-axis) variance is
            twice this value.

    VALIDITY. The result holds only when the Fresnel zone is small against the
    aperture, sqrt(L/k) << D (Ch. 6, text below Eq. (83), printed p. 200). This
    function does not gate on that condition.

    NEW HOME: `olb.turbulence.andrews.structure.angle_of_arrival_variance`. That
    function takes Cn2 and the path length directly, and it also gives the
    inner-scale and outer-scale branches of Ch. 6, Eq. (83).
    '''
    # Rebuild the path moment Cn2 * L from the Fried parameter, so that this
    # signature stays unchanged. r0 = (0.423 k^2 Cn2 L)^(-3/5), so
    # Cn2 L = r0^(-5/3) / (0.423 k^2). The Andrews function uses only the
    # product Cn2 * z, so z = 1 m carries it.
    k = 2.0 * np.pi / wavelength
    cn2_l = np.asarray(r0, dtype=float) ** (-5.0 / 3.0) / (0.423 * k ** 2)
    return _andrews_angle_of_arrival_variance(D, wavelength, 1.0, cn2_l)


if __name__ == '__main__':
    # Pure-physics self-check. Use plain numeric inputs. No scenario import.
    from ..beam import free_space_radius

    lam = 1550e-9
    w0 = 0.02
    L = 3e3
    hs = np.linspace(0.0, L, 200)

    def _variance(cn2, length=L, waist=w0):
        grid = np.linspace(0.0, length, 200)
        cn2_slant = np.full_like(grid, cn2)
        w_profile = free_space_radius(waist, grid, None, lam)
        return wander_arrival_angle_variance(length, cn2_slant, w_profile, grid)

    # A real, positive tilt variance for a turbulent path.
    v = _variance(1e-14)
    assert np.isfinite(v) and v > 0.0, v

    # Stronger turbulence gives a larger tilt variance.
    assert _variance(1e-13) > _variance(1e-15)

    # The tilt variance scales linearly with Cn2 (beam wander is linear in Cn2).
    ratio = _variance(2e-14) / _variance(1e-14)
    assert np.isclose(ratio, 2.0, rtol=1e-6), ratio

    # The aperture angle-of-arrival term. It falls as D^(-1/3) and it grows as
    # the Fried parameter falls.
    D_ap, r0_ap = 0.2, 0.1
    v_ap = aperture_arrival_angle_variance(D_ap, r0_ap, lam)
    assert v_ap > 0.0, v_ap
    assert aperture_arrival_angle_variance(0.4, r0_ap, lam) < v_ap
    assert aperture_arrival_angle_variance(D_ap, 0.05, lam) > v_ap
    # The gradient-tilt recast, Andrews Ch. 6, Eq. (84), printed p. 201.
    recast = 0.174 * (D_ap / r0_ap) ** (5.0 / 3.0) * (lam / D_ap) ** 2
    pct = abs(v_ap - recast) / recast * 100.0
    assert pct < 2.0, pct

    print(f"wander tilt variance (3 km, Cn2=1e-14) = {v:.3e} rad^2")
    print(f"aperture AoA tilt (D=0.2 m, r0=0.1 m) = {v_ap:.3e} rad^2  "
          f"per-axis 1-sigma = {np.sqrt(v_ap) * 1e6:.3f} urad  "
          f"(gradient tilt, {pct:.3f} % from the 0.174 recast)")
    print(f"  radial 1-sigma = {np.sqrt(v) * 1e6:.3f} urad  "
          f"per-axis 1-sigma = {np.sqrt(v / 2) * 1e6:.3f} urad")
    print("self-check passed")
