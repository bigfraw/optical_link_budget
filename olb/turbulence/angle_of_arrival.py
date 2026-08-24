'''
Received tip-tilt (angle of arrival) of a Gaussian beam over a turbulent path.

This module gives the received tip-tilt angle variance of a Gaussian beam. The
tip-tilt is the random slope of the arriving wavefront. A receive telescope
focuses the beam onto a fibre tip. A tip-tilt of angle theta moves the focal
spot by f*theta, with f the focal length. So the received tip-tilt drives the
fibre-coupling loss. The Term factories live in olb.models.coupling.

The module is pure physics. It imports only numpy and olb._deps. It does not
import the scenario, the terminal, or the results.

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
     This is DEFERRED. See aperture_arrival_angle_variance below.
'''

import numpy as np

from .._deps import beam_wander_variance


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
    Return the aperture angle-of-arrival tip-tilt variance. DEFERRED.

    This is the second, smaller aperture angle-of-arrival "corrugation" tilt. It
    is the classic plane-wave form. The explicit Andrews treatment is DEFERRED.
    The repo owner will specify it. Do not guess the coefficient.

    The working received tip-tilt used by the coupling Terms is the beam-wander
    term (wander_arrival_angle_variance) only. See docs/andrews-crosscheck.md
    batch 2.
    '''
    raise NotImplementedError(
        "Aperture angle-of-arrival tip-tilt is DEFERRED. The repo owner will "
        "specify the explicit Andrews and Phillips treatment. Do not guess the "
        "coefficient. The working received tip-tilt is the beam-wander term "
        "(wander_arrival_angle_variance) only."
    )


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

    # The aperture angle-of-arrival term is DEFERRED: it must raise.
    try:
        aperture_arrival_angle_variance(0.2, 0.1, lam)
        raise AssertionError("aperture AoA must be deferred (NotImplementedError)")
    except NotImplementedError:
        pass

    print(f"wander tilt variance (3 km, Cn2=1e-14) = {v:.3e} rad^2")
    print(f"  radial 1-sigma = {np.sqrt(v) * 1e6:.3f} urad  "
          f"per-axis 1-sigma = {np.sqrt(v / 2) * 1e6:.3f} urad")
    print("self-check passed")
