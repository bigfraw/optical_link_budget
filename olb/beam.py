'''
Deliberate transmit-beam divergence, as a virtual waist behind the aperture.

A collimated transmitter of aperture radius w0 has the diffraction-limited
far-field half-angle divergence lambda/(pi*w0). To send a wider beam on purpose
(for example to relax the pointing budget), recast the diverged transmitter as
an ordinary Gaussian beam. The beam starts from a virtual waist behind the
aperture:

    w_v = lambda / (pi * theta)
    d   = zR(w_v) * sqrt((w0/w_v)^2 - 1)        # virtual distance behind aperture

Its free-space radius at range z is then gaussz(w_v, d + z). No separate
curvature bookkeeping is needed; the ordinary Gaussian machinery does the rest.

theta cannot be below the aperture diffraction limit lambda/(pi*w0). theta equal
to that limit is the collimated case, which returns (w0, 0). None also means
collimated.

This module is pure physics. It imports only numpy and olb._deps. It does not
import the scenario, the results, or the assumptions.
'''

import numpy as np

from ._deps import gaussz, zR, w0_to_div


def virtual_waist(w0, divergence_rad=None, wavelength=1550e-9):
    '''
    Recast a diverged transmitter as a Gaussian beam from a virtual waist.

    Parameters:
        w0 : float
            Beam radius at the transmit aperture [m].
        divergence_rad : float, optional
            Far-field HALF-angle divergence [rad]. None (or the diffraction
            limit) means collimated.
        wavelength : float
            Wavelength [m].

    Returns:
        tuple
            (w_v, d) : the virtual waist radius [m] and its distance behind the
            aperture [m]. Collimated returns (w0, 0.0).
    '''
    theta_min = w0_to_div(w0, wavelength)
    if divergence_rad is None:
        return w0, 0.0
    if divergence_rad < theta_min * (1 - 1e-9):
        raise ValueError(
            f"divergence {divergence_rad * 1e6:.2f} urad is below the diffraction "
            f"limit {theta_min * 1e6:.2f} urad of a {w0} m aperture. Widen the "
            "aperture or accept the larger divergence."
        )
    if np.isclose(divergence_rad, theta_min):
        return w0, 0.0

    w_v = wavelength / (np.pi * divergence_rad)
    d = zR(w_v, wavelength) * np.sqrt((w0 / w_v) ** 2 - 1)
    return w_v, d


def free_space_radius(w0, z, divergence_rad=None, wavelength=1550e-9):
    '''
    Turbulence-free beam radius at range z, for a transmitter of aperture radius
    w0 and far-field half-angle divergence divergence_rad.

    Reduces exactly to gaussz(w0, z, wavelength) for a collimated beam.

    Parameters:
        w0 : float
            Beam radius at the transmit aperture [m].
        z : float or numpy.ndarray
            Range from the aperture [m].
        divergence_rad : float, optional
            Far-field half-angle divergence [rad]. None = collimated.
        wavelength : float
            Wavelength [m].

    Returns:
        float or numpy.ndarray
            Free-space beam radius [m].
    '''
    w_v, d = virtual_waist(w0, divergence_rad, wavelength)
    return gaussz(w_v, d + np.asarray(z, dtype=float), wavelength)


if __name__ == '__main__':
    w0 = 0.05
    lam = 1550e-9
    theta_min = w0_to_div(w0, lam)

    # Collimated request reproduces plain gaussz exactly.
    for z in (0.0, 1e3, 1e6):
        assert np.isclose(free_space_radius(w0, z, None, lam), gaussz(w0, z, lam))
    # theta equal to the diffraction limit is also the collimated case.
    assert np.isclose(free_space_radius(w0, 1e6, theta_min, lam), gaussz(w0, 1e6, lam))

    # The virtual-waist recast reproduces the aperture size and the far field.
    for factor in (1.0, 2.0, 5.0):
        theta = factor * theta_min
        assert np.isclose(free_space_radius(w0, 0.0, theta, lam), w0)
        assert np.isclose(free_space_radius(w0, 1e9, theta, lam) / 1e9, theta, rtol=1e-6)

    # A diverged beam is wider at range than a collimated one.
    assert free_space_radius(w0, 1e6, 5 * theta_min, lam) > free_space_radius(w0, 1e6, None, lam)

    # Sub-diffraction divergence is not a beam.
    try:
        free_space_radius(w0, 1e3, 0.5 * theta_min, lam)
        raise AssertionError("sub-diffraction divergence should raise ValueError")
    except ValueError:
        pass

    print(f"diffraction limit for w0={w0} m: {theta_min * 1e6:.2f} urad")
    print(f"collimated w at 600 km: {free_space_radius(w0, 600e3, None, lam):.2f} m")
    print(f"5x diverged w at 600 km: {free_space_radius(w0, 600e3, 5 * theta_min, lam):.2f} m")
    print("self-check passed")
