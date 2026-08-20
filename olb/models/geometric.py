'''
Geometric (free-space spreading) loss for an FSO link.

Assumes a Gaussian transmit beam and a circular receive aperture, optionally
with a central circular obscuration (e.g. a Cassegrain secondary mirror). The
fraction of transmitted power captured by the receiver is the Gaussian-beam
power enclosed between the obscuration radius and the aperture radius:

    eta = exp(-2*b_obs**2 / w(z)**2) - exp(-2*a_rx**2 / w(z)**2)

where a_rx is the receiver aperture radius, b_obs = obscuration_ratio*a_rx is
the obscuration radius, and w(z) is the 1/e^2 Gaussian beam radius at range z
(borrowed from _deps.gaussz). Loss is deterministic and reported as positive dB.
'''

import numpy as np

from .._deps import gaussz
from ..results import Term
from ..assumptions import Assumptions, BEAM_GAUSSIAN, REGIME_NA, SPECTRUM_NA


def geometric_loss_db(range_m, w0, rx_diameter, wavelength=1550e-9,
                      obscuration_ratio=0.0):
    '''
    Geometric spreading loss of a Gaussian beam into a circular aperture.

    Parameters:
        range_m : float or ndarray
            Slant range from transmitter to receiver [m].
        w0 : float
            Transmit beam waist radius [m].
        rx_diameter : float
            Receiver aperture diameter [m].
        wavelength : float
            Wavelength [m].
        obscuration_ratio : float
            Ratio of central obscuration diameter to rx_diameter (e.g. the
            secondary mirror of a Cassegrain telescope). 0 = unobscured.

    Returns:
        float or ndarray
            Geometric loss [dB], positive.
    '''
    w_z = gaussz(w0, range_m, wavelength)
    a_rx = rx_diameter / 2
    b_obs = obscuration_ratio * a_rx
    eta = np.exp(-2 * b_obs**2 / w_z**2) - np.exp(-2 * a_rx**2 / w_z**2)
    return -10 * np.log10(eta)


def geometric_loss_term(scenario, geometry):
    '''
    Deterministic geometric-loss Term for a Scenario over a geometry.

    Parameters:
        scenario : Scenario
            Provides link.tx_waist_m, rx_diameter_m, rx_obscuration_ratio,
            wavelength_m.
        geometry : geometry object
            Provides slant_range_m (float or ndarray).

    Returns:
        Term
            category="geometric", mean_db = the loss over the geometry.
    '''
    link = scenario.link
    loss = geometric_loss_db(
        geometry.slant_range_m,
        link.tx_waist_m,
        link.rx_diameter_m,
        wavelength=link.wavelength_m,
        obscuration_ratio=link.rx_obscuration_ratio,
    )
    return Term(name="geometric spreading", category="geometric", mean_db=loss,
                assumptions=Assumptions(
                    beam_type=BEAM_GAUSSIAN,
                    turbulence_regime=REGIME_NA,
                    spectrum=SPECTRUM_NA,
                    validity="Far-field Gaussian beam into a circular aperture. "
                             "Paraxial. No turbulence.",
                ))


if __name__ == '__main__':
    from ..geometry import CircularOrbit
    from ..scenario import Scenario, Link

    geom = CircularOrbit(altitude_m=550e3, elevation_deg=[30, 60, 90])

    # pure formula: known limits
    r = geom.slant_range_m
    assert np.all(geometric_loss_db(r, 0.02, 1e6) < 1e-6)     # rx >> beam -> ~0 dB
    assert np.all(geometric_loss_db(r, 0.02, 1e-6) > 100)     # rx <<< beam -> huge loss
    assert np.all(geometric_loss_db(r, 0.02, 0.7, obscuration_ratio=0.3)
                  > geometric_loss_db(r, 0.02, 0.7))          # obscuration adds loss

    # Term path
    scn = Scenario(link=Link(tx_waist_m=0.035, rx_diameter_m=0.7,
                             rx_obscuration_ratio=0.3, wavelength_m=1550e-9))
    term = geometric_loss_term(scn, geom)
    assert term.category == "geometric"
    assert np.shape(term.mean_db) == (3,)
    assert term.assumptions is not None

    for el, loss in zip(geom.elevation_deg, np.atleast_1d(term.mean_db)):
        print(f"elevation {el:>4.0f} deg  ->  geometric loss {loss:6.2f} dB")
    print("self-check passed")
