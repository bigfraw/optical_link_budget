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

from ..beam import free_space_radius
from ..results import Term
from ..assumptions import Assumptions, BEAM_GAUSSIAN, REGIME_NA, SPECTRUM_NA


def geometric_loss_db(range_m, w0, rx_diameter, wavelength=1550e-9,
                      obscuration_ratio=0.0, divergence_rad=None):
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
        divergence_rad : float, optional
            Transmit far-field half-angle divergence [rad]. None means
            collimated (the diffraction limit). A larger value spreads the beam
            more and adds loss.

    Returns:
        float or ndarray
            Geometric loss [dB], positive.
    '''
    w_z = free_space_radius(w0, range_m, divergence_rad, wavelength)
    a_rx = rx_diameter / 2
    b_obs = obscuration_ratio * a_rx
    eta = np.exp(-2 * b_obs**2 / w_z**2) - np.exp(-2 * a_rx**2 / w_z**2)
    return -10 * np.log10(eta)


def geometric_loss_term(scenario, geometry):
    '''
    Deterministic geometric-loss Term for a scenario over a geometry.

    Parameters:
        scenario : SpaceScenario or TerrestrialScenario
            Reads the transmit terminal (waist, divergence, wavelength) and the
            receive terminal (aperture, obscuration). See olb.scenario.
        geometry : geometry object
            Provides slant_range_m (float or ndarray).

    Returns:
        Term
            category="geometric", mean_db = the loss over the geometry.
    '''
    tx = scenario.tx_terminal
    rx = scenario.rx_terminal
    loss = geometric_loss_db(
        geometry.slant_range_m,
        tx.transmitter.waist_m,
        rx.aperture_m,
        wavelength=tx.wavelength_m,
        obscuration_ratio=rx.obscuration_ratio,
        divergence_rad=tx.transmitter.divergence_rad,
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
    from ..scenario import SpaceScenario
    from ..terminal import Terminal, Transmitter

    geom = CircularOrbit(altitude_m=550e3, elevation_deg=[30, 60, 90])

    # pure formula: known limits
    r = geom.slant_range_m
    assert np.all(geometric_loss_db(r, 0.02, 1e6) < 1e-6)     # rx >> beam -> ~0 dB
    assert np.all(geometric_loss_db(r, 0.02, 1e-6) > 100)     # rx <<< beam -> huge loss
    assert np.all(geometric_loss_db(r, 0.02, 0.7, obscuration_ratio=0.3)
                  > geometric_loss_db(r, 0.02, 0.7))          # obscuration adds loss

    # A diverged beam spreads more, so it costs more geometric loss than a
    # collimated beam of the same w0.
    from ..units import w0_to_div
    theta_min = w0_to_div(0.02, 1550e-9)
    assert np.all(geometric_loss_db(r, 0.02, 0.08, divergence_rad=5 * theta_min)
                  > geometric_loss_db(r, 0.02, 0.08))         # divergence adds loss
    assert np.allclose(geometric_loss_db(r, 0.02, 0.08, divergence_rad=theta_min),
                       geometric_loss_db(r, 0.02, 0.08))      # limit == collimated

    # Term path: uplink -> tx=ground (waist), rx=space (aperture, obscuration).
    scn = SpaceScenario(
        ground=Terminal(aperture_m=0.3, wavelength_m=1550e-9,
                        transmitter=Transmitter(waist_m=0.035)),
        space=Terminal(aperture_m=0.7, obscuration_ratio=0.3, wavelength_m=1550e-9),
        direction="uplink")
    term = geometric_loss_term(scn, geom)
    assert term.category == "geometric"
    assert np.shape(term.mean_db) == (3,)
    assert term.assumptions is not None

    for el, loss in zip(geom.elevation_deg, np.atleast_1d(term.mean_db)):
        print(f"elevation {el:>4.0f} deg  ->  geometric loss {loss:6.2f} dB")
    print("self-check passed")
