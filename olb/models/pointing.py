'''
Pointing / tracking jitter fade for an FSO link.

The system aims a Gaussian transmit beam (1/e^2 radius w_z at the receiver
range) with 2-D Gaussian pointing jitter of 1-sigma angle sigma_theta. The
radial pointing displacement at the receiver is r = sigma_r * |unit 2-D
Gaussian|, with sigma_r = sigma_theta * range. For the small-aperture, on-axis
Gaussian approximation, the collected-power fraction against boresight is

    h(r) = exp(-2 r^2 / w_z^2)   ->   loss_db = (20/ln10) * r^2 / w_z^2

The two jitter axes are i.i.d. Gaussian, so r^2 is exponential. The loss in dB
then has an exponential distribution:

    loss_db ~ Exponential(mean = (20/ln10) * 2 * sigma_r^2 / w_z^2)

This gives the closed-form quantile and sampler below. Assumption: small receive
aperture, on-axis Gaussian beam (no aperture averaging of the fade).
'''

import numpy as np

from .._deps import gaussz
from ..results import Term
from ..assumptions import Assumptions, BEAM_GAUSSIAN, REGIME_NA, SPECTRUM_NA

_K = 20.0 / np.log(10.0)   # dB per (r^2 / w_z^2)


def pointing_loss_mean_db(range_m, w0, sigma_theta_rad, wavelength=1550e-9):
    '''
    Expected pointing-jitter loss E[loss].

    Parameters:
        range_m : float or ndarray
            Slant range from transmitter to receiver [m].
        w0 : float
            Transmit beam waist radius [m].
        sigma_theta_rad : float
            1-sigma pointing (tracking) jitter angle [rad].
        wavelength : float
            Wavelength [m].

    Returns:
        float or ndarray
            Mean pointing loss [dB], positive.
    '''
    w_z = gaussz(w0, range_m, wavelength)
    sigma_r = sigma_theta_rad * np.asarray(range_m)
    return _K * 2.0 * sigma_r**2 / w_z**2


def pointing_loss_term(scenario, geometry):
    '''
    Pointing-jitter fade Term for a Scenario over a geometry.

    For zero jitter the term is deterministic 0 dB. For nonzero jitter the loss
    has an exponential distribution, so the Term gives all three views (mean,
    quantile, sampler).

    Parameters:
        scenario : Scenario
            Provides link.tx_waist_m, wavelength_m, pointing_jitter_rad.
        geometry : geometry object
            Provides slant_range_m (float or ndarray).

    Returns:
        Term
            category="pointing".
    '''
    link = scenario.link
    range_m = geometry.slant_range_m
    sigma_theta = link.pointing_jitter_rad

    def _assumptions():
        return Assumptions(
            beam_type=BEAM_GAUSSIAN,
            turbulence_regime=REGIME_NA,
            spectrum=SPECTRUM_NA,
            validity="Small receive aperture relative to the beam. On-axis "
                     "Gaussian beam. 2-D Gaussian jitter.",
        )

    if sigma_theta == 0:
        return Term(name="pointing jitter", category="pointing", mean_db=0.0,
                    note="no jitter", assumptions=_assumptions())

    mean = pointing_loss_mean_db(range_m, link.tx_waist_m, sigma_theta,
                                 wavelength=link.wavelength_m)
    shape = np.shape(mean)

    def quantile(p):
        return -mean * np.log(1.0 - p)   # inverse exponential CDF

    def sampler(n, rng):
        return rng.exponential(scale=mean, size=(n, *shape))

    a = _assumptions()
    w_z = gaussz(link.tx_waist_m, range_m, link.wavelength_m)
    a_rx = link.rx_diameter_m / 2.0
    if np.any(a_rx > 0.5 * w_z):
        a.flag("Receive aperture is not small relative to the beam; the on-axis "
               "approximation is weak.")

    return Term(name="pointing jitter", category="pointing",
                mean_db=mean, sampler=sampler, quantile=quantile,
                assumptions=a)


if __name__ == '__main__':
    from ..geometry import CircularOrbit
    from ..scenario import Scenario, Link

    geom = CircularOrbit(altitude_m=550e3, elevation_deg=[30, 60, 90])

    # zero jitter -> deterministic 0 dB term, no sampler/quantile
    z = pointing_loss_term(Scenario(link=Link(pointing_jitter_rad=0.0)), geom)
    assert z.mean_db == 0.0 and not z.stochastic and z.quantile is None
    assert z.assumptions is not None

    # larger jitter -> larger mean loss
    r = geom.slant_range_m
    m_small = pointing_loss_mean_db(r, 0.035, 5e-6)
    m_big = pointing_loss_mean_db(r, 0.035, 20e-6)
    assert np.all(m_big > m_small)

    # statistical term: sampled mean ~= mean_loss_db, quantile(0.99) > mean
    scn = Scenario(link=Link(tx_waist_m=0.035, wavelength_m=1550e-9,
                             pointing_jitter_rad=10e-6))
    term = pointing_loss_term(scn, geom)
    assert term.stochastic and term.quantile is not None
    assert term.assumptions is not None
    rng = np.random.default_rng(0)
    draws = term.sample_db(200_000, rng)
    assert draws.shape == (200_000, 3)
    rel = np.abs(draws.mean(axis=0) - term.mean_db) / term.mean_db
    assert np.all(rel < 0.02), rel
    assert np.all(term.quantile_db(0.99) > term.mean_db)

    for el, mu, q99 in zip(geom.elevation_deg, np.atleast_1d(term.mean_db),
                           np.atleast_1d(term.quantile_db(0.99))):
        print(f"elevation {el:>4.0f} deg  ->  mean {mu:6.3f} dB  "
              f"99% fade {q99:6.3f} dB")
    print("self-check passed")
