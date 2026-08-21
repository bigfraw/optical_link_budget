'''
Transmit Gaussian efficiency: the truncation loss at the launch aperture.

A real transmitter sends a Gaussian beam through a finite, possibly obscured,
circular aperture. The aperture clips the wings of the Gaussian and the central
obscuration blocks the middle. This is a fixed hardware loss. It does not depend
on range.

The efficiency is the on-axis far-field gain of the truncated beam, referenced
to the gain of the untruncated source Gaussian:

    alpha = a / w_T = (tx_aperture_m / 2) / tx_waist_m
    eta   = (exp(-alpha^2) - exp(-alpha^2 * Cr^2))^2

where a is the aperture radius, w_T is the Gaussian waist (1/e^2 radius) at the
aperture, and Cr is the linear central-obscuration ratio (obscuration diameter /
aperture diameter). The loss is -10*log10(eta), positive dB.

eta goes to 1 when the aperture is much wider than the beam (no truncation) and
to 0 when the aperture is much smaller than the beam (the wings carry the power).
The classic optimal truncation is near alpha = 1.12.

This is the corrected antenna-gain form. It has NO 2/alpha^2 prefactor; that
prefactor double-counts a normalisation the untruncated-source reference already
carries. The removal is validated against a numerical Fraunhofer propagation of a
truncated Gaussian (tn2_kepler test_gauss_prop).

This module is direction-agnostic. For an uplink the transmitter is the ground
station; for a downlink it is the satellite.
'''

import numpy as np

from ..results import Term
from ..assumptions import Assumptions, BEAM_GAUSSIAN, REGIME_NA, SPECTRUM_NA


def gaussian_efficiency(alpha, obscuration_ratio=0.0):
    '''
    On-axis truncation efficiency of a Gaussian beam in a circular aperture.

    Parameters:
        alpha : float or ndarray
            Ratio of aperture radius to Gaussian waist, a / w_T.
        obscuration_ratio : float
            Linear central-obscuration ratio Cr (obscuration diameter / aperture
            diameter). 0 = unobscured.

    Returns:
        float or ndarray
            Efficiency eta in (0, 1].
    '''
    a2 = np.asarray(alpha, dtype=float) ** 2
    return (np.exp(-a2) - np.exp(-a2 * obscuration_ratio ** 2)) ** 2


def tx_efficiency_loss_db(tx_aperture_m, tx_waist_m, obscuration_ratio=0.0):
    '''
    Transmit truncation loss of a Gaussian beam at a circular aperture.

    Parameters:
        tx_aperture_m : float
            Transmit aperture diameter [m].
        tx_waist_m : float
            Gaussian waist (1/e^2 radius) at the aperture [m].
        obscuration_ratio : float
            Linear central-obscuration ratio Cr. 0 = unobscured.

    Returns:
        float
            Truncation loss [dB], positive.
    '''
    alpha = (tx_aperture_m / 2) / tx_waist_m
    eta = gaussian_efficiency(alpha, obscuration_ratio)
    return -10 * np.log10(eta)


def tx_gaussian_efficiency_term(scenario, geometry=None):
    '''
    Transmit Gaussian-efficiency (truncation) Term for a Scenario.

    Range-independent, so the geometry is not read. The signature keeps the
    common f(scenario, geometry) -> Term shape.

    Parameters:
        scenario : Scenario
            Reads the transmit terminal aperture_m and obscuration_ratio and its
            Transmitter waist_m. See olb.scenario.
        geometry : object, optional
            Unused. Present for a uniform model signature.

    Returns:
        Term
            name="transmit Gaussian efficiency", category="system".
    '''
    tx = scenario.tx_terminal
    aperture_m = tx.aperture_m
    waist_m = tx.transmitter.waist_m
    obscuration_ratio = tx.obscuration_ratio
    alpha = (aperture_m / 2) / waist_m
    loss = tx_efficiency_loss_db(aperture_m, waist_m, obscuration_ratio)
    return Term(
        name="transmit Gaussian efficiency",
        category="system",
        mean_db=float(loss),
        note=f"aperture truncation, alpha={alpha:.3f}, "
             f"Cr={obscuration_ratio:g}",
        meta={"alpha": float(alpha),
              "eta": float(gaussian_efficiency(alpha, obscuration_ratio))},
        assumptions=Assumptions(
            beam_type=BEAM_GAUSSIAN,
            turbulence_regime=REGIME_NA,
            spectrum=SPECTRUM_NA,
            validity="On-axis far-field gain of a truncated Gaussian, referenced "
                     "to the untruncated source. Paraxial. No turbulence.",
        ),
    )


if __name__ == '__main__':
    from ..scenario import Scenario
    from ..terminal import Terminal, Transmitter

    # Limits: wide aperture -> no loss; narrow aperture -> large loss.
    assert gaussian_efficiency(10.0) > 1 - 1e-6              # aperture >> beam -> eta ~ 1
    assert tx_efficiency_loss_db(2.0, 0.1) < 1e-3           # aperture >> beam -> ~0 dB
    assert tx_efficiency_loss_db(0.01, 0.1) > 20            # aperture << beam -> big loss

    # Obscuration only adds loss.
    assert (tx_efficiency_loss_db(0.15, 0.12, obscuration_ratio=0.3)
            > tx_efficiency_loss_db(0.15, 0.12))

    # Corrected form: NO 2/alpha^2 prefactor. eta is bounded by 1, so the loss is
    # never a gain. The stale form (with 2/alpha^2) gives eta > 1 (negative loss)
    # for a small alpha; the corrected form does not.
    assert np.all(gaussian_efficiency(np.array([0.1, 0.3, 0.5, 1.0, 2.0])) <= 1.0)

    # TN-2 launch: Da=150 mm, w_T=0.8*Da, Cr=0.3 -> alpha=0.625.
    tn2 = tx_efficiency_loss_db(0.150, 0.8 * 0.150, obscuration_ratio=0.3)
    assert 10.0 < tn2 < 12.0, tn2                            # ~10.8 dB

    # Term path: uplink -> tx=ground carries the aperture, obscuration, waist.
    scn = Scenario(
        ground=Terminal(aperture_m=0.150, obscuration_ratio=0.3,
                        transmitter=Transmitter(waist_m=0.12)),
        space=Terminal(aperture_m=0.05),
        direction="uplink")
    term = tx_gaussian_efficiency_term(scn)
    assert term.category == "system"
    assert np.isscalar(term.mean_db) and term.mean_db > 0
    assert 0.0 < term.meta["eta"] <= 1.0
    assert term.assumptions is not None

    print(f"TN-2 transmit truncation loss: {tn2:.2f} dB")
    print(f"Term: {term.name}  {term.mean_db:.2f} dB  ({term.note})")
    print("self-check passed")
