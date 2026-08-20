'''
Retroreflected budget assembly: up-leg Terms plus down-leg Terms.

The ground station launches a beam up. The satellite retroreflector captures the
power. The reflected wavefront is flat at the satellite, so the return is a plane
wave that carries the captured power back to the ground receiver. The two legs
use independent turbulence. This module concatenates the up-leg Terms and the
down-leg Terms and adds a fixed retro-reflection Term.
'''

from dataclasses import replace

import numpy as np

from ..results import Budget, Term
from ..assumptions import (Assumptions, BEAM_PLANE_WAVE, REGIME_NA, SPECTRUM_NA)
from ..models.geometric import geometric_loss_term
from ..models.transmittance import atmospheric_loss_term, DEFAULT_TAU_ZENITH
from ..models.pointing import pointing_loss_term
from ..turbulence.profiles import default_cn2_profile
from .uplink import uplink_turbulence_term
from .downlink import downlink_scintillation_term


def retro_budget(scenario, geometry, *, turbulence=True, tau_zenith=None,
                 n_samples=3000, cn2_profile=None, retro_loss_db=0.0):
    '''
    Assemble the retroreflected budget: up-leg Terms plus down-leg Terms.

    The ground station launches a beam up. The satellite retroreflector captures
    the power. The reflected wavefront is flat at the satellite, so the return is
    a plane wave that carries the captured power back to the ground receiver. The
    two legs use independent turbulence. The retroreflector aperture is the
    hinge: it is the up-leg receive aperture and the down-leg transmit aperture.
    The losses are dB, so the Terms add.

    Parameters:
        scenario : Scenario
            The link case. link.retro_aperture_m must not be None.
        geometry : CircularOrbit or TLEPass
            The link geometry.
        turbulence : bool
            Add the up-leg coupled-flux turbulence Term when true.
        tau_zenith : float, optional
            Zenith optical depth. Defaults to transmittance.DEFAULT_TAU_ZENITH.
        n_samples : int
            Monte Carlo draws for the turbulence Term mean estimate.
        cn2_profile : numpy.ndarray, optional
            Explicit zenith Cn2 profile. Defaults to default_cn2_profile.
        retro_loss_db : float
            Fixed loss of the retroreflection [dB].

    Returns:
        Budget
            The budget with the original scenario set.
    '''
    link = scenario.link
    if link.retro_aperture_m is None:
        raise ValueError(
            "retro_budget needs the retroreflector aperture. "
            "Set scenario.link.retro_aperture_m."
        )
    retro_aperture_m = link.retro_aperture_m
    tau = DEFAULT_TAU_ZENITH if tau_zenith is None else tau_zenith

    # The retro aperture is the uplink receiver. A corner-cube retro has no
    # central obscuration.
    up_link = replace(link, direction="uplink",
                      rx_diameter_m=retro_aperture_m,
                      rx_obscuration_ratio=0.0)
    # The retro is the plane-wave transmitter. The Gaussian-equivalent waist is
    # half the aperture diameter. The ground stays the receiver. The retro is
    # passive, so there is no active pointing jitter on the return.
    down_link = replace(link, direction="downlink",
                        tx_waist_m=retro_aperture_m / 2.0,
                        pointing_jitter_rad=0.0)
    up_scn = replace(scenario, link=up_link)
    down_scn = replace(scenario, link=down_link)

    if cn2_profile is None:
        cn2_profile = default_cn2_profile(scenario.site)

    up_terms = [
        geometric_loss_term(up_scn, geometry),
        atmospheric_loss_term(up_scn, geometry, tau_zenith=tau),
        pointing_loss_term(up_scn, geometry),
    ]
    if turbulence:
        up_terms.append(uplink_turbulence_term(up_scn, geometry, n_samples=n_samples,
                                               cn2_profile=cn2_profile))
    for t in up_terms:
        t.name = "uplink " + t.name

    down_terms = [
        geometric_loss_term(down_scn, geometry),
        atmospheric_loss_term(down_scn, geometry, tau_zenith=tau),
        downlink_scintillation_term(down_scn, geometry),
    ]
    for t in down_terms:
        t.name = "downlink " + t.name

    retro_term = Term(
        name="retro reflection", category="system", mean_db=retro_loss_db,
        assumptions=Assumptions(
            beam_type=BEAM_PLANE_WAVE, turbulence_regime=REGIME_NA,
            spectrum=SPECTRUM_NA,
            validity="The reflected wavefront is flat at the satellite; the "
                     "return is a plane wave. The up and down legs use "
                     "independent turbulence (no reciprocity). The "
                     "retroreflector aperture is modelled as a Gaussian waist "
                     "of half the aperture diameter. The model does not include "
                     "velocity aberration or point-ahead loss on the return.",
        ),
    )

    return Budget(up_terms + down_terms + [retro_term], scenario=scenario)


if __name__ == '__main__':
    from ..scenario import Scenario, Link
    from ..geometry import CircularOrbit

    retro_scn = Scenario(
        link=Link(wavelength_m=1550e-9, tx_waist_m=0.06, retro_aperture_m=0.05,
                  rx_diameter_m=0.7, rx_obscuration_ratio=0.3,
                  pointing_jitter_rad=2e-6, tx_power_dbm=40,
                  rx_sensitivity_dbm=-50),
        altitude_m=1500e3,
    )
    retro_geom = CircularOrbit(altitude_m=1500e3, elevation_deg=45.0)

    no_aperture = replace(retro_scn, link=replace(retro_scn.link,
                                                  retro_aperture_m=None))
    try:
        retro_budget(no_aperture, retro_geom)
        raise AssertionError("retro_budget must reject a None retro aperture")
    except ValueError:
        pass

    retro = retro_budget(retro_scn, retro_geom)
    assert retro.to_frame().shape[0] == 8, retro.to_frame().shape
    frame = retro.to_frame()
    n_atmos = (frame["category"] == "atmospheric").sum()
    assert n_atmos == 2, n_atmos
    # The leg prefix must give "downlink scintillation", not the old double name.
    assert "downlink scintillation" in [t.name for t in retro.terms], \
        [t.name for t in retro.terms]
    retro_mc = retro.monte_carlo(2000, rng=np.random.default_rng(0),
                                 availabilities=(0.99,))
    retro_margin = retro_mc["margin_db"][0.99]
    assert np.isfinite(retro_margin), retro_margin
    af = retro.assumptions_frame()
    retro_row = af[af["name"] == "retro reflection"]
    assert not retro_row.empty, "retro reflection row missing"
    assert "plane wave" in retro_row.iloc[0]["validity"], retro_row.iloc[0]

    print("retro budget terms:")
    print(retro.to_frame().to_string(index=False))
    print(f"\nretro 45 deg 99% margin: {retro_margin:.2f} dB")
    print("self-check passed")
