'''
Retroreflected ground-to-space budget: retroreflection as a retransmission.

The crucial idea is that a retroreflector RE-TRANSMITS the beam. The ground
station launches a beam up. The satellite retroreflector captures the power over
its aperture. It then re-emits that captured power back down as a new beam. So a
retro link is an up-leg transmission followed by a down-leg transmission, with
the retro aperture as the hinge: it is the up-leg receive aperture and the
down-leg transmit aperture.

Because every Term is a dB loss and the losses add, the return power is the
retransmission chain:

    P_return = P_launch
             - (up-leg losses)     # fraction of launch power that hits the retro
             - retro_loss
             - (down-leg losses)   # fraction of re-emitted power that hits ground

The up-leg geometric Term already gives the fraction of the launch power that the
retro aperture captures. The down-leg geometric Term treats the retro as a fresh
Gaussian transmitter (waist = half the aperture diameter). So the chained dB sum
IS the retransmission model; no explicit power bookkeeping is needed.

SPACE ONLY. This retransmission picture holds for a ground-to-space link, where
the slant range is long, the return beam diverges far past the ground aperture,
and the two legs see independent turbulence. It does NOT hold for a short
terrestrial (horizontal-path) retro link, where the return does not fully
diverge and reciprocity couples the legs. That case needs a different module.
'''

from dataclasses import replace

import numpy as np

from ..results import Budget, Term
from ..assumptions import (Assumptions, BEAM_PLANE_WAVE, REGIME_NA, SPECTRUM_NA)
from ..models.geometric import geometric_loss_term
from ..models.gaussian_efficiency import (uniform_aperture_correction_db,
                                          tx_gaussian_efficiency_term)
from ..models.extinction import slant_extinction_term, DEFAULT_TAU_ZENITH
from ..models.pointing import pointing_loss_term
from ..terminal import Terminal, Transmitter
from ..turbulence.profiles import default_cn2_profile
from .uplink import uplink_turbulence_term, TX_TRUNCATION_MIN_DB
from .downlink import downlink_scintillation_term


def retro_space_budget(scenario, geometry, *, fidelity=1, turbulence=True,
                       tau_zenith=None, n_samples=3000, cn2_profile=None,
                       retro_loss_db=0.0, fast_params=None):
    '''
    Assemble the retroreflected ground-to-space budget as a retransmission.

    Retroreflection is modelled as a retransmission: the up-leg carries the
    launch power to the retro aperture, and the down-leg re-emits the captured
    power back to the ground receiver. The two legs use independent turbulence.
    The retroreflector aperture is the hinge: it is the up-leg receive aperture
    and the down-leg transmit aperture. The losses are dB, so the Terms add and
    the sum is the retransmission chain (see the module docstring).

    This is the SPACE model. It assumes a long slant range, a fully diverged
    return, and independent turbulence on the two legs. Do not use it for a short
    terrestrial retro link.

    Parameters:
        scenario : SpaceScenario
            The link case. The `space` Terminal is the passive retroreflector;
            its aperture_m is the retro aperture. The direction is "retro".
        geometry : CircularOrbit or TLEPass
            The link geometry.
        turbulence : bool
            Add the up-leg coupled-flux turbulence Term when true.
        tau_zenith : float, optional
            Zenith optical depth. Defaults to extinction.DEFAULT_TAU_ZENITH.
        n_samples : int
            Monte Carlo draws for the turbulence Term mean estimate.
        cn2_profile : numpy.ndarray, optional
            Explicit zenith Cn2 profile. Defaults to default_cn2_profile.
        retro_loss_db : float
            Fixed loss of the retroreflection [dB].
        fidelity : int
            The down-leg receive-coupling fidelity: 1 (the default, FAST modal
            overlap) or 0 (analytic mean-only). The UP-leg turbulence stays the
            coupled-flux Monte Carlo at either value (there is no analytic
            mean-only uncorrected uplink model, so the up-leg is fidelity 1
            regardless). fidelity=2 (wave optics) is NOT supported: the folded
            double pass shares its screens (the two legs are correlated), which
            needs its own design (see CLAUDE.md, the deferred folded double pass).
        fast_params : dict, optional
            Extra FAST parameters for the fidelity-1 down-leg coupling.

    Returns:
        Budget
            The budget with the original scenario set.

    Raises:
        ValueError
            If fidelity is not 0 or 1 (fidelity=2 is deferred for retro).
    '''
    if fidelity == 2:
        raise ValueError(
            "fidelity=2 (wave optics) is not supported for a retro link. The "
            "folded double pass shares its screens, so the two legs are "
            "correlated; that needs its own design. Use fidelity=0 or 1."
        )
    if fidelity not in (0, 1):
        raise ValueError(f"fidelity must be 0 or 1 for retro, got {fidelity!r}.")
    smf_fidelity = "fast" if fidelity == 1 else "mean"
    retro_aperture_m = scenario.space.aperture_m
    wavelength = scenario.space.wavelength_m
    tau = DEFAULT_TAU_ZENITH if tau_zenith is None else tau_zenith

    # The retro is the uplink receiver. A corner-cube retro has no central
    # obscuration. The up-leg keeps the ground transmit terminal.
    retro_rx = Terminal(aperture_m=retro_aperture_m, obscuration_ratio=0.0,
                        wavelength_m=wavelength)
    up_scn = replace(scenario, direction="uplink", space=retro_rx)
    # The retro re-transmits the captured power: it is the plane-wave transmitter
    # on the return. The Gaussian-equivalent waist is half the aperture diameter.
    # The retro is passive, so there is no active pointing jitter. The ground
    # stays the receiver.
    retro_tx = Terminal(aperture_m=retro_aperture_m, obscuration_ratio=0.0,
                        wavelength_m=wavelength,
                        transmitter=Transmitter(waist_m=retro_aperture_m / 2.0))
    down_scn = replace(scenario, direction="downlink", space=retro_tx)

    if cn2_profile is None:
        cn2_profile = default_cn2_profile(scenario.channel.site)

    up_terms = [
        geometric_loss_term(up_scn, geometry),
        slant_extinction_term(up_scn, geometry, tau_zenith=tau),
    ]
    # Up-leg pointing jitter folds into the coupled-flux turbulence Term (it
    # shares the beam-wander displacement), so add the standalone pointing Term
    # only when the turbulence Term is off. Adding both double-counts the jitter.
    # See olb.links.uplink.uplink_budget.
    if not turbulence:
        up_terms.append(pointing_loss_term(up_scn, geometry))
    # The launch aperture truncates the up-leg beam, exactly as uplink_budget
    # does. Opt-in: it fires only when the ground transmitter truncates the beam
    # by more than TX_TRUNCATION_MIN_DB. A bistatic ground reads its beam-director
    # aperture through the Transmitter override, not the receive-telescope
    # aperture (see olb.terminal.Transmitter).
    up_tx = up_scn.tx_terminal
    if up_tx.transmitter is not None:
        eff = tx_gaussian_efficiency_term(up_scn, geometry)
        if eff.mean_db > TX_TRUNCATION_MIN_DB:
            up_terms.append(eff)
    if turbulence:
        up_terms.append(uplink_turbulence_term(up_scn, geometry, n_samples=n_samples,
                                               cn2_profile=cn2_profile))
    for t in up_terms:
        t.name = "uplink " + t.name

    # The retro reflects the flat wavefront that fills its aperture, so the
    # return is a uniformly illuminated (top-hat) aperture, not a Gaussian. The
    # geometric Term models a Gaussian of waist aperture/2, which over-states the
    # on-axis gain by 2/(1-Cr^2). This fixed Term converts it to the top-hat.
    tophat_db = uniform_aperture_correction_db(retro_tx.obscuration_ratio)
    tophat_term = Term(
        name="top-hat correction", category="system", mean_db=tophat_db,
        note="retro return is a uniform aperture, not a Gaussian(waist=D/2)",
        assumptions=Assumptions(
            beam_type=BEAM_PLANE_WAVE, turbulence_regime=REGIME_NA,
            spectrum=SPECTRUM_NA,
            validity="Converts the Gaussian(waist=aperture/2) geometric model to "
                     "a uniformly illuminated (top-hat) aperture, via the "
                     "Gaussian aperture-illumination efficiency. Far-field, "
                     "on-axis (receive aperture much smaller than the return "
                     "lobe).",
        ),
    )
    down_terms = [
        geometric_loss_term(down_scn, geometry),
        slant_extinction_term(down_scn, geometry, tau_zenith=tau),
        tophat_term,
    ]
    # The return-leg receive term follows the ground receiver, exactly as
    # downlink_budget does: when the ground has a detector, the receive-coupling
    # Term owns the receive-side physics (SMF adds the fibre-coupling loss; an
    # Aperture reproduces the plain scintillation). Without a detector, fall back
    # to the standalone plane-wave scintillation.
    rx = down_scn.rx_terminal
    if rx is not None and rx.detector is not None:
        # Import here to break the downlink <-> coupling import cycle.
        from ..models.coupling import downlink_coupling_term
        down_terms.append(downlink_coupling_term(down_scn, geometry, n_samples=n_samples,
                                           smf_fidelity=smf_fidelity,
                                           fast_params=fast_params))
    else:
        down_terms.append(downlink_scintillation_term(down_scn, geometry))
    for t in down_terms:
        t.name = "downlink " + t.name

    retro_term = Term(
        name="retro reflection", category="system", mean_db=retro_loss_db,
        assumptions=Assumptions(
            beam_type=BEAM_PLANE_WAVE, turbulence_regime=REGIME_NA,
            spectrum=SPECTRUM_NA,
            validity="Retroreflection is modelled as a retransmission: the retro "
                     "captures the up-leg power and re-emits it down. This holds "
                     "for a ground-to-space link (long range, fully diverged "
                     "return, independent turbulence on the two legs), not for a "
                     "short terrestrial link. The reflected wavefront is flat at "
                     "the satellite; the return is a plane wave. The "
                     "retroreflector aperture is modelled as a Gaussian waist of "
                     "half the aperture diameter. The model does not include "
                     "velocity aberration or point-ahead loss on the return.",
        ),
    )

    return Budget(up_terms + down_terms + [retro_term], scenario=scenario)


if __name__ == '__main__':
    from ..scenario import SpaceScenario, Channel
    from ..terminal import Terminal, Transmitter, Aperture
    from ..geometry import CircularOrbit

    # The ground terminal transmits up and receives the return. The space
    # terminal is the passive retroreflector (aperture only). The ground is
    # bistatic: a small beam director (0.15 m) transmits the up-leg, so the launch
    # truncation reads the director aperture, not the 0.7 m receive telescope.
    retro_scn = SpaceScenario(
        ground=Terminal(aperture_m=0.7, obscuration_ratio=0.3, wavelength_m=1550e-9,
                        pointing_jitter_rad=0e-6,
                        transmitter=Transmitter(waist_m=0.06, power_dbm=40,
                                                aperture_m=0.15, obscuration_ratio=0.0),
                        detector=Aperture(sensitivity_dbm=-50)),
        space=Terminal(aperture_m=0.05, wavelength_m=1550e-9),
        direction="retro", channel=Channel(altitude_m=1500e3),
    )
    retro_geom = CircularOrbit(altitude_m=1500e3, elevation_deg=45.0)

    retro = retro_space_budget(retro_scn, retro_geom)
    # 9 terms: the up-leg carries the opt-in launch-truncation term because the
    # 0.15 m beam director truncates the 0.06 m waist beam. With turbulence on
    # there is NO standalone up-leg pointing Term (the jitter folds into the
    # coupled-flux turbulence Term).
    assert retro.to_frame().shape[0] == 9, retro.to_frame().shape
    assert not any(t.category == "pointing" for t in retro.terms)
    names = [t.name for t in retro.terms]
    assert "uplink transmit Gaussian efficiency" in names, names
    # The truncation reads the director aperture (0.15 m), not the receive
    # telescope (0.7 m): alpha = (0.15/2)/0.06 = 1.25, an unobscured director.
    eff = next(t for t in retro.terms if t.name == "uplink transmit Gaussian efficiency")
    assert abs(eff.meta["alpha"] - 1.25) < 1e-9, eff.meta["alpha"]
    frame = retro.to_frame()
    n_atmos = (frame["category"] == "atmospheric").sum()
    assert n_atmos == 2, n_atmos
    # The return leg carries the top-hat correction (+3.01 dB, unobscured retro).
    names = [t.name for t in retro.terms]
    assert "downlink top-hat correction" in names, names
    tophat = next(t for t in retro.terms if t.name == "downlink top-hat correction")
    assert abs(tophat.mean_db - 3.0103) < 1e-3, tophat.mean_db
    # The ground has an Aperture detector, so the return leg carries the
    # receive-coupling Term (leg-prefixed), not the standalone scintillation.
    assert "downlink receive coupling (aperture)" in [t.name for t in retro.terms], \
        [t.name for t in retro.terms]
    retro_mc = retro.monte_carlo(2000, rng=np.random.default_rng(0),
                                 availabilities=(0.99,))
    retro_margin = retro_mc["margin_db"][0.99]
    assert np.isfinite(retro_margin), retro_margin
    af = retro.assumptions_frame()
    retro_row = af[af["name"] == "retro reflection"]
    assert not retro_row.empty, "retro reflection row missing"
    assert "plane wave" in retro_row.iloc[0]["validity"], retro_row.iloc[0]
    assert "retransmission" in retro_row.iloc[0]["validity"], retro_row.iloc[0]

    print("retro (space) budget terms:")
    print(retro.to_frame().to_string(index=False))
    print(f"\nretro 45 deg 99% margin: {retro_margin:.2f} dB")
    print("self-check passed")
