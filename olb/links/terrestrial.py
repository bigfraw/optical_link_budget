'''
Terrestrial (horizontal-path) Terms and budget assembly.

This module builds the fidelity-zero terrestrial link budget: a ground-to-ground
horizontal path. It reuses the direction-agnostic Terms (geometric spreading,
pointing jitter, transmit truncation) and adds the horizontal Beer-Lambert
extinction Term. See olb.models.

A horizontal path differs from a space link in two ways that matter here. The
range is a constant path length, not a slant range that changes with elevation.
And the Gaussian-beam properties (waist, divergence, curvature) steer the
turbulence result strongly, because the whole path sits in the near field of a
finite beam. So the horizontal scintillation is NOT the plane-wave slant-path
model that the downlink uses. It needs the Gaussian-beam analytic forms.

That analytic scintillation Term is a RESERVED SLOT. It is not implemented here.
The terrestrial_scintillation_term below raises NotImplementedError and names
the Andrews equations that it needs. Until it is filled, build the budget with
scintillation=False (the default). The budget then holds the deterministic
Terms, which are exact and direction-agnostic.
'''

from ..results import Budget
from ..models.geometric import geometric_loss_term
from ..models.transmittance import terrestrial_extinction_term
from ..models.pointing import pointing_loss_term
from ..models.gaussian_efficiency import tx_gaussian_efficiency_term

# Below this launch-truncation loss the beam is an untruncated Gaussian, so the
# transmit Gaussian-efficiency term is skipped [dB]. Matches olb.links.uplink.
TX_TRUNCATION_MIN_DB = 1e-2


def terrestrial_scintillation_term(scenario, geometry):
    '''
    Reserved slot: the analytic Gaussian-beam horizontal scintillation Term.

    NOT IMPLEMENTED. A horizontal-path scintillation model is steered by the
    Gaussian-beam parameters at the receiver, so it is not the plane-wave
    slant-path model in olb.turbulence.scintillation. Do not substitute an
    invented formula. This slot needs the following analytic forms, all for a
    horizontal path of constant Cn2 (scenario.channel.cn2) and length L
    (geometry.slant_range_m), from Andrews and Phillips, Laser Beam Propagation
    through Random Media, 2nd ed. (2005):

      1. The plane-wave Rytov variance
             sigma_R^2 = 1.23 * Cn2 * k^(7/6) * L^(11/6)
         as the weak-fluctuation reference (Ch. 8).
      2. The Gaussian-beam on-axis (longitudinal) scintillation index
             sigma_I^2(0, L)
         in terms of the receiver-plane beam parameters Theta (and
         Theta_bar = 1 - Theta) and Lambda (Ch. 8). This is the term the beam
         divergence and curvature steer.
      3. The horizontal aperture-averaging factor for a circular receive
         aperture (Ch. 10), to reduce the point index to the flux index.

    With sigma_I^2 (and its aperture-averaged flux value) in hand, the lognormal
    fade faces (mean_db, quantile, sampler) follow the same closed form as the
    downlink lognormal Term (see olb.links.downlink._lognormal_term). Reuse that
    machinery; only the sigma^2 computation is missing.

    Raises:
        NotImplementedError
            Always. Fill this slot with the Andrews Gaussian-beam forms above.
    '''
    raise NotImplementedError(
        "terrestrial_scintillation_term is a reserved slot. A horizontal-path "
        "scintillation model needs the Andrews Gaussian-beam scintillation index "
        "sigma_I^2(0, L) (Andrews & Phillips 2nd ed. 2005, Ch. 8), the Rytov "
        "variance sigma_R^2 = 1.23*Cn2*k^(7/6)*L^(11/6), and the horizontal "
        "aperture-averaging factor (Ch. 10). See the docstring. Build the budget "
        "with scintillation=False until this slot is filled."
    )


def terrestrial_budget(scenario, geometry, *, scintillation=False):
    '''
    Assemble the terrestrial budget: geometric, horizontal extinction, pointing.

    The deterministic Terms are exact and direction-agnostic. The horizontal
    scintillation Term is a reserved slot (see terrestrial_scintillation_term),
    so scintillation defaults to False. Set it True only once that slot is
    filled; it raises NotImplementedError otherwise.

    Parameters:
        scenario : TerrestrialScenario
            A terrestrial link case. tx = near end, rx = far end. Its
            TerrestrialChannel carries path_length_m, attenuation_db_per_km, cn2.
        geometry : HorizontalPath
            The horizontal path (reads slant_range_m = path length).
        scintillation : bool
            Add the (pending) analytic scintillation Term when True. Raises
            NotImplementedError until the reserved slot is filled.

    Returns:
        Budget
            The budget with the scenario set.
    '''
    terms = [
        geometric_loss_term(scenario, geometry),
        terrestrial_extinction_term(scenario, geometry),
        pointing_loss_term(scenario, geometry),
    ]
    # The transmit Gaussian-efficiency (launch truncation) Term is opt-in. It
    # fires only when the transmit terminal has a Transmitter whose launch
    # aperture truncates the beam by more than TX_TRUNCATION_MIN_DB. A wide
    # aperture leaves the beam an untruncated Gaussian, so the Term is skipped.
    # Same rule as olb.links.uplink.
    tx = scenario.tx_terminal
    if tx.transmitter is not None:
        eff = tx_gaussian_efficiency_term(scenario, geometry)
        if eff.mean_db > TX_TRUNCATION_MIN_DB:
            terms.append(eff)
    # An SMF detector on the far terminal adds the fidelity-0 (mean-only) fibre-
    # coupling loss, using the horizontal Gaussian-beam r0 and the compensation
    # stack (tip-tilt, AO). It is MEAN-ONLY, so it locks the budget to fidelity 0
    # and the budget then refuses a fade margin. An Aperture (bucket) detector is
    # phase-insensitive; its turbulence penalty is the (reserved) scintillation
    # Term, so it is not added here. See olb.models.coupling.
    from ..terminal import SMF
    rx = scenario.rx_terminal
    if isinstance(rx.detector, SMF):
        # Lazy import breaks the terrestrial <-> coupling import cycle.
        from ..models.coupling import terrestrial_smf_coupling_term
        terms.append(terrestrial_smf_coupling_term(scenario, geometry))
    if scintillation:
        terms.append(terrestrial_scintillation_term(scenario, geometry))
    return Budget(terms, scenario=scenario)


if __name__ == '__main__':
    import numpy as np

    from ..scenario import TerrestrialScenario, TerrestrialChannel
    from ..geometry import HorizontalPath
    from ..terminal import Terminal, Transmitter, Aperture, SMF, TipTilt, AO

    def _terr(w0, L, *, divergence=None, power=None, jitter=0.0,
              near_aperture=0.1, near_obscuration=0.0, far_aperture=0.1,
              attenuation=0.5, sensitivity=None):
        '''Build a TerrestrialScenario: tx = near, rx = far.'''
        detector = None if sensitivity is None else Aperture(sensitivity_dbm=sensitivity)
        return TerrestrialScenario(
            near=Terminal(aperture_m=near_aperture, obscuration_ratio=near_obscuration,
                          wavelength_m=1550e-9, pointing_jitter_rad=jitter,
                          transmitter=Transmitter(waist_m=w0, power_dbm=power,
                                                  divergence_rad=divergence)),
            far=Terminal(aperture_m=far_aperture, wavelength_m=1550e-9,
                         detector=detector),
            channel=TerrestrialChannel(path_length_m=L, attenuation_db_per_km=attenuation))

    # A clean 5 km link. A wide near aperture (0.3 m for a 0.02 m waist) leaves
    # the beam untruncated, so the launch-truncation Term does not fire.
    scn = _terr(0.02, 5e3, power=30, jitter=5e-6, sensitivity=-40,
                near_aperture=0.3)
    geom = HorizontalPath(5e3)
    budget = terrestrial_budget(scn, geom)
    names = [t.name for t in budget.terms]
    assert names == ["geometric spreading", "atmospheric extinction (horizontal)",
                     "pointing jitter"], names
    # The extinction Term is exact: 5 km * 0.5 dB/km = 2.5 dB.
    ext = next(t for t in budget.terms if t.category == "atmospheric")
    assert np.isclose(ext.mean_db, 2.5), ext.mean_db

    # A narrow near aperture (0.02 m for a 0.02 m waist) truncates the beam, so
    # the launch-truncation Term fires.
    scn_ap = _terr(0.02, 5e3, power=30, near_aperture=0.02, near_obscuration=0.2)
    budget_ap = terrestrial_budget(scn_ap, geom)
    assert "transmit Gaussian efficiency" in [t.name for t in budget_ap.terms]
    assert budget_ap.total_loss_db() > budget.total_loss_db()

    # A longer path costs more geometric spread AND more extinction.
    long_budget = terrestrial_budget(_terr(0.02, 10e3, near_aperture=0.3),
                                     HorizontalPath(10e3))
    assert long_budget.total_loss_db() > budget.total_loss_db()

    # The deterministic budget has an analytic fade (every Term has a quantile).
    fade = budget.fade_margin_db(0.99)
    assert np.isfinite(fade)

    # --- SMF fibre coupling (fidelity-0, mean-only) -------------------------
    import warnings as _warnings

    def _smf(compensation=None, near_aperture=0.3, far_aperture=0.2, w0=0.02,
             cn2=1e-14):
        scn = TerrestrialScenario(
            near=Terminal(aperture_m=near_aperture, wavelength_m=1550e-9,
                          transmitter=Transmitter(waist_m=w0, power_dbm=30)),
            far=Terminal(aperture_m=far_aperture, wavelength_m=1550e-9,
                         detector=SMF(sensitivity_dbm=-40),
                         compensation=compensation or []),
            channel=TerrestrialChannel(path_length_m=3e3, attenuation_db_per_km=0.5,
                                       cn2=cn2))
        return terrestrial_budget(scn, HorizontalPath(3e3))

    smf_budget = _smf()
    with _warnings.catch_warnings():
        _warnings.simplefilter("ignore")
        assert "receive coupling (SMF)" in [t.name for t in smf_budget.terms]
        coupling = next(t for t in smf_budget.terms if t.category == "coupling")
    # It is mean-only, so the budget is fidelity 0 and refuses a fade margin.
    assert coupling.mean_only and not smf_budget.provides_fade
    try:
        smf_budget.fade_margin_db(0.99)
    except ValueError as e:
        assert "fidelity-0" in str(e) and "mean-only" in str(e)
    else:
        raise AssertionError("a mean-only budget must refuse fade_margin_db")
    # The mean total loss is still reported (that is the fidelity-0 deliverable).
    assert np.isfinite(smf_budget.total_loss_db())
    # Monte Carlo reports the mean but suppresses the fade.
    with _warnings.catch_warnings():
        _warnings.simplefilter("ignore")
        mc = smf_budget.monte_carlo(2000, np.random.default_rng(0))
    assert mc["fade_available"] is False and mc["fade_db"] is None
    assert mc["margin_db"] is None and np.isfinite(mc["mean_loss_db"])

    # Tip-tilt, then full AO, each buys back coupling (less loss than none).
    with _warnings.catch_warnings():
        _warnings.simplefilter("ignore")
        loss_none = next(t for t in _smf().terms if t.category == "coupling").mean_db
        loss_tt = next(t for t in _smf([TipTilt()]).terms
                       if t.category == "coupling").mean_db
        loss_ao = next(t for t in _smf([TipTilt(), AO(200)]).terms
                       if t.category == "coupling").mean_db
    assert loss_ao < loss_tt < loss_none, (loss_ao, loss_tt, loss_none)

    # A bucket (Aperture) detector does NOT add a coupling Term, and the budget
    # keeps its analytic fade (no mean-only term).
    assert budget.provides_fade
    assert np.isfinite(budget.fade_margin_db(0.99))

    # A path-length sweep broadcasts over the geometry shape.
    sweep = terrestrial_budget(_terr(0.02, np.array([2e3, 5e3, 10e3]),
                                     near_aperture=0.3),
                               HorizontalPath(np.array([2e3, 5e3, 10e3])))
    assert np.shape(sweep.total_loss_db()) == (3,)

    # The scintillation slot is reserved: it raises, and so does scintillation=True.
    try:
        terrestrial_scintillation_term(scn, geom)
    except NotImplementedError:
        pass
    else:
        raise AssertionError("scintillation slot must raise NotImplementedError")
    try:
        terrestrial_budget(scn, geom, scintillation=True)
    except NotImplementedError:
        pass
    else:
        raise AssertionError("scintillation=True must raise NotImplementedError")

    print(budget.to_frame().to_string(index=False))
    print(f"\ntotal (5 km, clean): {budget.total_loss_db():.2f} dB")
    print(f"99% analytic fade:   {fade:.2f} dB")
    print(f"with launch truncation: {budget_ap.total_loss_db():.2f} dB")
    print(f"SMF coupling (mean-only) loss: none={float(loss_none):.2f} dB  "
          f"tip-tilt={float(loss_tt):.2f} dB  AO200={float(loss_ao):.2f} dB "
          f"(fade margin refused: fidelity-0)")
    print("self-check passed")
