'''
Terrestrial (horizontal-path) Terms and budget assembly.

This module builds the terrestrial link budget: a ground-to-ground horizontal
path. It reuses the direction-agnostic Terms (geometric spreading, pointing
jitter, transmit truncation), adds the horizontal Beer-Lambert extinction Term,
and adds the horizontal Gaussian-beam scintillation Term.

A horizontal path differs from a space link in two ways that matter here. The
range is a constant path length, not a slant range that changes with elevation.
And the Gaussian-beam properties (waist, divergence, curvature) steer the
turbulence result strongly, because the whole path sits in the near field of a
finite beam. So the horizontal scintillation is NOT the plane-wave slant-path
model that the downlink uses. It uses the Gaussian-beam analytic form.

The scintillation Term reuses two pieces of physics that the package already
carries. The point scintillation index is the on-axis Gaussian beam-wave index
sigma2_I(0, L) of olb.turbulence.beam_scintillation (Dios et al., Applied Optics
43 (2004) 3866, Eq. 16), evaluated on a constant-Cn2 horizontal grid. The
aperture-averaging win is the Andrews weak-turbulence Kolmogorov factor of
olb.turbulence.scintillation.aperture_averaging_factor_weak (Andrews and
Phillips, 2nd ed. 2005, Ch. 10). The lognormal fade faces (mean_db, quantile,
sampler) follow the same closed form as the downlink lognormal Term (see
olb.links.downlink._lognormal_term).
'''

import warnings

import numpy as np
from scipy.stats import norm

from ..results import Budget, Term
from ..assumptions import (Assumptions, BEAM_GAUSSIAN, REGIME_WEAK,
                           SPECTRUM_KOLMOGOROV)
from ..models.geometric import geometric_loss_term
from ..models.transmittance import terrestrial_extinction_term
from ..models.pointing import pointing_loss_term
from ..models.gaussian_efficiency import tx_gaussian_efficiency_term
from ..turbulence.beam_scintillation import on_axis_scintillation_index
from ..turbulence.scintillation import (aperture_averaging_factor_weak,
                                        WEAK_FLUCTUATION_LIMIT)

# Below this launch-truncation loss the beam is an untruncated Gaussian, so the
# transmit Gaussian-efficiency term is skipped [dB]. Matches olb.links.uplink.
TX_TRUNCATION_MIN_DB = 1e-2

# Natural log of ten, for the dB conversions in the lognormal faces.
_LN10 = np.log(10.0)

# Points on the constant-Cn2 horizontal grid for the scintillation integral. The
# on-axis index integrates over the path, so a few hundred points converge it.
_SCINT_GRID_N = 400


def terrestrial_scintillation_term(scenario, geometry, *, n_grid=_SCINT_GRID_N):
    '''
    Build the horizontal Gaussian-beam scintillation Term for an aperture receiver.

    This Term gives the analytic lognormal turbulence fade of a power-in-bucket
    (Aperture) receiver on a horizontal path. It carries a real fade, so it has
    all three faces (mean_db, quantile, sampler).

    Physics (two pieces the package already carries, both cited in place):

      1. The on-axis Gaussian beam-wave scintillation index sigma2_I(0, L). This
         is the finite-beam point index. It sits between the plane-wave limit
         (large waist) and the spherical-wave limit (small waist), so the beam
         waist and the range steer it. See
         olb.turbulence.beam_scintillation.on_axis_scintillation_index (Dios et
         al., Applied Optics 43 (2004) 3866, Eq. 16). A horizontal path is a
         constant-Cn2 grid from the transmitter (z=0) to the receiver (z=L) at
         elevation 90 (sec=1), so the range L is the grid length.

      2. The aperture-averaging factor A that reduces the point index to the
         aperture-averaged flux index sigma2_P = A * sigma2_I. A telescope of
         diameter D averages the fine scintillation, so A falls as D grows. See
         olb.turbulence.scintillation.aperture_averaging_factor_weak (Andrews
         and Phillips, 2nd ed. 2005, Ch. 10, weak Kolmogorov, small inner scale).

    The flux index sigma2_P then feeds the same lognormal fade faces as the
    downlink lognormal Term (olb.links.downlink._lognormal_term). With
    sigma_l^2 = ln(1 + sigma2_P) the faces are
        mean_db     = (5/ln10) * sigma_l^2
        quantile(p) = -10*log10( exp(-sigma_l^2/2 + sigma_l * Phi_inv(1-p)) )
    Source of the lognormal irradiance PDF: Andrews and Phillips, 2nd ed. (2005),
    Ch. 5. Phi_inv is the inverse standard normal CDF.

    Validity: the lognormal model is a weak-fluctuation model. The Term carries
    sigma2_I in Term.meta and sets a weak_fluctuation_valid flag. It flags a
    violation and gives a warning when sigma2_I exceeds WEAK_FLUCTUATION_LIMIT.
    Above that limit use a gamma-gamma or a Monte Carlo model.

    Parameters:
        scenario : TerrestrialScenario
            tx = near (its Transmitter waist launches the beam); rx = far (its
            aperture diameter D and wavelength set the averaging). Reads the
            channel path length L and the constant Cn2.
        geometry : HorizontalPath
            Unused here (the path length and Cn2 come from the channel). Kept for
            the f(scenario, geometry) -> Term signature.
        n_grid : int
            Points on the constant-Cn2 path grid for the index integral.

    Returns:
        Term
            name="scintillation", category="turbulence". It has a real fade.

    Raises:
        ValueError
            If the near terminal has no Transmitter (no launch beam).
    '''
    tx = scenario.tx_terminal
    rx = scenario.rx_terminal
    if tx.transmitter is None:
        raise ValueError(
            "terrestrial scintillation needs a launch beam: set the near terminal "
            "transmitter = Transmitter(waist_m=...)."
        )
    w0 = tx.transmitter.waist_m
    D = rx.aperture_m
    wavelength = rx.wavelength_m
    L = float(scenario.channel.path_length_m)
    cn2 = float(scenario.channel.cn2)

    # Horizontal path: distance-from-transmitter grid, constant Cn2, elevation 90
    # (sec=1), L = the grid length. This gives the on-axis index sigma2_I(0, L).
    hs = np.linspace(0.0, L, int(n_grid))
    cn2_profile = np.full_like(hs, cn2)
    sigma2_I = float(on_axis_scintillation_index(
        hs, cn2_profile, w0, wavelength, elevation_deg=90.0, path_length_m=None))

    # TODO(pointing jitter): this Term is ON-AXIS only (r=0), so it does NOT yet
    # carry pointing/tracking jitter. When adding it, note the ANALYTIC-PATH
    # ASYMMETRY vs the uplink MC (olb.turbulence.coupled_flux, which folds jitter
    # into r=beta and gets everything for free because it works in ABSOLUTE
    # irradiance):
    #   1. Fold the jitter displacement into the OFF-AXIS radial index at r=beta
    #      (beam_scintillation.radial_scintillation_index), beta drawn from the
    #      jitter (+ beam-wander) 2-D variance. This adds the FLUCTUATION boost.
    #   2. You STILL need a SEPARATE mean-power loss term for the same jitter,
    #      because the off-axis sigma2_I is normalised to the LOCAL mean and does
    #      NOT carry the exp(-2*beta^2/W^2) mean drop. The MC path merges the two;
    #      the analytic path cannot. Adding both is NOT double-counting -- they
    #      are different statistical moments (variance vs mean).
    # This is the reason to converge the analytic and MC Dios paths rather than
    # patch one. See memory dios-scintillation-convergence / pointing-jitter-into-beta.

    # Aperture-averaging win: a larger D averages more, so A and sigma2_P fall.
    A = float(aperture_averaging_factor_weak(D, wavelength, L))
    sigma2_P = A * sigma2_I

    sigma_l2 = np.log(1.0 + sigma2_P)
    sigma_l = np.sqrt(sigma_l2)

    mean_db = (5.0 / _LN10) * sigma_l2

    def quantile(p):
        # Loss exceeded a fraction (1 - p) of the time. Phi_inv = norm.ppf.
        z = norm.ppf(1.0 - p)
        return -10.0 / _LN10 * (-sigma_l2 / 2.0 + sigma_l * z)

    def sampler(n, rng):
        P = rng.lognormal(mean=-sigma_l2 / 2.0, sigma=sigma_l, size=n)
        return -10.0 * np.log10(P)

    # The Rytov (lognormal) validity is a turbulence property. Test it with the
    # point index sigma2_I, not the aperture-averaged sigma2_P.
    valid = sigma2_I < WEAK_FLUCTUATION_LIMIT
    assumptions = Assumptions(
        beam_type=BEAM_GAUSSIAN,
        turbulence_regime=REGIME_WEAK,
        spectrum=SPECTRUM_KOLMOGOROV,
        validity="Weak fluctuation: sigma2_I < 0.25. The point index is the "
                 "on-axis Gaussian beam-wave index over the constant-Cn2 "
                 "horizontal path (Dios et al. 2004, Eq. 16). The weak "
                 "aperture-averaging factor assumes a uniform circular aperture "
                 "with no central obscuration and a small inner scale (Andrews "
                 "and Phillips 2005, Ch. 10).",
    )
    # The weak Kolmogorov averaging factor models a uniform circular aperture. A
    # central obscuration (Cassegrain secondary) breaks that. Flag the violation.
    if rx.obscuration_ratio > 0.0:
        assumptions.flag(
            f"The far aperture has a central obscuration "
            f"(ratio={rx.obscuration_ratio:.3f}); the weak aperture-averaging "
            "factor assumes a uniform circular aperture and does not model it."
        )
    if not valid:
        assumptions.flag(
            f"sigma2_I={sigma2_I:.3f} exceeds the weak-fluctuation limit "
            f"{WEAK_FLUCTUATION_LIMIT}; the lognormal fade is not trusted. Use "
            "gamma-gamma or Monte Carlo."
        )
        warnings.warn(
            f"horizontal Gaussian-beam scintillation index sigma2_I={sigma2_I:.3f} "
            f">= {WEAK_FLUCTUATION_LIMIT} -- the Rytov weak-fluctuation model is "
            "exceeded. The lognormal fade is not trusted. Use a gamma-gamma or a "
            "Monte Carlo model."
        )

    return Term(
        name="scintillation",
        category="turbulence",
        mean_db=float(mean_db),
        sampler=sampler,
        quantile=quantile,
        note="horizontal Gaussian-beam lognormal scintillation, aperture-averaged",
        meta={
            "model": "lognormal",
            "beam": "gaussian-horizontal",
            "sigma2_I": sigma2_I,
            "sigma2_P": float(sigma2_P),
            "aperture_averaging_factor": A,
            "weak_fluctuation_valid": bool(valid),
            "weak_fluctuation_limit": WEAK_FLUCTUATION_LIMIT,
        },
        assumptions=assumptions,
    )


def terrestrial_budget(scenario, geometry, *, scintillation=True, turbulence=True):
    '''
    Assemble the terrestrial budget: geometric, extinction, pointing, turbulence.

    The deterministic Terms (geometric spreading, horizontal extinction, pointing
    jitter) are exact and direction-agnostic. The turbulence effect depends on the
    receiver front end, and it mirrors the downlink rule:

      - An Aperture (bucket) detector, or no detector, gets the horizontal
        Gaussian-beam scintillation Term (terrestrial_scintillation_term). It is
        a real analytic fade. This is the natural default for a realistic
        aperture budget, so scintillation defaults to True.
      - An SMF detector gets the fidelity-0 mean-only fibre-coupling Term instead
        (olb.models.coupling.terrestrial_smf_coupling_term). The coupling loss IS
        the turbulence effect for the fibre, so the scintillation Term is NOT
        also added (no double-count). The mean-only coupling Term locks the budget
        to fidelity 0, so the budget then refuses a fade margin. When the SMF sets
        the coupling optics (focal_length_m and mode_field_radius_m), the budget
        also adds the receive tip-tilt walk-off fade Term (smf_walkoff_term).
      - An MMF detector gets the multimode-fibre coupling Term instead
        (olb.models.coupling.mmf_coupling_term): the geometric spot-in-core loss
        plus the tip-tilt walk-off fade. It replaces the scintillation Term (no
        double-count). The MMF Term has a real fade, so the budget keeps its fade.

    Set scintillation=False to drop the scintillation Term and keep only the
    deterministic Terms (for example to sweep an array path length, where the
    scalar-only scintillation Term does not broadcast; loop per distance instead).

    Parameters:
        scenario : TerrestrialScenario
            A terrestrial link case. tx = near end, rx = far end. Its
            TerrestrialChannel carries path_length_m, attenuation_db_per_km, cn2.
        geometry : HorizontalPath
            The horizontal path (reads slant_range_m = path length).
        scintillation : bool
            Add the horizontal Gaussian-beam scintillation Term for an aperture /
            no-detector receiver when True (the default). An SMF detector always
            replaces it with the coupling Term.
        turbulence : bool
            Master turbulence switch. When False, drop EVERY turbulence quantity:
            no scintillation Term, and the fibre-coupling Terms keep only their
            static parts. An SMF coupling Term becomes the static mode-match loss,
            an MMF Term keeps its spot-overfill loss, and the walk-off Term keeps
            only the receive mechanical jitter (the beam-wander tilt drops). The
            deterministic Terms (geometric, extinction, launch truncation) and the
            transmit pointing jitter stay. So a coupling budget with angular jitter
            still runs, only without turbulence.

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
    # The receive-side turbulence effect. An SMF detector on the far terminal
    # takes the fidelity-0 (mean-only) fibre-coupling Term, using the horizontal
    # Gaussian-beam r0 and the compensation stack. That coupling loss REPLACES the
    # scintillation Term (no double-count): it IS the turbulence effect for the
    # fibre. It is MEAN-ONLY, so it locks the budget to fidelity 0, and the budget
    # then refuses a fade margin. An Aperture (bucket) detector, or no detector,
    # is phase-insensitive; its turbulence penalty is the scintillation Term.
    from ..terminal import SMF, MMF
    rx = scenario.rx_terminal
    if isinstance(rx.detector, SMF):
        # Lazy import breaks the terrestrial <-> coupling import cycle.
        from ..models.coupling import (terrestrial_smf_coupling_term,
                                       smf_walkoff_term)
        # The receive tip-tilt walk-off fade fires when the fibre-coupling optics
        # are set (focal length + mode field radius). Without them a tip-tilt has
        # no focal-plane displacement, so the Term is skipped. When the walk-off
        # fires, it carries the tip-tilt coupling loss. So the coupling Term keeps
        # the HIGHER-ORDER residual only (drop_tiptilt=True). This stops the
        # tip-tilt from being counted two times.
        walkoff_on = (getattr(rx.detector, "optimal_focus", False)
                      or (rx.detector.focal_length_m is not None
                          and rx.detector.mode_field_radius_m is not None))
        terms.append(terrestrial_smf_coupling_term(scenario, geometry,
                                                   drop_tiptilt=walkoff_on,
                                                   turbulence=turbulence))
        if walkoff_on:
            terms.append(smf_walkoff_term(scenario, geometry,
                                          turbulence=turbulence))
    elif isinstance(rx.detector, MMF):
        # An MMF (light bucket) replaces the scintillation Term with the geometric
        # spot-in-core coupling plus the tip-tilt walk-off fade (no double-count).
        from ..models.coupling import mmf_coupling_term
        terms.append(mmf_coupling_term(scenario, geometry, turbulence=turbulence))
    elif scintillation and turbulence:
        terms.append(terrestrial_scintillation_term(scenario, geometry))
    return Budget(terms, scenario=scenario)


if __name__ == '__main__':
    from ..scenario import TerrestrialScenario, TerrestrialChannel
    from ..geometry import HorizontalPath
    from ..terminal import Terminal, Transmitter, Aperture, SMF, MMF, TipTilt, AO

    def _terr(w0, L, *, divergence=None, power=None, jitter=0.0,
              near_aperture=0.1, near_obscuration=0.0, far_aperture=0.1,
              far_obscuration=0.0, attenuation=0.5, sensitivity=None,
              cn2=3e-16):
        '''Build a TerrestrialScenario: tx = near, rx = far. Weak Cn2 by default.'''
        detector = None if sensitivity is None else Aperture(sensitivity_dbm=sensitivity)
        return TerrestrialScenario(
            near=Terminal(aperture_m=near_aperture, obscuration_ratio=near_obscuration,
                          wavelength_m=1550e-9, pointing_jitter_rad=jitter,
                          transmitter=Transmitter(waist_m=w0, power_dbm=power,
                                                  divergence_rad=divergence)),
            far=Terminal(aperture_m=far_aperture, obscuration_ratio=far_obscuration,
                         wavelength_m=1550e-9, detector=detector),
            channel=TerrestrialChannel(path_length_m=L, attenuation_db_per_km=attenuation,
                                       cn2=cn2))

    # A clean 5 km link. A wide near aperture (0.3 m for a 0.02 m waist) leaves
    # the beam untruncated, so the launch-truncation Term does not fire. The Cn2
    # default (3e-16) keeps the path in the weak regime.
    scn = _terr(0.02, 5e3, power=30, jitter=5e-6, sensitivity=-40,
                near_aperture=0.3)
    geom = HorizontalPath(5e3)
    budget = terrestrial_budget(scn, geom)
    names = [t.name for t in budget.terms]
    # The scintillation Term is now in the aperture budget (the default).
    assert names == ["geometric spreading", "atmospheric extinction (horizontal)",
                     "pointing jitter", "scintillation"], names
    # The extinction Term is exact: 5 km * 0.5 dB/km = 2.5 dB.
    ext = next(t for t in budget.terms if t.category == "atmospheric")
    assert np.isclose(ext.mean_db, 2.5), ext.mean_db

    # --- horizontal Gaussian-beam scintillation Term ------------------------
    scint = next(t for t in budget.terms if t.name == "scintillation")
    # sigma2_I > 0, and the aperture averages it down: A < 1, sigma2_P < sigma2_I.
    assert scint.meta["sigma2_I"] > 0.0, scint.meta["sigma2_I"]
    assert 0.0 < scint.meta["aperture_averaging_factor"] < 1.0
    assert scint.meta["sigma2_P"] < scint.meta["sigma2_I"]
    # This scenario stays weak, so the Term is valid (no flag, no warning).
    assert scint.meta["weak_fluctuation_valid"] and scint.assumptions.ok
    # It has a real fade: a working analytic 99% quantile deeper than the mean.
    q99_scint = scint.quantile_db(0.99)
    assert np.isfinite(q99_scint) and q99_scint > scint.mean_db, (q99_scint, scint.mean_db)

    # The aperture-averaging win: a larger receive aperture shrinks the flux index
    # and the fade. Sweep D and check both fall monotonically.
    D_sweep = [0.05, 0.1, 0.2, 0.4, 0.8]
    sig_P, fades = [], []
    for D in D_sweep:
        s = terrestrial_scintillation_term(
            _terr(0.02, 5e3, far_aperture=D), geom)
        sig_P.append(s.meta["sigma2_P"])
        fades.append(float(s.quantile_db(0.99)))
    assert all(np.diff(sig_P) < 0.0), sig_P        # flux index shrinks with D
    assert all(np.diff(fades) < 0.0), fades        # 99% fade shrinks with D

    # The weak-fluctuation flag trips at a strong Cn2 (and again at a long path).
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        strong = terrestrial_scintillation_term(
            _terr(0.02, 5e3, far_aperture=0.2, cn2=1e-13), geom)
        long_path = terrestrial_scintillation_term(
            _terr(0.02, 20e3, far_aperture=0.2, cn2=1e-14), HorizontalPath(20e3))
    assert not strong.meta["weak_fluctuation_valid"] and not strong.assumptions.ok
    assert not long_path.meta["weak_fluctuation_valid"] and not long_path.assumptions.ok
    # A central obscuration on the far aperture flags the circular-aperture filter.
    obsc = terrestrial_scintillation_term(
        _terr(0.02, 5e3, far_aperture=0.2, far_obscuration=0.3), geom)
    assert any("central obscuration" in v for v in obsc.assumptions.violations)

    # The aperture budget now carries the turbulence fade, so its 99% fade is
    # DEEPER than the pointing-only (scintillation-off) budget.
    budget_noscint = terrestrial_budget(scn, geom, scintillation=False)
    fade_with = budget.fade_margin_db(0.99)
    fade_without = budget_noscint.fade_margin_db(0.99)
    assert np.isfinite(fade_with) and np.isfinite(fade_without)
    assert fade_with > fade_without, (fade_with, fade_without)

    # A narrow near aperture (0.02 m for a 0.02 m waist) truncates the beam, so
    # the launch-truncation Term fires.
    scn_ap = _terr(0.02, 5e3, power=30, near_aperture=0.02, near_obscuration=0.2)
    budget_ap = terrestrial_budget(scn_ap, geom)
    assert "transmit Gaussian efficiency" in [t.name for t in budget_ap.terms]
    assert budget_ap.total_loss_db() > budget.total_loss_db()

    # A longer path costs more geometric spread AND more extinction. Keep the same
    # aperture and a weak Cn2 so the comparison is deterministic-plus-turbulence.
    long_budget = terrestrial_budget(_terr(0.02, 10e3, near_aperture=0.3,
                                           sensitivity=-40),
                                     HorizontalPath(10e3))
    assert long_budget.total_loss_db() > budget.total_loss_db()

    # The aperture budget has an analytic fade (every Term has a quantile).
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
    # An SMF detector replaces the scintillation Term (no double-count).
    assert "scintillation" not in [t.name for t in smf_budget.terms]
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

    # An SMF with the coupling optics set also adds the receive tip-tilt walk-off
    # Term. The budget stays fidelity-0 (mean-only coupling term still present).
    with _warnings.catch_warnings():
        _warnings.simplefilter("ignore")
        scn_opt = TerrestrialScenario(
            near=Terminal(aperture_m=0.3, wavelength_m=1550e-9,
                          transmitter=Transmitter(waist_m=0.02, power_dbm=30)),
            far=Terminal(aperture_m=0.2, wavelength_m=1550e-9,
                         pointing_jitter_rad=5e-6,
                         detector=SMF(focal_length_m=0.02,
                                      mode_field_radius_m=5.2e-6,
                                      sensitivity_dbm=-40)),
            channel=TerrestrialChannel(path_length_m=3e3, attenuation_db_per_km=0.5,
                                       cn2=1e-14))
        smf_opt_budget = terrestrial_budget(scn_opt, HorizontalPath(3e3))
    assert "SMF tip-tilt walk-off" in [t.name for t in smf_opt_budget.terms]
    assert not smf_opt_budget.provides_fade   # mean-only coupling term still locks it

    # --- MMF (multimode-fibre light bucket) ---------------------------------
    def _mmf(core_radius=25e-6, focal=None, jitter=5e-6, far_aperture=0.2,
             cn2=1e-15, optimal_focus=True):
        # optimal_focus fills the spot to the core, so a tip-tilt walks it off and
        # the Term carries a real fade. A weak Cn2 keeps the offset moderate.
        scn = TerrestrialScenario(
            near=Terminal(aperture_m=0.3, wavelength_m=1550e-9,
                          transmitter=Transmitter(waist_m=0.02, power_dbm=30)),
            far=Terminal(aperture_m=far_aperture, wavelength_m=1550e-9,
                         pointing_jitter_rad=jitter,
                         detector=MMF(core_radius_m=core_radius, focal_length_m=focal,
                                      optimal_focus=optimal_focus, sensitivity_dbm=-38)),
            channel=TerrestrialChannel(path_length_m=3e3, attenuation_db_per_km=0.5,
                                       cn2=cn2))
        return terrestrial_budget(scn, HorizontalPath(3e3))

    mmf_budget = _mmf()
    mmf_names = [t.name for t in mmf_budget.terms]
    # The MMF coupling Term replaces the scintillation Term (no double-count).
    assert "receive coupling (MMF)" in mmf_names and "scintillation" not in mmf_names
    mmf_term = next(t for t in mmf_budget.terms if t.category == "coupling")
    # The MMF Term has a real fade, so the budget keeps its fade margin.
    assert not mmf_term.mean_only and mmf_budget.provides_fade
    assert np.isfinite(mmf_budget.fade_margin_db(0.99))
    assert mmf_budget.fade_margin_db(0.99) > mmf_budget.total_loss_db()

    # A bucket (Aperture) detector adds the scintillation Term, and the budget
    # keeps its analytic fade (scintillation has a quantile; no mean-only term).
    assert budget.provides_fade
    assert np.isfinite(budget.fade_margin_db(0.99))

    # --- master turbulence switch (turbulence=False) ------------------------
    # An aperture budget with turbulence off drops the scintillation Term but keeps
    # the deterministic Terms and the transmit pointing jitter (still a real fade).
    off = terrestrial_budget(scn, geom, turbulence=False)
    off_names = [t.name for t in off.terms]
    assert "scintillation" not in off_names, off_names
    assert off.provides_fade and np.isfinite(off.fade_margin_db(0.99))
    # An SMF with the walk-off optics: turbulence off keeps the static mode-match
    # coupling (deterministic, NOT mean-only) plus the jitter walk-off fade. The
    # walk-off carries the receive jitter alone (no beam-wander tilt), so the
    # budget still reports a fade and the jitter still drives the coupling.
    smf_wo = SMF(focal_length_m=0.02, mode_field_radius_m=5.2e-6, sensitivity_dbm=-40)
    scn_off = TerrestrialScenario(
        near=Terminal(aperture_m=0.3, wavelength_m=1550e-9,
                      transmitter=Transmitter(waist_m=0.02, power_dbm=30)),
        far=Terminal(aperture_m=0.2, wavelength_m=1550e-9, pointing_jitter_rad=8e-6,
                     detector=smf_wo),
        channel=TerrestrialChannel(path_length_m=3e3, attenuation_db_per_km=0.5, cn2=1e-14))
    smf_off = terrestrial_budget(scn_off, HorizontalPath(3e3), turbulence=False)
    cpl_off = next(t for t in smf_off.terms if t.category == "coupling")
    wo_off = next(t for t in smf_off.terms if t.name == "SMF tip-tilt walk-off")
    assert cpl_off.meta["model"] == "static" and not cpl_off.mean_only
    assert wo_off.meta["sigma2_wander"] == 0.0 and wo_off.meta["sigma2_jitter"] > 0.0
    assert smf_off.provides_fade and wo_off.quantile_db(0.99) > wo_off.mean_db
    # The jitter drives the coupling fade even with turbulence off (the whole point).
    scn_calm = TerrestrialScenario(
        near=scn_off.near, far=Terminal(aperture_m=0.2, wavelength_m=1550e-9,
                                        pointing_jitter_rad=1e-6, detector=smf_wo),
        channel=scn_off.channel)
    wo_calm = next(t for t in terrestrial_budget(scn_calm, HorizontalPath(3e3),
                                                 turbulence=False).terms
                   if t.name == "SMF tip-tilt walk-off")
    assert wo_off.mean_db > wo_calm.mean_db, (wo_off.mean_db, wo_calm.mean_db)

    # A path-length sweep is scalar-only for the scintillation and coupling Terms,
    # so an array geometry must turn scintillation off (or loop per distance).
    sweep = terrestrial_budget(_terr(0.02, np.array([2e3, 5e3, 10e3]),
                                     near_aperture=0.3),
                               HorizontalPath(np.array([2e3, 5e3, 10e3])),
                               scintillation=False)
    assert np.shape(sweep.total_loss_db()) == (3,)

    print(budget.to_frame().to_string(index=False))
    print(f"\ntotal (5 km, clean): {budget.total_loss_db():.2f} dB")
    print(f"scintillation: sigma2_I={scint.meta['sigma2_I']:.4f} "
          f"A={scint.meta['aperture_averaging_factor']:.4f} "
          f"sigma2_P={scint.meta['sigma2_P']:.4f}")
    print(f"99% fade: turbulence-on={fade_with:.2f} dB  "
          f"pointing-only={fade_without:.2f} dB")
    print("aperture-averaging win (99% scintillation fade vs D):")
    for D, f in zip(D_sweep, fades):
        print(f"    D={D * 100:5.1f} cm -> {f:.3f} dB")
    print(f"with launch truncation: {budget_ap.total_loss_db():.2f} dB")
    print(f"SMF coupling (mean-only) loss: none={float(loss_none):.2f} dB  "
          f"tip-tilt={float(loss_tt):.2f} dB  AO200={float(loss_ao):.2f} dB "
          f"(fade margin refused: fidelity-0)")
    print("self-check passed")
