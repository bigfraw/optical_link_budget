'''
Downlink scintillation Term and budget assembly.

This module builds the downlink scintillation Term and the downlink budget. The
Term gives the analytic lognormal plane-wave fade for a LEO-to-ground link, or
the gamma-gamma fade when the turbulence is too strong for the lognormal model.
The pure scintillation physics lives in olb.turbulence.plane_wave_scintillation
and in the Andrews foundation layer, olb.turbulence.andrews.

The received power P is lognormal with E[P] = 1. With sigma_l^2 = ln(1 +
sigma2_P) the loss faces are
    mean_db     = (5/ln10) * sigma_l^2
    quantile(p) = -10*log10( exp(-sigma_l^2/2 + sigma_l * Phi_inv(1-p)) )
Source of the lognormal irradiance PDF: Andrews and Phillips, Laser Beam
Propagation through Random Media, 2nd ed. (2005), Ch. 5. Phi_inv is the inverse
standard normal CDF.

Validity: the lognormal Term is trusted while the point index sigma2_I stays
below the lognormal-PDF house rule LOGNORMAL_PDF_LIMIT (0.25); the gamma-gamma
Term carries a separate rytov_regime label from the REGIME gate (boundary
sigma_R^2 = 1). The lognormal Term warns when sigma2_I exceeds LOGNORMAL_PDF_LIMIT.
Above that limit use model="auto", which selects the gamma-gamma Term. The gamma-gamma Term holds for
every fluctuation strength: Andrews and Phillips, 2nd ed. (2005),
DOI 10.1117/3.626196, Ch. 9, Eqs. (137) and (138), printed p. 370, and Ch. 12,
Eq. (40), printed p. 497.
'''

import warnings

import numpy as np

from ..results import Budget
from ..assumptions import (trace_assumptions, BEAM_PLANE_WAVE, REGIME_STRONG,
                          REGIME_WEAK, SPECTRUM_KOLMOGOROV)
from ..models.fade import irradiance_fade_term
from ..models.geometric import geometric_loss_term
from ..models.gaussian_efficiency import tx_gaussian_efficiency_term
from .uplink import TX_TRUNCATION_MIN_DB
from ..models.extinction import slant_extinction_term, DEFAULT_TAU_ZENITH
from ..models.pointing import pointing_loss_term
from ..turbulence.andrews.distributions import (gamma_gamma_mean_log,
                                                gamma_gamma_params,
                                                gamma_gamma_quantile,
                                                gamma_gamma_rvs,
                                                gamma_gamma_scintillation_index,
                                                lognormal_params,
                                                lognormal_mean_log,
                                                lognormal_quantile,
                                                lognormal_rvs)
from ..turbulence.andrews.scintillation import (large_scale_log_variance,
                                                small_scale_log_variance,
                                                rytov_weak, LOGNORMAL_PDF_LIMIT,
                                                WEAK_REGIME_LIMIT)
from ..turbulence.profiles import DEFAULT_HS, default_cn2_profile
from ..turbulence.plane_wave_scintillation import (plane_wave_scintillation_index,
                                        aperture_averaged_scintillation_index)


def _lognormal_term(scenario, geometry, *, aperture_average, hs, cn2_profile):
    '''
    Build the analytic lognormal downlink scintillation Term.

    Give all three faces: mean_db, quantile(p), sampler(n, rng). Carry sigma2_I
    and the weak_fluctuation_valid flag in Term.meta. Warn when sigma2_I exceeds
    the weak-fluctuation limit.
    '''
    rx = scenario.rx_terminal
    wavelength = rx.wavelength_m
    D = rx.aperture_m
    elev = geometry.elevation_deg

    # Open the collection context around the PHYSICS CALLS only. Each decorated
    # function registers its own assumptions, so the Term inherits the union.
    with trace_assumptions() as trace:
        sigma2_I = plane_wave_scintillation_index(elev, wavelength, hs,
                                                  cn2_profile)
        if aperture_average:
            # Use the distributed-path aperture-averaging integral directly.
            sigma2_P = aperture_averaged_scintillation_index(D, elev, wavelength,
                                                             hs, cn2_profile)
    if aperture_average:
        A = np.asarray(sigma2_P) / np.asarray(sigma2_I)
    else:
        A = 1.0
        sigma2_P = sigma2_I

    # The lognormal irradiance model. lognormal_params turns the aperture-averaged
    # index into the log-irradiance variance; the three dB faces come from the ONE
    # shared adapter (olb.models.fade.irradiance_fade_term), the SAME path the
    # gamma-gamma Term uses (backlog I-2 / crosscheck TL-01..04).
    sigma_l2 = lognormal_params(sigma2_P)
    base_shape = np.shape(sigma_l2)

    # This is a LOGNORMAL Term, so its binding validity is the lognormal-PDF
    # house rule on the point index sigma2_I (NOT the aperture-averaged sigma2_P,
    # and NOT the regime boundary sigma_R^2 = 1, which is 4x looser -- see
    # Conflict C-05). Test the PDF shape with LOGNORMAL_PDF_LIMIT.
    valid = np.asarray(sigma2_I) < LOGNORMAL_PDF_LIMIT
    # The traced physics functions own the beam type, the regime, the spectrum,
    # the slant-geometry constraint, and the circular-aperture / no-obscuration
    # constraints; the merge inherits their union and any traced violation. State
    # the three headline fields explicitly (this is a plane-wave lognormal Term).
    assumptions = trace.merge(
        beam_type=BEAM_PLANE_WAVE,
        turbulence_regime=REGIME_WEAK,
        spectrum=SPECTRUM_KOLMOGOROV,
        validity="Lognormal fade PDF trusted for sigma2_I < LOGNORMAL_PDF_LIMIT "
                 "(0.25, the Ch. 11.3 optimistic-tail house rule, NOT the wider "
                 "regime boundary sigma_R^2 = 1). Aperture averaging uses the "
                 "distributed-path spectral integral over height and spatial "
                 "wavenumber. The aperture-averaging filter assumes a uniform "
                 "circular aperture with no central obscuration.",
    )
    # The central obscuration is a scenario-level fact the physics never sees (the
    # circular-aperture filter takes no obscuration ratio). A Cassegrain secondary
    # breaks the unobscured-aperture constraint, so flag it at the factory level.
    obscuration = rx.obscuration_ratio
    if aperture_average and obscuration > 0.0:
        assumptions.flag(
            f"The receive aperture has a central obscuration "
            f"(ratio={obscuration:.3f}); the circular-aperture averaging "
            "filter does not model it.",
            source="factory:links.downlink",
        )
    # The lognormal-PDF house rule (0.25 on sigma2_I) is a PDF-SHAPE decision that
    # the physics does not gate (the traced regime check is the WIDER sigma_R^2 = 1
    # boundary). So it stays a factory flag, source-tagged the same way.
    if not np.all(valid):
        worst = float(np.asarray(sigma2_I)[~valid].max())
        assumptions.flag(
            f"sigma2_I={worst:.3f} exceeds the lognormal-PDF house limit "
            f"{LOGNORMAL_PDF_LIMIT}; the lognormal fade tail is optimistic. Use "
            "gamma-gamma (model='auto') or Monte Carlo.",
            source="factory:links.downlink",
        )
        warnings.warn(
            f"plane-wave scintillation index sigma2_I={np.max(sigma2_I):.3f} >= "
            f"{LOGNORMAL_PDF_LIMIT} -- the lognormal fade tail is not trusted. "
            "Use the gamma-gamma model (model='auto') or the Monte Carlo model."
        )

    return irradiance_fade_term(
        "scintillation", "turbulence",
        mean_log=lognormal_mean_log(sigma_l2),
        quantile=lambda p: lognormal_quantile(p, sigma_l2),
        rvs=lambda n, rng: lognormal_rvs(n, sigma_l2, rng),
        note="plane-wave downlink lognormal scintillation, aperture-averaged",
        meta={
            "model": "lognormal",
            "sigma2_I": float(sigma2_I) if base_shape == () else np.asarray(sigma2_I),
            "sigma2_P": float(sigma2_P) if base_shape == () else np.asarray(sigma2_P),
            "aperture_averaging_factor": float(A) if np.ndim(A) == 0 else np.asarray(A),
            "weak_fluctuation_valid": bool(valid) if base_shape == () else valid,
            "weak_fluctuation_limit": LOGNORMAL_PDF_LIMIT,
        },
        assumptions=assumptions,
    )


def _gamma_gamma_term(scenario, geometry, *, aperture_average, hs, cn2_profile):
    '''
    Build the analytic gamma-gamma downlink scintillation Term.

    The gamma-gamma model covers the moderate-to-strong fluctuation regime. Use
    it when sigma2_I meets or passes the weak-fluctuation limit. The Term gives
    all three faces: mean_db, quantile(p), sampler(n, rng).

    The Term COMPOSES the Andrews foundation layer. It derives no new physics:
        sigma_R^2       the same slant plane-wave index the lognormal Term uses,
                        Ch. 12, Eq. (38), printed p. 495
        sigma_lnX^2     andrews.scintillation.large_scale_log_variance,
                        Ch. 9, Eq. (41), printed p. 335
        sigma_lnY^2     andrews.scintillation.small_scale_log_variance,
                        Ch. 9, Eq. (46), printed p. 336
        alpha, beta     andrews.distributions.gamma_gamma_params,
                        Ch. 9, Eq. (138), printed p. 370
        the dB faces    olb.models.fade.irradiance_fade_term
    Source: Andrews and Phillips, Laser Beam Propagation through Random Media,
    2nd ed. (2005), DOI 10.1117/3.626196.

    The scintillation index of the result, 1/alpha + 1/beta + 1/(alpha beta)
    (Ch. 9, Eq. (139), printed p. 371), is identically
    exp(sigma_lnX^2 + sigma_lnY^2) - 1, which is the book's own weak-to-strong
    downlink index, Ch. 12, Eq. (40), printed p. 497. So this Term and
    andrews.paths.downlink_scintillation_index(regime="strong") agree exactly.
    The module self-check measures that identity.

    POINT RECEIVER. The book gives NO aperture-averaged downlink index in the
    moderate-to-strong regime: Ch. 12, Eq. (39), printed p. 496, is a weak-theory
    form and Ch. 12, Eq. (40) is a point form, and the book prints no product of
    the two. So this Term models a POINT receiver. Its fade is deeper than the
    fade of a real aperture, which is the safe direction. The Assumptions record
    flags that when aperture_average is true.
    '''
    rx = scenario.rx_terminal
    # Open the collection context around the PHYSICS CALLS only.
    with trace_assumptions() as trace:
        sigma2_R = plane_wave_scintillation_index(geometry.elevation_deg,
                                                  rx.wavelength_m, hs, cn2_profile)
        if np.ndim(sigma2_R) != 0:
            raise NotImplementedError(
                "the gamma_gamma model takes a scalar elevation only. The "
                "gamma-gamma quantile and sampler carry one (alpha, beta) pair "
                "(olb.turbulence.andrews.distributions). Loop over the elevations."
            )
        sigma2_R = float(sigma2_R)
        sigma2_lnX = float(large_scale_log_variance(sigma2_R, wave='plane'))
        sigma2_lnY = float(small_scale_log_variance(sigma2_R, wave='plane'))
        alpha, beta = gamma_gamma_params(sigma2_lnX, sigma2_lnY)
        alpha, beta = float(alpha), float(beta)
        sigma2_I = float(gamma_gamma_scintillation_index(alpha, beta))

    # State the headline fields explicitly: this is the extended-Rytov strong
    # model, so REGIME_STRONG overrides the weak regime of the plane-wave feeder.
    assumptions = trace.merge(
        beam_type=BEAM_PLANE_WAVE,
        turbulence_regime=REGIME_STRONG,
        spectrum=SPECTRUM_KOLMOGOROV,
        validity="All fluctuation strengths. Andrews and Phillips, 2nd ed. "
                 "(2005), DOI 10.1117/3.626196, state that the extended Rytov "
                 "chain of Ch. 12, Eq. (40), printed p. 497, holds for every "
                 "value of the Rytov variance. The chain uses the Kolmogorov "
                 "spectrum only (Ch. 12, Eq. (15), printed p. 490), so it has "
                 "no inner scale and no outer scale. The receiver is a POINT.",
    )
    # The point-receiver limit is a scenario-level fact (the physics never sees
    # the receive aperture here), so flag it at the factory level.
    if aperture_average:
        assumptions.flag(
            "aperture averaging is NOT applied: the book gives no "
            "aperture-averaged downlink index in the moderate-to-strong regime "
            "(Ch. 12, Eq. (39) is weak theory, Ch. 12, Eq. (40) is a point "
            "form). This point-receiver fade is deeper than the true aperture "
            "fade.",
            source="factory:links.downlink",
        )

    return irradiance_fade_term(
        "scintillation", "turbulence",
        mean_log=gamma_gamma_mean_log(alpha, beta),
        quantile=lambda p: gamma_gamma_quantile(p, alpha, beta),
        rvs=lambda n, rng: gamma_gamma_rvs(n, alpha, beta, rng),
        note="plane-wave downlink gamma-gamma scintillation, point receiver",
        meta={
            "model": "gamma_gamma",
            "alpha": alpha,
            "beta": beta,
            "sigma2_lnX": sigma2_lnX,
            "sigma2_lnY": sigma2_lnY,
            "sigma2_I": sigma2_I,
            "sigma2_P": sigma2_I,
            "aperture_averaging_factor": 1.0,
            "sigma2_R": sigma2_R,
            # The gamma-gamma Term is valid at ALL strengths, so this flag is
            # purely informational: is the case actually in the weak REGIME? Test
            # the true Rytov variance against the REGIME boundary (1.0), NOT the
            # lognormal-PDF house rule (0.25) -- that was the factor-of-4 error of
            # Conflict C-05. sigma2_R here is a plane wave, so Lambda is None.
            "rytov_regime": rytov_weak(float(np.max(sigma2_R))),
            "weak_fluctuation_valid": bool(np.all(np.asarray(sigma2_R)
                                                  < WEAK_REGIME_LIMIT)),
            "weak_fluctuation_limit": WEAK_REGIME_LIMIT,
        },
        assumptions=assumptions,
    )


def _downlink_fidelity2_terms(scenario, geometry, wave, hs, cn2_profile,
                              turbulence=True):
    '''
    The fidelity-2 wave-optics Terms of a downlink.

    A space link cannot be simulated end to end (the sim propagates only the
    ~20 km atmosphere slab with a plane-wave input; the full slant path is
    absent), so the loss splits into:
      - a DETERMINISTIC geometric-loss Term (geometric spread + aperture capture
        + launch truncation). By default this is the ANALYTIC far-field Term
        (wave.vacuum is None), because a space link is far field and the wave
        vacuum run is slow and grid-noise-limited over the full slant range (see
        run_fidelity2). With vacuum="wave" it is the wave-optics vacuum Term
        (which also carries the vacuum fibre coupling for an SMF receiver).
      - one or two STOCHASTIC Terms (the slab turbulence penalty, vacuum-limit
        1.0).
    Together they replace the analytic geometric and scintillation/coupling Terms.
    `wave` is a Fidelity2Bundle from olb.models.waveoptics.run_fidelity2.

    An SMF receiver gets the composite fibre penalty (aperture-power penalty x
    ABSOLUTE slab fibre coupling smf_eta). In the wave-vacuum branch the static
    co-moving coupling is cancelled against the vacuum Term (vac_smf_db); in the
    analytic branch the analytic geometric Term carries no coupling, so the
    absolute smf_eta stands alone. An MMF (light-bucket) receiver gets TWO
    stochastic Terms: the aperture-power penalty, and one MMF coupling Term. The
    MMF coupling reads the ABSOLUTE core-capture (mmf_eta), so it holds the
    static encircled-energy floor, and NO vacuum baseline is subtracted. An
    Aperture / no-detector receiver gets the aperture-power penalty alone.

    With `turbulence` False, or with a VACUUM-ONLY bundle (turbulent None), the
    Term set keeps the DETERMINISTIC geometric Terms alone. A fibre receiver
    then needs a wave vacuum run, because the default analytic geometric Term
    carries no coupling number (see the raise below).
    '''
    from ..models.waveoptics import (waveoptics_vacuum_term,
                                     waveoptics_turbulence_term,
                                     waveoptics_mmf_coupling_term,
                                     waveoptics_vacuum_mmf_term)
    from ..terminal import SMF, MMF
    rx = scenario.rx_terminal
    is_smf = isinstance(rx.detector, SMF)
    is_mmf = isinstance(rx.detector, MMF)
    elev = np.asarray(geometry.elevation_deg, dtype=float)
    if elev.size == 1:
        # A one-element array IS one line of sight, so accept it as a scalar.
        elev = elev.reshape(())
    if elev.ndim != 0:
        raise ValueError(
            "the fidelity-2 downlink takes a scalar elevation (one range per "
            "record). Loop over the elevations and build one bundle each."
        )
    vacuum_only = (not turbulence) or wave.turbulent is None
    if wave.turbulent is None and turbulence:
        raise ValueError(
            "the `wave` bundle is vacuum-only (turbulent is None), but the "
            "budget asks for turbulence. Run "
            "olb.models.waveoptics.run_fidelity2 WITHOUT turbulence=False, or "
            "pass turbulence=False to the budget."
        )
    if vacuum_only and (is_smf or is_mmf) and wave.vacuum is None:
        raise ValueError(
            "a fibre receiver at fidelity=2 with turbulence=False has no "
            "coupling number: the space bundle skipped the wave vacuum run "
            "(vacuum='analytic'), and the analytic geometric Term carries no "
            "coupling. Run run_fidelity2(vacuum='wave', turbulence=False) for "
            "the static coupling, or use an Aperture receiver. The fidelity-0 "
            "analytic coupling is NOT wired in here."
        )
    # The scintillation index is a turbulence quantity, so read it only when a
    # stochastic Term needs it.
    sigma2_I = (None if vacuum_only else
                float(plane_wave_scintillation_index(
                    float(elev), rx.wavelength_m, hs, cn2_profile)))

    if wave.vacuum is None:
        # ANALYTIC geometric loss (the default for a space link). The link is far
        # field, so the analytic Term is exact and the wave vacuum run is skipped
        # (it is slow and grid-noise-limited over the full slant range; see
        # olb.models.waveoptics.run_fidelity2 and validation/vacuum_loss). The
        # launch-truncation Term is opt-in, the same rule as the analytic budget.
        geo = [geometric_loss_term(scenario, geometry)]
        eff = tx_gaussian_efficiency_term(scenario, geometry)
        if eff.mean_db > TX_TRUNCATION_MIN_DB:
            geo.append(eff)
        # No wave coupling baseline to cancel: the slab smf_eta is the ABSOLUTE
        # fibre coupling, so it stands alone in the SMF penalty below.
        vac_smf_db = 0.0
    else:
        # The wave-optics vacuum Term (opt-in for space, vacuum="wave").
        geo = [waveoptics_vacuum_term(wave.vacuum, include_smf=is_smf,
                                      beam_type=BEAM_PLANE_WAVE)]
        vac_smf_db = wave.vacuum.smf_coupling_db if is_smf else 0.0
    if vacuum_only:
        # VACUUM-ONLY: the deterministic Terms alone. An SMF coupling already
        # sits inside the wave vacuum Term (include_smf); an MMF light bucket
        # needs its own deterministic core-capture Term.
        if is_mmf:
            geo.append(waveoptics_vacuum_mmf_term(wave.vacuum, rx.detector,
                                                  rx.aperture_m,
                                                  beam_type=BEAM_PLANE_WAVE))
        return geo
    trials = wave.turbulent.trials
    if is_smf:
        coll = np.array([t.collected_power for t in trials], dtype=float)
        eta = np.array([t.smf_eta for t in trials], dtype=float)
        # collected_power is the aperture-power penalty (vacuum-limit 1.0), and
        # the slab smf_eta is the ABSOLUTE fibre coupling. In the wave-vacuum
        # branch the vacuum Term charges +smf_coupling_db, so -vac_smf_db cancels
        # it and the static coupling is counted once. In the analytic branch
        # vac_smf_db is 0 (the analytic geometric Term carries no coupling), so
        # the absolute smf_eta stands alone.
        loss_db = (-10.0 * np.log10(coll) - 10.0 * np.log10(eta) - vac_smf_db)
        pen = waveoptics_turbulence_term(
            wave.turbulent, loss_db=loss_db, beam_type=BEAM_PLANE_WAVE,
            sigma2_I=sigma2_I,
            note="downlink turbulence penalty (wave optics): aperture-power and "
                 "fibre-coupling loss relative to the vacuum baseline.")
        return geo + [pen]
    if is_mmf:
        # The light bucket: the aperture-power penalty (the bucket scintillation)
        # plus ONE MMF coupling evaluation on the turbulent field. mmf_eta is the
        # ABSOLUTE core-capture, so it holds the static encircled-energy floor. No
        # vacuum coupling baseline is subtracted (unlike the SMF composite): the
        # vacuum Term carries no coupling here (include_smf is False for an MMF rx).
        pen = waveoptics_turbulence_term(
            wave.turbulent, quantity="collected_power", beam_type=BEAM_PLANE_WAVE,
            sigma2_I=sigma2_I,
            note="downlink scintillation (wave optics): light-bucket aperture-power "
                 "penalty, vacuum-normalised.")
        cpl = waveoptics_mmf_coupling_term(
            wave.turbulent, beam_type=BEAM_PLANE_WAVE, sigma2_I=sigma2_I,
            note="downlink MMF coupling (wave optics): core-capture of the "
                 "turbulent focused spot, absolute (holds the static floor).")
        return geo + [pen, cpl]
    # An Aperture / no-detector receiver: the aperture-power penalty alone.
    pen = waveoptics_turbulence_term(
        wave.turbulent, quantity="collected_power", beam_type=BEAM_PLANE_WAVE,
        sigma2_I=sigma2_I,
        note="downlink scintillation (wave optics): aperture-power penalty, "
             "vacuum-normalised.")
    return geo + [pen]


def _auto_select(scenario, geometry, *, aperture_average, hs, cn2_profile):
    '''
    Selector layer: choose the downlink scintillation model from the regime.

    The selector reads the POINT scintillation index sigma2_I, which is the
    slant plane-wave Rytov variance. It returns:
        sigma2_I <  LOGNORMAL_PDF_LIMIT   the lognormal Term
        sigma2_I >= LOGNORMAL_PDF_LIMIT   the gamma-gamma Term

    The switch point is the lognormal-PDF house rule 0.25, DELIBERATELY tighter
    than the regime boundary sigma_R^2 = 1 (Andrews and Phillips, 2nd ed. (2005),
    DOI 10.1117/3.626196, Ch. 5, Eq. (15), printed p. 140; Ch. 12, Eq. (40),
    printed p. 497). This is NOT the factor-of-4 conflation: the switch is a
    PDF-fidelity decision, not a regime test. olb switches early because Ch. 11,
    Sec. 11.3, printed p. 451, says the lognormal tail is too thin, and this
    selector reports fade depths from that tail. The gamma-gamma chain of
    Ch. 12, Eq. (40) is valid at every fluctuation strength, so the early switch
    costs no validity. See `LOGNORMAL_PDF_LIMIT` and Conflict C-05 in
    docs/andrews-crosscheck.md.

    The gamma-gamma Term takes a scalar elevation only. For an elevation array
    that breaks the limit, the selector keeps the lognormal Term and warns.
    '''
    # Read the point index first, so that the lognormal Term is never built
    # (and never warns) for a case that goes to the gamma-gamma Term.
    sigma2_I = plane_wave_scintillation_index(geometry.elevation_deg,
                                              scenario.rx_terminal.wavelength_m,
                                              hs, cn2_profile)
    weak = np.all(np.asarray(sigma2_I) < LOGNORMAL_PDF_LIMIT)
    if not weak and np.ndim(sigma2_I) != 0:
        warnings.warn(
            "auto selector: sigma2_I exceeds the lognormal-PDF limit at one "
            "elevation or more, but the gamma-gamma Term takes a scalar "
            "elevation only. The lognormal Term is returned. Loop over the "
            "elevations to get the gamma-gamma fade."
        )
        weak = True
    if weak:
        return _lognormal_term(scenario, geometry,
                               aperture_average=aperture_average, hs=hs,
                               cn2_profile=cn2_profile)
    return _gamma_gamma_term(scenario, geometry,
                             aperture_average=aperture_average, hs=hs,
                             cn2_profile=cn2_profile)


_MODELS = {
    "lognormal": _lognormal_term,
    "gamma_gamma": _gamma_gamma_term,
    "auto": _auto_select,
}


def downlink_scintillation_term(scenario, geometry, *, model="lognormal",
                                aperture_average=True, hs=None, cn2_profile=None):
    '''
    Build the analytic downlink scintillation Term (fidelity 0/1, aperture).

    Dispatch to the requested model. The "lognormal" model is the analytic
    weak-fluctuation model. The "gamma_gamma" model is the analytic
    moderate-to-strong model. The "auto" model is the selector layer. The
    fidelity-2 wave-optics downlink is NOT a scintillation-model choice: it is the
    whole-path `fidelity=2` route of downlink_budget (two Terms).

    Parameters:
        scenario : SpaceScenario
            Reads the receive terminal wavelength and aperture. Reads the site to
            build the default Cn2 profile with get_c2n.
        geometry : CircularOrbit or TLEPass
            Reads elevation_deg. A scalar elevation gives a scalar Term. An
            elevation array gives a Term that broadcasts over that shape.
        model : str
            One of "lognormal", "gamma_gamma", "auto".
        aperture_average : bool
            Apply the plane-wave aperture-averaging factor when true.
        hs : numpy.ndarray, optional
            Heights above the ground station [m]. Defaults to DEFAULT_HS.
        cn2_profile : numpy.ndarray, optional
            Zenith Cn2(h) profile on the hs grid. Defaults to the site profile
            from get_c2n.

    Returns:
        Term
            name="scintillation", category="turbulence".

    Raises:
        ValueError
            If model is not a known name.
        NotImplementedError
            If model is "gamma_gamma" with an elevation array.
    '''
    if model not in _MODELS:
        raise ValueError(
            f"unknown model {model!r}. Use one of {sorted(_MODELS)}."
        )
    hs = DEFAULT_HS if hs is None else hs
    if cn2_profile is None:
        cn2_profile = default_cn2_profile(scenario.channel.site, hs)
    return _MODELS[model](scenario, geometry, aperture_average=aperture_average,
                          hs=hs, cn2_profile=cn2_profile)


def downlink_budget(scenario, geometry, *, fidelity=1, tau_zenith=None,
                    scintillation=True, turbulence=True, n_samples=2000,
                    fast_params=None, scint_model="lognormal", wave=None):
    '''
    Assemble the downlink budget at a chosen fidelity.

    `fidelity` is a WHOLE-PATH choice (see the README fidelity ladder):

      - fidelity=0 (analytic). An SMF detector gets the mean-only analytic
        fibre-coupling Term (downlink_coupling_term(smf_fidelity="mean")). An
        Aperture / no detector gets the analytic scintillation Term
        (`scint_model`).
      - fidelity=1 (the default, statistical). An SMF detector gets the FAST
        modal-overlap Term (downlink_coupling_term(smf_fidelity="fast"), needs
        fast-aosim). An Aperture / no detector uses the SAME analytic
        scintillation Term as fidelity 0 (the closed-form lognormal / gamma-gamma
        is the model of record and already carries a fade — tiers 0 and 1
        coincide for an aperture; only the SMF coupling model changes).
      - fidelity=2 (wave optics). The whole path is a field simulation. It gives
        a deterministic vacuum-optics Term (geometric spread + aperture capture +
        vacuum fibre coupling over the full slant range) and a stochastic
        turbulence Term (the slab penalty). An SMF or Aperture receiver gives TWO
        Terms; an MMF (light-bucket) receiver gives THREE (the vacuum Term carries
        no coupling, and the aperture-power penalty and the absolute MMF coupling
        are two Terms). They REPLACE the geometric and the scintillation /
        coupling Terms. Only the analytic extinction and pointing Terms stay. It
        needs a precomputed `wave` bundle (olb.models.waveoptics.run_fidelity2);
        the budget never runs the sim.

    The geometric, extinction, and pointing Terms are the deterministic backbone
    at fidelity 0/1. The downlink keeps a standalone pointing Term because it has
    no coupled-flux beam-wave machinery to fold jitter into (unlike the uplink).

    Parameters:
        scenario : SpaceScenario
            The link case.
        geometry : CircularOrbit or TLEPass
            The link geometry.
        fidelity : int
            0 (analytic), 1 (statistical, the default), or 2 (wave optics, needs
            `wave`).
        tau_zenith : float, optional
            Zenith optical depth. Defaults to extinction.DEFAULT_TAU_ZENITH.
        scintillation : bool
            Add the analytic scintillation Term for an aperture / no-detector
            receiver at fidelity 0/1 when true.
        turbulence : bool
            Master turbulence switch, at every fidelity. At fidelity 0/1, drop
            every turbulence quantity and keep the static parts. At fidelity 2,
            drop the wave-optics turbulence Term and the stochastic coupling
            Term, so the budget keeps the deterministic geometric Terms plus
            extinction and pointing. Pair it with
            olb.models.waveoptics.run_fidelity2(turbulence=False), which makes
            no screens and no trials; the EMPTY bundle (vacuum=None,
            turbulent=None) of a space link is accepted. A `wave` bundle is
            still required at fidelity 2, so the call shape is uniform.
        n_samples : int
            FAST Monte Carlo draws (NITER) for the fidelity-1 SMF coupling.
        fast_params : dict, optional
            Extra FAST parameters for the fidelity-1 SMF coupling.
        scint_model : str
            The analytic scintillation MODEL for an APERTURE receiver:
            "lognormal" (the default), "gamma_gamma", or "auto". It is not a
            fidelity axis; it applies at fidelity 0/1 only.
        wave : Fidelity2Bundle, list, or Campaign, optional
            The precomputed wave-optics record for fidelity=2: a
            Fidelity2Bundle, a list of them, or a Campaign. Run it with
            olb.models.waveoptics.run_fidelity2, or store it with
            olb.waveoptics.turbulence.Campaign and pass the campaign itself
            (olb.models.waveoptics.resolve_wave turns it into the bundle).

    A receive detector is opt-in. A BUCKET receiver -- no detector or a plain
    Aperture -- gets the aperture-averaged scintillation Term (with the
    scint_model selector); None and Aperture() are the SAME bucket, so there is
    no "no detector" special case. An SMF detector instead gets the
    receive-coupling Term, which owns the receive-side turbulence physics and
    REPLACES the scintillation Term. An MMF or Camera receiver raises at
    fidelity 0/1.

    Returns:
        Budget
            The budget with the scenario set.

    Raises:
        ValueError
            If fidelity is not 0/1/2, or if fidelity=2 without a `wave` bundle.
    '''
    if fidelity not in (0, 1, 2):
        raise ValueError(f"fidelity must be 0, 1, or 2, got {fidelity!r}.")
    tau = DEFAULT_TAU_ZENITH if tau_zenith is None else tau_zenith

    if fidelity == 2:
        # A Campaign is a wave record too: turn it into the bundle it holds.
        from ..models.waveoptics import resolve_wave
        wave = resolve_wave(wave)
        if wave is None:
            raise ValueError(
                "fidelity=2 needs a precomputed `wave` bundle. Run "
                "olb.models.waveoptics.run_fidelity2(scenario, geometry, ...) and "
                "pass it as wave. The budget does not run the split-step "
                "propagation implicitly."
            )
        # The two wave-optics Terms replace geometric and scintillation/coupling.
        # Only extinction (absorption) and pointing (mechanical jitter) stay.
        hs = DEFAULT_HS
        cn2_profile = default_cn2_profile(scenario.channel.site, hs)
        terms = [
            slant_extinction_term(scenario, geometry, tau_zenith=tau),
            pointing_loss_term(scenario, geometry),
        ]
        terms += _downlink_fidelity2_terms(scenario, geometry, wave, hs,
                                           cn2_profile, turbulence=turbulence)
        return Budget(terms, scenario=scenario)

    # fidelity 0/1: the analytic backbone plus the receive-side turbulence Term.
    smf_fidelity = "fast" if fidelity == 1 else "mean"
    terms = [
        geometric_loss_term(scenario, geometry),
        slant_extinction_term(scenario, geometry, tau_zenith=tau),
        pointing_loss_term(scenario, geometry),
    ]
    terminal = getattr(scenario, "rx_terminal", None)
    detector = terminal.detector if terminal is not None else None
    # A bucket receiver -- no detector or a plain Aperture -- is phase-insensitive,
    # so its turbulence penalty is the aperture-averaged scintillation Term (with
    # the scint_model selector). None and Aperture() are the SAME bucket: there is
    # no "no detector" case, that is just Aperture(). An SMF/MMF/Camera detector
    # takes the receive-coupling Term instead (SMF couples; MMF and Camera raise).
    from ..terminal import Aperture
    if detector is not None and not isinstance(detector, Aperture):
        # Import here to break the downlink <-> coupling import cycle.
        from ..models.coupling import downlink_coupling_term
        terms.append(downlink_coupling_term(scenario, geometry, n_samples=n_samples,
                                      smf_fidelity=smf_fidelity,
                                      fast_params=fast_params,
                                      turbulence=turbulence))
    elif scintillation and turbulence:
        terms.append(downlink_scintillation_term(scenario, geometry,
                                                 model=scint_model,
                                                 aperture_average=True))
    return Budget(terms, scenario=scenario)


if __name__ == '__main__':
    from ..scenario import SpaceScenario, Channel
    from ..geometry import CircularOrbit
    from ..terminal import Terminal, Transmitter, Aperture, SMF, TipTilt, AO

    lam = 1550e-9

    def _dl(ground, *, jitter=0.0, power=None):
        '''Build a downlink SpaceScenario: tx=space (satellite waist 0.035), rx=ground.'''
        space = Terminal(aperture_m=0.05, wavelength_m=lam,
                         pointing_jitter_rad=jitter,
                         transmitter=Transmitter(waist_m=0.035, power_dbm=power))
        return SpaceScenario(ground=ground, space=space, direction="downlink",
                        channel=Channel(altitude_m=600e3))

    scenario = _dl(Terminal(aperture_m=0.7, wavelength_m=lam))
    hs = DEFAULT_HS
    cn2 = default_cn2_profile(scenario.channel.site, hs)

    # The lognormal Term gives a working quantile deeper than the mean loss.
    geom = CircularOrbit(600e3, 30.0)
    term = downlink_scintillation_term(scenario, geom, cn2_profile=cn2)
    assert term.name == "scintillation", term.name
    q99 = term.quantile_db(0.99)
    assert q99 is not None and np.isfinite(q99) and q99 > term.mean_db, (q99, term.mean_db)
    assert term.assumptions is not None

    # WP3a: the lognormal factory now INHERITS its assumptions from the traced
    # physics. The Term carries the plane-wave scintillation source and (when
    # aperture-averaged) the aperture-averaging source; the provenance is not
    # empty, so Budget.check()'s untraced guard stays quiet.
    prov = term.assumptions.provenance
    assert prov, "the lognormal Term must carry traced provenance"
    assert any("plane_wave_scintillation_index" in s for s in prov), prov
    assert any("aperture_averaged_scintillation_index" in s for s in prov), prov

    # A 20 deg elevation breaks the weak-fluctuation limit; 60 deg does not.
    import warnings as _w
    with _w.catch_warnings():
        _w.simplefilter("ignore")
        t20 = downlink_scintillation_term(scenario, CircularOrbit(600e3, 20.0),
                                          cn2_profile=cn2)
        t60 = downlink_scintillation_term(scenario, CircularOrbit(600e3, 60.0),
                                          cn2_profile=cn2)
    assert not t20.assumptions.ok
    assert t60.assumptions.ok

    # The sampled mean matches mean_db for large n.
    rng = np.random.default_rng(0)
    sampled = term.sample_db(200_000, rng)
    assert sampled.shape == (200_000,)
    assert abs(sampled.mean() - term.mean_db) < 0.02, (sampled.mean(), term.mean_db)

    # Fidelity 2 is now a WHOLE-PATH budget route, not a scintillation model.
    # "montecarlo" is gone from downlink_scintillation_term.
    try:
        downlink_scintillation_term(scenario, geom, model="montecarlo",
                                    cn2_profile=cn2)
    except ValueError as e:
        assert "unknown model" in str(e)
    else:
        raise AssertionError("model='montecarlo' must raise unknown model")

    # Guards (no run): bad fidelity, and fidelity=2 with no bundle.
    for bad in (3, "fast"):
        try:
            downlink_budget(scenario, geom, fidelity=bad)
        except ValueError as e:
            assert "fidelity must be 0, 1, or 2" in str(e)
        else:
            raise AssertionError(f"fidelity={bad!r} must raise")
    try:
        downlink_budget(scenario, geom, fidelity=2)
    except ValueError as e:
        assert "needs a precomputed `wave` bundle" in str(e)
    else:
        raise AssertionError("fidelity=2 without a bundle must raise")
    # The default aperture budget is UNCHANGED: fidelity=1 aperture = lognormal.
    def_scint = next(t for t in downlink_budget(scenario, geom).terms
                     if t.category == "turbulence")
    assert def_scint.meta["model"] == "lognormal"

    # A real fidelity-2 aperture downlink (skip if aotools absent). The DEFAULT
    # geometric loss is ANALYTIC (a space link is far field, so wave.vacuum is
    # None), so the budget shows the analytic "geometric spreading" Term plus the
    # wave-optics turbulence Term, and a real fade.
    from ..models.waveoptics import run_fidelity2
    try:
        import warnings as _w2
        with _w2.catch_warnings():
            _w2.simplefilter("ignore")
            f2_bundle = run_fidelity2(scenario, CircularOrbit(600e3, 30.0),
                                      preset="rapid", n_trials=16, seed=5,
                                      hs=hs, cn2_profile=cn2, progress=False)
            f2 = downlink_budget(scenario, CircularOrbit(600e3, 30.0), fidelity=2,
                                 wave=f2_bundle)
    except ImportError:
        print("aotools not installed; skipping the downlink fidelity-2 run.")
        f2 = None
    if f2 is not None:
        assert f2_bundle.vacuum is None, "space defaults to the analytic vacuum"
        geo = next(t for t in f2.terms if t.name == "geometric spreading")
        turb = next(t for t in f2.terms if t.meta.get("model") == "waveoptics")
        assert geo.category == "geometric" and not geo.stochastic and turb.stochastic
        assert not any(t.meta.get("model") == "waveoptics-vacuum" for t in f2.terms)
        with _w2.catch_warnings():
            _w2.simplefilter("ignore")
            assert f2.provides_fade and np.isfinite(f2.fade_margin_db(0.9))
        print(f"downlink fidelity 2 (600 km, 30 deg, rapid, 16 trials): analytic "
              f"geometry {geo.mean_db:.2f} dB + turbulence {turb.mean_db:.3f} dB")

        # The wave vacuum stays available as an OPT-IN (vacuum="wave"): the
        # budget then shows the wave-optics vacuum Term instead of the analytic
        # geometric Term.
        with _w2.catch_warnings():
            _w2.simplefilter("ignore")
            w_bundle = run_fidelity2(scenario, CircularOrbit(600e3, 30.0),
                                     preset="rapid", n_trials=16, seed=5, hs=hs,
                                     cn2_profile=cn2, progress=False, vacuum="wave")
            fw = downlink_budget(scenario, CircularOrbit(600e3, 30.0), fidelity=2,
                                 wave=w_bundle)
        assert w_bundle.vacuum is not None
        vacw = next(t for t in fw.terms
                    if t.meta.get("model") == "waveoptics-vacuum")
        assert not vacw.stochastic
        print(f"  opt-in wave vacuum: {vacw.mean_db:.2f} dB "
              f"(vs analytic {geo.mean_db:.2f} dB)")

    # A fidelity-2 MMF (light-bucket) downlink gives THREE loss Terms beside the
    # extinction and pointing Terms: the DEFAULT analytic geometric spreading,
    # the aperture-power scintillation, and the MMF coupling. The MMF coupling
    # reads the ABSOLUTE mmf_eta, so it holds the static floor and no vacuum
    # baseline is subtracted.
    from ..terminal import MMF
    scn_mmf = _dl(Terminal(aperture_m=0.5, wavelength_m=lam,
                           detector=MMF(core_radius_m=25e-6, optimal_focus=True,
                                        numerical_aperture=0.2,
                                        sensitivity_dbm=-110)),
                  power=30)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            mmf_bundle = run_fidelity2(
                scn_mmf, CircularOrbit(600e3, 30.0), preset="rapid", n_trials=16,
                seed=5, hs=hs, progress=False,
                cn2_profile=default_cn2_profile(scn_mmf.channel.site, hs))
            f2_mmf = downlink_budget(scn_mmf, CircularOrbit(600e3, 30.0),
                                     fidelity=2, wave=mmf_bundle)
    except ImportError:
        f2_mmf = None
    if f2_mmf is not None:
        cpl_m = next(t for t in f2_mmf.terms if t.name == "receive coupling (MMF)")
        geo_m = next(t for t in f2_mmf.terms if t.name == "geometric spreading")
        # The analytic geometric Term carries NO coupling; no vacuum baseline.
        assert mmf_bundle.vacuum is None and geo_m.category == "geometric"
        assert not any(t.meta.get("model") == "waveoptics-vacuum"
                       for t in f2_mmf.terms)
        assert cpl_m.category == "coupling" and cpl_m.stochastic and not cpl_m.mean_only
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            assert f2_mmf.provides_fade and np.isfinite(f2_mmf.fade_margin_db(0.9))
        print(f"downlink fidelity 2 MMF (600 km, 30 deg): analytic geometry "
              f"{geo_m.mean_db:.2f} dB + MMF coupling {cpl_m.mean_db:.2f} dB")

    # --- fidelity-2 master turbulence switch (no simulation needed) ----------
    # The EMPTY bundle (vacuum=None, turbulent=None) of a space link needs no
    # run at all, so this case is cheap. The budget then shows the analytic
    # deterministic Terms alone.
    from ..models.waveoptics import Fidelity2Bundle
    empty = Fidelity2Bundle(vacuum=None, turbulent=None)
    f2_off = downlink_budget(scenario, CircularOrbit(600e3, 30.0), fidelity=2,
                             wave=empty, turbulence=False)
    off_names = [t.name for t in f2_off.terms]
    off_cats = [t.category for t in f2_off.terms]
    assert "turbulence" not in off_cats and "coupling" not in off_cats, off_cats
    assert "geometric spreading" in off_names, off_names
    assert "atmospheric" in off_cats and "pointing" in off_cats, off_cats
    # An empty bundle with turbulence=True raises a helpful error.
    try:
        downlink_budget(scenario, CircularOrbit(600e3, 30.0), fidelity=2,
                        wave=empty)
    except ValueError as e:
        assert "vacuum-only" in str(e), str(e)
    else:
        raise AssertionError("an empty bundle with turbulence must raise")
    # A fibre receiver has no vacuum coupling number in the analytic-default
    # bundle, so it raises and names the fix.
    scn_smf_off = _dl(Terminal(aperture_m=0.7, wavelength_m=lam,
                               detector=SMF(sensitivity_dbm=-40)))
    try:
        downlink_budget(scn_smf_off, CircularOrbit(600e3, 30.0), fidelity=2,
                        wave=empty, turbulence=False)
    except ValueError as e:
        assert "vacuum='wave'" in str(e), str(e)
    else:
        raise AssertionError("an SMF receiver with no wave vacuum must raise")
    # A `wave` bundle is still REQUIRED, so the call shape stays uniform.
    try:
        downlink_budget(scenario, CircularOrbit(600e3, 30.0), fidelity=2,
                        turbulence=False)
    except ValueError as e:
        assert "needs a precomputed `wave` bundle" in str(e)
    else:
        raise AssertionError("fidelity=2 turbulence=False still needs a bundle")
    print(f"downlink fidelity 2, turbulence=False (600 km, 30 deg): "
          f"{f2_off.total_loss_db():.2f} dB, terms {off_names}")

    # --- gamma-gamma Term ---------------------------------------------------
    # A 15 deg elevation is a strong case: the point index passes the house
    # limit, so the auto selector must route to the gamma-gamma Term.
    geom15 = CircularOrbit(600e3, 15.0)
    gg = downlink_scintillation_term(scenario, geom15, model="auto",
                                     cn2_profile=cn2)
    assert gg.meta["model"] == "gamma_gamma", gg.meta["model"]
    # WP3a: the gamma-gamma factory inherits the traced extended-Rytov chain.
    gg_prov = gg.assumptions.provenance
    assert gg_prov, "the gamma-gamma Term must carry traced provenance"
    assert any("plane_wave_scintillation_index" in s for s in gg_prov), gg_prov
    assert any("large_scale_log_variance" in s for s in gg_prov), gg_prov
    assert any("small_scale_log_variance" in s for s in gg_prov), gg_prov
    # The case passed the lognormal-PDF switch (sigma2_I >= LOGNORMAL_PDF_LIMIT),
    # so it routed to gamma-gamma; sigma2_R is comfortably past 0.25 too.
    assert gg.meta["sigma2_I"] >= LOGNORMAL_PDF_LIMIT, gg.meta["sigma2_I"]
    gg99 = gg.quantile_db(0.99)
    assert gg99 > gg.mean_db > 0.0, (gg99, gg.mean_db)
    assert not gg.mean_only
    # The Term reproduces the book weak-to-strong index, Ch. 12, Eq. (40),
    # printed p. 497, which andrews.paths gives on the same slant path.
    from ..turbulence.andrews.paths import downlink_scintillation_index
    # The two feeders differ only by the ground datum: andrews.paths integrates
    # (h - h0)^(5/6) and plane_wave_scintillation.py integrates h^(5/6), and
    # DEFAULT_HS starts at h0 = 1 m. So this is a MEASUREMENT, not an identity.
    book = downlink_scintillation_index(hs, cn2, lam, 15.0, regime="strong")
    gg_gap = abs(gg.meta["sigma2_I"] / book - 1.0)
    assert gg_gap < 1e-2, (gg.meta["sigma2_I"], book)
    # The sampled mean of the dB loss matches mean_db.
    gg_draws = gg.sample_db(200_000, np.random.default_rng(3))
    assert abs(gg_draws.mean() - gg.mean_db) < 0.02, (gg_draws.mean(), gg.mean_db)
    # A weak case still routes to the lognormal Term.
    assert downlink_scintillation_term(scenario, CircularOrbit(600e3, 60.0),
                                       model="auto",
                                       cn2_profile=cn2).meta["model"] == "lognormal"
    # The gamma-gamma fade is deeper than the lognormal fade of the same case.
    with _w.catch_warnings():
        _w.simplefilter("ignore")
        ln15 = downlink_scintillation_term(scenario, geom15, model="lognormal",
                                           aperture_average=False,
                                           cn2_profile=cn2)
    assert gg99 > ln15.quantile_db(0.99), (gg99, ln15.quantile_db(0.99))
    # An elevation array is refused by the gamma-gamma model itself.
    try:
        downlink_scintillation_term(scenario,
                                    CircularOrbit(600e3, np.array([15.0, 20.0])),
                                    model="gamma_gamma", cn2_profile=cn2)
    except NotImplementedError:
        pass
    else:
        raise AssertionError("gamma_gamma must refuse an elevation array")

    # An elevation sweep broadcasts over the geometry shape.
    sweep = CircularOrbit(600e3, np.array([20.0, 45.0, 90.0]))
    sweep_term = downlink_scintillation_term(scenario, sweep, cn2_profile=cn2)
    assert np.shape(sweep_term.mean_db) == (3,)
    assert sweep_term.sample_db(1000, rng).shape == (1000, 3)

    # A central obscuration flags a filter-validity violation at 60 deg (weak
    # fluctuation holds). A zero obscuration does not flag it.
    scen_obs = _dl(Terminal(aperture_m=0.7, wavelength_m=lam,
                            obscuration_ratio=0.3))
    cn2_obs = default_cn2_profile(scen_obs.channel.site, hs)
    t_obs = downlink_scintillation_term(scen_obs, CircularOrbit(600e3, 60.0),
                                        cn2_profile=cn2_obs)
    assert not t_obs.assumptions.ok, "obscuration=0.3 must flag a violation"
    assert t60.assumptions.ok, "obscuration=0.0 must not flag a violation"

    # --- downlink budget self-check -----------------------------------------
    budget_scn = _dl(Terminal(aperture_m=0.7, wavelength_m=lam),
                     jitter=2e-6, power=40)
    down = downlink_budget(budget_scn, CircularOrbit(altitude_m=600e3,
                                                    elevation_deg=60.0))
    assert down.to_frame().shape[0] == 4, down.to_frame().shape
    # Every downlink Term has a closed-form quantile, so analytic fade works.
    down_fade = down.fade_margin_db(0.99)
    assert np.isfinite(down_fade), down_fade
    # The plain downlink budget shows the base Term name "scintillation".
    assert "scintillation" in [t.name for t in down.terms], [t.name for t in down.terms]
    # WP3a: the demo budget's turbulence Term carries traced provenance, so
    # Budget.check() reports NO untraced-guard entry.
    guard = down.check(warn=False)
    assert not any("did not open" in reason for _, reason in guard), guard

    # --- opt-in receive terminal --------------------------------------------
    # The no-detector budget is unchanged: 4 terms, base scintillation name.
    assert scenario.rx_terminal.detector is None
    assert down.to_frame().shape[0] == 4
    assert "scintillation" in [t.name for t in down.terms]

    geom60 = CircularOrbit(altitude_m=600e3, elevation_deg=60.0)

    # Aperture detector: a bucket is the SAME as no detector, so it gets the
    # scintillation Term, byte-for-byte identical to the plain (no-detector) total.
    scn_ap = _dl(Terminal(aperture_m=0.7, wavelength_m=lam,
                          detector=Aperture(sensitivity_dbm=-40)),
                 jitter=2e-6, power=40)
    down_ap = downlink_budget(scn_ap, geom60)
    assert down_ap.to_frame().shape[0] == 4                 # same count as plain
    assert "scintillation" in [t.name for t in down_ap.terms]
    assert np.isclose(down_ap.total_loss_db(), down.total_loss_db()), (
        down_ap.total_loss_db(), down.total_loss_db())     # parity

    # SMF detector, no AO: the coupling loss deepens the total over the aperture.
    scn_smf = _dl(Terminal(aperture_m=0.7, wavelength_m=lam,
                           detector=SMF(sensitivity_dbm=-40)),
                  jitter=2e-6, power=40)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        down_smf = downlink_budget(scn_smf, geom60)
    assert "receive coupling (SMF)" in [t.name for t in down_smf.terms]
    assert down_smf.total_loss_db() > down_ap.total_loss_db()

    # SMF detector with tip-tilt + AO: the coupling loss shrinks toward the
    # aperture case.
    scn_ao = _dl(Terminal(aperture_m=0.7, wavelength_m=lam,
                          detector=SMF(sensitivity_dbm=-40),
                          compensation=[TipTilt(), AO(200)]),
                 jitter=2e-6, power=40)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        down_ao = downlink_budget(scn_ao, geom60)
    assert down_ao.total_loss_db() < down_smf.total_loss_db()
    # The SMF budget still has a closed-form analytic fade (coupling + fade).
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        assert np.isfinite(down_ao.fade_margin_db(0.99))

    # --- master turbulence switch (turbulence=False) ------------------------
    # An Aperture detector is a bucket (the SAME as no detector), so turbulence off
    # drops the scintillation Term entirely. The budget keeps its fade from the
    # pointing jitter.
    ap_off = downlink_budget(scn_ap, geom60, turbulence=False)
    assert not any(t.category in ("coupling", "turbulence") for t in ap_off.terms), \
        [t.name for t in ap_off.terms]
    assert ap_off.provides_fade and np.isfinite(ap_off.fade_margin_db(0.99))
    # An SMF detector with turbulence off keeps only the static mode-match loss
    # (deterministic, NOT mean-only), so the budget still provides a fade and the
    # off total is cheaper than the on total. No FAST run happens.
    smf_off = downlink_budget(scn_smf, geom60, turbulence=False)
    cpl_off = next(t for t in smf_off.terms if t.category == "coupling")
    assert cpl_off.meta["model"] == "static" and not cpl_off.mean_only
    assert np.isclose(cpl_off.mean_db, -10.0 * np.log10(0.8145))
    assert smf_off.provides_fade and smf_off.total_loss_db() < down_smf.total_loss_db()

    print(f"mean loss = {term.mean_db:.4f} dB   99% fade = {q99:.4f} dB")
    print(f"sampled mean = {sampled.mean():.4f} dB (n=200000)")
    print(f"auto @15deg: model={gg.meta['model']}  "
          f"sigma2_R={gg.meta['sigma2_R']:.4f}  "
          f"alpha={gg.meta['alpha']:.3f}  beta={gg.meta['beta']:.3f}  "
          f"sigma2_I={gg.meta['sigma2_I']:.4f}  "
          f"(book Ch. 12, Eq. (40) = {book:.4f}, {gg_gap * 100:+.3f} %)")
    print(f"gamma-gamma @15deg: mean loss = {gg.mean_db:.4f} dB   "
          f"99% fade = {gg99:.4f} dB   "
          f"(lognormal point 99% fade = {ln15.quantile_db(0.99):.4f} dB)")
    print(down.to_frame().to_string(index=False))
    print(f"terminal totals @60deg: aperture={down_ap.total_loss_db():.3f} dB  "
          f"SMF no-AO={down_smf.total_loss_db():.3f} dB  "
          f"SMF+AO200={down_ao.total_loss_db():.3f} dB")
    print("self-check passed")
