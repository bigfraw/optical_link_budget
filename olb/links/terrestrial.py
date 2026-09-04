'''
Terrestrial (horizontal-path) Terms and budget assembly.

This module builds the terrestrial link budget: a ground-to-ground horizontal
path. It reuses the shared Terms (geometric spreading, pointing
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
sigma2_I(0, L) of olb.turbulence.beam_wave_scintillation (Dios et al., Applied
Optics 43 (2004) 3866, Eq. 16), evaluated on a constant-Cn2 horizontal grid. The
aperture-averaging win is the Andrews weak-turbulence Kolmogorov factor of
olb.turbulence.plane_wave_scintillation.aperture_averaging_factor_weak (Andrews
and Phillips, 2nd ed. 2005, Ch. 10). The lognormal fade faces (mean_db, quantile,
sampler) follow the same closed form as the downlink lognormal Term (see
olb.links.downlink._lognormal_term).
'''

import warnings

import numpy as np

from ..results import Budget
from ..assumptions import (trace_assumptions, BEAM_GAUSSIAN, REGIME_WEAK,
                           SPECTRUM_KOLMOGOROV)
from ..models.geometric import geometric_loss_term
from ..models.extinction import terrestrial_extinction_term
from ..models.pointing import pointing_loss_term
from ..models.gaussian_efficiency import tx_gaussian_efficiency_term
from ..turbulence.beam_wave_scintillation import on_axis_scintillation_index
from ..turbulence.plane_wave_scintillation import aperture_averaging_factor_weak
from ..turbulence.andrews.beam import beam_params
from ..turbulence.andrews.scintillation import (rytov_weak, rytov_variance,
                                        LOGNORMAL_PDF_LIMIT,
                                        RYTOV_CONFIDENT_WEAK, WEAK_REGIME_LIMIT)
from ..turbulence.andrews.distributions import (lognormal_params,
                                                lognormal_mean_log,
                                                lognormal_quantile,
                                                lognormal_rvs)
from ..models.fade import irradiance_fade_term

# Below this launch-truncation loss the beam is an untruncated Gaussian, so the
# transmit Gaussian-efficiency term is skipped [dB]. Matches olb.links.uplink.
TX_TRUNCATION_MIN_DB = 1e-2

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
         olb.turbulence.beam_wave_scintillation.on_axis_scintillation_index (Dios et
         al., Applied Optics 43 (2004) 3866, Eq. 16). A horizontal path is a
         constant-Cn2 grid from the transmitter (z=0) to the receiver (z=L) at
         elevation 90 (sec=1), so the range L is the grid length.

      2. The aperture-averaging factor A that reduces the point index to the
         aperture-averaged flux index sigma2_P = A * sigma2_I. A telescope of
         diameter D averages the fine scintillation, so A falls as D grows. See
         olb.turbulence.plane_wave_scintillation.aperture_averaging_factor_weak (Andrews
         and Phillips, 2nd ed. 2005, Ch. 10, weak Kolmogorov, small inner scale).

    The flux index sigma2_P then feeds the same lognormal fade faces as the
    downlink lognormal Term (olb.links.downlink._lognormal_term). With
    sigma_l^2 = ln(1 + sigma2_P) the faces are
        mean_db     = (5/ln10) * sigma_l^2
        quantile(p) = -10*log10( exp(-sigma_l^2/2 + sigma_l * Phi_inv(1-p)) )
    Source of the lognormal irradiance PDF: Andrews and Phillips, 2nd ed. (2005),
    Ch. 5. Phi_inv is the inverse standard normal CDF.

    Validity: TWO separate weak-fluctuation tests, kept distinct (Conflict C-05,
    TL-05). The Term carries a `rytov_regime` label ("weak"/"soft"/"hard") from
    the beam-aware regime gate on the Rytov variance (both Ch. 5, Eq. (16)
    conditions, so a focused beam is caught), and a `weak_fluctuation_valid` flag
    from the tighter lognormal-PDF house rule sigma2_I < LOGNORMAL_PDF_LIMIT
    (0.25). A "soft" regime gives a soft warning; a "hard" regime or an invalid
    PDF flags the assumptions and warns. Above the weak regime use the fidelity-2
    Monte Carlo; above the PDF limit use gamma-gamma or Monte Carlo.

    Parameters:
        scenario : TerrestrialScenario
            The tx terminal's Transmitter waist launches the beam; the rx
            terminal's aperture diameter D and wavelength set the averaging.
            The scenario `direction` maps the roles (forward: tx = near,
            rx = far; reverse swaps them). Reads the channel path length L and
            the constant Cn2.
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

    # Open the collection context around the PHYSICS CALLS only. Each decorated
    # function registers its own assumptions (beam type, regime, spectrum, and
    # its constraints and checks), so the Term inherits the union. A strong path
    # trips the Dios beam-wave reliability check (sigma_chi^2 >= 0.6) inside
    # on_axis_scintillation_index automatically, so the hard-regime violation is
    # traced, not hand-built.
    with trace_assumptions() as trace:
        sigma2_I = float(on_axis_scintillation_index(
            hs, cn2_profile, w0, wavelength, elevation_deg=90.0, path_length_m=None))
        # Aperture-averaging win: a larger D averages more, so A and sigma2_P fall.
        A = float(aperture_averaging_factor_weak(D, wavelength, L))
        # Lambda (receiver-plane beam parameter, collimated launch) for the
        # regime label below.
        Lambda = float(beam_params(w0, wavelength, L).lam)

    # TODO(pointing jitter): this Term is ON-AXIS only (r=0), so it does NOT yet
    # carry pointing/tracking jitter. When adding it, note the ANALYTIC-PATH
    # ASYMMETRY vs the uplink MC (olb.turbulence.uplink_flux, which folds jitter
    # into r=beta and gets everything for free because it works in ABSOLUTE
    # irradiance):
    #   1. Fold the jitter displacement into the OFF-AXIS radial index at r=beta
    #      (beam_wave_scintillation.radial_scintillation_index), beta drawn from the
    #      jitter (+ beam-wander) 2-D variance. This adds the FLUCTUATION boost.
    #   2. You STILL need a SEPARATE mean-power loss term for the same jitter,
    #      because the off-axis sigma2_I is normalised to the LOCAL mean and does
    #      NOT carry the exp(-2*beta^2/W^2) mean drop. The MC path merges the two;
    #      the analytic path cannot. Adding both is NOT double-counting -- they
    #      are different statistical moments (variance vs mean).
    # This is the reason to converge the analytic and MC Dios paths rather than
    # patch one. See memory dios-scintillation-convergence / pointing-jitter-into-beta.

    sigma2_P = A * sigma2_I

    # The lognormal irradiance model. lognormal_params turns the aperture-averaged
    # index into the log-irradiance variance; the three dB faces come from the ONE
    # shared adapter (olb.models.fade.irradiance_fade_term), the SAME path the
    # downlink and gamma-gamma Terms use (backlog I-2 / crosscheck TL-01..04).
    sigma_l2 = lognormal_params(sigma2_P)

    # TWO separate weak-fluctuation tests (see olb.turbulence.andrews.scintillation
    # and Conflict C-05 / TL-05 in docs/andrews-crosscheck.md):
    #
    #  1. The REGIME gate: is the analytic Rytov index valid at all? This is a
    #     turbulence-strength test on the Rytov variance sigma2_R, and for a
    #     Gaussian beam it needs BOTH sigma_R^2 < 1 AND sigma_R^2 Lambda^(5/6) < 1
    #     (Andrews and Phillips 2005, Ch. 5, Eq. (16), printed p. 140). rytov_weak
    #     reads Lambda (the receiver-plane beam parameter, collimated launch) so a
    #     focused beam trips the gate a plane-wave test would pass. It returns
    #     "weak" (firm), "soft" (canonical weak, a soft warning) or "hard"
    #     (leaving weak, a hard warning).
    #  2. The lognormal-PDF house rule: is the fade PDF SHAPE trusted? This is a
    #     tighter test on the index sigma2_I < LOGNORMAL_PDF_LIMIT (0.25, from the
    #     optimistic lognormal tail, Ch. 11.3, printed p. 451), kept SEPARATE.
    # The plane-wave Rytov variance sets the regime gate. Use the canonical
    # andrews form (this is outside the trace above, so it adds no provenance).
    sigma2_R = float(rytov_variance(wavelength, L, cn2, wave='plane'))
    regime = rytov_weak(sigma2_R, Lambda)
    pdf_valid = bool(sigma2_I < LOGNORMAL_PDF_LIMIT)
    # The traced physics functions own the beam type, the regime, the spectrum,
    # and the circular-aperture / no-obscuration constraints; the merge inherits
    # their union and the traced Dios reliability violation. State the three
    # headline fields explicitly (this is a Gaussian-beam lognormal Term).
    assumptions = trace.merge(
        beam_type=BEAM_GAUSSIAN,
        turbulence_regime=REGIME_WEAK,
        spectrum=SPECTRUM_KOLMOGOROV,
        validity="Weak fluctuation, TWO conditions. Regime: the Gaussian-beam "
                 "gate sigma_R^2 < 1 AND sigma_R^2 Lambda^(5/6) < 1 (Andrews and "
                 "Phillips 2005, Ch. 5, Eq. (16)). PDF shape: sigma2_I < 0.25, a "
                 "house rule (Ch. 11.3). The point index is the on-axis Gaussian "
                 "beam-wave index over the constant-Cn2 horizontal path (Dios et "
                 "al. 2004, Eq. 16). The weak aperture-averaging factor assumes a "
                 "uniform circular aperture with no central obscuration and a "
                 "small inner scale (Andrews and Phillips 2005, Ch. 10).",
    )
    # The weak Kolmogorov averaging factor models a uniform circular aperture. A
    # central obscuration (Cassegrain secondary) breaks that. It is a
    # scenario-level fact the physics never sees, so flag it at the factory level.
    if rx.obscuration_ratio > 0.0:
        assumptions.flag(
            f"The far aperture has a central obscuration "
            f"(ratio={rx.obscuration_ratio:.3f}); the weak aperture-averaging "
            "factor assumes a uniform circular aperture and does not model it.",
            source="factory:links.terrestrial",
        )
    if regime == 'soft':
        warnings.warn(
            f"horizontal Gaussian-beam Rytov variance sigma2_R={sigma2_R:.3f} "
            f"(Lambda={Lambda:.3f}) is past the confident-weak value "
            f"{RYTOV_CONFIDENT_WEAK} but within the book weak limit "
            f"{WEAK_REGIME_LIMIT}; the analytic index is usable, but check it "
            "against a Monte Carlo near the top of the band."
        )
    elif regime == 'hard':
        # The hard-tier regime violation is now OWNED by the physics: the traced
        # on_axis_scintillation_index check flags the Dios beam-wave index as
        # unreliable (sigma_chi^2 >= 0.6). So the factory keeps only the warning
        # here (verbatim); the assumptions violation arrives through the trace.
        warnings.warn(
            f"horizontal Gaussian-beam Rytov variance sigma2_R={sigma2_R:.3f} "
            f"(Lambda={Lambda:.3f}) leaves the weak regime -- the analytic index "
            "is not trusted. Use the fidelity-2 Monte Carlo."
        )
    if not pdf_valid:
        # The lognormal-PDF house rule (0.25 on sigma2_I) is a PDF-SHAPE decision
        # that the physics does not gate (the traced regime check is the wider
        # Dios reliability bound). So it stays a factory flag, source-tagged.
        assumptions.flag(
            f"sigma2_I={sigma2_I:.3f} exceeds the lognormal-PDF house limit "
            f"{LOGNORMAL_PDF_LIMIT}; the lognormal fade tail is not trusted. Use "
            "gamma-gamma or Monte Carlo.",
            source="factory:links.terrestrial",
        )
        warnings.warn(
            f"horizontal Gaussian-beam scintillation index sigma2_I={sigma2_I:.3f} "
            f">= {LOGNORMAL_PDF_LIMIT} -- the lognormal fade tail is optimistic. "
            "Use a gamma-gamma or a Monte Carlo model."
        )

    return irradiance_fade_term(
        "scintillation", "turbulence",
        mean_log=lognormal_mean_log(sigma_l2),
        quantile=lambda p: lognormal_quantile(p, sigma_l2),
        rvs=lambda n, rng: lognormal_rvs(n, sigma_l2, rng),
        note="horizontal Gaussian-beam lognormal scintillation, aperture-averaged",
        meta={
            "model": "lognormal",
            "beam": "gaussian-horizontal",
            "sigma2_I": sigma2_I,
            "sigma2_P": float(sigma2_P),
            "aperture_averaging_factor": A,
            "sigma2_R": sigma2_R,
            "Lambda": Lambda,
            "rytov_regime": regime,
            # weak_fluctuation_valid keeps its meaning: the lognormal PDF SHAPE
            # is trusted (sigma2_I < LOGNORMAL_PDF_LIMIT). The regime label above
            # is the separate index-validity test.
            "weak_fluctuation_valid": pdf_valid,
            "weak_fluctuation_limit": LOGNORMAL_PDF_LIMIT,
        },
        assumptions=assumptions,
    )


def _terrestrial_fidelity2_terms(scenario, geometry, wave, turbulence=True):
    '''
    The fidelity-2 wave-optics Terms of a terrestrial link.

    A terrestrial link is fully simulated end to end on ONE flat grid, so the
    FULL launch-to-detector loss splits cleanly into:
      - a DETERMINISTIC vacuum-optics Term (the no-turbulence loss: launch
        truncation + geometric spread + aperture capture + vacuum fibre coupling);
      - a STOCHASTIC turbulence Term (the fade).
    Their sum reconstructs the direct launch-to-detector turbulent loss exactly
    (the vacuum baselines cancel; see olb.models.waveoptics). Together
    they replace the analytic geometric, launch-truncation, scintillation, and
    coupling Terms. `wave` is a Fidelity2Bundle from
    olb.models.waveoptics.run_fidelity2.

    An SMF receiver gets the composite fibre penalty (aperture capture x fibre
    coupling). An MMF receiver gets the aperture-power penalty PLUS the
    light-bucket core-coupling Term (waveoptics_mmf_coupling_term). That coupling
    is the ABSOLUTE core capture (relative to the COLLECTED power), so it does not
    double-count the aperture capture, and it already holds the detector defocus,
    because the run folds it into mmf_eta (see olb.waveoptics.mmf and
    olb.waveoptics.turbulence.run). An Aperture (bucket) receiver gets the
    aperture-power penalty only.

    With `turbulence` False, or with a VACUUM-ONLY bundle (turbulent None), the
    Term set is DETERMINISTIC: the vacuum-optics Term alone, plus the vacuum MMF
    core-capture Term for an MMF receiver. No stochastic Term is built.
    '''
    from ..models.waveoptics import (waveoptics_vacuum_term,
                                     waveoptics_turbulence_term,
                                     waveoptics_mmf_coupling_term,
                                     waveoptics_vacuum_mmf_term)
    from ..waveoptics.field import Power
    from ..terminal import SMF, MMF
    tx = scenario.tx_terminal
    rx = scenario.rx_terminal
    is_smf = isinstance(rx.detector, SMF)
    is_mmf = isinstance(rx.detector, MMF)

    if wave.turbulent is None and turbulence:
        raise ValueError(
            "the `wave` bundle is vacuum-only (turbulent is None), but the "
            "budget asks for turbulence. Run "
            "olb.models.waveoptics.run_fidelity2 WITHOUT turbulence=False, or "
            "pass turbulence=False to the budget."
        )

    vac = waveoptics_vacuum_term(wave.vacuum, include_smf=is_smf)
    if not turbulence or wave.turbulent is None:
        # VACUUM-ONLY: the deterministic Terms alone. The SMF coupling is
        # already inside the vacuum Term (include_smf); the MMF light bucket
        # needs its own deterministic core-capture Term, computed on the
        # receive-clipped vacuum field.
        if is_mmf:
            return [vac, waveoptics_vacuum_mmf_term(wave.vacuum, rx.detector,
                                                    rx.aperture_m)]
        return [vac]

    trials = wave.turbulent.trials
    coll = np.array([t.collected_power for t in trials], dtype=float)
    # The vacuum aperture fraction on the SAME grid: collected / after-tx-clip,
    # matching the terrestrial collected_power normalisation.
    vac_coll = float(Power(wave.vacuum.stages[3][1])
                     / Power(wave.vacuum.stages[1][1]))
    if is_smf:
        eta = np.array([t.smf_eta for t in trials], dtype=float)
        vac_smf_eta = 10.0 ** (-wave.vacuum.smf_coupling_db / 10.0)
        loss_db = (-10.0 * np.log10(coll / vac_coll)
                   - 10.0 * np.log10(eta / vac_smf_eta))
        note = ("terrestrial turbulence penalty (wave optics): aperture-power and "
                "fibre-coupling loss relative to the vacuum baseline.")
    else:
        loss_db = -10.0 * np.log10(coll / vac_coll)
        note = ("terrestrial turbulence penalty (wave optics): aperture-power "
                "loss relative to the vacuum baseline.")

    L = float(scenario.channel.path_length_m)
    hs = np.linspace(0.0, L, _SCINT_GRID_N)
    cn2_profile = np.full_like(hs, float(scenario.channel.cn2))
    sigma2_I = float(on_axis_scintillation_index(
        hs, cn2_profile, tx.transmitter.waist_m, rx.wavelength_m,
        elevation_deg=90.0, path_length_m=None))
    pen = waveoptics_turbulence_term(
        wave.turbulent, loss_db=loss_db, beam_type=BEAM_GAUSSIAN,
        sigma2_I=sigma2_I, note=note)
    terms = [vac, pen]
    if is_mmf:
        # The light-bucket core coupling (absolute, with fade). It is the fraction
        # of the COLLECTED power that enters the core, so it composes with the
        # vacuum and aperture-penalty Terms (launch -> collected) without a
        # double-count. It already carries the detector defocus (spot growth),
        # because the run folds it into mmf_eta.
        mmf_term = waveoptics_mmf_coupling_term(
            wave.turbulent, beam_type=BEAM_GAUSSIAN, sigma2_I=sigma2_I,
            note="terrestrial MMF light-bucket coupling (wave optics): absolute "
                 "core capture relative to the collected power, with the detector "
                 "defocus.")
        terms.append(mmf_term)
    return terms


def terrestrial_budget(scenario, geometry, *, fidelity=0, scintillation=True,
                       turbulence=True, wave=None):
    '''
    Assemble the terrestrial budget at a chosen fidelity.

    `fidelity` is a WHOLE-PATH choice (see the README fidelity ladder):

      - fidelity=0 (the default, analytic). The deterministic Terms (geometric
        spreading, horizontal extinction, pointing jitter) are exact. The
        receive-side turbulence effect depends on the front end:
          * an Aperture (bucket) or no detector gets the horizontal Gaussian-beam
            scintillation Term (terrestrial_scintillation_term, a real analytic
            fade);
          * an SMF detector gets the mean-only fibre-coupling Term
            (terrestrial_smf_coupling_term) plus the tip-tilt walk-off fade when
            the coupling optics are set. The mean-only Term locks the budget to
            fidelity 0 (it then refuses a fade margin);
          * an MMF detector gets the multimode spot-in-core coupling plus the
            walk-off fade (a real fade).
        The coupling / scintillation Term REPLACES the standalone scintillation
        for a fibre receiver (no double-count).
      - fidelity=1 is UNAVAILABLE for a terrestrial link and raises. FAST is a
        far-field plane-wave-source model; a near-field finite Gaussian beam
        needs the split-step model of fidelity 2 (see backlog 1-1).
      - fidelity=2 (wave optics). The whole path is a field simulation. It gives
        TWO Terms: a deterministic vacuum-optics Term (launch truncation +
        geometric spread + aperture capture + vacuum fibre coupling) and a
        stochastic turbulence Term (the fade). Together they REPLACE the
        geometric, launch-truncation, scintillation, and coupling Terms. Only the
        analytic extinction (molecular absorption, never in the field sim) and
        pointing (mechanical jitter) Terms stay. It needs a precomputed `wave`
        bundle (olb.models.waveoptics.run_fidelity2); the budget never runs the
        split-step propagation itself.

    Set scintillation=False to drop the scintillation Term at fidelity 0 and keep
    only the deterministic Terms (for example to sweep an array path length, where
    the scalar-only scintillation Term does not broadcast; loop per distance).

    Parameters:
        scenario : TerrestrialScenario
            A terrestrial link case. The `direction` maps the roles (forward:
            tx = near, rx = far; reverse swaps them). Its TerrestrialChannel
            carries path_length_m, attenuation_db_per_km, cn2.
        geometry : HorizontalPath
            The horizontal path (reads slant_range_m = path length).
        fidelity : int
            0 (analytic, the default), 1 (unavailable, raises), or 2 (wave optics,
            needs `wave`).
        scintillation : bool
            Add the fidelity-0 scintillation Term for an aperture / no-detector
            receiver when True (the default).
        turbulence : bool
            Master turbulence switch, at fidelity 0 AND fidelity 2. When False,
            drop EVERY turbulence quantity. At fidelity 0: no scintillation
            Term, and the fibre-coupling Terms keep only their static parts. At
            fidelity 2: no wave-optics turbulence Term and no stochastic
            coupling Term, so the budget shows the deterministic vacuum-optics
            Term alone (plus the vacuum MMF core capture for an MMF receiver).
            The deterministic Terms (geometric, extinction, launch truncation)
            and the transmit pointing jitter stay at both rungs. At fidelity 2
            pair it with olb.models.waveoptics.run_fidelity2(turbulence=False),
            which makes no screens and no trials.
        wave : Fidelity2Bundle, list, or Campaign, optional
            The precomputed wave-optics record for fidelity=2: a
            Fidelity2Bundle, a list of them, or a Campaign. Run it with
            olb.models.waveoptics.run_fidelity2, or store it with
            olb.waveoptics.turbulence.Campaign and pass the campaign itself
            (olb.models.waveoptics.resolve_wave turns it into the bundle).

    Returns:
        Budget
            The budget with the scenario set.

    Raises:
        ValueError
            If fidelity is not 0/1/2, if fidelity=1 (unavailable for terrestrial),
            if fidelity=2 without a `wave` bundle, or if fidelity=2 with
            turbulence=True and a VACUUM-ONLY bundle.
    '''
    if fidelity not in (0, 1, 2):
        raise ValueError(f"fidelity must be 0, 1, or 2, got {fidelity!r}.")
    if fidelity == 1:
        raise ValueError(
            "fidelity=1 is unavailable for a terrestrial link. FAST is a "
            "far-field plane-wave-source model; a near-field finite Gaussian "
            "beam needs the split-step model of fidelity 2. Use fidelity=0 "
            "(analytic) or fidelity=2 (wave optics)."
        )
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
        # The two wave-optics Terms replace geometric, launch truncation,
        # scintillation, and coupling. Only extinction (absorption) and pointing
        # (mechanical jitter) stay analytic.
        terms = [
            terrestrial_extinction_term(scenario, geometry),
            pointing_loss_term(scenario, geometry),
        ]
        terms += _terrestrial_fidelity2_terms(scenario, geometry, wave,
                                              turbulence=turbulence)
        return Budget(terms, scenario=scenario)

    # fidelity 0: the analytic budget.
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
                                       terrestrial_smf_walkoff_term)
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
            terms.append(terrestrial_smf_walkoff_term(scenario, geometry,
                                          turbulence=turbulence))
    elif isinstance(rx.detector, MMF):
        # An MMF (light bucket) replaces the scintillation Term with the geometric
        # spot-in-core coupling plus the tip-tilt walk-off fade (no double-count).
        from ..models.coupling import terrestrial_mmf_coupling_term
        terms.append(terrestrial_mmf_coupling_term(scenario, geometry, turbulence=turbulence))
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
    assert q99_scint is not None and np.isfinite(q99_scint) and q99_scint > scint.mean_db, \
        (q99_scint, scint.mean_db)

    # --- parity with the retired inline lognormal faces (backlog I-2) --------
    # The three dB faces now come from olb.models.fade.irradiance_fade_term, not
    # the old inline formula. Rebuild the retired faces and assert a match, so a
    # future change to the shared adapter cannot silently move the numbers. The
    # mean and the sampler are BYTE identical; the quantile matches to machine
    # precision (the adapter takes -10 log10(exp(x)), the retired code took
    # -10 x / ln10 directly).
    from scipy.stats import norm as _norm
    _ln10 = np.log(10.0)
    _sl2 = np.log(1.0 + scint.meta["sigma2_P"])
    _sl = np.sqrt(_sl2)
    assert scint.mean_db == (5.0 / _ln10) * _sl2, scint.mean_db
    for _p in (0.01, 0.5, 0.99):
        _old_q = -10.0 / _ln10 * (-_sl2 / 2.0 + _sl * _norm.ppf(1.0 - _p))
        assert abs(scint.quantile_db(_p) - _old_q) < 1e-12, \
            (_p, scint.quantile_db(_p), _old_q)
    _a = scint.sample_db(20_000, np.random.default_rng(4321))
    _b = -10.0 * np.log10(np.random.default_rng(4321).lognormal(
        mean=-_sl2 / 2.0, sigma=_sl, size=20_000))
    assert np.max(np.abs(_a - _b)) == 0.0
    print(f"[parity] terrestrial lognormal faces match retired inline "
          f"(mean {scint.mean_db:.6f} dB, byte-identical sampler)")

    # The aperture-averaging win: a larger receive aperture shrinks the flux index
    # and the fade. Sweep D and check both fall monotonically.
    D_sweep = [0.05, 0.1, 0.2, 0.4, 0.8]
    sig_P, fades = [], []
    for D in D_sweep:
        s = terrestrial_scintillation_term(
            _terr(0.02, 5e3, far_aperture=D), geom)
        sig_P.append(s.meta["sigma2_P"])
        q = s.quantile_db(0.99)
        fades.append(float(q) if q is not None else float("nan"))
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

    # --- WP3c: the scintillation Term inherits traced provenance --------------
    # The weak default Term names its physics sources.
    prov = scint.assumptions.provenance
    assert prov, "the scintillation Term must carry traced provenance"
    assert any("on_axis_scintillation_index" in s for s in prov), prov
    assert any("aperture_averaging_factor_weak" in s for s in prov), prov
    # The hard-tier regime flag MIGRATED to the physics: the strong and long-path
    # cases now read not-ok through the TRACED Dios reliability check, not a
    # hand-built factory flag. Assert membership, not the exact count.
    assert any("on_axis_scintillation_index" in v and "unreliable" in v
               for v in strong.assumptions.violations), strong.assumptions.violations
    assert any("on_axis_scintillation_index" in v and "unreliable" in v
               for v in long_path.assumptions.violations), long_path.assumptions.violations
    # Budget.check() on the default aperture budget reports NO untraced-guard
    # entry (the turbulence-category scintillation Term now carries provenance).
    guard = [(n, r) for n, r in budget.check(warn=False)
             if "did not open the assumption collection context" in r]
    assert guard == [], guard

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

    # --- fidelity=2 whole-path wave optics -----------------------------------
    def _smf_scn():
        return TerrestrialScenario(
            near=Terminal(aperture_m=0.3, wavelength_m=1550e-9,
                          transmitter=Transmitter(waist_m=0.02, power_dbm=30)),
            far=Terminal(aperture_m=0.2, wavelength_m=1550e-9,
                         detector=SMF(sensitivity_dbm=-40)),
            channel=TerrestrialChannel(path_length_m=3e3, attenuation_db_per_km=0.5,
                                       cn2=1e-14))

    wo_scn = _smf_scn()
    # Guards need NO run: bad fidelity, fidelity=1 (unavailable), fidelity=2 with
    # no bundle.
    for bad in (3, -1, "mean"):
        try:
            terrestrial_budget(wo_scn, HorizontalPath(3e3), fidelity=bad)
        except ValueError as e:
            assert "fidelity must be 0, 1, or 2" in str(e)
        else:
            raise AssertionError(f"fidelity={bad!r} must raise")
    try:
        terrestrial_budget(wo_scn, HorizontalPath(3e3), fidelity=1)
    except ValueError as e:
        assert "unavailable for a terrestrial link" in str(e)
    else:
        raise AssertionError("fidelity=1 must raise for terrestrial")
    try:
        terrestrial_budget(wo_scn, HorizontalPath(3e3), fidelity=2)
    except ValueError as e:
        assert "needs a precomputed `wave` bundle" in str(e)
    else:
        raise AssertionError("fidelity=2 without a bundle must raise")
    # The default is UNCHANGED: fidelity=0 SMF stays mean-only (fade locked).
    with _warnings.catch_warnings():
        _warnings.simplefilter("ignore")
        assert not terrestrial_budget(wo_scn, HorizontalPath(3e3)).provides_fade

    # The real fidelity-2 build needs one run_fidelity2 (skip if aotools absent).
    from ..models.waveoptics import run_fidelity2
    try:
        with _warnings.catch_warnings():
            _warnings.simplefilter("ignore")
            bundle = run_fidelity2(wo_scn, HorizontalPath(3e3), preset="rapid",
                                   n_trials=16, seed=3)
            f2 = terrestrial_budget(wo_scn, HorizontalPath(3e3), fidelity=2,
                                    wave=bundle)
    except ImportError:
        print("aotools not installed; skipping the terrestrial fidelity-2 run.")
        f2 = None
    if f2 is not None:
        names = [t.name for t in f2.terms]
        cats = [t.category for t in f2.terms]
        # Two wave-optics Terms; no analytic geometric/scintillation/coupling.
        vac = next(t for t in f2.terms if t.meta.get("model") == "waveoptics-vacuum")
        turb = next(t for t in f2.terms if t.meta.get("model") == "waveoptics")
        assert not vac.stochastic and turb.stochastic and not turb.mean_only
        assert "geometric spreading" not in names and "scintillation" not in names
        assert "receive coupling (SMF)" not in names
        assert "SMF tip-tilt walk-off" not in names
        # Extinction and pointing stay analytic.
        assert "atmospheric" in cats and "pointing" in cats
        # The fidelity-0 lock is gone: a real fade margin.
        assert f2.provides_fade and np.isfinite(f2.fade_margin_db(0.9))
        print(f"terrestrial fidelity 2 (3 km, rapid, 16 trials): vacuum "
              f"{vac.mean_db:.2f} dB + turbulence {turb.mean_db:.2f} dB, "
              f"total {f2.total_loss_db():.2f} dB")

    # --- fidelity-2 master turbulence switch ---------------------------------
    # turbulence=False at fidelity 2 gives a VACUUM-ONLY Term set: the
    # deterministic vacuum-optics Term plus extinction and pointing. The bundle
    # comes from run_fidelity2(turbulence=False), which makes no trials.
    try:
        with _warnings.catch_warnings():
            _warnings.simplefilter("ignore")
            vac_bundle = run_fidelity2(wo_scn, HorizontalPath(3e3),
                                       preset="rapid", turbulence=False,
                                       progress=False)
            f2_off = terrestrial_budget(wo_scn, HorizontalPath(3e3), fidelity=2,
                                        wave=vac_bundle, turbulence=False)
    except ImportError:
        f2_off = None
    if f2_off is not None:
        assert vac_bundle.turbulent is None
        off_cats = [t.category for t in f2_off.terms]
        assert "turbulence" not in off_cats, off_cats
        assert "coupling" not in off_cats, off_cats
        # The deterministic backbone stays: the vacuum Term, extinction, pointing.
        assert any(t.meta.get("model") == "waveoptics-vacuum" for t in f2_off.terms)
        assert "atmospheric" in off_cats and "pointing" in off_cats
        # A vacuum-only bundle with turbulence=True raises a helpful error.
        try:
            terrestrial_budget(wo_scn, HorizontalPath(3e3), fidelity=2,
                               wave=vac_bundle)
        except ValueError as e:
            assert "vacuum-only" in str(e), str(e)
        else:
            raise AssertionError("a vacuum-only bundle with turbulence must raise")
        # turbulence=False with the FULL bundle is allowed: it drops the
        # stochastic Terms, so the Term set matches the vacuum-only bundle.
        if f2 is not None:
            full_off = terrestrial_budget(wo_scn, HorizontalPath(3e3),
                                          fidelity=2, wave=bundle,
                                          turbulence=False)
            assert ([t.name for t in full_off.terms]
                    == [t.name for t in f2_off.terms])
        print(f"terrestrial fidelity 2, turbulence=False (3 km): "
              f"{f2_off.total_loss_db():.2f} dB, "
              f"terms {[t.name for t in f2_off.terms]}")

    # --- fidelity-2 MMF light bucket, and the non-focal-plane detector -------
    # The MMF core coupling is now routed into the fidelity-2 budget. The
    # detector defocus grows the spot, which folds into the per-trial mmf_eta in
    # the run, so the coupling loss grows off the true focus.
    def _f2_mmf(defocus_m=0.0):
        s = TerrestrialScenario(
            near=Terminal(aperture_m=0.3, wavelength_m=1550e-9,
                          transmitter=Transmitter(waist_m=0.02, power_dbm=30)),
            far=Terminal(aperture_m=0.2, wavelength_m=1550e-9,
                         pointing_jitter_rad=5e-6,
                         detector=MMF(core_radius_m=25e-6, optimal_focus=True,
                                      defocus_m=defocus_m, sensitivity_dbm=-38)),
            channel=TerrestrialChannel(path_length_m=3e3, attenuation_db_per_km=0.5,
                                       cn2=1e-14))
        with _warnings.catch_warnings():
            _warnings.simplefilter("ignore")
            wave = run_fidelity2(s, HorizontalPath(3e3), preset="rapid",
                                 n_trials=24, seed=3)
            bud = terrestrial_budget(s, HorizontalPath(3e3), fidelity=2, wave=wave)
        mmf = next(t for t in bud.terms if t.name == "receive coupling (MMF)")
        return bud, mmf

    try:
        b_focus, m_focus = _f2_mmf(defocus_m=0.0)
        # A NEGATIVE defocus moves the detector AWAY from the true focus. The
        # received beam is a diverging Gaussian, so its true focus sits BEYOND the
        # focal plane (at +dz_curv, some 6 mm here). So a positive defocus of a few
        # mm moves TOWARD the true focus and wins power back. See
        # olb.models.coupling.terrestrial and validation/defocus.
        b_defoc, m_defoc = _f2_mmf(defocus_m=-2e-3)
    except ImportError:
        b_focus = None
    if b_focus is not None:
        # The MMF coupling Term is now in the fidelity-2 budget, with a real fade.
        assert m_focus.category == "coupling" and m_focus.stochastic
        assert b_focus.provides_fade
        # A defocus away from the true focus grows the spot, so the core captures
        # less (more loss).
        assert m_defoc.mean_db > m_focus.mean_db, (m_defoc.mean_db, m_focus.mean_db)
        print(f"terrestrial fidelity 2 MMF (25 um core): focus {m_focus.mean_db:.2f} "
              f"dB -> defocus -2 mm {m_defoc.mean_db:.2f} dB")

        # An MMF receiver with turbulence=False keeps ONE deterministic
        # core-capture Term, computed on the receive-clipped vacuum field.
        s_mmf = TerrestrialScenario(
            near=Terminal(aperture_m=0.3, wavelength_m=1550e-9,
                          transmitter=Transmitter(waist_m=0.02, power_dbm=30)),
            far=Terminal(aperture_m=0.2, wavelength_m=1550e-9,
                         pointing_jitter_rad=5e-6,
                         detector=MMF(core_radius_m=25e-6, optimal_focus=True,
                                      sensitivity_dbm=-38)),
            channel=TerrestrialChannel(path_length_m=3e3,
                                       attenuation_db_per_km=0.5, cn2=1e-14))
        with _warnings.catch_warnings():
            _warnings.simplefilter("ignore")
            mmf_vac_bundle = run_fidelity2(s_mmf, HorizontalPath(3e3),
                                           preset="rapid", turbulence=False,
                                           progress=False)
            mmf_off = terrestrial_budget(s_mmf, HorizontalPath(3e3), fidelity=2,
                                         wave=mmf_vac_bundle, turbulence=False)
        mmf_vac_term = next(t for t in mmf_off.terms if t.category == "coupling")
        assert not mmf_vac_term.stochastic and not mmf_vac_term.mean_only
        assert mmf_vac_term.meta["model"] == "waveoptics-vacuum"
        assert not any(t.category == "turbulence" for t in mmf_off.terms)
        # The vacuum core capture is BETTER than the turbulent mean (no fade).
        assert mmf_vac_term.mean_db < m_focus.mean_db, (mmf_vac_term.mean_db,
                                                        m_focus.mean_db)
        print(f"terrestrial fidelity 2 MMF, turbulence=False: core capture "
              f"{mmf_vac_term.mean_db:.2f} dB (turbulent mean "
              f"{m_focus.mean_db:.2f} dB)")

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
    wo_off_q99 = wo_off.quantile_db(0.99)
    assert smf_off.provides_fade and wo_off_q99 is not None and wo_off_q99 > wo_off.mean_db
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
