'''
Receive-coupling Term for a downlink receive terminal.

This module builds the ONE receive-coupling Term of a downlink receive terminal.
The Compensation stack and the Detector are one physical chain: the residual
wavefront that the Compensation leaves sets the coupling into the Detector. So
this module emits a single Term, not two. The shared single-mode-fibre physics
lives in olb.models.coupling._common; the phase physics lives in
olb.turbulence.ao.

Two detector front ends:

  Aperture (bucket): a power-in-bucket detector integrates the intensity. It is
  phase-insensitive, so the Compensation stack does not change its coupling. The
  turbulence penalty is aperture-averaged scintillation. This term reuses the
  downlink lognormal scintillation Term unchanged, so it reproduces the plain
  downlink behaviour.

  SMF (single-mode fibre): the fibre couples only the field that matches the
  fibre mode. Two fidelities:
    smf_fidelity="fast" (the default, and the ONLY statistical model): the true
        LP01 modal overlap from FAST. It gives the mean, the quantile, and the
        fade. See olb.models.fast.
    smf_fidelity="mean": a cheap analytic MEAN-ONLY estimate, for when only the
        expected coupling loss is wanted (no fade). The coupling efficiency eta
        falls with the residual phase variance sigma^2_res that the Compensation
        stack leaves. Two limits of one overlap physics (see _common):
          small residual: extended Marechal, eta = eta_max * exp(-sigma^2_res).
          large residual: the Dikmelik-Davidson uncorrected coupling curve.
        The mean loss is -10*log10(eta). The Term is DETERMINISTIC: it carries no
        fade (no sampler, no quantile). For the fade, use smf_fidelity="fast".

How the SMF Term avoids double-counting the geometric aperture capture:
  The geometric spreading Term already carries the free-space spread and the
  receive-aperture power-in-bucket capture. The SMF eta is a MULTIPLICATIVE
  fibre-coupling efficiency on the aperture-collected field. So the SMF Term adds
  only the coupling loss -10*log10(eta); it does NOT re-add the aperture capture.
  In the budget the receive-coupling Term REPLACES the standalone scintillation
  Term (the turbulence effect is now the coupling loss and fade), and the
  geometric spreading Term stays. See olb.links.downlink.downlink_budget.

Sources:
  Noll residual variance: R. J. Noll, JOSA 66(3), 207 (1976). See
  olb.turbulence.ao.
  SMF coupling / Marechal: extended Marechal approximation.
  Uncorrected SMF coupling against D/r0: Y. Dikmelik and F. M. Davidson, Appl.
  Opt. 44(23), 4946-4952 (2005), DOI 10.1364/AO.44.004946.
'''

import numpy as np

from ...results import Term
from ...assumptions import (Assumptions, BEAM_PLANE_WAVE, BEAM_GAUSSIAN,
                            REGIME_WEAK, REGIME_NA, SPECTRUM_KOLMOGOROV,
                            SPECTRUM_NA)
from ...terminal import Aperture, SMF
from ...turbulence.profiles import DEFAULT_HS, default_cn2_profile
from ...turbulence.ao import (plane_wave_fried_parameter_profile,
                             apply_compensation)
from ...links.downlink import downlink_scintillation_term
from ._common import (_smf_eta_max, _smf_coupling_efficiency, _effective_dr0,
                     _smf_static_term, SMF_DEEP_TURBULENCE_DR0)


def _aperture_term(scenario, geometry, *, hs, cn2_profile):
    '''
    Build the receive-coupling Term for an aperture (bucket) detector.

    Reuse the downlink lognormal aperture-averaged scintillation Term unchanged.
    A bucket detector is phase-insensitive, so the compensation stack does not
    change it. Rename the Term to the coupling category, so the total is parity
    with the plain downlink scintillation.
    '''
    term = downlink_scintillation_term(scenario, geometry, model="lognormal",
                                       aperture_average=True, hs=hs,
                                       cn2_profile=cn2_profile)
    term.name = "receive coupling (aperture)"
    term.category = "coupling"
    term.note = "bucket detector: aperture-averaged lognormal scintillation"
    return term


def _smf_mean_term(scenario, geometry, *, hs, cn2_profile, turbulence=True):
    '''
    Build the MEAN-ONLY receive-coupling Term for a single-mode-fibre detector.

    Get r0 from the plane-wave Fried parameter. Get the residual phase variance
    from the compensation stack. Convert it to the coupling efficiency eta, then
    to the mean coupling loss -10*log10(eta). The Term is DETERMINISTIC: it has no
    sampler and no quantile, so it carries no fade. For the fade (and the true
    modal overlap) use smf_fidelity="fast" (olb.models.fast).

    With turbulence=False the residual wavefront error drops to zero, so the Term
    is the static mode-match loss only (_smf_static_term).
    '''
    terminal = scenario.rx_terminal
    detector = terminal.detector
    wavelength = terminal.wavelength_m
    D = terminal.aperture_m
    eta_max = _smf_eta_max(detector, D, wavelength)

    if not turbulence:
        return _smf_static_term(eta_max)

    elev = geometry.elevation_deg
    r0 = plane_wave_fried_parameter_profile(cn2_profile, hs, wavelength, elev)
    residual = apply_compensation(terminal.compensation, D, r0)
    sigma2_res = residual.variance
    eta = _smf_coupling_efficiency(sigma2_res, eta_max)
    coupling_loss = -10.0 * np.log10(eta)     # positive dB, scalar or per-elevation
    dr0_eff = _effective_dr0(sigma2_res)
    base_shape = np.shape(coupling_loss)

    assumptions = Assumptions(
        beam_type=BEAM_PLANE_WAVE,
        turbulence_regime=REGIME_WEAK,
        spectrum=SPECTRUM_KOLMOGOROV,
        validity="MEAN-ONLY single-mode-fibre coupling loss. Extended Marechal for "
                 "a small residual; Dikmelik-Davidson uncorrected coupling for a "
                 "large residual. The Dikmelik-Davidson coupling assumes a uniform "
                 "circular aperture with no central obscuration. This Term is "
                 "DETERMINISTIC: it gives the expected coupling loss only and models "
                 "NO fade. For the fade and the true modal overlap use "
                 "smf_fidelity='fast'.",
    )
    # MEAN-ONLY. The Term carries no coupling fade. Always flag it, so a fade
    # margin (for example a 99% link margin) is never read off this Term.
    assumptions.flag(
        "Mean-only SMF coupling: this Term is the expected coupling loss and models "
        "NO fade (no sampler, no quantile). Every fade margin read from it is wrong. "
        "Use smf_fidelity='fast' for the statistical (fidelity-1) coupling."
    )
    # Dikmelik-Davidson assumes a uniform circular aperture. A central obscuration
    # on the receive aperture breaks the coupling curve. Flag the violation.
    if terminal.obscuration_ratio > 0.0:
        assumptions.flag(
            f"The receive aperture has a central obscuration "
            f"(ratio={terminal.obscuration_ratio:.3f}); the Dikmelik-Davidson "
            "coupling curve assumes a uniform circular aperture and does not "
            "model it."
        )
    # Flag deep turbulence, where the practical coupling curve is extrapolated.
    dr0_max = float(np.max(dr0_eff))
    if dr0_max > SMF_DEEP_TURBULENCE_DR0:
        assumptions.flag(
            f"effective D/r0={dr0_max:.1f} exceeds {SMF_DEEP_TURBULENCE_DR0:.0f}; "
            "the practical uncorrected coupling curve is extrapolated. Use "
            "smf_fidelity='fast'."
        )

    note = (f"SMF coupling (mean-only), eta_max={eta_max:g}, "
            f"n_comp_modes={residual.n_modes}")
    return Term(
        name="receive coupling (SMF)",
        category="coupling",
        mean_db=float(coupling_loss) if base_shape == () else coupling_loss,
        sampler=None,       # deterministic: mean-only, no fade
        quantile=None,
        note=note,
        meta={
            "detector": "SMF",
            "model": "mean-only",
            "eta": float(eta) if base_shape == () else np.asarray(eta),
            "coupling_loss_db": float(coupling_loss) if base_shape == () else np.asarray(coupling_loss),
            "sigma2_res": float(sigma2_res) if np.ndim(sigma2_res) == 0 else np.asarray(sigma2_res),
            "effective_D_over_r0": float(dr0_eff) if np.ndim(dr0_eff) == 0 else np.asarray(dr0_eff),
            "r0_m": float(r0) if np.ndim(r0) == 0 else np.asarray(r0),
            "n_comp_modes": residual.n_modes,
        },
        assumptions=assumptions,
        mean_only=True,     # fidelity-0: expected coupling loss, no fade
    )


def downlink_coupling_term(scenario, geometry, *, hs=None, cn2_profile=None,
                           n_samples=2000, smf_fidelity="fast", fast_params=None,
                           turbulence=True):
    '''
    Build the ONE receive-coupling Term of a downlink receive terminal.

    Read scenario.rx_terminal. Dispatch on the detector type. An Aperture detector
    reuses the downlink aperture-averaged scintillation, so it is parity with the
    plain downlink. An SMF detector picks the fidelity with smf_fidelity:

    - "fast" (the default): the fidelity-1 true LP01 modal overlap from FAST. It
      gives the mean, the quantile, and the fade. Needs fast-aosim. See
      olb.models.fast.
    - "mean": a cheap analytic MEAN-ONLY estimate (extended-Marechal / Dikmelik-
      Davidson from the residual wavefront). The Term is DETERMINISTIC and models
      no fade. Use it when only the expected coupling loss is wanted. See
      _smf_mean_term.

    Parameters:
        scenario : SpaceScenario or TerrestrialScenario
            Reads rx_terminal, link.wavelength_m, and the site Cn2 profile.
        geometry : CircularOrbit or TLEPass
            Reads elevation_deg. A scalar elevation gives a scalar Term.
        hs : numpy.ndarray, optional
            Heights above the ground station [m]. Defaults to DEFAULT_HS.
        cn2_profile : numpy.ndarray, optional
            Zenith Cn2(h) profile. Defaults to the site profile.
        n_samples : int
            FAST Monte Carlo draws (NITER) for smf_fidelity="fast". Ignored for an
            Aperture detector and for smf_fidelity="mean".
        smf_fidelity : str
            The SMF coupling model: "fast" (default, fidelity-1 true modal overlap,
            needs fast-aosim) or "mean" (analytic mean-only, no fade).
        fast_params : dict, optional
            Extra FAST parameters, passed through when smf_fidelity="fast".
        turbulence : bool
            When False, drop every turbulence quantity: an Aperture detector has
            0 dB coupling (no scintillation) and an SMF detector keeps only the
            static mode-match loss (no residual wavefront error, no FAST run). The
            static parts survive; the fade parts that come from turbulence do not.

    Returns:
        Term
            category="coupling".

    Raises:
        ValueError
            If rx_terminal is None or has no detector, or the detector type is
            unknown, or smf_fidelity is unknown.
    '''
    terminal = getattr(scenario, "rx_terminal", None)
    if terminal is None or terminal.detector is None:
        raise ValueError(
            "downlink_coupling_term needs a scenario.rx_terminal with a detector. "
            "Set scenario.rx_terminal = Terminal(..., detector=Aperture() or SMF())."
        )
    # Turbulence off: no scintillation, no residual wavefront, no FAST engine. An
    # Aperture bucket has no static coupling loss (0 dB); an SMF keeps its static
    # mode-match floor. Resolve this before the Cn2 profile, so no turbulence
    # quantity is even built.
    if not turbulence:
        detector = terminal.detector
        if isinstance(detector, Aperture):
            return Term(
                name="receive coupling (aperture)",
                category="coupling",
                mean_db=0.0,
                note="turbulence off: no scintillation",
                meta={"detector": "Aperture", "model": "static"},
                assumptions=Assumptions(
                    beam_type=BEAM_GAUSSIAN,
                    turbulence_regime=REGIME_NA,
                    spectrum=SPECTRUM_NA,
                    validity="Turbulence off: a bucket detector has no static "
                             "coupling loss and no scintillation, so 0 dB."),
            )
        if isinstance(detector, SMF):
            eta_max = _smf_eta_max(detector, terminal.aperture_m,
                                   terminal.wavelength_m)
            return _smf_static_term(eta_max)
        raise ValueError(
            f"unknown detector {type(detector).__name__!r}. Use Aperture or SMF."
        )
    hs = DEFAULT_HS if hs is None else hs
    if cn2_profile is None:
        cn2_profile = default_cn2_profile(scenario.channel.site, hs)

    detector = terminal.detector
    if isinstance(detector, Aperture):
        return _aperture_term(scenario, geometry, hs=hs, cn2_profile=cn2_profile)
    if isinstance(detector, SMF):
        if smf_fidelity == "fast":
            # Fidelity-1: the true LP01 modal overlap from FAST. Lazy import keeps
            # the fast-aosim dependency optional.
            from ..fast import smf_fast_term
            return smf_fast_term(scenario, geometry, hs=hs, cn2_profile=cn2_profile,
                                 n_samples=n_samples, fast_params=fast_params)
        if smf_fidelity == "mean":
            return _smf_mean_term(scenario, geometry, hs=hs, cn2_profile=cn2_profile)
        raise ValueError(
            f"unknown smf_fidelity {smf_fidelity!r}. Use 'fast' or 'mean'."
        )
    raise ValueError(
        f"unknown detector {type(detector).__name__!r}. Use Aperture or SMF."
    )


if __name__ == '__main__':
    import warnings

    from ...scenario import SpaceScenario, Channel
    from ...geometry import CircularOrbit
    from ...terminal import Terminal, Transmitter, TipTilt, AO

    lam = 1550e-9
    hs = DEFAULT_HS
    geom = CircularOrbit(600e3, 60.0)

    def _downlink(ground):
        '''Build a downlink SpaceScenario: tx=space (satellite), rx=ground.'''
        return SpaceScenario(
            ground=ground,
            space=Terminal(aperture_m=0.05, wavelength_m=lam,
                           transmitter=Transmitter(waist_m=0.035)),
            direction="downlink", channel=Channel(altitude_m=600e3))

    # The SMF tests use the analytic mean-only fidelity (the "fast" default needs
    # fast-aosim, which is optional). The Aperture detector ignores smf_fidelity.
    def build(scn):
        return downlink_coupling_term(
            scn, geom, hs=hs, smf_fidelity="mean",
            cn2_profile=default_cn2_profile(scn.channel.site, hs))

    # --- Aperture detector: parity with the plain downlink scintillation ---
    scn_ap = _downlink(Terminal(aperture_m=0.7, wavelength_m=lam,
                                detector=Aperture()))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        t_ap = build(scn_ap)
        t_scint = downlink_scintillation_term(
            scn_ap, geom, cn2_profile=default_cn2_profile(scn_ap.channel.site, hs))
    assert t_ap.category == "coupling"
    assert np.isclose(t_ap.mean_db, t_scint.mean_db)    # exact parity
    q_ap, q_scint = t_ap.quantile_db(0.99), t_scint.quantile_db(0.99)
    assert q_ap is not None and q_scint is not None
    assert np.isclose(q_ap, q_scint)

    # --- SMF, no AO: mean-only Marechal / Dikmelik coupling loss ------------
    scn_smf = _downlink(Terminal(aperture_m=0.7, wavelength_m=lam, detector=SMF()))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        t_smf = build(scn_smf)
    assert t_smf.name == "receive coupling (SMF)"
    assert t_smf.meta["model"] == "mean-only"
    assert t_smf.quantile is None and not t_smf.stochastic   # deterministic
    assert t_smf.meta["coupling_loss_db"] > 0.0
    assert any("Mean-only" in v for v in t_smf.assumptions.violations)
    assert not t_smf.assumptions.ok

    # --- SMF with AO: AO buys back coupling (less mean loss than no AO) -----
    scn_ao = _downlink(Terminal(aperture_m=0.7, wavelength_m=lam, detector=SMF(),
                                compensation=[TipTilt(), AO(200)]))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        t_ao = build(scn_ao)
    assert t_ao.meta["model"] == "mean-only"
    assert t_ao.meta["coupling_loss_db"] > 0.0
    assert float(t_ao.mean_db) < float(t_smf.mean_db), (t_ao.mean_db, t_smf.mean_db)

    # --- Dikmelik-Davidson circular-aperture assumption --------------------
    scn_obsc = _downlink(Terminal(aperture_m=0.7, obscuration_ratio=0.3,
                                  wavelength_m=lam, detector=SMF(),
                                  compensation=[TipTilt(), AO(200)]))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        t_obsc = build(scn_obsc)
    assert any("uniform circular aperture" in v for v in t_obsc.assumptions.violations), \
        t_obsc.assumptions.violations
    assert not any("uniform circular aperture" in v for v in t_ao.assumptions.violations)

    # --- Mean-only: deterministic, no fade ---------------------------------
    rng = np.random.default_rng(0)
    for t in (t_smf, t_ao, t_obsc):
        assert any("Mean-only" in v for v in t.assumptions.violations), \
            t.assumptions.violations
        assert not t.assumptions.ok
        assert not t.stochastic and t.quantile is None
        q99 = t.quantile_db(0.99)   # deterministic term returns its mean, never None
        assert q99 is not None and np.isclose(q99, t.mean_db)
        draws = t.sample_db(1000, rng)
        assert draws.shape == (1000,) and np.allclose(draws, t.mean_db)

    # An elevation sweep broadcasts the mean-only term.
    sweep = CircularOrbit(600e3, np.array([40.0, 60.0, 90.0]))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        t_sweep = downlink_coupling_term(
            scn_ao, sweep, hs=hs, smf_fidelity="mean",
            cn2_profile=default_cn2_profile(scn_ao.channel.site, hs))
    assert np.shape(t_sweep.mean_db) == (3,)
    assert t_sweep.sample_db(100, rng).shape == (100, 3)

    # An unknown fidelity is refused.
    try:
        downlink_coupling_term(scn_smf, geom, hs=hs, smf_fidelity="reciprocity",
                               cn2_profile=default_cn2_profile(scn_smf.channel.site, hs))
        raise AssertionError("unknown smf_fidelity must raise")
    except ValueError as e:
        assert "smf_fidelity" in str(e)

    # --- turbulence=False: static coupling, no turbulence quantity ---------
    ap_off = downlink_coupling_term(scn_ap, geom, turbulence=False)
    assert ap_off.mean_db == 0.0 and ap_off.meta["model"] == "static"
    smf_off = downlink_coupling_term(scn_smf, geom, turbulence=False)
    assert smf_off.meta["model"] == "static" and not smf_off.mean_only
    assert not smf_off.stochastic and smf_off.quantile is None
    assert np.isclose(smf_off.mean_db, -10.0 * np.log10(0.8145))

    print(f"aperture coupling mean = {float(t_ap.mean_db):.4f} dB "
          f"(= scintillation {float(t_scint.mean_db):.4f} dB)")
    print(f"SMF no-AO (mean-only):   coupling loss {float(t_smf.mean_db):.2f} dB  "
          f"eta={t_smf.meta['eta']:.4f}")
    print(f"SMF +AO200 (mean-only):  eta={t_ao.meta['eta']:.4f}  "
          f"coupling loss={t_ao.meta['coupling_loss_db']:.2f} dB")
    print("coupling downlink self-check passed")
