'''
Receive-coupling Term for a downlink terminal (approach A).

This module builds the ONE receive-coupling Term of a downlink receive terminal.
The Compensation stack and the Detector are one physical chain: the residual
wavefront that the Compensation leaves sets the coupling into the Detector. So
this module emits a single Term, not two. The pure physics lives in
olb.turbulence.ao and olb.turbulence.scintillation.

Two detector front ends:

  Aperture (bucket): a power-in-bucket detector integrates the intensity. It is
  phase-insensitive, so the Compensation stack does not change its coupling. The
  turbulence penalty is aperture-averaged scintillation. This term reuses the
  downlink lognormal scintillation Term unchanged, so it reproduces the plain
  downlink behaviour.

  SMF (single-mode fibre): the fibre couples only the field that matches the
  fibre mode. The coupling efficiency eta falls with the residual phase variance
  sigma^2_res that the Compensation stack leaves. Two limits of one overlap
  physics:
    small residual (sigma^2_res < SMF_SMALL_RESIDUAL_LIMIT): extended Marechal,
        eta = eta_max * exp(-sigma^2_res).
    large residual (sigma^2_res >= the limit): the Dikmelik-Davidson uncorrected
        single-mode-fibre coupling curve against D/r0. See _smf_large_residual.
  The mean loss is -10*log10(eta). The Term also carries a fade: it adds the
  aperture-averaged lognormal scintillation fade on top of the constant coupling
  loss. This is the fidelity-0 fade model. Fidelity 1 would sample the
  residual-phase PSD (ResidualWavefront.psd) instead.

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
  SMF coupling / Marechal: V. W. S. Chan and others; extended Marechal
  approximation.
  Uncorrected SMF coupling against D/r0: Y. Dikmelik and F. M. Davidson,
  "Fiber-coupling efficiency for free-space optical communication through
  atmospheric turbulence," Appl. Opt. 44(23), 4946-4952 (2005).
'''

import numpy as np

from dataclasses import replace

from ..results import Term
from ..assumptions import (Assumptions, BEAM_PLANE_WAVE, BEAM_GAUSSIAN,
                           REGIME_WEAK, SPECTRUM_KOLMOGOROV)
from ..terminal import Aperture, SMF, TipTilt, AO, Transmitter
from ..turbulence.profiles import DEFAULT_HS, default_cn2_profile
from ..turbulence.ao import (plane_wave_fried_parameter, apply_compensation,
                            NOLL_PISTON)
from ..links.downlink import downlink_scintillation_term

# Optimal single-mode-fibre coupling. The back-projected fibre mode is a Gaussian.
# The coupling into an unobscured circular aperture is a maximum when the ratio of
# the aperture radius to the aperture-plane Gaussian 1/e^2 radius is
# BETA_SMF_OPT = 1.12, which gives eta_max = 0.8145. So the fibre-mode waist is
# w = (aperture/2)/BETA_SMF_OPT. Source: S. Shaklan and F. Roddier, "Coupling
# starlight into single-mode fiber optics," Appl. Opt. 27(11), 2334 (1988);
# C. Ruilier, SPIE 3350 (1998).
BETA_SMF_OPT = 1.12


def smf_mode_waist(aperture_m):
    '''
    Aperture-plane Gaussian waist of the back-projected single-mode-fibre mode.

    Parameters:
        aperture_m : float
            Receive aperture diameter [m].

    Returns:
        float
            The 1/e^2 Gaussian radius [m] at the aperture (see BETA_SMF_OPT).
    '''
    return (aperture_m / 2.0) / BETA_SMF_OPT

# Residual-variance threshold [rad^2] that selects the SMF coupling limit. Below
# it the extended Marechal approximation holds. Above it the beam is far from a
# flat wavefront, so use the Dikmelik-Davidson uncorrected coupling curve.
SMF_SMALL_RESIDUAL_LIMIT = 1.0

# Effective-D/r0 bound above which even the practical uncorrected coupling curve
# is extrapolated deep turbulence. Flag it. Use the exact Dikmelik-Davidson
# integral (fidelity 1) above this.
SMF_DEEP_TURBULENCE_DR0 = 10.0


def _effective_dr0(sigma2_res):
    '''
    Return the effective D/r0 that carries the residual phase variance.

    Invert the piston-removed Noll relation sigma^2 = NOLL_PISTON*(D/r0)^(5/3).
    For an uncorrected aperture this returns the physical D/r0 exactly. For a
    partly corrected aperture it returns the residual-equivalent D/r0. So the
    uncorrected coupling curve reads the residual, not the raw, turbulence.
    '''
    return (np.asarray(sigma2_res) / NOLL_PISTON) ** (3.0 / 5.0)


def _smf_large_residual(sigma2_res, eta_max):
    '''
    Return the SMF coupling efficiency in the large-residual limit.

    Use the practical uncorrected single-mode-fibre coupling evaluation against
    the effective D/r0. It holds the Dikmelik-Davidson limits: eta -> eta_max as
    D/r0 -> 0, and eta ~ eta_max * (r0/D)^2 for a large D/r0 (only ~ (r0/D)^2
    coherent cells couple into the one fibre mode). The exact Dikmelik-Davidson
    double integral is the fidelity-1 upgrade.

    formula:
        eta = eta_max * [ 1 + (D_eff/r0)^(5/3) ]^(-6/5)
            = eta_max * [ 1 + sigma^2_res / NOLL_PISTON ]^(-6/5)
    The exponent -6/5 gives the large-aperture asymptote eta ~ (r0/D)^2.
    Source: Dikmelik and Davidson, Appl. Opt. 44(23), 2005 (limits and curve).
    '''
    return eta_max * (1.0 + np.asarray(sigma2_res) / NOLL_PISTON) ** (-6.0 / 5.0)


def _smf_coupling_efficiency(sigma2_res, eta_max):
    '''
    Return the SMF coupling efficiency eta from the residual phase variance.

    Select the limit by SMF_SMALL_RESIDUAL_LIMIT. Small residual uses the
    extended Marechal approximation. Large residual uses the uncorrected
    Dikmelik-Davidson coupling curve. The two are limits of one overlap physics.
    They cross over near sigma^2_res = 1 rad^2, where Marechal gives
    eta_max/e = 0.37*eta_max and the Dikmelik-Davidson branch gives about
    0.44*eta_max.

    Parameters:
        sigma2_res : float or numpy.ndarray
            Residual phase variance [rad^2].
        eta_max : float
            Maximum fibre-to-aperture mode match (flat wavefront).

    Returns:
        numpy.ndarray
            Coupling efficiency eta in (0, eta_max].
    '''
    sigma2 = np.asarray(sigma2_res, dtype=float)
    marechal = eta_max * np.exp(-sigma2)
    large = _smf_large_residual(sigma2, eta_max)
    return np.where(sigma2 < SMF_SMALL_RESIDUAL_LIMIT, marechal, large)


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


def _smf_term(scenario, geometry, *, hs, cn2_profile):
    '''
    Build the receive-coupling Term for a single-mode-fibre detector.

    Get r0 from the plane-wave Fried parameter. Get the residual phase variance
    from the compensation stack. Convert it to the coupling efficiency eta, then
    to the mean coupling loss. Add the aperture-averaged lognormal scintillation
    fade on top (the fidelity-0 fade model).
    '''
    terminal = scenario.rx_terminal
    detector = terminal.detector
    wavelength = terminal.wavelength_m
    D = terminal.aperture_m
    elev = geometry.elevation_deg

    r0 = plane_wave_fried_parameter(cn2_profile, hs, wavelength, elev)
    residual = apply_compensation(terminal.compensation, D, r0)
    sigma2_res = residual.variance

    eta = _smf_coupling_efficiency(sigma2_res, detector.eta_max)
    coupling_loss = -10.0 * np.log10(eta)     # positive dB, scalar or per-elevation
    dr0_eff = _effective_dr0(sigma2_res)

    # The fade reuses the downlink aperture-averaged lognormal scintillation. The
    # coupling loss is a constant offset on the fluctuating scintillation. So the
    # Term keeps a closed-form quantile and a sampler.
    scint = downlink_scintillation_term(scenario, geometry, model="lognormal",
                                        aperture_average=True, hs=hs,
                                        cn2_profile=cn2_profile)
    offset = coupling_loss
    base_shape = np.shape(offset)

    mean_db = np.asarray(scint.mean_db) + offset

    def quantile(p):
        return scint.quantile(p) + offset

    def sampler(n, rng):
        return scint.sampler(n, rng) + offset

    assumptions = Assumptions(
        beam_type=BEAM_PLANE_WAVE,
        turbulence_regime=REGIME_WEAK,
        spectrum=SPECTRUM_KOLMOGOROV,
        validity="Single-mode-fibre coupling. Extended Marechal for a small "
                 "residual; Dikmelik-Davidson uncorrected coupling for a large "
                 "residual. The Dikmelik-Davidson coupling assumes a uniform "
                 "circular aperture with no central obscuration. FIDELITY-0 FADE: "
                 "the fade uses the aperture-averaged lognormal scintillation as a "
                 "stand-in for the residual-coupling fluctuation, so every fade "
                 "margin EXCLUDES the residual fibre-coupling fade and is "
                 "optimistic on the fibre side. Weak fluctuation: sigma2_I < 0.25.",
    )
    # Carry over any scintillation-side flag (weak-fluctuation, obscuration).
    if scint.assumptions is not None:
        for reason in scint.assumptions.violations:
            assumptions.flag(reason)
    # FIDELITY-0 FADE. The fade is aperture-averaged scintillation shifted by a
    # CONSTANT mean coupling loss (eta at the mean residual variance). The true
    # time-varying fibre-coupling fluctuation -- the residual wavefront coupling
    # into the single fibre mode (AO/tip-tilt residual jitter, speckle) -- is NOT
    # modelled. So the quantile and every fade margin (for example the 99% link
    # margin) EXCLUDE the residual fibre-coupling fade and are OPTIMISTIC on the
    # fibre side. Always flag it. The fidelity-1 upgrade samples the residual-
    # phase PSD (olb.turbulence.ao.ResidualWavefront.psd).
    assumptions.flag(
        "Fidelity-0 fade: the fibre-coupling loss is a constant offset on the "
        "aperture-averaged scintillation fade. The 99% (and every) fade margin "
        "EXCLUDES the residual fibre-coupling fluctuation (speckle / AO-residual "
        "jitter into the single mode), so the margin is optimistic on the fibre "
        "side. Fidelity 1 samples ResidualWavefront.psd (olb.turbulence.ao)."
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
            "the practical uncorrected coupling curve is extrapolated. Use the "
            "exact Dikmelik-Davidson integral (fidelity 1)."
        )

    note = (f"SMF coupling, eta_max={detector.eta_max:g}, "
            f"n_comp_modes={residual.n_modes}")
    return Term(
        name="receive coupling (SMF)",
        category="coupling",
        mean_db=float(mean_db) if base_shape == () else mean_db,
        sampler=sampler,
        quantile=quantile,
        note=note,
        meta={
            "detector": "SMF",
            "eta": float(eta) if base_shape == () else np.asarray(eta),
            "coupling_loss_db": float(offset) if base_shape == () else np.asarray(offset),
            "sigma2_res": float(sigma2_res) if np.ndim(sigma2_res) == 0 else np.asarray(sigma2_res),
            "effective_D_over_r0": float(dr0_eff) if np.ndim(dr0_eff) == 0 else np.asarray(dr0_eff),
            "r0_m": float(r0) if np.ndim(r0) == 0 else np.asarray(r0),
            "n_comp_modes": residual.n_modes,
        },
        assumptions=assumptions,
    )


def _smf_reciprocity_term(scenario, geometry, *, hs, cn2_profile, n_samples):
    '''
    Build the SMF receive-coupling Term with the reciprocity Strehl proxy.

    Approach (no adaptive optics; tip-tilt at most). By Shapiro reciprocity the
    downlink fibre-coupling Strehl equals the on-axis Strehl of the back-projected
    fibre-mode Gaussian launched up the same turbulent column. The Dios coupled-
    flux (see olb.links.uplink) already gives that fluctuating on-axis Strehl
    S(t). The coupling efficiency is

        eta(t) = eta_max * S(t),

    where eta_max is the static (flat-wavefront) mode match and S(t) carries the
    turbulence coupling fade -- most importantly the beam-wander (angle-of-arrival)
    tip-tilt fade that the fidelity-0 model drops. The fibre mode is a Gaussian of
    aperture-plane waist w = smf_mode_waist(aperture) (see BETA_SMF_OPT).

    IMPORTANT. This is a Strehl PROXY, NOT a true LP01 modal overlap. The on-axis
    Strehl and the coherent mode overlap agree well for a tip-tilt-dominated fade,
    but not for scintillation and high-order aberration. So the Term is best for
    the no-AO / tip-tilt regime. The Monte-Carlo-only fidelity-1 upgrade is the
    true Gaussian-mode overlap (for example the FAST tool). Because the coupled-
    flux is Monte-Carlo-only, this Term gives a sampler and sets quantile=None.

    Parameters:
        scenario : Scenario
            The downlink case. Reads rx_terminal (aperture, wavelength, detector,
            compensation) and the site Cn2 profile.
        geometry : CircularOrbit or TLEPass
            The link geometry (elevation, slant range).
        hs : numpy.ndarray
            Height grid [m].
        cn2_profile : numpy.ndarray
            Zenith Cn2(h) profile matching hs (fast-free, from default_cn2_profile).
        n_samples : int
            Monte Carlo draws for the reciprocal coupled-flux mean estimate.

    Returns:
        Term
            name="receive coupling (SMF)", category="coupling", Monte-Carlo-only.
    '''
    # Import here to break the coupling <-> uplink import cycle.
    from ..links.uplink import uplink_turbulence_term

    rx = scenario.rx_terminal
    detector = rx.detector
    eta_max = detector.eta_max
    w_fib = smf_mode_waist(rx.aperture_m)
    floor_db = -10.0 * np.log10(eta_max)   # static mode-match loss, a constant

    # The reciprocal uplink: launch the back-projected fibre-mode Gaussian from
    # the receive aperture, collimated, up the same column. Only the transmit
    # waist matters to the coupled-flux; the satellite receiver is on-axis.
    recip_ground = replace(rx, transmitter=Transmitter(waist_m=w_fib))
    recip_scn = replace(scenario, direction="uplink", ground=recip_ground)
    turb = uplink_turbulence_term(recip_scn, geometry, n_samples=n_samples,
                                  cn2_profile=cn2_profile)

    base_shape = np.shape(turb.mean_db)
    mean_db = np.asarray(turb.mean_db) + floor_db

    def sampler(n, rng):
        return turb.sampler(n, rng) + floor_db

    assumptions = Assumptions(
        beam_type=BEAM_GAUSSIAN,
        turbulence_regime=REGIME_WEAK,
        spectrum=SPECTRUM_KOLMOGOROV,
        validity="Reciprocity Strehl proxy for single-mode-fibre coupling: "
                 "eta = eta_max * S(t), where S(t) is the on-axis Strehl of the "
                 "back-projected fibre-mode Gaussian launched up the same column "
                 "(Dios coupled-flux, olb.links.uplink). The fibre mode is a "
                 f"Gaussian of aperture-plane waist w = (aperture/2)/{BETA_SMF_OPT} "
                 "(Shaklan-Roddier optimum, eta_max=0.8145). Weak fluctuation: "
                 "sigma2_x < 0.6.",
    )
    # STREHL PROXY, not a true modal overlap. Always flag it.
    assumptions.flag(
        "Strehl proxy: the fade is eta_max * on-axis Strehl, NOT a coherent LP01 "
        "modal overlap. It captures the tip-tilt / angle-of-arrival coupling fade "
        "but approximates scintillation and high-order coupling. The fidelity-1 "
        "upgrade is the true Gaussian-mode overlap (for example the FAST tool)."
    )
    # Reciprocity assumes the up and down columns coincide. Always flag it.
    assumptions.flag(
        "Reciprocity assumes the up-leg and down-leg columns coincide. For a LEO "
        "the point-ahead angle can exceed the isoplanatic angle, so the coupling "
        "statistics are approximate."
    )
    # Tip-tilt receiver: the reciprocity fade still carries the full wander.
    if any(isinstance(c, TipTilt) for c in rx.compensation):
        assumptions.flag(
            "Tip-tilt correction present: the reciprocity fade still includes the "
            "full angle-of-arrival (beam-wander) fluctuation, which a tip-tilt "
            "tracker partly removes. So the fade is conservative for a tip-tilt "
            "receiver. Wander removal is not yet modelled."
        )
    # Carry the weak-fluctuation flag from the reciprocal up-leg.
    wf = turb.meta.get("weak_fluctuation_valid")
    if wf is not None and not np.all(wf):
        assumptions.flag(
            "Rytov weak-fluctuation limit exceeded on the reciprocal up-leg; "
            "the scintillation approaches saturation and the Strehl is not "
            "trustworthy."
        )

    note = (f"SMF reciprocity Strehl proxy, w_fibre={w_fib:.3g} m, "
            f"eta_max={eta_max:g}")
    return Term(
        name="receive coupling (SMF)",
        category="coupling",
        mean_db=float(mean_db) if base_shape == () else mean_db,
        sampler=sampler,
        quantile=None,   # MC-only: coupled-flux has no closed form -> monte_carlo()
        note=note,
        meta={
            "detector": "SMF",
            "model": "reciprocity-strehl",
            "eta_max": float(eta_max),
            "w_fibre_mode_m": float(w_fib),
            "beta_smf_opt": BETA_SMF_OPT,
            "floor_db": float(floor_db),
            "n_samples": n_samples,
            "sigma2_x": turb.meta.get("sigma2_x"),
        },
        assumptions=assumptions,
    )


def rx_coupling_term(scenario, geometry, *, hs=None, cn2_profile=None,
                     n_samples=2000, smf_fidelity="reciprocity", fast_params=None):
    '''
    Build the ONE receive-coupling Term of a downlink receive terminal.

    Read scenario.rx_terminal. Dispatch on the detector type and the correction.
    An Aperture detector reuses the downlink aperture-averaged scintillation, so it
    is parity with the plain downlink. An SMF detector splits by the correction:

    - No adaptive optics (empty stack or tip-tilt only): the reciprocity Strehl
      proxy. It reuses the Dios coupled-flux to carry the tip-tilt / angle-of-
      arrival coupling fade. Monte-Carlo-only. See _smf_reciprocity_term.
    - Adaptive optics present: the extended-Marechal / Dikmelik-Davidson coupling
      from the residual wavefront, with the aperture-averaged scintillation as the
      fidelity-0 fade. See _smf_term.

    Parameters:
        scenario : Scenario
            Reads rx_terminal, link.wavelength_m, and the site Cn2 profile.
        geometry : CircularOrbit or TLEPass
            Reads elevation_deg. A scalar elevation gives a scalar Term.
        hs : numpy.ndarray, optional
            Heights above the ground station [m]. Defaults to DEFAULT_HS.
        cn2_profile : numpy.ndarray, optional
            Zenith Cn2(h) profile. Defaults to the site profile.
        n_samples : int
            Monte Carlo draws for the reciprocity Strehl proxy (SMF, no AO) and
            for the FAST fidelity-1 term (its NITER).
        smf_fidelity : str
            The SMF coupling model. "reciprocity" (default) uses the Strehl proxy
            for a no-AO fibre and the Dikmelik/Marechal model for an AO fibre.
            "fast" uses the FAST fidelity-1 true modal overlap (needs fast-aosim).
        fast_params : dict, optional
            Extra FAST parameters, passed through when smf_fidelity="fast".

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
            "rx_coupling_term needs a scenario.rx_terminal with a detector. "
            "Set scenario.rx_terminal = Terminal(..., detector=Aperture() or SMF())."
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
            from .coupling_fast import smf_fast_term
            return smf_fast_term(scenario, geometry, hs=hs, cn2_profile=cn2_profile,
                                 n_samples=n_samples, fast_params=fast_params)
        if smf_fidelity != "reciprocity":
            raise ValueError(
                f"unknown smf_fidelity {smf_fidelity!r}. Use 'reciprocity' or 'fast'."
            )
        # No AO -> reciprocity Strehl proxy (carries the tip-tilt fade). AO
        # present -> extended-Marechal / Dikmelik coupling (fidelity-0 fade).
        has_ao = any(isinstance(c, AO) for c in terminal.compensation)
        if has_ao:
            return _smf_term(scenario, geometry, hs=hs, cn2_profile=cn2_profile)
        return _smf_reciprocity_term(scenario, geometry, hs=hs,
                                     cn2_profile=cn2_profile, n_samples=n_samples)
    raise ValueError(
        f"unknown detector {type(detector).__name__!r}. Use Aperture or SMF."
    )


if __name__ == '__main__':
    import warnings

    from ..scenario import Scenario, Channel
    from ..geometry import CircularOrbit
    from ..terminal import Terminal, Transmitter, TipTilt, AO

    lam = 1550e-9
    hs = DEFAULT_HS
    geom = CircularOrbit(600e3, 60.0)

    def _downlink(ground):
        '''Build a downlink Scenario: tx=space (satellite), rx=ground.'''
        return Scenario(
            ground=ground,
            space=Terminal(aperture_m=0.05, wavelength_m=lam,
                           transmitter=Transmitter(waist_m=0.035)),
            direction="downlink", channel=Channel(altitude_m=600e3))

    # --- SMF efficiency limits --------------------------------------------
    # Small residual -> Marechal. Large residual -> Dikmelik-Davidson branch.
    eta0 = _smf_coupling_efficiency(0.0, 0.8145)
    assert np.isclose(eta0, 0.8145)                     # flat wavefront -> eta_max
    assert _smf_coupling_efficiency(0.5, 0.8145) < 0.8145
    # eta falls as the residual rises, in both branches.
    grid = np.array([0.0, 0.3, 0.9, 2.0, 8.0, 30.0])
    etas = _smf_coupling_efficiency(grid, 0.8145)
    assert np.all(np.diff(etas) < 0.0), etas
    # Large-aperture asymptote: eta ~ (r0/D)^2, so eta ~ sigma2_res^(-6/5).
    e1 = _smf_large_residual(2000.0, 0.8145)
    e2 = _smf_large_residual(4000.0, 0.8145)
    slope = np.log(e2 / e1) / np.log(4000.0 / 2000.0)
    assert abs(slope - (-6.0 / 5.0)) < 1e-3, slope

    def build(scn):
        return rx_coupling_term(scn, geom, hs=hs, n_samples=800,
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
    assert np.isclose(t_ap.quantile_db(0.99), t_scint.quantile_db(0.99))

    # --- SMF, no AO: reciprocity Strehl proxy (Monte-Carlo-only) -----------
    # A 0.7 m fibre receiver with no correction. The reciprocity path launches the
    # back-projected fibre mode (waist = (aperture/2)/1.12) up the same column and
    # reads the on-axis Strehl.
    scn_smf = _downlink(Terminal(aperture_m=0.7, wavelength_m=lam, detector=SMF()))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        t_smf = build(scn_smf)
    assert t_smf.name == "receive coupling (SMF)"
    assert t_smf.meta["model"] == "reciprocity-strehl"
    assert t_smf.quantile is None and t_smf.stochastic      # Monte-Carlo-only
    # Fibre-mode sizing relative to the aperture (deliberate).
    assert np.isclose(t_smf.meta["w_fibre_mode_m"], (0.7 / 2) / BETA_SMF_OPT)
    assert t_smf.meta["beta_smf_opt"] == BETA_SMF_OPT
    # The mean loss is at least the static mode-match floor.
    assert t_smf.mean_db >= t_smf.meta["floor_db"] - 1e-9
    # Strehl-proxy and reciprocity caveats are always flagged, so it is never ok.
    assert any("Strehl proxy" in v for v in t_smf.assumptions.violations)
    assert any("point-ahead" in v for v in t_smf.assumptions.violations)
    assert not t_smf.assumptions.ok

    # --- SMF, tip-tilt only: still reciprocity, plus a conservative flag ----
    scn_tt = _downlink(Terminal(aperture_m=0.7, wavelength_m=lam, detector=SMF(),
                                compensation=[TipTilt()]))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        t_tt = build(scn_tt)
    assert t_tt.meta["model"] == "reciprocity-strehl"
    assert any("Tip-tilt correction present" in v for v in t_tt.assumptions.violations)

    # --- SMF with AO: extended-Marechal / Dikmelik coupling (fidelity-0) ----
    scn_ao = _downlink(Terminal(aperture_m=0.7, wavelength_m=lam, detector=SMF(),
                                compensation=[TipTilt(), AO(200)]))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        t_ao = build(scn_ao)
    assert t_ao.meta.get("model") is None                  # the Marechal/Dikmelik path
    assert t_ao.meta["coupling_loss_db"] > 0.0
    # AO buys back coupling: the AO mean loss is far below the no-AO reciprocity mean.
    assert float(t_ao.mean_db) < float(t_smf.mean_db), (t_ao.mean_db, t_smf.mean_db)

    # --- Dikmelik-Davidson circular-aperture assumption (AO path) ----------
    # A central obscuration on the receive aperture flags a violation naming the
    # coupling curve; an unobscured aperture does not. This flag lives on the AO
    # (Marechal/Dikmelik) path.
    scn_obsc = _downlink(Terminal(aperture_m=0.7, obscuration_ratio=0.3,
                                  wavelength_m=lam, detector=SMF(),
                                  compensation=[TipTilt(), AO(200)]))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        t_obsc = build(scn_obsc)
    assert any("uniform circular aperture" in v for v in t_obsc.assumptions.violations), \
        t_obsc.assumptions.violations
    assert not any("uniform circular aperture" in v for v in t_ao.assumptions.violations)

    # --- Fidelity-0 fade: flagged on the AO (Marechal/Dikmelik) path -------
    # An AO SMF term flags that the fade excludes the residual fibre-coupling
    # fluctuation, so it is never fully ok.
    for t in (t_ao, t_obsc):
        assert any("Fidelity-0 fade" in v for v in t.assumptions.violations), \
            t.assumptions.violations
        assert not t.assumptions.ok

    # The AO SMF term is stochastic with a closed-form quantile deeper than the mean.
    rng = np.random.default_rng(0)
    assert t_ao.stochastic and t_ao.quantile is not None
    assert t_ao.quantile_db(0.99) > t_ao.mean_db
    draws = t_ao.sample_db(50_000, rng)
    assert abs(draws.mean() - t_ao.mean_db) < 0.05, (draws.mean(), t_ao.mean_db)

    # The reciprocity SMF term samples too (Monte-Carlo-only, no quantile).
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r_draws = t_smf.sample_db(1000, rng)
    assert r_draws.shape == (1000,) and np.all(np.isfinite(r_draws))

    # An elevation sweep broadcasts on the AO path.
    sweep = CircularOrbit(600e3, np.array([40.0, 60.0, 90.0]))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        t_sweep = rx_coupling_term(scn_ao, sweep, hs=hs,
                                   cn2_profile=default_cn2_profile(scn_ao.channel.site, hs))
    assert np.shape(t_sweep.mean_db) == (3,)
    assert t_sweep.sample_db(100, rng).shape == (100, 3)

    print(f"aperture coupling mean = {float(t_ap.mean_db):.4f} dB "
          f"(= scintillation {float(t_scint.mean_db):.4f} dB)")
    print(f"SMF no-AO (reciprocity): mean {float(t_smf.mean_db):.2f} dB  "
          f"w_fibre={t_smf.meta['w_fibre_mode_m']:.3f} m  "
          f"floor={t_smf.meta['floor_db']:.2f} dB")
    print(f"SMF +AO200 (Dikmelik):   eta={t_ao.meta['eta']:.4f}  "
          f"coupling loss={t_ao.meta['coupling_loss_db']:.2f} dB  "
          f"sigma2_res={t_ao.meta['sigma2_res']:.4f}")
    print("self-check passed")
