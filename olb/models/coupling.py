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
  fibre mode. Two fidelities:
    smf_fidelity="fast" (the default, and the ONLY statistical model): the true
        LP01 modal overlap from FAST. It gives the mean, the quantile, and the
        fade. See olb.models.coupling_fast.
    smf_fidelity="mean": a cheap analytic MEAN-ONLY estimate, for when only the
        expected coupling loss is wanted (no fade). The coupling efficiency eta
        falls with the residual phase variance sigma^2_res that the Compensation
        stack leaves. Two limits of one overlap physics:
          small residual (sigma^2_res < SMF_SMALL_RESIDUAL_LIMIT): extended
              Marechal, eta = eta_max * exp(-sigma^2_res).
          large residual (sigma^2_res >= the limit): the Dikmelik-Davidson
              uncorrected single-mode-fibre coupling curve against D/r0. See
              _smf_large_residual.
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
  SMF coupling / Marechal: V. W. S. Chan and others; extended Marechal
  approximation.
  Uncorrected SMF coupling against D/r0: Y. Dikmelik and F. M. Davidson,
  "Fiber-coupling efficiency for free-space optical communication through
  atmospheric turbulence," Appl. Opt. 44(23), 4946-4952 (2005).
'''

import numpy as np

from ..results import Term
from ..assumptions import (Assumptions, BEAM_PLANE_WAVE, BEAM_GAUSSIAN,
                           REGIME_WEAK, REGIME_NA, SPECTRUM_KOLMOGOROV,
                           SPECTRUM_NA)
from ..terminal import Aperture, SMF, MMF, TipTilt, AO
from ..beam import free_space_radius
from ..turbulence.profiles import DEFAULT_HS, default_cn2_profile
from ..turbulence.ao import (plane_wave_fried_parameter, apply_compensation,
                            NOLL_PISTON)
from ..turbulence.gaussian_fried import gaussian_fried_parameter_profile
from ..turbulence.angle_of_arrival import wander_arrival_angle_variance
from ..links.downlink import downlink_scintillation_term

# dB per (r^2 / w^2), from the exponential Gaussian power falloff. Same constant
# as olb.models.pointing (20/ln10). Source: Andrews and Phillips, 2nd ed. (2005),
# DOI 10.1117/3.626196 (Gaussian power in a circular region).
_K = 20.0 / np.log(10.0)

# Residual-variance threshold [rad^2] that selects the SMF coupling limit. Below
# it the extended Marechal approximation holds. Above it the beam is far from a
# flat wavefront, so use the Dikmelik-Davidson uncorrected coupling curve.
SMF_SMALL_RESIDUAL_LIMIT = 1.0

# Effective-D/r0 bound above which even the practical uncorrected coupling curve
# is extrapolated deep turbulence. Flag it. Use the exact Dikmelik-Davidson
# integral (fidelity 1) above this.
SMF_DEEP_TURBULENCE_DR0 = 10.0

# Optimal single-mode-fibre coupling parameter a = w_m/w_s. It maximises the
# mode-overlap eta_max(a)=0.8145. Source: Shaklan and Roddier, Appl. Opt. 27
# (1988) 2334, DOI 10.1364/AO.27.002334.
SMF_OPTIMAL_A = 1.12
# Default mode field radius [m] for optimal_focus, from SMF-28 at 1550 nm (mode
# field diameter 10.4 um).
SMF28_MODE_FIELD_RADIUS_M = 5.2e-6


def _smf_optics(detector, D, wavelength):
    '''
    Resolve the (focal_length, mode_field_radius) of a single-mode-fibre detector.

    Honour SMF.optimal_focus. When it is True the model assumes the optimal
    coupling parameter a=SMF_OPTIMAL_A and derives the focal length from the mode
    field radius and the aperture:
        f = pi*(D/2)*w_m / (lambda*SMF_OPTIMAL_A).
    A missing mode field radius defaults to the SMF-28 value. An explicit
    focal_length_m always wins. Source: Shaklan and Roddier, Appl. Opt. 27 (1988)
    2334, DOI 10.1364/AO.27.002334 (the a parameter).

    Returns:
        tuple
            (focal_length_m, mode_field_radius_m). An entry is None when the
            detector does not set it and optimal_focus is False.
    '''
    w_m = detector.mode_field_radius_m
    f = detector.focal_length_m
    if getattr(detector, "optimal_focus", False):
        if w_m is None:
            w_m = SMF28_MODE_FIELD_RADIUS_M
        if f is None:
            # Optimal focus: pick f so a=SMF_OPTIMAL_A (the eta_max peak).
            f = np.pi * (D / 2.0) * w_m / (wavelength * SMF_OPTIMAL_A)
    return f, w_m


def _mmf_focal_length(detector, D, wavelength):
    '''
    Resolve the focal length of a multimode-fibre detector, honouring optimal_focus.

    An explicit focal_length_m always wins. With optimal_focus the model matches
    the spot to the core: it picks f so the spot radius is the core radius over
    SMF_OPTIMAL_A, that is a_core/w_s = SMF_OPTIMAL_A. So
        f = pi*(D/2)*a_core / (lambda*SMF_OPTIMAL_A).
    This is a geometric spot-to-core match (about 92% static capture), NOT a
    mode-overlap optimum: a shorter f captures more, up to the practical numerical
    aperture. Source of the a parameter: Shaklan and Roddier, Appl. Opt. 27 (1988)
    2334, DOI 10.1364/AO.27.002334.

    Returns:
        float or None
            The focal length [m], or None when the detector sets neither an
            explicit focal length nor optimal_focus.
    '''
    if detector.focal_length_m is not None:
        return detector.focal_length_m
    if getattr(detector, "optimal_focus", False):
        # NOTE: for now, optimal focus is based on SMF28 core size
        core_radius = 5.2e-6
        return np.pi * (D / 2.0) * core_radius / (wavelength * SMF_OPTIMAL_A)
    return None


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


def smf_eta_max_from_a(a):
    '''
    Return the flat-wavefront single-mode-fibre coupling eta_max from a.

    The coupling parameter is a = w_m/w_s = pi*(D/2)*w_m/(lambda*f), with w_m the
    fibre mode field radius, w_s the focal spot radius, D the aperture diameter,
    and f the focal length. The mode-overlap of an unobscured, uniformly
    illuminated circular aperture with a Gaussian fibre mode is

        eta_max(a) = 2 * [ (1 - exp(-a^2)) / a ]^2.

    It peaks at eta_max=0.8145 near a=1.12, and it falls on both sides. Source:
    Shaklan and Roddier, Applied Optics 27, 2334 (1988), DOI
    10.1364/AO.27.002334 (also Ruilier, Proc. SPIE 3350, 1998, DOI
    10.1117/12.317094; and Dikmelik and Davidson, Applied Optics 44(23), 4946
    (2005), DOI 10.1364/AO.44.004946).

    Parameters:
        a : float or numpy.ndarray
            The coupling parameter a = pi*(D/2)*w_m/(lambda*f).

    Returns:
        float or numpy.ndarray
            The flat-wavefront coupling efficiency eta_max in (0, 0.8145].
    '''
    a = np.asarray(a, dtype=float)
    return 2.0 * ((1.0 - np.exp(-a ** 2)) / a) ** 2


def _smf_eta_max(detector, D, wavelength):
    '''
    Return the flat-wavefront eta_max for a single-mode-fibre detector.

    When the detector has no focal length, use the eta_max field (today's
    behaviour). When it has a focal length, derive a and eta_max(a) from the
    optics. A focal length needs a mode field radius. See smf_eta_max_from_a.
    With optimal_focus and no explicit focal length, a=SMF_OPTIMAL_A, so
    eta_max is the peak 0.8145.
    '''
    if getattr(detector, "optimal_focus", False) and detector.focal_length_m is None:
        return float(smf_eta_max_from_a(SMF_OPTIMAL_A))
    if detector.focal_length_m is None:
        return detector.eta_max
    if detector.mode_field_radius_m is None:
        raise ValueError(
            "SMF.focal_length_m needs SMF.mode_field_radius_m to derive the "
            "coupling parameter a."
        )
    a = (np.pi * (D / 2.0) * detector.mode_field_radius_m
         / (wavelength * detector.focal_length_m))
    return float(smf_eta_max_from_a(a))


# --- received tip-tilt and the focal-plane walk-off fade ---------------------

def _received_tiptilt_variance(scenario, *, n_grid, turbulence=True):
    '''
    Return the radial (2-axis) received tip-tilt variance at the fibre [rad^2].

    Sum two contributions:
      A. The beam-wander arrival tilt over the horizontal path. Reuse
         olb.turbulence.angle_of_arrival.wander_arrival_angle_variance. A receive
         tip-tilt or AO stage tracks it out, so the code gates it: a TipTilt or an
         AO stage in the rx compensation stack removes the wander tilt (residual
         ~ 0). Without a tip-tilt corrector, the full wander tilt reaches the
         fibre.
         TODO (deferred): a finite-bandwidth tracking loop leaves a residual
         tilt. There is no tracking-bandwidth model yet, so the correction is
         all-or-nothing.
      B. The mechanical pointing (tracking) jitter of the RECEIVE terminal. A
         per-axis 1-sigma jitter sigma adds 2*sigma^2 to the radial (2-axis)
         variance.

    The rx mechanical jitter billed here is the RECEIVE terminal jitter. It is
    NOT the transmit-side pointing loss (olb.models.pointing uses the TRANSMIT
    terminal against the far aperture). So the two do not double-count.

    The beam-wander arrival tilt (A) is a turbulence quantity; the mechanical
    jitter (B) is not. With turbulence=False the wander tilt drops to zero and only
    the receive mechanical jitter reaches the fibre, so the fibre walk-off fade is
    the jitter alone. See the master turbulence switch on the link budgets.

    Returns:
        tuple
            (sigma2_theta_radial, meta) : the radial variance [rad^2] and a meta
            dict with the split and the tracking flag.
    '''
    rx = scenario.rx_terminal

    if turbulence:
        tx = scenario.tx_terminal
        if tx.transmitter is None:
            raise ValueError(
                "the received tip-tilt needs a launch beam: set the near terminal "
                "transmitter = Transmitter(waist_m=...)."
            )
        w0 = tx.transmitter.waist_m
        divergence = tx.transmitter.divergence_rad
        wavelength = rx.wavelength_m
        L = float(scenario.channel.path_length_m)
        cn2 = float(scenario.channel.cn2)

        hs = np.linspace(0.0, L, int(n_grid))
        cn2_slant = np.full_like(hs, cn2)
        w_profile = free_space_radius(w0, hs, divergence, wavelength)
        sigma2_wander = wander_arrival_angle_variance(L, cn2_slant, w_profile, hs)

        # A tip-tilt (or AO, which includes tilt) tracking loop removes the wander
        # arrival tilt. The all-or-nothing gate is the ponytail model (no bandwidth).
        tracks = any(isinstance(s, (TipTilt, AO)) for s in rx.compensation)
        if tracks:
            sigma2_wander = 0.0
    else:
        # Turbulence off: no beam-wander arrival tilt, only the mechanical jitter.
        sigma2_wander = 0.0
        tracks = False

    sigma2_jitter = 2.0 * rx.pointing_jitter_rad ** 2   # per-axis -> 2 axes
    sigma2_total = sigma2_wander + sigma2_jitter
    meta = {
        "sigma2_theta_radial": sigma2_total,
        "sigma2_wander": sigma2_wander,
        "sigma2_jitter": sigma2_jitter,
        "wander_tracked": bool(tracks),
    }
    return sigma2_total, meta


def _walkoff_faces(f, w_eff, sigma2_theta_radial):
    '''
    Return the (mean_db, quantile, sampler) of the focal-plane walk-off fade.

    A received tip-tilt moves the focal spot on the fibre tip by f*theta. The
    radial displacement dx has two i.i.d. Gaussian axes, so dx^2 is exponential.
    The captured-power fraction against the fibre feature of size w_eff is
        h(dx) = exp(-2 dx^2 / w_eff^2)   ->   loss_db = _K * 2 dx^2 / w_eff^2.
    So the loss in dB is exponential:
        loss_db ~ Exponential(mean = _K * f^2 * sigma2_theta_radial / w_eff^2).
    This is the same closed form as olb.models.pointing. w_eff is the effective
    coupling scale that the caller gives. For a single-mode fibre it is
    sqrt(w_s^2+w_m^2), the two-Gaussian overlap of the focal spot (radius w_s) and
    the fibre mode (radius w_m). For a multimode fibre it is the core radius.
    Source: Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196 (Gaussian
    power falloff and 2-D Gaussian jitter); two-Gaussian mode overlap, Shaklan and
    Roddier, Appl. Opt. 27 (1988) 2334, DOI 10.1364/AO.27.002334.
    '''
    # Per-axis focal displacement variance: sigma_dx^2 = f^2 * (radial/2).
    mean = _K * f ** 2 * float(sigma2_theta_radial) / w_eff ** 2

    def quantile(p):
        return -mean * np.log(1.0 - p)      # inverse exponential CDF

    def sampler(n, rng):
        if mean <= 0.0:
            return np.zeros(n)
        return rng.exponential(scale=mean, size=n)

    return mean, quantile, sampler


def _smf_static_term(eta_max):
    '''
    Turbulence-off single-mode-fibre coupling Term: the static mode-match loss.

    With turbulence off there is no residual wavefront error, so eta = eta_max and
    the coupling loss is the fixed mode-match floor -10*log10(eta_max). This is a
    real static optical loss, NOT a turbulence quantity, so the Term is
    DETERMINISTIC but NOT mean-only: it does not lock the budget out of a fade, so a
    jitter walk-off Term can still carry the fade. Source of eta_max: Shaklan and
    Roddier, Applied Optics 27, 2334 (1988), DOI 10.1364/AO.27.002334.
    '''
    coupling_loss = -10.0 * np.log10(eta_max)
    assumptions = Assumptions(
        beam_type=BEAM_GAUSSIAN,
        turbulence_regime=REGIME_NA,
        spectrum=SPECTRUM_NA,
        validity="Turbulence off: static single-mode-fibre mode-match coupling only "
                 "(eta_max). No residual wavefront error and no fade.",
    )
    return Term(
        name="receive coupling (SMF)",
        category="coupling",
        mean_db=float(coupling_loss),
        note=f"SMF static coupling (turbulence off), eta_max={eta_max:g}",
        meta={"detector": "SMF", "model": "static", "eta_max": float(eta_max),
              "eta": float(eta_max), "coupling_loss_db": float(coupling_loss)},
        assumptions=assumptions,
    )


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
    modal overlap) use smf_fidelity="fast" (olb.models.coupling_fast).

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
    r0 = plane_wave_fried_parameter(cn2_profile, hs, wavelength, elev)
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


def rx_coupling_term(scenario, geometry, *, hs=None, cn2_profile=None,
                     n_samples=2000, smf_fidelity="fast", fast_params=None,
                     turbulence=True):
    '''
    Build the ONE receive-coupling Term of a downlink receive terminal.

    Read scenario.rx_terminal. Dispatch on the detector type. An Aperture detector
    reuses the downlink aperture-averaged scintillation, so it is parity with the
    plain downlink. An SMF detector picks the fidelity with smf_fidelity:

    - "fast" (the default): the fidelity-1 true LP01 modal overlap from FAST. It
      gives the mean, the quantile, and the fade. Needs fast-aosim. See
      olb.models.coupling_fast.
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
            "rx_coupling_term needs a scenario.rx_terminal with a detector. "
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
            from .coupling_fast import smf_fast_term
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


def terrestrial_smf_coupling_term(scenario, geometry, *, n_grid=64,
                                  drop_tiptilt=False, turbulence=True):
    '''
    Fidelity-0 single-mode-fibre coupling Term for a terrestrial link.

    Compute the MEAN fibre-coupling loss of the received Gaussian beam at the far
    aperture. Get r0 from the horizontal Gaussian-beam Fried parameter (weak
    turbulence). Get the residual phase variance from the compensation stack
    (tip-tilt, AO). Convert it to the coupling efficiency eta, then to the mean
    coupling loss -10*log10(eta). The Term is MEAN-ONLY: it has no fade, so it
    locks the budget to fidelity 0.

    IMPORTANT (effective-r0 weak-turbulence approximation, well flagged): the Noll
    residual and the Dikmelik-Davidson coupling are PLANE-WAVE, Kolmogorov,
    phase-only forms. This Term evaluates them at the GAUSSIAN-beam r0. That
    substitution holds ONLY in weak turbulence, where the wavefront is
    phase-dominated. It ignores the beam-wave amplitude scintillation, the beam
    wander, and the near-field beam curvature that a finite horizontal beam
    carries. For those a fidelity-2 split-step beam-propagation model is needed
    (FAST does not model a near-field finite beam). See the Term assumptions.

    Parameters:
        scenario : TerrestrialScenario
            tx = near (its Transmitter waist launches the beam); rx = far (its
            SMF detector, aperture, obscuration, and compensation stack).
        geometry : HorizontalPath
            Unused here (path length and Cn2 come from the channel). Kept for the
            f(scenario, geometry) -> Term signature.
        n_grid : int
            Points on the constant-Cn2 path grid for the r0 integral.
        drop_tiptilt : bool
            Remove the tip-tilt (Noll modes 2 and 3) from the residual phase
            variance. Set True when the budget also adds the receive tip-tilt
            walk-off Term (smf_walkoff_term). The walk-off Term then carries the
            tip-tilt coupling loss, and this Term keeps the HIGHER-ORDER residual
            only. So the tip-tilt is not counted two times.
        turbulence : bool
            When False, drop the residual wavefront error (no r0, no Fried
            parameter), so the Term is the static mode-match loss only
            (_smf_static_term). It is then DETERMINISTIC but NOT mean-only, so the
            budget can still report a fade from the jitter walk-off Term.

    Returns:
        Term
            name="receive coupling (SMF)", category="coupling", mean_only=True
            (turbulence on) or a static deterministic Term (turbulence off).
    '''
    rx = scenario.rx_terminal
    detector = rx.detector
    if not isinstance(detector, SMF):
        raise ValueError(
            "terrestrial_smf_coupling_term needs an SMF detector on the far "
            "terminal. Set far = Terminal(..., detector=SMF())."
        )
    D = rx.aperture_m
    wavelength = rx.wavelength_m
    eta_max = _smf_eta_max(detector, D, wavelength)

    if not turbulence:
        return _smf_static_term(eta_max)

    tx = scenario.tx_terminal
    if tx.transmitter is None:
        raise ValueError(
            "terrestrial SMF coupling needs a launch beam: set the near terminal "
            "transmitter = Transmitter(waist_m=...)."
        )
    w0 = tx.transmitter.waist_m
    L = float(scenario.channel.path_length_m)
    cn2 = float(scenario.channel.cn2)

    # Horizontal Gaussian-beam Fried parameter over the constant-Cn2 path.
    hs = np.linspace(0.0, L, int(n_grid))
    cn2_profile = np.full_like(hs, cn2)
    r0 = gaussian_fried_parameter_profile(hs, cn2_profile, w0, wavelength,
                                          path='terrestrial')

    # When the receive tip-tilt walk-off Term is active, it carries the tip-tilt
    # coupling loss. So this coupling Term keeps the HIGHER-ORDER residual only.
    # A virtual TipTilt removes the tip-tilt (Noll modes 2 and 3). The
    # best-correcting stage wins, so a stack that already corrects more than
    # tip-tilt does not change. Source: Noll 1976 (residual variance by mode).
    stack = [*rx.compensation, TipTilt()] if drop_tiptilt else rx.compensation
    residual = apply_compensation(stack, D, r0)
    sigma2_res = residual.variance
    eta = _smf_coupling_efficiency(sigma2_res, eta_max)
    coupling_loss = -10.0 * np.log10(eta)
    dr0_eff = _effective_dr0(sigma2_res)

    assumptions = Assumptions(
        beam_type=BEAM_GAUSSIAN,
        turbulence_regime=REGIME_WEAK,
        spectrum=SPECTRUM_KOLMOGOROV,
        validity="Fidelity-0 MEAN-ONLY single-mode-fibre coupling for a horizontal "
                 "Gaussian beam. eta comes from the residual phase variance "
                 "(extended Marechal for a small residual, Dikmelik-Davidson "
                 "uncorrected coupling for a large one), evaluated at the "
                 "horizontal Gaussian-beam r0. This Term is DETERMINISTIC and "
                 "models NO fade.",
    )
    assumptions.flag(
        "Mean-only SMF coupling: this Term is the expected coupling loss and models "
        "NO fade (no sampler, no quantile). It locks the budget to fidelity 0, so "
        "no budget fade margin is reported."
    )
    assumptions.flag(
        "Effective-r0 weak-turbulence approximation: the Noll residual and the "
        "Dikmelik-Davidson coupling are plane-wave, Kolmogorov, phase-only forms "
        "evaluated at the Gaussian-beam r0. Valid only in weak turbulence; it "
        "ignores beam-wave amplitude scintillation, beam wander, and near-field "
        "curvature. A fidelity-2 split-step beam-propagation model is needed for "
        "those."
    )
    if rx.obscuration_ratio > 0.0:
        assumptions.flag(
            f"The far aperture has a central obscuration "
            f"(ratio={rx.obscuration_ratio:.3f}); the Dikmelik-Davidson coupling "
            "curve assumes a uniform circular aperture and does not model it."
        )
    dr0_max = float(np.max(dr0_eff))
    if dr0_max > SMF_DEEP_TURBULENCE_DR0:
        assumptions.flag(
            f"effective D/r0={dr0_max:.1f} exceeds {SMF_DEEP_TURBULENCE_DR0:.0f}; "
            "the practical uncorrected coupling curve is extrapolated."
        )

    base_shape = np.shape(coupling_loss)
    note = (f"terrestrial SMF coupling (mean-only), eta_max={eta_max:g}, "
            f"n_comp_modes={residual.n_modes}, r0={float(r0) * 100:.1f} cm")
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
            "beam": "gaussian-horizontal",
            "eta": float(eta) if base_shape == () else np.asarray(eta),
            "coupling_loss_db": float(coupling_loss) if base_shape == () else np.asarray(coupling_loss),
            "sigma2_res": float(sigma2_res) if np.ndim(sigma2_res) == 0 else np.asarray(sigma2_res),
            "effective_D_over_r0": float(dr0_eff) if np.ndim(dr0_eff) == 0 else np.asarray(dr0_eff),
            "r0_m": float(r0) if np.ndim(r0) == 0 else np.asarray(r0),
            "n_comp_modes": residual.n_modes,
        },
        assumptions=assumptions,
        mean_only=True,
    )


def smf_walkoff_term(scenario, geometry, *, n_grid=64, turbulence=True):
    '''
    Build the single-mode-fibre tip-tilt walk-off fade Term (terrestrial).

    This is the RECEIVE-side fibre walk-off. It is NOT the transmit pointing loss
    (olb.models.pointing, transmit beam against the far aperture). The received
    tip-tilt (beam wander + receive mechanical jitter) moves the focal spot on the
    fibre tip by f*theta. The focal spot radius w_s and the fibre mode radius w_m
    set the tolerated displacement together, through the effective scale
    w_eff=sqrt(w_s^2+w_m^2) (the two-Gaussian overlap). The fade is exponential in
    dB (see _walkoff_faces).

    The received tip-tilt uses _received_tiptilt_variance: a receive tip-tilt or
    AO stage tracks out the wander tilt, and the receive mechanical jitter always
    reaches the fibre. This Term carries a real fade (mean, quantile, sampler).
    With turbulence=False the beam-wander tilt drops, so the walk-off fade is the
    receive mechanical jitter alone.

    Note (single-mode-fibre subtlety): at a fixed coupling parameter a, the focal
    length f cancels in the tilt-to-coupling response, because w_m scales with f.
    This Term reads the PHYSICAL w_m and f, so it holds for the fixed optics that
    the terminal specifies.

    Parameters:
        scenario : TerrestrialScenario
            tx = near (its Transmitter launches the beam); rx = far (its SMF
            detector with focal_length_m and mode_field_radius_m, its
            pointing_jitter_rad, and its compensation stack).
        geometry : HorizontalPath
            Unused (the path length and Cn2 come from the channel). Kept for the
            f(scenario, geometry) -> Term signature.
        n_grid : int
            Points on the constant-Cn2 path grid for the wander integral.

    Returns:
        Term
            name="SMF tip-tilt walk-off", category="pointing". It has a real fade.
    '''
    rx = scenario.rx_terminal
    detector = rx.detector
    if not isinstance(detector, SMF):
        raise ValueError("smf_walkoff_term needs an SMF detector on the far terminal.")
    f, w_m = _smf_optics(detector, rx.aperture_m, rx.wavelength_m)
    if f is None or w_m is None:
        raise ValueError(
            "smf_walkoff_term needs the fibre-coupling optics to map a tip-tilt to "
            "a focal-plane displacement. Set SMF.focal_length_m and "
            "SMF.mode_field_radius_m, or set SMF.optimal_focus=True to derive the "
            "optimal focal length from the mode field radius."
        )
    # Effective focal-plane coupling scale for a lateral spot offset. The captured
    # power is the overlap integral of two Gaussians: the focal spot (1/e^2 radius
    # w_s) and the fibre LP01 mode (radius w_m). The offset dependence is
    # exp(-2*dx^2/(w_s^2+w_m^2)), so w_eff=sqrt(w_s^2+w_m^2). The spot radius is
    # w_s=lambda*f/(pi*(D/2)). A point spot (w_s->0) gives w_eff=w_m; that
    # underestimates the tolerated offset and gives too large a loss. Source:
    # two-Gaussian mode overlap, Shaklan and Roddier, Appl. Opt. 27 (1988) 2334,
    # DOI 10.1364/AO.27.002334.
    w_s = rx.wavelength_m * f / (np.pi * rx.aperture_m / 2.0)
    w_eff = np.sqrt(w_s ** 2 + w_m ** 2)
    sigma2_theta, meta = _received_tiptilt_variance(scenario, n_grid=n_grid,
                                                    turbulence=turbulence)
    mean, quantile, sampler = _walkoff_faces(f, w_eff, sigma2_theta)

    assumptions = Assumptions(
        beam_type=BEAM_GAUSSIAN,
        turbulence_regime=REGIME_WEAK,
        spectrum=SPECTRUM_KOLMOGOROV,
        validity="Receive tip-tilt walk-off of the focal spot against the "
                 "single-mode fibre mode. The wander arrival tilt is the weak "
                 "beam-wander model (Dios et al. 2004). A receive tip-tilt or AO "
                 "stage tracks it out (all-or-nothing, no bandwidth). The receive "
                 "mechanical jitter adds to it. The fade is exponential in dB. The "
                 "coupling loss has no upper limit: more tilt gives more loss. Only "
                 "the weak-turbulence validity of the tilt variance limits the Term.",
    )
    return Term(
        name="SMF tip-tilt walk-off",
        category="pointing",
        mean_db=float(mean),
        sampler=sampler,
        quantile=quantile,
        note=f"SMF walk-off, w_eff={w_eff * 1e6:.1f} um "
             f"(w_s={w_s * 1e6:.1f}, w_m={w_m * 1e6:.1f}), f={f * 1e3:.1f} mm, "
             f"wander_tracked={meta['wander_tracked']}",
        meta={**meta, "detector": "SMF", "focal_length_m": f,
              "mode_field_radius_m": w_m, "spot_radius_m": float(w_s),
              "w_eff_m": float(w_eff)},
        assumptions=assumptions,
    )


def mmf_coupling_term(scenario, geometry, *, n_grid=64, turbulence=True):
    '''
    Build the multimode-fibre coupling Term (terrestrial).

    A multimode fibre is a light bucket of core radius a_core in the fibre plane.
    The focal spot has the diffraction radius w_s = lambda*f / (pi*(D/2)). The
    coupling is the overlap of the focal-plane intensity with the core disk. The
    Term splits the loss into two parts:

      (i) STATIC spot-overfill loss. The on-axis Gaussian focal spot puts a
          fraction of its power inside the core:
              eta_static = 1 - exp(-2 a_core^2 / w_s^2).
          A spot smaller than the core gives eta_static ~ 1 (little loss). A spot
          larger than the core overfills it and loses power. Source: Andrews and
          Phillips, 2nd ed. (2005), DOI 10.1117/3.626196 (Gaussian power in a
          circular region; the Gaussian is the standard approximation to the
          Airy focal spot, with the 1/e^2 radius w_s = lambda*f/(pi*(D/2))).

      (ii) WALK-OFF fluctuation. The received tip-tilt moves the spot by f*theta.
           The spot walks off the core over the core radius scale. The fade is
           exponential in dB, with the core radius as the tolerated displacement
           (see _walkoff_faces). This carries a real fade.

    Anti-double-count: this Term is a MULTIPLICATIVE fibre-coupling efficiency on
    the aperture-collected field. It adds only -10log10(eta_mmf). It does NOT
    re-add the aperture capture (the geometric spreading Term already carries the
    free-space spread and the aperture capture). This is the same pattern the SMF
    Term uses. The receive mechanical jitter here uses the RECEIVE terminal, so it
    does not double-count the transmit pointing Term.

    With turbulence=False the static spot-overfill loss (i) stays (it is a fixed
    optical loss, not turbulence), and the walk-off (ii) keeps only the receive
    mechanical jitter (the beam-wander tilt drops).

    Parameters:
        scenario : TerrestrialScenario
            tx = near (its Transmitter launches the beam); rx = far (its MMF
            detector, aperture, pointing_jitter_rad, and compensation stack).
        geometry : HorizontalPath
            Unused (the path length and Cn2 come from the channel). Kept for the
            f(scenario, geometry) -> Term signature.
        n_grid : int
            Points on the constant-Cn2 path grid for the wander integral.

    Returns:
        Term
            name="receive coupling (MMF)", category="coupling". It has a real fade.
    '''
    rx = scenario.rx_terminal
    detector = rx.detector
    if not isinstance(detector, MMF):
        raise ValueError("mmf_coupling_term needs an MMF detector on the far terminal.")
    D = rx.aperture_m
    wavelength = rx.wavelength_m
    f = _mmf_focal_length(detector, D, wavelength)
    if f is None:
        raise ValueError(
            "mmf_coupling_term needs a focal length to map a tip-tilt to a "
            "focal-plane displacement. Set MMF.focal_length_m, or set "
            "MMF.optimal_focus=True to match the spot to the core."
        )
    a_core = detector.core_radius_m

    # Diffraction focal spot radius (1/e^2), Gaussian approximation to the Airy.
    w_s = wavelength * f / (np.pi * (D / 2.0))
    eta_static = 1.0 - np.exp(-2.0 * a_core ** 2 / w_s ** 2)
    static_db = -10.0 * np.log10(eta_static)

    sigma2_theta, meta = _received_tiptilt_variance(scenario, n_grid=n_grid,
                                                    turbulence=turbulence)
    walk_mean, walk_quantile, walk_sampler = _walkoff_faces(f, a_core, sigma2_theta)

    mean_db = static_db + walk_mean

    def quantile(p):
        return static_db + walk_quantile(p)

    def sampler(n, rng):
        return static_db + walk_sampler(n, rng)

    assumptions = Assumptions(
        beam_type=BEAM_GAUSSIAN,
        turbulence_regime=REGIME_WEAK,
        spectrum=SPECTRUM_KOLMOGOROV,
        validity="Multimode-fibre coupling: a Gaussian focal spot (radius w_s = "
                 "lambda*f/(pi*(D/2))) overlaps a core disk of radius a_core. The "
                 "static loss is the on-axis spot power outside the core. The "
                 "walk-off is the receive tip-tilt (weak beam wander, tracked by a "
                 "tip-tilt or AO stage; plus the receive mechanical jitter). The "
                 "spot model assumes a uniform, unobscured circular aperture.",
    )
    if rx.obscuration_ratio > 0.0:
        assumptions.flag(
            f"The far aperture has a central obscuration "
            f"(ratio={rx.obscuration_ratio:.3f}); the Gaussian focal-spot model "
            "assumes a uniform circular aperture and does not model it."
        )
    return Term(
        name="receive coupling (MMF)",
        category="coupling",
        mean_db=float(mean_db),
        sampler=sampler,
        quantile=quantile,
        note=f"MMF coupling, a_core={a_core * 1e6:.1f} um, "
             f"w_s={w_s * 1e6:.1f} um, static={static_db:.2f} dB",
        meta={**meta, "detector": "MMF", "core_radius_m": a_core,
              "focal_length_m": f, "spot_radius_m": float(w_s),
              "eta_static": float(eta_static), "static_loss_db": float(static_db),
              "walkoff_mean_db": float(walk_mean)},
        assumptions=assumptions,
    )


if __name__ == '__main__':
    import warnings

    from ..scenario import SpaceScenario, Channel
    from ..geometry import CircularOrbit
    from ..terminal import Terminal, Transmitter, TipTilt, AO

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

    # The SMF tests use the analytic mean-only fidelity (the "fast" default needs
    # fast-aosim, which is optional). The Aperture detector ignores smf_fidelity.
    def build(scn):
        return rx_coupling_term(scn, geom, hs=hs, smf_fidelity="mean",
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

    # --- SMF, no AO: mean-only Marechal / Dikmelik coupling loss ------------
    # A 0.7 m fibre receiver with no correction. Deterministic: mean loss only.
    scn_smf = _downlink(Terminal(aperture_m=0.7, wavelength_m=lam, detector=SMF()))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        t_smf = build(scn_smf)
    assert t_smf.name == "receive coupling (SMF)"
    assert t_smf.meta["model"] == "mean-only"
    assert t_smf.quantile is None and not t_smf.stochastic   # deterministic
    assert t_smf.meta["coupling_loss_db"] > 0.0
    # Mean-only caveat is always flagged, so it is never ok.
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
    # A central obscuration on the receive aperture flags a violation naming the
    # coupling curve; an unobscured aperture does not.
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
    # Every SMF mean-only term flags that it carries no fade, so it is never ok.
    rng = np.random.default_rng(0)
    for t in (t_smf, t_ao, t_obsc):
        assert any("Mean-only" in v for v in t.assumptions.violations), \
            t.assumptions.violations
        assert not t.assumptions.ok
        # Deterministic: the quantile is the constant mean, and samples broadcast it.
        assert not t.stochastic and t.quantile is None
        assert np.isclose(t.quantile_db(0.99), t.mean_db)
        draws = t.sample_db(1000, rng)
        assert draws.shape == (1000,) and np.allclose(draws, t.mean_db)

    # An elevation sweep broadcasts the mean-only term.
    sweep = CircularOrbit(600e3, np.array([40.0, 60.0, 90.0]))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        t_sweep = rx_coupling_term(scn_ao, sweep, hs=hs, smf_fidelity="mean",
                                   cn2_profile=default_cn2_profile(scn_ao.channel.site, hs))
    assert np.shape(t_sweep.mean_db) == (3,)
    assert t_sweep.sample_db(100, rng).shape == (100, 3)

    # An unknown fidelity is refused.
    try:
        rx_coupling_term(scn_smf, geom, hs=hs, smf_fidelity="reciprocity",
                         cn2_profile=default_cn2_profile(scn_smf.channel.site, hs))
        raise AssertionError("unknown smf_fidelity must raise")
    except ValueError as e:
        assert "smf_fidelity" in str(e)

    print(f"aperture coupling mean = {float(t_ap.mean_db):.4f} dB "
          f"(= scintillation {float(t_scint.mean_db):.4f} dB)")
    print(f"SMF no-AO (mean-only):   coupling loss {float(t_smf.mean_db):.2f} dB  "
          f"eta={t_smf.meta['eta']:.4f}")
    print(f"SMF +AO200 (mean-only):  eta={t_ao.meta['eta']:.4f}  "
          f"coupling loss={t_ao.meta['coupling_loss_db']:.2f} dB  "
          f"sigma2_res={t_ao.meta['sigma2_res']:.4f}")

    # --- Step 2: eta_max from the coupling parameter a ---------------------
    # The overlap curve peaks at 0.8145 near a=1.12 and falls on both sides.
    assert np.isclose(smf_eta_max_from_a(1.12), 0.8145, atol=1e-3), smf_eta_max_from_a(1.12)
    assert smf_eta_max_from_a(0.5) < 0.8145 and smf_eta_max_from_a(2.5) < 0.8145
    a_peak = np.linspace(0.9, 1.4, 51)
    assert np.isclose(a_peak[np.argmax(smf_eta_max_from_a(a_peak))], 1.12, atol=0.03)

    # f=None keeps the eta_max field exactly (today's behaviour).
    assert _smf_eta_max(SMF(), 0.2, lam) == 0.8145
    # Set f and w_m to sit on the optimum a=1.12: eta_max ~ 0.8145.
    D_test, wm = 0.2, 5.2e-6
    f_opt = np.pi * (D_test / 2.0) * wm / (lam * 1.12)   # a = 1.12
    assert np.isclose(_smf_eta_max(SMF(focal_length_m=f_opt, mode_field_radius_m=wm),
                                   D_test, lam), 0.8145, atol=1e-3)
    # An off-optimum focal length (a far from 1.12) gives a lower eta_max.
    eta_off = _smf_eta_max(SMF(focal_length_m=4.0 * f_opt, mode_field_radius_m=wm),
                           D_test, lam)
    assert eta_off < 0.8145, eta_off
    # A focal length with no mode field radius is refused.
    try:
        _smf_eta_max(SMF(focal_length_m=0.02), D_test, lam)
        raise AssertionError("focal length without mode field radius must raise")
    except ValueError as e:
        assert "mode_field_radius_m" in str(e)

    # --- Steps 3 and 4: received tip-tilt, SMF walk-off, MMF coupling ------
    from ..scenario import TerrestrialScenario, TerrestrialChannel
    from ..geometry import HorizontalPath
    from ..terminal import Aperture

    def _terr(detector, *, jitter=0.0, compensation=None, cn2=1e-14,
              far_aperture=0.2, w0=0.02, L=3e3):
        return TerrestrialScenario(
            near=Terminal(aperture_m=0.3, wavelength_m=lam,
                          transmitter=Transmitter(waist_m=w0, power_dbm=30)),
            far=Terminal(aperture_m=far_aperture, wavelength_m=lam,
                         pointing_jitter_rad=jitter, detector=detector,
                         compensation=compensation or []),
            channel=TerrestrialChannel(path_length_m=L, attenuation_db_per_km=0.5,
                                       cn2=cn2))

    hpath = HorizontalPath(3e3)

    # Received tip-tilt: a tip-tilt stage tracks out the wander, jitter remains.
    v_full, _ = _received_tiptilt_variance(_terr(Aperture(), jitter=5e-6), n_grid=64)
    v_track, m_track = _received_tiptilt_variance(
        _terr(Aperture(), jitter=5e-6, compensation=[TipTilt()]), n_grid=64)
    assert v_full > v_track, (v_full, v_track)         # tracking removes the wander
    assert m_track["wander_tracked"] and m_track["sigma2_wander"] == 0.0
    assert np.isclose(v_track, 2.0 * 5e-6 ** 2)        # only the jitter remains

    # SMF walk-off: a real fade. A larger jitter gives a deeper walk-off loss.
    smf_opt = SMF(focal_length_m=f_opt, mode_field_radius_m=wm, sensitivity_dbm=-40)
    wo_small = smf_walkoff_term(_terr(smf_opt, jitter=2e-6), hpath)
    wo_big = smf_walkoff_term(_terr(smf_opt, jitter=10e-6), hpath)
    assert wo_small.stochastic and wo_small.quantile is not None
    assert wo_big.mean_db > wo_small.mean_db
    assert wo_big.quantile_db(0.99) > wo_big.mean_db   # exponential tail
    rng = np.random.default_rng(0)
    draws = wo_big.sample_db(50_000, rng)
    assert np.abs(draws.mean() - wo_big.mean_db) / wo_big.mean_db < 0.03
    # optimal_focus: derive f from a=1.12 (the far aperture is D_test=0.2). The
    # derived f equals the explicit optimal f_opt, so the walk-off matches, and
    # eta_max is the 0.8145 peak. It also works from the MFD alone (SMF-28 default).
    smf_focus = SMF(mode_field_radius_m=wm, optimal_focus=True, sensitivity_dbm=-40)
    wo_focus = smf_walkoff_term(_terr(smf_focus, jitter=10e-6), hpath)
    assert np.isclose(wo_focus.meta["focal_length_m"], f_opt)
    assert np.isclose(wo_focus.mean_db, wo_big.mean_db)
    assert np.isclose(_smf_eta_max(smf_focus, D_test, lam), 0.8145, atol=1e-4)
    wo_default = smf_walkoff_term(_terr(SMF(optimal_focus=True), jitter=10e-6), hpath)
    assert np.isclose(wo_default.mean_db, wo_big.mean_db)
    # An SMF without the optics cannot map a tip-tilt to a displacement.
    try:
        smf_walkoff_term(_terr(SMF()), hpath)
        raise AssertionError("SMF walk-off without optics must raise")
    except ValueError:
        pass

    # MMF coupling: static overfill + walk-off, with a real fade.
    mmf = MMF(core_radius_m=25e-6, focal_length_m=0.05, sensitivity_dbm=-38)
    t_mmf = mmf_coupling_term(_terr(mmf, jitter=5e-6), hpath)
    assert t_mmf.name == "receive coupling (MMF)" and t_mmf.category == "coupling"
    assert t_mmf.stochastic and t_mmf.quantile is not None and not t_mmf.mean_only
    assert 0.0 < t_mmf.meta["eta_static"] <= 1.0
    assert t_mmf.meta["static_loss_db"] >= 0.0
    assert t_mmf.quantile_db(0.99) > t_mmf.mean_db     # walk-off adds a fade
    # A smaller core overfills more (larger static loss) and walks off sooner.
    t_small_core = mmf_coupling_term(
        _terr(MMF(core_radius_m=5e-6, focal_length_m=0.05), jitter=5e-6), hpath)
    t_big_core = mmf_coupling_term(
        _terr(MMF(core_radius_m=50e-6, focal_length_m=0.05), jitter=5e-6), hpath)
    assert t_small_core.mean_db > t_big_core.mean_db, (t_small_core.mean_db,
                                                       t_big_core.mean_db)
    # A larger jitter deepens the MMF walk-off loss.
    t_mmf_calm = mmf_coupling_term(_terr(mmf, jitter=1e-6), hpath)
    assert t_mmf.mean_db > t_mmf_calm.mean_db
    # MMF optimal_focus: derive f to match the spot to the core (a_core/w_s=1.12).
    # It matches an explicit derived focal length, and gives about 92% static
    # capture. An MMF with no focal length and no optimal_focus is refused.
    a_core = 25e-6
    f_mmf = np.pi * (D_test / 2.0) * a_core / (lam * 1.12)
    mmf_focus = MMF(core_radius_m=a_core, optimal_focus=True, sensitivity_dbm=-38)
    t_focus = mmf_coupling_term(_terr(mmf_focus, jitter=5e-6), hpath)
    t_explicit = mmf_coupling_term(
        _terr(MMF(core_radius_m=a_core, focal_length_m=f_mmf), jitter=5e-6), hpath)
    assert np.isclose(t_focus.meta["focal_length_m"], f_mmf)
    assert np.isclose(t_focus.mean_db, t_explicit.mean_db)
    assert np.isclose(t_focus.meta["eta_static"], 1.0 - np.exp(-2.0 * 1.12 ** 2), atol=1e-3)
    try:
        mmf_coupling_term(_terr(MMF(core_radius_m=a_core)), hpath)
        raise AssertionError("MMF without a focal length must raise")
    except ValueError:
        pass
    # --- turbulence=False: static coupling + jitter, no turbulence quantity --
    # Aperture: 0 dB (no scintillation). SMF: static mode-match loss only, and
    # NOT mean-only, so it does not lock a budget out of a jitter fade.
    ap_off = rx_coupling_term(scn_ap, geom, turbulence=False)
    assert ap_off.mean_db == 0.0 and ap_off.meta["model"] == "static"
    smf_off = rx_coupling_term(scn_smf, geom, turbulence=False)
    assert smf_off.meta["model"] == "static" and not smf_off.mean_only
    assert not smf_off.stochastic and smf_off.quantile is None
    assert np.isclose(smf_off.mean_db, -10.0 * np.log10(0.8145))
    # The received tip-tilt with turbulence off is the mechanical jitter alone.
    v_off, m_off = _received_tiptilt_variance(_terr(Aperture(), jitter=5e-6),
                                              n_grid=64, turbulence=False)
    assert m_off["sigma2_wander"] == 0.0 and np.isclose(v_off, 2.0 * 5e-6 ** 2)
    # The SMF walk-off then carries the jitter fade with NO turbulence.
    wo_off = smf_walkoff_term(_terr(smf_opt, jitter=10e-6), hpath, turbulence=False)
    assert wo_off.meta["sigma2_wander"] == 0.0 and wo_off.meta["sigma2_jitter"] > 0.0
    assert wo_off.quantile_db(0.99) > wo_off.mean_db
    # An SMF static term needs no launch beam (turbulence off skips the r0 path).
    terr_off = terrestrial_smf_coupling_term(
        TerrestrialScenario(
            near=Terminal(aperture_m=0.3, wavelength_m=lam),   # no transmitter
            far=Terminal(aperture_m=0.2, wavelength_m=lam, detector=SMF()),
            channel=TerrestrialChannel(path_length_m=3e3, attenuation_db_per_km=0.5,
                                       cn2=1e-14)),
        hpath, turbulence=False)
    assert terr_off.meta["model"] == "static" and not terr_off.mean_only

    print(f"SMF walk-off (10 urad): mean {wo_big.mean_db:.3f} dB  "
          f"99% {float(wo_big.quantile_db(0.99)):.3f} dB")
    print(f"MMF (25 um core, 5 urad): static {t_mmf.meta['static_loss_db']:.3f} dB  "
          f"mean {t_mmf.mean_db:.3f} dB  99% {float(t_mmf.quantile_db(0.99)):.3f} dB")
    print("self-check passed")
