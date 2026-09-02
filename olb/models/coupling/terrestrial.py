'''
Receive-coupling Terms for a terrestrial (horizontal) link.

This module builds the receive-coupling Terms of a terrestrial far terminal: the
single-mode-fibre mean coupling, the single-mode-fibre tip-tilt walk-off fade, and
the multimode-fibre (light-bucket) coupling. The shared single-mode-fibre coupling
physics lives in olb.models.coupling._common; this module reads it. The received
tip-tilt (beam wander plus receive mechanical jitter) lives here, because only the
terrestrial Terms use it.

RECEIVED CURVATURE (the olb convention). A terrestrial received beam is a
DIVERGING Gaussian, not a plane wave. So the true focus behind a coupling lens of
focal length f is at z = f + dz_curv, with dz_curv = f^2/(R_rx - f). Every Term
here ALWAYS charges that curvature defocus, at the ACTUAL fibre plane: the
detector sits at z = f + defocus_m, so its distance from the true focus is
dz_eff = defocus_m - dz_curv. The optimal_focus flag keeps its meaning (a
focal-LENGTH rule) and never moves the detector. Use curvature_focus_shift to get
dz_curv and set detector.defocus_m to it for a coupler aligned on the true focus.

Sources:
  Marechal / Dikmelik-Davidson coupling: see olb.models.coupling._common.
  eta_max(a): Shaklan and Roddier, Appl. Opt. 27, 2334 (1988), DOI
  10.1364/AO.27.002334.
  Beam-wander arrival tilt: Dios et al. 2004 (see olb.turbulence.angle_of_arrival).
  Off-axis Gaussian encircled energy (Marcum Q): Marcum, RAND RM-753 (1950).
  Thin-lens focus shift of a spherical Gaussian input: S. A. Self, "Focusing of
  spherical Gaussian beams," Appl. Opt. 22, 658 (1983), DOI 10.1364/AO.22.000658.
  Aberrated single-mode coupling: Ruilier and Cassaing, JOSA A 18, 143 (2001),
  DOI 10.1364/JOSAA.18.000143.
'''

import numpy as np

from ...results import Term
from ...assumptions import (trace_assumptions, BEAM_GAUSSIAN, REGIME_WEAK,
                            SPECTRUM_KOLMOGOROV)
from ...terminal import SMF, MMF, TipTilt, AO
from ...beam import (free_space_radius, launch_curvature, gaussz,
                     phase_front_radius)
from ...turbulence.gaussian_fried import gaussian_fried_parameter_profile
from ...turbulence.angle_of_arrival import wander_arrival_angle_variance
from ...turbulence.ao import apply_compensation
from ...turbulence.andrews.beam import beam_params
from ...turbulence.andrews.scintillation import (rytov_weak, rytov_variance,
                                                 WEAK_REGIME_LIMIT)
from ._common import (_smf_eta_max, _smf_coupling_efficiency, _effective_dr0,
                     _smf_static_term, smf_eta_defocused, SMF_OPTIMAL_A,
                     SMF_DEEP_TURBULENCE_DR0)

# dB per (r^2 / w^2), from the exponential Gaussian power falloff. Same constant
# as olb.models.pointing (20/ln10). Source: Andrews and Phillips, 2nd ed. (2005),
# DOI 10.1117/3.626196 (Gaussian power in a circular region).
_K = 20.0 / np.log(10.0)

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
        return np.pi * (D / 2.0) * detector.core_radius_m / (wavelength * SMF_OPTIMAL_A)
    return None


# --- the received-curvature focus shift --------------------------------------

def _received_curvature(scenario, f):
    '''
    Return the received phase-front radius and the focus shift it causes.

    A terrestrial received beam is a DIVERGING Gaussian, not a plane wave. Its
    phase-front radius at the receive aperture is R_rx (see
    olb.beam.phase_front_radius; Andrews and Phillips, 2nd ed. (2005),
    DOI 10.1117/3.626196, Ch. 4, Eqs. (7) and (8)). A thin lens of focal length f
    images that diverging input BEYOND its focal plane, at

        z = f + dz_curv,        dz_curv = f^2 / (R_rx - f).

    Source: S. A. Self, "Focusing of spherical Gaussian beams," Appl. Opt. 22,
    658 (1983), DOI 10.1364/AO.22.000658. For a space link R_rx is enormous, so
    dz_curv is about zero; for a short horizontal link it is millimetres, which
    is many focal Rayleigh ranges.

    The olb CONVENTION: the curvature defocus is ALWAYS charged, at the ACTUAL
    fibre plane. The detector sits at z = f + defocus_m, so its distance from the
    TRUE focus is

        dz_eff = defocus_m - dz_curv.

    The optimal_focus flag is a focal-LENGTH rule; it never moves the detector.
    A user who wants a tracked (aligned-at-true-focus) coupler sets
    detector.defocus_m = curvature_focus_shift(scenario).

    A scenario with no launch beam (no Transmitter) gives (inf, 0.0): the model
    cannot know the received curvature, so it charges nothing and the caller
    flags it.

    Parameters:
        scenario : TerrestrialScenario
        f : float
            The resolved focal length of the receive coupling optic [m].

    Returns:
        tuple
            (R_rx, dz_curv) in m. R_rx is positive for a diverging front.
    '''
    tx = scenario.tx_terminal
    if tx.transmitter is None:
        return float('inf'), 0.0
    R_rx = float(phase_front_radius(tx.transmitter.waist_m,
                                    float(scenario.channel.path_length_m),
                                    tx.transmitter.divergence_rad,
                                    scenario.rx_terminal.wavelength_m))
    return R_rx, float(f ** 2 / (R_rx - f))


def curvature_focus_shift(scenario):
    '''
    Return the received-curvature focus shift dz_curv of the RECEIVE optics [m].

    This is the PUBLIC helper. The true focus of the received diverging beam sits
    at z = f + dz_curv, not at z = f (see _received_curvature). So a coupler that
    is aligned on the true focus has

        detector.defocus_m = curvature_focus_shift(scenario).

    Set that value on the detector to model a TRACKED (aligned) coupler. Leave
    defocus_m at 0.0 to model a fibre at the nominal focal plane, which pays the
    full curvature defocus.

    Parameters:
        scenario : TerrestrialScenario
            Its rx terminal must carry an SMF or MMF detector with resolvable
            coupling optics.

    Returns:
        float
            dz_curv [m]. Positive: the true focus is BEYOND the focal plane.

    Raises:
        ValueError
            The receive detector is not a fibre, or its focal length cannot be
            resolved.
    '''
    rx = scenario.rx_terminal
    detector = rx.detector
    if isinstance(detector, SMF):
        f, _ = _smf_optics(detector, rx.aperture_m, rx.wavelength_m)
    elif isinstance(detector, MMF):
        f = _mmf_focal_length(detector, rx.aperture_m, rx.wavelength_m)
    else:
        raise ValueError(
            "curvature_focus_shift needs an SMF or MMF detector on the receive "
            "terminal.")
    if f is None:
        raise ValueError(
            "curvature_focus_shift needs the coupling focal length. Set "
            "focal_length_m, or set optimal_focus=True.")
    return _received_curvature(scenario, f)[1]


def _flag_curvature(assumptions, R_rx, f, dz_curv, defocus_m):
    '''
    Add the received-curvature assumption flags to a Term.

    Two facts the physics functions never see:
      - the thin-lens focus shift is a SMALL-shift geometry. It needs the input
        phase front well outside the focal length. R_rx <= 2*f breaks that (the
        image runs away), so flag it but still compute.
      - a scenario with no launch beam gives no received curvature, so the model
        charges NONE. Say so.
    '''
    if not np.isfinite(R_rx):
        assumptions.flag(
            "No launch beam on the transmit terminal, so the received phase-front "
            "curvature is unknown and NO curvature defocus is charged. The "
            "coupling loss is therefore OPTIMISTIC for a short horizontal link.",
            source="factory:models.coupling.terrestrial",
        )
        return
    if R_rx <= 2.0 * f:
        assumptions.flag(
            f"The received phase-front radius R_rx={R_rx:.3g} m is not large "
            f"against the focal length f={f:.3g} m (R_rx <= 2f). The thin-lens "
            "focus shift dz_curv = f^2/(R_rx - f) (S. A. Self, Appl. Opt. 22, 658 "
            "(1983), DOI 10.1364/AO.22.000658) is then no longer a small shift, "
            "so the image geometry of this model breaks. The value is still "
            "computed.",
            source="factory:models.coupling.terrestrial",
        )


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


def _spot_offset_sigma(f, dz, sigma2_theta_radial):
    '''
    Return the per-axis 1-sigma spot offset at the detector plane [m].

    The detector sits at z = f + dz, with dz the defocus distance (dz=0 at focus).
    A received arrival tilt theta moves the spot centre off the optical axis by the
    ray-optics chief-ray lever of a thin lens:
        d_spot = (f+dz)*theta,
    where theta is the received arrival tilt (radial, 2-axis). At the focus (dz=0)
    the lever is f, the focal-plane displacement; off focus the longer lever arm
    (f+dz) moves the spot more. theta has two independent Gaussian axes, so the
    per-axis variance is:
        sigma_d^2 = (f+dz)^2 * (sigma2_theta_radial/2).
    sigma2_theta_radial is the RADIAL (2-axis) tilt variance. Source: geometric
    (ray-optics chief-ray) of a thin lens; Gaussian relations, Andrews and
    Phillips, 2nd ed. (2005), DOI 10.1117/3.626196, Ch. 4.
    '''
    return np.sqrt((f + dz) ** 2 * float(sigma2_theta_radial) / 2.0)


def _walkoff_faces(w_eff, sigma_d):
    '''
    Return the (mean_db, quantile, sampler) of the walk-off fade at the detector.

    A received tip-tilt moves the spot on the detector by a per-axis 1-sigma
    sigma_d (see _spot_offset_sigma). The radial
    displacement dx has two i.i.d. Gaussian axes, so dx^2 is exponential. The
    captured-power fraction against the coupling feature of size w_eff is
        h(dx) = exp(-2 dx^2 / w_eff^2)   ->   loss_db = _K * dx^2 / w_eff^2.
    So the loss in dB is exponential:
        loss_db ~ Exponential(mean = _K * 2 * sigma_d^2 / w_eff^2).
    This is the same closed form as olb.models.pointing. w_eff is the effective
    coupling scale that the caller gives. For a single-mode fibre it is
    sqrt(w_det^2+w_m^2), the two-Gaussian overlap of the (defocused) spot (radius
    w_det) and the fibre mode (radius w_m). At focus w_det=w_s, the diffraction
    spot radius. Source: Andrews and Phillips, 2nd ed. (2005),
    DOI 10.1117/3.626196 (Gaussian power falloff and 2-D Gaussian jitter);
    two-Gaussian mode overlap, Shaklan and Roddier, Appl. Opt. 27 (1988) 2334,
    DOI 10.1364/AO.27.002334.
    '''
    # Radial displacement variance E[dx^2] = 2*sigma_d^2 (two i.i.d. axes).
    mean = _K * 2.0 * float(sigma_d) ** 2 / w_eff ** 2

    def quantile(p):
        return -mean * np.log(1.0 - p)      # inverse exponential CDF

    def sampler(n, rng):
        if mean <= 0.0:
            return np.zeros(n)
        return rng.exponential(scale=mean, size=n)

    return mean, quantile, sampler


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
            The tx terminal's Transmitter waist launches the beam; the rx
            terminal supplies the SMF detector, aperture, obscuration, and
            compensation stack. The scenario `direction` maps the roles
            (forward: tx = near, rx = far; reverse swaps them).
        geometry : HorizontalPath
            Unused here (path length and Cn2 come from the channel). Kept for the
            f(scenario, geometry) -> Term signature.
        n_grid : int
            Points on the constant-Cn2 path grid for the r0 integral.
        drop_tiptilt : bool
            Remove the tip-tilt (Noll modes 2 and 3) from the residual phase
            variance. Set True when the budget also adds the receive tip-tilt
            walk-off Term (terrestrial_smf_walkoff_term). The walk-off Term then
            carries the tip-tilt coupling loss, and this Term keeps the
            HIGHER-ORDER residual only. So the tip-tilt is not counted two times.
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

    # The flat-wavefront mode match, corrected for the RECEIVED CURVATURE. The
    # received beam is a diverging Gaussian, so its true focus is at
    # z = f + dz_curv, not at z = f (see _received_curvature). The fibre sits at
    # z = f + defocus_m, so it is dz_eff = defocus_m - dz_curv from the true
    # focus, and the pupil then carries a quadratic (defocus) phase of edge
    # coefficient c. The closed form smf_eta_defocused(a, c) holds that overlap,
    # and it reduces to eta_max(a) at c = 0 (Shaklan and Roddier, DOI
    # 10.1364/AO.27.002334; Ruilier and Cassaing, DOI 10.1364/JOSAA.18.000143).
    # Without resolvable optics (no focal length) there is no a and no c, so keep
    # the plain eta_max field and flag that the curvature is NOT modelled.
    f_rx, w_m_rx = _smf_optics(detector, D, wavelength)
    curvature_flags = []
    if f_rx is not None and w_m_rx is not None:
        a_rx = np.pi * (D / 2.0) * w_m_rx / (wavelength * f_rx)
        R_rx, dz_curv = _received_curvature(scenario, f_rx)
        dz_eff = float(detector.defocus_m) - dz_curv
        c_rx = np.pi * dz_eff * (D / 2.0) ** 2 / (wavelength * f_rx ** 2)
        eta_flat = float(smf_eta_defocused(a_rx, 0.0))
        eta_max = float(smf_eta_defocused(a_rx, c_rx))
        curvature_flags.append(('curvature', R_rx, f_rx, dz_curv,
                                float(detector.defocus_m)))
    else:
        a_rx = R_rx = dz_curv = dz_eff = c_rx = None
        eta_max = _smf_eta_max(detector, D, wavelength)
        eta_flat = eta_max
        curvature_flags.append(('no-optics', None, None, None, None))
    curvature_meta = {
        "received_curvature_m": R_rx,
        "curvature_defocus_m": dz_curv,
        "defocus_m": float(detector.defocus_m),
        "dz_eff_m": dz_eff,
        "defocus_coefficient_rad": c_rx,
        "eta_flat_wavefront": eta_flat,
        "eta_max": float(eta_max),
        "curvature_penalty_db": (float(-10.0 * np.log10(eta_max / eta_flat))
                                 if eta_flat > 0.0 else None),
    }

    def _add_curvature_flags(assumptions):
        '''Attach the curvature flags of this Term (shared by both branches).'''
        kind = curvature_flags[0][0]
        if kind == 'no-optics':
            assumptions.flag(
                "The single-mode-fibre detector has no resolvable coupling optics "
                "(no focal length), so the RECEIVED-CURVATURE defocus penalty is "
                "NOT modelled here: eta_max is the flat-wavefront mode match. Set "
                "SMF.focal_length_m and SMF.mode_field_radius_m, or "
                "SMF.optimal_focus=True, to charge it.",
                source="factory:models.coupling.terrestrial",
            )
        else:
            _flag_curvature(assumptions, R_rx, f_rx, dz_curv,
                            float(detector.defocus_m))

    if not turbulence:
        # Curvature is STATIC optics, not turbulence, so the turbulence-off Term
        # charges it too. Pass the already-aberrated eta to the shared static Term
        # (its signature stays as the downlink uses it), then add the meta and the
        # curvature flags here.
        term = _smf_static_term(eta_max)
        term.meta.update(curvature_meta)
        _add_curvature_flags(term.assumptions)
        return term

    tx = scenario.tx_terminal
    if tx.transmitter is None:
        raise ValueError(
            "terrestrial SMF coupling needs a launch beam: set the near terminal "
            "transmitter = Transmitter(waist_m=...)."
        )
    w0 = tx.transmitter.waist_m
    L = float(scenario.channel.path_length_m)
    cn2 = float(scenario.channel.cn2)

    # Horizontal Gaussian-beam Fried parameter over the constant-Cn2 path. The
    # launch curvature f0 of a deliberately diverged beam enters the beam
    # parameters (Theta0 = 1 - L/f0), so a diverged beam gets its own r0.
    # This closes olb Gap 3. See olb.beam.launch_curvature.
    f0 = launch_curvature(w0, tx.transmitter.divergence_rad, wavelength)
    hs = np.linspace(0.0, L, int(n_grid))
    cn2_profile = np.full_like(hs, cn2)

    # Open the collection context around the PHYSICS CALLS only. The Gaussian
    # Fried parameter and the compensation chain register their own assumptions
    # (beam type, weak regime, spectrum, launch curvature, path weight, and the
    # extended-Marechal residual check), so the Term inherits the union. A strong
    # path trips the Marechal residual check inside apply_compensation
    # automatically.
    with trace_assumptions() as trace:
        r0 = gaussian_fried_parameter_profile(hs, cn2_profile, w0, wavelength,
                                              path='terrestrial', f0=f0)
        # When the receive tip-tilt walk-off Term is active, it carries the
        # tip-tilt coupling loss. So this coupling Term keeps the HIGHER-ORDER
        # residual only. A virtual TipTilt removes the tip-tilt (Noll modes 2 and
        # 3). The best-correcting stage wins, so a stack that already corrects
        # more than tip-tilt does not change. Source: Noll 1976 (residual
        # variance by mode).
        stack = [*rx.compensation, TipTilt()] if drop_tiptilt else rx.compensation
        residual = apply_compensation(stack, D, r0)
    sigma2_res = residual.variance
    eta = _smf_coupling_efficiency(sigma2_res, eta_max)
    coupling_loss = -10.0 * np.log10(eta)
    dr0_eff = _effective_dr0(sigma2_res)

    # The traced physics own the beam type, the weak regime, the spectrum, and
    # the launch-curvature / path-weight / Marechal constraints; the merge
    # inherits their union and any traced violation. State the three headline
    # fields explicitly (this is a Gaussian-beam mean-only coupling Term).
    assumptions = trace.merge(
        beam_type=BEAM_GAUSSIAN,
        turbulence_regime=REGIME_WEAK,
        spectrum=SPECTRUM_KOLMOGOROV,
        validity="Fidelity-0 MEAN-ONLY single-mode-fibre coupling for a horizontal "
                 "Gaussian beam. eta comes from the residual phase variance "
                 "(extended Marechal for a small residual, Dikmelik-Davidson "
                 "uncorrected coupling for a large one), evaluated at the "
                 "horizontal Gaussian-beam r0. The flat-wavefront mode match is "
                 "replaced by the DEFOCUS-ABERRATED closed form "
                 "eta(a, c) = 2 a^2 |(1-exp(-(a^2-ic)))/(a^2-ic)|^2 (Shaklan and "
                 "Roddier, DOI 10.1364/AO.27.002334; Ruilier and Cassaing, "
                 "DOI 10.1364/JOSAA.18.000143), with "
                 "c = pi*dz_eff*(D/2)^2/(lambda*f^2). The received curvature "
                 "defocus is ALWAYS charged at the ACTUAL fibre plane: "
                 "dz_eff = defocus_m - dz_curv, dz_curv = f^2/(R_rx - f) "
                 "(S. A. Self, Appl. Opt. 22, 658 (1983), "
                 "DOI 10.1364/AO.22.000658), R_rx the received phase-front radius "
                 "(Andrews and Phillips 2005, Ch. 4, DOI 10.1117/3.626196). "
                 "optimal_focus is a focal-LENGTH rule; it never moves the "
                 "detector. This Term is DETERMINISTIC and models NO fade.",
    )
    _add_curvature_flags(assumptions)
    # Scenario-level facts the physics never sees stay factory flags, source
    # tagged so they carry the same "[source] reason" prefix as a traced check.
    assumptions.flag(
        "Mean-only SMF coupling: this Term is the expected coupling loss and models "
        "NO fade (no sampler, no quantile). It locks the budget to fidelity 0, so "
        "no budget fade margin is reported.",
        source="factory:models.coupling.terrestrial",
    )
    assumptions.flag(
        "Effective-r0 weak-turbulence approximation: the Noll residual and the "
        "Dikmelik-Davidson coupling are plane-wave, Kolmogorov, phase-only forms "
        "evaluated at the Gaussian-beam r0. Valid only in weak turbulence; it "
        "ignores beam-wave amplitude scintillation, beam wander, and near-field "
        "curvature. A fidelity-2 split-step beam-propagation model is needed for "
        "those.",
        source="factory:models.coupling.terrestrial",
    )
    if rx.obscuration_ratio > 0.0:
        assumptions.flag(
            f"The far aperture has a central obscuration "
            f"(ratio={rx.obscuration_ratio:.3f}); the Dikmelik-Davidson coupling "
            "curve assumes a uniform circular aperture and does not model it.",
            source="factory:models.coupling.terrestrial",
        )
    dr0_max = float(np.max(dr0_eff))
    if dr0_max > SMF_DEEP_TURBULENCE_DR0:
        assumptions.flag(
            f"effective D/r0={dr0_max:.1f} exceeds {SMF_DEEP_TURBULENCE_DR0:.0f}; "
            "the practical uncorrected coupling curve is extrapolated.",
            source="factory:models.coupling.terrestrial",
        )

    base_shape = np.shape(coupling_loss)
    note = (f"terrestrial SMF coupling (mean-only), eta_max={eta_max:g}, "
            f"n_comp_modes={residual.n_modes}, r0={float(r0) * 100:.1f} cm"
            + (f", curvature defocus {dz_curv * 1e3:.3f} mm "
               f"({curvature_meta['curvature_penalty_db']:.2f} dB)"
               if dz_curv else ""))
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
            "f0_m": float(f0),
            "n_comp_modes": residual.n_modes,
            **curvature_meta,
        },
        assumptions=assumptions,
        mean_only=True,
    )


def terrestrial_smf_walkoff_term(scenario, geometry, *, n_grid=64, turbulence=True):
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
            The tx terminal's Transmitter launches the beam; the rx terminal
            supplies the SMF detector with focal_length_m and
            mode_field_radius_m, its pointing_jitter_rad, and its compensation
            stack. The scenario `direction` maps the roles (forward: tx = near,
            rx = far; reverse swaps them).
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
        raise ValueError(
            "terrestrial_smf_walkoff_term needs an SMF detector on the far terminal.")
    f, w_m = _smf_optics(detector, rx.aperture_m, rx.wavelength_m)
    if f is None or w_m is None:
        raise ValueError(
            "terrestrial_smf_walkoff_term needs the fibre-coupling optics to map a "
            "tip-tilt to a focal-plane displacement. Set SMF.focal_length_m and "
            "SMF.mode_field_radius_m, or set SMF.optimal_focus=True to derive the "
            "optimal focal length from the mode field radius."
        )
    dz = float(detector.defocus_m)
    # Effective focal-plane coupling scale for a lateral spot offset. The captured
    # power is the overlap integral of two Gaussians: the spot (1/e^2 radius w_det)
    # and the fibre LP01 mode (radius w_m). The offset dependence is
    # exp(-2*dx^2/(w_det^2+w_m^2)), so w_eff=sqrt(w_det^2+w_m^2). The focal spot
    # radius is w_s=lambda*f/(pi*(D/2)). Off focus, the detector sits at z=f+dz, so
    # the spot grows to w_det=gaussz(w_s, dz): the Gaussian beam radius versus the
    # distance from the waist (Andrews and Phillips 2005, Ch. 4,
    # DOI 10.1117/3.626196). At focus (dz=0) w_det=w_s, so this reduces to the
    # focal-plane case exactly. A point spot (w_det->0) gives w_eff=w_m; that
    # underestimates the tolerated offset and gives too large a loss. Source:
    # two-Gaussian mode overlap, Shaklan and Roddier, Appl. Opt. 27 (1988) 2334,
    # DOI 10.1364/AO.27.002334.
    # The spot grows with the distance from the TRUE focus, not from the nominal
    # focal plane. The received beam is a diverging Gaussian, so its true focus is
    # at z = f + dz_curv, and the fibre is dz_eff = dz - dz_curv away from it (see
    # _received_curvature; S. A. Self, Appl. Opt. 22, 658 (1983), DOI
    # 10.1364/AO.22.000658). The chief-ray LEVERS below keep the PHYSICAL dz,
    # because the detector position, not the focus position, sets them.
    R_rx, dz_curv = _received_curvature(scenario, f)
    dz_eff = dz - dz_curv
    w_s = rx.wavelength_m * f / (np.pi * rx.aperture_m / 2.0)
    w_det = gaussz(w_s, dz_eff, rx.wavelength_m)
    w_eff = np.sqrt(w_det ** 2 + w_m ** 2)
    # Open the collection context around the PHYSICS CALL only. The received
    # tip-tilt reads the decorated beam-wander arrival-tilt kernel (through
    # _received_tiptilt_variance), so the Term inherits its assumptions and
    # carries traced provenance.
    with trace_assumptions() as trace:
        sigma2_theta, meta = _received_tiptilt_variance(scenario, n_grid=n_grid,
                                                        turbulence=turbulence)
    # Per-axis spot offset at the detector plane: the (f+dz) tilt lever (see
    # _spot_offset_sigma). At focus (dz=0) this reduces to f*sqrt(sigma2_theta/2),
    # the focal-plane displacement.
    sigma_d = _spot_offset_sigma(f, dz, sigma2_theta)
    mean, quantile, sampler = _walkoff_faces(w_eff, sigma_d)

    assumptions = trace.merge(
        beam_type=BEAM_GAUSSIAN,
        turbulence_regime=REGIME_WEAK,
        spectrum=SPECTRUM_KOLMOGOROV,
        validity="Receive tip-tilt walk-off of the spot against the single-mode "
                 "fibre mode. The wander arrival tilt is the weak beam-wander "
                 "model (Dios et al. 2004). A receive tip-tilt or AO stage tracks "
                 "it out (all-or-nothing, no bandwidth). The receive mechanical "
                 "jitter adds to it. The fade is exponential in dB. The coupling "
                 "loss has no upper limit: more tilt gives more loss. The spot "
                 "grows to w_det=gaussz(w_s, dz_eff), with dz_eff = "
                 "defocus_m - dz_curv the distance from the TRUE focus of the "
                 "received diverging beam (dz_curv = f^2/(R_rx - f); S. A. Self, "
                 "Appl. Opt. 22, 658 (1983), DOI 10.1364/AO.22.000658). The "
                 "displacement uses the ray-optics chief-ray lever "
                 "(f+dz)*theta with the PHYSICAL dz (Andrews and "
                 "Phillips 2005, Ch. 4, DOI 10.1117/3.626196). The walk-off "
                 "RESPONSE to a displacement is GEOMETRIC ONLY. Only the "
                 "weak-turbulence validity of the tilt variance limits the Term.",
    )
    _flag_curvature(assumptions, R_rx, f, dz_curv, dz)
    # The MEAN modal quadratic-phase penalty of a defocused single-mode fibre is
    # now MODELLED, in the coupling Term (terrestrial_smf_coupling_term, the
    # aberrated closed form eta(a, c)). What stays geometric here is the WALK-OFF
    # RESPONSE: how the coupling falls as the spot slides off the mode. This Term
    # models that as a two-Gaussian overlap of the defocused spot and the fibre
    # mode, so it ignores how the defocus phase itself reshapes the modal overlap
    # against a displacement. So this walk-off fade stays OPTIMISTIC off focus.
    # Flag it loudly, source tagged. The trigger stays the DELIBERATE defocus_m,
    # so a plain focal-plane terrestrial link does not raise it: there the mean
    # curvature penalty is fully modelled in the coupling Term, and only the
    # second-order walk-off response stays geometric.
    if dz != 0.0:
        assumptions.flag(
            f"Off-true-focus single-mode-fibre walk-off (dz_eff="
            f"{dz_eff * 1e3:.3f} mm, defocus_m={dz * 1e3:.3f} mm, dz_curv="
            f"{dz_curv * 1e3:.3f} mm): the MEAN modal quadratic-phase penalty is "
            "modelled in the coupling Term (receive coupling (SMF)), but the "
            "walk-off DISPLACEMENT response here is GEOMETRIC ONLY (spot growth "
            "plus a two-Gaussian overlap). It does not model how the defocus phase "
            "reshapes the modal overlap against a displacement, so this walk-off "
            "fade is OPTIMISTIC. For the full modal treatment use Ruilier and "
            "Cassaing, JOSA A 18 (2001) 143, DOI 10.1364/JOSAA.18.000143, or an "
            "MMF (light bucket) or fidelity-2 model.",
            source="factory:models.coupling.terrestrial",
        )
    # Close the known gap (the walk-off Term declared REGIME_WEAK but never
    # flagged). The beam-wander arrival tilt is a WEAK-turbulence model, but the
    # vendored kernel (olb.turbulence.coupled_flux.beam_wander_variance) carries
    # NO runtime weak-regime check, so the trace alone does not flag a strong
    # path. Flag it here as a scenario-level regime fact (source tagged), only
    # when the wander tilt actually contributes (turbulence on and not tracked
    # out, so sigma2_wander > 0). The gate is the same beam-aware Rytov test the
    # scintillation Term uses (Andrews and Phillips 2005, Ch. 5, Eq. (16), DOI
    # 10.1117/3.626196).
    if meta["sigma2_wander"] > 0.0:
        L = float(scenario.channel.path_length_m)
        cn2 = float(scenario.channel.cn2)
        w0 = scenario.tx_terminal.transmitter.waist_m
        # Read the plane-wave Rytov variance from the andrews kernel, do not
        # re-derive it: sigma_R^2 = 1.23 Cn2 k^(7/6) L^(11/6), Andrews and
        # Phillips 2005, Ch. 8, Eq. (10), printed p. 262, DOI 10.1117/3.626196.
        sigma2_R = float(rytov_variance(rx.wavelength_m, L, cn2))
        Lambda = float(beam_params(w0, rx.wavelength_m, L).lam)
        if rytov_weak(sigma2_R, Lambda) == 'hard':
            assumptions.flag(
                f"sigma2_R={sigma2_R:.3f} (Lambda={Lambda:.3f}) meets or exceeds "
                f"the Gaussian-beam weak limit {WEAK_REGIME_LIMIT}; the beam-wander "
                "arrival-tilt model is not trusted. Use the fidelity-2 Monte Carlo.",
                source="factory:models.coupling.terrestrial",
            )
    return Term(
        name="SMF tip-tilt walk-off",
        category="pointing",
        mean_db=float(mean),
        sampler=sampler,
        quantile=quantile,
        note=f"SMF walk-off, w_eff={w_eff * 1e6:.1f} um "
             f"(w_det={w_det * 1e6:.1f}, w_m={w_m * 1e6:.1f}), f={f * 1e3:.1f} mm, "
             f"defocus={dz * 1e3:.2f} mm, dz_eff={dz_eff * 1e3:.2f} mm, "
             f"wander_tracked={meta['wander_tracked']}",
        meta={**meta, "detector": "SMF", "focal_length_m": f,
              "mode_field_radius_m": w_m, "spot_radius_m": float(w_s),
              "spot_radius_detector_m": float(w_det), "w_eff_m": float(w_eff),
              "defocus_m": dz,
              "received_curvature_m": float(R_rx),
              "curvature_defocus_m": float(dz_curv), "dz_eff_m": float(dz_eff),
              "spot_offset_1sigma_m": float(sigma_d)},
        assumptions=assumptions,
    )


def _mmf_encircled_efficiency(offset_m, w_s, a_core):
    '''
    Return the fraction of the focal-spot power inside the fibre core.

    A multimode fibre is a light bucket. The core is a HARD disk of radius a_core
    in the fibre plane. The focal spot is a Gaussian of 1/e^2 intensity radius w_s.
    When the spot centre is at a radial offset offset_m from the core centre, the
    coupled power is the ENCIRCLED ENERGY of the displaced Gaussian inside the
    disk. So the fibre collects ALL the spot power that lands inside the core. It
    is NOT a mode overlap.

    The Gaussian intensity is a 2-D normal with a per-axis standard deviation
    sigma = w_s/2. So the encircled power is the non-central chi-square CDF:
        eta = P(R <= a_core) = ncx2.cdf( (a_core/sigma)^2, df=2, nc=(offset/sigma)^2 ).
    This equals 1 - Q1(2*offset/w_s, 2*a_core/w_s) with Q1 the Marcum Q-function.
    At offset=0 it reduces to 1 - exp(-2*a_core^2/w_s^2), the on-axis encircled
    energy. Source: Marcum, RAND RM-753 (1950); encircled energy of an off-axis
    Gaussian.

    Parameters:
        offset_m : float or numpy.ndarray
            Radial offset of the spot centre from the core centre [m].
        w_s : float
            Focal spot 1/e^2 intensity radius [m].
        a_core : float
            Core radius [m].

    Returns:
        numpy.ndarray
            Coupled power fraction eta in (0, 1].
    '''
    from scipy.stats import ncx2
    sigma = w_s / 2.0
    nc = (np.asarray(offset_m, dtype=float) / sigma) ** 2
    return ncx2.cdf((a_core / sigma) ** 2, df=2, nc=nc)


def terrestrial_mmf_coupling_term(scenario, geometry, *, n_grid=64, turbulence=True):
    '''
    Build the multimode-fibre coupling Term (terrestrial).

    A multimode fibre is a light bucket of core radius a_core in the fibre plane.
    The focal spot has the diffraction radius w_s = lambda*f / (pi*(D/2)). The
    coupling is the ENCIRCLED ENERGY of the focal spot inside the hard core disk.
    A light bucket collects ALL the spot power that lands inside the core. It is
    NOT a mode overlap.

    The received tip-tilt moves the spot centre off the core by dx = f*theta. So
    the coupled power is the encircled energy of the DISPLACED Gaussian spot inside
    the core:
        eta(dx) = 1 - Q1(2*dx/w_s, 2*a_core/w_s),
    with Q1 the Marcum Q-function (see _mmf_encircled_efficiency). At dx=0 this
    reduces to the on-axis encircled energy eta_static = 1 - exp(-2 a_core^2/w_s^2).
    A small spot deep inside the core loses nothing until it nears the edge (a
    flat-top acceptance). A spot larger than the core overfills it and loses power
    from dx=0. The tilt dx has two i.i.d. Gaussian axes, so its radial magnitude is
    Rayleigh. The Term averages the loss over that Rayleigh offset, and it carries
    a real fade.

    Anti-double-count: this Term is a MULTIPLICATIVE fibre-coupling efficiency on
    the aperture-collected field. It adds only -10log10(eta_mmf). It does NOT
    re-add the aperture capture (the geometric spreading Term already carries the
    free-space spread and the aperture capture). This is the same pattern the SMF
    Term uses. The receive mechanical jitter here uses the RECEIVE terminal, so it
    does not double-count the transmit pointing Term.

    With turbulence=False the received tip-tilt keeps only the receive mechanical
    jitter (the beam-wander tilt drops). So the coupling is the encircled energy
    with only the jitter offset. The on-axis static loss stays, because it is a
    fixed optical loss, not turbulence.

    Parameters:
        scenario : TerrestrialScenario
            The tx terminal's Transmitter launches the beam; the rx terminal
            supplies the MMF detector, aperture, pointing_jitter_rad, and
            compensation stack. The scenario `direction` maps the roles
            (forward: tx = near, rx = far; reverse swaps them).
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
        raise ValueError(
            "terrestrial_mmf_coupling_term needs an MMF detector on the far terminal.")
    D = rx.aperture_m
    wavelength = rx.wavelength_m
    f = _mmf_focal_length(detector, D, wavelength)
    if f is None:
        raise ValueError(
            "terrestrial_mmf_coupling_term needs a focal length to map a tip-tilt "
            "to a focal-plane displacement. Set MMF.focal_length_m, or set "
            "MMF.optimal_focus=True to match the spot to the core."
        )
    a_core = detector.core_radius_m
    dz = float(detector.defocus_m)

    # Diffraction focal spot radius (1/e^2), Gaussian approximation to the Airy.
    w_s = wavelength * f / (np.pi * (D / 2.0))
    # Defocused spot radius. The received beam is a DIVERGING Gaussian, so its
    # true focus is at z = f + dz_curv, not at z = f (see _received_curvature;
    # S. A. Self, Appl. Opt. 22, 658 (1983), DOI 10.1364/AO.22.000658). The
    # detector sits at z = f + dz, so it is dz_eff = dz - dz_curv from the true
    # focus, and the spot grows to w_det = gaussz(w_s, dz_eff) (Gaussian beam
    # radius versus distance from the waist, Andrews and Phillips 2005, Ch. 4,
    # DOI 10.1117/3.626196). At large |dz_eff| this tends to the geometric blur
    # (D/2)*|dz_eff|/f. The coupling scale (the eta_static and the offset
    # averaging) uses w_det, so a spot away from the true focus spills more from
    # the core: this is the mean defocus loss. The curvature defocus is ALWAYS
    # charged; set defocus_m = curvature_focus_shift(scenario) for a coupler
    # aligned on the true focus.
    R_rx, dz_curv = _received_curvature(scenario, f)
    dz_eff = dz - dz_curv
    w_det = gaussz(w_s, dz_eff, wavelength)

    # Angular acceptance gate (numerical aperture). The focusing optic makes a cone
    # of half-angle NA_optic = (D/2)/f. A ray from aperture radius rho focuses at
    # angle rho/f, so only rays from rho <= f*NA_fibre stay within the acceptance
    # cone theta_a = arcsin(NA). For a uniform aperture the guided POWER fraction is
    # (f*NA/(D/2))^2 = (NA/NA_optic)^2, clipped to 1. This is a flat multiplicative
    # loss on the coupled power; when NA_optic <= NA it is 1 (no loss). It is the
    # etendue penalty a core-radius-only bucket misses. Source: Snyder and Love,
    # Optical Waveguide Theory (1983), DOI 10.1007/978-1-4613-2813-1.
    na_optic = (D / 2.0) / f
    na_fibre = detector.numerical_aperture
    if na_fibre is not None:
        na_factor = float(min(1.0, (na_fibre / na_optic) ** 2))
    else:
        na_factor = 1.0

    # The received tip-tilt moves the spot centre off the core. The per-axis spot
    # offset has a 1-sigma of sigma_d = sqrt((f+dz)^2*(sigma2_theta/2)): the
    # ray-optics chief-ray lever (see _spot_offset_sigma). At focus (dz=0) this
    # reduces to f*sqrt(sigma2_theta/2), the focal-plane displacement. The coupled
    # power is the encircled energy of the displaced Gaussian spot inside the hard
    # core (a light bucket, NOT a mode overlap). See _mmf_encircled_efficiency.
    #
    # Open the collection context around the PHYSICS CALL only. The received
    # tip-tilt reads the decorated beam-wander arrival-tilt kernel, so the Term
    # inherits its assumptions and carries traced provenance.
    with trace_assumptions() as trace:
        sigma2_theta, meta = _received_tiptilt_variance(scenario, n_grid=n_grid,
                                                        turbulence=turbulence)
    sigma_d = _spot_offset_sigma(f, dz, sigma2_theta)   # spot offset [m]
    eta_static = na_factor * float(_mmf_encircled_efficiency(0.0, w_det, a_core))
    static_db = -10.0 * np.log10(eta_static)

    def _loss_db(offset_m):
        eta = np.clip(na_factor * _mmf_encircled_efficiency(offset_m, w_det, a_core),
                      1e-300, None)
        return -10.0 * np.log10(eta)

    # Mean over the Rayleigh offset distribution (two i.i.d. Gaussian axes). A
    # quadrature over the radial offset gives the expected loss. The 8-sigma grid
    # covers the Rayleigh tail (the pdf there is exp(-32), negligible).
    if sigma_d > 0.0:
        dd = np.linspace(0.0, 8.0 * sigma_d, 2000)
        rayleigh = dd / sigma_d ** 2 * np.exp(-dd ** 2 / (2.0 * sigma_d ** 2))
        mean_db = float(np.trapezoid(_loss_db(dd) * rayleigh, dd))
    else:
        mean_db = static_db

    def quantile(p):
        # The loss rises with the offset, and the radial offset is Rayleigh, so the
        # p-quantile of the loss is the loss at the p-quantile of the offset.
        d_p = sigma_d * np.sqrt(-2.0 * np.log(1.0 - p))
        return float(_loss_db(d_p))

    def sampler(n, rng):
        dx = rng.normal(0.0, sigma_d, n)
        dy = rng.normal(0.0, sigma_d, n)
        return _loss_db(np.sqrt(dx ** 2 + dy ** 2))

    # The traced beam-wander kernel owns the beam type, the weak regime, and the
    # spectrum; the merge inherits its union and provenance. State the three
    # headline fields explicitly (this is a Gaussian-beam light-bucket Term).
    assumptions = trace.merge(
        beam_type=BEAM_GAUSSIAN,
        turbulence_regime=REGIME_WEAK,
        spectrum=SPECTRUM_KOLMOGOROV,
        validity="Multimode-fibre coupling (light bucket): the coupled power is the "
                 "encircled energy of the Gaussian spot (1/e^2 radius w_det = "
                 "gaussz(w_s, dz), w_s = lambda*f/(pi*(D/2))) inside the hard core "
                 "disk of radius a_core, displaced by the received tip-tilt. It is the "
                 "non-central chi-square CDF (Marcum Q), so it collects ALL the spot "
                 "power inside the core, not a mode overlap. The DEFOCUS model is "
                 "geometric. The received beam is a diverging Gaussian, so its true "
                 "focus is at z = f + dz_curv with dz_curv = f^2/(R_rx - f) "
                 "(S. A. Self, Appl. Opt. 22, 658 (1983), DOI 10.1364/AO.22.000658), "
                 "and R_rx is the received phase-front radius (Andrews and Phillips "
                 "2005, Ch. 4, DOI 10.1117/3.626196). This curvature defocus is "
                 "ALWAYS charged: the detector sits at z = f + defocus_m, so the "
                 "spot grows over dz_eff = defocus_m - dz_curv, while the "
                 "displacement uses the ray-optics chief-ray lever "
                 "(f+dz)*theta with the PHYSICAL dz (Andrews and "
                 "Phillips 2005, Ch. 4, DOI 10.1117/3.626196; ray-optics chief-ray "
                 "of a thin lens). optimal_focus is a focal-LENGTH rule; it never "
                 "moves the detector. The "
                 "tip-tilt is the weak beam wander (tracked by a tip-tilt or AO "
                 "stage) plus the receive mechanical jitter. The spot model assumes "
                 "a uniform, unobscured "
                 "circular aperture. The numerical-aperture gate (when set) is a "
                 "flat power-transmission factor, NOT a re-truncated aperture: it "
                 "does not re-broaden the focal spot.",
    )
    # Scenario-level facts the physics never sees stay factory flags, source
    # tagged so they carry the same "[source] reason" prefix as a traced check.
    _flag_curvature(assumptions, R_rx, f, dz_curv, dz)
    if na_fibre is not None and na_optic > na_fibre:
        assumptions.flag(
            f"The focusing cone NA_optic={na_optic:.3f} exceeds the fibre "
            f"NA={na_fibre:.3f}; the fibre does not guide the steep rays. The "
            f"angular gate cuts the coupled power by {(-10.0 * np.log10(na_factor)):.2f} "
            "dB. Shorten nothing further, or use a larger-NA fibre.",
            source="factory:models.coupling.terrestrial",
        )
        assumptions.flag(
            "The NA gate (NA/NA_optic)^2 is an aperture-AREA fraction, so it assumes "
            "a uniformly illuminated pupil that fills the aperture. Here D is the "
            "ILLUMINATED diameter, not only the mechanical stop. An underfilled or "
            "Gaussian-apodized pupil carries less power in the steep marginal rays, "
            "so its true NA loss is SMALLER. Thus this gate is CONSERVATIVE "
            "(pessimistic) for an underfilled pupil, and exact for a filled uniform "
            "one. A distant-source receive beam is near uniform, so the assumption "
            "holds for a receive aperture.",
            source="factory:models.coupling.terrestrial",
        )
    if rx.obscuration_ratio > 0.0:
        assumptions.flag(
            f"The far aperture has a central obscuration "
            f"(ratio={rx.obscuration_ratio:.3f}); the Gaussian focal-spot model "
            "assumes a uniform circular aperture and does not model it.",
            source="factory:models.coupling.terrestrial",
        )
    return Term(
        name="receive coupling (MMF)",
        category="coupling",
        mean_db=float(mean_db),
        sampler=sampler,
        quantile=quantile,
        note=f"MMF coupling, a_core={a_core * 1e6:.1f} um, "
             f"w_det={w_det * 1e6:.1f} um (w_s={w_s * 1e6:.1f}), "
             f"defocus={dz * 1e3:.2f} mm, dz_curv={dz_curv * 1e3:.2f} mm, "
             f"dz_eff={dz_eff * 1e3:.2f} mm, static={static_db:.2f} dB"
             + (f", NA gate {(-10.0 * np.log10(na_factor)):.2f} dB"
                if na_factor < 1.0 else ""),
        meta={**meta, "detector": "MMF", "core_radius_m": a_core,
              "focal_length_m": f, "spot_radius_m": float(w_s),
              "spot_radius_detector_m": float(w_det),
              "defocus_m": dz,
              "received_curvature_m": float(R_rx),
              "curvature_defocus_m": float(dz_curv), "dz_eff_m": float(dz_eff),
              "eta_static": eta_static, "static_loss_db": float(static_db),
              "walkoff_mean_db": float(mean_db - static_db),
              "spot_offset_1sigma_m": float(sigma_d),
              "offset_tilt_1sigma_m": float((f + dz) * np.sqrt(sigma2_theta / 2.0)),
              "numerical_aperture": na_fibre, "na_optic": float(na_optic),
              "na_factor": na_factor,
              "na_gate_loss_db": float(-10.0 * np.log10(na_factor))},
        assumptions=assumptions,
    )


if __name__ == '__main__':
    import warnings

    from ...scenario import TerrestrialScenario, TerrestrialChannel
    from ...geometry import HorizontalPath
    from ...terminal import Terminal, Transmitter, Aperture, SMF, MMF, TipTilt

    lam = 1550e-9
    D_test, wm = 0.2, 5.2e-6
    f_opt = np.pi * (D_test / 2.0) * wm / (lam * 1.12)   # a = 1.12

    def _terr(detector, *, jitter=0.0, compensation=None, cn2=1e-14,
              far_aperture=0.2, w0=0.02, L=3e3, divergence=None):
        return TerrestrialScenario(
            near=Terminal(aperture_m=0.3, wavelength_m=lam,
                          transmitter=Transmitter(waist_m=w0, power_dbm=30,
                                                  divergence_rad=divergence)),
            far=Terminal(aperture_m=far_aperture, wavelength_m=lam,
                         pointing_jitter_rad=jitter,
                         detector=detector, compensation=compensation or []),
            channel=TerrestrialChannel(path_length_m=L, attenuation_db_per_km=0.5,
                                       cn2=cn2))

    hpath = HorizontalPath(3e3)

    # --- Received tip-tilt: a tip-tilt stage tracks out the wander -----------
    v_full, _ = _received_tiptilt_variance(_terr(Aperture(), jitter=5e-6), n_grid=64)
    v_track, m_track = _received_tiptilt_variance(
        _terr(Aperture(), jitter=5e-6, compensation=[TipTilt()]), n_grid=64)
    assert v_full > v_track, (v_full, v_track)         # tracking removes the wander
    assert m_track["wander_tracked"] and m_track["sigma2_wander"] == 0.0
    assert np.isclose(v_track, 2.0 * 5e-6 ** 2)        # only the jitter remains

    # --- SMF walk-off: a real fade. A larger jitter deepens the loss ---------
    smf_opt = SMF(focal_length_m=f_opt, mode_field_radius_m=wm, sensitivity_dbm=-40)
    wo_small = terrestrial_smf_walkoff_term(_terr(smf_opt, jitter=2e-6), hpath)
    wo_big = terrestrial_smf_walkoff_term(_terr(smf_opt, jitter=10e-6), hpath)
    assert wo_small.stochastic and wo_small.quantile is not None
    assert wo_big.mean_db > wo_small.mean_db
    wo_big_q = wo_big.quantile_db(0.99)
    assert wo_big_q is not None and wo_big_q > wo_big.mean_db   # exponential tail
    rng = np.random.default_rng(0)
    draws = wo_big.sample_db(50_000, rng)
    assert np.abs(draws.mean() - wo_big.mean_db) / wo_big.mean_db < 0.03
    # optimal_focus: derive f from a=1.12 (the far aperture is D_test=0.2).
    smf_focus = SMF(mode_field_radius_m=wm, optimal_focus=True, sensitivity_dbm=-40)
    wo_focus = terrestrial_smf_walkoff_term(_terr(smf_focus, jitter=10e-6), hpath)
    assert np.isclose(wo_focus.meta["focal_length_m"], f_opt)
    assert np.isclose(wo_focus.mean_db, wo_big.mean_db)
    wo_default = terrestrial_smf_walkoff_term(
        _terr(SMF(optimal_focus=True), jitter=10e-6), hpath)
    assert np.isclose(wo_default.mean_db, wo_big.mean_db)
    # An SMF without the optics cannot map a tip-tilt to a displacement.
    try:
        terrestrial_smf_walkoff_term(_terr(SMF()), hpath)
        raise AssertionError("SMF walk-off without optics must raise")
    except ValueError:
        pass

    # --- MMF coupling: encircled energy of the displaced spot in the core ----
    a_core = 25e-6
    # On-axis reduces to the encircled energy 1 - exp(-2 a_core^2/w_s^2).
    assert np.isclose(float(_mmf_encircled_efficiency(0.0, 5e-6, a_core)),
                      1.0 - np.exp(-2.0 * a_core ** 2 / (5e-6) ** 2))
    # Flat-top acceptance: a small spot deep inside the core loses almost nothing.
    assert float(_mmf_encircled_efficiency(0.1 * a_core, 1e-6, a_core)) > 0.999
    # The spot centre at the core edge collects about half the power (about 3 dB).
    e_edge = float(_mmf_encircled_efficiency(a_core, 1e-6, a_core))
    assert 0.45 < e_edge < 0.55, e_edge
    # The coupling falls monotonically as the offset grows.
    offs = np.array([0.0, 0.5, 1.0, 1.5, 2.0]) * a_core
    etas = _mmf_encircled_efficiency(offs, 8e-6, a_core)
    assert np.all(np.diff(etas) < 0.0), etas

    # An MMF Term at optimal focus, so a tip-tilt walks it off. A real fade.
    mmf = MMF(core_radius_m=a_core, optimal_focus=True, sensitivity_dbm=-38)
    t_mmf = terrestrial_mmf_coupling_term(_terr(mmf, jitter=10e-6, cn2=1e-15), hpath)
    assert t_mmf.name == "receive coupling (MMF)" and t_mmf.category == "coupling"
    assert t_mmf.stochastic and t_mmf.quantile is not None and not t_mmf.mean_only
    assert 0.0 < t_mmf.meta["eta_static"] <= 1.0 and t_mmf.meta["static_loss_db"] >= 0.0
    t_mmf_q = t_mmf.quantile_db(0.99)
    assert t_mmf_q is not None and t_mmf_q > t_mmf.mean_db  # walk-off adds a fade
    t_mmf_calm = terrestrial_mmf_coupling_term(
        _terr(mmf, jitter=2e-6, cn2=1e-15), hpath)
    assert t_mmf.mean_db > t_mmf_calm.mean_db               # more jitter, more loss
    # optimal_focus derives f to match the spot to the core (a_core/w_s=1.12).
    f_mmf = np.pi * (D_test / 2.0) * a_core / (lam * 1.12)
    t_explicit = terrestrial_mmf_coupling_term(
        _terr(MMF(core_radius_m=a_core, focal_length_m=f_mmf), jitter=10e-6,
              cn2=1e-15), hpath)
    assert np.isclose(t_mmf.meta["focal_length_m"], f_mmf)
    assert np.isclose(t_mmf.mean_db, t_explicit.mean_db)
    # The received curvature is ALWAYS charged, so a fibre AT the focal plane is
    # dz_curv away from the true focus and does NOT reach the flat-wavefront
    # a_core/w_s=1.12 value. Move the fibre to the true focus and it does.
    assert t_mmf.meta["curvature_defocus_m"] > 0.0
    t_tracked = terrestrial_mmf_coupling_term(
        _terr(MMF(core_radius_m=a_core, focal_length_m=f_mmf,
                  defocus_m=t_mmf.meta["curvature_defocus_m"]),
              jitter=10e-6, cn2=1e-15), hpath)
    assert np.isclose(t_tracked.meta["dz_eff_m"], 0.0, atol=1e-12)
    assert np.isclose(t_tracked.meta["eta_static"], 1.0 - np.exp(-2.0 * 1.12 ** 2),
                      atol=1e-3), t_tracked.meta["eta_static"]
    assert t_tracked.meta["static_loss_db"] < t_mmf.meta["static_loss_db"]
    # curvature_focus_shift is the public route to that value.
    assert np.isclose(
        curvature_focus_shift(_terr(MMF(core_radius_m=a_core,
                                        focal_length_m=f_mmf))),
        t_mmf.meta["curvature_defocus_m"])
    assert t_mmf.meta["numerical_aperture"] is None and t_mmf.meta["na_factor"] == 1.0
    # Numerical-aperture angular gate.
    na_optic = t_mmf.meta["na_optic"]
    t_wide_na = terrestrial_mmf_coupling_term(
        _terr(MMF(core_radius_m=a_core, focal_length_m=f_mmf,
                  numerical_aperture=10.0 * na_optic), jitter=10e-6, cn2=1e-15), hpath)
    assert np.isclose(t_wide_na.meta["na_factor"], 1.0)
    assert np.isclose(t_wide_na.mean_db, t_mmf.mean_db)     # no gate: unchanged
    t_tight_na = terrestrial_mmf_coupling_term(
        _terr(MMF(core_radius_m=a_core, focal_length_m=f_mmf,
                  numerical_aperture=0.5 * na_optic), jitter=10e-6, cn2=1e-15), hpath)
    assert np.isclose(t_tight_na.meta["na_factor"], 0.25)   # (0.5)^2
    assert np.isclose(t_tight_na.mean_db - t_mmf.mean_db, -10.0 * np.log10(0.25), atol=1e-6)
    assert any("does not guide the steep rays" in v
               for v in t_tight_na.assumptions.violations)
    # An MMF with no focal length and no optimal_focus is refused.
    try:
        terrestrial_mmf_coupling_term(_terr(MMF(core_radius_m=a_core)), hpath)
        raise AssertionError("MMF without a focal length must raise")
    except ValueError:
        pass

    # --- Gap 3: the launch curvature reaches the Fried parameter -------------
    smf_plain = SMF(sensitivity_dbm=-40)
    c_coll = terrestrial_smf_coupling_term(_terr(smf_plain), hpath)
    theta_min_terr = lam / (np.pi * 0.02)
    c_div = terrestrial_smf_coupling_term(
        _terr(smf_plain, divergence=5 * theta_min_terr), hpath)
    assert np.isinf(c_coll.meta["f0_m"]) and c_div.meta["f0_m"] < 0.0
    # A diverged beam is more spherical (Theta0 > 1), so the transmitter-referred
    # path weight drops: r0 grows, the residual falls, the loss falls.
    assert c_div.meta["r0_m"] > c_coll.meta["r0_m"], (
        c_div.meta["r0_m"], c_coll.meta["r0_m"])
    assert c_div.mean_db < c_coll.mean_db, (c_div.mean_db, c_coll.mean_db)

    # --- turbulence=False: static coupling + jitter, no turbulence quantity --
    v_off, m_off = _received_tiptilt_variance(_terr(Aperture(), jitter=5e-6),
                                              n_grid=64, turbulence=False)
    assert m_off["sigma2_wander"] == 0.0 and np.isclose(v_off, 2.0 * 5e-6 ** 2)
    # The SMF walk-off then carries the jitter fade with NO turbulence.
    wo_off = terrestrial_smf_walkoff_term(
        _terr(smf_opt, jitter=10e-6), hpath, turbulence=False)
    assert wo_off.meta["sigma2_wander"] == 0.0 and wo_off.meta["sigma2_jitter"] > 0.0
    wo_off_q = wo_off.quantile_db(0.99)
    assert wo_off_q is not None and wo_off_q > wo_off.mean_db
    # An SMF static term needs no launch beam (turbulence off skips the r0 path).
    terr_off = terrestrial_smf_coupling_term(
        TerrestrialScenario(
            near=Terminal(aperture_m=0.3, wavelength_m=lam),   # no transmitter
            far=Terminal(aperture_m=0.2, wavelength_m=lam, detector=SMF()),
            channel=TerrestrialChannel(path_length_m=3e3, attenuation_db_per_km=0.5,
                                       cn2=1e-14)),
        hpath, turbulence=False)
    assert terr_off.meta["model"] == "static" and not terr_off.mean_only

    # --- WP3c: traced provenance, the walk-off gap closure, and the guard ----
    from ...results import Budget

    # (1) The coupling and walk-off Terms now carry traced provenance naming
    #     their physics sources.
    smf_plain2 = SMF(sensitivity_dbm=-40)
    cpl = terrestrial_smf_coupling_term(_terr(smf_plain2, cn2=1e-15), hpath)
    assert cpl.assumptions.provenance, "SMF coupling must carry traced provenance"
    assert any("gaussian_fried_parameter_profile" in s
               for s in cpl.assumptions.provenance), cpl.assumptions.provenance
    assert any("apply_compensation" in s
               for s in cpl.assumptions.provenance), cpl.assumptions.provenance
    wo_prov = terrestrial_smf_walkoff_term(_terr(smf_opt, jitter=10e-6, cn2=1e-15),
                                           hpath)
    assert any("wander_arrival_angle_variance" in s
               for s in wo_prov.assumptions.provenance), wo_prov.assumptions.provenance
    mmf_prov = terrestrial_mmf_coupling_term(_terr(mmf, jitter=10e-6, cn2=1e-15),
                                             hpath)
    assert mmf_prov.assumptions.provenance, "MMF coupling must carry provenance"

    # (2) THE GAP CLOSED. The walk-off Term declared REGIME_WEAK but NEVER
    #     flagged. In strong turbulence it now reads not-ok; in genuinely weak
    #     turbulence it stays ok. No spurious field-region violation reaches it:
    #     the walk-off reads wander_arrival_angle_variance, NOT the z=1.0-carrier
    #     aperture_arrival_angle_variance, so the delegate Fresnel-zone check is
    #     never on the trace.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        wo_weak = terrestrial_smf_walkoff_term(
            _terr(smf_opt, jitter=10e-6, cn2=1e-16), hpath)
        wo_strong = terrestrial_smf_walkoff_term(
            _terr(smf_opt, jitter=10e-6, cn2=1e-13), hpath)
    assert wo_weak.assumptions.ok, wo_weak.assumptions.violations
    assert not wo_strong.assumptions.ok, \
        "the walk-off Term must flag in strong turbulence (gap closed)"
    assert any("beam-wander arrival-tilt model is not trusted" in v
               for v in wo_strong.assumptions.violations), \
        wo_strong.assumptions.violations
    assert not any("field-region" in v.lower() or "Fresnel" in v
                   for v in wo_strong.assumptions.violations), \
        "no spurious field-region violation must reach the walk-off Term"

    # (3) Budget.check() reports NO untraced-guard entry for the coupling Terms.
    guard_reason = "did not open the assumption collection context"
    for term in (cpl, mmf_prov):
        found = [(n, r) for n, r in Budget([term]).check(warn=False)
                 if guard_reason in r]
        assert found == [], found

    # --- Non-focal-plane (defocus) sensing -----------------------------------
    a_core_d = 25e-6
    f_mmf_d = np.pi * (D_test / 2.0) * a_core_d / (lam * 1.12)     # optimal focus
    w_s_d = lam * f_mmf_d / (np.pi * (D_test / 2.0))
    dz_big = 5e-3                                                  # 5 mm defocus

    def _mmf_det(dz=0.0):
        return MMF(core_radius_m=a_core_d, focal_length_m=f_mmf_d, defocus_m=dz,
                   sensitivity_dbm=-38)

    # (a) dz=0 reproduces the existing focal-plane result EXACTLY (bit-identical).
    m_focus = terrestrial_mmf_coupling_term(
        _terr(_mmf_det(0.0), jitter=10e-6, cn2=1e-15), hpath)
    m_focus_ref = terrestrial_mmf_coupling_term(
        _terr(MMF(core_radius_m=a_core_d, focal_length_m=f_mmf_d), jitter=10e-6,
              cn2=1e-15), hpath)
    assert m_focus.mean_db == m_focus_ref.mean_db, (m_focus.mean_db, m_focus_ref.mean_db)
    assert m_focus.meta["defocus_m"] == 0.0

    # (b) At large dz the spot grows and spills from the core, so the loss changes.
    m_dz = terrestrial_mmf_coupling_term(
        _terr(_mmf_det(dz_big), jitter=10e-6, cn2=1e-15), hpath)

    # (c) w_det grows with |dz_eff| and equals gaussz(w_s, dz_eff), with
    #     dz_eff = defocus_m - dz_curv the distance from the TRUE focus. The
    #     curvature defocus is always charged, so a fibre at the focal plane
    #     (defocus_m=0) already carries the spot growth of |dz_curv|.
    dzc = m_focus.meta["curvature_defocus_m"]
    assert dzc > 0.0                       # a collimated 3 km launch diverges
    assert np.isclose(m_focus.meta["dz_eff_m"], -dzc)
    assert np.isclose(m_focus.meta["spot_radius_detector_m"],
                      gaussz(w_s_d, -dzc, lam))
    assert np.isclose(m_dz.meta["spot_radius_detector_m"],
                      gaussz(w_s_d, dz_big - dzc, lam))
    # The spot is symmetric about the TRUE focus, not about the focal plane.
    delta = 3e-3
    m_plus = terrestrial_mmf_coupling_term(
        _terr(_mmf_det(dzc + delta), jitter=10e-6, cn2=1e-15), hpath)
    m_minus = terrestrial_mmf_coupling_term(
        _terr(_mmf_det(dzc - delta), jitter=10e-6, cn2=1e-15), hpath)
    assert np.isclose(m_plus.meta["spot_radius_detector_m"],
                      m_minus.meta["spot_radius_detector_m"])
    # Far from the true focus the spot tends to the geometric blur
    # (D/2)*|dz_eff|/f.
    dz_huge = 0.05
    m_huge = terrestrial_mmf_coupling_term(
        _terr(MMF(core_radius_m=a_core_d, focal_length_m=f_mmf_d, defocus_m=dz_huge),
              jitter=1e-9, cn2=1e-18), hpath)
    blur = (D_test / 2.0) * abs(m_huge.meta["dz_eff_m"]) / f_mmf_d
    assert np.isclose(m_huge.meta["spot_radius_detector_m"], blur, rtol=0.02), (
        m_huge.meta["spot_radius_detector_m"], blur)

    # (d) The mean loss grows with |dz_eff| (the spot spills from the core). This
    #     is the static-spill effect, so read it in the low-jitter limit where the
    #     walk-off does not dominate. (With a large tip-tilt a bigger defocused
    #     spot is more robust to walk-off, so the TOTAL mean can be non-monotonic;
    #     the static spill always grows.) The BEST plane is the true focus, so
    #     moving the fibre from the focal plane TOWARD dz_curv wins power back.
    assert m_dz.meta["static_loss_db"] < m_focus.meta["static_loss_db"]
    m_focus_calm = terrestrial_mmf_coupling_term(
        _terr(_mmf_det(0.0), jitter=1e-7, cn2=1e-18), hpath)
    m_track_calm = terrestrial_mmf_coupling_term(
        _terr(_mmf_det(dzc), jitter=1e-7, cn2=1e-18), hpath)
    m_far_calm = terrestrial_mmf_coupling_term(
        _terr(_mmf_det(2.0 * dzc), jitter=1e-7, cn2=1e-18), hpath)
    assert m_track_calm.mean_db < m_focus_calm.mean_db, (m_track_calm.mean_db,
                                                         m_focus_calm.mean_db)
    assert m_far_calm.mean_db > m_track_calm.mean_db
    # The tracked (true-focus) MMF Term recovers the near-flat-wavefront static
    # coupling of the a_core/w_s=1.12 spot match, about 0.37 dB.
    assert np.isclose(m_track_calm.meta["static_loss_db"], 0.37, atol=0.02), \
        m_track_calm.meta["static_loss_db"]

    # --- SMF walk-off defocus (geometric only, with a loud flag) --------------
    smf_d = SMF(focal_length_m=f_opt, mode_field_radius_m=wm, sensitivity_dbm=-40)
    smf_dz = SMF(focal_length_m=f_opt, mode_field_radius_m=wm, defocus_m=dz_big,
                 sensitivity_dbm=-40)
    # dz=0 is unchanged versus the plain focal-plane SMF walk-off.
    wo_focus0 = terrestrial_smf_walkoff_term(_terr(smf_d, jitter=10e-6), hpath)
    assert wo_focus0.mean_db == wo_big.mean_db and wo_focus0.meta["defocus_m"] == 0.0
    assert not any("response here is GEOMETRIC ONLY" in v
                   for v in wo_focus0.assumptions.violations)
    # Off focus the spot grows (a larger w_det), so the walk-off response changes.
    wo_dz = terrestrial_smf_walkoff_term(_terr(smf_dz, jitter=10e-6), hpath)
    assert (wo_dz.meta["spot_radius_detector_m"]
            > wo_focus0.meta["spot_radius_detector_m"])
    # The loud SMF defocus flag fires when dz != 0 (optimistic geometric model).
    assert any("response here is GEOMETRIC ONLY" in v
               for v in wo_dz.assumptions.violations), wo_dz.assumptions.violations
    assert not wo_dz.assumptions.ok
    assert np.isclose(wo_dz.meta["spot_radius_detector_m"], gaussz(
        lam * f_opt / (np.pi * D_test / 2.0), wo_dz.meta["dz_eff_m"], lam))
    # The walk-off spot too grows from the TRUE focus, so the fibre AT the focal
    # plane already carries the curvature defocus.
    assert wo_focus0.meta["curvature_defocus_m"] > 0.0
    assert np.isclose(wo_focus0.meta["dz_eff_m"],
                      -wo_focus0.meta["curvature_defocus_m"])

    # --- the validation/defocus report scenario ------------------------------
    # lambda 1550 nm, L = 5 km, collimated w0 = 0.02 m, D = 0.2 m, core 25 um,
    # f = 4.5242 m. See validation/defocus/fidelity2_mmf_coupling_gap.md.
    f_rep = np.pi * (D_test / 2.0) * 25e-6 / (lam * 1.12)
    rep_kw = dict(jitter=1e-9, cn2=1e-18, w0=0.02, L=5e3, far_aperture=D_test)
    rep_mmf = terrestrial_mmf_coupling_term(
        _terr(MMF(core_radius_m=25e-6, focal_length_m=f_rep, sensitivity_dbm=-38),
              **rep_kw), hpath)
    assert np.isclose(rep_mmf.meta["received_curvature_m"], 5131.5, rtol=1e-3), \
        rep_mmf.meta["received_curvature_m"]
    assert np.isclose(rep_mmf.meta["curvature_defocus_m"], 3.99e-3, rtol=2e-3), \
        rep_mmf.meta["curvature_defocus_m"]
    # The fibre AT the focal plane reads about 8.5 dB (the analytic Gaussian-spot
    # value; fidelity 2 reads about 7.1 dB, the known 2-W1 Airy-versus-Gaussian
    # gap).
    assert np.isclose(rep_mmf.meta["static_loss_db"], 8.54, atol=0.1), \
        rep_mmf.meta["static_loss_db"]
    # The SMF closed form at the same optics: a = 1.12, w_m = 25 um, so
    # c = -3.95 rad and eta = 0.215 (6.68 dB).
    rep_smf = terrestrial_smf_coupling_term(
        _terr(SMF(focal_length_m=f_rep, mode_field_radius_m=25e-6,
                  sensitivity_dbm=-40), **rep_kw), hpath)
    assert np.isclose(rep_smf.meta["defocus_coefficient_rad"], -3.95, atol=0.01), \
        rep_smf.meta["defocus_coefficient_rad"]
    assert np.isclose(rep_smf.meta["eta_max"], 0.215, atol=5e-4), \
        rep_smf.meta["eta_max"]
    assert np.isclose(rep_smf.meta["curvature_penalty_db"], 5.79, atol=0.02), \
        rep_smf.meta["curvature_penalty_db"]
    # A fibre moved to the TRUE focus pays no curvature penalty at all.
    rep_smf_tracked = terrestrial_smf_coupling_term(
        _terr(SMF(focal_length_m=f_rep, mode_field_radius_m=25e-6,
                  defocus_m=rep_mmf.meta["curvature_defocus_m"],
                  sensitivity_dbm=-40), **rep_kw), hpath)
    assert np.isclose(rep_smf_tracked.meta["curvature_penalty_db"], 0.0, atol=1e-9)
    assert np.isclose(rep_smf_tracked.meta["eta_max"], 0.8145, atol=1e-3)

    wo_big_q99 = wo_big.quantile_db(0.99)
    t_mmf_q99 = t_mmf.quantile_db(0.99)
    print(f"SMF walk-off (10 urad): mean {wo_big.mean_db:.3f} dB  "
          f"99% {float(wo_big_q99) if wo_big_q99 is not None else float('nan'):.3f} dB")
    print(f"MMF (25 um core, 5 urad): static {t_mmf.meta['static_loss_db']:.3f} dB  "
          f"mean {t_mmf.mean_db:.3f} dB  "
          f"99% {float(t_mmf_q99) if t_mmf_q99 is not None else float('nan'):.3f} dB")
    print(f"report scenario (5 km, D=0.2 m, f={f_rep:.4f} m): "
          f"R_rx={rep_mmf.meta['received_curvature_m']:.1f} m, "
          f"dz_curv={rep_mmf.meta['curvature_defocus_m'] * 1e3:.3f} mm")
    print(f"  MMF fibre at f: static {rep_mmf.meta['static_loss_db']:.2f} dB")
    print(f"  SMF fibre at f: eta {rep_smf.meta['eta_max']:.3f} "
          f"({-10 * np.log10(rep_smf.meta['eta_max']):.2f} dB, curvature penalty "
          f"{rep_smf.meta['curvature_penalty_db']:.2f} dB)")
    print("coupling terrestrial self-check passed")
