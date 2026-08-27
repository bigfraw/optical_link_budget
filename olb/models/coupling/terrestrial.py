'''
Receive-coupling Terms for a terrestrial (horizontal) link.

This module builds the receive-coupling Terms of a terrestrial far terminal: the
single-mode-fibre mean coupling, the single-mode-fibre tip-tilt walk-off fade, and
the multimode-fibre (light-bucket) coupling. The shared single-mode-fibre coupling
physics lives in olb.models.coupling._common; this module reads it. The received
tip-tilt (beam wander plus receive mechanical jitter) lives here, because only the
terrestrial Terms use it.

Sources:
  Marechal / Dikmelik-Davidson coupling: see olb.models.coupling._common.
  eta_max(a): Shaklan and Roddier, Appl. Opt. 27, 2334 (1988), DOI
  10.1364/AO.27.002334.
  Beam-wander arrival tilt: Dios et al. 2004 (see olb.turbulence.angle_of_arrival).
  Off-axis Gaussian encircled energy (Marcum Q): Marcum, RAND RM-753 (1950).
'''

import numpy as np

from ...results import Term
from ...assumptions import (Assumptions, BEAM_GAUSSIAN, REGIME_WEAK,
                            SPECTRUM_KOLMOGOROV)
from ...terminal import SMF, MMF, TipTilt, AO
from ...beam import free_space_radius, launch_curvature
from ...turbulence.gaussian_fried import gaussian_fried_parameter_profile
from ...turbulence.angle_of_arrival import wander_arrival_angle_variance
from ...turbulence.ao import apply_compensation
from ._common import (_smf_eta_max, _smf_coupling_efficiency, _effective_dr0,
                     _smf_static_term, SMF_OPTIMAL_A, SMF_DEEP_TURBULENCE_DR0)

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

    # Horizontal Gaussian-beam Fried parameter over the constant-Cn2 path. The
    # launch curvature f0 of a deliberately diverged beam enters the beam
    # parameters (Theta0 = 1 - L/f0), so a diverged beam gets its own r0.
    # This closes olb Gap 3. See olb.beam.launch_curvature.
    f0 = launch_curvature(w0, tx.transmitter.divergence_rad, wavelength)
    hs = np.linspace(0.0, L, int(n_grid))
    cn2_profile = np.full_like(hs, cn2)
    r0 = gaussian_fried_parameter_profile(hs, cn2_profile, w0, wavelength,
                                          path='terrestrial', f0=f0)

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
            "f0_m": float(f0),
            "n_comp_modes": residual.n_modes,
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

    # Diffraction focal spot radius (1/e^2), Gaussian approximation to the Airy.
    w_s = wavelength * f / (np.pi * (D / 2.0))

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

    # The received tip-tilt moves the spot centre off the core by dx = f*theta.
    # theta is the radial (2-axis) tip-tilt, so the per-axis spot offset has a
    # 1-sigma of sigma_d = f*sqrt(sigma2_theta/2). The coupled power is the
    # encircled energy of the displaced Gaussian spot inside the hard core (a
    # light bucket, NOT a mode overlap). See _mmf_encircled_efficiency.
    sigma2_theta, meta = _received_tiptilt_variance(scenario, n_grid=n_grid,
                                                    turbulence=turbulence)
    sigma_d = f * np.sqrt(sigma2_theta / 2.0)          # per-axis spot offset [m]
    eta_static = na_factor * float(_mmf_encircled_efficiency(0.0, w_s, a_core))
    static_db = -10.0 * np.log10(eta_static)

    def _loss_db(offset_m):
        eta = np.clip(na_factor * _mmf_encircled_efficiency(offset_m, w_s, a_core),
                      1e-300, None)
        return -10.0 * np.log10(eta)

    # Mean over the Rayleigh offset distribution (two i.i.d. Gaussian axes). A
    # quadrature over the radial offset gives the expected loss. The 8-sigma grid
    # covers the Rayleigh tail (the pdf there is exp(-32), negligible).
    if sigma_d > 0.0:
        dd = np.linspace(0.0, 8.0 * sigma_d, 2000)
        rayleigh = dd / sigma_d ** 2 * np.exp(-dd ** 2 / (2.0 * sigma_d ** 2))
        mean_db = float(np.trapz(_loss_db(dd) * rayleigh, dd))
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

    assumptions = Assumptions(
        beam_type=BEAM_GAUSSIAN,
        turbulence_regime=REGIME_WEAK,
        spectrum=SPECTRUM_KOLMOGOROV,
        validity="Multimode-fibre coupling (light bucket): the coupled power is the "
                 "encircled energy of the Gaussian focal spot (1/e^2 radius w_s = "
                 "lambda*f/(pi*(D/2))) inside the hard core disk of radius a_core, "
                 "offset by the received tip-tilt dx = f*theta. It is the "
                 "non-central chi-square CDF (Marcum Q), so it collects ALL the spot "
                 "power inside the core, not a mode overlap. The tip-tilt is the "
                 "weak beam wander (tracked by a tip-tilt or AO stage) plus the "
                 "receive mechanical jitter. The spot model assumes a uniform, "
                 "unobscured circular aperture. The numerical-aperture gate (when "
                 "set) is a flat power-transmission factor, NOT a re-truncated "
                 "aperture: it does not re-broaden the focal spot.",
    )
    if na_fibre is not None and na_optic > na_fibre:
        assumptions.flag(
            f"The focusing cone NA_optic={na_optic:.3f} exceeds the fibre "
            f"NA={na_fibre:.3f}; the fibre does not guide the steep rays. The "
            f"angular gate cuts the coupled power by {(-10.0 * np.log10(na_factor)):.2f} "
            "dB. Shorten nothing further, or use a larger-NA fibre."
        )
        assumptions.flag(
            "The NA gate (NA/NA_optic)^2 is an aperture-AREA fraction, so it assumes "
            "a uniformly illuminated pupil that fills the aperture. Here D is the "
            "ILLUMINATED diameter, not only the mechanical stop. An underfilled or "
            "Gaussian-apodized pupil carries less power in the steep marginal rays, "
            "so its true NA loss is SMALLER. Thus this gate is CONSERVATIVE "
            "(pessimistic) for an underfilled pupil, and exact for a filled uniform "
            "one. A distant-source receive beam is near uniform, so the assumption "
            "holds for a receive aperture."
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
             f"w_s={w_s * 1e6:.1f} um, static={static_db:.2f} dB"
             + (f", NA gate {(-10.0 * np.log10(na_factor)):.2f} dB"
                if na_factor < 1.0 else ""),
        meta={**meta, "detector": "MMF", "core_radius_m": a_core,
              "focal_length_m": f, "spot_radius_m": float(w_s),
              "eta_static": eta_static, "static_loss_db": float(static_db),
              "walkoff_mean_db": float(mean_db - static_db),
              "spot_offset_1sigma_m": float(sigma_d),
              "numerical_aperture": na_fibre, "na_optic": float(na_optic),
              "na_factor": na_factor,
              "na_gate_loss_db": float(-10.0 * np.log10(na_factor))},
        assumptions=assumptions,
    )


if __name__ == '__main__':
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
                         pointing_jitter_rad=jitter, detector=detector,
                         compensation=compensation or []),
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
    assert np.isclose(t_mmf.meta["eta_static"], 1.0 - np.exp(-2.0 * 1.12 ** 2), atol=1e-3)
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

    wo_big_q99 = wo_big.quantile_db(0.99)
    t_mmf_q99 = t_mmf.quantile_db(0.99)
    print(f"SMF walk-off (10 urad): mean {wo_big.mean_db:.3f} dB  "
          f"99% {float(wo_big_q99) if wo_big_q99 is not None else float('nan'):.3f} dB")
    print(f"MMF (25 um core, 5 urad): static {t_mmf.meta['static_loss_db']:.3f} dB  "
          f"mean {t_mmf.mean_db:.3f} dB  "
          f"99% {float(t_mmf_q99) if t_mmf_q99 is not None else float('nan'):.3f} dB")
    print("coupling terrestrial self-check passed")
