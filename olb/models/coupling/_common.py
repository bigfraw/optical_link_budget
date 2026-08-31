'''
Shared single-mode-fibre coupling physics, independent of the link.

Both the downlink coupling Term and the terrestrial coupling Term read the same
core single-mode-fibre physics. This module holds that shared set: the coupling
efficiency against the residual phase variance, the flat-wavefront eta_max from
the coupling parameter a, and the turbulence-off static mode-match Term. The
downlink module and the terrestrial module import from here. This is intra-package
reuse; it does not import another model family.

Two limits of one overlap physics set the coupling efficiency:
  small residual (sigma^2_res < SMF_SMALL_RESIDUAL_LIMIT): extended Marechal,
      eta = eta_max * exp(-sigma^2_res).
  large residual (sigma^2_res >= the limit): the Dikmelik-Davidson uncorrected
      single-mode-fibre coupling curve against the effective D/r0.

Sources:
  Noll residual variance: R. J. Noll, JOSA 66(3), 207 (1976). See
  olb.turbulence.ao.
  SMF coupling / Marechal: extended Marechal approximation.
  Uncorrected SMF coupling against D/r0: Y. Dikmelik and F. M. Davidson,
  "Fiber-coupling efficiency for free-space optical communication through
  atmospheric turbulence," Appl. Opt. 44(23), 4946-4952 (2005), DOI
  10.1364/AO.44.004946.
  eta_max(a): Shaklan and Roddier, Appl. Opt. 27, 2334 (1988), DOI
  10.1364/AO.27.002334.
'''

import numpy as np

from ...results import Term
from ...assumptions import (Assumptions, BEAM_GAUSSIAN, REGIME_NA, SPECTRUM_NA)
from ...turbulence.ao import NOLL_PISTON

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

    It peaks at eta_max=0.8145 near a=1.12, and it falls on both sides.

    Assumes a UNIFORMLY illuminated aperture and a FLAT (best-focus) wavefront.
    So it holds when the received spot overfills the aperture (the aperture reads
    the near-flat centre of the beam) and the receiver focuses for the incoming
    curvature (the optimal_focus case). A near-field terrestrial link inside the
    Rayleigh range can break BOTH: the received Gaussian tapers across the
    aperture, and the wavefront is curved. Curvature is removable by refocus. The
    residual taper error runs SAFE, because a Gaussian-into-Gaussian overlap can
    exceed the 0.8145 top-hat value, so this constant is then conservative. A
    curvature-aware, illumination-aware eta_max is the Gap-3 upgrade (see
    CLAUDE.md and olb.models.coupling.terrestrial). Source:
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


def smf_eta_defocused(a, c):
    '''
    Return the single-mode-fibre coupling eta with a DEFOCUS aberration.

    smf_eta_max_from_a assumes a FLAT (best-focus) wavefront. Put the fibre tip a
    distance dz_eff away from the true focus and the uniformly illuminated pupil
    carries a quadratic (defocus) phase. The mode-overlap integral of that pupil
    with the Gaussian fibre mode stays closed form: the a^2 of the Gaussian
    weight becomes the COMPLEX a^2 - i*c, so

        eta(a, c) = 2 a^2 | (1 - exp(-(a^2 - i c))) / (a^2 - i c) |^2,
        c = pi * dz_eff * (D/2)^2 / (lambda * f^2)   [rad, the edge defocus
                                                      aberration coefficient].

    At c = 0 this reduces EXACTLY to eta_max(a) = 2*((1-exp(-a^2))/a)^2. The
    result depends on |c| only, so the direction of the defocus does not matter.

    Sources: Shaklan and Roddier, Appl. Opt. 27, 2334 (1988), DOI
    10.1364/AO.27.002334 (the a parameter and the flat-wavefront overlap);
    Ruilier and Cassaing, JOSA A 18, 143 (2001), DOI 10.1364/JOSAA.18.000143
    (single-mode coupling with an aberrated pupil).

    Parameters:
        a : float or numpy.ndarray
            The coupling parameter a = pi*(D/2)*w_m/(lambda*f).
        c : float or numpy.ndarray
            The defocus aberration coefficient [rad] (see above). 0.0 is the
            flat-wavefront case.

    Returns:
        float or numpy.ndarray
            The coupling efficiency in (0, 0.8145].
    '''
    a = np.asarray(a, dtype=float)
    z = a ** 2 - 1j * np.asarray(c, dtype=float)
    return 2.0 * a ** 2 * np.abs((1.0 - np.exp(-z)) / z) ** 2


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
    # This Term is category "coupling" but reads NO decorated turbulence physics
    # (it is a deterministic mode-match loss, not a fade). So it opens no
    # collection context. It self-declares an "untraced: static optics"
    # provenance, so the Budget.check() untraced-Term guard stays quiet (the
    # guard flags a coupling Term with EMPTY provenance; see olb.results).
    assumptions = Assumptions(
        beam_type=BEAM_GAUSSIAN,
        turbulence_regime=REGIME_NA,
        spectrum=SPECTRUM_NA,
        validity="Turbulence off: static single-mode-fibre mode-match coupling only "
                 "(eta_max). No residual wavefront error and no fade.",
        provenance=["untraced: static optics"],
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


if __name__ == '__main__':
    import warnings

    from ...terminal import SMF

    lam = 1550e-9

    # --- SMF efficiency limits: small residual -> Marechal, large -> Dikmelik ---
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

    # --- eta_max from the coupling parameter a ------------------------------
    # The overlap curve peaks at 0.8145 near a=1.12 and falls on both sides.
    assert np.isclose(smf_eta_max_from_a(1.12), 0.8145, atol=1e-3), smf_eta_max_from_a(1.12)
    assert smf_eta_max_from_a(0.5) < 0.8145 and smf_eta_max_from_a(2.5) < 0.8145
    a_peak = np.linspace(0.9, 1.4, 51)
    assert np.isclose(a_peak[np.argmax(smf_eta_max_from_a(a_peak))], 1.12, atol=0.03)

    # --- the defocus-aberrated closed form ----------------------------------
    # c = 0 reduces EXACTLY to the flat-wavefront eta_max(a).
    for a_t in (0.5, 1.12, 2.5):
        assert np.isclose(smf_eta_defocused(a_t, 0.0), smf_eta_max_from_a(a_t),
                          rtol=1e-12), a_t
    # A defocus can only lose power, and more defocus loses more.
    cs = np.array([0.0, 1.0, 2.0, 4.0, 8.0])
    etas_c = smf_eta_defocused(1.12, cs)
    assert np.all(np.diff(etas_c) < 0.0), etas_c
    # The sign of c does not matter (|...|^2 of a conjugate pair).
    assert np.isclose(smf_eta_defocused(1.12, 3.95),
                      smf_eta_defocused(1.12, -3.95))
    # The report scenario (validation/defocus): a = 1.12, c = -3.95 rad from the
    # received-curvature defocus of a 5 km collimated terrestrial link.
    assert np.isclose(float(smf_eta_defocused(1.12, -3.95)), 0.215, atol=5e-4), \
        float(smf_eta_defocused(1.12, -3.95))

    # --- _smf_eta_max from the detector optics ------------------------------
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

    # --- static Term: deterministic, NOT mean-only, static mode-match floor --
    st = _smf_static_term(0.8145)
    assert np.isclose(st.mean_db, -10.0 * np.log10(0.8145))
    assert st.quantile is None and not st.mean_only and st.meta["model"] == "static"
    # The static Term is category "coupling" with no traced physics, so it
    # self-declares the "untraced: static optics" provenance that satisfies the
    # Budget.check() untraced-Term guard (olb.results).
    from ...results import Budget
    assert st.assumptions.provenance == ["untraced: static optics"], \
        st.assumptions.provenance
    with warnings.catch_warnings():
        warnings.simplefilter("error")   # the guard must NOT warn on this Term
        guard = [(n, r) for n, r in Budget([st]).check(warn=True)
                 if "did not open the assumption collection context" in r]
    assert guard == [], guard

    # --- _effective_dr0 inverts the piston-removed Noll relation ------------
    assert np.isclose(_effective_dr0(NOLL_PISTON), 1.0)

    print(f"eta_max(a=1.12) = {float(smf_eta_max_from_a(1.12)):.4f}")
    print(f"static SMF coupling floor = {float(st.mean_db):.3f} dB")
    print("coupling _common self-check passed")
