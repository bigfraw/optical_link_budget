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

from ..results import Term
from ..assumptions import (Assumptions, BEAM_PLANE_WAVE, REGIME_WEAK,
                           SPECTRUM_KOLMOGOROV)
from ..terminal import Aperture, SMF
from ..turbulence.profiles import DEFAULT_HS, default_cn2_profile
from ..turbulence.ao import (plane_wave_fried_parameter, apply_compensation,
                            NOLL_PISTON)
from ..links.downlink import downlink_scintillation_term

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
                 "residual. The fade uses the aperture-averaged lognormal "
                 "scintillation as a stand-in for the residual-coupling "
                 "fluctuation. Weak fluctuation: sigma2_I < 0.25.",
    )
    # Carry over any scintillation-side flag (weak-fluctuation, obscuration).
    if scint.assumptions is not None:
        for reason in scint.assumptions.violations:
            assumptions.flag(reason)
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


def rx_coupling_term(scenario, geometry, *, hs=None, cn2_profile=None):
    '''
    Build the ONE receive-coupling Term of a downlink receive terminal.

    Read scenario.rx_terminal. Dispatch on the detector type. An Aperture
    detector reuses the downlink aperture-averaged scintillation, so it is parity
    with the plain downlink. An SMF detector computes the fibre-coupling loss from
    the residual wavefront that the compensation stack leaves, and adds the
    scintillation fade.

    Parameters:
        scenario : Scenario
            Reads rx_terminal, link.wavelength_m, and the site Cn2 profile.
        geometry : CircularOrbit or TLEPass
            Reads elevation_deg. A scalar elevation gives a scalar Term.
        hs : numpy.ndarray, optional
            Heights above the ground station [m]. Defaults to DEFAULT_HS.
        cn2_profile : numpy.ndarray, optional
            Zenith Cn2(h) profile. Defaults to the site profile.

    Returns:
        Term
            category="coupling".

    Raises:
        ValueError
            If rx_terminal is None or has no detector, or the detector type is
            unknown.
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
        return _smf_term(scenario, geometry, hs=hs, cn2_profile=cn2_profile)
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
        return rx_coupling_term(scn, geom, hs=hs,
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

    # --- SMF, no correction: large coupling loss ---------------------------
    scn_smf = _downlink(Terminal(aperture_m=0.7, wavelength_m=lam, detector=SMF()))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        t_smf = build(scn_smf)
    assert t_smf.name == "receive coupling (SMF)"
    assert t_smf.meta["coupling_loss_db"] > 3.0         # 0.7 m, no AO -> big loss

    # --- SMF with AO: much lower coupling loss ------------------------------
    scn_ao = _downlink(Terminal(aperture_m=0.7, wavelength_m=lam, detector=SMF(),
                                compensation=[TipTilt(), AO(200)]))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        t_ao = build(scn_ao)
    assert t_ao.meta["coupling_loss_db"] < t_smf.meta["coupling_loss_db"]
    assert t_ao.meta["eta"] > t_smf.meta["eta"]

    # The SMF term is stochastic with a closed-form quantile deeper than the mean.
    rng = np.random.default_rng(0)
    assert t_ao.stochastic and t_ao.quantile is not None
    assert t_ao.quantile_db(0.99) > t_ao.mean_db
    draws = t_ao.sample_db(50_000, rng)
    assert abs(draws.mean() - t_ao.mean_db) < 0.05, (draws.mean(), t_ao.mean_db)

    # An elevation sweep broadcasts.
    sweep = CircularOrbit(600e3, np.array([40.0, 60.0, 90.0]))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        t_sweep = rx_coupling_term(scn_ao, sweep, hs=hs,
                                   cn2_profile=default_cn2_profile(scn_ao.channel.site, hs))
    assert np.shape(t_sweep.mean_db) == (3,)
    assert t_sweep.sample_db(100, rng).shape == (100, 3)

    print(f"aperture coupling mean = {float(t_ap.mean_db):.4f} dB "
          f"(= scintillation {float(t_scint.mean_db):.4f} dB)")
    print(f"SMF no-AO:  eta={t_smf.meta['eta']:.4f}  "
          f"coupling loss={t_smf.meta['coupling_loss_db']:.2f} dB  "
          f"D/r0_eff={t_smf.meta['effective_D_over_r0']:.2f}")
    print(f"SMF +AO200: eta={t_ao.meta['eta']:.4f}  "
          f"coupling loss={t_ao.meta['coupling_loss_db']:.2f} dB  "
          f"sigma2_res={t_ao.meta['sigma2_res']:.4f}")
    print("self-check passed")
