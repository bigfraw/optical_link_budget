'''
Fidelity-1 FAST models: downlink fibre coupling, and the pre-compensated uplink.

This module holds two FAST-driven Terms. `smf_fast_term` gives the downlink
single-mode-fibre coupling. `uplink_fast_term` gives the turbulence penalty of a
pre-compensated (downlink-beacon plus adaptive-optics) uplink. The two share the
FAST loader, the Cn2 layering, and the adaptive-optics mapping below.

The analytic mean-only (Dikmelik/Marechal) model does NOT compute the true fibre-
mode overlap, and it carries no fade. This module does both. It drives FAST (the
`fast-aosim` package) to get the downlink single-mode-fibre coupling as the
coherent overlap of the turbulent aperture field with the back-projected fibre
mode:

    eta(t) = |integral (Aperture * M_fibre) * exp(chi + i*phi) dA| ^ 2
             / |integral (Aperture * M_fibre) dA| ^ 2

FAST propagates Monte-Carlo phase screens (phase) with an aperture-averaged
log-normal scintillation (log-amplitude), and forms the mode overlap directly. So
this is a true modal-coupling metric, not a Strehl proxy. See
olb.models.coupling.downlink for the lower-fidelity models and the fidelity ladder.

FAST is an OPTIONAL dependency (GPLv3). This module imports it lazily. Install it
with `pip install fast-aosim`.

Mapping FAST -> the Term:
    floor_db = -link_budget['smf_coupling']    # static mode-match loss (aperture
                                               # + obscuration; ~0.9 dB unobscured)
    loss(t)  = floor_db - result.dB_rel(t)     # dB_rel is the turbulence penalty
The static floor lives in FAST's diffraction limit, not in the normalised power
result.dB_rel, so the two add. The Term stores the per-sample loss, so it gives an
empirical mean, an empirical quantile, and a resampling sampler.

FIRST-CUT LIMITS (flagged in the Term assumptions):
    - Scalar elevation only. An elevation array needs one FAST run per elevation.
    - Point-ahead is off (DTHETA=0): the up/down anisoplanatism of a moving
      satellite is not modelled. Pass fast_params={'DTHETA': [x, y]} in arcsec.
    - NOAO low-order (tilt) accuracy depends on the FAST grid (NPXLS). The auto
      grid may undersample tilt. Pass fast_params={'NPXLS': ...} for production.
'''

import logging

import numpy as np

from ..results import Term
from ..assumptions import (BEAM_GAUSSIAN, REGIME_WEAK,
                            SPECTRUM_KOLMOGOROV, SPECTRUM_VON_KARMAN,
                            trace_assumptions)
from ..terminal import TipTilt, AO
from ..turbulence.profiles import DEFAULT_HS, default_cn2_profile
from ..turbulence.plane_wave_scintillation import plane_wave_scintillation_index
from ..turbulence.andrews.scintillation import rytov_weak, LOGNORMAL_PDF_LIMIT


def _load_fast():
    '''Import fast-aosim lazily, with a helpful error when it is absent.'''
    try:
        import fast
        return fast
    except ImportError as e:
        raise ImportError(
            "the fidelity-1 SMF coupling needs the `fast-aosim` package. "
            "Run `pip install fast-aosim`, or use the analytic mean-only model "
            "(downlink_coupling_term smf_fidelity='mean'), which carries no fade."
        ) from e


def _cn2_layers(cn2_zenith, hs):
    '''
    Integrated Cn2 per layer [m^1/3] from a zenith Cn2(h) profile.

    FAST wants the integrated Cn2 per layer at ZENITH (it applies the airmass
    itself). So multiply the profile by the layer thickness.
    '''
    return np.asarray(cn2_zenith, float) * np.gradient(np.asarray(hs, float))


def _ao_params(compensation):
    '''
    Map the olb compensation stack to FAST adaptive-optics parameters.

    olb AO(n_modes) removes the first n_modes Zernike (Noll) modes, and TipTilt
    removes the first three (piston, tip, tilt). FAST corrects modally up to ZMAX
    Zernikes when MODAL is set. So an AO stage maps to AO_MODE="AO", MODAL=True,
    ZMAX=max(n_modes); a tip-tilt-only stack maps to AO_MODE="TT" (Zmax=3); an
    empty stack maps to AO_MODE="NOAO".

    Returns:
        dict
            The AO-related FAST parameters (AO_MODE, and MODAL/ZMAX for AO).
    '''
    ao_modes = [c.n_modes for c in compensation if isinstance(c, AO)]
    if ao_modes:
        return {"AO_MODE": "AO", "MODAL": True, "ZMAX": int(max(ao_modes))}
    if any(isinstance(c, TipTilt) for c in compensation):
        return {"AO_MODE": "TT"}
    return {"AO_MODE": "NOAO"}


def _spectrum_label(params):
    '''
    Read the spectrum label and the scale note from the RESOLVED FAST scales.

    FAST is a von Karman engine (it always calls turb_powerspectrum_vonKarman).
    olb sets the scales itself (L0=inf, l0=1e-6), so the spectrum is Kolmogorov
    by our own choice, not inherited from the FAST conf.py. A finite L0 (or a
    large l0) from fast_params makes it a true von Karman spectrum. So read the
    label from the resolved scales, not from a fixed constant. The two spectra
    are Ch. 3 of Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196.

    Parameters:
        params : dict
            The resolved FAST parameters. Reads L0 and l0.

    Returns:
        tuple
            (spectrum, scale_note, L0, l0).
    '''
    L0 = float(params["L0"])
    l0 = float(params["l0"])
    kolmogorov = np.isinf(L0) and l0 <= 1e-3
    spectrum = SPECTRUM_KOLMOGOROV if kolmogorov else SPECTRUM_VON_KARMAN
    scale_note = (
        "The outer scale is infinite and the inner scale is 1 um (the Kolmogorov "
        "limit); pass fast_params={'L0': ...} [m] or {'l0': ...} for a finite "
        "von Karman scale." if kolmogorov else
        f"von Karman spectrum with outer scale L0={L0:g} m and inner scale "
        f"l0={l0:g} m.")
    return spectrum, scale_note, L0, l0


def smf_fast_term(scenario, geometry, *, hs=None, cn2_profile=None,
                  n_samples=1000, fast_params=None):
    '''
    Fidelity-1 SMF receive-coupling Term for a downlink, computed by FAST.

    Parameters:
        scenario : SpaceScenario or TerrestrialScenario
            The downlink case. Reads rx_terminal (aperture, obscuration,
            wavelength, compensation), the site (Cn2, wind), and the orbit.
        geometry : CircularOrbit or TLEPass
            The link geometry. SCALAR elevation only in this first cut.
        hs : numpy.ndarray, optional
            Zenith height grid [m]. Defaults to DEFAULT_HS.
        cn2_profile : numpy.ndarray, optional
            Zenith Cn2(h) profile matching hs. Defaults to the site profile.
        n_samples : int
            FAST Monte Carlo draws (NITER).
        fast_params : dict, optional
            Extra FAST parameters, merged last (for example {'NPXLS': 128,
            'DTHETA': [4, 0], 'L0': 25.0}). Overrides the mapped defaults.

    Returns:
        Term
            name="receive coupling (SMF)", category="coupling". Carries an
            empirical mean, quantile, and sampler from the FAST draws.

    Raises:
        ImportError
            If fast-aosim is not installed.
        ValueError
            If the geometry elevation is not scalar.
    '''
    fast = _load_fast()

    elev = np.asarray(geometry.elevation_deg, dtype=float)
    if elev.ndim != 0:
        raise ValueError(
            "smf_fast_term takes a scalar elevation in this first cut. Loop over "
            "elevations and build one Term each."
        )
    elev = float(elev)

    rx = scenario.rx_terminal
    D = rx.aperture_m
    obsc_ratio = rx.obscuration_ratio
    hs = DEFAULT_HS if hs is None else hs
    if cn2_profile is None:
        cn2_profile = default_cn2_profile(scenario.channel.site, hs)
    hs = np.asarray(hs, dtype=float)
    cn2_layer = _cn2_layers(cn2_profile, hs)
    wind = float(scenario.channel.site.wind_rms_m_s)
    ao_params = _ao_params(rx.compensation)   # AO_MODE, and MODAL/ZMAX for AO

    params = dict(
        WVL=rx.wavelength_m,
        D_GROUND=D,
        OBSC_GROUND=obsc_ratio * D,          # FAST wants the obscuration DIAMETER
        SMF=True,
        PROP_DIR="down",                     # SMF coupling is a receive-side quantity
        W0="opt",                            # FAST optimises the fibre-mode size
        H_SAT=scenario.channel.altitude_m,
        ZENITH_ANGLE=90.0 - elev,
        H_TURB=hs,
        CN2_TURB=cn2_layer,
        WIND_SPD=wind * np.ones_like(hs),
        WIND_DIR=np.zeros_like(hs),
        DTHETA=[0, 0],                       # first cut: no point-ahead
        SUBHARM=True,                        # capture low-order tilt (see below)
        L0=np.inf,                           # outer scale: infinite -> Kolmogorov
        l0=1e-6,                             # inner scale: 1 um, below any optical scale
        NITER=int(n_samples),
        LOGLEVEL="ERROR",
    )
    params.update(ao_params)                 # AO_MODE (+ MODAL/ZMAX from n_modes)
    if fast_params:
        params.update(fast_params)           # a finite L0/l0 here makes it von Karman

    # Quiet the FAST logger (it logs each init step at INFO) and its tqdm progress
    # bar. run() wraps the chunk loop in tqdm with no disable argument, so replace
    # the name FAST resolves (fast.fast.tqdm) with a pass-through for the run.
    fast_logger = logging.getLogger("fast")
    old_level = fast_logger.level
    old_tqdm = fast.fast.tqdm
    fast_logger.setLevel(logging.ERROR)
    fast.fast.tqdm = lambda iterable=None, *a, **k: iterable
    try:
        sim = fast.Fast(params)
        floor_db = -float(sim.compute_link_budget()["smf_coupling"])
        result = sim.run()
    finally:
        fast_logger.setLevel(old_level)
        fast.fast.tqdm = old_tqdm

    # dB_rel is the turbulence coupling penalty (10*log10 of power normalised to
    # the diffraction limit). The static mode-match floor is separate, so add it.
    loss_db = floor_db - np.asarray(result.dB_rel, dtype=float)
    mean_db = float(loss_db.mean())

    # AMPLITUDE-regime check. FAST models the phase with real Monte-Carlo screens
    # (the phase-driven modal-coupling fade is fidelity-1), but it models the log-
    # AMPLITUDE as an aperture-averaged log-normal, which only holds in the weak
    # fluctuation regime. The regime is set by the plane-wave scintillation index
    # sigma2_I (NOT by result.scintillation_index, which is the COUPLED-power index
    # -- for a no-AO fibre that is dominated by phase speckle and is routinely > 1
    # even when the amplitude is weak). So flag on the amplitude sigma2_I only.
    # AMPLITUDE-regime physics. Trace it so the weak-regime check on
    # plane_wave_scintillation_index registers automatically. FAST models the
    # phase with real Monte-Carlo screens (the phase-driven modal-coupling fade is
    # fidelity-1), but it models the log-AMPLITUDE as an aperture-averaged
    # log-normal, which only holds in the weak fluctuation regime. The regime is
    # set by the plane-wave scintillation index sigma2_I (NOT by
    # result.scintillation_index, which is the COUPLED-power index -- for a no-AO
    # fibre that is dominated by phase speckle and is routinely > 1 even when the
    # amplitude is weak).
    with trace_assumptions() as trace:
        sigma2_I_amp = float(plane_wave_scintillation_index(
            elev, rx.wavelength_m, hs, cn2_profile))

    def quantile(p):
        return float(np.quantile(loss_db, p))

    def sampler(n, rng):
        return rng.choice(loss_db, size=n, replace=True)

    # The spectrum label follows the RESOLVED scales (see _spectrum_label).
    spectrum, scale_note, L0, l0 = _spectrum_label(params)
    assumptions = trace.merge(
        beam_type=BEAM_GAUSSIAN,
        turbulence_regime=REGIME_WEAK,
        spectrum=spectrum,
        validity="Fidelity-1 single-mode-fibre coupling: the true LP01 Gaussian-"
                 "mode overlap under turbulence, from FAST (fast-aosim). The "
                 "received field is the aperture times the fibre mode, propagated "
                 "through Monte-Carlo phase screens with an aperture-averaged log-"
                 "normal scintillation. eta_max and the fibre-mode size come from "
                 "FAST (W0='opt'). Low-order tilt is captured with subharmonics "
                 "(SUBHARM=True); without them the small auto grid undersamples the "
                 "tilt and understates the loss by several dB. " + scale_note +
                 " Weak-to-moderate fluctuation (log-normal scintillation).",
    )
    # The core FAST coupling is an EXTERNAL Monte Carlo (fast-aosim), not an
    # @assumes physics function, so the trace cannot see it. Declare it honestly,
    # so the provenance is not empty and the untraced guard is satisfied.
    assumptions.provenance.append("untraced: fast-aosim")
    # First-cut limitation: no point-ahead. A scenario-level fact the traced
    # physics never sees, so it stays a factory flag. Point-ahead anisoplanatism
    # only matters for a PRE-COMPENSATED beam (the corrected modes decorrelate
    # over the point-ahead angle); a plain receive coupling never uses it, so the
    # flag only trips when the scenario carries a precompensation source.
    if getattr(scenario, "precompensation", None) is not None:
        assumptions.flag(
            "Point-ahead is off (DTHETA=0): the up-leg / down-leg anisoplanatism "
            "of a moving satellite is not modelled. Pass fast_params={'DTHETA': "
            "[x, y]} in arcsec to include it.",
            source="factory:models.fast"
        )
    # AMPLITUDE saturation: only here does FAST's log-normal scintillation break.
    # A large coupled-power fade (deep 99% tail) does NOT trip this -- that fade is
    # phase-driven and modelled correctly by the screens.
    # TWO distinct amplitude tests, kept separate (Conflict C-05): the REGIME
    # boundary (does the amplitude log-normal saturate? sigma_R^2 = 1) and the
    # tighter lognormal-PDF house rule (is the amplitude tail optimistic?
    # sigma2_I = 0.25). The amplitude index is a plane wave, so Lambda is None.
    amp_regime = rytov_weak(sigma2_I_amp)
    # The REGIME hard-flag (sigma2_I >= WEAK_REGIME_LIMIT) is now produced by the
    # traced weak-regime check on plane_wave_scintillation_index (an equivalent
    # violation), so no factory flag is built for it. The tighter lognormal-PDF
    # caution (>= 0.25 but still in the weak regime) has no traced check, so it
    # stays a factory flag.
    if amp_regime != 'hard' and sigma2_I_amp >= LOGNORMAL_PDF_LIMIT:
        assumptions.flag(
            f"Plane-wave amplitude scintillation sigma2_I={sigma2_I_amp:.2f} is past "
            f"the lognormal-PDF limit {LOGNORMAL_PDF_LIMIT} (but within the weak "
            "regime); the amplitude log-normal TAIL is optimistic. The phase-driven "
            "coupling fade dominates and is modelled by the screens.",
            source="factory:models.fast"
        )

    zmax = params.get("ZMAX")
    ao_desc = params["AO_MODE"] + (f"(ZMAX={zmax})" if zmax is not None else "")
    note = (f"FAST fidelity-1 modal coupling, AO={ao_desc}, "
            f"floor={floor_db:.2f} dB, NITER={params['NITER']}")
    return Term(
        name="receive coupling (SMF)",
        category="coupling",
        mean_db=mean_db,
        sampler=sampler,
        quantile=quantile,
        note=note,
        meta={
            "detector": "SMF",
            "model": "fast-modal",
            "ao_mode": params["AO_MODE"],
            "zmax": params.get("ZMAX"),
            "floor_db": floor_db,
            "subharmonics": bool(params.get("SUBHARM", False)),
            "npxls": int(getattr(sim, "Npxls", 0)),
            "coupled_scintillation_index": float(result.scintillation_index),
            "amplitude_sigma2_I": sigma2_I_amp,
            "amplitude_rytov_regime": amp_regime,
            "amplitude_regime_weak": amp_regime != 'hard',
            "r0_los_m": float(getattr(sim, "r0_los", np.nan)),
            "L0_m": L0,
            "l0_m": l0,
            "n_samples": int(params["NITER"]),
        },
        assumptions=assumptions,
    )


# Arcsec per radian (180 / pi * 3600). FAST takes the point-ahead offset DTHETA
# in arcsec, and olb keeps every angle in radians.
ARCSEC_PER_RAD = 206265.0


def uplink_fast_term(scenario, geometry, *, hs=None, cn2_profile=None,
                     n_samples=1000, fast_params=None):
    '''
    Fidelity-1 turbulence Term for a PRE-COMPENSATED uplink, computed by FAST.

    This is the model of record for a downlink-beacon plus adaptive-optics
    uplink. It replaces the two analytic phase Terms of olb.links.uplink, which
    give the mean Strehl loss only and carry no scintillation and no fade.

    HOW FAST GIVES AN UPLINK NUMBER. The FAST Monte Carlo is direction-agnostic.
    `compute_detector` overlaps the ground-pupil field (a hard circular aperture
    times a Gaussian mode) with the residual phase and a log-normal
    log-amplitude, and it normalises to the vacuum overlap. By reciprocity that
    normalised overlap is the uplink on-axis flux at the satellite, when the
    launch mode is that Gaussian and the residual phase is the pre-compensation
    error at the pupil. Source: J. H. Shapiro, JOSA 61(4), 492 (1971),
    DOI 10.1364/JOSA.61.000492. The FAST parameter PROP_DIR changes the analytic
    `compute_link_budget` only, which olb does not read here. FAST is the
    published tool for this case: O. J. D. Farley and others, Opt. Express
    30(13), 23050 (2022), DOI 10.1364/OE.458659.

    POINT-AHEAD. DTHETA is the point-ahead offset in arcsec, [x, y].
    `ao_power_spectra.G_AO_PAOLA` turns it into a per-layer displacement
    h * theta / 206265, and it builds the anisoplanatism-plus-servo residual
    filter on the CORRECTED (low-order) mask only. So DTHETA models exactly the
    point-ahead decorrelation of the corrected modes.

    WHAT THIS TERM DOES NOT HOLD. `result.dB_rel` is the received power in dB
    relative to the diffraction limit (no turbulence), so this Term is the PURE
    turbulence penalty. It holds no static loss. Do not double-count: the
    free-space spread and the receive capture are in the geometric Term, the
    launch truncation is in tx_gaussian_efficiency_term, and the mechanical
    tracking jitter is in the standalone pointing Term.

    PHASE-ONLY PRE-COMPENSATION. The adaptive-optics filter touches the phase
    only. The log-amplitude stays an aperture-and-mode-filtered log-normal, so
    the amplitude scintillation is NOT corrected. That is the physics: a
    wavefront correction does not remove the amplitude fluctuation. The
    log-normal holds in the weak-fluctuation regime only, so this Term gates on
    the plane-wave scintillation index, the same as smf_fast_term.

    Parameters:
        scenario : SpaceScenario
            The uplink case. Reads tx_terminal (the ground terminal: aperture,
            obscuration, wavelength, Transmitter waist, compensation),
            rx_terminal (the satellite aperture), the site (Cn2, wind), and the
            orbit altitude.
        geometry : CircularOrbit or TLEPass
            The link geometry. SCALAR elevation only. Reads elevation_deg and
            point_ahead_rad.
        hs : numpy.ndarray, optional
            Zenith height grid [m]. Defaults to DEFAULT_HS.
        cn2_profile : numpy.ndarray, optional
            Zenith Cn2(h) profile matching hs. Defaults to the site profile.
        n_samples : int
            FAST Monte Carlo draws (NITER).
        fast_params : dict, optional
            Extra FAST parameters, merged last (for example {'NPXLS': 128,
            'DTHETA': [0, 0], 'L0': 25.0}). Overrides the mapped defaults.

    Returns:
        Term
            name="turbulence (pre-compensated, FAST)", category="turbulence".
            Carries an empirical mean, quantile, and sampler from the draws.

    Raises:
        ImportError
            If fast-aosim is not installed.
        ValueError
            If the scenario is not an uplink, if the elevation is not scalar, if
            the ground terminal has no Transmitter waist, or if the ground
            terminal has no adaptive-optics stage.
    '''
    fast = _load_fast()

    if getattr(scenario, "direction", None) != "uplink":
        raise ValueError(
            "uplink_fast_term takes an uplink SpaceScenario. Set "
            "direction='uplink', or use smf_fast_term for the downlink."
        )

    elev = np.asarray(geometry.elevation_deg, dtype=float)
    if elev.ndim != 0:
        raise ValueError(
            "uplink_fast_term takes a scalar elevation in this first cut. Loop "
            "over elevations and build one Term each."
        )
    elev = float(elev)

    tx = scenario.tx_terminal
    if tx.transmitter is None or tx.transmitter.waist_m is None:
        raise ValueError(
            "uplink_fast_term needs a launch beam. Give the ground terminal a "
            "Transmitter with waist_m."
        )
    if not any(isinstance(c, AO) for c in tx.compensation):
        raise ValueError(
            "uplink_fast_term models a PRE-COMPENSATED uplink, so the ground "
            "terminal must have an AO(n_modes) stage in its compensation "
            "stack. With no adaptive optics the FAST run degenerates to an "
            "uncorrected beam that FAST does not model correctly for an "
            "uplink. Use olb.links.uplink.uplink_turbulence_term (the "
            "coupled-flux Monte Carlo) for the uncorrected route."
        )

    # Bistatic override (see olb.terminal.Transmitter): the Transmitter values
    # win when they are set, else the owning Terminal values apply.
    D = (tx.transmitter.aperture_m if tx.transmitter.aperture_m is not None
         else tx.aperture_m)
    obsc_ratio = (tx.transmitter.obscuration_ratio
                  if tx.transmitter.obscuration_ratio is not None
                  else tx.obscuration_ratio)

    rx = scenario.rx_terminal
    hs = DEFAULT_HS if hs is None else hs
    if cn2_profile is None:
        cn2_profile = default_cn2_profile(scenario.channel.site, hs)
    hs = np.asarray(hs, dtype=float)
    cn2_layer = _cn2_layers(cn2_profile, hs)
    wind = float(scenario.channel.site.wind_rms_m_s)
    ao_params = _ao_params(tx.compensation)   # AO_MODE, MODAL, ZMAX

    params = dict(
        WVL=tx.wavelength_m,
        D_GROUND=D,
        OBSC_GROUND=obsc_ratio * D,          # FAST wants the obscuration DIAMETER
        W0=float(tx.transmitter.waist_m),    # the numeric launch waist, not "opt"
        PROP_DIR="up",                       # analytic budget only; the MC is agnostic
        SMF=False,                           # the satellite is not a fibre receiver
        D_SAT=rx.aperture_m,
        OBSC_SAT=rx.obscuration_ratio * rx.aperture_m,
        H_SAT=scenario.channel.altitude_m,
        ZENITH_ANGLE=90.0 - elev,
        H_TURB=hs,
        CN2_TURB=cn2_layer,
        WIND_SPD=wind * np.ones_like(hs),
        WIND_DIR=np.zeros_like(hs),
        # Point-ahead: rad -> arcsec, along x. This drives the decorrelation of
        # the corrected modes (G_AO_PAOLA, see the docstring).
        DTHETA=[float(geometry.point_ahead_rad) * ARCSEC_PER_RAD, 0.0],
        SUBHARM=True,                        # capture low-order tilt
        L0=np.inf,                           # outer scale: infinite -> Kolmogorov
        l0=1e-6,                             # inner scale: 1 um, below any optical scale
        NITER=int(n_samples),
        LOGLEVEL="ERROR",
    )
    params.update(ao_params)
    if fast_params:
        params.update(fast_params)           # a finite L0/l0 here makes it von Karman

    # Quiet the FAST logger and its tqdm progress bar, the same as smf_fast_term.
    fast_logger = logging.getLogger("fast")
    old_level = fast_logger.level
    old_tqdm = fast.fast.tqdm
    fast_logger.setLevel(logging.ERROR)
    fast.fast.tqdm = lambda iterable=None, *a, **k: iterable
    try:
        sim = fast.Fast(params)
        result = sim.run()
    finally:
        fast_logger.setLevel(old_level)
        fast.fast.tqdm = old_tqdm

    # dB_rel is negative for a turbulence penalty, and olb loss is positive.
    # No static floor is added: this Term is the PURE turbulence penalty.
    loss_db = -np.asarray(result.dB_rel, dtype=float)
    mean_db = float(loss_db.mean())

    # AMPLITUDE-regime physics. The adaptive-optics filter does not touch the
    # log-amplitude, and that log-normal holds in the weak-fluctuation regime
    # only. Gate on the plane-wave index, NOT on result.scintillation_index
    # (which is the COUPLED-power index and is phase-dominated). Trace it so the
    # weak-regime check on plane_wave_scintillation_index registers automatically.
    with trace_assumptions() as trace:
        sigma2_I_amp = float(plane_wave_scintillation_index(
            elev, tx.wavelength_m, hs, cn2_profile))

    def quantile(p):
        return float(np.quantile(loss_db, p))

    def sampler(n, rng):
        return rng.choice(loss_db, size=n, replace=True)

    spectrum, scale_note, L0, l0 = _spectrum_label(params)
    assumptions = trace.merge(
        beam_type=BEAM_GAUSSIAN,
        turbulence_regime=REGIME_WEAK,
        spectrum=spectrum,
        validity="Fidelity-1 pre-compensated uplink turbulence penalty from FAST "
                 "(fast-aosim). Source: O. J. D. Farley and others, Opt. Express "
                 "30(13), 23050 (2022), DOI 10.1364/OE.458659. The FAST Monte "
                 "Carlo overlaps the ground-pupil field with the residual phase "
                 "and a log-normal log-amplitude, and it normalises to the vacuum "
                 "overlap. By reciprocity that overlap is the uplink on-axis flux "
                 "at the satellite. Source: J. H. Shapiro, JOSA 61(4), 492 "
                 "(1971), DOI 10.1364/JOSA.61.000492. The point-ahead offset "
                 "DTHETA decorrelates the CORRECTED modes only "
                 "(ao_power_spectra.G_AO_PAOLA). "
                 "PHASE-ONLY PRE-COMPENSATION: the adaptive-optics filter does "
                 "not touch the log-amplitude, so the amplitude scintillation "
                 "stays uncorrected. That log-normal holds in the weak "
                 "fluctuation regime only. "
                 "FAST servo and wavefront-sensor defaults in force: DSUBAP=0.02 "
                 "m, TLOOP=0.001 s, TEXP=0.001 s, ALIAS=True, NOISE=0. Pass "
                 "fast_params to change them. " + scale_note + " "
                 "This Term is the PURE turbulence penalty. It holds no static "
                 "loss. The free-space spread and the receive capture are in the "
                 "geometric Term, the launch truncation is in "
                 "tx_gaussian_efficiency_term, and the mechanical tracking "
                 "jitter is in the standalone pointing Term. Do not "
                 "double-count them.",
    )
    # The core FAST coupling is an EXTERNAL Monte Carlo (fast-aosim), not an
    # @assumes physics function, so the trace cannot see it. Declare it honestly,
    # so the provenance is not empty and the untraced guard is satisfied.
    assumptions.provenance.append("untraced: fast-aosim")
    # TWO distinct amplitude tests, kept separate (Conflict C-05): the REGIME
    # boundary (does the amplitude log-normal saturate? sigma_R^2 = 1) and the
    # tighter lognormal-PDF house rule (is the amplitude tail optimistic?
    # sigma2_I = 0.25). The amplitude index is a plane wave, so Lambda is None.
    amp_regime = rytov_weak(sigma2_I_amp)
    # The REGIME hard-flag (sigma2_I >= WEAK_REGIME_LIMIT) is now produced by the
    # traced weak-regime check on plane_wave_scintillation_index (an equivalent
    # violation), so no factory flag is built for it. The tighter lognormal-PDF
    # caution (>= 0.25 but still in the weak regime) has no traced check, so it
    # stays a factory flag.
    if amp_regime != 'hard' and sigma2_I_amp >= LOGNORMAL_PDF_LIMIT:
        assumptions.flag(
            f"Plane-wave amplitude scintillation sigma2_I={sigma2_I_amp:.2f} is past "
            f"the lognormal-PDF limit {LOGNORMAL_PDF_LIMIT} (but within the weak "
            "regime); the amplitude log-normal TAIL is optimistic. The phase-driven "
            "fade dominates and is modelled by the screens.",
            source="factory:models.fast"
        )

    zmax = params.get("ZMAX")
    dtheta_arcsec = float(params["DTHETA"][0])
    note = (f"FAST pre-compensated uplink, AO(ZMAX={zmax}), "
            f"point-ahead {dtheta_arcsec:.2f} arcsec, NITER={params['NITER']}")
    return Term(
        name="turbulence (pre-compensated, FAST)",
        category="turbulence",
        mean_db=mean_db,
        sampler=sampler,
        quantile=quantile,
        note=note,
        meta={
            "model": "fast-precomp-uplink",
            "ao_mode": params["AO_MODE"],
            "zmax": zmax,
            "dtheta_arcsec": dtheta_arcsec,
            "coupled_scintillation_index": float(result.scintillation_index),
            "amplitude_sigma2_I": sigma2_I_amp,
            "amplitude_rytov_regime": amp_regime,
            "amplitude_regime_weak": amp_regime != 'hard',
            "r0_los_m": float(getattr(sim, "r0_los", np.nan)),
            "L0_m": L0,
            "l0_m": l0,
            "n_samples": int(params["NITER"]),
        },
        assumptions=assumptions,
    )


if __name__ == '__main__':
    import warnings

    from ..scenario import SpaceScenario, Channel, Site
    from ..geometry import CircularOrbit
    from ..terminal import Terminal, Transmitter, SMF, TipTilt, AO

    try:
        _load_fast()
    except ImportError as e:
        print(f"fast-aosim not installed ({e.__class__.__name__}); skipping the "
              "fidelity-1 self-check.")
        raise SystemExit

    lam = 1550e-9

    def _downlink(ground):
        return SpaceScenario(
            ground=ground,
            space=Terminal(aperture_m=0.05, wavelength_m=lam,
                           transmitter=Transmitter(waist_m=0.035)),
            direction="downlink",
            channel=Channel(site=Site(cn2_ground=1.7e-14), altitude_m=1500e3))

    geom = CircularOrbit(1500e3, elevation_deg=30.0)
    scn = _downlink(Terminal(aperture_m=0.7, wavelength_m=lam, detector=SMF()))

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        term = smf_fast_term(scn, geom, n_samples=300)

    assert term.name == "receive coupling (SMF)"
    assert term.meta["model"] == "fast-modal"
    assert term.category == "coupling"
    # The mean loss is at least the static mode-match floor.
    assert term.mean_db >= term.meta["floor_db"] - 1e-6
    # Stochastic with an empirical quantile deeper than the mean.
    rng = np.random.default_rng(0)
    assert term.stochastic and term.quantile is not None
    q99 = term.quantile_db(0.99)
    assert q99 is not None and q99 > term.mean_db
    draws = term.sample_db(5000, rng)
    assert draws.shape == (5000,) and np.all(np.isfinite(draws))
    # Subharmonics on (captures the tilt); point-ahead caveat flagged.
    assert term.meta["subharmonics"] is True
    assert "subharmonics" in term.assumptions.validity

    # Spectrum label follows the RESOLVED scales, not a fixed constant. The
    # default (L0=inf, l0=1e-6) is the Kolmogorov limit; a finite L0 flips the
    # label to von Karman and names the scale.
    from ..assumptions import SPECTRUM_KOLMOGOROV, SPECTRUM_VON_KARMAN
    assert term.assumptions.spectrum == SPECTRUM_KOLMOGOROV
    assert np.isinf(term.meta["L0_m"]) and term.meta["l0_m"] == 1e-6
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        term_vk = smf_fast_term(scn, geom, n_samples=200,
                                fast_params={"L0": 25.0})
    assert term_vk.assumptions.spectrum == SPECTRUM_VON_KARMAN
    assert term_vk.meta["L0_m"] == 25.0
    assert "von Karman" in term_vk.assumptions.validity
    # No precompensation on a downlink, so the point-ahead caveat stays silent.
    assert not any("Point-ahead" in v for v in term.assumptions.violations)
    # No correction -> NOAO, no ZMAX.
    assert term.meta["ao_mode"] == "NOAO" and term.meta["zmax"] is None

    # WP3d: the traced olb-side physics (the plane-wave amplitude index) is in the
    # provenance, and the external FAST Monte Carlo self-declares.
    from ..results import Budget
    prov = term.assumptions.provenance
    assert any("plane_wave_scintillation_index" in s for s in prov), prov
    assert "untraced: fast-aosim" in prov, prov
    # The demo (coupling category) Term passes the untraced guard.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        checked = Budget([term]).check(warn=False)
    assert not any("did not open the assumption collection context" in reason
                   for _name, reason in checked), checked

    # The migrated REGIME hard-flag: a 10 deg downlink drives the plane-wave
    # amplitude index past the weak boundary (sigma2_I ~ 1.6), so the TRACED
    # weak-regime check on plane_wave_scintillation_index fires. The Term is
    # not ok, and the violation carries the physics source prefix.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        term_low = smf_fast_term(scn, CircularOrbit(1500e3, 10.0), n_samples=100)
    assert not term_low.assumptions.ok
    assert any("plane_wave_scintillation_index" in v and "weak" in v
               for v in term_low.assumptions.violations), \
        term_low.assumptions.violations

    # AO order is honoured: AO(n_modes) -> FAST MODAL ZMAX=n_modes, and more
    # corrected modes give less coupling loss. Tip-tilt sits between NOAO and AO.
    def build(comp, n=400):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return smf_fast_term(_downlink(Terminal(
                aperture_m=0.7, wavelength_m=lam, detector=SMF(),
                compensation=comp)), geom, n_samples=n)

    t_tt = build([TipTilt()])
    t_ao20 = build([TipTilt(), AO(20)])
    t_ao6 = build([TipTilt(), AO(6)])
    assert t_tt.meta["ao_mode"] == "TT" and t_tt.meta["zmax"] is None
    assert t_ao20.meta["ao_mode"] == "AO" and t_ao20.meta["zmax"] == 20
    assert t_ao6.meta["zmax"] == 6                        # n_modes maps to ZMAX
    # More corrected modes -> less loss. NOAO > TT > AO(6) > AO(20).
    assert term.mean_db > t_tt.mean_db > t_ao6.mean_db > t_ao20.mean_db, (
        term.mean_db, t_tt.mean_db, t_ao6.mean_db, t_ao20.mean_db)

    # An elevation array is refused in this first cut.
    try:
        smf_fast_term(scn, CircularOrbit(1500e3, np.array([30.0, 60.0])))
        raise AssertionError("array elevation must raise")
    except ValueError:
        pass

    print(f"FAST fidelity-1 SMF (0.7 m, 30 deg): "
          f"NOAO {term.mean_db:.2f} dB | TT {t_tt.mean_db:.2f} dB | "
          f"AO(6) {t_ao6.mean_db:.2f} dB | AO(20) {t_ao20.mean_db:.2f} dB")

    # --- pre-compensated uplink (uplink_fast_term) ---------------------------
    from ..scenario import DownlinkBeacon

    def _uplink(comp, aperture_m=1.5):
        return SpaceScenario(
            ground=Terminal(aperture_m=aperture_m, wavelength_m=lam,
                            transmitter=Transmitter(waist_m=0.2),
                            compensation=comp),
            space=Terminal(aperture_m=0.05, wavelength_m=lam),
            direction="uplink",
            channel=Channel(site=Site(cn2_ground=1.7e-14), altitude_m=600e3),
            precompensation=DownlinkBeacon())

    up_geom = CircularOrbit(600e3, elevation_deg=60.0)
    up_scn = _uplink([TipTilt(), AO(60)])

    def up_term(scn, n=400, **kw):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return uplink_fast_term(scn, up_geom, n_samples=n, **kw)

    u60 = up_term(up_scn)
    assert u60.name == "turbulence (pre-compensated, FAST)"
    assert u60.category == "turbulence" and u60.meta["model"] == "fast-precomp-uplink"
    # Stochastic: a 99% quantile is deeper than the mean.
    assert u60.stochastic and u60.quantile is not None and not u60.mean_only
    u_q99 = u60.quantile_db(0.99)
    assert u_q99 is not None and u_q99 > u60.mean_db, (u_q99, u60.mean_db)
    assert np.all(np.isfinite(u60.sample_db(2000, np.random.default_rng(1))))
    # The orbit gives a real point-ahead offset.
    assert u60.meta["dtheta_arcsec"] > 0, u60.meta["dtheta_arcsec"]
    assert u60.meta["zmax"] == 60 and u60.meta["ao_mode"] == "AO"

    # WP3d: the pre-compensated uplink Term carries the same traced provenance and
    # the external FAST self-declaration, and passes the untraced guard.
    up_prov = u60.assumptions.provenance
    assert any("plane_wave_scintillation_index" in s for s in up_prov), up_prov
    assert "untraced: fast-aosim" in up_prov, up_prov
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        up_checked = Budget([u60]).check(warn=False)
    assert not any("did not open the assumption collection context" in reason
                   for _name, reason in up_checked), up_checked

    # More corrected modes -> less mean loss.
    u6 = up_term(_uplink([TipTilt(), AO(6)]))
    assert u6.mean_db > u60.mean_db, (u6.mean_db, u60.mean_db)

    # Point-ahead costs loss: a zero offset gives less mean loss than the orbit.
    u_nopa = up_term(up_scn, fast_params={"DTHETA": [0, 0]})
    assert u_nopa.meta["dtheta_arcsec"] == 0.0
    assert u_nopa.mean_db < u60.mean_db, (u_nopa.mean_db, u60.mean_db)

    # A downlink scenario is refused.
    try:
        uplink_fast_term(scn, up_geom, n_samples=100)
        raise AssertionError("a downlink scenario must raise")
    except ValueError:
        pass

    # A ground terminal with no AO stage is refused.
    try:
        uplink_fast_term(_uplink([TipTilt()]), up_geom, n_samples=100)
        raise AssertionError("no AO stage must raise")
    except ValueError:
        pass

    print(f"FAST pre-compensated uplink (1.5 m, 60 deg, "
          f"{u60.meta['dtheta_arcsec']:.2f} arcsec point-ahead): "
          f"AO(6) {u6.mean_db:.2f} dB | AO(60) {u60.mean_db:.2f} dB | "
          f"AO(60) no point-ahead {u_nopa.mean_db:.2f} dB "
          f"(99% fade {u_q99:.2f} dB)")
    print("self-check passed")
