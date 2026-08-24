'''
Fidelity-1 single-mode-fibre coupling from FAST (the true LP01 modal overlap).

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

from ...results import Term
from ...assumptions import (Assumptions, BEAM_GAUSSIAN, REGIME_WEAK,
                            SPECTRUM_KOLMOGOROV, SPECTRUM_VON_KARMAN)
from ...terminal import TipTilt, AO
from ...turbulence.profiles import DEFAULT_HS, default_cn2_profile
from ...turbulence.plane_wave_scintillation import (plane_wave_scintillation_index,
                                                   WEAK_FLUCTUATION_LIMIT)


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
    sigma2_I_amp = float(plane_wave_scintillation_index(
        elev, rx.wavelength_m, hs, cn2_profile))

    def quantile(p):
        return float(np.quantile(loss_db, p))

    def sampler(n, rng):
        return rng.choice(loss_db, size=n, replace=True)

    # FAST is a von Karman engine (it always calls turb_powerspectrum_vonKarman).
    # olb sets the scales itself (L0=inf, l0=1e-6 above), so the spectrum is
    # Kolmogorov by our own choice, not inherited from FAST's conf.py. A finite
    # L0 (or a large l0) from fast_params makes it a true von Karman spectrum, so
    # read the label from the RESOLVED scales, not from a fixed constant.
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
    assumptions = Assumptions(
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
    # First-cut limitation: no point-ahead.
    assumptions.flag(
        "Point-ahead is off (DTHETA=0): the up-leg / down-leg anisoplanatism of a "
        "moving satellite is not modelled. Pass fast_params={'DTHETA': [x, y]} in "
        "arcsec to include it."
    )
    # AMPLITUDE saturation: only here does FAST's log-normal scintillation break.
    # A large coupled-power fade (deep 99% tail) does NOT trip this -- that fade is
    # phase-driven and modelled correctly by the screens.
    if sigma2_I_amp > WEAK_FLUCTUATION_LIMIT:
        assumptions.flag(
            f"Plane-wave amplitude scintillation sigma2_I={sigma2_I_amp:.2f} exceeds "
            f"the weak-fluctuation limit {WEAK_FLUCTUATION_LIMIT}; FAST's log-normal "
            "scintillation (the amplitude part) departs from data in saturation. The "
            "phase-driven coupling fade is still modelled by the screens, but the "
            "amplitude contribution to the fade tail is not trustworthy. Raise the "
            "elevation, or note the amplitude regime."
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
            "amplitude_regime_weak": bool(sigma2_I_amp <= WEAK_FLUCTUATION_LIMIT),
            "r0_los_m": float(getattr(sim, "r0_los", np.nan)),
            "L0_m": L0,
            "l0_m": l0,
            "n_samples": int(params["NITER"]),
        },
        assumptions=assumptions,
    )


if __name__ == '__main__':
    import warnings

    from ...scenario import SpaceScenario, Channel, Site
    from ...geometry import CircularOrbit
    from ...terminal import Terminal, Transmitter, SMF, TipTilt, AO

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
    assert term.quantile_db(0.99) > term.mean_db
    draws = term.sample_db(5000, rng)
    assert draws.shape == (5000,) and np.all(np.isfinite(draws))
    # Subharmonics on (captures the tilt); point-ahead caveat flagged.
    assert term.meta["subharmonics"] is True
    assert "subharmonics" in term.assumptions.validity

    # Spectrum label follows the RESOLVED scales, not a fixed constant. The
    # default (L0=inf, l0=1e-6) is the Kolmogorov limit; a finite L0 flips the
    # label to von Karman and names the scale.
    from ...assumptions import SPECTRUM_KOLMOGOROV, SPECTRUM_VON_KARMAN
    assert term.assumptions.spectrum == SPECTRUM_KOLMOGOROV
    assert np.isinf(term.meta["L0_m"]) and term.meta["l0_m"] == 1e-6
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        term_vk = smf_fast_term(scn, geom, n_samples=200,
                                fast_params={"L0": 25.0})
    assert term_vk.assumptions.spectrum == SPECTRUM_VON_KARMAN
    assert term_vk.meta["L0_m"] == 25.0
    assert "von Karman" in term_vk.assumptions.validity
    assert any("Point-ahead" in v for v in term.assumptions.violations)
    # No correction -> NOAO, no ZMAX.
    assert term.meta["ao_mode"] == "NOAO" and term.meta["zmax"] is None

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
    print("self-check passed")
