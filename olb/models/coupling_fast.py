'''
Fidelity-1 single-mode-fibre coupling from FAST (the true LP01 modal overlap).

The fidelity-0 (Dikmelik/Marechal) and the reciprocity Strehl-proxy models do NOT
compute the true fibre-mode overlap. This module does. It drives FAST (the
`fast-aosim` package) to get the downlink single-mode-fibre coupling as the
coherent overlap of the turbulent aperture field with the back-projected fibre
mode:

    eta(t) = |integral (Aperture * M_fibre) * exp(chi + i*phi) dA| ^ 2
             / |integral (Aperture * M_fibre) dA| ^ 2

FAST propagates Monte-Carlo phase screens (phase) with an aperture-averaged
log-normal scintillation (log-amplitude), and forms the mode overlap directly. So
this is a true modal-coupling metric, not a Strehl proxy. See olb.models.coupling
for the lower-fidelity models and the fidelity ladder.

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
from ..assumptions import (Assumptions, BEAM_GAUSSIAN, REGIME_WEAK,
                           SPECTRUM_KOLMOGOROV)
from ..terminal import TipTilt, AO
from ..turbulence.profiles import DEFAULT_HS, default_cn2_profile


def _load_fast():
    '''Import fast-aosim lazily, with a helpful error when it is absent.'''
    try:
        import fast
        return fast
    except ImportError as e:
        raise ImportError(
            "the fidelity-1 SMF coupling needs the `fast-aosim` package. "
            "Run `pip install fast-aosim`, or use a lower fidelity "
            "(rx_coupling_term smf_fidelity='reciprocity')."
        ) from e


def _cn2_layers(cn2_zenith, hs):
    '''
    Integrated Cn2 per layer [m^1/3] from a zenith Cn2(h) profile.

    FAST wants the integrated Cn2 per layer at ZENITH (it applies the airmass
    itself). So multiply the profile by the layer thickness.
    '''
    return np.asarray(cn2_zenith, float) * np.gradient(np.asarray(hs, float))


def _ao_mode(compensation):
    '''Map the olb compensation stack to a FAST AO_MODE (first cut: NOAO/TT/AO).'''
    if any(isinstance(c, AO) for c in compensation):
        return "AO"
    if any(isinstance(c, TipTilt) for c in compensation):
        return "TT"
    return "NOAO"


def smf_fast_term(scenario, geometry, *, hs=None, cn2_profile=None,
                  n_samples=1000, fast_params=None):
    '''
    Fidelity-1 SMF receive-coupling Term for a downlink, computed by FAST.

    Parameters:
        scenario : Scenario
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
    ao_mode = _ao_mode(rx.compensation)

    params = dict(
        WVL=rx.wavelength_m,
        D_GROUND=D,
        OBSC_GROUND=obsc_ratio * D,          # FAST wants the obscuration DIAMETER
        SMF=True,
        PROP_DIR="down",                     # SMF coupling is a receive-side quantity
        AO_MODE=ao_mode,
        W0="opt",                            # FAST optimises the fibre-mode size
        H_SAT=scenario.channel.altitude_m,
        ZENITH_ANGLE=90.0 - elev,
        H_TURB=hs,
        CN2_TURB=cn2_layer,
        WIND_SPD=wind * np.ones_like(hs),
        WIND_DIR=np.zeros_like(hs),
        DTHETA=[0, 0],                       # first cut: no point-ahead
        NITER=int(n_samples),
        LOGLEVEL="ERROR",
    )
    if fast_params:
        params.update(fast_params)

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

    def quantile(p):
        return float(np.quantile(loss_db, p))

    def sampler(n, rng):
        return rng.choice(loss_db, size=n, replace=True)

    assumptions = Assumptions(
        beam_type=BEAM_GAUSSIAN,
        turbulence_regime=REGIME_WEAK,
        spectrum=SPECTRUM_KOLMOGOROV,
        validity="Fidelity-1 single-mode-fibre coupling: the true LP01 Gaussian-"
                 "mode overlap under turbulence, from FAST (fast-aosim). The "
                 "received field is the aperture times the fibre mode, propagated "
                 "through Monte-Carlo phase screens with an aperture-averaged log-"
                 "normal scintillation. eta_max and the fibre-mode size come from "
                 "FAST (W0='opt'). Weak-to-moderate fluctuation (log-normal "
                 "scintillation).",
    )
    # First-cut limitation: no point-ahead.
    assumptions.flag(
        "Point-ahead is off (DTHETA=0): the up-leg / down-leg anisoplanatism of a "
        "moving satellite is not modelled. Pass fast_params={'DTHETA': [x, y]} in "
        "arcsec to include it."
    )
    # NOAO tilt sampling depends on the FAST grid.
    if ao_mode == "NOAO":
        assumptions.flag(
            "NOAO low-order (tilt) accuracy depends on the FAST grid (NPXLS); the "
            "auto grid may undersample tilt. Pass fast_params={'NPXLS': ...} for "
            "production accuracy."
        )

    note = (f"FAST fidelity-1 modal coupling, AO_MODE={ao_mode}, "
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
            "ao_mode": ao_mode,
            "floor_db": floor_db,
            "scintillation_index": float(result.scintillation_index),
            "r0_los_m": float(getattr(sim, "r0_los", np.nan)),
            "n_samples": int(params["NITER"]),
        },
        assumptions=assumptions,
    )


if __name__ == '__main__':
    import warnings

    from ..scenario import Scenario, Channel, Site
    from ..geometry import CircularOrbit
    from ..terminal import Terminal, Transmitter, SMF

    try:
        _load_fast()
    except ImportError as e:
        print(f"fast-aosim not installed ({e.__class__.__name__}); skipping the "
              "fidelity-1 self-check.")
        raise SystemExit

    lam = 1550e-9

    def _downlink(ground):
        return Scenario(
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
    # Point-ahead and NOAO-grid caveats are flagged.
    assert any("Point-ahead" in v for v in term.assumptions.violations)
    assert any("NOAO" in v for v in term.assumptions.violations)

    # An elevation array is refused in this first cut.
    try:
        smf_fast_term(scn, CircularOrbit(1500e3, np.array([30.0, 60.0])))
        raise AssertionError("array elevation must raise")
    except ValueError:
        pass

    print(f"FAST fidelity-1 SMF (0.7 m, no AO, 30 deg): mean {term.mean_db:.2f} dB, "
          f"99% {term.quantile_db(0.99):.2f} dB, floor {term.meta['floor_db']:.2f} dB, "
          f"S.I. {term.meta['scintillation_index']:.2f}")
    print("self-check passed")
