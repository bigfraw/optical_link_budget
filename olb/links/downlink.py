'''
Downlink scintillation Term and budget assembly.

This module builds the downlink scintillation Term and the downlink budget. The
Term gives the analytic lognormal plane-wave fade for a LEO-to-ground link. The
pure scintillation physics lives in olb.turbulence.scintillation.

The received power P is lognormal with E[P] = 1. With sigma_l^2 = ln(1 +
sigma2_P) the loss faces are
    mean_db     = (5/ln10) * sigma_l^2
    quantile(p) = -10*log10( exp(-sigma_l^2/2 + sigma_l * Phi_inv(1-p)) )
Source of the lognormal irradiance PDF: Andrews and Phillips, Laser Beam
Propagation through Random Media, 2nd ed. (2005), Ch. 5. Phi_inv is the inverse
standard normal CDF.

Validity: the Rytov (lognormal) model is a weak-fluctuation model. The code
carries sigma2_I in Term.meta and sets a weak_fluctuation_valid flag. It gives a
warning when sigma2_I exceeds WEAK_FLUCTUATION_LIMIT. Above that limit use the
gamma-gamma or the Monte Carlo model.
'''

import warnings

import numpy as np
from scipy.stats import norm

from ..results import Budget, Term
from ..assumptions import (Assumptions, BEAM_PLANE_WAVE, REGIME_WEAK,
                          SPECTRUM_KOLMOGOROV)
from ..models.geometric import geometric_loss_term
from ..models.transmittance import atmospheric_loss_term, DEFAULT_TAU_ZENITH
from ..models.pointing import pointing_loss_term
from ..turbulence.profiles import DEFAULT_HS, default_cn2_profile
from ..turbulence.scintillation import (plane_wave_scintillation_index,
                                        aperture_averaged_scintillation_index,
                                        WEAK_FLUCTUATION_LIMIT)

_LN10 = np.log(10.0)


def _lognormal_term(scenario, geometry, *, aperture_average, hs, cn2_profile):
    '''
    Build the analytic lognormal downlink scintillation Term.

    Give all three faces: mean_db, quantile(p), sampler(n, rng). Carry sigma2_I
    and the weak_fluctuation_valid flag in Term.meta. Warn when sigma2_I exceeds
    the weak-fluctuation limit.
    '''
    link = scenario.link
    wavelength = link.wavelength_m
    D = link.rx_diameter_m
    elev = geometry.elevation_deg

    sigma2_I = plane_wave_scintillation_index(elev, wavelength, hs, cn2_profile)
    if aperture_average:
        # Use the distributed-path aperture-averaging integral directly.
        sigma2_P = aperture_averaged_scintillation_index(D, elev, wavelength,
                                                         hs, cn2_profile)
        A = np.asarray(sigma2_P) / np.asarray(sigma2_I)
    else:
        A = 1.0
        sigma2_P = sigma2_I

    sigma_l2 = np.log(1.0 + sigma2_P)
    sigma_l = np.sqrt(sigma_l2)
    base_shape = np.shape(sigma_l2)

    mean_db = (5.0 / _LN10) * sigma_l2

    def quantile(p):
        # Loss exceeded a fraction (1 - p) of the time. Phi_inv = norm.ppf.
        z = norm.ppf(1.0 - p)
        return -10.0 / _LN10 * (-sigma_l2 / 2.0 + sigma_l * z)

    def sampler(n, rng):
        P = rng.lognormal(mean=-sigma_l2 / 2.0, sigma=sigma_l,
                          size=(n, *base_shape))
        return -10.0 * np.log10(P)

    # The Rytov approximation validity is a turbulence property. Test it with
    # the point sigma2_I, not the aperture-averaged sigma2_P.
    valid = np.asarray(sigma2_I) < WEAK_FLUCTUATION_LIMIT
    assumptions = Assumptions(
        beam_type=BEAM_PLANE_WAVE,
        turbulence_regime=REGIME_WEAK,
        spectrum=SPECTRUM_KOLMOGOROV,
        validity="Weak fluctuation: sigma2_I < 0.25. Aperture averaging uses "
                 "the distributed-path spectral integral over height and "
                 "spatial wavenumber. The aperture-averaging filter assumes a "
                 "uniform circular aperture with no central obscuration.",
    )
    # The circular-aperture filter models an unobscured aperture. A central
    # obscuration (Cassegrain secondary) breaks that. Flag the violation.
    obscuration = getattr(scenario.link, "rx_obscuration_ratio", 0.0)
    if aperture_average and obscuration > 0.0:
        assumptions.flag(
            f"The receive aperture has a central obscuration "
            f"(ratio={obscuration:.3f}); the circular-aperture averaging "
            "filter does not model it."
        )
    if not np.all(valid):
        worst = float(np.asarray(sigma2_I)[~valid].max())
        assumptions.flag(
            f"sigma2_I={worst:.3f} exceeds the weak-fluctuation limit 0.25; "
            "use gamma-gamma or Monte Carlo."
        )
        warnings.warn(
            f"plane-wave scintillation index sigma2_I={np.max(sigma2_I):.3f} >= "
            f"{WEAK_FLUCTUATION_LIMIT} -- the Rytov weak-fluctuation model is "
            "exceeded. The lognormal fade is not trusted. Use the gamma-gamma "
            "model or the Monte Carlo model."
        )

    return Term(
        name="scintillation",
        category="turbulence",
        mean_db=float(mean_db) if base_shape == () else mean_db,
        sampler=sampler,
        quantile=quantile,
        note="plane-wave downlink lognormal scintillation, aperture-averaged",
        meta={
            "model": "lognormal",
            "sigma2_I": float(sigma2_I) if base_shape == () else np.asarray(sigma2_I),
            "sigma2_P": float(sigma2_P) if base_shape == () else np.asarray(sigma2_P),
            "aperture_averaging_factor": float(A) if np.ndim(A) == 0 else np.asarray(A),
            "weak_fluctuation_valid": bool(valid) if base_shape == () else valid,
            "weak_fluctuation_limit": WEAK_FLUCTUATION_LIMIT,
        },
        assumptions=assumptions,
    )


def _gamma_gamma_term(scenario, geometry, *, aperture_average, hs, cn2_profile):
    '''
    Reserved slot: the analytic gamma-gamma downlink scintillation Term.

    The gamma-gamma model covers the moderate-to-strong fluctuation regime. Use
    it when sigma2_I exceeds the weak-fluctuation limit. This slot is not
    implemented yet.
    '''
    raise NotImplementedError(
        "the gamma_gamma model is a reserved slot. It will hold the analytic "
        "gamma-gamma downlink scintillation Term for the moderate-to-strong "
        "fluctuation regime. Use model='lognormal' for now."
    )


def _montecarlo_term(scenario, geometry, *, aperture_average, hs, cn2_profile):
    '''
    Reserved slot: the Monte Carlo downlink scintillation Term.

    The Monte Carlo model uses phase-screen field propagation. It lives in the
    heavier Monte Carlo area of the code. This slot is not implemented yet.
    '''
    raise NotImplementedError(
        "the montecarlo model is a reserved slot. It will hold the phase-screen "
        "field-propagation downlink scintillation Term, in the heavier Monte "
        "Carlo area. Use model='lognormal' for now."
    )


def _auto_select(scenario, geometry, *, aperture_average, hs, cn2_profile):
    '''
    Selector layer: choose the downlink scintillation model from the regime.

    This layer chooses the best model from the fluctuation regime. For now it
    returns the lognormal Term. It warns when sigma2_I exceeds the
    weak-fluctuation limit. In that case it should use the gamma-gamma model or
    the Monte Carlo model once those slots exist.
    '''
    term = _lognormal_term(scenario, geometry, aperture_average=aperture_average,
                           hs=hs, cn2_profile=cn2_profile)
    if not np.all(term.meta["weak_fluctuation_valid"]):
        warnings.warn(
            "auto selector: sigma2_I exceeds the weak-fluctuation limit. The "
            "lognormal Term is returned now. The selector should use the "
            "gamma-gamma model or the Monte Carlo model once those slots exist."
        )
    return term


_MODELS = {
    "lognormal": _lognormal_term,
    "gamma_gamma": _gamma_gamma_term,
    "montecarlo": _montecarlo_term,
    "auto": _auto_select,
}


def downlink_scintillation_term(scenario, geometry, *, model="lognormal",
                                aperture_average=True, hs=None, cn2_profile=None):
    '''
    Build the downlink scintillation Term.

    Dispatch to the requested model. The "lognormal" model is the analytic
    weak-fluctuation model. The "auto" model is the selector layer.

    Parameters:
        scenario : Scenario
            Reads link.wavelength_m and link.rx_diameter_m. Reads the site to
            build the default Cn2 profile with get_c2n.
        geometry : CircularOrbit or TLEPass
            Reads elevation_deg. A scalar elevation gives a scalar Term. An
            elevation array gives a Term that broadcasts over that shape.
        model : str
            One of "lognormal", "gamma_gamma", "montecarlo", "auto".
        aperture_average : bool
            Apply the plane-wave aperture-averaging factor when true.
        hs : numpy.ndarray, optional
            Heights above the ground station [m]. Defaults to DEFAULT_HS.
        cn2_profile : numpy.ndarray, optional
            Zenith Cn2(h) profile on the hs grid. Defaults to the site profile
            from get_c2n.

    Returns:
        Term
            name="scintillation", category="turbulence".

    Raises:
        ValueError
            If model is not a known name.
        NotImplementedError
            If model is a reserved slot ("gamma_gamma" or "montecarlo").
    '''
    if model not in _MODELS:
        raise ValueError(
            f"unknown model {model!r}. Use one of {sorted(_MODELS)}."
        )
    hs = DEFAULT_HS if hs is None else hs
    if cn2_profile is None:
        cn2_profile = default_cn2_profile(scenario.site, hs)
    return _MODELS[model](scenario, geometry, aperture_average=aperture_average,
                          hs=hs, cn2_profile=cn2_profile)


def downlink_budget(scenario, geometry, *, tau_zenith=None, scintillation=True):
    '''
    Assemble the downlink budget: geometric, atmospheric, pointing, scintillation.

    The coupled-flux Term models a ground-launched uplink beam. It does not
    apply to the downlink. The downlink scintillation Term now exists. It gives
    the analytic lognormal plane-wave fade. Every downlink Term has a closed-form
    quantile, so the downlink budget supports analytic fade. The budget also
    supports Monte Carlo.

    Parameters:
        scenario : Scenario
            The link case.
        geometry : CircularOrbit or TLEPass
            The link geometry.
        tau_zenith : float, optional
            Zenith optical depth. Defaults to transmittance.DEFAULT_TAU_ZENITH.
        scintillation : bool
            Add the lognormal downlink scintillation Term when true.

    Returns:
        Budget
            The budget with the scenario set.
    '''
    tau = DEFAULT_TAU_ZENITH if tau_zenith is None else tau_zenith
    terms = [
        geometric_loss_term(scenario, geometry),
        atmospheric_loss_term(scenario, geometry, tau_zenith=tau),
        pointing_loss_term(scenario, geometry),
    ]
    if scintillation:
        terms.append(downlink_scintillation_term(scenario, geometry,
                                                 model="lognormal",
                                                 aperture_average=True))
    return Budget(terms, scenario=scenario)


if __name__ == '__main__':
    from ..scenario import Scenario, Link
    from ..geometry import CircularOrbit

    scenario = Scenario(
        link=Link(wavelength_m=1550e-9, direction="downlink", tx_waist_m=0.035,
                  rx_diameter_m=0.7),
        altitude_m=600e3,
    )
    hs = DEFAULT_HS
    cn2 = default_cn2_profile(scenario.site, hs)

    # The lognormal Term gives a working quantile deeper than the mean loss.
    geom = CircularOrbit(600e3, 30.0)
    term = downlink_scintillation_term(scenario, geom, cn2_profile=cn2)
    assert term.name == "scintillation", term.name
    q99 = term.quantile_db(0.99)
    assert np.isfinite(q99) and q99 > term.mean_db, (q99, term.mean_db)
    assert term.assumptions is not None

    # A 20 deg elevation breaks the weak-fluctuation limit; 60 deg does not.
    import warnings as _w
    with _w.catch_warnings():
        _w.simplefilter("ignore")
        t20 = downlink_scintillation_term(scenario, CircularOrbit(600e3, 20.0),
                                          cn2_profile=cn2)
        t60 = downlink_scintillation_term(scenario, CircularOrbit(600e3, 60.0),
                                          cn2_profile=cn2)
    assert not t20.assumptions.ok
    assert t60.assumptions.ok

    # The sampled mean matches mean_db for large n.
    rng = np.random.default_rng(0)
    sampled = term.sample_db(200_000, rng)
    assert sampled.shape == (200_000,)
    assert abs(sampled.mean() - term.mean_db) < 0.02, (sampled.mean(), term.mean_db)

    # The reserved slots raise NotImplementedError.
    for slot in ("gamma_gamma", "montecarlo"):
        try:
            downlink_scintillation_term(scenario, geom, model=slot, cn2_profile=cn2)
        except NotImplementedError:
            pass
        else:
            raise AssertionError(f"model={slot!r} must raise NotImplementedError")

    # An elevation sweep broadcasts over the geometry shape.
    sweep = CircularOrbit(600e3, np.array([20.0, 45.0, 90.0]))
    sweep_term = downlink_scintillation_term(scenario, sweep, cn2_profile=cn2)
    assert np.shape(sweep_term.mean_db) == (3,)
    assert sweep_term.sample_db(1000, rng).shape == (1000, 3)

    # A central obscuration flags a filter-validity violation at 60 deg (weak
    # fluctuation holds). A zero obscuration does not flag it.
    scen_obs = Scenario(
        link=Link(wavelength_m=1550e-9, direction="downlink", tx_waist_m=0.035,
                  rx_diameter_m=0.7, rx_obscuration_ratio=0.3),
        altitude_m=600e3,
    )
    cn2_obs = default_cn2_profile(scen_obs.site, hs)
    t_obs = downlink_scintillation_term(scen_obs, CircularOrbit(600e3, 60.0),
                                        cn2_profile=cn2_obs)
    assert not t_obs.assumptions.ok, "obscuration=0.3 must flag a violation"
    assert t60.assumptions.ok, "obscuration=0.0 must not flag a violation"

    # --- downlink budget self-check -----------------------------------------
    budget_scn = Scenario(
        link=Link(wavelength_m=1550e-9, direction="downlink", tx_waist_m=0.035,
                  rx_diameter_m=0.7, tx_power_dbm=40, rx_sensitivity_dbm=-40,
                  pointing_jitter_rad=2e-6),
        altitude_m=600e3,
    )
    down = downlink_budget(budget_scn, CircularOrbit(altitude_m=600e3,
                                                    elevation_deg=60.0))
    assert down.to_frame().shape[0] == 4, down.to_frame().shape
    # Every downlink Term has a closed-form quantile, so analytic fade works.
    down_fade = down.fade_margin_db(0.99)
    assert np.isfinite(down_fade), down_fade
    # The plain downlink budget shows the base Term name "scintillation".
    assert "scintillation" in [t.name for t in down.terms], [t.name for t in down.terms]

    print(f"mean loss = {term.mean_db:.4f} dB   99% fade = {q99:.4f} dB")
    print(f"sampled mean = {sampled.mean():.4f} dB (n=200000)")
    print(down.to_frame().to_string(index=False))
    print("self-check passed")
