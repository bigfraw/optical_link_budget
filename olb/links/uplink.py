'''
Uplink Terms and budget assembly (beam wander + scintillation).

This module gives the coupled-flux turbulence Term for a ground-launched uplink
beam, and the uplink budget that assembles it with the geometric, atmospheric,
and pointing Terms.

The turbulence Term wraps the coupled-flux Monte Carlo. It estimates the
turbulence-induced fade at the satellite receiver. It captures beam wander and
scintillation together. There is no closed form for the coupled fade, so this is
a MONTE-CARLO-ONLY Term. It gives a sampler and sets ``quantile=None``. This
value tells the budget to evaluate the term with ``monte_carlo()``, not with the
analytic fade sum.
'''

import numpy as np

from ..results import Budget, Term
from ..assumptions import (Assumptions, BEAM_GAUSSIAN, REGIME_WEAK,
                          SPECTRUM_KOLMOGOROV)
from ..models.geometric import geometric_loss_term
from ..models.transmittance import atmospheric_loss_term, DEFAULT_TAU_ZENITH
from ..models.pointing import pointing_loss_term
from ..models.gaussian_efficiency import tx_gaussian_efficiency_term
from ..turbulence.coupled_flux import _flux_result, WEAK_FLUCTUATION_LIMIT
from ..turbulence.profiles import DEFAULT_HS, default_cn2_profile

# Below this launch-truncation loss the beam is an untruncated Gaussian, so the
# transmit Gaussian-efficiency term is skipped [dB].
TX_TRUNCATION_MIN_DB = 1e-2


def uplink_turbulence_term(scenario, geometry, n_samples=3000, n_apertures=1,
                           hs=None, cn2_profile=None):
    '''
    Monte-Carlo turbulence Term (uplink beam wander + scintillation).

    MC-only: it gives a sampler and sets quantile=None, so the budget evaluates
    it with monte_carlo(). The code fills ``mean_db`` from a representative draw
    at construction, so the budget table still has a value.

    Parameters:
        scenario : Scenario
            Reads the transmit terminal (waist w0, divergence, wavelength) and
            site.cn2_ground (passed as the HV57 ground scale hv57_A).
        geometry : CircularOrbit or TLEPass
            Reads elevation_deg and slant_range_m. Scalar elevation -> the
            sampler returns shape (n,); an elevation array -> shape (n, E),
            evaluated with one MC per elevation (expensive).
        n_samples : int
            MC draws used for the construction-time mean estimate.
        n_apertures : int
            Independent on-axis samples averaged per receiver aperture
            (receive-side aperture averaging).
        hs : numpy.ndarray, optional
            Turbulence altitude grid [m]. Defaults to DEFAULT_HS.
        cn2_profile : numpy.ndarray, optional
            Explicit Cn2(h) profile at zenith matching ``hs``. If None the
            kernel builds an HV57 profile (requires the `fast` package).

    Returns:
        Term
            name="turbulence (coupled-flux)", category="turbulence".
    '''
    hs = DEFAULT_HS if hs is None else hs
    tx = scenario.tx_terminal
    w0 = tx.transmitter.waist_m
    wavelength = tx.wavelength_m
    hv57_A = scenario.channel.site.cn2_ground
    divergence_rad = tx.transmitter.divergence_rad

    elev = np.atleast_1d(np.asarray(geometry.elevation_deg, dtype=float))
    ranges = np.atleast_1d(np.asarray(geometry.slant_range_m, dtype=float))
    scalar = np.ndim(geometry.elevation_deg) == 0

    # Representative draw per elevation -> table mean + validity metadata.
    reps = [_flux_result(w0, e, r, wavelength, hs, cn2_profile, hv57_A,
                         n_samples, n_apertures, divergence_rad=divergence_rad)
            for e, r in zip(elev, ranges)]
    mean_db = np.array([-10 * np.log10(np.mean(rep["Is_summed"])) for rep in reps])
    sigma2_x = np.array([rep["sigma2_x_mean"] for rep in reps])
    valid = np.array([rep["weak_fluctuation_valid"] for rep in reps])

    assumptions = Assumptions(
        beam_type=BEAM_GAUSSIAN,
        turbulence_regime=REGIME_WEAK,
        spectrum=SPECTRUM_KOLMOGOROV,
        validity="Rytov weak fluctuation: sigma2_x < 0.6 (WEAK_FLUCTUATION_LIMIT). "
                 "Divergence enters the beam broadening AND the scintillation "
                 "index (through the diverged receiver-plane Lambda and Theta).",
    )
    if not np.all(valid):
        worst = float(sigma2_x[~valid].max())
        assumptions.flag(
            f"sigma2_x={worst:.2f} exceeds the weak-fluctuation limit "
            f"{WEAK_FLUCTUATION_LIMIT}; scintillation approaches saturation."
        )

    def sampler(n, rng):
        # rng bridge: coupled_flux_montecarlo draws from numpy's GLOBAL RNG
        # (np.random), not a passed Generator, so seed the global RNG from the
        # budget's seeded `rng` to keep the draw reproducible.
        np.random.seed(int(rng.integers(0, 2 ** 32 - 1)))
        cols = [-10 * np.log10(
                    _flux_result(w0, e, r, wavelength, hs, cn2_profile, hv57_A,
                                 n, n_apertures, divergence_rad=divergence_rad)["Is_summed"])
                for e, r in zip(elev, ranges)]   # one MC per elevation (expensive)
        return cols[0] if scalar else np.stack(cols, axis=1)

    return Term(
        name="turbulence (coupled-flux)",
        category="turbulence",
        mean_db=float(mean_db[0]) if scalar else mean_db,
        sampler=sampler,
        quantile=None,   # MC-only: no closed form -> budget must monte_carlo()
        note="uplink beam wander + scintillation, coupled-flux Monte Carlo",
        meta={
            "weak_fluctuation_valid": bool(valid[0]) if scalar else valid,
            "sigma2_x": float(sigma2_x[0]) if scalar else sigma2_x,
            "weak_fluctuation_limit": WEAK_FLUCTUATION_LIMIT,
            "w_diffraction_limited": reps[0]["w_diffraction_limited"] if scalar
                else np.array([rep["w_diffraction_limited"] for rep in reps]),
            "w_st": reps[0]["w_st"] if scalar
                else np.array([rep["w_st"] for rep in reps]),
            "n_apertures": n_apertures,
        },
        assumptions=assumptions,
    )


def uplink_budget(scenario, geometry, *, turbulence=True, tau_zenith=None,
                  n_samples=3000, cn2_profile=None):
    '''
    Assemble the uplink budget: geometric, atmospheric, pointing, turbulence.

    Add the coupled-flux turbulence Term when `turbulence` is true. Build a
    default Cn2 profile when `cn2_profile` is None, so the budget runs without
    the `fast` package.

    Parameters:
        scenario : Scenario
            The link case.
        geometry : CircularOrbit or TLEPass
            The link geometry.
        turbulence : bool
            Add the coupled-flux turbulence Term when true.
        tau_zenith : float, optional
            Zenith optical depth. Defaults to transmittance.DEFAULT_TAU_ZENITH.
        n_samples : int
            Monte Carlo draws for the turbulence Term mean estimate.
        cn2_profile : numpy.ndarray, optional
            Explicit zenith Cn2 profile. Defaults to default_cn2_profile.

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
    # The transmit Gaussian-efficiency term is opt-in. It fires only when the
    # transmit terminal has a Transmitter and its launch aperture truncates the
    # beam by more than TX_TRUNCATION_MIN_DB. A wide aperture leaves the beam an
    # untruncated Gaussian, so the term is skipped.
    tx = scenario.tx_terminal
    if tx.transmitter is not None:
        eff = tx_gaussian_efficiency_term(scenario, geometry)
        if eff.mean_db > TX_TRUNCATION_MIN_DB:
            terms.append(eff)
    if turbulence:
        if cn2_profile is None:
            cn2_profile = default_cn2_profile(scenario.channel.site)
        terms.append(uplink_turbulence_term(scenario, geometry, n_samples=n_samples,
                                            cn2_profile=cn2_profile))
    return Budget(terms, scenario=scenario)


if __name__ == '__main__':
    from ..scenario import Scenario, Channel
    from ..geometry import CircularOrbit
    from ..terminal import Terminal, Transmitter, Aperture

    def _uplink(w0, *, divergence=None, power=None, jitter=0.0,
                ground_aperture=0.5, ground_obscuration=0.0,
                space_aperture=0.05, sensitivity=None):
        '''Build an uplink Scenario: tx=ground, rx=space (satellite).'''
        detector = None if sensitivity is None else Aperture(sensitivity_dbm=sensitivity)
        return Scenario(
            ground=Terminal(aperture_m=ground_aperture, obscuration_ratio=ground_obscuration,
                            wavelength_m=1550e-9, pointing_jitter_rad=jitter,
                            transmitter=Transmitter(waist_m=w0, power_dbm=power,
                                                    divergence_rad=divergence)),
            space=Terminal(aperture_m=space_aperture, wavelength_m=1550e-9,
                           detector=detector),
            direction="uplink", channel=Channel(altitude_m=600e3))

    scenario = _uplink(0.1)
    rng = np.random.default_rng(0)

    # Is the `fast` package available? Try a build without an explicit profile.
    try:
        uplink_turbulence_term(scenario, CircularOrbit(600e3, 55.0), n_samples=500)
        fast_available = True
        cn2 = None
    except ImportError as e:
        fast_available = False
        cn2 = 1e-16 * np.ones_like(DEFAULT_HS)   # moderate substitute profile so the check still runs
        print(f"`fast` unavailable ({e.__class__.__name__}); using explicit cn2_profile.")

    for elevation_deg in [30,60,90]:
        geom = CircularOrbit(600e3, float(elevation_deg))
        term = uplink_turbulence_term(scenario, geom, n_samples=2000, cn2_profile=cn2)
        samples = term.sample_db(3000, rng)
        fade_99 = np.percentile(samples, 99)

        print('=' * 40)
        print(f"Elevation: {elevation_deg} deg")
        print(f"Mean turbulence loss: {term.mean_db:.2f} dB")
        print(f"99% fade:             {fade_99:.2f} dB")
        print(f"sigma2_x={term.meta['sigma2_x']:.3f} "
              f"(weak_fluctuation_valid={term.meta['weak_fluctuation_valid']})")

        assert samples.shape == (3000,)
        assert np.all(np.isfinite(samples))
        assert term.quantile_db(0.99) is None   # MC-only: no closed form
        assert fade_99 > term.mean_db           # a 99% fade is deeper than the mean loss

    print('\n' + '=' * 40)
    # weak_fluctuation_valid must follow the threshold: negligible Cn2 -> valid,
    # strong Cn2 -> invalid. This does not depend on whether the sweep above is
    # inside the trusted regime.
    weak_cn2 = 1e-18 * np.ones_like(DEFAULT_HS)    # negligible turbulence
    strong_cn2 = 1e-15 * np.ones_like(DEFAULT_HS)   # strong, sigma2_x finite but past the limit
    geom = CircularOrbit(600e3, 90.0)
    valid_term = uplink_turbulence_term(scenario, geom, n_samples=1000, cn2_profile=weak_cn2)
    invalid_term = uplink_turbulence_term(scenario, geom, n_samples=1000, cn2_profile=strong_cn2)
    assert valid_term.meta["weak_fluctuation_valid"] is True
    assert invalid_term.meta["weak_fluctuation_valid"] is False
    assert valid_term.assumptions is not None
    assert valid_term.assumptions.ok            # weak Cn2 -> no violation
    assert not invalid_term.assumptions.ok      # strong Cn2 -> violation flagged

    print(f"weak Cn2  -> weak_fluctuation_valid={valid_term.meta['weak_fluctuation_valid']}")
    print(f"strong Cn2 -> weak_fluctuation_valid={invalid_term.meta['weak_fluctuation_valid']}")

    # Divergence: it now enters the beam broadening AND the scintillation index.
    # A diverged beam is wider and more spherical-wave-like, so it both dilutes
    # the broadening loss and scintillates less. Neither link raises a
    # divergence-specific violation, because the model no longer approximates it.
    from .._deps import w0_to_div
    tx0 = scenario.tx_terminal
    theta_min = w0_to_div(tx0.transmitter.waist_m, tx0.wavelength_m)
    div_scn = _uplink(0.1, divergence=5 * theta_min)
    moderate_cn2 = 1e-16 * np.ones_like(DEFAULT_HS)
    np.random.seed(0)
    coll_term = uplink_turbulence_term(scenario, geom, n_samples=4000, cn2_profile=moderate_cn2)
    np.random.seed(0)
    div_term = uplink_turbulence_term(div_scn, geom, n_samples=4000, cn2_profile=moderate_cn2)
    assert not any("Divergence" in v for v in div_term.assumptions.violations)
    # Diverging widens the free-space baseline and dilutes the turbulence loss.
    assert div_term.meta["w_diffraction_limited"] > coll_term.meta["w_diffraction_limited"]
    assert div_term.mean_db < coll_term.mean_db, (div_term.mean_db, coll_term.mean_db)
    # The diverged beam scintillates less (lower log-amplitude variance).
    assert div_term.meta["sigma2_x"] < coll_term.meta["sigma2_x"], (
        div_term.meta["sigma2_x"], coll_term.meta["sigma2_x"])
    print(f"collimated sigma2_x={coll_term.meta['sigma2_x']:.4f}, "
          f"diverged sigma2_x={div_term.meta['sigma2_x']:.4f}")

    # --- uplink budget self-check -------------------------------------------
    # A wide launch aperture (1.5 m for a 0.2 m waist) leaves the beam untruncated,
    # so the transmit Gaussian-efficiency term does not fire.
    budget_scn = _uplink(0.2, power=40, jitter=2e-6, sensitivity=-40,
                         ground_aperture=1.5)
    budget_geom = CircularOrbit(altitude_m=600e3, elevation_deg=60.0)
    up = uplink_budget(budget_scn, budget_geom,
                       cn2_profile=default_cn2_profile(budget_scn.channel.site))
    assert up.to_frame().shape[0] == 4, up.to_frame().shape
    up_mc = up.monte_carlo(2000, rng=np.random.default_rng(0), availabilities=(0.99,))
    up_margin = up_mc["margin_db"][0.99]
    assert np.isfinite(up_margin), up_margin

    # A narrow launch aperture (0.15 m for a 0.2 m waist) truncates the beam, so
    # the transmit Gaussian-efficiency term fires.
    ap_scn = _uplink(0.2, power=40, jitter=2e-6, sensitivity=-40,
                     ground_aperture=0.15, ground_obscuration=0.3)
    up_ap = uplink_budget(ap_scn, budget_geom,
                          cn2_profile=default_cn2_profile(ap_scn.channel.site))
    assert up_ap.to_frame().shape[0] == 5, up_ap.to_frame().shape
    eff = next(t for t in up_ap.terms if t.category == "system")
    assert eff.mean_db > 0                       # truncation is a loss
    assert up_ap.total_loss_db() > up.total_loss_db()   # aperture truncation costs margin

    print('\n' + '=' * 40)
    print(up.to_frame().to_string(index=False))
    print(f"\nuplink 60 deg 99% margin: {up_margin:.2f} dB")
    print(f"with aperture: +{eff.mean_db:.2f} dB transmit truncation ({eff.note})")
    print(f"fast_available={fast_available}")
    print("self-check passed.")
