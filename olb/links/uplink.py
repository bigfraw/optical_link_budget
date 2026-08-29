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
from ..assumptions import (Assumptions, BEAM_GAUSSIAN, BEAM_PLANE_WAVE,
                          REGIME_WEAK, SPECTRUM_KOLMOGOROV)
from ..models.geometric import geometric_loss_term
from ..models.extinction import slant_extinction_term, DEFAULT_TAU_ZENITH
from ..models.pointing import pointing_loss_term
from ..models.gaussian_efficiency import tx_gaussian_efficiency_term
from ..turbulence.anisoplanatism import (anisoplanatic_phase_variance,
                                         max_radial_order)
from ..turbulence.ao import plane_wave_fried_parameter_profile, apply_compensation
from ..turbulence.uplink_flux import _flux_result, WEAK_FLUCTUATION_LIMIT
from ..turbulence.plane_wave_scintillation import plane_wave_scintillation_index
from ..turbulence.profiles import DEFAULT_HS, default_cn2_profile
from ..terminal import AO
from ..scenario import DownlinkBeacon, LaserGuideStar

# Below this launch-truncation loss the beam is an untruncated Gaussian, so the
# transmit Gaussian-efficiency term is skipped [dB].
TX_TRUNCATION_MIN_DB = 1e-2

# Above this residual phase variance [rad^2] the extended Marechal mean
# eta = exp(-sigma2) departs from the true on-axis mean. The real far field
# breaks into a speckled core plus a halo, and the exponential decays faster
# than the real core, so the Term overstates the loss. Source: T. S. Ross,
# Appl. Opt. 48(10), 1812 (2009), DOI 10.1364/AO.48.001812.
MARECHAL_SIGMA2_MAX = 1.0


def _flag_marechal(assumptions, sigma2):
    '''Flag a residual phase variance past the extended-Marechal limit.'''
    worst = float(np.max(sigma2))
    if worst > MARECHAL_SIGMA2_MAX:
        assumptions.flag(
            f"MARECHAL LIMIT: the residual phase variance sigma2={worst:.2f} "
            f"rad^2 is more than {MARECHAL_SIGMA2_MAX:g} rad^2. The extended "
            "Marechal mean eta = exp(-sigma2) is a small-residual form. Past "
            "this limit the exponential decays faster than the true on-axis "
            "mean, so the Term overstates the loss. Source: T. S. Ross, "
            "Appl. Opt. 48(10), 1812 (2009), DOI 10.1364/AO.48.001812."
        )


def uplink_turbulence_term(scenario, geometry, n_samples=3000, n_apertures=1,
                           hs=None, cn2_profile=None):
    '''
    Monte-Carlo turbulence Term (uplink beam wander + scintillation).

    MC-only: it gives a sampler and sets quantile=None, so the budget evaluates
    it with monte_carlo(). The code fills ``mean_db`` from a representative draw
    at construction, so the budget table still has a value.

    Parameters:
        scenario : SpaceScenario
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
    # Mechanical pointing jitter folds into the beam-wander displacement inside
    # the coupled-flux MC, so this Term now carries BOTH the turbulence wander
    # and the tracking jitter. uplink_budget therefore drops the standalone
    # pointing-loss Term when turbulence is on (adding both double-counts it).
    sigma_theta = tx.pointing_jitter_rad

    elev = np.atleast_1d(np.asarray(geometry.elevation_deg, dtype=float))
    ranges = np.atleast_1d(np.asarray(geometry.slant_range_m, dtype=float))
    scalar = np.ndim(geometry.elevation_deg) == 0

    # Representative draw per elevation -> table mean + validity metadata.
    reps = [_flux_result(w0, e, r, wavelength, hs, cn2_profile, hv57_A,
                         n_samples, n_apertures, divergence_rad=divergence_rad,
                         sigma_theta_rad=sigma_theta)
            for e, r in zip(elev, ranges)]
    mean_db = np.array([-10 * np.log10(np.mean(rep["Is_summed"])) for rep in reps])
    sigma2_x = np.array([rep["sigma2_x_mean"] for rep in reps])
    valid = np.array([rep["weak_fluctuation_valid"] for rep in reps])

    assumptions = Assumptions(
        beam_type=BEAM_GAUSSIAN,
        turbulence_regime=REGIME_WEAK,
        spectrum=SPECTRUM_KOLMOGOROV,
        validity="Rytov weak fluctuation: sigma2_x < 0.25 (WEAK_FLUCTUATION_LIMIT). "
                 "Divergence enters the beam broadening AND the scintillation "
                 "index (through the diverged receiver-plane Lambda and Theta). "
                 "Mechanical pointing jitter folds into the beam-wander "
                 "displacement, so this Term carries the tracking-jitter loss and "
                 "fade too (no separate uplink pointing Term). "
                 "The Dios coupled-flux analysis assumes an untruncated Gaussian "
                 "launch beam, so it does not model a central obscuration on the "
                 "launch aperture. The obscuration MEAN loss is carried elsewhere "
                 "(tx_gaussian_efficiency_term). The size of the obscuration effect "
                 "on this fade is UNRESOLVED (see the investigation note).",
    )
    if not np.all(valid):
        worst = float(sigma2_x[~valid].max())
        assumptions.flag(
            f"sigma2_x={worst:.2f} exceeds the weak-fluctuation limit "
            f"{WEAK_FLUCTUATION_LIMIT}; scintillation approaches saturation."
        )
    # Dios assumes an untruncated Gaussian launch beam. A central obscuration on
    # the launch aperture (the Transmitter override, else the Terminal value)
    # breaks that. Flag the violation.
    tx_obsc = (tx.transmitter.obscuration_ratio
               if tx.transmitter.obscuration_ratio is not None
               else tx.obscuration_ratio)
    if tx_obsc > 0.0:
        assumptions.flag(
            f"The launch aperture has a central obscuration (ratio={tx_obsc:.3f}); "
            "the Dios coupled-flux analysis assumes an untruncated Gaussian beam "
            "and does not model it."
        )

    def sampler(n, rng):
        # rng bridge: the coupled-flux kernels (coupled_flux_sample) draw from
        # numpy's GLOBAL RNG (np.random), not a passed Generator, so seed the
        # global RNG from the budget's seeded `rng` to keep the draw reproducible.
        np.random.seed(int(rng.integers(0, 2 ** 32 - 1)))
        cols = [-10 * np.log10(
                    _flux_result(w0, e, r, wavelength, hs, cn2_profile, hv57_A,
                                 n, n_apertures, divergence_rad=divergence_rad,
                                 sigma_theta_rad=sigma_theta)["Is_summed"])
                for e, r in zip(elev, ranges)]   # one MC per elevation (expensive)
        return cols[0] if scalar else np.stack(cols, axis=1)

    return Term(
        name="turbulence (coupled-flux)",
        category="turbulence",
        mean_db=float(mean_db[0]) if scalar else mean_db,
        sampler=sampler,
        quantile=None,   # MC-only: no closed form -> budget must monte_carlo()
        note="uplink beam wander + jitter + scintillation, coupled-flux Monte Carlo",
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


def uplink_point_ahead_term(scenario, geometry, hs=None, cn2_profile=None,
                            max_order='auto'):
    '''
    Point-ahead anisoplanatism Term (uplink pre-compensation residual).

    This Term is the error of a downlink-beacon uplink pre-compensation. The
    terminal senses the turbulence on the downlink beam and applies the conjugate
    to the uplink beam. The up and down paths share the same turbulence
    (reciprocity), so the downlink phase gives the uplink correction. But the
    downlink arrives from where the satellite was, and the uplink goes to where
    the satellite will be. The two directions differ by the point-ahead angle.

    The correction removes the part of each Zernike order that stays correlated
    across that angle. The DECORRELATION residual stays. The error is the sum of
    that residual over the corrected orders 2..max_order, with the piston and the
    two tilts left out (a separate tracking loop points the beam). The residual
    per order is 2 sigma_n^2 (1 - rho_n): it is small for a well-correlated low
    order and it saturates at twice the mode variance for a fully decorrelated
    order. So the loss grows with the adaptive-optics order, up to an
    infinite-order limit. This is NOT a penalty for correcting. It is the part of
    the turbulence that the two directions do not share. See
    anisoplanatic_phase_variance and Fig. 2 of Stone et al. (1994).

    The Term is MEAN-ONLY: it gives the expected loss and no fade. It has no
    sampler and no quantile, because the phase variance is a steady-state
    ensemble value with no time-domain draw.

    Parameters:
        scenario : SpaceScenario
            Reads the transmit terminal (aperture, wavelength, compensation) and
            site (through the default Cn2 profile).
        geometry : CircularOrbit or TLEPass
            Reads elevation_deg and point_ahead_rad. Scalar elevation -> a scalar
            mean_db; an elevation array -> one value per elevation.
        hs : numpy.ndarray, optional
            Turbulence altitude grid [m]. Defaults to DEFAULT_HS.
        cn2_profile : numpy.ndarray, optional
            Explicit Cn2(h) profile at zenith matching ``hs``. If None the code
            builds it with default_cn2_profile from the site.
        max_order : "auto", int, or None
            Highest Zernike radial order that the correction touches. "auto" (the
            default) reads it from the transmit terminal. An AO(n_modes) stage
            gives max_radial_order(n_modes). No AO stage gives None (the ideal
            infinite-order limit). An int forces that order. None forces the
            infinite-order limit.

    Returns:
        Term
            name="point-ahead anisoplanatism", category="anisoplanatism",
            mean_only=True.
    '''
    hs = DEFAULT_HS if hs is None else hs
    # Resolve the profile here, because a caller can use this Term alone.
    if cn2_profile is None:
        cn2_profile = default_cn2_profile(scenario.channel.site, hs)
    tx = scenario.tx_terminal
    D = tx.aperture_m
    wavelength = tx.wavelength_m

    # Map the compensation stack to a corrected radial order. The largest AO
    # stage sets it. With no AO stage, fall back to the infinite-order limit.
    if max_order == 'auto':
        ao_modes = [c.n_modes for c in tx.compensation if isinstance(c, AO)]
        max_order = max_radial_order(max(ao_modes)) if ao_modes else None

    elev = np.atleast_1d(np.asarray(geometry.elevation_deg, dtype=float))
    theta = np.broadcast_to(
        np.asarray(geometry.point_ahead_rad, dtype=float), elev.shape)
    scalar = np.ndim(geometry.elevation_deg) == 0

    # anisoplanatic_phase_variance takes ONE angle and ONE elevation at a time,
    # because hs is already a grid. So loop over the elevations.
    sigma2 = np.array([
        anisoplanatic_phase_variance(D, t, hs, cn2_profile, wavelength,
                                     remove='piston_tilt', max_order=max_order,
                                     elevation_deg=e)
        for e, t in zip(elev, theta)])
    # Extended Marechal: eta = exp(-sigma2), so the loss is -10*log10(eta).
    # Source: V. W. S. Chan and others; extended Marechal approximation.
    # Derivation and validity: T. S. Ross, Appl. Opt. 48(10), 1812 (2009),
    # DOI 10.1364/AO.48.001812. The same relation is in olb.models.coupling.
    loss_db = (10.0 / np.log(10.0)) * sigma2

    order_note = "all orders" if max_order is None else f"orders 2..{max_order}"
    theta_urad = theta * 1e6
    note = ("point-ahead anisoplanatism, " + order_note + ", theta="
            + (f"{theta_urad[0]:.2f}" if scalar
               else f"{theta_urad.min():.2f}-{theta_urad.max():.2f}") + " urad")

    assumptions = Assumptions(
        beam_type=BEAM_GAUSSIAN,
        turbulence_regime=REGIME_WEAK,
        spectrum=SPECTRUM_KOLMOGOROV,
        validity="Decorrelation residual of a downlink-beacon uplink "
                 "pre-compensation. The terminal senses the turbulence on the "
                 "downlink beam and applies the conjugate to the uplink beam. The "
                 "correction removes the part of each Zernike order that stays "
                 "correlated across the point-ahead angle. The error is the "
                 "decorrelation residual summed over the corrected orders "
                 + order_note + ", with the piston and the two tilts left out (a "
                 "separate tracking loop points the beam). The residual per order "
                 "is 2 sigma_n^2 (1 - rho_n). It grows with the corrected order. "
                 "Source: Stone et al. (1994), DOI 10.1364/JOSAA.11.000347. "
                 "The phase variance becomes a loss with the extended Marechal "
                 "approximation, the same relation as in olb.models.coupling "
                 "(V. W. S. Chan and others; extended Marechal approximation). "
                 "This Term is the anisoplanatic part only. The companion "
                 "uplink_fitting_term gives the uncorrected high-order error. "
                 "NEITHER Term models the scintillation, so it must not be added "
                 "to a full uncorrected turbulence Term; the two stand in for the "
                 "corrected turbulence error. "
                 "The point-ahead angle comes from geometry.point_ahead_rad, thus "
                 "from my_analysis_modules.satellite.SatellitePass."
                 "point_ahead_angle(). That function uses the simple form "
                 "2 * v_orbit * sin(elevation) / c. It does not use the more "
                 "general form 2 * omega_line_of_sight * slant_range / c, and its "
                 "source file gives no citation. This is a limit of the input "
                 "accuracy. This Term does not correct it. "
                 "This Term gives the mean loss only. It models no fade.",
    )
    # BIG LIMITATION: the pre-compensated uplink model is phase-only and
    # mean-only. Adaptive optics corrects the phase; it does not remove the
    # amplitude scintillation. No trustworthy analytic model exists for the
    # scintillation of a pre-compensated beam (decision 2026-08-27): the
    # correction decorrelates over the point-ahead angle mode by mode, and a
    # decorrelated correction reshapes the beam, so the analytic normalisation
    # breaks. The model of record is the fidelity-1 FAST Monte Carlo.
    assumptions.flag(
        "NO SCINTILLATION, NO FADE: the pre-compensated uplink budget models "
        "the phase (wavefront) only, and this Term gives the mean Strehl loss "
        "only. Adaptive optics does not remove the amplitude scintillation, "
        "and no trustworthy analytic model exists for the scintillation of a "
        "pre-compensated beam. The model of record is the fidelity-1 FAST "
        "Monte Carlo with the point-ahead offset (olb.models.fast, "
        "uplink_fast_term). Call uplink_budget with fidelity=1 "
        "to get it. Do not read a fade for a pre-compensated uplink from this "
        "budget."
    )
    _flag_marechal(assumptions, sigma2)
    # With no adaptive-optics stage the correction order is unknown, so the Term
    # uses the infinite-order limit. That is an UPPER bound of the true error.
    if max_order is None:
        assumptions.flag(
            "No adaptive-optics stage sets the corrected order, so this Term "
            "uses the infinite-order limit. This is an upper bound. Pass an "
            "AO(n_modes) stage, or set max_order, for the true adaptive-optics "
            "order."
        )

    return Term(
        name="point-ahead anisoplanatism",
        category="anisoplanatism",
        mean_db=float(loss_db[0]) if scalar else loss_db,
        note=note,
        meta={
            "theta_paa_rad": float(theta[0]) if scalar else np.asarray(theta),
            "sigma2_rad2": float(sigma2[0]) if scalar else sigma2,
            "max_order": max_order,
        },
        assumptions=assumptions,
        mean_only=True,   # fidelity-0: expected residual only, no fade (see results.Budget)
    )


def uplink_fitting_term(scenario, geometry, hs=None, cn2_profile=None):
    '''
    Adaptive-optics fitting-error Term for the uplink (uncorrected high orders).

    The adaptive optics corrects the low Zernike orders. The high orders stay
    uncorrected. This Term gives the loss of that uncorrected wavefront error. It
    is the companion of uplink_point_ahead_term: the point-ahead Term gives the
    decorrelation residual of the CORRECTED orders, and this Term gives the full
    error of the UNCORRECTED orders. The two mode sets do not overlap, so the two
    Terms add.

    The residual is the Noll variance after the correction:
        sigma^2 = c * (D / r0)^(5/3)
    with the Noll coefficient c set by the compensation stack (see
    olb.turbulence.ao). An empty stack gives c = NOLL_PISTON, so the Term is then
    the total uncorrected phase variance (the piston removed). r0 is the
    plane-wave Fried parameter at the ground aperture. By reciprocity the up and
    down paths share the ground-aperture phase, so the plane-wave (downlink) r0
    sets the sensed and the corrected wavefront.

    The variance becomes a loss with the extended Marechal approximation, the same
    relation as in olb.models.coupling. The Term is MEAN-ONLY: it has no sampler
    and no quantile.

    Parameters:
        scenario : SpaceScenario
            Reads the transmit terminal (aperture, wavelength, compensation) and
            site (through the default Cn2 profile).
        geometry : CircularOrbit or TLEPass
            Reads elevation_deg. Scalar elevation -> a scalar mean_db; an
            elevation array -> one value per elevation.
        hs : numpy.ndarray, optional
            Turbulence altitude grid [m]. Defaults to DEFAULT_HS.
        cn2_profile : numpy.ndarray, optional
            Explicit Cn2(h) profile at zenith matching ``hs``. If None the code
            builds it with default_cn2_profile from the site.

    Returns:
        Term
            name="AO fitting error", category="fitting", mean_only=True.
    '''
    hs = DEFAULT_HS if hs is None else hs
    if cn2_profile is None:
        cn2_profile = default_cn2_profile(scenario.channel.site, hs)
    tx = scenario.tx_terminal
    D = tx.aperture_m
    wavelength = tx.wavelength_m
    elev = geometry.elevation_deg
    scalar = np.ndim(elev) == 0

    # Reciprocity: the ground-aperture phase is common to the up and down paths,
    # so the plane-wave (downlink) r0 sets the sensed and corrected wavefront.
    r0 = plane_wave_fried_parameter_profile(cn2_profile, hs, wavelength, elev)
    residual = apply_compensation(tx.compensation, D, r0)
    sigma2_fit = np.asarray(residual.variance, dtype=float)
    # Extended Marechal: eta = exp(-sigma2), so the loss is -10*log10(eta).
    # Source: V. W. S. Chan and others; extended Marechal approximation.
    # Derivation and validity: T. S. Ross, Appl. Opt. 48(10), 1812 (2009),
    # DOI 10.1364/AO.48.001812. The same relation is in olb.models.coupling.
    loss_db = (10.0 / np.log(10.0)) * sigma2_fit

    note = (f"AO fitting error, {residual.n_modes} modes corrected, "
            f"Noll c={residual.coefficient:.4f}")
    assumptions = Assumptions(
        beam_type=BEAM_PLANE_WAVE,
        turbulence_regime=REGIME_WEAK,
        spectrum=SPECTRUM_KOLMOGOROV,
        validity="Uncorrected high-order wavefront error of the uplink adaptive "
                 "optics. It is the Noll residual sigma^2 = c (D/r0)^(5/3) after "
                 "the correction, with c set by the compensation stack. An empty "
                 "stack gives the total uncorrected phase variance (piston "
                 "removed). Source: R. J. Noll, JOSA 66(3), 207 (1976), "
                 "DOI 10.1364/JOSA.66.000207. r0 is the plane-wave Fried "
                 "parameter at the ground aperture; by reciprocity it is common "
                 "to the up and down paths. The variance becomes a loss with the "
                 "extended Marechal approximation (olb.models.coupling; "
                 "V. W. S. Chan and others). It holds for a small residual "
                 "(sigma2 <= MARECHAL_SIGMA2_MAX; a larger residual gets a "
                 "flag). This Term gives the mean loss only. It models no "
                 "fade and no scintillation.",
    )
    # BIG LIMITATION: the pre-compensated uplink model is phase-only and
    # mean-only. See uplink_point_ahead_term for the full decision note.
    assumptions.flag(
        "NO SCINTILLATION, NO FADE: the pre-compensated uplink budget models "
        "the phase (wavefront) only, and this Term gives the mean Strehl loss "
        "only. Adaptive optics does not remove the amplitude scintillation, "
        "and no trustworthy analytic model exists for the scintillation of a "
        "pre-compensated beam. The model of record is the fidelity-1 FAST "
        "Monte Carlo with the point-ahead offset (olb.models.fast, "
        "uplink_fast_term). Call uplink_budget with fidelity=1 "
        "to get it. Do not read a fade for a pre-compensated uplink from this "
        "budget."
    )
    _flag_marechal(assumptions, sigma2_fit)
    return Term(
        name="AO fitting error",
        category="fitting",
        mean_db=float(loss_db) if scalar else loss_db,
        note=note,
        meta={
            "sigma2_fit_rad2": float(sigma2_fit) if scalar else sigma2_fit,
            "r0_m": float(r0) if scalar else np.asarray(r0),
            "n_modes": int(residual.n_modes),
            "noll_c": float(residual.coefficient),
        },
        assumptions=assumptions,
        mean_only=True,   # fidelity-0: expected residual only, no fade
    )


def _uplink_fidelity2_terms(scenario, geometry, wave, hs, cn2_profile):
    '''
    The two fidelity-2 wave-optics Terms of an UNCORRECTED uplink.

    A space uplink cannot be simulated end to end (the turbulent runner
    propagates only the ~20 km slab as a downlink and reads the uplink through
    reciprocity; the full slant range is absent), so the loss splits into:
      - a DETERMINISTIC vacuum-optics Term (the full no-turbulence loss over the
        slant range: launch truncation + geometric spread + satellite-aperture
        capture, from the co-moving-grid vacuum run);
      - a STOCHASTIC turbulence Term from the reciprocity overlap eta_turb
        (Shapiro, DOI 10.1364/JOSA.61.000492), the pure turbulence penalty.
    Together they replace the analytic geometric, launch-truncation, and
    coupled-flux Terms. The tracking jitter stays in the standalone pointing Term
    (the reciprocity overlap holds no jitter). `wave` is a Fidelity2Bundle from
    olb.models.waveoptics.run_fidelity2.

    The reciprocity route reads the SAME screens up and down, so it does NOT model
    an adaptive-optics correction or the point-ahead decorrelation. So it fits the
    uncorrected uplink only, not a pre-compensated one.
    '''
    from ..models.waveoptics import (waveoptics_vacuum_term,
                                     waveoptics_turbulence_term)
    elev = np.asarray(geometry.elevation_deg, dtype=float)
    if elev.ndim != 0:
        raise ValueError(
            "the fidelity-2 uplink takes a scalar elevation (one range per "
            "record). Loop over the elevations and build one bundle each."
        )
    sigma2_I = float(plane_wave_scintillation_index(
        float(elev), scenario.tx_terminal.wavelength_m, hs, cn2_profile))
    vac = waveoptics_vacuum_term(wave.vacuum, include_smf=False,
                                 beam_type=BEAM_GAUSSIAN)
    pen = waveoptics_turbulence_term(
        wave.turbulent, quantity="eta_turb", beam_type=BEAM_GAUSSIAN,
        sigma2_I=sigma2_I,
        note="uplink turbulence penalty by reciprocity (wave optics). PURE "
             "turbulence penalty: the free-space spread, the launch truncation, "
             "and the satellite-aperture capture are in the vacuum-optics Term, "
             "and the tracking jitter is in the pointing Term.")
    return [vac, pen]


def uplink_budget(scenario, geometry, *, fidelity=1, turbulence=True,
                  tau_zenith=None, n_samples=3000, cn2_profile=None, wave=None):
    '''
    Assemble the uplink budget at a chosen fidelity.

    `fidelity` is a WHOLE-PATH choice (see the README fidelity ladder). The routes
    also depend on the pre-compensation source (SpaceScenario.precompensation):

      UNCORRECTED (no source, or a tip-tilt-only beacon):
        fidelity=0 : raises. There is no analytic mean-only model for an
            uncorrected uplink (beam wander plus scintillation has no closed
            form).
        fidelity=1 (the default): the coupled-flux Monte Carlo Term
            (uplink_turbulence_term, beam wander + scintillation). It carries the
            tracking jitter too, so there is no standalone pointing Term.
        fidelity=2 : two wave-optics Terms -- a deterministic vacuum-optics Term
            (launch truncation + geometric spread + satellite-aperture capture)
            and a stochastic turbulence Term from the reciprocity overlap
            eta_turb. They replace the geometric, launch-truncation, and
            coupled-flux Terms. The reciprocity overlap holds no jitter, so the
            standalone pointing Term stays.
      PRE-COMPENSATED (DownlinkBeacon with an AO stage):
        fidelity=0 : the AO error budget -- two analytic mean-only phase Terms:
            the AO fitting error (uplink_fitting_term) and the point-ahead
            anisoplanatism (uplink_point_ahead_term). PHASE-ONLY and MEAN-ONLY:
            no scintillation, no fade (decision 2026-08-27, the pre-compensated
            beam has no trustworthy analytic scintillation form). The Terms flag
            it, so Budget.check() warns.
        fidelity=1 (the default): ONE FAST Monte-Carlo Term
            (uplink_fast_term) with the residual phase, the point-ahead
            decorrelation, and the uncorrected log-amplitude, by reciprocity. It
            carries a real fade. It needs the optional `fast-aosim` package. The
            standalone pointing Term carries the tracking jitter.
        fidelity=2 : raises. The reciprocity screens carry no adaptive-optics
            correction or point-ahead decorrelation, so wave optics does not model
            a pre-compensated uplink. Use fidelity=1 (FAST).
      LaserGuideStar: not implemented yet. It raises NotImplementedError.

    Pointing jitter: the coupled-flux Term (uncorrected fidelity 1) carries the
    tracking jitter itself. Every other route leaves it to a standalone pointing
    Term, so the jitter is never lost and never double-counted.

    Parameters:
        scenario : SpaceScenario
            The link case. Reads `precompensation` for the uplink correction.
        geometry : CircularOrbit or TLEPass
            The link geometry.
        fidelity : int
            0 (analytic), 1 (statistical, the default), or 2 (wave optics, needs
            `wave` and an uncorrected uplink).
        turbulence : bool
            Add the turbulence Term when true (fidelity 0/1). When false the
            budget is geometric-only.
        tau_zenith : float, optional
            Zenith optical depth. Defaults to extinction.DEFAULT_TAU_ZENITH.
        n_samples : int
            Monte Carlo draws for the coupled-flux Term (fidelity 1 uncorrected)
            and the FAST draw count (fidelity 1 pre-compensated).
        cn2_profile : numpy.ndarray, optional
            Explicit zenith Cn2 profile. Defaults to default_cn2_profile.
        wave : Fidelity2Bundle, optional
            The precomputed wave-optics records for fidelity=2. Run it with
            olb.models.waveoptics.run_fidelity2.

    Returns:
        Budget
            The budget with the scenario set.

    Raises:
        NotImplementedError
            If the scenario uses a LaserGuideStar pre-compensation source.
        ValueError
            If fidelity is not 0/1/2; if fidelity=0 for an uncorrected uplink; if
            fidelity=2 for a pre-compensated uplink or without a `wave` bundle.
        ImportError
            If fidelity=1 pre-compensated and `fast-aosim` is not installed.
    '''
    if fidelity not in (0, 1, 2):
        raise ValueError(f"fidelity must be 0, 1, or 2, got {fidelity!r}.")
    tau = DEFAULT_TAU_ZENITH if tau_zenith is None else tau_zenith
    tx = scenario.tx_terminal

    # Resolve the pre-compensation source. A laser guide star is not modelled
    # yet. A downlink beacon with an AO stage pre-compensates the uplink.
    pc = scenario.precompensation
    if isinstance(pc, LaserGuideStar):
        raise NotImplementedError(
            "laser-guide-star pre-compensation is not modelled yet. Its focal "
            "(cone) anisoplanatism differs from the downlink-beacon point-ahead "
            "anisoplanatism. Use a DownlinkBeacon or no pre-compensation for now."
        )
    precomp = (isinstance(pc, DownlinkBeacon)
               and any(isinstance(c, AO) for c in tx.compensation))

    if fidelity == 2:
        if precomp:
            raise ValueError(
                "fidelity=2 does not model a PRE-COMPENSATED uplink. The "
                "reciprocity screens carry no adaptive-optics correction or "
                "point-ahead decorrelation. Use fidelity=1 (FAST) for the "
                "pre-compensated case, or remove the pre-compensation source."
            )
        if wave is None:
            raise ValueError(
                "fidelity=2 needs a precomputed `wave` bundle. Run "
                "olb.models.waveoptics.run_fidelity2(scenario, geometry, ...) and "
                "pass it as wave. The budget does not run the split-step "
                "propagation implicitly."
            )
        if cn2_profile is None:
            cn2_profile = default_cn2_profile(scenario.channel.site)
        # The vacuum-optics Term owns the geometric spread and launch truncation;
        # the reciprocity Term holds no jitter, so the pointing Term stays.
        terms = [
            slant_extinction_term(scenario, geometry, tau_zenith=tau),
            pointing_loss_term(scenario, geometry),
        ]
        terms += _uplink_fidelity2_terms(scenario, geometry, wave, DEFAULT_HS,
                                         cn2_profile)
        return Budget(terms, scenario=scenario)

    # fidelity 0/1: the analytic backbone.
    terms = [
        geometric_loss_term(scenario, geometry),
        slant_extinction_term(scenario, geometry, tau_zenith=tau),
    ]
    # The coupled-flux Term (uncorrected fidelity 1) carries the jitter itself; on
    # every other route the standalone pointing Term carries it.
    if not turbulence or precomp:
        terms.append(pointing_loss_term(scenario, geometry))
    # The transmit Gaussian-efficiency (launch truncation) Term is opt-in.
    if tx.transmitter is not None:
        eff = tx_gaussian_efficiency_term(scenario, geometry)
        if eff.mean_db > TX_TRUNCATION_MIN_DB:
            terms.append(eff)
    if turbulence:
        if cn2_profile is None:
            cn2_profile = default_cn2_profile(scenario.channel.site)
        if precomp and fidelity == 1:
            # Fidelity 1, pre-compensated: ONE FAST Monte-Carlo Term (residual
            # phase + point-ahead + log-amplitude, by reciprocity). Import here to
            # keep the `fast-aosim` dependency optional.
            from ..models.fast import uplink_fast_term
            terms.append(uplink_fast_term(scenario, geometry,
                                          cn2_profile=cn2_profile,
                                          n_samples=n_samples))
        elif precomp:
            # Fidelity 0, pre-compensated: the AO error budget (fitting error +
            # point-ahead anisoplanatism), phase-only and mean-only.
            terms.append(uplink_fitting_term(scenario, geometry,
                                             cn2_profile=cn2_profile))
            terms.append(uplink_point_ahead_term(scenario, geometry,
                                                 cn2_profile=cn2_profile))
        elif fidelity == 0:
            raise ValueError(
                "fidelity=0 has no analytic mean-only model for an UNCORRECTED "
                "uplink (beam wander plus scintillation has no closed form). Use "
                "fidelity=1 (coupled-flux Monte Carlo) or fidelity=2 (wave "
                "optics)."
            )
        else:
            # Fidelity 1, uncorrected: the coupled-flux Monte Carlo. It carries
            # the tracking jitter, so there is no standalone pointing Term.
            terms.append(uplink_turbulence_term(scenario, geometry,
                                                n_samples=n_samples,
                                                cn2_profile=cn2_profile))
    return Budget(terms, scenario=scenario)


if __name__ == '__main__':
    from ..scenario import SpaceScenario, Channel
    from ..geometry import CircularOrbit
    from ..terminal import Terminal, Transmitter, Aperture, TipTilt, AO

    def _uplink(w0, *, divergence=None, power=None, jitter=0.0,
                ground_aperture=0.5, ground_obscuration=0.0,
                space_aperture=0.05, sensitivity=None, compensation=None,
                precompensation=None):
        '''Build an uplink SpaceScenario: tx=ground, rx=space (satellite).'''
        detector = None if sensitivity is None else Aperture(sensitivity_dbm=sensitivity)
        return SpaceScenario(
            ground=Terminal(aperture_m=ground_aperture, obscuration_ratio=ground_obscuration,
                            wavelength_m=1550e-9, pointing_jitter_rad=jitter,
                            transmitter=Transmitter(waist_m=w0, power_dbm=power,
                                                    divergence_rad=divergence),
                            compensation=compensation or []),
            space=Terminal(aperture_m=space_aperture, wavelength_m=1550e-9,
                           detector=detector),
            direction="uplink", channel=Channel(altitude_m=600e3),
            precompensation=precompensation)

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
    from ..units import w0_to_div
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

    # Dios assumes an untruncated Gaussian launch beam. A central obscuration on
    # the launch aperture flags a violation; a clean launch does not.
    obsc_scn = _uplink(0.1, ground_aperture=0.15, ground_obscuration=0.3)
    obsc_term = uplink_turbulence_term(obsc_scn, geom, n_samples=500, cn2_profile=weak_cn2)
    assert any("untruncated Gaussian" in v for v in obsc_term.assumptions.violations), \
        obsc_term.assumptions.violations
    clean_term = uplink_turbulence_term(scenario, geom, n_samples=500, cn2_profile=weak_cn2)
    assert not any("untruncated Gaussian" in v for v in clean_term.assumptions.violations)

    # --- uplink budget self-check -------------------------------------------
    # A wide launch aperture (1.5 m for a 0.2 m waist) leaves the beam untruncated,
    # so the transmit Gaussian-efficiency term does not fire.
    budget_scn = _uplink(0.2, power=40, jitter=2e-6, sensitivity=-40,
                         ground_aperture=1.5)
    budget_geom = CircularOrbit(altitude_m=600e3, elevation_deg=60.0)
    up = uplink_budget(budget_scn, budget_geom,
                       cn2_profile=default_cn2_profile(budget_scn.channel.site))
    # With turbulence on there is NO separate pointing Term: geometric,
    # atmospheric, turbulence. The jitter lives inside the turbulence Term.
    assert up.to_frame().shape[0] == 3, up.to_frame().shape
    assert not any(t.category == "pointing" for t in up.terms)
    up_mc = up.monte_carlo(2000, rng=np.random.default_rng(0), availabilities=(0.99,))
    up_margin = up_mc["margin_db"][0.99]
    assert np.isfinite(up_margin), up_margin

    # Jitter is not lost: a larger tracking jitter deepens the turbulence Term's
    # loss, and so costs budget margin, WITHOUT any standalone pointing Term.
    calm_scn = _uplink(0.2, power=40, jitter=0.0, sensitivity=-40, ground_aperture=1.5)
    up_calm = uplink_budget(calm_scn, budget_geom,
                            cn2_profile=default_cn2_profile(calm_scn.channel.site))
    jitt_scn = _uplink(0.2, power=40, jitter=2e-6, sensitivity=-40, ground_aperture=1.5)
    up_jitt = uplink_budget(jitt_scn, budget_geom,
                            cn2_profile=default_cn2_profile(jitt_scn.channel.site))
    turb_calm = next(t for t in up_calm.terms if t.category == "turbulence")
    turb_jitt = next(t for t in up_jitt.terms if t.category == "turbulence")
    assert np.isfinite(turb_jitt.mean_db) and np.isfinite(turb_calm.mean_db)
    assert turb_jitt.mean_db > turb_calm.mean_db, (turb_jitt.mean_db, turb_calm.mean_db)

    # With turbulence OFF the standalone pointing Term returns, so jitter is
    # never silently dropped.
    up_noturb = uplink_budget(budget_scn, budget_geom, turbulence=False)
    assert any(t.category == "pointing" for t in up_noturb.terms)

    # A narrow launch aperture (0.15 m for a 0.2 m waist) truncates the beam, so
    # the transmit Gaussian-efficiency term fires.
    ap_scn = _uplink(0.2, power=40, jitter=2e-6, sensitivity=-40,
                     ground_aperture=0.15, ground_obscuration=0.3)
    up_ap = uplink_budget(ap_scn, budget_geom,
                          cn2_profile=default_cn2_profile(ap_scn.channel.site))
    assert up_ap.to_frame().shape[0] == 4, up_ap.to_frame().shape
    eff = next(t for t in up_ap.terms if t.category == "system")
    assert eff.mean_db > 0                       # truncation is a loss
    assert up_ap.total_loss_db() > up.total_loss_db()   # aperture truncation costs margin

    # --- point-ahead anisoplanatism self-check -------------------------------
    pa_geom = CircularOrbit(altitude_m=600e3, elevation_deg=60.0)
    ao_scn = _uplink(0.2, power=40, jitter=2e-6, sensitivity=-40,
                     ground_aperture=1.5, compensation=[TipTilt(), AO(60)])
    pa_term = uplink_point_ahead_term(ao_scn, pa_geom)
    assert np.isfinite(pa_term.mean_db) and pa_term.mean_db > 0, pa_term.mean_db
    # The Term is mean-only: no sampler and no quantile, so it carries no fade.
    assert not pa_term.stochastic and pa_term.quantile is None
    assert pa_term.meta["theta_paa_rad"] > 0
    assert pa_term.category == "anisoplanatism"
    # AO(60) fills radial orders up to n=9 (55 modes through order 9), so the
    # Term corrects orders 2..9, not the infinite-order limit.
    assert pa_term.meta["max_order"] == max_radial_order(60) == 9
    # A real AO order is set, so no upper-bound flag. But the scintillation gap
    # always flags, so the Term is not "ok".
    assert not any("upper bound" in v for v in pa_term.assumptions.violations)
    assert any("NO SCINTILLATION" in v for v in pa_term.assumptions.violations)

    # The error grows with the corrected AO order, up to the infinite-order limit
    # (Fig. 2 of Stone et al. 1994). More corrected orders inject more error.
    v6 = uplink_point_ahead_term(
        _uplink(0.2, ground_aperture=1.5, compensation=[AO(6)]), pa_geom).mean_db
    v60 = pa_term.mean_db
    v_inf = uplink_point_ahead_term(ao_scn, pa_geom, max_order=None).mean_db
    assert v6 < v60 < v_inf, (v6, v60, v_inf)

    # The error grows with the aperture at a fixed order (Eq. 29 of Stone et al.).
    small_scn = _uplink(0.2, ground_aperture=0.3, compensation=[AO(60)])
    large_scn = _uplink(0.2, ground_aperture=1.0, compensation=[AO(60)])
    small_term = uplink_point_ahead_term(small_scn, pa_geom)
    large_term = uplink_point_ahead_term(large_scn, pa_geom)
    assert large_term.mean_db > small_term.mean_db, (large_term.mean_db,
                                                     small_term.mean_db)
    # A small residual gets no extended-Marechal flag.
    assert small_term.meta["sigma2_rad2"] <= MARECHAL_SIGMA2_MAX, \
        small_term.meta["sigma2_rad2"]
    assert not any("MARECHAL" in v for v in small_term.assumptions.violations)

    # A direct call with no AO stage falls back to the infinite-order upper bound.
    tt_scn = _uplink(0.2, power=40, jitter=2e-6, sensitivity=-40,
                     ground_aperture=1.5, compensation=[TipTilt()])
    tt_term = uplink_point_ahead_term(tt_scn, pa_geom)
    assert tt_term.meta["max_order"] is None
    assert any("upper bound" in v for v in tt_term.assumptions.violations), \
        tt_term.assumptions.violations

    # --- AO fitting error (uncorrected high orders, Noll) --------------------
    # Correcting more modes leaves less uncorrected fitting error. An empty stack
    # is the total uncorrected phase variance.
    fit_none = uplink_fitting_term(
        _uplink(0.2, ground_aperture=1.5, compensation=[]), pa_geom)
    fit_ao6 = uplink_fitting_term(
        _uplink(0.2, ground_aperture=1.5, compensation=[AO(6)]), pa_geom)
    fit_ao60 = uplink_fitting_term(ao_scn, pa_geom)
    assert not fit_none.stochastic and fit_none.quantile is None   # mean-only
    assert fit_none.mean_db > fit_ao6.mean_db > fit_ao60.mean_db > 0, (
        fit_none.mean_db, fit_ao6.mean_db, fit_ao60.mean_db)
    assert fit_none.category == "fitting"
    # The fitting Term flags the missing scintillation too.
    assert any("NO SCINTILLATION" in v for v in fit_ao60.assumptions.violations)
    # The uncompensated stack leaves a residual far past the extended-Marechal
    # limit, so the Term flags it (T. S. Ross, DOI 10.1364/AO.48.001812).
    assert fit_none.meta["sigma2_fit_rad2"] > MARECHAL_SIGMA2_MAX
    assert any("MARECHAL" in v for v in fit_none.assumptions.violations), \
        fit_none.assumptions.violations

    # --- source-driven budget dispatch ---------------------------------------
    pa_cn2 = default_cn2_profile(ao_scn.channel.site)

    # No source: the uplink is uncorrected. The coupled-flux turbulence Term is
    # present and there is no anisoplanatism Term.
    uncorr = uplink_budget(ao_scn, pa_geom, n_samples=500, cn2_profile=pa_cn2)
    assert any(t.category == "turbulence" for t in uncorr.terms)
    assert not any(t.category == "anisoplanatism" for t in uncorr.terms)

    # DownlinkBeacon + AO: the uplink is pre-compensated. The coupled-flux Term is
    # REPLACED by the AO error budget -- the fitting error (uncorrected orders)
    # plus the point-ahead anisoplanatism (corrected orders). A standalone
    # pointing Term carries the jitter that the coupled-flux Term used to hold.
    beacon_scn = _uplink(0.2, power=40, jitter=2e-6, sensitivity=-40,
                         ground_aperture=1.5, compensation=[TipTilt(), AO(60)],
                         precompensation=DownlinkBeacon())
    precomp = uplink_budget(beacon_scn, pa_geom, cn2_profile=pa_cn2,
                            fidelity=0)
    assert any(t.category == "anisoplanatism" for t in precomp.terms)
    assert any(t.category == "fitting" for t in precomp.terms)          # Noll piece
    assert not any(t.category == "turbulence" for t in precomp.terms)   # replaced
    assert any(t.category == "pointing" for t in precomp.terms)         # jitter kept

    # The anisoplanatism Term is mean-only, so the budget locks to fidelity 0.
    pa_row = next(t for t in precomp.terms if t.category == "anisoplanatism")
    assert pa_row.mean_db > 0
    assert pa_row.mean_only and not precomp.provides_fade
    try:
        precomp.fade_margin_db(0.99)
    except ValueError as e:
        assert "fidelity-0" in str(e) and "mean-only" in str(e)
    else:
        raise AssertionError("a mean-only budget must refuse fade_margin_db")
    assert np.isfinite(precomp.total_loss_db())   # the mean total is still reported
    # The missing scintillation is a flagged violation, so Budget.check() warns.
    assert any("NO SCINTILLATION" in reason
               for _, reason in precomp.check(warn=False))

    # fidelity takes 0, 1, or 2 only.
    try:
        uplink_budget(beacon_scn, pa_geom, cn2_profile=pa_cn2, fidelity=3)
    except ValueError as e:
        assert "fidelity must be 0, 1, or 2" in str(e)
    else:
        raise AssertionError("an unknown fidelity must raise")

    # --- fidelity-1 pre-compensated uplink (fidelity=1) ----------------------
    # This needs the optional `fast-aosim` package, so skip it when it is
    # absent. ONE FAST Monte-Carlo Term replaces the two analytic phase Terms,
    # and it carries a real fade. The pointing Term stays, because the FAST Term
    # holds no mechanical jitter.
    try:
        precomp_fast = uplink_budget(beacon_scn, pa_geom, cn2_profile=pa_cn2,
                                     n_samples=400, fidelity=1)
    except ImportError as e:
        precomp_fast = None
        print(f"`fast-aosim` unavailable ({e.__class__.__name__}); skipping the "
              "fidelity=1 pre-compensated check.")
    if precomp_fast is not None:
        turb_terms = [t for t in precomp_fast.terms if t.category == "turbulence"]
        assert len(turb_terms) == 1, [t.name for t in turb_terms]
        assert not any(t.category == "anisoplanatism" for t in precomp_fast.terms)
        assert not any(t.category == "fitting" for t in precomp_fast.terms)
        assert any(t.category == "pointing" for t in precomp_fast.terms)
        assert precomp_fast.provides_fade
        fast_margin = precomp_fast.fade_margin_db(0.99)
        assert np.isfinite(fast_margin), fast_margin
        print(f"\npre-compensated uplink, FAST (D=1.5 m, 60 deg): turbulence "
              f"{turb_terms[0].mean_db:.2f} dB mean, 99% budget margin "
              f"{fast_margin:.2f} dB")
        print(f"  ({turb_terms[0].note})")

    # DownlinkBeacon with only a tip-tilt stage corrects no order above the tilt,
    # so the uplink stays uncorrected (coupled flux), no anisoplanatism Term.
    tt_beacon_scn = _uplink(0.2, power=40, jitter=2e-6, sensitivity=-40,
                            ground_aperture=1.5, compensation=[TipTilt()],
                            precompensation=DownlinkBeacon())
    tt_beacon = uplink_budget(tt_beacon_scn, pa_geom, n_samples=500,
                              cn2_profile=pa_cn2)
    assert any(t.category == "turbulence" for t in tt_beacon.terms)
    assert not any(t.category == "anisoplanatism" for t in tt_beacon.terms)

    # LaserGuideStar: not modelled yet, so the budget raises.
    lgs_scn = _uplink(0.2, ground_aperture=1.5, compensation=[AO(60)],
                      precompensation=LaserGuideStar())
    try:
        uplink_budget(lgs_scn, pa_geom, cn2_profile=pa_cn2)
    except NotImplementedError:
        pass
    else:
        raise AssertionError("a laser-guide-star source must raise")

    # --- fidelity-2 wave-optics uplink ---------------------------------------
    # Guards (no run): fidelity=0 uncorrected raises; fidelity=2 pre-compensated
    # raises; fidelity=2 with no bundle raises.
    try:
        uplink_budget(budget_scn, budget_geom, fidelity=0,
                      cn2_profile=default_cn2_profile(budget_scn.channel.site))
    except ValueError as e:
        assert "no analytic mean-only model for an UNCORRECTED" in str(e)
    else:
        raise AssertionError("fidelity=0 uncorrected must raise")
    try:
        uplink_budget(beacon_scn, pa_geom, cn2_profile=pa_cn2, fidelity=2,
                      wave=object())
    except ValueError as e:
        assert "PRE-COMPENSATED" in str(e)
    else:
        raise AssertionError("fidelity=2 pre-compensated must raise")
    try:
        uplink_budget(budget_scn, budget_geom, fidelity=2)
    except ValueError as e:
        assert "needs a precomputed `wave` bundle" in str(e)
    else:
        raise AssertionError("fidelity=2 without a bundle must raise")
    # The default route is UNCHANGED: fidelity=1 uncorrected is the coupled-flux
    # Monte Carlo, not wave optics.
    default_turb = next(t for t in up.terms if t.category == "turbulence")
    assert default_turb.meta.get("model") != "waveoptics"
    assert default_turb.name == "turbulence (coupled-flux)"

    # A real fidelity-2 uncorrected uplink (skip if aotools absent): vacuum +
    # reciprocity Terms, the standalone pointing Term kept, a real fade.
    from ..models.waveoptics import run_fidelity2
    try:
        wo_bundle = run_fidelity2(
            budget_scn, budget_geom, preset="rapid", n_trials=16, seed=9,
            cn2_profile=default_cn2_profile(budget_scn.channel.site))
        wo_up = uplink_budget(
            budget_scn, budget_geom, fidelity=2, wave=wo_bundle,
            cn2_profile=default_cn2_profile(budget_scn.channel.site))
    except ImportError:
        wo_up = None
        print("aotools not installed; skipping the uplink fidelity-2 run.")
    if wo_up is not None:
        vac = next(t for t in wo_up.terms if t.meta.get("model") == "waveoptics-vacuum")
        turb = next(t for t in wo_up.terms if t.meta.get("model") == "waveoptics")
        assert not vac.stochastic and turb.stochastic
        assert "geometric spreading" not in [t.name for t in wo_up.terms]
        # The reciprocity Term holds no jitter, so the pointing Term fires.
        assert any(t.category == "pointing" for t in wo_up.terms)
        assert wo_up.provides_fade and np.isfinite(wo_up.fade_margin_db(0.9))
        print(f"uplink fidelity 2 (600 km, 60 deg, rapid, 16 trials): vacuum "
              f"{vac.mean_db:.2f} dB + turbulence {turb.mean_db:.2f} dB")

    print('\n' + '=' * 40)
    print(f"point-ahead angle: {pa_term.meta['theta_paa_rad'] * 1e6:.2f} urad, "
          f"sigma2={pa_term.meta['sigma2_rad2']:.3f} rad^2")
    print(f"point-ahead loss (D=1.5 m, 60 deg): AO(6) {v6:.2f} dB  |  "
          f"AO(60) {v60:.2f} dB  |  ideal {v_inf:.2f} dB")
    print(f"  ({pa_term.note})")

    # AO error budget (D=1.5 m, 60 deg): the corrected wavefront is the fitting
    # error (uncorrected orders, Noll) plus the point-ahead anisoplanatism
    # (corrected orders, Stone). More modes -> less fitting, more anisoplanatism.
    print("\nAO error budget (D=1.5 m, 60 deg): fitting + point-ahead = total")
    for j in (6, 20, 60, 200):
        scn_j = _uplink(0.2, ground_aperture=1.5, compensation=[AO(j)])
        f_db = uplink_fitting_term(scn_j, pa_geom).mean_db
        a_db = uplink_point_ahead_term(scn_j, pa_geom).mean_db
        print(f"  AO({j:>3}): fitting {f_db:5.2f} dB + point-ahead {a_db:5.2f} dB"
              f" = {f_db + a_db:5.2f} dB")

    print('\n' + '=' * 40)
    print(up.to_frame().to_string(index=False))
    print(f"\nuplink 60 deg 99% margin: {up_margin:.2f} dB")
    print(f"with aperture: +{eff.mean_db:.2f} dB transmit truncation ({eff.note})")
    print(f"fast_available={fast_available}")
    print("self-check passed.")
