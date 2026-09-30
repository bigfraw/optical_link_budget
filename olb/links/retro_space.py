'''
Retroreflected ground-to-space budget: retroreflection as a retransmission.

The crucial idea is that a retroreflector RE-TRANSMITS the beam. The ground
station launches a beam up. The satellite retroreflector captures the power over
its aperture. It then re-emits that captured power back down as a new beam. So a
retro link is an up-leg transmission followed by a down-leg transmission, with
the retro aperture as the hinge: it is the up-leg receive aperture and the
down-leg transmit aperture.

Because every Term is a dB loss and the losses add, the return power is the
retransmission chain:

    P_return = P_launch
             - (up-leg losses)     # fraction of launch power that hits the retro
             - retro_loss
             - (down-leg losses)   # fraction of re-emitted power that hits ground

The up-leg geometric Term already gives the fraction of the launch power that the
retro aperture captures. The down-leg geometric Term treats the retro as a fresh
Gaussian transmitter (waist = half the aperture diameter). So the chained dB sum
IS the retransmission model; no explicit power bookkeeping is needed.

SPACE ONLY. This retransmission picture holds for a ground-to-space link, where
the slant range is long, the return beam diverges far past the ground aperture,
and the two legs see independent turbulence. It does NOT hold for a short
terrestrial (horizontal-path) retro link, where the return does not fully
diverge and reciprocity couples the legs. That case needs a different module.
'''

from dataclasses import replace

import numpy as np

from ..results import Budget, Term
from ..assumptions import (Assumptions, merge_assumptions, BEAM_GAUSSIAN,
                          BEAM_PLANE_WAVE, REGIME_NA, SPECTRUM_NA)
from ..models.geometric import geometric_loss_term
from ..models.gaussian_efficiency import (uniform_aperture_correction_db,
                                          tx_gaussian_efficiency_term)
from ..models.extinction import slant_extinction_term, DEFAULT_TAU_ZENITH
from ..models.pointing import pointing_loss_term
from scipy.special import j1

from ..terminal import Terminal, Transmitter, SpoiledCornerCube
from ..turbulence.profiles import default_cn2_profile
from .uplink import uplink_turbulence_term, TX_TRUNCATION_MIN_DB
from .downlink import downlink_scintillation_term

# The first zero of J1: past it, the station sits outside the Airy core.
AIRY_FIRST_NULL_X = 3.8317


def retro_velocity_aberration_term(scenario, geometry, point_ahead_rad="geometry"):
    '''
    The loss of the return lobe that the velocity aberration moves off the station.

    A corner cube returns the beam along the incoming ray in ITS frame. The
    satellite moves, so in the ground frame the return leaves at the
    aberration angle theta = 2 v_perp / c. That is the point-ahead angle
    (J. J. Degnan, Geodynamics Series 25, 133 (1993), DOI 10.1029/GD025p0133).
    The far field of the unobscured cube aperture D is the Airy pattern, so the
    station reads the fraction [2 J1(x) / x]^2 of the on-axis peak, with
    x = pi D theta / lambda (M. Born and E. Wolf, Principles of Optics, 7th ed.,
    Sec. 8.5.2, DOI 10.1017/CBO9781139644181). The on-axis spread of the lobe is
    the down-leg geometric Term; this Term is the offset only.

    Parameters:
        scenario : SpaceScenario
            The retro link. `space.aperture_m` is the cube aperture.
        geometry : CircularOrbit
            Gives point_ahead_rad when point_ahead_rad="geometry".
        point_ahead_rad : "geometry" or float
            The aberration angle [rad]. A float overrides the geometry (the PAA
            case); 0 gives a 0 dB Term.

    Returns:
        Term
            Category "geometric", deterministic.
    '''
    theta = (geometry.point_ahead_rad if isinstance(point_ahead_rad, str)
             and point_ahead_rad == "geometry" else point_ahead_rad)
    theta = np.asarray(theta, dtype=float)
    D, lam = scenario.space.aperture_m, scenario.space.wavelength_m
    x = np.pi * D * theta / lam
    xs = np.where(x > 0, x, 1.0)
    airy = np.where(x > 0, (2 * j1(xs) / xs) ** 2, 1.0)   # -> 1 as x -> 0
    with np.errstate(divide="ignore"):
        loss_db = -10 * np.log10(airy)
    assumptions = Assumptions(
        beam_type=BEAM_PLANE_WAVE, turbulence_regime=REGIME_NA,
        spectrum=SPECTRUM_NA,
        validity="A standard corner cube with an Airy far field (a circular, "
                 "unobscured aperture, a flat incident wavefront, no dihedral "
                 "offset, no polarisation split). The aberration angle is the "
                 "point-ahead angle 2 v_perp / c. Monostatic: the offset is "
                 "measured from the transmitter.",
    )
    if np.any(x >= AIRY_FIRST_NULL_X):
        assumptions.flag(
            f"PAST THE AIRY NULL: x = pi D theta / lambda = {float(np.max(x)):.2f} "
            f">= {AIRY_FIRST_NULL_X}, so the station sits outside the Airy core. "
            "The real cube pattern (hexagonal pupil, polarisation) differs most "
            "there, so the loss is indicative only. A spoiled cube is the fix.")
    return Term(name="velocity aberration", category="geometric",
                mean_db=loss_db if loss_db.ndim else float(loss_db),
                note="return lobe offset by the velocity aberration 2 v_perp / c",
                meta={"point_ahead_rad": theta, "x": x},
                assumptions=assumptions)


def retro_space_budget(scenario, geometry, *, fidelity=1, turbulence=True,
                       tau_zenith=None, n_samples=3000, cn2_profile=None,
                       retro_loss_db=0.0, fast_params=None,
                       point_ahead_rad="geometry", wave=None):
    '''
    Assemble the retroreflected ground-to-space budget as a retransmission.

    Retroreflection is modelled as a retransmission: the up-leg carries the
    launch power to the retro aperture, and the down-leg re-emits the captured
    power back to the ground receiver. At fidelity 0 and 1 the two legs use
    independent turbulence; fidelity 2 correlates them (see `fidelity`).
    The retroreflector aperture is the hinge: it is the up-leg receive aperture
    and the down-leg transmit aperture. The losses are dB, so the Terms add and
    the sum is the retransmission chain (see the module docstring).

    This is the SPACE model. It assumes a long slant range, a fully diverged
    return, and independent turbulence on the two legs. Do not use it for a short
    terrestrial retro link.

    Parameters:
        scenario : SpaceScenario
            The link case. The `space` Terminal is the passive retroreflector;
            its aperture_m is the retro aperture. The direction is "retro".
        geometry : CircularOrbit or TLEPass
            The link geometry.
        turbulence : bool
            Add the up-leg coupled-flux turbulence Term when true.
        tau_zenith : float, optional
            Zenith optical depth. Defaults to extinction.DEFAULT_TAU_ZENITH.
        n_samples : int
            Monte Carlo draws for the turbulence Term mean estimate.
        cn2_profile : numpy.ndarray, optional
            Explicit zenith Cn2 profile. Defaults to default_cn2_profile.
        retro_loss_db : float
            Fixed loss of the retroreflection [dB].
        fidelity : int
            0, 1 (the default) or 2. At 0 and 1 it is the down-leg
            receive-coupling fidelity: 1 is the FAST modal overlap, 0 the
            analytic mean-only coupling. The UP-leg turbulence stays the
            coupled-flux Monte Carlo at either value (there is no analytic
            mean-only uncorrected uplink model). The two legs are INDEPENDENT
            there, which under-reads the fade of a fibre return (the Term
            flags it; validation/retro_bracket). At 2 (wave optics) ONE
            stochastic Term holds BOTH legs of each trial, from a retro wave
            record (`wave`): the up leg and the return of one pulse share one
            line, and the wind moves the air between them over 2R/c
            (olb.waveoptics.turbulence.run.retro_sensing_geometry). Fidelity 2
            takes a bucket (None or Aperture) or an SMF receiver.
        wave : Fidelity2Bundle or Campaign, optional
            The retro wave record of fidelity 2, from
            olb.models.waveoptics.run_fidelity2 on THIS retro scenario, or a
            retro Campaign. The space geometric loss is analytic, so the
            bundle holds no vacuum record.
        fast_params : dict, optional
            Extra FAST parameters for the fidelity-1 down-leg coupling.
        point_ahead_rad : "geometry" or float
            The velocity-aberration angle of the return [rad] (the PAA case).
            "geometry" reads geometry.point_ahead_rad; a float overrides it.
            The kind of retro is `scenario.space.retroreflector` (None is the
            standard CornerCube; a SpoiledCornerCube raises, not built).

    Returns:
        Budget
            The budget with the original scenario set.

    Raises:
        ValueError
            If fidelity is not 0, 1 or 2, or a fidelity-2 call has no retro
            wave record, a vacuum record, or more than one line of sight.
        NotImplementedError
            If the retroreflector is a SpoiledCornerCube, or a fidelity-2
            receiver is an MMF or a Camera.
    '''
    if isinstance(scenario.space.retroreflector, SpoiledCornerCube):
        raise NotImplementedError(
            "a SpoiledCornerCube (dihedral-angle offset) is not built. Its "
            "return is six beams on a ring, not one Airy lobe. Use CornerCube.")
    if fidelity not in (0, 1, 2):
        raise ValueError(f"fidelity must be 0, 1 or 2 for retro, got "
                         f"{fidelity!r}.")
    smf_fidelity = "fast" if fidelity == 1 else "mean"
    retro_aperture_m = scenario.space.aperture_m
    wavelength = scenario.space.wavelength_m
    tau = DEFAULT_TAU_ZENITH if tau_zenith is None else tau_zenith

    # The retro is the uplink receiver. A corner-cube retro has no central
    # obscuration. The up-leg keeps the ground transmit terminal.
    retro_rx = Terminal(aperture_m=retro_aperture_m, obscuration_ratio=0.0,
                        wavelength_m=wavelength)
    up_scn = replace(scenario, direction="uplink", space=retro_rx)
    # The retro re-transmits the captured power: it is the plane-wave transmitter
    # on the return. The Gaussian-equivalent waist is half the aperture diameter.
    # The retro is passive, so there is no active pointing jitter. The ground
    # stays the receiver.
    retro_tx = Terminal(aperture_m=retro_aperture_m, obscuration_ratio=0.0,
                        wavelength_m=wavelength,
                        transmitter=Transmitter(waist_m=retro_aperture_m / 2.0))
    down_scn = replace(scenario, direction="downlink", space=retro_tx)

    if cn2_profile is None:
        cn2_profile = default_cn2_profile(scenario.channel.site)

    up_terms = [
        geometric_loss_term(up_scn, geometry),
        slant_extinction_term(up_scn, geometry, tau_zenith=tau),
    ]
    # Up-leg pointing jitter folds into the coupled-flux turbulence Term (it
    # shares the beam-wander displacement), so add the standalone pointing Term
    # only when the turbulence Term is off. Adding both double-counts the jitter.
    # See olb.links.uplink.uplink_budget. The wave-optics Term holds no jitter,
    # so fidelity 2 always keeps it.
    if not turbulence or fidelity == 2:
        up_terms.append(pointing_loss_term(up_scn, geometry))
    # The launch aperture truncates the up-leg beam, exactly as uplink_budget
    # does. Opt-in: it fires only when the ground transmitter truncates the beam
    # by more than TX_TRUNCATION_MIN_DB. A bistatic ground reads its beam-director
    # aperture through the Transmitter override, not the receive-telescope
    # aperture (see olb.terminal.Transmitter).
    up_tx = up_scn.tx_terminal
    if up_tx.transmitter is not None:
        eff = tx_gaussian_efficiency_term(up_scn, geometry)
        if eff.mean_db > TX_TRUNCATION_MIN_DB:
            up_terms.append(eff)
    if turbulence and fidelity != 2:
        up_terms.append(uplink_turbulence_term(up_scn, geometry, n_samples=n_samples,
                                               cn2_profile=cn2_profile))
    for t in up_terms:
        t.name = "uplink " + t.name

    # The retro reflects the flat wavefront that fills its aperture, so the
    # return is a uniformly illuminated (top-hat) aperture, not a Gaussian. The
    # geometric Term models a Gaussian of waist aperture/2, which over-states the
    # on-axis gain by 2/(1-Cr^2). This fixed Term converts it to the top-hat.
    tophat_db = uniform_aperture_correction_db(retro_tx.obscuration_ratio)
    tophat_term = Term(
        name="top-hat correction", category="system", mean_db=tophat_db,
        note="retro return is a uniform aperture, not a Gaussian(waist=D/2)",
        assumptions=Assumptions(
            beam_type=BEAM_PLANE_WAVE, turbulence_regime=REGIME_NA,
            spectrum=SPECTRUM_NA,
            validity="Converts the Gaussian(waist=aperture/2) geometric model to "
                     "a uniformly illuminated (top-hat) aperture, via the "
                     "Gaussian aperture-illumination efficiency. Far-field, "
                     "on-axis (receive aperture much smaller than the return "
                     "lobe).",
        ),
    )
    down_terms = [
        geometric_loss_term(down_scn, geometry),
        slant_extinction_term(down_scn, geometry, tau_zenith=tau),
        tophat_term,
        retro_velocity_aberration_term(scenario, geometry, point_ahead_rad),
    ]
    # The return-leg receive term follows the ground receiver, exactly as
    # downlink_budget does: a bucket receiver (no detector or a plain Aperture) is
    # phase-insensitive and gets the standalone plane-wave scintillation. None and
    # Aperture() are the SAME bucket. An SMF detector adds the fibre-coupling loss.
    rx = down_scn.rx_terminal
    detector = rx.detector if rx is not None else None
    from ..terminal import Aperture, SMF
    wave_terms = []
    if fidelity == 2:
        wave_terms = _retro_fidelity2_terms(scenario, geometry, wave,
                                            turbulence)
    elif detector is not None and not isinstance(detector, Aperture):
        # Import here to break the downlink <-> coupling import cycle.
        from ..models.coupling import downlink_coupling_term
        cpl = downlink_coupling_term(down_scn, geometry, n_samples=n_samples,
                                     smf_fidelity=smf_fidelity,
                                     fast_params=fast_params)
        if isinstance(detector, SMF) and turbulence:
            cpl.assumptions.flag(
                "INDEPENDENT LEGS: the up leg and a FIBRE return fade TOGETHER "
                "(they share one line through the air), and this rung draws "
                "them apart, so it under-reads the fade (2.5 to 9 dB at p5 on "
                "the hero 0.7 m case, validation/retro_bracket). Use "
                "fidelity=2 for the fade.")
        down_terms.append(cpl)
    else:
        down_terms.append(downlink_scintillation_term(down_scn, geometry))
    for t in down_terms:
        t.name = "downlink " + t.name

    retro_term = Term(
        name="retro reflection", category="system", mean_db=retro_loss_db,
        assumptions=Assumptions(
            beam_type=BEAM_PLANE_WAVE, turbulence_regime=REGIME_NA,
            spectrum=SPECTRUM_NA,
            validity="Retroreflection is modelled as a retransmission: the retro "
                     "captures the up-leg power and re-emits it down. This holds "
                     "for a ground-to-space link (long range, fully diverged "
                     "return, independent turbulence on the two legs), not for a "
                     "short terrestrial link. The reflected wavefront is flat at "
                     "the satellite; the return is a plane wave. The "
                     "retroreflector aperture is modelled as a Gaussian waist of "
                     "half the aperture diameter. The velocity aberration of "
                     "the return is its own Term.",
        ),
    )

    return Budget(up_terms + down_terms + wave_terms + [retro_term],
                  scenario=scenario)


def _retro_fidelity2_terms(scenario, geometry, wave, turbulence):
    """Give the ONE stochastic Term of a fidelity-2 retro budget.

    Each trial of a retro wave record holds both legs of one pulse: the up leg
    as the reciprocity overlap eta_turb (on the wind-shifted window; Shapiro,
    DOI 10.1364/JOSA.61.000492) and the return as the ground receiver
    (collected_power, and smf_eta for a fibre). The loss of the trial is the
    SUM of the two leg losses in dB, so the Term keeps their correlation. Both
    scalars are vacuum-normalised, and smf_eta is the ABSOLUTE fibre coupling,
    exactly as downlink_budget reads it; the analytic geometric Terms carry the
    rest.

    Returns:
        A list with the Term, or an empty list when turbulence is off.
    """
    from ..models.waveoptics import resolve_wave, waveoptics_turbulence_term
    from ..terminal import Aperture, SMF
    detector = scenario.ground.detector
    if detector is not None and not isinstance(detector, (Aperture, SMF)):
        raise NotImplementedError(
            f"the fidelity-2 retro takes a bucket or an SMF receiver, not "
            f"{type(detector).__name__}.")
    if np.asarray(geometry.elevation_deg, dtype=float).size != 1:
        raise ValueError("the fidelity-2 retro takes ONE line of sight. Loop "
                         "over the elevations (olb.sweep).")
    if not turbulence:
        if isinstance(detector, SMF):
            raise ValueError(
                "a fibre receiver at fidelity=2 with turbulence=False has no "
                "coupling number (the space geometric loss is analytic). Use "
                "fidelity 0 or 1, or an Aperture receiver.")
        return []
    wave = resolve_wave(wave)
    rec = None if wave is None else wave.turbulent
    if rec is None or getattr(rec, "retro_wind_ground_m_s", None) is None:
        raise ValueError(
            "fidelity=2 needs a RETRO wave record: run "
            "olb.models.waveoptics.run_fidelity2 on this retro scenario, or "
            "pass a retro Campaign as wave=.")
    if wave.vacuum is not None:
        raise ValueError("the fidelity-2 retro reads the ANALYTIC geometric "
                         "loss; pass a bundle with vacuum=None.")
    up = np.array([t.eta_turb for t in rec.trials], dtype=float)
    ret = np.array([t.collected_power for t in rec.trials], dtype=float)
    if isinstance(detector, SMF):
        ret = ret * np.array([t.smf_eta for t in rec.trials], dtype=float)
    term = waveoptics_turbulence_term(
        rec, loss_db=-10.0 * np.log10(up * ret), beam_type=BEAM_GAUSSIAN,
        name="retro turbulence (wave optics)",
        L0_m=scenario.channel.site.outer_scale_m,
        note="both legs of each trial: the up-leg overlap on the wind-shifted "
             "window times the return at the ground receiver, vacuum-"
             "normalised.",
        meta_extra={"retro_wind_ground_m_s": rec.retro_wind_ground_m_s})
    for reason in (
            "FROZEN FLOW OVER THE ROUND TRIP: the up leg reads the air that the "
            "Bufton wind moved in 2R/c (integer pixels, one wind direction, no "
            "boiling; Taylor, DOI 10.1098/rspa.1938.0032).",
            "NO ENHANCED BACKSCATTER: the model is two passes of one "
            "atmosphere, not the coupled double passage of Andrews and "
            "Phillips Ch. 13 (DOI 10.1117/3.626196).",
            "POINT-SOURCE CUBE: the up leg is the on-axis overlap, so the cube "
            "must be much smaller than the uplink coherence width at the "
            "satellite."):
        term.assumptions.flag(reason, source="retro_space_budget")
    return [term]


if __name__ == '__main__':
    from ..scenario import SpaceScenario, Channel
    from ..terminal import Terminal, Transmitter, Aperture, SMF
    from ..geometry import CircularOrbit

    # The ground terminal transmits up and receives the return. The space
    # terminal is the passive retroreflector (aperture only). The ground is
    # bistatic: a small beam director (0.15 m) transmits the up-leg, so the launch
    # truncation reads the director aperture, not the 0.7 m receive telescope.
    retro_scn = SpaceScenario(
        ground=Terminal(aperture_m=0.7, obscuration_ratio=0.3, wavelength_m=1550e-9,
                        pointing_jitter_rad=0e-6,
                        transmitter=Transmitter(waist_m=0.06, power_dbm=40,
                                                aperture_m=0.15, obscuration_ratio=0.3),
                        detector=Aperture(sensitivity_dbm=-50)),
        space=Terminal(aperture_m=0.05, wavelength_m=1550e-9),
        direction="retro", channel=Channel(altitude_m=1500e3),
    )
    elevation = 30.0
    retro_geom = CircularOrbit(altitude_m=1500e3, elevation_deg=elevation)

    retro = retro_space_budget(retro_scn, retro_geom)
    # 10 terms (9 + the velocity aberration): the up-leg carries the opt-in launch-truncation term because the
    # 0.15 m beam director truncates the 0.06 m waist beam. With turbulence on
    # there is NO standalone up-leg pointing Term (the jitter folds into the
    # coupled-flux turbulence Term).
    assert retro.to_frame().shape[0] == 10, retro.to_frame().shape
    assert not any(t.category == "pointing" for t in retro.terms)
    names = [t.name for t in retro.terms]
    assert "uplink transmit Gaussian efficiency" in names, names
    # The truncation reads the director aperture (0.15 m), not the receive
    # telescope (0.7 m): alpha = (0.15/2)/0.06 = 1.25, an unobscured director.
    eff = next(t for t in retro.terms if t.name == "uplink transmit Gaussian efficiency")
    assert abs(eff.meta["alpha"] - 1.25) < 1e-9, eff.meta["alpha"]
    frame = retro.to_frame()
    n_atmos = (frame["category"] == "atmospheric").sum()
    assert n_atmos == 2, n_atmos
    # The return leg carries the top-hat correction (+3.01 dB, unobscured retro).
    names = [t.name for t in retro.terms]
    assert "downlink top-hat correction" in names, names
    tophat = next(t for t in retro.terms if t.name == "downlink top-hat correction")
    assert abs(tophat.mean_db - 3.0103) < 1e-3, tophat.mean_db
    # The ground has an Aperture detector, a bucket, so the return leg carries the
    # leg-prefixed scintillation Term (a bucket is the SAME as no detector).
    assert "downlink scintillation" in [t.name for t in retro.terms], \
        [t.name for t in retro.terms]
    retro_mc = retro.monte_carlo(2000, rng=np.random.default_rng(0),
                                 availabilities=(0.99,))
    retro_fade = retro_mc["fade_db"][0.99]
    retro_mean = retro_mc["mean_loss_db"]
    assert np.isfinite(retro_fade), retro_fade
    af = retro.assumptions_frame()
    retro_row = af[af["name"] == "retro reflection"]
    assert not retro_row.empty, "retro reflection row missing"
    assert "plane wave" in retro_row.iloc[0]["validity"], retro_row.iloc[0]
    assert "retransmission" in retro_row.iloc[0]["validity"], retro_row.iloc[0]

    # --- assumption trace wiring (WP3b) -------------------------------------
    # retro recomposes finished Terms; it opens NO trace of its own. The up-leg
    # turbulence Term comes from the wired uplink factory, so it carries traced
    # provenance, and the untraced guard must not flag it.
    up_turb = next(t for t in retro.terms
                   if t.name == "uplink turbulence (coupled-flux)")
    assert up_turb.assumptions.provenance, "the up-leg turbulence Term needs provenance"
    assert any("uplink_flux._flux_result" in s
               for s in up_turb.assumptions.provenance), up_turb.assumptions.provenance
    guard = [name for name, reason in retro.check(warn=False)
             if "did not open" in reason]
    # This WP owns the up-leg turbulence Term. (The down-leg scintillation Term
    # is a sibling WP, so do not assert on it here.)
    assert "uplink turbulence (coupled-flux)" not in guard, guard
    # merge_assumptions recomposes finished records with NO trace of its own: the
    # folded up+down record carries the union of the two provenances.
    down_cpl = next(t for t in retro.terms
                    if t.name == "downlink scintillation")
    folded = merge_assumptions(up_turb.assumptions, down_cpl.assumptions)
    assert set(up_turb.assumptions.provenance) <= set(folded.provenance)

    # --- velocity aberration (the PAA case) ----------------------------------
    va = next(t for t in retro.terms if t.name == "downlink velocity aberration")
    x = va.meta["x"]
    assert abs(va.mean_db + 10 * np.log10((2 * j1(x) / x) ** 2)) < 1e-12
    # A 5 cm cube at 1500 km, 30 deg sits inside the Airy core (x ~ 2.4, the
    # D_opt point), so no flag; x past the null is flagged.
    assert 0 < x < AIRY_FIRST_NULL_X and not va.assumptions.violations, x
    far = retro_velocity_aberration_term(retro_scn, retro_geom, 50e-6)
    assert far.meta["x"] > AIRY_FIRST_NULL_X and far.assumptions.violations
    # The PAA override: 0 is a 0 dB Term; a 1 cm cube at 50 urad is ~1.14 dB.
    zero = retro_velocity_aberration_term(retro_scn, retro_geom, 0.0)
    assert zero.mean_db == 0.0 and not zero.assumptions.violations
    small = replace(retro_scn, space=Terminal(aperture_m=0.01))
    one_cm = retro_velocity_aberration_term(small, retro_geom, 50e-6)
    assert abs(one_cm.mean_db - 1.14) < 0.005, one_cm.mean_db
    # The spoiled kind is a stub.
    from ..terminal import SpoiledCornerCube as _Spoiled
    try:
        retro_space_budget(replace(retro_scn, space=Terminal(
            aperture_m=0.05, retroreflector=_Spoiled(5e-6))), retro_geom)
        raise AssertionError("a SpoiledCornerCube must raise")
    except NotImplementedError:
        pass
    print(f"velocity aberration: theta = {va.meta['point_ahead_rad']*1e6:.1f} urad, "
          f"x = {x:.2f}, loss = {va.mean_db:.2f} dB")

    # --- fidelity 2: ONE Term holds both legs of each trial -----------------
    from ..models.waveoptics import run_fidelity2
    wave = run_fidelity2(retro_scn, retro_geom, n_trials=8, seed=5,
                         preset="rapid", progress=False)
    assert wave.vacuum is None and wave.turbulent.retro_wind_ground_m_s == 10.0
    f2 = retro_space_budget(retro_scn, retro_geom, fidelity=2, wave=wave)
    names2 = [t.name for t in f2.terms]
    for want in ("retro turbulence (wave optics)", "uplink pointing",
                 "downlink velocity aberration", "downlink top-hat correction"):
        assert any(want in n for n in names2), (want, names2)
    assert not any("coupled-flux" in n or "scintillation" in n
                   for n in names2), names2
    tw = next(t for t in f2.terms if t.name == "retro turbulence (wave optics)")
    legs = np.array([t.eta_turb * t.collected_power
                     for t in wave.turbulent.trials])
    assert abs(tw.mean_db - np.mean(-10 * np.log10(legs))) < 1e-9, tw.mean_db
    assert any("FROZEN FLOW" in v for v in tw.assumptions.violations)
    # A record of another direction is not a retro record.
    up_wave = run_fidelity2(replace(retro_scn, direction="uplink"), retro_geom,
                            n_trials=2, seed=5, preset="rapid", progress=False)
    try:
        retro_space_budget(retro_scn, retro_geom, fidelity=2, wave=up_wave)
        raise AssertionError("an uplink record must not feed a retro budget")
    except ValueError:
        pass
    # Fidelity 0/1 draw the legs apart: a fibre return says so.
    smf_scn = replace(retro_scn, ground=replace(retro_scn.ground,
                                                detector=SMF()))
    f0 = retro_space_budget(smf_scn, retro_geom, fidelity=0)
    cpl0 = next(t for t in f0.terms if t.category == "coupling")
    assert any("INDEPENDENT LEGS" in v for v in cpl0.assumptions.violations)
    print(f"retro fidelity 2 (8 trials): turbulence {tw.mean_db:.2f} dB mean, "
          f"total {f2.to_frame()['mean_db'].sum():.2f} dB")

    print("retro (space) assumptions:")
    retro.check()
    print("retro (space) budget terms:")
    print(retro.to_frame().to_string(index=False))
    print(f"\nretro {elevation} deg 99% fade: {retro_fade:.2f} dB")
    print(f"retro {elevation} deg mean: {retro_mean:.2f} dB")
    print("self-check passed")
