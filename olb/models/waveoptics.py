'''
Fidelity-2 wave-optics Terms: turn a turbulent field propagation into a Term.

This module is the bridge from the split-step wave-optics layer
(olb.waveoptics.turbulence) to the budget. It has two parts.

`run_waveoptics` runs the turbulent propagation ONE time. It sizes the grid and
plans the screens once (turbulent_grid), then it shares that grid and plan across
every trial. It gives back the minimal scalar record TurbWaveResult.

`waveoptics_turbulence_term` turns that record into a Term. It NEVER runs a
simulation: it reads the per-trial scalars of a record that a caller already
computed. So the expensive propagation runs once, not once for each budget build.
The reduction mirrors the fidelity-1 FAST factory (olb.models.fast): the
per-trial loss makes an empirical mean, an empirical quantile, and a resampling
sampler.

WHICH SCALAR. Each trial holds three scalars, and the direction sets which one is
the loss (see olb.waveoptics.turbulence.run.TurbTrial):
  - smf_eta         the single-mode-fibre coupling efficiency (a terrestrial or a
                    downlink fibre receiver). It is the ABSOLUTE efficiency, so it
                    already holds the static mode-match floor.
  - collected_power the power in the receive aperture. The SPACE case normalises
                    it to the vacuum baseline, so it is the pure turbulence
                    penalty (vacuum limit 1.0). Do NOT use the terrestrial
                    collected_power as a turbulence Term: it holds the geometric
                    spread, which the geometric Term already carries.
  - eta_turb        the uplink reciprocity overlap against the free-space
                    baseline. It is the pure turbulence penalty.
Each scalar is a ratio, so the loss is -10*log10(ratio). NO floor is added: the
three ratios are absolute or vacuum-normalised. (This differs from FAST, where
result.dB_rel is a penalty relative to the diffraction limit and the static floor
adds on top.)

SNAPSHOT-ONLY. The trials are independent atmosphere snapshots with no time axis.
The Term gives a fade-DEPTH distribution, not a fade RATE or a fade DURATION. See
olb.waveoptics.turbulence.temporal for the planned frozen-flow extension.

Sources:
- Schmidt, Numerical Simulation of Optical Wave Propagation, DOI 10.1117/3.866274,
  Ch. 9. The split-step propagation that gives the trials.
- Shapiro, DOI 10.1364/JOSA.61.000492. The uplink-downlink reciprocity of
  eta_turb.
'''

import warnings
from dataclasses import dataclass

import numpy as np

from ..results import Term, EmpiricalSampler
from ..assumptions import (Assumptions, BEAM_GAUSSIAN, BEAM_PLANE_WAVE,
                            REGIME_WEAK, REGIME_STRONG, REGIME_NA,
                            SPECTRUM_KOLMOGOROV, SPECTRUM_VON_KARMAN, SPECTRUM_NA)
from ..turbulence.plane_wave_scintillation import WEAK_FLUCTUATION_LIMIT

# Each quantity fixes the default Term name and category.
_QUANTITY_SPEC = {
    "smf_eta": ("receive coupling (SMF)", "coupling"),
    "collected_power": ("scintillation", "turbulence"),
    "eta_turb": ("turbulence (wave optics, reciprocity)", "turbulence"),
}


@dataclass(frozen=True)
class Fidelity2Bundle:
    '''The two wave-optics records a fidelity-2 budget needs.

    A fidelity-2 budget shows TWO Terms: a deterministic vacuum-optics Term (the
    full no-turbulence loss launch to detector) and a stochastic turbulence Term
    (the fade). They come from two runs, bundled here:

        vacuum    : the WaveResult of one no-turbulence propagation
                    (olb.waveoptics.run.propagate_scenario). It gives the full
                    geometric spread, the launch truncation, the aperture
                    capture, and the vacuum fibre coupling.
        turbulent : the TurbWaveResult of the split-step Monte Carlo
                    (olb.waveoptics.turbulence.propagate_turbulent_scenario). It
                    gives the turbulence penalty.

    Build one with run_fidelity2, which runs both on the correct grids.
    '''
    vacuum: object
    turbulent: object


def run_waveoptics(scenario, geometry, *, n_trials=200, preset="standard",
                   seed=None, threader=None, grid=None, plan=None, hs=None,
                   cn2_profile=None, L0_m=np.inf, subharmonics=True):
    '''
    Run the turbulent split-step propagation ONE time.

    This is the one place the expensive simulation runs. It sizes the grid and
    plans the screens once (turbulent_grid), then it shares that grid and plan
    across every trial, so a longer run repeats a shorter run's atmosphere. Give
    the result to `waveoptics_turbulence_term`, which never runs a simulation.

    Parameters:
        scenario : SpaceScenario or TerrestrialScenario
            The link case. The direction sets which scalar is the loss (see the
            module docstring).
        geometry : CircularOrbit, TLEPass, or HorizontalPath
            The link geometry. It must resolve to ONE range (scalar elevation).
        n_trials : int
            The number of independent atmosphere snapshots.
        preset : str or QualityPreset
            The sampling quality: "reference", "standard", or "rapid".
        seed : int or None
            The seed of the trial set. The same seed repeats the same set.
        threader : olb.waveoptics.Threader or None
            Run the trials across threads when set (the FFT releases the GIL).
        grid, plan : optional
            A precomputed grid and screen plan. Give both or neither. When None,
            turbulent_grid sizes them from the scenario.
        hs, cn2_profile : numpy.ndarray, optional
            The height grid and the zenith Cn2 profile (a space link only).
        L0_m : float
            The turbulence outer scale [m]. Infinite is the Kolmogorov limit.
        subharmonics : bool
            Add the aotools subharmonic low-frequency screen content.

    Returns:
        TurbWaveResult
            The minimal scalar record (olb.waveoptics.turbulence.run).
    '''
    from ..waveoptics.turbulence import (propagate_turbulent_scenario,
                                          turbulent_grid)
    if grid is None or plan is None:
        grid, plan, _ = turbulent_grid(scenario, geometry, preset=preset, hs=hs,
                                       cn2_profile=cn2_profile, L0_m=L0_m)
    return propagate_turbulent_scenario(
        scenario, geometry, n_trials=n_trials, seed=seed, preset=preset,
        grid=grid, plan=plan, hs=hs, cn2_profile=cn2_profile, L0_m=L0_m,
        subharmonics=subharmonics, threader=threader)


def waveoptics_turbulence_term(result, *, quantity=None, loss_db=None,
                               beam_type=BEAM_GAUSSIAN, name=None, category=None,
                               spectrum=None, sigma2_I=None, L0_m=np.inf,
                               note=None, meta_extra=None):
    '''
    Build a fidelity-2 Term from a turbulent wave-optics record.

    This never runs a simulation. It reduces the per-trial loss of a computed
    TurbWaveResult to the three Term faces, the same way the FAST factory does: an
    empirical mean, an empirical quantile, and a resampling sampler. The Term
    carries a real fade (it is not mean-only).

    Give `quantity` (read one trial scalar) OR `loss_db` (a per-trial loss array
    the caller already computed, for a COMPOSITE penalty such as collected_power
    times smf_eta). At least one is required; with `loss_db`, `quantity` only
    labels the meta.

    Parameters:
        result : TurbWaveResult
            The record from run_waveoptics (or propagate_turbulent_scenario). It
            is read for the trial count, the preset, and the seed metadata, even
            when `loss_db` gives the losses.
        quantity : str, optional
            The per-trial scalar to read: "smf_eta", "collected_power", or
            "eta_turb". The module docstring says which one fits each direction.
        loss_db : array-like, optional
            A per-trial loss in dB (positive = loss), one value per trial. Use it
            for a composite penalty the caller computes from more than one scalar.
            When given, `quantity` only labels the meta; `name` and `category`
            default to a generic turbulence Term.
        beam_type : str
            The launch beam type for the assumptions (BEAM_GAUSSIAN for a launched
            Gaussian, BEAM_PLANE_WAVE for a downlink slab).
        name, category : str, optional
            Override the default Term name and category of the quantity.
        spectrum : str, optional
            The turbulence spectrum label. When None, it follows L0_m: Kolmogorov
            for an infinite outer scale, von Karman for a finite one. It must
            match the sibling turbulence Terms, because Budget.check() allows one
            spectrum per budget.
        sigma2_I : float, optional
            The plane-wave scintillation index, for the weak/strong regime flag
            and the meta. It is not on the minimal record, so a caller passes it.
        L0_m : float
            The outer scale [m] that the run used, for the spectrum label and the
            meta. The minimal record does not carry it, so a caller passes it.
        note : str, optional
            The Term note. A default names the model, the quantity, and the count.
        meta_extra : dict, optional
            Extra meta merged into the Term meta.

    Returns:
        Term
            A three-face fidelity-2 Term (mean_only=False).

    Raises:
        ValueError
            If neither quantity nor loss_db is given, if quantity is unknown, if
            the record has no trials, or if the chosen scalar is None for every
            trial (the direction does not set it).
    '''
    if quantity is None and loss_db is None:
        raise ValueError("give quantity, or loss_db (with an optional quantity "
                         "label).")

    if loss_db is not None:
        # The caller computed a per-trial COMPOSITE loss (dB) and hands it in.
        losses = np.asarray(loss_db, dtype=float).ravel()
        if losses.size == 0:
            raise ValueError("loss_db is empty.")
        quantity_label = quantity if quantity is not None else "composite"
        name = "turbulence (wave optics)" if name is None else name
        category = "turbulence" if category is None else category
        loss_desc = "a composite per-trial loss (dB)"
    else:
        if quantity not in _QUANTITY_SPEC:
            raise ValueError(
                f"unknown quantity {quantity!r}. Use one of "
                f"{sorted(_QUANTITY_SPEC)}, or pass loss_db."
            )
        default_name, default_category = _QUANTITY_SPEC[quantity]
        name = default_name if name is None else name
        category = default_category if category is None else category
        raw = [getattr(t, quantity) for t in result.trials]
        if len(raw) == 0:
            raise ValueError("the TurbWaveResult holds no trials.")
        if any(v is None for v in raw):
            raise ValueError(_none_message(quantity))
        # Each scalar is a power ratio (absolute or vacuum-normalised), so the
        # loss is -10*log10(ratio). No static floor is added.
        losses = -10.0 * np.log10(np.asarray(raw, dtype=float))
        quantity_label = quantity
        loss_desc = f"-10*log10({quantity})"

    # The EmpiricalSampler owns the mean, the empirical quantile, the bootstrap
    # sampler, AND the tail-adequacy check (it warns when a quantile reads too few
    # tail samples).
    sampler = EmpiricalSampler(losses)
    mean_db = sampler.mean_db
    n_trials = sampler.n

    if spectrum is None:
        spectrum = SPECTRUM_KOLMOGOROV if np.isinf(L0_m) else SPECTRUM_VON_KARMAN
    scale_note = (
        "The outer scale is infinite (the Kolmogorov limit)." if np.isinf(L0_m)
        else f"The outer scale is L0={L0_m:g} m (von Karman).")

    weak = sigma2_I is None or float(sigma2_I) <= WEAK_FLUCTUATION_LIMIT
    regime = REGIME_WEAK if weak else REGIME_STRONG
    assumptions = Assumptions(
        beam_type=beam_type,
        turbulence_regime=regime,
        spectrum=spectrum,
        validity="Fidelity-2 wave optics: a split-step field propagation through "
                 "independent atmosphere snapshots (olb.waveoptics.turbulence; "
                 "Schmidt 2010, DOI 10.1117/3.866274, Ch. 9). The Term is the "
                 f"per-trial loss {loss_desc}, with no added floor, "
                 "reduced to an empirical mean, quantile, and sampler. " +
                 scale_note +
                 " SNAPSHOT-ONLY: the trials have no time axis, so this Term is a "
                 "fade-DEPTH distribution, not a fade RATE or a fade DURATION.",
    )
    # SNAPSHOT-ONLY is a real limit, so flag it always.
    assumptions.flag(
        "SNAPSHOT-ONLY: the wave-optics trials are independent atmosphere "
        "snapshots with no time axis. This Term gives the fade depth "
        "distribution only, not the fade rate or the fade duration. See "
        "olb.waveoptics.turbulence.temporal for the planned frozen-flow axis."
    )
    # The split-step solver stays valid in strong turbulence. So a strong index is
    # NOT an accuracy break: it only means the deep tail needs more trials.
    if not weak:
        assumptions.flag(
            f"Plane-wave scintillation index sigma2_I={float(sigma2_I):.2f} is "
            f"past the weak-fluctuation limit {WEAK_FLUCTUATION_LIMIT}. The "
            "split-step field solver stays valid, but the fade deepens, so the "
            "deep-tail quantile needs a large n_trials to converge."
        )

    if note is None:
        note = (f"wave-optics {quantity_label}, {n_trials} snapshots, "
                f"preset={result.preset}")
    meta = {
        "model": "waveoptics",
        "quantity": quantity_label,
        "preset": result.preset,
        "n_trials": n_trials,
        "seed_entropy": result.seed_entropy,
        "sigma2_I": None if sigma2_I is None else float(sigma2_I),
        "L0_m": float(L0_m),
        "weak_fluctuation_valid": bool(weak),
        "weak_fluctuation_limit": WEAK_FLUCTUATION_LIMIT,
    }
    if meta_extra:
        meta.update(meta_extra)

    return Term(
        name=name,
        category=category,
        mean_db=mean_db,
        sampler=sampler,             # EmpiricalSampler: bootstrap resample
        quantile=sampler.quantile,   # empirical quantile, warns when under-sampled
        note=note,
        meta=meta,
        assumptions=assumptions,
    )


def waveoptics_smf_coupling_term(result, **kwargs):
    '''
    Build the fidelity-2 turbulent SMF-coupling Term (category "coupling").

    This is the coupling-category face of the wave-optics record: the per-trial
    single-mode-fibre efficiency, reduced to the three Term faces. It is the
    fidelity-2 companion of the fidelity-1 smf_fast_term. It is a thin wrapper on
    waveoptics_turbulence_term with quantity="smf_eta", named so a caller finds
    the coupling Term without knowing the quantity flag; the coupling package
    re-exports it. See that function for the keyword arguments.

    Parameters:
        result : TurbWaveResult
            The record from run_waveoptics. The receive terminal must be an SMF,
            or the smf_eta scalar is None for every trial and this raises.
        **kwargs :
            Passed to waveoptics_turbulence_term (beam_type, name, spectrum,
            sigma2_I, L0_m, note, meta_extra). `quantity` is fixed to "smf_eta".

    Returns:
        Term
            A three-face coupling Term (mean_only=False).
    '''
    return waveoptics_turbulence_term(result, quantity="smf_eta", **kwargs)


def _none_message(quantity):
    '''The helpful error when the chosen scalar is None for the record.'''
    if quantity == "eta_turb":
        return ("quantity='eta_turb' is None: the record is not an uplink. "
                "eta_turb is the uplink reciprocity overlap; a downlink or a "
                "terrestrial record does not set it. Use 'collected_power' for a "
                "downlink aperture, or 'smf_eta' for a fibre receiver.")
    if quantity == "smf_eta":
        return ("quantity='smf_eta' is None: the receive terminal has no SMF "
                "detector, so the run computed no fibre coupling. Give the rx "
                "terminal detector=SMF(...), or use 'collected_power'.")
    return (f"quantity={quantity!r} is None for the record. The direction does "
            "not set it.")


def _full_vacuum_loss_db(result):
    '''
    The full no-turbulence loss launch to collected power, from a WaveResult.

    This is the exact launch-to-collected power ratio (stage 0 to stage 3), NOT
    the sum tx_truncation_db + geometric_loss_db, which skips the propagation
    grid-tail loss (the scaled and Fresnel routes lose a few percent of the power
    at the receive plane). See olb.waveoptics.run.
    '''
    from ..waveoptics.field import Power
    p_launch = Power(result.stages[0][1])       # "launch"
    p_collected = Power(result.stages[3][1])    # "after rx clip"
    return -10.0 * np.log10(float(p_collected) / float(p_launch))


def waveoptics_vacuum_term(result, *, include_smf=None, name=None,
                           category="geometric", beam_type=BEAM_GAUSSIAN,
                           note=None, meta_extra=None):
    '''
    Build the DETERMINISTIC vacuum-optics Term from a no-turbulence WaveResult.

    This is the fidelity-2 companion of the analytic geometric Term. It is the
    full no-turbulence loss from launch to detector: the launch truncation, the
    free-space diffraction spread, the receive-aperture capture, and (for an SMF
    receiver) the vacuum fibre coupling. The wave-optics field solve computes all
    of these directly, so at fidelity 2 this ONE Term replaces the analytic
    geometric_loss_term, tx_gaussian_efficiency_term, and the static coupling.

    It carries NO fade (one propagation, no turbulence), so it has no sampler and
    no quantile. It is NOT mean_only: a deterministic Term does not lock the
    budget to fidelity 0 (that lock is for a fluctuating quantity modelled by its
    mean alone). Its quantile is its mean, constant across availability.

    Parameters:
        result : WaveResult
            The record of one no-turbulence propagation
            (olb.waveoptics.run.propagate_scenario), or Fidelity2Bundle.vacuum.
        include_smf : bool, optional
            Add the vacuum SMF coupling loss (result.smf_coupling_db). None (the
            default) adds it when the record has it (an SMF receiver).
        name, category : str, optional
            Override the default Term name and category ("geometric").
        beam_type : str
            The launch beam type for the assumptions.
        note : str, optional
            The Term note.
        meta_extra : dict, optional
            Extra meta merged into the Term meta.

    Returns:
        Term
            A deterministic geometric-category Term (mean_db only).
    '''
    full_db = _full_vacuum_loss_db(result)
    smf_db = result.smf_coupling_db
    use_smf = (smf_db is not None) if include_smf is None else bool(include_smf)
    if use_smf and smf_db is None:
        raise ValueError(
            "include_smf=True, but the vacuum WaveResult has no smf_coupling_db "
            "(the receive terminal is not an SMF). Use include_smf=False or give "
            "the rx terminal an SMF detector."
        )
    mean_db = full_db + (float(smf_db) if use_smf else 0.0)

    if name is None:
        name = ("vacuum optics (geometry+truncation+coupling)" if use_smf
                else "vacuum optics (geometry+truncation)")
    if note is None:
        note = (f"no-turbulence field loss, {result.propagator}, "
                f"full {full_db:.3f} dB"
                + (f" + SMF {float(smf_db):.3f} dB" if use_smf else ""))
    # A deterministic vacuum loss carries no turbulence, so it declares no
    # regime and no spectrum. That keeps it neutral to the one-spectrum
    # Budget.check() rule (it sits beside the turbulence Term).
    assumptions = Assumptions(
        beam_type=beam_type,
        turbulence_regime=REGIME_NA,
        spectrum=SPECTRUM_NA,
        validity="Fidelity-2 vacuum optics: the full NO-TURBULENCE loss from "
                 "launch to detector, from one field propagation "
                 "(olb.waveoptics.run.propagate_scenario). It holds the launch "
                 "truncation, the free-space diffraction spread, the receive "
                 "aperture capture, and the vacuum fibre coupling. It replaces "
                 "the analytic geometric, launch-truncation, and static-coupling "
                 "Terms at fidelity 2. It carries no fade.",
    )
    meta = {
        "model": "waveoptics-vacuum",
        "propagator": result.propagator,
        "full_db": full_db,
        "smf_coupling_db": None if smf_db is None else float(smf_db),
        "include_smf": use_smf,
    }
    if meta_extra:
        meta.update(meta_extra)
    return Term(name=name, category=category, mean_db=mean_db, note=note,
                meta=meta, assumptions=assumptions)


def run_fidelity2(scenario, geometry, *, n_trials=200, preset="standard",
                  seed=None, threader=None, hs=None, cn2_profile=None,
                  L0_m=np.inf, subharmonics=True):
    '''
    Run BOTH wave-optics propagations a fidelity-2 budget needs, ONE time each.

    A fidelity-2 budget shows two Terms, from two runs:
      - the TURBULENT split-step Monte Carlo (the fade), sized by turbulent_grid;
      - the VACUUM no-turbulence propagation (the full geometric loss).

    The grids differ by direction, and getting them right is the whole point:
      - SPACE: the turbulent runner propagates only the ~20 km atmosphere slab
        with a plane-wave input on a flat grid; it holds NONE of the full-path
        (600 km) diffraction, and its outputs are vacuum-limit-1.0 penalties. So
        the vacuum run uses its OWN co-moving grid over the full slant range, and
        the two Terms add cleanly.
      - TERRESTRIAL: the turbulent runner propagates the real launch beam over
        the full path, so its collected power holds the geometric spread. The
        vacuum run therefore uses the SAME flat grid, so the turbulence penalty
        (turbulent / vacuum on that grid) is exact.

    The budget never runs a simulation; it consumes the Fidelity2Bundle this
    returns. Give the same seed to repeat the atmosphere.

    Parameters:
        scenario, geometry : the link case and geometry (one range only).
        n_trials, preset, seed, threader, hs, cn2_profile, L0_m, subharmonics :
            passed to the turbulent run (see run_waveoptics).

    Returns:
        Fidelity2Bundle
            The vacuum WaveResult and the turbulent TurbWaveResult.
    '''
    from ..waveoptics.turbulence import (propagate_turbulent_scenario,
                                          turbulent_grid)
    from ..waveoptics.run import propagate_scenario
    grid, plan, _ = turbulent_grid(scenario, geometry, preset=preset, hs=hs,
                                   cn2_profile=cn2_profile, L0_m=L0_m)
    turbulent = propagate_turbulent_scenario(
        scenario, geometry, n_trials=n_trials, seed=seed, preset=preset,
        grid=grid, plan=plan, hs=hs, cn2_profile=cn2_profile, L0_m=L0_m,
        subharmonics=subharmonics, threader=threader)
    if hasattr(scenario, "direction"):
        # SPACE: the vacuum run needs its own co-moving grid over the full path.
        vacuum = propagate_scenario(scenario, geometry)
    else:
        # TERRESTRIAL: the vacuum run shares the flat turbulent grid, so the
        # turbulence penalty (turbulent / vacuum) is exact.
        vacuum = propagate_scenario(scenario, geometry, grid=grid)
    return Fidelity2Bundle(vacuum=vacuum, turbulent=turbulent)


if __name__ == '__main__':
    from ..waveoptics.turbulence.run import TurbTrial, TurbWaveResult

    # --- Part A: the reducer, with synthetic trials (no simulation) -----------
    # A hand-made record proves the reduction bit-for-bit and the guards, fast
    # and deterministically, without running the split-step layer.
    etas = np.linspace(0.3, 0.7, 200)
    trials = [TurbTrial(collected_power=float(e), smf_eta=float(e),
                        eta_turb=None, seed_key=(0, i), wall_time_s=0.0)
              for i, e in enumerate(etas)]
    rec = TurbWaveResult(trials=trials, grid=None, plan=None, report=None,
                         preset="rapid", seed_entropy=1234)

    term = waveoptics_turbulence_term(rec, quantity="smf_eta", sigma2_I=0.1)
    assert term.name == "receive coupling (SMF)" and term.category == "coupling"
    assert term.meta["model"] == "waveoptics" and term.meta["n_trials"] == 200
    assert not term.mean_only and term.stochastic and term.quantile is not None
    # The coupling wrapper is exactly quantity="smf_eta".
    wrap = waveoptics_smf_coupling_term(rec, sigma2_I=0.1)
    assert wrap.category == "coupling" and wrap.mean_db == term.mean_db
    # The reduction is EXACTLY the empirical loss statistics.
    loss = -10.0 * np.log10(etas)
    assert term.mean_db == float(loss.mean())
    assert term.quantile_db(0.9) == float(np.quantile(loss, 0.9))
    assert term.quantile_db(0.9) > term.mean_db
    # The sampler bootstraps the loss samples (draws are exact loss values).
    draws = term.sample_db(5000, np.random.default_rng(0))
    assert draws.shape == (5000,) and np.all(np.isfinite(draws))
    assert set(np.unique(draws)).issubset(set(loss.tolist()))
    # Two builds from one record give identical faces (deterministic reduction).
    term2 = waveoptics_turbulence_term(rec, quantity="smf_eta", sigma2_I=0.1)
    assert term2.mean_db == term.mean_db
    assert term2.quantile_db(0.75) == term.quantile_db(0.75)
    # The weak index gives a weak regime; the snapshot-only flag always fires.
    assert term.meta["weak_fluctuation_valid"]
    assert any("SNAPSHOT-ONLY" in v for v in term.assumptions.violations)
    # A strong index flags the regime but stays valid (not an accuracy break).
    strong = waveoptics_turbulence_term(rec, quantity="smf_eta", sigma2_I=2.0)
    assert strong.assumptions.turbulence_regime == REGIME_STRONG
    assert any("split-step field solver stays valid" in v
               for v in strong.assumptions.violations)

    # Spectrum follows the outer scale.
    assert term.assumptions.spectrum == SPECTRUM_KOLMOGOROV
    vk = waveoptics_turbulence_term(rec, quantity="smf_eta", L0_m=20.0)
    assert vk.assumptions.spectrum == SPECTRUM_VON_KARMAN and vk.meta["L0_m"] == 20.0

    # Guards: a None scalar, an unknown quantity, and an empty record all raise.
    try:
        waveoptics_turbulence_term(rec, quantity="eta_turb")   # None on a downlink
        raise AssertionError("eta_turb=None must raise")
    except ValueError as e:
        assert "not an uplink" in str(e)
    try:
        waveoptics_turbulence_term(rec, quantity="mystery")
        raise AssertionError("an unknown quantity must raise")
    except ValueError:
        pass
    empty = TurbWaveResult(trials=[], grid=None, plan=None, report=None,
                           preset="rapid", seed_entropy=0)
    try:
        waveoptics_turbulence_term(empty, quantity="smf_eta")
        raise AssertionError("an empty record must raise")
    except ValueError:
        pass
    # The tail adequacy is a PROPERTY OF THE SAMPLER (results.EmpiricalSampler),
    # not of any caller. 200 trials populate the 0.9 tail (20 samples), but the
    # 0.999 tail (0.2 samples) is under-sampled and the quantile warns.
    assert not term.sampler.undersampled(0.9)
    assert term.sampler.undersampled(0.999)
    assert term.sampler.n == 200
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        term.quantile_db(0.999)
    assert any("UNDER-SAMPLED" in str(x.message) for x in w)

    # The loss_db path takes a per-trial COMPOSITE loss directly (for a penalty
    # built from more than one scalar). It reduces to the same three faces.
    comp_loss = -10.0 * np.log10(etas) - 10.0 * np.log10(etas)   # e.g. capture x coupling
    comp = waveoptics_turbulence_term(rec, loss_db=comp_loss, quantity="smf_eta")
    assert comp.name == "turbulence (wave optics)" and comp.category == "turbulence"
    assert comp.mean_db == float(comp_loss.mean())
    assert comp.meta["quantity"] == "smf_eta"      # quantity only labels here
    # At least one of quantity / loss_db is required.
    try:
        waveoptics_turbulence_term(rec)
        raise AssertionError("neither quantity nor loss_db must raise")
    except ValueError as e:
        assert "quantity" in str(e)
    print("Part A (reducer + guards): passed")

    # --- Part B: one real terrestrial run (skips if aotools is absent) --------
    from ..scenario import TerrestrialScenario, TerrestrialChannel
    from ..geometry import HorizontalPath
    from ..terminal import Terminal, Transmitter, SMF

    scn = TerrestrialScenario(
        near=Terminal(aperture_m=0.2, wavelength_m=1550e-9,
                      transmitter=Transmitter(waist_m=0.02, power_dbm=30)),
        far=Terminal(aperture_m=0.2, wavelength_m=1550e-9,
                     detector=SMF(sensitivity_dbm=-40)),
        channel=TerrestrialChannel(path_length_m=3e3, attenuation_db_per_km=0.5,
                                   cn2=5e-15))
    geom = HorizontalPath(3e3)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            bundle = run_fidelity2(scn, geom, preset="rapid", n_trials=16, seed=7)
    except ImportError as e:
        print(f"aotools not installed ({e.__class__.__name__}); skipping the "
              "real-run part of the self-check.")
        print("self-check passed (Part A)")
        raise SystemExit

    # The deterministic vacuum-optics Term: the full no-turbulence loss.
    vac = waveoptics_vacuum_term(bundle.vacuum)
    assert vac.category == "geometric" and vac.meta["model"] == "waveoptics-vacuum"
    assert not vac.stochastic and vac.quantile is None and not vac.mean_only
    assert vac.meta["smf_coupling_db"] is not None      # an SMF receiver
    # A deterministic Term's quantile is its mean (constant across availability),
    # so it does NOT lock the budget to fidelity 0.
    assert vac.quantile_db(0.99) == vac.mean_db

    # The terrestrial SMF turbulence penalty (composite), vacuum on the SAME grid.
    tr = bundle.turbulent.trials
    coll = np.array([t.collected_power for t in tr])
    eta = np.array([t.smf_eta for t in tr])
    from ..waveoptics.field import Power
    vac_coll = float(Power(bundle.vacuum.stages[3][1])
                     / Power(bundle.vacuum.stages[1][1]))   # collected / after tx clip
    vac_smf_eta = 10.0 ** (-bundle.vacuum.smf_coupling_db / 10.0)
    penalty_loss = -10.0 * np.log10(coll / vac_coll) - 10.0 * np.log10(eta / vac_smf_eta)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pen = waveoptics_turbulence_term(
            bundle.turbulent, loss_db=penalty_loss, quantity="smf_eta",
            name="turbulence (wave optics)")
    assert pen.stochastic and not pen.mean_only and pen.meta["n_trials"] == 16

    # RECONSTRUCTION: vacuum Term + penalty == the direct launch->fibre loss of
    # each trial, tx_trunc - 10log10(collected_i * smf_i). The vac_coll and
    # vac_smf baselines cancel exactly.
    tx_trunc = bundle.vacuum.tx_truncation_db
    direct = tx_trunc - 10.0 * np.log10(coll * eta)
    reconstructed = vac.mean_db + penalty_loss
    assert np.allclose(reconstructed, direct, atol=1e-9), \
        (reconstructed[:3], direct[:3])
    print(f"terrestrial SMF fidelity 2 (3 km, rapid, 16 trials): "
          f"vacuum {vac.mean_db:.2f} dB + turbulence {pen.mean_db:.2f} dB")
    print("self-check passed")
