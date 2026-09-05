"""The fidelity-2 field against FAST and the analytic model: the space-downlink
SMF coupling-loss gap (backlog 2-W1, 2-AO).

THE QUESTION. Does the fidelity-2 wave-optics field read LESS single-mode-fibre
(SMF) coupling loss than the fidelity-1 FAST model on the space downlink, and by
how much? Backlog 2-W1 claims the field reads 0.7 to 2.9 dB LESS loss than FAST.
This study measures that gap, elevation by elevation, like for like.

THE THREE MODELS, per elevation. All three are UNCORRECTED (no adaptive optics,
no tip-tilt). This is the only fair comparison, because the fidelity-2 layer
applies NO correction (backlog 2-AO).
  FAST      the fidelity-1 model (olb.models.fast.smf_fast_term). The ground
            terminal has an EMPTY compensation stack, so FAST runs AO_MODE=NOAO
            automatically. term.mean_db is the mean of the per-sample coupling
            loss in dB.
  field     the fidelity-2 split-step field, run through a Campaign process pool
            (olb.waveoptics.turbulence.campaign.Campaign, workers processes; it
            saturates the machine, unlike a serial/threaded run_fidelity2).
            Each trial gives collected_power and smf_eta. The composite loss is
            -10 log10(collected_power * smf_eta). For a space link the
            collected_power is vacuum-normalised to 1.0 and smf_eta holds the
            static mode-match floor, so this composite matches FAST's floor plus
            turbulence normalisation.
  analytic  the fidelity-0 model (olb.models.coupling.downlink
            .downlink_coupling_term, smf_fidelity="mean"). It is deterministic.
            The empty compensation stack makes it the uncorrected
            Dikmelik-Davidson coupling curve.

THE HEADLINE. The gap = FAST mean_db - field mean_db, per elevation. A positive
gap means the field reads LESS loss than FAST, which is the 2-W1 claim.

PARITY (or the comparison is a confound).
  - Same Cn2. All three integrate the site Hufnagel-Valley profile by default
    (cn2=None). The study passes no custom profile to one model alone.
  - Same outer scale. The study sets a FINITE von Karman outer scale L0 = 25 m
    on BOTH FAST (fast_params L0) and the field (run_fidelity2 L0_m). A finite L0
    is required: with L0 = inf the largest tilt cells are capped by the GRID SIDE,
    so the SMF loss depends on the grid size (owner decision 2026-09-05; olb
    outer-scale study, backlog 2-P5). The analytic term is L0-agnostic (Kolmogorov
    r0), so it is an approximate reference here, not a full parity match.
  - mean-of-dB on every side. FAST and the field both take the mean of the
    per-sample or per-trial dB. The analytic is one dB value. The study does NOT
    mix in a -10 log10(mean eta).
  - defocus = 0. A 500 km downlink is far field, so the received curvature is
    about zero and the SMF defocus_m stays 0.0. This differs from the terrestrial
    near-field case, where the received curvature drives a real defocus (backlog
    0-P11, 2-W2).

THE FAST NOAO GRID GUARD. FAST's auto grid can undersample the low-order tilt in
NOAO mode and understate the loss by several dB (see olb.models.fast, the SUBHARM
and NPXLS notes). So the study runs an NPXLS convergence check at one mid
elevation BEFORE the sweep. It picks the smallest NPXLS whose mean_db is within a
tolerance of the largest NPXLS, and it uses that pinned NPXLS for every
elevation. SUBHARM stays True throughout.

THE ELEVATION LOOP is required. FAST and run_fidelity2 take a scalar elevation
only, so the study loops one line of sight at a time (backlog I-1).

VALIDATION ONLY. The study reads the production layer. It changes no olb module.

CAVEAT. This study certifies the UNCORRECTED rung only (backlog 2-AO). Fidelity 2
applies no tip-tilt removal and no adaptive optics, so the fidelity-1 NOAO run is
the only like-for-like FAST comparison. This is NOT a reference for an
AO-corrected link.

Sources:
- O. J. D. Farley and others, Opt. Express 30(13), 23050 (2022),
  DOI 10.1364/OE.458659. The FAST method (fidelity 1).
- Y. Dikmelik and F. M. Davidson, Appl. Opt. 44(23), 4946 (2005),
  DOI 10.1364/AO.44.004946. The analytic uncorrected SMF coupling curve.
- C. Ruilier, Proc. SPIE 3350, 319 (1998), DOI 10.1117/12.317094. The 0.8145
  mode-match limit of an unobscured circular aperture (see
  olb.models.coupling._common and olb.waveoptics.smf).
- Schmidt, Numerical Simulation of Optical Wave Propagation with Examples in
  MATLAB (2010), DOI 10.1117/3.866274, Ch. 9. The split-step field solve
  (fidelity 2).
- Andrews and Phillips, Laser Beam Propagation through Random Media, 2nd ed.
  (2005), DOI 10.1117/3.626196. Ch. 12 the Hufnagel-Valley profile; Ch. 8 the
  Rytov variance.
- J. H. Shapiro, JOSA 61(4), 492 (1971), DOI 10.1364/JOSA.61.000492. The
  reciprocity used inside the field solve.

Run it from the repository root:

    python -m validation.waveoptics_vs_fast.waveoptics_vs_fast
    python -m validation.waveoptics_vs_fast.waveoptics_vs_fast --elevations 30 60
    python -m validation.waveoptics_vs_fast.waveoptics_vs_fast --analyse-only
"""

import argparse
import json
import os
import time
import warnings

import numpy as np

# numpy 2.4 removed the np.trapz alias, and FAST still calls it. Restore the name
# here, in the validation script only. Do not change production code.
if not hasattr(np, "trapz"):
    np.trapz = np.trapezoid

from olb.geometry import CircularOrbit
from olb.scenario import Channel, Site, SpaceScenario
from olb.terminal import SMF, Terminal, Transmitter
from olb.models.coupling.downlink import downlink_coupling_term
from olb.models.fast import smf_fast_term
from olb.waveoptics.turbulence.campaign import Campaign

HERE = os.path.dirname(os.path.abspath(__file__))
# The campaign store (gitignored). An env var overrides it.
CAMPAIGN_ROOT = os.environ.get("OLB_WVF_ROOT", os.path.join(HERE, "campaigns"))

LAM = 1550e-9
ALT_M = 500e3

# The scenario Cn2 and wind. They match tail_convergence.py, so the two studies
# use the SAME hero downlink.
CN2_GROUND = 1.7e-14
WIND_RMS = 21.0

# The default seed of the fidelity-2 run.
SEED = 20260905

# The outer scale of the screens, in m. A FINITE value is REQUIRED: with L0 = inf
# the largest tilt cells are capped by the grid side, so the SMF loss depends on
# the grid size (owner decision 2026-09-05; olb outer-scale study, backlog 2-P5,
# 0-W4). The SAME L0 goes to the field AND FAST for parity; the analytic term is
# L0-agnostic (Kolmogorov r0).
OUTER_SCALE_M = 25.0

# The default sweep. FAST and run_fidelity2 take a scalar elevation only, so the
# study loops (backlog I-1).
ELEVATIONS = (20.0, 30.0, 60.0, 90.0)

# The default trial and sample counts.
N_TRIALS = 1000        # fidelity-2 trials
N_SAMPLES = 1000       # FAST Monte Carlo draws

# The FAST NOAO grid guard.
NPXLS_SET = (128, 256, 512)
NPXLS_TOL = 0.15       # dB, the convergence tolerance of the pinned NPXLS
NPXLS_CAL_ELEVATION = 30.0

# The exceedance probabilities of the fade. pX is the loss that the link EXCEEDS
# X percent of the time, so it is the (1 - X) quantile of the loss.
EXCEEDANCE = (0.50, 0.05, 0.01)

# The bootstrap of every statistic: the resample count and the interval.
N_BOOT = 1000
BOOT_INTERVAL = 0.68
BOOT_SEED = 7

# The 2-W1 claim: the field reads 0.7 to 2.9 dB LESS loss than FAST.
GAP_CLAIM_LOW = 0.7
GAP_CLAIM_HIGH = 2.9


# ---------------------------------------------------------------------------
# The scenario
# ---------------------------------------------------------------------------

def scenario_and_geometry(elevation_deg):
    """Build the downlink scenario and the orbit of one elevation.

    The case is the hero downlink of tail_convergence.py: a 1550 nm
    space-to-ground link to a 700 mm ground telescope with a single-mode-fibre
    receiver, from a 100 mm space terminal at 500 km.

    The ground terminal has an EMPTY compensation stack, so every model runs
    uncorrected. The SMF defocus_m stays 0.0, because a 500 km downlink is far
    field and the received curvature is about zero.

    Args:
        elevation_deg: the elevation of the line of sight, in deg.

    Returns:
        The pair (SpaceScenario, CircularOrbit). The orbit holds a length-1
        elevation array, which is ONE line of sight.
    """
    site = Site(cn2_ground=CN2_GROUND, wind_rms_m_s=WIND_RMS)
    channel = Channel(site=site, altitude_m=ALT_M)
    ground = Terminal(aperture_m=0.7, wavelength_m=LAM,
                      pointing_jitter_rad=2e-6,
                      detector=SMF(sensitivity_dbm=-45.0))
    space = Terminal(aperture_m=0.10, wavelength_m=LAM,
                     pointing_jitter_rad=1e-6,
                     transmitter=Transmitter(waist_m=0.04, power_dbm=30.0))
    scn = SpaceScenario(ground=ground, space=space, direction="downlink",
                        channel=channel)
    geom = CircularOrbit(altitude_m=ALT_M, elevation_deg=[float(elevation_deg)])
    return scn, geom


def _scalar_orbit(elevation_deg):
    """Give the orbit of ONE scalar elevation. FAST refuses an array."""
    return CircularOrbit(altitude_m=ALT_M, elevation_deg=float(elevation_deg))


def _confirm_no_ao(scn):
    """Assert the ground terminal has an empty compensation stack.

    The uncorrected comparison needs no adaptive optics and no tip-tilt. This
    check makes the parity mechanical.
    """
    stack = scn.rx_terminal.compensation
    assert not stack, ("the ground terminal must have an empty compensation "
                       "stack for the uncorrected comparison; found %r" % (stack,))
    assert scn.rx_terminal.detector.defocus_m == 0.0, \
        "the SMF defocus_m must be 0.0 for the far-field downlink"


# ---------------------------------------------------------------------------
# The statistics. Copied from validation/tail_convergence/tail_convergence.py.
# ---------------------------------------------------------------------------

def _quantiles(loss_db, probs=EXCEEDANCE):
    """Give the loss that the link EXCEEDS each probability.

    A fade is a LARGE loss, so the loss exceeded a fraction q of the time is the
    (1 - q) quantile of the loss sample. Source: the same helper in
    validation/tail_convergence/tail_convergence.py.

    Args:
        loss_db: the per-trial loss, in dB.
        probs:   the exceedance probabilities.

    Returns:
        A float array, one value for each probability.
    """
    x = np.asarray(loss_db, dtype=float)
    return np.quantile(x, [1.0 - q for q in probs])


def _bootstrap(x, fn, n_boot=N_BOOT, seed=BOOT_SEED):
    """Give the bootstrap half-width of a statistic.

    The function resamples the values with replacement, it takes the statistic
    of each resample, and it gives the half-width of the central BOOT_INTERVAL
    band. The rng seed is FIXED, so a rerun gives the same interval. Source: the
    same helper in validation/tail_convergence/tail_convergence.py.

    Args:
        x:      the sample, a 1-D array.
        fn:     a callable fn(sample) -> a float or a float array.
        n_boot: the number of resamples.
        seed:   the rng seed.

    Returns:
        The half-width. It has the shape of fn(x).
    """
    x = np.asarray(x, dtype=float)
    rng = np.random.default_rng(seed)
    draws = np.array([np.atleast_1d(fn(x[rng.integers(0, x.size, x.size)]))
                      for _ in range(n_boot)])
    lo = np.quantile(draws, 0.5 - BOOT_INTERVAL / 2.0, axis=0)
    hi = np.quantile(draws, 0.5 + BOOT_INTERVAL / 2.0, axis=0)
    half = np.atleast_1d((hi - lo) / 2.0)
    return half if half.size > 1 else float(half[0])


def _by_p(values):
    """Key an exceedance array by its pX label."""
    return {f"p{int(p * 100)}": float(v)
            for p, v in zip(EXCEEDANCE, np.atleast_1d(values))}


# ---------------------------------------------------------------------------
# The three models
# ---------------------------------------------------------------------------

def _fast_term(scn, elevation_deg, n_samples, npxls, l0):
    """Build the fidelity-1 FAST SMF-coupling Term of one elevation.

    FAST runs AO_MODE=NOAO, because the ground terminal has no compensation
    stack. SUBHARM stays True (the FAST default), so the low-order tilt is
    captured. The outer scale L0 matches the field, for parity. Source:
    olb.models.fast.smf_fast_term.

    Args:
        scn:           the downlink SpaceScenario.
        elevation_deg: the scalar elevation, in deg.
        n_samples:     the FAST Monte Carlo draws.
        npxls:         the pinned FAST grid size (NPXLS).
        l0:            the outer scale, in m (matches the field).

    Returns:
        The FAST Term.
    """
    geom = _scalar_orbit(elevation_deg)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return smf_fast_term(scn, geom, n_samples=int(n_samples),
                             fast_params={"NPXLS": int(npxls), "L0": float(l0)})


def _fast_row(scn, elevation_deg, n_samples, npxls, l0):
    """Measure the fidelity-1 FAST model at one elevation.

    The FAST Term holds an empirical loss distribution. The mean and the p5/p1
    exceedance quantiles come from the Term directly. The mean bootstrap bar
    comes from a resample of the Term's empirical distribution, which is an
    estimate of the sampling uncertainty of the mean.

    Args:
        scn:           the downlink SpaceScenario.
        elevation_deg: the scalar elevation, in deg.
        n_samples:     the FAST Monte Carlo draws.
        npxls:         the pinned FAST grid size.
        l0:            the outer scale, in m (matches the field).

    Returns:
        A dict of the FAST measurement.
    """
    term = _fast_term(scn, elevation_deg, n_samples, npxls, l0)
    mean_db = float(term.mean_db)
    # The exceedance quantiles from the Term (the true empirical quantiles of the
    # run). quantile_db(p) is the p quantile, so pX (exceeded X percent) is the
    # (1 - X) quantile.
    quant = {f"p{int(p * 100)}": float(term.quantile_db(1.0 - p))
             for p in EXCEEDANCE}
    # The mean bootstrap. term.sample_db draws from the empirical loss
    # distribution, so this bar estimates the sampling uncertainty of the mean.
    rng = np.random.default_rng(BOOT_SEED)
    draws = np.asarray(term.sample_db(int(n_samples), rng), dtype=float).ravel()
    mean_half = _bootstrap(draws, lambda s: float(s.mean()))
    return {
        "elevation_deg": float(elevation_deg),
        "mean_db": mean_db,
        "mean_db_half": float(mean_half),
        "quantiles_db": quant,
        "npxls": int(npxls),
        "n_samples": int(n_samples),
        "ao_mode": term.meta.get("ao_mode"),
        "zmax": term.meta.get("zmax"),
        "floor_db": float(term.meta.get("floor_db")),
        "r0_los_m": float(term.meta.get("r0_los_m")),
        "amplitude_sigma2_I": float(term.meta.get("amplitude_sigma2_I")),
    }


def _field_row(scn, geom, elevation_deg, n_trials, seed, preset, l0, workers,
               block_size, root, mode="process"):
    """Measure the fidelity-2 field at one elevation.

    Three execution modes for the SAME physics, so the study can compare the
    parallel routes:
      process (default) -- a `Campaign` runs the trials across a WARM
        ProcessPoolExecutor (workers processes, one block per process). No GIL,
        so it SATURATES the machine, and it is resumable.
      thread            -- bare `run_fidelity2` with a `Threader(workers)`. The
        threads share memory; the FFT releases the GIL, but the Python work
        between FFTs does not, so a thread run is GIL-capped (about 0.35
        efficiency at 16) and does not saturate.
      serial            -- `run_fidelity2` with no pool, one trial at a time.
    The numbers agree to Monte-Carlo noise across the three; only the wall time
    and the CPU use differ. Sources:
    olb.waveoptics.turbulence.campaign.Campaign (process);
    olb.models.waveoptics.run_fidelity2 with olb.waveoptics.Threader (thread).

    The composite per-trial loss is -10 log10(collected_power * smf_eta). This
    is the SAME composite that tail_convergence.py uses: the collected_power is
    vacuum-normalised to 1.0 for a space link, and smf_eta holds the static
    mode-match floor. So the composite matches FAST's floor plus turbulence
    normalisation.

    Args:
        scn:           the downlink SpaceScenario.
        geom:          the length-1 orbit of this elevation (one line of sight).
        elevation_deg: the elevation, in deg (for the record).
        n_trials:      the fidelity-2 trials.
        seed:          the run seed.
        preset:        the sampling preset ("standard").
        l0:            the outer scale, in m (matches FAST).
        workers:       the campaign process-pool size.
        block_size:    the trials in one block file.
        root:          the campaign root directory of this elevation (process
                       mode only).
        mode:          "process" (Campaign pool), "thread" (Threader) or
                       "serial".

    Returns:
        The pair (a dict of the field measurement, the per-trial loss array).
    """
    # The outer scale L0 matches FAST, for parity (a finite L0 makes the SMF loss
    # grid-independent; backlog 2-P5).
    if mode == "process":
        # The Campaign takes the length-1 array orbit directly (it does not hit
        # the run_fidelity2 recap).
        # The stored campaigns are DOUBLE precision (2026-09-05).
        camp = Campaign(scn, geom, root, seed=int(seed), preset=preset,
                        block_size=int(block_size), L0_m=float(l0),
                        precision="double")
        camp.run(int(n_trials), workers=int(workers), progress=True)
        trials = camp.load(int(n_trials), fields=False).trials
    elif mode in ("thread", "serial"):
        from olb.models.waveoptics import run_fidelity2
        from olb.waveoptics.threader import Threader
        threader = Threader(max_workers=int(workers)) if mode == "thread" else None
        # A SCALAR orbit avoids the run_fidelity2 recap bug on a length-1 array;
        # progress=False skips the recap entirely.
        geom_scalar = _scalar_orbit(elevation_deg)
        bundle = run_fidelity2(scn, geom_scalar, n_trials=int(n_trials),
                               preset=preset, seed=int(seed), L0_m=float(l0),
                               threader=threader, progress=False)
        trials = bundle.turbulent.trials
    else:
        raise ValueError(f"unknown field mode {mode!r} "
                         "(use process, thread or serial)")
    power = np.array([t.collected_power for t in trials], dtype=float)
    eta = np.array([t.smf_eta for t in trials], dtype=float)
    loss = -10.0 * np.log10(power * eta)
    quant = _quantiles(loss)
    row = {
        "elevation_deg": float(elevation_deg),
        "mean_db": float(loss.mean()),
        "mean_db_half": _bootstrap(loss, lambda s: float(s.mean())),
        "quantiles_db": _by_p(quant),
        "quantiles_half_db": _by_p(_bootstrap(loss, _quantiles)),
        "n_trials": int(loss.size),
        "seed": int(seed),
        "preset": preset,
        "field_mode": mode,
    }
    return row, loss


def _analytic_loss(scn, geom):
    """Give the fidelity-0 analytic SMF coupling loss of one elevation.

    downlink_coupling_term with smf_fidelity="mean" is deterministic. The empty
    compensation stack makes it the uncorrected Dikmelik-Davidson coupling curve.
    Source: olb.models.coupling.downlink.downlink_coupling_term.

    Args:
        scn:  the downlink SpaceScenario.
        geom: the length-1 orbit of this elevation.

    Returns:
        The coupling loss in dB, a float.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        term = downlink_coupling_term(scn, geom, smf_fidelity="mean")
    return float(np.ravel(np.asarray(term.mean_db, dtype=float))[0])


# ---------------------------------------------------------------------------
# The FAST NOAO NPXLS convergence guard
# ---------------------------------------------------------------------------

def npxls_convergence(scn, elevation_deg, npxls_set, n_samples, tol, l0):
    """Pick the FAST NPXLS grid at one mid elevation.

    FAST's auto grid can undersample the low-order tilt in NOAO mode. So the
    study runs FAST at each NPXLS of the set with SUBHARM True. It picks the
    smallest NPXLS whose mean_db is within `tol` of the largest NPXLS.

    Args:
        scn:           the downlink SpaceScenario.
        elevation_deg: the calibration elevation, in deg.
        npxls_set:     the NPXLS values to test.
        n_samples:     the FAST Monte Carlo draws.
        tol:           the convergence tolerance, in dB.
        l0:            the outer scale, in m (matches the field).

    Returns:
        The pair (the pinned NPXLS, the table rows).
    """
    npxls_set = sorted(int(n) for n in npxls_set)
    rows = []
    for npxls in npxls_set:
        term = _fast_term(scn, elevation_deg, n_samples, npxls, l0)
        rows.append({"npxls": npxls, "mean_db": float(term.mean_db),
                     "r0_los_m": float(term.meta.get("r0_los_m"))})
    reference = rows[-1]["mean_db"]        # the largest NPXLS
    pinned = npxls_set[-1]
    for r in rows:
        if abs(r["mean_db"] - reference) <= tol:
            pinned = r["npxls"]
            break
    return pinned, rows


# ---------------------------------------------------------------------------
# The figures
# ---------------------------------------------------------------------------

def _plot(rows, survival, fig_dir):
    """Draw the figures. Give the list of the written paths."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib is not installed; the figures are skipped.")
        return []

    os.makedirs(fig_dir, exist_ok=True)
    written = []
    elev = [r["elevation_deg"] for r in rows]

    # ---- 1. the mean loss of the three models ----
    fig, ax = plt.subplots(figsize=(7.5, 5.0))
    ax.errorbar(elev, [r["fast_mean_db"] for r in rows],
                yerr=[r["fast_mean_half_db"] for r in rows],
                marker="o", capsize=3, label="FAST (fidelity 1, NOAO)")
    ax.errorbar(elev, [r["field_mean_db"] for r in rows],
                yerr=[r["field_mean_half_db"] for r in rows],
                marker="s", capsize=3, label="field (fidelity 2)")
    ax.plot(elev, [r["analytic_db"] for r in rows], "^--",
            label="analytic (fidelity 0, uncorrected)")
    ax.set_xlabel("elevation [deg]")
    ax.set_ylabel("SMF coupling loss [dB]")
    ax.set_title("The mean SMF coupling loss of the three models")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    path = os.path.join(fig_dir, "mean_loss_vs_elevation.png")
    fig.savefig(path, dpi=140)
    plt.close(fig)
    written.append(path)

    # ---- 2. the gap: FAST minus field ----
    fig, ax = plt.subplots(figsize=(7.5, 5.0))
    ax.errorbar(elev, [r["gap_db"] for r in rows],
                yerr=[r["gap_half_db"] for r in rows],
                marker="o", capsize=3, label="gap = FAST - field")
    ax.axhspan(GAP_CLAIM_LOW, GAP_CLAIM_HIGH, alpha=0.15, color="green",
               label=f"2-W1 claim {GAP_CLAIM_LOW} to {GAP_CLAIM_HIGH} dB")
    ax.axhline(0.0, color="k", lw=0.8)
    ax.set_xlabel("elevation [deg]")
    ax.set_ylabel("gap [dB] (positive = field reads less loss)")
    ax.set_title("The FAST minus field SMF coupling-loss gap (2-W1)")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    path = os.path.join(fig_dir, "gap_vs_elevation.png")
    fig.savefig(path, dpi=140)
    plt.close(fig)
    written.append(path)

    # ---- 3. the p5 and p1 fade of FAST and the field ----
    fig, ax = plt.subplots(figsize=(7.5, 5.0))
    for key, style in (("p5", "o-"), ("p1", "s--")):
        ax.plot(elev, [r["fast_quantiles_db"][key] for r in rows], style,
                label=f"FAST {key}")
        ax.plot(elev, [r["field_quantiles_db"][key] for r in rows], style,
                label=f"field {key}")
    ax.set_xlabel("elevation [deg]")
    ax.set_ylabel("SMF coupling loss [dB]")
    ax.set_title("The p5 and p1 fade of FAST and the field")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    path = os.path.join(fig_dir, "fade_vs_elevation.png")
    fig.savefig(path, dpi=140)
    plt.close(fig)
    written.append(path)

    # ---- 4. the survival functions of one elevation ----
    if survival:
        fig, ax = plt.subplots(figsize=(7.5, 5.0))
        for name, loss in survival.items():
            s = np.sort(np.asarray(loss, dtype=float))
            frac = 1.0 - np.arange(s.size) / s.size
            ax.semilogy(s, frac, label=name)
        ax.set_xlabel("SMF coupling loss [dB]")
        ax.set_ylabel("fraction of trials with a larger loss")
        ax.set_title("The field survival function (the mid elevation)")
        ax.grid(alpha=0.3, which="both")
        ax.legend(fontsize=8)
        fig.tight_layout()
        path = os.path.join(fig_dir, "field_survival.png")
        fig.savefig(path, dpi=140)
        plt.close(fig)
        written.append(path)
    return written


# ---------------------------------------------------------------------------
# The driver
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--elevations", nargs="+", type=float,
                    default=list(ELEVATIONS),
                    help="the elevations of the sweep [deg]")
    ap.add_argument("--n-trials", type=int, default=N_TRIALS,
                    help="the fidelity-2 trials of each elevation")
    ap.add_argument("--n-samples", type=int, default=N_SAMPLES,
                    help="the FAST Monte Carlo draws of each elevation")
    ap.add_argument("--npxls-set", nargs="+", type=int, default=list(NPXLS_SET),
                    help="the FAST NPXLS values of the convergence guard")
    ap.add_argument("--npxls-tol", type=float, default=NPXLS_TOL,
                    help="the NPXLS convergence tolerance [dB]")
    ap.add_argument("--npxls-cal-elevation", type=float,
                    default=NPXLS_CAL_ELEVATION,
                    help="the elevation of the NPXLS convergence guard [deg]")
    ap.add_argument("--seed", type=int, default=SEED,
                    help="the fidelity-2 run seed")
    ap.add_argument("--workers", type=int, default=16,
                    help="the campaign process-pool size for the field (16 on a "
                         "32-core box; processes saturate, threads do not)")
    ap.add_argument("--block-size", type=int, default=50,
                    help="the trials in one campaign block file")
    ap.add_argument("--field-mode", choices=("process", "thread", "serial"),
                    default="process",
                    help="how to run the field trials: process (Campaign pool, "
                         "saturates), thread (Threader, GIL-capped), or serial")
    ap.add_argument("--L0", type=float, default=OUTER_SCALE_M,
                    help="the outer scale [m], threaded to the field AND FAST "
                         "(default 25; 'inf' is the grid-dependent Kolmogorov "
                         "limit and is NOT recommended)")
    ap.add_argument("--analyse-only", action="store_true",
                    help="skip the runs and re-read the stored results")
    args = ap.parse_args()

    # Tag the outputs by outer scale and field mode, so an inf run and a 25 m
    # run, or a process run and a thread run, do not clobber each other.
    l0_lab = "inf" if not np.isfinite(args.L0) else f"{args.L0:g}m"
    tag = f"L0{l0_lab}_{args.field_mode}"
    log_path = os.path.join(HERE, f"waveoptics_vs_fast_{tag}.log")
    json_path = os.path.join(HERE, f"waveoptics_vs_fast_{tag}_results.json")
    # The log is written LINE BY LINE, so a killed run keeps its log.
    with open(log_path, "w", encoding="utf-8"):
        pass

    def say(text=""):
        print(text, flush=True)
        with open(log_path, "a", encoding="utf-8") as fh:
            fh.write(text + "\n")

    say("The fidelity-2 field against FAST: the space-downlink SMF gap "
        "(backlog 2-W1, 2-AO)")
    say(f"elevations    : {args.elevations} deg")
    say(f"trials/field  : {args.n_trials}      samples/FAST: {args.n_samples}      "
        f"seed: {args.seed}")
    say(f"field mode    : {args.field_mode}   workers: {args.workers}   "
        f"(process = Campaign pool, saturates; thread = Threader, GIL-capped)")
    say(f"scenario      : downlink, 1550 nm, 500 km, 700 mm ground SMF, "
        f"100 mm space terminal, Cn2_0={CN2_GROUND:.1e}")
    say("correction    : NONE on every model. Fidelity 2 applies no tip-tilt and "
        "no AO (backlog 2-AO), so the FAST NOAO run is the like-for-like pair.")
    say(f"outer scale   : L0 = {args.L0:g} m on the field AND FAST (a finite L0 "
        f"makes the SMF loss grid-independent; backlog 2-P5). The analytic term "
        f"is L0-agnostic.")
    say("parity        : same Cn2 (site HV), same outer scale on field and FAST, "
        "mean-of-dB on FAST and the field, SMF defocus_m=0 (far field).")
    say(f"mode          : {'ANALYSE ONLY' if args.analyse_only else 'RUN'}")
    if args.n_trials < 1000:
        say("WARNING: under 1000 trials, so the p1 tail is UNDER-SAMPLED. Read "
            "the mean and p5 rows only.")
    say()

    if args.analyse_only:
        if not os.path.exists(json_path):
            say(f"no stored results at {json_path}. Run without --analyse-only "
                "first.")
            return
        with open(json_path, "r", encoding="utf-8") as fh:
            stamp = json.load(fh)
        rows = stamp["rows"]
        npxls_rows = stamp.get("npxls_convergence", [])
        pinned = stamp.get("pinned_npxls")
        survival = {}
        _report(say, rows, npxls_rows, pinned, args)
        figs = _plot(rows, survival, os.path.join(HERE, "figures"))
        for f in figs:
            say(f"wrote {f}")
        say(f"read  {json_path}")
        say(f"wrote {log_path}")
        return

    # ---- the FAST NOAO NPXLS convergence guard ----
    say("THE FAST NOAO NPXLS CONVERGENCE GUARD")
    say(f"  calibration elevation {args.npxls_cal_elevation:.1f} deg, "
        f"NPXLS {args.npxls_set}, tolerance {args.npxls_tol} dB, SUBHARM True")
    scn_cal, _ = scenario_and_geometry(args.npxls_cal_elevation)
    _confirm_no_ao(scn_cal)
    pinned, npxls_rows = npxls_convergence(scn_cal, args.npxls_cal_elevation,
                                           args.npxls_set, args.n_samples,
                                           args.npxls_tol, args.L0)
    say(f"  {'NPXLS':>7s}{'mean_db':>10s}{'r0_los[cm]':>12s}")
    for r in npxls_rows:
        say(f"  {r['npxls']:7d}{r['mean_db']:10.3f}{r['r0_los_m'] * 100:12.2f}")
    say(f"  PINNED NPXLS = {pinned} (the smallest within {args.npxls_tol} dB of "
        f"NPXLS {npxls_rows[-1]['npxls']})")
    say()

    # ---- the elevation sweep ----
    rows = []
    survival = {}
    mid_index = len(args.elevations) // 2
    for i, el in enumerate(args.elevations):
        say(f"ELEVATION {el:.1f} deg")
        scn, geom = scenario_and_geometry(el)
        _confirm_no_ao(scn)

        t0 = time.perf_counter()
        fast = _fast_row(scn, el, args.n_samples, pinned, args.L0)
        t_fast = time.perf_counter() - t0
        say(f"  FAST    : mean {fast['mean_db']:.3f} +-{fast['mean_db_half']:.3f} "
            f"dB, p5 {fast['quantiles_db']['p5']:.3f}, "
            f"p1 {fast['quantiles_db']['p1']:.3f} dB, "
            f"r0_los {fast['r0_los_m'] * 100:.2f} cm, "
            f"ao {fast['ao_mode']}, floor {fast['floor_db']:.3f} dB "
            f"({t_fast:.1f} s)")

        t0 = time.perf_counter()
        el_root = os.path.join(CAMPAIGN_ROOT, l0_lab, f"el{el:g}")
        field, loss = _field_row(scn, geom, el, args.n_trials, args.seed,
                                 "standard", args.L0, args.workers,
                                 args.block_size, el_root, mode=args.field_mode)
        t_field = time.perf_counter() - t0
        say(f"  field   : mean {field['mean_db']:.3f} "
            f"+-{field['mean_db_half']:.3f} dB, "
            f"p5 {field['quantiles_db']['p5']:.3f}, "
            f"p1 {field['quantiles_db']['p1']:.3f} dB "
            f"({t_field:.1f} s)")

        analytic = _analytic_loss(scn, geom)
        say(f"  analytic: mean {analytic:.3f} dB (deterministic, uncorrected)")

        gap = fast["mean_db"] - field["mean_db"]
        gap_half = float(np.hypot(fast["mean_db_half"], field["mean_db_half"]))
        field_vs_analytic = field["mean_db"] - analytic
        say(f"  GAP     : FAST - field = {gap:+.3f} +-{gap_half:.3f} dB   "
            f"field - analytic = {field_vs_analytic:+.3f} dB")
        say()

        rows.append({
            "elevation_deg": float(el),
            "r0_los_m": fast["r0_los_m"],
            "fast_mean_db": fast["mean_db"],
            "fast_mean_half_db": fast["mean_db_half"],
            "fast_quantiles_db": fast["quantiles_db"],
            "fast_ao_mode": fast["ao_mode"],
            "fast_zmax": fast["zmax"],
            "fast_floor_db": fast["floor_db"],
            "fast_amplitude_sigma2_I": fast["amplitude_sigma2_I"],
            "field_mean_db": field["mean_db"],
            "field_mean_half_db": field["mean_db_half"],
            "field_quantiles_db": field["quantiles_db"],
            "field_quantiles_half_db": field["quantiles_half_db"],
            "field_n_trials": field["n_trials"],
            "analytic_db": analytic,
            "gap_db": float(gap),
            "gap_half_db": gap_half,
            "field_minus_analytic_db": float(field_vs_analytic),
            "fast_seconds": float(t_fast),
            "field_seconds": float(t_field),
        })
        if i == mid_index:
            survival["field (fidelity 2)"] = loss

    # ---- the tables and the verdict ----
    _report(say, rows, npxls_rows, pinned, args)

    stamp = {
        "study": "waveoptics_vs_fast",
        "backlog": ["2-W1", "2-AO"],
        "elevations_deg": list(args.elevations),
        "n_trials": args.n_trials, "n_samples": args.n_samples,
        "workers": args.workers, "block_size": args.block_size,
        "seed": args.seed, "wavelength_m": LAM, "altitude_m": ALT_M,
        "outer_scale_m": float(args.L0),
        "cn2_ground": CN2_GROUND, "wind_rms_m_s": WIND_RMS,
        "exceedance": list(EXCEEDANCE), "n_bootstrap": N_BOOT,
        "bootstrap_interval": BOOT_INTERVAL,
        "pinned_npxls": int(pinned),
        "npxls_tol_db": args.npxls_tol,
        "npxls_cal_elevation_deg": args.npxls_cal_elevation,
        "npxls_convergence": npxls_rows,
        "gap_claim_db": [GAP_CLAIM_LOW, GAP_CLAIM_HIGH],
        "rows": rows,
    }
    with open(json_path, "w", encoding="utf-8") as fh:
        json.dump(stamp, fh, indent=2)
    figs = _plot(rows, survival, os.path.join(HERE, "figures"))
    print(f"\nwrote {json_path}")
    for f in figs:
        print(f"wrote {f}")
    print(f"wrote {log_path}")


def _report(say, rows, npxls_rows, pinned, args):
    """Print the tables and the verdict from the measured rows."""
    if pinned is not None:
        say(f"PINNED FAST NPXLS = {pinned}")
        say()

    # ---- the per-elevation table ----
    say("THE PER-ELEVATION TABLE (SMF coupling loss, dB; pX = loss EXCEEDED X %)")
    hdr = (f"{'elev':>6s}{'r0[cm]':>8s}{'FASTmean':>18s}{'fieldmean':>18s}"
           f"{'analytic':>9s}{'GAP':>16s}{'f-an':>8s}"
           f"{'FASTp5':>8s}{'FASTp1':>8s}{'fldp5':>8s}{'fldp1':>8s}")
    say(hdr)
    say("-" * len(hdr))
    for r in rows:
        say(f"{r['elevation_deg']:6.1f}{r['r0_los_m'] * 100:8.2f}"
            f"{r['fast_mean_db']:11.3f}+-{r['fast_mean_half_db']:4.3f}"
            f"{r['field_mean_db']:11.3f}+-{r['field_mean_half_db']:4.3f}"
            f"{r['analytic_db']:9.3f}"
            f"{r['gap_db']:9.3f}+-{r['gap_half_db']:4.3f}"
            f"{r['field_minus_analytic_db']:8.3f}"
            f"{r['fast_quantiles_db']['p5']:8.2f}"
            f"{r['fast_quantiles_db']['p1']:8.2f}"
            f"{r['field_quantiles_db']['p5']:8.2f}"
            f"{r['field_quantiles_db']['p1']:8.2f}")
    say("  GAP = FASTmean - fieldmean (positive = field reads LESS loss, the "
        "2-W1 claim). f-an = fieldmean - analytic.")
    say("  FASTp5/p1 and fldp5/p1 are the loss EXCEEDED 5 and 1 percent of the "
        "time. The +- values are 68 % bootstrap half-widths.")
    say()

    # ---- the verdict ----
    say("THE VERDICT")
    gaps = np.array([r["gap_db"] for r in rows], dtype=float)
    bars = np.array([r["gap_half_db"] for r in rows], dtype=float)
    elevs = np.array([r["elevation_deg"] for r in rows], dtype=float)
    gmin, gmax = float(gaps.min()), float(gaps.max())
    say(f"  The gap runs from {gmin:+.2f} to {gmax:+.2f} dB across "
        f"{elevs.min():.0f} to {elevs.max():.0f} deg.")

    # Does the field read consistently LESS loss than FAST (a positive gap by
    # more than the bootstrap bar)?
    resolved_less = bool(np.all(gaps > bars))
    any_less = bool(np.all(gaps > 0.0))
    if resolved_less:
        say("  The field reads LESS loss than FAST at EVERY elevation, by more "
            "than the bootstrap bar. This CONFIRMS the 2-W1 direction.")
    elif any_less:
        say("  The field reads less loss than FAST at every elevation, but not "
            "always past the bootstrap bar. The 2-W1 direction holds, weakly.")
    else:
        say("  The field does NOT read consistently less loss than FAST. The "
            "2-W1 direction is NOT reproduced here.")

    # The trend with elevation.
    if gaps.size >= 2:
        trend = "rises" if gaps[-1] > gaps[0] else "falls"
        say(f"  The gap {trend} from {gaps[0]:+.2f} dB at {elevs[0]:.0f} deg to "
            f"{gaps[-1]:+.2f} dB at {elevs[-1]:.0f} deg (low elevation = more "
            "turbulence).")

    # Against the backlog claim band.
    in_band = int(np.sum((gaps >= GAP_CLAIM_LOW) & (gaps <= GAP_CLAIM_HIGH)))
    say(f"  {in_band} of {gaps.size} elevations sit inside the 2-W1 claim band "
        f"{GAP_CLAIM_LOW} to {GAP_CLAIM_HIGH} dB.")

    say("  DEFOCUS is about zero for the 500 km downlink (far field), so the "
        "residual gap is NOT a defocus effect. It is the Airy-versus-Gaussian "
        "spot shape of the field against the FAST mode overlap. That is an owner "
        "interpretation question. This study does NOT resolve it.")
    say("  CAVEAT: this certifies the UNCORRECTED rung only (backlog 2-AO). "
        "Fidelity 2 applies no correction, so this is NOT a reference for an "
        "AO-corrected link.")
    say()


if __name__ == '__main__':
    main()
