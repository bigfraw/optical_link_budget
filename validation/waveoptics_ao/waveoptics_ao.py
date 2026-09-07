"""The perfect-AO validation of the fidelity-2 wave-optics layer (backlog 2-AO).

WHAT IT MEASURES. The fidelity-2 layer now removes the first N Noll modes of the
wavefront over the receive aperture (`olb.waveoptics.compensation`). This study
measures that chain against the book, and it measures what the correction does
to the single-mode-fibre (SMF) fade of the hero downlink.

THE FIVE STUDIES.
  V0  the trust gate of the summed-screen source. It compares the summed screen
      phase at the aperture against the phase of the receive FIELD, through the
      wrapped-gradient slopes. A space link senses the summed screen phase, so
      this measurement is the validity envelope of that source.
  V2  the physics cross-check. It measures the residual phase variance over the
      aperture after the removal of J = 1, 3, 10, 21 and 35 Noll modes, and it
      compares that with Noll's Delta_J (D / r0)^(5/3).
  V1  the fade. It measures the SMF loss of the hero downlink with no
      correction, with TipTilt, with AO(10) and with AO(21). It compares the
      result with the tracked fidelity-1 FAST Term.
  V3  method (a) against method (c). It corrects ONE stored campaign two ways,
      from the stored screen phase and from the field slopes, and it compares
      the coupling efficiency trial for trial.
  V4  the terrestrial sanity check. It applies the slope correction post hoc to
      two stored 2-TC terrestrial campaigns.
  V5  method (a) against method (c) on the TERRESTRIAL link. It is V3 and V0 on
      a horizontal path. It makes its OWN terrestrial campaigns, because a
      2-TC campaign holds no screen phase.

THE HERO DOWNLINK. A 1550 nm space-to-ground link to a 700 mm ground telescope
with an SMF receiver, from a 100 mm space terminal at 500 km. It is the case of
`validation/tail_convergence/` and `validation/waveoptics_vs_fast/`.

THE OUTER SCALE IS 25 m ON EVERY RUN. The owner fixed it on 2026-09-05 (backlog
2-P5): L0 = inf is grid-dependent, and it is pessimistic on the fibre fade tail.
Noll's law assumes the Kolmogorov spectrum (L0 = inf), so V2 reports a LOW low-
order residual, and it measures that outer-scale effect with a control.

PERFECT AO, AND NOTHING MORE. The correction is an ideal modal fit of one
snapshot. There is no wavefront-sensor noise, no finite subaperture, no
aliasing, no servo lag, no branch point and no anisoplanatism. So every
corrected number here is the UPPER BOUND of the benefit of a corrector.

Sources:
- R. J. Noll, "Zernike polynomials and atmospheric turbulence," J. Opt. Soc.
  Am. 66(3), 207-211 (1976), DOI 10.1364/JOSA.66.000207. Table I (the mode
  order) and Table IV (the residual coefficients Delta_J).
- Schmidt, Numerical Simulation of Optical Wave Propagation (2010),
  DOI 10.1117/3.866274, Ch. 9. The split-step method of the trials.
- Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196. Ch. 12 the
  Hufnagel-Valley profile; Ch. 8 the Rytov variance; Ch. 6 the Fried parameter.
- O. J. D. Farley and others, Opt. Express 30(13), 23050 (2022),
  DOI 10.1364/OE.458659. The FAST method (fidelity 1).
- D. L. Fried, J. Opt. Soc. Am. 56, 1372 (1966), DOI 10.1364/JOSA.56.001372.
  The Fried parameter r0.

Run it from the repository root. The campaigns take the CUDA backend:

    python -m validation.waveoptics_ao.waveoptics_ao --study v0 --n-trials 200
    python -m validation.waveoptics_ao.waveoptics_ao --study v2 --n-trials 200
    python -m validation.waveoptics_ao.waveoptics_ao --study v1 --n-trials 1000
    python -m validation.waveoptics_ao.waveoptics_ao --study v3 --n-trials 200
    python -m validation.waveoptics_ao.waveoptics_ao --study v4 --n-trials 500
    python -m validation.waveoptics_ao.waveoptics_ao --study v5 --n-trials 500

Add `--analyse-only` to read what is stored and to compute no trial. Add
`--fft-backend numpy` to run on a host without a CUDA device. V4 always runs on
the host, because it reopens the stored 2-TC campaigns of `numpy`.
"""

import argparse
import json
import os
import time
import warnings

import numpy as np

from olb.geometry import CircularOrbit
from olb.scenario import Channel, Site, SpaceScenario
from olb.terminal import AO, SMF, Terminal, TipTilt, Transmitter
from olb.turbulence.ao import (NOLL_PISTON, NOLL_TIPTILT, _AO_EXP,
                               _NOLL_AO_CONST,
                               plane_wave_fried_parameter_profile)
from olb.turbulence.profiles import DEFAULT_HS, default_cn2_profile
from olb.waveoptics.compensation import (ApertureModes, circle, max_abs_step,
                                         wrapped_gradient)
from olb.waveoptics.turbulence import Campaign
# _rebuilt_fields is PRIVATE. A campaign gives no public Field of a stored
# trial (the open gap of CLAUDE.md, "ONE campaign gap is OPEN"), and
# examples/waveoptics/camera_tracking.py reads the same helper.
from olb.waveoptics.turbulence.run import SLOPE_STEP_WARN_RAD, _rebuilt_fields

HERE = os.path.dirname(os.path.abspath(__file__))
ENV_ROOT = "OLB_WAVEOPTICS_AO_ROOT"

# ---------------------------------------------------------------------------
# The case
# ---------------------------------------------------------------------------

LAM = 1550e-9
ALTITUDE_M = 500e3
GROUND_APERTURE_M = 0.70
PATCH_RADIUS_M = 0.35           # the radius of the 0.7 m aperture
SEED = 20260907
PRESET = "standard"

# The outer scale of the screens, in m. The owner fixed it on 2026-09-05
# (backlog 2-P5). L0 = inf is grid-dependent, and it is pessimistic on the
# fibre fade tail.
L0_M = 25.0

# The exceedance probabilities of every fade table. pX is the loss that the
# link EXCEEDS X percent of the time.
EXCEEDANCE = (0.5, 0.1, 0.05, 0.01)
N_BOOT = 400
BOOT_INTERVAL = 0.68
BOOT_SEED = 12345

# The compensation stacks under test. The name is the campaign directory.
STACKS = {
    "base": (),
    "tiptilt": (TipTilt(),),
    "ao10": (AO(n_modes=10),),
    "ao21": (AO(n_modes=21),),
}
STACK_ORDER = ("base", "tiptilt", "ao10", "ao21")

# The mode counts of the V0 fit and of the V2 ladder.
V0_MODES = 21
V2_MODES = (1, 3, 10, 21, 35)

# The FAST grid guard. FAST can undersample the low-order tilt, so the study
# picks the smallest NPXLS whose mean loss is within the tolerance of the
# largest one. Copied from validation/waveoptics_vs_fast/.
NPXLS_SET = (128, 256, 512)
NPXLS_TOL = 0.15


def noll_delta(j):
    """Give the Noll residual coefficient Delta_J of J removed modes.

    The value is the coefficient of sigma^2 = Delta_J (D / r0)^(5/3), the
    residual phase variance over a circular aperture after the removal of the
    first J Noll modes. J = 1 (piston only) and J = 3 (tip-tilt) take the
    printed Table IV values. A larger J takes the large-J asymptotic form
    0.2944 J^(-sqrt(3)/2), which agrees with Table IV to better than 2 percent
    at J = 10 and at J = 21.

    The three constants come from `olb.turbulence.ao`, so this study measures
    the values that the analytic ladder uses.

    Source: R. J. Noll, J. Opt. Soc. Am. 66(3), 207 (1976),
    DOI 10.1364/JOSA.66.000207, Table IV, printed p. 210.

    Args:
        j: the number of removed Noll modes.

    Returns:
        The coefficient Delta_J, as a float.
    """
    if int(j) == 1:
        return float(NOLL_PISTON)
    if int(j) == 3:
        return float(NOLL_TIPTILT)
    return float(_NOLL_AO_CONST * float(j) ** _AO_EXP)


def hero_scenario(elevation_deg, compensation=()):
    """Build the hero downlink and the orbit of one elevation.

    Args:
        elevation_deg: the elevation of the line of sight, in deg.
        compensation:  the compensation stack of the GROUND terminal.

    Returns:
        The pair (SpaceScenario, CircularOrbit).
    """
    site = Site(cn2_ground=1.7e-14, wind_rms_m_s=21.0)
    channel = Channel(site=site, altitude_m=ALTITUDE_M)
    ground = Terminal(aperture_m=GROUND_APERTURE_M, wavelength_m=LAM,
                      pointing_jitter_rad=2e-6,
                      detector=SMF(sensitivity_dbm=-45.0),
                      compensation=list(compensation))
    space = Terminal(aperture_m=0.10, wavelength_m=LAM,
                     pointing_jitter_rad=1e-6,
                     transmitter=Transmitter(waist_m=0.04, power_dbm=30.0))
    scn = SpaceScenario(ground=ground, space=space, direction="downlink",
                        channel=channel)
    geom = CircularOrbit(altitude_m=ALTITUDE_M,
                         elevation_deg=[float(elevation_deg)])
    return scn, geom


def campaigns_root():
    """Give the parent directory of every campaign of this study."""
    return os.environ.get(ENV_ROOT) or os.path.join(HERE, "campaigns")


def campaign_of(elevation_deg, stack_name, args):
    """Build (or reopen) the campaign of one elevation and one stack.

    Every campaign stores the field patch of the 0.7 m aperture AND the summed
    screen phase, so any study reads any campaign post hoc.

    Args:
        elevation_deg: the elevation, in deg.
        stack_name:    a key of STACKS.
        args:          the parsed command line.

    Returns:
        The pair (Campaign, the sizer warning texts).
    """
    scn, geom = hero_scenario(elevation_deg)
    root = os.path.join(campaigns_root(), f"el{elevation_deg:02.0f}",
                        stack_name)
    stack = STACKS[stack_name]
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        camp = Campaign(scn, geom, root, seed=SEED, preset=PRESET,
                        block_size=int(args.block_size),
                        patch_radius_m=PATCH_RADIUS_M, L0_m=L0_M,
                        precision="single", fft_backend=args.fft_backend,
                        compensation=(list(stack) or None),
                        store_screen_phase=True)
    return camp, sorted({str(w.message) for w in caught})


def ensure_trials(camp, n_trials, args, say):
    """Run the missing trials of a campaign, and report the wall time.

    Args:
        camp:     the Campaign.
        n_trials: the wanted trial count.
        args:     the parsed command line.
        say:      the log function.

    Returns:
        A dict of the timing, or None when the call computed no trial.
    """
    if args.analyse_only:
        say(f"  ANALYSE ONLY: {camp.n_stored} trials on disk in "
            f"{camp.root_dir}")
        return None
    missing = max(0, int(n_trials) - int(camp.n_stored))
    if missing == 0:
        say(f"  {camp.n_stored} trials already on disk. No trial computed.")
        return None
    t0 = time.perf_counter()
    # The CUDA backend runs ONE process for ONE device, so workers stays None.
    n_done = camp.run(int(n_trials), workers=args.workers, progress=False)
    wall = time.perf_counter() - t0
    say(f"  {n_done} trials on disk, {wall:.1f} s wall for {missing} new "
        f"trials ({wall / max(missing, 1):.3f} s/trial)")
    return {"wall_s": float(wall), "n_new": int(missing),
            "n_stored": int(n_done),
            "s_per_trial": float(wall / max(missing, 1))}


# ---------------------------------------------------------------------------
# The shared statistics
# ---------------------------------------------------------------------------

def _quantiles(loss_db, probs=EXCEEDANCE):
    """Give the loss that the link EXCEEDS each probability.

    A fade is a LARGE loss, so the loss exceeded a fraction q of the time is
    the (1 - q) quantile of the loss sample.
    """
    x = np.asarray(loss_db, dtype=float)
    return np.quantile(x, [1.0 - q for q in probs])


def _bootstrap(x, fn, n_boot=N_BOOT, seed=BOOT_SEED):
    """Give the bootstrap half-width of a statistic.

    The function resamples the trials with replacement, it takes the statistic
    of each resample, and it gives the half-width of the central BOOT_INTERVAL
    band. The seed is FIXED, so a rerun gives the same interval.
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


def _fade_row(loss_db):
    """Give the mean and the exceedance quantiles of a loss sample."""
    loss = np.asarray(loss_db, dtype=float)
    return {
        "n": int(loss.size),
        "mean_db": float(loss.mean()),
        "mean_db_half": _bootstrap(loss, lambda s: float(s.mean())),
        "quantiles_db": _by_p(_quantiles(loss)),
        "quantiles_half_db": _by_p(_bootstrap(loss, _quantiles)),
    }


def _gain_and_rms(reference, other):
    """Give the regression gain and the RMS difference of two coefficient sets.

    A single ratio is noisy, because one coefficient can be near zero. So the
    gain is the least-squares slope of `other` against `reference`, and the
    difference is an RMS over the whole set.

    Args:
        reference: the reference values, a 1-D array.
        other:     the compared values, the same shape.

    Returns:
        The triple (gain, the RMS difference, the RMS of the reference).
    """
    a = np.asarray(reference, dtype=float).ravel()
    b = np.asarray(other, dtype=float).ravel()
    rms = float(np.sqrt((a ** 2).mean()))
    diff = float(np.sqrt(((b - a) ** 2).mean()))
    gain = float((b @ a) / (a @ a))
    return gain, diff, rms


def aperture_modes(camp, n_modes):
    """Build the projector of the receive aperture of a campaign.

    The mask is the aperture of the clip terminal on the stored patch grid.
    This is the SAME mask that the runner and `recouple_compensated` build.

    Args:
        camp:    the Campaign.
        n_modes: the number of Noll modes of the fit.

    Returns:
        The ApertureModes object.
    """
    patch = camp.patch
    rx = camp.scenario.rx_terminal
    mask = circle(patch.n, rx.aperture_m / patch.pixel_m,
                  rx.obscuration_ratio)
    return ApertureModes(n_modes, patch.n, mask)


def _screen_map(result, row, patch):
    """Give the summed screen phase of one trial as a full-grid map."""
    phase = np.zeros(patch.n * patch.n, dtype=np.float64)
    phase[patch.indices] = result.screen_phase[row]
    return phase.reshape(patch.n, patch.n)


def _dir_bytes(path):
    """Give the total byte count of the files in one directory."""
    if not os.path.isdir(path):
        return 0
    return int(sum(os.path.getsize(os.path.join(path, f))
                   for f in os.listdir(path)
                   if os.path.isfile(os.path.join(path, f))))


# ---------------------------------------------------------------------------
# V0: the trust gate of the summed-screen source
# ---------------------------------------------------------------------------

def v0_row(camp, n_trials):
    """Compare the summed screen phase against the receive-field phase.

    THE QUESTION. A space link senses the SUMMED SCREEN PHASE (method (c)),
    because the downlink slab starts from a plane wave. Does that sum equal the
    phase that arrives at the aperture? The field carries the diffraction
    between the screens; the sum does not.

    THE MEASUREMENT. Both sources go through the SAME 21-mode projector. The
    screen source uses `estimate` on the unwrapped sum. The field source uses
    `estimate_from_slopes` on the wrapped-gradient slopes of the stored complex
    field. The study reports:
      - the tilt gain and the tilt RMS difference (Noll j = 2 and j = 3),
      - the 21-mode gain and RMS difference (j = 2 to 21; piston has no slope),
      - the RMS residual of the summed screen phase after the removal of the
        FIELD's own 21-mode fit, in rad.

    Args:
        camp:     the Campaign. It must hold the fields and the screen phase.
        n_trials: the number of trials to read.

    Returns:
        The pair (a result dict, a dict of the raw coefficient arrays).
    """
    result = camp.load(int(n_trials), fields=True)
    if result.screen_phase is None:
        raise ValueError(f"{camp.root_dir} holds no stored screen phase.")
    patch = camp.patch
    modes = aperture_modes(camp, V0_MODES)
    rx = camp.scenario.rx_terminal

    screen_c, slope_c = [], []
    residual, screen_rms, steps = [], [], []
    for row, array in _rebuilt_fields(result, rx.aperture_m, None):
        phase = _screen_map(result, row, patch)
        screen_c.append(modes.estimate(phase))
        sx, sy = wrapped_gradient(array)
        slope_c.append(modes.estimate_from_slopes(sx, sy))
        # The step statistic reads the APERTURE pixels only, so the dark
        # pixels outside the mask do not enter it.
        mx, my = wrapped_gradient(array, mask=modes.mask)
        steps.append(max_abs_step(mx, my))
        # The residual of the summed screen phase against the FIELD's own fit.
        # The variance removes the mean, so the piston never enters it.
        residual.append(np.sqrt(modes.residual_variance(phase, slope_c[-1])))
        inside = phase.ravel()[modes.indices]
        screen_rms.append(float(np.std(inside)))

    screen_c = np.asarray(screen_c)
    slope_c = np.asarray(slope_c)
    tilt_gain, tilt_diff, tilt_rms = _gain_and_rms(screen_c[:, 1:3],
                                                   slope_c[:, 1:3])
    all_gain, all_diff, all_rms = _gain_and_rms(screen_c[:, 1:],
                                                slope_c[:, 1:])
    residual = np.asarray(residual)
    screen_rms = np.asarray(screen_rms)
    row = {
        "elevation_deg": float(np.min(np.asarray(
            camp.geometry.elevation_deg, dtype=float))),
        "n_trials": int(screen_c.shape[0]),
        "n_modes": int(V0_MODES),
        "sigma2_R": float(np.sum(camp.plan.sigma2_r)),
        "r0_total_m": float(camp.plan.r0_total_m),
        "d_over_r0": float(rx.aperture_m / camp.plan.r0_total_m),
        "grid_n": int(camp.grid.n),
        "pixel_mm": float(camp.grid.pixel_m * 1e3),
        "aperture_px": float(rx.aperture_m / camp.grid.pixel_m),
        "tilt_gain": tilt_gain,
        "tilt_rms_diff_rad": tilt_diff,
        "tilt_rms_rad": tilt_rms,
        "tilt_rel_diff": tilt_diff / tilt_rms,
        "modes_gain": all_gain,
        "modes_rms_diff_rad": all_diff,
        "modes_rms_rad": all_rms,
        "modes_rel_diff": all_diff / all_rms,
        "screen_rms_rad": float(screen_rms.mean()),
        "residual_rms_rad": float(residual.mean()),
        "residual_over_screen": float((residual / screen_rms).mean()),
        "max_step_rad_mean": float(np.mean(steps)),
        "max_step_rad_worst": float(np.max(steps)),
    }
    return row, {"screen": screen_c, "slope": slope_c}


# ---------------------------------------------------------------------------
# V2: the Noll cross-check of the modal chain
# ---------------------------------------------------------------------------

def v2_row(camp, n_trials, modes_list=V2_MODES):
    """Measure the residual phase variance against Noll's Table IV.

    THE MEASUREMENT. For each J of `modes_list` the function fits the first J
    Noll modes to the summed screen phase over the aperture, and it takes the
    variance of what is left. The mean over the trials goes against
    Delta_J (D / r0)^(5/3).

    THE FRIED PARAMETER. D / r0 uses `plan.r0_total_m`, the PLANE-WAVE
    composite Fried parameter of the screen plan
    (`olb.waveoptics.turbulence.sampling._composite_r0`, Fried,
    DOI 10.1364/JOSA.56.001372). The caller cross-checks it against
    `olb.turbulence.ao.plane_wave_fried_parameter_profile`.

    Source: Noll 1976, DOI 10.1364/JOSA.66.000207, Table IV, printed p. 210.

    Args:
        camp:       the Campaign. It must hold the screen phase.
        n_trials:   the number of trials to read.
        modes_list: the removed-mode counts.

    Returns:
        A result dict.
    """
    # fields=True is necessary: `Campaign.load` gates the stored screen phase
    # on the same flag as the stored field patch.
    result = camp.load(int(n_trials), fields=True)
    if result.screen_phase is None:
        raise ValueError(f"{camp.root_dir} holds no stored screen phase.")
    patch = camp.patch
    rx = camp.scenario.rx_terminal
    projectors = {j: aperture_modes(camp, j) for j in modes_list}
    totals = {j: [] for j in modes_list}
    for row in range(result.screen_phase.shape[0]):
        phase = _screen_map(result, row, patch)
        for j, mo in projectors.items():
            totals[j].append(mo.residual_variance(phase, mo.estimate(phase)))

    r0 = float(camp.plan.r0_total_m)
    d_over_r0 = float(rx.aperture_m / r0)
    scale = d_over_r0 ** (5.0 / 3.0)
    rows = []
    for j in modes_list:
        meas = np.asarray(totals[j], dtype=float)
        book = noll_delta(j) * scale
        rows.append({
            "J": int(j),
            "measured_rad2": float(meas.mean()),
            "measured_half": _bootstrap(meas, lambda s: float(s.mean())),
            "noll_rad2": float(book),
            "ratio": float(meas.mean() / book),
            "delta_J": float(noll_delta(j)),
        })
    return {
        "elevation_deg": float(np.min(np.asarray(
            camp.geometry.elevation_deg, dtype=float))),
        "n_trials": int(result.screen_phase.shape[0]),
        "r0_total_m": r0,
        "d_over_r0": d_over_r0,
        "L0_m": L0_M,
        "rows": rows,
    }


def r0_cross_check(elevation_deg, camp):
    """Cross-check the plan Fried parameter against the analytic profile.

    `plan.r0_total_m` adds the screen Fried parameters of the plan. The
    analytic value integrates the Hufnagel-Valley profile over the slant path.
    The two must agree, because both are the plane-wave slant r0 of the same
    atmosphere.

    Sources: Fried, DOI 10.1364/JOSA.56.001372; Andrews and Phillips,
    DOI 10.1117/3.626196, Ch. 12, printed p. 481 (the airmass).

    Args:
        elevation_deg: the elevation, in deg.
        camp:          the Campaign.

    Returns:
        A dict of the two values and their ratio.
    """
    site = camp.scenario.channel.site
    profile = default_cn2_profile(site, DEFAULT_HS)
    analytic = plane_wave_fried_parameter_profile(profile, DEFAULT_HS, LAM,
                                                  float(elevation_deg))
    plan_r0 = float(camp.plan.r0_total_m)
    return {"plan_r0_m": plan_r0, "analytic_r0_m": float(analytic),
            "ratio": plan_r0 / float(analytic)}


def outer_scale_control(camp, n_seeds, modes_list=V2_MODES):
    """Measure the Noll ratio of ONE screen at L0 = 25 m and at L0 = inf.

    WHY. Noll's law assumes the Kolmogorov spectrum, which has an INFINITE
    outer scale. The campaign screens run at L0 = 25 m, so they hold less
    low-order power, and the measured J = 1 and J = 3 residuals must read LOW.
    This control isolates that effect: it makes ONE screen of the composite r0
    on the SAME grid, at both outer scales, and it fits the same ladder. It
    runs no propagation, so it separates the SPECTRUM from the split step.

    Source of the von Karman spectrum: Andrews and Phillips,
    DOI 10.1117/3.626196, Ch. 3.

    Args:
        camp:       the Campaign (it gives the grid and the composite r0).
        n_seeds:    the number of screens of each outer scale.
        modes_list: the removed-mode counts.

    Returns:
        A dict L0 label -> {J: ratio}.
    """
    from olb.waveoptics.turbulence.screens import ScreenFactory

    patch = camp.patch
    rx = camp.scenario.rx_terminal
    r0 = float(camp.plan.r0_total_m)
    scale = (rx.aperture_m / r0) ** (5.0 / 3.0)
    projectors = {j: aperture_modes(camp, j) for j in modes_list}
    out = {}
    for label, l0 in (("25", L0_M), ("inf", np.inf)):
        factory = ScreenFactory(patch.n, camp.grid.pixel_m, L0_m=l0)
        rng = np.random.default_rng(BOOT_SEED)
        totals = {j: 0.0 for j in modes_list}
        for _ in range(int(n_seeds)):
            screen = factory.make(r0, rng)
            for j, mo in projectors.items():
                totals[j] += mo.residual_variance(screen, mo.estimate(screen))
        out[label] = {int(j): float(totals[j] / n_seeds
                                    / (noll_delta(j) * scale))
                      for j in modes_list}
    return out


# ---------------------------------------------------------------------------
# V1: the fade of the corrected hero downlink
# ---------------------------------------------------------------------------

def v1_rows(camp):
    """Give the SMF fade and the bucket fade of one campaign.

    THE TWO QUANTITIES.
      - the SMF coupling loss -10 log10(smf_eta). The correction acts here.
      - the collected (bucket) power loss -10 log10(collected_power). The
        correction is a PHASE factor, so it must NOT move this number.
    The study also reports the composite -10 log10(power * eta), which is the
    fidelity-2 downlink SMF quantity of `olb.models.waveoptics`.

    Args:
        camp: the Campaign.

    Returns:
        A dict of the three fade rows.
    """
    result = camp.load(fields=False)
    power = np.array([t.collected_power for t in result.trials], dtype=float)
    eta = np.array([t.smf_eta for t in result.trials], dtype=float)
    return {
        "n_trials": int(eta.size),
        "n_modes_corrected": int(camp.n_modes_corrected),
        "smf": _fade_row(-10.0 * np.log10(eta)),
        "bucket": _fade_row(-10.0 * np.log10(power)),
        "composite": _fade_row(-10.0 * np.log10(power * eta)),
        "mean_eta": float(eta.mean()),
        "mean_power": float(power.mean()),
        "disk_bytes": _dir_bytes(camp.root_dir),
    }


def fast_row(elevation_deg, stack, n_samples, npxls_set, tol, say):
    """Measure the fidelity-1 FAST SMF Term of one stack.

    FAST reads the compensation stack of the receive terminal:
    an empty stack gives AO_MODE = NOAO, a TipTilt stage gives TT, and an AO(n)
    stage gives modal AO with ZMAX = n (`olb.models.fast._ao_params`). The
    outer scale matches the field, for parity.

    A GRID GUARD RUNS FIRST. FAST can undersample the low-order tilt, so the
    function tries each NPXLS and it takes the smallest one within `tol` dB of
    the largest one. Copied from validation/waveoptics_vs_fast/.

    Source: O. J. D. Farley and others, DOI 10.1364/OE.458659.

    Args:
        elevation_deg: the elevation, in deg.
        stack:         the compensation stack.
        n_samples:     the FAST draws.
        npxls_set:     the NPXLS values of the guard.
        tol:           the guard tolerance, in dB.
        say:           the log function.

    Returns:
        A dict of the FAST measurement, or None when fast-aosim is missing.
    """
    from olb.models.fast import smf_fast_term

    scn, _ = hero_scenario(elevation_deg, compensation=stack)
    geom = CircularOrbit(altitude_m=ALTITUDE_M,
                         elevation_deg=float(elevation_deg))

    def term_of(npxls):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return smf_fast_term(scn, geom, n_samples=int(n_samples),
                                 fast_params={"NPXLS": int(npxls),
                                              "L0": float(L0_M)})

    guard = []
    for npxls in npxls_set:
        t = term_of(npxls)
        guard.append({"npxls": int(npxls), "mean_db": float(t.mean_db)})
    reference = guard[-1]["mean_db"]
    pinned = next((g for g in guard
                   if abs(g["mean_db"] - reference) <= tol), guard[-1])
    say(f"    FAST grid guard: "
        + ", ".join(f"{g['npxls']}:{g['mean_db']:.3f}" for g in guard)
        + f"  -> NPXLS {pinned['npxls']}")
    term = term_of(pinned["npxls"])
    quant = {f"p{int(p * 100)}": float(term.quantile_db(1.0 - p))
             for p in EXCEEDANCE}
    return {
        "elevation_deg": float(elevation_deg),
        "mean_db": float(term.mean_db),
        "quantiles_db": quant,
        "npxls": int(pinned["npxls"]),
        "guard": guard,
        "n_samples": int(n_samples),
        "ao_mode": term.meta.get("ao_mode"),
        "zmax": term.meta.get("zmax"),
        "floor_db": float(term.meta.get("floor_db")),
        "r0_los_m": float(term.meta.get("r0_los_m")),
    }


def posthoc_equals_inrun(base_camp, comp_camp, stack, n_trials):
    """Check that the post-hoc correction equals the in-run correction.

    The stored patch of an UNCORRECTED campaign holds the raw field, so
    `recouple_compensated` must reproduce the coupling efficiency of a campaign
    that corrected the SAME seeds inside the run. The two campaigns share the
    seed and the plan, so trial k is the same atmosphere in both.

    Args:
        base_camp: the uncorrected Campaign with the stored screen phase.
        comp_camp: the in-run corrected Campaign of the same stack.
        stack:     the compensation stack.
        n_trials:  the number of trials to compare.

    Returns:
        A dict of the worst and the RMS relative difference.
    """
    n = int(min(n_trials, base_camp.n_stored, comp_camp.n_stored))
    detector = base_camp.scenario.rx_terminal.detector
    post = base_camp.recouple_compensated(list(stack), detector, n_trials=n,
                                          source="screens")
    inrun = np.array([t.smf_eta for t in
                      comp_camp.load(n, fields=False).trials], dtype=float)
    rel = np.abs(post - inrun) / inrun
    return {"n": n, "max_rel": float(rel.max()),
            "rms_rel": float(np.sqrt((rel ** 2).mean())),
            "mean_eta_posthoc": float(post.mean()),
            "mean_eta_inrun": float(inrun.mean())}


# ---------------------------------------------------------------------------
# V3: method (a) against method (c)
# ---------------------------------------------------------------------------

def v3_row(camp, stack_name, n_trials):
    """Correct one stored campaign from the screens and from the slopes.

    Both routes read the SAME stored field. The screen route reads the stored
    summed screen phase (method (c)); the slope route reads the wrapped
    gradient of the field (method (a)). They must agree on a SPACE link, where
    both are valid. The check qualifies the slope route before a terrestrial
    link trusts it alone.

    Args:
        camp:       the uncorrected Campaign with the stored screen phase.
        stack_name: "tiptilt" or "ao21".
        n_trials:   the number of trials.

    Returns:
        A dict of the agreement.
    """
    stack = list(STACKS[stack_name])
    detector = camp.scenario.rx_terminal.detector
    n = int(min(n_trials, camp.n_stored))
    screens = camp.recouple_compensated(stack, detector, n_trials=n,
                                        source="screens")
    slopes = camp.recouple_compensated(stack, detector, n_trials=n,
                                       source="slopes")
    rel = (slopes - screens) / screens
    d_db = -10.0 * np.log10(slopes / screens)
    return {
        "stack": stack_name,
        "n": n,
        "mean_eta_screens": float(screens.mean()),
        "mean_eta_slopes": float(slopes.mean()),
        "rel_rms": float(np.sqrt((rel ** 2).mean())),
        "rel_max": float(np.abs(rel).max()),
        "mean_db_screens": float((-10.0 * np.log10(screens)).mean()),
        "mean_db_slopes": float((-10.0 * np.log10(slopes)).mean()),
        "db_rms_diff": float(np.sqrt((d_db ** 2).mean())),
        "db_max_diff": float(np.abs(d_db).max()),
    }


# ---------------------------------------------------------------------------
# V4: the terrestrial sanity check
# ---------------------------------------------------------------------------

def v4_row(path_m, cn2, preset, n_trials, stacks, say):
    """Apply the slope correction post hoc to a stored 2-TC campaign.

    The function REOPENS the campaign of `validation/terrestrial_campaigns/`.
    It copies the construction of that module exactly, because the fingerprint
    must match. It computes NO trial: it reads the stored blocks only.

    The terrestrial source is the SLOPES, because a horizontal link is near
    field and the summed screen phase over the pupil stops matching the
    arriving phase.

    Args:
        path_m:   the path length, in m.
        cn2:      the Cn2 of the path, in m^-2/3.
        preset:   the sampling preset of the stored campaign.
        n_trials: the number of trials to read.
        stacks:   the stack names to apply.
        say:      the log function.

    Returns:
        A result dict.
    """
    from validation.terrestrial_campaigns.run_campaigns import (
        L0_M as TC_L0_M, PATCH_RADIUS_M as TC_PATCH_M, PRECISION as TC_PREC,
        SEED as TC_SEED, build_scenario, campaigns_root as tc_root, cell_tag)

    scn, geom = build_scenario(path_m, cn2, "collimated")
    root = os.path.join(tc_root(), cell_tag(path_m, cn2, preset,
                                            "collimated", False))
    camp = Campaign(scn, geom, root, seed=TC_SEED, preset=preset,
                    block_size=50, patch_radius_m=TC_PATCH_M, L0_m=TC_L0_M,
                    precision=TC_PREC)
    n = int(min(n_trials, camp.n_stored))
    say(f"  {cell_tag(path_m, cn2, preset, 'collimated', False)}: "
        f"{camp.n_stored} trials on disk, {n} read, grid {camp.grid.n} px, "
        f"{camp.plan.z_m.size} screens, r0_total "
        f"{camp.plan.r0_total_m * 100:.2f} cm")

    detector = camp.scenario.rx_terminal.detector
    result = camp.load(n, fields=True)
    eta = np.array([t.smf_eta for t in result.trials], dtype=float)

    # The sampling limit of method (a): the wrapped gradient aliases above pi.
    modes = aperture_modes(camp, 21)
    steps = []
    for _, array in _rebuilt_fields(result, camp.scenario.rx_terminal.aperture_m,
                                    None):
        sx, sy = wrapped_gradient(array, mask=modes.mask)
        steps.append(max_abs_step(sx, sy))
    steps = np.asarray(steps)

    out = {
        "cell": cell_tag(path_m, cn2, preset, "collimated", False),
        "path_m": float(path_m), "cn2": float(cn2), "preset": preset,
        "n_trials": int(n),
        "grid_n": int(camp.grid.n),
        "n_screens": int(camp.plan.z_m.size),
        "r0_total_m": float(camp.plan.r0_total_m),
        "step_warn_rad": float(SLOPE_STEP_WARN_RAD),
        "step_mean_rad": float(steps.mean()),
        "step_max_rad": float(steps.max()),
        "frac_over_warn": float((steps > SLOPE_STEP_WARN_RAD).mean()),
        "untracked": _fade_row(-10.0 * np.log10(eta)),
        "tracked": {},
    }
    for name in stacks:
        corrected = camp.recouple_compensated(list(STACKS[name]), detector,
                                              n_trials=n, source="slopes")
        out["tracked"][name] = _fade_row(-10.0 * np.log10(corrected))
    return out


# ---------------------------------------------------------------------------
# V5: method (a) against method (c) on the terrestrial link
# ---------------------------------------------------------------------------

# The V5 cells, as "<path>km:<cn2>" tokens. Three path lengths at one Cn2, and
# one stronger cell at the shortest path.
V5_CELLS = ("2km:3e-15", "5km:3e-15", "10km:3e-15", "2km:1e-14")

# The stacks of the V5 comparison.
V5_STACKS = ("tiptilt", "ao21")


def parse_v5_cell(token):
    """Turn a V5 cell token into the pair (path_m, cn2).

    Args:
        token: a token such as "5km:3e-15".

    Returns:
        The pair (the path length in m, the Cn2 in m^-2/3).

    Raises:
        ValueError: the token has no ":" separator, or a bad number.
    """
    if ":" not in token:
        raise ValueError(f"a V5 cell token needs a colon, got {token!r}. "
                         f"Use a token such as 5km:3e-15.")
    path_text, cn2_text = token.split(":", 1)
    return float(path_text.rstrip("kmKM")) * 1e3, float(cn2_text)


def v5_campaign(path_m, cn2, preset, args):
    """Build (or reopen) the terrestrial campaign of one V5 cell.

    THE CELL IS THE 2-TC CELL. The scenario, the seed, the outer scale, the
    patch radius and the precision come from
    `validation/terrestrial_campaigns/run_campaigns.py`, so this campaign is
    the 2-TC cell plus two settings: it stores the summed screen phase, and it
    runs on the CUDA backend. Each setting enters the campaign fingerprint, so
    the store gets its own root under `validation/waveoptics_ao/campaigns/`.

    Args:
        path_m: the horizontal path length, in m.
        cn2:    the Cn2 of the path, in m^-2/3.
        preset: the sampling preset name.
        args:   the parsed command line.

    Returns:
        The pair (Campaign, the sizer warning texts).
    """
    from validation.terrestrial_campaigns.run_campaigns import (
        L0_M as TC_L0_M, PATCH_RADIUS_M as TC_PATCH_M, PRECISION as TC_PREC,
        SEED as TC_SEED, build_scenario, cell_tag)

    scn, geom = build_scenario(path_m, cn2, "collimated")
    root = os.path.join(campaigns_root(),
                        "terr_" + cell_tag(path_m, cn2, preset, "collimated",
                                           False))
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        camp = Campaign(scn, geom, root, seed=TC_SEED, preset=preset,
                        block_size=int(args.block_size),
                        patch_radius_m=TC_PATCH_M, L0_m=TC_L0_M,
                        precision=TC_PREC, fft_backend=args.fft_backend,
                        store_screen_phase=True)
    return camp, sorted({str(w.message) for w in caught})


def v5_row(camp, n_trials):
    """Compare the two sensing sources on ONE terrestrial campaign.

    THE QUESTION. On a horizontal near-field link the summed screen phase is
    not guaranteed to be the arriving wavefront: the screens sit at different
    ranges, the beam footprint changes along the path, and the field
    diffracts between the screens. The field slopes are the arriving
    wavefront, while the phase step for each pixel stays under the stencil
    limit. So this row measures how far the two sources move apart.

    THE FOUR PARTS.
      (i)   the V3 comparison of the coupling efficiency, for each stack.
      (ii)  the V0 comparison of the modal coefficients, at 21 modes.
      (iii) the step statistic of the slope stencil, so the reader sees where
            the slope reference itself is past its limit.
      (iv)  the SMF fade quantiles of the untracked field and of both sources.

    Args:
        camp:     the Campaign. It must hold the fields and the screen phase.
        n_trials: the number of trials to read.

    Returns:
        A result dict.
    """
    n = int(min(n_trials, camp.n_stored))
    result = camp.load(n, fields=True)
    if result.screen_phase is None:
        raise ValueError(f"{camp.root_dir} holds no stored screen phase.")
    patch = camp.patch
    rx = camp.scenario.rx_terminal
    detector = rx.detector
    modes = aperture_modes(camp, V0_MODES)

    # (ii) and (iii): one pass over the stored fields.
    screen_c, slope_c, steps, residual, screen_rms = [], [], [], [], []
    for row, array in _rebuilt_fields(result, rx.aperture_m, None):
        phase = _screen_map(result, row, patch)
        screen_c.append(modes.estimate(phase))
        sx, sy = wrapped_gradient(array)
        slope_c.append(modes.estimate_from_slopes(sx, sy))
        mx, my = wrapped_gradient(array, mask=modes.mask)
        steps.append(max_abs_step(mx, my))
        residual.append(np.sqrt(modes.residual_variance(phase, slope_c[-1])))
        screen_rms.append(float(np.std(phase.ravel()[modes.indices])))
    screen_c = np.asarray(screen_c)
    slope_c = np.asarray(slope_c)
    steps = np.asarray(steps)
    residual = np.asarray(residual)
    screen_rms = np.asarray(screen_rms)
    tilt_gain, tilt_diff, tilt_rms = _gain_and_rms(screen_c[:, 1:3],
                                                   slope_c[:, 1:3])
    all_gain, all_diff, all_rms = _gain_and_rms(screen_c[:, 1:],
                                                slope_c[:, 1:])

    # (i) and (iv): the coupling efficiency of each stack and each source.
    eta = {"untracked": np.array([t.smf_eta for t in result.trials],
                                 dtype=float)}
    compare, fades = [], {"untracked": _fade_row(
        -10.0 * np.log10(eta["untracked"]))}
    for name in V5_STACKS:
        for source in ("screens", "slopes"):
            eta[f"{name}_{source}"] = camp.recouple_compensated(
                list(STACKS[name]), detector, n_trials=n, source=source)
            fades[f"{name}_{source}"] = _fade_row(
                -10.0 * np.log10(eta[f"{name}_{source}"]))
        a = eta[f"{name}_screens"]
        b = eta[f"{name}_slopes"]
        rel = (b - a) / a
        d_db = -10.0 * np.log10(b / a)
        compare.append({
            "stack": name, "n": n,
            "mean_eta_screens": float(a.mean()),
            "mean_eta_slopes": float(b.mean()),
            "rel_rms": float(np.sqrt((rel ** 2).mean())),
            "rel_max": float(np.abs(rel).max()),
            "mean_db_screens": float((-10.0 * np.log10(a)).mean()),
            "mean_db_slopes": float((-10.0 * np.log10(b)).mean()),
            "db_rms_diff": float(np.sqrt((d_db ** 2).mean())),
            "db_max_diff": float(np.abs(d_db).max()),
        })

    return {
        "n_trials": n,
        "n_modes": int(V0_MODES),
        "grid_n": int(camp.grid.n),
        "pixel_mm": float(camp.grid.pixel_m * 1e3),
        "n_screens": int(camp.plan.z_m.size),
        "sigma2_R": float(np.sum(camp.plan.sigma2_r)),
        "r0_total_m": float(camp.plan.r0_total_m),
        "d_over_r0": float(rx.aperture_m / camp.plan.r0_total_m),
        "aperture_px": float(rx.aperture_m / camp.grid.pixel_m),
        "tilt_gain": tilt_gain,
        "tilt_rms_diff_rad": tilt_diff,
        "tilt_rms_rad": tilt_rms,
        "tilt_rel_diff": tilt_diff / tilt_rms,
        "modes_gain": all_gain,
        "modes_rms_diff_rad": all_diff,
        "modes_rms_rad": all_rms,
        "modes_rel_diff": all_diff / all_rms,
        "screen_rms_rad": float(screen_rms.mean()),
        "residual_rms_rad": float(residual.mean()),
        "residual_over_screen": float((residual / screen_rms).mean()),
        "step_warn_rad": float(SLOPE_STEP_WARN_RAD),
        "step_mean_rad": float(steps.mean()),
        "step_max_rad": float(steps.max()),
        "frac_over_warn": float((steps > SLOPE_STEP_WARN_RAD).mean()),
        "compare": compare,
        "fades": fades,
        "disk_bytes": _dir_bytes(camp.root_dir),
    }


# ---------------------------------------------------------------------------
# The figures
# ---------------------------------------------------------------------------

def _pyplot():
    """Give the matplotlib pyplot module with the Agg backend, or None."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib is not installed; the figures are skipped.",
              flush=True)
        return None
    return plt


def _fig_dir():
    """Give the figure directory, and make it."""
    path = os.path.join(HERE, "figures")
    os.makedirs(path, exist_ok=True)
    return path


def plot_v0(coeffs, rows):
    """Draw the coefficient agreement of the two sensing sources."""
    plt = _pyplot()
    if plt is None:
        return []
    fig, axes = plt.subplots(1, len(coeffs), figsize=(5.5 * len(coeffs), 5.0),
                             squeeze=False)
    for ax, (el, c) in zip(axes[0], sorted(coeffs.items())):
        row = next(r for r in rows if abs(r["elevation_deg"] - el) < 1e-6)
        ax.plot(c["screen"][:, 1:3].ravel(), c["slope"][:, 1:3].ravel(),
                ".", ms=3, alpha=0.5, label="tilt (j = 2, 3)")
        ax.plot(c["screen"][:, 3:].ravel(), c["slope"][:, 3:].ravel(),
                ".", ms=3, alpha=0.5, label="j = 4 to 21")
        lim = float(np.abs(c["screen"][:, 1:]).max()) * 1.05
        ax.plot([-lim, lim], [-lim, lim], "k-", lw=0.8, label="y = x")
        ax.set_xlabel("summed screen coefficient [rad]")
        ax.set_ylabel("field slope coefficient [rad]")
        ax.set_title(f"{el:.0f} deg: tilt gain {row['tilt_gain']:.3f}")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
    fig.tight_layout()
    path = os.path.join(_fig_dir(), "v0_source_agreement.png")
    fig.savefig(path, dpi=140)
    plt.close(fig)
    return [path]


def plot_v2(results):
    """Draw the measured residual variance against Noll's law."""
    plt = _pyplot()
    if plt is None:
        return []
    fig, ax = plt.subplots(figsize=(7.5, 5.0))
    for res in results:
        j = [r["J"] for r in res["rows"]]
        ax.loglog(j, [r["measured_rad2"] for r in res["rows"]], "o-",
                  label=f"measured, {res['elevation_deg']:.0f} deg")
        ax.loglog(j, [r["noll_rad2"] for r in res["rows"]], "--",
                  label=f"Noll, {res['elevation_deg']:.0f} deg")
    ax.set_xlabel("removed Noll modes J")
    ax.set_ylabel("residual phase variance [rad^2]")
    ax.set_title("The modal residual against Noll (L0 = 25 m screens)")
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=8)
    fig.tight_layout()
    path = os.path.join(_fig_dir(), "v2_noll_residual.png")
    fig.savefig(path, dpi=140)
    plt.close(fig)
    return [path]


def plot_v1(survival):
    """Draw the SMF loss survival curve of each stack."""
    plt = _pyplot()
    if plt is None:
        return []
    fig, axes = plt.subplots(1, len(survival), figsize=(6.0 * len(survival),
                                                        5.0), squeeze=False)
    for ax, (el, curves) in zip(axes[0], sorted(survival.items())):
        for name in STACK_ORDER:
            if name not in curves:
                continue
            s = np.sort(np.asarray(curves[name], dtype=float))
            frac = 1.0 - np.arange(s.size) / s.size
            ax.semilogy(s, frac, label=name)
        ax.set_xlabel("SMF coupling loss [dB]")
        ax.set_ylabel("fraction of trials with a larger loss")
        ax.set_title(f"{el:.0f} deg")
        ax.grid(alpha=0.3, which="both")
        ax.legend(fontsize=8)
    fig.tight_layout()
    path = os.path.join(_fig_dir(), "v1_smf_survival.png")
    fig.savefig(path, dpi=140)
    plt.close(fig)
    return [path]


# ---------------------------------------------------------------------------
# The drivers
# ---------------------------------------------------------------------------

def _log_maker(name):
    """Make the log file of one study, and give its `say` function."""
    path = os.path.join(HERE, f"waveoptics_ao_{name}.log")
    with open(path, "w", encoding="utf-8"):
        pass

    def say(text=""):
        print(text, flush=True)
        with open(path, "a", encoding="utf-8") as fh:
            fh.write(text + "\n")

    return say, path


def _write_results(name, payload):
    """Write the JSON record of one study."""
    path = os.path.join(HERE, f"waveoptics_ao_{name}_results.json")
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=1)
    return path


def _header(say, args, title):
    """Print the common header of a study."""
    say(title)
    say(f"date          : 2026-09-07      seed: {SEED}      preset: {PRESET}")
    say(f"outer scale   : L0 = {L0_M:g} m (the owner rule, backlog 2-P5)")
    say(f"precision     : single      fft backend: {args.fft_backend}      "
        f"workers: {args.workers}")
    say("scenario      : downlink, 1550 nm, 500 km, 700 mm ground SMF, "
        "100 mm space terminal")
    say("caveat        : PERFECT AO. It is an ideal modal fit of one snapshot: "
        "no sensor noise, no lag, no anisoplanatism.")
    say()


def run_v0(args):
    """Run V0, the trust gate of the summed-screen source."""
    say, log_path = _log_maker("v0")
    _header(say, args, "V0: the summed screen phase against the field phase "
                       "(backlog 2-AO)")
    rows, coeffs, timing = [], {}, {}
    for el in args.elevations:
        say(f"ELEVATION {el:.0f} deg")
        camp, warns = campaign_of(el, "base", args)
        for w in warns:
            say(f"  sizer warning: {w}")
        timing[str(el)] = ensure_trials(camp, args.n_trials, args, say)
        row, c = v0_row(camp, min(args.n_trials, camp.n_stored))
        rows.append(row)
        coeffs[float(el)] = c
        say(f"  grid {row['grid_n']} px, {row['pixel_mm']:.2f} mm pixel, "
            f"aperture {row['aperture_px']:.0f} px, "
            f"sigma2_R {row['sigma2_R']:.4f}, D/r0 {row['d_over_r0']:.2f}")
        say()

    hdr = (f"{'elev':>6s}{'sigma2_R':>10s}{'D/r0':>7s}{'tilt gain':>11s}"
           f"{'tilt dRMS':>11s}{'tilt RMS':>10s}{'21 gain':>9s}"
           f"{'21 dRMS':>9s}{'resid':>8s}{'res/scr':>9s}{'step':>7s}")
    say("THE V0 TABLE (the coefficients are in rad)")
    say(hdr)
    say("-" * len(hdr))
    for r in rows:
        say(f"{r['elevation_deg']:6.0f}{r['sigma2_R']:10.4f}"
            f"{r['d_over_r0']:7.2f}{r['tilt_gain']:11.4f}"
            f"{r['tilt_rms_diff_rad']:11.3f}{r['tilt_rms_rad']:10.3f}"
            f"{r['modes_gain']:9.4f}{r['modes_rms_diff_rad']:9.3f}"
            f"{r['residual_rms_rad']:8.3f}{r['residual_over_screen']:9.3f}"
            f"{r['max_step_rad_worst']:7.2f}")
    say("  tilt gain = the slope tilt regressed on the screen tilt "
        "(1.0 = they agree).")
    say("  tilt dRMS = the RMS difference of the two tilt sets; tilt RMS = the "
        "RMS screen tilt.")
    say("  resid = the RMS of the summed screen phase after the removal of the "
        "FIELD's 21-mode fit.")
    say("  res/scr = that residual divided by the RMS screen phase over the "
        "aperture.")
    say("  step = the worst wrapped phase step of the field, in rad per pixel "
        f"(the warn level is {SLOPE_STEP_WARN_RAD:g}).")
    say()

    figs = plot_v0(coeffs, rows) if args.figures else []
    path = _write_results("v0", {"study": "V0", "date": "2026-09-07",
                                 "rows": rows, "timing": timing,
                                 "figures": figs})
    say(f"wrote {path}")
    say(f"wrote {log_path}")
    for f in figs:
        say(f"wrote {f}")


def run_v2(args):
    """Run V2, the Noll cross-check of the modal chain."""
    say, log_path = _log_maker("v2")
    _header(say, args, "V2: the residual phase variance against Noll's "
                       "Table IV (backlog 2-AO)")
    say("Noll: sigma^2 = Delta_J (D / r0)^(5/3). Source: Noll 1976, "
        "DOI 10.1364/JOSA.66.000207, Table IV, printed p. 210.")
    say("Delta_1 and Delta_3 are the printed values; a larger J takes the "
        "asymptotic 0.2944 J^(-sqrt(3)/2).")
    say()

    results, checks, controls = [], {}, {}
    for el in args.elevations:
        say(f"ELEVATION {el:.0f} deg")
        camp, warns = campaign_of(el, "base", args)
        for w in warns:
            say(f"  sizer warning: {w}")
        ensure_trials(camp, args.n_trials, args, say)
        res = v2_row(camp, min(args.n_trials, camp.n_stored))
        results.append(res)
        checks[str(el)] = r0_cross_check(el, camp)
        say(f"  r0: plan {checks[str(el)]['plan_r0_m'] * 100:.2f} cm, "
            f"analytic profile "
            f"{checks[str(el)]['analytic_r0_m'] * 100:.2f} cm, "
            f"ratio {checks[str(el)]['ratio']:.3f}")
        say(f"  D/r0 = {res['d_over_r0']:.2f}, "
            f"{res['n_trials']} trials")
        say(f"  {'J':>4}{'measured [rad^2]':>20}{'Noll [rad^2]':>16}"
            f"{'ratio':>9}")
        for r in res["rows"]:
            say(f"  {r['J']:>4}{r['measured_rad2']:>15.4f}"
                f" +-{r['measured_half']:<4.3f}{r['noll_rad2']:>16.4f}"
                f"{r['ratio']:>9.3f}")
        if args.control_seeds > 0:
            controls[str(el)] = outer_scale_control(camp, args.control_seeds)
            say(f"  THE OUTER-SCALE CONTROL ({args.control_seeds} single "
                f"screens on the same grid, r0 = "
                f"{camp.plan.r0_total_m * 100:.2f} cm):")
            say(f"  {'J':>4}{'ratio at L0 = 25 m':>22}"
                f"{'ratio at L0 = inf':>20}")
            for j in V2_MODES:
                say(f"  {j:>4}{controls[str(el)]['25'][j]:>22.3f}"
                    f"{controls[str(el)]['inf'][j]:>20.3f}")
        say()

    figs = plot_v2(results) if args.figures else []
    path = _write_results("v2", {"study": "V2", "date": "2026-09-07",
                                 "results": results, "r0_check": checks,
                                 "outer_scale_control": controls,
                                 "figures": figs})
    say(f"wrote {path}")
    say(f"wrote {log_path}")
    for f in figs:
        say(f"wrote {f}")


def run_v1(args):
    """Run V1, the fade of the corrected hero downlink."""
    say, log_path = _log_maker("v1")
    _header(say, args, "V1: the SMF fade of the corrected hero downlink "
                       "(backlog 2-AO)")
    say("The stacks: base (no correction), tiptilt (3 Noll modes), ao10 (10), "
        "ao21 (21).")
    say("pX is the loss EXCEEDED X percent of the time. A LARGE loss is a deep "
        "fade.")
    say()

    out = {"study": "V1", "date": "2026-09-07", "elevations": {},
           "fast": {}, "posthoc_check": {}, "timing": {}}
    survival = {}
    for el in args.elevations:
        say(f"ELEVATION {el:.0f} deg")
        camps, rows, curves = {}, {}, {}
        for name in STACK_ORDER:
            camp, warns = campaign_of(el, name, args)
            camps[name] = camp
            for w in warns:
                say(f"  {name}: sizer warning: {w}")
            say(f"  RUN {name} ({camp.n_modes_corrected} Noll modes)")
            t = ensure_trials(camp, args.n_trials, args, say)
            out["timing"][f"el{el:.0f}_{name}"] = t
            rows[name] = v1_rows(camp)
            eta = np.array([tr.smf_eta for tr in
                            camp.load(fields=False).trials], dtype=float)
            curves[name] = -10.0 * np.log10(eta)
        out["elevations"][str(el)] = rows
        survival[float(el)] = curves

        say()
        say(f"  THE SMF COUPLING LOSS -10 log10(smf_eta) [dB], "
            f"{rows['base']['n_trials']} trials")
        hdr = (f"  {'stack':<9s}{'modes':>6s}{'mean':>8s}{'p50':>8s}"
               f"{'p10':>8s}{'p5':>8s}{'p1':>8s}{'d_p5':>7s}{'d_p1':>7s}"
               f"{'gain@p5':>9s}")
        say(hdr)
        say("  " + "-" * (len(hdr) - 2))
        base_q = rows["base"]["smf"]["quantiles_db"]
        for name in STACK_ORDER:
            q = rows[name]["smf"]["quantiles_db"]
            say(f"  {name:<9s}{rows[name]['n_modes_corrected']:>6d}"
                f"{rows[name]['smf']['mean_db']:>8.2f}{q['p50']:>8.2f}"
                f"{q['p10']:>8.2f}{q['p5']:>8.2f}{q['p1']:>8.2f}"
                f"{q['p5'] - q['p50']:>7.2f}{q['p1'] - q['p50']:>7.2f}"
                f"{base_q['p5'] - q['p5']:>9.2f}")
        say("  d_p5 and d_p1 are the FADE DEPTH pX - p50. gain@p5 is the p5 "
            "improvement against the base.")
        say()
        say(f"  THE BUCKET LOSS -10 log10(collected_power) [dB]. It MUST NOT "
            f"move: the correction is a phase factor.")
        say(f"  {'stack':<9s}{'mean':>8s}{'p50':>8s}{'p5':>8s}{'p1':>8s}")
        for name in STACK_ORDER:
            q = rows[name]["bucket"]["quantiles_db"]
            say(f"  {name:<9s}{rows[name]['bucket']['mean_db']:>8.3f}"
                f"{q['p50']:>8.3f}{q['p5']:>8.3f}{q['p1']:>8.3f}")
        say()

        # The post-hoc against the in-run check.
        chk = posthoc_equals_inrun(camps["base"], camps["tiptilt"],
                                   STACKS["tiptilt"],
                                   min(args.check_trials, args.n_trials))
        out["posthoc_check"][str(el)] = chk
        say(f"  POST-HOC AGAINST IN-RUN (tiptilt, {chk['n']} trials): "
            f"worst relative difference {chk['max_rel']:.2e}, "
            f"RMS {chk['rms_rel']:.2e}")
        say()

        # The fidelity-1 FAST comparison.
        if not args.no_fast:
            say("  THE FIDELITY-1 FAST COMPARISON (the same L0, the same "
                "stack)")
            fast_rows = {}
            for name in STACK_ORDER:
                try:
                    fast_rows[name] = fast_row(el, STACKS[name],
                                               args.fast_samples, NPXLS_SET,
                                               NPXLS_TOL, say)
                except ImportError as exc:
                    say(f"    fast-aosim is not available: {exc}")
                    break
                except Exception as exc:            # noqa: BLE001
                    say(f"    FAST failed for {name}: "
                        f"{type(exc).__name__}: {exc}")
                    continue
            out["fast"][str(el)] = fast_rows
            if fast_rows:
                say(f"  {'stack':<9s}{'ao_mode':>9s}{'zmax':>6s}"
                    f"{'FAST mean':>11s}{'field mean':>12s}{'gap':>8s}"
                    f"{'FAST p5':>9s}{'field p5':>10s}")
                for name in STACK_ORDER:
                    if name not in fast_rows:
                        continue
                    f = fast_rows[name]
                    # The field composite matches the FAST quantity: FAST
                    # reports the SMF coupling loss of the whole receive path.
                    fld = rows[name]["composite"]
                    say(f"  {name:<9s}{str(f['ao_mode']):>9s}"
                        f"{str(f['zmax']):>6s}{f['mean_db']:>11.2f}"
                        f"{fld['mean_db']:>12.2f}"
                        f"{f['mean_db'] - fld['mean_db']:>8.2f}"
                        f"{f['quantiles_db']['p5']:>9.2f}"
                        f"{fld['quantiles_db']['p5']:>10.2f}")
                say("  gap = FAST minus field, on the composite "
                    "-10 log10(power * eta). A positive gap means the field "
                    "reads LESS loss.")
            say()

    figs = plot_v1(survival) if args.figures else []
    out["figures"] = figs
    path = _write_results("v1", out)
    say(f"wrote {path}")
    say(f"wrote {log_path}")
    for f in figs:
        say(f"wrote {f}")


def run_v3(args):
    """Run V3, method (a) against method (c) on the space link."""
    say, log_path = _log_maker("v3")
    _header(say, args, "V3: the slope source against the summed-screen source "
                       "(backlog 2-AO)")
    say("Both routes correct the SAME stored field of the SAME uncorrected "
        "campaign. Only the sensing source changes.")
    say()
    rows = []
    for el in args.elevations:
        say(f"ELEVATION {el:.0f} deg")
        camp, _ = campaign_of(el, "base", args)
        ensure_trials(camp, args.n_trials, args, say)
        for name in ("tiptilt", "ao21"):
            row = v3_row(camp, name, min(args.n_trials, camp.n_stored))
            row["elevation_deg"] = float(el)
            rows.append(row)
        say()

    hdr = (f"{'elev':>6s}{'stack':>9s}{'n':>6s}{'eta screens':>13s}"
           f"{'eta slopes':>12s}{'rel RMS':>10s}{'rel max':>10s}"
           f"{'dB RMS':>9s}{'dB max':>9s}")
    say("THE V3 TABLE")
    say(hdr)
    say("-" * len(hdr))
    for r in rows:
        say(f"{r['elevation_deg']:6.0f}{r['stack']:>9s}{r['n']:6d}"
            f"{r['mean_eta_screens']:13.5f}{r['mean_eta_slopes']:12.5f}"
            f"{r['rel_rms']:10.4f}{r['rel_max']:10.4f}"
            f"{r['db_rms_diff']:9.3f}{r['db_max_diff']:9.3f}")
    say("  rel = (slopes - screens) / screens on the coupling efficiency of "
        "each trial.")
    say("  dB = -10 log10(slopes / screens), the per-trial coupling "
        "difference.")
    say()
    path = _write_results("v3", {"study": "V3", "date": "2026-09-07",
                                 "rows": rows})
    say(f"wrote {path}")
    say(f"wrote {log_path}")


def run_v4(args):
    """Run V4, the terrestrial sanity check."""
    say, log_path = _log_maker("v4")
    say("V4: the post-hoc slope correction of two stored terrestrial "
        "campaigns (backlog 2-AO, 1-9)")
    say("date          : 2026-09-07")
    say("source        : validation/terrestrial_campaigns/ (2-TC), seed "
        "20260906, L0 = 25 m, single precision")
    say("method        : the SLOPES of the stored field. A horizontal link is "
        "near field, so the summed screen phase is not the arriving phase.")
    say("caveat        : SANITY ONLY. It is a perfect snapshot fit with no "
        "sensor noise and no lag. The numbers feed backlog 1-9.")
    say()
    rows = []
    for path_m in args.paths_km:
        rows.append(v4_row(float(path_m) * 1e3, args.cn2, args.preset,
                           args.n_trials, ("tiptilt", "ao21"), say))
        say()

    hdr = (f"{'cell':<24s}{'n':>6s}  {'kind':<10s}{'mean':>8s}{'p50':>8s}"
           f"{'p10':>8s}{'p5':>8s}{'p1':>8s}{'d_p5':>7s}")
    say("THE V4 TABLE: the SMF coupling loss -10 log10(smf_eta) [dB]")
    say(hdr)
    say("-" * len(hdr))
    for r in rows:
        for kind, fade in [("untracked", r["untracked"])] + \
                [(k, v) for k, v in r["tracked"].items()]:
            q = fade["quantiles_db"]
            say(f"{r['cell']:<24s}{fade['n']:6d}  {kind:<10s}"
                f"{fade['mean_db']:8.2f}{q['p50']:8.2f}{q['p10']:8.2f}"
                f"{q['p5']:8.2f}{q['p1']:8.2f}{q['p5'] - q['p50']:7.2f}")
    say()
    say("THE SAMPLING LIMIT OF THE SLOPE METHOD")
    say(f"{'cell':<24s}{'mean step':>11s}{'worst step':>12s}"
        f"{'over warn':>11s}")
    for r in rows:
        say(f"{r['cell']:<24s}{r['step_mean_rad']:11.3f}"
            f"{r['step_max_rad']:12.3f}{r['frac_over_warn'] * 100:10.1f}%")
    say(f"  step = the worst wrapped phase difference of one trial, in rad per "
        f"pixel. The warn level is {SLOPE_STEP_WARN_RAD:g}, and the method "
        f"aliases above pi.")
    say()
    path = _write_results("v4", {"study": "V4", "date": "2026-09-07",
                                 "rows": rows})
    say(f"wrote {path}")
    say(f"wrote {log_path}")


def v5_verdict(row):
    """Give the plain-language verdict of one V5 cell.

    THE RULE. The field slopes are the arriving wavefront while the phase step
    for each pixel stays under the stencil limit. The summed screens are exact
    only in the geometric-optics limit, which a near-field horizontal path does
    not hold. So the verdict reads the step statistic first.

    Args:
        row: a `v5_row` result dict.

    Returns:
        A list of text lines.
    """
    over = row["frac_over_warn"] * 100.0
    lines = []
    if over <= 1.0:
        lines.append(
            f"    TRUST THE SLOPES. The worst phase step is "
            f"{row['step_max_rad']:.2f} rad for each pixel, and only "
            f"{over:.1f} percent of the trials go over the "
            f"{row['step_warn_rad']:g} rad warn level. The slope source "
            f"measures the wavefront that ARRIVES.")
    else:
        lines.append(
            f"    THE SLOPE SOURCE IS AT ITS LIMIT. {over:.1f} percent of the "
            f"trials carry a step above the {row['step_warn_rad']:g} rad warn "
            f"level, and the worst step is {row['step_max_rad']:.2f} rad. The "
            f"wrapped gradient aliases above pi, so read this cell with care.")
    lines.append(
        f"    The summed screens are NOT the arriving phase here: the path is "
        f"near field, sigma2_R is {row['sigma2_R']:.3f}, and the screens hold "
        f"no diffraction between the planes.")
    lines.append(
        f"    The 21-mode gain of the slopes against the screens is "
        f"{row['modes_gain']:.3f}, and the RMS difference is "
        f"{row['modes_rms_diff_rad']:.3f} rad against an RMS screen "
        f"coefficient of {row['modes_rms_rad']:.3f} rad "
        f"({row['modes_rel_diff'] * 100:.0f} percent).")
    return lines


def run_v5(args):
    """Run V5, the two sensing sources on the terrestrial link."""
    say, log_path = _log_maker("v5")
    say("V5: the summed-screen source against the field slopes on the "
        "TERRESTRIAL link (backlog 2-AO)")
    say("date          : 2026-09-07")
    say(f"cells         : {', '.join(args.v5_cells)} at the "
        f"{args.preset} preset")
    say("scenario      : the 2-TC cell. A collimated 5 mm launch, a 10 cm SMF "
        "receiver, L0 = 25 m, single precision, seed 20260906.")
    say("difference    : this campaign STORES the summed screen phase, so it "
        "gets its own root. A 2-TC store holds no screen phase.")
    say(f"precision     : single      fft backend: {args.fft_backend}")
    say("physics       : the field slopes are the arriving wavefront while "
        "the phase step for each pixel stays under the stencil limit.")
    say("                The summed screens are exact only in the "
        "geometric-optics limit, and a horizontal path is near field.")
    say("caveat        : PERFECT AO. It is an ideal modal fit of one "
        "snapshot: no sensor noise, no lag, no anisoplanatism.")
    say()

    rows, timing = [], {}
    for token in args.v5_cells:
        path_m, cn2 = parse_v5_cell(token)
        say(f"CELL {token}")
        camp, warns = v5_campaign(path_m, cn2, args.preset, args)
        for w in warns:
            say(f"  sizer warning: {w}")
        say(f"  root {camp.root_dir}")
        say(f"  grid {camp.grid.n} px, "
            f"{camp.grid.pixel_m * 1e3:.2f} mm pixel, "
            f"{camp.plan.z_m.size} screens, r0_total "
            f"{camp.plan.r0_total_m * 100:.2f} cm, "
            f"sigma2_R {np.sum(camp.plan.sigma2_r):.4f}")
        timing[token] = ensure_trials(camp, args.n_trials, args, say)
        row = v5_row(camp, args.n_trials)
        row.update({"cell": token, "path_m": float(path_m),
                    "cn2": float(cn2), "preset": args.preset})
        rows.append(row)
        say(f"  {row['n_trials']} trials read, "
            f"{row['disk_bytes'] / 2 ** 20:.0f} MB on disk")
        say()

    hdr = (f"{'cell':<12s}{'stack':>9s}{'n':>6s}{'eta screens':>13s}"
           f"{'eta slopes':>12s}{'rel RMS':>10s}{'rel max':>10s}"
           f"{'dB RMS':>9s}{'dB max':>9s}")
    say("THE V5 COUPLING TABLE (the V3 measurement on a horizontal path)")
    say(hdr)
    say("-" * len(hdr))
    for r in rows:
        for c in r["compare"]:
            say(f"{r['cell']:<12s}{c['stack']:>9s}{c['n']:6d}"
                f"{c['mean_eta_screens']:13.5f}{c['mean_eta_slopes']:12.5f}"
                f"{c['rel_rms']:10.4f}{c['rel_max']:10.4f}"
                f"{c['db_rms_diff']:9.3f}{c['db_max_diff']:9.3f}")
    say("  rel = (slopes - screens) / screens on the coupling efficiency of "
        "each trial.")
    say("  dB = -10 log10(slopes / screens), the per-trial coupling "
        "difference.")
    say()

    hdr = (f"{'cell':<12s}{'sigma2_R':>10s}{'D/r0':>7s}{'tilt gain':>11s}"
           f"{'tilt dRMS':>11s}{'tilt RMS':>10s}{'21 gain':>9s}"
           f"{'21 dRMS':>9s}{'21 RMS':>9s}{'resid':>8s}{'res/scr':>9s}")
    say("THE V5 COEFFICIENT TABLE (the V0 measurement, the values are in rad)")
    say(hdr)
    say("-" * len(hdr))
    for r in rows:
        say(f"{r['cell']:<12s}{r['sigma2_R']:10.4f}{r['d_over_r0']:7.2f}"
            f"{r['tilt_gain']:11.4f}{r['tilt_rms_diff_rad']:11.3f}"
            f"{r['tilt_rms_rad']:10.3f}{r['modes_gain']:9.4f}"
            f"{r['modes_rms_diff_rad']:9.3f}{r['modes_rms_rad']:9.3f}"
            f"{r['residual_rms_rad']:8.3f}{r['residual_over_screen']:9.3f}")
    say("  tilt gain = the slope tilt regressed on the screen tilt "
        "(1.0 = they agree).")
    say("  resid = the RMS of the summed screen phase after the removal of "
        "the FIELD's 21-mode fit.")
    say("  res/scr = that residual divided by the RMS screen phase over the "
        "aperture.")
    say()

    hdr = (f"{'cell':<12s}{'mean step':>11s}{'worst step':>12s}"
           f"{'over warn':>11s}{'grid':>7s}{'px mm':>8s}")
    say("THE SAMPLING LIMIT OF THE SLOPE METHOD")
    say(hdr)
    say("-" * len(hdr))
    for r in rows:
        say(f"{r['cell']:<12s}{r['step_mean_rad']:11.3f}"
            f"{r['step_max_rad']:12.3f}{r['frac_over_warn'] * 100:10.1f}%"
            f"{r['grid_n']:7d}{r['pixel_mm']:8.2f}")
    say(f"  step = the worst wrapped phase difference of one trial, in rad "
        f"for each pixel. The warn level is {SLOPE_STEP_WARN_RAD:g}, and the "
        f"method aliases above pi.")
    say()

    kinds = ["untracked"] + [f"{s}_{src}" for s in V5_STACKS
                             for src in ("screens", "slopes")]
    hdr = (f"{'cell':<12s}  {'kind':<18s}{'mean':>8s}{'p50':>8s}{'p10':>8s}"
           f"{'p5':>8s}{'p1':>8s}{'d_p5':>7s}")
    say("THE SMF FADE: the coupling loss -10 log10(smf_eta) [dB]")
    say(hdr)
    say("-" * len(hdr))
    for r in rows:
        for kind in kinds:
            fade = r["fades"][kind]
            q = fade["quantiles_db"]
            say(f"{r['cell']:<12s}  {kind:<18s}{fade['mean_db']:8.2f}"
                f"{q['p50']:8.2f}{q['p10']:8.2f}{q['p5']:8.2f}{q['p1']:8.2f}"
                f"{q['p5'] - q['p50']:7.2f}")
    say("  pX is the loss EXCEEDED X percent of the time. d_p5 is the fade "
        "depth p5 - p50.")
    say()

    say("THE VERDICT OF EACH CELL")
    for r in rows:
        say(f"  {r['cell']}:")
        for line in v5_verdict(r):
            say(line)
        say()

    path = _write_results("v5", {"study": "V5", "date": "2026-09-07",
                                 "preset": args.preset, "rows": rows,
                                 "timing": timing})
    say(f"wrote {path}")
    say(f"wrote {log_path}")


STUDIES = {"v0": run_v0, "v1": run_v1, "v2": run_v2, "v3": run_v3,
           "v4": run_v4, "v5": run_v5}


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--study", nargs="+", default=["v0"],
                    help="the studies to run, from v0 v1 v2 v3 v4 v5")
    ap.add_argument("--elevations", nargs="+", type=float, default=[30.0, 20.0],
                    help="the elevations of the space studies [deg]")
    ap.add_argument("--n-trials", type=int, default=200,
                    help="the trials of each campaign")
    ap.add_argument("--check-trials", type=int, default=200,
                    help="the trials of the post-hoc against in-run check")
    ap.add_argument("--block-size", type=int, default=50,
                    help="the trials in one block file")
    ap.add_argument("--workers", default=None,
                    help="the campaign pool size. Leave it out for the CUDA "
                         "backend: one device runs one stream.")
    ap.add_argument("--fft-backend", default="cupy",
                    choices=["numpy", "scipy", "cupy"],
                    help="the FFT backend of the campaigns")
    ap.add_argument("--analyse-only", action="store_true",
                    help="read what is stored and compute no trial")
    ap.add_argument("--no-figures", dest="figures", action="store_false",
                    help="skip the figures")
    ap.add_argument("--control-seeds", type=int, default=200,
                    help="the single screens of the V2 outer-scale control")
    ap.add_argument("--fast-samples", type=int, default=1000,
                    help="the FAST Monte Carlo draws of the V1 comparison")
    ap.add_argument("--no-fast", action="store_true",
                    help="skip the fidelity-1 FAST comparison")
    ap.add_argument("--paths-km", nargs="+", type=float, default=[2.0, 10.0],
                    help="the V4 terrestrial path lengths [km]")
    ap.add_argument("--cn2", type=float, default=3e-15,
                    help="the V4 terrestrial Cn2 [m^-2/3]")
    ap.add_argument("--preset", default="rapid",
                    help="the preset of the V4 and V5 terrestrial campaigns")
    ap.add_argument("--v5-cells", nargs="+", default=list(V5_CELLS),
                    help="the V5 terrestrial cells, as <path>km:<cn2> tokens, "
                         "for example 5km:3e-15")
    args = ap.parse_args()
    if args.workers not in (None, "auto"):
        args.workers = int(args.workers)

    for name in args.study:
        if name not in STUDIES:
            raise SystemExit(f"unknown study {name!r}. Use one of "
                             f"{sorted(STUDIES)}.")
    for name in args.study:
        STUDIES[name](args)


if __name__ == '__main__':
    main()
