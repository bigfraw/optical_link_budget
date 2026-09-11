"""The point-ahead validation of the fidelity-2 uplink (backlog 2-P4).

WHAT IT MEASURES. The fidelity-2 runner now propagates the slab a second time
for each POINT-AHEAD angle
(`olb.waveoptics.turbulence.run.propagate_turbulent_scenario`,
`point_ahead_rad`). The ground terminal senses the downlink beacon, it applies
the conjugate wavefront to the uplink beam, and the uplink goes to where the
satellite WILL BE. The two directions read the SAME atmosphere through a
LATERALLY SHIFTED window of each screen. The drop of the reciprocity overlap
between the two directions is the point-ahead anisoplanatism.

THE FOUR STUDIES.
  p0  the identities. The stored angle 0.0 must equal the beacon-direction
      overlap bit for bit; the post-hoc read of the stored planes must equal
      an in-run corrected campaign; and the regeneration route must give the
      stored column back.
  p1  the penalty against the other two rungs. It gives the mean loss of each
      angle, and it compares the anisoplanatic penalty L(theta) - L(0) with
      the fidelity-1 FAST Term and with the analytic Stone Term.
  p2  the fade. It gives p50, p5 and p1 of the per-trial loss at each angle
      and each stack.
  p3  the screen count. It repeats the 30 deg penalty on plans of 5, 9, 15 and
      25 screens and on a ground-split plan, at the PINNED production grid.

THE CASE. The hero uplink of `validation/anisoplanatism_screens/common.py`: a
1550 nm uplink from a 700 mm ground terminal with a full-aperture launch
(waist 0.35 m) to a 100 mm space terminal at 500 km, with a `DownlinkBeacon`
pre-compensation source. It is the hero downlink of
`validation/waveoptics_ao/` turned around.

THE TILT STAYS IN (owner decision, 2026-09-11). The terminal senses the
DOWNLINK beacon tilt and the steering mirror adds the point-ahead offset
geometrically, so the uplink has no tilt reference of its own. The runner
corrects the modes of the stack and it keeps the tilt in the error; FAST keeps
the piston and the tilt in its modal mask; the Stone Term takes
`remove='piston'`. So the three rungs read the SAME mode set, except the
piston, which changes no overlap integral.

THE OUTER SCALE IS 25 m ON EVERY RUN (the owner rule of 2026-09-05, backlog
2-P5), on the screens, in FAST (`fast_params={"L0": 25.0}`) and in the Stone
Term (`L0_m=25.0`).

PERFECT AO, AND NOTHING MORE. The correction is an ideal modal fit of one
snapshot: no wavefront-sensor noise, no finite subaperture, no aliasing, no
servo lag and no branch point. So every corrected number here is the UPPER
BOUND of the benefit of a corrector, and the point-ahead penalty is the
CLEANEST possible one.

Sources:
- J. Stone, P. H. Hu, S. P. Mills and S. Ma, "Anisoplanatic effects in
  finite-aperture optical systems," J. Opt. Soc. Am. A 11(1), 347-357 (1994),
  DOI 10.1364/JOSAA.11.000347. The analytic angular phase variance.
- J. H. Shapiro, "Reciprocity of the turbulent atmosphere," J. Opt. Soc. Am.
  61(4), 492-495 (1971), DOI 10.1364/JOSA.61.000492. The uplink overlap.
- R. J. Noll, J. Opt. Soc. Am. 66(3), 207-211 (1976),
  DOI 10.1364/JOSA.66.000207, Table I. The mode order of the correction.
- O. J. D. Farley and others, Opt. Express 30(13), 23050 (2022),
  DOI 10.1364/OE.458659. The FAST method (fidelity 1).
- Schmidt, Numerical Simulation of Optical Wave Propagation (2010),
  DOI 10.1117/3.866274, Ch. 9. The split-step method of the trials.
- Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196. Ch. 12 the
  Hufnagel-Valley profile, Ch. 8 the Rytov variance.
- T. S. Ross, Appl. Opt. 48(10), 1812 (2009), DOI 10.1364/AO.48.001812. The
  extended Marechal mapping of the Stone Term.

Run it from the repository root. The campaigns take the CUDA backend:

    python -m validation.waveoptics_pointahead.waveoptics_pointahead \\
        --study p0 p1 p2 p3

Add `--analyse-only` to read what is stored and to compute no trial. Add
`--dry-run` to print the campaign list and to exit. Add `--fft-backend numpy`
to run on a host with no CUDA device, and `--no-fast` when `fast-aosim` is not
installed.
"""

import argparse
import json
import os
import warnings
from dataclasses import replace

import numpy as np

from olb.geometry import CircularOrbit
from olb.terminal import AO
from olb.waveoptics.turbulence import Campaign
# The case and the screen plans of the point-ahead family. They are IMPORTED,
# not copied: `validation/anisoplanatism_screens/` owns the hero uplink.
from validation.anisoplanatism_screens.common import (ARCSEC, L0_M,
                                                      hero_uplink, plan_set)
# The driver helpers of the perfect-AO study. They are IMPORTED too: the fade
# row, the bootstrap, the trial runner and the FAST grid guard are the same.
from validation.waveoptics_ao.waveoptics_ao import (BOOT_INTERVAL, BOOT_SEED,
                                                    EXCEEDANCE, NPXLS_SET,
                                                    NPXLS_TOL, N_BOOT, STACKS,
                                                    STACK_ORDER, _dir_bytes,
                                                    _fade_row, _pyplot,
                                                    ensure_trials)

HERE = os.path.dirname(os.path.abspath(__file__))
ENV_ROOT = "OLB_WAVEOPTICS_POINTAHEAD_ROOT"

# ---------------------------------------------------------------------------
# The case
# ---------------------------------------------------------------------------

SEED = 20260907
PRESET = "standard"

# The stored patch radius, in m. It is half of the 0.7 m ground aperture with
# the 1.5x margin of the campaign default (PATCH_MARGIN_FACTOR).
PATCH_RADIUS_M = 0.525

# The FIXED point-ahead angles, in arcsec. The geometric angle of the
# elevation goes in front of them, so every campaign holds five angles:
# 0 (the beacon direction), the geometric angle, and these three.
FIXED_ARCSEC = (2.0, 5.0, 10.0)

# The screen counts of p3, and the sub-screen count of the ground-split plan.
COUNTS = (5, 9, 15, 25)
SPLIT_N = 4
P3_PLANS = ("n5", "n9", "n15", "n25", f"ground_split_x{SPLIT_N}")

# The pass band of the p1 comparison, in dB. It is the like-for-like tolerance
# of docs/physics.md Section 9l.
FAST_PASS_DB = 0.5

# The measured GPU rate of the 2-AO campaigns: 8200 trials in about 11 minutes
# at 512 px and 9 screens, so 12.4 trials/s for ONE split step per trial. A
# point-ahead trial runs 1 + n_angles split steps.
GPU_TRIALS_PER_S = 12.4


def campaigns_root():
    """Give the parent directory of every campaign of this study."""
    return os.environ.get(ENV_ROOT) or os.path.join(HERE, "campaigns")


def geometric_angle(elevation_deg):
    """Give the point-ahead angle of one elevation, in rad.

    The value comes from `olb.geometry.CircularOrbit.point_ahead_rad`.
    """
    _scn, geom = hero_uplink(elevation_deg)
    return float(np.atleast_1d(geom.point_ahead_rad).ravel()[0])


def angles_of(elevation_deg):
    """Give the point-ahead angles of one elevation, in rad.

    The order is the order of every stored column: 0.0 (the beacon
    direction), the geometric angle of the elevation, then the fixed angles.

    Args:
        elevation_deg: the elevation, in deg.

    Returns:
        A tuple of floats, in rad.
    """
    return (0.0, geometric_angle(elevation_deg)) + tuple(
        a * ARCSEC for a in FIXED_ARCSEC)


def angle_labels(angles):
    """Give the short label of each angle, for a table column."""
    out = []
    for i, a in enumerate(angles):
        if a == 0.0:
            out.append("beacon")
        elif i == 1:
            out.append(f"geom {a / ARCSEC:.2f}\"")
        else:
            out.append(f"{a / ARCSEC:.0f}\"")
    return out


class AngleGeometry:
    """A geometry of ONE elevation and ONE chosen point-ahead angle.

    WHY IT EXISTS. `olb.links.uplink.uplink_point_ahead_term` reads the angle
    from `geometry.point_ahead_rad`, and that property of a `CircularOrbit` is
    the GEOMETRIC angle of the elevation. This study sweeps the angle at a
    FIXED elevation, so it hands the Term a small object that carries the two
    attributes the Term reads. The FAST Term needs no such object: its angle
    goes through `fast_params={"DTHETA": [arcsec, 0]}`, which the Term merges
    last.

    Attributes:
        elevation_deg:   the scalar elevation, in deg.
        point_ahead_rad: the scalar point-ahead angle, in rad.
    """

    def __init__(self, elevation_deg, point_ahead_rad):
        self.elevation_deg = float(elevation_deg)
        self.point_ahead_rad = float(point_ahead_rad)


def stacked_scenario(elevation_deg, stack):
    """Give the hero uplink with one compensation stack on the ground.

    The CAMPAIGN never uses this scenario: a campaign takes the plain
    scenario and it declares the stack through `Campaign(compensation=...)`,
    so every stack of one elevation shares the same atmosphere key. The
    ANALYTIC Terms read the stack from the terminal, so they take this copy.

    Args:
        elevation_deg: the elevation, in deg.
        stack:         a sequence of TipTilt and AO stages.

    Returns:
        The SpaceScenario.
    """
    scn, _geom = hero_uplink(elevation_deg)
    return replace(scn, ground=replace(scn.ground, compensation=list(stack)))


def campaign_of(elevation_deg, stack_name, args, plan_name=None,
                grid=None, plan=None, root=None):
    """Build (or reopen) the point-ahead campaign of one cell.

    Every campaign stores the field patch AND the summed screen phase of the
    beacon plane and of every point-ahead plane, so a post-hoc read answers
    any compensation stack with no new propagation.

    Args:
        elevation_deg: the elevation, in deg.
        stack_name:    a key of STACKS. "base" runs no in-run correction.
        args:          the parsed command line.
        plan_name:     the plan label of a p3 campaign, or None.
        grid:          a caller GridSpec, or None for the sizer.
        plan:          a caller ScreenPlan, or None for the sizer.
        root:          the store directory, or None for the default name.

    Returns:
        The pair (Campaign, the sizer warning texts).
    """
    scn, geom = hero_uplink(elevation_deg)
    stack = STACKS[stack_name]
    if root is None:
        tag = stack_name if plan_name is None else f"{plan_name}_{stack_name}"
        root = os.path.join(campaigns_root(), f"el{elevation_deg:02.0f}", tag)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        camp = Campaign(scn, geom, root, seed=SEED, preset=PRESET,
                        block_size=int(args.block_size),
                        patch_radius_m=PATCH_RADIUS_M, L0_m=L0_M,
                        precision="single", fft_backend=args.fft_backend,
                        compensation=(list(stack) or None),
                        store_screen_phase=True,
                        point_ahead_rad=angles_of(elevation_deg),
                        screen_margin_m=None, grid=grid, plan=plan)
    return camp, sorted({str(w.message) for w in caught})


def p3_campaigns(elevation_deg, args):
    """Give the pinned-grid campaigns of the p3 screen-count sweep.

    THE GRID IS PINNED. Each campaign takes the PRODUCTION grid and one
    override plan, so the sweep moves the screen count only. The recipe comes
    from `validation/terrestrial_screen_count/screen_count_sweep.py` through
    `validation.anisoplanatism_screens.common.plan_set`.

    THE n9 PLAN IS THE PRODUCTION PLAN, and it still gets its OWN store: the
    campaign fingerprint holds the repr of the caller grid and plan
    (`olb.waveoptics.turbulence.fingerprint.cache_key`), so a campaign that
    PASSES the production plan has a different key from one that lets the
    sizer build it. The function asserts the two plans agree.

    Args:
        elevation_deg: the elevation, in deg.
        args:          the parsed command line.

    Returns:
        A dict plan name -> (Campaign, the sizer warnings, the ScreenPlan).
    """
    scn, geom = hero_uplink(elevation_deg)
    grid, plans = plan_set(scn, geom, COUNTS, split_n=SPLIT_N)
    assert np.allclose(plans["n9"].z_m, plans["production"].z_m), \
        (plans["n9"].z_m, plans["production"].z_m)
    out = {}
    for name in P3_PLANS:
        camp, warns = campaign_of(elevation_deg, "base", args,
                                  plan_name=name, grid=grid,
                                  plan=plans[name])
        out[name] = (camp, warns, plans[name])
    return out


# ---------------------------------------------------------------------------
# The shared statistics
# ---------------------------------------------------------------------------

def pa_etas(camp, stack_name, n_trials, workers=None):
    """Give the point-ahead overlap of every trial and every angle.

    The read runs on the STORED planes, so it costs NO propagation. The
    sensing source follows the channel family, so a space campaign senses the
    stored summed screen phase. See `Campaign.recouple_point_ahead`.

    Args:
        camp:       the Campaign.
        stack_name: a key of STACKS. "base" applies no correction.
        n_trials:   the number of trials.
        workers:    None, an int or "auto", for the post-hoc pool.

    Returns:
        A float array of the shape (n_trials, n_angles).
    """
    stack = list(STACKS[stack_name])
    return camp.recouple_point_ahead(stack or None, n_trials=int(n_trials),
                                     workers=workers)


def _boot_rows(values, fn, n_boot=N_BOOT, seed=BOOT_SEED):
    """Give the bootstrap half-width of a statistic over the TRIAL rows.

    The function resamples the ROWS of a (n_trials, n_angles) array with
    replacement, so every column of one resample holds the same trials. That
    keeps the pairing between the beacon column and the point-ahead column,
    which a penalty difference needs. The half-width is the half of the
    central BOOT_INTERVAL band, and the seed is FIXED.

    Args:
        values: an array whose FIRST axis is the trial.
        fn:     a callable fn(resampled values) -> a float or a 1-D array.
        n_boot: the resample count.
        seed:   the fixed seed.

    Returns:
        A float or a 1-D array, the same shape as fn gives.
    """
    x = np.asarray(values, dtype=float)
    rng = np.random.default_rng(seed)
    draws = np.array([np.atleast_1d(fn(x[rng.integers(0, x.shape[0],
                                                      x.shape[0])]))
                      for _ in range(n_boot)])
    lo = np.quantile(draws, 0.5 - BOOT_INTERVAL / 2.0, axis=0)
    hi = np.quantile(draws, 0.5 + BOOT_INTERVAL / 2.0, axis=0)
    half = np.atleast_1d((hi - lo) / 2.0)
    return half if half.size > 1 else float(half[0])


def penalty_row(eta):
    """Give the mean loss and the anisoplanatic penalty of every angle.

    THE PENALTY is L(theta) - L(0), the extra loss that the point-ahead angle
    adds to the beacon direction. It is positive dB, because the two
    directions share less and less of the correction as the angle grows.

    Args:
        eta: the (n_trials, n_angles) overlap array.

    Returns:
        A dict of the loss, the penalty and their bootstrap half-widths.
    """
    eta = np.asarray(eta, dtype=float)

    def _loss(x):
        return -10.0 * np.log10(x.mean(axis=0))

    def _pen(x):
        v = _loss(x)
        return v - v[0]

    return {
        "n_trials": int(eta.shape[0]),
        "loss_db": [float(v) for v in _loss(eta)],
        "loss_half_db": [float(v) for v in
                         np.atleast_1d(_boot_rows(eta, _loss))],
        "penalty_db": [float(v) for v in _pen(eta)],
        "penalty_half_db": [float(v) for v in
                            np.atleast_1d(_boot_rows(eta, _pen))],
        "mean_eta": [float(v) for v in eta.mean(axis=0)],
    }


def fade_row(eta):
    """Give the per-trial fade of every angle.

    The per-trial loss is -10 log10(eta) of one trial, so the quantiles are
    the FADE of the link. pX is the loss the link EXCEEDS X percent of the
    time.

    Args:
        eta: the (n_trials, n_angles) overlap array.

    Returns:
        A dict with one `_fade_row` for each angle, and the p5 penalty.
    """
    eta = np.asarray(eta, dtype=float)
    loss = -10.0 * np.log10(eta)
    rows = [_fade_row(loss[:, i]) for i in range(loss.shape[1])]

    def _p5(x):
        return np.quantile(-10.0 * np.log10(x), 0.95, axis=0)

    def _dp5(x):
        v = _p5(x)
        return v - v[0]

    p5 = np.atleast_1d(_p5(eta))
    return {
        "angles": rows,
        "p5_db": [float(v) for v in p5],
        "p5_penalty_db": [float(v) for v in (p5 - p5[0])],
        "p5_penalty_half_db": [float(v) for v in
                               np.atleast_1d(_boot_rows(eta, _dp5))],
    }


# ---------------------------------------------------------------------------
# The other two rungs
# ---------------------------------------------------------------------------

def fast_penalties(elevation_deg, stack_name, angles, n_samples, say):
    """Measure the fidelity-1 FAST penalty of one stack at every angle.

    THE ANGLE SWEEP GOES THROUGH `fast_params`. `uplink_fast_term` builds
    DTHETA from `geometry.point_ahead_rad`, and it merges `fast_params` LAST,
    so `fast_params={"DTHETA": [arcsec, 0.0]}` sets the angle of the run. The
    geometry keeps its own elevation.

    THE GRID GUARD RUNS FIRST, at DTHETA = 0. FAST can undersample the
    low-order tilt, so the function takes the smallest NPXLS within
    NPXLS_TOL dB of the largest one. It copies the guard of
    `validation/waveoptics_ao/`.

    EVERY STACK HAS A FAST ROW. An empty stack is the NOAO launch and a
    tip-tilt stack is the TT launch (`olb.models.fast._ao_params`), so the
    `base` and the `tiptilt` rungs sit next to the AO rungs on the same grid.

    Source: O. J. D. Farley and others, DOI 10.1364/OE.458659.

    Args:
        elevation_deg: the elevation, in deg.
        stack_name:    a key of STACKS.
        angles:        the point-ahead angles, in rad.
        n_samples:     the FAST Monte Carlo draws.
        say:           the log function.

    Returns:
        A dict of the FAST measurement.
    """
    from olb.models.fast import uplink_fast_term

    scn = stacked_scenario(elevation_deg, STACKS[stack_name])
    geom = CircularOrbit(altitude_m=scn.channel.altitude_m,
                         elevation_deg=float(elevation_deg))

    def term_of(npxls, theta_rad):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return uplink_fast_term(
                scn, geom, n_samples=int(n_samples),
                fast_params={"NPXLS": int(npxls), "L0": float(L0_M),
                             "DTHETA": [float(theta_rad) / ARCSEC, 0.0]})

    guard = []
    for npxls in NPXLS_SET:
        t = term_of(npxls, 0.0)
        guard.append({"npxls": int(npxls), "mean_db": float(t.mean_db)})
    reference = guard[-1]["mean_db"]
    pinned = next((g for g in guard
                   if abs(g["mean_db"] - reference) <= NPXLS_TOL), guard[-1])
    say(f"    FAST grid guard ({stack_name}): "
        + ", ".join(f"{g['npxls']}:{g['mean_db']:.3f}" for g in guard)
        + f"  -> NPXLS {pinned['npxls']}")

    loss, quant, meta = [], [], {}
    for theta in angles:
        term = term_of(pinned["npxls"], theta)
        loss.append(float(term.mean_db))
        quant.append({f"p{int(p * 100)}": float(term.quantile_db(1.0 - p))
                      for p in EXCEEDANCE})
        meta = {"ao_mode": term.meta.get("ao_mode"),
                "zmax": term.meta.get("zmax")}
    return {
        "stack": stack_name,
        "elevation_deg": float(elevation_deg),
        "npxls": int(pinned["npxls"]),
        "guard": guard,
        "n_samples": int(n_samples),
        "loss_db": loss,
        "penalty_db": [v - loss[0] for v in loss],
        "quantiles_db": quant,
        **meta,
    }


def stone_penalties(elevation_deg, stack_name, angles):
    """Measure the analytic Stone penalty of one stack at every angle.

    The Term gives the anisoplanatic phase variance of the point-ahead angle
    over the corrected orders, and it maps it to a loss with the extended
    Marechal relation loss_db = (10 / ln 10) sigma^2 (T. S. Ross, Appl. Opt.
    48(10), 1812 (2009), DOI 10.1364/AO.48.001812). Source of the variance:
    Stone et al. (1994), DOI 10.1364/JOSAA.11.000347.

    THE MODE SET. `remove='piston'` keeps the tilt, the convention of this
    design. `max_order='auto'` reads the AO stage of the terminal; a TipTilt
    stack has no AO stage, so the caller passes `max_order=1` (the tilt is
    radial order 1).

    THE BASE STACK HAS NO STONE VALUE. With no corrected mode there is no
    decorrelation residual, so the anisoplanatic penalty is ZERO by
    definition and the caller skips it.

    Args:
        elevation_deg: the elevation, in deg.
        stack_name:    a key of STACKS, not "base".
        angles:        the point-ahead angles, in rad.

    Returns:
        A dict of the Stone measurement.
    """
    from olb.links.uplink import uplink_point_ahead_term

    scn = stacked_scenario(elevation_deg, STACKS[stack_name])
    max_order = 1 if stack_name == "tiptilt" else "auto"
    loss = []
    for theta in angles:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            term = uplink_point_ahead_term(
                scn, AngleGeometry(elevation_deg, theta), remove="piston",
                max_order=max_order, L0_m=L0_M)
        loss.append(float(np.atleast_1d(term.mean_db).ravel()[0]))
    return {
        "stack": stack_name,
        "elevation_deg": float(elevation_deg),
        "max_order": (None if max_order == "auto" else int(max_order)),
        "remove": "piston",
        "L0_m": float(L0_M),
        "loss_db": loss,
        "penalty_db": [v - loss[0] for v in loss],
    }


# ---------------------------------------------------------------------------
# The drivers
# ---------------------------------------------------------------------------

def _log_maker(name):
    """Make the log file of one study, and give its `say` function."""
    path = os.path.join(HERE, f"waveoptics_pointahead_{name}.log")
    with open(path, "w", encoding="utf-8"):
        pass

    def say(text=""):
        print(text, flush=True)
        with open(path, "a", encoding="utf-8") as fh:
            fh.write(text + "\n")

    return say, path


def _write_results(name, payload):
    """Write the JSON record of one study."""
    path = os.path.join(HERE, f"waveoptics_pointahead_{name}_results.json")
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=1)
    return path


def _fig_dir():
    """Give the figure directory of this study, and make it."""
    path = os.path.join(HERE, "figures")
    os.makedirs(path, exist_ok=True)
    return path


def _header(say, args, title):
    """Print the common header of a study."""
    say(title)
    say(f"date          : 2026-09-11      seed: {SEED}      preset: {PRESET}")
    say(f"outer scale   : L0 = {L0_M:g} m (the owner rule, backlog 2-P5)")
    say(f"precision     : single      fft backend: {args.fft_backend}      "
        f"workers: {args.workers}")
    say("scenario      : uplink, 1550 nm, 500 km, 700 mm ground terminal, "
        "0.35 m launch waist, 100 mm space terminal, DownlinkBeacon")
    say("mode set      : the TILT STAYS IN. The runner keeps it, FAST keeps "
        "it, and the Stone Term takes remove='piston'.")
    say("caveat        : PERFECT AO. It is an ideal modal fit of one "
        "snapshot: no sensor noise, no lag, no aliasing.")
    say()


def _campaign_line(name, camp, n_trials):
    """Give one line of the campaign list, and the stored byte estimate."""
    n_planes = 1 + len(camp.point_ahead_rad)
    pixels = int(camp.patch.indices.size)
    # complex64 field (8 B) plus float32 screen phase (4 B) for each plane.
    per_trial = pixels * n_planes * 12
    total = per_trial * int(n_trials)
    seconds = int(n_trials) * n_planes / GPU_TRIALS_PER_S
    text = (f"  {name:<22s}{int(n_trials):>7d}{n_planes:>8d}"
            f"{camp.grid.n:>7d}{camp.screen_n:>9d}"
            f"{camp.plan.z_m.size:>8d}{total / 2 ** 20:>10.0f}"
            f"{seconds / 60.0:>9.1f}")
    return text, total, seconds


def run_dry(args):
    """Print the campaign list of the asked studies, and compute nothing."""
    say, log_path = _log_maker("dry")
    _header(say, args, "THE CAMPAIGN LIST (a dry run: no trial is computed)")
    say(f"  {'campaign':<22s}{'trials':>7s}{'planes':>8s}{'grid':>7s}"
        f"{'screen':>9s}{'screens':>8s}{'MB':>10s}{'GPU min':>9s}")
    say("  " + "-" * 80)
    total_bytes, total_s = 0, 0.0
    seen = set()
    for el in args.elevations:
        for name in ("base", "ao10"):
            camp, _ = campaign_of(el, name, args)
            if camp.root_dir in seen:
                continue
            seen.add(camp.root_dir)
            # The ao10 campaign is the p0 reference only, so it runs to the
            # check count, not to the full trial count.
            wanted = args.n_trials if name == "base" else args.check_trials
            text, nbytes, seconds = _campaign_line(f"el{el:.0f}/{name}", camp,
                                                   wanted)
            say(text)
            total_bytes += nbytes
            total_s += seconds
    if "p3" in args.study:
        for el in [e for e in args.elevations if abs(e - 30.0) < 1e-6]:
            for plan_name, (camp, _w, _p) in p3_campaigns(el, args).items():
                if camp.root_dir in seen:
                    continue
                seen.add(camp.root_dir)
                text, nbytes, seconds = _campaign_line(
                    f"el{el:.0f}/{plan_name}", camp, args.count_trials)
                say(text)
                total_bytes += nbytes
                total_s += seconds
    say("  " + "-" * 80)
    say(f"  {'TOTAL':<22s}{'':>7s}{'':>8s}{'':>7s}{'':>9s}{'':>8s}"
        f"{total_bytes / 2 ** 20:>10.0f}{total_s / 60.0:>9.1f}")
    say()
    say("  planes = 1 + the angle count. Each plane is ONE split step and one "
        "stored field, so the cost scales with it.")
    say("  screen = the oversize screen side in px. The point-ahead window "
        "slides across it by whole pixels.")
    say(f"  GPU min uses the 2-AO rate of {GPU_TRIALS_PER_S:g} trials/s for "
        "ONE split step at 512 px and 9 screens.")
    say(f"wrote {log_path}")


def run_p0(args):
    """Run p0, the identities of the point-ahead record."""
    say, log_path = _log_maker("p0")
    _header(say, args, "p0: the identities of the point-ahead record "
                       "(backlog 2-P4)")
    say("THREE CHECKS, on each elevation:")
    say("  1. the stored angle 0.0 equals the beacon overlap eta_turb, bit "
        "for bit.")
    say("  2. the post-hoc read of the UNCORRECTED planes equals an IN-RUN "
        "AO(10) campaign of the same seeds.")
    say("  3. the regeneration route at a STORED angle gives the stored "
        "column back.")
    say("  Check 3 runs the regeneration on the backend of the campaign, so "
        "the gate is BIT identity. The host (numpy)")
    say("  read of a CUDA campaign is printed next to it as a report: it "
        "agrees at the float32 rounding level only.")
    say()

    rows = {}
    for el in args.elevations:
        say(f"ELEVATION {el:.0f} deg")
        angles = angles_of(el)
        base, warns = campaign_of(el, "base", args)
        for w in warns:
            say(f"  sizer warning: {w}")
        say(f"  grid {base.grid.n} px, {base.grid.pixel_m * 1e3:.2f} mm "
            f"pixel, screen {base.screen_n} px, margin "
            f"{base.screen_margin_m:.3f} m, {base.plan.z_m.size} screens")
        say(f"  angles: " + ", ".join(
            f"{lab} = {a * 1e6:.3f} urad" for lab, a
            in zip(angle_labels(angles), angles)))
        ensure_trials(base, args.n_trials, args, say)
        ao10, _ = campaign_of(el, "ao10", args)
        ensure_trials(ao10, args.check_trials, args, say)

        n_check = int(min(args.check_trials, base.n_stored, ao10.n_stored))
        result = base.load(args.n_trials, fields=False)
        stored = np.array([t.eta_turb_pa for t in result.trials], dtype=float)
        beacon = np.array([t.eta_turb for t in result.trials], dtype=float)
        walls = np.array([t.wall_time_s for t in result.trials], dtype=float)
        same = bool(np.array_equal(stored[:, 0], beacon))

        # Check 2: the stored UNCORRECTED planes, corrected post hoc.
        post = base.recouple_point_ahead([AO(n_modes=10)], n_trials=n_check,
                                         workers=args.posthoc_workers)
        inrun = np.array([t.eta_turb_pa for t in
                          ao10.load(n_check, fields=False).trials], dtype=float)
        rel_post = float(np.abs(post / inrun - 1.0).max())

        # Check 3: the regeneration route at the STORED 5 arcsec angle.
        i5 = len(angles) - 1 - FIXED_ARCSEC[::-1].index(5.0)
        # The SAME backend as the campaign gives the stored column back bit
        # for bit. The host read is a report only (the float32 rounding of
        # two FFT libraries).
        regen = base.point_ahead([angles[i5]], None, n_trials=n_check,
                                 fft_backend=args.fft_backend)
        rel_regen = float(np.abs(regen[:, 0] / stored[:n_check, i5]
                                 - 1.0).max())
        bit_regen = bool(np.array_equal(regen[:, 0], stored[:n_check, i5]))
        rel_host = None
        if args.fft_backend != "numpy":
            host = base.point_ahead([angles[i5]], None, n_trials=n_check,
                                    fft_backend="numpy")
            rel_host = float(np.abs(host[:, 0] / stored[:n_check, i5]
                                    - 1.0).max())

        rows[str(el)] = {
            "elevation_deg": float(el),
            "angles_rad": [float(a) for a in angles],
            "n_trials": int(stored.shape[0]),
            "n_check": n_check,
            "grid_n": int(base.grid.n),
            "screen_n": int(base.screen_n),
            "screen_margin_m": float(base.screen_margin_m),
            "n_screens": int(base.plan.z_m.size),
            "beacon_identity": same,
            "posthoc_max_rel": rel_post,
            "regenerate_max_rel": rel_regen,
            "regenerate_bit_identical": bit_regen,
            "regenerate_host_max_rel": rel_host,
            "wall_s_per_trial": float(walls.mean()),
            "fft_backend": args.fft_backend,
            "disk_bytes": _dir_bytes(base.root_dir),
        }
        say(f"  1. the beacon identity eta_turb_pa[:, 0] == eta_turb: "
            f"{'PASS (bit for bit)' if same else 'FAIL'}")
        say(f"  2. the post-hoc AO(10) against the in-run AO(10), {n_check} "
            f"trials: worst relative error {rel_post:.2e} "
            f"({'PASS' if rel_post < 1e-6 else 'FAIL'}, the gate is 1e-6)")
        say(f"  3. the regeneration at 5 arcsec on the {args.fft_backend} "
            f"backend, {n_check} trials: worst relative error "
            f"{rel_regen:.2e} ({'PASS, bit identical' if bit_regen else 'FAIL'}"
            f", the gate is bit identity)"
            + (f"; the host read differs by {rel_host:.2e}"
               if rel_host is not None else ""))
        say(f"  the trial wall time is {walls.mean():.3f} s for "
            f"{1 + len(angles)} split steps "
            f"({walls.mean() / (1 + len(angles)):.3f} s for each pass)")
        say(f"  the screen side is {base.screen_n} px against a {base.grid.n} "
            f"px propagation grid "
            f"({base.screen_n / base.grid.n:.2f}x)")
        say()

    path = _write_results("p0", {"study": "p0", "date": "2026-09-11",
                                 "rows": rows})
    say(f"wrote {path}")
    say(f"wrote {log_path}")


def run_p1(args):
    """Run p1, the penalty against FAST and against Stone."""
    say, log_path = _log_maker("p1")
    _header(say, args, "p1: the point-ahead penalty against FAST and Stone "
                       "(backlog 2-P4)")
    say("THE FIELD PENALTY is L(theta) - L(0) with L = -10 log10(mean "
        "eta_pa). Loss is positive dB.")
    say("THE FAST PENALTY is the same difference of the fidelity-1 Term, at "
        "the same stack and the same outer scale.")
    say("THE STONE PENALTY is the analytic (10 / ln10) sigma^2 of the "
        "corrected orders. It SATURATES past sigma^2 = 1 rad^2 (the extended")
    say("Marechal limit), so read it as a report, not as a gate.")
    say(f"THE PASS BAND of the field against FAST is {FAST_PASS_DB:g} dB, "
        "the like-for-like tolerance of docs/physics.md Section 9l.")
    say()

    out = {"study": "p1", "date": "2026-09-11", "elevations": {},
           "pass_band_db": FAST_PASS_DB}
    for el in args.elevations:
        say(f"ELEVATION {el:.0f} deg")
        angles = angles_of(el)
        labels = angle_labels(angles)
        base, warns = campaign_of(el, "base", args)
        for w in warns:
            say(f"  sizer warning: {w}")
        ensure_trials(base, args.n_trials, args, say)
        n = int(min(args.n_trials, base.n_stored))

        field = {}
        for name in STACK_ORDER:
            field[name] = penalty_row(pa_etas(base, name, n,
                                              args.posthoc_workers))
        say()
        say(f"  THE FIELD LOSS -10 log10(mean eta_pa) [dB], {n} trials")
        hdr = ("  " + f"{'stack':<9s}" +
               "".join(f"{lab:>12s}" for lab in labels))
        say(hdr)
        say("  " + "-" * (len(hdr) - 2))
        for name in STACK_ORDER:
            say(f"  {name:<9s}" + "".join(
                f"{v:>12.3f}" for v in field[name]["loss_db"]))
        say("  The beacon column IS the 2-AO perfect pre-compensation loss. "
            "The point-ahead columns degrade it.")
        say()
        say(f"  THE ANISOPLANATIC PENALTY L(theta) - L(0) [dB]")
        say(hdr)
        say("  " + "-" * (len(hdr) - 2))
        for name in STACK_ORDER:
            say(f"  {name:<9s}" + "".join(
                f"{v:>12.3f}" for v in field[name]["penalty_db"]))
        say()

        fast_rows, stone_rows, compare = {}, {}, []
        for name in STACK_ORDER:
            if name != "base":
                stone_rows[name] = stone_penalties(el, name, angles)
            if args.no_fast:
                continue
            try:
                fast_rows[name] = fast_penalties(el, name, angles,
                                                 args.fast_samples, say)
            except ImportError as exc:
                say(f"    fast-aosim is not available: {exc}")
                args.no_fast = True
            except Exception as exc:            # noqa: BLE001
                say(f"    FAST failed for {name}: "
                    f"{type(exc).__name__}: {exc}")

        say()
        say("  THE THREE RUNGS, on the anisoplanatic penalty [dB]")
        say(f"  {'stack':<9s}{'angle':>12s}{'field':>9s}{'+-':>7s}"
            f"{'FAST':>9s}{'gap':>8s}{'band':>7s}{'Stone':>9s}")
        say("  " + "-" * 68)
        for name in STACK_ORDER:
            for i, lab in enumerate(labels):
                if i == 0:
                    continue
                fld = field[name]["penalty_db"][i]
                half = field[name]["penalty_half_db"][i]
                f_row = fast_rows.get(name)
                s_row = stone_rows.get(name)
                fast_v = None if f_row is None else f_row["penalty_db"][i]
                stone_v = None if s_row is None else s_row["penalty_db"][i]
                gap = None if fast_v is None else fast_v - fld
                verdict = ("" if gap is None
                           else ("PASS" if abs(gap) <= FAST_PASS_DB
                                 else "OVER"))
                say(f"  {name:<9s}{lab:>12s}{fld:>9.3f}{half:>7.3f}"
                    + (f"{fast_v:>9.3f}{gap:>8.3f}{verdict:>7s}"
                       if fast_v is not None else f"{'-':>9s}{'-':>8s}"
                                                  f"{'-':>7s}")
                    + (f"{stone_v:>9.3f}" if stone_v is not None
                       else f"{'-':>9s}"))
                compare.append({
                    "stack": name, "angle_label": lab,
                    "angle_rad": float(angles[i]),
                    "field_db": fld, "field_half_db": half,
                    "fast_db": fast_v, "gap_db": gap, "verdict": verdict,
                    "stone_db": stone_v})
        say("  gap = FAST minus the field. The band column reads PASS while "
            f"|gap| <= {FAST_PASS_DB:g} dB.")
        say("  Stone has no `base` row: with no corrected mode there is no "
            "decorrelation residual, so the penalty is zero by definition.")
        say()

        out["elevations"][str(el)] = {
            "angles_rad": [float(a) for a in angles],
            "angle_labels": labels, "n_trials": n,
            "field": field, "fast": fast_rows, "stone": stone_rows,
            "compare": compare}

    figs = plot_p1(out) if args.figures else []
    out["figures"] = figs
    path = _write_results("p1", out)
    say(f"wrote {path}")
    say(f"wrote {log_path}")
    for f in figs:
        say(f"wrote {f}")


def run_p2(args):
    """Run p2, the fade of the point-ahead uplink."""
    say, log_path = _log_maker("p2")
    _header(say, args, "p2: the point-ahead FADE (backlog 2-P4)")
    say("pX is the loss the link EXCEEDS X percent of the time. A LARGE loss "
        "is a deep fade.")
    say("The p5 penalty is p5(theta) - p5(0): the extra deep-fade loss of the "
        "point-ahead angle.")
    say("The +- bar is the 68 percent bootstrap half-width over the trial "
        "rows; the 2-sigma gate is twice it.")
    say()

    out = {"study": "p2", "date": "2026-09-11", "elevations": {}}
    for el in args.elevations:
        say(f"ELEVATION {el:.0f} deg")
        angles = angles_of(el)
        labels = angle_labels(angles)
        base, warns = campaign_of(el, "base", args)
        for w in warns:
            say(f"  sizer warning: {w}")
        ensure_trials(base, args.n_trials, args, say)
        n = int(min(args.n_trials, base.n_stored))

        rows = {}
        for name in STACK_ORDER:
            rows[name] = fade_row(pa_etas(base, name, n,
                                          args.posthoc_workers))
        say()
        say(f"  THE PER-TRIAL LOSS -10 log10(eta_pa) [dB], {n} trials")
        say(f"  {'stack':<9s}{'angle':>12s}{'mean':>8s}{'p50':>8s}"
            f"{'p10':>8s}{'p5':>8s}{'p1':>8s}{'d p5':>8s}{'+-':>7s}")
        say("  " + "-" * 68)
        for name in STACK_ORDER:
            for i, lab in enumerate(labels):
                fade = rows[name]["angles"][i]
                q = fade["quantiles_db"]
                say(f"  {name:<9s}{lab:>12s}{fade['mean_db']:>8.2f}"
                    f"{q['p50']:>8.2f}{q['p10']:>8.2f}{q['p5']:>8.2f}"
                    f"{q['p1']:>8.2f}"
                    f"{rows[name]['p5_penalty_db'][i]:>8.2f}"
                    f"{rows[name]['p5_penalty_half_db'][i]:>7.2f}")
            say()
        say("  d p5 = p5(theta) - p5(beacon), the p5 anisoplanatic penalty.")
        say()
        out["elevations"][str(el)] = {
            "angles_rad": [float(a) for a in angles],
            "angle_labels": labels, "n_trials": n, "stacks": rows}

    figs = plot_p2(out) if args.figures else []
    out["figures"] = figs
    path = _write_results("p2", out)
    say(f"wrote {path}")
    say(f"wrote {log_path}")
    for f in figs:
        say(f"wrote {f}")


def run_p3(args):
    """Run p3, the screen-count sweep of the point-ahead penalty."""
    say, log_path = _log_maker("p3")
    _header(say, args, "p3: the SCREEN COUNT of the point-ahead penalty "
                       "(backlog 2-P4)")
    say("THE GRID IS PINNED at the production grid, so the sweep moves the "
        "screen count only.")
    say("THE PASS RULE (the 9i criterion): the penalty is FLAT from 9 screens "
        "up, inside the 2-sigma bootstrap band.")
    say("The n9 plan IS the production plan, and it still gets its own store: "
        "the campaign fingerprint reads the repr of the")
    say("caller grid and plan, so a PASSED production plan keys differently "
        "from a SIZED one. The script asserts the two plans agree.")
    say()

    elevations = [e for e in args.elevations if abs(e - 30.0) < 1e-6] or [30.0]
    out = {"study": "p3", "date": "2026-09-11", "elevations": {}}
    for el in elevations:
        say(f"ELEVATION {el:.0f} deg")
        angles = angles_of(el)
        labels = angle_labels(angles)
        # The reported angles: the 5 arcsec and the 10 arcsec columns.
        want = [len(angles) - 1 - FIXED_ARCSEC[::-1].index(a)
                for a in (5.0, 10.0)]
        cells = p3_campaigns(el, args)
        rows = {}
        for plan_name, (camp, warns, plan) in cells.items():
            for w in warns:
                say(f"  {plan_name}: sizer warning: {w}")
            say(f"  {plan_name}: {plan.z_m.size} screens, r0_total "
                f"{plan.r0_total_m * 100:.2f} cm, grid {camp.grid.n} px, "
                f"screen {camp.screen_n} px")
            ensure_trials(camp, args.count_trials, args, say)
            n = int(min(args.count_trials, camp.n_stored))
            per_stack = {}
            for name in STACK_ORDER:
                eta = pa_etas(camp, name, n, args.posthoc_workers)
                per_stack[name] = {"mean": penalty_row(eta),
                                   "fade": fade_row(eta)}
            rows[plan_name] = {"n_screens": int(plan.z_m.size),
                               "n_trials": n,
                               "r0_total_m": float(plan.r0_total_m),
                               "stacks": per_stack}
        say()

        for i in want:
            say(f"  THE PENALTY AT {labels[i]} [dB], "
                f"{rows[P3_PLANS[0]]['n_trials']} trials for each plan")
            say(f"  {'plan':<18s}{'screens':>8s}" + "".join(
                f"{name:>18s}" for name in STACK_ORDER))
            say("  " + "-" * (26 + 18 * len(STACK_ORDER)))
            for plan_name in P3_PLANS:
                r = rows[plan_name]
                line = f"  {plan_name:<18s}{r['n_screens']:>8d}"
                for name in STACK_ORDER:
                    m = r["stacks"][name]["mean"]
                    line += (f"{m['penalty_db'][i]:>11.3f}"
                             f"+-{m['penalty_half_db'][i]:<5.3f}")
                say(line)
            say()
            say(f"  THE p5 PENALTY AT {labels[i]} [dB]")
            say(f"  {'plan':<18s}{'screens':>8s}" + "".join(
                f"{name:>18s}" for name in STACK_ORDER))
            say("  " + "-" * (26 + 18 * len(STACK_ORDER)))
            for plan_name in P3_PLANS:
                r = rows[plan_name]
                line = f"  {plan_name:<18s}{r['n_screens']:>8d}"
                for name in STACK_ORDER:
                    f = r["stacks"][name]["fade"]
                    line += (f"{f['p5_penalty_db'][i]:>11.3f}"
                             f"+-{f['p5_penalty_half_db'][i]:<5.3f}")
                say(line)
            say()

        # The 9i pass rule, against the 25-screen plan.
        verdicts = []
        for i in want:
            for name in STACK_ORDER:
                ref = rows["n25"]["stacks"][name]["mean"]
                for plan_name in ("n9", "n15", "n25"):
                    m = rows[plan_name]["stacks"][name]["mean"]
                    band = 2.0 * (m["penalty_half_db"][i]
                                  + ref["penalty_half_db"][i])
                    diff = m["penalty_db"][i] - ref["penalty_db"][i]
                    verdicts.append({
                        "angle_label": labels[i], "stack": name,
                        "plan": plan_name, "diff_db": float(diff),
                        "band_db": float(band),
                        "flat": bool(abs(diff) <= band)})
        n_flat = sum(1 for v in verdicts if v["flat"])
        say(f"  THE 9i PASS RULE against the 25-screen plan: "
            f"{n_flat} of {len(verdicts)} readings are flat inside the "
            f"2-sigma band.")
        for v in verdicts:
            if not v["flat"]:
                say(f"    NOT FLAT: {v['plan']}, {v['stack']}, "
                    f"{v['angle_label']}: {v['diff_db']:+.3f} dB against a "
                    f"{v['band_db']:.3f} dB band.")
        say()
        out["elevations"][str(el)] = {
            "angles_rad": [float(a) for a in angles],
            "angle_labels": labels, "reported": want,
            "plans": rows, "verdicts": verdicts}

    path = _write_results("p3", out)
    say(f"wrote {path}")
    say(f"wrote {log_path}")


# ---------------------------------------------------------------------------
# The figures
# ---------------------------------------------------------------------------

def plot_p1(out):
    """Draw the mean penalty against the angle, one panel for each elevation."""
    plt = _pyplot()
    if plt is None:
        return []
    els = sorted(out["elevations"], key=float)
    fig, axes = plt.subplots(1, len(els), figsize=(6.0 * len(els), 5.0),
                             squeeze=False)
    for ax, key in zip(axes[0], els):
        cell = out["elevations"][key]
        x = np.asarray(cell["angles_rad"], dtype=float) / ARCSEC
        for name in STACK_ORDER:
            ax.plot(x, cell["field"][name]["penalty_db"], "o-",
                    label=f"field {name}")
            if name in cell["fast"]:
                ax.plot(x, cell["fast"][name]["penalty_db"], "s--",
                        label=f"FAST {name}")
            if name in cell["stone"]:
                ax.plot(x, cell["stone"][name]["penalty_db"], ":",
                        label=f"Stone {name}")
        ax.set_xlabel("point-ahead angle [arcsec]")
        ax.set_ylabel("anisoplanatic penalty [dB]")
        ax.set_title(f"{float(key):.0f} deg")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7)
    fig.tight_layout()
    path = os.path.join(_fig_dir(), "p1_penalty_vs_angle.png")
    fig.savefig(path, dpi=140)
    plt.close(fig)
    return [path]


def plot_p2(out):
    """Draw the p5 fade against the angle, one panel for each elevation."""
    plt = _pyplot()
    if plt is None:
        return []
    els = sorted(out["elevations"], key=float)
    fig, axes = plt.subplots(1, len(els), figsize=(6.0 * len(els), 5.0),
                             squeeze=False)
    for ax, key in zip(axes[0], els):
        cell = out["elevations"][key]
        x = np.asarray(cell["angles_rad"], dtype=float) / ARCSEC
        for name in STACK_ORDER:
            ax.plot(x, cell["stacks"][name]["p5_db"], "o-", label=name)
        ax.set_xlabel("point-ahead angle [arcsec]")
        ax.set_ylabel("p5 loss [dB]")
        ax.set_title(f"{float(key):.0f} deg")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
    fig.tight_layout()
    path = os.path.join(_fig_dir(), "p2_p5_vs_angle.png")
    fig.savefig(path, dpi=140)
    plt.close(fig)
    return [path]


STUDIES = {"p0": run_p0, "p1": run_p1, "p2": run_p2, "p3": run_p3}


def main():
    global PRESET
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--study", nargs="+", default=["p0", "p1", "p2", "p3"],
                    help="the studies to run, from p0 p1 p2 p3")
    ap.add_argument("--elevations", nargs="+", type=float,
                    default=[30.0, 20.0],
                    help="the elevations of the campaigns [deg]")
    ap.add_argument("--n-trials", type=int, default=1000,
                    help="the trials of each p0, p1 and p2 campaign")
    ap.add_argument("--count-trials", type=int, default=400,
                    help="the trials of each p3 screen-count campaign")
    ap.add_argument("--check-trials", type=int, default=200,
                    help="the trials of the p0 identity checks")
    ap.add_argument("--block-size", type=int, default=50,
                    help="the trials in one block file")
    ap.add_argument("--workers", default=None,
                    help="the campaign pool size. Leave it out for the CUDA "
                         "backend: one device runs one stream.")
    ap.add_argument("--posthoc-workers", default=None,
                    help="the pool size of the post-hoc reads. None runs "
                         "them in this process.")
    ap.add_argument("--fft-backend", default="cupy",
                    choices=["numpy", "scipy", "cupy"],
                    help="the FFT backend of the campaigns")
    ap.add_argument("--preset", default=PRESET,
                    help="the sampling preset of the campaigns")
    ap.add_argument("--analyse-only", action="store_true",
                    help="read what is stored and compute no trial")
    ap.add_argument("--dry-run", action="store_true",
                    help="print the campaign list and exit")
    ap.add_argument("--no-figures", dest="figures", action="store_false",
                    help="skip the figures")
    ap.add_argument("--fast-samples", type=int, default=1000,
                    help="the FAST Monte Carlo draws of the p1 comparison")
    ap.add_argument("--no-fast", action="store_true",
                    help="skip the fidelity-1 FAST comparison")
    args = ap.parse_args()
    for name in ("workers", "posthoc_workers"):
        value = getattr(args, name)
        if value not in (None, "auto"):
            setattr(args, name, int(value))
    PRESET = args.preset

    for name in args.study:
        if name not in STUDIES:
            raise SystemExit(f"unknown study {name!r}. Use one of "
                             f"{sorted(STUDIES)}.")
    if args.dry_run:
        run_dry(args)
        return
    for name in args.study:
        STUDIES[name](args)


if __name__ == '__main__':
    main()
