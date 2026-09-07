"""The terrestrial SCREEN-COUNT convergence sweep (backlog 2-TC1).

THE QUESTION. The turbulent planner takes its screen count from the Schmidt
per-screen cap: N = max(min_screens, ceil(sigma_R^2 / 0.4)). At the 10 km,
Cn2 = 1e-14 standard cell that cap asks for 35 screens. The cap is a per-screen
THIN-SCREEN VALIDITY rule (Schmidt, Numerical Simulation of Optical Wave
Propagation (2010), DOI 10.1117/3.866274, Listing 9.5, lines 37 and 38, printed
p. 175, which the book credits to Martin and Flatte, DOI 10.1364/AO.27.002111).
It is NOT a convergence result, and olb has no terrestrial convergence sweep. So
the count that a uniform horizontal path really needs is UNMEASURED. This study
measures it.

THE CELL IS A CHOICE. `--path-km` and `--cn2` select any backbone cell. The
default is the 10 km, Cn2 = 1e-14 cell, where the cap binds at 35 screens. The
second cell of record is 5 km at Cn2 = 3e-15, which sits at the `min_screens`
floor of 9 screens, so it asks the OTHER half of the question: is the FLOOR
enough? The reference screen count is never hard-coded: the script asks the
sizer for the count of the UNMODIFIED standard preset of the selected cell.

THE METHOD. The study PINS the grid and it moves the screens only. Each case is
one `Campaign` with a CALLER plan: the sweep asks `turbulent_grid` for the same
equal-Rytov-weight cut of `_plan_terrestrial` with the cap turned off
(`sigma2_r_screen_max = 1e9`) and the floor set to the wanted count
(`min_screens = n`), so the plan holds EXACTLY n screens of the production
shape. The grid comes from the UNMODIFIED standard sizing, and the script
asserts that the override does not move it.

THE REFERENCE is the standard-preset campaign of the terrestrial backbone run
of the SAME cell, for example
`validation/terrestrial_campaigns/campaigns/L10km_cn21e-14_standard`. This
script NEVER runs it: it reopens that store with the same `Campaign` settings
that `run_campaigns.make_campaign` built it with, and it reads the FIRST
`--n-trials` trials. Those trials are bit-identical to a shorter run of the same
seed, because the runner seeds trial k off (entropy, k). See
`olb.waveoptics.turbulence.campaign`.

THE REFERENCE-COUNT BIT-IDENTITY CHECK (on by default, `--no-ref-check` skips
it) PROVES the pipeline before any conclusion. It builds the override campaign
at n = the reference screen count, into its own root, it runs
`--ref-check-trials` trials there, and it asserts that the `collected_power`
and the `smf_eta` of those trials equal the reference values EXACTLY. It first
asserts that the override plan equals the production plan (`z_m` and `r0_m`),
because the `min_screens = n` route and the production cap route could in
principle place the screens differently. So a PASS says that the grid pinning,
the plan override, the speed opt-ins and the campaign reopen all reproduce the
production path.

THE CELL is exactly a backbone cell: a horizontal path of `--path-km` at
`--cn2`, a 5 mm waist from a 100 mm terminal, a 100 mm receive aperture with a
single-mode fibre, the standard preset, L0 = 25 m, single precision, seed
20260906. The cell constants and the scenario builder are IMPORTED from
`validation.terrestrial_campaigns.run_campaigns`, so the two studies cannot
drift apart.

THE MEASURED QUANTITIES. Each case reports four quantities, and each one gives
a scintillation index and the fade quantiles p10, p5 and p1 of
-10 log10(x / <x>), with a 1000-resample bootstrap 68 percent interval:
  P 10 cm     the collected power of the 100 mm receive bucket.
  point       the irradiance of the CENTRE PIXEL of the receive field. Nothing
              averages a point, so it is the sharpest probe of the screen plan.
              The dry-run table gives the pixel of the grid.
  smf_eta     the single-mode-fibre coupling efficiency of that aperture.
  P 5 cm      the collected power of a 50 mm bucket, a post-hoc `recollect` of
              the SAME stored fields.
The index is the normalised variance var(x)/mean(x)^2. Source: Andrews and
Phillips, Laser Beam Propagation through Random Media, 2nd ed. (2005),
DOI 10.1117/3.626196, Ch. 8. olb holds the definition in
`olb.turbulence.andrews.scintillation`.

THE TWO VERDICTS. Each count gets both, because they answer different
questions.
  RESOLUTION (inside the reference noise). A count is CONVERGED when, for every
  one of those four quantities, the p5 and the p1 delta against the reference
  sits inside the bootstrap half-width of the reference, AND the index ratio
  sits inside 1 +/- 0.10. Else it is NOT CONVERGED, and the line names the
  quantity that fails. This says whether the run can SEE a difference, not
  whether the difference matters.
  TOLERANCE (`--tolerance-db`, default 1.0 dB). A count PASSES the tolerance
  when every |p5| and |p1| delta of the 10 cm bucket AND of the SMF coupling is
  under that many dB. This says whether the difference matters to a link
  budget.
Each count is a Monte Carlo estimate too, so the noise on a DELTA is about
sqrt(2) times the reference bootstrap half-width. The reference bar table
prints both numbers.

NO CORRECTION. Fidelity 2 models no tip-tilt removal and no adaptive optics
(backlog 2-AO), so every fibre number here holds the full beam wander.

VALIDATION ONLY. The script reads the production layer. It changes no olb
module.

Sources:
- Schmidt, DOI 10.1117/3.866274, Ch. 9. The split-step method, the per-screen
  cap of Listing 9.5, printed p. 175, and the layer moment rule Eq. (9.65),
  printed p. 164.
- Andrews and Phillips, DOI 10.1117/3.626196. Ch. 8, Eq. (20), printed p. 264:
  the plane-wave Rytov variance and its path weight (the equal-weight cut).
- Martin and Flatte, DOI 10.1364/AO.27.002111. The per-screen strength cap and
  the pixel-per-coherence-length rule.
- Fried, DOI 10.1364/JOSA.56.001372. The Fried parameter r0.
- The analysis idioms (the quantile bootstrap, the centre pixel, the figure
  layout) come from validation/tail_convergence/tail_convergence.py and
  validation/outer_scale_tail/outer_scale_tail.py.

Run it from the repository root:

    python -m validation.terrestrial_screen_count.screen_count_sweep --dry-run
    python -m validation.terrestrial_screen_count.screen_count_sweep --workers 12
    python -m validation.terrestrial_screen_count.screen_count_sweep --analyse-only
    python -m validation.terrestrial_screen_count.screen_count_sweep \
        --fft-backend scipy --screen-generator olb-lean
    python -m validation.terrestrial_screen_count.screen_count_sweep \
        --path-km 5 --cn2 3e-15 --n-trials 2000 \
        --fft-backend scipy --screen-generator olb-lean

THE FILE NAMES carry the cell tag, for example
`screen_count_sweep_L5km_cn23e-15_results.json`. The DEFAULT cell (10 km,
Cn2 = 1e-14) keeps the plain names `screen_count_sweep_results.json`,
`screen_count_sweep.log` and `figures/screen_count_sweep.png`, so the stored
record of that cell stays where it is.

THE TWO SPEED OPT-INS. `--fft-backend scipy` and `--screen-generator olb-lean`
go to BOTH the reference reopen and every override campaign, and each one
ENTERS the campaign fingerprint. So they give every root the suffix `_scipy`,
`_lean` or `_scipy_lean`, and the reference they name is the opt-in backbone
cell, for example `L10km_cn21e-14_standard_scipy_lean`. A reopen with the wrong
settings RAISES a fingerprint mismatch. See
validation/terrestrial_campaigns/README.md.
"""

import argparse
import dataclasses
import json
import os
import time
import warnings

import numpy as np

from olb.waveoptics.turbulence.campaign import Campaign
from olb.waveoptics.turbulence.sampling import PRESETS, turbulent_grid
from validation.terrestrial_campaigns import run_campaigns as backbone

HERE = os.path.dirname(os.path.abspath(__file__))

# The environment override of the campaign store of THIS study. The reference
# store has its own override, the one that run_campaigns uses.
ENV_ROOT = "OLB_SCREEN_COUNT_ROOT"

# The DEFAULT cell. It is a backbone cell, and the constants come from that
# module. `--path-km` and `--cn2` select another one.
PATH_M = 10e3
CN2 = 1e-14
PRESET = "standard"
LAUNCH = "collimated"

# The screen counts of the sweep, cheapest first.
DEFAULT_COUNTS = (5, 10, 15, 20)

# The default trials of each count, and of the reference slice.
DEFAULT_TRIALS = 2000

# The default trials of the reference-count bit-identity check.
REF_CHECK_TRIALS = 100

# The default dB tolerance of the second verdict. See the module docstring.
TOLERANCE_DB = 1.0

# The quantities that the dB tolerance verdict reads. They are the two
# RECEIVER kinds of a real terrestrial link: a bucket and a fibre.
TOLERANCE_QUANTITIES = ("P10cm", "smf_eta")

# The cost model of the dry run. Each anchor is (s/trial, screens, workers) of
# a MEASURED backbone run on bigfraw, and a trial is close to linear in the
# screen count, so the projection scales by n / screens and by workers / the
# asked workers.
#   10 km, 1e-14: the backbone smoke run, about 6.0 s of WALL time for one
#                 trial at 35 screens on an 8-worker pool (the pool time, not
#                 the time of one process).
#   5 km, 3e-15:  the backbone full run, 1.12 s for one trial at 9 screens on
#                 a 12-worker pool with the two speed opt-ins.
# A cell that is not in the table takes the 5 km anchor, the nearest
# measurement, and the projection is then a rough guide only.
COST_ANCHORS = {
    (10e3, 1e-14): (6.0, 35, 8),
    (5e3, 3e-15): (1.12, 9, 12),
}
FALLBACK_ANCHOR = (1.12, 9, 12)

# The receive diameters that the analysis reads. The 100 mm value is the
# propagated aperture; the 50 mm value is a post-hoc crop of the stored field.
APERTURE_M = 0.10
SMALL_APERTURE_M = 0.05

# The exceedance probabilities. pX is the fade that the link EXCEEDS X percent
# of the time, so it is the (100 - X) percentile of the fade.
EXCEEDANCE = (0.10, 0.05, 0.01)

# The bootstrap of every quantile and every index.
N_BOOT = 1000
BOOT_INTERVAL = 0.68
BOOT_SEED = 7

# The convergence gate on the index ratio. See the module docstring.
INDEX_TOLERANCE = 0.10

# The four measured quantities, in the print order.
QUANTITIES = ("P10cm", "point", "smf_eta", "P5cm")
QUANTITY_DOC = {
    "P10cm": "the collected power of the 100 mm receive bucket",
    "point": "the irradiance of the CENTRE PIXEL of the receive grid",
    "smf_eta": "the single-mode-fibre coupling efficiency",
    "P5cm": "the collected power of a 50 mm bucket (a post-hoc recollect)",
}


# ---------------------------------------------------------------------------
# The campaigns
# ---------------------------------------------------------------------------

def campaigns_root():
    """Give the parent directory of the override campaigns of this study."""
    return os.environ.get(ENV_ROOT) or os.path.join(HERE, "campaigns")


def case_tag(n_screens, path_m=PATH_M, cn2=CN2,
             fft_backend=backbone.FFT_BACKEND,
             screen_generator=backbone.SCREEN_GENERATOR):
    """Give the directory name of one override campaign.

    The name carries the CELL tag of the backbone, then the screen count, then
    the opt-in suffix, for example `L5km_cn23e-15_standard_n5_scipy_lean`. The
    two speed opt-ins enter the campaign fingerprint, so they give the root the
    same suffix that the backbone gives its own roots.

    Args:
        n_screens:        the screen count.
        path_m:           the horizontal path length, in m.
        cn2:              the Cn2 of the path, in m^-2/3.
        fft_backend:      "numpy" or the "scipy" opt-in.
        screen_generator: "olb" or the "olb-lean" opt-in.

    Returns:
        The directory name.
    """
    return (backbone.cell_tag(path_m, cn2, PRESET, LAUNCH, False)
            + f"_n{int(n_screens)}"
            + backbone.settings_suffix(fft_backend, screen_generator))


def file_stem(path_m=PATH_M, cn2=CN2):
    """Give the cell part of an output file name.

    The DEFAULT cell gives an EMPTY string, so the log, the JSON record and
    the figure of the 10 km cell keep the names they already have on disk.
    Every other cell gives `_L<path>km_cn2<cn2>`, the same cell text that the
    backbone `cell_tag` builds.

    Args:
        path_m: the horizontal path length, in m.
        cn2:    the Cn2 of the path, in m^-2/3.

    Returns:
        The string.
    """
    if (float(path_m), float(cn2)) == (PATH_M, CN2):
        return ""
    return f"_L{path_m / 1e3:g}km_cn2{cn2:.0e}"


def reference_screen_count(path_m=PATH_M, cn2=CN2):
    """Give the screen count that the PRODUCTION sizer picks for one cell.

    This is the count of the UNMODIFIED standard preset, so it is the count of
    the backbone campaign that the sweep compares against. The sizer runs no
    trial, so the call is cheap.

    Args:
        path_m: the horizontal path length, in m.
        cn2:    the Cn2 of the path, in m^-2/3.

    Returns:
        The screen count as an int.
    """
    scn, geom = backbone.build_scenario(path_m, cn2, LAUNCH)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        _, plan, _ = turbulent_grid(scn, geom, preset=PRESET,
                                    L0_m=backbone.L0_M)
    return int(plan.z_m.size)


def override_plan(n_screens, path_m=PATH_M, cn2=CN2):
    """Build the grid and the EXACT n-screen plan of one override case.

    THE OVERRIDE. `_plan_terrestrial` takes the count
    N = max(min_screens, ceil(sigma_R^2 / sigma2_r_screen_max)) and it then
    cuts the path into N slabs of equal plane-wave Rytov weight (Andrews and
    Phillips, DOI 10.1117/3.626196, Ch. 8, Eq. (20), printed p. 264). A guard
    loop can raise N by one when the midpoint of the last slab overshoots the
    cap. So a preset copy with the cap effectively removed
    (sigma2_r_screen_max = 1e9) and the floor set to n gives EXACTLY n screens
    of the production shape, and the guard loop never fires.

    THE GRID stays the production grid: the function returns the grid of the
    UNMODIFIED standard sizing, and it reports the grid that the override
    itself would size, so the caller can check that the two agree.

    Args:
        n_screens: the wanted screen count.
        path_m:    the horizontal path length, in m.
        cn2:       the Cn2 of the path, in m^-2/3.

    Returns:
        A dict with the keys grid, plan, grid_override, warnings.

    Raises:
        AssertionError: the plan does not hold exactly n_screens screens.
    """
    scn, geom = backbone.build_scenario(path_m, cn2, LAUNCH)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        grid, _, _ = turbulent_grid(scn, geom, preset=PRESET,
                                    L0_m=backbone.L0_M)
        preset = dataclasses.replace(PRESETS[PRESET],
                                     name=f"{PRESET}_n{int(n_screens)}",
                                     sigma2_r_screen_max=1e9,
                                     min_screens=int(n_screens))
        grid_over, plan, _ = turbulent_grid(scn, geom, preset=preset,
                                            L0_m=backbone.L0_M)
    assert plan.z_m.size == int(n_screens), (plan.z_m.size, n_screens)
    return {"grid": grid, "plan": plan, "grid_override": grid_over,
            "warnings": sorted({str(w.message) for w in caught})}


def make_override_campaign(n_screens, block_size, path_m=PATH_M, cn2=CN2,
                           fft_backend=backbone.FFT_BACKEND,
                           screen_generator=backbone.SCREEN_GENERATOR):
    """Open (or make) the campaign of one screen count.

    The campaign takes the PINNED grid and the caller plan. Both enter the
    fingerprint, so each count is its own store.

    Args:
        n_screens:        the screen count.
        block_size:       the trials in one block.
        path_m:           the horizontal path length, in m.
        cn2:              the Cn2 of the path, in m^-2/3.
        fft_backend:      "numpy" or the "scipy" opt-in.
        screen_generator: "olb" or the "olb-lean" opt-in.

    Returns:
        The pair (Campaign, the dict of `override_plan`).
    """
    scn, geom = backbone.build_scenario(path_m, cn2, LAUNCH)
    spec = override_plan(n_screens, path_m, cn2)
    root = os.path.join(campaigns_root(),
                        case_tag(n_screens, path_m, cn2, fft_backend,
                                 screen_generator))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        camp = Campaign(scn, geom, root, seed=backbone.SEED, preset=PRESET,
                        block_size=int(block_size),
                        patch_radius_m=backbone.PATCH_RADIUS_M,
                        L0_m=backbone.L0_M, precision=backbone.PRECISION,
                        fft_backend=fft_backend,
                        screen_generator=screen_generator,
                        grid=spec["grid"], plan=spec["plan"])
    return camp, spec


def reference_campaign(block_size, path_m=PATH_M, cn2=CN2,
                       fft_backend=backbone.FFT_BACKEND,
                       screen_generator=backbone.SCREEN_GENERATOR):
    """Reopen the backbone campaign of one cell. It is NEVER run here.

    The call goes through `run_campaigns.make_campaign`, so the settings are
    the settings the store was built with and the fingerprint matches. A
    mismatch raises inside `Campaign`.

    Args:
        block_size:       the block size of the reference store (50 in the
                          backbone run).
        path_m:           the horizontal path length, in m.
        cn2:              the Cn2 of the path, in m^-2/3.
        fft_backend:      "numpy" or the "scipy" opt-in. It must match the
                          settings the reference store was made with.
        screen_generator: "olb" or the "olb-lean" opt-in. The same rule.

    Returns:
        The Campaign.
    """
    camp, _, _ = backbone.make_campaign(path_m, cn2, PRESET, LAUNCH,
                                        int(block_size), False, fft_backend,
                                        screen_generator)
    return camp


def reference_check(ref_camp, override_camp, n_trials, workers, say):
    """Prove that the override pipeline reproduces the production path.

    THE TEST. The override campaign at n = the reference screen count must
    give EXACTLY the reference trials. The function first compares the two
    PLANS (`z_m` and `r0_m`), because the `min_screens = n` route and the
    production cap route could in principle place the screens differently. It
    then runs `n_trials` override trials and compares the `collected_power`
    and the `smf_eta` of trials 0 to n_trials - 1 with `==`.

    A PASS says that the grid pinning, the plan override, the two speed opt-ins
    and the campaign reopen all reproduce the production path, so a delta at
    another count is a SCREEN-COUNT effect and nothing else.

    Args:
        ref_camp:      the reopened backbone Campaign.
        override_camp: the override Campaign at the reference screen count.
        n_trials:      the trials of the check.
        workers:       the process-pool size.
        say:           the log function.

    Returns:
        A JSON-ready dict.

    Raises:
        AssertionError: a stored value differs from the reference value.
    """
    out = {"n_trials": int(n_trials),
           "root": override_camp.root_dir,
           "n_screens": int(override_camp.plan.z_m.size)}
    same_z = np.array_equal(override_camp.plan.z_m, ref_camp.plan.z_m)
    same_r0 = np.array_equal(override_camp.plan.r0_m, ref_camp.plan.r0_m)
    out["plans_equal"] = bool(same_z and same_r0)
    if not out["plans_equal"]:
        say("  REF CHECK: the override plan DIFFERS from the production plan. "
            "The min_screens route and the cap route placed the screens "
            "differently, so the sweep is NOT like for like.")
        say(f"    z_m equal  : {same_z}")
        say(f"    r0_m equal : {same_r0}")
        say(f"    override z_m : {np.asarray(override_camp.plan.z_m)}")
        say(f"    reference z_m: {np.asarray(ref_camp.plan.z_m)}")
        out["passed"] = False
        return out
    say(f"  REF CHECK: the plans agree ({out['n_screens']} screens). "
        f"Running {n_trials} trials into {override_camp.root_dir}")
    override_camp.run(int(n_trials), workers=workers, progress=True)
    n_use = int(min(n_trials, override_camp.n_stored, ref_camp.n_stored))
    a = override_camp.load(n_use, fields=False).trials
    b = ref_camp.load(n_use, fields=False).trials
    got_p = np.array([t.collected_power for t in a], dtype=float)
    ref_p = np.array([t.collected_power for t in b], dtype=float)
    got_e = np.array([t.smf_eta for t in a], dtype=float)
    ref_e = np.array([t.smf_eta for t in b], dtype=float)
    out["n_compared"] = n_use
    out["max_abs_diff_collected_power"] = float(np.max(np.abs(got_p - ref_p)))
    out["max_abs_diff_smf_eta"] = float(np.max(np.abs(got_e - ref_e)))
    out["power_identical"] = bool(np.all(got_p == ref_p))
    out["eta_identical"] = bool(np.all(got_e == ref_e))
    out["passed"] = bool(out["power_identical"] and out["eta_identical"])
    say(f"    compared {n_use} trials: max |d collected_power| "
        f"{out['max_abs_diff_collected_power']:.3e}, max |d smf_eta| "
        f"{out['max_abs_diff_smf_eta']:.3e}")
    say(f"    bit-identical: collected_power {out['power_identical']}, "
        f"smf_eta {out['eta_identical']}")
    assert out["passed"], (
        "the override campaign at the reference screen count does NOT "
        "reproduce the reference trials. max |d collected_power| = "
        f"{out['max_abs_diff_collected_power']:.3e}, max |d smf_eta| = "
        f"{out['max_abs_diff_smf_eta']:.3e}. The sweep pipeline (the grid "
        "pinning, the plan override, the speed opt-ins or the campaign "
        "reopen) does not match the production path, so no count comparison "
        "is trustworthy.")
    return out


# ---------------------------------------------------------------------------
# The statistics (the idioms of validation/tail_convergence)
# ---------------------------------------------------------------------------

def _quantiles(fade_db, probs=EXCEEDANCE):
    """Give the fade that the link EXCEEDS each probability.

    A fade is a LARGE loss, so the fade exceeded a fraction q of the time is
    the (1 - q) quantile of the sample. Source:
    validation/tail_convergence/tail_convergence.py.

    Args:
        fade_db: the per-trial fade, in dB.
        probs:   the exceedance probabilities.

    Returns:
        A float array, one value for each probability.
    """
    x = np.asarray(fade_db, dtype=float)
    return np.quantile(x, [1.0 - q for q in probs])


def _bootstrap(x, fn, n_boot=N_BOOT, seed=BOOT_SEED):
    """Give the bootstrap half-width of a statistic.

    The function resamples the trials with replacement, it takes the statistic
    of each resample, and it gives the half-width of the central
    BOOT_INTERVAL band. The rng seed is FIXED, so a rerun of the analysis gives
    the same interval. Source:
    validation/tail_convergence/tail_convergence.py.

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


def _index(x):
    """Give the normalised variance var(x)/mean(x)^2 of a sample.

    This is the scintillation index sigma^2 = <x^2>/<x>^2 - 1. Source: Andrews
    and Phillips, DOI 10.1117/3.626196, Ch. 8. olb holds the definition in
    olb.turbulence.andrews.scintillation. Source of the helper:
    validation/tail_convergence/tail_convergence.py.
    """
    a = np.asarray(x, dtype=float)
    return float(a.var() / a.mean() ** 2)


def _centre_irradiance(result):
    """Give the irradiance of the CENTRE PIXEL of each stored trial.

    The centre pixel is the flat index (n // 2) * n + (n // 2), the axis pixel
    of the CircAperture convention that
    olb.waveoptics.turbulence.run._field_patch uses. Source:
    validation/outer_scale_tail/outer_scale_tail.py and
    validation/tail_convergence/tail_convergence.py.

    Args:
        result: a TurbWaveResult that holds the fields and the patch.

    Returns:
        A float array of |E|^2, one value for each trial.

    Raises:
        AssertionError: the centre pixel is outside the stored patch.
    """
    patch = result.patch
    n = int(patch.n)
    flat = (n // 2) * n + (n // 2)
    where = np.searchsorted(patch.indices, flat)
    assert where < patch.indices.size and patch.indices[where] == flat, \
        "the centre pixel is not inside the stored patch"
    return np.abs(result.fields[:, where]) ** 2


def _dir_bytes(path):
    """Give the total byte count of the files in one directory."""
    return int(sum(os.path.getsize(os.path.join(path, f))
                   for f in os.listdir(path)
                   if os.path.isfile(os.path.join(path, f))))


def measure_quantity(x):
    """Measure one quantity: its mean, its index and its fade quantiles.

    The fade is -10 log10(x / <x>), so it is POSITIVE for a value under the
    mean (loss is positive dB) and it needs no vacuum reference.

    Args:
        x: the per-trial values, a 1-D array of positive numbers.

    Returns:
        A JSON-ready dict.
    """
    a = np.asarray(x, dtype=float)
    fade = -10.0 * np.log10(a / a.mean())
    q = _quantiles(fade)
    half = np.atleast_1d(_bootstrap(fade, _quantiles))
    keys = [f"p{int(p * 100)}" for p in EXCEEDANCE]
    return {
        "n": int(a.size),
        "mean": float(a.mean()),
        "mean_loss_db": float(-10.0 * np.log10(a.mean())),
        "index": _index(a),
        "index_half": float(_bootstrap(a, _index)),
        "fades_db": {k: float(v) for k, v in zip(keys, q)},
        "fades_half_db": {k: float(v) for k, v in zip(keys, half)},
    }


def measure(camp, label, n_trials):
    """Measure one campaign: the four quantities and the run record.

    Args:
        camp:     the Campaign.
        label:    the case label.
        n_trials: the number of trials to read.

    Returns:
        A JSON-ready dict.
    """
    result = camp.load(n_trials, fields=True)
    power = np.array([t.collected_power for t in result.trials], dtype=float)
    eta = np.array([t.smf_eta for t in result.trials], dtype=float)
    point = _centre_irradiance(result)
    small = camp.recollect(aperture_m=SMALL_APERTURE_M, n_trials=n_trials)
    values = {"P10cm": power, "point": point, "smf_eta": eta, "P5cm": small}
    return {
        "case": label,
        "n_trials": int(power.size),
        "n_screens": int(camp.plan.z_m.size),
        "grid_n": int(camp.grid.n),
        "grid_size_m": float(camp.grid.size_m),
        "grid_pixel_m": float(camp.grid.size_m / camp.grid.n),
        "sigma2_r_max": float(camp.plan.sigma2_r.max()),
        "sigma2_r_sum": float(camp.plan.sigma2_r.sum()),
        "r0_total_m": float(camp.plan.r0_total_m),
        "plan_z_m": camp.plan.z_m.tolist(),
        "plan_sigma2_r": camp.plan.sigma2_r.tolist(),
        "root": camp.root_dir,
        "n_stored": int(camp.n_stored),
        "disk_bytes": _dir_bytes(camp.root_dir),
        "mean_wall_time_s": float(np.mean([t.wall_time_s
                                           for t in result.trials])),
        "quantities": {k: measure_quantity(v) for k, v in values.items()},
    }


# ---------------------------------------------------------------------------
# The comparison
# ---------------------------------------------------------------------------

def compare(row, ref, tolerance_db=TOLERANCE_DB):
    """Compare one case against the reference, and give the two verdicts.

    THE RESOLUTION VERDICT. A quantity PASSES when the p5 and the p1 delta
    both sit inside the bootstrap half-width of the REFERENCE, and the index
    ratio sits inside 1 +/- INDEX_TOLERANCE.

    THE TOLERANCE VERDICT. The case PASSES when every |p5| and |p1| delta of
    the TOLERANCE_QUANTITIES (the 10 cm bucket and the SMF coupling) is under
    `tolerance_db`.

    THE DELTA NOISE. The half-width above is the noise of the reference alone.
    The case is a Monte Carlo estimate of the same size, so the noise on the
    DELTA is about sqrt(2) times that, and the function reports it as
    `fade_delta_half_db`.

    Args:
        row:          the measured dict of the case.
        ref:          the measured dict of the reference.
        tolerance_db: the dB tolerance of the second verdict.

    Returns:
        A JSON-ready dict with the per-quantity deltas and the two verdicts.
    """
    out, failed, over = {}, [], []
    for name in QUANTITIES:
        a, b = ref["quantities"][name], row["quantities"][name]
        item = {
            "mean_delta_db": float(b["mean_loss_db"] - a["mean_loss_db"]),
            "index": b["index"], "index_ref": a["index"],
            "index_ratio": (float(b["index"] / a["index"])
                            if a["index"] > 0 else float("nan")),
            "fade_delta_db": {}, "fade_ref_half_db": {},
            "fade_delta_half_db": {}, "fade_inside": {},
        }
        for key in ("p10", "p5", "p1"):
            d = b["fades_db"][key] - a["fades_db"][key]
            half = a["fades_half_db"][key]
            item["fade_delta_db"][key] = float(d)
            item["fade_ref_half_db"][key] = float(half)
            # The two estimates are independent samples of the same size, so
            # the noise on their difference is about sqrt(2) times the noise
            # on one of them.
            item["fade_delta_half_db"][key] = float(np.sqrt(2.0) * half)
            item["fade_inside"][key] = bool(abs(d) <= half)
        bad = [k for k in ("p5", "p1") if not item["fade_inside"][k]]
        ratio = item["index_ratio"]
        if not np.isfinite(ratio) or abs(ratio - 1.0) > INDEX_TOLERANCE:
            bad.append("index")
        item["passed"] = not bad
        if bad:
            failed.append(f"{name}: " + ", ".join(bad))
        if name in TOLERANCE_QUANTITIES:
            wide = [k for k in ("p5", "p1")
                    if abs(item["fade_delta_db"][k]) > tolerance_db]
            item["within_tolerance"] = not wide
            if wide:
                over.append(f"{name}: " + ", ".join(
                    f"{k} {item['fade_delta_db'][k]:+.2f} dB" for k in wide))
        out[name] = item
    return {"quantities": out, "converged": not failed, "failed": failed,
            "tolerance_db": float(tolerance_db),
            "within_tolerance": not over, "over_tolerance": over}


# ---------------------------------------------------------------------------
# The printing
# ---------------------------------------------------------------------------

def make_say(path):
    """Give a print function that also appends to a log file, line by line."""
    def say(text=""):
        print(text, flush=True)
        with open(path, "a", encoding="utf-8") as fh:
            fh.write(text + "\n")
    return say


def print_table(rows, say):
    """Print a table with one space-padded column for each field."""
    width = [max(len(r[i]) for r in rows) for i in range(len(rows[0]))]
    for k, row in enumerate(rows):
        say("  " + "  ".join(v.rjust(width[i]) if i else v.ljust(width[i])
                             for i, v in enumerate(row)))
        if k == 0:
            say("  " + "  ".join("-" * w for w in width))


def grid_held(spec, camp):
    """Say whether the override sizing keeps the production grid.

    The pixel count must be EQUAL. The side is a float that the sizer rebuilds
    from the screen plan, so two plans of the same cell can differ in the last
    bits; the test takes a relative tolerance of 1e-9, far under any physical
    difference.

    Args:
        spec: the dict of `override_plan`.
        camp: the Campaign, which holds the PINNED production grid.

    Returns:
        True when the two grids agree.
    """
    over = spec["grid_override"]
    return bool(over.n == camp.grid.n
                and np.isclose(over.size_m, camp.grid.size_m, rtol=1e-9,
                               atol=0.0))


def cost_anchor(path_m, cn2):
    """Give the (s/trial, screens, workers) cost anchor of one cell.

    A cell that COST_ANCHORS does not name takes the 5 km anchor, the nearest
    measurement, so the projection is a rough guide only.

    Args:
        path_m: the horizontal path length, in m.
        cn2:    the Cn2 of the path, in m^-2/3.

    Returns:
        The triple.
    """
    return COST_ANCHORS.get((float(path_m), float(cn2)), FALLBACK_ANCHOR)


def dry_run_rows(specs, workers, path_m=PATH_M, cn2=CN2):
    """Give the dry-run table: the grid, the screens and the projected cost.

    The cost model is the measured anchor of the cell, scaled linearly by the
    screen count, and divided by the worker count. See COST_ANCHORS.

    Args:
        specs:   a list of dicts with the keys label, n_screens, trials, camp,
                 spec.
        workers: the process count of the run.
        path_m:  the horizontal path length, in m.
        cn2:     the Cn2 of the path, in m^-2/3.

    Returns:
        A list of string lists, the header first.
    """
    anchor_s, anchor_screens, anchor_workers = cost_anchor(path_m, cn2)
    rows = [["case", "screens", "n px", "side m", "px mm", "max sigma2_r",
             "grid held", "trials", "s/trial",
             "h at %d workers" % workers, "root"]]
    for s in specs:
        camp, spec = s["camp"], s["spec"]
        s_per = anchor_s * camp.plan.z_m.size / anchor_screens
        same = grid_held(spec, camp)
        trials = int(s["trials"])
        rows.append([
            s["label"],
            f"{camp.plan.z_m.size:d}",
            f"{camp.grid.n:d}",
            f"{camp.grid.size_m:.3f}",
            f"{camp.grid.size_m / camp.grid.n * 1e3:.2f}",
            f"{camp.plan.sigma2_r.max():.4f}",
            "yes" if same else "NO",
            f"{trials:d}",
            f"{s_per:.2f}",
            f"{trials * s_per * anchor_workers / max(workers, 1) / 3600.0:.2f}",
            # The directory name carries the opt-in suffix, so the dry run
            # shows WHICH store each count writes into.
            os.path.basename(camp.root_dir),
        ])
    return rows


def print_results(rows, deltas, say, ref_screens):
    """Print the ONE results table: the reference first, then each count.

    Args:
        rows:        the measured dicts, the reference first.
        deltas:      a dict case -> the compare() dict. The reference has none.
        say:         the log function.
        ref_screens: the screen count of the reference.
    """
    head = (f"{'case':<8s}{'scr':>4s}{'n':>6s}{'s2r_max':>9s}")
    for name in QUANTITIES:
        head += f"{name + ' s2':>17s}{name + ' p5':>17s}{name + ' p1':>17s}"
    say(head)
    say("-" * len(head))
    for row in rows:
        d = deltas.get(row["case"])
        line = (f"{row['case']:<8s}{row['n_screens']:4d}{row['n_trials']:6d}"
                f"{row['sigma2_r_max']:9.4f}")
        for name in QUANTITIES:
            q = row["quantities"][name]
            if d is None:
                line += (f"{q['index']:17.4f}"
                         f"{q['fades_db']['p5']:17.2f}"
                         f"{q['fades_db']['p1']:17.2f}")
            else:
                item = d["quantities"][name]
                line += (f"{q['index']:10.4f}(x{item['index_ratio']:4.2f})"
                         f"{q['fades_db']['p5']:10.2f}"
                         f"({item['fade_delta_db']['p5']:+5.2f})"
                         f"{q['fades_db']['p1']:10.2f}"
                         f"({item['fade_delta_db']['p1']:+5.2f})")
        say(line)
    say("  scr = the screen count. s2 = the scintillation index "
        "var(x)/mean(x)^2, with the ratio to the reference in brackets.")
    say(f"  p5 and p1 = the fade -10 log10(x/<x>) EXCEEDED 5 and 1 percent of "
        f"the time [dB], with the delta against the {ref_screens}-screen "
        f"reference.")
    for name in QUANTITIES:
        say(f"  {name:<8s} {QUANTITY_DOC[name]}")


def print_reference_bars(ref, say):
    """Print the bootstrap half-widths that the verdict tests against.

    It prints the REFERENCE bar and, next to it, the sqrt(2)-scaled bar of a
    DELTA: each count is a Monte Carlo estimate of the same size, so the noise
    on the difference of the two is about sqrt(2) times the noise on one.
    """
    root2 = float(np.sqrt(2.0))
    say(f"THE REFERENCE BOOTSTRAP HALF-WIDTHS ({BOOT_INTERVAL * 100:.0f} "
        f"percent, {N_BOOT} resamples). A delta inside these bars is not "
        f"resolved. The 'delta' bar is sqrt(2) times the reference bar, "
        f"because a count is an estimate of the same size.")
    for name in QUANTITIES:
        q = ref["quantities"][name]
        say(f"  {name:<8s} p10 +-{q['fades_half_db']['p10']:.2f}  "
            f"p5 +-{q['fades_half_db']['p5']:.2f}  "
            f"p1 +-{q['fades_half_db']['p1']:.2f} dB   "
            f"index {q['index']:.4f} +-{q['index_half']:.4f}   "
            f"mean loss {q['mean_loss_db']:.3f} dB")
        say(f"  {'':<8s} delta bar (x sqrt 2): "
            f"p10 +-{root2 * q['fades_half_db']['p10']:.2f}  "
            f"p5 +-{root2 * q['fades_half_db']['p5']:.2f}  "
            f"p1 +-{root2 * q['fades_half_db']['p1']:.2f} dB")


# ---------------------------------------------------------------------------
# The figure
# ---------------------------------------------------------------------------

def plot(rows, ref, fig_dir, ref_screens, path_m=PATH_M, cn2=CN2):
    """Draw the sweep figure. Give the list of the written paths.

    The figure holds two rows of three panels, one column for each of the
    three propagated quantities (the 100 mm power, the centre pixel and the
    fibre coupling). The top row is the index against the screen count, the
    bottom row the p5 and the p1 fade. The reference is a horizontal band of
    its bootstrap half-width.

    Args:
        rows:        the measured dicts of the sweep counts.
        ref:         the measured dict of the reference.
        fig_dir:     the output directory.
        ref_screens: the screen count of the reference.
        path_m:      the horizontal path length, in m.
        cn2:         the Cn2 of the path, in m^-2/3.

    Returns:
        A list of the written paths.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib is not installed; the figure is skipped.")
        return []

    os.makedirs(fig_dir, exist_ok=True)
    names = ("P10cm", "point", "smf_eta")
    counts = [r["n_screens"] for r in rows]
    fig, axes = plt.subplots(2, 3, figsize=(13.5, 7.5), sharex=True)
    for col, name in enumerate(names):
        q_ref = ref["quantities"][name]
        # ---- the index ----
        ax = axes[0][col]
        ax.errorbar(counts, [r["quantities"][name]["index"] for r in rows],
                    yerr=[r["quantities"][name]["index_half"] for r in rows],
                    marker="o", capsize=3, color="C0")
        ax.axhspan(q_ref["index"] - q_ref["index_half"],
                   q_ref["index"] + q_ref["index_half"],
                   color="0.7", alpha=0.5,
                   label=f"{ref_screens} screens")
        ax.axhline(q_ref["index"], color="0.3", lw=1)
        ax.set_title(name)
        ax.set_ylabel("index var/mean^2" if col == 0 else "")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
        # ---- the fades ----
        ax = axes[1][col]
        for key, colour in (("p5", "C1"), ("p1", "C3")):
            ax.errorbar(counts,
                        [r["quantities"][name]["fades_db"][key] for r in rows],
                        yerr=[r["quantities"][name]["fades_half_db"][key]
                              for r in rows],
                        marker="o", capsize=3, color=colour, label=key)
            ax.axhspan(q_ref["fades_db"][key] - q_ref["fades_half_db"][key],
                       q_ref["fades_db"][key] + q_ref["fades_half_db"][key],
                       color=colour, alpha=0.15)
            ax.axhline(q_ref["fades_db"][key], color=colour, lw=1, ls="--")
        ax.set_xlabel("screens")
        ax.set_ylabel("fade [dB]" if col == 0 else "")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
    fig.suptitle(f"The terrestrial screen-count sweep, {path_m / 1e3:g} km, "
                 f"Cn2 = {cn2:.0e}, {PRESET} preset. The band is the "
                 f"{ref_screens}-screen reference.")
    fig.tight_layout()
    path = os.path.join(fig_dir,
                        f"screen_count_sweep{file_stem(path_m, cn2)}.png")
    fig.savefig(path, dpi=140)
    plt.close(fig)
    return [path]


# ---------------------------------------------------------------------------
# The driver
# ---------------------------------------------------------------------------

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--path-km", type=float, default=PATH_M / 1e3,
                    help=f"the horizontal path length of the cell, in km "
                         f"(default {PATH_M / 1e3:g})")
    ap.add_argument("--cn2", type=float, default=CN2,
                    help=f"the Cn2 of the cell, in m^-2/3 (default {CN2:.0e})")
    ap.add_argument("--counts", nargs="+", type=int,
                    default=list(DEFAULT_COUNTS),
                    help="the screen counts of the sweep")
    ap.add_argument("--n-trials", type=int, default=DEFAULT_TRIALS,
                    help=f"the trials of each count, and of the reference "
                         f"slice (default {DEFAULT_TRIALS})")
    ap.add_argument("--tolerance-db", type=float, default=TOLERANCE_DB,
                    help=f"the dB tolerance of the second verdict (default "
                         f"{TOLERANCE_DB:g}). A count PASSES when every |p5| "
                         f"and |p1| delta of "
                         f"{' and '.join(TOLERANCE_QUANTITIES)} is under it.")
    ap.add_argument("--no-ref-check", action="store_true",
                    help="skip the reference-count bit-identity check")
    ap.add_argument("--ref-check-trials", type=int, default=REF_CHECK_TRIALS,
                    help=f"the trials of that check (default "
                         f"{REF_CHECK_TRIALS})")
    ap.add_argument("--workers", type=int, default=12,
                    help="the process-pool size of Campaign.run")
    ap.add_argument("--block-size", type=int, default=50,
                    help="the trials in one block of the sweep campaigns")
    ap.add_argument("--ref-block-size", type=int, default=backbone.FULL_BLOCK,
                    help="the block size the REFERENCE store was built with")
    ap.add_argument("--analyse-only", action="store_true",
                    help="skip the runs and read what is stored")
    ap.add_argument("--run-only", action="store_true",
                    help="store the trials and stop; no analysis, no figure. "
                         "Use it for a timing stage or a chained run, and "
                         "analyse ONE time at the end with --analyse-only.")
    ap.add_argument("--dry-run", action="store_true",
                    help="size every count, print the table, and run nothing")
    ap.add_argument("--fft-backend", choices=("numpy", "scipy"),
                    default=backbone.FFT_BACKEND,
                    help=f"the Forvard transform backend (default "
                         f"{backbone.FFT_BACKEND}). It goes to the REFERENCE "
                         f"reopen and to every override campaign, and 'scipy' "
                         f"gives every root the suffix _scipy.")
    ap.add_argument("--screen-generator", choices=("olb", "olb-lean"),
                    default=backbone.SCREEN_GENERATOR,
                    help=f"the phase-screen generator (default "
                         f"{backbone.SCREEN_GENERATOR}). The same rule, with "
                         f"the suffix _lean.")
    args = ap.parse_args(argv)

    counts = sorted({int(c) for c in args.counts})      # cheapest first
    path_m = float(args.path_km) * 1e3
    cn2 = float(args.cn2)
    stem = file_stem(path_m, cn2)
    ref_screens = reference_screen_count(path_m, cn2)
    ref_label = f"ref{ref_screens}"
    log_path = os.path.join(HERE, f"screen_count_sweep{stem}.log")
    if not args.dry_run:
        with open(log_path, "w", encoding="utf-8"):
            pass
        say = make_say(log_path)
    else:
        say = print

    say("The terrestrial screen-count convergence sweep (backlog 2-TC1), "
        f"{time.strftime('%Y-%m-%d %H:%M:%S')}")
    say(f"cell          : {path_m / 1e3:g} km, Cn2 = {cn2:.0e}, "
        f"{PRESET} preset, L0 = {backbone.L0_M:g} m, "
        f"{backbone.PRECISION} precision, seed {backbone.SEED}")
    say(f"counts        : {counts}   (the reference is {ref_screens} screens, "
        f"the count the production sizer picks, and it is NEVER run here)")
    say(f"trials        : {args.n_trials} for each case")
    say(f"ref check     : "
        f"{'OFF' if args.no_ref_check else str(args.ref_check_trials) + ' trials at n=' + str(ref_screens)}")
    say(f"tolerance     : {args.tolerance_db:g} dB on "
        f"{' and '.join(TOLERANCE_QUANTITIES)}")
    say(f"store         : {campaigns_root()}")
    say(f"speed opt-ins : fft {args.fft_backend}, screens "
        f"{args.screen_generator}, root suffix "
        f"{backbone.settings_suffix(args.fft_backend, args.screen_generator) or '(none)'}")
    say(f"mode          : "
        f"{'DRY RUN' if args.dry_run else ('ANALYSE ONLY' if args.analyse_only else 'RUN')}")
    say("correction    : NONE. Fidelity 2 applies no tip-tilt removal and no "
        "AO (backlog 2-AO).")
    say()

    # ---- size every case ----
    def size(n, label, trials):
        """Size one override campaign, and check that the grid is held."""
        camp, spec = make_override_campaign(n, args.block_size, path_m, cn2,
                                            args.fft_backend,
                                            args.screen_generator)
        assert grid_held(spec, camp), (
            f"the {n}-screen override moved the grid "
            f"({spec['grid_override'].n} px, {spec['grid_override'].size_m} m) "
            f"against the standard sizing ({camp.grid.n} px, "
            f"{camp.grid.size_m} m). The study pins the grid, so it passes the "
            f"standard grid explicitly; this assert only reports the case.")
        for text in spec["warnings"]:
            say(f"  n={n}: sizer warning: {text}")
        return {"n_screens": n, "label": label, "trials": int(trials),
                "camp": camp, "spec": spec}

    specs = [size(n, f"n{n}", args.n_trials) for n in counts]
    ref_spec = (None if args.no_ref_check
                else size(ref_screens, f"n{ref_screens} ref-check",
                          args.ref_check_trials))

    print_table(dry_run_rows(specs + ([ref_spec] if ref_spec else []),
                             args.workers, path_m, cn2), say)
    say()
    if args.dry_run:
        say("dry run: nothing ran. Drop --dry-run to store the trials.")
        return

    # ---- the reference store, and the bit-identity check ----
    ref_camp = None
    try:
        ref_camp = reference_campaign(args.ref_block_size, path_m, cn2,
                                      args.fft_backend, args.screen_generator)
    except (ValueError, OSError) as exc:
        say(f"WARNING: the reference campaign did not open: {exc}")

    ref_check = None
    if ref_spec is not None and ref_camp is not None and ref_camp.n_stored:
        # --analyse-only runs nothing, so the check then reads the trials the
        # override store already holds. `Campaign.run` of a count it already
        # holds stores nothing, so the call is a no-op there.
        n_check = min(args.ref_check_trials, ref_camp.n_stored)
        if args.analyse_only:
            n_check = min(n_check, ref_spec["camp"].n_stored)
        if n_check:
            ref_check = reference_check(ref_camp, ref_spec["camp"], n_check,
                                        args.workers, say)
        else:
            say("  REF CHECK: --analyse-only, and the check store holds no "
                "trial. The bit-identity check is SKIPPED.")
        say()
    elif ref_spec is not None:
        say("  REF CHECK: the reference store holds no trial, so the "
            "bit-identity check is SKIPPED.")
        say()

    # ---- the runs, cheapest first ----
    timing = {}
    if not args.analyse_only:
        for s in specs:
            camp, n = s["camp"], s["n_screens"]
            say(f"RUN n={n}: {args.n_trials} trials into {camp.root_dir}")
            n_before = int(camp.n_stored)
            t0 = time.perf_counter()
            n_done = camp.run(args.n_trials, workers=args.workers,
                              progress=True)
            wall = time.perf_counter() - t0
            # The rate counts the NEW trials only: a resumed store already
            # holds the old blocks (the same rule as campaign_resources.py).
            n_new = int(n_done) - n_before
            timing[str(n)] = {"wall_s": float(wall), "n_stored": int(n_done),
                              "n_new": n_new, "workers": args.workers,
                              "s_per_trial_new": (float(wall / n_new)
                                                  if n_new > 0 else None)}
            say(f"  n={n}: {n_done} trials on disk, {n_new} new in "
                f"{wall:.1f} s wall, "
                f"{(wall / n_new) if n_new > 0 else float('nan'):.2f} "
                f"s/trial (the new trials of this call)")
            say()

    if args.run_only:
        # A timing stage or a chained run: keep the store and the run lines,
        # and skip the analysis (it rebuilds the full grid for every stored
        # trial, minutes of one core).
        run_path = os.path.join(HERE, f"screen_count_sweep{stem}_runs.json")
        record = []
        if os.path.exists(run_path):
            with open(run_path, encoding="utf-8") as fh:
                record = json.load(fh)
        record.append({"time": time.strftime("%Y-%m-%d %H:%M:%S"),
                       "counts": counts, "n_trials": args.n_trials,
                       "workers": args.workers, "block_size": args.block_size,
                       "fft_backend": args.fft_backend,
                       "screen_generator": args.screen_generator,
                       "timing": timing})
        with open(run_path, "w", encoding="utf-8") as fh:
            json.dump(record, fh, indent=2)
        say(f"run only: the trials are stored. wrote {run_path}")
        return

    # ---- the analysis ----
    rows, deltas = [], {}
    ref_row = None
    if ref_camp is not None:
        n_ref = min(args.n_trials, ref_camp.n_stored)
        if n_ref == 0:
            say(f"WARNING: the reference store {ref_camp.root_dir} holds no "
                f"trial. The comparison is SKIPPED.")
        else:
            if n_ref < args.n_trials:
                say(f"WARNING: the reference holds {n_ref} trials, under the "
                    f"{args.n_trials} asked for. Every case is compared on "
                    f"{n_ref} trials.")
            ref_row = measure(ref_camp, ref_label, n_ref)
            rows.append(ref_row)

    n_common = ref_row["n_trials"] if ref_row else args.n_trials
    sweep_rows = []
    for s in specs:
        camp, n = s["camp"], s["n_screens"]
        if camp.n_stored == 0:
            say(f"  n={n}: no stored trial. Skipped.")
            continue
        n_use = min(n_common, camp.n_stored)
        if n_use < n_common:
            say(f"  n={n}: only {n_use} trials on disk, under the "
                f"{n_common} of the reference. The comparison is UNMATCHED.")
        row = measure(camp, f"n{n}", n_use)
        row["run_wall_s"] = timing.get(str(n), {}).get("wall_s")
        if ref_row is not None:
            deltas[row["case"]] = compare(row, ref_row, args.tolerance_db)
        rows.append(row)
        sweep_rows.append(row)
    say()

    if not sweep_rows:
        say("no case holds a trial: there is nothing to compare.")
    else:
        say("THE RESULTS (the reference row first)")
        print_results(rows, deltas, say, ref_screens)
        say()
        if ref_row is not None:
            print_reference_bars(ref_row, say)
            say()
            say("THE VERDICTS (two lines for each count)")
            say(f"  1 RESOLUTION (inside the reference noise). CONVERGED "
                f"means: every p5 and p1 delta sits inside the reference "
                f"bootstrap half-width, AND every index ratio sits inside "
                f"1 +/- {INDEX_TOLERANCE:.2f}.")
            say(f"  2 TOLERANCE. WITHIN means: every |p5| and |p1| delta of "
                f"{' and '.join(TOLERANCE_QUANTITIES)} is under "
                f"{args.tolerance_db:g} dB.")
            for row in sweep_rows:
                d = deltas[row["case"]]
                if d["converged"]:
                    say(f"  n={row['n_screens']:<3d} resolution: CONVERGED "
                        f"against {ref_screens} screens.")
                else:
                    say(f"  n={row['n_screens']:<3d} resolution: NOT "
                        f"CONVERGED: " + "; ".join(d["failed"]))
                if d["within_tolerance"]:
                    say(f"  {'':<6s} tolerance : WITHIN "
                        f"{args.tolerance_db:g} dB.")
                else:
                    say(f"  {'':<6s} tolerance : OVER {args.tolerance_db:g} "
                        f"dB: " + "; ".join(d["over_tolerance"]))
        say()

    stamp = {
        "study": "terrestrial_screen_count",
        "backlog": ["2-TC1"],
        "cell": {"path_length_m": path_m, "cn2_m_m23": cn2, "preset": PRESET,
                 "launch": LAUNCH, "L0_m": backbone.L0_M,
                 "precision": backbone.PRECISION, "seed": backbone.SEED,
                 "patch_radius_m": backbone.PATCH_RADIUS_M,
                 "wavelength_m": backbone.LAM, "waist_m": backbone.WAIST_M,
                 "aperture_m": APERTURE_M,
                 "small_aperture_m": SMALL_APERTURE_M},
        "counts": counts, "reference_screens": ref_screens,
        "n_trials": args.n_trials, "workers": args.workers,
        "block_size": args.block_size,
        "fft_backend": args.fft_backend,
        "screen_generator": args.screen_generator,
        "analyse_only": bool(args.analyse_only),
        "exceedance": list(EXCEEDANCE), "n_bootstrap": N_BOOT,
        "bootstrap_interval": BOOT_INTERVAL,
        "index_tolerance": INDEX_TOLERANCE,
        "tolerance_db": float(args.tolerance_db),
        "tolerance_quantities": list(TOLERANCE_QUANTITIES),
        "reference_check": ref_check,
        "timing": timing, "cases": rows, "comparison": deltas,
    }
    json_path = os.path.join(HERE,
                             f"screen_count_sweep{stem}_results.json")
    with open(json_path, "w", encoding="utf-8") as fh:
        json.dump(stamp, fh, indent=2)
    figs = (plot(sweep_rows, ref_row, os.path.join(HERE, "figures"),
                 ref_screens, path_m, cn2)
            if (sweep_rows and ref_row is not None) else [])
    say(f"wrote {json_path}")
    for f in figs:
        say(f"wrote {f}")
    say(f"wrote {log_path}")


if __name__ == '__main__':
    main()
