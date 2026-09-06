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

THE METHOD. The study PINS the grid and it moves the screens only. Each case is
one `Campaign` with a CALLER plan: the sweep asks `turbulent_grid` for the same
equal-Rytov-weight cut of `_plan_terrestrial` with the cap turned off
(`sigma2_r_screen_max = 1e9`) and the floor set to the wanted count
(`min_screens = n`), so the plan holds EXACTLY n screens of the production
shape. The grid comes from the UNMODIFIED standard sizing, and the script
asserts that the override does not move it.

THE REFERENCE is the 35-screen campaign of the terrestrial backbone run,
`validation/terrestrial_campaigns/campaigns/L10km_cn21e-14_standard`. This
script NEVER runs it: it reopens that store with the same `Campaign` settings
that `run_campaigns.make_campaign` built it with, and it reads the FIRST
`--n-trials` trials. Those trials are bit-identical to a shorter run of the same
seed, because the runner seeds trial k off (entropy, k). See
`olb.waveoptics.turbulence.campaign`.

THE CELL is exactly the backbone cell: a 10 km horizontal path at
Cn2 = 1e-14 m^(-2/3), a 5 mm waist from a 100 mm terminal, a 100 mm receive
aperture with a single-mode fibre, the standard preset, L0 = 25 m, single
precision, seed 20260906. The cell constants and the scenario builder are
IMPORTED from `validation.terrestrial_campaigns.run_campaigns`, so the two
studies cannot drift apart.

THE MEASURED QUANTITIES. Each case reports four quantities, and each one gives
a scintillation index and the fade quantiles p10, p5 and p1 of
-10 log10(x / <x>), with a 1000-resample bootstrap 68 percent interval:
  P 10 cm     the collected power of the 100 mm receive bucket.
  point       the irradiance of the CENTRE PIXEL of the receive field. Nothing
              averages a point, so it is the sharpest probe of the screen plan.
              The pixel is 4.48 mm on this grid.
  smf_eta     the single-mode-fibre coupling efficiency of that aperture.
  P 5 cm      the collected power of a 50 mm bucket, a post-hoc `recollect` of
              the SAME stored fields.
The index is the normalised variance var(x)/mean(x)^2. Source: Andrews and
Phillips, Laser Beam Propagation through Random Media, 2nd ed. (2005),
DOI 10.1117/3.626196, Ch. 8. olb holds the definition in
`olb.turbulence.andrews.scintillation`.

THE VERDICT. A count is CONVERGED when, for every one of those four
quantities, the p5 and the p1 delta against the 35-screen reference sits inside
the bootstrap half-width of the reference, AND the index ratio sits inside
1 +/- 0.10. Else it is NOT CONVERGED, and the line names the quantity that
fails.

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

THE TWO SPEED OPT-INS. `--fft-backend scipy` and `--screen-generator olb-lean`
go to BOTH the reference reopen and every override campaign, and each one
ENTERS the campaign fingerprint. So they give every root the suffix `_scipy`,
`_lean` or `_scipy_lean`, and the reference they name is the opt-in backbone
cell `L10km_cn21e-14_standard_scipy_lean`. A reopen with the wrong settings
RAISES a fingerprint mismatch. See validation/terrestrial_campaigns/README.md.
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

# The cell. It is the backbone cell, and the constants come from that module.
PATH_M = 10e3
CN2 = 1e-14
PRESET = "standard"
LAUNCH = "collimated"

# The reference campaign: the 35-screen backbone run. This script never runs it.
REF_SCREENS = 35
REF_LABEL = "ref35"

# The screen counts of the sweep, cheapest first.
DEFAULT_COUNTS = (5, 10, 15, 20)

# The cost model of the dry run. The backbone smoke run measured about 6.0 s
# of WALL time for one trial of this cell at 35 screens on an 8-worker pool
# (the pool time, not the time of one process), and a trial is close to linear
# in the screen count, so the projection scales by n / 35 and by 8 / workers.
REF_S_PER_TRIAL = 6.0
REF_WORKERS = 8

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
    "point": "the irradiance of the CENTRE PIXEL (4.48 mm on this grid)",
    "smf_eta": "the single-mode-fibre coupling efficiency",
    "P5cm": "the collected power of a 50 mm bucket (a post-hoc recollect)",
}


# ---------------------------------------------------------------------------
# The campaigns
# ---------------------------------------------------------------------------

def campaigns_root():
    """Give the parent directory of the override campaigns of this study."""
    return os.environ.get(ENV_ROOT) or os.path.join(HERE, "campaigns")


def case_tag(n_screens, fft_backend=backbone.FFT_BACKEND,
             screen_generator=backbone.SCREEN_GENERATOR):
    """Give the directory name of one override campaign.

    The two speed opt-ins enter the campaign fingerprint, so they give the
    root the same suffix that the backbone gives its own roots.

    Args:
        n_screens:        the screen count.
        fft_backend:      "numpy" or the "scipy" opt-in.
        screen_generator: "olb" or the "olb-lean" opt-in.

    Returns:
        The directory name.
    """
    return (backbone.cell_tag(PATH_M, CN2, PRESET, LAUNCH, False)
            + f"_n{int(n_screens)}"
            + backbone.settings_suffix(fft_backend, screen_generator))


def override_plan(n_screens):
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

    Returns:
        A dict with the keys grid, plan, grid_override, warnings.

    Raises:
        AssertionError: the plan does not hold exactly n_screens screens.
    """
    scn, geom = backbone.build_scenario(PATH_M, CN2, LAUNCH)
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


def make_override_campaign(n_screens, block_size,
                           fft_backend=backbone.FFT_BACKEND,
                           screen_generator=backbone.SCREEN_GENERATOR):
    """Open (or make) the campaign of one screen count.

    The campaign takes the PINNED grid and the caller plan. Both enter the
    fingerprint, so each count is its own store.

    Args:
        n_screens:        the screen count.
        block_size:       the trials in one block.
        fft_backend:      "numpy" or the "scipy" opt-in.
        screen_generator: "olb" or the "olb-lean" opt-in.

    Returns:
        The pair (Campaign, the dict of `override_plan`).
    """
    scn, geom = backbone.build_scenario(PATH_M, CN2, LAUNCH)
    spec = override_plan(n_screens)
    root = os.path.join(campaigns_root(),
                        case_tag(n_screens, fft_backend, screen_generator))
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


def reference_campaign(block_size, fft_backend=backbone.FFT_BACKEND,
                       screen_generator=backbone.SCREEN_GENERATOR):
    """Reopen the 35-screen backbone campaign. It is NEVER run here.

    The call goes through `run_campaigns.make_campaign`, so the settings are
    the settings the store was built with and the fingerprint matches. A
    mismatch raises inside `Campaign`.

    Args:
        block_size:       the block size of the reference store (50 in the
                          backbone run).
        fft_backend:      "numpy" or the "scipy" opt-in. It must match the
                          settings the reference store was made with.
        screen_generator: "olb" or the "olb-lean" opt-in. The same rule.

    Returns:
        The Campaign.
    """
    camp, _, _ = backbone.make_campaign(PATH_M, CN2, PRESET, LAUNCH,
                                        int(block_size), False, fft_backend,
                                        screen_generator)
    return camp


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

def compare(row, ref):
    """Compare one case against the reference, and give the verdict.

    A quantity PASSES when the p5 and the p1 delta both sit inside the
    bootstrap half-width of the REFERENCE, and the index ratio sits inside
    1 +/- INDEX_TOLERANCE.

    Args:
        row: the measured dict of the case.
        ref: the measured dict of the 35-screen reference.

    Returns:
        A JSON-ready dict with the per-quantity deltas and the verdict.
    """
    out, failed = {}, []
    for name in QUANTITIES:
        a, b = ref["quantities"][name], row["quantities"][name]
        item = {
            "mean_delta_db": float(b["mean_loss_db"] - a["mean_loss_db"]),
            "index": b["index"], "index_ref": a["index"],
            "index_ratio": (float(b["index"] / a["index"])
                            if a["index"] > 0 else float("nan")),
            "fade_delta_db": {}, "fade_ref_half_db": {}, "fade_inside": {},
        }
        for key in ("p10", "p5", "p1"):
            d = b["fades_db"][key] - a["fades_db"][key]
            half = a["fades_half_db"][key]
            item["fade_delta_db"][key] = float(d)
            item["fade_ref_half_db"][key] = float(half)
            item["fade_inside"][key] = bool(abs(d) <= half)
        bad = [k for k in ("p5", "p1") if not item["fade_inside"][k]]
        ratio = item["index_ratio"]
        if not np.isfinite(ratio) or abs(ratio - 1.0) > INDEX_TOLERANCE:
            bad.append("index")
        item["passed"] = not bad
        if bad:
            failed.append(f"{name}: " + ", ".join(bad))
        out[name] = item
    return {"quantities": out, "converged": not failed, "failed": failed}


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


def dry_run_rows(specs, n_trials, workers):
    """Give the dry-run table: the grid, the screens and the projected cost.

    The cost model is REF_S_PER_TRIAL at REF_SCREENS screens, scaled linearly
    by the screen count, and divided by the worker count.

    Args:
        specs:    a list of dicts with the keys n_screens, camp, spec.
        n_trials: the trials of each case.
        workers:  the process count of the run.

    Returns:
        A list of string lists, the header first.
    """
    rows = [["case", "screens", "n px", "side m", "px mm", "max sigma2_r",
             "grid held", "s/trial", "h at %d workers" % workers, "root"]]
    for s in specs:
        camp, spec = s["camp"], s["spec"]
        s_per = REF_S_PER_TRIAL * camp.plan.z_m.size / REF_SCREENS
        same = ((spec["grid_override"].n, spec["grid_override"].size_m)
                == (camp.grid.n, camp.grid.size_m))
        rows.append([
            f"n{int(s['n_screens'])}",
            f"{camp.plan.z_m.size:d}",
            f"{camp.grid.n:d}",
            f"{camp.grid.size_m:.3f}",
            f"{camp.grid.size_m / camp.grid.n * 1e3:.2f}",
            f"{camp.plan.sigma2_r.max():.4f}",
            "yes" if same else "NO",
            f"{s_per:.2f}",
            f"{n_trials * s_per * REF_WORKERS / max(workers, 1) / 3600.0:.2f}",
            # The directory name carries the opt-in suffix, so the dry run
            # shows WHICH store each count writes into.
            os.path.basename(camp.root_dir),
        ])
    return rows


def print_results(rows, deltas, say):
    """Print the ONE results table: the reference first, then each count.

    Args:
        rows:   the measured dicts, the reference first.
        deltas: a dict case -> the compare() dict. The reference has none.
        say:    the log function.
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
    say("  p5 and p1 = the fade -10 log10(x/<x>) EXCEEDED 5 and 1 percent of "
        "the time [dB], with the delta against the 35-screen reference.")
    for name in QUANTITIES:
        say(f"  {name:<8s} {QUANTITY_DOC[name]}")


def print_reference_bars(ref, say):
    """Print the bootstrap half-widths that the verdict tests against."""
    say(f"THE REFERENCE BOOTSTRAP HALF-WIDTHS ({BOOT_INTERVAL * 100:.0f} "
        f"percent, {N_BOOT} resamples). A delta inside these bars is not "
        f"resolved.")
    for name in QUANTITIES:
        q = ref["quantities"][name]
        say(f"  {name:<8s} p10 +-{q['fades_half_db']['p10']:.2f}  "
            f"p5 +-{q['fades_half_db']['p5']:.2f}  "
            f"p1 +-{q['fades_half_db']['p1']:.2f} dB   "
            f"index {q['index']:.4f} +-{q['index_half']:.4f}   "
            f"mean loss {q['mean_loss_db']:.3f} dB")


# ---------------------------------------------------------------------------
# The figure
# ---------------------------------------------------------------------------

def plot(rows, ref, fig_dir):
    """Draw the sweep figure. Give the list of the written paths.

    The figure holds two rows of three panels, one column for each of the
    three propagated quantities (the 100 mm power, the centre pixel and the
    fibre coupling). The top row is the index against the screen count, the
    bottom row the p5 and the p1 fade. The reference is a horizontal band of
    its bootstrap half-width.

    Args:
        rows:    the measured dicts of the sweep counts.
        ref:     the measured dict of the reference.
        fig_dir: the output directory.

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
                   label=f"{REF_SCREENS} screens")
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
    fig.suptitle("The terrestrial screen-count sweep, 10 km, Cn2 = 1e-14, "
                 f"standard preset. The band is the {REF_SCREENS}-screen "
                 f"reference.")
    fig.tight_layout()
    path = os.path.join(fig_dir, "screen_count_sweep.png")
    fig.savefig(path, dpi=140)
    plt.close(fig)
    return [path]


# ---------------------------------------------------------------------------
# The driver
# ---------------------------------------------------------------------------

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--counts", nargs="+", type=int,
                    default=list(DEFAULT_COUNTS),
                    help="the screen counts of the sweep")
    ap.add_argument("--n-trials", type=int, default=1000,
                    help="the trials of each count, and of the reference")
    ap.add_argument("--workers", type=int, default=12,
                    help="the process-pool size of Campaign.run")
    ap.add_argument("--block-size", type=int, default=50,
                    help="the trials in one block of the sweep campaigns")
    ap.add_argument("--ref-block-size", type=int, default=backbone.FULL_BLOCK,
                    help="the block size the REFERENCE store was built with")
    ap.add_argument("--analyse-only", action="store_true",
                    help="skip the runs and read what is stored")
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
    log_path = os.path.join(HERE, "screen_count_sweep.log")
    if not args.dry_run:
        with open(log_path, "w", encoding="utf-8"):
            pass
        say = make_say(log_path)
    else:
        say = print

    say("The terrestrial screen-count convergence sweep (backlog 2-TC1), "
        f"{time.strftime('%Y-%m-%d %H:%M:%S')}")
    say(f"cell          : {PATH_M / 1e3:g} km, Cn2 = {CN2:.0e}, "
        f"{PRESET} preset, L0 = {backbone.L0_M:g} m, "
        f"{backbone.PRECISION} precision, seed {backbone.SEED}")
    say(f"counts        : {counts}   (the reference is {REF_SCREENS} screens, "
        f"and it is NEVER run here)")
    say(f"trials        : {args.n_trials} for each case")
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
    specs = []
    for n in counts:
        camp, spec = make_override_campaign(n, args.block_size,
                                            args.fft_backend,
                                            args.screen_generator)
        same = ((spec["grid_override"].n, spec["grid_override"].size_m)
                == (camp.grid.n, camp.grid.size_m))
        assert same, (
            f"the {n}-screen override moved the grid "
            f"({spec['grid_override'].n} px, {spec['grid_override'].size_m} m) "
            f"against the standard sizing ({camp.grid.n} px, "
            f"{camp.grid.size_m} m). The study pins the grid, so it passes the "
            f"standard grid explicitly; this assert only reports the case.")
        specs.append({"n_screens": n, "camp": camp, "spec": spec})
        for text in spec["warnings"]:
            say(f"  n={n}: sizer warning: {text}")

    print_table(dry_run_rows(specs, args.n_trials, args.workers), say)
    say()
    if args.dry_run:
        say("dry run: nothing ran. Drop --dry-run to store the trials.")
        return

    # ---- the runs, cheapest first ----
    timing = {}
    if not args.analyse_only:
        for s in specs:
            camp, n = s["camp"], s["n_screens"]
            say(f"RUN n={n}: {args.n_trials} trials into {camp.root_dir}")
            t0 = time.perf_counter()
            n_done = camp.run(args.n_trials, workers=args.workers,
                              progress=True)
            wall = time.perf_counter() - t0
            timing[str(n)] = {"wall_s": float(wall), "n_stored": int(n_done)}
            say(f"  n={n}: {n_done} trials on disk, {wall:.1f} s wall, "
                f"{wall / max(n_done, 1):.2f} s/trial (this call)")
            say()

    # ---- the analysis ----
    rows, deltas = [], {}
    ref_row = None
    try:
        ref_camp = reference_campaign(args.ref_block_size, args.fft_backend,
                                      args.screen_generator)
        n_ref = min(args.n_trials, ref_camp.n_stored)
        if n_ref == 0:
            say(f"WARNING: the reference store {ref_camp.root_dir} holds no "
                f"trial. The comparison is SKIPPED.")
        else:
            if n_ref < args.n_trials:
                say(f"WARNING: the reference holds {n_ref} trials, under the "
                    f"{args.n_trials} asked for. Every case is compared on "
                    f"{n_ref} trials.")
            ref_row = measure(ref_camp, REF_LABEL, n_ref)
            rows.append(ref_row)
    except (ValueError, OSError) as exc:
        say(f"WARNING: the reference campaign did not open: {exc}")

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
            deltas[row["case"]] = compare(row, ref_row)
        rows.append(row)
        sweep_rows.append(row)
    say()

    if not sweep_rows:
        say("no case holds a trial: there is nothing to compare.")
    else:
        say("THE RESULTS (the reference row first)")
        print_results(rows, deltas, say)
        say()
        if ref_row is not None:
            print_reference_bars(ref_row, say)
            say()
            say("THE VERDICT")
            say(f"  CONVERGED means: every p5 and p1 delta sits inside the "
                f"reference bootstrap half-width, AND every index ratio sits "
                f"inside 1 +/- {INDEX_TOLERANCE:.2f}.")
            for row in sweep_rows:
                d = deltas[row["case"]]
                if d["converged"]:
                    say(f"  n={row['n_screens']:<3d} CONVERGED against "
                        f"{REF_SCREENS} screens.")
                else:
                    say(f"  n={row['n_screens']:<3d} NOT CONVERGED: "
                        + "; ".join(d["failed"]))
        say()

    stamp = {
        "study": "terrestrial_screen_count",
        "backlog": ["2-TC1"],
        "cell": {"path_length_m": PATH_M, "cn2_m_m23": CN2, "preset": PRESET,
                 "launch": LAUNCH, "L0_m": backbone.L0_M,
                 "precision": backbone.PRECISION, "seed": backbone.SEED,
                 "patch_radius_m": backbone.PATCH_RADIUS_M,
                 "wavelength_m": backbone.LAM, "waist_m": backbone.WAIST_M,
                 "aperture_m": APERTURE_M,
                 "small_aperture_m": SMALL_APERTURE_M},
        "counts": counts, "reference_screens": REF_SCREENS,
        "n_trials": args.n_trials, "workers": args.workers,
        "block_size": args.block_size,
        "fft_backend": args.fft_backend,
        "screen_generator": args.screen_generator,
        "analyse_only": bool(args.analyse_only),
        "exceedance": list(EXCEEDANCE), "n_bootstrap": N_BOOT,
        "bootstrap_interval": BOOT_INTERVAL,
        "index_tolerance": INDEX_TOLERANCE,
        "timing": timing, "cases": rows, "comparison": deltas,
    }
    json_path = os.path.join(HERE, "screen_count_sweep_results.json")
    with open(json_path, "w", encoding="utf-8") as fh:
        json.dump(stamp, fh, indent=2)
    figs = (plot(sweep_rows, ref_row, os.path.join(HERE, "figures"))
            if (sweep_rows and ref_row is not None) else [])
    say(f"wrote {json_path}")
    for f in figs:
        say(f"wrote {f}")
    say(f"wrote {log_path}")


if __name__ == '__main__':
    main()
