"""Check that the two speed OPT-INS do not change the physics.

THE QUESTION. The terrestrial campaigns run with two settings of record: the
numpy transform backend (`fft_backend="numpy"`) and the default screen
generator (`screen_generator="olb"`). Two speed opt-ins exist: the scipy
transform backend and the lean screen generator (`"olb-lean"`). Each opt-in
enters the campaign fingerprint, so an opt-in run writes its OWN store, with
the root suffix `_scipy_lean`. `validation/memory_cut/screen_generator_lean.py`
measured the two on ONE trial. This script measures them on a WHOLE campaign.

THE TEST. It opens the default store `L5km_cn23e-15_standard` and the opt-in
store `L5km_cn23e-15_standard_scipy_lean` through
`run_campaigns.make_campaign`. The two hold the SAME seed, the SAME grid and
the SAME screen plan, so trial k of one store is trial k of the other. The
script compares them trial by trial. It reports:

1. THE TRIALS. The largest and the median relative difference of
   `collected_power` and of `smf_eta`.
2. THE STATISTICS. The mean loss and the p5 and p1 fades of each store.

The expectation is a difference at the rounding level of single precision,
about 1e-6, and statistics that agree to much better than the Monte Carlo
error.

THE DATA IS NOT LOCAL. The campaign stores live on the desktop machine. Run
this script where the stores are, or point the environment variable
OLB_TERRESTRIAL_CAMPAIGNS_ROOT at them. Use `--dry-run` to see the two roots
without opening them.

Run it from the repository root:

    python -m validation.terrestrial_campaigns.optin_crosscheck --dry-run
    python -m validation.terrestrial_campaigns.optin_crosscheck
"""

import argparse
import json
import os
import time

import numpy as np

from .run_campaigns import (FULL_BLOCK, campaigns_root, make_campaign,
                            make_say, print_table)

HERE = os.path.dirname(os.path.abspath(__file__))

# The cell of the check. It is the one cell that holds BOTH a default store
# and an opt-in store.
CELL_PATH_M = 5e3
CELL_CN2 = 3e-15
CELL_PRESET = "standard"
CELL_LAUNCH = "collimated"

# The two settings pairs: the defaults of record, then the two opt-ins.
DEFAULT_SETTINGS = {"fft_backend": "numpy", "screen_generator": "olb"}
OPTIN_SETTINGS = {"fft_backend": "scipy", "screen_generator": "olb-lean"}

# The fade quantiles of the report. The availability targets of the study.
QUANTILES = (0.05, 0.01)


def open_campaign(settings, block_size):
    """Open one campaign of the check cell.

    The construction only reads the manifest. It runs no trial.

    Args:
        settings:   a dict with the keys fft_backend and screen_generator.
        block_size: the trials in one block.

    Returns:
        The Campaign.
    """
    camp, _, _ = make_campaign(CELL_PATH_M, CELL_CN2, CELL_PRESET, CELL_LAUNCH,
                               block_size, smoke=False, **settings)
    return camp


def loss_db(eta):
    """Give the loss in dB of a power fraction. Loss is positive dB.

    L = -10 log10(eta). This is the house sign rule of olb.

    Args:
        eta: an array of power fractions.

    Returns:
        An array of dB values.
    """
    eta = np.asarray(eta, dtype=float)
    return -10.0 * np.log10(eta)


def trial_arrays(camp, n_trials):
    """Give the per-trial scalars of a campaign.

    Args:
        camp:     the Campaign.
        n_trials: the number of trials to load.

    Returns:
        The pair (collected_power, smf_eta) as float arrays. A missing
        smf_eta gives an array of NaN.
    """
    res = camp.load(n_trials, fields=False)
    power = np.array([t.collected_power for t in res.trials], dtype=float)
    eta = np.array([np.nan if t.smf_eta is None else t.smf_eta
                    for t in res.trials], dtype=float)
    return power, eta


def relative_difference(a, b):
    """Give the largest and the median relative difference of two arrays.

    The reference is the first array, so the difference is |b / a - 1|.

    Args:
        a: the reference array.
        b: the comparison array.

    Returns:
        A dict with the keys max and median.
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    keep = np.isfinite(a) & np.isfinite(b) & (a != 0.0)
    if not np.any(keep):
        return {"max": None, "median": None, "n": 0}
    rel = np.abs(b[keep] / a[keep] - 1.0)
    return {"max": float(rel.max()), "median": float(np.median(rel)),
            "n": int(keep.sum())}


def statistics(power, eta):
    """Give the mean loss and the fade quantiles of one campaign.

    The total power fraction of a trial is the collected power times the
    fibre coupling efficiency. The loss is -10 log10 of that fraction. The
    fade DEPTH of a quantile is the loss at that quantile minus the mean
    loss, so it is positive dB of extra loss.

    Args:
        power: the collected-power array.
        eta:   the fibre coupling array. An array of NaN drops the fibre part.

    Returns:
        A JSON-ready dict.
    """
    total = power * eta if np.all(np.isfinite(eta)) else power
    keep = np.isfinite(total) & (total > 0.0)
    total = total[keep]
    mean_loss = float(loss_db(total.mean()))
    out = {"n_trials": int(total.size),
           "mean_eta": float(total.mean()),
           "mean_loss_db": mean_loss,
           "fibre_included": bool(np.all(np.isfinite(eta)))}
    for q in QUANTILES:
        # The q quantile of the POWER is the (1 - q) quantile of the loss: a
        # low power is a deep fade.
        low = float(np.quantile(total, q))
        out[f"p{int(q * 100)}_loss_db"] = float(loss_db(low))
        out[f"p{int(q * 100)}_fade_db"] = float(loss_db(low)) - mean_loss
    return out


def summary_rows(stats_a, stats_b):
    """Give the statistics table rows, the header first."""
    head = ["store", "trials", "mean eta", "mean loss dB"]
    for q in QUANTILES:
        head += [f"p{int(q * 100)} loss dB", f"p{int(q * 100)} fade dB"]
    rows = [head]
    for name, s in (("default", stats_a), ("scipy+lean", stats_b)):
        row = [name, f"{s['n_trials']:d}", f"{s['mean_eta']:.6e}",
               f"{s['mean_loss_db']:.4f}"]
        for q in QUANTILES:
            row += [f"{s[f'p{int(q * 100)}_loss_db']:.4f}",
                    f"{s[f'p{int(q * 100)}_fade_db']:+.4f}"]
        rows.append(row)
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--n-trials", type=int, default=None,
                    help="the trials to compare. The default takes every "
                         "trial that BOTH stores hold.")
    ap.add_argument("--block-size", type=int, default=FULL_BLOCK,
                    help=f"the trials in one block (default {FULL_BLOCK}). It "
                         f"must match the block size of the stores.")
    ap.add_argument("--dry-run", action="store_true",
                    help="print the two campaign roots and exit")
    args = ap.parse_args(argv)

    log_paths = () if args.dry_run else (os.path.join(HERE,
                                                      "optin_crosscheck.log"),)
    say = make_say(*log_paths)

    camp_a = open_campaign(DEFAULT_SETTINGS, args.block_size)
    camp_b = open_campaign(OPTIN_SETTINGS, args.block_size)

    say(f"opt-in cross-check, {time.strftime('%Y-%m-%d %H:%M:%S')}")
    say(f"cell          : {CELL_PATH_M / 1e3:g} km, Cn2 {CELL_CN2:.0e}, "
        f"{CELL_PRESET} preset, {CELL_LAUNCH} launch")
    say(f"store parent  : {campaigns_root()}")
    say(f"default store : {camp_a.root_dir}")
    say(f"opt-in store  : {camp_b.root_dir}")
    say(f"default       : fft {DEFAULT_SETTINGS['fft_backend']}, screens "
        f"{DEFAULT_SETTINGS['screen_generator']}")
    say(f"opt-in        : fft {OPTIN_SETTINGS['fft_backend']}, screens "
        f"{OPTIN_SETTINGS['screen_generator']}")
    say()

    if args.dry_run:
        say("dry run: no store opened. Drop --dry-run to compare the trials.")
        return

    n_a, n_b = camp_a.n_stored, camp_b.n_stored
    say(f"stored trials : default {n_a}, opt-in {n_b}")
    n = min(n_a, n_b)
    if args.n_trials is not None:
        n = min(n, args.n_trials)
    if n < 1:
        raise SystemExit(
            "optin_crosscheck: one store holds no trial. Run this script on "
            "the machine that holds the campaigns, or set "
            "OLB_TERRESTRIAL_CAMPAIGNS_ROOT.")
    say(f"compared      : {n} trials, trial for trial")
    say()

    power_a, eta_a = trial_arrays(camp_a, n)
    power_b, eta_b = trial_arrays(camp_b, n)

    diff = {"collected_power": relative_difference(power_a, power_b),
            "smf_eta": relative_difference(eta_a, eta_b)}
    for name, d in diff.items():
        if d["max"] is None:
            say(f"  {name:16s}: no comparable trial")
        else:
            say(f"  {name:16s}: max rel diff {d['max']:.3e}, "
                f"median {d['median']:.3e}, over {d['n']} trials")
    say()

    stats_a = statistics(power_a, eta_a)
    stats_b = statistics(power_b, eta_b)
    print_table(summary_rows(stats_a, stats_b), say)

    out = {"cell": {"path_length_m": CELL_PATH_M, "cn2_m_m23": CELL_CN2,
                    "preset": CELL_PRESET, "launch": CELL_LAUNCH},
           "default": {"root": camp_a.root_dir, "settings": DEFAULT_SETTINGS,
                       "n_stored": int(n_a), "fingerprint": camp_a.fingerprint,
                       "statistics": stats_a},
           "optin": {"root": camp_b.root_dir, "settings": OPTIN_SETTINGS,
                     "n_stored": int(n_b), "fingerprint": camp_b.fingerprint,
                     "statistics": stats_b},
           "n_compared": int(n),
           "relative_difference": diff}
    path = os.path.join(HERE, "optin_crosscheck_results.json")
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=2)
    say()
    say(f"outputs       : {path}")


if __name__ == "__main__":
    main()
