"""Compare the tilt-sensing routes of a terrestrial TipTilt correction.

THE QUESTION. A terrestrial fidelity-2 tip-tilt correction senses the tilt of
the stored receive field. The default route is the WRAPPED-GRADIENT SLOPES
(`olb.waveoptics.compensation.slopes`). That route differences adjacent pixels,
so it ALIASES where the local phase step passes pi per pixel. On a strong path
(the 10 km 2-TC cells) many trials pass that step, so the slope-sensed tilt is
not trustworthy there.

THE NEW ROUTE. The GRADIENT (centroid) tilt reads the far-field intensity
centroid (`olb.waveoptics.compensation.gtilt`). It reads the complex field, so
it never unwraps a phase and it does not alias a local step: it holds the
overall tilt where the slopes alias. The "auto" route is the slopes with a
per-trial fall-back to the G-tilt where the slope step passes pi.

WHAT THIS SCRIPT MEASURES. For each stored 2-TC terrestrial campaign it runs
FOUR post-hoc passes over the 10 cm single-mode-fibre receiver, with NO new
propagation:
  - uncorrected                (recouple)
  - TipTilt sensed by slopes   (recouple_compensated, source="slopes")
  - TipTilt sensed by G-tilt   (recouple_compensated, source="gtilt")
  - TipTilt sensed by auto     (recouple_compensated, source="auto")
It also reads the per-trial aliasing step. The decisive metric is the FIBRE
COUPLING FADE: the deep-fade (p5, p1) coupling loss of each route. The route
that centres the focused spot better simply couples better, so this sidesteps
the G-tilt-against-Z-tilt question. The mechanism table then bins the trials by
the aliasing step and shows the two routes AGREE below pi and DIVERGE above it.

THE DATA. The campaigns are the terrestrial 2-TC backbone
(validation/terrestrial_campaigns/): 2 / 5 / 10 km paths, Cn2 = 3e-15 / 1e-14,
a collimated 5 mm launch into a 10 cm SMF receiver, L0 = 25 m, single
precision. The stored complex field patch of each trial holds the UNCORRECTED
field, so every route here is a post-hoc read of the SAME atmosphere. The store
lives on bigfraw under the campaigns root; run this script there.

Sources:
- G. A. Tyler, "Bandwidth considerations for tracking through turbulence,"
  J. Opt. Soc. Am. A 11(1), 358-367 (1994), DOI 10.1364/JOSAA.11.000358. The
  Z-tilt / G-tilt split.
- R. J. Noll, J. Opt. Soc. Am. 66(3), 207-211 (1976),
  DOI 10.1364/JOSA.66.000207. The Zernike (tip-tilt) modes and their order.

Run it from the repository root:

    python -m validation.gtilt_sensing.gtilt_vs_slopes --dry-run
    python -m validation.gtilt_sensing.gtilt_vs_slopes --workers 8
    python -m validation.gtilt_sensing.gtilt_vs_slopes --cells 10km:1e-14:standard
"""

import argparse
import json
import math
import os
import time

import numpy as np

from olb.terminal import SMF, TipTilt
from olb.waveoptics.compensation import circle, max_abs_step, wrapped_gradient

from validation.terrestrial_campaigns.run_campaigns import (
    RX_APERTURE_M, campaigns_root, make_campaign)

HERE = os.path.dirname(os.path.abspath(__file__))

# The aliasing threshold. The wrapped-gradient slope aliases above pi per
# pixel; the runner warns above this value. Source:
# olb.waveoptics.turbulence.run.SLOPE_STEP_WARN_RAD.
STEP_WARN_RAD = 2.8

# The receive fibre of the 2-TC study.
DETECTOR = SMF(optimal_focus=True, sensitivity_dbm=-40)
STACK = [TipTilt()]


class _StepFn:
    """The per-trial aliasing step, over the receive aperture.

    The step is the largest wrapped phase difference between two adjacent
    pixels inside the aperture. It is the quantity the runner tests, but
    measured over the mask so the zero pixels outside the stored patch disc do
    not read as a step. It must pickle, so it is a top-level class.
    """

    def __init__(self, aperture_m, obscuration_ratio=0.0):
        self.aperture_m = float(aperture_m)
        self.obscuration_ratio = float(obscuration_ratio)

    def __call__(self, rec):
        mask = rec.context.get("step_mask")
        if mask is None:
            side = rec.patch.crop().side
            mask = circle(side, self.aperture_m / rec.patch.pixel_m,
                          self.obscuration_ratio)
            rec.context["step_mask"] = mask
        sx, sy = wrapped_gradient(rec.array, mask=mask)
        return max_abs_step(sx, sy)


def open_campaign(path_m, cn2, preset, launch, fft_backend, screen_generator,
                  block_size):
    """Reopen one stored 2-TC campaign from its manifest.

    The function builds the campaign the SAME way the 2-TC runner did, so the
    root and the fingerprint match and it reopens the store. It runs no trial.

    Args:
        path_m:           the path length, in m.
        cn2:              the Cn2, in m^-2/3.
        preset:           the quality preset name.
        launch:           "collimated" or "diverged".
        fft_backend:      the stored FFT backend.
        screen_generator: the stored screen generator.
        block_size:       the stored block size.

    Returns:
        The Campaign, or None when the store holds no trial.
    """
    camp, _, _ = make_campaign(path_m, cn2, preset, launch, block_size, False,
                               fft_backend, screen_generator)
    return camp if camp.n_stored > 0 else None


def _settings_from_name(name):
    """Read the two speed opt-ins from a campaign directory name.

    The directory name carries the opt-in suffix (_scipy, _lean or
    _scipy_lean), so it is the source of truth for the FFT backend and the
    screen generator, even when an older cell.json holds no such key. See
    validation.terrestrial_campaigns.run_campaigns.settings_suffix.

    Args:
        name: the campaign directory name.

    Returns:
        The pair (fft_backend, screen_generator).
    """
    fft = "scipy" if "_scipy" in name else "numpy"
    gen = "olb-lean" if "_lean" in name else "olb"
    return fft, gen


def discover(root, launch):
    """Find the finished campaigns under a root, from their cell.json.

    Args:
        root:   the campaigns parent directory.
        launch: keep only campaigns of this launch.

    Returns:
        A list of dicts with the reopen parameters and the cell metadata, one
        for each directory that holds a cell.json and at least one block.
    """
    out = []
    if not os.path.isdir(root):
        return out
    for name in sorted(os.listdir(root)):
        cell_json = os.path.join(root, name, "cell.json")
        if name.endswith("_smoke") or not os.path.isfile(cell_json):
            continue
        with open(cell_json, encoding="utf-8") as fh:
            rec = json.load(fh)
        if rec.get("launch", "collimated") != launch:
            continue
        # The directory-name suffix is the source of truth for the opt-ins,
        # because an older cell.json may hold no fft_backend key.
        fft, gen = _settings_from_name(name)
        camp_meta = rec.get("campaign", {})
        out.append({
            "name": name,
            "cell": rec["cell"],
            "path_m": rec["path_length_m"],
            "cn2": rec["cn2_m_m23"],
            "preset": rec["preset"],
            "launch": rec.get("launch", "collimated"),
            "sigma2_R": rec["sigma2_R_plane"],
            "r0_m": rec["r0_m"],
            "grid_n": rec["grid"]["n"],
            "screens": rec["screens"]["n"],
            "fft_backend": camp_meta.get("fft_backend", fft),
            "screen_generator": camp_meta.get("screen_generator", gen),
            "block_size": camp_meta.get("block_size", 50),
        })
    return out


def loss_db(eta):
    """Give the coupling loss in dB (positive dB) of a coupling fraction."""
    eta = np.asarray(eta, dtype=float)
    eta = np.where(eta > 0.0, eta, np.nan)
    return -10.0 * np.log10(eta)


def quantile_losses(eta):
    """Give the coupling loss at the mean and the deep-fade quantiles.

    Loss is POSITIVE dB. A DEEP fade is a LOW coupling, so it is a HIGH loss,
    at a LOW percentile of the coupling. The report gives the mean loss and the
    loss at the 50th, 5th and 1st percentiles of the coupling.

    Args:
        eta: the coupling fraction of each trial.

    Returns:
        A dict of the mean and the p50, p5 and p1 losses, in dB.
    """
    eta = np.asarray(eta, dtype=float)
    finite = eta[np.isfinite(eta) & (eta > 0.0)]
    if finite.size == 0:
        return {"mean_db": float("nan"), "p50_db": float("nan"),
                "p5_db": float("nan"), "p1_db": float("nan"),
                "mean_eta": float("nan")}
    return {
        "mean_eta": float(finite.mean()),
        "mean_db": float(-10.0 * math.log10(finite.mean())),
        "p50_db": float(-10.0 * math.log10(np.percentile(finite, 50))),
        "p5_db": float(-10.0 * math.log10(np.percentile(finite, 5))),
        "p1_db": float(-10.0 * math.log10(np.percentile(finite, 1))),
    }


def analyse_cell(spec, args, say):
    """Run the four routes and the step metric on one campaign.

    Args:
        spec: the discovery dict of the cell.
        args: the parsed arguments.
        say:  the log function.

    Returns:
        A JSON-ready result dict.
    """
    camp = open_campaign(spec["path_m"], spec["cn2"], spec["preset"],
                         spec["launch"], spec["fft_backend"],
                         spec["screen_generator"], spec["block_size"])
    if camp is None:
        say(f"  {spec['cell']}: no stored trial, skipped")
        return None
    a = RX_APERTURE_M
    n = args.n_trials
    w = args.workers
    say(f"--- {spec['cell']} ({spec['name']}) ---")
    say(f"  sigma_R^2 {spec['sigma2_R']:.2f}, r0 {spec['r0_m'] * 100:.2f} cm, "
        f"grid {spec['grid_n']} px, {spec['screens']} screens, "
        f"{camp.n_stored} trials on disk")

    t0 = time.time()
    step = camp.map_trials(_StepFn(a), n_trials=n, workers=w, fields=True,
                           screen_phase=False)
    eta_u = camp.recouple(DETECTOR, aperture_m=a, n_trials=n, workers=w)
    eta_s = camp.recouple_compensated(STACK, DETECTOR, aperture_m=a,
                                      n_trials=n, source="slopes", workers=w)
    eta_g = camp.recouple_compensated(STACK, DETECTOR, aperture_m=a,
                                      n_trials=n, source="gtilt", workers=w)
    wall = time.time() - t0

    routes = {"uncorrected": eta_u, "slopes": eta_s, "gtilt": eta_g}
    losses = {k: quantile_losses(v) for k, v in routes.items()}

    # The aliasing census. A wrapped step CANNOT pass pi (the wrap caps it), so
    # the useful proxy is the fraction of trials whose largest step sits near
    # pi, above STEP_WARN_RAD.
    n_trials = int(step.size)
    frac_warn = float(np.mean(step > STEP_WARN_RAD))

    # THE MECHANISM. Bin the trials by the aliasing proxy. In the low-step bin
    # the slope fit is trustworthy; in the near-pi bin the slope step is at its
    # wrap limit. The metric is the MEDIAN coupling of the bin, so one deep
    # outlier does not move it.
    below = step <= STEP_WARN_RAD
    above = ~below
    bins = {}
    for label, sel in (("step<=warn", below), ("step>warn", above)):
        if not np.any(sel):
            bins[label] = {"n": 0}
            continue
        bins[label] = {
            "n": int(sel.sum()),
            "median_eta_uncorrected": float(np.nanmedian(eta_u[sel])),
            "median_eta_slopes": float(np.nanmedian(eta_s[sel])),
            "median_eta_gtilt": float(np.nanmedian(eta_g[sel])),
        }

    say(f"  aliasing   : max step {step.max():.2f} rad/px, "
        f"{100 * frac_warn:.1f}% of trials past {STEP_WARN_RAD} (near pi)")
    say(f"  {'route':<12}{'mean dB':>9}{'p50 dB':>9}{'p5 dB':>9}{'p1 dB':>9}")
    for name in ("uncorrected", "slopes", "gtilt"):
        L = losses[name]
        say(f"  {name:<12}{L['mean_db']:>9.2f}{L['p50_db']:>9.2f}"
            f"{L['p5_db']:>9.2f}{L['p1_db']:>9.2f}")
    d_g = losses["slopes"]["p5_db"] - losses["gtilt"]["p5_db"]
    say(f"  p5 gain over slopes: gtilt {d_g:+.2f} dB "
        f"(positive = G-tilt fades less)")
    for label in ("step<=warn", "step>warn"):
        b = bins[label]
        if b["n"] == 0:
            say(f"  {label:<10}: no trial")
            continue
        say(f"  {label:<10} ({b['n']:5d} trials): median eta "
            f"slopes {b['median_eta_slopes']:.4f}, "
            f"gtilt {b['median_eta_gtilt']:.4f}")
    say(f"  ({wall:.1f} s for post-hoc passes)")
    say()

    if args.figures:
        _plot_cdf(spec, routes, say)

    return {
        "cell": spec["cell"], "name": spec["name"],
        "path_m": spec["path_m"], "cn2": spec["cn2"], "preset": spec["preset"],
        "sigma2_R": spec["sigma2_R"], "r0_m": spec["r0_m"],
        "grid_n": spec["grid_n"], "screens": spec["screens"],
        "n_trials": n_trials, "aperture_m": a,
        "max_step_rad": float(step.max()),
        "frac_step_gt_warn": frac_warn, "frac_step_gt_pi": frac_pi,
        "losses": losses, "bins": bins, "wall_s": wall,
    }


def _plot_cdf(spec, routes, say):
    """Plot the coupling-loss CDF of the four routes for one cell."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:                       # pragma: no cover
        say(f"  (no figure: {exc})")
        return
    fig, ax = plt.subplots(figsize=(6.0, 4.2))
    for name, eta in routes.items():
        L = np.sort(loss_db(eta))
        L = L[np.isfinite(L)]
        if L.size == 0:
            continue
        y = np.linspace(0.0, 1.0, L.size)
        ax.plot(L, y, label=name, lw=1.6)
    ax.set_xlabel("SMF coupling loss (dB)")
    ax.set_ylabel("cumulative fraction of trials")
    ax.set_title(f"{spec['cell']}  (sigma_R^2 = {spec['sigma2_R']:.2f})")
    ax.grid(alpha=0.3)
    ax.legend()
    path = os.path.join(HERE, "figures", f"cdf_{spec['name']}.png")
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)
    say(f"  figure     : {path}")


def make_say(path=None):
    """Give a print function that also appends to a log file."""
    def say(text=""):
        print(text, flush=True)
        if path:
            with open(path, "a", encoding="utf-8") as fh:
                fh.write(text + "\n")
    return say


def parse_workers(text):
    if text == "auto":
        return "auto"
    return int(text)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--cells", nargs="*", default=None,
                    help="the cells to analyse, as path:cn2:preset tokens. The "
                         "default is every finished campaign on disk.")
    ap.add_argument("--launch", choices=("collimated", "diverged"),
                    default="collimated")
    ap.add_argument("--n-trials", type=int, default=None,
                    help="the trials to read from each campaign (default all)")
    ap.add_argument("--workers", type=parse_workers, default=None,
                    help="the post-hoc process-pool size: an int, 'auto', or "
                         "none for this process")
    ap.add_argument("--figures", action="store_true",
                    help="write a coupling-loss CDF for each cell")
    ap.add_argument("--dry-run", action="store_true",
                    help="list the campaigns on disk and analyse nothing")
    args = ap.parse_args(argv)

    root = campaigns_root()
    log_path = None if args.dry_run else os.path.join(HERE, "gtilt_vs_slopes.log")
    say = make_say(log_path)
    say(f"gtilt vs slopes, {time.strftime('%Y-%m-%d %H:%M:%S')}")
    say(f"store         : {root}")
    say(f"detector      : 10 cm SMF, TipTilt() correction")
    say()

    specs = discover(root, args.launch)
    if args.cells:
        want = set(args.cells)
        specs = [s for s in specs if s["cell"] in want]
    if not specs:
        say("no finished campaign found. Check the campaigns root, or run "
            "validation/terrestrial_campaigns/run_campaigns.py first.")
        return
    say(f"found {len(specs)} campaign(s): "
        + ", ".join(s["cell"] for s in specs))
    say()
    if args.dry_run:
        for s in specs:
            say(f"  {s['cell']:<24} {s['name']:<36} "
                f"sigma_R^2 {s['sigma2_R']:.2f}, {s['grid_n']} px")
        say()
        say("dry run: nothing analysed. Drop --dry-run to run the routes.")
        return

    results = [r for r in (analyse_cell(s, args, say) for s in specs)
               if r is not None]
    out_path = os.path.join(HERE, "gtilt_vs_slopes_results.json")
    with open(out_path, "w", encoding="utf-8") as fh:
        json.dump(results, fh, indent=2)
    say(f"outputs       : {out_path}")


if __name__ == "__main__":
    main()
