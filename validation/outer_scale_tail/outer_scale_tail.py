"""The fidelity-2 SMF outer-scale AND preset fade-tail study (2-P5, 2-I3, 0-W4).

THIS STUDY NOW SERVES TWO BACKLOG ITEMS.
  1. 2-P5: the OUTER-SCALE measurement. Does the assumed outer scale L0 bias the
     fidelity-2 single-mode-fibre (SMF) fade tail of a downlink, and by how much?
  2. 2-I3: the RAPID-PRESET decision. Is the rapid preset good to use as the
     shipped default? Backlog 2-I3 flagged rapid as the probable default but
     asked to "run the catalogue before the switch". This study runs it.

THE OUTER-SCALE QUESTION (2-P5). The default screens run L0_m = inf. But three
subharmonic levels reach only about 27 times the grid side (about 95 m at 30 deg),
so the screens behave like an EFFECTIVE outer scale of 30 to 100 m, not infinity
(backlog 2-P5). A finite outer scale removes the largest turbulence cells. Those
cells carry the tip and the tilt of the wave. A fibre pays the tilt hard, so the
outer scale can move the SMF fade tail.

WHY THE OUTER SCALE MATTERS. The tail sets the link availability margin. Backlog
2-P5 says the L0 = inf default claims an outer scale a small grid cannot hold, and
it ESTIMATES the SMF p5 fade moves by "of the order of 2 dB" with the outer-scale
choice. That estimate came from the tilt share of the fade, not a measurement.
This study turns the estimate into a measured number. The measured number then
feeds 0-W4: an explicit site L0 threaded to the screens AND to the analytic tilt
Terms.

THE TWO DIMENSIONS: a CONFIG crossed with an OUTER SCALE (a 2 by 2 by default).
  CONFIG. Each config has its OWN grid and its OWN screen plan.
    ref   the well-resolved REFERENCE: the pinned standard 15-screen grid and
          plan (n = 1024 at 30 deg). Campaign preset "standard". It is the 2-P5
          deliverable.
    rapid the RAPID preset AS SHIPPED: its own coarser grid and its 5-screen
          plan, sized by turbulent_grid with PRESETS["rapid"]. Campaign preset
          "rapid". It is the 2-I3 candidate default.
  OUTER SCALE. L0 = inf (the current default) and L0 = 25 m (a finite site
    value), as a CLI list.

THE MATCHED-SEED RULE (critical). WITHIN a config the grid, the plan and the seed
are sized ONCE and shared by BOTH of that config's L0 campaigns. So within a
config only the outer-scale filter of the screens differs: the screen PLACEMENT
and the Fried radius are identical, and a per-quantile difference between the two
L0 cases is the outer scale ALONE. That is a clean matched-seed pair, and the ref
config pair is THE 2-P5 answer. ACROSS configs the grid, the pixel and the screen
count differ, so a rapid-against-ref comparison is a QUANTILE comparison, NOT a
matched-seed one; the report labels it that way.

NO CORRECTION. The fidelity-2 layer applies no tip-tilt removal and no adaptive
optics (backlog 2-AO). So the SMF coupled power here is the RAW uncorrected
atmosphere. The SMF tail is where the outer scale bites most, because the fibre
pays the full received tilt. A tracked (aligned) terminal removes the tilt, so it
cares less about the outer scale. The verdict of this study applies to the
uncorrected case only.

THE MEASURED QUANTITIES. Each trial gives one composite loss

    loss_db = -10 log10(collected_power * smf_eta)

which is the fidelity-2 downlink SMF quantity of olb.models.waveoptics. The grid
is pinned WITHIN a config, so the vacuum reference is one constant that cancels in
a comparison between the two L0 cases of that config. The study reports p50, p10,
p5 and p1 of that loss, where pX is the loss the link EXCEEDS X percent of the
time. It also reports the POINT irradiance fade

    point_loss_db = -10 log10(I / <I>)

from the CENTRE PIXEL of the stored receive-field patch. A point receiver is the
sharpest probe of the tilt there is, because nothing averages it. A point compared
ACROSS configs is pixel-limited, because the two grids average a point over a
different pixel; the report flags that.

VALIDATION ONLY. The script reads the production layer. It changes no olb module.
The scenario, the pinned grid rule and the statistics helpers are copied from
validation/tail_convergence/tail_convergence.py (the same hero downlink), with
attribution at each copy.

Sources:
- Andrews and Phillips, Laser Beam Propagation through Random Media, 2nd ed.
  (2005), DOI 10.1117/3.626196. Ch. 8, the scintillation index is the
  normalised irradiance variance. Ch. 12, Eq. (14), printed p. 482, the slant
  airmass sec(zeta). The Kolmogorov and von Karman outer-scale spectra.
- Schmidt, Numerical Simulation of Optical Wave Propagation with Examples in
  MATLAB (2010), DOI 10.1117/3.866274, Ch. 9, Eq. (9.50), printed p. 159, the
  von Karman phase power spectrum with the outer-scale term.
- Fried, DOI 10.1364/JOSA.56.001372. The Fried parameter r0.

Run it from the repository root:

    python -m validation.outer_scale_tail.outer_scale_tail
    python -m validation.outer_scale_tail.outer_scale_tail --L0 inf 25
    python -m validation.outer_scale_tail.outer_scale_tail --configs ref rapid
    python -m validation.outer_scale_tail.outer_scale_tail --configs ref
    python -m validation.outer_scale_tail.outer_scale_tail --analyse-only
"""

import argparse
import json
import os
import time
import tracemalloc
import warnings
from dataclasses import replace

import numpy as np

from olb.geometry import CircularOrbit
from olb.scenario import Channel, Site, SpaceScenario
from olb.terminal import SMF, Terminal, Transmitter
from olb.turbulence.andrews.paths import sec_zeta
from olb.waveoptics.turbulence.campaign import Campaign
from olb.waveoptics.turbulence.sampling import PRESETS, turbulent_grid

HERE = os.path.dirname(os.path.abspath(__file__))
CAMPAIGN_ROOT = os.environ.get("OLB_OUTER_SCALE_ROOT",
                               os.path.join(HERE, "campaigns"))

# The wavelength and the base seed. The seed is SHARED across the L0 cases of one
# config, so the two atmospheres of a config differ only in the outer-scale
# filter.
LAM = 1550e-9
SEED = 20260904

# The pinned grid of the REFERENCE config comes from the sizer of this screen
# count. 15 screens gives the finer grid (n = 1024 at 30 deg). Both L0 cases of
# the ref config run on this ONE grid and this ONE 15-screen plan. Source: the
# pinned-grid rule of validation/tail_convergence/tail_convergence.py (PIN_N).
PIN_N = 15

# The config names this study knows. "ref" is the well-resolved reference (2-P5),
# "rapid" is the shipped rapid preset as a default candidate (2-I3).
ALL_CONFIGS = ("ref", "rapid")

# The exceedance probabilities. pX is the loss that the link EXCEEDS X percent
# of the time, so it is the (100 - X) percentile of the loss.
EXCEEDANCE = (0.50, 0.10, 0.05, 0.01)

# The bootstrap of every quantile: the resample count, the interval and the
# fixed rng seed. Source: validation/tail_convergence/tail_convergence.py.
N_BOOT = 1000
BOOT_INTERVAL = 0.68
BOOT_SEED = 7

# The trial counts of the growth check (how a tail estimate settles).
GROWTH_N = (100, 200, 400, 600, 800, 1000)

# The subharmonic factor: three subharmonic levels reach about 27 times the grid
# side. Source: backlog 2-P5 (three levels, 3^3 = 27).
SUBHARMONIC_FACTOR = 27


# ---------------------------------------------------------------------------
# The scenario (copied verbatim from tail_convergence.py; same hero downlink)
# ---------------------------------------------------------------------------

def scenario_and_geometry(elevation_deg):
    """Build the downlink scenario and the orbit of one elevation.

    The case is the presentation downlink: a 1550 nm space-to-ground link to a
    700 mm ground telescope with a single-mode-fibre receiver, from a 100 mm
    space terminal at 500 km.

    Source: validation/tail_convergence/tail_convergence.py, copied verbatim.

    Args:
        elevation_deg: the elevation of the line of sight, in deg.

    Returns:
        The pair (SpaceScenario, CircularOrbit).
    """
    site = Site(cn2_ground=1.7e-14, wind_rms_m_s=21.0)
    channel = Channel(site=site, altitude_m=500e3)
    ground = Terminal(aperture_m=0.7, wavelength_m=LAM,
                      pointing_jitter_rad=2e-6,
                      detector=SMF(sensitivity_dbm=-45.0))
    space = Terminal(aperture_m=0.10, wavelength_m=LAM,
                     pointing_jitter_rad=1e-6,
                     transmitter=Transmitter(waist_m=0.04, power_dbm=30.0))
    scn = SpaceScenario(ground=ground, space=space, direction="downlink",
                        channel=channel)
    geom = CircularOrbit(altitude_m=500e3, elevation_deg=[float(elevation_deg)])
    return scn, geom


def shared_grid_plan(scn, geom):
    """Size the pinned grid and plan of the REFERENCE config.

    The function pins the grid the same way the tail-convergence study does: it
    sizes the standard preset with a screen floor of PIN_N, which gives the
    finer grid. It returns the grid AND the plan. Both L0 campaigns of the ref
    config take this SAME grid and this SAME plan, so only the outer-scale filter
    of the screens differs. The sizer runs with L0_m = inf; the outer scale does
    not move the screen PLACEMENT or the Fried radius, only the screen spectrum.

    Source: the _sizer helper of
    validation/tail_convergence/tail_convergence.py.

    Args:
        scn:  the SpaceScenario.
        geom: the CircularOrbit.

    Returns:
        The triple (GridSpec, ScreenPlan, the sorted warning texts).
    """
    preset = replace(PRESETS["standard"], name=f"std{int(PIN_N)}",
                     min_screens=int(PIN_N))
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        grid, plan, _ = turbulent_grid(scn, geom, preset=preset)
    return grid, plan, sorted({str(w.message) for w in caught})


def rapid_grid_plan(scn, geom):
    """Size the grid and plan of the RAPID config AS SHIPPED.

    The function sizes the rapid preset end to end, so the grid AND the 5-screen
    plan are exactly what a preset="rapid" user gets. Both L0 campaigns of the
    rapid config take this SAME grid and this SAME plan.

    Args:
        scn:  the SpaceScenario.
        geom: the CircularOrbit.

    Returns:
        The triple (GridSpec, ScreenPlan, the sorted warning texts).
    """
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        grid, plan, _ = turbulent_grid(scn, geom, preset=PRESETS["rapid"])
    return grid, plan, sorted({str(w.message) for w in caught})


def build_configs(scn, geom, names):
    """Build the (grid, plan, preset, doc, warnings) of each named config.

    Each config sizes its OWN grid and its OWN screen plan ONCE. The (grid, plan)
    pair is then shared by every L0 campaign of that config, so within a config
    only the outer-scale filter changes.

    Args:
        scn:   the SpaceScenario.
        geom:  the CircularOrbit.
        names: the config names, from ALL_CONFIGS.

    Returns:
        A dict name -> {"grid", "plan", "preset", "doc", "warnings"}.
    """
    out = {}
    for name in names:
        if name == "ref":
            grid, plan, warns = shared_grid_plan(scn, geom)
            out[name] = {
                "grid": grid, "plan": plan, "preset": "standard",
                "warnings": warns,
                "doc": "the well-resolved REFERENCE: the pinned standard "
                       f"{PIN_N}-screen grid and plan (the 2-P5 deliverable).",
            }
        elif name == "rapid":
            grid, plan, warns = rapid_grid_plan(scn, geom)
            out[name] = {
                "grid": grid, "plan": plan, "preset": "rapid",
                "warnings": warns,
                "doc": "the RAPID preset AS SHIPPED: its own coarser grid and "
                       "5-screen plan (the 2-I3 default candidate).",
            }
        else:
            raise ValueError(f"unknown config {name!r}. Use one of "
                             f"{ALL_CONFIGS}.")
    return out


# ---------------------------------------------------------------------------
# The L0 command-line values
# ---------------------------------------------------------------------------

def parse_l0(token):
    """Turn a command-line token into an outer-scale value, in m.

    The tokens "inf", "np.inf" and "infinity" all give numpy inf. Any other
    token is read as a float.

    Args:
        token: one command-line string.

    Returns:
        The outer scale as a float, in m.
    """
    low = str(token).strip().lower()
    if low in ("inf", "np.inf", "infinity", "infinite"):
        return float(np.inf)
    return float(low)


def l0_label(value):
    """Give a short, file-safe label for an outer-scale value."""
    return "inf" if not np.isfinite(value) else f"{value:g}m"


def l0_text(value):
    """Give a short, human label for an outer-scale value."""
    return "inf" if not np.isfinite(value) else f"{value:g} m"


# ---------------------------------------------------------------------------
# The statistics (copied from tail_convergence.py; small pure helpers)
# ---------------------------------------------------------------------------

def _quantiles(loss_db, probs=EXCEEDANCE):
    """Give the loss that the link EXCEEDS each probability.

    A fade is a LARGE loss, so the loss exceeded a fraction q of the time is the
    (1 - q) quantile of the loss sample.

    Source: validation/tail_convergence/tail_convergence.py.

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

    The function resamples the trials with replacement, it takes the statistic
    of each resample, and it gives the half-width of the central BOOT_INTERVAL
    band. The rng seed is FIXED, so a rerun of the analysis gives the same
    interval.

    Source: validation/tail_convergence/tail_convergence.py.

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
    and Phillips, DOI 10.1117/3.626196, Ch. 8 (the scintillation index is the
    normalised irradiance variance). olb holds the definition in
    olb.turbulence.andrews.scintillation. Copied from
    validation/tail_convergence/tail_convergence.py.
    """
    a = np.asarray(x, dtype=float)
    return float(a.var() / a.mean() ** 2)


def _centre_irradiance(result):
    """Give the irradiance of the CENTRE PIXEL of each stored trial.

    The centre pixel is the flat index (n // 2) * n + (n // 2), the axis pixel
    of the CircAperture convention. Source:
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
    """Give the total byte count of the files in one directory.

    Source: validation/tail_convergence/tail_convergence.py.
    """
    return int(sum(os.path.getsize(os.path.join(path, f))
                   for f in os.listdir(path)
                   if os.path.isfile(os.path.join(path, f))))


# ---------------------------------------------------------------------------
# The analysis of one (config, L0) case
# ---------------------------------------------------------------------------

def analyse(camp, config, label, doc, n_trials, growth=False):
    """Measure one (config, L0) campaign.

    The body follows the analyse helper of
    validation/tail_convergence/tail_convergence.py: it loads the stored
    fields, it forms the composite SMF loss and the point fade, and it gives the
    quantiles with a bootstrap.

    Args:
        camp:     the Campaign.
        config:   the config name ("ref" or "rapid").
        label:    the file-safe L0 label.
        doc:      the one-line human text of the case.
        n_trials: the number of trials to read.
        growth:   True adds the growth-with-n table.

    Returns:
        The triple (a result dict, the SMF loss array, the point loss array).
    """
    tracemalloc.start()
    result = camp.load(n_trials, fields=True)
    fields_bytes = int(result.fields.nbytes)
    peak_bytes = int(tracemalloc.get_traced_memory()[1])
    tracemalloc.stop()

    power = np.array([t.collected_power for t in result.trials], dtype=float)
    eta = np.array([t.smf_eta for t in result.trials], dtype=float)
    # The composite fidelity-2 downlink SMF loss. The vacuum reference is one
    # constant for the pinned grid of a config, so it cancels between the two L0
    # cases of that config. See olb.models.waveoptics.
    loss = -10.0 * np.log10(power * eta)

    q = _quantiles(loss)
    q_half = _bootstrap(loss, _quantiles)
    point = _centre_irradiance(result)
    # The point fade, against the mean irradiance. Nothing averages a point, so
    # this is the sharpest tail of the study.
    point_loss = -10.0 * np.log10(point / point.mean())
    pq = _quantiles(point_loss)
    pq_half = np.atleast_1d(_bootstrap(point_loss, _quantiles))

    def _by_p(values):
        return {f"p{int(p * 100)}": float(v)
                for p, v in zip(EXCEEDANCE, np.atleast_1d(values))}

    out = {
        "config": config,
        "case": label,
        "doc": doc,
        "L0_m": (None if not np.isfinite(camp.L0_m) else float(camp.L0_m)),
        "preset": camp.preset,
        "n_trials": int(loss.size),
        "n_screens": int(camp.plan.z_m.size),
        "grid_n": int(camp.grid.n),
        "grid_size_m": float(camp.grid.size_m),
        "grid_pixel_m": float(camp.grid.pixel_m),
        "subharmonic_reach_m": float(SUBHARMONIC_FACTOR * camp.grid.size_m),
        "bottom_cn2_share": float(camp.plan.cn2_int_m13[-1]
                                  / camp.plan.cn2_int_m13.sum()),
        "r0_total_m": float(camp.plan.r0_total_m),
        "mean_db": float(loss.mean()),
        "mean_db_half": _bootstrap(loss, lambda s: float(s.mean())),
        "quantiles_db": {f"p{int(p * 100)}": float(v)
                         for p, v in zip(EXCEEDANCE, q)},
        "quantiles_half_db": {f"p{int(p * 100)}": float(v)
                              for p, v in zip(EXCEEDANCE, np.atleast_1d(q_half))},
        "fade_depth_db": {f"p{int(p * 100)}": float(v - q[0])
                          for p, v in zip(EXCEEDANCE, q)},
        "n_past_p1": int((loss > q[-1]).sum()),
        "point_quantiles_db": _by_p(pq),
        "point_quantiles_half_db": _by_p(pq_half),
        "sigma2_P": _index(power),
        "sigma2_P_half": _bootstrap(power, _index),
        "sigma2_I_point": _index(point),
        "sigma2_I_point_half": _bootstrap(point, _index),
        "fields_bytes": fields_bytes,
        "load_peak_bytes": peak_bytes,
        "disk_bytes": _dir_bytes(camp.root_dir),
        "n_stored": int(camp.n_stored),
        "mean_wall_time_s": float(np.mean([t.wall_time_s
                                           for t in result.trials])),
    }
    if growth:
        rows = []
        for n in GROWTH_N:
            if n > loss.size:
                continue
            gq = _quantiles(loss[:n])
            gh = np.atleast_1d(_bootstrap(loss[:n], _quantiles))
            rows.append({"n": int(n),
                         "quantiles_db": {f"p{int(p * 100)}": float(v)
                                          for p, v in zip(EXCEEDANCE, gq)},
                         "half_db": {f"p{int(p * 100)}": float(v)
                                     for p, v in zip(EXCEEDANCE, gh)}})
        out["growth"] = rows
    return out, loss, point_loss


# ---------------------------------------------------------------------------
# The figures
# ---------------------------------------------------------------------------

def _plot(rows, growth_rows, survival, ref_title, config_order, fig_dir):
    """Draw the figures. Give the list of the written paths.

    Args:
        rows:         the per-case result dicts.
        growth_rows:  the growth-with-n rows of the reference case.
        survival:     a dict of two survival panels; each maps a
                      "config/L0" label to a loss array.
        ref_title:    the human label of the reference (config, L0) case.
        config_order: the config names, in report order.
        fig_dir:      the output directory.

    Returns:
        The list of the written paths.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib is not installed; the figures are skipped.")
        return []

    os.makedirs(fig_dir, exist_ok=True)
    written = []
    names = [f"p{int(p * 100)}" for p in EXCEEDANCE]
    cfgs = [c for c in config_order if any(r["config"] == c for r in rows)]

    # ---- 1. the empirical survival functions: one curve per (config, L0) ----
    fig, axes = plt.subplots(1, 2, figsize=(12.0, 5.0))
    for ax, (title, curves) in zip(axes, survival.items()):
        for name, loss in curves.items():
            s = np.sort(np.asarray(loss, dtype=float))
            frac = 1.0 - np.arange(s.size) / s.size
            ax.semilogy(s, frac, label=name)
        ax.set_xlabel("loss [dB]")
        ax.set_ylabel("fraction of trials with a larger loss")
        ax.set_title(title)
        ax.grid(alpha=0.3, which="both")
        ax.legend(fontsize=8)
    fig.tight_layout()
    path = os.path.join(fig_dir, "survival.png")
    fig.savefig(path, dpi=140)
    plt.close(fig)
    written.append(path)

    # ---- 2. the SMF and the point quantiles against L0, one series per config-
    # The x axis is the outer scale in metres, with L0 = inf drawn at the right.
    def _x_of(r):
        return 1e4 if r["L0_m"] is None else float(r["L0_m"])

    fig, axes = plt.subplots(1, 2, figsize=(12.0, 5.0))
    for ax, (title, qkey, hkey) in zip(
            axes,
            (("The SMF loss quantiles", "quantiles_db", "quantiles_half_db"),
             ("The point fade quantiles", "point_quantiles_db",
              "point_quantiles_half_db"))):
        for cfg in cfgs:
            cfg_rows = sorted([r for r in rows if r["config"] == cfg],
                              key=_x_of)
            if not cfg_rows:
                continue
            for key in names:
                ax.errorbar([_x_of(r) for r in cfg_rows],
                            [r[qkey][key] for r in cfg_rows],
                            yerr=[r[hkey][key] for r in cfg_rows],
                            marker="o", capsize=3, label=f"{cfg} {key}")
        ax.set_xscale("log")
        ax.set_xlabel("outer scale L0 [m] (inf drawn at 1e4)")
        ax.set_ylabel("loss [dB]")
        ax.set_title(title)
        ax.grid(alpha=0.3, which="both")
        ax.legend(fontsize=7, ncol=len(cfgs))
    fig.tight_layout()
    path = os.path.join(fig_dir, "quantiles_vs_L0.png")
    fig.savefig(path, dpi=140)
    plt.close(fig)
    written.append(path)

    # ---- 3. the growth with the trial count (the reference case) ----
    if growth_rows:
        fig, ax = plt.subplots(figsize=(7.5, 5.0))
        for key in names:
            ax.errorbar([g["n"] for g in growth_rows],
                        [g["quantiles_db"][key] for g in growth_rows],
                        yerr=[g["half_db"][key] for g in growth_rows],
                        marker="o", capsize=3, label=key)
        ax.set_xlabel("trials")
        ax.set_ylabel("loss [dB]")
        ax.set_title(f"How the tail estimate settles ({ref_title})")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
        fig.tight_layout()
        path = os.path.join(fig_dir, "growth_with_trials.png")
        fig.savefig(path, dpi=140)
        plt.close(fig)
        written.append(path)
    return written


# ---------------------------------------------------------------------------
# The driver
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--elevation", type=float, default=30.0,
                    help="the elevation of the line of sight [deg]")
    ap.add_argument("--n-trials", type=int, default=1000,
                    help="the trials of each (config, L0) case")
    ap.add_argument("--workers", type=int, default=8,
                    help="the processes of the campaign pool (8, not 16: a "
                         "16-worker 1024 px run has run out of memory)")
    ap.add_argument("--block-size", type=int, default=50,
                    help="the trials in one block file")
    ap.add_argument("--L0", nargs="+", default=["inf", "25"],
                    help="the outer-scale values [m]; 'inf' gives infinity")
    ap.add_argument("--configs", nargs="+", default=list(ALL_CONFIGS),
                    help=f"the configs to run, from {list(ALL_CONFIGS)}")
    ap.add_argument("--analyse-only", action="store_true",
                    help="skip the runs and read what is stored")
    args = ap.parse_args()

    el = float(args.elevation)
    tag = f"el{el:02.0f}"

    l0_values = [parse_l0(t) for t in args.L0]
    labels = [l0_label(v) for v in l0_values]
    if len(set(labels)) != len(labels):
        raise ValueError(f"the L0 values give a repeated label: {labels}")

    config_names = list(args.configs)
    for c in config_names:
        if c not in ALL_CONFIGS:
            raise ValueError(f"unknown config {c!r}. Use one of {ALL_CONFIGS}.")
    if len(set(config_names)) != len(config_names):
        raise ValueError(f"a config is repeated: {config_names}")

    text_of = {lab: l0_text(v) for lab, v in zip(labels, l0_values)}
    l0_of = {lab: v for lab, v in zip(labels, l0_values)}

    log_path = os.path.join(HERE, f"outer_scale_tail_{tag}.log")
    # The log is written LINE BY LINE, so a killed run keeps its log.
    with open(log_path, "w", encoding="utf-8"):
        pass

    def say(text=""):
        print(text, flush=True)
        with open(log_path, "a", encoding="utf-8") as fh:
            fh.write(text + "\n")

    scn, geom = scenario_and_geometry(el)
    configs = build_configs(scn, geom, config_names)

    say("The fidelity-2 SMF outer-scale AND preset fade-tail study "
        "(2-P5, 2-I3, 0-W4)")
    say(f"elevation     : {el:.1f} deg      trials/case: {args.n_trials}")
    say(f"seed          : {SEED}      block size: {args.block_size}      "
        f"workers: {args.workers}")
    say(f"scenario      : downlink, 1550 nm, 500 km, 700 mm ground SMF, "
        f"100 mm space terminal")
    say(f"configs       : " + ", ".join(config_names))
    say(f"L0 values [m] : " + ", ".join(l0_text(v) for v in l0_values))
    say(f"cases         : {len(config_names)} configs x {len(labels)} outer "
        f"scales = {len(config_names) * len(labels)} campaigns")
    say(f"mode          : {'ANALYSE ONLY' if args.analyse_only else 'RUN'}")
    say("method        : WITHIN a config, a MATCHED-SEED L0 pair (one grid, one "
        "plan, one seed; only the outer-scale filter changes).")
    say("              : ACROSS configs, a QUANTILE comparison (the grids, the "
        "pixels and the screen counts differ), NOT matched-seed.")
    say("correction    : NONE. Fidelity 2 applies no tip-tilt removal and no "
        "AO (backlog 2-AO); the SMF power is the raw atmosphere.")
    say("primary       : the SMF coupled-power fade and the POINT irradiance "
        "fade. The aperture index s2_P is a footnote (strongly averaged).")
    if args.n_trials < 1000:
        say("WARNING: under 1000 trials, so the p1 tail is UNDER-SAMPLED. Read "
            "the p10 and p5 rows only.")
    say()

    # ---- the grid and the plan of each config ----
    say("THE CONFIGURATIONS (each config has its OWN grid and screen plan)")
    elev_sec = float(sec_zeta(el))
    for cfg in config_names:
        spec = configs[cfg]
        grid, plan = spec["grid"], spec["plan"]
        heights = (plan.z_total_m - plan.z_m) / elev_sec
        share = plan.cn2_int_m13 / plan.cn2_int_m13.sum()
        reach = SUBHARMONIC_FACTOR * grid.size_m
        say(f"  {cfg} ({spec['preset']} preset): {spec['doc']}")
        say(f"    grid {grid.n} px, {grid.size_m:.3f} m "
            f"({grid.pixel_m * 1e3:.2f} mm px); {plan.z_m.size} screens; "
            f"r0_total {plan.r0_total_m * 100:.2f} cm")
        say(f"    screen h [m] : " + " ".join(f"{v:.0f}" for v in heights))
        say(f"    Cn2 share    : " + " ".join(f"{v:.3f}" for v in share))
        # The grid side sets the largest turbulence cell the screens can hold.
        say(f"    subharmonic reach ~ {reach:.0f} m ({SUBHARMONIC_FACTOR} x the "
            f"{grid.size_m:.1f} m grid side; backlog 2-P5)")
        for w in spec["warnings"]:
            say(f"    sizer warning: {w}")
        say()
    if "ref" in config_names and "rapid" in config_names:
        r_reach = SUBHARMONIC_FACTOR * configs["ref"]["grid"].size_m
        p_reach = SUBHARMONIC_FACTOR * configs["rapid"]["grid"].size_m
        smaller = "SMALLER than" if p_reach < r_reach else "NOT smaller than"
        say(f"  NOTE: the rapid subharmonic reach ({p_reach:.0f} m) is "
            f"{smaller} the reference ({r_reach:.0f} m), so rapid holds a "
            f"smaller effective outer scale (backlog 2-P5).")
        say()

    # ---- one campaign for each (config, L0), on that config's grid and plan ---
    doc_of = {}
    camps = {}
    for cfg in config_names:
        spec = configs[cfg]
        for lab in labels:
            key = (cfg, lab)
            doc_of[key] = (f"config {cfg}, matched-seed L0 = {text_of[lab]}")
            root = os.path.join(CAMPAIGN_ROOT, tag, cfg, lab)
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                camp = Campaign(scn, geom, root, seed=SEED,
                                preset=spec["preset"],
                                block_size=args.block_size,
                                grid=spec["grid"], plan=spec["plan"],
                                L0_m=l0_of[lab])
            camps[key] = camp
            for w in sorted({str(x.message) for x in caught}):
                say(f"  {cfg}/{lab}: campaign warning: {w}")
    say()

    # The reference (config, L0) case: it carries the growth check and is the
    # 2-P5 deliverable. The reference is the ref config at L0 = inf when present.
    ref_config = "ref" if "ref" in config_names else config_names[0]
    ref_l0 = "inf" if "inf" in labels else labels[0]
    ref_key = (ref_config, ref_l0)
    ref_title = f"{ref_config}, L0 = {text_of[ref_l0]}"

    # ---- the runs ----
    timing = {}
    if not args.analyse_only:
        for cfg in config_names:
            for lab in labels:
                key = (cfg, lab)
                camp = camps[key]
                say(f"RUN {cfg}/{lab}: {args.n_trials} trials")
                t0 = time.perf_counter()
                n_done = camp.run(args.n_trials, workers=args.workers,
                                  progress=True)
                wall = time.perf_counter() - t0
                timing[f"{cfg}/{lab}"] = {"wall_s": float(wall),
                                          "n_stored": int(n_done)}
                per = wall / max(n_done, 1)
                say(f"  {cfg}/{lab}: {n_done} trials on disk, {wall:.1f} s "
                    f"wall, {per:.2f} s/trial (this call), "
                    f"{_dir_bytes(camp.root_dir) / 1e6:.1f} MB on disk")
                # A CHECKPOINT: one partial result as soon as the case is done,
                # so a long run is not a black box. The scalars load in a moment.
                quick = camp.load(n_done, fields=False)
                ql = -10.0 * np.log10(np.array(
                    [t.collected_power * t.smf_eta for t in quick.trials]))
                qq = _quantiles(ql)
                say(f"  CHECKPOINT {cfg}/{lab}: SMF loss mean {ql.mean():.2f} "
                    f"dB, p50 {qq[0]:.2f}, p10 {qq[1]:.2f}, p5 {qq[2]:.2f}, "
                    f"p1 {qq[3]:.2f} dB (n = {ql.size})")
                say()

    # ---- the analysis ----
    rows, growth_rows = [], []
    row_of = {}
    survival = {"the SMF loss": {}, "the point fade": {}}
    for cfg in config_names:
        for lab in labels:
            key = (cfg, lab)
            camp = camps[key]
            if camp.n_stored == 0:
                say(f"  {cfg}/{lab}: no stored trial. Skipped.")
                continue
            n = min(args.n_trials, camp.n_stored)
            want_growth = key == ref_key
            row, loss, point_loss = analyse(camp, cfg, lab, doc_of[key], n,
                                            growth=want_growth)
            row["run_wall_s"] = timing.get(f"{cfg}/{lab}", {}).get("wall_s")
            if want_growth and row.get("growth"):
                growth_rows = row["growth"]
                row["growth_case"] = True
            rows.append(row)
            row_of[key] = row
            curve = f"{cfg}/{lab}"
            survival["the SMF loss"][curve] = loss
            survival["the point fade"][curve] = point_loss

    # ---- the per-case quantile table ----
    say("THE PER-CASE QUANTILE TABLE (pX = the loss EXCEEDED X percent "
        "of the time)")
    hdr = (f"{'cfg':>7s}{'L0':>6s}{'scr':>4s}{'grid':>6s}{'mean':>8s}"
           f"{'p50':>8s}{'p10':>15s}{'p5':>15s}{'p1':>15s}{'d_p5':>7s}"
           f"{'d_p1':>7s}{'pt_p5':>14s}{'pt_p1':>14s}{'s2_I':>8s}{'s2_P':>8s}"
           f"{'s/tr':>7s}{'MB':>7s}")
    say(hdr)
    say("-" * len(hdr))
    for r in rows:
        q, h = r["quantiles_db"], r["quantiles_half_db"]
        pq, ph = r["point_quantiles_db"], r["point_quantiles_half_db"]
        say(f"{r['config']:>7s}{r['case']:>6s}{r['n_screens']:4d}"
            f"{r['grid_n']:6d}{r['mean_db']:8.2f}{q['p50']:8.2f}"
            f"{q['p10']:9.2f}+-{h['p10']:4.2f}"
            f"{q['p5']:9.2f}+-{h['p5']:4.2f}"
            f"{q['p1']:9.2f}+-{h['p1']:4.2f}"
            f"{r['fade_depth_db']['p5']:7.2f}{r['fade_depth_db']['p1']:7.2f}"
            f"{pq['p5']:8.2f}+-{ph['p5']:4.2f}"
            f"{pq['p1']:8.2f}+-{ph['p1']:4.2f}"
            f"{r['sigma2_I_point']:8.3f}{r['sigma2_P']:8.3f}"
            f"{r['mean_wall_time_s']:7.2f}{r['disk_bytes'] / 1e6:7.1f}")
    say("  cfg = the config (ref = well-resolved reference, rapid = shipped "
        "rapid preset). L0 = the outer scale [m] ('inf' is infinite).")
    say("  scr = screens. p10, p5, p1 = the SMF loss EXCEEDED 10, 5, 1 percent "
        "of the time. d_p5 and d_p1 are the FADE DEPTH pX - p50 [dB].")
    say("  pt_p5, pt_p1 = the POINT irradiance fade -10 log10(I/<I>) exceeded "
        "5 and 1 percent of the time [dB]. s2_I = the point index.")
    say("  s2_P = the aperture index of the 700 mm bucket. It is strongly "
        "averaged and INSENSITIVE: read it as a footnote only.")
    say("  s/tr = the mean stored wall time of one trial [s]. The +- values "
        f"are {BOOT_INTERVAL * 100:.0f} % bootstrap half-widths "
        f"({N_BOOT} resamples).")
    say()

    # ---- the growth check ----
    if growth_rows:
        say(f"HOW THE TAIL ESTIMATE SETTLES WITH THE TRIAL COUNT ({ref_title})")
        say(f"{'n':>6s}{'p50':>9s}{'p10':>15s}{'p5':>15s}{'p1':>15s}")
        for g in growth_rows:
            q, h = g["quantiles_db"], g["half_db"]
            say(f"{g['n']:6d}{q['p50']:9.2f}"
                f"{q['p10']:9.2f}+-{h['p10']:4.2f}"
                f"{q['p5']:9.2f}+-{h['p5']:4.2f}"
                f"{q['p1']:9.2f}+-{h['p1']:4.2f}")
        say()

    # ---- the verdict, in three labelled parts ----
    def _delta(a, b, qkey, hkey, key):
        """Give (b - a, the combined bootstrap bar, the sigma count)."""
        d = b[qkey][key] - a[qkey][key]
        bar = float(np.hypot(a[hkey][key], b[hkey][key]))
        return float(d), bar, (abs(d) / bar if bar > 0 else float("inf"))

    # The two primary quantities, in the same words.
    faces = (("SMF", "quantiles_db", "quantiles_half_db"),
             ("POINT", "point_quantiles_db", "point_quantiles_half_db"))
    verdict = {"outer_scale_within_config": {},
               "rapid_vs_ref_matched_L0": {},
               "bottom_line": {}}

    say("THE VERDICT")
    say()

    # ---- (A) the OUTER-SCALE effect, WITHIN each config (matched-seed) ----
    say("(A) OUTER-SCALE EFFECT, WITHIN each config (matched-seed, clean)")
    say("    delta = loss(this L0) - loss(the reference L0 = inf). A NEGATIVE "
        "delta means the smaller outer scale gives LESS fade (less loss).")
    if "inf" not in labels:
        say("    no L0 = inf reference in the L0 list; part (A) is skipped.")
    else:
        for cfg in config_names:
            ref_row = row_of.get((cfg, "inf"))
            if ref_row is None:
                say(f"  {cfg}: no stored L0 = inf case; skipped.")
                continue
            tag_2p5 = " (THE 2-P5 DELIVERABLE)" if cfg == "ref" else ""
            say(f"  config {cfg}{tag_2p5}:")
            cfg_block = {}
            any_other = False
            for lab in labels:
                if lab == "inf":
                    continue
                other = row_of.get((cfg, lab))
                if other is None:
                    continue
                any_other = True
                pair_block = {}
                for face, qkey, hkey in faces:
                    for k in ("p5", "p1"):
                        d, bar, sig = _delta(ref_row, other, qkey, hkey, k)
                        moves = sig > 1.0
                        pair_block[f"{face}_{k}"] = {
                            "delta_db": d, "bar_db": bar, "sigmas": float(sig),
                            "resolved": bool(moves)}
                        say(f"    {face} {k}: L0 inf -> {text_of[lab]} moves the "
                            f"fade by {d:+.2f} dB, against a combined bootstrap "
                            f"bar of {bar:.2f} dB ({sig:.1f} sigma). "
                            f"{'IT MOVES' if moves else 'NO RESOLVED MOVE'}.")
                cfg_block[lab] = pair_block
            if not any_other:
                say("    only the L0 = inf case is stored; nothing to compare.")
            verdict["outer_scale_within_config"][cfg] = cfg_block
    say()

    # ---- (B) RAPID vs REFERENCE at matched L0 (quantile comparison) ----
    say("(B) RAPID vs REFERENCE at matched L0 (QUANTILE comparison, NOT "
        "matched-seed)")
    say("    delta = loss(rapid) - loss(ref) at the same outer scale. The two "
        "grids and screen counts differ, so this is a quantile comparison.")
    if "ref" not in config_names or "rapid" not in config_names:
        say("    part (B) needs BOTH the ref and the rapid configs; skipped.")
    else:
        for lab in labels:
            r_ref = row_of.get(("ref", lab))
            r_rap = row_of.get(("rapid", lab))
            if r_ref is None or r_rap is None:
                say(f"    L0 = {text_of[lab]}: a config has no stored trial; "
                    f"skipped.")
                continue
            block = {}
            for face, qkey, hkey in faces:
                for k in ("p5", "p1"):
                    d, bar, sig = _delta(r_ref, r_rap, qkey, hkey, k)
                    within = abs(d) <= bar
                    block[f"{face}_{k}"] = {
                        "delta_db": d, "bar_db": bar, "sigmas": float(sig),
                        "within_bar": bool(within)}
                    # A POINT compared across configs is pixel-limited: the two
                    # grids average a point over a different pixel. See the
                    # grid-effect note of tail_convergence.py.
                    note = ""
                    if face == "POINT":
                        note = (" PIXEL-LIMITED: the grids average a point over "
                                f"a different pixel "
                                f"({r_rap['grid_pixel_m'] * 1e3:.2f} against "
                                f"{r_ref['grid_pixel_m'] * 1e3:.2f} mm), so this "
                                "is NOT a physics difference.")
                    say(f"    {face} {k} at L0 = {text_of[lab]}: rapid - ref = "
                        f"{d:+.2f} dB against a bar of {bar:.2f} dB "
                        f"({sig:.1f} sigma).{note}")
            verdict["rapid_vs_ref_matched_L0"][lab] = block
    say()

    # ---- (C) the BOTTOM LINE for the rapid-as-default decision ----
    say("(C) BOTTOM LINE for the rapid-as-default decision")
    rule = ("rapid is supported when |rapid - ref| at SMF p5 is within the "
            "combined bootstrap bar at EVERY tested outer scale")
    say(f"    RULE: {rule}.")
    if "ref" not in config_names or "rapid" not in config_names:
        say("    the decision needs BOTH the ref and the rapid configs; "
            "skipped.")
        verdict["bottom_line"] = {"rule": rule, "per_L0": {},
                                  "rapid_supported": None}
    else:
        per_l0 = {}
        supported = True
        parts = []
        for lab in labels:
            r_ref = row_of.get(("ref", lab))
            r_rap = row_of.get(("rapid", lab))
            if r_ref is None or r_rap is None:
                supported = None
                continue
            d, bar, sig = _delta(r_ref, r_rap, "quantiles_db",
                                 "quantiles_half_db", "p5")
            within = abs(d) <= bar
            per_l0[lab] = {"smf_p5_delta_db": d, "bar_db": bar,
                           "sigmas": float(sig), "within_bar": bool(within)}
            if supported is not None and not within:
                supported = False
            parts.append(f"at L0 = {text_of[lab]} by {d:+.2f} dB")
            say(f"    At L0 = {text_of[lab]}, rapid differs from the reference "
                f"by {d:+.2f} dB at SMF p5 (combined bootstrap bar {bar:.2f} "
                f"dB); rapid {'IS' if within else 'IS NOT'} inside the "
                f"reference spread.")
        if supported is None:
            say("    a case has no stored trial; no rapid-as-default verdict.")
            word = "UNDETERMINED"
        else:
            word = "IS" if supported else "IS NOT"
            joined = "; ".join(parts) if parts else "no case"
            say(f"    Rapid differs from the reference at SMF p5 {joined}. "
                f"Rapid {word} inside the reference bootstrap spread at every "
                f"tested outer scale, so rapid {word} supported as the shipped "
                f"default on this scenario.")
        verdict["bottom_line"] = {"rule": rule, "per_L0": per_l0,
                                  "rapid_supported": (None if supported is None
                                                      else bool(supported))}
    say()

    stamp = {
        "study": "outer_scale_tail",
        "backlog": ["2-P5", "2-I3", "0-W4"],
        "elevation_deg": el, "n_trials": args.n_trials, "seed": SEED,
        "block_size": args.block_size, "workers": args.workers,
        "wavelength_m": LAM,
        "configs": {
            cfg: {"preset": configs[cfg]["preset"],
                  "grid": {"n": int(configs[cfg]["grid"].n),
                           "size_m": float(configs[cfg]["grid"].size_m),
                           "pixel_m": float(configs[cfg]["grid"].pixel_m)},
                  "n_screens": int(configs[cfg]["plan"].z_m.size),
                  "r0_total_m": float(configs[cfg]["plan"].r0_total_m),
                  "subharmonic_reach_m": float(
                      SUBHARMONIC_FACTOR * configs[cfg]["grid"].size_m)}
            for cfg in config_names},
        "L0_values_m": [None if not np.isfinite(v) else float(v)
                        for v in l0_values],
        "reference_config": ref_config, "reference_L0": ref_l0,
        "exceedance": list(EXCEEDANCE), "n_bootstrap": N_BOOT,
        "bootstrap_interval": BOOT_INTERVAL,
        "analyse_only": bool(args.analyse_only),
        "timing": timing, "verdict": verdict, "cases": rows,
    }
    json_path = os.path.join(HERE, f"outer_scale_tail_{tag}_results.json")
    with open(json_path, "w", encoding="utf-8") as fh:
        json.dump(stamp, fh, indent=2)
    figs = _plot(rows, growth_rows, survival, ref_title, config_names,
                 os.path.join(HERE, "figures"))
    print(f"\nwrote {json_path}")
    for f in figs:
        print(f"wrote {f}")
    print(f"wrote {log_path}")


if __name__ == '__main__':
    main()
