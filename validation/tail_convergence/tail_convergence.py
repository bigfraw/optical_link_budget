"""The fidelity-2 SMF fade-tail convergence study (backlog 2-I2T and 2-N6).

THE QUESTION. Does the deep single-mode-fibre (SMF) fade tail of a fidelity-2
downlink CONVERGE as the near-ground Cn2 is resolved with more, thinner phase
screens? The MEAN is already validated flat (docs/schmidt-crosscheck.md, WP7).
The TAIL is the open risk, because the tail sets the link availability margin.
A post-WP7 matched-seed re-test saw about 2 dB of p5 movement between two screen
placements, but at 200 trials that movement stayed inside the Monte-Carlo noise
(docs/schmidt-crosscheck.md, the post-WP7 measurement). This study resolves it.

THE METHOD. The study PINS the grid and it moves the screens only. The shipped
sizer widens and refines the grid when the screen count grows, so a naive sweep
of `min_screens` moves two variables at one time. Each case here is one
`Campaign`: one scenario, one geometry, one seed, one grid, one screen plan.

THE CASES.
  rapid    the shipped RAPID preset: its own grid and its 5-screen plan. It is
           what a preset="rapid" user gets, end to end.
  prod     the shipped default. It is the production baseline of the budgets.
  pin05    5 screens (the rapid floor) on the PINNED grid. Against rapid it
           isolates the rapid GRID; against pin09 the rapid COUNT.
  pin07    7 screens (the WP7 converged count) on the pinned grid.
  pin09    the default screen count on the PINNED grid. prod against pin09
           isolates the GRID effect.
  pin15    the same grid, 15 screens.
  pin25    the same grid, 25 screens.
  pin40    the same grid, 40 screens.
  gnd09x4  the same grid and the pin09 plan, with the BOTTOM screen split into
           four sub-screens of equal integrated Cn2. It isolates the near-ground
           resolution alone, at an almost unchanged screen placement above.

THE MEASURED QUANTITIES. Each trial gives one composite loss

    loss_db = -10 log10(collected_power * smf_eta)

which is the fidelity-2 downlink SMF quantity of olb.models.waveoptics (the
vacuum references cancel in a comparison between cases, because the grid is
pinned). The study reports p50, p10, p5 and p1 of that loss, where pX is the
loss that the link EXCEEDS X percent of the time. So p5 is the loss exceeded
5 percent of the time, and p1 is the 1 percent worst loss. It also reports the
FADE DEPTH pX - p50, so two cases with a different vacuum coupling still
compare.

NO CORRECTION. The fidelity-2 layer applies no tip-tilt removal and no adaptive
optics (backlog 2-AO), and it has no correction stage that a switch could turn
on. So the SMF coupled power here is the RAW uncorrected atmosphere, and the
verdict of this study applies to that case only.

THE SECOND PRIMARY QUANTITY IS THE POINT IRRADIANCE. Each stored trial gives
the irradiance I of the CENTRE PIXEL of the receive field. The study reports
the POINT fade

    point_loss_db = -10 log10(I / <I>)

with the same p50, p10, p5, p1 and the same bootstrap. A point receiver is the
sharpest probe of the near-ground resolution there is: nothing averages it, and
a fibre pays the point figure, not the aperture figure (the receiver-kind
question of backlog 2-I3).

It also reports two scintillation indices, both the normalised irradiance
variance sigma^2 = <x^2>/<x>^2 - 1 = var(x)/mean(x)^2. Source: Andrews and
Phillips, Laser Beam Propagation through Random Media, 2nd ed. (2005),
DOI 10.1117/3.626196, Ch. 8 (the scintillation index is the normalised
irradiance variance). olb holds that definition in
olb.turbulence.andrews.scintillation.
  sigma2_I_point  on the centre-pixel irradiance (the POINT index).
  sigma2_P        on the collected power (the APERTURE index). It is a FOOTNOTE
                  only: a 700 mm bucket averages the irradiance strongly, so
                  this index is of the order 1e-3 and it is INSENSITIVE to the
                  screen plan. A flat sigma2_P proves nothing about the tail.

THE 2-N6 MEASUREMENTS. The study also records, for each case, the wall time of
the run, the seconds per trial, the bytes on disk and the stored trial count,
and it repeats the quantiles at a growing trial count (100 to 1000) to show how
a tail estimate settles. Those are the large-campaign numbers of backlog 2-N6.

VALIDATION ONLY. The script reads the production layer. It changes no olb
module. It uses the private planner helpers of
olb.waveoptics.turbulence.sampling to build the near-ground refinement; that is
acceptable in a validation script, and the split is ASSERTED to conserve the
parent plan (the integrated Cn2, the Rytov variance and r0_total).

Sources:
- Andrews and Phillips, Laser Beam Propagation through Random Media, 2nd ed.
  (2005), DOI 10.1117/3.626196. Ch. 8, Eq. (20), printed p. 264: the plane-wave
  Rytov variance and its path weight. Ch. 12, Eq. (14), printed p. 482: the
  slant airmass sec(zeta). Ch. 12, Eq. (23): the composite Fried parameter.
- Fried, DOI 10.1364/JOSA.56.001372. The Fried parameter r0.
- Schmidt, Numerical Simulation of Optical Wave Propagation with Examples in
  MATLAB (2010), DOI 10.1117/3.866274, Ch. 9. The split-step method and the
  layer moment rule Eq. (9.65), printed p. 164.

Run it from the repository root:

    python -m validation.tail_convergence.tail_convergence
    python -m validation.tail_convergence.tail_convergence --elevation 20
    python -m validation.tail_convergence.tail_convergence --analyse-only
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
from olb.turbulence.andrews.beam import wavenumber
from olb.turbulence.andrews.paths import sec_zeta
from olb.turbulence.profiles import get_c2n
from olb.waveoptics.turbulence.campaign import Campaign
from olb.waveoptics.turbulence.sampling import (DEFAULT_H_TOP_M, PRESETS,
                                                ScreenPlan, _composite_r0,
                                                _cumtrapz, _integration_heights,
                                                _screen_rytov, turbulent_grid)
from olb.waveoptics.turbulence.screens import screen_r0

HERE = os.path.dirname(os.path.abspath(__file__))
CAMPAIGN_ROOT = os.environ.get("OLB_TAIL_ROOT", os.path.join(HERE, "campaigns"))

LAM = 1550e-9
SEED = 20260904

# The pinned grid comes from the sizer of THIS screen count. 15 screens gives
# the finer grid (n = 1024 at 30 deg), so every case of the sweep runs on a grid
# that is at least as good as its own sizer would ask for.
PIN_N = 15

# The screen counts of the pinned sweep.
PIN_COUNTS = (9, 15, 25, 40)

# The number of sub-screens of the near-ground refinement case.
GROUND_SPLIT = 4

# The exceedance probabilities. pX is the loss that the link EXCEEDS X percent
# of the time, so it is the (100 - X) percentile of the loss.
EXCEEDANCE = (0.50, 0.10, 0.05, 0.01)

# The bootstrap of every quantile: the resample count and the interval.
N_BOOT = 1000
BOOT_INTERVAL = 0.68
BOOT_SEED = 7

# The trial counts of the growth check (backlog 2-N6).
GROWTH_N = (100, 200, 400, 600, 800, 1000)

ALL_CASES = ("rapid", "prod", "pin05", "pin07", "pin09", "pin15", "pin25",
             "pin40", "gnd09x4")

# One line for each case. It says what the case ISOLATES.
CASE_DOC = {
    "rapid": "the shipped RAPID preset: its own grid AND its 5-screen plan. It "
             "is what a preset='rapid' user gets.",
    "prod": "the shipped default: the standard sizer picks the grid AND the "
            "plan. It is the production baseline.",
    "pin05": "5 screens (the rapid floor) on the pinned grid. Against rapid it "
             "isolates the rapid GRID; against pin09 the rapid COUNT.",
    "pin07": "7 screens (the WP7 converged count) on the pinned grid.",
    "pin09": "9 screens on the pinned grid. Against prod it isolates the GRID.",
    "pin15": "15 screens on the pinned grid. It isolates the SCREEN COUNT.",
    "pin25": "25 screens on the pinned grid. It isolates the SCREEN COUNT.",
    "pin40": "40 screens on the pinned grid. It isolates the SCREEN COUNT.",
    "gnd09x4": f"the pin09 plan with the BOTTOM screen split into "
               f"{GROUND_SPLIT} equal-Cn2 sub-screens. It isolates the "
               f"NEAR-GROUND resolution alone.",
}


# ---------------------------------------------------------------------------
# The scenario
# ---------------------------------------------------------------------------

def scenario_and_geometry(elevation_deg):
    """Build the downlink scenario and the orbit of one elevation.

    The case is the presentation downlink: a 1550 nm space-to-ground link to a
    700 mm ground telescope with a single-mode-fibre receiver, from a 100 mm
    space terminal at 500 km.

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


def _site_cn2(scn):
    """Give the site Hufnagel-Valley zenith Cn2 callable cn2(h).

    It repeats the default of olb.waveoptics.turbulence.sampling._plan_space,
    so the study integrates the SAME atmosphere that the planner integrates.
    Source: Andrews and Phillips, DOI 10.1117/3.626196, Ch. 12, Eq. (1),
    printed p. 481, through olb.turbulence.profiles.get_c2n.
    """
    site = scn.channel.site

    def cn2(h):
        return get_c2n(h, site.wind_rms_m_s, site.cn2_ground)

    return cn2


def _sizer(scn, geom, min_screens):
    """Size the grid and plan the screens with a changed screen floor.

    Args:
        scn:         the SpaceScenario.
        geom:        the CircularOrbit.
        min_screens: the screen-count floor of the preset copy.

    Returns:
        The triple (GridSpec, ScreenPlan, the warning texts).
    """
    preset = replace(PRESETS["standard"], name=f"std{int(min_screens)}",
                     min_screens=int(min_screens))
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        grid, plan, _ = turbulent_grid(scn, geom, preset=preset)
    return grid, plan, sorted({str(w.message) for w in caught})


# ---------------------------------------------------------------------------
# The near-ground refinement
# ---------------------------------------------------------------------------

def split_bottom_screen(scn, geom, plan, n_sub):
    """Split the BOTTOM screen of a space plan into n_sub sub-screens.

    THE BOTTOM SCREEN is the last element of `plan.z_m`, the largest distance
    from the input plane, so it is the LOWEST slab of the atmosphere. It holds
    the ground boundary layer, and it carries most of the integrated Cn2.

    THE METHOD.
    1. Integrate the slant Cn2 density m(h) = Cn2(h) sec(zeta) over the same
       internal height grid that the planner uses. The airmass sec(zeta) is
       Andrews and Phillips, DOI 10.1117/3.626196, Ch. 12, Eq. (14), printed
       p. 482.
    2. Invert the cumulative integral to find the top height h1 of the bottom
       slab: INT_0^h1 m dh equals the bottom screen `cn2_int_m13`.
    3. Cut [0, h1] at equal cumulative Cn2, so each sub-slab holds the same
       integrated Cn2.
    4. Put each sub-screen at the Cn2-weighted centroid of its sub-slab, the
       placement rule of _plan_space_continuous.
    5. Take the Rytov contribution of each sub-slab from the same Rytov density
       that the planner integrates, w(h) = 2.25 k^(7/6) Cn2(h) sec (h sec)^(5/6)
       (olb.waveoptics.turbulence.sampling._screen_rytov; Andrews and Phillips,
       Ch. 8, Eq. (20), printed p. 264). The distance to the ground receiver is
       h sec(zeta).
    6. Take r0 of each sub-screen from screen_r0 (Fried,
       DOI 10.1364/JOSA.56.001372), and r0_total from _composite_r0.

    The function ASSERTS that the split conserves the parent plan.

    Args:
        scn:   the SpaceScenario.
        geom:  the CircularOrbit.
        plan:  the parent ScreenPlan.
        n_sub: the number of sub-screens.

    Returns:
        A new ScreenPlan.
    """
    lam = scn.tx_terminal.wavelength_m
    k = wavenumber(lam)
    h_top = DEFAULT_H_TOP_M
    elevation = float(np.min(np.asarray(geom.elevation_deg, dtype=float)))
    sec = float(sec_zeta(elevation))

    h = _integration_heights(h_top)
    cn2_h = np.asarray(_site_cn2(scn)(h), dtype=float)
    z_of_h = (h_top - h) * sec
    m_dens = cn2_h * sec
    w_dens = _screen_rytov(k, cn2_h * sec, np.maximum(h * sec, 0.0))
    m_cum = _cumtrapz(m_dens, h)
    w_cum = _cumtrapz(w_dens, h)
    zm_cum = _cumtrapz(z_of_h * m_dens, h)

    # Step 2: the top height of the bottom slab.
    m_bottom = float(plan.cn2_int_m13[-1])
    h1 = float(np.interp(m_bottom, m_cum, h))

    # Step 3: the equal-Cn2 cut of [0, h1].
    targets = m_bottom * np.arange(1, n_sub) / n_sub
    h_edges = np.concatenate(([0.0], np.interp(targets, m_cum, h), [h1]))

    def _seg(cum, lo, hi):
        return float(np.interp(hi, h, cum) - np.interp(lo, h, cum))

    # Steps 4 and 5. The sub-slabs run from the ground up, so the reverse order
    # gives an ASCENDING z, the order that the split step wants.
    sub_z, sub_m, sub_s2 = [], [], []
    for j in range(n_sub - 1, -1, -1):
        lo, hi = h_edges[j], h_edges[j + 1]
        m_slab = _seg(m_cum, lo, hi)
        sub_m.append(m_slab)
        sub_z.append(_seg(zm_cum, lo, hi) / m_slab)
        sub_s2.append(_seg(w_cum, lo, hi))
    sub_z = np.asarray(sub_z)
    sub_m = np.asarray(sub_m)
    sub_s2 = np.asarray(sub_s2)

    z_new = np.concatenate((plan.z_m[:-1], sub_z))
    m_new = np.concatenate((plan.cn2_int_m13[:-1], sub_m))
    s2_new = np.concatenate((plan.sigma2_r[:-1], sub_s2))
    r0_new = screen_r0(m_new, lam)
    r0_total_new = _composite_r0(r0_new)

    # The conservation checks. The split must move the SCREEN COUNT only.
    assert abs(m_new.sum() / plan.cn2_int_m13.sum() - 1.0) < 1e-3, \
        (m_new.sum(), plan.cn2_int_m13.sum())
    assert abs(s2_new.sum() / plan.sigma2_r.sum() - 1.0) < 1e-3, \
        (s2_new.sum(), plan.sigma2_r.sum())
    assert abs(r0_total_new / plan.r0_total_m - 1.0) < 1e-6, \
        (r0_total_new, plan.r0_total_m)
    assert np.all(np.diff(z_new) > 0.0), z_new

    return ScreenPlan(z_m=z_new, cn2_int_m13=m_new, r0_m=r0_new,
                      sigma2_r=s2_new, z_total_m=plan.z_total_m,
                      r0_total_m=r0_total_new, direction=plan.direction)


# ---------------------------------------------------------------------------
# The cases
# ---------------------------------------------------------------------------

def build_cases(scn, geom, names):
    """Build the (grid, plan) pair of each named case.

    Args:
        scn:   the SpaceScenario.
        geom:  the CircularOrbit.
        names: the case names.

    Returns:
        A dict name -> {"grid", "plan", "warnings", "preset"}. A None grid and
        a None plan mean the case lets the Campaign size itself from the preset
        (the `prod` and the `rapid` cases).
    """
    pinned_grid, _, pin_warns = _sizer(scn, geom, PIN_N)
    out = {}
    for name in names:
        if name == "prod":
            out[name] = {"grid": None, "plan": None, "warnings": []}
        elif name == "rapid":
            out[name] = {"grid": None, "plan": None, "warnings": [],
                         "preset": "rapid"}
        elif name.startswith("pin"):
            _, plan, warns = _sizer(scn, geom, int(name[3:]))
            out[name] = {"grid": pinned_grid, "plan": plan,
                         "warnings": warns + pin_warns}
        elif name == "gnd09x4":
            _, plan9, warns = _sizer(scn, geom, 9)
            out[name] = {"grid": pinned_grid,
                         "plan": split_bottom_screen(scn, geom, plan9,
                                                     GROUND_SPLIT),
                         "warnings": warns + pin_warns}
        else:
            raise ValueError(f"unknown case {name!r}. Use one of {ALL_CASES}.")
    return out


def _bottom_share(plan):
    """Give the Cn2 share of the LOWEST screen of a plan."""
    return float(plan.cn2_int_m13[-1] / plan.cn2_int_m13.sum())


def _bottom_height_m(plan, geom):
    """Give the height above the ground of the LOWEST screen, in m."""
    elevation = float(np.min(np.asarray(geom.elevation_deg, dtype=float)))
    sec = float(sec_zeta(elevation))
    return float((plan.z_total_m - plan.z_m[-1]) / sec)


# ---------------------------------------------------------------------------
# The statistics
# ---------------------------------------------------------------------------

def _quantiles(loss_db, probs=EXCEEDANCE):
    """Give the loss that the link EXCEEDS each probability.

    A fade is a LARGE loss, so the loss exceeded a fraction q of the time is the
    (1 - q) quantile of the loss sample.

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
    of each resample, and it gives the half-width of the central
    BOOT_INTERVAL band. The rng seed is FIXED, so a rerun of the analysis gives
    the same interval.

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
    olb.turbulence.andrews.scintillation.
    """
    a = np.asarray(x, dtype=float)
    return float(a.var() / a.mean() ** 2)


def _centre_irradiance(result):
    """Give the irradiance of the CENTRE PIXEL of each stored trial.

    The centre pixel is the flat index (n // 2) * n + (n // 2), the axis pixel
    of the CircAperture convention that olb.waveoptics.turbulence.run
    _field_patch uses.

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


# ---------------------------------------------------------------------------
# The analysis of one case
# ---------------------------------------------------------------------------

def analyse(camp, name, n_trials, growth=False):
    """Measure one campaign.

    Args:
        camp:     the Campaign.
        name:     the case name.
        n_trials: the number of trials to read.
        growth:   True adds the growth-with-n table (backlog 2-N6).

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
    # The composite fidelity-2 downlink SMF loss. The budget Term takes the
    # same product and it divides by the vacuum references, which are one
    # constant for one grid. See olb.models.waveoptics.
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
        "case": name,
        "doc": CASE_DOC[name],
        "preset": camp.preset,
        "n_trials": int(loss.size),
        "n_screens": int(camp.plan.z_m.size),
        "grid_n": int(camp.grid.n),
        "grid_size_m": float(camp.grid.size_m),
        "grid_pixel_m": float(camp.grid.pixel_m),
        "bottom_screen_h_m": None,          # the caller fills it (it needs the geometry)
        "bottom_cn2_share": _bottom_share(camp.plan),
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

def _plot(rows, growth_rows, survival, fig_dir):
    """Draw the three figures. Give the list of the written paths."""
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

    # ---- 1. the quantiles against the screen count ----
    fig, ax = plt.subplots(figsize=(7.5, 5.0))
    pinned = [r for r in rows if r["case"].startswith("pin")]
    pinned.sort(key=lambda r: r["n_screens"])
    for key in names:
        if pinned:
            ax.errorbar([r["n_screens"] for r in pinned],
                        [r["quantiles_db"][key] for r in pinned],
                        yerr=[r["quantiles_half_db"][key] for r in pinned],
                        marker="o", capsize=3, label=f"{key} (pinned grid)")
    for r in rows:
        if r["case"] == "prod":
            ax.plot([r["n_screens"]] * len(names),
                    [r["quantiles_db"][k] for k in names], "kx",
                    markersize=9, label="prod (own grid)")
        if r["case"] == "gnd09x4":
            ax.plot([r["n_screens"]] * len(names),
                    [r["quantiles_db"][k] for k in names], "r+",
                    markersize=11, label="gnd09x4 (ground split)")
    ax.set_xlabel("screens")
    ax.set_ylabel("loss [dB]")
    ax.set_title("The SMF fade quantiles against the screen count")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    path = os.path.join(fig_dir, "tail_vs_screens.png")
    fig.savefig(path, dpi=140)
    plt.close(fig)
    written.append(path)

    # ---- 2. the empirical survival functions: the SMF loss, the point fade ----
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

    # ---- 3. the growth with the trial count ----
    if growth_rows:
        fig, ax = plt.subplots(figsize=(7.5, 5.0))
        for key in names:
            ax.errorbar([g["n"] for g in growth_rows],
                        [g["quantiles_db"][key] for g in growth_rows],
                        yerr=[g["half_db"][key] for g in growth_rows],
                        marker="o", capsize=3, label=key)
        ax.set_xlabel("trials")
        ax.set_ylabel("loss [dB]")
        ax.set_title("How the tail estimate settles with the trial count")
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
                    help="the trials of each case")
    ap.add_argument("--workers", type=int, default=16,
                    help="the processes of the campaign pool")
    ap.add_argument("--block-size", type=int, default=50,
                    help="the trials in one block file")
    ap.add_argument("--cases", nargs="+", default=list(ALL_CASES),
                    help=f"the cases to run, from {list(ALL_CASES)}")
    ap.add_argument("--analyse-only", action="store_true",
                    help="skip the runs and read what is stored")
    args = ap.parse_args()

    el = float(args.elevation)
    tag = f"el{el:02.0f}"
    log_path = os.path.join(HERE, f"tail_convergence_{tag}.log")
    # The log is written LINE BY LINE, so a killed run keeps its log.
    with open(log_path, "w", encoding="utf-8"):
        pass

    def say(text=""):
        print(text, flush=True)
        with open(log_path, "a", encoding="utf-8") as fh:
            fh.write(text + "\n")

    scn, geom = scenario_and_geometry(el)
    cases = build_cases(scn, geom, args.cases)

    say("The fidelity-2 SMF fade-tail convergence study (backlog 2-I2T, 2-N6)")
    say(f"elevation     : {el:.1f} deg      trials/case: {args.n_trials}")
    say(f"seed          : {SEED}      preset: standard (rapid for the rapid "
        f"case)      block size: {args.block_size}      "
        f"workers: {args.workers}")
    say(f"scenario      : downlink, 1550 nm, 500 km, 700 mm ground SMF, "
        f"100 mm space terminal")
    say(f"mode          : {'ANALYSE ONLY' if args.analyse_only else 'RUN'}")
    say("correction    : NONE. Fidelity 2 applies no tip-tilt removal and no "
        "AO (backlog 2-AO); the SMF power is the raw atmosphere.")
    say("primary       : the SMF coupled-power fade and the POINT irradiance "
        "fade. The aperture index s2_P is a footnote (strongly averaged).")
    if args.n_trials < 1000:
        say("WARNING: under 1000 trials, so the p1 tail is UNDER-SAMPLED. Read "
            "the p10 and p5 rows only.")
    say()

    say("THE CASES")
    for name in args.cases:
        say(f"  {name:<9s} {CASE_DOC[name]}")
    say()

    # ---- the grid and the plan of each case ----
    camps = {}
    for name in args.cases:
        root = os.path.join(CAMPAIGN_ROOT, tag, name)
        spec = cases[name]
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            # The stored campaigns are DOUBLE precision (2026-09-04).
            camp = Campaign(scn, geom, root, seed=SEED,
                            precision="double",
                            preset=spec.get("preset", "standard"),
                            block_size=args.block_size,
                            grid=spec["grid"], plan=spec["plan"])
        camps[name] = camp
        for w in sorted({str(x.message) for x in caught}):
            say(f"  {name}: campaign warning: {w}")
        for w in spec["warnings"]:
            say(f"  {name}: sizer warning: {w}")

        plan = camp.plan
        elev_sec = float(sec_zeta(el))
        heights = (plan.z_total_m - plan.z_m) / elev_sec
        share = plan.cn2_int_m13 / plan.cn2_int_m13.sum()
        say(f"  {name}: grid {camp.grid.n} px, {camp.grid.size_m:.3f} m "
            f"({camp.grid.pixel_m * 1e3:.2f} mm px); "
            f"{plan.z_m.size} screens; r0_total {plan.r0_total_m * 100:.2f} cm")
        say(f"    screen h [m] : "
            + " ".join(f"{v:.0f}" for v in heights))
        say(f"    Cn2 share    : "
            + " ".join(f"{v:.3f}" for v in share))
        say(f"    sigma2_r     : "
            + " ".join(f"{v:.4f}" for v in plan.sigma2_r))
        say()

    # ---- the runs ----
    timing = {}
    if not args.analyse_only:
        for name in args.cases:
            camp = camps[name]
            say(f"RUN {name}: {args.n_trials} trials")
            t0 = time.perf_counter()
            n_done = camp.run(args.n_trials, workers=args.workers,
                              progress=True)
            wall = time.perf_counter() - t0
            timing[name] = {"wall_s": float(wall), "n_stored": int(n_done)}
            per = wall / max(n_done, 1)
            say(f"  {name}: {n_done} trials on disk, {wall:.1f} s wall, "
                f"{per:.2f} s/trial (this call), "
                f"{_dir_bytes(camp.root_dir) / 1e6:.1f} MB on disk")
            # A CHECKPOINT: one partial result as soon as the case is done, so
            # a long run is not a black box. The scalars load in a moment.
            quick = camp.load(n_done, fields=False)
            ql = -10.0 * np.log10(np.array(
                [t.collected_power * t.smf_eta for t in quick.trials]))
            qq = _quantiles(ql)
            say(f"  CHECKPOINT {name}: SMF loss mean {ql.mean():.2f} dB, "
                f"p50 {qq[0]:.2f}, p10 {qq[1]:.2f}, p5 {qq[2]:.2f}, "
                f"p1 {qq[3]:.2f} dB (n = {ql.size})")
            say()

    # ---- the analysis ----
    rows, growth_rows = [], []
    survival = {"the SMF loss": {}, "the point fade": {}}
    # The growth check runs on pin09 when it is present, else on the first case.
    growth_case = "pin09" if "pin09" in args.cases else args.cases[0]
    for name in args.cases:
        camp = camps[name]
        if camp.n_stored == 0:
            say(f"  {name}: no stored trial. Skipped.")
            continue
        n = min(args.n_trials, camp.n_stored)
        want_growth = name == growth_case
        row, loss, point_loss = analyse(camp, name, n, growth=want_growth)
        row["bottom_screen_h_m"] = _bottom_height_m(camp.plan, geom)
        row["run_wall_s"] = timing.get(name, {}).get("wall_s")
        if want_growth and row.get("growth"):
            growth_rows = row["growth"]
            row["growth_case"] = True
        rows.append(row)
        survival["the SMF loss"][name] = loss
        survival["the point fade"][name] = point_loss

    # ---- the convergence table ----
    say("THE CONVERGENCE TABLE (pX = the loss EXCEEDED X percent of the time)")
    hdr = (f"{'case':<9s}{'scr':>4s}{'grid':>6s}{'h_bot':>8s}{'share':>7s}"
           f"{'mean':>8s}{'p50':>8s}{'p10':>15s}{'p5':>15s}{'p1':>15s}"
           f"{'d_p5':>7s}{'d_p1':>7s}"
           f"{'pt_p5':>14s}{'pt_p1':>14s}{'s2_I':>8s}{'s2_P':>8s}"
           f"{'s/tr':>7s}{'MB':>7s}")
    say(hdr)
    say("-" * len(hdr))
    order = [c for c in ALL_CASES if c in {r["case"] for r in rows}]
    for name in order:
        r = next(x for x in rows if x["case"] == name)
        q, h = r["quantiles_db"], r["quantiles_half_db"]
        pq, ph = r["point_quantiles_db"], r["point_quantiles_half_db"]
        say(f"{r['case']:<9s}{r['n_screens']:4d}{r['grid_n']:6d}"
            f"{r['bottom_screen_h_m']:8.0f}{r['bottom_cn2_share']:7.3f}"
            f"{r['mean_db']:8.2f}{q['p50']:8.2f}"
            f"{q['p10']:9.2f}+-{h['p10']:4.2f}"
            f"{q['p5']:9.2f}+-{h['p5']:4.2f}"
            f"{q['p1']:9.2f}+-{h['p1']:4.2f}"
            f"{r['fade_depth_db']['p5']:7.2f}{r['fade_depth_db']['p1']:7.2f}"
            f"{pq['p5']:8.2f}+-{ph['p5']:4.2f}"
            f"{pq['p1']:8.2f}+-{ph['p1']:4.2f}"
            f"{r['sigma2_I_point']:8.3f}{r['sigma2_P']:8.3f}"
            f"{r['mean_wall_time_s']:7.2f}{r['disk_bytes'] / 1e6:7.1f}")
    say("  scr = screens, h_bot = the height of the lowest screen [m], "
        "share = its Cn2 share.")
    say("  p10, p5, p1 = the SMF loss EXCEEDED 10, 5, 1 percent of the time. "
        "d_p5 and d_p1 are the FADE DEPTH pX - p50 [dB].")
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
        say("HOW THE TAIL ESTIMATE SETTLES WITH THE TRIAL COUNT (backlog 2-N6)")
        say(f"{'n':>6s}{'p50':>9s}{'p10':>15s}{'p5':>15s}{'p1':>15s}")
        for g in growth_rows:
            q, h = g["quantiles_db"], g["half_db"]
            say(f"{g['n']:6d}{q['p50']:9.2f}"
                f"{q['p10']:9.2f}+-{h['p10']:4.2f}"
                f"{q['p5']:9.2f}+-{h['p5']:4.2f}"
                f"{q['p1']:9.2f}+-{h['p1']:4.2f}")
        say()

    # ---- the verdict ----
    say("THE VERDICT")
    pinned = sorted([r for r in rows if r["case"].startswith("pin")],
                    key=lambda r: r["n_screens"])
    verdict = {}

    def _delta(a, b, qkey, hkey, key):
        """Give (b - a, the combined bootstrap bar, the sigma count)."""
        d = b[qkey][key] - a[qkey][key]
        bar = float(np.hypot(a[hkey][key], b[hkey][key]))
        return float(d), bar, (abs(d) / bar if bar > 0 else float("inf"))

    # The two primary quantities, in the same words.
    faces = (("SMF", "quantiles_db", "quantiles_half_db"),
             ("POINT", "point_quantiles_db", "point_quantiles_half_db"))
    if len(pinned) >= 2:
        first, last = pinned[0], pinned[-1]
        for label, qkey, hkey in faces:
            for key in ("p5", "p1"):
                d, bar, sigmas = _delta(first, last, qkey, hkey, key)
                moves = sigmas > 1.0
                verdict[f"{label}_{key}"] = {
                    "delta_db": d, "bar_db": bar, "sigmas": float(sigmas),
                    "resolved": bool(moves),
                    "from_screens": first["n_screens"],
                    "to_screens": last["n_screens"]}
                say(f"  {label} {key}: {first['n_screens']} -> "
                    f"{last['n_screens']} screens moves the loss by {d:+.2f} "
                    f"dB, against a combined bootstrap bar of {bar:.2f} dB "
                    f"({sigmas:.1f} sigma). "
                    f"{'IT MOVES' if moves else 'NO RESOLVED MOVE'}.")
            spread = {k: float(max(r[qkey][k] for r in pinned)
                               - min(r[qkey][k] for r in pinned))
                      for k in ("p50", "p10", "p5", "p1")}
            verdict[f"{label}_spread_db"] = spread
            say(f"  {label} full spread across the pinned counts [dB]: "
                + ", ".join(f"{k} {v:.2f}" for k, v in spread.items()))
    gnd = next((r for r in rows if r["case"] == "gnd09x4"), None)
    base = next((r for r in rows if r["case"] == "pin09"), None)
    if gnd and base:
        for label, qkey, hkey in faces:
            for key in ("p5", "p1"):
                d, bar, sigmas = _delta(base, gnd, qkey, hkey, key)
                verdict[f"{label}_ground_split_{key}"] = {
                    "delta_db": d, "bar_db": bar, "sigmas": float(sigmas)}
                say(f"  {label} {key}: the near-ground split moves it by "
                    f"{d:+.2f} dB against a bar of {bar:.2f} dB "
                    f"({sigmas:.1f} sigma).")
    prod = next((r for r in rows if r["case"] == "prod"), None)
    if prod and base:
        for label, qkey, hkey in faces:
            d, bar, sigmas = _delta(base, prod, qkey, hkey, "p5")
            verdict[f"{label}_grid_effect_p5"] = {
                "delta_db": d, "bar_db": bar, "sigmas": float(sigmas)}
            # A "point" is one PIXEL, and the two grids have a different pixel.
            # The 2026-09-04 re-bin check showed the two fields agree at a
            # MATCHED averaging area, so a point difference here is pixel
            # averaging, not physics. Read the point across the pinned series.
            note = ("" if label == "SMF" else
                    " PIXEL-LIMITED: the two grids average a point over a "
                    f"different pixel ({prod['grid_pixel_m'] * 1e3:.2f} against "
                    f"{base['grid_pixel_m'] * 1e3:.2f} mm), so this is NOT a "
                    "physics difference.")
            say(f"  {label} p5: the GRID (prod against pin09) moves it by "
                f"{d:+.2f} dB against a bar of {bar:.2f} dB "
                f"({sigmas:.1f} sigma).{note}")
    say()

    stamp = {
        "study": "tail_convergence",
        "backlog": ["2-I2T", "2-N6"],
        "elevation_deg": el, "n_trials": args.n_trials, "seed": SEED,
        "preset": "standard", "block_size": args.block_size,
        "workers": args.workers, "wavelength_m": LAM,
        "exceedance": list(EXCEEDANCE), "n_bootstrap": N_BOOT,
        "bootstrap_interval": BOOT_INTERVAL,
        "analyse_only": bool(args.analyse_only),
        "timing": timing, "verdict": verdict, "cases": rows,
    }
    json_path = os.path.join(HERE, f"tail_convergence_{tag}_results.json")
    with open(json_path, "w", encoding="utf-8") as fh:
        json.dump(stamp, fh, indent=2)
    figs = _plot(rows, growth_rows, survival, os.path.join(HERE, "figures"))
    print(f"\nwrote {json_path}")
    for f in figs:
        print(f"wrote {f}")
    print(f"wrote {log_path}")


if __name__ == '__main__':
    main()
