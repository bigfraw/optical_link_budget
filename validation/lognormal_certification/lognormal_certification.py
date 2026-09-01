"""Certify the analytic aperture-averaged lognormal power draw (backlog 1-6).

THE QUESTION. In WEAK turbulence the cheap analytic route takes the
aperture-averaged index sigma2_P = A sigma2_I and it draws the received power
from a lognormal (olb.links.terrestrial.terrestrial_scintillation_term, and the
same closed form in olb.links.downlink._lognormal_term). But an aperture
integrates a CORRELATED lognormal field, and a sum of lognormals is not a
lognormal: as the diameter D grows, the power moves to a Gaussian and the tails
get thin. A finite Gaussian beam adds beam wander, which puts a POINTING tail on
the distribution that one index cannot describe.

THE TEST. Run the fidelity-2 split-step Monte Carlo of
olb.waveoptics.turbulence, histogram the aperture-collected power, and measure
the deep-fade quantiles against the analytic lognormal. Each case reports:

  - the analytic index sigma2_P against the measured index var(P)/mean(P)^2, so
    an INDEX error and a SHAPE error do not mix;
  - the fade quantiles (10 %, 5 %, and 1 % of the time) of the sim against the
    analytic lognormal AND against a lognormal REFIT to the measured index. The
    refit leg is the pure SHAPE test;
  - the skew of ln P. A lognormal power gives a Gaussian ln P, so its skew is
    0. A negative skew is the drift to a Gaussian POWER (thin fade tail); a
    strong negative skew is the pointing tail.

THE INDEX ERROR IS SPLIT IN TWO (added 2026-09-01). sigma2_P = A sigma2_I holds
TWO analytic models, and a sigma2_P comparison alone cannot say which one is
wrong: the POINT index sigma2_I (the Dios on-axis Gaussian beam-wave form, which
reads the launch waist only) and the aperture-averaging FILTER A (the Churnside
weak plane-wave fit). So the script measures the point index too:

  - ONE propagation for each trial gives the whole aperture sweep. The runner
    entry point olb.waveoptics.turbulence.run.propagate_turbulent_field returns
    the complex receive-plane field of one snapshot, and the script clips that
    ONE field at the point estimator AND at every receive diameter. So every
    aperture reads the SAME atmosphere, and the cost does not grow with the
    number of diameters.
  - THE POINT ESTIMATOR is the mean irradiance inside a small on-axis disc of
    diameter D_POINT_M = 8 mm, taken with the SAME clip that the apertures use.
    A single centre PIXEL is the literal point value, but the pixel is only 0.9
    to 1.8 mm on the grids here, so a one-pixel estimate rides on grid-scale
    noise. The 8 mm disc is 0.14 of the Fresnel scale sqrt(lambda L) = 5.6 cm,
    and the analytic filter says it averages the index down by about 2 % only
    (the script prints that number, a_point_analytic). The CAVEAT is the pixel: the
    grid holds no irradiance structure below one pixel, so the estimator is a
    lower bound on a true point index. The Fresnel scale is 30 to 60 pixels
    across here, so that bias is small.
  - The table then shows FOUR columns: sigma2_I analytic against sigma2_I
    measured (the POINT-INDEX error) and A analytic against
    A_eff = sigma2_P_sim / sigma2_I_sim (the FILTER error).

THE BEAM-FILLING FLAG (added 2026-09-01). The Churnside filter A is a
PLANE-WAVE fit: it assumes the aperture reads a piece of a much wider field. A
receive aperture that holds most of the beam breaks that assumption (backlog
2-N2), and its A comparison tests the aperture-holds-the-beam regime, not the
filter. So the script gives each case the vacuum received beam radius w(L) and
the captured power fraction

    eta_fill = 1 - exp(-2 (D/2)^2 / w(L)^2),

and it FLAGS a case with eta_fill past FILL_LIMIT = 0.5 as
BEAM-FILLING-LIMITED, in the log, in the results JSON and in the figure.

THE ABSOLUTE-IMPACT COLUMNS (added 2026-09-01). Where the averaging is heavy,
a LARGE relative index error moves a SMALL number of dB. So each case prints
the measured fade spread (the standard deviation of the loss in dB) next to the
5 % fade depth, and the index ratio sits beside them.

The sweep holds ONE horizontal path and it moves the receive diameter from
point-like to strong averaging (D/rho_0 from about 0.2 to 8), for a collimated
and for a diverged launch.

NOTE. Fidelity 2 models NO tip-tilt correction (backlog 2-AO), so the sim holds
the FULL beam wander. That is part of what this script tests: the analytic route
carries no pointing tail at all.

VALIDATION ONLY. The script reads the production layer and it changes no olb
module. It does NOT extend TurbTrial or TurbWaveResult (backlog 2-I1): it works
through the public single-snapshot field entry point instead.

Sources:
- Andrews and Phillips, Laser Beam Propagation through Random Media, 2nd ed.
  (2005), DOI 10.1117/3.626196. Ch. 4 (the Gaussian free-space beam radius and
  irradiance profile, Eqs. (7) and (8)), Ch. 5 (the lognormal irradiance PDF and
  the weak-fluctuation conditions), Ch. 6 (the coherence radius rho_0 and the
  Fried parameter r_0 = 2.1 rho_0), Ch. 10 (the weak aperture-averaging factor).
- Dios et al., Applied Optics 43 (2004) 3866, DOI 10.1364/AO.43.003866, Eq. 16.
  The on-axis Gaussian beam-wave scintillation index of the analytic leg.
- Schmidt, Numerical Simulation of Optical Wave Propagation, (2010),
  DOI 10.1117/3.866274, Ch. 9. The split-step method of the fidelity-2 leg.

Run it from the repository root:

    python -m validation.lognormal_certification.lognormal_certification
    python -m validation.lognormal_certification.lognormal_certification --full
"""

import argparse
import json
import os
import time
import warnings

import numpy as np
from scipy.stats import norm, skew

from olb.beam import gaussz, virtual_waist
from olb.geometry import HorizontalPath
from olb.links.terrestrial import terrestrial_scintillation_term
from olb.scenario import TerrestrialChannel, TerrestrialScenario
from olb.terminal import Aperture, Terminal, Transmitter
from olb.turbulence.andrews.scintillation import rytov_weak
from olb.turbulence.andrews.structure import coherence_radius, fried_parameter
from olb.turbulence.plane_wave_scintillation import aperture_averaging_factor_weak
from olb.waveoptics.field import Power
from olb.waveoptics.run import _clip
from olb.waveoptics.threader import Threader
from olb.waveoptics.turbulence.run import (_start_field,
                                           propagate_turbulent_field,
                                           propagate_turbulent_scenario)
from olb.waveoptics.turbulence.sampling import PRESETS, turbulent_grid

HERE = os.path.dirname(os.path.abspath(__file__))
LAM = 1550e-9
_LN10 = np.log(10.0)

# The fixed horizontal path. Cn2 = 3e-15 over 2 km gives sigma_R^2 near 0.2, so
# the path is FIRMLY weak (the rytov_weak "weak" tier, below 0.3). That is the
# band the analytic route claims, so it is the band to certify.
PATH_M = 2000.0
CN2 = 3e-15

# The launch. A 5 mm waist puts a beam radius near 0.2 m at 2 km, so the
# aperture sweep runs from far inside the beam to the beam size.
WAIST_M = 5e-3
# The diverged launch. The diffraction divergence of the 5 mm waist is about
# 1.0e-4 rad, so 2.0e-4 rad is a beam opened by about a factor of two.
DIVERGENCE_RAD = 2.0e-4

# The receive diameters, from point-like to strong averaging.
DIAMETERS_M = (0.01, 0.05, 0.15, 0.40)

# The point-index estimator: the mean irradiance in an on-axis disc of this
# diameter. See the module docstring (THE POINT ESTIMATOR).
D_POINT_M = 8e-3

# A case whose receive aperture holds more than this fraction of the beam power
# is BEAM-FILLING-LIMITED: the plane-wave aperture-averaging fit does not claim
# it (backlog 2-N2).
FILL_LIMIT = 0.5

# The fade exceedance probabilities to report. A fade at q is the loss that the
# link exceeds a fraction q of the time.
EXCEEDANCE = (0.10, 0.05, 0.01)

# A disagreement past this value at the 5 % fade is NOTABLE (the verdict rule).
NOTABLE_DB = 0.5

# The number of trials of the matched-seed cross-check against the old,
# one-run-for-each-aperture path.
N_MATCH_CHECK = 3


def _scenario(diameter_m, divergence_rad):
    """Build the terrestrial case: near = launch, far = the bucket receiver."""
    return TerrestrialScenario(
        near=Terminal(aperture_m=0.10, wavelength_m=LAM,
                      transmitter=Transmitter(waist_m=WAIST_M,
                                              divergence_rad=divergence_rad)),
        far=Terminal(aperture_m=diameter_m, wavelength_m=LAM,
                     detector=Aperture(sensitivity_dbm=-40)),
        channel=TerrestrialChannel(path_length_m=PATH_M,
                                   attenuation_db_per_km=0.0, cn2=CN2))


def _beam_radius_m(divergence_rad):
    """Give the VACUUM received beam radius w(L) of the launch, in m.

    A deliberately diverged transmitter is an ordinary Gaussian beam from a
    virtual waist w_v at the distance d behind the aperture (olb.beam), so the
    free-space radius at the range L is w(L) = gaussz(w_v, d + L). Source:
    Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196, Ch. 4, Eqs. (7)
    and (8), printed p. 87.
    """
    w_v, offset = virtual_waist(WAIST_M, divergence_rad, LAM)
    return float(gaussz(w_v, offset + PATH_M, LAM))


def _fill_fraction(diameter_m, w_L):
    """Give the vacuum power fraction that a diameter D catches, eta_fill.

    A Gaussian beam of 1/e^2 radius w carries the irradiance
    I(r) = (2/(pi w^2)) exp(-2 r^2 / w^2). The radial integral to r = D/2 is

        eta_fill = 1 - exp(-2 (D/2)^2 / w^2).

    Source of the irradiance profile: Andrews and Phillips, 2nd ed. (2005),
    DOI 10.1117/3.626196, Ch. 4, Eq. (8), printed p. 87. The same closed form
    is the on-axis limit of olb.models.coupling.terrestrial
    _mmf_encircled_efficiency.
    """
    return float(1.0 - np.exp(-2.0 * (diameter_m / 2.0) ** 2 / w_L ** 2))


def _lognormal_fades(sigma2_P):
    """Give the lognormal fade loss [dB] at each EXCEEDANCE probability.

    The analytic route holds sigma_l^2 = ln(1 + sigma2_P) and it draws the
    power P from a lognormal of unit mean. The loss that the link exceeds a
    fraction q of the time is

        fade(q) = -10/ln10 * ( -sigma_l^2/2 + sigma_l * Phi^-1(q) )

    with Phi^-1 the inverse standard normal CDF. This is the quantile face of
    olb.links.terrestrial.terrestrial_scintillation_term, written on the
    exceedance axis. Source of the lognormal irradiance PDF: Andrews and
    Phillips, 2nd ed. (2005), DOI 10.1117/3.626196, Ch. 5.
    """
    s2 = np.log(1.0 + sigma2_P)
    s = np.sqrt(s2)
    return {q: float(-10.0 / _LN10 * (-s2 / 2.0 + s * norm.ppf(q)))
            for q in EXCEEDANCE}


def _lognormal_mean_db(sigma2_P):
    """Give the mean loss [dB] of the unit-mean lognormal, (5/ln10) sigma_l^2."""
    return float((5.0 / _LN10) * np.log(1.0 + sigma2_P))


def _index(samples):
    """Give the normalised variance var(x)/mean(x)^2 of a sample set."""
    x = np.asarray(samples, dtype=float)
    return float(x.var() / x.mean() ** 2)


def _run_geometry(name, divergence_rad, n_trials, preset, seed, threader):
    """Measure ONE launch geometry: every aperture on the SAME trials.

    The runner propagates one snapshot for each trial and it clips that ONE
    field at the point estimator and at every receive diameter. So the whole
    aperture sweep costs one propagation set, and every aperture reads the same
    atmosphere.

    Returns:
        A tuple (cases, shared). cases is a list of one result dict for each
        diameter; shared holds the grid, the timing and the matched-seed check.
    """
    geom = HorizontalPath(PATH_M)
    p = PRESETS[preset] if isinstance(preset, str) else preset

    # THE GRID IS THE GRID OF THE LARGEST APERTURE. The sizer widens the grid
    # with the receive aperture, so the largest diameter gives the grid that
    # serves every diameter of the sweep.
    scn_ref = _scenario(max(DIAMETERS_M), divergence_rad)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        grid, plan, report = turbulent_grid(scn_ref, geom, preset=p)
    grid_warnings = sorted({str(w.message) for w in caught})

    # The reference power: the launched field after the transmit clip. It is
    # the SAME reference that TurbTrial.collected_power uses, so the matched
    # check below compares like with like.
    p_reference = Power(_start_field(scn_ref, grid, LAM, is_space=False))
    diameters = (D_POINT_M,) + DIAMETERS_M

    def one_trial(k):
        """Propagate snapshot k and read every aperture of the sweep."""
        F_rx, _, _ = propagate_turbulent_field(scn_ref, geom, seed=seed,
                                               trial=k, preset=p, grid=grid,
                                               plan=plan)
        return [float(Power(_clip(F_rx, D, 0.0)) / p_reference)
                for D in diameters]

    t0 = time.perf_counter()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        if threader is None or threader.max_workers == 1:
            rows = [one_trial(k) for k in range(n_trials)]
        else:
            rows = threader.map(one_trial, range(n_trials))
    sim_warnings = sorted({str(w.message) for w in caught})
    wall_s = time.perf_counter() - t0
    rows = np.asarray(rows, dtype=float)

    # THE MATCHED-SEED CHECK. The old path ran propagate_turbulent_scenario one
    # time for each aperture. The same seed and the same grid give the same
    # atmosphere, so the collected power of trial k of that runner must equal
    # the clip of the field of trial k here.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ref_run = propagate_turbulent_scenario(
            scn_ref, geom, n_trials=min(N_MATCH_CHECK, n_trials), seed=seed,
            preset=p, grid=grid, plan=plan)
    ref_power = np.array([t.collected_power for t in ref_run.trials])
    match = rows[:ref_power.size, -1]
    match_max_rel = float(np.max(np.abs(match / ref_power - 1.0)))

    point = rows[:, 0]
    sigma2_I_sim = _index(point)
    w_L = _beam_radius_m(divergence_rad)

    # The coherence radius rho_0 and the Fried parameter r_0 = 2.1 rho_0 of the
    # horizontal path. Source: Andrews and Phillips, 2nd ed. (2005),
    # DOI 10.1117/3.626196, Ch. 6, Eq. (64) and the text below it, printed
    # p. 194. A plane wave is the reference scale of the aperture-averaging
    # discussion, so the sweep is reported on the plane-wave rho_0.
    rho0 = float(coherence_radius(LAM, PATH_M, CN2, wave='plane'))
    r0 = float(fried_parameter(rho0))
    a_point = float(aperture_averaging_factor_weak(D_POINT_M, LAM, PATH_M))

    cases = []
    for j, D in enumerate(DIAMETERS_M, start=1):
        scn = _scenario(D, divergence_rad)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            term = terrestrial_scintillation_term(scn, geom)
        analytic_warnings = [str(w.message) for w in caught]

        sigma2_I = float(term.meta["sigma2_I"])
        sigma2_P = float(term.meta["sigma2_P"])
        sigma2_R = float(term.meta["sigma2_R"])
        Lambda = float(term.meta["Lambda"])
        A_analytic = float(term.meta["aperture_averaging_factor"])

        power = rows[:, j]
        # Normalise by the measured mean. The analytic fade faces describe a
        # unit-mean power, so the SHAPE comparison must divide out the
        # deterministic geometric loss that collected_power also holds.
        loss_db = -10.0 * np.log10(power / power.mean())
        sigma2_meas = _index(power)

        fade_analytic = _lognormal_fades(sigma2_P)
        fade_refit = _lognormal_fades(sigma2_meas)
        fade_emp = {q: float(np.quantile(loss_db, 1.0 - q)) for q in EXCEEDANCE}

        eta_fill = _fill_fraction(D, w_L)
        cases.append({
            "name": name,
            "diameter_m": D,
            "divergence_rad": divergence_rad,
            "n_trials": int(n_trials),
            "preset": p.name,
            "seed": int(seed),
            "grid_n": int(grid.n),
            "grid_size_m": float(grid.size_m),
            "n_screens": int(plan.z_m.size),
            "rho0_m": rho0,
            "r0_m": r0,
            "D_over_rho0": D / rho0,
            "D_over_r0": D / r0,
            "sigma2_R": sigma2_R,
            "Lambda": Lambda,
            "rytov_regime": rytov_weak(sigma2_R, Lambda),
            # --- the split index legs ---
            "beam_radius_m": w_L,
            "eta_fill": eta_fill,
            "beam_filling_limited": bool(eta_fill > FILL_LIMIT),
            "d_point_m": D_POINT_M,
            "a_point_analytic": a_point,
            "sigma2_I_analytic": sigma2_I,
            "sigma2_I_measured": sigma2_I_sim,
            "sigma2_I_ratio": sigma2_I_sim / sigma2_I,
            "aperture_averaging_factor": A_analytic,
            "a_effective": sigma2_meas / sigma2_I_sim,
            "a_ratio": (sigma2_meas / sigma2_I_sim) / A_analytic,
            # --- the whole-route index ---
            "sigma2_P_analytic": sigma2_P,
            "sigma2_P_measured": sigma2_meas,
            "sigma2_ratio": sigma2_meas / sigma2_P,
            # --- the absolute impact ---
            "loss_std_db": float(loss_db.std()),
            "mean_db_analytic": _lognormal_mean_db(sigma2_P),
            "mean_db_measured": float(loss_db.mean()),
            "skew_ln_power": float(skew(np.log(power))),
            "fade_analytic_db": {str(q): fade_analytic[q] for q in EXCEEDANCE},
            "fade_refit_db": {str(q): fade_refit[q] for q in EXCEEDANCE},
            "fade_measured_db": {str(q): fade_emp[q] for q in EXCEEDANCE},
            "delta_analytic_db": {str(q): fade_emp[q] - fade_analytic[q]
                                  for q in EXCEEDANCE},
            "delta_refit_db": {str(q): fade_emp[q] - fade_refit[q]
                               for q in EXCEEDANCE},
            "analytic_warnings": analytic_warnings,
            "sim_warnings": sim_warnings,
            "loss_db": loss_db.tolist(),
        })

    shared = {
        "name": name,
        "wall_s": wall_s,
        "grid_n": int(grid.n),
        "grid_size_m": float(grid.size_m),
        "grid_pixel_m": float(grid.pixel_m),
        "n_screens": int(plan.z_m.size),
        "beam_radius_m": w_L,
        "sigma2_I_measured": sigma2_I_sim,
        "a_point_analytic": a_point,
        "match_trials": int(ref_power.size),
        "match_max_rel_diff": match_max_rel,
        "grid_warnings": grid_warnings,
        "sim_warnings": sim_warnings,
        "report_step_over_limit_max": (None if report is None
                                       else float(report.step_over_limit_max)),
    }
    return cases, shared


def _plot(cases, path):
    """Draw one panel for each case: the measured PDF against the lognormal."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib is not installed; the PNG is skipped.")
        return None
    n = len(cases)
    ncol = 4
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.2 * ncol, 3.4 * nrow),
                             squeeze=False)
    for ax, c in zip(axes.ravel(), cases):
        loss = np.array(c["loss_db"])
        flagged = c["beam_filling_limited"]
        ax.hist(loss, bins=40, density=True,
                color="0.85" if flagged else "0.75", label="fidelity 2")
        # The analytic lognormal PDF on the loss axis. P is lognormal of unit
        # mean, so ln P is normal with mean -sigma_l^2/2 and sigma sigma_l, and
        # loss = -10 log10 P is normal with the signs turned over.
        s2 = np.log(1.0 + c["sigma2_P_analytic"])
        s = np.sqrt(s2)
        mu_db = 10.0 / _LN10 * (s2 / 2.0)
        sd_db = 10.0 / _LN10 * s
        x = np.linspace(loss.min(), loss.max(), 400)
        ax.plot(x, norm.pdf(x, mu_db, sd_db), "b-",
                label="lognormal (analytic sigma2_P)")
        s2r = np.log(1.0 + c["sigma2_P_measured"])
        ax.plot(x, norm.pdf(x, 10.0 / _LN10 * (s2r / 2.0),
                            10.0 / _LN10 * np.sqrt(s2r)), "r--",
                label="lognormal (refit)")
        for q, style in zip(EXCEEDANCE, (":", "-.", "-")):
            ax.axvline(c["fade_measured_db"][str(q)], color="k", ls=style, lw=1)
            ax.axvline(c["fade_analytic_db"][str(q)], color="b", ls=style, lw=1)
        title = (f"{c['name']}\nD={c['diameter_m'] * 100:.0f} cm, "
                 f"D/rho0={c['D_over_rho0']:.2f}, "
                 f"fill={c['eta_fill'] * 100:.0f} %")
        ax.set_title(title, fontsize=9,
                     color="0.45" if flagged else "black")
        if flagged:
            # A flagged case tests the aperture-holds-the-beam regime, not the
            # plane-wave averaging filter. Grey the panel so the reader sees it.
            ax.set_facecolor("0.94")
            ax.text(0.98, 0.95, "BEAM-FILLING-LIMITED", transform=ax.transAxes,
                    ha="right", va="top", fontsize=7, color="0.35")
        ax.set_xlabel("loss [dB]")
        ax.set_ylabel("pdf")
        ax.tick_params(labelsize=8)
    for ax in axes.ravel()[n:]:
        ax.axis("off")
    axes[0][0].legend(fontsize=7)
    fig.suptitle("Aperture-averaged power: fidelity 2 against the analytic "
                 "lognormal (black = measured fade, blue = analytic fade; "
                 "a grey panel is beam-filling-limited)", fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)
    return path


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--trials", type=int, default=None,
                    help="trials for each case (default 150 quick, 1500 full)")
    ap.add_argument("--preset", default=None,
                    help="the sampling preset (default rapid quick, standard "
                         "full)")
    ap.add_argument("--seed", type=int, default=20260901)
    ap.add_argument("--workers", type=int, default=None,
                    help="threads for the trials (default min(16, cores))")
    ap.add_argument("--full", action="store_true",
                    help="the long run: enough trials for the 1 %% fade tail")
    args = ap.parse_args()

    n_trials = args.trials or (1500 if args.full else 150)
    preset = args.preset or ("standard" if args.full else "rapid")
    threader = Threader(max_workers=args.workers)

    lines = []

    def say(text=""):
        print(text)
        lines.append(text)

    say("Aperture-averaged lognormal certification (backlog 1-6)")
    say(f"mode          : {'FULL' if args.full else 'QUICK'}")
    say(f"trials/case   : {n_trials}    preset: {preset}    "
        f"seed: {args.seed}    workers: {threader.max_workers}")
    say(f"path          : {PATH_M / 1e3:.1f} km horizontal, Cn2 = {CN2:.1e}, "
        f"lambda = {LAM * 1e9:.0f} nm")
    say(f"launch        : waist {WAIST_M * 1e3:.1f} mm, collimated and "
        f"diverged ({DIVERGENCE_RAD:.1e} rad)")
    say(f"point index   : mean irradiance in a {D_POINT_M * 1e3:.0f} mm on-axis "
        f"disc (the analytic filter averages it by "
        f"{(1 - aperture_averaging_factor_weak(D_POINT_M, LAM, PATH_M)) * 100:.2f} %)")
    if n_trials < 1000:
        say("WARNING: under 1000 trials, so the 1 % fade tail is UNDER-SAMPLED. "
            "Read the 10 % and 5 % rows only, and use --full for the deep tail.")
    say("NOTE: fidelity 2 models NO tip-tilt correction (backlog 2-AO), so the "
        "sim holds the FULL beam wander. The pointing tail is part of the test.")
    say()

    geometries = (("collimated", None), ("diverged", DIVERGENCE_RAD))
    cases = []
    shared = []
    for gname, div in geometries:
        t0 = time.perf_counter()
        gcases, gshared = _run_geometry(gname, div, n_trials, preset,
                                        args.seed, threader)
        cases.extend(gcases)
        shared.append(gshared)
        say(f"  ran {gname:<10s} {n_trials} trials, every aperture on the same "
            f"fields: grid {gshared['grid_n']}px/"
            f"{gshared['grid_size_m']:.2f} m ({gshared['grid_pixel_m'] * 1e3:.2f} "
            f"mm px), {gshared['n_screens']} screens, w(L) = "
            f"{gshared['beam_radius_m'] * 100:.1f} cm  "
            f"[{time.perf_counter() - t0:.1f} s]")
        say(f"      matched-seed check against propagate_turbulent_scenario "
            f"({gshared['match_trials']} trials, D = "
            f"{max(DIAMETERS_M) * 100:.0f} cm): max relative difference "
            f"{gshared['match_max_rel_diff']:.2e}")
        for w in gshared["sim_warnings"]:
            say(f"      sim warning: {w}")
    say()

    match_ok = all(s["match_max_rel_diff"] < 1e-9 for s in shared)
    say(f"REFACTOR CHECK: the one-propagation sweep reproduces the old "
        f"one-run-for-each-aperture path: {'PASS' if match_ok else 'FAIL'}")
    say()

    say("INDEX SPLIT: the POINT index and the AVERAGING FILTER, apart")
    hdr = (f"{'case':<11s}{'D [cm]':>8s}{'D/rho0':>8s}{'s2I an':>9s}"
           f"{'s2I sim':>9s}{'ratio':>7s}{'A an':>9s}{'A eff':>9s}"
           f"{'ratio':>7s}{'fill':>7s}  flag")
    say(hdr)
    say("-" * (len(hdr) + 8))
    for c in cases:
        say(f"{c['name']:<11s}{c['diameter_m'] * 100:8.1f}"
            f"{c['D_over_rho0']:8.2f}{c['sigma2_I_analytic']:9.4f}"
            f"{c['sigma2_I_measured']:9.4f}{c['sigma2_I_ratio']:7.2f}"
            f"{c['aperture_averaging_factor']:9.4f}{c['a_effective']:9.4f}"
            f"{c['a_ratio']:7.2f}{c['eta_fill']:7.2f}  "
            f"{'BEAM-FILLING-LIMITED' if c['beam_filling_limited'] else ''}")
    say("  s2I sim is one number for each launch: the point estimator does not "
        "read the receive aperture. A eff = sigma2_P_sim / sigma2_I_sim.")
    say("  A flagged case holds more than "
        f"{FILL_LIMIT * 100:.0f} % of the beam power, so its A comparison "
        "tests the aperture-holds-the-beam regime (backlog 2-N2), NOT the "
        "plane-wave filter.")
    say()

    say("WHOLE-ROUTE INDEX AND ABSOLUTE IMPACT")
    hdr3 = (f"{'case':<11s}{'D [cm]':>8s}{'s2P an':>10s}{'s2P sim':>10s}"
            f"{'ratio':>7s}{'sd loss':>9s}{'p5 sim':>9s}{'p5 an':>9s}"
            f"{'delta':>8s}{'skew lnP':>10s}  flag")
    say(hdr3)
    say("-" * (len(hdr3) + 8))
    for c in cases:
        say(f"{c['name']:<11s}{c['diameter_m'] * 100:8.1f}"
            f"{c['sigma2_P_analytic']:10.4f}{c['sigma2_P_measured']:10.4f}"
            f"{c['sigma2_ratio']:7.2f}{c['loss_std_db']:9.3f}"
            f"{c['fade_measured_db']['0.05']:9.3f}"
            f"{c['fade_analytic_db']['0.05']:9.3f}"
            f"{c['delta_analytic_db']['0.05']:8.3f}"
            f"{c['skew_ln_power']:10.2f}  "
            f"{'FILL' if c['beam_filling_limited'] else ''}")
    say("  sd loss, p5 and delta are dB. A large index ratio with a small "
        "sd loss is a big RELATIVE error that moves few dB.")
    say()

    say("FADE QUANTILES [dB], loss exceeded a fraction q of the time")
    hdr2 = (f"{'case':<11s}{'D [cm]':>8s}{'q':>7s}{'sim':>8s}{'analytic':>10s}"
            f"{'delta':>8s}{'refit':>8s}{'d refit':>9s}")
    say(hdr2)
    say("-" * len(hdr2))
    for c in cases:
        for q in EXCEEDANCE:
            k = str(q)
            say(f"{c['name']:<11s}{c['diameter_m'] * 100:8.1f}{q:7.2f}"
                f"{c['fade_measured_db'][k]:8.3f}"
                f"{c['fade_analytic_db'][k]:10.3f}"
                f"{c['delta_analytic_db'][k]:8.3f}"
                f"{c['fade_refit_db'][k]:8.3f}"
                f"{c['delta_refit_db'][k]:9.3f}")
    say()

    # THE VERDICT. The rule of the backlog item: a disagreement past 0.5 dB at
    # the 5 % fade is notable. The shape leg (the refit lognormal) says whether
    # the LOGNORMAL FAMILY is wrong; the analytic leg says whether the whole
    # cheap route (index and shape together) is trustworthy.
    worst_a = max(abs(c["delta_analytic_db"]["0.05"]) for c in cases)
    worst_s = max(abs(c["delta_refit_db"]["0.05"]) for c in cases)
    notable = [c for c in cases
               if abs(c["delta_analytic_db"]["0.05"]) > NOTABLE_DB]
    filled = [c for c in cases if c["beam_filling_limited"]]
    verdict = "PASS" if not notable else "FAIL"
    say(f"VERDICT: {verdict}")
    say(f"  worst |delta| at the 5 % fade, whole analytic route : "
        f"{worst_a:.3f} dB")
    say(f"  worst |delta| at the 5 % fade, SHAPE only (refit)   : "
        f"{worst_s:.3f} dB")
    if filled:
        say(f"  BEAM-FILLING-LIMITED cases (eta_fill past "
            f"{FILL_LIMIT:.2f}); read their A columns as a REGIME note, not as "
            f"a filter error:")
        for c in filled:
            say(f"    {c['name']:<11s} D={c['diameter_m'] * 100:5.1f} cm  "
                f"eta_fill={c['eta_fill']:.2f}  "
                f"w(L)={c['beam_radius_m'] * 100:.1f} cm  "
                f"A ratio={c['a_ratio']:.2f}")
    if notable:
        say(f"  notable cases (past {NOTABLE_DB} dB at the 5 % fade):")
        for c in notable:
            say(f"    {c['name']:<11s} D={c['diameter_m'] * 100:5.1f} cm  "
                f"delta={c['delta_analytic_db']['0.05']:+.2f} dB  "
                f"sigma2 ratio={c['sigma2_ratio']:.2f}  "
                f"skew lnP={c['skew_ln_power']:+.2f}")
        say("  A FAIL BOUNDS the cheap route. Where the index ratio is near 1 "
            "and the shape leg is small, the lognormal FAMILY is adequate and "
            "the index is the fault. Read the INDEX SPLIT table to see whether "
            "the POINT index or the averaging FILTER carries that fault. Where "
            "the shape leg is large, use a composite (lognormal x pointing) "
            "model or a direct empirical sampler.")
    else:
        say("  The cheap analytic weak-turbulence route gives a trustworthy "
            "power distribution over the tested band.")
    if n_trials < 1000:
        say("  QUICK MODE: this verdict reads the 10 % and 5 % fades only. Run "
            "--full before you record a certification.")

    stamp = {
        "mode": "full" if args.full else "quick",
        "n_trials": n_trials, "preset": preset, "seed": args.seed,
        "workers": threader.max_workers,
        "path_m": PATH_M, "cn2": CN2, "wavelength_m": LAM,
        "waist_m": WAIST_M, "divergence_rad": DIVERGENCE_RAD,
        "d_point_m": D_POINT_M, "fill_limit": FILL_LIMIT,
        "exceedance": list(EXCEEDANCE), "notable_db": NOTABLE_DB,
        "verdict": verdict,
        "match_check_pass": match_ok,
        "worst_delta_p5_db": worst_a, "worst_shape_delta_p5_db": worst_s,
        "geometries": shared,
        "cases": cases,
    }
    json_path = os.path.join(HERE, "lognormal_certification_results.json")
    with open(json_path, "w") as fh:
        json.dump(stamp, fh, indent=2)
    # The house rule of validation/: a figure goes into the figures/ subfolder.
    fig_dir = os.path.join(HERE, "figures")
    os.makedirs(fig_dir, exist_ok=True)
    png = _plot(cases, os.path.join(fig_dir, "lognormal_certification.png"))
    log_path = os.path.join(HERE, "lognormal_certification.log")
    with open(log_path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"\nwrote {json_path}")
    if png:
        print(f"wrote {png}")
    print(f"wrote {log_path}")


if __name__ == '__main__':
    main()
