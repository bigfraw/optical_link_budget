'''
The fidelity-2 uplink through reciprocity, against the coupled-flux Monte Carlo.

The satellite of an uplink sits outside the atmosphere, so a field simulation
cannot put a receive aperture at the far end of a gridded path. It does not
need to. The turbulent atmosphere is RECIPROCAL: the flux that the satellite
sees is the overlap of the received DOWNLINK field with the ground transmit
mode. See Shapiro, DOI 10.1364/JOSA.61.000492. The runner propagates the
downlink slab and it reads that overlap as `eta_turb`, against the free-space
baseline. So `-10*log10(eta_turb)` is the uplink turbulence loss, on the same
free-space reference as the analytic Terms.

The fidelity-1 model of the same number is the Dios coupled-flux Monte Carlo of
`olb.turbulence.uplink_flux`, which `olb.links.uplink.uplink_turbulence_term`
calls. It draws a beam-wander offset, it evaluates the off-axis Dios
scintillation at that offset, and it rescales the flux onto the free-space
baseline. Its loss is `-10*log10(Is_summed)`.

The two are compared here at the zenith and at 30 degrees.

NO POINTING JITTER. The coupled-flux model folds a mechanical jitter into the
same wander variance beta2 that carries the turbulence, so a jitter would
change one model and not the other. Both terminals therefore carry
`pointing_jitter_rad=0.0`, and the comparison is turbulence-only.

WHAT TO EXPECT.

- The MEANS agree to about 1 dB in the weak regime. The zenith case is weak.
  The 30-degree case is NOT: the coupled-flux model prints
  `weak_fluctuation_valid = False` there, so its own author does not trust it.
  The script prints that flag next to the numbers.
- The TAILS differ by construction, and the script REPORTS the difference
  instead of testing it. The coupled-flux tail is a PARAMETRIC lognormal: one
  wander offset per draw, and a lognormal irradiance about the local mean. The
  wave-optics tail is the tail of a real field, so it carries the speckle of a
  broken wavefront, the deep destructive interference of the whole aperture,
  and the correlation between the wander and the scintillation. A field Monte
  Carlo goes deeper in the tail than a parametric lognormal does.

The figure goes to `examples/waveoptics/turbulent_uplink_reciprocity.png`.

The layer builds NO Term and it changes NO budget. See the README.

Run from the repo root:
    python -m examples.waveoptics.turbulent_uplink_reciprocity
'''

import dataclasses
import time
import warnings

import matplotlib.pyplot as plt
import numpy as np

from olb import CircularOrbit, Terminal, Transmitter
from olb.scenario import Channel, SpaceScenario
from olb.turbulence.profiles import DEFAULT_HS, default_cn2_profile
from olb.turbulence.uplink_flux import _flux_result
from olb.waveoptics import Threader
from olb.waveoptics.turbulence import (propagate_turbulent_field,
                                       propagate_turbulent_scenario,
                                       turbulent_grid)

WAVELENGTH_M = 1550e-9
ALTITUDE_M = 600e3
WAIST_M = 0.10              # the transmit beam radius at the launch plane
GROUND_APERTURE_M = 0.50    # 5 waists: the launch clip takes 3e-5 of the power
SPACE_APERTURE_M = 0.30
ELEVATIONS_DEG = (90.0, 30.0)
PRESET = "rapid"
N_TRIALS = 200
BLOCK = 50                  # trials for each progress line
SEED = 20260826
FLUX_SAMPLES = 3000         # coupled-flux Monte Carlo draws

# The trials are independent, so they run across threads. None takes one worker
# for each core. See olb.waveoptics.Threader.
THREADER = Threader()

PNG = "examples/waveoptics/turbulent_uplink_reciprocity.png"
FIELD_PNG = "examples/waveoptics/turbulent_uplink_reciprocity_field.png"


def build_scenario():
    '''Build the uplink scenario: a 100 mm waist from a 500 mm ground terminal.

    The pointing jitter is zero on both terminals. The coupled-flux model folds
    a jitter into the wander variance and the reciprocity overlap does not, so
    a jitter would compare two different things.
    '''
    ground = Terminal(aperture_m=GROUND_APERTURE_M, wavelength_m=WAVELENGTH_M,
                      pointing_jitter_rad=0.0,
                      transmitter=Transmitter(waist_m=WAIST_M, power_dbm=30.0))
    space = Terminal(aperture_m=SPACE_APERTURE_M, wavelength_m=WAVELENGTH_M,
                     pointing_jitter_rad=0.0)
    return SpaceScenario(ground=ground, space=space, direction="uplink",
                         channel=Channel(altitude_m=ALTITUDE_M))


def run_blocks(scenario, geometry, label, *, n_trials, block, seed, **kwargs):
    '''Run the trials in blocks, and print one progress line for each block.

    The runner keys its seeds on the trial INDEX, so a second call with the
    same seed repeats the same trials. Each block therefore takes its own seed,
    seed + block index. The set stays repeatable.
    '''
    trials, t0 = [], time.time()
    result = None
    for i, start in enumerate(range(0, n_trials, block)):
        n = min(block, n_trials - start)
        result = propagate_turbulent_scenario(scenario, geometry, n_trials=n,
                                              seed=seed + i, **kwargs)
        trials += result.trials
        print(f"    {label}: {len(trials):4d} / {n_trials} trials, "
              f"{time.time() - t0:6.1f} s")
    return dataclasses.replace(result, trials=trials)


def stats(loss_db):
    '''Give the mean, the standard deviation, and the 99% fade of a loss set.

    The 99% fade is the loss that 99 percent of the samples stay below: the
    99th percentile of the loss in dB.
    '''
    loss = np.asarray(loss_db, dtype=float)
    return (float(loss.mean()), float(loss.std(ddof=1)),
            float(np.percentile(loss, 99)))


def draw(cases):
    '''Draw the two loss distributions on one axis for each elevation.

    The histograms are normalised to a density, so the two sample counts do not
    have to match. The MEANS sit on the same axis as vertical lines, and so do
    the two 99% fades: the tail difference is the point of the figure.
    '''
    fig, axes = plt.subplots(1, len(cases), figsize=(6.4 * len(cases), 5.2),
                             constrained_layout=True, squeeze=False)
    for ax, case in zip(axes[0], cases):
        wave, flux = case["wave_db"], case["flux_db"]
        edges = np.linspace(min(wave.min(), flux.min()) - 1.0,
                            max(np.percentile(wave, 99.5),
                                np.percentile(flux, 99.5)) + 1.0, 40)
        ax.hist(flux, bins=edges, density=True, color="tab:red", alpha=0.45,
                label=f"fidelity 1, coupled flux, {flux.size} draws")
        ax.hist(wave, bins=edges, density=True, color="tab:blue", alpha=0.45,
                label=f"fidelity 2, field, {wave.size} snapshots")
        ax.axvline(case["flux_mean"], color="tab:red", linewidth=2.0,
                   label=f"fidelity 1 mean, {case['flux_mean']:.2f} dB")
        ax.axvline(case["wave_mean"], color="tab:blue", linewidth=2.0,
                   label=f"fidelity 2 mean, {case['wave_mean']:.2f} dB")
        ax.axvline(case["flux_q99"], color="tab:red", linestyle="--",
                   label=f"fidelity 1 99% fade, {case['flux_q99']:.2f} dB")
        ax.axvline(case["wave_q99"], color="tab:blue", linestyle="--",
                   label=f"fidelity 2 99% fade, {case['wave_q99']:.2f} dB")
        weak = ("weak" if case["weak_valid"]
                else "PAST the weak-fluctuation limit")
        ax.set_xlabel("uplink turbulence loss, dB")
        ax.set_ylabel("probability density")
        ax.set_title(f"{case['elevation_deg']:.0f} deg elevation\n"
                     f"coupled flux says sigma2_x = "
                     f"{case['sigma2_x']:.3f}, {weak}", fontsize=10)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)

    fig.suptitle("The uplink turbulence loss, two models: the reciprocity "
                 "overlap of a propagated field,\nand the Dios coupled-flux "
                 f"Monte Carlo. {ALTITUDE_M * 1e-3:.0f} km, waist "
                 f"{WAIST_M * 1e3:.0f} mm, {WAVELENGTH_M * 1e9:.0f} nm, "
                 "no pointing jitter", fontsize=12)
    fig.savefig(PNG, dpi=150)
    plt.close(fig)


def draw_field(F, aperture_m, obscuration, path, title):
    '''Draw the received-field amplitude and phase, with the aperture ring.

    The left panel is the amplitude |E| of one snapshot at the receive plane.
    The right panel is the phase; it is blanked where the intensity is below 2
    percent of the peak, so the phase of the dark background does not swamp the
    colour scale. The dashed ring is the receive aperture, and the dotted ring
    is the central obscuration when there is one. See Goodman, Introduction to
    Fourier Optics, ISBN 978-0974707723, for the amplitude and the phase of a
    scalar field.
    '''
    half_mm = F.siz / 2 * 1e3
    extent = [-half_mm, half_mm, -half_mm, half_mm]
    amp = np.abs(F.field)
    intensity = amp ** 2
    phase = np.where(intensity > 0.02 * intensity.max(),
                     np.angle(F.field), np.nan)

    fig, (a0, a1) = plt.subplots(1, 2, figsize=(11.6, 5.3),
                                 constrained_layout=True)
    im0 = a0.imshow(amp, extent=extent, origin="lower", cmap="inferno")
    fig.colorbar(im0, ax=a0, shrink=0.85, label="|E|, arb.")
    a0.set_title("amplitude |E|", fontsize=10)
    twilight = plt.cm.twilight.copy()
    twilight.set_bad("white")
    im1 = a1.imshow(phase, extent=extent, origin="lower", cmap=twilight,
                    vmin=-np.pi, vmax=np.pi)
    fig.colorbar(im1, ax=a1, shrink=0.85, label="phase, rad")
    a1.set_title("phase (blank below 2% intensity)", fontsize=10)

    r_mm = aperture_m / 2 * 1e3
    for ax in (a0, a1):
        ax.add_patch(plt.Circle((0, 0), r_mm, fill=False, edgecolor="cyan",
                                 linewidth=1.7, linestyle="--",
                                 label="receive aperture"))
        if obscuration > 0:
            ax.add_patch(plt.Circle((0, 0), obscuration * r_mm, fill=False,
                                    edgecolor="cyan", linewidth=1.4,
                                    linestyle=":"))
        ax.set_xlabel("x, mm")
        ax.set_ylabel("y, mm")
    a0.legend(loc="upper right", fontsize=8, framealpha=0.7)

    fig.suptitle(title, fontsize=12)
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main():
    t_start = time.time()
    scenario = build_scenario()
    hs = DEFAULT_HS
    cn2_profile = default_cn2_profile(scenario.channel.site, hs)

    print("=" * 78)
    print(f"{ALTITUDE_M * 1e-3:.0f} km uplink, the fidelity-2 reciprocity "
          f"overlap against the coupled-flux MC")
    print("=" * 78)
    print(f"  wavelength              {WAVELENGTH_M * 1e9:11.1f} nm")
    print(f"  transmit waist          {WAIST_M * 1e3:11.1f} mm")
    print(f"  ground aperture         {GROUND_APERTURE_M * 1e3:11.1f} mm "
          f"(alpha = {(GROUND_APERTURE_M / 2) / WAIST_M:.1f}, a clean launch)")
    print(f"  orbit                   {ALTITUDE_M * 1e-3:11.1f} km")
    print(f"  Cn2 profile             the default site, "
          f"ground {cn2_profile[0]:.3g} m^(-2/3)")
    print(f"  pointing jitter         {0.0:11.1f} urad, both terminals")
    print(f"  preset                  {PRESET:>11}")
    print(f"  snapshots per elevation {N_TRIALS:11d}")
    print(f"  coupled-flux draws      {FLUX_SAMPLES:11d}")
    print(f"  thread workers          {THREADER.max_workers:11d}")

    cases = []
    for elevation_deg in ELEVATIONS_DEG:
        orbit = CircularOrbit(ALTITUDE_M, elevation_deg=[elevation_deg])
        range_m = float(orbit.slant_range_m[0])
        with warnings.catch_warnings(record=True):
            warnings.simplefilter("always")
            grid, plan, report = turbulent_grid(scenario, orbit, preset=PRESET,
                                                hs=hs, cn2_profile=cn2_profile)
        print("")
        print(f"  elevation {elevation_deg:.0f} deg: slant range "
              f"{range_m * 1e-3:.1f} km; grid {grid.n} px, {grid.size_m:.3f} m; "
              f"{plan.z_m.size} screens; r0 = {plan.r0_total_m * 1e2:.2f} cm")
        print(f"    sampling: pixels per r0 {report.pixels_per_r0:.2f}, "
              f"grid margin {report.grid_margin:.2f}, Fresnel pixels "
              f"{report.fresnel_pixels_min:.2f}, "
              f"step/limit {report.step_over_limit_max:.3f}, "
              f"strongest screen {report.sigma2_r_screen_max:.4f}, "
              f"warnings {len(report.warnings)}")

        result = run_blocks(scenario, orbit, f"{elevation_deg:.0f} deg",
                            n_trials=N_TRIALS, block=BLOCK, seed=SEED,
                            preset=PRESET, hs=hs, cn2_profile=cn2_profile,
                            grid=grid, plan=plan, threader=THREADER)
        eta_turb = np.array([tr.eta_turb for tr in result.trials])
        times = np.array([tr.wall_time_s for tr in result.trials])
        wave_db = -10 * np.log10(eta_turb)

        # The SAME function that olb.links.uplink.uplink_turbulence_term calls.
        # It draws from the global numpy RNG, so seed that RNG to keep the run
        # repeatable. sigma_theta_rad stays 0: no jitter, see the docstring.
        np.random.seed(SEED % (2 ** 31))
        with warnings.catch_warnings(record=True):
            warnings.simplefilter("always")
            flux = _flux_result(WAIST_M, elevation_deg, range_m, WAVELENGTH_M,
                                hs, cn2_profile, 1.7e-14, FLUX_SAMPLES, 1,
                                sigma_theta_rad=0.0)
        flux_db = -10 * np.log10(flux["Is_summed"])

        wave_mean, wave_sigma, wave_q99 = stats(wave_db)
        flux_mean, flux_sigma, flux_q99 = stats(flux_db)
        cases.append({
            "elevation_deg": elevation_deg,
            "wave_db": wave_db, "flux_db": flux_db,
            "wave_mean": wave_mean, "wave_sigma": wave_sigma,
            "wave_q99": wave_q99,
            "flux_mean": flux_mean, "flux_sigma": flux_sigma,
            "flux_q99": flux_q99,
            "weak_valid": bool(flux["weak_fluctuation_valid"]),
            "sigma2_x": float(flux["sigma2_x_mean"]),
            "report": report,
            "mean_time_s": float(times.mean()),
            "min_time_s": float(times.min()),
            "max_time_s": float(times.max()),
        })

    # ---- the table ----
    print("")
    print("the uplink turbulence loss, in positive dB")
    print(f"{'elev':>6}{'model':>16}{'mean':>9}{'sigma':>9}{'99% fade':>11}"
          f"{'weak?':>8}")
    print("-" * 59)
    for case in cases:
        weak = "yes" if case["weak_valid"] else "NO"
        print(f"{case['elevation_deg']:>5.0f} {'fidelity 2 field':>16}"
              f"{case['wave_mean']:>9.2f}{case['wave_sigma']:>9.2f}"
              f"{case['wave_q99']:>11.2f}{'-':>8}")
        print(f"{'':>5} {'fidelity 1 flux':>16}"
              f"{case['flux_mean']:>9.2f}{case['flux_sigma']:>9.2f}"
              f"{case['flux_q99']:>11.2f}{weak:>8}")
        print(f"{'':>5} {'difference':>16}"
              f"{case['wave_mean'] - case['flux_mean']:>9.2f}"
              f"{case['wave_sigma'] - case['flux_sigma']:>9.2f}"
              f"{case['wave_q99'] - case['flux_q99']:>11.2f}{'':>8}")
    print("")
    gaps = [abs(c["wave_mean"] - c["flux_mean"]) for c in cases]
    print(f"  EXPECTED: the two MEANS inside 1 dB in the weak regime. "
          f"ACHIEVED: {min(gaps):.2f} to {max(gaps):.2f} dB.")
    print("  The 'weak?' column is weak_fluctuation_valid from the coupled-flux")
    print("  model. A NO there means that the fidelity-1 number is outside the")
    print("  regime of its own Rytov model, so the row is a report, not a test.")
    print("")
    print("  The TAILS are NOT compared. The two 99% fades come from two "
          "different")
    print("  constructions. The coupled flux draws ONE wander offset for each "
          "sample and")
    print("  puts a parametric lognormal irradiance about the local mean. The "
          "field")
    print("  simulation propagates a real wavefront, so its tail carries the "
          "speckle of")
    print("  a broken wavefront and the correlation between the wander and the")
    print("  scintillation. A field Monte Carlo therefore reaches deeper. "
          "Neither tail")
    print("  is validated here.")

    # ---- the sampling that the preset achieved ----
    print("")
    print("the wave-optics sampling report, ACHIEVED against the preset:")
    print(f"{'elev':>6}{'px per r0':>12}{'grid margin':>14}"
          f"{'Fresnel px':>13}{'step/limit':>13}{'warnings':>10}")
    print("-" * 68)
    for case in cases:
        rep = case["report"]
        print(f"{case['elevation_deg']:>5.0f} {rep.pixels_per_r0:>11.2f}"
              f"{rep.grid_margin:>14.2f}{rep.fresnel_pixels_min:>13.2f}"
              f"{rep.step_over_limit_max:>13.3f}{len(rep.warnings):>10d}")
    print("  The preset asks for 2 pixels per r0 and a grid margin of 2. A "
          "Fresnel pixel")
    print("  count of 2 or more is good. A step over the Forvard limit is cut "
          "into")
    print("  sub-steps by the split-step engine, so a value above 1 is not a "
          "fault.")

    # ---- the timing ----
    print("")
    print("per-trial wall time:")
    for case in cases:
        print(f"  {case['elevation_deg']:>3.0f} deg   "
              f"mean {case['mean_time_s']:6.3f} s   "
              f"min {case['min_time_s']:6.3f} s   "
              f"max {case['max_time_s']:6.3f} s")

    # ---- the phase screens ----
    # Each trial makes one fresh screen stack, so the run generates
    # screens-per-trial x trials x elevations screens.
    per_trial = int(plan.z_m.size)
    n_elev = len(ELEVATIONS_DEG)
    total_screens = per_trial * N_TRIALS * n_elev
    print("")
    print("phase screens created:")
    print(f"  {per_trial} per trial x {N_TRIALS} trials x {n_elev} elevations "
          f"= {total_screens} screens")

    draw(cases)

    # ---- the received field of one snapshot ----
    # The satellite sits outside the atmosphere, so the field to picture is the
    # DOWNLINK field at the ground: the uplink reads that same field through
    # reciprocity. Trial 0 of the lowest-elevation run, the strongest
    # turbulence.
    field_elev = ELEVATIONS_DEG[-1]
    F_rx, _, _ = propagate_turbulent_field(
        scenario, CircularOrbit(ALTITUDE_M, elevation_deg=[field_elev]),
        seed=SEED, trial=0, preset=PRESET, hs=hs, cn2_profile=cn2_profile)
    draw_field(F_rx, GROUND_APERTURE_M, 0.0, FIELD_PNG,
               f"Received field at the ground (downlink slab, read by "
               f"reciprocity), one snapshot, {field_elev:.0f} deg\n"
               f"{ALTITUDE_M * 1e-3:.0f} km uplink, ground aperture "
               f"{GROUND_APERTURE_M * 1e3:.0f} mm, {WAVELENGTH_M * 1e9:.0f} nm")

    print("")
    print(f"figure saved: {PNG}")
    print(f"field figure saved: {FIELD_PNG}")
    print(f"(elapsed {time.time() - t_start:.1f} s)")


if __name__ == "__main__":
    main()
