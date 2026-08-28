'''
The fidelity-2 turbulent split step on a 600 km downlink, against fidelity 0 and 1.

The script propagates a real complex field DOWN through the atmosphere slab of
a 600 km space link, one snapshot for each seed, at three elevations. It
compares two quantities with the models that the downlink budget already
carries:

- the aperture-averaged scintillation index of the collected power, against the
  fidelity-0 analytic plane-wave integral that `downlink_scintillation_term`
  calls (`olb.turbulence.plane_wave_scintillation`);
- the single-mode-fibre coupling loss, against the fidelity-1 FAST Term
  (`olb.models.fast.smf_fast_term`).

THE SLAB, NOT THE ORBIT. The satellite sits outside the atmosphere, so the
gridded path is the 20 km atmosphere slab only, and a unit PLANE WAVE enters at
the top of it. The 600 km of vacuum above the slab adds no turbulence. The
runner normalises the collected power against a zero-screen vacuum reference
through the SAME mask and the SAME hops, so `collected_power` has an exact
vacuum limit of 1.0 and it is a pure turbulence penalty.

WHAT TO EXPECT.

- The scintillation index agrees with the analytic aperture-averaged integral
  to about 20 percent. Both are the same physics on the same Cn2 profile.
- The fibre-coupling mean was expected to sit within 0.5 to 1 dB of FAST. It
  meets that at 30 degrees and it misses it at the zenith: the field simulation
  reads 0.7 dB less loss at 30 degrees and 2.9 dB less at the zenith. On the
  turbulence part alone (the total minus each static floor) the gap runs from
  0.2 dB to 2.0 dB. The script prints the gap and it names the candidates. It
  does not pick one, because an example is a demonstration, not a test.

WHY THE TWO FIBRE MODELS CAN DIFFER. They are two different models, not two
runs of one model:

- FAST samples the pupil on a very small grid (the script prints the achieved
  pixel count; it is 8 by default), and it models the log-AMPLITUDE as an
  aperture-averaged lognormal instead of propagating it. The field simulation
  propagates the amplitude and it samples the 500 mm pupil across about 90
  pixels.
- FAST optimises the fibre-mode radius (W0="opt"). The wave-optics overlap
  takes the fixed Ruilier radius D / 2.24, which is the optimum of an
  UNOBSCURED pupil. This case has a 0.3 central obscuration.
- The two static mode-match floors therefore differ, and the script prints both
  so that the reader can compare the TURBULENCE part alone.

The figure goes to `examples/waveoptics/figures/turbulent_downlink.png`.

The layer builds NO Term and it changes NO budget. See the README.

Run from the repo root:
    python -m examples.waveoptics.turbulent_downlink
'''

import dataclasses
import time
import warnings

import matplotlib.pyplot as plt
import numpy as np

from olb import SMF, CircularOrbit, Terminal, Transmitter
from olb.models.fast import _load_fast, smf_fast_term
from olb.scenario import Channel, SpaceScenario
from olb.turbulence.plane_wave_scintillation import \
    aperture_averaged_scintillation_index
from olb.turbulence.profiles import DEFAULT_HS, default_cn2_profile
from olb.waveoptics import Threader
from olb.waveoptics.turbulence import (propagate_turbulent_field,
                                       propagate_turbulent_scenario,
                                       turbulent_grid)

WAVELENGTH_M = 1550e-9
ALTITUDE_M = 600e3
GROUND_APERTURE_M = 0.50
GROUND_OBSCURATION = 0.30
SPACE_APERTURE_M = 0.10
SPACE_WAIST_M = 0.05
ELEVATIONS_DEG = (30.0, 60.0, 90.0)
PRESET = "rapid"
N_TRIALS = 70
BLOCK = 35                  # trials for each progress line
SEED = 20260826
FAST_SAMPLES = 20000

# The Cn2 of a still atmosphere, for the static mode-match floor. 1e-24 is four
# decades below any real site, so one trial gives the vacuum coupling.
QUIET_CN2 = 1e-24

# The trials are independent, so they run across threads. None takes one worker
# for each core. See olb.waveoptics.Threader.
THREADER = Threader()

PNG = "examples/waveoptics/figures/turbulent_downlink.png"
FIELD_PNG = "examples/waveoptics/figures/turbulent_downlink_field.png"


def build_scenario():
    '''Build the downlink scenario: a 500 mm obscured fibre receiver on the ground.

    The pointing jitter is zero on both terminals, so the comparison carries
    the TURBULENCE only.
    '''
    ground = Terminal(aperture_m=GROUND_APERTURE_M,
                      obscuration_ratio=GROUND_OBSCURATION,
                      wavelength_m=WAVELENGTH_M, pointing_jitter_rad=0.0,
                      detector=SMF(sensitivity_dbm=-110.0))
    space = Terminal(aperture_m=SPACE_APERTURE_M, wavelength_m=WAVELENGTH_M,
                     pointing_jitter_rad=0.0,
                     transmitter=Transmitter(waist_m=SPACE_WAIST_M,
                                             power_dbm=30.0))
    return SpaceScenario(ground=ground, space=space, direction="downlink",
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


def verdict(ratio, low, high):
    '''Give PASS when the ratio sits inside the band, and CHECK otherwise.'''
    return "PASS" if low <= ratio <= high else "CHECK"


def loss_db(values):
    '''Give the mean loss in positive dB of a set of efficiency samples.'''
    return float(-10 * np.log10(np.mean(values)))


def static_floor_db(scenario, hs):
    '''Give the still-atmosphere fibre-coupling loss, in dB.

    One trial with a negligible Cn2 gives the STATIC mode-match floor of the
    wave-optics overlap: the loss of a flat wavefront on the obscured pupil
    against the fixed Ruilier fibre mode. The turbulence part of any other
    number is that number minus this floor.
    '''
    quiet = np.full(hs.size, QUIET_CN2)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        result = propagate_turbulent_scenario(
            scenario, CircularOrbit(ALTITUDE_M, elevation_deg=[90.0]),
            n_trials=1, seed=1, preset=PRESET, hs=hs, cn2_profile=quiet)
    return -10 * np.log10(result.trials[0].smf_eta)


def draw(rows, floor_wave_db, floor_fast_db, has_fast):
    '''Draw the scintillation index and the coupling loss against the elevation.

    The left panel puts the measured aperture index on the analytic curve. The
    right panel puts the measured fibre-coupling loss, with the standard error
    of its mean as the error bar, on the FAST mean and the FAST 99% fade. The
    two static floors sit on the same axis, so the reader sees how much of each
    total is the mode match and how much is the turbulence.
    '''
    fig, axes = plt.subplots(1, 2, figsize=(12.8, 5.2),
                             constrained_layout=True)

    elev = np.array([r["elevation_deg"] for r in rows])
    fine = np.linspace(min(elev) - 3, 92, 120)

    axes[0].semilogy(fine, rows[0]["analytic_curve"](fine), color="tab:red",
                     linewidth=2.0,
                     label="fidelity 0, the aperture-averaged integral")
    axes[0].semilogy(elev, [r["sigma2_wave"] for r in rows], "o",
                     color="tab:blue", markersize=9,
                     label=f"fidelity 2, {N_TRIALS} snapshots each")
    axes[0].set_xlabel("elevation, deg")
    axes[0].set_ylabel("sigma2_P, the aperture flux index")
    axes[0].set_title(f"1. the bucket scintillation index\n"
                      f"D = {GROUND_APERTURE_M * 1e3:.0f} mm, obscuration "
                      f"{GROUND_OBSCURATION:g}", fontsize=10)
    axes[0].grid(alpha=0.3, which="both")
    axes[0].legend(fontsize=8)

    axes[1].errorbar(elev, [r["smf_wave_db"] for r in rows],
                     yerr=[r["smf_wave_sem_db"] for r in rows], fmt="o",
                     color="tab:blue", markersize=9, capsize=4,
                     label="fidelity 2, the mean coupling loss")
    axes[1].plot(elev, [r["q99_wave_db"] for r in rows], "v",
                 color="tab:blue", markersize=8, alpha=0.6,
                 label="fidelity 2, the 99% fade")
    axes[1].axhline(floor_wave_db, color="tab:blue", linestyle=":",
                    label=f"fidelity 2 static floor, {floor_wave_db:.2f} dB")
    if has_fast:
        axes[1].plot(elev, [r["smf_fast_db"] for r in rows], "s--",
                     color="tab:red", markersize=8,
                     label="fidelity 1 FAST, the mean coupling loss")
        axes[1].plot(elev, [r["q99_fast_db"] for r in rows], "v--",
                     color="tab:red", markersize=8, alpha=0.6,
                     label="fidelity 1 FAST, the 99% fade")
        axes[1].axhline(floor_fast_db, color="tab:red", linestyle=":",
                        label=f"FAST static floor, {floor_fast_db:.2f} dB")
    axes[1].invert_yaxis()
    axes[1].set_xlabel("elevation, deg")
    axes[1].set_ylabel("fibre-coupling loss, dB")
    axes[1].set_title("2. the single-mode-fibre coupling\n"
                      "two models, not two runs of one model", fontsize=10)
    axes[1].grid(alpha=0.3)
    # The band between the static floors and the mean losses is empty, so the
    # legend goes there and it covers no marker.
    axes[1].legend(fontsize=7, loc="upper center", ncol=2)

    fig.suptitle("Fidelity-2 turbulent split step against fidelity 0 and 1, "
                 f"{ALTITUDE_M * 1e-3:.0f} km downlink\n"
                 f"ground aperture {GROUND_APERTURE_M * 1e3:.0f} mm, "
                 f"obscuration {GROUND_OBSCURATION:g}, single-mode fibre, "
                 f"{WAVELENGTH_M * 1e9:.0f} nm, default site", fontsize=12)
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

    try:
        _load_fast()
        has_fast = True
    except ImportError as exc:
        has_fast = False
        print(f"NOTE: fast-aosim is absent ({exc.__class__.__name__}). The "
              f"fidelity-1 columns are skipped.")

    print("=" * 78)
    print(f"{ALTITUDE_M * 1e-3:.0f} km downlink, fidelity-2 turbulent split "
          f"step against fidelity 0 and 1")
    print("=" * 78)
    print(f"  wavelength              {WAVELENGTH_M * 1e9:11.1f} nm")
    print(f"  ground aperture         {GROUND_APERTURE_M * 1e3:11.1f} mm, "
          f"obscuration {GROUND_OBSCURATION:g}, SMF")
    print(f"  orbit                   {ALTITUDE_M * 1e-3:11.1f} km")
    print(f"  Cn2 profile             the default site, "
          f"ground {cn2_profile[0]:.3g} m^(-2/3)")
    print(f"  preset                  {PRESET:>11}")
    print(f"  snapshots per elevation {N_TRIALS:11d}")
    print(f"  thread workers          {THREADER.max_workers:11d}")

    floor_wave_db = static_floor_db(scenario, hs)
    print(f"  static mode-match floor {floor_wave_db:11.3f} dB "
          f"(fidelity 2, one still-atmosphere trial)")

    rows = []
    floor_fast_db = float("nan")
    for elevation_deg in ELEVATIONS_DEG:
        orbit = CircularOrbit(ALTITUDE_M, elevation_deg=[elevation_deg])
        with warnings.catch_warnings(record=True):
            warnings.simplefilter("always")
            grid, plan, report = turbulent_grid(scenario, orbit, preset=PRESET,
                                                hs=hs, cn2_profile=cn2_profile)
        print("")
        print(f"  elevation {elevation_deg:.0f} deg: grid {grid.n} px, "
              f"{grid.size_m:.3f} m; {plan.z_m.size} screens; "
              f"r0 = {plan.r0_total_m * 1e2:.2f} cm; "
              f"slab = {plan.z_total_m * 1e-3:.1f} km")
        print(f"    sampling: pixels per r0 {report.pixels_per_r0:.2f}, "
              f"grid margin {report.grid_margin:.2f}, Fresnel pixels "
              f"{report.fresnel_pixels_min:.2f}, "
              f"step/limit {report.step_over_limit_max:.3f}, "
              f"warnings {len(report.warnings)}")
        result = run_blocks(scenario, orbit, f"{elevation_deg:.0f} deg",
                            n_trials=N_TRIALS, block=BLOCK, seed=SEED,
                            preset=PRESET, hs=hs, cn2_profile=cn2_profile,
                            grid=grid, plan=plan, threader=THREADER)

        power = np.array([tr.collected_power for tr in result.trials])
        smf_eta = np.array([tr.smf_eta for tr in result.trials])
        times = np.array([tr.wall_time_s for tr in result.trials])

        sigma2_wave = float(power.var() / power.mean() ** 2)
        # The SAME kernel that olb.links.downlink.downlink_scintillation_term
        # calls for its aperture-averaged flux index. See Andrews and Phillips,
        # DOI 10.1117/3.626196, Ch. 10 and Ch. 12.
        sigma2_analytic = float(aperture_averaged_scintillation_index(
            GROUND_APERTURE_M, elevation_deg, WAVELENGTH_M, hs, cn2_profile))

        row = {
            "elevation_deg": elevation_deg,
            "sigma2_wave": sigma2_wave,
            "sigma2_analytic": sigma2_analytic,
            "smf_wave_db": loss_db(smf_eta),
            # The standard error of the mean efficiency, in dB.
            "smf_wave_sem_db": float(10.0 / np.log(10.0)
                                     * smf_eta.std(ddof=1)
                                     / np.sqrt(smf_eta.size) / smf_eta.mean()),
            "q99_wave_db": float(-10 * np.log10(np.percentile(smf_eta, 1))),
            "mean_time_s": float(times.mean()),
            "min_time_s": float(times.min()),
            "max_time_s": float(times.max()),
            "analytic_curve": lambda e: aperture_averaged_scintillation_index(
                GROUND_APERTURE_M, e, WAVELENGTH_M, hs, cn2_profile),
        }
        if has_fast:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                term = smf_fast_term(scenario,
                                     CircularOrbit(ALTITUDE_M,
                                                   elevation_deg=elevation_deg),
                                     hs=hs, cn2_profile=cn2_profile,
                                     n_samples=FAST_SAMPLES)
            floor_fast_db = float(term.meta["floor_db"])
            row["smf_fast_db"] = float(term.mean_db)
            row["q99_fast_db"] = float(term.quantile_db(0.99))
            row["fast_npxls"] = int(term.meta["npxls"])
            row["fast_weak"] = bool(term.meta["amplitude_regime_weak"])
        rows.append(row)

    # ---- table 1: the bucket scintillation ----
    print("")
    print("table 1: the aperture-averaged scintillation index of the collected "
          "power")
    print(f"{'elev':>6}{'fidelity 2':>13}{'fidelity 0':>13}{'ratio':>9}"
          f"{'':>8}")
    print("-" * 49)
    for row in rows:
        ratio = row["sigma2_wave"] / row["sigma2_analytic"]
        print(f"{row['elevation_deg']:>5.0f}{'':>1}{row['sigma2_wave']:>13.5f}"
              f"{row['sigma2_analytic']:>13.5f}{ratio:>9.3f}   "
              f"{verdict(ratio, 0.70, 1.35)}")
    print(f"  The MC error of a variance from {N_TRIALS} snapshots is about "
          f"{100 * np.sqrt(2.0 / N_TRIALS):.0f} percent.")

    # ---- table 2: the fibre coupling ----
    print("")
    print("table 2: the single-mode-fibre coupling loss, in positive dB")
    if has_fast:
        print(f"{'elev':>6}{'fid 2':>9}{'FAST':>9}{'diff':>8}"
              f"{'fid 2 turb':>12}{'FAST turb':>11}{'diff':>8}"
              f"{'fid 2 99%':>11}{'FAST 99%':>10}")
        print("-" * 84)
        for row in rows:
            turb_wave = row["smf_wave_db"] - floor_wave_db
            turb_fast = row["smf_fast_db"] - floor_fast_db
            print(f"{row['elevation_deg']:>5.0f}{'':>1}"
                  f"{row['smf_wave_db']:>9.2f}{row['smf_fast_db']:>9.2f}"
                  f"{row['smf_wave_db'] - row['smf_fast_db']:>8.2f}"
                  f"{turb_wave:>12.2f}{turb_fast:>11.2f}"
                  f"{turb_wave - turb_fast:>8.2f}"
                  f"{row['q99_wave_db']:>11.2f}{row['q99_fast_db']:>10.2f}")
        print(f"  static floors: fidelity 2 {floor_wave_db:.2f} dB, "
              f"FAST {floor_fast_db:.2f} dB. The 'turb' columns are the total")
        print("  minus the floor, so they hold the turbulence penalty alone.")
        print(f"  FAST samples the pupil across {rows[0]['fast_npxls']} pixels; "
              f"the field simulation uses about "
              f"{GROUND_APERTURE_M / (grid.size_m / grid.n):.0f}.")
        gaps = [abs(r["smf_wave_db"] - r["smf_fast_db"]) for r in rows]
        turb_gaps = [abs((r["smf_wave_db"] - floor_wave_db)
                         - (r["smf_fast_db"] - floor_fast_db)) for r in rows]
        print("")
        print(f"  EXPECTED: the two means inside 1 dB. ACHIEVED: {min(gaps):.1f} "
              f"to {max(gaps):.1f} dB on the")
        print(f"  totals, and {min(turb_gaps):.1f} to {max(turb_gaps):.1f} dB "
              f"on the turbulence part alone. The gap is")
        print("  smallest at the LOW elevation, where the turbulence "
              "dominates, and largest")
        print("  at the zenith, where the static floor and the fibre-mode "
              "choice matter most.")
        print("  The two models differ in the pupil sampling, in the amplitude "
              "treatment")
        print("  (FAST puts an aperture-averaged lognormal on the amplitude; "
              "the field")
        print("  simulation propagates it), and in the fibre-mode radius (FAST "
              "optimises it;")
        print("  the wave-optics overlap takes the fixed Ruilier D / 2.24 of "
              "an UNOBSCURED")
        print("  pupil). The script does not pick one. Which model is right is "
              "an owner")
        print("  decision, and it is the reason that no fidelity-2 Term is "
              "wired.")
    else:
        print(f"{'elev':>6}{'fid 2':>9}{'fid 2 turb':>12}{'fid 2 99%':>11}")
        print("-" * 38)
        for row in rows:
            print(f"{row['elevation_deg']:>5.0f}{'':>1}"
                  f"{row['smf_wave_db']:>9.2f}"
                  f"{row['smf_wave_db'] - floor_wave_db:>12.2f}"
                  f"{row['q99_wave_db']:>11.2f}")
        print("  fast-aosim is absent, so the fidelity-1 columns are missing.")

    # ---- the timing ----
    print("")
    print("per-trial wall time:")
    for row in rows:
        print(f"  {row['elevation_deg']:>3.0f} deg   "
              f"mean {row['mean_time_s']:6.3f} s   "
              f"min {row['min_time_s']:6.3f} s   "
              f"max {row['max_time_s']:6.3f} s")

    # ---- the phase screens ----
    # Each trial makes one fresh screen stack, so the run generates
    # screens-per-trial x trials x elevations screens, plus the one static-floor
    # trial.
    per_trial = int(plan.z_m.size)
    n_elev = len(ELEVATIONS_DEG)
    total_screens = per_trial * N_TRIALS * n_elev + per_trial
    print("")
    print("phase screens created:")
    print(f"  {per_trial} per trial x {N_TRIALS} trials x {n_elev} elevations "
          f"+ {per_trial} (static floor) = {total_screens} screens")

    draw(rows, floor_wave_db, floor_fast_db, has_fast)

    # ---- the received field of one snapshot ----
    # Trial 0 of the lowest-elevation run, the strongest turbulence. The seed
    # and the trial index match that run, so this IS its first snapshot.
    field_elev = ELEVATIONS_DEG[0]
    F_rx, _, _ = propagate_turbulent_field(
        scenario, CircularOrbit(ALTITUDE_M, elevation_deg=[field_elev]),
        seed=SEED, trial=0, preset=PRESET, hs=hs, cn2_profile=cn2_profile)
    draw_field(F_rx, GROUND_APERTURE_M, GROUND_OBSCURATION, FIELD_PNG,
               f"Received downlink field at the ground, one snapshot, "
               f"{field_elev:.0f} deg elevation\n"
               f"{ALTITUDE_M * 1e-3:.0f} km, ground aperture "
               f"{GROUND_APERTURE_M * 1e3:.0f} mm, obscuration "
               f"{GROUND_OBSCURATION:g}, {WAVELENGTH_M * 1e9:.0f} nm")

    print("")
    print(f"figure saved: {PNG}")
    print(f"field figure saved: {FIELD_PNG}")
    print(f"(elapsed {time.time() - t_start:.1f} s)")


if __name__ == "__main__":
    main()
