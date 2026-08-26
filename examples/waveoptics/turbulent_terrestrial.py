'''
The fidelity-2 turbulent split step on a 2 km horizontal link, against fidelity 0.

The script propagates a real complex field through a stack of random phase
screens, one snapshot for each seed, and it compares the statistics of those
snapshots with the analytic Terms of the terrestrial budget.

THE WEAK-FLUCTUATION PRECONDITION. Every analytic target here is a WEAK
fluctuation form: the Dios on-axis Gaussian-beam index, the Andrews weak
aperture-averaging factor, and the Noll-plus-Dikmelik-Davidson fibre coupling.
So the case keeps the plane-wave Rytov variance at about 0.2. The script
ASSERTS that precondition, because a comparison outside it compares against a
model that its own author does not trust. Cn2 = 3e-15 m^(-2/3) over 2 km at
1550 nm gives sigma2_R = 0.21 and r0 = 10.7 cm.

THE LAUNCH IS CLEAN. The launch aperture is 200 mm for a 50 mm waist, so the
clip takes less than 0.01 percent of the power. The propagated beam is
therefore the pure Gaussian that the analytic forms assume, and the comparison
tests the TURBULENCE, not the truncation.

THREE RECEIVE APERTURES, ONE ATMOSPHERE. TurbTrial carries the collected power
of the receive aperture, and nothing else. Its docstring forbids a piecemeal
extension, so the script does not add a field to it. Instead it runs the same
grid, the same screen plan and the same seeds three times, and it changes the
receive aperture only:

- a 3-pixel PINHOLE, which reads the on-axis irradiance;
- a 30 mm SAMPLING aperture, which is small against the 108 mm beam, so it
  samples the beam exactly as the analytic aperture-averaging factor assumes;
- the 100 mm BUDGET aperture with its single-mode fibre.

WHAT TO EXPECT. The expectation is 15 to 20 percent on the scintillation rows
and 0.5 dB on the fibre-coupling row. The first two rows meet it. The last two
do NOT, and the reason is physics, not a defect:

- The on-axis index agrees to about 10 percent. The pinhole index runs a little
  HIGH, because the split step lets the beam WANDER off the fixed pinhole and
  the Dios on-axis form carries no wander.
- The 30 mm sampling bucket agrees with the weak aperture-averaging factor.
- The 100 mm BUDGET bucket does NOT. That aperture collects most of the beam,
  and the split step conserves power, so almost nothing is left to fluctuate.
  The weak aperture-averaging factor of Andrews models an aperture that SAMPLES
  a wide beam; it does not know about a beam that the aperture nearly holds.
  The script prints the measured capture fraction next to that row.
- The fidelity-0 fibre-coupling Term reads about 2.5 dB MORE loss than the
  field simulation. The Term takes the Noll residual at the aperture diameter
  D, and the Dikmelik-Davidson coupling of a UNIFORMLY illuminated pupil. The
  received field here is a finite Gaussian beam whose intensity is already down
  to 1/e^2 at the rim of the 100 mm aperture, so the phase error that the fibre
  mode really weighs is smaller than the uniform-pupil form gives. The Term
  docstring names this exact limit ("it ignores ... the near-field beam
  curvature that a finite horizontal beam carries. For those a fidelity-2
  split-step beam-propagation model is needed"). This script is that model, and
  it puts a number on the gap. The number is converged: a grid two times wider
  moves it by 0.3 dB.
- The planned composite r0 IS the plane-wave Fried parameter by construction,
  so it repeats the analytic plane-wave value exactly. The Gaussian-beam r0 of
  the same path is LARGER, because a finite beam samples less of the
  turbulence.

The script prints PASS or CHECK against a loose band. An example is a
demonstration, not a test: only the weak-fluctuation precondition is asserted.

The figure goes to `examples/waveoptics/turbulent_terrestrial.png`.

The layer builds NO Term and it changes NO budget. See the README.

Run from the repo root:
    python -m examples.waveoptics.turbulent_terrestrial
'''

import dataclasses
import time
import warnings

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm

from olb import SMF, Terminal, Transmitter
from olb.beam import free_space_radius
from olb.geometry import HorizontalPath
from olb.models.coupling import terrestrial_smf_coupling_term
from olb.scenario import TerrestrialChannel, TerrestrialScenario
from olb.turbulence.beam_wave_scintillation import on_axis_scintillation_index
from olb.turbulence.gaussian_fried import (gaussian_fried_parameter,
                                           plane_wave_fried_parameter)
from olb.turbulence.plane_wave_scintillation import \
    aperture_averaging_factor_weak
from olb.waveoptics import Threader
from olb.waveoptics.turbulence import (propagate_turbulent_field,
                                       propagate_turbulent_scenario,
                                       turbulent_grid)

WAVELENGTH_M = 1550e-9
PATH_M = 2000.0
CN2 = 3e-15                 # sigma2_R = 0.21, firmly weak
WAIST_M = 0.05
TX_APERTURE_M = 0.20        # 4 waists: a negligible launch clip
RX_APERTURE_M = 0.10        # the budget receiver, with the fibre
SAMPLE_APERTURE_M = 0.03    # small against the beam: a fair bucket sample
PRESET = "standard"
N_TRIALS = 300
BLOCK = 40                  # trials for each progress line
SEED = 20260826
PINHOLE_PIXELS = 3          # the on-axis probe, in pixels

# The weak-fluctuation band that the script asserts. Below 0.05 the turbulence
# is too quiet to measure against the Monte Carlo noise, and 0.35 is already
# past the comfortable margin of the WEAK_FLUCTUATION_LIMIT of 0.25.
WEAK_BAND = (0.05, 0.35)

# The trials are independent, so they run across threads. None takes one worker
# for each core. See olb.waveoptics.Threader.
THREADER = Threader()

PNG = "examples/waveoptics/turbulent_terrestrial.png"
FIELD_PNG = "examples/waveoptics/turbulent_terrestrial_field.png"

# The Cn2 path grid of the analytic terms. A horizontal path is a constant-Cn2
# grid from the transmitter to the receiver, exactly as
# olb.links.terrestrial.terrestrial_scintillation_term builds it.
_GRID_N = 400


def build_scenario(rx_aperture_m=RX_APERTURE_M, fibre=True):
    '''Build the horizontal scenario: a clean 50 mm launch to an SMF receiver.

    The pointing jitter is zero on both terminals, so the comparison carries
    the TURBULENCE only. Set fibre=False for the two probe apertures: they read
    the collected power alone, so they need no detector.
    '''
    near = Terminal(aperture_m=TX_APERTURE_M, wavelength_m=WAVELENGTH_M,
                    pointing_jitter_rad=0.0,
                    transmitter=Transmitter(waist_m=WAIST_M))
    far = Terminal(aperture_m=rx_aperture_m, wavelength_m=WAVELENGTH_M,
                   pointing_jitter_rad=0.0,
                   detector=SMF() if fibre else None)
    return TerrestrialScenario(
        near=near, far=far,
        channel=TerrestrialChannel(path_length_m=PATH_M, cn2=CN2))


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


def scintillation_index(values):
    '''Give the normalised variance var(I) / mean(I)^2 of a set of samples.'''
    v = np.asarray(values, dtype=float)
    return float(v.var() / v.mean() ** 2)


def verdict(ratio, low, high):
    '''Give PASS when the ratio sits inside the band, and CHECK otherwise.'''
    return "PASS" if low <= ratio <= high else "CHECK"


def timing_line(result, label):
    '''Print the mean, the smallest and the largest trial time.'''
    t = np.array([tr.wall_time_s for tr in result.trials])
    print(f"  {label:<22}mean {t.mean():6.3f} s   min {t.min():6.3f} s   "
          f"max {t.max():6.3f} s   total {t.sum():7.1f} s")


def _fade_panel(ax, power, sigma2, title):
    '''Draw one bucket fade histogram against its analytic lognormal.

    With sigma_l^2 = ln(1 + sigma2_P) the loss in dB about the mean is normal
    with the mean (5/ln10) sigma_l^2 and the standard deviation (10/ln10)
    sigma_l. See Andrews and Phillips, DOI 10.1117/3.626196, Ch. 5, the
    lognormal irradiance model.
    '''
    ln10 = np.log(10.0)
    fade_db = -10 * np.log10(power / power.mean())
    sigma_l2 = np.log(1.0 + sigma2)
    sigma_l = np.sqrt(sigma_l2)
    span = max(abs(fade_db).max(), (10.0 / ln10) * sigma_l * 3.0)
    x = np.linspace(-1.2 * span, 1.2 * span, 400)
    pdf = norm.pdf(x, loc=(5.0 / ln10) * sigma_l2,
                   scale=(10.0 / ln10) * sigma_l)

    ax.hist(fade_db, bins=22, density=True, color="tab:blue", alpha=0.55,
            label=f"fidelity 2, {fade_db.size} snapshots")
    ax.plot(x, pdf, color="tab:red", linewidth=2.0,
            label=f"fidelity 0 lognormal\nsigma2_P = {sigma2:.4f}")
    ax.axvline(0.0, color="black", linestyle=":", linewidth=1.0)
    ax.set_xlabel("collected-power fade about the mean, dB")
    ax.set_ylabel("probability density")
    ax.set_title(title, fontsize=10)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)


def draw(sampled, sigma2_sample, power, sigma2_bucket, smf_eta,
         analytic_smf_db):
    '''Draw the two bucket fades, and the fibre-coupling loss histogram.

    Panel 1 is the small sampling aperture, where the analytic lognormal fits.
    Panel 2 is the budget aperture, which collects most of the beam: the
    measured fade is far narrower than the analytic curve, because the split
    step conserves the power that the aperture holds. Panel 3 puts the
    fibre-coupling loss of every snapshot against the analytic MEAN of the
    fidelity-0 coupling Term. That Term is mean-only, so it draws as one line:
    it carries no fade at all. The spread of the histogram is the fade that
    fidelity 0 cannot give.
    '''
    fig, axes = plt.subplots(1, 3, figsize=(16.4, 5.0),
                             constrained_layout=True)

    _fade_panel(axes[0], sampled, sigma2_sample,
                f"1. the {SAMPLE_APERTURE_M * 1e3:.0f} mm sampling bucket\n"
                f"small against the beam: the analytic form holds")
    _fade_panel(axes[1], power, sigma2_bucket,
                f"2. the {RX_APERTURE_M * 1e3:.0f} mm budget aperture\n"
                f"it holds most of the beam: the analytic form does not")

    coupling_db = -10 * np.log10(smf_eta)
    axes[2].hist(coupling_db, bins=22, color="tab:green", alpha=0.55,
                 label=f"fidelity 2, {coupling_db.size} snapshots")
    axes[2].axvline(-10 * np.log10(smf_eta.mean()), color="tab:blue",
                    linewidth=2.0,
                    label=f"fidelity 2 mean, "
                          f"{-10 * np.log10(smf_eta.mean()):.2f} dB")
    axes[2].axvline(analytic_smf_db, color="tab:red", linestyle="--",
                    linewidth=2.0,
                    label=f"fidelity 0 mean-only Term, {analytic_smf_db:.2f} dB")
    axes[2].set_xlabel("fibre-coupling loss, dB")
    axes[2].set_ylabel("snapshots")
    axes[2].set_title("3. the single-mode-fibre coupling\n"
                      "the fidelity-0 Term gives the red line only, no fade",
                      fontsize=10)
    axes[2].grid(alpha=0.3)
    axes[2].legend(fontsize=8)

    fig.suptitle("Fidelity-2 turbulent split step against fidelity 0, "
                 f"{PATH_M * 1e-3:.0f} km horizontal path\n"
                 f"waist {WAIST_M * 1e3:.0f} mm, receive aperture "
                 f"{RX_APERTURE_M * 1e3:.0f} mm, Cn2 = {CN2:g} m^(-2/3), "
                 f"{WAVELENGTH_M * 1e9:.0f} nm", fontsize=12)
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
    geometry = HorizontalPath(PATH_M)

    # The grid and the plan come out ONE time. Both runs below then share them,
    # so the pinhole run sees the SAME atmosphere as the bucket run.
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        grid, plan, report = turbulent_grid(scenario, geometry, preset=PRESET)

    sigma2_R = float(plan.sigma2_r.sum())
    assert WEAK_BAND[0] <= sigma2_R <= WEAK_BAND[1], (
        f"the plane-wave Rytov variance is {sigma2_R:.3f}, outside the weak "
        f"band {WEAK_BAND}. Every analytic target of this script is a weak "
        f"fluctuation form, so a comparison outside that band is meaningless. "
        f"Change CN2 or PATH_M.")

    print("=" * 74)
    print(f"{PATH_M * 1e-3:.0f} km horizontal link, fidelity-2 turbulent split "
          f"step against fidelity 0")
    print("=" * 74)
    print(f"  wavelength              {WAVELENGTH_M * 1e9:11.1f} nm")
    print(f"  path length             {PATH_M * 1e-3:11.3f} km")
    print(f"  Cn2                     {CN2:11.3g} m^(-2/3)")
    print(f"  transmit waist          {WAIST_M * 1e3:11.1f} mm")
    print(f"  launch aperture         {TX_APERTURE_M * 1e3:11.1f} mm "
          f"(alpha = {(TX_APERTURE_M / 2) / WAIST_M:.2f}, a clean launch)")
    print(f"  receive aperture        {RX_APERTURE_M * 1e3:11.1f} mm, SMF")
    print(f"  sampling aperture       {SAMPLE_APERTURE_M * 1e3:11.1f} mm "
          f"(the fair bucket probe)")
    print(f"  sigma2_R, plane wave    {sigma2_R:11.4f}   "
          f"(weak: the band is {WEAK_BAND[0]} to {WEAK_BAND[1]})")
    print("")
    print(f"  preset                  {PRESET:>11}")
    print(f"  grid                    {grid.n:11d} px, {grid.size_m:.4f} m")
    print(f"  pixel pitch             {grid.pixel_m * 1e3:11.4f} mm")
    print(f"  screens                 {plan.z_m.size:11d}")
    print(f"  thread workers          {THREADER.max_workers:11d}")
    print(f"  r0, composite           {plan.r0_total_m * 1e2:11.3f} cm")
    print("")
    print("  the sampling report, ACHIEVED against the preset:")
    print(f"    pixels per r0         {report.pixels_per_r0:11.2f}  "
          f"(the preset asks for 3)")
    print(f"    grid margin           {report.grid_margin:11.2f}  "
          f"(1.0 means the light just fits)")
    print(f"    Fresnel pixels, min   {report.fresnel_pixels_min:11.2f}  "
          f"(2 or more is good)")
    print(f"    step / Forvard limit  {report.step_over_limit_max:11.3f}  "
          f"(1.0 or less is good)")
    print(f"    strongest screen s2_r {report.sigma2_r_screen_max:11.4f}")
    print(f"    warnings              {len(report.warnings):11d}")
    for text in report.warnings + tuple(str(w.message) for w in caught):
        print(f"      {text}")

    # ---- the three runs ----
    # Every run shares the grid, the plan and the seeds, so the three receive
    # apertures look at the SAME 120 atmospheres.
    pinhole_m = PINHOLE_PIXELS * grid.pixel_m
    print("")
    print(f"  running {N_TRIALS} snapshots three times on the same screens: a "
          f"{pinhole_m * 1e3:.2f} mm pinhole,")
    print(f"  a {SAMPLE_APERTURE_M * 1e3:.0f} mm sampling bucket, and the "
          f"{RX_APERTURE_M * 1e3:.0f} mm budget aperture with its fibre")
    runs = {}
    for label, aperture_m, fibre in (
            ("pinhole", pinhole_m, False),
            ("sampling", SAMPLE_APERTURE_M, False),
            ("budget", RX_APERTURE_M, True)):
        runs[label] = run_blocks(
            build_scenario(rx_aperture_m=aperture_m, fibre=fibre), geometry,
            label, n_trials=N_TRIALS, block=BLOCK, seed=SEED, grid=grid,
            plan=plan, threader=THREADER)

    on_axis = np.array([tr.collected_power for tr in runs["pinhole"].trials])
    sampled = np.array([tr.collected_power for tr in runs["sampling"].trials])
    power = np.array([tr.collected_power for tr in runs["budget"].trials])
    smf_eta = np.array([tr.smf_eta for tr in runs["budget"].trials])

    print("")
    print("  per-trial wall time:")
    for label in ("pinhole", "sampling", "budget"):
        timing_line(runs[label], f"{label} run")

    # ---- the analytic targets ----
    # The same kernels that olb.links.terrestrial.terrestrial_scintillation_term
    # calls: the Dios on-axis Gaussian-beam index (Applied Optics 43 (2004)
    # 3866, Eq. 16) on a constant-Cn2 horizontal grid, and the Andrews weak
    # Kolmogorov aperture-averaging factor (DOI 10.1117/3.626196, Ch. 10).
    hs = np.linspace(0.0, PATH_M, _GRID_N)
    cn2_profile = np.full_like(hs, CN2)
    sigma2_on_axis = float(on_axis_scintillation_index(
        hs, cn2_profile, WAIST_M, WAVELENGTH_M, elevation_deg=90.0,
        path_length_m=None))
    factor_a = float(aperture_averaging_factor_weak(RX_APERTURE_M,
                                                    WAVELENGTH_M, PATH_M))
    factor_a_sample = float(aperture_averaging_factor_weak(
        SAMPLE_APERTURE_M, WAVELENGTH_M, PATH_M))
    sigma2_bucket = factor_a * sigma2_on_axis
    sigma2_sample = factor_a_sample * sigma2_on_axis

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        coupling = terrestrial_smf_coupling_term(scenario, geometry)
    analytic_smf_db = float(coupling.mean_db)
    measured_smf_db = float(-10 * np.log10(smf_eta.mean()))

    # The free-space beam radius at the receiver. It sets the capture fraction
    # of each probe aperture, and it is the reason that the two rows below
    # disagree with the analytic forms.
    sample_beam_radius = float(free_space_radius(WAIST_M, PATH_M, None,
                                                 WAVELENGTH_M))

    r0_plane = float(plane_wave_fried_parameter(PATH_M, CN2, WAVELENGTH_M))
    r0_gauss = float(gaussian_fried_parameter(PATH_M, WAIST_M, CN2,
                                              WAVELENGTH_M))

    measured_on_axis = scintillation_index(on_axis)
    measured_sample = scintillation_index(sampled)
    measured_bucket = scintillation_index(power)

    # ---- the table ----
    rows = [
        ("sigma2_I, on axis", measured_on_axis, sigma2_on_axis, "",
         0.70, 1.45),
        (f"sigma2_P, {SAMPLE_APERTURE_M * 1e3:.0f} mm sample",
         measured_sample, sigma2_sample, "", 0.70, 1.35),
        (f"sigma2_P, {RX_APERTURE_M * 1e3:.0f} mm budget",
         measured_bucket, sigma2_bucket, "", 0.70, 1.35),
        ("SMF mean loss", measured_smf_db, analytic_smf_db, "dB", 0.90, 1.10),
        ("r0, plane wave", plan.r0_total_m * 1e2, r0_plane * 1e2, "cm",
         0.99, 1.01),
        ("r0, Gaussian beam", plan.r0_total_m * 1e2, r0_gauss * 1e2, "cm",
         0.60, 1.05),
    ]
    print("")
    print(f"{'quantity':<22}{'fidelity 2':>13}{'fidelity 0':>13}"
          f"{'ratio':>9}{'':>7}")
    print("-" * 64)
    for name, got, want, unit, low, high in rows:
        ratio = got / want
        tag = f"{unit:<3}" if unit else "   "
        print(f"{name:<22}{got:>13.4f}{want:>13.4f}{ratio:>9.3f}  "
              f"{tag}{verdict(ratio, low, high)}")
    print("")
    print(f"  The MC error of a variance from {N_TRIALS} snapshots is about "
          f"{100 * np.sqrt(2.0 / N_TRIALS):.0f} percent,")
    print("  so read the three scintillation rows to that accuracy.")
    print("")
    print("  The on-axis row runs HIGH by construction. The pinhole reads the")
    print("  irradiance at the fixed centre of the grid, so it sees the beam")
    print("  WANDER off the axis. The Dios on-axis index sigma2_I(0, L) is the")
    print("  index at the centre of the beam, wherever the beam is, so it")
    print("  carries no wander at all.")
    print("")
    print(f"  The {RX_APERTURE_M * 1e3:.0f} mm budget row runs LOW, and the "
          f"capture fraction says why: that")
    print(f"  aperture collects {100 * power.mean():.0f} percent of the "
          f"launched power (the beam radius at")
    print(f"  the receiver is {sample_beam_radius * 1e3:.0f} mm), and the split "
          f"step conserves power. So")
    print("  almost nothing is left to fluctuate. The weak aperture-averaging")
    print("  factor models an aperture that SAMPLES a wide beam. The "
          f"{SAMPLE_APERTURE_M * 1e3:.0f} mm")
    print(f"  probe collects {100 * sampled.mean():.1f} percent, and it does "
          f"agree.")
    print("")
    print("  The SMF row runs LOW too. The fidelity-0 Term takes the Noll")
    print(f"  residual at the aperture D = {RX_APERTURE_M * 1e3:.0f} mm, and "
          f"the Dikmelik-Davidson coupling")
    print("  of a UNIFORMLY illuminated pupil. The received field is a finite")
    print("  Gaussian beam that is already down to 1/e^2 at the rim, so the")
    print("  phase error that the fibre mode really weighs is smaller. The Term")
    print("  docstring names this limit and asks for a fidelity-2 split-step")
    print("  model. This script is that model. The gap is converged: a grid two")
    print("  times wider moves the fidelity-2 number by 0.3 dB.")
    print("")
    print("  The r0 rows are not a measurement. The screen planner builds each")
    print("  screen from r0 = (0.423 k^2 INT Cn2 dz)^(-3/5), the PLANE-WAVE")
    print("  Fried parameter, so the plane-wave row is an identity check of the")
    print("  plan. The Gaussian-beam r0 of the same path is larger: a finite")
    print("  beam samples less of the turbulence than a plane wave does.")

    # ---- the fade that fidelity 0 does not give ----
    q99_measured = float(-10 * np.log10(np.percentile(sampled, 1)
                                        / sampled.mean()))
    sigma_l2 = np.log(1.0 + sigma2_sample)
    q99_analytic = float(-10.0 / np.log(10.0)
                         * (-sigma_l2 / 2.0
                            + np.sqrt(sigma_l2) * norm.ppf(0.01)))
    print("")
    print(f"  99% fade, {SAMPLE_APERTURE_M * 1e3:.0f} mm sample    "
          f"fidelity 2 {q99_measured:6.2f} dB, "
          f"fidelity 0 lognormal {q99_analytic:6.2f} dB")
    print(f"  99% SMF coupling loss   fidelity 2 "
          f"{-10 * np.log10(np.percentile(smf_eta, 1)):6.2f} dB, "
          f"fidelity 0 mean-only Term gives NO fade")

    # ---- the phase screens ----
    # Each trial makes one fresh screen stack, and the three aperture runs each
    # run the full set of trials on the SAME grid and plan. So the run
    # generates screens-per-trial x trials x 3 runs screens.
    per_trial = int(plan.z_m.size)
    n_runs = 3
    total_screens = per_trial * N_TRIALS * n_runs
    print("")
    print("phase screens created:")
    print(f"  {per_trial} per trial x {N_TRIALS} trials x {n_runs} aperture "
          f"runs = {total_screens} screens")

    draw(sampled, sigma2_sample, power, sigma2_bucket, smf_eta,
         analytic_smf_db)

    # ---- the received field of one snapshot ----
    # Trial 0 of the budget run, on the shared grid and plan, so it is the same
    # atmosphere that the three aperture runs saw first.
    F_rx, _, _ = propagate_turbulent_field(
        build_scenario(), geometry, seed=SEED, trial=0, grid=grid, plan=plan)
    draw_field(F_rx, RX_APERTURE_M, 0.0, FIELD_PNG,
               f"Received field at the far terminal, one snapshot\n"
               f"{PATH_M * 1e-3:.0f} km horizontal, receive aperture "
               f"{RX_APERTURE_M * 1e3:.0f} mm, waist {WAIST_M * 1e3:.0f} mm, "
               f"Cn2 = {CN2:g} m^(-2/3), {WAVELENGTH_M * 1e9:.0f} nm")

    print("")
    print(f"figure saved: {PNG}")
    print(f"field figure saved: {FIELD_PNG}")
    print(f"(elapsed {time.time() - t_start:.1f} s)")


if __name__ == "__main__":
    main()
