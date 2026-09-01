r"""The fidelity-2 profiling baseline: where does one turbulent trial spend its
time?

This script splits ONE turbulent split-step trial into its parts, for a set of
representative cases. It gates the rest of the speed plan (see
docs/waveoptics-efficiency-plan.md): the later tasks attack the biggest slice
first. It changes NO production code. It reproduces the run_one body of
olb.waveoptics.turbulence.run with a timer around each part, and it imports the
private helpers the same way the run.py self-check does.

The trial has three parts:

1. Screen generation. For each screen of the stack, the time splits into the
   base Fourier-series draw (aotools ft_phase_screen) and the subharmonic
   addition (the extra time of aotools ft_sh_phase_screen). See
   olb.waveoptics.turbulence.screens.phase_screen.
2. The split step. The time splits into the free-space hops (Forvard, two FFTs
   each) and the screen and mask multiplies. The script derives the hop and
   sub-step count from the plan (z_m, z_total_m, forvard_max_z) and it reports
   the seconds per FFT.
3. The scalar reads: the receive clip, Power, coupling_efficiency (an SMF
   detector on the receive terminal turns this path on), and the uplink overlap
   for the uplink case.

The script uses the seeds that propagate_turbulent_scenario would use for
seed=0, trial 0, so the numbers repeat. It times the whole trial REPS times and
it reports the median of each part, because the wall clock is noisy; the
physics numbers come from the deterministic trial 0.

RESULTS (2026-08-29, this machine; see profile_baseline_run.log):
- SCREEN GENERATION dominates EVERY case, at 80 to 84 percent of the trial. The
  split step (the Forvard hops plus the screen and mask multiplies) is 15 to 19
  percent, and the scalar reads are about 1 percent.
- Inside screen generation, the SUBHARMONIC addition is about 7 times the base
  Fourier draw. So the subharmonics alone are about 70 percent of the whole
  trial. That is the biggest slice, and it is what P1 must attack (the
  27-mode subharmonic sum builds a full N x N complex exponential per mode in a
  Python loop; the separable outer-product win of the plan removes that).
- The share is flat across the cases:
  * terrestrial 2 km standard   (256 px,  9 screens): screen gen 83%
  * space downlink 30 deg rapid (256 px,  5 screens): screen gen 82%
  * space downlink 30 deg std   (512 px,  9 screens): screen gen 84%
  * space downlink 30 deg ref  (1024 px, 15 screens): screen gen 80%
  * space uplink 30 deg std     (512 px,  9 screens): screen gen 83%
- A Forvard hop costs about 3 times the raw fft2 of its size, because it rebuilds
  the N x N transfer function (cos and sin) on every call. That is a smaller,
  secondary target (relevant to the P3 batched-FFT idea).
See the printed SUMMARY block and profile_baseline_results.json for the numbers.

Run from the repository root:
    python -m validation.waveoptics_speed.profile_baseline
"""

import json
import os
import platform
import time
import warnings
from statistics import median

import numpy as np
import scipy.fft as _sfft

from olb.geometry import CircularOrbit, HorizontalPath
from olb.scenario import (Channel, SpaceScenario, TerrestrialChannel,
                          TerrestrialScenario)
from olb.terminal import SMF, Terminal, Transmitter
from olb.waveoptics.field import Begin, Power
from olb.waveoptics.grid import forvard_max_z
from olb.waveoptics.propagators import Forvard
from olb.waveoptics.run import _clip, _launch_aperture, _normalised_gauss
from olb.waveoptics.smf import coupling_efficiency
from olb.waveoptics.sources import GaussBeam
from olb.beam import virtual_waist
# The private helpers of the trial runner, imported the way the run.py
# self-check does.
from olb.waveoptics.turbulence.run import (_ground_transmit_mode,
                                           _resolve_seed, _screen_seed)
from olb.waveoptics.turbulence.sampling import PRESETS, turbulent_grid
from olb.waveoptics.turbulence.screens import phase_screen
from olb.waveoptics.turbulence.splitstep import (_apply_mask, _substeps,
                                                 super_gaussian_boundary)
from olb.waveoptics.turbulence.screens import Screen

LAM = 1550e-9
SEED = 0                 # propagate_turbulent_scenario(seed=0), trial 0
TRIAL = 0
REPS = 3                 # the whole trial runs REPS times; report the median
FFT_SIZES = (512, 1024, 2048, 4096)
FFT_REPS = 5


# --------------------------------------------------------------------------
# the cases
# --------------------------------------------------------------------------
def terrestrial_case():
    """The 2 km, Cn2 = 5e-15 horizontal case, standard preset, SMF receiver.

    The terminals match the sampling.py self-check case 1, plus an SMF detector
    on the receiver so the coupling read runs.
    """
    scn = TerrestrialScenario(
        near=Terminal(aperture_m=0.10, wavelength_m=LAM,
                      transmitter=Transmitter(waist_m=0.02)),
        far=Terminal(aperture_m=0.20, wavelength_m=LAM, detector=SMF()),
        channel=TerrestrialChannel(path_length_m=2000.0, cn2=5e-15))
    return scn, HorizontalPath(2000.0), None, None


def space_case(direction):
    """A 30 deg, 600 km space case, ground with a transmitter AND an SMF.

    The ground carries a Transmitter (the uplink overlap reads its mode) and an
    SMF detector (the coupling read runs). The terminals match the sampling.py
    self-check space case.
    """
    ground = Terminal(aperture_m=0.50, wavelength_m=LAM,
                      transmitter=Transmitter(waist_m=0.05), detector=SMF())
    scn = SpaceScenario(ground=ground,
                        space=Terminal(aperture_m=0.30, wavelength_m=LAM),
                        direction=direction, channel=Channel(altitude_m=600e3))
    geom = CircularOrbit(altitude_m=600e3, elevation_deg=[30.0])
    from olb.turbulence.profiles import DEFAULT_HS, default_cn2_profile
    hs = DEFAULT_HS
    cn2 = default_cn2_profile(scn.channel.site, hs)
    return scn, geom, hs, cn2


CASES = [
    ("terrestrial 2km standard", "standard", terrestrial_case),
    ("space downlink 30deg rapid", "rapid",
     lambda: space_case("downlink")),
    ("space downlink 30deg standard", "standard",
     lambda: space_case("downlink")),
    ("space downlink 30deg reference", "reference",
     lambda: space_case("downlink")),
    ("space uplink 30deg standard", "standard",
     lambda: space_case("uplink")),
]


# --------------------------------------------------------------------------
# the timed trial
# --------------------------------------------------------------------------
def timed_screen_stack(plan, grid, entropy):
    """Make the screen stack the way run_one does, and time it.

    Each screen is ONE aotools ft_sh_phase_screen draw (subharmonics on), the
    exact call of run_one. So this time IS the real screen-generation cost of
    the trial. The base-vs-subharmonic split comes from measure_base_fft, which
    times ft_phase_screen (subharmonics off) separately; the subharmonic part is
    then the difference (ft_sh time - ft time).

    Returns:
        (stack, sh_s) with the total ft_sh time over the stack, in seconds.
    """
    stack, sh_s = [], 0.0
    for j in range(int(plan.z_m.size)):
        s = _screen_seed(entropy, TRIAL, j)
        t0 = time.perf_counter()
        full = phase_screen(plan.r0_m[j], grid.n, grid.pixel_m, L0_m=np.inf,
                            seed=s, subharmonics=True)
        sh_s += time.perf_counter() - t0
        stack.append(full)
    return stack, sh_s


def measure_base_fft(plan, grid, entropy):
    """Time the base Fourier draw (aotools ft_phase_screen, subharmonics off).

    This is a DIAGNOSTIC pass, outside the trial timer. The base time splits the
    ft_sh cost into the base FFT part and the subharmonic addition (the
    difference). It does not go into the trial total, because run_one draws
    ft_sh only.

    Returns:
        The total base-draw time over the stack, in seconds.
    """
    base_s = 0.0
    for j in range(int(plan.z_m.size)):
        s = _screen_seed(entropy, TRIAL, j)
        t0 = time.perf_counter()
        phase_screen(plan.r0_m[j], grid.n, grid.pixel_m, L0_m=np.inf,
                     seed=s, subharmonics=False)
        base_s += time.perf_counter() - t0
    return base_s


def timed_split_step(F_start, plan, stack, mask, max_step_m):
    """Reproduce splitstep.split_step with per-part timers.

    The function mirrors the split_step body: hop, screen, mask, ... , final
    hop. It counts the Forvard calls (each is two FFTs) and it sums the hop
    time, the screen time and the mask time.

    Returns:
        (F_out, hop_s, screen_s, mask_s, n_forvard, n_substeps).
    """
    z = np.asarray(plan.z_m, dtype=float).ravel()
    z_total = float(plan.z_total_m)
    hop_s = screen_s = mask_s = 0.0
    n_forvard = n_substeps = 0

    def hop(F, gap_m):
        nonlocal hop_s, mask_s, n_forvard, n_substeps
        for dz in _substeps(gap_m, max_step_m):
            t0 = time.perf_counter()
            F = Forvard(F, dz)
            hop_s += time.perf_counter() - t0
            n_forvard += 1
            n_substeps += 1
            t0 = time.perf_counter()
            F = _apply_mask(F, mask)
            mask_s += time.perf_counter() - t0
        return F

    from olb.waveoptics.field import Field
    F = Field.copy(F_start)
    here = 0.0
    for zi, scr in zip(z, stack):
        F = hop(F, zi - here)
        t0 = time.perf_counter()
        F = Screen(F, scr)
        screen_s += time.perf_counter() - t0
        t0 = time.perf_counter()
        F = _apply_mask(F, mask)
        mask_s += time.perf_counter() - t0
        here = zi
    F = hop(F, z_total - here)
    return F, hop_s, screen_s, mask_s, n_forvard, n_substeps


def run_case(name, preset_name, builder):
    """Profile one case. Return a result dict."""
    scn, geom, hs, cn2 = builder()
    p = PRESETS[preset_name]
    is_space = hasattr(scn, "ground")

    warns = []
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        grid, plan, report = turbulent_grid(scn, geom, preset=p, hs=hs,
                                            cn2_profile=cn2)
    warns = [str(w.message) for w in caught]

    lam = scn.tx_terminal.wavelength_m
    rx = scn.ground if is_space else scn.rx_terminal
    mask = super_gaussian_boundary(grid.n, p.boundary_width_frac)
    max_step_m = grid.n * grid.pixel_m ** 2 / lam

    # ---- the read-only setup (computed once, timed for interest) ----
    t0 = time.perf_counter()
    psi_tx = o_vac = None
    if is_space:
        F_plane = Begin(grid.size_m, lam, grid.n)
        flat = np.zeros((grid.n, grid.n))
        F_vac, *_ = timed_split_step(F_plane, plan, [flat] * int(plan.z_m.size),
                                     mask, max_step_m)
        p_reference = Power(_clip(F_vac, rx.aperture_m, rx.obscuration_ratio))
        if scn.direction == "uplink":
            psi_tx = _ground_transmit_mode(scn.ground, grid)
            o_vac = float(np.abs((F_vac.field * np.conj(psi_tx)).sum()) ** 2)
        F_in = None
    else:
        tx = scn.tx_terminal
        t = tx.transmitter
        w_v, offset = virtual_waist(t.waist_m, t.divergence_rad, lam)
        F0 = _normalised_gauss(GaussBeam(Begin(grid.size_m, lam, grid.n), w_v))
        if offset > 0:
            from olb.waveoptics.propagators import GForvard
            F0 = GForvard(F0, offset)
        F_in = _clip(F0, *_launch_aperture(tx))
        p_reference = Power(F_in)
    setup_s = time.perf_counter() - t0

    entropy = _resolve_seed(SEED)
    n_screens = int(plan.z_m.size)

    # ---- the base-FFT diagnostic (outside the trial timer) ----
    base_list = [measure_base_fft(plan, grid, entropy) for _ in range(REPS)]

    # ---- the timed trial, REPS times, median of each part ----
    parts = {k: [] for k in ("scr", "hop", "screen", "mask",
                             "clip", "power", "coupling", "overlap", "total")}
    n_forvard = n_substeps = 0
    physics = {}
    for _ in range(REPS):
        t_trial = time.perf_counter()
        stack, scr_s = timed_screen_stack(plan, grid, entropy)
        F_start = (Begin(grid.size_m, lam, grid.n) if is_space else F_in)
        (F_rx, hop_s, screen_s, mask_s,
         n_forvard, n_substeps) = timed_split_step(F_start, plan, stack, mask,
                                                   max_step_m)

        t0 = time.perf_counter()
        collected = _clip(F_rx, rx.aperture_m, rx.obscuration_ratio)
        clip_s = time.perf_counter() - t0
        t0 = time.perf_counter()
        cp = float(Power(collected) / p_reference)
        power_s = time.perf_counter() - t0
        t0 = time.perf_counter()
        smf_eta = (float(coupling_efficiency(collected, rx.aperture_m))
                   if isinstance(rx.detector, SMF) else None)
        coupling_s = time.perf_counter() - t0
        overlap_s = 0.0
        eta_turb = None
        if is_space and scn.direction == "uplink":
            t0 = time.perf_counter()
            o = float(np.abs((F_rx.field * np.conj(psi_tx)).sum()) ** 2)
            eta_turb = o / o_vac
            overlap_s = time.perf_counter() - t0
        total_s = time.perf_counter() - t_trial

        parts["scr"].append(scr_s)
        parts["hop"].append(hop_s)
        parts["screen"].append(screen_s)
        parts["mask"].append(mask_s)
        parts["clip"].append(clip_s)
        parts["power"].append(power_s)
        parts["coupling"].append(coupling_s)
        parts["overlap"].append(overlap_s)
        parts["total"].append(total_s)
        physics = {"collected_power": cp, "smf_eta": smf_eta,
                   "eta_turb": eta_turb}

    med = {k: median(v) for k, v in parts.items()}
    # screen generation IS the ft_sh cost of the trial. The base FFT part is the
    # diagnostic ft draw; the subharmonic part is the difference (clamped at 0
    # against timer noise, because ft_sh always does the base draw too).
    screen_gen_s = med["scr"]
    base_fft_s = median(base_list)
    subharmonic_s = max(0.0, screen_gen_s - base_fft_s)
    propagation_s = med["hop"] + med["screen"] + med["mask"]
    scalar_s = med["clip"] + med["power"] + med["coupling"] + med["overlap"]
    n_ffts = 2 * n_forvard
    s_per_fft = med["hop"] / n_ffts if n_ffts else 0.0
    denom = screen_gen_s + propagation_s + scalar_s

    return {
        "name": name,
        "preset": preset_name,
        "grid": {"n": grid.n, "side_m": grid.size_m, "pixel_m": grid.pixel_m},
        "plan": {"n_screens": n_screens, "z_total_m": plan.z_total_m,
                 "r0_total_m": plan.r0_total_m,
                 "n_substeps": n_substeps, "n_forvard": n_forvard,
                 "n_ffts": n_ffts},
        "report": {"pixels_per_r0": report.pixels_per_r0,
                   "grid_margin": report.grid_margin,
                   "fresnel_pixels_min": report.fresnel_pixels_min,
                   "step_over_limit_max": report.step_over_limit_max,
                   "sigma2_r_screen_max": report.sigma2_r_screen_max,
                   "n_clamped": report.n_clamped},
        "physics": physics,
        "timing_s": {
            "screen_base_fft": base_fft_s, "screen_subharmonic": subharmonic_s,
            "screen_gen_total": screen_gen_s,
            "prop_hops": med["hop"], "prop_screen_apply": med["screen"],
            "prop_mask_apply": med["mask"], "prop_total": propagation_s,
            "scalar_clip": med["clip"], "scalar_power": med["power"],
            "scalar_coupling": med["coupling"], "scalar_overlap": med["overlap"],
            "scalar_total": scalar_s,
            "s_per_fft": s_per_fft, "setup_s": setup_s,
            "trial_total": med["total"]},
        "share_pct": {
            "screen_gen": 100.0 * screen_gen_s / denom,
            "propagation": 100.0 * propagation_s / denom,
            "scalar_reads": 100.0 * scalar_s / denom},
        "warnings": warns,
    }


# --------------------------------------------------------------------------
# the raw FFT microbenchmark
# --------------------------------------------------------------------------
def fft_microbench():
    """Time numpy.fft.fft2 and scipy.fft.fft2(workers=1) on complex128."""
    rng = np.random.default_rng(0)
    out = {}
    for n in FFT_SIZES:
        a = (rng.standard_normal((n, n))
             + 1j * rng.standard_normal((n, n))).astype(np.complex128)
        np.fft.fft2(a)                       # warm up
        tn = min(_time(lambda: np.fft.fft2(a)) for _ in range(FFT_REPS))
        _sfft.fft2(a, workers=1)
        ts = min(_time(lambda: _sfft.fft2(a, workers=1))
                 for _ in range(FFT_REPS))
        out[str(n)] = {"numpy_fft2_s": tn, "scipy_fft2_w1_s": ts}
    return out


def _time(fn):
    t0 = time.perf_counter()
    fn()
    return time.perf_counter() - t0


# --------------------------------------------------------------------------
# the driver
# --------------------------------------------------------------------------
def main():
    t_start = time.time()
    env = {
        "numpy": np.__version__,
        "scipy": __import__("scipy").__version__,
        "aotools": __import__("aotools").__version__,
        "cpu": platform.processor(),
        "machine": platform.machine(),
        "cores": os.cpu_count(),
        "platform": platform.platform(),
        "reps": REPS,
    }
    print("environment:")
    for k, v in env.items():
        print(f"  {k:22s} {v}")
    print("")

    print("raw FFT microbenchmark (complex128, best of "
          f"{FFT_REPS}):")
    fft = fft_microbench()
    print(f"  {'n':>6}{'numpy fft2 [ms]':>18}{'scipy fft2 w1 [ms]':>20}")
    for n in FFT_SIZES:
        row = fft[str(n)]
        print(f"  {n:>6}{row['numpy_fft2_s'] * 1e3:>18.3f}"
              f"{row['scipy_fft2_w1_s'] * 1e3:>20.3f}")
    print("")

    results = []
    for name, preset_name, builder in CASES:
        print(f"profiling: {name} ...", flush=True)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            res = run_case(name, preset_name, builder)
        results.append(res)

    # ---- the per-case table ----
    print("")
    print("per-case breakdown (median of the trial, shares in percent):")
    hdr = (f"  {'case':<34}{'n':>6}{'scr':>4}{'ffts':>6}"
           f"{'trial[s]':>10}{'scr%':>7}{'prop%':>7}{'scal%':>7}")
    print(hdr)
    for r in results:
        print(f"  {r['name']:<34}{r['grid']['n']:>6}"
              f"{r['plan']['n_screens']:>4}{r['plan']['n_ffts']:>6}"
              f"{r['timing_s']['trial_total']:>10.3f}"
              f"{r['share_pct']['screen_gen']:>7.1f}"
              f"{r['share_pct']['propagation']:>7.1f}"
              f"{r['share_pct']['scalar_reads']:>7.1f}")
    print("")

    # ---- screen generation split ----
    print("screen generation split (base FFT vs subharmonics, seconds):")
    print(f"  {'case':<34}{'base[s]':>10}{'subh[s]':>10}{'sub/base':>10}")
    for r in results:
        b = r['timing_s']['screen_base_fft']
        s = r['timing_s']['screen_subharmonic']
        print(f"  {r['name']:<34}{b:>10.4f}{s:>10.4f}"
              f"{(s / b if b else 0):>10.2f}")
    print("")

    # ---- propagation detail ----
    print("propagation detail:")
    print(f"  {'case':<34}{'hops[s]':>10}{'s/fft[ms]':>11}"
          f"{'substeps':>9}{'scr+mask[s]':>13}")
    for r in results:
        t = r['timing_s']
        print(f"  {r['name']:<34}{t['prop_hops']:>10.3f}"
              f"{t['s_per_fft'] * 1e3:>11.3f}"
              f"{r['plan']['n_substeps']:>9}"
              f"{t['prop_screen_apply'] + t['prop_mask_apply']:>13.4f}")
    print("")

    # ---- the SUMMARY block: the biggest slice per case ----
    print("SUMMARY (the biggest slice per case):")
    for r in results:
        sh = r['share_pct']
        biggest = max(sh, key=sh.get)
        label = {"screen_gen": "screen generation",
                 "propagation": "the split step (Forvard hops)",
                 "scalar_reads": "the scalar reads"}[biggest]
        print(f"  {r['name']}: {label} takes {sh[biggest]:.0f}% "
              f"({r['grid']['n']} px, {r['plan']['n_screens']} screens, "
              f"{r['plan']['n_ffts']} FFTs).")
    print("")

    out = {"environment": env, "fft_microbench": fft, "cases": results}
    path = os.path.join(os.path.dirname(__file__),
                        "profile_baseline_results.json")
    with open(path, "w") as f:
        json.dump(out, f, indent=1)
    print(f"wrote {os.path.relpath(path)}")
    print(f"(elapsed {time.time() - t_start:.1f} s)")


if __name__ == '__main__':
    main()
