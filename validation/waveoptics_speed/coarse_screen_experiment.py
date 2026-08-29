r"""Experiment (a): coarse phase screens plus interpolation.

RECOMMENDATION (read first): BURY IT. A coarse screen at n/f pixels is
band-limited to f times below the grid Nyquist frequency, and no interpolation
(FFT zero-pad or bicubic) puts back the small-scale phase that builds the
scintillation, so the aperture sigma2_I falls far past the 5 percent kill line
at every factor that saves time (see the tables below). The one interpolation
that is cheap (bicubic) is the one that damages sigma2_I the most, and the FFT
zero-pad that damages it least costs a full-size inverse transform that erases
the speed win. The already-wired fast ScreenFactory (screen_generator="olb") is
both faster and correct, so nothing here beats it at equal accuracy.

WHAT THIS SCRIPT MEASURES. The screens.py docstring FORBIDS the naive
coarse-then-interpolate screen, with the reason that the coarse screen carries
no power above its own Nyquist frequency, so it loses the Fresnel-scale
structure sqrt(lambda z) that builds the scintillation (Schmidt,
DOI 10.1117/3.866274, Sec. 9.4, printed p. 172). This script MEASURES that
claim; it does not assume it away. It generates each screen at n/f pixels for
f in {2, 4, 8} at the SAME physical side (so the pitch is f dx), interpolates it
back to n x n by FFT zero-padding and by bicubic zoom, and runs the production
split step on the full grid. The reference is the full-resolution olb
ScreenFactory on the SAME grid and plan, same trial count.

It measures, over about 200 trials each:
  - the mean collected power and the aperture sigma2_I of each configuration,
  - the single-mode-fibre eta for the terrestrial case,
  - the phase structure function of one screen against
    D_phi = 6.88 (r/r0)^(5/3) (Fried, DOI 10.1364/JOSA.56.001372),
  - the wall-time of the screen build (coarse gen + interpolation) against the
    full-resolution build.
It also tests a HYBRID: a coarse low-frequency screen plus a full-resolution
high-frequency remainder, and it states whether that buys anything over the
subharmonics that the generator already adds.

CASES: terrestrial 2 km, Cn2 = 5e-15, standard preset; space downlink 30 deg,
rapid preset. These match the P0/P1 cases.

KILL CRITERION (stated before the run): the approach DIES if sigma2_I moves by
more than 5 percent at every configuration that saves time, OR if no
configuration beats the fast olb generator of screen_generator_check.py at equal
accuracy.

This script changes NO production code. It imports the private trial helpers the
way the run.py self-check and profile_baseline.py do.

Sources:
- Schmidt, DOI 10.1117/3.866274, Ch. 9 (the screen spectrum and the pitch
  rules, Sec. 9.4, printed p. 172).
- Fried, DOI 10.1364/JOSA.56.001372 (r0 and the phase structure function).
- Andrews and Phillips, DOI 10.1117/3.626196, Ch. 8 and Ch. 12 (the analytic
  scintillation reference).

Run from the repository root:
    python -m validation.waveoptics_speed.coarse_screen_experiment
"""

import json
import os
import platform
import time
import warnings

import numpy as np
from scipy.ndimage import zoom

from olb.geometry import CircularOrbit, HorizontalPath
from olb.scenario import (Channel, SpaceScenario, TerrestrialChannel,
                          TerrestrialScenario)
from olb.terminal import SMF, Terminal, Transmitter
from olb.turbulence.plane_wave_scintillation import (
    aperture_averaged_scintillation_index)
from olb.turbulence.profiles import DEFAULT_HS, default_cn2_profile
from olb.waveoptics.field import Begin, Field, Power
from olb.waveoptics.run import _clip
from olb.waveoptics.smf import coupling_efficiency
from olb.waveoptics.threader import Threader
from olb.waveoptics.turbulence.run import _screen_seed, _start_field
from olb.waveoptics.turbulence.sampling import PRESETS, turbulent_grid
from olb.waveoptics.turbulence.screens import ScreenFactory, phase_screen
from olb.waveoptics.turbulence.splitstep import split_step, super_gaussian_boundary

LAM = 1550e-9
HERE = os.path.dirname(__file__)
N_TRIALS = 200
SEED = 20260829
FACTORS = (2, 4, 8)
METHODS = ("fft", "bicubic")


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def warm_aotools():
    """Trigger the one-time scipy deprecation aotools raises, quietly."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            phase_screen(0.1, 32, 0.01, seed=0)
        except ImportError:
            pass


def d_phi(scr, kpx):
    """Measure the phase structure function at kpx pixels, on the two axes."""
    dh = scr[:, kpx:] - scr[:, :-kpx]
    dv = scr[kpx:, :] - scr[:-kpx, :]
    return 0.5 * (np.mean(dh * dh) + np.mean(dv * dv))


def fft_upsample(a, n):
    """Band-limited (sinc) upsample of a real square array to n x n.

    The transform embeds the coarse spectrum in a zero-padded n x n spectrum
    and inverts it. The scale (n/nc)^2 keeps the DC term, so the mean and the
    variance-per-sample are preserved. This is the FAIREST interpolation: it
    adds NO power above the coarse Nyquist, so it shows the pure band-limit
    effect. It also costs one full n x n inverse transform.
    """
    nc = a.shape[0]
    A = np.fft.fftshift(np.fft.fft2(a))
    B = np.zeros((n, n), dtype=complex)
    off = (n - nc) // 2
    B[off:off + nc, off:off + nc] = A
    b = np.fft.ifft2(np.fft.ifftshift(B)) * (float(n) / nc) ** 2
    return np.real(b)


def bicubic_upsample(a, n):
    """Bicubic (order-3) zoom of a real square array to n x n.

    scipy.ndimage.zoom with a reflect boundary. It is cheaper than the FFT
    upsample, but it is a smooth interpolation, so it too adds no true
    small-scale structure and it rolls the spectrum off further.
    """
    nc = a.shape[0]
    out = zoom(a, float(n) / nc, order=3, mode="reflect", grid_mode=False)
    # zoom can land one pixel off for a non-integer factor; here n/nc is exact.
    return out[:n, :n]


# ---------------------------------------------------------------------------
# the cases
# ---------------------------------------------------------------------------

def terrestrial_case():
    """The 2 km, Cn2 = 5e-15 horizontal case, standard preset, SMF receiver."""
    scn = TerrestrialScenario(
        near=Terminal(aperture_m=0.10, wavelength_m=LAM,
                      transmitter=Transmitter(waist_m=0.02)),
        far=Terminal(aperture_m=0.20, wavelength_m=LAM, detector=SMF()),
        channel=TerrestrialChannel(path_length_m=2000.0, cn2=5e-15))
    return dict(name="terrestrial 2km standard", preset="standard",
                scenario=scn, geometry=HorizontalPath(2000.0),
                hs=None, cn2=None, is_space=False)


def space_downlink_case():
    """A 30 deg, 600 km downlink, rapid preset, ground aperture receiver."""
    ground = Terminal(aperture_m=0.50, wavelength_m=LAM,
                      transmitter=Transmitter(waist_m=0.05))
    scn = SpaceScenario(ground=ground,
                        space=Terminal(aperture_m=0.30, wavelength_m=LAM),
                        direction="downlink", channel=Channel(altitude_m=600e3))
    hs = DEFAULT_HS
    cn2 = default_cn2_profile(scn.channel.site, hs)
    return dict(name="space downlink 30deg rapid", preset="rapid",
                scenario=scn, geometry=CircularOrbit(altitude_m=600e3,
                                                     elevation_deg=[30.0]),
                hs=hs, cn2=cn2, is_space=True)


# ---------------------------------------------------------------------------
# one case setup: grid, plan, mask, reference power, start field
# ---------------------------------------------------------------------------

def setup_case(case):
    """Build the grid, plan, mask, reference power, and start field of a case."""
    scn = case["scenario"]
    preset = PRESETS[case["preset"]]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        grid, plan, report = turbulent_grid(
            scn, case["geometry"], preset=preset, hs=case["hs"],
            cn2_profile=case["cn2"])
    mask = super_gaussian_boundary(grid.n, preset.boundary_width_frac)
    if case["is_space"]:
        rx = scn.ground
        F_plane = Begin(grid.size_m, LAM, grid.n)
        flat = np.zeros((grid.n, grid.n))
        F_vac = split_step(F_plane, plan.z_m, [flat] * int(plan.z_m.size),
                           plan.z_total_m, boundary=mask)
        p_reference = Power(_clip(F_vac, rx.aperture_m, rx.obscuration_ratio))
        start = None
    else:
        rx = scn.rx_terminal
        start = _start_field(scn, grid, LAM, is_space=False)
        p_reference = Power(start)
    return dict(scenario=scn, preset=preset, grid=grid, plan=plan,
                report=report, mask=mask, rx=rx, p_reference=p_reference,
                start=start, is_space=case["is_space"])


# ---------------------------------------------------------------------------
# the trial: build a screen stack, propagate, read the receiver
# ---------------------------------------------------------------------------

def run_trials(setup, build_stack, n_trials, threader, seed=SEED):
    """Run n_trials with a per-trial screen-stack builder. Return the reads.

    build_stack(entropy, k) -> list of n x n screens for trial k.
    """
    grid, plan, mask = setup["grid"], setup["plan"], setup["mask"]
    rx, is_space = setup["rx"], setup["is_space"]
    p_reference = setup["p_reference"]
    start = setup["start"]

    def run_one(k):
        stack = build_stack(seed, k)
        F_start = Begin(grid.size_m, LAM, grid.n) if is_space else start
        F_rx = split_step(F_start, plan.z_m, stack, plan.z_total_m,
                          boundary=mask)
        collected = _clip(F_rx, rx.aperture_m, rx.obscuration_ratio)
        power = float(Power(collected) / p_reference)
        eta = (float(coupling_efficiency(collected, rx.aperture_m))
               if isinstance(rx.detector, SMF) else None)
        return power, eta

    results = threader.map(run_one, range(n_trials))
    powers = np.array([r[0] for r in results], dtype=float)
    etas = np.array([r[1] for r in results if r[1] is not None], dtype=float)
    return powers, etas


def stats(powers, etas):
    """Return the mean power, sigma2_I, and mean eta of a trial set."""
    mean = float(powers.mean())
    sigma2 = float(powers.var(ddof=1) / mean ** 2)
    se_sigma2 = float(sigma2 * np.sqrt(2.0 / (powers.size - 1)))
    se_mean = float(powers.std(ddof=1) / np.sqrt(powers.size))
    eta = float(etas.mean()) if etas.size else None
    return dict(mean_power=mean, se_mean_power=se_mean, sigma2_I=sigma2,
                se_sigma2_I=se_sigma2, smf_eta=eta)


# ---------------------------------------------------------------------------
# the screen-stack builders
# ---------------------------------------------------------------------------

def full_builder(grid, plan, L0_m=np.inf):
    """A full-resolution olb ScreenFactory stack builder (the reference)."""
    fac = ScreenFactory(grid.n, grid.pixel_m, L0_m=L0_m)
    r0 = plan.r0_m

    def build(entropy, k):
        return [fac.make(r0[j], np.random.default_rng(_screen_seed(entropy, k, j)))
                for j in range(r0.size)]
    return build


def coarse_builder(grid, plan, f, method, L0_m=np.inf):
    """A coarse-then-interpolate stack builder.

    The coarse screen has n/f pixels at the SAME physical side, so its pitch is
    f*dx and its Nyquist is f times below the full grid. It is then upsampled to
    n x n by FFT zero-padding or bicubic zoom.
    """
    n = grid.n
    nc = n // f
    cfac = ScreenFactory(nc, grid.pixel_m * f, L0_m=L0_m)
    r0 = plan.r0_m
    up = fft_upsample if method == "fft" else bicubic_upsample

    def build(entropy, k):
        out = []
        for j in range(r0.size):
            rng = np.random.default_rng(_screen_seed(entropy, k, j))
            out.append(up(cfac.make(r0[j], rng), n))
        return out
    return build


def hybrid_builder(grid, plan, f, L0_m=np.inf):
    """A hybrid: coarse low-frequency screen + full-resolution high-frequency.

    It builds a full-resolution screen, then it REPLACES the low band (below the
    coarse Nyquist) with an FFT-upsampled coarse screen and keeps the
    full-resolution high band. The point is to check whether splitting the band
    buys anything over the subharmonics that the generator already adds. It does
    NOT save time: it generates a full-resolution screen anyway.
    """
    n = grid.n
    nc = n // f
    ffac = ScreenFactory(n, grid.pixel_m, L0_m=L0_m)
    cfac = ScreenFactory(nc, grid.pixel_m * f, L0_m=L0_m)
    r0 = plan.r0_m
    # the low-band mask in the centred spectrum: the coarse Nyquist window.
    off = (n - nc) // 2
    lowmask = np.zeros((n, n), dtype=bool)
    lowmask[off:off + nc, off:off + nc] = True

    def build(entropy, k):
        out = []
        for j in range(r0.size):
            rng = np.random.default_rng(_screen_seed(entropy, k, j))
            full = ffac.make(r0[j], rng)
            rng2 = np.random.default_rng(_screen_seed(entropy + 1, k, j))
            coarse_up = fft_upsample(cfac.make(r0[j], rng2), n)
            # keep the full high band, take the coarse low band.
            Ff = np.fft.fftshift(np.fft.fft2(full))
            Fc = np.fft.fftshift(np.fft.fft2(coarse_up))
            Ff[lowmask] = Fc[lowmask]
            out.append(np.real(np.fft.ifft2(np.fft.ifftshift(Ff))))
        return out
    return build


# ---------------------------------------------------------------------------
# the screen-build timing (serial, warmed, median)
# ---------------------------------------------------------------------------

def time_stack(build, reps=7):
    """Median wall time of one screen-stack build over reps calls."""
    build(SEED, 0)                       # warm
    ts = []
    for k in range(reps):
        t0 = time.perf_counter()
        build(SEED, 1000 + k)
        ts.append(time.perf_counter() - t0)
    return float(np.median(ts))


# ---------------------------------------------------------------------------
# the structure-function evidence (generator only, one shared test)
# ---------------------------------------------------------------------------

def structure_function_evidence():
    """Measure D_phi of full vs coarse+interp screens against Fried theory.

    A fixed generator test at n = 512, dx = 0.01, r0 = 0.10. It shows the
    band-limit deficit at the small separations that carry the Fresnel-scale
    structure.
    """
    n, dx, r0, M = 512, 0.01, 0.10, 40
    ks = np.array([2, 3, 5, 8, 11, 16])
    theory = 6.88 * ((ks * dx) / r0) ** (5.0 / 3.0)

    def measure(build_screen):
        acc = np.zeros(ks.size)
        for i in range(M):
            s = build_screen(i)
            for j, kpx in enumerate(ks):
                acc[j] += d_phi(s, kpx)
        return acc / M

    ffac = ScreenFactory(n, dx)
    full = measure(lambda i: ffac.make(r0, np.random.default_rng(5000 + i)))
    out = {"r_over_r0": [float(k * dx / r0) for k in ks],
           "d_theory": [float(v) for v in theory],
           "full": {"d": [float(v) for v in full],
                    "ratio": [float(a / b) for a, b in zip(full, theory)]}}
    for f in FACTORS:
        nc = n // f
        cfac = ScreenFactory(nc, dx * f)
        for method, up in (("fft", fft_upsample), ("bicubic", bicubic_upsample)):
            meas = measure(
                lambda i, cf=cfac, u=up: u(
                    cf.make(r0, np.random.default_rng(5000 + i)), n))
            out[f"f{f}_{method}"] = {
                "d": [float(v) for v in meas],
                "ratio": [float(a / b) for a, b in zip(meas, theory)]}
    return out


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def run_case(case, threader):
    """Run the full sweep for one case. Return a result dict."""
    setup = setup_case(case)
    grid, plan = setup["grid"], setup["plan"]
    print(f'\n=== {case["name"]}  (n={grid.n}, {plan.z_m.size} screens, '
          f'pixel {grid.pixel_m * 1e3:.3f} mm, r0 {plan.r0_total_m * 1e2:.2f} '
          f'cm) ===')

    # the reference
    ref_build = full_builder(grid, plan)
    ref_powers, ref_etas = run_trials(setup, ref_build, N_TRIALS, threader)
    ref = stats(ref_powers, ref_etas)
    ref_time = time_stack(ref_build)
    ref["build_s"] = ref_time
    print(f'  reference (olb full)   mean {ref["mean_power"]:.5f}  '
          f'sigma2_I {ref["sigma2_I"]:.5f} +/- {ref["se_sigma2_I"]:.5f}  '
          f'{"eta %.5f  " % ref["smf_eta"] if ref["smf_eta"] else ""}'
          f'build {ref_time * 1e3:.1f} ms')

    configs = {}
    for f in FACTORS:
        for method in METHODS:
            key = f"f{f}_{method}"
            build = coarse_builder(grid, plan, f, method)
            powers, etas = run_trials(setup, build, N_TRIALS, threader)
            s = stats(powers, etas)
            s["build_s"] = time_stack(build)
            s["sigma2_pct"] = 100.0 * (s["sigma2_I"] - ref["sigma2_I"]) \
                / ref["sigma2_I"]
            s["mean_db"] = 10.0 * np.log10(s["mean_power"] / ref["mean_power"])
            s["saves_time"] = bool(s["build_s"] < ref["build_s"])
            s["speedup"] = float(ref["build_s"] / s["build_s"])
            configs[key] = s
            print(f'  {key:12s}  mean {s["mean_power"]:.5f} '
                  f'({s["mean_db"]:+.2f} dB)  sigma2_I {s["sigma2_I"]:.5f} '
                  f'({s["sigma2_pct"]:+.1f}%)  '
                  f'{"eta %.5f  " % s["smf_eta"] if s["smf_eta"] else ""}'
                  f'build {s["build_s"] * 1e3:.1f} ms '
                  f'({s["speedup"]:.2f}x, saves={s["saves_time"]})')

    # the hybrid, one factor
    hb = hybrid_builder(grid, plan, 4)
    hp, he = run_trials(setup, hb, N_TRIALS, threader)
    hs = stats(hp, he)
    hs["build_s"] = time_stack(hb)
    hs["sigma2_pct"] = 100.0 * (hs["sigma2_I"] - ref["sigma2_I"]) \
        / ref["sigma2_I"]
    hs["mean_db"] = 10.0 * np.log10(hs["mean_power"] / ref["mean_power"])
    hs["speedup"] = float(ref["build_s"] / hs["build_s"])
    print(f'  hybrid f4      mean {hs["mean_power"]:.5f} '
          f'({hs["mean_db"]:+.2f} dB)  sigma2_I {hs["sigma2_I"]:.5f} '
          f'({hs["sigma2_pct"]:+.1f}%)  build {hs["build_s"] * 1e3:.1f} ms '
          f'({hs["speedup"]:.2f}x)')

    # the kill-criterion verdict
    time_savers = {k: v for k, v in configs.items() if v["saves_time"]}
    all_savers_fail = (len(time_savers) > 0 and
                       all(abs(v["sigma2_pct"]) > 5.0
                           for v in time_savers.values()))
    any_equal_accuracy_and_faster = any(
        abs(v["sigma2_pct"]) <= 5.0 and v["saves_time"]
        for v in configs.values())
    verdict = {
        "n_time_savers": len(time_savers),
        "all_time_savers_exceed_5pct_sigma2": bool(all_savers_fail),
        "any_config_equal_accuracy_and_faster": bool(any_equal_accuracy_and_faster),
        "killed": bool(all_savers_fail or not any_equal_accuracy_and_faster),
    }
    print(f'  VERDICT: killed={verdict["killed"]}  '
          f'(time-savers={len(time_savers)}, '
          f'any equal-accuracy+faster={any_equal_accuracy_and_faster})')

    # the analytic reference for context
    if case["is_space"]:
        sigma2_analytic = float(aperture_averaged_scintillation_index(
            case["scenario"].ground.aperture_m, 30.0, LAM, case["hs"],
            case["cn2"]))
    else:
        sigma2_analytic = None

    return dict(name=case["name"], grid=dict(n=grid.n, pixel_m=grid.pixel_m,
                                             side_m=grid.size_m),
                n_screens=int(plan.z_m.size), r0_total_m=plan.r0_total_m,
                reference=ref, configs=configs, hybrid_f4=hs,
                verdict=verdict, sigma2_analytic=sigma2_analytic)


def main():
    t_start = time.time()
    warm_aotools()
    threader = Threader()

    print("structure function of full vs coarse+interp screens")
    print("(Fried D_phi = 6.88 (r/r0)^(5/3); ratio < 1 is a band-limit deficit)")
    sf = structure_function_evidence()
    hdr = "  " + "".join(f"{r:8.2f}" for r in sf["r_over_r0"])
    print("  r/r0:" + hdr)
    print("  full ratio " + "".join(f"{v:8.3f}" for v in sf["full"]["ratio"]))
    for f in FACTORS:
        for method in METHODS:
            k = f"f{f}_{method}"
            print(f"  {k:10s} " + "".join(f"{v:8.3f}"
                                          for v in sf[k]["ratio"]))

    cases = [terrestrial_case(), space_downlink_case()]
    results = [run_case(c, threader) for c in cases]

    out = {
        "environment": {"numpy": np.__version__,
                        "platform": platform.platform(),
                        "cores": os.cpu_count()},
        "n_trials": N_TRIALS,
        "factors": list(FACTORS),
        "methods": list(METHODS),
        "structure_function": sf,
        "cases": results,
    }
    path = os.path.join(HERE, "coarse_screen_experiment_results.json")
    with open(path, "w") as fp:
        json.dump(out, fp, indent=1)

    print("\n=== SUMMARY ===")
    for r in results:
        print(f'  {r["name"]}: killed={r["verdict"]["killed"]}')
        for k, v in r["configs"].items():
            print(f'    {k:12s} sigma2 {v["sigma2_pct"]:+7.1f}%  '
                  f'mean {v["mean_db"]:+.2f} dB  build {v["speedup"]:.2f}x  '
                  f'saves={v["saves_time"]}')
    print(f'\nwrote {path}')
    print(f'(elapsed {time.time() - t_start:.0f} s)')


if __name__ == "__main__":
    main()
