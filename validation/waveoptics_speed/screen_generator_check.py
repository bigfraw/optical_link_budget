'''
Check the fast olb ScreenFactory against the aotools phase_screen baseline.

THE POINT. P0 (profile_baseline.py) shows that screen generation is 80 to 84%
of a fidelity-2 trial, and the subharmonic addition is about 70% of a trial.
screens.py now holds a fast, cached generator, ScreenFactory, beside the
aotools-backed phase_screen. This script is the evidence for the three
acceptance items:

  1. STRUCTURE FUNCTION. The new generator passes the structure-function bounds
     of the screens.py self-check: D_phi/theory inside [0.85, 1.02] over
     r/r0 = 0.3 to 1.6, the subharmonics are necessary, and two screens add as
     r0^(-5/3). The tight tests live in the module self-check; this script
     repeats the D_phi ratios as printed evidence.

  2. STATISTICAL EQUIVALENCE. 200 trials of the space downlink at 30 deg, rapid
     preset, aotools vs olb. The mean collected power and the aperture sigma2_I
     agree inside the Monte-Carlo error bars. The two generators draw DIFFERENT
     atmospheres for the same seed, so this is an ensemble comparison, not a
     screen-by-screen one.

  3. SPEED. Seconds per screen, aotools vs olb, for n = 512 to 4096, and
     seconds per full stack for the P0 cases. The batch make_stack uses the
     two-screens-per-transform pairing.

The script also measures the float32 (complex64) error against float64.

Every screen equation cites Schmidt (2010), DOI 10.1117/3.866274, Ch. 9, and
Lane, Glindemann and Dainty, DOI 10.1088/0959-7174/2/3/003 (the subharmonics).
The r0 law is Fried, DOI 10.1364/JOSA.56.001372.

Run from the repository root:
    python -m validation.waveoptics_speed.screen_generator_check
'''

import json
import os
import platform
import time
import warnings

import numpy as np

from olb.geometry import CircularOrbit
from olb.scenario import Channel, SpaceScenario
from olb.terminal import Terminal, Transmitter
from olb.turbulence.plane_wave_scintillation import (
    aperture_averaged_scintillation_index)
from olb.turbulence.profiles import DEFAULT_HS, default_cn2_profile
from olb.waveoptics.threader import Threader
from olb.waveoptics.turbulence.run import propagate_turbulent_scenario
from olb.waveoptics.turbulence.screens import ScreenFactory, phase_screen

LAM = 1550e-9
HERE = os.path.dirname(__file__)


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def d_phi(scr, kpx):
    '''Measure the phase structure function at kpx pixels, on the two axes.'''
    dh = scr[:, kpx:] - scr[:, :-kpx]
    dv = scr[kpx:, :] - scr[:-kpx, :]
    return 0.5 * (np.mean(dh * dh) + np.mean(dv * dv))


def warm_aotools():
    '''Trigger the one-time scipy deprecation that aotools raises, quietly.

    aotools 1.0.7 reads the deprecated scipy.ndimage.interpolation namespace on
    its first call. The deprecation fires ONCE per process. Trigger it here so
    it does not land in a later timed or recorded block.
    '''
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        phase_screen(0.1, 64, 0.01, seed=0)


# ---------------------------------------------------------------------------
# 1. the structure-function evidence
# ---------------------------------------------------------------------------

def structure_function_evidence():
    '''Return the D_phi ratios of the olb generator, and the two side checks.'''
    n, dx, r0, M = 512, 0.01, 0.10, 40
    ks = np.array([3, 5, 8, 11, 16])          # r/r0 = 0.3 to 1.6.
    fac = ScreenFactory(n, dx)
    acc = np.zeros(ks.size)
    for i in range(M):
        s = fac.make(r0, np.random.default_rng(5000 + i))
        for j, kpx in enumerate(ks):
            acc[j] += d_phi(s, kpx)
    acc /= M
    theory = 6.88 * ((ks * dx) / r0) ** (5.0 / 3.0)
    ratio = acc / theory

    # the subharmonics are necessary.
    n2, dx2, r0_2, M2, kbig = 256, 0.01, 0.10, 30, 32
    fac_on = ScreenFactory(n2, dx2, subharmonics=True)
    fac_off = ScreenFactory(n2, dx2, subharmonics=False)
    d_on = d_off = 0.0
    for i in range(M2):
        d_on += d_phi(fac_on.make(r0_2, np.random.default_rng(6000 + i)), kbig)
        d_off += d_phi(fac_off.make(r0_2, np.random.default_rng(6000 + i)),
                       kbig)
    d_on /= M2
    d_off /= M2

    # two screens add as r0^(-5/3), read back at r/r0 = 0.8.
    def fit_r0(scr_list, kpx):
        d = np.mean([d_phi(s, kpx) for s in scr_list])
        return (kpx * dx2) * (6.88 / d) ** 0.6

    r0a, r0b, kfit, M3 = 0.10, 0.16, 8, 30
    sa, sb, ssum = [], [], []
    for i in range(M3):
        a1 = fac_on.make(r0a, np.random.default_rng(7000 + i))
        b1 = fac_on.make(r0b, np.random.default_rng(8000 + i))
        sa.append(a1)
        sb.append(b1)
        ssum.append(a1 + b1)
    fa, fb, fab = fit_r0(sa, kfit), fit_r0(sb, kfit), fit_r0(ssum, kfit)
    add_law = fab / (fa ** (-5 / 3) + fb ** (-5 / 3)) ** -0.6

    return {
        'r_over_r0': [float(k * dx / r0) for k in ks],
        'd_measured': [float(v) for v in acc],
        'd_theory': [float(v) for v in theory],
        'ratio': [float(v) for v in ratio],
        'ratio_in_band': bool(np.all(ratio > 0.85) and np.all(ratio < 1.02)),
        'subharmonic_on': float(d_on),
        'subharmonic_off': float(d_off),
        'subharmonic_deficit': float(1 - d_off / d_on),
        'two_screen_add_law': float(add_law),
    }


# ---------------------------------------------------------------------------
# 2. the statistical-equivalence evidence
# ---------------------------------------------------------------------------

def _stat(power):
    '''Return (mean, se_mean, sigma2_I, se_sigma2) for a power sample.'''
    p = np.asarray(power, dtype=float)
    n = p.size
    mean = float(p.mean())
    se_mean = float(p.std(ddof=1) / np.sqrt(n))
    var = float(p.var(ddof=1))
    sigma2 = var / mean ** 2
    # The standard error of a variance estimate of a near-normal sample is
    # var * sqrt(2/(n-1)). The mean error feeds through, but it is small here.
    se_sigma2 = float(sigma2 * np.sqrt(2.0 / (n - 1)))
    return mean, se_mean, sigma2, se_sigma2


def statistical_equivalence(n_trials=200, seed=2024):
    '''Compare the two generators on 200 downlink trials at 30 deg, rapid.'''
    hs = DEFAULT_HS
    orbit = CircularOrbit(altitude_m=600e3, elevation_deg=[30.0])
    ground = Terminal(aperture_m=0.40, wavelength_m=LAM,
                      transmitter=Transmitter(waist_m=0.06))
    down = SpaceScenario(ground=ground,
                         space=Terminal(aperture_m=0.30, wavelength_m=LAM),
                         direction='downlink', channel=Channel())
    threader = Threader()

    out = {}
    powers = {}
    for gen in ('aotools', 'olb'):
        t0 = time.perf_counter()
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            res = propagate_turbulent_scenario(
                down, orbit, n_trials=n_trials, seed=seed, preset='rapid',
                threader=threader, screen_generator=gen)
        wall = time.perf_counter() - t0
        power = np.array([tr.collected_power for tr in res.trials])
        powers[gen] = power
        mean, se_mean, sigma2, se_sigma2 = _stat(power)
        out[gen] = {'mean_power': mean, 'se_mean_power': se_mean,
                    'sigma2_I': sigma2, 'se_sigma2_I': se_sigma2,
                    'wall_s': float(wall)}

    # The analytic reference, for context.
    cn2 = default_cn2_profile(down.channel.site, hs)
    sigma2_theory = float(aperture_averaged_scintillation_index(
        ground.aperture_m, 30.0, LAM, hs, cn2))

    d_mean = out['aotools']['mean_power'] - out['olb']['mean_power']
    se_mean = np.hypot(out['aotools']['se_mean_power'],
                       out['olb']['se_mean_power'])
    d_sig = out['aotools']['sigma2_I'] - out['olb']['sigma2_I']
    se_sig = np.hypot(out['aotools']['se_sigma2_I'], out['olb']['se_sigma2_I'])
    out['comparison'] = {
        'n_trials': n_trials,
        'sigma2_I_analytic': sigma2_theory,
        'delta_mean_power': float(d_mean),
        'se_delta_mean_power': float(se_mean),
        'mean_power_sigmas_apart': float(abs(d_mean) / se_mean),
        'delta_sigma2_I': float(d_sig),
        'se_delta_sigma2_I': float(se_sig),
        'sigma2_I_sigmas_apart': float(abs(d_sig) / se_sig),
        'mean_power_agrees_2sigma': bool(abs(d_mean) < 2.0 * se_mean),
        'sigma2_I_agrees_2sigma': bool(abs(d_sig) < 2.0 * se_sig),
    }
    return out


# ---------------------------------------------------------------------------
# 3. the speed table
# ---------------------------------------------------------------------------

def _time_call(fn, reps):
    '''Return the mean wall time of fn over reps calls, after one warm call.'''
    fn()
    t0 = time.perf_counter()
    for _ in range(reps):
        fn()
    return (time.perf_counter() - t0) / reps


def speed_per_screen(ns=(512, 1024, 2048, 4096), dx=0.01, r0=0.10):
    '''Seconds per screen, aotools vs olb make, over a range of grid sizes.'''
    rows = []
    for n in ns:
        reps = max(3, int(round(6 * (512 / n) ** 2)))
        # aotools: one full spectrum build for each screen.
        rng_seed = [0]

        def aot():
            rng_seed[0] += 1
            return phase_screen(r0, n, dx, seed=rng_seed[0])

        # olb: build the factory once, then time make only.
        t_build = time.perf_counter()
        fac = ScreenFactory(n, dx)
        build_s = time.perf_counter() - t_build
        rng = np.random.default_rng(1)

        def olb_make():
            return fac.make(r0, rng)

        # olb make_stack: the pairing, amortised per screen over a stack.
        rng2 = np.random.default_rng(2)

        def olb_stack():
            return fac.make_stack([r0] * 8, rng2)

        t_aot = _time_call(aot, reps)
        t_make = _time_call(olb_make, reps)
        t_stack = _time_call(olb_stack, max(3, reps // 2)) / 8.0
        rows.append({
            'n': int(n),
            'aotools_s': float(t_aot),
            'olb_make_s': float(t_make),
            'olb_make_stack_s_per_screen': float(t_stack),
            'olb_factory_build_s': float(build_s),
            'speedup_make': float(t_aot / t_make),
            'speedup_stack': float(t_aot / t_stack),
        })
    return rows


def speed_per_stack(cases, dx=0.01, r0=0.10):
    '''Seconds for a whole screen stack, for the P0 cases.

    The timing is independent of the r0 values, so one r0 stands for the whole
    stack. The grid n and the screen count come from the P0 cases.
    '''
    rows = []
    for c in cases:
        n = int(c['grid']['n'])
        m = int(c['plan']['n_screens'])
        reps = max(3, int(round(4 * (256 / n) ** 2)))

        def aot():
            return [phase_screen(r0, n, dx, seed=1000 + j) for j in range(m)]

        fac = ScreenFactory(n, dx)

        def olb_loop():
            rng = np.random.default_rng(3)
            return [fac.make(r0, rng) for _ in range(m)]

        def olb_stack():
            return fac.make_stack([r0] * m, np.random.default_rng(4))

        rows.append({
            'name': c['name'],
            'n': n,
            'n_screens': m,
            'aotools_stack_s': float(_time_call(aot, reps)),
            'olb_make_loop_s': float(_time_call(olb_loop, reps)),
            'olb_make_stack_s': float(_time_call(olb_stack, reps)),
        })
    for r in rows:
        r['speedup_stack'] = float(r['aotools_stack_s'] / r['olb_make_stack_s'])
    return rows


# ---------------------------------------------------------------------------
# 4. the float32 error
# ---------------------------------------------------------------------------

def float32_error(n=512, dx=0.01, r0=0.10, M=20):
    '''Measure the complex64/float32 screen error against float64.'''
    fac64 = ScreenFactory(n, dx, dtype=np.float64)
    fac32 = ScreenFactory(n, dx, dtype=np.float32)
    rels = []
    for i in range(M):
        s64 = fac64.make(r0, np.random.default_rng(200 + i))
        s32 = fac32.make(r0, np.random.default_rng(200 + i))
        rels.append(float(np.sqrt(np.mean((s32 - s64) ** 2)) / np.std(s64)))
    return {'rel_rms_mean': float(np.mean(rels)),
            'rel_rms_max': float(np.max(rels)),
            'note': ('numpy.fft upcasts to complex128 internally, so float32 '
                     'saves memory, not FFT time; the error is the float32 '
                     'storage of the filter and the output only')}


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main():
    t_start = time.time()
    warm_aotools()

    # Load the P0 cases, to keep the per-stack table in sync with the baseline.
    p0_path = os.path.join(HERE, 'profile_baseline_results.json')
    with open(p0_path) as f:
        p0 = json.load(f)
    p0_cases = p0['cases']

    print('1. structure function of the olb generator')
    sf = structure_function_evidence()
    for r, dm, th, rt in zip(sf['r_over_r0'], sf['d_measured'],
                             sf['d_theory'], sf['ratio']):
        print(f'   r/r0 {r:4.2f}  D {dm:8.3f}  theory {th:8.3f}  '
              f'ratio {rt:6.3f}')
    print(f'   ratio in [0.85, 1.02]: {sf["ratio_in_band"]}')
    print(f'   subharmonic on {sf["subharmonic_on"]:.3f}, '
          f'off {sf["subharmonic_off"]:.3f}, '
          f'deficit {sf["subharmonic_deficit"]:.3f}')
    print(f'   two-screen add law (target 1.0): {sf["two_screen_add_law"]:.4f}')
    print('')

    print('2. statistical equivalence, 200 downlink trials at 30 deg, rapid')
    se = statistical_equivalence()
    for gen in ('aotools', 'olb'):
        g = se[gen]
        print(f'   {gen:8s}  mean power {g["mean_power"]:.5f} '
              f'+/- {g["se_mean_power"]:.5f}   '
              f'sigma2_I {g["sigma2_I"]:.5f} +/- {g["se_sigma2_I"]:.5f}   '
              f'({g["wall_s"]:.1f} s)')
    cmp = se['comparison']
    print(f'   analytic sigma2_I {cmp["sigma2_I_analytic"]:.5f}')
    print(f'   delta mean power  {cmp["delta_mean_power"]:+.5f} '
          f'({cmp["mean_power_sigmas_apart"]:.2f} sigma, '
          f'agrees {cmp["mean_power_agrees_2sigma"]})')
    print(f'   delta sigma2_I    {cmp["delta_sigma2_I"]:+.5f} '
          f'({cmp["sigma2_I_sigmas_apart"]:.2f} sigma, '
          f'agrees {cmp["sigma2_I_agrees_2sigma"]})')
    print('')

    print('3a. seconds per screen, aotools vs olb')
    per_screen = speed_per_screen()
    print('    n     aotools     olb make   stack/scr   x make   x stack')
    for r in per_screen:
        print(f'   {r["n"]:5d}  {r["aotools_s"]:9.4f}  {r["olb_make_s"]:9.4f}  '
              f'{r["olb_make_stack_s_per_screen"]:9.4f}  '
              f'{r["speedup_make"]:6.1f}x  {r["speedup_stack"]:6.1f}x')
    print('')

    print('3b. seconds per full stack, the P0 cases')
    per_stack = speed_per_stack(p0_cases)
    print('    case                            n   scr   aotools   olb loop  '
          'olb stack   x')
    for r in per_stack:
        print(f'   {r["name"]:28s}  {r["n"]:5d} {r["n_screens"]:4d}  '
              f'{r["aotools_stack_s"]:8.4f}  {r["olb_make_loop_s"]:8.4f}  '
              f'{r["olb_make_stack_s"]:8.4f}  {r["speedup_stack"]:5.1f}x')
    print('')

    print('4. float32 (complex64) error against float64')
    f32 = float32_error()
    print(f'   rel rms mean {f32["rel_rms_mean"]:.2e}, '
          f'max {f32["rel_rms_max"]:.2e}')
    print(f'   {f32["note"]}')
    print('')

    out = {
        'environment': {
            'numpy': np.__version__,
            'platform': platform.platform(),
            'cores': os.cpu_count(),
        },
        'structure_function': sf,
        'statistical_equivalence': se,
        'speed_per_screen': per_screen,
        'speed_per_stack': per_stack,
        'float32_error': f32,
    }
    path = os.path.join(HERE, 'screen_generator_check_results.json')
    with open(path, 'w') as f:
        json.dump(out, f, indent=1)
    print(f'wrote {path}')
    print(f'(elapsed {time.time() - t_start:.0f} s)')


if __name__ == '__main__':
    main()
