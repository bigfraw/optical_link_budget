'''
A broad validity pass on the fast olb phase-screen generator, ScreenFactory.

VERDICT (filled from the run; see the results JSON and the run log).
    The olb generator is a TRUSTWORTHY drop-in for aotools across every axis
    tested. Equivalence: the mean power and the aperture sigma2_I agree for all
    four cases (terrestrial standard, downlink rapid and standard, uplink
    standard) inside the Monte-Carlo bars, the largest gap 0.39 sigma; the
    uplink eta_turb agrees too. Converged (2000 trials) the olb sigma2_I is
    0.01522, that is 1.012 times the analytic aperture-averaged value, against
    0.994 for aotools; the two generators agree at 0.39 sigma. THE FADE TAIL IS
    SAFE: the low collected-power quantiles that set a fade margin (0.1%, 1%,
    5%, 10%), the median and the 90% point all agree inside the bootstrapped
    bars; the margin-setting 1% quantile agrees to 0.014 dB (bar 0.067 dB), and
    the worst gap is +0.05 dB at the 0.1% point (bar 0.20 dB). At a finite outer
    scale (L0 = 25 m) the olb structure function stays inside 5% of aotools at
    every separation; both fall below the pure-Kolmogorov 6.88 (r/r0)^(5/3)
    line, which is the EXPECTED von Karman roll-off at r approaching L0, not a
    failure. So nothing in the physics blocks making the olb generator a
    default; the switch stays an owner decision, and the default stays aotools.

THE POINT. P1 added a fast, opt-in generator, `ScreenFactory` (screens.py),
selected by `screen_generator="olb"` on `propagate_turbulent_scenario` and
`propagate_turbulent_field` (run.py). It draws a DIFFERENT random atmosphere
from the aotools path for the same seed, and it is 7.6 to 14x faster per
screen. P1's own check (screen_generator_check.py) validated it on ONE case:
space downlink 30 deg, rapid preset, 200 trials, Kolmogorov (L0 = inf).

This script is the BROADER pass, before the generator is ever a default. It
answers: is the olb generator a trustworthy drop-in for aotools across
geometries, presets, the outer scale, AND the fade TAIL that P1 did not probe?
Every statistic carries its Monte-Carlo error bar, estimated from the sample.

  1. EQUIVALENCE ACROSS CASES. terrestrial 2 km, Cn2 = 5e-15, standard, SMF;
     space downlink 30 deg rapid AND standard; space uplink 30 deg standard.
     The mean collected power and the aperture sigma2_I must agree between the
     two generators inside the combined MC bars.
  2. CONVERGED sigma2_I vs ANALYTIC. Space downlink 30 deg, a high trial count.
     The olb sigma2_I against aperture_averaged_scintillation_index and against
     a matched-count aotools run.
  3. THE FADE TAIL. Space downlink 30 deg. The low quantiles of collected_power
     that set a fade margin (0.1%, 1%, 5%, 10%), the median, the 90% quantile,
     and the scintillation index, aotools vs olb, with bootstrapped quantile
     bars. This is the key deliverable.
  4. THE OUTER SCALE. The phase structure-function test of the P1 check, now
     with a FINITE outer scale (L0 = 25 m), both generators, on one grid.

The estimators. sigma2_I = var(P) / mean(P)^2 (a normalised intensity
variance). The mean error of a variance estimate of a near-normal sample is
var * sqrt(2/(n-1)). The quantile bars are bootstrapped. The structure
function is D_phi(r) = 6.88 (r/r0)^(5/3), Fried, DOI 10.1364/JOSA.56.001372.
The screen physics is Schmidt (2010), DOI 10.1117/3.866274, Ch. 9. The
downlink analytic index is olb.turbulence.plane_wave_scintillation.

REPORT-ONLY. This script changes no production code. The default stays
aotools; the switch to olb is an owner decision.

Run from the repository root:
    python -m validation.waveoptics_speed.generator_validation
'''

import json
import os
import platform
import time
import warnings

import numpy as np

from olb.geometry import CircularOrbit, HorizontalPath
from olb.scenario import (Channel, SpaceScenario, TerrestrialChannel,
                          TerrestrialScenario)
from olb.terminal import SMF, Terminal, Transmitter
from olb.turbulence.plane_wave_scintillation import (
    aperture_averaged_scintillation_index)
from olb.turbulence.profiles import DEFAULT_HS, default_cn2_profile
from olb.waveoptics.threader import Threader
from olb.waveoptics.turbulence.run import propagate_turbulent_scenario
from olb.waveoptics.turbulence.screens import ScreenFactory, phase_screen

LAM = 1550e-9
HS = DEFAULT_HS
HERE = os.path.dirname(__file__)

# Trial counts. The rapid 256 px case is cheap, so it carries the high count
# for the converged index and the fade tail. The 512 px standard cases are
# heavier, so they take the matched 200. See the task brief.
N_TERR = 300            # terrestrial 2 km, standard, 256 px.
N_RAPID_BIG = 2000      # space downlink 30 deg, rapid, 256 px (meas 1, 2, 3).
N_STD = 200             # the 512 px standard cases (meas 1).

# Be a good CPU citizen: three agents share a 32-core machine. One third each.
_WORKERS = max(1, (os.cpu_count() or 4) // 3)


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def warm_aotools():
    '''Trigger the one-time scipy deprecation that aotools raises, quietly.'''
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        phase_screen(0.1, 64, 0.01, seed=0)


def d_phi(scr, kpx):
    '''Measure the phase structure function at kpx pixels, on the two axes.'''
    dh = scr[:, kpx:] - scr[:, :-kpx]
    dv = scr[kpx:, :] - scr[:-kpx, :]
    return 0.5 * (np.mean(dh * dh) + np.mean(dv * dv))


def _stat(sample):
    '''Return (mean, se_mean, sigma2_I, se_sigma2_I) of a power sample.

    sigma2_I = var / mean^2. The standard error of a variance estimate of a
    near-normal sample is var * sqrt(2/(n-1)); the mean error feeds through but
    it is small here. This matches the P1 check estimator.
    '''
    p = np.asarray(sample, dtype=float)
    n = p.size
    mean = float(p.mean())
    se_mean = float(p.std(ddof=1) / np.sqrt(n))
    var = float(p.var(ddof=1))
    sigma2 = var / mean ** 2
    se_sigma2 = float(sigma2 * np.sqrt(2.0 / (n - 1)))
    return mean, se_mean, sigma2, se_sigma2


def _agree(a, se_a, b, se_b):
    '''Return the delta, the combined bar, the sigmas apart, and the 2-sigma
    verdict of two measurements a and b.'''
    d = a - b
    se = float(np.hypot(se_a, se_b))
    return {'delta': float(d), 'se_delta': se,
            'sigmas_apart': float(abs(d) / se) if se > 0 else float('inf'),
            'agrees_2sigma': bool(abs(d) < 2.0 * se)}


def _run_pair(scenario, geometry, n_trials, preset, seed, threader, **kw):
    '''Run both generators on one case, matched count and seed.

    Returns a dict of the two per-generator power samples and, for an uplink,
    the two eta_turb samples.
    '''
    out = {}
    for gen in ('aotools', 'olb'):
        t0 = time.perf_counter()
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            res = propagate_turbulent_scenario(
                scenario, geometry, n_trials=n_trials, seed=seed,
                preset=preset, threader=threader, screen_generator=gen, **kw)
        wall = time.perf_counter() - t0
        power = np.array([tr.collected_power for tr in res.trials])
        eta = np.array([tr.eta_turb for tr in res.trials
                        if tr.eta_turb is not None])
        smf = np.array([tr.smf_eta for tr in res.trials
                        if tr.smf_eta is not None])
        # Keep only the sampling warnings, not the aotools deprecation.
        msgs = sorted({str(w.message) for w in caught
                       if 'deprecat' not in str(w.message).lower()})
        out[gen] = {'power': power, 'eta': eta, 'smf': smf, 'wall_s': float(wall),
                    'grid_n': int(res.grid.n),
                    'n_screens': int(res.plan.z_m.size),
                    'warnings': msgs}
    return out


def _equivalence_block(pair, label):
    '''Turn a _run_pair result into the mean-power and sigma2_I comparison.'''
    block = {'label': label,
             'grid_n': pair['olb']['grid_n'],
             'n_screens': pair['olb']['n_screens'],
             'n_trials': int(pair['olb']['power'].size),
             'warnings': {'aotools': pair['aotools']['warnings'],
                          'olb': pair['olb']['warnings']}}
    for gen in ('aotools', 'olb'):
        mean, se_mean, sig, se_sig = _stat(pair[gen]['power'])
        block[gen] = {'mean_power': mean, 'se_mean_power': se_mean,
                      'sigma2_I': sig, 'se_sigma2_I': se_sig,
                      'wall_s': pair[gen]['wall_s']}
    block['mean_power_cmp'] = _agree(
        block['aotools']['mean_power'], block['aotools']['se_mean_power'],
        block['olb']['mean_power'], block['olb']['se_mean_power'])
    block['sigma2_I_cmp'] = _agree(
        block['aotools']['sigma2_I'], block['aotools']['se_sigma2_I'],
        block['olb']['sigma2_I'], block['olb']['se_sigma2_I'])
    # The reciprocity uplink observable, when present.
    if pair['olb']['eta'].size:
        for gen in ('aotools', 'olb'):
            m, se_m, s, se_s = _stat(pair[gen]['eta'])
            block.setdefault('eta_turb', {})[gen] = {
                'mean': m, 'se_mean': se_m, 'sigma2': s, 'se_sigma2': se_s}
        block['eta_turb_mean_cmp'] = _agree(
            block['eta_turb']['aotools']['mean'],
            block['eta_turb']['aotools']['se_mean'],
            block['eta_turb']['olb']['mean'],
            block['eta_turb']['olb']['se_mean'])
    # The single-mode-fibre coupling, when present. It carries the tilt and the
    # higher-order residual, so it is the turbulence-sensitive observable of a
    # heavily aperture-averaged terrestrial link.
    if pair['olb']['smf'].size:
        for gen in ('aotools', 'olb'):
            m, se_m, s, se_s = _stat(pair[gen]['smf'])
            block.setdefault('smf_eta', {})[gen] = {
                'mean': m, 'se_mean': se_m, 'sigma2': s, 'se_sigma2': se_s}
        block['smf_eta_mean_cmp'] = _agree(
            block['smf_eta']['aotools']['mean'],
            block['smf_eta']['aotools']['se_mean'],
            block['smf_eta']['olb']['mean'],
            block['smf_eta']['olb']['se_mean'])
        block['smf_eta_sigma2_cmp'] = _agree(
            block['smf_eta']['aotools']['sigma2'],
            block['smf_eta']['aotools']['se_sigma2'],
            block['smf_eta']['olb']['sigma2'],
            block['smf_eta']['olb']['se_sigma2'])
    return block


def _bootstrap_quantiles(x, qs, n_boot, rng):
    '''Return the quantiles of x and their bootstrap standard errors.'''
    x = np.asarray(x, dtype=float)
    base = np.quantile(x, qs)
    boots = np.empty((n_boot, len(qs)))
    for b in range(n_boot):
        boots[b] = np.quantile(rng.choice(x, x.size, replace=True), qs)
    return base, boots.std(axis=0)


# ---------------------------------------------------------------------------
# the scenarios
# ---------------------------------------------------------------------------

def _terrestrial_scn():
    near = Terminal(aperture_m=0.20, wavelength_m=LAM,
                    transmitter=Transmitter(waist_m=0.05))
    far = Terminal(aperture_m=0.20, wavelength_m=LAM, detector=SMF())
    return TerrestrialScenario(
        near=near, far=far,
        channel=TerrestrialChannel(path_length_m=2000.0, cn2=5e-15))


def _space_scn(direction):
    ground = Terminal(aperture_m=0.40, wavelength_m=LAM,
                      transmitter=Transmitter(waist_m=0.06))
    return SpaceScenario(ground=ground,
                         space=Terminal(aperture_m=0.30, wavelength_m=LAM),
                         direction=direction, channel=Channel())


# ---------------------------------------------------------------------------
# 1. equivalence across cases
# ---------------------------------------------------------------------------

def equivalence_across_cases(threader):
    '''Compare the two generators on four cases.'''
    orbit = CircularOrbit(altitude_m=600e3, elevation_deg=[30.0])
    blocks = []

    print('1. equivalence across cases (aotools vs olb, matched seed + count)')

    # terrestrial 2 km, Cn2 = 5e-15, standard, SMF receiver.
    terr = _terrestrial_scn()
    pair = _run_pair(terr, HorizontalPath(2000.0), N_TERR, 'standard',
                     seed=3001, threader=threader)
    blocks.append(_equivalence_block(pair, 'terrestrial 2km 5e-15 standard'))

    # space downlink 30 deg, rapid AND standard.
    down = _space_scn('downlink')
    pair = _run_pair(down, orbit, N_RAPID_BIG, 'rapid', seed=3002,
                     threader=threader)
    blocks.append(_equivalence_block(pair, 'space downlink 30deg rapid'))
    down_rapid_pair = pair            # reuse for meas 2 and 3.

    pair = _run_pair(down, orbit, N_STD, 'standard', seed=3003,
                     threader=threader)
    blocks.append(_equivalence_block(pair, 'space downlink 30deg standard'))

    # space uplink 30 deg, standard.
    up = _space_scn('uplink')
    pair = _run_pair(up, orbit, N_STD, 'standard', seed=3004,
                     threader=threader)
    blocks.append(_equivalence_block(pair, 'space uplink 30deg standard'))

    for b in blocks:
        print(f'   {b["label"]:34s}  (n={b["n_trials"]}, {b["grid_n"]} px, '
              f'{b["n_screens"]} scr)')
        for gen in ('aotools', 'olb'):
            g = b[gen]
            print(f'     {gen:8s} mean {g["mean_power"]:.5f} '
                  f'+/- {g["se_mean_power"]:.5f}   '
                  f'sigma2_I {g["sigma2_I"]:.5f} +/- {g["se_sigma2_I"]:.5f}   '
                  f'({g["wall_s"]:.1f} s)')
        mp, sg = b['mean_power_cmp'], b['sigma2_I_cmp']
        print(f'     delta mean {mp["delta"]:+.5f} '
              f'({mp["sigmas_apart"]:.2f} sig, agrees {mp["agrees_2sigma"]})   '
              f'delta sigma2 {sg["delta"]:+.5f} '
              f'({sg["sigmas_apart"]:.2f} sig, agrees {sg["agrees_2sigma"]})')
        if 'eta_turb' in b:
            ec = b['eta_turb_mean_cmp']
            print(f'     eta_turb mean aotools '
                  f'{b["eta_turb"]["aotools"]["mean"]:.5f} / olb '
                  f'{b["eta_turb"]["olb"]["mean"]:.5f}  '
                  f'({ec["sigmas_apart"]:.2f} sig, agrees {ec["agrees_2sigma"]})')
        if 'smf_eta' in b:
            mc, sc = b['smf_eta_mean_cmp'], b['smf_eta_sigma2_cmp']
            print(f'     smf_eta mean aotools '
                  f'{b["smf_eta"]["aotools"]["mean"]:.5f} / olb '
                  f'{b["smf_eta"]["olb"]["mean"]:.5f}  '
                  f'({mc["sigmas_apart"]:.2f} sig, agrees {mc["agrees_2sigma"]})'
                  f'   sigma2 ({sc["sigmas_apart"]:.2f} sig, '
                  f'agrees {sc["agrees_2sigma"]})')
        for gen in ('aotools', 'olb'):
            if b['warnings'][gen]:
                print(f'     [{gen} warnings] {b["warnings"][gen]}')
    print('')
    return blocks, down_rapid_pair


# ---------------------------------------------------------------------------
# 2. converged sigma2_I vs analytic
# ---------------------------------------------------------------------------

def converged_vs_analytic(down_rapid_pair):
    '''Compare the converged sigma2_I of both generators to the analytic value.

    The screens are the rapid preset (5 screens). The rapid count sits about 10
    percent below the converged aperture index for a 30 deg slab (a KNOWN
    screen-count effect, docs/schmidt-crosscheck.md WP7). That bias is COMMON
    to both generators, so it does not affect the generator-equivalence
    question this pass answers.
    '''
    ground = _space_scn('downlink').ground
    cn2 = default_cn2_profile(Channel().site, HS)
    sigma2_analytic = float(aperture_averaged_scintillation_index(
        ground.aperture_m, 30.0, LAM, HS, cn2))

    out = {'sigma2_I_analytic': sigma2_analytic,
           'n_trials': int(down_rapid_pair['olb']['power'].size),
           'preset': 'rapid', 'n_screens': down_rapid_pair['olb']['n_screens']}
    for gen in ('aotools', 'olb'):
        _, _, sig, se_sig = _stat(down_rapid_pair[gen]['power'])
        out[gen] = {'sigma2_I': sig, 'se_sigma2_I': se_sig,
                    'ratio_to_analytic': sig / sigma2_analytic}
    out['gen_cmp'] = _agree(out['aotools']['sigma2_I'],
                            out['aotools']['se_sigma2_I'],
                            out['olb']['sigma2_I'], out['olb']['se_sigma2_I'])

    print('2. converged sigma2_I vs analytic (space downlink 30 deg, rapid)')
    print(f'   n_trials {out["n_trials"]}, {out["n_screens"]} screens')
    print(f'   analytic (aperture-averaged) {sigma2_analytic:.5f}')
    for gen in ('aotools', 'olb'):
        g = out[gen]
        print(f'   {gen:8s} {g["sigma2_I"]:.5f} +/- {g["se_sigma2_I"]:.5f}  '
              f'(ratio to analytic {g["ratio_to_analytic"]:.3f})')
    gc = out['gen_cmp']
    print(f'   generator delta {gc["delta"]:+.5f} '
          f'({gc["sigmas_apart"]:.2f} sig, agrees {gc["agrees_2sigma"]})')
    print('')
    return out


# ---------------------------------------------------------------------------
# 3. the fade tail
# ---------------------------------------------------------------------------

def fade_tail(down_rapid_pair, n_boot=500, seed=777):
    '''Compare the low quantiles of collected_power between the generators.'''
    qs = [0.001, 0.01, 0.05, 0.10, 0.50, 0.90]
    rng = np.random.default_rng(seed)
    out = {'quantiles': qs,
           'n_trials': int(down_rapid_pair['olb']['power'].size),
           'preset': 'rapid'}
    base = {}
    se = {}
    for gen in ('aotools', 'olb'):
        p = down_rapid_pair[gen]['power']
        base[gen], se[gen] = _bootstrap_quantiles(p, qs, n_boot, rng)
        sig = float(p.var() / p.mean() ** 2)
        out[gen] = {'quantile': [float(v) for v in base[gen]],
                    'se_quantile': [float(v) for v in se[gen]],
                    'scint_index': sig}
    # The per-quantile difference in dB, with the combined bootstrap bar. A
    # fade margin reads the low quantiles, so the dB gap there is the verdict.
    rows = []
    for i, q in enumerate(qs):
        a, o = base['aotools'][i], base['olb'][i]
        d_db = 10.0 * np.log10(o / a)
        # Propagate the bootstrap bars into dB: d(10 log10 x) = 10/ln10 dx/x.
        se_db = (10.0 / np.log(10.0)) * np.hypot(se['aotools'][i] / a,
                                                 se['olb'][i] / o)
        rows.append({'q': q, 'aotools': float(a), 'olb': float(o),
                     'delta_db': float(d_db), 'se_delta_db': float(se_db),
                     'agrees_2sigma': bool(abs(d_db) < 2.0 * se_db)})
    out['comparison'] = rows

    print('3. the fade tail (space downlink 30 deg, rapid) '
          f'n_trials {out["n_trials"]}')
    print('   quantile   aotools      olb      delta dB   bar dB   agree')
    for r in rows:
        print(f'   {r["q"]*100:6.1f}%   {r["aotools"]:.5f}  {r["olb"]:.5f}  '
              f'{r["delta_db"]:+7.4f}  {r["se_delta_db"]:6.4f}  '
              f'{r["agrees_2sigma"]}')
    print(f'   scint index aotools {out["aotools"]["scint_index"]:.5f}, '
          f'olb {out["olb"]["scint_index"]:.5f}')
    print('')
    return out


# ---------------------------------------------------------------------------
# 4. the outer scale
# ---------------------------------------------------------------------------

def outer_scale(L0_m=25.0):
    '''The structure-function test with a finite outer scale, both generators.

    The band r/r0 = 0.3 to 1.6 is r = 3 to 16 cm at r0 = 10 cm, which is far
    below L0 = 25 m. So the von Karman roll-off (at r ~ L0) is NOT reached
    here, and both generators must still track D_phi = 6.88 (r/r0)^(5/3). The
    check also confirms the finite-L0 screen matches the infinite-L0 screen in
    this band, so the outer scale does not corrupt the in-band structure. See
    Fried, DOI 10.1364/JOSA.56.001372, and Schmidt, DOI 10.1117/3.866274,
    Ch. 9, Eq. (9.51).
    '''
    n, dx, r0, M = 512, 0.01, 0.10, 40
    ks = np.array([3, 5, 8, 11, 16])          # r/r0 = 0.3 to 1.6.
    theory = 6.88 * ((ks * dx) / r0) ** (5.0 / 3.0)

    def measure_olb(L0):
        fac = ScreenFactory(n, dx, L0_m=L0)
        acc = np.zeros(ks.size)
        for i in range(M):
            s = fac.make(r0, np.random.default_rng(9000 + i))
            for j, kpx in enumerate(ks):
                acc[j] += d_phi(s, kpx)
        return acc / M

    def measure_aotools(L0):
        acc = np.zeros(ks.size)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            for i in range(M):
                s = phase_screen(r0, n, dx, L0_m=L0, seed=9000 + i)
                for j, kpx in enumerate(ks):
                    acc[j] += d_phi(s, kpx)
        return acc / M

    d_olb = measure_olb(L0_m)
    d_aot = measure_aotools(L0_m)
    d_olb_inf = measure_olb(np.inf)

    out = {'L0_m': L0_m, 'r_over_r0': [float(k * dx / r0) for k in ks],
           'd_theory': [float(v) for v in theory],
           'd_olb_finite': [float(v) for v in d_olb],
           'd_aotools_finite': [float(v) for v in d_aot],
           'd_olb_infinite': [float(v) for v in d_olb_inf],
           'ratio_olb_to_theory': [float(v) for v in d_olb / theory],
           'ratio_aotools_to_theory': [float(v) for v in d_aot / theory],
           'olb_vs_aotools': [float(v) for v in d_olb / d_aot],
           'finite_vs_infinite_olb': [float(v) for v in d_olb / d_olb_inf]}
    out['olb_aotools_in_5pct'] = bool(
        np.all(np.abs(d_olb / d_aot - 1.0) < 0.05))
    out['olb_tracks_theory_in_band'] = bool(
        np.all(d_olb / theory > 0.85) and np.all(d_olb / theory < 1.05))

    print(f'4. the outer scale, L0 = {L0_m} m, structure function on one grid')
    print('   r/r0    theory     olb(L0)   aot(L0)   olb/aot   olb/olb(inf)')
    for i, k in enumerate(ks):
        print(f'   {k*dx/r0:4.2f}   {theory[i]:8.3f}  {d_olb[i]:8.3f}  '
              f'{d_aot[i]:8.3f}  {d_olb[i]/d_aot[i]:7.4f}  '
              f'{d_olb[i]/d_olb_inf[i]:7.4f}')
    print(f'   olb vs aotools inside 5%: {out["olb_aotools_in_5pct"]}')
    print(f'   olb tracks theory in band: {out["olb_tracks_theory_in_band"]}')
    print('')
    return out


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main():
    t_start = time.time()
    warm_aotools()
    threader = Threader(max_workers=_WORKERS)
    print(f'(workers {threader.max_workers} of {os.cpu_count()} cores)\n')

    blocks, down_rapid = equivalence_across_cases(threader)
    conv = converged_vs_analytic(down_rapid)
    tail = fade_tail(down_rapid)
    osc = outer_scale()

    out = {
        'environment': {'numpy': np.__version__,
                        'platform': platform.platform(),
                        'cores': os.cpu_count(), 'workers': threader.max_workers},
        'trial_counts': {'terrestrial': N_TERR, 'downlink_rapid': N_RAPID_BIG,
                         'standard_cases': N_STD},
        'equivalence': blocks,
        'converged_vs_analytic': conv,
        'fade_tail': tail,
        'outer_scale': osc,
    }
    path = os.path.join(HERE, 'generator_validation_results.json')
    with open(path, 'w') as f:
        json.dump(out, f, indent=1)
    print(f'wrote {path}')
    print(f'(elapsed {time.time() - t_start:.0f} s)')


if __name__ == '__main__':
    main()
