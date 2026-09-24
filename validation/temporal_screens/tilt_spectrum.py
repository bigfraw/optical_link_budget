'''
Gate (d) of the temporal strip screens: the tilt spectrum against Greenwood.

A record of frames only holds a TIME AXIS if the phase it carries changes at
the right RATE. This script measures that rate. It runs ONE frozen-flow record
of the hero downlink, it reads the SUMMED SCREEN PHASE that each frame stores,
it fits the Z-tilt of the 0.7 m ground pupil in every frame, and it takes the
power spectrum of that time series.

WHAT THE SPECTRUM MUST SHOW. The one-axis tilt of a pupil under frozen flow has
two power laws. Below the pupil corner frequency the spectrum falls slowly, as
f^(-2/3); above it the pupil averages the small scales away and the spectrum
falls steeply. The break between the two is the CORNER. The script fits both
laws, it crosses them to find the corner, and it compares that corner with two
reference frequencies:

  1. THE GREENWOOD FREQUENCY, `olb.turbulence.andrews.temporal
     .greenwood_frequency`, fed with the SAME per-layer velocity array that the
     strips use (`StripPlan.v_x_m_s` and `v_y_m_s`). It is NOT fed its default
     Bufton wind profile: that profile carries its own slew term, and this
     record already holds the slew in the layer velocity, so the default value
     reads high. The script prints both, so the reader sees the size of that
     trap.
  2. THE PUPIL CORNER `0.3 v / D` of the FASTEST layer. The Greenwood frequency
     is a PHASE quantity: it is set by the smallest scale of the wavefront, and
     it does not know the pupil diameter. The TILT is a pupil quantity, so its
     corner sits well below f_G, on the V/D scale. Source: Tyler, J. Opt. Soc.
     Am. A 11(1), pp. 358 to 367 (1994), DOI 10.1364/JOSAA.11.000358. The
     script also prints the Tyler tilt-tracking frequency f_T and the corner of
     the slowest layer, because a five-layer stack whose speeds span 11 to
     131 m/s has FIVE corners and its break is a blend of them.

The band on the corner is 30 percent on the ratio, and the band on the
low-frequency exponent is 0.30 on the exponent itself.

WHAT THIS GATE FOUND (2026-09-13, 4000 frames at 0.5 ms, `rapid`, 30 deg). The
low-frequency exponent reads -0.56 against the -2/3 = -0.667 law, so the
frozen-flow axis carries the right tilt spectrum. The break sits at 38 Hz. That
is 0.06 of f_G = 642 Hz, so f_G is NOT the tilt corner, and it is 0.68 of the
56 Hz corner of the fastest layer, which is where a blend of five layers must
sit. So the two corner bands FAIL as written, and the reason is physics, not
the strips: the tilt of a 0.7 m pupil breaks on the V/D scale.

THE PRESET IS `rapid`, so the record runs locally in a few minutes. The
velocity model and the sample step do not depend on the preset.

Sources:
- Taylor, Proc. R. Soc. Lond. A 164, pp. 476 to 490 (1938),
  DOI 10.1098/rspa.1938.0032. The frozen-flow hypothesis.
- Greenwood, J. Opt. Soc. Am. 67(3), pp. 390 to 393 (1977),
  DOI 10.1364/JOSA.67.000390. The Greenwood frequency.
- Tyler, J. Opt. Soc. Am. A 11(1), pp. 358 to 367 (1994),
  DOI 10.1364/JOSAA.11.000358. The tilt power spectrum and its corner.
- Andrews and Phillips, Laser Beam Propagation through Random Media, 2nd ed.,
  DOI 10.1117/3.626196, Ch. 14, Eq. (38), printed p. 622 (the Greenwood
  frequency), and Ch. 12, Eqs. (2) and (3), printed p. 481 (the Bufton wind).
- Noll, J. Opt. Soc. Am. 66(3), pp. 207 to 211 (1976),
  DOI 10.1364/JOSA.66.000207. The Zernike tilt modes Z2 and Z3.

Outputs, next to this script:
    data/tilt_series.csv     the tilt of every frame
    data/tilt_psd.csv        the averaged power spectrum
    figures/tilt_spectrum.png the spectrum with the two fits and the two
                              reference frequencies

Run from the repo root:
    python -m validation.temporal_screens.tilt_spectrum
'''

import argparse
import csv
import os
import shutil
import time
import warnings

import numpy as np

from olb.turbulence.andrews.temporal import greenwood_frequency
from olb.turbulence.profiles import default_cn2_profile
from olb.waveoptics.compensation import ApertureModes, circle
from olb.waveoptics.turbulence.run import propagate_turbulent_scenario
from olb.waveoptics.turbulence.temporal import TemporalSpec, strip_plan
from validation.waveoptics_ao.waveoptics_ao import hero_scenario

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, 'data')
FIGS = os.path.join(HERE, 'figures')

# ---- the case ----
ELEVATION_DEG = 30.0
PRESET = 'rapid'
PRECISION = 'single'
SEED = 20260913
DT_S = 5e-4                  # the sample step, in s.
PATCH_RADIUS_M = 0.35        # half the hero ground aperture.
N_SEGMENTS = 2               # the segments of the averaged periodogram.
LOW_BAND_HZ = (1.0, 5.0)     # the fit band of the f^(-2/3) law.
HIGH_BAND_HZ = (60.0, 400.0)  # the fit band of the steep law.
SLOPE_TOL = 0.30             # the band on the low-frequency exponent.
CORNER_TOL = 0.30            # the band on the corner frequency.

BANDS = []


def band(name, value, low, high):
    '''Record one PASS or FAIL line, and print it.'''
    ok = low <= value <= high
    line = (f"  {'PASS' if ok else 'FAIL'}  {name:46s} {value:9.4f}  "
            f"band [{low:.4f}, {high:.4f}]")
    BANDS.append(line)
    print(line)
    return ok


def layer_speeds(sp):
    '''Give the speed of every layer of a StripPlan, in m/s.'''
    return np.hypot(np.asarray(sp.v_x_m_s, dtype=float),
                    np.asarray(sp.v_y_m_s, dtype=float))


def layer_heights(plan, elevation_deg):
    '''Give the altitude of every screen of a downlink plan, in m.

    A space plan counts z from the TOP of the slab, so the distance of a screen
    from the ground is z_total - z, and the altitude is that distance times
    sin(elevation). This is the rule of
    olb.waveoptics.turbulence.run._ground_distance.
    '''
    z_g = float(plan.z_total_m) - np.asarray(plan.z_m, dtype=float)
    return z_g * np.sin(np.deg2rad(float(elevation_deg)))


def tilt_series(result, aperture_m):
    '''Give the two Noll tilt coefficients of every frame, in rad.

    The runner stores the SUMMED SCREEN PHASE of each frame at the patch
    pixels. A space slab starts from a plane wave, so that sum IS the sensed
    wavefront (the runner uses the same source for its own correction). The fit
    is the least-squares projection on the Noll modes Z2 and Z3 over the ground
    pupil. Source: Noll, DOI 10.1364/JOSA.66.000207, Table I.

    Args:
        result:     the TurbWaveResult of the record.
        aperture_m: the ground aperture diameter, in m.

    Returns:
        A float array of the shape (n_frames, 2): the Z2 and the Z3
        coefficient of each frame.
    '''
    grid, patch = result.grid, result.patch
    modes = ApertureModes(3, grid.n,
                          circle(grid.n, aperture_m / grid.pixel_m, 0.0))
    full = np.zeros(int(grid.n) ** 2, dtype=float)
    out = np.empty((result.screen_phase.shape[0], 2))
    for k, row in enumerate(result.screen_phase):
        full[patch.indices] = row
        coeffs = modes.estimate(full[modes.indices])
        out[k] = coeffs[1], coeffs[2]
    return out


def averaged_psd(series, dt_s, n_segments=N_SEGMENTS):
    '''Give the one-sided power spectrum of a series, averaged over segments.

    The segments are disjoint, each one is detrended and multiplied by a Hann
    window, and the window loss is divided out. This is the plain Welch
    estimate, written out because scipy.signal is not a dependency of olb.

    Args:
        series:      the time series.
        dt_s:        the sample step, in s.
        n_segments:  the number of disjoint segments.

    Returns:
        The pair (frequency [Hz], power spectral density [unit^2/Hz]). The
        zero frequency is dropped.
    '''
    x = np.asarray(series, dtype=float)
    m = int(x.size // int(n_segments))
    win = np.hanning(m)
    scale = dt_s * (win ** 2).sum()
    acc = None
    for s in range(int(n_segments)):
        seg = x[s * m:(s + 1) * m]
        seg = seg - seg.mean()
        spec = np.abs(np.fft.rfft(seg * win)) ** 2 / scale
        spec[1:-1] *= 2.0                      # the one-sided fold.
        acc = spec if acc is None else acc + spec
    freq = np.fft.rfftfreq(m, d=dt_s)
    return freq[1:], (acc / int(n_segments))[1:]


def fit_power_law(freq, psd, lo_hz, hi_hz):
    '''Fit log(psd) = a + b log(freq) over one band. Give (b, 10^a).'''
    keep = (freq >= lo_hz) & (freq <= hi_hz) & (psd > 0.0)
    b, a = np.polyfit(np.log10(freq[keep]), np.log10(psd[keep]), 1)
    return float(b), float(10.0 ** a)


def write_series_csv(times, tilt):
    '''Write the tilt of every frame. Give the path back.'''
    os.makedirs(DATA, exist_ok=True)
    path = os.path.join(DATA, 'tilt_series.csv')
    with open(path, 'w', newline='', encoding='utf-8') as fh:
        w = csv.writer(fh)
        w.writerow(['t_s', 'a2_rad', 'a3_rad'])
        for t, (a2, a3) in zip(times, tilt):
            w.writerow([t, a2, a3])
    return path


def write_psd_csv(freq, psd):
    '''Write the averaged spectrum. Give the path back.'''
    os.makedirs(DATA, exist_ok=True)
    path = os.path.join(DATA, 'tilt_psd.csv')
    with open(path, 'w', newline='', encoding='utf-8') as fh:
        w = csv.writer(fh)
        w.writerow(['f_hz', 'psd_rad2_per_hz'])
        for f, p in zip(freq, psd):
            w.writerow([f, p])
    return path


def draw_psd(freq, psd, low_fit, high_fit, corner_hz, f_g, f_fast):
    '''Draw the spectrum with the two fits. Give the path back.'''
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    os.makedirs(FIGS, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7.5, 5.0))
    ax.loglog(freq, psd, lw=0.8, color='0.4', label='measured')
    for (slope, amp), (lo, hi), style, name in (
            (low_fit, LOW_BAND_HZ, '-', 'low-f fit'),
            (high_fit, HIGH_BAND_HZ, '--', 'high-f fit')):
        f = np.geomspace(lo, hi, 32)
        ax.loglog(f, amp * f ** slope, style, lw=2.0,
                  label=f'{name}, exponent {slope:.2f}')
    for value, name, colour in ((corner_hz, 'corner', 'C2'),
                                (f_fast, 'fastest 0.3 v/D', 'C1'),
                                (f_g, 'Greenwood f_G', 'C3')):
        ax.axvline(value, color=colour, ls=':', lw=1.5,
                   label=f'{name} = {value:.1f} Hz')
    ax.set_xlabel('frequency [Hz]')
    ax.set_ylabel('tilt power spectral density [rad^2/Hz]')
    ax.set_title(f'gate (d): the Z-tilt spectrum of one record, hero downlink '
                 f'{ELEVATION_DEG:.0f} deg')
    ax.grid(True, which='both', alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    path = os.path.join(FIGS, 'tilt_spectrum.png')
    fig.savefig(path, dpi=140)
    plt.close(fig)
    return path


def main():
    '''Run the gate, print the bands, and write the outputs.'''
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--frames', type=int, default=4000,
                        help='the frames of the record (default 4000).')
    args = parser.parse_args()

    t_start = time.time()
    scn, geom = hero_scenario(ELEVATION_DEG)
    strip_root = os.path.join(HERE, 'strips_tilt')
    shutil.rmtree(strip_root, ignore_errors=True)
    spec = TemporalSpec(dt_s=DT_S, n_frames=args.frames, strip_dir=strip_root)

    print(f'gate (d), the hero downlink at {ELEVATION_DEG:.0f} deg, '
          f'{PRESET} preset, {PRECISION} precision, seed {SEED}')
    print(f'  frames                  {args.frames:11d} at '
          f'{DT_S * 1e3:.2f} ms ({(args.frames - 1) * DT_S:.2f} s)')
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        try:
            res = propagate_turbulent_scenario(
                scn, geom, n_trials=args.frames, seed=SEED, preset=PRESET,
                precision=PRECISION, temporal=spec,
                patch_radius_m=PATCH_RADIUS_M, store_screen_phase=True)
            sp = strip_plan(res.plan, res.grid, spec, geom, 25.0)
        finally:
            shutil.rmtree(strip_root, ignore_errors=True)
    run_s = time.time() - t_start

    aperture_m = scn.ground.aperture_m
    speeds = layer_speeds(sp)
    heights = layer_heights(res.plan, ELEVATION_DEG)
    print(f'  grid                    {res.grid.n:11d} px, '
          f'{res.grid.pixel_m * 1e3:.2f} mm')
    print(f'  screens                 {res.plan.z_m.size:11d}')
    print(f'  record run time         {run_s:11.1f} s')
    print('')
    print('  the layer velocity of the strips (slew + Bufton, ws = 0)')
    print(f"  {'layer':>6}{'h [m]':>12}{'v [m/s]':>10}{'r0 [cm]':>10}")
    for j, (h, v, r0) in enumerate(zip(heights, speeds, res.plan.r0_m)):
        print(f'  {j:>6}{h:>12.0f}{v:>10.2f}{r0 * 1e2:>10.2f}')

    # THE REFERENCE FREQUENCIES. The Greenwood integral reads the site Cn2 on a
    # fine height grid, and the SAME velocity model that the strips use,
    # interpolated onto that grid.
    hs = np.concatenate(([0.0], np.geomspace(1.0, float(heights.max()), 512)))
    cn2 = default_cn2_profile(scn.channel.site, hs)
    wind = np.interp(hs, heights[::-1], speeds[::-1])
    f_g = float(greenwood_frequency(hs, cn2, scn.ground.wavelength_m,
                                    elevation_deg=ELEVATION_DEG,
                                    wind_profile=wind))
    f_g_default = float(greenwood_frequency(hs, cn2, scn.ground.wavelength_m,
                                            elevation_deg=ELEVATION_DEG))
    # THE TYLER TILT FREQUENCY, the tracking bandwidth of the Z-tilt of a
    # pupil of the diameter D:
    #     f_T = 0.331 D^(-1/6) lambda^(-1) [INT Cn2(z) v(z)^2 dz]^(1/2)
    # Source: Tyler, DOI 10.1364/JOSAA.11.000358. The plan already holds the
    # SLANT integrated Cn2 of each screen, so the path integral is the sum over
    # the screens with the strip speed of each one.
    w = np.asarray(res.plan.cn2_int_m13, dtype=float)
    f_tyler = float(0.331 * aperture_m ** (-1.0 / 6.0)
                    / scn.ground.wavelength_m
                    * np.sqrt((w * speeds ** 2).sum()))
    v_eff = float((w * speeds ** (5.0 / 3.0)).sum() / w.sum()) ** (3.0 / 5.0)
    # The pupil corner of the SLOWEST and of the FASTEST layer. Each layer
    # breaks at its own 0.3 v / D, so the two numbers bracket the break of the
    # whole stack. Source: Tyler, DOI 10.1364/JOSAA.11.000358 (the V/D scale).
    f_slow = 0.3 * float(speeds.min()) / aperture_m
    f_fast = 0.3 * float(speeds.max()) / aperture_m

    print('')
    print(f'  Greenwood f_G, strip velocity  {f_g:10.1f} Hz')
    print(f'  Greenwood f_G, default wind    {f_g_default:10.1f} Hz '
          f'({f_g_default / f_g:.1f}x, it carries its own slew)')
    print(f'  effective layer speed          {v_eff:10.1f} m/s')
    print(f'  Tyler tilt frequency f_T       {f_tyler:10.1f} Hz')
    print(f'  pupil corner, slowest layer    {f_slow:10.1f} Hz')
    print(f'  pupil corner, fastest layer    {f_fast:10.1f} Hz')

    tilt = tilt_series(res, aperture_m)
    times = np.arange(tilt.shape[0]) * DT_S
    freq, psd_x = averaged_psd(tilt[:, 0], DT_S)
    _f, psd_y = averaged_psd(tilt[:, 1], DT_S)
    psd = 0.5 * (psd_x + psd_y)          # the two axes of one tilt.

    low_fit = fit_power_law(freq, psd, *LOW_BAND_HZ)
    high_fit = fit_power_law(freq, psd, *HIGH_BAND_HZ)
    # The corner is where the two power laws cross.
    corner_hz = float(10.0 ** (np.log10(low_fit[1] / high_fit[1])
                               / (high_fit[0] - low_fit[0])))

    print('')
    print('1. the two power laws of the tilt spectrum')
    print(f'  {LOW_BAND_HZ[0]:.0f} to {LOW_BAND_HZ[1]:.0f} Hz     '
          f'exponent {low_fit[0]:8.3f}  (the law is -2/3 = -0.667)')
    print(f'  {HIGH_BAND_HZ[0]:.0f} to {HIGH_BAND_HZ[1]:.0f} Hz  '
          f'exponent {high_fit[0]:8.3f}')
    band('low-frequency exponent', low_fit[0],
         -2.0 / 3.0 - SLOPE_TOL, -2.0 / 3.0 + SLOPE_TOL)

    print('')
    print('2. the corner of the spectrum against the two references')
    print(f'  measured corner                {corner_hz:10.1f} Hz')
    band('corner / Greenwood f_G', corner_hz / f_g,
         1.0 - CORNER_TOL, 1.0 + CORNER_TOL)
    band('corner / fastest-layer 0.3 v/D', corner_hz / f_fast,
         1.0 - CORNER_TOL, 1.0 + CORNER_TOL)

    print('')
    print(f'  file saved: {write_series_csv(times, tilt)}')
    print(f'  file saved: {write_psd_csv(freq, psd)}')
    print(f'  figure saved: '
          f'{draw_psd(freq, psd, low_fit, high_fit, corner_hz, f_g, f_fast)}')
    print('    Caption: the averaged Z-tilt power spectrum of one frozen-flow '
          'record, with')
    print('    the two fitted power laws, their crossing, and the two '
          'reference frequencies.')

    failed = [line for line in BANDS if 'FAIL' in line]
    print('')
    print(f'  {len(BANDS) - len(failed)} of {len(BANDS)} pass bands hold')
    for line in failed:
        print(line)
    print('')
    print(f'(elapsed {time.time() - t_start:.1f} s)')


if __name__ == '__main__':
    main()
