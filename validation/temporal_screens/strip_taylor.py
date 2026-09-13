'''
Gate (b) of the temporal strip screens: FROZEN FLOW on one strip.

A record crops its frames out of one long strip at the integer pixel offset of
the absolute position v*k*dt (`olb/waveoptics/turbulence/temporal.py`). This
script measures that the time axis which comes out of that crop is the frozen
flow of Taylor, DOI 10.1098/rspa.1938.0032.

THREE MEASUREMENTS.

  1. TAYLOR. The temporal structure function D(tau) of ONE pixel of the frame
     must equal the SPATIAL structure function D(v*tau) of the same strip. The
     band is 2 standard errors. The residual is the integer-pixel rounding of
     the offset, which is at most half a pixel.
  2. THE SPECTRUM. The temporal phase power spectrum of one pixel must follow
     f^(-8/3) over the inertial decade. Andrews and Phillips,
     DOI 10.1117/3.626196, Ch. 12, give the frozen-flow time axis; the -8/3 law
     is the Kolmogorov phase spectrum f^(-11/3) after the integration over the
     axis across the wind. The script prints the fitted exponent and its 95
     percent confidence interval.
  3. THE SEAM. A Fourier screen is periodic, so the strip carries a pad of
     `pad_outer_scales` outer scales, and the record must not walk back into
     its own start. The test is the correlation of two windows that sit the
     WHOLE travel apart. ONE such pair says nothing: a 1.28 m window of a
     screen with a 25 m outer scale is almost a plane, so after the mean of
     each window goes, the correlation is little more than the sign agreement
     of two tilts, and it lands near +1 or -1 at random. So the script takes
     the MEAN over every disjoint pair of the strip, and it asks that the mean
     is zero inside 2 standard errors. It also checks, exactly, that the last
     window ends before the pad.

Sources:
- Taylor, The spectrum of turbulence, DOI 10.1098/rspa.1938.0032. Frozen flow.
- Andrews and Phillips, Laser Beam Propagation through Random Media, 2nd ed.,
  DOI 10.1117/3.626196, Ch. 12, Eqs. (2) and (3), printed p. 481 (the Bufton
  wind), and Ch. 12 (the temporal spectra of a frozen atmosphere).
- Greenwood, DOI 10.1364/JOSA.67.000390. The Greenwood frequency: the time
  scale that the step dt must resolve.
- Schmidt, DOI 10.1117/3.866274, Ch. 9, Eqs. (9.78) to (9.81), printed pp. 166
  to 169. The screen that the strip is made of.
- Assemat and Wilson, DOI 10.1364/OE.14.000988, Eq. (5). The von Karman
  covariance that bounds the seam.

Outputs, next to this script:
    data/taylor_dphi.csv       D(tau) against the spatial D(v tau)
    data/taylor_psd.csv        the temporal power spectrum
    figures/strip_taylor.png   the two panels

Run from the repo root:
    python -m validation.temporal_screens.strip_taylor
'''

import os
import shutil
import time
import types

import numpy as np
from scipy.signal import welch

from olb.waveoptics.turbulence.temporal import (TemporalSpec, build_strips,
                                                frame_offsets, frame_stack,
                                                open_strips, strip_paths,
                                                strip_plan)

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, 'data')
FIGS = os.path.join(HERE, 'figures')
STRIPS = os.path.join(HERE, 'strips')

# ---- the case ----
N = 128                      # the frame side, in pixels.
DX = 0.01                    # the pitch, in m.
R0 = 0.10                    # the Fried parameter of the layer, in m.
L0 = 25.0                    # the outer scale, in m.
DT = 0.5e-3                  # the step, in s.
N_FRAMES = 4096
WIND = 50.0                  # the ground wind Vg, in m/s.
SEED = 20260913
LAGS = np.array([4, 8, 16, 32, 64, 128])      # frames.
FIT_BAND_HZ = (10.0, 200.0)                   # the inertial decade of the fit.
SEAM_MAX = 0.01

BANDS = []


def band(name, value, low, high):
    '''Record one PASS or FAIL line, and print it.'''
    ok = low <= value <= high
    line = (f"  {'PASS' if ok else 'FAIL'}  {name:46s} {value:10.4f}  "
            f"band [{low:.4f}, {high:.4f}]")
    BANDS.append(line)
    print(line)
    return ok


def build():
    '''Make the one-layer strip of the record, and open it.'''
    # A one-layer horizontal plan: the module reads four attributes only, so a
    # namespace is enough and this script needs no scenario.
    plan = types.SimpleNamespace(z_m=np.array([0.0]), z_total_m=0.0,
                                 direction='terrestrial',
                                 r0_m=np.array([R0]))
    grid = types.SimpleNamespace(n=N, pixel_m=DX)
    geometry = types.SimpleNamespace(elevation_deg=90.0, slew_deg_s=0.0)
    spec = TemporalSpec(dt_s=DT, n_frames=N_FRAMES, strip_dir=STRIPS,
                        record=0, wind_ground_m_s=WIND, wind_dir_deg=0.0,
                        slew_rad_s=0.0, pad_outer_scales=2.0)
    sp = strip_plan(plan, grid, spec, geometry, L0)
    paths = strip_paths(spec, 1)
    build_strips(sp, plan.r0_m, paths, L0, True, [SEED])
    return sp, paths


def sample(sp, strips):
    '''Give the time series of one column of the frame, one row at a time.

    The result is an (n_frames, n) array: the column n//2 of each frame. The
    rows are independent samples of the same process, so the ensemble is n
    series and not one.
    '''
    out = np.empty((sp.n_frames, sp.n), dtype=np.float32)
    for k in range(sp.n_frames):
        out[k] = next(iter(frame_stack(strips, sp, k)))[:, sp.n // 2]
    return out


def spatial_dphi(strip, lag_px):
    '''Give D(r) along x of the strip at one integer pixel lag.'''
    k = int(lag_px)
    diff = np.asarray(strip[:, k:], dtype=np.float64) - strip[:, :-k]
    return float(np.mean(diff * diff)), float(np.std(diff * diff)
                                              / np.sqrt(diff.size))


def temporal_dphi(series, lag):
    '''Give D(tau) of the time series at one frame lag, and its error.'''
    diff = series[lag:, :].astype(np.float64) - series[:-lag, :]
    sq = diff * diff
    # The samples along the time axis are correlated, so the effective count
    # is the number of INDEPENDENT rows times the travel in coherence lengths.
    n_eff = series.shape[1] * max(1.0, (series.shape[0] - lag) * lag
                                  / series.shape[0])
    return float(np.mean(sq)), float(np.std(sq) / np.sqrt(n_eff))


def fit_slope(f, psd, low, high):
    '''Fit log(psd) = a + b log(f) over a band, and give b and its 95 percent
    confidence half width.

    The fit is a closed-form least squares on two vectors. It does NOT call a
    LAPACK routine, so it does not depend on the BLAS of the environment.
    '''
    keep = (f >= low) & (f <= high) & (psd > 0.0)
    x = np.log10(f[keep])
    y = np.log10(psd[keep])
    xm, ym = x.mean(), y.mean()
    sxx = float(np.sum((x - xm) ** 2))
    slope = float(np.sum((x - xm) * (y - ym)) / sxx)
    resid = y - (ym + slope * (x - xm))
    s2 = float(np.sum(resid ** 2) / (x.size - 2))
    return slope, 1.96 * np.sqrt(s2 / sxx), int(x.size)


def write_csv(name, header, rows):
    '''Write one CSV next to this script.'''
    os.makedirs(DATA, exist_ok=True)
    path = os.path.join(DATA, name)
    with open(path, 'w', encoding='utf-8') as fh:
        fh.write(header + '\n')
        for row in rows:
            fh.write(','.join(f'{v:.6g}' for v in row) + '\n')
    return path


def draw(taus, d_time, d_space, f, psd, slope, half):
    '''Draw the Taylor panel and the spectrum panel.'''
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    os.makedirs(FIGS, exist_ok=True)
    path = os.path.join(FIGS, 'strip_taylor.png')
    fig, axes = plt.subplots(1, 2, figsize=(10.0, 4.0))
    axes[0].loglog(taus * 1e3, d_time, 'o-', label='temporal D(tau)')
    axes[0].loglog(taus * 1e3, d_space, 's--', label='spatial D(v tau)')
    axes[0].set_xlabel('lag tau [ms]')
    axes[0].set_ylabel('D [rad^2]')
    axes[0].set_title('Taylor frozen flow')
    axes[0].legend()
    axes[0].grid(alpha=0.3)

    keep = f > 0
    axes[1].loglog(f[keep], psd[keep], lw=0.8, label='measured')
    ref = (f[keep] / FIT_BAND_HZ[0]) ** (-8.0 / 3.0)
    scale = np.interp(FIT_BAND_HZ[0], f[keep], psd[keep])
    axes[1].loglog(f[keep], ref * scale, 'k--', label='f^(-8/3)')
    axes[1].axvspan(*FIT_BAND_HZ, color='0.85', zorder=0)
    axes[1].set_xlabel('frequency [Hz]')
    axes[1].set_ylabel('phase PSD [rad^2 / Hz]')
    axes[1].set_title(f'fitted exponent {slope:.3f} +/- {half:.3f}')
    axes[1].legend()
    axes[1].grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)
    return path


def main():
    '''Run the three measurements and print the pass bands.'''
    t_start = time.time()
    sp, paths = build()
    v = float(sp.v_x_m_s[0])
    print(f'gate (b): frozen flow on one strip, n = {N}, dx = '
          f'{DX * 1e2:.0f} cm, r0 = {R0 * 1e2:.0f} cm, L0 = {L0:.0f} m')
    print(f'  speed {v:.2f} m/s, dt {DT * 1e3:.2f} ms, {N_FRAMES} frames, '
          f'strip {sp.shape[0]} px')
    print(f'  travel {v * (N_FRAMES - 1) * DT:.1f} m, seam pad '
          f'{sp.pad * DX:.1f} m')
    print('')

    strips = open_strips(paths)
    series = sample(sp, strips)

    # ---- 1. Taylor ----
    print('1. the temporal D(tau) against the spatial D(v tau)')
    rows, d_time, d_space, taus = [], [], [], []
    for lag in LAGS:
        lag_px = int(np.rint(v * lag * DT / DX))
        dt_val, dt_se = temporal_dphi(series, int(lag))
        ds_val, ds_se = spatial_dphi(strips[0], lag_px)
        d_time.append(dt_val)
        d_space.append(ds_val)
        taus.append(lag * DT)
        tol = 2.0 * np.hypot(dt_se, ds_se) / ds_val
        rows.append((lag * DT, lag_px * DX, dt_val, dt_se, ds_val, ds_se))
        band(f'D(tau) / D(v tau) at tau = {lag * DT * 1e3:6.1f} ms',
             dt_val / ds_val, 1.0 - tol, 1.0 + tol)
    taus = np.array(taus)
    d_time = np.array(d_time)
    d_space = np.array(d_space)

    # ---- 2. the spectrum ----
    print('')
    print('2. the temporal power spectrum')
    f, pxx = welch(series.astype(np.float64), fs=1.0 / DT, axis=0,
                   nperseg=1024, detrend='constant')
    psd = pxx.mean(axis=1)
    slope, half, n_fit = fit_slope(f, psd, *FIT_BAND_HZ)
    print(f'  fitted exponent over {FIT_BAND_HZ[0]:.0f} to '
          f'{FIT_BAND_HZ[1]:.0f} Hz ({n_fit} bins): {slope:.3f} +/- {half:.3f}'
          f' (95 percent)')
    print(f'  the law is -8/3 = {-8.0 / 3.0:.3f}')
    band('temporal PSD exponent', slope, -8.0 / 3.0 - max(half, 0.10),
         -8.0 / 3.0 + max(half, 0.10))

    # ---- 3. the seam ----
    print('')
    print('3. the seam')
    ny, nx = sp.shape[0]
    travel = frame_offsets(sp, N_FRAMES - 1)[0][1]
    print(f'  the last window sits at column {travel}, of a strip of {nx}, '
          f'and the pad is {sp.pad} px')
    band('the last window ends before the pad',
         float(travel + sp.n), 0.0, float(nx - sp.pad))
    pairs = []
    for start in range(0, nx - travel - sp.n, sp.n):
        a = np.asarray(strips[0][:, start:start + sp.n], dtype=np.float64)
        b = np.asarray(strips[0][:, start + travel:start + travel + sp.n],
                       dtype=np.float64)
        pairs.append(float(np.corrcoef(a.ravel(), b.ravel())[0, 1]))
    pairs = np.array(pairs)
    se = float(pairs.std(ddof=1) / np.sqrt(pairs.size))
    print(f'  {pairs.size} disjoint pairs at the travel separation: mean '
          f'{pairs.mean():+.4f}, standard error {se:.4f}, rms '
          f'{np.sqrt((pairs ** 2).mean()):.4f}')
    band('mean corr at the travel separation', float(pairs.mean()),
         -max(2.0 * se, SEAM_MAX), max(2.0 * se, SEAM_MAX))

    # ---- the files ----
    print('')
    print('  file saved: ' + write_csv(
        'taylor_dphi.csv',
        'tau_s,r_m,d_temporal,d_temporal_se,d_spatial,d_spatial_se', rows))
    print('  file saved: ' + write_csv(
        'taylor_psd.csv', 'frequency_hz,psd_rad2_per_hz',
        list(zip(f, psd))))
    print(f'  figure saved: {draw(taus, d_time, d_space, f, psd, slope, half)}')
    print('    Caption: LEFT, the temporal structure function of one pixel of '
          'the frame')
    print('    against the spatial structure function of the same strip at the '
          'travel')
    print('    distance. RIGHT, the temporal phase power spectrum with the '
          'f^(-8/3) law;')
    print('    the grey band is the band of the fit.')

    del strips
    shutil.rmtree(STRIPS, ignore_errors=True)

    failed = [line for line in BANDS if 'FAIL' in line]
    print('')
    print(f'  {len(BANDS) - len(failed)} of {len(BANDS)} pass bands hold')
    for line in failed:
        print(line)
    print('')
    print(f'(elapsed {time.time() - t_start:.1f} s)')


if __name__ == '__main__':
    main()
