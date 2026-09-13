'''
Are the THREE strip routes the same screen? (backlog 2-P1b item 10.)

THREE ROUTES give one frame of one layer, and each one holds a RECTANGULAR
Fourier screen of a different SHORT axis:

  A. ALONG TRACK, the route of record. The wind runs along x, so the walk has
     no y part and the short axis is EXACTLY n (`strip_plan` adds no y pad
     when v_y is zero).
  B. The axis-aligned BOX of a 2-D walk. The short axis is n plus the
     perpendicular travel plus the seam pad, so a crosswind makes it many
     times n.
  C. The ROTATED thin strip. The long axis runs along the resultant velocity,
     and the short axis is the rotation footprint n(|cos|+|sin|) plus the crop
     margin on each side.

`rotated_strip_gate.py` measured C against the UNROTATED cut of the SAME crop,
which is the rotation alone. It left the REAL question open: the three routes
read 0.57 to 1.69 on the same field statistics, so ARE they the same screen?

THE HYPOTHESIS UNDER TEST. The factor is the PISTON of a window. A Fourier
screen holds no power below its own fundamental 1/(side), so a THIN strip loses
the outer-scale part of the spectrum on its short axis. That part is almost
pure PISTON over a receive aperture of 0.7 m, and a piston changes no image.
So the RAW window variance must differ between the routes, and every
PISTON-FREE quantity (the structure function, the tilt, the piston-removed
aperture variance) must agree.

THE FIVE MEASUREMENTS, route against route and against the analytic von Karman
law of Assemat and Wilson, DOI 10.1364/OE.14.000988, Eq. (5):

  1. D(r) along x and along y, at the pixel lags 1, 2 and 4 and at 3 cm to the
     aperture diameter.
  2. The Zernike tilt variance over the receive aperture, on each axis, against
     the Noll filter (`validation.screens.helpers.tilt_filter_variance`), and
     the gradient tilt against the Andrews and Phillips filter.
  3. The PISTON-REMOVED phase variance over the aperture, AND the RAW window
     variance over the whole frame. The second is the quantity that carries
     the factor.
  4. The 2-D power spectrum in radial bands, in particular the LOWEST band
     (the outer-scale content that a short axis loses) and the OUTER ring near
     Nyquist (the corner modes that a rotation loses).
  5. The PREDICTION of `captured_fraction`: the share of the tilt variance that
     a screen of the short-axis side can hold, under a sharp low-frequency
     cutoff. It is a LOWER bound, because the subharmonics put part of the band
     back.

THE WINDOWS COME FROM THE WHOLE STRIP. A spatial statistic of a frame does not
know WHERE along the strip that frame sits, so the script reads `N_WINDOWS`
windows spread over the whole long axis instead of the first few frames of a
record. Each window is built the way its route builds a frame: routes A and B
take the plain slice, and route C takes the padded crop, turns it back with
`rotate_fourier` and cuts the central n. So the strips stay the PRODUCTION
shapes, and the windows are far enough apart to be nearly independent.

THE ERROR BAR IS OVER SEEDS. Two windows of one strip 13 m apart are not
independent at a 25 m outer scale. The script averages the windows of a seed
first and it takes the standard error over the SEEDS. Every route sees the same
seed list.

THE BANDS. 3 percent on D(r) at the pixel lags, 5 percent everywhere else, each
against route B (the box), which holds the most of the spectrum. The script
also prints 2 standard errors of each difference, because the tilt is a noisy
statistic.

Outputs, next to this script:
    data/route_equivalence.csv     every route, every case, every quantity

Run from the repository root:

    python -m validation.temporal_screens.route_equivalence
    python -m validation.temporal_screens.route_equivalence --quick

Sources:
- Assemat and Wilson, Opt. Express 14(3), pp. 988 to 999 (2006),
  DOI 10.1364/OE.14.000988, Eq. (5). The closed-form von Karman phase
  covariance, which gives D(r) = 2 [B(0) - B(r)].
- Schmidt, Numerical Simulation of Optical Wave Propagation with Examples in
  MATLAB, DOI 10.1117/3.866274, Ch. 3, Eq. (3.16), printed p. 48 (the
  covariance and the structure function); Ch. 9, Eqs. (9.78) to (9.81),
  printed pp. 166 to 169 (the Fourier screen and its subharmonics).
- Noll, J. Opt. Soc. Am. 66(3), pp. 207 to 211 (1976),
  DOI 10.1364/JOSA.66.000207. Eq. (8), p. 208, the Zernike tilt filter, and
  Table IV, p. 209, the residual after piston.
- Andrews and Phillips, Laser Beam Propagation through Random Media, 2nd ed.,
  DOI 10.1117/3.626196, Ch. 6, Eqs. (80) to (84), printed pp. 200 and 201. The
  gradient tilt.
- Taylor, DOI 10.1098/rspa.1938.0032. The frozen flow of the moving crop.
- Unser, Thevenaz and Yaroslavsky, IEEE Trans. Image Process. 4, 1371 (1995),
  DOI 10.1109/83.469963. The three-shear rotation of route C.
'''

import argparse
import dataclasses
import os
import shutil
import tempfile
import time
import types

import numpy as np

from olb.turbulence.profiles import v_wind
from olb.waveoptics.schmidt.turbulence import von_karman_phase_psd
from olb.waveoptics.turbulence.temporal import (TemporalSpec, build_strips,
                                                open_strips, rotate_fourier,
                                                strip_paths, strip_plan)
from validation.screens.helpers import (captured_fraction, gradient_tilt,
                                        gtilt_filter_variance, pupil_mask,
                                        tilt_filter_variance,
                                        vk_covariance_closed, zernike_tilt)
from validation.temporal_screens.rect_factory import (DX, L0, R0, dphi_axes,
                                                      theory_dphi)

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, 'data')

# ---- the case ----
N = 256                      # the frame side, in pixels (2.56 m at DX).
H_M = 10000.0                # the layer altitude, in m.
VG = 10.0                    # the ground wind of the Bufton profile, in m/s.
# The record is NOMINAL. The windows come from the whole strip, so the travel
# only has to be non-zero. The seam pad then sets the strip length.
TRAVEL_M = 0.05
N_FRAMES = 2
PAD_OS = 2.0                 # the PRODUCTION seam pad, in outer scales.
ROT_MARGIN = 0.5             # the crop margin of route C, in units of n.
SEEDS = tuple(20260913 + i for i in range(64))
N_WINDOWS = 8                # windows along each axis of a strip.
THETA_CASES = (('hero fast layer', 6.5), ('45 deg', 45.0))
LAM = 1550e-9
PUPIL_M = 0.7                # the hero ground aperture.
LAGS_PX = (1, 2, 4)
R_M = (0.03, 0.07, 0.175, 0.35, 0.70)     # 3 cm to the aperture diameter.
BAND_SMALL = 0.03            # the band of D(r) at the pixel lags.
BAND = 0.05                  # the band everywhere else.
# The radial bands of measurement 4, as fractions of the Nyquist radius.
PSD_BANDS = ((0.0, 0.05), (0.05, 0.2), (0.2, 0.5), (0.5, 0.8), (0.8, 1.0))
REFERENCE = 'box'            # the route that every ratio divides by.

BANDS = []


def band(name, value, low, high):
    '''Record one PASS or FAIL line, and print it.'''
    ok = bool(low <= value <= high)
    line = (f"  {'PASS' if ok else 'FAIL'}  {name:56s} {value:9.4f}  "
            f"band [{low:.4f}, {high:.4f}]")
    BANDS.append(line)
    print(line)
    return ok


# ---------------------------------------------------------------------------
# The three routes
# ---------------------------------------------------------------------------

def route_plans(theta_deg, n_frames=N_FRAMES, travel_m=TRAVEL_M):
    '''Give the spec and the StripPlan of the three routes at one angle.

    The wind blows along y (`wind_dir_deg = 90`), so the y speed is the Bufton
    wind at the layer and the x speed is the slew. The slew rate that puts the
    resultant at `theta_deg` follows from tan(theta) = v_y / v_x.

    ROUTE A carries the SAME resultant SPEED along x alone: its slew rate
    carries the whole speed less the along-track wind, so the three routes hold
    strips of the same LENGTH in metres of travel and they differ in the SHORT
    axis only.

    Args:
        theta_deg: the angle of the resultant velocity, in deg, from the x
                   axis. It must be inside (0, 90).
        n_frames:  the frames of the record.
        travel_m:  the travel of the record along the resultant, in m.

    Returns:
        A tuple (routes, speed_m_s, dt_s). `routes` maps the route name to the
        pair (TemporalSpec, StripPlan), in the order A, B, C.
    '''
    v_y = float(v_wind(H_M, ws=0.0, Vg=VG))
    omega = v_y / (H_M * np.tan(np.deg2rad(theta_deg)))
    speed = v_y / np.sin(np.deg2rad(theta_deg))
    dt = travel_m / ((int(n_frames) - 1) * speed)
    plan = types.SimpleNamespace(z_m=np.array([0.0]), z_total_m=H_M,
                                 direction='down', r0_m=np.array([R0]))
    grid = types.SimpleNamespace(n=N, pixel_m=DX)
    geom = types.SimpleNamespace(elevation_deg=90.0,
                                 slew_deg_s=np.rad2deg(omega))
    # ROUTE A: the whole resultant speed along x, so v_y is exactly zero and
    # `strip_plan` gives a short axis of exactly n.
    geom_trk = types.SimpleNamespace(
        elevation_deg=90.0, slew_deg_s=np.rad2deg((speed - v_y) / H_M))
    common = dict(dt_s=dt, n_frames=int(n_frames), strip_dir='unused',
                  wind_ground_m_s=VG, pad_outer_scales=PAD_OS)
    trk = TemporalSpec(wind_dir_deg=0.0, **common)
    box = TemporalSpec(wind_dir_deg=90.0, **common)
    rot = TemporalSpec(wind_dir_deg=90.0, rotated=True,
                       rot_margin=ROT_MARGIN, **common)
    routes = {
        'along track': (trk, strip_plan(plan, grid, trk, geom_trk, L0)),
        'box': (box, strip_plan(plan, grid, box, geom, L0)),
        'rotated': (rot, strip_plan(plan, grid, rot, geom, L0)),
    }
    return routes, float(speed), float(dt)


def windows_of(spec, sp, seed, root, record, n_windows):
    '''Build one strip and give `n_windows` frames spread over its long axis.

    Each window is built the way the route builds a frame. Routes A and B take
    the plain n by n slice, which is what `frame_stack` yields for them. Route
    C takes the padded crop of side `m_crop`, turns it back by the layer angle
    with the three-shear rotation (Unser, Thevenaz and Yaroslavsky,
    DOI 10.1109/83.469963) and cuts the central n by n, which is the body of
    `frame_stack` for that route.

    THE OFFSET IS NOT THE FRAME TIME. A spatial statistic of a window does not
    know where that window sits, and the windows of the first few frames of a
    record overlap. So the offsets spread over the WHOLE long axis, which
    includes the seam pad.

    Args:
        spec:      the TemporalSpec. It names the strip file.
        sp:        the StripPlan.
        seed:      the integer seed of the layer.
        root:      the directory of the strip files.
        record:    the record index. It names the strip file, so two routes of
                   one seed must not share it.
        n_windows: the windows along each axis. A WIDE strip gives a grid of
                   n_windows by n_windows, and a THIN one gives a single row.

    Returns:
        A list of N by N float arrays.
    '''
    spec = dataclasses.replace(spec, strip_dir=root, record=int(record))
    paths = strip_paths(spec, 1)
    build_strips(sp, [R0], paths, L0, True, [int(seed)],
                 table_dtype=np.float32)
    strips = open_strips(paths)
    ny, nx = sp.shape[0]
    side = sp.n if sp.m_crop is None else int(sp.m_crop[0])
    low = (side - sp.n) // 2
    # A WIDE strip (the box) gives windows on BOTH axes, so one build gives
    # more samples. A THIN strip gives one row of windows.
    oxs = np.linspace(0, nx - side, int(n_windows)).astype(int)
    oys = (np.linspace(0, ny - side, int(n_windows)).astype(int)
           if ny - side >= side else np.array([(ny - side) // 2]))
    out = []
    for oy in oys:
        for ox in oxs:
            patch = np.array(strips[0][oy:oy + side, ox:ox + side],
                             dtype=float)
            if sp.m_crop is None:
                out.append(patch)
            else:
                out.append(np.asarray(rotate_fourier(
                    patch, -float(sp.theta[0])))[low:low + sp.n,
                                                 low:low + sp.n])
    del strips
    for path in paths:
        os.remove(path)
    return out


# ---------------------------------------------------------------------------
# The measurements
# ---------------------------------------------------------------------------

def psd_band_power(frame, band_masks, taper, plane):
    '''Give the mean 2-D power spectrum inside each radial band.

    THE PISTON AND THE TILT GO FIRST, AND THE WINDOW IS TAPERED. A raw square
    window of a turbulent phase is not periodic, so its transform carries a
    LEAKAGE tail that falls only as f^(-2) and that scales with the total
    window variance. A thin strip and a wide box hold different amounts of
    that variance, so the raw band power would read the leakage of the piston
    and the tilt, not the turbulence of the band. The least-squares plane
    removal takes the two largest terms away, and the separable Hann taper
    (Schmidt, DOI 10.1117/3.866274, Ch. 3) makes the rest of the edge smooth.

    Args:
        frame:       one n by n phase frame.
        band_masks:  one boolean mask for each band, on the shifted k grid.
        taper:       the separable Hann window, of the shape of the frame.
        plane:       the design matrix (1, x, y) of the whole frame.

    Returns:
        A float array, one mean power for each band.
    '''
    flat = frame.ravel()
    coeff = np.linalg.lstsq(plane, flat, rcond=None)[0]
    residual = (flat - plane @ coeff).reshape(frame.shape) * taper
    power = np.abs(np.fft.fftshift(np.fft.fft2(residual))) ** 2
    return np.array([float(np.mean(power[m])) for m in band_masks])


def piston_free_variance(frame, mask):
    '''Give the phase variance over the aperture, with the piston removed.

    The piston of a window is the mean of the phase over the aperture. It
    delays every ray by the same amount, so it moves no image and it couples
    no differently.

    Args:
        frame: one n by n phase frame, in rad.
        mask:  the boolean aperture mask.

    Returns:
        The variance, in rad^2.
    '''
    values = frame[mask]
    return float(np.mean((values - values.mean()) ** 2))


def piston_free_law(mask, dx_m, r0_m, l0_m):
    '''Give the analytic piston-removed phase variance over an aperture.

    formula:
        Var(phi - mean(phi)) = (1 / Np^2) SUM_i SUM_j D(|r_i - r_j|) / 2
    Source: the covariance identity D(r) = 2 [B(0) - B(r)] of Schmidt,
    DOI 10.1117/3.866274, Ch. 3, Eq. (3.16), printed p. 48, applied to
    Var = B(0) - (1/Np^2) SUM SUM B(r_ij). The covariance B is the closed form
    of Assemat and Wilson, DOI 10.1364/OE.14.000988, Eq. (5).

    THE PAIR COUNT comes from the autocorrelation of the mask, which is one
    transform. The mask is padded to twice its side, so the correlation does
    not wrap.

    Args:
        mask:  the boolean aperture mask.
        dx_m:  the pixel pitch, in m.
        r0_m:  the Fried parameter, in m.
        l0_m:  the outer scale, in m.

    Returns:
        The variance, in rad^2.
    '''
    m = np.asarray(mask, dtype=float)
    n = m.shape[0]
    padded = np.zeros((2 * n, 2 * n))
    padded[:n, :n] = m
    count = np.fft.irfft2(np.abs(np.fft.rfft2(padded)) ** 2, s=(2 * n, 2 * n))
    lag = np.fft.fftfreq(2 * n, d=1.0 / (2 * n))
    lx, ly = np.meshgrid(lag, lag)
    r = np.hypot(lx, ly) * float(dx_m)
    b0 = float(vk_covariance_closed(0.0, r0_m, l0_m)[0])
    d_of_r = 2.0 * (b0 - vk_covariance_closed(r.ravel(),
                                              r0_m, l0_m).reshape(r.shape))
    weight = count / count.sum()
    return float(np.sum(weight * d_of_r / 2.0))


def frame_row(frame, mask, ks, band_masks, taper, plane):
    '''Give every scalar of ONE frame, as a flat dict.

    Args:
        frame:      one n by n phase frame, in rad.
        mask:       the boolean aperture mask.
        ks:         the pixel lags of the structure function.
        band_masks: the radial band masks of the power spectrum.
        taper:      the Hann window of `psd_band_power`.
        plane:      the design matrix of `psd_band_power`.

    Returns:
        A dict of the name of each quantity against its value.
    '''
    dy, dx = dphi_axes(frame, ks)
    _, _, a2, a3 = zernike_tilt(frame, mask, DX, LAM)
    gx, gy = gradient_tilt(frame, mask, DX, LAM)
    row = {'var_raw': float(np.var(frame)),
           'var_piston_free': piston_free_variance(frame, mask),
           'ztilt_x': a2 ** 2, 'ztilt_y': a3 ** 2,
           'gtilt_x': gx ** 2, 'gtilt_y': gy ** 2}
    for i, k in enumerate(ks):
        row[f'dphi_y_{int(k)}'] = float(dy[i])
        row[f'dphi_x_{int(k)}'] = float(dx[i])
    for i, power in enumerate(psd_band_power(frame, band_masks, taper,
                                             plane)):
        row[f'psd_{i}'] = float(power)
    return row


def measure_route(spec, sp, seeds, root, record0, n_windows, mask, ks,
                  band_masks, taper, plane):
    '''Give the mean and the standard error of every quantity of one route.

    The windows of ONE seed are averaged FIRST, because they are windows of
    one strip and the outer scale is 25 m, so they are not independent. The
    standard error then runs over the SEEDS.

    Args:
        spec:       the TemporalSpec of the route.
        sp:         the StripPlan of the route.
        seeds:      the integer seed of each record.
        root:       the directory of the strip files.
        record0:    the first record index of this route.
        n_windows:  the windows of each strip.
        mask:       the aperture mask.
        ks:         the pixel lags.
        band_masks: the radial band masks.
        taper:      the Hann window of `psd_band_power`.
        plane:      the design matrix of `psd_band_power`.

    Returns:
        A tuple (mean, standard_error), two dicts of the same keys.
    '''
    per_seed = []
    for s, seed in enumerate(seeds):
        rows = [frame_row(f, mask, ks, band_masks, taper, plane)
                for f in windows_of(spec, sp, seed, root, record0 + s,
                                    n_windows)]
        per_seed.append({k: float(np.mean([r[k] for r in rows]))
                         for k in rows[0]})
    keys = list(per_seed[0])
    mean, err = {}, {}
    for k in keys:
        values = np.array([r[k] for r in per_seed], dtype=float)
        mean[k] = float(values.mean())
        err[k] = float(values.std(ddof=1) / np.sqrt(values.size))
    return mean, err


# ---------------------------------------------------------------------------
# The report
# ---------------------------------------------------------------------------

def ratio_line(name, got, got_se, ref, ref_se, tol):
    '''Print one route ratio with its 2 SE bar, and record the PASS band.'''
    ratio = got / ref
    two_se = 2.0 * ratio * np.hypot(got_se / got if got else 0.0,
                                    ref_se / ref if ref else 0.0)
    print(f'      {name:34s}{got:12.5g}{ratio:10.4f}  +/- {two_se:7.4f}')
    return band(name, ratio, 1.0 - tol, 1.0 + tol), two_se


def run_case(case_name, theta_deg, args, mask, ks, band_masks, taper,
             plane, rows):
    '''Run the three routes of one case and print every table.'''
    routes, speed, dt = route_plans(theta_deg, args.frames, args.travel)
    print('')
    print(f'CASE {case_name}, resultant {theta_deg:.1f} deg, '
          f'|v| = {speed:.1f} m/s, dt = {dt * 1e3:.3f} ms')
    print(f'  {"route":14s}{"strip [px]":>22s}{"short axis [m]":>16s}'
          f'{"Mpx":>10s}')
    for j, (name, (_, sp)) in enumerate(routes.items()):
        ny, nx = sp.shape[0]
        short = (sp.m_crop[0] if sp.m_crop is not None else ny) * DX
        print(f'  {name:14s}{str((ny, nx)):>22s}{short:16.2f}'
              f'{ny * nx / 1e6:10.2f}')

    got = {}
    root = tempfile.mkdtemp(prefix='olb_routeeq_')
    try:
        for j, (name, (spec, sp)) in enumerate(routes.items()):
            t0 = time.time()
            got[name] = measure_route(spec, sp, args.seeds, root,
                                      1000 * (j + 1), args.windows, mask, ks,
                                      band_masks, taper, plane)
            print(f'  route {name:14s} done in {time.time() - t0:6.1f} s')
    finally:
        shutil.rmtree(root, ignore_errors=True)

    ref_mean, ref_se = got[REFERENCE]
    others = [n for n in got if n != REFERENCE]
    theory = theory_dphi(np.asarray(ks) * DX)
    psd = lambda f: von_karman_phase_psd(f, R0, L0)
    law_ztilt = tilt_filter_variance(psd, PUPIL_M)
    law_gtilt = gtilt_filter_variance(psd, PUPIL_M, LAM)
    law_pf = piston_free_law(mask, DX, R0, L0)

    print('')
    print('  1. the structure function D(r), against the box and the law')
    print(f'      {"quantity":34s}{"value":>12s}{"/ box":>10s}     '
          f'{"2 SE":>7s}')
    for i, k in enumerate(ks):
        r_m = k * DX
        tol = BAND_SMALL if k in LAGS_PX else BAND
        for axis in ('y', 'x'):
            key = f'dphi_{axis}_{int(k)}'
            print(f'      {"box   D_" + axis + f"({r_m:.3f} m) / law":34s}'
                  f'{ref_mean[key]:12.5g}'
                  f'{ref_mean[key] / theory[i]:10.4f}')
            for name in others:
                ok, _ = ratio_line(f'{case_name[:12]} {name} D_{axis}'
                                   f'({r_m:.3f} m)',
                                   got[name][0][key], got[name][1][key],
                                   ref_mean[key], ref_se[key], tol)
                rows.append((case_name, name, f'dphi_{axis}_{r_m:.3f}',
                             got[name][0][key], got[name][1][key],
                             ref_mean[key], theory[i]))

    print('')
    print('  2. the tilt over the 0.70 m aperture')
    print(f'      the Noll law <a^2> = {law_ztilt:.4f} rad^2, the G-tilt law '
          f'<alpha^2> = {law_gtilt:.4e} rad^2')
    for key, law, label in (('ztilt_x', law_ztilt, 'Z-tilt x'),
                            ('ztilt_y', law_ztilt, 'Z-tilt y'),
                            ('gtilt_x', law_gtilt, 'G-tilt x'),
                            ('gtilt_y', law_gtilt, 'G-tilt y')):
        print(f'      {"box   " + label + " / law":34s}{ref_mean[key]:12.5g}'
              f'{ref_mean[key] / law:10.4f}')
        for name in others:
            ratio_line(f'{case_name[:12]} {name} {label}',
                       got[name][0][key], got[name][1][key],
                       ref_mean[key], ref_se[key], BAND)
            rows.append((case_name, name, label, got[name][0][key],
                         got[name][1][key], ref_mean[key], law))

    print('')
    print('  3. the piston-removed aperture variance and the RAW window '
          'variance')
    print(f'      the piston-free law over the 0.70 m aperture = '
          f'{law_pf:.4f} rad^2')
    for key, law, label in (('var_piston_free', law_pf, 'piston-free var'),
                            ('var_raw', np.nan, 'raw window var')):
        print(f'      {"box   " + label:34s}{ref_mean[key]:12.5g}'
              f'{ref_mean[key] / law:10.4f}')
        for name in others:
            ratio_line(f'{case_name[:12]} {name} {label}',
                       got[name][0][key], got[name][1][key],
                       ref_mean[key], ref_se[key], BAND)
            rows.append((case_name, name, label, got[name][0][key],
                         got[name][1][key], ref_mean[key], law))

    print('')
    print('  4. the 2-D power spectrum in radial bands (of the Nyquist '
          'radius)')
    for i, (lo, hi) in enumerate(PSD_BANDS):
        key = f'psd_{i}'
        label = f'PSD {lo:.2f} to {hi:.2f}'
        print(f'      {"box   " + label:34s}{ref_mean[key]:12.5g}')
        for name in others:
            ratio_line(f'{case_name[:12]} {name} {label}',
                       got[name][0][key], got[name][1][key],
                       ref_mean[key], ref_se[key], BAND)
            rows.append((case_name, name, label, got[name][0][key],
                         got[name][1][key], ref_mean[key], np.nan))

    print('')
    print('  5. the sharp-cutoff PREDICTION of the tilt deficit of a short '
          'axis')
    print(f'      {"route":14s}{"short axis [m]":>16s}{"captured":>12s}'
          f'{"measured Z-tilt / box":>24s}')
    for name, (spec, sp) in routes.items():
        ny = sp.shape[0][0]
        side = (sp.m_crop[0] if sp.m_crop is not None else ny) * DX
        share = captured_fraction(psd, PUPIL_M, 1.0 / side)
        measured = ((got[name][0]['ztilt_x'] + got[name][0]['ztilt_y'])
                    / (ref_mean['ztilt_x'] + ref_mean['ztilt_y']))
        print(f'      {name:14s}{side:16.2f}{share:12.4f}{measured:24.4f}')
        rows.append((case_name, name, 'captured_fraction', share, 0.0,
                     np.nan, np.nan))
    return got


def write_csv(rows):
    '''Write every measured quantity. Give the path back.'''
    os.makedirs(DATA, exist_ok=True)
    path = os.path.join(DATA, 'route_equivalence.csv')
    with open(path, 'w', encoding='utf-8') as fh:
        fh.write('case,route,quantity,value,standard_error,box,law\n')
        for r in rows:
            fh.write(f'{r[0]},{r[1]},{r[2]},{r[3]:.8g},{r[4]:.8g},'
                     f'{r[5]:.8g},{r[6]:.8g}\n')
    return path


def _self_check():
    '''Check the two new estimators against hand cases. It runs every time.'''
    # A pure PISTON has no piston-free variance and a large raw variance.
    mask = pupil_mask(32, 0.1, 2.0)
    assert piston_free_variance(np.full((32, 32), 3.0), mask) < 1e-24
    # A ramp phi = a x over a mask: the piston-free variance is a^2 <x^2>.
    axis = (np.arange(32) - 16) * 0.1
    xx = np.meshgrid(axis, axis)[0]
    got = piston_free_variance(2.0 * xx, mask)
    want = 4.0 * float(np.var(xx[mask]))
    assert abs(got / want - 1.0) < 1e-12, (got, want)
    # The piston-free LAW of a tiny aperture goes to zero, and it grows with
    # the aperture.
    small = piston_free_law(pupil_mask(64, 0.01, 0.05), 0.01, R0, L0)
    large = piston_free_law(pupil_mask(64, 0.01, 0.60), 0.01, R0, L0)
    assert 0.0 < small < large, (small, large)
    # The law of a 0.6 m aperture must sit near the Noll piston residual
    # 1.0299 (D/r0)^(5/3) (Noll, DOI 10.1364/JOSA.66.000207, Table IV, p. 209).
    # The von Karman outer scale makes the measured value LOWER, never higher.
    noll = 1.0299 * (0.60 / R0) ** (5.0 / 3.0)
    assert 0.3 * noll < large < noll, (large, noll)


def main():
    '''Run the three routes of both cases and print every band.'''
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument('--quick', action='store_true',
                        help='three seeds and one case, for a smoke run.')
    parser.add_argument('--seeds', type=int, default=len(SEEDS),
                        help='the records of each route.')
    parser.add_argument('--frames', type=int, default=N_FRAMES,
                        help='the frames of one record.')
    parser.add_argument('--travel', type=float, default=TRAVEL_M,
                        help='the travel of one record [m].')
    parser.add_argument('--windows', type=int, default=N_WINDOWS,
                        help='the windows of each strip.')
    args = parser.parse_args()
    cases = THETA_CASES
    if args.quick:
        args.seeds = 4
        cases = THETA_CASES[1:]
    args.seeds = list(SEEDS[:int(args.seeds)])

    _self_check()
    t_start = time.time()
    print(f'the three strip routes, n = {N}, dx = {DX * 1e2:.0f} cm, '
          f'r0 = {R0 * 1e2:.0f} cm, L0 = {L0:.0f} m, pupil '
          f'{PUPIL_M * 1e2:.0f} cm')
    print(f'  {len(args.seeds)} seeds x {args.windows} windows, seam pad '
          f'{PAD_OS:.1f} L0, rotation margin {ROT_MARGIN:.2f} n')
    print(f'  every ratio divides by the {REFERENCE!r} route')

    ks = np.array(list(LAGS_PX) + [int(round(r / DX)) for r in R_M])
    mask = pupil_mask(N, DX, PUPIL_M)
    f_axis = np.fft.fftshift(np.fft.fftfreq(N))
    radius = np.hypot(f_axis[None, :], f_axis[:, None]) / f_axis.max()
    band_masks = [(radius >= lo) & (radius <= hi) for lo, hi in PSD_BANDS]
    hann = np.hanning(N)
    taper = hann[:, None] * hann[None, :]
    grid_axis = (np.arange(N) - N // 2) * DX
    gx, gy = np.meshgrid(grid_axis, grid_axis)
    plane = np.column_stack([np.ones(N * N), gx.ravel(), gy.ravel()])

    rows = []
    for case_name, theta_deg in cases:
        run_case(case_name, theta_deg, args, mask, ks, band_masks, taper,
                 plane, rows)

    print('')
    print(f'  file saved: {write_csv(rows)}')
    failed = [line for line in BANDS if 'FAIL' in line]
    print('')
    print(f'  {len(BANDS) - len(failed)} of {len(BANDS)} pass bands hold')
    for line in failed:
        print(line)
    print('')
    print(f'(elapsed {time.time() - t_start:.1f} s)')


if __name__ == '__main__':
    main()
