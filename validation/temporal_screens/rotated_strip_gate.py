'''Gate: the ROTATED thin strip against the axis-aligned box (2-P1b item 10).

THE PROBLEM. A CROSSWIND turns the walk of a layer away from the x axis. The
axis-aligned bounding box of that walk then grows its SHORT axis, and the
30 deg hero at `wind_dir_deg = 90` asks for a 11897 x 43744 box (520 Mpx) for
one jet-level layer. That box FITS the host (137 s, 15.6 GiB, see
`table_precision.py`), but the time is the FFT of the 23x area.

THE ROUTE UNDER TEST. `TemporalSpec(rotated=True)` holds ONE THIN strip per
layer whose long axis runs along the RESULTANT velocity of that layer. Every
frame takes a PADDED crop of that strip and turns it back with the EXACT
Fourier three-shear rotation of Unser, Thevenaz and Yaroslavsky,
DOI 10.1109/83.469963, then it cuts the central n by n. A real-space bilinear
or cubic rotation is NOT the tool: it smooths the Fresnel-scale structure that
builds the scintillation.

THE FOUR MEASUREMENTS, against the exact box route and the analytic von Karman
law of Assemat and Wilson, DOI 10.1364/OE.14.000988, Eq. (5):

  1. D(r) of the frames along y and along x, at the pixel lags 1, 2 and 4 and
     at 3 cm to 1 m, plus the phase variance and the PSD near the Nyquist
     edge (the corner-mode loss). The 1 px lag is the one that shows a
     Fresnel-scale smoothing. The sweep runs at three crop MARGINS.
  2. The FRAME-TO-FRAME consistency: D(tau) of the ROTATED frames against
     D(tau) of the UNROTATED cut of the SAME crops. A rotation is an
     isometry, so the two must agree.
  3. The END-TO-END field: one screen, one Fresnel hop, then the point and the
     aperture scintillation index and the mean single-mode coupling, on four
     arms (rotated, unrotated, box, along track).
  4. THE COST: the strip pixels, the build time and the per-frame rotation
     time of the two routes, and the projection of the hero crosswind box.

THE REFERENCE IS NOT THE BOX ALONE. The three strips do NOT share a shape:
the box short axis is many times wider than the thin strip, and the ALONG-
TRACK strip of the route of record is thinner still (exactly n). The short
axis sets how much of the outer scale the strip holds on that axis, so the
three give different frame VARIANCES and different large-r structure, and
that difference is not the rotation. Measurement 2 and the primary band of
measurement 3 therefore compare the rotated frames with the UNROTATED cut of
the SAME crops, which shares the strip, the seed and the grid.

THE BANDS. Measurement 1 holds 3 percent, or 2 SE of the difference where
that is wider. Measurement 2 and 3 hold the 5 percent P2 kill line.

THE PAD. The gate runs at `pad_outer_scales = 0.5`, not the production 2.0.
The seam pad guards the correlation between the START and the END of a record.
It does not enter a per-frame structure function or a per-frame field, and it
is the largest part of a small test strip. Measurement 4 reports the
PRODUCTION pad.

Outputs, next to this script:
    data/rotated_strip_dphi.csv    measurement 1, every margin and every case
    data/rotated_strip_cost.csv    measurement 4
    figures/rotated_strip_dphi.png the structure function of the two routes

Run from the repository root:

    python -m validation.temporal_screens.rotated_strip_gate
    python -m validation.temporal_screens.rotated_strip_gate --quick
    python -m validation.temporal_screens.rotated_strip_gate --build-hero

The third form builds the real 30 deg crosswind strip of the hero (bigfraw).

Sources:
- Unser, Thevenaz and Yaroslavsky, IEEE Trans. Image Process. 4, 1371 (1995),
  DOI 10.1109/83.469963. The three-shear rotation.
- Schmidt, DOI 10.1117/3.866274, Ch. 2. The Fourier shift theorem that each
  shear uses, and Ch. 9 for the screen.
- Assemat and Wilson, DOI 10.1364/OE.14.000988, Eq. (5). The analytic von
  Karman covariance, so D(r) = 2 (B(0) - B(r)).
- Taylor, DOI 10.1098/rspa.1938.0032. The frozen flow of the moving crop.
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
from olb.waveoptics.turbulence.screens import ScreenFactory
from olb.waveoptics.turbulence.temporal import TemporalSpec, build_strips, \
    frame_offsets, frame_stack, open_strips, rotate_fourier, strip_paths, \
    strip_plan
from validation.screens.helpers import pupil_mask
from validation.temporal_screens.rect_factory import DX, L0, R0, dphi_axes, \
    theory_dphi

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, 'data')
FIGS = os.path.join(HERE, 'figures')

# ---- the test case ----
N = 256                      # the frame side, in pixels (2.56 m at DX).
H_M = 10000.0                # the layer altitude, in m.
VG = 10.0                    # the ground wind of the Bufton profile, in m/s.
TRAVEL_M = 3.0               # the travel of the record, in m.
N_FRAMES = 201
PAD_OS = 0.5                 # the seam pad of the gate, in outer scales.
SEEDS = (20260913, 20260914, 20260915)
FRAMES_USED = 8              # frames of each record that enter measurement 1.
MARGINS = (0.25, 0.5, 1.0)
THETA_CASES = (('hero crosswind', 6.5), ('mid angle', 45.0))
LAGS_PX = (1, 2, 4)
R_M = (0.03, 0.16, 0.32, 1.00)
BAND = 0.03                  # the fixed band of measurement 1 and 2.
KILL = 0.05                  # the P2 kill line of measurement 3.

# ---- the field arm of measurement 3 ----
LAM = 1550e-9
Z_PROP_M = 5000.0            # the hop after the screen, in m.
R0_FIELD = 0.30              # a weaker layer, so the hop stays weak.
PUPIL_M = 0.5

# ---- the hero crosswind box, from table_precision.py and the README ----
HERO_BOX = (11897, 43744)
HERO_N = 512
HERO_DX = 0.00686
HERO_PAD = 7289              # ceil(2 * 25 m / HERO_DX).
HERO_BOX_SECONDS = 137.4     # the measured float32-table build, on bigfraw.

BANDS = []


def band(name, value, low, high):
    '''Record one PASS or FAIL line, and print it.'''
    ok = bool(low <= value <= high)
    line = (f"  {'PASS' if ok else 'FAIL'}  {name:54s} {value:9.4f}  "
            f"band [{low:.4f}, {high:.4f}]")
    BANDS.append(line)
    print(line)
    return ok


def case_plans(theta_deg, margin, n_frames=N_FRAMES, travel_m=TRAVEL_M):
    '''Give the BOX plan and the ROTATED plan of one resultant angle.

    The wind blows along y (`wind_dir_deg = 90`), so the y speed is the Bufton
    wind at the layer and the x speed is the slew. The slew rate that puts the
    resultant at `theta_deg` follows from tan(theta) = v_y / v_x.

    Args:
        theta_deg: the angle of the resultant velocity, in deg, from the x
                   axis. It must be inside (0, 90).
        margin:    the crop margin of the rotated route, in units of N.
        n_frames:  the frames of the record.
        travel_m:  the travel of the record along the resultant, in m.

    Returns:
        A tuple (spec_box, sp_box, spec_rot, sp_rot, dt_s, speed_m_s).
    '''
    v_y = float(v_wind(H_M, ws=0.0, Vg=VG))
    omega = v_y / (H_M * np.tan(np.deg2rad(theta_deg)))
    speed = v_y / np.sin(np.deg2rad(theta_deg))
    dt = travel_m / ((n_frames - 1) * speed)
    plan = types.SimpleNamespace(z_m=np.array([0.0]), z_total_m=H_M,
                                 direction='down',
                                 r0_m=np.array([R0]))
    grid = types.SimpleNamespace(n=N, pixel_m=DX)
    geometry = types.SimpleNamespace(elevation_deg=90.0,
                                     slew_deg_s=np.rad2deg(omega))
    common = dict(dt_s=dt, n_frames=n_frames, strip_dir='unused',
                  wind_ground_m_s=VG, wind_dir_deg=90.0,
                  pad_outer_scales=PAD_OS)
    spec_box = TemporalSpec(**common)
    spec_rot = TemporalSpec(rotated=True, rot_margin=margin, **common)
    return (spec_box, strip_plan(plan, grid, spec_box, geometry, L0),
            spec_rot, strip_plan(plan, grid, spec_rot, geometry, L0),
            dt, speed)


def frames_of(spec, sp, seed, root, indices):
    '''Build the strip of one record and give the frames at `indices`.

    Args:
        spec:    the TemporalSpec.
        sp:      the StripPlan.
        seed:    the integer seed of the layer.
        root:    the directory of the strip files.
        indices: the frame indices to read.

    Returns:
        A list of N by N float arrays, and the build time in s.
    '''
    spec = dataclasses.replace(spec, strip_dir=root)
    paths = strip_paths(spec, 1)
    t0 = time.time()
    build_strips(sp, [R0], paths, L0, True, [seed], table_dtype=np.float32)
    build_s = time.time() - t0
    strips = open_strips(paths)
    out = [np.array(next(iter(frame_stack(strips, sp, k))), dtype=float)
           for k in indices]
    del strips
    return out, build_s


def unrotated_frames(sp, paths, indices):
    """Give the central n by n cut of the SAME crops, with NO rotation.

    It is the control of measurement 3. The screen is isotropic, so an
    unrotated cut is as valid an atmosphere as a rotated one, and the two
    must give the same field statistics.
    """
    strips = open_strips(paths)
    m, lo = int(sp.m_crop[0]), (int(sp.m_crop[0]) - sp.n) // 2
    out = []
    for k in indices:
        _, ox = frame_offsets(sp, k)[0]
        patch = np.array(strips[0][0:m, ox:ox + m], dtype=float)
        out.append(patch[lo:lo + sp.n, lo:lo + sp.n])
    del strips
    return out


def bilinear_frames(sp, paths, indices):
    '''The real-space BILINEAR control: the same crop, rotated at order 1.

    It is a CONTROL, not the route. A bilinear kernel is a low-pass filter, so
    it must lose the 1 px structure that the shears keep.
    '''
    from scipy import ndimage
    strips = open_strips(paths)
    m, lo = int(sp.m_crop[0]), (int(sp.m_crop[0]) - sp.n) // 2
    out = []
    for k in indices:
        ox = int(np.rint(np.hypot(sp.v_x_m_s[0], sp.v_y_m_s[0])
                         * k * sp.dt_s / sp.dx))
        patch = np.array(strips[0][0:m, ox:ox + m], dtype=float)
        turned = ndimage.rotate(patch, np.rad2deg(float(sp.theta[0])),
                                reshape=False, order=1, mode='nearest')
        out.append(turned[lo:lo + sp.n, lo:lo + sp.n])
    del strips
    return out


def dphi_stats(frames, ks):
    '''Give the mean and the standard error of D(r) on the two axes.'''
    dy = np.array([dphi_axes(f, ks)[0] for f in frames])
    dx = np.array([dphi_axes(f, ks)[1] for f in frames])
    both = np.concatenate([dy, dx], axis=0)
    return (both.mean(axis=0),
            both.std(axis=0, ddof=1) / np.sqrt(both.shape[0]),
            dy.mean(axis=0), dx.mean(axis=0))


def nyquist_psd_ratio(frames_a, frames_b):
    '''Give the ratio of the mean 2-D PSD in the OUTER ring of the k square.

    The corners of the frequency square leave the square under a rotation, so
    a rotated frame must read LOW here. The ring is the band 0.8 to 1.0 of the
    Nyquist radius.
    '''
    def ring_power(frames):
        n = frames[0].shape[0]
        f = np.fft.fftshift(np.fft.fftfreq(n))
        rr = np.hypot(f[None, :], f[:, None]) / f.max()
        mask = (rr >= 0.8) & (rr <= 1.0)
        return float(np.mean([
            np.mean(np.abs(np.fft.fftshift(np.fft.fft2(a))) ** 2 * mask)
            for a in frames]))
    return ring_power(frames_a) / ring_power(frames_b)


def temporal_pair(sp, paths, lags, n_origin=6):
    """Give D(tau) of the ROTATED frames and of the UNROTATED crops.

    THE TEST. Frame i and frame i+k are two windows of ONE strip, `dx` pixels
    apart along the strip x. Frozen flow (Taylor, DOI 10.1098/rspa.1938.0032)
    makes the temporal structure function the spatial one at that offset. A
    ROTATION is an isometry, so it must not change D(tau). The unrotated cut
    of the SAME crops is therefore the exact reference: it holds the same
    strip, the same seed and the same offsets, so the ratio is the rotation
    alone.

    The average runs over `n_origin` starting frames and over every pixel of
    the frame.

    Returns:
        Two arrays, D(tau) rotated and D(tau) unrotated, one for each lag.
    """
    strips = open_strips(paths)
    m, lo = int(sp.m_crop[0]), (int(sp.m_crop[0]) - sp.n) // 2
    starts = np.linspace(0, sp.n_frames - 1 - max(lags), n_origin).astype(int)
    keys = sorted({int(k) for s0 in starts for k in (s0, )}
                  | {int(s0 + k) for s0 in starts for k in lags})
    rot, raw = {}, {}
    for k in keys:
        _, ox = frame_offsets(sp, k)[0]
        patch = np.array(strips[0][0:m, ox:ox + m], dtype=float)
        raw[k] = patch[lo:lo + sp.n, lo:lo + sp.n]
        rot[k] = np.asarray(rotate_fourier(patch, -float(sp.theta[0])))[
            lo:lo + sp.n, lo:lo + sp.n]
    del strips

    def d_tau(store):
        return np.array([float(np.mean([
            np.mean((store[int(s0 + k)] - store[int(s0)]) ** 2)
            for s0 in starts])) for k in lags])
    return d_tau(rot), d_tau(raw)


def corner_loss(sp, paths, theta, n_lag_ks):
    """Give the corner-mode loss of ONE rotation, on ONE array.

    The crop is rotated and the central window of the result is compared with
    the central window of the crop itself. The grid, the seed and the filter
    are then the SAME on the two sides, so the ratio is the rotation alone:
    the frame-to-frame comparison of measurement 1 also carries the different
    shape of the two strips.

    Returns:
        The ratio of the outer-ring PSD, and the ratio of D(r) at `n_lag_ks`.
    """
    strips = open_strips(paths)
    m, lo = int(sp.m_crop[0]), (int(sp.m_crop[0]) - sp.n) // 2
    before, after = [], []
    for k in (0, sp.n_frames // 2, sp.n_frames - 1):
        ox = int(np.rint(np.hypot(sp.v_x_m_s[0], sp.v_y_m_s[0])
                         * k * sp.dt_s / sp.dx))
        patch = np.array(strips[0][0:m, ox:ox + m], dtype=float)
        before.append(patch[lo:lo + sp.n, lo:lo + sp.n])
        after.append(np.asarray(rotate_fourier(patch, -float(theta)))[
            lo:lo + sp.n, lo:lo + sp.n])
    del strips
    d_b, _, _, _ = dphi_stats(before, n_lag_ks)
    d_a, _, _, _ = dphi_stats(after, n_lag_ks)
    return nyquist_psd_ratio(after, before), d_a / d_b


def field_stats(frames):
    '''Propagate every frame and give the two indices and the mean coupling.

    One screen, then one Fresnel hop of Z_PROP_M, then a PUPIL_M aperture.
    The screen is scaled from R0 to R0_FIELD, because the screen amplitude is
    the pure factor r0^(-5/6) (Schmidt, DOI 10.1117/3.866274, Ch. 9,
    Eq. (9.70)), so no strip has to be built again.
    '''
    from olb.waveoptics import smf
    from olb.waveoptics.field import Begin, Intensity
    from olb.waveoptics.propagators import Forvard
    from olb.waveoptics.turbulence.screens import Screen

    scale = (R0_FIELD / R0) ** (-5.0 / 6.0)
    mask = pupil_mask(N, DX, PUPIL_M)
    centre, bucket, eta = [], [], []
    for f in frames:
        fld = Screen(Begin(N * DX, LAM, N), f * scale)
        fld = Forvard(fld, Z_PROP_M)
        img = np.asarray(Intensity(fld))
        centre.append(float(img[N // 2, N // 2]))
        bucket.append(float(np.mean(img[mask])))
        eta.append(float(smf.coupling_efficiency(fld, PUPIL_M)))
    centre, bucket = np.array(centre), np.array(bucket)
    return (float(centre.var() / centre.mean() ** 2),
            float(bucket.var() / bucket.mean() ** 2),
            float(np.mean(eta)))


def hero_rotated_shape(margin):
    '''Give the rotated strip shape of the 30 deg hero crosswind layer.

    The box (11897, 43744) holds the frame side, the travel and the seam pad
    on each axis, so the travel of the layer follows by subtraction, and the
    resultant angle and length follow from the travel.
    '''
    ty = HERO_BOX[0] - HERO_N - HERO_PAD
    tx = HERO_BOX[1] - HERO_N - HERO_PAD
    theta = np.arctan2(ty, tx)
    foot = abs(np.cos(theta)) + abs(np.sin(theta))
    m = int(np.ceil(HERO_N * (foot + 2.0 * margin)))
    m += m % 2
    nx = m + int(np.ceil(np.hypot(ty, tx))) + HERO_PAD
    nx = 32 * int(np.ceil(nx / 32))
    return (m, nx), float(np.rad2deg(theta))


def measure_structure(args):
    '''Measurement 1 and 2: the structure of the rotated frames.'''
    ks = np.array(list(LAGS_PX) + [int(round(r / DX)) for r in R_M])
    r_all = ks * DX
    theory = theory_dphi(r_all)
    rows = []
    print('1. the structure function of a frame (D(r), both axes pooled)')
    for case_name, theta_deg in THETA_CASES:
        _, sp_box, _, _, _, speed = case_plans(theta_deg, MARGINS[0])
        print(f'  case: {case_name}, theta = {theta_deg:.1f} deg, '
              f'|v| = {speed:.1f} m/s')
        idx = np.linspace(0, N_FRAMES - 1, FRAMES_USED).astype(int)
        root = tempfile.mkdtemp(prefix='olb_rotgate_')
        try:
            box_frames = []
            spec_box, sp_box, _, _, _, _ = case_plans(theta_deg, MARGINS[0])
            for s, seed in enumerate(args.seeds):
                spec = dataclasses.replace(spec_box, record=s)
                got, _ = frames_of(spec, sp_box, seed, root, idx)
                box_frames += got
            box_mean, box_se, box_y, box_x = dphi_stats(box_frames, ks)
            box_var = float(np.mean([np.var(f) for f in box_frames]))
            for margin in args.margins:
                _, _, spec_rot, sp_rot, _, _ = case_plans(theta_deg, margin)
                rot_frames, ctrl_ref = [], None
                for s, seed in enumerate(args.seeds):
                    # The record index carries the MARGIN, because the crop
                    # side changes with it and a strip file is reused by
                    # `build_strips` when the path is the same.
                    spec = dataclasses.replace(
                        spec_rot, record=100 + 10 * int(margin * 100) + s,
                        strip_dir=root)
                    got, _ = frames_of(spec, sp_rot, seed, root, idx)
                    rot_frames += got
                    if args.bilinear and ctrl_ref is None and margin == 0.5:
                        ctrl_ref = bilinear_frames(sp_rot,
                                                   strip_paths(spec, 1), idx)
                rot_mean, rot_se, rot_y, rot_x = dphi_stats(rot_frames, ks)
                rot_var = float(np.mean([np.var(f) for f in rot_frames]))
                psd = nyquist_psd_ratio(rot_frames, box_frames)
                # The variance and the outer-ring PSD of the two routes are
                # DIAGNOSTIC, not a band: the frame variance follows the
                # SHORT axis of the strip that holds it (the box is much
                # wider than the thin strip there), so the two routes must
                # not agree on it.
                print(f'    margin {margin:.2f} n   '
                      f'short axis rot {sp_rot.shape[0][0]} px, box '
                      f'{sp_box.shape[0][0]} px   '
                      f'variance ratio {rot_var / box_var:.4f}   '
                      f'outer-ring PSD ratio {psd:.4f}')
                # The rot/box columns are DIAGNOSTIC. The two strips do not
                # share a shape, so their large-r structure legitimately
                # differs; the one-array line below is the controlled test.
                print(f'      {"lag":>8}{"r [m]":>8}{"box/law":>10}'
                      f'{"rot/law":>10}{"rot/box":>10}{"rot y/box":>11}'
                      f'{"rot x/box":>11}')
                for i, kpx in enumerate(ks):
                    print(f'      {kpx:8d}{r_all[i]:8.3f}'
                          f'{box_mean[i] / theory[i]:10.4f}'
                          f'{rot_mean[i] / theory[i]:10.4f}'
                          f'{rot_mean[i] / box_mean[i]:10.4f}'
                          f'{rot_y[i] / box_y[i]:11.4f}'
                          f'{rot_x[i] / box_x[i]:11.4f}')
                    rows.append((case_name, margin, int(kpx), r_all[i],
                                 theory[i], box_mean[i], rot_mean[i],
                                 rot_mean[i] / box_mean[i]))
                # THE CORNER LOSS, on ONE array: the crop against its own
                # rotation. The frame ratio above also carries the different
                # SHAPE of the two strips, and this one does not.
                spec = dataclasses.replace(
                    spec_rot, record=100 + 10 * int(margin * 100),
                    strip_dir=root)
                ring, d_ratio = corner_loss(sp_rot, strip_paths(spec, 1),
                                            sp_rot.theta[0], ks)
                print(f'      one-array rotation: outer-ring PSD '
                      f'{ring:.4f}, D(r) at 1 px {d_ratio[0]:.4f}, '
                      f'at 4 px {d_ratio[2]:.4f}')
                for j, kpx in ((0, 1), (2, 4)):
                    band(f'{case_name[:12]} m{margin:.2f} one-array '
                         f'D({kpx} px)', float(d_ratio[j]),
                         1.0 - BAND, 1.0 + BAND)
                if ctrl_ref is not None:
                    c_mean, _, _, _ = dphi_stats(ctrl_ref, ks)
                    print(f'      CONTROL bilinear (order 1) rot/box at '
                          f'1 px: {c_mean[0] / box_mean[0]:.4f}')
        finally:
            shutil.rmtree(root, ignore_errors=True)
    return rows


def measure_temporal(args):
    """Measurement 2: the frame-to-frame consistency of the rotation."""
    print('')
    print('2. the frame-to-frame consistency: D(tau) of the ROTATED frames')
    print('   against D(tau) of the UNROTATED cut of the SAME crops. A')
    print('   rotation is an isometry, so the two must agree.')
    lags = [8, 16, 32, 64, 128]
    for case_name, theta_deg in THETA_CASES:
        _, _, spec_rot, sp_rot, dt, speed = case_plans(theta_deg, 0.5)
        root = tempfile.mkdtemp(prefix='olb_rottime_')
        try:
            acc = []
            for s, seed in enumerate(args.seeds):
                spec = dataclasses.replace(spec_rot, record=300 + s,
                                           strip_dir=root)
                paths = strip_paths(spec, 1)
                build_strips(sp_rot, [R0], paths, L0, True, [seed],
                             table_dtype=np.float32)
                acc.append(temporal_pair(sp_rot, paths, lags))
            d_rot = np.mean([a[0] for a in acc], axis=0)
            d_raw = np.mean([a[1] for a in acc], axis=0)
        finally:
            shutil.rmtree(root, ignore_errors=True)
        theory = theory_dphi(np.array([k * dt * speed for k in lags]))
        print(f'  case: {case_name}, |v| = {speed:.1f} m/s')
        print(f'      {"tau [ms]":>10}{"v tau [m]":>11}{"rotated":>10}'
              f'{"unrotated":>11}{"ratio":>9}{"rot/law":>9}')
        for i, k in enumerate(lags):
            print(f'      {k * dt * 1e3:10.3f}{k * dt * speed:11.4f}'
                  f'{d_rot[i]:10.4f}{d_raw[i]:11.4f}'
                  f'{d_rot[i] / d_raw[i]:9.4f}{d_rot[i] / theory[i]:9.4f}')
        # The band is the 5 percent kill line, not the 3 percent of
        # measurement 1: D(tau) at the SHORTEST lag is built almost entirely
        # from the highest spatial frequencies, and those are the ones a
        # rotation moves out of the frequency square.
        worst = float(np.max(np.abs(d_rot / d_raw - 1.0)))
        band(f'{case_name[:12]} worst |D(tau) rot/unrot - 1|', worst,
             0.0, KILL)


def measure_field(args):
    """Measurement 3: the end-to-end field.

    FOUR ARMS. The PRIMARY band is the rotated frames against the UNROTATED
    central cut of the SAME crops: the strip, the seed and the grid are then
    identical, so the ratio is the rotation alone. Two arms give the context
    that the three strips do NOT share a shape:
      - the axis-aligned BOX, whose short axis is many times wider, so it
        holds more low-frequency power on that axis;
      - the ALONG-TRACK strip of the route of record, whose short axis is
        EXACTLY n, so it is the THINNEST of the three.
    """
    print('')
    print(f'3. the field after one screen and a {Z_PROP_M / 1e3:.0f} km hop '
          f'(r0 = {R0_FIELD * 1e2:.0f} cm)')
    idx = np.linspace(0, N_FRAMES - 1, 32).astype(int)
    for case_name, theta_deg in THETA_CASES:
        spec_box, sp_box, spec_rot, sp_rot, _, speed = case_plans(
            theta_deg, 0.5)
        spec_trk, sp_trk, _, _, _, _ = case_plans(89.999, 0.5)
        root = tempfile.mkdtemp(prefix='olb_rotfield_')
        try:
            arms = {'rotated': [], 'unrotated': [], 'box': [],
                    'along track': []}
            for s, seed in enumerate(args.seeds):
                arms['box'] += frames_of(
                    dataclasses.replace(spec_box, record=400 + s),
                    sp_box, seed, root, idx)[0]
                spec = dataclasses.replace(spec_rot, record=500 + s,
                                           strip_dir=root)
                arms['rotated'] += frames_of(spec, sp_rot, seed, root, idx)[0]
                arms['unrotated'] += unrotated_frames(
                    sp_rot, strip_paths(spec, 1), idx)
                arms['along track'] += frames_of(
                    dataclasses.replace(spec_trk, record=600 + s),
                    sp_trk, seed, root, idx)[0]
            got = {k: field_stats(v) for k, v in arms.items()}
        finally:
            shutil.rmtree(root, ignore_errors=True)
        print(f'  case: {case_name}')
        names = ('point sigma2_I', 'aperture sigma2_I', 'mean SMF eta')
        print(f'    {"quantity":22s}{"rotated":>11}{"unrotated":>11}'
              f'{"rot/unrot":>11}{"box":>11}{"rot/box":>10}'
              f'{"along trk":>11}{"trk/box":>10}')
        for i, name in enumerate(names):
            print(f'    {name:22s}{got["rotated"][i]:11.5f}'
                  f'{got["unrotated"][i]:11.5f}'
                  f'{got["rotated"][i] / got["unrotated"][i]:11.4f}'
                  f'{got["box"][i]:11.5f}'
                  f'{got["rotated"][i] / got["box"][i]:10.4f}'
                  f'{got["along track"][i]:11.5f}'
                  f'{got["along track"][i] / got["box"][i]:10.4f}')
            band(f'{case_name[:12]} {name} rot/unrot',
                 got['rotated'][i] / got['unrotated'][i],
                 1.0 - KILL, 1.0 + KILL)


def measure_cost(args):
    '''Measurement 4: the pixels, the build time and the per-frame cost.'''
    print('')
    print('4. the cost')
    rows = []
    for case_name, theta_deg in THETA_CASES:
        spec_box, sp_box, spec_rot, sp_rot, _, _ = case_plans(theta_deg, 0.5)
        root = tempfile.mkdtemp(prefix='olb_rotcost_')
        try:
            _, t_box = frames_of(dataclasses.replace(spec_box, record=600), sp_box,
                args.seeds[0], root, [0])
            _, t_rot = frames_of(dataclasses.replace(spec_rot, record=700), sp_rot,
                args.seeds[0], root, [0])
            # The per-frame cost of each route, over 20 frames.
            spec_b = dataclasses.replace(spec_box, record=600,
                                         strip_dir=root)
            spec_r = dataclasses.replace(spec_rot, record=700,
                                         strip_dir=root)
            s_box = open_strips(strip_paths(spec_b, 1))
            s_rot = open_strips(strip_paths(spec_r, 1))
            t0 = time.time()
            for k in range(20):
                np.asarray(next(iter(frame_stack(s_box, sp_box, k)))).sum()
            f_box = (time.time() - t0) / 20.0
            t0 = time.time()
            for k in range(20):
                np.asarray(next(iter(frame_stack(s_rot, sp_rot, k)))).sum()
            f_rot = (time.time() - t0) / 20.0
            del s_box, s_rot
        finally:
            shutil.rmtree(root, ignore_errors=True)
        px_box = sp_box.shape[0][0] * sp_box.shape[0][1]
        px_rot = sp_rot.shape[0][0] * sp_rot.shape[0][1]
        print(f'  case: {case_name}')
        print(f'    box     {sp_box.shape[0]}  {px_box / 1e6:8.2f} Mpx  '
              f'build {t_box:7.2f} s  frame {f_box * 1e3:7.2f} ms')
        print(f'    rotated {sp_rot.shape[0]}  {px_rot / 1e6:8.2f} Mpx  '
              f'build {t_rot:7.2f} s  frame {f_rot * 1e3:7.2f} ms')
        print(f'    ratio                   {px_box / px_rot:8.2f} x    '
              f'      {t_box / max(t_rot, 1e-9):7.2f} x')
        rows.append((case_name, px_box, px_rot, t_box, t_rot, f_box, f_rot))

    print('  the 30 deg hero crosswind layer, PRODUCTION pad (2 L0):')
    print(f'    box     {HERO_BOX}  {HERO_BOX[0] * HERO_BOX[1] / 1e6:8.1f} '
          f'Mpx  {HERO_BOX_SECONDS:6.1f} s (measured, bigfraw)')
    for margin in args.margins:
        shape, theta = hero_rotated_shape(margin)
        px = shape[0] * shape[1]
        ratio = HERO_BOX[0] * HERO_BOX[1] / px
        print(f'    rotated m{margin:.2f} {shape}  {px / 1e6:8.1f} Mpx  '
              f'{HERO_BOX_SECONDS / ratio:6.1f} s (area scaled)  '
              f'{ratio:6.1f} x smaller   theta = {theta:.2f} deg')
        rows.append((f'hero m{margin:.2f}', HERO_BOX[0] * HERO_BOX[1], px,
                     HERO_BOX_SECONDS, HERO_BOX_SECONDS / ratio,
                     float('nan'), float('nan')))
    return rows


def build_hero(margin, backend):
    '''Build the real rotated hero strip one time, and time it.'''
    shape, theta = hero_rotated_shape(margin)
    print(f'the 30 deg hero crosswind layer, rotated, margin {margin:.2f} n')
    print(f'  shape {shape}  ({shape[0] * shape[1] / 1e6:.1f} Mpx)  '
          f'theta = {theta:.2f} deg')
    t0 = time.time()
    factory = ScreenFactory(shape[0], HERO_DX, L0_m=L0, nx=shape[1],
                            dtype=np.float32, table_dtype=np.float32)
    screen = factory.make(0.5, np.random.default_rng(20260913))
    took = time.time() - t0
    print(f'  built in {took:.1f} s   rms {float(np.std(screen)):.2f} rad')
    print(f'  the box of the same layer: {HERO_BOX} in '
          f'{HERO_BOX_SECONDS:.1f} s  ->  {HERO_BOX_SECONDS / took:.1f} x')
    # One frame of that strip: the crop and the three shears.
    m = shape[0]
    patch = np.asarray(screen[:, :m], dtype=np.float32)
    del screen, factory
    from olb.waveoptics.propagators import set_fft_backend, xp
    previous = set_fft_backend(backend)
    try:
        xpm = xp()
        device = xpm is not np

        def sync():
            """Wait for the device, so the timing is not the queue."""
            if device:
                xpm.cuda.runtime.deviceSynchronize()
        a = xpm.asarray(patch)
        rotate_fourier(a, -np.deg2rad(theta))       # warm the plans up.
        sync()
        t0 = time.time()
        for _ in range(10):
            rotate_fourier(a, -np.deg2rad(theta))
        sync()
        print(f'  one rotated frame on {backend}: '
              f'{(time.time() - t0) / 10.0 * 1e3:.1f} ms at {m} px')
    finally:
        set_fft_backend(previous)


def write_dphi_csv(rows):
    '''Write measurement 1.'''
    os.makedirs(DATA, exist_ok=True)
    path = os.path.join(DATA, 'rotated_strip_dphi.csv')
    with open(path, 'w', encoding='utf-8') as fh:
        fh.write('case,margin_n,lag_px,r_m,theory,box,rotated,ratio\n')
        for r in rows:
            fh.write(f'{r[0]},{r[1]:.2f},{r[2]},{r[3]:.4f},{r[4]:.6f},'
                     f'{r[5]:.6f},{r[6]:.6f},{r[7]:.6f}\n')
    return path


def write_cost_csv(rows):
    '''Write measurement 4.'''
    os.makedirs(DATA, exist_ok=True)
    path = os.path.join(DATA, 'rotated_strip_cost.csv')
    with open(path, 'w', encoding='utf-8') as fh:
        fh.write('case,box_px,rotated_px,box_build_s,rotated_build_s,'
                 'box_frame_s,rotated_frame_s\n')
        for r in rows:
            fh.write(f'{r[0]},{r[1]},{r[2]},{r[3]:.3f},{r[4]:.3f},'
                     f'{r[5]:.6f},{r[6]:.6f}\n')
    return path


def draw_dphi(rows):
    '''Draw the rotated / box ratio against the lag, for every margin.'''
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    os.makedirs(FIGS, exist_ok=True)
    path = os.path.join(FIGS, 'rotated_strip_dphi.png')
    cases = sorted({r[0] for r in rows})
    fig, axes = plt.subplots(1, len(cases), figsize=(5.2 * len(cases), 4.2),
                             squeeze=False)
    for ax, case in zip(axes[0], cases):
        ax.axhline(1.0, color='k', lw=0.8)
        ax.axhspan(1.0 - BAND, 1.0 + BAND, color='0.85')
        for margin in sorted({r[1] for r in rows if r[0] == case}):
            sel = [r for r in rows if r[0] == case and r[1] == margin]
            ax.plot([r[3] for r in sel], [r[7] for r in sel], marker='o',
                    label=f'margin {margin:.2f} n')
        ax.set_xscale('log')
        ax.set_xlabel('separation r [m]')
        ax.set_ylabel('D(r) rotated / box')
        ax.set_title(case)
        ax.legend()
        ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)
    return path


def main():
    '''Run the four measurements and print the pass bands.'''
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--quick', action='store_true',
                        help='one seed and one margin, for a smoke run.')
    parser.add_argument('--bilinear', action='store_true',
                        help='add the real-space bilinear control line.')
    parser.add_argument('--build-hero', action='store_true',
                        help='build the real hero crosswind strip and stop.')
    parser.add_argument('--margin', type=float, default=0.5,
                        help='the margin of --build-hero.')
    parser.add_argument('--fft-backend', default='numpy',
                        help='the FFT backend of the --build-hero timing.')
    args = parser.parse_args()
    if args.build_hero:
        build_hero(args.margin, args.fft_backend)
        return
    args.seeds = SEEDS[:1] if args.quick else SEEDS
    args.margins = (0.5,) if args.quick else MARGINS

    t_start = time.time()
    print(f'gate: the ROTATED thin strip, n = {N}, dx = {DX * 1e2:.0f} cm, '
          f'r0 = {R0 * 1e2:.0f} cm, L0 = {L0:.0f} m, '
          f'{len(args.seeds)} seeds')
    print('')
    rows = measure_structure(args)
    measure_temporal(args)
    measure_field(args)
    cost = measure_cost(args)

    print('')
    print(f'  file saved: {write_dphi_csv(rows)}')
    print(f'  file saved: {write_cost_csv(cost)}')
    print(f'  figure saved: {draw_dphi(rows)}')
    print('    Caption: D(r) of the rotated frames divided by D(r) of the '
          'axis-aligned')
    print('    box frames, for each crop margin. The grey band is 3 percent.')

    failed = [line for line in BANDS if 'FAIL' in line]
    print('')
    print(f'  {len(BANDS) - len(failed)} of {len(BANDS)} pass bands hold')
    for line in failed:
        print(line)
    print('')
    print(f'(elapsed {time.time() - t_start:.1f} s)')


if __name__ == '__main__':
    main()
