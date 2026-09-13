'''
Gate (a) of the temporal strip screens: the RECTANGULAR ScreenFactory.

A frozen-flow record crops its frames out of ONE long STRIP screen (see
`olb/waveoptics/turbulence/temporal.py`). So a crop of a strip must hold the
same turbulence as a square screen of the same side. This script measures that.

THREE MEASUREMENTS.

  1. THE BIT-IDENTITY. `ScreenFactory(n, dx, nx=n)` must give the screen of
     `ScreenFactory(n, dx)`, bit for bit, in both precisions and at both outer
     scales. Every stored campaign depends on it.
  2. THE STRUCTURE FUNCTION. A (512, 4096) strip against a 512 square ensemble,
     on the y axis AND on the x axis, and against the analytic von Karman law.
     The band is 2 standard errors of the difference.
  3. THE TILT. The Z-tilt of a 1 m pupil, on 32 crops of the strips against 32
     square screens. The tilt is the mode that a band-limited screen loses
     first, so it is the hardest of the three.

WHY THE STRIP NEEDS ITS OWN LOW-FREQUENCY RULE. The subharmonics of the factory
keep the spacing of the SHORT axis, and the factory cuts the matching
low-frequency box out of the main grid. Without that cut a strip reads D(r)
along y 8 to 23 percent LOW and D(r) along x 1 to 12 percent HIGH: the near-DC
columns of the long axis hold a large power, and the main grid gives all of it
to the frequency fy = 0, which carries no y structure at all.

Sources:
- Schmidt, Numerical Simulation of Optical Wave Propagation with Examples in
  MATLAB, DOI 10.1117/3.866274. Ch. 9, Eq. (9.44), printed p. 160 (the
  Kolmogorov structure function); Ch. 9, Eqs. (9.78) to (9.81), printed
  pp. 166 to 169 (the Fourier screen and its subharmonics).
- Assemat and Wilson, Opt. Express 14(3), pp. 988 to 999 (2006),
  DOI 10.1364/OE.14.000988, Eq. (5). The closed-form von Karman phase
  covariance, which gives the analytic D(r) = 2 (B(0) - B(r)).
- Noll, J. Opt. Soc. Am. 66(3), pp. 207 to 211 (1976),
  DOI 10.1364/JOSA.66.000207. The Zernike tilt.
- Lane, Glindemann and Dainty, DOI 10.1088/0959-7174/2/3/003. The subharmonic
  method.

Outputs, next to this script:
    data/rect_dphi.csv      the structure function of every arm
    data/rect_tilt.csv      the Z-tilt variance of every arm
    figures/rect_dphi.png   the structure function against the analytic law

Run from the repo root:
    python -m validation.temporal_screens.rect_factory
'''

import os
import time

import numpy as np

from olb.waveoptics.turbulence.screens import ScreenFactory
from validation.screens.helpers import pupil_mask, vk_covariance_closed, \
    zernike_tilt

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, 'data')
FIGS = os.path.join(HERE, 'figures')

# ---- the case ----
N = 512                      # the short side, in pixels.
NX = 4096                    # the long side, in pixels.
DX = 0.01                    # the pitch, in m.
R0 = 0.10                    # the Fried parameter, in m.
L0 = 25.0                    # the operating outer scale, in m.
LAM = 1550e-9
M_STRIP = 16                 # strips in the ensemble.
M_SQUARE = 32                # square screens in the ensemble.
SEED_STRIP = 20260913
SEED_SQUARE = 20260914
R_LIST = np.array([0.03, 0.08, 0.16, 0.32, 0.64, 1.00])   # up to 1 m.
PUPIL_M = 1.0
N_CROPS = 32

BANDS = []


def band(name, value, low, high):
    '''Record one PASS or FAIL line, and print it.'''
    ok = low <= value <= high
    line = (f"  {'PASS' if ok else 'FAIL'}  {name:52s} {value:9.4f}  "
            f"band [{low:.4f}, {high:.4f}]")
    BANDS.append(line)
    print(line)
    return ok


def theory_dphi(r_m):
    '''Give the analytic von Karman D(r) = 2 (B(0) - B(r)).'''
    b = vk_covariance_closed(np.concatenate(([0.0], np.asarray(r_m))), R0, L0)
    return 2.0 * (b[0] - b[1:])


def dphi_axes(screen, ks):
    '''Give D(r) along y and along x of one screen, at the pixel lags ks.'''
    dy = [float(np.mean((screen[k:, :] - screen[:-k, :]) ** 2)) for k in ks]
    dx = [float(np.mean((screen[:, k:] - screen[:, :-k]) ** 2)) for k in ks]
    return np.array(dy), np.array(dx)


def mean_se(samples):
    '''Give the mean and the standard error of a (M, K) sample array.'''
    a = np.asarray(samples, dtype=float)
    return a.mean(axis=0), a.std(axis=0, ddof=1) / np.sqrt(a.shape[0])


def check_bit_identity():
    '''Measurement 1: nx = n must not move one bit.'''
    print('1. the bit-identity of nx = n')
    ok = True
    for dtype in (np.float64, np.float32):
        for l0 in (L0, np.inf):
            square = ScreenFactory(128, DX, L0_m=l0, dtype=dtype).make(
                R0, np.random.default_rng(5))
            rect = ScreenFactory(128, DX, L0_m=l0, dtype=dtype, nx=128).make(
                R0, np.random.default_rng(5))
            same = bool(np.array_equal(square, rect))
            ok = ok and same
            print(f"  {'PASS' if same else 'FAIL'}  {dtype.__name__:8s} "
                  f"L0 = {l0!s:5s}  bit-identical {same}")
    BANDS.append(f"  {'PASS' if ok else 'FAIL'}  nx = n is bit-identical")
    return ok


def measure_dphi():
    '''Measurement 2: the structure function of the strip and of the square.'''
    ks = np.rint(R_LIST / DX).astype(int)
    strip_y, strip_x, square = [], [], []
    fac_strip = ScreenFactory(N, DX, L0_m=L0, nx=NX)
    fac_square = ScreenFactory(N, DX, L0_m=L0)
    crops = []
    for i in range(M_STRIP):
        s = fac_strip.make(R0, np.random.default_rng(SEED_STRIP + i))
        dy, dx = dphi_axes(s, ks)
        strip_y.append(dy)
        strip_x.append(dx)
        if len(crops) < N_CROPS:
            # Disjoint square crops, so the tilt samples are independent.
            step = NX // (N_CROPS // M_STRIP + 1)
            for j in range(N_CROPS // M_STRIP):
                crops.append(np.array(s[:, j * step:j * step + N]))
        del s
    for i in range(M_SQUARE):
        s = fac_square.make(R0, np.random.default_rng(SEED_SQUARE + i))
        dy, dx = dphi_axes(s, ks)
        square.append(0.5 * (dy + dx))
        del s
    return (mean_se(strip_y), mean_se(strip_x), mean_se(square), crops)


def measure_tilt(crops):
    '''Measurement 3: the Z-tilt variance of a 1 m pupil.'''
    mask = pupil_mask(N, DX, PUPIL_M)
    fac_square = ScreenFactory(N, DX, L0_m=L0)
    a_strip, a_square = [], []
    for crop in crops:
        _, _, a2, a3 = zernike_tilt(crop, mask, DX, LAM)
        a_strip.append([a2 ** 2, a3 ** 2])
    for i in range(len(crops)):
        s = fac_square.make(R0, np.random.default_rng(SEED_SQUARE + 500 + i))
        _, _, a2, a3 = zernike_tilt(s, mask, DX, LAM)
        a_square.append([a2 ** 2, a3 ** 2])
    # The two axes are two samples of one variance, so flatten them.
    return (mean_se(np.array(a_strip).reshape(-1, 1)),
            mean_se(np.array(a_square).reshape(-1, 1)))


def write_dphi_csv(theory, sy, sx, sq):
    '''Write the structure function of every arm.'''
    os.makedirs(DATA, exist_ok=True)
    path = os.path.join(DATA, 'rect_dphi.csv')
    with open(path, 'w', encoding='utf-8') as fh:
        fh.write('r_m,theory,strip_y,strip_y_se,strip_x,strip_x_se,'
                 'square,square_se\n')
        for i, r in enumerate(R_LIST):
            fh.write(f'{r:.4f},{theory[i]:.6f},{sy[0][i]:.6f},{sy[1][i]:.6f},'
                     f'{sx[0][i]:.6f},{sx[1][i]:.6f},{sq[0][i]:.6f},'
                     f'{sq[1][i]:.6f}\n')
    return path


def write_tilt_csv(strip, square):
    '''Write the Z-tilt variance of the two arms.'''
    os.makedirs(DATA, exist_ok=True)
    path = os.path.join(DATA, 'rect_tilt.csv')
    with open(path, 'w', encoding='utf-8') as fh:
        fh.write('arm,a_squared_rad2,standard_error\n')
        fh.write(f'strip_crop,{strip[0][0]:.6f},{strip[1][0]:.6f}\n')
        fh.write(f'square,{square[0][0]:.6f},{square[1][0]:.6f}\n')
    return path


def draw_dphi(theory, sy, sx, sq):
    '''Draw the structure function of every arm, divided by the theory.'''
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    os.makedirs(FIGS, exist_ok=True)
    path = os.path.join(FIGS, 'rect_dphi.png')
    fig, ax = plt.subplots(figsize=(7.0, 4.2))
    ax.axhline(1.0, color='k', lw=0.8)
    for (mean, se), label, colour in ((sy, 'strip, along y', 'tab:blue'),
                                      (sx, 'strip, along x', 'tab:orange'),
                                      (sq, 'square 512', 'tab:green')):
        ax.errorbar(R_LIST, mean / theory, yerr=2.0 * se / theory,
                    marker='o', capsize=3, label=label, color=colour)
    ax.set_xscale('log')
    ax.set_xlabel('separation r [m]')
    ax.set_ylabel('D(r) / analytic von Karman')
    ax.set_title(f'({N}, {NX}) strip against a {N} square, '
                 f'r0 = {R0 * 1e2:.0f} cm, L0 = {L0:.0f} m')
    ax.legend()
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)
    return path


def main():
    '''Run the three measurements and print the pass bands.'''
    t_start = time.time()
    print(f'gate (a): the rectangular ScreenFactory, ({N}, {NX}) at '
          f'dx = {DX * 1e2:.0f} cm, r0 = {R0 * 1e2:.0f} cm, L0 = {L0:.0f} m')
    print('')
    check_bit_identity()

    print('')
    print('2. the structure function (2 SE of the difference)')
    theory = theory_dphi(R_LIST)
    sy, sx, sq, crops = measure_dphi()
    print(f'{"r [m]":>8}{"theory":>10}{"strip y":>10}{"strip x":>10}'
          f'{"square":>10}{"y/sq":>8}{"x/sq":>8}{"y/th":>8}{"x/th":>8}')
    for i, r in enumerate(R_LIST):
        print(f'{r:8.3f}{theory[i]:10.3f}{sy[0][i]:10.3f}{sx[0][i]:10.3f}'
              f'{sq[0][i]:10.3f}{sy[0][i] / sq[0][i]:8.3f}'
              f'{sx[0][i] / sq[0][i]:8.3f}{sy[0][i] / theory[i]:8.3f}'
              f'{sx[0][i] / theory[i]:8.3f}')
    print('')
    for i, r in enumerate(R_LIST):
        for (mean, se), axis in ((sy, 'y'), (sx, 'x')):
            tol = 2.0 * np.hypot(se[i], sq[1][i]) / sq[0][i]
            band(f'strip {axis} / square at r = {r:.2f} m',
                 mean[i] / sq[0][i], 1.0 - tol, 1.0 + tol)

    print('')
    print('3. the Z-tilt of a 1 m pupil (2 SE of the difference)')
    t_strip, t_square = measure_tilt(crops)
    tol = 2.0 * np.hypot(t_strip[1][0], t_square[1][0]) / t_square[0][0]
    print(f'  strip crops  <a^2> = {t_strip[0][0]:.4f} +/- '
          f'{t_strip[1][0]:.4f} rad^2 ({len(crops)} crops, two axes)')
    print(f'  square       <a^2> = {t_square[0][0]:.4f} +/- '
          f'{t_square[1][0]:.4f} rad^2')
    band('strip Z-tilt / square Z-tilt', t_strip[0][0] / t_square[0][0],
         1.0 - tol, 1.0 + tol)

    print('')
    print(f'  file saved: {write_dphi_csv(theory, sy, sx, sq)}')
    print(f'  file saved: {write_tilt_csv(t_strip, t_square)}')
    print(f'  figure saved: {draw_dphi(theory, sy, sx, sq)}')
    print('    Caption: the ensemble structure function of the strip on each '
          'axis and of')
    print('    the square screen, divided by the analytic von Karman law. The '
          'bars are')
    print('    two standard errors.')

    failed = [line for line in BANDS if 'FAIL' in line]
    print('')
    print(f'  {len(BANDS) - len(failed)} of {len(BANDS)} pass bands hold')
    for line in failed:
        print(line)
    print('')
    print(f'(elapsed {time.time() - t_start:.1f} s)')


if __name__ == '__main__':
    main()
