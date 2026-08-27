'''
The Schmidt phase screens against the aotools screens of the production layer.

`olb/waveoptics/schmidt/turbulence.py` builds a Kolmogorov phase screen from
the book equations. `olb/waveoptics/turbulence/screens.py` gets the same object
from `aotools`, an outside package. The two are independent implementations of
one method, so each one validates the other. This script measures how far apart
they are, and it puts a number on the ONE defect that both of them carry.

THE STORY IS THE LOW-FREQUENCY DEFICIT. A Fourier screen holds no power below
the grid fundamental 1/(N dx), so its structure function sits BELOW the
Kolmogorov law of Ch. 9, Eq. (9.44), printed p. 160, and the gap grows with the
separation (Ch. 9, Fig. 9.3, printed p. 169). The subharmonic method of Lane,
Glindemann and Dainty (DOI 10.1088/0959-7174/2/3/003) adds three levels of low
frequency back. It lifts the curve, and it does NOT close the gap. The book
says so itself (Ch. 9, text above Sec. 9.4, printed p. 172: the match is close,
not exact).

THE THREE THINGS THIS SCRIPT DOES.

  1. ONE SCREEN, THREE WAYS. `ft_phase_screen`, `ft_sh_phase_screen` and the
     production `phase_screen`, on one r0, one N and one pitch, on a shared
     colour scale. The subharmonic pair carries a visible large-scale tilt and
     curvature that the plain Fourier screen does not.
  2. THE ENSEMBLE STRUCTURE FUNCTION. The theory curve, and the three
     generators, measured with `schmidt.fourier.structure_function` (Ch. 3,
     Eqs. (3.19) to (3.25), printed pp. 49 and 50) over a seeded ensemble.
  3. THE PER-SCREEN STRENGTH RULE, AND THE FACTOR OF FOUR. The book caps one
     screen's share of the LOG-AMPLITUDE variance sigma_chi^2 at rmax = 0.1
     (Ch. 9, Listing 9.5, lines 37 and 38, printed p. 175). The production
     planner caps one screen's share of the plane-wave RYTOV variance
     sigma_R^2 at `sigma2_r_screen_max`. Those are not the same quantity:
     sigma_R^2 = 4 sigma_chi^2, and the two constants in front are 2.25 and
     0.563. The script proves that bridge from the code, and then it puts a
     real `turbulent_grid` screen plan against BOTH caps.

THE MEASURED NUMBERS (they print below, and they are stable across seeds):
  - the schmidt subharmonic screen reaches 0.88 to 0.93 of the theory over
    r/r0 = 0.3 to 1.6;
  - the aotools screen reaches 0.91 to 0.94 there, which is 1 to 3 percent
    ABOVE the book generator. The two agree well inside the band;
  - the plain Fourier screen reaches 0.69 to 0.82 in the band, and it falls to
    0.47 at r/r0 = 8;
  - BOTH subharmonic generators fall away past the band: 0.80 and 0.86 at
    r/r0 = 8. The deficit is real, and no subharmonic level removes it.

ONE ESTIMATOR FOR ALL THREE. Every curve comes from
`schmidt.fourier.structure_function`, so the comparison carries no estimator
difference. The production self-check of
`olb/waveoptics/turbulence/screens.py` uses a direct shifted-difference
estimator instead, and it reads slightly different numbers for the same
screens. Compare the numbers of this script with each other, not with that
self-check.

This script changes NO olb module. It reads the production layer only.

Figures:
    examples/schmidt/screens_examples.png    the three screens, one colour scale
    examples/schmidt/screens_structure.png   the ensemble structure function
    examples/schmidt/screens_strength.png    the per-screen strength bars

Sources:
- Schmidt, Numerical Simulation of Optical Wave Propagation with Examples in
  MATLAB, DOI 10.1117/3.866274. Ch. 3, Eqs. (3.19) to (3.25), printed pp. 49
  and 50 (the structure-function estimator); Ch. 9, Eq. (9.44), printed p. 160
  (D(r) = 6.88 (r/r0)^(5/3)); Ch. 9, Eq. (9.52), printed p. 161 (the Kolmogorov
  phase PSD); Ch. 9, Eqs. (9.78) to (9.80), printed p. 167 (the Fourier
  screen); Ch. 9, Eq. (9.81), printed p. 169 (the subharmonics); Ch. 9,
  Eq. (9.70), printed p. 165 (the layer r0); Ch. 9, Eqs. (9.71) and (9.72),
  printed p. 165 (the composite r0); Ch. 9, Eqs. (9.73) and (9.74), printed
  p. 165 (the per-screen log-amplitude share); Ch. 9, Eq. (9.75), printed
  p. 165, and Listing 9.5, printed p. 175 (the screen-strength solve and the
  rmax cap); Ch. 9, Eq. (9.64) and the text below it, printed p. 163 (the weak
  threshold sigma_chi^2 < 0.25).
- McGlamery, J. Opt. Soc. Am. 57(3), pp. 293 to 297 (1967),
  DOI 10.1364/JOSA.57.000293. The Fourier phase screen.
- Lane, Glindemann and Dainty, Waves in Random Media 2, pp. 209 to 224 (1992),
  DOI 10.1088/0959-7174/2/3/003. The subharmonic method.
- Martin and Flatte, Appl. Opt. 27(11), pp. 2111 to 2126 (1988),
  DOI 10.1364/AO.27.002111. The per-screen strength guideline.
- Fried, DOI 10.1364/JOSA.56.001372. The Fried parameter r0.
- Andrews and Phillips, Laser Beam Propagation through Random Media, 2nd ed.,
  DOI 10.1117/3.626196, Ch. 8, Eq. (20). The plane-wave Rytov variance that the
  production planner integrates.

Run from the repo root:
    python -m examples.schmidt.screens_and_turbulence
'''

import time
import warnings

import matplotlib.pyplot as plt
import numpy as np

from olb.geometry import CircularOrbit
from olb.scenario import Channel, SpaceScenario
from olb.terminal import Terminal, Transmitter
from olb.waveoptics.schmidt.fourier import structure_function
from olb.waveoptics.schmidt.turbulence import (RMAX, composite_r0,
                                               ft_phase_screen,
                                               ft_sh_phase_screen,
                                               kolmogorov_structure_function,
                                               max_screen_strength,
                                               screen_r0 as book_screen_r0,
                                               screen_rytov_share,
                                               screen_strengths)
from olb.waveoptics.turbulence import phase_screen, turbulent_grid
from olb.waveoptics.turbulence.sampling import PRESETS, _screen_rytov

WAVELENGTH_M = 1550e-9

# ---- the screen ensemble ----
SCREEN_N = 512
SCREEN_PITCH_M = 0.01           # a 5.12 m side, that is 51 r0
SCREEN_R0_M = 0.10
N_SCREENS = 30
PUPIL_RADIUS_M = 1.2            # a guard band, so the circular estimate holds
SEED = 5000

# The separations that the table reads, in units of r0. The band that the
# production self-check of olb/waveoptics/turbulence/screens.py measures is
# r/r0 = 0.3 to 1.6, so the same band gets the stated tolerance here.
PROBES_R0 = np.array([0.3, 0.5, 0.8, 1.2, 1.6, 3.2, 8.0])
BAND_MAX_R0 = 1.6
BAND_TOL = (0.80, 1.02)

# ---- the space case that gives the real screen plan ----
SPACE_ALTITUDE_M = 600e3
SPACE_ELEVATION_DEG = 30.0
SPACE_GROUND_APERTURE_M = 0.50
SPACE_SAT_APERTURE_M = 0.30
SPACE_WAIST_M = 0.05
PRESET = "standard"

EXAMPLES_PNG = "examples/schmidt/screens_examples.png"
STRUCTURE_PNG = "examples/schmidt/screens_structure.png"
STRENGTH_PNG = "examples/schmidt/screens_strength.png"


def make_screens(index):
    '''Build one screen from each of the three generators.

    The two book generators SHARE a seeded generator per call, so the plain
    Fourier screen and the subharmonic screen of the same index carry the same
    high-frequency draw. The production screen takes its own integer seed,
    because aotools builds its own generator. The comparison is between
    ENSEMBLE means, so the draws do not have to match.

    The scales: the book generators run at L0 = infinite and l0 = 0, and the
    production generator runs at L0 = 1e6 m and l0 = 1 um. Over the sampled
    band, 1/(3^3 N dx) to 1/(2 dx), the two spectra agree to better than 1e-6,
    so both are Kolmogorov here. See Schmidt (2010), DOI 10.1117/3.866274,
    Ch. 9, Eqs. (9.51) and (9.52), printed p. 161.
    '''
    return {
        "schmidt, Fourier only": ft_phase_screen(
            SCREEN_R0_M, SCREEN_N, SCREEN_PITCH_M,
            rng=np.random.default_rng(SEED + index)),
        "schmidt, subharmonic": ft_sh_phase_screen(
            SCREEN_R0_M, SCREEN_N, SCREEN_PITCH_M,
            rng=np.random.default_rng(SEED + index)),
        "aotools, subharmonic": phase_screen(
            SCREEN_R0_M, SCREEN_N, SCREEN_PITCH_M, seed=SEED + index),
    }


def run_ensemble():
    '''Average the structure function of each generator over the ensemble.

    The estimator is `schmidt.fourier.structure_function`, Ch. 3, Eqs. (3.19)
    to (3.25), printed pp. 49 and 50. It makes the correlation CIRCULAR, so the
    window leaves a guard band at the grid edge.
    '''
    axis = (np.arange(SCREEN_N) - SCREEN_N // 2) * SCREEN_PITCH_M
    xx, yy = np.meshgrid(axis, axis)
    pupil = (np.hypot(xx, yy) <= PUPIL_RADIUS_M).astype(float)

    keys = list(make_screens(0))
    total = {k: np.zeros((SCREEN_N, SCREEN_N)) for k in keys}
    for i in range(N_SCREENS):
        for key, screen in make_screens(i).items():
            total[key] += structure_function(screen, pupil, SCREEN_PITCH_M)
    for key in keys:
        total[key] /= N_SCREENS

    # Read the x axis of the mean structure function.
    index = (SCREEN_N // 2
             + np.round(PROBES_R0 * SCREEN_R0_M / SCREEN_PITCH_M).astype(int))
    separation = (index - SCREEN_N // 2) * SCREEN_PITCH_M
    theory = kolmogorov_structure_function(separation, SCREEN_R0_M)
    rows = {k: total[k][SCREEN_N // 2, index] for k in keys}
    return total, index, separation, theory, rows


def space_scenario():
    '''Build the 600 km downlink case, so the screen plan is a real one.'''
    return SpaceScenario(
        ground=Terminal(aperture_m=SPACE_GROUND_APERTURE_M,
                        wavelength_m=WAVELENGTH_M, pointing_jitter_rad=0.0),
        space=Terminal(aperture_m=SPACE_SAT_APERTURE_M,
                       wavelength_m=WAVELENGTH_M, pointing_jitter_rad=0.0,
                       transmitter=Transmitter(waist_m=SPACE_WAIST_M)),
        direction="downlink", channel=Channel(altitude_m=SPACE_ALTITUDE_M))


def bridge_constants():
    '''Prove the factor of four between the two per-screen cap arithmetics.

    The book share of the LOG-AMPLITUDE variance is Ch. 9, Eq. (9.73), printed
    p. 165, for a plane wave:

        d(sigma_chi^2)_i = 1.33 k^(-5/6) z_to_rx^(5/6) r0_i^(-5/3)

    Put the layer r0 of Ch. 9, Eq. (9.70), printed p. 165, into it,
    r0_i^(-5/3) = 0.423 k^2 (INT Cn2 dz)_i, and the share becomes

        d(sigma_chi^2)_i = 1.33 * 0.423 k^(7/6) (INT Cn2 dz)_i z_to_rx^(5/6)
                         = 0.563 k^(7/6) (INT Cn2 dz)_i z_to_rx^(5/6)

    The production `_screen_rytov` of olb/waveoptics/turbulence/sampling.py
    writes the SAME form with the constant 2.25, from the plane-wave Rytov
    variance of Andrews and Phillips, DOI 10.1117/3.626196, Ch. 8, Eq. (20).
    And 2.25 / 0.563 = 4, which is the standard relation
    sigma_R^2 = 4 sigma_chi^2 (Ch. 9, text below Eq. (9.64), printed p. 163).

    This function measures that ratio from the two live code paths, so the
    bridge is a MEASUREMENT, not a claim.
    '''
    k = 2.0 * np.pi / WAVELENGTH_M
    cn2_int, z_to_rx = 1e-14, 5000.0
    r0_layer = float(book_screen_r0(cn2_int, WAVELENGTH_M))
    # alpha = 0, so the plane-wave weight (1 - alpha)^(5/6) is exactly 1 and
    # the distance carries the whole path factor.
    book = float(screen_rytov_share(r0_layer, 0.0, z_to_rx, WAVELENGTH_M,
                                    wave='plane'))
    prod = float(_screen_rytov(k, cn2_int, z_to_rx))
    hand_book = 1.33 * 0.423 * k ** (7.0 / 6.0) * cn2_int * z_to_rx ** (5 / 6)
    return book, prod, hand_book, prod / book


# ---------------------------------------------------------------------------
# The figures
# ---------------------------------------------------------------------------

def draw_examples(screens):
    '''Draw one screen from each generator, on one colour scale.'''
    span = max(float(np.abs(s).max()) for s in screens.values())
    fig, axes = plt.subplots(1, 3, figsize=(17.0, 5.4),
                             constrained_layout=True)
    half = SCREEN_N * SCREEN_PITCH_M / 2
    for ax, (name, screen) in zip(axes, screens.items()):
        image = ax.imshow(screen, extent=[-half, half, -half, half],
                          origin="lower", cmap="RdBu_r", vmin=-span, vmax=span)
        ax.set_title(f"{name}\npeak-to-peak "
                     f"{screen.max() - screen.min():.1f} rad", fontsize=10)
        ax.set_xlabel("x, m")
        ax.set_ylabel("y, m")
        fig.colorbar(image, ax=ax, shrink=0.82, label="phase, rad")
    fig.suptitle(f"One phase screen, three generators. r0 = "
                 f"{SCREEN_R0_M * 1e2:.0f} cm, {SCREEN_N} px at "
                 f"{SCREEN_PITCH_M * 1e2:.0f} cm, "
                 f"{WAVELENGTH_M * 1e9:.0f} nm.\n"
                 f"The two subharmonic screens carry the large-scale tilt and "
                 f"curvature that the Fourier screen is missing.", fontsize=12)
    fig.savefig(EXAMPLES_PNG, dpi=150)
    plt.close(fig)


def draw_structure(total):
    '''Draw the ensemble structure function against the Kolmogorov law.'''
    axis = (np.arange(SCREEN_N) - SCREEN_N // 2) * SCREEN_PITCH_M
    keep = (axis > 0) & (axis <= 2.0 * PUPIL_RADIUS_M)
    r_over_r0 = axis[keep] / SCREEN_R0_M
    exact = kolmogorov_structure_function(axis[keep], SCREEN_R0_M)

    colours = {"schmidt, Fourier only": "tab:green",
               "schmidt, subharmonic": "tab:blue",
               "aotools, subharmonic": "tab:red"}

    fig, (a0, a1) = plt.subplots(1, 2, figsize=(13.6, 5.6),
                                 constrained_layout=True)
    a0.loglog(r_over_r0, exact, color="black", linewidth=2.4,
              label="theory, 6.88 (r/r0)^(5/3), Ch. 9, Eq. (9.44), p. 160")
    for name, colour in colours.items():
        a0.loglog(r_over_r0, total[name][SCREEN_N // 2, keep], color=colour,
                  linewidth=1.8, label=name)
    a0.set_xlabel("separation r / r0")
    a0.set_ylabel("D(r), rad^2")
    a0.set_title(f"the ensemble structure function, {N_SCREENS} screens",
                 fontsize=10)
    a0.grid(alpha=0.3, which="both")
    a0.legend(fontsize=8, loc="upper left")

    a1.axhline(1.0, color="black", linewidth=2.4, label="theory")
    for name, colour in colours.items():
        a1.semilogx(r_over_r0,
                    total[name][SCREEN_N // 2, keep] / exact, color=colour,
                    linewidth=1.8, label=name)
    a1.axvspan(PROBES_R0[0], BAND_MAX_R0, color="0.85", zorder=0,
               label=f"the stated band, {PROBES_R0[0]} to {BAND_MAX_R0} r0")
    a1.set_ylim(0.0, 1.15)
    a1.set_xlabel("separation r / r0")
    a1.set_ylabel("measured D(r) / theory")
    a1.set_title("the same, as a ratio: the low-frequency deficit",
                 fontsize=10)
    a1.grid(alpha=0.3)
    a1.legend(fontsize=8, loc="lower left")

    fig.suptitle("The low-frequency deficit of a Fourier phase screen. The "
                 "subharmonics lift the curve; they do not close the gap.",
                 fontsize=12)
    fig.savefig(STRUCTURE_PNG, dpi=150)
    plt.close(fig)


def draw_strength(plan, preset, sigma2_chi, book_r0):
    '''Draw the per-screen shares against the two caps, and the two r0 sets.'''
    index = np.arange(plan.z_m.size)
    fig, (a0, a1) = plt.subplots(1, 2, figsize=(14.4, 5.6),
                                 constrained_layout=True)

    a0.bar(index, sigma2_chi, color="tab:blue",
           label="the share of each screen, sigma_chi^2")
    a0.axhline(RMAX, color="tab:red", linewidth=2.0, linestyle="--",
               label=f"book cap, rmax = {RMAX} (Listing 9.5, p. 175)")
    a0.axhline(preset.sigma2_r_screen_max / 4.0, color="tab:green",
               linewidth=2.0, linestyle=":",
               label=f"production cap, {preset.sigma2_r_screen_max} of "
                     f"sigma_R^2 = {preset.sigma2_r_screen_max / 4:.3f} of "
                     f"sigma_chi^2")
    a0.set_yscale("log")
    a0.set_xlabel("screen index, from the top of the slab")
    a0.set_ylabel("log-amplitude share, sigma_chi^2")
    a0.set_title(f"the per-screen strength of a real plan\n"
                 f"{plan.z_m.size} screens, total sigma_chi^2 = "
                 f"{sigma2_chi.sum():.4f}", fontsize=10)
    a0.grid(alpha=0.3, axis="y")
    a0.legend(fontsize=8, loc="lower right")

    width = 0.4
    a1.bar(index - width / 2, plan.r0_m * 1e2, width, color="tab:blue",
           label="turbulent_grid, from the Cn2 profile")
    a1.bar(index + width / 2, book_r0 * 1e2, width, color="tab:orange",
           label="screen_strengths, the Ch. 9, Eq. (9.75) solve")
    a1.set_yscale("log")
    a1.set_xlabel("screen index, from the top of the slab")
    a1.set_ylabel("screen r0, cm")
    a1.set_title("the two ways to give a screen its strength\n"
                 "a Cn2 integral against a two-equation least-squares solve",
                 fontsize=10)
    a1.grid(alpha=0.3, axis="y")
    a1.legend(fontsize=8)

    fig.suptitle(f"The per-screen strength rule. Space downlink slab at "
                 f"{SPACE_ELEVATION_DEG:.0f} deg, {PRESET} preset.\n"
                 f"The two caps differ by exactly 4, because "
                 f"sigma_R^2 = 4 sigma_chi^2.", fontsize=12)
    fig.savefig(STRENGTH_PNG, dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------

def main():
    t_start = time.time()

    print("=" * 78)
    print("The Schmidt phase screens against the aotools screens of olb")
    print("=" * 78)
    print(f"  wavelength              {WAVELENGTH_M * 1e9:11.1f} nm")
    print(f"  screen grid             {SCREEN_N:11d} px at "
          f"{SCREEN_PITCH_M * 1e2:.1f} cm, side "
          f"{SCREEN_N * SCREEN_PITCH_M:.2f} m")
    print(f"  r0                      {SCREEN_R0_M * 1e2:11.1f} cm "
          f"({SCREEN_R0_M / SCREEN_PITCH_M:.0f} px, the side is "
          f"{SCREEN_N * SCREEN_PITCH_M / SCREEN_R0_M:.0f} r0)")
    print(f"  ensemble                {N_SCREENS:11d} screens per generator")
    print(f"  structure-function pupil{PUPIL_RADIUS_M:11.2f} m radius "
          f"(a guard band for the circular estimate)")

    t0 = time.time()
    total, index, separation, theory, rows = run_ensemble()
    t_ensemble = time.time() - t0

    print("")
    print("1. ONE SCREEN, THREE WAYS")
    example = make_screens(0)
    for name, screen in example.items():
        print(f"  {name:<24}rms {screen.std():7.3f} rad   "
              f"peak-to-peak {screen.max() - screen.min():7.3f} rad")
    print("  The subharmonic screens carry more rms, and the excess is all at")
    print("  low frequency: it is the tilt and the curvature that the plain")
    print("  Fourier screen cannot hold (Ch. 9, text below Listing 9.2,")
    print("  printed p. 167).")

    print("")
    print("2. THE ENSEMBLE STRUCTURE FUNCTION, against Ch. 9, Eq. (9.44), "
          "p. 160")
    print("")
    keys = list(rows)
    head = (f"  {'r/r0':>6}{'D theory':>11}"
            + "".join(f"{k.split(', ')[1][:9]:>11}{'ratio':>8}" for k in keys))
    print(f"  {'':>6}{'':>11}"
          + "".join(f"{k.split(',')[0][:9]:>19}" for k in keys))
    print(head)
    for j, probe in enumerate(PROBES_R0):
        line = f"  {probe:>6.1f}{theory[j]:>11.2f}"
        for key in keys:
            line += f"{rows[key][j]:>11.2f}{rows[key][j] / theory[j]:>8.3f}"
        print(line)

    band = PROBES_R0 <= BAND_MAX_R0
    ratio_sh = rows["schmidt, subharmonic"] / theory
    ratio_ao = rows["aotools, subharmonic"] / theory
    ratio_ft = rows["schmidt, Fourier only"] / theory
    print("")
    print(f"  the stated band is r/r0 = {PROBES_R0[0]} to {BAND_MAX_R0}, "
          f"tolerance {BAND_TOL[0]} to {BAND_TOL[1]}")
    print(f"  schmidt subharmonic in the band  "
          f"{ratio_sh[band].min():.3f} to {ratio_sh[band].max():.3f}")
    print(f"  aotools subharmonic in the band  "
          f"{ratio_ao[band].min():.3f} to {ratio_ao[band].max():.3f}")
    gap = 100.0 * ((ratio_ao / ratio_sh)[band] - 1.0)
    print(f"  aotools against schmidt          "
          f"{gap.min():+.1f} to {gap.max():+.1f} percent")
    print(f"  Fourier only, everywhere         "
          f"{ratio_ft.min():.3f} to {ratio_ft.max():.3f}")
    print(f"  (the ensemble took {t_ensemble:.1f} s)")

    # The asserts. They repeat the tolerance philosophy of the two self-checks
    # that this script cross-checks: assert a BAND, not a point.
    assert np.all(ratio_sh[band] > BAND_TOL[0]), ratio_sh[band]
    assert np.all(ratio_sh[band] < BAND_TOL[1]), ratio_sh[band]
    assert np.all(ratio_ao[band] > BAND_TOL[0]), ratio_ao[band]
    assert np.all(ratio_ao[band] < BAND_TOL[1]), ratio_ao[band]
    # The plain Fourier screen is LOW everywhere, and the deficit grows with
    # the separation (Ch. 9, Fig. 9.3, printed p. 169).
    assert np.all(ratio_ft < ratio_sh), (ratio_ft, ratio_sh)
    assert ratio_ft[-1] < ratio_ft[0], ratio_ft
    # Neither subharmonic screen closes the gap at a large separation.
    assert ratio_sh[-1] < 0.90, ratio_sh[-1]
    assert ratio_ao[-1] < 0.90, ratio_ao[-1]
    print("  The two subharmonic generators agree inside the band, and BOTH")
    print("  fall away past it. That is the book's own statement (Ch. 9, text")
    print("  above Sec. 9.4, printed p. 172): the match is close, not exact.")

    # -----------------------------------------------------------------
    print("")
    print("3. THE PER-SCREEN STRENGTH RULE, AND THE FACTOR OF FOUR")
    book, prod, hand_book, ratio = bridge_constants()
    print("")
    print("  One layer, Cn2 dz = 1e-14 m^(1/3), 5 km from the receiver, "
          "plane wave:")
    print(f"    schmidt.screen_rytov_share (sigma_chi^2)    {book:.6e}")
    print(f"    1.33 * 0.423 k^(7/6) Cn2dz z^(5/6), by hand {hand_book:.6e}")
    print(f"    production _screen_rytov   (sigma_R^2)      {prod:.6e}")
    print(f"    ratio, production over book                 {ratio:.6f}")
    assert abs(book / hand_book - 1.0) < 1e-9, (book, hand_book)
    # The two printed constants are ROUNDED: the book prints 1.33 and 0.423,
    # and Andrews and Phillips print 2.25. So the ratio is 4 to three decimal
    # places, not to machine precision.
    assert abs(ratio / 4.0 - 1.0) < 1e-3, ratio
    print("  The two constants are 2.25 and 1.33 * 0.423 = 0.5626, and their")
    print(f"  ratio is {ratio:.5f}. That is 4 to the rounding of the two "
          f"PRINTED")
    print("  constants. The two layers cap DIFFERENT quantities:")
    print("  the book caps the log-amplitude variance sigma_chi^2 of Ch. 9,")
    print("  Eq. (9.64), printed p. 163, and the production planner caps the")
    print("  plane-wave Rytov variance sigma_R^2 = 4 sigma_chi^2.")

    space = space_scenario()
    orbit = CircularOrbit(altitude_m=SPACE_ALTITUDE_M,
                          elevation_deg=[SPACE_ELEVATION_DEG])
    with warnings.catch_warnings(record=True):
        warnings.simplefilter("always")
        grid, plan, report = turbulent_grid(space, orbit, preset=PRESET)
    preset = PRESETS[PRESET]
    sigma2_chi = plan.sigma2_r / 4.0
    alpha = plan.z_m / plan.z_total_m

    print("")
    print(f"  A REAL PLAN. The downlink slab at {SPACE_ELEVATION_DEG:.0f} deg, "
          f"{PRESET} preset:")
    print(f"    slab length           {plan.z_total_m * 1e-3:11.3f} km")
    print(f"    screens               {plan.z_m.size:11d}")
    print(f"    r0 of the whole path  {plan.r0_total_m * 1e2:11.3f} cm")
    print(f"    sigma_R^2 total       {plan.sigma2_r.sum():11.4f}")
    print(f"    sigma_chi^2 total     {sigma2_chi.sum():11.4f}   "
          f"(weak below 0.25, Ch. 9, p. 163)")
    print(f"    strongest screen      {sigma2_chi.max():11.4f} of "
          f"sigma_chi^2")
    print(f"    book cap rmax         {RMAX:11.4f}")
    print(f"    production cap        "
          f"{preset.sigma2_r_screen_max / 4.0:11.4f}   "
          f"({preset.sigma2_r_screen_max} of sigma_R^2)")
    print(f"    the production cap is "
          f"{RMAX / (preset.sigma2_r_screen_max / 4.0):.1f} times STRICTER "
          f"than the book cap")
    assert sigma2_chi.max() <= RMAX, sigma2_chi.max()
    assert sigma2_chi.sum() < 0.25, sigma2_chi.sum()

    # The composite r0 of the two layers must be the same number, because the
    # two add the screens the same way: Ch. 9, Eq. (9.71), printed p. 165.
    r0_book = composite_r0(plan.r0_m, wave='plane')
    print(f"    composite r0, book Eq. (9.71) {r0_book * 1e2:9.4f} cm")
    print(f"    composite r0, turbulent_grid  {plan.r0_total_m * 1e2:9.4f} cm")
    print(f"    ratio                         "
          f"{r0_book / plan.r0_total_m:9.6f}")
    assert abs(r0_book / plan.r0_total_m - 1.0) < 1e-9, r0_book

    # The book's own screen-strength solve, on the SAME screen positions.
    solved = screen_strengths(alpha, plan.r0_total_m, sigma2_chi.sum(),
                              plan.z_total_m, WAVELENGTH_M)
    cap_x2 = max_screen_strength(alpha, plan.z_total_m, WAVELENGTH_M)
    print("")
    print("  THE OTHER ROUTE TO A SCREEN STRENGTH. Ch. 9, Eq. (9.75), printed")
    print("  p. 165, does not read a Cn2 profile at all. It solves a")
    print("  two-equation bounded least-squares problem for r0_i, so that the")
    print("  layering reproduces the path r0 and the path sigma_chi^2. The")
    print("  system is UNDERDETERMINED, so the answer is not the profile.")
    print("")
    print(f"    {'i':>3}{'alpha':>8}{'r0 plan, cm':>14}{'r0 solved, cm':>16}"
          f"{'sigma_chi^2':>14}{'cap on r0, cm':>15}")
    # An infinite bound on r0_i^(-5/3) is NO constraint, so the smallest
    # allowed r0 is zero there. `max_screen_strength` returns that infinity for
    # a screen whose spherical-wave path weight is zero, that is at alpha = 0
    # (Schmidt (2010), DOI 10.1117/3.866274, Ch. 9, Eq. (9.74), printed
    # p. 165).
    cap_r0 = np.where(np.isfinite(cap_x2), cap_x2 ** (-3.0 / 5.0), 0.0)
    for i in range(plan.z_m.size):
        cap_text = ("none" if cap_r0[i] == 0.0
                    else f"{cap_r0[i] * 1e2:.3f}")
        print(f"    {i:>3}{alpha[i]:>8.4f}{plan.r0_m[i] * 1e2:>14.3f}"
              f"{solved[i] * 1e2:>16.3f}{sigma2_chi[i]:>14.5f}"
              f"{cap_text:>15}")
    print("    The `cap on r0` column is `max_screen_strength` inverted: it is")
    print("    the SMALLEST r0 that a screen at that alpha may have and still")
    print("    obey rmax. A screen weaker than the cap has a larger r0. The")
    print("    first screen sits at alpha = 0, where the spherical-wave path")
    print("    weight of Ch. 9, Eq. (9.74), printed p. 165, is zero, so it")
    print("    carries no cap at all.")
    print("    THE SOLVED COLUMN IS A SPHERICAL-WAVE ANSWER (Ch. 9, Eq. (9.74),")
    print("    printed p. 165) on a PLANE-WAVE plan, so read the METHOD, not")
    print("    the agreement: the solve puts nearly all of the strength into a")
    print("    few screens, and the profile spreads it over the whole slab.")
    # Every screen of the plan must obey the book cap on r0_i.
    assert np.all(plan.r0_m >= cap_r0 - 1e-12), (plan.r0_m, cap_r0)

    draw_examples(example)
    draw_structure(total)
    draw_strength(plan, preset, sigma2_chi, solved)

    print("")
    print(f"figure saved: {EXAMPLES_PNG}")
    print("  Caption: one screen from each generator on one colour scale. The "
          "subharmonic")
    print("  pair carries the large-scale tilt that the plain Fourier screen "
          "cannot hold.")
    print(f"figure saved: {STRUCTURE_PNG}")
    print("  Caption: the ensemble structure function against the Kolmogorov "
          "law, and the")
    print("  same as a ratio. The grey band is the stated tolerance band; both "
          "subharmonic")
    print("  generators sit inside it and both fall away past it.")
    print(f"figure saved: {STRENGTH_PNG}")
    print("  Caption: the per-screen log-amplitude share of a real screen plan "
          "against the")
    print("  two caps, and the two ways to give a screen its strength.")
    print("")
    print(f"(elapsed {time.time() - t_start:.1f} s)")


if __name__ == "__main__":
    main()
