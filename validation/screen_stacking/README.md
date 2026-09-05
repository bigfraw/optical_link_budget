# The phase-screen stacking test

Does a STACK of phase screens hold the statistics of ONE screen of the same
composite `r0`? Phase only, no propagation. It is the generator half of the
tail-convergence study (`validation/tail_convergence/`, backlog 2-I2T), and it
adds a measurement to backlog 2-N2.

## The question

The FFT-plus-subharmonic screen generator is close to Kolmogorov, but not
faithful: it misses the low-frequency (tip and tilt) band below the grid
fundamental (`validation/screens/FINDINGS.md`, Q1: the three-level subharmonic
screen holds 0.75 to 0.78 of the Noll Z-tilt variance). The owner's hypothesis
(2026-09-04): when a split-step plan divides the SAME turbulence among MORE
screens, the stack loses MORE of that band than one screen would, so a
many-screen plan reads a little LESS effective turbulence and a little LESS
fibre-coupling loss. If that holds, a screen-count sweep carries a generator
bias.

The tail-convergence study measures the PRODUCT of the screens and the Fresnel
propagation between them. This test isolates the SCREENS.

## The method

Each configuration is a list of per-screen `r0` values taken from the plans of
the tail-convergence study at 30 deg. The script draws 100 stacks with the
production `ScreenFactory` on the pinned 1024 px grid (3.43 mm pixel), SUMS the
screens of a stack into one phase map, and measures three things against the
Kolmogorov theory of the COMPOSITE `r0` (`r0c = (SUM r0_i^(-5/3))^(-3/5)`,
Andrews and Phillips, DOI 10.1117/3.626196, Ch. 12, Eq. (23)):

1. The phase structure function `D(r)` at six separations, averaged over the
   whole grid. Theory: `D(r) = 6.88 (r/r0)^(5/3)`, Fried,
   DOI 10.1364/JOSA.56.001372.
2. `Delta1`, the piston-removed phase variance over a 0.7 m aperture. Theory:
   `1.0299 (D/r0)^(5/3)`, Noll, DOI 10.1364/JOSA.66.000207, Table IV.
3. `Delta3`, the tip-tilt-removed phase variance over the same aperture.
   Theory: `0.134 (D/r0)^(5/3)`, the same table.

Every number is a RATIO to its theory; 1.0 is a faithful screen. The 25
disjoint 0.7 m apertures that tile the grid interior give the aperture
statistics of one draw, and the standard error comes from the spread of the
per-draw means over the 100 draws, so the correlation between the apertures of
one draw does not fake a small bar.

Run it from the repository root:

    python -m validation.screen_stacking.screen_stacking [--draws 100]

It imports the case builder of `validation.tail_convergence.tail_convergence`
for the plans, so the two studies read the same screens.

## The results (2026-09-04, 100 draws for each configuration)

| configuration | screens | r0c [cm] | D(0.05) | D(0.1) | D(0.2) | D(0.35) | D(0.7) | D(1.4) | Delta1 | Delta3 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| ground layer, 1 screen (the `pin09` bottom screen) | 1 | 15.38 | 0.954 | 0.880 | 0.851 | 0.829 | 0.783 | 0.721 | 0.823 +-0.033 | 1.005 +-0.010 |
| ground layer, 4 sub-screens (`gnd09x4`) | 4 | 15.38 | 0.962 | 0.888 | 0.861 | 0.841 | 0.793 | 0.716 | 0.839 +-0.036 | 0.996 +-0.011 |
| whole plan, 5 screens (`pin05`) | 5 | 12.73 | 0.935 | 0.857 | 0.821 | 0.790 | 0.729 | 0.638 | 0.772 +-0.024 | 1.022 +-0.010 |
| whole plan, 9 screens (`pin09`) | 9 | 12.73 | 0.931 | 0.852 | 0.816 | 0.788 | 0.739 | 0.680 | 0.765 +-0.027 | 1.014 +-0.010 |
| whole plan, 25 screens (`pin25`) | 25 | 12.73 | 0.920 | 0.838 | 0.797 | 0.761 | 0.697 | 0.592 | 0.748 +-0.025 | 1.005 +-0.010 |

## The results at a FINITE outer scale (2026-09-04, `--L0 25`, 100 draws)

The screens are drawn with a von Karman outer scale of 25 m AND the reference
is the von Karman theory at 25 m (the Noll piston and tilt filters integrated
over the von Karman PSD, and the closed-form covariance of Assemat and Wilson).
The reference route reproduces Fried and Noll at L0 = 1e7 m inside 1 percent
(the script asserts it).

| configuration | screens | D(0.05) | D(0.1) | D(0.2) | D(0.35) | D(0.7) | D(1.4) | Delta1 | Delta3 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| ground layer, 1 screen | 1 | 1.048 | 0.993 | 0.993 | 1.002 | 0.994 | 0.960 | 1.000 +-0.026 | 0.996 +-0.010 |
| ground layer, 4 sub-screens | 4 | 1.049 | 0.994 | 0.993 | 1.002 | 0.990 | 0.932 | 1.004 +-0.027 | 0.987 +-0.011 |
| whole plan, 5 screens | 5 | 1.046 | 0.991 | 0.989 | 0.993 | 0.979 | 0.927 | 0.974 +-0.022 | 1.013 +-0.010 |
| whole plan, 9 screens | 9 | 1.033 | 0.975 | 0.969 | 0.969 | 0.961 | 0.946 | 0.946 +-0.026 | 1.005 +-0.010 |
| whole plan, 25 screens | 25 | 1.029 | 0.969 | 0.960 | 0.956 | 0.939 | 0.862 | 0.942 +-0.026 | 0.996 +-0.010 |

## The verdict

0. **The L0 = inf "deficit" is the outer scale the grid cannot hold, not a
   broken generator.** Three subharmonic levels reach 27 x the grid side, so
   the 3.51 m grid holds scales up to about 95 m and nothing beyond, while the
   Kolmogorov reference assumes an INFINITE outer scale. The measured L0 = inf
   screens match the von Karman theory at L0 = 95 m almost exactly (9-screen
   plan: `Delta1` 0.765 measured against 0.765 theory; D(0.7 m) 0.739 against
   0.712; D(1.4 m) 0.680 against 0.637). Asked for L0 = 25 m and judged
   against it, the screens read `Delta1` 1.000 +-0.026 for one screen and D(r)
   inside 1 to 5 percent of the theory up to r = 0.7 m. So the production
   default `L0_m = inf` CLAIMS an outer scale it does not deliver. The
   physics behind the choice is large: the von Karman theory (the Noll piston
   filter over the von Karman PSD) gives a piston-removed aperture variance of
   0.630 of the Kolmogorov value at L0 = 25 m and 0.765 at L0 = 95 m, for
   D = 0.7 m. The fibre pays that tilt, so the fibre fade tail moves with the
   choice by a MEASURED 2.5 dB (30 deg) to 2.8 dB (20 deg) at p5
   (2026-09-05, `validation/outer_scale_tail/`, a matched-seed `L0 = inf`
   against `L0 = 25 m` pair; owner decision: run fidelity 2 at a fixed
   `L0 = 25 m`). That is backlog 2-P5, HIGH.
1. **Against L0 = inf the missing power is ALL tip and tilt, at every count.** `Delta3` (the
   tilt-removed variance) reads 1.00 +-0.01 in every configuration, and
   `Delta1` (the piston-removed variance) reads 0.75 to 0.84. So the generator
   holds the high orders exactly and misses 16 to 25 percent of the aperture
   phase variance, all of it in the tip-tilt band. This is the 2-N2 deficit of
   `validation/screens/` measured again on the production grid, and the two
   agree (0.75 to 0.78 there).
2. **Stacking barely changes it.** At L0 = inf: the ground layer as 1 or 4
   screens gives `Delta1` 0.823 -> 0.839 (0.3 sigma) and the whole plan at
   5 -> 25 screens 0.772 -> 0.748 (0.7 sigma). At L0 = 25 m, where the offset
   of item 0 is gone, a MILD real drift shows: the ground layer as 1 or 4
   screens is a null (1.000 -> 1.004, 0.1 sigma), but the whole plan reads
   0.974 / 0.946 / 0.942 at 5 / 9 / 25 screens against 1.000 for one screen,
   about 5 percent low at 1 to 2 sigma, with D(r) 3 to 6 percent low at
   r >= 0.35 m for the 9- and 25-screen stacks. So the owner's hypothesis
   holds at the FEW-PERCENT level for a many-screen plan and not at all for
   the ground split. Five percent of the tilt variance is a fraction of a dB
   on an SMF fade, against the 2 dB p5 trend of the tail-convergence study,
   so the count trend of that study still comes from the PROPAGATION and
   placement side of the product, not from the screens.
3. **The structure function at large separation.** At L0 = inf `D(1.4 m)`
   reads 0.59 to 0.72 of Fried: that is item 0 again (the L0 = 95 m theory
   gives 0.64). At L0 = 25 m it reads 0.86 to 0.96, the 25-screen stack the
   lowest. A separation of 1.4 m is 0.4 of the grid side, where a periodic
   FFT screen is at its weakest, so read D(1.4 m) as a bound.
4. **What matters is item 0.** Every fidelity-2 SMF fade today runs on the
   tilt of an effective 95 m outer scale while the code claims infinity. The
   fibre overlap pays the tilt, so the fidelity-2 SMF tail depends on a
   physics choice (the site L0) that no input sets. The work is backlog 2-P5:
   an explicit site L0 threaded to the screens and to the analytic tilt Terms,
   and a sizer check that raises the subharmonic depth (a factor 3 a level)
   when L0 exceeds 3^n times the grid side, which bites for a SMALL aperture.
   The known screen limits stay (2-N2, S-27; `validation/screens/FINDINGS.md`).

## Files

| File | Purpose |
| --- | --- |
| [screen_stacking.py](screen_stacking.py) | The test. `--L0 <m>` draws the screens with a von Karman outer scale and judges them against the von Karman theory at that L0; the default is infinite (Fried / Noll). It writes `screen_stacking[_L0<m>]_results.json` and the matching `.log`. |
| [screen_stacking_results.json](screen_stacking_results.json) | The L0 = inf ratios with their standard errors, and the verdict. |
| [screen_stacking_L025_results.json](screen_stacking_L025_results.json) | The same at L0 = 25 m. |
| [screen_stacking.log](screen_stacking.log), [screen_stacking_L025.log](screen_stacking_L025.log) | The printed tables. |

This study is VALIDATION ONLY. It reads the production generator and it
changes no `olb` module.
