# Phase-screen low-frequency study: findings

## What this study is

A phase screen on a finite grid holds no power below its grid fundamental.
That missing low-frequency band is the tip and the tilt. This study measures
how much tilt each screen route holds, and it compares each route against the
analytic value. The study answers six questions. It also settles two open
rows in `docs/schmidt-crosscheck.md`, and it gives a design recommendation for
`olb/waveoptics/turbulence/temporal.py`.

The study is VALIDATION ONLY. It reads the production layer. It changes no olb
module.

Five scripts live in `validation/screens/`:

- `helpers.py` — the analytic truths and the shared estimators.
- `oversize_crop.py` — arm 1, the Fourier and the oversize-and-crop screens.
- `infinite_screen_stats.py` — arm 2, the spatial statistics of the extruded
  screens.
- `extrusion_stationarity.py` — arm 3, the drift test of the extrusion.
- `n_columns_sweep.py` — arm 4, the `n_columns` sweep of the extrusion (Q6).

Common parameters: lambda = 1550 nm, r0 = 0.10 m, dx = 0.01 m, reference grid
512 pixels (side 5.12 m), pupil D = 1.0 m (D/r0 = 10), von Karman L0 = 25 m
where the spectrum is finite. The first run of the screens came from aotools
1.0.7. The rerun of 2026-09-08 came from aotools
`0.1.dev477+g8f38ff7c1.d20260605` (dist-info 1.0.8). Read the rerun section
before you use a number of the extruded screens.

The analytic truths are the Noll per-axis Z-tilt filter integral (Noll,
DOI 10.1364/JOSA.66.000207), the Andrews G-tilt filter (Andrews and Phillips,
DOI 10.1117/3.626196, Ch. 6, Eq. (84), printed p. 201), the von Karman
covariance in a float64 closed form (Assemat and Wilson,
DOI 10.1364/OE.14.000988, Eq. (5)) and as a Hankel transform of the phase power
spectral density (Schmidt, DOI 10.1117/3.866274, Ch. 9), and the structure
function D(r) = 6.88 (r/r0)^(5/3) (Schmidt, DOI 10.1117/3.866274, Ch. 9,
printed p. 160).

## How to reproduce

Run each command from the repository root. The first four commands take about
15 minutes in total. The fifth command is the Q6 sweep, and it takes about 15
minutes on its own.

```
python -m validation.screens.helpers
python -m validation.screens.oversize_crop
python -m validation.screens.infinite_screen_stats
python -m validation.screens.extrusion_stationarity
python -m validation.screens.n_columns_sweep
python -m validation.screens.subharmonic_start
python -m validation.screens.kept_start
```

The sixth command is the Q7 start check, about 4.5 minutes. The seventh is
the Q8 kept-start check, about 29 minutes at 32 seeds.

## Q1 — Does oversize-and-crop beat subharmonic augmentation for tip/tilt?

**Answer: no for a pure Kolmogorov spectrum, and yes for a von Karman spectrum
with a finite outer scale.**

The metric is the Z-tilt angle variance over the 1.0 m pupil, divided by the
Noll filter integral (DOI 10.1364/JOSA.66.000207). The error bars are 2
standard errors.

Kolmogorov (`oversize_tilt.csv`):

| arm | Z-tilt ratio | 2 SE |
|---|---|---|
| book plain Fourier, K = 1 | 0.406 | 0.081 |
| aotools three-level subharmonic, K = 1 | 0.752 | 0.150 |
| book three-level subharmonic, K = 1 | 0.782 | 0.156 |
| oversize-and-crop, K = 2 | 0.547 | 0.109 |
| oversize-and-crop, K = 4 | 0.532 | 0.137 |
| oversize-and-crop, K = 8 | 0.840 | 0.266 |

The K = 8 oversize arm reaches the subharmonic level, but its error bar is
wide. The sharp-cutoff capture model predicts 0.773 for a subharmonic reach of
1/(27 x 5.12 m). So the two cures buy the same band, and neither one wins.

The structure function agrees (`oversize_dphi.csv`). Across r/r0 = 0.5 to 20
the K = 8 arm falls from 0.912 to 0.707 of the theory, and the book subharmonic
arm falls from 0.908 to 0.709. The two curves track each other at every
separation.

Von Karman, L0 = 25 m (`oversize_tilt.csv`):

| arm | Z-tilt ratio | 2 SE |
|---|---|---|
| book plain Fourier, K = 1 | 0.879 | 0.227 |
| aotools three-level subharmonic, K = 1 | 0.915 | 0.236 |
| book three-level subharmonic, K = 1 | 1.035 | 0.267 |
| oversize-and-crop, K = 2 | 0.939 | 0.242 |
| oversize-and-crop, K = 4 | 0.848 | 0.268 |
| oversize-and-crop, K = 8 | 1.050 | 0.383 |

The K = 2 arm (side 10.24 m) already sits at the analytic value inside the
error bar. The K = 8 arm holds D(r)/theory between 0.995 and 1.007 for r/r0
between 1 and 10. The subharmonic arms also reach about 1.0 there. So an
oversize screen is a full cure for a finite outer scale, and the subharmonic
route is an equal cure at a lower cost.

## Q2 — At what oversize factor does the tilt converge?

**Answer: it converges from K = 2 with a finite outer scale, and it never
converges at a practical factor without one.**

Von Karman, L0 = 25 m: the tilt sits inside the noise of the analytic value
from K = 2, which is a screen side of 10.24 m. The convergence is complete at a
side of about 41 m. State the rule in outer scales: the screen side must reach
about L0/2 for a first agreement, and about 1.6 L0 for a full agreement.

Pure Kolmogorov: the captured share grows as about K^(1/3). The sharp-cutoff
model predicts 0.32, 0.46, 0.57 and 0.66 for K = 1, 2, 4 and 8. The pupil tilt
variance of a Kolmogorov spectrum is finite, but the screen side that captures
it grows without bound.

The practical rule for olb: select the outer scale from the physics first.
Then size the screen side to at least one to two outer scales.

## Q3 — Is the infinite-screen covariance wrong?

**Answer: the formula is correct, the screen variance is spin-up limited, and
the extrusion axis carries a real defect.**

The outside claim mixes three levels. The study separates them.

### Level 1 — the formula: CORRECT

`aotools.turbulence.turb.phase_covariance` implements Eq. (5) of Assemat and
Wilson (DOI 10.1364/OE.14.000988). It matches the float64 closed form to
3.9e-7 at L0 = 25 m, and to 1.2e-6 at L0 = 2.56 m. The float32 cast and the
1e-40 offset are harmless. The closed route and the Hankel route differ by a
flat 0.454 percent. That offset is the rounded printed constant 0.023 of the
von Karman power spectral density (Schmidt, DOI 10.1117/3.866274, Ch. 9).

### Level 2 — the screen variance: NOT wrong, but spin-up limited

At L0 = 2.56 m a 512-row spin-up covers 2.0 outer scales. The raw variance then
lands to 0.08 percent of the theory (`infinite_covariance.csv`: B(0) = 19.178
against 19.193 rad^2). The transverse-axis covariance holds inside 2.4 percent
of B(0) at every lag.

At L0 = 25 m the same 512 rows cover only 0.2 of one outer scale. The row-lag
record then shows a flat 12 percent deficit at every lag
(`rowlag_covariance.csv`: C(0) = 756.3 against a theory of 860.3 rad^2). A
single 512-row frame is worse, because it also loses the frame piston
(`infinite_covariance.csv`: B(0) = 322.5 against 856.3 rad^2). A spin-up ladder
over 128, 512 and 2048 rows shows the piston at 0.29 of the theory at 512 rows,
and at 0.72 at 2048 rows. The piston still climbs there.

The deficit is missing outer-scale power. It is not a wrong kernel.

### Level 3 — the extrusion axis: GENUINELY DEFECTIVE

(REVISED in Q6.)

The 2-column Markov recursion (`n_columns=2`) over-correlates its own
direction. Three measurements show it:

- At L0 = 2.56 m the extrusion axis holds 1.498 rad^2 at a 2.56 m lag. The
  theory gives 0.085 rad^2, and the transverse axis reads -0.104 rad^2
  (`infinite_covariance.csv`).
- The D(r) anisotropy reaches 20 to 30 percent of the theory at separations of
  1.0 to 2.5 m. The pooled trend statistic over 16 runs gives |t| = 3.07 and
  |t| = 2.58.
- The normalised row-lag correlation rho(k) carries a growing excess:
  +0.002 at 0.64 m, +0.023 at 1.28 m, and +0.095 at 2.56 m
  (`rowlag_covariance.csv`).

This bias is stationary. It does not drift. It is the substance behind the
outside advice. So the advice is right about the extrusion, and wrong about the
formula.

## Q4 — Does the extrusion drift after the spin-up?

**Answer: no. The extrusion is stationary after the spin-up.**

The test fits a straight line to each metric against the window index over
windows 1 to 4, which are 512 to 4096 cumulative rows. All five trend tests
pass at |t| < 1.1. The script prints `stationary after spin-up: yes`.

Window 0 is the initial frame. The aotools source confirms why it sits low:
`make_initial_screen` calls `ft_phase_screen`, which is the plain Fourier
screen with no subharmonics. The step from window 0 to window 1 lifts D(r) by
50 percent on average, and it lifts the Z-tilt variance by 24 percent
(`stationarity_windows.csv`). That step is the spin-up, and the study expects
it.

## Q5 — The recommendation for `olb/waveoptics/turbulence/temporal.py`

The module is a `NotImplementedError` stub today. Its plan reads
`PhaseScreenVonKarman` plus `add_row`. The study gives five points.

1. **A finite outer scale is mandatory.** The class requires it. A pure
   Kolmogorov tilt never converges on any finite screen (Q2).

2. **Specify the spin-up in outer scales, not in rows.** Discard at least
   2 L0 / dx rows before the first frame. At L0 = 25 m and dx = 0.01 m that is
   5000 rows. A shorter spin-up carries a variance deficit of 10 percent or
   more (Q3, level 2).

3. **The extrusion-axis over-correlation is a real defect for a temporal
   model.** Under frozen flow a row lag is a time lag (Taylor,
   DOI 10.1098/rspa.1938.0032). So the defect smooths the TEMPORAL axis
   exactly. A temporal spectrum reads too slow, and a fade duration reads too
   benign. Raise `n_columns` to reduce it. That trade of memory against
   accuracy is unverified here.

4. **Prefer the shifted large screen.** This study validates the alternative.
   Pre-compute one large Fourier screen, and shift a crop window across it for
   each frame. At L0 = 25 m a K = 8 screen (side 41 m, 4096 by 4096) reproduces
   the tilt and D(r) to the theory inside the noise (Q1). It carries no axis
   anisotropy by construction. It costs one fast Fourier transform at the start
   and one slice per frame. It supports a subpixel shift through Fourier
   interpolation. A wraparound reuse must respect one full screen traverse. The
   memory cost is 134 MB per layer for a 4096 by 4096 float64 array. Weigh that
   cost against the extrusion defect.

   **Recommendation: make the shifted large screen the default design for
   `temporal.py`.** Keep the extruded class as a fallback for a memory limit
   only. Use it with a raised `n_columns` and with a measured re-check.

5. **The snapshot layer needs no change.** The subharmonic route stays correct
   for snapshot statistics. Production runs `subharmonics=True` today, and that
   route matches the K = 8 oversize screen on every measured metric (Q1). This
   study indicates no change to `olb/waveoptics/turbulence/screens.py`.

**Addendum, 2026-09-08. The fixed-build numbers of points 3 and 4.** Q6
measured `n_columns` from 2 to 32. Point 3 stays true, but its cure changes.
The control of the over-correlation is the FRAME WIDTH in outer scales, not
`n_columns`. At a frame side of 2 L0 the extrusion axis reads the theory at
`n_columns = 2`. At a frame side of 0.2 L0 the excess is +0.228 on rho at
`n_columns = 2`, and `n_columns = 8` reduces it only to +0.177. So a raise of
`n_columns` reduces the defect by 1.3 to 3.8 times, and it never removes it.
The production grid is 1024 pixels at about 5 mm, which is a side of about
5 m. At L0 = 25 m that side is 0.2 L0, and it is the SAME regime as the
N = 512 / L0 = 25 m cell of Q6. So an extruded screen on the production grid
keeps a large extrusion-axis over-correlation at every practical `n_columns`.
Point 4 keeps its recommendation, and Q6 makes it stronger. The shifted large
screen has a side of 1.6 L0 or more by construction, so it sits in the regime
where the statistics hold. Make the shifted large screen the default design
for `temporal.py`.

## Side findings

**S-27 settlement.** At N = 512, with shared draws and with both estimators,
the aotools `ft_sh_phase_screen` reads ABOVE the book generator. The excess is
+4.4 to +9.5 percent for the fast Fourier transform estimator, and +7.1 to
+15.0 percent for the direct estimator (`s27_settlement.csv`). At r/r0 = 8 the
direct estimator reads 0.707 for the book generator and 0.813 for aotools. So
aotools is CLOSER to the theory there. This result supports the numbers in the
`screens.py` docstring. It CONTRADICTS the gap S-27 row and the forward-map row
in `docs/schmidt-crosscheck.md`, whose run used N = 256. The difference is a
grid-size effect. The owner must correct those two rows. This study does not
edit them.

**The aotools shared-seed quirk.** `ft_sh_phase_screen` reuses the integer seed
for the subharmonic draws and for the high-frequency draws. So the 9
low-frequency Gaussians duplicate the first 9 high-frequency draws. The
measured tilt correlation is 0.087 +- 0.158 over 80 samples, and the variance
shift stays well inside 2 standard errors (`seed_quirk.csv`). The quirk is
cosmetic. It gives no bias.

**The printed G-tilt constant is rounded.** The book prints 0.174 (Andrews and
Phillips, DOI 10.1117/3.626196, Ch. 6, Eq. (84), printed p. 201). The filter
integral gives 0.1698, which matches Sasiela
(DOI 10.1007/978-3-642-59022-0). The printed structure constant 6.88 is 6.8839
exactly.

**The capture model under-estimates a plain Fourier screen.** A sharp-cutoff
capture model assumes that the screen holds nothing below its grid fundamental.
The measured Z-tilt variance is 1.15 times that prediction on average over the
eight crop arms (range 0.89 to 1.47; the plain Kolmogorov Fourier screen at
K = 1 reads 1.27 times). The cause is the cells next to the direct-current cell,
which lump the sub-fundamental band. So `helpers.captured_fraction` is a slight
under-estimate. One pass band failed for exactly that reason: the book plain
Fourier Kolmogorov arm reads 0.406 against a band top of 0.400.

**Read the structure-function comparison inside 80 percent of the window.** The
fast Fourier transform estimator makes the correlation circular. It collapses
in its last bin, where the separation equals the mask diameter. The
`infinite_dphi.csv` row at r = 2.560 m shows that collapse (ratio 0.344 against
about 0.75 at the next separation).

## Rerun on the fixed kernel (2026-09-08)

### Why the study ran again

The first run used aotools 1.0.7. That build cast the separations of
`turbulence.turb.phase_covariance` to float32. The stencil covariance
`cov_zz` has a condition number of 7e7 at N = 512, dx = 1 cm and L0 = 25 m.
Float32 has an epsilon of 6e-8. So the inverse was noise-limited, and the A
matrix of the extrusion was noise-limited with it. Upstream fixed the cast to
float64 in pull request 111 (2026-04-02). Two upstream issues describe the
symptom: issue 107 (the noise grows on one side of the screen after many
`add_row` calls) and issue 109 (the variance, so r0, is not preserved).

- https://github.com/AOtools/aotools/pull/111
- https://github.com/AOtools/aotools/issues/107
- https://github.com/AOtools/aotools/issues/109

The installed build carries the fix. The four commands of "How to reproduce"
ran again on it. All four completed. The 1.0.7 tables stay under
`data/*_aotools107.csv` for the comparison.

### What did not move

Four tables are BIT-IDENTICAL to the 1.0.7 run: `oversize_tilt.csv`,
`oversize_dphi.csv`, `s27_settlement.csv` and `seed_quirk.csv`. Those tables do
not read the extruded screen. So Q1, Q2, the S-27 settlement and the seed
quirk stand without a change. Window 0 of `stationarity_windows.csv` is also
identical, because it is the plain Fourier start.

### What moved

| quantity | aotools 1.0.7 | fixed build | theory |
|---|---|---|---|
| extrusion-axis B at a 2.56 m lag, L0 = 2.56 m, rad^2 | 1.498 | 1.509 | 0.085 |
| transverse-axis B at the same lag, rad^2 | -0.104 | -0.104 | 0.085 |
| D(r) anisotropy (Dx - Dy) / D_theory at 1.00 m | +0.20, \|t\| = 3.07 | +0.195, \|t\| = 3.96 | 0 |
| D(r) anisotropy (Dx - Dy) / D_theory at 2.50 m | +0.30, \|t\| = 2.58 | +0.257, \|t\| = 2.90 | 0 |
| row-lag rho excess at 0.64 m | +0.002 | +0.008 | 0 |
| row-lag rho excess at 1.28 m | +0.023 | +0.037 | 0 |
| row-lag rho excess at 2.56 m | +0.095 | +0.121 | 0 |
| row-lag C(0), rad^2 | 756.3 | 783.7 | 860.3 |
| screen variance B(0), L0 = 2.56 m, rad^2 | 19.178 | 19.148 | 19.193 |
| single-frame B(0), L0 = 25 m, rad^2 | 322.5 | 305.0 | 856.3 |
| extruded-screen Z-tilt ratio | 0.916 | 0.946 | 1.000 |
| extruded-screen G-tilt ratio | 0.904 | 0.961 | 1.000 |
| spin-up step of D(r), window 0 to 1 | +50 percent | +43 percent | - |
| spin-up step of the Z-tilt variance, window 0 to 1 | +24 percent | -1.5 percent | - |

### The verdict on Q3 Level 3

**The Q3 Level 3 result SURVIVES. The extrusion axis stays GENUINELY
DEFECTIVE.** The fixed kernel makes the defect slightly LARGER, not smaller.
All three measurements hold:

1. The extrusion axis holds 1.509 rad^2 at a 2.56 m lag against a theory of
   0.085 rad^2. The transverse axis reads -0.104 rad^2, which is the same
   value as the 1.0.7 run.
2. The D(r) anisotropy stays at 20 to 26 percent at 1.0 to 2.5 m. Both trend
   statistics grew: \|t\| = 3.96 and \|t\| = 2.90.
3. The row-lag correlation excess grew at every lag. It reads +0.121 at
   2.56 m, against +0.095 before.

So the float32 kernel was not the cause of the over-correlation. The Q5
point 3 advice keeps its evidence, and the Q5 point 4 recommendation stands.

### The verdict on Q4

**The Q4 result SURVIVES. The extrusion is stationary after the spin-up.** All
five trend tests pass, and the largest statistic is \|t\| = 0.58. That is
better than the \|t\| < 1.1 of the 1.0.7 run. The script prints
`stationary after spin-up: yes`.

One number of Q4 changed sign. The spin-up step of the Z-tilt variance reads
-1.5 percent, where 1.0.7 read +24 percent. That metric is noisy: its 2
standard errors are about 0.25 on a ratio near 1.0, and the window-to-window
scatter of the ratio runs from 0.60 to 1.08. So read the D(r) step, which
stays a clear lift of +43 percent. Do not read the Z-tilt step as a trend.

### The gain of the fix

The fixed kernel moves three quantities toward the theory. The row-lag C(0)
deficit falls from 12.1 to 8.9 percent. The Z-tilt ratio rises from 0.916 to
0.946. The G-tilt ratio rises from 0.904 to 0.961. So the fix corrects a
VARIANCE deficit. It does not correct the extrusion-axis MEMORY.

## Q6 — What does `n_columns` buy? (2026-09-08)

**Answer: it buys very little. The frame width in outer scales controls the
extrusion-axis error, not the stencil depth.**

### The method

`validation/screens/n_columns_sweep.py` extrudes
`PhaseScreenVonKarman` at `n_columns` in {2, 4, 8, 16, 32}, at L0 = 2.56 m and
25 m, on a 128-pixel grid, plus a 512-pixel grid at `n_columns` in {2, 8}. Each
cell runs 8 seeds, it discards 2 L0 / dx rows of spin-up, and it records 20
L0 / dx rows (10000 rows at N = 512, L0 = 25 m). Every metric compares against
the closed-form von Karman covariance of Assemat and Wilson
(DOI 10.1364/OE.14.000988, Eq. (5)), with no mean subtraction.

### The band table

The frame side is N dx. The column `side / L0` is that side in outer scales.
Metric 1 is the extrusion-axis rho(k), band 0.02 past 0.5 L0. Metric 2 is the
transverse rho(k), the control, same band. Metric 3 is the structure-function
anisotropy D_ext / D_trans, band 5 percent. Metric 4 is the chunk-mean variance
of the row piston, band 30 percent.

| N | L0, m | side / L0 | n_col | 1. rho ext | 2. rho trans | 3. anisotropy | 4. piston |
|---|---|---|---|---|---|---|---|
| 128 | 2.56 | 0.50 | 2 | FAIL +0.102 | FAIL +0.039 | FAIL 16.6 % | PASS 1.146 |
| 128 | 2.56 | 0.50 | 4 | FAIL +0.102 | FAIL +0.039 | FAIL 14.4 % | PASS 1.149 |
| 128 | 2.56 | 0.50 | 8 | FAIL +0.097 | FAIL +0.039 | FAIL 14.2 % | PASS 1.068 |
| 128 | 2.56 | 0.50 | 16 | FAIL +0.094 | FAIL +0.039 | FAIL 14.2 % | PASS 1.023 |
| 128 | 2.56 | 0.50 | 32 | FAIL +0.093 | FAIL +0.040 | FAIL 14.3 % | PASS 1.005 |
| 128 | 25.0 | 0.051 | 2 | FAIL +0.297 | pass +0.017 | FAIL 112.8 % | FAIL 1.754 |
| 128 | 25.0 | 0.051 | 4 | FAIL +0.268 | pass +0.014 | FAIL 68.2 % | FAIL 1.696 |
| 128 | 25.0 | 0.051 | 8 | FAIL +0.206 | pass -0.006 | FAIL 55.7 % | FAIL 1.517 |
| 128 | 25.0 | 0.051 | 16 | FAIL +0.135 | pass -0.007 | FAIL 59.9 % | FAIL 1.384 |
| 128 | 25.0 | 0.051 | 32 | FAIL +0.078 | pass -0.008 | FAIL 61.8 % | PASS 1.268 |
| 512 | 2.56 | 2.00 | 2 | FAIL +0.042 | pass +0.005 | FAIL 6.7 % | PASS 0.993 |
| 512 | 2.56 | 2.00 | 8 | FAIL -0.037 | pass +0.003 | PASS 3.4 % | PASS 0.959 |
| 512 | 25.0 | 0.205 | 2 | FAIL +0.228 | FAIL -0.072 | FAIL 28.5 % | PASS 1.079 |
| 512 | 25.0 | 0.205 | 8 | FAIL +0.177 | FAIL -0.063 | FAIL 15.1 % | PASS 0.966 |

21 of 61 bands pass. The guard band passes at every cell: no `lstsq` fallback
runs, and the minimum eigenvalue of BB^T holds at 1.117e-2 everywhere. So
errors E3 and E4 of the plan stay untriggered.

Read the metric 1 column against its error bar. The 2 standard errors of the
worst lag run from 0.019 to 0.168. At L0 = 2.56 m the worst deviation is 1.5 to
2.5 standard errors, and it is a small number. At L0 = 25 m the deviation is a
large number at every `n_columns`.

Two columns need a caution. Metric 2 fails at L0 = 2.56 m on the 128-pixel
grid, and it fails on the 512-pixel grid at L0 = 25 m. In both cases the
largest transverse lag equals half the grid side, where the estimator reads few
pairs. That failure is a grid-width effect, not a recursion effect, because the
recursion does not touch the transverse axis. Metric 3 at L0 = 25 m divides two
very small structure functions (D of order 0.002 to 0.09 against a variance of
856 rad^2), so its ratio is noisy at small r.

### The reading: frame width, not stencil depth

Sort the cells by `side / L0`, and the metric 1 column falls into order:

| side / L0 | worst rho excess, n_col 2 | worst rho excess, best n_col |
|---|---|---|
| 2.00 (N 512, L0 2.56 m) | +0.042 | -0.037 at 8 |
| 0.50 (N 128, L0 2.56 m) | +0.102 | +0.093 at 32 |
| 0.205 (N 512, L0 25 m) | +0.228 | +0.177 at 8 |
| 0.051 (N 128, L0 25 m) | +0.297 | +0.078 at 32 |

The stencil depth is the weak lever. At L0 = 2.56 m a raise from 2 to 32
columns moves the excess by 0.009 only. At L0 = 25 m it moves it by 2 to 4
times, and it still leaves 4 to 40 times the band. The frame width is the
strong lever: the same `n_columns = 2` reads +0.042 at 2 L0 of frame and +0.297
at 0.05 L0 of frame. The row piston gives the same reading. Its chunk-mean
variance is 1.75 times the theory at 25 m and `n_columns = 2` on the 128-pixel
grid, and it is still 1.27 times at `n_columns = 32`; on the 512-pixel grid it
reads 1.08 at `n_columns = 2`. So a wider frame cures the piston, and the
stencil depth only reduces it.

### The Phase 2 window verdict

Phase 2 asked one question: is the excess a short-record artefact? It is NOT.
The script cut one 100 L0 record into windows of 2, 5, 10, 20 and 100 outer
scales, and it read rho at a lag of one outer scale in each window.

- At L0 = 2.56 m the apparent rho holds the theory (0.0044) at EVERY window.
  It runs -0.020 to +0.012 at `n_columns = 2`, and -0.029 to -0.003 at
  `n_columns = 16`. Both longest windows pass the band.
- At L0 = 25 m the excess GROWS with the window. At `n_columns = 2` it climbs
  from +0.092 at 2 L0 to +0.227 at 100 L0. At `n_columns = 16` it climbs from
  +0.022 to +0.079. Both longest windows fail the band.

An artefact of a short record shrinks as the record grows. This excess grows.
So it is a real slow mode of the recursion, and a longer record does not cure
it. The 2026-09-08 probe that read the theory on a 16384-row record was a
sampling accident.

### The cost

| n_columns | cov_zz size | cov_zz, MB | setup, s (N 128) | add_row, ms (N 128) | setup, s (N 512) | add_row, ms (N 512) |
|---|---|---|---|---|---|---|
| 2 | 2 N | 0.5 / 8.0 | 0.03 to 0.64 | 0.013 to 0.015 | 0.60 to 0.77 | 0.91 |
| 4 | 4 N | 2.0 | 0.10 to 0.12 | 0.017 | not run | not run |
| 8 | 8 N | 8.0 / 128 | 0.27 to 0.33 | 0.023 | 5.18 to 6.93 | 1.11 |
| 16 | 16 N | 32.0 | 0.96 to 1.23 | 0.038 to 0.041 | not run | not run |
| 32 | 32 N | 128.0 | 4.10 to 4.93 | 0.097 to 0.104 | not run, 2.1 GB | not run |

`cov_zz` holds (`n_columns` N)^2 doubles, so it grows as the square of the
stencil depth. At `n_columns = 32` and N = 512 it would be 16384 by 16384, or
2.1 GB, and its Cholesky factor would cost minutes. The run time of one
`add_row` is the mild cost: at N = 512 it grows from 0.91 to 1.11 ms between
`n_columns` 2 and 8, which is 22 percent. The setup is the hard cost: it grows
from 0.6 to 6.9 s over the same step, and it grows as the cube of the stencil
depth.

### The reconciliation with Q3 Level 3

Q3 Level 3 reports, at L0 = 2.56 m and N = 512 with `n_columns = 2`, an
extrusion-axis covariance of 1.509 rad^2 at a 2.56 m lag, against a theory of
0.085 rad^2, with the transverse axis at -0.104 rad^2. The sweep reads the SAME
cell at rho = -0.001 against a theory of 0.005 at a 2.49 m lag
(`ncol_rowlag.csv`). A direct check settled the difference (6 seeds, N = 512,
L0 = 2.56 m, `n_columns = 2`; one record per seed, both estimators applied to
the SAME rows):

| estimator | rows | B at a 2.56 m lag, rad^2 | 2 SE |
|---|---|---|---|
| `infinite_screen_stats.axis_covariance`, extrusion | the 512-row spin-up frame | +1.288 | 3.546 |
| `extrusion_stationarity.row_lag_sums` | the SAME 512 rows | +1.288 | 3.546 |
| `axis_covariance`, transverse | the same frame | -0.761 | 3.646 |
| `axis_covariance`, extrusion | a late frame, no transient rows | -1.403 | 2.003 |
| `row_lag_sums`, rho | the 5120-row settled record | +0.0016 | 0.0464 |

**NEITHER estimator carries a bias. The two agree to every printed digit on the
same rows.** The two numbers differ because of the SAMPLE, not the estimator. A
lag of 256 rows on a 512-row frame averages 256 row pairs of ONE realisation,
so one frame gives a standard deviation of about 4.3 rad^2 at that lag. Sixteen
frames then give 2 standard errors of about 2.2 rad^2. The reported 1.509 rad^2
sits 0.65 of that distance from the theory, and the transverse -0.104 rad^2 is
the same noise. `infinite_covariance.csv` carries no error column, so the study
read that noise as a signal. The long record gives 5120 rows and a 2 SE of
0.046 on rho, which is 100 times tighter, and it reads the theory.

**So the L0 = 2.56 m part of Q3 Level 3 is an ARTEFACT of the sample size. It
is not a defect of the extrusion.** The other two measurements of Q3 Level 3
SURVIVE, and the sweep reproduces them, but both came from L0 = 25 m on a
512-pixel grid. That is a frame side of 0.205 L0. The D(r) anisotropy of
`extrusion_stationarity` (+0.195 at 1.0 m, |t| = 3.96; +0.257 at 2.5 m,
|t| = 2.90) and the row-lag rho excess of `rowlag_covariance.csv` (+0.121 at
2.56 m) both sit in that cell, and the sweep reads 28.5 percent of anisotropy
and +0.228 of rho excess there. **The defect is REAL, and it is the
frame-width-against-outer-scale effect. It shows at L0 = 25 m, and it does not
show at L0 = 2.56 m on a frame of 2 outer scales.**

### The verdict on the smallest acceptable `n_columns`

Select `n_columns` from the frame side in outer scales, not from a fixed rule.
At a frame side of 2 L0 or more, `n_columns = 2` is enough: the extrusion axis
holds inside 0.042 on rho, the piston reads 0.99 of the theory, and the
anisotropy holds inside 7 percent, so the paper's own claim of two rows is
correct there (Assemat and Wilson, DOI 10.1364/OE.14.000988). Take
`n_columns = 8` if the anisotropy must hold inside 5 percent; it costs 128 MB
of `cov_zz` and 7 s of setup at N = 512, and 22 percent of the `add_row` time.
At a frame side near 0.5 L0, `n_columns` buys nothing on the extrusion axis
(+0.102 at 2 columns against +0.093 at 32), so keep 2 columns and accept a rho
excess of about 0.10; raise it to 16 only when the row piston matters, because
that raise moves the piston from 1.15 to 1.02 of the theory. At a frame side of
0.2 L0 or less — the production regime — NO `n_columns` is acceptable. The best
measured cell still fails by four times the band, the cost of the next step is
2.1 GB of `cov_zz`, and the correct cure is a wider frame or the shifted large
screen of Q5 point 4.

## Q7 — Does a subharmonic initial frame remove the spin-up? (2026-09-08)

**Answer: yes for every quantity the optics sees. The residual "deficit" of a
subharmonic start is the frame piston, which the generators zero on purpose
and which does nothing to a propagation.**

The stock `make_initial_screen` starts the extrusion from a plain Fourier
frame (`ft_phase_screen`, no subharmonics). The recursion is Markov in the
last `n_columns` rows (Assemat and Wilson, DOI 10.1364/OE.14.000988), so the
whole future depends on the initial frame only through those rows. Arm 5
(`subharmonic_start.py`) subclasses the class and swaps the initial frame for
a three-level subharmonic screen (Lane, Glindemann and Dainty,
DOI 10.1088/0959-7174/2/3/003), from the olb production generator
(`olb.waveoptics.turbulence.screens.phase_screen`) and from aotools
(`ft_sh_phase_screen`). N = 512, dx = 1 cm, r0 = 10 cm, L0 = 25 m,
`n_columns = 2`, 32 seeds, frames read at row 0, 512 and 2048.

| start | rows | piston-removed variance / theory | Z-tilt over 1 m / Noll |
|---|---|---|---|
| plain | 0 | 0.317 | 0.713 |
| plain | 512 | 0.587 | 0.840 |
| plain | 2048 | 0.854 | 1.074 |
| olb subharmonic | 0 | 0.895 | 1.084 |
| olb subharmonic | 512 | 0.935 | 1.055 |
| olb subharmonic | 2048 | 0.725 | 0.952 |
| aotools subharmonic | 0 | 1.029 | 1.044 |
| aotools subharmonic | 512 | 0.740 | 0.874 |
| aotools subharmonic | 2048 | 0.970 | 0.876 |

The Z-tilt standard error is about 18 percent. The band is 0.20.

Three readings:

1. **The optical spin-up is gone at row 0.** Both subharmonic starts read the
   piston-removed variance and the 1 m Z-tilt inside the band on the first
   frame, against 0.32 and 0.71 for the plain start. So a worker can use row
   1 onward, with no discarded rows.
2. **The frame piston is zero by construction.** `ft_sh_phase_screen` and
   `ScreenFactory` both subtract the mean of the low-frequency part. So the
   raw frame variance of a subharmonic start reads 0.27 to 0.31 of B(0), which
   is B(0) minus the piston theory (856 minus 598 rad^2). That is a report, not
   a defect: a frame-wide constant phase does nothing to a propagation.
3. **The 2 L0 spin-up rule of Q5 point 2 was set by the piston, not by the
   optics.** With the plain start the piston-removed variance and the Z-tilt
   reach 0.85 and 1.07 at 2048 rows, which is 0.2 L0 at 25 m. So even the
   plain start needs about 0.2 L0 / dx rows for the optical quantities, ten
   times fewer than the rule, and the subharmonic start needs none.

**Not settled: whether the recursion KEEPS the correct start on a narrow
frame.** The subharmonic rows at 512 and 2048 scatter between 0.73 and 0.97
on the variance and 0.87 and 1.06 on the tilt, and one of the four "kept"
bands fails (olb at 2048 rows). Q6 says the recursion has its own stationary
level on a 0.2 L0 frame, so a drift from the correct start toward that level
is expected, and 32 seeds do not resolve its rate. A long record with the
subharmonic start, read with the Q6 row-lag estimator, is the test.

## Q8 — Is the subharmonic start kept on a 0.2 L0 frame? (2026-09-08)

**Answer: the two starts are indistinguishable on the extrusion axis, and
both carry the Q6 over-correlation in full. The start does not set the slow
modes of the record; the recursion on the narrow frame does.**

Arm 6 (`kept_start.py`) runs the Q6 production cell (N = 512, L0 = 25 m, the
frame is 0.2 L0) with NO spin-up from both starts, `n_columns` 2 and 8, 32
seeds, 10000 rows (4 L0) each, and the Q6 row-lag estimator over the whole
record. The error bar is 2 standard errors over the seeds. A first run with
8 seeds and no error bar (2026-09-08, `kept_start.log`) read the subharmonic
start ABOVE the plain start by 0.19 to 0.20; the 32-seed run shows that was
noise.

| `n_columns` | start | rho at 0.5 L0 | excess | 2 SE | Q6 plain, spun up |
|---|---|---|---|---|---|
| 2 | plain | 0.355 | +0.271 | 0.085 | +0.228 |
| 2 | olb subharmonic | 0.341 | +0.257 | 0.083 | +0.228 |
| 8 | plain | 0.262 | +0.178 | 0.078 | +0.177 |
| 8 | olb subharmonic | 0.272 | +0.189 | 0.078 | +0.177 |

The theory is 0.084. The difference between the starts, subharmonic minus
plain, with its 2 SE:

| `n_columns` | lag 0.1 L0 | lag 0.5 L0 | lag 1 L0 |
|---|---|---|---|
| 2 | -0.021 +- 0.048 | -0.015 +- 0.119 | -0.068 +- 0.141 |
| 8 | -0.015 +- 0.050 | +0.010 +- 0.110 | -0.049 +- 0.130 |

Every difference includes zero. The subharmonic start against the theory
reads +0.062 +- 0.036 (2 columns) and +0.031 +- 0.036 (8 columns) at 0.1 L0,
and +0.257 +- 0.083 and +0.189 +- 0.078 at 0.5 L0, so the excess is
significant at both lags for 2 columns and at 0.5 L0 for 8 columns.

The checkpoints (`kept_start_frames.csv`, the piston-removed variance and the
1 m Z-tilt at rows 0, 512, 2048, 5000 and 10000, 2 SE of 0.07 to 0.21 on the
variance and 0.24 to 0.46 on the tilt): the plain start reads 0.41 +- 0.07 on
the variance at row 0 and climbs to 0.87 to 1.10 by 5000 rows; the
subharmonic start reads 0.74 +- 0.10 at row 0 and stays between 0.71 and
0.89. The Z-tilt reads 0.7 to 1.2 everywhere, inside its own error bar, so
this run cannot separate the two starts on the tilt.

Three readings:

1. **The Q6 defect is confirmed from a second start, with an error bar.**
   No spin-up, either start, both `n_columns`: the extrusion-axis excess at
   0.5 L0 is +0.19 to +0.27 with a 2 SE of 0.08, and it matches the spun-up
   Q6 values inside 0.05. The 8-column excess is smaller than the 2-column
   one on both starts.
2. **The start is not kept, and it is not carried either.** The extrusion
   axis reads the same from a plain start and from a subharmonic start over a
   4 L0 record, so the recursion imposes its own axis statistics within the
   record whatever the first frame held. A parallel worker gains the first
   frame from a subharmonic start (Q7), and nothing on the extrusion axis.
3. **The plain-start variance recovers by about 0.2 L0 of rows, the
   subharmonic start needs none.** That repeats the Q7 reading with an error
   bar on the variance: the spin-up of the optical quantities is short, and
   it is removed by the subharmonic start.

## Products

| file | what it holds |
|---|---|
| `oversize_tilt.csv` | The Z-tilt and G-tilt variance of every arm, with the predicted capture. |
| `oversize_tilt.png` | The Z-tilt ratio against the oversize factor. |
| `oversize_dphi.csv` | The structure function of every arm over 12 separations. |
| `oversize_dphi.png` | The structure-function ratio against the separation. |
| `s27_settlement.csv` | The four generators against two estimators, at five separations. |
| `seed_quirk.csv` | The aotools shared-seed measurement over 80 samples. |
| `infinite_covariance.csv` | B(r) per axis, at two outer scales, against the theory. |
| `infinite_covariance.png` | The same covariance, plus a residual panel. |
| `infinite_dphi.csv` | The structure function of the two extruded screen classes. |
| `infinite_dphi.png` | The structure-function ratio of the two classes. |
| `infinite_tilt.csv` | The pooled Z-tilt and G-tilt variance of the extruded screens. |
| `stationarity_windows.csv` | The window table: D(r) and the Z-tilt against the cumulative row count. |
| `rowlag_covariance.csv` | The row-lag covariance C(k) over 257 lags, against the theory. |
| `stationarity.png` | The window-drift panel and the row-lag panel. |
| `n_columns_sweep.py` | Q6. The sweep script, Phases 1 and 2 of `PLAN_n_columns.md`. It runs in about 15 minutes. |
| `ncol_rowlag.csv` | Q6. rho(k) per axis, per cell, for every `n_columns`. |
| `ncol_rowlag.png` | Q6. rho(k) against the row lag, one panel per (L0, N). |
| `ncol_anisotropy.png` | Q6. D_ext(r) / D_trans(r) against the separation. |
| `ncol_piston.csv` | Q6. The chunk-mean variance and the periodogram of the row piston. |
| `ncol_piston.png` | Q6. The same two piston quantities. |
| `ncol_cost.csv` | Q6. The setup time, the `add_row` time, the size of `cov_zz`, and the guards. |
| `ncol_cost.png` | Q6. The same three cost quantities against `n_columns`. |
| `ncol_window.csv` | Q6, Phase 2. The apparent rho(1 L0) against the analysis window. |
| `ncol_window.png` | Q6, Phase 2. The same window sweep. |
| `n_columns_sweep.log` | Q6. The full stdout of the sweep, with every PASS and FAIL line. |
| `*_aotools107.csv` | The nine tables of the first run, on aotools 1.0.7. Compare a rerun table against its twin here. |
| `subharmonic_start.py` | Q7. The subharmonic-start check, three starts, 32 seeds. It runs in about 4.5 minutes. |
| `subharmonic_start.csv` | Q7. start, rows, metric, measured, theory, ratio. |
| `subharmonic_start.png` | Q7. The four ratios against the row count. |
| `kept_start.py` | Q8. The two starts with no spin-up on the production cell, 32 seeds, 10000 rows, 2 SE on every number. About 29 minutes. |
| `kept_start_rho.csv` | Q8. rho(k) on the extrusion axis per start and `n_columns`. |
| `kept_start_frames.csv` | Q8. The checkpoint ratios. |
| `kept_start.png` | Q8. rho(k) and the two checkpoint ratios. |

## Proposed row for validation/README.md

The owner adds this row. This study edits no tracked file.

```
| [screens/](screens/) | The low-frequency phase-screen study. It compares the Fourier, subharmonic, oversized-and-cropped, and extruded (infinite) screens against the analytic tilt and structure-function values. See [screens/FINDINGS.md](screens/FINDINGS.md). |
```
