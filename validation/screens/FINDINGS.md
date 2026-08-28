# Phase-screen low-frequency study: findings

## What this study is

A phase screen on a finite grid holds no power below its grid fundamental.
That missing low-frequency band is the tip and the tilt. This study measures
how much tilt each screen route holds, and it compares each route against the
analytic value. The study answers five questions. It also settles two open
rows in `docs/schmidt-crosscheck.md`, and it gives a design recommendation for
`olb/waveoptics/turbulence/temporal.py`.

The study is VALIDATION ONLY. It reads the production layer. It changes no olb
module.

Four scripts live in `validation/screens/`:

- `helpers.py` — the analytic truths and the shared estimators.
- `oversize_crop.py` — arm 1, the Fourier and the oversize-and-crop screens.
- `infinite_screen_stats.py` — arm 2, the spatial statistics of the extruded
  screens.
- `extrusion_stationarity.py` — arm 3, the drift test of the extrusion.

Common parameters: lambda = 1550 nm, r0 = 0.10 m, dx = 0.01 m, reference grid
512 pixels (side 5.12 m), pupil D = 1.0 m (D/r0 = 10), von Karman L0 = 25 m
where the spectrum is finite. The screens come from aotools 1.0.7.

The analytic truths are the Noll per-axis Z-tilt filter integral (Noll,
DOI 10.1364/JOSA.66.000207), the Andrews G-tilt filter (Andrews and Phillips,
DOI 10.1117/3.626196, Ch. 6, Eq. (84), printed p. 201), the von Karman
covariance in a float64 closed form (Assemat and Wilson,
DOI 10.1364/OE.14.000988, Eq. (5)) and as a Hankel transform of the phase power
spectral density (Schmidt, DOI 10.1117/3.866274, Ch. 9), and the structure
function D(r) = 6.88 (r/r0)^(5/3) (Schmidt, DOI 10.1117/3.866274, Ch. 9,
printed p. 160).

## How to reproduce

Run each command from the repository root. The four commands take about 15
minutes in total.

```
python -m validation.screens.helpers
python -m validation.screens.oversize_crop
python -m validation.screens.infinite_screen_stats
python -m validation.screens.extrusion_stationarity
```

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

## Proposed row for validation/README.md

The owner adds this row. This study edits no tracked file.

```
| [screens/](screens/) | The low-frequency phase-screen study. It compares the Fourier, subharmonic, oversized-and-cropped, and extruded (infinite) screens against the analytic tilt and structure-function values. See [screens/FINDINGS.md](screens/FINDINGS.md). |
```
