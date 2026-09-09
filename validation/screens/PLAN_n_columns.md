# Plan: the `n_columns` exploration of the aotools infinite phase screen

Written 2026-09-08. This plan is a study of
`aotools.turbulence.infinitephasescreen.PhaseScreenVonKarman`. It answers one
question: can the extruded screen give an acceptable von Karman phase
covariance, and at what `n_columns`? It is VALIDATION ONLY. It changes no olb
module.

## 1. What changed since FINDINGS.md

`FINDINGS.md` (2026-08-28) ran on aotools 1.0.7. The installed build is now
`0.1.dev477+g8f38ff7c1.d20260605` (dist-info 1.0.8). Between the two, upstream
merged PR 111 (2026-04-02): `turb.phase_covariance` cast its input to FLOAT32,
and now casts to FLOAT64. Two open upstream issues describe the symptom of the
float32 kernel: issue 107 (noise grows on one side of the screen after many
`add_row` calls) and issue 109 (the variance, so r0, is not preserved). The
mechanism is the conditioning: at N = 512, dx = 1 cm, L0 = 25 m the stencil
covariance `cov_zz` has a condition number of 7e7, and float32 (eps 6e-8)
cannot invert it. So the A matrix of 1.0.7 was noise-limited on exactly the
study's parameters.

CONSEQUENCE: the FINDINGS Level 3 result ("the extrusion axis is GENUINELY
DEFECTIVE" at `n_columns = 2`) was measured with the broken kernel. It must be
rerun before `n_columns` is varied. A first probe on the fixed build (N = 128,
L0 = 2.56 m, 16384 rows, one seed) reads the row-lag correlation at one outer
scale as 0.006 against a theory of 0.004, where FINDINGS read a large excess.
That is one seed, and it is not yet a result.

Sources:
- https://github.com/AOtools/aotools/issues/107
- https://github.com/AOtools/aotools/issues/109
- https://github.com/AOtools/aotools/pull/111

## 2. Genuine errors found in the code (installed build)

| # | Where | Finding | Verdict |
|---|---|---|---|
| E1 | `turb.phase_covariance`, `r = numpy.float32(r)` (1.0.7) | Float32 covariance into a 7e7-conditioned inverse. Fixed upstream by PR 111. | REAL, FIXED in the installed build. Explains issues 107 and 109. |
| E2 | `turb.phase_covariance`, `r += 1e-40` | `numpy.float64(arr)` returns the SAME array (shares memory, verified 2026-09-08), so the function adds 1e-40 to the caller's array in place (`self.seperations`). | REAL, HARMLESS (1e-40 on separations of order 1e-2). Report upstream. |
| E3 | `PhaseScreen.makeAMatrix`, the `lstsq(rcond=1e-8)` fallback | A `print`, not a warning, and a truncated pseudo-inverse that changes the statistics silently. | NOT TRIGGERED at the study's parameters (Cholesky succeeds at N = 512, L0 = 25 m, `n_columns` 2 to 8; probed 2026-09-08). Guard it in the sweep: assert no fallback. |
| E4 | `PhaseScreen.makeBMatrix`, `svd` then `sqrt(W)` | A negative eigenvalue of BB^T becomes +abs. | NOT TRIGGERED: BB^T is positive (minimum eigenvalue 1.1e-2 in every probe). Log the minimum eigenvalue in the sweep. |
| E5 | `PhaseScreenKolmogorov.get_new_row`, the reference point | `A.(Z - ref) + ref` is exact only when the rows of A sum to 1. A is built from the raw covariance, not the difference covariance of Fried 2008. | HEURISTIC, not a bug for olb: the rows of A sum to 0.996 to 1.000 (probed). olb plans the VonKarman class only. |
| E6 | `make_initial_screen` | Plain `ft_phase_screen`, no subharmonics; the first frame is low in the outer-scale band. | KNOWN (FINDINGS Q4). The spin-up cures it. Not a bug. |
| E7 | `PhaseScreenVonKarman` docstring, "nCols" | The parameter is `n_columns`. | COSMETIC. |

No error in the geometry: X sits at row -1, the stencil at rows 0 to
`n_columns - 1`, `add_row` prepends, so the stencil of the next step is the
newest rows. Consistent.

## 3. Restore the figures and the tables

The PNGs of `validation/screens/` were NEVER committed (only the CSVs were,
and commit ea597e3 untracked them). So there is nothing to restore from git
for the figures; they come back from a rerun. The CSVs of the 1.0.7 run come
back from git:

```
git show ea597e3^:validation/screens/data/oversize_tilt.csv > validation/screens/data/oversize_tilt_aotools107.csv
```

(nine files, the names in FINDINGS "Products"). Keep them under `data/` with
the `_aotools107` suffix so the fixed-build rerun can be compared row by row.
The scripts write to `figures/` and `data/`; create both directories first.
`*.png` is tracked by `.gitignore`, `*.csv` is not: commit the figures, do not
force-add the tables.

## 4. The exploration, in phases

### Phase 0. Rerun FINDINGS on the fixed build (about 15 min)

Run the four commands of FINDINGS "How to reproduce" unchanged. Compare every
table against the restored 1.0.7 CSV. The expected outcome: Q1, Q2, S-27 and
the seed quirk do not move (they do not read the infinite screen); Q3 Level 3
and Q4 MOVE. Write the diff into FINDINGS as a new section "Rerun on the fixed
kernel". If Level 3 no longer shows the excess, the Q5 point 3 advice ("raise
`n_columns`") loses its evidence and point 4 (the shifted large screen) must be
re-weighed on cost only.

### Phase 1. The `n_columns` sweep (the new script `n_columns_sweep.py`)

The grid:
- `n_columns` in {2, 4, 8, 16, 32}.
- L0 in {2.56 m, 25 m} at r0 = 0.10 m, dx = 0.01 m (the FINDINGS parameters).
- N = 128 for the long records (cheap), N = 512 at `n_columns` in {2, 8} as
  the production-size check.
- 8 seeds for each cell. Spin-up 2 L0 / dx rows (512 at 2.56 m; 5000 at 25 m,
  FINDINGS Q5 point 2). Record 20 L0 / dx rows after the spin-up (5120 at
  2.56 m; 50000 at 25 m, only at N = 128).

The metrics, all against the analytic von Karman covariance
(`helpers.vk_covariance_closed`):
1. The extrusion-axis correlation rho(k) at lags 1 px to 4 L0. The band:
   |rho_meas - rho_theory| < 0.02 at every lag past 0.5 L0, pooled over the
   seeds (2 SE).
2. The transverse-axis covariance at the same lags (the control; it must hold
   as in FINDINGS).
3. The structure-function anisotropy D_ext(r) / D_trans(r) at r = 0.1 to 2 L0.
   The band: inside 5 percent.
4. The row-piston series m(t) = the mean of each new row: its variance against
   the transverse-averaged theory, and its power spectrum against the 1-D von
   Karman spectrum of a row average. This is the quantity the near-unit-root
   mode of the recursion controls. A single-seed probe (2026-09-08) showed
   20 m chunk means that wander between -0.45 and 1.9 rad while the theory
   covariance at 20 m is zero. The chunk-mean variance across seeds is the
   statistic.
5. The screen variance B(0) / theory, and the Z-tilt over D = 1 m
   (`helpers.zernike_tilt`), pooled.
6. The cost: setup wall time, the size of `cov_zz` ((`n_columns` N)^2
   doubles: 2 MB at 128 x 16, 512 MB at 512 x 16, 2 GB at 512 x 32), and the
   `add_row` time.
7. The guards of Section 2: no lstsq fallback (capture stdout), the minimum
   eigenvalue of BB^T, the condition number of `cov_zz`.

The output: `data/ncol_rowlag.csv`, `data/ncol_piston.csv`,
`data/ncol_cost.csv`; `figures/ncol_rowlag.png` (rho(k) per `n_columns`, one
panel per L0, the theory in black), `figures/ncol_anisotropy.png`,
`figures/ncol_piston.png` (the chunk-mean variance and the spectrum),
`figures/ncol_cost.png`.

### Phase 2. The record-length and piston question

The 2026-09-08 probes disagree between a 6144-row record (a flat rho floor of
0.13 to 0.16 at every `n_columns`) and a 16384-row record (rho at 1 L0 reads
the theory). So the excess FINDINGS measured on 512-row windows may be a
record-length effect and not a Markov-memory effect. Test it directly: at
`n_columns = 2` and 16, one long record of 100 L0 / dx rows, and compute
rho(k) in windows of 2, 5, 10, 20 and 100 L0. Plot the apparent rho(1 L0)
against the window length. If it converges to the theory with the window, the
"over-correlation" is the wander of the row piston across a short window, and
the cure is the spin-up plus a long record, not `n_columns`.

If the piston wander is real and slow (the chunk-mean variance exceeds the
theory), test one anchor: subtract the reference point the Kolmogorov class
uses (Fried 2008), or condition on the row mean. Measure it, do not adopt it.

### Phase 3. The decision for `olb/waveoptics/turbulence/temporal.py`

Restate FINDINGS Q5 with the fixed-build numbers. The two candidates:
- the extruded screen at the smallest `n_columns` that passes the Phase 1
  bands, with its spin-up rule in outer scales and its setup cost;
- the shifted large screen (K = 8 at L0 = 25 m: 4096 x 4096, 134 MB per
  layer).

The comparison is cost per frame and memory per layer at the production grid
(1024 px, dx from `GridSpec`), for the same pass bands. That is the owner's
decision. The plan does not build `temporal.py`.

### Phase 4. The speed evaluation of the time axis (added 2026-09-08, after the sweep)

THE QUESTION. Which screen route gives the cheapest CORRECT frozen-flow time
axis for a tracked LEO pass? The routes are (a) the extruded screen at the
`n_columns` Phase 1 certified, (b) the shifted large screen of FINDINGS Q5
point 4, and (c) the DOWN-SAMPLED extrusion: an extruded screen at a pitch set
by the layer's OWN r0, upsampled to the propagation grid, so that a sub-r0
move is a shift of the screen the layer already holds, and `add_row` runs only
when the accumulated move reaches one coarse pixel.

WHY (c) IS NOT THE BURIED P2 IDEA. P2 (`validation/waveoptics_speed/
coarse_screen_experiment.py`) coarsened EVERY screen by one uniform factor and
measured the SNAPSHOT scintillation. It died on its own kill line: the coarse
screen loses the Fresnel-scale structure that builds sigma2_I, and the FFT
zero-pad that hurts least erased the speed win. Route (c) is per LAYER, and
its purpose is the TIME axis, not the snapshot. The layers that move fastest
are the high, weak ones with the largest r0, so the pitch that each layer
needs is the one to test. The kill line is the same (sigma2_I inside 5
percent, the structure function at the fine pitch) plus the temporal one.

THE PRODUCTION NUMBERS (probed 2026-09-08 with `turbulent_grid`, preset
`standard`, 500 km orbit, 0.7 m ground aperture, 1550 nm):

| elevation | grid | pixel | side | screens | r0 total | slew |
|---|---|---|---|---|---|---|
| 30 deg | 512 px | 6.9 mm | 3.5 m | 9 | 12.7 cm | 0.44 deg/s |
| 90 deg | 512 px | 5.4 mm | 2.8 m | 9 | 19.3 cm | 0.87 deg/s |

Per screen at 30 deg (height, the layer r0, the apparent translation
`omega_slew * z_slant`, and the rows per second at the FULL pitch against the
rows at a pitch of r0_layer / 4):

| height | r0 layer | move | rows/s at 6.9 mm | rows/s at r0/4 |
|---|---|---|---|---|
| 16.2 km | 2.25 m | 247 m/s | 36,000 | 440 |
| 9.7 km | 1.74 m | 148 m/s | 21,000 | 340 |
| 4.1 km | 1.13 m | 63 m/s | 9,200 | 220 |
| 2.0 km | 0.79 m | 30 m/s | 4,400 | 150 |
| 0.6 km | 0.42 m | 8.7 m/s | 1,300 | 80 |
| 0.08 km | 0.15 m | 1.2 m/s | 180 | 30 |

The Bufton wind (`olb.turbulence.profiles.v_wind`, 10 to 30 m/s) adds to the
move as a vector, so the numbers above are the floor of the top layers and
the bulk of the ground layer. The side of the production grid is 0.14 L0 at
L0 = 25 m, which is the regime where the Phase 1 sweep shows the extrusion
over-correlates along its axis at every `n_columns` tested. So route (a) has
to carry that excess (about 0.2 in rho at 0.5 L0 at N = 512, `n_columns`
8), and Phase 4 must say what it costs on the TEMPORAL spectrum.

THE EXPERIMENT (`validation/screens/temporal_speed.py`, validation only):
1. One layer at a time, at the 30 deg plan: for each of the 9 screens, the
   three routes give a 2 s record of frames at 1 kHz (the frame rate is a
   parameter). Route (b) uses one K = 8 Fourier screen at L0 = 25 m and a
   Fourier-interpolated sub-pixel shift. Route (c) sweeps the coarse pitch
   dx_c / r0_layer in {1/2, 1/4, 1/8, 1/16} and the upsample by FFT zero-pad
   ONLY (P2 showed bicubic damages sigma2_I most); the shift between rows is
   the Fourier shift of the coarse screen, which is cheap because it runs on
   the coarse grid.
2. The pass bands for each layer route: the phase structure function of one
   upsampled frame against 6.88 (r / r0)^(5/3) at r from the fine pitch to
   the Fresnel scale sqrt(lambda z) of that layer, inside 5 percent; the
   temporal phase structure function D(tau) of one pixel against the
   frozen-flow theory D(v tau) (Taylor, DOI 10.1098/rspa.1938.0032), inside 5
   percent up to tau = r0_layer / v; and the Greenwood-like temporal spectrum
   slope -8/3 (Andrews and Phillips, DOI 10.1117/3.626196, Ch. 12).
3. The end-to-end check: the full 9-screen split step over the 2 s record on
   each route, against the snapshot statistics of the campaign (the same
   grid, the same seed rule): the aperture sigma2_I and the SMF eta mean
   inside 5 percent (the P2 kill line), and the fade DURATION at the 5
   percent level within 10 percent between routes (b) and (c).
4. The cost: wall time per frame per layer, and per frame end to end, plus
   the memory per layer (route (b) 134 MB at 4096 px double; route (c)
   the coarse screen plus its (n_columns N_c)^2 covariance, which at
   dx_c = r0/4 and the 3.5 m side is N_c of order 2 to 24 px).

THE KILL LINES, stated before the run: route (c) DIES if no pitch coarser
than the fine pitch passes bands 2 and 3 on the top five layers; route (a)
DIES if its temporal spectrum on the extrusion axis reads more than 10
percent slow at the fade-relevant frequencies (the Q3 Level 3 defect is
exactly a slower temporal axis under Taylor); route (b) is the reference
and cannot die, it can only cost too much.

## 5. Products and the record

- `validation/screens/n_columns_sweep.py` (new; it reuses `helpers.py` and
  the `row_lag_sums` of `extrusion_stationarity.py`, no new physics).
- The rerun figures of the three existing scripts, under `figures/`.
- FINDINGS.md: a new "Rerun on the fixed kernel" section, a new "Q6 -
  `n_columns`" section, and the aotools version line corrected.
- An upstream note for E2 (the in-place mutation) on the aotools tracker,
  owner's call.
- `docs/schmidt-crosscheck.md` rows S-27 and the forward map: the FINDINGS
  side finding still asks the owner to correct them.
- Phase 4 (not started): `validation/screens/temporal_speed.py`, its tables
  under `data/temporal_*.csv`, its figures under `figures/temporal_*.png`,
  and a "Q7 - the time axis" section in FINDINGS. Its verdict feeds the
  design of `olb/waveoptics/turbulence/temporal.py`, which stays a stub
  until the owner decides.
