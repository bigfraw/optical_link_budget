# The fidelity-2 SMF fade-tail convergence study

Backlog items **2-I2T** (the tail-convergence question) and **2-N6** (the
large-campaign measurements).

## The question

Does the deep single-mode-fibre (SMF) fade tail of a fidelity-2 downlink
CONVERGE as the near-ground `Cn2` is resolved with more, thinner phase screens?

The MEAN of a fidelity-2 space budget is already validated flat against the
screen count (`docs/schmidt-crosscheck.md`, work package 7). The TAIL is the
open risk, because the tail sets the link availability margin. The post-WP7
matched-seed re-test saw about 2 dB of p5 movement between two screen
placements, but at 200 trials that movement stayed inside the Monte-Carlo
noise. This study resolves the movement above the noise.

## The method

The study PINS the grid and it moves the screens only. The shipped sizer widens
and refines the grid when the screen count grows, so a naive sweep of
`min_screens` moves two variables at one time. Each case is one
`olb.waveoptics.turbulence.Campaign`: one scenario, one geometry, one seed, one
grid and one screen plan. A campaign is resumable and idempotent, so a killed
run continues where it stopped.

Each trial gives one composite loss

    loss_db = -10 log10(collected_power * smf_eta)

which is the fidelity-2 downlink SMF quantity of `olb.models.waveoptics`. The
vacuum references are one constant for one grid, so they cancel in a comparison
between the pinned cases.

`pX` is the loss that the link EXCEEDS X percent of the time. So `p5` is the
loss exceeded 5 percent of the time, and `p1` is the 1 percent worst loss. The
table also gives the FADE DEPTH `pX - p50`, so two cases with a different vacuum
coupling still compare.

**No correction.** The fidelity-2 layer applies no tip-tilt removal and no
adaptive optics (backlog 2-AO), and it has no correction stage that a switch
could turn on. So the SMF coupled power here is the RAW uncorrected atmosphere,
and the verdict applies to that case only.

**The second primary quantity is the POINT irradiance.** Each stored trial
gives the irradiance `I` of the centre pixel of the receive field, and the
study reports the point fade `-10 log10(I / <I>)` with the same quantiles and
the same bootstrap. Nothing averages a point, so it is the sharpest probe of
the near-ground resolution, and a fibre pays the point figure (backlog 2-I3).

**The aperture index `sigma2_P` is a footnote only.** A 700 mm bucket averages
the irradiance strongly, so that index is of the order 1e-3 and it is
insensitive to the screen plan. A flat `sigma2_P` proves nothing about the
tail.

## The scenario

A 1550 nm downlink. `Site(cn2_ground=1.7e-14, wind_rms_m_s=21.0)`,
`Channel(altitude_m=500e3)`. The ground terminal has a 700 mm aperture, 2 urad
of pointing jitter and an SMF detector. The space terminal has a 100 mm
aperture, 1 urad of jitter and a 40 mm waist transmitter at 30 dBm. The geometry
is a 500 km circular orbit at one elevation (30 deg by default; 20 deg is the
low-elevation run).

## The cases

| case | what it isolates |
| --- | --- |
| `rapid` | the shipped RAPID preset: its own grid AND its 5-screen plan. It is what a `preset="rapid"` user gets, end to end. |
| `prod` | the shipped default. The standard sizer picks the grid AND the plan. It is the production baseline of the budgets. |
| `pin05` | 5 screens (the rapid floor) on the pinned grid. Against `rapid` it isolates the rapid GRID; against `pin09` the rapid COUNT. |
| `pin07` | 7 screens (the WP7 converged count) on the pinned grid. |
| `pin09` | 9 screens on the pinned grid. Against `prod` it isolates the GRID. |
| `pin15` | 15 screens on the pinned grid. It isolates the SCREEN COUNT. |
| `pin25` | 25 screens on the pinned grid. It isolates the SCREEN COUNT. |
| `pin40` | 40 screens on the pinned grid. It isolates the SCREEN COUNT. |
| `gnd09x4` | the `pin09` plan with the BOTTOM screen split into four sub-screens of equal integrated `Cn2`. It isolates the NEAR-GROUND resolution alone. |

The pinned grid is the grid that the sizer gives for `min_screens=15`. At 30 deg
that is 1024 px on a 3.514 m side (3.43 mm pixel); the shipped default gives
512 px (6.86 mm pixel).

The `gnd09x4` split inverts the cumulative slant `Cn2` integral to recover the
top height of the bottom slab, cuts that slab at equal cumulative `Cn2`, and
puts each sub-screen at the `Cn2`-weighted centroid of its sub-slab. The script
ASSERTS that the split conserves the parent integrated `Cn2`, the parent Rytov
variance and `r0_total`.

## How to run it

Run it from the repository root:

    python -m validation.tail_convergence.tail_convergence \
        --elevation 30 --n-trials 1000 --workers 16 --block-size 50

Then repeat it at the low elevation:

    python -m validation.tail_convergence.tail_convergence --elevation 20

Options:

- `--cases prod pin09 ...` selects the cases. The default is all nine.
- `OLB_TAIL_ROOT=<dir>` (an environment variable) moves the campaign root, for
  a smoke test that must not touch the study store.
- `--analyse-only` skips every run and reads what is stored. A partial campaign
  is fine: the analysis reads `n_stored`.

## How to resume

A `Campaign` computes only the MISSING blocks, so a killed run needs no special
handling: run the same command again. The campaign directories are

    validation/tail_convergence/campaigns/el{elevation}/{case}/

and they are gitignored. A campaign directory holds ONE physics case, so do not
change `--block-size` or the seed between two runs of the same directory; the
manifest check raises if you do.

## The outputs

- `tail_convergence_el{elevation}_results.json` — every measured number.
- `tail_convergence_el{elevation}.log` — the run log, with the progress lines
  and the wall time of each run.
- `figures/tail_vs_screens.png` — the quantiles against the screen count.
- `figures/survival.png` — the empirical survival function of every case.
- `figures/growth_with_trials.png` — how the tail estimate settles with the
  trial count (backlog 2-N6).

## The results (2026-09-04, 30 deg, 1000 trials for each case)

`pX` is the SMF loss exceeded X percent of the time, in dB. `pt pX` is the
POINT irradiance fade `-10 log10(I/<I>)` exceeded X percent of the time. The
`+-` values are 68 percent bootstrap half-widths (1000 resamples). `pin40` and
the 20 deg run were NOT run (owner decision, 2026-09-04): the 30 deg series
below answers the question, and the rapid end of the count axis was judged
more useful than a fifth pinned count.

| case | screens | grid [px] | pixel [mm] | bottom h [m] | bottom Cn2 share | mean | p50 | p10 | p5 | p1 | depth p5 | depth p1 | pt p5 | pt p1 | sigma2_I point | sigma2_P | s/trial | MB |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `rapid` | 5 | 256 | 10.31 | 140 | 0.853 | 18.15 | 17.14 | 26.66 +-0.64 | 29.90 +-0.38 | 38.11 +-1.70 | 12.76 | 20.97 | 3.83 +-0.16 | 5.75 +-0.37 | 0.214 | 0.005 | 0.15 | 28.9 |
| `prod` | 9 | 512 | 6.86 | 80 | 0.730 | 18.71 | 17.96 | 27.47 +-0.64 | 31.47 +-0.65 | 37.64 +-1.55 | 13.51 | 19.68 | 3.87 +-0.18 | 6.11 +-0.75 | 0.239 | 0.005 | 2.02 | 65.5 |
| `pin05` | 5 | 1024 | 3.43 | 140 | 0.853 | 18.74 | 17.90 | 28.00 +-0.65 | 31.51 +-0.70 | 38.74 +-0.95 | 13.61 | 20.83 | 4.28 +-0.16 | 5.86 +-0.34 | 0.239 | 0.005 | 3.90 | 261.7 |
| `pin07` | 7 | 1024 | 3.43 | 99 | 0.793 | 18.57 | 17.69 | 27.07 +-0.33 | 30.95 +-0.72 | 38.03 +-1.29 | 13.26 | 20.34 | 4.35 +-0.18 | 5.92 +-0.25 | 0.243 | 0.005 | 5.43 | 261.7 |
| `pin09` | 9 | 1024 | 3.43 | 80 | 0.730 | 18.59 | 17.59 | 27.49 +-0.43 | 31.20 +-1.08 | 37.71 +-1.46 | 13.61 | 20.11 | 4.29 +-0.15 | 6.17 +-0.23 | 0.259 | 0.005 | 9.63 | 261.7 |
| `pin15` | 15 | 1024 | 3.43 | 55 | 0.593 | 18.51 | 17.71 | 27.06 +-0.47 | 30.45 +-0.59 | 38.02 +-2.12 | 12.73 | 20.31 | 3.97 +-0.11 | 5.81 +-0.25 | 0.234 | 0.005 | 9.29 | 261.7 |
| `pin25` | 25 | 1024 | 3.43 | 39 | 0.470 | 18.22 | 17.37 | 26.73 +-0.31 | 29.04 +-0.68 | 37.60 +-1.35 | 11.67 | 20.23 | 4.18 +-0.25 | 5.88 +-0.26 | 0.232 | 0.005 | 14.37 | 261.7 |
| `gnd09x4` | 12 | 1024 | 3.43 | 13 | 0.182 | 18.19 | 17.45 | 26.58 +-0.32 | 29.79 +-0.62 | 36.75 +-1.11 | 12.34 | 19.30 | 4.19 +-0.21 | 5.79 +-0.21 | 0.255 | 0.005 | 12.65 | 261.7 |

`s/trial` is the mean stored wall time of one trial inside its worker (the
wall time of the whole campaign is in the 2-N6 section). `sigma2_P` is the
700 mm bucket index; it is 0.005 in every case and it proves nothing.

### The verdict

1. **The grid is a null.** `prod` (512 px) against `pin09` (1024 px), the same
   9-screen plan: every SMF quantile agrees inside its bar (p5 +0.27 dB,
   0.2 sigma; mean 18.71 against 18.59). The shipped 6.86 mm pixel does not
   under-resolve the SMF coupling. The point columns of the two DO differ
   (index 0.239 against 0.259), and a re-bin check on the stored fields showed
   why: a "point" is one pixel, and at a MATCHED averaging area the two grids
   agree exactly (6.86 mm: 0.239 against 0.246; 13.7 mm: 0.225 against 0.226;
   27.4 mm: 0.194 against 0.195). The index falls about 5 percent per pixel
   doubling, because the bottom screen puts real irradiance structure at its
   ~1.6 cm Fresnel scale. So the point comparison is clean ACROSS the pinned
   series only (one 3.43 mm pixel), and the point column of `prod` and `rapid`
   is not comparable to it.
2. **The SMF mid-tail moves with the screen count, in the SAFE direction, and
   it is not converged at 25 screens.** On the pinned grid, p5 falls
   31.51 -> 30.95 -> 31.20 -> 30.45 -> 29.04 dB from 5 to 25 screens
   (-2.47 dB, 2.5 sigma), p10 falls 1.3 dB, the mean and p50 fall 0.5 dB.
   MORE screens give LESS fade. The near-ground split `gnd09x4` (the 9-screen
   plan with its 80 m ground screen cut into four at 174 / 89 / 44 / 13 m)
   lands on the 25-screen line (p5 29.79, -1.41 dB against `pin09`, 1.1 sigma),
   so the live variable is the resolution of the GROUND layer, as backlog 2-I2
   supposed. The trend past 25 screens is UNRESOLVED (`pin40` not run).
3. **The deep tail does not move.** p1 reads 38.74 / 38.03 / 37.71 / 38.02 /
   37.60 dB across the pinned series (spread 1.14 dB, 0.7 sigma against the
   +-1.5 dB bars) and 36.75 for the ground split. So the 99 percent
   availability margin of the default plan is unchanged by the count; the
   95 percent margin is about 2 dB PESSIMISTIC at 9 screens against 25.
4. **The point irradiance does not move at all.** Across the pinned series the
   point p5 spreads 0.37 dB and the point p1 0.36 dB (0.1 to 0.3 sigma), and
   the point index sits at 0.23 to 0.26 at every count (the analytic plane-wave
   Rytov variance of this slab is 0.22). So the count effect is NOT
   scintillation. It is the PHASE statistics that the fibre overlap pays (the
   tip, the tilt and the low orders of the ground layer), which a pixel does
   not see and a bucket receiver would not see either.
5. **The screen generator is not the cause.** The separate phase-only study
   `validation/screen_stacking/` stacks the same per-screen r0 lists with no
   propagation: every configuration holds 0.75 to 0.84 of the Noll
   piston-removed aperture variance and 1.00 of the tilt-removed variance, so
   the generator misses the SAME tip-tilt fraction at every count (the 2-N2
   deficit), the ground layer as 1 or 4 screens agrees inside 0.3 sigma, and
   the whole plan at 5 -> 25 screens changes by -0.024 +-0.034 (0.7 sigma).
   That 20 percent "deficit" is the OUTER SCALE: the grid holds scales up to
   27 x its side (95 m) and the L0 = inf reference has no limit; the same
   screens match a von Karman L0 = 95 m theory exactly, and asked for
   L0 = 25 m they read 1.00 +-0.03 for one screen and 0.94 to 0.97 for a
   5- to 25-screen plan (a mild, few-percent stacking drift). The count trend
   of item 2 therefore comes from the PROPAGATION and placement side of the
   product, not from the screens; the outer scale itself is backlog 2-P5
   (HIGH), worth an estimated (not measured) 2 dB at p5 of this SMF tail.
6. **The rapid rung holds.** `rapid` as shipped (256 px, 10.31 mm, 5 screens)
   reads p5 29.90 and p1 38.11, inside the spread of the standard series, at
   1/30 of the cost (24 s for 1000 trials). Against the same 5-screen plan on
   the fine grid (`pin05`) its p5 is 1.6 dB lower (2 sigma) and its mean
   0.6 dB lower, so the 10.3 mm pixel starts to soften the SMF coupling where
   the 6.86 mm pixel does not; its p1 agrees.
7. **The post-WP7 hint is reversed.** That 200-trial re-test read MORE fade
   with thin near-pupil screens. This matched-grid, 1000-trial series reads
   LESS fade with a finer ground layer, at every quantile. Trust this one.

WHAT THIS MEANS FOR THE DEFAULT. The 9-screen `standard` plan is converged for
the mean, the deep tail and every point quantity; it is about 2 dB pessimistic
at p5 and 1 dB at p10 against a 25-screen plan, and the true converged p5 is
not known (it may sit below 29 dB). No re-tiering is forced: the error is on
the safe side and the 99 percent margin is unaffected. The number belongs in
the 2-I3 preset revision as the SMF (fibre) row of the catalogue.

SCOPE. One elevation (30 deg), one site profile (HV5/7 at the hero site), one
uncorrected 0.7 m SMF receiver, no tip-tilt and no AO (backlog 2-AO). The
statement is "the default plan is converged for THIS link to the numbers
above", not "fidelity 2 is right".

### The 2-N6 numbers

- **Wall time.** At 512 px (`prod`): 1000 trials in 172 s on 16 workers
  (0.17 s/trial). At 256 px (`rapid`): 24 s. At 1024 px with 9 screens: 781 s
  on 16 workers (0.78 s/trial); with 5 / 7 / 15 / 25 screens on 8 workers:
  577 / 794 / 1351 / 1333 s (0.58 / 0.79 / 1.35 / 1.33 s/trial; `pin25`
  resumed from 8 stored blocks, so its wall time covers 12 blocks). The cost
  is close to linear in the screen count, because the screen generation
  dominates a 1024 px trial.
- **Disk.** 261.7 MB for 1000 trials at 1024 px (the field patch of a 0.7 m
  aperture at 3.43 mm, complex64), 65.5 MB at 512 px, 28.9 MB at 256 px. The
  eight campaigns hold 1.7 GB.
- **Memory, and a kill.** A 16-worker pool at 1024 px with 15 screens ran out
  of memory (about 1 GB for each worker against 9 GB free on the machine), and
  a worker `MemoryError` ended the `run` call. Every finished block stayed on
  disk. A second call with 8 workers resumed from those blocks with no rerun,
  and the three complete cases returned in 0.0 s with their checkpoint lines.
  So the resume after a kill WORKS, and `Campaign.run` has no memory guard: the
  caller sizes the pool. The load of one 1000-trial case with its fields holds
  262 MB of fields; the scalars alone are negligible.
- **How the tail settles with the trial count** (`pin09`):

  | n | p50 | p10 | p5 | p1 |
  | --- | --- | --- | --- | --- |
  | 100 | 17.42 | 27.96 +-1.93 | 29.94 +-1.01 | 34.93 +-3.27 |
  | 200 | 17.42 | 26.75 +-1.47 | 29.94 +-1.47 | 37.63 +-3.17 |
  | 400 | 17.89 | 26.88 +-0.74 | 29.94 +-1.47 | 37.63 +-1.65 |
  | 600 | 17.87 | 27.18 +-0.58 | 31.00 +-1.26 | 40.57 +-1.10 |
  | 800 | 17.63 | 27.24 +-0.41 | 30.70 +-1.20 | 40.04 +-1.48 |
  | 1000 | 17.59 | 27.49 +-0.43 | 31.20 +-1.08 | 37.71 +-1.46 |

  The p1 estimate wanders over 5 dB between 100 and 800 trials and settles to
  +-1.5 dB at 1000; p5 settles to +-1 dB. So 1000 trials resolve a p5 shift of
  about 1.5 dB and a p1 shift of about 2 dB, and the ten-samples-past-p1 rule
  of `EmpiricalSampler` is the right floor for a 99 percent statement.

## Notes

This study is VALIDATION ONLY. It reads the production layer and it changes no
`olb` module. It uses the private planner helpers of
`olb.waveoptics.turbulence.sampling` to build the near-ground refinement.
