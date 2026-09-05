# The fidelity-2 speed plan

Status: DONE, 2026-08-29. P0 to P4 all ran; see Section 8 for the findings and
the execution design. The evidence lives in `validation/waveoptics_speed/`.
Sections 1 to 7 below are the original work plan, kept as written.

This document is a work plan. It holds two work packages and five
ready-to-paste subagent prompts (Section 7). Each prompt is self-contained:
paste it into one Opus subagent with no other context. Run P0 first. P1, P2,
and P3 are independent after P0. Run P4 last, because it reads the results of
the others.

## 1. Goal

Fidelity-2 trials are slow, and the verification campaigns ahead run
thousands of them. Two tracks attack the cost:

- **Work package A — per-trial efficiency.** Make one trial cheaper: faster
  screen generation, and a grid that follows the beam.
- **Work package B — execution planning.** Run the trials better: measure
  what is truly parallel, compare threads with processes, and decide what to
  cache and how a campaign executes.

The rule for both tracks: **measure first**. No production change lands
before the profile says where the time goes, and no change lands without a
statistical acceptance test.

## 2. Ground rules (these go into every task)

- All documentation, comments, and commit messages use ASD-STE100 Simplified
  Technical English. See `CONVENTIONS.md`.
- Every equation cites its source by DOI.
- Ponytail: the laziest solution that works. Borrow the shared kernels. No
  speculative abstraction.
- **The schmidt layer is now a resource (owner decision, 2026-08-29).** The
  old rule said `olb/waveoptics/schmidt/` is validation-only. For this speed
  work the rule is RELAXED: read that layer, and borrow or wire its kernels
  where they help. The LightPipes production core still keeps its bodies.
- aotools stays import-only (LGPL-3.0). Do not copy its code. A new
  generator written from the Schmidt book equations is allowed in
  production.
- The seed contract is a design constraint: for one seed, trial k is
  bit-identical, and the trial count does not change it (see
  `olb/waveoptics/turbulence/run.py`, `_screen_seed`). A new generator
  gives DIFFERENT draws. That is accepted, but it must be a VERSIONED,
  opt-in choice, never a silent replacement.
- `TurbTrial` stays a frozen scalar record (owner decision). Do not extend
  it.
- Validation scripts follow the `validation/` pattern: one script, one
  results JSON, one run log.

## 3. The baseline facts

Facts read from the code on 2026-08-29:

- One trial = one screen stack (`phase_screen` per screen) + one
  `split_step` + the scalar reads (`_clip`, `Power`, `coupling_efficiency`,
  `mmf_coupling_efficiency`, the uplink overlap). See
  `olb/waveoptics/turbulence/run.py`, `run_one`.
- `phase_screen` wraps aotools `ft_sh_phase_screen`. Per screen that is one
  full FFT draw PLUS 3 subharmonic levels x 9 modes, and each mode builds a
  full N x N complex exponential in a Python loop. Three structural wins are
  visible in the aotools source:
  1. The sqrt-PSD filter depends on (N, dx, L0, l0) only. r0 enters as the
     scalar factor r0^(-5/6). So the filter caches per grid.
  2. Each subharmonic mode exp(2 pi i (fx x + fy y)) is a separable outer
     product of two N-vectors. The 27-mode sum collapses to small matrix
     products.
  3. One complex Gaussian draw + one FFT gives TWO independent real screens
     (the real part and the imaginary part).
- `Forvard` uses single-threaded `numpy.fft` (pocketfft). It releases the
  GIL. Each hop is 2 full-grid FFTs. A space standard run has 10+ hops.
- The turbulent grid is sized ONCE for the whole path
  (`olb/waveoptics/turbulence/sampling.py`). The space slab input is a
  plane wave that FILLS the grid.
- `Threader` (threads) exists and gives a real speed-up, but no scaling
  study exists. Processes were never measured. The platform is Windows
  (spawn, not fork).
- `olb/waveoptics/schmidt/` already holds: a book-cited screen generator
  that takes an `rng` (`turbulence.ft_phase_screen`, `subharmonic_screen`,
  Listings 9.2 and 9.3); the per-plane-pitch propagation chain
  (`fresnel.partial_propagations`, Ch. 8, Eq. (8.18)); and the sampling
  constraints 1 to 4 (`sampling.constraint1_max_delta2` ...
  `constraint4_min_n`, `check_sampling`, `partial_grid_spacing`,
  `partial_max_step`, `partial_plane_count`).
- `olb/waveoptics/turbulence/screens.py` documents WHY a naive
  coarse-then-interpolate screen is wrong: the coarse screen is band-limited
  above its own Nyquist frequency, so it loses the Fresnel-scale structure
  that builds the scintillation (Schmidt, DOI 10.1117/3.866274, Sec. 9.4,
  printed p. 172). The interpolation experiment (P2) must measure against
  that documented failure mode, not assume it away.

## 4. Work package A — per-trial efficiency

- **A0 (prompt P0). The profiling baseline.** Split one trial's wall time
  into its parts, across representative cases. This gates everything: the
  later tasks attack the biggest slice first.
- **A1 (prompt P1). The fast screen generator.** An olb generator built on
  the schmidt book kernels, with the three structural wins of Section 3.
  Opt-in, versioned, statistically validated.
- **A2 (prompt P2). The two grid ideas, measured honestly.**
  (a) Coarse screens + interpolation, with a stated kill criterion.
  (b) A grid that follows the beam: beam-sized screens on the flat grid,
  and the per-plane-pitch chain of gap S-14, built on the schmidt
  machinery. Report only; no production wiring.

## 5. Work package B — execution and caching

- **B1 (prompt P3). The parallel scaling study.** The Threader scaling
  curve, a process-pool comparison, and a batched-FFT prototype, with a
  memory budget. Output: a recommendation per case.
- **B2 (prompt P4). Caching and the campaign design.** Decide what to
  cache at which level, then write the execution design for a fidelity-2
  campaign, and build the smallest useful cache.

## 6. Order and gates

1. **P0 first.** Its JSON gates the other tasks.
2. **P1, P2, P3 in any order after P0.** They are independent.
3. **P4 last.** It reads the results of P0, P1, and P3.

Acceptance is STATISTICAL everywhere: a change to the generator or the grid
changes the random draws, so bit-exact comparison is impossible. The tests
compare converged statistics (mean collected power, the aperture
scintillation index sigma2_I, smf_eta) inside the Monte-Carlo error bars of
the trial count.

## 7. The subagent prompts

Each block below is one complete prompt. Paste it into one Opus subagent.

### P0 — the profiling baseline (run this first)

```text
You work in the repository D:\repos\optical_link_budget (package `olb`).
Model guidance: be the laziest engineer that still works correctly. Borrow
the shared kernels. No speculative abstraction. All documentation, comments,
and commit messages use ASD-STE100 Simplified Technical English (see
CONVENTIONS.md). Every equation cites its source by DOI. Run modules with
`python -m ...` from the repository root.

TASK. Build a profiling baseline for the fidelity-2 turbulent wave-optics
trials, so later speed work attacks the biggest cost first.

READ FIRST:
- olb/waveoptics/turbulence/run.py   (run_one: the trial body)
- olb/waveoptics/turbulence/screens.py (phase_screen wraps aotools)
- olb/waveoptics/turbulence/splitstep.py (split_step, hop structure)
- olb/waveoptics/turbulence/sampling.py (turbulent_grid, PRESETS)
- olb/waveoptics/propagators.py (Forvard: 2 numpy FFTs per hop)
- validation/ (the script + results-JSON + run-log pattern)

BUILD: validation/waveoptics_speed/profile_baseline.py

The script measures, for ONE trial of each case, the wall time of:
1. Screen generation, split into (a) the base FFT screen
   (aotools ft_phase_screen) and (b) the subharmonic addition (the
   difference against ft_sh_phase_screen). Time each screen of the stack.
2. The split step, split into hops. Derive the hop and sub-step count from
   the plan (z_m, z_total_m, forvard_max_z) and report seconds per FFT.
3. The scalar reads: the receive clip, Power, coupling_efficiency (put an
   SMF detector on the receive terminal so this path runs), and the uplink
   overlap for the uplink case.
Reproduce the run_one body inside the script with timers around each part.
Import the private helpers the way the run.py self-check does. Use the same
seeds as propagate_turbulent_scenario would (seed=0, trial 0), so the
numbers are repeatable.

CASES (build them the way the sampling.py and run.py self-checks do):
- Terrestrial: 2 km, Cn2 = 5e-15, standard preset, SMF receiver.
- Space downlink at 30 deg elevation, 600 km orbit: rapid, standard, AND
  reference presets.
- Space uplink at 30 deg, standard preset (adds the overlap read).
Record for each case: grid n, side, pixel, screen count, sub-step count,
and the achieved report numbers.

ALSO measure raw FFT cost: time numpy.fft.fft2 on complex128 arrays of
n = 512, 1024, 2048, 4096, and the same with scipy.fft.fft2(workers=1) for
reference. Record numpy and scipy versions, CPU name, and core count.

OUTPUT:
- validation/waveoptics_speed/profile_baseline_results.json (all numbers).
- A printed table: per case, the share of screen generation (FFT part /
  subharmonic part), propagation, and scalar reads.
- validation/waveoptics_speed/profile_baseline_run.log (the run output).
- A short RESULTS section at the top of the script docstring: the biggest
  slice per case, in one sentence each.

DO NOT change any production code. The deliverable is the script, the JSON,
the log, and the docstring summary.
```

### P1 — the fast screen generator

```text
You work in the repository D:\repos\optical_link_budget (package `olb`).
Model guidance: be the laziest engineer that still works correctly. Borrow
the shared kernels. No speculative abstraction. All documentation, comments,
and commit messages use ASD-STE100 Simplified Technical English (see
CONVENTIONS.md). Every equation cites its source by DOI. Run modules with
`python -m ...` from the repository root.

TASK. Add a FAST phase-screen generator to
olb/waveoptics/turbulence/screens.py, beside the existing aotools-backed
phase_screen. Opt-in only: the default path stays aotools, so old runs stay
reproducible. The new generator gives DIFFERENT random draws; that is
accepted and documented.

READ FIRST:
- olb/waveoptics/turbulence/screens.py (the current wrapper, its docstring
  rules, and its self-check: the structure-function tests are the template
  for yours)
- olb/waveoptics/schmidt/turbulence.py (ft_phase_screen,
  subharmonic_screen: a book-cited generator that takes an rng. OWNER
  DECISION 2026-08-29: the schmidt layer may be borrowed and wired for this
  speed work; the old validation-only rule is relaxed here.)
- olb/waveoptics/turbulence/run.py (_screen_seed: the seed contract)
- validation/waveoptics_speed/profile_baseline_results.json (which part of
  screen generation dominates; attack that first)

LICENCE RULE: aotools is LGPL-3.0 and stays import-only. Do not copy its
code. The new generator derives from the Schmidt book equations and the
olb schmidt layer: the modified von Karman phase PSD is Schmidt,
DOI 10.1117/3.866274, Ch. 9, Eqs. (9.51) and (9.52), printed p. 161; the
Fourier-series screen is Eqs. (9.78) to (9.80), printed pp. 166-167; the
subharmonics are Eq. (9.81), printed p. 169, from Lane, Glindemann and
Dainty, DOI 10.1088/0959-7174/2/3/003.

THE THREE STRUCTURAL WINS (verified in the aotools source):
1. The sqrt-PSD filter depends on (n, dx, L0, l0) only; r0 enters as the
   scalar factor r0^(-5/6). Compute the filter ONCE per grid and scale it
   per screen.
2. Each subharmonic mode exp(2 pi i (fx x + fy y)) is a separable outer
   product of two length-n vectors. Precompute the per-level basis vectors
   once per grid and build the 27-mode sum with small matrix products, not
   27 full-grid complex exponentials.
3. Optional: one complex Gaussian draw + one FFT gives TWO independent real
   screens (real and imaginary parts). Offer a batch call that fills a
   whole trial stack and uses the pairing internally. Keep the seed
   deterministic: one SeedSequence-derived rng per stack.

SHAPE (ponytail): the laziest good API is a small factory, for example
`ScreenFactory(n, pixel_m, L0_m=inf, l0_m=1e-6)` with
`make(r0_m, rng) -> screen` and `make_stack(r0_m_array, rng) -> list`.
Cache the filter and the subharmonic basis inside it. Add an optional
dtype switch (complex64 internally, float32 out) and MEASURE its error
before you recommend it. Then wire an OPT-IN argument into
propagate_turbulent_scenario and propagate_turbulent_field (for example
screen_generator="aotools" (default) | "olb"), threading the same
_screen_seed integers into a numpy default_rng per screen. Document in the
docstring that the two generators give different draws for the same seed.

ACCEPTANCE (all must pass):
1. The new generator passes structure-function tests equal to the screens.py
   self-check cases 2, 3, and 5 (D_phi ratio inside [0.85, 1.02] over
   r/r0 = 0.3 to 1.6; the subharmonics are necessary; two screens add as
   r0^(-5/3)). Add these to the module self-check, gated the same way.
2. Statistical equivalence: 200 trials of the space downlink at 30 deg,
   rapid preset, old vs new generator (different seeds are fine). The mean
   collected power and the aperture sigma2_I agree inside the Monte-Carlo
   error bars (estimate the error of each statistic from the sample; state
   the bars in the output).
3. A speed table: seconds per screen, old vs new, n = 512 to 4096, and
   seconds per full stack for the P0 cases.
Save the evidence as validation/waveoptics_speed/screen_generator_check.py
with a results JSON and a run log, following the validation/ pattern.

DO NOT silently change any default. DO NOT extend TurbTrial.
```

### P2 — coarse screens and the beam-following grid (measure, do not wire)

```text
You work in the repository D:\repos\optical_link_budget (package `olb`).
Model guidance: be the laziest engineer that still works correctly. Borrow
the shared kernels. No speculative abstraction. All documentation, comments,
and commit messages use ASD-STE100 Simplified Technical English (see
CONVENTIONS.md). Every equation cites its source by DOI. Run modules with
`python -m ...` from the repository root.

TASK. Two measured experiments on the fidelity-2 grid cost. Both are
REPORT-ONLY: no production wiring. Each has a kill criterion, and a negative
result is a valid deliverable.

READ FIRST:
- olb/waveoptics/turbulence/screens.py — its docstring FORBIDS the naive
  coarse-then-interpolate screen, with the reason: the coarse screen is
  band-limited above its own Nyquist frequency, so it loses the
  Fresnel-scale structure sqrt(lambda z) that builds the scintillation
  (Schmidt, DOI 10.1117/3.866274, Sec. 9.4, printed p. 172). Experiment (a)
  MEASURES that claim; it does not assume it away, and it does not ignore
  it.
- olb/waveoptics/turbulence/sampling.py (the one-grid sizer and its rules)
- olb/waveoptics/turbulence/run.py and splitstep.py (the trial body)
- olb/waveoptics/schmidt/fresnel.py (partial_propagations,
  two_step_fresnel: the book per-plane-pitch chain, Ch. 8, Eq. (8.18))
- olb/waveoptics/schmidt/sampling.py (constraint1_max_delta2 ...
  constraint4_min_n, check_sampling, partial_grid_spacing,
  partial_max_step, partial_plane_count)
OWNER DECISION 2026-08-29: the schmidt layer may be borrowed and wired for
this speed work; the old validation-only rule is relaxed here.

EXPERIMENT (a): coarse screens + interpolation.
Script: validation/waveoptics_speed/coarse_screen_experiment.py.
Generate each screen at n/f pixels for f in {2, 4, 8}, interpolate to n
(test both FFT zero-padding and bicubic), and run the split step on the
full grid. Reference: the standard full-resolution run, same preset, same
trial count. Cases: terrestrial 2 km Cn2 5e-15 standard, space downlink
30 deg rapid. Measure, with about 200 trials each: mean collected power,
aperture sigma2_I, smf_eta (terrestrial), the phase structure function of
one screen against D_phi = 6.88 (r/r0)^(5/3) (Fried,
DOI 10.1364/JOSA.56.001372), and the wall-time saving. Also test the
hybrid: coarse low-frequency part + full-resolution high-frequency part,
and state whether it buys anything the subharmonics do not already give.
KILL CRITERION (state it in the script docstring before you run): the
approach dies if sigma2_I moves by more than 5 percent at every
configuration that saves time, or if no configuration beats the fast
generator of screen_generator_check.py (if that exists) at equal accuracy.

EXPERIMENT (b): a grid that follows the beam.
Script: validation/waveoptics_speed/beam_grid_experiment.py.
Step 1: VERIFY the claim that the space slab gains nothing from beam-sized
screens, because its input plane wave fills the grid. Measure the
irradiance support at each screen plane of a space 30 deg standard run and
report it.
Step 2 (terrestrial): the launch beam is small near the transmitter.
Prototype beam-extent-sized screens: at each plane, generate an m x m
screen at the SAME pixel pitch that covers the beam support plus a guard,
embed it in the full grid (zero phase outside), and measure speed vs error
against the full-screen reference (same statistics as experiment (a)).
This saves screen generation only, not FFT cost; say so in the report.
Step 3 (the real prize, gap S-14): a per-plane-pitch chain. Build a
prototype split step on schmidt fresnel.partial_propagations, with the
pitch schedule from sampling.partial_grid_spacing and the plane geometry
checked by constraints 1 to 4 (check_sampling). The screens then live on
per-plane pitches (screen_r0 stays valid; the pitch enters phase_screen
directly). Compare the total pixel-operations count and the wall time
against the one-flat-grid reference at equal accuracy, for the terrestrial
case AND the space slab.
KILL CRITERION: a variant dies if it cannot hold mean power inside 0.1 dB
and sigma2_I inside 10 percent of the reference while saving time.

OUTPUT: the two scripts, one results JSON and one run log each
(validation/ pattern), and a RECOMMENDATION section at the top of each
script docstring: wire it, or bury it, and why, in three sentences or
fewer. DO NOT change production code. DO NOT extend TurbTrial.
```

### P3 — the parallel scaling study

```text
You work in the repository D:\repos\optical_link_budget (package `olb`).
Model guidance: be the laziest engineer that still works correctly. Borrow
the shared kernels. No speculative abstraction. All documentation, comments,
and commit messages use ASD-STE100 Simplified Technical English (see
CONVENTIONS.md). Every equation cites its source by DOI. Run modules with
`python -m ...` from the repository root. The platform is Windows 11
(spawn, not fork).

TASK. Measure how fidelity-2 turbulent trials scale across workers, and
recommend an execution mode per case. REPORT-ONLY: no production change
except (optionally) a measured default for Threader max_workers, and only
if the data supports one.

READ FIRST:
- olb/waveoptics/threader.py (the existing thread pool; pocketfft releases
  the GIL, so threads already help)
- olb/waveoptics/turbulence/run.py (propagate_turbulent_scenario takes
  threader=; run_one is the unit of work; the seed contract makes any
  execution order give identical trials)
- validation/waveoptics_speed/profile_baseline_results.json (the per-part
  costs; interpret your curves against it)

BUILD: validation/waveoptics_speed/scaling_study.py

MEASUREMENTS (cases: terrestrial 2 km standard; space downlink 30 deg
rapid AND standard; about 32 trials per point, fixed seed):
1. THREADS. The Threader scaling curve: workers = 1, 2, 4, ... up to the
   core count. Report trials/second and parallel efficiency. Find where
   the curve saturates and say whether memory bandwidth is the plausible
   limit (compare against the raw fft2 numbers of the P0 JSON).
2. PROCESSES. A ProcessPoolExecutor variant: a module-level worker function
   that rebuilds the read-only setup once per process (an initializer) and
   runs trials by index. Measure the spawn + pickle overhead separately
   from the steady-state rate. Guard everything under
   if __name__ == "__main__". Check that the scenario dataclasses pickle;
   report if they do not.
3. BATCHED FFT. A prototype batched split step: stack B trials into
   (B, N, N) complex arrays and run the hops with numpy fft2 over the last
   two axes (numpy vectorises the leading axis), and with
   scipy.fft.fft2(..., axes=(-2, -1), workers=W). Generate the screens per
   hop on the fly to cap memory. Report trials/second against B and W,
   and a memory budget: fields are N^2 complex128 = 64 MB at N = 2048, and
   the budget must stay honest about the screen arrays too.
4. ONE cross-check: any parallel mode must give the SAME trial statistics
   as the serial loop for the same seed (the seed contract). Assert it.

OUTPUT: validation/waveoptics_speed/scaling_study_results.json, a run log,
and a RECOMMENDATION table printed and repeated in the script docstring:
per case, the best mode, the worker count, and the speed-up over one
worker. If the data supports a better Threader default on this machine,
say so in one sentence; do not change it yourself.
```

### P4 — caching and the campaign execution design (run last)

```text
You work in the repository D:\repos\optical_link_budget (package `olb`).
Model guidance: be the laziest engineer that still works correctly. Borrow
the shared kernels. No speculative abstraction. All documentation, comments,
and commit messages use ASD-STE100 Simplified Technical English (see
CONVENTIONS.md). Every equation cites its source by DOI. Run modules with
`python -m ...` from the repository root.

TASK. Decide what the fidelity-2 layer caches, design how a verification
campaign executes, and build the SMALLEST useful cache. This task runs
LAST: it reads the results of the earlier speed tasks.

READ FIRST:
- validation/waveoptics_speed/*.json and the run logs (the measured costs
  and the recommended execution modes)
- olb/waveoptics/turbulence/run.py (the vacuum baseline is recomputed on
  every propagate_turbulent_scenario call; the seed contract: trial k is
  bit-identical for one seed, and a longer run repeats the trials of a
  shorter run — this makes cached runs EXTENDABLE)
- olb/models/waveoptics.py (run_fidelity2, Fidelity2Bundle: the budget
  precompute pattern; the budget never runs the sim)
- docs/waveoptics-efficiency-plan.md (this plan; you append to it)

DECIDE, with the measured numbers, for each cache level:
1. Per-grid filter and subharmonic basis: in-process, inside the screen
   factory (P1 likely built this; confirm and say so).
2. The vacuum baseline per (scenario, grid, plan): it costs about one
   trial per call. Cache it in-process keyed on the grid + plan, and
   measure the saving for a repeated-call session.
3. Whole TurbWaveResult runs on disk, keyed by (a scenario content hash,
   seed_entropy, preset, generator version). The seed contract makes this
   cache EXTENDABLE: a stored 200-trial run serves any request for the
   first 200 trials, and a request for 500 computes only trials 200..499.
   Design the key and the storage (a JSON of scalars is tiny; say why the
   fields suffice).
4. Full screen stacks on disk: size it honestly (for example 9 screens x
   200 trials x 2048^2 float64 is about 57 GB) and rule on it.

THEN write the campaign execution design as a new section 8 ("Findings and
the execution design") of docs/waveoptics-efficiency-plan.md: how a
verification campaign runs fidelity-2 — precompute the Fidelity2Bundle
once, reuse it across budgets, pick the execution mode from the P3 table,
extend cached runs instead of recomputing, and store the evidence in the
validation/ pattern. Keep it under two pages.

BUILD the smallest piece that pays: by the numbers above, that is most
likely cache level 3 (the extendable disk cache of scalar results) plus
level 2 (the in-process vacuum baseline). Put the disk cache in a small
module (for example olb/waveoptics/turbulence/cache.py) with a load-or-run
entry point around propagate_turbulent_scenario. Opt-in, off by default.
Include a module self-check. DO NOT extend TurbTrial. DO NOT change any
budget default. Any wiring into the budgets stays owner-gated; say so in
the docstring.
```

## 8. Findings and the execution design

Status: DONE, 2026-08-29 (P4). This section reads the results of P0, P1, and
P3, decides what the layer caches, and records the campaign execution design.
The P4 disk cache (`cache.py`) was built on 2026-08-29 and RETIRED on
2026-09-04: `Campaign` (`olb/waveoptics/turbulence/campaign.py`, level 5 below)
replaces it, and its content key lives on in
`olb/waveoptics/turbulence/fingerprint.py`. The measured numbers of level 3
stay in this section as the record. No budget calls either module.

### 8.1 Two defaults changed (commit e8c7f77) — campaign-relevant facts

- **The default screen generator is now "olb".** `propagate_turbulent_scenario`
  and `propagate_turbulent_field` build screens with the cached `ScreenFactory`
  (P1), not aotools. "aotools" is the opt-in reference path. The two draw
  DIFFERENT atmospheres for the same seed, so the generator name is part of
  every cache key.
- **`Threader` defaults to `min(16, cores)` workers.** The scaling study (P3)
  saturates near 8 to 16 workers; more workers only add contention.

### 8.2 The five cache levels

| Level | What | Where | Measured cost / saving | Ruling |
|---|---|---|---|---|
| 1 | sqrt-PSD filter + subharmonic basis | in-process, `ScreenFactory` | build once per grid (n=2048: 0.19 s); then 7.5x per screen, 12x per stack vs aotools | BUILT (P1). Confirmed: `__init__` caches `_filt`, `_sub_filt`, `_E`; `make` scales by r0^(-5/6). |
| 2 | vacuum baseline per call | in-process | one trial per call (~0.35 s space standard propagation); a repeated request hits the session memo and skips the whole call | BUILT as a whole-result memo. See note. |
| 3 | whole `TurbWaveResult` runs | disk (JSON scalars) | 100-trial run 9.2 s; disk hit 21 ms (~400x); grow 100->150 computes only the new 50 (4.6 s, = a fresh 50, not a fresh 150); 150-trial file ~10 KB | RETIRED 2026-09-04. Level 5 replaces it; the key moved to `fingerprint.py`. |
| 4 | full screen stacks | disk | 9 screens x 200 trials x 2048^2 float64 ~= 57 GB for one case | RULED OUT. The screens regenerate in ~0.17 s per stack (P1); disk for tens of GB buys nothing. |
| 5 | the campaign store: the trial scalars PLUS the masked receive-field patch | disk (one `.npz` for each block, plus a JSON manifest) | a 1 m patch at a 5 mm pitch is ~320 KB for each trial, so ~3.2 GB for 10,000 trials; a new receive aperture, obscuration, detector or defocus then costs NO propagation (`recouple`, `recollect`) | BUILT 2026-09-04. `campaign.py`. It supersedes level 3. |

**Level 5 replaces level 3.** `olb/waveoptics/turbulence/campaign.py` holds
the blocks of ONE seeded native run (the runner's new `start_index`), so a
campaign slice is bit-identical to a native run, and it keeps the receive field.
`cache.py` seeded each block from a sub-seed and stored no field. It is DELETED
(2026-09-04, backlog 2-W4 DONE), with its self-check and
`validation/waveoptics_speed/cache_check.py`. Its content key moved unchanged to
`fingerprint.py`.

**Level 2, honestly.** The space vacuum baseline is recomputed INSIDE
`propagate_turbulent_scenario` (`olb/waveoptics/turbulence/run.py`). A
store around that runner cannot inject a precomputed baseline without editing
it, which is out of scope. So a `Campaign` does not cache the baseline alone.
Instead it sizes the grid and the plan ONE time (the manifest keeps them, so a
resumed campaign never re-sizes), and a `load` of stored blocks makes no call
at all. Only the cheap baseline split-step repeats for each new block, not the
grid sizing. A dedicated baseline cache would need a `baseline=` argument in
the runner, an owner-gated change, not taken here.

### 8.3 The content key (`fingerprint.py`), and the retired storage (level 3)

- **The key** is `cache_key(...)` in `olb/waveoptics/turbulence/fingerprint.py`,
  a SHA-256 of everything that changes a trial: `repr(scenario)` (dataclass
  repr, stable), a canonical geometry signature (the object has no stable
  repr), the preset, the base seed, the screen generator, `KEY_VERSION`, `L0_m`,
  the subharmonic switch, the Cn2 inputs (a callable is sampled at a fixed set
  of heights), the block size, and any caller grid/plan. A change to any input
  gives a new key. `Campaign` stores it in the manifest, and a mismatch raises.
  The key was born in `cache.py` and it moved here unchanged on 2026-09-04.
- **The retired storage (the record).** `cache.py` stored a run as a small
  JSON of five scalars for each trial (`collected_power`, `smf_eta`,
  `eta_turb`, `mmf_eta`, `wall_time_s`; ~10 KB for 150 trials), in fixed
  blocks that each ran as a self-contained runner call seeded from
  `SeedSequence(base_seed, spawn_key=(CACHE_VERSION, b))`. A disk hit was ~400x
  faster than a cold run, and a grow computed only the missing blocks (measured:
  grow 100->150 = a fresh 50). Its blocks were i.i.d. snapshots, but NOT the
  trials of one native single-seed run. The runner's `start_index` (2026-09-04)
  removed that limit, and `Campaign` (level 5) keeps the native seeding and the
  field, so the cache is deleted.

### 8.4 How a verification campaign runs fidelity 2

1. **Precompute the bundle once per scenario.** Build the `Fidelity2Bundle`
   with `olb.models.waveoptics.run_fidelity2` (the vacuum WaveResult plus the
   turbulent `TurbWaveResult`). The budgets never run the sim; they consume the
   bundle. Reuse ONE bundle across every budget of that scenario (uplink,
   downlink, SMF, MMF, and every availability quantile) — the reducer is cheap.
2. **Store and extend the turbulent run.** For a long or repeated campaign,
   fill the turbulent record through a `Campaign` with a fixed integer seed:
   `Campaign.run(n)` grows the stored blocks and computes only the missing
   ones, and `Campaign.load(n)` gives the `TurbWaveResult` for the reducer with
   no new call. A deeper tail (more trials for a deep-fade quantile) is one
   more `run(n)`. Keep the vacuum run as its own single deterministic call (it
   is cheap; it is not stored).
3. **Pick the execution mode from the P3 table.** For the "olb" generator, the
   P3 reading was that PROCESSES win by about 1.4x over threads. The fair rerun
   of 2026-09-04 REVISES that number to a wall-time TIE for one run; see
   Section 8.6. Use a `ProcessPoolExecutor` with
   8 to 16 workers: terrestrial and space rapid saturate at 16 (~17 trials/s),
   space standard at 8 (~5.9 trials/s). Threads (`Threader`, default 16) are the
   built-in, no-pickle fallback. Beyond 16 workers the rate falls (spawn and
   memory-bandwidth limits; the raw fft2 numbers in the P0 JSON confirm the
   band is the ceiling). The cache and any execution mode give the SAME trial
   statistics (the seed contract).
4. **Store the evidence in the validation/ pattern.** One script, one results
   JSON, one run log. `validation/waveoptics_speed/cache_check.py` measured the
   retired cache (a memo hit, a disk hit, a JSON file size), and it was DELETED
   with it (2026-09-04): those measurements do not map onto `Campaign`. The
   campaign self-check (`python -m olb.waveoptics.turbulence.campaign`) proves
   the grow-only-missing-block behaviour and the bit-identity with a native
   run, and `examples/waveoptics/campaign_demo.py` shows a resumed campaign.

### 8.5 Owner-gated follow-ups

- Whether a `Campaign` ever backs a budget by default. A budget already READS
  one: `wave=campaign` is a fidelity-2 wave record, and every Monte Carlo
  example uses it. No DEFAULT budget reads one, because fidelity 2 is opt-in.
- A start-index in the runner, for a single-seed tail extension that is
  bit-identical to a native run: DONE (2026-09-04, `start_index`; `Campaign`
  uses it).
- An in-process baseline cache, if the runner grows a `baseline=` argument.

### 8.6 P3 fair rerun (2026-09-04)

Script: `validation/waveoptics_speed/fair_scaling_rerun.py`. Log:
`fair_scaling_rerun_run.log`. The rerun corrects four defects of P3: it runs
200 trials for each point (P3 ran 32, so one straggler trial moved a point), it
pins the BLAS thread count to 1 in BOTH routes, it gives each process a chunk of
n/W trials instead of one trial for each task, and it counts the pool SPAWN
inside the wall time. It adds a pure-`fft2` process ceiling with NO Python work
between the transforms, as the machine limit.

| Case (grid) | Route | Best W | Wall, 200 trials | Steady ratio |
|---|---|---|---|---|
| space downlink 30 deg, standard (N=512) | threads | 16 | 31.92 s | — |
| space downlink 30 deg, standard (N=512) | processes | 16 | 32.20 s (spawn 4.22 s, steady 27.98 s) | 1.14x |
| terrestrial 2 km, standard (N=256) | threads | 8 | 7.77 s | — |
| terrestrial 2 km, standard (N=256) | processes | 12 | 7.50 s (spawn 2.57 s, steady 4.93 s) | 1.74x |

Processes over threads on WALL time: 0.99x on the space case and 1.04x on the
terrestrial case. So the two routes TIE for one run. Processes win 1.14x and
1.74x in STEADY state, and the Windows pool spawn (2.5 to 4.4 s) cancels that
win.

**The fft2 ceiling is the machine, not the GIL.** With no Python between the
transforms, the process speed-up plateaus at 5 to 6x on N=512 (5.36x at W=8,
5.77x at W=24) and at 8x on N=256 (8.40x at W=12). The threads reach 4.23x and
3.82x. So the plateau comes from the memory bandwidth and the hybrid P/E cores
of the i9-14900HX, and the GIL costs only on the small grid.

**The BLAS pin does nothing here.** `threadpoolctl` finds no BLAS pool in this
numpy build, and an alternating A/B of the environment pin against no pin is
flat inside the noise. A wired pin was tried and reverted.

**Noise.** The timing noise on this laptop is +/-10 percent, and it drifts down
over a session. The owner sees the same swing when the window moves to the
foreground (Windows 11 throttles the power of a background process on a hybrid
CPU). So only alternating A/B pairs count as evidence.

**The decision.** There is NO automatic selector between threads and processes:
a selector between two options that tie inside the noise is not worth a branch.
Threads stay the default of ONE run. A process pool pays only when it stays WARM
across many blocks, and that is the `Campaign`
(`olb/waveoptics/turbulence/campaign.py`, `Campaign.run(n, workers=W)`): ONE
pool for the whole call, and the blocks run serially inside each process, so the
parallelism lives at one level only.
