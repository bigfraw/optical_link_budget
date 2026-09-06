# terrestrial_screen_count/

The terrestrial SCREEN-COUNT convergence sweep (backlog 2-TC1).

## Purpose

The turbulent planner takes its screen count from the Schmidt per-screen cap:

    N = max(min_screens, ceil(sigma_R^2 / sigma2_r_screen_max))

At the 10 km, `Cn2 = 1e-14` standard cell that cap asks for **35 screens**. The
owner reads that as far more than a uniform horizontal path needs.

The cap is a per-screen **thin-screen validity rule**, not a convergence result:
Schmidt, *Numerical Simulation of Optical Wave Propagation* (2010),
DOI 10.1117/3.866274, Listing 9.5, lines 37 and 38, printed p. 175, which the
book credits to Martin and Flatte, DOI 10.1364/AO.27.002111. The book caps the
log-amplitude variance at `rmax = 0.1`, and `sigma_R^2 = 4 sigma_chi^2`, so the
olb field `sigma2_r_screen_max = 0.4` IS that book cap. The book gives no
statement that a plan under the cap has converged statistics, and olb has no
terrestrial convergence sweep (WP7 measured a SLANT slab). So the count that a
horizontal path really needs is UNMEASURED. This study measures it.

## The cell is a CHOICE

`--path-km` and `--cn2` select any backbone cell of
`validation/terrestrial_campaigns/`. The constants and the scenario builder are
IMPORTED from that module, so the two studies cannot drift apart.

The reference screen count is NEVER hard-coded: the script asks the sizer for
the count of the UNMODIFIED `standard` preset of the selected cell, and it
labels the reference row `ref<count>`.

| Item | 10 km cell (the default) | 5 km cell |
| --- | --- | --- |
| path | 10 km horizontal, `Cn2 = 1e-14 m^(-2/3)` | 5 km horizontal, `Cn2 = 3e-15 m^(-2/3)` |
| reference count | 35 screens (the CAP binds) | 9 screens (the `min_screens` FLOOR binds) |
| grid | 2048 px, 9.174 m side, 4.48 mm pixel | 2048 px, 3.704 m side, 1.81 mm pixel |
| backbone root | `L10km_cn21e-14_standard_scipy_lean` | `L5km_cn23e-15_standard_scipy_lean` |
| reference trials | 2000 | 2000 |

The two cells ask the two halves of the question: is the CAP too high, and is
the FLOOR enough?

Everything else is common to both: no absorption; a 100 mm transmit terminal
with a 5 mm collimated waist at 1550 nm; a 100 mm receive aperture with
`SMF(optimal_focus=True)`; the `standard` preset; `L0 = 25 m` (the owner
decision of backlog 2-P5); `single` precision; seed `20260906`; patch radius
5 cm; block size 50.

## The file names

Every output name and every override root carries the CELL tag, so two cells
never overwrite each other:

| Cell | Override root | JSON | Log | Figure |
| --- | --- | --- | --- | --- |
| 10 km (default) | `L10km_cn21e-14_standard_n{n}_scipy_lean` | `screen_count_sweep_results.json` | `screen_count_sweep.log` | `figures/screen_count_sweep.png` |
| 5 km | `L5km_cn23e-15_standard_n{n}_scipy_lean` | `screen_count_sweep_L5km_cn23e-15_results.json` | `screen_count_sweep_L5km_cn23e-15.log` | `figures/screen_count_sweep_L5km_cn23e-15.png` |

The naming rule gives the DEFAULT cell an EMPTY tag, so the stored record of
the 10 km cell keeps the plain names it already has on disk and it stays
reproducible. Every other cell takes `_L<path>km_cn2<cn2>`. The override roots
always carried the cell tag, because they are built from the backbone
`cell_tag`, so those names do not move either.

## The override mechanism

`_plan_terrestrial` picks the count from the cap and then cuts the path into
slabs of EQUAL plane-wave Rytov weight (Andrews and Phillips, 2nd ed. (2005),
DOI 10.1117/3.626196, Ch. 8, Eq. (20), printed p. 264). A guard loop can raise
the count by one when the midpoint of the last slab overshoots the cap.

To force EXACTLY `n` screens of the SAME shape, the script asks
`turbulent_grid` for a preset copy with the cap effectively removed and the
floor set to `n`:

    dataclasses.replace(PRESETS["standard"], sigma2_r_screen_max=1e9,
                        min_screens=n)

The floor then binds, the guard loop never fires, and the script ASSERTS
`plan.z_m.size == n`.

The GRID comes from the UNMODIFIED `standard` sizing, and it is passed to
`Campaign(grid=..., plan=...)` explicitly, so the study moves the screens only.
The script also sizes the grid that the modified preset would pick and it
ASSERTS the two agree; they do at every count of both cells (the pixel count
clamps at `n_max = 2048` in every case), and the table column `grid held` says
so. The pixel count must be EQUAL and the side is compared with a relative
tolerance of `1e-9`, because the sizer rebuilds that float from the screen
plan and two plans of one cell can differ in the last bits.

Both `grid` and `plan` enter the campaign fingerprint, so each count is its own
store: `campaigns/L10km_cn21e-14_standard_n{n}/` (environment override
`OLB_SCREEN_COUNT_ROOT`).

## The reference-count bit-identity check

ON by default; `--no-ref-check` skips it, `--ref-check-trials` sets its size
(default 100). It runs BEFORE the sweep, so a broken pipeline stops the study
before it costs an hour.

The check builds the override campaign at `n` = the reference screen count,
into its OWN root (`..._n9_...` for the 5 km cell, `..._n35_...` for the
10 km cell), runs `--ref-check-trials` trials there, and asserts that the
`collected_power` and the `smf_eta` of those trials equal the reference values
EXACTLY (`==`). It first asserts that the override plan equals the production
plan (`z_m` and `r0_m`), because the `min_screens = n` route and the production
cap route could in principle place the screens differently; a difference gives
a clear message and no run. The maximum absolute difference is printed either
way.

A PASS proves that the grid pinning, the plan override, the two speed opt-ins
and the campaign reopen all reproduce the production path, so a delta at
another count is a SCREEN-COUNT effect and nothing else.

## The two speed opt-ins

`--fft-backend scipy` and `--screen-generator olb-lean` are the same two
opt-ins that the backbone runner takes, and the same suffix rule holds. They
go to BOTH the reference reopen and every override campaign, so the roots
become `campaigns/L10km_cn21e-14_standard_n{n}_scipy_lean/`. The `--dry-run`
table shows the root of each count in its last column.

Each opt-in ENTERS the campaign fingerprint, so a reopen with the wrong
settings RAISES a fingerprint mismatch. The settings must match the store the
reference was made with. The physics agrees at the rounding level of single
precision, about 6e-7 in the power and in the coupling efficiency
(`validation/memory_cut/`). See
`validation/terrestrial_campaigns/README.md` for the owner decision of
2026-09-06.

## The reference

The standard-preset campaign of the backbone run of the SAME cell, which on
`bigfraw` is the OPT-IN root
`validation/terrestrial_campaigns/campaigns/L10km_cn21e-14_standard_scipy_lean`
or `..._L5km_cn23e-15_standard_scipy_lean` (2000 trials, seed 20260906;
environment override `OLB_TERRESTRIAL_CAMPAIGNS_ROOT`). So the sweep must run
with
`--fft-backend scipy --screen-generator olb-lean`, and every override campaign
then takes the same settings, which keeps the comparison like for like. This
script **never runs the reference**: it reopens that store through
`run_campaigns.make_campaign`, so the settings and the fingerprint match, and
it reads the FIRST `--n-trials` trials. Those trials are bit-identical to a
shorter run of the same seed, because the runner seeds trial k off
`(entropy, k)`.

## How to run

From the repository root, with the environment python:

    # size every count, print the table and the projected cost, run nothing
    python -m validation.terrestrial_screen_count.screen_count_sweep --dry-run

    # the same for the 5 km cell
    python -m validation.terrestrial_screen_count.screen_count_sweep \
        --dry-run --path-km 5 --cn2 3e-15 \
        --fft-backend scipy --screen-generator olb-lean

    # the full sweep: 5, 10, 15 and 20 screens, 2000 trials each
    python -m validation.terrestrial_screen_count.screen_count_sweep \
        --workers 12 --block-size 50 \
        --fft-backend scipy --screen-generator olb-lean

    # read what is stored, and redo the analysis and the figure
    python -m validation.terrestrial_screen_count.screen_count_sweep \
        --analyse-only

The counts run CHEAPEST FIRST, and `Campaign.run` skips a stored block, so the
sweep is resumable after a kill. The priority boost of
`olb.waveoptics.priority` is automatic inside `Campaign.run`.

The DRY RUN projects the cost from a MEASURED backbone anchor of the cell: 6.0
s for one trial at 35 screens on 8 workers for the 10 km cell (the backbone
smoke run), and 1.12 s for one trial at 9 screens on 12 workers with the two
opt-ins for the 5 km cell (the backbone full run). A trial is close to linear
in the screen count, so the projection scales by `n / screens` and by the
worker count. A cell that has no anchor takes the 5 km one, and its projection
is a rough guide only.

### The 5 km launch line on bigfraw

    $cmd = 'cmd /c "cd /d D:\repos\optical_link_budget && C:\Users\alexf\anaconda3\envs\olb\python.exe -u -m validation.terrestrial_screen_count.screen_count_sweep --path-km 5 --cn2 3e-15 --counts 5 10 15 20 --n-trials 2000 --workers 12 --fft-backend scipy --screen-generator olb-lean > validation\terrestrial_screen_count\sweep_launch_L5km_cn23e-15.log 2>&1"'
    Invoke-CimMethod -ClassName Win32_Process -MethodName Create -Arguments @{CommandLine=$cmd}

The WMI launch detaches the process from the ssh session, so the run survives
a disconnect.

### The 10 km bigfraw queue note

The sweep is QUEUED behind the terrestrial backbone run on `bigfraw`: it is
launched by a `Wait-Process` wrapper that waits for the backbone python process
and then starts this module. The launch command:

    $cmd = 'powershell -NoProfile -Command "Wait-Process -Id <backbone PID> -ErrorAction SilentlyContinue; cd D:\repos\optical_link_budget; & C:\Users\alexf\anaconda3\envs\olb\python.exe -u -m validation.terrestrial_screen_count.screen_count_sweep --workers 12 --fft-backend scipy --screen-generator olb-lean > validation\terrestrial_screen_count\sweep_launch.log 2>&1"'
    Invoke-CimMethod -ClassName Win32_Process -MethodName Create -Arguments @{CommandLine=$cmd}

`<backbone PID>` is the process id that the backbone WMI launch returned
(3980 on 2026-09-06). The wrapper sleeps until that process ends, then it
runs the sweep in the same detached way.

## The two verdicts

Each case measures four quantities from the SAME stored trials:

| Quantity | What it is |
| --- | --- |
| `P10cm` | the collected power of the 100 mm receive bucket |
| `point` | the irradiance of the CENTRE PIXEL of the receive grid |
| `smf_eta` | the single-mode-fibre coupling efficiency |
| `P5cm` | the collected power of a 50 mm bucket, a post-hoc `recollect` |

For each one the script gives the scintillation index `var(x)/mean(x)^2`
(Andrews and Phillips, DOI 10.1117/3.626196, Ch. 8) and the fades p10, p5 and
p1 of `-10 log10(x / <x>)`, each with a 1000-resample bootstrap 68 percent
interval. Loss is positive dB, so a fade is a positive number.

Each count gets BOTH verdicts, because they answer different questions.

**1 RESOLUTION (inside the reference noise).** A count is CONVERGED when, for
every one of the four quantities:

- the p5 delta against the reference sits inside the bootstrap half-width of
  the reference, AND
- the p1 delta does too, AND
- the index ratio sits inside `1 +/- 0.10`.

Else the count is NOT CONVERGED, and the verdict line names the quantity and
the failing test. This says whether the run can SEE a difference, not whether
the difference matters.

**2 TOLERANCE (`--tolerance-db`, default 1.0 dB).** A count is WITHIN the
tolerance when every `|p5|` and `|p1|` delta of the 10 cm bucket AND of the SMF
coupling is under that many dB. Those two are the two RECEIVER kinds of a real
terrestrial link. This says whether the difference matters to a link budget.

The bootstrap bar is the noise of the REFERENCE alone. Each count is an
estimate of the same size, so the noise on a DELTA is about `sqrt(2)` times
that, and the reference-bar block prints both numbers, one line each.

A NOT CONVERGED verdict at 4 trials means nothing: the bootstrap bars are then
wider than the whole effect. Read the verdicts at 2000 trials.

The `point` quantity is the sharpest test, because nothing averages a point.
The `P10cm` bucket holds much of the beam on this path, so its index is small
and it is the least sensitive.

## Outputs

| File | What it holds |
| --- | --- |
| `screen_count_sweep_results.json` | every number, plus the `z_m` and the `sigma2_r` of each plan, and the `reference_check` record |
| `screen_count_sweep.log` | the printed lines, written line by line |
| `figures/screen_count_sweep.png` | the index and the p5 / p1 fade against the screen count, for `P10cm`, `point` and `smf_eta`. The reference is a horizontal band. |

A cell other than the default takes the same three names with the cell tag in
them. See "The file names" above.

## The worker plateau after the memory cut (2026-09-06)

Measured on bigfraw on ONE fixed configuration, the 5 km / Cn2 = 3e-15
standard reference plan (9 screens, 2048 px, the opt-ins `scipy` +
`olb-lean`), one campaign in `campaigns_workers/` with blocks of 10 (so the
pool is never short of blocks), grown by 300 trials at each worker count.
The record is `workers_L5km_cn23e-15.log` and `resources_L5km_cn23e-15.csv`
(a 30 s sampler: utility, processor time, RAM, python count, working set).

| workers | wall for 300 trials | s/trial | utility mean | RAM used | python working set |
|---|---|---|---|---|---|
| 8 | 348 s | 1.16 | 62 % | 20 GB | 5.0 GB |
| 12 | 329 s | 1.10 | 83 % | 22 GB | 6.9 GB |
| 16 | 318 s | 1.06 | 102 % | 24 GB | 9.9 GB |
| 20 | 318 s | 1.06 | 103 % | 26 GB | 12.1 GB |

The curve is FLAT: 8 to 20 workers spans 9 percent, and 20 equals 16. The
2048 px pool is memory-bandwidth bound after the cut as before it, so more
processes add utility, RAM and commit, not trials. THE SETTING OF RECORD IS
12 WORKERS for a 2048 px cell. The 24 and 28 stages were NOT run, on purpose:
each worker COMMITS about 2.2 GB (its allocation high-water mark, kept on
the heap) while it touches 0.6 GB, and the commit limit of the box is 61 GB
(32.5 GB RAM + a 28.7 GB page file) with about 15 GB taken by the rest of
the desktop, so 24 workers would have exhausted the commit charge and failed
allocations system-wide. A larger page file raises that limit at no runtime
cost (the reserved pages are never written), but the plateau says it is not
worth it. Two gotchas of the measurement: the sweep script's "s/trial (this
call)" line divided by the trials on disk, so a resumed stage read low (fixed:
the rate now counts the new trials), and the script's per-call analysis
rebuilt the full grid for every stored trial between stages, minutes of one
core that a Task Manager reading caught as "5 percent" (fixed: `--run-only`).

## Results, the 10 km cell

Run on bigfraw, 2026-09-06, queued behind the backbone run: counts 5 / 10 /
15 at 1000 trials each (12 workers, the opt-ins `scipy` + `olb-lean`, roots
`L10km_cn21e-14_standard_n{5,10,15}_scipy_lean`), against the first 1000
trials of the 35-screen backbone cell `L10km_cn21e-14_standard_scipy_lean`.
The owner stopped the sweep before the 20-screen count. The record is
`screen_count_sweep.log`, `screen_count_sweep_results.json`,
`figures/screen_count_sweep.png`, and the two launch logs.

| count | 10 cm bucket p5 / p1 vs 35 [dB] | bucket index ratio | 5 cm bucket index ratio | centre-pixel index ratio | verdict |
|---|---|---|---|---|---|
| 5 | -0.62 / -0.71 | 0.90 | 0.70 | 0.44 | NOT CONVERGED |
| 10 | -0.53 / -0.65 | 1.06 | 0.84 | 0.58 | NOT CONVERGED |
| 15 | -0.46 / -0.25 | 0.89 | 0.65 | 0.40 | NOT CONVERGED |

Reference bootstrap half-widths (68 percent): bucket p5 0.24 dB, p1 0.31 dB,
index 8 percent; 5 cm bucket index 21 percent; centre pixel index 42 percent.

VERDICT, BY RECEIVER KIND (the owner's reading, 2026-09-06). The bars above
are the bootstrap noise of the REFERENCE alone; each count is a 1000-trial
estimate too, so the noise on a DELTA is about sqrt(2) larger: +-0.34 dB at
the bucket p5, +-1.05 dB at the SMF p5, +-3.4 dB at the SMF p1. Read this way:

- SMF receiver (the link of interest): the p5 deltas -0.49 / +0.68 / +0.31 dB
  and the p1 deltas -0.73 / +0.88 / +0.04 dB carry no sign pattern and sit
  inside the noise. From 5 to 35 screens there is NO detectable screen-count
  effect on the fibre fade at 1000 trials. The fibre fade in saturation is
  tilt and low-order phase, the large scales that few screens already carry;
  the count changes the small-scale amplitude structure. So 10 screens is
  defensible for a fibre-coupled terrestrial link on this cell, and the cap
  count buys nothing measurable. The resolution is about +-1 dB at p5, so this
  rules out a large effect, not a small one; +-0.5 dB needs about four times
  the trials.
- Bucket receiver: a small, consistently OPTIMISTIC bias. The 10 cm bucket p5
  reads 0.62 / 0.53 / 0.46 dB less fade at 5 / 10 / 15 screens (1.4 to 1.8
  sigma each, the same sign every time), the 5 cm bucket index 0.70 / 0.84 /
  0.65 of the reference, the centre-pixel index 0.44 / 0.58 / 0.40. Within a
  1 dB tolerance 10 or 15 screens pass; a small sensor or a point statistic
  pays more.
- Not tested: whether 35 itself is converged (a 50 to 70 screen run, about
  1 h). The "CONVERGED" rule in the log means "inside the reference noise",
  a resolution statement, not a tolerance.

That reading is the reason the script now prints BOTH verdicts and the
`sqrt(2)` delta bar: the owner had to do both by hand for the table above.

The five other standard cells of the backbone sit at the `min_screens` floor
(9 to 11), so the cap bound only this cell. The count rule should key on the
receiver kind and a dB tolerance (2-I3), not on the cap alone.

## Results, the 5 km cell

NOT RUN YET. This cell sits at the `min_screens` FLOOR of 9 screens, so it asks
the other half of the question: is the floor enough? The plan is the counts
5 / 10 / 15 / 20 at 2000 trials each, against ALL 2000 trials of the backbone
reference `L5km_cn23e-15_standard_scipy_lean`, with the reference-count
bit-identity check at n = 9 first. The projected cost is about 3.5 h of pool
time at 12 workers (0.35 + 0.69 + 1.04 + 1.38 h), plus about 2 minutes for the
check. The 2000 trials, against the 1000 of the 10 km run, are what make the
p1 fade resolve.
