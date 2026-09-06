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

## The cell

Exactly the backbone cell of `validation/terrestrial_campaigns/`. The constants
and the scenario builder are IMPORTED from that module, so the two studies
cannot drift apart.

| Item | Value |
| --- | --- |
| path | 10 km horizontal, `Cn2 = 1e-14 m^(-2/3)`, no absorption |
| transmit | 100 mm terminal, 5 mm waist, collimated, 1550 nm |
| receive | 100 mm aperture, `SMF(optimal_focus=True)` |
| preset | `standard` (grid 2048 px, 9.174 m side, 4.48 mm pixel) |
| outer scale | `L0 = 25 m` (the owner decision of backlog 2-P5) |
| precision | `single`; seed `20260906`; patch radius 5 cm; block size 50 |

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
ASSERTS the two agree; they do at every count of this cell (the pixel count
clamps at `n_max = 2048` in both cases), and the table column `grid held` says
so.

Both `grid` and `plan` enter the campaign fingerprint, so each count is its own
store: `campaigns/L10km_cn21e-14_standard_n{n}/` (environment override
`OLB_SCREEN_COUNT_ROOT`).

## The reference

The 35-screen campaign of the backbone run,
`validation/terrestrial_campaigns/campaigns/L10km_cn21e-14_standard`
(2000 trials, seed 20260906; environment override
`OLB_TERRESTRIAL_CAMPAIGNS_ROOT`). This script **never runs it**: it reopens
that store through `run_campaigns.make_campaign`, so the settings and the
fingerprint match, and it reads the FIRST `--n-trials` trials. Those trials are
bit-identical to a shorter run of the same seed, because the runner seeds trial
k off `(entropy, k)`.

## How to run

From the repository root, with the environment python:

    # size every count, print the table and the projected cost, run nothing
    python -m validation.terrestrial_screen_count.screen_count_sweep --dry-run

    # the full sweep: 5, 10, 15 and 20 screens, 1000 trials each
    python -m validation.terrestrial_screen_count.screen_count_sweep \
        --workers 12 --block-size 50

    # read what is stored, and redo the analysis and the figure
    python -m validation.terrestrial_screen_count.screen_count_sweep \
        --analyse-only

The counts run CHEAPEST FIRST, and `Campaign.run` skips a stored block, so the
sweep is resumable after a kill. The priority boost of
`olb.waveoptics.priority` is automatic inside `Campaign.run`.

### The bigfraw queue note

The sweep is QUEUED behind the terrestrial backbone run on `bigfraw`: it is
launched by a `Wait-Process` wrapper that waits for the backbone python process
and then starts this module. The launch command:

    $cmd = 'powershell -NoProfile -Command "Wait-Process -Id <backbone PID> -ErrorAction SilentlyContinue; cd D:\repos\optical_link_budget; & C:\Users\alexf\anaconda3\envs\olb\python.exe -u -m validation.terrestrial_screen_count.screen_count_sweep --workers 12 > validation\terrestrial_screen_count\sweep_launch.log 2>&1"'
    Invoke-CimMethod -ClassName Win32_Process -MethodName Create -Arguments @{CommandLine=$cmd}

`<backbone PID>` is the process id that the backbone WMI launch returned
(3980 on 2026-09-06). The wrapper sleeps until that process ends, then it
runs the sweep in the same detached way.

## What CONVERGED means

Each case measures four quantities from the SAME stored trials:

| Quantity | What it is |
| --- | --- |
| `P10cm` | the collected power of the 100 mm receive bucket |
| `point` | the irradiance of the CENTRE PIXEL (the pixel is 4.48 mm here) |
| `smf_eta` | the single-mode-fibre coupling efficiency |
| `P5cm` | the collected power of a 50 mm bucket, a post-hoc `recollect` |

For each one the script gives the scintillation index `var(x)/mean(x)^2`
(Andrews and Phillips, DOI 10.1117/3.626196, Ch. 8) and the fades p10, p5 and
p1 of `-10 log10(x / <x>)`, each with a 1000-resample bootstrap 68 percent
interval. Loss is positive dB, so a fade is a positive number.

A count is **CONVERGED** when, for every one of the four quantities:

- the p5 delta against the 35-screen reference sits inside the bootstrap
  half-width of the reference, AND
- the p1 delta does too, AND
- the index ratio sits inside `1 +/- 0.10`.

Else the count is **NOT CONVERGED**, and the verdict line names the quantity
and the failing test. A NOT CONVERGED verdict at 4 trials means nothing: the
bootstrap bars are then wider than the whole effect. Read the verdict at 1000
trials.

The `point` quantity is the sharpest test, because nothing averages a point.
The `P10cm` bucket holds much of the beam on this path, so its index is small
and it is the least sensitive.

## Outputs

| File | What it holds |
| --- | --- |
| `screen_count_sweep_results.json` | every number, plus the `z_m` and the `sigma2_r` of each plan |
| `screen_count_sweep.log` | the printed lines, written line by line |
| `figures/screen_count_sweep.png` | the index and the p5 / p1 fade against the screen count, for `P10cm`, `point` and `smf_eta`. The reference is a horizontal band. |

## Results

**Not yet run.** The sweep is queued behind the terrestrial backbone campaign
on `bigfraw`. Fill this section from `screen_count_sweep.log` when it ends.

If 10 screens sit on the 35-screen line inside the Monte Carlo error, the cap
is over-conservative on a uniform horizontal path, and a terrestrial preset
rule (a count from a convergence table, not from the cap) replaces it. See
backlog 2-TC1 and 2-I3.
