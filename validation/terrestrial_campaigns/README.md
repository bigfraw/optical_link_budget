# The terrestrial fidelity-2 campaign set

`run_campaigns.py` stores a backbone dataset of terrestrial (horizontal)
fidelity-2 snapshots, and it reports what the dataset costs. It is a RUNNER
only: it stores the trials and it measures the time, the memory and the disk.
It does NO distribution analysis. The owner deferred every analysis question;
they are listed at the end of this file.

## Why the dataset

The stored trials serve four open questions:

- **backlog 1-8, gate (b).** Is the terrestrial fade lognormal over this band?
  The band runs from firmly weak to strongly saturated.
- **The `rapid` preset.** Is it usable for a terrestrial link? Each case runs
  at `rapid` and at `standard`, on the same seed, so the two are comparable.
- **backlog 2-N2, beam filling.** A receive aperture that holds much of the
  beam breaks the plane-wave aperture-averaging fit. `cell.json` gives the
  vacuum beam radius w(L) and the captured fraction eta_fill of each aperture.
- **The single-mode-fibre coupling distribution** of a horizontal link, which
  no analytic Term of olb models as a fade today.

## The cells

Three path lengths, two turbulence strengths, two quality presets: twelve
campaigns. The numbers below come from `--dry-run` (collimated launch,
L0 = 25 m, seed 20260906, single precision).

| cell | preset | sigma_R^2 | rho_0 | w(L) | grid | side | pixel | screens |
|---|---|---|---|---|---|---|---|---|
| 2 km, 3e-15 | rapid | 0.213 | 5.07 cm | 19.7 cm | 1024 px | 0.945 m | 0.92 mm | 5 |
| 2 km, 1e-14 | rapid | 0.709 | 2.46 cm | 19.7 cm | 1024 px | 1.022 m | 1.00 mm | 5 |
| 5 km, 3e-15 | rapid | 1.142 | 2.93 cm | 49.3 cm | 1024 px | 2.494 m | 2.44 mm | 5 |
| 5 km, 1e-14 | rapid | 3.806 | 1.42 cm | 49.3 cm | 1024 px | 2.858 m | 2.79 mm | 10 |
| 10 km, 3e-15 | rapid | 4.069 | 1.93 cm | 98.7 cm | 1024 px | 5.331 m | 5.21 mm | 11 |
| 10 km, 1e-14 | rapid | 13.564 | 0.94 cm | 98.7 cm | 1024 px | 6.529 m | 6.38 mm | 35 |
| 2 km, 3e-15 | standard | 0.213 | 5.07 cm | 19.7 cm | 2048 px | 1.425 m | 0.70 mm | 9 |
| 2 km, 1e-14 | standard | 0.709 | 2.46 cm | 19.7 cm | 2048 px | 1.509 m | 0.74 mm | 9 |
| 5 km, 3e-15 | standard | 1.142 | 2.93 cm | 49.3 cm | 2048 px | 3.704 m | 1.81 mm | 9 |
| 5 km, 1e-14 | standard | 3.806 | 1.42 cm | 49.3 cm | 2048 px | 4.112 m | 2.01 mm | 10 |
| 10 km, 3e-15 | standard | 4.069 | 1.93 cm | 98.7 cm | 2048 px | 7.788 m | 3.80 mm | 11 |
| 10 km, 1e-14 | standard | 13.564 | 0.94 cm | 98.7 cm | 2048 px | 9.174 m | 4.48 mm | 35 |

`sigma_R^2` is the plane-wave Rytov variance (Andrews and Phillips, 2nd ed.
(2005), DOI 10.1117/3.626196, Ch. 8, Eq. (20), printed p. 264), `rho_0` the
coherence radius (Ch. 6) and `w(L)` the vacuum received beam radius (Ch. 4,
Eqs. (7) and (8), printed p. 87). The screen counts are the counts the shipped
planner gives TODAY. Do not read them as fixed: the planner changes, and the
script reads the count from the plan. Always run `--dry-run` for the current
numbers.

The link is the same in every cell: 1550 nm, a 5 mm launch waist through a
10 cm transmit aperture, a 10 cm receive aperture with a single-mode fibre at
`optimal_focus` (the Shaklan and Roddier coupling parameter a = 1.12,
DOI 10.1364/AO.27.002334), no extinction, and the fixed outer scale
L0 = 25 m (the owner decision of 2026-09-05, backlog 2-P5).

`--launch diverged` repeats a cell with the launch opened to 2.0e-4 rad, about
two times the diffraction divergence of the 5 mm waist. It is a spot check, and
it gets its own campaign roots (`..._diverged`).

## What a trial stores

Each trial holds:

- `collected_power`: the power inside the 10 cm receive aperture, as a
  fraction of the launched power. It holds the geometric spread too, because a
  terrestrial trial divides by the power after the transmit clip.
- `smf_eta`: the single-mode-fibre coupling efficiency of that 10 cm aperture,
  UNTRACKED, at `defocus_m = 0`.
- the complex64 receive-plane field on a disc of radius 5 cm, BEFORE the
  aperture clip. That radius equals the clip radius of the 10 cm aperture.

## What is post hoc

The stored field makes these a `recouple` or a `recollect` of the campaign,
with no new propagation:

- the 5 cm bucket (`camp.recollect(aperture_m=0.05)`) and the 5 cm fibre
  (`camp.recouple(SMF(optimal_focus=True), aperture_m=0.05)`);
- a central obscuration, another focal length, another mode field radius,
  another defocus, or another detector type;
- the TRACKED (aligned) coupler: set
  `detector.defocus_m = cell["apertures"][k]["curvature_focus_shift_m"]`, the
  focus shift that the received curvature causes (S. A. Self, Appl. Opt. 22
  (1983) 658, DOI 10.1364/AO.22.000658). `cell.json` holds that value for each
  aperture;
- the POINT scintillation index: read the CENTRE PIXEL of the stored patch, as
  `validation/outer_scale_tail/outer_scale_tail.py` (`_centre_irradiance`,
  line 365) does. CAUTION: at 10 km the pixel is 3.8 to 6.4 mm, so the "point"
  is a coarse point. The pixel of each cell is in `cell.json`.

**NO tip-tilt correction.** Fidelity 2 models no adaptive optics and no
tracking (backlog 2-AO), so every fibre number here is untracked and it holds
the full beam wander.

## The clamped grids

The 5 mm launch waist gives a beam radius of 49 cm at 5 km and 99 cm at 10 km.
The grid side must hold that beam plus the scattering cone, so the pixel count
the sizer wants passes the preset `n_max`. The sizer then KEEPS the side and it
takes a coarser pixel, and it says so. Five of the twelve cells carry that
warning today, for example:

    turbulent_grid: the pixel count wants 4096, but n_max is 1024. The grid
    keeps its side and takes a coarse pixel.

Every warning is stored in `cell.json` under `sizer.warnings`, next to
`sizer.pixels_per_r0`, `sizer.fresnel_pixels_min` and
`sizer.step_over_limit_max`. Read them before you read a long-path result.

## The storage

The stored disc is the patch, so the cost falls with the grid pitch. Measured
from the sized campaigns, for 2000 trials in each cell:

| cell | patch pixels | kB per trial | MB per 2000 trials |
|---|---|---|---|
| 2 km, 3e-15, standard | 16237 | 127 | 248 |
| 2 km, 1e-14, standard | 14441 | 113 | 220 |
| 2 km, 3e-15, rapid | 9233 | 72 | 141 |
| 2 km, 1e-14, rapid | 7869 | 61 | 120 |
| 5 km (four cells) | 1005 to 2393 | 8 to 19 | 15 to 37 |
| 10 km (four cells) | 193 to 545 | 1.5 to 4.3 | 3 to 8 |

The whole twelve-cell set at 2000 trials is about 0.9 GB. The 2 km standard
cells hold almost all of it, because their pixel is the smallest.

## How to run it

**Size everything and run nothing.** It writes `cell.json` next to each
campaign root and it prints the table:

    python -m validation.terrestrial_campaigns.run_campaigns --dry-run

**The smoke pass.** A few trials for each cell, in separate `_smoke` roots. It
reports the seconds for each NEW trial, the projected hours for 2000 trials at
that worker count, the peak python working set from a background poll, and the
sizer warnings. It also cross-checks the post-hoc `recouple` against the stored
`smf_eta`, and it checks that the 5 cm crops are finite and physical. It writes
`smoke_results.json` and `smoke.log`:

    python -m validation.terrestrial_campaigns.run_campaigns --smoke --workers 8

**The full run.** It stores 2000 trials for each cell in blocks of 50, and it
writes `run_<cell>.log` for each cell plus the shared `run_campaigns.log`:

    python -m validation.terrestrial_campaigns.run_campaigns --workers 8 --block-size 50

The run RESUMES: a block that already sits on disk is not recomputed, so a
killed run loses at most the blocks in flight. The script prints the campaign
root and the trials already on disk before it runs anything.

Other switches: `--cells 2km:3e-15:rapid 5km:1e-14:standard` picks the cells,
`--order table` keeps the table order (the default `cheap-first` runs the small
grids first), `--launch diverged` opens the launch, `--workers auto` lets the
campaign size its own pool (`--workers` takes an integer or `auto`), and
`OLB_TERRESTRIAL_CAMPAIGNS_ROOT` moves the whole store off this folder.

## The two speed opt-ins

Two settings make a trial faster. Each one is OFF by default:

| flag | default | the opt-in | the root suffix |
|---|---|---|---|
| `--fft-backend` | `numpy` | `scipy` | `_scipy` |
| `--screen-generator` | `olb` | `olb-lean` | `_lean` |

Both together give the suffix `_scipy_lean`, so the 5 km standard cell writes
into `campaigns/L5km_cn23e-15_standard_scipy_lean/` and it logs to
`run_L5km_cn23e-15_standard_scipy_lean.log`. `cell.json` records both settings
under `campaign.fft_backend` and `campaign.screen_generator`, next to the
campaign fingerprint it already holds. The `--dry-run` table shows the root of
each cell in its last column, so the suffix is visible before a run.

**THE FINGERPRINT RULE.** Each opt-in ENTERS the campaign fingerprint. So a
`Campaign` that reopens an existing root with different settings RAISES:

    ValueError: the campaign in ...\L2km_cn23e-15_rapid_scipy_lean was made
    with fft_backend='scipy', and this Campaign asks for fft_backend='numpy'.
    A stored campaign is ONE physics case. Use a new directory, or match the
    stored settings.

A campaign that takes an opt-in therefore needs a NEW directory. It cannot
resume a store that the default settings made.

**THE PHYSICS AGREES.** The opt-ins keep the same random stream and the same
equations. Measured in `validation/memory_cut/`, they move the collected power
and the coupling efficiency by about 6e-7 in relative terms, the rounding
level of single precision, the same level that the double-to-single switch
measured (`validation/precision/`).

**THE OWNER DECISION (2026-09-06).** The REMAINING standard cells run with
BOTH opt-ins on, so the 5 km and the 10 km standard cells carry
`_scipy_lean`. The other cells do NOT:

- the eight finished cells keep their old directories and their old manifests:
  all six `rapid` cells, and the two 2 km `standard` cells;
- the partial `L5km_cn23e-15_standard` (24 blocks, the default settings) stays
  on disk untouched as a cross-check. It is NOT resumed.

To run the remaining standard cells:

    python -m validation.terrestrial_campaigns.run_campaigns \
        --cells 5km:3e-15:standard 5km:1e-14:standard \
                10km:3e-15:standard 10km:1e-14:standard \
        --fft-backend scipy --screen-generator olb-lean \
        --workers 8 --block-size 50

## Run it on bigfraw

The desktop runs a windowless ssh process throttled by default. The rules are
the rules of `validation/campaign_resources/README.md`; this is the same recipe
for this module.

1. Chain remote commands with `;`, never `&&`. The login shell of
   `ssh desktop` is Windows PowerShell 5.1, where `&&` is a parse error.
2. Ship the code from a PowerShell prompt, not from a bash shell.
3. Launch through WMI, so the run outlives the ssh session. Redirect the
   stdout to a DIFFERENT file than `smoke.log` / `run_campaigns.log`: the
   script opens those itself, and a redirect onto the same file makes it
   fail with a permission error at the first line.

   ```
   $cmd = 'cmd /c "cd /d D:\repos\optical_link_budget && C:\Users\alexf\anaconda3\envs\olb\python.exe -u -m validation.terrestrial_campaigns.run_campaigns --smoke --workers 8 > validation\terrestrial_campaigns\smoke_launch.log 2>&1"'
   Invoke-CimMethod -ClassName Win32_Process -MethodName Create -Arguments @{CommandLine=$cmd}
   ```

   and the full form:

   ```
   $cmd = 'cmd /c "cd /d D:\repos\optical_link_budget && C:\Users\alexf\anaconda3\envs\olb\python.exe -u -m validation.terrestrial_campaigns.run_campaigns --workers 8 --block-size 50 > validation\terrestrial_campaigns\run_launch.log 2>&1"'
   Invoke-CimMethod -ClassName Win32_Process -MethodName Create -Arguments @{CommandLine=$cmd}
   ```

4. The priority boost is automatic. `Campaign.run(boost=True)` is the default,
   so it boosts the parent and every pool worker. The script adds nothing.
5. Give the pool enough blocks. The effective process count is
   `min(workers, ceil(n_trials / block_size))`. 2000 trials in blocks of 50
   gives 40 blocks, which is enough for any worker count up to 40. The script
   warns when the block count falls below the worker count.
6. 8 workers is the chosen default on bigfraw. Judge the load with
   `% Processor Utility` next to `% Processor Time`.
7. Kill orphan python processes before a relaunch. The campaign resumes from
   the blocks on disk.

## Deferred: the analysis

The owner deferred all of this. This script must NOT do it:

- the fade distribution: lognormal against gamma-gamma, over the sigma_R^2
  band from 0.21 to 13.6;
- the index split: the point index against the aperture-averaged index, and
  the effective averaging factor;
- the fade quantiles (p10, p5, p1) and their bootstrap intervals;
- the `rapid` against `standard` verdict;
- the fidelity-2 fibre coupling against `terrestrial_smf_coupling_term` and
  the walk-off Term `terrestrial_smf_walkoff_term`;
- the tracked-focus `recouple` at `curvature_focus_shift`.

## Results

Not yet run.
