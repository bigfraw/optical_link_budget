# HANDOVER: the diverged fidelity-2 uplink (2026-09-23)

> STATUS (later on 2026-09-23): DONE. The blocker of Section 3c is RESOLVED
> and the plan of Section 4 ran WITHOUT the Dios arm (owner decision). See
> [../uplink_divergence/README.md](../uplink_divergence/README.md).

Read this first in the next session. It records the state of the work, the
findings, the open problem and the next plan. The detail of each finding is in
[README.md](README.md).

## 1. The state of the repository

- Branch: `uplink-divergence`, made from LOCAL `main` (2defbf5). Local `main`
  is 18 commits AHEAD of `origin/main` (the `waveoptics-pointahead` merge is
  not pushed). The branch is NOT pushed.
- NOTHING IS COMMITTED. The working tree holds:
  - `olb/waveoptics/turbulence/run.py` and `campaign.py`: the conjugate fix
    (section 3).
  - `CLAUDE.md`, `docs/physics.md`, `docs/api-waveoptics.md`: the overlap
    formula without the conjugate.
  - `validation/README.md`: the `divergence_sampling/` index entry.
  - `validation/divergence_sampling/` (untracked): the study script, its
    README, `figures/`, `diagnostics/` (three rough scripts) and this file.
- `validation/temporal_screens/` shows as untracked on this branch. It is the
  leftover data of the `waveoptics-temporal` branch (campaigns, logs). Do NOT
  commit it and do NOT delete it.
- The `waveoptics-temporal` branch is untouched.
- The first job of the next session: ask the owner, then commit the fix and the
  study on `uplink-divergence` (two commits: the fix with its docs, then the
  validation study).

## 2. The question that started it

Can a fidelity-2 uplink take a DELIBERATE beam divergence? The owner asked how
the divergence enters, because the downlink slab starts from a plane wave.

The answer: the downlink slab NEVER sees a divergence. The uplink divergence
enters in two places, both from `Transmitter.divergence_rad` of the GROUND
terminal (there is no separate fidelity-2 setting):

1. The analytic geometric Term (`geometric_loss_term`). It is correct at every
   divergence.
2. The ground transmit mode `psi_tx` (`_ground_transmit_mode`): a Gaussian at
   the virtual waist, moved to the aperture with `GForvard`, so it carries the
   launch curvature. It enters the Shapiro reciprocity overlap
   (DOI 10.1364/JOSA.61.000492), `eta_turb = |sum(F_turb psi_tx)|^2 /
   |sum(F_vac psi_tx)|^2`.

## 3. Findings

### 3a. The sampling limit (DONE, documented)

The space sizer does not read the launch curvature. Past about
`lambda / (2 dx)` (113 urad on the 256 px, 6.83 mm test grid) `psi_tx` aliases
into a lattice of false bowl centres (spacing `lambda R / dx`). The vacuum
overlap then errs by -27 to +16 dB, with either sign. Below the limit the error
is under 0.28 dB. A sizer guard (`dx <= lambda / (4 theta)`) is NOT BUILT.
Figures: `figures/1_phase_cut.png`, `2_phase_map.png`, `3_overlap_error.png`.

### 3b. The wrong conjugate (FIXED, not committed)

The runner used `|sum(F conj(psi_tx))|^2`, the mode-match form of two fields
that travel the SAME way. Reciprocity puts the up-going launch against the
down-coming field and the Green's function is symmetric, so the correct form
has NO conjugate. The optimum launch is `psi_tx = conj(F)` (phase conjugation).

- Proof: direct propagation (launch, propagate, read the axis) against both
  forms. In vacuum with a lens f = +2L the truth is 3.76, no-conjugate 3.71,
  conjugate 0.45 (f = -2L is the mirror). Through one phase screen the
  no-conjugate overlap follows the truth at a constant ratio; the conjugate is
  5x to 35x off for a diverged launch.
- Effect: a curved `psi_tx` read the OPPOSITE curvature (a diverged launch was
  a converging beam that focused `R = w / theta` in front of the telescope). A
  REAL `psi_tx` (every collimated launch) is bit-identical, so every earlier
  study and campaign stands.
- The fix: `reciprocity_overlap(E, psi)` in `olb/waveoptics/turbulence/run.py`,
  used by every host overlap site (the in-run tail, the vacuum baseline, the
  point-ahead passes, `point_ahead_overlap`, `_PointAheadRunner`, and
  `campaign.py`); the device tail drops the conjugate too. The module
  self-check tests it against direct propagation.
- Checks that PASS: `python -m olb.waveoptics.turbulence.run` (about 50 s),
  `python -m olb.waveoptics.turbulence.campaign`, `python -m
  olb.models.waveoptics`.
- NOT changed: `validation/waveoptics_speed/profile_baseline.py` and
  `validation/uplink_sigma2i/uplink_farfield_reciprocity.py` still write the
  conjugate. Both use a collimated (real) launch, so their results stand; fix
  them if they are reused with a curved launch.

### 3c. A diverged launch is NOT grid-converged (OPEN, the blocker)

After the fix, the mean `eta_turb` moves with the pixel count at a fixed grid
side and plan (0.3 m aperture, 0.1 m waist, 600 km, 30 deg, rapid, 64 trials,
seed 7):

| divergence | 256 px | 512 px | 1024 px |
| --- | --- | --- | --- |
| none | 0.270 +/- 0.028 | 0.222 +/- 0.025 | - |
| 50 urad | 0.563 +/- 0.049 | 0.718 +/- 0.046 | 1.795 +/- 0.176 |
| 100 urad | 1.090 +/- 0.089 | 0.774 +/- 0.039 | 0.768 +/- 0.060 |

A mean above 1 is not physical: the far field of this launch falls
monotonically, so turbulent blurring can only LOWER the mean. The expected
long-term value (the Kolmogorov plane-wave coherence, r0 = 0.127 m) is 0.79 at
50 urad and 0.96 at 100 urad (`diagnostics/farfield_expectation.py`).

THE PROBABLE CAUSE: the vacuum baseline is not a plane wave. The slab starts
from a plane wave that fills the grid, and the absorbing edge makes Fresnel
rings over the 20 km slab (the Fresnel scale is about 0.18 m against a 1.75 m
grid). Inside the aperture `|F_vac|` spans 0.22 to 1.32 at 256 px, 0.78 to 1.30
at 512 px and 0.80 to 1.47 at 1024 px, with 0.5 to 0.7 rad of phase ripple
(`diagnostics/vacuum_flatness.py`). Rings are circular chirps, and so is a
curved `psi_tx`, so the overlap picks up the grid-dependent ring pattern. The
turbulent numerator scrambles the rings and the denominator keeps them. A
collimated `psi_tx` has no chirp and is much less sensitive (this agrees with
the earlier FAST parity of the collimated uplink). This is NOT PROVEN yet: no
fix has been tried.

## 4. The next plan (owner request): fidelity 1 (Dios) against fidelity 2

The owner wants to compare the POWER DISTRIBUTIONS of an uplink across the
launch divergence, a 15 cm ground aperture with NO central obscuration, and to
run the heavy part on the desktop (`ssh desktop`, bigfraw, the cupy GPU).

CAUTION given to the owner: a Dios comparison CANNOT prove reciprocity. Dios is
a different model (analytic, weak fluctuation, untruncated Gaussian, no outer
scale; `uplink_sigma2i/` found it over-reads sigma2_I 2x to 7x for a FILLED
launch). Only a direct upward propagation through the SAME screens tests
reciprocity. Arm A below does that.

- THE CONFIGURATION (proposed): D = 0.15 m, no obscuration, 1550 nm, HV5/7,
  a point receiver on the satellite. Waist 37.5 mm (D = 4 w: a negligible clip,
  so the field sees the same Gaussian as Dios; theta_min = 13.2 urad).
  Divergence: collimated, 2, 4, 8, 16 x theta_min (26, 53, 105, 210 urad).
  Elevations 30 and 60 deg. L0 = 25 m on the field.
- GATE 0 (FIRST, locally): make the field converge for a curved launch.
  Candidates: (a) a wider grid side so the rings leave the aperture, (b) a
  wide tapered (super-Gaussian flat-top) start field, (c) a normalisation on
  an ideal plane-wave baseline. Judge on the vacuum flatness in the aperture
  and on the mean `eta_turb` against the pixel count at 50 / 100 / 210 urad.
  Use `diagnostics/grid_convergence.py` and `diagnostics/vacuum_flatness.py` as
  the starting point.
- ARM A, the direct reciprocity check: rebuild the production loop by hand
  (assert `numpy.array_equal` against the public runner, as
  `validation/receiver_cone_clip/` does), propagate `psi_tx` UP through the
  SAME screens in reverse order plus the vacuum hop to the satellite, and
  compare the on-axis irradiance with the overlap, trial by trial. Keep it to
  <= 53 urad: the diverged beam at the slab top needs about 7000 px above that.
- ARM B, the fidelity-2 campaigns (bigfraw, cupy): one `Campaign` for each
  elevation WITH the field patch stored. The divergence enters through
  `psi_tx` only, so a validation helper rebuilds `psi_tx` for each divergence
  and re-reads eta post hoc: one set of runs, paired atmospheres across the
  divergence. 2000 trials for each elevation. The pixel must be
  `<= lambda / (4 theta_max)` = 1.8 mm.
- ARM C, fidelity 1: `uplink_turbulence_term` (the coupled-flux Monte Carlo),
  20k samples for each cell, locally. It takes the divergence through
  `launch_curvature` (the scintillation) and the `w_free` override (the
  broadening).
- THE METRICS, turbulence row only (the geometric row is analytic in both):
  the mean loss, sigma2_I, p1 / p5 / p50, CDF overlays, a two-sample KS
  statistic, and the trend against the divergence. Report, not fix: the outer
  scale (Dios is Kolmogorov, so its wander is larger on a collimated launch),
  the Dios off-axis sigma2_I, the weak-regime limit at 30 deg. An optional
  L0 = 1000 m field arm brackets the outer scale.
- THE OUTPUT: `validation/uplink_divergence_dios/` (scripts, a git-ignored
  `campaigns/`, figures, a README with the verdict).

OWNER DECISIONS STILL OPEN:
1. The waist: 37.5 mm (Dios-clean) or a filled launch (about 55 mm, more
   realistic, Dios known to be poor).
2. The orbit altitude (500 km matches the campaigns, 600 km the uplink
   self-checks) and the elevations.
3. Gate 0 first, locally, before bigfraw (recommended: yes).

## 5. Practical notes

- Python: `C:\Users\alexf\anaconda3\envs\olb\python.exe` (not on PATH). Run
  every script from the repository root as a module.
- bigfraw: `ssh desktop`, PowerShell 5.1 (chain with `;`, never `&&`), launch
  long runs through WMI `Win32_Process Create`, call
  `boost_process_priority()` in the parent AND each pool worker, the cupy venv
  `olb-gpu-venv`, pull data back as ONE zip. The desktop repo drifts; see the
  memory notes before a switch.
- `run_fidelity2` refuses `grid=`; use `propagate_turbulent_scenario(...,
  grid=g, plan=plan)` (both together) for a convergence study.
- Never force-add generated `.log` / `.json` / `.npz` files.
