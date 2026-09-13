# The screen plan against the point-ahead anisoplanatism

Date 2026-09-11, branch `waveoptics-pointahead`, backlog 2-P4. Both studies
run on a laptop CPU. No campaign, no propagation.

## The purpose

A pre-compensated uplink senses the downlink beacon and launches ahead of it
by the point-ahead angle. In the plane-parallel screen model the two paths
cross screen j at a lateral offset `theta * z_g,j`, where `z_g,j` is the
ground distance of the screen. The part of the sensed phase that the two
offsets do not share is the anisoplanatic residual, and Stone et al. (1994)
give its variance as a continuous integral over the Cn2 profile.

The fidelity-2 runner replaces that integral with a sum over its screens. Two
questions follow, and the owner suspected that the answer to the first one is
"more screens near the ground, where the offset is small":

1. Does the production screen plan (equal Rytov weight, 9 screens at the
   standard preset) reproduce the continuous Stone variance, and does a plan
   cut by the anisoplanatic weight do better?
2. Does a screen drawn once and read through two integer-pixel windows give
   the Stone variance when the two windows are summed and differenced?

THE CASE. The hero of `validation/waveoptics_ao/` turned around: a 1550 nm
uplink from a 700 mm ground SMF terminal to a 100 mm space terminal at 500 km,
HV5/7 site, `L0 = 25 m` on the screens where the study says so, the standard
preset grid (512 px at every elevation). The ground transmitter fills the
aperture (waist 0.35 m). `common.py` holds the case and the plan builders.

## THE CONVENTION: the tilt stays in

The mode set of record is `remove = "piston"` (owner decision, 2026-09-11).
The terminal senses the DOWNLINK beacon tilt, and the steering mirror adds the
point-ahead offset geometrically. So the terminal holds no uplink tilt
reference, and the uplink pays the FULL tilt anisoplanatism. The old
convention `remove = "piston_tilt"` assumed a separate uplink tilt loop, and
that loop does not exist in this design.

The tilt is the biggest single part of the error. At 30 deg, D = 0.7 m and
10 arcsec the continuous Stone variance goes from 0.703 rad^2 with the tilt
out to 1.708 rad^2 with the tilt in at `L0 = inf` (a factor 2.43), and from
0.702 to 1.680 rad^2 at `L0 = 25 m` (a factor 2.39).

The library function `anisoplanatic_phase_variance` keeps its general
`remove` argument, and the production Term
`olb.links.uplink.uplink_point_ahead_term` now defaults to `remove="piston"`.
Study 1 takes `--remove` (default `piston`), so the old rows come back with
`--remove piston_tilt`.

## The outer scale in Stone

Stone writes the Kolmogorov spectrum. The screens of record carry a finite
outer scale, so the two must agree on it. The substitution is one line. The
kernel `u^(-8/3)` of Eq. (36) comes from the polar element `kappa d kappa`
times the spectrum `kappa^(-11/3)`, with `u = kappa R`:

    Kolmogorov:  u^(-8/3) = u (u^2)^(-11/6)
    von Karman:  u (u^2 + u0^2)^(-11/6),   u0 = kappa0 R = 2 pi R / L0

with `Phi(kappa) = 0.033 Cn2 (kappa^2 + kappa0^2)^(-11/6)` and
`kappa0 = 2 pi / L0` (Andrews and Phillips, 2nd ed. (2005),
DOI 10.1117/3.626196, Ch. 3, Eq. (20), printed p. 68). Every prefactor stays
the same, and `L0 = inf` gives Stone back exactly.
`olb.turbulence.anisoplanatism._inner_integral` takes `u0`, and
`anisoplanatic_phase_variance` takes `L0`. The isoplanatic angle stays
Kolmogorov only: a finite outer scale breaks the `(theta/theta0)^(5/3)` power
law that defines `theta0`.

THE EFFECT IS SMALL, and that is the physics. The anisoplanatic variance is a
DIFFERENCE of two wavefronts, so it already cancels every scale much larger
than the offset `theta z_g`. At 30 deg, D = 0.7 m and 10 arcsec, `L0 = 25 m`
takes 1.6 percent off the piston-removed variance and 0.1 percent off the
piston-and-tilt-removed one. The tilt holds the large scales, so it is the
mode the outer scale touches most, but the whole effect stays under 2 percent
in the owner regime.

Study 2 at `L0 = 25 m` is the NUMERICAL VALIDATION of that kernel: the screens
carry the 25 m outer scale and the Stone reference carries it too, so the
ratio is a pass test at both outer scales.

## The run lines

    C:\Users\alexf\anaconda3\envs\olb\python.exe -m validation.anisoplanatism_screens.anisoplanatism_screens
    C:\Users\alexf\anaconda3\envs\olb\python.exe -m validation.anisoplanatism_screens.anisoplanatism_screens --L0 25
    C:\Users\alexf\anaconda3\envs\olb\python.exe -m validation.anisoplanatism_screens.phase_only_shift
    C:\Users\alexf\anaconda3\envs\olb\python.exe -m validation.anisoplanatism_screens.phase_only_shift --L0 25

The first script takes 1 to 5 minutes (the finite outer scale needs one Stone
table for each aperture, so `--L0 25` is the slower one). The second takes
about 15 minutes at the default 100 draws; `--n-draws 10` is the smoke run.
Each script writes its log, its results JSON and its figure next to itself,
and study 1 also writes `production_table[_L0<m>].md`.

## Study 1: the discrete Stone sum of a plan (`anisoplanatism_screens.py`)

The script evaluates the delta-layer form of Stone Eq. (29) on a screen plan:

    sigma^2 = 2 (2 pi)^(8/3) C_A k0^2 R^(5/3) SUM_j m_j I(z_g,j theta / R)

where `m_j` is the slant Cn2 integral of screen j (`ScreenPlan.cn2_int_m13`,
which already carries the airmass factor), `R = D/2`, and `I(beta)` is the
Stone inner integral of `olb/turbulence/anisoplanatism.py` with the piston
removed and the tilt kept. The reference is the continuous integral of the
same module on the fine height grid of the planner (ratio 0.9997 against the
library function, at both outer scales).

The sweep: elevations 20, 30, 60, 90 deg; apertures 0.4, 0.7, 1.0 m; angles
2, 5, 10 arcsec plus the geometry angle; `max_order` 3 (AO of 10 Noll modes,
orders 1 to 3 with the tilt in), 5 (21 modes, orders 1 to 5) and None (every
order). Plans: the production plan, the equal-Rytov override at 5, 9, 15, 25
screens, the ground split (the lowest screen cut in 4), and two
equal-anisoplanatic-weight families (the slab cut at equal shares of
`Cn2 sec (h sec)^2`, the screen at the Cn2 centroid or at the matched
centroid). The small-beta exponent of `I` is still 2.000 with the tilt kept
(the piston-removed weight `M(u)` still goes as `u^2` at a small `u`), so the
plan family does not change.

The production plan at D = 0.7 m, `max_order` None, `L0 = inf`:

| elevation | angle [arcsec] | continuous [rad^2] | discrete [rad^2] | ratio |
| --- | --- | --- | --- | --- |
| 20 | 10.0 | 3.167 | 3.146 | 0.994 |
| 30 | 10.0 | 1.708 | 1.697 | 0.994 |
| 60 | 10.0 | 0.6695 | 0.6653 | 0.994 |
| 90 | 10.0 | 0.5174 | 0.5142 | 0.994 |
| 30 | 5.24 (geometry) | 1.075 | 1.068 | 0.994 |

The same cells at `L0 = 25 m` (the variance falls by 1 to 2 percent, the ratio
does not move):

| elevation | angle [arcsec] | continuous [rad^2] | discrete [rad^2] | ratio |
| --- | --- | --- | --- | --- |
| 20 | 10.0 | 3.097 | 3.077 | 0.994 |
| 30 | 10.0 | 1.680 | 1.669 | 0.994 |
| 60 | 10.0 | 0.6627 | 0.6585 | 0.994 |
| 90 | 10.0 | 0.5128 | 0.5096 | 0.994 |
| 30 | 5.24 (geometry) | 1.065 | 1.058 | 0.994 |

The per-screen shares at 30 deg, 10 arcsec, D = 0.7 m (`L0 = inf`) are much
FLATTER than they were with the tilt out, because the tilt is a large-scale
mode and the high screens carry the large offsets: the top screen (16.2 km,
offset 1.57 m) holds 9.7 percent, the screen at 4.1 km holds the most at
16.2 percent, and the two lowest screens (0.57 km and 0.08 km, offsets 8.0 px
and 1.1 px) hold 10.3 percent together. With the tilt out the same two lowest
screens held 21.4 percent.

**VERDICT: NO PLANNER CHANGE.** Over the 135 production cells of the owner
regime (D at most 1 m, angle at most 10 arcsec, elevation at least 20 deg) the
9-screen ratio holds between 0.987 and 0.995 at `L0 = inf`, and between 0.987
and 0.995 at `L0 = 25 m`. The worst cell of both is 20 deg, D = 0.4 m,
10 arcsec, `max_order` 3 at 0.987. The equal-anisoplanatic-weight families do
NOT improve it: a slab that is thick in height breaks the small-angle limit
that the centroid assumes, and the worst owner-regime ratio of that family is
3.25 at 5 screens. The integer-pixel rounding of the offsets moves the
production ratio by at most +0.035 (20 deg, D = 0.4 m, 5 arcsec). The
suspicion that the ground layers need more screens is not supported: the
lowest screen holds 2.2 percent of the variance at an offset of one pixel, and
the ground split moves nothing.

The full 144-row table is `production_table.md` (`L0 = inf`) and
`production_table_L025.md` (`L0 = 25 m`). The script writes it on every run.

![discrete against continuous](figures/discrete_vs_continuous.png)

## Study 2: the shifted screen windows (`phase_only_shift.py`)

The script draws each screen ONE time on an oversize grid with the production
`ScreenFactory` (the same L0, subharmonics and precision as the runner), reads
the beacon window `scr[0:n, 0:n]` and the uplink window `scr[0:n, s_j:s_j+n]`
with `s_j = round(theta z_g,j / dx)`, sums the difference over the plan, and
fits 36 Noll modes over 16 disjoint 0.7 m apertures (25 of 0.4 m) for each of
100 draws. Three quantities, all with the TILT IN: the piston-removed
variance, and the Noll bands 2 to 10 and 2 to 21 (the `max_order` 3 and 5 of
Stone). The reference is the delta-layer sum of study 1 at the SAME rounded
offsets and at the SAME outer scale, and the continuous integral as a second
column.

At `L0 = inf`, D = 0.7 m, the production count (9 screens):

| angle [arcsec] | quantity | measured [rad^2] | Stone [rad^2] | ratio +- 2 sigma |
| --- | --- | --- | --- | --- |
| 2 | piston removed | 0.415 | 0.421 | 0.987 +- 0.025 |
| 5 | piston removed | 1.058 | 1.046 | 1.012 +- 0.038 |
| 10 | piston removed | 1.716 | 1.688 | 1.017 +- 0.041 |
| 10 | band 2 to 10 | 1.419 | 1.385 | 1.024 +- 0.050 |
| 10 | band 2 to 21 | 1.530 | 1.498 | 1.021 +- 0.047 |

The same cells at `L0 = 25 m`:

| angle [arcsec] | quantity | measured [rad^2] | Stone [rad^2] | ratio +- 2 sigma |
| --- | --- | --- | --- | --- |
| 2 | piston removed | 0.413 | 0.419 | 0.987 +- 0.025 |
| 5 | piston removed | 1.048 | 1.036 | 1.012 +- 0.038 |
| 10 | piston removed | 1.687 | 1.660 | 1.016 +- 0.040 |
| 10 | band 2 to 10 | 1.390 | 1.358 | 1.024 +- 0.050 |
| 10 | band 2 to 21 | 1.501 | 1.470 | 1.021 +- 0.046 |

**VERDICT: the window rule holds at BOTH outer scales.** 52 of 54 cells (5, 9,
15 screens, three angles, two apertures, three quantities) sit inside the 0.95
to 1.05 band or inside their own 2-sigma bar, at `L0 = inf` and at
`L0 = 25 m`. The two cells outside are the SAME two in both runs, and they
miss by a hair: 15 screens, 10 arcsec, D = 0.7 m, band 2 to 10 at
1.055 +- 0.054, and band 2 to 21 at 1.051 +- 0.050. The ratio holds between
0.959 and 1.055 over every cell.

**THE OUTER-SCALE KERNEL IS VALIDATED.** Between the two runs the MEASURED
variance of the screens falls by 1.7 percent (1.716 to 1.687 rad^2 on the
headline cell) and the von Karman Stone reference falls by the SAME 1.7
percent (1.688 to 1.660 rad^2), so the ratio moves by 0.001. The kernel of
item 1 tracks the screens. The grid holds the 25 m scale with room to spare:
no cell raised a subharmonic-reach warning (768 px at 6.86 mm is a 5.3 m side,
and the three subharmonic levels reach 3^3 times that).

The sub-pixel arm shows why the rounding matters at the bottom: with the two
lowest screens forced to 0 px the piston-removed variance at 10 arcsec falls
from 1.716 to 1.533 rad^2 at `L0 = inf` (the true offsets are 8.0 and 1.1 px),
so the runner must keep the rounded offset of every screen, and it does.

![shifted windows against Stone](figures/phase_only_shift.png)

## The caveats

- The tilt is a LARGE-SCALE quantity, so it carries much more sample spread
  than the tilt-removed residual of the old convention. The 2-sigma bar of a
  100-draw cell is therefore wider than it was.
- Both studies are phase only. The propagated check (diffraction between the
  screens, the reciprocity overlap, the AO correction) is the campaign study
  `validation/waveoptics_pointahead/`.
- The plane-parallel offset `theta z_g` ignores the Earth curvature. At
  20 deg and 58 km of slant path the error is under 1 percent of the offset.
- The screens carry a finite subharmonic reach. `ScreenFactory` raises the
  level count when `L0` is larger than the reach of the grid
  (`screens.required_n_sub_levels`), and the log of each cell prints any
  warning it gets.

## The files

| File | Purpose |
| --- | --- |
| `common.py` | The hero uplink case, the plan builders (production, override, ground split, equal anisoplanatic weight) and the log helpers. |
| `anisoplanatism_screens.py` | Study 1. Takes `--remove` and `--L0`. Writes `anisoplanatism_screens[_L0<m>].log`, `_results.json`, `production_table[_L0<m>].md`, `figures/discrete_vs_continuous[_L0<m>].png`. |
| `phase_only_shift.py` | Study 2. Takes `--L0`. Writes `phase_only_shift[_L0<m>].log`, `_results.json`, `figures/phase_only_shift[_L0<m>].png`. |

## Sources

- J. Stone, P. H. Hu, S. P. Mills and S. Ma, J. Opt. Soc. Am. A 11(1),
  347-357 (1994), DOI 10.1364/JOSAA.11.000347. Eqs. (27), (29), (36), (43).
- R. J. Noll, J. Opt. Soc. Am. 66, 207 (1976), DOI 10.1364/JOSA.66.000207.
- L. C. Andrews and R. L. Phillips, 2nd ed. (2005), DOI 10.1117/3.626196.
  Ch. 12, Eq. (14), the plane-parallel airmass. Ch. 3, Eq. (20), printed
  p. 68, the von Karman spectrum of the outer-scale kernel.
