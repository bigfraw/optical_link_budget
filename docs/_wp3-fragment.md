# WP3 fragment — sampling constraints (Schmidt Ch. 6.3-6.4, Ch. 7, Ch. 8)

A later step merges this fragment into `docs/schmidt-crosscheck.md`. Do not
edit the tracker directly while the parallel work packages run.

Book: J. D. Schmidt, "Numerical Simulation of Optical Wave Propagation with
Examples in MATLAB", SPIE Press (2010). DOI 10.1117/3.866274.

Module: `olb/waveoptics/schmidt/sampling.py`.

## Table 1 — existing olb code that implements a Ch. 6 to Ch. 8 relation

| olb code (file:line) | What it computes | Book equation | Agreement |
|---|---|---|---|
| `olb/waveoptics/grid.py:209-223` `forvard_max_z` | `z_max = N dx^2 / lambda` | Ch. 7, Eq. (7.59), printed p. 127, inverted with delta2 = delta1 | EXACT. But the docstring cites "Ch. 6". The rule is constraint 4 of Ch. 7. `schmidt/sampling.py` `angular_spectrum_max_z` reproduces it, and the self-check proves the two agree. |
| `olb/waveoptics/grid.py:105-110` FLAT EXTENT rule, `size = guard * 2 * r_max` | The grid side of a vacuum propagation | No book equation. Ch. 7, Eq. (7.18), printed p. 120, gives `D_grid >= (D_illum + D2)/2` instead | DIFFERENT RULE. The book sizes the grid from the illuminated area and the region of interest, and it ALLOWS the wrapped light outside D2. olb puts a fixed margin around the beam. See S-003. |
| `olb/waveoptics/grid.py:118-124` and `:166-172` RESOLUTION rule, `dx <= feature/(P/2)` | The pixel pitch | No book equation. The book picks the spacing per example: Ch. 7, Listing 7.1, printed p. 124 ("at least 50 grid pts across ap"); Ch. 8, text, printed p. 144 ("at least 30 grid points across D1 and D2") | COMPATIBLE IN FORM, COARSER IN VALUE. See Table 3. |
| `olb/waveoptics/grid.py:169` `n_wanted = 2 ** ceil(log2(size/dx))` | The pixel count | Ch. 7, Listing 7.1 line 11, printed p. 124, and Listing 7.2 line 13, printed p. 128 | EXACT. The book rounds N up to the next power of two for the FFT. |
| `olb/waveoptics/grid.py:179-185` the range warning | Warns when `z > forvard_max_z` | Ch. 7, Eq. (7.59), printed p. 127 | EXACT rule, wrong citation (Ch. 6). |
| `olb/waveoptics/turbulence/sampling.py:370-380` and `:441-443` the extent rule | `side = [guard*2*r_beam + 2*(lambda/r0)*z] / (1 - b)` | Ch. 9 (the scattering cone). Ch. 8, Sec. 8.1, printed p. 134, gives the absorbing band only | OUT OF SCOPE for WP3, except the `(1 - b)` divisor. See Table 3. |
| `olb/waveoptics/turbulence/sampling.py:382-386` and `:451-455` the pixel rule | `dx <= min(r0/P_r0, sqrt(lambda z_i)/2, feature/(P/2))` | Ch. 9 (r0), Andrews Ch. 8 (the Fresnel scale). Ch. 7 gives none of the three | NO Ch. 7 CONTENT. The turbulent sizer never evaluates constraints 1 to 4. See S-002. |
| `olb/waveoptics/turbulence/sampling.py:457-458` `n = clamp(2**ceil(log2(side/dx)), 256, n_max)` | The pixel count | Ch. 7, Listing 7.2 line 13, printed p. 128 (the power-of-two step only) | The power-of-two step is EXACT. The `[256, n_max]` clamp has no book source. See Table 3. |
| `olb/waveoptics/turbulence/sampling.py:469` `step_over_limit_max = max(gap)/forvard_max_z` | Reports the worst split-step length against constraint 4 | Ch. 8, Eq. (8.24), printed p. 144 | SAME IDEA, DIFFERENT ROUTE. The book SETS the step count from Eq. (8.24). olb sets it from the Cn2 profile and only REPORTS the ratio. See S-006. |
| `olb/waveoptics/turbulence/splitstep.py:29` `super_gaussian_boundary(n, width_frac=0.125, power=8)` | The absorbing boundary | Ch. 8, Eq. (8.1), printed p. 134; Listing 8.1, printed p. 142; Fig. 8.1, printed p. 134 | FORM AGREES (a super-Gaussian, exponent > 2). The NUMBERS are a different parameterisation. See Table 3 and S-007. |

## Table 2 — gaps

| Id | Gap | Book source | Note |
|---|---|---|---|
| S-001 | No TWO-STEP Fresnel propagator. `olb/waveoptics/propagators.py` has `Forvard`, `Fresnel` and `GForvard`, and `lenses.py` has the co-moving route. None of them frees the magnification with a second Fresnel integral. | Ch. 6, Sec. 6.3.2, Eqs. (6.24) to (6.29), printed pp. 94 and 95 | `schmidt/sampling.py` `two_step_planes` gives the two intermediate-plane geometries. A propagator that uses them is not built. |
| S-002 | No budget, sizer, or runner calls a sampling checker. `GridSpec.for_scenario` and `turbulent_grid` warn on their OWN rules only. | Ch. 7, Sec. 7.3.3, printed p. 127 | Wire `check_sampling` into `GridSpec.for_scenario` and into `SamplingReport`. It never raises, so it cannot break an existing run. |
| S-003 | Constraints 1 and 2 are implemented NOWHERE in olb. The observation-region extent D2 never enters a grid decision. | Ch. 7, Eqs. (7.14) and (7.20), printed pp. 119 and 120 | olb uses the receive aperture as a FEATURE (a pixel rule), never as D2 (an extent rule). |
| S-004 | Constraint 3 is not checked. The transmit-beam curvature R never reaches a grid rule. | Ch. 7, Eq. (7.53), printed p. 126 | This is the same missing curvature thread as Gap 3 of the Andrews cross-check. |
| S-005 | The Fresnel-integral MINIMUM distance is not checked. `propagators.py:136` `Fresnel` has no near-distance guard, so a short call aliases silently. | Ch. 7, Eqs. (7.41) and (7.42), printed p. 123 | `fresnel_min_distance` gives the bound. |
| S-006 | No vacuum partial-propagation planner. A long vacuum link takes the co-moving lens route instead of a chain of angular-spectrum steps. | Ch. 8, Eqs. (8.23) and (8.24), printed pp. 143 and 144 | `partial_max_step` and `partial_plane_count` give the count. The turbulent planner sets its count from Cn2, not from Eq. (8.24). |
| S-007 | The absorbing-boundary width carries no Schmidt number. | Ch. 8, Listing 8.1, printed p. 142; Fig. 8.1, printed p. 134 | `absorbing_boundary_sigma` gives the book value. See Table 3. |
| S-008 | `forvard_max_z` cites "Ch. 6". The rule is Ch. 7, Eq. (7.59), printed p. 127. | Ch. 7, Eq. (7.59), printed p. 127 | A one-line docstring fix. |

## Table 3 — book values for the seeded constants

| olb constant (file:line) | olb value | Book value | Verdict |
|---|---|---|---|
| `grid.py:98` `guard=4.0` | The grid half-side is 4 beam radii | The book gives NO guard factor. Ch. 7, Eq. (7.18), printed p. 120, sizes the grid as `(D_illum + D2)/2`, so the wrapped light gets exactly half way around and stops at the edge of D2 | UNCITED, and it is a DIFFERENT PHILOSOPHY. The book tolerates aliasing outside the region of interest; olb forbids it everywhere. The olb rule is stricter for a wide beam and it is silent about D2. |
| `grid.py:98` `pixels_per_feature=16` | 16 pixels across the smallest hard edge | Ch. 7, Listing 7.1, printed p. 124: 50 points across D1. Ch. 8, text, printed p. 144: 30 points across D1 and across D2 | UNCITED, and COARSER than both worked examples by a factor of 2 to 3. The book treats the number as a per-problem choice, not a constant. |
| `turbulence/sampling.py:56` `PIXELS_PER_FEATURE = 8` | 8 pixels across the smallest hard edge | The same two book choices, 50 and 30 | UNCITED, coarser again by a factor of 4 to 6. |
| `grid.py:36` `N_MIN = 256` | A floor on the pixel count | The book has NO floor. N comes from the constraints and then rounds to the next power of two. The two worked examples land at N = 128 (Ch. 7, printed p. 123) and N = 512 (Ch. 7, printed p. 128); the Ch. 8 example uses N = 128 (printed p. 144) | UNCITED. Note that all three book examples sit at or BELOW the olb floor, so the floor never binds on a book-sized problem. It is a convenience, not physics. |
| `grid.py:209` `forvard_max_z = N dx^2 / lambda` | The angular-spectrum range limit | Ch. 7, Eq. (7.59), printed p. 127, with delta2 = delta1 | EXACT. Only the citation is wrong (S-008). |
| `turbulence/splitstep.py:29` `width_frac=0.125`, `power=8` | The super-Gaussian absorber | Ch. 8, Listing 8.1, printed p. 142: `sigma = 0.47 N` grid points with the exponent 16 on r. Fig. 8.1, printed p. 134: `sigma = 0.45 L`, n = 16 | FORM AGREES (Eq. (8.1) needs n > 2). The parameterisations differ: the book scales sigma by the FULL side (0.47 N, so 0.94 of the half-side), olb by the half-side band width (0.125). The exponent differs too: 16 in the book against 8 in olb, so the olb taper is softer. |

## WP3 note

### Built

`olb/waveoptics/schmidt/sampling.py`, plus a one-line placeholder
`olb/waveoptics/schmidt/__init__.py` (the package directory did not exist).

Small pure functions, numpy only, no olb import:

- Ch. 7.1 and 7.2, the band limit and the geometry: `nyquist_max_angle`
  (Eq. (7.7)), `geometric_max_angle` (Eqs. (7.8), (7.9), (7.12)),
  `illuminated_diameter` (Eq. (7.16)).
- The four numbered constraints: `constraint1_max_delta2` (Eq. (7.14)),
  `constraint2_min_n` (Eq. (7.20)), `constraint3_delta2_window` (Eq. (7.53)),
  `constraint4_min_n` (Eq. (7.59)), and `constraint3_is_slack` (Eq. (7.60)).
- The local-frequency analysis the constraints come from:
  `local_spatial_frequency_source` (Eqs. (7.37), (7.39), (7.51)) and
  `local_spatial_frequency_transfer` (Eqs. (7.55), (7.57)).
- Per-kernel: `one_step_delta2` (Eq. (7.21)), `one_step_min_n` (Eq. (7.25)),
  `fresnel_min_distance` (Eqs. (7.41), (7.42)), `two_step_planes`
  (Eqs. (6.24) to (6.29)), `angular_spectrum_max_z` (Eq. (7.59) inverted).
- Ch. 8: `partial_grid_spacing` (Table 8.2), `partial_max_step` (Eq. (8.24)),
  `partial_plane_count` (text below Eq. (8.24)), `absorbing_boundary_sigma`
  (Listing 8.1 and Fig. 8.1).
- `check_sampling` returns five `Rule` tuples
  (name, satisfied, bound, actual, citation). It NEVER raises and it never
  warns. The caller decides.

The self-check reproduces the book's own worked numbers:

| Book place | Quantity | Book | Module |
|---|---|---|---|
| Ch. 7, Eq. (7.43), printed p. 123 | N_min, one step | 66 | 65.79 |
| Ch. 7, Eq. (7.43), printed p. 123 | N used | 128 | 128 |
| Ch. 7, Eq. (7.43), printed p. 123 | delta2 | 97.7 um | 97.66 um |
| Ch. 7, Eq. (7.42), printed p. 123 | z_min | 8 cm | 8.00 cm |
| Ch. 7, Sec. 7.3.2, printed p. 127 | log2 N, constraint 4 | 8.55 | 8.55 |
| Ch. 7, Sec. 7.3.2, printed p. 128 | log2 N, constraint 2 | 8.51 | 8.51 |
| Ch. 7, Sec. 7.3.2, printed p. 128 | N used | 512 | 512 |
| Ch. 8, Eq. (8.25), printed p. 144 | Delta_z max | 0.567 m | 0.569 m |
| Ch. 8, text, printed p. 144 | planes n | 5 | 5 |
| Ch. 6, Table 6.2, printed p. 95 | two-step planes, m = 2, 1, 1/2 | 1/3, 2/3, -1, 2 and 1/2, 1/2, inf and 2/3, 1/3, 2, -1 | identical |

The self-check also proves the DERIVATIONS close, not only the numbers: the
Nyquist rule on Eq. (7.39) at the source edge gives back Eq. (7.42); the same
rule on Eq. (7.51) gives back the constraint-3 upper bound; the same rule on
Eq. (7.57) gives back constraint 4; Eq. (7.31) reproduces Eq. (7.25) exactly;
and Eq. (7.18) on `illuminated_diameter` reproduces constraint 2 exactly.

Run: `python -m olb.waveoptics.schmidt.sampling`.

### Decisions

- ONE checker, a list of five plain namedtuples. No severity levels, no
  fixer, no auto-sizer. The module measures; the caller acts.
- `check_sampling` returns ALL five rows for every call, and each row carries
  its citation. It does not take a `method` argument, because the citation
  already tells the caller which kernel a row governs (rows 1 and 2 are
  geometry and hold for all three kernels; rows 3 and 4 are the
  angular-spectrum kernel; row 5 is the two Fresnel-integral kernels).
- The per-kernel assumption sets live in the docstrings of `one_step_delta2`,
  `two_step_planes` and `angular_spectrum_max_z`, one kernel per docstring.
  Each names the Fresnel-approximation validity, what the kernel fixes or
  frees about the grid spacing, when it aliases, and which constraint governs
  it.
- `local_spatial_frequency_source` takes ONE optional `m`. `m=None` gives the
  Fresnel-integral phase curvature 1/z + 1/R (Eq. (7.39)); a value of `m`
  gives the angular-spectrum curvature (1 - m)/z + 1/R (Eq. (7.51)). Two
  equations, one function, because only the curvature differs.
- No MATLAB listing is ported. Every function is written from the printed
  equations.
- The module does NOT change any existing sizer. Wiring it is S-002, an owner
  decision, because a wired check moves no numbers but it will print warnings
  on grids that run today.

### The book would not give

- NO guard factor and NO margin philosophy that matches `grid.py`. The book
  lets the wrapped light come half way around the grid, up to the edge of D2
  (Eq. (7.18), printed p. 120). It never asks for empty space around the beam.
  So `guard=4.0` cannot be justified from Ch. 7; only replaced by constraint 2.
- NO fixed pixels-per-feature number and NO N floor. The book picks 50 points
  across the aperture in one example and 30 in another, and it calls the whole
  analysis "a guideline ... not unbreakable rules" (Ch. 7, Sec. 7.3.3, printed
  p. 129).
- NO obscured or annular aperture rule. D1 and D2 are plain extents.
- NO turbulence in Ch. 7 or Ch. 8. The screen rules are Ch. 9, so the
  `min_screens` and `pixels_per_r0` question stays open after this work package.
- TWO arithmetic slips in the book's own worked numbers, both reproduced above
  and both harmless:
  1. Eq. (8.25), printed p. 144, prints `(66.7 um)^2 * 128 / 1 um = 0.567 m`.
     The arithmetic gives 0.569 m. The plane count (5) is the same either way.
  2. The Ch. 8 example, printed p. 144, reads "at least N = 2^7 = 128 grid
     points are required" off the Fig. 8.5 contour plot. Constraint 2 with
     delta1 = 66.7 um and delta_n = 133 um gives N >= 142.8, which is 2^7.16.
     (The book prints D1 = 2 mm and the 30-point choice. It does not print D2;
     D2 = 4 mm follows from delta_n = 133 um times 30 points.)
     The book then uses N = 128 for Eq. (8.25). So the printed example
     VIOLATES its own constraint 2 by 11 percent.
- NO two-sided form of Eq. (7.41). The printed bound is one-sided. A
  converging source with `lambda R < D1 delta1` has no valid distance at all;
  the module returns `math.inf` there, which is an olb decision, not a book
  statement.
