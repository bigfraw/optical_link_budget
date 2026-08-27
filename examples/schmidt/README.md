# examples/schmidt — the Schmidt numerical foundation layer

Three runnable scripts. Each one prints a labelled table of numbers that the
package computed, and saves its figures next to the script. No script invents a
number, and no script changes an `olb` module.

Every equation comes from one book:

> J. D. Schmidt, *Numerical Simulation of Optical Wave Propagation with
> Examples in MATLAB*, SPIE Press Monograph PM199 (2010).
> DOI: 10.1117/3.866274

Each module docstring names the chapter, the equation number and the printed
page of the method that it shows. The package itself is
`olb/waveoptics/schmidt/`. Its documentation is `docs/physics.md` Section 8 and
`docs/api-waveoptics.md` Section 10. The equation-by-equation audit is
`docs/schmidt-crosscheck.md`.

The Andrews and Phillips suite in `examples/andrews/` is the twin of this one.
Andrews owns the ANALYTIC value of a quantity. Schmidt owns the SIMULATION
rule: the transforms, the propagators, the sampling constraints, the absorbing
boundary, and the phase screens.

## How to run

From the repository root:

    python -m examples.schmidt.propagator_kernels

The `-m` form is the house idiom of `examples/`. It puts the repository root on
the import path, so `import olb` works. Run each script the same way, with its
own module name. The screen script needs `aotools`, because it reads the
production generator: `pip install aotools`, or `pip install olb[screens]`.

## The scripts

| Script | What it shows |
| --- | --- |
| `propagator_kernels.py` | The book kernels against the production LightPipes propagators, in three tiers. It bridges the two convention gaps first: the piston phase `exp(i k z)` that the production propagators keep and the book listings drop, and the quadrature. Figures: `propagator_kernels_cuts.png`, `propagator_kernels_diffs.png`. |
| `sampling_and_edges.py` | Part A is a gallery of deliberate sampling failures, each one paired with the grid that obeys the rule: the transfer-function range limit, the transfer phase that does the folding, and the edge wrap of a partial-propagation chain with and without the absorber. Part B runs `check_sampling` and `properly_sampled_checklist` on the REAL production grids, with a citation on every row. Figures: `sampling_artefacts.png`, `sampling_artefacts_edges.png`, `sampling_absorbers.png`. |
| `screens_and_turbulence.py` | One screen by three generators on one colour scale; the ensemble structure function of the three against Eq. (9.44); and the per-screen strength rule with the factor-4 bridge between the two variance conventions. Figures: `screens_examples.png`, `screens_structure.png`, `screens_strength.png`. |

## The measured agreements

The numbers below print when the scripts run. They are stable across seeds.

**The propagators.**

- Tier (a), the SAME algorithm: `schmidt.fresnel.angular_spectrum` at m = 1
  against the production `Forvard`. The two agree to about 1e-10 over the whole
  grid.
- Tier (b), the same integral by another quadrature: the one-step and two-step
  Fresnel kernels against the production `Fresnel`. In the interior of the grid
  the agreement is 6e-4 for a soft Gaussian and 1.5e-2 for a hard truncation.
  The production `Fresnel` convolves on a doubled grid and differences four
  shifted copies, so its artefact sits in the outer band, where the field is
  below 1e-8 of the peak.
- Tier (c), the co-moving grid: `two_step_fresnel` with a free output pitch
  against the production `Lens -> LensFresnel -> Convert` recipe, at a
  magnification of 247. The agreement is 1.7e-3 soft and 2.3e-2 hard.

**The screens.**

- The book subharmonic generator reaches 0.88 to 0.93 of the theory of
  Eq. (9.44) over `r/r0 = 0.3` to 1.6.
- The `aotools` generator of the production layer reads 1 to 3 percent ABOVE
  the book generator there. The two agree well inside the band.
- The plain Fourier screen reaches 0.69 to 0.82 in the band, and it falls to
  0.47 at `r/r0 = 8`. Both subharmonic generators fall away past the band too,
  to 0.80 and 0.86. The deficit is real, and no subharmonic level removes it.
- The bridge between the two per-screen conventions measures 3.9994. The book
  caps the LOG-AMPLITUDE variance `sigma_chi^2` at `rmax = 0.1` (Listing 9.5,
  printed p. 175). The production planner caps the plane-wave RYTOV variance
  `sigma_R^2`, and `sigma_R^2 = 4 sigma_chi^2`. So the book cap is 0.4 on the
  production number, and the presets 0.05 / 0.10 / 0.25 are 8x / 4x / 1.6x
  stricter than the book.

## Wiring status

### LIVE: nothing

The Schmidt layer is VALIDATION ONLY, by owner decision. No budget, no Term, no
sizer and no runner reads it. `olb/waveoptics/schmidt/__init__.py` exports
nothing. The production wave-optics code stays the LightPipes port, and its
bodies do not change.

Work package 6 wrote the book citations INTO the production modules. That
retrofit changed docstrings and comments only, and it moved no number. It
corrected one mis-citation: `olb/waveoptics/grid.py` `forvard_max_z` cited
"Ch. 6", and the rule is constraint 4, Ch. 7, Eq. (7.59), printed p. 127, at
m = 1 (also Ch. 8, Eq. (8.24), printed p. 144).

### AVAILABLE: the checker functions

Any future work can call these. Each one measures and returns; none warns, and
none raises on a broken rule.

- `schmidt.sampling.check_sampling` — the five vacuum rows of Ch. 7, each with
  its citation. Wiring it into `GridSpec.for_scenario` is gap S-06.
- `schmidt.turbulence.properly_sampled_checklist` — the Sec. 9.5 rows: the same
  constraints with the turbulence-blurred extents of Eqs. (9.84) and (9.85),
  plus the two pitch rules of Sec. 9.4.
- `schmidt.turbulence.moment_error`, `profile_moments` and `layer_moments` —
  the layer moment rule of Eq. (9.65). It is the principled replacement for the
  `_merge_layers` bail-out, and it gives a screen-count floor of 4.
- `schmidt.turbulence.screen_strengths` and `max_screen_strength` — the
  constrained solve for the screen strengths under the `rmax` cap.
- `schmidt.sampling.fresnel_min_distance` — the minimum distance that the
  production `Fresnel` states in words and does not check (gap S-09).
- `schmidt.fresnel.partial_propagations` — the book chain with a per-plane
  pitch. The production split step holds one flat pitch (gap S-14).

### The notable conflicts

Both are recorded, not changed.

- **The absorber shape.** The production `super_gaussian_boundary` is RADIAL,
  of power 8, exactly 1.0 inside 0.875 of the half-side, and `exp(-1)` at the
  middle of an edge. The book runs a SEPARABLE mask of power 16 with one
  half-width of `0.47 N` pixels (Listing 8.1, printed p. 142; Listing 9.7,
  printed p. 179), which is 0.068 at the middle of an edge. So the book absorbs
  about 5 times harder there. Eq. (8.1), printed p. 134, allows any power above
  2 and gives no half-width, so the production power of 8 is inside the
  equation. The book itself records that Flatte and others used 8. Gap S-15.
- **The split-step flat pitch.** The production `split_step` calls `Forvard` at
  m = 1 on ONE flat grid for every step. Ch. 8, Eq. (8.18), printed p. 139,
  gives each partial propagation its own pitch, from the linear rule of
  Eq. (8.8), printed p. 136. So the production grid cannot grow with a
  diverging beam: the beam must fit the SOURCE grid at the receiver. Gap S-14.

### The book errata found

Three printed numbers do not reproduce. All three are recorded in
`docs/schmidt-crosscheck.md`, and none of them changes an `olb` value.

- **Eq. (8.25), printed p. 144.** The book prints
  `(66.7 um)^2 * 128 / 1 um = 0.567 m`. The arithmetic gives 0.569 m. The plane
  count of 5 is the same either way.
- **The Ch. 8 example, printed p. 144.** The text reads "at least N = 2^7 = 128
  grid points are required" off a contour plot. Constraint 2 with
  `delta1 = 66.7 um` and `delta_n = 133 um` asks for `N >= 142.8`, that is
  2^7.16. So the printed example breaks its own constraint 2 by 11 percent.
- **`r0_sw = 17.7 cm`, Sec. 9.5.1, printed p. 176.** Listing 9.5 with the same
  inputs gives 12.66 cm, and Problem 2, printed p. 183, confirms the
  `(3/8)^(-3/5)` factor. The printed `sigma_chi,sw^2 = 0.436` DOES reproduce
  (0.4365), so the odd number is `r0_sw` alone. Do not calibrate anything
  against 17.7 cm. Gap S-28.
