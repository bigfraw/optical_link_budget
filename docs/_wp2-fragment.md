# WP2 fragment — Fresnel propagators (Ch. 6) and partial propagation (Ch. 8)

A later step merges these rows into `docs/schmidt-crosscheck.md`. Do not edit
that file from here.

Source of every equation:
    J. D. Schmidt, "Numerical Simulation of Optical Wave Propagation with
    Examples in MATLAB", SPIE Press Monograph PM199 (2010).
    DOI: 10.1117/3.866274

---

## Table 1 — forward map (olb code to book), WP2 rows

| olb id | location (file:line) | quantity | book eq | printed p | pdf p | status | note |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `Forvard` | `olb/waveoptics/propagators.py:57` | The angular-spectrum propagator, m = 1 | (6.31), (6.32) | 95 | 108 | checked | The transfer function in the docstring, `exp(i k z) exp(-i pi lam z f^2)`, is EXACTLY Eq. (6.32). olb KEEPS the piston factor `exp(i k z)`; `schmidt.fresnel.angular_spectrum` drops it, as Listing 6.5, printed p. 102, does. So the two differ by one constant phase. The irradiance is the same. |
| `Forvard` | `olb/waveoptics/propagators.py:121` | The transform pair | (2.6), (2.9) | 16, 17 | 29, 30 | checked | Duplicate of the WP1 row. A bare `fft2`/`ifft2` with a sign-alternation trick in place of `fftshift`. The missing `dx^2` and `(N df)^2` cancel inside one propagator. |
| `Forvard` grid rule | `olb/waveoptics/propagators.py:69` | "give the grid a side of about 8 times the largest beam radius" | (7.14), (7.20) | 119, 120 | 132, 133 | conflict | The olb rule is a rule of thumb with no source. Constraints 1 and 2 give the real bound, and both need D1, D2 and z. See gap S-08. |
| `Fresnel` | `olb/waveoptics/propagators.py:136` | The CONVOLUTION form of the Fresnel integral | (6.6) | 88 | 101 | checked | olb convolves on a grid of twice the side and integrates the kernel over one pixel in closed form (C and S). The book does NOT do that: it multiplies by the analytic transfer function of Eq. (6.49), printed p. 99. The two solve the same Eq. (6.6). The pixel integral is a REFINEMENT of the book, not a departure from it. |
| `Fresnel` minimum distance | `olb/waveoptics/propagators.py:149` | "z comparable with, or less than, the aperture" | (7.41), (7.42) | 123 | 136 | checked | The docstring cites Ch. 7 for the minimum distance. Eqs. (7.41) and (7.42) give the number: `z >= D1 dx1 R /(lam R - D1 dx1)`, and `z >= D1 dx1 / lam` for a flat source. olb states the rule in words only. See gap S-09. |
| `Fresnel` doubled grid | `olb/waveoptics/propagators.py:139` | The zero-padded convolution | — | — | — | no book equation | The book has NO zero-padded convolution. It controls the wrap with the ABSORBING BOUNDARY of Ch. 8, Eq. (8.1), printed p. 134, and with the grid-size constraint of Eq. (7.20). Two different cures for one problem. |
| `GForvard` | `olb/waveoptics/propagators.py:253` | The analytic ABCD Gaussian route | (6.70), (6.77), (6.81) | 103, 104 | 116, 117 | partial | The book gives the ray matrices and the generalized Huygens-Fresnel (ABCD) INTEGRAL. It does NOT give the closed-form q-parameter transform `q2 = (A q1 + B)/(C q1 + D)` that `_ABCD` uses. olb cites Siegman, ISBN 978-0935702118. Eq. (6.77) assumes azimuthal symmetry, which a pure Gaussian meets. |
| `Lens` | `olb/waveoptics/lenses.py:70` | The thin lens as a quadratic phase | (6.76) | 104 | 117 | checked | Eq. (6.76) is the thin-lens ray matrix. The phase form is the operator `Q[-1/f, r]` of Eq. (6.7), printed p. 89. The book states the lens phase delay at Sec. 6.5, printed p. 104. |
| `LensFresnel` | `olb/waveoptics/lenses.py:173` | The co-moving (spherical) grid | — | — | — | NO book equation | The book NAMES the Coles and Rubio angular-grid method (Ch. 6, text, printed p. 87) and then does NOT develop it. Schmidt's own answer to the same problem is the SCALING PARAMETER m of Eq. (6.65), printed p. 100, on a FLAT grid. The olb `Lens -> LensFresnel -> Convert` recipe has no Schmidt equation to check against. See gap S-05. |
| `Convert` | `olb/waveoptics/lenses.py:233` | Return from spherical to a flat grid | — | — | — | NO book equation | Same as `LensFresnel`. The book never leaves the flat grid. |
| `LensFresnel` magnification | `olb/waveoptics/lenses.py:187` | `fA = z/(m - 1)`, grid scale `(f - z)/f` | (6.24), (6.52) to (6.54) | 94, 99, 100 | 107, 112, 113 | partial | The olb `m` and the book `m` mean the SAME thing: the ratio of the two grid pitches. But olb gets it from a virtual lens and the book gets it from a free parameter in the exponent. The two routes are not the same algorithm. |
| `split_step` | `olb/waveoptics/turbulence/splitstep.py:94` | The partial-propagation loop | (8.18) | 139 | 152 | conflict | olb calls `Forvard` (m = 1) on ONE flat grid for every step. Eq. (8.18) gives each step its OWN pitch, from the linear rule of Eq. (8.8), printed p. 136, and its own magnification m_i. So olb cannot grow the grid with the beam. See gap S-06. |
| `super_gaussian_boundary` | `olb/waveoptics/turbulence/splitstep.py:29` | The absorbing boundary | (8.1) | 134 | 147 | conflict | The FORM matches Eq. (8.1), but the parameterisation and the shape do not match Listing 8.1. See Table 3 below and gap S-07. |
| `forvard_max_z` | `olb/waveoptics/grid.py:209` | `z_max = N dx^2 / lam` | (7.59), (8.24) | 127, 144 | 140, 157 | checked | This IS constraint 4 of Eq. (7.59), rearranged for m = 1: `N >= lam z/(dx1 dx2)` with `dx1 = dx2 = dx` gives `z <= N dx^2/lam`. It is also Eq. (8.24), the step cap, with `min(dx1, dxn) = dx`. Fill the Table 3 row: the constant is DERIVED, not a guess. |

---

## Table 2 — gaps and suggestions, WP2 rows

| gap id | book section | book eq | capability | target module | priority |
| --- | --- | --- | --- | --- | --- |
| S-05 | 6.4, 6.5 | (6.65) | The SCALED flat-grid propagator, as the alternative to the olb `Lens -> LensFresnel -> Convert` co-moving route. `schmidt.fresnel.angular_spectrum(..., dx2=)` now gives it. The two routes solve the same problem, and NO test compares them. Compare them on the 600 km uplink of `olb/waveoptics/run.py`. | an example in WP5 | high |
| S-06 | 8.3 | (8.18), (8.8) | The PER-PLANE grid pitch in the turbulent split step. `olb/waveoptics/turbulence/splitstep.py` holds one flat pitch for the whole path, so a diverging uplink beam must fit the SOURCE grid at the RECEIVER. `schmidt.fresnel.partial_propagations` now shows the book's linear pitch rule. Wiring it into `splitstep.py` is an owner decision, because it moves every turbulent number. | `olb/waveoptics/turbulence/splitstep.py` | high |
| S-07 | 8.1 | (8.1), Listing 8.1 | The absorbing boundary of olb is a DIFFERENT shape from the book's. olb: power 8, a taper band of 0.125 of the half-side, so the mask is exactly 1.0 out to 0.875 of the half-side and `exp(-1)` at the middle of an edge. Book: power 16, sigma = 0.47 N pixels, so the mask is 0.99999 at 0.2 N and 0.0678 at the middle of an edge. The book's boundary bites HARDER at the edge and it has no flat-then-taper break. Decide which one olb keeps. | `olb/waveoptics/turbulence/splitstep.py` | medium |
| S-08 | 7.2, 7.3 | (7.14), (7.20) | Constraints 1 and 2 as CODE. `olb/waveoptics/grid.py` sizes the grid from a `guard` of 4 and a `pixels_per_feature` of 16, with no citation. The book gives two inequalities in D1, D2, z, dx1 and dx2. | `olb/waveoptics/schmidt/sampling.py` (WP3), then `olb/waveoptics/grid.py` | high |
| S-09 | 7.3.1.2 | (7.41), (7.42) | The MINIMUM one-step distance as a number. `olb/waveoptics/propagators.py:149` states the rule in words: "z comparable with, or less than, the size of the aperture". The book gives `z >= D1 dx1 / lam` for a flat source. | `olb/waveoptics/schmidt/sampling.py` (WP3) | medium |
| S-10 | 6.6 | (6.82), (6.89), (6.92) | The MODEL point source. A true delta has infinite bandwidth, so the book replaces it with a sinc that gives the wanted windowed target field. olb has no point source: `olb/waveoptics/sources.py` gives `GaussBeam`, `PlaneWave`, `CircAperture` and `CircScreen` only. A retro link or a beacon may need one. | `olb/waveoptics/schmidt/` (a later work package) | low |
| S-11 | 6.5 | (6.77), (6.80), (6.81) | The general ABCD propagator for a NON-Gaussian field. `GForvard` handles a pure Gaussian only, and it raises on any other field. Eq. (6.81) gives the ABCD transfer function, which works on any field. | not needed yet | low |

---

## Table 3 — constants ledger, WP2 rows

The book values for the absorbing boundary. Both rows were `flagged` with an
empty book column.

| olb constant | olb value | location | book quantity | book value | book eq | printed p | pdf p | status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `super_gaussian_boundary` `power` | 8 | `olb/waveoptics/turbulence/splitstep.py:29` | The super-Gaussian exponent n of Eq. (8.1) | **16** | (8.1), Listing 8.1 line 12 | 134, 142 | 147, 155 | CONFLICT. The book RUNS `sg = exp(-nsq.^8/w^16)`, which is `exp(-(r/w)^16)`, so n = 16. Figure 8.1, printed p. 134, also plots n = 16. Eq. (8.1) itself only needs n > 2, so the olb value of 8 is ALLOWED by the equation but it is not the book's number. The book also records that Flatte and others used n = 8 (Ch. 8, text, printed p. 134), so 8 has a source in the literature, not in Schmidt's own runs. |
| `super_gaussian_boundary` `width_frac` | 0.125 | `olb/waveoptics/turbulence/splitstep.py:29` | The half-width sigma of Eq. (8.1) | **0.47 N pixels** (Listing 8.1); **0.45 L** (Fig. 8.1) | (8.1), Listing 8.1 line 11 | 134, 142 | 147, 155 | CONFLICT of PARAMETERISATION. The book states one half-width sigma in PIXELS, measured from the centre. olb states a taper BAND width as a fraction of the half-side, with a hard flat region inside it. The two cannot be converted. Book at the middle of an edge: `exp(-(0.5/0.47)^16) = 0.0678`. olb at the middle of an edge: `exp(-1) = 0.368`. The book absorbs about 5 times harder there. |
| `forvard_max_z` | `z_max = N dx^2 / lam` | `olb/waveoptics/grid.py:209` | Constraint 4 with m = 1; the step cap | `N >= lam z/(dx1 dx2)`; `Delta_z_max = min(dx1,dxn)^2 N/lam` | (7.59), (8.24) | 127, 144 | 140, 157 | RESOLVED. The olb formula is the book formula. Change the status from `flagged` to `cited`. |

Note for the Fig. 8.1 caption: it prints `sigma = 0.45 L` and `n = 16`, while
Listing 8.1 prints `w = 0.47*N`. L and N are the same length in pixel units, so
the two numbers differ by 4%. `schmidt.fresnel.super_gaussian_absorber` takes
0.47 as the default, because that is the value the book RUNS.

---

## WP2 note

### Built

`olb/waveoptics/schmidt/fresnel.py`. The module gives five names. It imports
numpy and `schmidt.fourier` only.

- `one_step_fresnel(Uin, wavelength, dx1, z) -> (Uout, dx2)` — Ch. 6, Sec.
  6.3.1. The transform form is Eq. (6.5), printed p. 88. The operator chain is
  Eq. (6.15), printed p. 90. The FIXED observation pitch
  `dx2 = lambda z /(N dx1)` is Eq. (6.16), printed p. 90.
- `two_step_fresnel(Uin, wavelength, dx1, dx2, z) -> Uout` — Ch. 6, Sec. 6.3.2.
  The operator chain is Eq. (6.18), printed p. 93. The pitch chain is
  Eqs. (6.19) to (6.21), printed p. 94. The scaling parameter is Eq. (6.24) and
  the intermediate plane is Eq. (6.25), both printed p. 94.
- `angular_spectrum(Uin, wavelength, dx, z, dx2=None) -> Uout` — Ch. 6,
  Sec. 6.4. With `dx2=None` it is the baseline form, Eqs. (6.31) and (6.32),
  printed p. 95. With `dx2` it is the SCALED form, Eq. (6.65), printed p. 100.
- `super_gaussian_absorber(n, sigma_frac=0.47, power=16)` — Ch. 8, Eq. (8.1),
  printed p. 134, with the reference values of Listing 8.1, printed p. 142.
- `partial_propagations(Uin, wavelength, dx1, dxn, z_planes, absorber=None)` —
  Ch. 8, Sec. 8.3. The general chain is Eq. (8.18), printed p. 139. The linear
  pitch rule is Eq. (8.8), printed p. 136. The cancellation of the middle
  quadratic phases is Eqs. (8.14) to (8.16), printed p. 138.

Every kernel docstring carries a VALIDITY paragraph. It names four things: the
Fresnel and paraxial condition (Ch. 1, Eqs. (1.49), (1.50) and (1.57), printed
pp. 8 and 10; Ch. 6, text, printed p. 87), what the kernel fixes about the
output pitch, why the kernel aliases, and WHICH Chapter 7 constraint governs it
by equation number. No constraint is implemented here. `sampling.py` (WP3) owns
the tests.

### Self-check numbers

`python -m olb.waveoptics.schmidt.fresnel`. The common geometry is N = 1024,
lambda = 1 um, dx1 = 80 um, z = 10 m, W0 = 3 mm. The one-step pitch of
Eq. (6.16) is then 122.070 um, so m = 1.5259 and the three Chapter 6 kernels
share one output grid.

| check | measured | target |
| --- | --- | --- |
| `one_step_fresnel` against the Gaussian closed form | 4.857e-16 | 1e-12 |
| `two_step_fresnel` against the Gaussian closed form | 4.090e-16 | 1e-12 |
| `angular_spectrum`, scaled (m = 1.5259), against the closed form | 5.405e-16 | 1e-12 |
| `angular_spectrum`, baseline (m = 1, z = 4 m), against the closed form | 4.623e-16 | 1e-12 |
| one-step against two-step | 5.005e-16 | 1e-12 |
| one-step against angular spectrum | 6.867e-16 | 1e-12 |
| two-step against angular spectrum | 6.867e-16 | 1e-12 |
| `partial_propagations`, 6 planes, against one angular-spectrum step | 9.800e-16 | 1e-12 |
| `partial_propagations` against the closed form | 7.900e-16 | 1e-12 |
| the absorber does not touch the beam | 2.822e-15 | 1e-9 |

The absorber values: centre = 1.000000 exactly; the smallest value inside a
pixel radius of 0.2 N is 0.999999; the value at the middle of an edge is
0.067796, which equals the book value `exp(-(0.5/0.47)^16) = 0.067796` to 1e-12;
the mask never grows along a radius.

Every error is at the floor of double precision. That is the correct answer, not
a loose test: a Gaussian on a well-sampled grid has a spectrum far inside the
band, so the discrete Fresnel transform has no aliasing left to make an error.

The self-check also prints the step cap of Eq. (8.24) for the common geometry:
`Delta_z_max = 6.554 m`, so `n >= 3` planes. The run used 6.

### Decisions

- **The piston phase `exp(i k z)` is DROPPED from every kernel.** The book's own
  Listings 6.1, 6.3 and 6.5 (printed pp. 91, 96 and 102) drop it, and Listing
  8.1 (printed p. 142) drops it too. The factor is constant across a plane, so
  it changes no irradiance and no relative phase. Dropping it in all four
  kernels is what lets them agree with each other to 1e-16. NOTE: olb
  `Forvard` KEEPS `exp(i k z)`. Any future comparison must remove one or add the
  other.
- **The SCALED angular spectrum is BUILT, not skipped.** The task allowed the
  baseline `m = 1` form only. Eq. (6.65) reduces to Eqs. (6.31) and (6.32) when
  m = 1, so building the general form and defaulting `dx2=None` costs four lines
  and gives both. The tracker calls Eq. (6.50) "the workhorse of Chs. 7 to 9",
  and gap S-05 needs it to compare against the olb co-moving route.
- **`partial_propagations` takes `(dx1, dxn, z_planes)`, not a per-plane pitch
  list.** Eq. (8.8), printed p. 136, DERIVES the per-plane pitch from the two
  end pitches and the fractional distance. A caller cannot choose an arbitrary
  pitch per plane: Eq. (8.15), printed p. 138, needs the linear rule for the
  middle quadratic phases to cancel. So a `dx_planes` argument would let a
  caller break the algorithm. The signature follows Listing 8.1.
- **`z_planes` holds the FULL list, and the first value is normally 0.** The
  book's MATLAB argument starts at the second plane and line 14 of Listing 8.1
  prepends a zero. The full list is clearer at a call site.
- **The absorber goes on EVERY plane after the first, the observation plane
  included.** Eq. (8.18), printed p. 139, has `A[r_i+1]` for i = 1 to n-1, and
  Listing 8.1 line 38 applies `sg` inside the loop. `absorber=None` gives the
  pure vacuum result, which is what the self-check compares.
- **`two_step_fresnel` REFUSES m = 1.** Eq. (6.25), printed p. 94, then puts the
  intermediate plane at infinity. Table 6.2, printed p. 95, records that case.
  Use `angular_spectrum` for m = 1.
- **`two_step_fresnel` takes the MINUS branch of Eq. (6.25).** Listing 6.3,
  line 13, printed p. 96, uses `Dz1 = Dz/(1 - m)`. Eq. (6.30), printed p. 95,
  proves that the plus branch gives the same magnitude of m. A step distance may
  be negative, and the code does not refuse it.
- **This module holds NO sampling test.** `sampling.py` (WP3) owns the four
  constraints. The docstrings name the governing constraint by equation number
  and stop there. A caller that breaks a constraint gets an aliased result and
  no warning. The module docstring says so.
- **The Tukey window of Eq. (8.2), printed p. 134, is NOT built.** Nothing needs
  it.
- **The point source of Sec. 6.6 is NOT built.** See gap S-10.
- **The ABCD route of Sec. 6.5 is NOT built.** `GForvard` and `Lens` already
  cover the olb need. See gap S-11.

### The book would not give

- **A numerical threshold for the Fresnel approximation.** Chapter 1 defines the
  paraxial approximation as `cos alpha ~ 1` and `cos beta ~ 1` (Eqs. (1.49) and
  (1.50), printed p. 8), and Chapter 6 states that the approximation "is a very
  good one" for parallel planes (text, printed p. 87). Neither gives an
  inequality in D, z and lambda. The docstrings say that the book gives none.
- **A sampling constraint set for TWO-STEP propagation.** Section 7.3.1 analyses
  ONE step (Sec. 7.3.1.1, printed p. 120). Chapter 7 never revisits Sec. 6.3.2.
  The `two_step_fresnel` docstring tells the caller to apply the one-step rules
  to each of the two steps, with the pairs `(dx1, dxa, z1)` and `(dxa, dx2,
  z2)`. That is our reading, not a printed rule.
- **A Gaussian-beam TEST CASE.** The book DOES give the closed form: Ch. 1,
  Eqs. (1.53) to (1.56), printed p. 9. But it never uses it as a numerical
  target. Throughout Chs. 6 to 8 it compares against Fresnel diffraction from a
  SQUARE aperture (Ch. 1, Eq. (1.60), printed p. 11). The self-check uses the
  Gaussian instead, because a Gaussian has no truncation error on a large grid,
  so it isolates the kernel from the grid. Andrews and Phillips, 2nd ed. (2005),
  DOI 10.1117/3.626196, Ch. 4, Eqs. (37) and (38), printed p. 93, print the same
  solution.
- **A tolerance for any comparison.** The book compares its figures by eye
  ("the comparison is very close", printed p. 92). The self-check sets its own
  targets.
- **A rule that ties the absorber half-width to the region of interest.** The
  values 0.47 N and 16 are the book's own numbers in Listing 8.1. Section 8.1
  gives the reason for a super-Gaussian in words only: "we must be careful not
  to alter light in the central region of the grid" (printed p. 134).
- **Any equation for a spherical or co-moving grid.** Chapter 6 NAMES the Coles
  and Rubio angular-grid method (text, printed p. 87) and then does not develop
  it. The olb `LensFresnel` and `Convert` route has no Schmidt equation to check
  against. Schmidt's own answer to the same problem is the scaling parameter m
  on a flat grid. See gap S-05.
