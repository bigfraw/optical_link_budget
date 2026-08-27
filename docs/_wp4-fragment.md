# WP4 fragment — Chapter 9, turbulence and screens

This file is a FRAGMENT. Merge it into `docs/schmidt-crosscheck.md` in WP6. Do
not edit the tracker from WP4.

Built module: `olb/waveoptics/schmidt/turbulence.py`.

Source of every equation below:

    Schmidt (2010), DOI 10.1117/3.866274, Ch. N, Eq. (nn), printed p. NNN

---

# Table 1 — forward map (olb code to Chapter 9)

| olb id | location (file:line) | quantity | book eq | printed p | pdf p | status | note |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `phase_screen` PSD | `olb/waveoptics/turbulence/screens.py:81` | modified von Karman PHASE PSD, `0.023 r0^(-5/3) exp(-(f/fm)^2) / (f^2+f0^2)^(11/6)` | (9.51), (9.52) | 161 | 174 | checked | The expression is the book's, and `fm = 5.92/(2 pi l0)`, `f0 = 1/L0` are the book's too. olb defaults `l0 = 1e-6 m` and `L0 = 1e6 m`, not 0 and infinity, so the numbers are the Kolmogorov ones to 12 digits. The constant agrees: `0.49 (2 pi)^(-5/3) = 0.02290` against the printed 0.023 (0.42% apart). |
| `phase_screen` FT draw | `olb/waveoptics/turbulence/screens.py:81` (aotools `ft_phase_screen`) | the Fourier-series screen | (9.78)–(9.80), Listing 9.2 | 167 | 180 | checked | The `ft_phase_screen` of the new module and the aotools one give the SAME mean structure function to 3 digits: ratio to Eq. (9.44) is 0.636 / 0.540 / 0.422 / 0.347 / 0.286 at r/r0 = 1 / 2 / 4 / 6 / 8 (24 screens, N = 256, r0 = 10 px, direct-difference estimator). So the aotools FT screen IS Listing 9.2. |
| `phase_screen(subharmonics=True)` | `olb/waveoptics/turbulence/screens.py:81` (aotools `ft_sh_phase_screen`) | the subharmonic low-frequency screen | (9.81), Listing 9.3 | 169, 170 | 182, 183 | checked | The two do NOT agree. On the run above, the book form gives 0.863 / 0.826 / 0.783 / 0.760 / 0.741 of theory, and aotools gives 0.824 / 0.778 / 0.725 / 0.692 / 0.661. The book form is 5 to 12% closer to Eq. (9.44) at every separation. See Table 2, row S-12. |
| `screen_r0` | `olb/waveoptics/turbulence/screens.py:56` | `r0_i = (0.423 k^2 Cn2_i dz_i)^(-3/5)` | (9.70) | 165 | 178 | checked | Exact match, constant included. olb cites Fried and Andrews Ch. 12; the book credits Roggemann et al., DOI 10.1364/AO.34.004037. Add the Schmidt citation. |
| `_composite_r0` | `olb/waveoptics/turbulence/sampling.py:198` | `r0 = (SUM r0_i^(-5/3))^(-3/5)` | (9.71) | 165 | 178 | checked | Exact match. It is the PLANE-wave composite. The book also gives the spherical one, Eq. (9.72), which olb has no name for. |
| `_screen_rytov` | `olb/waveoptics/turbulence/sampling.py:175` | one screen's path weight | (9.63), (9.73) | 163, 165 | 176, 178 | checked | olb computes `2.25 k^(7/6) (INT Cn2 dz) (z - z_i)^(5/6)`, which is the plane-wave RYTOV variance `sigma_R^2`. The book's per-screen quantity is the LOG-AMPLITUDE variance `sigma_chi^2`, constant 0.563. The ratio is `2.25/0.563 = 3.997`. Measured in the self-check: 3.9994. See Table 3. |
| `sigma2_r_screen_max` | `olb/waveoptics/turbulence/sampling.py:107` | the per-screen cap | Listing 9.5, lines 37, 38 | 175 | 188 | checked | The book caps `sigma_chi^2` at `rmax = 0.1`. olb caps `sigma_R^2 = 4 sigma_chi^2` at 0.05 / 0.10 / 0.25. See Table 3 for the factor analysis. |
| extent rule, the scattering cone | `olb/waveoptics/turbulence/sampling.py:442` | `2 (lambda/r0) z` added to the grid side | (9.84), (9.85) | 173 | 186 | checked | The added term is `c lambda dz / r0` with `c = 2`, which is the book's low value. Listing 9.6, line 2, printed p. 177, uses `c = 2` too. The book states `c = 2` holds 97% of the light and `c = 4` holds 99% (text below Eq. (9.85), printed p. 173). BUT olb adds the blur to the grid SIDE. The book adds it to D1' and D2' and then feeds constraints 1 to 3. Different route, same constant. |
| pixel rule, `pixels_per_r0` | `olb/waveoptics/turbulence/sampling.py:451` | `dx <= r0_total / pixels_per_r0` | Sec. 9.4 text | 172 | 185 | checked | The book gives the rule of Johnston and Lane, DOI 10.1364/AO.39.004761: pick the pitch at which the phase step between two adjacent samples stays below pi for more than 99.7% of the draws. With Eq. (9.44) that reads `3 sqrt(6.88 (dx/r0)^(5/3)) <= pi`, so `dx <= 0.332 r0`, that is **3.01 pixels per r0**. olb's `standard` preset value 3 lands on it. |
| pixel rule, the Fresnel scale | `olb/waveoptics/turbulence/sampling.py:454` | `dx <= sqrt(lambda z)/2` | Sec. 9.4 text | 172 | 185 | checked | olb ALREADY has the book's scintillation pitch rule, exactly. It cites Andrews Ch. 8 for it. The rule is Schmidt Sec. 9.4, printed p. 172, from Johnston and Lane. The tracker glossary row `sqrt(lambda z) (172)` says "No olb name ... uses `pixels_per_r0` only". That row is WRONG. Fix it in WP6. |
| `_merge_layers` | `olb/waveoptics/turbulence/sampling.py:209` | where the screens go, and what each carries | (9.65) | 164 | 177 | gap | olb groups adjacent Cn2 layers under the Rytov cap, and it BAILS OUT to one screen per layer when the merge undershoots `min_screens`. It matches NO moment of Eq. (9.65). See Table 2, row S-07, and the WP7 verdict. |
| `turbulent_grid` | `olb/waveoptics/turbulence/sampling.py:366` | the grid sizer | (9.86)–(9.88) | 173, 174 | 186, 187 | gap | olb applies NONE of the three turbulent geometry constraints. It sizes the side from a beam-plus-cone rule and the pixel from `r0` and the Fresnel scale, then it rounds N up to a power of two. See Table 2, row S-06. |
| `super_gaussian_boundary` | `olb/waveoptics/turbulence/splitstep.py:29` | the absorbing boundary | (8.1); Listing 9.7, line 19 | 134, 179 | 147, 192 | checked | Eq. (8.1) gives the SHAPE `exp(-(r/sigma)^n)`, `n > 2`, and no numbers. Listing 9.7 gives the numbers, and they are not olb's. See Table 3. |
| `split_step` max hop | `olb/waveoptics/turbulence/splitstep.py:170` | `max_step = N dx^2 / lambda` | (9.89) | 174 | 187 | checked | Exact match. It repeats Ch. 8, Eq. (8.24), printed p. 144. olb cites Ch. 6; the turbulent statement is Eq. (9.89). |
| `split_step` loop | `olb/waveoptics/turbulence/splitstep.py:94` | the split-step chain | (9.1)–(9.3) | 150 | 163 | checked | olb hops to a screen, applies the screen, and hops on. The book applies the screen AT each partial-propagation plane, Eq. (9.3), printed p. 150. The two agree when the screens sit at the plane positions. olb's screens sit at slab CENTRES, so the two differ by half a slab. The book does not treat that case. |
| `min_screens` | `olb/waveoptics/turbulence/sampling.py:107` | the screen-count floor | — | — | — | gap | Chapter 9 gives NO such floor. See Table 3 and the WP7 verdict. |

# Table 2 — gaps and suggestions

Numbering continues from row S-04 of the tracker.

| gap id | book section | book eq | capability | target module | priority |
| --- | --- | --- | --- | --- | --- |
| S-05 | 9.4 | Sec. 9.4 text, printed p. 172 | The Johnston and Lane PHASE pitch rule, `dx <= 0.332 r0`. olb's `pixels_per_r0` is a bare preset integer with a Martin and Flatte citation and no derivation. The book's prose plus Eq. (9.44) give the number. Built here as `phase_pitch_max`. | `olb/waveoptics/turbulence/sampling.py` | medium |
| S-06 | 9.4, 9.5.2 | (9.86), (9.87), (9.88) | The three turbulent geometry constraints, and the blurred extents D1', D2' of Eqs. (9.84) and (9.85). olb checks none of them, so a bad pitch pair gives no warning. Built here as `constraint1_pitch_max`, `constraint2_n_min`, `constraint3_pitch_range`, `blurred_extent`. | `olb/waveoptics/turbulence/sampling.py`, or a validation example | high |
| S-07 | 9.2.5 | (9.65) | The layered-atmosphere MOMENT rule for `0 <= m <= 7`. It is the only screen-placement rule that Chapter 9 gives. `_merge_layers` satisfies no part of it. Built here as `profile_moments`, `layer_moments`, `moment_error`. | `olb/waveoptics/turbulence/sampling.py` | high (the WP7 hinge) |
| S-08 | 9.2.5, 9.5.1 | (9.75), Listing 9.5 | The constrained least-squares solve for the screen `r0` values from a target `r0_sw` and `sigma_chi,sw^2`. olb never solves for a screen strength; it takes the Cn2 layers as given. Built here as `screen_strengths` and `max_screen_strength`. | `olb/waveoptics/turbulence/sampling.py` | medium |
| S-09 | 9.5.5 | (9.32), (9.44) | The observation-plane coherence factor as the end-to-end verification of a turbulent run. It is row S-02 of the tracker seen from Chapter 9. `properly_sampled_checklist` names it as an advisory step; nothing measures it. | a later work package | medium |
| S-10 | 9.4 | — | `QualityPreset.fresnel_weight_min` (`olb/waveoptics/turbulence/sampling.py:450`) exempts a weak screen from the Fresnel pitch rule. Chapter 9 states NO such exemption; Sec. 9.4 applies the rule to every step. The exemption is a real cost saver, so keep it, but mark it as an olb rule, not a book rule. | `olb/waveoptics/turbulence/sampling.py` | low |
| S-11 | 9.4 | — | Documentation only. The tracker glossary row for `sqrt(lambda z)` says olb has no such rule. `sampling.py:454` has it exactly. Correct the row. | `docs/schmidt-crosscheck.md` | low |
| S-12 | 9.3 | (9.81), Listing 9.3 | The aotools subharmonic screen reads 5 to 12% lower than the book's Listing 9.3 form across `r/r0 = 1` to 8. Both are below theory. Decide whether to keep aotools, to pass the book form from `schmidt.turbulence`, or to move to Johansson and Gavel, DOI 10.1117/12.177254, which the book calls the closest match (Ch. 9, text above Sec. 9.4, printed p. 172). | `olb/waveoptics/turbulence/screens.py` | medium |
| S-13 | 9.5.1 | Listing 9.5, lines 15 to 18 | A BOOK ERROR to record. Sec. 9.5.1, printed p. 176, prints `r0_sw = 17.7 cm` for the 50 km example. Listing 9.5 with the same inputs gives **12.66 cm**, and Problem 2, printed p. 183, confirms the `(3/8)^(-3/5)` factor. The printed `sigma_chi,sw^2 = 0.436` DOES reproduce (0.4365). So the printed `r0_sw` is the odd number. Do not calibrate anything against 17.7 cm. | none (a note) | — |

# Table 3 — constants ledger, the Chapter 9 rows

| olb constant | olb value | location | book quantity | book value | book eq | printed p | pdf p | status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `QualityPreset.min_screens` | 15 / 9 / 5 | `.../sampling.py:107` | a screen-count FLOOR | **the book would not give** | — | — | — | checked, no source |
| `QualityPreset.sigma2_r_screen_max` | 0.05 / 0.10 / 0.25 | `.../sampling.py:107` | `rmax`, the per-screen cap | 0.1 on `sigma_chi^2`, so **0.4 on `sigma_R^2`** | Listing 9.5, lines 37, 38 | 175 | 188 | checked, olb is 1.6x to 8x stricter |
| `QualityPreset.pixels_per_r0` | 4 / 3 / 2 | `.../sampling.py:107` | the phase pitch rule | **3.01** pixels per r0 | Sec. 9.4 text with (9.44) | 172, 160 | 185, 173 | checked, `standard` = 3 matches |
| `QualityPreset.guard` | 4 / 3 / 2 | `.../sampling.py:107` | the beam-radius margin | the book gives no such margin; it sizes from D1', D2' and constraint 2 | (9.84)–(9.87) | 173, 174 | 186, 187 | checked, different route |
| the scattering-cone factor 2 | 2 (hard-coded) | `.../sampling.py:442` | `c`, the blur sensitivity | **2 to 8; c = 2 holds 97%, c = 4 holds 99%**; Listing 9.6 uses 2 | (9.84), (9.85) | 173 | 186 | checked, olb sits at the book's low end |
| `super_gaussian_boundary` `power` | 8 | `.../splitstep.py:29` | the super-Gaussian order | **16** | Listing 9.7, line 19 | 179 | 192 | checked, olb is softer |
| `super_gaussian_boundary` `width_frac` | 0.125 | `.../splitstep.py:29` | the boundary half-width | `0.47 N dx`, that is 0.94 of the HALF-SIDE at the `exp(-1)` point | Listing 9.7, line 19 | 179 | 192 | checked, see the note below |
| `QualityPreset.boundary_width_frac` | 0.125 / 0.125 / 0.10 | `.../sampling.py:107` | same as above | same as above | Listing 9.7, line 19 | 179 | 192 | checked, no book value in this parameterisation |
| `QualityPreset.fresnel_weight_min` | 0.005 / 0.02 / 0.05 | `.../sampling.py:107` | a weak-screen exemption | **the book would not give** | — | — | — | checked, no source |
| `MAX_SCREENS` | 500 | `.../sampling.py:60` | a screen-count cap | **the book would not give** | — | — | — | checked, no source |

## The `rmax` versus `sigma2_r_screen_max` factor

Listing 9.5, lines 37 and 38, printed p. 175, read

    rmax = 0.1;
    x2 = rmax/1.33*(k/Dz)^(5/6) ./ A(2,:);

where `x` holds `r0_i^(-5/3)` and row 2 of `A` holds
`alpha^(5/6) (1-alpha)^(5/6)`. Multiply both sides by `A(2,i)`:

    1.33 k^(-5/6) z^(5/6) r0_i^(-5/3) alpha^(5/6) (1-alpha)^(5/6) <= 0.1

The left side is one term of Eq. (9.74), printed p. 165, which is the
SPHERICAL-WAVE LOG-AMPLITUDE variance `sigma_chi,sw^2` of Eq. (9.64), printed
p. 163. So `rmax = 0.1` bounds `sigma_chi^2`, NOT the Rytov variance. The
book's own text calls it "the overall Rytov number" (Sec. 9.5.1, printed
p. 176), and that phrase is loose.

olb's `_screen_rytov` (`sampling.py:175`) computes
`2.25 k^(7/6) (INT Cn2 dz) (z - z_i)^(5/6)`. Substituting Eq. (9.70) into the
book's plane-wave term, Eq. (9.73), gives
`0.563 k^(7/6) (INT Cn2 dz) (z - z_i)^(5/6)`. The ratio is

    2.25 / 0.563 = 3.997

The self-check measures 3.9994. So **olb's per-screen number is
`sigma_R^2 = 4 sigma_chi^2`, and the book's cap of 0.1 on `sigma_chi^2` is a
cap of 0.4 on olb's number.**

| preset | olb cap on `sigma_R^2` | the same as a cap on `sigma_chi^2` | against the book's 0.1 |
| --- | --- | --- | --- |
| `reference` | 0.05 | 0.0125 | 8x stricter |
| `standard` | 0.10 | 0.025 | 4x stricter |
| `rapid` | 0.25 | 0.0625 | 1.6x stricter |

Two more differences, both small:

- The book weights the screen for a SPHERICAL wave,
  `alpha^(5/6) (1-alpha)^(5/6)`. olb weights it for a PLANE wave,
  `(1-alpha)^(5/6)`. For a downlink slab the source is far away, so the
  plane-wave weight is the right one, and it is the LARGER of the two, so
  olb's choice stays conservative.
- The book applies the cap as an OPTIMISER bound while it solves for the
  screen `r0` values. olb applies it as a merge rule on a fixed Cn2 profile.

**Verdict on this constant: olb is conservative, and it is not wrong.** No
change is forced. If a run is too slow, `rapid` at 0.25 is still 1.6x inside
the book's guideline, and 0.4 is the book value.

## The absorbing boundary constants

Listing 9.7, line 19, printed p. 179, reads

    sg = exp(-(x1/(0.47*N*d1)).^16) .* exp(-(y1/(0.47*N*d1)).^16);

so the mask is SEPARABLE in x and y, the order is 16, and the `exp(-1)` point
sits at `x = 0.47 L`, that is 0.94 of the half-side. The mask first falls below
0.99 at 0.705 of the half-side, and it is 0.068 at the middle of an edge.

`olb.super_gaussian_boundary` is RADIAL, of order 8, exactly 1.0 inside 0.875
of the half-side, and `exp(-1)` at the edge. The two shapes are not the same
family, so the numbers do not map one to one. Recorded, not changed. Eq. (8.1),
printed p. 134, allows any order above 2 and gives no `sigma`.

---

# WP4 note

## Built

`olb/waveoptics/schmidt/turbulence.py`, twenty names, each with its chapter,
equation number and printed page:

**The spectra (Secs. 9.2.3, 9.3).** `phase_psd` (the one shared expression),
`kolmogorov_phase_psd` (9.49), (9.52), `von_karman_phase_psd` (9.50),
`modified_von_karman_phase_psd` (9.51), `kolmogorov_structure_function`
(9.44).

**The screens (Sec. 9.3).** `ft_phase_screen` (9.78) to (9.80) with Listing
9.2, `subharmonic_screen` (9.81) with Listing 9.3, `ft_sh_phase_screen` (the
sum).

**The per-screen bound (Sec. 9.2.5, Listing 9.5).** `screen_rytov_share`
(9.73), (9.74), `max_screen_strength` (Listing 9.5, lines 37 to 39),
`screen_strengths` (9.75), `composite_r0` (9.71), (9.72), `screen_r0` (9.70),
and the constants `RMAX = 0.1` and `WEAK_SIGMA2_CHI = 0.25`.

**The layer rule (Sec. 9.2.5).** `profile_moments`, `layer_moments`,
`moment_error`, all Eq. (9.65).

**The sampling bounds (Sec. 9.4).** `fresnel_pitch_max`, `phase_pitch_max`,
`blurred_extent` (9.84), (9.85), `constraint1_pitch_max` (9.86),
`constraint2_n_min` (9.87), `constraint3_pitch_range` (9.88),
`max_partial_step` (9.89), `min_planes` (9.90).

**The procedure (Sec. 9.5).** `properly_sampled_checklist`, which returns one
`(rule, satisfied, bound, actual, citation)` tuple per step. Its arguments are
plain numbers, and their names match `GridSpec` (`n`, `pixel_m`, `size_m`) and
`ScreenPlan` (`z_m`, `r0_m`, `r0_total_m`, `z_total_m`) one to one. It imports
no olb module outside `schmidt`.

## Measured

Self-check numbers, from `python -m olb.waveoptics.schmidt.turbulence`
(7.5 s):

- The modified von Karman PSD reduces to Kolmogorov for `L0 = inf`, `l0 = 0`
  to a relative error of 0.0. The Kolmogorov branch equals
  `0.023 r0^(-5/3) f^(-11/3)` to 4e-16. The angular constant converts:
  `0.49 (2 pi)^(-5/3) = 0.02290`, 0.42% from the printed 0.023.
- The mean structure function of 24 screens, `N = 512`, `r0 = 10` px, through
  `schmidt.fourier.structure_function` with a 1.2 m pupil, against Eq. (9.44):

  | r/r0 | subharmonic ratio | FT-only ratio |
  | --- | --- | --- |
  | 0.3 | 0.908 | 0.822 |
  | 0.5 | 0.898 | 0.797 |
  | 0.8 | 0.885 | 0.765 |
  | 1.2 | 0.870 | 0.733 |
  | 1.6 | 0.857 | 0.706 |
  | 3.2 | 0.815 | 0.625 |
  | 8.0 | 0.763 | 0.505 |

  The stated band is `r/r0 = 0.3` to 1.6, and the tolerance there is 0.85 to
  1.02. That is the band and the tolerance of the self-check of
  `olb/waveoptics/turbulence/screens.py`. The subharmonic screen lands inside
  it; the FT-only screen is below 0.85 everywhere and it falls to 0.505 at
  `r/r0 = 8`. The subharmonics do NOT close the gap at a large separation.
- Moment matching of a uniform 50 km slab, Eq. (9.65), `m = 0` to 7:
  4 screens at the 4-point Gauss-Legendre nodes match every moment to 1.2e-8
  (the trapezium error of the reference profile). 11 uniformly spaced screens
  of equal strength, which is the layering of the book's own worked example,
  miss `m = 2` by 5.0% and `m = 7` by 31.5%.
- The Sec. 9.5.1 example: `sigma_chi,sw^2 = 0.4365` against the printed 0.436.
  `r0_sw = 12.66 cm` against the printed 17.7 cm; see Table 2, row S-13. The
  11-screen `screen_strengths` solve returns `r0_sw` to 1.2e-5 and
  `sigma_chi,sw^2` exactly, and its
  largest screen share is 0.0745 against the cap of 0.1.
- The factor between the book's per-screen quantity and olb's: 3.9994.
- The Sec. 9.5.2 example: constraint 2 asks for `N >= 355.2`, and the book
  picks 512 because "the required number of grid points is more than 2^8"
  (printed p. 177). `min_planes` returns 2, and the book says two. Both match.
- The Johnston and Lane phase pitch rule gives 3.01 pixels per r0.

## Decisions

- ONE expression carries all three spectra. `kolmogorov_phase_psd`,
  `von_karman_phase_psd` and `modified_von_karman_phase_psd` are one-line
  wrappers over `phase_psd`. The book itself derives the three the same way
  (Eqs. (9.49) to (9.51), printed p. 161).
- `phase_psd` returns infinity at `f = 0` when `L0` is infinite, because the
  divergence is real physics. The two screen generators zero that sample, as
  Listing 9.2, line 16, and Listing 9.3, line 26, do.
- The subharmonic part is its OWN function, `subharmonic_screen`, and
  `ft_sh_phase_screen` sums the two. That is the structure of Listings 9.2 and
  9.3. A caller can measure the two parts apart, which the self-check does.
- `screen_strengths` calls `scipy.optimize.lsq_linear`, which solves the same
  bounded linear least-squares problem as the book's `fmincon`. The MATLAB
  listing is not ported.
- The moment rule is a CHECKER (`moment_error`), not a solver. Chapter 9 gives
  no solver for Eq. (9.65); it states the equality and then, at Sec. 9.5.5,
  printed p. 182, tells the reader to "adjust the values of z_i and dz_i
  attempting to match turbulence moments". WP7 owns the adjustment.
- `properly_sampled_checklist` exempts a screen of zero path weight from the
  scintillation pitch rule, because such a screen adds no scintillation. For a
  spherical wave those are the screens at `alpha = 0` and `alpha = 1`. The
  book states no exemption; it follows from Eq. (9.74).
- Steps 9.5.3, 9.5.4 and 9.5.5 come back as ADVISORY rows, with
  `satisfied = None`. They are procedures, not inequalities.
- Constraint 3 is not exempted. Ch. 7, Eq. (7.60), printed p. 129, exempts it
  when `1 + dz/R < D2/D1`. That belongs to WP3.

## The book would not give

- **A minimum screen count.** See the WP7 verdict below.
- **A tolerance for the structure function.** The book compares Fig. 9.3 and
  Fig. 9.9 by eye and calls the match "close". The self-check sets its own
  band, and it borrows the band and the tolerance from
  `olb/waveoptics/turbulence/screens.py` so that the two files agree.
- **An equation for the phase pitch rule.** Sec. 9.4, printed p. 172, states
  it in prose only: "phase differences less than pi in adjacent grid points
  occur more than 99.7% of the time". The algebra that turns that into
  `dx <= 0.332 r0` is ours. The 99.7% is a 3-sigma reading of a Gaussian
  phase difference, and the variance is Eq. (9.44).
- **An equation number for the scintillation pitch rule.** Sec. 9.4, printed
  p. 172, gives `sqrt(lambda z)/2` in prose, with no equation number.
- **A solver for Eq. (9.65).** See above.
- **A rule for a screen at a slab CENTRE.** The book puts one screen at each
  partial-propagation plane. olb puts a screen at the Cn2-weighted centre of a
  merged slab. Chapter 9 does not treat that placement.
- **A temporal axis.** Sec. 9.5.4, printed p. 179, states the frozen-flow
  method in prose and points to the Greenwood frequency. No equation, no code.

## WP7 GATE VERDICT — what Chapter 9 does and does not justify

The question for WP7 is: does Schmidt justify `QualityPreset.min_screens`
(15 / 9 / 5), and does Eq. (9.65) give a principled replacement for the
`_merge_layers` bail-out?

**1. Chapter 9 justifies NO screen-count floor. This is now settled.** Three
pieces of text bear on it, and none of them is a derivation:

- Eq. (9.90), printed p. 174, `n_min = ceil(dz / dz_max) + 1`, is a SAMPLING
  floor only. It comes from Constraint 4, which is a rule of the FFT method,
  not of the atmosphere. On the book's own 50 km example it returns 2.
- Sec. 9.2.5, printed p. 165, says "Using a typical number of phase screens,
  like 5-10, there are 10-20 unknown parameters". The "5-10" counts the
  UNKNOWNS of the underdetermined system of Eq. (9.75). It is not a floor, and
  the book gives no reason for it.
- Sec. 9.5.2, printed p. 177, says "the minimum number of planes is two, so we
  could use just one propagation. However, we use ten propagations (11 planes)
  to represent the atmosphere properly." The book gives NO formula, no
  criterion, and no convergence study for the 11.

So the 15 / 9 / 5 integers cannot be sourced to Schmidt, and neither can any
other integer. **The `min_screens` field stays uncited after WP4.** Two routes
remain open, and WP7 must pick one:

- (a) Delete the floor and let Eq. (9.65) set the count. See point 2.
- (b) Keep a floor and justify it with a CONVERGENCE SWEEP inside olb: hold the
  path fixed, sweep the screen count, and find where the measured scintillation
  index and the coherence factor stop moving. That is olb's own evidence, not
  the book's, and the docstring must say so.

**2. Eq. (9.65) IS a principled replacement for `_merge_layers`, and it is a
STRONGER rule than the one olb uses now.** The equation is

    INTEGRAL Cn2(z) z^m dz = SUM_i Cn2_i z_i^m dz_i,   0 <= m <= 7

It fixes both the screen POSITIONS and the screen STRENGTHS at once, and it is
the only screen-placement rule in the chapter. Three consequences for WP7:

- **It gives a real lower bound on the screen count: 4.** A layering with `n`
  screens has `2n` free numbers, and Eq. (9.65) is 8 equations. So `n = 4` is
  the smallest set that CAN satisfy it. The self-check shows that 4 screens at
  the 4-point Gauss-Legendre nodes match all 8 moments of a uniform slab
  EXACTLY (error 1.2e-8). This is a moment-matching bound, NOT a
  scintillation-fidelity bound: 4 screens match `r0`, `theta_0` and
  `sigma_chi^2`, and they say nothing about the irradiance PDF. It is
  nevertheless the first number in this whole area that follows from the book.
- **It decouples the screen count from the profile sampling, which is exactly
  the bug that `_merge_layers` has.** The moments of the CONTINUOUS profile do
  not depend on how finely `hs` samples it. So a 20-layer `DEFAULT_HS` and a
  200-layer real profile give the same target moments, and therefore the same
  screen count. That removes the "200 layers gives 200 screens" failure that
  `CLAUDE.md` records.
- **It condemns the layering that olb produces today, and the book's own worked
  example too.** 11 uniformly spaced equal screens on a uniform slab miss
  moment 2 by 5% and moment 7 by 31%. `_merge_layers` groups by Rytov weight,
  which is a `(1-alpha)^(5/6)` weighting, so it matches moment 0 (the Cn2
  integral) and nothing else.

**3. What Chapter 9 still does not settle for WP7.** Eq. (9.65) constrains the
layering, but it does not pick ONE layering:

- The chapter gives no solver, and no tolerance on the moment error.
- It gives no guidance on what happens when the moment-matched positions
  violate the per-screen `rmax` cap of Listing 9.5. The two rules can conflict
  on a strong path, and the book's own Listing 9.5 sidesteps this by FIXING the
  positions first and then solving only for the strengths under the cap. That
  is a defensible route for olb too, and it needs no moment machinery: pick the
  positions from the sampling rules, then solve Eq. (9.75).
- A generalized Gauss quadrature with weight `Cn2(z)` would satisfy Eq. (9.65)
  exactly with 4 nodes for ANY profile. That is the clean generalisation of
  the self-check case. The book does not name it, and it is NOT built here.
  Flag it as the leading candidate for WP7, not as a decision.

**Practical recommendation for WP7, on the evidence above.** Replace the
`_merge_layers` bail-out with a rule that (i) picks a screen count from the
larger of `min_planes` (Eq. (9.90)) and the per-screen `rmax` cap, with a hard
floor of 4 from the moment count; (ii) places and weights the screens to
minimise `moment_error`; and (iii) reports the achieved moment error in
`SamplingReport`, because Chapter 9 gives the equality and no tolerance. Do
NOT keep 15 / 9 / 5 with a Schmidt citation. The book does not support it.
