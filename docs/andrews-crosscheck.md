# Andrews cross-check tracker

This file tracks the sections of Andrews and Phillips that we must check our
equations against, or bring into our analytical modules. The user flags the
sections in batches. Each entry keeps a status.

Source book: Andrews and Phillips, *Laser Beam Propagation through Random
Media*, 2nd ed. (SPIE Press, 2005). DOI: 10.1117/3.626196. The equation numbers
are per chapter, so "Ch. 6, Eq. (108)" means equation (108) of Chapter 6. The
end of each chapter collates its equations. Use that collation as the quick
index.

## Status keys

- **flagged** — the user marked the section. We have not checked it yet.
- **checked** — we compared the section to the code. See the note for the
  result.
- **incorporated** — the section is now in a module. See the note for the
  module and the `docs/physics.md` section.

## Policy: strong fluctuation regime

For the strong fluctuation regime, we use numerical simulation, not an analytic
extrapolation of weak-fluctuation theory. So each flagged section that gives a
weak-fluctuation form must also record its regime limit (for example the
Rytov-variance bound). When a scenario passes that limit, the analytic Term is
not valid, and the budget must fall back to a numerical path. See also
`olb/assumptions.py`, which each Term uses to declare its turbulence regime.

## Batch 1 — terrestrial and slant-path Gaussian beam

| Section | Equations | Topic | Maps to | Status |
| --- | --- | --- | --- | --- |
| 6.8.1 | Ch. 6, Eq. (108); Eq. (40), (45) from 6.3.1 | Mean irradiance. Turbulence-induced beam spread for an arbitrary refractive-index structure function. | Beam spread / effective long-term beam radius. `olb/beam.py`, `olb/turbulence/gaussian_fried.py`, `docs/physics.md` §1 and §5e. | incorporated in part — `olb/turbulence/andrews/wander.py` holds W_LT (Ch. 6, Eq. (86)) and W_ST (Ch. 6, Eq. (100)), and `andrews/beam.py` holds the curvature-general beam parameters. Answered by Table 1 rows KR-12 and KR-13 (Ch. 6, Eqs. (109) and (111), both `exact`) and GF-01 to GF-04. Eq. (108) itself stays a gap: Table 2 row G-42. |
| 6.6.1 | Ch. 6, Eq. (88) | Beam wander. | Beam-wander variance that folds into the coupled-flux wander term and the Dios off-axis model. `olb/turbulence/uplink_flux.py`, `olb/turbulence/beam_wave_scintillation.py`, `docs/physics.md` §5c and §5d. | incorporated — `olb/turbulence/andrews/wander.py` holds the book beam-wander variance (Ch. 6, Eqs. (93) to (99)), the pointing error (Ch. 8, Eqs. (36) to (38)) and their slant forms (Ch. 12, Eqs. (50) and (53)). Answered by Table 1 row KR-04 (`other source`: the kernel constant is 2.07, from Dios Eq. (11); the book gives 7.25). C-01 is CLOSED on 2026-08-25: the kernel copies Dios correctly, so 2.07 stays. Belmonte 2000 prints the same 2.07 form, traces it to the Yura / Mironov-Nosov image-motion level arm, and validates it against a split-step simulation; the Andrews 7.25 is the beam-wave spectral-filter derivation and over-counts by 3.50. See the C-01 closure block and Table 2 rows G-38, G-44, G-49. |
| 6.8 (general) | slant-path extension | Extension to slant paths for an arbitrary Cn2. | Slant-path generalisation of the Gaussian-beam turbulence Terms. Ties to the CLAUDE.md "Next task" (curvature past the collimated case). | incorporated — `olb/turbulence/andrews/paths.py` holds the Chapter 12 slant forms (the path moments mu_0 to mu_3, the uplink and downlink indices, the coherence radius, the isoplanatic angle), and `andrews/beam.py` takes any input curvature f0, which closes olb gap 3 at the physics layer. Answered by Table 1 rows GF-15 to GF-21, of which GF-18 is `wrong` (mirrored path weight). See Conflicts C-02 and Table 2 rows G-43 and G-130. |
| 6.7 | temporal spectra | Temporal spectra of the beam parameters. | The planned temporal-vs-snapshot option. See the temporal-statistics side-step. | incorporated in part — `olb/turbulence/andrews/temporal.py` holds the irradiance temporal spectra (Ch. 8, Sec. 8.5; Ch. 9, Sec. 9.8), the quasi-frequency, the Greenwood frequency and the coherence time, and `andrews/distributions.py` holds the fade rate and the mean fade time. No Table 1 row exists, because the olb Terms still have no temporal axis. Row G-40 is closed; row G-41 (the Ch. 6.7 mutual-coherence temporal spectrum) stays open; see also G-11, G-67, G-99, G-115 and G-149. |
| 8.2 | scintillation index | Scintillation index for a tracked and an untracked Gaussian beam. Find the restrictions of weak-fluctuation theory for this case. | Scintillation index Terms; the weak/strong regime limit that sets when we switch to a numerical path. `olb/turbulence/plane_wave_scintillation.py`, `olb/turbulence/beam_wave_scintillation.py`, `docs/physics.md` §5b and §5d. | incorporated — `olb/turbulence/andrews/scintillation.py` holds the tracked and untracked Gaussian-beam index with its weak and all-regime branches (Ch. 8, Sec. 8.2; Ch. 9, Sec. 9.6), and `andrews/paths.py` holds the slant-path pair (Ch. 12, Eqs. (54) to (61)). Answered by Table 1 rows BW-08 to BW-16. UF-01 is fixed. PW-01 keeps 0.25 as a labelled house rule, and `olb/links/downlink.py` now uses that limit as the switch point to the gamma-gamma Term. TL-05 stays `wrong`. See Conflicts C-05 and Table 2 rows G-20, G-50 to G-53, G-84. |

## Batch 2 — Gaussian-beam angle of arrival / aperture tip-tilt

| Section | Equations | Topic | Maps to | Status |
| --- | --- | --- | --- | --- |
| (owner to specify) | (owner to specify) | Aperture angle-of-arrival "corrugation" tip-tilt of a Gaussian beam (the classic plane-wave form ~0.182*(D/r0)^(5/3)*(lambda/D)^2). | The second, smaller received tip-tilt term. `olb/turbulence/angle_of_arrival.py` `aperture_arrival_angle_variance` (a stub that raises NotImplementedError). The working received tip-tilt is the beam-wander term only. | incorporated — the owner DECIDED C-04 on the GRADIENT tilt. `olb/turbulence/andrews/structure.angle_of_arrival_variance` gives the book form, Ch. 6, Eq. (84), printed 201, with the inner-scale and outer-scale branches of Eq. (83), and the stub `aperture_arrival_angle_variance` now calls it with the same signature. The value is 0.174 (D/r0)^(5/3)(lambda/D)^2 per axis, NOT the Noll Zernike 0.182. Both docstrings say so. Note that `olb/turbulence/ao.py` still uses the Noll coefficients, so olb holds BOTH tilt conventions: a caller that adds the two must say which one it means. See Conflicts C-04 and Table 2 rows G-34 to G-36. |

### Notes for batch 2

- The working received tip-tilt used by the coupling Terms is the beam-wander
  arrival tilt, `wander_arrival_angle_variance` (Dios et al. 2004, DOI
  10.1364/AO.43.003866). It reuses the `beam_wander_variance` kernel. Its
  "radial (2-axis)" docstring is CONFIRMED correct on 2026-08-25 by Dios
  Eqs. (9) and (10), printed p. 3868. See Conflicts C-03.
- The aperture angle-of-arrival "corrugation" term is a separate, smaller
  contribution. `aperture_arrival_angle_variance` in
  `olb/turbulence/angle_of_arrival.py` now delegates to
  `andrews.structure.angle_of_arrival_variance` (the gradient-tilt form, C-04),
  so it no longer raises. It still feeds NO coupling Term (backlog 0-W3), so the
  received tip-tilt stays a lower bound.

### Notes for batch 1

- Section 6.8.1 gives the mean irradiance for an arbitrary structure function.
  Eq. (40) and Eq. (45) come from 6.3.1 and give the base Gaussian-beam wave
  parameters that Eq. (108) uses. So check Eq. (40) and Eq. (45) first.
- Section 6.8 extends the beam-spread and mean-irradiance forms to a slant path
  with an arbitrary Cn2. This is useful for the slant-path generalisation.
- Section 8.2 needs the weak-fluctuation limit written down next to the
  equation. When a terrestrial or slant scenario is in strong turbulence, the
  analytic scintillation index is not valid, and we use a numerical path (see
  the policy above).

---

# Merged cross-reference matrix

Seven readers (R1 to R7) read the book against the code. R1 read Ch. 3 and
Appendix III, R2 read Ch. 4 and Ch. 5, R3 read Ch. 6 and Ch. 7, R4 read Ch. 8,
R5 read Ch. 9, R6 read Ch. 10 and Ch. 11, R7 read Ch. 12 to Ch. 14. This
section merges their three tables into one matrix.

How to read the tables:

- `olb id` is a stable key `XX-nn`. The prefix names the module: `BM` beam.py,
  `DL` links/downlink.py, `RT` links/retro_space.py, `TL` links/terrestrial.py,
  `PT` models/pointing.py, `RS` results.py, `AA` angle_of_arrival.py,
  `AN` anisoplanatism.py, `AO` ao.py, `BW` beam_wave_scintillation.py,
  `GF` gaussian_fried.py, `PW` plane_wave_scintillation.py, `PR` profiles.py,
  `UF` uplink_flux.py, `KR` the external kernel in `my_analysis_modules`.
- `status` is `exact`, `reduction` (the olb form adds assumptions),
  `approximate` (a different fit to the same quantity), `other source` (a
  non-Andrews paper, with its DOI in the note), `unmatched` (not found) or
  `wrong` (the code disagrees with the book).
- Where two readers gave one location different statuses, the table keeps the
  STRICTER status and the Conflicts section holds the disagreement.
- A reader who could not find an equation inside the pages assigned to that
  reader marked it `unmatched`. Where a second reader then found the same
  equation in another chapter, the merged row carries the positive result and
  the note records the resolution.
- `printed p` is the arabic page in the book. `pdf p` is the page in
  `REFS/9780819478320.pdf`. Printed page N = PDF page N - 25.

## Findings summary

Ranked. The physics faults come first, then the citation faults.

**Physics faults (`wrong` rows)**

1. **PW-09 and KR-24 — the point plane-wave scintillation index has three
   wrong constants.** The code uses 0.54, 1.22 and 0.509; the book uses 0.49,
   1.11 and 0.51. Confirmed by 3 readers: Ch. 9, Eq. (47) printed 336; Ch. 12,
   Eqs. (40) and (93) printed 497 and 522; App. III Table VII(b) printed 769.
   Two faults follow. The weak limit becomes 1.049 sigma_R^2 instead of
   1.000 sigma_R^2, a 4.9 % error. And the module disagrees with itself,
   because `aperture_averaged_index_andrews` reduces at d = 0 to the correct
   0.49/1.11/0.51/0.69. The same constants sit in the shared kernel
   (`general_atmospherics.py:23`), which is the origin. One fix must cover both.
   — the olb half is fixed on andrews-foundation (PW-09: 0.49, 1.11 and 0.51
   are now in `plane_wave_scintillation_index_closed`; the constant 0.69 was
   already correct). The kernel half (KR-24) is STILL OPEN, because
   `my_analysis_modules` is a different repository.
2. **KR-18 — a misplaced parenthesis in the on-axis Dios integrand**
   (`coupled_flux.py:288`). The code computes
   `A^(5/6) [1 - (1+ratio^2)^(5/12)] cos(...)`; Dios Eq. (16) and Andrews
   Ch. 8, Eq. (17) need the cosine to multiply ONLY the second term. The error
   is non-zero for every finite Gaussian beam, and the olb uplink turbulence
   Term inherits it through `coupled_flux_sample`. The olb twin at
   `beam_wave_scintillation.py:137` is correct. — FIXED, and the fix is
   CONFIRMED against Dios Eq. (16), printed p. 3869, on 2026-08-25.
3. **KR-04 — the beam-wander variance constant is 2.07 against a book value of
   7.25.** The integrand is identical. R3 re-derived 7.25 from Ch. 6, Eqs. (88)
   and (89). This is a cross-source conflict (Dios against Andrews), not a
   plain bug. See Conflicts C-01. — CLOSED 2026-08-25: the kernel copies Dios
   Eq. (11) correctly, so 2.07 STAYS. Belmonte 2000 (Dios ref 23) prints the
   same 2.07 form (his Eq. (21)), traces it to the Yura / Mironov-Nosov image-
   motion level arm, and validates it against a split-step simulation (his
   Figs. 11 and 12) — a second confirmation after Dios Fig. 3. The Andrews 7.25
   is the beam-wave spectral-filter derivation; it over-counts by 3.50 against
   both simulations. See the C-01 closure block below.
4. **KR-20 — a second mean-irradiance weight on an already-normalised index**
   (`coupled_flux.py:379`). Andrews Ch. 8, Eqs. (9) and (15) normalise the
   scintillation index by the square of the LOCAL mean irradiance, so the extra
   `I_off^2` factor suppresses the index one more time. — NOT A FAULT, and the
   2026-08 removal is REVERSED. Dios Eq. (25), printed p. 3870, prints that
   weight, because Dios re-normalises the index to the beam CENTRE before his
   Eq. (26). The kernel implements Dios, so the weight is back.
5. **KR-01 and GF-18 — the spherical path weight is mirrored.** The code uses
   `((L-z)/L)^(5/3)`; Andrews Ch. 6, Eqs. (115) and (116) use `(z/L)^(5/3)`.
   The geometry resolves this: see Conflicts C-02. Do not "fix" it blind.
6. **UF-01 — RESOLVED: `WEAK_FLUCTUATION_LIMIT` was 0.6 (warned 2.4 times too
   late), now 0.25.**
   Confirmed by 2 readers (Ch. 8 text after Eq. (23) printed 264; Ch. 12,
   Eqs. (40) and (93) printed 497 and 522). The book limit sigma_R^2 < 1 gives
   sigma_x^2 < 0.25. The uplink budget can return numbers the book calls
   untrustworthy with no warning. — fixed on andrews-foundation. The limit is
   now 0.25, with the Ch. 8 (printed pp. 264-265) citation.
7. **PW-01 — `WEAK_FLUCTUATION_LIMIT = 0.25` is not a book number.** Confirmed
   by 4 readers. The book states sigma_R^2 < 1. The code errs in the safe
   direction, but the comment cites the book for a number the book does not
   give. See Conflicts C-05. — fixed on andrews-foundation. The value stays at
   0.25, but the comment now calls it a house rule that is stricter than the
   book, and gives the Ch. 11.3 (printed p. 451) lognormal-tail justification.

   > **STALE (UF-01, PW-01, and the constants-ledger rows below):** the name
   > `WEAK_FLUCTUATION_LIMIT` no longer exists in the code. The 0.25 house rule
   > has ONE canonical definition, `LOGNORMAL_PDF_LIMIT = 0.25` in
   > `andrews/scintillation.py`. The old file and line references
   > (`plane_wave_scintillation.py:45`, `uplink_flux.py:72`) are historical. The
   > findings above stay as the audit record; only the symbol name moved.
8. **TL-05 — the terrestrial weak gate uses a plane-wave threshold on a
   Gaussian beam.** Andrews Ch. 5, Eq. (16) printed 140 needs BOTH
   sigma_R^2 < 1 AND sigma_R^2 Lambda^(5/6) < 1. A focused or a strongly
   diffracted terrestrial beam can pass a gate it must fail.
9. **KR-05 — the long-term waist adds twice the wander variance.** Andrews
   Ch. 6, Eq. (100) printed 205 has the factor 1. See Conflicts C-03. — NOT A
   FAULT. Dios Eq. (1), printed p. 3867, prints the factor 2 on a RADIAL
   `<beta^2>` (Dios Eq. (10)). No code change. The net difference from Andrews
   is 1.75, not 3.50, because the factor 2 and the constant 2.07 partly cancel.

**Citation faults**

1. **PW-02 — a wrong equation number in four places.**
   `plane_wave_scintillation.py` cites "Ch. 12, Eq. (12.44)" in the module
   docstring (twice), the function docstring and `_scintillation_integral`.
   Ch. 12, Eq. (44) printed 498 is the downlink irradiance covariance. The
   correct number is Ch. 12, Eq. (38), printed 495 (repeated as Eq. (92),
   printed 522). The formula itself is exact. — fixed on andrews-foundation.
   All 4 places now cite Ch. 12, Eq. (38), printed p. 495.
2. **PW-12, PW-13 and PW-14 — three aperture-averaging factors credit "Andrews
   and Phillips Ch. 10" for Churnside constants.** The constants 1.07, 2.21,
   0.908 and 0.162 are nowhere in the 809-page book. They come from Churnside,
   Appl. Opt. 30 (1991) 1982, DOI 10.1364/AO.30.001982, which Andrews cites as
   Ch. 10 reference [12] but never reproduces. The book's own weak fit is
   Ch. 10, Eq. (61) printed 412: A = [1 + 1.062 k D_G^2/(4 L)]^(-7/6), with the
   7/6 OUTSIDE the bracket. The two fits differ by up to 12 %. — the
   ATTRIBUTION is fixed on andrews-foundation: all three docstrings now credit
   Churnside and name the different Andrews counterpart. The functional forms
   are unchanged, and a later work package supersedes them.
3. **AO-07 — `ao.py:151` credits "Andrews Ch. 3"** for the residual phase power
   spectral density 0.023 r0^(-5/3) f^(-11/3). That is Noll, JOSA 66 (1976)
   207, DOI 10.1364/JOSA.66.000207. The nearest Andrews statement, Ch. 14,
   Eq. (88) printed 635, is a different (von Karman, geometrical-optics) form.
4. **PW-01 — the 0.25 comment cites "Andrews Ch. 5"** for a threshold the book
   does not state. See physics fault 7.
5. **PW-05 — the docstring cites Ch. 12** for the aperture-averaged double
   integral. Ch. 9, Eq. (25) printed 333 and Ch. 10, Eq. (59) printed 412 are
   the closer citations.
6. **GF-06 and GF-07 — NOT a fault.** R3 called the "Ch. 9" docstring citation
   of the effective beam parameters a mis-citation, because the source is
   Ch. 7, Eq. (58). R5 then found Ch. 9, Eqs. (85) and (86) printed 349, which
   state the same result identically. The citation is valid, only less
   specific. Recorded so that nobody "fixes" it.

**Measurements**

1. **GAP 9 is now MEASURED and CLOSED.** The Andrews closed form Ch. 8, Eq. (23)
   printed 264 is now built, as
   `olb/turbulence/andrews/scintillation.py:beam_rytov_variance`. The self-check
   of that module compares it with the Dios path integral
   `beam_wave_scintillation.on_axis_scintillation_index` over one homogeneous
   horizontal path (lambda = 1550 nm, L = 2000 m, Cn2 = 3e-16, sigma_R^2 =
   0.0213, w0 = 5 cm). The measured difference is **+3.06 %** for a collimated
   beam and **-1.57 %** for a divergent beam with f0 = -1000 m. So the two
   models agree inside the 15 % gate. The residual is the book rounding of the
   constants 3.86 and 0.40, which gives 0.998 and 0.3996 in the plane-wave and
   spherical-wave limits.

## Table 1 — forward map

One row per equation. Sorted by olb file path, then by line number.

| olb id | symbol | location | quantity | book section | book equation | printed p | pdf p | reduction | status | note |
|---|---|---|---|---|---|---|---|---|---|---|
| BM-01 | w_v = lambda/(pi theta) | olb/beam.py:59 | virtual waist radius that gives a wanted far-field half-angle divergence | 4.4.1 | Ch. 4, Eq. (37) | 93 | 118 | collimated | reduction | The far-field limit of W = W0 sqrt(Theta0^2+Lambda0^2) with Theta0 = 1 and Lambda0 >> 1. The book states no explicit divergence formula in Ch. 4. |
| BM-02 | d = zR(w_v) sqrt((w0/w_v)^2-1) | olb/beam.py:60 | distance of the virtual waist behind the transmit aperture | 4.4.1 | Ch. 4, Eq. (37) | 93 | 118 | collimated | exact | Algebraic inversion of Eq. (37) with Theta0 = 1: W/W0 = sqrt(1+Lambda0^2), Lambda0 = d/zR. |
| BM-03 | zR = 0.5 k w^2 | olb/beam.py:60 (kernel `zR`) | Rayleigh range | 4.5.2 | Ch. 4, Sec. 4.5.2 text (the distance at which Lambda0 = 1) | 98 | 123 | collimated | exact | Andrews writes zR = 0.5 k W0^2 in running text, not as a numbered equation. |
| BM-04 | W(z) = w_v sqrt(1+(2z/(k w_v^2))^2) | olb/beam.py:86 (kernel `gaussz`) | free-space beam radius at range z | 4.4.1 | Ch. 4, Eq. (37); repeated as Eq. (139) | 93 (119) | 118 (144) | collimated | exact | The olb recast reaches a divergent beam through a shifted collimated beam, not through Theta0 > 1. |
| DL-01 | sigma_l^2 = ln(1+sigma2_P) | olb/links/downlink.py:64 | log-irradiance variance from the scintillation index | 9.2.3 | Ch. 9, Eq. (10) | 329 | 354 | weak | exact | Confirmed by 2 readers. Ch. 5, Eq. (95) printed 157 gives the same result as sigma_I^2 = exp(4 sigma_chi^2) - 1, with sigma_lnI^2 = 4 sigma_chi^2 from Eqs. (91)-(92). |
| DL-02 | mean_db = (5/ln10) sigma_l^2 | olb/links/downlink.py:68 | mean dB loss of a unit-mean lognormal | 9.11 | Ch. 9, Eq. (158) | 384 | 409 | weak | exact | Confirmed by 2 readers. A dB consequence of the lognormal PDF (Ch. 5, Eq. (93) printed 156), not a book constant. |
| DL-03 | quantile(p) = -(10/ln10)(-sigma_l^2/2 + sigma_l Phi_inv(1-p)) | olb/links/downlink.py:73 | fade depth not exceeded a fraction p of the time | 9.11 | Ch. 9, Eq. (158) | 384 | 409 | weak | exact | Confirmed by 2 readers. The book warns on the same page that the lognormal tail near the origin is too thin. That tail is the fade-critical region this face reports. See also Ch. 5, Eq. (93) and Sec. 5.10 printed 167. |
| DL-04 | mean-log offset = -sigma_l^2/2 | olb/links/downlink.py:73,76 | E[I] = 1 normalisation of the lognormal | 9.11 | Ch. 9, Eq. (158) | 384 | 409 | weak | exact | R2 marked this `unmatched` inside Ch. 4 and Ch. 5, because Ch. 5, Eq. (93) leaves the mean log-amplitude free. R5 then found the explicit -sigma_I^2/2 offset in Ch. 9, Eq. (158). Range-limited non-find, resolved. |
| DL-05 | `_gamma_gamma_term` | olb/links/downlink.py | gamma-gamma downlink scintillation Term, all fluctuation strengths | 9.10 | Ch. 9, Eqs. (137), (138), (139), (140); Ch. 12, Eq. (40) | 370-371 (497) | 395-396 (522) | point receiver | incorporated | Was `unmatched` (a reserved slot that raised). WP7 built it. The Term composes `andrews.scintillation.large_scale_log_variance` and `small_scale_log_variance` (Ch. 9, Eqs. (41) and (46)) with `andrews.distributions.gamma_gamma_params` (Eq. (138)), then turns the model into the three dB faces through `olb/models/fade.py`. Its index 1/alpha + 1/beta + 1/(alpha beta) (Eq. (139)) is identically the book weak-to-strong index of Ch. 12, Eq. (40). `model="auto"` selects it at sigma2_I >= 0.25. POINT receiver only: the book gives no aperture-averaged downlink index in this regime. See Table 2 rows G-102 to G-105. |
| RT-01 | top-hat correction | olb/links/retro_space.py:137 | Gaussian(waist = D/2) to uniform-aperture conversion | - | - | - | - | unobscured; on-axis | unmatched | Not in Ch. 12, Ch. 13.7 or Ch. 14. The source is `olb/models/gaussian_efficiency.py`, which no reader was assigned. |
| RT-02 | independent-turbulence assumption | olb/links/retro_space.py:169 | down-leg turbulence drawn independently of the up-leg | 13.7.4 | Ch. 13, Eqs. (147)-(159) | 581-584 | 606-609 | - | unmatched | Ch. 13, Eq. (146) printed 581 DOES support the module for spatial coherence: the reflected coherence radius equals a one-way spherical wave. But the monostatic scintillation index of 13.7.4 is not the sum of two independent legs. The module docstring already limits itself to a long slant path, so the gap is declared, not hidden. |
| RT-03 | P_return chain | olb/links/retro_space.py:190 | retroreflected return as an up-leg plus a down-leg dB sum | 13.7 | Ch. 13, Eqs. (132)-(142) | 577-580 | 602-605 | - | unmatched | The book models the double passage as one coupled problem and gives a backscatter amplification factor, which olb has no Term for. See Table 2 rows G-152 to G-156. |
| TL-01 | sigma_l^2 = ln(1+sigma2_P) | olb/links/terrestrial.py:152 | log-irradiance variance from the scintillation index | 5.7.2 | Ch. 5, Eq. (95), inverted | 157 | 182 | weak | exact | Duplicate of DL-01 (olb gap 10). |
| TL-02 | mean_db = (5/ln10) sigma_l^2 | olb/links/terrestrial.py:155 | mean dB loss of a unit-mean lognormal | 5.7.2 | Ch. 5, Eq. (93) | 156 | 181 | weak | exact | Duplicate of DL-02 (olb gap 10). |
| TL-03 | quantile(p) | olb/links/terrestrial.py:160 | fade depth not exceeded a fraction p of the time | 5.7.2 | Ch. 5, Eq. (93) | 156 | 181 | weak | exact | Duplicate of DL-03 (olb gap 10). |
| TL-04 | mean-log offset = -sigma_l^2/2 | olb/links/terrestrial.py:160,163 | E[I] = 1 normalisation | 9.11 | Ch. 9, Eq. (158) | 384 | 409 | weak | exact | Duplicate of DL-04 (olb gap 10). |
| TL-05 | validity test sigma2_I < 0.25 on a Gaussian beam | olb/links/terrestrial.py | weak-fluctuation gate for a beam wave | 5.2.2 | Ch. 5, Eq. (16) | 140 | 165 | on-axis | FIXED 2026-08-29 | Andrews printed 140 says the plane-wave criterion is not adequate for a Gaussian beam, and needs BOTH sigma_R^2 < 1 AND sigma_R^2 Lambda^(5/6) < 1, with Lambda = 2L/(k W^2). RESOLVED: the terrestrial Term now calls the shared `rytov_weak(sigma2_R, Lambda)` (olb/turbulence/andrews/scintillation.py), which applies both conditions via the binding strength sigma2_R * max(1, Lambda^(5/6)), so a focused beam is caught. The lognormal-PDF house rule is now the distinct LOGNORMAL_PDF_LIMIT on sigma2_I. |
| PT-01 | h(r) = exp(-2 r^2/w_z^2) | olb/models/pointing.py:10 (docstring), used at :54 | on-axis Gaussian roll-off of collected power against boresight | 11.4.1 | Ch. 11, Eqs. (40) and (41) | 459 | 484 | on-axis; unobscured | exact | The book writes the same free-space Gaussian profile in the detector plane. |
| PT-02 | E[loss] = (20/ln10) 2 sigma_r^2/w_z^2 | olb/models/pointing.py:54 | mean pointing-jitter loss in dB | - | - | - | - | on-axis; unobscured | unmatched | Ch. 11.3 states at printed 451 that it neglects any possible pointing error. Ch. 12 DOES treat pointing error, but through a beam-wander pointing variance folded into the off-axis scintillation index (Ch. 12, Eqs. (53) and (54), printed 503), not as an exponential dB fade. See Table 2 rows G-140 to G-142. |
| PT-03 | loss_db quantile = -mean ln(1-p) | olb/models/pointing.py:101 | inverse exponential CDF of the pointing loss | - | - | - | - | on-axis | unmatched | Ch. 11.3.1 builds the fade quantile from a lognormal (Eq. (24)) or a gamma-gamma (Eq. (26)) irradiance PDF, printed 451-452. No exponential-in-dB fade law appears in Ch. 10, Ch. 11 or Ch. 12. |
| PT-04 | loss_db ~ Exponential(mean) sampler | olb/models/pointing.py:104 | Monte Carlo draw of the pointing loss | - | - | - | - | on-axis | unmatched | Same as PT-03. |
| RS-01 | q(p) = mean_db for a deterministic term | olb/results.py:62 | per-term quantile face | 11.3.1 | Ch. 11, Eq. (23) | 451 | 476 | none | exact | A degenerate PDF integrates to a step, which is consistent with Eq. (23). |
| RS-02 | total_loss_db = SUM of mean_db | olb/results.py:130 | deterministic budget total in dB | - | - | - | - | none | unmatched | Ch. 10 and Ch. 11 never add dB losses. They carry the linear mean irradiance and the flux variance separately. Ch. 11, Eq. (28) printed 452 is the only "total" the book forms, and it is a variance. |
| RS-03 | total fade = SUM over terms of q_t(p) | olb/results.py:223 (sum at :249-258) | analytic fade margin of the whole budget | 11.3.1 | Ch. 11, Eq. (23) | 451 | 476 | none | approximate | The book computes the fade from ONE joint irradiance PDF. A sum of per-term p-quantiles is not a book form. It is a conservative upper bound, as the olb docstring itself states. |
| RS-04 | fade level at availability p | olb/results.py:296 (monte_carlo percentile) | joint loss level not exceeded a fraction p of the time | 11.3.1 | Ch. 11, Eqs. (23) and (25) | 451 | 476 | none | exact | Numerical inverse of Pr(I <= I_T). The book fade threshold F_T = 10 log10(mean I / I_T) [dB] is the same dB quantity olb reports. |
| AA-01 | sigma2_theta = <r_c^2>/L^2 | olb/turbulence/angle_of_arrival.py:60 | beam-wander arrival tilt | 6.6.2 | Ch. 6, Eq. (94) and the text below it | 204 | 229 | Kolmogorov; l0=0; L0=inf | exact | Book: "<r_c^2> = L^2 <beta_a^2> by equating W_G = W0". The MAPPING is exact. The INPUT <r_c^2> is not: it comes from the kernel row KR-04, which is `wrong`. |
| AA-02 | <beta_a^2> | olb/turbulence/angle_of_arrival.py:75 | aperture angle-of-arrival tilt variance | 6.5 | Ch. 6, Eq. (84); definition Eq. (82); scales Eq. (83) | 201 (200) | 226 (225) | Kolmogorov; l0=0; L0=inf; plane wave | incorporated | Was `unmatched` (a deferred stub that raised). The owner decided C-04 on the GRADIENT tilt, and the stub now delegates to `andrews.structure.angle_of_arrival_variance`, which also carries the inner-scale and outer-scale branches of Eq. (83). THE BOOK FORM: <beta_a^2> = 2.91 Cn2 L (2 W_G)^(-1/3) for 2 W_G >> l0. It is ONE AXIS. With D = 2 W_G and r0 = (0.423 k^2 Cn2 L)^(-3/5) it becomes 0.174 (D/r0)^(5/3)(lambda/D)^2, NOT the 0.182 named in batch 2. The slant-path version is Ch. 12, Eq. (28) printed 492. See Conflicts C-04. |
| AN-01 | p_n(u) = 4 (n+1)^2 [J_{n+1}(u)/u]^2 | olb/turbulence/anisoplanatism.py:117 | Zernike radial-order weight | 14.5.3 | Ch. 14, Eq. (86) | 634 | 659 | unobscured; plane wave | exact | STRONG INDEPENDENT CONFIRMATION. Summing the Andrews per-mode filter over the n+1 azimuthal modes of radial order n and averaging over the azimuth gives exactly the Stone p_n. |
| AN-02 | M(u) = 1 - p_0(u) | olb/turbulence/anisoplanatism.py:153 | piston-removed modal weight | 14.5.4 | Ch. 14, Eq. (89) | 635 | 660 | unobscured; plane wave | exact | The Andrews piston-removed integrand carries {1 - [2 J1(kappa D/2)/(kappa D/2)]^2}. |
| AN-03 | (n+1)(n+2)/2 | olb/turbulence/anisoplanatism.py:225 | Noll mode count through radial order n | 14.5.2 | Ch. 14, Table 14.1 | 630 | 655 | - | exact | The counting rule implied by the Noll ordering that Andrews tabulates through (m, n) = (3, 3). |
| AN-04 | theta0, C1 = 2.914381 | olb/turbulence/anisoplanatism.py:268 | classical isoplanatic angle | 12.4.3 | Ch. 12, Eq. (30) | 493 | 518 | plane wave; Kolmogorov; L0=inf | exact | olb multiplies the integral by airmass^(8/3) then raises it to -3/5, which gives cos^(8/5)(zeta). Identical to the book. Repeated as Ch. 14, Eqs. (37) and (97). The Stone 2.914381 is the unrounded 2.91. |
| AN-05 | sigma^2 = 2 (2 pi)^(8/3) C_A k^2 R^(5/3) INT Cn2 I(beta) dh | olb/turbulence/anisoplanatism.py:358 | finite-aperture anisoplanatic phase variance | - | - | - | - | plane wave; Kolmogorov; L0=inf | other source | Stone et al., JOSA A 11 (1994) 347, DOI 10.1364/JOSAA.11.000347, Eqs. (29) and (36). Andrews has NO finite-aperture anisoplanatism: Ch. 12.4.3 printed 493 and Ch. 14.3.4 printed 622 give only the zero-aperture law. This is a genuine olb capability beyond the book. |
| AN-06 | sigma^2 = (theta/theta0)^(5/3) | olb/turbulence/anisoplanatism.py:392 | classical anisoplanatic phase variance | 14.3.4 | Ch. 14, Eq. (36) | 622 | 647 | plane wave; Kolmogorov; L0=inf | exact | Andrews writes it as the MTF exp[-(theta/theta0)^(5/3)]; the exponent is the variance. |
| AO-01 | NOLL_PISTON = 1.0299 | olb/turbulence/ao.py:39, applied at :178 | piston-removed phase variance, c (D/r0)^(5/3) | 14.5.4 | Ch. 14, Eq. (90) | 635 | 660 | plane wave; Kolmogorov; L0=inf; unobscured | exact | Andrews (after Sasiela) gives 1.03; the Andrews ABCD route, Eq. (91) on the same page, gives 1.02. The Noll 1.0299 sits inside that pair. |
| AO-02 | NOLL_TIPTILT = 0.134 | olb/turbulence/ao.py:40, applied at :178 | piston-and-tilt-removed phase variance | 14.5.4 | Ch. 14, Eq. (94) | 636 | 661 | plane wave; Kolmogorov; L0=inf; unobscured | exact | Andrews gives 0.13, formed as Eq. (90) 1.03 minus Eq. (93) 0.90. The Andrews tilt-only value 0.90 is NOT in olb. See Table 2 row G-165. |
| AO-03 | 0.2944 J^(-sqrt(3)/2) | olb/turbulence/ao.py:41, applied at :97 | large-J adaptive-optics residual coefficient | - | - | - | - | plane wave; Kolmogorov; L0=inf; unobscured | other source | Noll, JOSA 66 (1976) 207, DOI 10.1364/JOSA.66.000207. Ch. 14.5.4 stops at piston and tilt and gives no large-J asymptote. Correctly cited in the code. |
| AO-04 | sec(zeta) | olb/turbulence/ao.py:67 | airmass 1/sin(elevation) on the Cn2 integral | 12.4.1 | Ch. 12, Eq. (23) | 492 | 517 | flat atmosphere | exact | R3 marked this `unmatched` inside Ch. 6, because Sec. 6.8 says only that Cn2 becomes Cn2(h) and gives no slant scale factor. R7 resolved it: Ch. 12, Eq. (23) carries the sec(zeta) explicitly. |
| AO-05 | r0 = (0.423 k^2 sec(zeta) INT Cn2 dh)^(-3/5) | olb/turbulence/ao.py:69 | plane-wave Fried parameter | 12.4.1 | Ch. 12, Eq. (23) | 492 | 517 | plane wave; Kolmogorov; l0=0; L0=inf | exact | Confirmed by 2 readers. Also Ch. 6, Eq. (113) printed 209 and Ch. 14, Eqs. (25) printed 617 and (95) printed 636. The book prints 0.42; 0.423 is the unrounded 1.46/2.1^(5/3) = 0.4240. |
| AO-06 | f_c = sqrt(n_modes)/(2 D) | olb/turbulence/ao.py:147 | adaptive-optics correction cutoff | - | - | - | - | unobscured | unmatched | No counterpart in Ch. 12.4, Ch. 14.3 or Ch. 14.5. It is a heuristic. |
| AO-07 | Phi_phi(f) = 0.023 r0^(-5/3) f^(-11/3) | olb/turbulence/ao.py:151 | residual phase power spectral density | - | - | - | - | Kolmogorov; L0=inf | other source | Noll 1976, DOI 10.1364/JOSA.66.000207. CITATION FAULT: the docstring credits "Andrews Ch. 3". The nearest Andrews statement, Ch. 14, Eq. (88) printed 635, is a different (von Karman, geometrical-optics) form 1.83 (L0/r0)^(5/3). |
| BW-01 | `_KOLMOGOROV = 0.033`, Phi_n = 0.033 Cn2 kappa^(-11/3) | olb/turbulence/beam_wave_scintillation.py:55 | Kolmogorov refractive-index spectrum carried by the Dios integrals | 3.3.1 | Ch. 3, Eq. (18) | 67 | 92 | Kolmogorov; l0=0; L0=inf | exact | The code integrates a Cn2(h) grid, so the exact form used is the path-varying Eq. (26) printed 72. Restated as Ch. 6, Eq. (110) printed 208, Ch. 8, Eq. (120) printed 299 and Ch. 12, Eq. (15) printed 490. |
| BW-02 | Theta0 = 1 - L/F0 | olb/turbulence/beam_wave_scintillation.py:84 | input-plane curvature parameter | 8.1 | Ch. 8, Eq. (5) | 261 | 286 | - | exact | Confirmed by 2 readers. Also App. III Table VI footer printed 768 and Ch. 12, Eq. (8) printed 488. The code default f0 = inf gives Theta0 = 1, but the argument is threaded, so this call site is NOT collimated-only. |
| BW-03 | Lambda0 = 2 L/(k W0^2) | olb/turbulence/beam_wave_scintillation.py:85 | input-plane Gaussian beam parameter | 8.1 | Ch. 8, Eq. (5) | 261 | 286 | - | exact | Confirmed by 2 readers. Also App. III Table VI footer printed 768. |
| BW-04 | Theta = Theta0/(Theta0^2+Lambda0^2) | olb/turbulence/beam_wave_scintillation.py:87 | output-plane curvature parameter | 8.1 | Ch. 8, Eq. (6) | 261 | 286 | - | exact | Confirmed by 2 readers. |
| BW-05 | Lambda = Lambda0/(Theta0^2+Lambda0^2) | olb/turbulence/beam_wave_scintillation.py:88 | output-plane beam parameter | 8.1 | Ch. 8, Eq. (6) | 261 | 286 | - | exact | Confirmed by 2 readers. |
| BW-06 | W^2(L) = 2 L/(k Lambda) | olb/turbulence/beam_wave_scintillation.py:89 | receiver-plane beam radius squared | 8.1 | Ch. 8, Eq. (6) | 261 | 286 | - | exact | The code inverts the table definition Lambda = 2L/(k W^2). |
| BW-07 | sec(zeta) over an altitude grid | olb/turbulence/beam_wave_scintillation.py:91-92, 140 | slant-path weighting of Cn2, with the (L-z)/L transmitter weight | 8.7.1 | Ch. 8, Eqs. (115) and (118) | 299 | 324 | Kolmogorov; l0=0; L0=inf; weak | reduction | The book integrates Cn2(z) along the path element dz. The code integrates over altitude dh and multiplies by sec(zeta). Equivalent only for a straight slant path with no refraction. |
| BW-08 | A(z) = (Lambda L/k)((L-z)/L)^2 | olb/turbulence/beam_wave_scintillation.py:94 | Gaussian attenuation scale in the spectral integrand | 8.2 | Ch. 8, Eqs. (14) and (17) exponent | 262-263 | 287-288 | Kolmogorov; l0=0; L0=inf; weak | other source | Dios et al., DOI 10.1364/AO.43.003866, Eq. (17). It is the Andrews exponential kernel with xi = 1 - z/L after the kappa integration. The Andrews closed-form twin is App. III Table IX(a) printed 771. |
| BW-09 | B(z) = (L/k)((L-z)/L)(Theta + (1-Theta) z/L) | olb/turbulence/beam_wave_scintillation.py:95 | phase term of the spectral integrand | 8.2 | Ch. 8, Eq. (17) cosine argument | 263 | 288 | Kolmogorov; l0=0; L0=inf; weak | other source | Dios Eq. (18), DOI 10.1364/AO.43.003866. Same phase term as the Andrews longitudinal integrand. |
| BW-10 | on_axis_scintillation_index | olb/turbulence/beam_wave_scintillation.py:137-141 | longitudinal (on-axis) scintillation index | 8.2 | Ch. 8, Eq. (17); closed form Eq. (19); slant form Eq. (118) | 263 (299) | 288 (324) | Kolmogorov; l0=0; L0=inf; weak | other source | Dios Eq. (16), DOI 10.1364/AO.43.003866. It is the Andrews Eq. (17) integral written over z with a variable Cn2. The module self-check reproduces the Andrews Eq. (20) plane and spherical limits to 2-3 %. The Andrews closed form (App. III Table IX(a) printed 771, and Ch. 8, Eq. (23) printed 264) is NOT compared anywhere in olb: olb gap 9. |
| BW-11 | pred = 4.42 sigma_R^2 Lambda^(5/6)(r/W)^2 | olb/turbulence/beam_wave_scintillation.py:163, 231 | small-r radial approximation | 8.2 | Ch. 8, Eq. (22) | 264 | 289 | Kolmogorov; l0=0; L0=inf; weak | exact | The book states r < W. |
| BW-12 | radial_scintillation_index | olb/turbulence/beam_wave_scintillation.py:169-172 | radial (off-axis) component | 8.2 | Ch. 8, Eq. (16); closed form Eq. (18) | 263 | 288 | Kolmogorov; l0=0; L0=inf; weak | other source | Dios Eq. (20). Andrews Eq. (18) reads 2.64 sigma_R^2 Lambda^(5/6)[1 - 1F1(-5/6;1;2r^2/W^2)]. The code has (1F1 - 1) times Gamma(-5/6) < 0, so the sign agrees. The Summary prints the constant as 2.65 (Eq. (129), printed 302). |
| BW-13 | gaussian_scintillation_index = longitudinal + radial | olb/turbulence/beam_wave_scintillation.py:192 | total Gaussian-beam scintillation index | 8.2 | Ch. 8, Eq. (15); closed form Eq. (21); slant form Eq. (117) | 263-264 (299) | 288-289 (324) | Kolmogorov; l0=0; L0=inf; weak | exact | Dios Eq. (13) is the Andrews Eq. (15) split. |
| BW-14 | sigma2_R = 1.23 Cn2 k^(7/6) L^(11/6) | olb/turbulence/beam_wave_scintillation.py:208 | plane-wave Rytov variance | 8.2 | Ch. 8, Eq. (20) | 264 | 289 | Kolmogorov; l0=0; L0=inf; weak; plane wave; homogeneous Cn2 | exact | |
| BW-15 | 0.404 sigma2_R | olb/turbulence/beam_wave_scintillation.py:216 | spherical-wave limit beta0^2 | 8.2 | Ch. 8, Eq. (20) | 264 | 289 | Kolmogorov; l0=0; L0=inf; weak; homogeneous Cn2 | exact | The book gives 0.4 exactly; the code asserts 0.404 with rtol 3e-2. Ch. 9, Eqs. (63) and (64) printed 341 give the same pair. |
| BW-16 | sigma2_I approximately 4 sigma2_chi | olb/turbulence/beam_wave_scintillation.py:242, 261 | link between the index and the log-amplitude variance | 8.2 | Ch. 8, Eq. (13) | 262 | 287 | weak | exact | Book: sigma_I^2 = exp(4 sigma_x^2) - 1, which is approximately 4 sigma_x^2. |
| BW-17 | Hufnagel-Valley Cn2(h) | olb/turbulence/beam_wave_scintillation.py:246-249 | Cn2 profile used in the self-check | 12.2.1 | Ch. 12, Eq. (1) | 481 | 506 | - | exact | R4 marked this `unmatched` inside Ch. 8, which has no Cn2 profile. R7 resolved it: every constant, exponent and scale height matches Ch. 12, Eq. (1). The defaults w = 21 m/s and A = 1.7e-14 are the book H-V5/7. |
| GF-01 | COLLIMATED_THETA0 = 1 | olb/turbulence/gaussian_fried.py:32 | input-plane curvature parameter, fixed at 1 | 6.2.1 | Ch. 6, Eq. (6) | 183 | 208 | collimated | reduction | Confirmed by 3 readers. Also Ch. 4, Eq. (33) printed 92 and the Ch. 9 worked example printed 384. Andrews keeps Theta0 = 1 - L/F0 general: Theta0 = 1 collimated, less than 1 convergent, more than 1 divergent. The profile form at line 310 already carries the general 1 - L/f0, and the terrestrial call site now passes f0 (Gap 3 wired, 2026-08-27). Only this single-path form keeps the fixed Theta0, and no budget calls it. |
| GF-02 | Lambda0 = 2 z/(k w0^2) | olb/turbulence/gaussian_fried.py:49 | input-plane Fresnel ratio | 6.2.1 | Ch. 6, Eq. (6) | 183 | 208 | - | exact | Confirmed by 3 readers. Also Ch. 4, Eqs. (33) and (136) printed 92 and 118. |
| GF-03 | Lambda = Lambda0/(Theta0^2+Lambda0^2) | olb/turbulence/gaussian_fried.py:63 | output-plane diffraction parameter | 6.2.1 | Ch. 6, Eq. (7) | 183 | 208 | collimated | reduction | Confirmed by 3 readers. The form is exact; only the hard-wired Theta0 = 1 reduces it. Also Ch. 4, Eqs. (44) and (138) printed 95. |
| GF-04 | Theta = Theta0/(Theta0^2+Lambda0^2) | olb/turbulence/gaussian_fried.py:64 | output-plane refraction parameter | 6.2.1 | Ch. 6, Eq. (7) | 183 | 208 | collimated | reduction | Confirmed by 3 readers. Andrews Ch. 4, Eq. (45) adds the complement Theta_bar = 1 - Theta, which olb computes only in the profile form. |
| GF-05 | sigma_R = (1.23 Cn2 k^(7/6) z^(11/6))^0.5 | olb/turbulence/gaussian_fried.py:78 | plane-wave Rytov standard deviation | 6.9 | Ch. 6, Eq. (119) | 210 | 235 | Kolmogorov; l0=0; L0=inf; plane wave; homogeneous Cn2 | exact | Confirmed by 3 readers. Also Ch. 5, Eq. (15) printed 140, Ch. 7, Eq. (1) printed 230 and Ch. 9 Sec. 9.2 printed 323. Duplicate of PW-07 and KR-23: olb gap 10. |
| GF-06 | Theta_e = (Theta - 0.81 s^(12/5) Lambda)/(1 + 1.63 s^(12/5) Lambda) | olb/turbulence/gaussian_fried.py:98 | strong-turbulence effective curvature | 7.4 | Ch. 7, Eq. (58) | 242 | 267 | Kolmogorov | exact | Confirmed by 2 readers. The book writes (Theta + 2 q Lambda/3)/(1 + 4 q Lambda/3) with q = 1.22 s^(12/5), so 2q/3 = 0.813 and 4q/3 = 1.627. Ch. 9, Eq. (85) printed 349 restates it identically, so the olb docstring citation of "Ch. 9" is VALID, not a fault. See Conflicts C-08. |
| GF-07 | Lambda_e = Lambda/(1 + 1.63 s^(12/5) Lambda) | olb/turbulence/gaussian_fried.py:99 | strong-turbulence effective Fresnel ratio | 7.4 | Ch. 7, Eq. (58) | 242 | 267 | Kolmogorov | exact | Confirmed by 2 readers. Ch. 9, Eq. (86) printed 349, restated as Eq. (150) printed 382, gives Lambda_e = 2L/(k W_LT^2). |
| GF-08 | a = (1 -/+ abs(Theta)^(8/3))/(1 - Theta) | olb/turbulence/gaussian_fried.py:114 | a-factor of the Gaussian-beam structure function | 7.4.1 | Ch. 7, Eq. (60) | 243 | 268 | - | exact | olb feeds Theta_e, so Ch. 7, Eq. (60) is the exact match. The free-space form is Ch. 6, Eq. (55) printed 192. R5 marked it `unmatched` because Ch. 9 only refers to Sec. 7.4.1; range-limited non-find, resolved. |
| GF-09 | rho0_e = (8/(3(a_e + 0.62 Lambda_e^(11/6))))^(3/5) | olb/turbulence/gaussian_fried.py:129 | beam coherence ratio | 7.4.1 | Ch. 7, Eq. (59) lower | 243 | 268 | Kolmogorov; l0=0; L0=inf | exact | Ch. 6, Eq. (79) printed 199 gives the same 0.62; the Ch. 6 summary Eq. (132) printed 213 prints 0.618. R3 verified it numerically. R5 range-limited non-find, resolved. |
| GF-10 | rho_pl = (1.46 Cn2 k^2 z)^(-3/5) | olb/turbulence/gaussian_fried.py:141 | plane-wave coherence radius | 6.4.1 | Ch. 6, Eq. (64) lower; Eq. (130) | 194 (213) | 219 (238) | Kolmogorov; l0=0; L0=inf; plane wave; homogeneous Cn2 | exact | Duplicate of PW-08: olb gap 10. App. III Table IV printed 767 is the direct table, and it also gives the inner-scale branches 1.64 and 1.87 that olb lacks. |
| GF-11 | r0 = 2.1 rho0 | olb/turbulence/gaussian_fried.py:155 | plane-wave Fried parameter | 6.4.1 | Ch. 6, text below Eq. (64) | 194 | 219 | - | exact | Book: "atmospheric coherence width r0 = 2.1 rho0". Also the App. III note under Table IV printed 767. |
| GF-12 | (8/3)^(3/5) | olb/turbulence/gaussian_fried.py:163 | spherical-over-plane Fried ratio | 6.4.2 | Ch. 6, Eq. (71) | 196 | 221 | homogeneous Cn2 | approximate | STATUS DISAGREEMENT, see Conflicts C-07. R3: the book gives rho_sp = (0.55 Cn2 k^2 L)^(-3/5); the exact 3/8 weight gives 0.5475, so 0.55 is rounded, and olb 1.7963 differs from the book 1.7913 by 0.3 %. R5: the Ch. 9 worked example printed 384 gives r0 = (0.16 Cn2 k^2 L)^(-3/5), and 0.423 x 3/8 = 0.1586, which the book rounds to 0.16, so the olb constant is confirmed. |
| GF-13 | r0_sph | olb/turbulence/gaussian_fried.py:180 | spherical-wave Fried parameter | 6.4.2 | Ch. 6, Eq. (71) | 196 | 221 | homogeneous Cn2; Kolmogorov; l0=0; L0=inf | reduction | |
| GF-14 | r0_gauss = 2.1 rho0_e rho_pl | olb/turbulence/gaussian_fried.py:207 | Gaussian-beam Fried parameter, single path | 7.4.1 | Ch. 7, Eq. (59) lower, times 2.1 | 243 | 268 | collimated; Kolmogorov; l0=0; L0=inf; homogeneous Cn2 | reduction | The book Eq. (59) gives the 1/e radius. olb applies the 2.1 Fried factor of Ch. 6 printed 194. |
| GF-15 | Theta0 = 1 - L/f0 (profile form) | olb/turbulence/gaussian_fried.py:310 | input curvature, general | 6.2.1 | Ch. 6, Eq. (6) | 183 | 208 | - | exact | The profile form is already general in f0. The single-path form (GF-01) is not. |
| GF-16 | Lambda0 = 2L/(k w0^2) (profile form) | olb/turbulence/gaussian_fried.py:311 | input Fresnel ratio | 6.2.1 | Ch. 6, Eq. (6) | 183 | 208 | - | exact | |
| GF-17 | free-space Theta, Lambda in the profile form | olb/turbulence/gaussian_fried.py:313 | output beam parameters | 9.6.1 | Ch. 9, Eqs. (85) and (86) | 349 | 374 | weak | reduction | Confirmed by 2 readers. olb uses the free-space Theta and Lambda, not Theta_e and Lambda_e. The replacement is coded at lines 98-99 (GF-06, GF-07) but NO caller uses it. This is the deliberate deferral marked "ponytail" at line 308, and olb gap 4. |
| GF-18 | mu1 = (Theta + Theta_bar xi)^(5/3), xi = (L-z)/L | olb/turbulence/gaussian_fried.py:319 | wave-structure-function path weight | 6.8.2 | Ch. 6, Eq. (115) | 209 | 234 | Kolmogorov; l0=0; L0=inf; weak | wrong | The book weight is (Theta + Theta_bar z/L)^(5/3) with z measured FROM THE TRANSMITTER, so the spherical limit is (z/L)^(5/3) and the Cn2 near the RECEIVER carries the weight (Eq. (116) and the text below it, printed 209). olb uses the MIRROR weight, which matches Dios et al., Appl. Opt. 43 (2004) 3866, Eq. (3), the transmitter-plane (reciprocal) coherence radius of an uplink beam. See Conflicts C-02. Do not "fix" it blind. |
| GF-19 | 0.62 Lambda^(11/6) INT Cn2 xi^(5/3) | olb/turbulence/gaussian_fried.py:320,323 | mu2 term of the profile coherence radius | 6.8.2 | Ch. 6, Eq. (115) | 209 | 234 | Kolmogorov; l0=0; L0=inf | exact | The coefficient 0.62 agrees. The same mirror-weight caution as GF-18 applies. |
| GF-20 | rho0 = (1.46 k^2 sec(zeta)(mu1 + 0.62 Lambda^(11/6) mu2))^(-3/5) | olb/turbulence/gaussian_fried.py:322 | profile coherence radius | 6.8.2 | Ch. 6, Eq. (115) | 209 | 234 | Kolmogorov; l0=0; L0=inf; weak | reduction | The book Eq. (115) has no sec(zeta); olb adds the airmass. The airmass itself is book-supported: Ch. 12, Eq. (23) printed 492. See AO-04. |
| GF-21 | r0 = 2.1 rho0 (profile form) | olb/turbulence/gaussian_fried.py:324 | profile Fried parameter | 6.4.1 | Ch. 6, text below Eq. (64) | 194 | 219 | - | exact | |
| PW-01 | WEAK_FLUCTUATION_LIMIT = 0.25 | olb/turbulence/plane_wave_scintillation.py:45 | weak-fluctuation validity gate | 5.2.2 | Ch. 5, Eq. (15) and the text after it | 140 | 165 | weak | wrong | Confirmed by 4 readers that 0.25 is NOT a book number. The comment cites "Andrews Ch. 5, sigma2_R less than about 0.25". Andrews printed 140 says weak is sigma_R^2 < 1 and moderate is sigma_R^2 about 1; Ch. 10, Eq. (61) printed 412 and Ch. 12, Eq. (40) printed 497 say the same. R6 argues 0.25 is a defensible house rule, because Ch. 11.3 printed 451 says the lognormal tail is optimistic. The code errs safe; the CITATION is false. See Conflicts C-05. |
| PW-02 | sigma2_I = 2.25 k^(7/6) sec^(11/6) INT Cn2 h^(5/6) dh | olb/turbulence/plane_wave_scintillation.py:80 | slant-path plane-wave point scintillation index | 12.5.2 | Ch. 12, Eq. (38) | 495 | 520 | plane wave; Kolmogorov; l0=0; L0=inf; weak | exact | CITATION FAULT. The code says "Eq. (12.44)" in 4 places (module docstring twice, function docstring, `_scintillation_integral` docstring). Ch. 12, Eq. (44) printed 498 is the downlink irradiance covariance. The correct number is Ch. 12, Eq. (38) printed 495, repeated as Eq. (92) printed 522. The formula agrees exactly. |
| PW-03 | [2 J1(x)/x]^2, x = kappa D/2 | olb/turbulence/plane_wave_scintillation.py:108 | circular-aperture averaging filter | 10.3.2 | Ch. 10, Eq. (59); Ch. 10, Eq. (54) | 412 (410) | 437 (435) | unobscured | approximate | STATUS DISAGREEMENT, see Conflicts C-06. R6: Andrews Eq. (59) uses the SOFT Gaussian aperture exp(-D_G^2 kappa^2/16) with D_G^2 = 8 W_G^2; the Airy form is the Fourier transform of the hard circular MTF in the brackets of Eq. (54); the two agree in the limits but not in between. R7: the identical function IS printed, as the piston Zernike filter, Ch. 14, Eq. (86) with m = n = 0, printed 634; Andrews uses it in the other role, as {1 - [2J1/x]^2} in Eq. (89). The obscured (annular) aperture is in neither. |
| PW-04 | Phi_n = 0.033 kappa^(-11/3) | olb/turbulence/plane_wave_scintillation.py:138 | Kolmogorov refractive-index spectrum, Cn2 factored out | 9.2.2 | Ch. 9, Eq. (3) with G = 1 | 327 | 352 | Kolmogorov; l0=0; L0=inf | reduction | Confirmed by 3 readers. Eq. (3) is Phi_n,e = 0.033 Cn2 kappa^(-11/3) G(kappa,l0,L0); the olb code sets G = 1, so it drops the inner-scale and outer-scale filters. The bare form is Ch. 3, Eq. (18) printed 67 and Ch. 12, Eq. (15) printed 490. Ch. 3 warns that extending Eq. (18) to all kappa can make integrals diverge, and the olb kappa grid has no cut. |
| PW-05 | 8 pi^2 k^2 sec INT INT kappa Phi (1 - cos(kappa^2 z/k)) F(kappa) | olb/turbulence/plane_wave_scintillation.py:149 | aperture-averaged flux scintillation index, weak fluctuation | 10.3.2 | Ch. 10, Eq. (59) | 412 | 437 | plane wave; Kolmogorov; l0=0; L0=inf; weak; unobscured | reduction | Confirmed by 3 readers. The same double integral. olb replaces the normalised path variable with a Cn2(h)-weighted height integral, which is the standard slant-path generalisation; the aperture filter differs (see PW-03). The point form is Ch. 9, Eq. (25) printed 333; the slant form is Ch. 12, Eqs. (16) printed 491 and (75) printed 514. CITATION: the docstring cites Ch. 12; Ch. 9, Eq. (25) is the closer match. Andrews computes the same quantity in closed form for a hard aperture by the ABCD route in Ch. 12, Eq. (39) printed 496; the two are not compared in olb. |
| PW-06 | A = sigma2_I(D)/sigma2_I(0) | olb/turbulence/plane_wave_scintillation.py:185 | aperture-averaging factor, definition | 10.3 | Ch. 10, Eq. (56) | 410 | 435 | unobscured; plane wave | exact | The book defines A from the normalised pupil-plane covariance; olb takes the ratio of the same wavenumber integral with and without the filter. Equivalent definition. Ch. 12, Eq. (39) printed 496 states that the D_G = 0 limit reduces to the point index, which is the consistency this ratio relies on. |
| PW-07 | sigma_1 = (1.23 Cn2 k^(7/6) L^(11/6))^0.5 | olb/turbulence/plane_wave_scintillation.py:251 | plane-wave Rytov standard deviation, single path | 5.2.2 | Ch. 5, Eq. (15) | 140 | 165 | Kolmogorov; plane wave; weak; homogeneous Cn2 | exact | Confirmed by 6 readers. Also Ch. 6, Eq. (119) printed 210; Ch. 8, Eq. (20) printed 264; Ch. 9 Sec. 9.2 printed 323; Ch. 10, text at Eq. (60) printed 412; App. III Table VII(a) footer printed 769. Ch. 13 writes the same quantity as sigma_1^2. Triplicated in olb: olb gap 10. |
| PW-08 | rho_c = (1.46 Cn2 k^2 L)^(-3/5) | olb/turbulence/plane_wave_scintillation.py:264 | plane-wave spatial coherence radius | 6.4.1 | Ch. 6, Eq. (64) lower; Eq. (130) | 194 (213) | 219 (238) | Kolmogorov; plane wave; l0=0; L0=inf; homogeneous Cn2 | exact | Confirmed by 2 readers. The code cites "Ch. 6", which is right. App. III Table IV printed 767 is the direct table and adds the inner-scale branches 1.64 and 1.87 that olb lacks. The Ch. 12 uplink form Eq. (27) printed 492 carries the same 1.46; the Ch. 12 downlink form Eq. (22) printed 491 uses 1.45. R5 and R6 marked it `unmatched` inside Ch. 9 and Ch. 10, which only USE the radius. |
| PW-09 | sigma_I^2 = exp[0.54 s^2/(1+1.22 s^(12/5))^(7/6) + 0.509 s^2/(1+0.69 s^(12/5))^(5/6)] - 1 | olb/turbulence/plane_wave_scintillation.py:284-285 | point plane-wave scintillation index, all regimes | 9.4.1 | Ch. 9, Eq. (47) | 336 | 361 | Kolmogorov; plane wave; l0=0; L0=inf | wrong | CONFIRMED BY 3 READERS. The book gives 0.49, 1.11 and 0.51: Ch. 9, Eq. (47) printed 336; Ch. 12, Eqs. (40) printed 497 and (93) printed 522; App. III Table VII(b) printed 769. A full-text search of the book finds NO 0.54 and NO 0.509. The 1.22 belongs one level down, in q = L/(k rho0^2) = 1.22 sigma_R^(12/5) printed 342, which enters eta_X as 0.35 x 1.22/0.38 = 1.11. Two faults follow. (a) The weak limit becomes 1.049 sigma_R^2 instead of 1.000: measured +4.7 % at sigma_R = 0.3 and -0.5 % at sigma_R = 5. (b) The module disagrees with itself, because PW-11 at d = 0 gives 0.49/1.11/0.51/0.69. FIX: set 0.49, 1.11 and 0.51 here AND in the kernel (KR-24). |
| PW-10 | d = (k D^2/(4 L))^0.5 | olb/turbulence/plane_wave_scintillation.py:297 | aperture parameter | 10.3.2 | Ch. 10, Eq. (68) | 413 | 438 | - | exact | The book also relates it to the soft aperture by Omega_G = 16 L/(k D_G^2) = 4/d^2 with D_G^2 = 8 W_G^2. Ch. 12, Eq. (39) printed 496 scales its own aperture variable with k D_G^2/(16 L), a different definition. |
| PW-11 | sigma_I^2(D) = exp[0.49 s^2/(1+0.65 d^2+1.11 s^(12/5))^(7/6) + 0.51 s^2 (1+0.69 s^(12/5))^(-5/6)/(1+0.90 d^2+0.62 d^2 s^(12/5))] - 1 | olb/turbulence/plane_wave_scintillation.py:319-321 | aperture-averaged plane-wave index, weak to strong | 10.3.2 | Ch. 10, Eq. (69) | 413 | 438 | Kolmogorov; plane wave; l0=0; L0=inf; unobscured | exact | Every constant matches: 0.49, 0.65 d^2, 1.11, 0.51, 0.69, 0.90 d^2, 0.62 d^2. The book valid range is 0 to infinity in sigma_R^2. The d = 0 limit is exactly Ch. 9, Eq. (47) and Ch. 12, Eq. (40), which is what condemns PW-09. |
| PW-12 | A = [1 + 1.07 d^(7/3)]^(-1) | olb/turbulence/plane_wave_scintillation.py:336 | weak-turbulence aperture-averaging factor, small inner scale | - | not in the book | - | - | Kolmogorov; l0=0; L0=inf; weak; plane wave; unobscured | other source | Churnside, Appl. Opt. 30 (1991) 1982, DOI 10.1364/AO.30.001982, the book Ch. 10 reference [12] printed 438. CITATION FAULT: the docstring credits "Andrews and Phillips Ch. 10". The book form is DIFFERENT: Ch. 10, Eq. (61) printed 412 is A = [1 + 1.062 (k D_G^2/(4 L))]^(-7/6), with the 7/6 OUTSIDE the bracket and the constant 1.062. Same D^(-7/3) large-aperture asymptote; the two differ by up to 12 %, largest near d^2 = 1. The constant 1.07 is nowhere in the 809-page book. |
| PW-13 | A = [1 + 2.21 (D/l0)^(7/3)]^(-1) | olb/turbulence/plane_wave_scintillation.py:353 | weak-turbulence aperture-averaging factor, large inner scale | - | not in the book | - | - | weak; plane wave; unobscured | other source | Churnside 1991, DOI 10.1364/AO.30.001982. CITATION FAULT: the docstring credits "Andrews and Phillips Ch. 10". Andrews handles a finite inner scale only through the two-scale forms Ch. 10, Eqs. (62)-(68). The constant 2.21 is nowhere in the book. |
| PW-14 | A (strong, Churnside two-term, 0.908 and 0.162) | olb/turbulence/plane_wave_scintillation.py:377-379 | strong-turbulence aperture-averaging factor, small inner scale | - | not in the book | - | - | strong; plane wave; unobscured | other source | Churnside 1991, DOI 10.1364/AO.30.001982. Andrews cites Churnside and plots his curve in Figs. 10.11-10.12 printed 418, but never prints his formula. The Andrews replacement for the same job is Ch. 10, Eq. (69) (see PW-11). NOTE: this function calls PW-09, so the wrong constants propagate into A_strong. |
| PR-01 | `DEFAULT_HS = logspace(1 m, 20 km, 20)` | olb/turbulence/profiles.py:13 | turbulence altitude grid | - | - | - | - | - | unmatched | A numerical grid choice, not a book equation. The grid stops at 20 km, so it truncates the high-altitude Cn2 layer that the Hufnagel-Valley model puts near 10 km. |
| PR-02 | `default_cn2_profile` calls `get_c2n(hs, wind_rms, cn2_ground)` | olb/turbulence/profiles.py:34 | zenith Cn2(h) profile | 12.2.1 | Ch. 12, Eq. (1) | 481 | 506 | - | exact | R1 marked this `unmatched` inside Ch. 3, because Sec. 3.3.4 defers the Cn2(h) models to Sec. 12.2.1. R7 resolved it against the kernel (see KR-27). |
| UF-01 | WEAK_FLUCTUATION_LIMIT = 0.25 | olb/turbulence/uplink_flux.py:72 | validity ceiling on the log-amplitude variance sigma2_x | 8.2 | Ch. 8, text below Eq. (23) | 264-265 | 289-290 | weak | reduction | RESOLVED. The book defines weak fluctuation by sigma_R^2 < 1 (Ch. 8 printed 264; Ch. 12, Eqs. (40) printed 497 and (93) printed 522). With sigma_I^2 = 4 sigma_x^2 (Ch. 8, Eq. (13)) the book limit is sigma_x^2 < 0.25. The code PREVIOUSLY allowed 0.6 (sigma_R^2 = 2.4, so the warning fired 2.4 times too late); it now uses 0.25, the exact book-derived limit. The 0.6 came from the Dios split-step comparison, not from Andrews. The 0.6 in the `beam_wave_scintillation.py:25` docstring is a DIFFERENT, correct statement: the Dios model's own agreement range against a split-step reference, not this gate. |
| UF-02 | TODO closed form 3.86 sigma_R^2 {0.40[(1+2Theta)^2+4Lambda^2]^(5/12) cos[(5/6)atan((1+2Theta)/(2Lambda))] - (11/16)Lambda^(5/6)} | olb/turbulence/uplink_flux.py:92-95 | validation target for the Dios integrator | 8.2 | Ch. 8, Eq. (23); repeated as Summary Eq. (130) | 264 (303) | 289 (328) | Kolmogorov; l0=0; L0=inf; weak; on-axis; homogeneous Cn2 | exact | FOUND. The comment reproduces Eq. (23) verbatim (its longitudinal half; the comment drops the leading 4.42 sigma_R^2 Lambda^(5/6) r^2/W^2 radial term, which is 0 on axis). The book states Eq. (23) holds "in the case of a collimated OR divergent beam", so it IS the correct validation target for the diverged feed of olb gap 9. Also Ch. 9, Eqs. (92), (93) and (148). |
| UF-03 | Lambda0 = 2L/(k W0^2) | olb/turbulence/uplink_flux.py:98 | input-plane Fresnel parameter | 8.1 | Ch. 8, Eq. (5) | 261 | 286 | - | exact | Confirmed by 2 readers. Also Ch. 12, Eq. (8) printed 488. |
| UF-04 | Theta0 = 1 - L/F0, with F0 < 0 for a diverging wavefront | olb/turbulence/uplink_flux.py:100-103 | input-plane curvature parameter | 8.1 | Ch. 8, Eq. (5) | 261 | 286 | - | exact | Confirmed by 2 readers. Also Ch. 12, Eq. (8) printed 488. The book states the plane limit Theta = 1, Lambda = 0 and the spherical limit Theta = Lambda = 0 on the same page, so Theta0 > 1 for a divergent beam agrees. |
| UF-05 | f0 = -(d + zR^2/d) | olb/turbulence/uplink_flux.py:102 | phase-front radius of a Gaussian beam at distance d from its waist | 4.4.1 | Ch. 4, Eqs. (38), (50), (66) | 93, 96, 101 | 118, 121, 126 | - | unmatched | The free-space R(z) = z + zR^2/z is in Ch. 4, but no reader verified this exact expression against a book equation. See Table 2 row G-15. |
| UF-06 | Theta = Theta0/(Theta0^2 + Lambda0^2) | olb/turbulence/uplink_flux.py:104 | output-plane curvature parameter | 8.1 | Ch. 8, Eq. (6) | 261 | 286 | - | exact | Confirmed by 2 readers. Also Ch. 12, Eq. (9) printed 489. |
| UF-07 | Z0 = L/sqrt(1/Theta - 1) | olb/turbulence/uplink_flux.py:105 | effective Rayleigh range | 8.1 | Ch. 8, Eq. (6) | 261 | 286 | - | exact | Algebraic inverse of the kernel parameterisation Theta = [1 + (L/Z0)^2]^-1, so the diverged Theta of Eq. (6) reaches the kernel through Z0. It works around the collimated-only kernel row KR-15. |
| UF-08 | beta2 += 2 (sigma_theta L)^2 | olb/turbulence/uplink_flux.py:183 | mechanical jitter folded into the beam-wander displacement variance | 12.6.2 | Ch. 12, Eqs. (50), (51), (53) | 502-503 | 527-528 | weak; homogeneous Cn2 | approximate | Confirmed by 2 readers. Andrews keeps the wander displacement <r_c^2> (Eqs. (50) and (51)) and the pointing-error variance sigma_pe^2 (Eq. (53)) as SEPARATE quantities that share one integral, then feeds alpha_pe into the untracked index Eq. (54) and the rms wander angle into the tracked index Eq. (57). Andrews NEVER adds a mechanical tracking jitter into <r_c^2>. The olb variance addition is a defensible extension, but it is not the book construction and it has no Andrews citation. See Conflicts C-09. |
| UF-09 | w_lt | olb/turbulence/uplink_flux.py:184 | uplink long-term waist from the wander variance | 12.6.1 | Ch. 12, Eq. (48); repeated as Eq. (97) | 500 (523) | 525 (548) | weak | approximate | Andrews builds W_LT from reciprocity, W_LT = W[1 + (D0/r0)^(5/3)]^(1/2) with D0^2 = 8 W0^2, or from Eq. (97) with 4.35 mu2u. olb adds the wander variance beta2 to the short-term waist instead. |
| UF-10 | betax, betay, beta | olb/turbulence/uplink_flux.py:190-192 | per-axis normal draws of variance beta2/2, Rayleigh radius | 8.3 | Ch. 8, Eqs. (32) and (33) | 271-272 | 296-297 | homogeneous Cn2 | exact | Consistent with <r_c^2> being the two-dimensional displacement variance. Also Ch. 12, Eq. (50) printed 502. |
| UF-11 | Is_summed x (w_free/w_st)^2 | olb/turbulence/uplink_flux.py:209 | rescale the flux from the short-term waist onto the free-space baseline | 12.5.1 | Ch. 12, Eq. (32); repeated as Eq. (66) | 494 (510) | 519 (535) | on-axis | exact | Confirmed by 2 readers. The on-axis mean irradiance is W0^2/W_LT^2; the rescale is that ratio referred to the free-space width. Also Ch. 8, Eq. (32) printed 271. |
| KR-01 | r0s = (0.42 k^2 INT Cn2 ((L-z)/L)^(5/3) dz)^(-3/5) | my_analysis_modules/coupled_flux.py:46 | spherical-wave coherence diameter | 6.8.2 | Ch. 6, Eq. (116) | 209 | 234 | Kolmogorov; l0=0; L0=inf | wrong | The coefficient agrees (1.46 for rho becomes 0.423 for r0). The WEIGHT IS MIRRORED: the book uses (z/L)^(5/3) with z from the transmitter, and says the receiver end carries the weight. The kernel weights the transmitter end. It matches Dios, DOI 10.1364/AO.43.003866, Eq. (3). Same caution as GF-18. See Conflicts C-02. |
| KR-02 | W^2 = W0^2 (1 + L^2/Z0^2) | my_analysis_modules/coupled_flux.py:83 | free-space beam radius squared | 6.6.3 | Ch. 6, Eq. (86), the W^2 factor | 202 | 227 | collimated | exact | |
| KR-03 | 2[(4.2 L/(k r0s))(1 - 0.26 (r0s/W0)^(1/3))]^2 | my_analysis_modules/coupled_flux.py:86 | short-term waist turbulence term | 6.6.3 | Ch. 6, Eq. (101) | 206 | 231 | - | other source | The book short-term radius is W sqrt(1 + 1.33 s_R^2 Lambda^(5/6)[1 - 0.66(Lambda0^2/(1+Lambda0^2))^(1/6)]), a Rytov-variance form. The 4.2/0.26 form is the Yura and Fried coherence-radius form used by Dios, DOI 10.1364/AO.43.003866. Neither 4.2 nor 0.26 is in Ch. 6 or Ch. 7. |
| KR-04 | <beta^2> = 2.07 INT Cn2 (L-z)^2 W(z)^(-1/3) dz | my_analysis_modules/coupled_flux.py:113 | beam-wander displacement variance | 6.6.1 | Ch. 6, Eq. (93) with kappa0 = 0; Eqs. (117)-(118) | 203 (209-210) | 228 (234-235) | Kolmogorov; l0=0; L0=inf; frozen atmosphere | other source | ADJUDICATED 2026-08-25 against the Dios paper itself. THE KERNEL COPIES DIOS CORRECTLY: Dios Eq. (11), printed p. 3868, prints the constant 2.07, the weight (L-z)^2 and the factor [1/W_s(z)]^(1/3). Dios Eq. (10) makes <beta^2> the RADIAL variance. So the kernel is not a mis-copy, and the status moves from `wrong` to `other source`. THE INTEGRAND IS IDENTICAL, THE COEFFICIENT IS NOT: the book gives 7.25 (R3 re-derived it from Eq. (88) with the filter of Eq. (89): 8 pi^2 (0.033)(1/2) Gamma(1/6) = 7.252), and both quantities are radial, so the 3.50 gap stands as a source-against-source difference. Dios validates his Eq. (11) against a split-step (FFT-BPM) simulation of the same uplink in his Fig. 3, printed p. 3871, and the two agree. See Conflicts C-01. |
| KR-05 | W_LT = sqrt(W_ST^2 + 2 <beta^2>) | my_analysis_modules/coupled_flux.py:127 | long-term waist | 6.6.3 | Ch. 6, Eq. (100) | 205 | 230 | - | other source | ADJUDICATED 2026-08-25. The kernel copies Dios Eq. (1), printed p. 3867 (repeated as Eq. (29), printed p. 3870), which prints W_LT^2 = W_ST^2 + 2<beta^2> with a RADIAL <beta^2> (Dios Eq. (10)). So the factor 2 is the paper's own factor, NOT a per-axis conversion, and `olb/turbulence/angle_of_arrival.py:57` is right to call the quantity radial. The book gives W_LT^2 = W_ST^2 + <r_c^2>, factor 1 on a radial <r_c^2>. With each source's own constant the wander part of W_LT^2 is 1.38 Cn2 L^3 W0^(-1/3) by Dios and 2.42 by Andrews, so the two rules differ by 1.75, not 3.50. See Conflicts C-03. |
| KR-06 | 2(4 L/(k r0s))^2 | my_analysis_modules/coupled_flux.py:167 | long-term waist, collimated | 6.6.3 | Ch. 6, Eq. (86); Eq. (124) upper | 202 (212) | 227 (237) | collimated | other source | Dios Eq. (2), DOI 10.1364/AO.43.003866. The book form W sqrt(1 + 1.33 s_R^2 Lambda^(5/6)) is linear in Cn2; the r0s form goes as Cn2^(6/5). They are the weak-Rytov and the coherence-limited far-field limits, not the same equation. |
| KR-07 | G_u goes to 1.33 s_R^2 Lambda^(5/6) | my_analysis_modules/coupled_flux.py:200 | uniform-Cn2 limit of the spreading integral | 6.6.1 | Ch. 6, Eq. (86) | 202 | 227 | homogeneous Cn2; Kolmogorov | exact | 4.35 x 3/8 / 1.23 = 1.326; the book prints 1.33. |
| KR-08 | Lambda = 2L/(k W^2(L)) | my_analysis_modules/coupled_flux.py:232 | output-plane beam parameter | 6.2.1 | Ch. 6, Eq. (7) | 183 | 208 | - | exact | |
| KR-09 | (1 - z/L)^(5/3), z from the transmitter | my_analysis_modules/coupled_flux.py:234 | path weight of the spreading integral | 6.8.1 | Ch. 6, Eq. (109) | 208 | 233 | - | exact | CORRECT orientation here. Contrast KR-01, which mirrors it. Also Ch. 8, Eq. (118) printed 299. |
| KR-10 | Phi_n = 0.033 Cn2 kappa^(-11/3) | my_analysis_modules/coupled_flux.py:235 | Kolmogorov spectrum | 6.8.1 | Ch. 6, Eq. (110) | 208 | 233 | Kolmogorov; l0=0; L0=inf | exact | |
| KR-11 | 4 pi^2 (0.033)(-0.5 Gamma(-5/6)) = 4.3508 | my_analysis_modules/coupled_flux.py:235 | closed form of the spectral integral | 6.8.1 | Ch. 6, Eq. (109) lower | 208 | 233 | Kolmogorov | exact | R3 evaluated it: Gamma(-5/6) = -6.6796, so the product is 4.3508. The book prints 4.35. |
| KR-12 | G_u = 4 pi^2 k^2 INT INT kappa Phi_n {1 - exp[-(Lambda L kappa^2/k)(1-z/L)^2]} | my_analysis_modules/coupled_flux.py:236 | turbulence spreading integral (the T integral) | 6.8.1 | Ch. 6, Eq. (109) upper | 208 | 233 | Kolmogorov | exact | R4 marked the same lines `unmatched` inside Ch. 8, because Ch. 8, Eq. (32) states only the split and refers the derivation to Sec. 6.6. R3 resolved it. |
| KR-13 | W_LT = W sqrt(1 + G_u) | my_analysis_modules/coupled_flux.py:238 | long-term beam waist | 6.8.1 | Ch. 6, Eq. (111) | 208 | 233 | Kolmogorov | exact | |
| KR-14 | Lambda = 2L/(k wL^2) | my_analysis_modules/coupled_flux.py:249-250 | output-plane beam parameter | 8.1 | Ch. 8, Eq. (6) | 261 | 286 | - | exact | |
| KR-15 | Theta = [1 + (L/Z0)^2]^-1 | my_analysis_modules/coupled_flux.py:253-254 | output-plane curvature parameter | 8.1 | Ch. 8, Eq. (6) | 261 | 286 | collimated | reduction | Only the collimated receiver Theta. olb works around it with the effective Z0 of UF-07. |
| KR-16 | _A | my_analysis_modules/coupled_flux.py:257-260 | A(z), the Gaussian attenuation scale | 8.2 | Ch. 8, Eqs. (14) and (17) exponent | 262-263 | 287-288 | Kolmogorov; l0=0; L0=inf; weak | other source | Dios Eq. (17), DOI 10.1364/AO.43.003866. Duplicate of BW-08: olb gap 10. |
| KR-17 | _B | my_analysis_modules/coupled_flux.py:263-266 | B(z), the phase term | 8.2 | Ch. 8, Eq. (17) cosine argument | 263 | 288 | Kolmogorov; l0=0; L0=inf; weak | other source | Dios Eq. (18), DOI 10.1364/AO.43.003866. Duplicate of BW-09. |
| KR-18 | on_axis_scintillation_index | my_analysis_modules/coupled_flux.py:288 | longitudinal scintillation index integrand | 8.2 | Ch. 8, Eq. (17); closed forms Eqs. (19) and (23) | 263-264 | 288-289 | Kolmogorov; l0=0; L0=inf; weak | exact | FIXED 2026-08, and the fix is CONFIRMED against the Dios paper on 2026-08-25: Eq. (16), printed p. 3869, closes the large bracket after the cosine, so the cosine multiplies only the second term. The kernel now agrees with Dios Eq. (16) and with Andrews Ch. 8, Eq. (17). The fault, for the record: MISPLACED PARENTHESIS. The code computed A^(5/6)[1 - (1+ratio^2)^(5/12)] cos((5/6)atan(ratio)). Dios Eq. (16) and Andrews Eq. (17) need A^(5/6)[1 - (1+ratio^2)^(5/12) cos((5/6)atan(ratio))], that is the cosine multiplies ONLY the second term. The error is minus A^(5/6)[1 - cos((5/6)atan(B/A))] per unit path, so it vanishes only where A goes to zero, and it is non-zero for every finite Gaussian beam. `uplink_flux._flux_result` reaches this function through `coupled_flux_sample` (line 377), so the olb uplink turbulence Term inherits it. The olb twin at `beam_wave_scintillation.py:137-139` is CORRECT. |
| KR-19 | off_axis_scintillation_index | my_analysis_modules/coupled_flux.py:310-316 | radial scintillation index | 8.2 | Ch. 8, Eqs. (16) and (18) | 263 | 288 | Kolmogorov; l0=0; L0=inf; weak | other source | Dios Eq. (20). Agrees with the olb twin BW-12. |
| KR-20 | sigma2_gauss = (sigma2_on + sigma2_off) x I_off^2 | my_analysis_modules/coupled_flux.py:379 | total index at the wander offset | 8.2 | Ch. 8, Eqs. (9), (15) and (117) | 261, 263 (299) | 286, 288 (324) | Kolmogorov; weak | other source | ADJUDICATED 2026-08-25, and the earlier patch is REVERSED. The weight IS in Dios: Eq. (25), printed p. 3870, gives sigma2_I,Gb = (sigma2_I + sigma2_I,r) <I>^2, with <I> the Eq. (24) mean irradiance at the wander offset. Dios says why in the text above Eq. (24): Eqs. (13), (16) and (20) normalise the index to the LOCAL mean irradiance, and Eq. (25) re-normalises it to the mean at the BEAM CENTRE, which is the normalisation that Eq. (26) needs. Section 5, step (c)(ii) tells the reader to use Eq. (25) at this point. An agent removed the weight in 2026-08 on the Andrews reading (Ch. 8, Eqs. (9) and (15) keep the local normalisation), which made the kernel disagree with the paper it cites. The weight is back. Andrews and Dios normalise differently; that is a source difference, not a bug. |
| KR-21 | sigma2_x = 0.25 ln(1 + sigma2_gauss) | my_analysis_modules/coupled_flux.py:382 | log-amplitude variance from the index | 8.2 | Ch. 8, Eq. (13) | 262 | 287 | weak | exact | Exact inversion of sigma_I^2 = exp(4 sigma_x^2) - 1. |
| KR-22 | xi drawn from N(-sigma2_x, sqrt(sigma2_x)) | my_analysis_modules/coupled_flux.py:385 | lognormal log-amplitude draw with a conserved mean irradiance | 9.11 | Ch. 9, Eq. (158) | 384 | 409 | weak | unmatched | Ch. 8 gives only the variance, not the PDF. The lognormal PDF with the mean offset is Ch. 9, Eq. (158) printed 384 and Ch. 12, Eqs. (65)-(66) printed 510, but no reader verified this exact kernel line against them. |
| KR-23 | rytov = (1.23 Cn2 k^(7/6) L^(11/6))^0.5 | my_analysis_modules/general_atmospherics.py:22 | plane-wave Rytov standard deviation | 5.2.2 | Ch. 5, Eq. (15) | 140 | 165 | Kolmogorov; plane wave; weak; homogeneous Cn2 | exact | Confirmed by 3 readers. Also Ch. 6, Eq. (119) printed 210 and App. III Table VII(a) footer printed 769. olb duplicates it two more times instead of importing it through `olb/_deps.py`: olb gap 10. |
| KR-24 | scint_pl = exp[0.54 r^2/(1+1.22 r^(12/5))^(7/6) + 0.509 r^2/(1+0.69 r^(12/5))^(5/6)] - 1 | my_analysis_modules/general_atmospherics.py:23 | plane-wave scintillation index | 9.4.1 | Ch. 9, Eq. (47) | 336 | 361 | Kolmogorov; plane wave; l0=0; L0=inf | wrong | CONFIRMED BY 2 READERS. Byte-for-byte the same constants as PW-09 and the same disagreement with the book (Ch. 12, Eqs. (40) and (93); App. III Table VII(b)). The kernel is the ORIGIN of the error; olb copied it. One fix must cover both. |
| KR-25 | r0 = (0.423 k^2 Cn2 L)^(-3/5) | my_analysis_modules/general_atmospherics.py:24 | plane-wave Fried parameter | 6.4.1 | Ch. 6, Eq. (64) lower, times 2.1 | 194 | 219 | Kolmogorov; plane wave; l0=0; L0=inf; homogeneous Cn2 | exact | Confirmed by 3 readers. Also Ch. 12, Eq. (23) printed 492, Ch. 14, Eq. (25) printed 617 and App. III Table IV with its note printed 767. 2.1^(-5/3) x 1.46 = 0.4236; the book prints 0.42. Third copy of the same equation (with GF-11 and AO-05): olb gap 10. |
| KR-26 | v(h) = (slew rate) h + Vg + 30 exp(-((h-9400)/4800)^2) | my_analysis_modules/general_atmospherics.py:41 | Bufton slew-plus-wind speed profile | 12.2.1 | Ch. 12, Eq. (3) | 481 | 506 | - | exact | R1 marked this `unmatched` inside Sec. 3.4, which holds only the Taylor hypothesis and has no wind profile. R7 resolved it: all three constants match. The code `ws` is in deg/s and is converted with deg2rad, so it is the Andrews slew rate in rad/s. |
| KR-27 | Cn2(h) = 0.00594 (v/27)^2 (1e-5 h)^10 exp(-h/1000) + 2.7e-16 exp(-h/1500) + A exp(-h/100) | my_analysis_modules/general_atmospherics.py:56-58 | Hufnagel-Valley Cn2 altitude profile | 12.2.1 | Ch. 12, Eq. (1) | 481 | 506 | - | exact | R1 marked this `unmatched` inside Ch. 3, which defers to Sec. 12.2.1. R7 resolved it: EVERY constant, exponent and scale height matches term by term. The defaults w = 21 m/s and A = 1.7e-14 are the book H-V5/7. |
| KR-28 | W_LT = W (1 + 1.33 sigma_R^2 Lambda^(5/6))^0.5 | my_analysis_modules/general_atmospherics.py:108 | long-term spot size | 6.6.1 | Ch. 6, Eq. (86) | 202 | 227 | weak; collimated | exact | R1 marked this `unmatched` inside Ch. 3 and App. III. R3 resolved it: Ch. 6, Eq. (86) is the source of the 1.33. |
| KR-29 | Lambda0 = 2L/(k W0^2); Lambda; W | my_analysis_modules/general_atmospherics.py:135, 152-155 | Gaussian beam parameters | 6.2.1 | Ch. 6, Eqs. (6) and (7); App. III Table VI footer | 183 (768) | 208 (793) | collimated | reduction | The kernel FIXES Theta0 = 1 at line 135 and HARDCODES the waist W0 = 0.125 m at line 152, so the caller's waist input is ignored. olb threads f0, so olb is the more general of the two. See olb gap 3. |
| KR-30 | W_LT = W (1 + (0.35/r0)^(5/3))^(3/5) | my_analysis_modules/general_atmospherics.py:162 | long-term spot size from r0 | - | - | - | - | - | unmatched | The length 0.35 m is undocumented. Not an Andrews form. |
| KR-31 | rho_meas = 2.1 x (blob radius at 1/e) | my_analysis_modules/general_atmospherics.py:441 | Fried parameter from a measured coherence function | 6.4.1 | Ch. 6, text below Eq. (64) | 194 | 219 | Kolmogorov; plane wave | exact | Also the App. III note under Table IV printed 767: "Fried's parameter is related by r0 = 2.1 rho_pl". |
| KR-32 | r0(lambda2) = r0(lambda1)(lambda2/lambda1)^(6/5) | my_analysis_modules/general_atmospherics.py:508 | Fried parameter wavelength scaling | 6.4.1 | Ch. 6, Eq. (64); App. III Table IV | 194 (767) | 219 (792) | Kolmogorov; plane wave | exact | Follows from r0 proportional to k^(-6/5). |
| KR-33 | r0_net = (sum r0_i^(-5/3))^(-3/5) | my_analysis_modules/general_atmospherics.py:526 | net Fried parameter of N layers | - | - | - | - | Kolmogorov | unmatched | It follows from the additive Cn2 path integral, but the book does not state it. |
| KR-34 | phase-screen amplitude proportional to (f^2)^(-11/12) | my_analysis_modules/general_atmospherics.py:542, 588 | Kolmogorov phase-screen Fourier amplitude | 3.3.1 | Ch. 3, Eq. (18) | 67 | 92 | Kolmogorov; l0=0; L0=inf | exact | The amplitude is the square root of a two-dimensional kappa^(-11/3) power spectrum, so the exponent is -11/6 = 2 x (-11/12). |
| KR-35 | Roddier random-tilt scale 3.30 | my_analysis_modules/general_atmospherics.py:609 | low-frequency tilt injection in a phase screen | - | - | - | - | Kolmogorov | other source | Roddier, Progress in Optics 19 (1981) 281. |
| KR-36 | var = 1.0299 (D/r0)^(5/3); var_notilt = 0.134 (D/r0)^(5/3) | my_analysis_modules/general_atmospherics.py:785-786 | phase variance over a circular aperture, with and without tilt | 14.5.4 | Ch. 14, Eqs. (90) and (94) | 635-636 | 660-661 | Kolmogorov; unobscured | exact | R1 marked this `other source` (Noll, JOSA 66 (1976) 207, DOI 10.1364/JOSA.66.000207) inside Ch. 3 and App. III. R7 resolved it: the book prints 1.03 and 0.13 for the same two quantities. The Noll unrounded values sit inside the Andrews pair (Eq. (91) gives 1.02 by the ABCD route). Duplicate of AO-01 and AO-02: olb gap 10. |

## Table 2 — reverse map

Book capabilities that olb lacks. Sorted by book section. `closes olb gap` uses
the numbered list: 1 aperture angle-of-arrival tilt feeds no Term (0-W3); 2 NO-SCINTILLATION
pre-compensated uplink; 3 collimated-only Fried parameter; 4 unused
strong-regime effective beam parameters; 5 missing gamma-gamma; 6 no inner or
outer scale; 7 no temporal statistics; 8 annular aperture unmodelled; 9 Dios
diverged feed unvalidated against the Andrews closed form; 10 duplicate physics
copies. `target module` names the planned home in `olb/turbulence/andrews/`.

| gap id | book section | book equation | printed p | pdf p | capability | target module | closes olb gap | priority |
|---|---|---|---|---|---|---|---|---|
| G-01 | 3.2.3 | Ch. 3, Eqs. (11) and (14) | 64-65 | 89-90 | Cn2 from measured temperature: n = 1 + 79e-6 P/T and Cn2 = (79e-6 P/T^2)^2 CT2. Converts site weather data into Cn2. | paths.py | - | P3 (site data conversion, not link physics) |
| G-02 | 3.2.3 | Ch. 3, Eq. (13) | 64 | 89 | Inner scale from the dissipation rate, l0 = 7.4 (nu^3/epsilon)^(1/4). Gives a physical default l0 for the new spectrum models. | spectra.py | 6 | P3 |
| G-03 | 3.3 | Ch. 3, Eqs. (15)-(17) | 66 | 91 | The spectrum, covariance and structure-function transform triple, Dn(R) = 8 pi INT kappa^2 Phi_n(kappa)[1 - sin(kappa R)/(kappa R)] dkappa. Lets any new spectrum give its own structure function. | structure.py | 6 | P2 |
| G-04 | 3.3.1 | Ch. 3, Eq. (18); restated as Eq. (28) | 67 | 92 | An explicit Phi_n(kappa, Cn2) Kolmogorov spectrum function. olb hardcodes 0.033 kappa^(-11/3) in two modules with no shared definition. | spectra.py | 10 | P1 |
| G-05 | 3.3.2 | Ch. 3, Eq. (19); restated as Eq. (29) | 67 | 92 | Tatarskii spectrum, 0.033 Cn2 kappa^(-11/3) exp(-kappa^2/km^2) with km = 5.92/l0. Gives inner-scale truncation. | spectra.py | 6 | P1 |
| G-06 | 3.3.2 | Ch. 3, Eq. (20); restated as Eq. (30) | 68 | 93 | von Karman and modified von Karman spectrum, 0.033 Cn2 exp(-kappa^2/km^2)/(kappa^2+k0^2)^(11/6) with k0 = 2 pi/L0. Gives an outer scale, and both scales together. Restated as Ch. 8, Eq. (25) printed 265 with km = 5.92/l0 and k0 = 2 pi/L0. | spectra.py | 6 | P1 |
| G-07 | 3.3.2 | Ch. 3, Eq. (21) | 68 | 93 | Exponential spectrum, 0.033 Cn2 kappa^(-11/3)[1 - exp(-kappa^2/k0^2)] with k0 = C0/L0. The book uses C0 = 8 pi in the Ch. 9 scintillation model and C0 = 4 pi to approximate von Karman. | spectra.py | 6 | P2 |
| G-08 | 3.3.3 | Ch. 3, Eq. (22); restated as Eq. (31) | 69 | 94 | Modified atmospheric (Andrews-Hill) spectrum with the high-wavenumber bump, [1 + 1.802(kappa/kl) - 0.254(kappa/kl)^(7/6)] exp(-kappa^2/kl^2)/(kappa^2+k0^2)^(11/6) with kl = 3.3/l0. The book states this is the only listed model that has the bump, and that the bump matters for scintillation. Ch. 9, Eqs. (3), (5) and (6) printed 327-328, restated as Eqs. (145)-(147) printed 380-381, give the same filters as f(kappa l0) and g(kappa L0) with kappa_0 = 8 pi/L0. | spectra.py | 6 | P1 |
| G-09 | 3.3.3 | Ch. 3, Eq. (23) | 69 | 94 | Modified atmospheric spectrum in the alternative outer-scale form, with k0 = 4 pi/L0 (or 2 pi/L0, or 8 pi/L0). Needed to match the Ch. 9 scintillation results, which use their own convention. | spectra.py | 6 | P2 |
| G-10 | 3.3.3 | Ch. 3, Eqs. (24) and (25) | 71 | 96 | Refractive-index structure function with a finite inner scale, from the Tatarskii and the modified spectrum, using the confluent hypergeometric function. Shows the quadratic roll-off that olb cannot represent. | structure.py | 6 | P3 (reference check only; olb needs the wave structure function, not Dn) |
| G-11 | 3.4 | Ch. 3, Eq. (27) | 73 | 98 | The Taylor frozen-turbulence hypothesis, u(R, t+t') = u(R - V_perp t', t). This is the entry point for every temporal quantity (Greenwood frequency, fade duration, temporal covariance). olb has no temporal axis at all. | temporal.py | 7 | P1 |
| G-12 | 3.4 | text, printed 72-73 | 72 | 97 | The two atmospheric time scales: advection time L0/V_perp, about 1 s, and eddy turnover, about 10 s. Also the stated failure mode of the Taylor hypothesis when V_perp is small or the wind is along the line of sight. Bounds the validity of any temporal model and of the snapshot assumption olb makes now. | temporal.py | 7 | P2 |
| G-13 | 4.4.1 | Ch. 4, Eq. (33) | 92 | 117 | Input curvature parameter Theta0 = 1 - z/F0 for a convergent (Theta0 < 1) or a divergent (Theta0 > 1) beam. The single missing input that makes the Fried parameter curvature-general. | beam.py | 3 | P1 |
| G-14 | 4.4.1 | Ch. 4, Eqs. (37) and (49); Eq. (139) | 93, 96 (119) | 118, 121 (144) | Spot radius W = W0 sqrt(Theta0^2+Lambda0^2) = W0/sqrt(Theta^2+Lambda^2) for a general Theta0. An independent cross-check of `olb.beam.free_space_radius`, which reaches the same answer by a virtual waist. | beam.py | 9 | P1 |
| G-15 | 4.4.1, 4.4.2 | Ch. 4, Eqs. (38), (50), (66); Eq. (140) | 93, 96, 101 | 118, 121, 126 | Phase-front radius of curvature F at the receiver, in three equivalent forms. Needed to know the wavefront that arrives at a focusing telescope, and to chain to the fibre focal plane. Also closes the unmatched row UF-05. | beam.py | - | P2 |
| G-16 | 4.4.2 | Ch. 4, Eqs. (44), (45) and (47) | 95 | 120 | Output beam parameters Theta, Theta_bar = 1 - Theta and Lambda from a GENERAL Theta0, plus the receiver-plane identities Theta = 1 - z/F and Lambda = 2z/(k W^2). | beam.py | 3 | P1 |
| G-17 | 4.4.2 | Ch. 4, Eqs. (41), (51), (52); Eq. (141) | 94, 96 (119) | 119, 121 (144) | On-axis and off-axis free-space irradiance, (Theta^2+Lambda^2) exp(-2 r^2/W^2). Out of scope: olb already gets the same number from the spot radius in the geometric Term. | out of scope | - | P3 |
| G-18 | 4.5, 4.5.1 | Ch. 4, Eqs. (55), (57), (58), (60); Eqs. (142)-(144) | 96, 98 (119) | 121, 123 (144) | Focusing parameter Theta_f = 2 F0/(k W0^2), spot at the geometric focus, waist distance and waist radius. This is the closed-form replacement for the olb virtual-waist recast, and it handles a CONVERGENT beam, which the olb recast cannot. | beam.py | 3 | P1 |
| G-19 | 4.9.1 | Ch. 4, Eqs. (149)-(151); Eq. (152) | 121-122 | 146-147 | Three-plane beam parameters through a "Gaussian lens" (a thin lens plus a finite aperture stop, with Theta_G = 2 L1/(k W_G^2)). Gives the beam parameters AT the fibre focal plane, aperture truncation included, in place of the present ad-hoc focal-spot size in the coupling models. | beam.py | - | P2 |
| G-20 | 5.2.2 | Ch. 5, Eq. (16) | 140 | 165 | The correct Gaussian-beam weak-fluctuation gate: sigma_R^2 < 1 AND sigma_R^2 Lambda^(5/6) < 1, with Lambda = 2L/(k W^2). Fixes TL-05, where a beam-wave Term is gated by a plane-wave threshold. | scintillation.py | - | P1 |
| G-21 | 5.2.2 | Ch. 5, Eq. (17) | 141 | 166 | Spectrum-independent weak gate, q = L/(k rho_pl^2) < 1 and q Lambda < 1. Works for ANY refractive-index spectrum, so it stays valid after an inner or an outer scale is added. | scintillation.py | 6 | P2 |
| G-22 | 5.7.1 | Ch. 5, Eq. (87) | 155 | 180 | Modified Rician (Rice-Nakagami) irradiance PDF from the first Born approximation. DO NOT BUILD: Andrews printed 155 reports that the theoretical moments lie below measured data and that the model is not suitable for irradiance fluctuations. Record the finding, not the code. | distributions.py | - | P3 |
| G-23 | 5.7.1, 5.7.2 | Ch. 5, Eqs. (88) and (94) | 155, 157 | 180, 182 | Normalised-moment tests for the lognormal and the Rician. A cheap validation face for any fitted irradiance PDF, and the standard way to say which PDF fits the data. | distributions.py | - | P2 |
| G-24 | 5.7.2 | Ch. 5, Eq. (93) | 156 | 181 | DONE (2026-09-04, backlog I-2). The lognormal PDF has one home: `andrews.distributions` `lognormal_params`/`_mean_log`/`_quantile`/`_rvs`, turned into the three Term faces by the one adapter `models/fade.py` `irradiance_fade_term`. Both `downlink._lognormal_term` and `terrestrial_scintillation_term` now build through it (was the inline DL-01..04 / TL-01..04 copies). | distributions.py | 10 | P2 |
| G-25 | 5.7.2 | Ch. 5, Eq. (95) | 157 | 182 | Forward map sigma_I^2 = exp(4 sigma_chi^2) - 1. olb uses only the inverse. The forward direction is needed to feed a log-amplitude phase screen or an adaptive-optics residual into a scintillation index. | distributions.py | 2 | P2 |
| G-26 | 5.9.1 | Ch. 5, Eqs. (103), (104), (105) | 160-161 | 185-186 | Two-scale split n1 = n_LS + n_SS and Phi_n = Phi_n,LS + Phi_n,SS, with the large and the small scales uncorrelated. The structural assumption under G-27. | scintillation.py | 4 | P2 |
| G-27 | 5.9.2 | Ch. 5, Eq. (108) | 163 | 188 | Effective spectrum Phi_n,e = Phi_n G_X + Phi_n G_Y with the large-scale and small-scale spatial filters. This is the parent theory of every strong-regime "effective" form olb already uses. Building it makes those forms derivable, not copied. | scintillation.py | 4 | P1 |
| G-28 | 5.9.3 | Ch. 5, Sec. 5.9.3 (a list, not numbered) | 165 | 190 | The six special scale sizes that set the filter cutoffs: inner scale, coherence radius, first Fresnel zone, beam radius, scattering disk, outer scale. Andrews states which set matters for which quantity: beam wander uses W and L0; scintillation uses the Fresnel zone, rho0 and l0. This is the entry point for an inner and an outer scale. | scintillation.py, spectra.py | 6 | P1 |
| G-29 | 5.10 | Ch. 5, text after Eq. (117) | 167 | 192 | Recorded caveat, no code: Andrews states the irradiance is not truly lognormal, because the second-order Rytov term is not Gaussian, and that the lognormal does not fit simulation data well in the tails. The olb quantile face lives exactly in that tail. Add this to the Term validity string. | out of scope | - | P2 |
| G-30 | 6.4.1 | Ch. 6, Eqs. (62) and (63); App. III Table I | 194 (765) | 219 (790) | Plane-wave wave structure function with an inner and an outer scale, and its two asymptotic regimes. The Kolmogorov row, D = 2.914 Cn2 k^2 L r^(5/3), is the base of every coherence radius in olb. | structure.py | 6 | P1 |
| G-31 | 6.4.2 | Ch. 6, Eqs. (69) and (70); App. III Table II | 195 (765) | 220 (790) | Spherical-wave wave structure function with an inner scale, Kolmogorov row D = 1.093 Cn2 k^2 L r^(5/3). Needed for a short terrestrial link fed by a point source. | structure.py | 6 | P2 |
| G-32 | 6.4.3 | Ch. 6, Eqs. (75), (76), (77); App. III Table III | 196-197 (766) | 221-222 (791) | Gaussian-beam wave structure function, general, with the inner-scale and outer-scale factor. The Kolmogorov row carries a = (1 - Theta^(8/3))/(1 - Theta), which is the general-curvature result olb needs. | structure.py | 3, 6 | P1 |
| G-33 | 6.4.3 | Ch. 6, Eq. (78) upper; App. III Table IV | 198 (767) | 223 (792) | Coherence radius when r0 is much less than l0, and the plane-wave inner-scale branches rho_pl = (1.64 Cn2 k^2 L l0^(-1/3))^(-1/2) (von Karman) and 1.87 (modified). olb has only the Kolmogorov branch, so it overestimates rho_pl in strong turbulence with a real inner scale. | structure.py | 6 | P1 |
| G-34 | 6.5 | Ch. 6, Eq. (84); definition Eq. (82) | 201 (200) | 226 (225) | Aperture angle-of-arrival variance, <beta_a^2> = 2.91 Cn2 L (2 W_G)^(-1/3), one axis, plane wave, Kolmogorov. NOW IMPLEMENTED in `structure.angle_of_arrival_variance`; `angle_of_arrival.aperture_arrival_angle_variance` delegates to it (AA-02 closed at module level; no Term consumes it yet, 0-W3). The slant-path version is Ch. 12, Eq. (28) printed 492, repeated as Eq. (90) printed 522: 2.91 mu0 sec(zeta)(2 W_G)^(-1/3), uplink and downlink. | structure.py, aperture.py | 1 | P1 |
| G-35 | 6.5 | Ch. 6, Eq. (83) | 200 | 225 | Angle of arrival with an inner and an outer scale: 1.64 Cn2 L l0^(-1/3)[1 - 0.72(kappa0 l0)^(1/3)] for a small aperture, and the [1 - 0.81(2 kappa0 W_G)^(1/3)] outer-scale reduction for a large one. | structure.py | 1, 6 | P2 |
| G-36 | 6.5 | Ch. 6, Eq. (85) | 201 | 226 | Rms image jitter equals the focal length times the rms angle of arrival. This is the focal-spot displacement that the fibre-coupling Terms need. | structure.py | 1 | P1 |
| G-37 | 6.6.1 | Ch. 6, Eqs. (86) and (87) | 202 | 227 | Partition of the spreading integral into a small-scale and a large-scale part, which separates diffraction, beam breathing and beam wander. It is the derivation that justifies not counting tip-tilt two times. | beam.py | 10 | P3 |
| G-38 | 6.6.1, 6.6.2 | Ch. 6, Eqs. (88)-(99) | 203-204 | 228-229 | Beam wander from first principles: the large-scale filter exp(-kappa^2 W^2(z)), the general integral Eq. (93) with the coefficient 7.25, the collimated (2.42) and focused (2.72) closed forms, and the outer-scale forms Eqs. (97) and (99). Replaces the kernel copy that carries 2.07 (KR-04). | wander.py | 10 | P1 |
| G-39 | 6.6.3 | Ch. 6, Eqs. (100) and (101) | 205-206 | 230-231 | Short-term beam radius from W_LT and the wander variance, in Rytov-variance form. Gives an independent check on the kernel 4.2/0.26 form (KR-03) and settles the factor in KR-05. | beam.py | 10 | P2 |
| G-40 | 6.7 | Ch. 6, Eqs. (102)-(105) | 206 | 231 | Angular spectrum of the mutual coherence function; the Gaussian approximation exp(-theta^2/(4 theta_c^2)) with theta_c = 1/(k rho_pl); spectral width 2 theta_c = lambda/(pi rho_pl). | temporal.py | 7 | P2 |
| G-41 | 6.7 | Ch. 6, Eqs. (106) and (107) | 207 | 232 | Temporal frequency spectrum from the mutual coherence function and Taylor frozen turbulence, exp(-omega^2/(4 omega_c^2)) with omega_c = V_perp/rho_pl. Gives the fade bandwidth and the adaptive-optics or tracking servo bandwidth. | temporal.py | 7 | P1 |
| G-42 | 6.8.1 | Ch. 6, Eqs. (108), (109), (111) | 208 | 233 | Slant-path mean irradiance, the spreading integral, and W_LT = W sqrt(1 + 4.35 k^(7/6) Lambda^(5/6) L^(5/6) INT Cn2(z)(1-z/L)^(5/3) dz). The kernel already has this (KR-12, KR-13); it needs an olb-native home. | paths.py | 9, 10 | P1 |
| G-43 | 6.8.2 | Ch. 6, Eqs. (114) and (115); App. III Table VI | 209 (768) | 234 (793) | Slant-path Gaussian-beam wave structure function and coherence radius, general in Theta and Lambda, so general in curvature: r0 = [8/(3(a + 0.618 Lambda^(11/6)))]^(3/5)(1.46 Cn2 k^2 L)^(-3/5) with Theta0 = 1 - L/F0 kept general. This is exactly the CLAUDE.md "next task" target, and it also gives the inner-scale branch. | paths.py, structure.py | 3 | P1 |
| G-44 | 6.8.3 | Ch. 6, Eqs. (117) and (118) | 209-210 | 234-235 | Slant-path beam wander, general in Theta0 and with a height-dependent outer-scale wavenumber. | wander.py | 1, 6 | P1 |
| G-45 | 6.9 | Ch. 6, Eq. (120) | 210 | 235 | Mean (coherent) field attenuation, exp(-0.39 Cn2 k^2 L kappa0^(-5/3)). Needed only for a coherent-detection budget. | out of scope | - | P3 |
| G-46 | 7.3, 7.4.1 | Ch. 7, text at Eq. (53) and Fig. 7.2 | 241, 244 | 266, 269 | Validation benchmarks: r0 goes to 2.27 rho_pl (exact, Belen'kii and Mironov) and 2.11 rho_pl from the effective-parameter form, as q goes to infinity. Gives a unit test for GF-09. | out of scope | 4 | P2 |
| G-47 | 7.4 | Ch. 7, Eq. (57); App. III Table VI note; Ch. 9, Eqs. (85) and (86) | 242 (768, 349) | 267 (793, 374) | WIRE UP the effective beam parameters. Andrews gives W_LT = W sqrt(1 + 4 q Lambda/3), the effective curvature F_LT, and the rule to substitute Theta_e and Lambda_e for Theta and Lambda in EVERY strong-turbulence result. olb already codes Theta_e and Lambda_e at `gaussian_fried.py:98-99`, but NO caller uses them (GF-17). Feed them to `gaussian_fried_parameter_profile` and to the beam-wave scintillation, so strong turbulence spreads the beam. | beam.py | 4 | P1 |
| G-48 | 7.4.1 | Ch. 7, Eq. (59) upper | 243 | 268 | Effective-parameter coherence radius in the regime where r0 is much less than l0. | structure.py | 4, 6 | P3 |
| G-49 | 7.4.2 | Ch. 7, Eqs. (62) and (63) | 245-246 | 270-271 | Beam wander valid from weak to strong fluctuation: the filter uses W_LT(z), so the 7.25 integral gains the [1 + 1.63 s^(12/5) Lambda0 (1-xi)^(16/5)] factor. Validated against experiment in Figs. 7.5 and 7.6. | wander.py | 4 | P1 |
| G-50 | 8.2 | Ch. 8, Eqs. (14), (16), (17) | 262-263 | 287-288 | General spectral-integral form of the scintillation index for ANY Phi_n, not only Kolmogorov. | scintillation.py | 6 | P2 |
| G-51 | 8.2 | Ch. 8, Eqs. (18) and (19) | 263 | 288 | Closed-form radial (2.64 with 1F1) and longitudinal (3.86 with 2F1) components, with no numerical path integral. | scintillation.py | 9 | P2 |
| G-52 | 8.2 | Ch. 8, Eq. (21) | 264 | 289 | Exact Kolmogorov Gaussian-beam scintillation index in 2F1 and 1F1. The reference the Dios numerical integrator must reproduce. | scintillation.py | 9 | P1 |
| G-53 | 8.2 | Ch. 8, Eq. (23); App. III Table IX(a); Ch. 9, Eq. (102) | 264 (771, 352) | 289 (796, 377) | Simple algebraic Gaussian-beam scintillation index, sigma_I^2(r,L) = 4.42 sR^2 Lambda^(5/6) r^2/W^2 + 3.86 sR^2{0.40[(1+2Theta)^2+4Lambda^2]^(5/12) cos[(5/6) arctan((1+2Theta)/(2Lambda))] - (11/16)Lambda^(5/6)}. The book states it is valid for a collimated OR divergent beam, so it is exactly the closed form the `uplink_flux` TODO names (UF-02) and the independent check of the Dios path. Ch. 9, Eq. (102) gives the all-regime on-axis version. | scintillation.py | 9 | P1 |
| G-54 | 8.2.1 | Ch. 8, Eq. (29) | 266 | 291 | Plane-wave inner-scale scintillation index on the von Karman spectrum, in Qm = 35.05 L/(k l0^2). | scintillation.py | 6 | P1 |
| G-55 | 8.2.1 | Ch. 8, Eqs. (30) and (31) | 267 | 292 | Gaussian-beam scintillation index with BOTH an inner and an outer scale (von Karman), radial and longitudinal parts. | scintillation.py | 6 | P1 |
| G-56 | 8.3 | Ch. 8, Eq. (32) | 271 | 296 | Cited decomposition W_LT^2 = W^2 + W_TSS^2 + <r_c^2>. olb takes this from the external kernel with no book citation. | wander.py | 10 | P2 |
| G-57 | 8.3 | Ch. 8, Eq. (33) | 272 | 297 | Rms beam wander of a collimated beam, 0.69 (lambda L/(2W0))(2W0/r0)^(5/6), with r0 = (0.16 Cn2 k^2 L)^(-3/5). | wander.py | 2 | P1 |
| G-58 | 8.3 | Ch. 8, Eqs. (34) and (35) | 272 | 297 | Beam-jitter spatial-frequency filter with kappa_r = Cr/r0, Cr about 2 pi. Separates hot-spot dancing from beam jitter. | wander.py | 2 | P2 |
| G-59 | 8.3 | Ch. 8, Eqs. (36)-(38); App. III Table IX(b) footer | 273-274 (772) | 298-299 (797) | Rms wander-induced effective pointing error, collimated (0.48) and focused (0.54) closed forms in 2W0/r0, with the beam-wander displacement variance sigma_pe^2 = 7.25 Cn2 L^3 W0^(-1/3) INT(...) built in. This separates jitter from scintillation cleanly. | wander.py | 2 | P1 |
| G-60 | 8.3 | Ch. 8, Eq. (39) | 274 | 299 | Asymptotic pointing error for a small and a large 2W0/r0. Both limits drive the pointing error to zero. | wander.py | 2 | P2 |
| G-61 | 8.3.1 | Ch. 8, Eq. (40) | 274 | 299 | Untracked-beam longitudinal scintillation index: the conventional Rytov longitudinal term PLUS 4.42 sR^2 Lambda^(5/6) spe^2/W^2 from the wander-induced pointing error. App. III Table IX(b) printed 772 collates the same result. | scintillation.py | 2 | P1 |
| G-62 | 8.3.1 | Ch. 8, Eqs. (41) and (42) | 275 | 300 | Untracked-beam index across the whole beam profile, with a unit-step floor inside the pointing-error radius. | scintillation.py | 2 | P1 |
| G-63 | 8.3.2 | Ch. 8, Eqs. (43) and (44) | 275-276 | 300-301 | Tracked-beam scintillation index: the plain Rytov longitudinal term, and a radial term shifted by the rms wander displacement, with a unit step. This is the saturating replacement for the non-saturating Dios model. | scintillation.py | 2 | P1 |
| G-64 | 8.4 | Ch. 8, Eq. (50); Summary Eq. (133) | 280 (304) | 305 (329) | Gaussian-beam irradiance covariance function, with delta_t = 0.67 - 0.17 Theta. | aperture.py | 8 | P2 |
| G-65 | 8.4 | Ch. 8, Eqs. (48), (49), (51), (52) | 280-281 | 305-306 | Plane-wave and spherical-wave irradiance covariance, exact and approximate. | aperture.py | 8 | P2 |
| G-66 | 8.4 | Ch. 8, Eq. (53); Summary text | 281 (304) | 306 (329) | Irradiance correlation width: about 1.7 sqrt(L/k) plane, 3 sqrt(L/k) spherical, sqrt(L/k) collimated with Lambda0 about 1. Sets the point-receiver limit above which aperture averaging starts. | aperture.py | 8 | P2 |
| G-67 | 8.5 | Ch. 8, Eqs. (54) and (55) | 282 | 307 | Fourier pair between the temporal irradiance covariance and the power spectral density, under the Taylor frozen-turbulence hypothesis. | temporal.py | 7 | P1 |
| G-68 | 8.5.1 | Ch. 8, Eqs. (56) and (57) | 282-283 | 307-308 | Plane-wave temporal irradiance power spectrum and the Fresnel frequency omega_t = V_perp/sqrt(L/k): flat below omega_t, omega^(-8/3) above. | temporal.py | 7 | P1 |
| G-69 | 8.5.2 | Ch. 8, Eqs. (58) and (59) | 284 | 309 | Spherical-wave temporal covariance and power spectrum. | temporal.py | 7 | P2 |
| G-70 | 8.5.3 | Ch. 8, Eqs. (62)-(65) | 285 | 310 | Gaussian-beam longitudinal temporal covariance and power spectrum. | temporal.py | 7 | P1 |
| G-71 | 8.5.3 | Ch. 8, Eq. (66) and the radial spectrum that follows | 286-287 | 311-312 | Off-axis (radial) temporal covariance and power spectrum term. | temporal.py | 7 | P2 |
| G-72 | 8.6.1 | Ch. 8, Eqs. (75), (77), (80) | 289-291 | 314-316 | Phase variance: plane-wave geometrical-optics form with the outer scale, the diffraction correction, and the Gaussian-beam generalisation. | structure.py | 6 | P2 |
| G-73 | 8.6.2 | Ch. 8, Eqs. (91)-(93) | 292 | 317 | Phase structure function with an inner and an outer scale, plus the inertial-range limit. | structure.py | 6 | P2 |
| G-74 | 8.6.3 | Ch. 8, Eqs. (102), (104), (106) | 295-296 | 320-321 | Phase spatial covariance function, plane wave and Gaussian beam, with the diffraction correction. | structure.py | 6 | P3 (no coherent receiver in olb yet) |
| G-75 | 8.6.4 | Ch. 8, Eqs. (110), (113), (114) | 297-298 | 322-323 | Longitudinal phase temporal power spectrum. | temporal.py | 7 | P2 |
| G-76 | 8.7.1 | Ch. 8, Eq. (115) | 299 | 324 | Slant-path plane-wave scintillation, 2.25 k^(7/6) L^(5/6) INT Cn2(z)(1-z/L)^(5/6) dz: the UPLINK weighting, heavy near the transmitter. | paths.py | - | P1 |
| G-77 | 8.7.1 | Ch. 8, Eq. (116) | 299 | 324 | Slant-path spherical-wave scintillation, 2.25 k^(7/6) INT Cn2(z) z^(5/6)(1-z/L)^(5/6) dz. | paths.py | - | P1 |
| G-78 | 8.7.1 | Ch. 8, Eq. (118) | 299 | 324 | Slant-path UNTRACKED Gaussian-beam total scintillation index (coefficients 14.50 and 8.70), with a variable Cn2(z). A direct book-cited replacement for the Dios integrator on the uplink. | paths.py | 2 | P1 |
| G-79 | 8.7.1 | Ch. 8, Eq. (119) | 300 | 325 | Slant-path TRACKED Gaussian-beam total scintillation index. | paths.py | 2 | P1 |
| G-80 | 8.7.1 | Ch. 8, Eqs. (121)-(123); Summary Eq. (134) | 300-301 (304) | 325-326 (329) | Slant-path irradiance covariance, Gaussian beam, plane wave and spherical wave. | paths.py | 8 | P2 |
| G-81 | 8.7.2 | Ch. 8, Eqs. (124) and (125) | 301 | 326 | Slant-path phase variance, 0.78 k^2 k0^(-5/3) INT Cn2(z) dz, von Karman outer scale. | paths.py | 6 | P2 |
| G-82 | 9.2.1 | Ch. 9, Eq. (2) | 325 | 350 | The three controlling scale sizes: coherence radius, Fresnel zone and scattering disk. A cheap regime diagnostic for `olb/assumptions.py`, better than the single 0.25 threshold olb uses now. | scintillation.py | - | P2 |
| G-83 | 9.2.3 | Ch. 9, Eqs. (9), (11), (12) | 329 | 354 | The large-scale times small-scale modulation frame itself: I = X Y, sigma_I^2 = sigma_X^2 + sigma_Y^2 + sigma_X^2 sigma_Y^2, and sigma_I^2 = exp(sigma_lnX^2 + sigma_lnY^2) - 1. Every model below hangs off this one frame, so code it once as the shared structure. | scintillation.py | 5, 10 | P1 |
| G-84 | 9.3.1, 9.3.2 | Ch. 9, Eqs. (18), (19), (23), (24) | 331-332 | 356-357 | Saturation-regime asymptotes: 1 + 0.86/sigma_R^(4/5) (plane), 1 + 2.73/sigma_R^(4/5) (spherical), 1 + (0.86 + 1.87 Theta)/sigma_R^(4/5) (beam), plus the inner-scale version. Use them as REGRESSION TESTS on any strong-regime model, the way the book uses them in Fig. 9.5. | scintillation.py | 4 | P2 |
| G-85 | 9.4.1 | Ch. 9, Eqs. (41) and (46) | 335-336 | 360-361 | The plane-wave large-scale and small-scale log-irradiance variances that FEED the gamma-gamma parameters: sigma_lnX^2 = 0.49 sigma_R^2/(1+1.11 sigma_R^(12/5))^(7/6) and sigma_lnY^2 = 0.51 sigma_R^2/(1+0.69 sigma_R^(12/5))^(5/6). olb has only their exponentiated sum, and with wrong constants (PW-09). | scintillation.py | 5 | P1 |
| G-86 | 9.4.1 | Ch. 9, Eq. (47); Ch. 12, Eqs. (40) and (93); App. III Table VII(b) | 336 (497, 522, 769) | 361 (522, 547, 794) | The CORRECT plane-wave all-regime scintillation index, with 0.49, 1.11, 0.51 and 0.69. This REPLACES the wrong constants at `plane_wave_scintillation.py:284` (PW-09) and `general_atmospherics.py:23` (KR-24). A straight constant fix, confirmed by 3 readers in 3 chapters. | scintillation.py | 10 | P1 |
| G-87 | 9.4.2 | Ch. 9, Eq. (48) | 338 | 363 | The weak-fluctuation plane-wave index on the MODIFIED spectrum, with Q_l = 10.89 L/(k l0^2). It is the input that Eq. (60) needs. | scintillation.py | 6 | P1 |
| G-88 | 9.4.2 | Ch. 9, Eqs. (55), (56), (60), (61) | 339-340 | 364-365 | Plane-wave scintillation index WITH a finite inner scale and a finite outer scale. The book shows (Fig. 9.6) that a 1 m outer scale, typical near the ground, changes the answer strongly once sigma_R > 2. A terrestrial link near the ground sits exactly there. | scintillation.py | 6 | P1 |
| G-89 | 9.5 | Ch. 9, Eqs. (63) and (64) | 341 | 366 | The spherical-wave Rytov variance beta0^2 = 0.4 sigma_R^2, and the pair sigma_R^2 = 2.5 beta0^2. olb mixes plane and spherical Rytov variances without this conversion. | scintillation.py | 10 | P2 |
| G-90 | 9.5.1 | Ch. 9, Eqs. (69) and (72) | 342 | 367 | The SPHERICAL-wave large-scale and small-scale log-irradiance variances, 0.20 sigma_R^2/(1+0.19 sigma_R^(12/5))^(7/6) and 0.20 sigma_R^2/(1+0.23 sigma_R^(12/5))^(5/6). Needed for a gamma-gamma retro or point-source link. | scintillation.py | 5 | P1 |
| G-91 | 9.5.2 | Ch. 9, Eqs. (75), (78), (80), (82), (83) | 343-345 | 368-370 | Spherical-wave scintillation index with an inner and an outer scale, plus its weak-regime input. The book validates it against 1200 m horizontal data (Sec. 9.5.3) and shows the infinite-outer-scale curve MISSES the data while a 0.6 m outer scale fits. | scintillation.py | 6 | P2 |
| G-92 | 9.6.1 | Ch. 9, Eqs. (89), (90), (91) | 350 | 375 | Beam-wander pointing-error variance and the tracked/untracked radial scintillation split. The book says an untracked collimated or divergent horizontal beam has a negligible wander effect, but a CONVERGENT beam does not. olb has no such test. | wander.py | - | P2 |
| G-93 | 9.6.2 | Ch. 9, Eqs. (97) and (101) | 352 | 377 | The GAUSSIAN-BEAM large-scale and small-scale log variances, 0.49 sigma_B^2/(1+0.56(1+Theta) sigma_B^(12/5))^(7/6) and 0.51 sigma_B^2/(1+0.69 sigma_B^(12/5))^(5/6). This closes gamma-gamma for the uplink and the terrestrial beam, where olb now has only a lognormal. | scintillation.py | 5 | P1 |
| G-94 | 9.6.2, 9.11 | Ch. 9, Eqs. (92), (93); restated as Eq. (148) | 351 (381) | 376 (406) | sigma_B^2, the weak-fluctuation Gaussian-beam Rytov variance, in exact hypergeometric form and in the simple 3.86/0.40/(11/16) approximation for a collimated or divergent beam. It is the input to G-93 and to the on-axis all-regime index. | scintillation.py | 9 | P1 |
| G-95 | 9.6.3 | Ch. 9, Eqs. (104), (105), (108), (109), (110), (112), (114) | 354-356 | 379-381 | Gaussian-beam scintillation index with an inner and an outer scale, longitudinal and radial. This is the full untracked-beam model. | scintillation.py | 6 | P2 |
| G-96 | 9.7.1, 9.11 | Ch. 9, Eqs. (119), (120); restated as Eqs. (154) and (155) | 362 (383) | 387 (408) | Plane-wave irradiance covariance in the strong regime, with eta_X = 2.61/(1+1.11 sigma_R^(12/5)) and eta_Y = 3(1+0.69 sigma_R^(12/5)), and the normalised covariance that defines the correlation width. The two-scale shape (a short core plus a long tail) is the physical basis of temporal statistics. | temporal.py | 7 | P1 |
| G-97 | 9.7.1 | Ch. 9, Eqs. (121) and (122) | 363 | 388 | Plane-wave irradiance covariance WITH an inner and an outer scale. | temporal.py | 6, 7 | P2 |
| G-98 | 9.7.2 | Ch. 9, Eqs. (123), (124), (125) | 363-364 | 388-389 | Spherical-wave irradiance covariance, zero inner scale. The book states that no Gaussian-beam covariance has been computed, so that case stays out of scope. | temporal.py | 7 | P2 |
| G-99 | 9.8, 9.11 | Ch. 9, Eqs. (126), (127), (128); restated as Eqs. (156) and (157) | 365 (383) | 390 (408) | Temporal covariance from Taylor frozen turbulence, and the irradiance power spectrum scaled by the transition frequency omega_t = V_perp/sqrt(L/k). This is the direct answer to the "temporal-statistics side-step": it gives a fade RATE and a fade DURATION, not only a fade depth. | temporal.py | 7 | P1 |
| G-100 | 9.9.1 | Ch. 9, Eqs. (131) and (132) | 368 | 393 | K distribution and its index sigma_I^2 = 1 + 2/alpha. It is the alpha = 1 or beta = 1 special case of gamma-gamma, so it is nearly free once G-102 exists. Valid only for sigma_I^2 > 1. | distributions.py | 5 | P2 |
| G-101 | 9.9.2 | Ch. 9, Eq. (133) | 369 | 394 | Lognormal-Rician (Beckmann) PDF. The book flags two problems: no closed form, and no way to tie its parameters to atmospheric conditions. So build it AFTER gamma-gamma, and only as a comparison model. | distributions.py | 5 | P2 |
| G-102 | 9.10 | Ch. 9, Eq. (137); restated as Eq. (159); Ch. 12, Eq. (67) | 370 (384, 510) | 395 (409, 535) | Gamma-gamma irradiance PDF, valid in ALL fluctuation regimes, unlike the lognormal. Closes the reserved slot DL-05. | distributions.py | 5 | P1 |
| G-103 | 9.10 | Ch. 9, Eq. (138); restated as Eq. (160); Ch. 12, Eq. (68) | 370 (384, 511) | 395 (409, 536) | The two gamma-gamma parameters straight from atmospheric physics: alpha = 1/(exp(sigma_lnX^2) - 1) and beta = 1/(exp(sigma_lnY^2) - 1). No fitted parameter. This is why the book prefers gamma-gamma. | distributions.py | 5 | P1 |
| G-104 | 9.10 | Ch. 9, Eq. (139) | 371 | 396 | Consistency identity sigma_I^2 = 1/alpha + 1/beta + 1/(alpha beta). Use it as the unit test that ties the gamma-gamma parameters back to the scintillation index. | distributions.py | 5 | P1 |
| G-105 | 9.10 | Ch. 9, Eq. (140) | 371 | 396 | Closed-form gamma-gamma cumulative distribution as a sum of two 1F2 hypergeometric functions. The book states that the CDF, not the PDF, is what a link budget needs, because it gives the fade probability directly. | distributions.py | 5 | P1 |
| G-106 | 9.11 | Ch. 9, Eq. (158) | 384 | 409 | Move the lognormal PDF out of `olb/links/downlink.py` and `olb/links/terrestrial.py` into a shared distributions module, so lognormal, K and gamma-gamma share one interface and one fade-quantile path. | distributions.py | 10 | P2 |
| G-107 | 10.2.4 | Ch. 10, Eq. (53) | 409 | 434 | Exact spherical-wave aperture-averaging factor in closed form, a 2F1 hypergeometric. A cheap analytic cross-check for the olb numerical integral (PW-06). | aperture.py | 10 | P2 |
| G-108 | 10.3 (whole book searched) | none | - | - | NEGATIVE RESULT, recorded: the book gives NO aperture filter, MTF or flux variance for an ANNULAR (centrally obscured) RECEIVE aperture. Sections 10.3.1-10.3.6 use only a soft Gaussian aperture with D_G^2 = 8 W_G^2, or the unobscured circular MTF of Eq. (54). A full-text search of all 809 PDF pages returns "annular" only in Ch. 17 (an annular TRANSMIT beam, printed 720-729) and "obscur" only there plus one Ch. 14 aside. olb gap 8 CANNOT be closed from Andrews and Phillips; it needs another source. | out of scope | 8 | P2 |
| G-109 | 10.3.2 | Ch. 10, Eqs. (60) and (61) | 412 | 437 | The book's own weak plane-wave Kolmogorov flux variance Eq. (60), and its 7 %-error aperture-averaging fit Eq. (61), A = [1 + 1.062 (k D_G^2/(4 L))]^(-7/6). olb ships a different algebraic fit that is not in the book (PW-12). | aperture.py | 10 | P2 |
| G-110 | 10.3.2 | Ch. 10, Eqs. (62)-(68) | 412-413 | 437-438 | Plane-wave aperture-averaged flux variance WITH a finite inner scale and a finite outer scale (Q_l = 10.89 L/(k l0^2), Q_0 = 64 pi^2 L/(k L0^2)). The olb Eq. (69) implementation (PW-11) is the zero-inner-scale, infinite-outer-scale limit only. | aperture.py | 6 | P2 |
| G-111 | 10.3.3 | Ch. 10, Eqs. (70)-(77) | 415-416 | 440-441 | Spherical-wave aperture-averaged flux variance, the two-scale Eq. (71) and the zero-scale Eq. (77) with 0.49 / 0.18 d^2 / 0.56 / 0.51 / 0.69 / 0.90 d^2 / 0.62 d^2. Validated against Churnside's measured data to about 15 %. olb has no spherical-wave aperture averaging. | aperture.py | 6 | P2 |
| G-112 | 10.3.5 | Ch. 10, Eq. (78) | 419 | 444 | The weak-fluctuation Gaussian-beam flux-variance double integral, which predicts that the flux variance goes to ZERO when the collecting-lens radius equals the incident beam radius. The olb plane-wave model can never reproduce this finite-beam capture effect. | aperture.py | 9 | P1 |
| G-113 | 10.3.5 | Ch. 10, Eqs. (79)-(86) | 419-420 | 444-445 | Gaussian-beam aperture-averaged flux variance WITH a finite inner scale and a finite outer scale. | aperture.py | 6 | P2 |
| G-114 | 10.3.5 | Ch. 10, Eqs. (87)-(90) | 420 | 445 | Gaussian-beam aperture-averaged flux variance for zero inner scale and infinite outer scale, in the beam parameters and Omega_G = 16 L/(k D_G^2). Valid weak to strong. olb has NO Gaussian-beam aperture averaging at all; it uses the plane-wave form for every receiver. | aperture.py | 9 | P1 |
| G-115 | 10.3.6 | Ch. 10, Eqs. (91)-(97) | 421-422 | 446-447 | Aperture-averaged TEMPORAL covariance and temporal spectrum of the irradiance under Taylor frozen turbulence, with the closed forms (94) and (95), the Fresnel frequency, and the spectrum by Eqs. (96) and (97). olb has no temporal model of any kind. | temporal.py | 7 | P1 |
| G-116 | 11.2.2 | Ch. 11, Eqs. (12) and (15) | 447-448 | 472-473 | Rice level-crossing rate: expected threshold crossings per second from the joint PDF of a process and its time derivative, and the quasi-frequency n0. The general machinery behind every fade-rate quantity. | temporal.py | 7 | P1 |
| G-117 | 11.3.1 | Ch. 11, Eqs. (23), (24), (25); Ch. 12, Eqs. (69)-(71) | 451 (511) | 476 (536) | Probability of fade for a LOGNORMAL irradiance in closed form, Pr = 0.5{1 + erf[(0.5 sigma_I^2 - 0.23 F_T)/(sigma_I sqrt2)]}, with the fade threshold parameter F_T in dB. SIGN CORRECTED by WP1: the reader read "+" from a two-column scan interleave; the minus follows from Ch. 12, Eq. (71) and the F_T = 0 check (Pr must exceed 0.5). Implemented and asserted in andrews/distributions.py. olb previously reached the same number only by a Monte Carlo percentile (RS-04). | distributions.py | 7 | P1 |
| G-118 | 11.3.1 | Ch. 11, Eqs. (21), (26), (27), (28) | 450-452 | 475-477 | Gamma-gamma irradiance PDF in the DETECTOR plane and its closed-form cumulative fade probability, with alpha and beta set from the APERTURE-AVERAGED log variances. The book states the lognormal is optimistic in the deep-fade tail, which is exactly the region olb reports. | distributions.py | 5 | P1 |
| G-119 | 11.3.2 | Ch. 11, Eqs. (33), (34), (35); Ch. 12, Eq. (72) | 455-456 (513) | 480-481 (538) | Expected number of fades per second for a lognormal irradiance, with the quasi-frequency n0 from the second derivative of the covariance. | temporal.py | 7 | P1 |
| G-120 | 11.3.2 | Ch. 11, Eqs. (36), (37), (38); Ch. 12, Eq. (74) | 456 (514) | 481 (539) | Expected number of fades for a GAMMA-GAMMA irradiance, a modified-Bessel closed form, with its own quasi-frequency. | temporal.py, distributions.py | 7 | P1 |
| G-121 | 11.3.3 | Ch. 11, Eq. (39); Ch. 12, Eqs. (78) and (79) | 456 (515) | 481 (540) | Mean fade time, the probability of fade divided by the fade rate: the average seconds the link stays below threshold. olb reports a fade DEPTH but never a fade DURATION. | temporal.py | 7 | P1 |
| G-122 | 11.2, 11.4 | Ch. 11, Eqs. (6), (9)-(11), (15)-(18), (69)-(74) | 445-464, 471-472 | 470-489, 496-497 | Receiver electrical chain: shot-noise-limited SNR, probability of detection and false alarm, threshold-to-noise ratio, false-alarm rate, mean SNR in turbulence, conditional fade probability and on-off-keying bit error rate. Out of scope: olb models the optical channel only and has no modem or detector-noise model. | out of scope | - | P3 |
| G-123 | 11.5 | Ch. 11, Secs. 11.5.1 and 11.5.2 | 465-470 | 490-495 | Spatial diversity: aperture averaging by an ARRAY of small receivers, and the resulting bit error rate. Out of scope: olb assumes a single receive aperture, and the error-rate half needs a modem model. | out of scope | - | P3 |
| G-124 | 12.2.1 | Ch. 12, Eq. (2) | 481 | 506 | Rms pseudowind computed from V(h) over 5-20 km. | paths.py | - | P2 |
| G-125 | 12.2.1 | Ch. 12, Eqs. (4) and (5) | 482 | 507 | SLC day and SLC night piecewise Cn2 profiles. | paths.py | - | P2 |
| G-126 | 12.2.2 | Ch. 12, Eq. (6) | 483 | 508 | SCIDAR outer-scale altitude profile, L0(h) = 4/[1 + ((h-8500)/2500)^2], capped at 4 m. | paths.py | 6 | P1 |
| G-127 | 12.2.2 | Ch. 12, Eq. (7) | 483 | 508 | Second SCIDAR outer-scale profile, L0(h) = 5/[1 + ((h-7500)/2500)^2], capped at 5 m. | paths.py | 6 | P1 |
| G-128 | 12.3.3 | Ch. 12, text and Fig. 12.4 | 488 | 513 | Point-ahead angle theta_p = 2V/c, and the statement that it usually exceeds the isoplanatic angle. | paths.py | 2 | P2 |
| G-129 | 12.4.1 | Ch. 12, Eqs. (17)-(19) | 491 | 516 | Downlink Gaussian-beam wave structure function, before the plane-wave reduction. | structure.py | - | P2 |
| G-130 | 12.4.1 | Ch. 12, Eqs. (24)-(27) | 492 | 517 | Uplink wave structure function and uplink spatial coherence radius, with the Lambda^(11/6) term. This is the book's own uplink path weighting, and it is what settles the mirrored-weight question of GF-18 and KR-01. | structure.py | 3 | P1 |
| G-131 | 12.4.3 | Ch. 12, Eq. (29) | 493 | 518 | Isoplanatic angle for an upward-propagating Gaussian-beam wave, not the spherical-wave limit. | paths.py | - | P2 |
| G-132 | 12.5.1 | Ch. 12, Eq. (33) | 494 | 519 | Downlink long-term spot size, for low-elevation paths where the second path moment cannot be dropped. | beam.py | - | P2 |
| G-133 | 12.5.2 | Ch. 12, Eqs. (36) and (37) | 495 | 520 | Downlink Gaussian-beam scintillation index including the off-axis radial term. | scintillation.py | - | P2 |
| G-134 | 12.5.2 | Ch. 12, Eq. (39) | 496 | 521 | Downlink irradiance flux variance for a HARD aperture D_G, in closed form by the ABCD route. A direct analytic replacement for the olb numerical aperture integral (PW-05). | aperture.py | 8 | P1 |
| G-135 | 12.5.3 | Ch. 12, Eqs. (41)-(43) | 498 | 523 | Downlink irradiance covariance and the correlation width rc = sqrt(45e3 sec(zeta)/k), which sets when aperture averaging starts. | scintillation.py | 8 | P2 |
| G-136 | 12.5.3 | Ch. 12, Eqs. (44)-(47) | 498-499 | 523-524 | Downlink large-scale and small-scale log variances with eta_X and eta_Y: the direct inputs to gamma-gamma for the slant path. | distributions.py | 5 | P1 |
| G-137 | 12.6.1 | Ch. 12, Eqs. (48) and (49) | 500-501 | 525-526 | Uplink long-term spot size and Strehl ratio by reciprocity, with D0^2 = 8 W0^2, both branches of D0/r0. Settles the approximate row UF-09. | beam.py | - | P2 |
| G-138 | 12.6.2 | Ch. 12, Eq. (50) | 502 | 527 | Uplink beam-wander displacement variance with the explicit outer-scale wavenumber kappa0(h). | wander.py | 6 | P1 |
| G-139 | 12.6.2 | Ch. 12, Eqs. (51) and (52) | 502 | 527 | Simplified wander variance, 0.54 (H-h0)^2 sec^2(zeta)(lambda/(2W0))^2 (2W0/r0)^(5/3)[1 - (kappa0^2 W0^2/(1+kappa0^2 W0^2))^(1/6)]. | wander.py | 6 | P1 |
| G-140 | 12.6.3 | Ch. 12, Eq. (53) | 503 | 528 | Uplink pointing-error variance with the scaling constant Cr. This is the book's own home for the jitter fold of UF-08, and the book form that PT-02 lacks. | wander.py | - | P1 |
| G-141 | 12.6.3 | Ch. 12, Eq. (54); repeated as Eq. (100) | 503 (524) | 528 (549) | Uplink UNTRACKED longitudinal scintillation index, 5.95 (H-h0)^2 sec^2(zeta)(2W0/r0)^(5/3)(alpha_pe/W)^2 plus the Rytov term. | scintillation.py | 2 | P1 |
| G-142 | 12.6.3 | Ch. 12, Eq. (56) | 504 | 529 | Uplink untracked OFF-AXIS scintillation index, with the unit step at a radial angle beyond the pointing error. | scintillation.py | 2 | P1 |
| G-143 | 12.6.3 | Ch. 12, Eqs. (57) and (58) | 504 | 529 | Uplink TRACKED scintillation index and its Rytov variance. A tilt-tracked uplink still scintillates. NOTE (2026-08-27): this is NOT the pre-compensated residual — that case has no analytic form, by decision (backlog 0-W1). | scintillation.py | 2 | P1 |
| G-144 | 12.6.4 | Ch. 12, Eqs. (59)-(61); repeated as Eqs. (99) and (100) | 506 (524) | 531 (549) | Uplink strong-fluctuation scintillation index, tracked and untracked, valid for all beam Rytov variances. | scintillation.py | 4 | P1 |
| G-145 | 12.6.4 | Ch. 12, Eq. (62) | 509 | 534 | Uplink spherical-wave weak scintillation index, 2.25 k^(7/6)(H-h0)^(5/6) sec^(11/6)(zeta) INT Cn2 (1-x)^(5/6) x^(5/6) dh. | scintillation.py | - | P2 |
| G-146 | 12.6.5 | Ch. 12, Eqs. (63) and (64) | 509-510 | 534-535 | Uplink irradiance covariance function. Out of scope: Andrews states the uplink correlation width is tens of metres, so any satellite receiver is a point receiver and aperture averaging never applies. | out of scope | - | P3 |
| G-147 | 12.7 | Ch. 12, Eqs. (65) and (66) | 510 | 535 | Lognormal irradiance PDF for a Gaussian-beam wave, with the radial mean-irradiance roll-off exp(-2 r^2/W_LT^2) inside the log. Closes the unmatched kernel row KR-22. | distributions.py | - | P2 |
| G-148 | 12.7.1 | Ch. 12, text after Eq. (71) and the Fig. 12.20 discussion | 513 | 538 | Aperture-averaged gamma-gamma parameters, alpha = [0.49 sigma_I^2(D_G)]^(-1) and beta = [0.51 sigma_I^2(D_G)]^(-1). | distributions.py | 5 | P1 |
| G-149 | 12.7.2 | Ch. 12, Eqs. (73), (75), (76), (77) | 514 | 539 | Quasi-frequency n0, the temporal irradiance covariance under Taylor frozen flow, and the second derivative at zero lag with the wind-weighted path moment. | temporal.py | 7 | P1 |
| G-150 | 12.8.1 | Ch. 12, Eqs. (80) and (81) | 517 | 542 | UPLINK gamma-gamma parameters: the large-scale variance folds the pointing-error term in, and the small-scale variance is the rest. | distributions.py | 5 | P1 |
| G-151 | 12.8.3 | Ch. 12, Eqs. (82)-(84) | 519 | 544 | Uplink second derivative of the temporal covariance, with the two uplink path moments. | temporal.py | 7 | P1 |
| G-152 | 13.7.1 | Ch. 13, Eqs. (132) and (133) | 577 | 602 | Backscatter amplification factor, weak fluctuations. Out of scope: olb models retro as a retransmission and no user need is stated. | out of scope | - | P3 |
| G-153 | 13.7.2 | Ch. 13, Eqs. (134)-(142) | 578-580 | 603-605 | Backscatter amplification factor, strong fluctuations. Out of scope: same reason. | out of scope | - | P3 |
| G-154 | 13.7.3 | Ch. 13, Eq. (146) | 581 | 606 | Point-target spatial coherence radius of the reflected wave. Out of scope: it already agrees with the one-way spherical wave that olb assumes. | out of scope | - | P3 |
| G-155 | 13.7.4 | Ch. 13, Eqs. (147)-(159) | 581-584 | 606-609 | Point-target covariance and scintillation index, weak fluctuations, including the monostatic-against-bistatic correlation term. Out of scope: same reason. It is the book evidence behind the unmatched row RT-02. | out of scope | - | P3 |
| G-156 | 13.7.5 | Ch. 13, Eqs. (161) onward | 585 | 610 | Point-target scintillation index, strong fluctuations. Out of scope: same reason. | out of scope | - | P3 |
| G-157 | 14.3.2 | Ch. 14, Eqs. (26)-(28) | 617-618 | 642-643 | Long-term turbulence MTF, exp[-3.44 (lambda F n/r0)^(5/3)], and the total system MTF for a Gaussian and for a hard aperture. | aperture.py | - | P2 |
| G-158 | 14.3.3 | Ch. 14, Eqs. (33)-(35) | 621 | 646 | Short-term (tilt-removed) point-spread function and MTF, 0.28 (D_G/r0)^(5/3), with the plane-wave and spherical-wave branches. This is the Andrews route to a tip-tilt-corrected Strehl ratio. | aperture.py | 1 | P2 |
| G-159 | 14.3.4 | Ch. 14, Eq. (38) | 622 | 647 | Greenwood time constant, tau0 = [2.91 k^2 INT Cn2(z) V^(5/3)(z) dz]^(-3/5), the wind-weighted profile integral. | temporal.py | 7 | P1 |
| G-160 | 14.3.4 | Ch. 14, Eq. (39) | 622 | 647 | Greenwood time constant for a constant wind, tau0 = 0.32 r0/V_perp, and the Greenwood frequency. | temporal.py | 7 | P1 |
| G-161 | 14.3.5 | Ch. 14, Eqs. (40), (42), (43); Ch. 10, Eqs. (48)-(50) | 623 (408) | 648 (433) | Strehl ratio against D/r0: the weak form 1/[1+(D/r0)^(5/3)], the all-conditions form [1+(D/r0)^(5/3)]^(-6/5), and the Sasiela asymptotic series (r0/D)^2 + 0.6159 (r0/D)^3 + 0.0500 (r0/D)^5. olb has no Strehl-ratio output. | aperture.py | - | P1 |
| G-162 | 14.5.3 | Ch. 14, Eqs. (85) and (86) | 634 | 659 | GENERAL Zernike aperture filter functions for any azimuthal and radial order, not only piston and tilt. olb hard-codes orders 0 and 1 (AN-01, AN-02). | aperture.py | 8 | P1 |
| G-163 | 14.5.3 | Ch. 14, Eq. (87) | 634 | 659 | Gaussian-beam version of the Zernike filters, by the replacement kappa -> gamma kappa D_G/2 with gamma = 1 - (Theta + i Lambda)(1 - z/L). Turns the plane-wave anisoplanatism and the Noll variance into beam-wave forms. | aperture.py | 3 | P1 |
| G-164 | 14.5.4 | Ch. 14, Eq. (88) | 635 | 660 | Geometrical-optics phase variance, 1.83 (L0/r0)^(5/3), under the von Karman spectrum. It shows that the UNFILTERED phase variance diverges without an outer scale. | structure.py | 6 | P2 |
| G-165 | 14.5.4 | Ch. 14, Eq. (93) | 635 | 660 | Tilt-only aperture-averaged phase variance, 0.90 (D_G/r0)^(5/3). olb has the piston-removed and the piston-plus-tilt-removed variances, but not tilt alone. | aperture.py | 1 | P2 |
| G-166 | App. III, Table V | Table V, three spectrum rows | 767 | 792 | Spherical-wave coherence radius, rho_sp = (0.55 Cn2 k^2 L)^(-3/5), plus the inner-scale branches 0.55 and 0.62. | structure.py | 6 | P2 |

## Table 3 — constants ledger

Every literal physics number, deduplicated by (value, location) and sorted by
file then line. The last block holds book constants that olb does NOT have; it
is the numeric face of Table 2.

| value | location | role | book equation | printed p | pdf p | agrees |
|---|---|---|---|---|---|---|
| 1e-9 | olb/beam.py:50 | numerical tolerance on the diffraction-limit test | not physics | - | - | not found (no book equation applies) |
| pi | olb/beam.py:59 | Gaussian far-field divergence, theta = lambda/(pi w) = 2/(k w) | Ch. 4, Eq. (37), far-field limit | 93 | 118 | yes |
| 0.5 | olb/beam.py:60 (kernel zR = 0.5 k w^2) | Rayleigh range | Ch. 4, Sec. 4.5.2 text | 98 | 123 | yes |
| 5/ln10 | olb/links/downlink.py:68 | dB form of the mean lognormal loss | Ch. 9, Eq. (158); Ch. 5, Eq. (93) | 384, 156 | 409, 181 | yes (a dB conversion, not a book constant) |
| 10/ln10 | olb/links/downlink.py:73 | dB form of the lognormal quantile | Ch. 9, Eq. (158); Ch. 5, Eq. (93) | 384, 156 | 409, 181 | yes |
| 0.5 (in -sigma_l^2/2) | olb/links/downlink.py:73,76 | E[I] = 1 normalisation of the lognormal | Ch. 9, Eq. (158) | 384 | 409 | yes. R2 could not find it in Ch. 4 or Ch. 5, because Ch. 5, Eq. (93) leaves the mean log-amplitude free; R5 found it explicitly in Ch. 9, Eq. (158). |
| 0.8145 | olb/links/downlink.py:452 (self-check) | single-mode-fibre static mode-match efficiency | not a book quantity (coupling) | - | - | not found |
| 5/ln10 | olb/links/terrestrial.py:155 | dB form of the mean lognormal loss | Ch. 5, Eq. (93) | 156 | 181 | yes |
| 10/ln10 | olb/links/terrestrial.py:160 | dB form of the lognormal quantile | Ch. 5, Eq. (93) | 156 | 181 | yes |
| 0.5 (in -sigma_l^2/2) | olb/links/terrestrial.py:160,163 | E[I] = 1 normalisation | Ch. 9, Eq. (158) | 384 | 409 | yes |
| 2 (in exp(-2 r^2/w_z^2)) | olb/models/pointing.py:10 | the 1/e^2 Gaussian beam convention | Ch. 11, Eqs. (40) and (41) | 459 | 484 | yes |
| 20/ln10 = 8.6859 | olb/models/pointing.py:27 | dB per (r^2/w_z^2) | Ch. 11, Eq. (25); Eq. (26) uses 0.23 = ln(10)/10 | 451-452 | 476-477 | yes (unit convention only) |
| 2 (in 2 sigma_r^2/w_z^2) | olb/models/pointing.py:54 | two independent jitter axes, exponential mean of r^2 | none | - | - | not found. Ch. 11 explicitly neglects pointing error (printed 451). Ch. 12, Eq. (53) printed 503 is the book's own pointing-error construction, and it is different. |
| 0.0096932 | olb/turbulence/anisoplanatism.py:80 | Stone structure-function constant | Stone 1994 Eq. (14), DOI 10.1364/JOSAA.11.000347 | - | - | not found alone; its product with the other two factors reproduces the book 2.91 |
| 1.11833 | olb/turbulence/anisoplanatism.py:84 | the HJ1(8/3, 0, 1) magnitude | Stone 1994 note 15, DOI 10.1364/JOSAA.11.000347 | - | - | not found alone |
| (2 pi)^(8/3) | olb/turbulence/anisoplanatism.py:87 | spatial-frequency scale of Stone Eq. (29) | Stone 1994, DOI 10.1364/JOSAA.11.000347 | - | - | not found alone |
| 4 and (n+1)^2 | olb/turbulence/anisoplanatism.py:117 | Zernike radial-order weight p_n(u) | Ch. 14, Eq. (86) | 634 | 659 | yes; the Andrews per-mode filter summed over the n+1 azimuthal modes of order n gives exactly 4 (n+1)^2 [J_{n+1}/kappa]^2 |
| 8/3 (in u^(-8/3)) | olb/turbulence/anisoplanatism.py:189 | exponent of the Stone spatial-frequency integral | Stone 1994 Eq. (36), DOI 10.1364/JOSAA.11.000347 | - | - | not found |
| (n+1)(n+2)/2 | olb/turbulence/anisoplanatism.py:225 | Noll mode count through radial order n | Ch. 14, Table 14.1 | 630 | 655 | yes |
| 5/3 | olb/turbulence/anisoplanatism.py:266, :392 | the h^(5/3) weight and the (theta/theta0)^(5/3) law | Ch. 12, Eq. (30); Ch. 14, Eqs. (36) and (37) | 493, 622 | 518, 647 | yes |
| 2.914381 | olb/turbulence/anisoplanatism.py:267 | isoplanatic-angle constant | Ch. 12, Eq. (30); Ch. 14, Eqs. (37) and (97) | 493, 622, 636 | 518, 647, 661 | yes; the book prints the rounded 2.91 |
| 1.0299 | olb/turbulence/ao.py:39 | piston-removed phase-variance coefficient | Ch. 14, Eq. (90) | 635 | 660 | yes; the book gives 1.03, and 1.02 by the ABCD route (Eq. (91), same page) |
| 0.134 | olb/turbulence/ao.py:40 | piston-and-tilt-removed coefficient | Ch. 14, Eq. (94) | 636 | 661 | yes; the book gives 0.13 |
| 0.2944 | olb/turbulence/ao.py:41 | large-J adaptive-optics residual prefactor | Noll 1976, DOI 10.1364/JOSA.66.000207 | - | - | not found; Ch. 14.5.4 stops at piston and tilt |
| sqrt(3)/2 | olb/turbulence/ao.py:42 | large-J residual exponent | Noll 1976, DOI 10.1364/JOSA.66.000207 | - | - | not found |
| sec(zeta) | olb/turbulence/ao.py:67 | airmass on the Cn2 integral | Ch. 12, Eq. (23) | 492 | 517 | yes; Ch. 6 Sec. 6.8 gives no slant scale factor, so R3 could not confirm it there |
| 0.423 | olb/turbulence/ao.py:69 | plane-wave Fried-parameter prefactor | Ch. 12, Eq. (23); Ch. 6, Eq. (113); Ch. 14, Eqs. (25) and (95) | 492, 209, 617, 636 | 517, 234, 642, 661 | yes; the book prints the rounded 0.42, and 1.46/2.1^(5/3) = 0.4240 |
| 2.1 (implicit at ao.py:69) | olb/turbulence/ao.py:69 | ratio r0 = 2.1 rho0 | Ch. 12, Eq. (89); Ch. 6, text below Eq. (64); Ch. 14, text at Eq. (24) | 522, 194, 617 | 547, 219, 642 | yes |
| 2.0 (in f_c = sqrt(N)/(2 D)) | olb/turbulence/ao.py:147 | adaptive-optics cutoff spatial frequency | none | - | - | not found; a heuristic with no Andrews counterpart |
| 0.023 | olb/turbulence/ao.py:151 | Kolmogorov residual phase PSD prefactor | Noll 1976, DOI 10.1364/JOSA.66.000207 | - | - | not found. CITATION FAULT: the docstring credits "Andrews Ch. 3". The nearest Andrews form, Ch. 14, Eq. (88) printed 635, is a different von Karman geometrical-optics statement. |
| 5/3 | olb/turbulence/ao.py:178 | the (D/r0)^(5/3) scaling of every residual variance | Ch. 14, Eqs. (90), (93), (94) | 635-636 | 660-661 | yes |
| 0.6 (docstring) | olb/turbulence/beam_wave_scintillation.py:25 | Dios split-step agreement ceiling on sigma2_chi | Ch. 8, text after Eq. (23) | 264-265 | 289-290 | no; Andrews puts the weak limit at sigma_I^2 < 1, that is sigma_chi^2 < 0.25 by Eq. (13). 0.6 is 2.4 times looser and is a Dios empirical bound. |
| 0.033 | olb/turbulence/beam_wave_scintillation.py:55 | Kolmogorov spectrum leading constant | Ch. 3, Eq. (18); Ch. 6, Eq. (110); Ch. 8, Eqs. (25) and (120); Ch. 12, Eq. (15) | 67, 208, 265, 299, 490 | 92, 233, 290, 324, 515 | yes |
| Gamma(-5/6) = -6.6865 | olb/turbulence/beam_wave_scintillation.py:56 | Kolmogorov kappa-integration factor in the Dios prefactor | not a literal in the book | 264 | 289 | yes, indirectly: the combination 4 pi^2 k^2 Gamma(-5/6) 0.033 reproduces the Andrews Eq. (20) plane-wave 1.23 and spherical-wave 0.4 limits to 2-3 % in the module self-check |
| 1 (Theta0 for f0 = inf) | olb/turbulence/beam_wave_scintillation.py:84 | collimated-beam curvature parameter | Ch. 8, Eq. (5); App. III Table VI footer | 261, 768 | 286, 793 | yes for the collimated case. The argument is threaded here, so this call site is not collimated-only. |
| 2.0 | olb/turbulence/beam_wave_scintillation.py:85, :89, :171 | Lambda0 = 2L/(k W0^2); W^2 = 2L/(k Lambda); the 1F1 argument 2 r^2/W^2 | Ch. 8, Eqs. (5), (6), (18) | 261, 263 | 286, 288 | yes |
| 5/6 | olb/turbulence/beam_wave_scintillation.py:137, :169 | Kolmogorov exponent on A and on Lambda | Ch. 8, Eqs. (18), (19), (22), (23) | 263-264 | 288-289 | yes |
| 5/12 | olb/turbulence/beam_wave_scintillation.py:138 | exponent on (A^2 + B^2), half of 5/6 | Ch. 8, Eq. (23) | 264 | 289 | yes |
| 4 pi^2 | olb/turbulence/beam_wave_scintillation.py:140, :170 | prefactor of the scintillation path integral | Ch. 8, Eq. (14) uses 8 pi^2 k^2 L | 262 | 287 | yes; the book prefactor is over the normalised path variable, the Dios form is over z with the Gamma factor already applied |
| -5/6 | olb/turbulence/beam_wave_scintillation.py:172 | first parameter of 1F1(-5/6, 1, 2 r^2/W^2) | Ch. 8, Eq. (18) | 263 | 288 | yes |
| 1.23 | olb/turbulence/beam_wave_scintillation.py:208 | plane-wave Rytov variance constant | Ch. 8, Eq. (20) | 264 | 289 | yes |
| 0.404 | olb/turbulence/beam_wave_scintillation.py:216 | spherical-wave limit factor on sigma_R^2 | Ch. 8, Eq. (20); Ch. 9, Eqs. (63) and (64) | 264, 341 | 289, 366 | yes; the book value is 0.4 exactly and the self-check tolerance is 3 % |
| 4.42 | olb/turbulence/beam_wave_scintillation.py:231 | small-r radial scintillation slope | Ch. 8, Eq. (22); Ch. 9, Eqs. (87), (88), (113) | 264, 349-356 | 289, 374-381 | yes |
| 0.00594, 27.0, 1e-5, power 10, 1000, 2.7e-16, 1500, 1.7e-14, 100, 21.0 | olb/turbulence/beam_wave_scintillation.py:246-249 | Hufnagel-Valley Cn2(h) in the self-check | Ch. 12, Eq. (1) and the text after Eq. (3) | 481 | 506 | yes; R7 checked every constant term by term. R4 could not, because Ch. 8 has no Cn2 profile. |
| 1.0 (COLLIMATED_THETA0) | olb/turbulence/gaussian_fried.py:32 | input curvature parameter Theta0 | Ch. 6, Eq. (6); Ch. 4, Eq. (33) with F0 to infinity; Ch. 9 Sec. 9.12 | 183, 92, 384 | 208, 117, 409 | yes for a collimated beam ONLY. The book keeps Theta0 general; olb fixes it. This is olb gap 3. |
| 2.0 | olb/turbulence/gaussian_fried.py:49 | numerator of Lambda0 = 2 z/(k w0^2) | Ch. 6, Eq. (6); Ch. 4, Eq. (33) | 183, 92 | 208, 117 | yes |
| 1.23 | olb/turbulence/gaussian_fried.py:78 | Rytov-variance coefficient | Ch. 6, Eq. (119); Ch. 5, Eq. (15) | 210, 140 | 235, 165 | yes |
| 7/6, 11/6 | olb/turbulence/gaussian_fried.py:78-79 | Rytov-variance exponents on k and z | Ch. 6, Eq. (119); Ch. 5, Eq. (15) | 210, 140 | 235, 165 | yes |
| 12/5 | olb/turbulence/gaussian_fried.py:96 | exponent on sigma_R in q | Ch. 7, Eq. (58); Ch. 9, Eqs. (85) and (86) | 242, 349 | 267, 374 | yes |
| 1.63 | olb/turbulence/gaussian_fried.py:97 | Theta_e and Lambda_e denominator, 4q/3 | Ch. 7, Eq. (58); Ch. 9, Eqs. (85), (86), (150) | 242, 349, 382 | 267, 374, 407 | yes |
| 0.81 | olb/turbulence/gaussian_fried.py:98 | Theta_e numerator, 2q/3 with q = 1.22 sigma_R^(12/5) | Ch. 7, Eq. (58); Ch. 9, Eq. (85) | 242, 349 | 267, 374 | yes |
| 8/3 | olb/turbulence/gaussian_fried.py:113 | exponent in the a-factor | Ch. 6, Eq. (55); Ch. 7, Eq. (60) | 192, 243 | 217, 268 | yes |
| 8/3 (as 8.0/3.0) | olb/turbulence/gaussian_fried.py:129 | numerator 8 over 3 in rho0_e | Ch. 7, Eq. (59) lower; Ch. 6, Eq. (79) | 243, 199 | 268, 224 | yes |
| 0.62 | olb/turbulence/gaussian_fried.py:129 | Lambda_e^(11/6) coefficient in rho0_e | Ch. 7, Eq. (59); Ch. 6, Eqs. (74) and (79) | 243, 196, 199 | 268, 221, 224 | yes; the Ch. 6 summary Eq. (132) printed 213 prints 0.618 |
| 11/6 | olb/turbulence/gaussian_fried.py:129 | exponent on Lambda_e | Ch. 7, Eq. (59) | 243 | 268 | yes |
| 1.46 | olb/turbulence/gaussian_fried.py:141 | plane-wave coherence radius coefficient | Ch. 6, Eqs. (64) and (130); App. III Table IV | 194, 213, 767 | 219, 238, 792 | yes |
| -3/5, 3/5 | olb/turbulence/gaussian_fried.py:142, :163, :323 | coherence-radius exponent | Ch. 6, Eq. (64); Ch. 9 Sec. 9.12 | 194, 384 | 219, 409 | yes |
| 2.1 | olb/turbulence/gaussian_fried.py:155, :207, :324 | Fried factor r0 = 2.1 rho0 | Ch. 6, text below Eq. (64); App. III note under Table IV | 194, 767 | 219, 792 | yes |
| (8/3)^(3/5) = 1.7963 | olb/turbulence/gaussian_fried.py:163 | spherical-over-plane Fried ratio | Ch. 6, Eq. (71) against Eq. (64); Ch. 9 Sec. 9.12 | 196, 194, 384 | 221, 219, 409 | DISPUTED, see Conflicts C-07. R3: the book 0.55/1.46 gives 1.7913, a 0.3 % difference, because 0.55 is a rounded 0.5475. R5: the Ch. 9 worked example r0 = (0.16 Cn2 k^2 L)^(-3/5) confirms it, because 0.423 x 3/8 = 0.1586 rounds to 0.16. |
| 0.42 | olb/turbulence/gaussian_fried.py:174 (docstring, cites Dios Eq. (3)) | spherical-wave Fried coefficient | Ch. 6, Eqs. (71) and (116), times 2.1 | 196, 209 | 221, 234 | yes; 1.46/2.1^(5/3) = 0.4240 |
| 5/3 | olb/turbulence/gaussian_fried.py:319, :320 | wave-structure-function path-weight exponent | Ch. 6, Eqs. (112), (115), (116); Ch. 9 Sec. 9.3.1 text | 209, 330 | 234, 355 | yes for the exponent, but see GF-18: olb mirrors the weight DIRECTION |
| 1.46 | olb/turbulence/gaussian_fried.py:322 | profile coherence coefficient | Ch. 6, Eq. (115) | 209 | 234 | yes |
| 0.62 | olb/turbulence/gaussian_fried.py:323 | mu2 coefficient in the profile form | Ch. 6, Eq. (115) | 209 | 234 | yes |
| 0.423 | olb/turbulence/gaussian_fried.py:358, :359, :403 (self-check) | plane-wave Fried reference | Ch. 6, Eq. (64), times 2.1 | 194 | 219 | yes |
| 3/8 | olb/turbulence/gaussian_fried.py:360 (self-check) | spherical path weight for a constant Cn2 | Ch. 6, Eq. (71) | 196 | 221 | yes |
| 0.25 | olb/turbulence/plane_wave_scintillation.py:45 | WEAK_FLUCTUATION_LIMIT on sigma2_I | Ch. 5, Eq. (15) and the text after it; Ch. 10, Eq. (61) condition; Ch. 12, Eq. (40) | 140, 412, 497 | 165, 437, 522 | NO. Confirmed by 4 readers. The book states weak is sigma_R^2 < 1. The value 0.25 is nowhere in the book. The code errs safe (R6 notes Ch. 11.3 printed 451 supports a tighter bound), but the comment cites the book for a number the book does not give. See Conflicts C-05. |
| 1e-2, 1e4, 2000 | olb/turbulence/plane_wave_scintillation.py:89 | kappa integration grid limits and point count | - | - | - | not found (a numerical choice). The grid has no inner-scale or outer-scale cut, which matches the Ch. 3, Eq. (18) warning that this can make integrals diverge. |
| 5/6 | olb/turbulence/plane_wave_scintillation.py:79 | Cn2(h) h^(5/6) slant-path weight | Ch. 12, Eq. (38) | 495 | 520 | yes |
| 2.25 | olb/turbulence/plane_wave_scintillation.py:80 | slant plane-wave scintillation-index prefactor | Ch. 12, Eq. (38); repeated as Eq. (92) | 495, 522 | 520, 547 | yes (but the code cites Eq. (44), which is wrong; see PW-02) |
| 7/6, 11/6 | olb/turbulence/plane_wave_scintillation.py:80 | wavenumber and airmass exponents | Ch. 12, Eq. (38) | 495 | 520 | yes |
| 2 (in x = kappa D/2) | olb/turbulence/plane_wave_scintillation.py:105-108 | half-diameter in the aperture filter argument | Ch. 10, Eq. (54); Ch. 14, Eq. (86) with m = n = 0 | 410, 634 | 435, 659 | yes |
| 0.033 | olb/turbulence/plane_wave_scintillation.py:138 | Kolmogorov spectrum constant | Ch. 3, Eq. (18); Ch. 9, Eq. (3); Ch. 12, Eq. (15) | 67, 327, 490 | 92, 352, 515 | yes |
| 11/3 | olb/turbulence/plane_wave_scintillation.py:138 | Kolmogorov spectral slope | Ch. 3, Eq. (18); Ch. 9, Eq. (3); Ch. 12, Eq. (15) | 67, 327, 490 | 92, 352, 515 | yes |
| 8 pi^2 | olb/turbulence/plane_wave_scintillation.py:149 | leading constant of the Rytov double integral | Ch. 9, Eq. (25); Ch. 10, Eq. (59); Ch. 12, Eqs. (16) and (75) | 333, 412, 491, 514 | 358, 437, 516, 539 | yes |
| 1.23 | olb/turbulence/plane_wave_scintillation.py:251 | Rytov-variance coefficient | Ch. 5, Eq. (15); Ch. 6, Eq. (119); Ch. 10, text at Eq. (60); App. III Table VII(a) footer | 140, 210, 412, 769 | 165, 235, 437, 794 | yes |
| 7/6, 11/6 | olb/turbulence/plane_wave_scintillation.py:251-252 | Rytov-variance exponents | Ch. 5, Eq. (15); Ch. 10, Eq. (60) text | 140, 412 | 165, 437 | yes |
| 1.46 | olb/turbulence/plane_wave_scintillation.py:264 | plane-wave coherence radius constant | Ch. 6, Eqs. (64) and (130); Ch. 12, Eq. (27); App. III Table IV | 194, 213, 492, 767 | 219, 238, 517, 792 | yes; the Ch. 12 DOWNLINK form Eq. (22) printed 491 uses 1.45 for the same structure |
| 12/5 | olb/turbulence/plane_wave_scintillation.py:283, :317 | exponent on sigma_R in both denominators | Ch. 9, Eqs. (40), (45), (47); Ch. 12, Eq. (40) | 335-336, 497 | 360-361, 522 | yes |
| 0.54 | olb/turbulence/plane_wave_scintillation.py:284 | large-scale term of the point plane-wave index | Ch. 9, Eq. (47); Ch. 12, Eqs. (40) and (93); App. III Table VII(b) | 336, 497, 522, 769 | 361, 522, 547, 794 | NO. Confirmed by 3 readers: the book gives 0.49. A full-text search of the book finds no 0.54 scintillation constant. With 0.509 the weak limit becomes 1.049 sigma_R^2 instead of 1.000. |
| 1.22 | olb/turbulence/plane_wave_scintillation.py:284 | large-scale denominator of the point plane-wave index | Ch. 9, Eq. (47); Ch. 12, Eqs. (40) and (93) | 336, 497, 522 | 361, 522, 547 | NO. The book gives 1.11. The value 1.22 belongs one level down, in q = L/(k rho0^2) = 1.22 sigma_R^(12/5) printed 342, which enters eta_X as 0.35 x 1.22/0.38 = 1.11. |
| 0.509 | olb/turbulence/plane_wave_scintillation.py:285 | small-scale term of the point plane-wave index | Ch. 9, Eq. (47); Ch. 12, Eqs. (40) and (93); App. III Table VII(b) | 336, 497, 522, 769 | 361, 522, 547, 794 | NO. The book gives 0.51, and the same file uses 0.51 at line 320. |
| 0.69 | olb/turbulence/plane_wave_scintillation.py:285, :320 | small-scale denominator constant | Ch. 9, Eqs. (45)-(47); Ch. 10, Eq. (69); Ch. 12, Eq. (40) | 336, 413, 497 | 361, 438, 522 | yes; eta_Y = 3(1 + 0.69 sigma_R^(12/5)) |
| 4 (in d^2 = k D^2/(4 L)) | olb/turbulence/plane_wave_scintillation.py:297 | aperture-parameter denominator | Ch. 10, Eq. (68) | 413 | 438 | yes; note Ch. 12, Eq. (39) printed 496 uses k D_G^2/(16 L) for its own variable |
| 0.49 | olb/turbulence/plane_wave_scintillation.py:319 | large-scale numerator of the aperture-averaged index | Ch. 10, Eq. (69); Ch. 9, Eq. (47) | 413, 336 | 438, 361 | yes |
| 0.65 | olb/turbulence/plane_wave_scintillation.py:319 | d^2 term in the large-scale denominator | Ch. 10, Eq. (69) | 413 | 438 | yes (R6 verified it in Ch. 10; R7's text search of PDF 400-440 did not return it) |
| 1.11 | olb/turbulence/plane_wave_scintillation.py:319 | sigma_R^(12/5) coefficient in the large-scale denominator | Ch. 10, Eq. (69); Ch. 9, Eqs. (40) and (47) | 413, 335-336 | 438, 360-361 | yes |
| 0.51 | olb/turbulence/plane_wave_scintillation.py:320 | small-scale numerator of the aperture-averaged index | Ch. 10, Eq. (69); Ch. 9, Eqs. (46) and (47) | 413, 336 | 438, 361 | yes |
| 0.90 | olb/turbulence/plane_wave_scintillation.py:321 | d^2 term in the small-scale denominator | Ch. 10, Eq. (69) | 413 | 438 | yes |
| 0.62 | olb/turbulence/plane_wave_scintillation.py:321 | d^2 sigma_R^(12/5) cross term | Ch. 10, Eq. (69) | 413 | 438 | yes |
| 7/6, 5/6 | olb/turbulence/plane_wave_scintillation.py:319-321 | exponents of the aperture-averaged index | Ch. 10, Eq. (69) | 413 | 438 | yes |
| 1.07 | olb/turbulence/plane_wave_scintillation.py:336 | weak aperture-averaging fit constant | Ch. 10, Eq. (61) uses 1.062 in a DIFFERENT algebraic form | 412 | 437 | NO. 1.07 is nowhere in the 809-page book. Eq. (61) reads A = [1 + 1.062 (k D_G^2/(4 L))]^(-7/6), with the 7/6 outside the bracket. The two fits differ by up to 12 %. The source is Churnside 1991, DOI 10.1364/AO.30.001982. |
| 7/6 (placement in the 1.07 fit) | olb/turbulence/plane_wave_scintillation.py:336 | exponent placement | Ch. 10, Eq. (61) | 412 | 437 | no; the book applies -7/6 to the whole bracket, olb applies +7/6 to d^2 inside and -1 outside |
| 2.21 | olb/turbulence/plane_wave_scintillation.py:353 | large-inner-scale aperture-averaging constant | not in the book | - | - | NO; Churnside 1991, DOI 10.1364/AO.30.001982, cited by the book as Ch. 10 reference [12] printed 438 but never reproduced |
| 7/3 | olb/turbulence/plane_wave_scintillation.py:353, :379 | large-aperture asymptote exponent | not in the book | - | - | not found; the same Churnside source |
| 0.908 | olb/turbulence/plane_wave_scintillation.py:377 | strong-regime coherence-scale term | not in the book | - | - | NO; Churnside 1991, DOI 10.1364/AO.30.001982 |
| 0.162 | olb/turbulence/plane_wave_scintillation.py:379 | strong-regime scattering-disk term | not in the book | - | - | NO; Churnside 1991, DOI 10.1364/AO.30.001982 |
| 1, 20e3, 20 | olb/turbulence/profiles.py:13 | altitude grid start, end and point count | - | - | - | not found (a numerical grid choice). The grid stops at 20 km, so it truncates the Hufnagel-Valley high-altitude layer near 10 km. |
| 0.6 | olb/turbulence/uplink_flux.py:48 | WEAK_FLUCTUATION_LIMIT on the log-amplitude variance | Ch. 8, text after Eq. (23); Ch. 12, Eqs. (40) and (93) | 264-265, 497, 522 | 289-290, 522, 547 | NO. Confirmed by 2 readers. The book limit sigma_R^2 < 1 means sigma_x^2 < 0.25, so the warning fires about 2.4 times too late. See Conflicts C-05. |
| 3.86 | olb/turbulence/uplink_flux.py:93 (comment) | longitudinal scintillation coefficient in the validation target | Ch. 8, Eq. (23); Summary Eq. (130); Ch. 9, Eqs. (92), (93), (148) | 264, 303, 351, 381 | 289, 328, 376, 406 | yes |
| 0.40 | olb/turbulence/uplink_flux.py:93 (comment) | coefficient on [(1+2Theta)^2+4Lambda^2]^(5/12) | Ch. 8, Eq. (23) | 264 | 289 | yes |
| 1 + 2 Theta, 2 Lambda | olb/turbulence/uplink_flux.py:93-94 (comment) | arguments of the power and the arctangent | Ch. 8, Eq. (23) | 264 | 289 | yes |
| 11/16 | olb/turbulence/uplink_flux.py:94 (comment) | coefficient on Lambda^(5/6) | Ch. 8, Eq. (23) | 264 | 289 | yes |
| 2 (in Lambda0 = 2L/(k W0^2)) | olb/turbulence/uplink_flux.py:98 | input-plane Fresnel parameter | Ch. 8, Eq. (5); Ch. 12, Eq. (8) | 261, 488 | 286, 513 | yes |
| 1 (in Theta0 = 1 - L/F0) | olb/turbulence/uplink_flux.py:103 | input-plane curvature parameter | Ch. 8, Eq. (5); Ch. 12, Eq. (8) | 261, 488 | 286, 513 | yes |
| Theta0/(Theta0^2 + Lambda0^2) | olb/turbulence/uplink_flux.py:104 | output-plane curvature parameter | Ch. 8, Eq. (6); Ch. 12, Eq. (9) | 261, 489 | 286, 514 | yes |
| 2.0 (in beta2 += 2 (sigma_theta L)^2) | olb/turbulence/uplink_flux.py:183 | mechanical jitter folded into the two-dimensional wander variance | Ch. 8, Eq. (32); Ch. 12, Eqs. (50) and (53) | 271, 502-503 | 296, 527-528 | DISPUTED, see Conflicts C-09. R4: yes as a per-axis to two-dimensional conversion. R7: not found, because Andrews keeps the wander variance and the pointing-error variance separate and NEVER sums a tracking jitter into the wander variance. |
| 0.5 (in sqrt(0.5 beta2)) | olb/turbulence/uplink_flux.py:190-191 | per-axis split of the two-dimensional wander variance | Ch. 8, Eq. (32); Ch. 12, Eq. (50) | 271, 502 | 296, 527 | yes; the wander displacement is the two-dimensional variance |
| 21.0 (w), 1.7e-14 (A) | olb/turbulence/uplink_flux.py:244 | Hufnagel-Valley 5/7 defaults | Ch. 12, text after Eq. (3) | 481 | 506 | yes |
| 0.42 | my_analysis_modules/coupled_flux.py:46 | spherical-wave Fried coefficient | Ch. 6, Eq. (116), times 2.1 | 209 | 234 | yes for the coefficient only; the path weight is MIRRORED (KR-01) |
| 4.2 | my_analysis_modules/coupled_flux.py:86 | short-term waist turbulence factor | none in Ch. 6 or Ch. 7 | - | - | not found (the Yura and Fried form used by Dios) |
| 0.26 | my_analysis_modules/coupled_flux.py:86 | short-term waist wander-removal factor | none in Ch. 6 or Ch. 7 | - | - | not found; the book Eq. (101) printed 206 uses 0.66 in a different parameterisation |
| 2.07 | my_analysis_modules/coupled_flux.py:113 | beam-wander variance coefficient | Ch. 6, Eqs. (93), (117), (118) give 7.25 | 203, 209, 210 | 228, 234, 235 | NO. Same integrand, coefficient 3.50 times low as a radial variance. R3 re-derived 8 pi^2 (0.033)(0.5) Gamma(1/6) = 7.252 from Eqs. (88) and (89). The book 2.07 (Ch. 9, eta_Y) is an unrelated quantity. See Conflicts C-01. |
| 2 (factor) | my_analysis_modules/coupled_flux.py:127 | W_LT^2 = W_ST^2 + 2 beta^2 | Ch. 6, Eq. (100) has the factor 1 | 205 | 230 | no; right only if beta^2 is a PER-AXIS variance, but `olb/turbulence/angle_of_arrival.py:57` calls the same number RADIAL. See Conflicts C-03. |
| 4.0 | my_analysis_modules/coupled_flux.py:167 | long-term waist turbulence term | none in Ch. 6 or Ch. 7 | - | - | not found (Dios Eq. (2), DOI 10.1364/AO.43.003866) |
| 1.33 | my_analysis_modules/coupled_flux.py:200 (docstring) | uniform-Cn2 limit of the spreading integral | Ch. 6, Eq. (86) | 202 | 227 | yes; 4.35 x 3/8 / 1.23 = 1.326 |
| 1.23 | my_analysis_modules/coupled_flux.py:200 (docstring) | Rytov variance | Ch. 6, Eq. (119) | 210 | 235 | yes |
| 5/3 | my_analysis_modules/coupled_flux.py:234 | the (1 - z/L) path-weight exponent | Ch. 6, Eq. (109); Ch. 8, Eq. (118) | 208, 299 | 233, 324 | yes |
| 0.033 | my_analysis_modules/coupled_flux.py:235, :291, :315 | Kolmogorov spectral constant | Ch. 6, Eq. (110); Ch. 8, Eq. (120) | 208, 299 | 233, 324 | yes |
| 4 pi^2 (0.033)(-0.5 Gamma(-5/6)) = 4.3508 | my_analysis_modules/coupled_flux.py:235 | slant-path spreading coefficient | Ch. 6, Eq. (109) | 208 | 233 | yes; the book prints 4.35 |
| 5/6 | my_analysis_modules/coupled_flux.py:236 | exponent on Lambda L/k | Ch. 6, Eq. (109) | 208 | 233 | yes |
| 5/6, 5/12 | my_analysis_modules/coupled_flux.py:288 | Kolmogorov exponents of the longitudinal integrand | Ch. 8, Eqs. (19) and (23) | 263-264 | 288-289 | yes for the exponents; the PARENTHESIS placement is wrong (KR-18) |
| -5/6 | my_analysis_modules/coupled_flux.py:312 | first parameter of 1F1 | Ch. 8, Eq. (18) | 263 | 288 | yes |
| 2 (in exp(2 xi) and exp(-2 beta^2/w^2)) | my_analysis_modules/coupled_flux.py:336 | irradiance from the log amplitude, and the Gaussian on-axis roll-off | Ch. 8, Eq. (13) | 262 | 287 | yes |
| 0.25 | my_analysis_modules/coupled_flux.py:382 | sigma2_x = 0.25 ln(1 + sigma2_I) | Ch. 8, Eq. (13) | 262 | 287 | yes |
| 1.23 | my_analysis_modules/general_atmospherics.py:22 | Rytov-variance coefficient | Ch. 5, Eq. (15); Ch. 6, Eq. (119); App. III Table VII(a) footer | 140, 210, 769 | 165, 235, 794 | yes |
| 7/6, 11/6 | my_analysis_modules/general_atmospherics.py:22 | Rytov-variance exponents | Ch. 5, Eq. (15) | 140 | 165 | yes |
| 0.54, 1.22, 0.509 | my_analysis_modules/general_atmospherics.py:23 | plane-wave scintillation index constants | Ch. 9, Eq. (47); Ch. 12, Eqs. (40) and (93); App. III Table VII(b) | 336, 497, 522, 769 | 361, 522, 547, 794 | NO. The book gives 0.49, 1.11 and 0.51. THIS IS THE ORIGIN of the fault; olb copied it into `plane_wave_scintillation.py:284`. One fix must cover both. |
| 0.69 | my_analysis_modules/general_atmospherics.py:23 | small-scale denominator constant | Ch. 9, Eq. (47); Ch. 12, Eq. (40) | 336, 497 | 361, 522 | yes |
| 0.423 | my_analysis_modules/general_atmospherics.py:24 | plane-wave Fried-parameter constant | Ch. 6, Eq. (64) times 2.1; Ch. 12, Eq. (23); Ch. 14, Eq. (25); App. III Table IV plus its note | 194, 492, 617, 767 | 219, 517, 642, 792 | yes; 2.1^(-5/3) x 1.46 = 0.4236, and the book prints 0.42 |
| Vg = 10 (default) | my_analysis_modules/general_atmospherics.py:26 | ground wind speed | Ch. 12, Eq. (3) | 481 | 506 | not found; Andrews leaves the ground wind a free site parameter and gives no default |
| slew rate = 1 deg/s (default) | my_analysis_modules/general_atmospherics.py:26 | satellite slew rate | Ch. 12, Eq. (3) | 481 | 506 | not found; Andrews leaves the slew rate free and gives no default |
| 30, 9400, 4800 | my_analysis_modules/general_atmospherics.py:41 | Bufton tropopause jet amplitude, centre altitude and width | Ch. 12, Eq. (3) | 481 | 506 | yes; all three match. R1 could not confirm this in Sec. 3.4, which has no wind profile. |
| 21 | my_analysis_modules/general_atmospherics.py:44 | default rms wind speed for the Hufnagel-Valley model | Ch. 12, text after Eq. (3) | 481 | 506 | yes; this is the book H-V5/7 |
| 0.00594, 27, 1e-5, power 10, 1000 | my_analysis_modules/general_atmospherics.py:56 | Hufnagel-Valley high-altitude wind-driven term | Ch. 12, Eq. (1) | 481 | 506 | yes |
| 2.7e-16, 1500 | my_analysis_modules/general_atmospherics.py:57 | Hufnagel-Valley middle term | Ch. 12, Eq. (1) | 481 | 506 | yes |
| 1.7e-14, 100 | my_analysis_modules/general_atmospherics.py:58 | Hufnagel-Valley ground term and its scale height | Ch. 12, Eq. (1) | 481 | 506 | yes |
| 1.33 | my_analysis_modules/general_atmospherics.py:108 | long-term spot-size growth constant | Ch. 6, Eq. (86) | 202 | 227 | yes; R1 could not confirm it inside Ch. 3 or App. III |
| 1 (Theta0 hardcoded) | my_analysis_modules/general_atmospherics.py:135 | collimated-beam curvature parameter | Ch. 6, Eq. (6); App. III Table VI footer | 183, 768 | 208, 793 | yes for the collimated case only. The kernel HARDCODES it; olb threads f0, so olb is the more general of the two. |
| 2 (in Lambda0 = 2L/(k W0^2)) | my_analysis_modules/general_atmospherics.py:152 | Gaussian beam parameter definition | Ch. 6, Eq. (6); App. III Table VI footer | 183, 768 | 208, 793 | yes |
| 0.125 | my_analysis_modules/general_atmospherics.py:152 | beam waist hardcoded inside `wLT` | - | - | - | not found (a hardcoded input, not a book constant; it silently overrides the caller) |
| 0.35 | my_analysis_modules/general_atmospherics.py:162 | undocumented length in the r0-based long-term spot size | - | - | - | not found |
| 2.1 | my_analysis_modules/general_atmospherics.py:441 | ratio of the Fried parameter to the coherence radius | Ch. 6, text below Eq. (64); App. III note under Table IV | 194, 767 | 219, 792 | yes |
| 6/5 | my_analysis_modules/general_atmospherics.py:508 | Fried parameter wavelength-scaling exponent | Ch. 6, Eq. (64); App. III Table IV | 194, 767 | 219, 792 | yes |
| -5/3, -3/5 | my_analysis_modules/general_atmospherics.py:526 | net Fried parameter of N layers | - | - | - | not found; it follows from the additive Cn2 path integral but the book does not state it |
| -11/12 | my_analysis_modules/general_atmospherics.py:542, :588 | Kolmogorov phase-screen Fourier amplitude exponent | Ch. 3, Eq. (18) | 67 | 92 | yes; the square root of a two-dimensional kappa^(-11/3) spectrum |
| 3.30 | my_analysis_modules/general_atmospherics.py:609 | Roddier random-tilt scale in a phase screen | Roddier, Progress in Optics 19 (1981) 281 | - | - | not found in Andrews |
| 1.0299 | my_analysis_modules/general_atmospherics.py:785 | total phase variance over a circular aperture | Ch. 14, Eq. (90); Noll 1976, DOI 10.1364/JOSA.66.000207 | 635 | 660 | yes; the book prints 1.03 (and 1.02 by the ABCD route) |
| 0.134 | my_analysis_modules/general_atmospherics.py:786 | tilt-removed phase variance over a circular aperture | Ch. 14, Eq. (94); Noll 1976, DOI 10.1364/JOSA.66.000207 | 636 | 661 | yes; the book prints 0.13 |
| 5.92 | ABSENT from olb | Tatarskii inner-scale wavenumber, km = 5.92/l0 | Ch. 3, Eq. (19), restated as Eq. (29); Ch. 8, Eq. (25) | 67, 265 | 92, 290 | not found in olb. The book fixes 5.92 so that the structure function takes the quadratic form of Eq. (13). See olb gap 6. |
| 3.3 | ABSENT from olb | modified atmospheric inner-scale wavenumber, kl = 3.3/l0 | Ch. 3, Eq. (22), restated as Eq. (31); Ch. 9, Eqs. (5) and (146) | 69, 328, 381 | 94, 353, 406 | not found in olb. See olb gap 6. |
| 1.802, 0.254 | ABSENT from olb | high-wavenumber bump terms of the modified atmospheric spectrum | Ch. 3, Eq. (22); Ch. 9, Eqs. (5) and (146) | 69, 328 | 94, 353 | not found in olb. See olb gap 6. |
| 2 pi, 4 pi, 8 pi (in k0 = C0/L0) | ABSENT from olb | outer-scale wavenumber convention | Ch. 3, Eqs. (20), (21), (23); Ch. 9, Eqs. (6) and (147) | 68-69, 328 | 93-94, 353 | not found in olb. The book uses 2 pi/L0 or 1/L0 for von Karman, and 4 pi/L0 or 8 pi/L0 for the modified spectrum. Any new spectra.py must make the convention explicit. |
| 2.914, 1.093 | ABSENT from olb | plane-wave and spherical-wave Kolmogorov structure-function constants | App. III Tables I and II; Ch. 9 Sec. 9.3.1 text prints 2.91 | 765, 330 | 790, 355 | not found in olb |
| 1.64, 1.87 | ABSENT from olb | plane-wave inner-scale coherence radius constants, von Karman and modified | App. III Table IV; Ch. 6, Eq. (83) | 767, 200 | 792, 225 | not found in olb. See olb gap 6. |
| 0.55, 0.62 | ABSENT from olb | spherical-wave coherence radius constants | App. III Table V | 767 | 792 | not found in olb |
| 0.618, 8/3 (in the a-factor) | ABSENT from the general path | Gaussian-beam structure function and coherence radius, a = (1 - Theta^(8/3))/(1 - Theta) | App. III Tables III and VI; Ch. 6, Eq. (55) | 766, 768, 192 | 791, 793, 217 | present in olb only through the collimated path. See olb gap 3. |
| 4.42, 3.86, 0.40, 11/16 | ABSENT from the code body | Andrews closed-form Gaussian-beam scintillation index | App. III Table IX(a); Ch. 8, Eq. (23) | 771, 264 | 796, 289 | present only in the `uplink_flux.py:93-94` COMMENT, never executed. See olb gap 9. |
| 7.25 | ABSENT from olb | beam-wander displacement variance coefficient | Ch. 6, Eqs. (93) and (117); App. III Table IX(b) footer | 203, 209, 772 | 228, 234, 797 | not found in olb; the kernel uses 2.07. See Conflicts C-01. |
| 2.42, 2.72 | ABSENT from olb | collimated and focused beam-wander closed forms | Ch. 6, Eqs. (95), (96), (127), (128) | 204, 212 | 229, 237 | not found in olb |
| 0.66 | ABSENT from olb | wander-removal factor of the book short-term radius | Ch. 6, Eq. (101) | 206 | 231 | not found in olb; the kernel uses the Dios 4.2/0.26 form instead |
| 2.91 | `structure.angle_of_arrival_variance` | one-axis aperture angle-of-arrival variance | Ch. 6, Eqs. (84) and (133); Ch. 12, Eqs. (28) and (90) | 201, 213, 492, 522 | 226, 238, 517, 547 | IMPLEMENTED; `angle_of_arrival.aperture_arrival_angle_variance` delegates to it (no Term reads it yet, 0-W3) |
| 0.182 | ABSENT from olb; named in the batch-2 flag as the target | one-axis tilt variance as 0.182 (D/r0)^(5/3)(lambda/D)^2 | none in Andrews | - | - | NO. The Andrews Eq. (84) route converts to 0.174: 2.91/(0.423 x 4 pi^2) = 0.1743. Andrews gives the GRADIENT tilt; 0.182 is the ZERNIKE tilt from Noll. R3 searched Ch. 6 and Ch. 7 (printed 179-255) and found no 0.182. See Conflicts C-04. |
| 0.81 (angle of arrival) | ABSENT from olb | outer-scale reduction of the angle of arrival | Ch. 6, Eqs. (83) and (133) | 200, 213 | 225, 238 | not found in olb |
| 1.062 | ABSENT from olb | the book's own weak aperture-averaging constant | Ch. 10, Eq. (61) | 412 | 437 | not found in olb; olb ships the Churnside 1.07 fit instead |
| 0.16 | ABSENT from olb as a literal | spherical-wave Fried radius constant | Ch. 10, text after Eq. (47); Ch. 9 Sec. 9.12; Ch. 8, Eq. (33) | 408, 384, 272 | 433, 409, 297 | not found as a literal; olb reaches it as 0.423 x 3/8 = 0.1586 |
| 0.18, 0.56 | ABSENT from olb | spherical-wave counterparts of the olb 0.65 and 1.11 | Ch. 10, Eq. (77); Ch. 9, Eqs. (96), (97), (102) | 416, 351-352 | 441, 376-377 | not found in olb; olb has no spherical-wave aperture averaging |
| 10.89, 64 pi^2 | ABSENT from olb | inner-scale and outer-scale parameters Q_l and Q_0 | Ch. 10, Eq. (68); Ch. 9, Eqs. (24) and (48) | 413, 332, 338 | 438, 357, 363 | not found in olb; olb has no inner or outer scale |
| 35.05 | ABSENT from olb | von Karman inner-scale parameter Qm | Ch. 8, Eq. (29) | 266 | 291 | not found in olb |
| 8 (in D_G^2 = 8 W_G^2) | ABSENT from olb | hard-to-soft aperture conversion | Ch. 10, text after Eq. (57); Ch. 12, Eq. (48) | 411, 500 | 436, 525 | not found in olb; olb uses a hard aperture only |
| 0.23 | ABSENT from olb | ln(10)/10, the dB-to-natural-log factor in the fade cumulative | Ch. 11, Eqs. (24), (26), (34); Ch. 12, Eqs. (69)-(71) | 451-455, 511 | 476-480, 536 | not found in olb; olb never uses the book closed-form fade cumulative |
| 2.61 | ABSENT from olb | eta_X numerator, plane wave | Ch. 9, Eqs. (40) and (155) | 335, 383 | 360, 408 | not found in olb |
| 0.20, 0.19, 0.23 (spherical) | ABSENT from olb | spherical-wave large-scale and small-scale constants | Ch. 9, Eqs. (69) and (72) | 342 | 367 | not found in olb |
| 0.86, 2.73, 1.87 | ABSENT from olb | saturation-regime asymptote constants | Ch. 9, Eqs. (18), (19), (23) | 331-332 | 356-357 | not found in olb |
| 2.64 (2.65 in the Summary) | ABSENT from olb as a literal | closed-form radial scintillation constant | Ch. 8, Eqs. (18) and (129) | 263, 302 | 288, 327 | not found; olb computes the same quantity numerically through the Dios integral |
| 0.48, 0.54 (pointing error) | ABSENT from olb | collimated and focused closed forms of the wander-induced pointing error | Ch. 8, Eqs. (36)-(38); Ch. 12, Eqs. (51) and (52) | 273-274, 502 | 298-299, 527 | not found in olb |
| 5.95 | ABSENT from olb | uplink untracked longitudinal scintillation index constant | Ch. 12, Eq. (54), repeated as Eq. (100) | 503, 524 | 528, 549 | not found in olb |
| 14.50, 8.70 | ABSENT from olb | slant-path untracked Gaussian-beam index coefficients | Ch. 8, Eq. (118) | 299 | 324 | not found in olb |
| 45e3 | ABSENT from olb | downlink correlation width, rc = sqrt(45e3 sec(zeta)/k) | Ch. 12, Eqs. (41)-(43) | 498 | 523 | not found in olb |
| 0.67 - 0.17 Theta | ABSENT from olb | the delta_t factor in the Gaussian-beam irradiance covariance | Ch. 8, Eq. (50); Summary Eq. (133) | 280, 304 | 305, 329 | not found in olb |
| 1.7, 3 | ABSENT from olb | irradiance correlation width factors, plane and spherical | Ch. 8, Eq. (53) | 281 | 306 | not found in olb |
| 3.503, 3.63 | ABSENT from olb | downlink and uplink second derivatives of the temporal covariance | Ch. 12, Eqs. (75)-(77) and (82)-(84) | 514, 519 | 539, 544 | not found in olb |
| 550 Hz | ABSENT from olb | the book's nominal quasi-frequency for its worked fade-time figures | Ch. 11, Sec. 11.3.3 text | 457 | 482 | not found in olb; olb has no temporal model |
| 3.44 | ABSENT from olb | long-term turbulence MTF constant | Ch. 14, Eqs. (26)-(28) | 617-618 | 642-643 | not found in olb |
| 0.28 | ABSENT from olb | short-term (tilt-removed) MTF constant | Ch. 14, Eqs. (33)-(35) | 621 | 646 | not found in olb |
| 0.32 | ABSENT from olb | Greenwood time constant for a constant wind, tau0 = 0.32 r0/V_perp | Ch. 14, Eq. (39) | 622 | 647 | not found in olb |
| 0.90 (tilt only) | ABSENT from olb | tilt-only aperture-averaged phase variance | Ch. 14, Eq. (93) | 635 | 660 | not found in olb; olb has piston-removed (1.0299) and piston-plus-tilt-removed (0.134) but not tilt alone |
| 0.6159, 0.0500 | ABSENT from olb | Sasiela Strehl-ratio asymptotic series | Ch. 14, Eq. (43) | 623 | 648 | not found in olb |
| 1.83 | ABSENT from olb | geometrical-optics phase variance, 1.83 (L0/r0)^(5/3) | Ch. 14, Eq. (88) | 635 | 660 | not found in olb |
| 2.27, 2.11 | ABSENT from olb | strong-regime benchmarks for r0 over the coherence radius | Ch. 7, text at Eq. (53) and Fig. 7.2 | 241, 244 | 266, 269 | not found in olb; they are the unit test that GF-09 needs |
| 4 and 5 (SCIDAR L0 caps) | ABSENT from olb | outer-scale altitude-profile caps | Ch. 12, Eqs. (6) and (7) | 483 | 508 | not found in olb |
| 1.22 (as q) | ABSENT from olb as a named quantity | strong-fluctuation parameter q = L/(k rho_pl^2) = 1.22 sigma_R^(12/5) | Ch. 6, Eq. (122); Ch. 7, Eq. (58); Ch. 9 printed 342 | 211, 242, 342 | 236, 267, 367 | present only folded into 0.81 and 1.63 at `gaussian_fried.py:97-98`. The SAME 1.22 appears WRONGLY at `plane_wave_scintillation.py:284`, one level too high. |

## Conflicts

Rows where readers, or the sources themselves, disagree. Each conflict needs an
owner decision. A conflict is NOT a plain bug: do not "fix" one of these until
the decision is made.

| id | subject | the two positions | evidence | recommendation |
|---|---|---|---|---|
| C-01 **RESOLVED 2026-08-25** | Beam-wander variance constant. `beam_wander_variance` at `coupled_flux.py:113` (Table 1 row KR-04), which feeds AA-01 and the whole uplink wander path. | (a) The kernel uses 2.07 and cites Dios et al. 2004, DOI 10.1364/AO.43.003866. (b) Andrews Ch. 6, Eqs. (93), (117) and (118) give 7.25 for an IDENTICAL integrand. | R3 re-derived 7.25 from Ch. 6, Eq. (88) with the filter of Eq. (89): 8 pi^2 (0.033)(1/2) Gamma(1/6) = 7.252. The book value enters W_LT^2 = W_ST^2 + <r_c^2> (Eq. (100)) and R3 verified it against the book Worked Example 2, printed 215. A full-book search for 2.07 returns only eta_Y = 2.07 sigma_R^(12/5) in Ch. 9, an unrelated quantity. The kernel is 3.50 times low as a radial variance, 1.75 times low as a per-axis variance. UPDATE (kernel-patch agent, 2026-08-25): the per-axis escape is RULED OUT. Andrews Eq. (93) reduces through Eq. (94) to the standard 2.42 Cn2 L^3 W0^(-1/3); the kernel 2.07 gives 0.69, a factor 3.50 under EITHER axis convention (a per-axis reading needs 3.63, not 2.07). The Dios paper is paywalled, so its own definition stays unverified. The kernel keeps 2.07 with a comparison comment. RESOLUTION (2026-08-25, the owner supplied the paper, `REFS/Dios et al. - 2004 ...pdf`): the kernel does NOT mis-copy Dios. Dios Eq. (11), printed p. 3868, prints exactly `<beta^2> = 2.07 INT_0^L Cn2(z)(L-z)^2 [1/W_s(z)]^(1/3) dz`. Dios does not derive it; he takes it from his reference 23, Belmonte, Applied Optics 39, 5426 (2000), DOI 10.1364/AO.39.005426, and he prints no filter function for it. Dios cites Andrews for the long-term spread (his Eqs. (4)-(6)) and for the scintillation index (his Eqs. (14)-(20)), but NOT for the beam wander, so the paper makes no attempt to reconcile 2.07 with 7.25. | **DECIDED: KEEP 2.07. The kernel is faithful to its citation and gets no change to the constant.** The gap is a true source-against-source difference, Belmonte against Andrews, and it is 3.50 under a consistent RADIAL convention on both sides (Dios Eq. (10) makes `<beta^2>` radial; Andrews `<r_c^2>` is radial too). Two reasons to keep 2.07. (1) Dios Fig. 3, printed p. 3871, plots his Eq. (11) against a split-step FFT-BPM wave-optics simulation of the same uplink, and the two agree; the text, printed p. 3870, calls the analytic and the numerical results "similar but slightly different". A factor of 3.50 would be plain on that figure. (2) The whole olb uplink chain (wander offset, off-axis index, Eq. (23) irradiance) is Dios, so mixing in an Andrews constant would break the internal consistency of the model. NEXT STEP to close the gap for good: read Belmonte 2000, DOI 10.1364/AO.39.005426, which is the only document that can say which filter gives 2.07. Until then `olb/turbulence/andrews/wander.py` stays the independent Andrews measurement, and the two numbers stay side by side. |
| C-02 | Spherical and Gaussian-beam path weight direction. `gaussian_fried.py:319` (GF-18) and `coupled_flux.py:46` (KR-01) use ((L-z)/L)^(5/3); Andrews Ch. 6, Eqs. (115) and (116) use (z/L)^(5/3). | (a) R3: the book measures z FROM THE TRANSMITTER and states below Eq. (116) that the Cn2 near the RECEIVER carries the weight, so olb uses the MIRROR weight. (b) The olb form matches Dios, DOI 10.1364/AO.43.003866, Eq. (3), which gives the TRANSMITTER-plane (reciprocal) coherence radius of an UPLINK beam. | THE GEOMETRY RESOLVES IT. R4 found Ch. 8, Eq. (115) printed 299: the slant UPLINK plane-wave weighting is (1 - z/L)^(5/6), heavy near the TRANSMITTER. R7 found the book's own uplink wave structure function and coherence radius, Ch. 12, Eqs. (24)-(27) printed 492, which is a separate uplink derivation from the downlink pair Eqs. (17)-(22). And R3 confirmed that the same kernel weights the SPREADING integral correctly at `coupled_flux.py:234` (KR-09), where the book Ch. 6, Eq. (109) also measures from the transmitter. So Ch. 6, Eqs. (115) and (116) give the DOWNLINK (receiver-referred) coherence radius, and the olb weight is the UPLINK (transmitter-referred) one. | NOT A BUG, a plane-of-reference difference. DO NOT flip the weight. Do two things instead. (1) Document the reference plane in the GF-18 and KR-01 docstrings, and say which link direction each serves. (2) Build Table 2 row G-130 (Ch. 12, Eqs. (24)-(27)) and check the olb uplink form against the book's own uplink form, not against the downlink Eq. (115). |
| C-03 **RESOLVED 2026-08-25** | Long-term waist factor. `coupled_flux.py:127` (KR-05) computes W_LT^2 = W_ST^2 + 2 <beta^2>; Andrews Ch. 6, Eq. (100) printed 205 has the factor 1 on the wander variance. | (a) The kernel factor 2 is correct if <beta^2> is a PER-AXIS variance, because the book <r_c^2> is the two-dimensional (radial) displacement variance. (b) `olb/turbulence/angle_of_arrival.py:57` documents the SAME number as the RADIAL (two-axis) variance, in which case the factor must be 1. | The two readings cannot both hold. R4 independently read the same quantity as two-dimensional at `uplink_flux.py:183` and :190-191 (Table 3 rows), where olb doubles a per-axis jitter and then halves it again per Cartesian axis. So olb itself uses BOTH conventions in two places. RESOLUTION (2026-08-25, from the paper): position (a) is WRONG and position (b) is RIGHT. Dios Eq. (9), printed p. 3868, gives `beta = sqrt(beta_x^2 + beta_y^2)`, and Dios Eq. (10) gives `<beta_x^2> = <beta_y^2> = 0.5<beta^2>`. So `<beta^2>` is the RADIAL (two-axis) variance. Dios Eq. (1), printed p. 3867, then prints `W_LT^2(z) = W_ST^2(z) + 2<beta^2>` on that radial quantity, and Dios repeats the same combination in his error term, Eq. (29), printed p. 3870. | **DECIDED: the convention is RADIAL, and the kernel factor 2 is Dios Eq. (1). No code change.** The factor 2 is not a per-axis to radial conversion; it is the paper's own factor. Every olb site is already consistent with the radial reading: `angle_of_arrival.py:57` calls it radial (correct), `uplink_flux.py:195-197` draws each Cartesian axis with the variance 0.5*beta2 (Dios Eq. (10), correct), and `uplink_flux.py:188` folds a per-axis jitter as 2*(sigma_theta L)^2 into the radial total (correct arithmetic, see C-09). The kernel docstrings of `beam_wander_variance` and `long_term_beam_waist` now state the radial convention and cite Dios Eqs. (1), (9), (10) and (11). RESIDUAL, for the record: Dios adds 2 x 2.07 and Andrews Eq. (100) adds 1 x 7.25, so the wander part of W_LT^2 differs by 1.75, not 3.50. The factor 2 and the constant 2.07 partly cancel. That 1.75 is the same number the WP5 table calls "the factor-2 escape route"; it is not an escape route, it is the true residual. |
| C-04 | Angle-of-arrival tilt coefficient. The batch-2 flag names 0.182 (D/r0)^(5/3)(lambda/D)^2; the Andrews route gives 0.174. | (a) Andrews Ch. 6, Eq. (84) printed 201 gives <beta_a^2> = 2.91 Cn2 L (2 W_G)^(-1/3), one axis. With D = 2 W_G and r0 = (0.423 k^2 Cn2 L)^(-3/5) this is 2.91/(0.423 x 4 pi^2) = 0.1743. (b) 0.182 is the Noll ZERNIKE tilt, a different tilt definition. | R3 searched Ch. 6 and Ch. 7 (printed 179-255) and found no 0.182 anywhere. Andrews Eq. (82) defines the tilt as a phase DIFFERENCE across the pupil (a gradient tilt); Noll defines it as the Zernike Z2/Z3 coefficient. | A DEFINITION CHOICE FOR THE OWNER, not an error in either source. Pick the tilt definition FIRST, then take the matching coefficient: gradient tilt gives 0.174 with Ch. 6, Eq. (84); Zernike tilt gives 0.182 with Noll 1976. Note that `olb/turbulence/ao.py` already uses the NOLL convention (1.0299 and 0.134), so the Noll 0.182 is the consistent choice inside olb. Record the choice in the `aperture_arrival_angle_variance` docstring before you fill the stub. |
| C-05 | The two weak-fluctuation limits: 0.25 on sigma2_I at `plane_wave_scintillation.py:45` (PW-01) and 0.6 on sigma2_x at `uplink_flux.py:48` (UF-01). | (a) R2 marked 0.25 `wrong` (a false citation); R5 and R7 marked it `unmatched` (not a book number); R6 marked it `approximate` (a defensible tighter house rule); R1 could not find it in range. (b) R4 and R7 both marked 0.6 `wrong`: it is 2.4 times looser than the book. | The book states ONE weak criterion in five places: sigma_R^2 < 1 (Ch. 5 printed 140; Ch. 8 printed 264; Ch. 9 printed 323; Ch. 10, Eq. (61) printed 412; Ch. 12, Eq. (40) printed 497). Ch. 5, Eq. (16) printed 140 adds that a Gaussian beam ALSO needs sigma_R^2 Lambda^(5/6) < 1. With sigma_I^2 = 4 sigma_x^2 (Ch. 8, Eq. (13)), sigma_R^2 = 1 maps to sigma_x^2 = 0.25. So the two olb limits are inconsistent with EACH OTHER: 0.25 on the index is 4 times tighter than the book, and 0.6 on the log amplitude is 2.4 times looser. Ch. 11.3 printed 451 does say the lognormal tail is optimistic against simulation, which argues for a tighter bound than the book's. | ONE RECOMMENDATION, three parts. (1) Write the BOOK criterion once, as a shared helper: sigma_R^2 < 1, plus sigma_R^2 Lambda^(5/6) < 1 for a beam wave (Table 2 row G-20). Cite Ch. 5, Eq. (16). (2) Keep a SEPARATE, tighter house limit for the LOGNORMAL Term only, and label it a house rule, not an Andrews number, with the Ch. 11.3 printed 451 tail argument as its justification. (3) Lower the uplink 0.6 to 0.25 in sigma2_x units, so it matches the book sigma_R^2 < 1. Do (3) first: it is the only one of the three that lets an untrustworthy number leave the budget with no warning. |
| C-06 | Aperture filter form. `plane_wave_scintillation.py:108` (PW-03) uses the hard Airy filter [2 J1(x)/x]^2. | (a) R6 marked it `approximate`: Andrews Ch. 10, Eq. (59) printed 412 uses a SOFT Gaussian aperture, exp(-D_G^2 kappa^2/16) with D_G^2 = 8 W_G^2; the Airy form is the Fourier transform of the hard circular MTF inside Eq. (54) printed 410; the two agree in the limits but not in between. (b) R7 marked it `exact`: the identical function IS printed in the book, as the piston Zernike filter, Ch. 14, Eq. (86) with m = n = 0, printed 634. | Both readings are correct about different pages. Andrews prints the function (Ch. 14) but does NOT use it for aperture averaging (Ch. 10 uses the soft aperture). The merged table keeps the stricter `approximate`. | No code change is implied. State in the docstring that the filter is the hard circular MTF, cite Ch. 14, Eq. (86) for the function, and note that the Andrews aperture-averaging chain (Ch. 10, Eqs. (57)-(69)) uses a soft Gaussian aperture instead. If a numerical comparison is wanted, build Table 2 row G-134 (Ch. 12, Eq. (39)), which IS a hard-aperture closed form. |
| C-07 | The (8/3)^(3/5) spherical-over-plane Fried ratio at `gaussian_fried.py:163` (GF-12). | (a) R3 marked it `approximate`: the book Ch. 6, Eq. (71) printed 196 gives rho_sp = (0.55 Cn2 k^2 L)^(-3/5), and the exact 3/8 weight gives 0.5475, so the book 0.55 is ROUNDED; olb 1.7963 against the book 1.7913 is a 0.3 % difference. (b) R5 marked it `exact`: the Ch. 9 worked example printed 384 gives r0 = (0.16 Cn2 k^2 L)^(-3/5), and 0.423 x 3/8 = 0.1586, which the book rounds to 0.16, so the book confirms the olb constant. | Both readers computed the same numbers. The disagreement is only about whether a rounded book constant makes the exact olb constant "approximate". | Keep the exact (8/3)^(3/5). It is the analytic value; the book prints a rounded one. Note the 0.3 % in the docstring so that a future numerical comparison against a book figure does not read as a fault. No code change. |
| C-08 | The "Ch. 9" docstring citation of the effective beam parameters at `gaussian_fried.py:96-99` (GF-06, GF-07). | (a) R3 called it a mis-citation, because the source is Ch. 7, Eq. (58) printed 242. (b) R5 found Ch. 9, Eqs. (85) and (86) printed 349, restated as Eq. (150) printed 382, which state the same result identically. | Ch. 9 Sec. 9.6.1 explicitly names Sec. 7.4.1 as the home of the derivation, then restates the result. | NOT a citation fault. Leave the docstring alone, or make it more specific by naming both: Ch. 7, Eq. (58) for the derivation and Ch. 9, Eqs. (85)-(86) for the restatement. Recorded here so that nobody "fixes" it. |
| C-09 **CONFIRMED 2026-08-25** | The jitter fold at `uplink_flux.py:183` (UF-08), beta2 += 2 (sigma_theta L)^2. | (a) R4 marked the factor 2 `yes`: Ch. 8, Eq. (32) printed 271 treats the wander displacement as the two-dimensional variance, so a per-axis jitter variance doubles correctly. (b) R7 marked it `not found`: Andrews keeps the wander displacement (Ch. 12, Eqs. (50) and (51) printed 502) and the pointing-error variance (Eq. (53) printed 503) as SEPARATE quantities that share one integral, and feeds them into the untracked index Eq. (54) and the tracked index Eq. (57). Andrews NEVER adds a mechanical tracking jitter into the wander variance. | The two readers agree on the ARITHMETIC (per axis to two dimensions) and disagree on the CONSTRUCTION (whether the book supports adding a mechanical jitter into the wander variance at all). DIOS CHECKED 2026-08-25: the paper does NOT treat mechanical jitter. It models atmospheric beam wander only; the words jitter, tracking and pointing error do not occur in it. So Dios gives the fold no support either. | Both are right, and the paper does not change that. The factor 2 is correct arithmetic, and it is now on a FIRM convention: Dios Eq. (10) makes beta2 the radial variance, so a per-axis jitter variance (sigma_theta L)^2 does add as 2*(sigma_theta L)^2. The construction is an olb extension with NO citation in either source. Keep the code, and label the extension in the docstring. The memory note "pointing jitter into beta" already records why olb does it this way (to avoid double-counting a stacked pointing Term). Build Table 2 row G-140 (Ch. 12, Eq. (53)) to give the book's own route, then compare the two numerically before you decide whether to change anything. |
| C-10 | Range-limited non-finds, as a class. 14 Table 1 rows were marked `unmatched` by one reader and then found by a second reader in another chapter. | The merged rows carry the POSITIVE result, and each note names the reader who could not find it and why. | Examples: AO-04 and PR-02 (Ch. 6 against Ch. 12), BW-17 and KR-26 and KR-27 (Ch. 3 and Ch. 8 against Ch. 12), DL-04 (Ch. 5 against Ch. 9), GF-08 and GF-09 (Ch. 9 against Ch. 7), KR-12 (Ch. 8 against Ch. 6), KR-28 (Ch. 3 against Ch. 6), KR-36 (App. III against Ch. 14). | No action. This is a property of the split reading, not a disagreement. Listed so that a future reader does not re-open them. |

### C-01/C-03 measurements (WP5)

`olb/turbulence/andrews/wander.py` is an INDEPENDENT implementation of the
Andrews beam-wander chain: Ch. 6, Eq. (93) (general), Eq. (94) (infinite outer
scale), Eq. (100) (short-term radius), Ch. 8, Eq. (36) (pointing error), and
Ch. 12, Eqs. (50) and (53) (slant path). It changes NOTHING in the kernel
`coupled_flux.py` and nothing in the olb Dios path. It only measures.

The module first reproduces the book's own numbers. Ch. 6 Worked Example 2
(printed p. 215) gives sqrt(<r_c^2>) = 3.35 cm; the module gives 3.3492 cm, a
0.02 % difference. The same example gives W_LT = 6.52 cm (module 6.5195 cm,
0.01 %) and W_ST = 5.59 cm (module 5.5934 cm, 0.06 %). Worked Example 4
(printed p. 216) gives 1.81 cm collimated and 1.90 cm for a beam focused at
900 m over a 1 km path; the module gives 1.811 cm and 1.902 cm.

Two cases. TERRESTRIAL: lambda 1550 nm, L 2000 m, Cn2 3e-16, W0 5 cm,
collimated. UPLINK: the defaults of the `olb/turbulence/uplink_flux.py`
self-check, so W0 1 m, lambda 1550 nm, slant range 600 km, zenith, the
`DEFAULT_HS` grid, and the HV5/7 profile `get_c2n(hs, 21.0, 1.7e-14)`.

| case | quantity | value |
|---|---|---|
| terrestrial | Andrews Eq. (93)/(94) `<r_c^2>` | 1.57436e-05 m^2 (rms 3.968 mm) |
| terrestrial | kernel 2.07, free-space W(z) | 4.48369e-06 m^2 |
| terrestrial | kernel 2.07, geometrical-optics W(z) = W0 | 4.49508e-06 m^2 |
| terrestrial | Eq. (94) reduced form 2.42 Cn2 L^3 W0^(-1/3) | 1.57436e-05 m^2 |
| terrestrial | **ratio Andrews / kernel (free-space W(z))** | **3.5113** |
| terrestrial | **ratio Andrews / kernel (same GOM W(z))** | **3.5024** |
| terrestrial | ratio Eq. (94) reduced / Andrews general | 1.0000 |
| terrestrial | ratio Eq. (94) reduced / kernel (GOM) | 3.5024 |
| uplink | Andrews Ch. 12, Eq. (50) slant `<r_c^2>` | 6.0264 m^2 (rms 2.45 m, 4.09 urad) |
| uplink | kernel 2.07, free-space W(z) | 1.72064 m^2 |
| uplink | **ratio Andrews / kernel** | **3.5024** |
| uplink | homogeneous 2.42 surrogate with Cn2 = mu0/L | 2.01517 m^2 |
| uplink | ratio surrogate / Andrews slant | 0.3344 |
| terrestrial | W_ST from the kernel, the shared input of the next three rows | 0.0537708 m |
| terrestrial | W_LT by Andrews Eq. (100), factor 1 on a RADIAL `<r_c^2>` | 0.053917 m |
| terrestrial | W_LT by the kernel, factor 2 on a PER-AXIS `<beta^2>` | 0.0538541 m |
| terrestrial | **ratio W_LT Andrews / kernel** | **1.00117** |
| both | 7.25 / (2 x 2.07), the factor-2 escape route | 1.7512 |

Reading of the table.

1. C-01. The constant ratio is 3.5024 = 7.25 / 2.07 EXACTLY, on both cases and
   on both path geometries, once the two sides read the same beam-radius
   profile W(z). The 3.5113 row differs from 3.5024 only because the kernel
   takes the true free-space (diffracting) W(z) while Andrews Eq. (93) takes
   the refractive (geometrical-optics) W(z) = W0 |Theta0 + Theta0_bar xi|. So
   the integrand is the same and the whole gap is the leading constant. The
   Eq. (94) reduced form and the Eq. (93) general form agree to 1 part in 1e12
   for a collimated beam, so the reduced form is not a separate position.
2. C-03. Feeding both combination rules the SAME short-term waist, the
   long-term waists differ by only 0.12 %. That is a COINCIDENCE of this case,
   not agreement: the kernel doubles a variance that is 3.50 times too small
   (2 x 2.07 = 4.14 against 7.25, a residual factor 1.7512), and W_LT is
   dominated by W_ST here, so the wander term is a small correction either way.
   The coincidence disappears for a small transmitter or a strong path, where
   the wander term dominates W_LT.
3. NO convention reading reconciles the two. A radial reading of the kernel
   leaves it 3.50 times low. A per-axis reading of the kernel leaves it 1.75
   times low, and a per-axis reading also contradicts
   `olb/turbulence/angle_of_arrival.py:57` and the Eq. (100) factor 1.
   SUPERSEDED 2026-08-25: the paper is now read (see the next block). Dios
   prints 2.07 and defines `<beta^2>` as RADIAL, so the kernel copies Dios
   correctly and the 3.50 gap is real. But the last row of the table above,
   1.7512, is NOT an escape route: it is the true residual difference in the
   long-term waist, because Dios adds 2<beta^2> where Andrews adds 1<r_c^2>.
4. Uplink note. The homogeneous 2.42 form with an equivalent constant
   Cn2 = mu0/L understates the slant answer by a factor 3 (row ratio 0.3344).
   That is not a conflict. On an uplink all the turbulence sits in the first
   20 km of a 600 km path, so the transmitter weight xi^2 is close to 1 over
   the whole turbulent layer instead of averaging to 1/3. Do not use the
   homogeneous reduction on a slant path.

Not built here, and named so that nobody looks for them: Ch. 8, Eqs. (40),
(41), (43) and (44), and Ch. 12, Eqs. (54), (56) and (57), put sqrt(<r_c^2>)
and sigma_pe into the tracked and untracked scintillation index. Those live in
`olb/turbulence/andrews/scintillation.py`, whose `scintillation_index` already
takes `wander_rms_m` and `pointing_error_m`. `wander.py` supplies exactly those
two quantities, as variances; take the square root at the call site.

### C-01, C-03 and C-09 adjudication against the Dios paper (2026-08-25)

The owner supplied the primary source, `REFS/Dios et al. - 2004 - Scintillation
and beam-wander analysis in an optical ground station-satellite uplink.pdf`
(Applied Optics 43 (19) 3866, DOI 10.1364/AO.43.003866). The kernel
`my_analysis_modules/coupled_flux.py` was then read equation by equation
against it. This block records the result.

**What the paper says**

| Dios equation | printed p. | statement | kernel function | verdict |
|---|---|---|---|---|
| Eq. (1) | 3867 | `W_LT^2(z) = W_ST^2(z) + 2<beta^2>` | `long_term_beam_waist` | exact copy |
| Eq. (2) | 3868 | collimated long-term waist | `long_term_beam_waist_collimated` | exact copy |
| Eq. (3) | 3868 | `r0,s` with the `((L-z)/L)^(5/3)` weight, stated "for the uplink" | `spherical_wave_coherence_diameter` | exact copy; also settles C-02 |
| Eqs. (4)-(6) | 3868 | second-order Rytov long-term waist | `long_term_beam_waist_rytov` | exact copy |
| Eq. (7) | 3868 | Yura short-term waist, 4.2 and 0.26 | `short_term_beam_waist` | exact copy |
| Eq. (9) | 3868 | `beta = sqrt(beta_x^2 + beta_y^2)` | - | `<beta^2>` is RADIAL |
| Eq. (10) | 3868 | `<beta_x^2> = <beta_y^2> = 0.5<beta^2>` | `uplink_flux.py:195-197` | exact copy |
| Eq. (11) | 3868 | `<beta^2> = 2.07 INT Cn2(z)(L-z)^2 [1/W_s(z)]^(1/3) dz` | `beam_wander_variance` | exact copy |
| Eqs. (15), (17), (18) | 3869 | `Lambda`, `Theta`, `A(z)`, `B(z)` | `_lambda_function`, `_theta_function`, `_A`, `_B` | exact copy |
| Eq. (16) | 3869 | on-axis index; the bracket closes AFTER the cosine | `on_axis_scintillation_index` | exact copy after the 2026-08 fix |
| Eq. (20) | 3869 | off-axis index with `1F1(-5/6, 1, 2r^2/W^2(L))` | `off_axis_scintillation_index` | exact copy |
| Eq. (23) | 3869 | `I = exp(2 chi) exp(-2 beta^2/W_ST^2(L))` | `on_axis_irradiance` | exact copy |
| Eq. (24) | 3869 | `<I(r,L)> = exp(-2 r^2/W_LT^2(L))` | `mean_off_axis_irradiance` | exact copy |
| Eq. (25) | 3870 | `sigma2_I,Gb = (sigma2_I + sigma2_I,r) <I>^2` | `coupled_flux_sample` | weight RESTORED 2026-08-25 |
| Eq. (26) | 3870 | `<chi> = -sigma_chi^2`, `sigma_chi^2 = 0.25 ln(1+sigma2_I,Gb)` | `coupled_flux_sample` | exact copy |

**C-01 — DECIDED: keep 2.07.** The kernel is a faithful copy of Dios Eq. (11),
in the constant, in the `(L-z)^2` weight and in the `[1/W_s(z)]^(1/3)` factor.
Dios does not derive Eq. (11); he takes it from his reference 23, Belmonte,
Applied Optics 39, 5426 (2000), DOI 10.1364/AO.39.005426, and he prints no
filter function for it. He cites Andrews for the long-term spread and for the
scintillation index, but NOT for the beam wander, so the paper makes no attempt
to reconcile 2.07 with the Andrews 7.25. The gap is therefore a genuine
source-against-source difference and it is 3.50 on BOTH sides of a radial
convention. Two reasons keep 2.07 in the kernel. (1) Dios Fig. 3, printed
p. 3871, plots Eq. (11) against a split-step (FFT-BPM) wave-optics simulation
of the same uplink, and the two agree closely; the text on p. 3870 calls the
analytic and the numerical results "similar but slightly different". A factor
of 3.50 would be plain on that log plot. (2) The olb uplink chain is Dios end
to end, so an Andrews constant inside it would break the model's internal
consistency. To close the gap for good, read Belmonte 2000 next: it is the only
document that can say which filter gives 2.07.

**C-03 — DECIDED: RADIAL, and the factor 2 is Dios's own.** Dios Eq. (10) makes
`<beta^2>` the radial (two-axis) variance, and Dios Eq. (1) then puts the factor
2 on that radial quantity. So position (a) of the C-03 row is wrong: the factor
2 is NOT a per-axis to radial conversion. Position (b) is right, and every olb
site already agrees with it. No code change. The residual difference from
Andrews Eq. (100) in the long-term waist is 1.75, not 3.50, because 2 x 2.07 =
4.14 against 7.25.

**C-09 — CONFIRMED: the jitter fold is an olb extension.** The paper models
atmospheric beam wander only. Mechanical jitter, tracking and pointing error do
not appear in it. So neither source supports adding a mechanical jitter into
the wander variance. The ARITHMETIC of `uplink_flux.py:188` is now on a firm
footing, because Dios Eq. (10) fixes the radial convention. Keep the code and
keep the extension labelled.

**Kernel edits made (2026-08-25)**

`my_analysis_modules/coupled_flux.py` is untracked in its own repository, so it
was edited in place and gets no commit.

1. `coupled_flux_sample` — the Dios Eq. (25) mean-irradiance weight is
   RESTORED. This reverses the 2026-08 KR-20 patch, which was made on the
   Andrews reading before the paper was available.
2. `beam_wander_variance` — the docstring now cites Dios Eq. (11) and Belmonte
   2000, states the RADIAL convention from Dios Eqs. (9) and (10), keeps the
   Andrews 3.50 comparison as a known difference, and gives the Fig. 3
   justification for keeping 2.07. The constant is unchanged.
3. `long_term_beam_waist` — the docstring now cites Dios Eqs. (1) and (29) and
   says the factor 2 is the paper's own factor on a radial variance. The code
   is unchanged.
4. `spherical_wave_coherence_diameter` — the docstring now cites Dios Eq. (3)
   and records the C-02 plane-of-reference decision. The code is unchanged.
5. `coupled_flux_sample` — the `wL` and `wL_lt` parameter descriptions are
   corrected. `wL` is the free-space width `W(L)` of Dios Eq. (15), not a
   short-term waist.

**Measured effect of restoring the Eq. (25) weight**

The weight `<I>^2 = exp(-4 beta^2/W_LT^2)` is 1 on axis and falls with the
wander offset, so it lowers the index and the fade. Measured on the
`olb/turbulence/uplink_flux.py` self-check case (W0 1 m, 1550 nm, 600 km,
zenith, HV5/7 `get_c2n(hs, 21.0, 1.7e-14)`), for which `<beta^2>` = 1.7206 m^2:

| wander offset `beta` | index without the weight | index with Eq. (25) | ratio |
|---|---|---|---|
| 0 (on axis) | 8.505 | 8.505 | 1.000 |
| 0.928 m (one axis, 1 sigma) | 40.96 | 34.50 | 0.842 |
| 1.312 m (radial rms) | 81.82 | 58.04 | 0.709 |
| 2.623 m (2 x radial rms) | 14641 | 3706 | 0.253 |

The self-checks stay green. `python my_analysis_modules/coupled_flux.py` passes
its cross-validation and its demo. `python -m olb.turbulence.uplink_flux`
passes, with these shifts (the runs are unseeded Monte Carlo, so a part of each
shift is sampling noise):

| self-check line | before | after |
|---|---|---|
| weak Cn2 `sigma2_x` | 0.0186 | 0.0177 |
| strong Cn2 `sigma2_x` | 6.2978 | 6.2418 |
| no jitter, loss | 12.073 dB | 12.038 dB |
| no jitter, 99 % fade | 46.974 dB | 42.855 dB |
| 5 urad jitter, loss | 18.131 dB | 18.055 dB |
| 5 urad jitter, 99 % fade | 438.345 dB | 422.473 dB |
| collimated `sigma2_x` | 0.9839 | 0.8959 |
| 5x diverged `sigma2_x` | 0.5139 | 0.4549 |

The deep-tail fade moves most, which is correct: the weight is smallest at a
large wander offset, which is where the deep fades come from.

### C-01 — CLOSED against Belmonte 2000 (2026-08-25)

The owner supplied the primary source that the Dios adjudication asked for:
Belmonte, "Feasibility study for the simulation of beam propagation:
consideration of coherent lidar performance", Applied Optics 39 (30) 5426
(2000), DOI 10.1364/AO.39.005426. It is Dios reference 23, so it is the origin
of the 2.07 constant in the olb uplink chain. The paper CLOSES C-01. Keep 2.07.

**The paper prints the exact kernel form.** Belmonte Eq. (21), printed p. 5435,
is `<rho^2(z)> = 2.07 INT_0^L Cn2(z)(L-z)^2 W_S(z)^(-1/3) dz`. It is identical to
Dios Eq. (11) and to `my_analysis_modules/coupled_flux.py:beam_wander_variance`,
in the constant 2.07, in the `(L-z)^2` moment arm and in the short-term-radius
weight `W_S(z)^(-1/3)`. The line below Eq. (21) makes it RADIAL: each of the two
orthogonal components has the standard deviation `<rho^2>^(1/2)/sqrt(2)`, so the
per-axis variance is `<rho^2>/2`. This is the same radial convention as the
Andrews `<r_c^2>` and as C-03.

**The two constants are two derivations of the SAME radial quantity, not two
quantities.** Both reduce, for a collimated homogeneous Kolmogorov path, to
`C Cn2 L^3 W0^(-1/3)` with the same `INT_0^1 xi^2 dxi = 1/3` (substitute
`xi = 1 - z/L`, so `(L-z)^2 = L^2 xi^2`). So the ratio is purely the leading
constant: Belmonte reduces to 0.69 (`= 2.07/3`), Andrews reduces to 2.42
(`= 7.25/3`), and 2.42/0.69 = 3.50 exactly. No radial-versus-per-axis factor and
no integrand weight closes the gap. `olb/turbulence/andrews/wander.py` asserts
this ratio (7.25/2.07) in its self-check.

**Provenance of 2.07: the image-motion level arm, not the beam-wave filter.**
Belmonte derives Eq. (21) by analogy with the image-motion (angle-of-arrival)
problem, "the level arm weighted proportionally to the strength of the turbulence
along the path" (printed p. 5435). His sources are Yura, JOSA 63, 567 (1973)
(ref 46, the short-term spread and the level arm), Hufnagel, The Infrared
Handbook (1978) (ref 47), and Mironov and Nosov, JOSA 64, 516 (1977) (ref 48).
So 2.07 is the Yura / Mironov-Nosov centroid-tilt constant. Andrews 7.25 comes
from the beam-wave large-scale spectral filter of Ch. 6, Eqs. (88), (89) and
(93). The Andrews filter keeps ALL large-scale refraction of the beam; the
image-motion route keeps the centroid tilt (G-tilt) only, which is a subset. So
the image-motion constant is the smaller one. Belmonte cites the Andrews group
(refs 49, 53, 56) for the spread and the coherence, but NOT for the wander, and
he prints no filter for Eq. (21). So neither source reconciles the constant on
paper.

**The tie-break is simulation, and it backs 2.07.** Belmonte measures the true
centroid displacement directly from a split-step (phase-screen) wave-optics
simulation, through his Eq. (20), the intensity first moment. His Figs. 11 and
12, printed p. 5435, compare Eq. (21) with that simulation over 2- and 10-um
beams: the level-arm form matches, and it runs slightly LOW at the longest
ranges, which the paper attributes to scintillation left out of the level-arm
analysis (printed p. 5435). A constant 3.50 times larger would sit far above the
simulated points on those plots. His Section 4 conclusion, printed p. 5442,
states that beam wander "is properly considered by the available analytical
approach when turbulence is weak to moderate"; only strong-turbulence beam
break-up saturates the deflection, which is the regime olb sends to a numerical
path anyway. This is a SECOND independent split-step confirmation of 2.07, after
Dios Fig. 3. So the empirical evidence is unambiguous even though the analytic
derivations disagree.

**Consequence for olb.** The kernel and the whole Dios uplink chain keep 2.07;
no code changes. `olb/turbulence/andrews/wander.py` is the book-faithful Andrews
7.25 form and it stays as built, because it is the Andrews layer, but NO budget
uses it for the uplink wander (the coupled-flux path is Dios end to end). Any
future caller that wants a simulation-validated beam wander must use the Dios
kernel route (2.07), not the Andrews `beam_wander_variance` (7.25), and the
wander.py docstring already carries this warning. C-01 is resolved.

MACHINE-READABLE RECORD: `C01_WANDER`, a `kind="conflict"` `Constraint` on
`coupled_flux.beam_wander_variance` in `olb/turbulence/coupled_flux.py`. It cites
both sides (Belmonte 2.07 against Andrews 7.25). A traced Term that reads the Dios
wander kernel inherits it, and `constraints_frame()` lists it.

---

## WP3 notes — spectra, structure functions and aperture averaging

Work package WP3 built `olb/turbulence/andrews/spectra.py`, `structure.py` and
`aperture.py`, and it filled the inner-scale and outer-scale branches of
`andrews/scintillation.py`. This block records what was built, what was
delegated, and what the book would not give.

### Built

- `spectra.py` — Ch. 3, Eqs. (18) to (23), printed pp. 67 to 69. The five
  models `kolmogorov`, `tatarskii`, `von_karman`, `exponential` and
  `modified_atmospheric`, plus the plain dict `SPECTRA`. On the outer-scale
  constant the book is explicit and inconsistent on purpose: Eq. (20), printed
  p. 68, and Eq. (22), printed p. 69, print k0 = 2 pi/L0 "or sometimes
  k0 = 1/L0"; Eq. (23), printed p. 69, prints k0 = 4 pi/L0 "or 2 pi/L0, or
  8 pi/L0"; the Ch. 9 scintillation model uses C0 = 8 pi. Each function names
  its default and takes the constant as a keyword. Closes G-04 to G-09.
- `structure.py` — the wave structure function and the coherence radius of
  Appendix III, Tables I to VI, printed pp. 765 to 768, with Ch. 6, Secs. 6.4
  and 6.5. Also `fried_parameter` (r_0 = 2.1 rho_0, Ch. 6, text below Eq. (64),
  printed p. 194), `angle_of_arrival_variance` (Ch. 6, Eqs. (83) and (84),
  printed pp. 200 and 201) and `rms_image_jitter` (Ch. 6, Eq. (85), printed
  p. 201). Closes G-30 to G-36 and the coherence-radius half of G-43.
- `aperture.py` — Ch. 10, Sec. 10.3. Plane weak Eq. (60) and the Eq. (61) fit,
  spherical weak Eq. (53), plane strong Eq. (69) and the two-scale Eqs. (62) to
  (68), spherical strong Eq. (77) and the two-scale Eqs. (71) to (76), Gaussian
  strong Eqs. (87) to (90). Closes G-107, G-109, G-110, G-111, G-112 and G-114.
- `scintillation.py` — `two_scale_parameters`, `weak_two_scale_index` (Ch. 9,
  Eqs. (48), (75) and (104)) and the two-scale large-scale and small-scale log
  variances for the plane and the spherical wave. Closes G-54, G-87, G-88 and
  G-91.

### Decisions recorded

- **C-04 is DECIDED by the owner: the GRADIENT tilt.**
  `structure.angle_of_arrival_variance` returns
  <beta_a^2> = 2.91 Cn2 L D^(-1/3) = 0.174 (D/r_0)^(5/3)(lambda/D)^2 per axis,
  Ch. 6, Eq. (84), printed p. 201. The Noll Zernike tilt 0.182 is NOT what it
  returns, and both docstrings say so. The deferred stub
  `angle_of_arrival.aperture_arrival_angle_variance` now calls it, with the
  same signature. Note that `olb/turbulence/ao.py` still uses the NOLL
  coefficients 1.0299 and 0.134, so the package now holds BOTH tilt
  conventions. A caller that adds the two must say which one it means.
  MACHINE-READABLE RECORD: `_C04_TILT_CONFLICT`, a `kind="conflict"` `Constraint`
  on `ao.apply_compensation` in `olb/turbulence/ao.py`. A traced Term inherits it
  and `constraints_frame()` lists it.
- **C-06 is honoured.** `aperture.py` uses the book's own SOFT Gaussian
  aperture, D_G^2 = 8 W_G^2 (Ch. 10, text below Eq. (57), printed p. 411). It
  does NOT reuse the olb Airy filter. The module docstring states the
  hard-against-soft difference and points to Ch. 14, Eq. (86), printed p. 634.
  MACHINE-READABLE RECORD: `RECEIVER_AIRY_CONFLICT`, a `kind="conflict"`
  `Constraint` in `olb/turbulence/andrews/aperture.py`.
- **G-108 is confirmed.** The `aperture.py` docstring records that the book has
  NO annular receive aperture, so olb gap 8 needs another source.
- **C-07 is honoured.** `gaussian_fried.spherical_wave_fried_parameter` keeps
  the exact (8/3)^(3/5). The docstring now says the book row 0.55 gives 1.7913
  in place of 1.7963, a 0.3 % difference. The self-check MEASURES that 0.273 %.
  MACHINE-READABLE RECORD: `C07_SPHERICAL_RATIO`, a `kind="conflict"`
  `Constraint` in `olb/turbulence/gaussian_fried.py`.
- **C-02 is honoured.** `structure.py` takes one path length and one scalar
  Cn2, so it makes no path integral and picks no reference plane. Its docstring
  says so and points at C-02.
  MACHINE-READABLE RECORD: `PATH_WEIGHT_CONFLICT`, a `kind="conflict"`
  `Constraint` in `olb/turbulence/andrews/structure.py`.

### Reading resolved by a limit, not by the scan

- Ch. 9, Eq. (108), printed p. 355, and Ch. 10, Eq. (80), printed p. 419, carry
  the factor (1/3 - Theta/2 + Theta^2/5). The scan cannot tell Theta from
  Theta-bar. The factor is the integral INT_0^1 xi^2 (1 - Theta_bar xi)^2 dxi of
  Eq. (107), printed p. 354, so it is **Theta-bar**. The plane-wave limit then
  gives 0.49/3 = 0.163, which is the 0.16 of Eq. (55); the spherical-wave limit
  gives 0.49/30 = 0.0163, which is the 0.04 beta_0^2 of Eq. (72). Both check.
- Ch. 10, Eq. (53), printed p. 409, carries two signs the scan cannot resolve:
  the sign inside the hypergeometric argument and the sign of the 11/16 term.
  Only one of the four combinations gives both A(0) = 1 and a factor that falls
  to zero. `aperture._weak_spherical_factor` uses that one and says so. The
  coefficients 9.66 and 11/16 are as printed.

### The book would not give (owner action needed)

1. **Appendix III, Table III, printed p. 766 — the Gaussian row of the MODIFIED
   spectrum.** The two Lambda-only bump terms read as
   0.438 (Lambda Q_l)^(1/6) and 0.056 (Lambda Q_l)^(1/6). Those fall only as
   Lambda^(1/6), which breaks the plane-wave reduction by 2.3 %. Ch. 6, text
   below Eq. (77), printed p. 197, states the Gaussian row MUST reduce to the
   plane row. So the numerators are wrong as read, and they are NOT guessed.
   `structure.wave_structure_function` raises NotImplementedError for
   wave="gaussian" with spectrum="modified". Every other cell of Tables I to III
   is built and its plane and spherical reductions are measured.
2. **Ch. 9, Eq. (109), printed p. 355, and Ch. 10, Eq. (84), printed p. 420 —
   the Gaussian-beam eta_X of the two-scale STRONG theory.** No reading
   recovered from the scan gives both the plane-wave value 2.61 (Ch. 9,
   Eq. (54), printed p. 339) and the spherical-wave value 8.56 (Ch. 10,
   Eq. (74), printed p. 415) in the two limits. So the Gaussian two-scale strong
   branch is NOT built, in `scintillation.large_scale_log_variance` and in
   `aperture.averaged_index`. Both raise NotImplementedError with the citation.
   The WEAK Gaussian two-scale index, Ch. 9, Eq. (104), printed p. 354, IS
   built, and its three limits are measured inside 1 %.
3. **Ch. 10, Eq. (78), printed p. 419 — the weak Gaussian-beam flux-variance
   double integral.** The book prints no closed form. `averaged_index` raises
   for wave="gaussian" with regime="weak" and points at the all-regime chain.

### Measured findings

- **The two-scale large-scale log variance does NOT reduce to the Kolmogorov
  one as l0 goes to zero.** Ch. 9, Eq. (54), printed p. 339, states its
  substitution L/(k rho_0^2) = 1.02 sigma_R^2 Q_l^(1/6) for the case
  rho_0 << l0 ONLY. So the two branches agree only where
  0.45 sigma_R^2 Q_l^(1/6) equals 1.11 sigma_R^(12/5). Measured at
  sigma_R^2 = 7.1, lambda = 1550 nm, L = 2 km, the ratio two-scale over
  Kolmogorov runs 0.85 (l0 = 1 mm), 1.23 (3 mm), 1.43 (5 mm), 1.67 (10 mm) and
  0.004 (l0 = 1 nm). The book's claim below Eq. (61), printed p. 413, that the
  chain "reduces" as Q_l goes to infinity is loose. Treat the two-scale branch
  as a moderate-to-strong model WITH a real inner scale, not as a superset. The
  outer-scale term DOES vanish exactly as L0 grows.
- **The Andrews weak plane aperture-averaging models against the Churnside fit
  that olb ships** (lambda = 1550 nm, L = 2 km, sigma_R^2 = 0.021):

  | d | Ch. 10, Eq. (60) exact | Ch. 10, Eq. (61) fit | Churnside 1.07 | exact/Churnside |
  |---|---|---|---|---|
  | 0.5 | 0.72192 | 0.75979 | 0.82487 | 0.875 |
  | 1.0 | 0.43413 | 0.42986 | 0.48309 | 0.899 |
  | 2.0 | 0.14785 | 0.14455 | 0.15643 | 0.945 |
  | 5.0 | 0.01927 | 0.02089 | 0.02139 | 0.901 |

  The Churnside fit is optimistic by 5 % to 13 % against the book's own exact
  Eq. (60) over that range. The Eq. (61) fit tracks Eq. (60) inside 8 %, which
  agrees with the book's own "less than 7 % error" claim at printed p. 412.
- **The Gaussian Eq. (88) chain is an independent fit.** It does not reduce
  exactly to the plane Eq. (69) or to the spherical Eq. (77). Measured ratios at
  sigma_R^2 = 0.71: plane limit 1.04 (d = 0.5), 1.10 (d = 1), 1.17 (d = 2);
  spherical limit 0.99, 0.97, 0.91. That is the size of the book's own fitting
  error, and it confirms the reading of Eq. (88).
- **`ao.plane_wave_fried_parameter_profile` changed by -0.135 %.** It used the
  Fried 1966 constant 0.423. It now delegates to the Andrews chain
  r_0 = 2.1 (1.46 Cn2 k^2 L)^(-3/5), which is the equivalent of 0.4240. The book
  itself prints the rounded 0.42 at Ch. 12, Eq. (23), printed p. 492. Measured
  at 60 deg elevation, 1550 nm, HV57: 17.2915 cm against 17.3149 cm.

### Delegations (names, signatures and docstrings kept)

| old home | new home |
|---|---|
| `plane_wave_scintillation.coherence_radius` | `andrews.structure.coherence_radius` |
| `plane_wave_scintillation.aperture_averaged_index_andrews` | `andrews.aperture.averaged_index` |
| `gaussian_fried.plane_wave_coherence_radius` | `andrews.structure.coherence_radius` |
| `gaussian_fried.plane_wave_fried_parameter` | `andrews.structure.fried_parameter` |
| `gaussian_fried.spherical_wave_fried_parameter` | through `plane_wave_fried_parameter`, exact (8/3)^(3/5) kept |
| `ao.plane_wave_fried_parameter_profile` | `andrews.structure.fried_parameter` |
| `angle_of_arrival.aperture_arrival_angle_variance` | `andrews.structure.angle_of_arrival_variance` |

The Churnside trio `aperture_averaging_factor_weak`, `_weak_inner` and
`_strong` KEEPS its own bodies, because those constants are Churnside 1991,
DOI 10.1364/AO.30.001982, not Andrews. Each docstring now names its book-form
alternative in `andrews.aperture`.

---

## WP6 notes — slant paths and the satellite link

Work package WP6 built `olb/turbulence/andrews/paths.py`, the Chapter 12
slant-path module. This block records what was built, what it measures, and the
two readings that the book leaves open.

### Built

- `sec_zeta` — Ch. 12, Eq. (14) geometry, printed p. 490. The module constant
  `ZENITH_LIMIT_DEG = 60.0` carries the book's own bound on the weak-fluctuation
  slant results: Ch. 12, Sec. 12.1, printed p. 478 ("less than 60 deg in most
  cases but may be restricted to zenith angles less than 45 deg in cases where
  ground-level Cn2 is large"), repeated at Ch. 12, Sec. 12.9, printed p. 521.
  The book gives NO Earth-curvature correction: it writes H = h0 + L cos(zeta)
  and puts sec(zeta) in front of each path integral. So the geometry is
  plane-parallel by construction, and the 45 to 60 deg bound is the only limit
  the book states.
- `hufnagel_valley` and `bufton_wind` — Ch. 12, Eqs. (1) and (3), printed
  p. 481. Both are THIN CITED WRAPPERS on the shared kernel `get_c2n` and
  `v_wind`, which reader R7 verified exact. No second implementation exists.
- `rms_wind` — Ch. 12, Eq. (2), printed p. 481. Closes G-124.
- `outer_scale_profile` plus the plain dict `OUTER_SCALE_MODELS` — Ch. 12,
  Eqs. (6) and (7), printed p. 483. Closes G-126 and G-127.
- `mu` — the path moments mu_0 to mu_3, both directions: Ch. 12, Eqs. (18),
  (19), (21), (25), (26), (37) and (55), printed pp. 491 to 503, restated as
  Eqs. (85) to (88) and (94) to (96), printed pp. 522 and 523.
- `downlink_scintillation_index` — Ch. 12, Eq. (38) (point, weak), Eq. (39)
  (hard aperture D_G, weak) and Eq. (40) (point, weak to strong), printed
  pp. 495 to 497. Closes G-134.
- `uplink_scintillation_index` — Ch. 12, Eqs. (54), (56), (57) and (58) (weak)
  and Eqs. (59) to (61) (weak to strong), printed pp. 503 to 506. It composes
  `wander.beam_wander_variance_slant` and `wander.pointing_error_variance_slant`
  instead of repeating them. Closes G-141, G-142, G-143 and G-144.
- `uplink_coherence_radius` — Ch. 12, Eqs. (24) to (27), printed p. 492. This is
  row G-130.
- `isoplanatic_angle` — Ch. 12, Eq. (29) (Gaussian beam) and Eq. (30)
  (spherical wave), printed p. 493. Closes G-131.
- `point_ahead_angle` — Ch. 12, Sec. 12.3.3, printed p. 488. Closes G-128.

### G-130 RESOLVED, and C-02 with it

The book's own uplink coherence radius, Ch. 12, Eq. (27), and the kernel
`spherical_wave_coherence_diameter` (`coupled_flux.py:21`) are DIFFERENT
QUANTITIES, not two readings of one quantity. Measured at 1.06 um, 60 deg
elevation, GEO (H = 38.5e3 km), W0 = 2 cm collimated, H-V5/7:

| quantity | value |
| --- | --- |
| book Ch. 12, Eq. (27), rho_0 of the uplink wave AT THE SATELLITE | 907.489 m |
| kernel `spherical_wave_coherence_diameter` | 0.11263 m |
| ratio kernel / Eq. (27) | 1.241e-04 |
| book Ch. 12, Eq. (23), GROUND Fried parameter r0 | 0.11262 m |
| **ratio kernel / Eq. (23)** | **1.000025** |
| ratio kernel / (2.1 x downlink Eq. (22) rho_0) | 1.001528 |

Reading of the table.

1. The kernel is Andrews Ch. 12, Eq. (23) to 2.5 parts in 1e5. The kernel
   weight ((L-z)/L)^(5/3) is 1 to within (h/L)^(5/3) over the whole turbulent
   layer of a satellite path, so the weight drops out and the kernel reduces to
   the flat mu_0 integral of Eq. (23). Eq. (23) is exactly the r0 that the book
   ITSELF feeds into the uplink beam-wander and pointing-error equations (50),
   (51) and (53), printed pp. 502 and 503. So the kernel computes the right
   number for the job it does.
2. Eq. (27) is the coherence radius IN THE SATELLITE PLANE. It weights the
   turbulence by the distance from the ground transmitter, which makes it
   hundreds of metres. The book states that below Eq. (27), printed p. 492:
   "the spatial coherence radius at the satellite will be many times larger than
   the probable size of the satellite". A 907 m result is that statement.
3. So C-02 is CLOSED. GF-18 and KR-01 are NOT mirrored. The kernel weight is
   the transmitter-referred one, and the book's own uplink chain uses the same
   reference plane. DO NOT flip the weight. The only action left on C-02 is the
   docstring note, which `paths.py:uplink_coherence_radius` now carries.

### A reading the book leaves open: mu_1

Ch. 12, Eq. (18) (downlink) and Eq. (25) (uplink) both PRINT the mu_1 bracket
with the plain height fraction (h - h0)/(H - h0). Read that way the two
equations are identical and the downlink coherence radius at the ground comes
out near 900 m, which is absurd. `paths.py` uses
|Theta + Theta_bar (1 - xi)|^(5/3) instead, that is the weight of the distance
FROM THE TRANSMITTER, with xi from Eq. (14). Three book facts fix that reading:

1. Ch. 6, Eq. (115), printed p. 209, gives the same moment on a general slant
   path as INT Cn2(z) |Theta + Theta_bar z/L|^(5/3) dz with z from the
   transmitter, and Ch. 6, Eq. (116) confirms it at Theta = Lambda = 0.
2. The text below Ch. 12, Eq. (19), printed p. 491, states mu_1d = mu_0 for a
   downlink from space. Only the (1 - xi) reading gives that.
3. The text below Ch. 12, Eq. (27), printed p. 492, states the uplink coherence
   radius at the satellite is huge. Only the (1 - xi) reading gives that.

Measured, same case as the table above: the (1 - xi) reading gives
mu_1d = 2.2339e-12, which equals mu_0 = 2.2340e-12 as fact 2 needs; the literal
reading gives 1.9277e-19. The book's own Worked Example 2, printed p. 525,
prints mu_1d = 1.98e-19, which is the LITERAL reading and which contradicts
facts 2 and 3. Recorded so that nobody "fixes" the module towards the worked
example. The module self-check prints both numbers.

### Measured findings

1. **The Andrews slant uplink longitudinal index and the Dios path integral
   agree to 0.02 %.** Case: 1.06 um, 600 km slant range, 60 deg elevation,
   W0 = 10 cm collimated, H-V5/7 on a 4001-point 0 to 20 km grid. Andrews
   Ch. 12, Eq. (58) gives 1.117061e-02; `beam_wave_scintillation.
   on_axis_scintillation_index` (Dios, DOI 10.1364/AO.43.003866) gives
   1.117327e-02, a difference of -0.02 %. This is the gap-9 twin on a real
   slant path, and it is tighter than the +3.06 % of the horizontal case,
   because a satellite path is close to the spherical-wave limit where the two
   forms share their leading constant. The UNTRACKED index of Ch. 12, Eq. (54)
   is 2.980385e-01, that is 27 times the Rytov value. The whole difference is
   the beam-wander pointing term, which is exactly the effect that Ch. 12.6.3
   exists to model (see Figs. 12.13 and 12.14, printed pp. 505 and 507). So the
   Dios route in olb reports the TRACKED index, and any untracked uplink budget
   that uses it alone understates the on-axis scintillation by more than an
   order of magnitude for a 10 cm beam.
2. **Ch. 12, Eq. (39) against the olb hard-Airy numerical integral.** This is
   Conflict C-06 measured. Same GEO case, downlink, 1.06 um, 60 deg elevation.
   No assert: the two use different aperture filters.

   | D_G [m] | Eq. (39) | olb Airy integral | ratio |
   | --- | --- | --- | --- |
   | 0.05 | 7.364280e-02 | 7.245216e-02 | 1.0164 |
   | 0.20 | 1.606343e-02 | 1.405438e-02 | 1.1429 |
   | 1.00 | 4.344872e-04 | 4.359956e-04 | 0.9965 |

   The two agree at the small-aperture and the large-aperture end and part by
   14 % in the middle, which is where the filter shape matters. Eq. (39) is a
   CLOSED FORM and needs no wavenumber grid, so it should supersede the
   numerical integral in `plane_wave_scintillation` when a later work package
   moves the downlink Term over.
3. **Eq. (39) reduces to Eq. (38) at D_G = 0 to +0.077 %.** The residual is the
   book rounding: 8.70 cos(5 pi/12) = 2.2517 and the book prints 2.25.
4. **Ch. 12, Eq. (29) reduces to Eq. (30) exactly** (measured 2.2e-16) in the
   spherical-wave limit, and the Andrews route agrees with the Stone route
   already in `anisoplanatism.isoplanatic_angle` to 9.0e-04. The whole
   difference is the constant, 2.91 against Stone's 2.914381:
   (2.914381/2.91)^(-3/5) = 0.99910.
5. **The book's own numbers reproduce.** Ch. 12, Worked Examples 1 and 2,
   printed pp. 524 and 525: mu_0 = 2.2340e-12 (book 2.24e-12), mu_3u =
   3.6382e-17 (book 3.70e-17), W = 750.0 m (book 750 m), r0 = 11.26 cm (book
   11.24 cm), uplink tracked 0.0688 (book 0.07), uplink untracked 0.0928 (book
   0.095), downlink on axis 0.1260 (book 0.13), isoplanatic angle 13.73 urad
   (book 13.5). The untracked value reproduces with the wander module default
   C_r = 2 pi, which is the value that makes the book's own worked example come
   out; the book leaves C_r free (Ch. 12, text below Eq. (53), printed p. 503)
   and uses 3.86/r0 in Fig. 12.13 and pi/r0 in Figs. 12.14 to 12.17.

### Bearing on olb gap 2 (the NO-SCINTILLATION corrected uplink)

`olb/links/uplink.py` flags its beacon-plus-adaptive-optics budget NO
SCINTILLATION, because the phase-only error budget drops the intensity
fluctuation that the coupled-flux Term carried. Ch. 12, Eqs. (57) to (60) give
the TRACKED index: a tracked uplink beam still scintillates by sigma_B_u^2 on
axis. Tracking removes the WANDER term, not the RYTOV term. On the case above
that index is 1.1e-02, not zero. WP6 read `uplink_scintillation_index(...,
tracked=True)` as the floor the corrected budget must carry. The owner
REJECTED that reading (2026-08-27, backlog 0-W1): the tracked form removes
the wander fully, which is a perfect tilt correction, the beacon decorrelates
from the uplink path over the point-ahead angle mode by mode, and a
decorrelated correction reshapes the beam, so the form is OPTIMISTIC for a
pre-compensated uplink and is not a bound in either direction. No analytic
Term will model the pre-compensated scintillation. The model of record is the
fidelity-1 FAST Monte Carlo with the point-ahead offset, wired 2026-08-27 as
`uplink_fast_term` (`uplink_budget(fidelity=1)`, the default for a
pre-compensated scenario).

### Refused, and named so that nobody looks for it

- An inner scale or an outer scale on ANY Ch. 12 slant scintillation form.
  Chapter 12 uses the Kolmogorov spectrum only (Ch. 12, Eq. (15), printed
  p. 490). Both `downlink_scintillation_index` and `uplink_scintillation_index`
  raise NotImplementedError on `l0` or `L0` and name
  `andrews.scintillation.weak_two_scale_index` as the single-path route. No
  coefficient is guessed.
- An aperture-averaged downlink index in the STRONG regime. Eq. (39) is a weak
  form and Eq. (40) is a point form; the book gives no product of the two.
- The SLC day and night profiles, Ch. 12, Eqs. (4) and (5), printed p. 482
  (row G-125). Not built: no caller needs them, and `hufnagel_valley` covers
  every olb site today.
- The full Gaussian-beam DOWNLINK index, Ch. 12, Eqs. (36) and (37), printed
  p. 495 (row G-133). Not built: the book itself reduces the downlink to the
  plane wave (Ch. 12, text below Eq. (21), printed p. 491). Compose
  `mu(..., order=3, direction="downlink")` with a beam if it is ever needed.

---

# WP4 — temporal statistics and the two fade-rate faces

This block records the work package that built
`olb/turbulence/andrews/temporal.py` and filled the two fade-rate stubs of
`olb/turbulence/andrews/distributions.py`. It closes the Table 2 rows G-11,
G-12, G-67 to G-71, G-96, G-99, G-115, G-116, G-119 to G-121, G-149, G-159 and
G-160. Row G-41 (Ch. 6.7, the mutual-coherence temporal spectrum) and rows G-75,
G-97, G-98, G-151 stay open.

## What is now built

| Function | Book source | printed p |
|---|---|---|
| `taylor_wavenumber` | Ch. 3, Eq. (27) | 73 |
| `fresnel_frequency` | Ch. 8, text below Eq. (57) | 283 |
| `irradiance_temporal_spectrum`, weak | Ch. 8, Eq. (65), with Eqs. (57) and (59) as its limits | 285, 283, 284 |
| `irradiance_temporal_spectrum`, strong | Ch. 10, Eqs. (93)-(97); at D = 0, Ch. 9, Eqs. (126)-(128) | 421-422, 365 |
| `quasi_frequency` | Ch. 11, Eqs. (14), (15), (35), (38); Ch. 12, Eq. (73) | 448, 456, 514 |
| `greenwood_frequency`, `coherence_time` | Ch. 14, Eqs. (38) and (39); text | 622, 623 |
| `expected_number_of_fades` | Ch. 11, Eqs. (34) and (37); Ch. 12, Eqs. (72) and (74) | 455-456, 513-514 |
| `mean_fade_time` | Ch. 11, Eq. (39); Ch. 12, Eqs. (78) and (79) | 456, 515 |

## Findings

1. **Ch. 8, Eq. (65) prints one factor in the wrong place.** The printed
   equation puts `(omega/omega_t)^(-8/3)` in front of BOTH spectral groups. Its
   own plane-wave and spherical-wave limits, Ch. 8, Eqs. (57) and (59), printed
   pp. 283 and 284, put it in front of the FIRST group only. The code follows
   Eqs. (57) and (59). With the printed placement the second group vanishes at
   high frequency, and the constants 0.72 and 0.24 of Eqs. (57) and (59) have
   nothing to cancel. With the code reading, Eq. (65) gives Eq. (57) at
   Theta = 1, Lambda = 0 (measured amplitude 6.949 against the printed 6.95) and
   Eq. (59) at Theta = 0, Lambda = 0 (measured 5.445 against the printed 5.47).

2. **Ch. 8, Eq. (64) leaves a branch cut open.** Eq. (65) writes the second
   group as `0.29 i^(4/3) a_j^(-4/3)`, and `a_j` is complex. The printed form
   does not name the branch, and the principal branch does NOT reproduce
   Eq. (57). The code writes every argument through `q_j = 1/(4 i a_j)`, which
   is real and positive for a plane wave (1/2) and a spherical wave (2/9),
   exactly as Eqs. (57) and (59) print the arguments.

3. **The printed constants 0.72 and 0.24 make the spectrum go NEGATIVE.** Each
   of the two spectral groups carries a 1/omega tail. The tails must cancel,
   because the book states at printed p. 283 that the spectrum decays as
   omega^(-8/3). The exact coefficient is
   `C q^(4/3)` with `C = -Gamma(-1/3) Gamma(11/6)/[Gamma(1/2) Gamma(7/3)] =
   1.810729`, which gives 0.7186 and 0.2437. The book rounds them to two
   figures. With the printed 0.72 the residual tail takes the plane-wave
   spectrum below zero above about 100 Fresnel frequencies (measured: -3.06e-10
   at 10 kHz on the module test case), and the spectral integral no longer
   converges. The code uses the derived C, and the self-check prints both.

4. **The quasi-frequency nu0 has NO upper limit of its own.** The book defines
   nu0 by the second spectral moment (Ch. 12, Eq. (73), printed p. 514). With a
   Kolmogorov spectrum and a zero inner scale the spectrum decays as
   omega^(-8/3), so the moment integrand decays as f^(-2/3) and b_2 grows as
   f_max^(1/3). The measured growth is x1.49 per decade of band. This is why the
   book sets nu0 to a fixed 550 Hz for its figures (printed pp. 457 and 514)
   instead of computing it. Any olb caller MUST set the band from the detector
   bandwidth or from an inner scale. The book gives no temporal spectrum with an
   inner scale, so olb cannot close that gap from this source.

5. **The book gives NO numeric worked example for the fade rate or the fade
   time.** Ch. 11.7 Example 1, printed pp. 472-473, stops at the probability of
   fade. Ch. 11 Problem 6, printed p. 474, asks for both at nu0 = 100 Hz and
   prints no answer. Ch. 12.10 has none either. So the two new faces are checked
   against the book's own internal identities instead: the rate equals nu0 at
   the threshold where 0.23 F_T equals sigma_l2/2 (printed p. 448), the
   gamma-gamma rate matches the printed Eq. (37) to 2e-15, and
   Pr(fade) = <n> <t> holds to machine precision for all three models.

6. **Ch. 9, Eq. (126) prints 0.50 where Ch. 9, Eq. (46) and Ch. 10, Eq. (95)
   print 0.51** for the small-scale log-irradiance amplitude. The code uses
   0.51, so that the strong temporal covariance at zero lag matches the
   aperture-averaged index of Ch. 10, Eq. (69), printed p. 413. The residual
   difference measured at zero lag is +0.58 %, which comes from the small-scale
   limit `(x)^(5/12) K_(5/6)(sqrt(x)) -> 0.5 Gamma(5/6) 2^(5/6) = 1.0056`, not 1.

7. **The weak and the strong temporal spectra carry the same POWER but not the
   same SHAPE.** At sigma_R^2 = 7e-4 the two integrals agree to 0.30 %, but the
   ratio of the two spectra is 1.374 at 1 Hz, 0.848 at 100 Hz and 0.422 at
   1 kHz on the module test case. Ch. 8.5 comes from the Rytov covariance;
   Ch. 9.8 comes from the two-scale extended-Rytov covariance. A caller that
   needs a spectral SHAPE must pick one and say which.

## Refused, and named so that nobody looks for it

- A temporal spectrum with a finite inner scale or outer scale, in ANY regime.
  Ch. 9.8, printed p. 364, states "We will also ignore the effects of a finite
  inner scale and outer scale". Ch. 10, printed p. 425, states only that a scale
  changes the peak values and does not shift the peak position. No closed form
  is printed, so `irradiance_temporal_spectrum` raises NotImplementedError on
  `l0` or `L0`.
- A strong-regime spherical wave or Gaussian beam. Ch. 9.8, printed p. 364,
  limits the analysis to a plane wave, and Ch. 9, printed p. 364, states that no
  Gaussian-beam covariance has been computed.
- A weak-only aperture-averaged temporal spectrum. Ch. 10.3.6, printed
  pp. 421-422, gives the all-regime form only, so `D` needs `regime="strong"`.
- The off-axis (radial) weak temporal spectrum, Ch. 8, Eqs. (66) and (67),
  printed pp. 286-287 (row G-71). Not built: no olb caller reads an off-axis
  spectrum today. The longitudinal part is built.

---

# WP7 — the wire-in

This block records the last work package of the Andrews foundation layer. WP7
wrote no new physics. It connected the nine modules to the link budgets and
synchronised the documentation.

## What was wired

1. **`olb/turbulence/andrews/__init__.py` exports all nine modules.** WP1 exported
   `distributions.py` only. The file is now a flat name list with no logic, so
   `from olb.turbulence.andrews import <name>` reaches every public function,
   and `python -m olb.turbulence.andrews.<module>` still runs each self-check.
2. **The gamma-gamma downlink Term is real.** `_gamma_gamma_term` in
   `olb/links/downlink.py` no longer raises. It composes the slant plane-wave
   Rytov variance (Ch. 12, Eq. (38), printed p. 495), the two log variances
   (Ch. 9, Eqs. (41) and (46), printed pp. 335 and 336), the gamma-gamma
   parameters (Ch. 9, Eq. (138), printed p. 370) and the Term adapter
   `olb/models/fade.py`. Table 1 row DL-05 is now `incorporated`.
3. **`_auto_select` switches model.** Below sigma2_I = 0.25 it returns the
   lognormal Term; at or above it, the gamma-gamma Term. The old warning that
   said only the lognormal model exists is gone. The switch point is the house
   limit of Conflict C-05, not the book limit sigma_R^2 = 1: the book bound is 4
   times looser, and the gamma-gamma chain of Ch. 12, Eq. (40), printed p. 497,
   is valid at every strength, so the early switch costs no validity. The
   gamma-gamma Term takes a scalar elevation only, because the quantile and the
   sampler carry one (alpha, beta) pair; an elevation array is refused, and the
   selector keeps the lognormal Term and warns.
4. **`olb/assumptions.py` names five spectra.** `SPECTRUM_TATARSKII`,
   `SPECTRUM_EXPONENTIAL` and `SPECTRUM_MODIFIED` join the two that were there,
   one for each model in `andrews/spectra.py` (Ch. 3, Eqs. (18) to (23)). No
   Term uses the new three yet.

## Measured

The downlink self-check now exercises the strong path. A 0.7 m ground aperture
at 1550 nm, a 600 km orbit, 15 deg elevation, the default H-V Cn2 profile:

| Quantity | Value |
|---|---|
| point sigma_R^2 (Ch. 12, Eq. (38)) | 0.7567 |
| alpha, beta (Ch. 9, Eq. (138)) | 4.850, 3.143 |
| sigma_I^2 (Ch. 9, Eq. (139)) | 0.5899 |
| the same by Ch. 12, Eq. (40) through `andrews.paths` | 0.5892 (+0.114 %) |
| gamma-gamma mean loss | 1.1901 dB |
| gamma-gamma 99 % fade | 10.0710 dB |
| lognormal point 99 % fade, same case | 8.8074 dB |

Two notes on that table. First, the +0.114 % is a DATUM difference, not a
physics difference: `andrews/paths.py` integrates (h - h0)^(5/6) and
`plane_wave_scintillation.py` integrates h^(5/6), and `DEFAULT_HS` starts at
h0 = 1 m. Second, the gamma-gamma fade is 1.26 dB DEEPER than the lognormal fade
of the same point index. That is the tail that Ch. 11, Sec. 11.3, printed
p. 451, says the lognormal model misses.

## Still open after WP7

- **The pre-compensated uplink carries NO SCINTILLATION, by DECISION.** WP6
  recorded `andrews/paths.uplink_scintillation_index(tracked=True)` as the
  floor of the residual. The owner rejected that reading on 2026-08-27: the
  form is optimistic for a pre-compensated beam, and no trustworthy analytic
  form exists for that case. The budget stays phase-only and mean-only, with
  loud flags. The model of record is the fidelity-1 FAST route, wired
  2026-08-27 as `uplink_fast_term`. See backlog 0-W1 and the "Bearing on olb
  gap 2" note above.
- **The gamma-gamma downlink Term models a POINT receiver.** The book gives no
  aperture-averaged downlink index in the moderate-to-strong regime. The Term
  flags that through its `Assumptions` record. Its fade is deeper than the true
  aperture fade, which is the safe direction.
- **`downlink_budget` still asks for `model="lognormal"`.** WP7 did not change
  the budget default, because the gamma-gamma Term drops the aperture averaging
  and so it would change the budget total by several dB at a low elevation. The
  owner must choose that. Call `downlink_scintillation_term(model="auto")`
  directly for the strong case.
- **The annular receive aperture (olb gap 8) needs a source that is not this
  book.** See the WP3 notes.
- **Conflict C-01 is CLOSED (2026-08-25) with Belmonte 2000.** The paper prints
  the same 2.07 wander form, gives its image-motion provenance, and validates it
  against a split-step simulation. Keep 2.07 in the Dios kernel; the Andrews
  7.25 form over-counts by 3.50. See the C-01 closure block.
- **The curvature-general Fried parameter is WIRED (2026-08-27).** The
  terrestrial fibre-coupling call site now passes f0 from the transmitter
  divergence (through `olb.beam.launch_curvature`), so a diverged beam drives
  its own r0. The single-path `gaussian_fried.gaussian_fried_parameter` still
  keeps its collimated signature (a tidy-up; the budgets use the profile form,
  which is general in f0). See `docs/physics.md` Section 5e.
- **TL-05 code half FIXED (2026-08-29)**: the terrestrial weak gate is now the
  shared beam-aware `rytov_weak(sigma2_R, Lambda)` of
  `olb/turbulence/andrews/scintillation.py`, which applies BOTH Ch. 5, Eq. (16)
  conditions (binding strength `sigma2_R * max(1, Lambda**(5/6))`), so a focused
  beam is caught. The same helper serves the uplink coupled-flux path (Dios edge
  `UPLINK_SIGMA2X_LIMIT = 0.6`). The lognormal-PDF house rule is now the distinct
  `LOGNORMAL_PDF_LIMIT = 0.25` on sigma2_I. This closes the code half of
  Conflict C-05 recommendation parts (1) and (3), and the beam half of TL-05.
  Remaining: the downlink (plane-wave) and `fast.py` limits are untouched
  follow-ups (both individually correct, no longer conflated).
