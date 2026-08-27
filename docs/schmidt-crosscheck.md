# Schmidt cross-check tracker

This file tracks the sections of Schmidt that we must check the olb wave-optics
layer against, or bring into it. It is the Schmidt twin of
`docs/andrews-crosscheck.md`. Andrews and Phillips gives the ANALYTIC physics
(fidelity 0 and 1). Schmidt gives the NUMERICAL method (fidelity 2): the
propagators, the sampling constraints, the absorbing boundary, and the phase
screens. So this tracker maps `olb/waveoptics/` and
`olb/waveoptics/turbulence/`, not the analytic Terms.

## Source

Jason D. Schmidt, *Numerical Simulation of Optical Wave Propagation with
Examples in MATLAB*, SPIE Press Monograph PM199, 2010. DOI 10.1117/3.866274.

The citation format for this book, everywhere in this repository, is:

    Schmidt (2010), DOI 10.1117/3.866274, Ch. N, Eq. (nn), printed p. NNN

The equation numbers are per chapter, and the book prints the chapter number in
the equation number. So "Ch. 7, Eq. (7.14)" is equation 14 of Chapter 7. Write
the full printed form `(7.14)`, not `(14)`. The book gives NO end-of-chapter
equation collation. Use the section index below as the index.

### Page rule

    printed p. N = PDF p. N + 13

The offset is +13. Four checks give the same value:

| printed p | content | PDF p |
| --- | --- | --- |
| 1 | Chapter 1 first page | 14 |
| 87 | Chapter 6 first page | 100 |
| 152 | Section 9.2.1 | 165 |
| 187 | Appendix B first page | 200 |

The PDF holds 212 pages. The front matter (roman pages i to xi) fills PDF pages
1 to 13.

### Access

The Read tool FALSELY refuses this PDF as password-protected. Use pdftotext
from the Bash tool, and stream to stdout:

    pdftotext -f <first> -l <last> "C:\Users\alexf\Zotero\storage\74VRZ7IT\Schmidt - 2010 - Numerical Simulation of Optical Wave Propagation with Examples in MATLAB.pdf" -

Give PDF page numbers to `-f` and `-l`, not printed page numbers.

## Status keys

- **flagged** — the owner marked the section. We have not checked it yet.
- **checked** — we compared the section to the code. See the note for the
  result.
- **incorporated** — the section is now in a module. See the note for the
  module and the `docs/physics.md` section.

## Policy: the book is the numerical authority

Where Schmidt and Andrews both give a quantity, Andrews owns the ANALYTIC
value and Schmidt owns the SIMULATION rule. Example: the Fried parameter comes
from Andrews Ch. 12; the per-screen `r0i` split and the screen count come from
Schmidt Ch. 9. A conflict between the two is a real finding. Record it in
Table 2, do not average the two.

---

# Chapter and section index

`printed` is the arabic page in the book. `pdf` is the page in the Zotero PDF.

## Chapter 1 — Foundations of Scalar Diffraction Theory (printed 1–13, pdf 14–26)

| Section | Title | printed | pdf |
| --- | --- | --- | --- |
| 1.1 | Basics of Classical Electrodynamics | 1 | 14 |
| 1.1.1 | Sources of electric and magnetic fields | 2 | 15 |
| 1.1.2 | Electric and magnetic fields | 2 | 15 |
| 1.2 | Simple Traveling-Wave Solutions to Maxwell's Equations | 5 | 18 |
| 1.2.1 | Obtaining a wave equation | 5 | 18 |
| 1.2.2 | Simple traveling-wave fields | 7 | 20 |
| 1.3 | Scalar Diffraction Theory | 9 | 22 |
| 1.4 | Problems | 12 | 25 |

## Chapter 2 — Digital Fourier Transforms (printed 15–38, pdf 28–51)

| Section | Title | printed | pdf |
| --- | --- | --- | --- |
| 2.1 | Basics of Digital Fourier Transforms | 15 | 28 |
| 2.1.1 | Fourier transforms: from analytic to numerical | 15 | 28 |
| 2.1.2 | Inverse Fourier transforms: from analytic to numerical | 17 | 30 |
| 2.1.3 | Discrete Fourier transforms in software | 18 | 31 |
| 2.2 | Sampling Pure-Frequency Functions | 21 | 34 |
| 2.3 | Discrete vs Continuous Fourier Transforms | 23 | 36 |
| 2.4 | Alleviating Effects of Discretization | 26 | 39 |
| 2.5 | Three Case Studies in Transforming Signals | 30 | 43 |
| 2.5.1 | Sinc signals | 30 | 43 |
| 2.5.2 | Gaussian signals | 31 | 44 |
| 2.5.3 | Gaussian signals with quadratic phase | 33 | 46 |
| 2.6 | Two-Dimensional Discrete Fourier Transforms | 35 | 48 |
| 2.7 | Problems | 37 | 50 |

## Chapter 3 — Simple Computations Using Fourier Transforms (printed 39–53, pdf 52–66)

| Section | Title | printed | pdf |
| --- | --- | --- | --- |
| 3.1 | Convolution | 39 | 52 |
| 3.2 | Correlation | 43 | 56 |
| 3.3 | Structure Functions | 47 | 60 |
| 3.4 | Derivatives | 50 | 63 |
| 3.5 | Problems | 53 | 66 |

## Chapter 4 — Fraunhofer Diffraction and Lenses (printed 55–64, pdf 68–77)

| Section | Title | printed | pdf |
| --- | --- | --- | --- |
| 4.1 | Fraunhofer Diffraction | 55 | 68 |
| 4.2 | Fourier-Transforming Properties of Lenses | 58 | 71 |
| 4.2.1 | Object against the lens | 59 | 72 |
| 4.2.2 | Object before the lens | 59 | 72 |
| 4.2.3 | Object behind the lens | 61 | 74 |
| 4.3 | Problems | 64 | 77 |

## Chapter 5 — Imaging Systems and Aberrations (printed 65–85, pdf 78–98)

| Section | Title | printed | pdf |
| --- | --- | --- | --- |
| 5.1 | Aberrations | 65 | 78 |
| 5.1.1 | Seidel aberrations | 66 | 79 |
| 5.1.2 | Zernike circle polynomials | 66 | 79 |
| 5.1.2.1 | Decomposition and mode removal | 73 | 86 |
| 5.1.2.2 | RMS wavefront aberration | 75 | 88 |
| 5.2 | Impulse Response and Transfer Function of Imaging Systems | 77 | 90 |
| 5.2.1 | Coherent imaging | 77 | 90 |
| 5.2.2 | Incoherent imaging | 79 | 92 |
| 5.2.3 | Strehl ratio | 82 | 95 |
| 5.3 | Problems | 84 | 97 |

## Chapter 6 — Fresnel Diffraction in Vacuum (printed 87–113, pdf 100–126)

| Section | Title | printed | pdf |
| --- | --- | --- | --- |
| 6.1 | Different Forms of the Fresnel Diffraction Integral | 88 | 101 |
| 6.2 | Operator Notation | 89 | 102 |
| 6.3 | Fresnel-Integral Computation | 90 | 103 |
| 6.3.1 | One-step propagation | 90 | 103 |
| 6.3.2 | Two-step propagation | 92 | 105 |
| 6.4 | Angular-Spectrum Propagation | 95 | 108 |
| 6.5 | Simple Optical Systems | 102 | 115 |
| 6.6 | Point Sources | 107 | 120 |
| 6.7 | Problems | 113 | 126 |

### Key equations of Chapter 6

| Eq. | printed | pdf | Gloss |
| --- | --- | --- | --- |
| (6.1) | 88 | 101 | The Fresnel diffraction integral. The base form. |
| (6.5) | 88 | 101 | The FT form. One Fourier transform gives the observation field. |
| (6.6) | 88 | 101 | The convolution form. The kernel is the free-space amplitude spread function. |
| (6.7)–(6.11) | 89 | 102 | The operator definitions: Q quadratic phase, V scale, F and F^-1 transforms, R the Fresnel integral. |
| (6.12) | 89 | 102 | Q2, the quadratic-phase operator in the FREQUENCY domain. |
| (6.13) | 90 | 103 | Q2 written as Q. A convenience, used in Sec. 6.4. |
| (6.15) | 90 | 103 | One-step propagation as an operator chain: Q V F Q. |
| (6.16) | 90 | 103 | The FIXED observation pitch of one step: Delta2 = lambda z / (N Delta1). ASSUMPTION: one step gives NO control of the output pitch. |
| (6.18) | 93 | 106 | Two-step propagation as an operator chain. It frees the output pitch. |
| (6.19)–(6.21) | 94 | 107 | The intermediate-plane and observation-plane pitches of the two steps. |
| (6.24) | 94 | 107 | The scaling parameter m as the ratio of the two step distances. |
| (6.25) | 94 | 107 | The intermediate-plane location. Two solutions give the same m. |
| (6.30) | 95 | 108 | The proof that both solutions give abs(m). |
| (6.31) | 95 | 108 | The angular-spectrum form. Two transforms and a transfer function. |
| (6.32) | 95 | 108 | The free-space transfer function H(f) = exp(ikz) exp(-i pi lambda z f^2). ASSUMPTION: paraxial, parallel planes. |
| (6.38) | 98 | 111 | The rearrangement of the quadratic exponent that introduces m. |
| (6.41) | 98 | 111 | The scaled source field U', with the (1-m) quadratic phase folded in. |
| (6.47) | 99 | 112 | The impulse response of the scaled convolution. |
| (6.49) | 99 | 112 | The amplitude transfer function of the scaled convolution. |
| (6.50) | 99 | 112 | The SCALED angular-spectrum operator chain. This is the workhorse of Chs. 7–9. |
| (6.51)–(6.54) | 99–100 | 112–113 | The pitch chain: Delta_f1 = 1/(N Delta1), and Delta2 = m Delta1. |
| (6.55), (6.56) | 100 | 113 | Two pitch identities. Chapter 8 uses them. |
| (6.67) | 101 | 114 | The compact plus-or-minus m angular-spectrum form. Both signs are valid. |
| (6.70)–(6.72) | 103 | 116 | The ray matrices: the general ABCD, ray transfer, and refraction. |
| (6.76) | 104 | 117 | The thin-lens ray matrix. |
| (6.77) | 104 | 117 | The generalized Huygens-Fresnel (ABCD) integral. ASSUMPTION: azimuthal symmetry only. |
| (6.80) | 104 | 117 | The ABCD integral as a convolution. |
| (6.81) | 104 | 117 | The ABCD transfer function. |
| (6.82) | 107 | 120 | The true point source as a Dirac delta. It has infinite bandwidth. |
| (6.85) | 108 | 121 | The observation field of a point source. A paraxial spherical wave. |
| (6.86) | 108 | 121 | The WINDOWED target field. The window must be larger than the detector and smaller than the grid. |
| (6.89) | 110 | 123 | The model point source that gives that windowed target field. |
| (6.92) | 110 | 123 | The sinc model point source, for a square window of width D. |

## Chapter 7 — Sampling Requirements for Fresnel Diffraction (printed 115–131, pdf 128–144)

| Section | Title | printed | pdf |
| --- | --- | --- | --- |
| 7.1 | Imposing a Band Limit | 115 | 128 |
| 7.2 | Propagation Geometry | 117 | 130 |
| 7.3 | Validity of Propagation Methods | 120 | 133 |
| 7.3.1 | Fresnel-integral propagation | 120 | 133 |
| 7.3.1.1 | One step, fixed observation-plane grid spacing | 120 | 133 |
| 7.3.1.2 | Avoiding aliasing | 121 | 134 |
| 7.3.2 | Angular-spectrum propagation | 124 | 137 |
| 7.3.3 | General guidelines | 128 | 141 |
| 7.4 | Problems | 130 | 143 |

### Key equations of Chapter 7

| Eq. | printed | pdf | Gloss |
| --- | --- | --- | --- |
| (7.1) | 116 | 129 | The Nyquist criterion: Delta <= 1/(2 f_max). |
| (7.3)–(7.5) | 116 | 129 | The plane wave and its direction cosines. The spatial-frequency spectrum IS the plane-wave spectrum. |
| (7.6) | 116 | 129 | Delta1 <= lambda / (2 theta_max). It ties the pitch to the ray angle. |
| (7.7) | 117 | 130 | The inverse: theta_max = lambda / (2 Delta1). |
| (7.8) | 118 | 131 | theta_edges = (D1 + D2) / (2 Delta_z). The source-edge to receiver-edge angle. |
| (7.9) | 118 | 131 | theta_k, the tilt of the virtual spherical wave at the source edge. |
| (7.13) | 118 | 131 | The combined geometric limit before the algebra. |
| **(7.14)** | 119 | 132 | **CONSTRAINT 1**: Delta2 <= -(D2/D1) Delta1 + lambda Delta_z / D1. It samples every ray that reaches the region of interest. |
| (7.15), (7.16) | 119 | 132 | D_illum, the diameter that the source illuminates in the observation plane. |
| (7.17), (7.18) | 119–120 | 132–133 | The grid must be at least the mean of D_illum and D2. Then the wrap-around stops at the aperture edge. |
| **(7.20)** | 120 | 133 | **CONSTRAINT 2**: N >= D1/(2 Delta1) + D2/(2 Delta2) + lambda Delta_z / (2 Delta1 Delta2). |
| (7.21) | 120 | 133 | The fixed one-step pitch, repeated from (6.16). |
| (7.25) | 121 | 134 | The minimum N for one-step Fresnel propagation. It diverges as lambda Delta_z approaches D2 Delta1. |
| (7.31) | 121 | 134 | Constraint 2 reduces to the SAME bound for one step. A consistency check. |
| (7.32) | 122 | 135 | The model source: an apodized amplitude times a parabolic wavefront of radius R. R < 0 diverges, R > 0 converges. |
| (7.36) | 122 | 135 | The one-step chain with the source curvature folded into one Q. |
| (7.37) | 122 | 135 | The local spatial frequency: f_loc = grad(phase) / (2 pi). |
| (7.40) | 123 | 136 | The Nyquist test on the source quadratic phase, at the aperture edge x = D1/2. |
| (7.41), (7.42) | 123 | 136 | The MINIMUM one-step distance. Finite R: Delta_z >= D1 Delta1 R / (lambda R - D1 Delta1). Infinite R: Delta_z >= D1 Delta1 / lambda. The book calls this a guideline, not a rule. |
| (7.45) | 125 | 138 | The angular-spectrum chain with the source curvature folded in. |
| (7.46), (7.47) | 125 | 138 | The TWO quadratic phase factors that the angular-spectrum method must sample. |
| (7.52) | 126 | 139 | The Nyquist test on the first factor. |
| **(7.53)** | 126 | 139 | **CONSTRAINT 3**: (1 + Delta_z/R) Delta1 - lambda Delta_z / D1 <= Delta2 <= (1 + Delta_z/R) Delta1 + lambda Delta_z / D1. |
| (7.58) | 126 | 139 | The Nyquist test on the transfer function, at the edge of the frequency grid. |
| **(7.59)** | 127 | 140 | **CONSTRAINT 4**: N >= lambda Delta_z / (Delta1 Delta2). It depends on the METHOD, not the geometry. Chapter 8 relaxes it. |
| (constraint set) | 127 | 140 | The four constraints listed together. Solve them in the (Delta1, Delta2) plane. |
| (7.60) | 129 | 142 | Constraint 3 does not apply when 1 + Delta_z/R < D2/D1. The geometric beam then stays inside D2. |

## Chapter 8 — Relaxed Sampling Constraints with Partial Propagations (printed 133–147, pdf 146–160)

| Section | Title | printed | pdf |
| --- | --- | --- | --- |
| 8.1 | Absorbing Boundaries | 134 | 147 |
| 8.2 | Two Partial Propagations | 135 | 148 |
| 8.3 | Arbitrary Number of Partial Propagations | 138 | 151 |
| 8.4 | Sampling for Multiple Partial Propagations | 139 | 152 |
| 8.5 | Problems | 146 | 159 |

### Key equations of Chapter 8

| Eq. | printed | pdf | Gloss |
| --- | --- | --- | --- |
| (8.1) | 134 | 147 | The SUPER-GAUSSIAN absorbing boundary: g = exp(-(r/sigma)^n), n > 2. |
| (8.2) | 134 | 147 | The Tukey (cosine-taper) window. The alternative boundary. |
| (8.3)–(8.7) | 136 | 149 | The similar-triangle proof of the linear pitch rule. |
| (8.8) | 136 | 149 | The LINEAR pitch rule: Delta_i = (1 - alpha_i) Delta1 + alpha_i Delta_n, with alpha_i = z_i / Delta_z. |
| (8.9) | 137 | 150 | Two partial propagations, with the absorbing-boundary operator A in the middle plane. |
| (8.14), (8.15) | 138 | 151 | The two middle-plane quadratic phases are inverse. Their product is 1, so both drop out. |
| (8.16) | 138 | 151 | The reduced two-propagation chain. One quadratic phase at each END only. |
| (8.18) | 139 | 152 | The GENERAL form: a product over n-1 partial propagations, each with its own boundary, transform pair, and Q2. |
| (8.19)–(8.22) | 139, 143 | 152, 156 | Constraint 3 is UNCHANGED by partial propagation. The intermediate pitches cancel. |
| (8.23) | 143 | 156 | Constraint 4 per partial step: N >= lambda Delta_z_i / (Delta_i Delta_i+1). |
| **(8.24)** | 144 | 157 | **The step-length cap**: Delta_z_max = min(Delta1, Delta_n)^2 N / lambda. This is the rule that partial propagation buys. |
| (8.24 text) | 144 | 157 | The plane count: n >= ceil(Delta_z / Delta_z_max) + 1. More planes are always allowed. |
| (8.25) | 144 | 157 | The worked example: 66.7 um pitch, N = 128, lambda = 1 um gives Delta_z_max = 0.567 m, so n = 5 planes. |
| (procedure) | 144 | 157 | The three-step recipe: (1) pick N, Delta1, Delta_n from constraints 1 and 2; (2) get Delta_z_max and n from (8.24); (3) more steps are safe. |
| (Listing 8.1) | 142 | 155 | The reference boundary in code: `w = 0.47*N`, `sg = exp(-nsq.^8/w^16)`. That is a super-Gaussian of POWER 16 in r, with the half-width at 0.47 N pixels. |

## Chapter 9 — Propagation through Atmospheric Turbulence (printed 149–183, pdf 162–196)

| Section | Title | printed | pdf |
| --- | --- | --- | --- |
| 9.1 | Split-Step Beam Propagation Method | 149 | 162 |
| 9.2 | Refractive Properties of Atmospheric Turbulence | 150 | 163 |
| 9.2.1 | Kolmogorov theory of turbulence | 152 | 165 |
| 9.2.2 | Optical propagation through turbulence | 156 | 169 |
| 9.2.3 | Optical parameters of the atmosphere | 157 | 170 |
| 9.2.4 | Layered atmosphere model | 164 | 177 |
| 9.2.5 | Theory | 164 | 177 |
| 9.3 | Monte-Carlo Phase Screens | 166 | 179 |
| 9.4 | Sampling Constraints | 172 | 185 |
| 9.5 | Executing a Properly Sampled Simulation | 174 | 187 |
| 9.5.1 | Determine propagation geometry and turbulence conditions | 174 | 187 |
| 9.5.2 | Analyze the sampling constraints | 176 | 189 |
| 9.5.3 | Perform a vacuum simulation | 178 | 191 |
| 9.5.4 | Perform the turbulent simulations | 179 | 192 |
| 9.5.5 | Verify the output | 180 | 193 |
| 9.6 | Conclusion | 182 | 195 |
| 9.7 | Problems | 183 | 196 |

### Key equations of Chapter 9

| Eq. | printed | pdf | Gloss |
| --- | --- | --- | --- |
| (9.1) | 150 | 163 | The SPLIT STEP: propagate a half gap, apply the phase, propagate the other half. ASSUMPTION: the index deviation is small. |
| (9.2) | 150 | 163 | The refraction operator T = exp(-i psi), with psi the phase along the gap. |
| (9.3) | 150 | 163 | The full turbulent chain. It is (8.18) with T in place of the boundary operator. Set T = 1 to get vacuum. |
| (9.5), (9.7) | 153 | 166 | The velocity and potential-temperature structure functions. The 2/3 power law. |
| (9.13) | 154 | 167 | The refractive-index structure function. Two branches, split at l0. |
| (9.14) | 154 | 167 | Cn2 from the temperature structure constant. |
| (9.16) | 155 | 168 | The KOLMOGOROV PSD: Phi_n = 0.033 Cn2 kappa^(-11/3). Valid for 1/L0 << kappa << 1/l0. |
| (9.17) | 155 | 168 | The von Karman PSD. It adds the outer scale. |
| (9.18) | 155 | 168 | The MODIFIED von Karman PSD. It adds both scales. kappa_m = 5.92/l0, kappa_0 = 2 pi / L0. |
| (9.19) | 156 | 169 | The Taylor frozen-flow hypothesis. It converts space statistics to time statistics. |
| (9.31) | 158 | 171 | The coherence factor (modulus of the complex coherence factor). |
| (9.32), (9.33) | 158 | 171 | The wave structure function D = D_chi + D_phi = -2 ln(mu). |
| (9.35) | 158 | 171 | The mean atmospheric MTF from the phase structure function. |
| (9.38) | 159 | 172 | The plane-wave coherence factor for Kolmogorov turbulence. |
| (9.41) | 159 | 172 | The r0 definitions: D(r0) = 6.88 rad^2, and r0 = 2.1 rho_0. |
| (9.42) | 159 | 172 | r0 for a PLANE wave: (0.423 k^2 integral Cn2 dz)^(-3/5). |
| (9.43) | 159 | 172 | r0 for a SPHERICAL wave. The (z/Z)^(5/3) path weight. |
| (9.44) | 160 | 173 | D(r) = 6.88 (r/r0)^(5/3). The Kolmogorov phase structure function. Use it to verify a screen. |
| (9.45) | 160 | 173 | The von Karman structure function. Finite outer scale. |
| (9.46), (9.47) | 160 | 173 | The modified von Karman structure function, and the Andrews algebraic fit with under 2% error. |
| (9.48) | 160 | 173 | Phi_phi = 2 pi^2 k^2 z Phi_n. The link from the index PSD to the phase PSD. ASSUMPTION: plane wave, weak turbulence. |
| (9.49)–(9.51) | 160–161 | 173–174 | The phase PSDs in ANGULAR frequency: 0.49 r0^(-5/3) times the three spectrum shapes. |
| **(9.52)** | 161 | 174 | The phase PSD in ORDINARY frequency: Phi_phi(f) = 0.023 r0^(-5/3) f^(-11/3). This is the screen PSD that the code uses. |
| (9.53), (9.54) | 161 | 174 | The Fried MTF, long and short exposure. |
| (9.55) | 161 | 174 | The exposure switch: 0 long, 1 short without scintillation, 1/2 short with scintillation. |
| (9.57) | 162 | 175 | The Andrews Strehl fit: S = [1 + (D/r0)^(5/3)]^(-6/5). |
| (9.58) | 162 | 175 | The Sasiela Strehl polynomial. Accurate for D/r0 > 2. |
| (9.60), (9.61) | 163 | 176 | The isoplanatic angle: D_phi(theta0) = 1 rad^2, and its Cn2 integral. |
| (9.62) | 163 | 176 | The log-amplitude variance definition. |
| (9.63), (9.64) | 163 | 176 | The plane-wave and spherical-wave log-amplitude variance. **Weak fluctuations are sigma_chi^2 < 0.25.** Rytov theory holds only there. |
| **(9.65)** | 164 | 177 | The LAYER MOMENT MATCH: the layered Cn2 must match the continuous profile for the moments 0 <= m <= 7. This is the rule that decides where screens go. |
| (9.66)–(9.69) | 164–165 | 177–178 | The discrete-sum forms of r0 (plane and spherical) and sigma_chi^2 (plane and spherical). |
| **(9.70)** | 165 | 178 | The PER-SCREEN r0: r0_i = (0.423 k^2 Cn2_i Delta_z_i)^(-3/5). ASSUMPTION: it is the PLANE-wave r0, so the layer must be thin. |
| (9.71)–(9.74) | 165 | 178 | The composite r0 and Rytov quantities written from the screen r0 values. |
| (9.75) | 165 | 178 | The 5-screen matrix system. Row 1 holds alpha^(5/3), row 2 holds alpha^(5/6) (1-alpha)^(5/6). |
| (9.76)–(9.78) | 166–167 | 179–180 | The Fourier-series representation of a phase screen. |
| (9.79), (9.80) | 167 | 180 | The Fourier-coefficient variance: mean of abs(c)^2 = Phi_phi(f) Delta_fx Delta_fy = Phi_phi / (Lx Ly). |
| (9.81) | 169 | 182 | The SUBHARMONIC low-frequency screen. Np = 3 grids, a 3x3 frequency set, spacing 1/(3^p L). |
| (constraint text) | 172 | 185 | Johnston and Lane: sample the SCINTILLATION too. The scale is the Fresnel length sqrt(lambda z), so use half of it as a pitch cap. |
| **(9.84), (9.85)** | 173 | 186 | The turbulence-blurred apertures: D1' = D1 + c lambda Delta_z / r0_rev, D2' = D2 + c lambda Delta_z / r0. **c is 2 to 8. c = 2 holds 97% of the light, c = 4 holds 99%.** |
| (9.86)–(9.88) | 173–174 | 186–187 | Constraints 1, 2 and 3 REPEATED with D1' and D2'. Constraint 4 does not change, because it is a method rule. |
| (9.89), (9.90) | 174 | 187 | Delta_z_max = min(Delta1, Delta_n)^2 N / lambda, and n_min = ceil(Delta_z / Delta_z_max) + 1. |
| (Listing 9.5) | 175 | 188 | The per-screen Rytov CAP in code: `rmax = 0.1`. The book credits Martin and Flatte. It is the upper bound on one screen's share of the Rytov number. |
| (Sec. 9.5.2 text) | 177 | 190 | The 50 km example needs only 2 planes by the sampling rules. The book uses 11 planes "to represent the atmosphere properly". It gives NO formula for that floor. |

### What Chapter 9 does NOT give

- **No explicit minimum screen count.** The book gives an UPPER bound on one
  screen's Rytov share (`rmax = 0.1`, Listing 9.5 line 37, printed p. 175) and
  it says "a typical number of phase screens, like 5–10" (Sec. 9.2.5, printed
  p. 165). It gives no derived lower floor. This matters for WP7. See Table 3.
- **No annular (obscured) aperture treatment.** Same gap as the Andrews
  tracker.
- **No temporal axis in the worked example.** Section 9.5.4 (printed p. 179)
  states the frozen-flow method in prose and points to the Greenwood frequency,
  but it gives no equation and no code.
- **No spherical (co-moving) grid.** Section 6.4 names the Coles and Rubio
  angular-grid method and then does NOT develop it. The olb `LensFresnel` and
  `Convert` route has no Schmidt equation to check against.

## Appendices and back matter

| Part | Title | printed | pdf |
| --- | --- | --- | --- |
| A | Function Definitions | 185 | 198 |
| B | MATLAB Code Listings | 187 | 200 |
| — | References | 189 | 202 |
| — | Index | 195 | 208 |

---

# Glossary — book symbol to olb name

`printed p` gives where the book defines the symbol. An empty olb column means
that olb has no name for it yet.

| Book symbol (printed p) | Meaning, and the olb name |
| --- | --- |
| Delta1 (89) | Source-plane grid pitch. `GridSpec.pixel_m`; the `d1` argument in `olb/waveoptics/propagators.py`. |
| Delta2, Delta_n (89, 140) | Observation-plane grid pitch. olb keeps ONE pitch for the whole flat run, so `GridSpec.pixel_m` covers both. The scaled route changes the pitch inside `LensFresnel`. |
| Delta_i (140) | The pitch of the ith plane, from the linear rule (8.8). olb has no per-plane pitch: `split_step` holds one flat grid. |
| Delta_f1 (99) | Spatial-frequency pitch, 1/(N Delta1). Internal to the FFT calls. |
| N (91) | Grid points per side. `GridSpec.n`; the floor `N_MIN = 256` and the cap `QualityPreset.n_max` in `olb/waveoptics/grid.py` and `.../turbulence/sampling.py`. |
| L = N Delta (121) | Grid side length. `GridSpec.size_m`. |
| D1 (117) | Source maximum extent. Built by `_features` in `olb/waveoptics/grid.py` from the waist and the transmit aperture. |
| D2 (117) | Observation-plane region of interest. The receive `aperture_m`. |
| D1', D2' (173) | The turbulence-blurred extents. olb has no direct name. `side_core` in `sampling.turbulent_grid` plays the same role: `guard * 2 * r_beam + 2 * (lambda / r0_total) * z`. |
| c (173) | The blur sensitivity, 2 to 8. olb hard-codes the equivalent factor as 2. |
| m = Delta2/Delta1 (89) | The scaling parameter. `beam_magnification` in `olb/waveoptics/grid.py`; the flag `GridSpec.scaled`. |
| Delta_z, z (89) | Step distance, and total distance. `geometry.slant_range_m`; `z_total_m` in `ScreenPlan`. |
| alpha_i = z_i / Delta_z (140) | Fractional distance along the path. `ScreenPlan.z_m / ScreenPlan.z_total_m`. |
| n (planes) (139) | Plane count. `ScreenPlan.z_m.size`; the preset floor `QualityPreset.min_screens`. |
| R (122) | Source wavefront radius. `Field.curvature` in `olb/waveoptics/field.py`. |
| theta_max (117) | Maximum ray angle the grid can carry. No olb name. |
| Q, Q2, V, F, R (89) | The propagation operators. olb hides them inside `Forvard`, `Fresnel`, `GForvard`, `Lens`, `LensFresnel` and `Convert`. |
| T[z_i, z_i+1] (150) | The refraction (phase screen) operator. `Screen(Fin, phase_rad)` in `olb/waveoptics/turbulence/screens.py`. |
| g_sg (134) | The super-Gaussian absorbing boundary. `super_gaussian_boundary(n, width_frac, power)` in `.../turbulence/splitstep.py`. |
| sigma, n (window) (134) | The boundary half-width and power. olb uses `width_frac` and `power`. NOTE: olb states the width as a fraction of the half-side; the book states it in pixels. |
| Cn2 (154) | Refractive-index structure parameter. `olb/turbulence/profiles.py`; `ScreenPlan.cn2_int_m13` holds the INTEGRATED value per screen. |
| l0, L0 (153) | Inner and outer scale. `l0_m` and `L0_m` in `phase_screen`. |
| kappa_m = 5.92/l0, kappa_0 = 2 pi/L0 (155) | The two spectrum corner frequencies. Internal to the screen generator. |
| Phi_n(kappa) (155) | Refractive-index PSD. Held by aotools inside `phase_screen`. |
| Phi_phi(f) (161) | Phase PSD, 0.023 r0^(-5/3) f^(-11/3). The screen PSD. |
| r0 (159) | Fried parameter of the whole path. `ScreenPlan.r0_total_m`; `_composite_r0` in `sampling.py`. |
| r0_i (165) | Per-screen Fried parameter. `screen_r0(cn2_integral_m13, wavelength_m)`; `ScreenPlan.r0_m`. |
| rho_0 (159) | Coherence radius, r0 = 2.1 rho_0. No olb name. |
| theta_0 (163) | Isoplanatic angle. `olb/turbulence/andrews/paths.py`. |
| sigma_chi^2 (163) | LOG-AMPLITUDE variance. Weak below 0.25. **CAUTION**: `ScreenPlan.sigma2_r` and `QualityPreset.sigma2_r_screen_max` hold the plane-wave RYTOV variance, which is 4 sigma_chi^2. Do not compare the two numbers directly. |
| mu(r) (158) | Coherence factor. The verification target of Sec. 9.5.5. No olb name. |
| D(r), D_phi(r) (158) | Wave and phase structure functions. No olb name. |
| H(f) (158) | Mean atmospheric MTF. No olb name. |
| sqrt(lambda z) (172) | Fresnel length, the scintillation scale. olb has no NAME for it, but it HAS the rule: `olb/waveoptics/turbulence/sampling.py:454` caps the pitch at `sqrt(lambda z_i)/2`. olb cites Andrews Ch. 8 for it; the same rule is Schmidt Sec. 9.4, printed p. 172, from Johnston and Lane. So the olb pitch rule uses BOTH `pixels_per_r0` and the Fresnel scale. See Table 2, row S-26. |
| rect (185) | Appendix A. No olb equivalent. |
| tri (185) | Appendix A. No olb equivalent. |
| sinc (185) | Appendix A. Used by the model point source. No olb equivalent. |
| comb (185) | Appendix A. No olb equivalent. |
| circ (185) | Appendix A. `CircAperture` and `CircScreen` in `olb/waveoptics/sources.py`. |
| jinc (185) | Appendix A. No olb equivalent. |

---

# Table 1 — forward map (olb code to book)

| olb id | location (file:line) | quantity | book eq | printed p | pdf p | status | note |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `Power`, `Normal` | `olb/waveoptics/field.py:181` | Riemann-sum scaling, `sum(I) * dx^2` | (2.3), (2.32) | 15, 36 | 28, 48 | checked | The same rule that puts `dx^2` on `ft2`. The two agree. The book cites it for the transform; olb cites Goodman for the power. No change needed. |
| `Forvard`, `Fresnel` | `olb/waveoptics/propagators.py:121` | UNSCALED `fft2`/`ifft2` | (2.6), (2.9) | 16, 17 | 29, 30 | checked | The LightPipes propagators use a bare `fft2` with a sign-alternation trick in place of `fftshift`. The forward and the inverse transform cancel inside one propagator, so the missing `dx^2` and `(N df)^2` cancel too. The result is correct, but the intermediate spectrum carries NO physical scaling. Do NOT mix `schmidt.fourier.ft2` with those intermediates. |
| `forvard_max_z` | `olb/waveoptics/grid.py:209` | `z_max = N dx^2 / lambda` | (7.59) | 127 | 140 | checked | EXACT. The rule is constraint 4 of Ch. 7, inverted with Delta2 = Delta1. The docstring cites "Ch. 6", which is the wrong chapter. `schmidt/sampling.py` `angular_spectrum_max_z` reproduces it, and the self-check proves the two agree. See gap S-12. |
| the flat EXTENT rule | `olb/waveoptics/grid.py:105` | `size = guard * 2 * r_max`, the grid side of a vacuum propagation | (7.18) | 120 | 133 | conflict | DIFFERENT RULE. The book sizes the grid from the illuminated area and the region of interest, `D_grid >= (D_illum + D2)/2`, and it ALLOWS the wrapped light outside D2. olb puts a fixed margin around the beam. See gap S-07. |
| the RESOLUTION rule | `olb/waveoptics/grid.py:118`, `:166` | `dx <= feature/(P/2)`, the pixel pitch | Listing 7.1; Ch. 8 text | 124, 144 | 137, 157 | partial | The book gives no equation. It picks the spacing per example: "at least 50 grid pts across ap" (Listing 7.1), and at least 30 grid points across D1 and across D2 (Ch. 8 text). COMPATIBLE IN FORM, COARSER IN VALUE. See Table 3. |
| `n_wanted` | `olb/waveoptics/grid.py:169` | `n = 2 ** ceil(log2(size/dx))`, the pixel count | Listing 7.1 line 11; Listing 7.2 line 13 | 124, 128 | 137, 141 | checked | EXACT. The book rounds N up to the next power of two for the FFT. |
| the range warning | `olb/waveoptics/grid.py:179` | warns when `z > forvard_max_z` | (7.59) | 127 | 140 | checked | EXACT rule, wrong citation. The docstring says Ch. 6. See gap S-12. |
| the turbulent extent rule | `olb/waveoptics/turbulence/sampling.py:370`, `:441` | `side = [guard*2*r_beam + 2*(lambda/r0)*z] / (1 - b)` | Ch. 9 (the scattering cone); (8.1) for the band | 134 | 147 | out of scope | Chapters 7 and 8 give no such rule, except the `(1 - b)` absorbing-band divisor of Sec. 8.1. Chapter 9 owns the rest of the row. See Table 3. |
| the turbulent pixel rule | `olb/waveoptics/turbulence/sampling.py:382`, `:451` | `dx <= min(r0/P_r0, sqrt(lambda z_i)/2, feature/(P/2))` | — | — | — | gap | NO Ch. 7 CONTENT. The two scale rules come from Ch. 9 (r0) and from the Fresnel scale. The turbulent sizer never evaluates constraints 1 to 4. See gap S-06. |
| the pixel-count clamp | `olb/waveoptics/turbulence/sampling.py:457` | `n = clamp(2**ceil(log2(side/dx)), 256, n_max)` | Listing 7.2 line 13 | 128 | 141 | partial | The power-of-two step is EXACT. The `[256, n_max]` clamp has no book source. See Table 3. |
| `step_over_limit_max` | `olb/waveoptics/turbulence/sampling.py:469` | `max(gap)/forvard_max_z`, the worst split-step length against constraint 4 | (8.24) | 144 | 157 | partial | SAME IDEA, DIFFERENT ROUTE. The book SETS the plane count from Eq. (8.24). olb sets it from the Cn2 profile and only REPORTS the ratio. See gap S-10. |
| `super_gaussian_boundary` | `olb/waveoptics/turbulence/splitstep.py:29` | the absorbing boundary | (8.1); Listing 8.1; Fig. 8.1 | 134, 142 | 147, 155 | conflict | FORM AGREES (a super-Gaussian of exponent above 2). The NUMBERS are a different parameterisation. See Table 3 and gap S-11. |
| `Forvard` | `olb/waveoptics/propagators.py:57` | The angular-spectrum propagator, m = 1 | (6.31), (6.32) | 95 | 108 | checked | The transfer function in the docstring, `exp(i k z) exp(-i pi lam z f^2)`, is EXACTLY Eq. (6.32). olb KEEPS the piston factor `exp(i k z)`; `schmidt.fresnel.angular_spectrum` drops it, as Listing 6.5, printed p. 102, does. So the two differ by one constant phase. The irradiance is the same. |
| `Forvard` | `olb/waveoptics/propagators.py:121` | The transform pair | (2.6), (2.9) | 16, 17 | 29, 30 | checked | Duplicate of the WP1 row. A bare `fft2`/`ifft2` with a sign-alternation trick in place of `fftshift`. The missing `dx^2` and `(N df)^2` cancel inside one propagator. |
| `Forvard` grid rule | `olb/waveoptics/propagators.py:69` | "give the grid a side of about 8 times the largest beam radius" | (7.14), (7.20) | 119, 120 | 132, 133 | conflict | The olb rule is a rule of thumb with no source. Constraints 1 and 2 give the real bound, and both need D1, D2 and z. See gap S-16. |
| `Fresnel` | `olb/waveoptics/propagators.py:136` | The CONVOLUTION form of the Fresnel integral | (6.6) | 88 | 101 | checked | olb convolves on a grid of twice the side and integrates the kernel over one pixel in closed form (C and S). The book does NOT do that: it multiplies by the analytic transfer function of Eq. (6.49), printed p. 99. The two solve the same Eq. (6.6). The pixel integral is a REFINEMENT of the book, not a departure from it. |
| `Fresnel` minimum distance | `olb/waveoptics/propagators.py:149` | "z comparable with, or less than, the aperture" | (7.41), (7.42) | 123 | 136 | checked | The docstring cites Ch. 7 for the minimum distance. Eqs. (7.41) and (7.42) give the number: `z >= D1 dx1 R /(lam R - D1 dx1)`, and `z >= D1 dx1 / lam` for a flat source. olb states the rule in words only. See gap S-17. |
| `Fresnel` doubled grid | `olb/waveoptics/propagators.py:139` | The zero-padded convolution | — | — | — | no book equation | The book has NO zero-padded convolution. It controls the wrap with the ABSORBING BOUNDARY of Ch. 8, Eq. (8.1), printed p. 134, and with the grid-size constraint of Eq. (7.20). Two different cures for one problem. |
| `GForvard` | `olb/waveoptics/propagators.py:253` | The analytic ABCD Gaussian route | (6.70), (6.77), (6.81) | 103, 104 | 116, 117 | partial | The book gives the ray matrices and the generalized Huygens-Fresnel (ABCD) INTEGRAL. It does NOT give the closed-form q-parameter transform `q2 = (A q1 + B)/(C q1 + D)` that `_ABCD` uses. olb cites Siegman, ISBN 978-0935702118. Eq. (6.77) assumes azimuthal symmetry, which a pure Gaussian meets. |
| `Lens` | `olb/waveoptics/lenses.py:70` | The thin lens as a quadratic phase | (6.76) | 104 | 117 | checked | Eq. (6.76) is the thin-lens ray matrix. The phase form is the operator `Q[-1/f, r]` of Eq. (6.7), printed p. 89. The book states the lens phase delay at Sec. 6.5, printed p. 104. |
| `LensFresnel` | `olb/waveoptics/lenses.py:173` | The co-moving (spherical) grid | — | — | — | NO book equation | The book NAMES the Coles and Rubio angular-grid method (Ch. 6, text, printed p. 87) and then does NOT develop it. Schmidt's own answer to the same problem is the SCALING PARAMETER m of Eq. (6.65), printed p. 100, on a FLAT grid. The olb `Lens -> LensFresnel -> Convert` recipe has no Schmidt equation to check against. See gap S-13. |
| `Convert` | `olb/waveoptics/lenses.py:233` | Return from spherical to a flat grid | — | — | — | NO book equation | Same as `LensFresnel`. The book never leaves the flat grid. |
| `LensFresnel` magnification | `olb/waveoptics/lenses.py:187` | `fA = z/(m - 1)`, grid scale `(f - z)/f` | (6.24), (6.52) to (6.54) | 94, 99, 100 | 107, 112, 113 | partial | The olb `m` and the book `m` mean the SAME thing: the ratio of the two grid pitches. But olb gets it from a virtual lens and the book gets it from a free parameter in the exponent. The two routes are not the same algorithm. |
| `split_step` | `olb/waveoptics/turbulence/splitstep.py:94` | The partial-propagation loop | (8.18) | 139 | 152 | conflict | olb calls `Forvard` (m = 1) on ONE flat grid for every step. Eq. (8.18) gives each step its OWN pitch, from the linear rule of Eq. (8.8), printed p. 136, and its own magnification m_i. So olb cannot grow the grid with the beam. See gap S-14. |
| `super_gaussian_boundary` | `olb/waveoptics/turbulence/splitstep.py:29` | The absorbing boundary | (8.1) | 134 | 147 | conflict | The FORM matches Eq. (8.1), but the parameterisation and the shape do not match Listing 8.1. See Table 3 and gap S-15. |
| `forvard_max_z` | `olb/waveoptics/grid.py:209` | `z_max = N dx^2 / lam` | (7.59), (8.24) | 127, 144 | 140, 157 | checked | This IS constraint 4 of Eq. (7.59), rearranged for m = 1: `N >= lam z/(dx1 dx2)` with `dx1 = dx2 = dx` gives `z <= N dx^2/lam`. It is also Eq. (8.24), the step cap, with `min(dx1, dxn) = dx`. The Table 3 row is now filled: the constant is DERIVED, not a guess. |
| `phase_screen` PSD | `olb/waveoptics/turbulence/screens.py:81` | modified von Karman PHASE PSD, `0.023 r0^(-5/3) exp(-(f/fm)^2) / (f^2+f0^2)^(11/6)` | (9.51), (9.52) | 161 | 174 | checked | The expression is the book's, and `fm = 5.92/(2 pi l0)`, `f0 = 1/L0` are the book's too. olb defaults `l0 = 1e-6 m` and `L0 = 1e6 m`, not 0 and infinity, so the numbers are the Kolmogorov ones to 12 digits. The constant agrees: `0.49 (2 pi)^(-5/3) = 0.02290` against the printed 0.023 (0.42% apart). |
| `phase_screen` FT draw | `olb/waveoptics/turbulence/screens.py:81` (aotools `ft_phase_screen`) | the Fourier-series screen | (9.78)–(9.80), Listing 9.2 | 167 | 180 | checked | The `ft_phase_screen` of the new module and the aotools one give the SAME mean structure function to 3 digits: the ratio to Eq. (9.44) is 0.636 / 0.540 / 0.422 / 0.347 / 0.286 at r/r0 = 1 / 2 / 4 / 6 / 8 (24 screens, N = 256, r0 = 10 px, direct-difference estimator). So the aotools FT screen IS Listing 9.2. |
| `phase_screen(subharmonics=True)` | `olb/waveoptics/turbulence/screens.py:81` (aotools `ft_sh_phase_screen`) | the subharmonic low-frequency screen | (9.81), Listing 9.3 | 169, 170 | 182, 183 | conflict | The two do NOT agree. On the run above, the book form gives 0.863 / 0.826 / 0.783 / 0.760 / 0.741 of theory, and aotools gives 0.824 / 0.778 / 0.725 / 0.692 / 0.661. The book form is 5 to 12% closer to Eq. (9.44) at every separation. See Table 2, row S-27. |
| `screen_r0` | `olb/waveoptics/turbulence/screens.py:56` | `r0_i = (0.423 k^2 Cn2_i dz_i)^(-3/5)` | (9.70) | 165 | 178 | checked | Exact match, the constant included. olb cites Fried and Andrews Ch. 12; the book credits Roggemann et al., DOI 10.1364/AO.34.004037. Add the Schmidt citation. |
| `_composite_r0` | `olb/waveoptics/turbulence/sampling.py:198` | `r0 = (SUM r0_i^(-5/3))^(-3/5)` | (9.71) | 165 | 178 | checked | Exact match. It is the PLANE-wave composite. The book also gives the spherical one, Eq. (9.72), which olb has no name for. |
| `_screen_rytov` | `olb/waveoptics/turbulence/sampling.py:175` | one screen's path weight | (9.63), (9.73) | 163, 165 | 176, 178 | checked | olb computes `2.25 k^(7/6) (INT Cn2 dz) (z - z_i)^(5/6)`, which is the plane-wave RYTOV variance `sigma_R^2`. The book's per-screen quantity is the LOG-AMPLITUDE variance `sigma_chi^2`, constant 0.563. The ratio is `2.25/0.563 = 3.997`. The self-check measures 3.9994. See Table 3. |
| `sigma2_r_screen_max` | `olb/waveoptics/turbulence/sampling.py:107` | the per-screen cap | Listing 9.5, lines 37, 38 | 175 | 188 | checked | The book caps `sigma_chi^2` at `rmax = 0.1`. olb caps `sigma_R^2 = 4 sigma_chi^2` at 0.05 / 0.10 / 0.25. See Table 3 for the factor analysis. |
| the extent rule, the scattering cone | `olb/waveoptics/turbulence/sampling.py:442` | `2 (lambda/r0) z` added to the grid side | (9.84), (9.85) | 173 | 186 | checked | The added term is `c lambda dz / r0` with `c = 2`, which is the book's low value. Listing 9.6, line 2, printed p. 177, uses `c = 2` too. The book states that `c = 2` holds 97% of the light and `c = 4` holds 99% (text below Eq. (9.85), printed p. 173). BUT olb adds the blur to the grid SIDE. The book adds it to D1' and D2' and then feeds constraints 1 to 3. Different route, same constant. |
| the pixel rule, `pixels_per_r0` | `olb/waveoptics/turbulence/sampling.py:451` | `dx <= r0_total / pixels_per_r0` | Sec. 9.4 text | 172 | 185 | checked | The book gives the rule of Johnston and Lane, DOI 10.1364/AO.39.004761: pick the pitch at which the phase step between two adjacent samples stays below pi for more than 99.7% of the draws. With Eq. (9.44) that reads `3 sqrt(6.88 (dx/r0)^(5/3)) <= pi`, so `dx <= 0.332 r0`, that is **3.01 pixels per r0**. The olb `standard` preset value 3 lands on it. |
| the pixel rule, the Fresnel scale | `olb/waveoptics/turbulence/sampling.py:454` | `dx <= sqrt(lambda z)/2` | Sec. 9.4 text | 172 | 185 | checked | olb ALREADY has the book's scintillation pitch rule, exactly. It cites Andrews Ch. 8 for it. The rule is Schmidt Sec. 9.4, printed p. 172, from Johnston and Lane. The tracker glossary row `sqrt(lambda z) (172)` said that olb has no such rule. That row was WRONG, and it is now corrected. See Table 2, row S-26. |
| `_merge_layers` | `olb/waveoptics/turbulence/sampling.py:311` | where the screens go, and what each carries | (9.65) | 164 | 177 | partial | olb groups adjacent Cn2 layers under the Rytov cap, and WP7 replaced the bail-out: a weak path now clamps to EXACTLY `min_screens` equal-weight groups. It does not SOLVE Eq. (9.65), but the Cn2-weighted centroid holds all 8 moments of the default profile inside 1 percent; the module self-check measures that. See Table 2, row S-22, and the WP7 note. |
| `turbulent_grid` | `olb/waveoptics/turbulence/sampling.py:366` | the grid sizer | (9.86)–(9.88) | 173, 174 | 186, 187 | gap | olb applies NONE of the three turbulent geometry constraints. It sizes the side from a beam-plus-cone rule and the pixel from `r0` and the Fresnel scale, then it rounds N up to a power of two. See Table 2, row S-21. |
| `super_gaussian_boundary` | `olb/waveoptics/turbulence/splitstep.py:29` | the absorbing boundary | (8.1); Listing 9.7, line 19 | 134, 179 | 147, 192 | checked | Eq. (8.1) gives the SHAPE `exp(-(r/sigma)^n)`, `n > 2`, and no numbers. Listing 9.7 gives the numbers, and they are not olb's. See Table 3. |
| `split_step` max hop | `olb/waveoptics/turbulence/splitstep.py:170` | `max_step = N dx^2 / lambda` | (9.89) | 174 | 187 | checked | Exact match. It repeats Ch. 8, Eq. (8.24), printed p. 144. olb cites Ch. 6; the turbulent statement is Eq. (9.89). |
| `split_step` loop | `olb/waveoptics/turbulence/splitstep.py:94` | the split-step chain | (9.1)–(9.3) | 150 | 163 | checked | olb hops to a screen, applies the screen, and hops on. The book applies the screen AT each partial-propagation plane, Eq. (9.3), printed p. 150. The two agree when the screens sit at the plane positions. The olb screens sit at slab CENTRES, so the two differ by half a slab. The book does not treat that case. |
| `min_screens` | `olb/waveoptics/turbulence/sampling.py:110` | the screen-count floor | — | — | — | olb rule | Chapter 9 gives NO such floor. WP7 kept 15 / 9 / 5 and it re-sourced them to an olb convergence sweep, with the moment floor of 4 as the absolute lower bound. See Table 3 and the WP7 note. |

# Table 2 — gaps and suggestions

| gap id | book section | book eq | capability | target module | priority |
| --- | --- | --- | --- | --- | --- |
| S-01 | 3.3 | (3.15), (3.17), (3.25) | The structure function of a phase screen. olb generates screens but never verifies one against D(r) = 6.88 (r/r0)^(5/3), Eq. (9.44). | `olb/waveoptics/schmidt/fourier.py` (built), then a screen check in WP4 | high |
| S-02 | 3.2 | (3.11), (3.14) | The windowed auto-correlation, which gives the coherence factor mu(r) of Sec. 9.5.5. It is the OTHER verification of a turbulent run. | a later work package | medium |
| S-03 | 2.5.2, 2.5.3 | (2.27), (2.31) | The p-fraction bandwidth of a Gaussian, with and without a quadratic phase. It gives a grid pitch from the beam alone. `GridSpec.for_scenario` uses a `pixels_per_feature` guess instead. | `olb/waveoptics/grid.py` | low |
| S-04 | 3.4 | (3.26) | The derivative by transform. The book states it is not used again. olb has no wavefront sensor. Do NOT build it. | none | none |
| S-05 | 6.3.2 | (6.24)–(6.29) | No TWO-STEP Fresnel propagator. `olb/waveoptics/propagators.py` has `Forvard`, `Fresnel` and `GForvard`, and `lenses.py` has the co-moving route. None of them frees the magnification with a second Fresnel integral. `schmidt/sampling.py` `two_step_planes` gives the two intermediate-plane geometries. A propagator that uses them is not built. | `olb/waveoptics/propagators.py` | — |
| S-06 | 7.3.3 | — | No budget, sizer, or runner calls a sampling checker. `GridSpec.for_scenario` and `turbulent_grid` warn on their OWN rules only. Wire `check_sampling` into `GridSpec.for_scenario` and into `SamplingReport`. It never raises, so it cannot break an existing run. | `olb/waveoptics/grid.py`, `olb/waveoptics/turbulence/sampling.py` | — |
| S-07 | 7.2, 7.3 | (7.14), (7.20) | Constraints 1 and 2 are implemented NOWHERE in olb. The observation-region extent D2 never enters a grid decision. olb uses the receive aperture as a FEATURE (a pixel rule), never as D2 (an extent rule). | `olb/waveoptics/grid.py` | — |
| S-08 | 7.3.2 | (7.53) | Constraint 3 is not checked. The transmit-beam curvature R never reaches a grid rule. This is the same missing curvature thread as Gap 3 of the Andrews cross-check. | `olb/waveoptics/grid.py` | — |
| S-09 | 7.3.1.2 | (7.41), (7.42) | The Fresnel-integral MINIMUM distance is not checked. `olb/waveoptics/propagators.py:136` `Fresnel` has no near-distance guard, so a short call aliases silently. `fresnel_min_distance` gives the bound. | `olb/waveoptics/propagators.py` | — |
| S-10 | 8.4 | (8.23), (8.24) | No vacuum partial-propagation planner. A long vacuum link takes the co-moving lens route instead of a chain of angular-spectrum steps. `partial_max_step` and `partial_plane_count` give the count. The turbulent planner sets its count from Cn2, not from Eq. (8.24). | `olb/waveoptics/run.py` | — |
| S-11 | 8.1 | Listing 8.1; Fig. 8.1 | The absorbing-boundary width carries no Schmidt number. `absorbing_boundary_sigma` gives the book value. See Table 3. | `olb/waveoptics/turbulence/splitstep.py` | — |
| S-12 | 7.3.2 | (7.59) | `forvard_max_z` cites "Ch. 6". The rule is Ch. 7, Eq. (7.59), printed p. 127. A one-line docstring fix. | `olb/waveoptics/grid.py` | — |
| S-13 | 6.4, 6.5 | (6.65) | The SCALED flat-grid propagator, as the alternative to the olb `Lens -> LensFresnel -> Convert` co-moving route. `schmidt.fresnel.angular_spectrum(..., dx2=)` now gives it. The two routes solve the same problem, and NO test compares them. Compare them on the 600 km uplink of `olb/waveoptics/run.py`. | an example in WP5 | high |
| S-14 | 8.3 | (8.18), (8.8) | The PER-PLANE grid pitch in the turbulent split step. `olb/waveoptics/turbulence/splitstep.py` holds one flat pitch for the whole path, so a diverging uplink beam must fit the SOURCE grid at the RECEIVER. `schmidt.fresnel.partial_propagations` now shows the book's linear pitch rule. Wiring it into `splitstep.py` is an owner decision, because it moves every turbulent number. | `olb/waveoptics/turbulence/splitstep.py` | high |
| S-15 | 8.1 | (8.1), Listing 8.1 | The absorbing boundary of olb is a DIFFERENT shape from the book's. olb: power 8, a taper band of 0.125 of the half-side, so the mask is exactly 1.0 out to 0.875 of the half-side and `exp(-1)` at the middle of an edge. Book: power 16, sigma = 0.47 N pixels, so the mask is 0.99999 at 0.2 N and 0.0678 at the middle of an edge. The book's boundary bites HARDER at the edge and it has no flat-then-taper break. Decide which one olb keeps. | `olb/waveoptics/turbulence/splitstep.py` | medium |
| S-16 | 7.2, 7.3 | (7.14), (7.20) | Constraints 1 and 2 as CODE. `olb/waveoptics/grid.py` sizes the grid from a `guard` of 4 and a `pixels_per_feature` of 16, with no citation. The book gives two inequalities in D1, D2, z, dx1 and dx2. | `olb/waveoptics/schmidt/sampling.py` (built in WP3), then `olb/waveoptics/grid.py` | high |
| S-17 | 7.3.1.2 | (7.41), (7.42) | The MINIMUM one-step distance as a number. `olb/waveoptics/propagators.py:149` states the rule in words: "z comparable with, or less than, the size of the aperture". The book gives `z >= D1 dx1 / lam` for a flat source. | `olb/waveoptics/schmidt/sampling.py` (built in WP3) | medium |
| S-18 | 6.6 | (6.82), (6.89), (6.92) | The MODEL point source. A true delta has infinite bandwidth, so the book replaces it with a sinc that gives the wanted windowed target field. olb has no point source: `olb/waveoptics/sources.py` gives `GaussBeam`, `PlaneWave`, `CircAperture` and `CircScreen` only. A retro link or a beacon may need one. | `olb/waveoptics/schmidt/` (a later work package) | low |
| S-19 | 6.5 | (6.77), (6.80), (6.81) | The general ABCD propagator for a NON-Gaussian field. `GForvard` handles a pure Gaussian only, and it raises on any other field. Eq. (6.81) gives the ABCD transfer function, which works on any field. | not needed yet | low |
| S-20 | 9.4 | Sec. 9.4 text, printed p. 172 | The Johnston and Lane PHASE pitch rule, `dx <= 0.332 r0`. The olb `pixels_per_r0` is a bare preset integer with a Martin and Flatte citation and no derivation. The book's prose plus Eq. (9.44) give the number. Built as `phase_pitch_max`. | `olb/waveoptics/turbulence/sampling.py` | medium |
| S-21 | 9.4, 9.5.2 | (9.86), (9.87), (9.88) | The three turbulent geometry constraints, and the blurred extents D1', D2' of Eqs. (9.84) and (9.85). olb checks none of them, so a bad pitch pair gives no warning. Built as `constraint1_pitch_max`, `constraint2_n_min`, `constraint3_pitch_range`, `blurred_extent`. | `olb/waveoptics/turbulence/sampling.py`, or a validation example | high |
| S-22 | 9.2.5 | (9.65) | The layered-atmosphere MOMENT rule for `0 <= m <= 7`. It is the only screen-placement rule that Chapter 9 gives. Built as `profile_moments`, `layer_moments`, `moment_error`. WP7 MEASURED it against the production grouping: the error stays inside 1 percent for every moment. The planner still does not SOLVE the rule, and `SamplingReport` does not carry the error. | `olb/waveoptics/turbulence/sampling.py` | low (WP7 measured it) |
| S-23 | 9.2.5, 9.5.1 | (9.75), Listing 9.5 | The constrained least-squares solve for the screen `r0` values from a target `r0_sw` and `sigma_chi,sw^2`. olb never solves for a screen strength; it takes the Cn2 layers as given. Built as `screen_strengths` and `max_screen_strength`. | `olb/waveoptics/turbulence/sampling.py` | medium |
| S-24 | 9.5.5 | (9.32), (9.44) | The observation-plane coherence factor as the end-to-end verification of a turbulent run. It is row S-02 of this table seen from Chapter 9. `properly_sampled_checklist` names it as an advisory step; nothing measures it. | a later work package | medium |
| S-25 | 9.4 | — | `QualityPreset.fresnel_weight_min` (`olb/waveoptics/turbulence/sampling.py:450`) exempts a weak screen from the Fresnel pitch rule. Chapter 9 states NO such exemption; Sec. 9.4 applies the rule to every step. The exemption is a real cost saver, so keep it, but mark it as an olb rule, not a book rule. | `olb/waveoptics/turbulence/sampling.py` | low |
| S-26 | 9.4 | — | Documentation only. The tracker glossary row for `sqrt(lambda z)` said that olb has no such rule. `olb/waveoptics/turbulence/sampling.py:454` has it exactly. DONE: the glossary row is corrected. | `docs/schmidt-crosscheck.md` | low (done) |
| S-27 | 9.3 | (9.81), Listing 9.3 | The aotools subharmonic screen reads 5 to 12% lower than the book's Listing 9.3 form across `r/r0 = 1` to 8. Both are below theory. Decide whether to keep aotools, to pass the book form from `schmidt.turbulence`, or to move to Johansson and Gavel, DOI 10.1117/12.177254, which the book calls the closest match (Ch. 9, text above Sec. 9.4, printed p. 172). | `olb/waveoptics/turbulence/screens.py` | medium |
| S-28 | 9.5.1 | Listing 9.5, lines 15 to 18 | A BOOK ERROR to record. Sec. 9.5.1, printed p. 176, prints `r0_sw = 17.7 cm` for the 50 km example. Listing 9.5 with the same inputs gives **12.66 cm**, and Problem 2, printed p. 183, confirms the `(3/8)^(-3/5)` factor. The printed `sigma_chi,sw^2 = 0.436` DOES reproduce (0.4365). So the printed `r0_sw` is the odd number. Do not calibrate anything against 17.7 cm. | none (a note) | — |

# Table 3 — constants ledger

Each row holds an olb constant that carries NO citation today. Fill the book
columns as the cross-check proceeds. An empty book column means that we have
not yet found a source, or that the book gives none.

The `book value` column holds the number that the book prints, when it prints
one. An empty book column means that we have not yet found a source, or that the
book gives none.

| olb constant | olb value | location | book quantity | book value | book eq | printed p | pdf p | status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `QualityPreset.min_screens` | 15 / 9 / 5 | `olb/waveoptics/turbulence/sampling.py:110` | a screen-count FLOOR | **the book would not give** | — | — | — | checked, no book source. WP7 sourced it to an olb convergence sweep, and it kept the values. See the WP7 note. |
| `QualityPreset.sigma2_r_screen_max` | 0.05 / 0.10 / 0.25 | `olb/waveoptics/turbulence/sampling.py:105` | `rmax`, the per-screen cap | 0.1 on `sigma_chi^2`, so **0.4 on `sigma_R^2`** | Listing 9.5, lines 37, 38 | 175 | 188 | checked, olb is 1.6x to 8x stricter. See the factor note below. |
| `QualityPreset.boundary_width_frac` | 0.125 / 0.125 / 0.10 | `olb/waveoptics/turbulence/sampling.py:105` | the boundary half-width | `0.47 N dx` (Listing 9.7); `0.45 L` (Fig. 8.1) | (8.1); Listing 9.7, line 19 | 134, 179 | 147, 192 | checked, the book gives no value in this parameterisation |
| `super_gaussian_boundary` `power` | 8 | `olb/waveoptics/turbulence/splitstep.py:29` | the super-Gaussian exponent n of Eq. (8.1) | **16** | (8.1); Listing 8.1, line 12; Listing 9.7, line 19 | 134, 142, 179 | 147, 155, 192 | CONFLICT. The book RUNS `sg = exp(-nsq.^8/w^16)`, which is `exp(-(r/w)^16)`, so n = 16. Figure 8.1, printed p. 134, also plots n = 16. Eq. (8.1) itself only needs n > 2, so the olb value of 8 is ALLOWED by the equation but it is not the book's number. The book also records that Flatte and others used n = 8 (Ch. 8, text, printed p. 134), so 8 has a source in the literature, not in Schmidt's own runs. |
| `super_gaussian_boundary` `width_frac` | 0.125 | `olb/waveoptics/turbulence/splitstep.py:29` | the half-width sigma of Eq. (8.1) | **0.47 N pixels** (Listing 8.1); **0.45 L** (Fig. 8.1) | (8.1); Listing 8.1, line 11 | 134, 142 | 147, 155 | CONFLICT of PARAMETERISATION. The book states one half-width sigma in PIXELS, measured from the centre. olb states a taper BAND width as a fraction of the half-side, with a hard flat region inside it. The two cannot be converted. Book at the middle of an edge: `exp(-(0.5/0.47)^16) = 0.0678`. olb at the middle of an edge: `exp(-1) = 0.368`. The book absorbs about 5 times harder there. |
| `QualityPreset.guard` | 4 / 3 / 2 | `olb/waveoptics/turbulence/sampling.py:105` | the beam-radius margin | the book gives no such margin; it sizes from D1', D2' and constraint 2 | (9.84)–(9.87) | 173, 174 | 186, 187 | checked, different route |
| `QualityPreset.pixels_per_r0` | 4 / 3 / 2 | `olb/waveoptics/turbulence/sampling.py:105` | the phase pitch rule | **3.01** pixels per r0 | Sec. 9.4 text with (9.44) | 172, 160 | 185, 173 | checked, `standard` = 3 matches |
| `QualityPreset.fresnel_weight_min` | 0.005 / 0.02 / 0.05 | `olb/waveoptics/turbulence/sampling.py:105` | a weak-screen exemption | **the book would not give** | — | — | — | checked, no source. See gap S-25. |
| `GridSpec.for_scenario` `guard` | 4.0 | `olb/waveoptics/grid.py:98` | the grid half-side margin, 4 beam radii | the book gives NO guard factor; Eq. (7.18) sizes the grid as `(D_illum + D2)/2` | (7.18) | 120 | 133 | UNCITED, and a DIFFERENT PHILOSOPHY. The book tolerates aliasing outside the region of interest; olb forbids it everywhere. The olb rule is stricter for a wide beam and it is silent about D2. See gap S-07. |
| `GridSpec.for_scenario` `pixels_per_feature` | 16 | `olb/waveoptics/grid.py:98` | points across the smallest hard edge | **50** across D1 (Listing 7.1); **30** across D1 and D2 (Ch. 8 text) | Listing 7.1; Ch. 8 text | 124, 144 | 137, 157 | UNCITED, and COARSER than both worked examples by a factor of 2 to 3. The book treats the number as a per-problem choice, not a constant. |
| `N_MIN` | 256 | `olb/waveoptics/grid.py:36` | a floor on the pixel count | the book has NO floor; its three worked examples land at **128**, **512** and **128** | Ch. 7 text; Ch. 8 text | 123, 128, 144 | 136, 141, 157 | UNCITED. All three book examples sit at or BELOW the olb floor, so the floor never binds on a book-sized problem. It is a convenience, not physics. |
| `forvard_max_z` | z_max = N dx^2 / lambda | `olb/waveoptics/grid.py:209` | constraint 4 with m = 1; the step cap | `N >= lam z/(dx1 dx2)`; `Delta_z_max = min(dx1,dxn)^2 N/lam` | (7.59), (8.24) | 127, 144 | 140, 157 | cited, checked. RESOLVED by WP2 and WP3: the olb formula IS the book formula. The constant is DERIVED, not a guess. Only the docstring chapter is wrong. See gap S-12. |
| `PIXELS_PER_FEATURE` | 8 | `olb/waveoptics/turbulence/sampling.py:56` | points across the smallest hard edge | **50** and **30**, as above | Listing 7.1; Ch. 8 text | 124, 144 | 137, 157 | UNCITED, coarser again by a factor of 4 to 6. |
| the scattering-cone factor | 2 (hard-coded) | `olb/waveoptics/turbulence/sampling.py:442` | `c`, the blur sensitivity | **2 to 8; c = 2 holds 97%, c = 4 holds 99%**; Listing 9.6 uses 2 | (9.84), (9.85) | 173 | 186 | checked, olb sits at the book's low end |
| `MAX_SCREENS` | 500 | `olb/waveoptics/turbulence/sampling.py:60` | a screen-count cap | **the book would not give** | — | — | — | checked, no source |

## The Fig. 8.1 and Listing 8.1 half-width

The Fig. 8.1 caption prints `sigma = 0.45 L` and `n = 16`, while Listing 8.1
prints `w = 0.47*N`. L and N are the same length in pixel units, so the two
numbers differ by 4%. `schmidt.fresnel.super_gaussian_absorber` takes 0.47 as
the default, because that is the value the book RUNS.

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

The olb `_screen_rytov` (`sampling.py:175`) computes
`2.25 k^(7/6) (INT Cn2 dz) (z - z_i)^(5/6)`. Substitute Eq. (9.70) into the
book's plane-wave term, Eq. (9.73), and it gives
`0.563 k^(7/6) (INT Cn2 dz) (z - z_i)^(5/6)`. The ratio is

    2.25 / 0.563 = 3.997

The self-check measures 3.9994. So **the olb per-screen number is
`sigma_R^2 = 4 sigma_chi^2`, and the book's cap of 0.1 on `sigma_chi^2` is a
cap of 0.4 on the olb number.**

| preset | olb cap on `sigma_R^2` | the same as a cap on `sigma_chi^2` | against the book's 0.1 |
| --- | --- | --- | --- |
| `reference` | 0.05 | 0.0125 | 8x stricter |
| `standard` | 0.10 | 0.025 | 4x stricter |
| `rapid` | 0.25 | 0.0625 | 1.6x stricter |

Two more differences, both small:

- The book weights the screen for a SPHERICAL wave,
  `alpha^(5/6) (1-alpha)^(5/6)`. olb weights it for a PLANE wave,
  `(1-alpha)^(5/6)`. For a downlink slab the source is far away, so the
  plane-wave weight is the right one, and it is the LARGER of the two. The olb
  choice stays conservative.
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

The olb `super_gaussian_boundary` is RADIAL, of order 8, exactly 1.0 inside
0.875 of the half-side, and `exp(-1)` at the edge. The two shapes are not the
same family, so the numbers do not map one to one. Recorded, not changed.
Eq. (8.1), printed p. 134, allows any order above 2 and gives no `sigma`.

---

# Work-package notes

## WP1 — Fourier foundations (Chs. 2, 3)

**Built.** The sub-package `olb/waveoptics/schmidt/` and its first module
`fourier.py`. The module gives four names:

- `ft2(g, dx)` — Ch. 2, Eq. (2.6), printed p. 16, with the two-dimensional
  scaling of Sec. 2.6, printed p. 36.
- `ift2(G, df)` — Ch. 2, Eq. (2.9), printed p. 17, with the two-dimensional
  scaling of Sec. 2.6, printed p. 37.
- `freq_pitch(n, dx)` — df = 1/(N dx). Ch. 2, text below Eq. (2.3), printed
  p. 16; Ch. 6, Eq. (6.51), printed p. 99.
- `structure_function(ph, mask, dx)` — Ch. 3, Eqs. (3.15), (3.17), (3.19) to
  (3.25), printed pp. 47 to 50.

**Decisions.**

- The shift pair is `fftshift(fft2(ifftshift(g)))`. The book prints
  `fftshift(fft2(fftshift(g)))` (Listing 2.5, printed p. 36). The two agree for
  an even grid count, which is the only count the book discusses (Sec. 2.1.3,
  printed p. 18). The form above is also correct for an odd count.
- `structure_function` uses the INVERSE transform of Listing 3.7, printed p. 48,
  not the forward transform of Eq. (3.25), printed p. 50. The bracket is real
  and even in f, so the two give the same result.
- `structure_function` does NOT multiply the result by the mask, and it needs no
  extra `dx^2`. The book's Listing 3.8, printed p. 48, divides its example by
  `delta^2`. A ramp check to a relative error of 1.8e-15 shows that the plain
  Eqs. (3.17) and (3.25) need no such factor. The listing's factor belongs to
  its example, not to the equation.
- The result is left unmasked. The caller selects the separation range where the
  overlap area A(dr) is not zero.
- Convolution (Sec. 3.1), correlation (Sec. 3.2), and the derivative (Sec. 3.4)
  are NOT built. Nothing needs them yet. See Table 2, rows S-02 and S-04.

**The book would not give.**

- The two-dimensional Gaussian transform pair. Section 2.5.2, printed pp. 31 and
  32, prints the ONE-dimensional pair, Eqs. (2.23) and (2.24), and the
  one-dimensional bandwidth, Eq. (2.27). The self-check builds the
  two-dimensional pair from Eq. (2.32), printed p. 36, and it writes the
  constants from the transform convention of Eq. (2.1), printed p. 15.
- A tolerance for the structure function. The book compares its figures by eye.
  The self-check sets its own targets.

## WP2 — Fresnel propagators (Ch. 6)

**Built.** `olb/waveoptics/schmidt/fresnel.py`. The module gives five names. It
imports numpy and `schmidt.fourier` only.

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

**Self-check numbers.**

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

The absorber values: the centre is 1.000000 exactly; the smallest value inside a
pixel radius of 0.2 N is 0.999999; the value at the middle of an edge is
0.067796, which equals the book value `exp(-(0.5/0.47)^16) = 0.067796` to 1e-12;
the mask never grows along a radius.

Every error is at the floor of double precision. That is the correct answer, not
a loose test: a Gaussian on a well-sampled grid has a spectrum far inside the
band, so the discrete Fresnel transform has no aliasing left to make an error.

The self-check also prints the step cap of Eq. (8.24) for the common geometry:
`Delta_z_max = 6.554 m`, so `n >= 3` planes. The run used 6.

**Decisions.**

- **The piston phase `exp(i k z)` is DROPPED from every kernel.** The book's own
  Listings 6.1, 6.3 and 6.5 (printed pp. 91, 96 and 102) drop it, and Listing
  8.1 (printed p. 142) drops it too. The factor is constant across a plane, so
  it changes no irradiance and no relative phase. To drop it in all four
  kernels is what lets them agree with each other to 1e-16. NOTE: the olb
  `Forvard` KEEPS `exp(i k z)`. Any future comparison must remove one or add the
  other.
- **The SCALED angular spectrum is BUILT, not skipped.** The task allowed the
  baseline `m = 1` form only. Eq. (6.65) reduces to Eqs. (6.31) and (6.32) when
  m = 1, so the general form with `dx2=None` as the default costs four lines
  and gives both. The tracker calls Eq. (6.50) "the workhorse of Chs. 7 to 9",
  and gap S-13 needs it to compare against the olb co-moving route.
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
- **The point source of Sec. 6.6 is NOT built.** See gap S-18.
- **The ABCD route of Sec. 6.5 is NOT built.** `GForvard` and `Lens` already
  cover the olb need. See gap S-19.

**The book would not give.**

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
  target. Through Chs. 6 to 8 it compares against Fresnel diffraction from a
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
  on a flat grid. See gap S-13.

## WP3 — Sampling constraints (Chs. 7, 8)

**Built.** `olb/waveoptics/schmidt/sampling.py`, plus a one-line placeholder
`olb/waveoptics/schmidt/__init__.py` (the package directory did not exist).

Small pure functions, numpy only, no olb import:

- Ch. 7.1 and 7.2, the band limit and the geometry: `nyquist_max_angle`
  (Eq. (7.7)), `geometric_max_angle` (Eqs. (7.8), (7.9), (7.12)),
  `illuminated_diameter` (Eq. (7.16)).
- The four numbered constraints: `constraint1_max_delta2` (Eq. (7.14)),
  `constraint2_min_n` (Eq. (7.20)), `constraint3_delta2_window` (Eq. (7.53)),
  `constraint4_min_n` (Eq. (7.59)), and `constraint3_is_slack` (Eq. (7.60)).
- The local-frequency analysis that the constraints come from:
  `local_spatial_frequency_source` (Eqs. (7.37), (7.39), (7.51)) and
  `local_spatial_frequency_transfer` (Eqs. (7.55), (7.57)).
- Per kernel: `one_step_delta2` (Eq. (7.21)), `one_step_min_n` (Eq. (7.25)),
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

The self-check also proves that the DERIVATIONS close, not only the numbers: the
Nyquist rule on Eq. (7.39) at the source edge gives back Eq. (7.42); the same
rule on Eq. (7.51) gives back the constraint-3 upper bound; the same rule on
Eq. (7.57) gives back constraint 4; Eq. (7.31) reproduces Eq. (7.25) exactly;
and Eq. (7.18) on `illuminated_diameter` reproduces constraint 2 exactly.

Run: `python -m olb.waveoptics.schmidt.sampling`.

**Decisions.**

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
- The module does NOT change any existing sizer. To wire it is gap S-06, an
  owner decision, because a wired check moves no numbers but it will print
  warnings on grids that run today.

**The book would not give.**

- NO guard factor and NO margin philosophy that matches `grid.py`. The book
  lets the wrapped light come half way around the grid, up to the edge of D2
  (Eq. (7.18), printed p. 120). It never asks for empty space around the beam.
  So `guard=4.0` cannot be justified from Ch. 7; it can only be replaced by
  constraint 2.
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

## WP4 — Turbulence and screens (Ch. 9)

**Built.** `olb/waveoptics/schmidt/turbulence.py`, twenty names, each with its
chapter, equation number and printed page:

- **The spectra (Secs. 9.2.3, 9.3).** `phase_psd` (the one shared expression),
  `kolmogorov_phase_psd` (9.49), (9.52), `von_karman_phase_psd` (9.50),
  `modified_von_karman_phase_psd` (9.51), `kolmogorov_structure_function`
  (9.44).
- **The screens (Sec. 9.3).** `ft_phase_screen` (9.78) to (9.80) with Listing
  9.2, `subharmonic_screen` (9.81) with Listing 9.3, `ft_sh_phase_screen` (the
  sum).
- **The per-screen bound (Sec. 9.2.5, Listing 9.5).** `screen_rytov_share`
  (9.73), (9.74), `max_screen_strength` (Listing 9.5, lines 37 to 39),
  `screen_strengths` (9.75), `composite_r0` (9.71), (9.72), `screen_r0` (9.70),
  and the constants `RMAX = 0.1` and `WEAK_SIGMA2_CHI = 0.25`.
- **The layer rule (Sec. 9.2.5).** `profile_moments`, `layer_moments`,
  `moment_error`, all Eq. (9.65).
- **The sampling bounds (Sec. 9.4).** `fresnel_pitch_max`, `phase_pitch_max`,
  `blurred_extent` (9.84), (9.85), `constraint1_pitch_max` (9.86),
  `constraint2_n_min` (9.87), `constraint3_pitch_range` (9.88),
  `max_partial_step` (9.89), `min_planes` (9.90).
- **The procedure (Sec. 9.5).** `properly_sampled_checklist`, which returns one
  `(rule, satisfied, bound, actual, citation)` tuple per step. Its arguments are
  plain numbers, and their names match `GridSpec` (`n`, `pixel_m`, `size_m`) and
  `ScreenPlan` (`z_m`, `r0_m`, `r0_total_m`, `z_total_m`) one to one. It imports
  no olb module outside `schmidt`.

**Measured.** Self-check numbers, from
`python -m olb.waveoptics.schmidt.turbulence` (7.5 s):

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
  `r0_sw = 12.66 cm` against the printed 17.7 cm; see Table 2, row S-28. The
  11-screen `screen_strengths` solve returns `r0_sw` to 1.2e-5 and
  `sigma_chi,sw^2` exactly, and its largest screen share is 0.0745 against the
  cap of 0.1.
- The factor between the book's per-screen quantity and the olb one: 3.9994.
- The Sec. 9.5.2 example: constraint 2 asks for `N >= 355.2`, and the book
  picks 512 because "the required number of grid points is more than 2^8"
  (printed p. 177). `min_planes` returns 2, and the book says two. Both match.
- The Johnston and Lane phase pitch rule gives 3.01 pixels per r0.

**Decisions.**

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

**The book would not give.**

- **A minimum screen count.** See the WP7 gate verdict below.
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

## WP5 — Examples

## WP6 — Retrofit and documentation

**Done.** Two parts, and neither one moved a number.

Part 1, the citation retrofit. Table 1 gave the map; the equation numbers now
sit in the production modules, in the docstrings and the comments only. The
bodies did not change, and the six self-checks
(`olb.waveoptics.grid`, `.propagators`, `.lenses`, `.turbulence.sampling`,
`.turbulence.splitstep`, `.turbulence.screens`) print the same numbers as
before. The changes:

- `grid.py` — `forvard_max_z` cited "Ch. 6". It is now Ch. 7, Eq. (7.59),
  printed p. 127, constraint 4 at m = 1, and Ch. 8, Eq. (8.24), printed p. 144.
  Gap S-12 is CLOSED. The guard, the pixels-per-feature value and `N_MIN` now
  say that they are olb rules, against Eq. (7.18) (gaps S-07 and S-16).
- `propagators.py` — `Forvard` cites Eqs. (6.31) and (6.32), and it records
  that it KEEPS the piston phase that the book listings drop. `Fresnel` cites
  the convolution form, Eq. (6.6), and it records that the pixel-integrated
  kernel is a refinement that the book does not use (the book transforms it,
  Eq. (6.49)). The doubled grid has no book equation. The minimum distance now
  carries Eqs. (7.41) and (7.42).
- `lenses.py` — the Ch. 7 hint is GONE. The module now records that the book
  names the Coles and Rubio angular grid and does not develop it, and that the
  book answer is the scaling parameter m of Eq. (6.65) on a flat grid (gap
  S-13). `Lens` keeps a real citation: Eqs. (6.7) and (6.76).
- `turbulence/sampling.py` — the `sqrt(lambda z)/2` pitch rule gains Sec. 9.4,
  printed p. 172; `pixels_per_r0` gains the 3.01 derivation;
  `sigma2_r_screen_max` gains the `rmax = 0.1` relation and the factor of 4;
  the scattering cone gains Eqs. (9.84) and (9.85) at c = 2;
  `fresnel_weight_min` is marked an olb rule; `_merge_layers` gains the
  Eq. (9.65) pointer.
- `turbulence/splitstep.py` — the absorber cites Eq. (8.1) and records the
  conflict with the book values (power 16, sigma 0.47 N). The sub-step note
  cites Eqs. (8.19) to (8.22). The screen placement at a slab centre is
  recorded against Eq. (9.3).
- `turbulence/screens.py` — the PSD cites Eq. (9.51); the screen draw cites
  Eqs. (9.78) to (9.80); the subharmonics cite Eq. (9.81), with Lane et al.,
  DOI 10.1088/0959-7174/2/3/003, and Johansson and Gavel,
  DOI 10.1117/12.177254. `screen_r0` cites Eq. (9.70), and `Screen` cites
  Eq. (9.2).

Part 2, the documentation. `docs/physics.md` gains Section 8, the twin of the
Andrews Section 5h. `docs/api-waveoptics.md` gains Section 10, the public
functions of the four modules, and the `QualityPreset` table gains the
`min_screens` caveat and the factor-4 note. `examples/schmidt/README.md` is
new. `docs/examples.md` gains the suite. `CLAUDE.md` records the layer and the
new `min_screens` evidence. `docs/backlog.md` gains the item 2-S1, a pointer to
Table 2.

## WP7 — The `min_screens` revision

### WP7 GATE VERDICT — what Chapter 9 does and does not justify

This verdict comes from WP4. It decides the `min_screens` revision.

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
  index and the coherence factor stop to move. That is the olb evidence, not
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
  200-layer real profile give the same target moments, and thus the same
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
  on a strong path, and the book's own Listing 9.5 avoids this: it FIXES the
  positions first and then solves only for the strengths under the cap. That
  is a defensible route for olb too, and it needs no moment machinery: pick the
  positions from the sampling rules, then solve Eq. (9.75).
- A generalized Gauss quadrature with the weight `Cn2(z)` would satisfy
  Eq. (9.65) exactly with 4 nodes for ANY profile. That is the clean
  generalisation of the self-check case. The book does not name it, and it is
  NOT built here. Flag it as the leading candidate for WP7, not as a decision.

**Practical recommendation for WP7, on the evidence above.** Replace the
`_merge_layers` bail-out with a rule that (i) picks a screen count from the
larger of `min_planes` (Eq. (9.90)) and the per-screen `rmax` cap, with a hard
floor of 4 from the moment count; (ii) places and weights the screens to
minimise `moment_error`; and (iii) reports the achieved moment error in
`SamplingReport`, because Chapter 9 gives the equality and no tolerance. Do
NOT keep 15 / 9 / 5 with a Schmidt citation. The book does not support it.

### WP7 — Built / Decisions / Measured

**Built.** `olb/waveoptics/turbulence/sampling.py`:

- `_equal_weight_groups(weights, n_groups)` — a new helper. It cuts the `Cn2`
  layers into EXACTLY `n_groups` contiguous groups of equal Rytov weight. The
  cut goes BEFORE the layer that passes the target share, so a layer that is
  heavier than one share keeps its own group; the strongest group is then the
  strongest LAYER, which is the least that any grouping can give.
- `_merge_layers` — the bail-out is gone. A natural merge that undershoots
  `min_screens` now calls `_equal_weight_groups`. Each screen keeps its
  `Cn2`-weighted centroid, and the Rytov cap still RAISES the count above the
  floor on a strong path.
- `_plan_space` — a new warning when the profile has fewer layers than
  `min_screens`. The planner does not split a layer, so the plan keeps the
  layer count and it says so.
- The self-check gains case 6 (the floor is a preset choice, not a profile
  artefact) and case 7 (the moment error of the grouping, through a
  validation-only import of `olb.waveoptics.schmidt.turbulence`). Case 5 moved
  to a 10 degree path, because the merged bottom group of a high-elevation plan
  is no longer exempt from the Fresnel pixel rule.

**Decisions.**

1. **The floor stays, and it keeps 15 / 9 / 5.** The gate verdict left two
   routes. WP7 took route (b), the convergence sweep. The docstring, the
   `PRESETS` comment, `docs/api-waveoptics.md` and `docs/physics.md` all name
   the source as an olb convergence sweep, NOT Schmidt.
2. **The absolute lower bound is 4.** It is the moment count of Eq. (9.65),
   printed p. 164: 8 equations against 2 free numbers for each screen. The
   self-check asserts that no preset goes under it.
3. **The planner does NOT solve Eq. (9.65).** The equal-weight grouping plus
   the `Cn2`-weighted centroid is the lazy route, and the measurement below
   shows that it is enough for the default profile. A generalized Gauss
   quadrature with the weight `Cn2(z)` stays the clean generalisation, and it
   stays unbuilt.
4. **`SamplingReport` does not carry the moment error.** The gate verdict asked
   for it. WP7 measures it in the self-check instead, because no caller reads
   it. Add the field when a caller needs it.

**Measured 1: the convergence sweep.** The sweep holds the GRID fixed (the grid
of the largest count), it moves the screen count only, and it runs 200 snapshots
for each count with one fixed seed. The metrics are the mean collected power in
dB and the aperture scintillation index `sigma2_P` of the collected power. The
Monte Carlo error of a variance from 200 snapshots is about 10 percent, so read
the `sigma2` columns to that accuracy. The space rows use a 60-layer `Cn2`
profile, so that 26 screens is reachable; the count is now independent of the
layer count, so the finer profile is free.

Case A, the 600 km downlink slab at 90 degrees, `rapid` preset, grid 256 px:

| screens | mean power, dB | delta, dB | `sigma2_P` | delta, % |
|---|---|---|---|---|
| 3 | 0.018 | 0.028 | 0.00149 | +3.9 |
| 5 | 0.001 | 0.011 | 0.00142 | -1.5 |
| 7 | 0.018 | 0.027 | 0.00145 | +0.7 |
| 9 | 0.021 | 0.031 | 0.00128 | -11.1 |
| 12 | 0.018 | 0.028 | 0.00131 | -8.9 |
| 15 | -0.008 | 0.002 | 0.00137 | -4.4 |
| 20 | -0.009 | 0.001 | 0.00143 | -0.5 |
| 26 | -0.010 | 0.000 | 0.00144 | 0.0 |

Case B, the same slab at 30 degrees, `rapid` preset, grid 256 px. THIS IS THE
BINDING CASE:

| screens | mean power, dB | delta, dB | `sigma2_P` | delta, % |
|---|---|---|---|---|
| 3 | 0.022 | 0.109 | 0.00846 | -19.0 |
| 5 | 0.021 | 0.108 | 0.00944 | -9.7 |
| 7 | -0.020 | 0.067 | 0.01069 | +2.3 |
| 9 | -0.011 | 0.076 | 0.01022 | -2.2 |
| 12 | -0.011 | 0.076 | 0.01049 | +0.4 |
| 15 | -0.083 | 0.004 | 0.01047 | +0.2 |
| 20 | -0.087 | 0.000 | 0.01037 | -0.8 |
| 26 | -0.087 | 0.000 | 0.01045 | 0.0 |

Case C, a 2 km horizontal path at `Cn2 = 3e-15` into a 30 mm sampling bucket,
`standard` preset, grid 256 px. The Rytov cap sets 4 as the smallest count that
the planner builds here:

| screens | mean power, dB | delta, dB | `sigma2_P` | delta, % |
|---|---|---|---|---|
| 4 | 8.431 | -0.042 | 0.07850 | -4.6 |
| 5 | 8.463 | -0.009 | 0.07992 | -2.9 |
| 7 | 8.461 | -0.011 | 0.08168 | -0.7 |
| 9 | 8.465 | -0.007 | 0.08411 | +2.2 |
| 12 | 8.447 | -0.026 | 0.08585 | +4.3 |
| 15 | 8.442 | -0.031 | 0.08663 | +5.3 |
| 20 | 8.461 | -0.012 | 0.08509 | +3.4 |
| 26 | 8.473 | 0.000 | 0.08230 | 0.0 |

**The verdict.** The tolerance is 0.1 dB on the mean power and 5 percent on the
index. The mean power meets it at every count in every case. The index meets it
at every count of case A and case C, because those two paths are weak and
homogeneous. Case B, the SLANT path, is the one that binds: the index is 19
percent low at 3 screens, 10 percent low at 5, and flat from 7 up. So 7 is the
smallest converged count, and the preset ladder follows it with one step of
conservatism:

- `reference` 15 — two steps above the converged count. CONFIRMED.
- `standard` 9 — one step above the converged count. CONFIRMED.
- `rapid` 5 — one step BELOW it, and the stated compromise: the slant index runs
  about 10 percent low, and the mean power holds inside 0.11 dB. It stays above
  the moment floor of 4. CONFIRMED as a deliberate trade.

**Measured 2: the moment error of the grouping.** The `sampling.py` self-check
compares the production grouping with `schmidt.turbulence.moment_error` on the
default site profile at the zenith, 200 layers, `standard` preset, 9 screens.
Every moment `0 <= m <= 7` holds inside 0.15 percent, against 0.04 percent for
one screen for each layer. So the `Cn2`-weighted centroid grouping satisfies
Eq. (9.65) in practice without solving it. CAVEAT: this metric is not sensitive
to the near-ground layers, because they all sit at almost the same distance from
the source.

**Measured 3: the three turbulent examples, re-run.** The screen counts and the
agreement numbers moved, and the doc statements were updated with them.

| Script | Preset | Screens, before | Screens, after | Agreement, before | Agreement, after |
|---|---|---|---|---|---|
| `turbulent_terrestrial.py` | `standard` | 9 | 9 | fidelity-0 SMF Term about 2.5 dB above the field | 2.3 dB (4.61 dB against 2.32 dB). The horizontal planner takes no layer list, so WP7 did not move it. |
| `turbulent_downlink.py` | `rapid` | 20 | 5 | index ratios 1.10 to 1.27; FAST 0.7 to 2.9 dB above the field | index ratios 1.01, 1.19, 1.28; FAST 2.7 to 3.9 dB above the field, and 1.8 to 3.0 dB on the turbulence part alone |
| `turbulent_uplink_reciprocity.py` | `rapid` | 20 | 5 | means 0.18 dB (90 deg), 0.54 dB (30 deg) | 0.19 dB (90 deg), 1.05 dB (30 deg) |

The two space scripts now run in about one minute each, against several minutes
before. No sampling warning appears in any of the three runs.

**One effect to watch.** The downlink fibre-coupling loss moved by about 1.2 dB
at 30 degrees when the plan went from 20 screens to 5. The sweep shows that a
properly grouped plan does NOT move that number with the count (the `smf` column
of case B scatters about 1.1 dB with no trend, at 12.6 to 13.8 dB). So the shift comes from the LAYERING,
not from the count: the old one-screen-per-layer plan put screens 2 m and 3 m
from the pupil, and the new bottom group sits at the `Cn2`-weighted centroid of
the ground layers, about 150 m away. Which layering is right for a near-pupil
`Cn2` spike is an open question. It is a candidate for the Gauss-quadrature
route of decision 3.
