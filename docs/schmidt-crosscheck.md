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
| sqrt(lambda z) (172) | Fresnel length, the scintillation scale. No olb name. The olb pitch rule uses `pixels_per_r0` only. |
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

# Table 2 — gaps and suggestions

| gap id | book section | book eq | capability | target module | priority |
| --- | --- | --- | --- | --- | --- |

# Table 3 — constants ledger

Each row holds an olb constant that carries NO citation today. Fill the book
columns as the cross-check proceeds. An empty book column means that we have
not yet found a source, or that the book gives none.

| olb constant | olb value | location | book quantity | book eq | printed p | pdf p | status |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `QualityPreset.min_screens` | 15 / 9 / 5 | `olb/waveoptics/turbulence/sampling.py:105` | | | | | flagged |
| `QualityPreset.sigma2_r_screen_max` | 0.05 / 0.10 / 0.25 | `olb/waveoptics/turbulence/sampling.py:105` | | | | | flagged |
| `QualityPreset.boundary_width_frac` | 0.125 / 0.125 / 0.10 | `olb/waveoptics/turbulence/sampling.py:105` | | | | | flagged |
| `super_gaussian_boundary` `power` | 8 | `olb/waveoptics/turbulence/splitstep.py:29` | | | | | flagged |
| `super_gaussian_boundary` `width_frac` | 0.125 | `olb/waveoptics/turbulence/splitstep.py:29` | | | | | flagged |
| `QualityPreset.guard` | 4 / 3 / 2 | `olb/waveoptics/turbulence/sampling.py:105` | | | | | flagged |
| `QualityPreset.pixels_per_r0` | 4 / 3 / 2 | `olb/waveoptics/turbulence/sampling.py:105` | | | | | flagged |
| `QualityPreset.fresnel_weight_min` | 0.005 / 0.02 / 0.05 | `olb/waveoptics/turbulence/sampling.py:105` | | | | | flagged |
| `GridSpec.for_scenario` `guard` | 4.0 | `olb/waveoptics/grid.py:98` | | | | | flagged |
| `GridSpec.for_scenario` `pixels_per_feature` | 16 | `olb/waveoptics/grid.py:98` | | | | | flagged |
| `N_MIN` | 256 | `olb/waveoptics/grid.py:36` | | | | | flagged |
| `forvard_max_z` | z_max = N dx^2 / lambda | `olb/waveoptics/grid.py:209` | | | | | flagged |

---

# Work-package notes

## WP1 — Fourier foundations (Chs. 2, 3)

## WP2 — Fresnel propagators (Ch. 6)

## WP3 — Sampling constraints (Chs. 7, 8)

## WP4 — Turbulence and screens (Ch. 9)

## WP5 — Examples

## WP6 — Retrofit and documentation

## WP7 — The `min_screens` revision
