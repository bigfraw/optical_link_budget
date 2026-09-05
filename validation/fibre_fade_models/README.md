# Analytic models of the fibre-coupled power fade (a plan, backlog 1-9)

STATUS (2026-09-05): a PLAN and a design note. No script exists yet. The owner
wants to review it before any code is written. This folder is the home of the
work when it starts.

This note is for an agent that picks up the terrestrial fidelity-1 rung
(backlog 1-9 in `docs/backlog.md`). It records what the repo already shows
about the fibre-coupled fade, the candidate analytic models, the order to try
them, and the temporal extension the owner likes. Read `CLAUDE.md`, the 1-9
entry of `docs/backlog.md`, and the two studies below before you start.

It is VALIDATION ONLY: nothing here changes an `olb` module until step 4.

## 1. The question that started it

Is the point scintillation index close to the scintillation index of the
uncorrected fibre-coupled power? The stored results say NO. The two fades are
far apart, and they do not move together.

## 2. What the repo already holds

Three fidelity-2 studies store per-case statistics of the hero downlink
(1550 nm, a 700 mm ground terminal with an uncorrected SMF, a 100 mm space
terminal at 500 km, `cn2_ground = 1.7e-14`, wind 21 m/s, 1000 trials for
each case):

- `validation/tail_convergence/tail_convergence_el30_results.json`
- `validation/outer_scale_tail/outer_scale_tail_el30_results.json`
- `validation/outer_scale_tail/outer_scale_tail_el20_results.json`

Each case gives `sigma2_I_point` (the centre-pixel irradiance index),
`sigma2_P` (the 700 mm BUCKET index, of order 0.005, a footnote), the SMF
loss quantiles (`quantiles_db`) and the point fade quantiles
(`point_quantiles_db`). No stored JSON gives the fibre-coupled power INDEX
itself. The raw `Campaign` directories are git-ignored and live on the
owner's desktop; each trial there holds `smf_eta`, `collected_power`, and the
receive-field patch, so the direct index is a few lines of post-hoc analysis
with no rerun.

The comparison, with the spread as p5 minus p50 and p1 minus p50:

| case | point index | point p5 | SMF p5 | point p1 | SMF p1 |
| --- | --- | --- | --- | --- | --- |
| 30 deg, 15 screens, 1024 px, L0 = inf | 0.23 | 3.5 dB | 12.7 dB | 5.4 dB | 20.3 dB |
| 30 deg, 15 screens, 1024 px, L0 = 25 m | 0.23 | 3.6 dB | 11.3 dB | 5.2 dB | 18.2 dB |
| 30 deg, 25 screens, 1024 px, L0 = inf | 0.23 | 3.8 dB | 11.7 dB | 5.5 dB | 20.2 dB |
| 30 deg, rapid (5 screens, 256 px) | 0.21 | 3.5 dB | 12.8 dB | 5.4 dB | 21.0 dB |
| 20 deg, 15 screens, 1024 px, L0 = inf | 0.49 | 5.4 dB | 13.7 dB | 8.0 dB | 20.2 dB |
| 20 deg, 15 screens, 1024 px, L0 = 25 m | 0.48 | 5.4 dB | 11.8 dB | 7.9 dB | 18.3 dB |

Three readings:

1. The point fade agrees with its own lognormal index. A lognormal with
   sigma2_I = 0.23 has a p5 spread of 3.6 dB (Andrews and Phillips, 2nd ed.,
   DOI 10.1117/3.626196, Ch. 5, Eq. (95), printed 157, and Ch. 11, Eq. (24),
   printed 451). The data reads 3.5 dB.
2. The SMF fade is three to four times wider, and it does NOT follow the point
   index. From 30 deg to 20 deg the point index doubles, but the SMF spread
   stays at 11 to 14 dB.
3. The SMF spread sits at the SPECKLE limit. A negative-exponential
   distribution (index 1.0, the many-cell limit of a coherent sum) has a p5
   spread of 11.3 dB and a p1 spread of 18.4 dB. The stored SMF spreads are
   11.3 to 13.7 dB at p5 and 18.2 to 21.0 dB at p1. The excess above the
   speckle limit is the tilt, which is what the outer scale moves (2-P5).

So at D/r0 of about 5.5 the uncorrected fibre-coupled power is PHASE
dominated, and the point index tells you nothing about it. A fibre fade model
must model the phase.

## 3. The candidate analytic models

The owner's decision: run the fidelity-2 sweep FIRST, then fit every
candidate post hoc against the stored trials. None of the four below is a
derived closed form of the fibre-coupled power. Each is a candidate family
with a proposed parameter map, and the sweep judges it.

### 3.1 Lognormal times Rician (the candidate Term)

Split the coupled power into two factors: the bucket power times the coupling
efficiency, `P_smf = P_bucket * eta`.

- `P_bucket` is the aperture-averaged lognormal that 1-6 certified on one
  weak path (`validation/lognormal_certification/`).
- `eta` is the squared modulus of the overlap of the pupil field with the
  fibre mode. Write the pupil field as a coherent part plus a random part.
  The coherent amplitude is the Strehl amplitude `exp(-sigma_phi^2 / 2)`;
  the random part is a complex Gaussian that holds the rest of the power.
  Then `eta` is Rician with `K = exp(-sigma_phi^2) / (1 - exp(-sigma_phi^2))`.
- `sigma_phi^2` is the Noll residual: `1.03 (D/r0)^(5/3)` untracked, or
  `0.134 (D/r0)^(5/3)` after the tip-tilt walk-off Term takes the tilt
  (Noll 1976, DOI 10.1364/JOSA.66.000207, Table IV).

The product is the lognormal-Rician FAMILY that `olb/turbulence/andrews/
distributions.py` already holds (Ch. 9, Sec. 9.9.2, Eq. (133), printed 369).
BE CLEAR about what that is: the book gives it as a heuristic PDF for the
irradiance at a POINT receiver (Churnside and Clifford,
DOI 10.1364/JOSAA.4.001923). It derives nothing about a fibre. The fibre
overlap has the same algebraic shape, so the family is reusable, but the
parameter map above is olb's, not the book's. The olb copy is the PDF only
(no CDF, no quantile, no sampler), so the Term needs a sampler.

Three approximations sit inside it, and the sweep must test each:

- The Gaussian random part holds only at large D/r0 (many coherence cells).
  At small D/r0 a few low-order modes carry the phase, and the overlap
  fluctuation is not Gaussian. A terrestrial link often sits in that
  transition.
- The bucket power and `eta` are taken as independent. On a strong path the
  amplitude and the phase correlate, so this is a weak-regime model. Flag it
  at `LOGNORMAL_PDF_LIMIT = 0.25`, the house rule.
- The Strehl amplitude is the extended Marechal form, valid below about
  `MARECHAL_SIGMA2_MAX = 1.0` rad^2 (T. S. Ross, DOI 10.1364/AO.48.001812).

A free check: at D/r0 = 5.5 the K factor is near zero, the Rician becomes
negative-exponential, and that is exactly the 11 to 14 dB p5 spread of
Section 2.

### 3.2 The Zernike Monte Carlo (the cheap reference)

The von Karman phase spectrum filtered by the pupil IS the Noll Zernike
covariance matrix, mode by mode (Noll 1976, DOI 10.1364/JOSA.66.000207). So:
draw a Gaussian Zernike vector from that covariance, build the pupil phase,
take the numeric overlap with the fibre mode (`olb/waveoptics/smf.py`
`coupling_efficiency` on a small vacuum grid), and multiply by ONE lognormal
amplitude scalar. That is a Monte Carlo with NO propagation, so thousands of
trials run in seconds. It is what FAST does in the far field, but olb owns
the near-field Gaussian weighting, which FAST lacks (it is far-field only).
Use it to calibrate 3.1 and to find the D/r0 where the Rician K map breaks.
Take `r0` from `gaussian_fried_parameter_profile`, not the plane-wave r0,
because the Gaussian illumination changes the effective D/r0.

### 3.3 The coupled-power variance from the coherence function

The MEAN coupling is the Dikmelik and Davidson integral of the fibre mode
against the mutual coherence function `exp(-3.44 (r/r0)^(5/3))`
(DOI 10.1364/AO.44.004946). The second moment needs the FOURTH-order
coherence; under the Gaussian-phase assumption it follows from the structure
function at four points. That gives an analytic index of the coupled power
with no draw, and then a family selector can read it, the way the downlink
`model="auto"` does. It is the most elegant route and the most derivation
work. Not a first step.

### 3.4 The receiver-kind split

The bucket and the fibre are DIFFERENT random variables (the 1-9 entry says
why). Fit them separately: the bucket against the lognormal, the gamma-gamma
and the K distribution; the fibre against the same three plus 3.1. Split the
fibre fade into its tilt part and its higher-order part (`spot_metrics` gives
the centroid), because the tilt is what a tracked terminal removes (2-AO) and
what the outer scale sets (2-P5).

## 4. The temporal extension (the owner likes this one)

The fade DEPTH work above has a temporal twin. It matters because a budget
needs the fade rate and the fade duration as well as the depth (the
`andrews/temporal.py` pieces exist and no Term reads them).

A first guess, "convolve the bucket temporal spectrum with the tilt
spectrum", is WRONG, and the owner said why: the coupling is a nonlinear
function of ALL the modes. Every Zernike mode moves the overlap, with a
weight set by the fibre mode: tilt walks the spot off the core, defocus grows
it, then astigmatism and coma, and so on. So the spectrum of `eta` is the
spectrum of a nonlinear functional of the whole phase vector, not the tilt
spectrum.

Each mode has its OWN temporal spectrum. Conan, Rousset and Madec (1995,
DOI 10.1364/JOSAA.12.001559) give the temporal power spectrum of each Zernike
coefficient under frozen flow: the tilt falls as `f^(-2/3)` at low frequency;
every higher order is FLAT at low frequency with a knee near
`0.3 (n + 1) v / D` that rises with the radial order n; all of them roll off
as `f^(-17/3)` past the knee. So the higher orders put their power at HIGHER
frequencies than the tilt, and a tilt-only model under-counts the fast fades.

The one part that survives: the product of two independent processes gives
the convolution of their spectra. `P_bucket * eta` is that product, so the
bucket spectrum (Andrews Ch. 8) convolved with the TRUE `eta` spectrum is the
right assembly. The gap is the `eta` spectrum itself. Two ways to get it:

- **Fidelity 2, the truth.** The frozen-flow axis is the `temporal.py` stub
  in `olb/waveoptics/turbulence/`. Translate each screen at its layer wind
  (the Bufton profile in `olb/turbulence/profiles.py`), propagate every time
  step, and store the per-step coupled power. That gives the `eta` spectrum,
  the fade rate and the fade duration with no model. It is the same
  `Campaign` machinery with a time index inside a trial, so the post-hoc
  analysis is the same as the depth study.
- **Fidelity 1, the Zernike temporal Monte Carlo.** Draw a TIME SERIES of
  Zernike coefficients from the Conan spectra, push each sample through the
  numeric overlap, and read the `eta` spectrum from the result. It is the
  temporal twin of 3.2: exact in the phase-only limit, seconds to run, and
  the natural fit target once a fidelity-2 series exists.

## 5. The order of work

Do NOT write a script before the owner reviews this note.

1. The fidelity-2 sweep of 1-9 (the desktop, 8 to 12 workers, about 20 min
   for each 1000-trial campaign at 512 px; see
   `validation/campaign_resources/README.md`): a stronger Cn2, a longer path,
   a focused launch, at the operating `L0 = 25 m`, with the receive-field
   patch stored so every detector is a post-hoc `recouple`.
2. The post-hoc fit script: for each campaign, form `P_bucket` (`recollect`)
   and `P_smf` (`recouple`) trial by trial, and fit every family of Section 3
   (a moment fit and a maximum-likelihood fit). Report the p10, p5, p1 fade
   against the empirical quantiles with the bootstrap bars the existing
   studies use, the direct fibre-coupled index, and the tilt/higher-order
   split.
3. The Zernike Monte Carlo of 3.2 as the cheap cross-check at every D/r0.
4. Decide the family for each receiver kind, and wire
   `terrestrial_budget(fidelity=1)` through `olb/models/fade.py` (step 4 of
   the 1-9 entry). Add a sampler for the lognormal-Rician if 3.1 wins.
5. The temporal extension of Section 4, once the depth rung is wired.

## 6. Conventions that apply

- ASD-STE100 in every docstring, comment and commit message
  (`CONVENTIONS.md`).
- Every equation cites a DOI.
- Loss is positive dB.
- A study keeps its scripts, results JSON and logs at its top level, and its
  figures in `figures/`. Run it as `python -m validation.fibre_fade_models.<script>`.
- A campaign directory is git-ignored; add `validation/fibre_fade_models/campaigns/`
  to `.gitignore` when the first script lands.
