# G-tilt tilt-sensing against the wrapped slopes

This study asks one question: for a terrestrial fidelity-2 TipTilt correction,
does the GRADIENT (centroid) tilt sense the tilt better than the
wrapped-gradient slopes?

## The two routes

- **slopes** — the shipped route. It differences the WRAPPED phase between two
  adjacent pixels and fits the modes by least squares
  (`olb.waveoptics.compensation.slopes`). It aliases where the local phase step
  passes pi.
- **gtilt** — the far-field intensity centroid, which is the intensity-weighted
  mean phase gradient (`olb.waveoptics.compensation.gtilt`). It reads the
  COMPLEX field, so it never unwraps a phase and it does not alias a local step.
  It is tilt-only. Source: G. A. Tyler, J. Opt. Soc. Am. A 11, 358 (1994),
  DOI 10.1364/JOSAA.11.000358.

## The method

`gtilt_vs_slopes.py` reopens each stored terrestrial 2-TC campaign
(`validation/terrestrial_campaigns/`) and runs three POST-HOC coupling passes
over the 10 cm single-mode fibre, with NO new propagation: uncorrected, a
TipTilt correction sensed by the slopes, and the same correction sensed by
G-tilt. The judge is the deep-fade (p5) coupling loss.

    python -m validation.gtilt_sensing.gtilt_vs_slopes --workers 8 --figures

## The result: G-tilt never wins

The p5 coupling-loss GAIN of G-tilt over the slopes (a positive number means
G-tilt fades less), 2000 trials for each cell, a 10 cm SMF receiver:

| sigma_R^2 | grid | p5 gain of gtilt over slopes |
|-----------|------|------------------------------|
| 0.21      | 1024 | -0.04 dB                     |
| 0.71      | 1024 / 2048 | -0.08 / -0.12 dB      |
| 1.14      | 1024 / 2048 | -0.18 dB              |
| 4.07      | 1024 / 2048 | -0.36 / -0.44 dB      |
| 13.56     | 1024 / 2048 | -3.65 / -0.50 dB      |

G-tilt is worse at every strength, and it gets worse as the turbulence gets
stronger. It is WORST on the coarsest grid (1024 px, the most aliasing) — the
opposite of what an anti-aliasing sensor would do.

## Why the slopes win

1. **The slope fit is robust to aliasing for the LOW-order tilt.** It is a
   least-squares fit over every in-aperture pixel pair, so the aliased pairs
   average out as noise. The aliasing damage that item 2-AO / V5 saw was to the
   HIGHER-order `AO(21)` reconstruction, which G-tilt cannot address.
2. **The far-field centroid is speckle-noisy in scintillation.** At
   sigma_R^2 >~ 4 the received field is broken into speckle, so its centroid
   (still the exact mean gradient) is not the tilt that best couples a fibre.

## A second route that also lost: the spectral phase gradient

A spectral gradient of the complex field, `Im(conj(E) dE)/|E|^2` with
`dE = IFFT(i k FFT(E))`, is the exact gradient for a band-limited PERIODIC
field. On a real turbulent field it fails: the field on a finite grid is not
periodic, so the FFT derivative rings and reads the gradient with about 40 %
RMS error, worse than the local wrapped-slope even with no aliasing. It has no
winning regime (the local slope is exact below pi, and above pi the field is
under-sampled and non-periodic), so it is not kept.

## The lesson

For modal sensing on a finite, non-periodic, cropped turbulent field, a LOCAL
least-squares operator (the wrapped slopes) beats a GLOBAL FFT operator (the
centroid, the spectral gradient). `gtilt` stays in the code as an opt-in source
only, because it is a faithful model of a real focal-plane centroid tracker.
