# Aperture-averaged lognormal certification (backlog 1-6)

Does the cheap analytic weak-turbulence route give a TRUSTWORTHY received-power
distribution?

The analytic route (`olb.links.terrestrial.terrestrial_scintillation_term`, and
the same closed form in `olb.links.downlink._lognormal_term`) takes the
aperture-averaged index `sigma2_P = A sigma2_I` and it draws the power from a
lognormal. But an aperture integrates a CORRELATED lognormal field, and a sum of
lognormals is not a lognormal: as the diameter D grows, the power moves to a
Gaussian and the fade tail gets thin. A finite Gaussian beam adds beam wander,
which puts a POINTING tail on the distribution that one index cannot describe.

The script measures that against the fidelity-2 split-step Monte Carlo
(`olb.waveoptics.turbulence`), which solves the field and makes no distribution
assumption.

## Method

One fixed horizontal path (2 km, `Cn2 = 3e-15`, 1550 nm) that stays FIRMLY weak
(`sigma_R^2 = 0.21`, the `rytov_weak` "weak" tier). Two launches (collimated and
diverged) x four receive diameters (1, 5, 15, 40 cm), so `D/rho_0` runs from 0.20
(point-like) to 7.89 (strong averaging). `rho_0` is the plane-wave coherence
radius of Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196, Ch. 6.

Each case reports FIVE separate things, so an INDEX error, a FILTER error and a
SHAPE error do not mix:

1. the analytic `sigma2_P` against the measured `var(P)/mean(P)^2`;
2. **the index split (added 2026-09-01).** `sigma2_P = A sigma2_I` holds TWO
   analytic models, so a `sigma2_P` comparison alone cannot say which one is
   wrong. The script measures the POINT index too, and it prints four columns:
   the analytic `sigma2_I` (the Dios on-axis Gaussian beam-wave form) against the
   measured point `sigma2_I`, and the analytic `A` (the Churnside weak
   plane-wave fit) against `A_eff = sigma2_P_sim / sigma2_I_sim`;
3. **the beam-filling check (added 2026-09-01).** The vacuum received beam
   radius `w(L)` and the captured power fraction
   `eta_fill = 1 - exp(-2 (D/2)^2 / w(L)^2)`. A case past `eta_fill = 0.5` is
   FLAGGED as BEAM-FILLING-LIMITED (see below);
4. **the absolute impact (added 2026-09-01).** The measured fade spread (the
   standard deviation of the loss in dB) beside the 5 % fade depth, so the table
   shows both "how wrong" and "how much it matters";
5. the fade quantiles (the loss exceeded 10 %, 5 % and 1 % of the time) of the
   sim against the analytic lognormal, AND against a lognormal REFIT to the
   MEASURED index. The refit leg is the pure SHAPE test; plus the skew of
   `ln P`. A lognormal power gives a Gaussian `ln P`, so the skew is 0. A
   negative skew is the drift to a Gaussian POWER (a thin fade tail).

Fidelity 2 models NO tip-tilt correction (backlog 2-AO), so the sim holds the
FULL beam wander. That is part of what the script tests.

### One propagation for the whole aperture sweep

The script runs `propagate_turbulent_field` (the public single-snapshot entry
point of `olb.waveoptics.turbulence.run`), which gives back the complex
receive-plane field of one trial. It then clips that ONE field at the point
estimator AND at every receive diameter. So:

- every aperture reads the SAME atmosphere, and
- the cost does not grow with the number of diameters (the quick run fell from
  about 9 minutes to about 2.5 minutes).

The grid is the grid of the LARGEST diameter, because the sizer widens the grid
with the receive aperture. The script does NOT extend `TurbTrial` or
`TurbWaveResult` (backlog 2-I1). It cross-checks itself against the old
one-run-for-each-aperture path: three matched seeds on the 40 cm aperture, where
the two paths share a grid, agree BIT FOR BIT (`max relative difference 0.0`).
The small apertures move by under 1 % against the earlier reading, because they
now sit on the 40 cm grid and not on their own.

### The point-index estimator

The point index is the mean irradiance inside a small ON-AXIS disc of diameter
8 mm, taken with the same clip that the apertures use. A single centre PIXEL is
the literal point value, but the pixel is only 0.93 mm (collimated) to 1.82 mm
(diverged) on these grids, so a one-pixel estimate rides on grid-scale noise.
The 8 mm disc is 0.14 of the Fresnel scale `sqrt(lambda L) = 5.6 cm`, and the
analytic filter says it averages the index down by 1.9 % only.

CAVEAT: the grid holds no irradiance structure below one pixel, so the estimator
is a LOWER BOUND on a true point index. The Fresnel scale is 30 to 60 pixels
across here, so that bias is small. The consistency check is the 1 cm aperture,
where the measured `A_eff` (0.986) sits within 2 % of the analytic `A` (0.968).

### The beam-filling flag

The Churnside factor `A` is a PLANE-WAVE fit: it assumes the aperture reads a
piece of a much wider field. A receive aperture that holds most of the beam
breaks that assumption (backlog 2-N2), and then its `A` comparison tests the
aperture-holds-the-beam regime, not the filter. The flag travels into the log,
the results JSON (`beam_filling_limited`) and the figure (a grey panel).

## Run

    python -m validation.lognormal_certification.lognormal_certification
    python -m validation.lognormal_certification.lognormal_certification --full

QUICK (the default) is 150 trials for each LAUNCH on the `rapid` preset, about
2.5 minutes on 16 threads (one propagation set serves all four diameters). The
1 % fade is UNDER-SAMPLED at that trial count, and
the script prints a warning that says so; read the 10 % and the 5 % rows only.
`--full` (1500 trials, the `standard` preset) is the run for the deep tail, and
it must run before a certification is recorded.

Outputs, next to the script: `lognormal_certification.log` and
`lognormal_certification_results.json` (it holds the per-trial loss of each
case). The figure `figures/lognormal_certification.png` shows the measured
histogram against the analytic and the refit lognormal, with the fade lines.

## Verdict

A PASS certifies that the cheap analytic weak-turbulence calculation gives a
trustworthy power distribution. A FAIL BOUNDS where it may be used, and it points
to a composite (lognormal x pointing) model or to a direct empirical sampler.
A disagreement past 0.5 dB at the 5 % fade is NOTABLE.

### First reading, QUICK mode (2026-09-01, 150 trials, `rapid` preset)

PASS, but READ THE LEGS SEPARATELY. The worst 5 % disagreement of the whole
analytic route is 0.30 dB, and of the SHAPE leg alone 0.13 dB.

| case | D [cm] | `D/rho_0` | `s2I` an | `s2I` sim | ratio | `A` an | `A_eff` | ratio | fill | flag |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| collimated | 1 | 0.20 | 0.0744 | 0.0615 | 0.83 | 0.9681 | 0.9864 | 1.02 | 0.00 | |
| collimated | 5 | 0.99 | 0.0744 | 0.0615 | 0.83 | 0.4149 | 0.5817 | 1.40 | 0.03 | |
| collimated | 15 | 2.96 | 0.0744 | 0.0615 | 0.83 | 0.0518 | 0.1335 | 2.58 | 0.25 | |
| collimated | 40 | 7.89 | 0.0744 | 0.0615 | 0.83 | 0.0055 | 0.0029 | 0.53 | 0.87 | BEAM-FILLING-LIMITED |
| diverged | 1 | 0.20 | 0.0744 | 0.0676 | 0.91 | 0.9681 | 0.9864 | 1.02 | 0.00 | |
| diverged | 5 | 0.99 | 0.0744 | 0.0676 | 0.91 | 0.4149 | 0.5713 | 1.38 | 0.01 | |
| diverged | 15 | 2.96 | 0.0744 | 0.0676 | 0.91 | 0.0518 | 0.1514 | 2.92 | 0.07 | |
| diverged | 40 | 7.89 | 0.0744 | 0.0676 | 0.91 | 0.0055 | 0.0140 | 2.54 | 0.39 | |

| case | D [cm] | `s2P` an | `s2P` sim | ratio | sd loss [dB] | p5 sim [dB] | p5 an [dB] | delta [dB] |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| collimated | 1 | 0.0720 | 0.0607 | 0.84 | 1.071 | 1.736 | 2.035 | -0.299 |
| collimated | 5 | 0.0309 | 0.0358 | 1.16 | 0.816 | 1.342 | 1.312 | +0.030 |
| collimated | 15 | 0.0039 | 0.0082 | 2.13 | 0.396 | 0.684 | 0.451 | +0.233 |
| collimated | 40 | 0.0004 | 0.0002 | 0.43 | 0.058 | 0.118 | 0.146 | -0.027 |
| diverged | 1 | 0.0720 | 0.0667 | 0.93 | 1.074 | 1.882 | 2.035 | -0.153 |
| diverged | 5 | 0.0309 | 0.0386 | 1.25 | 0.849 | 1.453 | 1.312 | +0.141 |
| diverged | 15 | 0.0039 | 0.0102 | 2.65 | 0.439 | 0.709 | 0.451 | +0.257 |
| diverged | 40 | 0.0004 | 0.0009 | 2.31 | 0.134 | 0.237 | 0.146 | +0.091 |

- **The lognormal FAMILY holds.** With the index refit to the measured value,
  every case agrees inside 0.13 dB at the 5 % fade and inside 0.31 dB at the 1 %
  fade, at every `D/rho_0` up to 7.9. The skew of `ln P` stays between -0.38 and
  +0.19, with no trend against D. So the drift to a Gaussian power is NOT visible
  in this weak band, and the beam wander of the sim does not build a pointing
  tail that the lognormal cannot hold.
- **The FILTER `A` is the fault, not the POINT index.** The split says it
  clearly. The analytic point `sigma2_I` reads 1.10 (diverged) to 1.21
  (collimated) times HIGH, an error of 10 to 20 %. The Churnside filter
  OVER-AVERAGES much more: `A_eff/A` is 1.4 at `D/rho_0 = 1`, 2.6 to 2.9 at
  `D/rho_0 = 3`, and it stays at 2.5 out at `D/rho_0 = 7.9` for the launch that
  does not fill the aperture. So the `sigma2_P` error near `D/rho_0 = 3` is
  almost all filter, and the point index partly HIDES it (a high point index and
  an over-averaging filter pull in opposite directions).
- **The D = 40 cm collimated column is BEAM-FILLING-LIMITED. Do not read it as a
  filter error.** The collimated received beam radius is `w(L) = 19.7 cm`, so
  the 40 cm aperture catches `eta_fill = 0.87` of the beam. It is close to a
  total-power measurement, and total power fluctuates little, so `A_eff` falls
  BELOW the analytic `A` (ratio 0.53) and the whole route reads 2.3 times HIGH.
  The diverged launch at the same diameter has `w(L) = 40.4 cm` and
  `eta_fill = 0.39`, and it behaves like every other unfilled case (`A` ratio
  2.54). The reversal is the fill fraction, not the diameter. This is
  backlog 2-N2 measured.
- The diverged launch reads a bigger point index than the collimated launch
  (0.0676 against 0.0615), but the analytic Term gives BOTH the same number: the
  on-axis Gaussian beam-wave index takes the waist only, and no divergence. That
  is a known blind spot of the analytic leg.
- **The absolute impact is small in this band.** The fade spread falls from
  1.07 dB at D = 1 cm to 0.06 to 0.13 dB at D = 40 cm. The WORST relative index
  error (2.9 times, diverged, D = 15 cm) moves the 5 % fade by 0.26 dB only,
  because the averaging has already taken the spread to 0.44 dB. So a large
  relative index error does NOT make a large budget error here.

The quick-mode numbers above are kept for history; the certification of record
is the FULL run below.

### Certification of record, FULL run (2026-09-01, 1500 trials, `standard` preset)

PASS. Worst 5 % disagreement of the whole analytic route 0.289 dB, of the SHAPE
leg alone 0.128 dB; at the 1 % fade 0.413 dB and 0.210 dB. Both worst cases are
the diverged 15 cm receiver. See `lognormal_certification.log` for the full
tables. What the full run changes against the quick reading:

- **The point-index bias mostly disappears.** The analytic `sigma2_I` reads only
  3 to 4 % high (`sim/an` 0.96 collimated, 0.97 diverged). The 10 to 20 % bias of
  the quick run was a `rapid`-preset artifact, not a model error. The diverged
  launch still reads a slightly larger index than the collimated one, and the
  analytic Term still gives both the same number (the waist-only blind spot),
  but at this path the difference is 1 %.
- **The filter fault stands, milder.** `A_eff/A` is about 1.2 at `D/rho_0 = 1`
  and 1.8 (collimated) to 2.5 (diverged) at `D/rho_0 = 3`, still 2.4 at 7.9 for
  the unfilled launch. The Churnside `A` remains the weak leg of the route.
- **The shape holds to the deep tail.** At 1500 trials the 1 % fade is sampled
  (about 15 trials past it); the refit lognormal agrees inside 0.21 dB there in
  every case, and the skew of `ln P` sits near -0.2 with no trend against D.
- **The beam-filling flag confirms.** Collimated 40 cm: `eta_fill = 0.87`,
  `A` ratio 0.40, flagged. Diverged 40 cm (`eta_fill = 0.39`) behaves like the
  unfilled cases.

So the DISTRIBUTION SHAPE is CERTIFIED for this weak band, to the 1 % fade, and
the open question is NAMED: the weak aperture-averaging factor `A` over
`1 <= D/rho_0 <= 8`, below the beam-filling limit. Still untested (the 1-8
gates): a focused launch, a stronger `Cn2`, and a longer path.
