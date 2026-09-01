# RESOLVED: fidelity-1 vs fidelity-2 uplink scintillation index

Date: 2026-08-28. Status: RESOLVED AND CONFIRMED (Section 6). This note
replaces the 2026-08-28 unresolved handoff.

## 1. Verdict

- The two indices AGREE for an UNDERFILLED launch (beam wander much smaller
  than the far-field spot): fidelity 1 sits 1.2x to 1.3x above fidelity 2,
  inside the joint error bars.
- For a FILLED launch (wander comparable to the far-field spot,
  beta_rms / w_L near 1), **fidelity 1 OVER-PREDICTS sigma2_I by 2x to 7x**.
- The cause is the sampled OFF-AXIS term of the Dios model
  (DOI 10.1364/AO.43.003866, Eq. (20)): an unsaturated weak-fluctuation
  (Rytov) expression, evaluated at the sampled wander displacement. When the
  wander reaches the beam edge, that expression is far past its validity and
  it has no saturation. The paper's own validation figure (Fig. 5) shows the
  same departure: its FFT-BPM reference falls BELOW the semianalytic curve at
  large W0 ("the appearance of the saturation effect").
- Fidelity 2 is the trustworthy leg in the filled regime. Its result is
  converged against the grid, the screen count, and the subharmonic
  (low-frequency tilt) content. The mode-matched reciprocity comparison is
  valid: the earlier "different quantities" excuse stays wrong, and the
  earlier suspicion that the fidelity-2 grid truncates the wander is now
  REJECTED by measurement (Section 4).

## 2. The experiment

`validation/uplink_sigma2i/uplink_farfield_reciprocity.py`. Three upgrades over the earlier
(void) comparison:

1. **Mode-matched.** The launch aperture is 1.0 m, so the clip on the
   transmit Gaussian is < 2e-7 of the power. Fidelity 1 (pure Gaussian of
   waist w0) and fidelity 2 (`_ground_transmit_mode`) see the same mode. The
   mode is built from `Transmitter.waist_m`, never from the aperture.
2. **The far-field MAP, not one point.** A satellite offset x is exactly a
   tilt exp(-i k x.r / L) on the downlink plane wave (Shapiro reciprocity,
   DOI 10.1364/JOSA.61.000492). So one zoom FFT (chirp-z) of
   E_down * conj(psi_tx) gives the uplink flux at EVERY satellite offset:
   the full instantaneous spot. The f = 0 bin equals the runner's eta_turb
   bit for bit (asserted). Per snapshot the script reads the on-axis flux,
   the instantaneous beam-centre offset (the wander beta), and the flux at
   that centre (the tracked flux). So each INGREDIENT of the Dios model is
   measured on its own.
3. **Realistic turbulence.** The default HV57 site profile (and a 0.3x scale
   that keeps fidelity 1 weak-valid for the filled beam), not the near-vacuum
   Cn2 = 1e-15 of the void comparison.

## 3. The numbers (600 km, 60 deg, 1550 nm, DEFAULT_HS, HV57 site profile)

250 fidelity-2 snapshots per case; 40000 fidelity-1 samples. Full data:
`uplink_farfield_reciprocity_results.json` (key case + variants; the JSON of
the first full run was lost to a restart, the log
`uplink_farfield_reciprocity_run.log` holds every headline number).

| case | fid-1 sigma2_I | fid-2 sigma2_I | ratio | fid-1 weak_valid |
|------|---------------|----------------|-------|------------------|
| w0=0.06, HV x0.3 | 0.052 | 0.042 | 1.2 | True |
| w0=0.06, HV x1.0 | 0.342 | 0.255 | 1.3 | True |
| w0=0.18, HV x0.3 | 0.938 | 0.412 | 2.3 | True (sigma2_x = 0.167) |
| w0=0.18, HV x1.0 | 9.017 | 1.245 | 7.2 | **False** (sigma2_x = 0.639) |

Decomposition of the clean filled case (w0 = 0.18, HV x0.3): fidelity 1's
own wander-only contribution is 0.085 and its on-axis index is 0.004, so
about 0.85 of its 0.938 comes from the lognormal driven by
sigma2_off(beta). The measured beam-frame (tracked) index of fidelity 2 at
the same point is **0.018**: the real fluctuation on the displaced axis is
small, and the off-axis Rytov term inflates it by an order of magnitude.

## 4. Ingredient findings

- **Short-term waist.** The Dios w_st matches the measured beam-frame width
  to a few percent in every weak-valid case (5.81 vs 5.92; 2.66 vs 2.40;
  6.68 vs 6.52 m). It drifts only in the flagged-invalid case. The beam
  spread machinery is NOT the problem.
- **Beam wander.** The measured satellite-plane wander variance is
  1.8x to 2.0x the Dios/Belmonte 2.07 form and 0.55x the Andrews 7.25 form,
  CONSISTENTLY across all four cases (ratios 2.02, 1.93, 2.00, 1.81). It is
  CONVERGED: the wide_x2 variant (doubled grid side at a fixed pixel, which
  doubles the subharmonic reach) does not raise it. Dropping the
  subharmonics loses ~30% of it, so the subharmonics are necessary and
  sufficient. This bears on Conflict C-01 (docs/andrews-crosscheck.md):
  in THIS simulation, neither constant is exact, and 2.07 is LOW by ~2x.
  (Caveat: the beam centre is a windowed centroid on the speckled spot;
  the definition can bias the absolute value by tens of percent, not 2x.)
- **Convergence variants** (filled, HV x0.3; 120 snapshots each):

  | variant | sigma2_I | <beta^2> | w_lt |
  |---------|----------|----------|------|
  | standard (base, 250 snapshots) | 0.412 | 2.48 m^2 | 3.12 m |
  | reference preset (15 screens, 2048 px) | 0.390 | 2.33 m^2 | 3.09 m |
  | wide x2 (fixed pixel) | 0.359 | 2.12 m^2 | 3.05 m |
  | no subharmonics | 0.254 | 1.71 m^2 | 2.96 m |

  The top three rows agree inside the ~15% Monte-Carlo noise. So the
  fidelity-2 answer is stable, and the fidelity-1 excess is not a sampling
  artefact of the wave optics.
- **Third leg (Andrews Ch. 12, DOI 10.1117/3.626196).** The untracked
  analytic index (0.136 / 0.374 for the filled cases) sits BELOW the
  fidelity-2 measurement (0.41 / 1.25). No analytic leg matches the wave
  optics in the filled regime; the wave optics sits between Andrews
  untracked and Dios.

## 5. The Fig. 5 replication, and a slant-geometry defect

`validation/uplink_sigma2i/dios_fig5_replication.py` replicates Dios et al. 2004, Fig. 5
(GEO uplink, 0.84 um, sigma2_chi against W0, elevations 90 and 30 deg) with
the vendored kernels. Estimator: sigma2_chi = var(ln I)/4 on the raw
samples, both legs.

- **The vendored fidelity-1 kernels are FAITHFUL to the paper.** The curve
  overlays the figure: 90 deg plateau 0.0298 (figure ~0.028), 30 deg plateau
  0.0938 after the slant correction below (figure ~0.095), the 30 deg curve
  crosses sigma2_chi = 1 at W0 ~ 0.068 m and the 90 deg curve reaches 0.96
  at W0 = 0.1 m, as printed.
- **A real slant-geometry defect in the wrapper.**
  `olb.turbulence.uplink_flux._flux_result` puts the airmass on Cn2 but
  keeps the VERTICAL height grid as the path coordinate, so the Dios path
  weights A(z), B(z) read z = h instead of z = h sec(zeta). The on-axis
  index then scales as sec^1 instead of sec^(11/6) in the small-w0 limit:
  -40% at 30 deg elevation, -13% at 60 deg (negligible in Section 3, where
  the on-axis term is tiny). The exact fix, validated against the figure:
  give the kernels the slant-mapped grid hs*sec with the ZENITH profile.
  The off-axis and wander integrals are insensitive (< 1%). A separate
  owner-gated fix task is running for this.
- **Validity flag algebra, for the record.** sigma2_I ~= 4 sigma2_chi
  (I = exp(2 chi)), so the Rytov criterion sigma_R^2 < 1 is
  sigma2_chi < 0.25 -- which is what `WEAK_FLUCTUATION_LIMIT` implements.

## 6. CONFIRMED: the fidelity-2 GEO leg of the replication

Fidelity 2 on the exact Fig. 5 case (GEO, 0.84 um, 90 and 30 deg,
120 snapshots per point; sigma2_chi = var(ln eta)/4 on the raw samples):

| point | fid-1 (slant-corr.) | fid-2 | Fig. 5 FFT-BPM (read off) |
|-------|--------------------:|------:|---------------------------|
| 90 deg, W0 = 0.01 | 0.023 | 0.026 | ~0.028 |
| 90 deg, W0 = 0.03 | ~0.10 | 0.171 | above the solid line |
| 90 deg, W0 = 0.10 | 0.956 | 0.680 | ~0.5 to 0.7 (saturated) |
| 30 deg, W0 = 0.01 | 0.075 | 0.109 | ~0.105 |
| 30 deg, W0 = 0.03 | 0.260 | 0.692 | above the solid line |
| 30 deg, W0 = 0.10 | 2.459 | 0.646 | ~0.5 to 0.7 (saturated) |

The fidelity-2 points reproduce the paper's FFT-BPM behaviour at every
station: they sit ON the reference at small W0 (30 deg, W0 = 0.01:
0.109 vs the ~0.105 asterisk), they OVERSHOOT the weak-theory curve in the
mid range (the focusing regime, visible in the figure as asterisks above
the solid line), and they SATURATE near sigma2_chi ~ 0.65 at large W0
while the fidelity-1 curve climbs without limit (0.96 / 2.46). The GEO
wander ratio repeats the LEO finding: sim <beta^2> = 1.8x to 2.3x the Dios
form at every point. Plot: `figures/dios_fig5_replication.png`. So the fidelity-2
method is CONSISTENT with the exact published simulation that Dios
validated against, and Section 1 stands confirmed.

## 7. Proposals (owner-gated; none applied)

1. **A better fidelity-1 validity gate.** The current flag averages the
   per-sample sigma2_x, and the mean hides the broken tail: the filled
   HV x0.3 case reads "valid" (0.167 < 0.25) while over-predicting 2.3x.
   Gate on the typical wander radius instead: flag when
   sigma2_x(beta_rms) >= 0.25, or simply when beta_rms / w_L >= ~0.5.
   One line in `_flux_result`; changes a validity flag, so owner-gated.
2. **The slant-coordinate fix** (Section 5; separate task running).
3. **Model-of-record note.** For an uncorrected uplink with a filled launch
   (beta_rms / w_L over ~0.5), prefer fidelity 2; treat the fidelity-1 fade
   as pessimistic there. The docs (physics.md Section 5c, uplink docstrings)
   should say so once the pending run confirms.

## 8. What this unblocks

The original goal (the obscuration effect on uplink scintillation) can now
proceed against fidelity 2 ALONE: compare the obscured against the
unobscured psi_tx through the same far-field-map machinery. The fidelity-1
comparison is not a valid reference in the filled regime, and an obscured
launch is by construction a filled launch. The obscuration MEAN loss result
stands unchanged (tx_gaussian_efficiency_term, ~0.5 to 2.4 dB against the
vacuum far field).

## 9. Files

- `validation/uplink_sigma2i/uplink_farfield_reciprocity.py` -- the main experiment
  (four cases + variants; `--variants` reruns the key case lean).
- `validation/uplink_sigma2i/uplink_farfield_reciprocity_run.log` -- all four cases.
- `validation/uplink_sigma2i/uplink_farfield_variants_run.log`,
  `uplink_farfield_reciprocity_results_variants.json` -- key case +
  variants (full JSON incl. the sigma2_I(r) profiles and eta samples).
- `validation/uplink_sigma2i/dios_fig5_replication.py`, `dios_fig5_replication_run.log`,
  `dios_fig5_replication_results.json` -- the Fig. 5 replication.
- Earlier artifacts: `uplink_obscuration_dios_vs_waveoptics.py` (the MEAN
  panel stands; the scintillation panel remains void),
  `uplink_obscuration_farfield.py` (fine, vacuum only).
