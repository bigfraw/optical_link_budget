# uplink_divergence: the fidelity-2 uplink with a deliberate divergence

The follow-up of [../divergence_sampling/](../divergence_sampling/README.md)
(2026-09-23, branch `uplink-divergence`). That study found that a diverged
fidelity-2 uplink does not converge with the pixel count. This study finds the
cause, fixes it for this validation, checks the fix against a direct upward
propagation, and gives the power distributions up to 200 urad. The fidelity-1
(Dios) comparison of the handover plan was DROPPED by the owner.

## The configuration

`config.py`: a 15 cm ground aperture with no obscuration, a FILLED launch
(waist 55 mm, so D = 2.73 w), 1550 nm, a 500 km orbit, elevations 30 and
60 deg, the site HV5/7 atmosphere, L0 = 25 m. The divergences are 2, 4, 8 and
16 times `theta_min = lambda / (pi w) = 8.97 urad` (17.9, 35.9, 71.8,
143.5 urad) and 200 urad, plus the collimated launch.

## Verdicts

1. **THE CAUSE IS THE PLANE-WAVE START.** The slab starts from a plane wave
   that fills the grid. The absorbing mask cuts it at every hop, and the hop
   count follows the pixel count (`max_step = N dx^2 / lambda`), so the
   Fresnel-ring pattern in the aperture changes with the pixel count. A curved
   `psi_tx` is a chirp and it reads the rings. Against the exact vacuum
   answer, the plane-wave overlap errs by -16 to +10 dB and does not converge
   from 256 to 2048 px (`gate0.py vacuum`).
2. **A WIDE GAUSSIAN START FIXES IT.** A Gaussian of waist 0.15 to 0.2 x the
   grid side makes no rings, and its masked split-step vacuum matches the
   closed form (GForvard) within 0.02 dB at every pixel count.
3. **ARM A: THE GAUSSIAN START MATCHES A DIRECT UPWARD PROPAGATION.** The truth
   launches `psi_tx` UP through the same screens on a wide grid and reads the
   satellite axis with the Fresnel sum. The hand downward loop equals the
   production runner bit for bit (step 0). 200 trials for each case:

   | start | mean eta vs truth | correlation | rms dB per trial |
   | --- | --- | --- | --- |
   | Gaussian, 0.3 x side | within 3 % | 0.989 to 0.998 | 0.13 to 0.34 |
   | Gaussian, 0.2 x side | within 5 % | 0.97 to 0.995 | 0.23 to 0.49 |
   | plane wave (of record) | up to 1.6 dB high | down to -0.05 | 0.5 to 5.0 |

   (collimated, 2x and 4x; 30 and 60 deg on the rapid grid, and 60 deg at 2x
   the pixel count, which agrees, so the fix converges.) 0.3 x side is the
   choice. The mean eta at 4x is 0.98 to 1.02 in the TRUTH, so a mean near 1
   is physical. Arm A stops at 4x: above that the upward beam needs a wide grid
   of tens of thousands of pixels.
4. **THE OVERLAP HAS NO CONJUGATE, AND THE LAUNCH MODE DIVERGES.** Rechecked on
   owner pushback: `psi_tx` grows as it propagates (37 to 62 mm rms over 1 km
   at 36 urad), and through one strong screen the no-conjugate overlap reads
   0.28 / 0.40 against a direct-propagation truth of 0.32 / 0.55 (36 / 144
   urad), where the conjugate reads 0.048 / 0.024.
5. **ARM B: THE TURBULENCE FADE IS STABLE UP TO 200 URAD.** 2000 trials for
   each elevation (standard preset, 1024 px at 30 deg, 512 px at 60 deg,
   `dx <= lambda / (4 theta_max)`), cupy on bigfraw, the divergences read post
   hoc from ONE campaign, so they share their atmospheres. Past 4x the
   spread holds: the p5 fade is about 5 dB (30 deg) and 3 dB (60 deg) from
   36 to 200 urad. The MEAN loss oscillates with the divergence (30 deg: 0.25 /
   0.99 / 0.45 / 0.93 dB at 36 / 72 / 144 / 200 urad), and a screen-free
   Kolmogorov far-field estimate (`farfield_check.py`) gives the same pattern
   (0.10 / 0.86 / 0.10 / 0.78 dB). So the oscillation is physics, below.
6. **THE HARD-CLIP FAR-FIELD RIPPLE IS PHYSICAL.** The clip at 1.36 w leaves
   16 % of the peak amplitude at the edge. On the axis the far field is the
   sum of the beam centre and the edge-diffraction term, and their relative
   phase `k a^2 / (2 R)` turns with the divergence, so the on-axis vacuum power
   swings by the factor `|1 - q exp(i phi)|^2`, `q = exp(-a^2/w^2) = 0.156`:
   -1.47 to +1.26 dB, period about 30 urad. The turbulence averages the
   ripple out, so `eta` (turbulent / vacuum ON THE AXIS) ripples the other way.
   It needs a sharp round edge and an exact divergence, so it is an idealised
   limit.
7. **A BUDGET BUG: THE TRUNCATION TERM IGNORES THE DIVERGENCE.**
   `tx_gaussian_efficiency_term` charges a FIXED 1.47 dB, the on-axis clip
   penalty of a COLLIMATED beam, at every divergence. The exact clipped
   launch (`total_loss.py`) matches the analytic geometric Terms to 0.01 dB
   when collimated, and it reads 0.7 to 2.7 dB LESS loss for every diverged
   case (the ripple around the smooth curve, mean -0.1 dB). So a diverged
   fidelity-2 budget over-reads the total by that much. NOT FIXED.

## The total loss (the consistent route: exact clipped vacuum + turbulence)

| p50 / p99 (dB) | coll | 36 urad | 72 urad | 144 urad | 200 urad |
| --- | --- | --- | --- | --- | --- |
| 30 deg | 51.7 / 64.6 | 60.9 / 66.5 | 65.9 / 71.2 | 72.1 / 77.7 | 74.6 / 80.1 |
| 60 deg | 46.3 / 53.4 | 56.3 / 59.7 | 61.1 / 64.5 | 67.5 / 71.3 | 69.8 / 73.2 |

The receiver is a 5 cm point-like aperture on the satellite. The budget route
(analytic geometric + truncation + turbulence) is in `figures/total_loss_cdf.png`
next to it.

## Open items

- The Gaussian start is a VALIDATION PATCH (`gate0.start_field` replaces
  `run.Begin`). The production runner still starts from the plane wave, and
  the campaign fingerprint does not record the start. A production option
  (default off, append-only fingerprint) is the next step, and whether it
  becomes the default for a diverged uplink is an OWNER decision.
- The truncation Term fix (verdict 7).
- The space sizer does not read the launch curvature (the
  `dx <= lambda / (4 theta)` guard of `divergence_sampling/` is NOT BUILT).
- 8x and above rest on the pixel rule and the trend, not on Arm A.

## Files

| File | Purpose |
| --- | --- |
| [config.py](config.py) | The shared scenario, geometry and divergence list. |
| [gate0.py](gate0.py) | The start-field patch, the vacuum convergence test and a turbulent pixel sweep. |
| [arm_a.py](arm_a.py) | The direct upward-propagation check (step 0 plus a process pool; run on bigfraw). |
| [arm_b.py](arm_b.py) | The campaigns (`run`, cupy on bigfraw) and the post-hoc read (`read`); writes `figures/arm_b_ccdf.png`. |
| [farfield_check.py](farfield_check.py) | The screen-free expected mean loss against the divergence. |
| [total_loss.py](total_loss.py) | The total-loss CDFs by the two routes; writes `figures/total_loss_cdf.png`. |

The campaigns (`campaigns/`, about 0.4 GB) and the `.npy` / `.json` outputs
stay on bigfraw under `D:/repos/optical_link_budget/validation/uplink_divergence/`;
they are git-ignored.

Run each script from the repository root as a module, for example
`-m validation.uplink_divergence.arm_b read`, with the olb env interpreter
(the GPU venv for `arm_b`).
