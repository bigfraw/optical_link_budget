# Next: tip-tilt at both ends of the 45 deg link

Status: PLAN, 2026-10-20. Nothing is built. Do the steps in order; each step
can stop the work.

## The question

The 45 deg terminal (`ellipse45.py`, README section "A 1 inch terminal seen at
45 deg") has deep fades at 10 km: p99 is 9 dB (1e-15) and 14 dB (5e-15 and
1e-14) more loss than the mean. Does tip-tilt at both ends remove a large part
of them?

ONE MEASUREMENT, TWO USES. Each terminal senses the tip-tilt of the beam that
it RECEIVES; that is the only tilt it can act on, and it is not zero. The
same angle theta drives the receive mirror and the transmit mirror. What
differs is the lever arm:

- **On the fibre tip** the spot moves by f theta. At f = 70 mm and tens of
  urad, that is approximately 1 to 2 um against a 52.5 um core radius, and
  the mean eta is already 0.98. So the receive correction changes the
  coupling by little: the tilt is real, but it does not move the spot
  appreciably on the core.
- **At the far end** the same angle on the launch moves the beam by L theta.
  At L = 10 km and tens of urad, that is 0.1 to 0.3 m against a 0.53 m beam
  radius: the wander. The transmit mirror points the launch along the
  conjugate of the sensed tilt. By reciprocity (Shapiro,
  DOI 10.1364/JOSA.61.000492) this is a tip-tilt PRE-COMPENSATION, and it
  removes the wander. Expected gain: real, size unknown.
- Order of magnitude of theta: the arrival-angle variance of a 1 inch
  aperture on a 10 km path is approximately 2.91 Cn2 L D^(-1/3) for a plane
  wave, less for a spherical wave (Andrews and Phillips, Ch. 6,
  DOI 10.1117/3.626196): approximately 5 to 10 urad rms at 1e-15 and 20 to
  30 urad at 1e-14. Step 1 measures it.
- A horizontal link is the CLEAN case: the terminals do not move (no
  point-ahead), and the two ends have the same aperture, so the beacon and
  the launch see the same pupil.

## Step 1: receive tip-tilt from the stored campaigns (no new trial)

Read a perfect receive tip-tilt (3 Noll modes, Noll DOI 10.1364/JOSA.66.000207)
from the stored crops of the 10 km campaigns, as `run.CropEta` did for the
old families.

- The circle receiver: `_PostCorrector(patch, [TipTilt()], 0.0254, 0.0,
  "gtilt")` (and "slopes" as the cross-check), as in `run.CropEta`.
- The ellipse receiver: `_PostCorrector` builds a CIRCULAR mask. Use
  `olb.waveoptics.compensation.modal.ApertureModes(3, N, ellipse_mask)` (the
  least-squares fit takes any mask), or remove the G-tilt (the mean wrapped
  phase gradient) over the ellipse.
- Add the result as two more columns to `EllipseRead`; the read is a few
  seconds for each cell. It runs on bigfraw (the npz blocks are there).
- Also store the sensed tilt angle of each trial (rad, two axes). Its rms
  times f is the spot motion on the tip, and its rms times L is the beam
  motion that the same angle gives at the far end. That one number connects
  step 1 to step 3.
- **The wander from the tilt, statistically.** The forward-beam wander comes
  from the turbulence near the TRANSMITTER (weight (1 - z/L)^2), and the tilt
  that the far end senses comes from the turbulence near the RECEIVER. So the
  far-end tilt does not follow the forward wander trial by trial; the
  NEAR-end tilt (of the reverse beam) does, by reciprocity. But the path has a
  constant Cn2 and the two ends are the same, so the two ends have the same
  tilt statistics: L x rms(far-end tilt) estimates the rms wander. Compare it
  with `beam_wander_variance` (Ch. 6, DOI 10.1117/3.626196, at L0 = 25 m).
  At 5e-15 and 1e-14 (Rytov 6.7 and 13) the beam breaks up, so the centroid
  wander is less clear. Step 1b measures the wander directly
  (`rx_moments`), and step 2b compares it with L x the near-end tilt, trial
  by trial.
- **Expected:** the TT eta moves the mean and p99 by less than 0.2 dB,
  because f theta is small against the core. That closes the RECEIVE
  correction as a coupling fix; it does not close the tilt, which step 3
  uses.

## Step 1b: collect the wander of the stored trials (data collection)

The stored patch is 19 mm in radius, and the 10 km beam is 0.53 m in radius,
so the stored fields cannot give the beam centroid (the speckle sets the
intensity inside the patch, not the slope of the displaced Gaussian). The
seed regenerates each atmosphere bit for bit on the same backend, so
RE-PROPAGATE the stored trials and keep the receive-plane moments of each
trial. They then pair, trial for trial, with the collected power that is
already stored.

**What to collect, for each trial** (one file `wander_ell45.npz` in each
campaign folder, one row for each trial):

| Value | What it is |
| --- | --- |
| `cx`, `cy` | the power centroid of the receive-plane beam, in m |
| `w_st` | the short-term beam radius, sqrt(2 <\|r - c\|^2>) about the centroid, in m (a Gaussian exp(-2 r^2 / W^2) gives <r^2> = W^2 / 2) |
| `p_grid` | the power on the grid (the absorbing boundary takes the rest) |
| `collected_power` | the bucket power again, for the identity check |

The moments read a window of radius 3 W(L) about the axis (W(L) the vacuum
beam radius), so the light that the turbulence scatters to the grid edge
does not pull the centroid.

**Which cells** (all are on bigfraw):

| Path | Cells | GPU time (one forward pass, as the first run) |
| --- | --- | --- |
| 10 km | ellipse launch at 1e-15, 5e-15, 1e-14; circle launch at 1e-15 | approx. 57 + 72 + 102 + 57 min = 4.8 h |
| 1 km | all six (circle and ellipse at the three Cn2) | approx. 20 min |

**How.** Add one opt-in runner option, `rx_moments=False` on
`propagate_turbulent_scenario`, that puts (cx, cy, w_st, p_grid) on each
`TurbTrial` (step 2b and step 3 need the same numbers, so build them one
time). The receiver is an MMF, so the receive field is on the host already
(the MMF host tail); the moments are some sums. Default off, so every
stored record and key stays the same. Do NOT reopen the stored campaigns
with it (their manifests do not know the option): call
`propagate_turbulent_scenario` directly with the campaign scenario, grid,
plan, `seed=SEED`, `fft_backend="cupy"`, and `start_index` block by block,
the same arguments that `Campaign.run` gives (`Campaign._runner_kwargs()`).
A validation script `wander_split.py` in this folder does that and writes
the npz. Launch it on bigfraw through WMI with a log, as `ellipse45.py`.

**The identity gate (first, on one block).** The regenerated
`collected_power` must equal the stored one for each trial (the same
backend, so to the float32 rounding at most). If it does not, the
atmospheres are not the same and the pairing is invalid: stop.

**The analysis** (local, from the npz and the stored power):

- r = |c|, and the wander loss of each trial, L_w = 10 log10(e) 2 r^2 /
  w_st^2 dB (the aperture is small against w_st, so the bucket reads the
  local Gaussian irradiance; Andrews and Phillips Ch. 6 and 8,
  DOI 10.1117/3.626196).
- The beam-shape loss, from w_st against the vacuum W(L): the short-term
  spread.
- The remainder: the total turbulence loss less the two, which is the
  scintillation (the speckle on the aperture).
- **Figures** (`figures/wander_split_<case>.png`): the collected power (dB)
  against r / w_st with the model curve; the CDFs of the total fade, the
  wander part and the remainder; and the share of the p95 and p99 fade that
  the wander explains, against Cn2, for 1 km and 10 km.
- Compare rms(r) with `beam_wander_variance` (L0 = 25 m) and with
  L x rms(the far-end tilt) of step 1.

**What it answers:** how much of each fade is the wander (which the transmit
tip-tilt of step 3 can remove) and how much is the speckle (which it cannot).
That is the upper bound of the step 3 gain, from the trials that exist, at
the cost of one forward pass.

## Step 2: an analytic estimate of the wander part (minutes, local)

Before any runner change, estimate how much of the tail the wander gives.

- The wander variance: `olb.turbulence.andrews.wander.beam_wander_variance`
  and `pointing_error_variance` (Ch. 6, DOI 10.1117/3.626196), at L0 = 25 m,
  for the 9.3 mm launch at 10 km.
- The tracked and the untracked Gaussian-beam index:
  `olb.turbulence.andrews.scintillation.scintillation_index(..., wave=
  "gaussian", tracked=True|False)` (Ch. 8).
- Put both through the bucket over the 1 inch aperture and compare the p99 of
  the two.
- **Caution:** the Rytov variance is 1.3, 6.7 and 13 at the three Cn2, so
  5e-15 and 1e-14 are outside the weak-regime model. The estimate is a
  guide for 1e-15 only and an order of magnitude past it. It decides whether
  step 3 is worth the cost; it does not replace step 3.

## Step 2b: the tilt at the two ends, trial by trial (the diagnostic)

Before the correction, measure what the two ends see at the SAME instant.
This needs the reverse pass of step 3 but NO correction, so it is the first
thing to build, and it tests the reverse pass and the sign.

For each trial, on the same screens, propagate both directions UNCORRECTED and
store four values (two axes each):

| Value | What it is |
| --- | --- |
| `theta_far` | the tilt of the FORWARD beam over the far aperture (what the far terminal senses) |
| `theta_near` | the tilt of the REVERSE beam over the near aperture (what the near terminal senses) |
| `c_far` | the centroid of the forward beam at the far end (its wander) |
| `c_near` | the centroid of the reverse beam at the near end |

Questions, each a correlation over the trials:

1. **Is the same-instant tilt the same at the two ends?** corr(theta_near,
   theta_far). Expected: partial, because each end senses mostly its own
   end of the path.
2. **Which tilt follows the wander?** corr(L theta_near, c_far) against
   corr(L theta_far, c_far). Expected: the near-end tilt follows the forward
   wander (reciprocity, Shapiro DOI 10.1364/JOSA.61.000492); the far-end tilt
   does not.
3. **Are the statistics the same?** rms(theta_near) = rms(theta_far) and
   rms(c_far) = rms(c_near) for a uniform Cn2 and identical terminals (the
   mirror symmetry of the path).
4. **The sign (a gate).** The sign of the slope in question 2 tells which
   sign of the correction removes the wander (step 3, item 3).

**The weighting, made visible (cheap and instructive).** Run the same with
ONE screen and move it along the path, z/L = 0.05 to 0.95. Plot the four
variances and the two correlations against z/L. The wander must grow
towards the transmitter, and each tilt towards its own receiver (the path
weights of Andrews and Phillips Ch. 6, DOI 10.1117/3.626196). One screen at a
fixed r0 and a small grid make it seconds per trial.

**Figures** (`figures/tilt_two_ends_<case>.png`, one per path and Cn2):

1. **Scatter, theta_near against theta_far,** x axis and y axis in two
   panels, in urad, with the 1:1 line, the Pearson r with its bootstrap SE,
   and the two rms values in the title. This is THE plot of the asymmetry: a
   cloud on the 1:1 line means the two ends agree; a round cloud means they do
   not.
2. **Scatter, c_far against L theta_near and against L theta_far,** in m,
   side by side, with r on each. The left must line up (reciprocity); the
   right must not.
3. **Single-screen sweep:** r(theta_near, theta_far), r(L theta_near, c_far)
   and r(L theta_far, c_far) against the screen position z/L, with the
   bootstrap band; and the four variances on a second axis.
4. **Summary:** the three correlations against Cn2 for 1 km and 10 km on one
   plot, so the fall in strong turbulence is visible.

Use the G-tilt for the figures and print the Z-tilt r next to it.

**Where to run it.** First at 1 km (1024 px, weak turbulence, minutes for
1000 trials), where the centroid is clear. Then the three 10 km Cn2, where
the beam breaks up at 5e-15 and 1e-14 and the centroid is a less clear
wander.

**Code.** The reverse pass is item 1 of step 3. The tilt is the G-tilt (the
power-weighted mean phase gradient, so a centroid in the far field) over
each aperture, plus the Z-tilt (the least-squares plane, `ApertureModes`) as
the cross-check; the two differ in strong turbulence. The centroid is one sum
over the full receive grid. A validation script (no package change) can do
this if it calls `split_step` directly with the seeded screens of a stored
campaign; that keeps the runner untouched until step 3.

## Step 3: terrestrial transmit tip-tilt in the runner (the real test)

The runner pre-compensates a SPACE uplink only (`compensation="terminal"`
reads the screens of the downlink slab). A terrestrial path needs its own
reverse pass. The stored campaigns cannot answer this post hoc: the stored
patch is only 1.5 apertures wide, so it does not hold the beam centroid, and
the launch must change.

Proposed option (default OFF, so every stored key stays valid; the
fingerprint tail is append-only):

    propagate_turbulent_scenario(..., tx_compensation=None | [TipTilt()])

For each trial, on the SAME screens:

1. **Reverse pass (the beacon).** Launch the FAR terminal's beam (its own
   Transmitter, the same 45 deg ellipse) from z = L to z = 0: the screens in
   REVERSE order at the mirrored positions L - z_j.
2. **Sense.** Fit the tip-tilt of the beacon field over the NEAR launch
   aperture (the ellipse mask; the terrestrial slope source, as the
   receive-side compensation already does). This IS the near terminal's
   receive tilt: the same measurement that step 1 reads at the far end.
3. **Pre-compensate.** Remove that tilt from the near launch field: apply the
   CONJUGATE of the sensed tilt plane (the retraced wave goes back along the
   beacon path).
4. **Forward pass.** Propagate the corrected launch as now, and store the
   same scalars and patch.

Design points to decide:

- **Screen memory.** `split_step` takes the screens one at a time (a
  generator), so the reverse order needs either (a) a second draw from the
  seeds (the host draw is the cost of a trial today, so about 2x the trial
  time) or (b) all screens held on the device (35 x 4096^2 float32 =
  2.3 GB, which fits the 8 GB RTX 4070). Recommend (b) with (a) as the
  fallback.
- **Sign check (a gate).** One strong screen next to the near terminal: the
  pre-compensated beam must land centred at the far end (its centroid near
  zero, where the uncorrected one wanders). Add it to the runner self-check.
- **Measure the wander directly.** Turn on `rx_moments` (step 1b): the
  centroid with and without the correction is the direct test.
- **Flags.** PERFECT TT (no servo lag, no sensor noise, a snapshot), as the
  space perfect-AO route flags it. A real loop at a crosswind of some m/s
  gets less.

Cost on the bigfraw GPU, 1000 trials, 10 km, ellipse at both ends: today
57 / 72 / 102 min for the three Cn2. With the reverse pass, about 2x: approx.
2 / 2.4 / 3.4 h, 8 h for the three. The 1 km cells add some minutes.

## What the result answers

- The p99 and p95 of the ellipse terminal with TT at both ends, against the
  numbers of the README section (43.13 / 48.97 / 49.68 dB at p99).
- How much of the 9 to 14 dB tail is wander, and how much is the speckle that
  no tip-tilt removes (a tracked beam still scintillates).
- If the gain is large: whether a real tracking loop (bandwidth against the
  tilt frequency of the path) can get it. That is a temporal question (the
  frozen-flow axis is built for the space downlink only, backlog 2-P1b).
