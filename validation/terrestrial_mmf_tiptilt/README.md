# 10 km terrestrial MMF link: 1 inch and 2 inch terminals

Status: 2026-10-07. Fidelity-2 split-step wave optics, 1000 independent
snapshots for each case, 2048 px grid, on the bigfraw GPU.

## Summary

A 2 inch terminal decreases the mean launch-to-fibre loss by 14.5 to 17 dB
on a 10 km, 1550 nm horizontal link. It does not decrease the deep fades.

- **The collected power gives the gain.** The 36.3 mm launch has a radius of
  0.27 m at 10 km. The 9.1 mm launch has a radius of 1.09 m. Thus the 2 inch
  aperture collects approximately 16 dB more power.
- **The fibre coupling is not the problem when the focal length agrees with
  the aperture.** The 105 um core catches all of the turbulent spot in each
  trial (`figures/fibre_tip_*.png`). The fibre NA sets the coupling loss:
  1.2 dB at f = 50.24 mm on 1 inch, 0.1 dB at f = 200.96 mm on 2 inch.
- **The fade is a bucket fade.** p95 is 5 to 11 dB more loss than the mean.
  The power that gets to the aperture sets it, not the coupling. Receive
  tip-tilt gives less than 0.1 dB at each strength.
- **The 2 inch tail is shorter at 1e-15 and 5e-15 and longer at 1e-14.** No
  aperture averages the large speckle. At 1e-14 the narrower 2 inch beam
  wanders by approximately its own radius (1.06 W), and this adds deep fades.

## Setup

| | 1 inch family | 2 inch family |
| --- | --- | --- |
| Aperture, both ends | 25.4 mm | 50.8 mm |
| Launch 1/e^2 diameter | 4 x 2.27 = 9.08 mm | 16 x 2.27 = 36.32 mm |
| Beam radius at 10 km (vacuum) | 1.09 m | 0.27 m |
| Coupling focal length (read) | 4 x 12.56 = 50.24 mm | 16 x 12.56 = 200.96 mm |
| NA limit D / (2 NA) | 57.7 mm | 115 mm |
| Grid pixel / px across aperture | 3.9 to 4.8 mm / 5 | 1.2 to 2.6 mm / 20 to 42 |
| Campaign folders | `cn2<Cn2>_standard_focal_n2048` | `cn2<Cn2>_standard_focal_n2048_ap2in` |

- **Link:** 10 km horizontal path, 1550 nm, constant Cn2 of 1e-15, 5e-15 and
  1e-14 m^-2/3, outer scale 25 m, no absorption. A collimated launch. The
  core is at the nominal focal plane.
- **Fibre:** 105 um core multimode, NA 0.22. The focal length is a design
  value. The script does not optimise it.
- **Coupling read:** post hoc from the stored receive field of each trial
  (`run.CropEta`), so the same atmospheres give each focal length and each
  tip-tilt case. The read focuses the stored CROP, not the full grid, because
  the auto upsample of a short focal length makes a 2048 px grid 32768 px.
  Tip-tilt is a perfect 3-mode correction (Noll, DOI 10.1364/JOSA.66.000207),
  sensed from the field slopes or from the far-field centroid (G-tilt).
- **Loss** is launched power to fibre power, in positive dB. It includes the
  beam spread.
- The in-run focal length (the campaign key) is not the read focal length:
  1.21 m for the 1 inch family, 50.24 mm for the 2 inch family. It sets only
  the in-run scalar `mmf_eta`.

## Results

Launch-to-fibre loss in dB, no tip-tilt.

| Cn2 | Aperture | Mean | Median | p95 | p99 | Mean eta |
| --- | --- | --- | --- | --- | --- | --- |
| 1e-15 | 1 inch | 36.64 | 37.56 | 43.04 | 44.96 | 0.745 |
| 1e-15 | 2 inch | 19.55 | 20.12 | 24.87 | 26.80 | 0.975 |
| 5e-15 | 1 inch | 37.55 | 40.40 | 48.38 | 50.52 | 0.748 |
| 5e-15 | 2 inch | 21.67 | 23.75 | 31.43 | 34.51 | 0.973 |
| 1e-14 | 1 inch | 37.68 | 40.97 | 48.13 | 49.92 | 0.713 |
| 1e-14 | 2 inch | 23.04 | 26.26 | 34.33 | 37.32 | 0.969 |

The bootstrap standard errors are 0.1 to 0.2 dB on the mean and 0.3 to
0.7 dB on p99. The full statistics, with the two tip-tilt cases, are in
`results_standard_focal_n2048_f50mm.json` and
`results_standard_focal_n2048_ap2in_f201mm.json`. The 1 inch 4096 px check
is in `results_standard_focal_n4096_f50mm.json`.

### The geometric loss

The geometric loss is the vacuum loss from the launched power to the power
in the receive aperture (`geometric_split.py`, which writes a `geometric`
block into each results JSON). It is deterministic, so each quantile moves by
the same constant. It has three parts:

- **Spread:** the Gaussian beam into the aperture, 1 - exp(-2 a^2 / w(L)^2)
  (DOI 10.1117/3.626196, Ch. 4).
- **Truncation:** the launch aperture clips the 2 inch launch at
  a / w = 1.4, which decreases the on-axis far-field gain (Klein and Degnan,
  DOI 10.1364/AO.13.002134). The value is referenced to the power after the
  clip, as the campaign is. The 1 inch launch is not clipped.
- **Grid aperture:** the receive mask on the campaign grid is a set of whole
  pixels. At 5 px across the 1 inch aperture its area is 10 percent wrong,
  and the sign changes with the grid pixel, thus with Cn2. This is a sampling
  error of the run, not physics.

| Cn2 | Aperture | Spread | Truncation | Grid aperture | Total geometric | Bucket turbulence (mean) | Mean - geometric | Median - geometric | p95 - geometric | p99 - geometric |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1e-15 | 1 inch | 35.64 | 0.00 | -0.40 | 35.24 | 0.15 | 1.41 | 2.32 | 7.81 | 9.72 |
| 1e-15 | 2 inch | 17.63 | 1.24 | 0.00 | 18.87 | 0.57 | 0.68 | 1.25 | 6.00 | 7.93 |
| 5e-15 | 1 inch | 35.64 | 0.00 | +0.32 | 35.96 | 0.39 | 1.59 | 4.44 | 12.42 | 14.56 |
| 5e-15 | 2 inch | 17.63 | 1.24 | +0.03 | 18.90 | 2.66 | 2.77 | 4.85 | 12.53 | 15.61 |
| 1e-14 | 1 inch | 35.64 | 0.00 | +0.22 | 35.86 | 0.47 | 1.82 | 5.11 | 12.27 | 14.06 |
| 1e-14 | 2 inch | 17.63 | 1.24 | +0.06 | 18.93 | 4.00 | 4.11 | 7.34 | 15.41 | 18.40 |

All values are dB, no tip-tilt. "Bucket turbulence" is the mean collected
power less the geometric loss. "- geometric" is the launch-to-fibre loss less
the geometric loss, thus the turbulence and the fibre coupling together. A
wave-optics vacuum run (`olb.waveoptics.run.propagate_scenario`) gives
18.92 dB for the 2 inch family, 0.06 dB from the analytic total.

- **The 2 inch gain is geometric.** The geometric loss is 16.4 to 17.0 dB
  less. At 1e-14 the 2 inch turbulence costs 2.3 dB more on the mean and
  4.3 dB more on p99, so the net gain is 14.6 dB.
- **The 1 inch bucket turbulence is small on the mean** (0.15 to 0.47 dB):
  the 1.09 m beam is wide, so wander does not move the power off the
  aperture. The 2 inch beam loses 0.6 to 4.0 dB to the wander.
- **The grid aperture is the 2048 against 4096 px bias.** At 4096 px the
  term is -0.07 and -0.05 dB, and the mean less the geometric loss is 1.35
  and 1.36 dB, against 1.41 and 1.59 dB at 2048 px. Thus most of the 0.3 to
  0.6 dB grid bias in the 1 inch mean comes from the mask area.

The split of the loss shows its source. At 1e-15 on 1 inch, the collected
power gives 6.3 dB of the 6.4 dB between the mean and p95, and the coupling
gives 0.2 dB (eta from 0.66 to 0.79). On 1 inch the NA accepts the inner
76 percent of the aperture area. Receive tip-tilt changes no statistic by
more than 0.4 dB. The 0.4 dB is the slope sensor at 1e-14 on 1 inch, where
that sensor is not reliable on a 5 px pupil.

## The fade tail

Two effects set the tail: speckle that is too large for each aperture to
average, and beam wander, which is large for the narrower 2 inch beam.

| Cn2 (m^-2/3) | rho0, small speckle | L/(k rho0), large speckle | rms wander / beam radius, 1 -> 2 inch | collected-power SI, 1 -> 2 inch | p99 - mean, 1 -> 2 inch |
| --- | --- | --- | --- | --- | --- |
| 1e-15 | 6.7 cm | approx. Fresnel, 5.0 cm | 0.11 -> 0.33 | 0.55 -> 0.35 | 8.1 -> 7.3 dB |
| 5e-15 | 2.6 cm | 9.6 cm | 0.24 -> 0.75 | 1.90 -> 1.32 | 12.9 -> 12.7 dB |
| 1e-14 | 1.7 cm | 14.6 cm | 0.34 -> 1.06 | 2.21 -> 2.55 | 11.6 -> 14.2 dB |

rho0 is the spherical-wave coherence radius, rho0 = (0.545 Cn2 k^2 L)^(-3/5),
and L/(k rho0) is the large scale of the strong-fluctuation irradiance
(Andrews and Phillips, 2nd ed., DOI 10.1117/3.626196, Ch. 6 and Ch. 9). The
wander is `olb.turbulence.andrews.wander.beam_wander_variance` (Ch. 6) with
L0 = 25 m.

- **The averaging is small.** The Fresnel scale sqrt(L/k) is 5.0 cm. At
  1e-15 a 2 inch aperture holds approximately one speckle, and a 1 inch
  aperture holds a part of one. Two times the aperture decreases the
  scintillation index by only 37 percent.
- **No aperture averages the large speckle.** The Rytov variance is
  approximately 1.3, 6.7 and 13 at the three Cn2 (DOI 10.1117/3.626196,
  Ch. 5). In strong turbulence the irradiance has a small scale near rho0,
  which the 2 inch aperture partly averages, and a large scale of 10 to
  15 cm, which is wider than the two apertures. The deep fades are on the
  large scale.
- **Beam wander has an effect on the narrow beam.** The rms wander is
  approximately the same in metres (0.29 m and 0.37 m at 1e-14). That is
  1.06 beam radii for 2 inch and 0.34 for 1 inch. Thus the 2 inch receiver
  is frequently on the edge of the Gaussian, which gives a slow, deep fade
  that no aperture averages.
- **Result:** at 1e-15 and 5e-15 the averaging is stronger and the 2 inch
  tail is shorter. At 1e-14 the wander is stronger and the 2 inch p99 is
  2.6 dB more loss.

The wander split is analytic. The stored field patch is only 1.5 apertures
wide, so these runs cannot measure the beam centroid.

## Implications and next steps

For the smallest power fluctuations, launch a small beam (wide at the
receiver, so that the wander is a small fraction of it) and collect with a
large aperture. The aperture recovers the mean power.

- **Aperture size:** to average the large speckle, D must be more than
  approximately L/(k rho0), 4 to 6 inch at 5e-15 and 1e-14.
- **Fibre limit:** the multimode fibre accepts a turbulent spot only while
  D (lambda/rho0) <= 2 a NA (the etendue of a step-index fibre, Snyder and
  Love, DOI 10.1007/978-1-4613-2813-1). Thus D <= 0.25 m at 1e-14 and
  approximately 1 m at 1e-15 for this core and NA. The focal length must also
  increase as f >= D / (2 NA).
- **Wander:** transmit tip-tilt or more than one transmit beam removes the
  wander fade, and keeps the mean power. Receive tip-tilt does not help.
- **Proposed run:** keep the 9.08 mm launch and size the grid and the stored
  patch for a 6 inch receiver (`Campaign(sizing_aperture_m=...)`). Then 1, 2,
  4 and 6 inch all come from one run, post hoc, at 2048 px and the three Cn2.

Limits: snapshots only (the fade depth, not the duration); a perfect,
noise-free tip-tilt; 2048 px against 4096 px moves the 1 inch mean by 0.3 to
0.6 dB, and the sign changes with Cn2 (most of it is the pixel area of the
receive mask, see "The geometric loss"); the 2 inch family has no 4096 px
check.

## A 1 inch terminal seen at 45 deg

Status: 2026-10-20. `ellipse45.py`, 1000 snapshots for each case, on the
bigfraw GPU.

### Setup

- **Terminal:** a 1 inch aperture at both ends, seen at 45 deg (for example
  through a window or a fold). The aperture is then an ellipse of
  25.4 x 17.96 mm (the minor axis is 25.4 cos 45 mm), on TRANSMIT and on
  RECEIVE.
- **Launch:** a collimated Gaussian of 1/e^2 diameter 4 x 4.65 mm = 18.6 mm
  (the collimator beam through the 4x telescope). The circle clips 0.11 dB of
  it, the ellipse 0.39 dB.
- **Fibre:** 105 um core, NA 0.22, at f = 4 x 17.5 mm = 70 mm. The NA gate
  rho <= f NA = 15.4 mm is outside the aperture, so it is open and the fibre
  sees all of the pupil. The core is at the nominal focal plane.
- **Paths:** 10 km and 1 km, Cn2 1e-15, 5e-15 and 1e-14, L0 = 25 m.
- **The transmit ellipse** is the new `Transmitter(tilt_deg=45)`: the
  fidelity-2 launch clips the projected ellipse (`olb.waveoptics.run.
  _launch_clip`). An unset tilt leaves every campaign key alone. The analytic
  truncation Terms refuse a tilt.
- **The receive ellipse** is a post-hoc read of the stored patch
  (`ellipse45.EllipseRead`). The masks are built on the upsampled pupil, so
  the whole-pixel area error of "The geometric loss" does not enter. Each
  campaign gives the circle and the ellipse receiver from the same trials.
- **Loss** is from the telescope output (before the launch clip) to the
  fibre, in positive dB.
- **Grid:** the sized side, with 4096 px at 10 km (1.0 to 1.6 mm, 11 to 18 px
  across the minor axis) and 1024 px at 1 km (0.4 mm). A vacuum run at 2x the
  pixels moves the loss by 0.13 dB or less at 10 km and 0.04 dB or less at
  1 km (`ellipse45.py --vacuum-check`).
- **Not run:** the circle launch at 10 km for 5e-15 and 1e-14 (the owner
  stopped them for the ellipse cases).

### Results

Telescope-to-fibre loss in dB. The bootstrap standard error is 0.01 to
0.03 dB on the 1 km mean, 0.09 to 0.25 dB on the 10 km mean, and 0.3 to
0.8 dB on the 10 km p99.

| Path | Cn2 | TX -> RX | Mean | Median | p95 | p99 | Mean eta |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 1 km | 1e-15 | circle -> circle | 11.15 | 11.15 | 11.57 | 11.73 | 0.984 |
| 1 km | 1e-15 | ellipse -> ellipse | 13.99 | 13.99 | 14.44 | 14.60 | 0.980 |
| 1 km | 5e-15 | circle -> circle | 11.20 | 11.23 | 12.20 | 12.56 | 0.984 |
| 1 km | 5e-15 | ellipse -> ellipse | 14.07 | 14.10 | 15.15 | 15.51 | 0.980 |
| 1 km | 1e-14 | circle -> circle | 11.30 | 11.35 | 12.75 | 13.31 | 0.984 |
| 1 km | 1e-14 | ellipse -> ellipse | 14.13 | 14.19 | 15.72 | 16.25 | 0.980 |
| 10 km | 1e-15 | circle -> circle | 31.21 | 31.92 | 37.44 | 40.00 | 0.983 |
| 10 km | 1e-15 | ellipse -> ellipse | 34.08 | 34.82 | 40.61 | 43.13 | 0.980 |
| 10 km | 5e-15 | ellipse -> ellipse | 34.70 | 37.45 | 46.17 | 48.97 | 0.978 |
| 10 km | 1e-14 | ellipse -> ellipse | 35.32 | 38.58 | 46.57 | 49.68 | 0.978 |

The mixed cases (one end circle, one end ellipse) are in
`results_ell45_standard.json` (not tracked).

- **The 45 deg terminal costs 2.8 to 2.9 dB on the mean,** at each path and
  each Cn2. The receive ellipse gives 1.48 to 1.51 dB, which is the area
  factor cos 45 (1.51 dB). The transmit ellipse gives 1.34 to 1.37 dB: 0.28 dB
  of more clip, and the rest because the narrow axis spreads the far field.
  At 10 km and 1e-15 the p99 cost is 3.1 dB, a little more than the mean.
- **The fibre does not set the loss.** Mean eta is 0.978 to 0.984 at each
  case, and the NA gate is open. The turbulence on the mean is in the
  collected power.
- **At 1 km the fades are small:** p99 is 0.6 to 2.1 dB more loss than the
  mean. The Rytov variance is approximately 0.02 to 0.2 (weak).
- **At 10 km the fades are deep:** p99 is 9 dB (1e-15) and 14 dB (5e-15 and
  1e-14) more loss than the mean. Against the vacuum value (33.89 dB, ellipse
  at both ends), the turbulence adds 0.2, 0.8 and 1.4 dB to the mean.
- **This design is better than the 1 inch family of "Results".** With the
  circle at both ends at 1e-15 the mean is 5.4 dB less loss (31.21 against
  36.64 dB) and p99 is 5.0 dB less. The wider launch gives a narrower beam at
  10 km (0.53 m against 1.09 m radius), and f = 70 mm keeps the NA gate open
  (eta 0.98 against 0.745). With the ellipse at both ends the mean is still
  2.4 to 2.9 dB less loss than that circular family, and p99 is 0.2 to
  1.8 dB less.

Limits: snapshots only; no tip-tilt read (at f = 70 mm the tilt moves the spot
a few um against a 52.5 um core radius); a perfect, flat 45 deg projection
(no window, no fold aberration, no polarisation); the vacuum check covers the
circle receiver only.

`prelim_ellipse.py` draws the vacuum picture of the case
(`figures/prelim_ellipse45_L10km.png`, `figures/prelim_ellipse45_L1km.png`).
At 10 km its receive-pupil panel samples the field at 15.5 mm, so the cross in
that panel is the interpolation, not the beam; the beam is flat to 0.1 percent
over the aperture.

## Reproduce

The commands run from the repository root on bigfraw (the GPU venv). The
campaigns ran with the cupy backend, and the backend is part of the campaign
key, so `--fft-backend cupy` is necessary to open them. With the trials on
disk, the run makes no new trial and only reads the coupling.

```bash
python -m validation.terrestrial_mmf_tiptilt.run --n 2048 --fft-backend cupy --workers 8
```

```bash
python -m validation.terrestrial_mmf_tiptilt.run --n 2048 --fft-backend cupy --workers 8 --aperture-in 2
```

```bash
python -m validation.terrestrial_mmf_tiptilt.fibre_tip --n 2048
```

The 45 deg terminal (all cells, then the 10 km ellipse launch only, as run):

```bash
python -m validation.terrestrial_mmf_tiptilt.ellipse45 --length-km 1 10 --fft-backend cupy --workers 8
```

```bash
python -m validation.terrestrial_mmf_tiptilt.ellipse45 --length-km 10 --tx ellipse --fft-backend cupy --workers 8
```

The geometric split reads only the results JSON files, so it runs locally:

```bash
"C:\Users\alexf\anaconda3\envs\olb\python.exe" -m validation.terrestrial_mmf_tiptilt.geometric_split
```

The campaign data (`campaigns/`, about 170 npz blocks) stays on bigfraw and
git does not track it. The other scripts of this folder are the earlier
studies: `smoke.py` and `analyse_smoke.py` (the sampling smoke runs),
`compare_grids.py` (2048 against 4096 px), `focal_sweep.py` (the 1 inch
focal-length sweep) and `rekey.py` (the campaign re-key after the MMF
`optimal_focus` removal).
