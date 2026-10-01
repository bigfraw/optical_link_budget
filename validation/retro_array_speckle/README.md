# The target speckle of a multi-cube retroreflector array

Date: 2026-10-01. Branch `spherical-point-ahead`. Machine: the laptop.
Status: a DESIGN RECORD. Nothing is simulated yet. The numbers below are
ORDER-OF-MAGNITUDE ESTIMATES, not results, except the `lageos_pass.py` output.

## The purpose

olb models a space retro link with ONE corner cube
(`olb.links.retro_space`). A real geodetic satellite carries MANY cubes:
LAGEOS-1 (NORAD 8820) carries 426 (422 fused silica, 4 germanium) on a 60 cm
sphere. This study records three things:

1. Why the single-cube velocity-aberration Term does not apply to such an
   array.
2. The mechanism that a single cube cannot show: the light of the cubes
   INTERFERES at the station, so the return FADES with no atmosphere at all
   (target speckle).
3. The simulation that measures that fade, and how it combines with the
   fidelity-2 retro record.

## The trigger: LAGEOS-1 over Yarragadee, 2026-10-01

`lageos_pass.py` propagates the CelesTrak TLE (epoch 2026-09-30 12:56 UTC)
from MOBLAS-5 (29.0464 S, 115.3467 E, 244 m). The TLE matches the RC700 mount
readout (RA 18h17m54.8s, Dec -49d33m28s, J2000) at 05:19:50 UTC to 8.3 arcmin,
at 34.4 deg elevation and 7450 km range.

| | 532 nm | 1064 nm |
| --- | --- | --- |
| Point-ahead angle at the readout | 35.6 urad (7.35 arcsec) | same |
| Point-ahead angle over the pass | 6.7 to 8.0 arcsec | same |
| `x = pi D theta / lambda`, one 3.81 cm cube | 8.02 | 4.01 |
| Airy loss of one UNSPOILED cube, at the readout | 24.6 dB | 29.3 dB |
| Airy loss over the pass | 23.8 to 32.1 dB | 21.0 to 79.4 dB |

The first Airy null is at x = 3.83. At 532 nm the station sits in the outer
rings of the lobe; at 1064 nm it sits on the null, so a small angle change
moves the loss by tens of dB. These numbers are NOT the LAGEOS loss. The
LAGEOS cubes have a dihedral-angle offset of about 1.25 arcsec (SPOILED
cubes), which turns the lobe into a ring of about the aberration angle. The
retro budget refuses a `SpoiledCornerCube` for this reason.

## The mechanisms

### 1. One cube is one small transmitter

A cube returns the light along the incoming ray in its own frame. Seen from
the station, the return lobe is offset by the aberration angle
`theta = 2 v_perp / c` (Degnan, DOI 10.1029/GD025p0133). Two properties of
cube i set what it gives the station:

- The STRENGTH, the cube cross section `sigma_i` toward the station. It
  falls with the incidence angle (the effective area of a tilted cube
  shrinks, and a cube stops returning past a cut-off angle), and it is the
  far-field pattern of the cube read at the offset `theta`. For an unspoiled
  cube the on-axis peak is `sigma = 4 pi A^2 / lambda^2` (Degnan,
  DOI 10.1029/GD025p0133). For one 3.81 cm cube at 532 nm that is
  5.8e7 m^2.
- The PHASE. The light makes the trip twice, so

      phi_i = 2 k (n . r_i) + psi_i

  with `k = 2 pi / lambda`, `n` the unit vector to the station, `r_i` the
  effective reflection point of the cube (inside the glass, behind the
  face), and `psi_i` the fixed phase of the reflections inside the cube.

### 2. The array factor: add the fields, not the powers

The field at the station is the sum of the cube fields. The cross section at
one instant is

    sigma_now = | sum_i sqrt(sigma_i) exp(i phi_i) |^2

That sum is the ARRAY FACTOR. Two equal cubes show the effect: in phase they
give 4 sigma_1, in anti-phase 0, and 2 sigma_1 on average (the two-slit
pattern).

### 3. The phases are effectively random: fully developed speckle

A path change of `lambda / 4` (0.13 um at 532 nm) moves the round-trip phase
by pi. The motion of the satellite moves the cube ranges by much more, so
with about 20 to 30 active cubes the sum is a random walk of phasors. Then
(Goodman, DOI 10.1364/JOSA.66.001145):

- The MEAN of `sigma_now` is `sum_i sigma_i`. That is the published lidar
  cross section (the ILRS quotes about 7e6 m^2 for LAGEOS; CHECK this value).
  So the lidar cross section is the MEAN of this fade, and it holds no fade.
- `sigma_now` is EXPONENTIAL: P(below 10 % of the mean) = 1 - exp(-0.1) =
  9.5 %, and P(below 1 %) = 1.0 %. The fades go to true nulls.
- One dominant cube (a cube that looks straight at the station) makes the
  law RICIAN, with shallower fades.

A cross-check of the order: 5.8e7 m^2 per cube, times the 24 dB Airy loss,
times about 25 active cubes, gives a few 1e6 m^2. That is the order of the
published value.

### 4. The fade is ONE number across the receive aperture

Across the ground, the phase differences change and the speckle repeats at
about `lambda R / D_act`, with `D_act` the projected size of the active cap.
For an active half-angle of 40 deg, `D_act = 2 R_sphere sin 40 = 0.39 m`, and
the speckle at 7450 km and 532 nm is about 10 m across. A 1 m telescope sits
inside one speckle, so aperture averaging does not help.

### 5. The time scale

The speckle changes fully when the view direction turns by about
`lambda / (2 D_act)` = 0.69 urad. With the attitude fixed, the inertial line
of sight turns at `v_perp / R` = 0.72 mrad/s at the readout, so the speckle
time is about 1 ms. A spin makes it shorter. The spin state of LAGEOS-1 is a
simulation INPUT; do not assume a value.

### 6. What reduces the fade

| Effect | What it does | Status for LAGEOS |
| --- | --- | --- |
| Polarization | Uncoated (TIR) cubes depolarize; a polarization-blind detector adds two independent speckles: gamma of order 2, P(below 10 %) = 1.75 %, P(below 1 %) = 0.02 % (Goodman, DOI 10.1364/JOSA.66.001145) | Applies (uncoated fused silica). The real depolarization of each cube is an input |
| Integration time T | Averages `M = T / tau_c` speckles; the contrast falls as `1 / sqrt(M)` | Ranging over seconds: averaged. Comms bits (ns): no effect |
| Laser linewidth | A linewidth above `c / (round-trip depth spread)` decorrelates the cubes | The active cap is about `R (1 - cos 40) = 7 cm` deep, 14 cm round trip, so a linewidth above about 2 GHz averages. Comms lasers (kHz to MHz) do NOT |
| Several wavelengths | Each wavelength speckles independently | A design option |
| Aperture averaging | Needs a telescope near the speckle size (10 m) | No effect |

### 7. The pulse case (ranging)

A ranging pulse of 10 to 50 ps is 3 to 15 mm long. That is much shorter than
the 14 cm round-trip depth spread, so the returns of cubes at different
depths do not overlap in time and their POWERS add. Only the cubes inside
about one pulse length of the same range interfere. The return is a
stretched pulse made of a few speckling groups: a weaker fade, but a changed
pulse shape, which moves the range centroid. That is the LAGEOS "target
signature" and its centre-of-mass correction of about 250 mm (Otsubo and
Appleby, DOI 10.1029/2002JB002209).

## Why ONE wave-optics run covers all 426 cubes

The atmosphere is the SAME for every cube on both legs:

- UP LEG. The turbulent patches of the beam at the satellite are of order
  `lambda R / r0` (about 40 m at 7450 km and r0 = 10 cm). The array is
  0.6 m. Every cube receives the same up-leg power.
- RETURN LEG. The array spans 0.6 m / 7450 km = 0.08 urad from the ground.
  The isoplanatic angle is a few urad (Fried, DOI 10.1364/JOSA.72.000052).
  Every cube's light crosses the same turbulence.

So for each trial the field at the receive aperture FACTORIZES:

    U_rx(x) = U_point(x) * AF(t),   AF(t) = sum_i sqrt(sigma_i) exp(i phi_i(t))

`U_point` is the turbulent return of ONE point source, which is the existing
fidelity-2 retro record (`direction="retro"`, with the Shapiro reciprocity
overlap for the up leg, DOI 10.1364/JOSA.61.000492). `AF` needs no
propagation: 426 complex multiply-adds per time step.

The factorization FAILS only when reflectors are far enough apart to see
different turbulence: tens of metres on the return leg at this range
(isoplanatic angle times R). Then each group takes its own window of the
SAME drawn atmosphere (the point-ahead shifted-window route), not a new
atmosphere. No single-array spacecraft comes near that.

## The simulation plan

### Stage A: the array factor, no atmosphere (CPU, seconds)

Inputs:

- A cube catalogue: the position, normal, type and dihedral offset of each
  cube. Source for LAGEOS: Fitzmaurice et al., NASA TP-1062 (1977), no DOI.
- A cube far-field model: the strength `sigma_i` against the incidence angle
  and the offset `theta`. First cut: the effective area of a tilted cube
  times the Airy lobe. Real model: the FFT of the six-sector cube pupil with
  the dihedral phase of each sector (the spoiled ring). Arnold, SAO Special
  Report 382 (1979), no DOI, gives the method.
- The attitude: the spin axis and the rate.
- The line of sight from a `TLEPass` (the 2026-10-01 pass is in
  `lageos_pass.py`).
- The wavelength, and the mode: coherent (CW), polarization-blind (two
  independent sums), or pulsed (power sum over range bins of one pulse
  length).

Outputs: the time series of `sigma_now`, its mean, its CDF, and its
correlation time.

Gates:

| Gate | Test | Pass |
| --- | --- | --- |
| A1 | ONE cube | No speckle: `sigma_now` is the cube pattern, constant at a fixed geometry |
| A2 | N equal cubes, random phases | The CDF is exponential (Kolmogorov-Smirnov) |
| A3 | Pulse mode with a long pulse | Equals the coherent mode; with a short pulse the mean equals `sum sigma_i` with a weak fade |
| A4 | LAGEOS mean `sigma_now` | The order of the published cross section (the exact ILRS value to be checked) |
| A5 | LAGEOS pulse mode | The range centroid gives the published centre-of-mass correction, about 250 mm (Otsubo and Appleby) |

### Stage B: the time scale

Sweep the attitude (no spin, slow spin, fast spin) over the 2026-10-01 pass.
Record the correlation time of `sigma_now` against the 1 ms estimate of
mechanism 5. This sets the sample step of stage A and tells whether a comms
receiver sees frozen fades.

### Stage C: combine with the fidelity-2 retro record

For each trial of a retro campaign, multiply the power of the existing
record by an independent draw of `|AF|^2 / mean`. The target speckle and the
turbulence come from different physics, so the draws are independent for
SNAPSHOT statistics. Report the combined p1 / p5 / mean against the
atmosphere alone. A JOINT time series needs the retro time axis, which the
runner refuses today (a temporal retro is open in the CLAUDE.md retro item).

### Stage D (only if needed): separated arrays

For a spacecraft with arrays tens of metres apart, one shifted return window
for each array group in the same trial. Not needed for LAGEOS.

## Open decisions (owner)

1. The priority: CW (comms or beacon, full speckle) or pulsed (ranging,
   target signature).
2. The cube far-field model: the first cut (area times Airy) or the real
   spoiled six-sector pupil.
3. The cube catalogue and the spin state of LAGEOS-1: which source.
4. Where the array factor lives in olb if it becomes a Term. The
   `LidarCrossSection(sigma_m2)` retroreflector (the mean only, fidelity 0
   and 1) is BUILT (2026-10-01, `olb.links.retro_space.retro_cross_section_term`).

## References to check

These DOIs come from memory and are NOT yet checked against the papers:
Goodman 1976 (10.1364/JOSA.66.001145), Fried 1982 (10.1364/JOSA.72.000052),
Otsubo and Appleby 2003 (10.1029/2002JB002209). The Degnan, Shapiro and Born
and Wolf DOIs are the ones that the package already uses (the Degnan and the
Born and Wolf DOIs are also unchecked, see the CLAUDE.md retro item). The
ILRS LAGEOS cross section (about 7e6 m^2) and the 1.25 arcsec dihedral
offset are also from memory.

## The files

| File | Purpose |
| --- | --- |
| [lageos_pass.py](lageos_pass.py) | The 2026-10-01 LAGEOS-1 pass over Yarragadee: the match to the mount readout, the point-ahead angle, and the single unspoiled cube Airy loss at 532 and 1064 nm. It needs `skyfield` only. |
