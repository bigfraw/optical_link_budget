# The retro bracket: do the two legs of a retro link fade together?

Date: 2026-09-29. Branch `retro-velocity-aberration`. Machine: the laptop
(`bigfraw` was offline), on the CPU with `fft_backend="scipy"`.

## The purpose

A fidelity-2 SPACE retro link has an UP LEG (ground to the corner cube) and a
RETURN LEG (corner cube to ground). The cheapest fidelity-2 route pairs the two
legs from DIFFERENT trials of one downlink campaign, so the legs see
INDEPENDENT atmospheres and no point-ahead record is needed. This study asks
if that route is safe. It measures the retro fade for three pairings of the
two legs on the SAME stored trials, with no new propagation for the read:

- SAME: the up leg and the return leg of trial k read the same atmosphere.
- PA: the up leg reads the point-ahead window of trial k and the return leg
  reads the beacon window of trial k (angles 0, geometric, 2, 5 and 10 arcsec).
- INDEP: the legs read different atmospheres. The retro fade in dB is the sum
  of the two leg fades, so the script uses all n^2 pairs of the two measured
  distributions.

## The legs

- UP LEG: the Shapiro reciprocity overlap of the stored ground field with the
  ground transmit mode, over its own vacuum overlap
  (`Campaign.uplink_overlaps`; DOI 10.1364/JOSA.61.000492). A small corner cube
  is a point receiver in the uplink far field, so the overlap IS the power it
  catches. Two launches on the 0.7 m ground aperture: the hero full-aperture
  waist (0.35 m, monostatic, near the fibre mode) and a small 0.06 m waist.
- RETURN LEG: the beacon-window field at the ground, into a 0.7 m single-mode
  fibre (`Campaign.recouple(SMF())`) and a 0.7 m bucket (`Campaign.recollect`).

## The case

The hero uplink of `validation/anisoplanatism_screens/common.py`
(`hero_uplink`): 1550 nm, 500 km, site `cn2_ground = 1.7e-14`,
`wind_rms = 21 m/s`, `L0 = 25 m`, seed 20260907, preset `standard`, single
precision, no correction. The campaigns are NEW point-ahead `base` campaigns
under `campaigns/` (gitignored), built by `waveoptics_pointahead.campaign_of`
with block size 25. They are NOT the bigfraw campaigns (the key holds the FFT
backend and the block size).

## The run line

From the repository root:

```
$env:OLB_WAVEOPTICS_POINTAHEAD_ROOT = "validation/retro_bracket/campaigns"
python -m validation.retro_bracket.retro_bracket --elevations 30 20 `
    --n-trials 1000 --block-size 25 --fft-backend scipy --preset standard `
    --n-boot 200 --run-workers 4
```

`--run-workers` computes the missing trials first; leave it out to read only.
The run took 467 s (30 deg) and 368 s (20 deg) for the trials with 4 workers,
plus the read and 200 bootstraps. The script writes
`retro_bracket_el30_el20.log` and `.json` (both gitignored), so the table
below is the record.

## The result

Pairing minus INDEP, in dB (positive = more loss), with the paired-bootstrap
standard error. `rho` is the correlation of the two leg losses in dB. 1000
trials for each elevation.

| launch, return | el | pairing | rho | mean | p5 | p1 |
| --- | --- | --- | --- | --- | --- | --- |
| hero 0.35 m, SMF | 30 | SAME | +0.99 | -3.84 +/- 0.14 | +7.22 +/- 0.83 | +11.81 +/- 2.69 |
| | 30 | PA 5.24" (geom) | +0.56 | -2.92 +/- 0.24 | +4.05 +/- 0.70 | +7.09 +/- 2.01 |
| | 30 | PA 10" | +0.44 | -2.33 +/- 0.27 | +3.63 +/- 0.83 | +3.68 +/- 1.66 |
| | 20 | SAME | +0.99 | -3.86 +/- 0.17 | +8.77 +/- 0.94 | +9.78 +/- 1.93 |
| | 20 | PA 3.58" (geom) | +0.45 | -2.26 +/- 0.25 | +3.11 +/- 0.68 | +3.59 +/- 1.59 |
| | 20 | PA 10" | +0.34 | -2.04 +/- 0.28 | +2.48 +/- 0.73 | +1.20 +/- 1.57 |
| small 0.06 m, SMF | 30 | SAME | +0.39 | -1.27 +/- 0.07 | +2.80 +/- 0.67 | +4.49 +/- 1.32 |
| | 30 | PA 5.24" (geom) | +0.36 | -1.05 +/- 0.10 | +2.51 +/- 0.83 | +2.93 +/- 3.04 |
| | 20 | SAME | +0.36 | -1.54 +/- 0.14 | +2.46 +/- 0.62 | +2.63 +/- 2.14 |
| | 20 | PA 3.58" (geom) | +0.27 | -1.29 +/- 0.19 | +2.99 +/- 0.64 | +2.22 +/- 2.04 |
| hero 0.35 m, bucket | 30 | SAME | +0.05 | -0.01 +/- 0.01 | -0.24 +/- 0.19 | -0.21 +/- 0.30 |
| | 20 | SAME | +0.08 | -0.04 +/- 0.02 | -0.05 +/- 0.20 | -0.15 +/- 0.39 |
| small 0.06 m, bucket | 30 | SAME | +0.03 | -0.02 +/- 0.01 | -0.31 +/- 0.18 | -0.06 +/- 0.32 |
| | 20 | SAME | +0.09 | -0.06 +/- 0.02 | +0.26 +/- 0.20 | +0.43 +/- 0.48 |

The INDEP reference (absolute, dB; mean / p5 / p1): hero SMF 27.21 / 49.85 /
58.78 (30 deg) and 31.46 / 54.04 / 62.91 (20 deg); hero bucket 13.26 / 28.26
/ 35.95 and 15.40 / 30.69 / 36.49. The bucket values have no vacuum reference
(they are normalised to the mean), so only their differences carry meaning.
Every bucket PA column also sits inside about 1.5 SE of zero. The other PA
angles are in the log.

THE IDENTITY: the PA 0" column equals the SAME column on every row (the script
asserts that angle 0 is the beacon window).

## The verdict

1. BUCKET RETURN: the legs do not correlate (rho about 0 in every pairing), so
   the INDEP route is SAFE.
2. SMF RETURN: the legs fade TOGETHER. INDEP under-reads the p5 fade by 2.5 to
   9 dB, depending on the pairing. The mean moves the other way (1 to 4 dB
   LESS loss), so a link sized on the mean hides the effect. The hero launch
   correlates more than the small launch, because its mode is near the fibre
   mode (reciprocity: SAME gives rho = 0.99).

## THE GEOMETRY CAVEAT (found after the run)

The PA pairing is probably the WRONG model of a retro link. One piece of light
goes up along the line to the satellite position P1 at the time it arrives,
and the return that reaches the station comes back from P1 along the SAME
line. The point-ahead lead and the turn of the line of sight during the flight
are the same angle (2 v_perp / c), so they cancel. The two passes of one pulse
share one line; only the WIND moves the air between them, over the round trip
dt = 2R/c (about 6.1 ms at 30 deg and 8.0 ms at 20 deg for 500 km). The PA
pairing puts the up leg of one pulse against the return of a pulse sent 2R/c
EARLIER, which is the pre-compensated-uplink geometry (2-P4), not the retro.

So the SAME column is the dt -> 0 limit of the right geometry, and the owner
agreed the argument (2026-09-30). The runner now models it (`direction=
"retro"`, below).

## The retro run (2026-09-30): the right geometry, measured

The runner takes `direction="retro"`
(`olb.waveoptics.turbulence.run.retro_sensing_geometry`). Each trial draws ONE
oversize atmosphere: the RETURN leg reads the unshifted window, and the UP leg
reads a window shifted at each screen by the Bufton wind drift over the round
trip, `|V(h)| * 2R/c` (Vg = 10 m/s, one wind direction, integer pixels). The
margin is 24.2 cm at 30 deg and 31.7 cm at 20 deg (the screen grows from 512
to 576 px). The script mode `--retro` builds these campaigns and compares the
RETRO pairing (the same trial) with INDEP (all n^2 pairs of the same record):

```
python -m validation.retro_bracket.retro_bracket --retro --elevations 30 20 `
    --n-trials 1000 --block-size 25 --fft-backend scipy --preset standard `
    --n-boot 200 --run-workers 4
```

The run took 150 s (30 deg) and 170 s (20 deg) for the trials with 4 workers.
The hero launch (0.35 m waist) only, 1000 trials for each elevation. RETRO
minus INDEP, in dB:

| return | el | rho | mean | p5 | p1 |
| --- | --- | --- | --- | --- | --- |
| SMF | 30 | +0.73 | -3.54 +/- 0.16 | +6.48 +/- 0.49 | +6.12 +/- 0.97 |
| SMF | 20 | +0.59 | -3.36 +/- 0.19 | +2.85 +/- 0.48 | +3.36 +/- 1.96 |
| bucket | 30 | +0.01 | -0.01 +/- 0.01 | -0.06 +/- 0.16 | -0.51 +/- 0.27 |
| bucket | 20 | +0.08 | -0.05 +/- 0.02 | +0.06 +/- 0.21 | +0.01 +/- 0.29 |

INDEP reference (mean / p5 / p1, dB): SMF 27.29 / 50.47 / 59.57 (30 deg) and
31.42 / 52.09 / 59.43 (20 deg).

THE READING:

1. BUCKET: no correlation, as in every pairing above. Independent legs are
   safe for a bucket return.
2. SMF: the wind decorrelates the legs only in part. At 30 deg (dt = 6.1 ms)
   the RETRO p5 penalty is +6.5 dB, near SAME (+7.2) and above the PA pairing
   (+4.1). At 20 deg (dt = 8.0 ms, a larger drift) it is +2.9 dB, near the PA
   pairing (+3.1) and far below SAME (+8.8). So the drift over the round trip
   matters, and no fixed pairing stands in for it. The fidelity-2 retro budget
   reads this record.
3. The mean moves the other way (3.4 to 3.5 dB LESS loss than INDEP), so a link
   sized on the mean hides the fade.

## Open

- A tip-tilt receive loop on the return leg: it removes part of the shared
  tilt, so the SMF gap can shrink.
- The sensitivity to Vg and to the wind direction of each layer (the model
  takes one direction for all layers).
- The small launch (0.06 m) on a retro campaign (the up leg of a retro record
  reads the scenario launch only).
- A rerun on bigfraw at more trials.
