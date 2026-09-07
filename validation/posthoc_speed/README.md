# The post-hoc read of a stored campaign: the CROP

Date: 2026-09-07. Script: `posthoc_speed.py`. Record: `posthoc_speed.json`.

## The problem

A stored fidelity-2 trial holds the pixels of the patch DISC only. The old
read-back scattered those pixels into the FULL grid, and every later step then
swept the zero padding, which is most of the grid. At 1024 px the patch of the
0.7 m hero downlink is a 271 px square, so the padding is 14.3 times the useful
pixels. The old read also rebuilt the aperture mask, the fibre mode, the modal
basis and the slope reconstructor for EACH trial.

## The change

The read now works on the square CROP that just holds the disc, and it builds
the clip mask, the fibre mode, the modal basis and the slope reconstructor ONE
time for a call.

**THE CROP RULE: PUPIL-plane quantities on the crop, FOCAL-plane quantities on
the padded grid.** The crop keeps the pixel pitch and the centre pixel of the
grid, so the clip, the collected power, the single-mode overlap and the modal
fit read the SAME pixels. An `MMF` or a `Camera` FOCUSES the field, and the
focal-plane pixel scale reads the grid EXTENT, so those pad the crop back.

## The case

The 30 deg hero downlink, 500 km, a 0.7 m ground aperture with a 0.3
obscuration, an SMF receiver, `rapid` preset, `L0 = 25 m`, single precision,
seed 20260907. The grid is PINNED at 1024 px (2.638 m, 5 screens). The store
holds 256 trials in blocks of 8, with the summed screen phase. The machine is
the laptop, and the page cache is warm.

## 1. The agreement, crop against full grid

The worst relative difference over the 256 trials:

| Quantity | Worst relative difference |
|---|---|
| `recollect` | 4.4e-16 |
| `recouple`, SMF | 3.1e-15 |
| `recouple`, MMF | 0 (bit-identical) |
| `recouple_compensated`, `source="screens"` | 2.9e-15 |
| `recouple_compensated`, `source="slopes"` | 1.1e-15 |

The pupil quantities agree at the float rounding level: the crop and the full
grid hold the same pixels, and a sum over the crop adds in a different ORDER
only. The MMF route pads the crop back to the full grid before it focuses, so
it is bit-identical.

## 2. The time of one trial, one core

| Read | Full grid (s) | Crop (s) | Speed-up |
|---|---|---|---|
| `recollect` | 0.017 | 0.001 | 11.7x |
| `recouple`, SMF | 0.035 | 0.014 | 2.5x |
| `recouple_compensated`, `source="screens"` | 0.084 | 0.018 | 4.7x |
| `recouple_compensated`, `source="slopes"` | 0.565 | 0.048 | 11.9x |

The slope route wins the most, because it pays a full-grid wrapped gradient, a
modal fit and a modal apply for each trial, and it builds a slope reconstructor
one time. The SMF row wins the least, because the crop read is now bound by the
npz block read and not by the arithmetic.

## 3. The block-parallel `map_trials`, 4 workers

The whole 256-trial read, on the crop:

| Read | Serial (s) | Pool of 4 (s) | Speed-up |
|---|---|---|---|
| `recouple`, SMF | 3.64 | 6.20 | 0.59x |
| `recouple_compensated`, `source="slopes"` | 12.22 | 9.51 | 1.28x |

The pool gives the SAME numbers, value for value.

**THE HONEST READING: a process pool does not pay for a light read.** The
Windows pool spawn costs 2.5 to 4.4 s (the same cost the 2026-09-04 fair rerun
measured for `Campaign.run`), and each worker must import olb again. The SMF
read of 256 trials is 3.6 s of work, so the spawn is more than the work. The
compensated slope read is 12 s of work, and the pool wins 1.28x on this laptop.

So use `workers=` for a HEAVY read of a LARGE campaign (a compensated read of
thousands of trials on a many-core machine), and leave the default `None` for
everything else. The crop is the change that always pays.
