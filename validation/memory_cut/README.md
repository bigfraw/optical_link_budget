# The memory cut of a fidelity-2 trial, and the two speed opt-ins

Date: 2026-09-06. Machine: a 4-core Linux cloud container with 15 GB (NOT
bigfraw, NOT the laptop). The wall times below are for that box only; the
ratios carry over, the seconds do not.

Two scripts. Run each one from the repository root on a QUIET machine.

```
python -m validation.memory_cut.memory_cut_check
python -m validation.memory_cut.screen_generator_lean
```

Each one writes a JSON of its numbers next to itself.

## 1. `memory_cut_check.py`: the bit-identical cut

Three claims of the 2026-09-06 change, checked mechanically:

1. **Lazy screens.** `split_step` reads its screens one at a time. A
   generator and a list of the same screens give the same field, bit for bit,
   and a wrong count raises.
2. **The Forvard cache.** A warm cache gives the same field as a cold one, and
   a trial with the cache OFF (`FORVARD_CACHE_BYTES = 0`, the old behaviour)
   gives the same collected power and SMF eta as a trial with it ON, bit for
   bit.
3. **The pool sizer.** `worker_memory_bytes` sits above the measured peak of
   one serial trial, and `auto_workers` obeys both limits.

The measurement, one serial trial of the 30 deg hero downlink (0.7 m SMF
ground terminal, 500 km, L0 = 25 m), pinned at 1024 px:

| precision | screens | cache | s/trial | peak MiB | sizer estimate MiB |
|---|---|---|---|---|---|
| single | 9 | off | 3.61 | 156 | |
| single | 9 | on | 2.58 | 156 | 410 |
| double | 15 | off | 8.09 | 200 | |
| double | 15 | on | 6.64 | 200 | 740 |

The cache cuts the time by 18 to 29 percent. The peak does not move with the
cache: it is set by the FFT work arrays. The lazy screens set the peak, and
the earlier session measurement against the pre-change code read 321 to
201 MiB (double, 15 screens) and 193 to 157 MiB (single, 9 screens).

## 2. `screen_generator_lean.py`: the two opt-ins

Neither one is bit-identical to the default of record, so neither one is a
default. Both enter a campaign fingerprint, so a stored campaign never mixes
with them.

**The lean screen body** (`screen_generator="olb-lean"`,
`ScreenFactory(lean=True)`). The same physics, the same random stream, one
third fewer full-grid passes: the normals go straight into the complex grid,
the filter multiply and the inverse transform run in place through scipy.fft,
and the two centring shifts become sign flips (the double-shift identity
Forvard uses).

| check | result |
|---|---|
| draw agreement, same seed, float64 | 2e-16 to 4e-16 relative rms |
| draw agreement, same seed, float32 | 1e-7 to 2e-7 relative rms |
| structure-function r0 fit, 40 screens, 1024 px float32 | default 0.1111 +- 0.0010 m, lean 0.1111 +- 0.0010 m |

| grid | dtype | default ms / MiB | lean ms / MiB |
|---|---|---|---|
| 1024 | float32 | 86 / 48 | 60 / 28 |
| 1024 | float64 | 137 / 80 | 87 / 56 |
| 2048 | float32 | 448 / 192 | 299 / 112 |
| 2048 | float64 | 770 / 320 | 588 / 224 |

**The scipy FFT backend** (`fft_backend="scipy"`,
`olb.waveoptics.propagators.set_fft_backend`). `Forvard` runs its two
transforms through `scipy.fft` with `overwrite_x=True`. The raw 1024 px
complex64 `fft2` on this box: numpy 69.5 ms, scipy 15.6 ms, scipy in place
13.0 ms. The numpy transform of numpy 2.4 is not vectorised for complex64
here; scipy's pocketfft is.

**The whole trial**, single precision, 1024 px, 9 screens, best of three:

| backend | generator | s/trial | gain | power / eta from the default |
|---|---|---|---|---|
| numpy | olb (the default) | 2.74 | | 0 / 0 |
| scipy | olb | 2.19 | 1.25x | 5.6e-7 / 6.1e-7 |
| numpy | olb-lean | 2.57 | 1.07x | 1.1e-7 / 1.5e-7 |
| scipy | olb-lean | 1.99 | 1.38x | 5.6e-7 / 2.3e-7 |

The agreement is the rounding of single precision, the same level the
double-to-single switch measured (`validation/precision/`).

## What the numbers say

- The Forvard cache and the lazy screens are free: bit-identical, on by
  default.
- The scipy backend is the larger of the two opt-ins, and it also moves
  fewer bytes (the transform writes into its input). The lean generator adds
  on top. Together they give 1.38x on one core, on top of the 1.3x of the
  cache.
- On a memory-bandwidth-bound pool the gain in trials per second can be
  larger than the one-core gain, because both opt-ins cut the bytes each
  trial moves. That is NOT measured: measure it on bigfraw with
  `validation/campaign_resources/ --workers auto` when the box is free, once
  with the defaults and once with `fft_backend="scipy"` and
  `screen_generator="olb-lean"` (the script does not take those flags yet).
- Whether either opt-in becomes a DEFAULT is an owner decision, because each
  one changes every seeded fidelity-2 number at the rounding level.
