# The single-against-double precision check of the fidelity-2 layer

`precision_check.py` runs the SAME turbulent trials two times, one time in
double precision and one time in single precision, and it measures the
difference.

- Single precision: a complex64 field, float32 phase screens. It is the
  DEFAULT of every runner since 2026-09-05 (owner decision):
  `precision="single"` on `propagate_turbulent_scenario`, `Campaign`,
  `run_waveoptics` and `run_fidelity2`.
- Double precision: a complex128 field, float64 phase screens. It is
  `precision="double"`, and it reproduces every run made before that date bit
  for bit.

The seed is the same in the two runs, so trial k sees the same atmosphere.

## Why single precision exists

A `Campaign` on the desktop is memory-bandwidth bound: twelve processes
saturate the memory channels on a 512 px grid. Half the bytes for each element
is the one change with real leverage there. `numpy.fft` keeps complex64, so no
FFT library changes.

## The case

It is the case of `validation/campaign_resources`: a 700 mm ground aperture
with a 30 percent central obscuration, an SMF detector, 1550 nm, a 500 km orbit
at 30 deg, and the fixed outer scale L0 = 25 m. The preset is `rapid`, so the
check is short.

## Run it

```
python -m validation.precision.precision_check
python -m validation.precision.precision_check --n-trials 100
python -m validation.precision.precision_check --preset standard
```

It writes `precision_check.json` next to itself.

## Read the result

The script prints, and the JSON holds:

- `max_d_collected_power`: the largest relative difference of the collected
  power over the trials.
- `max_d_smf_eta`: the same for the fibre coupling efficiency.
- `max_d_field_rms`: the largest rms difference of the receive-plane field,
  over the peak amplitude of the double-precision field.
- `speed_ratio_double_over_single`: the wall time of the double run over the
  wall time of the single run, on the threaded runner.

EXPECT about 1e-5 or better on the three difference columns. A number worse
than 1e-3 means a kernel lost precision, and the script prints a warning and
exits 1. Look at the phase wrap of `olb.waveoptics.propagators.Forvard` first:
it keeps `Bus` in double precision on purpose, because that subtraction takes
the fractional part of a number in the thousands.

## The caution

A single-precision campaign is a DIFFERENT record. Its trials are not
bit-identical to a double-precision run of the same seed, and the arithmetic
carries about 7 digits, not 16. `precision` enters the campaign fingerprint, so
a single-precision store never mixes with a double-precision store. Run this
check before a budget reads a single-precision campaign.

## Measured, 2026-09-05

8 trials, `rapid` preset, 30 deg, a 256 px grid of 2.638 m, 5 screens, on the
laptop:

| quantity | value |
|----------|-------|
| max relative difference, collected power | 4.6e-07 |
| max relative difference, SMF eta | 3.4e-06 |
| max relative rms, receive field | 1.3e-06 |
| wall time, double | 0.41 s |
| wall time, single | 0.31 s |
| speed ratio (double / single) | 1.36x |

So single precision changes the physics by parts in a million, far below the
Monte Carlo spread of a fade statistic. The wall-time ratio of this small case
is noisy; measure it again on the target grid and the target machine. On the
desktop at 12 workers and 512 px the gain is 1.32x (11.2 against 8.5
trials/s, `validation/campaign_resources/`). A `precision="double"` run is
unchanged, bit for bit.

20 trials, `standard` preset, 30 deg, a 512 px grid of 3.514 m, 9 screens, on
the laptop (`precision_check_standard.json`):

| quantity | value |
|----------|-------|
| max relative difference, collected power | 5.6e-07 |
| max relative difference, SMF eta | 2.4e-06 |
| max relative rms, receive field | 1.5e-06 |
| wall time, double | 5.35 s |
| wall time, single | 4.02 s |
| speed ratio (double / single) | 1.33x |

The agreement holds at the production grid and screen count.
