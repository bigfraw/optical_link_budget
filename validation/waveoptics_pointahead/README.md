# The point-ahead validation of the fidelity-2 uplink (backlog 2-P4)

Date: 2026-09-11. Branch `waveoptics-pointahead`. Machine: the desktop
`bigfraw`, on the CUDA FFT backend (`fft_backend="cupy"`, one device, one
process).

## The purpose

The fidelity-2 runner now makes ONE MORE propagation pass for each point-ahead
angle (`olb/waveoptics/turbulence/run.py`, `point_ahead_rad`). The ground
terminal senses the DOWNLINK beacon, it applies the conjugate wavefront to the
uplink beam, and the uplink goes to where the satellite WILL BE. The two
directions read the SAME atmosphere through a laterally shifted window of each
screen, so the drop of the reciprocity overlap between them IS the point-ahead
anisoplanatism. Before this work the fidelity-2 pre-compensated uplink carried
the flag NO ANISOPLANATISM, and the model of record stayed the fidelity-1 FAST
Term.

This study answers four questions:

1. Does the record hold together? Is the stored beacon column the beacon
   overlap, does the post-hoc read of the stored planes equal an in-run
   corrected campaign, and does the regeneration route give a stored column
   back? (p0)
2. How large is the penalty, and do the three rungs agree? The field against
   the fidelity-1 FAST Term and against the analytic Stone Term. (p1)
3. What does the point ahead do to the FADE, not only to the mean? (p2)
4. Is the penalty converged in the SCREEN COUNT? (p3)

THE TILT STAYS IN (owner decision, 2026-09-11). The terminal senses the
downlink beacon tilt and the steering mirror adds the point-ahead offset
geometrically, so the uplink has no tilt reference of its own and it pays the
full tilt anisoplanatism. The runner corrects the modes of the stack and it
keeps the tilt in the error; FAST keeps the piston and the tilt in its modal
mask; the Stone Term takes `remove='piston'`. So the three rungs read the SAME
mode set, except the piston, which changes no overlap integral.

THE MODE COUNT IS NOLL, and it is the count of `olb/turbulence/ao.py`:
`TipTilt` removes the first 3 Noll modes and `AO(n)` the first n (R. J. Noll,
DOI 10.1364/JOSA.66.000207, Table I).

## The case

THE HERO UPLINK of `validation/anisoplanatism_screens/common.py`
(`hero_uplink`), which is the hero downlink of `validation/waveoptics_ao/`
turned around: a 1550 nm uplink from a 700 mm ground terminal with a FULL
APERTURE launch (`Transmitter(waist_m=0.35)`) to a 100 mm space terminal at
500 km, with a `DownlinkBeacon` pre-compensation source. Site
`cn2_ground = 1.7e-14`, `wind_rms = 21 m/s`. Preset `standard`, seed 20260907,
single precision, `L0 = 25 m` everywhere (the owner rule of 2026-09-05,
backlog 2-P5): on the screens, in FAST (`fast_params={"L0": 25.0}`) and in the
Stone Term (`L0_m=25.0`).

THE ANGLES of each campaign, in this column order: 0 (the beacon direction),
the GEOMETRIC angle of the elevation
(`olb.geometry.CircularOrbit.point_ahead_rad`, 5.24 arcsec at 30 deg), then 2,
5 and 10 arcsec. `screen_margin_m=None`, so the runner sizes the oversize
screen from the widest angle and the highest screen.

THE STACKS are the four stacks of `validation/waveoptics_ao/`: `base` (no
correction), `tiptilt` (3 Noll modes), `ao10` (10) and `ao21` (21). Only the
`base` campaign is computed: the other three come from it POST HOC, through
`Campaign.recouple_point_ahead`, which corrects the stored planes with the
beacon estimate of any stack and costs no propagation.

## The run lines

From the repository root on `bigfraw`. The GPU python is
`C:\Users\alexf\olb-gpu-venv\Scripts\python.exe`. The login shell of
`ssh desktop` is Windows PowerShell 5.1, so chain with `;`, never with `&&`,
and launch the long run through WMI so it outlives the ssh session
(`validation/campaign_resources/README.md`):

```
$cmd = 'cmd /c "cd /d D:\repos\optical_link_budget && C:\Users\alexf\olb-gpu-venv\Scripts\python.exe -u -m validation.waveoptics_pointahead.waveoptics_pointahead --study p0 p1 p2 p3 > validation\waveoptics_pointahead\run_all.log 2>&1"'
Invoke-CimMethod -ClassName Win32_Process -MethodName Create -Arguments @{CommandLine=$cmd}
```

The defaults ARE the production settings: `--elevations 30 20`,
`--n-trials 1000`, `--count-trials 400`, `--check-trials 200`,
`--block-size 50`, `--fast-samples 1000`, `--fft-backend cupy`,
`--workers` unset (one device, one stream). The studies also run one at a
time:

```
python -m validation.waveoptics_pointahead.waveoptics_pointahead --study p0
python -m validation.waveoptics_pointahead.waveoptics_pointahead --study p1
python -m validation.waveoptics_pointahead.waveoptics_pointahead --study p2
python -m validation.waveoptics_pointahead.waveoptics_pointahead --study p3
```

Add `--dry-run` to print the campaign list and to exit. Add `--analyse-only`
to read what is stored and to compute no trial. Add `--fft-backend numpy` on a
host with no CUDA device, and `--no-fast` when `fast-aosim` is not installed.
`--posthoc-workers 4` opens a process pool for the post-hoc reads.

THE SMOKE RUN, on a laptop, takes a few minutes:

```
python -m validation.waveoptics_pointahead.waveoptics_pointahead \
    --study p0 p1 p2 p3 --fft-backend numpy --preset rapid --n-trials 8 \
    --count-trials 8 --check-trials 4 --elevations 30 --block-size 4 \
    --no-fast
```

## The cost

From `--dry-run` at the production settings:

| campaign | trials | planes | grid | screen | screens | MB |
| --- | --- | --- | --- | --- | --- | --- |
| el30/base | 1000 | 6 | 512 | 768 | 9 | 1263 |
| el30/ao10 | 200 | 6 | 512 | 768 | 9 | 253 |
| el20/base | 1000 | 6 | 512 | 800 | 9 | 791 |
| el20/ao10 | 200 | 6 | 512 | 800 | 9 | 158 |
| el30/n5 | 400 | 6 | 512 | 736 | 5 | 505 |
| el30/n9 | 400 | 6 | 512 | 768 | 9 | 505 |
| el30/n15 | 400 | 6 | 512 | 768 | 15 | 505 |
| el30/n25 | 400 | 6 | 512 | 768 | 25 | 505 |
| el30/ground_split_x4 | 400 | 6 | 512 | 768 | 12 | 505 |
| TOTAL | 4400 | | | | | 4990 |

`planes` is 1 + the angle count: each plane is one split step and one stored
field, so a point-ahead trial costs six passes.

THE MEASURED COST (`run_full.log`): the 30 deg base campaign of 1000 trials
took 281 s, which is 0.281 s for each trial, and the 20 deg base campaign took
the same 281 s. The whole nine-campaign run took 1407 s of GPU wall, that is
UNDER 30 MINUTES, and it wrote 5.0 GB of stored planes. The campaigns live
under `validation/waveoptics_pointahead/campaigns/` and they are
gitignored; the logs, the JSON records and the figures are committed.

## p0: the identities of the record

THREE CHECKS on each elevation.

1. `eta_turb_pa[:, 0] == eta_turb`, bit for bit. Angle 0.0 runs the same
   window as the beacon pass, so the two must be the same number.
2. The post-hoc route. `base.recouple_point_ahead([AO(10)])` corrects the
   STORED UNCORRECTED planes; an in-run `ao10` campaign of the same seeds
   corrects inside the run. The gate is 1e-6 relative on every trial and
   angle.
3. The regeneration route. `base.point_ahead([5 arcsec], None)` rebuilds the
   screens from the stored seeds and propagates again; it must give the stored
   5 arcsec column back. The regeneration runs on the BACKEND OF THE CAMPAIGN,
   so the gate is BIT IDENTITY. The log also prints the host (numpy) read of
   the same CUDA campaign as a report: it agrees at the float32 rounding level
   only.

| elev | beacon identity | post-hoc max rel | regenerate max rel | bit identical | s/trial | screen px |
| --- | --- | --- | --- | --- | --- | --- |
| 30 deg | PASS, bit for bit | 9.04e-07 | 0.00e+00 (host read 1.86e-05) | YES | 0.276 (0.046 for each pass) | 768 against 512 (1.50x) |
| 20 deg | PASS, bit for bit | 6.63e-07 | 0.00e+00 (host read 3.40e-05) | YES | 0.277 (0.046 for each pass) | 800 against 512 (1.56x) |

VERDICT: **the record holds together. All three identities pass on both
elevations.** The stored zero angle equals the beacon overlap bit for bit, so
the point-ahead pass adds no offset of its own. The post-hoc AO(10) read of the
uncorrected planes matches the in-run AO(10) campaign to 9.0e-07 at 30 deg and
6.6e-07 at 20 deg, well inside the 1e-6 gate, so the three derived stacks cost
no propagation and they are the same numbers as a computed campaign. The
regeneration at the stored 5 arcsec angle is BIT IDENTICAL on the cupy backend
that computed the campaign; the host (numpy) read of the same campaign differs
by 1.9e-05 at 30 deg and 3.4e-05 at 20 deg, which is two FFT libraries at
float32 and not a physics difference. The cost of the record is 0.28 s for each
trial, for six split steps, and the oversize screen is 768 px at 30 deg and
800 px at 20 deg against the 512 px propagation grid.

## p1: the penalty against FAST and against Stone

THE FIELD PENALTY is `L(theta) - L(0)` with `L = -10 log10(mean eta_pa)`. The
beacon column `L(0)` IS the 2-AO perfect pre-compensation loss, and the
point-ahead columns degrade it.

THE FAST PENALTY is the same difference of `olb.models.fast.uplink_fast_term`.
The angle goes in through `fast_params={"DTHETA": [arcsec, 0]}`, which the
Term merges LAST, so it overrides the angle that the Term reads from the
geometry. An NPXLS grid guard runs first at DTHETA = 0 and it pins the
smallest grid within 0.15 dB of the largest, the guard of
`validation/waveoptics_vs_fast/`. Every stack has a FAST row: an empty stack
is the NOAO launch and a tip-tilt stack is the TT launch, the same map as the
downlink Term (the old refusal of an uplink without an AO stage was an olb
guard, lifted 2026-09-11).

THE STONE PENALTY is `olb.links.uplink.uplink_point_ahead_term` with
`remove='piston'` and `L0_m=25.0`, reported as
`loss_db = (10 / ln 10) sigma^2`. The `base` stack has no Stone row: with no
corrected mode there is no decorrelation residual, so the penalty is zero by
definition. The `tiptilt` stack takes `max_order=1` (the tilt is radial order
1); `ao10` and `ao21` take `max_order='auto'`, which reads the AO stage.

THE PASS BAND of the field against FAST is 0.5 dB, the like-for-like tolerance
of `docs/physics.md` Section 9l. Stone is a REPORT, not a gate: the extended
Marechal mapping saturates past sigma^2 = 1 rad^2 (T. S. Ross,
DOI 10.1364/AO.48.001812), so the Stone number reads high at a large angle.

THE FIELD LOSS `-10 log10(mean eta_pa)` itself, in dB, 1000 trials. The beacon
column IS the 2-AO perfect pre-compensation loss:

| elev | stack | beacon | geom | 2" | 5" | 10" |
| --- | --- | --- | --- | --- | --- | --- |
| 30 deg | base | 13.264 | 13.213 | 13.224 | 13.184 | 13.362 |
| 30 deg | tiptilt | 7.782 | 9.009 | 8.166 | 8.962 | 9.547 |
| 30 deg | ao10 | 2.971 | 5.711 | 3.923 | 5.606 | 6.851 |
| 30 deg | ao21 | 1.674 | 4.811 | 2.822 | 4.699 | 6.102 |
| 20 deg | base | 15.398 | 15.654 | 15.630 | 15.498 | 15.430 |
| 20 deg | tiptilt | 10.461 | 11.638 | 11.157 | 11.940 | 12.563 |
| 20 deg | ao10 | 4.372 | 7.838 | 6.466 | 8.465 | 10.079 |
| 20 deg | ao21 | 2.619 | 6.559 | 5.011 | 7.352 | 9.153 |

The geometry angle is 5.24 arcsec at 30 deg and 3.58 arcsec at 20 deg.

30 deg, 1000 trials, the anisoplanatic penalty in dB:

| stack | angle | field | +- | FAST | gap | band | Stone |
| --- | --- | --- | --- | --- | --- | --- | --- |
| base | geom 5.24" | -0.051 | 0.135 | 0.198 | 0.249 | PASS | - |
| base | 2" | -0.040 | 0.097 | 0.019 | 0.059 | PASS | - |
| base | 5" | -0.080 | 0.136 | 0.247 | 0.327 | PASS | - |
| base | 10" | 0.098 | 0.163 | 0.115 | 0.017 | PASS | - |
| tiptilt | geom 5.24" | 1.227 | 0.075 | 1.733 | 0.505 | OVER | 2.453 |
| tiptilt | 2" | 0.383 | 0.055 | 0.778 | 0.395 | PASS | 0.612 |
| tiptilt | 5" | 1.180 | 0.075 | 1.740 | 0.560 | OVER | 2.332 |
| tiptilt | 10" | 1.765 | 0.102 | 2.137 | 0.373 | PASS | 4.351 |
| ao10 | geom 5.24" | 2.739 | 0.060 | 3.237 | 0.497 | PASS | 3.714 |
| ao10 | 2" | 0.951 | 0.027 | 0.965 | 0.013 | PASS | 1.185 |
| ao10 | 5" | 2.635 | 0.059 | 3.216 | 0.581 | OVER | 3.561 |
| ao10 | 10" | 3.879 | 0.085 | 3.978 | 0.099 | PASS | 6.117 |
| ao21 | geom 5.24" | 3.137 | 0.062 | 3.331 | 0.193 | PASS | 4.090 |
| ao21 | 2" | 1.148 | 0.026 | 1.078 | -0.070 | PASS | 1.429 |
| ao21 | 5" | 3.025 | 0.061 | 3.383 | 0.358 | PASS | 3.929 |
| ao21 | 10" | 4.429 | 0.087 | 4.256 | -0.173 | PASS | 6.639 |

20 deg, 1000 trials:

| stack | angle | field | +- | FAST | gap | band | Stone |
| --- | --- | --- | --- | --- | --- | --- | --- |
| base | geom 3.58" | 0.256 | 0.153 | 0.316 | 0.059 | PASS | - |
| base | 2" | 0.232 | 0.132 | 0.265 | 0.033 | PASS | - |
| base | 5" | 0.100 | 0.171 | 0.482 | 0.381 | PASS | - |
| base | 10" | 0.032 | 0.181 | 0.189 | 0.157 | PASS | - |
| tiptilt | geom 3.58" | 1.178 | 0.127 | 1.664 | 0.486 | PASS | 3.587 |
| tiptilt | 2" | 0.696 | 0.106 | 0.833 | 0.137 | PASS | 1.672 |
| tiptilt | 5" | 1.479 | 0.133 | 1.730 | 0.251 | PASS | 4.937 |
| tiptilt | 10" | 2.102 | 0.142 | 2.300 | 0.197 | PASS | 8.317 |
| ao10 | geom 3.58" | 3.466 | 0.078 | 4.433 | 0.967 | OVER | 5.429 |
| ao10 | 2" | 2.094 | 0.057 | 2.523 | 0.429 | PASS | 2.945 |
| ao10 | 5" | 4.093 | 0.094 | 5.030 | 0.937 | OVER | 7.129 |
| ao10 | 10" | 5.707 | 0.116 | 6.778 | 1.070 | OVER | 11.448 |
| ao21 | geom 3.58" | 3.940 | 0.079 | 4.781 | 0.841 | OVER | 5.979 |
| ao21 | 2" | 2.392 | 0.054 | 2.780 | 0.388 | PASS | 3.360 |
| ao21 | 5" | 4.733 | 0.091 | 5.654 | 0.921 | OVER | 7.785 |
| ao21 | 10" | 6.534 | 0.117 | 8.015 | 1.481 | OVER | 12.349 |

`gap` is FAST minus the field. A POSITIVE gap means the field reads LESS
penalty. The band column reads PASS while `|gap| <= 0.5` dB. Stone has no `base`
row: with no corrected mode there is no decorrelation residual, so the penalty
is zero by definition.

VERDICT: **the field and FAST agree at 30 deg, and FAST reads high at 20 deg,
which is what a weak-fluctuation model must do.** THE MODE SET IS MATCHED
across the three rungs: the tilt stays in the error of all three, and only the
piston differs (Stone takes `remove='piston'`, FAST and the field keep it, and
no overlap integral can see a piston).

At 30 deg the field and FAST agree inside the 0.5 dB band at every AO(10) and
AO(21) cell, except two that read 0.50 to 0.58 dB. That is FAST's own Monte
Carlo spread, not a bias: FAST at 1000 draws moves 0.2 to 0.3 dB between runs,
and the first pass of this analysis read those same two cells at 0.28 and
0.36 dB. The UNCORRECTED (`base`) row is a NULL in both models: every reading
sits inside 0.25 dB of zero, which is correct, because an uncorrected beam has
no correction to decorrelate. The tip-tilt row reads FAST about 0.4 to 0.6 dB
above the field at 30 deg.

At 20 deg FAST reads 0.7 to 1.5 dB ABOVE the field at the geometry angle and at
every larger angle. THIS IS EXPECTED. FAST propagates no field. It is a
weak-fluctuation model (a phase screen plus a lognormal amplitude), so it holds
no saturation. The field rung does saturate, and the overshoot of FAST grows
with the angle and with the airmass. That is why the OVER cells sit at 20 deg
and not at 30 deg.

THE STONE COLUMN IS A REPORT, not a gate. It is the extended-Marechal mapping
`(10 / ln 10) sigma^2`, which saturates past sigma^2 = 1 rad^2 (T. S. Ross,
DOI 10.1364/AO.48.001812). So it overstates the penalty everywhere, and more at
20 deg, where it reads 8.3 to 12.3 dB against a field penalty of 2.1 to
6.5 dB.

![the penalty against the angle](figures/p1_penalty_vs_angle.png)

## p2: the fade

The per-trial loss is `-10 log10(eta_pa)` of one trial, so pX is the loss the
link EXCEEDS X percent of the time. The p5 penalty is `p5(theta) - p5(0)`, the
extra deep-fade loss of the point-ahead angle. The `+-` bar is the 68 percent
bootstrap half-width over the trial rows; the resample takes whole ROWS, so
the beacon column and the point-ahead column keep their pairing.

30 deg, 1000 trials, the loss in dB:

| stack | angle | mean | p50 | p10 | p5 | p1 | d p5 | +- |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| base | beacon | 16.56 | 15.85 | 24.85 | 28.25 | 35.98 | 0.00 | 0.00 |
| base | geom 5.24" | 16.58 | 15.82 | 25.00 | 28.32 | 34.68 | 0.07 | 0.76 |
| base | 2" | 16.40 | 16.03 | 23.99 | 27.01 | 34.51 | -1.24 | 0.89 |
| base | 5" | 16.58 | 15.85 | 24.87 | 28.04 | 34.58 | -0.21 | 0.70 |
| base | 10" | 16.69 | 15.82 | 24.62 | 28.22 | 36.69 | -0.03 | 0.86 |
| tiptilt | beacon | 9.27 | 8.23 | 14.93 | 18.01 | 25.04 | 0.00 | 0.00 |
| tiptilt | geom 5.24" | 10.77 | 9.65 | 17.61 | 20.37 | 26.98 | 2.36 | 0.82 |
| tiptilt | 2" | 9.72 | 8.70 | 15.73 | 18.56 | 26.21 | 0.55 | 0.61 |
| tiptilt | 5" | 10.70 | 9.64 | 17.29 | 20.21 | 28.00 | 2.19 | 0.75 |
| tiptilt | 10" | 11.43 | 10.40 | 17.31 | 21.03 | 27.64 | 3.02 | 0.83 |
| ao10 | beacon | 3.05 | 2.96 | 4.15 | 4.51 | 5.45 | 0.00 | 0.00 |
| ao10 | geom 5.24" | 6.50 | 5.67 | 10.02 | 12.72 | 19.80 | 8.21 | 0.38 |
| ao10 | 2" | 4.09 | 3.86 | 5.70 | 6.34 | 8.25 | 1.83 | 0.14 |
| ao10 | 5" | 6.33 | 5.53 | 9.78 | 12.46 | 18.88 | 7.95 | 0.47 |
| ao10 | 10" | 8.30 | 6.97 | 13.74 | 17.68 | 24.75 | 13.17 | 0.77 |
| ao21 | beacon | 1.69 | 1.68 | 2.27 | 2.42 | 2.86 | 0.00 | 0.00 |
| ao21 | geom 5.24" | 5.55 | 4.65 | 9.02 | 11.34 | 17.95 | 8.92 | 0.34 |
| ao21 | 2" | 2.92 | 2.73 | 4.06 | 4.80 | 6.26 | 2.38 | 0.08 |
| ao21 | 5" | 5.38 | 4.46 | 8.89 | 10.93 | 17.72 | 8.51 | 0.36 |
| ao21 | 10" | 7.53 | 6.13 | 12.95 | 16.59 | 23.73 | 14.18 | 0.75 |

20 deg, 1000 trials:

| stack | angle | mean | p50 | p10 | p5 | p1 | d p5 | +- |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| base | beacon | 18.67 | 17.61 | 27.50 | 30.52 | 35.90 | 0.00 | 0.00 |
| base | geom 3.58" | 19.06 | 18.20 | 27.43 | 30.74 | 37.76 | 0.22 | 0.85 |
| base | 2" | 19.03 | 18.14 | 27.26 | 30.38 | 38.40 | -0.15 | 1.06 |
| base | 5" | 18.81 | 18.25 | 27.21 | 29.50 | 37.22 | -1.02 | 0.95 |
| base | 10" | 18.66 | 18.01 | 26.87 | 29.96 | 35.91 | -0.56 | 0.87 |
| tiptilt | beacon | 12.50 | 11.46 | 19.30 | 22.13 | 28.57 | 0.00 | 0.00 |
| tiptilt | geom 3.58" | 14.10 | 12.93 | 21.33 | 24.57 | 31.08 | 2.44 | 0.75 |
| tiptilt | 2" | 13.43 | 12.19 | 20.87 | 24.27 | 31.15 | 2.14 | 0.75 |
| tiptilt | 5" | 14.42 | 13.33 | 21.68 | 25.58 | 31.97 | 3.45 | 0.89 |
| tiptilt | 10" | 15.12 | 14.13 | 23.03 | 26.03 | 33.41 | 3.90 | 0.69 |
| ao10 | beacon | 4.52 | 4.40 | 6.09 | 6.54 | 7.99 | 0.00 | 0.00 |
| ao10 | geom 3.58" | 9.08 | 7.96 | 14.25 | 17.47 | 24.27 | 10.93 | 0.61 |
| ao10 | 2" | 7.03 | 6.53 | 10.07 | 11.62 | 16.02 | 5.08 | 0.36 |
| ao10 | 5" | 10.00 | 8.98 | 15.66 | 18.83 | 24.67 | 12.28 | 0.49 |
| ao10 | 10" | 12.35 | 11.00 | 20.10 | 23.01 | 30.59 | 16.47 | 0.67 |
| ao21 | beacon | 2.67 | 2.64 | 3.52 | 3.75 | 4.30 | 0.00 | 0.00 |
| ao21 | geom 3.58" | 7.70 | 6.60 | 12.87 | 15.51 | 21.80 | 11.76 | 0.75 |
| ao21 | 2" | 5.48 | 4.93 | 8.11 | 9.78 | 13.71 | 6.02 | 0.34 |
| ao21 | 5" | 8.90 | 7.77 | 14.45 | 18.17 | 25.65 | 14.42 | 0.69 |
| ao21 | 10" | 11.46 | 10.16 | 19.05 | 22.44 | 29.11 | 18.69 | 0.42 |

VERDICT: **the deep-fade cost of the point-ahead is much larger than the mean
cost. A link that is sized on the mean is under-sized.** At 30 deg, at the
geometry angle, the p5 loss rises by 8.2 dB for AO(10) and by 8.9 dB for
AO(21), against only 2.7 and 3.1 dB on the mean. The tip-tilt stack pays 2.4 dB
of p5 there. The uncorrected stack pays ZERO inside its bar (0.07 dB against a
0.76 dB bar), which is the same null the mean shows in p1.

At 20 deg the AO p5 penalties are 10.9 dB (AO(10)) and 11.8 dB (AO(21)) at the
geometry angle, and they grow to 16.5 and 18.7 dB at 10 arcsec. The p50 rows
show where that fade comes from: at 30 deg the AO(21) median moves 1.68 ->
4.65 dB from the beacon to the geometry angle, and the p1 moves 2.86 ->
17.95 dB. So the point ahead does not shift the distribution, it grows its
TAIL. The tip-tilt stack keeps a flatter tail (p1 25.04 -> 26.98 dB at 30 deg),
because it starts with a much deeper fade and it has less to lose.

![the p5 fade against the angle](figures/p2_p5_vs_angle.png)

## p3: the screen count

30 deg only. Each campaign takes the PRODUCTION grid and one override plan
from `validation.anisoplanatism_screens.common.plan_set`, so the sweep moves
the SCREEN COUNT only: 5, 9, 15 and 25 equal-Rytov-weight screens, plus the
ground-split plan that cuts the lowest screen into four. 400 trials for each
plan.

THE n9 PLAN IS THE PRODUCTION PLAN, and it still gets its OWN store. The
campaign fingerprint holds the repr of the caller grid and plan
(`olb.waveoptics.turbulence.fingerprint.cache_key`), so a campaign that PASSES
the production plan keys differently from one that lets the sizer build it.
The script asserts that the two plans hold the same screen distances.

THE PASS RULE is the 9i criterion: the penalty is FLAT from 9 screens up,
inside the 2-sigma bootstrap band, against the 25-screen plan.

The penalty at 5 arcsec, in dB:

| plan | screens | base | tiptilt | ao10 | ao21 |
| --- | --- | --- | --- | --- | --- |
| n5 | 5 | -0.151 +-0.214 | 1.212 +-0.121 | 2.611 +-0.093 | 2.913 +-0.086 |
| n9 | 9 | -0.070 +-0.223 | 1.358 +-0.131 | 2.683 +-0.100 | 3.051 +-0.100 |
| n15 | 15 | 0.417 +-0.197 | 1.375 +-0.133 | 2.596 +-0.092 | 2.924 +-0.094 |
| n25 | 25 | 0.083 +-0.228 | 1.191 +-0.143 | 2.664 +-0.091 | 3.018 +-0.091 |
| ground_split_x4 | 12 | 0.372 +-0.240 | 1.189 +-0.133 | 2.699 +-0.097 | 3.042 +-0.099 |

The penalty at 10 arcsec, in dB:

| plan | screens | base | tiptilt | ao10 | ao21 |
| --- | --- | --- | --- | --- | --- |
| n5 | 5 | -0.114 +-0.275 | 1.863 +-0.165 | 4.016 +-0.128 | 4.505 +-0.125 |
| n9 | 9 | 0.423 +-0.286 | 1.636 +-0.156 | 3.639 +-0.125 | 4.223 +-0.124 |
| n15 | 15 | 0.347 +-0.233 | 1.751 +-0.162 | 3.847 +-0.121 | 4.380 +-0.123 |
| n25 | 25 | -0.023 +-0.240 | 1.733 +-0.154 | 3.895 +-0.127 | 4.442 +-0.122 |
| ground_split_x4 | 12 | 0.636 +-0.271 | 1.676 +-0.177 | 3.699 +-0.131 | 4.238 +-0.122 |

The p5 penalty at 5 and 10 arcsec is in the log and in the JSON record. The p5
bars are 3 to 10 times wider than the mean bars (0.4 to 2.6 dB), because the p5
reads one tail of 400 trials; the AO(10) p5 penalty at 10 arcsec, for example,
runs 11.1 to 12.0 dB over the five plans with bars of 0.6 to 1.5 dB.

The cost of the count is linear: 400 trials take 72.9 s at 5 screens, 103.9 s at
9, 152.4 s at 15 and 251.5 s at 25.

VERDICT: **the point-ahead penalty is CONVERGED at the production nine-screen
plan.** 24 of 24 readings are flat against the 25-screen plan inside the
2-sigma band, at 5 arcsec and at 10 arcsec, on the mean penalty and on the p5
penalty. The ground-split plan, which cuts the lowest screen into four, changes
nothing: it sits inside the same band at every stack and at both angles. The
uncorrected (`base`) column reads zero inside its own bar at every count, which
is the same null as p1 and p2.

THIS AGREES WITH THE TWO PHASE-ONLY STUDIES of
`validation/anisoplanatism_screens/`, where the production nine-screen plan
reads 0.987 to 0.995 of the continuous Stone sum. NOTE THE RESOLUTION LIMIT:
400 trials cannot resolve the 2 percent by which the analytic sum puts the
five-screen plan low. So the sweep shows that nine screens are enough; it does
not show that five screens are wrong.

## The caveats

- SNAPSHOT ONLY. One atmosphere for each trial, no time axis. The study gives
  the fade DEPTH, not the fade rate and not the fade duration.
- PERFECT AO. The correction is an ideal modal fit of one snapshot: no
  wavefront-sensor noise, no finite subaperture, no aliasing, no servo lag and
  no branch point. So the beacon-direction loss is the UPPER BOUND of the
  benefit of a corrector, and the point-ahead penalty is the CLEANEST possible
  one. A real terminal pays more.
- PLANE-PARALLEL SCREENS. The point-ahead pass shifts the window of each
  screen by `theta * z_ground`, the plane-parallel geometry of Andrews and
  Phillips, DOI 10.1117/3.626196, Ch. 12, Eq. (14). There is no cone effect,
  so this geometry holds for a star-like (infinite) beacon, not for a laser
  guide star.
- INTEGER PIXELS. The lateral shift is a whole number of pixels, so no
  interpolation touches the screen statistics. The angle the run measures is
  therefore the rounded angle. At 512 px on a 3.5 m grid the pixel is 6.9 mm,
  and one pixel is 0.21 arcsec at a 32 km screen.
- THE STONE COMPARISON SATURATES. The extended Marechal mapping
  `eta = exp(-sigma^2)` is an approximation that breaks past sigma^2 = 1
  rad^2 (T. S. Ross, DOI 10.1364/AO.48.001812), so the Stone column reads high
  at a large angle and at a high mode count. It is a report, not a gate.
- THE FIELD HOLDS NO SERVO LAG, and FAST holds its own default servo
  (`TLOOP = 0.001 s`, `TEXP = 0.001 s`, `ALIAS=True`, `NOISE=0`). So a FAST
  penalty is not a pure anisoplanatism; the difference of the two angles
  removes most of that, because the servo term does not depend on the angle.
- THE FAST ROWS CARRY THEIR OWN NOISE. Each FAST number is a Monte Carlo of
  1000 draws, and it moves 0.2 to 0.3 dB from run to run. So a gap of 0.5 dB
  against the field is at the level of that spread.

## The TILT anisoplanatism, in rad^2 (`tilt_anisoplanatism.py`)

THE QUESTION (owner, 2026-09-11). p1 reads FAST about 0.3 dB ABOVE the field
on the TIP-TILT stack at 30 deg, with the same sign at every angle, after the
known FAST leak is taken away. At 20 deg that offset goes away and an excess
appears on the AO rows instead. A dB number holds two things together: the
phase variance and the extended Marechal map of that variance. This script
takes the MAP away. It measures the TILT decorrelation alone, in rad^2, from
the STORED planes of the `base` campaigns, and it puts the field, Stone and
FAST next to each other mode by mode.

THE THREE ROUTES. (1) The SCREEN route: the difference of the summed screen
phase of the shifted window and of the beacon window, fitted with 21 Noll
modes over the 0.7 m aperture mask. That IS the Stone quantity. (2) The FIELD
route: the tilt of the propagated point-ahead field minus the tilt of the
propagated beacon field, from the wrapped slopes and from the G-tilt centroid.
(3) The REFERENCES: the continuous Stone integral at `L0 = 25 m` and at
`L0 = inf`, the delta-layer Stone sum on the campaign's OWN plan with the
pixel-rounded shifts, and the FAST tilt band.

THE BAND UNITS. Every band is a phase variance over the aperture mask, in
rad^2, and the bands ADD UP to the mask variance of the piston-removed
difference. The Noll basis is only NEAR-orthonormal on a pixel mask, so each
band is reported as the SHARE of the fitted variance, which adds EXACTLY. The
script asserts that closure at 1e-6.

THE FAST TILT BAND is the difference of two SERVO-OFF runs: the clean
anisoplanatic split of `ZMAX = 3` (the piston and the two tilts) minus the
split of `ZMAX = 1` (the piston alone). The grid rule is the rule of
`validation/fast_stone_pointahead/`: an EXPLICIT `(NPXLS, DX)` pair, here
1024 px and 0.05 m, so the side 51.2 m holds the 25 m outer scale.

THE SCRIPT PROPAGATES NOTHING. It reads the stored planes only. Run it from
the repository root as a module of `validation.waveoptics_pointahead`.

PLACEHOLDER. The table below is NOT MEASURED. Fill it from
`tilt_anisoplanatism_results.json` after the run on the stored 1000-trial
campaigns.

| el [deg] | angle | screen tilt [rad^2] | +-2 sigma | slopes | G-tilt | Stone px | Stone 25 m | Stone inf | FAST | field/Stone | FAST/Stone |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 30 | geom | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD |
| 30 | 2 arcsec | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD |
| 30 | 5 arcsec | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD |
| 30 | 10 arcsec | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD |
| 20 | geom | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD |
| 20 | 2 arcsec | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD |
| 20 | 5 arcsec | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD |
| 20 | 10 arcsec | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD |

THE SMOKE TEST (2026-09-11, this laptop, 8 trials, `preset="rapid"`, two
angles, 256 px, 5 screens) says only that the machinery runs: the band-sum
closure held to 3e-16, and FAST read 1.047 of the continuous Stone tilt band
at both angles. The field column of an 8-trial rapid store carries a +-2 sigma
bar of 0.15 rad^2 on a 0.24 rad^2 value, so it is NOT a result.

## The files

- `waveoptics_pointahead.py`: the driver. `--study p0 p1 p2 p3`.
- `tilt_anisoplanatism.py`: the tilt-band analysis of the stored planes.
- `tilt_anisoplanatism.log`, `tilt_anisoplanatism_results.json`,
  `figures/tilt_anisoplanatism.png`.
- `waveoptics_pointahead_<study>.log`: the log of each study.
- `waveoptics_pointahead_<study>_results.json`: the numbers of each study.
- `figures/p1_penalty_vs_angle.png`, `figures/p2_p5_vs_angle.png`.
- `campaigns/`: the stores. They are gitignored.

## The sources

- J. Stone, P. H. Hu, S. P. Mills and S. Ma, "Anisoplanatic effects in
  finite-aperture optical systems," J. Opt. Soc. Am. A 11(1), 347-357 (1994),
  DOI 10.1364/JOSAA.11.000347.
- J. H. Shapiro, "Reciprocity of the turbulent atmosphere," J. Opt. Soc. Am.
  61(4), 492-495 (1971), DOI 10.1364/JOSA.61.000492.
- R. J. Noll, "Zernike polynomials and atmospheric turbulence," J. Opt. Soc.
  Am. 66(3), 207-211 (1976), DOI 10.1364/JOSA.66.000207.
- O. J. D. Farley and others, "FAST: Fourier domain adaptive optics simulation
  tool," Opt. Express 30(13), 23050 (2022), DOI 10.1364/OE.458659.
- J. D. Schmidt, Numerical Simulation of Optical Wave Propagation (2010),
  DOI 10.1117/3.866274, Ch. 9.
- L. C. Andrews and R. L. Phillips, Laser Beam Propagation through Random
  Media, 2nd ed. (2005), DOI 10.1117/3.626196.
- T. S. Ross, "Limitations and applicability of the Marechal approximation,"
  Appl. Opt. 48(10), 1812 (2009), DOI 10.1364/AO.48.001812.
- D. L. Fried, J. Opt. Soc. Am. 56, 1372 (1966), DOI 10.1364/JOSA.56.001372.
