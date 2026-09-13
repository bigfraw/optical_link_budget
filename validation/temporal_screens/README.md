# temporal_screens — the strip screens of the frozen-flow time axis

The fidelity-2 layer gives SNAPSHOTS. A fade RATE and a fade DURATION need a
screen stack that MOVES. `olb/waveoptics/turbulence/temporal.py` makes that
stack: each layer gets ONE oversized RECTANGULAR screen, a STRIP, and a frame
is a crop of the strip at an integer pixel offset (Taylor frozen flow,
DOI 10.1098/rspa.1938.0032).

These five scripts are gates (a) to (e) of that build. Run each from the
repository root.

| script | gate | what it measures |
|---|---|---|
| `rect_factory.py` | (a) | `ScreenFactory(nx=...)`: `nx = n` is bit-identical; the strip holds the structure function and the Z-tilt of a square screen |
| `strip_taylor.py` | (b) | one pixel of the moving frame: Taylor, the -8/3 temporal spectrum, and the seam |
| `frame0_parity.py` | (c) | frame 0 of a record against a drawn snapshot: the aperture scintillation index and the mean fibre coupling |
| `tilt_spectrum.py` | (d) | the Z-tilt spectrum of one record: the -2/3 law and the corner, against the Greenwood frequency |
| `hero_temporal.py` | (e) | the fade RATE and the fade DURATION at the 5 percent level, from 8 records of 2 s at 30 and 20 deg. SMOKE DONE 2026-09-13 (3.6 s); the FULL run DONE 2026-09-13 on the bigfraw GPU (727 s) |
| `rotated_strip_gate.py` | — | the ROTATED thin strip of a crosswind against the axis-aligned box and against the unrotated cut of the same crop (backlog 2-P1b item 10) |
| `route_equivalence.py` | — | are the THREE strip routes (along track, box, rotated) the SAME screen? D(r), the tilt, the piston-free and the raw variance, and the power spectrum in radial bands, over 64 seeds |
| `rotation_taper.py` | — | six cures for the seam ringing of the three shears, against the rotation of the full periodic screen. It picks the per-shear taper (V4) and the roll-off rule of `TemporalSpec.rot_margin` |
| `rotated_record_parity.py` | — | one 4000-frame ROTATED record against 4000 independent snapshots: the scintillation index, the mean power and the 5 and 1 percent fade of three receivers |
| `record_plots.py` | — | the PICTURES of record 0 at 30 deg: the power against time over 2 s, and a 0.1 s animation of the phase, the aperture intensity, the fibre tip and the fibre power around the deepest fade. It only READS the record. DONE 2026-09-13 (55 s): the deepest fibre fade is 40.3 dB under the median at t = 96.5 ms |
| `record_ao_plots.py` | — | the SAME record under four PERFECT-AO stacks (none, TipTilt, AO(10), AO(50)): the four fibre power series over 2 s, and the same 0.1 s animation with one ROW for each stack. It only READS the record and it corrects each stored field post hoc. DONE 2026-09-13 (113 s, the SLOPES fallback) |

```
python -m validation.temporal_screens.rect_factory
python -m validation.temporal_screens.strip_taylor
python -m validation.temporal_screens.frame0_parity
python -m validation.temporal_screens.tilt_spectrum
python -m validation.temporal_screens.hero_temporal --smoke
python -m validation.temporal_screens.record_plots
python -m validation.temporal_screens.record_ao_plots
python -m validation.temporal_screens.record_ao_plots --source screens
```

## Gate (e), the run lines

`hero_temporal.py` keeps ONE campaign for each (elevation, record) under
`campaigns/el<deg>/r<record>` (not tracked), and the frames of a record ARE the
trials of its campaign. It DELETES the strips of a complete record, because the
seed rebuilds them; `--keep-strips` holds them. `--analyse` reads what is stored
and computes no frame. The environment variable `OLB_TEMPORAL_ROOT` moves the
campaigns to another disk.

The SMOKE run holds 30 deg, one record, 200 frames, the `rapid` preset and the
host backend. It ran end to end on the laptop on 2026-09-13 in **3.6 s** (0.015
s per frame, a 256 px grid, 5 screens): 3 fibre events and 2 bucket events, both
INFO because 0.1 s of record cannot reach the 20 events that a usable bar needs.
It wrote `campaigns/smoke/el30/r0/` (one block file, the manifest, the patch
indices), `data/hero_temporal_30.csv` and `figures/hero_temporal_30.png`.

The FULL run takes the CUDA device of `bigfraw`. The GPU python is
`C:\Users\alexf\olb-gpu-venv\Scripts\python.exe`. The login shell of
`ssh desktop` is Windows PowerShell 5.1, so chain with `;`, never with `&&`, and
launch the long run through WMI so it outlives the ssh session
(`validation/campaign_resources/README.md`):

```
$cmd = 'cmd /c "cd /d D:\repos\optical_link_budget && C:\Users\alexf\olb-gpu-venv\Scripts\python.exe -u -m validation.temporal_screens.hero_temporal > validation\temporal_screens\hero_temporal.log 2>&1"'
Invoke-CimMethod -ClassName Win32_Process -MethodName Create -Arguments @{CommandLine=$cmd}
```

`--store-screen-phase` keeps the summed screen phase of every frame. The
perfect-AO read of a SPACE link senses that phase (`record_ao_plots.py`). It is
OFF by default, and the flag enters the campaign fingerprint, so a record that
holds the phase is a NEW campaign and the default run stays byte-identical.

The defaults ARE the production settings: `--elevations 30 20`, `--records 8`,
`--frames 4000`, `--dt 5e-4`, `--preset standard`, `--block-size 500`,
`--fft-backend cupy`, `--workers` unset (one device, one stream). The block size
is clamped to the frame count, and it must DIVIDE it, so a campaign never runs a
frame past the travel that its strips hold.

Each script prints PASS or FAIL bands, writes its CSVs to `data/` (not
tracked) and its figures to `figures/` (tracked). `strip_taylor.py` writes its
strip to `strips/` and removes it at the end. `frame0_parity.py` and
`tilt_spectrum.py` do the same under `strips_frame0/` and `strips_tilt/`.

## Results, 2026-09-13

Gate (a), a (512, 4096) strip at dx = 1 cm, r0 = 10 cm, L0 = 25 m, against a
512 square ensemble and against the analytic von Karman law. 14 of 14 bands
hold.

| r [m] | D(r) along y / law | D(r) along x / law | square / law |
|---|---|---|---|
| 0.03 | 0.989 | 0.980 | 0.993 |
| 0.16 | 0.997 | 0.980 | 1.009 |
| 0.32 | 0.994 | 0.970 | 1.016 |
| 1.00 | 0.991 | 0.947 | 1.021 |

The Z-tilt of a 1 m pupil over 32 crops reads 0.89 of the square value, inside
the 2 SE band (the tilt is a noisy statistic, so that band is wide).

Gate (b), one layer at 50.6 m/s, dt = 0.5 ms, 4096 frames, a (128, 15520)
strip. 9 of 9 bands hold.

- The temporal D(tau) equals the spatial D(v tau) of the same strip to 2.4
  percent, at every lag from 2 to 64 ms.
- The temporal spectrum fits an exponent of **-2.643 +/- 0.054** (95 percent)
  over 10 to 200 Hz, against the -8/3 = -2.667 law.
- The last window ends 22 pixels before the seam pad, and the mean correlation
  of 40 disjoint window pairs at the travel separation is +0.065 +/- 0.078.

Gate (c), the hero downlink at 30 deg, `rapid` preset, single precision, seed
20260913, 256 trials of each arm. 3 of 4 bands hold.

| statistic | snapshot | frame 0 | ratio | 2 SE of the difference |
|---|---|---|---|---|
| aperture sigma2_I | 0.0042 | 0.0039 | 0.939 | 0.251 |
| mean SMF eta | 0.03940 | 0.04095 | 1.039 | 0.215 |

Both ratios sit well inside the 2 SE bar, so **a strip crop is the same
atmosphere as a drawn screen stack**. The fixed 5 percent band is narrower than
the error bar at 256 trials, so the index misses it at 0.939; that band is a
target, not a test, at this trial count. A 128-trial run of the same seed stream
read 0.64 on the index, which was the tail of the index estimator and not a
deficit.

Gate (d), one record of 4000 frames at 0.5 ms (2 s), the same case. 1 of 3 bands
holds, and the two failures are physics, not the strips.

| quantity | value |
|---|---|
| low-frequency exponent, 1 to 5 Hz | **-0.559** (the law is -2/3 = -0.667) |
| high-frequency exponent, 60 to 400 Hz | -5.26 |
| measured corner | 38.2 Hz |
| Greenwood f_G, strip velocity | 641.5 Hz |
| Greenwood f_G, default Bufton wind | 1095.3 Hz (1.7x: it carries its own slew) |
| Tyler tilt frequency f_T | 15.2 Hz |
| pupil corner 0.3 v/D, slowest layer (11.8 m/s) | 5.1 Hz |
| pupil corner 0.3 v/D, fastest layer (131.3 m/s) | 56.3 Hz |

The **-2/3 law holds**, so the frozen-flow axis carries the right tilt spectrum.
The corner does NOT sit at f_G, and it must not: f_G is a PHASE quantity that
does not know the pupil, and the tilt of a 0.7 m pupil breaks on the V/D scale
(Tyler, DOI 10.1364/JOSAA.11.000358). The measured 38 Hz sits between the
corner of the slowest layer and the corner of the fastest one, which is where
the blend of a five-layer stack must sit.

Gate (e), the hero downlink, `standard` preset, single precision, seed
20260913, the cupy backend on bigfraw. It ran 8 records of 4000 frames at
dt = 0.5 ms for each elevation, which is 16 s of record for each elevation. It
took **727 s** in total, 0.010 to 0.011 s per frame at 512 px. The full log is
`data/hero_temporal.log` (not tracked), the series are
`data/hero_temporal_{30,20}.csv` and the figures are
`figures/hero_temporal_{30,20}.png`.

The 5 percent fade level is the 5th percentile of the pooled series. An event
is a maximal run of frames under that level. The rate bar is Poisson and the
duration bar is 2 SE.

| case | level below the median | events | rate [1/s] | mean duration [ms] | median [ms] | longest [ms] |
|---|---|---|---|---|---|---|
| 30 deg SMF | 12.60 dB | 596 | 37.25 +/- 1.53 | 1.342 +/- 0.089 | 1.000 | 9.5 |
| 30 deg bucket | 0.49 dB | 407 | 25.44 +/- 1.26 | 1.966 +/- 0.143 | 1.500 +/- 0.071 | 7.0 |
| 20 deg SMF | 12.33 dB | 562 | 35.12 +/- 1.48 | 1.423 +/- 0.100 | 1.000 | 11.5 |
| 20 deg bucket | 0.85 dB | 291 | 18.19 +/- 1.07 | 2.749 +/- 0.215 | 2.500 +/- 0.464 | 9.5 |

Both elevations hold the 20-event band, so both PASS. The record grid is
512 px, 6.86 mm at 30 deg and 8.68 mm at 20 deg. The Greenwood frequency at the
strip velocities is 658.7 Hz at 30 deg and 684.4 Hz at 20 deg, and tau0 is
1.518 ms and 1.461 ms; those are context, not a test (see gate (d)).

TWO CAUTIONS.

1. The median SMF event lasts 1.0 ms, which is 2 frames. So dt = 0.5 ms only
   just resolves the SMF fade duration. A finer dt is the next study.
2. These are the FIRST temporal fade numbers of the package, and no reference
   model checks them. `olb/turbulence/andrews/temporal.py` holds an analytic
   fade rate and fade duration that have no external check (backlog 0-N6), so a
   temporal fade reference model is the other next study.

## The perfect-AO pictures, 2026-09-13

`record_ao_plots.py` reads record 0 at 30 deg and it corrects every stored
field four ways: no correction, TipTilt (3 Noll modes), AO(10) and AO(50). The
correction is the post-hoc perfect-AO read of the package
(`Campaign.recouple_compensated`), so it is an ideal modal fit of one snapshot
with no sensor noise, no servo lag and no anisoplanatism. Every corrected
number is therefore an UPPER BOUND. Source: Noll,
DOI 10.1364/JOSA.66.000207 (the mode order and the count).

THE SENSING SOURCE. A SPACE link senses the SUMMED SCREEN PHASE. Record 0 was
run BEFORE the `--store-screen-phase` flag existed, so it holds none, and the
run below took the SLOPES fallback. The slope route reads a high-order
correction 4 to 6 percent LOW on a space link, so the table below is
INDICATIVE. Run one record with `--store-screen-phase` and read it with
`--source screens` for the route of record.

The run took **113 s** on the laptop (4000 frames, 3 corrected passes of 3.0 to
4.4 s each, then a 200-frame animation). The columns are dB over the
UNCORRECTED median. A fade event is a maximal run of frames under a level. The
"shared" columns use the UNCORRECTED 5 percent level, so every stack counts the
SAME fades; the "own" columns use the 5 percent level of that stack.

| stack | Noll modes | mean eta | median [dB] | 5 percent [dB] | shared N | shared mean [ms] | own N | own mean [ms] |
|---|---|---|---|---|---|---|---|---|
| none | 0 | 0.0413 | 0.00 | -11.72 | 77 | 1.299 | 77 | 1.299 |
| TipTilt | 3 | 0.1445 | +6.98 | -2.10 | 6 | 1.250 | 33 | 3.030 |
| AO(10) | 10 | 0.4029 | +11.89 | +9.78 | 0 | — | 28 | 3.571 |
| AO(50) | 50 | 0.6502 | +13.91 | +13.31 | 0 | — | 53 | 1.887 |

THE DEEP FADES GO FIRST. The tip-tilt corrector removes 71 of the 77 deep
fades of the record, and both AO stacks remove all 77. The 5 percent level
rises by 9.6 dB (TipTilt), 21.5 dB (AO(10)) and 25.0 dB (AO(50)), which is much
more than the median gain, so the correction cuts the TAIL harder than the
bulk. That matches the campaign result of `validation/waveoptics_ao/` (9.4 /
22.1 / 24.3 dB of p5 at 30 deg), and the two agree although this record senses
the slopes and that campaign senses the screens.

The figures are `figures/record_ao_timeseries_30.png` (0.25 MB) and
`figures/record_ao_field_30.gif` (14.7 MB, 200 frames at 20 fps, dpi 60). The
animation rows share the colour scale of each column, so the phase flattens and
the fibre-tip spot sharpens down the rows.

## What the two structure gates found

**The strip needs its own low-frequency rule.** The row `fy = 0` of a
rectangular grid carries the power of every frequency with |fy| < dfy/2, and
the main grid gives all of it to `fy = 0`, which holds NO structure along y. On
a square grid that row is thin and the 3 by 3 subharmonics repair it. On a
strip the row is the widest part of the spectrum, because dfy is many times
dfx. Two wrong routes were measured before the right one:

1. **Per-axis subharmonics** (a 3 by 3 set at dfy/3^p by dfx/3^p). D(r) along y
   read 8 to 23 percent LOW and D(r) along x read 1 to 12 percent HIGH.
2. **The square subharmonic box** (keep the spacing of the short axis and cut
   the matching box out of the main grid). D(r) came right on both axes, but
   the low frequencies along x became a COMB of a few tones with a period of
   `3^P * side_y` = 34.6 m, so a 103 m record walked back into its own start.

The route of record is the **band**: the `fy = 0` row is cut out and drawn
again as 2P+1 sub-rows in y, each keeping the FINE x sampling of the strip.
Both axes then hold the analytic law inside 1 to 5 percent, and the long axis
carries a continuum.

**One window pair is not a seam test.** A 1.28 m window of a screen with a 25 m
outer scale is almost a plane. After the mean of each window goes, its
correlation with another window is little more than the sign agreement of two
tilts, so it lands near +1 or -1 at random, even for two windows 2 m apart.
Only the MEAN over many pairs says anything.

## The crosswind box: the memory of one strip build, 2026-09-13

`table_precision.py` measures the two memory levers of `ScreenFactory` on a
strip build (backlog 2-P1b item 10). A crosswind grows the SHORT axis of the
axis-aligned box, and the 30 deg hero at `wind_dir_deg = 90` asks for a
11897 x 43744 box (520 Mpx) for one jet-level layer. Before the change the
build drew the noise through two float64 grids and two complex128 grids (48
bytes per pixel, 25 GB), and it built the filter through float64 `meshgrid`
tables, so the box swapped and took 170 s.

Two changes are BIT-IDENTICAL on every route of record (the old and the new
module were run side by side on eight seeded cases, square and strip, float64
and float32): the noise is written straight into the real and the imaginary
part of one complex grid, and the transform chain rebinds one name, so one grid
dies as the next is born. One change is an OPT-IN: `table_dtype=np.float32`
builds the filter tables in single precision.

| route, 520 Mpx box on bigfraw | wall | peak working set |
|---|---|---|
| before (float64 meshgrid, four-grid draw) | 170 s, swapping | above the 32 GB box |
| float64 tables, new draw and chain | 151.5 s | 19.99 GiB |
| `table_dtype=np.float32` | 137.4 s | 15.64 GiB |

The float32 tables against the float64 tables on the same seed (512 x 4096,
rms 22 rad): the maximum phase difference is 1.14e-5 rad, the relative rms
difference 8.0e-8, and D(r) along y and along x agree to six digits at every
lag from 3 cm to 1 m. GATE PASS. So the box FITS the host now, and the
generation time is the FFT of the 23x area, not the memory: a 9-layer crosswind
record still costs about 20 min of strip build against 22.7 s along track. The
route that cuts the TIME is the exact Fourier three-shear rotation of a thin
strip (queued; a real-space bilinear or cubic rotation is NOT the tool, it
smooths the Fresnel scale). The speed-only mapping (every layer along the slew
axis at its resultant speed) is REJECTED: it holds for a rotation-invariant
scalar receiver only, and the downlink phase field in time depends on each
layer's velocity DIRECTION.

## The rotated thin strip, 2026-09-13

`rotated_strip_gate.py` measures the ROTATED thin strip (backlog 2-P1b item
10, route (b)). `TemporalSpec(rotated=True, rot_margin=0.5)` holds ONE thin
strip for each layer along that layer's RESULTANT velocity, and each frame
takes a PADDED crop and turns it back with `rotate_fourier`, the exact Fourier
three-shear rotation (Unser, Thevenaz and Yaroslavsky, DOI 10.1109/83.469963).
The opt-in is OFF by default, so the along-track route does not move.

THE REFERENCE IS NOT THE BOX ALONE. The three strips do NOT share a shape. The
box short axis is 1540 to 1719 px in the gate case, the rotated thin strip is
412 to 876 px, and the ALONG-TRACK strip of the route of record is `n` = 256
px. The short axis sets how much of the outer scale a strip holds on that axis,
so the three give different frame variances and different large-r structure,
and that difference is NOT the rotation. The along-track route, which is the
route of record, differs from the box by 0.57 to 1.69 on the same field
statistics. So the controlled test is the rotated frame against the UNROTATED
cut of the SAME crop: one strip, one seed, one grid.

The case is `n = 256`, `dx = 1 cm`, `r0 = 10 cm`, `L0 = 25 m`, 3 seeds, 201
frames over 3 m of travel, `pad_outer_scales = 0.5`. 19 of 20 bands hold.

Gate 1, the rotation of ONE array, D(r) after / before. The band is 3 percent.

| case | margin | D(1 px) | D(4 px) | outer-ring PSD |
|---|---|---|---|---|
| 6.5 deg | 0.25 n | 1.014 | 0.996 | 1.083 |
| 6.5 deg | 0.50 n | 1.008 | 1.006 | 0.852 |
| 6.5 deg | 1.00 n | 0.999 | 0.998 | 1.072 |
| 45 deg | 0.25 n | **1.040 FAIL** | 1.026 | 1.175 |
| 45 deg | 0.50 n | 1.011 | 1.017 | 1.055 |
| 45 deg | 1.00 n | 0.999 | 0.993 | 1.109 |

THE MARGIN ANSWER: 0.25 n is NOT enough at 45 deg, 0.5 n holds inside 2
percent, and 1.0 n is clean inside 0.8 percent. The CONTROL, a real-space
bilinear (order 1) rotation of the same crop, reads D(1 px) at 0.815 (6.5 deg)
and 0.851 (45 deg): it loses 15 to 18 percent of the smallest-scale structure.
That is why bilinear is not the tool.

Gate 2, the frame-to-frame consistency: D(tau) of the rotated frames against
D(tau) of the unrotated cut of the SAME crops, at `rot_margin = 0.5`. The band
is 5 percent, because D(tau) at the shortest lag is built almost entirely from
the highest spatial frequencies.

| case | worst \|ratio - 1\| | at the shortest lag |
|---|---|---|
| 6.5 deg | 0.010 | 0.992 |
| 45 deg | 0.032 | 0.976 |

Gate 3, one screen and a 5 km hop, `r0 = 30 cm`, a 0.5 m pupil, 96 frames.
The band is the 5 percent P2 kill line, on rotated / unrotated.

| quantity | 6.5 deg rot/unrot | 45 deg rot/unrot | 45 deg rot/box | 45 deg along track / box |
|---|---|---|---|---|
| point sigma2_I | 1.010 | 1.001 | 0.616 | 0.690 |
| aperture sigma2_I | 1.000 | 0.992 | 0.655 | 0.676 |
| mean SMF eta | 1.000 | 1.000 | 1.066 | 1.483 |

The rotation changes NOTHING that the field sees. The rot/box and the
along-track/box columns show the strip-SHAPE effect, and the ALONG-TRACK route
of record moves as far from the box as the rotated one does.

Gate 4, the cost. The 30 deg hero crosswind layer, PRODUCTION pad (2 L0), on
bigfraw:

| route | shape | pixels | build |
|---|---|---|---|
| box | (11897, 43744) | 520.4 Mpx | 137.4 s |
| rotated, 0.25 n | (824, 44320) | 36.5 Mpx | (14.3x smaller) |
| rotated, 0.50 n | (1080, 44576) | 48.1 Mpx | **6.7 s (20.4x)** |
| rotated, 1.00 n | (1592, 45088) | 71.8 Mpx | (7.3x smaller) |

One rotated frame of that 1080 px crop costs **6.7 ms on the cupy backend**
and 285 ms on numpy. So a 9-layer crosswind record trades about 20 min of box
build for about 63 s of strip build plus about 60 ms of rotation per frame.
At 4000 frames that is 240 s of rotation, which is 6x the propagation of one
512 px frame. A margin rule by ANGLE (0.25 n holds at the small angles of the
fast layers) and a per-layer crop side would cut it.

NOT DONE: the runner and `Campaign` do not read `rotated` yet. The gate drives
`strip_plan`, `build_strips` and `frame_stack` directly.

```
python -m validation.temporal_screens.rotated_strip_gate
python -m validation.temporal_screens.rotated_strip_gate --quick --bilinear
python -m validation.temporal_screens.rotated_strip_gate --build-hero --fft-backend cupy
```

## Are the three routes the same screen? 2026-09-13

`route_equivalence.py` answers the question the gate above left open. THREE
routes give one frame of one layer, and each one holds a rectangular Fourier
screen of a different SHORT axis: A, the ALONG-TRACK strip of record (short
axis exactly `n`); B, the axis-aligned BOX of the 2-D walk (short axis `n` plus
the perpendicular travel plus the seam pad); C, the ROTATED thin strip (short
axis `n(|cos| + |sin|)` plus the crop margin).

The script reads 8 windows across the WHOLE long axis of each strip, for 64
SEEDS, and it averages the windows of a seed first, so the standard error runs
over independent atmospheres. The case is `n = 256`, `dx = 1 cm`, `r0 = 10 cm`,
`L0 = 25 m`, a 0.70 m receive aperture, the PRODUCTION seam pad of 2 L0 and the
0.5 n rotation margin. The strips are (256, 5280) = 2.56 m short axis,
(5257, 5280) = 52.6 m, and (540, 5568) = 5.40 m at 6.5 deg or (620, 5632) =
6.20 m at 45 deg. 87 of 108 bands hold, and the ratios that miss are listed
below with their bars. The reference of every ratio is the BOX.

**THE ANSWER: YES, inside a few percent.** The three routes hold the analytic
von Karman law and each other to 5 percent on every quantity a receiver sees.
There is NO factor of 0.57 to 1.69 in the screens.

| quantity (route / box) | along track, 6.5 deg | rotated, 6.5 deg | along track, 45 deg | rotated, 45 deg | 2 SE |
|---|---|---|---|---|---|
| D_y(1 px) | 1.0022 | 1.0000 | 1.0021 | 0.9931 | 0.009 |
| D_x(1 px) | 0.9831 | 1.0032 | 0.9829 | 0.9999 | 0.009 |
| D_y(0.18 m) | 1.0059 | 0.9898 | 1.0059 | 0.9865 | 0.027 |
| D_x(0.18 m) | 0.9499 | 0.9855 | 0.9492 | 0.9807 | 0.027 |
| D_y(0.70 m) | 0.9936 | 0.9651 | 1.0020 | 0.9325 | 0.050 |
| D_x(0.70 m) | 0.9112 | 0.9723 | 0.9055 | 0.9671 | 0.052 |
| Z-tilt, both axes | 0.9352 | 0.9268 | 0.9493 | 0.9501 | 0.08 |
| piston-free aperture variance | 0.9529 | 0.9418 | 0.9638 | 0.9550 | 0.08 |
| RAW window variance | 0.9241 | 0.9299 | 0.9296 | 0.9179 | 0.06 |
| PSD, 0.00 to 0.05 Nyquist | 1.0043 | 1.0117 | 0.9977 | 0.9836 | 0.08 |
| PSD, 0.20 to 0.50 Nyquist | 0.9989 | 0.9980 | 1.0001 | 0.9996 | 0.004 |
| PSD, 0.80 to 1.00 Nyquist | 1.0004 | 1.0045 | 1.0004 | **0.9846** | 0.002 |

The BOX itself holds the analytic law: D(r) / law is 0.92 at the 1 px lag (the
band limit of any discrete screen) and 0.98 to 1.02 from 2 px to 0.70 m, and
its Z-tilt is 0.99 to 1.05 of the Noll filter value.

WHAT IS REAL, and what a receiver would see:

1. **The Z-tilt of a thin route is 5 to 7 percent LOW** (0.93 to 0.95 of the
   box, 2 SE 0.08). That is 0.2 to 0.3 dB of tilt POWER, not a factor. The
   gradient tilt says the same.
2. **D(r) along the LONG axis of the along-track strip reads 5 to 9 percent low
   at 0.18 to 0.70 m** (2 SE 0.03 to 0.05). Gate (a) already reported that for
   a (512, 4096) strip (x / law 0.947 at 1 m); this case has an aspect ratio of
   21, so it is a little larger. The ROTATED route splits that deficit between
   the two axes, because it turns the strip onto the lab axes.
3. **The RAW window variance of a thin route is 7 to 8 percent low** (2 SE
   0.06). That is the only quantity the short axis really moves, and it is a
   few percent, not a factor.
4. **The rotation loses 1.5 percent of the outer Nyquist ring at 45 deg**
   (0.9846 +/- 0.0018) and nothing measurable at 6.5 deg (1.0045 +/- 0.0021).
   That IS the corner-mode loss the three-shear rotation predicts: the corners
   of the frequency square leave the square under a rotation, and the loss goes
   with the angle. Every other band agrees to 0.5 percent.

WHERE THE 0.57 TO 1.69 CAME FROM. It is SAMPLING, not the strip shape. Gate 3
of `rotated_strip_gate.py` reads the scintillation index of 96 propagated
frames drawn from THREE records, so it holds three independent atmospheres, and
that gate prints no error bar on its rot/box and along-track/box columns. An
index of three atmospheres carries a standard error of tens of percent. The
per-frame ensemble here uses 64 independent seeds and it finds the screens
agree. A second reason: the arm that gate 3 labels "along track" is
`case_plans(89.999)`, a BOX whose x travel is almost zero, so its short axis is
`n` plus the seam pad (1506 px = 15 m), not the `n` = 256 px (2.56 m) of the
route of record. So that column measured neither the along-track route nor the
rotation.

THE SHARP-CUTOFF MODEL IS A LOWER BOUND, as
`validation.screens.helpers.captured_fraction` says. It predicts that a 2.56 m
screen holds only 0.42 of the Z-tilt variance of a 0.70 m pupil, and the
measured value is 0.94. The subharmonic BAND rule of the rectangular
`ScreenFactory` (see the section above) puts nearly the whole low-frequency
band back.

```
python -m validation.temporal_screens.route_equivalence
python -m validation.temporal_screens.route_equivalence --quick
```

## The rotated record against independent snapshots (step 3), 2026-09-13

`rotated_record_parity.py` runs ONE 4000-frame ROTATED record of the hero
downlink at 30 deg with a pure CROSSWIND (`wind_dir_deg = 90`,
`rotated=True`, the STEP 3 plan-wide `rot_margin = 0.5`) and 4000
INDEPENDENT snapshot trials on
the same grid, the same 9-screen plan, the `standard` preset, single
precision and the same seed. A frozen-flow record only correlates its frames
in TIME, so the two arms must hold the same one-point statistics. It is gate
(c) again, on the whole rotated record instead of frame 0 of an along-track
one.

The bar is a bootstrap of 400 resamples. Arm R takes a moving BLOCK bootstrap
of 40 frames (20 ms; Kunsch, DOI 10.1214/aos/1176347265), because its frames
are correlated. Arm S takes the ordinary bootstrap. The band is 5 percent on
an index and 0.3 dB on a mean or a fade quantile. The fade quantile is in dB
UNDER the median of its own arm, so it carries no absolute level.

| receiver | statistic | record | snapshot | comparison | 2 SE | verdict |
|---|---|---|---|---|---|---|
| point (3 cm) | sigma2_I | 0.1772 | 0.1870 | 0.9475 | 0.0996 | FAIL, inside 1 SE |
| point | mean [dB] | -32.143 | -32.145 | +0.001 | 0.130 | PASS |
| point | p5 [dB] | -3.275 | -3.178 | -0.097 | 0.217 | PASS |
| point | p1 [dB] | -5.014 | -4.659 | -0.355 | 0.468 | FAIL, inside 1 SE |
| bucket (0.7 m) | sigma2_I | 0.00474 | 0.00453 | 1.0466 | 0.1621 | PASS |
| bucket | mean [dB] | +0.017 | 0.000 | +0.017 | 0.034 | PASS |
| bucket | p5 [dB] | -0.499 | -0.495 | -0.004 | 0.055 | PASS |
| bucket | p1 [dB] | -0.696 | -0.728 | +0.031 | 0.074 | PASS |
| SMF | sigma2_I | 1.4455 | 1.5427 | 0.9370 | 0.2143 | FAIL, inside 1 SE |
| SMF | mean [dB] | -14.183 | -14.054 | -0.129 | 0.814 | PASS |
| SMF | p5 [dB] | -12.249 | -12.710 | +0.461 | 1.520 | FAIL, inside 1 SE |
| SMF | p1 [dB] | -19.579 | -19.332 | -0.247 | 3.041 | PASS |

8 of 12 bands hold, and **EVERY band that misses sits inside ONE standard
error of its own bootstrap**. The BUCKET, which has the tightest bars, holds
all four. So a rotated frozen-flow record is the same atmosphere as an
ensemble of independent snapshots. The fixed bands are a target that 4000
correlated frames of a heavy-tailed fibre statistic cannot resolve: 2 s of
record holds only about 75 independent atmospheres at the tilt time scale.

THE COST, on the bigfraw GPU (cupy), 512 px, 9 layers, 4000 frames:

| item | rotated record | along-track record (gate (e)) |
|---|---|---|
| strip pixels, 9 layers | 293.6 Mpx (1120 MB float32) | 140.2 Mpx |
| the axis-aligned BOX of the same crosswind | 3265.5 Mpx (12.5 GB) | not applicable |
| strip build | about 77 s (293.6 Mpx at the measured 3.79 Mpx/s) | 22.7 s |
| wall of the whole record | **745.4 s (0.186 s/frame)** | about 44 s (0.011 s/frame) |
| the 4000 independent snapshots, same grid | 79.1 s (0.020 s/frame) | — |
| peak working set | 2.23 GiB | — |

So the rotated route holds the crosswind in 1120 MB where the exact box needs
12.5 GB (11.1x), and it builds in about 77 s where the box needs about 860 s.
It pays for that at the FRAME: 0.186 s against 0.011 s, which is 17x, because
each frame crops a 1234 px patch of every layer, uploads it and turns it with
three shears. The crop side is ONE number for the whole plan and the widest
footprint of the plan sets it (a 40 deg layer here), so a PER-LAYER crop side
and a margin rule by ANGLE are the obvious cuts (STEP 4 below does both). A
crosswind record is
therefore about 12 minutes against about 15 minutes for the box build alone,
and it needs a tenth of the disk.

```
python -m validation.temporal_screens.rotated_record_parity --smoke
python -m validation.temporal_screens.rotated_record_parity
python -m validation.temporal_screens.rotated_record_parity --analyse
```

## Step 4: the taper, the per-layer crop and the device strips, 2026-09-13

Step 3 left the rotated route at 0.186 s for one frame against 0.011 s for an
along-track frame. TWO causes were named: the crop side was ONE number for the
whole plan (1234 px, set by the widest footprint), and every frame uploaded
about 55 MB of host crops to the device. Step 4 removes both, and it adds the
taper that makes a SMALL crop safe.

### The taper (`rotation_taper.py`)

A shear reads each line as PERIODIC, so the STEP between the two ends of a
crop rings inward with a sinc tail. Six cures were measured against the
rotation of the FULL periodic screen (1024 px, `n = 256`, `r0 = 10 cm`,
`L0 = 25 m`, 4 seeds, 6.5 and 45 deg). THE WINNER IS V4: a 1-D raised-cosine
taper on the ends of every row or column, applied BEFORE EACH of the three
shear transforms. A single 2-D window (V2, V3) does NOT work, because the two
later shears undo it. The flat top of the taper is the rotation footprint
`n(|cos t| + |sin t|)`, so the kept frame never reads a tapered pixel.

| variant | 6.5 deg, edge ring | 45 deg, edge ring | smallest margin that passes |
|---|---|---|---|
| V0 plain (the control) | 1.79e-2 rad at 0.5 n | 1.89e-2 rad at 0.5 n | none of the tested |
| V4 per-shear taper | 5.0e-3 rad at 0.5 n | 6.6e-3 rad at 0.5 n | 0.125 n and 0.25 n |

`rotate_fourier(patch, theta, taper=flat_px)` is that taper, and `taper=None`
(the default) keeps the plain call, which the validation scripts use as the
control. The FRAME route always passes the footprint.

### The crop side of a layer

`TemporalSpec.rot_margin` now means the ROLL-OFF PAST the geometric rotation
footprint, on each side, and the crop side of LAYER j is

    m_j = n (|cos t_j| + |sin t_j|) + 2 * rot_margin * n.

The default is 0.07. A sweep of the roll-off at `n = 256`, 4 seeds, against
the same reference (the gate is an outer-16 px ring rms under 1e-2 rad and
every `|D(r) ratio - 1|` under 1 percent at 1, 4 and 16 px):

| angle | roll-off 0.04 | roll-off 0.07 | roll-off 0.10 |
|---|---|---|---|
| 3 deg | 7.97e-3 PASS | 7.37e-3 PASS | 7.03e-3 PASS |
| 6.5 deg | 7.03e-3 PASS | 6.44e-3 PASS | 6.10e-3 PASS |
| 20 deg | 6.65e-3 PASS | 6.41e-3 PASS | 6.22e-3 PASS |
| 45 deg | 7.41e-3 PASS | 7.25e-3 PASS | 7.11e-3 PASS |

The ring error barely falls past 0.04, because what is left is the CORNER-mode
loss of the rotation, not the seam. So 0.07 carries a safety factor of about
two on the roll-off and it stays cheap. A layer that walks nearly along x now
asks for a crop of about 1.14 n where the old plan-wide 0.5 n margin asked
for 2.41 n, and the cost of a frame goes as the crop AREA.

### The strips stay on the device

`open_strips(paths, on_device=True)` uploads the strips ONE time, and the
runner does that for the ROTATED route under the `"cupy"` backend. A frame
then crops and turns them where the field is and it uploads nothing. The call
falls back to the host memory map, with a warning, when the strips need more
than half the free device memory. The along-track route keeps the memory map,
unchanged: its frame is a plain slice, so it has nothing to gain.

### The cost and the parity, the SAME hero record

The 4000-frame hero record of the section above, rerun on the bigfraw GPU
(cupy, 512 px, 9 layers, 30 deg, `wind_dir_deg = 90`, `standard`, single
precision, seed 20260913). The wall of each arm holds the strip build.

| item | step 3 | step 4 | change |
|---|---|---|---|
| crop side | 1234 px, ONE for the plan | 614 to 794 px, one for EACH layer | 0.25 to 0.41 of the area |
| strips on disk | 1120 MB | 643 MB | 1.74x smaller |
| peak working set | 2.23 GiB | 1.55 GiB | 1.44x smaller |
| wall, 4000 frames | 745.4 s (**0.186 s/frame**) | 217.2 s (**0.0543 s/frame**) | **3.43x faster** |
| an ALONG-TRACK record, same grid | 0.011 s/frame | 0.011 s/frame | the rotation now costs 4.9x, not 17x |
| 4000 INDEPENDENT snapshots | 0.020 s/frame | 0.020 s/frame | the rotation now costs 2.7x, not 9.3x |
| pass bands | 8 of 12 | **9 of 12** | — |

The per-layer crop sides of the hero, with the resultant angle of each layer:

| layer | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|---|---|---|
| angle [deg] | 7 | 15 | 22 | 28 | 32 | 31 | 40 | 69 | 87 |
| crop side [px] | 640 | 696 | 740 | 766 | 778 | 776 | 794 | 736 | 614 |

THE PARITY HOLDS, and it is BETTER than step 3. The POINT receiver now holds
all four bands (it held two), the BUCKET holds all four again, and the SMF
holds one.

| receiver | statistic | record | snapshot | comparison | 2 SE | verdict |
|---|---|---|---|---|---|---|
| point (3 cm) | sigma2_I | 0.1802 | 0.1870 | 0.9637 | 0.0944 | PASS |
| point | mean [dB] | -32.055 | -32.145 | +0.090 | 0.123 | PASS |
| point | p5 [dB] | -3.144 | -3.178 | +0.034 | 0.261 | PASS |
| point | p1 [dB] | -4.838 | -4.659 | -0.179 | 0.590 | PASS |
| bucket (0.7 m) | sigma2_I | 0.00432 | 0.00453 | 0.9538 | 0.1175 | PASS |
| bucket | mean [dB] | +0.006 | 0.000 | +0.006 | 0.031 | PASS |
| bucket | p5 [dB] | -0.479 | -0.495 | +0.016 | 0.055 | PASS |
| bucket | p1 [dB] | -0.765 | -0.728 | -0.037 | 0.104 | PASS |
| SMF | sigma2_I | 2.128 | 1.543 | 1.3796 | 0.3660 | FAIL, about 1.8 SE |
| SMF | mean [dB] | -14.967 | -14.054 | -0.913 | 1.041 | FAIL, inside 1 SE |
| SMF | p5 [dB] | -12.261 | -12.710 | +0.449 | 1.361 | FAIL, inside 1 SE |
| SMF | p1 [dB] | -19.452 | -19.332 | -0.120 | 2.784 | PASS |

THE CAVEAT. The smaller crop also makes the strip SHORT AXIS smaller: 8.5 m
before, 4.2 to 5.4 m now. That is the range `route_equivalence.py` already
measured against the exact box (D(r), the tilt, the piston-free variance and
every radial spectrum band inside a few percent), and the ALONG-TRACK route of
record runs at 3.5 m, so the new sides sit between the two validated routes.
The SMF scintillation INDEX is the noisiest statistic of the set: 2 s of
record holds only about 75 independent atmospheres at the tilt time scale, and
the index of a heavy-tailed fibre series is a fourth-moment estimator. Its
mean moved 0.91 dB, which sits inside its own bootstrap bar.

## Sources

- Taylor, DOI 10.1098/rspa.1938.0032. The frozen flow hypothesis.
- Andrews and Phillips, DOI 10.1117/3.626196, Ch. 12, Eqs. (2) and (3),
  printed p. 481. The Bufton wind with the slew term.
- Schmidt, DOI 10.1117/3.866274, Ch. 9, Eqs. (9.78) to (9.81), printed
  pp. 166 to 169. The Fourier screen and its subharmonics.
- Lane, Glindemann and Dainty, DOI 10.1088/0959-7174/2/3/003. The subharmonic
  method.
- Assemat and Wilson, DOI 10.1364/OE.14.000988, Eq. (5). The von Karman
  covariance that bounds the seam, and the extrusion that this route replaces.
- Noll, DOI 10.1364/JOSA.66.000207. The Zernike tilt.
- Greenwood, DOI 10.1364/JOSA.67.000390. The Greenwood frequency.
- Tyler, DOI 10.1364/JOSAA.11.000358. The tilt power spectrum and its V/D
  corner.
- Unser, Thevenaz and Yaroslavsky, DOI 10.1109/83.469963. The three-shear
  rotation of the rotated thin strip.
- Kunsch, Ann. Statist. 17(3), pp. 1217 to 1241 (1989),
  DOI 10.1214/aos/1176347265. The moving block bootstrap of a correlated
  series, which `rotated_record_parity.py` uses on the record arm.
- Andrews and Phillips, DOI 10.1117/3.626196, Ch. 14, Eq. (38), printed p. 622.
  The Greenwood frequency of a slant path, which `greenwood_frequency` gives.

## The low-frequency band of a strip

`figures/strip_band.png` shows the frequency plane near zero. A square
grid has one missing cell at f = 0, and the Lane subharmonics fill it. A
strip has an fx spacing of 1/Lx, much finer than its fy spacing of 1/Ly, so
the fy = 0 row is a wide slab of thin cells, each sampled at fy = 0 only.
Per-axis subharmonics fill the origin box only, and D(r) along y reads 8 to
23 percent low. The landed rule zeroes that row in the main filter and
redraws the slab as 2P+1 sub-rows in fy (the 1-D Lane partition), each at
every fx. See the `ScreenFactory` docstring and gate (a).
