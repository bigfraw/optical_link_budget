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
| `hero_temporal.py` | (e) | the fade RATE and the fade DURATION at the 5 percent level, from 8 records of 2 s at 30 and 20 deg. SMOKE DONE 2026-09-13 (3.6 s); the full run is bigfraw time |

```
python -m validation.temporal_screens.rect_factory
python -m validation.temporal_screens.strip_taylor
python -m validation.temporal_screens.frame0_parity
python -m validation.temporal_screens.tilt_spectrum
python -m validation.temporal_screens.hero_temporal --smoke
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
- Andrews and Phillips, DOI 10.1117/3.626196, Ch. 14, Eq. (38), printed p. 622.
  The Greenwood frequency of a slant path, which `greenwood_frequency` gives.
