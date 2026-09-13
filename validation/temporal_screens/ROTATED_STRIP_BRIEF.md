# Brief: the rotated strip for a crosswind (frozen-flow time axis)

A handoff for a code agent. Repo `optical_link_budget`, branch
`waveoptics-temporal`. Backlog item 2-P1b (10). Write all code and prose in
ASD-STE100 Simplified Technical English. Cite every equation by DOI. Follow the
ponytail rule: the laziest solution that works, reuse the kernels, no
speculative abstraction.

## The context

The fidelity-2 time axis makes one oversized Fourier "strip" screen for each
turbulence layer, then it scrolls an `n x n` crop window across the strip, one
frame at each step (`olb/waveoptics/turbulence/temporal.py`). The velocity of a
layer is a 2-D vector: the satellite slew along `+x` plus the Bufton wind at
`wind_dir_deg`.

    v_x = omega_slew * altitude + v_buf * cos(wind_dir)
    v_y =                          v_buf * sin(wind_dir)

Today the strip is the AXIS-ALIGNED BOUNDING BOX of that 2-D walk:

    nx = n + ceil(|v_x| * T / dx) + pad
    ny = n + ceil(|v_y| * T / dx) + pad   (the +pad only when v_y != 0)

`pad = ceil(pad_outer_scales * L0 / dx)`, default `pad_outer_scales = 2.0`
(the torus seam pad). `nx` rounds up to a multiple of `STRIP_ALIGN = 32`.

So a CROSSWIND grows the SHORT axis and the box balloons. The default
(`wind_dir_deg = 0`, along track) keeps the short axis at `n + pad` and is
cheap.

## The problem, measured (bigfraw, 30 deg hero, standard preset)

Generation of the box is the cost, not the memory or the per-frame crop.

| case | layers built | size | generation time |
|---|---|---|---|
| along-track (thin), `wind_dir_deg = 0` | all 9 | 0.52 GiB | 22.7 s |
| perp crosswind (box), `wind_dir_deg = 90` | 1 of 9 | 1.94 GiB | 170 s |

One perpendicular box layer (11897 x 43744) took 170 s and it maxed the host
RAM before layer 1. The full 9-box build is about 25 min plus the blow-up. Two
causes: the box FFT is about 23x the area of the thin strip, and the
`ScreenFactory` filter builds float64 `meshgrid` arrays of the full 2-D size
(about 3.9 GiB each at that shape). Strips build on the HOST (scipy backend),
one layer at a time, then memory-map.

The per-frame crop is cheap: a GPU bilinear rotate of a 512 px frame is 0.071 ms
(0.28 s per layer per 4000-frame record); cubic is 0.390 ms (1.56 s). The
propagation is about 40 s per record. So the rotation route wins on GENERATION
time and RAM. That is the case for a crosswind.

## Two routes to build

- **(a) rotate once onto a big grid.** FFT the thin strip, rotate it ONE time
  onto a big axis-aligned grid, then slice with integer offsets as now. Saves
  the big FFT, keeps the big memory, one interpolation pass over the whole
  screen.
- **(b) thin strip, rotate each frame.** Hold the thin strip along the
  resultant velocity, and rotate each frame's `n x n` crop. Saves the FFT AND
  the memory, one interpolation per crop. Use `cupyx.scipy.ndimage`.

Both routes need a ROTATION-CORNER PAD. A grid-aligned `n x n` crop rotated by
theta has a strip footprint of `n(|cos theta| + |sin theta|)`, so the SHORT
axis needs `n(cos theta + sin theta) + seam_pad`. That is an extra about
`0.41 n` at 45 deg, and it is negligible at the small angles of the fast
layers. The along-track default needs NONE of this, because its box IS the thin
strip.

## The open question: fidelity, not speed

Interpolating a phase screen SMOOTHS its high-frequency (Fresnel-scale)
structure. That structure builds the scintillation, and it is the reason the
coarse screen died in P2 (`validation/waveoptics_speed/`). So a rotated route
MUST pass a gate before it is trusted:

1. The phase structure function `D(r)` of a rotated crop against the
   bounding-box strip and the von Karman law, inside a few percent on BOTH
   axes (reuse the estimators in `validation/screens/helpers.py` and the
   gate-(a) script `validation/temporal_screens/rect_factory.py`).
2. The aperture scintillation index `sigma2_I` of a short end-to-end record
   built by the rotated route against the bounding-box route, inside 5 percent
   (the P2 kill line). Reuse `validation/temporal_screens/frame0_parity.py` for
   the snapshot side.

Bilinear (`order=1`) likely fails these; cubic (`order=3`) is safer but still
smooths. Measure both. A route that fails the gate is dead, whatever it saves.

## The task

1. Add a `wind_dir_deg` sweep to a scratch script and confirm the generation
   numbers above (the along-track plan is the reference).
2. Prototype route (b) first (it saves the most): a rotated `frame_stack`
   variant that holds one thin strip per layer along the resultant and rotates
   each crop on the device. Add the corner pad to `strip_plan`.
3. Run the two fidelity gates. Report `D(r)` and `sigma2_I` for bilinear and
   cubic against the bounding-box route.
4. If a route passes, wire it behind an OPT-IN flag on `TemporalSpec` (default
   OFF, so the along-track path is bit-identical), thread it through the runner
   and `Campaign`, and enter it in the fingerprint tail only when non-default.
5. If both routes fail the gate, record the fidelity numbers and close the item
   as NO. A wrong screen is worse than a slow one.

## A separate, smaller lever

The `ScreenFactory` filter uses float64 `meshgrid` arrays of the full 2-D size.
On any large strip (crosswind or not) that wastes RAM and time. A broadcast
float32 build would cut the peak RAM and the generation time for EVERY strip.
It is bit-changing at the rounding level, so it is its own opt-in study, not
part of the rotation work.

## How to run

- Python on this laptop: `C:\Users\alexf\anaconda3\envs\olb\python.exe`, run
  modules as `python -m ...` from the repo root.
- Heavy or GPU work goes to bigfraw: `ssh desktop` (the shell is PowerShell
  5.1, so chain with `;`, never `&&`). The GPU venv is
  `C:\Users\alexf\olb-gpu-venv\Scripts\python.exe`. Launch a long run detached
  through WMI `Win32_Process` Create, not `Start-Process`. Pull data back as
  ONE zip (`Compress-Archive` then `scp`), never `scp -r`.
- The temporal API: `strip_plan(plan, grid, spec, geometry, L0_m)`,
  `frame_offsets(sp, k)`, `frame_stack(strips, sp, k)`,
  `build_strips(sp, r0_m, paths, L0_m, subharmonics, seeds, dtype)`,
  `open_strips(paths)`, and `ScreenFactory(ny, dx, nx=nx, L0_m=..., dtype=...)`
  for the rectangular strip.
