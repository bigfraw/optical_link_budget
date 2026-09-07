# GPU FFT microbenchmark (backlog 2-N8)

A quick look at what a CUDA device gives the fidelity-2 split-step trial,
before any olb code changes. The script is `gpu_microbench.py`; the record
is `gpu_microbench_results.json`.

## Why

A trial is FFTs: 40 Forvard transforms at 1024 px (60 at 2048 px) plus one
screen transform for each of the 9 screens, and the 12-worker process pool
on bigfraw is memory-bandwidth bound (see
`validation/terrestrial_screen_count/`). A GPU has several times the memory
bandwidth of the host DDR5, so it is the one lever left on the plateau.

## What the script measures

1. **Raw transforms.** One fft2 at 512 to 4096 px, complex64 and
   complex128: numpy, scipy (`workers=1`, in place) and cupy (cuFFT).
2. **The real CPU trial.** One serial trial of the 30 deg hero downlink
   (0.7 m SMF, 500 km) at a pinned pixel count, single precision, the
   `standard` preset, L0 = 25 m, with the numpy and the scipy backends. A
   counter on the Forvard hooks gives the FFT count of the plan.
3. **A synthetic GPU trial.** The same operations replayed on the device:
   the hops (sign, fft2, transfer function, ifft2, piston, sign, mask), the
   screen multiplies, and one filtered-white-noise ifft2 plus three
   subharmonic levels for each screen, then one download of the field. It
   has no Python glue, so it is a lower bound on a real GPU trial.

The white noise of the screens is drawn three ways, because that draw is
the design question (see below): `host` (numpy PCG64, serial, uploaded),
`threads` (numpy PCG64, one thread per screen, drawn for the next trial
while the GPU runs this one) and `device` (cuRAND).

## The random-stream question

`ScreenFactory.make` draws one complex n x n grid of double normals for
each screen from numpy PCG64, seeded by (entropy, trial, screen). The whole
campaign design hangs off that seed chain: a block is a bit-identical slice
of one run, a resume never re-draws, the seed rebuilds any screen, and an
opt-in is certified by one relative-rms number against the default path.

That draw is 19 million normals for a 1024 px trial and 75 million at
2048 px, about 0.2 s and 0.9 s serial on one core, which is 3 to 5 times
the GPU work it feeds. So a GPU trial that keeps the CPU stream and draws
serially gains nothing at 2048 px.

DECISION (owner, 2026-09-06): keep the CPU stream and hide the draw. Each
screen owns its seed, so the order of the draws does not matter, and numpy
releases the GIL on a big draw, so one thread per screen draws the next
trial while the GPU runs this one. A device draw (cuRAND) is the fallback
if the threaded draw still outpaces the GPU; it is a different atmosphere
for the same seed, so it gives up the trial-for-trial cross-check. A
float32 host draw is out: it consumes the generator differently and changes
the stream.

## The first run (2026-09-06, CONTENDED)

The run on bigfraw (RTX 4070 Laptop, 8 GB, cupy 14.2) happened while the
terrestrial screen-count sims held the CPU, so the CPU numbers are high
and the CPU-to-GPU ratios overstate the win. The GPU was idle, so its
numbers stand. That run replayed 5 noise grids, not 9; the script now
draws 9.

| complex64 fft2 | numpy (contended) | scipy (contended) | cupy |
|---|---|---|---|
| 1024 px | 100 ms | 14.3 ms | 0.12 ms |
| 2048 px | 458 ms | 80 ms | 0.53 ms |

| synthetic GPU trial | 1024 px | 2048 px |
|---|---|---|
| propagation | 18 ms | 96 ms |
| screens, device noise | 22 ms | 200 ms |
| screens, host noise (5 grids, contended) | 340 ms | 1255 ms |
| download | 4 ms | 16 ms |

The complex128 cuFFT at 2048 px read 41 ms, slower than 4096 px at 17 ms.
Recheck; double is not the target (the 4070 is FP64-starved).

## The clean run (2026-09-07, pool idle)

Same machine, no sims running, 9 noise grids per trial. The CPU baselines
are the serial single-precision trial of the hero downlink.

| complex64 fft2 | numpy | scipy | cupy | cupy / scipy |
|---|---|---|---|---|
| 1024 px | 30.6 ms | 4.2 ms | 0.074 ms | 57x |
| 2048 px | 154 ms | 22.0 ms | 0.61 ms | 36x |

| one trial | 1024 px (40 FFTs) | 2048 px (60 FFTs) |
|---|---|---|
| CPU serial, numpy backend | 1.94 s | 12.3 s |
| CPU serial, scipy backend | 1.52 s | 9.0 s |
| GPU, host noise serial | 310 ms | 1173 ms |
| GPU, host noise threaded (plan of record) | 66 ms | 349 ms |
| GPU, device noise (cuRAND) | 21 ms | 284 ms |
| host draw alone, 9 threads | 40 ms | 164 ms |
| field download | 1.5 ms | 5.5 ms |

The threaded host draw HIDES the draw at 2048 px (349 against 284 ms, a
1.2x cost for keeping the CPU stream) and only PARTLY at 1024 px (66
against 21 ms: the GPU work is 20 ms, shorter than the 40 ms threaded
draw). Against the 12-worker pool of record (1.10 s/trial at 2048 px,
`validation/terrestrial_screen_count/`), one GPU stream with the threaded
draw is about 3x at 2048 px; at 1024 px about 4x, from the serial number
and the measured pool plateau.

The "gpu work" column of the script's `host` and `threads` rows includes
the draw, because the un-pipelined call draws inside it; the
`s_per_trial` row is the pipelined number and the one to read.

Odd numbers to know: the complex64 cuFFT at 4096 px read 41.5 ms, slower
than complex128 at 16.7 ms (a cuFFT plan quirk at that size; not the
target grid), and the earlier complex128 2048 px anomaly did NOT repeat
(3.9 ms).

VERDICT: the plan of record holds. Build the `fft_backend="cupy"` opt-in
with the threaded host draw; the cuRAND fallback is not needed at the
grids in use.

## The run command

Run it with the pool idle:

    ssh desktop "cd D:\repos\optical_link_budget; C:\Users\alexf\olb-gpu-venv\Scripts\python.exe validation\gpu_fft\gpu_microbench.py --n 1024 2048 --trials 3"

The venv `C:\Users\alexf\olb-gpu-venv` on bigfraw is a
`--system-site-packages` venv over the `olb` conda env with `cupy-cuda12x`
and the nvidia CUDA wheels (`nvidia-cufft-cu12`, `nvidia-curand-cu12`,
`nvidia-cuda-runtime-cu12`, `nvidia-cuda-nvrtc-cu12`, `nvidia-cublas-cu12`)
added. The `olb` env is untouched. cupy 14 does not bundle the CUDA
libraries, so the nvidia wheels are required; no CUDA toolkit is needed.

## The design that follows, if the clean run holds

A third `fft_backend="cupy"` next to `"numpy"` and `"scipy"`, an opt-in
and a different fingerprint. `propagators.py` (`_fft2`/`_ifft2`, the
Forvard factor cache, `xp.array` not `np.array`), `splitstep.py` (the mask
uploaded once) and `screens.py` (the filter and the subharmonics through
`xp`, the rng on the host) move to the device; the runner downloads the
receive field once after `split_step`, and the clip, the coupling and the
patch store stay CPU code. A GPU campaign runs one process, not the pool.
