"""GPU FFT microbenchmark for the fidelity-2 trial (backlog 2-N8).

The script measures three things on one machine, and it writes them to
JSON:

1. RAW TRANSFORMS. One complex64 and one complex128 fft2 at 512, 1024,
   2048 and 4096 px: numpy, scipy (workers=1, in place) and cupy (cuFFT).
   The cupy time includes the device synchronise, not the upload.

2. THE REAL CPU TRIAL. One serial trial of the 30 deg hero downlink at a
   pinned pixel count, single precision, the standard preset, L0 = 25 m,
   with the numpy backend and with the scipy backend. A counter on the
   Forvard transforms gives the number of FFTs the plan makes, so the
   synthetic GPU trial below replays the SAME work.

3. THE SYNTHETIC GPU TRIAL. On the device, the same sequence of
   operations a trial makes: for each hop the sign multiply, fft2, the
   transfer-function multiply, ifft2, the piston multiply, the sign
   multiply and the boundary mask; for each screen the exp(i phi)
   multiply; for each screen ONE filtered-white-noise ifft2 (the real part
   kept, as ScreenFactory.make does) plus the three subharmonic levels of
   nine outer products; one download of the receive field. The white noise
   is drawn three ways: on the HOST by numpy PCG64, serial, and uploaded
   (the CPU random stream kept); on the host in one thread per screen,
   overlapped with the GPU work of the previous trial (the stream kept, the
   draw hidden; owner decision 2026-09-06, the plan of record); on the
   DEVICE by cuRAND (a different stream). The estimate is a lower bound on
   the real GPU trial, because it has no Python glue.

Run on the machine with the CUDA device:

    python validation/gpu_fft/gpu_microbench.py --n 1024 2048 --trials 3

Documentation uses ASD-STE100 Simplified Technical English.
"""

import argparse
import json
import os
import platform
import sys
import time

import numpy as np
from scipy import fft as sfft

from olb.geometry import CircularOrbit
from olb.scenario import Channel, Site, SpaceScenario
from olb.terminal import SMF, Terminal, Transmitter
from olb.waveoptics import propagators as P
from olb.waveoptics.grid import GridSpec
from olb.waveoptics.priority import boost_process_priority
from olb.waveoptics.turbulence.run import propagate_turbulent_scenario
from olb.waveoptics.turbulence.sampling import turbulent_grid

try:
    import cupy as cp
except ImportError:                       # the CPU parts still run
    cp = None


def hero_scenario():
    """The 30 deg hero downlink: 0.7 m SMF ground terminal, 500 km."""
    ground = Terminal(aperture_m=0.7, obscuration_ratio=0.3,
                      wavelength_m=1550e-9, pointing_jitter_rad=1e-6,
                      detector=SMF(sensitivity_dbm=-40.0, focal_length_m=2.0,
                                   mode_field_radius_m=5.2e-6))
    space = Terminal(aperture_m=0.1, obscuration_ratio=0.0,
                     wavelength_m=1550e-9, pointing_jitter_rad=1e-6,
                     transmitter=Transmitter(waist_m=0.03, power_dbm=30.0,
                                             divergence_rad=None))
    scn = SpaceScenario(ground=ground, space=space,
                        channel=Channel(site=Site(), altitude_m=500e3),
                        direction="downlink")
    return scn, CircularOrbit(altitude_m=500e3, elevation_deg=30.0)


def pinned(scn, geo, preset, n, L0_m=25.0):
    """Size with the preset, then pin the pixel count only."""
    g, plan, _ = turbulent_grid(scn, geo, preset=preset, L0_m=L0_m)
    return GridSpec(size_m=g.size_m, n=int(n), scaled=g.scaled), plan


def best_of(fn, reps):
    """The best wall time of `reps` calls, in s."""
    best = np.inf
    for _ in range(reps):
        t0 = time.perf_counter()
        fn()
        best = min(best, time.perf_counter() - t0)
    return best


# ---------------------------------------------------------------- part 1

def raw_transforms(sizes, reps):
    out = {}
    for n in sizes:
        row = {}
        for cdt, tag in ((np.complex64, "c64"), (np.complex128, "c128")):
            a = (np.random.default_rng(0).standard_normal((n, n))
                 + 0j).astype(cdt)
            row[f"numpy_{tag}_s"] = best_of(lambda: np.fft.fft2(a), reps)
            b = a.copy()
            row[f"scipy_{tag}_s"] = best_of(
                lambda: sfft.fft2(b, overwrite_x=True, workers=1), reps)
            if cp is not None:
                d = cp.asarray(a)
                cp.fft.fft2(d)                     # plan and warm up
                cp.cuda.Device().synchronize()

                def gpu():
                    cp.fft.fft2(d)
                    cp.cuda.Device().synchronize()
                row[f"cupy_{tag}_s"] = best_of(gpu, reps)
        out[str(n)] = row
        print(f"  fft2 {n:5d} px: " + "  ".join(
            f"{k[:-2]} {v * 1e3:8.3f} ms" for k, v in row.items()))
    return out


# ---------------------------------------------------------------- part 2

class FftCounter:
    """Count the Forvard transforms of one run through the backend hooks."""

    def __init__(self):
        self.n = 0

    def __enter__(self):
        self._f, self._i = P._fft2, P._ifft2

        def f(a):
            self.n += 1
            return self._f(a)

        def i(a):
            self.n += 1
            return self._i(a)
        P._fft2, P._ifft2 = f, i
        return self

    def __exit__(self, *exc):
        P._fft2, P._ifft2 = self._f, self._i


def cpu_trial(scn, geo, grid, plan, trials, backend):
    kw = dict(seed=7, preset="standard", threader=None, precision="single",
              L0_m=25.0, grid=grid, plan=plan, fft_backend=backend)
    with FftCounter() as c:
        propagate_turbulent_scenario(scn, geo, n_trials=1, **kw)   # warm up
    n_fft = c.n
    t0 = time.perf_counter()
    res = propagate_turbulent_scenario(scn, geo, n_trials=trials, **kw)
    dt = (time.perf_counter() - t0) / trials
    return dt, n_fft, float(res.trials[0].collected_power)


# ---------------------------------------------------------------- part 3

def gpu_trial(n, n_fft, n_screens, reps, n_sub_levels=3, noise="host"):
    """Replay the trial's operations on the device. Gives s per trial and a
    breakdown: propagation, screens, download.

    noise: "host" draws the white noise serially with numpy PCG64 and
    uploads it; "threads" draws it in one thread per screen for the NEXT
    trial while the GPU runs this one (numpy releases the GIL on a big
    draw, and each screen owns its seed, so the order does not matter);
    "device" draws it with cuRAND.
    """
    if cp is None:
        return None
    from concurrent.futures import ThreadPoolExecutor
    rng = np.random.default_rng(0)
    cdt, rdt = cp.complex64, cp.float32
    n_hops = n_fft // 2
    field = cp.asarray((rng.standard_normal((n, n)) + 0j).astype(np.complex64))
    iiij = cp.asarray(rng.choice([-1.0, 1.0], (n, n)).astype(np.float32))
    CC = cp.exp(1j * cp.asarray(rng.standard_normal((n, n)).astype(np.float32)))
    CC = CC.astype(cdt)
    mask = cp.asarray(rng.random((n, n)).astype(np.float32))
    filt = cp.asarray(rng.random((n, n)).astype(np.float32))
    piston = cp.complex64(np.exp(1j * 0.3))
    xs = cp.asarray(np.linspace(-1, 1, n).astype(np.float32))
    ang = cp.exp(1j * cp.asarray(rng.random((n_sub_levels, 9, 2)) * 6.283,
                                dtype=cp.float32))
    pool = (ThreadPoolExecutor(max_workers=n_screens + 1)   # +1: the outer submit
            if noise == "threads" else None)

    def host_draw(seed):
        """One screen's coefficient grid, as ScreenFactory._base_pair draws
        it: two n x n double normals, cast to complex64."""
        r = np.random.default_rng(seed)
        g = r.standard_normal((n, n)) + 1j * r.standard_normal((n, n))
        return g.astype(np.complex64)

    def draw_all(trial):
        seeds = [(trial, j) for j in range(n_screens)]
        if noise == "device":
            return None
        if noise == "threads":
            return list(pool.map(host_draw, [hash(sd) & 0xffffffff for sd in seeds]))
        return [host_draw(hash(sd) & 0xffffffff) for sd in seeds]

    def one_screen(cn_host):
        if cn_host is None:
            cn = (cp.random.standard_normal((n, n), dtype=rdt)
                  + 1j * cp.random.standard_normal((n, n), dtype=rdt))
        else:
            cn = cp.asarray(cn_host)
        cn *= filt
        s1 = cp.fft.ifft2(cn).real.copy()
        for lev in range(n_sub_levels):
            for q in range(9):
                ph = cp.exp(1j * 0.01 * (lev + 1) * xs)
                s1 += (ph[:, None] * ph[None, :] * ang[lev, q, 0]).real
        return s1

    def gpu_part(noise_grids):
        """The device work of one trial, given the host noise (or None)."""
        f = field.copy()
        grids = noise_grids if noise_grids is not None else [None] * n_screens
        screens = [one_screen(g) for g in grids]
        every = max(1, n_hops // max(1, n_screens))
        for h in range(n_hops):
            f *= iiij
            f = cp.fft.fft2(f)
            f *= CC
            f = cp.fft.ifft2(f)
            f *= piston
            f *= iiij
            f *= mask
            if h % every == 0 and screens:
                f *= cp.exp(1j * screens.pop())
        cp.cuda.Device().synchronize()
        return f

    # warm up: plans, pool threads
    gpu_part(draw_all(-1))

    # the pipelined loop: draw trial k+1 (threads) while the GPU does trial k
    def run_pipeline(k_trials):
        if noise == "threads":
            nxt = pool.submit(draw_all, 0)
            for k in range(k_trials):
                grids = nxt.result()
                nxt = pool.submit(draw_all, k + 1)
                f = gpu_part(grids)
            nxt.result()
        else:
            for k in range(k_trials):
                f = gpu_part(draw_all(k))
        return f

    t = best_of(lambda: run_pipeline(reps), 1) / reps
    t_draw = (0.0 if noise == "device" else
              best_of(lambda: draw_all(0), reps))
    t_gpu = best_of(lambda: gpu_part(draw_all(0)), reps)
    f = gpu_part(draw_all(0))
    t_download = best_of(lambda: f.get(), reps)
    if pool is not None:
        pool.shutdown()
    return {"s_per_trial": t, "host_draw_s": t_draw, "gpu_only_s": t_gpu,
            "download_s": t_download, "n_hops": n_hops,
            "n_screens": n_screens, "noise": noise}


# ---------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--n", type=int, nargs="+", default=[1024, 2048])
    ap.add_argument("--trials", type=int, default=3)
    ap.add_argument("--reps", type=int, default=5)
    ap.add_argument("--out", default=os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "gpu_microbench_results.json"))
    args = ap.parse_args()
    boost_process_priority()

    out = {"environment": {
        "python": sys.version.split()[0], "numpy": np.__version__,
        "cupy": None if cp is None else cp.__version__,
        "cpu": platform.processor(), "cores": os.cpu_count(),
        "gpu": (None if cp is None else
                cp.cuda.runtime.getDeviceProperties(0)["name"].decode()),
        "platform": platform.platform()}}
    print("environment:", json.dumps(out["environment"]))

    print("\n1. raw transforms (best of reps)")
    out["raw"] = raw_transforms([512, 1024, 2048, 4096], args.reps)

    scn, geo = hero_scenario()
    out["cases"] = []
    for n in args.n:
        grid, plan = pinned(scn, geo, "standard", n)
        n_screens = int(plan.z_m.size)
        print(f"\n2. real CPU trial at {n} px, {n_screens} screens, single")
        row = {"n": n, "n_screens": n_screens, "cpu": {}}
        for backend in ("numpy", "scipy"):
            dt, n_fft, cpow = cpu_trial(scn, geo, grid, plan, args.trials, backend)
            row["cpu"][backend] = {"s_per_trial": dt, "n_forvard_ffts": n_fft,
                                   "collected_power": cpow}
            print(f"   {backend:6s}: {dt:7.3f} s/trial, {n_fft} Forvard FFTs")
        n_fft = row["cpu"]["numpy"]["n_forvard_ffts"]
        print(f"3. synthetic GPU trial at {n} px ({n_fft} FFTs, {n_screens} screens)")
        for noise in ("host", "threads", "device"):
            g = gpu_trial(n, n_fft, n_screens, args.reps, noise=noise)
            if g is None:
                print("   no cupy")
                break
            row[f"gpu_noise_{noise}"] = g
            print(f"   noise {noise:7s}: {g['s_per_trial'] * 1e3:8.1f} ms/trial "
                  f"(gpu work {g['gpu_only_s'] * 1e3:.1f}, host draw "
                  f"{g['host_draw_s'] * 1e3:.1f}, download "
                  f"{g['download_s'] * 1e3:.1f} ms)")
            row[f"speedup_{noise}_vs_numpy_serial"] = (
                row["cpu"]["numpy"]["s_per_trial"] / g["s_per_trial"])
            row[f"speedup_{noise}_vs_scipy_serial"] = (
                row["cpu"]["scipy"]["s_per_trial"] / g["s_per_trial"])
            print(f"      over one serial CPU trial: "
                  f"{row[f'speedup_{noise}_vs_numpy_serial']:.0f}x numpy, "
                  f"{row[f'speedup_{noise}_vs_scipy_serial']:.0f}x scipy")
        out["cases"].append(row)

    with open(args.out, "w") as fh:
        json.dump(out, fh, indent=2)
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
