"""The certification of the cupy FFT backend (backlog 2-N8, milestone 1).

The script measures the "cupy" FFT backend against the "numpy" backend of
record, on the 30 deg hero downlink of `gpu_microbench.py`. It answers three
questions:

(a) THE ATMOSPHERE. Does `ScreenFactory.make` give the SAME screen on the
    device and on the host for the same seed? The white noise is a numpy
    PCG64 draw on the host in both cases, and only the filter multiply and
    the transform move, so the two screens must agree at the rounding level
    of float32 (about 1e-7 relative rms).

(b) THE PROPAGATION. Does `split_step` give the same receive field on the
    device and on the host, for the same start field and the SAME screens?
    The script gives the host screens to both runs, so the difference is the
    transforms only. The expectation is about 1e-6 relative rms in single
    precision, and a collected power that agrees to a small fraction of a
    dB. Loss is positive dB, so a positive dB difference means the device
    collected LESS power.

(c) THE SPEED. One 1024 px `split_step` under the three backends.

The script writes `cupy_backend_check.json` next to itself, and it prints a
table. Run it on the machine with the CUDA device:

    cd D:\\repos\\optical_link_budget
    C:\\Users\\alexf\\olb-gpu-venv\\Scripts\\python.exe validation\\gpu_fft\\cupy_backend_check.py

Documentation uses ASD-STE100 Simplified Technical English.
"""

import argparse
import json
import os
import platform
import sys
import time

import numpy as np

from olb.geometry import CircularOrbit
from olb.scenario import Channel, Site, SpaceScenario
from olb.terminal import SMF, Terminal, Transmitter
from olb.waveoptics.field import Begin, Power, to_host
from olb.waveoptics.grid import GridSpec
from olb.waveoptics.priority import boost_process_priority
from olb.waveoptics.propagators import (clear_forvard_cache, set_fft_backend)
from olb.waveoptics.turbulence.sampling import turbulent_grid
from olb.waveoptics.turbulence.screens import ScreenFactory
from olb.waveoptics.turbulence.splitstep import (split_step,
                                                 super_gaussian_boundary)

try:
    import cupy as cp
except ImportError:                       # the host parts still run
    cp = None

LAM_M = 1550e-9
SEED = 20260907


def hero_scenario():
    """The 30 deg hero downlink: 0.7 m SMF ground terminal, 500 km.

    It is a copy of `gpu_microbench.hero_scenario`, so the two scripts
    measure the same link.
    """
    ground = Terminal(aperture_m=0.7, obscuration_ratio=0.3,
                      wavelength_m=LAM_M, pointing_jitter_rad=1e-6,
                      detector=SMF(sensitivity_dbm=-40.0, focal_length_m=2.0,
                                   mode_field_radius_m=5.2e-6))
    space = Terminal(aperture_m=0.1, obscuration_ratio=0.0,
                     wavelength_m=LAM_M, pointing_jitter_rad=1e-6,
                     transmitter=Transmitter(waist_m=0.03, power_dbm=30.0,
                                             divergence_rad=None))
    scn = SpaceScenario(ground=ground, space=space,
                        channel=Channel(site=Site(), altitude_m=500e3),
                        direction="downlink")
    return scn, CircularOrbit(altitude_m=500e3, elevation_deg=30.0)


def pinned(scn, geo, preset, n, L0_m=25.0):
    """Size with the preset, then pin the pixel count only.

    It is a copy of `gpu_microbench.pinned`.
    """
    g, plan, _ = turbulent_grid(scn, geo, preset=preset, L0_m=L0_m)
    return GridSpec(size_m=g.size_m, n=int(n), scaled=g.scaled), plan


def rel_rms(a, b):
    """Give the rms difference of two arrays over the rms of the reference."""
    a = np.asarray(a)
    b = np.asarray(b)
    return float(np.sqrt(np.mean(np.abs(a - b) ** 2)
                         / np.mean(np.abs(b) ** 2)))


def sync():
    """Wait for the device. It is a no-op with no cupy."""
    if cp is not None:
        cp.cuda.Device().synchronize()


def best_of(fn, reps):
    """Give the best wall time of `reps` calls, in s."""
    best = np.inf
    for _ in range(reps):
        t0 = time.perf_counter()
        fn()
        sync()
        best = min(best, time.perf_counter() - t0)
    return best


def make_screens(backend, n, pixel_m, r0_m, L0_m=25.0):
    """Make one screen for each r0 under one FFT backend.

    The factory reads the backend in __init__, so the function sets the
    backend first. The seeds and the draw order do not change with the
    backend.

    Returns:
        A list of screens. They are device arrays under the "cupy" backend.
    """
    previous = set_fft_backend(backend)
    try:
        factory = ScreenFactory(n, pixel_m, L0_m=L0_m, dtype=np.float32)
        return [factory.make(float(r0), np.random.default_rng(SEED + 100 + j))
                for j, r0 in enumerate(r0_m)]
    finally:
        set_fft_backend(previous)


def run_split_step(backend, grid, plan, screens, mask):
    """Run one split step under one FFT backend. Give the host field.

    The screens are HOST arrays, so the two backends see exactly the same
    atmosphere and the difference is the transforms only.
    """
    previous = set_fft_backend(backend)
    try:
        F0 = Begin(grid.size_m, LAM_M, grid.n, dtype=np.complex64)
        F_rx = split_step(F0, plan.z_m, list(screens), plan.z_total_m,
                          boundary=mask)
        power = Power(F_rx)
        return to_host(F_rx).field.copy(), float(power)
    finally:
        set_fft_backend(previous)


def time_split_step(backend, grid, plan, screens, mask, reps):
    """Give the best wall time of one split step, in s."""
    previous = set_fft_backend(backend)
    try:
        F0 = Begin(grid.size_m, LAM_M, grid.n, dtype=np.complex64)
        scr = list(screens)

        def call():
            split_step(F0, plan.z_m, list(scr), plan.z_total_m, boundary=mask)
        call()                              # warm the caches and the plans
        sync()
        return best_of(call, reps)
    finally:
        set_fft_backend(previous)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--n", type=int, nargs="+", default=[512, 1024])
    ap.add_argument("--reps", type=int, default=3)
    ap.add_argument("--out", default=os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "cupy_backend_check.json"))
    args = ap.parse_args()
    boost_process_priority()

    out = {"environment": {
        "python": sys.version.split()[0], "numpy": np.__version__,
        "cupy": None if cp is None else cp.__version__,
        "cores": os.cpu_count(),
        "gpu": (None if cp is None else
                cp.cuda.runtime.getDeviceProperties(0)["name"].decode()),
        "platform": platform.platform()},
        "seed": SEED, "precision": "single", "L0_m": 25.0,
        "preset": "standard"}
    print("environment:", json.dumps(out["environment"]))
    if cp is None:
        print("cupy is absent. The script needs the CUDA device. Stop.")
        return 1

    scn, geo = hero_scenario()
    out["cases"] = []
    for n in args.n:
        grid, plan = pinned(scn, geo, "standard", n)
        n_screens = int(plan.z_m.size)
        row = {"n": n, "n_screens": n_screens,
               "grid_size_m": float(grid.size_m)}
        print(f"\n{n} px, {n_screens} screens, side {grid.size_m:.2f} m")

        # ---- (a) the atmosphere ----
        clear_forvard_cache()
        s_host = make_screens("numpy", n, grid.pixel_m, plan.r0_m)
        s_dev = make_screens("cupy", n, grid.pixel_m, plan.r0_m)
        errs = [rel_rms(cp.asnumpy(d), h) for h, d in zip(s_host, s_dev)]
        row["screen_rel_rms"] = errs
        row["screen_rel_rms_max"] = float(max(errs))
        row["screen_rms_rad"] = float(np.sqrt(np.mean(s_host[0] ** 2)))
        print(f"  (a) screen rel rms, {n_screens} screens: "
              f"max {max(errs):.2e}, min {min(errs):.2e}")

        # ---- (b) the propagation ----
        mask = super_gaussian_boundary(n)
        f_host, p_host = run_split_step("numpy", grid, plan, s_host, mask)
        f_dev, p_dev = run_split_step("cupy", grid, plan, s_host, mask)
        row["field_rel_rms"] = rel_rms(f_dev, f_host)
        row["power_numpy"] = p_host
        row["power_cupy"] = p_dev
        # Loss is positive dB: a positive number means the device collected
        # less power than the host.
        row["power_diff_db"] = float(-10.0 * np.log10(p_dev / p_host))
        print(f"  (b) receive field rel rms {row['field_rel_rms']:.2e}, "
              f"power difference {row['power_diff_db']:+.2e} dB")

        # ---- (c) the speed ----
        row["split_step_s"] = {}
        for backend in ("numpy", "scipy", "cupy"):
            scr = s_dev if backend == "cupy" else s_host
            dt = time_split_step(backend, grid, plan, scr, mask, args.reps)
            row["split_step_s"][backend] = dt
            print(f"  (c) split_step {backend:5s} {dt * 1e3:9.2f} ms")
        base = row["split_step_s"]["numpy"]
        row["speedup_vs_numpy"] = {
            k: base / v for k, v in row["split_step_s"].items()}
        out["cases"].append(row)

    with open(args.out, "w") as fh:
        json.dump(out, fh, indent=2)

    print("\n| n px | screens | screen rel rms | field rel rms | dP dB | "
          "numpy ms | scipy ms | cupy ms | cupy speed-up |")
    print("|---|---|---|---|---|---|---|---|---|")
    for r in out["cases"]:
        t = r["split_step_s"]
        print(f"| {r['n']} | {r['n_screens']} | "
              f"{r['screen_rel_rms_max']:.1e} | {r['field_rel_rms']:.1e} | "
              f"{r['power_diff_db']:+.1e} | {t['numpy'] * 1e3:.1f} | "
              f"{t['scipy'] * 1e3:.1f} | {t['cupy'] * 1e3:.1f} | "
              f"{r['speedup_vs_numpy']['cupy']:.1f}x |")
    print(f"\nwrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
