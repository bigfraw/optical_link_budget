"""Verify the 2026-09-06 memory cut of a fidelity-2 trial, and measure it.

Three claims, each one checked here:

1. The LAZY SCREENS. `split_step` reads its screens one at a time, so a
   generator keeps ONE screen in memory. A generator and a list of the same
   screens give the SAME field, bit for bit.
2. The FORVARD CACHE. `Forvard` keeps its sign pattern and its wrapped
   transfer function for each distinct hop. A warm cache gives the SAME field
   as a cold one, bit for bit, and a trial with the cache OFF
   (`FORVARD_CACHE_BYTES = 0`, the pre-2026-09-06 behaviour) gives the same
   numbers as a trial with the cache ON.
3. The POOL SIZER. `worker_memory_bytes` sits ABOVE the measured peak of one
   serial trial, so a pool that `auto_workers` sizes does not run out of
   memory, and `auto_workers` obeys both limits.

The script then measures the wall time and the traced peak memory of one
serial trial with the cache OFF and ON, at 1024 px, in both precisions. It
writes `memory_cut_check.json` next to itself.

Run it from the repository root, on the laptop or in a container. DO NOT run
it on a busy machine: the timing is the point.

    python -m validation.memory_cut.memory_cut_check
    python -m validation.memory_cut.memory_cut_check --n 512 --trials 2
"""
import argparse
import json
import os
import time
import tracemalloc
import warnings

import numpy as np

from olb.geometry import CircularOrbit
from olb.scenario import Channel, Site, SpaceScenario
from olb.terminal import SMF, Terminal, Transmitter
from olb.waveoptics import propagators as P
from olb.waveoptics.field import Begin
from olb.waveoptics.grid import GridSpec
from olb.waveoptics.resources import auto_workers, worker_memory_bytes
from olb.waveoptics.turbulence.run import propagate_turbulent_scenario
from olb.waveoptics.turbulence.sampling import turbulent_grid
from olb.waveoptics.turbulence.screens import ScreenFactory
from olb.waveoptics.turbulence.splitstep import split_step, super_gaussian_boundary

HERE = os.path.dirname(os.path.abspath(__file__))


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


def check_lazy_screens():
    """Claim 1: a generator and a list give the same field."""
    n, dx, lam = 256, 0.01, 1550e-9
    fac = ScreenFactory(n, dx, L0_m=25.0, dtype=np.float32)
    screens = [fac.make(0.05, np.random.default_rng(k)) for k in range(4)]
    z = np.array([100.0, 300.0, 600.0, 900.0])
    mask = super_gaussian_boundary(n)
    F0 = Begin(n * dx, lam, n, dtype=np.complex64)
    a = split_step(F0, z, list(screens), 1000.0, boundary=mask)
    b = split_step(F0, z, (s for s in screens), 1000.0, boundary=mask)
    assert np.array_equal(a.field, b.field), "generator differs from list"
    for bad in ([screens[0]], screens + [screens[0]]):
        try:
            split_step(F0, z, iter(bad), 1000.0, boundary=mask)
            raise AssertionError("a wrong screen count must raise")
        except ValueError:
            pass
    print("1. lazy screens: a generator equals a list, bit for bit")


def check_forvard_cache():
    """Claim 2: warm equals cold, and cache off equals cache on."""
    n, lam = 256, 1550e-9
    F0 = Begin(2.0, lam, n, dtype=np.complex64)
    F0.field = F0.field * super_gaussian_boundary(n, 0.3)
    P.clear_forvard_cache()
    cold = P.Forvard(F0, 1234.5)
    warm = P.Forvard(F0, 1234.5)
    assert P.forvard_cache_bytes() > 0
    assert np.array_equal(cold.field, warm.field)
    saved = P.FORVARD_CACHE_BYTES
    try:
        P.FORVARD_CACHE_BYTES = 0
        P.clear_forvard_cache()
        off = P.Forvard(F0, 1234.5)
        assert P.forvard_cache_bytes() == 0
    finally:
        P.FORVARD_CACHE_BYTES = saved
    assert np.array_equal(cold.field, off.field)
    print("2. Forvard cache: warm equals cold equals off, bit for bit")


def check_sizer(n, precision, peak_bytes, n_hops):
    """Claim 3: the estimate sits above the measured peak."""
    est = worker_memory_bytes(n, precision, block_size=1, patch_pixels=0,
                              n_hops=n_hops)
    assert est > peak_bytes, (est, peak_bytes)
    assert auto_workers(2 ** 30, free_bytes=2 ** 30, cores=8) == (1, "memory")
    assert auto_workers(2 ** 20, free_bytes=2 ** 40, cores=8) == (7, "cpu")
    return est


def measure(scn, geo, grid, plan, precision, trials, cache_on):
    """One serial run: the mean wall time and the traced peak of a trial."""
    saved = P.FORVARD_CACHE_BYTES
    P.FORVARD_CACHE_BYTES = saved if cache_on else 0
    P.clear_forvard_cache()
    try:
        kw = dict(seed=7, preset="standard", threader=None,
                  precision=precision, L0_m=25.0, grid=grid, plan=plan)
        propagate_turbulent_scenario(scn, geo, n_trials=1, **kw)   # warm up
        tracemalloc.start()
        t0 = time.perf_counter()
        res = propagate_turbulent_scenario(scn, geo, n_trials=trials, **kw)
        dt = (time.perf_counter() - t0) / trials
        _, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
    finally:
        P.FORVARD_CACHE_BYTES = saved
        P.clear_forvard_cache()
    return dt, peak, res


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--n", type=int, default=1024, help="the pixel count")
    ap.add_argument("--trials", type=int, default=3)
    args = ap.parse_args()
    warnings.simplefilter("ignore")

    check_lazy_screens()
    check_forvard_cache()

    scn, geo = hero_scenario()
    out = {"grid_n": args.n, "trials": args.trials,
           "cores": os.cpu_count(), "cases": []}
    print(f"3. one serial trial at {args.n} px ({os.cpu_count()} cores)")
    for precision, preset in (("single", "standard"), ("double", "reference")):
        grid, plan = pinned(scn, geo, preset, args.n)
        row = {"precision": precision, "preset": preset,
               "n_screens": int(plan.z_m.size)}
        results = {}
        for cache_on in (False, True):
            dt, peak, res = measure(scn, geo, grid, plan, precision,
                                    args.trials, cache_on)
            tag = "cache_on" if cache_on else "cache_off"
            row[tag] = {"s_per_trial": dt, "peak_mib": peak / 2 ** 20}
            results[tag] = res
            print(f"   {precision:6s} {plan.z_m.size:2d} screens, cache "
                  f"{'on ' if cache_on else 'off'}: {dt:6.2f} s/trial, "
                  f"peak {peak / 2 ** 20:6.0f} MiB")
        same = all(a.collected_power == b.collected_power
                   and a.smf_eta == b.smf_eta
                   for a, b in zip(results["cache_off"].trials,
                                   results["cache_on"].trials))
        assert same, "the cache changed a number"
        est = check_sizer(args.n, precision, row["cache_on"]["peak_mib"]
                          * 2 ** 20, int(plan.z_m.size) + 1)
        row["worker_estimate_mib"] = est / 2 ** 20
        print(f"          cache on equals cache off, bit for bit; "
              f"the sizer estimates {est / 2 ** 20:.0f} MiB per worker")
        out["cases"].append(row)

    path = os.path.join(HERE, "memory_cut_check.json")
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=2)
    print(f"wrote {path}")
    print("memory_cut_check passed")


if __name__ == "__main__":
    main()
