"""Measure the two 2026-09-06 OPT-INS against the defaults of record: the
lean screen generator (`screen_generator="olb-lean"`) and the scipy FFT
backend (`fft_backend="scipy"`).

Neither one is bit-identical to the default, so neither one is a default.
This script gives the numbers that a decision needs:

1. THE DRAWS. The lean body reads the SAME random stream as the default
   body. For one seed, the two screens agree to the rounding of the output
   type: float64 to about 1e-16, float32 to about 1e-7 relative rms.
2. THE STATISTICS. The structure-function r0 fit of 40 lean screens equals
   the fit of 40 default screens inside the standard error.
3. THE SCREEN COST. The wall time and the traced peak memory of one screen,
   both bodies, both types, at 1024 and 2048 px.
4. THE TRANSFORM. The raw 2-d transform of a 1024 px complex64 array,
   numpy against scipy.
5. THE TRIAL. One serial single-precision trial of the 30 deg hero downlink
   at 1024 px, the four combinations of generator and backend, interleaved
   and best of three, and the agreement of the collected power and the SMF
   eta between each combination and the default.

It writes `screen_generator_lean.json` next to itself.

Run it from the repository root, on the laptop or in a container. DO NOT run
it on a busy machine: the timing is the point.

    python -m validation.memory_cut.screen_generator_lean
"""
import argparse
import json
import os
import time
import tracemalloc
import warnings

import numpy as np
from scipy import fft as sfft

from olb.waveoptics.turbulence.run import propagate_turbulent_scenario
from olb.waveoptics.turbulence.screens import ScreenFactory

from .memory_cut_check import hero_scenario, pinned

HERE = os.path.dirname(os.path.abspath(__file__))


def structure_fn_r0(scr, dx, seps=(2, 4, 8, 16)):
    """Fit r0 from D(r) = 6.88 (r / r0)^(5/3) at a few small separations.

    Andrews and Phillips, DOI 10.1117/3.626196, Ch. 3, Eq. (75), printed
    p. 82 (the Kolmogorov phase structure function).
    """
    D = np.array([np.mean((scr[:, k:] - scr[:, :-k]).astype(np.float64) ** 2)
                  for k in seps])
    r = np.array(seps) * dx
    return float(np.mean(r / (D / 6.88) ** 0.6))


def draws(out):
    print("1. the draws, same seed, r0 = 0.1 m, L0 = 25 m")
    rows = []
    for n, dx in ((512, 6.86e-3), (1024, 3.43e-3)):
        for dt in (np.float64, np.float32):
            a = ScreenFactory(n, dx, L0_m=25.0, dtype=dt).make(
                0.1, np.random.default_rng(5))
            b = ScreenFactory(n, dx, L0_m=25.0, dtype=dt, lean=True).make(
                0.1, np.random.default_rng(5))
            rel = float(np.sqrt(np.mean((a - b) ** 2) / np.mean(a ** 2)))
            rows.append({"n": n, "dtype": dt.__name__, "rel_rms": rel})
            print(f"   {n:5d} px {dt.__name__:8s}: rel rms {rel:.2e}")
            assert rel < (1e-12 if dt is np.float64 else 1e-5)
    out["draws"] = rows


def statistics(out):
    print("2. the structure-function r0 of 40 screens, 1024 px float32")
    n, dx = 1024, 3.43e-3
    rows = {}
    for lean in (False, True):
        f = ScreenFactory(n, dx, L0_m=25.0, dtype=np.float32, lean=lean)
        r0 = [structure_fn_r0(f.make(0.1, np.random.default_rng(100 + k)), dx)
              for k in range(40)]
        rows["lean" if lean else "default"] = {
            "r0_fit_m": float(np.mean(r0)),
            "sem_m": float(np.std(r0) / np.sqrt(len(r0)))}
        print(f"   {'lean   ' if lean else 'default'}: r0 fit "
              f"{np.mean(r0):.4f} +- {np.std(r0) / np.sqrt(len(r0)):.4f} m")
    gap = abs(rows["lean"]["r0_fit_m"] - rows["default"]["r0_fit_m"])
    assert gap < 3 * rows["default"]["sem_m"], gap
    out["statistics"] = rows


def screen_cost(out, count=30):
    print(f"3. one screen, {count} draws, r0 = 0.1 m, L0 = 25 m")
    rows = []
    for n, dx in ((1024, 3.43e-3), (2048, 1.715e-3)):
        for dt in (np.float32, np.float64):
            for lean in (False, True):
                f = ScreenFactory(n, dx, L0_m=25.0, dtype=dt, lean=lean)
                f.make(0.1, np.random.default_rng(0))
                tracemalloc.start()
                t0 = time.perf_counter()
                for k in range(count):
                    f.make(0.1, np.random.default_rng(k))
                dt_s = (time.perf_counter() - t0) / count
                _, peak = tracemalloc.get_traced_memory()
                tracemalloc.stop()
                rows.append({"n": n, "dtype": dt.__name__, "lean": lean,
                             "ms_per_screen": dt_s * 1e3,
                             "peak_mib": peak / 2 ** 20})
                print(f"   {n:5d} px {dt.__name__:8s} "
                      f"{'lean   ' if lean else 'default'}: "
                      f"{dt_s * 1e3:7.1f} ms, peak {peak / 2 ** 20:6.1f} MiB")
    out["screen_cost"] = rows


def transform(out):
    print("4. the raw 1024 px complex64 fft2")
    a = (np.random.default_rng(0).standard_normal((1024, 1024))
         + 0j).astype(np.complex64)
    rows = {}
    for name, f in (("numpy", np.fft.fft2),
                    ("scipy", lambda x: sfft.fft2(x, workers=1)),
                    ("scipy_overwrite",
                     lambda x: sfft.fft2(x.copy(), overwrite_x=True,
                                         workers=1))):
        f(a)
        t0 = time.perf_counter()
        for _ in range(10):
            f(a)
        rows[name] = (time.perf_counter() - t0) / 10 * 1e3
        print(f"   {name:16s} {rows[name]:6.1f} ms")
    out["transform_ms"] = rows


def trial(out, n, reps):
    print(f"5. one serial single-precision trial at {n} px, best of {reps}")
    scn, geo = hero_scenario()
    grid, plan = pinned(scn, geo, "standard", n)
    combos = [("numpy", "olb"), ("scipy", "olb"), ("numpy", "olb-lean"),
              ("scipy", "olb-lean")]

    def run(backend, gen, trials=3):
        return propagate_turbulent_scenario(
            scn, geo, n_trials=trials, seed=7, preset="standard",
            threader=None, precision="single", L0_m=25.0, grid=grid,
            plan=plan, screen_generator=gen, fft_backend=backend)

    ref = run("numpy", "olb", 1)
    times = {c: [] for c in combos}
    for _ in range(reps):
        for c in combos:
            t0 = time.perf_counter()
            run(*c)
            times[c].append((time.perf_counter() - t0) / 3)
    rows = []
    for c in combos:
        got = run(*c, trials=1)
        d_p = abs(got.trials[0].collected_power
                  / ref.trials[0].collected_power - 1.0)
        d_e = abs(got.trials[0].smf_eta / ref.trials[0].smf_eta - 1.0)
        rows.append({"fft_backend": c[0], "screen_generator": c[1],
                     "s_per_trial": min(times[c]),
                     "rel_diff_power": d_p, "rel_diff_eta": d_e})
        print(f"   {c[0]:5s} + {c[1]:8s}: {min(times[c]):5.2f} s/trial, "
              f"power {d_p:.1e}, eta {d_e:.1e} from the default")
        assert d_p < 1e-5 and d_e < 1e-5
    out["trial"] = {"n": n, "n_screens": int(plan.z_m.size), "rows": rows}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--n", type=int, default=1024)
    ap.add_argument("--reps", type=int, default=3)
    args = ap.parse_args()
    warnings.simplefilter("ignore")
    out = {"cores": os.cpu_count()}
    draws(out)
    statistics(out)
    screen_cost(out)
    transform(out)
    trial(out, args.n, args.reps)
    path = os.path.join(HERE, "screen_generator_lean.json")
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=2)
    print(f"wrote {path}")
    print("screen_generator_lean passed")


if __name__ == "__main__":
    main()
