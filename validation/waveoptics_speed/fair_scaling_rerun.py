r"""A FAIR re-measurement of threads against processes for fidelity-2 trials.

WHY. The first study (scaling_study.py) has four fairness defects. This script
keeps its cases, its setup and its trial function (it IMPORTS them, it does not
copy them) and it corrects the four:

  1. TRIAL COUNT. The old study ran 32 trials for each point, so W = 16 gave 2
     trials for each worker. That measures the tail, not the steady state. This
     script runs N_TRIALS = 200 for each point.
  2. BLAS THREADS. This script sets OPENBLAS_NUM_THREADS, MKL_NUM_THREADS and
     OMP_NUM_THREADS to 1 BEFORE numpy is imported. A spawned child inherits the
     environment, so BOTH routes run BLAS with one thread. This is DELIBERATE:
     it stops each route from over-subscribing the cores, and it makes the
     thread route and the process route comparable.
  3. CHUNKING. The process route maps with chunksize = ceil(n_trials / W), so
     each worker gets ONE contiguous block (for example 10 processes with 20
     trials each). That is the work shape the owner wants to test.
  4. SPAWN. The process route reports THREE numbers: spawn_s (the pool creation
     and the initializer, forced with one warm task for each worker), steady_s
     (the 200 trials), and wall_s = spawn_s + steady_s. The thread route reports
     wall_s only, because a thread pool starts in milliseconds. The comparison
     uses WALL time, and it also shows the steady-state ratio.

A FIFTH measurement gives the machine ceiling: a pure fft2 process sweep with NO
Python trial logic. Each of W processes runs K in-place np.fft.fft2 calls on a
complex128 (N, N) array, N = the grid size of the case. It shows where the
memory bandwidth of the machine plateaus with no GIL and no Python work.

REPORT ONLY. This script changes NO production code under olb/.

Run from the repository root:
    python validation/waveoptics_speed/fair_scaling_rerun.py
    python validation/waveoptics_speed/fair_scaling_rerun.py --smoke --out x.json
"""

import os

# The BLAS pin MUST come before numpy is imported. See defect 2 above.
BLAS_ENV = {"OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1",
            "OMP_NUM_THREADS": "1"}
os.environ.update(BLAS_ENV)

import json
import math
import platform
import sys
import time
import warnings
from concurrent.futures import ProcessPoolExecutor

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__)))))

import numpy as np

from validation.waveoptics_speed.scaling_study import (SEED, _attach_builder,
                                                       _proc_init, _proc_trial,
                                                       _proc_warm, build_setup,
                                                       thread_scaling,
                                                       warm_aotools)

HERE = os.path.dirname(os.path.abspath(__file__))
GENERATOR = "olb"                       # the production default generator
CASE_NAMES = ["space downlink 30deg standard", "terrestrial 2km standard"]
N_TRIALS = 200
BASELINE_TRIALS = 40                    # the cheap one-worker rate baseline
WORKERS = [4, 8, 12, 16, 24]            # W = 1 comes from the serial baseline
FFT_WORKERS = [1, 4, 8, 12, 16, 24]
K_FFT = 20


# ---------------------------------------------------------------------------
# 5. the pure-fft2 bandwidth ceiling. The worker is module level, for spawn.
# ---------------------------------------------------------------------------

def fft_worker(args):
    """Run k in-place fft2 calls on a complex128 (n, n) array. Time inside."""
    n, k = args
    rng = np.random.default_rng(7)
    a = (rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n)))
    np.fft.fft2(a)                      # one warm-up transform, not timed
    t0 = time.perf_counter()
    for _ in range(k):
        a = np.fft.fft2(a)
    return time.perf_counter() - t0


def fft_ceiling(n, workers, k):
    """Sweep the process count on raw fft2. Report the aggregate fft2 rate."""
    rows = []
    base = None
    for w in workers:
        with ProcessPoolExecutor(max_workers=w) as pool:
            list(pool.map(_proc_warm, range(w)))     # force every worker up
            times = list(pool.map(fft_worker, [(n, k)] * w, chunksize=1))
        # The workers time themselves, so the rate holds no dispatch cost. The
        # slowest worker sets the aggregate rate.
        wall = max(times)
        rate = w * k / wall
        base = rate if base is None else base
        rows.append({"workers": w, "wall_s": wall,
                     "worker_s_mean": float(np.mean(times)),
                     "fft2_per_s": rate,
                     "speedup": rate / base, "efficiency": (rate / base) / w})
    return rows


# ---------------------------------------------------------------------------
# 2. processes, with one contiguous block for each worker
# ---------------------------------------------------------------------------

def process_scaling_fair(name, generator, workers, n_trials):
    """Sweep the process count. Report spawn, steady state and wall time."""
    rows = []
    for w in workers:
        t0 = time.perf_counter()
        pool = ProcessPoolExecutor(max_workers=w, initializer=_proc_init,
                                   initargs=(name, generator))
        list(pool.map(_proc_warm, range(w)))      # force spawn + initializer
        spawn_s = time.perf_counter() - t0

        chunk = max(1, math.ceil(n_trials / w))
        t0 = time.perf_counter()
        got = list(pool.map(_proc_trial, range(n_trials), chunksize=chunk))
        steady_s = time.perf_counter() - t0
        pool.shutdown(wait=True)

        rows.append({"workers": w, "chunksize": chunk, "spawn_s": spawn_s,
                     "steady_s": steady_s, "wall_s": spawn_s + steady_s,
                     "steady_trials_per_s": n_trials / steady_s,
                     "wall_trials_per_s": n_trials / (spawn_s + steady_s),
                     "collected_first": float(got[0])})
    return rows


# ---------------------------------------------------------------------------
# the plateau read: the smallest W inside 5% of the best rate
# ---------------------------------------------------------------------------

def _plateau(rows, key):
    best = max(r[key] for r in rows)
    for r in rows:
        if r[key] >= 0.95 * best:
            return r["workers"], best
    return rows[-1]["workers"], best


def _add_efficiency(rows, key, base_rate):
    for r in rows:
        r["speedup_over_serial"] = r[key] / base_rate
        r["efficiency"] = r["speedup_over_serial"] / r["workers"]
    return rows


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main():
    smoke = "--smoke" in sys.argv
    out_path = os.path.join(HERE, "fair_scaling_rerun_results.json")
    if "--out" in sys.argv:
        out_path = sys.argv[sys.argv.index("--out") + 1]

    n_trials = 8 if smoke else N_TRIALS
    base_trials = 4 if smoke else BASELINE_TRIALS
    workers = [1, 4] if smoke else WORKERS
    fft_workers = [1, 4] if smoke else FFT_WORKERS
    k_fft = 3 if smoke else K_FFT

    t_start = time.time()
    warm_aotools()
    cpu = os.cpu_count() or 1
    env = {
        "numpy": np.__version__,
        "scipy": __import__("scipy").__version__,
        "cpu": platform.processor(),
        "cores": cpu,
        "platform": platform.platform(),
        "blas_env": {k: os.environ[k] for k in BLAS_ENV},
        "generator": GENERATOR,
        "n_trials": n_trials,
        "baseline_trials": base_trials,
        "workers": workers,
        "fft_workers": fft_workers,
        "k_fft": k_fft,
        "smoke": smoke,
        "seed": SEED,
        "baseline_note": ("W=1 is a SEPARATE serial baseline of "
                          f"{base_trials} trials, extrapolated to a rate, so "
                          "the 200-trial sweeps do not pay a one-worker run"),
    }
    print("environment:")
    for k, v in env.items():
        print(f"  {k:18s} {v}")
    print("")

    results = {}
    for name in CASE_NAMES:
        print(f"case {name}: baseline ...", flush=True)
        # The one-worker baseline: the production runner with one thread.
        base = thread_scaling(name, GENERATOR, [1],
                              n_trials=base_trials)[0]["trials_per_s"]

        print(f"case {name}: threads ...", flush=True)
        thr = thread_scaling(name, GENERATOR, workers, n_trials=n_trials)
        _add_efficiency(thr, "trials_per_s", base)

        print(f"case {name}: processes ...", flush=True)
        prc = process_scaling_fair(name, GENERATOR, workers, n_trials)
        _add_efficiency(prc, "steady_trials_per_s", base)

        setup = build_setup(name, GENERATOR)
        n_grid = int(setup["grid"].n)
        print(f"case {name}: fft2 ceiling (N={n_grid}) ...", flush=True)
        fft = fft_ceiling(n_grid, fft_workers, k_fft)

        results[name] = {"grid_n": n_grid, "baseline_trials_per_s": base,
                         "threads": thr, "processes": prc, "fft": fft}

    # ---- the tables ----
    for name, r in results.items():
        print("")
        print("=" * 74)
        print(f"{name}   (N={r['grid_n']}, {n_trials} trials, generator "
              f"{GENERATOR})")
        print(f"  serial baseline: {r['baseline_trials_per_s']:.3f} trials/s "
              f"({base_trials} trials)")
        print("=" * 74)
        print(f"(a) THREADS   {'W':>4}{'wall_s':>10}{'trials/s':>11}"
              f"{'speedup':>9}{'eff':>7}")
        for w in r["threads"]:
            print(f"              {w['workers']:>4}{w['wall_s']:>10.2f}"
                  f"{w['trials_per_s']:>11.2f}"
                  f"{w['speedup_over_serial']:>9.2f}{w['efficiency']:>7.2f}")
        print(f"(b) PROCESSES {'W':>4}{'chunk':>7}{'spawn_s':>9}"
              f"{'steady_s':>10}{'wall_s':>9}{'steady/s':>10}{'eff':>7}")
        for w in r["processes"]:
            print(f"              {w['workers']:>4}{w['chunksize']:>7}"
                  f"{w['spawn_s']:>9.2f}{w['steady_s']:>10.2f}"
                  f"{w['wall_s']:>9.2f}{w['steady_trials_per_s']:>10.2f}"
                  f"{w['efficiency']:>7.2f}")
        print(f"(c) FFT2      {'W':>4}{'wall_s':>10}{'fft2/s':>11}"
              f"{'speedup':>9}{'eff':>7}")
        for w in r["fft"]:
            print(f"              {w['workers']:>4}{w['wall_s']:>10.2f}"
                  f"{w['fft2_per_s']:>11.2f}{w['speedup']:>9.2f}"
                  f"{w['efficiency']:>7.2f}")

    # ---- the recommendation ----
    print("")
    print("=" * 74)
    print("RECOMMENDATION  (wall time decides; the steady ratio is shown too)")
    print("=" * 74)
    reco = {}
    for name, r in results.items():
        t_best = min(r["threads"], key=lambda x: x["wall_s"])
        p_best = min(r["processes"], key=lambda x: x["wall_s"])
        route = "processes" if p_best["wall_s"] < t_best["wall_s"] else "threads"
        best = p_best if route == "processes" else t_best
        t_steady = max(x["trials_per_s"] for x in r["threads"])
        p_steady = max(x["steady_trials_per_s"] for x in r["processes"])
        p_plateau, _ = _plateau(r["processes"], "steady_trials_per_s")
        f_plateau, _ = _plateau(r["fft"], "fft2_per_s")
        entry = {
            "best_route": route, "best_workers": best["workers"],
            "best_wall_s": best["wall_s"],
            "thread_best_wall_s": t_best["wall_s"],
            "thread_best_workers": t_best["workers"],
            "process_best_wall_s": p_best["wall_s"],
            "process_best_workers": p_best["workers"],
            "process_over_thread_wall": t_best["wall_s"] / p_best["wall_s"],
            "process_over_thread_steady": p_steady / t_steady,
            "process_plateau_w": p_plateau, "fft_plateau_w": f_plateau,
            "plateau_matches": p_plateau == f_plateau,
        }
        reco[name] = entry
        print(f"  {name}")
        print(f"     best wall: {route} at W={best['workers']} "
              f"({best['wall_s']:.2f} s for {n_trials} trials)")
        print(f"     threads best W={t_best['workers']} "
              f"{t_best['wall_s']:.2f} s | processes best "
              f"W={p_best['workers']} {p_best['wall_s']:.2f} s")
        print(f"     processes/threads: wall "
              f"{entry['process_over_thread_wall']:.2f}x, steady "
              f"{entry['process_over_thread_steady']:.2f}x "
              f"(>1 means processes win)")
        print(f"     plateau W: processes {p_plateau}, fft2 ceiling "
              f"{f_plateau} -> "
              f"{'MATCH' if entry['plateau_matches'] else 'DIFFERENT'}")

    out = {"environment": env, "cases": results, "recommendation": reco}
    with open(out_path, "w") as f:
        json.dump(out, f, indent=1)
    print("")
    print(f"wrote {out_path}")
    print(f"(elapsed {time.time() - t_start:.1f} s)")


if __name__ == "__main__":
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        main()
