"""The end-to-end check of the cupy backend (backlog 2-N8, milestone 2).

Milestone 1 certified the kernels (`cupy_backend_check.py`). This script
checks the ONE knob a caller has, `fft_backend`, through the runner and
through `Campaign`. It answers two questions:

(a) THE RUNNER. Does `propagate_turbulent_scenario` give the same trials
    under "cupy" as under "numpy"? The two runs share the seed, the grid and
    the screen plan, and the white noise stays a host numpy PCG64 draw, so
    the two see the SAME atmosphere and they must agree at the rounding
    level of float32 (about 1e-6 to 1e-5 relative). The script reports the
    trial-for-trial difference of `collected_power` and `smf_eta`, and the
    wall time of one trial for each backend.

(b) THE CAMPAIGN. How fast is ONE GPU stream against the 12-worker CPU pool
    of record? The script stores the SAME campaign two times, under "cupy"
    (one process, whatever `workers` says) and under "numpy" with a pool,
    and it reports the trials for each second of each. The two stores have
    DIFFERENT keys, because the backend enters the fingerprint, so the
    script compares the stored scalars by trial index.

Run it on the machine with the CUDA device, with the pool idle:

    cd D:\\repos\\optical_link_budget
    C:\\Users\\alexf\\olb-gpu-venv\\Scripts\\python.exe validation\\gpu_fft\\cupy_campaign_check.py

The campaigns go under `validation/gpu_fft/_campaigns/` (git-ignored). Give
`--fresh` to delete them first.

Documentation uses ASD-STE100 Simplified Technical English.
"""

import argparse
import json
import os
import platform
import shutil
import sys
import time

import numpy as np

from olb.geometry import CircularOrbit
from olb.scenario import Channel, Site, SpaceScenario
from olb.terminal import SMF, Terminal, Transmitter
from olb.waveoptics.grid import GridSpec
from olb.waveoptics.priority import boost_process_priority
from olb.waveoptics.turbulence.campaign import Campaign
from olb.waveoptics.turbulence.run import propagate_turbulent_scenario
from olb.waveoptics.turbulence.sampling import turbulent_grid

try:
    import cupy as cp
except ImportError:                       # the host parts still run
    cp = None

LAM_M = 1550e-9
SEED = 7
L0_M = 25.0
HERE = os.path.dirname(os.path.abspath(__file__))


def hero_scenario():
    """The 30 deg hero downlink: 0.7 m SMF ground terminal, 500 km.

    It is a copy of `gpu_microbench.hero_scenario`, so every script of this
    directory measures the same link.
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


def pinned(scn, geo, preset, n, L0_m=L0_M):
    """Size with the preset, then pin the pixel count only."""
    g, plan, _ = turbulent_grid(scn, geo, preset=preset, L0_m=L0_m)
    return GridSpec(size_m=g.size_m, n=int(n), scaled=g.scaled), plan


def max_rel(a, b):
    """The largest relative difference of two arrays, as a float."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    keep = np.isfinite(a) & np.isfinite(b) & (b != 0.0)
    if not keep.any():
        return float("nan")
    return float(np.max(np.abs(a[keep] / b[keep] - 1.0)))


# ---------------------------------------------------------------- part (a)

def runner_case(scn, geo, n, n_trials):
    """Run the same trials under numpy and under cupy, and compare them."""
    grid, plan = pinned(scn, geo, "standard", n)
    row = {"n": n, "n_screens": int(plan.z_m.size), "n_trials": n_trials}
    got = {}
    for backend in ("numpy", "cupy"):
        t0 = time.perf_counter()
        res = propagate_turbulent_scenario(
            scn, geo, n_trials=n_trials, seed=SEED, preset="standard",
            grid=grid, plan=plan, L0_m=L0_M, precision="single",
            fft_backend=backend)
        dt = time.perf_counter() - t0
        got[backend] = res
        row[f"s_per_trial_{backend}"] = dt / n_trials
        print(f"  {backend:5s}: {dt / n_trials * 1e3:8.1f} ms per trial")
    for name in ("collected_power", "smf_eta"):
        a = [getattr(t, name) for t in got["cupy"].trials]
        b = [getattr(t, name) for t in got["numpy"].trials]
        row[f"max_rel_{name}"] = max_rel(a, b)
        row[f"{name}_numpy"] = [float(v) for v in b]
        row[f"{name}_cupy"] = [float(v) for v in a]
        print(f"  max relative difference, {name:16s} "
              f"{row[f'max_rel_{name}']:.2e}")
    row["speedup"] = row["s_per_trial_numpy"] / row["s_per_trial_cupy"]
    return row


# ---------------------------------------------------------------- part (b)

def campaign_case(scn, geo, n, n_trials, block_size, workers, root, fresh):
    """Store the same campaign under cupy and under numpy, and time both."""
    grid, plan = pinned(scn, geo, "standard", n)
    n_blocks = -(-n_trials // block_size)
    if n_blocks < workers:
        print(f"  CAUTION: {n_blocks} blocks for {workers} workers. The pool "
              "leaves workers idle, so the CPU number reads low.")
    row = {"n": n, "n_screens": int(plan.z_m.size), "n_trials": n_trials,
           "block_size": block_size, "n_blocks": n_blocks,
           "workers_numpy": workers}
    stores = {}
    for backend in ("cupy", "numpy"):
        path = os.path.join(root, f"n{n}_{backend}")
        if fresh and os.path.isdir(path):
            shutil.rmtree(path)
        camp = Campaign(scn, geo, path, seed=SEED, preset="standard",
                        block_size=block_size, grid=grid, plan=plan,
                        L0_m=L0_M, precision="single", fft_backend=backend)
        t0 = time.perf_counter()
        # The cupy campaign ignores `workers` and says so. The numpy campaign
        # opens the pool of record.
        camp.run(n_trials, workers=(None if backend == "cupy" else workers),
                 progress=True)
        dt = time.perf_counter() - t0
        row[f"s_total_{backend}"] = dt
        row[f"trials_per_s_{backend}"] = n_trials / dt
        row[f"key_{backend}"] = camp.fingerprint
        stores[backend] = camp.load(n_trials, fields=False)
        print(f"  {backend:5s}: {n_trials / dt:7.2f} trials/s "
              f"({dt:.1f} s for {n_trials})")
    assert row["key_cupy"] != row["key_numpy"], "the two keys must differ"
    for name in ("collected_power", "smf_eta"):
        a = [getattr(t, name) for t in stores["cupy"].trials]
        b = [getattr(t, name) for t in stores["numpy"].trials]
        row[f"max_rel_{name}"] = max_rel(a, b)
        print(f"  max relative difference, {name:16s} "
              f"{row[f'max_rel_{name}']:.2e}")
    row["speedup"] = row["trials_per_s_cupy"] / row["trials_per_s_numpy"]
    return row


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--trials", type=int, default=10,
                    help="the trials of part (a). The one-time host setup "
                         "(the vacuum baseline) divides over them, so a "
                         "small count reads high.")
    ap.add_argument("--campaign-trials", type=int, default=120,
                    help="the trials of the 1024 px campaign pair")
    ap.add_argument("--block", type=int, default=10,
                    help="the block size of the 1024 px campaign pair. THE "
                         "BLOCK COUNT MUST REACH THE POOL SIZE: a pool with "
                         "more workers than blocks leaves workers idle, and "
                         "the CPU number then reads low.")
    ap.add_argument("--workers", type=int, default=12,
                    help="the numpy pool size (the pool of record)")
    ap.add_argument("--big-trials", type=int, default=48,
                    help="the trials of the 2048 px campaign pair")
    ap.add_argument("--big-block", type=int, default=4,
                    help="the block size of the 2048 px campaign pair")
    ap.add_argument("--big-budget-s", type=float, default=180.0,
                    help="run the 2048 px pair only when the 1024 px pair "
                         "took less than this")
    ap.add_argument("--root", default=os.path.join(HERE, "_campaigns"))
    ap.add_argument("--fresh", action="store_true",
                    help="delete the stored campaigns first")
    ap.add_argument("--out", default=os.path.join(
        HERE, "cupy_campaign_check.json"))
    args = ap.parse_args()

    boost_process_priority()
    out = {"environment": {
        "python": sys.version.split()[0],
        "numpy": np.__version__,
        "cupy": None if cp is None else cp.__version__,
        "cores": os.cpu_count(),
        "gpu": (None if cp is None else
                cp.cuda.runtime.getDeviceProperties(0)["name"].decode()),
        "platform": platform.platform()},
        "seed": SEED, "precision": "single", "L0_m": L0_M,
        "preset": "standard"}
    print("environment:", json.dumps(out["environment"]))
    if cp is None:
        print("cupy is absent. The script needs the CUDA device. Stop.")
        return 1
    os.makedirs(args.root, exist_ok=True)

    scn, geo = hero_scenario()

    print("\n(a) the runner, 1024 px, numpy against cupy")
    out["runner"] = runner_case(scn, geo, 1024, args.trials)

    print("\n(b) the campaign pair, 1024 px")
    t0 = time.perf_counter()
    out["campaigns"] = [campaign_case(
        scn, geo, 1024, args.campaign_trials, args.block, args.workers,
        args.root, args.fresh)]
    pair_s = time.perf_counter() - t0

    if pair_s < args.big_budget_s:
        print(f"\n(c) the campaign pair, 2048 px "
              f"(the 1024 px pair took {pair_s:.0f} s)")
        out["campaigns"].append(campaign_case(
            scn, geo, 2048, args.big_trials, args.big_block, args.workers,
            args.root, args.fresh))
    else:
        print(f"\n(c) skipped: the 1024 px pair took {pair_s:.0f} s, more "
              f"than the {args.big_budget_s:.0f} s budget.")

    with open(args.out, "w") as fh:
        json.dump(out, fh, indent=2)

    r = out["runner"]
    print("\n(a) the runner, one trial")
    print("| n px | screens | trials | numpy ms | cupy ms | speed-up | "
          "max rel dP | max rel eta |")
    print("|---|---|---|---|---|---|---|---|")
    print(f"| {r['n']} | {r['n_screens']} | {r['n_trials']} | "
          f"{r['s_per_trial_numpy'] * 1e3:.0f} | "
          f"{r['s_per_trial_cupy'] * 1e3:.0f} | {r['speedup']:.1f}x | "
          f"{r['max_rel_collected_power']:.1e} | "
          f"{r['max_rel_smf_eta']:.1e} |")

    print("\n(b) and (c) one GPU stream against the CPU pool")
    print("| n px | trials | cupy trials/s | numpy pool trials/s | workers | "
          "speed-up | max rel dP | max rel eta |")
    print("|---|---|---|---|---|---|---|---|")
    for c in out["campaigns"]:
        print(f"| {c['n']} | {c['n_trials']} | "
              f"{c['trials_per_s_cupy']:.2f} | "
              f"{c['trials_per_s_numpy']:.2f} | {c['workers_numpy']} | "
              f"{c['speedup']:.1f}x | "
              f"{c['max_rel_collected_power']:.1e} | "
              f"{c['max_rel_smf_eta']:.1e} |")
    print(f"\nwrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
