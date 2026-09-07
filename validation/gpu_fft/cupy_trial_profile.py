"""Where the real cupy trial spends its time (backlog 2-N8).

THE QUESTION. The synthetic microbench (`gpu_microbench.py`) predicted 66 ms
for one 1024 px trial with the threaded host draw, and 349 ms at 2048 px. The
real `fft_backend="cupy"` trial measures 333 ms at 1024 px
(`cupy_campaign_check.json`). So about 260 ms of a 1024 px trial is host work
that the device never sees, and the microbench does not hold it. This script
finds that work.

HOW. It runs the REAL `propagate_turbulent_scenario`, and it wraps the stages
of the trial loop with timers. The wrap is a monkeypatch of the module
globals that `run_one` reads, so the measured code IS the shipped code. The
timers are in this script and NEVER in the package.

Each timer reads `cupy.cuda.Device().synchronize()` before it stops, so a
device stage carries its own device time and not the time of the next
synchronising call.

THE STAGES.

- `draw_wait`      the main thread waiting for the host noise of this trial.
                   A wait near the whole draw time says the pipeline does not
                   overlap.
- `noise_upload`   the host-to-device copy of the noise of one screen.
- `screen_filter`  the filter multiply and the inverse transform of one screen.
- `screen_sub`     the subharmonic levels of one screen.
- `forvard`        one Forvard hop (the sign, the two transforms, the transfer
                   function).
- `screen_apply`   the phase screen into the field (the upload, if any, and
                   the multiply).
- `mask`           the absorbing boundary multiply.
- `to_host`        the ONE download of the receive field.
- `patch_store`    the masked patch write.
- `clip`           the receive-aperture clip.
- `power`          the collected power.
- `detector_eta`   the fibre coupling (the fibre mode build plus the overlap).

It also counts the Forvard factor cache misses, so a cache that thrashes is
visible as a count and not as a guess.

The script writes `cupy_trial_profile.json`, and it prints the table. Pass
`--tag before` or `--tag after` to keep the two runs of one fix in the SAME
file: the record is a dictionary of tags.

Run it on a machine with a CUDA device and the GPU venv:

    ssh desktop "cd D:\\repos\\optical_link_budget; \
        C:\\Users\\alexf\\olb-gpu-venv\\Scripts\\python.exe \
        validation\\gpu_fft\\cupy_trial_profile.py --tag before"
"""

import argparse
import cProfile
import io
import json
import os
import pstats
import sys
import time
from collections import defaultdict

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__)))))

from gpu_microbench import hero_scenario, pinned          # noqa: E402

from olb.waveoptics import propagators as _prop           # noqa: E402
from olb.waveoptics.turbulence import run as _run         # noqa: E402
from olb.waveoptics.turbulence import screens as _screens  # noqa: E402
from olb.waveoptics.turbulence import splitstep as _ss    # noqa: E402
from olb.waveoptics.turbulence.run import (               # noqa: E402
    propagate_turbulent_scenario)

try:
    import cupy as cp
except ImportError:                                       # pragma: no cover
    cp = None

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "cupy_trial_profile.json")


# ------------------------------------------------------------------ timers

class Timers:
    """Collect the total time and the call count of each named stage."""

    def __init__(self):
        self.total = defaultdict(float)
        self.count = defaultdict(int)
        self.on = False

    def reset(self):
        self.total.clear()
        self.count.clear()

    def sync(self):
        """Wait for the device, so a device stage owns its own time."""
        if cp is not None:
            cp.cuda.Device().synchronize()

    def wrap(self, name, fn, sync=True):
        """Give a function that runs `fn` and adds its time to `name`."""
        def timed(*a, **kw):
            if not self.on:
                return fn(*a, **kw)
            t0 = time.perf_counter()
            out = fn(*a, **kw)
            if sync:
                self.sync()
            self.total[name] += time.perf_counter() - t0
            self.count[name] += 1
            return out
        timed.__name__ = getattr(fn, "__name__", name)
        return timed


T = Timers()


class TimedFuture:
    """A future whose `result()` time counts as `draw_wait`."""

    def __init__(self, inner):
        self._inner = inner

    def result(self, *a, **kw):
        t0 = time.perf_counter()
        out = self._inner.result(*a, **kw)
        T.total["draw_wait"] += time.perf_counter() - t0
        T.count["draw_wait"] += 1
        return out


def install(n_px):
    """Monkeypatch the stages of the trial loop, and give the undo function.

    The patch replaces module globals only. Nothing in the package changes,
    and `restore` puts every name back.
    """
    saved = []

    def patch(mod, name, new):
        saved.append((mod, name, getattr(mod, name)))
        setattr(mod, name, new)

    # ---- the tail of run_one, all host code today ----
    patch(_run, "to_host", T.wrap("to_host", _run.to_host))
    patch(_run, "_clip", T.wrap("clip", _run._clip))
    patch(_run, "Power", T.wrap("power", _run.Power))
    patch(_run, "_detector_eta", T.wrap("detector_eta", _run._detector_eta))

    # ---- the propagation, hop by hop ----
    patch(_ss, "Forvard", T.wrap("forvard", _ss.Forvard))
    patch(_ss, "Screen", T.wrap("screen_apply", _ss.Screen))
    patch(_ss, "_apply_mask", T.wrap("mask", _ss._apply_mask))

    # ---- the screen build, in three parts ----
    factory_cls = _screens.ScreenFactory
    base_from = factory_cls._base_pair_from
    sub_from = factory_cls._subharmonic_from

    def timed_base(self, g):
        if not T.on:
            return base_from(self, g)
        t0 = time.perf_counter()
        if self._device:
            g = self._xp.asarray(g)
            T.sync()
            T.total["noise_upload"] += time.perf_counter() - t0
            T.count["noise_upload"] += 1
            # The dtype of the uploaded noise is the finding of check (iv).
            T.total["_noise_bytes"] = float(g.nbytes)
        t1 = time.perf_counter()
        cn = g * self._filt
        full = self._ift_series(cn)
        out = self._xp.real(full), self._xp.imag(full)
        T.sync()
        T.total["screen_filter"] += time.perf_counter() - t1
        T.count["screen_filter"] += 1
        return out

    patch(factory_cls, "_base_pair_from", timed_base)
    patch(factory_cls, "_subharmonic_from",
          T.wrap("screen_sub", sub_from))

    # ---- the patch store, timed through the fields array write ----
    # run_one writes `fields[k - start_index] = ...`. The write is one line of
    # run_one, so it cannot be patched from here. It is measured separately
    # by the no-patch / patch pair of runs (see `measure`).

    # ---- the draw pipeline ----
    import concurrent.futures as _cf
    real_pool = _cf.ThreadPoolExecutor

    class TimedPool(real_pool):
        def submit(self, fn, *a, **kw):
            return TimedFuture(real_pool.submit(self, fn, *a, **kw))

    saved.append((_cf, "ThreadPoolExecutor", real_pool))
    _cf.ThreadPoolExecutor = TimedPool

    # ---- the Forvard factor cache misses ----
    real_factors = _prop._forvard_factors
    state = {"calls": 0, "misses": 0}

    def counted(N, size, lam, z, cdtype):
        state["calls"] += 1
        cdt = np.dtype(cdtype)
        key = (int(N), float(size), float(lam), float(z), cdt.str,
               _prop._cache_place())
        if key not in _prop._forvard_cache:
            state["misses"] += 1
        return real_factors(N, size, lam, z, cdtype)

    saved.append((_prop, "_forvard_factors", real_factors))
    _prop._forvard_factors = counted

    def restore():
        for mod, name, old in reversed(saved):
            setattr(mod, name, old)

    return restore, state


# ------------------------------------------------------------- the measure

def measure(n_px, trials, backend, patch_radius_m):
    """Run the hero downlink and give the per-stage breakdown of one trial.

    Args:
        n_px:           the pixel count of one side.
        trials:         the number of timed trials, after one warm-up.
        backend:        "cupy" or "numpy".
        patch_radius_m: the stored patch radius, in m, or None.

    Returns:
        A dictionary of the record.
    """
    scn, geo = hero_scenario()
    grid, plan = pinned(scn, geo, "standard", n_px)
    common = dict(preset="standard", grid=grid, plan=plan, L0_m=25.0,
                  precision="single", fft_backend=backend, seed=7,
                  patch_radius_m=patch_radius_m)

    restore, state = install(n_px)
    try:
        T.on = False
        # The warm-up fills the Forvard cache and the cuFFT plan cache, so the
        # timed runs measure the steady state.
        propagate_turbulent_scenario(scn, geo, n_trials=1, **common)

        def timed_run(n):
            """Give the wall time, the stage totals and the miss count of a
            run of n trials. Each call HOLDS the one-time setup of the runner
            (the vacuum baseline, the transmit mode), so the difference of two
            runs is the marginal cost of one trial."""
            state["calls"] = state["misses"] = 0
            T.reset()
            T.on = True
            t0 = time.perf_counter()
            out = propagate_turbulent_scenario(scn, geo, n_trials=n, **common)
            T.sync()
            wall = time.perf_counter() - t0
            T.on = False
            return (wall, dict(T.total), dict(T.count),
                    state["calls"], state["misses"], out)

        # THE SETUP IS SUBTRACTED. The runner makes a VACUUM baseline before
        # the trials, and that baseline runs the same split step. So a single
        # run over-counts every stage. The difference of a run of n trials and
        # a run of one trial, over n - 1, is the cost of ONE trial and it
        # holds no setup.
        w1, tot1, cnt1, calls1, miss1, _ = timed_run(1)
        wn, totn, cntn, callsn, missn, res = timed_run(trials + 1)
    finally:
        restore()

    d = float(trials)
    noise1 = tot1.pop("_noise_bytes", 0.0)
    noise_bytes = totn.pop("_noise_bytes", noise1)
    keys = set(tot1) | set(totn)
    stages = {}
    for k in keys:
        ms = 1e3 * (totn.get(k, 0.0) - tot1.get(k, 0.0)) / d
        stages[k] = {"ms_per_trial": ms,
                     "calls_per_trial": (cntn.get(k, 0) - cnt1.get(k, 0)) / d}
    stages = dict(sorted(stages.items(),
                         key=lambda kv: -kv[1]["ms_per_trial"]))
    named = sum(v["ms_per_trial"] for v in stages.values())
    wall_ms = 1e3 * (wn - w1) / d
    row = {
        "n_px": n_px,
        "backend": backend,
        "trials": trials,
        "n_screens": int(plan.z_m.size),
        "patch_radius_m": patch_radius_m,
        # A run of ONE trial is the SETUP plus one trial. The setup is the
        # vacuum baseline, the transmit mode and the start field, and a
        # Campaign pays it at EVERY block. So a small block hides a large
        # fixed cost, and this number says how large.
        "ms_setup_plus_one_trial": 1e3 * w1,
        "ms_setup": 1e3 * w1 - wall_ms,
        "ms_per_trial_wall": wall_ms,
        "ms_per_trial_named": named,
        "ms_per_trial_unattributed": wall_ms - named,
        "forvard_factor_calls_per_trial": (callsn - calls1) / d,
        "forvard_factor_misses_per_trial": (missn - miss1) / d,
        "forvard_cache_bytes": _prop.forvard_cache_bytes(),
        "uploaded_noise_bytes": noise_bytes,
        "stages": stages,
        "collected_power_0": float(res.trials[0].collected_power),
        "smf_eta_0": float(res.trials[0].smf_eta),
    }
    return row


def show(row):
    """Print one record as a table."""
    print(f"\n--- {row['backend']} {row['n_px']} px, "
          f"{row['n_screens']} screens, {row['trials']} trials ---")
    print(f"  wall            {row['ms_per_trial_wall']:9.1f} ms/trial")
    print(f"  setup, one time {row['ms_setup']:9.1f} ms per RUN")
    print(f"  named stages    {row['ms_per_trial_named']:9.1f} ms/trial")
    print(f"  unattributed    {row['ms_per_trial_unattributed']:9.1f} ms/trial")
    print(f"  factor misses   {row['forvard_factor_misses_per_trial']:9.1f} "
          f"of {row['forvard_factor_calls_per_trial']:.0f} calls/trial")
    print(f"  cache holds     {row['forvard_cache_bytes'] / 2**20:9.1f} MiB")
    print(f"  noise upload    {row['uploaded_noise_bytes'] / 2**20:9.2f} MiB "
          "per screen")
    print(f"  {'stage':<16}{'ms/trial':>10}{'calls':>9}")
    for name, v in row["stages"].items():
        print(f"  {name:<16}{v['ms_per_trial']:10.2f}"
              f"{v['calls_per_trial']:9.1f}")


def profile_calls(n_px, trials):
    """Run cProfile over the trials and give the top 15 by cumulative time."""
    scn, geo = hero_scenario()
    grid, plan = pinned(scn, geo, "standard", n_px)
    kw = dict(preset="standard", grid=grid, plan=plan, L0_m=25.0,
              precision="single", fft_backend="cupy", seed=7)
    propagate_turbulent_scenario(scn, geo, n_trials=1, **kw)   # warm up
    pr = cProfile.Profile()
    pr.enable()
    propagate_turbulent_scenario(scn, geo, n_trials=trials, **kw)
    pr.disable()
    buf = io.StringIO()
    pstats.Stats(pr, stream=buf).sort_stats("cumulative").print_stats(15)
    text = buf.getvalue()
    print(text)
    return text


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tag", default="before",
                    help="the key of this record in the JSON file")
    ap.add_argument("--n", type=int, nargs="+", default=[1024, 2048])
    ap.add_argument("--trials", type=int, default=5)
    ap.add_argument("--no-cprofile", action="store_true")
    args = ap.parse_args()

    if cp is None:
        print("cupy is absent. This script needs a CUDA device.")
        return 1
    dev = cp.cuda.Device()
    print(f"cupy {cp.__version__}, device {dev.id}, "
          f"free {dev.mem_info[0] / 2**30:.2f} GiB of "
          f"{dev.mem_info[1] / 2**30:.2f} GiB")

    rows = []
    for n_px in args.n:
        # A patch of 0.35 m is the receive aperture radius of the hero
        # terminal, so the store is the store a real campaign makes.
        for radius in (None, 0.35):
            rows.append(measure(n_px, args.trials, "cupy", radius))
            show(rows[-1])

    prof = None
    if not args.no_cprofile:
        prof = profile_calls(args.n[0], args.trials)

    record = {}
    if os.path.exists(OUT):
        with open(OUT) as fh:
            record = json.load(fh)
    # THE ROWS MERGE. A tag collects the rows of every run that carries it,
    # keyed by the grid and the patch, so a second call with another --n adds
    # a row and it does not drop the first one. A repeat of the same case
    # replaces its own row.
    tag = record.setdefault(args.tag, {})
    kept = {(r["n_px"], r["patch_radius_m"]): r for r in tag.get("rows", [])}
    kept.update({(r["n_px"], r["patch_radius_m"]): r for r in rows})
    tag.update({
        "when": time.strftime("%Y-%m-%d %H:%M:%S"),
        "cupy_version": cp.__version__,
        "forvard_cache_bytes_limit": _prop.FORVARD_CACHE_BYTES,
        "rows": [kept[k] for k in sorted(kept, key=lambda k: (k[0], k[1] or 0))],
    })
    if prof is not None:
        tag["cprofile_top15"] = prof
    with open(OUT, "w") as fh:
        json.dump(record, fh, indent=1)
    print(f"\nwrote {OUT} under the tag {args.tag!r}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
