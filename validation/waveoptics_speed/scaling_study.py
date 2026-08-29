r"""How fidelity-2 turbulent trials scale across workers: threads, processes,
and a batched split step.

THE POINT. P0 (profile_baseline.py) shows screen generation is 80 to 84% of a
serial trial. P1 (screen_generator_check.py) adds a fast, cached olb screen
generator (ScreenFactory) that is 7 to 14 times faster per screen, so it SHIFTS
the cost balance toward the split step. This script measures how one whole trial
scales when many trials run at the same time, on THREE routes, for BOTH screen
generators:

  1. THREADS. The production path: propagate_turbulent_scenario with a Threader.
     The pocketfft transform releases the GIL, so threads help. This script
     sweeps the worker count and reports trials per second and the parallel
     efficiency, then it says where the curve saturates.
  2. PROCESSES. A ProcessPoolExecutor prototype with a module-level worker and a
     once-per-process initializer that rebuilds the read-only setup. It measures
     the spawn-plus-initialise overhead apart from the steady-state rate, and it
     checks that the scenario dataclasses pickle.
  3. BATCHED FFT. A prototype that stacks B trials into a (B, N, N) array and
     runs the free-space hops with one transform over the last two axes (numpy
     vectorises the leading axis; scipy adds a workers pool). It reports trials
     per second against B and W, and an honest memory budget.

It also cross-checks the seed contract: each parallel route gives the SAME trial
result as the serial loop for the same seed.

REPORT ONLY. This script changes NO production code. It reads the production
split-step layer and it builds its process and batched prototypes on top.

RESULTS AND RECOMMENDATION. See scaling_study_results.json and the printed
RECOMMENDATION block. The one-line summary is regenerated on each run and pasted
back here by the owner; until then read the JSON. The measured table has the
form (per case, per generator): the best mode, its worker or batch count, and
the speed-up over one worker.

RECOMMENDATION (2026-08-29, 32 cores, numpy 1.26.4, scipy 1.11.4,
aotools 1.0.7, Windows 11; see scaling_study_run.log). Best mode, worker/batch
count, speed-up over one worker:

  case | generator                         mode   count   x1
  terrestrial 2km standard | aotools    processes     8   6.1
  terrestrial 2km standard | olb        processes    16   4.5
  space downlink 30deg rapid | aotools    threads    16   3.5
  space downlink 30deg rapid | olb      processes    16   6.5
  space downlink 30deg standard| aotools processes   16   9.9
  space downlink 30deg standard| olb    processes     8   5.0

READ OF THE NUMBERS.
- PROCESSES beat THREADS on 5 of the 6 (case, generator) points, by 1.35x to
  1.79x on the best rate. The GIL is the reason: the split step is more than the
  FFT (screen exp, mask, clip, and the aotools screen draw all run in Python),
  so many threads contend on the interpreter while separate processes do not.
  The one tie is the rapid aotools case (0.96x), where the tiny grid (N=256)
  makes the process steady-state win too small to pay for the spawn.
- THREADS saturate at W = 8 to 16 and the efficiency has fallen to about 0.35
  by W = 16 and about 0.15 by W = 32. The plateau near half the cores is the
  memory-bandwidth wall the P0 raw fft2 numbers predict (fft2 is bandwidth bound
  at N >= 1024), COMPOUNDED by the GIL on the non-FFT Python work.
- PROCESSES pay a real, measured spawn+init cost on Windows spawn: 2.4 to 8.8 s
  for the pool, rising with the worker count. It is a ONE-TIME cost. It only
  pays off past a few hundred trials; below that, threads win outright because
  they spawn in milliseconds.
- The olb generator lifts the ONE-WORKER rate 2 to 8x (it removes ~70% of the
  serial trial, per P1), so the absolute best rates come with olb. It also makes
  each trial more FFT-bound, which is why olb saturates slightly earlier.
- BATCHED FFT is NOT competitive here. On the N=512 standard case its best is
  2.82 trials/s (olb, scipy workers=32), against 5.9 trials/s for processes and
  4.3 for threads. scipy workers help the batch (numpy 1.75 -> scipy_w32 2.82),
  but a single batched split step still runs the screen and mask work
  single-process, and the (B, N, N) FFT is not faster per trial than B separate
  FFTs across processes. Batching wins memory locality, not wall time, on this
  machine. Keep it as a GPU/low-core fallback, not the default.

A BETTER THREADER DEFAULT? The Threader default is one thread per core (32).
The data says that OVER-subscribes: on every case the thread rate peaks at
W = 8 to 16 and then falls or flattens (efficiency 0.15 to 0.35 at W = 32). A
default of about half the cores (min(16, cores) or so) would give nearly the
peak rate at a quarter of the threads. This is a RECOMMENDATION only; this
script changes nothing.

Sources:
- Schmidt, Numerical Simulation of Optical Wave Propagation with Examples in
  MATLAB, DOI 10.1117/3.866274, Ch. 9. The split-step method and the absorbing
  boundary that the batched prototype reproduces.
- Andrews and Phillips, Laser Beam Propagation through Random Media, 2nd ed.,
  DOI 10.1117/3.626196. The scintillation the trials measure.

Run from the repository root:
    python -m validation.waveoptics_speed.scaling_study
"""

import json
import os
import pickle
import platform
import time
import warnings
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import scipy.fft as _sfft

from olb.geometry import CircularOrbit, HorizontalPath
from olb.scenario import (Channel, SpaceScenario, TerrestrialChannel,
                          TerrestrialScenario)
from olb.terminal import SMF, Terminal, Transmitter
from olb.turbulence.profiles import DEFAULT_HS, default_cn2_profile
from olb.waveoptics.field import Begin, Field, Power
from olb.waveoptics.run import _clip, _launch_aperture, _normalised_gauss
from olb.waveoptics.sources import GaussBeam
from olb.beam import virtual_waist
from olb.waveoptics.threader import Threader
from olb.waveoptics.turbulence.run import (_resolve_seed, _screen_builder,
                                           _screen_seed,
                                           propagate_turbulent_scenario)
from olb.waveoptics.turbulence.sampling import PRESETS, turbulent_grid
from olb.waveoptics.turbulence.screens import phase_screen
from olb.waveoptics.turbulence.splitstep import (_substeps,
                                                 super_gaussian_boundary)

LAM = 1550e-9
HERE = os.path.dirname(__file__)
SEED = 2026
TRIALS = 32                       # about 32 trials per point (a good citizen)
GENERATORS = ("aotools", "olb")


def _worker_list():
    """The worker sweep: 1, 2, 4, ... up to the core count."""
    cores = os.cpu_count() or 1
    ws, w = [], 1
    while w < cores:
        ws.append(w)
        w *= 2
    ws.append(cores)
    return sorted(set(ws))


# ---------------------------------------------------------------------------
# the cases
# ---------------------------------------------------------------------------
# A case is a light key. The builders rebuild the scenario and the geometry from
# the key, so a fresh worker process needs no pickled scenario.

def _terrestrial_case():
    scn = TerrestrialScenario(
        near=Terminal(aperture_m=0.10, wavelength_m=LAM,
                      transmitter=Transmitter(waist_m=0.02)),
        far=Terminal(aperture_m=0.20, wavelength_m=LAM, detector=SMF()),
        channel=TerrestrialChannel(path_length_m=2000.0, cn2=5e-15))
    return scn, HorizontalPath(2000.0), None, None


def _space_case(direction):
    ground = Terminal(aperture_m=0.50, wavelength_m=LAM,
                      transmitter=Transmitter(waist_m=0.05), detector=SMF())
    scn = SpaceScenario(ground=ground,
                        space=Terminal(aperture_m=0.30, wavelength_m=LAM),
                        direction=direction, channel=Channel(altitude_m=600e3))
    geom = CircularOrbit(altitude_m=600e3, elevation_deg=[30.0])
    hs = DEFAULT_HS
    cn2 = default_cn2_profile(scn.channel.site, hs)
    return scn, geom, hs, cn2


# name -> (preset, builder). The builders match the P0 profile cases.
CASES = {
    "terrestrial 2km standard": ("standard", _terrestrial_case),
    "space downlink 30deg rapid": ("rapid", lambda: _space_case("downlink")),
    "space downlink 30deg standard": ("standard",
                                      lambda: _space_case("downlink")),
}
SCALING_CASES = list(CASES)          # every case feeds the thread/process sweep
BATCH_CASE = "space downlink 30deg standard"   # the batched prototype case


def _build_case(name):
    """Rebuild (scenario, geometry, hs, cn2, preset_name) from a case name."""
    preset_name, builder = CASES[name]
    scn, geom, hs, cn2 = builder()
    return scn, geom, hs, cn2, preset_name


def warm_aotools():
    """Trigger the one-time scipy deprecation aotools raises, quietly.

    aotools 1.0.7 reads the deprecated scipy.ndimage.interpolation namespace on
    its first call, once per process. Trigger it here so it does not land in a
    timed block.
    """
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            phase_screen(0.1, 32, 0.01, seed=0)
    except ImportError:
        pass


# ---------------------------------------------------------------------------
# the read-only setup and one trial, shared by the serial, thread reference,
# process and batched routes
# ---------------------------------------------------------------------------
# This reproduces the collected_power face of run_one in
# olb.waveoptics.turbulence.run. The cross-check proves it is bit-identical to
# the production propagate_turbulent_scenario, so the process and batched routes
# inherit that guarantee.

def build_setup(name, generator):
    """Build the read-only setup a trial needs, for one case and generator.

    Returns a dict with the grid, the plan, the boundary mask, the reference
    power, the receive aperture, and a screen builder. It is the once-per-run
    (and once-per-process) work.
    """
    scn, geom, hs, cn2, preset_name = _build_case(name)
    p = PRESETS[preset_name]
    is_space = hasattr(scn, "direction")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        grid, plan, report = turbulent_grid(scn, geom, preset=p, hs=hs,
                                            cn2_profile=cn2)
    lam = scn.tx_terminal.wavelength_m
    rx = scn.ground if is_space else scn.rx_terminal
    mask = super_gaussian_boundary(grid.n, p.boundary_width_frac)
    max_step_m = grid.n * grid.pixel_m ** 2 / lam

    if is_space:
        # THE VACUUM BASELINE: the same plane wave and the same hops with flat
        # screens. The turbulent collected_power then has a vacuum limit of 1.0.
        from olb.waveoptics.turbulence.splitstep import split_step
        F_plane = Begin(grid.size_m, lam, grid.n)
        flat = np.zeros((grid.n, grid.n))
        F_vac = split_step(F_plane, plan.z_m, [flat] * int(plan.z_m.size),
                           plan.z_total_m, boundary=mask)
        p_reference = Power(_clip(F_vac, rx.aperture_m, rx.obscuration_ratio))
        F_in = None
    else:
        tx = scn.tx_terminal
        t = tx.transmitter
        w_v, offset = virtual_waist(t.waist_m, t.divergence_rad, lam)
        F0 = _normalised_gauss(GaussBeam(Begin(grid.size_m, lam, grid.n), w_v))
        if offset > 0:
            from olb.waveoptics.propagators import GForvard
            F0 = GForvard(F0, offset)
        F_in = _clip(F0, *_launch_aperture(tx))
        p_reference = Power(F_in)

    build_screen = _screen_builder(generator, grid, np.inf, True)
    return {
        "name": name, "generator": generator, "is_space": is_space,
        "grid": grid, "plan": plan, "mask": mask, "max_step_m": max_step_m,
        "p_reference": float(p_reference), "lam": lam, "F_in": F_in,
        "aperture_m": rx.aperture_m, "obscuration": rx.obscuration_ratio,
        "n_screens": int(plan.z_m.size),
        "report": report,
    }


def trial_collected_power(setup, entropy, k):
    """Run trial k and return its collected_power. It mirrors run_one."""
    from olb.waveoptics.turbulence.splitstep import split_step
    plan = setup["plan"]
    grid = setup["grid"]
    build = setup["build_screen"]
    stack = [build(_screen_seed(entropy, k, j), plan.r0_m[j])
             for j in range(setup["n_screens"])]
    F_start = (Begin(grid.size_m, setup["lam"], grid.n) if setup["is_space"]
               else setup["F_in"])
    F_rx = split_step(F_start, plan.z_m, stack, plan.z_total_m,
                      boundary=setup["mask"])
    collected = _clip(F_rx, setup["aperture_m"], setup["obscuration"])
    return float(Power(collected) / setup["p_reference"])


# build_setup does not store the screen builder (a closure does not pickle and a
# process rebuilds its own). Attach it after the fact for the in-process routes.
def _attach_builder(setup):
    setup["build_screen"] = _screen_builder(setup["generator"], setup["grid"],
                                            np.inf, True)
    return setup


def serial_collected(setup, entropy, n_trials):
    """The serial reference: collected_power of trials 0..n-1."""
    return np.array([trial_collected_power(setup, entropy, k)
                     for k in range(n_trials)])


# ---------------------------------------------------------------------------
# 1. THREADS: the production path
# ---------------------------------------------------------------------------

def thread_scaling(name, generator, workers, n_trials=TRIALS):
    """Sweep the Threader worker count on the production runner."""
    scn, geom, hs, cn2, preset_name = _build_case(name)
    rows = []
    base_rate = None
    for w in workers:
        threader = Threader(max_workers=w)
        t0 = time.perf_counter()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            propagate_turbulent_scenario(
                scn, geom, n_trials=n_trials, seed=SEED, preset=preset_name,
                hs=hs, cn2_profile=cn2, threader=threader,
                screen_generator=generator)
        wall = time.perf_counter() - t0
        rate = n_trials / wall
        if base_rate is None:
            base_rate = rate
        rows.append({"workers": w, "wall_s": wall, "trials_per_s": rate,
                     "speedup": rate / base_rate,
                     "efficiency": (rate / base_rate) / w})
    return rows


# ---------------------------------------------------------------------------
# 2. PROCESSES: a module-level worker with a once-per-process initializer
# ---------------------------------------------------------------------------

_WORKER = {}                       # per-process global: the read-only setup


def _proc_init(name, generator):
    """Initialise a worker process: warm aotools, build the setup once."""
    warm_aotools()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        _WORKER["setup"] = _attach_builder(build_setup(name, generator))
    _WORKER["entropy"] = _resolve_seed(SEED)


def _proc_trial(k):
    """Run one trial by index in a worker process."""
    return trial_collected_power(_WORKER["setup"], _WORKER["entropy"], k)


def _proc_warm(_):
    """A tiny task that forces a worker to spawn and initialise."""
    return os.getpid()


def process_scaling(name, generator, workers, n_trials=TRIALS):
    """Sweep the ProcessPool worker count. Split spawn from steady state."""
    rows = []
    base_rate = None
    for w in workers:
        # ---- spawn + initialise: force every worker up with tiny tasks ----
        t0 = time.perf_counter()
        pool = ProcessPoolExecutor(max_workers=w, initializer=_proc_init,
                                   initargs=(name, generator))
        # One warm task per worker (and a few extra) forces all to start.
        list(pool.map(_proc_warm, range(max(w, 4))))
        spawn_s = time.perf_counter() - t0

        # ---- steady state: the real trials ----
        t0 = time.perf_counter()
        got = list(pool.map(_proc_trial, range(n_trials)))
        steady_s = time.perf_counter() - t0
        pool.shutdown(wait=True)

        rate = n_trials / steady_s
        if base_rate is None:
            base_rate = rate
        rows.append({"workers": w, "spawn_init_s": spawn_s,
                     "steady_s": steady_s, "trials_per_s": rate,
                     "speedup": rate / base_rate,
                     "efficiency": (rate / base_rate) / w,
                     "collected_first": float(got[0])})
    return rows


def scenario_pickles():
    """Check that the scenario dataclasses of every case pickle, round-trip."""
    out = {}
    for name in CASES:
        scn, geom, _, _, _ = _build_case(name)
        entry = {"scenario": True, "geometry": True, "error": None}
        try:
            s2 = pickle.loads(pickle.dumps(scn))
            g2 = pickle.loads(pickle.dumps(geom))
            entry["scenario"] = s2.tx_terminal.aperture_m == scn.tx_terminal.aperture_m
            entry["geometry"] = np.allclose(np.atleast_1d(g2.slant_range_m),
                                            np.atleast_1d(geom.slant_range_m))
        except Exception as exc:              # report, do not raise
            entry["scenario"] = entry["geometry"] = False
            entry["error"] = f"{type(exc).__name__}: {exc}"
        out[name] = entry
    return out


# ---------------------------------------------------------------------------
# 3. BATCHED FFT: stack B trials into (B, N, N) and run the hops in one call
# ---------------------------------------------------------------------------

def _iiij(N):
    """The alternating sign pattern that replaces a double fftshift."""
    ii = np.ones((N,), dtype=float)
    ii[1::2] = -1.0
    return np.outer(ii, ii)


def _transfer(dz, lam, size, N):
    """The Forvard transfer function CC and the piston, for one hop.

    This is exactly the array olb.waveoptics.propagators.Forvard builds for a
    non-negative z (Schmidt, DOI 10.1117/3.866274, Ch. 6, Eq. (6.32), p. 95).
    The batched hop reuses it across the whole (B, N, N) stack.
    """
    _2pi = 2.0 * 3.141592654
    z = abs(dz)
    kz = _2pi / lam * z
    z1 = z * lam / 2.0
    No2 = int(N / 2)
    SW = np.arange(-No2, N - No2) / size
    SW = SW * SW
    SSW = SW.reshape((-1, 1)) + SW
    Bus = z1 * SSW
    Ir = Bus.astype(int)
    Abus = _2pi * (Ir - Bus)
    CC = np.cos(Abus) + 1j * np.sin(Abus)
    return CC, complex(np.cos(kz), np.sin(kz))


def _batched_forvard(arr, dz, lam, size, N, iiij, backend, workers):
    """Propagate a (B, N, N) stack one hop, the batched Forvard.

    It matches Forvard for dz >= 0: sign flip, fft2, transfer, ifft2, piston,
    sign flip. numpy vectorises the leading axis; scipy adds a workers pool.
    """
    if dz == 0.0:
        return arr
    CC, piston = _transfer(dz, lam, size, N)
    arr = arr * iiij
    if backend == "numpy":
        arr = np.fft.fft2(arr, axes=(-2, -1))
        arr *= CC
        arr = np.fft.ifft2(arr, axes=(-2, -1))
    else:
        arr = _sfft.fft2(arr, axes=(-2, -1), workers=workers)
        arr *= CC
        arr = _sfft.ifft2(arr, axes=(-2, -1), workers=workers)
    arr = arr * piston
    arr *= iiij
    return arr


def _batch_screens(setup, entropy, ks, j):
    """Make one screen for plane j of each trial in ks, stacked (B, N, N).

    The screens are made on the fly, one plane at a time, to cap memory. They
    use the SAME per-(trial, screen) seed as the serial run, so a batched trial
    is bit-comparable with the serial trial of the same index.
    """
    plan = setup["plan"]
    grid = setup["grid"]
    r0 = plan.r0_m[j]
    build = setup["build_screen"]
    return np.stack([build(_screen_seed(entropy, k, j), r0) for k in ks])


def batched_split_step(setup, entropy, ks, backend, workers):
    """Run a batch of trials through the split step in (B, N, N) arrays.

    It mirrors olb.waveoptics.turbulence.splitstep.split_step: hop with
    sub-steps and a mask after each, apply a screen, mask, and a final hop. The
    only change is the leading batch axis. Returns the (B, N, N) receive field.
    """
    plan = setup["plan"]
    grid = setup["grid"]
    N = grid.n
    lam = setup["lam"]
    size = grid.size_m
    mask = setup["mask"]
    max_step_m = setup["max_step_m"]
    iiij = _iiij(N)
    B = len(ks)

    if setup["is_space"]:
        start = np.ones((N, N), dtype=np.complex128)
    else:
        start = setup["F_in"].field
    arr = np.broadcast_to(start, (B, N, N)).astype(np.complex128)

    def hop(a, gap_m):
        for dz in _substeps(gap_m, max_step_m):
            a = _batched_forvard(a, dz, lam, size, N, iiij, backend, workers)
            a = a * mask
        return a

    here = 0.0
    z = np.asarray(plan.z_m, dtype=float).ravel()
    for j, zi in enumerate(z):
        arr = hop(arr, zi - here)
        screens = _batch_screens(setup, entropy, ks, j)
        arr = arr * np.exp(1j * screens)
        arr = arr * mask
        here = zi
    arr = hop(arr, plan.z_total_m - here)
    return arr


def _batched_collected(setup, arr):
    """Collected_power of each slice of a batched receive field."""
    template = Begin(setup["grid"].size_m, setup["lam"], setup["grid"].n)
    out = np.empty(arr.shape[0])
    for b in range(arr.shape[0]):
        F = Field.copy(template)
        F.field = arr[b]
        F._IsGauss = False
        out[b] = float(Power(_clip(F, setup["aperture_m"],
                                   setup["obscuration"])) / setup["p_reference"])
    return out


def batched_scaling(name, b_list, generator="olb"):
    """Sweep the batch size B and the FFT backend for one case.

    trials_per_s is B / (the split-step time of the batch). The screens are made
    on the fly per hop, so the memory holds one field stack and one screen stack.
    """
    setup = _attach_builder(build_setup(name, generator))
    entropy = _resolve_seed(SEED)
    N = setup["grid"].n
    cpu = os.cpu_count() or 1
    backends = [("numpy", 1), ("scipy_w1", 1), (f"scipy_w{cpu}", cpu)]
    rows = []
    for B in b_list:
        ks = list(range(B))
        row = {"B": B}
        for label, w in backends:
            backend = "numpy" if label == "numpy" else "scipy"
            # warm once (screen cache, fft plan), then time.
            _ = batched_split_step(setup, entropy, ks[:1], backend, w)
            t0 = time.perf_counter()
            arr = batched_split_step(setup, entropy, ks, backend, w)
            wall = time.perf_counter() - t0
            row[label + "_s"] = wall
            row[label + "_trials_per_s"] = B / wall
        rows.append(row)
    # a memory budget, honest about the screen stack too.
    field_mb = N * N * 16 / 1e6                  # one complex128 field
    screen_mb = N * N * 8 / 1e6                  # one float64 screen
    mem = {
        "grid_n": N, "n_screens": setup["n_screens"],
        "field_mb_each": field_mb, "screen_mb_each": screen_mb,
        "note": ("peak ~ B*(field + exp(screen)) + B*screen while a plane is "
                 "applied; screens are made one plane at a time"),
        "peak_mb_per_B": 2.0 * field_mb + screen_mb,
        "example_N2048_B32_gb": 32 * (2.0 * (2048 ** 2 * 16 / 1e6)
                                      + (2048 ** 2 * 8 / 1e6)) / 1e3,
    }
    return {"case": name, "generator": generator, "grid_n": N, "rows": rows,
            "memory": mem}


# ---------------------------------------------------------------------------
# 4. the seed-contract cross-checks
# ---------------------------------------------------------------------------

def crosscheck(name, generator, n=8):
    """Every route must give the SAME trials as the serial loop for one seed."""
    scn, geom, hs, cn2, preset_name = _build_case(name)
    entropy = _resolve_seed(SEED)
    setup = _attach_builder(build_setup(name, generator))

    # (a) the replicated serial must equal the PRODUCTION serial, bit for bit.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        prod = propagate_turbulent_scenario(
            scn, geom, n_trials=n, seed=SEED, preset=preset_name, hs=hs,
            cn2_profile=cn2, screen_generator=generator)
    prod_cp = np.array([t.collected_power for t in prod.trials])
    repl_cp = serial_collected(setup, entropy, n)
    d_repl = float(np.max(np.abs(repl_cp - prod_cp)))

    # (b) the production Threader must equal the production serial, bit for bit.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        thr = propagate_turbulent_scenario(
            scn, geom, n_trials=n, seed=SEED, preset=preset_name, hs=hs,
            cn2_profile=cn2, threader=Threader(max_workers=4),
            screen_generator=generator)
    thr_cp = np.array([t.collected_power for t in thr.trials])
    d_thread = float(np.max(np.abs(thr_cp - prod_cp)))

    # (c) the ProcessPool must equal the production serial, bit for bit.
    with ProcessPoolExecutor(max_workers=2, initializer=_proc_init,
                             initargs=(name, generator)) as pool:
        proc_cp = np.array(list(pool.map(_proc_trial, range(n))))
    d_proc = float(np.max(np.abs(proc_cp - prod_cp)))

    # (d) the batched split step must match the serial, to FFT round-off. It
    # uses the same per-(trial, screen) seeds, so this is a tight tolerance.
    arr = batched_split_step(setup, entropy, list(range(n)), "numpy", 1)
    batch_cp = _batched_collected(setup, arr)
    d_batch = float(np.max(np.abs(batch_cp - prod_cp)
                           / np.maximum(prod_cp, 1e-30)))

    return {
        "name": name, "generator": generator, "n": n,
        "serial_repl_vs_production_max_abs": d_repl,
        "thread_vs_serial_max_abs": d_thread,
        "process_vs_serial_max_abs": d_proc,
        "batched_vs_serial_max_rel": d_batch,
        "serial_repl_exact": d_repl == 0.0,
        "thread_exact": d_thread == 0.0,
        "process_exact": d_proc == 0.0,
        "batched_within_1e-9": d_batch < 1e-9,
    }


# ---------------------------------------------------------------------------
# the recommendation
# ---------------------------------------------------------------------------

def _best_mode(thread_rows, proc_rows, batch_res, generator):
    """Pick the fastest route for one (case, generator)."""
    best = {"mode": "threads", "count": thread_rows[0]["workers"],
            "trials_per_s": thread_rows[0]["trials_per_s"]}
    for r in thread_rows:
        if r["trials_per_s"] > best["trials_per_s"]:
            best = {"mode": "threads", "count": r["workers"],
                    "trials_per_s": r["trials_per_s"]}
    for r in proc_rows:
        if r["trials_per_s"] > best["trials_per_s"]:
            best = {"mode": "processes", "count": r["workers"],
                    "trials_per_s": r["trials_per_s"]}
    if batch_res is not None and batch_res["generator"] == generator:
        cpu = os.cpu_count() or 1
        key = f"scipy_w{cpu}_trials_per_s"
        for r in batch_res["rows"]:
            best_batch = max(r.get("numpy_trials_per_s", 0),
                             r.get(key, 0),
                             r.get("scipy_w1_trials_per_s", 0))
            if best_batch > best["trials_per_s"]:
                best = {"mode": f"batched(B={r['B']})", "count": r["B"],
                        "trials_per_s": best_batch}
    one_worker = thread_rows[0]["trials_per_s"]
    best["speedup_over_1"] = best["trials_per_s"] / one_worker
    best["one_worker_trials_per_s"] = one_worker
    return best


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main():
    t_start = time.time()
    warm_aotools()
    workers = _worker_list()
    cpu = os.cpu_count() or 1

    env = {
        "numpy": np.__version__,
        "scipy": __import__("scipy").__version__,
        "aotools": __import__("aotools").__version__,
        "cpu": platform.processor(),
        "cores": cpu,
        "platform": platform.platform(),
        "trials_per_point": TRIALS,
        "worker_list": workers,
    }
    print("environment:")
    for k, v in env.items():
        print(f"  {k:20s} {v}")
    print("")

    # ---- the raw fft2 memory-bandwidth reference (from P0, re-read) ----
    p0_path = os.path.join(HERE, "profile_baseline_results.json")
    fft_ref = {}
    if os.path.exists(p0_path):
        with open(p0_path) as f:
            fft_ref = json.load(f).get("fft_microbench", {})

    # ---- the scenario-pickle check ----
    print("scenario pickle round-trip:")
    picks = scenario_pickles()
    for name, e in picks.items():
        ok = e["scenario"] and e["geometry"]
        print(f"  {name:34s} {'OK' if ok else 'FAIL ' + str(e['error'])}")
    print("")

    # ---- 1 + 2: threads and processes, every case, both generators ----
    threads = {}
    procs = {}
    for name in SCALING_CASES:
        for gen in GENERATORS:
            key = f"{name} | {gen}"
            print(f"threads:   {key} ...", flush=True)
            threads[key] = thread_scaling(name, gen, workers)
            print(f"processes: {key} ...", flush=True)
            procs[key] = process_scaling(name, gen, workers)

    # ---- 3: the batched prototype, one case, both generators ----
    batched = {}
    b_list = [b for b in (1, 2, 4, 8, 16, 32) if b <= 32]
    for gen in GENERATORS:
        print(f"batched:   {BATCH_CASE} | {gen} ...", flush=True)
        batched[gen] = batched_scaling(BATCH_CASE, b_list, generator=gen)

    # ---- 4: the seed-contract cross-check ----
    print("cross-check (seed contract) ...", flush=True)
    checks = {}
    for name in SCALING_CASES:
        for gen in GENERATORS:
            checks[f"{name} | {gen}"] = crosscheck(name, gen)

    # ---- the printed tables ----
    print("")
    print("=" * 78)
    print("1. THREADS  (trials/s, speed-up over 1 worker, parallel efficiency)")
    print("=" * 78)
    for key, rows in threads.items():
        print(f"  {key}")
        print(f"     {'W':>4}{'trials/s':>12}{'speedup':>10}{'eff':>8}")
        for r in rows:
            print(f"     {r['workers']:>4}{r['trials_per_s']:>12.2f}"
                  f"{r['speedup']:>10.2f}{r['efficiency']:>8.2f}")

    print("")
    print("=" * 78)
    print("2. PROCESSES  (spawn+init apart from the steady-state rate)")
    print("=" * 78)
    for key, rows in procs.items():
        print(f"  {key}")
        print(f"     {'W':>4}{'spawn_s':>10}{'steady/s':>12}"
              f"{'speedup':>10}{'eff':>8}")
        for r in rows:
            print(f"     {r['workers']:>4}{r['spawn_init_s']:>10.2f}"
                  f"{r['trials_per_s']:>12.2f}{r['speedup']:>10.2f}"
                  f"{r['efficiency']:>8.2f}")

    print("")
    print("=" * 78)
    print(f"3. BATCHED FFT  ({BATCH_CASE}, N={batched['olb']['grid_n']})")
    print("=" * 78)
    for gen in GENERATORS:
        br = batched[gen]
        print(f"  generator {gen}: trials/s")
        print(f"     {'B':>4}{'numpy':>12}{'scipy_w1':>12}"
              f"{'scipy_w' + str(cpu):>12}")
        for r in br["rows"]:
            print(f"     {r['B']:>4}{r['numpy_trials_per_s']:>12.2f}"
                  f"{r['scipy_w1_trials_per_s']:>12.2f}"
                  f"{r[f'scipy_w{cpu}_trials_per_s']:>12.2f}")
    mem = batched["olb"]["memory"]
    print(f"  memory: one field {mem['field_mb_each']:.2f} MB, one screen "
          f"{mem['screen_mb_each']:.2f} MB at N={mem['grid_n']}; "
          f"peak ~ {mem['peak_mb_per_B']:.2f} MB per trial in the batch.")
    print(f"          at N=2048, B=32 the batch needs about "
          f"{mem['example_N2048_B32_gb']:.2f} GB.")

    print("")
    print("=" * 78)
    print("4. SEED-CONTRACT CROSS-CHECK  (max deviation from the serial loop)")
    print("=" * 78)
    all_ok = True
    for key, c in checks.items():
        ok = (c["thread_exact"] and c["process_exact"]
              and c["batched_within_1e-9"] and c["serial_repl_exact"])
        all_ok = all_ok and ok
        print(f"  {key}")
        print(f"     serial-repl vs prod {c['serial_repl_vs_production_max_abs']:.2e}"
              f"   thread {c['thread_vs_serial_max_abs']:.2e}"
              f"   process {c['process_vs_serial_max_abs']:.2e}"
              f"   batched(rel) {c['batched_vs_serial_max_rel']:.2e}"
              f"   {'PASS' if ok else 'FAIL'}")
    # A hard gate: the parallel routes must reproduce the serial loop.
    assert all_ok, "a parallel route broke the seed contract"

    # ---- the RECOMMENDATION table ----
    print("")
    print("=" * 78)
    print("RECOMMENDATION  (best mode, worker/batch count, speed-up over 1)")
    print("=" * 78)
    recommendation = {}
    proc_beats_thread = {}
    for name in SCALING_CASES:
        for gen in GENERATORS:
            key = f"{name} | {gen}"
            best = _best_mode(threads[key], procs[key],
                              batched.get(gen), gen)
            recommendation[key] = best
            # does the best process rate beat the best thread rate?
            best_thr = max(r["trials_per_s"] for r in threads[key])
            best_prc = max(r["trials_per_s"] for r in procs[key])
            proc_beats_thread[key] = {
                "best_thread_trials_per_s": best_thr,
                "best_process_trials_per_s": best_prc,
                "process_over_thread": best_prc / best_thr}
    print(f"  {'case | generator':<40}{'mode':>16}{'count':>7}{'x1':>7}")
    for key, best in recommendation.items():
        print(f"  {key:<40}{best['mode']:>16}{best['count']:>7}"
              f"{best['speedup_over_1']:>7.1f}")
    print("")
    print("  processes vs threads (best rate ratio, >1 means processes win):")
    for key, d in proc_beats_thread.items():
        print(f"     {key:<40}{d['process_over_thread']:>6.2f}x")

    # a one-line thread-default read, from the aotools standard downlink curve.
    default_key = f"space downlink 30deg standard | aotools"
    trows = threads.get(default_key, [])
    if trows:
        knee = max(trows, key=lambda r: r["trials_per_s"])
        eff_at_cpu = trows[-1]["efficiency"]
        print("")
        print(f"  Threader default read ({default_key}): peak trials/s at "
              f"W={knee['workers']} (efficiency {knee['efficiency']:.2f}); "
              f"at W={trows[-1]['workers']} efficiency is {eff_at_cpu:.2f}.")

    out = {
        "environment": env,
        "fft_microbench_ref_from_p0": fft_ref,
        "scenario_pickles": picks,
        "threads": threads,
        "processes": procs,
        "batched": batched,
        "crosschecks": checks,
        "recommendation": recommendation,
        "process_vs_thread": proc_beats_thread,
    }
    path = os.path.join(HERE, "scaling_study_results.json")
    with open(path, "w") as f:
        json.dump(out, f, indent=1)
    print("")
    print(f"wrote {os.path.relpath(path)}")
    print(f"(elapsed {time.time() - t_start:.1f} s)")


if __name__ == "__main__":
    main()
