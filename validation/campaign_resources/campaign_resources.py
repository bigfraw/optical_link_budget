"""A 4000-trial fidelity-2 downlink Campaign with a CPU and RAM monitor.

WHAT IT DOES. It runs ONE fidelity-2 space downlink through the production
`Campaign` store (olb.waveoptics.turbulence.campaign): a 700 mm ground
aperture with a 30 percent central obscuration, an SMF detector, 1550 nm, a
500 km orbit at 30 deg, the standard preset, and the fixed outer scale
L0 = 25 m (the owner decision of 2026-09-05). The trials go to disk in blocks,
so a killed run resumes with no rerun.

WHY THE MONITOR. The user wants proof that the parallel work is real. A
Campaign has ONE level of parallelism: `workers=W` opens one warm process pool
and runs each block SERIALLY inside its process; `workers=None` runs the
blocks one after the other, each block THREADED inside (the Threader of the
runner). A background thread of this script samples the machine every
`--sample-s` seconds and writes a CSV:

    t_s, cpu_pct, cores_busy, ram_used_mb, ram_avail_mb, n_python, python_rss_mb, cpu_mhz

- `cpu_pct` and `cores_busy` come from the Windows `GetSystemTimes` deltas
  (idle, kernel, user; system wide, all logical cores). `cores_busy` is
  the busy fraction times the logical core count, so it reads directly
  against `workers`.
- `ram_used_mb` and `ram_avail_mb` come from `GlobalMemoryStatusEx`.
- `n_python` and `python_rss_mb` come from `tasklist`: the count of
  python.exe processes and the sum of their working sets. With `workers=W`
  the count must read W + 1 (the pool plus this parent) while blocks run.

- `cpu_mhz` is the mean current clock of the logical cores
  (`CallNtPowerInformation`). CAUTION: on the i9-14900HX it reads a
  constant nominal 1900 MHz, so it is NOT a throttling probe there. Judge
  the clocks with the `% Processor Utility` counter (frequency weighted,
  the Task Manager number) against `% Processor Time`: utility above
  time means the cores boost, utility below time means they throttle.

THE WORKER BOOST. Windows does not pass the power-throttling opt-out of the
parent to a spawned child, so the script wraps the Campaign pool initializer
and every worker calls `boost_process_priority()` itself one time.

No psutil: the env does not carry it, so the monitor is ctypes and one
Windows tool only. Off Windows the CPU and RAM columns read NaN.

THE VERIFICATION RULE (docs/api-waveoptics.md Section 9g). The effective
process count is `min(workers, ceil(n_trials / block_size))`, so keep
`block_size <= n_trials / workers`. The default 4000 trials in blocks of 50
gives 80 blocks, enough for any worker count up to 80.

Usage (on the desktop, after the sync script):

    python -m validation.campaign_resources.campaign_resources
    python -m validation.campaign_resources.campaign_resources --workers 24
    python -m validation.campaign_resources.campaign_resources --threads
    python -m validation.campaign_resources.campaign_resources --workers 12 --grid-n 1024
    python -m validation.campaign_resources.campaign_resources --smoke

`--smoke` is a small local check: 8 trials, blocks of 2, 2 workers, rapid.
"""

import argparse
import csv
import json
import os
import subprocess
import sys
import threading
import time

import numpy as np

from olb.geometry import CircularOrbit
from olb.scenario import Channel, Site, SpaceScenario
from olb.terminal import SMF, Terminal, Transmitter
from olb.waveoptics.grid import GridSpec
from olb.waveoptics.turbulence import Campaign
from olb.waveoptics.turbulence.sampling import turbulent_grid
from olb.waveoptics.turbulence import campaign as _campaign_module

HERE = os.path.dirname(os.path.abspath(__file__))
LAM = 1550e-9
SEED = 20260905
L0_M = 25.0


def boost_process_priority(high=False):
    """Windows only: raise the priority class of THIS process and opt it out of
    power throttling (EcoQoS). Pure ctypes, no dependency. No-op off Windows,
    never raises.

    The default is ABOVE_NORMAL. HIGH (`high=True`) puts a 16-worker pool above
    sshd and the VS Code server, and a Remote-SSH connection then times out
    (seen 2026-09-05). The EcoQoS opt-out, not the class, is what stops the
    throttling, so ABOVE_NORMAL keeps the speed and keeps the machine usable."""
    if not sys.platform.startswith("win"):
        return
    try:
        import ctypes
        from ctypes import wintypes
        k32 = ctypes.WinDLL("kernel32", use_last_error=True)
        k32.GetCurrentProcess.restype = wintypes.HANDLE
        k32.SetPriorityClass.argtypes = [wintypes.HANDLE, wintypes.DWORD]
        k32.SetPriorityClass.restype = wintypes.BOOL
        hproc = k32.GetCurrentProcess()
        HIGH_PRIORITY_CLASS = 0x00000080
        ABOVE_NORMAL_PRIORITY_CLASS = 0x00008000
        k32.SetPriorityClass(hproc, HIGH_PRIORITY_CLASS if high else ABOVE_NORMAL_PRIORITY_CLASS)

        class _PPTS(ctypes.Structure):
            _fields_ = [("Version", wintypes.ULONG),
                        ("ControlMask", wintypes.ULONG),
                        ("StateMask", wintypes.ULONG)]
        k32.SetProcessInformation.argtypes = [wintypes.HANDLE, ctypes.c_int,
                                              ctypes.c_void_p, wintypes.DWORD]
        k32.SetProcessInformation.restype = wintypes.BOOL
        st = _PPTS(1, 0x1, 0)
        k32.SetProcessInformation(hproc, 4, ctypes.byref(st), ctypes.sizeof(st))
    except Exception:
        pass


_ORIGINAL_INIT_WORKER = _campaign_module._init_worker


def _boosted_init_worker(payload):
    """The pool initializer: boost THIS worker, then run the Campaign one.

    Windows does not pass the power-throttling opt-out to a spawned child, so
    a worker of the pool runs parked and downclocked unless it opts out
    itself. This wrapper runs in each worker one time.
    """
    boost_process_priority()
    _ORIGINAL_INIT_WORKER(payload)


# ---------------------------------------------------------------------------
# The scenario
# ---------------------------------------------------------------------------

def scenario_and_geometry(elevation_deg, aperture_m=0.7, obscuration=0.3):
    """Build the downlink scenario: an obscured ground SMF receiver.

    Args:
        elevation_deg: the elevation of the one line of sight, in deg.
        aperture_m:    the ground receive aperture diameter, in m.
        obscuration:   the central obscuration ratio of that aperture.

    Returns:
        The pair (SpaceScenario, CircularOrbit).
    """
    site = Site(cn2_ground=1.7e-14, wind_rms_m_s=21.0)
    channel = Channel(site=site, altitude_m=500e3)
    ground = Terminal(aperture_m=aperture_m, obscuration_ratio=obscuration,
                      wavelength_m=LAM, pointing_jitter_rad=2e-6,
                      detector=SMF(sensitivity_dbm=-45.0))
    space = Terminal(aperture_m=0.10, wavelength_m=LAM,
                     pointing_jitter_rad=1e-6,
                     transmitter=Transmitter(waist_m=0.04, power_dbm=30.0))
    scn = SpaceScenario(ground=ground, space=space, direction="downlink",
                        channel=channel)
    geom = CircularOrbit(altitude_m=500e3, elevation_deg=[float(elevation_deg)])
    return scn, geom


# ---------------------------------------------------------------------------
# The resource monitor
# ---------------------------------------------------------------------------

def _system_times():
    """Give (idle, kernel, user) in 100 ns ticks from GetSystemTimes, or None."""
    if not sys.platform.startswith("win"):
        return None
    import ctypes
    from ctypes import wintypes
    k32 = ctypes.WinDLL("kernel32")
    idle, kern, user = wintypes.FILETIME(), wintypes.FILETIME(), wintypes.FILETIME()
    if not k32.GetSystemTimes(ctypes.byref(idle), ctypes.byref(kern),
                              ctypes.byref(user)):
        return None

    def tick(ft):
        return (ft.dwHighDateTime << 32) | ft.dwLowDateTime
    # The kernel time INCLUDES the idle time on Windows.
    return tick(idle), tick(kern), tick(user)


def _memory_mb():
    """Give (used_mb, avail_mb) from GlobalMemoryStatusEx, or (nan, nan)."""
    if not sys.platform.startswith("win"):
        return np.nan, np.nan
    import ctypes
    from ctypes import wintypes

    class MEMORYSTATUSEX(ctypes.Structure):
        _fields_ = [("dwLength", wintypes.DWORD),
                    ("dwMemoryLoad", wintypes.DWORD),
                    ("ullTotalPhys", ctypes.c_ulonglong),
                    ("ullAvailPhys", ctypes.c_ulonglong),
                    ("ullTotalPageFile", ctypes.c_ulonglong),
                    ("ullAvailPageFile", ctypes.c_ulonglong),
                    ("ullTotalVirtual", ctypes.c_ulonglong),
                    ("ullAvailVirtual", ctypes.c_ulonglong),
                    ("ullAvailExtendedVirtual", ctypes.c_ulonglong)]
    st = MEMORYSTATUSEX()
    st.dwLength = ctypes.sizeof(st)
    ctypes.WinDLL("kernel32").GlobalMemoryStatusEx(ctypes.byref(st))
    total, avail = st.ullTotalPhys / 2**20, st.ullAvailPhys / 2**20
    return total - avail, avail


def _cpu_mhz():
    """Give the mean current clock of the logical cores in MHz, or nan.

    It calls CallNtPowerInformation(ProcessorInformation). A throttled
    (parked, downclocked) run reads a low clock while the busy fraction
    reads high, so this column separates the two.
    """
    if not sys.platform.startswith("win"):
        return np.nan
    import ctypes
    from ctypes import wintypes

    class PPI(ctypes.Structure):
        _fields_ = [("Number", wintypes.ULONG), ("MaxMhz", wintypes.ULONG),
                    ("CurrentMhz", wintypes.ULONG), ("MhzLimit", wintypes.ULONG),
                    ("MaxIdleState", wintypes.ULONG),
                    ("CurrentIdleState", wintypes.ULONG)]
    n = os.cpu_count() or 1
    buf = (PPI * n)()
    rc = ctypes.WinDLL("powrprof").CallNtPowerInformation(
        11, None, 0, ctypes.byref(buf), ctypes.sizeof(buf))
    if rc != 0:
        return np.nan
    return float(np.mean([p.CurrentMhz for p in buf]))


def _python_processes():
    """Give (count, rss_mb) of every python.exe from tasklist, or (nan, nan)."""
    if not sys.platform.startswith("win"):
        return np.nan, np.nan
    try:
        out = subprocess.run(
            ["tasklist", "/FI", "IMAGENAME eq python.exe", "/FO", "CSV", "/NH"],
            capture_output=True, text=True, timeout=10).stdout
    except Exception:
        return np.nan, np.nan
    n, rss = 0, 0.0
    for row in csv.reader(out.splitlines()):
        if len(row) >= 5 and row[0].lower() == "python.exe":
            n += 1
            digits = "".join(ch for ch in row[4] if ch.isdigit())
            rss += float(digits or 0) / 1024
    return n, rss


class ResourceMonitor:
    """A background sampler of the system CPU, the RAM and the python processes.

    Use it as a context manager. It writes one CSV row for each sample, so a
    killed run keeps its samples.
    """

    def __init__(self, csv_path, sample_s=2.0):
        self.csv_path = csv_path
        self.sample_s = float(sample_s)
        self.ncpu = os.cpu_count() or 1
        self.rows = []
        self._stop = threading.Event()
        self._thread = threading.Thread(target=self._loop, daemon=True)

    def __enter__(self):
        with open(self.csv_path, "w", newline="", encoding="utf-8") as fh:
            csv.writer(fh).writerow(
                ["t_s", "cpu_pct", "cores_busy", "ram_used_mb", "ram_avail_mb",
                 "n_python", "python_rss_mb", "cpu_mhz"])
        self._thread.start()
        return self

    def __exit__(self, *exc):
        self._stop.set()
        self._thread.join(timeout=self.sample_s * 3)

    def _loop(self):
        t0 = time.time()
        prev = _system_times()
        while not self._stop.wait(self.sample_s):
            now = _system_times()
            if prev is None or now is None:
                cpu = np.nan
            else:
                d_idle = now[0] - prev[0]
                d_total = (now[1] - prev[1]) + (now[2] - prev[2])
                cpu = 100.0 * (1.0 - d_idle / d_total) if d_total > 0 else np.nan
            prev = now
            used, avail = _memory_mb()
            n_py, rss = _python_processes()
            def r(x, nd=0):
                # A NaN sample stays NaN: round() of a NaN raises.
                return x if not np.isfinite(x) else round(x, nd)
            row = [round(time.time() - t0, 1), r(cpu, 1),
                   r(cpu / 100.0 * self.ncpu, 2), r(used), r(avail),
                   n_py, r(rss), r(_cpu_mhz())]
            self.rows.append(row)
            with open(self.csv_path, "a", newline="", encoding="utf-8") as fh:
                csv.writer(fh).writerow(row)

    def summary(self):
        """Give the peak and the mean of each column over the run."""
        if not self.rows:
            return {}
        a = np.array(self.rows, dtype=float)
        names = ["cpu_pct", "cores_busy", "ram_used_mb", "ram_avail_mb",
                 "n_python", "python_rss_mb", "cpu_mhz"]
        out = {"n_samples": int(a.shape[0]), "logical_cores": self.ncpu}
        for i, name in enumerate(names, start=1):
            col = a[:, i]
            out[name] = {"mean": float(np.nanmean(col)),
                         "peak": float(np.nanmax(col))}
        return out


def plot(rows, ncpu, workers, path):
    """Draw the CPU, the RAM and the python process count against time."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    a = np.array(rows, dtype=float)
    t = a[:, 0]
    fig, ax = plt.subplots(3, 1, figsize=(9, 8), sharex=True)
    ax[0].plot(t, a[:, 2], label="cores busy (system)")
    if workers is not None:
        ax[0].axhline(workers, color="k", ls="--", label=f"workers = {workers}")
    ax[0].axhline(ncpu, color="gray", ls=":", label=f"logical cores = {ncpu}")
    ax[0].set_ylabel("busy cores")
    ax[0].legend(loc="lower right")
    ax[1].plot(t, a[:, 3] / 1024, label="RAM used (system)")
    ax[1].plot(t, a[:, 6] / 1024, label="python.exe RSS (sum)")
    ax[1].set_ylabel("GB")
    ax[1].legend(loc="lower right")
    ax[2].plot(t, a[:, 5], label="python.exe count")
    if workers is not None:
        ax[2].axhline(workers + 1, color="k", ls="--",
                      label=f"workers + parent = {workers + 1}")
    ax[2].set_ylabel("python.exe count")
    ax2 = ax[2].twinx()
    ax2.plot(t, a[:, 7], color="C3", label="mean clock [MHz]")
    ax2.set_ylabel("MHz")
    ax[2].legend(loc="lower left")
    ax2.legend(loc="lower right")
    ax[2].set_xlabel("time [s]")
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--n-trials", type=int, default=4000)
    ap.add_argument("--workers", type=int, default=16,
                    help="the process-pool size; --threads ignores it")
    ap.add_argument("--threads", action="store_true",
                    help="workers=None: serial blocks, threaded inside")
    ap.add_argument("--block-size", type=int, default=50)
    ap.add_argument("--elevation", type=float, default=30.0)
    ap.add_argument("--preset", default="standard")
    ap.add_argument("--sample-s", type=float, default=2.0)
    ap.add_argument("--grid-n", type=int, default=None,
                    help="pin the grid pixel count (the side and the screen "
                         "plan stay as the preset sizer gives them)")
    ap.add_argument("--root", default=None,
                    help="the campaign directory (default: next to this file)")
    ap.add_argument("--smoke", action="store_true",
                    help="8 trials, blocks of 2, 2 workers, rapid preset")
    args = ap.parse_args(argv)
    if args.smoke:
        args.n_trials, args.block_size, args.workers = 8, 2, 2
        args.preset = "rapid"

    boost_process_priority()
    _campaign_module._init_worker = _boosted_init_worker
    workers = None if args.threads else args.workers
    n_blocks = -(-args.n_trials // args.block_size)
    if workers is not None and n_blocks < workers:
        print(f"WARNING: only {n_blocks} blocks for {workers} workers; the pool "
              f"can use at most {n_blocks} processes. Lower --block-size.")

    grid_tag = "" if args.grid_n is None else f"_n{args.grid_n}"
    tag = (f"el{args.elevation:.0f}_{args.preset}{grid_tag}_"
           + ("threads" if workers is None else f"w{workers}")
           + ("_smoke" if args.smoke else ""))
    root = args.root or os.path.join(
        HERE, f"campaign_{args.preset}_el{args.elevation:.0f}{grid_tag}"
        + ("_smoke" if args.smoke else ""))
    csv_path = os.path.join(HERE, f"resources_{tag}.csv")
    json_path = os.path.join(HERE, f"resources_{tag}.json")
    png_path = os.path.join(HERE, f"resources_{tag}.png")

    scn, geom = scenario_and_geometry(args.elevation)
    grid = plan = None
    if args.grid_n is not None:
        # Size with the preset, then pin the pixel count only. The side and
        # the screen plan are the ones the preset sizer gives, so the pinned
        # run differs from the default run in the pixel pitch alone. Both
        # enter the campaign fingerprint.
        sized_grid, plan, _ = turbulent_grid(scn, geom, preset=args.preset,
                                             L0_m=L0_M)
        grid = GridSpec(size_m=sized_grid.size_m, n=int(args.grid_n),
                        scaled=sized_grid.scaled)
    camp = Campaign(scn, geom, root, seed=SEED, preset=args.preset,
                    block_size=args.block_size, L0_m=L0_M, grid=grid,
                    plan=plan)
    print(f"campaign      : {root}")
    print(f"scenario      : downlink 1550 nm, 500 km at {args.elevation:.0f} deg, "
          f"ground {scn.rx_terminal.aperture_m * 1e3:.0f} mm / "
          f"{scn.rx_terminal.obscuration_ratio:.0%} obscured, SMF")
    print(f"grid          : {camp.grid.n} px, {camp.grid.size_m:.3f} m, "
          f"{camp.plan.z_m.size} screens, L0 = {L0_M} m, preset {args.preset}")
    print(f"trials        : {args.n_trials} in blocks of {args.block_size} "
          f"({n_blocks} blocks), {camp.n_stored} already on disk")
    print("parallelism   : "
          + ("workers=None (serial blocks, threaded inside)" if workers is None
             else f"workers={workers} (one process pool, serial blocks inside)"))
    print(f"logical cores : {os.cpu_count()}      sample every {args.sample_s} s")

    n_before = camp.n_stored
    t0 = time.time()
    with ResourceMonitor(csv_path, sample_s=args.sample_s) as mon:
        stored = camp.run(args.n_trials, workers=workers, progress=True)
    wall = time.time() - t0

    summary = mon.summary()
    n_new = int(stored) - int(n_before)
    summary.update({"n_trials": int(stored), "n_new_trials": n_new,
                    "wall_s": wall,
                    # The rate counts the NEW trials only: a resumed campaign
                    # already holds the old blocks.
                    "trials_per_s": n_new / wall if wall > 0 else None,
                    "workers": workers, "block_size": args.block_size,
                    "grid_n": camp.grid.n, "n_screens": int(camp.plan.z_m.size),
                    "L0_m": L0_M, "preset": args.preset,
                    "campaign_root": root, "csv": csv_path})
    with open(json_path, "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)
    if mon.rows:
        plot(mon.rows, mon.ncpu, workers, png_path)

    print()
    print(f"done          : {stored} trials on disk, {n_new} new in {wall:.1f} s "
          f"({wall / max(n_new, 1):.3f} s/trial)")
    if summary.get("cores_busy"):
        print(f"cores busy    : mean {summary['cores_busy']['mean']:.1f}, "
              f"peak {summary['cores_busy']['peak']:.1f} of {mon.ncpu}")
        print(f"python procs  : peak {summary['n_python']['peak']:.0f}, "
              f"RSS peak {summary['python_rss_mb']['peak'] / 1024:.1f} GB")
        print(f"RAM used      : peak {summary['ram_used_mb']['peak'] / 1024:.1f} GB")
        print(f"clock         : mean {summary['cpu_mhz']['mean']:.0f} MHz")
    print(f"outputs       : {csv_path}\n                {json_path}\n"
          f"                {png_path}")


if __name__ == "__main__":
    main()
