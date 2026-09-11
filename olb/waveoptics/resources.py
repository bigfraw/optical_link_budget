"""The machine resources of a wave-optics run: the free memory, the memory
one pool worker needs, and the worker count that fills the machine.

The fidelity-2 Monte Carlo is bound by memory: by the BANDWIDTH in steady
state, and by the CAPACITY when a large pool runs a large grid (a 16-worker
pool at 1024 px and 15 screens ran out of memory on 2026-09-04, see
validation/tail_convergence/README.md). So the pool must know two numbers
before it starts: the free memory of the box, and the memory of one worker.
This module gives both, and the worker count that keeps a target fraction of
the cores busy without going past the free memory.

The module imports the standard library and numpy only.
"""
import os
import sys

import numpy as np

# The steady memory of one Python process that has imported numpy, scipy and
# olb, before it allocates a grid. Measured on Windows and Linux at 120 to
# 150 MB; the margin covers the allocator.
WORKER_BASE_BYTES = 160 * 2 ** 20

# The safety factor on the per-worker estimate. The Windows allocator keeps
# a high-water mark, and a transient can stack on another one.
WORKER_SAFETY = 1.25


def free_memory_bytes():
    """Give the memory the machine can give to new processes now, in bytes.

    On Linux it is `MemAvailable` of /proc/meminfo. On Windows it is
    `ullAvailPhys` of GlobalMemoryStatusEx. On macOS it is the free plus the
    inactive pages of `vm_stat`.

    Returns:
        An int, or None when the probe is not available.
    """
    try:
        if sys.platform.startswith("linux"):
            with open("/proc/meminfo") as fh:
                for line in fh:
                    if line.startswith("MemAvailable:"):
                        return int(line.split()[1]) * 1024
            return None
        if sys.platform.startswith("win"):
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
            if not ctypes.WinDLL("kernel32").GlobalMemoryStatusEx(
                    ctypes.byref(st)):
                return None
            return int(st.ullAvailPhys)
        if sys.platform == "darwin":
            import subprocess
            out = subprocess.run(["vm_stat"], capture_output=True,
                                 text=True, check=True).stdout
            page = 4096
            pages = 0
            for line in out.splitlines():
                if "page size of" in line:
                    page = int(line.split("page size of")[1].split()[0])
                for tag in ("Pages free:", "Pages inactive:"):
                    if line.startswith(tag):
                        pages += int(line.split(":")[1].strip().rstrip("."))
            return pages * page
    except (OSError, ValueError, ImportError):
        return None
    return None


def worker_memory_bytes(n, precision="single", block_size=1, patch_pixels=0,
                        n_hops=None, screen_n=None, n_screens=0):
    """Estimate the peak memory of ONE pool worker of a turbulent run.

    The estimate counts the arrays that a trial holds at its peak, on an
    N x N grid:

    - the field and its copies in the split step: the input, the hop copy,
      the FFT work array and its two FFT results, and the mask product,
      about 6 complex arrays;
    - ONE phase screen (the runner makes them one at a time) plus the
      generator transients, about 4 complex arrays at the screen precision;
    - the boundary mask and the sign pattern, 2 real arrays;
    - the Forvard factor cache: one transfer function for each distinct hop
      (a split-step plan makes one distinct hop for each screen gap, so
      `n_hops` is the screen count plus one), capped at the cache bound;
    - the stored field patch of the block, two times (the block array and
      its pickled copy on the way back to the parent);
    - the KEPT SCREEN NOISE of a POINT-AHEAD trial: that trial keeps the white
      noise of every screen, so each pass rebuilds the same oversize screen on
      its own window. The term is `n_screens * screen_n^2` complex128 values
      (the draw is always double, see ScreenFactory.draw). It is 0 by default,
      so a run with no point-ahead pass reads the number it always read;
    - the base of the interpreter with numpy, scipy and olb.

    The whole is scaled by WORKER_SAFETY. The estimate sits above the measured
    working set on purpose: a pool that starts too many workers is killed,
    and a pool that starts one too few loses a few percent.

    Args:
        n:            the pixel count of one grid side.
        precision:    "single" or "double".
        block_size:   the trials of one block.
        patch_pixels: the pixel count of the stored patch, 0 for none.
        n_hops:       the distinct hop lengths of the plan. None takes the
                      full cache bound.
        screen_n:     the pixel count of one drawn screen side. None (the
                      default) takes `n`, the screen of record.
        n_screens:    the screens whose noise a POINT-AHEAD trial keeps. 0 (the
                      default) keeps the old estimate.

    Returns:
        An int, in bytes.
    """
    from .propagators import FORVARD_CACHE_BYTES
    c = 8 if precision == "single" else 16
    r = c // 2
    n2 = int(n) * int(n)
    arrays = n2 * (10 * c + 2 * r)
    cache_bytes = FORVARD_CACHE_BYTES
    if n_hops is not None:
        cache_bytes = min(cache_bytes, int(n_hops) * n2 * c)
    patch = 2 * int(block_size) * int(patch_pixels) * 8      # complex64.
    # The kept noise of one POINT-AHEAD trial, in complex128 (ScreenFactory
    # draws in double). It is 0 when no trial keeps it.
    side = int(n if screen_n is None else screen_n)
    noise = int(n_screens) * side * side * 16
    total = WORKER_BASE_BYTES + arrays + int(cache_bytes) + patch + noise
    return int(np.ceil(total * WORKER_SAFETY))


def auto_workers(per_worker_bytes, *, cpu_fraction=0.9, memory_fraction=0.9,
                 free_bytes=None, cores=None):
    """Give the worker count that fills the machine.

    The count is the smaller of two limits: the cores, at `cpu_fraction`
    of the logical core count (the rest is left to the parent, the disk
    writes and the shell), and the memory, at `memory_fraction` of the free
    memory divided by the memory of one worker. The result is at least 1.

    The CPU limit is a TARGET, not a measurement. A memory-bandwidth-bound
    grid can plateau under it (12 workers of 32 threads on bigfraw at 512 px
    before the 2026-09-06 cache and lazy screens; measure again with
    validation/campaign_resources/ after a change). The memory limit is a
    hard one: past it a worker dies with MemoryError.

    Args:
        per_worker_bytes: the memory of one worker, from
                          worker_memory_bytes().
        cpu_fraction:     the fraction of the logical cores to use.
        memory_fraction:  the fraction of the free memory to use.
        free_bytes:       the free memory. None probes the machine; a
                          machine with no probe takes the CPU limit only.
        cores:            the logical core count. None reads os.cpu_count().

    Returns:
        The pair (workers, reason): an int and the string "cpu" or
        "memory" that names the binding limit.
    """
    if cores is None:
        cores = os.cpu_count() or 1
    by_cpu = max(1, int(np.floor(cpu_fraction * cores)))
    if free_bytes is None:
        free_bytes = free_memory_bytes()
    if free_bytes is None or per_worker_bytes <= 0:
        return by_cpu, "cpu"
    by_mem = max(1, int(np.floor(memory_fraction * free_bytes
                                 / per_worker_bytes)))
    if by_mem < by_cpu:
        return by_mem, "memory"
    return by_cpu, "cpu"


if __name__ == "__main__":
    free = free_memory_bytes()
    print(f"free memory: {free / 2 ** 30:.2f} GiB" if free else
          "free memory: no probe")
    for n in (256, 512, 1024, 2048):
        for prec in ("single", "double"):
            w = worker_memory_bytes(n, prec, block_size=50,
                                    patch_pixels=int(np.pi * (n / 10) ** 2),
                                    n_hops=10)
            k, why = auto_workers(w)
            print(f"  {n:5d} px {prec:6s}: {w / 2 ** 20:7.0f} MiB per worker"
                  f" -> {k:3d} workers ({why})")
    # Checks. The memory limit binds when the memory is small, the CPU limit
    # when it is large, and the count never goes under 1.
    assert auto_workers(2 ** 30, free_bytes=2 ** 30, cores=8) == (1, "memory")
    assert auto_workers(2 ** 20, free_bytes=2 ** 40, cores=8) == (7, "cpu")
    assert auto_workers(2 ** 40, free_bytes=2 ** 30, cores=8)[0] == 1
    assert auto_workers(2 ** 20, free_bytes=None, cores=8) == (7, "cpu")
    assert worker_memory_bytes(1024, "single") < worker_memory_bytes(
        1024, "double")
    # The POINT-AHEAD noise term is OPT-IN: the default reads the old number,
    # and a kept stack of oversize screens raises it.
    assert worker_memory_bytes(512, "single") == worker_memory_bytes(
        512, "single", n_screens=0), "the default must not move"
    assert worker_memory_bytes(512, "single", screen_n=576, n_screens=9) > \
        worker_memory_bytes(512, "single")
    print("resources self-check ok")
