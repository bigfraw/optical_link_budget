'''Gate: the float32 filter tables of a strip against the float64 tables.

THE QUESTION. `ScreenFactory` builds its win-1 filter (the radial frequency
grid and the von Karman PSD, Schmidt, DOI 10.1117/3.866274, Ch. 9, Eq. (9.51),
printed p. 161) as full-size tables. In float64 a crosswind box of 5e8 pixels
holds 4 GB in each table and the build swaps. `table_dtype=numpy.float32`
builds the tables in single precision. Is the screen the same?

THE MEASUREMENTS.

1. The filter of the default route equals the old meshgrid formula bit for
   bit, so the default screen of record did not move.
2. The same seed builds one strip two times, with float64 tables and with
   float32 tables. Report the maximum phase difference, the relative rms
   difference, and the structure function D(r) along y and along x of each
   route against the analytic von Karman law and against each other.
3. The cost of one large strip build on each route: the wall time, the
   tracemalloc peak and the process peak working set, each in a fresh
   process.

THE GATES. The phase difference must stay under 1e-3 rad (the float32 noise
floor of the transform is about 1e-5 rad on a 15 rad rms screen; the physics
reads 0.1 rad differences), and the two D(r) must agree inside 1e-3 on both
axes. A rounding-level change that passes them is a memory lever, not a
physics change.

Run from the repository root:

    python -m validation.temporal_screens.table_precision
    python -m validation.temporal_screens.table_precision --big 11897 43744

The second form measures the real 30 deg crosswind box (bigfraw).
'''
import argparse
import ctypes
import json
import os
import subprocess
import sys
import time
import tracemalloc

import numpy as np

from olb.waveoptics.propagators import set_fft_backend
from olb.waveoptics.turbulence.screens import ScreenFactory, _phase_psd_unit
from validation.temporal_screens.rect_factory import DX, L0, R0, R_LIST, \
    dphi_axes, theory_dphi

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, 'data')

MAX_DPHI_RAD = 1e-3
MAX_DPHI_RATIO = 1e-3


def peak_working_set_bytes():
    '''Give the peak working set of this process, in bytes (Windows only).'''
    if os.name != 'nt':
        return float('nan')

    class Counters(ctypes.Structure):
        _fields_ = [('cb', ctypes.c_uint32), ('PageFaultCount', ctypes.c_uint32),
                    ('PeakWorkingSetSize', ctypes.c_size_t),
                    ('WorkingSetSize', ctypes.c_size_t),
                    ('QuotaPeakPagedPoolUsage', ctypes.c_size_t),
                    ('QuotaPagedPoolUsage', ctypes.c_size_t),
                    ('QuotaPeakNonPagedPoolUsage', ctypes.c_size_t),
                    ('QuotaNonPagedPoolUsage', ctypes.c_size_t),
                    ('PagefileUsage', ctypes.c_size_t),
                    ('PeakPagefileUsage', ctypes.c_size_t)]
    c = Counters()
    c.cb = ctypes.sizeof(c)
    ctypes.windll.kernel32.GetCurrentProcess.restype = ctypes.c_void_p
    handle = ctypes.windll.kernel32.GetCurrentProcess()
    fn = getattr(ctypes.windll.kernel32, 'K32GetProcessMemoryInfo', None)
    if fn is None:
        fn = ctypes.windll.psapi.GetProcessMemoryInfo
    fn.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_uint32]
    fn(handle, ctypes.byref(c), c.cb)
    return float(c.PeakWorkingSetSize)


def build(ny, nx, table_dtype, seed):
    '''Build one float32 strip on the scipy backend, as `build_strips` does.'''
    previous = set_fft_backend('scipy')
    try:
        fac = ScreenFactory(ny, DX, L0_m=L0, nx=nx, dtype=np.float32,
                            table_dtype=table_dtype)
        return fac.make(R0, np.random.default_rng(seed))
    finally:
        set_fft_backend(previous)


def check_filter_identity():
    '''Measurement 1: the default filter equals the old meshgrid formula.'''
    print('1. the default (float64) filter against the old meshgrid formula')
    ok = True
    for n, nx in ((128, None), (128, 512)):
        fac = ScreenFactory(n, DX, L0_m=L0, nx=nx, dtype=np.float32)
        nxc = n if nx is None else nx
        dfy, dfx = 1.0 / (n * DX), 1.0 / (nxc * DX)
        df = dfy if dfy == dfx else np.sqrt(dfy * dfx)
        FX, FY = np.meshgrid((np.arange(nxc) - nxc // 2) * dfx,
                             (np.arange(n) - n // 2) * dfy)
        psd = _phase_psd_unit(np.hypot(FX, FY), L0, fac.l0_m)
        psd[n // 2, nxc // 2] = 0.0
        if nx is not None:
            psd[n // 2, :] = 0.0
        old = (np.sqrt(psd) * df).astype(np.float32)
        same = bool(np.array_equal(old, np.asarray(fac._filt)))
        ok = ok and same
        print(f"  {'PASS' if same else 'FAIL'}  shape {(n, nxc)!s:12s} "
              f"bit-identical {same}")
    return ok


def compare(ny, nx, seeds):
    '''Measurement 2: float32 tables against float64 tables, same seed.'''
    print(f'2. float32 tables against float64 tables, strip {(ny, nx)}, '
          f'{len(seeds)} seeds')
    ks = np.rint(R_LIST / DX).astype(int)
    ks = ks[ks < ny]
    theory = theory_dphi(ks * DX)
    max_d, rel = [], []
    d64, d32 = [], []
    for seed in seeds:
        s64 = build(ny, nx, None, seed)
        s32 = build(ny, nx, np.float32, seed)
        diff = s32.astype(np.float64) - s64.astype(np.float64)
        max_d.append(float(np.max(np.abs(diff))))
        rel.append(float(np.sqrt(np.mean(diff ** 2)) / np.std(s64)))
        d64.append(np.concatenate(dphi_axes(s64, ks)))
        d32.append(np.concatenate(dphi_axes(s32, ks)))
    d64 = np.mean(d64, axis=0)
    d32 = np.mean(d32, axis=0)
    ratio = d32 / d64
    print(f"  screen rms {float(np.std(s64)):.3f} rad; max |dphi| "
          f"{max(max_d):.3e} rad; rel rms {max(rel):.3e}")
    print(f"  {'r [m]':>8s} {'D64 y':>9s} {'D32 y':>9s} {'D64 x':>9s} "
          f"{'D32 x':>9s} {'theory':>9s} {'ratio y':>10s} {'ratio x':>10s}")
    m = len(ks)
    for i, k in enumerate(ks):
        print(f"  {k * DX:8.3f} {d64[i]:9.4f} {d32[i]:9.4f} {d64[m + i]:9.4f} "
              f"{d32[m + i]:9.4f} {theory[i]:9.4f} {ratio[i]:10.6f} "
              f"{ratio[m + i]:10.6f}")
    ok_phase = max(max_d) < MAX_DPHI_RAD
    ok_ratio = float(np.max(np.abs(ratio - 1.0))) < MAX_DPHI_RATIO
    print(f"  {'PASS' if ok_phase else 'FAIL'}  max |dphi| < {MAX_DPHI_RAD:g} rad")
    print(f"  {'PASS' if ok_ratio else 'FAIL'}  |D32/D64 - 1| < "
          f"{MAX_DPHI_RATIO:g} on both axes")
    return ok_phase and ok_ratio, {
        'ny': ny, 'nx': nx, 'seeds': list(seeds), 'max_dphi_rad': max(max_d),
        'rel_rms': max(rel), 'r_m': (ks * DX).tolist(),
        'D64_y': d64[:m].tolist(), 'D32_y': d32[:m].tolist(),
        'D64_x': d64[m:].tolist(), 'D32_x': d32[m:].tolist(),
        'theory': theory.tolist()}


def cost_child(ny, nx, table):
    '''The body of one cost measurement, in its own process.'''
    tracemalloc.start()
    t0 = time.perf_counter()
    strip = build(ny, nx, None if table == 'float64' else np.float32, 1)
    wall = time.perf_counter() - t0
    _, peak = tracemalloc.get_traced_memory()
    print(json.dumps({'table': table, 'ny': ny, 'nx': nx, 'wall_s': wall,
                      'tracemalloc_peak_gib': peak / 2 ** 30,
                      'peak_working_set_gib': peak_working_set_bytes() / 2 ** 30,
                      'rms': float(np.std(strip))}))


def cost(ny, nx, tables):
    '''Measurement 3: the build cost of each route, one fresh process each.'''
    print(f'3. the cost of one strip build, {(ny, nx)} = {ny * nx / 1e6:.0f} Mpx')
    rows = []
    for table in tables:
        out = subprocess.run(
            [sys.executable, '-m', 'validation.temporal_screens.table_precision',
             '--child', table, str(ny), str(nx)],
            capture_output=True, text=True, check=True)
        row = json.loads(out.stdout.strip().splitlines()[-1])
        rows.append(row)
        print(f"  {table:8s} wall {row['wall_s']:8.1f} s  tracemalloc peak "
              f"{row['tracemalloc_peak_gib']:6.2f} GiB  working-set peak "
              f"{row['peak_working_set_gib']:6.2f} GiB  rms {row['rms']:.3f}")
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--ny', type=int, default=512)
    parser.add_argument('--nx', type=int, default=4096)
    parser.add_argument('--seeds', type=int, default=4)
    parser.add_argument('--big', type=int, nargs=2, metavar=('NY', 'NX'),
                        default=(2048, 16384),
                        help='the strip shape of the cost measurement')
    parser.add_argument('--tables', nargs='+', default=['float64', 'float32'],
                        help='the routes of the cost measurement')
    parser.add_argument('--no-cost', dest='cost', action='store_false')
    parser.add_argument('--child', nargs=3, metavar=('TABLE', 'NY', 'NX'),
                        help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.child:
        cost_child(int(args.child[1]), int(args.child[2]), args.child[0])
        return
    ok1 = check_filter_identity()
    ok2, record = compare(args.ny, args.nx, list(range(1, args.seeds + 1)))
    record['filter_identity'] = ok1
    record['gate'] = ok1 and ok2
    if args.cost:
        record['cost'] = cost(args.big[0], args.big[1], args.tables)
    os.makedirs(DATA, exist_ok=True)
    path = os.path.join(DATA, f'table_precision_{args.big[0]}x{args.big[1]}.json')
    with open(path, 'w', encoding='utf-8') as handle:
        json.dump(record, handle, indent=1)
    print(f"{'GATE PASS' if record['gate'] else 'GATE FAIL'}; record {path}")
    if not record['gate']:
        sys.exit(1)


if __name__ == '__main__':
    main()
