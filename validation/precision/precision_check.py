"""A matched-seed single-against-double precision check of the fidelity-2 layer.

WHAT IT DOES. It runs the SAME turbulent trials two times, one time in double
precision (complex128 field, float64 screens) and one time in single precision
(complex64 field, float32 screens). The seed is the same, so trial k of the two
runs sees the SAME atmosphere. Then it compares, trial for trial:

- `collected_power`, the power inside the receive aperture,
- `smf_eta`, the single-mode-fibre coupling efficiency,
- the stored receive-plane field, as a relative rms over the patch.

It also gives the wall-time ratio of the two runs, on the threaded runner
(`workers=None`).

WHY. A `Campaign` on the desktop is memory-bandwidth bound. Half the bytes for
each element is the one change with real leverage there. So the question is
whether single precision changes the PHYSICS. This script answers it with
numbers.

THE CASE. It is the case of `validation/campaign_resources`: a 700 mm ground
aperture with a 30 percent central obscuration, an SMF detector, 1550 nm, a
500 km orbit at 30 deg, and the fixed outer scale L0 = 25 m. The preset is
`rapid`, so the check is short.

EXPECT about 1e-5 or better. A difference worse than 1e-3 means a kernel lost
precision. Look at the phase wrap of `olb.waveoptics.propagators.Forvard`
first: it keeps `Bus` in double precision on purpose.

Usage:

    python -m validation.precision.precision_check
    python -m validation.precision.precision_check --n-trials 100
    python -m validation.precision.precision_check --preset standard

It writes `precision_check.json` next to this file.
"""

import argparse
import json
import os
import sys
import time

import numpy as np

from olb.waveoptics.threader import Threader
from olb.waveoptics.turbulence.run import propagate_turbulent_scenario
from validation.campaign_resources.campaign_resources import (
    boost_process_priority, scenario_and_geometry)

HERE = os.path.dirname(os.path.abspath(__file__))
L0_M = 25.0
SEED = 20260905


def _rel(a, b):
    """Give the relative difference of two values, or NaN when b is zero."""
    if a is None or b is None or b == 0.0:
        return np.nan
    return abs(a / b - 1.0)


def run(n_trials, preset, elevation_deg, seed):
    """Run the double and the single case, and measure the difference.

    Args:
        n_trials:      the number of trials of each run.
        preset:        the name of a preset in sampling.PRESETS.
        elevation_deg: the elevation of the one line of sight, in deg.
        seed:          the integer base seed. The two runs share it.

    Returns:
        A dict of the results.
    """
    scn, geom = scenario_and_geometry(elevation_deg)
    patch_m = scn.rx_terminal.aperture_m / 2.0
    kw = dict(n_trials=n_trials, seed=seed, preset=preset, L0_m=L0_M,
              patch_radius_m=patch_m)

    out = {}
    for name, precision in (("double", "double"), ("single", "single")):
        t0 = time.perf_counter()
        out[name] = propagate_turbulent_scenario(
            scn, geom, threader=Threader(), precision=precision, **kw)
        out[name + "_s"] = time.perf_counter() - t0

    d, s = out["double"], out["single"]
    assert s.fields.dtype == np.complex64, s.fields.dtype

    rows = []
    for k, (a, b) in enumerate(zip(s.trials, d.trials)):
        # The stored patch is complex64 in BOTH runs, so this rms compares the
        # propagated field, not the store.
        fa, fb = s.fields[k].astype(np.complex128), d.fields[k].astype(
            np.complex128)
        peak = np.abs(fb).max()
        field_rms = float(np.sqrt(np.mean(np.abs(fa - fb) ** 2)) / peak)
        rows.append({
            "trial": k,
            "collected_power_double": float(b.collected_power),
            "collected_power_single": float(a.collected_power),
            "d_collected_power": float(_rel(a.collected_power,
                                            b.collected_power)),
            "smf_eta_double": None if b.smf_eta is None else float(b.smf_eta),
            "smf_eta_single": None if a.smf_eta is None else float(a.smf_eta),
            "d_smf_eta": float(_rel(a.smf_eta, b.smf_eta)),
            "d_field_rms": field_rms,
        })

    def worst(key):
        """Give the largest value of one column."""
        return float(np.nanmax([r[key] for r in rows]))

    return {
        "n_trials": int(n_trials), "preset": preset,
        "elevation_deg": float(elevation_deg), "seed": int(seed),
        "L0_m": L0_M,
        "grid_n": int(d.grid.n), "grid_size_m": float(d.grid.size_m),
        "n_screens": int(d.plan.z_m.size),
        "patch_pixels": int(d.patch.indices.size),
        "wall_double_s": out["double_s"], "wall_single_s": out["single_s"],
        "speed_ratio_double_over_single": out["double_s"] / out["single_s"],
        "max_d_collected_power": worst("d_collected_power"),
        "max_d_smf_eta": worst("d_smf_eta"),
        "max_d_field_rms": worst("d_field_rms"),
        "trials": rows,
    }


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--n-trials", type=int, default=20)
    ap.add_argument("--preset", default="rapid")
    ap.add_argument("--elevation", type=float, default=30.0)
    ap.add_argument("--seed", type=int, default=SEED)
    ap.add_argument("--out", default=None,
                    help="the JSON path (default: next to this file)")
    args = ap.parse_args(argv)

    boost_process_priority()
    res = run(args.n_trials, args.preset, args.elevation, args.seed)
    path = args.out or os.path.join(HERE, "precision_check.json")
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(res, fh, indent=2)

    print(f"case          : downlink {args.elevation:.0f} deg, "
          f"{args.preset} preset, L0 = {L0_M} m, seed {args.seed}")
    print(f"grid          : {res['grid_n']} px, {res['grid_size_m']:.3f} m, "
          f"{res['n_screens']} screens, patch {res['patch_pixels']} px")
    print(f"trials        : {res['n_trials']}")
    print()
    print("  k   collected power (single)   rel diff    smf eta (single)"
          "   rel diff    field rms")
    for r in res["trials"]:
        eta = "        none" if r["smf_eta_single"] is None else \
            f"{r['smf_eta_single']:12.6e}"
        print(f"{r['trial']:3d}   {r['collected_power_single']:22.9f}   "
              f"{r['d_collected_power']:8.2e}   {eta}   "
              f"{r['d_smf_eta']:8.2e}   {r['d_field_rms']:8.2e}")
    print()
    print(f"max rel diff, collected power : {res['max_d_collected_power']:.3e}")
    print(f"max rel diff, SMF eta         : {res['max_d_smf_eta']:.3e}")
    print(f"max rel rms, receive field    : {res['max_d_field_rms']:.3e}")
    print()
    print(f"wall time, double             : {res['wall_double_s']:.2f} s")
    print(f"wall time, single             : {res['wall_single_s']:.2f} s")
    print(f"speed ratio (double / single) : "
          f"{res['speed_ratio_double_over_single']:.2f}x")
    print(f"output                        : {path}")

    if res["max_d_collected_power"] > 1e-3 or res["max_d_smf_eta"] > 1e-3:
        print()
        print("WARNING: a difference worse than 1e-3. A kernel lost precision. "
              "Read the module docstring.")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
