"""Smoke run of the 10 km terrestrial MMF campaigns: 10 trials for each Cn2.

The same scenario, seed and preset as run.py, but each trial stores a WIDE
field patch (radius 1.5 m, wider than the 1.09 m vacuum beam radius), so the
received beam, its speckle and its tilt can be looked at post hoc
(analyse_smoke.py). The roots are separate from the production roots, because
the patch radius and the block size enter the manifest.

    python -m validation.terrestrial_mmf_tiptilt.smoke --fft-backend cupy
"""

import argparse
import os

from olb.waveoptics.priority import boost_process_priority
from olb.waveoptics.turbulence import Campaign
from validation.terrestrial_mmf_tiptilt.run import (CN2, HERE, SEED,
                                                    build_scenario,
                                                    refined_grid)

PATCH_M = 1.5


def smoke_campaign(cn2, n, fft_backend="numpy", n_px=None):
    """n_px=None takes the sized grid; an int keeps its side and plan and
    refines the pixel only."""
    scn, geo = build_scenario(cn2)
    root = os.path.join(HERE, "campaigns", f"smoke_cn2{cn2:.0e}"
                        + (f"_n{n_px}" if n_px else ""))
    grid, plan = (refined_grid(scn, geo, "standard", n_px) if n_px
                  else (None, None))
    return Campaign(scn, geo, root, seed=SEED, preset="standard",
                    block_size=n, patch_radius_m=PATCH_M, grid=grid,
                    plan=plan, fft_backend=fft_backend)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--trials", type=int, default=10)
    ap.add_argument("--fft-backend", default="numpy")
    ap.add_argument("--n", type=int, default=None,
                    help="pixel count at the sized side (None: sized grid)")
    args = ap.parse_args()
    boost_process_priority()
    for cn2 in CN2:
        camp = smoke_campaign(cn2, args.trials, args.fft_backend, args.n)
        print(f"Cn2 {cn2:.0e}: grid {camp.grid.n} px, "
              f"dx {1e3 * camp.grid.size_m / camp.grid.n:.2f} mm, "
              f"{camp.plan.z_m.size} screens", flush=True)
        camp.run(args.trials, progress=True)


if __name__ == "__main__":
    main()
