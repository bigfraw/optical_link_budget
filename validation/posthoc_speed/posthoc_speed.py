"""Measure the 2026-09-07 post-hoc read of a stored fidelity-2 campaign.

THE CHANGE. A stored trial holds the pixels of the patch disc only. The old
read-back scattered those pixels into the FULL grid, and every later step then
swept the zero padding, which is most of the grid. The new read works on the
SQUARE CROP that just holds the disc, and it builds the clip mask, the fibre
mode, the modal basis and the slope reconstructor ONE time for a call.

THE CROP RULE: PUPIL-plane quantities on the crop, FOCAL-plane quantities on
the padded grid. The collected power and the single-mode overlap are masked
sums over the pupil, so the crop gives them exactly. A multimode light bucket
and a camera focus the field, and the focal-plane pixel scale reads the grid
EXTENT, so those pad the crop back to the full grid.

THE SCRIPT. It stores a small campaign of the 30 deg hero downlink at 1024 px,
then it measures:

1. The AGREEMENT. `recollect`, `recouple` (SMF and MMF) and
   `recouple_compensated` (both sensing sources) on the crop against the full
   grid. The pupil quantities agree at the float rounding level, and the padded
   MMF route is bit-identical.
2. The TIME of one trial, full grid against crop, for `recouple` and for
   `recouple_compensated(source="slopes")`.
3. The BLOCK-PARALLEL speed-up of `Campaign.map_trials` with a process pool.

It writes `posthoc_speed.json` next to itself. Run it from the repository root
on an idle machine: the timing is the point.

    python -m validation.posthoc_speed.posthoc_speed
    python -m validation.posthoc_speed.posthoc_speed --n 512 --trials 8

Sources: Goodman, ISBN 978-0974707723 (the overlap integral); Ruilier and
Cassaing, DOI 10.1364/JOSAA.18.000143 (the fibre-mode match); Noll,
DOI 10.1364/JOSA.66.000207 (the modal basis).
"""
import argparse
import json
import os
import time
import warnings

import numpy as np

from olb.geometry import CircularOrbit
from olb.scenario import Channel, Site, SpaceScenario
from olb.terminal import MMF, SMF, TipTilt, Terminal, Transmitter
from olb.waveoptics.grid import GridSpec
from olb.waveoptics.turbulence import Campaign
from olb.waveoptics.turbulence.sampling import turbulent_grid

HERE = os.path.dirname(os.path.abspath(__file__))
SEED = 20260907
L0_M = 25.0


def hero_scenario():
    """Give the 30 deg hero downlink: a 0.7 m SMF ground terminal, 500 km.

    Returns:
        The pair (scenario, geometry).
    """
    ground = Terminal(aperture_m=0.7, obscuration_ratio=0.3,
                      wavelength_m=1550e-9, pointing_jitter_rad=1e-6,
                      detector=SMF(sensitivity_dbm=-40.0, focal_length_m=2.0,
                                   mode_field_radius_m=5.2e-6))
    space = Terminal(aperture_m=0.1, obscuration_ratio=0.0,
                     wavelength_m=1550e-9, pointing_jitter_rad=1e-6,
                     transmitter=Transmitter(waist_m=0.03, power_dbm=30.0,
                                             divergence_rad=None))
    scn = SpaceScenario(ground=ground, space=space,
                        channel=Channel(site=Site(), altitude_m=500e3),
                        direction="downlink")
    return scn, CircularOrbit(altitude_m=500e3, elevation_deg=30.0)


def timed(fn, repeat=1):
    """Give the best wall time of a callable, in seconds.

    Args:
        fn:     the callable. It takes no argument.
        repeat: the number of calls. The function keeps the shortest.

    Returns:
        The pair (seconds, the value of the last call).
    """
    best, out = None, None
    for _ in range(int(repeat)):
        t0 = time.perf_counter()
        out = fn()
        dt = time.perf_counter() - t0
        best = dt if best is None else min(best, dt)
    return best, out


def worst_rel(a, b):
    """Give the worst relative difference of two arrays."""
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    return float(np.abs(a / b - 1.0).max())


def main():
    """Store the campaign, check the agreement, and measure the time."""
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--n", type=int, default=1024,
                   help="the pixel count of one grid side")
    p.add_argument("--trials", type=int, default=32,
                   help="the number of stored trials")
    p.add_argument("--block-size", type=int, default=8,
                   help="the number of trials in one block")
    p.add_argument("--workers", type=int, default=4,
                   help="the pool size of the parallel measurement")
    args = p.parse_args()

    scn, geo = hero_scenario()
    grid, plan = turbulent_grid(scn, geo, preset="rapid", L0_m=L0_M)[:2]
    grid = GridSpec(size_m=grid.size_m, n=int(args.n), scaled=grid.scaled)
    root = os.path.join(HERE, "campaigns", f"hero_{args.n}px")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        camp = Campaign(scn, geo, root, seed=SEED, preset="rapid",
                        block_size=args.block_size, grid=grid, plan=plan,
                        L0_m=L0_M, store_screen_phase=True)
        print(f"campaign {camp.root_dir}")
        print(f"  grid {camp.grid.n} px, {camp.grid.size_m:.3f} m; "
              f"{camp.plan.z_m.size} screens")
        t0 = time.perf_counter()
        camp.run(args.trials, workers=args.workers)
        print(f"  {camp.n_stored} trials stored "
              f"({time.perf_counter() - t0:.1f} s)")

    crop = camp.patch.crop()
    ratio = (camp.patch.n / crop.side) ** 2
    print(f"  patch radius {camp.patch.radius_m * 100:.1f} cm, "
          f"{camp.patch.indices.size} pixels")
    print(f"  crop side {crop.side} px against the {camp.patch.n} px grid "
          f"({ratio:.1f}x fewer pixels)")

    ground = scn.ground
    mmf = MMF(sensitivity_dbm=-40.0, core_radius_m=25e-6, focal_length_m=2.0)
    stack = [TipTilt()]
    out = {"n": int(args.n), "trials": int(camp.n_stored),
           "block_size": int(args.block_size), "workers": int(args.workers),
           "crop_side": int(crop.side), "grid_n": int(camp.patch.n),
           "pixel_ratio": ratio}

    # ---- 1. the agreement of the crop and the full grid ----
    print("\nagreement, crop against full grid (worst relative difference):")
    checks = {}
    pw_crop = camp.recollect()
    pw_full = camp.recollect(compact=False)
    checks["recollect"] = worst_rel(pw_crop, pw_full)
    eta_crop = camp.recouple(ground.detector)
    eta_full = camp.recouple(ground.detector, compact=False)
    checks["recouple SMF"] = worst_rel(eta_crop, eta_full)
    mmf_crop = camp.recouple(mmf)
    mmf_full = camp.recouple(mmf, compact=False)
    checks["recouple MMF"] = worst_rel(mmf_crop, mmf_full)
    assert np.array_equal(mmf_crop, mmf_full), "the padded MMF must be exact"
    for src in ("screens", "slopes"):
        a = camp.recouple_compensated(stack, ground.detector, source=src)
        b = camp.recouple_compensated(stack, ground.detector, source=src,
                                      compact=False)
        checks[f"recouple_compensated {src}"] = worst_rel(a, b)
    for name, value in checks.items():
        print(f"  {name:32s} {value:9.2e}")
        assert value < 1e-6, (name, value)
    print("  recouple MMF is bit-identical (the crop pads back)")
    out["agreement"] = checks

    # ---- 2. the time of one trial ----
    print("\nseconds for one trial, one core:")
    rows = {}
    for name, call in (
            ("recouple SMF",
             lambda c: camp.recouple(ground.detector, compact=c)),
            ("recouple_compensated slopes",
             lambda c: camp.recouple_compensated(stack, ground.detector,
                                                 source="slopes", compact=c)),
            ("recouple_compensated screens",
             lambda c: camp.recouple_compensated(stack, ground.detector,
                                                 source="screens", compact=c)),
            ("recollect",
             lambda c: camp.recollect(compact=c))):
        t_full, _ = timed(lambda: call(False))
        t_crop, _ = timed(lambda: call(True))
        n = camp.n_stored
        rows[name] = {"full_s": t_full / n, "crop_s": t_crop / n,
                      "speedup": t_full / t_crop}
        print(f"  {name:30s} full {t_full / n:7.3f}  crop {t_crop / n:7.3f}  "
              f"{t_full / t_crop:5.1f}x")
    out["per_trial"] = rows

    # ---- 3. the block-parallel speed-up ----
    print(f"\nblock-parallel map_trials, {args.workers} workers:")
    par = {}
    for name, fn in (("recouple SMF",
                      lambda w: camp.recouple(ground.detector, workers=w)),
                     ("recouple_compensated slopes",
                      lambda w: camp.recouple_compensated(
                          stack, ground.detector, source="slopes",
                          workers=w))):
        t_one, a = timed(lambda: fn(None))
        t_pool, b = timed(lambda: fn(args.workers))
        assert np.array_equal(a, b), (name, a, b)
        par[name] = {"serial_s": t_one, "pool_s": t_pool,
                     "speedup": t_one / t_pool}
        print(f"  {name:30s} serial {t_one:7.2f} s  pool {t_pool:7.2f} s  "
              f"{t_one / t_pool:5.2f}x")
    print("  the pool gives the same numbers, value for value")
    out["parallel"] = par

    path = os.path.join(HERE, "posthoc_speed.json")
    with open(path, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nwrote {path}")


if __name__ == '__main__':
    main()
