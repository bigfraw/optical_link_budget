"""Compare the screen memory of the two fidelity-2 retro routes.

A retro trial needs the up leg and the return of one pulse, which share one
line and are 2R/c apart in time. Two routes give that:

- THE SHIFT ROUTE (built, `direction="retro"`): one oversize draw for each
  trial, and the up leg reads a window shifted by the wind drift
  |V(h)| * 2R/c at each screen (`olb.waveoptics.turbulence.run
  .retro_sensing_geometry`). The screen grows by the largest shift only.
- THE STRIP ROUTE (not built): a frozen-flow record of two frames at
  dt = 2R/c with `slew_rad_s=0.0` (`TemporalSpec`). Each layer is a strip with
  a seam pad of 2 L0 on the long axis, and each trial needs its own record
  (an independent atmosphere).

The script sizes both on the production grid and plan of the hero case (a
0.7 m ground terminal, `standard` preset, L0 = 25 m, Bufton Vg = 10 m/s) at two
altitudes and two elevations. It propagates nothing.

Run from the repository root:

    python -m validation.retro_bracket.retro_screen_sizes
"""

import shutil
import tempfile
import warnings
from dataclasses import replace

import numpy as np

from olb.geometry import CircularOrbit
from olb.waveoptics.turbulence import Campaign
from olb.waveoptics.turbulence.run import (_screen_draw_n,
                                           retro_sensing_geometry)
from olb.waveoptics.turbulence.temporal import TemporalSpec, strip_plan
from validation.anisoplanatism_screens.common import hero_uplink

L0_M = 25.0
VG_M_S = 10.0


def main():
    warnings.simplefilter("ignore")
    print(f"{'alt':>6} {'el':>4} {'2R/c':>8} {'max shift':>14} "
          f"{'shift route':>24} {'strip route (2 frames)':>26}")
    for alt in (500e3, 1000e3):
        for el in (30.0, 20.0):
            scn, _ = hero_uplink(el)
            scn = replace(scn, direction="retro", precompensation=None,
                          channel=replace(scn.channel, altitude_m=alt))
            geom = CircularOrbit(altitude_m=alt, elevation_deg=[el])
            root = tempfile.mkdtemp(prefix="olb_retro_sizes_")
            camp = Campaign(scn, geom, root, seed=1, preset="standard",
                            L0_m=L0_M, retro_wind_ground_m_s=VG_M_S)
            grid, plan = camp.grid, camp.plan
            dt = 2 * float(np.ravel(geom.slant_range_m)[0]) / 299792458.0
            shift = retro_sensing_geometry(plan, geom, VG_M_S).shift_m
            n, n_screens = grid.n, int(np.size(plan.z_m))
            n_draw = _screen_draw_n(n, shift.max(), grid.pixel_m)
            plain_mb = n_screens * n * n * 4 / 2 ** 20
            shift_mb = n_screens * n_draw * n_draw * 4 / 2 ** 20
            spec = TemporalSpec(dt_s=dt, n_frames=2, slew_rad_s=0.0,
                                strip_dir=root, wind_ground_m_s=VG_M_S)
            sp = strip_plan(plan, grid, spec, geom, L0_M)
            strip_px = sum(a * b for a, b in sp.shape)
            strip_mb = strip_px * 4 / 2 ** 20
            shutil.rmtree(root, ignore_errors=True)
            print(f"{alt / 1e3:4.0f}km {el:4.0f} {dt * 1e3:6.2f}ms "
                  f"{shift.max() * 100:5.1f} cm {shift.max() / grid.pixel_m:3.0f} px "
                  f"{n_draw:4d} px {shift_mb:4.0f} MB ({shift_mb / plain_mb:.2f}x) "
                  f"{max(b for _, b in sp.shape):6d} px {strip_mb:5.0f} MB "
                  f"({strip_mb / plain_mb:.1f}x)")
            print(f"{'':>30}shift by screen [cm]: "
                  + " ".join(f"{v * 100:.1f}" for v in shift))
    print(f"(plain stack: {n_screens} screens of {n} px, {plain_mb:.0f} MB; "
          f"the strip seam pad is 2 L0 = {2 * L0_M:.0f} m)")


if __name__ == "__main__":
    main()
