"""Step 1: the vacuum launch of the 10 km MMF link on flat grids and on the
Schmidt two-pitch chain.

THE QUESTION. The turbulent sizer clamps the 10 km terrestrial grid at
`n_max`, so the 4.54 mm launch waist gets about one pixel. Does the Schmidt
two-pitch chain (a fine pitch at the launch, a coarse pitch at the receiver)
hold the received field at a pixel count where the flat grid does not?

THE METHOD. VACUUM ONLY (no screens), so the truth is the analytic Gaussian
beam: the 25.4 mm launch aperture clips the 4.54 mm waist at 2.8 w0, which
removes exp(-15.6) of the power, so the clipped launch IS the Gaussian to
1e-7. The received field at z = L is

    E(r) = (w0 / w) exp(-r^2 / w^2) exp(-i k r^2 / (2R)) x (a constant phase)

(Siegman, Lasers (1986), Ch. 17; the same law as `olb.beam.gaussz`; Andrews
and Phillips, DOI 10.1117/3.626196, Ch. 4, Eqs. (7) and (8), printed p. 87).

- FLAT: the production loop, `split_step` with zero screens at the planes of
  the production plan, the production boundary, on the grid side of
  `turbulent_grid` at n_max = 2048 / 4096 / 8192.
- TWO-PITCH: `olb.waveoptics.schmidt.fresnel.partial_propagations` with the
  book absorber, the pitches and N from Schmidt Ch. 9, Eqs. (9.84) to (9.89),
  printed pp. 173 and 174 (DOI 10.1117/3.866274), and at least the
  production plane count.

THE METRICS, inside the stored-patch disc D2 = 38.1 mm (1.5 x the aperture
radius) and the 25.4 mm bucket:
- the field RMS error against the analytic field, relative, after the
  removal of one global phase;
- the bucket power against the analytic field summed over the SAME pixels
  (so the pixelised disc area cancels), in dB.

Run it from the repository root:

    python -m validation.terrestrial_sampling.vacuum_check
"""

import math
import time

import numpy as np

from olb.waveoptics.schmidt.fresnel import (partial_propagations,
                                            super_gaussian_absorber)
from olb.waveoptics.schmidt.turbulence import (blurred_extent,
                                               constraint1_pitch_max,
                                               constraint2_n_min,
                                               max_partial_step)
from olb.waveoptics.turbulence.sampling import PRESETS, turbulent_grid
from olb.waveoptics.turbulence.splitstep import (split_step,
                                                 super_gaussian_boundary)
from validation.terrestrial_mmf_tiptilt.run import (APERTURE_M, LAM, PATH_M,
                                                    WAIST_M, build_scenario)

import dataclasses

K = 2 * math.pi / LAM
ZR = math.pi * WAIST_M ** 2 / LAM
W_L = WAIST_M * math.sqrt(1 + (PATH_M / ZR) ** 2)
R_L = PATH_M * (1 + (ZR / PATH_M) ** 2)
D2 = 1.5 * APERTURE_M
CN2 = 1e-14                        # the strongest cell sets the side


def coords(n, dx):
    x = (np.arange(n) - n // 2) * dx
    return np.meshgrid(x, x)


def launch(n, dx, dtype=np.complex128):
    """The unit-power Gaussian clipped by the launch aperture."""
    X, Y = coords(n, dx)
    r2 = X ** 2 + Y ** 2
    U = np.sqrt(2 / math.pi) / WAIST_M * np.exp(-r2 / WAIST_M ** 2)
    U = U * (r2 <= (APERTURE_M / 2) ** 2)
    return U.astype(dtype)


def score(U, dx):
    """The field error and the bucket error of one received field."""
    X, Y = coords(U.shape[0], dx)
    r2 = X ** 2 + Y ** 2
    ref = (math.sqrt(2 / math.pi) / W_L * np.exp(-r2 / W_L ** 2)
           * np.exp(-1j * K * r2 / (2 * R_L)))
    disc = r2 <= (D2 / 2) ** 2
    best = None
    for E in (U, np.conj(U)):     # the propagators differ in sign convention
        e, r = E[disc], ref[disc]
        e = e * np.exp(-1j * np.angle(np.vdot(r, e)))
        rms = float(np.linalg.norm(e - r) / np.linalg.norm(r))
        best = rms if best is None else min(best, rms)
    bucket = r2 <= (APERTURE_M / 2) ** 2
    p = float((np.abs(U[bucket]) ** 2).sum() * dx ** 2)
    # The analytic field on the SAME pixels, so the pixelised disc area
    # cancels and only the field error is left.
    p_true = float((np.abs(ref[bucket]) ** 2).sum() * dx ** 2)
    return best, 10 * math.log10(p / p_true)


def flat(n_max):
    scn, geo = build_scenario(CN2)
    preset = dataclasses.replace(PRESETS["standard"], n_max=n_max)
    grid, plan = turbulent_grid(scn, geo, preset=preset, L0_m=25.0)[:2]
    from olb.waveoptics.field import Begin
    F = Begin(grid.size_m, LAM, grid.n, dtype=np.complex64)
    F.field = launch(grid.n, grid.pixel_m, np.complex64)
    zero = np.zeros((grid.n, grid.n), np.float32)
    mask = super_gaussian_boundary(grid.n, preset.boundary_width_frac)
    out = split_step(F, plan.z_m, (zero for _ in plan.z_m), PATH_M,
                     boundary=mask)
    return grid.n, grid.pixel_m, grid.pixel_m, out.field, plan.z_m.size


def two_pitch(n, dx1, dxn, n_planes):
    step = max_partial_step(dx1, dxn, n, LAM)
    n_planes = max(n_planes, math.ceil(PATH_M / step) + 1)
    z = np.linspace(0.0, PATH_M, n_planes)
    U = partial_propagations(launch(n, dx1), LAM, dx1, dxn, z,
                             absorber=super_gaussian_absorber(n))
    return U, n_planes


def schmidt_n(dx1, dxn, c):
    r0s = (0.423 * K ** 2 * CN2 * PATH_M * 3 / 8) ** -0.6
    d1 = blurred_extent(APERTURE_M, LAM, PATH_M, r0s, c)
    d2 = blurred_extent(D2, LAM, PATH_M, r0s, c)
    ok = dxn <= constraint1_pitch_max(dx1, d1, d2, LAM, PATH_M)
    n = constraint2_n_min(dx1, dxn, d1, d2, LAM, PATH_M)
    return ok, 2 ** math.ceil(math.log2(n))


def main():
    p_true = 1 - math.exp(-2 * (APERTURE_M / 2) ** 2 / W_L ** 2)
    print(f"w(L) = {W_L:.4f} m, R(L) = {R_L:.0f} m, "
          f"analytic bucket = {10 * math.log10(p_true):.3f} dB")
    print(f"{'route':34s} {'N':>5s} {'dx1 mm':>7s} {'dxn mm':>7s} "
          f"{'px/w0':>6s} {'planes':>6s} {'field RMS':>10s} {'bucket dB':>10s} "
          f"{'s':>6s}")
    rows = []
    for n_max in (2048, 4096, 8192):
        t = time.perf_counter()
        n, dx1, dxn, U, nscreen = flat(n_max)
        rms, db = score(U, dxn)
        rows.append((f"flat production n_max={n_max}", n, dx1, dxn,
                     nscreen + 1, rms, db, time.perf_counter() - t))
    for dx1, dxn, c in ((2e-3, 3e-3, 4), (2e-3, 4e-3, 2), (1.5e-3, 4e-3, 4),
                        (1e-3, 3e-3, 4), (1e-3, 6e-3, 2)):
        ok, n = schmidt_n(dx1, dxn, c)
        t = time.perf_counter()
        U, planes = two_pitch(n, dx1, dxn, 36)
        rms, db = score(U, dxn)
        rows.append((f"two-pitch c={c}{'' if ok else ' (C1 FAIL)'}", n, dx1,
                     dxn, planes, rms, db, time.perf_counter() - t))
    for name, n, dx1, dxn, planes, rms, db, s in rows:
        print(f"{name:34s} {n:5d} {dx1 * 1e3:7.2f} {dxn * 1e3:7.2f} "
              f"{WAIST_M / dx1:6.2f} {planes:6d} {rms:10.2e} {db:+10.4f} "
              f"{s:6.1f}")


if __name__ == "__main__":
    main()
