"""The gradient (centroid) tilt of a field, from the far-field centroid.

THE IDEA. The first moment of the far-field intensity is the intensity-weighted
mean phase gradient over the aperture. This mean gradient is the GRADIENT TILT,
or G-tilt: the tilt that a centroid tracker measures, and the tilt that centres
the focused spot. Source: G. A. Tyler, "Bandwidth considerations for tracking
through turbulence," J. Opt. Soc. Am. A 11(1), 358-367 (1994),
DOI 10.1364/JOSAA.11.000358 (the Z-tilt / G-tilt split, and G-tilt as the
aperture average of the gradient).

WHY NOT THE SLOPES. The wrapped-gradient slope route
(`olb.waveoptics.compensation.slopes`) differences ADJACENT pixels, so it
aliases wherever the LOCAL phase step passes pi, even when the OVERALL tilt is
small. Strong high-order structure and branch points break the slope fit that
way. The far-field centroid integrates the whole aperture, so a local strong
gradient scatters the far-field spot but does not move its centroid. So G-tilt
reads the overall tilt where the slope fit aliases. The one shared limit is
Nyquist: an overall tilt above pi per pixel wraps the spot in the FFT window,
exactly as it aliases a slope. `far_field_tilt` reports that margin so a caller
can flag it.

THE MEASUREMENT. Take the FFT of the aperture-masked field. The peak (the
centroid) of the far-field intensity sits at the spatial frequency of the mean
tilt (the Fourier shift theorem: a linear phase in the near field is a shift in
the far field). Its offset from the zero-frequency pixel is the mean gradient:

    g_x = 2 pi * (q_centroid - c) / N          [rad per pixel, along columns]
    g_y = 2 pi * (p_centroid - c) / N          [rad per pixel, along rows]

with c = int(N / 2), the zero-frequency pixel after `fftshift`, and (p, q) the
row and column of the intensity centroid. The hard aperture mask adds a
symmetric diffraction ring, which does not move the centroid.

THE PIXEL CONVENTION IS THE CONVENTION OF THE LAYER. The axis is the pixel
(int(N/2), int(N/2)), the convention of `olb.waveoptics.compensation.zernike`
and of `olb.waveoptics.sources.CircAperture`. `fftshift` puts the
zero-frequency pixel there too, so the tilt and the Zernike modes share one
origin.
"""

import numpy as np


def far_field_tilt(field, mask=None):
    """Give the G-tilt of a field, in radians per pixel.

    The value is the intensity-weighted mean phase gradient over the aperture,
    from the first moment of the far-field intensity. See the module docstring
    and Tyler 1994, DOI 10.1364/JOSAA.11.000358.

    Args:
        field: an (N, N) complex array. A real phase map is NOT accepted; the
               method reads the complex field, so it never unwraps a phase.
        mask:  an (N, N) bool array of the aperture, or None. The field is set
               to zero outside the mask before the transform, so the tilt is
               the aperture average. None uses the whole grid.

    Returns:
        The tuple (g_x, g_y, margin). g_x is the mean gradient along the
        columns (the x, or Noll j = 2, direction) and g_y is the mean gradient
        along the rows (the y, or Noll j = 3, direction), both in radians per
        pixel. `margin` is the centroid distance from the zero-frequency pixel
        as a fraction of the Nyquist half-grid: it approaches 1.0 as the tilt
        approaches pi per pixel and the far-field spot nears the FFT edge. A
        caller flags a value near 1.0, the G-tilt twin of an aliased slope.
    """
    E = np.asarray(field)
    if not np.iscomplexobj(E):
        raise ValueError(
            "far_field_tilt: the input must be a complex field, not a phase "
            "map. The method reads the field, so it never unwraps a phase.")
    if mask is not None:
        E = E * np.asarray(mask, dtype=bool)
    N = E.shape[0]
    c = int(N / 2)
    spec = np.abs(np.fft.fftshift(np.fft.fft2(E))) ** 2
    total = float(spec.sum())
    if total <= 0.0:
        return 0.0, 0.0, 0.0
    offs = np.arange(N, dtype=np.float64) - c
    q = float(spec.sum(axis=0) @ offs) / total       # column centroid (x)
    p = float(spec.sum(axis=1) @ offs) / total       # row centroid (y)
    g_x = 2.0 * np.pi * q / N
    g_y = 2.0 * np.pi * p / N
    margin = float(np.hypot(p, q)) / (0.5 * N)
    return g_x, g_y, margin


if __name__ == '__main__':
    from .zernike import circle, mask_radius_px, zernike_j
    from .slopes import SlopeReconstructor, max_abs_step, wrapped_gradient

    rng = np.random.default_rng(20260908)
    N = 256
    D = 200
    mask = circle(N, D)
    R = mask_radius_px(mask)
    c = int(N / 2)
    Y, X = np.mgrid[:N, :N]

    # 1. A PURE TILT ROUND-TRIP. A ramp of a known gradient must read back with
    #    that gradient, and the derived Noll coefficient must remove it. This
    #    catches any sign or normalisation error in the FFT convention.
    gx0, gy0 = 0.37, -0.21                            # rad per pixel
    ramp = gx0 * (X - c) + gy0 * (Y - c)
    field = np.exp(1j * ramp)
    gx, gy, margin = far_field_tilt(field, mask)
    print(f"pure tilt: asked g = ({gx0:+.3f}, {gy0:+.3f}) rad/px, "
          f"read ({gx:+.3f}, {gy:+.3f}), margin {margin:.3f}")
    assert abs(gx - gx0) < 2e-3 and abs(gy - gy0) < 2e-3, (gx, gy)

    # The Noll map: Z2 = 2 (X - c) / R, so d Z2 / d col = 2 / R per pixel, and
    # a2 = g_x R / 2. Rebuild the ramp from the coefficients and check the
    # residual gradient is near zero over the aperture.
    a2, a3 = gx * R / 2.0, gy * R / 2.0
    rebuilt = a2 * zernike_j(N, 2, radius_px=R) + a3 * zernike_j(N, 3, radius_px=R)
    corrected = field * np.exp(-1j * rebuilt)
    gxr, gyr, _ = far_field_tilt(corrected, mask)
    print(f"           after removal g = ({gxr:+.4f}, {gyr:+.4f}) rad/px")
    assert abs(gxr) < 5e-3 and abs(gyr) < 5e-3, (gxr, gyr)

    # 2. THE GROUND TRUTH. For a unit-amplitude field the intensity is flat, so
    #    the G-tilt is exactly the plain spatial mean of the phase gradient over
    #    the aperture. `far_field_tilt` must reproduce that mean on a REAL
    #    Kolmogorov screen, whatever the higher-order content. The screen phase
    #    is an unwrapped array, so the mean gradient is the reference, computed
    #    over an eroded mask to keep the edge out.
    from ..turbulence.screens import ScreenFactory

    def true_mean_gradient(phase, m):
        """The plain mean of the phase gradient over the mask, rad per pixel."""
        gy_map, gx_map = np.gradient(phase)
        inner = m & np.roll(m, 1, 0) & np.roll(m, -1, 0) & \
            np.roll(m, 1, 1) & np.roll(m, -1, 1)
        return float(gx_map[inner].mean()), float(gy_map[inner].mean())

    factory = ScreenFactory(N, 1.0 / D, L0_m=np.inf)     # the aperture is 1 m
    seed_rng = np.random.default_rng(11)
    print(f"{'D/r0':>5} {'step[rad/px]':>13} "
          f"{'G-tilt gx':>11} {'true gx':>9} {'slope gx':>9}")
    for d_over_r0 in (5.0, 12.0):
        rec = SlopeReconstructor(21, N, mask)
        worst_gtilt, worst_slope = 0.0, 0.0
        for _ in range(15):
            screen = factory.make(1.0 / d_over_r0, seed_rng)
            field = np.exp(1j * screen)
            tx, ty = true_mean_gradient(screen, mask)
            gx, gy, _ = far_field_tilt(field, mask)
            sx, sy = wrapped_gradient(field, mask=mask)
            step = max_abs_step(sx, sy)
            slope_gx = 2.0 * rec.estimate(sx, sy)[1] / R
            worst_gtilt = max(worst_gtilt, abs(gx - tx), abs(gy - ty))
            worst_slope = max(worst_slope, abs(slope_gx - tx))
        print(f"{d_over_r0:>5.0f} {step:>13.2f} {gx:>11.4f} "
              f"{tx:>9.4f} {slope_gx:>9.4f}")
        # G-tilt matches the true mean gradient at EVERY strength. The masked
        # FFT adds a symmetric diffraction ring only, so the centroid holds.
        assert worst_gtilt < 0.03, (d_over_r0, worst_gtilt)
    # The strong screen aliases the slope step (near pi), and there the slope
    # tilt drifts from the true mean gradient by more than G-tilt does. This is
    # the mechanism; the comparison script measures it on the campaign fields.
    print(f"  (strong screen: worst G-tilt error {worst_gtilt:.4f}, "
          f"worst slope error {worst_slope:.4f} rad/px)")
    assert worst_slope > worst_gtilt, (worst_slope, worst_gtilt)

    # 3. A ZERO FIELD gives a zero tilt, not a divide error.
    z = far_field_tilt(np.zeros((N, N), dtype=complex), mask)
    assert z == (0.0, 0.0, 0.0), z

    # 4. A phase map (real input) is refused.
    try:
        far_field_tilt(ramp, mask)
    except ValueError:
        pass
    else:
        raise AssertionError("a real phase map must raise")

    print("self-check passed")
