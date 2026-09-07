"""The wrapped-gradient slopes, and the modal fit of those slopes.

THE PROBLEM. A stored receive-plane field holds the phase modulo 2 pi. A modal
fit of a wrapped phase map is wrong, because the wraps are large jumps that no
low-order Zernike holds. A 2D phase unwrap is slow and fragile.

THE ANSWER (method (a) of the plan). Take the phase DIFFERENCE between two
adjacent pixels, and rewrap that difference into (-pi, pi]. The difference of a
complex field is one product:

    sx[i, k] = angle( E[i, k+1] * conj(E[i, k]) )

The rewrapped difference is the TRUE local gradient wherever the true phase
step per pixel is below pi. Above pi the measurement aliases and the fit is
wrong. Use `max_abs_step` to test a grid. See the plan section 9, and Schmidt
(2010), DOI 10.1117/3.866274, Ch. 9, on the sampling of a phase screen.

THE FIT. `SlopeReconstructor` builds the slope-influence matrix: it rasters
each Noll mode, then it differences the raster with the SAME finite-difference
stencil. So the model and the measurement carry the same discretisation error,
and the two cancel. The module does NOT use the analytic Zernike derivative.

HOLD ENOUGH MODES. The slope fit and the direct phase fit are two different
least-squares metrics. An unfitted mode aliases into the fitted modes
differently in each metric. The self-check measures this: on a Kolmogorov
screen a 6-mode fit reads a tilt 6 % below the direct phase fit, and a 21-mode
fit agrees to 0.1 %. So do not use a very short mode set with the slope route.

PISTON IS UNOBSERVABLE. Piston (Noll j = 1) has zero slope everywhere, so its
column of the influence matrix is zero. The fit returns 0 for j = 1. A piston
is a constant phase, and it does not change a coupling efficiency.

The Zernike order is Noll's. Source: R. J. Noll, J. Opt. Soc. Am. 66, 207
(1976), DOI 10.1364/JOSA.66.000207.
"""

import numpy as np

from .zernike import mask_radius_px, zernike_j


def _wrap(d):
    """Rewrap an angle difference into (-pi, pi]."""
    return -(np.mod(-d + np.pi, 2.0 * np.pi) - np.pi)


def wrapped_gradient(field_or_phase, mask=None):
    """Give the wrapped phase difference between adjacent pixels.

    The function accepts a complex field or a real phase map. For a complex
    field it takes the angle of the product of a pixel and the conjugate of
    its neighbour, which is the cheapest form. For a real phase map it wraps
    the plain difference.

    THE VALIDITY LIMIT. The true phase step between two adjacent pixels must
    be below pi. Above pi the wrapped difference aliases. Test a grid with
    `max_abs_step`.

    Args:
        field_or_phase: an (N, N) complex field, or an (N, N) real phase map
                        in radians.
        mask:           an (N, N) bool array, or None. A pair with one pixel
                        outside the mask reads 0.0.

    Returns:
        The tuple (sx, sy). sx has the shape (N, N - 1) and it holds the
        difference along x. sy has the shape (N - 1, N) and it holds the
        difference along y. The unit is radians per pixel.
    """
    a = np.asarray(field_or_phase)
    if np.iscomplexobj(a):
        sx = np.angle(a[:, 1:] * np.conj(a[:, :-1]))
        sy = np.angle(a[1:, :] * np.conj(a[:-1, :]))
    else:
        sx = _wrap(a[:, 1:] - a[:, :-1])
        sy = _wrap(a[1:, :] - a[:-1, :])
    if mask is not None:
        mask = np.asarray(mask, dtype=bool)
        sx = np.where(mask[:, 1:] & mask[:, :-1], sx, 0.0)
        sy = np.where(mask[1:, :] & mask[:-1, :], sy, 0.0)
    return np.asarray(sx, dtype=np.float64), np.asarray(sy, dtype=np.float64)


def max_abs_step(sx, sy):
    """Give the largest absolute phase step of a slope pair, in radians.

    A value near pi means the grid is too coarse for the turbulence strength.
    The wrapped difference then aliases, and the modal fit is wrong. Flag such
    a grid. See the plan section 9.

    Args:
        sx: the x slopes of `wrapped_gradient`.
        sy: the y slopes of `wrapped_gradient`.

    Returns:
        The largest absolute value, in radians per pixel, as a float.
    """
    return float(max(np.abs(sx).max(), np.abs(sy).max()))


class SlopeReconstructor:
    """The least-squares fit of Noll coefficients to wrapped-gradient slopes.

    The object builds the slope-influence matrix one time, then it fits any
    number of measurements. Build it one time for one aperture mask, and keep
    it.

    Attributes:
        n_modes:  the number of Noll modes of the fit.
        n:        the number of pixels along one side.
        mask:     the (N, N) bool aperture mask.
        pair_x:   the (N, N - 1) bool mask of the used x pairs.
        pair_y:   the (N - 1, N) bool mask of the used y pairs.
        influence: the (n_pairs, n_modes) slope-influence matrix.
    """

    def __init__(self, n_modes, N, mask, radius_px=None):
        """Build the slope-influence matrix and its pseudo-inverse.

        Args:
            n_modes:   the number of Noll modes. Mode 1 is piston, and it is
                       unobservable from slopes.
            N:         the number of pixels along one side.
            mask:      an (N, N) bool array of the aperture.
            radius_px: the radius that maps to rho = 1, in pixels. The
                       default is the radius of the outermost in-mask pixel.
        """
        mask = np.asarray(mask, dtype=bool)
        if mask.shape != (int(N), int(N)):
            raise ValueError(
                f"SlopeReconstructor: the mask shape {mask.shape} does not "
                f"match the grid ({N}, {N}).")
        self.n_modes = int(n_modes)
        self.n = int(N)
        self.mask = mask
        self.radius_px = (mask_radius_px(mask) if radius_px is None
                          else float(radius_px))
        self.pair_x = mask[:, 1:] & mask[:, :-1]
        self.pair_y = mask[1:, :] & mask[:-1, :]
        n_pairs = int(self.pair_x.sum() + self.pair_y.sum())
        infl = np.empty((n_pairs, self.n_modes), dtype=np.float64)
        for k in range(self.n_modes):
            Z = zernike_j(self.n, k + 1, radius_px=self.radius_px)
            # Difference the raster with the stencil of `wrapped_gradient`.
            gx, gy = wrapped_gradient(Z, mask=mask)
            infl[:, k] = np.concatenate((gx[self.pair_x], gy[self.pair_y]))
        self.influence = infl
        # The discrete modes are only near-orthogonal, and the piston column
        # is zero. So use the pseudo-inverse.
        self._recon = np.linalg.pinv(infl)

    def measure(self, sx, sy):
        """Select the used pairs of a slope measurement.

        Args:
            sx: the x slopes of `wrapped_gradient`.
            sy: the y slopes of `wrapped_gradient`.

        Returns:
            A 1d array of the used slope values.
        """
        return np.concatenate((np.asarray(sx)[self.pair_x],
                               np.asarray(sy)[self.pair_y]))

    def estimate(self, sx, sy):
        """Fit the Noll coefficients to one slope measurement.

        Args:
            sx: the x slopes of `wrapped_gradient`.
            sy: the y slopes of `wrapped_gradient`.

        Returns:
            A (n_modes,) array of the Noll coefficients, in radians. The
            piston coefficient (j = 1) is 0.0, because piston has no slope.
        """
        coeffs = self._recon @ self.measure(sx, sy)
        coeffs[0] = 0.0
        return coeffs


if __name__ == '__main__':
    from .zernike import circle

    rng = np.random.default_rng(20260907)

    # 1. A smooth synthetic phase. The slope route must give the coefficients
    #    back, because the model and the measurement share the stencil.
    N = 128
    D = 100
    mask = circle(N, D)
    radius = mask_radius_px(mask)
    n_modes = 10
    truth = np.zeros(n_modes)
    truth[1:] = rng.normal(size=n_modes - 1) * 2.0     # a few radians
    phase = np.zeros((N, N))
    for k in range(n_modes):
        phase += truth[k] * zernike_j(N, k + 1, radius_px=radius)
    field = np.exp(1j * phase)

    rec = SlopeReconstructor(n_modes, N, mask)
    sx, sy = wrapped_gradient(field, mask=mask)
    got = rec.estimate(sx, sy)
    err = float(np.abs(got[1:] - truth[1:]).max())
    print(f"smooth phase, worst coefficient error {err:.2e} rad "
          f"(max step {max_abs_step(sx, sy):.3f} rad/px)")
    assert err < 1e-6, err

    # 2. A wrapped Kolmogorov screen. The slope route and the direct phase fit
    #    must agree on the tilts. This is the plan V3 check in miniature.
    from ..turbulence.screens import ScreenFactory
    from .modal import ApertureModes

    N = 256
    D = 128
    mask = circle(N, D)
    pixel_m = 1.0 / D                      # the aperture is 1 m
    r0_m = 1.0 / 5.0                       # D / r0 = 5
    factory = ScreenFactory(N, pixel_m, L0_m=np.inf)

    def tilt_pair(n_modes, n_seeds=20):
        """Give the direct tilts and the slope tilts of n_seeds screens."""
        modes = ApertureModes(n_modes, N, mask)
        rec = SlopeReconstructor(n_modes, N, mask)
        seed_rng = np.random.default_rng(7)
        direct, slope, step = [], [], 0.0
        for _ in range(n_seeds):
            screen = factory.make(r0_m, seed_rng)
            direct.append(modes.estimate(screen)[1:3])
            sx, sy = wrapped_gradient(np.exp(1j * screen), mask=mask)
            step = max(step, max_abs_step(sx, sy))
            slope.append(rec.estimate(sx, sy)[1:3])
        return (np.asarray(direct).ravel(), np.asarray(slope).ravel(), step)

    def compare(direct, slope):
        """Give the gain and the relative RMS difference of two tilt sets.

        A single ratio is noisy, because one tilt coefficient can be near
        zero. So measure the difference against the RMS tilt of the set.
        """
        rms = float(np.sqrt((direct ** 2).mean()))
        rel = float(np.sqrt(((slope - direct) ** 2).mean()) / rms)
        gain = float((slope @ direct) / (direct @ direct))
        return gain, rel, rms

    print(f"Kolmogorov D/r0 = 5, 40 tilt coefficients, N = {N}, D = {D} px")
    print(f"{'modes':>6} {'gain':>9} {'rel. RMS diff':>15}")
    for n_modes in (6, 21):
        direct, slope, step = tilt_pair(n_modes)
        gain, rel, rms = compare(direct, slope)
        print(f"{n_modes:>6} {gain:>9.5f} {rel * 100:>14.3f} %")
        if n_modes == 21:
            # THE FIT MUST HOLD ENOUGH MODES. A 6-mode fit reads a tilt gain
            # of 0.94, because the slope metric aliases the unfitted coma
            # (j = 7, 8) into the tilt differently from the phase metric. At
            # 21 modes the two routes agree to 0.1 %. This is the plan V3
            # check in miniature.
            assert abs(gain - 1.0) < 0.01, gain
            assert rel < 0.01, rel
    print(f"  (max phase step {step:.3f} rad/px, RMS tilt {rms:.3f} rad)")

    # 3. A deliberately coarse grid. A ramp of 4 rad per pixel aliases, and
    #    the step check reads it.
    N = 64
    mask = circle(N, 50)
    c = int(N / 2)
    Y, X = np.mgrid[:N, :N]
    ramp = 4.0 * (X - c)
    sx_c, sy_c = wrapped_gradient(np.exp(1j * ramp), mask=mask)
    step = max_abs_step(sx_c, sy_c)
    print(f"coarse ramp of 4.000 rad/px reads {step:.3f} rad/px "
          f"(aliased, the true step is above pi)")
    assert step > 2.0
    assert abs(step - 4.0) > 1.0            # the measurement is wrong

    print("self-check passed")
