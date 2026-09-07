"""The perfect-AO modal projector over one receive aperture.

`ApertureModes` holds the Noll Zernike basis over an aperture mask, and the
least-squares reconstructor of that basis. It gives four operations:

    estimate(phase)          -> the Noll coefficients of a sensed phase
    estimate_from_slopes()   -> the same coefficients from wrapped slopes
    reconstruct(coeffs)      -> the phase map of those coefficients
    apply(field, coeffs, s)  -> field * exp(s * 1j * reconstruct(coeffs))

THE SOURCE AND THE TARGET ARE SEPARATE. `estimate` takes an ARBITRARY sensing
phase. It does not assume that the phase comes from the field that `apply`
corrects. A receive-side correction is
`apply(field, estimate(source), -1)`. The sign +1 gives the conjugate
pre-distortion of an uplink pre-compensation. The pre-compensation itself is
not built here.

PERFECT AO. The fit is an ideal modal fit. There is no wavefront-sensor noise,
no finite subaperture, no aliasing, no servo lag and no branch point. So the
result is the UPPER BOUND of the benefit of a corrector of that mode count.

THE FIT IS OVER THE APERTURE. A real corrector sees one pupil wavefront, not
each turbulent layer. So the fit runs over the receive aperture mask, one
snapshot at a time.

The residual of a fit that removes the first J Noll modes follows Noll:

    sigma^2 = Delta_J * (D / r0)^(5/3)      [rad^2]

Source: R. J. Noll, "Zernike polynomials and atmospheric turbulence,"
J. Opt. Soc. Am. 66(3), 207-211 (1976), DOI 10.1364/JOSA.66.000207, Table IV,
printed p. 210. The self-check of this module measures that law.
"""

import numpy as np

from .slopes import SlopeReconstructor, wrapped_gradient
from .zernike import mask_radius_px, zernike_basis


def modes_from_stack(stack):
    """Give the number of Noll modes that a compensation stack removes.

    The best-correcting stage wins. A TipTilt stage removes the first 3 Noll
    modes. An AO(n) stage removes the first n Noll modes. Piston (j = 1) is
    part of that count. This is the count rule of
    `olb.turbulence.ao._noll_coefficient_and_modes`, so the analytic ladder
    and the wave-optics ladder count the same modes.

    The function reads the CLASS NAME of a stage. So this package does not
    import `olb.terminal`, and the one-way dependency holds.

    Args:
        stack: the ordered compensation stack of a Terminal. It may be empty.

    Returns:
        The number of removed Noll modes. An empty stack gives 0.

    Raises:
        ValueError: the stack holds an unknown stage.
    """
    n_modes = 0
    for stage in (stack or ()):
        name = type(stage).__name__
        if name == "TipTilt":
            n_modes = max(n_modes, 3)
        elif name == "AO":
            n_modes = max(n_modes, int(stage.n_modes))
        else:
            raise ValueError(
                f"modes_from_stack: unknown compensation stage {name!r}. Use "
                "TipTilt or AO.")
    return n_modes


class ApertureModes:
    """The Noll modal basis and the reconstructor of one receive aperture.

    Build the object one time for one (n_modes, N, mask). Keep it, and reuse
    it for every trial of a campaign. The build cost is the pseudo-inverse of
    a (n_pix, n_modes) matrix.

    Attributes:
        n_modes:   the number of Noll modes. Mode 1 is piston.
        n:         the number of pixels along one side.
        mask:      the (N, N) bool aperture mask.
        indices:   the flat indices of the in-mask pixels, in C order.
        basis:     the (n_pix, n_modes) mode matrix at those pixels.
        radius_px: the radius that maps to rho = 1, in pixels.
    """

    def __init__(self, n_modes, N, mask, radius_px=None):
        """Build the basis and the least-squares reconstructor.

        Args:
            n_modes:   the number of Noll modes to fit. Mode 1 is piston.
            N:         the number of pixels along one side.
            mask:      an (N, N) bool array of the aperture. Make it with
                       `olb.waveoptics.compensation.zernike.circle`.
            radius_px: the radius that maps to rho = 1, in pixels. The
                       default is the radius of the outermost in-mask pixel,
                       so the modes are normalised over the aperture.

        Raises:
            ValueError: n_modes is less than 1, or the mask shape is wrong.
        """
        n_modes = int(n_modes)
        if n_modes < 1:
            raise ValueError(
                f"ApertureModes: n_modes must be 1 or more, got {n_modes}.")
        mask = np.asarray(mask, dtype=bool)
        if mask.shape != (int(N), int(N)):
            raise ValueError(
                f"ApertureModes: the mask shape {mask.shape} does not match "
                f"the grid ({N}, {N}).")
        self.n_modes = n_modes
        self.n = int(N)
        self.mask = mask
        self.indices = np.flatnonzero(mask.ravel())
        self.radius_px = (mask_radius_px(mask) if radius_px is None
                          else float(radius_px))
        self.basis = zernike_basis(self.n, n_modes, mask,
                                   radius_px=self.radius_px)
        # The discrete modes are only near-orthogonal, and an obscured mask
        # breaks the orthogonality outright. So use the pseudo-inverse.
        self._recon = np.linalg.pinv(self.basis)
        self._slopes = None

    @property
    def n_pix(self):
        """The number of in-mask pixels."""
        return int(self.indices.size)

    def _flat(self, phase):
        """Give the in-mask values of a phase map or a flat vector."""
        a = np.asarray(phase)
        if a.ndim == 2:
            return a.ravel()[self.indices]
        if a.ndim == 1 and a.size == self.indices.size:
            return a
        raise ValueError(
            f"ApertureModes: the phase must be an ({self.n}, {self.n}) map or "
            f"a flat vector of {self.indices.size} in-mask values, got the "
            f"shape {a.shape}.")

    def estimate(self, phase):
        """Fit the Noll coefficients to a sensing phase.

        The phase is ARBITRARY. It is the summed screen phase of a space
        link, the unwrapped phase of a field, or any other sensed wavefront.
        The phase must NOT be wrapped. Use `estimate_from_slopes` for a
        wrapped field.

        Args:
            phase: an (N, N) real phase map, or the flat in-mask vector, in
                   radians.

        Returns:
            A (n_modes,) array of the Noll coefficients, in radians.
        """
        return self._recon @ self._flat(phase)

    def estimate_from_slopes(self, sx, sy):
        """Fit the Noll coefficients to wrapped-gradient slopes.

        The call builds the slope reconstructor one time, then it keeps it.
        The piston coefficient is 0.0, because piston has no slope. See
        `olb.waveoptics.compensation.slopes`.

        Args:
            sx: the x slopes of `wrapped_gradient`.
            sy: the y slopes of `wrapped_gradient`.

        Returns:
            A (n_modes,) array of the Noll coefficients, in radians.
        """
        if self._slopes is None:
            self._slopes = SlopeReconstructor(self.n_modes, self.n, self.mask,
                                              radius_px=self.radius_px)
        return self._slopes.estimate(sx, sy)

    def reconstruct(self, coeffs):
        """Give the phase map of a set of Noll coefficients.

        The map is zero outside the aperture mask.

        Args:
            coeffs: a (n_modes,) array of the Noll coefficients, in radians.

        Returns:
            An (N, N) float64 phase map, in radians.
        """
        coeffs = np.asarray(coeffs, dtype=np.float64)
        if coeffs.shape != (self.n_modes,):
            raise ValueError(
                f"ApertureModes.reconstruct: the coefficients must have the "
                f"shape ({self.n_modes},), got {coeffs.shape}.")
        out = np.zeros(self.n * self.n, dtype=np.float64)
        out[self.indices] = self.basis @ coeffs
        return out.reshape(self.n, self.n)

    def apply(self, field, coeffs, sign=-1):
        """Multiply a field by the phase screen of a set of coefficients.

        The result is `field * exp(sign * 1j * reconstruct(coeffs))`.

        Use sign = -1 to REMOVE the sensed phase. That is a receive-side
        correction. Use sign = +1 to ADD the conjugate phase. That is the
        pre-distortion of an uplink pre-compensation.

        The amplitude does not change. A real corrector cannot fix the
        scintillation.

        Args:
            field:  an (N, N) complex array.
            coeffs: a (n_modes,) array of the Noll coefficients, in radians.
            sign:   -1 to remove the phase, +1 to add it.

        Returns:
            An (N, N) complex array of the dtype of the input field.

        Raises:
            ValueError: the sign is not -1 or +1.
        """
        if sign not in (-1, 1, -1.0, 1.0):
            raise ValueError(
                f"ApertureModes.apply: the sign must be -1 or +1, got {sign}.")
        field = np.asarray(field)
        phase = self.reconstruct(coeffs)
        out = field * np.exp((sign * 1j) * phase)
        return out.astype(field.dtype, copy=False)

    def residual_variance(self, phase, coeffs):
        """Give the phase variance that a fit leaves over the aperture.

        The value is the variance of (phase - reconstruct(coeffs)) over the
        in-mask pixels. The variance removes the mean, so the piston is
        removed in every case. Noll's Table IV reports the same quantity.
        Source: Noll 1976, DOI 10.1364/JOSA.66.000207, Table IV, printed
        p. 210.

        Args:
            phase:  an (N, N) real phase map, or the flat in-mask vector.
            coeffs: a (n_modes,) array of the Noll coefficients.

        Returns:
            The residual phase variance, in rad^2, as a float.
        """
        res = self._flat(phase) - self._flat(self.reconstruct(coeffs))
        return float(np.var(res))


if __name__ == '__main__':
    import time

    from ..turbulence.screens import ScreenFactory
    from .zernike import circle

    # The Noll residual coefficients Delta_J of Table IV, printed p. 210.
    # Source: Noll 1976, DOI 10.1364/JOSA.66.000207.
    NOLL_DELTA = {1: 1.0299, 3: 0.134, 10: 0.0401, 21: 0.0208}

    N = 256
    D_px = 128
    D_m = 1.0
    pixel_m = D_m / D_px
    d_over_r0 = 6.0
    r0_m = D_m / d_over_r0
    n_seeds = 200

    mask = circle(N, D_px)
    projectors = {J: ApertureModes(J, N, mask) for J in NOLL_DELTA}

    # V2 in miniature: the residual variance of the removed-mode count must
    # follow Noll's Table IV. The screens are Kolmogorov (L0 = inf), because
    # Noll assumes the Kolmogorov spectrum.
    t0 = time.perf_counter()
    factory = ScreenFactory(N, pixel_m, L0_m=np.inf)
    rng = np.random.default_rng(20260907)
    totals = {J: 0.0 for J in NOLL_DELTA}
    for _ in range(n_seeds):
        screen = factory.make(r0_m, rng)
        for J, mo in projectors.items():
            totals[J] += mo.residual_variance(screen, mo.estimate(screen))
    build_s = time.perf_counter() - t0

    scale = d_over_r0 ** (5.0 / 3.0)
    print(f"V2 Noll residual check: N = {N}, D = {D_px} px, "
          f"D/r0 = {d_over_r0:g}, {n_seeds} seeds, L0 = inf")
    print(f"{'J':>4} {'measured [rad^2]':>18} {'Noll [rad^2]':>14} "
          f"{'ratio':>8}")
    ratios = {}
    for J in sorted(NOLL_DELTA):
        meas = totals[J] / n_seeds
        noll = NOLL_DELTA[J] * scale
        ratios[J] = meas / noll
        print(f"{J:>4} {meas:>18.4f} {noll:>14.4f} {ratios[J]:>8.3f}")
    print(f"({n_seeds} screens and 4 fits took {build_s:.1f} s)")

    # THE J = 1 GAP IS THE OUTER SCALE, NOT THE FIT. The screens ask for
    # L0 = inf, but the three subharmonic levels reach 27 times the grid side
    # only. So the screens hold a finite outer scale, and they miss part of
    # the tilt. Noll's Delta_1 is a pure Kolmogorov value, so the measured
    # J = 1 residual reads LOW. See CLAUDE.md, backlog 2-P5. The J = 3 and
    # higher values remove the tilt, so they do not see that deficit and they
    # agree with the book.
    #
    # The test below shows the cause: it holds the aperture and the r0, and it
    # makes the GRID wider. A wider grid gives a larger effective outer scale,
    # so the J = 1 ratio must move toward 1.0.
    wide_ratios = []
    for grid_n in (256, 512):
        wf = ScreenFactory(grid_n, pixel_m, L0_m=np.inf)
        wm = circle(grid_n, D_px)
        wp = ApertureModes(1, grid_n, wm)
        wrng = np.random.default_rng(4242)
        tot = 0.0
        for _ in range(50):
            s = wf.make(r0_m, wrng)
            tot += wp.residual_variance(s, wp.estimate(s))
        wide_ratios.append(tot / 50 / (NOLL_DELTA[1] * scale))
        print(f"J = 1 ratio at the grid {grid_n} px "
              f"(grid side {grid_n * pixel_m:.1f} m): {wide_ratios[-1]:.3f}")
    assert wide_ratios[1] > wide_ratios[0], wide_ratios

    assert 0.60 <= ratios[1] <= 1.10, ratios[1]
    for J in (3, 10, 21):
        assert 0.85 <= ratios[J] <= 1.15, (J, ratios[J])

    # Idempotence: a second fit of the corrected phase gives zero.
    mo = projectors[10]
    screen = factory.make(r0_m, rng)
    c1 = mo.estimate(screen)
    left = screen - mo.reconstruct(c1)
    c2 = mo.estimate(left)
    print(f"idempotence: worst refit coefficient {np.abs(c2).max():.2e} rad "
          f"(first fit {np.abs(c1).max():.3f} rad)")
    assert np.abs(c2).max() < 1e-9 * max(1.0, float(np.abs(c1).max()) * 1e9)
    assert np.abs(c2).max() < 1e-8

    # apply(-1) then apply(+1) gives the field back.
    field = (np.exp(1j * screen) * mask).astype(np.complex64)
    back = mo.apply(mo.apply(field, c1, -1), c1, +1)
    assert back.dtype == np.complex64
    err = float(np.abs(back - field).max())
    print(f"apply(-1) then apply(+1): worst field error {err:.2e}")
    assert err < 1e-5

    # The correction removes the fitted phase from the field.
    corrected = mo.apply(field, c1, -1)
    inside = mask & (np.abs(field) > 0)
    # The comparison runs on the complex phasor, because the removed phase is
    # several radians and an angle difference wraps.
    d = corrected[inside] * np.conj(field[inside])
    want = np.exp(-1j * mo.reconstruct(c1)[inside])
    assert np.abs(d - want).max() < 1e-4

    # The stack map matches the count rule of olb.turbulence.ao.
    class TipTilt:
        pass

    class AO:
        def __init__(self, n_modes):
            self.n_modes = n_modes

    assert modes_from_stack([]) == 0
    assert modes_from_stack([TipTilt()]) == 3
    assert modes_from_stack([TipTilt(), AO(20)]) == 20
    assert modes_from_stack([AO(2), TipTilt()]) == 3
    try:
        modes_from_stack([object()])
    except ValueError:
        pass
    else:
        raise AssertionError("an unknown stage must raise")

    # The build cost of a realistic receive aperture.
    t0 = time.perf_counter()
    big = ApertureModes(60, 256, circle(256, 140))
    print(f"build ApertureModes(60, N = 256, D = 140 px): "
          f"{time.perf_counter() - t0:.2f} s, {big.n_pix} pixels")

    print("self-check passed")
