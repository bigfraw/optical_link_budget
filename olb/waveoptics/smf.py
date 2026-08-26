"""The single-mode-fibre pupil mode, and the coupling efficiency.

A single-mode fibre accepts one field shape only. The back-propagated fibre
mode is that shape in the pupil plane. It is a fundamental Gaussian. The
coupling efficiency is the normalised overlap of the received field with that
mode.

The module transcribes two helpers from the shared kernel repository
(my_analysis_modules: lightpipes_atmospherics.smf, coupling_efficiency, and
the overlap and power kernels of general_atmospherics). The transcription keeps
the package self-contained: it imports numpy and the local field module only.

Sources:
- Ruilier, A study of degraded light coupling into single-mode fibers,
  DOI 10.1117/12.317094. The mode radius that gives the largest overlap with a
  uniformly illuminated circular pupil, and the 0.8145 maximum.
- Goodman, Introduction to Fourier Optics, ISBN 978-0974707723. The overlap
  integral of two scalar fields.
- LightPipes manual, https://opticspy.github.io/lightpipes/manual.html. The
  lineage of the source helpers.
"""

import numpy as np

from .field import Begin, Intensity, SubIntensity
from .sources import GaussBeam

# The best ratio of the pupil diameter to the fibre-mode radius. The source
# kernel writes the same number as w0 = 2.24 / d and then uses 1 / w0 as the
# radius. Ruilier, DOI 10.1117/12.317094: the overlap is largest when the
# coupling parameter is a = pi*(D/2)*w_m/(lambda*f) = 1.12, and the pupil-plane
# mode radius is then D / 2.24.
MODE_RADIUS_RATIO = 2.24


def smf_mode(grid_size_m, wavelength_m, n, aperture_m):
    """Make the back-propagated fibre mode of a single-mode fibre.

    The mode is a fundamental Gaussian of the radius aperture_m / 2.24. See
    Ruilier, DOI 10.1117/12.317094. The function scales the mode so that the
    sum of the intensity over the grid is 1.0. So the overlap needs no other
    normalisation.

    Args:
        grid_size_m:   the physical side of the square grid, in m.
        wavelength_m:  the wavelength, in m.
        n:             the number of pixels along one side.
        aperture_m:    the pupil DIAMETER, in m.

    Returns:
        A Field that holds the mode.
    """
    w_mode = aperture_m / MODE_RADIUS_RATIO
    F = GaussBeam(Begin(grid_size_m, wavelength_m, n), w_mode)
    I = Intensity(F)
    return SubIntensity(F, I / I.sum())


def coupling_efficiency(field, aperture_m, mask=None):
    """Calculate the power fraction that couples into a single-mode fibre.

    eta = |sum(E * conj(M))|^2 / sum(|E|^2), with M the fibre mode of
    smf_mode(). The mode carries sum(|M|^2) = 1, so eta stays between 0 and 1.
    See Goodman, ISBN 978-0974707723 (the overlap integral), and Ruilier,
    DOI 10.1117/12.317094 (the fibre-mode match).

    The source kernel computes the overlap with an FFT correlation and reads
    the zero-shift sample. Parseval makes that sample equal to the plain inner
    product, so this transcription uses the inner product.

    Args:
        field:      the received Field. The grid, the wavelength and the pixel
                    count of the mode come from this field.
        aperture_m: the pupil DIAMETER, in m.
        mask:       an optional N x N array. The function multiplies the field
                    with the mask before the overlap. None applies no mask.

    Returns:
        The coupling efficiency, a float between 0 and 1.

    Note:
        The source kernel takes a LIST of fields and gives the ratio of the
        summed numerator to the summed denominator. This package propagates one
        field at a time, so the transcription takes one field and gives one
        float. Loop in the caller for a set of realisations.
    """
    E = field.field
    if mask is not None:
        E = E * mask
    M = smf_mode(field.siz, field.lam, field.N, aperture_m).field
    denominator = (np.abs(E) ** 2).sum()
    if denominator == 0.0:
        raise ValueError('coupling_efficiency: the field carries no power')
    numerator = np.abs((E * np.conj(M)).sum()) ** 2
    return float(numerator / denominator)


if __name__ == '__main__':
    from .field import Power
    from .sources import CircAperture, PlaneWave

    lam = 1550e-9
    D = 0.1                     # the pupil diameter
    size = 3 * D                # the grid holds the mode wings
    N = 512

    # The mode carries a unit intensity sum, and it is a real Gaussian.
    mode = smf_mode(size, lam, N, D)
    assert abs((np.abs(mode.field) ** 2).sum() - 1.0) < 1e-12
    c = N // 2
    r = mode.xvalues[c + N // 8]
    w_read = np.sqrt(r * r / np.log(abs(mode.field[c, c])
                                    / abs(mode.field[c, c + N // 8])))
    assert abs(w_read - D / MODE_RADIUS_RATIO) < 1e-6

    # A flat top-hat pupil couples at the Ruilier maximum, 0.8145.
    # Ruilier, DOI 10.1117/12.317094.
    flat = PlaneWave(Begin(size, lam, N), D)
    eta_flat = coupling_efficiency(flat, D)
    assert abs(eta_flat - 0.8145) < 0.01, eta_flat

    # The mode itself couples at 1.0. That is the physical limit.
    eta_mode = coupling_efficiency(mode, D)
    assert abs(eta_mode - 1.0) < 1e-9, eta_mode

    # A tilt across the pupil breaks the match, so the coupling falls.
    Y, X = flat.mgrid_cartesian
    tilted = SubIntensity(flat, Intensity(flat))     # a copy with the phase kept
    tilted.field = flat.field * np.exp(1j * 2 * np.pi * X / D)
    eta_tilt = coupling_efficiency(tilted, D)
    assert eta_tilt < 0.2 * eta_flat, eta_tilt

    # A mask that blocks the pupil centre also lowers the coupling.
    blocked = np.ones((N, N))
    blocked[flat.mgrid_Rsquared <= (0.3 * D / 2) ** 2] = 0.0
    eta_mask = coupling_efficiency(flat, D, mask=blocked)
    assert eta_mask < eta_flat, (eta_mask, eta_flat)

    print(f"pupil diameter D        {D * 1e3:9.2f} mm")
    print(f"mode radius D/2.24      {D / MODE_RADIUS_RATIO * 1e3:9.4f} mm")
    print(f"mode radius, read back  {w_read * 1e3:9.4f} mm")
    print("")
    print("coupling efficiency:")
    print(f"  top-hat pupil         {eta_flat:9.4f} (Ruilier max 0.8145)")
    print(f"  the mode itself       {eta_mode:9.4f}")
    print(f"  one-wave tilt         {eta_tilt:9.4f}")
    print(f"  obscured top hat      {eta_mask:9.4f}")
    print(f"  masked pupil power    {Power(flat):9.3e} W")
    print("self-check passed")
