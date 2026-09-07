"""The Noll Zernike modes and the circular aperture mask, in pure numpy.

The module gives three things:
- `noll_to_nm(j)`, the Noll single index j to the radial order n and the
  azimuthal order m.
- `zernike_j(N, j)`, one Noll mode on an N x N grid.
- `circle(N, d_px)`, a hard circular mask, with an optional central
  obscuration.

THE ORDER IS NOLL, NOT ANSI/OSA. The index starts at j = 1 (piston). j = 2 is
the x tilt and j = 3 is the y tilt. j = 4 is defocus. An EVEN j takes the
cosine term and an ODD j takes the sine term. This is the order of
R. J. Noll, "Zernike polynomials and atmospheric turbulence," J. Opt. Soc.
Am. 66(3), 207-211 (1976), DOI 10.1364/JOSA.66.000207, Table I, and it is the
order that `olb.turbulence.ao` counts: a TipTilt stage removes the first 3
Noll modes, and an AO(n) stage removes the first n Noll modes.

The normalisation is Noll's. Each mode has unit variance over the unit disc:

    Z_j = sqrt(n + 1) R_n^0(rho)                             for m = 0
    Z_j = sqrt(2 (n + 1)) R_n^m(rho) cos(m theta)            for m > 0
    Z_j = sqrt(2 (n + 1)) R_n^|m|(rho) sin(|m| theta)        for m < 0

    R_n^m(rho) = SUM_{s=0}^{(n-m)/2} (-1)^s (n-s)! rho^(n-2s)
                 / ( s! ((n+m)/2 - s)! ((n-m)/2 - s)! )

Source: Noll 1976, DOI 10.1364/JOSA.66.000207, Eqs. (2) to (4), printed
p. 208.

THE PIXEL CONVENTION IS THE CONVENTION OF THE LAYER. The axis is the pixel
(int(N/2), int(N/2)), and a mask keeps the pixels with r^2 <= radius^2. This is
the convention of `olb.waveoptics.sources.CircAperture` and of
`olb.waveoptics.turbulence.run._field_patch`. Do not change it. A stored field
patch and a mode raster must sit on the same pixels.
"""

from math import factorial

import numpy as np


def noll_to_nm(j):
    """Give the radial order n and the azimuthal order m of the Noll index j.

    The map is Noll's Table I. Source: Noll 1976,
    DOI 10.1364/JOSA.66.000207, Table I, printed p. 208. A positive m is the
    cosine term. A negative m is the sine term.

    Args:
        j: the Noll index. It starts at 1 (piston).

    Returns:
        The tuple (n, m).

    Raises:
        ValueError: j is less than 1.
    """
    j = int(j)
    if j < 1:
        raise ValueError(f"noll_to_nm: the Noll index j must be 1 or more, got {j}.")
    n = 0
    j1 = j - 1
    while j1 > n:
        n += 1
        j1 -= n
    # Noll's rule for the sign of m: an even j takes the cosine term.
    m = (-1) ** j * ((n % 2) + 2 * int((j1 + ((n + 1) % 2)) / 2))
    return n, int(m)


def _radial(n, m, rho):
    """Give the radial polynomial R_n^m(rho).

    Source: Noll 1976, DOI 10.1364/JOSA.66.000207, Eq. (4), printed p. 208.

    Args:
        n:   the radial order.
        m:   the absolute azimuthal order.
        rho: the normalised radius, an ndarray.

    Returns:
        An ndarray of the same shape as rho.
    """
    out = np.zeros_like(rho)
    for s in range((n - m) // 2 + 1):
        c = ((-1) ** s * factorial(n - s)
             / (factorial(s) * factorial((n + m) // 2 - s)
                * factorial((n - m) // 2 - s)))
        out += c * rho ** (n - 2 * s)
    return out


def _coordinates(N, radius_px):
    """Give the normalised radius and the angle of each pixel.

    The axis is the pixel (int(N/2), int(N/2)). See the module docstring.

    Args:
        N:         the number of pixels along one side.
        radius_px: the radius that maps to rho = 1, in pixels.

    Returns:
        The tuple (rho, theta).
    """
    c = int(N / 2)
    Y, X = np.mgrid[:N, :N]
    x = (X - c) / float(radius_px)
    y = (Y - c) / float(radius_px)
    return np.sqrt(x * x + y * y), np.arctan2(y, x)


def zernike_j(N, j, mask=None, radius_px=None):
    """Give the raster of the Noll mode j on an N x N grid.

    The mode is zero outside the unit circle. It is also zero outside the
    mask, when the caller gives a mask.

    Args:
        N:         the number of pixels along one side.
        j:         the Noll index. It starts at 1 (piston).
        mask:      an (N, N) bool array, or None. The mode is zero outside it.
        radius_px: the radius that maps to rho = 1, in pixels. The default is
                   N/2, the circle that the grid inscribes. Give the radius of
                   the aperture when the aperture is smaller than the grid.

    Returns:
        An (N, N) float64 array.
    """
    N = int(N)
    if radius_px is None:
        radius_px = N / 2.0
    n, m = noll_to_nm(j)
    rho, theta = _coordinates(N, radius_px)
    am = abs(m)
    out = _radial(n, am, rho)
    if m > 0:
        out = out * (np.sqrt(2.0 * (n + 1)) * np.cos(am * theta))
    elif m < 0:
        out = out * (np.sqrt(2.0 * (n + 1)) * np.sin(am * theta))
    else:
        out = out * np.sqrt(n + 1.0)
    out[rho > 1.0] = 0.0
    if mask is not None:
        out = np.where(mask, out, 0.0)
    return out


def zernike_basis(N, n_modes, mask, radius_px=None):
    """Give the Noll modes j = 1 .. n_modes at the in-mask pixels.

    Args:
        N:         the number of pixels along one side.
        n_modes:   the number of Noll modes. Mode 1 is piston.
        mask:      an (N, N) bool array of the aperture.
        radius_px: the radius that maps to rho = 1, in pixels. The default is
                   N/2.

    Returns:
        A (n_pix, n_modes) float64 array. The row order is the C order of the
        in-mask pixels of the grid.
    """
    mask = np.asarray(mask, dtype=bool)
    idx = np.flatnonzero(mask.ravel())
    out = np.empty((idx.size, int(n_modes)), dtype=np.float64)
    for k in range(int(n_modes)):
        out[:, k] = zernike_j(N, k + 1, radius_px=radius_px).ravel()[idx]
    return out


def circle(N, d_px, obscuration_ratio=0.0):
    """Give a hard circular mask of the diameter d_px pixels.

    The mask keeps the pixels with r^2 <= (d_px/2)^2, around the axis pixel
    (int(N/2), int(N/2)). This is the rule of
    `olb.waveoptics.sources.CircAperture`. An obscuration ratio cuts the
    central disc of the diameter obscuration_ratio * d_px.

    Args:
        N:                 the number of pixels along one side.
        d_px:              the diameter of the aperture, in pixels.
        obscuration_ratio: the ratio of the central obscuration diameter to
                           d_px. Use 0.0 for no obscuration.

    Returns:
        An (N, N) bool array.

    Raises:
        ValueError: the obscuration ratio is outside [0, 1).
    """
    N = int(N)
    if not 0.0 <= obscuration_ratio < 1.0:
        raise ValueError(
            "circle: the obscuration ratio must be 0.0 or more and less than "
            f"1.0, got {obscuration_ratio}.")
    c = int(N / 2)
    Y, X = np.mgrid[:N, :N]
    r2 = (X - c) ** 2.0 + (Y - c) ** 2.0
    out = r2 <= (d_px / 2.0) ** 2
    if obscuration_ratio > 0.0:
        out &= r2 > (obscuration_ratio * d_px / 2.0) ** 2
    return out


def mask_radius_px(mask):
    """Give the radius of the outermost in-mask pixel, in pixels.

    The radius comes from the axis pixel (int(N/2), int(N/2)). A modal fit
    uses it to normalise the Zernike modes over the aperture.

    Args:
        mask: an (N, N) bool array.

    Returns:
        The radius, in pixels, as a float.

    Raises:
        ValueError: the mask is empty.
    """
    mask = np.asarray(mask, dtype=bool)
    N = mask.shape[0]
    c = int(N / 2)
    Y, X = np.mgrid[:N, :N]
    r2 = (X - c) ** 2.0 + (Y - c) ** 2.0
    if not mask.any():
        raise ValueError("mask_radius_px: the mask holds no pixel.")
    return float(np.sqrt(r2[mask].max()))


if __name__ == '__main__':
    # 1. The index map against Noll's Table I, j = 1 to 15.
    # Source: Noll 1976, DOI 10.1364/JOSA.66.000207, Table I, printed p. 208.
    table = {1: (0, 0), 2: (1, 1), 3: (1, -1), 4: (2, 0), 5: (2, -2),
             6: (2, 2), 7: (3, -1), 8: (3, 1), 9: (3, -3), 10: (3, 3),
             11: (4, 0), 12: (4, 2), 13: (4, -2), 14: (4, 4), 15: (4, -4)}
    for j, nm in table.items():
        assert noll_to_nm(j) == nm, (j, noll_to_nm(j), nm)
    print("Noll index map j = 1 to 15: OK")

    # 2. The normalisation gives unit variance over the disc.
    N = 512
    disc = circle(N, N)
    npix = int(disc.sum())
    worst = 0.0
    for j in range(1, 16):
        Z = zernike_j(N, j, mask=disc)
        v = float((Z[disc] ** 2).mean())
        worst = max(worst, abs(v - 1.0))
    assert worst < 1e-2, worst
    print(f"unit variance over the disc, worst error {worst:.2e} (N = {N})")

    # 3. Two different modes are near-orthogonal on the discrete grid.
    B = zernike_basis(N, 15, disc)
    G = (B.T @ B) / npix
    off = np.abs(G - np.diag(np.diag(G))).max()
    assert off < 1e-2, off
    print(f"largest off-diagonal overlap {off:.2e}")

    # 4. The pixel count of a circle follows the area, and an obscuration
    #    removes the inner area.
    d = 200
    full = circle(N, d)
    ann = circle(N, d, obscuration_ratio=0.3)
    area = np.pi * (d / 2.0) ** 2
    assert abs(full.sum() - area) / area < 1e-2, full.sum()
    expect = area * (1.0 - 0.3 ** 2)
    assert abs(ann.sum() - expect) / expect < 1e-2, ann.sum()
    print(f"circle d = {d} px: {int(full.sum())} pixels (area {area:.0f})")
    print(f"obscured 0.3:      {int(ann.sum())} pixels (area {expect:.0f})")

    # 5. The axis pixel is the pixel int(N/2). Tip is odd about that pixel.
    c = int(N / 2)
    tip = zernike_j(N, 2, mask=disc)
    assert tip[c, c] == 0.0
    assert tip[c, c + 10] > 0.0 and tip[c, c - 10] < 0.0
    tilt = zernike_j(N, 3, mask=disc)
    assert tilt[c + 10, c] > 0.0 and tilt[c - 10, c] < 0.0
    print("axis pixel and tilt signs: OK")

    print("self-check passed")
