"""The beam sources and the hard apertures.

Ported and trimmed from LightPipes (https://github.com/opticspy/lightpipes),
BSD-3-Clause. See LIGHTPIPES_LICENSE.txt in this package.

The module is pure physics. It imports numpy only. It keeps the fundamental
Gaussian mode, the circular plane wave, the circular aperture and the
circular screen. It drops the Hermite-Gauss, Laguerre-Gauss and doughnut
modes, and it drops the tilt.

Sources:
- Goodman, Introduction to Fourier Optics, ISBN 978-0974707723. The scalar
  field, the hard aperture and the Babinet complement.
- LightPipes manual, https://opticspy.github.io/lightpipes/manual.html.
  The implementation lineage.
"""

import numpy as np

from .field import Field


def CircAperture(Fin, R, x_shift=0.0, y_shift=0.0):
    """Put a circular aperture in the field.

    The aperture is a hard mask. It sets the field to zero outside the
    radius R. See Goodman, ISBN 978-0974707723 (the circular pupil).

    Args:
        Fin:     the input field.
        R:       the radius of the aperture, in m.
        x_shift: the shift of the centre in x, in m.
        y_shift: the shift of the centre in y, in m.

    Returns:
        A new Field.
    """
    Fout = Field.copy(Fin)
    Y, X = Fout.mgrid_cartesian
    Y = Y - y_shift
    X = X - x_shift
    dist_sq = X**2 + Y**2       # squared, so no sqrt is necessary
    Fout.field[dist_sq > R**2] = 0.0
    Fout._IsGauss = False
    return Fout


def CircScreen(Fin, R, x_shift=0.0, y_shift=0.0):
    """Put a circular screen in the field.

    The screen is the Babinet complement of the aperture. It sets the
    field to zero inside the radius R. See Goodman,
    ISBN 978-0974707723 (Babinet's principle).

    Args:
        Fin:     the input field.
        R:       the radius of the screen, in m.
        x_shift: the shift of the centre in x, in m.
        y_shift: the shift of the centre in y, in m.

    Returns:
        A new Field.
    """
    Fout = Field.copy(Fin)
    Y, X = Fout.mgrid_cartesian
    Y = Y - y_shift
    X = X - x_shift
    dist_sq = X**2 + Y**2
    Fout.field[dist_sq <= R**2] = 0.0
    Fout._IsGauss = False
    return Fout


def PlaneWave(Fin, w, x_shift=0.0, y_shift=0.0):
    """Make a circular plane wave of the diameter w.

    The plane wave is the input field behind a circular aperture of the
    radius w/2. Begin() gives a uniform amplitude of 1.0, so a plane wave
    on a fresh field is uniform. See Goodman, ISBN 978-0974707723.

    Args:
        Fin:     the input field.
        w:       the diameter of the plane wave, in m.
        x_shift: the shift of the centre in x, in m.
        y_shift: the shift of the centre in y, in m.

    Returns:
        A new Field.
    """
    Fout = Field.copy(Fin)
    Fout = CircAperture(Fout, w / 2, x_shift, y_shift)
    Fout._IsGauss = False
    return Fout


def GaussBeam(Fin, w0, x_shift=0.0, y_shift=0.0):
    """Make a fundamental Gaussian beam in its waist.

    E(x,y) = exp(-((x-dx)^2 + (y-dy)^2) / w0^2)

    The amplitude is 1.0 on the axis. w0 is the 1/e radius of the
    amplitude, so it is the 1/e^2 radius of the intensity. See Siegman,
    Lasers, ISBN 978-0935702118, p. 642, and Goodman,
    ISBN 978-0974707723.

    The function keeps the Gaussian bookkeeping of LightPipes when the
    shift is zero. Then GForvard can propagate the beam with the ABCD
    matrix.

    Args:
        Fin:     the input field.
        w0:      the waist radius, in m.
        x_shift: the shift of the waist centre in x, in m.
        y_shift: the shift of the waist centre in y, in m.

    Returns:
        A new Field.
    """
    Fout = Field.copy(Fin)
    Y, X = Fout.mgrid_cartesian
    Y = Y - y_shift
    X = X - x_shift
    Fout.field = np.exp(-(X * X + Y * Y) / (w0 * w0))

    if x_shift == 0.0 and y_shift == 0.0:
        # The ABCD route needs the complex beam parameter q = z - i*z_R,
        # with z_R = pi*w0^2/lam. See Siegman, ISBN 978-0935702118.
        Fout._IsGauss = True
        Fout._q = -1j * np.pi * w0 * w0 / Fout.lam
        Fout._w0 = w0
        Fout._z = 0.0
        Fout._A = 1.0
    else:
        Fout._IsGauss = False
    return Fout


if __name__ == '__main__':
    from .field import Begin, Normal, Power

    lam = 1550e-9
    N = 512
    w0 = 5e-3
    size = 16 * w0          # the grid is 8 waist radii each side

    F0 = Begin(size, lam, N)

    # Normal() gives a power of 1.0 for each source.
    for name, F in (("Gauss", GaussBeam(F0, w0)),
                    ("plane wave", PlaneWave(F0, 8 * w0)),
                    ("uniform", F0)):
        assert abs(Power(Normal(F)) - 1.0) < 1e-12, name

    # An aperture of 3 waist radii removes almost no Gaussian power.
    FG = Normal(GaussBeam(F0, w0))
    p_clip = Power(CircAperture(FG, 3 * w0))
    assert abs(1.0 - p_clip) < 1e-6

    # A uniform wave keeps the area fraction of the mask.
    FP = Normal(F0)                             # uniform, power 1.0
    R = 0.3 * size
    frac = Power(CircAperture(FP, R))
    expected = np.pi * R**2 / size**2
    assert abs(frac - expected) < 1e-3
    # PlaneWave(w) is the same mask, with the diameter w = 2R.
    assert abs(Power(PlaneWave(FP, 2 * R)) - frac) < 1e-15

    # The screen and the aperture are complements.
    total = Power(FG)
    assert abs(Power(CircScreen(FG, 2 * w0))
               + Power(CircAperture(FG, 2 * w0)) - total) < 1e-15

    # A shift moves the centre of mass.
    Y, X = F0.mgrid_cartesian
    I = np.abs(GaussBeam(F0, w0, x_shift=2 * w0).field)**2
    xc = (X * I).sum() / I.sum()
    assert abs(xc - 2 * w0) < 1e-6

    print(f"waist radius            {w0 * 1e3:8.3f} mm")
    print(f"grid side               {size * 1e3:8.3f} mm")
    print(f"Gauss power in 3 w0     {p_clip:8.6f}")
    print(f"plane wave area frac    {frac:8.5f} (exact {expected:8.5f})")
    print(f"shifted Gauss centroid  {xc * 1e3:8.3f} mm (set {2 * w0 * 1e3:.3f} mm)")
    print("self-check passed")
