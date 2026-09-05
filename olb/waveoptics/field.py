"""The complex field on a square grid, and the field diagnostics.

Ported and trimmed from LightPipes (https://github.com/opticspy/lightpipes),
BSD-3-Clause. See LIGHTPIPES_LICENSE.txt in this package.

The module is pure physics. It imports numpy only. A field holds the complex
amplitude on an N x N square grid, the physical side of that grid, and the
wavelength. The grid is zero-centred: the pixel (N/2, N/2) is the axis.

Sources:
- Goodman, Introduction to Fourier Optics, ISBN 978-0974707723. Scalar
  diffraction theory, the intensity and the power of a scalar field.
- LightPipes manual, https://opticspy.github.io/lightpipes/manual.html.
  The implementation lineage and the grid convention.
"""

import copy as _copy

import numpy as np


class Field:
    """A scalar complex field on a square, zero-centred grid.

    THE PRECISION. The field holds complex128 by default. A caller asks for
    complex64 with dtype=numpy.complex64. The `field` setter casts each new
    array to that type, so the field keeps its precision after any operation.
    Single precision halves the bytes of each element. A large Monte Carlo is
    memory-bandwidth bound, so it runs faster. See the `precision` argument of
    olb.waveoptics.turbulence.run.propagate_turbulent_scenario.

    Attributes:
        field: the N x N complex amplitude array (complex128 or complex64).
        siz:   the physical side of the grid, in m.
        lam:   the wavelength, in m.
        N:     the number of pixels along one side.
        dx:    the distance between two pixels, in m.
    """

    @classmethod
    def begin(cls, grid_size, wavelength, N, dtype=np.complex128):
        """Make a new field. The amplitude is 1.0 at each pixel."""
        return cls(None, grid_size, wavelength, N, dtype=dtype)

    @classmethod
    def copy(cls, Fin):
        """Make a deep copy of a field. The copy shares no array."""
        return _copy.deepcopy(Fin)

    @classmethod
    def shallowcopy(cls, Fin):
        """Make a shallow copy of a field. The copy shares the array."""
        return _copy.copy(Fin)

    def __init__(self, Fin=None, grid_size=1.0, wavelength=1.0, N=0,
                 dtype=np.complex128):
        """Private. Use Begin() or the class methods.

        Raises:
            ValueError: dtype is not numpy.complex64 or numpy.complex128.
        """
        if np.dtype(dtype) not in (np.dtype(np.complex64),
                                   np.dtype(np.complex128)):
            raise ValueError(
                f"Field: dtype must be numpy.complex64 or numpy.complex128, "
                f"not {dtype!r}.")
        self._dtype = np.dtype(dtype).type
        if Fin is None:
            if not N:
                raise ValueError('Cannot create zero size field (N=0)')
            Fin = np.ones((N, N), dtype=self._dtype)
        else:
            Fin = np.asarray(Fin, dtype=self._dtype)
        self._field = Fin
        self._lam = wavelength
        self._siz = grid_size
        # The curvature of the spherical coordinate system, in 1/m. It is
        # 0.0 for a normal (flat) grid. LensForvard() and LensFresnel() set
        # it. Convert() removes it. See the LightPipes manual,
        # https://opticspy.github.io/lightpipes/manual.html, "Spherical
        # coordinates".
        self._curvature = 0.0
        # The Gaussian bookkeeping. GForvard reads these values.
        # q is the complex beam parameter, Siegman/Goodman convention:
        # q = z - i*z_R, with z_R = pi*w0^2/lam.
        self._IsGauss = False
        self._w0 = 0.2 * grid_size
        self._q = -1j * np.pi * self._w0 * self._w0 / wavelength
        self._z = 0.0
        self._A = 1.0

    # ---- the grid ----

    @property
    def grid_size(self):
        """The physical side of the grid, in m."""
        return self._siz

    @grid_size.setter
    def grid_size(self, value):
        self._siz = value

    siz = grid_size

    @property
    def wavelength(self):
        """The wavelength, in m."""
        return self._lam

    @wavelength.setter
    def wavelength(self, value):
        self._lam = value

    lam = wavelength

    @property
    def N(self):
        """The number of pixels along one side. The grid is square."""
        return self._field.shape[0]

    @property
    def dx(self):
        """The distance between two pixels, in m."""
        return self.siz / self.N

    @property
    def field(self):
        """The complex amplitude, an N x N array."""
        return self._field

    @field.setter
    def field(self, value):
        self._field = np.asarray(value, dtype=self._dtype)

    @property
    def xvalues(self):
        """The x coordinate of each pixel centre, a 1d array in m.

        The convention comes from matplotlib.pyplot.imshow: a positive
        shift in x is to the right. For an even N the centre pixel moves
        one step to the right.
        """
        w = self.N
        cx = int(w / 2)
        return self.dx * np.arange(-cx, (w - cx))

    @property
    def yvalues(self):
        """The y coordinate of each pixel centre, a 1d array in m.

        A positive shift in y is down.
        """
        h = self.N
        cy = int(h / 2)
        return self.dx * np.arange(-cy, (h - cy))

    @property
    def mgrid_cartesian(self):
        """The (Y, X) coordinate mesh of the grid, in m."""
        h, w = self.N, self.N
        cy, cx = int(h / 2), int(w / 2)
        Y, X = np.mgrid[:h, :w]
        return ((Y - cy) * self.dx, (X - cx) * self.dx)

    @property
    def mgrid_Rsquared(self):
        """The squared radius R^2 of each pixel, in m^2."""
        Y, X = self.mgrid_cartesian
        return X**2 + Y**2

    @property
    def mgrid_R(self):
        """The radius R of each pixel, in m."""
        return np.sqrt(self.mgrid_Rsquared)


def Begin(size, labda, N, dtype=np.complex128):
    """Make a new field of N x N pixels. The amplitude is 1.0.

    Args:
        size:  the physical side of the grid, in m.
        labda: the wavelength, in m.
        N:     the number of pixels along one side.
        dtype: numpy.complex128 (the default) or numpy.complex64. The field
               keeps this precision after each operation.

    Returns:
        A new Field.
    """
    return Field.begin(size, labda, N, dtype=dtype)


def field_dtype(precision):
    """Give the complex type of one precision name.

    The wave-optics runners take a `precision` keyword. This function is the
    ONE place that turns that name into a numpy type.

    Args:
        precision: "double" (numpy.complex128) or "single" (numpy.complex64).

    Returns:
        numpy.complex128 or numpy.complex64.

    Raises:
        ValueError: the name is not "double" or "single".
    """
    if precision == "double":
        return np.complex128
    if precision == "single":
        return np.complex64
    raise ValueError(
        f"precision must be 'double' or 'single', not {precision!r}.")


def Power(Fin):
    """Calculate the total power of the field.

    P = sum(|E|^2) * dx^2. See Goodman, Introduction to Fourier Optics,
    ISBN 978-0974707723 (the irradiance of a scalar field).
    """
    I = np.abs(Fin.field)**2
    return I.sum() * Fin.dx**2


def Normal(Fin):
    """Scale the field to a total power of 1.0.

    E_out = E_in / sqrt(P). See Goodman, ISBN 978-0974707723.
    """
    Fabs = np.abs(Fin.field)**2
    Fabs *= Fin.dx**2
    Ptot = Fabs.sum()
    if Ptot == 0.0:
        raise ValueError('Error in Normal(Fin): Zero beam power!')
    Fout = Field.copy(Fin)
    Fout.field = Fout.field * np.sqrt(1 / Ptot)
    return Fout


def Intensity(Fin, flag=0):
    """Calculate the intensity of the field.

    I(x,y) = E(x,y) * conj(E(x,y)). See Goodman, ISBN 978-0974707723.

    Args:
        Fin:  the input field.
        flag: 0 gives no normalisation. 1 normalises to 1. 2 normalises
              to 255 for a bitmap.

    Returns:
        An N x N array of real numbers.
    """
    I = np.abs(Fin.field)**2
    if flag > 0:
        Imax = I.max()
        if Imax == 0.0:
            raise ValueError('Cannot normalize because of 0 beam power.')
        I = I / Imax
        if flag == 2:
            I = I * 255
    return I


def Phase(Fin):
    """Calculate the phase of the field, in radians.

    The phase is the argument of the complex amplitude. See Goodman,
    ISBN 978-0974707723. The result stays in the interval [-pi, pi].
    """
    return np.angle(Fin.field)


def SubIntensity(Fin, Intens):
    """Replace the intensity of the field. The phase does not change.

    Args:
        Fin:    the input field.
        Intens: an N x N array of real numbers, or a scalar.

    Returns:
        A new Field.
    """
    Fout = Field.copy(Fin)
    Intens = np.asarray(Intens)
    if Intens.shape != Fout.field.shape:
        raise ValueError('Intensity map has wrong shape')
    phi = np.angle(Fout.field)
    Efield = np.sqrt(Intens)
    Fout.field = Efield * np.exp(1j * phi)
    Fout._IsGauss = False
    return Fout


if __name__ == '__main__':
    size = 20e-3
    lam = 1550e-9
    N = 256

    F = Begin(size, lam, N)
    assert F.N == N
    assert F.field.dtype == np.complex128
    assert abs(F.dx - size / N) < 1e-18
    assert abs(Power(F) - size**2) < 1e-12      # amplitude 1.0 on the grid

    # The grid is square and zero-centred.
    assert F.xvalues.shape == (N,)
    assert abs(F.mgrid_R[N // 2, N // 2]) < 1e-18

    # A fresh field is on a flat grid, and a copy keeps the curvature.
    assert F._curvature == 0.0
    F._curvature = -1e-6
    assert Field.copy(F)._curvature == -1e-6
    assert Field.shallowcopy(F)._curvature == -1e-6
    F._curvature = 0.0

    # Normal() always gives a power of 1.0.
    FN = Normal(F)
    assert abs(Power(FN) - 1.0) < 1e-12

    # Intensity and Phase.
    I = Intensity(FN, flag=1)
    assert abs(I.max() - 1.0) < 1e-12
    assert np.allclose(Phase(FN), 0.0)

    # SubIntensity keeps the phase and replaces the intensity.
    F2 = SubIntensity(FN, np.ones((N, N)))
    assert abs(Power(F2) - size**2) < 1e-12

    # ---- the single-precision switch ----
    # A complex64 field keeps complex64 after each operation, because the
    # `field` setter casts every new array to the stored type.
    F32 = Begin(size, lam, N, dtype=np.complex64)
    assert F32.field.dtype == np.complex64
    assert Field.copy(F32)._dtype == np.complex64
    assert Field.shallowcopy(F32)._dtype == np.complex64
    assert Normal(F32).field.dtype == np.complex64
    assert SubIntensity(Normal(F32),
                        np.ones((N, N))).field.dtype == np.complex64
    # A float64 array that a caller assigns comes back as complex64.
    F32.field = np.ones((N, N), dtype=np.complex128)
    assert F32.field.dtype == np.complex64
    # The two precisions agree on the power to about 1e-7.
    assert abs(Power(Normal(F32)) - 1.0) < 1e-5
    # field_dtype maps the two precision names, and it refuses each other name.
    assert field_dtype('double') is np.complex128
    assert field_dtype('single') is np.complex64
    try:
        field_dtype('half')
        raise AssertionError('field_dtype must refuse an unknown name')
    except ValueError as exc:
        assert 'single' in str(exc), str(exc)
    # An unknown dtype raises.
    try:
        Begin(size, lam, N, dtype=np.float64)
        raise AssertionError('Begin must refuse a real dtype')
    except ValueError as exc:
        assert 'complex64' in str(exc), str(exc)

    print(f"grid side          {size * 1e3:8.3f} mm")
    print(f"pixels per side    {N:8d}")
    print(f"pixel pitch        {F.dx * 1e6:8.3f} um")
    print(f"power, Begin       {Power(F):8.3e} W")
    print(f"power, Normal      {Power(FN):8.6f} W")
    print("self-check passed")
