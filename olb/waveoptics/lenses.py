"""The thin lens, and the propagators in spherical (co-moving) coordinates.

Ported and trimmed from LightPipes (https://github.com/opticspy/lightpipes),
BSD-3-Clause. See LIGHTPIPES_LICENSE.txt in this package.

The module is pure physics. It imports numpy only. It gives four functions:

- Lens:        an ideal thin lens. It is a quadratic phase mask.
- LensForvard: the spectral propagator in spherical coordinates.
- LensFresnel: the convolution propagator in spherical coordinates.
- Convert:     the return from spherical coordinates to a flat grid.

WHY the spherical route exists. A flat grid keeps its side. A space link
makes the beam grow by a factor of 100 or more, so the flat grid must hold
the far-field beam AND resolve the launch aperture. That needs a pixel count
which no computer holds. The spherical route moves the grid with the beam:
LensFresnel() rescales the grid side by (f - z)/f, and it keeps the phase
curvature in the field attribute _curvature. The pixel count stays small.

THE RECIPE for a diverging link, in three calls. The f of LensForvard() and
LensFresnel() is the focal length of the COORDINATE SYSTEM, not of a lens in
the beam. The functions add no phase to the field. So the caller puts the
equal and opposite PHYSICAL lens in the beam first. The two together add no
optical power, and the link stays the same link.

    m  = the grid magnification you want. For a Gaussian, m = w(z)/w0.
    fA = z / (m - 1)            the physical lens. It converges.

    F = Lens(F, fA)             it holds the beam on the small grid
    F = LensFresnel(F, -fA, z)  the coordinates diverge by m
    F = Convert(F)              it comes back to a flat grid

The identity behind the recipe is the ABCD factorisation of free space:

    [[1, z], [0, 1]] = Scale(m) . Free(z/m) . Lens(fA),   m = 1 + z/fA

The propagator does the SHORT step z/m on the launch grid. Then it relabels
the grid with the side m*size, and it keeps the residual curvature. See
Schmidt, DOI 10.1117/3.866274, Ch. 7 (the scaled Fresnel propagator), and
the LightPipes manual,
https://opticspy.github.io/lightpipes/manual.html, "Spherical coordinates".

The flat-grid propagators refuse a field that carries a curvature. Call
Convert() before you go back to Forvard() or Fresnel().

Sources:
- Goodman, Introduction to Fourier Optics, ISBN 978-0974707723. The thin-lens
  quadratic phase, and the focal spot of a Gaussian beam.
- Schmidt, Numerical Simulation of Optical Wave Propagation with Examples in
  MATLAB, DOI 10.1117/3.866274, Ch. 7 (the scaled two-step Fresnel
  propagator). The co-moving grid is the same idea: one internal step at a
  shorter distance, plus a coordinate scale factor.
- Siegman, Lasers, ISBN 978-0935702118. The ABCD matrix of a thin lens.
- LightPipes manual, https://opticspy.github.io/lightpipes/manual.html.
  The implementation lineage.
"""

import numpy as np

from .field import Field
from .propagators import Forvard, Fresnel, _ABCD

# The legacy constants of the C++ LightPipes. The port keeps them, so the
# numbers agree with the reference package.
_2PI_LEGACY = 3.1415926 * 2
_LARGENUMBER = 10000000.
_TINY_NUMBER = 1.0e-100


def Lens(Fin, f, x_shift=0.0, y_shift=0.0):
    """Put an ideal thin lens in the field.

    The lens multiplies the field with a quadratic phase:

        E_out(x,y) = E_in(x,y) * exp(-i*k*((x-dx)^2 + (y-dy)^2) / (2f))

    See Goodman, ISBN 978-0974707723 (the thin-lens phase transformation).

    A pure Gaussian beam with no shift takes the analytic ABCD route, with
    the lens matrix [[1, 0], [-1/f, 1]]. See Siegman, ISBN 978-0935702118.
    Each other case takes the phase mask, and the mask clears the Gaussian
    flag.

    Args:
        Fin:     the input field.
        f:       the focal length, in m. A negative f diverges the beam.
        x_shift: the shift of the lens centre in x, in m.
        y_shift: the shift of the lens centre in y, in m.

    Returns:
        A new Field.
    """
    if Fin._IsGauss and x_shift == 0.0 and y_shift == 0.0:
        return _ABCD(Fin, [[1.0, 0.0], [-1.0 / f, 1.0]])
    Fout = Field.copy(Fin)
    k = _2PI_LEGACY / Fout.lam
    yy, xx = Fout.mgrid_cartesian
    xx = xx - x_shift
    yy = yy - y_shift
    fi = -k * (xx**2 + yy**2) / (2 * f)
    Fout.field *= np.exp(1j * fi)
    Fout._IsGauss = False
    return Fout


def _combine_focal_lengths(f, curvature, size, lam):
    """Combine the virtual lens with the curvature the field already has.

    Two thin lenses in one plane add their powers: 1/f_total = 1/f + 1/f1.
    f1 comes from the curvature of the coordinate system. A flat grid gets a
    very long f1, so the result is f. See Goodman, ISBN 978-0974707723.
    """
    if curvature != 0.:
        f1 = 1. / curvature
    else:
        f1 = _LARGENUMBER * size**2 / lam
    if (f + f1) != 0.:
        return (f * f1) / (f + f1)
    return _LARGENUMBER * size**2 / lam


def LensForvard(Fin, f, z):
    """Propagate the field a distance z in spherical coordinates, with FFT.

    The function puts a virtual lens of the focal length f in the plane of
    the field, and it propagates a distance z. The grid side scales by
    (f - z)/f, and the internal spectral step is z1 = -z*f/(z - f). The
    residual curvature -1/(z - f) goes into the field. See the LightPipes
    manual, https://opticspy.github.io/lightpipes/manual.html, "Spherical
    coordinates", and Schmidt, DOI 10.1117/3.866274, Ch. 7.

    Use LensFresnel() for a long link. The spectral step of this function
    keeps the periodic artefact of Forvard().

    Args:
        Fin: the input field.
        f:   the focal length of the coordinate system, in m. Use
             f = -z / (m - 1) for a grid magnification m > 1, and put the
             physical lens Lens(F, z/(m - 1)) in the beam first. See the
             module docstring.
        z:   the propagation distance, in m.

    Returns:
        A new Field, in spherical coordinates.
    """
    size = Fin.siz
    lam = Fin.lam
    f = _combine_focal_lengths(f, Fin._curvature, size, lam)

    if (z - f) == 0:
        z1 = _LARGENUMBER
    else:
        z1 = -z * f / (z - f)

    Fout = Forvard(Fin, z1)

    ampl_scale = (f - z) / f
    size *= ampl_scale
    Fout._curvature = -1. / (z - f)
    Fout.siz = size
    if z1 >= 0:
        Fout.field /= ampl_scale
    else:
        # A backward internal step also mirrors the grid.
        ftemp = np.zeros_like(Fout.field, dtype=complex)
        ftemp.flat[:] = Fout.field.flat[::-1]
        Fout.field = ftemp
        Fout.field /= ampl_scale
    Fout._IsGauss = False
    return Fout


def LensFresnel(Fin, f, z):
    """Propagate the field a distance z in spherical coordinates.

    The function puts a virtual lens of the focal length f in the plane of
    the field, and it propagates a distance z with the convolution method.
    The grid side scales by (f - z)/f. The internal step is
    z1 = -z*f/(z - f), which is much shorter than z for a diverging virtual
    lens. The residual curvature -1/(z - f) goes into the field, so the
    output is NOT on a flat grid. Call Convert() to come back.

    THE RECIPE. f is the focal length of the COORDINATE SYSTEM. The function
    adds no phase to the field, so put the opposite PHYSICAL lens in the beam
    first. For a grid magnification m > 1:

        fA = z / (m - 1)
        F = Convert(LensFresnel(Lens(F, fA), -fA, z))

    For a Gaussian beam of the waist w0, take m = w(z)/w0, with
    w(z) = w0*sqrt(1 + (z/zR)^2) and zR = pi*w0^2/lam. The grid then holds
    the same number of beam radii at each end of the link. See Siegman,
    ISBN 978-0935702118, and the module docstring.

    The amplitude divides by the scale factor, so the power does not change:
    the grid area grows by the same factor squared.

    Args:
        Fin: the input field.
        f:   the focal length of the coordinate system, in m. A negative f
             makes the grid grow.
        z:   the propagation distance, in m.

    Returns:
        A new Field, in spherical coordinates.

    Raises:
        ValueError: the plane at z is behind the focus of the virtual lens.
    """
    size = Fin.siz
    lam = Fin.lam

    if f == z:
        f += _TINY_NUMBER

    f = _combine_focal_lengths(f, Fin._curvature, size, lam)

    z1 = -z * f / (z - f)
    if z1 < 0:
        raise ValueError('LensFresnel: behind the focus')

    Fout = Fresnel(Fin, z1)

    ampl_scale = (f - z) / f
    size *= ampl_scale
    Fout.siz = size
    Fout._curvature = -1. / (z - f)
    Fout.field /= ampl_scale
    Fout._IsGauss = False
    return Fout


def Convert(Fin):
    """Convert the field from spherical coordinates to a flat grid.

    The function removes the residual curvature of the coordinate system.
    It multiplies the field with the quadratic phase of a lens of the focal
    length f = -1/curvature, and it sets the curvature to zero. The grid
    side and the amplitude do not change. See Goodman,
    ISBN 978-0974707723, and the LightPipes manual,
    https://opticspy.github.io/lightpipes/manual.html.

    A field that is already on a flat grid comes back unchanged.

    Args:
        Fin: the input field.

    Returns:
        A new Field, on a flat grid.
    """
    # The copy happens for each field, so the return type is the same.
    Fout = Field.copy(Fin)
    curvature = Fin._curvature
    if curvature == 0.:
        return Fout

    f = -1. / curvature
    k = _2PI_LEGACY / Fin.lam
    kf = k / (2 * f)
    Fout.field *= np.exp(1j * kf * Fout.mgrid_Rsquared)
    Fout._curvature = 0.0
    Fout._IsGauss = False
    return Fout


if __name__ == '__main__':
    from .field import Begin, Power
    from .grid import GridSpec, forvard_max_z
    from .propagators import GForvard
    from .sources import CircAperture, GaussBeam

    def read_w(F):
        """Read the Gaussian radius w from a cut along the x axis.

        ln|E(r)| = ln|E(0)| - r^2/w^2, so a straight-line fit of ln|E|
        against r^2 gives w. Only the bright pixels take part.
        """
        c = F.N // 2
        r = F.xvalues[c:]
        a = np.abs(F.field[c, c:])
        keep = a > 0.1 * a[0]
        slope = np.polyfit(r[keep]**2, np.log(a[keep]), 1)[0]
        return np.sqrt(-1.0 / slope)

    def bucket(F, R):
        """The power inside a circle of the radius R, normalised."""
        return Power(CircAperture(F, R)) / Power(F)

    lam = 1550e-9

    # ---- 1. a thin lens makes the analytic focal spot ----
    w0_l = 1e-3
    f_lens = 0.1
    N_l = 1024
    # The aperture at 3 w0 clears the Gaussian flag, so Lens() takes the
    # phase-mask route. It removes 1e-8 of the power only.
    F = CircAperture(GaussBeam(Begin(8 * w0_l, lam, N_l), w0_l), 3 * w0_l)
    Ffoc = Fresnel(Lens(F, f_lens), f_lens)
    w_foc = read_w(Ffoc)
    w_foc_analytic = lam * f_lens / (np.pi * w0_l)   # Goodman
    assert abs(w_foc - w_foc_analytic) / w_foc_analytic < 0.02, w_foc

    # The ABCD branch of Lens() gives the same focal spot.
    Fg = GForvard(Lens(GaussBeam(Begin(8 * w0_l, lam, N_l), w0_l), f_lens),
                  f_lens)
    assert abs(read_w(Fg) - w_foc_analytic) / w_foc_analytic < 1e-6

    # ---- 2. the 600 km link on a co-moving grid ----
    w0 = 50e-3
    z = 600e3
    N = 512
    # The launch side is 10 w0. At 8 w0 the Fresnel convolution keeps about
    # 1 percent of power error at the edge of the grid.
    size0 = 10 * w0
    zR = np.pi * w0**2 / lam
    m = np.sqrt(1 + (z / zR)**2)          # w(z)/w0, the grid magnification
    fA = z / (m - 1)                      # the physical lens of the recipe

    F0 = GaussBeam(Begin(size0, lam, N), w0)
    # The analytic reference comes from the ABCD route.
    FG = GForvard(F0, z)
    wz = np.sqrt(-lam / np.pi * (FG._q.imag + FG._q.real**2 / FG._q.imag))
    assert abs(wz - w0 * m) / wz < 1e-12

    Flens = Lens(F0, fA)
    Fs = LensFresnel(Flens, -fA, z)
    assert Fs._curvature != 0.0
    assert abs(Fs.siz / size0 - m) / m < 1e-6

    # The coordinate rescale of LensFresnel adds no power: the amplitude
    # divides by m and the grid area grows by m^2.
    assert abs(Power(Fs) / Power(Fresnel(Flens, z / m)) - 1.0) < 1e-6

    Fz = Convert(Fs)
    assert Fz._curvature == 0.0
    # Convert is a pure phase. It adds no power.
    assert abs(Power(Fz) / Power(Fs) - 1.0) < 1e-12

    w_sim = read_w(Fz)
    assert abs(w_sim - wz) / wz < 0.01, (w_sim, wz)
    b_sim = bucket(Fz, wz)
    b_gauss = 1 - np.exp(-2)              # the Gaussian bucket at r = w
    assert abs(b_sim - b_gauss) / b_gauss < 0.01, (b_sim, b_gauss)
    # The Fresnel convolution keeps the power to about 1 percent here. The
    # error comes from the grid edge, not from the coordinate rescale.
    p_ratio = Power(Fz) / Power(F0)
    assert abs(p_ratio - 1.0) < 0.01, p_ratio

    # The flat grid cannot do this link. The spectral limit is short, and a
    # flat grid that holds w(z) at this pitch needs a huge pixel count.
    z_flat = forvard_max_z(GridSpec(size_m=size0, n=N), lam)
    assert z_flat < z / 1000
    n_flat = 8 * wz / (size0 / N)
    assert n_flat > 40000

    # ---- 3. after Convert the flat-grid propagators accept the field ----
    # The guard is happy now. The test is that no ValueError comes back.
    Fflat = Fresnel(Fz, 1e3)
    assert Fflat.field.shape == Fz.field.shape
    try:
        Fresnel(Fs, 1e3)
        raise AssertionError('Fresnel must refuse a spherical field')
    except ValueError as exc:
        assert 'Convert' in str(exc)
    # Convert on a flat field changes nothing.
    assert np.allclose(Convert(Fz).field, Fz.field)

    # ---- 4. LensForvard does the same bookkeeping ----
    # Forvard is a short-range propagator, so this is the bookkeeping check
    # only: the grid magnification and the power.
    z_s = 300.0
    m_s = np.sqrt(1 + (z_s / zR)**2)
    fA_s = z_s / (m_s - 1)
    Fv = LensForvard(Lens(F0, fA_s), -fA_s, z_s)
    assert abs(Fv.siz / size0 - m_s) / m_s < 1e-6
    assert abs(Power(Fv) / Power(F0) - 1.0) < 1e-9
    w_v = read_w(Convert(Fv))
    assert abs(w_v - w0 * m_s) / (w0 * m_s) < 0.01, (w_v, w0 * m_s)

    print("1. thin lens, focal spot")
    print(f"   waist w0                {w0_l * 1e3:12.3f} mm")
    print(f"   focal length f          {f_lens * 1e3:12.3f} mm")
    print(f"   w_focus, analytic       {w_foc_analytic * 1e6:12.3f} um")
    print(f"   w_focus, Lens+Fresnel   {w_foc * 1e6:12.3f} um")
    print(f"   w_focus, Lens ABCD      {read_w(Fg) * 1e6:12.3f} um")
    print("")
    print("2. the 600 km link on a co-moving grid")
    print(f"   waist w0                {w0 * 1e3:12.3f} mm")
    print(f"   range z                 {z * 1e-3:12.3f} km")
    print(f"   Rayleigh range zR       {zR * 1e-3:12.3f} km")
    print(f"   grid magnification m    {m:12.3f}")
    print(f"   physical lens fA        {fA * 1e-3:12.3f} km")
    print(f"   coordinate lens -fA     {-fA * 1e-3:12.3f} km")
    print(f"   internal step z/m       {z / m * 1e-3:12.3f} km")
    print(f"   launch grid side        {size0:12.4f} m  ({N} pixels)")
    print(f"   final grid side         {Fs.siz:12.4f} m  ({N} pixels)")
    print(f"   w(z), GForvard          {wz:12.4f} m")
    print(f"   w(z), co-moving grid    {w_sim:12.4f} m")
    print(f"   error                   {abs(w_sim - wz) / wz * 100:12.3f} %")
    print(f"   bucket at r = w(z)      {b_sim:12.6f}  (Gauss {b_gauss:.6f})")
    print(f"   power ratio             {p_ratio:12.6f}")
    print("")
    print("   for contrast, the same flat grid:")
    print(f"   Forvard limit N dx^2/lam{z_flat:12.1f} m  (link is {z:.0f} m)")
    print(f"   pixels to hold 8 w(z)   {n_flat:12.0f}     (grid has {N})")
    print("")
    print("3. LensForvard, 300 m, the bookkeeping only")
    print(f"   grid magnification      {Fv.siz / size0:12.6f}  (set {m_s:.6f})")
    print(f"   power ratio             {Power(Fv) / Power(F0):12.9f}")
    print("self-check passed")
