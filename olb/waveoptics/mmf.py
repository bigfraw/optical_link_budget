"""The multimode-fibre (light-bucket) coupling reduction.

A multimode fibre is a light bucket. It has no single mode shape. The core is a
hard disk in the fibre plane. The fibre collects the power of the focused spot
that lands inside that disk. So the coupled fraction is the ENCIRCLED ENERGY of
the focal spot inside the core, after an angular acceptance gate.

This module is the MMF sibling of olb.waveoptics.smf. It focuses the received
pupil field to the focal plane with a Fraunhofer FFT, and it sums the focal
intensity inside the core disk. It is a LIGHT BUCKET, not a modal overlap.

The core is FIXED on the axis (the fibre points at boresight). The turbulent
field already carries the instantaneous angle-of-arrival tilt, so the focused
spot walks off the core on its own. So the fade is intrinsic here. This is
UNLIKE the analytic terrestrial MMF Term, which adds the received tilt by hand.
The receive MECHANICAL jitter is NOT in this efficiency: it is a separate
analytic Term in the budget. So this eta is a turbulence-only quantity.

The module stays self-contained: it imports numpy and the local field module
only. The caller passes plain floats.

Sources:
- Goodman, Introduction to Fourier Optics, ISBN 978-0974707723. The focal-plane
  field is the Fourier transform of the pupil field (Fraunhofer diffraction).
- Snyder and Love, Optical Waveguide Theory (1983),
  DOI 10.1007/978-1-4613-2813-1. The angular acceptance of a fibre: it guides
  a ray only when the ray angle is inside the acceptance cone of the numerical
  aperture.
"""

import warnings

import numpy as np


def focal_intensity(field, focal_length_m, numerical_aperture=None, mask=None):
    """Focus a pupil field to the focal plane and give the focal intensity.

    The focal-plane amplitude is the 2-D Fourier transform of the pupil field
    (Fraunhofer diffraction, Goodman, ISBN 978-0974707723). This helper does the
    focus alone: it applies the optional mask, it applies the optional
    numerical-aperture pupil gate, it focuses, and it returns the focal intensity
    with the focal pixel size. mmf_coupling_efficiency uses it, and a caller that
    wants the focal-plane image (for a picture) uses the same helper, so the FFT
    focus lives in one place.

    The numerical-aperture gate is a PUPIL amplitude mask. A ray from the pupil
    radius rho focuses at the angle rho/focal_length_m. So the fibre guides only
    the rays with rho <= focal_length_m*numerical_aperture. See Snyder and Love,
    DOI 10.1007/978-1-4613-2813-1. None applies no gate.

    Args:
        field:              the PUPIL-plane Field. The grid, the wavelength and
                            the pixel count come from this field.
        focal_length_m:     the focal length of the coupling optic, in m.
        numerical_aperture: the fibre numerical aperture. None applies no
                            angular gate.
        mask:               an optional N x N array. The function multiplies the
                            field with the mask before the focus. None applies no
                            mask.

    Returns:
        A tuple (If, dx_focal). If is the N x N focal intensity |A|^2. dx_focal
        is the focal pixel size, in m.

    Note:
        norm='ortho' keeps Parseval exact, so sum(If) equals the summed power of
        the gated pupil field.
    """
    E = field.field
    if mask is not None:
        E = E * mask

    # The numerical-aperture gate. A ray from the pupil radius rho focuses at
    # the angle rho/focal_length_m, so the fibre guides only rho <= f*NA. See
    # Snyder and Love, DOI 10.1007/978-1-4613-2813-1.
    Eg = E
    if numerical_aperture is not None:
        rho = np.sqrt(field.mgrid_Rsquared)
        rho_max = focal_length_m * numerical_aperture
        Eg = np.where(rho <= rho_max, E, 0.0)

    # Focus the pupil field to the focal plane. The focal-plane amplitude is the
    # 2-D Fourier transform of the pupil field (Fraunhofer diffraction, Goodman,
    # ISBN 978-0974707723). norm='ortho' keeps Parseval exact, so the sum of the
    # focal intensity equals the sum of the (gated) pupil intensity.
    A = np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(Eg), norm='ortho'))
    If = np.abs(A) ** 2

    # The focal pixel size. The focal position is x = lambda*f * spatial
    # frequency, and the FFT frequency spacing is 1/(N*dx_pupil) = 1/siz. So
    # dx_focal = lambda*f/siz. See Goodman, ISBN 978-0974707723.
    dx_focal = field.lam * focal_length_m / field.siz
    return If, dx_focal


def mmf_coupling_efficiency(field, aperture_m, core_radius_m, focal_length_m,
                            numerical_aperture=None, mask=None):
    """Calculate the power fraction that couples into a multimode fibre.

    eta = P_core / P_total. P_total is the collected pupil power. P_core is the
    focal power inside the hard core disk, after the numerical-aperture gate.
    The focus is a Fraunhofer FFT (Goodman, ISBN 978-0974707723). So this is a
    light bucket: it sums the encircled energy inside the core, NOT a mode
    overlap.

    The core is FIXED on the axis. The turbulent field carries the tilt, so the
    spot walks off the core on its own. The receive mechanical jitter is NOT
    here; it is a separate analytic Term in the budget. So this eta is
    turbulence-only.

    The numerical-aperture gate is a PUPIL amplitude mask. A ray from the pupil
    radius rho focuses at the angle rho/focal_length_m. So the fibre guides only
    the rays with rho <= focal_length_m*numerical_aperture. See Snyder and Love,
    DOI 10.1007/978-1-4613-2813-1. This gate on the field is the same physics as
    a flat (NA/NA_optic)^2 power factor, but it also captures how the gate
    reshapes (broadens) the focal spot, because it changes the field before the
    focus. None applies no gate (the old light-bucket-only behaviour).

    Args:
        field:              the received PUPIL-plane Field. The caller clips it
                            to the receive aperture first. The grid, the
                            wavelength and the pixel count come from this field.
        aperture_m:         the pupil DIAMETER, in m. Kept for the signature of
                            the coupling calls, and to match olb.waveoptics.smf.
                            The field is already clipped to the aperture, so the
                            value is not used in the integral.
        core_radius_m:      the fibre core radius, in m.
        focal_length_m:     the focal length of the coupling optic, in m.
        numerical_aperture: the fibre numerical aperture. None applies no
                            angular gate.
        mask:               an optional N x N array. The function multiplies the
                            field with the mask before the focus (the same mask
                            convention as olb.waveoptics.smf). None applies no
                            mask.

    Returns:
        The coupling efficiency, a float between 0 and 1. It is the fraction of
        the collected aperture power that enters the core, and it includes the
        numerical-aperture gate loss.

    Raises:
        ValueError: the field carries no power.

    Note:
        By Parseval with norm='ortho', sum(|A|^2) = sum(|Eg|^2) <= P_total, and
        P_core <= sum(|A|^2). So eta stays in [0, 1].
    """
    E = field.field
    if mask is not None:
        E = E * mask

    # The denominator is the total collected pupil power. It is the masked but
    # UNGATED field, so the numerical-aperture gate loss lowers the numerator
    # only, and it shows up in eta.
    p_total = float((np.abs(E) ** 2).sum())
    if p_total == 0.0:
        raise ValueError('mmf_coupling_efficiency: the field carries no power')

    # Focus the pupil field with the shared helper: it applies the mask and the
    # numerical-aperture gate, and it gives the focal intensity and pixel size.
    If, dx_focal = focal_intensity(field, focal_length_m,
                                   numerical_aperture=numerical_aperture,
                                   mask=mask)

    # The sampling guard. The core is resolved by core_radius_m/dx_focal focal
    # pixels along the radius. Below about 3 the disk integral is coarse.
    n_core = core_radius_m / dx_focal
    if n_core < 3.0:
        warnings.warn(
            f"mmf_coupling_efficiency: the core radius spans only {n_core:.2f} "
            f"focal pixels (dx_focal={dx_focal:.3e} m). The core disk integral "
            f"is coarse below about 3 pixels. Use a wider grid.")

    # The on-axis core disk. The FFT with fftshift puts the axis at pixel N//2,
    # so the focal grid is zero-centred on that pixel with the spacing dx_focal.
    N = field.N
    c = N // 2
    idx = np.arange(N) - c
    XX, YY = np.meshgrid(idx, idx)
    r2_focal = (XX ** 2 + YY ** 2) * dx_focal ** 2
    p_core = float(If[r2_focal <= core_radius_m ** 2].sum())

    return float(p_core / p_total)


if __name__ == '__main__':
    from .field import Begin, Field
    from .sources import CircAperture, PlaneWave

    lam = 1550e-9
    D = 0.1                     # the pupil diameter
    size = 8 * D                # the grid holds the focal wings
    N = 512
    f = 0.5                     # the focal length

    # The diffraction focal spot 1/e^2 radius, for the print lines.
    w_s = lam * f / (np.pi * D / 2.0)
    dx_focal = lam * f / size

    flat = PlaneWave(Begin(size, lam, N), D)

    # A very LARGE core (radius >> spot) captures almost all the power. A hard
    # top-hat pupil makes an Airy spot whose rings decay slowly (the encircled
    # energy goes to 1 as 1 - 2/(pi*v)), so the core must be many spot radii.
    a_large = 160e-6
    eta_large = mmf_coupling_efficiency(flat, D, a_large, f)
    assert eta_large > 0.98, eta_large

    # A very SMALL core (radius << spot) captures little.
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        eta_small = mmf_coupling_efficiency(flat, D, 1e-6, f)
    assert eta_small < 0.3 * eta_large, (eta_small, eta_large)

    # eta rises as the core radius grows (encircled energy is non-decreasing).
    cores = np.array([2, 4, 6, 8, 12]) * 1e-6
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')       # the 2 um core spans < 3 pixels
        etas = np.array([mmf_coupling_efficiency(flat, D, a, f) for a in cores])
    assert np.all(np.diff(etas) > 0.0), etas

    # The NA gate reduces eta. A tight NA (rho_max well inside D/2) cuts the
    # captured power. A generous NA (rho_max >= D/2) equals the ungated eta.
    a_mid = 40e-6
    na_optic = (D / 2.0) / f
    eta_none = mmf_coupling_efficiency(flat, D, a_mid, f)
    eta_tight = mmf_coupling_efficiency(flat, D, a_mid, f,
                                        numerical_aperture=0.3 * na_optic)
    eta_wide = mmf_coupling_efficiency(flat, D, a_mid, f,
                                       numerical_aperture=1.5 * na_optic)
    assert eta_tight < eta_none, (eta_tight, eta_none)
    assert np.isclose(eta_wide, eta_none), (eta_wide, eta_none)

    # A pupil-plane phase tilt walks the spot toward the core edge, so the
    # captured power falls. The tilt moves the spot by about 2 core radii.
    a_tilt = 10e-6
    eta_untilt = mmf_coupling_efficiency(flat, D, a_tilt, f)
    Y, X = flat.mgrid_cartesian
    scale = lam * f / (2.0 * a_tilt)          # focal shift = lambda*f/scale
    tilted = Field.copy(flat)
    tilted.field = flat.field * np.exp(1j * 2 * np.pi * X / scale)
    eta_tilt = mmf_coupling_efficiency(tilted, D, a_tilt, f)
    assert eta_tilt < eta_untilt, (eta_tilt, eta_untilt)

    # focal_intensity gives the same focal plane the efficiency uses. By Parseval
    # (norm='ortho') the summed focal intensity equals the gated pupil power, and
    # dx_focal is lambda*f/size. Goodman, ISBN 978-0974707723.
    If_flat, dx_f = focal_intensity(flat, f)
    assert np.isclose(dx_f, lam * f / size), (dx_f, lam * f / size)
    p_pupil = float((np.abs(flat.field) ** 2).sum())
    assert np.isclose(If_flat.sum(), p_pupil, rtol=1e-9), (If_flat.sum(), p_pupil)
    # The gate removes power, so the gated focal sum drops below the ungated one.
    If_gate, _ = focal_intensity(flat, f, numerical_aperture=0.3 * na_optic)
    assert If_gate.sum() < If_flat.sum(), (If_gate.sum(), If_flat.sum())

    # A zero-power field raises.
    zero = Field.copy(flat)
    zero.field = np.zeros((N, N))
    try:
        mmf_coupling_efficiency(zero, D, a_large, f)
        raise AssertionError('a zero-power field must raise')
    except ValueError:
        pass

    print(f"pupil diameter D        {D * 1e3:9.2f} mm")
    print(f"focal length f          {f * 1e3:9.1f} mm")
    print(f"spot radius w_s         {w_s * 1e6:9.3f} um")
    print(f"focal pixel dx_focal    {dx_focal * 1e6:9.3f} um")
    print("")
    print("coupling efficiency (light bucket):")
    print(f"  large core (160 um)   {eta_large:9.4f}")
    print(f"  small core (1 um)     {eta_small:9.4f}")
    print(f"  mid core, no gate     {eta_none:9.4f}")
    print(f"  mid core, tight NA    {eta_tight:9.4f}")
    print(f"  mid core, wide NA     {eta_wide:9.4f}")
    print(f"  10 um core, on axis   {eta_untilt:9.4f}")
    print(f"  10 um core, tilted    {eta_tilt:9.4f}")
    print("mmf coupling self-check passed")
