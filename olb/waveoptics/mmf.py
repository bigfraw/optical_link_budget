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


def focal_intensity(field, focal_length_m, numerical_aperture=None, mask=None,
                    defocus_m=0.0):
    """Focus a pupil field to the detector plane and give the intensity.

    The focal-plane amplitude is the 2-D Fourier transform of the pupil field
    (Fraunhofer diffraction, Goodman, ISBN 978-0974707723). A thin lens makes
    that transform, so this single FFT IS the physical propagation to the focal
    plane of an ideal lens. This helper does the focus alone: it applies the
    optional mask, it applies the optional numerical-aperture pupil gate, it
    applies the optional defocus, it focuses, and it returns the intensity with
    the focal pixel size. mmf_coupling_efficiency uses it, and a caller that
    wants the detector-plane image (for a picture) uses the same helper, so the
    FFT focus lives in one place.

    The numerical-aperture gate is a PUPIL amplitude mask. A ray from the pupil
    radius rho focuses at the angle rho/focal_length_m. So the fibre guides only
    the rays with rho <= focal_length_m*numerical_aperture. See Snyder and Love,
    DOI 10.1007/978-1-4613-2813-1. None applies no gate.

    A non-zero defocus moves the observation plane to z = f + defocus_m. A
    displaced plane is a QUADRATIC PHASE across the pupil,
    W(rho) = -pi*defocus_m*rho^2/(lambda*f^2) rad. The Fraunhofer transform of the
    pupil field times this phase is the physical field at the displaced plane
    (Goodman, ISBN 978-0974707723, defocus as a quadratic pupil aberration). The
    MINUS sign follows the phase convention of this port: a DIVERGING beam carries
    exp(+i*k*r^2/2R) (see olb.waveoptics.propagators.GForvard) and a lens applies
    exp(-i*k*r^2/2f) (see olb.waveoptics.lenses.Lens). So a plane BEYOND the focus
    (defocus_m > 0) needs a weaker lens, that is a positive residual pupil
    curvature radius, which is the minus sign here. So
    the spot grows. This route holds while the defocused spot stays inside the FFT
    window N*lambda*f/siz; the caller (mmf_coupling_efficiency) warns near that
    limit. defocus_m=0.0 is the focal plane (unchanged).

    Args:
        field:              the PUPIL-plane Field. The grid, the wavelength and
                            the pixel count come from this field.
        focal_length_m:     the focal length of the coupling optic, in m.
        numerical_aperture: the fibre numerical aperture. None applies no
                            angular gate.
        mask:               an optional N x N array. The function multiplies the
                            field with the mask before the focus. None applies no
                            mask.
        defocus_m:          the detector offset from the focal plane, in m
                            (z = f + defocus_m). 0.0 is the focal plane.

    Returns:
        A tuple (If, dx_focal). If is the N x N intensity |A|^2 at the detector
        plane. dx_focal is the focal pixel size, in m.

    Note:
        norm='ortho' keeps Parseval exact, so sum(If) equals the summed power of
        the gated pupil field. The defocus phase keeps the power.
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

    # Defocus. The detector plane is z = f + defocus_m. A displaced plane is a
    # quadratic pupil phase W(rho) = -pi*defocus_m*rho^2/(lambda*f^2) rad, and the
    # Fraunhofer transform then gives the physical field at that plane. Goodman,
    # ISBN 978-0974707723. The MINUS sign is the phase convention of this port
    # (diverging = +i*k*r^2/2R, lens = -i*k*r^2/2f; see propagators.GForvard and
    # lenses.Lens). defocus_m=0.0 leaves the focal-plane field unchanged.
    if defocus_m != 0.0:
        rho2 = field.mgrid_Rsquared
        Eg = Eg * np.exp(-1j * np.pi * defocus_m * rho2
                         / (field.lam * focal_length_m ** 2))

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
                            numerical_aperture=None, mask=None, defocus_m=0.0):
    """Calculate the power fraction that couples into a multimode fibre.

    eta = P_core / P_total. P_total is the collected pupil power. P_core is the
    detector-plane power inside the hard core disk, after the numerical-aperture
    gate. The focus is a Fraunhofer FFT (Goodman, ISBN 978-0974707723). So this
    is a light bucket: it sums the encircled energy inside the core, NOT a mode
    overlap.

    The turbulent field carries the tilt, so the spot walks off the core on its
    own. The receive mechanical jitter is NOT in the field; it is a separate
    analytic Term in the budget.

    defocus_m models a non-focal-plane detector (see the module and
    olb.models.coupling.terrestrial): it grows the spot. The detector plane is
    z = f + defocus_m, so the spot has the field-computed defocused shape (the
    AXIAL effect). See focal_intensity. It defaults to 0.0 (the focal plane).

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
                            It sets the defocused-spot window guard only.
        core_radius_m:      the fibre core radius, in m.
        focal_length_m:     the focal length of the coupling optic, in m.
        numerical_aperture: the fibre numerical aperture. None applies no
                            angular gate.
        mask:               an optional N x N array. The function multiplies the
                            field with the mask before the focus (the same mask
                            convention as olb.waveoptics.smf). None applies no
                            mask.
        defocus_m:          the detector offset from the focal plane, in m
                            (z = f + defocus_m). 0.0 is the focal plane. A
                            non-zero value grows the spot, so a fixed core
                            captures less.

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

    # Focus the pupil field with the shared helper: it applies the mask, the
    # numerical-aperture gate, and the defocus, and it gives the intensity and
    # pixel size at the detector plane.
    If, dx_focal = focal_intensity(field, focal_length_m,
                                   numerical_aperture=numerical_aperture,
                                   mask=mask, defocus_m=defocus_m)

    # The sampling guard. The core is resolved by core_radius_m/dx_focal focal
    # pixels along the radius. Below about 3 the disk integral is coarse.
    n_core = core_radius_m / dx_focal
    if n_core < 3.0:
        warnings.warn(
            f"mmf_coupling_efficiency: the core radius spans only {n_core:.2f} "
            f"focal pixels (dx_focal={dx_focal:.3e} m). The core disk integral "
            f"is coarse below about 3 pixels. Use a wider grid.")

    # The defocus window guard. The defocused spot has the geometric radius
    # (aperture_m/2)*|defocus_m|/focal_length_m. The FFT window half-width is
    # (N//2)*dx_focal. When the spot fills the window the FFT aliases, so the core
    # power is not trustworthy. This is the pupil-limit regime; use a physical
    # co-moving propagation there. Warn, do not raise.
    N = field.N
    if defocus_m != 0.0:
        reach = (aperture_m / 2.0) * abs(defocus_m) / focal_length_m
        half_window = (N // 2) * dx_focal
        if reach > 0.5 * half_window:
            warnings.warn(
                f"mmf_coupling_efficiency: the defocused spot reaches "
                f"{reach:.3e} m, more than half the focal window "
                f"{half_window:.3e} m; the FFT may alias. Use a wider grid, a "
                f"smaller defocus, or a physical co-moving propagation.")

    # The core disk, on the axis. The FFT with fftshift puts the axis at pixel
    # N//2, so the detector grid is zero-centred on that pixel with the spacing
    # dx_focal.
    c = N // 2
    idx = np.arange(N) - c
    XX, YY = np.meshgrid(idx, idx)
    x = XX * dx_focal
    y = YY * dx_focal
    r2_focal = x ** 2 + y ** 2
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

    # --- defocus (axial displacement) grows the spot ------------------------
    # Keep the defocused spot inside the FFT window: the geometric spot radius
    # (D/2)*dz/f must stay well below the window half-width. Here f=0.5, D=0.1,
    # so dz=1 mm gives a 100 um spot inside the 248 um half-window.
    a_def = 40e-6
    eta_focus = mmf_coupling_efficiency(flat, D, a_def, f)
    eta_defocus = mmf_coupling_efficiency(flat, D, a_def, f, defocus_m=0.5e-3)
    eta_more = mmf_coupling_efficiency(flat, D, a_def, f, defocus_m=1.0e-3)
    assert eta_defocus < eta_focus, (eta_defocus, eta_focus)
    assert eta_more < eta_defocus, (eta_more, eta_defocus)   # monotone spot growth
    # defocus_m=0.0 is exactly the focal-plane result (unchanged path).
    assert mmf_coupling_efficiency(flat, D, a_def, f, defocus_m=0.0) == eta_focus
    # The defocus phase keeps the total focal power (Parseval).
    If_def, _ = focal_intensity(flat, f, defocus_m=0.5e-3)
    assert np.isclose(If_def.sum(), (np.abs(flat.field) ** 2).sum(), rtol=1e-9)

    # The SIGN of defocus_m. A DIVERGING input (phase-front radius R > 0, so the
    # pupil carries exp(+i*k*rho^2/2R); see olb.waveoptics.propagators.GForvard)
    # focuses BEYOND the lens focal length, at z = f + f^2/(R-f) (thin-lens image
    # of a spherical input; S. A. Self, Appl. Opt. 22, 658 (1983),
    # DOI 10.1364/AO.22.000658). So the best coupling must sit at a POSITIVE
    # defocus_m. A flat pupil is symmetric in dz, so it cannot test the sign.
    k_test = 2.0 * np.pi / lam
    R_test = 200.0                                  # a diverging pupil, R > 0
    dz_true = f ** 2 / (R_test - f)
    curved = Field.copy(flat)
    curved.field = flat.field * np.exp(1j * k_test * flat.mgrid_Rsquared
                                       / (2.0 * R_test))
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        eta_true = mmf_coupling_efficiency(curved, D, a_def, f, defocus_m=dz_true)
        eta_at_f = mmf_coupling_efficiency(curved, D, a_def, f)
        eta_mirror = mmf_coupling_efficiency(curved, D, a_def, f,
                                             defocus_m=-dz_true)
    assert eta_true > eta_at_f > eta_mirror, (eta_true, eta_at_f, eta_mirror)

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
