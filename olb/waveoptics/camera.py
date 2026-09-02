"""The focal-plane camera: the focused spot on finite pixels.

A tracking camera does not see the continuous focal intensity. It sees the
POWER IN EACH PIXEL. This module puts that discretisation on a fidelity-2
snapshot: it focuses the received pupil field, and it sums the focal power into
square camera pixels. So a wave-optics run gives the quantities a tracking loop
measures: the spot size in pixels, the spot centroid, and the power that spills
off the sensor.

The module REUSES the shared focal helper olb.waveoptics.mmf.focal_intensity. So
the Fraunhofer FFT, the defocus phase and its sign convention, and the Parseval
normalisation live in ONE place. This module adds the pixel grid alone.

THE BINNING. Camera pixel j along one axis covers the interval
    [(j - n_pixels/2)*pitch, (j - n_pixels/2 + 1)*pitch).
The fine focal sample at x goes to the pixel floor(x/pitch + n_pixels/2). The
module sums the samples of each pixel. It is an exact summation, not an
interpolation, so it keeps the power of every sample it accepts. A sample that
falls off the sensor is dropped: that power is the sensor spill.

THE NORMALISATION CONTRACT. camera_image divides the binned image by the TOTAL
masked pupil power. So image.sum() is the FRACTION of the collected power that
lands on the sensor, and it is 1.0 or less. spot_metrics reports that sum as
on_sensor_fraction.

The module stays self-contained. It imports numpy and the local mmf helper only.
It builds NO Term, and no budget reads it.

Sources:
- Goodman, Introduction to Fourier Optics, ISBN 978-0974707723. The focal-plane
  field is the Fourier transform of the pupil field (Fraunhofer diffraction),
  and a defocus is a quadratic pupil phase.
- ISO 11146-1:2021, Lasers and laser-related equipment - Test methods for laser
  beam widths, divergence angles and beam propagation ratios. The beam width and
  the beam position come from the first and second moments of the irradiance.
"""

import warnings
from dataclasses import dataclass

import numpy as np

from .mmf import focal_intensity


def camera_image(field, focal_length_m, pixel_pitch_m, n_pixels,
                 defocus_m=0.0, mask=None):
    """Focus a pupil field and bin it onto a square camera sensor.

    The function focuses the PUPIL field with the shared helper
    olb.waveoptics.mmf.focal_intensity (a Fraunhofer FFT; Goodman,
    ISBN 978-0974707723). It then sums the fine focal intensity into the camera
    pixels. The sensor is square, of n_pixels along each side, and it is centred
    on the axis.

    The fine focal grid has the pixel size dx_focal = lambda*f/siz, and it is
    zero-centred on the index N//2 (the convention of olb.waveoptics.mmf). The
    camera pixel j covers [(j - n_pixels/2)*pitch, (j - n_pixels/2 + 1)*pitch),
    so a fine sample at x goes to the pixel floor(x/pitch + n_pixels/2). The
    binning is an exact summation of the accepted samples.

    A non-zero defocus_m puts the sensor at z = f + defocus_m. The helper applies
    the quadratic pupil phase, so the spot grows. See focal_intensity for the
    sign convention: a DIVERGING received beam focuses at a POSITIVE defocus_m.

    Args:
        field:          the PUPIL-plane Field. The caller clips it to the
                        receive aperture first. The grid, the wavelength and the
                        pixel count come from this field.
        focal_length_m: the focal length of the imaging optic, in m.
        pixel_pitch_m:  the centre-to-centre pixel size of the camera, in m.
        n_pixels:       the number of camera pixels along one side.
        defocus_m:      the sensor offset from the focal plane, in m
                        (z = f + defocus_m). 0.0 is the focal plane.
        mask:           an optional N x N array. The function multiplies the
                        field with the mask before the focus (the mask
                        convention of olb.waveoptics.mmf). None applies no mask.

    Returns:
        A tuple (image, extent_m). image is the n_pixels x n_pixels array of the
        power fraction in each pixel; the first index is y and the second is x.
        image.sum() is the fraction of the collected power on the sensor.
        extent_m is the HALF-SIDE of the sensor, n_pixels*pixel_pitch_m/2, in m.
        It gives the imshow extent [-extent_m, extent_m, -extent_m, extent_m].

    Raises:
        ValueError: the field carries no power.

    Note:
        focal_intensity uses norm='ortho', so the summed fine focal intensity
        equals the summed masked pupil power (Parseval). A sensor that covers the
        whole fine focal window therefore reads image.sum() = 1.0.

        The optical axis falls on a pixel BOUNDARY, because n_pixels is usually
        even. The fine sample at x = 0 goes to the pixel on the positive side. So
        a symmetric spot reads a small positive centroid, of the order of a
        quarter of a pixel. This is the true response of an even-pixel sensor, not
        an error. A tracking loop calibrates that offset out.
    """
    E = field.field
    if mask is not None:
        E = E * mask

    # The denominator is the total collected pupil power. So image.sum() is the
    # fraction of that power on the sensor.
    p_total = float((np.abs(E) ** 2).sum())
    if p_total == 0.0:
        raise ValueError('camera_image: the field carries no power')

    # Focus the pupil field with the shared helper. It applies the mask and the
    # defocus, and it gives the intensity and the pixel size at the sensor plane.
    # A camera has no angular acceptance gate, so no numerical aperture is given.
    If, dx_focal = focal_intensity(field, focal_length_m, mask=mask,
                                   defocus_m=defocus_m)

    # The sampling guard. One camera pixel must hold several fine samples, or the
    # per-pixel sum is coarse and the pixel edges move by a full fine sample.
    if pixel_pitch_m < 3.0 * dx_focal:
        warnings.warn(
            f"camera_image: one camera pixel ({pixel_pitch_m:.3e} m) spans only "
            f"{pixel_pitch_m / dx_focal:.2f} fine focal samples "
            f"(dx_focal={dx_focal:.3e} m). The binning is coarse below about 3. "
            f"Use a wider pupil grid.")

    # The window guard. The FFT gives the focal intensity inside the half-width
    # (N//2)*dx_focal only. A sensor larger than that window reads zero on its
    # outer pixels, and that zero is false.
    N = field.N
    half_sensor = 0.5 * n_pixels * pixel_pitch_m
    half_window = (N // 2) * dx_focal
    if half_sensor > half_window:
        warnings.warn(
            f"camera_image: the sensor half-side {half_sensor:.3e} m is larger "
            f"than the focal window half-width {half_window:.3e} m. The outer "
            f"pixels read zero falsely. Use a wider pupil grid or fewer pixels.")

    # The fine focal coordinate of each sample, zero-centred on the index N//2
    # (the convention of olb.waveoptics.mmf).
    x_fine = (np.arange(N) - N // 2) * dx_focal

    # The camera pixel index of each fine sample. Pixel j covers
    # [(j - n_pixels/2)*pitch, (j - n_pixels/2 + 1)*pitch), so the index is the
    # floor of x/pitch + n_pixels/2. Keep the samples that land on the sensor.
    j = np.floor(x_fine / pixel_pitch_m + 0.5 * n_pixels).astype(int)
    keep = (j >= 0) & (j < n_pixels)

    image = np.zeros((n_pixels, n_pixels), dtype=float)
    if keep.any():
        jk = j[keep]
        # The first index is y (the rows), the second is x (the columns). This is
        # the index order of the focal intensity array.
        np.add.at(image, (jk[:, None], jk[None, :]), If[np.ix_(keep, keep)])

    return image / p_total, half_sensor


@dataclass(frozen=True)
class SpotMetrics:
    """The measurements a tracking camera makes on one binned image.

    Attributes:
        centroid_x_m:       the first moment of the pixel powers along x, in m.
                            It is zero on the sensor centre. Divide it by the
                            focal length to get the arrival angle.
        centroid_y_m:       the first moment along y, in m.
        rms_radius_m:       the root of the second central moment,
                            sqrt(<(x-xc)^2 + (y-yc)^2>), in m. It is the
                            second-moment spot size (ISO 11146-1:2021).
        peak_ix:            the column index of the brightest pixel.
        peak_iy:            the row index of the brightest pixel.
        on_sensor_fraction: the summed image. With the camera_image
                            normalisation it is the fraction of the collected
                            power that lands on the sensor.
    """

    centroid_x_m: float
    centroid_y_m: float
    rms_radius_m: float
    peak_ix: int
    peak_iy: int
    on_sensor_fraction: float


def spot_metrics(image, pixel_pitch_m):
    """Measure the spot position and the spot size of a binned camera image.

    The centroid is the FIRST moment of the pixel powers over the pixel-centre
    coordinates. The pixel j has the centre (j - n_pixels/2 + 0.5)*pitch, the
    same zero-centred grid that camera_image bins onto. The rms radius is the
    root of the SECOND CENTRAL moment. The first and second moments of the
    irradiance are the standard beam position and beam width (ISO 11146-1:2021).

    A camera measures a TRUNCATED moment: the pixels hold the power on the sensor
    only. So a spot that spills off the sensor reads a smaller rms radius than
    the true beam. on_sensor_fraction reports that truncation.

    Args:
        image:         the n x n array from camera_image (or any non-negative
                       pixel-power array). The first index is y.
        pixel_pitch_m: the centre-to-centre pixel size, in m.

    Returns:
        A SpotMetrics record.

    Raises:
        ValueError: the image carries no power.
    """
    img = np.asarray(image, dtype=float)
    total = float(img.sum())
    if total <= 0.0:
        raise ValueError('spot_metrics: the image carries no power')

    n = img.shape[0]
    # The pixel centres, on the zero-centred grid of camera_image.
    c = (np.arange(n) - 0.5 * n + 0.5) * pixel_pitch_m
    X, Y = np.meshgrid(c, c)          # X on the columns, Y on the rows

    xc = float((img * X).sum() / total)
    yc = float((img * Y).sum() / total)

    r2 = (X - xc) ** 2 + (Y - yc) ** 2
    rms = float(np.sqrt((img * r2).sum() / total))

    iy, ix = np.unravel_index(int(np.argmax(img)), img.shape)
    return SpotMetrics(centroid_x_m=xc, centroid_y_m=yc, rms_radius_m=rms,
                       peak_ix=int(ix), peak_iy=int(iy),
                       on_sensor_fraction=total)


if __name__ == '__main__':
    from .field import Begin, Field
    from .sources import PlaneWave

    lam = 1550e-9
    D = 0.1                     # the pupil diameter
    size = 8 * D                # the grid holds the focal wings
    N = 512
    f = 0.5                     # the focal length

    dx_focal = lam * f / size
    w_s = lam * f / (np.pi * D / 2.0)      # the diffraction spot 1/e^2 radius

    flat = PlaneWave(Begin(size, lam, N), D)

    # --- a large, well-sampled sensor holds almost all the power -------------
    pitch = 4.0 * dx_focal                 # 4 fine samples per camera pixel
    n_px = 126                             # the sensor stays inside the window
    img, extent = camera_image(flat, f, pitch, n_px)
    m = spot_metrics(img, pitch)
    assert m.on_sensor_fraction > 0.98, m.on_sensor_fraction
    assert abs(extent - 0.5 * n_px * pitch) < 1e-18
    # The spot sits on the axis, inside one pixel.
    assert abs(m.centroid_x_m) < pitch and abs(m.centroid_y_m) < pitch
    # The second-moment radius is of the order of the diffraction spot radius.
    # The Airy rings carry power far out, so the truncated moment reads several
    # spot radii (ISO 11146-1:2021 measures the true irradiance moment).
    assert w_s < m.rms_radius_m < 20.0 * w_s, (m.rms_radius_m, w_s)
    # The brightest pixel is the centre pixel.
    assert m.peak_ix == n_px // 2 and m.peak_iy == n_px // 2

    # --- a pupil tilt moves the centroid: the tracking measurement -----------
    # A pupil phase exp(i*2*pi*X/scale) is a tilt. The Fraunhofer transform moves
    # the spot by lambda*f/scale (Goodman, ISBN 978-0974707723).
    shift = 20e-6
    scale = lam * f / shift
    _, X_pupil = flat.mgrid_cartesian
    tilted = Field.copy(flat)
    tilted.field = flat.field * np.exp(1j * 2 * np.pi * X_pupil / scale)
    img_t, _ = camera_image(tilted, f, pitch, n_px)
    m_t = spot_metrics(img_t, pitch)
    assert abs(m_t.centroid_x_m - shift) < pitch, (m_t.centroid_x_m, shift)
    assert abs(m_t.centroid_y_m) < pitch, m_t.centroid_y_m

    # --- defocus grows the spot ---------------------------------------------
    # The geometric spot radius (D/2)*dz/f stays inside the focal window: f=0.5,
    # D=0.1, so dz=1 mm gives a 100 um spot inside the 248 um half-window.
    m_d1 = spot_metrics(camera_image(flat, f, pitch, n_px, defocus_m=0.5e-3)[0],
                        pitch)
    m_d2 = spot_metrics(camera_image(flat, f, pitch, n_px, defocus_m=1.0e-3)[0],
                        pitch)
    assert m.rms_radius_m < m_d1.rms_radius_m < m_d2.rms_radius_m, (
        m.rms_radius_m, m_d1.rms_radius_m, m_d2.rms_radius_m)
    # defocus_m=0.0 is exactly the focused image (the unchanged path).
    img0, _ = camera_image(flat, f, pitch, n_px, defocus_m=0.0)
    assert np.array_equal(img0, img)

    # --- a tiny sensor loses power ------------------------------------------
    img_small, _ = camera_image(flat, f, pitch, 6)
    m_small = spot_metrics(img_small, pitch)
    assert m_small.on_sensor_fraction < 0.9 * m.on_sensor_fraction, (
        m_small.on_sensor_fraction, m.on_sensor_fraction)

    # --- the coarse-pitch guard warns ---------------------------------------
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        camera_image(flat, f, 2.0 * dx_focal, 32)
    assert any('fine focal samples' in str(w.message) for w in caught), caught

    # --- the binning keeps the power ----------------------------------------
    # A sensor larger than the fine window accepts every fine sample. By Parseval
    # (norm='ortho' in focal_intensity) the summed focal intensity equals the
    # pupil power, so the image sums to 1.0. The window guard fires; ignore it.
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        img_all, _ = camera_image(flat, f, pitch, 256)
    assert abs(img_all.sum() - 1.0) < 1e-9, img_all.sum()

    # --- a zero-power field raises ------------------------------------------
    zero = Field.copy(flat)
    zero.field = np.zeros((N, N))
    try:
        camera_image(zero, f, pitch, n_px)
        raise AssertionError('a zero-power field must raise')
    except ValueError:
        pass

    # --- an empty image raises ----------------------------------------------
    try:
        spot_metrics(np.zeros((8, 8)), pitch)
        raise AssertionError('an empty image must raise')
    except ValueError:
        pass

    print(f"pupil diameter D        {D * 1e3:9.2f} mm")
    print(f"focal length f          {f * 1e3:9.1f} mm")
    print(f"spot radius w_s         {w_s * 1e6:9.3f} um")
    print(f"focal pixel dx_focal    {dx_focal * 1e6:9.3f} um")
    print(f"camera pitch            {pitch * 1e6:9.3f} um "
          f"({pitch / dx_focal:.1f} fine samples)")
    print("")
    print("case                 centroid x  centroid y   rms radius   on sensor")
    print("                         um          um           um                ")
    for name, mm in (("large sensor      ", m),
                     ("tilted pupil      ", m_t),
                     ("defocus 0.5 mm    ", m_d1),
                     ("defocus 1.0 mm    ", m_d2),
                     ("6 x 6 pixels      ", m_small)):
        print(f"{name} {mm.centroid_x_m * 1e6:10.3f}  "
              f"{mm.centroid_y_m * 1e6:10.3f}  {mm.rms_radius_m * 1e6:10.3f}  "
              f"{mm.on_sensor_fraction:10.4f}")
    print("camera self-check passed")
