# API reference: the wave-optics layer

This page documents `olb.waveoptics`, the fidelity-2 field propagation layer. It
gives the exact signatures and defaults from the source.

The package propagates a scalar complex field through free space on a square
grid. There is NO turbulence in the core. The turbulent split-step layer is the
sub-package `olb.waveoptics.turbulence`, and Section 9 documents it. The core is
a trimmed port of LightPipes
(https://github.com/opticspy/lightpipes, BSD-3-Clause). See
[`LIGHTPIPES_LICENSE.txt`](../olb/waveoptics/LIGHTPIPES_LICENSE.txt) in the
package. The port keeps the LightPipes names and the LightPipes call order, so a
script from that package runs here with no change.

Import from the package root:

```python
from olb.waveoptics import Begin, GaussBeam, Fresnel, GridSpec, propagate_scenario
```

Status: the layer IS wired into the budgets as `fidelity=2` (2026-08-28). The
bridge is `olb.models.waveoptics`: `run_fidelity2` makes the records, and the
Term factories turn them into Terms. Section 11 documents that module. The
modules of THIS package still build no Term and read no budget; they give the
fields and the per-trial scalars, and the bridge above does the rest. The vacuum
core stays the no-turbulence validator for the near-field and far-field limits of
the analytic Terms, and the turbulent sub-package gives snapshot statistics
against the same Terms. The remaining owner gate is whether fidelity 2 ever
becomes a DEFAULT.

The core (`field.py`, `sources.py`, `propagators.py`, `lenses.py`, `smf.py`,
`mmf.py`, `camera.py`) imports numpy and scipy only, and `threader.py` imports
the standard library only. They import nothing from the rest of `olb`. Only
`grid.py` and `run.py` read a scenario. The turbulent sub-package keeps the same
tiers (see Section 9).

---

## 1. The field (`olb/waveoptics/field.py`)

### `Field`

The scalar complex field on a square, zero-centred grid. The pixel `(N/2, N/2)`
is the axis. Do not call the constructor. Use `Begin()`.

| Member | Kind | Meaning |
|---|---|---|
| `field` | array | The `N` x `N` complex amplitude (`complex128`). |
| `siz` / `grid_size` | float | The physical side of the grid, in m. |
| `lam` / `wavelength` | float | The wavelength, in m. |
| `N` | int | The number of pixels along one side. |
| `dx` | float | The distance between two pixels, in m. `dx = siz / N`. |
| `xvalues`, `yvalues` | array | The coordinate of each pixel centre, in m. |
| `mgrid_cartesian` | `(Y, X)` | The coordinate mesh of the grid, in m. |
| `mgrid_Rsquared` | array | The squared radius of each pixel, in m^2. |
| `mgrid_R` | array | The radius of each pixel, in m. |
| `Field.copy(Fin)` | classmethod | A deep copy. The copy shares no array. |
| `Field.shallowcopy(Fin)` | classmethod | A shallow copy. The copy shares the array. |

The field also carries a private `_curvature`, in 1/m. It is `0.0` on a flat
grid. `LensForvard()` and `LensFresnel()` set it, and `Convert()` removes it. A
copy keeps it. See Section 3a.

### The field functions

- `Begin(size, labda, N, dtype=np.complex128)` — make a new `Field` of `N` x `N`
  pixels. The amplitude is 1.0 at each pixel. `dtype` is `numpy.complex128`
  (the default) or `numpy.complex64`; another type raises `ValueError`. The
  `field` setter casts every assignment to that type, so a field keeps its
  precision through every propagator. Each propagator builds its work arrays
  in the precision of the field it receives. `Forvard` keeps the phase wrap
  `Ir - Bus` in float64 on purpose, because that subtraction takes the
  fractional part of a number in the thousands, and it casts only the finished
  transfer function. `field_dtype(precision)` maps the name `"double"` or
  `"single"` to the numpy type.
- `Power(Fin)` — the total power, `P = sum(|E|^2) * dx^2`.
- `Normal(Fin)` — a new field that carries a total power of 1.0.
- `Intensity(Fin, flag=0)` — the intensity array `|E|^2`. `flag=1` normalises to
  1. `flag=2` normalises to 255 for a bitmap.
- `Phase(Fin)` — the phase array, in radians, in the interval `[-pi, pi]`.
- `SubIntensity(Fin, Intens)` — a new field with a replaced intensity. The phase
  does not change.

---

## 2. The sources and the apertures (`olb/waveoptics/sources.py`)

- `GaussBeam(Fin, w0, x_shift=0.0, y_shift=0.0)` — the fundamental Gaussian mode
  in its waist. `w0` is the 1/e amplitude radius. A zero shift keeps the Gaussian
  bookkeeping, so `GForvard` can propagate the beam.
- `PlaneWave(Fin, w, x_shift=0.0, y_shift=0.0)` — a circular plane wave of the
  DIAMETER `w`.
- `CircAperture(Fin, R, x_shift=0.0, y_shift=0.0)` — a hard circular aperture of
  the RADIUS `R`. It sets the field to zero outside `R`.
- `CircScreen(Fin, R, x_shift=0.0, y_shift=0.0)` — the Babinet complement. It sets
  the field to zero inside `R`.

Each mask clears the Gaussian flag, so `GForvard` refuses the field after it.

The port drops the Hermite-Gauss, Laguerre-Gauss and doughnut modes, and it drops
the tilt.

---

## 3. The propagators (`olb/waveoptics/propagators.py`)

- `Forvard(Fin, z)` — the FFT angular-spectrum method. A negative `z` propagates
  back. The grid keeps its side and its pitch.
- `Fresnel(Fin, z)` — the convolution method on a doubled grid. A negative `z`
  raises `ValueError`.
- `GForvard(Fin, z)` — the analytic ABCD route for a pure Gaussian beam. A field
  that is not a pure Gaussian raises `ValueError`.

Each function returns a new `Field`. A zero `z` returns a copy.

The three functions work on a FLAT grid only. Each one raises `ValueError` when
the field carries a curvature, because the grid side is then a function of the
coordinate system. Call `Convert()` first. LightPipes prints a message and gives
the field back unchanged; this port raises, because a silent bad result is worse
than a stop.

### 3a. The lens propagators (`olb/waveoptics/lenses.py`)

These four functions hold the thin lens and the spherical (co-moving) coordinate
route. The module imports numpy only.

- `Lens(Fin, f, x_shift=0.0, y_shift=0.0)` — an ideal thin lens. It multiplies the
  field with the quadratic phase `exp(-i*k*((x-dx)^2 + (y-dy)^2) / (2f))`
  (Goodman, ISBN 978-0974707723). `f` is in m, and a negative `f` diverges the
  beam. A pure Gaussian with no shift takes the analytic ABCD route with the lens
  matrix `[[1, 0], [-1/f, 1]]` (Siegman, ISBN 978-0935702118). Each other case
  takes the phase mask, and the mask clears the Gaussian flag.
- `LensForvard(Fin, f, z)` — the spectral (FFT) propagator in spherical
  coordinates. It puts a virtual lens of the focal length `f` in the plane of the
  field and propagates a distance `z`. The internal step is `z1 = -z*f/(z - f)`.
- `LensFresnel(Fin, f, z)` — the convolution propagator in spherical coordinates.
  The same bookkeeping as `LensForvard`, with the `Fresnel` internal step. A plane
  behind the focus of the virtual lens raises `ValueError`.
- `Convert(Fin)` — the return to a flat grid. It multiplies the field with the
  quadratic phase of a lens of the focal length `f = -1/curvature`, and it sets
  the curvature to zero. The grid side and the amplitude do not change. A field
  that is already flat comes back unchanged.

For `LensForvard` and `LensFresnel` the grid side scales by `(f - z)/f`, the
amplitude divides by the same factor, so the power does not change. The residual
curvature `-1/(z - f)` goes into the field. `f` is the focal length of the
COORDINATE SYSTEM, not of a lens in the beam, and the two functions add no phase
to the field.

`Lens(F, f)` with a curved input field combines the two powers first:
`1/f_total = 1/f + 1/f1`, with `f1 = 1/curvature`.

### The co-moving recipe

Put the equal and opposite PHYSICAL lens in the beam before the coordinate lens.
The two together add no optical power, so the link stays the same link. For a
grid magnification `m > 1`:

```python
m = w_z / w0                          # for a Gaussian, m = w(z)/w0
fA = z / (m - 1)                      # the physical lens. It converges.
F = Lens(F, fA)                       # it holds the beam on the small grid
F = LensFresnel(F, -fA, z)            # the coordinates diverge by m
F = Convert(F)                        # it comes back to a flat grid
```

The identity behind the recipe is the ABCD factorisation of free space,
`[[1, z], [0, 1]] = Scale(m) . Free(z/m) . Lens(fA)` with `m = 1 + z/fA`. The
propagator does the SHORT step `z/m` on the launch grid, then it relabels the
grid with the side `m*size`. See Schmidt, DOI 10.1117/3.866274, Ch. 7 (the scaled
Fresnel propagator), and the LightPipes manual,
https://opticspy.github.io/lightpipes/manual.html, "Spherical coordinates".

The self-check propagates a 50 mm waist over 600 km on a 512-pixel grid. The beam
radius agrees with the analytic ABCD value to 0.4 percent
(`python -m olb.waveoptics.lenses`).

---

## 4. The fibre coupling (`olb/waveoptics/smf.py`)

- `smf_mode(grid_size_m, wavelength_m, n, aperture_m)` — the back-propagated
  single-mode-fibre mode in the pupil plane. It is a Gaussian of the radius
  `aperture_m / MODE_RADIUS_RATIO`. The intensity sum is 1.0.
- `coupling_efficiency(field, aperture_m, mask=None, defocus_m=0.0,
  focal_length_m=None)` — the power fraction that
  couples into the fibre, a float between 0 and 1. It is the normalised overlap
  of the field with the mode. `aperture_m` is the pupil DIAMETER. An optional
  `mask` multiplies the field first. A field with no power raises `ValueError`.
  `defocus_m` puts the fibre tip at `z = f + defocus_m`: the function then
  multiplies the field with the quadratic pupil phase of
  `mmf.defocus_phase` (Section 4a), so a non-zero `defocus_m` needs
  `focal_length_m` and raises `ValueError` without it. The closed form of this
  overlap is `olb.models.coupling.smf_eta_defocused`, and the module self-check
  matches it to four decimals. `defocus_m=0.0` (the default) is the
  focal-plane overlap, unchanged.
- `MODE_RADIUS_RATIO = 2.24` — the best ratio of the pupil diameter to the
  pupil-plane mode radius (Ruilier, DOI 10.1117/12.317094). A flat pupil then
  couples at the maximum of 0.8145.

The function takes ONE field and gives ONE float. Loop in the caller for a set of
realisations.

### 4a. The multimode-fibre coupling (`olb/waveoptics/mmf.py`)

A multimode fibre is a LIGHT BUCKET, not a mode overlap. The core is a hard disk
in the fibre plane. The fibre collects the power of the focused spot that lands
inside that disk. So the coupled fraction is the ENCIRCLED ENERGY of the focal
spot inside the core.

- `mmf_coupling_efficiency(field, aperture_m, core_radius_m, focal_length_m,
  numerical_aperture=None, mask=None, defocus_m=0.0)` — the
  power fraction that couples into a
  multimode fibre, a float between 0 and 1. It is `eta = P_core / P_total`.
  `P_total` is the total collected pupil power. `P_core` is the detector-plane
  power inside the core disk of the radius `core_radius_m`, after the
  numerical-aperture gate. It focuses the pupil field with a Fraunhofer FFT
  (Goodman, ISBN 978-0974707723). The
  turbulent tilt walks the focused spot off the core on its own; the fade is
  intrinsic. The receive MECHANICAL jitter is NOT in this efficiency (it is a
  separate analytic Term), so this eta is a turbulence-only quantity.
  `aperture_m` is the pupil DIAMETER; it sets the defocused-spot window guard
  only. An
  optional `mask` multiplies the field first. A field with no power raises
  `ValueError`. The function WARNS when the core spans fewer than about 3 focal
  pixels.
- `focal_intensity(field, focal_length_m, numerical_aperture=None, mask=None,
  defocus_m=0.0)` —
  the shared helper that both `mmf_coupling_efficiency` and the example scripts
  use. It returns the tuple `(If, dx_focal)`. `If` is the detector-plane
  intensity,
  `|fftshift(fft2(ifftshift(Eg), norm='ortho'))|^2`, and `norm='ortho'` keeps
  Parseval exact. `dx_focal = field.lam * focal_length_m / field.siz` is the
  focal pixel size, in m. It applies the mask, the numerical-aperture pupil
  gate, and the defocus before the focus.
- `defocus_phase(field, defocus_m, focal_length_m)` — the quadratic pupil phase
  of the plane `z = f + defocus_m`,
  `exp(-i*pi*defocus_m*rho^2/(lambda*f^2))` (Goodman,
  ISBN 978-0974707723). It is the ONE shared defocus factor: `focal_intensity`
  and `olb.waveoptics.smf.coupling_efficiency` both call it, so the two
  fidelity-2 coupling legs read one sign convention. A `focal_length_m` of
  `None` raises `ValueError`.

#### The non-focal-plane detector: `defocus_m`

`defocus_m` moves the detector off the focal plane. It defaults to the
focal-plane behaviour, so an old call is unchanged.

- `defocus_m` puts the observation plane at `z = f + defocus_m`. A displaced
  plane is a QUADRATIC PHASE across the pupil,
  `W(rho) = -pi*defocus_m*rho^2/(lambda*f^2)` rad, and the Fraunhofer transform
  of the phased pupil is the physical field at that plane (Goodman,
  ISBN 978-0974707723; defocus as a quadratic pupil aberration). The MINUS sign
  is the phase convention of this port: a diverging beam carries
  `exp(+i*k*r^2/2R)` (`propagators.GForvard`) and a lens applies
  `exp(-i*k*r^2/2f)` (`lenses.Lens`). So the sign is now RIGHT-WAY-ROUND: a
  DIVERGING received beam couples best at a POSITIVE `defocus_m`, because a thin
  lens images a diverging input BEYOND its focal plane, at
  `z = f + f^2/(R - f)` (S. A. Self, Appl. Opt. 22, 658 (1983),
  DOI 10.1364/AO.22.000658). The module self-check asserts that sign. The phase
  keeps the power (Parseval).

The FFT route holds while the defocused spot stays inside the window
`N*lambda*f/siz`. `mmf_coupling_efficiency` WARNS when the geometric spot radius
passes half the window half-width, where the FFT can alias. Use a wider grid, a
smaller defocus, or a physical co-moving propagation there.

The numerical-aperture gate is a PUPIL amplitude mask. A ray from the pupil
radius `rho` focuses at the angle `rho/focal_length_m`, so the fibre guides only
the rays with `rho <= focal_length_m * numerical_aperture` (Snyder and Love,
DOI 10.1007/978-1-4613-2813-1). `None` applies no gate.

A SMALL receive aperture in a beam-sized grid gives a SHORT focal length and a
small focal field of view (about `lambda*f/dx_pupil`), so a plot window must stay
inside it; the focal integral is then also grid-limited.

### 4b. The focal-plane camera (`olb/waveoptics/camera.py`)

A tracking camera does not see the continuous focal intensity. It sees the POWER
IN EACH PIXEL. This module puts that discretisation on one fidelity-2 snapshot:
it focuses the received pupil field, and it sums the focal power into square
camera pixels. So a wave-optics run gives the quantities a tracking loop
measures: the spot size, the spot centroid, and the power that spills off the
sensor.

The module REUSES `focal_intensity` of Section 4a, so the Fraunhofer FFT, the
defocus phase and its sign convention, and the Parseval normalisation stay in ONE
place. This module adds the pixel grid alone. It imports numpy and that helper
only.

**DIAGNOSTIC ONLY.** The module builds NO Term, and no budget reads it. The
`Camera` detector of `olb.terminal` has no coupling model, so a `Camera` arm of
the runner (Section 9d) holds `None`.

- `camera_image(field, focal_length_m, pixel_pitch_m, n_pixels, defocus_m=0.0,
  mask=None)` — focus a PUPIL field and bin it onto a square sensor. The caller
  clips the field to the receive aperture first. It returns the tuple
  `(image, extent_m)`. `image` is the `n_pixels` x `n_pixels` array of the power
  fraction in each pixel; the first index is y and the second is x.
  `image.sum()` is the fraction of the COLLECTED power on the sensor, because the
  function divides by the total masked pupil power. `extent_m` is the HALF-SIDE
  of the sensor, `n_pixels*pixel_pitch_m/2`, for the `imshow` extent. A field
  with no power raises `ValueError`.
- `SpotMetrics` — a frozen dataclass with `centroid_x_m`, `centroid_y_m`,
  `rms_radius_m`, `peak_ix`, `peak_iy` and `on_sensor_fraction`. Divide a
  centroid by the focal length to get the arrival angle.
- `spot_metrics(image, pixel_pitch_m)` — the first and the second central moments
  of a binned image, as a `SpotMetrics`. The moments are the standard beam
  position and beam width (ISO 11146-1:2021). An empty image raises `ValueError`.

**THE BINNING.** Camera pixel `j` along one axis covers
`[(j - n_pixels/2)*pitch, (j - n_pixels/2 + 1)*pitch)`, and the fine focal sample
at `x` goes to the pixel `floor(x/pitch + n_pixels/2)`. The function SUMS the
samples of each pixel. It is an exact summation, not an interpolation. A sample
that falls off the sensor is dropped, and that power is the sensor spill.

`camera_image` WARNS two times. It warns when one camera pixel spans fewer than
about 3 fine focal samples (`dx_focal = lambda*f/siz`), because the binning is
then coarse. It warns when the sensor half-side is larger than the focal window
half-width `(N//2)*dx_focal`, because the outer pixels then read zero falsely.
Use a wider pupil grid in both cases.

The optical axis falls on a pixel BOUNDARY when `n_pixels` is even, so a
symmetric spot reads a small positive centroid, of the order of a quarter of a
pixel. That is the true response of an even-pixel sensor. A tracking loop
calibrates the offset out.

---

## 5. The grid (`olb/waveoptics/grid.py`)

### `GridSpec(size_m, n, scaled=False)`

A frozen dataclass. It holds the numbers that a propagation needs.

| Field | Type | Unit | Meaning |
|---|---|---|---|
| `size_m` | float | m | The physical side of the square grid. It is the side at the LAUNCH plane when `scaled` is True. |
| `n` | int | — | The number of pixels along one side. |
| `scaled` | bool | — | True selects the scaled (co-moving) route of Section 3a. The grid then grows with the beam, and the side at the receive plane is `size_m` times the magnification. False keeps a flat grid. |
| `pixel_m` | float | m | A property. The pitch, `size_m / n`. |

### `GridSpec.for_scenario(scenario, geometry, guard=4.0, pixels_per_feature=16, n_max=4096)`

Derive a grid from a scenario and a geometry. The method tries the FLAT route
first, and it takes the SCALED (co-moving) route when the flat route cannot
resolve the apertures.

- The FLAT EXTENT rule: `size = guard * 2 * r_max`. `r_max` is the largest of the
  beam radius at the launch plane, the beam radius at the longest range, the
  transmit aperture radius, and the receive aperture radius.
- The SCALED EXTENT rule: the grid starts at the launch plane, so
  `size = guard * 2 * max(waist, transmit aperture radius)` only. The grid then
  grows with the beam by the magnification `m`.
- The RESOLUTION rule: the smallest feature gets `pixels_per_feature` pixels
  across it. The features are the transmit waist and the hard edges of the two
  apertures (each aperture radius, and each central obscuration radius). The
  scaled route measures a receive feature at the LAUNCH plane, so it divides that
  feature by `m`. The pixel count goes up to the next power of two. It stays in
  the interval `[256, n_max]`.
- The transmit aperture obeys the bistatic rule of
  `olb.models.gaussian_efficiency`.

A beam that does not grow (`m = 1`) has no scaled route, because the recipe needs
a lens of the focal length `z/(m - 1)`.

The method WARNS. It does not raise. It warns when NEITHER route resolves the
smallest feature under the `n_max` clamp, and when the range of a flat grid is
longer than `forvard_max_z()`. A transmit terminal with no `Transmitter` raises
`ValueError`.

### `beam_magnification(scenario, z)`

The grid magnification `m = w(z)/w(0)` of the transmit beam, a float. It is 1.0 at
`z = 0`. A deliberately diverged beam brings its own `w(z)`, because
`olb.beam.free_space_radius` reads the `Transmitter` divergence. The sizer and
`propagate_scenario()` both call this function, so the grid grows by exactly the
factor that the beam grows by.

### `forvard_max_z(grid, wavelength_m)`

The longest range that the grid samples well, `z_max = N * dx^2 / lambda`. Past
this range the quadratic phase of the transfer function turns faster than one
sample, so `Forvard` aliases. See Schmidt, DOI 10.1117/3.866274, Ch. 6.

`N_MIN = 256` is the smallest pixel count.

---

## 6. One end-to-end propagation (`olb/waveoptics/run.py`)

### `propagate_scenario(scenario, geometry, grid=None, precision="single")`

Propagate the transmit beam of a scenario to the receive aperture. The steps are:
the launch, the launch-aperture clip, the free-space propagation, the
receive-aperture clip, and the fibre coupling. `grid=None` derives the grid with
`GridSpec.for_scenario()`. A geometry that gives more than one range raises
`ValueError`. Loop in the caller for a sweep.

A deliberately diverged beam starts at a virtual waist behind the aperture (see
`olb.beam`). Then the beam has the asked-for radius in the aperture plane.

The propagation step selects one of three routes: the exact ABCD route, the flat
`Fresnel` convolution, or the co-moving lens recipe. `GridSpec.for_scenario()`
selects the grid route, and `grid.scaled` records it. See Section 7.

### `WaveResult`

A frozen dataclass. All the losses are positive dB.

| Field | Type | Meaning |
|---|---|---|
| `stages` | list | The `(label, Field)` pairs. The labels are `"launch"`, `"after tx clip"`, `"at rx plane"` and `"after rx clip"`. |
| `grid` | `GridSpec` | The grid that the propagation used. |
| `tx_truncation_db` | float | The power that the launch aperture takes. |
| `geometric_loss_db` | float | The power that the receive aperture does not collect. |
| `smf_coupling_db` | float or None | The single-mode-fibre coupling loss. It reads `SMF.defocus_m` (the fibre tip at `z = f + defocus_m`). `None` when the receive terminal has no `SMF` detector. |
| `propagator` | str | The name of the propagator that ran: `"GForvard"`, `"Fresnel"` or `"LensFresnel"`. |

`PURE_GAUSS_CLIP = 1e-6` is the dispatch threshold. See Section 7.

### Compare the TOTAL, not the two parts

The fidelity-2 numbers do NOT compare one to one with the fidelity-0 analytic
Terms, because the two fidelities cut the loss in two different places:

- `olb.models.gaussian_efficiency.tx_efficiency_loss_db` is an on-axis FAR-FIELD
  gain ratio. `WaveResult.tx_truncation_db` is a plain power ratio at the launch
  aperture.
- `olb.models.geometric.geometric_loss_db` is the power fraction of the
  UNtruncated Gaussian in the receive aperture. `WaveResult.geometric_loss_db` is
  the power fraction of the PROPAGATED field.

The product of the fidelity-0 pair is the collected power fraction in the far
field with a small receiver. So the launch-to-collected TOTAL is the comparable
quantity, and the split is not. On the 600 km space link of
`examples/waveoptics/space_farfield.py` the two totals agree to 0.011 dB, and the
two splits do not.

---

## 7. The propagator regimes and their limits

| Propagator | Method | Field | Distance limits | Grid limits |
|---|---|---|---|---|
| `GForvard` | Analytic ABCD (Siegman, ISBN 978-0935702118) | A pure Gaussian beam ONLY. It refuses a clipped field. | Any `z`. | None. There is no grid error, because the route is analytic. |
| `Fresnel` | Convolution on a doubled grid (Schmidt, DOI 10.1117/3.866274, Ch. 7) | Any field. | MINIMUM `z`. The result is not valid when `z` is comparable with, or less than, the size of the diffracting aperture. (LightPipes manual; Schmidt Ch. 7) | No periodic wrap: the doubled grid absorbs it. The cost is 8 times the memory. The field must be zero at the grid edges. |
| `Forvard` | FFT angular spectrum (Schmidt, DOI 10.1117/3.866274, Ch. 6) | Any field. | Any `z` down to zero. But EACH call needs `z < forvard_max_z = N*dx^2/lambda`. Past that limit the transfer function aliases. | The boundary is periodic. A beam that reaches the edge wraps to the opposite edge. Give the grid a side of about 8 times the largest beam radius. |
| `LensFresnel` + `Convert` | The `Fresnel` convolution in spherical coordinates (Schmidt, DOI 10.1117/3.866274, Ch. 7; the LightPipes manual) | Any field. | The internal step is `z/m`, so the minimum-distance limit of `Fresnel` applies to `z/m`, not to `z`. A plane behind the focus of the virtual lens raises `ValueError`. | The grid GROWS with the beam by `m`. So one small pixel count holds a launch aperture and a far-field beam that is 100 times wider. The output is spherical until `Convert()` runs. |
| `LensForvard` + `Convert` | The `Forvard` spectral method in spherical coordinates | Any field. | The internal step is `z1 = -z*f/(z - f)`, and that step keeps the `forvard_max_z` limit and the periodic artefact of `Forvard`. | The same co-moving grid as `LensFresnel`. Use `LensFresnel` for a long link. |

The three flat-grid propagators refuse a spherical field. The two lens
propagators leave the field spherical, so call `Convert()` before you go back to
`Forvard`, `Fresnel` or `GForvard`. For the recipe that drives the two lens
propagators, see Section 3a.

### The dispatch rule of `propagate_scenario`

`propagate_scenario` selects one of three routes. It reads the power that the
launch aperture takes, and then the grid route:

- The clip takes less than `PURE_GAUSS_CLIP = 1e-6` of the power: the field stays
  a pure Gaussian. The function propagates the UNCLIPPED launch field with
  `GForvard`, the exact route, at any range.
- The clip takes more, and `grid.scaled` is False: the field carries the aperture
  edge on a flat grid. The function propagates the CLIPPED field with `Fresnel`.
- The clip takes more, and `grid.scaled` is True: the function runs the co-moving
  recipe on the CLIPPED field,
  `Convert(LensFresnel(Lens(F, fA), -fA, z))` with `fA = z/(m - 1)` and
  `m = beam_magnification(scenario, z)`. `WaveResult.propagator` is then
  `"LensFresnel"`. This is the route for a long space link.

`propagate_scenario` does not call `Forvard` or `LensForvard`.

---

## 8. The LightPipes propagators that the port does not hold

Each of these serves a regime that the current layer does not need. Add one only
when its trigger comes.

| Function | Regime | The trigger to add it |
|---|---|---|
| `Forward` | The direct integral. The output grid is different from the input grid. | You must magnify or shrink the grid between two planes. |
| `Steps` | Propagation through a medium with an index term. It holds a built-in absorbing boundary. | Nothing now. The turbulent split-step layer of Section 9 does the same work with `Forvard` and its own absorbing mask. |
| `Interpol` | A regrid of the field: a new side, a new pixel count, a shift, or a rotation. | You must pass a field between two propagators that need different grids. |

---

## 9. The turbulent split-step layer (`olb/waveoptics/turbulence/`)

The turbulent split-step layer is BUILT and self-checked. It moves a complex
field along a path and it puts a random phase screen at each slab of that path.
It gives one SNAPSHOT of the atmosphere for each seed. It carries NO time axis.

Status: the sub-package builds NO Term itself, but its records ARE wired into the
budgets as `fidelity=2` through `olb.models.waveoptics` (Section 11). The
remaining owner decision is whether fidelity 2 ever becomes a DEFAULT. See
`examples/waveoptics/README.md`.

Import from the sub-package:

```python
from olb.waveoptics.turbulence import (propagate_turbulent_scenario,
                                       turbulent_grid)
```

The seven modules are:

| Module | What it holds |
|---|---|
| `screens.py` | The phase screens: the fast `ScreenFactory` and the `aotools` wrapper. |
| `splitstep.py` | The propagate-screen-propagate loop, and the boundary mask. |
| `sampling.py` | The turbulent grid sizer, and the screen-placement planner. |
| `run.py` | The trial runner: one snapshot for each seed. |
| `campaign.py` | A large set of trials on disk, stored as blocks. |
| `fingerprint.py` | The content key that names one campaign (`cache_key`). |
| `temporal.py` | The frozen-flow time axis. PLANNED, NOT BUILT. |

**Two screen generators.** The DEFAULT is the `olb` generator, the fast
`ScreenFactory` of Section 9a. It imports numpy and scipy only, so a default
turbulent run needs NO `aotools`. The opt-in `aotools` generator is the
reference path; the caller selects it with `screen_generator="aotools"`.
`aotools` is a DEPENDENCY for that path only, and `screens.py` imports it
lazily. `olb` does not copy it, because `aotools` is LGPL-3.0. Install it with
`pip install aotools`, or with the extra `pip install olb[screens]`. An absent
`aotools` raises `ImportError` with that text, and only when the aotools path
draws a screen. The two generators give DIFFERENT random draws for the same
seed; the statistics agree (Section 9a).

The import tiers follow the tiers of the vacuum package. `screens.py` and
`splitstep.py` read the wave-optics core only. `sampling.py` and `run.py` read
the rest of `olb` (a scenario, the `Cn2` profiles, the Andrews layer).
`temporal.py` imports numpy only.

### 9a. The phase screens (`olb/waveoptics/turbulence/screens.py`)

- `screen_r0(cn2_integral_m13, wavelength_m)` — the Fried parameter of one slab,
  `r0 = (0.423 * k^2 * INT Cn2 dz)^(-3/5)` with `k = 2*pi/lambda`. The caller
  gives the INTEGRAL of `Cn2` over the slab, in m^(1/3), and that integral
  carries any slant factor already. The function adds no geometry. The input can
  be a scalar or an array, and the type of the result follows the input. See
  Fried, DOI 10.1364/JOSA.56.001372, and Andrews and Phillips,
  DOI 10.1117/3.626196, Ch. 12, Eq. (23).
- `phase_screen(r0_m, n, pixel_m, L0_m=np.inf, l0_m=1e-6, seed=None,
  subharmonics=True)` — one `n` x `n` random screen, in radians. It holds the
  modified von Karman spectrum
  `PHI(f) = 0.023 r0^(-5/3) exp(-(f/fm)^2) / (f^2 + f0^2)^(11/6)`, with
  `fm = 5.92/(2*pi*l0)` and `f0 = 1/L0`. The defaults give the pure Kolmogorov
  spectrum. `subharmonics=True` adds three low-frequency levels. `seed` is an
  INTEGER, because `aotools` builds its own generator. See Schmidt,
  DOI 10.1117/3.866274, Ch. 9.
- `Screen(Fin, phase_rad)` — a new `Field` with `E_out = E_in * exp(i*phi)`. The
  screen is a thin, pure phase element, so the power does not change. It raises
  `ValueError` on a spherical field and on a wrong-shape phase array.
- `ScreenFactory(n, pixel_m, L0_m=np.inf, l0_m=1e-6, subharmonics=True,
  n_sub_levels=3, dtype=np.float64)` — the FAST screen generator, and the
  DEFAULT of the runner (`screen_generator="olb"`). It caches the sqrt-PSD
  filter and the separable subharmonic basis ONE time for the grid, then it
  scales them for each screen by the scalar `r0^(-5/6)`. `make(r0_m, rng)` gives
  one screen; `make_stack(r0_m_array, rng)` gives a whole stack and it takes two
  independent screens from one complex FFT. It imports numpy and scipy only (the
  modified von Karman spectrum is Schmidt, DOI 10.1117/3.866274, Ch. 9,
  Eqs. (9.51) and (9.52); the Fourier-series screen is Eqs. (9.78) to (9.80);
  the subharmonics are Eq. (9.81), from Lane, Glindemann and Dainty,
  DOI 10.1088/0959-7174/2/3/003). It is 8 to 14 times faster per screen than the
  `aotools` path. It draws a DIFFERENT random atmosphere from `aotools` for the
  same seed. The broad validity pass shows the two agree in the mean collected
  power, the aperture `sigma2_I`, and the fade tail; see
  `validation/waveoptics_speed/generator_validation.py`.

**Make the screen AT the propagation pitch.** `phase_screen` takes the pitch and
the pixel count of the propagation grid. Do not make a coarse screen and
interpolate it up: a coarse screen carries no power above its own Nyquist
frequency, so it loses the structure at the Fresnel scale `sqrt(lambda*z)`, and
that structure builds the scintillation. The result then follows the coarse grid,
not the atmosphere. That route is the documented anti-pattern.

**The screen is band-limited, so the read-back is biased.** A Fourier screen
holds no power above the grid Nyquist frequency and too little power below
`1/(n*dx)`. So the measured structure function `D_phi(r)` stays BELOW the theory
`6.88 (r/r0)^(5/3)`. The three subharmonic levels lift it, but they do not close
the gap: the self-check of `screens.py` measures a deficit inside 15 percent over
`r/r0` from 0.3 to 1.6. An `r0` that is read back from `D_phi` therefore carries
that same bias. Read the bias as a RATIO between two cases, not as an absolute
value.

**The low-frequency limit is the OUTER SCALE the grid holds.** Three
subharmonic levels extend the screen to about `27 * n * dx` (a factor 3 for
each level) and no further. So a screen drawn with `L0_m=np.inf` on a 3.5 m
grid is, in its measured statistics, a von Karman screen with an outer scale
of about 95 m: the stacking test of `validation/screen_stacking/` (2026-09-04)
reads a piston-removed aperture variance of 0.765 of the Noll value on a
0.7 m aperture, and the von Karman theory at `L0 = 95 m` gives 0.765; asked
for `L0_m=25` and judged against the theory at 25 m, one screen reads
1.00 +-0.03. The tilt-removed variance is 1.00 at every `L0`. So the default
`L0_m=np.inf` CLAIMS an outer scale the grid does not deliver, and the fibre
tilt that an SMF fade pays moves with the choice by an estimated (not
measured) 2 dB at p5. Give an
explicit `L0_m`, and keep `L0_m <= 27 * n * dx` (raise `n_sub_levels` of
`ScreenFactory` when it is not; a small aperture makes a small grid). This is
backlog 2-P5.

### 9b. The split-step engine (`olb/waveoptics/turbulence/splitstep.py`)

- `super_gaussian_boundary(n, width_frac=0.125, power=8)` — an `n` x `n`
  absorbing mask, real, in the interval `[0, 1]`. It is exactly 1.0 inside the
  radius `(1 - width_frac)` of the half-side, and it falls as `exp(-t^p)` with
  `t = (rho - r_flat)/width_frac` outside it. `rho` is the radius in units of the
  half-side, so the mask is `exp(-1)` at the middle of an edge and zero at the
  corners. It raises `ValueError` when `width_frac` is outside `(0, 1]` or
  `power` is not positive.
- `split_step(Fin, z_screens_m, screens, z_total_m, *, boundary=None,
  max_step_m=None)` — one hop to each screen, the screen, and a last hop to
  `z_total_m`. Each hop uses `Forvard`. A hop longer than `max_step_m` breaks
  into EQUAL sub-steps. The default limit is `max_step = N * dx^2 / lambda`, the
  same formula as `forvard_max_z()` (the module repeats it, because `grid.py`
  reads the rest of `olb`). It returns a new `Field` at `z_total_m`. It raises
  `ValueError` on a spherical field, on unsorted distances, on a distance outside
  `[0, z_total_m]`, on a screen count that does not match the distances, and on a
  wrong-shape screen or mask.

**THE MASK IS NECESSARY. The sub-steps alone remove NO aliasing.** The sampled
transfer function of one long step is the product of the sampled transfer
functions of the sub-steps, so a split hop gives the same array as one long hop.
The mask is the part that helps: it removes the energy at the edge of the grid
between two sub-steps, before the periodic propagator brings that energy back at
the opposite edge. See Schmidt, DOI 10.1117/3.866274, Ch. 9. Give a boundary from
`super_gaussian_boundary()` for any path that spreads the beam.
`propagate_turbulent_scenario()` always does so.

**The subharmonics fight the periodic propagator.** The subharmonic content of a
screen is not periodic on the grid, and `Forvard` is periodic. So a run with
subharmonics needs the absorbing mask. The self-check of `splitstep.py` drops the
subharmonics for its plane-wave Rytov case, because that case fills the whole
grid and uses no mask: with the subharmonics on and no mask, the wrap raises the
measured variance. With the mask on, keep the subharmonics: the tilt content
drives the beam wander, and the uplink overlap of Section 9d reads that wander.

### 9c. The turbulent grid sizer (`olb/waveoptics/turbulence/sampling.py`)

The vacuum sizer `GridSpec.for_scenario()` is not sufficient for a turbulent
path. Turbulence spreads the beam, so the grid needs a wider side, and it adds
coherence structure at the Fried scale `r0`, so the grid needs a finer pixel.

#### `turbulent_grid(scenario, geometry, *, preset="standard", cn2=None, hs=None, cn2_profile=None, h_top_m=None, L0_m=np.inf)`

It returns the tuple `(GridSpec, ScreenPlan, SamplingReport)`. The geometry gives
`slant_range_m` (terrestrial) or `elevation_deg` (space), and the sizer takes the
WORST case: the longest range, or the lowest elevation.

For a SPACE link the planner is CONTINUOUS by default (item 2-I2 step 1): with no
`hs`/`cn2_profile` it builds the site Hufnagel-Valley callable, INTEGRATES it, and
cuts the atmosphere slab into equal-Rytov-weight screens at the Cn2-weighted
centroid of each slab. `cn2` overrides that callable (a `cn2(h) -> Cn2` function,
vectorised over an ndarray); `h_top_m` sets the integration top (`None` takes 20
km). The result is grid-free: a finer internal integration grid does not move the
plan. Pass an explicit `hs` (with or without `cn2_profile`) to take the LEGACY
array planner instead; `DEFAULT_HS` is now the fallback for that array caller
ONLY. `cn2`, `hs`, `cn2_profile`, `h_top_m` are all space only. `L0_m` sits here
so that one call site holds all the turbulence options; the SIZER does not read
it, and the runner passes it to `phase_screen()`.

It raises `ValueError` on an unknown preset name, and on a terrestrial transmit
terminal with no `Transmitter`. It WARNS on a sampling problem. It does not
raise, because an honest warning is better than a silent bad answer. The report
holds the same texts and the ACHIEVED numbers.

**THE EXTENT RULE.** The grid holds the beam AND the light that the turbulence
scatters out of it:

    side = [guard * 2 * r_beam + 2 * (lambda / r0_total) * z] / (1 - b)

The first part is the vacuum extent rule of `GridSpec.for_scenario()`. The second
part is the scattering cone: turbulence scatters light through the angle
`lambda/r0`, and that light must stay off the edge of the periodic grid. The
divisor `(1 - b)` makes room for the absorbing band, where `b` is
`boundary_width_frac`. See Schmidt, DOI 10.1117/3.866274, Ch. 9.

**THE PIXEL RULE.** The pixel obeys three limits, and the smallest wins:

    dx <= r0_total / pixels_per_r0    the coherence structure
    dx <= sqrt(lambda z_i) / 2        the Fresnel scale of screen i
    dx <= feature / (P / 2)           the hard edges, P = PIXELS_PER_FEATURE

The first limit comes from Schmidt, DOI 10.1117/3.866274, Ch. 9, and from Martin
and Flatte, DOI 10.1364/AO.27.002111. The second limit keeps the irradiance
correlation width sampled; the width is the Fresnel scale of the distance from
the screen to the receiver (Andrews and Phillips, DOI 10.1117/3.626196, Ch. 8).
Only a screen that carries more than `fresnel_weight_min` of the total Rytov
variance must obey it. A weak screen close to the receiver is exempt, because it
adds almost no scintillation.

**THE PIXEL COUNT.** `n` is the next power of two of `side/dx`, inside the
interval `[256, n_max]`. A clamp does NOT shrink the side, because the extent is
physics. The pixel grows instead, and the report says so.

Two module constants set the rest:

- `PIXELS_PER_FEATURE = 8` — the pixels across the smallest hard edge. The vacuum
  sizer asks for 16, but a turbulent grid is much wider, so 16 pixels across a
  small aperture makes an impossible pixel count. Pass a manual `GridSpec` to the
  runner for a finer edge.
- `MAX_SCREENS = 500` — the largest screen count that the planner builds. A path
  that asks for more gets the cap and a warning.

#### `QualityPreset` and `PRESETS`

A frozen dataclass, one named set of sampling rules. `PRESETS` maps the name to
the preset. `turbulent_grid()` and `propagate_turbulent_scenario()` both take a
name or a `QualityPreset`.

| Field | `reference` | `standard` (default) | `rapid` | Meaning |
|---|---|---|---|---|
| `pixels_per_r0` | 4 | 3 | 2 | `dx <= r0_total / pixels_per_r0`. Martin and Flatte, DOI 10.1364/AO.27.002111. Schmidt, DOI 10.1117/3.866274, Sec. 9.4, printed p. 172, gives the same rule from Johnston and Lane, and with Eq. (9.44) it reads 3.01 pixels per r0. So `standard` lands on the book value. |
| `guard` | 4 | 3 | 2 | The grid half-side over the beam radius. The same meaning as the guard of `GridSpec.for_scenario`. |
| `n_max` | 4096 | 2048 | 1024 | The largest pixel count. |
| `sigma2_r_screen_max` | 0.05 | 0.10 | 0.25 | The largest plane-wave Rytov contribution of ONE screen. A stronger screen breaks the thin-screen approximation. The book cap is `rmax = 0.1` on the LOG-AMPLITUDE variance (Schmidt, DOI 10.1117/3.866274, Listing 9.5, printed p. 175), and `sigma_R^2 = 4 sigma_chi^2`, so the book cap is 0.4 on this field. The three presets are 8x / 4x / 1.6x stricter than the book. |
| `min_screens` | 15 | 9 | 5 | The smallest screen count. `_merge_layers` clamps a weak path UP to exactly this count, so the count follows the PRESET and not the layer count of the `Cn2` profile. THE SOURCE IS olb, NOT THE BOOK: Schmidt gives no screen-count floor, and these integers come from an olb convergence sweep. The aperture scintillation index of a 30 degree downlink slab is 19 percent low at 3 screens, 10 percent low at 5, and flat from 7 up. No preset may go under 4, the moment floor of Eq. (9.65), printed p. 164. See WP7 in [schmidt-crosscheck.md](schmidt-crosscheck.md). |
| `fresnel_weight_min` | 0.005 | 0.02 | 0.05 | The Rytov share above which a screen must obey the Fresnel-scale pixel rule. The exemption is an olb rule; Schmidt, Sec. 9.4, applies the rule to every step. |
| `boundary_width_frac` | 0.125 | 0.125 | 0.10 | The width of the absorbing band, as a fraction of the half-side. It goes to `super_gaussian_boundary()`. |

#### `ScreenPlan`

A frozen dataclass. Where the screens sit, and what each one carries. THE PLAN IS
A LIST OF SLABS.

| Field | Type | Meaning |
|---|---|---|
| `z_m` | array | The distance of each screen from the INPUT plane, in m. The values go up. |
| `cn2_int_m13` | array | The integrated `Cn2` of each screen, in m^(1/3). It carries the slant factor already. |
| `r0_m` | array | The Fried parameter of each screen, in m. |
| `sigma2_r` | array | The plane-wave Rytov contribution of each screen. |
| `z_total_m` | float | The length of the gridded path, in m. |
| `r0_total_m` | float | The composite Fried parameter of the whole path, `(SUM r0_i^(-5/3))^(-3/5)`. |
| `direction` | str | `"terrestrial"` or `"down"`. A space plan always propagates DOWN. |

The per-screen Rytov contribution is
`d(sigma_R^2) = 2.25 k^(7/6) (INT Cn2 dz) (z_to_rx)^(5/6)`. A constant `Cn2`
gives the familiar `1.23 Cn2 k^(7/6) L^(11/6)`, because `2.25 * (6/11) = 1.23`.
See Andrews and Phillips, DOI 10.1117/3.626196, Ch. 8, Eq. (20), and Ch. 12,
Eqs. (36) and (38).

Where the boundaries go:

- **Terrestrial.** The path is uniform, so the screens share it EQUALLY and each
  screen sits at the centre of its slab. The planner starts from the mean-share
  estimate `sigma2_total / sigma2_r_screen_max` and raises the count until the
  STRONGEST screen (the one farthest from the receiver) obeys the cap.
- **Space.** The layers of the `Cn2` profile come from
  `olb.turbulence.profiles.default_cn2_profile`, times the airmass `sec(zenith)`
  (Andrews and Phillips, DOI 10.1117/3.626196, Ch. 12, Eq. (14)). The planner
  MERGES adjacent layers until each group stays under the Rytov cap, and it puts
  the screen at the `Cn2`-weighted centre of the group. It only merges: it does
  not split one layer, because a profile gives no sub-layer structure. A single
  layer that is stronger than the cap keeps its own screen, and the sizer warns
  to ask for a finer `hs` grid.

#### `SamplingReport`

A frozen dataclass. What the grid ACHIEVES, against what the preset asks for.

| Field | Type | Meaning |
|---|---|---|
| `pixels_per_r0` | float | The achieved `r0_total / dx`. |
| `grid_margin` | float | The untouched interior half-side, divided by the beam radius plus the scattering cone. 1.0 means the light just fits. The preset `guard` is the target. |
| `fresnel_pixels_min` | float | The smallest achieved pixel count across a REQUIRED Fresnel scale `sqrt(lambda z)`. It is infinite when no screen passes `fresnel_weight_min`. 2 or more is good. |
| `step_over_limit_max` | float | The largest planned gap between two screens, divided by `forvard_max_z()`. 1.0 or less is good. The engine cuts a longer gap into sub-steps. |
| `sigma2_r_screen_max` | float | The largest per-screen Rytov contribution that the plan holds. |
| `n_clamped` | bool | True means the pixel count hit `n_max`. |
| `warnings` | tuple | The warning texts that the sizer sent. |

The sizer warns when the pixel count hits `n_max`, when the achieved
`pixels_per_r0` is below the preset value, when `fresnel_pixels_min` is below
2.0, when the strongest screen passes `sigma2_r_screen_max`, and when the plan
hits `MAX_SCREENS`.

The report gives the achieved numbers of ONE grid. The layer runs NO automatic
convergence check. To prove a case, run it again on a finer preset, or on a wider
manual grid, and compare. Martin and Flatte validated the method that way,
DOI 10.1364/JOSAA.7.000838. The example scripts of Section 9f do exactly this by
hand.

### 9d. The trial runner (`olb/waveoptics/turbulence/run.py`)

#### `propagate_turbulent_scenario(scenario, geometry, *, n_trials=1, seed=None, preset="standard", grid=None, plan=None, cn2=None, hs=None, cn2_profile=None, h_top_m=None, L0_m=np.inf, subharmonics=True, threader=None, screen_generator="olb", progress=False, detectors=None, start_index=0, patch_radius_m=None, precision="single")`

It runs a set of turbulent split-step trials for one scenario and it returns a
`TurbWaveResult`. Each trial makes a NEW screen stack and moves one field through
it. The trials are independent snapshots.

- `grid` and `plan` come together, or neither comes. `None` for both calls
  `turbulent_grid()`, and the result then carries the report.
- `cn2` (a continuous `cn2(h)` callable), `hs`/`cn2_profile` (the legacy array),
  and `h_top_m` are the space profile options; they pass through to
  `turbulent_grid()`. `None` everywhere takes the continuous site profile. See
  `turbulent_grid()` above.
- The geometry must give ONE range. More than one raises `ValueError`. Loop in
  the caller.
- A `"retro"` direction raises `NotImplementedError`.
- `subharmonics=True` is the value to keep: the tilt content drives the beam
  wander, and the uplink overlap reads that wander.
- `screen_generator` is `"olb"` (the default, the fast `ScreenFactory`) or
  `"aotools"` (the reference path). The two give DIFFERENT draws for the same
  seed; the statistics agree. Only `"aotools"` needs the `aotools` package. An
  unknown name raises `ValueError`. `propagate_turbulent_field()` takes the same
  argument, with the same default.
- `precision` is `"single"` (the DEFAULT since 2026-09-05: a complex64 field
  and float32 screens from `ScreenFactory(dtype=np.float32)`, or an aotools
  screen cast to float32) or `"double"` (a complex128 field and float64
  screens). Another value raises `ValueError`. WHY: a campaign on a many-core machine is
  memory-bandwidth bound, and single precision halves the bytes each FFT
  moves. Measured (2026-09-05, `validation/campaign_resources/`): 12 workers
  on a 512 px grid give 11.2 trials/s single against 8.5 double, 1.32x, with
  FEWER busy threads. The physics agrees to parts per million: the collected
  power to 6e-7, the SMF eta to 3e-6 and the receive field rms to 2e-6
  relative, at the rapid and the standard presets
  (`validation/precision/`). CAUTION: a single-precision run is a DIFFERENT
  record. Its trials are not bit-identical to a double run of the same seed.
  `precision="double"` reproduces every run made before 2026-09-05 bit for
  bit; the studies that reopen those stored campaigns pass it explicitly. The stored patch is
  complex64 in both modes. `recouple()` and `recollect()` rebuild the grid in
  complex128 whatever the mode. `propagate_turbulent_field()` takes the same
  argument, with the same default.
- `threader` is an optional `olb.waveoptics.Threader`. `None` runs the trials one
  by one. A `Threader` runs them across threads and it keeps the trial order; the
  FFT releases the GIL, so it gives a real speed-up. `Threader()` with no
  argument takes `min(16, cores)` workers: the scaling study
  (`validation/waveoptics_speed/scaling_study.py`) finds the thread rate
  saturates at 8 to 16 workers.
- `progress=True` shows a tqdm bar that advances one step for each finished
  trial. It needs the optional `tqdm` package; without `tqdm` the run goes on
  with no bar and a warning. `False` (the default) shows no bar. With a threader
  the bar advances in the finishing order, and the returned trials keep the trial
  order.
- `detectors` is an optional sequence of detector objects, the arms behind a
  receive beamsplitter. See the paragraph below. `None` (the default) keeps the
  single-detector record, bit for bit.
- `start_index` is the index of the FIRST trial. The run covers the trials
  `start_index .. start_index + n_trials - 1`, and each `TurbTrial.seed_key`
  holds the TRUE index. `0` (the default) is the old behaviour.
- `patch_radius_m` stores the receive-plane field on a disc of that radius, in
  m. `None` (the default) stores no field, and the record is bit for bit the old
  record. A float fills `TurbWaveResult.fields` and `TurbWaveResult.patch`. See
  the paragraph below.

**THE BLOCK CONTRACT (`start_index`).** The runner seeds trial `k` off
`(seed, k)`, so a block of trials is a SLICE of one long run. A run of `n`
trials therefore equals the concatenation of its blocks, trial for trial and bit
for bit: `(start_index=0, n_trials=200)` plus `(start_index=200, n_trials=300)`
gives exactly the trials of `(start_index=0, n_trials=500)`. So a campaign
computes its blocks in any order, and on any number of processes. Section 9g
builds on that contract.

**THE STORED FIELD (`patch_radius_m`).** The runner writes the UNCLIPPED
receive-plane field at the pixels inside the radius, one row for each trial. The
disc uses the SAME pixel-centre convention as
`olb.waveoptics.sources.CircAperture`, so a patch of the radius `D/2` holds
exactly the pixels that a `CircAperture` of the diameter `D` keeps. The rows are
`complex64`. A radius larger than half the grid side raises `ValueError`. The
memory is small: a 1 m patch at a 5 mm pixel pitch is about 200 x 200
`complex64`, which is about 320 KB for each trial. The scalars do not change.
Read a stored field back with `recouple()` and `recollect()` below.

**ONE RUN, MANY ARMS (`detectors`).** With `detectors`, each trial computes the
coupling efficiency of EVERY arm on the SAME clipped receive field, and it
records them in `TurbTrial.detector_etas`, in the argument order. So N arms cost
ONE run, not N: the field is already in memory, and each arm is one more cheap
focal-plane calculation on that same array. The `frac` of a detector is IGNORED
here. A beamsplitter scales the field of an arm by a constant, and every coupling
efficiency is power-normalised, so the split ratio does not change it; the
fraction is a separate fixed dB Term (`olb.models.splitter`). An `SMF` arm gives
the mode overlap of Section 4, an `MMF` arm gives the core capture of Section 4a
(with its `defocus_m`), an `Aperture` arm gives exactly `1.0` (the aperture
capture is already in `collected_power`), and a `Camera` arm gives `None`,
because a `Camera` has no coupling model. An unknown detector type raises
`ValueError`.

**THE BOUNDARY MASK IS ALWAYS ON.** The runner builds
`super_gaussian_boundary(grid.n, preset.boundary_width_frac)` and it gives that
mask to every hop and every screen. The sizer keeps the receive aperture inside
the untouched interior of the mask, and the runner CHECKS that: a receive
aperture radius at or past `(1 - boundary_width_frac) * size_m / 2` warns that
the collected power is too low.

**TWO CASES.**

- **TERRESTRIAL.** The runner launches the transmit beam of the near terminal
  with the launch recipe of `propagate_scenario()` (a Gaussian at the virtual
  waist, an offset propagation, and the launch-aperture clip), and it propagates
  that beam along the horizontal path. It imports those helpers, so the vacuum
  limit of this module IS the vacuum module.
- **SPACE.** The gridded path is the DOWNLINK atmosphere slab ONLY. The satellite
  sits outside the atmosphere, so a unit PLANE WAVE enters at the top of the
  slab, and the vacuum above the slab carries no turbulence. **The uplink
  direction is never propagated.** A downlink reads the collected power at the
  ground. An uplink reads the SAME field through reciprocity.

**THE VACUUM BASELINE (space).** A unit plane wave fills the grid, so the
absorbing mask acts as a soft aperture, and over a long slab that soft edge makes
strong Fresnel rings on the axis. Those rings are a property of the GRID, not of
the atmosphere. So the space reference is the SAME plane wave, along the SAME
hops, through the SAME mask, with FLAT screens. Then the vacuum limit of each
space output is exactly 1.0, and each number is a pure turbulence penalty.

**THE RECIPROCITY OVERLAP (uplink).** The turbulent atmosphere is reciprocal, so
the uplink flux at the satellite is the overlap of the received downlink field
with the ground transmit mode. See Shapiro, DOI 10.1364/JOSA.61.000492. The
transmit mode `psi_tx` is the launch recipe above, scaled so that
`sum(|psi_tx|^2) = 1.0`. Then

    eta_turb = |SUM E_rx conj(psi_tx)|^2 / |SUM E_vac conj(psi_tx)|^2

where `E_vac` is the zero-screen vacuum run through the same mask and the same
hops. So the vacuum limit of `eta_turb` is exactly 1.0, and
`-10*log10(eta_turb)` is the uplink turbulence loss on the free-space baseline.
That is the baseline of the `(w_free/w_st)^2` rescale of
`olb.turbulence.uplink_flux`, so the two numbers compare. Point-ahead
anisoplanatism is NOT modelled: the uplink and the downlink read the same
screens.

**THE SEED CONTRACT.** `seed` takes an int, a numpy `Generator`, or `None` for a
fresh entropy, and the runner resolves it to ONE integer. Each screen then draws
from `SeedSequence(entropy, spawn_key=(trial, screen))`. So trial `k` is
bit-identical for one entropy AND one `screen_generator`, and the trial count
does not change it: a longer run repeats the trials of a shorter run.
`TurbWaveResult.seed_entropy` gives the integer back, and each trial records its
own `seed_key`. A change of `screen_generator` keeps the same integer seeds but
gives a different screen, so the two generators are not bit-identical.

#### `TurbTrial`

A frozen dataclass. One atmosphere snapshot.

| Field | Type | Meaning |
|---|---|---|
| `collected_power` | float | The power inside the receive aperture, as a fraction of the input power. The terrestrial case divides by the launched power AFTER the transmit clip, so it holds the geometric spread too. The space case divides by the VACUUM baseline on the same grid, so it holds the turbulence penalty only, and its vacuum limit is 1.0. |
| `smf_eta` | float or None | The single-mode-fibre coupling efficiency, from `olb.waveoptics.smf.coupling_efficiency`. `None` when the receive terminal has no `SMF` detector. |
| `eta_turb` | float or None | The uplink reciprocity overlap ratio, against the free-space baseline. `None` for a downlink and for a terrestrial case. |
| `mmf_eta` | float or None | The multimode-fibre (light-bucket) coupling efficiency, from `olb.waveoptics.mmf`. It is the encircled energy of the focused spot inside the core; the turbulent tilt walks the spot off the core. It also holds the NON-FOCAL-PLANE detector: the runner reads `MMF.defocus_m` for the plane `z = f + defocus_m`. At the focal plane (`defocus_m = 0`) this is the plain focal-plane coupling. `None` when the receive terminal has no `MMF` detector. |
| `seed_key` | tuple | The pair `(seed_entropy, trial_index)`. |
| `wall_time_s` | float | The time of the trial, in s. It holds the screen generation and the propagation. |
| `detector_etas` | tuple or None | The coupling efficiency of each detector of the `detectors` argument, in that order. `None` on the default path, so a single-detector record does not change. A `Camera` arm holds `None`, and an `Aperture` arm holds `1.0`. |

#### `TurbWaveResult`

A frozen dataclass. The result of a set of trials.

| Field | Type | Meaning |
|---|---|---|
| `trials` | list | One `TurbTrial` for each trial. |
| `grid` | `GridSpec` | The grid that the trials used. |
| `plan` | `ScreenPlan` | The screen plan that the trials used. |
| `report` | `SamplingReport` or None | `None` when the caller gives its own grid and plan. |
| `preset` | str | The name of the quality preset. |
| `seed_entropy` | int | The integer that seeds every trial. Give it back to repeat the set. |
| `fields` | `np.ndarray` or None | The stored receive-plane field on the patch, a `complex64` array of the shape `(n_trials, n_patch)`. The row order is the trial order. `None` when the caller asks for no patch. |
| `patch` | `FieldPatch` or None | The `FieldPatch` of those columns, or `None`. |

**The record holds the per-trial SCALARS, and, when the caller asks for it, the
OPTIONAL masked receive field (`fields` and `patch`).** The field capture is an
owner decision of 2026-09-04: a large campaign must recouple a stored field to a
new detector, without a new propagation. The scalars stay exactly as they are,
and a budget never reads the fields.

#### `FieldPatch`

A frozen dataclass. The mask that selects the stored receive-plane pixels. The
patch is a disc at the centre of the grid, in the pixel-centre convention of
`olb.waveoptics.sources.CircAperture`.

| Field | Type | Meaning |
|---|---|---|
| `radius_m` | float | The radius of the disc, in m. |
| `n` | int | The number of grid pixels along one side. |
| `pixel_m` | float | The distance between two pixels, in m. |
| `indices` | `np.ndarray` | The flat indices of the disc pixels, an `int32` array. The order is the C order of the `n` x `n` grid. |

#### `recouple(result, detector, aperture_m, obscuration_ratio, lam, *, trials=None)`

It couples a STORED receive field into a detector, after the run. The function
rebuilds the receive-plane field of each stored trial, it clips that field at
the receive aperture, and it gives the coupling efficiency of the detector. So a
campaign tries a NEW detector, a new focal length or a new defocus with NO new
propagation. The physics is the physics of the run: the function calls the same
internal helper on the same clipped field.

The reconstruction keeps the FULL grid. It scatters the stored patch values into
a zero array of the whole grid, because a crop would change the zero padding,
and the focal-plane pixel scale of a fibre coupling reads that padding. So the
value equals the in-run value. The rebuild gives ONE grid at a time, so the
memory holds one field only.

- `detector` is an `SMF`, an `MMF`, an `Aperture`, a `Camera`, or `None`. A
  detector with no coupling model (a `Camera`, or `None`) gives `NaN`.
- `trials` is an optional sequence of trial row indices. `None` takes every
  stored trial.
- It returns a float array, one value for each selected trial.
- A result that holds no field raises `ValueError`. A receive aperture larger
  than the stored patch also raises `ValueError`.

#### `recollect(result, aperture_m, obscuration_ratio, *, trials=None)`

It gives the collected power of each STORED trial, in grid units. The value is
the power of the clipped rebuilt field, and it is NOT normalised: the runner
divides its `collected_power` by a vacuum reference, and this function does not
know that reference. So divide by your OWN reference, or take the RATIO of two
trials, which needs no reference. The rebuild is the full-grid rebuild of
`recouple()` above.

- `trials` is an optional sequence of trial row indices. `None` takes every
  stored trial.
- It returns a float array, one value for each selected trial.
- A result that holds no field raises `ValueError`. A receive aperture larger
  than the stored patch also raises `ValueError`.

### 9e. The stubs

Each of these raises. Each one is a deliberate deferral, not a defect.

| Name | Where | Why it is deferred |
|---|---|---|
| `TemporalScreens` | `temporal.py` | The frozen-flow time axis. The constructor and `step()` raise `NotImplementedError`. The class docstring holds the recorded design: one `aotools` `PhaseScreenVonKarman` for each layer, `add_row()` to extrude, and a layer drift velocity that is the vector sum of the Bufton wind (Andrews and Phillips, DOI 10.1117/3.626196, Ch. 12, Eqs. (2) and (3)) and the apparent translation `omega_slew * z_i` of a tracked satellite. See also Assemat, Wilson and Gendron, DOI 10.1364/OE.14.000988, and Taylor, DOI 10.1098/rspa.1938.0032. |
| `folded_terrestrial()` | `run.py` | The double pass of a corner-cube retroreflector. The two passes share the same screens, so they are correlated. That correlation is the physics of the link, and it needs its own design. |
| The `"retro"` direction | `run.py` | `propagate_turbulent_scenario()` raises `NotImplementedError`. The same correlated double pass. |
| A co-moving screen | `screens.py`, `splitstep.py` | `Screen()` and `split_step()` raise `ValueError` on a spherical field. The split step runs on a FLAT grid only. Call `Convert()` first. |

### 9f. The example scripts

The suite in `examples/waveoptics/` holds TEN scripts: three vacuum scripts
(`space_farfield.py`, `terrestrial_stages.py`, `grid_artefacts.py`), the three
turbulent scripts below, the budget demonstration `budget_wiring.py`, two
multimode-fibre demonstrations (`mmf_core_psf.py`,
`mmf_core_psf_terrestrial.py`), and the camera demonstration
`camera_tracking.py`.

These three put the turbulent layer against the analytic models that the budgets
already use. Each one runs for about four to five minutes.

- `turbulent_terrestrial.py` — a 2 km horizontal link at `Cn2 = 3e-15`, three
  receive apertures on the same screens and the same seeds.
- `turbulent_downlink.py` — a 600 km downlink into a 500 mm obscured fibre
  receiver, at 30, 60 and 90 degrees.
- `turbulent_uplink_reciprocity.py` — a 600 km uplink through the overlap of
  Section 9d, at the zenith and at 30 degrees.

For the guide to each of the ten scripts, see
[examples/waveoptics/README.md](../examples/waveoptics/README.md). See also
[examples.md](examples.md) for what each one prints and what it shows.

### 9g. The campaign store (`olb/waveoptics/turbulence/campaign.py`)

A fade statistic needs thousands of snapshots. One trial is expensive, so a
campaign must survive a stopped process, grow later, and give its fields back to
a NEW detector with no new propagation. `Campaign` is that store: it keeps the
trials on disk in fixed BLOCKS, it is resumable, and it runs its blocks on ONE
warm process pool.

Import it from the sub-package:

```python
from olb.waveoptics.turbulence import Campaign
```

#### `Campaign(scenario, geometry, root_dir, *, seed, preset="standard", block_size=100, patch_radius_m=None, sizing_aperture_m=None, grid=None, plan=None, cn2=None, hs=None, cn2_profile=None, h_top_m=None, L0_m=np.inf, subharmonics=True, screen_generator="olb", precision="single")`

It opens a campaign, or it makes a new one. A `Campaign` names ONE physics case:
one scenario, one geometry, one grid, one screen plan, one seed.

- `seed` is REQUIRED, and it must be an integer: a campaign grows over more than
  one session, so its trials must repeat. `None` or a numpy `Generator` raises
  `ValueError`.
- `block_size` is the number of trials in one block. Block `b` holds the trials
  `b*block_size .. (b+1)*block_size - 1` of ONE native run.
- `patch_radius_m` is the radius of the stored field disc, in m. `None` takes
  `sizing_aperture_m / 2` when a sizing aperture is given, else the receive
  aperture / 2.
- `sizing_aperture_m` is an optional LARGER receive aperture that sizes the
  grid. See the rule below.
- `grid` is an optional `GridSpec`. The plan still comes from the `Cn2` inputs.
- `plan` is an optional `ScreenPlan`. Give it WITH `grid` to hold the grid
  fixed and move the screens only (a convergence study; the sizer moves the
  grid with the screen count, so a naive `min_screens` sweep moves two
  things). Both enter the fingerprint, so a different plan is a different
  campaign, and a reopen that drops the plan raises. The sizer still runs when
  one of the two is `None`, and the given one overrides its half. See
  `validation/tail_convergence/`.
- `cn2`, `hs`, `cn2_profile`, `h_top_m`, `L0_m`, `subharmonics` and
  `screen_generator` pass to the sizer and the runner of Section 9d, with the
  same meanings and the same defaults.
- `precision` passes to the runner (Section 9d). It enters the manifest and
  the fingerprint ONLY when it is `"single"`, so every key and every manifest
  of a campaign stored before 2026-09-05 (all double) stays valid; a manifest
  with no `precision` key reads as `"double"`. A single-precision campaign is
  its own store, and a reopen with the other precision raises, so pass
  `precision="double"` to reopen an old store.

Attributes: `root_dir`, `scenario`, `geometry`, `seed`, `preset`, `block_size`,
`patch_radius_m`, `grid`, `plan`, `patch`, `fingerprint`.

**The fingerprint.** `fingerprint` is `cache_key(...)` from
`olb/waveoptics/turbulence/fingerprint.py`: one SHA-256 of everything that
changes a trial (the scenario repr, a canonical geometry signature, the preset,
the seed, the screen generator, the outer scale, the subharmonic switch, the
Cn2 inputs, the block size, a caller grid and plan, and `KEY_VERSION`). The manifest
stores it, and an existing campaign whose fingerprint does not match raises.
This key came from the P4 scalar cache (`cache.py`), which `Campaign` replaced
and which was RETIRED on 2026-09-04; the value of the key did not change, so an
existing manifest still matches.

#### `Campaign.run(n_trials, *, workers=None, progress=False)`

It computes and stores the MISSING blocks up to `n_trials` trials, and it
returns the number of trials on disk. The call rounds `n_trials` up to a whole
number of blocks. A block that already sits on disk is NOT recomputed.
`progress=True` prints one line for each finished block.

#### `Campaign.load(n_trials=None, *, fields=True)`

It assembles a `TurbWaveResult` from the stored blocks. The record is the record
of a native run: the trial order is the trial order, and `seed_key` holds the
TRUE trial index, so
`olb.models.waveoptics.waveoptics_turbulence_term` reads it unchanged.
`n_trials=None` takes every stored trial. `fields=False` leaves
`TurbWaveResult.fields` and `TurbWaveResult.patch` `None`, so a budget-only load
stays small.

#### `Campaign.recouple(detector, aperture_m=None, obscuration_ratio=None, n_trials=None)`

It couples the STORED fields into a detector, with no new propagation, and it
returns a float array of the coupling efficiency of each trial. `aperture_m` and
`obscuration_ratio` of `None` take the values of the scenario receive terminal.
The call STREAMS: it reads one block at a time, so ten thousand trials never sit
in RAM at the same time. The physics is `recouple()` of Section 9d.

#### `Campaign.recollect(aperture_m=None, obscuration_ratio=None, n_trials=None)`

It gives the collected power of each STORED trial, in grid units, as a float
array. The value is NOT normalised: it holds no vacuum reference, so take the
RATIO of two trials, or divide by your own reference. It streams the same way.

#### `Campaign.n_stored`

A property. The number of trials on disk, counted from block 0 with no gap.

#### The file layout

| File | What it holds |
|---|---|
| `block_{b:05d}.npz` | One block. The five per-trial scalar columns (`collected_power`, `smf_eta`, `mmf_eta`, `eta_turb`, `wall_time_s`; `NaN` marks a `None`) and the `complex64` `fields` rows. |
| `manifest.json` | The fingerprint, the seed, the preset, the block size, the patch radius, the sizing aperture, the screen generator, the outer scale, the subharmonics flag, the olb version, the scenario text, the grid, the plan and the patch shape. |
| `patch_indices.npy` | The flat pixel indices of the `FieldPatch`. |

A block file holds ONE block, and the parent writes it with an atomic replace.
So a stopped campaign keeps every finished block.

**The manifest rebuild.** An EXISTING `root_dir` is checked, and the grid and
the plan then come from the manifest, NOT from a new sizing call. So a resumed
campaign NEVER re-sizes, and the atmosphere of a new block is the atmosphere of
the old blocks.

**The mismatch rule.** The fingerprint, the seed, the preset, the block size,
the patch radius and the sizing aperture must match the stored manifest. A
different value raises `ValueError`, and the message names the field. A stored
campaign is ONE physics case: use a new directory, or match the stored settings.

#### `workers`: ONE level of parallelism

- `workers=None` (the default) runs the blocks ONE AFTER THE OTHER in this
  process, and each block threads inside with the default `Threader` of
  Section 9d.
- `workers=W` opens ONE `ProcessPoolExecutor` of `W` processes for the WHOLE
  `run` call. A module-level initializer fills the worker state ONE time for
  each process, so the scenario, the geometry, the `GridSpec` and the
  `ScreenPlan` cross the process boundary once, not once for each block. A block
  then runs SERIALLY inside its process.

Never both: threads inside processes over-subscribe the cores. The parent writes
each block file as soon as that block arrives, so a killed campaign keeps every
finished block.

**The measured facts behind the rule**
(`validation/waveoptics_speed/fair_scaling_rerun.py`, 2026-09-04). Threads and
processes TIE on wall time for ONE 200-trial run, because a Windows process pool
costs 2.5 to 4.4 s to spawn. Processes beat threads by 1.15x to 1.7x in steady
state. So a process pool pays ONLY when it stays WARM across many blocks, which
is exactly the campaign case.

**PICKLING.** `workers=W` sends the scenario, the geometry, the `GridSpec` and
the `ScreenPlan` to each process. Those are dataclasses, so they pickle. The
`cn2` callable is NOT sent: the parent plans the screens one time and the
workers get the finished plan, so a lambda `cn2` is safe here. The scenario and
the geometry must still be picklable objects at module level.

#### `sizing_aperture_m`: size the grid one time for a family

Each trial stores the receive-plane field BEFORE the receive-aperture clip.
Store that field at the LARGEST receive aperture of the family ONE time. Then a
smaller receive aperture, a central obscuration, a different detector, a
different focal length and a different defocus are all a POST-HOC crop of the
stored field, through `recouple()` and `recollect()`, with NO new propagation.

`sizing_aperture_m` serves that plan: the grid is sized on a COPY of the
scenario whose RECEIVE terminal carries the larger aperture, and the trials then
run on THAT grid with the ORIGINAL scenario. So one campaign covers every
receive aperture up to the sizing aperture.

This is EXACT for a SPACE downlink, because the propagated slab does not read
the receive terminal at all: the input is a plane wave, and the receive terminal
enters only at the clip. The SPACE uplink reciprocity overlap reads the SAME
field. A TERRESTRIAL path does read the TRANSMIT terminal, so only the receive
side is free: a changed launch, or a changed path, needs a rerun.

#### Budgets from a campaign

A `Campaign` IS a fidelity-2 wave record, so it goes straight into the `wave`
slot of a budget:

```python
from olb import multi_detector_budgets
from olb.links import downlink_budget

campaign = Campaign(scenario, orbit, root, seed=2024).run(1000, workers=4)

# one budget of the scenario receive path
budget = downlink_budget(scenario, orbit, fidelity=2, wave=campaign)

# one budget for each beamsplitter arm, from the SAME campaign
arms = [SMF(frac=0.7), MMF(core_radius_m=25e-6, focal_length_m=0.5)]
budgets = multi_detector_budgets(scenario, orbit, arms, fidelity=2,
                                 wave=campaign)
```

`downlink_budget`, `uplink_budget`, `terrestrial_budget` and
`multi_detector_budgets` all accept it. Each one calls
`olb.models.waveoptics.resolve_wave(wave, detectors=None)`, the ONE adapter: it
returns a `Fidelity2Bundle` or a list of them unchanged, and it turns a
`Campaign` into the record the caller needs.

The two EXPLICIT forms do the same work, and a budget takes their output the
same way:

- `campaign_bundle(campaign, *, n_trials=None, vacuum=None)` gives ONE
  `Fidelity2Bundle`. The stored trials already carry the `smf_eta` or the
  `mmf_eta` of the scenario detector, so it reads the scalars only: no field is
  loaded and no re-coupling runs.
- `campaign_bundles(campaign, detectors, *, n_trials=None, vacuum=None)` gives a
  LIST, one bundle for each arm, in the `detectors` order. It is the POST-HOC
  twin of `run_fidelity2(..., detectors=[...])`: the campaign is the shared
  Monte Carlo, and each arm re-couples the SAME stored fields into its own
  detector, with no new propagation. The `frac` of a detector is IGNORED here;
  it enters one time, as the fixed splitter Term of `multi_detector_budgets`.

`vacuum` is the selector of `run_fidelity2`: `None` takes "analytic" for a space
link and "wave" for a terrestrial link.

`recouple` and `recollect` are DIAGNOSTIC, in the way a `Camera` is diagnostic.
They answer a "what if" question about the stored field (a smaller aperture, a
different focal length, a detector that was never on the terminal). They are NOT
the budget path: a budget reads the campaign itself.

A fidelity-2 budget takes ONE line of sight, so `geometry` must give one
elevation. A one-element elevation array is accepted (it IS one line of sight);
a true multi-element array raises. See
`examples/waveoptics/campaign_demo.py`.

#### Why the store keeps fields, and not screens

A phase screen is NOT stored. The seed regenerates a screen bit-identically in
tens of milliseconds, and ten thousand trials of screens would take 200 to
300 GB. The masked receive field is much smaller:

| Quantity | Value |
|---|---|
| The stored type | `complex64` |
| A 1 m patch at a 5 mm pixel pitch | about 200 x 200 pixels |
| The disk of one trial | about 320 KB |
| The disk of 10,000 trials | about 3.2 GB |

---

## 10. The Schmidt foundation layer (`olb/waveoptics/schmidt/`)

The sub-package holds the numerical method of Schmidt (2010),
DOI 10.1117/3.866274, as pure book physics. It imports numpy and scipy only, it
imports nothing from the rest of `olb`, and it returns no decibels. Each
function names its chapter, its equation number and its printed page.

Status: the layer is VALIDATION ONLY. No budget, no Term, no sizer and no runner
consumes it, by owner decision. `olb/waveoptics/schmidt/__init__.py` exports
nothing, so import each module by name:

```python
from olb.waveoptics.schmidt.fresnel import angular_spectrum
from olb.waveoptics.schmidt.sampling import check_sampling
```

The physics is in [physics.md](physics.md) Section 8. The equation-by-equation
map, the gaps and the constants ledger are in
[schmidt-crosscheck.md](schmidt-crosscheck.md).

### 10a. The transforms (`fourier.py`)

| Name | What it gives |
|---|---|
| `freq_pitch(n, dx)` | The frequency pitch `df = 1/(N dx)`. Ch. 2, text below Eq. (2.3); Ch. 6, Eq. (6.51). |
| `ft2(g, dx)` | The scaled two-dimensional forward transform. Ch. 2, Eq. (2.6), with the 2-D scaling of Sec. 2.6. |
| `ift2(G, df)` | The scaled two-dimensional inverse transform. Ch. 2, Eq. (2.9). |
| `structure_function(ph, mask, dx)` | The structure function of a masked phase screen, by transform. Ch. 3, Eqs. (3.15) to (3.25). |

### 10b. The propagation kernels (`fresnel.py`)

Every kernel DROPS the piston factor `exp(i k z)`, as the book listings do. The
production `Forvard` keeps it. Bridge that before any phase comparison.

| Name | What it gives |
|---|---|
| `one_step_fresnel(Uin, wavelength, dx1, z)` | The one-step Fresnel transform, and its FIXED output pitch `lambda z/(N dx1)`. Ch. 6, Eqs. (6.5), (6.15), (6.16). |
| `two_step_fresnel(Uin, wavelength, dx1, dx2, z)` | Two partial Fresnel integrals with a FREE output pitch. It refuses m = 1. Ch. 6, Eqs. (6.18) to (6.25). |
| `angular_spectrum(Uin, wavelength, dx, z, dx2=None)` | The angular-spectrum propagator. `dx2=None` gives the baseline form, Eqs. (6.31) and (6.32); a value gives the SCALED form, Eq. (6.65). |
| `super_gaussian_absorber(n, sigma_frac=0.47, power=16)` | The book absorbing boundary. Ch. 8, Eq. (8.1), with the values of Listing 8.1. |
| `partial_propagations(Uin, wavelength, dx1, dxn, z_planes, absorber=None)` | The general partial-propagation chain, with the linear per-plane pitch rule. Ch. 8, Eqs. (8.8), (8.14) to (8.18). |

### 10c. The sampling constraints (`sampling.py`)

Small pure functions. Nothing warns and nothing raises on a broken rule: the
module measures, and the caller acts.

| Name | What it gives |
|---|---|
| `nyquist_max_angle(wavelength_m, delta)` | The largest ray angle the grid carries. Ch. 7, Eq. (7.7). |
| `geometric_max_angle(D1, D2, delta1, delta2, z)` | The angle the geometry demands. Ch. 7, Eqs. (7.8), (7.9), (7.12). |
| `constraint1_max_delta2(D1, D2, delta1, wavelength_m, z)` | Constraint 1. Ch. 7, Eq. (7.14). |
| `illuminated_diameter(D1, delta1, delta2, wavelength_m, z)` | D_illum, the diameter that the source illuminates. Ch. 7, Eq. (7.16). |
| `constraint2_min_n(D1, D2, delta1, delta2, wavelength_m, z)` | Constraint 2. Ch. 7, Eq. (7.20). |
| `local_spatial_frequency_source(x1, wavelength_m, z, R=inf, m=None)` | The local frequency of the source quadratic phase. `m=None` gives Eq. (7.39); a value gives Eq. (7.51). |
| `local_spatial_frequency_transfer(f1, wavelength_m, z, delta1, delta2)` | The local frequency of the transfer function. Ch. 7, Eqs. (7.55), (7.57). |
| `constraint3_delta2_window(D1, delta1, wavelength_m, z, R=inf)` | Constraint 3, as a window on delta2. Ch. 7, Eq. (7.53). |
| `constraint3_is_slack(D1, D2, z, R=inf)` | True when constraint 3 does not apply. Ch. 7, Eq. (7.60). |
| `constraint4_min_n(delta1, delta2, wavelength_m, z)` | Constraint 4. Ch. 7, Eq. (7.59). |
| `one_step_delta2(N, delta1, wavelength_m, z)` | The fixed one-step output pitch. Ch. 7, Eq. (7.21). |
| `one_step_min_n(D1, D2, delta1, wavelength_m, z)` | The minimum N of a one-step run. Ch. 7, Eq. (7.25). |
| `fresnel_min_distance(D1, delta1, wavelength_m, R=inf)` | The MINIMUM one-step distance. Ch. 7, Eqs. (7.41), (7.42). |
| `two_step_planes(z, m)` | The two intermediate-plane geometries of a two-step run. Ch. 6, Eqs. (6.24) to (6.29). |
| `angular_spectrum_max_z(N, delta1, wavelength_m)` | The range limit `N dx^2/lambda`. Ch. 7, Eq. (7.59), inverted. It reproduces `olb.waveoptics.grid.forvard_max_z`. |
| `partial_grid_spacing(delta1, delta_n, alpha)` | The per-plane pitch of a partial propagation. Ch. 8, Table 8.2. |
| `partial_max_step(N, delta1, delta_n, wavelength_m)` | The step cap. Ch. 8, Eq. (8.24). |
| `partial_plane_count(z, N, delta1, delta_n, wavelength_m)` | The plane count from that cap. Ch. 8, text below Eq. (8.24). |
| `absorbing_boundary_sigma(N, frac=0.47)` | The book absorber half-width, in pixels. Ch. 8, Listing 8.1; Fig. 8.1. |
| `check_sampling(D1, D2, delta1, delta2, N, wavelength_m, z, R=inf)` | Five `Rule(name, satisfied, bound, actual, citation)` tuples. Ch. 7. |

### 10d. The turbulence and the screens (`turbulence.py`)

| Name | What it gives |
|---|---|
| `phase_psd(f, r0_m, L0_m=inf, l0_m=0.0)` | The one shared phase PSD expression. Ch. 9, Eqs. (9.49) to (9.52). |
| `kolmogorov_phase_psd`, `von_karman_phase_psd`, `modified_von_karman_phase_psd` | The three named spectra, as wrappers over `phase_psd`. Ch. 9, Eqs. (9.49), (9.50), (9.51). |
| `kolmogorov_structure_function(r_m, r0_m)` | `D(r) = 6.88 (r/r0)^(5/3)`. Ch. 9, Eq. (9.44). |
| `ft_phase_screen(r0_m, n, dx_m, L0_m=inf, l0_m=0.0, rng=None)` | The Fourier-series screen. Ch. 9, Eqs. (9.78) to (9.80). |
| `subharmonic_screen(..., n_p=3)` | The low-frequency part alone. Ch. 9, Eq. (9.81). |
| `ft_sh_phase_screen(..., n_p=3)` | The sum of the two. Ch. 9, Listings 9.2 and 9.3. |
| `screen_r0(cn2_integral_m13, wavelength_m)` | The per-screen Fried parameter. Ch. 9, Eq. (9.70). |
| `composite_r0(r0_i_m, alpha_i=None, wave='plane')` | The composite Fried parameter. Ch. 9, Eqs. (9.71), (9.72). |
| `screen_rytov_share(...)` | One screen's share of the LOG-AMPLITUDE variance. Ch. 9, Eqs. (9.73), (9.74). |
| `max_screen_strength(...)`, `RMAX = 0.1` | The per-screen cap of Listing 9.5, lines 37 to 39. |
| `screen_strengths(...)` | The bounded least-squares solve for the screen strengths. Ch. 9, Eq. (9.75). |
| `profile_moments`, `layer_moments`, `moment_error` | The layer moment rule for `0 <= m <= 7`. Ch. 9, Eq. (9.65). |
| `fresnel_pitch_max(wavelength_m, z_m)` | The `sqrt(lambda z)/2` pitch cap. Sec. 9.4. |
| `phase_pitch_max(r0_m, max_step_rad=pi, sigmas=3.0)` | The Johnston and Lane phase pitch rule, `dx <= 0.332 r0`. Sec. 9.4 with Eq. (9.44). |
| `blurred_extent(d_m, wavelength_m, dz_m, r0_m, c=2.0)` | The turbulence-blurred extents D1', D2'. Ch. 9, Eqs. (9.84), (9.85). |
| `constraint1_pitch_max`, `constraint2_n_min`, `constraint3_pitch_range` | The three turbulent geometry constraints. Ch. 9, Eqs. (9.86) to (9.88). |
| `max_partial_step`, `min_planes` | The turbulent step cap and the plane count. Ch. 9, Eqs. (9.89), (9.90). |
| `WEAK_SIGMA2_CHI = 0.25` | The weak-fluctuation threshold on the log-amplitude variance. Ch. 9, text below Eq. (9.64). |
| `properly_sampled_checklist(...)` | One `(rule, satisfied, bound, actual, citation)` tuple per step of Sec. 9.5. An advisory step returns `satisfied = None`. |

### 10e. The example scripts

Three scripts in `examples/schmidt/` put this layer against the production
layer. Each one reads `olb` and changes nothing in it.

- `propagator_kernels.py` â€” the book kernels against `Forvard`, `Fresnel` and
  the co-moving `Lens -> LensFresnel -> Convert` recipe, in three tiers.
- `sampling_and_edges.py` â€” a gallery of deliberate sampling failures, then the
  rule checker on the real production grids.
- `screens_and_turbulence.py` â€” the book screen generators against the
  `aotools` generator and against Eq. (9.44), plus the factor-4 bridge between
  the two per-screen variance conventions.

See [examples/schmidt/README.md](../examples/schmidt/README.md) for the
measured numbers and the wiring status.

---

## 11. The fidelity-2 runner and the Terms (`olb/models/waveoptics.py`)

This module sits OUTSIDE `olb.waveoptics`. It is the bridge from the layer above
to the budget: it runs the propagations one time, and it turns the records into
Terms. A budget NEVER runs a simulation; the caller precomputes the records and
passes them in. See [api-budget.md](api-budget.md) for the budget side.

| Name | What it gives |
|---|---|
| `run_fidelity2(...)` | The runner. It gives the `Fidelity2Bundle` (or a list of them) that a fidelity-2 budget needs. |
| `run_waveoptics(...)` | The turbulent run alone. It gives a `TurbWaveResult`. |
| `Fidelity2Bundle` | The two records: `vacuum` and `turbulent`. |
| `waveoptics_turbulence_term(result, ...)` | The stochastic turbulence Term, from the per-trial scalars. |
| `waveoptics_smf_coupling_term(result, ...)` | The turbulent SMF-coupling face (`quantity="smf_eta"`). |
| `waveoptics_mmf_coupling_term(result, ...)` | The turbulent MMF-coupling face (`quantity="mmf_eta"`). |
| `waveoptics_vacuum_term(result, ...)` | The deterministic vacuum-optics Term (launch to detector, no fade). |
| `waveoptics_vacuum_mmf_term(vacuum_result, detector, aperture_m, ...)` | The deterministic vacuum MMF core-capture Term (no fade). |

### 11a. `run_fidelity2(scenario, geometry, *, n_trials=200, preset="standard", seed=None, threader=None, cn2=None, hs=None, cn2_profile=None, h_top_m=None, L0_m=np.inf, subharmonics=True, progress=True, vacuum=None, turbulence=True, detectors=None, precision="single")`

It runs the wave-optics propagations that a fidelity-2 budget needs, one time
each: the TURBULENT split-step Monte Carlo (the fade), and a no-turbulence
GEOMETRIC loss.

`vacuum` selects the source of the geometric loss. `None` (the default) takes
`"analytic"` for a space link and `"wave"` for a terrestrial link. `"analytic"`
makes NO wave vacuum run, and the bundle `vacuum` is then `None`: a ground-space
link is far field, so the analytic geometric Term is exact and cheap, and the
full-path field solve is slow and grid-noise-limited. `"wave"` makes the run;
a space link opts back in that way, for research or a cross-check. `"analytic"`
raises for a terrestrial link, because the near-field penalty needs the vacuum
baseline on the SAME flat grid. An unknown value raises `ValueError`.

**THE MASTER TURBULENCE SWITCH (`turbulence`).** `True` (the default) runs the
split-step Monte Carlo. `False` SKIPS it fully: no screens, and no trials. The
bundle is then VACUUM-ONLY, with `turbulent=None`, and the budget shows the
deterministic Terms alone. The switch MIRRORS the fidelity-0 master `turbulence`
switch of `olb.links.terrestrial.terrestrial_budget`, so the fidelity ladder
reads the same at each rung. Pass `turbulence=False` to the budget too.

- A TERRESTRIAL vacuum-only run still sizes the grid with `turbulent_grid()` and
  propagates on that SAME grid. So the vacuum Term does NOT move when the caller
  toggles the switch.
- A SPACE vacuum-only run skips the grid sizing too. With the default
  `vacuum="analytic"` the bundle is then EMPTY (`vacuum=None`,
  `turbulent=None`), and that is VALID: the budget shows the analytic
  deterministic Terms alone. Use `vacuum="wave"` to get the receive-plane field
  of the co-moving vacuum solve.

**THE BEAMSPLITTER ARMS (`detectors`).** `detectors` takes a sequence of detector
objects, the arms behind a receive beamsplitter. The split-step Monte Carlo then
runs ONE time, and every arm reads the SAME field (Section 9d), so N arms cost
one run. This is EXACT: a beamsplitter scales the field of an arm by a constant,
and every coupling efficiency is power-normalised, so the split ratio does not
change it (the fraction is a separate fixed dB Term, `olb.models.splitter`). The
function then returns a LIST of `Fidelity2Bundle`, one for each arm, in the
`detectors` order. Each arm carries its own efficiency on the `smf_eta` or the
`mmf_eta` face of the shared trials, so the Term factories read it with no
change.

The VACUUM run is PER ARM, because a vacuum record holds the fibre coupling of
its own detector. It is one deterministic propagation, so it is cheap for a
terrestrial link. A space link with the default `vacuum="analytic"` makes NO
vacuum run at all; `vacuum="wave"` makes one full-path solve for each arm, and
that is slow (about 14 s each).

`progress=True` (the default) prints a recap of the auto-chosen grid, the screen
plan and the sampling quality, then it shows a tqdm bar over the turbulent
trials. The one-time vacuum run has no bar. Pass `progress=False` for a quiet
run.

#### `Fidelity2Bundle`

A frozen dataclass. The two wave-optics records a fidelity-2 budget needs.

| Field | Type | Meaning |
|---|---|---|
| `vacuum` | `WaveResult` or None | The record of one no-turbulence propagation (Section 6). It gives the geometric spread, the launch truncation, the aperture capture, and the vacuum fibre coupling. It is `None` when the geometric loss is analytic (`vacuum="analytic"`, the space default). |
| `turbulent` | `TurbWaveResult` or None | The record of the split-step Monte Carlo (Section 9d). It gives the turbulence penalty. It is `None` for a VACUUM-ONLY bundle from `run_fidelity2(turbulence=False)`. |

### 11b. `waveoptics_vacuum_mmf_term(vacuum_result, detector, aperture_m, *, beam_type=BEAM_GAUSSIAN, name=None, note=None, meta_extra=None)`

The DETERMINISTIC vacuum MMF core-capture Term, of the category `"coupling"`. It
carries NO fade.

A vacuum-only bundle has no per-trial `mmf_eta`, because it makes no trials. But
it holds the receive-clipped field, so the light-bucket core capture is a direct
calculation on that field. This Term is that calculation: the fraction of the
COLLECTED power that enters the fibre core. It reads stage 3 of the record
(`"after rx clip"`), and it calls `mmf_coupling_efficiency` of Section 4a with
the `numerical_aperture` and the `defocus_m` of the detector. Because it is
relative to the COLLECTED power, it composes with the vacuum-optics Term
(`waveoptics_vacuum_term`, launch to collected power) with NO double-count.

It is the vacuum companion of `waveoptics_mmf_coupling_term`. The focal length
follows the SAME rule as the turbulent runner: an explicit `MMF.focal_length_m`
wins; else `MMF.optimal_focus` matches the spot to the core through the
`a = 1.12` spot-to-core parameter (Shaklan and Roddier, Appl. Opt. 27 (1988)
2334, DOI 10.1364/AO.27.002334); else it raises `ValueError`.

`aperture_m` is the receive aperture DIAMETER, in m. The Term meta holds
`mmf_eta`, `focal_length_m` and `defocus_m`.

**The SMF leg reads the defocus too** (backlog 2-W2, DONE 2026-09-04).
`olb.waveoptics.smf.coupling_efficiency` takes `defocus_m` and
`focal_length_m`, and it applies the SAME quadratic pupil phase through the
shared helper `olb.waveoptics.mmf.defocus_phase`. So a fidelity-2 SMF budget
charges `SMF.defocus_m`, with the same sign as the multimode leg: a diverging
received beam couples best at a POSITIVE `defocus_m`. The focal length comes
from `SMF.focal_length_m`, or from `SMF.optimal_focus`
(`f = pi*(D/2)*w_m/(lambda*1.12)`); a defocus with neither raises `ValueError`.
`defocus_m=0.0` (the default) keeps the old focal-plane overlap exactly.
