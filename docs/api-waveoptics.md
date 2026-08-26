# API reference: the wave-optics layer

This page documents `olb.waveoptics`, the fidelity-2 field propagation layer. It
gives the exact signatures and defaults from the source.

The package propagates a scalar complex field through free space on a square
grid. There is NO turbulence in it. The core is a trimmed port of LightPipes
(https://github.com/opticspy/lightpipes, BSD-3-Clause). See
[`LIGHTPIPES_LICENSE.txt`](../olb/waveoptics/LIGHTPIPES_LICENSE.txt) in the
package. The port keeps the LightPipes names and the LightPipes call order, so a
script from that package runs here with no change.

Import from the package root:

```python
from olb.waveoptics import Begin, GaussBeam, Fresnel, GridSpec, propagate_scenario
```

Status: the package builds NO Term and it changes NO budget. It is the
no-turbulence validator for the near-field and far-field limits of the analytic
Terms. A fidelity-2 Term is an owner-gated later step.

The core (`field.py`, `sources.py`, `propagators.py`, `smf.py`) imports numpy and
scipy only. It imports nothing from the rest of `olb`. Only `grid.py` and `run.py`
read a scenario.

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

### The field functions

- `Begin(size, labda, N)` — make a new `Field` of `N` x `N` pixels. The amplitude
  is 1.0 at each pixel.
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

---

## 4. The fibre coupling (`olb/waveoptics/smf.py`)

- `smf_mode(grid_size_m, wavelength_m, n, aperture_m)` — the back-propagated
  single-mode-fibre mode in the pupil plane. It is a Gaussian of the radius
  `aperture_m / MODE_RADIUS_RATIO`. The intensity sum is 1.0.
- `coupling_efficiency(field, aperture_m, mask=None)` — the power fraction that
  couples into the fibre, a float between 0 and 1. It is the normalised overlap
  of the field with the mode. `aperture_m` is the pupil DIAMETER. An optional
  `mask` multiplies the field first. A field with no power raises `ValueError`.
- `MODE_RADIUS_RATIO = 2.24` — the best ratio of the pupil diameter to the
  pupil-plane mode radius (Ruilier, DOI 10.1117/12.317094). A flat pupil then
  couples at the maximum of 0.8145.

The function takes ONE field and gives ONE float. Loop in the caller for a set of
realisations.

---

## 5. The grid (`olb/waveoptics/grid.py`)

### `GridSpec(size_m, n)`

A frozen dataclass. It holds the two numbers that a propagation needs.

| Field | Type | Unit | Meaning |
|---|---|---|---|
| `size_m` | float | m | The physical side of the square grid. |
| `n` | int | — | The number of pixels along one side. |
| `pixel_m` | float | m | A property. The pitch, `size_m / n`. |

### `GridSpec.for_scenario(scenario, geometry, guard=4.0, pixels_per_feature=16, n_max=4096)`

Derive a grid from a scenario and a geometry.

- The EXTENT rule: `size = guard * 2 * r_max`. `r_max` is the largest of the beam
  radius at the launch plane, the beam radius at the longest range, the transmit
  aperture radius, and the receive aperture radius.
- The RESOLUTION rule: the smallest feature gets `pixels_per_feature` pixels
  across it. The pixel count goes up to the next power of two. It stays in the
  interval `[256, n_max]`.
- The transmit aperture obeys the bistatic rule of
  `olb.models.gaussian_efficiency`.

The method WARNS. It does not raise. It warns when the `n_max` clamp leaves too
few pixels on the smallest feature, and when the range is longer than
`forvard_max_z()`. A transmit terminal with no `Transmitter` raises `ValueError`.

### `forvard_max_z(grid, wavelength_m)`

The longest range that the grid samples well, `z_max = N * dx^2 / lambda`. Past
this range the quadratic phase of the transfer function turns faster than one
sample, so `Forvard` aliases. See Schmidt, DOI 10.1117/3.866274, Ch. 6.

`N_MIN = 256` is the smallest pixel count.

---

## 6. One end-to-end propagation (`olb/waveoptics/run.py`)

### `propagate_scenario(scenario, geometry, grid=None)`

Propagate the transmit beam of a scenario to the receive aperture. The steps are:
the launch, the launch-aperture clip, the free-space propagation, the
receive-aperture clip, and the fibre coupling. `grid=None` derives the grid with
`GridSpec.for_scenario()`. A geometry that gives more than one range raises
`ValueError`. Loop in the caller for a sweep.

A deliberately diverged beam starts at a virtual waist behind the aperture (see
`olb.beam`). Then the beam has the asked-for radius in the aperture plane.

### `WaveResult`

A frozen dataclass. All the losses are positive dB.

| Field | Type | Meaning |
|---|---|---|
| `stages` | list | The `(label, Field)` pairs. The labels are `"launch"`, `"after tx clip"`, `"at rx plane"` and `"after rx clip"`. |
| `grid` | `GridSpec` | The grid that the propagation used. |
| `tx_truncation_db` | float | The power that the launch aperture takes. |
| `geometric_loss_db` | float | The power that the receive aperture does not collect. |
| `smf_coupling_db` | float or None | The single-mode-fibre coupling loss. `None` when the receive terminal has no `SMF` detector. |
| `propagator` | str | The name of the propagator that ran. |

`PURE_GAUSS_CLIP = 1e-6` is the dispatch threshold. See Section 7.

---

## 7. The propagator regimes and their limits

| Propagator | Method | Field | Distance limits | Grid limits |
|---|---|---|---|---|
| `GForvard` | Analytic ABCD (Siegman, ISBN 978-0935702118) | A pure Gaussian beam ONLY. It refuses a clipped field. | Any `z`. | None. There is no grid error, because the route is analytic. |
| `Fresnel` | Convolution on a doubled grid (Schmidt, DOI 10.1117/3.866274, Ch. 7) | Any field. | MINIMUM `z`. The result is not valid when `z` is comparable with, or less than, the size of the diffracting aperture. (LightPipes manual; Schmidt Ch. 7) | No periodic wrap: the doubled grid absorbs it. The cost is 8 times the memory. The field must be zero at the grid edges. |
| `Forvard` | FFT angular spectrum (Schmidt, DOI 10.1117/3.866274, Ch. 6) | Any field. | Any `z` down to zero. But EACH call needs `z < forvard_max_z = N*dx^2/lambda`. Past that limit the transfer function aliases. | The boundary is periodic. A beam that reaches the edge wraps to the opposite edge. Give the grid a side of about 8 times the largest beam radius. |

### The dispatch rule of `propagate_scenario`

`propagate_scenario` reads the power that the launch aperture takes:

- The clip takes less than `PURE_GAUSS_CLIP = 1e-6` of the power: the field stays
  a pure Gaussian. The function propagates the UNCLIPPED launch field with
  `GForvard`, the exact route.
- The clip takes more: the field carries the aperture edge. The function
  propagates the CLIPPED field with `Fresnel`.

`propagate_scenario` does not call `Forvard`. A space link is longer than
`forvard_max_z` on any practical grid, so `GridSpec.for_scenario()` warns there.

---

## 8. The LightPipes propagators that the port does not hold

Each of these serves a regime that the current layer does not need. Add one only
when its trigger comes.

| Function | Regime | The trigger to add it |
|---|---|---|
| `Forward` | The direct integral. The output grid is different from the input grid. | You must magnify or shrink the grid between two planes. |
| `LensForvard`, `LensFresnel` with `Convert` | The co-moving spherical grid. The grid follows the beam, so a large expansion stays sampled. | This is the route to a samplable space link. Add it when a fidelity-2 space case must run. |
| `Steps` | Propagation through a medium with an index term. It holds a built-in absorbing boundary. | The turbulent split-step layer, or an absorbing edge (see Section 9). |
| `Interpol` | A regrid of the field: a new side, a new pixel count, a shift, or a rotation. | You must pass a field between two propagators that need different grids. |

---

## 9. Sampling design for the split-step layer (not built)

The turbulent split-step layer is NOT built. This section states the constraints
that it must satisfy. It does not design the layer.

Three constraints apply, and two of them look like a conflict:

1. **The screen spacing must be SHORT.** Each slab between two phase screens must
   scatter weakly: its Rytov variance must stay small. A slab that is too long
   puts the diffraction of the slab into a single screen, so the layer loses the
   correct irradiance statistics. See Martin and Flatte,
   DOI 10.1364/AO.27.002111, and Schmidt, DOI 10.1117/3.866274, Ch. 9.
2. **`Fresnel` cannot take a short hop.** The convolution method has a minimum
   distance (see Section 7). So `Fresnel` is NOT the split-step propagator.
3. **`Forvard` is the split-step propagator.** It has no minimum distance, and
   its `forvard_max_z` limit is PER CALL, not per path. A short hop therefore
   satisfies the weak-scatter bound and the sampling bound at the same time. A
   long path is many short `Forvard` calls, each one below `forvard_max_z`.

The split-step layer must also add three things that the current layer does not
have:

- A pixel pitch that resolves the coherence radius `r0`. The screen carries no
  phase structure below one pixel.
- A grid extent that holds the scattered spread, not only the free-space beam
  radius. Turbulence widens the beam.
- An absorbing edge mask between the hops. `Forvard` is periodic, so the light
  that turbulence pushes to the edge comes back. A super-Gaussian mask removes
  it. LightPipes solves this inside `Steps`.

### 9a. Where to put the screens on a non-uniform path

The three constraints above set how SHORT a slab must be. They do not set WHERE
the screen boundaries go. A slant path is not uniform: the ground holds most of
the `Cn2`, and the high air holds little. So equal `Cn2` weight, NOT equal
distance, sets the boundaries.

- **Equal-strength partition.** Put the boundaries so that each slab holds the
  same turbulence integral: the same `integral of Cn2 dz`. This is the same
  quantity that sets the Fried parameter, `r0 = (0.423 k^2 integral Cn2 dz)^(-3/5)`.
  So equal `Cn2` weight is equal `r0^(-5/3)` per slab, and equal phase variance
  per screen. See Lane, Glindemann and Dainty, DOI 10.1088/0959-7174/2/3/003, and
  Coles, Filice, Frehlich and Yadlowsky, DOI 10.1364/AO.34.002089.
- **Why not equal distance.** Equal `dz` puts many screens in the thin high air,
  where they add nothing, and too few screens at the ground, where the weak-scatter
  bound breaks first. The screens bunch near the ground on an uplink.
- **Reuse the profile.** The `Cn2` integral is the one that
  `olb.turbulence.profiles.default_cn2_profile` already gives. The partition is
  the inverse of its running integral, so the split-step layer must NOT hold its
  own profile.

The screen count comes from the weak-scatter bound and the total strength:
`N >= sigma_R^2 / sigma_per_screen^2`, with `sigma_per_screen^2` about 0.1. See
Martin and Flatte, DOI 10.1364/AO.27.002111.

### 9b. The convergence check

The bounds above give a screen count that SHOULD work. The honest check is to
prove it. Run the propagation. Then double the screen count and run it again. If
the receiver scintillation index (or the coupling efficiency) moves less than the
tolerance, the first count was enough. Martin and Flatte validated the method this
way, DOI 10.1364/JOSAA.7.000838.

This check fits the layer style of this package: it WARNS, it does not raise. A
path that does not converge at a practical screen count is an honest warning, the
same as the `GridSpec` warnings in Section 5.
