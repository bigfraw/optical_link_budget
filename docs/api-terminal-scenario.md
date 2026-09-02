# API reference: the pure-data layer

This page documents the pure-data layer of the `olb` package. These classes
hold the inputs that the models read. They compute no physics. The models read
a scenario and a geometry, and they return Terms.

The data moves in one direction, from the inputs to the models. A scenario does
not import the models. A terminal does not import the models.

Three rules hold across this layer:

- ALL terminal hardware lives on a `Terminal`. A channel holds no hardware.
- A terminal parameter can only be set through a `Terminal`.
- Loss is positive dB. Gain is negative dB.

Modules: `olb.terminal`, `olb.scenario`, `olb.geometry`.

---

## 1. Terminal hardware (`olb.terminal`)

A `Terminal` groups a telescope aperture, an optional transmitter, an optional
detector, and an optional compensation stack. One `Terminal` serves both link
directions. The scenario resolves which terminal transmits and which receives.

### `Terminal`

One optical terminal: aperture, transmitter, compensation, and detector.

| Field | Type | Unit | Default | Meaning |
|---|---|---|---|---|
| `aperture_m` | float | m | (required) | Telescope aperture diameter. |
| `obscuration_ratio` | float | — | `0.0` | Central obscuration diameter divided by aperture diameter. 0 means unobscured. |
| `wavelength_m` | float | m | `1550e-9` | The terminal operating wavelength. |
| `pointing_jitter_rad` | float | rad | `0.0` | 1-sigma tracking jitter. 0 means no jitter. |
| `transmitter` | `Transmitter` or None | — | `None` | The transmit source. None means the terminal only receives. |
| `detector` | `Aperture`, `SMF`, `MMF`, or `Camera`, or None | — | `None` | The detector front end. None means no receive-coupling Term. |
| `compensation` | list of `TipTilt` or `AO` | — | `[]` | The ordered wavefront-compensation stack. It may be empty. |

Constraints:

- Each `Terminal` owns its own `compensation` list. Two terminals do not share
  one list. The default is a fresh empty list per terminal.
- An empty `compensation` stack leaves the piston-removed turbulence.
- A `Terminal` holds ONE detector. A receive path that feeds SEVERAL detectors
  behind a beamsplitter (for example a tracking `Camera` and a comms fibre)
  stays one detector per `Terminal`. Each detector carries its splitter fraction
  `frac`, and the budget helper `olb.multidetector` makes one `Terminal` for
  each arm.

### `Transmitter`

The transmit source of a terminal. The transmitter launches a Gaussian beam
through a launch aperture.

By default the launch aperture is the owning `Terminal` aperture. The launch
truncation reads the `Terminal` `aperture_m` and `obscuration_ratio`. This is a
MONOSTATIC terminal. One aperture transmits and receives.

For a BISTATIC terminal the transmit beam director is a different aperture from
the receive telescope. Set `aperture_m` on the `Transmitter`. The launch
truncation then reads the transmitter values. The `Terminal` `aperture_m` and
`obscuration_ratio` then describe the receive telescope only.

| Field | Type | Unit | Default | Meaning |
|---|---|---|---|---|
| `waist_m` | float | m | (required) | Transmit Gaussian waist (1/e^2 radius) at the aperture. |
| `power_dbm` | float or None | dBm | `None` | Launch power. None if only losses matter. |
| `m2` | float | — | `1.0` | Beam quality M^2 (>= 1). |
| `divergence_rad` | float or None | rad | `None` | Transmit far-field 1/e^2 half-angle divergence. None means collimated (the diffraction limit). |
| `aperture_m` | float or None | m | `None` | Transmit (beam director) aperture diameter. None means the transmitter shares the owning `Terminal` aperture (monostatic). |
| `obscuration_ratio` | float or None | — | `None` | Central obscuration ratio of the transmit aperture. None means the transmitter shares the owning `Terminal` `obscuration_ratio`. Set 0.0 for an unobscured beam director on a terminal whose receive telescope is obscured. |

Constraints:

- `m2` must be at least 1.
- A `None` for `aperture_m` or `obscuration_ratio` keeps the monostatic default.
  The value then comes from the owning `Terminal`.

### Detectors

A detector is one of four front ends: `Aperture`, `SMF`, `MMF`, or `Camera`.
Each carries an optional receive sensitivity.

Each detector also carries a splitter fraction `frac`. It is the fraction of the
received power that the beamsplitter sends to this detector (0 to 1). `None`
means "take the remainder", so a detector that is alone gets 1.0. The fraction
rule lives in `olb.models.splitter.resolve_fracs`, and the per-arm budgets are in
`olb.multidetector.multi_detector_budgets` (see `api-budget.md`). The field is
pure data: `olb.terminal` computes no physics.

#### `Aperture`

Power-in-bucket detector. An aperture (bucket) detector integrates the intensity
over the aperture. It is phase-insensitive. So a compensation stack does not
change its coupling. Use it for parity with the plain downlink budget.

| Field | Type | Unit | Default | Meaning |
|---|---|---|---|---|
| `sensitivity_dbm` | float or None | dBm | `None` | Required received power. None if only losses matter. |
| `frac` | float or None | — | `None` | The fraction of the received power that the beamsplitter sends to this detector (0 to 1). None means "take the remainder" (1.0 when the detector is alone). |

#### `SMF`

Single-mode-fibre detector. A single-mode fibre couples only the field that
matches the fibre mode. The coupling falls when turbulence distorts the
wavefront.

A focusing optic of focal length `f` puts the collected beam onto the fibre tip.
Set `focal_length_m` and `mode_field_radius_m` to derive `eta_max` from the
optics. The model finds the coupling parameter `a = pi*(D/2)*w_m/(lambda*f)` and
`eta_max(a)` from the mode-overlap curve. The peak `eta_max=0.8145` is near
`a=1.12` (see `physics.md` section 6c). With `focal_length_m` None, the model
uses the `eta_max` field.

A single-mode-fibre subtlety: at a fixed `a` the focal length cancels in the
tilt-to-coupling response. So `f` is a static knob only. It sets `eta_max`
through `a`. It does not change the angular sensitivity on its own.

| Field | Type | Unit | Default | Meaning |
|---|---|---|---|---|
| `eta_max` | float | — | `0.8145` | Maximum fibre-to-aperture mode match for an unobscured circular aperture with a flat wavefront. Used when `focal_length_m` is None. |
| `sensitivity_dbm` | float or None | dBm | `None` | Required received power. None if only losses matter. |
| `focal_length_m` | float or None | m | `None` | Focal length of the fibre-coupling optic. None keeps the `eta_max` field. A value needs `mode_field_radius_m`. |
| `mode_field_radius_m` | float or None | m | `None` | Fibre mode field RADIUS (about 5.2e-6 m for SMF-28 at 1550 nm). It sets the fibre mode size for `a` and for the walk-off Term. |
| `optimal_focus` | bool | — | `False` | Design the coupling optic for the best coupling (see below). |
| `defocus_m` | float | m | `0.0` | Detector offset from the design focus. The fibre tip sits at `z = f + defocus_m`. `0.0` puts it at the nominal focal plane. |
| `frac` | float or None | — | `None` | The fraction of the received power that the beamsplitter sends to this detector (0 to 1). None means "take the remainder" (1.0 when the detector is alone). |

`optimal_focus=True` assumes the optimal coupling parameter `a=1.12`, so
`eta_max=0.8145`, and derives the focal length from the mode field radius and the
aperture: `f = pi*(D/2)*w_m/(lambda*1.12)`. A None `mode_field_radius_m` uses the
SMF-28 value (5.2e-6 m). An explicit `focal_length_m` overrides the derived
value. A bare `SMF()` (this flag False) does not change: it stays mean-only, with
no walk-off Term.

`optimal_focus` is a focal-LENGTH rule only. It never moves the detector, so it
does not put the fibre at the true focus of a curved received beam. A terrestrial
received beam diverges, so its true focus is BEYOND the focal plane, at
`z = f + dz_curv` (see `physics.md` section 6a). The terrestrial coupling Terms
ALWAYS charge that curvature defocus at the actual fibre plane. To model a
TRACKED (aligned) coupler, set

```python
from olb.models.coupling import curvature_focus_shift
detector.defocus_m = curvature_focus_shift(scenario)
```

A space link has an enormous `R_rx`, so its `dz_curv` is about zero.

#### `MMF`

Multimode-fibre detector: a light bucket in the fibre plane. A multimode fibre
accepts all the light that the focusing optic puts inside its core. The core is a
disk of fixed radius `core_radius_m`. So the coupling is a geometric overlap of
the focal spot with the core, not a modal overlap. A received tip-tilt of angle
theta moves the spot by `f*theta`, so the spot walks off the core when the
tip-tilt is large (see `physics.md` section 6c).

Unlike a single-mode fibre, the focal length does NOT cancel. The core is a fixed
physical size, so it subtends the angular field of view `core_radius_m/f`. So the
focal length is a genuine free parameter.

The fibre also has an ANGULAR limit. The numerical aperture `NA = n*sin(theta_a)`
sets the largest ray angle the fibre guides. The focusing optic makes a cone of
half-angle `NA_optic = (D/2)/f`. When `NA_optic > NA` the fibre does not guide the
steep rays, so the coupled power falls by `min(1, (NA/NA_optic)^2)`. This is the
etendue penalty a core-radius-only bucket misses (see `physics.md` section 6c).

| Field | Type | Unit | Default | Meaning |
|---|---|---|---|---|
| `core_radius_m` | float | m | (required) | Core RADIUS of the multimode fibre in the fibre plane. |
| `focal_length_m` | float or None | m | `None` | Focal length of the fibre-coupling optic. None needs `optimal_focus=True`. |
| `numerical_aperture` | float or None | — | `None` | Fibre NA. None turns the angular gate OFF (spatial encircled energy only). A value gates the focusing cone by `min(1, (NA/NA_optic)^2)`. |
| `sensitivity_dbm` | float or None | dBm | `None` | Required received power. None if only losses matter. |
| `optimal_focus` | bool | — | `False` | Match the spot to the core (see below). |
| `defocus_m` | float | m | `0.0` | Detector offset from the design focus. The core sits at `z = f + defocus_m`. `0.0` puts it at the nominal focal plane. A detector away from the TRUE focus sees a larger spot, so the core captures less. |
| `frac` | float or None | — | `None` | The fraction of the received power that the beamsplitter sends to this detector (0 to 1). None means "take the remainder" (1.0 when the detector is alone). |

`optimal_focus=True` derives the focal length so the spot radius is the core
radius over 1.12 (the same `a=1.12` that a single-mode fibre uses):
`f = pi*(D/2)*core_radius_m/(lambda*1.12)`. This gives about 92% static capture.
It is a geometric spot-to-core match, NOT a mode-overlap optimum: a shorter focal
length captures more, but the angular limit (`numerical_aperture`) then gates the
extra capture. Set `focal_length_m` to override the derived value. As for an
`SMF`, `optimal_focus` never moves the detector: the received-curvature focus
shift is charged at the actual fibre plane, and
`curvature_focus_shift(scenario)` gives the `defocus_m` of an aligned coupler
(see `physics.md` section 6a).

Import it as `from olb import MMF` (a top-level export), or as
`from olb.terminal import MMF`.

#### `Camera`

Focal-plane array detector: a tracking and spot-diagnostic sensor. A camera
images the focal spot on a grid of square pixels. It measures the spot SHAPE and
the spot POSITION. So it is the sensor of a tracking loop, and it is the
diagnostic front end of a wave-optics study.

`pixel_pitch_m` is the centre-to-centre distance of two pixels (the pixel
scale). The sensor is square: its side is `n_pixels * pixel_pitch_m`. One pixel
subtends the angle `pixel_pitch_m / focal_length_m` on the sky, so the focal
length sets the plate scale. A measured centroid `x` maps to the arrival angle
`theta = x / focal_length_m`.

`focal_length_m` is the imaging optic. `defocus_m` puts the sensor at
`z = f + defocus_m`, so `0.0` is the focal plane. This is the same convention as
`SMF` and `MMF`. See `olb.waveoptics.camera` for the focal-plane
discretisation.

| Field | Type | Unit | Default | Meaning |
|---|---|---|---|---|
| `pixel_pitch_m` | float | m | (required) | Centre-to-centre distance of two pixels (the pixel scale). |
| `n_pixels` | int | — | (required) | Number of pixels along one side. The sensor is square. |
| `focal_length_m` | float or None | m | `None` | Focal length of the imaging optic. None means the plate scale is not set, so the angular scale cannot be computed. |
| `defocus_m` | float | m | `0.0` | Sensor offset from the focal plane. The sensor sits at `z = f + defocus_m`. `0.0` is the focal plane. A non-zero value grows the spot. |
| `sensitivity_dbm` | float or None | dBm | `None` | Required received power. None if only losses matter. |
| `frac` | float or None | — | `None` | The fraction of the received power that the beamsplitter sends to this detector (0 to 1). None means "take the remainder" (1.0 when the detector is alone). A tracking camera usually takes a small fraction, and the comms fibre takes the remainder. |

LIMIT, budgets. No budget builds a coupling Term for a `Camera` today. The
dispatch is not the same in each budget, so read this before you put a `Camera`
on a budgeted terminal:

- `terrestrial_budget` treats a `Camera` like an `Aperture` (a power-in-bucket
  receiver): it adds the scintillation Term and no coupling Term.
- `downlink_budget` at `fidelity=2` also treats a `Camera` like an `Aperture`.
- `downlink_budget` at fidelity 0 or 1 RAISES `ValueError` ("unknown
  detector"), because `olb.models.coupling.downlink` knows `Aperture` and `SMF`
  only.

So use a `Camera` for the wave-optics focal-plane tools, and use an `Aperture`
for a power budget.

Import it as `from olb import Camera` (a top-level export), or as
`from olb.terminal import Camera`.

### Compensation stack

The compensation stack is an ordered list of correction stages. The residual
wavefront that the stack leaves sets the coupling into the detector. The
compensation and the detector are one physical chain. So the model emits one
receive-coupling Term, not two.

#### `TipTilt`

Tip-tilt correction stage. It removes the first three Zernike modes (piston,
tip, tilt). It has no fields.

#### `AO`

Adaptive-optics correction stage. It removes the first `n_modes` Zernike modes.

| Field | Type | Unit | Default | Meaning |
|---|---|---|---|---|
| `n_modes` | int | — | `20` | Number of Zernike modes that the stage removes. The model uses the large-order Noll asymptotic residual. |

### Snippet: monostatic and bistatic terminals

```python
from olb.terminal import Terminal, Transmitter, Aperture, SMF, TipTilt, AO

# Monostatic terminal: one 0.15 m aperture transmits and receives.
# The launch truncation reads the Terminal aperture_m and obscuration_ratio.
monostatic = Terminal(
    aperture_m=0.15,
    obscuration_ratio=0.3,
    transmitter=Transmitter(waist_m=0.12, power_dbm=40.0),
)

# Bistatic terminal: a small 0.15 m beam director transmits.
# A large 0.7 m telescope receives. The Transmitter carries its own aperture,
# so the launch truncation does not read the receive-telescope aperture.
bistatic = Terminal(
    aperture_m=0.7,
    obscuration_ratio=0.3,
    transmitter=Transmitter(waist_m=0.06, aperture_m=0.15, obscuration_ratio=0.0),
    detector=SMF(sensitivity_dbm=-40.0),
    compensation=[TipTilt(), AO(n_modes=60)],
)
```

---

## 2. Scenario families (`olb.scenario`)

A scenario is a link case. There are two families, and they share one interface.
A link is either a space link (a ground station and a satellite) or a
terrestrial link (two ground stations on a horizontal path).

Both families expose the same thin interface that the models read:

- `scenario.tx_terminal` — the transmit terminal.
- `scenario.rx_terminal` — the receive terminal.
- `scenario.channel` — the propagation channel.

So no model changes between the two families.

### `SpaceScenario`

A space link case: a ground terminal, a space terminal, and a direction.

| Field | Type | Unit | Default | Meaning |
|---|---|---|---|---|
| `ground` | `Terminal` | — | (required) | The ground-station terminal. |
| `space` | `Terminal` | — | (required) | The satellite terminal. |
| `direction` | `"uplink"` or `"downlink"` or `"retro"` | — | `"uplink"` | The link direction. It sets the tx / rx roles. |
| `channel` | `Channel` | — | `Channel()` | The space propagation channel (site plus orbit altitude). |
| `availability_target` | float | — | `0.99` | Target link availability (0-1). |
| `precompensation` | `DownlinkBeacon` or `LaserGuideStar` or None | — | `None` | The uplink pre-compensation source. None means the uplink is uncorrected. A downlink or a retro scenario refuses the field at construction (`ValueError`). |

The `direction` sets the roles:

| direction | tx_terminal | rx_terminal |
|---|---|---|
| uplink | ground | space |
| downlink | space | ground |
| retro | ground | ground |

#### Pre-compensation sources

The `precompensation` field names what the ground terminal senses to build the
uplink correction. It applies to the uplink direction only. The uplink budget
reads it (see `api-budget.md`).

- `None` — no source. The uplink is uncorrected.
- `DownlinkBeacon()` — the ground terminal senses the downlink beam and applies
  the conjugate to the uplink. The up and down paths share the turbulence, but
  the two directions differ by the point-ahead angle, so the correction leaves
  the modal decorrelation residual. The satellite terminal needs a transmitter
  for the downlink beam. It has no fields.
- `LaserGuideStar(altitude_m=90e3)` — a ground-launched guide star. NOT
  IMPLEMENTED yet. A guide star at a finite altitude gives focal (cone)
  anisoplanatism, a different effect from the point-ahead angular
  anisoplanatism. The `uplink_budget` raises `NotImplementedError` for this
  source. The default `altitude_m` is a sodium-layer guide star.

### `TerrestrialScenario`

A terrestrial (horizontal-path) link case: a near terminal and a far terminal.
Both ends are on the ground, so the terminals are named for the path ends. The
link is one-way, but the path is reciprocal, so `direction` selects which end
transmits.

The terrestrial `direction` is a different type from the space `direction`,
because "terrestrial" is a channel family, not a tx/rx geometry.

| Field | Type | Unit | Default | Meaning |
|---|---|---|---|---|
| `near` | `Terminal` | — | (required) | The local end of the path. |
| `far` | `Terminal` | — | (required) | The remote end of the path. |
| `direction` | `"forward"` \| `"reverse"` | — | `"forward"` | The transmit end of the path. |
| `channel` | `TerrestrialChannel` | — | `TerrestrialChannel()` | The horizontal propagation channel. |
| `availability_target` | float | — | `0.99` | Target link availability (0-1). |

Role mapping:

| direction | `tx_terminal` | `rx_terminal` |
|---|---|---|
| `forward` | `near` | `far` |
| `reverse` | `far` | `near` |

The channel does not change with the direction, because a horizontal path is
the same in the two directions.

---

## 3. Channels (`olb.scenario`)

A channel is the propagation medium plus its geometry parameters. A channel
holds no terminal hardware.

### `Site`

Ground station location and atmosphere (the propagation medium).

| Field | Type | Unit | Default | Meaning |
|---|---|---|---|---|
| `lat_deg` | float | deg | `-29.0468` | Latitude (TN-2 Kepler OGS default). |
| `lon_deg` | float | deg | `115.3467` | Longitude. |
| `alt_m` | float | m | `269.0` | Station height. |
| `cn2_ground` | float | m^-2/3 | `1.7e-14` | Hufnagel-Valley ground-level Cn2 scale (HV57 A). |
| `wind_rms_m_s` | float | m/s | `21.0` | Bufton wind profile rms. |
| `clear_sky_probability` | float | — | `1.0` | Cloud-free-line-of-sight fraction (0-1). |

### `Channel`

A space propagation channel: the ground site plus the satellite orbit. The space
links read `altitude_m` for the analytic orbit geometry. They build a Cn2(h)
profile from the site.

| Field | Type | Unit | Default | Meaning |
|---|---|---|---|---|
| `site` | `Site` | — | `Site()` | The ground station location and atmosphere (the medium). |
| `altitude_m` | float | m | `600e3` | The satellite altitude, for the analytic orbit geometry. |

### `TerrestrialChannel`

A horizontal (terrestrial) propagation channel: a ground-to-ground path. The
path is horizontal, so there is no orbit altitude and no elevation angle. A
single scalar Cn2 describes the path.

| Field | Type | Unit | Default | Meaning |
|---|---|---|---|---|
| `site` | `Site` | — | `Site()` | The ground atmosphere along the path (the medium). |
| `path_length_m` | float | m | `1e3` | Horizontal path length L. |
| `attenuation_db_per_km` | float | dB/km | `0.5` | Clear-air / haze extinction coefficient. The Beer-Lambert loss is `attenuation_db_per_km * (L / 1000)`. |
| `cn2` | float | m^-2/3 | `1e-14` | Constant refractive-index structure parameter along the path. Read by the terrestrial scintillation Term and the fibre-coupling Terms. |

---

## 4. Geometry (`olb.geometry`)

The models and the budget read two arrays from a geometry object:

- `geom.elevation_deg` — elevation above the horizon [deg].
- `geom.slant_range_m` — ground-station to satellite range [m].

The source of the two arrays does not change the models or the budget. Select
the backend for the task.

### `CircularOrbit`

Analytic circular-orbit geometry over an elevation grid. It is vectorised over
the elevation grid. Use it for parameter sweeps and Monte Carlo.

Constructor: `CircularOrbit(altitude_m, elevation_deg)`.

| Parameter | Type | Unit | Meaning |
|---|---|---|---|
| `altitude_m` | float | m | Orbital altitude above the Earth's surface. |
| `elevation_deg` | float or array | deg | Elevation angle(s) above the horizon. Use an array to sweep. |

Provides to a model:

- `slant_range_m` — the ground-station to satellite range.
- `point_ahead_rad` — the point-ahead angle from the finite speed of light.
- `slew_deg_s` — the apparent line-of-sight slew rate.

### `HorizontalPath`

Horizontal (terrestrial) path geometry: a constant range, no elevation. The
range is the path length. It does not change with any elevation angle. So this
geometry exposes only `slant_range_m`. It has no `elevation_deg`.

Constructor: `HorizontalPath(path_length_m)`.

| Parameter | Type | Unit | Meaning |
|---|---|---|---|
| `path_length_m` | float or array | m | Horizontal path length L. Use an array to sweep the range. |

Provides to a model:

- `slant_range_m` — the path length.

### `TLEPass`

A real satellite pass from a TLE, propagated with skyfield. Use it to replay an
actual pass.

Constructor: `TLEPass(tle_line1, tle_line2, lat_deg, lon_deg, alt_m, times, name="")`.

| Parameter | Type | Unit | Meaning |
|---|---|---|---|
| `tle_line1`, `tle_line2` | str | — | The two-line element set. |
| `lat_deg`, `lon_deg` | float | deg | Ground station geodetic latitude / longitude. |
| `alt_m` | float | m | Ground station height above the WGS84 ellipsoid. |
| `times` | skyfield Time | — | Array of times to sample the pass at. |
| `name` | str | — | Satellite name (cosmetic). Default `""`. |

After construction, `elevation_deg`, `azimuth_deg`, and `slant_range_m` are
arrays over `times`. Elevation is negative when the satellite is below the
horizon. Use the mask `elevation_deg > 0` for the visible pass.

Provides to a model: `elevation_deg`, `azimuth_deg`, `slant_range_m`, and
`times`.

#### `TLEPass.from_window`

Class method. It builds a pass. It samples a time window at a fixed step.

`TLEPass.from_window(tle_line1, tle_line2, lat_deg, lon_deg, alt_m, start_utc, duration_s, step_s=1.0, name="")`.

| Parameter | Type | Unit | Meaning |
|---|---|---|---|
| `start_utc` | tuple | — | `(year, month, day, hour, minute, second)` UTC start. |
| `duration_s` | float | s | Window length. |
| `step_s` | float | s | Sample step. Default `1.0`. |

The remaining parameters match the constructor.
