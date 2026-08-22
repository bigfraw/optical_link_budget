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
| `detector` | `Aperture` or `SMF` or None | — | `None` | The detector front end. None means no receive-coupling Term. |
| `compensation` | list of `TipTilt` or `AO` | — | `[]` | The ordered wavefront-compensation stack. It may be empty. |

Constraints:

- Each `Terminal` owns its own `compensation` list. Two terminals do not share
  one list. The default is a fresh empty list per terminal.
- An empty `compensation` stack leaves the piston-removed turbulence.

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

A detector is one of two front ends. Each carries an optional receive
sensitivity.

#### `Aperture`

Power-in-bucket detector. An aperture (bucket) detector integrates the intensity
over the aperture. It is phase-insensitive. So a compensation stack does not
change its coupling. Use it for parity with the plain downlink budget.

| Field | Type | Unit | Default | Meaning |
|---|---|---|---|---|
| `sensitivity_dbm` | float or None | dBm | `None` | Required received power. None if only losses matter. |

#### `SMF`

Single-mode-fibre detector. A single-mode fibre couples only the field that
matches the fibre mode. The coupling falls when turbulence distorts the
wavefront.

| Field | Type | Unit | Default | Meaning |
|---|---|---|---|---|
| `eta_max` | float | — | `0.8145` | Maximum fibre-to-aperture mode match for an unobscured circular aperture with a flat wavefront. |
| `sensitivity_dbm` | float or None | dBm | `None` | Required received power. None if only losses matter. |

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

The `direction` sets the roles:

| direction | tx_terminal | rx_terminal |
|---|---|---|
| uplink | ground | space |
| downlink | space | ground |
| retro | ground | ground |

### `TerrestrialScenario`

A terrestrial (horizontal-path) link case: a near terminal and a far terminal.
Both ends are on the ground, so the terminals are named for the path ends. The
link is one-way: tx = near (the local end), rx = far (the remote end).

There is no `direction`. "Terrestrial" is a channel family, not a tx/rx
geometry.

| Field | Type | Unit | Default | Meaning |
|---|---|---|---|---|
| `near` | `Terminal` | — | (required) | The local (transmit) end of the path. |
| `far` | `Terminal` | — | (required) | The remote (receive) end of the path. |
| `channel` | `TerrestrialChannel` | — | `TerrestrialChannel()` | The horizontal propagation channel. |
| `availability_target` | float | — | `0.99` | Target link availability (0-1). |

Role mapping: `tx_terminal` is `near`, `rx_terminal` is `far`. The scenario has
no `direction` attribute.

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
| `cn2` | float | m^-2/3 | `1e-14` | Constant refractive-index structure parameter along the path. Read by the (pending) terrestrial scintillation term. |

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
