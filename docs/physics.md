# Physics reference

This document is the physics reference for the `olb` package. Each section
describes one phenomenon. Each section gives what the code models, the governing
relation, the key inputs and outputs, and the source. The sources are the ones
that the code cites. Where the code cites a source but gives no literal DOI
string, this document gives the same citation and names the file.

Conventions in this package:

- Loss is a positive dB. Gain is a negative dB.
- Each Term is a factory `f(scenario, geometry) -> Term`.
- A Term has three faces: a mean, a quantile, and a sampler.
- Every equation in the code cites its source next to the equation.

The phenomena are:

1. Geometry and spreading
2. Transmit beam (Gaussian truncation, deliberate divergence)
3. Atmospheric extinction (slant airmass, horizontal Beer-Lambert)
4. Pointing and jitter fade
5. Turbulence (Cn2, scintillation, coupled-flux MC, Dios beam scintillation,
   Gaussian Fried parameter, AO residual)
6. Fibre coupling (analytic mean, FAST statistical)
7. The fidelity-2 turbulent split step (phase screens, sampling rules,
   reciprocity)

---

## 1. Geometry and spreading

File: `olb/models/geometric.py`

### What the code models

The code models the free-space spreading loss of a Gaussian transmit beam into a
circular receive aperture. The aperture can carry a central circular obscuration,
for example a Cassegrain secondary mirror. The loss is deterministic. It does not
model turbulence.

### Governing relation

The captured fraction is the Gaussian-beam power between the obscuration radius
and the aperture radius:

    eta = exp(-2*b_obs^2 / w(z)^2) - exp(-2*a_rx^2 / w(z)^2)
    loss_db = -10*log10(eta)

Here `a_rx` is the receive aperture radius, `b_obs = obscuration_ratio*a_rx` is
the obscuration radius, and `w(z)` is the 1/e^2 Gaussian beam radius at range z.
The code gets `w(z)` from `free_space_radius` (see Section 2), so a deliberate
transmit divergence widens the beam and adds loss.

### Inputs and outputs

- Inputs: slant range, transmit waist `w0`, receive aperture diameter,
  wavelength, obscuration ratio, and optional transmit divergence.
- Output: geometric loss [dB], positive, over the geometry.

### Assumptions and limits

- Far-field Gaussian beam into a circular aperture.
- Paraxial beam.
- No turbulence.
- Measured validity: see Section 9f. A fidelity-2 space budget takes THIS
  analytic Term, because the wave-optics vacuum solve is grid-noise-limited over
  a slant path.

### Source

The Gaussian-beam power-in-bucket form is stated in the module. The module gives
no literal DOI (source cited in `olb/models/geometric.py`). The beam radius comes
from `free_space_radius` in `olb/beam.py`.

---

## 2. Transmit beam

### 2a. Gaussian truncation efficiency

File: `olb/models/gaussian_efficiency.py`

#### What the code models

A real transmitter sends a Gaussian beam through a finite, possibly obscured,
circular aperture. The aperture clips the wings of the Gaussian. The obscuration
blocks the middle. This is a fixed hardware loss. It does not depend on range.

#### Governing relation

The efficiency is the on-axis far-field gain of the truncated beam, referenced to
the untruncated source Gaussian:

    alpha = a / w_T = (tx_aperture_m / 2) / tx_waist_m
    eta   = (exp(-alpha^2) - exp(-alpha^2 * Cr^2))^2
    loss_db = -10*log10(eta)

Here `a` is the aperture radius, `w_T` is the Gaussian waist (1/e^2 radius) at the
aperture, and `Cr` is the linear central-obscuration ratio. The efficiency goes to
1 when the aperture is much wider than the beam. It goes to 0 when the aperture is
much smaller than the beam. The classic optimal truncation is near alpha = 1.12.

This is the corrected antenna-gain form. It has no `2/alpha^2` prefactor. That
prefactor double-counts a normalisation that the untruncated-source reference
already carries. The code validated the removal against a numerical Fraunhofer
propagation of a truncated Gaussian (the `tn2_kepler` `test_gauss_prop`).

#### The top-hat correction

A retroreflector reflects the flat wavefront that fills its aperture, so its return
is a top-hat, not a Gaussian. The correction that converts the
Gaussian(w0 = aperture/2) geometric model to a top-hat is:

    correction_db = 10*log10( 2 / (1 - Cr^2) )

For an unobscured corner-cube (Cr = 0) the correction is +3.01 dB.

#### Inputs and outputs

- Inputs: transmit aperture diameter, transmit waist `w_T`, obscuration ratio.
  A bistatic beam director overrides the aperture and obscuration.
- Output: truncation loss [dB], positive.

#### Assumptions and limits

- On-axis far-field gain of a truncated Gaussian, referenced to the untruncated
  source.
- Paraxial beam. No turbulence.
- The far-field form breaks in the near field. A hard-truncated beam that a
  receiver sees inside the Rayleigh range `zR = pi*w_T^2/lambda` gets an on-axis
  intensity that oscillates with range (Fresnel edge diffraction). The
  `olb/waveoptics/` fidelity-2 layer is the no-turbulence check for such a link.
  See [architecture.md](architecture.md) Section 1.
- The Term now FLAGS that break. `tx_gaussian_efficiency_term(scenario,
  geometry)` reads `geometry.slant_range_m` for the validity test only; the loss
  stays range-independent. When `alpha` is below the module constant
  `TRUNCATION_NEAR_FIELD_ALPHA = 1.5` AND a range sits inside `zR`, the Term sets
  an ACTIVE assumptions violation, so `Budget.check()` reports it. A widely open
  aperture (`alpha` at or above 1.5) keeps almost no power in the clipped wings,
  so the near-field ripple is negligible. `geometry=None` skips the test, which
  is the space-budget case: a space link is always far field.
- The error is NOT conservative. Unlike the geometric spreading Term (exact at
  every range through `gaussz`) this Term does not self-correct, and unlike the
  single-mode-fibre `eta_max` of Section 6a the true value can sit above OR below
  it. Verify a flagged link with a fidelity-2 no-turbulence field propagation.
- Measured validity: see Sections 9c (the obscured launch pupil) and 9f (the
  vacuum geometric loss).

#### Source

The corrected antenna-gain form is stated in the module and is validated against
the `tn2_kepler` propagation test. The module gives no literal DOI (source cited
in `olb/models/gaussian_efficiency.py`).

### 2b. Deliberate transmit divergence (virtual waist)

File: `olb/beam.py`

#### What the code models

The code sends a wider beam on purpose, for example to relax the pointing budget.
It recasts a deliberately diverged transmitter as an ordinary Gaussian beam that
starts from a virtual waist behind the aperture. Then the ordinary Gaussian
machinery does the rest. No separate curvature bookkeeping is needed.

#### Governing relation

For a transmitter of aperture radius `w0` and far-field half-angle divergence
`theta`:

    w_v = lambda / (pi * theta)                 virtual waist radius
    d   = zR(w_v) * sqrt((w0/w_v)^2 - 1)         distance behind the aperture
    w(z) = gaussz(w_v, d + z)                     free-space radius at range z

`theta` cannot be below the aperture diffraction limit `lambda/(pi*w0)`. A `theta`
equal to that limit is the collimated case, which returns `(w0, 0)`. `None` also
means collimated. The collimated result reduces exactly to `gaussz(w0, z)`.

#### Inputs and outputs

- Inputs: aperture radius `w0`, divergence half-angle, wavelength, range z.
- Output: the turbulence-free beam radius at range z.

#### Assumptions and limits

- The divergence must be at least the diffraction limit. A sub-diffraction
  request raises a ValueError.
- Pure paraxial Gaussian optics. No turbulence.

#### Source

The virtual-waist recast is ported from `tn2_kepler` `fso_spot_size.py`
(`virtual_waist`, `free_space_radius`, `spot_sizes`). The module gives no literal
DOI (source cited in `olb/beam.py`).

---

## 3. Atmospheric extinction

File: `olb/models/extinction.py`

### 3a. Slant airmass extinction (space link)

#### What the code models

The code models the clear-sky molecular and aerosol extinction along a slant path.
It uses one zenith optical depth that the airmass scales. It is not a MODTRAN
line-by-line model. It is a one-parameter slant attenuation.

#### Governing relation

    airmass(elevation) = 1 / sin(elevation)
    T       = exp(-tau_zenith * airmass(elevation))
    loss_db = -10*log10(T) = (10/ln10) * tau_zenith * airmass(elevation)

The loss is linear in the airmass. `tau_zenith` combines molecular (Rayleigh) and
aerosol extinction into one clear-sky number. The default is
`DEFAULT_TAU_ZENITH = 0.05`. It is the near-IR clear-sky zenith optical depth at
1550 nm for a good, dry, high site. It gives a zenith transmittance of
exp(-0.05) = 0.95, about 0.22 dB. At 1550 nm Rayleigh scattering is negligible
(tau ~ 0.002), so aerosol dominates.

#### Inputs and outputs

- Inputs: elevation angle [deg], zenith optical depth.
- Output: extinction loss [dB], positive, larger at low elevation.

#### Assumptions and limits

- Clear sky. One zenith optical depth. Plane-parallel airmass.
- The airmass model diverges at the horizon. Do not use elevation 0.
- The code flags an elevation below 5 deg, which breaks the plane-parallel model.

### 3b. Horizontal Beer-Lambert extinction (terrestrial link)

#### What the code models

The code models the extinction over a horizontal terrestrial path. A horizontal
path has a constant extinction per unit length, so the loss is linear in the path
length.

#### Governing relation

    loss_db = attenuation_db_per_km * (path_length_m / 1000)

The coefficient is quoted directly in dB/km, so no logarithm is needed.

#### Inputs and outputs

- Inputs: horizontal path length [m], extinction coefficient [dB/km].
- Output: extinction loss [dB], positive.

#### Assumptions and limits

- Horizontal path. One dB-per-km coefficient, constant along the path.
- The coefficient value is weather- and visibility-dependent (fog, haze, rain).
  The user sets it per site.

### Source

The airmass model is borrowed from `fso_spot_size.airmass` in the sibling TN-2
analysis repo. The Beer-Lambert forms are stated in the module. The module gives
no literal DOI (source cited in `olb/models/extinction.py`).

---

## 4. Pointing and jitter fade

File: `olb/models/pointing.py`

### What the code models

The system aims a Gaussian transmit beam with 2-D Gaussian pointing jitter of
1-sigma angle `sigma_theta`. The code models the fade that this jitter causes at
the receiver.

### Governing relation

The radial pointing displacement at the receiver is `r = sigma_r * |unit 2-D
Gaussian|`, with `sigma_r = sigma_theta * range`. For the small-aperture, on-axis
Gaussian approximation the collected-power fraction against boresight is:

    h(r) = exp(-2*r^2 / w_z^2)
    loss_db = (20/ln10) * r^2 / w_z^2

The two jitter axes are i.i.d. Gaussian, so `r^2` is exponential. The loss in dB
then has an exponential distribution:

    loss_db ~ Exponential(mean = (20/ln10) * 2 * sigma_r^2 / w_z^2)

This gives a closed-form mean, quantile, and sampler. The quantile is the inverse
exponential CDF `quantile(p) = -mean * ln(1 - p)`.

### Inputs and outputs

- Inputs: slant range, transmit waist, jitter angle `sigma_theta`, wavelength,
  optional divergence.
- Output: a pointing-fade Term. Zero jitter gives a deterministic 0 dB Term. A
  nonzero jitter gives an exponential mean, quantile, and sampler.

A wider (diverged) beam gives less pointing loss than a collimated beam of the
same `w0`.

### Assumptions and limits

- Small receive aperture relative to the beam.
- On-axis Gaussian beam. No aperture averaging of the fade.
- 2-D Gaussian jitter.
- The code flags a receive aperture that is not small relative to the beam
  (`a_rx > 0.5*w_z`), where the on-axis approximation is weak.

### Note on the uplink path

The uplink coupled-flux Monte Carlo (Section 5c) does NOT stack a separate
pointing Term. It folds the mechanical jitter into the beam-wander variance
`beta2`. Stacking both would double-count the jitter displacement. See Section 5c.

### Source

The exponential-fade derivation is stated in the module. The module gives no
literal DOI (source cited in `olb/models/pointing.py`).

---

## 5. Turbulence

### 5a. Cn2 profiles

File: `olb/turbulence/profiles.py`

#### What the code models

The code builds the default zenith Cn2(h) profile from the site parameters. It
also holds the default turbulence altitude grid `DEFAULT_HS`, which is 20 points
log-spaced from 1 m to 20 km.

#### Governing relation

    Cn2(h) = get_c2n(hs, site.wind_rms_m_s, site.cn2_ground)

`get_c2n` is the Hufnagel-Valley model, defined in `olb.turbulence.profiles`. It
reads the site RMS wind and the ground-level Cn2. Use this profile when the
optional `fast` package is not available. The `fast` HV57 path fails without
that package.

#### Inputs and outputs

- Inputs: the site (RMS wind, ground Cn2), the altitude grid.
- Output: the zenith Cn2(h) profile on the grid.

#### Source

The profile builder is `get_c2n` (the Hufnagel-Valley model, Andrews and
Phillips, DOI 10.1117/3.626196, Ch. 12, Eq. (1); the formula and its citation
are in `olb/turbulence/profiles.py`).

### 5b. Scintillation index and aperture averaging (downlink plane wave)

File: `olb/turbulence/plane_wave_scintillation.py`

#### What the code models

The satellite is far away, so the downlink source is a plane wave at the top of
the atmosphere. The code gives the plane-wave scintillation index and the
aperture-averaging integral for that plane wave down to the ground aperture.

#### Governing relations

Plane-wave point scintillation index, integrated over the Cn2 slant path:

    sigma2_I = 2.25 * k^(7/6) * sec(zeta)^(11/6) * INT Cn2(h) h^(5/6) dh

Here `k = 2*pi/lambda`, `sec(zeta) = 1/sin(elevation)`, and `h` is the height
above the ground station.

Aperture-averaged flux index, the distributed-path Rytov double integral over
height h and spatial wavenumber kappa:

    sigma2_P = 8*pi^2 * k^2 * sec(zeta)
        * INT_h Cn2(h) [ INT_kappa kappa * 0.033 * kappa^(-11/3)
        * (1 - cos(kappa^2 * h * sec(zeta) / k))
        * (2*J1(kappa*D/2) / (kappa*D/2))^2 dkappa ] dh

The aperture-averaging factor is `A = sigma2_P / sigma2_I`. It obeys
`0 < A <= 1`, `A -> 1` as `D -> 0`, and the large-aperture asymptote `A ~ D^(-7/3)`.

#### Closed-form single-path forms

The module also gives closed-form algebraic approximations for one path length L
and one scalar Cn2. They hold for any turbulence strength. They give a fast answer
without an integral. Use them for a horizontal path or a quick estimate:

    Rytov std:       sigma_1 = ( 1.23 Cn2 k^(7/6) L^(11/6) )^0.5
    coherence radius: rho_c  = ( 1.46 Cn2 k^2 L )^(-3/5)
    point index:     sigma_I^2 = exp[ 0.54 s^2/(1+1.22 s^(12/5))^(7/6)
                        + 0.509 s^2/(1+0.69 s^(12/5))^(5/6) ] - 1
    aperture index:  sigma_I^2(D) = exp[ 0.49 s^2/(1+0.65 d^2+1.11 s^(12/5))^(7/6)
                        + 0.51 s^2 (1+0.69 s^(12/5))^(-5/6)
                          / (1+0.90 d^2+0.62 d^2 s^(12/5)) ] - 1

with `s = sigma_1` and the aperture parameter `d = ( k D^2 / (4 L) )^0.5`. The
module also gives the weak, weak-large-inner-scale, and strong aperture-averaging
factors.

#### Inputs and outputs

- Inputs: elevation, wavelength, height grid, Cn2 profile, receive aperture
  diameter D.
- Output: the scintillation index and the aperture-averaging factor.

#### Assumptions and limits

- Plane wave, weak fluctuation, isotropic turbulence.
- The lognormal irradiance PDF is trusted for `sigma2_I` below about
  `LOGNORMAL_PDF_LIMIT = 0.25`. Above it, focusing and saturation make the
  lognormal model depart from data. That 0.25 is a HOUSE RULE, 4 times stricter
  than the book bound `sigma_R^2 < 1` (Ch. 5, Eq. (15), printed p. 140). It is
  kept, because Ch. 11, Sec. 11.3, printed p. 451, says the lognormal tail is
  optimistic, and the fade faces report that tail.
- `downlink_scintillation_term(..., model="auto")` reads the point index and
  switches at that limit: below it the lognormal Term, at or above it the
  gamma-gamma Term of Section 5h. The gamma-gamma chain holds at every
  fluctuation strength (Ch. 12, Eq. (40), printed p. 497). It models a POINT
  receiver, because the book gives no aperture-averaged downlink index in that
  regime.
- The aperture filter `(2*J1(x)/x)^2` assumes a uniform circular aperture with no
  central obscuration. An annular aperture is not modelled yet.
- The Kolmogorov spectrum has no inner scale and no outer scale.
- Measured validity: see Section 9e. The lognormal SHAPE is certified in the weak
  band, and the weak aperture-averaging factor `A` over-averages by 1.4 to 2.9
  times over `1 <= D/rho_0 <= 8`.

#### Source

- Andrews and Phillips, Laser Beam Propagation through Random Media, 2nd ed.
  (2005), DOI 10.1117/3.626196: Ch. 12, Eq. (38), printed p. 495 (repeated as
  Ch. 12, Eq. (92), printed p. 522) for the plane-wave index; Ch. 10 for the
  aperture-averaging filter and the closed-form factors; Ch. 5, 6, 9 for the
  single-path forms.
- Churnside, Applied Optics 30 (1991) 1982, for the strong-turbulence
  aperture-averaging factor.

### 5c. Uplink coupled-flux Monte Carlo

File: `olb/turbulence/uplink_flux.py`

#### What the code models

The code runs the short uplink loop of the Dios et al. coupled-flux Monte Carlo
for one elevation. It models beam wander plus scintillation for a beam that
propagates up to a LEO satellite. It rescales the flux to the free-space baseline,
so the result is additive with the geometric Term.

#### Governing points

- The coupled flux normalises the on-axis irradiance to the beam's own short-term
  waist `w_st`. The short-term waist already includes the turbulence spread
  (`w_st = sqrt(w_free^2 + turbulence_term)`). This self-normalisation removes the
  beam-broadening loss. For a 600 km, ~1 m-waist uplink the removed loss is about
  10 dB. So the code rescales the flux samples by `(w_free / w_st)^2` to put the
  samples back on the free-space baseline. Then the per-sample dB loss is
  `-10*log10(flux)`.
- Divergence enters everywhere the beam geometry matters: the broadening baseline,
  the short and long-term waists (through the `w_free` override), and the
  scintillation index (through the receiver-plane Lambda and Theta of the diverged
  beam). A diverged beam is larger and more spherical-wave-like, so it
  scintillates less.
- Pointing jitter folds into the wander variance. `beta2` is the total 2-D wander
  variance. A per-axis 1-sigma jitter angle `sigma_theta` gives a per-axis
  displacement `sigma_r = sigma_theta*L`, which adds `2*sigma_r^2` to the 2-D
  total: `beta2 = beta2 + 2*(sigma_theta*L)^2`. The combined offset feeds both the
  Gaussian power falloff and the off-axis scintillation. This replaces a standalone
  pointing Term on the uplink. Adding both would double-count the jitter.

#### Inputs and outputs

- Inputs: transmit waist `w0`, elevation, slant range, wavelength, height grid,
  Cn2 profile, sample counts, optional divergence, optional jitter angle.
- Output: rescaled flux samples, the short and long-term waists, the coherence
  diameter, the mean log-amplitude variance, and a weak-fluctuation flag.

#### Assumptions and limits

- The Rytov model is a weak-fluctuation model.
- When the mean log-amplitude variance `sigma2_x` exceeds the Dios reliability
  edge `UPLINK_SIGMA2X_LIMIT = 0.6`, the scintillation approaches saturation and
  the numbers are not trustworthy. The code carries the flag and warns.
  `_flux_result` now OWNS this limit as a traced hard-tier check (the uplink Term
  inherits its violation), while the two `warnings.warn` lines stay verbatim.
- The launch beam is an untruncated Gaussian of waist `w0`. The code models no
  launch aperture and no central obscuration, so the fade does not change with an
  obscured pupil. `uplink_turbulence_term` flags a set obscuration in
  `budget.check()`. The MEAN loss from a central obscuration is separate and IS
  carried: the launch-truncation Term (`tx_gaussian_efficiency_term`, Section 2a)
  reads the obscuration and matches the wave-optics far-field to about 2 dB. The
  size of the obscuration effect on the FADE is now MEASURED against fidelity 2
  alone: measured validity, see Section 9c.
- Measured validity of the fade itself: see Section 9a (a filled launch) and
  Section 9b (the vendored kernels and a slant-coordinate defect).
- The coupled-flux MC needs the `fast` package to build the HV57 Cn2 profile, or
  an explicit `cn2_profile`.

#### Source

Dios et al., Applied Optics 43 (2004) 3866, for the coupled-flux and wander-offset
mechanism. The kernels are vendored in `olb.turbulence.coupled_flux`.

### 5d. Dios on-axis and off-axis beam scintillation

File: `olb/turbulence/beam_wave_scintillation.py`

#### What the code models

The code gives the scintillation index of a collimated Gaussian beam that
propagates up through turbulence to a receiver. It integrates a real Cn2(h) slant
profile. It returns the flux scintillation index on axis and off axis. This is the
analytic twin of the coupled-flux MC.

#### Governing relations

The total index at a point a radius r off the beam axis is the sum of the on-axis
(longitudinal) term and the radial (off-axis) term:

    sigma2_I(r, L) = sigma2_I(0, L) + sigma2_Ir(r, L)          (Dios Eq. 13)
    A(z) = (Lambda L / k) ((L-z)/L)^2                          (Eq. 17)
    B(z) = (L / k) ((L-z)/L) (Theta + (1-Theta) z/L)           (Eq. 18)
    sigma2_I(0,L) = 4 pi^2 k^2 Gamma(-5/6) 0.033
        * INT Cn2(z) [ A^(5/6) - (A^2+B^2)^(5/12)
                       * cos( (5/6) arctan(B/A) ) ] dz         (Eq. 16)
    sigma2_Ir(r,L) = 4 pi^2 k^2 Gamma(-5/6) 0.033
        * ( 1F1(-5/6, 1, 2 r^2 / W^2(L)) - 1 )
        * INT Cn2(z) A^(5/6) dz                                (Eq. 20)

Here z runs from the ground (z = 0) to the receiver (z = L). The weight `(L-z)/L`
gives the near-ground turbulence the full weight and the near-receiver turbulence
almost none. That is the uplink weighting. In the plane-wave limit (large w0) the
on-axis term returns the plane-wave Rytov variance. In the spherical-wave limit
(small w0) it returns 0.40 of that. The beam-wander model feeds `r = beta`, the
wander offset.

The geometric range L is separate from the turbulence grid. For a satellite uplink
pass `path_length_m` = the full slant range, so the weight stays near 1 across the
thin turbulence layer (the far-field limit).

#### Inputs and outputs

- Inputs: height grid, Cn2 profile, transmit waist `w0`, wavelength, elevation,
  phase-front curvature `f0`, optional slant range.
- Output: the on-axis, radial, or total scintillation index.

#### Assumptions and limits

- Weak-to-moderate turbulence. The model has no saturation. Measured validity:
  see Sections 9a (the unsaturated off-axis term at a filled launch) and 9e (the
  on-axis index against the field).
- Dios reports good agreement with a split-step reference up to
  `sigma2_chi ~ 0.6`. Above that the true index saturates and the model overshoots.
- `on_axis_scintillation_index` now OWNS this limit through a runtime check
  (`@assumes`, the `DIOS_RELIABILITY` constraint). The check fires when the index
  it returns leaves the weak regime (`hard_limit = 4 * UPLINK_SIGMA2X_LIMIT = 2.4`
  on the Rytov axis). The terrestrial scintillation Term therefore inherits its
  hard-tier violation from THIS function: the flag migrated from the factory's old
  `sigma_R^2 >= 1.0` test to the beam-wave index axis (`sigma_I^2 >= 2.4`). The two
  coincide on the tested strong and long-path triggers. A narrow band between them
  can now read `ok` where the old axis flagged, which is defensible; the tighter
  lognormal-PDF house rule (`sigma_I^2 < 0.25`, a PDF-shape flag the factory keeps)
  backstops the common cases.

#### Note (convergence, not duplication)

This analytic Dios path and the coupled-flux MC (Section 5c) use the same Dios
equations. The jitter-into-`r` correction is present only in the MC path. A code
comment marks this as debt: the terrestrial scintillation slot must converge on one
implementation with the jitter folded into r, not make a third copy.

#### Source

Dios et al., Applied Optics 43 (2004) 3866, Eqs. (13)-(20).

### 5e. Gaussian-beam Fried parameter

File: `olb/turbulence/gaussian_fried.py`

#### What the code models

The code gives the Fried parameter r0 of a collimated Gaussian beam through
turbulence. It has a closed-form single-path form and a profile-integral form.
The profile form serves an uplink, a downlink, or a horizontal terrestrial link.

#### Governing relations

Constant-Cn2 single path:

    r0_gauss = 2.1 * rho0_e * rho_pl
    rho_pl   = ( 1.46 Cn2 k^2 z )^(-3/5)              plane-wave coherence radius
    Lambda0  = 2 z / (k w0^2),  Theta0 = 1           collimated input
    Lambda   = Lambda0 / (Lambda0^2 + Theta0^2)
    Theta    = Theta0  / (Theta0^2 + Lambda0^2)
    rho0_e   = ( 8 / (3 (a + 0.62 Lambda_e^(11/6))) )^(3/5)

The effective parameters `Theta_e`, `Lambda_e` hold the strong-turbulence beam
spread through the Rytov std `sigma_R`. The plane-wave Fried parameter is
`r0_pl = 2.1*rho_pl`, which equals the standard `(0.423 k^2 Cn2 L)^(-3/5)`. The
spherical-wave form on a constant-Cn2 horizontal path is
`r0_sph = (8/3)^(3/5) * r0_pl ~ 1.80 * r0_pl`.

Profile integral (nonuniform Cn2, with the wave-structure-function weight):

    xi     = (L - z) / L                              1 at tx, 0 at rx
    mu1    = INT Cn2 (Theta + Theta_bar xi)^(5/3) dh
    mu2    = INT Cn2 xi^(5/3) dh
    rho0   = ( 1.46 k^2 sec(zeta) (mu1 + 0.62 Lambda^(11/6) mu2) )^(-3/5)
    r0     = 2.1 * rho0

In the plane-wave limit (Theta -> 1) the weight is 1. In the spherical-wave limit
(Theta -> 0) the weight is `xi^(5/3)`, which matches Dios Eq. (3). The `path`
argument sets the transmitter end and the weight direction: `uplink`, `downlink`,
or `terrestrial`.

#### Inputs and outputs

- Inputs: path length or height grid, transmit waist `w0`, Cn2 (scalar or
  profile), wavelength, path type, elevation, curvature, optional slant range.
- Output: the Gaussian-beam Fried parameter r0 [m].

#### Assumptions and limits

> **FLAG — collimated-beam limitation.** The single-path form
> `gaussian_fried_parameter` holds ONLY for a collimated Gaussian beam. It fixes
> the input curvature `Theta0 = 1` through the module constant
> `COLLIMATED_THETA0 = 1.0`. It takes no curvature input, so a focused or a
> diverged beam gives a wrong r0. The profile form
> `gaussian_fried_parameter_profile` is NOT limited this way: it takes a
> phase-front radius of curvature `f0` and computes `Theta0 = 1 - L/f0`. Its
> default `f0 = inf` still gives a collimated beam. The one call site, the
> terrestrial single-mode-fibre coupling in `olb/models/coupling/terrestrial.py`,
> now passes `f0`: it reads the launch curvature from the transmitter divergence
> through `olb.beam.launch_curvature`, so a deliberately diverged beam gets its
> own r0 (olb Gap 3, closed 2026-08-27). `launch_curvature` shares one
> implementation with the Dios scintillation feed in
> `olb/turbulence/uplink_flux.py`, so the geometric, the scintillation, and the
> Fried-parameter Terms read the SAME f0.
>
> STILL OPEN: the single-path form `gaussian_fried_parameter` keeps the
> collimated signature. To generalise it, add a curvature argument and thread it
> through `output_beam_params`, the way the profile form already does. The
> profile form is the one the budgets use, so this is a tidy-up, not a gap.

- The profile form assumes a Gaussian beam in the weak turbulence regime with a
  Kolmogorov spectrum (no inner scale, no outer scale). It holds to about a Rytov
  variance below 1 (Dios reports good agreement to `sigma_chi^2 ~ 0.6`).
- The profile form uses the free-space (diffractive) Theta and Lambda, not the
  strong-turbulence effective parameters. So it does not capture the
  turbulence-driven beam spread in strong turbulence. This is a deliberate
  deferral to match the Dios weak regime.

#### Source

- Andrews and Phillips, 2nd ed. (2005): Ch. 6 (Gaussian-beam coherence, beam
  parameters), Ch. 9 (effective beam parameters), Sec. 12.4.1 (wave structure
  function).
- Dios et al., Applied Optics 43 (2004) 3866, Eq. (3), for the uplink weight.

### 5f. AO and tip-tilt residual wavefront

File: `olb/turbulence/ao.py`

#### What the code models

The satellite is far away, so the downlink source is a plane wave. The code gives
the plane-wave Fried parameter for that plane wave, and the residual phase
variance that a wavefront-compensation stack leaves.

#### Governing relations

Plane-wave downlink Fried parameter:

    r0 = ( 0.423 * k^2 * airmass * INT Cn2(h) dh )^(-3/5)

with `airmass = 1/sin(elevation)`.

Residual phase variance over an aperture of diameter D after Zernike correction:

    sigma^2 = c * (D/r0)^(5/3)     [rad^2]

The Noll coefficient c depends on the correction:

    piston removed only        c = 1.0299    (no correction)
    first 3 Zernikes removed    c = 0.134     (tip-tilt)
    first J Zernikes removed    c = 0.2944 * J^(-sqrt(3)/2)   (large-J AO)

The AO form is the large-order asymptotic of the Noll residual series. The code
also carries a fidelity-1 hook: a high-pass Kolmogorov residual phase PSD above the
AO correction cutoff `f_c = sqrt(n_modes) / (2 D)`. Fidelity 0 does not evaluate
the PSD.

#### Inputs and outputs

- Inputs: Cn2 profile, height grid, wavelength, elevation, the compensation stack
  (TipTilt, AO), and the aperture diameter.
- Output: the plane-wave r0, and a `ResidualWavefront` with the residual variance
  and the PSD hook. The best-correcting stage sets the residual.

#### Assumptions and limits

- Plane-wave downlink source (the 0.423 constant).
- The spherical-wave uplink uses a different weight, which lives in the
  coupled-flux and Gaussian-Fried modules.
- The AO coefficient is the large-order asymptotic.

#### Source

- Fried (1966), and Andrews and Phillips, 2nd ed. (2005), Ch. 3, for the Fried
  parameter.
- R. J. Noll, "Zernike polynomials and atmospheric turbulence," J. Opt. Soc. Am.
  66(3), 207-211 (1976), for the residual variance coefficients.
- Kolmogorov phase PSD, Andrews and Phillips, 2nd ed. (2005), Ch. 3, for the
  fidelity-1 PSD hook.

---

### 5g. Angular anisoplanatism (point-ahead pre-compensation)

File: `olb/turbulence/anisoplanatism.py`

#### What the code models

Two wavefronts come to one aperture from two directions. The turbulence gives
them different phase. The code gives the mean-square phase difference between
them. This is the angular anisoplanatic error. It sets the limit of a system
that senses the turbulence on one source and corrects a second source at a
small angle. The main use is uplink pre-compensation. The ground terminal
senses the downlink beam and pre-corrects the uplink. The uplink goes to the
point-ahead direction, so the two directions differ by the point-ahead angle.

#### Governing relations

Classical isoplanatic angle (Stone Eq. 27):

    theta0 = ( C1 * k0^2 * INT Cn2(h) S^(5/3) dh )^(-3/5),   S = h * airmass
    C1 = 2 (2 pi)^(8/3) C_A |HJ1(8/3,0,1)| = 2.914381

Finite-aperture phase variance (Stone Eq. 29):

    sigma^2 = 2 (2 pi)^(8/3) C_A k0^2 R^(5/3) airmass * INT Cn2(h) I(S theta/R) dh
    I(b)    = INT_0^inf u^(-8/3) (1 - J0(b u)) M(u) du            (Eq. 36)

with `R = D/2`. The modal weight `M(u)` selects the Zernike orders. Each radial
order n has the weight (Stone Eq. A11):

    p_n(u) = 4 (n+1)^2 ( J_{n+1}(u) / u )^2

The order n=0 is the piston. The order n=1 is the two tilts. The weights of all
orders sum to 1. So the piston-and-tilt-removed weight is `M(u) = 1 - p0 - p1`.
This is the residual of a perfect, infinite-order correction.

Each order's variance is the DECORRELATION residual between the two directions:

    <[a_n(I) - a_n(II)]^2> = 2 sigma_n^2 ( 1 - rho_n )

where rho_n is the correlation of order n across the angle. An adaptive-optics
system removes the correlated part. The decorrelated part stays. So the error
of a system that corrects the orders 2..max_order is the band weight (Stone
Eq. 43, the frequency-restricted variance):

    M(u) = p2(u) + p3(u) + ... + p_{max_order}(u)

The band grows with `max_order`. Each added order brings its own decorrelation
residual, small for a well-correlated low order and up to twice the mode
variance for a fully decorrelated one. It goes up to the infinite-order value.
This is NOT a penalty for correcting; it is the part of the turbulence that the
two directions do not share. The classical `(theta/theta0)^(5/3)` law (Stone
Eqs. 1, 26) keeps the piston and overpredicts the error, up to about one order
of magnitude for a small aperture (Stone Fig. 1).

#### Inputs and outputs

- Inputs: the aperture diameter D, the angle theta, the Cn2 profile, the height
  grid, the wavelength, `remove` (`"none"`, `"piston"`, or `"piston_tilt"`),
  `max_order` (None for all orders, or an integer for a finite AO), and the
  elevation.
- Output: the phase variance sigma^2 [rad^2].
- `max_radial_order(n_zernike_modes)` turns a Noll mode count into the highest
  complete radial order, so an `AO(n_modes)` stage maps to a `max_order`.
- The code integrates the Bessel kernel directly. It does not use the paper's
  hypergeometric 3F2 series (Stone Eqs. 31-32).

#### Assumptions and limits

- Kolmogorov spectrum, no inner or outer scale.
- The pure angular case: both sources are at infinity and the two beams share
  the aperture radius.
- The piston is always removed as optically harmless. The tilt is removed for a
  beam that a separate tracking loop points.
- This is a phase quantity. It carries no amplitude scintillation. The uplink
  pre-compensation budget that uses this Term (`uplink_budget(fidelity=0)` with a
  `DownlinkBeacon` source) is therefore phase-only and mean-only: no
  scintillation and no fade. That is a recorded decision (2026-08-27, backlog
  0-W1): no trustworthy analytic form exists for the scintillation of a
  pre-compensated beam. The model of record is the fidelity-1 FAST route:
  `uplink_fast_term` (olb/models/fast.py) overlaps the ground-pupil
  field with the adaptive-optics residual phase (the point-ahead decorrelation
  included) and a log-normal log-amplitude, and reads the uplink flux by
  reciprocity (Shapiro, DOI 10.1364/JOSA.61.000492; Farley and others,
  DOI 10.1364/OE.458659). `uplink_budget(fidelity=1)` (the default for a
  pre-compensated scenario) consumes it. See `api-budget.md`.

#### Source

- J. Stone, P. H. Hu, S. P. Mills and S. Ma, "Anisoplanatic effects in
  finite-aperture optical systems," J. Opt. Soc. Am. A 11(1), 347-357 (1994).
  DOI: 10.1364/JOSAA.11.000347.
- R. J. Noll, J. Opt. Soc. Am. 66(3), 207-211 (1976), DOI 10.1364/JOSA.66.000207,
  for the Zernike mode count in `max_radial_order`.

---

### 5h. The Andrews foundation layer

Package: `olb/turbulence/andrews/`. Adapter: `olb/models/fade.py`.

#### What the code models

The nine modules of `olb/turbulence/andrews/` hold the propagation physics of
one book:

- L. C. Andrews and R. L. Phillips, *Laser Beam Propagation through Random
  Media*, 2nd ed., SPIE Press (2005). DOI: 10.1117/3.626196.

Each function names its chapter, its equation number and its printed page. The
package holds physics only. No module in it imports a scenario, a terminal, a
Term or a link, and no function in it returns decibels. So the package is the
one place to read a book equation, and the link layer above it stays thin.

| Module | What it gives | Book |
| --- | --- | --- |
| `beam.py` | Gaussian-beam parameters Theta0, Lambda0, Theta, Lambda, W, for ANY input curvature f0, and the strong-turbulence effective pair | Ch. 4, Eqs. (33), (37), (44), (45); Ch. 7, Eq. (58) |
| `spectra.py` | The five refractive-index spectra: Kolmogorov, Tatarskii, von Karman, exponential, modified atmospheric | Ch. 3, Eqs. (18) to (23) |
| `structure.py` | Wave structure function, coherence radius, Fried parameter, angle-of-arrival variance, rms image jitter | Ch. 6, Secs. 6.4 and 6.5; App. III, Tables I to VI |
| `scintillation.py` | Rytov variance, the two-scale weak index, and the large-scale and small-scale log variances that feed the gamma-gamma parameters | Ch. 8, Eq. (20); Ch. 9, Eqs. (41), (46), (48), (97), (101) |
| `wander.py` | Beam-wander variance, pointing-error variance, long-term and short-term beam radius, and their slant-path forms | Ch. 6, Eqs. (88) and (100); Ch. 12, Eqs. (50) to (53) |
| `aperture.py` | Aperture averaging of the irradiance flux, weak and all-regime, on a homogeneous path | Ch. 10, Sec. 10.3 |
| `distributions.py` | The irradiance PDFs (lognormal, gamma-gamma, K, lognormal-Rician), each with parameters, mean log, CDF, quantile and sampler; the probability of fade, the fade rate and the mean fade time | Ch. 9, Secs. 9.9 to 9.11; Ch. 11, Sec. 11.3 |
| `temporal.py` | Temporal spectra of the irradiance, the quasi-frequency, the Greenwood frequency and the coherence time | Ch. 8, Sec. 8.5; Ch. 9, Sec. 9.8; Ch. 14, Eqs. (38) and (39) |
| `paths.py` | The slant path and the satellite link: the H-V profile, the Bufton wind, the path moments mu_0 to mu_3, the uplink and downlink scintillation indices, the coherence radius, the isoplanatic angle and the point-ahead angle | Ch. 12 |

#### One distribution drives all three Term faces

A Term has three faces: `mean_db`, `quantile(p)` and `sampler(n, rng)`. Before
the foundation layer, each link module built the three faces by hand from a
lognormal. Now one irradiance model gives all three, through one adapter,
`olb/models/fade.py`:

    loss_db     = -10 log10(I),   with E[I] = 1
    mean_db     = -(10/ln10) E[ln I]
    quantile(p) = -10 log10( I(1 - p) )
    sampler     = -10 log10( draws of I )

`irradiance_fade_term(name, category, mean_log=..., quantile=..., rvs=...)`
takes the three faces of any model in `distributions.py` and returns the Term.
So a new distribution needs no new dB code. The loss at availability `p` reads
the `1 - p` quantile of the irradiance, because a deeper loss is a smaller
irradiance.

The downlink uses that adapter for the gamma-gamma Term. See Section 5b for the
selector: `downlink_scintillation_term(..., model="auto")` returns the lognormal
Term while the point index stays below the house limit 0.25, and the gamma-gamma
Term at or above it. The gamma-gamma chain is valid at every fluctuation
strength (Ch. 12, Eq. (40), printed p. 497), so the early switch costs no
validity. Its scintillation index `1/alpha + 1/beta + 1/(alpha beta)` (Ch. 9,
Eq. (139), printed p. 371) is identically `exp(sigma_lnX^2 + sigma_lnY^2) - 1`,
which is the book's own weak-to-strong downlink index. So the Term composes the
book, it does not re-derive it.

#### The delegation pattern

An older module keeps its NAME and its SIGNATURE, and its body calls the
foundation layer. Nothing above changed, and the physics has one home:

| Old home | New home |
| --- | --- |
| `plane_wave_scintillation.coherence_radius` | `andrews.structure.coherence_radius` |
| `plane_wave_scintillation.aperture_averaged_index_andrews` | `andrews.aperture.averaged_index` |
| `gaussian_fried.plane_wave_coherence_radius` | `andrews.structure.coherence_radius` |
| `gaussian_fried.plane_wave_fried_parameter` | `andrews.structure.fried_parameter` |
| `gaussian_fried.output_beam_params`, `effective_beam_params`, `rytov_std` | `andrews.beam.beam_params`, `andrews.beam.effective_beam_params`, `andrews.scintillation.rytov_variance` |
| `ao.plane_wave_fried_parameter_profile` | `andrews.structure.fried_parameter` |
| `angle_of_arrival.aperture_arrival_angle_variance` | `andrews.structure.angle_of_arrival_variance` |

The Churnside aperture-averaging trio in `plane_wave_scintillation.py` keeps its
own body, because those constants are Churnside 1991, DOI 10.1364/AO.30.001982,
not Andrews. Each docstring names its book-form alternative in
`andrews.aperture`.

#### The functions own their assumptions

Every public physics function in `olb/turbulence/**` now carries an
`@assumes(...)` decorator that attaches a machine-readable record and optional
runtime checks (see [architecture.md](architecture.md) Section 5 and
[api-budget.md](api-budget.md) for the mechanism). A Term factory opens
`trace_assumptions()` around its physics calls, so the Term inherits the union of
every function's assumptions automatically. The prose limit and the runtime check
are now one object: for example `scintillation.rytov_variance` carries the
`regime` constraint with a check on the Rytov variance it returns, and
`paths.*` carries the `zenith` check (its first enforcement anywhere).

One spectrum nuance is worth a headline. `andrews/scintillation.py` sets a MODULE
default `spectrum=SPECTRUM_KOLMOGOROV`, because the mainline indices are
Kolmogorov (no inner or outer scale). The opt-in two-scale `l0`/`L0` branches
(`weak_two_scale_index`, and the large-scale and small-scale log variances)
compute on the MODIFIED atmospheric spectrum, Ch. 9, Eqs. (48), (75) and (104).
So a Term that reads a two-scale branch inherits the module's Kolmogorov label
unless the factory states the modified spectrum. A factory that turns on `l0` or
`L0` must set `spectrum=SPECTRUM_MODIFIED` on the merge; the worked example in the
module self-check shows this.

#### Documented refusals

The package refuses a form that the book does not print. It does not guess a
coefficient. Each refusal raises `NotImplementedError` with its citation.

- **The annular (centrally obscured) receive aperture.** A full-text search of
  all 809 pages finds no aperture filter, no modulation transfer function and no
  flux variance for an obscured receive aperture. Secs. 10.3.1 to 10.3.6 use the
  soft Gaussian aperture or the unobscured circular transfer function of Ch. 10,
  Eq. (54), printed p. 410. So this gap needs another source. The Terms that
  assume an unobscured aperture keep saying so through `olb/assumptions.py`.
- **The Gaussian-beam two-scale strong chain**, in
  `scintillation.large_scale_log_variance` and in `aperture.averaged_index`. It
  needs eta_X of Ch. 9, Eq. (109), printed p. 355 (repeated as Ch. 10, Eq. (84),
  printed p. 420). No reading recovered from the scan gives both the plane-wave
  value 2.61 and the spherical-wave value 8.56 in the two limits.
- **The Gaussian row of the modified spectrum**, App. III, Table III, printed
  p. 766, in `structure.wave_structure_function`. As read, the two bump terms
  fall only as Lambda^(1/6), which breaks the plane-wave reduction by 2.3 %.
  Ch. 6, text below Eq. (77), printed p. 197, states the row MUST reduce.
- **The weak Gaussian-beam flux variance**, Ch. 10, Eq. (78), printed p. 419.
  The book prints no closed form for that double integral.
- **An inner scale or an outer scale on any Ch. 12 slant form**, in
  `paths.downlink_scintillation_index` and `paths.uplink_scintillation_index`.
  Chapter 12 uses the Kolmogorov spectrum only (Ch. 12, Eq. (15), printed
  p. 490). Use `scintillation.weak_two_scale_index` on a single homogeneous
  path.
- **An aperture-averaged downlink index in the strong regime.** Ch. 12,
  Eq. (39), printed p. 496, is a weak form and Ch. 12, Eq. (40) is a point form;
  the book prints no product of the two. So the gamma-gamma downlink Term models
  a POINT receiver, and its `Assumptions` record says so. That fade is deeper
  than the true aperture fade, which is the safe direction.
- **A temporal spectrum with a finite inner or outer scale, in any regime**, in
  `temporal.irradiance_temporal_spectrum`. Ch. 9, Sec. 9.8, printed p. 364,
  states that it ignores both scales.
- **A strong-regime spherical wave or Gaussian beam in the temporal module.**
  Ch. 9, Sec. 9.8, printed p. 364, limits that analysis to a plane wave.

Two more limits are inherent, not refused. First, the slant forms use a
plane-parallel atmosphere with no Earth-curvature correction, and the book
limits the weak-fluctuation slant results to zenith angles that do not exceed 45
to 60 deg (`paths.ZENITH_LIMIT_DEG`). Second, the two-scale large-scale log
variance does NOT reduce to the Kolmogorov branch as the inner scale goes to
zero: Ch. 9, Eq. (54), printed p. 339, states its substitution for the case
rho_0 << l0 only. Treat that branch as a moderate-to-strong model WITH a real
inner scale, not as a superset.

#### Source

- L. C. Andrews and R. L. Phillips, *Laser Beam Propagation through Random
  Media*, 2nd ed., SPIE Press (2005). DOI: 10.1117/3.626196. Every equation in
  the package cites its chapter, equation number and printed page.
- The equation-by-equation audit, the conflicts and the measured findings are in
  [andrews-crosscheck.md](andrews-crosscheck.md).

---

## 6. Fibre coupling

### 6a. Analytic mean coupling (fidelity 0)

File: `olb/models/coupling/_common.py` (shared SMF physics), used by
`olb/models/coupling/downlink.py` (`downlink_coupling_term`) and
`olb/models/coupling/terrestrial.py` (`terrestrial_smf_coupling_term`)

#### What the code models

A single-mode fibre couples only the field that matches the fibre mode. The
analytic model gives the mean coupling efficiency `eta` from the residual phase
variance that the compensation stack leaves. The Term is deterministic. It carries
no fade.

An aperture (bucket) detector is phase-insensitive, so it reuses the downlink
aperture-averaged scintillation Term unchanged (Section 5b).

#### Governing relations

Two limits of one overlap physics, selected by `SMF_SMALL_RESIDUAL_LIMIT = 1.0`
rad^2:

    small residual (sigma2_res < 1): extended Marechal
        eta = eta_max * exp(-sigma2_res)
    large residual (sigma2_res >= 1): Dikmelik-Davidson uncorrected coupling
        eta = eta_max * [ 1 + sigma2_res / NOLL_PISTON ]^(-6/5)

The mean loss is `-10*log10(eta)`. The `-6/5` exponent gives the large-aperture
asymptote `eta ~ (r0/D)^2`, where only about `(r0/D)^2` coherent cells couple into
the one fibre mode. The two branches cross near `sigma2_res = 1` rad^2, where
Marechal gives `eta_max/e = 0.37*eta_max` and the Dikmelik-Davidson branch gives
about `0.44*eta_max`. The effective D/r0 comes from inverting the piston-removed
Noll relation `sigma^2 = NOLL_PISTON*(D/r0)^(5/3)`.

For a terrestrial link the code evaluates the same forms at the horizontal
Gaussian-beam r0 (Section 5e). This is the effective-r0 weak-turbulence
approximation.

#### The received curvature and the defocus penalty (terrestrial)

The factor `eta_max` above is the FLAT-wavefront mode match. A terrestrial
received beam is not flat: it is a diverging Gaussian, with the phase-front
radius `R_rx` at the receive aperture (`olb.beam.phase_front_radius`; Andrews and
Phillips, 2nd ed. (2005), DOI 10.1117/3.626196, Ch. 4, Eqs. (7) and (8), printed
p. 87). A thin lens of focal length `f` images that diverging input BEYOND its
focal plane:

    dz_curv = f^2 / (R_rx - f)   > 0

S. A. Self, "Focusing of spherical Gaussian beams," Appl. Opt. 22, 658 (1983),
DOI 10.1364/AO.22.000658. So the TRUE focus of a terrestrial received beam sits
at `z = f + dz_curv`, beyond the nominal focal plane, and NOT at `z = f`.

The olb convention: the terrestrial coupling Terms ALWAYS charge that curvature
defocus, at the ACTUAL fibre plane. The detector sits at `z = f + defocus_m`, so
its distance from the true focus is

    dz_eff = defocus_m - dz_curv

and every spot-size and aberration quantity uses `dz_eff`. `optimal_focus` keeps
its meaning: it is a focal-LENGTH rule (`a = 1.12`), and it never moves the
detector. To model a tracked (aligned) coupler, read `dz_curv` from the public
helper `olb.models.coupling.curvature_focus_shift(scenario)` and set
`detector.defocus_m` to it.

A fibre `dz_eff` from the true focus sees a quadratic (defocus) phase across the
pupil. The mode-overlap integral stays closed form: the `a^2` of the Gaussian
weight becomes the complex `a^2 - i*c`, so

    eta(a, c) = 2 a^2 | (1 - exp(-(a^2 - i c))) / (a^2 - i c) |^2
    c = pi * dz_eff * (D/2)^2 / (lambda * f^2)     [rad]

`smf_eta_defocused(a, c)` in `olb/models/coupling/_common.py` holds this. It
reduces EXACTLY to `eta_max(a)` at `c = 0`, and it depends on `|c|` only.
Sources: Shaklan and Roddier, Appl. Opt. 27, 2334 (1988),
DOI 10.1364/AO.27.002334 (the `a` parameter and the flat-wavefront overlap);
Ruilier and Cassaing, JOSA A 18, 143 (2001), DOI 10.1364/JOSAA.18.000143
(single-mode coupling with an aberrated pupil).

Both terrestrial branches use it, the turbulent one and the `turbulence=False`
one, because the curvature is STATIC optics, not turbulence. For a space link
`R_rx` is enormous, so `dz_curv` is about zero: the downlink Terms are unchanged.

Example (1550 nm, `L = 5 km`, collimated `w0 = 0.02 m`, `D = 0.2 m`,
`w_m = 25 um`, `f = 4.524 m`): `R_rx = 5131.5 m`, `dz_curv = +3.99 mm`. A fibre at
the focal plane has `c = -3.95 rad`, so `eta = 0.215` (6.68 dB), of which 5.79 dB
is the curvature penalty on the 0.8145 flat-wavefront value. A fibre moved to the
true focus pays no curvature penalty. See
`validation/defocus/fidelity2_mmf_coupling_gap.md`.

#### How the code avoids double-counting

The geometric Term already carries the free-space spread and the aperture
power-in-bucket capture. The SMF `eta` is a multiplicative fibre-coupling
efficiency on the aperture-collected field. So the SMF Term adds only the coupling
loss `-10*log10(eta)`. In the budget the receive-coupling Term replaces the
standalone scintillation Term.

#### Inputs and outputs

- Inputs: the receive terminal (aperture, obscuration, wavelength, compensation),
  the Cn2 profile, the geometry.
- Output: a coupling Term. The mean-only SMF Term is deterministic.

#### Assumptions and limits

- Mean-only. The Term carries no fade (no sampler, no quantile). A fade margin read
  from it is wrong. This is the fidelity-0 lock: the Term flags itself, so the
  budget reports no fade margin. Use `downlink_budget(fidelity=1)` (FAST) or
  `fidelity=2` (wave optics) for the fade.
- The Dikmelik-Davidson coupling assumes a uniform circular aperture with no
  central obscuration. The code flags an obscured receive aperture.
- The static factor `eta_max(a) = 2*[(1 - exp(-a^2))/a]^2` (Section 6c) assumes a
  UNIFORMLY illuminated aperture and a FLAT wavefront. The CURVATURE half of that
  limit is now MODELLED for a terrestrial link: the Term charges the received
  curvature through the defocus-aberrated form `eta(a, c)` above, always, at the
  actual fibre plane. What stays is the ILLUMINATION half: a near-field
  terrestrial link inside the Rayleigh range tapers the received Gaussian across
  the aperture. That residual taper error runs SAFE, because a
  Gaussian-into-Gaussian overlap can pass the 0.8145 top-hat value, so the
  constant is then CONSERVATIVE. An illumination-aware `eta_max` stays open (see
  backlog 0-P11). See `olb/models/coupling/_common.py`.
- An SMF detector with NO resolvable coupling optics (no `focal_length_m` and no
  `optimal_focus`) has no `a` and no `c`. It keeps the plain `eta_max` field, and
  the Term flags that the curvature penalty is NOT modelled there.
- A terrestrial scenario whose transmit terminal carries no `Transmitter` gives
  no received curvature, so the Term charges none and flags itself OPTIMISTIC.
- The thin-lens focus shift is a SMALL-shift geometry. The Term flags
  `R_rx <= 2*f`, where the image runs away, and it still computes the value.
- The code flags an effective D/r0 above `SMF_DEEP_TURBULENCE_DR0 = 10`, where the
  practical coupling curve is extrapolated.
- The terrestrial form adds the effective-r0 weak-turbulence caveat: it evaluates
  plane-wave, Kolmogorov, phase-only forms at the Gaussian-beam r0. It ignores
  beam-wave amplitude scintillation, beam wander, and near-field curvature.
- Measured validity: see Section 9d. With the received curvature charged, the
  terrestrial multimode Term reads about 1 to 1.5 dB more loss than the field, and
  that residual is the Airy-versus-Gaussian spot shape.

#### Source

- Y. Dikmelik and F. M. Davidson, "Fiber-coupling efficiency for free-space optical
  communication through atmospheric turbulence," Appl. Opt. 44(23), 4946-4952
  (2005), for the uncorrected coupling curve and its limits.
- Extended Marechal approximation (Chan and others), for the small-residual
  limit. Derivation and validity: T. S. Ross, Appl. Opt. 48(10), 1812 (2009),
  DOI 10.1364/AO.48.001812.
- Noll 1976, for the residual variance (Section 5f).
- S. A. Self, "Focusing of spherical Gaussian beams," Appl. Opt. 22, 658 (1983),
  DOI 10.1364/AO.22.000658, for the focus shift of a curved input.
- C. Ruilier and F. Cassaing, JOSA A 18, 143 (2001), DOI 10.1364/JOSAA.18.000143,
  for single-mode coupling through an aberrated pupil.
- Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196, Ch. 4, Eqs. (7) and
  (8), printed p. 87, for the received phase-front radius.

### 6b. FAST statistical coupling (fidelity 1)

File: `olb/models/fast.py`

#### What the code models

The code drives FAST (the `fast-aosim` package) to get the true LP01 modal overlap
under turbulence, with the fade. It forms the coherent overlap of the turbulent
aperture field with the back-projected fibre mode. This is a true modal-coupling
metric, not a Strehl proxy.

#### Governing relation

    eta(t) = |INT (Aperture * M_fibre) * exp(chi + i*phi) dA|^2
             / |INT (Aperture * M_fibre) dA|^2

FAST propagates Monte-Carlo phase screens (phase) with an aperture-averaged
log-normal scintillation (log-amplitude), and forms the mode overlap directly. The
Term maps FAST to a loss per sample:

    floor_db = -link_budget['smf_coupling']    static mode-match loss
    loss(t)  = floor_db - result.dB_rel(t)     dB_rel is the turbulence penalty

The static floor lives in FAST's diffraction limit, so the two add. The Term stores
the per-sample loss, so it gives an empirical mean, an empirical quantile, and a
resampling sampler.

The code maps the olb compensation stack to FAST AO parameters: an AO(n) stage maps
to `AO_MODE="AO"`, `MODAL=True`, `ZMAX=n`; a tip-tilt stack maps to `AO_MODE="TT"`;
an empty stack maps to `AO_MODE="NOAO"`. Subharmonics capture the low-order tilt.

#### Inputs and outputs

- Inputs: the receive terminal, the site (Cn2, wind), the orbit, the sample count
  (NITER), and optional extra FAST parameters.
- Output: a coupling Term with an empirical mean, quantile, and sampler.

#### Assumptions and limits

- Scalar elevation only in this first cut. An elevation array needs one FAST run
  per elevation.
- Point-ahead is off (`DTHETA=0`): the up-leg and down-leg anisoplanatism of a
  moving satellite is not modelled. The Term flags this only when the scenario
  carries a `precompensation` source, because the point-ahead decorrelation
  applies to a pre-compensated beam. A plain downlink receive coupling never
  uses it, so it gets no flag.
- The default scales (L0 = inf, l0 = 1 um) are the Kolmogorov limit. A finite L0 or
  a large l0 from `fast_params` makes it a von Karman spectrum, and the Term reads
  the label from the resolved scales.
- FAST models the phase with real Monte-Carlo screens (the phase-driven coupling
  fade is fidelity-1). It models the log-amplitude as an aperture-averaged
  log-normal, which holds only in the weak fluctuation regime. The code flags when
  the plane-wave amplitude `sigma2_I` exceeds the lognormal-PDF house rule
  `LOGNORMAL_PDF_LIMIT = 0.25` (the retired name `WEAK_FLUCTUATION_LIMIT` is gone;
  the one canonical 0.25 lives in `andrews/scintillation.py`). A deep
  coupled-power fade does not trip that flag, because that fade is phase-driven and
  modelled correctly.
- FAST is an optional dependency (GPLv3). The module imports it lazily.

#### Source

FAST (`fast-aosim`), the modal-overlap engine, GPLv3. The regime constant comes
from the plane-wave scintillation index (Section 5b; Andrews and Phillips, 2nd
ed. (2005), DOI 10.1117/3.626196, Ch. 12, Eq. (38), printed p. 495).

### 6c. Fibre tip-tilt walk-off (terrestrial SMF and MMF)

File: `olb/models/coupling/terrestrial.py`, `olb/turbulence/angle_of_arrival.py`

#### What the code models

A focusing optic of focal length `f` puts the collected beam onto the fibre tip.
A received tip-tilt of angle theta moves the focal spot by `f*theta`. When the
spot moves, less light couples into the fibre. This is the walk-off. The Term
carries a real fade. This is the terrestrial receive-side effect that the
mean-only coupling Term (Section 6a) cannot carry.

#### Governing relations

The static mode match sets the coupling `eta_max`, which the coupling Term of
Section 6a then corrects for the received curvature. For a
single-mode fibre with a mode field radius `w_m`, a spot radius
`w_s = lambda*f/(pi*(D/2))`, and a coupling parameter `a = w_m/w_s`:

    eta_max(a) = 2 * [ (1 - exp(-a^2)) / a ]^2

It peaks at `eta_max = 0.8145` near `a = 1.12`. This is the mode overlap of a
uniform circular aperture with a Gaussian fibre mode.

The received tip-tilt has two contributions:

- A. The beam-wander arrival tilt. The turbulence moves the beam centroid at the
  receiver by an offset `r_c`, so the beam arrives from an apparent direction
  `r_c/L`. The radial (2-axis) variance is `sigma2_theta = <r_c^2>/L^2`. The code
  reuses the beam-wander kernel, which integrates the beam width profile along the
  path (`wander_arrival_angle_variance`). A receive tip-tilt or AO stage tracks
  this tilt out (all-or-nothing, no bandwidth model).
- B. The receive mechanical jitter. A per-axis 1-sigma jitter `sigma` adds
  `2*sigma^2` to the radial variance.

The captured power under a lateral offset is the overlap of two Gaussians, the
spot at the detector (`w_det`) and the fibre mode (`w_m`):

    eta(dx)/eta_max = exp( -2 * dx^2 / w_eff^2 ),   w_eff^2 = w_det^2 + w_m^2

The two axes of the offset are i.i.d. Gaussian, so the loss in dB is exponential
with mean `(20/ln10) * 2 * sigma_d^2 / w_eff^2`, with `sigma_d` the per-axis
spot offset below. The loss grows without a limit as the tilt grows. At the true
focus `w_det = w_s` and `sigma_d = f*sqrt(sigma2_theta/2)`, which is the
focal-plane case.

A single-mode-fibre subtlety: at a fixed `a` the focal length cancels in this
mean, because `w_s` scales with `f`. So `f` sets `eta_max` through `a`, but it
does not change the angular sensitivity on its own.

#### The detector plane: defocus and the chief-ray levers

The detector does not have to sit at the focal plane. It sits at
`z = f + defocus_m` (`SMF.defocus_m`, `MMF.defocus_m`; `0.0` is the focal plane).
The received beam is a diverging Gaussian, so the TRUE focus is at
`z = f + dz_curv` with `dz_curv = f^2/(R_rx - f) > 0` (Section 6a; S. A. Self,
Appl. Opt. 22, 658 (1983), DOI 10.1364/AO.22.000658). So the detector is

    dz_eff = defocus_m - dz_curv

from the true focus. Two effects follow, and the Terms separate them:

- SPOT GROWTH (axial). The spot grows away from the TRUE focus:
  `w_det = gaussz(w_s, dz_eff)`, the Gaussian beam radius against the distance
  from a waist (Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196,
  Ch. 4). At large `|dz_eff|` it tends to the geometric blur
  `(D/2)*|dz_eff|/f`. The coupling scale uses `w_det`, so a spot away from the
  true focus spills more from the MMF core and matches the SMF mode less well.
- SPOT DISPLACEMENT (lateral). The spot centre moves by the ray-optics chief-ray
  geometry of a thin lens, with the PHYSICAL `defocus_m` (the detector position,
  not the focus position, sets the lever):

      d_spot = (f + dz)*theta,    dz = defocus_m

  `theta` is the received arrival tilt (beam wander plus receive mechanical
  jitter), a radial 2-axis Gaussian, so the per-axis variance is:

      sigma_d^2 = (f + dz)^2 * sigma2_theta/2

  At the focal plane (`dz = 0`) the lever is `f`; off focus the longer lever arm
  `(f + dz)` moves the spot more.

No double-count with the mean coupling Term: the tip-tilt appears once. When the
walk-off Term is active, the mean coupling Term (Section 6a) keeps the
higher-order residual only. A virtual tip-tilt removes the Noll tip-tilt from its
residual (`drop_tiptilt=True`). The walk-off Term owns the tip-tilt.

A multimode fibre is a light bucket: the core is a HARD disk of radius `a_core`.
It collects ALL the spot power inside the core, so the coupling is the encircled
energy of the displaced Gaussian spot, NOT a mode overlap. With the spot of
radius `w_det` at the offset `dx`:

    eta(dx) = 1 - Q1( 2*dx/w_det ,  2*a_core/w_det )   (Marcum Q-function)

At `dx = 0` this reduces to `eta_static = 1 - exp(-2*a_core^2/w_det^2)`. A small
spot deep inside the core loses nothing until it nears the edge (a flat-top
acceptance); at the core edge it collects about half the power (about 3 dB). The
Term averages the loss over the Rayleigh offset. This differs from a single-mode
fibre, whose acceptance is a Gaussian mode, not a hard disk.

`optimal_focus` derives the focal length. For a single-mode fibre it picks `f` so
`a = 1.12` (the eta_max peak). For a multimode fibre it matches the spot to the
core, `a_core/w_s = 1.12`, which gives about 92% static capture. It is a
focal-LENGTH rule ONLY: it never moves the detector, so it does not cancel the
curvature defocus. A fibre at the focal plane of a 5 km horizontal link keeps the
full `dz_curv` offset, and reads about 8.5 dB static MMF loss for the example
optics of Section 6a; the same fibre at the true focus reads 0.37 dB.

A multimode fibre also has an ANGULAR acceptance. The numerical aperture
`NA = n*sin(theta_a)` sets the largest ray angle the fibre guides. The focusing
optic makes a cone of half-angle `NA_optic = (D/2)/f`. A ray from aperture radius
`rho` focuses at angle `rho/f`, so only rays from `rho <= f*NA` stay within the
acceptance cone. For a uniform aperture the guided POWER fraction is

    eta_NA = min( 1 ,  (NA / NA_optic)^2 )         NA_optic = (D/2)/f

a flat multiplicative loss on the coupled power. The spot size and the cone angle
are locked by the diffraction invariant `w_s * NA_optic = lambda/pi`, so a shorter
focal length shrinks the spot (better spatial capture) but steepens the cone: once
`NA_optic > NA` the extra light is not guided. This is the etendue penalty a
core-radius-only bucket misses; the fibre mode capacity is `V^2/2` with
`V = (2*pi/lambda)*a_core*NA`. The gate is OFF when `MMF.numerical_aperture` is
None. It is a flat power-transmission factor only: it does not re-broaden the
focal spot (that would need a re-truncated aperture).

#### Inputs and outputs

- Inputs: the receive terminal (aperture, wavelength, fibre optics,
  `defocus_m`, `pointing_jitter_rad`, compensation stack),
  the transmit waist and divergence, the path length, and the constant Cn2.
- Output: an SMF walk-off Term (category `pointing`), or an MMF coupling Term
  (category `coupling`). Both carry a real fade.

#### Assumptions and limits

- The walk-off falloff uses a Gaussian fit to the spot, and the eta_max value
  uses the more exact Airy-to-Gaussian overlap. The two spot models differ, which
  is standard practice near the peak. Against the fidelity-2 field the Gaussian
  spot model is the residual light-bucket gap of about 1 dB (backlog 2-W1).
- The MEAN modal penalty of a defocused single-mode fibre IS modelled, in the
  coupling Term of Section 6a (the closed form `eta(a, c)`). What stays geometric
  is the walk-off DISPLACEMENT response of `terrestrial_smf_walkoff_term`: it
  overlaps the defocused spot with the fibre mode, so it does not model how the
  defocus phase reshapes the modal overlap against a displacement. That walk-off
  fade is therefore OPTIMISTIC off focus, and the Term flags it loudly whenever
  `defocus_m` is not zero. Use an MMF (a light bucket), fidelity 2, or the full
  modal treatment of Ruilier and Cassaing (DOI 10.1364/JOSAA.18.000143).
- The defocus model is geometric: the spot keeps its Gaussian shape and only
  grows and moves. The chief-ray levers are ray optics.
- Contribution B of the received tip-tilt is the beam-wander tilt (A) only. The
  aperture angle-of-arrival "corrugation" tilt is available but feeds no Term:
  `aperture_arrival_angle_variance` in `olb/turbulence/angle_of_arrival.py` now
  delegates to `andrews.structure.angle_of_arrival_variance` (the gradient-tilt
  form, C-04), but no coupling Term adds contribution B. So the received tip-tilt
  is a lower bound. See `docs/andrews-crosscheck.md` batch 2 and backlog 0-W3.
- The beam-wander tilt is a weak-fluctuation model, so `sigma2_theta` is valid in
  weak turbulence only. The walk-off mapping itself has no upper limit. The SMF
  walk-off Term declared the weak regime and never flagged; that gap is CLOSED
  (2026-09-04) by a FUNCTION-OWNED check. The vendored Dios wander kernel
  `coupled_flux.beam_wander_variance` now takes an optional `wavelength` and runs
  the shared beam-aware `rytov_weak` gate itself, `wander_arrival_angle_variance`
  passes the keyword on, and the Term inherits the violation through the trace.
  The kernel value does not change: with no wavelength the check does not run.
  The MMF coupling Term reads the same wander model and asks for the check
  too (owner decision, 2026-09-04): a strong path now flags it, where before
  it read ok.
- The MMF `optimal_focus` is a geometric spot-to-core match, not a mode-overlap
  optimum: a shorter focal length captures more, up to the numerical-aperture gate.
- The MMF numerical-aperture gate is a flat power-transmission factor. It does not
  re-broaden the focal spot, and the spot itself stays diffraction-limited (no
  turbulence blur), so the full mode-count saturation is not modelled.

#### Source

- Shaklan and Roddier, "Coupling starlight into single-mode fiber optics," Appl.
  Opt. 27(11), 2334-2338 (1988), DOI 10.1364/AO.27.002334, for `eta_max(a)`.
- Dios et al., Applied Optics 43 (2004) 3866, DOI 10.1364/AO.43.003866, for the
  beam-wander offset that gives the arrival tilt.
- Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196, for the Gaussian
  power falloff and the 2-D Gaussian jitter.
- J. I. Marcum, "A statistical theory of target detection by pulsed radar," RAND
  RM-753 (1950), for the encircled energy of an off-axis Gaussian (the Marcum Q
  function) that gives the multimode-fibre coupling.
- Snyder and Love, Optical Waveguide Theory (1983), DOI 10.1007/978-1-4613-2813-1,
  for the numerical aperture, the acceptance cone, and the V-number that set the
  multimode-fibre angular gate.
- S. A. Self, Appl. Opt. 22, 658 (1983), DOI 10.1364/AO.22.000658, for the focus
  shift `dz_curv` of the received curved wavefront.
- C. Ruilier and F. Cassaing, JOSA A 18, 143 (2001), DOI 10.1364/JOSAA.18.000143,
  for the aberrated single-mode coupling that the walk-off response omits.

---

## 7. The fidelity-2 turbulent split step

Files: `olb/waveoptics/turbulence/` (`screens.py`, `splitstep.py`, `sampling.py`,
`run.py`, `temporal.py`)

### What the code models

The code moves a real complex field along the path and it puts a random phase
screen at each slab of that path. Each seed gives one SNAPSHOT of the atmosphere.
There is no time axis (`temporal.py` is PLANNED, NOT BUILT). The layer builds NO
Term and it changes NO budget. It is the reference against which the analytic
Terms of Sections 5 and 6 are read. For the API, see
[api-waveoptics.md](api-waveoptics.md) Section 9.

### Governing relations

The split step is propagate, apply a thin phase screen, propagate again. See
Schmidt, *Numerical Simulation of Optical Wave Propagation with Examples in
MATLAB*, DOI 10.1117/3.866274, Ch. 9. Each screen is a thin pure phase element,
`E_out = E_in * exp(i*phi)`, and each hop uses the `Forvard` FFT angular
spectrum. Each screen carries the integrated `Cn2` of its slab through the Fried
parameter `r0 = (0.423 k^2 INT Cn2 dz)^(-3/5)` (Fried,
DOI 10.1364/JOSA.56.001372; Andrews and Phillips, DOI 10.1117/3.626196, Ch. 12,
Eq. (23)), and the screens add as `r0 = (SUM r0_i^(-5/3))^(-3/5)`. The screen
spectrum is the modified von Karman form of Schmidt Ch. 9.

Five sampling rules size the grid and the screen count. `sampling.py` states each
one with its source, it WARNS when the grid misses one, and the `SamplingReport`
gives the ACHIEVED number:

1. **The extent.** `side = [guard*2*r_beam + 2*(lambda/r0_total)*z] / (1 - b)`.
   The first part is the vacuum extent rule. The second part is the scattering
   cone: turbulence scatters light through the angle `lambda/r0`, and that light
   must stay off the edge of the periodic grid. The divisor makes room for the
   absorbing band of width `b`. Schmidt, DOI 10.1117/3.866274, Ch. 9.
2. **The coherence pixel.** `dx <= r0_total / pixels_per_r0`. Schmidt,
   DOI 10.1117/3.866274, Ch. 9, and Martin and Flatte,
   DOI 10.1364/AO.27.002111.
3. **The Fresnel pixel.** `dx <= sqrt(lambda z_i) / 2` for the distance from
   screen `i` to the receiver. It keeps the irradiance correlation width sampled.
   Andrews and Phillips, DOI 10.1117/3.626196, Ch. 8. A screen that carries less
   than `fresnel_weight_min` of the total Rytov variance is exempt.
4. **The per-screen strength.** The plane-wave Rytov contribution of one screen,
   `2.25 k^(7/6) (INT Cn2 dz) (z_to_rx)^(5/6)`, must stay under
   `sigma2_r_screen_max`. A stronger screen breaks the thin-screen approximation.
   Andrews and Phillips, DOI 10.1117/3.626196, Ch. 8, Eq. (20), and Ch. 12,
   Eqs. (36) and (38); Schmidt, DOI 10.1117/3.866274, Ch. 9.

5. **The screen-count floor.** A weak path passes rule 4 with one screen, but
   one screen gives phase only and no scintillation. So `_merge_layers` clamps
   the count UP to exactly `min_screens` contiguous `Cn2`-weighted groups. The
   count follows the PRESET, not the layer count of the `Cn2` profile: a
   20-layer profile and a 200-layer profile of the same atmosphere give the
   same screens. A profile that has fewer layers than `min_screens` warns and
   keeps its layers, because the planner does not split a layer. THE FLOOR IS
   olb EVIDENCE, NOT BOOK PHYSICS. Schmidt gives no screen-count floor. An olb
   convergence sweep sets 15 / 9 / 5: it holds the grid fixed and it moves the
   count only, and the aperture scintillation index of a 30 degree downlink
   slab is 19 percent low at 3 screens, 10 percent low at 5, and flat from 7
   up, while the mean collected power holds inside 0.11 dB everywhere. The
   absolute lower bound is 4, the moment count of Schmidt,
   DOI 10.1117/3.866274, Ch. 9, Eq. (9.65), printed p. 164. The grouping does
   not SOLVE that equation, but the `Cn2`-weighted centroid matches all 8
   moments of the default profile to better than 1 percent.

A sixth, weaker limit gives `PIXELS_PER_FEATURE = 8` pixels across the smallest
hard edge.

### The reciprocity route (uplink)

The satellite of a space link sits outside the atmosphere, so a space scenario
ALWAYS simulates the DOWNLINK slab: a unit plane wave enters at the top of the
atmosphere on a flat grid. The uplink direction is never propagated. The
turbulent atmosphere is reciprocal, so the uplink flux at the satellite is the
overlap of the received downlink field with the ground transmit mode. See
Shapiro, "Reciprocity of the turbulent atmosphere," DOI 10.1364/JOSA.61.000492.
The code reads

    eta_turb = |SUM E_rx conj(psi_tx)|^2 / |SUM E_vac conj(psi_tx)|^2

with `psi_tx` the normalised ground transmit mode and `E_vac` a zero-screen
vacuum run through the SAME mask and the SAME hops. So the vacuum limit is
exactly 1.0, and `-10*log10(eta_turb)` sits on the free-space baseline of the
analytic Terms. Point-ahead anisoplanatism is NOT modelled here: the uplink and
the downlink read the same screens.

### Two numerical gotchas

- **The subharmonics fight the periodic propagator.** The subharmonic content of
  a screen is not periodic on the grid, and `Forvard` is periodic. So a run with
  subharmonics needs the absorbing boundary mask, and
  `propagate_turbulent_scenario` always applies one. The sub-steps alone remove
  no aliasing: the sampled transfer function of one long step is the product of
  the sampled transfer functions of the sub-steps, so a split hop gives the same
  array as one long hop. Only the mask helps, because it removes the edge energy
  between two sub-steps. Without the mask the wrap RAISES the measured variance;
  the self-check of `olb/waveoptics/turbulence/splitstep.py` therefore drops the
  subharmonics for its unmasked plane-wave Rytov case. Schmidt,
  DOI 10.1117/3.866274, Ch. 9.
- **The structure-function read-back is biased.** A Fourier screen holds no power
  above the grid Nyquist frequency and too little power below `1/(n*dx)`, so the
  measured `D_phi(r)` stays BELOW the theory `6.88 (r/r0)^(5/3)`. The three
  subharmonic levels lift it but they do not close the gap: the self-check of
  `olb/waveoptics/turbulence/screens.py` measures a deficit inside 15 percent
  over `r/r0` from 0.3 to 1.6. An `r0` read back from `D_phi` carries the same
  bias, so read it as a RATIO between two cases, not as an absolute value. Fried,
  DOI 10.1364/JOSA.56.001372; Schmidt, DOI 10.1117/3.866274, Ch. 9.

A related rule sits in `screens.py`: make the screen AT the propagation pitch.
A coarse screen that an FFT interpolates up carries no power above its own
Nyquist frequency, so it loses the structure at the Fresnel scale
`sqrt(lambda*z)` that builds the scintillation. That route is the documented
anti-pattern.

### Inputs and outputs

- Inputs: a `SpaceScenario` or a `TerrestrialScenario`, a geometry with ONE
  range, a quality preset, an optional `Cn2` profile and height grid, an outer
  scale, and a seed.
- Output: a `TurbWaveResult` of independent snapshots. Each `TurbTrial` holds the
  collected power, the single-mode-fibre coupling efficiency, the uplink
  `eta_turb`, its seed key, and its wall time. No Term, and no budget change.

### Assumptions and limits

- SNAPSHOT only. There is no fade rate and no fade duration.
- The `"retro"` direction and `folded_terrestrial()` raise
  `NotImplementedError`. The two passes of a retroreflector share the same
  screens, so they are correlated, and that correlation needs its own design.
- The split step runs on a FLAT grid only. `Screen()` and `split_step()` raise on
  a spherical (co-moving) field.
- The DEFAULT random draw comes from the self-contained `ScreenFactory`
  (`screen_generator="olb"`), which imports numpy and scipy only. It derives
  from Schmidt, DOI 10.1117/3.866274, Ch. 9, and the subharmonics from Lane,
  Glindemann and Dainty, DOI 10.1088/0959-7174/2/3/003. The opt-in `aotools`
  generator (`screen_generator="aotools"`, LGPL-3.0, the optional extra
  `screens`) is the reference path; `olb` imports it lazily and does not copy it.
  The two give different draws for the same seed; the statistics agree.
- Measured validity: see Sections 9g (how much tilt a screen holds), 9h (the two
  generators against each other, the analytic index and the fade tail) and 9i
  (the screen-count floor).

### Source

- J. D. Schmidt, *Numerical Simulation of Optical Wave Propagation with Examples
  in MATLAB*, SPIE Press (2010), DOI 10.1117/3.866274, Ch. 6, Ch. 7 and Ch. 9:
  the split step, the Fourier screen, the subharmonics, the absorbing boundary,
  the grid rules, and the range limit `z_max = N dx^2 / lambda`.
- J. H. Shapiro, "Reciprocity of the turbulent atmosphere,"
  DOI 10.1364/JOSA.61.000492: the uplink overlap.
- D. L. Fried, DOI 10.1364/JOSA.56.001372: the Fried parameter and the phase
  structure function.
- J. M. Martin and S. M. Flatte, DOI 10.1364/AO.27.002111: the
  pixel-per-coherence-length rule. The convergence practice is
  DOI 10.1364/JOSAA.7.000838.
- Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196, Ch. 7, Eq. (57),
  Ch. 8, Eq. (20), and Ch. 12, Eqs. (14), (23), (36) and (38): the long-term beam
  radius, the Rytov variance, the slant secant and the path weights.

---

## 8. The Schmidt foundation layer

Package: `olb/waveoptics/schmidt/`. Examples: `examples/schmidt/`.

### What the code models

The four modules of `olb/waveoptics/schmidt/` hold the NUMERICAL method of one
book:

- J. D. Schmidt, *Numerical Simulation of Optical Wave Propagation with Examples
  in MATLAB*, SPIE Press Monograph PM199 (2010). DOI: 10.1117/3.866274.

Section 5h holds the twin of this layer for the ANALYTIC physics. The two books
divide the work. Andrews and Phillips owns the analytic value of a quantity.
Schmidt owns the simulation rule: the transforms, the propagators, the sampling
constraints, the absorbing boundary, and the phase screens. A conflict between
the two is a real finding, and `docs/schmidt-crosscheck.md` records it.

Each function names its chapter, its equation number and its printed page. The
citation format is `Schmidt (2010), DOI 10.1117/3.866274, Ch. N, Eq. (nn),
printed p. NNN`. The package holds physics only. It imports numpy and scipy
only, it imports nothing from the rest of `olb`, and it returns no decibels.

| Module | What it gives | Book chapters |
| --- | --- | --- |
| `fourier.py` | The scaled two-dimensional transform pair `ft2` and `ift2`, the frequency pitch `freq_pitch`, and the transform-domain structure function `structure_function` | Ch. 2, Eqs. (2.3), (2.6), (2.9), (2.32); Ch. 3, Eqs. (3.15) to (3.25) |
| `fresnel.py` | The four propagation kernels `one_step_fresnel`, `two_step_fresnel`, `angular_spectrum` (baseline and scaled) and `partial_propagations`, plus the book absorber `super_gaussian_absorber` | Ch. 6, Eqs. (6.5), (6.15), (6.16), (6.18) to (6.25), (6.31), (6.32), (6.65); Ch. 8, Eqs. (8.1), (8.8), (8.14) to (8.18) |
| `sampling.py` | The four numbered constraints, the local-frequency analysis they come from, the per-kernel bounds, the partial-propagation planner, and the `check_sampling` rule table | Ch. 7, Eqs. (7.7) to (7.60); Ch. 8, Eqs. (8.23), (8.24) |
| `turbulence.py` | The three phase PSDs, the Fourier and subharmonic screen generators, the per-screen strength rule and its `rmax` cap, the layer moment rule, the turbulent sampling bounds, and the `properly_sampled_checklist` | Ch. 9, Eqs. (9.44), (9.49) to (9.52), (9.63) to (9.65), (9.70) to (9.75), (9.78) to (9.81), (9.84) to (9.90) |

### The layer is VALIDATION only

No budget, no Term, no sizer and no runner consumes this layer. The production
wave-optics code stays the LightPipes port of Section 7, and it keeps its
bodies. The Schmidt layer measures that code from the outside, and it gives the
book number that a future revision can move to.

The retrofit of work package 6 wrote the book citations INTO the production
modules, docstrings and comments only. Two of those citations changed a claim:

- `olb/waveoptics/grid.py` `forvard_max_z` cited "Ch. 6". The rule is
  constraint 4, Ch. 7, Eq. (7.59), printed p. 127, at m = 1, and the same
  expression is the step cap of Ch. 8, Eq. (8.24), printed p. 144. The constant
  is DERIVED, not a guess.
- `olb/waveoptics/lenses.py` hinted at Ch. 7 for the co-moving grid. The book
  develops no such grid. See the refusals below.

### The three example scripts

Each script reads the production layer and the Schmidt layer, prints a labelled
table, and saves a figure. No script changes an `olb` module.

- `propagator_kernels.py` — the book kernels against the production
  propagators, in three tiers. Tier (a), `angular_spectrum` at m = 1 against
  `Forvard`, is one algorithm on one grid, and the two agree to about 1e-10.
  Tier (b), the one-step and two-step Fresnel kernels against the production
  `Fresnel`, crosses a quadrature: the interior agreement is 6e-4 for a soft
  Gaussian and 1.5e-2 for a hard truncation. Tier (c), the two-step kernel
  against the production `Lens -> LensFresnel -> Convert` recipe at a
  magnification of 247, gives 1.7e-3 soft and 2.3e-2 hard.
- `sampling_and_edges.py` — a gallery of deliberate sampling failures, each
  paired with the grid that obeys the rule, and then the rule checker on the
  real production grids. It also plots the two absorber shapes on one axes.
- `screens_and_turbulence.py` — the screen generators against Eq. (9.44). The
  book subharmonic generator reaches 0.88 to 0.93 of theory over
  `r/r0 = 0.3` to 1.6, and the `aotools` generator of the production layer
  reads 1 to 3 percent above it, so the two agree well inside the band. Both
  fall away past it, and no subharmonic level removes that deficit. The script
  also proves the factor-4 bridge from the live code: the production per-screen
  number is the plane-wave Rytov variance and the book cap `rmax = 0.1` is on
  the log-amplitude variance, and the measured ratio is 3.9994.

### Documented refusals and absences

The book gives less than a reader expects in four places. Each absence is
recorded, not filled with a guess.

- **No obscured (annular) aperture.** D1 and D2 are plain extents everywhere in
  Chs. 7 to 9. This is the same gap as the Andrews layer, Section 5h.
- **No co-moving (spherical) grid.** Ch. 6, text, printed p. 87, names the
  Coles and Rubio angular-grid method and then does not develop it. The book
  never leaves the flat grid; its own answer to the same problem is the scaling
  parameter m of Ch. 6, Eq. (6.65), printed p. 100. So the production
  `LensFresnel` and `Convert` route has NO book equation to check against.
- **No temporal axis.** Sec. 9.5.4, printed p. 179, states the frozen-flow
  method in prose and points to the Greenwood frequency. It gives no equation
  and no code.
- **No screen-count floor.** Ch. 9 gives an UPPER bound on one screen's share
  (`rmax = 0.1`, Listing 9.5, printed p. 175) and a worked example with 11
  planes and no criterion (Sec. 9.5.2, printed p. 177). Eq. (9.90), printed
  p. 174, is a sampling floor of the FFT method, not of the atmosphere. So the
  production `QualityPreset.min_screens` (15 / 9 / 5) cannot be sourced to
  Schmidt, and neither can any other integer. olb states the floor as olb
  evidence instead: an internal convergence sweep sets it, and the absolute
  lower bound of 4 comes from the layer moment rule, Eq. (9.65), printed
  p. 164, which gives 8 equations against 2 free numbers for each screen. See
  the sampling rule 5 above, and WP7 in
  [schmidt-crosscheck.md](schmidt-crosscheck.md).

### Source

- J. D. Schmidt, *Numerical Simulation of Optical Wave Propagation with Examples
  in MATLAB*, SPIE Press (2010). DOI: 10.1117/3.866274. Every equation in the
  package cites its chapter, equation number and printed page.
- The equation-by-equation forward map, the 28 gaps, the constants ledger and
  the work-package notes are in
  [schmidt-crosscheck.md](schmidt-crosscheck.md).

---

## 9. Measured validity: what the validation scripts certify

Sections 1 to 8 give what each model COMPUTES. This section gives where each
model HOLDS. Every entry below is a MEASUREMENT, not a derivation.

The reference is the fidelity-2 field solve of Section 7. It solves the field on
a grid, and it makes no beam assumption, no regime assumption and no
distribution assumption. So it is the in-repo reference for the analytic
(fidelity-0) and statistical (fidelity-1) models, in the band where the grid
itself is trustworthy. Entry 9f gives one case where it is NOT trustworthy, and
where the analytic Term is the reference instead.

Each entry gives the physics QUESTION, the model under test, the reference, the
measured outcome with its date and its script, and a one-line VERDICT. The
scripts live in [validation/](../validation/), and
[validation/README.md](../validation/README.md) indexes them.

**The standing rule.** A validation script that reaches a verdict adds an entry
here, or it updates the entry that it supersedes. A result that lives only in a
run log, a memory note or a backlog aside is not documented.

### 9a. Where does the Dios uplink scintillation index apply?

- **Question.** Does the fidelity-1 coupled-flux uplink give the correct
  irradiance fluctuation of an uncorrected ground-to-space beam?
- **Model under test.** Section 5c, the coupled-flux Monte Carlo, and its
  off-axis term of Section 5d. Dios et al., Applied Optics 43 (2004) 3866,
  DOI 10.1364/AO.43.003866, Eq. (20).
- **Reference.** The fidelity-2 split step read through the Shapiro reciprocity
  overlap, DOI 10.1364/JOSA.61.000492. One zoom transform gives the uplink flux
  at every satellite offset, so each INGREDIENT of the analytic model (the
  on-axis index, the wander, the short-term waist, the beam-frame index) is
  measured on its own.
- **Measured (2026-08-28).** 600 km, 60 deg elevation, 1550 nm, HV57 site
  profile. For an UNDERFILLED launch (`w0 = 0.06 m`) fidelity 1 sits 1.2 to 1.3
  times above fidelity 2, inside the joint error bars. For a FILLED launch
  (`w0 = 0.18 m`, the wander comparable with the far-field spot) fidelity 1
  reads 2.3 times high at 0.3 times the profile strength, and 7.2 times high at
  the full profile. The cause is measured, not guessed: the analytic beam-frame
  index at the same point is 0.018, against the 0.85 of the fidelity-1 total
  that its unsaturated off-axis Rytov term contributes. The short-term waist is
  correct to a few percent. The fidelity-2 answer is CONVERGED: the reference
  preset, a doubled grid side and the standard plan agree inside the 15 percent
  Monte-Carlo noise. The measured wander variance is 1.8 to 2.0 times the
  Dios/Belmonte 2.07 form and 0.55 times the Andrews 7.25 form, at every case.
  The current fidelity-1 validity flag is a MEAN over the samples, so it can
  read "valid" (0.167 against the 0.25 limit) while the model reads 2.3 times
  high.
- **VERDICT.** Use fidelity 1 for an uncorrected uplink with an UNDERFILLED
  launch (`beta_rms / w_L` well below about 0.5). Do not use it for a filled
  launch: it is pessimistic there by 2 to 7 times on `sigma2_I`, and its
  weak-fluctuation flag does not catch the case. Use fidelity 2 instead.
- **Script.** `validation/uplink_sigma2i/uplink_farfield_reciprocity.py`;
  write-up
  [validation/uplink_sigma2i/UPLINK_SIGMA2I_INVESTIGATION.md](../validation/uplink_sigma2i/UPLINK_SIGMA2I_INVESTIGATION.md).

### 9b. Are the vendored coupled-flux kernels faithful to the paper?

- **Question.** Did the vendored copy of the Dios kernels
  (`olb/turbulence/coupled_flux.py`) keep the published behaviour?
- **Model under test.** The fidelity-1 uplink kernels, on the exact published
  case: a GEO uplink at 0.84 um, the log-amplitude variance against the
  transmit waist, at 90 and 30 deg elevation.
- **Reference.** Dios et al., DOI 10.1364/AO.43.003866, Fig. 5, and its own
  FFT-BPM points. The fidelity-2 leg runs the same case.
- **Measured (2026-08-28).** The fidelity-1 curve overlays the figure: the
  90 deg plateau reads 0.0298 against the printed 0.028, and the 30 deg plateau
  reads 0.0938 against the printed 0.095. The fidelity-2 points reproduce the
  paper's FFT-BPM behaviour at every station. They sit on the reference at a
  small waist (0.109 against about 0.105), they overshoot the weak-theory curve
  in the focusing range, and they SATURATE near 0.65 at `W0 = 0.10 m` while the
  fidelity-1 curve climbs to 0.96 (90 deg) and 2.46 (30 deg). The run also found
  a real DEFECT in the olb wrapper `olb.turbulence.uplink_flux._flux_result`: it
  puts the airmass on `Cn2` but it keeps the vertical height grid as the path
  coordinate, so the on-axis index scales as `sec(zeta)` and not as
  `sec(zeta)^(11/6)`. That is 40 percent low at 30 deg elevation and 13 percent
  low at 60 deg. The exact fix, validated against the figure, is to give the
  kernels the slant-mapped grid with the zenith profile.
- **VERDICT.** The vendored kernels are faithful. Trust the fidelity-1 uplink
  ONLY where the analytic curve has not left its weak band (see 9a). Read the
  on-axis index at a low elevation with the slant-coordinate defect in mind,
  until that fix lands.
- **Script.** `validation/uplink_sigma2i/dios_fig5_replication.py`.

### 9c. Does the uplink model see a central obscuration?

- **Question.** A launch telescope carries a secondary mirror. Does the
  fidelity-1 uplink read the obscuration, in the MEAN and in the FADE?
- **Model under test.** The transmit truncation Term of Section 2a
  (`tx_gaussian_efficiency_term`), and the coupled-flux index of Section 5c.
- **Reference.** The fidelity-2 far-field map of 9a, with an annular launch
  pupil. The obscuration ratio `eps` sweeps from 0 to 0.8.
- **Measured (2026-08-28).** 600 km, 60 deg, 1550 nm. The MEAN is NOT a blind
  spot: the analytic truncation Term tracks the wave-optics far-field mean
  inside 2.4 dB over a 60 dB sweep, and it is slightly conservative. The FADE
  is a blind spot. The coupled-flux kernel reads the launch through ONE number,
  the waist `w0`, so its `sigma2_I` is flat in `eps` by construction. The true
  index RISES: 1.5 times the unobscured value at an obscuration radius of
  0.44 waists (a filled launch), and 4.6 times at 2.33 waists (a small launch
  beam, where the obscuration blocks the core).
- **VERDICT.** Use the analytic MEAN loss at any obscuration. Trust the
  fidelity-1 FADE only when the obscuration radius stays well below one transmit
  waist. Past that, the fidelity-1 index is optimistic, and only fidelity 2
  gives the rise.
- **Script.**
  `validation/uplink_sigma2i/uplink_obscuration_dios_vs_waveoptics.py`.

### 9d. Does the terrestrial coupling Term agree with the field?

- **Question.** A terrestrial fidelity-2 run read about 7 dB more multimode
  coupling loss than the analytic Term of Section 6a. Which model is wrong?
- **Model under test.** `terrestrial_mmf_coupling_term` and
  `terrestrial_smf_coupling_term` of Section 6a, against the fidelity-2
  focal-plane coupling of `olb/waveoptics/mmf.py`.
- **Reference.** The fidelity-2 VACUUM field (no turbulence), plus an
  independent one-dimensional Fresnel (Hankel) quadrature of the truncated
  curved pupil that runs with no olb code in the loop.
- **Measured (2026-08-31).** 1550 nm, 5 km path, collimated `w0 = 0.02 m`,
  `D = 0.2 m`, core 25 um, `f = 4.5242 m`. NEITHER model was wrong. The received
  beam is a diverging Gaussian of `R_rx = 5131 m`, so the true focus sits
  `dz_curv = f^2/(R_rx - f) = +3.99 mm` BEYOND the focal plane (S. A. Self,
  Appl. Opt. 22, 658 (1983), DOI 10.1364/AO.22.000658). The analytic Term
  assumed best focus and the field observed at `f`. A defocus scan on the same
  vacuum field peaks at `+4 mm`, exactly the predicted shift, and the quadrature
  confirms it to better than 0.1 dB over a +/-8 mm sweep. The scan also found a
  SIGN fault: `olb.waveoptics.mmf.focal_intensity` applied the pupil phase with
  the wrong sign, so every earlier `defocus_m` scan came out mirrored. With the
  curvature charged in the analytic Terms, the gap falls from about 7 dB to
  about 1 to 1.5 dB (analytic 8.54 dB against the field 7.08 dB at the focal
  plane; about 1 dB at the true focus).
  The residual is the Airy-versus-Gaussian spot shape: the truncated pupil makes
  an Airy pattern whose slow rings a Gaussian spot model omits.
- **VERDICT.** Use the terrestrial coupling Terms for a finite-path link only
  with the received curvature charged, which they now do. Expect them to read
  about 1 to 1.5 dB MORE loss than the field, because of the spot shape. A space link
  is unaffected: `R_rx` is enormous there, so `dz_curv` is about zero. The
  fidelity-2 single-mode leg now takes the defocus too (backlog 2-W2, DONE
  2026-09-04), so the aberrated single-mode closed form HAS a field reference:
  the `olb/waveoptics/smf.py` self-check matches `smf_eta_defocused(a=1.12, c)`
  to four decimals at c = 0, 1, 2 and 4.
- **Script.** `validation/defocus/defocus_sensing.py`; write-up
  [validation/defocus/fidelity2_mmf_coupling_gap.md](../validation/defocus/fidelity2_mmf_coupling_gap.md).

### 9e. Is the aperture-averaged lognormal power draw trustworthy?

- **Question.** In weak turbulence, does the cheap analytic route
  `sigma2_P = A sigma2_I` plus a lognormal draw give a trustworthy received-power
  distribution? An aperture integrates a CORRELATED field, and a sum of
  lognormals is not a lognormal.
- **Model under test.** The three legs apart: the point index of Section 5d
  (the Dios on-axis Gaussian beam-wave form), the weak aperture-averaging factor
  `A` of Section 5b (the Churnside plane-wave fit,
  DOI 10.1364/AO.30.001982), and the lognormal SHAPE.
- **Reference.** The fidelity-2 split step. Fidelity 2 models no tip-tilt
  correction, so it holds the FULL beam wander.
- **Measured (2026-09-01, FULL run: 1500 trials for each launch on the
  `standard` preset).** One horizontal path (2 km, `Cn2 = 3e-15`, 1550 nm) that
  stays firmly weak (`sigma_R^2 = 0.21`), two launches, and four receive
  diameters, so `D/rho_0` runs from 0.20 to 7.89.
  - The lognormal FAMILY holds, to the deep tail. With the index refit to the
    MEASURED value, every case agrees inside 0.13 dB at the 5 percent fade, and
    inside 0.21 dB at the 1 percent fade. The skew of `ln P` sits near -0.2 in
    every case with no trend against the diameter. So the drift to a Gaussian
    power is not visible in this band, and the full beam wander of the sim
    builds no pointing tail that the lognormal cannot hold.
  - The FILTER is the fault, not the point index. The analytic point `sigma2_I`
    reads 3 to 4 percent HIGH only (the 10 to 20 percent bias of the quick run
    was a `rapid`-preset artifact). The Churnside factor OVER-AVERAGES by about
    1.2 times at `D/rho_0 = 1`, and by 1.8 (collimated) to 2.5 (diverged)
    times at `D/rho_0 = 3`.
  - One column is BEAM-FILLING-LIMITED, and it is not a filter error. The
    collimated 40 cm case catches `eta_fill = 0.87` of the beam, so it measures
    almost the total power, and `A_eff` falls BELOW the analytic `A` (ratio
    0.40). The diverged launch at the same diameter has `eta_fill = 0.39` and
    it behaves like every other unfilled case. The reversal is the fill
    fraction, not the diameter (backlog 2-N2).
  - The ABSOLUTE impact is small in this band. The fade spread falls from
    1.15 dB at `D = 1 cm` to 0.05 to 0.13 dB at `D = 40 cm`, so the worst
    relative index error (2.5 times) moves the 5 percent fade by 0.29 dB.
    The worst disagreement of the whole analytic route is 0.29 dB at the
    5 percent fade and 0.41 dB at the 1 percent fade (both the diverged
    15 cm case).
- **VERDICT.** Use the analytic lognormal draw for a weak horizontal path with
  a receive aperture that holds much less than half of the beam. The
  DISTRIBUTION SHAPE is certified there, to the 1 percent fade. Read the
  aperture-averaged INDEX as approximate over `1 <= D/rho_0 <= 8`: it is
  over-averaged, and only the small absolute fade spread keeps the budget error
  under about 0.3 dB at the 5 percent fade (0.41 dB at 1 percent). Do not read
  the analytic `A` at all when the aperture holds most of the beam. Still
  untested: a focused launch, a stronger `Cn2`, and a longer path (the 1-8
  gates).
- **Script.**
  `validation/lognormal_certification/lognormal_certification.py`; write-up
  [validation/lognormal_certification/README.md](../validation/lognormal_certification/README.md).
  See backlog 1-6.

### 9f. Which model gives the fidelity-2 geometric loss?

- **Question.** A fidelity-2 budget can take its no-turbulence geometric loss
  from the analytic Terms of Sections 1 and 2a, or from a wave-optics vacuum
  field solve. Which one is correct?
- **Model under test.** The wave-optics vacuum loss (`_full_vacuum_loss_db`)
  against the analytic pair `geometric_loss_term + tx_gaussian_efficiency_term`.
  The comparison is apples to apples: both sum the same two parts, and both use
  an aperture receiver.
- **Reference.** Each model is the reference of the other, on the path where it
  is trustworthy.
- **Measured (2026-08-31).** TERRESTRIAL, where the grid is small and well
  resolved: the two agree inside about 0.15 dB over 0.3 to 10 km collimated, and
  inside 0.002 dB on a tightly focused short path. This VALIDATES the analytic
  geometric Term against wave optics. SPACE: the wave vacuum loss of one 600 km
  case was measured as the grid refined from 4096 to 7168 pixels. It SCATTERS
  around the stable analytic 52.24 dB by +4.27, -1.10, +0.45 and -0.77 dB, a
  spread of 5.4 dB, and it does not converge. The cause is resolution: a
  practical grid cannot resolve the mm-scale aperture edges over a 2000 km path.
  That run also costs about 14 s.
- **VERDICT.** Use the ANALYTIC geometric Term for a space fidelity-2 budget. A
  ground-space link is always far field (the Fraunhofer distance of a 0.1 m
  aperture at 1550 nm is about 6 km), so the analytic form is exact, and it is
  both cheaper and more trustworthy than the field solve. `run_fidelity2` does
  this by default (`vacuum="analytic"`). Keep the WAVE vacuum for a terrestrial
  link, because the terrestrial turbulence penalty is turbulent / vacuum on the
  SAME grid, so the wave vacuum is the exact baseline that cancels the grid.
- **Script.** `validation/vacuum_loss/vacuum_loss_validation.py`.

### 9g. How much of the low-frequency (tilt) band does a phase screen hold?

- **Question.** A screen on a finite grid holds no power below its grid
  fundamental. That missing band is the tip and the tilt. How much does each
  screen route lose?
- **Model under test.** The screen routes of Section 7: the plain Fourier
  screen, the subharmonic screen (the production route), an oversized-and-cropped
  screen, and the extruded (infinite) screen of the planned temporal layer.
- **Reference.** The Noll per-axis Z-tilt filter integral,
  DOI 10.1364/JOSA.66.000207; the Andrews G-tilt filter,
  DOI 10.1117/3.626196, Ch. 6, Eq. (84), printed p. 201; the von Karman
  covariance of Assemat and Wilson, DOI 10.1364/OE.14.000988, Eq. (5); and the
  structure function `D(r) = 6.88 (r/r0)^(5/3)`, DOI 10.1117/3.866274, Ch. 9.
- **Measured (2026-08-28).** 1550 nm, `r0 = 0.10 m`, `D/r0 = 10`. With a PURE
  Kolmogorov spectrum no practical screen holds the tilt: the subharmonic route
  reaches 0.75 to 0.78 of the Noll value, and an oversize factor of 8 reaches
  0.84 with a wide error bar. With a finite outer scale (`L0 = 25 m`) both cures
  reach the analytic value: the subharmonic route reads 0.92 to 1.04, and the
  oversize route converges from a screen side of about `L0/2`. The infinite-screen
  covariance FORMULA is correct (it matches a float64 closed form to 3.9e-7). The
  screen variance is spin-up limited: 512 rows cover 2.0 outer scales at
  `L0 = 2.56 m` and land inside 0.08 percent, but they cover 0.2 outer scales at
  `L0 = 25 m` and read 12 percent low at every lag. The EXTRUSION AXIS is
  genuinely defective: the two-column Markov recursion over-correlates its own
  direction (1.498 rad^2 at a 2.56 m lag, against a theory of 0.085 rad^2), and
  the anisotropy of `D(r)` reaches 20 to 30 percent. The bias is stationary and
  it does not drift.
- **VERDICT.** The production SNAPSHOT route needs no change: the subharmonic
  screen matches the oversized screen on every measured metric. Select a finite
  outer scale from the physics, then size the screen side to one or two outer
  scales. UPDATE (2026-09-04, Section 9k): the subharmonic screen already
  reaches 27 x its side, so the rule is `L0 <= 27 x side` with three levels
  (raise the level count when a small grid breaks it), not a side of one or
  two outer scales; and a screen drawn with `L0 = inf` IS, in its measured
  statistics, a von Karman screen at `L0 = 27 x side`. Do NOT build the temporal layer on the extruded screen: under frozen
  flow (Taylor, DOI 10.1098/rspa.1938.0032) a row lag IS a time lag, so the
  extrusion defect smooths the temporal axis exactly, and a fade duration reads
  too benign. Prefer a shifted large screen.
- **Script.** `validation/screens/oversize_crop.py`,
  `infinite_screen_stats.py`, `extrusion_stationarity.py`; write-up
  [validation/screens/FINDINGS.md](../validation/screens/FINDINGS.md).

### 9h. Is the fast olb screen generator equal to the aotools reference?

- **Question.** The default screen generator is now the self-contained olb
  `ScreenFactory`. It draws a DIFFERENT random atmosphere from the `aotools`
  reference for the same seed. Do the two give the same statistics, and the same
  FADE TAIL?
- **Model under test.** `ScreenFactory` (`screen_generator="olb"`) against
  `aotools` (`screen_generator="aotools"`).
- **Reference.** The `aotools` generator, plus the analytic aperture-averaged
  index of Section 5b and the structure function
  `D_phi(r) = 6.88 (r/r0)^(5/3)`, DOI 10.1364/JOSA.56.001372.
- **Measured (2026-08-29).** Four cases (terrestrial 2 km, space downlink at
  30 deg on two presets, space uplink at 30 deg) agree on the mean collected
  power and on the aperture `sigma2_I` inside the combined Monte-Carlo bars; the
  largest gap is 0.76 sigma. The single-mode coupling and the uplink reciprocity
  `eta_turb` agree too (0.02 to 0.34 sigma). Converged over 2000 downlink trials
  the olb index is 0.01522, which is 1.012 times the analytic 0.01505, against
  0.994 for `aotools`. The tail number is the deciding one: the 1 percent fade
  quantile differs by 0.014 dB against a 0.067 dB bootstrap bar, and every
  quantile from 0.1 to 90 percent agrees inside 2 sigma. ONE caveat is common to
  both generators and is not a generator difference: the phase structure function
  sits 10 to 25 percent below the pure Kolmogorov value, which is the finite-grid
  low-frequency deficit of 9g. The two generators stay inside 5 percent of each
  other under `L0 = 25 m`.
- **VERDICT.** Use the olb generator. It is a trustworthy drop-in, and it is 7.6
  to 14 times faster per screen. Use `screen_generator="aotools"` to reproduce an
  older run bit for bit, because the seeded draws differ.
- **Script.** `validation/waveoptics_speed/generator_validation.py`.

### 9i. How many phase screens does a turbulent run need?

- **Question.** `QualityPreset.min_screens` is 15 / 9 / 5 for the reference,
  standard and rapid presets. Where do those integers come from?
- **Model under test.** The screen planner of Section 7,
  `olb/waveoptics/turbulence/sampling.py`.
- **Reference.** An olb convergence sweep. Schmidt gives NO screen-count floor:
  Eq. (9.90), printed p. 174, is a sampling floor of the FFT method, and the
  worked example of Sec. 9.5.2, printed p. 177, gives 11 planes with no
  criterion. The only book bound is the layer MOMENT rule, Ch. 9, Eq. (9.65),
  printed p. 164, which gives 8 equations against 2 free numbers for each screen,
  so 4 screens is the absolute lower bound.
- **Measured (2026-08-27).** The sweep holds the GRID fixed, it moves the screen
  count only, and it runs 200 snapshots for each count. The tolerance is 0.1 dB
  on the mean power and 5 percent on the index. The mean power meets it at every
  count of every case. The BINDING case is the 600 km downlink slab at 30 deg:
  the aperture index reads 19 percent low at 3 screens, 10 percent low at 5, and
  it stays flat from 7 up. A zenith slab and a 2 km horizontal path meet the
  tolerance at every count, because both are weak and homogeneous. Separately,
  the production grouping matches all 8 moments of the default profile inside
  0.15 percent, so the `Cn2`-weighted centroid satisfies Eq. (9.65) in practice
  without a solve.
- **VERDICT.** Use `standard` (9) or `reference` (15): both sit on the converged
  plateau. `rapid` (5) is a stated compromise, one step below the converged
  count of 7: its slant index runs about 10 percent low, and its mean power holds
  inside 0.11 dB. No preset may go under 4. The screen count follows the PRESET
  and not the layer count, so a 20-layer and a 200-layer profile of the same
  atmosphere give the same plan. Section 9k extends this sweep to the FADE
  TAIL at 1000 trials: the deep tail and every point quantity are flat from 5
  screens up, and the SMF p5 softens about 2 dB from 9 to 25 screens.
- **Script.** The sweep tables are in
  [schmidt-crosscheck.md](schmidt-crosscheck.md), the WP7 note.

### 9j. Do FAST and the Stone law agree on the point-ahead residual?

- **Question.** The pre-compensated uplink holds two models of the point-ahead
  decorrelation residual: the fidelity-1 FAST route (`uplink_fast_term`, the
  model of record; the PAOLA filter `G_AO_PAOLA`, Farley,
  DOI 10.1364/OE.458659) and the fidelity-0 Stone modal law
  (`uplink_point_ahead_term`; Stone, DOI 10.1364/JOSAA.11.000347). They compute
  the same quantity. Do they give the same number? Backlog 1-5.
- **Model under test.** Both routes at once, at matched conditions: the same
  HV5/7 profile, aperture (1.5 m), corrected order, and point-ahead angle, with
  the FAST servo and sensor effects OFF (`TLOOP=0`, `TEXP=0`, zero wind,
  `ALIAS=False`, `NOISE=0`), where the PAOLA filter reduces exactly to the pure
  two-path kernel `2 - 2 cos(delta_r . kappa)`.
- **Reference.** Each route checks the other; an independent whole-plane polar
  quadrature of the FAST filter breaks the tie, because the FAST grid truncates
  the low-frequency band. The mode sets are matched: the FAST modal mask keeps
  the piston and the tilts, so its analytic partner is the Stone band with NO
  mode removed, not the production `piston_tilt` form.
- **Measured (2026-09-02, FULL run: grid 1024 x 0.01 m, 3000 draws).**
  - The mode-matched ratio reads 1.044 to 1.055 across the whole sweep
    (point-ahead 0.25x to 2x nominal, ZMAX 1 to 66, elevation 30 to 90 deg);
    the single-layer case reads 0.991. The uncorrected anchor (ZMAX = 0) is
    exactly zero on both routes: with no corrected mode, nothing decorrelates. The fitting sides agree to 0.6 percent
    (FAST band integral 0.3334 rad^2 against Noll 0.3354 rad^2 at 55 modes).
  - The production pairing reads 3.5x at the production point, and the whole
    factor is the mode set: 2.08 rad^2 of piston plus 0.41 rad^2 of tilt
    decorrelation that `uplink_point_ahead_term` removes by design.
  - Two FAST cautions. The shipped `sim.aniso_servo_error` leaks
    `mask (1 - mask)` of the uncorrected band (0.061 rad^2 at a ZERO
    point-ahead angle, where the truth is 0). And the FAST grid misses 29 to
    48 percent of the whole-plane Kolmogorov residual (no support below `df`;
    the integral converges as `kappa^(1/3)`); the shipped Term's auto grid
    (`df` = 3.11 rad/m) misses more. The missing scales sit far above the
    aperture, so their effect on the coupled flux is damped; that damping is
    not quantified. OPEN follow-up.
  - The Term level: the fidelity-1 Monte Carlo mean reads 0.6 to 1.7 dB below
    the fidelity-0 pair at all five operating points (3.04 against 3.79 dB at
    AO(60), 60 deg — the backlog first reading, reproduced). The attribution
    ladder decomposes the gap into the mode-set convention, the auto-grid
    truncation, and the Marechal-against-Monte-Carlo mapping.
- **VERDICT.** MATCH — both routes are validated. At matched mode sets the
  PAOLA spatial-frequency filter and the Stone Zernike projection agree to
  about 5 percent across the swept angles, orders, and elevations, and the
  fitting sides agree to under 1 percent. The Term-level spread is a
  composition of measured conventions, not a physics disagreement. The FAST
  Term keeps its model-of-record role at the swept operating points.
- **Script.** `validation/fast_stone_pointahead/fast_stone_pointahead.py`
  (`--full` for the run of record); write-up
  [validation/fast_stone_pointahead/README.md](../validation/fast_stone_pointahead/README.md).
  See backlog 1-5.

### 9k. Does the fidelity-2 fade tail converge with the screen count, and what outer scale do the screens hold?

- **Question.** Section 9i converged the MEAN and the aperture index. The fade
  TAIL (p5, p1) sets the availability margin, and the post-WP7 re-test saw a
  2 dB p5 hint at 200 trials. Does the deep SMF fade tail converge as the
  near-ground `Cn2` is resolved with more screens? And, because the owner
  suspected the screen generator (a stack of N screens might lose more of the
  low-frequency band than one screen), does the generator itself depend on
  the count? Backlog 2-I2T, 2-N6, 2-P5.
- **Model under test.** The split-step layer of Section 7 end to end
  (`propagate_turbulent_scenario` through `Campaign`), on the 30 deg hero
  downlink (1550 nm, 500 km, a 0.7 m ground telescope with an UNCORRECTED SMF,
  no tip-tilt and no AO), with the GRID PINNED at 1024 px / 3.43 mm, because
  the sizer moves the grid with the screen count. Then the screen generator
  alone (`ScreenFactory`, no propagation), stacking the same per-screen `r0`
  lists.
- **Reference.** The tail has no analytic reference; the cases check each
  other with 68 percent bootstrap bars at 1000 trials (about +-0.7 dB at p5,
  +-1.5 dB at p1). The point index checks against the plane-wave Rytov
  variance of the slab (0.22). The screens check against Fried
  `D(r) = 6.88 (r/r0)^(5/3)`, DOI 10.1364/JOSA.56.001372, the Noll
  piston-removed and tilt-removed aperture variances 1.0299 and 0.134
  `(D/r0)^(5/3)`, DOI 10.1364/JOSA.66.000207, Table IV, and, at a finite outer
  scale, the same Noll filters integrated over the von Karman PSD (Schmidt,
  DOI 10.1117/3.866274, Ch. 9, Eq. (9.50)) with the Assemat and Wilson
  covariance, DOI 10.1364/OE.14.000988, Eq. (5).
- **Measured (2026-09-04, eight campaigns of 1000 trials).** The grid is a
  null: 512 against 1024 px on the same 9-screen plan agree at every quantile
  (p5 +0.27 dB, 0.2 sigma). On the pinned grid the SMF p5 falls 31.5 -> 31.0
  -> 31.2 -> 30.5 -> 29.0 dB from 5 to 25 screens (2.5 sigma), p10 falls
  1.3 dB, and the 9-screen plan with its 80 m ground screen cut into four
  equal-`Cn2` sub-screens lands on the 25-screen line (p5 29.8). p1 does not
  move (37.6 to 38.7 dB, 0.7 sigma). The POINT irradiance (one pixel) does not
  move at any count (quantile spread 0.4 dB; the index 0.23 to 0.26), and a
  re-bin of the stored fields shows the 512 and 1024 px grids agree exactly at
  a matched averaging area. `rapid` as shipped (256 px, 5 screens) reads
  inside the standard spread at 1/30 of the cost. The generator, phase only:
  against Kolmogorov (`L0 = inf`) every stack holds 0.75 to 0.84 of the Noll
  piston-removed variance and 1.00 of the tilt-removed variance, and the
  count changes it by -0.024 +-0.034 (5 -> 25 screens) and +0.015 +-0.049
  (the ground layer, 1 -> 4 screens). That deficit is the OUTER SCALE: three
  subharmonic levels reach 27 x the grid side (95 m), and the same screens
  match the von Karman theory at `L0 = 95 m` exactly (0.765 against 0.765).
  Drawn with `L0 = 25 m` and judged against the theory at 25 m, one screen
  reads 1.000 +-0.026 and a 5 / 9 / 25-screen plan 0.974 / 0.946 / 0.942.
- **VERDICT.** The default 9-screen plan is converged for the mean, the deep
  tail (p1) and every point quantity; its SMF p5 is about 2 dB PESSIMISTIC
  against 25 screens and the converged p5 is not known (40 screens not run).
  The count effect is the PHASE the fibre overlap pays (the ground-layer tilt
  and low orders), not scintillation, and it is NOT the generator, which
  misses the same fraction at every count with a mild few-percent stacking
  drift at a finite `L0`. The production `L0_m = inf` claims an outer scale
  the grid does not hold; the screens ARE a von Karman `L0 = 27 x side`
  atmosphere. The von Karman theory (the Noll piston filter over the von
  Karman PSD) gives a piston-removed aperture variance of 0.630 of the
  Kolmogorov value at `L0 = 25 m` and 0.765 at `L0 = 95 m`, for `D = 0.7 m`,
  so the fibre tilt that the SMF tail pays moves with the choice by an amount
  of the order of 2 dB at p5 (an ESTIMATE from the tilt share of the fade,
  not a measurement). Choose an explicit site `L0` and keep
  `L0 <= 27 x side` (backlog 2-P5, HIGH). NOT run: 20 deg, 40 screens.
- **Script.** `validation/tail_convergence/tail_convergence.py` and
  `validation/screen_stacking/screen_stacking.py [--L0 25]`; write-ups
  [validation/tail_convergence/README.md](../validation/tail_convergence/README.md)
  and
  [validation/screen_stacking/README.md](../validation/screen_stacking/README.md);
  the tracker note is in [schmidt-crosscheck.md](schmidt-crosscheck.md) after
  the post-WP7 measurement. See backlog 2-I2T, 2-N6, 2-P5.

---

## Source summary

- Andrews and Phillips, Laser Beam Propagation through Random Media, 2nd ed.
  (2005): scintillation index and aperture averaging (5b), Gaussian-beam coherence
  and Fried parameter (5e), residual and PSD (5f).
- Dios et al., Applied Optics 43 (2004) 3866: uplink coupled-flux and beam
  scintillation (5c, 5d), uplink Fried weight (5e).
- Churnside, Applied Optics 30 (1991) 1982: strong-turbulence aperture averaging
  (5b).
- Stone, Hu, Mills and Ma, J. Opt. Soc. Am. A 11(1), 347-357 (1994), DOI
  10.1364/JOSAA.11.000347: angular anisoplanatism (5g).
- Noll, J. Opt. Soc. Am. 66(3), 207-211 (1976): AO and tip-tilt residual (5f, 6a),
  the Zernike mode count (5g).
- Fried (1966): the Fried parameter (5f).
- Dikmelik and Davidson, Appl. Opt. 44(23), 4946-4952 (2005): analytic SMF coupling
  (6a).
- FAST (`fast-aosim`): statistical SMF coupling (6b).
- Schmidt, DOI 10.1117/3.866274: the split step, the Fourier phase screen, the
  absorbing boundary, and the grid rules (7); the numerical foundation layer
  `olb/waveoptics/schmidt/` (8).
- Shapiro, DOI 10.1364/JOSA.61.000492: the uplink reciprocity overlap (7).
- Martin and Flatte, DOI 10.1364/AO.27.002111: the pixel-per-coherence-length
  rule (7).

Where a section gives no literal DOI, the code cites the source in words next to
the equation. This document repeats those citations and names the file.
