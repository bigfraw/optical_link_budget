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
  `WEAK_FLUCTUATION_LIMIT = 0.25`. Above it, focusing and saturation make the
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
- When the mean log-amplitude variance `sigma2_x` exceeds
  `WEAK_FLUCTUATION_LIMIT = 0.25`, the scintillation approaches saturation and the
  numbers are not trustworthy. The code carries the flag and warns.
- The launch beam is an untruncated Gaussian of waist `w0`. The code models no
  launch aperture and no central obscuration, so the fade does not change with an
  obscured pupil. `uplink_turbulence_term` flags a set obscuration in
  `budget.check()`. The MEAN loss from a central obscuration is separate and IS
  carried: the launch-truncation Term (`tx_gaussian_efficiency_term`, Section 2a)
  reads the obscuration and matches the wave-optics far-field to about 2 dB. The
  size of the obscuration effect on the FADE is UNRESOLVED: an earlier validation
  compared this index against the fidelity-2 reciprocity overlap, but those two do
  not agree even with no obscuration, so that comparison is void. See the
  investigation note.
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

- Weak-to-moderate turbulence. The model has no saturation.
- Dios reports good agreement with a split-step reference up to
  `sigma2_chi ~ 0.6`. Above that the true index saturates and the model overshoots.

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
  UNIFORMLY illuminated aperture and a FLAT (best-focus) wavefront. It holds when
  the received spot overfills the aperture and the receiver focuses for the
  incoming curvature. A near-field terrestrial link inside the Rayleigh range can
  break both: the received Gaussian tapers across the aperture, and the wavefront
  is curved. A refocus removes the curvature. The residual taper error runs SAFE,
  because a Gaussian-into-Gaussian overlap can pass the 0.8145 top-hat value, so
  the constant is then CONSERVATIVE. A curvature-aware, illumination-aware
  `eta_max` is the open Gap-3 upgrade. See `olb/models/coupling/_common.py`.
- The code flags an effective D/r0 above `SMF_DEEP_TURBULENCE_DR0 = 10`, where the
  practical coupling curve is extrapolated.
- The terrestrial form adds the effective-r0 weak-turbulence caveat: it evaluates
  plane-wave, Kolmogorov, phase-only forms at the Gaussian-beam r0. It ignores
  beam-wave amplitude scintillation, beam wander, and near-field curvature.

#### Source

- Y. Dikmelik and F. M. Davidson, "Fiber-coupling efficiency for free-space optical
  communication through atmospheric turbulence," Appl. Opt. 44(23), 4946-4952
  (2005), for the uncorrected coupling curve and its limits.
- Extended Marechal approximation (Chan and others), for the small-residual
  limit. Derivation and validity: T. S. Ross, Appl. Opt. 48(10), 1812 (2009),
  DOI 10.1364/AO.48.001812.
- Noll 1976, for the residual variance (Section 5f).

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
  moving satellite is not modelled. The Term flags this.
- The default scales (L0 = inf, l0 = 1 um) are the Kolmogorov limit. A finite L0 or
  a large l0 from `fast_params` makes it a von Karman spectrum, and the Term reads
  the label from the resolved scales.
- FAST models the phase with real Monte-Carlo screens (the phase-driven coupling
  fade is fidelity-1). It models the log-amplitude as an aperture-averaged
  log-normal, which holds only in the weak fluctuation regime. The code flags when
  the plane-wave amplitude `sigma2_I` exceeds `WEAK_FLUCTUATION_LIMIT = 0.25`. A
  deep coupled-power fade does not trip that flag, because that fade is phase-driven
  and modelled correctly.
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

The static mode match sets the flat-wavefront coupling `eta_max`. For a
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
focal spot (`w_s`) and the fibre mode (`w_m`):

    eta(dx)/eta_max = exp( -2 * dx^2 / w_eff^2 ),   w_eff^2 = w_s^2 + w_m^2

The two axes of the offset are i.i.d. Gaussian, so the loss in dB is exponential
with mean `(20/ln10) * f^2 * sigma2_theta / w_eff^2`. The loss grows without a
limit as the tilt grows.

A single-mode-fibre subtlety: at a fixed `a` the focal length cancels in this
mean, because `w_s` scales with `f`. So `f` sets `eta_max` through `a`, but it
does not change the angular sensitivity on its own.

No double-count with the mean coupling Term: the tip-tilt appears once. When the
walk-off Term is active, the mean coupling Term (Section 6a) keeps the
higher-order residual only. A virtual tip-tilt removes the Noll tip-tilt from its
residual (`drop_tiptilt=True`). The walk-off Term owns the tip-tilt.

A multimode fibre is a light bucket: the core is a HARD disk of radius `a_core`.
It collects ALL the spot power inside the core, so the coupling is the encircled
energy of the displaced Gaussian spot, NOT a mode overlap. With the spot centre
at offset `dx = f*theta`:

    eta(dx) = 1 - Q1( 2*dx/w_s ,  2*a_core/w_s )       (Marcum Q-function)

At `dx = 0` this reduces to `eta_static = 1 - exp(-2*a_core^2/w_s^2)`. A small
spot deep inside the core loses nothing until it nears the edge (a flat-top
acceptance); at the core edge it collects about half the power (about 3 dB). The
Term averages the loss over the Rayleigh offset. This differs from a single-mode
fibre, whose acceptance is a Gaussian mode, not a hard disk.

`optimal_focus` derives the focal length. For a single-mode fibre it picks `f` so
`a = 1.12` (the eta_max peak). For a multimode fibre it matches the spot to the
core, `a_core/w_s = 1.12`, which gives about 92% static capture.

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
  `pointing_jitter_rad`, compensation stack), the transmit waist, the path length,
  and the constant Cn2.
- Output: an SMF walk-off Term (category `pointing`), or an MMF coupling Term
  (category `coupling`). Both carry a real fade.

#### Assumptions and limits

- The walk-off falloff uses a Gaussian fit to the focal spot, and the eta_max
  value uses the more exact Airy-to-Gaussian overlap. The two spot models differ,
  which is standard practice near the peak.
- Contribution B of the received tip-tilt is the beam-wander tilt (A) only. The
  aperture angle-of-arrival "corrugation" tilt is available but feeds no Term:
  `aperture_arrival_angle_variance` in `olb/turbulence/angle_of_arrival.py` now
  delegates to `andrews.structure.angle_of_arrival_variance` (the gradient-tilt
  form, C-04), but no coupling Term adds contribution B. So the received tip-tilt
  is a lower bound. See `docs/andrews-crosscheck.md` batch 2 and backlog 0-W3.
- The beam-wander tilt is a weak-fluctuation model, so `sigma2_theta` is valid in
  weak turbulence only. The walk-off mapping itself has no upper limit.
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
- The random draw comes from `aotools` (LGPL-3.0). `olb` imports it as the
  optional extra `screens`. `olb` does not copy it.

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
