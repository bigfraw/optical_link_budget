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

File: `olb/models/transmittance.py`

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
no literal DOI (source cited in `olb/models/transmittance.py`).

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

`get_c2n` is the shared kernel from `olb._deps`. It reads the site RMS wind and
the ground-level Cn2. Use this profile when the optional `fast` package is not
available. The `fast` HV57 path fails without that package.

#### Inputs and outputs

- Inputs: the site (RMS wind, ground Cn2), the altitude grid.
- Output: the zenith Cn2(h) profile on the grid.

#### Source

The profile builder is the shared `get_c2n` kernel (source cited in
`olb/turbulence/profiles.py` and `olb/_deps.py`).

### 5b. Scintillation index and aperture averaging (downlink plane wave)

File: `olb/turbulence/scintillation.py`

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
  lognormal model depart from data.
- The aperture filter `(2*J1(x)/x)^2` assumes a uniform circular aperture with no
  central obscuration. An annular aperture is not modelled yet.
- The Kolmogorov spectrum has no inner scale and no outer scale.

#### Source

- Andrews and Phillips, Laser Beam Propagation through Random Media, 2nd ed.
  (2005): Eq. (12.44) for the plane-wave index; Ch. 10 for the aperture-averaging
  filter and the closed-form factors; Ch. 5, 6, 9 for the single-path forms.
- Churnside, Applied Optics 30 (1991) 1982, for the strong-turbulence
  aperture-averaging factor.

### 5c. Uplink coupled-flux Monte Carlo

File: `olb/turbulence/coupled_flux.py`

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
  `WEAK_FLUCTUATION_LIMIT = 0.6`, the scintillation approaches saturation and the
  numbers are not trustworthy. The code carries the flag and warns.
- The coupled-flux MC needs the `fast` package to build the HV57 Cn2 profile, or
  an explicit `cn2_profile`.

#### Source

Dios et al., Applied Optics 43 (2004) 3866, for the coupled-flux and wander-offset
mechanism. The shared kernels come from `olb._deps`.

### 5d. Dios on-axis and off-axis beam scintillation

File: `olb/turbulence/beam_scintillation.py`

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
> terrestrial single-mode-fibre coupling in `olb/models/coupling.py`, does not
> pass `f0`, so it uses the collimated default too.
>
> To generalise, stop pinning `Theta0` to 1. For the profile form, pass a finite
> `f0` (a diverged beam has `f0 < 0`, so `Theta0 > 1`; a focused beam has
> `0 < f0`, so `Theta0 < 1`). For the single-path form, add a curvature argument
> and thread it through `output_beam_params`, the way the profile form already
> does. This matters for a deliberately diverged uplink beam: the package already
> recasts that beam through a virtual waist in `olb/beam.py` for the geometric and
> the scintillation Terms, but the Fried-parameter feed still assumes a collimated
> input.

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
  pre-compensation budget that uses this Term (`uplink_budget` with a
  `DownlinkBeacon` source) is therefore phase-only. It MISSES the scintillation
  and understates the deep fade. This is a major known gap. See `api-budget.md`.

#### Source

- J. Stone, P. H. Hu, S. P. Mills and S. Ma, "Anisoplanatic effects in
  finite-aperture optical systems," J. Opt. Soc. Am. A 11(1), 347-357 (1994).
  DOI: 10.1364/JOSAA.11.000347.
- R. J. Noll, J. Opt. Soc. Am. 66(3), 207-211 (1976), DOI 10.1364/JOSA.66.000207,
  for the Zernike mode count in `max_radial_order`.

---

## 6. Fibre coupling

### 6a. Analytic mean coupling (fidelity 0)

File: `olb/models/coupling.py`

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
  budget reports no fade margin. Use `smf_fidelity="fast"` for the fade.
- The Dikmelik-Davidson coupling assumes a uniform circular aperture with no
  central obscuration. The code flags an obscured receive aperture.
- The code flags an effective D/r0 above `SMF_DEEP_TURBULENCE_DR0 = 10`, where the
  practical coupling curve is extrapolated.
- The terrestrial form adds the effective-r0 weak-turbulence caveat: it evaluates
  plane-wave, Kolmogorov, phase-only forms at the Gaussian-beam r0. It ignores
  beam-wave amplitude scintillation, beam wander, and near-field curvature.

#### Source

- Y. Dikmelik and F. M. Davidson, "Fiber-coupling efficiency for free-space optical
  communication through atmospheric turbulence," Appl. Opt. 44(23), 4946-4952
  (2005), for the uncorrected coupling curve and its limits.
- Extended Marechal approximation (Chan and others), for the small-residual limit.
- Noll 1976, for the residual variance (Section 5f).

### 6b. FAST statistical coupling (fidelity 1)

File: `olb/models/coupling_fast.py`

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
from the plane-wave scintillation index (Section 5b, Andrews and Phillips, 2nd ed.,
Eq. 12.44).

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

Where a section gives no literal DOI, the code cites the source in words next to
the equation. This document repeats those citations and names the file.
