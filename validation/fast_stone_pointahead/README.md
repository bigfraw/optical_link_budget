# FAST against Stone: the point-ahead anisoplanatism residual (backlog 1-5)

## Question

An uplink terminal senses the turbulence on a downlink beacon. It applies the
conjugate phase to the uplink beam. The satellite moves during the round trip,
so the uplink goes to a different angle than the beacon came from. That
point-ahead angle decorrelates the correction, and a residual phase variance
stays.

olb holds TWO models of that residual:

- **Fidelity 1 (the model of record).** FAST (`fast-aosim` 0.1.7) builds the
  residual phase power spectrum `G_AO * Phi_n` and it integrates the corrected
  band. olb calls it through `olb.models.fast.uplink_fast_term`.
  Source: O. J. D. Farley and others, Opt. Express 30(13), 23050 (2022),
  DOI 10.1364/OE.458659.
- **Fidelity 0.** The Stone (1994) finite-aperture modal law
  (`olb.turbulence.anisoplanatism.anisoplanatic_phase_variance`) sums the
  decorrelation residual of the Zernike radial orders 2 to `max_order`. olb
  calls it through `olb.links.uplink.uplink_point_ahead_term`.
  Source: J. Stone, P. H. Hu, S. P. Mills and S. Ma, J. Opt. Soc. Am. A 11(1),
  347 (1994), DOI 10.1364/JOSAA.11.000347.

The two compute the SAME quantity. Do they give the same number? Where they do
not, what is the reason: the mode set, the numerics, or the physics?

The companion questions are the fitting error (FAST against Noll,
DOI 10.1364/JOSA.66.000207) and the whole-Term difference in dB.

## Method

The script drives FAST with the servo OFF: `TLOOP = 0`, `TEXP = 0`, and a zero
wind field. The PAOLA filter (`fast.ao_power_spectra.G_AO_PAOLA`) then reduces
to the pure anisoplanatic kernel

```
G_aniso(kappa) = 2 - 2 cos(delta_r . kappa),   delta_r_i = theta * h_i
```

which is the physics that Stone integrates. So the two routes see the same
atmosphere, the same aperture, and the same angle.

### The clean band split

FAST corrects modally (`MODAL=True`, `ZMAX=n`), so its low-frequency mask is a
SOFT Zernike Fourier filter `sum_j |Zhat_j(kappa)|^2` (piston included, clipped
at 1), not a hard cut. The shipped attribute `sim.aniso_servo_error` integrates
`G_ao * mask = aniso * mask^2 + mask (1 - mask)`, so it mixes a piece of the
UNCORRECTED band into the anisoplanatic number. The script rebuilds the kernel
from the sim attributes and it integrates each band ONE time:

```
aniso_corr = INT 2 pi k^2 Phi_n (2 - 2 cos(delta_r . kappa)) mask  d^2 kappa
fit_corr   = INT 2 pi k^2 Phi_n (1 - mask)                        d^2 kappa
```

The rebuilt kernel satisfies `G_ao = kernel * mask + (1 - mask)` to machine
precision, and `aniso_corr + fit_corr = sim.phs_var`. The script asserts both.

### The mode sets must match

The FAST modal mask is `sum_{j=1..ZMAX} |Zhat_j(kappa)|^2`, and the sum starts at
the Noll index 1. So the FAST corrected band KEEPS the piston and the two tilts.
The Stone set that holds the same modes is `remove='none'` over the band
`0..max_order`. That is the MODE-MATCHED pair, and the script reads its ratio
`Q/S_n` for the verdict.

The production Term `uplink_point_ahead_term` uses `remove='piston_tilt'`,
because a separate tracking loop points the beam. The script prints that pairing
too (`F/S_pt`), but it is NOT mode matched, so it is not the physics test.

### The low-frequency truncation

The Kolmogorov anisoplanatic integrand goes as `kappa^(-2/3)` at low frequency,
so the integral converges only as `kappa^(1/3)`. A FAST grid of `N` pixels and
pitch `dx` carries no frequency below `df = 2 pi / (N dx)`, so it MISSES a real
part of the variance. The script measures that part. It integrates the same
physics three ways:

1. on the FAST grid (the number FAST reports);
2. with an INDEPENDENT polar quadrature over the support that the grid holds
   (the square domain up to `kappa_max = pi/dx`, without the direct-current
   cell below `df/2`);
3. with the same quadrature over the WHOLE plane, which is what Stone
   integrates. The ratio of 2 to 3 is the truncation.

**The gate needs an outer scale.** In the Kolmogorov limit no grid can converge
on this integral, so a 1 % gate there has no meaning. The GATED legs (stage A0
and the convergence stage) therefore carry a finite outer scale `L0 = 5 m`. That
puts the low-frequency knee at `k0 = 2 pi / L0 = 1.26 rad/m`, well inside the
grid, so the FAST value and the independent quadrature must then agree. Every
other stage keeps the Kolmogorov limit, and the script reports the truncation
instead of gating it.

The quadrature shares no code with the FAST integrators. It builds the von
Karman spectrum `0.033 Cn2 dh (kappa^2 + k0^2)^(-11/6)` (Andrews and Phillips,
2nd ed. (2005), DOI 10.1117/3.626196, Ch. 3) and it reads the mask from
`fast.ao_power_spectra.zernike_squared_filter`.

### The stages

| Stage | What it does |
| --- | --- |
| A0 | One Gaussian Cn2 layer at 10 km, at the zenith. It is the cleanest test of the kernel, the mask, and the quadrature. It carries the 1 % gate. |
| A1 | The HV5/7 profile at the production point (1.5 m, 60 deg, the 600 km orbit point-ahead angle, ZMAX = 55). It prints the full attribution table, and it isolates the piston term and the tilt term. |
| B | Three sweeps: the point-ahead angle (0 to 2 times nominal), the corrected order (ZMAX 3 to 66, plus the production ZMAX = 60), and the elevation (30, 60, 90 deg). |
| Convergence | The A1 case on three grids: the base grid, a finer `df`, and a higher `kappa_max`. |
| C | The whole-Term comparison in dB. It gives the fidelity-1 FAST Monte Carlo mean against the fidelity-0 analytic pair, with the attribution ladder between them (Stone, FAST servo off, FAST default servo, each as an extended Marechal dB, T. S. Ross, DOI 10.1364/AO.48.001812). |

### Gates and verdicts

Hard asserts (the NUMERICS):

- A0: the FAST grid integral against the independent quadrature, on the same
  support and with the gate outer scale, within 1 %;
- every constructed sim: the clean closure
  `|aniso_corr + fit_corr - phs_var| / phs_var` below 0.5 %;
- convergence: `aniso_corr` moves less than 1 % between the grids, and
  `fit_corr` moves less than 1 % between the grids of one `kappa_max`. A grid
  of a DIFFERENT `kappa_max` is reported, not gated: a finer pixel opens a
  higher frequency band, and the Kolmogorov fitting tail there is real
  variance;
- a zero point-ahead angle gives a zero residual on both sides.

The PHYSICS agreement is NOT an assert. The script reports the mode-matched
ratio `Q/S_n` and it writes a mechanical verdict for each stage into the results
JSON: within 10 % is `MATCH`, 10 to 40 % is `MEASURED DIFFERENCE`, and past 40 %
is `INVESTIGATE`. Stage C reads the ratio of the fidelity-1 dB to the fidelity-0
dB.

## Run it

From the repository root:

```
python -m validation.fast_stone_pointahead.fast_stone_pointahead
python -m validation.fast_stone_pointahead.fast_stone_pointahead --full
```

Quick mode uses the 512 / 0.02 m grid and 500 Monte Carlo draws. Full mode uses
the 1024 / 0.01 m grid, 3000 draws, and every sweep point. The script writes
`fast_stone_pointahead.log`, `fast_stone_pointahead_results.json`, and three
figures in `figures/`.

The script needs `fast-aosim` (an optional, GPLv3 dependency). It reads the
production layer, and it changes no olb module.

Two notes about the current environment (numpy 2.4.6):

- numpy 2.4 removed the `np.trapz` alias, and olb still calls it. The script
  restores that one name at the top, in the validation script only.
- the matplotlib 3.11.1 build in this environment faults natively inside the
  Agg renderer AFTER the FAST libraries load into the process. So the script
  writes the JSON and the log first, and then it renders the figures from that
  JSON in a SEPARATE clean process (`--figures-only`). A render fault costs
  the figures only, never the numbers.

## Verdict

**FULL run, 2026-09-02** (grid 1024 / 0.01 m, `df` = 0.614 rad/m, 3000 Monte
Carlo draws, runtime 547 s). Every numerics gate passes: the A0 single-layer
FAST integral agrees with the independent quadrature to 0.53 % (gate 1 %), the
clean closure holds to machine precision, the refined-grid convergence moves
0.15 % (gate 1 %), and a zero point-ahead angle gives exactly zero on both
sides.

**MATCH — both routes are validated at the physics level.** The mode-matched
ratio (the whole-plane quadrature of the FAST soft Zernike filter against the
Stone band with no mode removed) reads `Q/S_n` = 1.044 to 1.055 across the
FULL sweep: the point-ahead angle 0.25 to 2 times nominal, the corrected order
ZMAX 1 to 66 (plus the exact-zero uncorrected anchor at ZMAX = 0), and the
elevation 30 to 90 deg. The single-layer case reads 0.991. The residual rises
steeply over the first corrected orders and then flattens: the low orders
decorrelate hardest over the point-ahead angle. So the PAOLA spatial-frequency filter and the Stone Zernike projection
— the two finite-aperture treatments the backlog asked about — differ by at
most about 5 %, and the difference grows gently with the angle. The fitting
side agrees to 0.6 %: the FAST clean uncorrected-band integral gives
0.3334 rad^2 against the Noll residual 0.3354 rad^2 at 55 modes.

**The production pairing differs by the MODE SET, not by the physics.** The
FAST modal mask keeps the piston and the two tilts, and
`uplink_point_ahead_term` removes them (a tracking loop points the beam). At
the production point that convention difference is 2.08 rad^2 of piston
decorrelation plus 0.41 rad^2 of tilt decorrelation, so the raw pairing reads
3.5x. That factor is explained, and it is not an error in either route.

**Two measured cautions on the FAST side.**

1. **The soft-mask leakage.** The shipped `sim.aniso_servo_error` mixes
   `mask (1 - mask)` of the uncorrected band into the anisoplanatic number: at
   a ZERO point-ahead angle it reads 0.061 rad^2 where the true residual is
   exactly 0. The clean band split of this script removes it. The leakage is
   small at the production point (2 %), and it grows as the true residual
   falls.
2. **The low-frequency grid truncation.** The Kolmogorov anisoplanatic
   integral converges only as `kappa^(1/3)`, and the FAST grid holds no
   frequency below `df`. On the fine 1024 / 0.01 study grid the FAST number
   misses 29 to 48 % of the whole-plane variance (38.6 % at the production
   point). The shipped `uplink_fast_term` runs on the FAST auto grid
   (NPXLS = 202, `df` = 3.11 rad/m), which misses more. The missing band sits
   at scales far above the 1.5 m aperture, where the phase across the pupil is
   close to a piston, so its effect on the COUPLED FLUX is damped; this study
   does not quantify that damping. OPEN follow-up.

**Stage C, the Term level: MEASURED DIFFERENCE, decomposed.** The fidelity-1
Monte Carlo mean reads 0.6 to 1.7 dB BELOW the fidelity-0 analytic pair at
all five operating points (worst ratio 0.78). The backlog first reading is
reproduced: 3.04 dB against 3.79 dB at AO(60), 60 deg, 1.5 m. The attribution
ladder shows where the gap lives: the two sides hold different mode sets
(above), the default-servo FAST sim runs on the coarse auto grid, and the
analytic side maps rad^2 to dB through the extended Marechal relation while
the Monte Carlo measures the true mode overlap. The gap is a composition of
known, measured conventions, not an unexplained physics disagreement.

**Bottom line for backlog 1-5.** The comparison the backlog asked for is done,
and it VALIDATES both routes: at matched conditions the FAST kernel and the
Stone law agree to about 5 % across the swept point-ahead angles, corrected
orders, and elevations, and the fitting sides agree to under 1 %. The
Term-level spread is explained by convention, grid, and mapping choices, each
one measured here.
