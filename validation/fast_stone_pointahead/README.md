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
  decorrelation residual of the Zernike radial orders 1 to `max_order`. olb
  calls it through `olb.links.uplink.uplink_point_ahead_term`, which removes the
  PISTON only from 2026-09-11, so the tilt (order 1) stays in.
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

**The production Term keeps the TILT from 2026-09-11.**
`uplink_point_ahead_term` defaults to `remove="piston"`. The terminal senses the
tilt of the DOWNLINK beacon, and the steering mirror adds the point-ahead offset
geometrically. So the terminal holds no uplink tilt reference, and the uplink
pays the full tilt decorrelation. The old default `remove="piston_tilt"` assumed
a separate uplink tilt loop, and that loop does not exist in this design.

So the production pairing is FAST against the Stone band with the PISTON
removed (the column `F/S_pis`). The two sets then differ by the piston ALONE.
A piston is a constant phase over the pupil: it changes no overlap integral and
no far-field irradiance, so it is a bookkeeping difference and not a physics
difference. The script keeps the old piston-and-tilt column in the tables for
the record.

The Term also reads the site outer scale `Site.outer_scale_m` (25 m by
default). Stage C therefore gives FAST the same outer scale, through
`fast_params={"L0": 25.0}`, so the arm is like for like. The mode-matched
stages A and B keep the Kolmogorov limit, because the Stone paper writes it.

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
| C | The whole-Term comparison in dB, LIKE FOR LIKE at the site outer scale `L0` = 25 m on every leg. It gives the fidelity-1 FAST Monte Carlo mean against the fidelity-0 analytic pair, with the attribution ladder between them (Stone, FAST servo off, FAST default servo, each as an extended Marechal dB, T. S. Ross, DOI 10.1364/AO.48.001812). It also prints the mode-set table: the Stone value with no mode removed, with the piston removed, and with the piston and the tilt removed, at 25 m and in the Kolmogorov limit. |

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

**FULL run, 2026-09-11** (grid 1024 / 0.01 m, `df` = 0.614 rad/m, 3000 Monte
Carlo draws, runtime 604 s). Every numerics gate passes: the A0 single-layer
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

**The production pairing now differs by the PISTON alone.** From 2026-09-11 the
Term keeps the tilt, so both routes hold it. At the production point (1.5 m,
60 deg, 9.01 arcsec, ZMAX = 60), like for like at `L0` = 25 m:

| quantity | rad^2 |
| --- | --- |
| FAST clean split, keeps the piston and the tilt | 1.8088 |
| Stone, no mode removed (the mode-matched partner) | 1.7840 |
| Stone, piston removed (the shipped Term) | 0.9479 |
| the piston decorrelation | 0.8361 |
| the tilt decorrelation, which BOTH routes keep | 0.3902 |

So `F/S_pis` reads 1.91, and 0.836 of the 0.861 rad^2 gap is the piston. The
rest is the 1.4 % that separates FAST from the mode-matched Stone value on this
grid. A piston is a constant phase over the pupil, so it changes no overlap
integral and no far-field irradiance: the two routes agree on every mode that
the link can see.

**The outer scale cuts the piston, and it does not touch the tilt.** The same
point in the Kolmogorov limit reads 2.079 rad^2 of piston decorrelation and
0.407 rad^2 of tilt decorrelation. At `L0` = 25 m the piston falls to
0.836 rad^2 (a 60 % cut) and the tilt only to 0.390 rad^2 (a 4 % cut). The
piston difference between two directions lives in the largest scales, and the
finite outer scale removes them. So the shipped Term moves very little with the
outer scale (0.9679 to 0.9479 rad^2, 0.09 dB), while the mode-matched Stone
value falls from 3.047 to 1.784 rad^2.

That cut also closes the FAST grid truncation. At `L0` = 25 m the FAST study
grid holds the band: `F/S_n` = 1.014, against 0.638 in the Kolmogorov limit.

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
   misses 29 to 56 % of the whole-plane variance (38.6 % at the production
   point). A finite outer scale removes that caution: at `L0` = 25 m the same
   grid holds the band (stage C). The shipped `uplink_fast_term` runs on the
   FAST auto grid
   (NPXLS = 202, `df` = 3.11 rad/m), which misses more. The missing band sits
   at scales far above the 1.5 m aperture, where the phase across the pupil is
   close to a piston, so its effect on the COUPLED FLUX is damped; this study
   does not quantify that damping. OPEN follow-up.

**Stage C, the Term level: the gap GREW, and the reason is the tilt.** The
fidelity-1 Monte Carlo mean now reads 2.1 to 4.2 dB BELOW the fidelity-0
analytic pair at all five operating points (worst ratio 0.55, so the mechanical
verdict reads INVESTIGATE where it read MEASURED DIFFERENCE before).

| point | fid0 = point-ahead + fitting [dB] | fid1 Monte Carlo [dB] |
| --- | --- | --- |
| AO(60) @ 60 deg | 5.467 = 4.117 + 1.351 | 3.026 |
| AO(60) @ 30 deg | 9.470 = 7.130 + 2.340 | 5.317 |
| AO(21) @ 60 deg | 7.078 = 3.726 + 3.353 | 4.522 |
| AO(60) @ 90 deg | 4.735 = 3.565 + 1.170 | 2.655 |
| AO(10) @ 60 deg | 9.550 = 3.175 + 6.375 | 6.145 |

The fidelity-1 side did NOT move: 3.026 dB against the 3.04 dB of the
2026-09-02 run at AO(60), 60 deg, so the outer scale costs the FAST Monte Carlo
nothing. The fidelity-0 side moved from 3.79 to 5.47 dB, and the move is the
point-ahead Term alone: 2.435 dB (piston and tilt removed, Kolmogorov) to
4.117 dB (piston removed, `L0` = 25 m). The tilt adds 1.69 dB and the outer
scale gives 0.09 dB back.

The attribution ladder says where the remaining gap lives, and it is the LOW
FREQUENCY band. The FAST Monte Carlo runs on the shipped auto grid
(NPXLS = 202, `df` = 3.11 rad/m), which holds no scale above 2 m, so it drops
most of the tilt decorrelation that the analytic side now charges in full: the
default-servo residual reads 0.825 rad^2 (3.59 dB) against the 1.809 rad^2
(7.86 dB) of the same case on the fine study grid. The analytic side also maps
rad^2 to dB through the extended Marechal relation, which over-charges a
tilt-dominated error, because a tilt displaces the far-field beam and does not
scatter it. So the fidelity-0 pair is the PESSIMISTIC bound and the fidelity-1
Monte Carlo is the OPTIMISTIC one, and the tilt sits between them. This is an
OPEN follow-up, and it does not touch the mode-matched verdict above.

**Bottom line for backlog 1-5.** The comparison the backlog asked for is done,
and it VALIDATES both routes: at matched conditions the FAST kernel and the
Stone law agree to about 5 % across the swept point-ahead angles, corrected
orders, and elevations, and the fitting sides agree to under 1 %. With the
production mode set of 2026-09-11 the two routes differ by the piston only, and
no overlap integral sees a piston. The Term-level spread in dB is larger than
before, and it is the tilt: the analytic Term charges the full tilt
decorrelation, and the shipped FAST grid holds little of it.
