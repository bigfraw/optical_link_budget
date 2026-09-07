# receiver_cone_clip

**A converging absorbing boundary that follows the receiver cone. Backlog 2-P3,
route (b). VERDICT: NO. Neither cone rule holds the turbulent field at any
factor of the sweep, and where the UNMAPPED rule holds the budget-visible
numbers it gives no extent saving at 5 km and a NEGATIVE one at 10 km.**

## The claim under test

The owner flagged this on 2026-09-06. The turbulent sizer holds the WHOLE beam
plus the scatter cone all the way to the receiver. But only the light inside
the BACK-PROJECTED cone of the receive aperture can land on that aperture. If
that is true, an absorbing boundary that converges onto the cone leaves the
field INSIDE the aperture unchanged and it removes most of the grid extent. A
smaller extent is a finer pixel at the same pixel count, which is what a 10 km
terrestrial cell needs: the sizer there asks for a 9.2 m side and it clamps the
pixel count at 2048.

## The method

The grid does NOT change. One fixed grid, one fixed screen plan, and, for each
trial, ONE set of phase screens. Only the MASK differs between the two arms, so
every difference inside the receive aperture is the doing of the clip.

- The REFERENCE arm rebuilds the production `split_step` loop by hand and it
  ASSERTS `numpy.array_equal` against the public `propagate_turbulent_field`.
  The log prints that check.
- The CLIPPED arm runs the SAME loop on the SAME screens, with the edge mask
  multiplied by a converging super-Gaussian `exp(-(r/R(z))^16)` at every plane
  where the production loop applies the edge mask.
- A VACUUM arm (no screens) runs the same masks.

The script now offers TWO cone rules, with `L` the path, `D` the receive
diameter, `c` the swept factor, and `delta(z) = 3 lambda (L - z) / D` the ramp
allowance of both:

    MAPPED    (--cone-rule mapped)
    R(z) = c * [ D/2 + theta_s(z) * (L - z) ] * (z / L) + delta(z)
    theta_s(z) = lambda / r0_rest,  r0_rest = (0.423 k^2 Cn2 (L - z))^(-3/5)

    UNMAPPED  (--cone-rule unmapped, the default)
    R(z) = c * [ (D/2) * (z / L) + theta_s * (L - z) ] + delta(z)
    theta_s     = lambda / r0_full,  r0_full = (0.423 k^2 Cn2 L)^(-3/5)

    the guard (--guard on, off by default)
    R(z) -> max( R(z), 3 w(z) )

`r0` is the plane-wave Fried parameter (Andrews and Phillips, 2nd ed. (2005),
DOI 10.1117/3.626196, Ch. 6, Eq. (64)). The `z/L` factor is the far-field ray
map: past the Rayleigh range a ray at `x` at the plane `z` arrives at `x L / z`,
so it must obey `x <= (D/2)(z/L)` to land inside the aperture (Ch. 4, Eqs. (7)
and (8), printed p. 87). The MAPPED rule puts that map on the SCATTER term too.
The UNMAPPED rule does not, because a screen sends scattered light away at a
fresh angle and that light crosses `theta_s (L - z)` wherever `z` sits; the
UNMAPPED rule also reads `theta_s` from the WHOLE path, because the light at the
plane `z` can already carry the scatter of the path behind it. `3 lambda (L-z)/D`
is the ramp allowance, so the diffraction cone of the mask itself stays inside
the aperture cone. The mask family is the absorbing operator of Schmidt (2010),
DOI 10.1117/3.866274, Ch. 8, Eq. (8.18), printed p. 139, whose shape is
Eq. (8.1), printed p. 134.

**The guard is OFF by default.** On a 5 km path a 5 mm waist grows to
`w(L) = 0.49 m`, so `3 w(L) = 1.48 m`, which is LARGER than the 1.43 m half-side
of the sizer grid. With the guard on the mask therefore covers the whole grid at
the receive end, it removes 5e-4 of the launched power (the far tail near the
launch, not the cone), and the clip is VACUOUS: it can save no extent, so a PASS
there says nothing about the cone. `--guard on` keeps it for the record, and
`--cone-rule mapped --guard on` reproduces the first run below, row for row.

## The cells

`validation.terrestrial_campaigns.run_campaigns.build_scenario(L, 1e-14,
"collimated")`: 1550 nm, a 5 mm collimated waist, a 10 cm receive aperture with
an optimal-focus single-mode fibre. Preset `rapid`, `L0_m = 25.0`, precision
`single`, seed 20260906, 6 trials. Two path lengths:

| Cell | Grid | Side | Pixel | Screens | r0 of the path |
| --- | --- | --- | --- | --- | --- |
| 5 km, Cn2 = 1e-14 | 1024 px | 2.8576 m | 2.791 mm | 10 | 2.99 cm |
| 10 km, Cn2 = 1e-14 | 1024 px | 6.5288 m | 6.376 mm | 35 | 1.97 cm |

The 9.2 m side and the 2048 px clamp that motivated route (b) belong to the
10 km cell at the `standard` preset. At `rapid`, which every run here uses, the
sizer asks for 6.5288 m at 1024 px.

## The pass rule

For one factor `c`, EVERY trial must give, inside the 10 cm aperture: a
relative field RMS below 1e-3, a collected-power step below 0.01 dB, and a
fibre-efficiency step below 1e-3. 1e-3 sits well ABOVE the single-precision
floor (`validation/precision` measured a field RMS of 1.3e-6 between a
complex64 and a complex128 run of the same seed), and well BELOW anything a
budget can see (0.01 dB is under the 0.11 dB screen-count spread of WP7).

The script reports the smallest passing factor TWO ways, because the two say
different things:

- the FIELD rule: the 1e-3 field RMS alone. It asks whether the clipped field
  IS the reference field, so it is the honest test of the claim.
- the BUDGET-VISIBLE rule: the 0.01 dB power step and the 1e-3 efficiency step
  alone. That is all a Term reads, so it is the test of whether a budget could
  ever tell the two runs apart.

## Run 1: the MAPPED rule, 5 km (2026-09-06)

    TURBULENCE, the 3 w(z) guard ON (the rule as written), 6 trials
        c    removed     max RMS  max |dP| dB  max |d eta|  verdict
     1.00  5.284e-04   1.066e-03    3.995e-05    1.768e-06     FAIL
     1.50  5.259e-04   1.041e-03    4.147e-05    1.429e-06     FAIL
     2.00  5.250e-04   1.031e-03    4.185e-05    1.291e-06     FAIL
     3.00  5.150e-04   9.975e-04    4.604e-05    1.213e-06     PASS
     5.00  4.895e-04   9.267e-04    4.984e-05    8.355e-07     PASS

    VACUUM, the 3 w(z) guard ON, 1 trial
     every c: removed 7.212e-06, RMS 0, dP 0, d eta 0            PASS

    TURBULENCE, no guard (the cone alone), 6 trials
        c    removed     max RMS  max |dP| dB  max |d eta|  verdict
     1.00  9.913e-01   2.298e-01    9.128e-01    2.730e-02     FAIL
     1.50  9.782e-01   7.862e-02    2.841e-02    1.348e-03     FAIL
     2.00  9.594e-01   5.680e-02    1.546e-02    9.606e-04     FAIL
     3.00  8.939e-01   3.297e-02    6.628e-03    3.938e-04     FAIL
     5.00  7.195e-01   1.234e-02    1.242e-03    5.583e-05     FAIL

    VACUUM, no guard, 1 trial
        c    removed     max RMS  max |dP| dB  max |d eta|  verdict
     1.00  9.824e-01   1.742e-01    6.507e-01    4.505e-02     FAIL
     1.50  9.609e-01   3.589e-04    1.419e-03    1.177e-04     PASS
     2.00  9.316e-01   3.597e-06    1.438e-05    1.254e-06     PASS
     3.00  8.527e-01   0.000e+00    0.000e+00    0.000e+00     PASS
     5.00  6.433e-01   0.000e+00    0.000e+00    0.000e+00     PASS

The implied grid, if the mask radius set the side at every plane
(`side = 2 max_z R(z) / (1 - boundary_width_frac)`, against the 2.8576 m of the
sizer, and the pixel count that keeps the 2.791 mm pixel of today):

    guard    c     max_z R(z)   side      factor   n at the same pixel
    on       any     1.4802 m   3.2894 m   0.87    2048   (today 1024)
    off      1.0     0.2362 m   0.5248 m   5.44     256
    off      1.5     0.2555 m   0.5677 m   5.03     256
    off      2.0     0.2814 m   0.6253 m   4.57     256
    off      3.0     0.3406 m   0.7568 m   3.78     512
    off      5.0     0.4687 m   1.0415 m   2.74     512

## Run 2: the UNMAPPED rule (2026-09-06)

The same grid, the same plan, the same seed and the same screens. Only the cone
rule changes. The files are `receiver_cone_clip_unmapped_off_5km.*`,
`receiver_cone_clip_unmapped_on_5km.*` and
`receiver_cone_clip_unmapped_off_10km.*`.

### 5 km, no guard (the cone alone)

    TURBULENCE, 6 trials
        c    removed     max RMS  max |dP| dB  max |d eta|  centre ratio  verdict
     1.00  9.914e-01   2.176e-01    9.849e-01    2.589e-02      0.971473    FAIL
     1.50  9.781e-01   4.654e-02    1.122e-02    7.322e-04      1.009474    FAIL
     2.00  9.593e-01   3.028e-02    6.183e-03    2.913e-04      1.010767    FAIL
     3.00  8.937e-01   1.076e-02    1.246e-03    5.968e-05      1.008528    FAIL
     5.00  7.190e-01   1.236e-03    5.628e-05    3.655e-06      1.000336    FAIL

    VACUUM, 1 trial
        c    removed     max RMS  max |dP| dB  max |d eta|  centre ratio  verdict
     1.00  9.824e-01   1.742e-01    6.507e-01    4.505e-02      1.000000    FAIL
     1.50  9.609e-01   3.589e-04    1.419e-03    1.177e-04      1.000000    PASS
     2.00  9.316e-01   3.597e-06    1.438e-05    1.254e-06      1.000000    PASS
     3.00  8.527e-01   0.000e+00    0.000e+00    0.000e+00      1.000000    PASS
     5.00  6.433e-01   0.000e+00    0.000e+00    0.000e+00      1.000000    PASS

    smallest passing c, FIELD rule            none of 1 to 5
    smallest passing c, BUDGET-VISIBLE rule   c = 2

### 5 km, the 3 w(z) guard ON (the clip is vacuous, kept for the record)

    TURBULENCE, 6 trials
        c    removed     max RMS  max |dP| dB  max |d eta|  centre ratio  verdict
     1.00  5.069e-04   9.466e-04    4.794e-05    1.023e-06      1.000330    PASS
     1.50  4.865e-04   8.743e-04    5.174e-05    9.102e-07      1.000258    PASS
     2.00  4.835e-04   8.266e-04    3.576e-05    1.032e-06      1.000274    PASS
     3.00  4.550e-04   7.132e-04    2.435e-05    7.516e-07      1.000271    PASS
     5.00  3.856e-04   4.031e-04    9.512e-06    4.964e-07      1.000099    PASS

    VACUUM, 1 trial
     every c: removed 7.212e-06, RMS 0, dP 0, d eta 0             PASS

    smallest passing c, FIELD rule            c = 1
    smallest passing c, BUDGET-VISIBLE rule   c = 1

### 10 km, no guard (35 screens, 6.5288 m side, 6.376 mm pixel)

    TURBULENCE, 6 trials
        c    removed     max RMS  max |dP| dB  max |d eta|  centre ratio  verdict
     1.50  9.963e-01   7.764e-02    9.926e-02    5.831e-04      1.066173    FAIL
     2.00  9.931e-01   3.668e-02    5.931e-02    4.282e-04      1.027782    FAIL
     3.00  9.860e-01   1.607e-02    1.568e-02    3.613e-04      0.998386    FAIL
     5.00  9.477e-01   5.607e-03    3.502e-03    1.775e-04      1.003710    FAIL

    VACUUM, 1 trial
        c    removed     max RMS  max |dP| dB  max |d eta|  centre ratio  verdict
     1.50  9.916e-01   3.807e-04    1.477e-03    1.209e-04      1.000000    PASS
     2.00  9.851e-01   3.817e-06    1.490e-05    1.229e-06      1.000000    PASS
     3.00  9.668e-01   0.000e+00    0.000e+00    0.000e+00      1.000000    PASS
     5.00  9.108e-01   0.000e+00    0.000e+00    0.000e+00      1.000000    PASS

    smallest passing c, FIELD rule            none of 1.5 to 5
    smallest passing c, BUDGET-VISIBLE rule   c = 5

### The extent that each factor implies

`side = 2 max_z R(z) / (1 - boundary_width_frac)`, and the pixel count that
keeps the pixel of the sizer grid:

    cell / guard        c    max_z R(z)   side      sizer/side   n at the same pixel
    5 km, off         1.0      0.4919 m   1.0930 m     2.61        512
    5 km, off         1.5      0.6215 m   1.3812 m     2.07        512
    5 km, off         2.0      0.7512 m   1.6694 m     1.71       1024   (today 1024)
    5 km, off         3.0      1.0106 m   2.2457 m     1.27       1024
    5 km, off         5.0      1.5293 m   3.3985 m     0.84       2048
    5 km, on        1 to 3     1.4802 m   3.2894 m     0.87       2048
    5 km, on          5.0      1.5293 m   3.3985 m     0.84       2048
    10 km, off        1.5      1.6444 m   3.6541 m     1.79       1024   (today 1024)
    10 km, off        2.0      2.0375 m   4.5277 m     1.44       1024
    10 km, off        3.0      2.8237 m   6.2749 m     1.04       1024
    10 km, off        5.0      4.3962 m   9.7693 m     0.67       2048

## What it means

**1. The MAPPED rule as written saves nothing.** The `3 w(z)` guard is larger
than the converging cone from about 2.3 km onward, and at the receive plane it
is 1.48 m against a 1.43 m half-side. So `max_z R(z)` is 1.48 m at EVERY factor,
the implied side (3.29 m) is LARGER than the side the sizer already picks
(2.86 m), and the pixel count would go UP, not down. That is why the guard is
now off by default: a PASS with the guard on is a pass of a mask that touches
almost nothing.

**2. The ray-cone claim IS correct in VACUUM, under both rules.** With the guard
off the cone throws away 93 to 99.6 percent of the launched power and the vacuum
field inside the aperture is untouched from `c = 1.5` up (a field RMS of 3.6e-4
at `c = 1.5`, 3.6e-6 at `c = 2`, and exactly zero at `c = 3` and 5, at 5 km and
at 10 km alike). So the geometry of the argument holds: a collimated beam in the
far field really does send only the light inside `(D/2)(z/L)` onto the aperture.

**3. The UNMAPPED rule is much better than the MAPPED rule, and it still FAILS
the field rule at every factor of the sweep.** At 5 km the turbulent field RMS
falls from 0.218 to 1.24e-3 over `c = 1` to 5, against the mapped rule's 0.230
to 1.23e-2: a factor of TEN better at `c = 5`, and it misses the 1e-3 limit by
24 percent. At 10 km it falls from 7.8e-2 to 5.6e-3, six times the limit. So the
shape correction is real and it is not enough.

**4. The remaining failure is the TAIL of the scattered light, not the shape.**
`theta_s = lambda / r0` is a COHERENCE angle, an RMS width. The Kolmogorov
angular spectrum has a power-law wing, so a fixed multiple of `theta_s` never
holds every ray that reaches the aperture; it holds a quantile of them, and the
field RMS measures exactly the part it misses. The measured convergence says the
same thing: the field RMS falls as about `c^-4.2` at 5 km and `c^-2.1` at 10 km,
a power law, not a cut-off. Extrapolating from `c = 3` and 5, the field rule
needs about `c = 5.3` at 5 km and about `c = 11.5` at 10 km.

**5. So the saving is NOT real at the cell that needs it.** Read the extent
table against the sizer side:

| Cell | Rule that HOLDS THE FIELD | Rule that holds the BUDGET numbers |
| --- | --- | --- |
| 5 km, guard off | none of `c = 1` to 5; about `c = 5.3` extrapolated, side 3.55 m, 1.24 x the sizer's 2.86 m | `c = 2`, side 1.67 m, 1.71 x smaller than the sizer, but still 1024 px at the same pixel |
| 5 km, guard on | `c = 1`, but the clip is vacuous (it removes 5e-4 of the power); side 3.29 m, LARGER than the sizer | the same vacuous `c = 1` |
| 10 km, guard off | none of `c = 1.5` to 5; about `c = 11.5` extrapolated, side 21.2 m, 3.25 x LARGER than the sizer's 6.53 m | `c = 5`, side 9.77 m, 1.5 x LARGER than the sizer |

At 5 km the budget-visible factor `c = 2` does shrink the extent by 1.71 x, but
`2.8576 / 2.791 mm = 1024` px and `1.6694 / 2.791 mm = 598` px both round up to
1024, so the pixel count does not fall. At 10 km, the cell that motivated route
(b), EVERY factor that a budget cannot see already asks for MORE extent than the
sizer gives. The cone that holds the aperture field is wider than the beam.

**6. The verdict, in one line.** Route (b) is closed for both rules. The
UNMAPPED cone is the right shape and the honest rule, and on this hardware it
buys nothing: the field rule costs more extent than it saves, and the only
factors that save extent are the ones a budget can only just not see.

## What a reader should take from the budget-visible column

The budget-visible numbers pass at `c = 2` (5 km) and `c = 5` (10 km) while the
field is plainly wrong. That is not a licence to clip. The metric that a Term
reads is a power and a mode overlap, both of them integrals over the aperture,
and an integral hides a field error that a fade tail does not: the pass rule
uses SIX trials, so it measures the middle of the distribution, not p1 or p5.
A clip that moves the field by 3 percent RMS is a clip whose deep-fade
statistics are untested.

## Run it

    python -m validation.receiver_cone_clip.receiver_cone_clip
    python -m validation.receiver_cone_clip.receiver_cone_clip --guard on
    python -m validation.receiver_cone_clip.receiver_cone_clip         --cone-rule mapped --guard on
    python -m validation.receiver_cone_clip.receiver_cone_clip         --path-km 10 --cn2 1e-14 --factors 1.5 2 3 5

The default is `--cone-rule unmapped --guard off`. Every output file carries the
rule, the guard and the cell, so a run does not write on the record of another
run. About 90 s for a five-factor 5 km sweep and about 210 s for a four-factor
10 km sweep, on a laptop.

## The files

| File | Content |
| --- | --- |
| [receiver_cone_clip.py](receiver_cone_clip.py) | The script. It rebuilds the production `split_step` loop by hand and asserts `numpy.array_equal` against `propagate_turbulent_field`, then reruns the same screens with the converging mask. Options: `--cone-rule mapped|unmapped`, `--guard on|off`, `--trials`, `--factors`, `--preset`, `--seed`, `--path-km`, `--cn2`. |
| [receiver_cone_clip.log](receiver_cone_clip.log) | Run 1, the MAPPED rule at 5 km, both guards. The untagged name of the first run. |
| [receiver_cone_clip_results.json](receiver_cone_clip_results.json) | The per-trial record of run 1. |
| [receiver_cone_clip_mapped_on_5km.log](receiver_cone_clip_mapped_on_5km.log) | The one-row reproduction check of run 1 (`c = 3`), row for row the same. |
| [receiver_cone_clip_unmapped_off_5km.log](receiver_cone_clip_unmapped_off_5km.log) | Run 2a: the UNMAPPED rule, 5 km, no guard. |
| [receiver_cone_clip_unmapped_on_5km.log](receiver_cone_clip_unmapped_on_5km.log) | Run 2b: the UNMAPPED rule, 5 km, the guard on. |
| [receiver_cone_clip_unmapped_off_10km.log](receiver_cone_clip_unmapped_off_10km.log) | Run 2c: the UNMAPPED rule, 10 km, no guard. |
| `receiver_cone_clip_<rule>_<guard>_<cell>.json` | Every per-trial metric, the pass flags, the two smallest passing factors, and the `R(z)`, `w(z)` profiles of that run. |
| [figures/](figures/) | One figure for each run. Left: the field RMS against c. Middle: the power step against c. Right: `R(z)` against `w(z)` and the grid half-side. |
