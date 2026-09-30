# The velocity aberration of a corner-cube return

Date: 2026-09-30. Branch `retro-velocity-aberration`. Machine: the laptop.

## The purpose

This study records the physics and the numbers of
`olb.links.retro_space.retro_velocity_aberration_term`, the Term that the
retro budget charges on EVERY rung (0, 1 and 2).

## The physics

A corner cube sends the light back along the incoming ray IN ITS OWN FRAME:
in the frame of the satellite, the return goes to where the station WAS. The
station moves on in that frame while the light flies (2R/c), so the centre of
the return lobe misses the station by `2 v_perp R / c`. Seen from the station
that is the angle `theta = 2 v_perp / c`, the point-ahead angle of the
geometry (Degnan, DOI 10.1029/GD025p0133; see `../point_ahead_geometry/`). A
passive cube cannot lead the station the way an active terminal leads the
satellite, so the SAME angle that the uplink leads by is the angle the return
misses by.

The station still receives light, because the return is a diffraction lobe,
not a ray. For an unobscured cube of aperture D the lobe is the Airy pattern,
about `lambda / D` wide, and the light that reaches the station is the part of
that lobe that points at it. So the Term is the Airy fraction at the offset
(Born and Wolf, DOI 10.1017/CBO9781139644181, Sec. 8.5.2):

    loss_db = -10 log10( [2 J1(x) / x]^2 ),   x = pi D theta / lambda

`x` is the offset in units of the lobe width. The Term is the OFFSET only:
the down-leg geometric Term (with the top-hat correction) already gives the
on-axis spread of the lobe. The offset changes the POWER only. The light that
reaches the station comes from the cube (a point source at that range)
through the SAME air as the up leg, so the offset adds NO angle to the
turbulence path (see `../retro_bracket/README.md`).

## The run line

From the repository root:

```
python -m validation.retro_velocity_aberration.airy_aberration
```

## The numbers (1550 nm)

The loss against the offset:

| x | where the station sits | loss (dB) |
| --- | --- | --- |
| 0 | the lobe centre | 0.00 |
| 1.0 | about half-way out | 1.11 |
| 2.405 | the best cube size (below) | 7.30 |
| 3.832 | the first null | infinite (the script reads 118 dB of rounding) |
| 5.136 | the first bright ring | 17.57 |

Typical cubes:

| D | lambda / D (urad) | x at 50 urad | loss at 50 urad (dB) | x at 500 km, 30 deg | loss (dB) |
| --- | --- | --- | --- | --- | --- |
| 1 cm | 155.0 | 1.01 | 1.14 | 0.51 | 0.29 |
| 2 cm | 77.5 | 2.03 | 4.92 | 1.03 | 1.18 |
| 5 cm | 31.0 | 5.07 | 17.59 | 2.57 | 8.61 |
| 10 cm | 15.5 | 10.13 | 54.24 | 5.15 | 17.57 |

The 500 km, 30 deg angle is the `CircularOrbit` value, 25.4 urad (it reads
low, backlog 0-P18).

One day of real ISS passes (the TLE of `../point_ahead_geometry/`), elevation
above 10 deg:

| D | loss range (dB) | median (dB) |
| --- | --- | --- |
| 1 cm | 0.18 to 1.12 | 0.58 |
| 2 cm | 0.73 to 4.84 | 2.42 |
| 5 cm | 5.00 to 27.94 | 17.73 |

## The best cube size

A larger cube catches more of the up leg (D^2) and makes a brighter, narrower
lobe (D^2 on axis), but it aims that narrower lobe further off the station.
The return to the station scales as

    D^4 [2 J1(x) / x]^2  ~  [x J1(x)]^2

and `d(x J1(x))/dx = x J0(x)`, so the maximum is at the first zero of J0,
`x = 2.405` (the script finds 2.4048 numerically). So

    D_opt = 2.405 lambda / (pi theta) = 0.765 lambda / theta

which is 5.9 cm at 20 urad, 4.0 cm at 30 urad, and 2.4 cm at 50 urad (1550
nm). Past D_opt a larger cube returns LESS to the station.

## The spoiled cube (background, NOT BUILT)

A real retro array for a fast satellite often SPOILS the dihedral angles: each
angle between two faces is 90 deg plus a small offset delta (of the order of
1 to 3 arcsec). The return then splits into SIX beams on a ring about the
incoming ray, one for each order in which the light meets the three faces,
and the designer sets the ring radius to the aberration angle. Each beam comes
from about one sixth of the aperture, so each lobe is wider and fainter, but
the station sits near a lobe peak and not near a null. So a spoiled cube lets
a cube larger than D_opt work. The aberration changes size and turns in
direction along a pass, so the ring covers it only in part.

`olb.terminal.SpoiledCornerCube(dihedral_offset_rad)` is a STUB:
`retro_space_budget` raises `NotImplementedError` for it. The model needs the
full far-field pattern of a spoiled cube, not the Airy pattern. A source to
check before any code: R. F. Chang, D. G. Currie, C. O. Alley and M. E.
Pittman, J. Opt. Soc. Am. 61, 431 (1971), DOI 10.1364/JOSA.61.000431
(UNVERIFIED, like the ring-radius rule of about 3.3 n delta that this session
quoted from memory).

## Limits of the Term

- A circular, unobscured aperture. A real cube has a hexagonal pupil, and an
  uncoated (total internal reflection) cube splits the lobe by polarisation,
  so the real pattern is less smooth than the Airy disc, most of all near and
  past the first null. The Term flags `PAST THE AIRY NULL` at x >= 3.8317.
- A flat incident wavefront: the cube must be much smaller than the uplink
  coherence width at the satellite (the uplink speckle is metres wide there).
- A monostatic station: the offset is measured from the transmitter.
- The two DOIs above are NOT yet checked against the papers (open in
  CLAUDE.md).
