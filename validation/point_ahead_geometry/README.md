# The point-ahead angle: where it comes from and where olb uses it

Date: 2026-09-30. Branch `retro-velocity-aberration`. Machine: the laptop.

## The purpose

The point-ahead angle is `2 v_perp / c`: `v_perp` is the velocity of the
satellite relative to the station across the line of sight, in the inertial
frame, and the angle is how far the line of sight turns during one round trip
of the light (Degnan, DOI 10.1029/GD025p0133). This study records three
things:

1. The TLE calculator `olb.geometry.TLEPass.point_ahead_rad` (new on
   2026-09-30) and its checks.
2. The error of the old analytic `CircularOrbit` form (backlog 0-P18, DONE
   2026-09-30).
3. The map of every place olb reads the angle, and every place it must NOT.

## The run line

From the repository root:

```
python -m validation.point_ahead_geometry.tle_point_ahead
```

It needs `skyfield` only. The TLEs are in the script: the ISS TLE of the
skyfield documentation (2014-01-20, about 420 km) and a synthetic
geostationary satellite. The station is 51.5 N, 0.1 W for the ISS.

## 1. The TLE calculator

`TLEPass` takes `v_perp` from the skyfield GCRS relative velocity of the
satellite and the station. The station velocity holds the rotation of the
Earth, so a geostationary satellite has a point-ahead angle.

| Check | Result |
| --- | --- |
| THE DEFINITION: the inertial line-of-sight turn between t - R/c and t + R/c, 467 ISS samples | worst relative difference 4.3e-06 |
| GEO from the sub-point, against `2 (v_geo - omega_E R_E) / c` = 17.41 urad | 17.41 urad |
| GEO from 40 deg latitude (the quoted value is about 18 to 19 urad) | 18.13 urad (43.8 deg elevation), 18.49 urad (34.4 deg) |
| GEO from 60 deg latitude | 18.96 urad (22.0 deg elevation) |
| ISS, one day, above 5 deg (the quoted LEO value is about 50 urad, 10 arcsec) | 18.1 to 49.6 urad (3.7 to 10.2 arcsec) |

CAUTION FOR THE DEFINITION CHECK. A float64 Julian date resolves only about
40 us, and R/c is about 3 ms. A shift of the plain date reads the turn 1.3
percent wrong. The script shifts the split whole + fraction date of skyfield.

The `olb.geometry` self-check repeats the GEO closed form and the definition
check.

## 2. The CircularOrbit form (backlog 0-P18)

UP TO 2026-09-30 `CircularOrbit.point_ahead_rad` was `2 v_orb sin(el) / c`
(the column "old flat" below). Now it IS the overhead `cos(eta)` form, and
the script asserts that. The old form was a FLAT-EARTH
form of an overhead pass. For an overhead pass the velocity is horizontal at
the SATELLITE, so the exact value is `v_orb cos(eta)`, with the nadir angle
`sin(eta) = R_E cos(el) / (R_E + h)`. `sin(el)` equals `cos(eta)` only at the
zenith or when `h -> 0`.

h = 420 km, `v_orb` = 7.66 km/s, in urad:

| el (deg) | ISS TLE (min to max) | overhead `cos(eta)` | old flat `sin(el)` |
| --- | --- | --- | --- |
| 85 | 49.2 | 50.9 | 50.9 |
| 60 | 48.7 to 49.2 | 45.1 | 44.3 |
| 45 | 37.0 to 44.3 | 38.2 | 36.1 |
| 30 | 28.4 to 34.7 | 29.8 | 25.6 |
| 20 | 23.2 to 27.5 | 24.1 | 17.5 |
| 10 | 19.0 to 42.7 | 19.6 | 8.9 |

THE READING:

1. `sin(el)` reads LOW away from the zenith: 14 percent at 30 deg, 27
   percent at 20 deg, 55 percent at 10 deg against the overhead form.
2. ELEVATION IS NOT ENOUGH. At one elevation an off-track pass moves more
   across the line of sight than an overhead pass, so a real pass spreads
   (19 to 43 urad at 10 deg). No elevation-only form holds that; a `TLEPass`
   does.
3. The overhead `cos(eta)` form is NOT a strict lower bound. The rotation of
   the Earth adds or removes up to about 0.4 km/s, so a real pass can read
   below it (28.4 against 29.8 urad at 30 deg).

The fix is DONE (2026-09-30, option B of backlog 0-P18, branch
`spherical-point-ahead`): `CircularOrbit` gives the overhead `cos(eta)` form,
and the slew rate reads `v_perp / (L sin(el))`. It moved every reader of the
angle (the list below), and a stored point-ahead campaign with
`"geometry"` angles now makes a new key.

## 3. Where olb uses the angle

ONE SOURCE: `geometry.point_ahead_rad` (`CircularOrbit` or `TLEPass`). Each
reader also takes an explicit value.

| Link | Rung | Use | Override |
| --- | --- | --- | --- |
| Pre-compensated uplink (`DownlinkBeacon`) | 0 | `uplink_point_ahead_term`: the Stone anisoplanatic phase variance at theta | the geometry only |
| Pre-compensated uplink | 1 | `uplink_fast_term`: the FAST `DTHETA` in arcsec | `fast_params={"DTHETA": ...}` |
| Pre-compensated uplink | 2 | the runner `point_ahead_rad="auto"` resolves to `"geometry"`: one more pass for each trial on a window shifted by theta * z at each screen (`eta_turb_pa`); the budget selects the stored angle of its geometry | the runner `point_ahead_rad=` (a float or a list) |
| Retro | 0, 1, 2 | `retro_velocity_aberration_term`: the Airy loss of the return lobe that the cube sends theta off the station (see `../retro_velocity_aberration/`) | `retro_space_budget(aberration_rad=...)`; 0.0 turns it off |

WHERE THE ANGLE MUST NOT GO:

- An UNCORRECTED uplink and a downlink: nothing decorrelates, so the runner
  `"auto"` resolves to None.
- THE RETRO TURBULENCE. The up leg and the return of ONE pulse share ONE
  line: the station leads the satellite by theta, and the line of sight
  turns by the same theta while the light flies, so the two cancel. Only the
  WIND moves the air between the two passes (over 2R/c). So the retro runner
  REFUSES an explicit `point_ahead_rad`, and it shifts the up-leg window by
  the wind drift (see `../retro_bracket/`). The point-ahead geometry holds
  only for two beams at the SAME moment (the beacon that arrives now and the
  uplink that leaves now).
- A terrestrial path: both terminals stand still, and the runner refuses the
  angle.

`point_ahead_rad` in the runner, `Campaign` and `run_fidelity2` means the
shifted-window pass of the pre-compensated uplink. The retro budget argument
is `aberration_rad`, a different quantity from the same source.
