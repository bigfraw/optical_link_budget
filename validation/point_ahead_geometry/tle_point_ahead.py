"""Check the TLE point-ahead angle, and measure the CircularOrbit form against it.

The point-ahead angle is 2 v_perp / c, with v_perp the velocity of the
satellite relative to the station across the line of sight, in the inertial
frame (Degnan, DOI 10.1029/GD025p0133). `olb.geometry.TLEPass.point_ahead_rad`
takes v_perp from the skyfield GCRS velocity. This script checks it three ways
and compares the analytic `CircularOrbit` form `2 v_orb sin(el) / c`:

1. THE DEFINITION. The lead is the turn of the INERTIAL line of sight between
   t - R/c and t + R/c. The script propagates the satellite to those two times
   (with the split whole + fraction Julian date: a float64 date resolves only
   about 40 us, and R/c is about 3 ms).
2. GEO. From the sub-point the relative inertial velocity is
   v_geo - omega_E R_E, all of it across the line of sight.
3. THE OVERHEAD CLOSED FORM. For an overhead pass (Earth rotation off),
   v_perp = v_orb cos(eta), with the nadir angle
   sin(eta) = R_E cos(el) / (R_E + h).

Run from the repository root:

    python -m validation.point_ahead_geometry.tle_point_ahead
"""

import numpy as np
from skyfield.api import load, wgs84, EarthSatellite

from olb.geometry import (CircularOrbit, TLEPass, _C, _EARTH_MASS,
                          _EARTH_RADIUS, _GRAV_CONST)

ARCSEC = np.pi / 180 / 3600
ts = load.timescale()


def tle(line):
    """Put the mod-10 checksum on column 69 of a TLE line."""
    s = sum(int(c) if c.isdigit() else int(c == "-") for c in line[:68])
    return line[:68] + str(s % 10)


# The ISS TLE of the skyfield documentation (2014-01-20), about 420 km.
ISS = (tle("1 25544U 98067A   14020.93268519  .00009878  00000-0  18200-3 0  5082"),
       tle("2 25544  51.6498 109.4756 0003572  55.9686 274.8005 15.49815350868473"))
# A synthetic GEO: i = 0, e = 0, one sidereal day.
GEO = (tle("1 99999U 23001A   23001.50000000  .00000000  00000-0  00000+0 0  9990"),
       tle("2 99999   0.0000   0.0000 0000001   0.0000 100.0000  1.00273791    10"))
STATION = dict(lat_deg=51.5, lon_deg=-0.1, alt_m=50.0)
H_ISS = 420e3


def main():
    # ---- 1. the ISS over one day ----
    iss = TLEPass.from_window(*ISS, **STATION, start_utc=(2014, 1, 21, 0, 0, 0),
                              duration_s=86400, step_s=5.0, name="ISS")
    vis = iss.elevation_deg > 5.0
    el, pa = iss.elevation_deg[vis], iss.point_ahead_rad[vis]
    t, rng = iss.times[vis], iss.slant_range_m[vis]
    topo = (EarthSatellite(*ISS, "ISS", ts)
            - wgs84.latlon(STATION["lat_deg"], STATION["lon_deg"],
                           STATION["alt_m"]))
    dt = rng / _C / 86400.0
    r0 = topo.at(ts.tt_jd(t.whole, t.tt_fraction - dt)).position.m
    r1 = topo.at(ts.tt_jd(t.whole, t.tt_fraction + dt)).position.m
    turn = np.arccos(np.clip(np.sum(r0 * r1, axis=0) / (
        np.linalg.norm(r0, axis=0) * np.linalg.norm(r1, axis=0)), -1, 1))
    print(f"ISS, one day, {vis.sum()} samples above 5 deg "
          f"(station {STATION['lat_deg']} N, {STATION['lon_deg']} E)")
    print(f"  1. velocity route against the line-of-sight turn over 2R/c: "
          f"worst relative difference {np.max(np.abs(turn / pa - 1)):.1e}")
    print(f"  range over the visible passes: {pa.min() * 1e6:.1f} to "
          f"{pa.max() * 1e6:.1f} urad ({pa.min() / ARCSEC:.1f} to "
          f"{pa.max() / ARCSEC:.1f} arcsec)")

    # ---- 3. the CircularOrbit form and the overhead closed form ----
    v_orb = np.sqrt(_GRAV_CONST * _EARTH_MASS / (_EARTH_RADIUS + H_ISS))
    print(f"\n  3. by elevation (h = {H_ISS / 1e3:.0f} km, v_orb = "
          f"{v_orb / 1e3:.2f} km/s), in urad:")
    print(f"  {'el':>5} {'ISS TLE (min-max)':>18} {'overhead cos(eta)':>18} "
          f"{'CircularOrbit sin(el)':>22}")
    for e in (85, 60, 45, 30, 20, 10):
        sel = np.abs(el - e) < 1.0
        eta = np.arcsin(_EARTH_RADIUS * np.cos(np.radians(e))
                        / (_EARTH_RADIUS + H_ISS))
        closed = 2 * v_orb * np.cos(eta) / _C
        circ = float(np.ravel(CircularOrbit(
            altitude_m=H_ISS, elevation_deg=[e]).point_ahead_rad)[0])
        tle_txt = (f"{pa[sel].min() * 1e6:6.1f} - {pa[sel].max() * 1e6:5.1f}"
                   if sel.any() else "(no sample)")
        print(f"  {e:5d} {tle_txt:>18} {closed * 1e6:18.1f} {circ * 1e6:22.1f}")

    # ---- 2. GEO ----
    sub = wgs84.subpoint(EarthSatellite(*GEO, "GEO", ts).at(
        ts.utc(2023, 1, 1, 12)))
    omega_e = 7.2921159e-5
    a_geo = (_GRAV_CONST * _EARTH_MASS / omega_e ** 2) ** (1 / 3)
    closed = 2 * (omega_e * a_geo - omega_e * 6378137.0) / _C
    print(f"\n  2. GEO (closed form from the sub-point "
          f"2 (v_geo - omega_E R_E) / c = {closed * 1e6:.2f} urad):")
    for lat, dlon in ((0.0, 0.0), (40.0, 0.0), (40.0, 30.0), (60.0, 0.0)):
        g = TLEPass.from_window(*GEO, lat_deg=lat,
                                lon_deg=sub.longitude.degrees + dlon, alt_m=0.0,
                                start_utc=(2023, 1, 1, 12, 0, 0), duration_s=1.0)
        print(f"  station lat {lat:4.1f}, lon offset {dlon:4.1f}: elevation "
              f"{g.elevation_deg[0]:5.1f} deg, point ahead "
              f"{g.point_ahead_rad[0] * 1e6:6.2f} urad")


if __name__ == "__main__":
    main()
