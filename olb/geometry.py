'''
Link geometry, with two backends that you can exchange.

The models and the budget read only two arrays from a geometry object:

    geom.elevation_deg    # elevation above the horizon [deg]
    geom.slant_range_m    # ground-station -> satellite range [m]

The source of the two arrays does not change the models or the budget. The
source is an analytic circular orbit or a real TLE that skyfield propagates.
Select the backend for the task:

    CircularOrbit  -- analytic, vectorised over an elevation grid. Use it for
                      parameter sweeps and Monte Carlo. It gives a regular
                      elevation axis and high speed.
    TLEPass        -- a real pass of a real satellite, sampled in time with
                      skyfield. Use it to replay an actual Kepler pass.

Each backend also gives the extra quantities that it can compute at low cost
(point-ahead and slew for the analytic orbit; point-ahead, azimuth and times
for the TLE pass). Select the backend that gives the quantities that your code needs.
'''

import numpy as np

# Physical constants for the circular-orbit geometry.
_GRAV_CONST = 6.674e-11    # gravitational constant [m^3 kg^-1 s^-2]
_EARTH_MASS = 5.972e24     # mass of the Earth [kg]
_EARTH_RADIUS = 6371e3     # mean radius of the Earth [m]
_C = 2.998e8               # speed of light [m/s]


class Satellite:
    '''A satellite in a circular orbit at a given altitude.'''

    def __init__(self, altitude):
        '''
        Parameters:
            altitude : float
                The orbital altitude above the Earth surface [m].
        '''
        self.altitude = altitude
        self.orbital_speed = np.sqrt(
            _GRAV_CONST * _EARTH_MASS / (_EARTH_RADIUS + self.altitude))

    def angular_speed(self):
        '''Return the angular speed about the Earth centre [deg/s].'''
        return np.rad2deg(self.orbital_speed / (_EARTH_RADIUS + self.altitude))


class SatellitePass:
    '''A pass of a Satellite over a ground station at a given elevation.'''

    def __init__(self, satellite, elevation):
        '''
        Parameters:
            satellite : Satellite
                The satellite that the ground station tracks.
            elevation : float
                The elevation angle above the horizon [deg].
        '''
        self.elevation = elevation
        self.satellite = satellite

    def tangential_velocity(self):
        '''
        Return the satellite velocity across the line of sight [m/s].

        This is the orbital-speed component transverse to the line of sight, as
        the ground station sees it, for an OVERHEAD pass on a SPHERICAL Earth
        that does not turn. The velocity is horizontal at the SATELLITE, so
        v_perp = v_orb cos(eta), with the nadir angle eta from the sine rule of
        the centre-station-satellite triangle,
        sin(eta) = R_E cos(el) / (R_E + h) (Degnan, DOI 10.1029/GD025p0133).
        The flat-Earth form v_orb sin(el) is the h -> 0 limit; at 30 deg it
        read 16 percent low at 500 km and 30 percent low at 1500 km (backlog
        0-P18). An off-track pass and the
        rotation of the Earth move a real pass off this value; use a TLEPass.
        '''
        Re = _EARTH_RADIUS
        h = self.satellite.altitude
        sin_eta = Re * np.cos(np.radians(self.elevation)) / (Re + h)
        return self.satellite.orbital_speed * np.sqrt(1.0 - sin_eta ** 2)

    def slant_range(self):
        '''Return the slant range from the ground station to the satellite [m].'''
        Re = _EARTH_RADIUS
        h = self.satellite.altitude
        el = np.radians(self.elevation)
        return -Re * np.sin(el) + np.sqrt(
            Re ** 2 * np.sin(el) ** 2 + h ** 2 + 2 * Re * h)

    def point_ahead_angle(self):
        '''
        Return the point-ahead angle [rad].

        This is the angular lead that accounts for the finite speed of light
        over the round trip.
        '''
        return 2 * self.tangential_velocity() / _C

    def apparent_slew_rate(self):
        '''Return the apparent angular slew rate of the line of sight [deg/s].

        THE RATE IS REFERENCED TO THE ALTITUDE, not to the slant range. It is
        the coefficient of the `ws*h` term of the Bufton wind profile with h
        the ALTITUDE of a layer (Andrews and Phillips, DOI 10.1117/3.626196,
        Ch. 12, Eqs. (2) and (3), printed p. 481). So the apparent speed at a
        layer is this rate times the ALTITUDE of the layer, and NOT times its
        slant distance. The frozen-flow time axis reads it that way (see
        olb.waveoptics.turbulence.temporal.strip_plan).

        The line of sight turns at omega = v_perp / L (L the slant range), and
        a layer at altitude h sits at the slant distance h / sin(el) (a thin,
        plane-parallel atmosphere). So the layer speed is
        omega h / sin(el) = [v_perp / (L sin(el))] h, and this rate is
        v_perp / (L sin(el)). On a flat Earth L sin(el) = h_sat; on the sphere
        L sin(el) < h_sat (backlog 0-P18).
        '''
        return np.rad2deg(self.tangential_velocity()
                          / (self.slant_range() * np.sin(np.radians(self.elevation))))


class CircularOrbit:
    '''Analytic circular-orbit geometry over an elevation grid.'''

    def __init__(self, altitude_m, elevation_deg):
        '''
        Parameters:
            altitude_m : float
                Orbital altitude above the Earth's surface [m].
            elevation_deg : float or array
                Elevation angle(s) above the horizon [deg]. Use an array to
                sweep.
        '''
        self.altitude_m = altitude_m
        self.elevation_deg = np.asarray(elevation_deg, dtype=float)
        self._pass = SatellitePass(Satellite(altitude_m), self.elevation_deg)

    @property
    def slant_range_m(self):
        return self._pass.slant_range()

    @property
    def point_ahead_rad(self):
        '''Point-ahead angle [rad] from the finite speed of light.'''
        return self._pass.point_ahead_angle()

    @property
    def slew_deg_s(self):
        '''Apparent line-of-sight slew rate [deg/s].'''
        return self._pass.apparent_slew_rate()


class HorizontalPath:
    '''
    Horizontal (terrestrial) path geometry: a constant range, no elevation.

    A terrestrial link runs ground-to-ground along a horizontal path. The range
    is the path length and it does not change with any elevation angle. So this
    geometry exposes only slant_range_m; it has no elevation_deg. The
    range-only models (geometric spreading, pointing) read slant_range_m and
    work unchanged. The horizontal extinction and scintillation terms read the
    path length and the constant Cn2 from the TerrestrialChannel, not from an
    elevation, so they do not need an elevation here.
    '''

    def __init__(self, path_length_m):
        '''
        Parameters:
            path_length_m : float or array
                Horizontal path length L [m]. Use an array to sweep the range.
        '''
        self.path_length_m = np.asarray(path_length_m, dtype=float)

    @property
    def slant_range_m(self):
        return self.path_length_m


class TLEPass:
    '''A real satellite pass from a TLE, propagated with skyfield.'''

    def __init__(self, tle_line1, tle_line2, lat_deg, lon_deg, alt_m, times,
                 name=""):
        '''
        Parameters:
            tle_line1, tle_line2 : str
                The two-line element set.
            lat_deg, lon_deg : float
                Ground station geodetic latitude / longitude [deg].
            alt_m : float
                Ground station height above the WGS84 ellipsoid [m].
            times : skyfield Time
                Array of times to sample the pass at (see from_window).
            name : str
                Satellite name (cosmetic).

        After construction, elevation_deg / azimuth_deg / slant_range_m /
        point_ahead_rad are arrays over `times`. Elevation is negative when the
        satellite is below the horizon. Use the mask elevation_deg > 0 for the
        visible pass.
        '''
        from skyfield.api import load, wgs84, EarthSatellite
        ts = load.timescale()
        satellite = EarthSatellite(tle_line1, tle_line2, name, ts)
        observer = wgs84.latlon(lat_deg, lon_deg, alt_m)
        topocentric = (satellite - observer).at(times)
        alt, az, dist = topocentric.altaz()

        self.name = name
        self.times = times
        self.elevation_deg = alt.degrees
        self.azimuth_deg = az.degrees
        self.slant_range_m = dist.m
        # THE POINT-AHEAD ANGLE is 2 v_perp / c, with v_perp the velocity of the
        # satellite RELATIVE TO THE STATION across the line of sight, in the
        # INERTIAL (GCRS) frame: the station velocity holds the rotation of
        # the Earth, so a geostationary satellite still has an angle
        # (J. J. Degnan, Geodynamics Series 25, 133 (1993),
        # DOI 10.1029/GD025p0133). skyfield gives the relative position and
        # velocity of the topocentric vector in GCRS axes.
        r = topocentric.position.m
        v = topocentric.velocity.m_per_s
        r_hat = r / np.linalg.norm(r, axis=0)
        v_perp = v - np.sum(v * r_hat, axis=0) * r_hat
        self.point_ahead_rad = 2.0 * np.linalg.norm(v_perp, axis=0) / _C

    @classmethod
    def from_window(cls, tle_line1, tle_line2, lat_deg, lon_deg, alt_m,
                    start_utc, duration_s, step_s=1.0, name=""):
        '''
        Build a pass. Sample a time window at a fixed step.

        Parameters:
            start_utc : tuple
                (year, month, day, hour, minute, second) UTC start.
            duration_s : float
                Window length [s].
            step_s : float
                Sample step [s].
            (remaining parameters: see __init__)
        '''
        from skyfield.api import load
        ts = load.timescale()
        year, month, day, hour, minute, second = start_utc
        seconds = second + np.arange(0.0, duration_s, step_s)
        times = ts.utc(year, month, day, hour, minute, seconds)
        return cls(tle_line1, tle_line2, lat_deg, lon_deg, alt_m, times, name)


if __name__ == '__main__':
    # ---- the CircularOrbit overhead pass, against a direct in-plane orbit ----
    # Station at (0, R_E); satellite at the central angle g on a circle of
    # radius R_E + h. The line-of-sight angle, differenced in time, is omega;
    # v_perp = omega * L, and the rate referenced to the altitude is
    # omega / sin(el) (see apparent_slew_rate).
    for h in (420e3, 1500e3):
        sat = Satellite(h)
        n = sat.orbital_speed / (_EARTH_RADIUS + h)
        g = np.linspace(0.0, 0.3, 3001)
        dt = 1e-3

        def los(t):
            a = g + n * t
            return np.arctan2((_EARTH_RADIUS + h) * np.cos(a) - _EARTH_RADIUS,
                              (_EARTH_RADIUS + h) * np.sin(a))
        el_deg = np.rad2deg(los(0.0))
        omega = np.abs(los(dt) - los(-dt)) / (2 * dt)
        keep = el_deg > 10.0
        orb = CircularOrbit(h, el_deg[keep])
        assert np.allclose(orb.point_ahead_rad,
                           2 * omega[keep] * orb.slant_range_m / _C, rtol=1e-5)
        assert np.allclose(np.deg2rad(orb.slew_deg_s),
                           omega[keep] / np.sin(np.radians(el_deg[keep])), rtol=1e-5)
    # The zenith keeps the old flat value; 420 km, 30 deg reads the 29.8 urad
    # of validation/point_ahead_geometry (the flat form gave 25.6).
    assert np.isclose(CircularOrbit(420e3, 90.0).point_ahead_rad,
                      2 * Satellite(420e3).orbital_speed / _C)
    assert abs(CircularOrbit(420e3, 30.0).point_ahead_rad * 1e6 - 29.8) < 0.1
    print(f"CircularOrbit 420 km: {CircularOrbit(420e3, 30.0).point_ahead_rad * 1e6:.2f}"
          f" urad at 30 deg (overhead sphere)")

    # ---- the TLE point-ahead angle, against two independent routes ----
    from skyfield.api import load, wgs84, EarthSatellite
    ts = load.timescale()

    def _tle(line):
        '''Put the mod-10 checksum on column 69 of a TLE line.'''
        s = sum(int(c) if c.isdigit() else int(c == "-") for c in line[:68])
        return line[:68] + str(s % 10)

    # A synthetic GEO (i = 0, e = 0). Seen from its sub-point the relative
    # inertial velocity is v_geo - omega_E * R_E, all of it across the line of
    # sight, so the angle is 2 (v_geo - omega_E R_E) / c.
    GEO = (_tle("1 99999U 23001A   23001.50000000  .00000000  00000-0  00000+0 0  9990"),
           _tle("2 99999   0.0000   0.0000 0000001   0.0000 100.0000  1.00273791    10"))
    sub = wgs84.subpoint(EarthSatellite(*GEO, "GEO", ts).at(ts.utc(2023, 1, 1, 12)))
    geo = TLEPass.from_window(*GEO, lat_deg=0.0, lon_deg=sub.longitude.degrees,
                              alt_m=0.0, start_utc=(2023, 1, 1, 12, 0, 0),
                              duration_s=1.0)
    omega_e = 7.2921159e-5                       # sidereal rate [rad/s]
    a_geo = (_GRAV_CONST * _EARTH_MASS / omega_e ** 2) ** (1 / 3)
    want = 2 * (omega_e * a_geo - omega_e * 6378137.0) / _C
    assert abs(geo.point_ahead_rad[0] / want - 1) < 5e-3, (geo.point_ahead_rad, want)

    # ISS (the skyfield documentation TLE). The lead IS the turn of the
    # INERTIAL line of sight between t - R/c and t + R/c.
    ISS = (_tle("1 25544U 98067A   14020.93268519  .00009878  00000-0  18200-3 0  5082"),
           _tle("2 25544  51.6498 109.4756 0003572  55.9686 274.8005 15.49815350868473"))
    iss = TLEPass.from_window(*ISS, lat_deg=51.5, lon_deg=-0.1, alt_m=50.0,
                              start_utc=(2014, 1, 21, 0, 0, 0),
                              duration_s=86400, step_s=30.0)
    vis = iss.elevation_deg > 5.0
    t, rng = iss.times[vis], iss.slant_range_m[vis]
    topo = EarthSatellite(*ISS, "ISS", ts) - wgs84.latlon(51.5, -0.1, 50.0)
    dt = rng / _C / 86400.0
    r0 = topo.at(ts.tt_jd(t.whole, t.tt_fraction - dt)).position.m
    r1 = topo.at(ts.tt_jd(t.whole, t.tt_fraction + dt)).position.m
    turn = np.arccos(np.sum(r0 * r1, axis=0) / (
        np.linalg.norm(r0, axis=0) * np.linalg.norm(r1, axis=0)))
    assert np.allclose(turn, iss.point_ahead_rad[vis], rtol=1e-4)
    assert 40e-6 < iss.point_ahead_rad[vis].max() < 55e-6    # ~10 arcsec
    print(f"GEO point ahead {geo.point_ahead_rad[0] * 1e6:.2f} urad "
          f"(closed form {want * 1e6:.2f}); ISS {vis.sum()} samples, "
          f"{iss.point_ahead_rad[vis].min() * 1e6:.1f} to "
          f"{iss.point_ahead_rad[vis].max() * 1e6:.1f} urad")
    print("self-check passed")
