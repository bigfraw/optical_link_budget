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
(point-ahead and slew for the analytic orbit; azimuth and times for the TLE
pass). Select the backend that gives the quantities that your code needs.
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
        the ground station sees it.
        '''
        return self.satellite.orbital_speed * np.sin(np.radians(self.elevation))

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
        '''Return the apparent angular slew rate of the line of sight [deg/s].'''
        return np.rad2deg(self.tangential_velocity() / self.satellite.altitude)


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

        After construction, elevation_deg / azimuth_deg / slant_range_m are
        arrays over `times`. Elevation is negative when the satellite is below
        the horizon. Use the mask elevation_deg > 0 for the visible pass.
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
