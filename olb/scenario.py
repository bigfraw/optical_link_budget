'''
Scenario: the pure-data description of a link case.

This module does not compute values. These dataclasses hold the inputs that the
models read. A "case" is a value that you build, copy, change, and sweep. The
models read a Scenario and a geometry and return Terms. The Scenario does not
import the models. The data moves in one direction, from the inputs to the
models.
'''

from dataclasses import dataclass, field
from typing import Literal, Optional

Direction = Literal["uplink", "downlink"]


@dataclass
class Link:
    '''Optical link hardware and beam parameters.'''
    wavelength_m: float = 1550e-9
    direction: Direction = "uplink"
    tx_waist_m: float = 0.1              # transmit 1/e^2 beam waist w0 [m]
    tx_power_dbm: Optional[float] = None  # launch power [dBm]; None if only losses matter
    m2: float = 1.0                      # beam quality M^2 (>= 1)
    rx_diameter_m: float = 0.08          # receive aperture diameter [m]
    rx_obscuration_ratio: float = 0.0    # central obscuration / aperture diameter
    rx_sensitivity_dbm: Optional[float] = None  # required received power [dBm]
    pointing_jitter_rad: float = 0.0     # 1-sigma tracking jitter [rad]
    divergence_rad: Optional[float] = None  # transmit far-field 1/e^2 HALF-angle divergence [rad]; None = collimated (diffraction limit)
    retro_aperture_m: Optional[float] = None  # satellite retroreflector aperture diameter [m]; used by retro_budget


@dataclass
class Site:
    '''Ground station location and atmosphere.'''
    lat_deg: float = -29.0468           # TN-2 default (Kepler OGS)
    lon_deg: float = 115.3467
    alt_m: float = 269.0
    cn2_ground: float = 1.7e-14         # Hufnagel-Valley ground-level Cn2 scale (HV57 A)
    wind_rms_m_s: float = 21.0          # Bufton wind profile rms [m/s]
    clear_sky_probability: float = 1.0  # cloud-free-line-of-sight fraction (0-1)


@dataclass
class Scenario:
    '''A full link case: hardware + site + orbit.'''
    link: Link = field(default_factory=Link)
    site: Site = field(default_factory=Site)
    altitude_m: float = 600e3           # satellite altitude, for analytic geometry
    availability_target: float = 0.99   # target link availability (0-1)
