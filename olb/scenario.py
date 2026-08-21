'''
Scenario: the pure-data description of a link case.

This module does not compute values. These dataclasses hold the inputs that the
models read. A "case" is a value that you build, copy, change, and sweep. The
models read a Scenario and a geometry and return Terms. The Scenario does not
import the models. The data moves in one direction, from the inputs to the
models.

Data model: all terminal hardware lives on a Terminal (see olb.terminal). A
Scenario holds two terminals, `ground` and `space`, a `Channel` (the propagation
channel), and the link `direction`. The channel holds no hardware.

The models do not read `ground` or `space` directly. They read the resolved
roles `tx_terminal` and `rx_terminal`, which the `direction` sets:

    direction   tx_terminal   rx_terminal
    uplink      ground        space
    downlink    space         ground
    retro       ground        ground

So one Terminal serves both link directions, and a terminal parameter can only
be set through a Terminal.
'''

from dataclasses import dataclass, field
from typing import Literal

from .terminal import Terminal

Direction = Literal["uplink", "downlink", "retro"]


@dataclass
class Site:
    '''Ground station location and atmosphere (the propagation medium).'''
    lat_deg: float = -29.0468           # TN-2 default (Kepler OGS)
    lon_deg: float = 115.3467
    alt_m: float = 269.0
    cn2_ground: float = 1.7e-14         # Hufnagel-Valley ground-level Cn2 scale (HV57 A)
    wind_rms_m_s: float = 21.0          # Bufton wind profile rms [m/s]
    clear_sky_probability: float = 1.0  # cloud-free-line-of-sight fraction (0-1)


@dataclass
class Channel:
    '''
    The propagation channel: the ground site plus the satellite orbit.

    The channel holds no terminal hardware. It is the intended seam for a later
    terrestrial (horizontal-path) channel. Do not build that variant now.

    Parameters:
        site : Site
            The ground station location and atmosphere (the medium).
        altitude_m : float
            The satellite altitude [m], for the analytic orbit geometry.
    '''
    site: Site = field(default_factory=Site)
    altitude_m: float = 600e3


@dataclass
class Scenario:
    '''A full link case: two terminals + a Channel + direction.'''
    ground: Terminal
    space: Terminal
    direction: Direction = "uplink"
    channel: Channel = field(default_factory=Channel)
    availability_target: float = 0.99   # target link availability (0-1)

    @property
    def tx_terminal(self) -> Terminal:
        '''The transmit terminal for the link direction (see the module docstring).'''
        return self.space if self.direction == "downlink" else self.ground

    @property
    def rx_terminal(self) -> Terminal:
        '''The receive terminal for the link direction (see the module docstring).'''
        return self.space if self.direction == "uplink" else self.ground


if __name__ == '__main__':
    from .terminal import Transmitter, Aperture

    ground = Terminal(aperture_m=0.7, transmitter=Transmitter(waist_m=0.1))
    space = Terminal(aperture_m=0.05, detector=Aperture())

    up = Scenario(ground=ground, space=space, direction="uplink")
    assert up.tx_terminal is ground and up.rx_terminal is space

    down = Scenario(ground=ground, space=space, direction="downlink")
    assert down.tx_terminal is space and down.rx_terminal is ground

    retro = Scenario(ground=ground, space=space, direction="retro")
    assert retro.tx_terminal is ground and retro.rx_terminal is ground

    assert up.channel.altitude_m == 600e3 and up.availability_target == 0.99
    assert up.channel.site.cn2_ground == 1.7e-14
    print("Scenario:", up)
    print("self-check passed")
