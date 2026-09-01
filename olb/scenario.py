'''
Scenario families: the pure-data description of a link case.

This module does not compute values. These dataclasses hold the inputs that the
models read. A "case" is a value that you build, copy, change, and sweep. The
models read a scenario and a geometry and return Terms. The scenario does not
import the models. The data moves in one direction, from the inputs to the
models.

Two scenario families, one contract. A link is either a SPACE link (a ground
station and a satellite) or a TERRESTRIAL link (two ground stations on a
horizontal path). Each family names its two terminals for what they physically
are, so the field names never lie:

    SpaceScenario        ground, space   + Channel (site + orbit altitude)
    TerrestrialScenario  near,   far     + TerrestrialChannel (site + path)

Both families expose the SAME thin interface that the models read:

    scenario.tx_terminal   the transmit terminal
    scenario.rx_terminal   the receive terminal
    scenario.channel       the propagation channel

So no model changes between the two families. A SpaceScenario resolves the two
roles from its `direction`:

    direction   tx_terminal   rx_terminal
    uplink      ground        space
    downlink    space         ground
    retro       ground        ground

A TerrestrialScenario resolves its two roles from its own `direction`. A
horizontal path is reciprocal, so the channel is the same in the two
directions and only the role mapping changes:

    direction   tx_terminal   rx_terminal
    forward     near          far
    reverse     far           near

The terrestrial `direction` is a DIFFERENT type from the space `Direction`:
"terrestrial" is a channel family, not a tx/rx geometry, so the two families
do not share their direction names. All terminal hardware lives on a Terminal
(see olb.terminal); a channel holds no hardware.

A SpaceScenario also carries an optional pre-compensation source for the uplink
(see `precompensation`). It names what the ground terminal senses to build the
uplink correction: the downlink beam (DownlinkBeacon), a laser guide star
(LaserGuideStar, a placeholder), or nothing (None, an uncorrected uplink).
'''

from dataclasses import dataclass, field
from typing import Literal, Optional, Union

from .terminal import Terminal

Direction = Literal["uplink", "downlink", "retro"]
TerrestrialDirection = Literal["forward", "reverse"]


# --- Pre-compensation source (uplink only) ----------------------------------

@dataclass
class DownlinkBeacon:
    '''
    Uplink pre-compensation sensed from the satellite downlink beam.

    The ground terminal senses the turbulence on the downlink beam and applies
    the conjugate to the uplink beam. This is reciprocity-based pre-compensation:
    the up and down paths share the same turbulence, so the downlink phase gives
    the uplink correction. But the downlink arrives from where the satellite was,
    and the uplink must go to where the satellite will be. The two directions
    differ by the point-ahead angle. So the correction removes only the part of
    each mode that stays correlated across that angle. The modal DECORRELATION
    residual stays. See olb.links.uplink.uplink_point_ahead_term and
    olb.turbulence.anisoplanatism. The satellite (space) terminal needs a
    transmitter for the downlink beam.
    '''
    pass


@dataclass
class LaserGuideStar:
    '''
    Uplink pre-compensation sensed from a ground-launched laser guide star.

    NOT YET IMPLEMENTED. A laser guide star focuses at a finite altitude, so its
    light samples a cone of the turbulence, not the full column. This gives focal
    (cone) anisoplanatism. That is a different effect from the point-ahead angular
    anisoplanatism of the downlink beacon. This class is a placeholder for a later
    task.

    Parameters:
        altitude_m : float
            Height of the guide star above the ground [m]. The default is a
            sodium-layer guide star.
    '''
    altitude_m: float = 90e3


# A pre-compensation source is one of the reference sources above. None means no
# pre-compensation: the uplink is uncorrected.
PreCompensationSource = Union[DownlinkBeacon, LaserGuideStar]


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
    A space propagation channel: the ground site plus the satellite orbit.

    The channel holds no terminal hardware. The space links read `altitude_m`
    for the analytic orbit geometry and build a Cn2(h) profile from the site.

    Parameters:
        site : Site
            The ground station location and atmosphere (the medium).
        altitude_m : float
            The satellite altitude [m], for the analytic orbit geometry.
    '''
    site: Site = field(default_factory=Site)
    altitude_m: float = 600e3


@dataclass
class TerrestrialChannel:
    '''
    A horizontal (terrestrial) propagation channel: a ground-to-ground path.

    The terrestrial counterpart of Channel. It holds no terminal hardware. The
    path is horizontal, so there is no orbit altitude and no elevation angle.
    Turbulence along a horizontal path is roughly uniform, so a single scalar
    Cn2 describes it. The extinction is a plain Beer-Lambert loss over the path,
    quoted directly as a dB-per-km coefficient (weather- and visibility-
    dependent; the user sets it per site).

    Parameters:
        site : Site
            The ground atmosphere along the path (the medium).
        path_length_m : float
            Horizontal path length L [m] for the terrestrial link.
        attenuation_db_per_km : float
            Clear-air / haze extinction coefficient [dB/km]. The Beer-Lambert
            loss is attenuation_db_per_km * (L / 1000).
        cn2 : float
            Constant refractive-index structure parameter Cn2 [m^-2/3] along the
            path. A single value, because a horizontal path sees ~uniform
            turbulence. Read by the (pending) terrestrial scintillation term.
    '''
    site: Site = field(default_factory=Site)
    path_length_m: float = 1e3
    attenuation_db_per_km: float = 0.5
    cn2: float = 1e-14


@dataclass
class SpaceScenario:
    '''
    A space link case: a ground terminal + a space terminal + direction.

    The optional `precompensation` field names the source that drives the uplink
    pre-compensation. It applies to the UPLINK direction only. None means the
    uplink is uncorrected. A DownlinkBeacon senses the downlink beam, so it drives
    the point-ahead anisoplanatism (see olb.links.uplink). A LaserGuideStar is a
    placeholder for a later task. A downlink or a retro scenario refuses the
    field at construction, because no model reads it there and a silent ignore
    hides a user error.
    '''
    ground: Terminal
    space: Terminal
    direction: Direction = "uplink"
    channel: Channel = field(default_factory=Channel)
    availability_target: float = 0.99
    precompensation: Optional[PreCompensationSource] = None   # uplink only

    def __post_init__(self):
        # Refuse a pre-compensation source on a non-uplink link. No model reads
        # the field there, so a silent ignore hides a user error (backlog I-4).
        if self.precompensation is not None and self.direction != "uplink":
            raise ValueError(
                f"precompensation applies to the uplink direction only; "
                f"this scenario has direction={self.direction!r}")

    @property
    def tx_terminal(self) -> Terminal:
        '''The transmit terminal for the link direction (see the module docstring).'''
        return self.space if self.direction == "downlink" else self.ground

    @property
    def rx_terminal(self) -> Terminal:
        '''The receive terminal for the link direction (see the module docstring).'''
        return self.ground if self.direction in ("downlink", "retro") else self.space


@dataclass
class TerrestrialScenario:
    '''
    A terrestrial (horizontal-path) link case: a near terminal + a far terminal.

    Both ends are on the ground, so the terminals are named for the path ends,
    not ground/space. The link is one-way, but the path is reciprocal, so
    `direction` selects which end transmits: "forward" (the default) gives
    tx = near (the local end) and rx = far (the remote end); "reverse" swaps
    the two. The channel does not change, because a horizontal path is the same
    in the two directions. The models read tx_terminal / rx_terminal / channel,
    exactly as for a SpaceScenario.

    Parameters:
        near : Terminal
            The local end of the path.
        far : Terminal
            The remote end of the path.
        direction : "forward" | "reverse"
            The transmit end. "forward" transmits from near, "reverse" from far.
        channel : TerrestrialChannel
            The horizontal propagation channel (path length, attenuation, Cn2).
        availability_target : float
            Target link availability (0-1).
    '''
    near: Terminal
    far: Terminal
    direction: TerrestrialDirection = "forward"
    channel: TerrestrialChannel = field(default_factory=TerrestrialChannel)
    availability_target: float = 0.99

    @property
    def tx_terminal(self) -> Terminal:
        '''The transmit terminal: near on a forward link, far on a reverse link.'''
        return self.far if self.direction == "reverse" else self.near

    @property
    def rx_terminal(self) -> Terminal:
        '''The receive terminal: far on a forward link, near on a reverse link.'''
        return self.near if self.direction == "reverse" else self.far


if __name__ == '__main__':
    from .terminal import Transmitter, Aperture

    ground = Terminal(aperture_m=0.7, transmitter=Transmitter(waist_m=0.1))
    space = Terminal(aperture_m=0.05, detector=Aperture())

    # --- space family -------------------------------------------------------
    up = SpaceScenario(ground=ground, space=space, direction="uplink")
    assert up.tx_terminal is ground and up.rx_terminal is space

    down = SpaceScenario(ground=ground, space=space, direction="downlink")
    assert down.tx_terminal is space and down.rx_terminal is ground

    retro = SpaceScenario(ground=ground, space=space, direction="retro")
    assert retro.tx_terminal is ground and retro.rx_terminal is ground

    assert up.channel.altitude_m == 600e3 and up.availability_target == 0.99
    assert up.channel.site.cn2_ground == 1.7e-14

    # Pre-compensation source: None by default; a DownlinkBeacon or LaserGuideStar
    # names the uplink correction reference.
    assert up.precompensation is None
    beacon_up = SpaceScenario(ground=ground, space=space, direction="uplink",
                              precompensation=DownlinkBeacon())
    assert isinstance(beacon_up.precompensation, DownlinkBeacon)
    lgs_up = SpaceScenario(ground=ground, space=space, direction="uplink",
                           precompensation=LaserGuideStar())
    assert isinstance(lgs_up.precompensation, LaserGuideStar)
    assert lgs_up.precompensation.altitude_m == 90e3

    # A non-uplink scenario refuses a pre-compensation source (backlog I-4).
    for bad_direction in ("downlink", "retro"):
        try:
            SpaceScenario(ground=ground, space=space, direction=bad_direction,
                          precompensation=DownlinkBeacon())
        except ValueError:
            pass
        else:
            raise AssertionError(
                f"precompensation on {bad_direction} did not raise")

    # --- terrestrial family -------------------------------------------------
    near = Terminal(aperture_m=0.1, transmitter=Transmitter(waist_m=0.02))
    far = Terminal(aperture_m=0.1, detector=Aperture())
    terr = TerrestrialScenario(near=near, far=far,
                               channel=TerrestrialChannel(path_length_m=5e3,
                                                          attenuation_db_per_km=0.5,
                                                          cn2=1e-14))
    assert terr.tx_terminal is near and terr.rx_terminal is far
    assert terr.channel.path_length_m == 5e3 and terr.channel.cn2 == 1e-14
    assert terr.channel.attenuation_db_per_km == 0.5
    # The default direction is forward: tx = near, rx = far.
    assert terr.direction == "forward"
    # The reverse direction swaps the two roles and keeps the same channel.
    rev = TerrestrialScenario(near=near, far=far, direction="reverse",
                              channel=terr.channel)
    assert rev.tx_terminal is far and rev.rx_terminal is near
    assert rev.channel is terr.channel
    for scn in (up, terr):
        assert isinstance(scn.tx_terminal, Terminal)
        assert isinstance(scn.rx_terminal, Terminal)

    print("SpaceScenario:", up)
    print("TerrestrialScenario:", terr)
    print("self-check passed")
