'''
Terminal: the pure-data description of one optical terminal.

This module holds the hardware of one optical terminal as plain dataclasses. It
computes no physics. ALL terminal hardware lives here. The channel (see
olb.scenario) holds no hardware. A terminal parameter can only be set through a
Terminal object.

A Terminal groups a telescope aperture, an optional Transmitter, an optional
wavefront Compensation stack, and an optional Detector front end. The models read
a Terminal and return Terms. The Terminal does not import the models. The data
moves in one direction, from the inputs to the models.

One Terminal serves both link directions. On an uplink the ground Terminal
transmits and the space Terminal receives. On a downlink the roles swap. The
scenario resolves which Terminal transmits and which receives (a SpaceScenario
from its direction; a TerrestrialScenario from near/far). See olb.scenario.

Approach A: the Compensation stack and the Detector are one physical chain. The
residual wavefront that the Compensation leaves sets the coupling into the
Detector. So the model emits ONE receive-coupling Term, not two. See
olb.models.coupling.
'''

from dataclasses import dataclass, field
from typing import Optional, Union


# --- Transmitter ------------------------------------------------------------

@dataclass
class Transmitter:
    '''
    The transmit source of a terminal.

    The transmitter launches a Gaussian beam through a launch aperture. By
    default the launch aperture is the owning Terminal aperture: the launch
    truncation reads the Terminal aperture_m and obscuration_ratio. This is a
    MONOSTATIC terminal, where one aperture transmits and receives.

    For a BISTATIC terminal the transmit beam director is a different aperture
    from the receive telescope. Set aperture_m (and, if it applies,
    obscuration_ratio) on the Transmitter. The launch truncation then reads these
    values, and the Terminal aperture_m and obscuration_ratio describe the
    RECEIVE telescope only. A value of None keeps the monostatic default (the
    Terminal value).

    Parameters:
        waist_m : float
            Transmit Gaussian waist (1/e^2 radius) at the aperture [m].
        power_dbm : float, optional
            Launch power [dBm]. None if only losses matter.
        m2 : float
            Beam quality M^2 (>= 1).
        divergence_rad : float, optional
            Transmit far-field 1/e^2 HALF-angle divergence [rad]. None means
            collimated (the diffraction limit).
        aperture_m : float, optional
            Transmit (beam director) aperture diameter [m]. None means the
            transmitter shares the owning Terminal aperture (monostatic).
        obscuration_ratio : float, optional
            Central obscuration ratio of the transmit aperture. None means the
            transmitter shares the owning Terminal obscuration_ratio. Set 0.0 for
            an unobscured beam director on a terminal whose receive telescope is
            obscured.
    '''
    waist_m: float
    power_dbm: Optional[float] = None
    m2: float = 1.0
    divergence_rad: Optional[float] = None
    aperture_m: Optional[float] = None
    obscuration_ratio: Optional[float] = None


# --- Detector front ends ----------------------------------------------------

@dataclass
class Aperture:
    '''
    Power-in-bucket detector.

    An aperture (bucket) detector integrates the intensity over the aperture. It
    is phase-insensitive. So a Compensation stack does not change its coupling.
    Use it for parity with the plain downlink budget.

    Parameters:
        sensitivity_dbm : float, optional
            Required received power [dBm]. None if only losses matter.
    '''
    sensitivity_dbm: Optional[float] = None


@dataclass
class SMF:
    '''
    Single-mode-fibre detector.

    A single-mode fibre couples only the field that matches the fibre mode. The
    coupling falls when turbulence distorts the wavefront. eta_max is the maximum
    fibre-to-aperture mode match for an unobscured circular aperture with a flat
    wavefront.

    A focusing optic of focal length f puts the collected beam onto the fibre tip.
    Set focal_length_m and mode_field_radius_m to derive eta_max from the optics.
    The model then computes the coupling parameter
        a = pi*(D/2)*w_m / (lambda*f)
    and eta_max(a) from the mode-overlap curve. The peak eta_max=0.8145 is near
    a=1.12. See olb.models.coupling.smf_eta_max_from_a.

    IMPORTANT (single-mode-fibre subtlety): at a fixed a, the focal length f
    cancels in the tilt-to-coupling response. So for a single-mode fibre, f is a
    STATIC knob only: it sets eta_max through a. It does NOT change the angular
    sensitivity on its own. To change eta_max, change a. See the walk-off Term
    in olb.models.coupling.

    Parameters:
        eta_max : float
            Maximum fibre-to-aperture mode match (flat wavefront). Used when
            focal_length_m is None.
        sensitivity_dbm : float, optional
            Required received power [dBm]. None if only losses matter.
        focal_length_m : float, optional
            Focal length of the fibre-coupling optic [m]. None keeps the eta_max
            field unchanged (today's behaviour). A value needs mode_field_radius_m
            so a can be derived.
        mode_field_radius_m : float, optional
            Fibre mode field RADIUS [m] (about 5.2e-6 m for SMF-28 at 1550 nm).
            It sets the fibre mode size for a and for the walk-off Term.
        optimal_focus : bool
            Design the fibre-coupling optic for the best coupling. When True the
            model assumes the optimal coupling parameter a=1.12 (so eta_max=0.8145)
            and derives the focal length f from the mode field radius and the
            aperture: f = pi*(D/2)*w_m/(lambda*1.12). If mode_field_radius_m is
            None, it uses the SMF-28 value (5.2e-6 m). Set focal_length_m to
            override the derived value. A bare SMF() (this flag False) is
            unchanged: it stays mean-only, with no walk-off Term.
    '''
    eta_max: float = 0.8145
    sensitivity_dbm: Optional[float] = None
    focal_length_m: Optional[float] = None
    mode_field_radius_m: Optional[float] = None
    optimal_focus: bool = False


@dataclass
class MMF:
    '''
    Multimode-fibre detector (a light bucket in the fibre plane).

    A multimode fibre accepts all the light that the focusing optic puts inside
    its core. The core is a disk of fixed radius core_radius_m in the fibre plane.
    So the coupling is a geometric overlap of the focal spot with the core, not a
    modal overlap. The focal spot has the diffraction radius
        w_s = lambda*f / (pi*(D/2)).
    A received tip-tilt of angle theta moves the spot by f*theta. So the spot
    walks off the core when the tip-tilt is large.

    Unlike a single-mode fibre, nothing cancels the focal length here. The core is
    a fixed physical size, so the core subtends the angular field of view
    core_radius_m/f. So the focal length is a GENUINE free parameter: a longer f
    magnifies the spot motion and narrows the field of view.

    Parameters:
        core_radius_m : float
            Core RADIUS of the multimode fibre in the fibre plane [m].
        focal_length_m : float
            Focal length of the fibre-coupling optic [m].
        sensitivity_dbm : float, optional
            Required received power [dBm]. None if only losses matter.
    '''
    core_radius_m: float
    focal_length_m: float
    sensitivity_dbm: Optional[float] = None


# --- Wavefront compensation stages ------------------------------------------

@dataclass
class TipTilt:
    '''
    Tip-tilt correction stage.

    A tip-tilt stage removes the first three Zernike modes (piston, tip, tilt).
    '''
    pass


@dataclass
class AO:
    '''
    Adaptive-optics correction stage.

    An adaptive-optics stage removes the first n_modes Zernike modes. Use the
    large-order Noll asymptotic residual (see olb.turbulence.ao).
    '''
    n_modes: int = 20


# A detector is one of the front ends; a compensation stage is one of the
# correctors. These aliases give the Terminal fields a concrete type, so a
# type checker knows the members (e.g. detector.eta_max, detector.sensitivity_dbm).
Detector = Union[Aperture, SMF, MMF]
Compensation = Union[TipTilt, AO]


# --- The terminal -----------------------------------------------------------

@dataclass
class Terminal:
    '''
    One optical terminal: aperture + transmitter + compensation + detector.

    All terminal hardware lives on a Terminal. The channel holds no hardware.

    Parameters:
        aperture_m : float
            Telescope aperture diameter [m].
        obscuration_ratio : float
            Central obscuration diameter divided by aperture diameter. 0 means
            unobscured.
        wavelength_m : float
            The terminal operating wavelength [m].
        pointing_jitter_rad : float
            1-sigma tracking (pointing) jitter [rad]. 0 means no jitter.
        transmitter : Transmitter, optional
            The transmit source. None means the terminal only receives.
        detector : Aperture or SMF, optional
            The detector front end. None means no receive-coupling Term.
        compensation : list
            The ordered wavefront-compensation stack. It may be empty. An empty
            stack leaves the piston-removed turbulence.
    '''
    aperture_m: float
    obscuration_ratio: float = 0.0
    wavelength_m: float = 1550e-9
    pointing_jitter_rad: float = 0.0
    transmitter: Optional[Transmitter] = None
    detector: Optional[Detector] = None
    compensation: list[Compensation] = field(default_factory=list)


if __name__ == '__main__':
    # Pure-data self-check. No physics here.
    t = Terminal(aperture_m=0.7)
    assert t.detector is None and t.compensation == [] and t.obscuration_ratio == 0.0
    assert t.transmitter is None and t.wavelength_m == 1550e-9
    assert t.pointing_jitter_rad == 0.0

    # A transmit terminal carries a Transmitter; truncation reads the aperture.
    tx = Terminal(aperture_m=0.15, obscuration_ratio=0.3,
                  transmitter=Transmitter(waist_m=0.12, power_dbm=40.0))
    assert tx.transmitter.waist_m == 0.12 and tx.transmitter.power_dbm == 40.0
    assert tx.transmitter.divergence_rad is None and tx.transmitter.m2 == 1.0
    # Monostatic default: the transmitter borrows the Terminal launch aperture.
    assert tx.transmitter.aperture_m is None
    assert tx.transmitter.obscuration_ratio is None

    # A bistatic terminal: a small beam director transmits, a large telescope
    # receives. The Transmitter carries its own aperture, so the launch
    # truncation does not read the receive-telescope aperture.
    bistatic = Terminal(aperture_m=0.7, obscuration_ratio=0.3,
                        transmitter=Transmitter(waist_m=0.06, aperture_m=0.15,
                                                obscuration_ratio=0.0))
    assert bistatic.aperture_m == 0.7 and bistatic.obscuration_ratio == 0.3
    assert bistatic.transmitter.aperture_m == 0.15
    assert bistatic.transmitter.obscuration_ratio == 0.0

    # A detector carries the receive sensitivity.
    smf = Terminal(aperture_m=0.7, detector=SMF(sensitivity_dbm=-40.0),
                   compensation=[TipTilt(), AO(n_modes=60)])
    assert isinstance(smf.detector, SMF)
    assert smf.detector.eta_max == 0.8145 and smf.detector.sensitivity_dbm == -40.0
    assert isinstance(smf.compensation[0], TipTilt)
    assert smf.compensation[1].n_modes == 60
    assert Aperture(sensitivity_dbm=-45.0).sensitivity_dbm == -45.0

    # The SMF optics fields default to None, so today's behaviour is unchanged.
    assert smf.detector.focal_length_m is None
    assert smf.detector.mode_field_radius_m is None
    # A focal length and a mode field radius set the coupling optics.
    smf_opt = SMF(focal_length_m=0.02, mode_field_radius_m=5.2e-6)
    assert smf_opt.focal_length_m == 0.02 and smf_opt.mode_field_radius_m == 5.2e-6
    assert smf_opt.eta_max == 0.8145            # the default flat-wavefront match

    # An MMF (light bucket) carries a core radius and a focal length.
    mmf = Terminal(aperture_m=0.2, detector=MMF(core_radius_m=25e-6,
                                                focal_length_m=0.05,
                                                sensitivity_dbm=-38.0))
    assert isinstance(mmf.detector, MMF)
    assert mmf.detector.core_radius_m == 25e-6 and mmf.detector.focal_length_m == 0.05
    assert mmf.detector.sensitivity_dbm == -38.0

    # Two terminals do not share one compensation list (default_factory works).
    a = Terminal(aperture_m=0.5)
    b = Terminal(aperture_m=0.5)
    a.compensation.append(TipTilt())
    assert b.compensation == [], "each Terminal must own its own list"

    print("Terminal:", smf)
    print("self-check passed")
