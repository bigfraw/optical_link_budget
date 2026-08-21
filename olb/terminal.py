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
Scenario resolves which Terminal transmits and which receives from the link
direction. See olb.scenario.

Approach A: the Compensation stack and the Detector are one physical chain. The
residual wavefront that the Compensation leaves sets the coupling into the
Detector. So the model emits ONE receive-coupling Term, not two. See
olb.models.coupling.
'''

from dataclasses import dataclass, field
from typing import Optional


# --- Transmitter ------------------------------------------------------------

@dataclass
class Transmitter:
    '''
    The transmit source of a terminal.

    The transmitter launches a Gaussian beam through the owning Terminal
    aperture. The launch truncation reads the owning Terminal aperture_m and
    obscuration_ratio, so the transmitter does not repeat them.

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
    '''
    waist_m: float
    power_dbm: Optional[float] = None
    m2: float = 1.0
    divergence_rad: Optional[float] = None


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

    Parameters:
        eta_max : float
            Maximum fibre-to-aperture mode match (flat wavefront).
        sensitivity_dbm : float, optional
            Required received power [dBm]. None if only losses matter.
    '''
    eta_max: float = 0.8145
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
    detector: Optional[object] = None
    compensation: list = field(default_factory=list)


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

    # A detector carries the receive sensitivity.
    smf = Terminal(aperture_m=0.7, detector=SMF(sensitivity_dbm=-40.0),
                   compensation=[TipTilt(), AO(n_modes=60)])
    assert isinstance(smf.detector, SMF)
    assert smf.detector.eta_max == 0.8145 and smf.detector.sensitivity_dbm == -40.0
    assert isinstance(smf.compensation[0], TipTilt)
    assert smf.compensation[1].n_modes == 60
    assert Aperture(sensitivity_dbm=-45.0).sensitivity_dbm == -45.0

    # Two terminals do not share one compensation list (default_factory works).
    a = Terminal(aperture_m=0.5)
    b = Terminal(aperture_m=0.5)
    a.compensation.append(TipTilt())
    assert b.compensation == [], "each Terminal must own its own list"

    print("Terminal:", smf)
    print("self-check passed")
