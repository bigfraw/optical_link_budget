'''
Terminal: the pure-data description of a receive optical terminal.

This module holds the hardware of one optical terminal as plain dataclasses. It
computes no physics. A Terminal groups a telescope aperture, an optional
wavefront Compensation stack, and a Detector front end. The models read a
Terminal and return Terms. The Terminal does not import the models. The data
moves in one direction, from the inputs to the models.

This first slice builds the RECEIVE side of the downlink. The satellite
transmits. The ground station detects with a telescope, an optional adaptive
optics stage, and a single-mode-fibre or aperture front end.

Approach A: the Compensation stack and the Detector are one physical chain. The
residual wavefront that the Compensation leaves sets the coupling into the
Detector. So the model emits ONE receive-coupling Term, not two. See
olb.models.coupling.
'''

from dataclasses import dataclass, field
from typing import Optional


# --- Detector front ends ----------------------------------------------------

@dataclass
class Aperture:
    '''
    Power-in-bucket detector.

    An aperture (bucket) detector integrates the intensity over the aperture. It
    is phase-insensitive. So a Compensation stack does not change its coupling.
    Use it for parity with the plain downlink budget.
    '''
    pass


@dataclass
class SMF:
    '''
    Single-mode-fibre detector.

    A single-mode fibre couples only the field that matches the fibre mode. The
    coupling falls when turbulence distorts the wavefront. eta_max is the maximum
    fibre-to-aperture mode match for an unobscured circular aperture with a flat
    wavefront.
    '''
    eta_max: float = 0.8145


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
    A receive optical terminal: aperture + compensation stack + detector.

    Parameters:
        aperture_m : float
            Telescope aperture diameter [m].
        detector : Aperture or SMF, optional
            The detector front end. None means no receive-coupling Term.
        compensation : list
            The ordered wavefront-compensation stack. It may be empty. An empty
            stack leaves the piston-removed turbulence.
        obscuration_ratio : float
            Central obscuration diameter divided by aperture diameter. 0 means
            unobscured.
        wavelength_m : float, optional
            Wavelength [m]. None falls back to the link or scenario wavelength.
    '''
    aperture_m: float
    detector: Optional[object] = None
    compensation: list = field(default_factory=list)
    obscuration_ratio: float = 0.0
    wavelength_m: Optional[float] = None


if __name__ == '__main__':
    # Pure-data self-check. No physics here.
    t = Terminal(aperture_m=0.7)
    assert t.detector is None and t.compensation == [] and t.obscuration_ratio == 0.0

    smf = Terminal(aperture_m=0.7, detector=SMF(),
                   compensation=[TipTilt(), AO(n_modes=60)])
    assert isinstance(smf.detector, SMF)
    assert smf.detector.eta_max == 0.8145
    assert isinstance(smf.compensation[0], TipTilt)
    assert smf.compensation[1].n_modes == 60

    # Two terminals do not share one compensation list (default_factory works).
    a = Terminal(aperture_m=0.5)
    b = Terminal(aperture_m=0.5)
    a.compensation.append(TipTilt())
    assert b.compensation == [], "each Terminal must own its own list"

    print("Terminal:", smf)
    print("self-check passed")
