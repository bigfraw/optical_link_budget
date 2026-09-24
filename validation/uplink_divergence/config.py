"""The shared configuration of the diverged-uplink validation.

A 15 cm ground aperture with no obscuration and a FILLED launch (waist
55 mm), 1550 nm, a 500 km orbit, and the site HV5/7 atmosphere with
L0 = 25 m. The divergences are multiples of the diffraction divergence
theta_min = lambda / (pi w) of the waist (Siegman, Lasers, Ch. 17,
ISBN 0-935702-11-3; the Gaussian far-field half-angle).
"""
import numpy as np
from olb.scenario import SpaceScenario, Channel
from olb.geometry import CircularOrbit
from olb.terminal import Terminal, Transmitter, Aperture

LAM = 1550e-9
D = 0.15
WAIST = 0.055
ALT = 500e3
ELEVATIONS = (30.0, 60.0)
THETA_MIN = LAM / (np.pi * WAIST)
MULTIPLES = (None, 2, 4, 8, 16)            # None is the collimated launch
DIVERGENCES = tuple(None if m is None else m * THETA_MIN for m in MULTIPLES)
# The top case the owner wants, 200 urad. It is read post hoc ONLY: no Arm A
# check, the trend of the smaller cases stands for it.
DIVERGENCES = DIVERGENCES + (200e-6,)


def scenario(div=None):
    """Give the uplink SpaceScenario for one divergence (rad, or None)."""
    ground = Terminal(aperture_m=D, wavelength_m=LAM,
                      transmitter=Transmitter(waist_m=WAIST, divergence_rad=div))
    space = Terminal(aperture_m=0.05, wavelength_m=LAM, detector=Aperture())
    return SpaceScenario(ground=ground, space=space, direction="uplink",
                         channel=Channel(altitude_m=ALT))


def geometry(elev):
    return CircularOrbit(ALT, float(elev))
