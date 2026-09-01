"""Shared setup for the presentation plots.

This module builds ONE hero downlink scenario and a clean light matplotlib
style. Every plot script imports from here, so the story stays consistent
across the slides. Loss is positive dB.
"""
import os
import sys

import numpy as np

# Make olb importable when a script runs from any directory.
_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

from olb import (SpaceScenario, Site, Channel, CircularOrbit, Terminal,
                 Transmitter, Aperture, SMF, TipTilt, AO)
from olb.turbulence.profiles import DEFAULT_HS, default_cn2_profile

# ----------------------------------------------------------------------------
# Paths
# ----------------------------------------------------------------------------
HERE = os.path.dirname(os.path.abspath(__file__))
FIGDIR = os.path.join(HERE, "figures")
DATADIR = os.path.join(HERE, "data")
os.makedirs(FIGDIR, exist_ok=True)
os.makedirs(DATADIR, exist_ok=True)


def figpath(name):
    return os.path.join(FIGDIR, name)


def datapath(name):
    return os.path.join(DATADIR, name)


# ----------------------------------------------------------------------------
# Hero downlink parameters (one LEO pass, 1550 nm)
# ----------------------------------------------------------------------------
WAVELENGTH_M = 1550e-9
ALTITUDE_M = 500e3           # LEO orbit height
HERO_ELEVATION_DEG = 30.0    # the running example elevation
TX_POWER_DBM = 30.0          # 1 W satellite launch

SITE = Site(cn2_ground=1.7e-14, wind_rms_m_s=21.0)
HS = DEFAULT_HS
CN2 = default_cn2_profile(SITE, HS)


GROUND_APERTURE_M = 0.7      # ground receive telescope diameter


def _satellite_tx():
    """The satellite transmit terminal (small aperture, 1 W)."""
    return Terminal(aperture_m=0.10, wavelength_m=WAVELENGTH_M,
                    pointing_jitter_rad=1e-6,
                    transmitter=Transmitter(waist_m=0.04, power_dbm=TX_POWER_DBM))


def downlink():
    """The ONE hero downlink: a 0.7 m ground telescope coupling to single-mode
    fibre, UNCOMPENSATED (no adaptive optics).

    An uncompensated coherent receiver lets every fidelity model the SAME
    physics (turbulent fibre coupling), so the three rungs of the ladder sit on
    one axis. Wave optics (fidelity 2) applies no AO, so this is also the only
    honest way to put it beside the analytic and FAST models.
    """
    ground = Terminal(aperture_m=GROUND_APERTURE_M, wavelength_m=WAVELENGTH_M,
                      pointing_jitter_rad=2e-6,
                      detector=SMF(sensitivity_dbm=-45.0))
    return SpaceScenario(ground=ground, space=_satellite_tx(),
                         direction="downlink",
                         channel=Channel(site=SITE, altitude_m=ALTITUDE_M))


def geometry(elevation_deg=HERO_ELEVATION_DEG):
    return CircularOrbit(ALTITUDE_M, elevation_deg=elevation_deg)


# ----------------------------------------------------------------------------
# Light matplotlib style (standard matplotlib vibes, tidied)
# ----------------------------------------------------------------------------
# A restrained, colour-blind-safe palette.
BLUE = "#1f77b4"
ORANGE = "#ff7f0e"
GREEN = "#2ca02c"
RED = "#d62728"
PURPLE = "#9467bd"
GREY = "#6b6b6b"

# Per-fidelity colours, used everywhere the ladder appears.
FID_COLORS = {0: "#4c72b0", 1: "#dd8452", 2: "#55a868"}
FID_LABELS = {
    0: "Fidelity 0  ·  analytic",
    1: "Fidelity 1  ·  statistical",
    2: "Fidelity 2  ·  wave optics",
}


def use_style():
    import matplotlib as mpl
    mpl.rcParams.update({
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
        "savefig.dpi": 160,
        "savefig.bbox": "tight",
        "font.size": 11,
        "axes.titlesize": 13,
        "axes.titleweight": "bold",
        "axes.labelsize": 11,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.color": "#e6e6e6",
        "grid.linewidth": 0.8,
        "legend.frameon": False,
        "figure.titlesize": 15,
        "figure.titleweight": "bold",
    })
