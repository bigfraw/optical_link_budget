'''
optical_link_budget (olb)
=========================

Optical ground-to-space laser link budgets with atmospheric propagation, fade
statistics and Monte Carlo.

The main parts:
    Scenario  -- pure-data description of a link case (olb.scenario)
    Terminal  -- one optical terminal; ALL terminal hardware lives here
    geometry  -- CircularOrbit (analytic) or TLEPass (skyfield)
    Term      -- one budget line, with mean / analytic-quantile / sampler views
    Budget    -- a collection of Terms -> table, analytic fade, or Monte Carlo

A Scenario holds two terminals (ground and space), a Channel (the propagation
channel: the site plus the orbit), and the link direction. The channel holds no
hardware; a terminal parameter can only be set through a Terminal. The models in
olb.models turn a Scenario + geometry into Terms. Monte Carlo is not a separate
path. It is the Budget that asks each Term for samples, not means.
'''

from .scenario import Scenario, Site, Channel
from .geometry import CircularOrbit, TLEPass
from .beam import virtual_waist, free_space_radius
from .results import Term, Budget
from .assumptions import Assumptions
from .terminal import Terminal, Transmitter, Aperture, SMF, TipTilt, AO
from .links import uplink_budget, downlink_budget, retro_space_budget, retro_budget
from . import units

__all__ = [
    "Scenario", "Site", "Channel",
    "CircularOrbit", "TLEPass",
    "virtual_waist", "free_space_radius",
    "Term", "Budget",
    "Assumptions",
    "Terminal", "Transmitter", "Aperture", "SMF", "TipTilt", "AO",
    "uplink_budget", "downlink_budget", "retro_space_budget", "retro_budget",
    "units",
]
