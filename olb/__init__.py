'''
optical_link_budget (olb)
=========================

Optical ground-to-space laser link budgets with atmospheric propagation, fade
statistics and Monte Carlo.

The main parts:
    SpaceScenario / TerrestrialScenario  -- pure-data link cases (olb.scenario)
    Terminal  -- one optical terminal; ALL terminal hardware lives here
    geometry  -- CircularOrbit (analytic) or TLEPass (skyfield)
    Term      -- one budget line, with mean / analytic-quantile / sampler views
    Budget    -- a collection of Terms -> table, analytic fade, or Monte Carlo

A SpaceScenario holds a ground and a space terminal, a Channel (site + orbit),
and the link direction. A TerrestrialScenario holds a near and a far terminal
and a TerrestrialChannel (site + horizontal path). Both expose the same
tx_terminal / rx_terminal / channel interface, so the models are shared. The
channel holds no hardware; a terminal parameter can only be set through a
Terminal. The models
in olb.models turn a scenario + geometry into Terms. Monte Carlo is not a
separate path. It is the Budget that asks each Term for samples, not means.
'''

from .scenario import (SpaceScenario, TerrestrialScenario, Site,
                       Channel, TerrestrialChannel)
from .geometry import CircularOrbit, TLEPass, HorizontalPath
from .beam import virtual_waist, free_space_radius
from .results import Term, Budget
from .assumptions import Assumptions
from .terminal import Terminal, Transmitter, Aperture, SMF, TipTilt, AO
from .links import (uplink_budget, downlink_budget, retro_space_budget,
                    retro_budget, terrestrial_budget)
from . import units

__all__ = [
    "SpaceScenario", "TerrestrialScenario",
    "Site", "Channel", "TerrestrialChannel",
    "CircularOrbit", "TLEPass", "HorizontalPath",
    "virtual_waist", "free_space_radius",
    "Term", "Budget",
    "Assumptions",
    "Terminal", "Transmitter", "Aperture", "SMF", "TipTilt", "AO",
    "uplink_budget", "downlink_budget", "retro_space_budget", "retro_budget",
    "terrestrial_budget",
    "units",
]
