'''
Per-direction Terms and budget assembly for optical link budgets.

This package builds the uplink budget, the downlink budget, and the retro
budget. It imports the direction-agnostic model factories and the pure
turbulence physics.
'''

from .uplink import uplink_turbulence_term, uplink_budget
from .downlink import downlink_scintillation_term, downlink_budget
from .retro_space import retro_space_budget
from .retro import retro_budget   # backward-compatible alias for retro_space_budget

__all__ = [
    "uplink_turbulence_term", "uplink_budget",
    "downlink_scintillation_term", "downlink_budget",
    "retro_space_budget", "retro_budget",
]
