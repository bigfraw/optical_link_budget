'''
Per-direction Terms and budget assembly for optical link budgets.

This package builds the uplink budget, the downlink budget, and the retro
budget. It imports the model factories and the pure turbulence physics.
'''

from .uplink import uplink_turbulence_term, uplink_budget
from .downlink import downlink_scintillation_term, downlink_budget
from .retro_space import retro_space_budget
# Old name kept for backward compatibility. See retro_space.retro_space_budget.
retro_budget = retro_space_budget
from .terrestrial import terrestrial_scintillation_term, terrestrial_budget

__all__ = [
    "uplink_turbulence_term", "uplink_budget",
    "downlink_scintillation_term", "downlink_budget",
    "retro_space_budget", "retro_budget",
    "terrestrial_scintillation_term", "terrestrial_budget",
]
