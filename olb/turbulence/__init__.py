'''
Pure turbulence physics for optical link budgets.

This package holds the direction-agnostic turbulence physics. It imports only
numpy, scipy, and olb._deps. It does not import the results, the scenario, the
assumptions, the models, or the links.
'''

from .profiles import DEFAULT_HS, default_cn2_profile, get_c2n
from .scintillation import (plane_wave_scintillation_index,
                            aperture_averaging_factor,
                            aperture_averaged_scintillation_index)
from .coupled_flux import coupled_flux_montecarlo

__all__ = [
    "DEFAULT_HS", "default_cn2_profile", "get_c2n",
    "plane_wave_scintillation_index",
    "aperture_averaging_factor",
    "aperture_averaged_scintillation_index",
    "coupled_flux_montecarlo",
]
