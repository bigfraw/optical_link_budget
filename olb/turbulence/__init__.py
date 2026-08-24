'''
Pure turbulence physics for optical link budgets.

This package holds the turbulence physics kernels. Some kernels use a
link-specific simplification: the plane-wave form for the space-to-ground
downlink, and the beam-wave form for the uplink. It imports only numpy, scipy,
and olb._deps. It does not import the results, the scenario, the assumptions, the
models, or the links.
'''

from .profiles import DEFAULT_HS, default_cn2_profile, get_c2n
from .plane_wave_scintillation import (plane_wave_scintillation_index,
                            aperture_averaging_factor,
                            aperture_averaged_scintillation_index,
                            sigma1_rytov,
                            coherence_radius,
                            plane_wave_scintillation_index_closed,
                            aperture_averaged_index_andrews,
                            aperture_averaging_factor_weak,
                            aperture_averaging_factor_weak_inner,
                            aperture_averaging_factor_strong)
from .gaussian_fried import (gaussian_fried_parameter,
                             gaussian_fried_parameter_profile,
                             plane_wave_coherence_radius)
from .beam_wave_scintillation import (on_axis_scintillation_index,
                                 radial_scintillation_index,
                                 gaussian_scintillation_index)
from .uplink_flux import coupled_flux_montecarlo

__all__ = [
    "DEFAULT_HS", "default_cn2_profile", "get_c2n",
    "plane_wave_scintillation_index",
    "aperture_averaging_factor",
    "aperture_averaged_scintillation_index",
    "sigma1_rytov",
    "coherence_radius",
    "plane_wave_scintillation_index_closed",
    "aperture_averaged_index_andrews",
    "aperture_averaging_factor_weak",
    "aperture_averaging_factor_weak_inner",
    "aperture_averaging_factor_strong",
    "gaussian_fried_parameter",
    "gaussian_fried_parameter_profile",
    "plane_wave_coherence_radius",
    "on_axis_scintillation_index",
    "radial_scintillation_index",
    "gaussian_scintillation_index",
    "coupled_flux_montecarlo",
]
