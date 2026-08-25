'''
The Andrews and Phillips foundation layer: pure irradiance and turbulence physics.

Every module in this package follows the book:
    L. C. Andrews and R. L. Phillips, "Laser Beam Propagation through Random
    Media", 2nd ed., SPIE Press (2005). DOI: 10.1117/3.626196
Each function names its section, its equation number, and its printed page.

The package holds physics only. No module here imports a scenario, a terminal, a
Term, or a link. No function here returns decibels. The model layer turns a face
into a dB loss. See olb/models/fade.py.
'''

from .distributions import (
    MODELS,
    expected_number_of_fades,
    fade_threshold_irradiance,
    gamma_gamma_cdf,
    gamma_gamma_mean_log,
    gamma_gamma_params,
    gamma_gamma_pdf,
    gamma_gamma_quantile,
    gamma_gamma_rvs,
    gamma_gamma_scintillation_index,
    k_cdf,
    k_mean_log,
    k_params,
    k_pdf,
    k_quantile,
    k_rvs,
    k_scintillation_index,
    lognormal_cdf,
    lognormal_mean_log,
    lognormal_params,
    lognormal_quantile,
    lognormal_rician_pdf,
    lognormal_rvs,
    mean_fade_time,
    probability_of_fade,
)

__all__ = [
    "MODELS",
    "expected_number_of_fades",
    "fade_threshold_irradiance",
    "gamma_gamma_cdf",
    "gamma_gamma_mean_log",
    "gamma_gamma_params",
    "gamma_gamma_pdf",
    "gamma_gamma_quantile",
    "gamma_gamma_rvs",
    "gamma_gamma_scintillation_index",
    "k_cdf",
    "k_mean_log",
    "k_params",
    "k_pdf",
    "k_quantile",
    "k_rvs",
    "k_scintillation_index",
    "lognormal_cdf",
    "lognormal_mean_log",
    "lognormal_params",
    "lognormal_quantile",
    "lognormal_rician_pdf",
    "lognormal_rvs",
    "mean_fade_time",
    "probability_of_fade",
]
