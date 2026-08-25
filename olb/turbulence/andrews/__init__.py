'''
The Andrews and Phillips foundation layer: pure irradiance and turbulence physics.

Every module in this package follows the book:
    L. C. Andrews and R. L. Phillips, "Laser Beam Propagation through Random
    Media", 2nd ed., SPIE Press (2005). DOI: 10.1117/3.626196
Each function names its section, its equation number, and its printed page.

The package holds physics only. No module here imports a scenario, a terminal, a
Term, or a link. No function here returns decibels. The model layer turns a face
into a dB loss. See olb/models/fade.py.

The nine modules:
    aperture.py       aperture averaging of the irradiance flux (Ch. 10.3)
    beam.py           Gaussian-beam parameters, any input curvature (Ch. 4)
    distributions.py  irradiance PDFs and the fade faces (Ch. 9.9, Ch. 11.3)
    paths.py          slant paths and the satellite link (Ch. 12)
    scintillation.py  Rytov variance and the two log variances (Ch. 8, Ch. 9)
    spectra.py        refractive-index spectra (Ch. 3)
    structure.py      wave structure function, coherence radius (Ch. 6, App. III)
    temporal.py       temporal spectra and the fade rate (Ch. 8.5, Ch. 9.8)
    wander.py         beam wander and the pointing error (Ch. 6.6, Ch. 12.6)
'''

from .aperture import (
    averaged_index,
    averaging_factor,
    d_param,
    omega_g,
    plane_weak_averaging_fit,
)
from .beam import BeamParams, beam_params, effective_beam_params, wavenumber
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
from .paths import (
    OUTER_SCALE_MODELS,
    ZENITH_LIMIT_DEG,
    bufton_wind,
    downlink_scintillation_index,
    hufnagel_valley,
    isoplanatic_angle,
    mu,
    outer_scale_profile,
    point_ahead_angle,
    rms_wind,
    sec_zeta,
    uplink_coherence_radius,
    uplink_scintillation_index,
)
from .scintillation import (
    Q0_CONSTANT,
    QL_CONSTANT,
    WEAK_REGIME_LIMIT,
    beam_rytov_variance,
    large_scale_log_variance,
    rytov_variance,
    scintillation_index,
    small_scale_log_variance,
    two_scale_parameters,
    weak_two_scale_index,
)
from .spectra import (
    EXPONENTIAL_C0,
    KOLMOGOROV_CONSTANT,
    MODIFIED_EQ23_C0,
    MODIFIED_KL,
    SPECTRA,
    TATARSKII_KM,
    VON_KARMAN_C0,
    exponential,
    kolmogorov,
    modified_atmospheric,
    tatarskii,
    von_karman,
)
from .structure import (
    FRIED_OVER_RHO0,
    a_factor,
    angle_of_arrival_variance,
    coherence_radius,
    fried_parameter,
    rms_image_jitter,
    wave_structure_function,
)
from .temporal import (
    GREENWOOD_CONSTANT,
    GREENWOOD_R0_CONSTANT,
    coherence_time,
    fresnel_frequency,
    greenwood_frequency,
    irradiance_temporal_spectrum,
    quasi_frequency,
    taylor_wavenumber,
)
from .wander import (
    C0_DEFAULT,
    CR_DEFAULT,
    WANDER_CONSTANT,
    WANDER_CONSTANT_COLLIMATED,
    beam_wander_variance,
    beam_wander_variance_slant,
    long_term_beam_radius,
    plane_fried_parameter_slant,
    pointing_error_variance,
    pointing_error_variance_slant,
    short_term_beam_radius,
    spherical_fried_parameter,
)

__all__ = [
    # aperture.py
    "averaged_index",
    "averaging_factor",
    "d_param",
    "omega_g",
    "plane_weak_averaging_fit",
    # beam.py
    "BeamParams",
    "beam_params",
    "effective_beam_params",
    "wavenumber",
    # distributions.py
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
    # paths.py
    "OUTER_SCALE_MODELS",
    "ZENITH_LIMIT_DEG",
    "bufton_wind",
    "downlink_scintillation_index",
    "hufnagel_valley",
    "isoplanatic_angle",
    "mu",
    "outer_scale_profile",
    "point_ahead_angle",
    "rms_wind",
    "sec_zeta",
    "uplink_coherence_radius",
    "uplink_scintillation_index",
    # scintillation.py
    "Q0_CONSTANT",
    "QL_CONSTANT",
    "WEAK_REGIME_LIMIT",
    "beam_rytov_variance",
    "large_scale_log_variance",
    "rytov_variance",
    "scintillation_index",
    "small_scale_log_variance",
    "two_scale_parameters",
    "weak_two_scale_index",
    # spectra.py
    "EXPONENTIAL_C0",
    "KOLMOGOROV_CONSTANT",
    "MODIFIED_EQ23_C0",
    "MODIFIED_KL",
    "SPECTRA",
    "TATARSKII_KM",
    "VON_KARMAN_C0",
    "exponential",
    "kolmogorov",
    "modified_atmospheric",
    "tatarskii",
    "von_karman",
    # structure.py
    "FRIED_OVER_RHO0",
    "a_factor",
    "angle_of_arrival_variance",
    "coherence_radius",
    "fried_parameter",
    "rms_image_jitter",
    "wave_structure_function",
    # temporal.py
    "GREENWOOD_CONSTANT",
    "GREENWOOD_R0_CONSTANT",
    "coherence_time",
    "fresnel_frequency",
    "greenwood_frequency",
    "irradiance_temporal_spectrum",
    "quasi_frequency",
    "taylor_wavenumber",
    # wander.py
    "C0_DEFAULT",
    "CR_DEFAULT",
    "WANDER_CONSTANT",
    "WANDER_CONSTANT_COLLIMATED",
    "beam_wander_variance",
    "beam_wander_variance_slant",
    "long_term_beam_radius",
    "plane_fried_parameter_slant",
    "pointing_error_variance",
    "pointing_error_variance_slant",
    "short_term_beam_radius",
    "spherical_fried_parameter",
]
