'''
Coupled-flux Monte Carlo wrapper for a LEO uplink (beam wander + scintillation).

This module wraps the Dios et al. coupled-flux Monte Carlo
(``coupled_flux_montecarlo``). It runs the kernel for one elevation and rescales
the flux to the diffraction-limited baseline. The wrapper is pure. The Term
factory lives in olb.links.uplink.

Diffraction-limited rescaling (essential): ``coupled_flux`` normalises the
on-axis irradiance to the beam's own short-term waist ``w_st``. The short-term
waist already includes the turbulence spread beyond the diffraction limit
(w_st = sqrt(w_diff**2 + turbulence_term)). This self-normalisation removes the
beam-broadening loss from the result. For a 600 km, ~1 m-waist uplink the loss
is approximately 10 dB. This is a large error, not a small one. The code
rescales the flux samples by (w_diff / w_st)**2. This puts the samples back on
the diffraction-limited baseline, which is the same reference that the geometric
term uses. The two terms are then directly additive. The per-sample dB loss is
then -10*log10(flux).

Validity: the Rytov model is a weak-fluctuation model. When the mean
log-amplitude variance sigma2_x exceeds WEAK_FLUCTUATION_LIMIT, the
scintillation approaches saturation and the numbers are not trustworthy. The
code carries a ``weak_fluctuation_valid`` flag and sigma2_x in the result, and
it gives a warning.
'''

import warnings

import numpy as np

from .._deps import coupled_flux_montecarlo, gaussz
from .profiles import DEFAULT_HS

WEAK_FLUCTUATION_LIMIT = 0.6   # log-amplitude variance limit; above it the Rytov model is not valid (saturation)


def _flux_result(w0, elevation_deg, range_m, wavelength, hs, cn2_profile,
                 hv57_A, n_samples, n_apertures):
    '''
    Run the coupled-flux MC for one elevation and rescale to the
    diffraction-limited baseline (see module docstring).

    Returns:
        dict
            The kernel output with ``Is_summed`` rescaled, plus
            ``w_diffraction_limited`` and ``weak_fluctuation_valid`` added.
    '''
    try:
        result = coupled_flux_montecarlo(
            w0=w0, elevation=elevation_deg, labda=wavelength, L=range_m,
            hs=hs, n_samples=n_samples, n_apertures=n_apertures,
            cn2_profile=cn2_profile, hv57_A=hv57_A,
        )
    except (ImportError, ModuleNotFoundError) as e:
        # The kernel builds its HV57 Cn2 profile from the `fast` package when
        # cn2_profile is None. If fast is missing that import fails here.
        raise ImportError(
            "coupled_flux_montecarlo needs the `fast` package to build the "
            "HV57 Cn2 profile. Run `pip install fast-atmosphere`, or pass an "
            "explicit cn2_profile to uplink_turbulence_term()."
        ) from e

    w_diff = gaussz(w0, range_m, wavelength)
    result["Is_summed"] = result["Is_summed"] * (w_diff / result["w_st"]) ** 2
    result["w_diffraction_limited"] = w_diff
    result["weak_fluctuation_valid"] = bool(result["sigma2_x_mean"] < WEAK_FLUCTUATION_LIMIT)
    if not result["weak_fluctuation_valid"]:
        warnings.warn(
            f"log-amplitude variance sigma2_x={result['sigma2_x_mean']:.2f} >= "
            f"{WEAK_FLUCTUATION_LIMIT} -- scintillation approaching saturation, "
            "turbulence loss not trustworthy (Rytov weak-fluctuation model exceeded)."
        )
    return result


if __name__ == '__main__':
    # Pure-physics self-check. Use plain numeric inputs; this module must not
    # import the scenario or the geometry. Seed the global RNG once because the
    # kernel draws from np.random.
    w0 = 1.0
    lam = 1550e-9
    hs = DEFAULT_HS
    range_m = 600e3            # zenith slant range for a 600 km orbit
    np.random.seed(0)

    weak_cn2 = 1e-18 * np.ones_like(hs)     # negligible turbulence
    strong_cn2 = 1e-15 * np.ones_like(hs)   # strong, sigma2_x finite but past the limit

    # Rescaling: the flux uses the diffraction-limited waist as the baseline.
    r_weak = _flux_result(w0, 90.0, range_m, lam, hs, weak_cn2, 1.7e-14, 1000, 1)
    assert np.isclose(r_weak["w_diffraction_limited"], gaussz(w0, range_m, lam))
    assert np.all(np.isfinite(r_weak["Is_summed"]))

    # weak_fluctuation_valid follows the threshold: negligible Cn2 -> valid,
    # strong Cn2 -> invalid.
    assert r_weak["weak_fluctuation_valid"] is True
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r_strong = _flux_result(w0, 90.0, range_m, lam, hs, strong_cn2,
                                1.7e-14, 1000, 1)
    assert r_strong["weak_fluctuation_valid"] is False

    print(f"w_diffraction_limited = {r_weak['w_diffraction_limited']:.4f} m")
    print(f"weak Cn2   -> weak_fluctuation_valid={r_weak['weak_fluctuation_valid']} "
          f"sigma2_x={r_weak['sigma2_x_mean']:.4f}")
    print(f"strong Cn2 -> weak_fluctuation_valid={r_strong['weak_fluctuation_valid']} "
          f"sigma2_x={r_strong['sigma2_x_mean']:.4f}")
    print("self-check passed.")
