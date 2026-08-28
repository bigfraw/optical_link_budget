'''
Coupled-flux Monte Carlo wrapper for a LEO uplink (beam wander + scintillation).

This module reimplements the short uplink loop of the Dios et al. coupled-flux
Monte Carlo (the shared ``coupled_flux_montecarlo``) for one elevation, and
rescales the flux to the free-space baseline. The reimplementation lets a
deliberate transmit divergence enter the beam-broadening baseline only (through
the ``w_free`` override on the waists), without a divergence argument on the
shared kernel and without editing it. The wrapper is pure. The Term factory
lives in olb.links.uplink.

Free-space rescaling (essential): the coupled flux normalises the on-axis
irradiance to the beam's own short-term waist ``w_st``. The short-term waist
already includes the turbulence spread beyond the free-space beam
(w_st = sqrt(w_free**2 + turbulence_term)). This self-normalisation removes the
beam-broadening loss from the result. For a 600 km, ~1 m-waist uplink the loss
is approximately 10 dB. This is a large error, not a small one. The code
rescales the flux samples by (w_free / w_st)**2. This puts the samples back on
the free-space baseline (``w_free_div``, the diverged free-space width when a
divergence is set), which is the same reference that the geometric term uses.
The two terms are then directly additive. The per-sample dB loss is then
-10*log10(flux).

Divergence enters everywhere the beam geometry matters: the broadening baseline,
the short/long-term waists (through the ``w_free`` override), and the
scintillation index (through the receiver-plane Lambda and Theta of the diverged
beam -- see ``_scintillation_beam``). A diverged beam is larger and more
spherical-wave-like, so it scintillates less.

Validity: the Rytov model is a weak-fluctuation model. When the mean
log-amplitude variance sigma2_x exceeds WEAK_FLUCTUATION_LIMIT, the
scintillation approaches saturation and the numbers are not trustworthy. The
code carries a ``weak_fluctuation_valid`` flag and sigma2_x in the result, and
it gives a warning.
'''

import warnings

import numpy as np

from ..beam import gaussz, zR
from .._deps import (coupled_flux_montecarlo,
                     spherical_wave_coherence_diameter, short_term_beam_waist,
                     long_term_beam_waist, beam_wander_variance,
                     coupled_flux_sample, on_axis_irradiance)
from ..beam import free_space_radius, launch_curvature
from .profiles import DEFAULT_HS

# Log-amplitude variance limit. Above it the Rytov model is not valid
# (saturation). Andrews and Phillips, Laser Beam Propagation through Random
# Media, 2nd ed. (2005), DOI 10.1117/3.626196, Ch. 8, the text below Eq. (23),
# printed pp. 264-265, give the weak limit as sigma_R^2 < 1. With
# sigma_I^2 = 4 sigma_x^2 (Ch. 8, Eq. (13)) that limit is sigma_x^2 < 0.25.
WEAK_FLUCTUATION_LIMIT = 0.25


def _scintillation_beam(w0, L, wavelength, divergence_rad):
    '''
    Receiver-plane Gaussian-beam parameters for the Dios scintillation index,
    for a beam with transmit far-field half-angle divergence ``divergence_rad``.

    The Dios scintillation integrals read the beam only through two numbers at
    the receiver: the Fresnel parameter Lambda (through the free-space width wL)
    and the curvature parameter Theta (through an effective Rayleigh range Z0).
    A deliberately diverged beam is, exactly, a Gaussian beam from a virtual
    waist behind the aperture, so it has its own Lambda and Theta. This function
    returns the wL and the effective Z0 that carry them, so the diverged beam
    feeds the shared scintillation index through its ordinary ``wL``/``Z0``
    arguments -- no change to the shared kernel.

    Transmitter-plane parameters (aperture, z=0):
        Lambda0 = 2*L/(k*w0^2)                  # w0 is the physical aperture radius
        Theta0  = 1 - L/F0                       # F0 = phase-front radius at aperture
    A diverging wavefront has F0 < 0, so Theta0 > 1. The collimated case is
    F0 = infinity -> Theta0 = 1. Transform to the receiver plane:
        Theta = Theta0 / (Theta0^2 + Lambda0^2)
        Z0    = L / sqrt(1/Theta - 1)            # so that 1/(1+(L/Z0)^2) = Theta
    The collimated case returns Z0 = zR(w0) and wL = gaussz(w0, L) exactly.

    Parameters:
        w0 : float
            Physical beam radius at the transmit aperture [m].
        L : float
            Range from the aperture to the receiver [m].
        wavelength : float
            Wavelength [m].
        divergence_rad : float or None
            Transmit far-field half-angle divergence [rad]. None = collimated.

    Returns:
        tuple
            (wL, Z0) : free-space width at the receiver [m] and the effective
            Rayleigh range [m] that carry the receiver-plane Lambda and Theta.
    '''
    # ponytail: TODO validate this diverged (Theta, Lambda) feed against the
    # Andrews & Phillips closed-form on-axis Gaussian-beam scintillation index
    # (weak fluctuation, constant Cn2). Feed the Dios integrator a constant Cn2
    # over a horizontal-equivalent path and assert on_axis_scintillation_index
    # matches 3.86*sigma_R^2*{0.40[(1+2Theta)^2+4Lambda^2]^(5/12)*
    # cos[(5/6)atan((1+2Theta)/(2Lambda))] - (11/16)Lambda^(5/6)} to a few
    # percent, for the collimated AND the diverged beam. Not added yet.
    k = 2 * np.pi / wavelength
    lambda0 = 2 * L / (k * w0 ** 2)
    # launch_curvature gives f0 = inf (collimated) or f0 < 0 (diverging), in
    # this repo's Theta0 = 1 - L/f0 convention. One shared implementation.
    f0 = launch_curvature(w0, divergence_rad, wavelength)
    theta0 = 1.0 - L / f0                  # 1.0 collimated, > 1 diverging
    theta = theta0 / (theta0 ** 2 + lambda0 ** 2)
    z0_eff = L / np.sqrt(1.0 / theta - 1.0)
    wL = float(free_space_radius(w0, L, divergence_rad, wavelength))
    return wL, z0_eff


def _flux_result(w0, elevation_deg, range_m, wavelength, hs, cn2_profile,
                 hv57_A, n_samples, n_apertures, divergence_rad=None,
                 sigma_theta_rad=0.0):
    '''
    Run the coupled-flux MC for one elevation and rescale to the free-space
    baseline (see module docstring).

    This reimplements the short uplink loop of the shared
    ``coupled_flux_montecarlo``, so that deliberate transmit divergence enters
    the beam-broadening baseline only. The short- and long-term waists use the
    DIVERGED free-space width ``w_free_div`` (through the ``w_free`` override),
    and the rescale baseline is ``w_free_div`` too -- the same reference the
    diverged geometric term uses. The turbulence term then carries the pure
    turbulence broadening and does not double-count the divergence. The
    scintillation index also reads the diverged beam, through its receiver-plane
    Lambda and Theta (see ``_scintillation_beam``).

    Parameters:
        divergence_rad : float, optional
            Transmit far-field half-angle divergence [rad]. None = collimated.
        sigma_theta_rad : float
            Mechanical pointing (tracking) jitter, per-axis 1-sigma angle [rad].
            It adds to the SAME receiver-plane displacement as the turbulence
            beam wander, so it folds into the wander variance beta2 before the
            per-sample offset is drawn (see below). 0.0 = perfect tracking.

    Returns:
        dict
            ``Is_summed`` (rescaled), ``w_st``, ``w_lt``, ``r0s``,
            ``sigma2_x_mean``, ``w_diffraction_limited`` and
            ``weak_fluctuation_valid``.
    '''
    if cn2_profile is None:
        try:
            from fast import turbulence_models
            cn2_profile = turbulence_models.HV57(hs, A=hv57_A)
        except (ImportError, ModuleNotFoundError) as e:
            raise ImportError(
                "the coupled-flux MC needs the `fast` package to build the "
                "HV57 Cn2 profile. Run `pip install fast-aosim`, or pass an "
                "explicit cn2_profile to uplink_turbulence_term()."
            ) from e

    L = float(range_m)
    k = 2 * np.pi / wavelength
    Z0 = zR(w0, wavelength)
    airmass = 1.0 / np.sin(np.radians(elevation_deg))
    cn2_slant = np.asarray(cn2_profile, dtype=float) * airmass

    # Diverged free-space beam: the broadening baseline. w_free_div at the
    # receiver, and the profile ws_div along the path for the wander integral.
    w_free_div = float(free_space_radius(w0, L, divergence_rad, wavelength))
    ws_div = free_space_radius(w0, hs, divergence_rad, wavelength)
    # Receiver-plane beam that the scintillation index reads: the DIVERGED
    # free-space width and an effective Rayleigh range that carry the diverged
    # Lambda and Theta. The collimated case reduces to gaussz(w0, L) and zR(w0).
    wL_scint, Z0_scint = _scintillation_beam(w0, L, wavelength, divergence_rad)

    r0s = spherical_wave_coherence_diameter(k, L, cn2_slant, hs)
    # w_free override -> waists broaden relative to the DIVERGED free-space beam.
    w_st = short_term_beam_waist(w0, L, Z0, k, r0s, w_free=w_free_div)
    beta2 = beam_wander_variance(L, cn2_slant, ws_div, hs)
    # Mechanical pointing jitter shares the receiver-plane displacement with the
    # turbulence beam wander, as an independent 2-D Gaussian offset. So sum the
    # displacement variances here: the combined per-sample offset then feeds BOTH
    # the Gaussian power falloff (on_axis_irradiance) AND the off-axis Dios
    # scintillation (coupled_flux_sample). beta2 is the total 2-D variance
    # <r^2>; a per-axis 1-sigma jitter angle sigma_theta gives a per-axis
    # displacement sigma_r = sigma_theta*L, which adds 2*sigma_r^2 to the 2-D
    # total. This is variance addition of independent offsets, not new physics,
    # so no extra citation: the wander-offset mechanism itself is Dios (Applied
    # Optics 43 (2004) 3866). This replaces a standalone pointing-loss term on
    # the uplink; adding both would double-count the jitter displacement.
    beta2 = beta2 + 2.0 * (sigma_theta_rad * L) ** 2
    w_lt = long_term_beam_waist(w_st, beta2)

    xis = np.zeros(n_samples)
    betas = np.zeros(n_samples)
    sigma2_xs = np.zeros(n_samples)
    for i in range(n_samples):
        betax = np.random.normal(0, np.sqrt(0.5 * beta2), 1)
        betay = np.random.normal(0, np.sqrt(0.5 * beta2), 1)
        beta = np.sqrt(betax ** 2 + betay ** 2)
        # The scintillation index reads the DIVERGED beam through its
        # receiver-plane width and effective Rayleigh range (wL_scint, Z0_scint).
        # A diverged beam is larger and more spherical-wave-like, so it
        # scintillates less. The collimated case reduces to the ordinary values.
        xi, _, s2x, _, _, _ = coupled_flux_sample(
            beta, cn2_slant, Z0_scint, hs, L, k, wL_scint, w_lt)
        betas[i] = np.squeeze(beta)
        xis[i] = np.squeeze(xi)
        sigma2_xs[i] = np.squeeze(s2x)

    Is = on_axis_irradiance(betas, w_st, xis)
    n_blocks = Is.shape[0] // n_apertures
    Is_summed = np.mean(Is[: n_blocks * n_apertures].reshape(n_blocks, n_apertures), axis=1)

    # Rescale onto the diverged free-space baseline w_free_div, the same
    # reference the diverged geometric term uses -> the two terms stay additive.
    Is_summed = Is_summed * (w_free_div / w_st) ** 2
    sigma2_x_mean = float(np.mean(sigma2_xs))

    result = {
        "w_st": float(w_st),
        "w_lt": float(w_lt),
        "r0s": float(r0s),
        "sigma2_x_mean": sigma2_x_mean,
        "Is_summed": Is_summed,
        "w_diffraction_limited": w_free_div,
        "weak_fluctuation_valid": bool(sigma2_x_mean < WEAK_FLUCTUATION_LIMIT),
    }
    if not result["weak_fluctuation_valid"]:
        warnings.warn(
            f"log-amplitude variance sigma2_x={sigma2_x_mean:.2f} >= "
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

    # Divergence: it enters the broadening baseline. A diverged beam has a wider
    # free-space width and a wider short-term waist, and it dilutes the mean
    # turbulence loss (the flux is spread over a larger beam that turbulence
    # broadens less, in relative terms).
    from ..units import w0_to_div
    theta_min = w0_to_div(w0, lam)
    np.random.seed(0)
    r_coll = _flux_result(w0, 90.0, range_m, lam, hs, weak_cn2, 1.7e-14, 4000, 1)
    np.random.seed(0)
    r_div = _flux_result(w0, 90.0, range_m, lam, hs, weak_cn2, 1.7e-14, 4000, 1,
                         divergence_rad=5 * theta_min)
    assert r_div["w_diffraction_limited"] > r_coll["w_diffraction_limited"]
    assert r_div["w_st"] > r_coll["w_st"]
    loss_coll = -10 * np.log10(np.mean(r_coll["Is_summed"]))
    loss_div = -10 * np.log10(np.mean(r_div["Is_summed"]))
    assert loss_div < loss_coll, (loss_div, loss_coll)

    # Collimated request (None) must match the diffraction-limited free-space width.
    assert np.isclose(r_coll["w_diffraction_limited"], gaussz(w0, range_m, lam))

    # The scintillation index now reads the diverged beam. A diverged beam is
    # larger and more spherical-wave-like, so its log-amplitude variance
    # sigma2_x is LOWER than the collimated beam under the same turbulence. Use a
    # moderate profile so sigma2_x is well above numerical noise.
    moderate_cn2 = 1e-16 * np.ones_like(hs)
    np.random.seed(1)
    m_coll = _flux_result(w0, 90.0, range_m, lam, hs, moderate_cn2, 1.7e-14, 4000, 1)
    np.random.seed(1)
    m_div = _flux_result(w0, 90.0, range_m, lam, hs, moderate_cn2, 1.7e-14, 4000, 1,
                         divergence_rad=5 * theta_min)
    assert m_div["sigma2_x_mean"] < m_coll["sigma2_x_mean"], (
        m_div["sigma2_x_mean"], m_coll["sigma2_x_mean"])

    # The collimated scintillation must be UNCHANGED by the receiver-plane recast
    # -- _scintillation_beam reduces to (gaussz(w0,L), zR(w0)) for a collimated
    # beam. Cross-check the effective parameters against the plain values.
    wL_c, Z0_c = _scintillation_beam(w0, range_m, lam, None)
    assert np.isclose(wL_c, gaussz(w0, range_m, lam))
    assert np.isclose(Z0_c, zR(w0, lam))

    # Pointing jitter: it folds into the wander displacement, so a larger jitter
    # widens the offset distribution -> a deeper mean loss AND a deeper fade,
    # with no separate pointing term. Zero jitter reproduces the no-jitter run.
    np.random.seed(2)
    r_nojit = _flux_result(w0, 90.0, range_m, lam, hs, moderate_cn2, 1.7e-14,
                           8000, 1)
    np.random.seed(2)
    r_jit = _flux_result(w0, 90.0, range_m, lam, hs, moderate_cn2, 1.7e-14,
                         8000, 1, sigma_theta_rad=5e-6)
    loss_nojit = -10 * np.log10(np.mean(r_nojit["Is_summed"]))
    loss_jit = -10 * np.log10(np.mean(r_jit["Is_summed"]))
    fade99_nojit = -10 * np.log10(np.percentile(r_nojit["Is_summed"], 1))
    fade99_jit = -10 * np.log10(np.percentile(r_jit["Is_summed"], 1))
    assert loss_jit > loss_nojit, (loss_jit, loss_nojit)
    assert fade99_jit > fade99_nojit, (fade99_jit, fade99_nojit)
    # Zero jitter is a no-op: same displacement variance -> same result.
    np.random.seed(3)
    r_a = _flux_result(w0, 90.0, range_m, lam, hs, moderate_cn2, 1.7e-14, 2000, 1)
    np.random.seed(3)
    r_b = _flux_result(w0, 90.0, range_m, lam, hs, moderate_cn2, 1.7e-14, 2000, 1,
                       sigma_theta_rad=0.0)
    assert np.allclose(r_a["Is_summed"], r_b["Is_summed"])
    print(f"no jitter   -> loss {loss_nojit:.3f} dB, 99% fade {fade99_nojit:.3f} dB")
    print(f"5 urad jit  -> loss {loss_jit:.3f} dB, 99% fade {fade99_jit:.3f} dB")

    print(f"collimated -> w_free={r_coll['w_diffraction_limited']:.2f} m, "
          f"turbulence loss {loss_coll:.3f} dB, sigma2_x={m_coll['sigma2_x_mean']:.4f}")
    print(f"5x diverged -> w_free={r_div['w_diffraction_limited']:.2f} m, "
          f"turbulence loss {loss_div:.3f} dB, sigma2_x={m_div['sigma2_x_mean']:.4f}")
    print("self-check passed.")
