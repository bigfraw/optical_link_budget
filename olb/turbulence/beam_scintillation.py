'''
Gaussian-beam scintillation index for an uplink, on axis and off axis.

This module gives the scintillation index of a collimated Gaussian beam that
propagates up through turbulence to a receiver. It uses the beam-wave Rytov
theory of Andrews et al., in the uplink form of Dios et al., Applied Optics 43
(2004) 3866, Eqs. (13)-(20). The functions are pure. They integrate a real
Cn2(h) slant profile. They return the flux scintillation index (the normalised
irradiance variance).

The scintillation index is NOT the Fried parameter. The Fried parameter is a
coherence length. The scintillation index is an irradiance variance. Both use
the same beam parameters Theta and Lambda (see gaussian_fried.output_beam_params)
but through different path integrals. So this module shares the beam parameters
with gaussian_fried, not the Fried parameter itself.

The total index at a point that sits a radius r off the beam axis is
    sigma2_I(r, L) = sigma2_I(0, L) + sigma2_Ir(r, L)          (Dios Eq. 13)
The first term is the on-axis (longitudinal) index. The second term is the
radial (off-axis) addition. A point detector at the beam centre sees the
on-axis term only. The beam-wander model feeds r = beta (the wander offset).

Regime: weak-to-moderate turbulence. The model has no saturation. Dios reports
good agreement with a split-step beam-propagation reference up to sigma2_chi ~
0.6. Above that the true index saturates and this model overshoots.

Physics (Dios Eqs. 14-20, collimated beam launched from the ground):
    Theta, Lambda are the receiver-plane beam parameters (Dios Eq. 15).
    A(z) = (Lambda L / k) ((L-z)/L)^2                          (Eq. 17)
    B(z) = (L / k) ((L-z)/L) (Theta + (1-Theta) z/L)           (Eq. 18)
    sigma2_I(0,L) = 4 pi^2 k^2 Gamma(-5/6) 0.033
        * INT Cn2(z) [ A^(5/6) - (A^2+B^2)^(5/12)
                       * cos( (5/6) arctan(B/A) ) ] dz         (Eq. 16)
    sigma2_Ir(r,L) = 4 pi^2 k^2 Gamma(-5/6) 0.033
        * ( 1F1(-5/6, 1, 2 r^2 / W^2(L)) - 1 )
        * INT Cn2(z) A^(5/6) dz                                (Eq. 20)
    Here z runs along the path from the ground (z=0) to the receiver (z=L), and
    (L-z)/L maps to (H-h)/(H-h0) over the altitude grid. So the turbulence near
    the ground gets the full weight and the turbulence near the receiver gets
    almost none. That is the uplink weighting.
'''

import numpy as np
from scipy.special import gamma, hyp1f1

# ponytail: DEBT. This analytic Dios path duplicates the beam-wave scintillation
# physics inside the shared coupled-flux MC (olb.turbulence.coupled_flux, which
# now folds pointing jitter + beam wander into the off-axis radius r=beta). This
# analytic twin still takes a bare r and knows nothing of that correction. The
# terrestrial scintillation slot will want the off-axis form WITH jitter folded
# into r -- converge on ONE implementation, do not make a third copy. See the
# memory note dios-scintillation-convergence.

# Leading constant of the Kolmogorov spectrum, Phi_n = 0.033 Cn2 kappa^(-11/3).
_KOLMOGOROV = 0.033
_GAMMA_M56 = gamma(-5.0 / 6.0)   # Gamma(-5/6) ~ -6.6865, negative by design


def _beam_and_path(hs, w0, wavelength, elevation_deg, f0, path_length_m):
    '''
    Return the beam parameters and the path integrand pieces over the grid.

    Build the receiver-plane beam parameters (Theta, Lambda), the receiver beam
    width squared W2, and the arrays A(h) and B(h).

    The geometric range L to the receiver is separate from the turbulence grid.
    For a horizontal link the beam ends at the top of the grid, so L defaults to
    (H - h0) sec(zeta). For a satellite uplink the receiver is far above the
    turbulence, so pass path_length_m = the full slant range. Then the weight
    (L - z)/L stays near 1 across the thin turbulence layer, which is the correct
    far-field uplink limit.

    Source: Dios et al., Applied Optics 43 (2004) 3866, Eqs. (15), (17), (18).
    '''
    hs = np.asarray(hs, dtype=float)
    k = 2.0 * np.pi / wavelength
    sec_z = 1.0 / np.sin(np.radians(elevation_deg))
    h0, H = hs[0], hs[-1]
    if path_length_m is None:
        L = (H - h0) * sec_z                     # horizontal link: L = the grid
    else:
        L = float(path_length_m)                 # uplink: L = slant range

    theta0 = 1.0 - L / f0                       # 1.0 for a collimated beam
    lambda0 = 2.0 * L / (k * w0 ** 2)
    denom = theta0 ** 2 + lambda0 ** 2
    theta = theta0 / denom
    lam = lambda0 / denom
    w2 = 2.0 * L / (k * lam)                     # W^2(L), receiver beam width^2

    z = (hs - h0) * sec_z                         # distance along the path [m]
    u = (L - z) / L                               # (L-z)/L, ~1 at ground
    z_over_l = z / L
    a = (lam * L / k) * u ** 2                   # A(h)
    b = (L / k) * u * (theta + (1.0 - theta) * z_over_l)   # B(h)
    return k, sec_z, L, theta, lam, w2, a, b


def on_axis_scintillation_index(hs, cn2_profile, w0, wavelength,
                                elevation_deg=90.0, f0=np.inf,
                                path_length_m=None):
    '''
    Return the on-axis Gaussian-beam scintillation index sigma2_I(0, L).

    Parameters:
        hs : numpy.ndarray
            Altitudes above the ground station [m]. Ascending. hs[0] is the
            ground. hs[-1] is the top of the turbulence path.
        cn2_profile : numpy.ndarray
            Zenith Cn2(h) on the hs grid [m^-2/3].
        w0 : float
            Beam RADIUS (1/e field) at the transmitter [m].
        wavelength : float
            Optical wavelength [m].
        elevation_deg : float
            Elevation angle above the horizon [deg]. 90 is the zenith.
        f0 : float
            Phase-front radius of curvature at the exit aperture [m]. Use
            numpy.inf for a collimated beam.
        path_length_m : float, optional
            Slant range L to the receiver [m]. Leave as None for a horizontal
            link (L = the grid length). Set it to the full slant range for a
            satellite uplink, where the receiver sits far above the turbulence.

    Returns:
        float
            sigma2_I(0, L).

    See the module docstring for the formula (Dios Eq. 16). In the plane-wave
    limit (large w0) it returns the plane-wave Rytov variance. In the
    spherical-wave limit (small w0) it returns 0.40 of that.
    '''
    cn2 = np.asarray(cn2_profile, dtype=float)
    k, sec_z, L, theta, lam, w2, a, b = _beam_and_path(hs, w0, wavelength,
                                                       elevation_deg, f0,
                                                       path_length_m)
    integrand = cn2 * (a ** (5.0 / 6.0)
                       - (a ** 2 + b ** 2) ** (5.0 / 12.0)
                       * np.cos((5.0 / 6.0) * np.arctan2(b, a)))
    prefactor = 4.0 * np.pi ** 2 * k ** 2 * _GAMMA_M56 * _KOLMOGOROV * sec_z
    return prefactor * np.trapz(integrand, np.asarray(hs, dtype=float))


def radial_scintillation_index(r, hs, cn2_profile, w0, wavelength,
                               elevation_deg=90.0, f0=np.inf,
                               path_length_m=None):
    '''
    Return the off-axis (radial) scintillation addition sigma2_Ir(r, L).

    Add this term to the on-axis index for a point that sits a radius r off the
    beam axis. At r = 0 it is 0. It grows with r.

    Parameters:
        r : float or numpy.ndarray
            Off-axis radius in the receiver plane [m].
        (the rest match on_axis_scintillation_index, including path_length_m)

    Returns:
        float or numpy.ndarray
            sigma2_Ir(r, L), broadcast over the r shape.

    See the module docstring for the formula (Dios Eq. 20). For small r it grows
    as 4.42 * sigma2_R * Lambda^(5/6) * (r/W)^2.
    '''
    cn2 = np.asarray(cn2_profile, dtype=float)
    k, sec_z, L, theta, lam, w2, a, b = _beam_and_path(hs, w0, wavelength,
                                                       elevation_deg, f0,
                                                       path_length_m)
    a_integral = np.trapz(cn2 * a ** (5.0 / 6.0), np.asarray(hs, dtype=float))
    prefactor = 4.0 * np.pi ** 2 * k ** 2 * _GAMMA_M56 * _KOLMOGOROV * sec_z
    arg = 2.0 * np.asarray(r, dtype=float) ** 2 / w2
    return prefactor * (hyp1f1(-5.0 / 6.0, 1.0, arg) - 1.0) * a_integral


def gaussian_scintillation_index(r, hs, cn2_profile, w0, wavelength,
                                 elevation_deg=90.0, f0=np.inf,
                                 path_length_m=None):
    '''
    Return the total Gaussian-beam scintillation index sigma2_I(r, L).

    Sum the on-axis and radial terms (Dios Eq. 13). Pass r = 0 for a point
    detector at the beam centre.

    Returns:
        float or numpy.ndarray
            sigma2_I(r, L) = sigma2_I(0, L) + sigma2_Ir(r, L).
    '''
    on_axis = on_axis_scintillation_index(hs, cn2_profile, w0, wavelength,
                                          elevation_deg, f0, path_length_m)
    radial = radial_scintillation_index(r, hs, cn2_profile, w0, wavelength,
                                        elevation_deg, f0, path_length_m)
    return on_axis + radial


if __name__ == '__main__':
    # Pure-physics self-check. Use plain numeric inputs. This module must not
    # import the scenario, the geometry, or the links.
    lam = 1550e-9
    k = 2.0 * np.pi / lam

    # A horizontal constant-Cn2 path (elevation 90, flat profile). The plane-,
    # spherical-, and radial-limit constants are known in closed form, so this
    # validates the absolute normalisation of Eqs. 16 and 20.
    L = 2000.0
    hs = np.linspace(0.0, L, 400)
    cn2_val = 3e-16                      # weak: keeps sigma2_R well below 1
    cn2 = np.full_like(hs, cn2_val)
    sigma2_R = 1.23 * cn2_val * k ** (7.0 / 6.0) * L ** (11.0 / 6.0)

    # Plane-wave limit: a very wide beam. On-axis index -> sigma2_R exactly.
    s_plane = on_axis_scintillation_index(hs, cn2, 5.0, lam)
    assert np.isclose(s_plane, sigma2_R, rtol=2e-2), (s_plane, sigma2_R)

    # Spherical-wave limit: a near-point source. On-axis index -> 0.40 sigma2_R.
    s_sph = on_axis_scintillation_index(hs, cn2, 5e-4, lam)
    assert np.isclose(s_sph, 0.404 * sigma2_R, rtol=3e-2), (s_sph / sigma2_R)

    # A finite beam sits between the two limits.
    w0 = 0.05
    s0 = on_axis_scintillation_index(hs, cn2, w0, lam)
    assert s_sph < s0 < s_plane or s_plane < s0 < s_sph, (s_sph, s0, s_plane)
    assert s0 > 0.0

    # Radial term: 0 on axis, grows with r, and the small-r slope matches
    # 4.42 sigma2_R Lambda^(5/6) (r/W)^2.
    _, _, _, theta, lam_b, w2, _, _ = _beam_and_path(hs, w0, lam, 90.0, np.inf,
                                                     None)
    assert radial_scintillation_index(0.0, hs, cn2, w0, lam) == 0.0
    r_small = 1e-3 * np.sqrt(w2)
    s_r = radial_scintillation_index(r_small, hs, cn2, w0, lam)
    pred = 4.42 * sigma2_R * lam_b ** (5.0 / 6.0) * (r_small ** 2 / w2)
    assert np.isclose(s_r, pred, rtol=3e-2), (s_r, pred)
    assert radial_scintillation_index(0.5 * np.sqrt(w2), hs, cn2, w0, lam) > s_r

    # Total is on-axis plus radial.
    tot = gaussian_scintillation_index(r_small, hs, cn2, w0, lam)
    assert np.isclose(tot, s0 + s_r), (tot, s0 + s_r)

    # Cross-check against the paper: Dios's geostationary uplink. lambda=0.84 um,
    # slant range L = 36e6 m, atmosphere to 20 km, HV Cn2 (their Eq. 28, v=21,
    # A=1.7e-14), collimated. Fig. 5 puts the small-waist on-axis sigma2_chi near
    # 0.02-0.03, so sigma2_I(0,L) ~ 4 sigma2_chi ~ 0.1 for a 1 cm waist.
    lam_geo = 0.84e-6
    hs_geo = np.linspace(1.0, 20e3, 800)
    v, Acn = 21.0, 1.7e-14
    cn2_geo = (0.00594 * (v / 27.0) ** 2 * (hs_geo * 1e-5) ** 10
               * np.exp(-hs_geo / 1000.0)
               + 2.7e-16 * np.exp(-hs_geo / 1500.0)
               + Acn * np.exp(-hs_geo / 100.0))
    s_geo = on_axis_scintillation_index(hs_geo, cn2_geo, 0.01, lam_geo,
                                        path_length_m=36e6)
    assert 0.03 < s_geo < 0.3, s_geo   # right order for Dios Fig. 5

    print(f"sigma2_R (plane ref)      = {sigma2_R:.5f}")
    print(f"on-axis, plane limit      = {s_plane:.5f}  (target {sigma2_R:.5f})")
    print(f"on-axis, spherical limit  = {s_sph:.5f}  (target "
          f"{0.404 * sigma2_R:.5f})")
    print(f"on-axis, w0={w0} m        = {s0:.5f}")
    print(f"radial small-r            = {s_r:.3e}  (target {pred:.3e})")
    print(f"Dios GEO, w0=1cm          = {s_geo:.4f}  "
          f"(sigma2_chi ~ {s_geo / 4:.4f}, Fig.5 ~ 0.02-0.03)")
    print("self-check passed")
