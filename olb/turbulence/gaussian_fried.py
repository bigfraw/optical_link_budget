'''
Gaussian-beam Fried parameter for a turbulent path.

This module gives the Fried parameter r0 of a collimated Gaussian beam that
propagates through turbulence. It is a closed-form model. It uses one path length
and one scalar Cn2. It does not integrate a Cn2 profile. The functions are pure.
They take numeric values or numpy arrays and return the same shape.

The plane-wave Fried parameter scales the plane-wave coherence radius. A finite
Gaussian beam has a different coherence radius. The beam parameters Theta and
Lambda set the difference. Strong turbulence changes the beam parameters. The
effective parameters Theta_e and Lambda_e hold that change.

Physics (collimated Gaussian beam, isotropic turbulence):
    r0_gauss = 2.1 * rho0_e * rho_pl
    rho_pl is the plane-wave coherence radius. rho0_e is the beam coherence radius
    divided by the plane-wave coherence radius. Source: Andrews and Phillips,
    Laser Beam Propagation through Random Media, 2nd ed. (2005),
    DOI 10.1117/3.626196, Ch. 6, Eq. (64), printed p. 194 (Gaussian-beam
    coherence) and Ch. 7, Eq. (58), printed p. 242 (effective beam parameters in
    strong turbulence). The beam parameters Theta0, Lambda0, Theta, Theta_bar
    and Lambda are Ch. 4, Eqs. (33), (44), (45) and (47), printed pp. 92 and 95.

    A radius w0 gives the beam its input size. The input Fresnel ratio is
    Lambda0 = 2 z / (k w0^2), with k = 2*pi/lambda. The beam is collimated, so
    the input curvature Theta0 = 1. The output parameters are
        Lambda = Lambda0 / (Lambda0^2 + Theta0^2)
        Theta  = Theta0  / (Theta0^2 + Lambda0^2)

NEW HOME. The beam parameters and the Rytov variance now live in the Andrews
foundation package, `olb.turbulence.andrews.beam` and
`olb.turbulence.andrews.scintillation`. Those modules are general in the input
curvature f0. The three functions below keep their names and their signatures,
but they call the new home. Use the new home for new code.
'''

import numpy as np

from .andrews.beam import beam_params
from .andrews.beam import effective_beam_params as _andrews_effective
from .andrews.scintillation import rytov_variance as _andrews_rytov_variance

# Input curvature of a collimated beam. The waist sits at the transmitter.
# Source: Andrews and Phillips, Laser Beam Propagation through Random Media,
# 2nd ed. (2005), DOI 10.1117/3.626196, Ch. 4, Eq. (33), printed p. 92:
# Theta0 = 1 - z/F0, which is 1 when F0 is infinite.
COLLIMATED_THETA0 = 1.0


def _k(wavelength):
    '''Return the optical wavenumber k = 2*pi/lambda.'''
    return 2.0 * np.pi / wavelength


def input_fresnel_ratio(z, w0, wavelength):
    '''
    Return the input-plane Fresnel ratio Lambda0 of a collimated beam.

    formula:
        Lambda0 = 2 z / (k w0^2),   k = 2*pi/lambda
    Here w0 is the beam RADIUS at the transmitter [m]. Source: Andrews and
    Phillips, 2nd ed. (2005), DOI 10.1117/3.626196, Ch. 4, Eq. (33), printed
    p. 92.
    '''
    return 2.0 * z / (_k(wavelength) * np.asarray(w0, dtype=float) ** 2)


def output_beam_params(z, w0, wavelength):
    '''
    Return the output-plane beam parameters (Theta, Lambda) of a collimated beam.

    formula:
        Lambda = Lambda0 / (Lambda0^2 + Theta0^2)
        Theta  = Theta0  / (Theta0^2 + Lambda0^2),   Theta0 = 1
    Source: Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196, Ch. 4,
    Eqs. (33) and (44), printed pp. 92 and 95.

    NEW HOME: `olb.turbulence.andrews.beam.beam_params`. That function is
    general in the input curvature f0 and it also gives Theta_bar and W.
    '''
    bp = beam_params(w0, wavelength, z)
    return bp.theta, bp.lam


def rytov_std(z, cn2, wavelength):
    '''
    Return the plane-wave Rytov standard deviation sigma_R.

    formula:
        sigma_R = ( 1.23 Cn2 k^(7/6) z^(11/6) )^0.5,   k = 2*pi/lambda
    The Rytov variance is sigma_R^2. Source: Andrews and Phillips, 2nd ed.
    (2005), DOI 10.1117/3.626196, Ch. 8, Eq. (20), printed p. 264.

    NEW HOME: `olb.turbulence.andrews.scintillation.rytov_variance`. That
    function gives the VARIANCE, and it also gives the spherical-wave and the
    Gaussian-beam forms.
    '''
    return np.sqrt(_andrews_rytov_variance(wavelength, z, cn2, wave='plane'))


def effective_beam_params(z, w0, cn2, wavelength):
    '''
    Return the strong-turbulence effective beam parameters (Theta_e, Lambda_e).

    Strong turbulence spreads the beam. The effective parameters hold that spread.

    formula:
        Theta_e  = (Theta - 0.81 sigma_R^(12/5) Lambda)
                   / (1 + 1.63 sigma_R^(12/5) Lambda)
        Lambda_e = Lambda / (1 + 1.63 sigma_R^(12/5) Lambda)
    Source: Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196, Ch. 7,
    Eq. (58), printed p. 242. Ch. 9, Eqs. (85) and (86), printed p. 349, restate
    it.

    NEW HOME: `olb.turbulence.andrews.beam.effective_beam_params`. That function
    takes a BeamParams and the Rytov variance, and it also gives the long-term
    beam radius W_LT.
    '''
    bp = beam_params(w0, wavelength, z)
    sigma2_R = _andrews_rytov_variance(wavelength, z, cn2, wave='plane')
    bp_e = _andrews_effective(bp, sigma2_R)
    return bp_e.theta, bp_e.lam


def _a_factor(theta_e):
    '''
    Return the a-factor of the beam coherence radius.

    formula:
        a = (1 - Theta_e^(8/3)) / (1 - Theta_e),        Theta_e >= 0
        a = (1 + |Theta_e|^(8/3)) / (1 - Theta_e),      Theta_e < 0
    Source: Andrews and Phillips, 2nd ed. (2005), Ch. 6.
    '''
    theta_e = np.asarray(theta_e, dtype=float)
    power = np.abs(theta_e) ** (8.0 / 3.0)
    numerator = np.where(theta_e >= 0.0, 1.0 - power, 1.0 + power)
    return numerator / (1.0 - theta_e)


def beam_coherence_ratio(z, w0, cn2, wavelength):
    '''
    Return rho0_e: the beam coherence radius over the plane-wave coherence radius.

    formula:
        rho0_e = ( 8 / (3 (a + 0.62 Lambda_e^(11/6))) )^(3/5)
    with a from the effective Theta_e. Source: Andrews and Phillips, 2nd ed.
    (2005), Ch. 6.
    '''
    theta_e, lam_e = effective_beam_params(z, w0, cn2, wavelength)
    a = _a_factor(theta_e)
    return (8.0 / (3.0 * (a + 0.62 * lam_e ** (11.0 / 6.0)))) ** (3.0 / 5.0)


def plane_wave_coherence_radius(z, cn2, wavelength):
    '''
    Return the plane-wave coherence radius rho_pl for a single path.

    formula:
        rho_pl = ( 1.46 Cn2 k^2 z )^(-3/5),   k = 2*pi/lambda
    Source: Andrews and Phillips, 2nd ed. (2005), Ch. 6.
    '''
    k = _k(wavelength)
    return (1.46 * np.asarray(cn2, dtype=float) * k ** 2
            * np.asarray(z, dtype=float)) ** (-3.0 / 5.0)


def plane_wave_fried_parameter(z, cn2, wavelength):
    '''
    Return the plane-wave Fried parameter r0 for a single path.

    formula:
        r0_pl = 2.1 * rho_pl
    with rho_pl the plane-wave coherence radius. This equals the standard
    r0 = (0.423 k^2 Cn2 L)^(-3/5). Source: Andrews and Phillips, 2nd ed. (2005),
    Ch. 6.
    '''
    return 2.1 * plane_wave_coherence_radius(z, cn2, wavelength)


# Ratio of the spherical-wave Fried parameter to the plane-wave Fried parameter
# for a constant Cn2 over the path. The spherical wave weights the path by
# (z/L)^(5/3). For a constant Cn2 the weighted path is 3/8 of the plane-wave
# path, so r0_sph = (8/3)^(3/5) * r0_pl. Source: Andrews and Phillips, 2nd ed.
# (2005), Ch. 6.
_SPHERICAL_OVER_PLANE = (8.0 / 3.0) ** (3.0 / 5.0)


def spherical_wave_fried_parameter(z, cn2, wavelength):
    '''
    Return the spherical-wave Fried parameter r0 for a single path.

    Use this parameter for a point source on a HORIZONTAL path with a constant
    Cn2. Do not use it for an uplink. An uplink weights the turbulence by
    ((L-z)/L)^(5/3), because the beam is small near the ground and large near the
    satellite. See Dios et al., Applied Optics 43 (2004) 3866, Eq. (3):
        r0_s = ( 0.42 k^2 INT Cn2(z) ((L-z)/L)^(5/3) dz )^(-3/5)

    formula (this function, horizontal path):
        r0_sph = (8/3)^(3/5) * r0_pl  ~ 1.80 * r0_pl
    Source: Andrews and Phillips, 2nd ed. (2005), Ch. 6.
    '''
    return _SPHERICAL_OVER_PLANE * plane_wave_fried_parameter(z, cn2, wavelength)


def gaussian_fried_parameter(z, w0, cn2, wavelength):
    '''
    Return the Fried parameter r0 of a collimated Gaussian beam.

    Parameters:
        z : float or numpy.ndarray
            Path length [m].
        w0 : float or numpy.ndarray
            Beam RADIUS at the transmitter [m].
        cn2 : float or numpy.ndarray
            Refractive-index structure constant Cn2 [m^-2/3].
        wavelength : float
            Optical wavelength [m].

    Returns:
        float or numpy.ndarray
            r0_gauss [m].

    formula:
        r0_gauss = 2.1 * rho0_e * rho_pl
    Source: Andrews and Phillips, 2nd ed. (2005), Ch. 6.
    '''
    rho0_e = beam_coherence_ratio(z, w0, cn2, wavelength)
    rho_pl = plane_wave_coherence_radius(z, cn2, wavelength)
    return 2.1 * rho0_e * rho_pl


def gaussian_fried_parameter_profile(hs, cn2_profile, w0, wavelength,
                                     path='uplink', elevation_deg=90.0,
                                     f0=np.inf, path_length_m=None):
    '''
    Return the Gaussian-beam Fried parameter over a nonuniform Cn2 path.

    Integrate the real Cn2 profile with the wave-structure-function path weight.
    Use this parameter for an uplink, a downlink, or a horizontal terrestrial
    link. It is more accurate than the constant-Cn2 gaussian_fried_parameter.
    Feed it to the Dios uplink model in place of the spherical-wave coherence
    r0_s.

    ASSUMPTION (weak turbulence): the model assumes a Gaussian beam in the WEAK
    turbulence regime with a Kolmogorov spectrum (no inner scale, no outer
    scale). It holds for weak-to-moderate turbulence, about a Rytov variance
    sigma_R^2 below 1 (Dios reports good agreement to sigma_chi^2 ~ 0.6). Above
    that the coherence saturates and this parameter reads too small a loss.

    IMPORTANT: the beam parameters use the FREE-SPACE (diffractive) Theta and
    Lambda, NOT the strong-turbulence effective Theta_e and Lambda_e. So the
    model does NOT capture the turbulence-driven beam spread that shortens the
    coherence length in strong turbulence. To extend it past the weak regime,
    substitute effective_beam_params for the theta/lambda block below. This is
    a deliberate deferral, kept to match the Dios weak-turbulence regime.

    The geometric range L to the receiver is SEPARATE from the Cn2 grid. For a
    terrestrial or a downlink path the beam ends at the far edge of the grid, so
    L defaults to the grid length. For a satellite uplink the receiver sits far
    above the turbulence, so pass path_length_m = the full slant range. Then the
    beam parameters and the weight use that range, and the near-ground turbulence
    keeps the full weight, which is the correct far-field uplink limit.

    Parameters:
        hs : numpy.ndarray
            For "uplink"/"downlink": altitudes above the ground station [m],
            ascending, hs[0] the ground and hs[-1] the top of the path. For
            "terrestrial": distance along the horizontal path from the near
            terminal [m], ascending, hs[0] = 0.
        cn2_profile : numpy.ndarray
            Cn2 on the hs grid [m^-2/3]. Zenith Cn2(h) for a slant path.
        w0 : float
            Beam RADIUS (1/e field) at the transmitter [m].
        wavelength : float
            Optical wavelength [m].
        path : str
            "uplink", "downlink", or "terrestrial". It sets the transmitter end
            and the weight direction.
        elevation_deg : float
            Elevation angle above the horizon [deg]. 90 is the zenith. Ignored
            for "terrestrial".
        f0 : float
            Phase-front radius of curvature at the exit aperture [m]. Use
            numpy.inf for a collimated beam.
        path_length_m : float, optional
            Geometric range L to the receiver [m]. Leave as None to use the grid
            length. Set it to the full slant range for a satellite uplink.

    Returns:
        float
            r0 = 2.1 * rho0 [m].

    formula:
        z      = distance along the path from the transmitter [m]
        xi     = (L - z) / L         (1 at the transmitter, 0 at the receiver)
        Theta0 = 1 - L / f0,   Lambda0 = 2 L / (k w0^2),   k = 2*pi/lambda
        Theta  = Theta0 / (Theta0^2 + Lambda0^2),   Theta_bar = 1 - Theta
        Lambda = Lambda0 / (Theta0^2 + Lambda0^2)
        mu1 = INT Cn2 (Theta + Theta_bar xi)^(5/3) dh
        mu2 = INT Cn2 xi^(5/3) dh
        rho0 = ( 1.46 k^2 sec(zeta) (mu1 + 0.62 Lambda^(11/6) mu2) )^(-3/5)
        r0   = 2.1 * rho0
    In the plane-wave limit (Theta -> 1) the weight is 1 and r0 is the plane-wave
    Fried parameter. In the spherical-wave limit (Theta -> 0, Lambda -> 0) the
    weight is xi^(5/3): the turbulence at the transmitter carries the full
    weight, which matches Dios Eq. (3). Source: Andrews and Phillips, 2nd ed.
    (2005), DOI 10.1117/3.626196, Ch. 12, Eq. (23), printed p. 492 (the slant
    airmass), Ch. 6, Eq. (115), printed p. 209 (the path weight), and Ch. 4,
    Eqs. (33), (44), (45) and (47), printed pp. 92 and 95 (the beam parameters).

    PLANE OF REFERENCE: the weight xi = (L - z)/L is TRANSMITTER-referred, per
    Dios et al., Appl. Opt. 43 (2004) 3866, DOI 10.1364/AO.43.003866, Eq. (3).
    The book Ch. 6, Eq. (115) writes the mirror weight z/L, which is
    RECEIVER-referred. The difference is the plane of reference, not a fault.
    '''
    hs = np.asarray(hs, dtype=float)
    cn2 = np.asarray(cn2_profile, dtype=float)
    k = _k(wavelength)
    h0, H = hs[0], hs[-1]

    if path == 'uplink':
        sec_z = 1.0 / np.sin(np.radians(elevation_deg))   # sec(zenith)
        z = (hs - h0) * sec_z                             # from ground tx
    elif path == 'downlink':
        sec_z = 1.0 / np.sin(np.radians(elevation_deg))
        z = (H - hs) * sec_z                              # from space tx
    elif path == 'terrestrial':
        sec_z = 1.0                                       # horizontal, no airmass
        z = hs - h0                                       # from near tx
    else:
        raise ValueError("path must be 'uplink', 'downlink', or 'terrestrial'")

    turb_len = (H - h0) * sec_z                            # length of the grid
    L = turb_len if path_length_m is None else float(path_length_m)

    # ponytail: free-space (weak-turbulence) beam parameters. For strong
    # turbulence swap in effective_beam_params (Theta_e, Lambda_e). Deferred to
    # match the Dios weak-turbulence regime. See the ASSUMPTION note above.
    theta0 = 1.0 - L / f0                                  # 1.0 for collimated
    lambda0 = 2.0 * L / (k * w0 ** 2)
    denom = theta0 ** 2 + lambda0 ** 2
    theta = theta0 / denom
    theta_bar = 1.0 - theta
    lam = lambda0 / denom

    xi = (L - z) / L                                      # 1 at tx, 0 at rx

    mu1 = np.trapz(cn2 * (theta + theta_bar * xi) ** (5.0 / 3.0), hs)
    mu2 = np.trapz(cn2 * xi ** (5.0 / 3.0), hs)

    rho0 = (1.46 * k ** 2 * sec_z
            * (mu1 + 0.62 * lam ** (11.0 / 6.0) * mu2)) ** (-3.0 / 5.0)
    return 2.1 * rho0


if __name__ == '__main__':
    # Pure-physics self-check. Use plain numeric inputs. This module must not
    # import the scenario, the geometry, or the links.
    lam = 1550e-9
    z = 750.0
    w0 = 0.05          # beam radius [m]
    cn2 = 1e-14

    r0 = gaussian_fried_parameter(z, w0, cn2, lam)
    assert r0 > 0.0, r0

    # Stronger turbulence gives a smaller Fried parameter.
    r0_weak = gaussian_fried_parameter(z, w0, 1e-15, lam)
    r0_strong = gaussian_fried_parameter(z, w0, 1e-13, lam)
    assert r0_strong < r0_weak, (r0_strong, r0_weak)

    # A longer path gives a smaller Fried parameter.
    r0_far = gaussian_fried_parameter(2400.0, w0, cn2, lam)
    assert r0_far < r0, (r0_far, r0)

    # The array path gives the same values as the scalar path.
    zs = np.array([750.0, 2400.0])
    r0s = gaussian_fried_parameter(zs, w0, cn2, lam)
    assert np.isclose(r0s[0], r0), (r0s[0], r0)
    assert np.isclose(r0s[1], r0_far), (r0s[1], r0_far)

    # Profile version: check the plane- and spherical-wave limits over a
    # constant-Cn2 slant path (same test as the thesis fried_parameter).
    hs = np.linspace(0.0, 1000.0, 200)
    cn2_flat = np.full_like(hs, 1e-15)
    L = hs[-1] - hs[0]
    r0_plane_ref = (0.423 * (2 * np.pi / lam) ** 2 * 1e-15 * L) ** (-3.0 / 5.0)
    r0_sph_ref = (0.423 * (2 * np.pi / lam) ** 2 * 1e-15 * L
                  * (3.0 / 8.0)) ** (-3.0 / 5.0)
    for p in ('uplink', 'downlink'):
        # Wide beam -> plane-wave limit.
        r0_pl = gaussian_fried_parameter_profile(hs, cn2_flat, 1e3, lam,
                                                 path=p, elevation_deg=90.0)
        assert np.isclose(r0_pl, r0_plane_ref, rtol=1e-2), (p, r0_pl)
        # Point source -> spherical-wave limit.
        r0_sph = gaussian_fried_parameter_profile(hs, cn2_flat, 1e-4, lam,
                                                  path=p, elevation_deg=90.0)
        assert np.isclose(r0_sph, r0_sph_ref, rtol=1e-2), (p, r0_sph)
    # Spherical r0 is ~1.8x the plane-wave r0 (the (8/3)^(3/5) factor).
    assert r0_sph_ref > r0_plane_ref

    # Weighting direction (uplink). Put the same total turbulence near the ground
    # or near the top. The uplink is dominated by near-ground turbulence, so the
    # ground-heavy profile must give the SMALLER r0. This fails on the old
    # (1 - xi) weight, which inverts the direction.
    hs_w = np.linspace(1.0, 20e3, 400)
    ground_heavy = 1e-15 * np.exp(-hs_w / 500.0)
    top_heavy = 1e-15 * np.exp(-(hs_w[-1] - hs_w) / 500.0)
    ground_heavy *= np.trapz(top_heavy, hs_w) / np.trapz(ground_heavy, hs_w)
    r0_ground = gaussian_fried_parameter_profile(hs_w, ground_heavy, 5e-4, lam,
                                                 path='uplink')
    r0_top = gaussian_fried_parameter_profile(hs_w, top_heavy, 5e-4, lam,
                                               path='uplink')
    assert r0_ground < r0_top, (r0_ground, r0_top)

    # Far-field L-separation. A tiny waist to a satellite (L=36e6 m) over the
    # 20 km layer: xi ~ 1 across the layer, so the spherical weight collapses to
    # the plane-wave weight and r0 approaches the plane-wave r0 of that layer.
    hv = 1e-15 * np.exp(-hs_w / 1000.0)
    r0_geo = gaussian_fried_parameter_profile(hs_w, hv, 0.01, lam,
                                              path='uplink', path_length_m=36e6)
    r0_geo_plane = 2.1 * (1.46 * (2 * np.pi / lam) ** 2
                          * np.trapz(hv, hs_w)) ** (-3.0 / 5.0)
    assert np.isclose(r0_geo, r0_geo_plane, rtol=5e-2), (r0_geo, r0_geo_plane)

    # Terrestrial path: hs is horizontal distance, no airmass. Constant Cn2 over
    # 2 km, wide beam -> plane-wave r0 of that horizontal path.
    d = np.linspace(0.0, 2000.0, 200)
    cn2_t = np.full_like(d, 1e-15)
    r0_terr = gaussian_fried_parameter_profile(d, cn2_t, 1e3, lam,
                                               path='terrestrial')
    r0_terr_ref = (0.423 * (2 * np.pi / lam) ** 2 * 1e-15 * 2000.0) ** (-3.0 / 5.0)
    assert np.isclose(r0_terr, r0_terr_ref, rtol=1e-2), (r0_terr, r0_terr_ref)

    print(f"r0_gauss z=750 m  cn2=1e-14 = {r0 * 100:.2f} cm")
    print(f"r0_gauss z=750 m  cn2=1e-15 = {r0_weak * 100:.2f} cm")
    print(f"r0_gauss z=750 m  cn2=1e-13 = {r0_strong * 100:.2f} cm")
    print(f"r0_gauss z=2400 m cn2=1e-14 = {r0_far * 100:.2f} cm")
    print(f"profile plane-wave limit  r0 = {r0_plane_ref * 100:.2f} cm")
    print(f"profile spherical limit   r0 = {r0_sph_ref * 100:.2f} cm")
    print(f"uplink ground-heavy r0    = {r0_ground * 100:.2f} cm  "
          f"< top-heavy {r0_top * 100:.2f} cm")
    print(f"GEO uplink w0=1cm r0      = {r0_geo * 100:.2f} cm  "
          f"(plane-layer ref {r0_geo_plane * 100:.2f} cm)")
    print(f"terrestrial 2 km r0       = {r0_terr * 100:.2f} cm")
    print("self-check passed")
