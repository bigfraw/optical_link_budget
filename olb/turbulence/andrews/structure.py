'''
Wave structure function, coherence radius and angle of arrival (Andrews).

This module gives the second-order coherence quantities of a wave that crossed a
turbulent path:

- `wave_structure_function` gives D(r, L), the wave structure function. It is the
  exponent of the complex degree of coherence, DOC = exp(-D/2).
- `coherence_radius` gives rho_0, the separation at which D = 2.
- `fried_parameter` gives the atmospheric coherence width r_0 = 2.1 rho_0.
- `angle_of_arrival_variance` gives the tilt variance of the wavefront across a
  collecting lens.

Source of every equation:
    L. C. Andrews and R. L. Phillips, "Laser Beam Propagation through Random
    Media", 2nd ed., SPIE Press (2005). DOI: 10.1117/3.626196
Chapter 6, Secs. 6.4 and 6.5, printed pp. 192 to 201, with the closed forms
collated in Appendix III, Tables I to VI, printed pp. 765 to 768. Each function
names its section, its equation number, and its printed page.

THREE WAVE TYPES, THREE SPECTRA. Every function takes `wave` ("plane",
"spherical" or "gaussian") and `spectrum` ("kolmogorov", "von_karman" or
"modified"). The Kolmogorov branch has no inner scale and no outer scale. The
other two branches need `l0`, and they take an optional `L0`.

PLANE OF REFERENCE. This module takes ONE path length L and ONE scalar Cn2, so
it makes no path integral and it picks no reference plane. The book results are
printed for a horizontal path with a constant Cn2 (Appendix III, header text,
printed p. 765). A caller that integrates a Cn2 profile must choose the
reference plane itself. See Conflict C-02 in docs/andrews-crosscheck.md: the olb
uplink weight is transmitter-referred by design, and the book Ch. 6, Eq. (115)
weight is receiver-referred. Do not flip either one.

This module holds physics only. It imports numpy and sibling andrews modules. It
returns no decibels.
'''

import numpy as np

from ...assumptions import Constraint, module_assumptions
from .beam import wavenumber

_WAVES = ('plane', 'spherical', 'gaussian')
_SPECTRA = ('kolmogorov', 'von_karman', 'modified')

# Fried's atmospheric coherence width over the spatial coherence radius. Source:
# Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196, Ch. 6, text below
# Eq. (64), printed p. 194: "it is customary to define spatial coherence by the
# related atmospheric coherence width r0 = 2.1 rho_0". Appendix III, Table IV
# footnote, printed p. 767, repeats it. Ch. 12, text at Eq. (23), printed p. 492,
# and Ch. 14, Eq. (89), printed p. 635, repeat it again.
FRIED_OVER_RHO0 = 2.1

# Inner-scale wavenumber ratios of Appendix III, Table III footer, printed
# p. 766: Q_m = L km^2/k = 35.04 L/(k l0^2) for the von Karman spectrum, and
# Q_l = L kl^2/k = 10.89 L/(k l0^2) for the modified atmospheric spectrum.
_QM_CONSTANT = 35.04
_QL_CONSTANT = 10.89

# Outer-scale wavenumber k0 = 2*pi/L0. Source: Ch. 3, Eq. (20), printed p. 68.
_K0_CONSTANT = 2.0 * np.pi


# ---------------------------------------------------------------------------
# Function-owned assumptions (see olb/assumptions.py and
# docs/assumptions-refactor-plan.md).
#
# Two statements hold for EVERY public function in the file, so they are module
# defaults: the one-path homogeneity (PLANE OF REFERENCE), and the C-02
# reference-plane conflict (the olb uplink weight is transmitter-referred, the
# book weight is receiver-referred). The spectrum is NOT a module default,
# because a function takes "kolmogorov", "von_karman" or "modified" by argument.
# ---------------------------------------------------------------------------
PATH_HOMOGENEITY = Constraint(
    "path-homogeneity",
    "One path length L and one scalar Cn2. The book prints these forms for a "
    "horizontal path with a constant Cn2. No profile integral, no reference "
    "plane.",
    "10.1117/3.626196", "Appendix III, header text, printed p. 765")

PATH_WEIGHT_CONFLICT = Constraint(
    "conflict",
    "A caller that builds a slant integral must pick the reference plane. The "
    "olb uplink weight is transmitter-referred; the book weight is "
    "receiver-referred. Do not flip either one (Conflict C-02).",
    "10.1117/3.626196", "Ch. 6, Eq. (115), printed p. 199")

assumes = module_assumptions(
    constraints=(PATH_HOMOGENEITY, PATH_WEIGHT_CONFLICT))

MODIFIED_GAUSSIAN_NOT_BUILT = Constraint(
    "not-built",
    "The Gaussian wave structure function on the modified spectrum is not "
    "built; wave='gaussian' with spectrum='modified' raises.",
    "10.1117/3.626196", "Appendix III, Table III, printed p. 766")

OUTER_SCALE_INFINITE = Constraint(
    "spectrum",
    "The outer scale is infinite: the book prints the coherence-radius rows for "
    "k0 = 0 only. A finite L0 raises.",
    "10.1117/3.626196", "Appendix III, Tables IV to VI, printed pp. 767-768")

TILT_CONVENTION_G = Constraint(
    "tilt-convention",
    "The returned tilt is the Andrews gradient tilt (G-tilt), what a centroid "
    "tracker measures. It is NOT the Noll Zernike tilt. A caller that mixes "
    "conventions must say which tilt it means (Conflict C-04).",
    "10.1117/3.626196", "Ch. 6, Eq. (84), printed p. 201")


def _fresnel_zone_check(args, result):
    '''Return a reason when the Fresnel zone is not small against the lens.

    Eq. (83) needs sqrt(L/k) << D (Ch. 6, text below Eq. (83), printed p. 200).
    The check reports the worst case (the largest Fresnel scale, the smallest
    lens) and flags the boundary where the Fresnel scale reaches the lens size.
    No warning here.
    '''
    L = float(np.max(np.asarray(args['z'], dtype=float)))
    D = float(np.min(np.asarray(args['D'], dtype=float)))
    k = float(np.min(np.asarray(wavenumber(args['wavelength']), dtype=float)))
    fresnel = float(np.sqrt(L / k))
    if fresnel >= D:
        return (f"sqrt(L/k) = {fresnel:.4f} m is not << D = {D:.4f} m; the "
                f"small-Fresnel-zone condition for Eq. (83) fails.")
    return None


FRESNEL_ZONE_CONSTRAINT = Constraint(
    "field-region",
    "The Fresnel zone is small against the lens: sqrt(L/k) << D.",
    "10.1117/3.626196", "Ch. 6, text below Eq. (83), printed p. 200",
    check=_fresnel_zone_check)


def _check(wave, spectrum):
    '''Refuse a wave type or a spectrum that this module does not know.'''
    if wave not in _WAVES:
        raise ValueError(f'wave must be one of {_WAVES}, not {wave!r}')
    if spectrum not in _SPECTRA:
        raise ValueError(f'spectrum must be one of {_SPECTRA}, '
                         f'not {spectrum!r}')


def _scales(spectrum, l0, L0):
    '''Return (l0, k0) after the spectrum-dependent checks.'''
    if spectrum == 'kolmogorov':
        if l0 is not None or L0 is not None:
            raise ValueError('spectrum="kolmogorov" carries no inner scale and '
                             'no outer scale; use "von_karman" or "modified"')
        return None, 0.0
    if l0 is None:
        raise ValueError(f'spectrum={spectrum!r} needs l0')
    k0 = 0.0 if L0 is None else _K0_CONSTANT / float(L0)
    return float(l0), k0


def _need_beam(beam, wave):
    '''Return the beam parameters, or refuse if the caller gave none.'''
    if wave == 'gaussian' and beam is None:
        raise ValueError('wave="gaussian" needs beam=BeamParams(...)')
    return beam


def a_factor(theta):
    '''
    Return the a-factor of the Gaussian-beam wave structure function.

    formula:
        a = (1 - Theta^(8/3)) / (1 - Theta),     Theta >= 0
        a = (1 + |Theta|^(8/3)) / (1 - Theta),   Theta < 0
    Source: Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196, Ch. 6,
    Eq. (55), printed p. 192. Appendix III, Table III footer, printed p. 766,
    repeats it.

    The plane-wave limit Theta -> 1 gives a = 8/3. The spherical-wave limit
    Theta = 0 gives a = 1. Source: Ch. 6, text below Eq. (55), printed p. 192.
    '''
    theta = np.asarray(theta, dtype=float)
    power = np.abs(theta) ** (8.0 / 3.0)
    numerator = np.where(theta >= 0.0, 1.0 - power, 1.0 + power)
    # The Theta = 1 limit of (1 - Theta^(8/3))/(1 - Theta) is 8/3.
    near_one = np.abs(1.0 - theta) < 1e-8
    safe = np.where(near_one, 1.0, 1.0 - theta)
    return np.where(near_one, 8.0 / 3.0, numerator / safe)


def _theta_difference(theta, term):
    '''
    Return [term(1) - Theta^3 term(Theta^2)] / (1 - Theta), with the Theta = 1
    limit taken analytically.

    The Gaussian-beam rows of Appendix III, Tables I to III, printed pp. 765 and
    766, all carry this shape. `term` must be a function u -> (base + u*x)^(-p)
    of the SCALED separation only, so the limit is
        3 term(1) + 2 dterm/du(1).
    '''
    theta = np.asarray(theta, dtype=float)
    near_one = np.abs(1.0 - theta) < 1e-6
    safe_theta = np.where(near_one, 0.0, theta)
    ratio = (term(1.0) - safe_theta ** 3 * term(safe_theta ** 2)) \
        / np.where(near_one, 1.0, 1.0 - safe_theta)
    if not np.any(near_one):
        return ratio
    # L'Hopital at Theta = 1. A central difference on u is enough, because the
    # term is smooth.
    du = 1e-6
    dterm = (term(1.0 + du) - term(1.0 - du)) / (2.0 * du)
    return np.where(near_one, 3.0 * term(1.0) + 2.0 * dterm, ratio)


@assumes(MODIFIED_GAUSSIAN_NOT_BUILT)
def wave_structure_function(rho, wavelength, z, cn2, *, wave='plane',
                            spectrum='kolmogorov', l0=None, L0=None, beam=None):
    '''
    Return the wave structure function D(rho, L).

    D is the exponent of the complex degree of coherence: DOC = exp(-D/2). The
    coherence radius is the separation at which D = 2 (Ch. 6, text at Eq. (56),
    printed p. 193).

    Parameters:
        rho : float or numpy.ndarray
            Separation of the two observation points [m].
        wavelength : float or numpy.ndarray
            Optical wavelength [m].
        z : float or numpy.ndarray
            Path length L [m].
        cn2 : float or numpy.ndarray
            Refractive-index structure constant [m^-2/3], constant over the path.
        wave : str
            "plane", "spherical" or "gaussian".
        spectrum : str
            "kolmogorov", "von_karman" or "modified".
        l0 : float, optional
            Inner scale [m]. Required for "von_karman" and "modified".
        L0 : float, optional
            Outer scale [m]. None gives an infinite outer scale.
        beam : BeamParams, optional
            The beam parameters at the receiver. Required for "gaussian".

    Returns:
        float or numpy.ndarray
            D(rho, L), dimensionless.

    formula (Kolmogorov):
        plane      D = 2.914 Cn2 k^2 L rho^(5/3)
        spherical  D = 1.093 Cn2 k^2 L rho^(5/3)
        gaussian   D = 1.093 Cn2 k^(7/6) L^(11/6)
                       [ a (k rho^2/L)^(5/6) + 0.618 Lambda^(11/6) (k rho^2/L) ]
    Source: Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196,
    Appendix III, Table I (plane, printed p. 765), Table II (spherical, printed
    p. 765) and Table III (Gaussian, printed p. 766). Ch. 6, Eqs. (63), (70) and
    (74), printed pp. 194, 195 and 196, print the same three results.

    formula (von Karman, with Q_m = 35.04 L/(k l0^2) and k0 = 2*pi/L0):
        plane      D = 3.280 Cn2 k^2 L l0^(-1/3) rho^2
                       [ (1 + 2.033 rho^2/l0^2)^(-1/6) - 0.715 (k0 l0)^(1/3) ]
        spherical  D = 1.093 Cn2 k^2 L l0^(-1/3) rho^2
                       [ (1 + rho^2/l0^2)^(-1/6) - 0.715 (k0 l0)^(1/3) ]
        gaussian   D = 1.093 Cn2 k^2 L l0^(-1/3) rho^2
                       { Lambda^2 (1 + 0.52 Lambda Q_m)^(-1/6)
                         - 0.715 (1 + Theta + Theta^2 + Lambda^2)(k0 l0)^(1/3)
                         + [ (1 + 0.11 Lambda Q_m + rho^2/l0^2)^(-1/6)
                             - Theta^3 (1 + 0.11 Lambda Q_m
                                        + Theta^2 rho^2/l0^2)^(-1/6) ]
                           / (1 - Theta) }
    Source: Appendix III, Tables I to III, printed pp. 765 and 766. Ch. 6,
    Eqs. (62), (69), (75) and (76), printed pp. 194, 195, 196 and 197, print the
    same three results.

    formula (modified atmospheric, with Q_l = 10.89 L/(k l0^2)):
        plane      D = 2.700 Cn2 k^2 L l0^(-1/3) rho^2
                       [ (1 + 0.632 rho^2/l0^2)^(-1/6)
                         + 0.438 (1 + 0.442 rho^2/l0^2)^(-2/3)
                         - 0.056 (1 + 0.376 rho^2/l0^2)^(-3/4)
                         - 0.868 (k0 l0)^(1/3) ]
        spherical  D = 0.900 Cn2 k^2 L l0^(-1/3) rho^2
                       [ (1 + 0.311 rho^2/l0^2)^(-1/6)
                         + 0.438 (1 + 0.183 rho^2/l0^2)^(-2/3)
                         - 0.056 (1 + 0.149 rho^2/l0^2)^(-3/4)
                         - 0.868 (k0 l0)^(1/3) ]
        gaussian   the same three roll-off groups, each divided by (1 - Theta)
                   and each carrying the Theta^3 second copy, plus the three
                   Lambda terms and the outer-scale term
    Source: Appendix III, Tables I to III, printed pp. 765 and 766.

    READING NOTE on the Gaussian row of the modified spectrum. That row is the
    longest cell of Table III. Its three Lambda-only terms all go to zero in the
    plane-wave limit (Lambda = 0), in the spherical-wave limit (Lambda = 0), and
    in the zero-inner-scale limit (Q_l -> infinity). So those three limits, which
    this module checks, cannot confirm the three Lambda-only terms. The other
    terms are confirmed by all three limits.
    '''
    _check(wave, spectrum)
    l0_val, k0 = _scales(spectrum, l0, L0)
    bm = _need_beam(beam, wave)

    rho = np.asarray(rho, dtype=float)
    z = np.asarray(z, dtype=float)
    cn2 = np.asarray(cn2, dtype=float)
    k = wavenumber(wavelength)

    if spectrum == 'kolmogorov':
        if wave == 'plane':
            return 2.914 * cn2 * k ** 2 * z * rho ** (5.0 / 3.0)
        if wave == 'spherical':
            return 1.093 * cn2 * k ** 2 * z * rho ** (5.0 / 3.0)
        u = k * rho ** 2 / z
        return (1.093 * cn2 * k ** (7.0 / 6.0) * z ** (11.0 / 6.0)
                * (a_factor(bm.theta) * u ** (5.0 / 6.0)
                   + 0.618 * bm.lam ** (11.0 / 6.0) * u))

    scale = cn2 * k ** 2 * z * l0_val ** (-1.0 / 3.0) * rho ** 2
    x = (rho / l0_val) ** 2
    outer = (k0 * l0_val) ** (1.0 / 3.0)

    if spectrum == 'von_karman':
        qm = _QM_CONSTANT * z / (k * l0_val ** 2)
        if wave == 'plane':
            return 3.280 * scale * ((1.0 + 2.033 * x) ** (-1.0 / 6.0)
                                    - 0.715 * outer)
        if wave == 'spherical':
            return 1.093 * scale * ((1.0 + x) ** (-1.0 / 6.0) - 0.715 * outer)
        base = 1.0 + 0.11 * bm.lam * qm
        roll = _theta_difference(bm.theta,
                                 lambda u: (base + u * x) ** (-1.0 / 6.0))
        weight = 1.0 + bm.theta + bm.theta ** 2 + bm.lam ** 2
        return 1.093 * scale * (
            bm.lam ** 2 * (1.0 + 0.52 * bm.lam * qm) ** (-1.0 / 6.0)
            - 0.715 * weight * outer + roll)

    ql = _QL_CONSTANT * z / (k * l0_val ** 2)
    if wave == 'plane':
        return 2.700 * scale * ((1.0 + 0.632 * x) ** (-1.0 / 6.0)
                                + 0.438 * (1.0 + 0.442 * x) ** (-2.0 / 3.0)
                                - 0.056 * (1.0 + 0.376 * x) ** (-3.0 / 4.0)
                                - 0.868 * outer)
    if wave == 'spherical':
        return 0.900 * scale * ((1.0 + 0.311 * x) ** (-1.0 / 6.0)
                                + 0.438 * (1.0 + 0.183 * x) ** (-2.0 / 3.0)
                                - 0.056 * (1.0 + 0.149 * x) ** (-3.0 / 4.0)
                                - 0.868 * outer)

    raise NotImplementedError(
        'wave="gaussian" with spectrum="modified" is NOT built. Andrews and '
        'Phillips, 2nd ed. (2005), DOI 10.1117/3.626196, print the row in '
        'Appendix III, Table III, printed p. 766, but the two Lambda-only bump '
        'terms cannot be read unambiguously from the source PDF. As read, they '
        'carry the numerators 0.438 (Lambda Q_l)^(1/6) and 0.056 (Lambda '
        'Q_l)^(1/6), which fall only as Lambda^(1/6). That breaks the '
        'plane-wave reduction by 2.3 %, and Ch. 6, text below Eq. (77), '
        'printed p. 197, states the Gaussian row MUST reduce to the plane row '
        'exactly. The owner must read Table III, printed p. 766, and give the '
        'two numerators. Do not guess them. Use spectrum="von_karman" for a '
        'two-scale Gaussian wave structure function.')


@assumes(OUTER_SCALE_INFINITE)
def coherence_radius(wavelength, z, cn2, *, wave='plane',
                     spectrum='kolmogorov', l0=None, L0=None, beam=None,
                     branch='auto'):
    '''
    Return the spatial coherence radius rho_0 [m].

    rho_0 is the separation at which the wave structure function reaches 2, so
    the complex degree of coherence falls to 1/e (Ch. 6, text at Eq. (56),
    printed p. 193).

    Parameters:
        wavelength, z, cn2, wave, spectrum, l0, beam
            The same as `wave_structure_function`.
        L0 : None
            The book prints Tables IV to VI for k0 = 0 only, so this function
            takes an INFINITE outer scale. A value raises ValueError.
        branch : str
            "inertial" gives the l0 << rho_0 << L0 row. "inner" gives the
            rho_0 << l0 row. "auto" computes the inertial row first and takes
            the inner row when that value falls below l0.

    Returns:
        float or numpy.ndarray
            rho_0 [m].

    formula (l0 << rho_0 << L0, every spectrum):
        plane      rho_pl = (1.46 Cn2 k^2 L)^(-3/5)
        spherical  rho_sp = (0.55 Cn2 k^2 L)^(-3/5)
        gaussian   rho_0  = [ 8 / (3 (a + 0.618 Lambda^(11/6))) ]^(3/5)
                            (1.46 Cn2 k^2 L)^(-3/5)
    formula (rho_0 << l0):
        plane      (C Cn2 k^2 L l0^(-1/3))^(-1/2),  C = 1.64 von Karman,
                                                    C = 1.87 modified
        spherical  (C Cn2 k^2 L l0^(-1/3))^(-1/2),  C = 0.55 von Karman,
                                                    C = 0.62 modified
        gaussian   [ 3 / (1 + Theta + Theta^2 + Lambda^2) ]^(1/2)
                   (C Cn2 k^2 L l0^(-1/3))^(-1/2),  same C as the plane row
    Source: Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196,
    Appendix III, Table IV (plane, printed p. 767), Table V (spherical, printed
    p. 767) and Table VI (Gaussian, printed p. 768). Ch. 6, Eqs. (64), (71) and
    (78), printed pp. 194, 196 and 198, print the same rows.

    The Gaussian rows reduce to the plane rows at Theta = 1, Lambda = 0 and to
    the spherical rows at Theta = Lambda = 0 (Ch. 6, text below Eq. (78),
    printed p. 198). The module self-check measures both reductions.
    '''
    _check(wave, spectrum)
    if L0 is not None:
        raise ValueError(
            'coherence_radius takes an infinite outer scale. Andrews and '
            'Phillips, 2nd ed. (2005), DOI 10.1117/3.626196, print Tables IV '
            'to VI (printed pp. 767 and 768) for k0 = 0 only.')
    l0_val, _ = _scales(spectrum, l0, None)
    bm = _need_beam(beam, wave)

    z = np.asarray(z, dtype=float)
    cn2 = np.asarray(cn2, dtype=float)
    k = wavenumber(wavelength)
    ck2l = cn2 * k ** 2 * z

    inertial_c = 1.46 if wave != 'spherical' else 0.55
    inertial = (inertial_c * ck2l) ** (-3.0 / 5.0)
    if wave == 'gaussian':
        a = a_factor(bm.theta)
        inertial = inertial * (8.0 / (3.0 * (a + 0.618
                                             * bm.lam ** (11.0 / 6.0)))
                               ) ** (3.0 / 5.0)

    if spectrum == 'kolmogorov':
        if branch == 'inner':
            raise ValueError('spectrum="kolmogorov" has no inner-scale branch')
        return inertial

    if wave == 'spherical':
        inner_c = 0.55 if spectrum == 'von_karman' else 0.62
    else:
        inner_c = 1.64 if spectrum == 'von_karman' else 1.87
    inner = (inner_c * ck2l * l0_val ** (-1.0 / 3.0)) ** (-0.5)
    if wave == 'gaussian':
        weight = 1.0 + bm.theta + bm.theta ** 2 + bm.lam ** 2
        inner = inner * (3.0 / weight) ** 0.5

    if branch == 'inertial':
        return inertial
    if branch == 'inner':
        return inner
    if branch != 'auto':
        raise ValueError(f'branch must be "auto", "inertial" or "inner", '
                         f'not {branch!r}')
    return np.where(inertial < l0_val, inner, inertial)


def fried_parameter(rho0):
    '''
    Return Fried's atmospheric coherence width r_0 = 2.1 rho_0 [m].

    Parameters:
        rho0 : float or numpy.ndarray
            The spatial coherence radius [m], from `coherence_radius`.

    formula:
        r_0 = 2.1 rho_0
    Source: Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196, Ch. 6,
    text below Eq. (64), printed p. 194. Appendix III, Table IV footnote,
    printed p. 767, Ch. 12, text at Eq. (23), printed p. 492, and Ch. 14,
    Eq. (89), printed p. 635, repeat it.

    NOTE on the numeric constant. The chain r_0 = 2.1 (1.46 Cn2 k^2 L)^(-3/5)
    equals (0.4240 Cn2 k^2 L)^(-3/5). The book Ch. 12, Eq. (23), printed p. 492,
    prints the rounded form (0.42 sec(zeta) k^2 INT Cn2 dh)^(-3/5). The classic
    Fried 1966 constant, which olb used before, is 0.423. The three constants
    give r_0 values inside 0.3 % of each other.
    '''
    return FRIED_OVER_RHO0 * np.asarray(rho0, dtype=float)


@assumes(TILT_CONVENTION_G, FRESNEL_ZONE_CONSTRAINT)
def angle_of_arrival_variance(D, wavelength, z, cn2, *,
                              spectrum='kolmogorov', l0=None, L0=None,
                              radial=False):
    '''
    Return the angle-of-arrival (tilt) variance across a collecting lens [rad^2].

    TILT DEFINITION - THE OWNER MADE THIS CHOICE. This function returns the
    ANDREWS GRADIENT TILT (G-tilt). Andrews defines the tilt as the total phase
    difference across the pupil divided by the pupil width (Ch. 6, Eqs. (80) to
    (82), printed p. 200), which is what a CENTROID TRACKER measures. It is NOT
    the Noll Zernike tilt.

    The two definitions give two different coefficients:
        gradient tilt (this function, Andrews Ch. 6, Eq. (84), printed p. 201)
            <beta_a^2> = 0.174 (D/r_0)^(5/3) (lambda/D)^2  per axis
        Zernike tilt  (Noll, JOSA 66 (1976) 207, DOI 10.1364/JOSA.66.000207)
            <beta_a^2> = 0.182 (D/r_0)^(5/3) (lambda/D)^2  per axis
    A full-book search finds no 0.182 in Andrews and Phillips. See Conflict C-04
    in docs/andrews-crosscheck.md. Note that `olb/turbulence/ao.py` uses the NOLL
    convention (1.0299 and 0.134), so a caller that mixes the two must say which
    tilt it means.

    Parameters:
        D : float or numpy.ndarray
            Collecting-lens diameter [m]. The book writes it as 2 W_G.
        wavelength : float or numpy.ndarray
            Optical wavelength [m]. The Kolmogorov result does NOT depend on it.
            See the VALIDITY note below.
        z : float or numpy.ndarray
            Path length L [m].
        cn2 : float or numpy.ndarray
            Refractive-index structure constant [m^-2/3], constant over the path.
        spectrum : str
            "kolmogorov", "von_karman" or "modified". The two-scale branches use
            Eq. (83), which the book writes for the von Karman roll-off. The
            "modified" name is accepted and gives the same Eq. (83) numbers,
            because the book prints no separate modified row for Eq. (83).
        l0, L0 : float, optional
            Inner scale and outer scale [m].
        radial : bool
            False (the default) returns the PER-AXIS variance, which is what
            Eq. (84) prints. True returns the two-axis (radial) variance, which
            is twice the per-axis value.

    Returns:
        float or numpy.ndarray
            The tilt variance [rad^2].

    formula (Kolmogorov, 2 W_G >> l0):
        <beta_a^2> = 2.91 Cn2 L (2 W_G)^(-1/3)
    Source: Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196, Ch. 6,
    Eq. (84), printed p. 201, with the definition Eq. (82), printed p. 200. The
    slant-path version is Ch. 12, Eq. (28), printed p. 492, repeated as Ch. 12,
    Eq. (90), printed p. 522.

    formula (two scales, Ch. 6, Eq. (83), printed p. 201):
        <beta_a^2> = 1.64 Cn2 L l0^(-1/3) [1 - 0.72 (k0 l0)^(1/3)],  D << l0
        <beta_a^2> = 2.91 Cn2 L D^(-1/3)  [1 - 0.81 (k0 D)^(1/3)],   D >> l0
    with k0 = 2*pi/L0 and D = 2 W_G.

    THE RECAST to Fried units. Put Cn2 L = r_0^(-5/3)/(0.423 k^2) and
    k = 2*pi/lambda into Eq. (84):
        <beta_a^2> = 2.91 Cn2 L D^(-1/3)
                   = 2.91 r_0^(-5/3) lambda^2 D^(-1/3) / (0.423 * 4 pi^2)
                   = [2.91 / (0.423 * 4 pi^2)] (D/r_0)^(5/3) (lambda/D)^2
                   = 0.1743 (D/r_0)^(5/3) (lambda/D)^2.
    The module self-check measures that recast.

    VALIDITY. Eq. (83) is independent of the wavelength, but the book adds that
    this holds only when the Fresnel zone is small against the lens, that is
    sqrt(L/k) << D (Ch. 6, text below Eq. (83), printed p. 200). This function
    does not gate on that condition. The caller must test it.
    '''
    if spectrum not in _SPECTRA:
        raise ValueError(f'spectrum must be one of {_SPECTRA}, '
                         f'not {spectrum!r}')
    l0_val, k0 = _scales(spectrum, l0, L0)

    D = np.asarray(D, dtype=float)
    z = np.asarray(z, dtype=float)
    cn2 = np.asarray(cn2, dtype=float)
    cn2_l = cn2 * z

    large = 2.91 * cn2_l * D ** (-1.0 / 3.0)
    if spectrum == 'kolmogorov':
        out = large
    else:
        large = large * (1.0 - 0.81 * (k0 * D) ** (1.0 / 3.0))
        small = (1.64 * cn2_l * l0_val ** (-1.0 / 3.0)
                 * (1.0 - 0.72 * (k0 * l0_val) ** (1.0 / 3.0)))
        out = np.where(D < l0_val, small, large)
    return 2.0 * out if radial else out


def rms_image_jitter(focal_length_m, variance):
    '''
    Return the rms image displacement in the focal plane [m].

    formula:
        rms image jitter = f * sqrt(<beta_a^2>)
    Source: Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196, Ch. 6,
    Eq. (85), printed p. 201, with the text above it: "the rms image
    displacement is the rms angle of arrival multiplied by the focal length f of
    the collecting lens".

    Parameters:
        focal_length_m : float or numpy.ndarray
            Focal length of the collecting lens [m].
        variance : float or numpy.ndarray
            The tilt variance [rad^2], from `angle_of_arrival_variance`. Pass the
            per-axis variance for a per-axis jitter, and the radial variance for
            a radial jitter.
    '''
    return np.asarray(focal_length_m, dtype=float) * np.sqrt(
        np.asarray(variance, dtype=float))


if __name__ == '__main__':
    # ---------------- physics self-checks ----------------
    from .beam import beam_params

    lam_m = 1550e-9
    L = 2000.0
    cn2_ref = 1e-15
    k_ref = wavenumber(lam_m)
    l0_ref = 5e-3

    # The wave structure function grows with separation and with turbulence.
    d1 = wave_structure_function(0.01, lam_m, L, cn2_ref)
    d2 = wave_structure_function(0.02, lam_m, L, cn2_ref)
    assert d2 > d1 > 0.0
    assert wave_structure_function(0.01, lam_m, L, 2 * cn2_ref) > d1

    # The Kolmogorov exponent is 5/3.
    rr = np.array([0.01, 0.02, 0.04])
    slope = np.polyfit(np.log(rr),
                       np.log(wave_structure_function(rr, lam_m, L, cn2_ref)),
                       1)[0]
    assert abs(slope - 5.0 / 3.0) < 1e-10, slope

    # The spherical WSF is 1.093/2.914 of the plane WSF. Ch. 6, Eqs. (63)
    # and (70), printed pp. 194 and 195.
    ratio = (wave_structure_function(0.01, lam_m, L, cn2_ref, wave='spherical')
             / d1)
    assert abs(ratio - 1.093 / 2.914) < 1e-12, ratio

    # A Gaussian beam sits between the plane wave and the spherical wave.
    bp_plane = beam_params(50.0, lam_m, L)
    bp_sph = beam_params(1e-5, lam_m, L)
    bp = beam_params(0.05, lam_m, L)
    dg_plane = wave_structure_function(0.01, lam_m, L, cn2_ref, wave='gaussian',
                                       beam=bp_plane)
    dg_sph = wave_structure_function(0.01, lam_m, L, cn2_ref, wave='gaussian',
                                     beam=bp_sph)
    assert abs(dg_plane / d1 - 1.0) < 1e-3, dg_plane / d1
    ref_sph = wave_structure_function(0.01, lam_m, L, cn2_ref, wave='spherical')
    assert abs(dg_sph / ref_sph - 1.0) < 1e-3, dg_sph / ref_sph

    # a_factor: 8/3 in the plane limit, 1 in the spherical limit.
    assert abs(a_factor(1.0) - 8.0 / 3.0) < 1e-9
    assert abs(a_factor(0.0) - 1.0) < 1e-12

    # Below the inner scale the two-scale WSF turns quadratic in rho.
    tiny = np.array([1e-4, 2e-4, 4e-4])
    d_small = wave_structure_function(tiny, lam_m, L, cn2_ref,
                                      spectrum='von_karman', l0=l0_ref)
    slope_small = np.polyfit(np.log(tiny), np.log(d_small), 1)[0]
    assert abs(slope_small - 2.0) < 5e-3, slope_small

    # The Gaussian row of the modified spectrum is refused, not guessed.
    try:
        wave_structure_function(2e-4, lam_m, L, cn2_ref, wave='gaussian',
                                spectrum='modified', l0=l0_ref, beam=bp)
    except NotImplementedError:
        pass
    else:
        raise AssertionError('gaussian + modified must raise')

    # A finite outer scale reduces the WSF.
    with_outer = wave_structure_function(0.01, lam_m, L, cn2_ref,
                                         spectrum='von_karman', l0=l0_ref,
                                         L0=1.0)
    no_outer = wave_structure_function(0.01, lam_m, L, cn2_ref,
                                       spectrum='von_karman', l0=l0_ref)
    assert with_outer < no_outer, (with_outer, no_outer)

    # The coherence radius is the D = 2 point. The book row 1.46 is a rounded
    # 2.914/2 = 1.457, so the measured D is 1.996, not exactly 2.
    rho0 = coherence_radius(lam_m, L, cn2_ref)
    assert abs(wave_structure_function(rho0, lam_m, L, cn2_ref) - 2.0) < 5e-3

    # Stronger turbulence gives a smaller coherence radius.
    assert coherence_radius(lam_m, L, 10 * cn2_ref) < rho0

    # The spherical coherence radius is larger than the plane-wave one.
    rho_sp = coherence_radius(lam_m, L, cn2_ref, wave='spherical')
    assert rho_sp > rho0
    assert abs(rho_sp / rho0 - (1.46 / 0.55) ** 0.6) < 1e-12

    # The Gaussian coherence radius reduces to the plane and spherical rows.
    err_pl = abs(coherence_radius(lam_m, L, cn2_ref, wave='gaussian',
                                  beam=bp_plane) / rho0 - 1.0)
    g_sph = coherence_radius(lam_m, L, cn2_ref, wave='gaussian', beam=bp_sph)
    # The Gaussian row uses (1.46 ...)^(-3/5) with a = 1, so the spherical
    # comparison carries the 1.46 against 0.55 rounding of the book.
    err_sp = abs(g_sph / rho_sp - 1.0)
    assert err_pl < 1e-3, err_pl
    assert err_sp < 1e-2, err_sp

    # The Gaussian two-scale rows also reduce to the plane and spherical rows.
    for spec in ('von_karman',):
        pl = wave_structure_function(2e-4, lam_m, L, cn2_ref, spectrum=spec,
                                     l0=l0_ref, L0=2.0)
        gpl = wave_structure_function(2e-4, lam_m, L, cn2_ref, wave='gaussian',
                                      spectrum=spec, l0=l0_ref, L0=2.0,
                                      beam=bp_plane)
        sp = wave_structure_function(2e-4, lam_m, L, cn2_ref, wave='spherical',
                                     spectrum=spec, l0=l0_ref, L0=2.0)
        gsp = wave_structure_function(2e-4, lam_m, L, cn2_ref, wave='gaussian',
                                      spectrum=spec, l0=l0_ref, L0=2.0,
                                      beam=bp_sph)
        assert abs(gpl / pl - 1.0) < 2e-3, (spec, gpl / pl)
        assert abs(gsp / sp - 1.0) < 2e-3, (spec, gsp / sp)

    # The tilt variance falls as the lens grows, with the D^(-1/3) law.
    v1 = angle_of_arrival_variance(0.2, lam_m, L, cn2_ref)
    v2 = angle_of_arrival_variance(0.4, lam_m, L, cn2_ref)
    assert abs(v2 / v1 - 2.0 ** (-1.0 / 3.0)) < 1e-12
    # The radial variance is twice the per-axis variance.
    assert abs(angle_of_arrival_variance(0.2, lam_m, L, cn2_ref, radial=True)
               / v1 - 2.0) < 1e-12
    # A finite outer scale reduces the tilt.
    assert angle_of_arrival_variance(0.2, lam_m, L, cn2_ref,
                                     spectrum='von_karman', l0=l0_ref,
                                     L0=5.0) < v1
    # Image jitter is the focal length times the rms tilt.
    assert abs(rms_image_jitter(1.0, v1) - np.sqrt(v1)) < 1e-15

    # ---------------- REDUCTION checks ----------------
    from .. import ao
    from .. import gaussian_fried as gf
    from .. import plane_wave_scintillation as pws
    from ..profiles import DEFAULT_HS, get_c2n

    # 1. coherence_radius(plane) reproduces the old olb copy.
    old_rho = pws.coherence_radius(cn2_ref, lam_m, L)
    err_rho = abs(rho0 - old_rho) / old_rho
    assert err_rho < 1e-12, err_rho
    old_rho_gf = gf.plane_wave_coherence_radius(L, cn2_ref, lam_m)
    err_rho_gf = abs(rho0 - old_rho_gf) / old_rho_gf
    assert err_rho_gf < 1e-12, err_rho_gf
    print(f'REDUCTION coherence_radius(plane) : pws err = {err_rho:.3e}  '
          f'gaussian_fried err = {err_rho_gf:.3e}  (target 1e-12)')

    # 2. fried_parameter reproduces the old plane-wave Fried parameter.
    mine_r0 = fried_parameter(rho0)
    old_r0 = gf.plane_wave_fried_parameter(L, cn2_ref, lam_m)
    err_r0 = abs(mine_r0 - old_r0) / old_r0
    assert err_r0 < 1e-9, err_r0
    print(f'REDUCTION fried_parameter(plane) : rel err = {err_r0:.3e}  '
          f'(target 1e-9)')

    # 3. The spherical-wave Fried parameter. olb keeps the EXACT (8/3)^(3/5)
    # ratio (Conflict C-07), and the book row 0.55 is a rounded 0.5475. Report
    # both numbers, and assert only the olb wiring.
    old_r0_sp = gf.spherical_wave_fried_parameter(L, cn2_ref, lam_m)
    wired_sp = (8.0 / 3.0) ** 0.6 * mine_r0
    err_sp_wire = abs(wired_sp - old_r0_sp) / old_r0_sp
    assert err_sp_wire < 1e-9, err_sp_wire
    book_sp = fried_parameter(rho_sp)
    pct_sp = (book_sp - old_r0_sp) / old_r0_sp * 100.0
    print(f'REDUCTION spherical Fried : olb (8/3)^(3/5) wiring err = '
          f'{err_sp_wire:.3e} (target 1e-9); the book Table V row 0.55 '
          f'differs by {pct_sp:+.3f} % (Conflict C-07, keep the exact ratio)')

    # 4. The ao profile Fried parameter. The old copy used the Fried 1966
    # constant 0.423. The Andrews chain 2.1 (1.46 Cn2 k^2 L)^(-3/5) is the
    # equivalent of 0.4240. Report the measured shift.
    hs = DEFAULT_HS
    cn2_prof = get_c2n(hs, 21.0, 1.7e-14)
    airmass = 1.0 / np.sin(np.radians(60.0))
    moment = float(np.trapezoid(cn2_prof, hs)) * airmass
    mine_prof = fried_parameter(coherence_radius(lam_m, 1.0, moment))
    old_prof = ao.plane_wave_fried_parameter_profile(cn2_prof, hs, lam_m, 60.0)
    pct_prof = (mine_prof - old_prof) / old_prof * 100.0
    assert abs(pct_prof) < 0.2, pct_prof
    print(f'REDUCTION ao profile r0 : Andrews chain {mine_prof * 100:.4f} cm '
          f'against the old 0.423 copy {old_prof * 100:.4f} cm, '
          f'{pct_prof:+.4f} % (0.423 is Fried 1966, the book chain is 0.4240)')

    # 5. The gradient-tilt recast against 0.174 (D/r0)^(5/3)(lambda/D)^2.
    D_ref = 0.3
    r0_ref = fried_parameter(coherence_radius(lam_m, L, cn2_ref))
    mine_tilt = angle_of_arrival_variance(D_ref, lam_m, L, cn2_ref)
    recast = 0.174 * (D_ref / r0_ref) ** (5.0 / 3.0) * (lam_m / D_ref) ** 2
    pct_tilt = abs(mine_tilt - recast) / recast * 100.0
    assert pct_tilt < 2.0, pct_tilt
    noll = 0.182 * (D_ref / r0_ref) ** (5.0 / 3.0) * (lam_m / D_ref) ** 2
    print(f'REDUCTION gradient tilt Eq. (84) against the 0.174 recast : '
          f'{pct_tilt:.3f} % (target 2 %). The Noll Zernike tilt 0.182 gives '
          f'{noll / mine_tilt:.4f} times this value.')

    print(f'rho_0 = {rho0 * 100:.3f} cm   r_0 = {mine_r0 * 100:.3f} cm   '
          f'tilt rms = {np.sqrt(mine_tilt) * 1e6:.3f} urad')

    # ---------------- assumption annotations ----------------
    import warnings

    from ...assumptions import trace_assumptions

    mod = __name__

    # (1) value parity in and out of a collection context.
    ref_wsf = wave_structure_function(0.01, lam_m, L, cn2_ref)
    with trace_assumptions():
        traced_wsf = wave_structure_function(0.01, lam_m, L, cn2_ref)
    assert traced_wsf == ref_wsf, (traced_wsf, ref_wsf)

    # (2) registration: the expected sources and constraint kinds appear.
    with trace_assumptions() as trace:
        wave_structure_function(0.01, lam_m, L, cn2_ref)
        coherence_radius(lam_m, L, cn2_ref)
        angle_of_arrival_variance(0.3, lam_m, L, cn2_ref)
    for name in ('wave_structure_function', 'coherence_radius',
                 'angle_of_arrival_variance'):
        assert f'{mod}.{name}' in trace.records, name
    kinds = {c.kind for rec in trace.records.values() for c in rec.constraints}
    for expected in ('path-homogeneity', 'conflict', 'not-built', 'spectrum',
                     'tilt-convention', 'field-region'):
        assert expected in kinds, expected
    print('[assumes] structure sources and kinds register ok')

    # (3) an out-of-range call yields a source-prefixed violation, and the
    #     physics layer emits NO warning. A long path with a tiny lens breaks
    #     the sqrt(L/k) << D condition of Eq. (83).
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        with trace_assumptions() as trace:
            angle_of_arrival_variance(1e-3, lam_m, 200e3, cn2_ref)
    assert any(v.startswith(f'[{mod}.angle_of_arrival_variance]')
               for v in trace.violations), trace.violations
    assert len(caught) == 0, 'a check must not warn'
    print('[assumes] angle-of-arrival Fresnel-zone check fires ok')

    print('self-check passed')
