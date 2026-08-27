'''
Sampling constraints for a numerical Fresnel propagation.

This module answers one question: what must hold for the answer to be right?
A wave-optics propagation is a discrete Fourier transform. A discrete Fourier
transform aliases. The book turns the propagation geometry into four numbered
inequalities. This module gives each inequality as a small pure function, and
it gives one checker that measures a grid against all of them.

Source of every equation:
    J. D. Schmidt, "Numerical Simulation of Optical Wave Propagation with
    Examples in MATLAB", SPIE Press (2010). DOI: 10.1117/3.866274
Chapter 6, Secs. 6.3 and 6.4, printed pp. 90 to 102 (the three kernels);
Chapter 7, printed pp. 115 to 131 (the four constraints); Chapter 8, printed
pp. 133 to 147 (the partial-propagation rules). Each function names its
chapter, its equation number, and its printed page.

THE FOUR CONSTRAINTS (Ch. 7, Sec. 7.3.3, printed p. 127):

    1. delta2 <= -(D2/D1) delta1 + lambda z / D1
    2. N      >= D1/(2 delta1) + D2/(2 delta2) + lambda z/(2 delta1 delta2)
    3. (1 + z/R) delta1 - lambda z/D1 <= delta2
                                      <= (1 + z/R) delta1 + lambda z/D1
    4. N      >= lambda z / (delta1 delta2)

Constraints 1 and 2 come from the propagation GEOMETRY only. They hold for
every kernel. Constraints 3 and 4 stop the aliasing of two quadratic phase
factors, so they belong to the kernel that carries those factors.

SYMBOLS. `D1` is the maximum spatial extent of the source field [m]. `D2` is
the extent of the observation-plane region of interest [m]. `delta1` and
`delta2` are the grid spacings in the two planes [m]. `N` is the number of
grid points along one side. `z` is the propagation distance [m]. `R` is the
radius of curvature of the source wavefront [m]; R < 0 is a diverging beam and
R > 0 is a converging beam (Ch. 7, Eq. (7.32), printed p. 122). An infinite R
is a flat wavefront.

VALIDITY OF THE WHOLE MODULE. Every rule below assumes the Fresnel (paraxial)
approximation and a scalar field. The book states the rules are a GUIDELINE,
not an unbreakable law: "simulations must be approached carefully and
validated fully" (Ch. 7, Sec. 7.3.3, printed p. 129). A grid that just passes
a bound can still give a result that does not match theory.

This module holds physics only. It imports numpy only. It returns no decibels,
and it imports nothing from the rest of olb.
'''

import math
from collections import namedtuple

import numpy as np

BOOK = 'Schmidt (2010), DOI 10.1117/3.866274'


def _cite(chapter, equation, page):
    '''Build one citation string.'''
    return f'{BOOK}, Ch. {chapter}, Eq. ({equation}), printed p. {page}'


# ---------------------------------------------------------------------------
# Ch. 7.1 and 7.2: the band limit and the propagation geometry
# ---------------------------------------------------------------------------

def nyquist_max_angle(wavelength_m, delta):
    '''
    Return the largest ray angle that a grid spacing can hold.

    formula:
        theta_max = lambda / (2 delta)
    Source: Schmidt (2010), DOI 10.1117/3.866274, Ch. 7, Eq. (7.7), printed
    p. 117. Eq. (7.6), printed p. 116, is the same rule as a bound on delta.

    The book maps the spatial-frequency spectrum of a field onto its plane-wave
    spectrum through the direction cosines alpha = lambda fx, beta = lambda fy
    (Ch. 7, Eq. (7.5), printed p. 116). So the Nyquist bound delta <= 1/(2
    fmax) becomes a bound on the ray angle.

    VALIDITY. The angle is paraxial. The book uses the small-angle form, so the
    result is a slope, not an arc.

    Parameters:
        wavelength_m : float
            The optical wavelength [m].
        delta : float
            The grid spacing [m].

    Returns:
        The largest angle [rad].
    '''
    return wavelength_m / (2.0 * delta)


def geometric_max_angle(D1, D2, delta1, delta2, z):
    '''
    Return the largest ray angle that the propagation geometry makes.

    formula:
        theta_edges = (D1 + D2) / (2 z)
        theta_k     = (D1 / (2 z)) (delta2/delta1 - 1)
        theta_max   = theta_edges + theta_k
                    = (D1 delta2/delta1 + D2) / (2 z)
    Source: Schmidt (2010), DOI 10.1117/3.866274, Ch. 7, Eqs. (7.8), (7.9) and
    (7.12), printed p. 118.

    A point at the edge of the source must fully illuminate the observation
    region. theta_edges is the angle from the bottom edge of the source to the
    top edge of the observation region. theta_k is the tilt of the virtual
    spherical wave inside the Fresnel integral, which the ratio of the two grid
    sizes sets.

    VALIDITY. The angles are paraxial, and the geometry is one-dimensional. The
    book states the result generalises to two dimensions (Ch. 7, Sec. 7.2,
    printed p. 117). The wavelength is not in the formula.

    Parameters:
        D1, D2 : float
            The source extent and the observation-region extent [m].
        delta1, delta2 : float
            The grid spacings [m].
        z : float
            The propagation distance [m].

    Returns:
        The largest ray angle [rad].
    '''
    return (D1 * delta2 / delta1 + D2) / (2.0 * z)


def constraint1_max_delta2(D1, D2, delta1, wavelength_m, z):
    '''
    Return the largest observation-plane grid spacing (constraint 1).

    formula:
        delta2 <= -(D2/D1) delta1 + lambda z / D1
    Source: Schmidt (2010), DOI 10.1117/3.866274, Ch. 7, Eq. (7.14), printed
    p. 119. Eq. (7.13), printed p. 118, is the step before the algebra.

    The rule joins the geometric ray angle of Eq. (7.12) to the Nyquist angle
    of Eq. (7.7). The Nyquist angle takes the SOURCE grid spacing, so Eq.
    (7.13) reads (D1 delta2/delta1 + D2)/(2z) <= lambda/(2 delta1). It makes
    the grid sample the spatial bandwidth that reaches the observation-plane
    region of interest.

    VALIDITY. Geometry only. The rule holds for all three kernels. The bound
    goes negative when D2 delta1 > lambda z; no grid then satisfies the rule,
    and the geometry needs a smaller delta1 or a smaller D2.

    Parameters:
        D1, D2 : float
            The source extent and the observation-region extent [m].
        delta1 : float
            The source-plane grid spacing [m].
        wavelength_m : float
            The optical wavelength [m].
        z : float
            The propagation distance [m].

    Returns:
        The largest delta2 [m].
    '''
    return (wavelength_m * z - D2 * delta1) / D1


def illuminated_diameter(D1, delta1, delta2, wavelength_m, z):
    '''
    Return the diameter of the illuminated area in the observation plane.

    formula:
        D_illum = D1 delta2/delta1 + 2 theta_max z
                = D1 delta2/delta1 + lambda z / delta1
    Source: Schmidt (2010), DOI 10.1117/3.866274, Ch. 7, Eqs. (7.15) and
    (7.16), printed p. 119.

    The source scales up by the grid ratio delta2/delta1, and the maximum
    angular content of the grid spreads it by lambda z / delta1.

    VALIDITY. Paraxial. It is the input to the grid-extent rule of Eq. (7.18).

    Returns:
        The illuminated diameter [m].
    '''
    return D1 * delta2 / delta1 + wavelength_m * z / delta1


def constraint2_min_n(D1, D2, delta1, delta2, wavelength_m, z):
    '''
    Return the smallest number of grid points (constraint 2).

    formula:
        D_grid >= (D_illum + D2) / 2
        N      >= D1/(2 delta1) + D2/(2 delta2) + lambda z/(2 delta1 delta2)
    Source: Schmidt (2010), DOI 10.1117/3.866274, Ch. 7, Eqs. (7.18) and
    (7.20), printed p. 120.

    The book lets the light alias, but not into the region of interest. The
    light that leaves the grid wraps around to the other side. The grid extent
    must be the mean of the illuminated area and the observation region, so the
    wrapped light gets only half way around.

    VALIDITY. Geometry only. The rule holds for all three kernels. The book
    proves that constraint 4 is more restrictive than constraints 1 and 2
    together (Ch. 7, Sec. 7.3.3, printed p. 128, and Problem 6, printed
    p. 131), so an angular-spectrum grid that passes constraint 4 passes this
    one.

    Returns:
        The smallest N. It is a float, not an integer. Round it up, and the
        book rounds again to the next power of two for the FFT.
    '''
    return (D1 / (2.0 * delta1) + D2 / (2.0 * delta2)
            + wavelength_m * z / (2.0 * delta1 * delta2))


# ---------------------------------------------------------------------------
# Ch. 7.3: the local spatial frequency of the two quadratic phase factors
# ---------------------------------------------------------------------------

def local_spatial_frequency_source(x1, wavelength_m, z, R=math.inf, m=None):
    '''
    Return the local spatial frequency of the source-plane quadratic phase.

    formula:
        f_loc     = grad(phi) / (2 pi)
        phi       = (k/2) c |r1|^2
        f_loc,x   = c x1 / lambda
        c = 1/z + 1/R              Fresnel-integral kernels,     m is None
        c = (1 - m)/z + 1/R        angular-spectrum kernel,      m = d2/d1
    Source: Schmidt (2010), DOI 10.1117/3.866274, Ch. 7, Eq. (7.37), printed
    p. 122 (the definition), Eq. (7.39), printed p. 123 (the Fresnel-integral
    case) and Eq. (7.51), printed p. 126 (the angular-spectrum case).

    The book models the source as an apodised beam of extent D1 with a
    parabolic wavefront of radius R (Ch. 7, Eq. (7.32), printed p. 122). Each
    kernel then holds ONE source-plane quadratic phase factor. Its local
    frequency grows with x1, so it is largest at the edge of the source,
    x1 = D1/2. The Nyquist rule at that edge gives the Fresnel-integral minimum
    distance (m is None) or constraint 3 (m is set).

    VALIDITY. Lambert and Fraser show the bandwidth of the product of the
    aperture and the phase is set by the APERTURE for a very small aperture,
    and by the PHASE at the edge for a larger aperture (Ch. 7, text below
    Eq. (7.36), printed p. 122). The book takes the second case, and so does
    this module. A source much smaller than one Fresnel zone breaks it.

    Parameters:
        x1 : float or numpy.ndarray
            The source-plane coordinate [m].
        wavelength_m : float
            The optical wavelength [m].
        z : float
            The propagation distance [m].
        R : float
            The radius of curvature of the source wavefront [m]. R < 0
            diverges, R > 0 converges, and math.inf is flat.
        m : float, optional
            The scaling parameter delta2/delta1 of the angular-spectrum
            kernel. None gives the Fresnel-integral form.

    Returns:
        The local spatial frequency [1/m].
    '''
    curvature = (1.0 if m is None else 1.0 - m) / z + 1.0 / R
    return curvature * np.asarray(x1, dtype=float) / wavelength_m


def local_spatial_frequency_transfer(f1, wavelength_m, z, delta1, delta2):
    '''
    Return the local spatial frequency of the angular-spectrum transfer phase.

    formula:
        phi'      = 2 pi^2 (lambda z / (m k)) |f1|^2,   m = delta2/delta1
        f'_loc,x  = lambda delta1 z f1x / delta2
    Source: Schmidt (2010), DOI 10.1117/3.866274, Ch. 7, Eqs. (7.55) and
    (7.57), printed p. 126.

    This is the SECOND quadratic phase factor. It lives in the spatial-
    frequency domain, so its local frequency has units of length. It is largest
    at the edge of the frequency grid, f1x = 1/(2 delta1). Constraint 4 comes
    from this result.

    VALIDITY. The Fresnel-integral kernels have no such factor, so this result
    and constraint 4 belong to the angular-spectrum kernel only.

    Parameters:
        f1 : float or numpy.ndarray
            The source-plane spatial frequency [1/m].
        wavelength_m, z, delta1, delta2 : float
            The wavelength [m], the distance [m], and the two spacings [m].

    Returns:
        The local spatial frequency of the transfer phase [m].
    '''
    return wavelength_m * delta1 * z * np.asarray(f1, dtype=float) / delta2


def constraint3_delta2_window(D1, delta1, wavelength_m, z, R=math.inf):
    '''
    Return the window of observation-plane grid spacings (constraint 3).

    formula:
        (1 + z/R) delta1 - lambda z / D1 <= delta2
                                         <= (1 + z/R) delta1 + lambda z / D1
    Source: Schmidt (2010), DOI 10.1117/3.866274, Ch. 7, Eq. (7.53), printed
    p. 126. Ch. 8, Eq. (8.22), printed p. 143, proves the SAME window holds for
    a chain of partial propagations, with delta_n in the place of delta2.

    The rule applies the Nyquist criterion to the local frequency of the
    source-plane quadratic phase (Eq. (7.52), printed p. 126) at the edge of
    the source aperture.

    VALIDITY. Angular-spectrum kernel only. The window is two-sided, because
    the phase curvature (1 - m)/z + 1/R changes sign with m = delta2/delta1.

    Returns:
        A (lower, upper) tuple of delta2 bounds [m]. The lower bound is often
        negative, which means no lower bound is active.
    '''
    slope = 1.0 + z / R
    spread = wavelength_m * z / D1
    return (slope * delta1 - spread, slope * delta1 + spread)


def constraint3_is_slack(D1, D2, z, R=math.inf):
    '''
    Say whether constraint 3 can be ignored.

    formula:
        1 + z/R < D2/D1
    Source: Schmidt (2010), DOI 10.1117/3.866274, Ch. 7, Eq. (7.60), printed
    p. 129.

    The book compares the slopes and the intercepts of constraints 1 and 3 in
    the (delta1, delta2) plane. When this test passes, constraint 1 is the
    tighter upper bound and the lower bound of constraint 3 is not reached. The
    physical reading is that the geometric beam stays inside a region of
    diameter D2 (Ch. 7, text below Eq. (7.60), printed p. 129).

    VALIDITY. Angular-spectrum kernel only. The test says nothing about
    constraint 4.

    Returns:
        True when constraint 3 is not a factor.
    '''
    return (1.0 + z / R) < (D2 / D1)


def constraint4_min_n(delta1, delta2, wavelength_m, z):
    '''
    Return the smallest number of grid points (constraint 4).

    formula:
        lambda z / (2 delta2) <= N delta1 / 2
        N >= lambda z / (delta1 delta2)
    Source: Schmidt (2010), DOI 10.1117/3.866274, Ch. 7, Eqs. (7.58) and
    (7.59), printed pp. 126 and 127.

    The rule applies the Nyquist criterion to the transfer-function phase at
    the edge of the frequency grid.

    VALIDITY. Angular-spectrum kernel only. This constraint holds the
    propagation METHOD, not the geometry, so it is the one that a long
    propagation breaks first. The book states it is usually the culprit when
    the required N is prohibitively large, and Ch. 8 relaxes it with partial
    propagations (Ch. 8, text, printed p. 133). With delta2 = delta1 the rule
    inverts to the familiar range limit z <= N delta1^2 / lambda.

    Returns:
        The smallest N, as a float.
    '''
    return wavelength_m * z / (delta1 * delta2)


# ---------------------------------------------------------------------------
# Ch. 6.3 and 7.3.1: the two Fresnel-integral kernels
# ---------------------------------------------------------------------------

def one_step_delta2(N, delta1, wavelength_m, z):
    '''
    Return the observation-plane grid spacing that one-step Fresnel FIXES.

    formula:
        delta2 = lambda z / (N delta1)
    Source: Schmidt (2010), DOI 10.1117/3.866274, Ch. 7, Eq. (7.21), printed
    p. 120. Ch. 6, Eq. (6.16), printed p. 90, gives the same result.

    KERNEL ASSUMPTIONS, one-step Fresnel integral:
    - The Fresnel (paraxial) approximation holds. The kernel is one FT between
      two confocal spheres of radius z (Ch. 6, Sec. 6.3.1, printed p. 90).
    - The kernel FIXES delta2. The caller has NO freedom: the frequency-domain
      spacing is 1/(N delta1), and the FT maps f1 = r2/(lambda z). A finer
      observation grid needs a larger N or a different geometry.
    - It ALIASES when the source-plane quadratic phase of curvature 1/z + 1/R
      turns faster than one sample at the edge of the source. See
      `fresnel_min_distance`. So the kernel is valid for LONG propagations
      only (Ch. 7, Sec. 7.3, printed p. 120).
    - The GOVERNING constraint is Eq. (7.25) through `one_step_min_n`, which is
      constraints 1 and 2 with delta2 substituted. The book proves the two
      give the identical bound (Eq. (7.31), printed p. 121). Constraints 3 and
      4 do NOT apply, because the kernel has no transfer-function phase factor.

    Returns:
        The fixed delta2 [m].
    '''
    return wavelength_m * z / (N * delta1)


def one_step_min_n(D1, D2, delta1, wavelength_m, z):
    '''
    Return the smallest N for a one-step Fresnel-integral propagation.

    formula:
        N >= D1 lambda z / (delta1 (lambda z - D2 delta1))
    Source: Schmidt (2010), DOI 10.1117/3.866274, Ch. 7, Eq. (7.25), printed
    p. 121. Eq. (7.31), printed p. 121, gets the identical bound from
    constraint 2 instead of constraint 1.

    KERNEL ASSUMPTIONS: see `one_step_delta2`. This function is the governing
    constraint of that kernel.

    VALIDITY. The book states two properties of the bound (text below
    Eq. (7.31), printed p. 121): the geometry must have lambda z > D2 delta1,
    because N can only be positive; and N goes to infinity as lambda z
    approaches D2 delta1.

    Returns:
        The smallest N, as a float. It is math.inf when lambda z <= D2 delta1.
    '''
    den = delta1 * (wavelength_m * z - D2 * delta1)
    if den <= 0.0:
        return math.inf
    return D1 * wavelength_m * z / den


def fresnel_min_distance(D1, delta1, wavelength_m, R=math.inf):
    '''
    Return the shortest distance that a Fresnel-integral kernel may propagate.

    formula:
        z >= D1 delta1 R / (lambda R - D1 delta1)   for a finite R
        z >= D1 delta1 / lambda                     for an infinite R
    Source: Schmidt (2010), DOI 10.1117/3.866274, Ch. 7, Eqs. (7.41) and
    (7.42), printed p. 123. Eq. (7.40), printed p. 123, is the Nyquist step.

    The bound stops the aliasing of the source-plane quadratic phase factor.
    This is why the book says the Fresnel-integral kernels suit LONG
    propagations and the angular-spectrum kernel suits SHORT propagations
    (Ch. 7, Sec. 7.3, printed p. 120).

    VALIDITY. The book calls it "just a guideline" and warns that a simulation
    near the minimum distance may not match an analytic result (Ch. 7, text
    below Eq. (7.42), printed p. 123). The printed form is one-sided; a
    converging source with lambda R < D1 delta1 has no valid distance, and this
    function returns math.inf there.

    Returns:
        The shortest distance [m], or math.inf when no distance works.
    '''
    if math.isinf(R):
        return D1 * delta1 / wavelength_m
    den = wavelength_m * R - D1 * delta1
    if den == 0.0:
        return math.inf
    z_min = D1 * delta1 * R / den
    return z_min if z_min > 0.0 else math.inf


def two_step_planes(z, m):
    '''
    Return the two intermediate-plane geometries of two-step Fresnel.

    formula:
        m = delta2/delta1 = Delta_z2 / Delta_z1
        Delta_z1 = z / (1 -+ m),   Delta_z2 = z (-+ m) / (1 -+ m)
    Source: Schmidt (2010), DOI 10.1117/3.866274, Ch. 6, Eqs. (6.24) to (6.29),
    printed pp. 94 and 95. Eq. (6.30), printed p. 95, proves the ratio.
    Table 6.2, printed p. 95, prints the example values.

    KERNEL ASSUMPTIONS, two-step Fresnel integral:
    - The Fresnel (paraxial) approximation holds. The kernel is two one-step
      Fresnel integrals through an intermediate plane at z1a.
    - The kernel FREES the magnification. The caller picks m, and the position
      of the intermediate plane follows. This is the whole reason for the
      second FT (Ch. 6, Sec. 6.3.2, printed p. 92). The cost is one more FT.
    - There are always TWO solutions. The "minus" branch puts the intermediate
      plane OUTSIDE the source-to-observation span (Fig. 6.2, printed p. 93),
      and the "plus" branch puts it BETWEEN the two planes (Fig. 6.3, printed
      p. 94). At m = 1 the minus branch runs to infinity and the plus branch
      sits half way.
    - It ALIASES under the same source-phase rule as the one-step kernel, now
      applied to EACH step. Give `fresnel_min_distance` each of Delta_z1 and
      Delta_z2 with the grid spacing of that step.
    - The GOVERNING constraints are constraints 1 and 2 on the pair
      (delta1, delta2). Constraints 3 and 4 do NOT apply.

    VALIDITY. A negative Delta_z is a real geometry in the book: it is a
    backward step to a plane behind the source. The book keeps it. This module
    keeps it too, and the caller must decide whether the propagator accepts it.

    Parameters:
        z : float
            The source-to-observation distance [m].
        m : float
            The wanted magnification delta2/delta1. It must be positive.

    Returns:
        A dict with the keys "minus" and "plus". Each value is a
        (Delta_z1, Delta_z2) tuple [m], and the two add up to z.
    '''
    if m <= 0.0:
        raise ValueError(
            f'two_step_planes needs a positive magnification, not {m!r}. '
            f'{_cite(6, "6.24", 94)}')
    out = {}
    for name, sign in (('minus', -1.0), ('plus', 1.0)):
        den = 1.0 + sign * m
        if den == 0.0:
            out[name] = (math.inf, -math.inf)
            continue
        dz1 = z / den
        out[name] = (dz1, z * sign * m / den)
    return out


def angular_spectrum_max_z(N, delta1, wavelength_m):
    '''
    Return the longest angular-spectrum step for the delta2 = delta1 case.

    formula:
        z <= N delta1^2 / lambda        (constraint 4 with delta2 = delta1)
    Source: Schmidt (2010), DOI 10.1117/3.866274, Ch. 7, Eq. (7.59), printed
    p. 127, inverted for z.

    KERNEL ASSUMPTIONS, angular spectrum:
    - The Fresnel (paraxial) approximation holds. The kernel is the convolution
      form of the Fresnel integral: an FT, a multiply by the transfer function,
      and an inverse FT (Ch. 6, Sec. 6.4, printed p. 96).
    - The kernel LEAVES delta1 and delta2 free and independent. A plain
      convolution would force delta2 = delta1 (m = 1); the book adds the
      scaling parameter m to lift that (Ch. 6, text below Eq. (6.32), printed
      p. 97).
    - It ALIASES when the transfer-function phase turns faster than one sample
      at the edge of the frequency grid. That happens as z grows, so the kernel
      is valid for SHORT propagations only. The wrap-around then creeps in from
      the edge of the grid (Ch. 8, text, printed p. 133).
    - The GOVERNING constraint is constraint 4, `constraint4_min_n`. It is more
      restrictive than constraints 1 and 2 together (Ch. 7, Sec. 7.3.3, printed
      p. 128). Constraint 3 also applies, and `constraint3_is_slack` says when
      it does not bite.
    - Ch. 8 relaxes constraint 4 with a chain of partial propagations. See
      `partial_max_step`.

    VALIDITY. This function is the m = 1 special case only. With delta2 not
    equal to delta1, call `constraint4_min_n` instead.

    Returns:
        The longest well-sampled step [m].
    '''
    return N * delta1 ** 2 / wavelength_m


# ---------------------------------------------------------------------------
# Ch. 8: partial propagations
# ---------------------------------------------------------------------------

def partial_grid_spacing(delta1, delta_n, alpha):
    '''
    Return the grid spacing at a fractional distance along the path.

    formula:
        delta_i = (1 - alpha_i) delta1 + alpha_i delta_n,   alpha_i = z_i / z
    Source: Schmidt (2010), DOI 10.1117/3.866274, Ch. 8, Table 8.2, printed
    p. 140. Ch. 8, Eq. (8.8), printed p. 136, derives the two-plane case from
    similar triangles.

    The grid spacing of a partial-propagation chain is LINEAR in the distance.
    The chain therefore has one straight ray cone from the source grid to the
    observation grid, and every intermediate grid sits on that cone.

    VALIDITY. The rule is geometric. It holds for an expanding chain
    (delta_n > delta1) and for a contracting chain.

    Parameters:
        delta1, delta_n : float
            The spacings in the first and the last plane [m].
        alpha : float or numpy.ndarray
            The fractional distance z_i / z, in [0, 1].

    Returns:
        The grid spacing [m], with the shape of alpha.
    '''
    alpha = np.asarray(alpha, dtype=float)
    return (1.0 - alpha) * delta1 + alpha * delta_n


def partial_max_step(N, delta1, delta_n, wavelength_m):
    '''
    Return the longest partial-propagation step.

    formula:
        Delta_z_max = min(delta1, delta_n)^2 N / lambda
    Source: Schmidt (2010), DOI 10.1117/3.866274, Ch. 8, Eq. (8.24), printed
    p. 144.

    Constraint 4 holds for EACH step: N >= lambda Delta_z_i / (delta_i
    delta_i+1) (Ch. 8, Eq. (8.23), printed p. 143). The book cannot solve that
    chain, because the spacings depend on the step positions that the chain is
    trying to find. It takes the worst case instead: the smallest spacing of
    the chain is delta1 or delta_n, because `partial_grid_spacing` is linear.
    The rule then becomes one bound on the step length.

    VALIDITY. Angular-spectrum kernel only. The bound is conservative, so a
    shorter step always stays valid: "One can always use more partial
    propagations" (Ch. 8, Sec. 8.4, printed p. 144). The rule assumes
    constraints 1, 2 and 3 already fixed N, delta1 and delta_n.

    Returns:
        The longest step [m].
    '''
    return min(delta1, delta_n) ** 2 * N / wavelength_m


def partial_plane_count(z, N, delta1, delta_n, wavelength_m):
    '''
    Return the smallest number of propagation PLANES of a partial chain.

    formula:
        n = ceil(z / Delta_z_max) + 1        planes
        n - 1                                partial propagations
    Source: Schmidt (2010), DOI 10.1117/3.866274, Ch. 8, text below Eq. (8.24),
    printed p. 144.

    VALIDITY. Angular-spectrum kernel only. The count is a MINIMUM. The book
    states more planes always stay valid. The chain is mathematically the same
    as one full propagation; the gain is that an absorbing boundary at each
    plane removes the wrap-around all along the path (Ch. 8, text, printed
    p. 133).

    Returns:
        The number of planes, an integer. It is 2 or more.
    '''
    dz_max = partial_max_step(N, delta1, delta_n, wavelength_m)
    return max(2, int(math.ceil(z / dz_max)) + 1)


def absorbing_boundary_sigma(N, frac=0.47):
    '''
    Return the width of the absorbing band, in grid points.

    formula:
        sigma = frac * N        [grid points]
    Source: Schmidt (2010), DOI 10.1117/3.866274, Ch. 8, Listing 8.1, printed
    p. 142, uses frac = 0.47 with a super-Gaussian of exponent n = 16.
    Fig. 8.1, printed p. 134, plots sigma = 0.45 L with n = 16.

    This function gives the SIZING rule only. The shape of the absorber is
    Eq. (8.1), printed p. 134, and it belongs to another module.

    VALIDITY. The band must be close to unity in the middle of the grid and
    close to zero at the edge, because the absorber must not touch the light in
    the region of interest (Ch. 8, Sec. 8.1, printed p. 134). The book states
    the Hamming and the Bartlett windows do NOT meet that rule, and that the
    super-Gaussian (n > 2) and the Tukey window do. A `frac` near 0.5 puts the
    half-power point near the edge of the grid.

    Parameters:
        N : int
            The number of grid points along one side.
        frac : float
            The width as a fraction of N.

    Returns:
        The width sigma, in grid points.
    '''
    return frac * N


# ---------------------------------------------------------------------------
# The checker
# ---------------------------------------------------------------------------

Rule = namedtuple('Rule', 'name satisfied bound actual citation')
Rule.__doc__ = '''One measured sampling rule. `bound` is a float, or a
(lower, upper) tuple for the two-sided constraint 3.'''


def check_sampling(D1, D2, delta1, delta2, N, wavelength_m, z, R=math.inf):
    '''
    Measure one grid against every sampling rule of Ch. 7.

    The function NEVER raises on a violation, and it never warns. It returns
    the measurements, and the caller decides what to do.

    Source: Schmidt (2010), DOI 10.1117/3.866274, Ch. 7, Sec. 7.3.3, printed
    p. 127 (the four constraints), and Ch. 7, Eqs. (7.41) and (7.42), printed
    p. 123 (the fifth row).

    Rows 1 and 2 are GEOMETRY. They hold for all three kernels. Rows 3 and 4
    hold for the angular-spectrum kernel only. Row 5 holds for the two
    Fresnel-integral kernels only. Each row carries its citation, so the caller
    can select the rows of its kernel.

    Parameters:
        D1, D2 : float
            The source extent and the observation-region extent [m].
        delta1, delta2 : float
            The grid spacings [m]. For a partial-propagation chain, pass
            delta_n as delta2 (Ch. 8, Sec. 8.4, printed p. 143).
        N : int
            The number of grid points along one side.
        wavelength_m : float
            The optical wavelength [m].
        z : float
            The propagation distance [m].
        R : float
            The radius of curvature of the source wavefront [m].

    Returns:
        A list of five Rule tuples.
    '''
    d2_max = constraint1_max_delta2(D1, D2, delta1, wavelength_m, z)
    n2 = constraint2_min_n(D1, D2, delta1, delta2, wavelength_m, z)
    lo, hi = constraint3_delta2_window(D1, delta1, wavelength_m, z, R)
    n4 = constraint4_min_n(delta1, delta2, wavelength_m, z)
    z_min = fresnel_min_distance(D1, delta1, wavelength_m, R)
    return [
        Rule('constraint 1: observation grid spacing', delta2 <= d2_max,
             d2_max, delta2, _cite(7, '7.14', 119)),
        Rule('constraint 2: grid points, geometry', N >= n2, n2, float(N),
             _cite(7, '7.20', 120)),
        Rule('constraint 3: source quadratic phase', lo <= delta2 <= hi,
             (lo, hi), delta2, _cite(7, '7.53', 126)),
        Rule('constraint 4: transfer-function phase', N >= n4, n4, float(N),
             _cite(7, '7.59', 127)),
        Rule('Fresnel-integral minimum distance', z >= z_min, z_min, z,
             _cite(7, '7.41', 123)),
    ]


if __name__ == '__main__':
    # ---------------- the book's worked example, Ch. 7.3.1 ----------------
    # One step of Fresnel-integral propagation from a square aperture.
    # Source: Ch. 7, Eq. (7.43), printed p. 123, and Listing 7.1, printed
    # p. 124.
    D1, D2, lam, z = 2e-3, 3e-3, 1e-6, 0.5
    d1 = D1 / 50.0                                   # 50 points across D1

    n_min = one_step_min_n(D1, D2, d1, lam, z)
    assert 65.0 < n_min < 66.0, n_min               # the book says 66 points
    n_use = 2 ** math.ceil(math.log2(n_min))
    assert n_use == 128, n_use                       # the book says 128
    d2 = one_step_delta2(n_use, d1, lam, z)
    assert abs(d2 - 97.7e-6) < 0.1e-6, d2            # the book says 97.7 um
    z_min = fresnel_min_distance(D1, d1, lam)
    assert abs(z_min - 0.08) < 1e-12, z_min          # the book says 8 cm
    assert z > 2.0 * z_min                           # "much farther"
    # The derivation closes: at z_min the local frequency of Eq. (7.39) at the
    # edge of the source is exactly the Nyquist limit 1/(2 delta1).
    f_edge_fr = float(local_spatial_frequency_source(D1 / 2.0, lam, z_min))
    assert abs(f_edge_fr / (1.0 / (2.0 * d1)) - 1.0) < 1e-12, f_edge_fr

    print('BOOK Ch. 7, Eq. (7.43), printed p. 123 (one-step Fresnel):')
    print(f'  delta1                {d1 * 1e6:9.2f} um   (book 40 um)')
    print(f'  N_min, Eq. (7.25)     {n_min:9.2f}      (book 66)')
    print(f'  N used                {n_use:9d}      (book 128)')
    print(f'  delta2, Eq. (7.21)    {d2 * 1e6:9.2f} um   (book 97.7 um)')
    print(f'  z_min, Eq. (7.42)     {z_min * 100:9.2f} cm   (book 8 cm)')

    # Eq. (7.25) and Eq. (7.31) are the same bound. Check the second route:
    # constraint 2 with delta2 substituted must give the identical number.
    n_via2 = constraint2_min_n(D1, D2, d1, one_step_delta2(n_min, d1, lam, z),
                               lam, z)
    assert abs(n_via2 / n_min - 1.0) < 1e-12, (n_via2, n_min)
    print(f'  Eq. (7.31) route      {n_via2:9.2f}      (identical to (7.25))')

    # ---------------- the book's worked example, Ch. 7.3.2 ----------------
    # The angular-spectrum method. Source: Ch. 7, text and Listing 7.2,
    # printed pp. 127 and 128.
    D1a, D2a, za, lama = 2e-3, 4e-3, 0.1, 1e-6
    d1a, d2a = 9.4848e-6, 28.1212e-6

    n4 = constraint4_min_n(d1a, d2a, lama, za)
    n2 = constraint2_min_n(D1a, D2a, d1a, d2a, lama, za)
    assert abs(math.log2(n4) - 8.55) < 0.01, math.log2(n4)   # book 2^8.55
    assert abs(math.log2(n2) - 8.51) < 0.01, math.log2(n2)   # book 2^8.51
    assert n4 > n2                                   # book: 4 is the tighter
    n_use_a = 2 ** math.ceil(math.log2(n4))
    assert n_use_a == 512, n_use_a                   # the book says 512

    # The book says constraint 1 is more restrictive than constraint 3.
    d2_max = constraint1_max_delta2(D1a, D2a, d1a, lama, za)
    lo3, hi3 = constraint3_delta2_window(D1a, d1a, lama, za)
    assert d2_max < hi3, (d2_max, hi3)
    assert d2a <= d2_max and lo3 <= d2a <= hi3
    assert constraint3_is_slack(D1a, D2a, za)        # 1 + z/R = 1 < D2/D1 = 2

    print('')
    print('BOOK Ch. 7, Sec. 7.3.2, printed p. 127 (angular spectrum):')
    print(f'  log2 N, constraint 4  {math.log2(n4):9.2f}      (book 8.55)')
    print(f'  log2 N, constraint 2  {math.log2(n2):9.2f}      (book 8.51)')
    print(f'  N used                {n_use_a:9d}      (book 512)')
    print(f'  delta2 max, cons. 1   {d2_max * 1e6:9.2f} um   (delta2 = '
          f'{d2a * 1e6:.2f} um)')
    print(f'  delta2 max, cons. 3   {hi3 * 1e6:9.2f} um   (book: 1 is the '
          f'tighter)')

    # The derivation closes: the Nyquist rule on the local spatial frequency of
    # Eq. (7.51) at x1 = D1/2 gives back the constraint-3 upper bound.
    f_edge = local_spatial_frequency_source(D1a / 2.0, lama, za, m=hi3 / d1a)
    assert abs(abs(float(f_edge)) / (1.0 / (2.0 * d1a)) - 1.0) < 1e-12, f_edge
    # The same rule on Eq. (7.57) at the edge of the frequency grid gives back
    # constraint 4. The frequency-domain Nyquist limit is N delta1 / 2.
    f_t = float(local_spatial_frequency_transfer(1.0 / (2.0 * d1a), lama, za,
                                                 d1a, d2a))
    assert abs(f_t - lama * za / (2.0 * d2a)) < 1e-18
    assert abs(f_t / (d1a / 2.0) / n4 - 1.0) < 1e-12         # -> Eq. (7.59)
    print('  Nyquist on Eqs. (7.51) and (7.57) reproduces constraints 3 and 4')

    # The band-limit chain of Sec. 7.1 and 7.2: at the constraint-1 boundary
    # the geometric ray angle of Eq. (7.12) equals the Nyquist angle of the
    # SOURCE grid spacing. That equality is Eq. (7.13), printed p. 118.
    theta_geo = geometric_max_angle(D1a, D2a, d1a, d2_max, za)
    theta_nyq = nyquist_max_angle(lama, d1a)
    assert abs(theta_geo / theta_nyq - 1.0) < 1e-12, (theta_geo, theta_nyq)
    print(f'  Eq. (7.12) angle = Eq. (7.7) angle at the constraint-1 edge: '
          f'{theta_geo * 1e3:.4f} mrad')

    # The illuminated area drives constraint 2 through Eq. (7.18).
    d_illum = illuminated_diameter(D1a, d1a, d2a, lama, za)
    assert abs(((d_illum + D2a) / 2.0) / d2a - n2) / n2 < 1e-12
    print(f'  D_illum, Eq. (7.16)   {d_illum * 1e3:9.3f} mm  -> constraint 2')

    # ---------------- the book's worked example, Ch. 8.4 ----------------
    # A plane wave from a 2 mm square aperture, 2 m to the sensor.
    # Source: Ch. 8, Eq. (8.25) and text, printed p. 144.
    D1b, D2b, lamb, zb = 2e-3, 4e-3, 1e-6, 2.0
    d1b, dnb = 66.7e-6, 133e-6                       # 30 points across each
    Nb = 128

    dz_max = partial_max_step(Nb, d1b, dnb, lamb)
    assert abs(dz_max / 0.567 - 1.0) < 0.01, dz_max  # the book says 0.567 m
    n_planes = partial_plane_count(zb, Nb, d1b, dnb, lamb)
    assert n_planes == 5, n_planes                   # the book says 5

    # Table 8.2: the spacing is linear in the fractional distance.
    alphas = np.linspace(0.0, 1.0, n_planes)
    spac = partial_grid_spacing(d1b, dnb, alphas)
    assert abs(spac[0] - d1b) < 1e-18 and abs(spac[-1] - dnb) < 1e-18
    assert np.all(np.diff(spac) > 0.0)
    # Ch. 8, Eq. (8.22): constraint 3 is unchanged by the chain.
    lo8, hi8 = constraint3_delta2_window(D1b, d1b, lamb, zb)
    assert lo8 <= dnb <= hi8

    n2b = constraint2_min_n(D1b, D2b, d1b, dnb, lamb, zb)
    sigma = absorbing_boundary_sigma(Nb)
    assert abs(sigma - 0.47 * Nb) < 1e-12

    print('')
    print('BOOK Ch. 8, Eq. (8.25), printed p. 144 (partial propagations):')
    print(f'  delta1                {d1b * 1e6:9.1f} um   (book 66.7 um)')
    print(f'  delta_n               {dnb * 1e6:9.1f} um   (book 133 um)')
    print(f'  Delta_z max           {dz_max:9.3f} m    (book 0.567 m)')
    print(f'  planes n              {n_planes:9d}      (book 5)')
    print(f'  N_min, constraint 2   {n2b:9.1f}      (log2 = '
          f'{math.log2(n2b):.2f}; the book reads 2^7 off Fig. 8.5)')
    print(f'  absorber sigma        {sigma:9.1f} px   (book 0.47 N)')

    # ---------------- Table 6.2, the two-step planes ----------------
    # Source: Ch. 6, Table 6.2, printed p. 95.
    zt = 1.0
    for m_val, want in ((2.0, {'plus': (1 / 3, 2 / 3), 'minus': (-1.0, 2.0)}),
                        (0.5, {'plus': (2 / 3, 1 / 3), 'minus': (2.0, -1.0)})):
        got = two_step_planes(zt, m_val)
        for branch, pair in want.items():
            assert abs(got[branch][0] - pair[0]) < 1e-12, (m_val, branch, got)
            assert abs(got[branch][1] - pair[1]) < 1e-12, (m_val, branch, got)
            assert abs(sum(got[branch]) - zt) < 1e-12
    half = two_step_planes(zt, 1.0)
    assert abs(half['plus'][0] - 0.5) < 1e-12
    assert math.isinf(half['minus'][0])              # the book: infinitely far
    print('')
    print('BOOK Ch. 6, Table 6.2, printed p. 95 (two-step planes, z = 1):')
    for m_val in (2.0, 1.0, 0.5):
        p = two_step_planes(zt, m_val)
        print(f'  m = {m_val:4.2f}  plus  ({p["plus"][0]:+.4f}, '
              f'{p["plus"][1]:+.4f})   minus ({p["minus"][0]:+.4f}, '
              f'{p["minus"][1]:+.4f})')

    # The angular-spectrum m = 1 case is constraint 4 inverted.
    z_as = angular_spectrum_max_z(512, 10e-6, 1e-6)
    assert abs(constraint4_min_n(10e-6, 10e-6, 1e-6, z_as) - 512) < 1e-9
    print(f'  m = 1 range limit, N = 512, delta = 10 um: {z_as:.4f} m')

    # ---------------- the checker ----------------
    # A PASSING grid: the Ch. 7.3.2 angular-spectrum example.
    good = check_sampling(D1a, D2a, d1a, d2a, 512, lama, za)
    assert len(good) == 5
    for rule in good:
        assert rule.satisfied, rule

    # A FAILING grid: the same spacings taken 20 times as far.
    bad = check_sampling(D1a, D2a, d1a, d2a, 512, lama, 20.0 * za)
    assert not bad[1].satisfied and not bad[3].satisfied
    assert bad[0].satisfied                          # geometry row 1 relaxes
    assert all(BOOK in rule.citation for rule in good + bad)

    print('')
    print('check_sampling, the Ch. 7.3.2 grid (N = 512, z = 0.1 m):')
    for rule in good:
        mark = 'OK  ' if rule.satisfied else 'FAIL'
        bound = (f'({rule.bound[0]:.3e}, {rule.bound[1]:.3e})'
                 if isinstance(rule.bound, tuple) else f'{rule.bound:.4g}')
        print(f'  {mark} {rule.name:<42s} bound {bound:>24s}  '
              f'actual {rule.actual:.4g}')
    print('check_sampling, the same grid at z = 2 m:')
    for rule in bad:
        mark = 'OK  ' if rule.satisfied else 'FAIL'
        bound = (f'({rule.bound[0]:.3e}, {rule.bound[1]:.3e})'
                 if isinstance(rule.bound, tuple) else f'{rule.bound:.4g}')
        print(f'  {mark} {rule.name:<42s} bound {bound:>24s}  '
              f'actual {rule.actual:.4g}')

    # A geometry with no valid N: lambda z <= D2 delta1.
    assert math.isinf(one_step_min_n(2e-3, 3e-3, 1e-3, 1e-6, 0.5))
    # A converging source that no distance samples.
    assert math.isinf(fresnel_min_distance(2e-3, 40e-6, 1e-6, R=0.01))
    try:
        two_step_planes(1.0, -1.0)
    except ValueError:
        pass
    else:
        raise AssertionError('two_step_planes must refuse a negative m')

    print('')
    print('self-check passed')
