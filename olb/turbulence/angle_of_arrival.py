'''
Received tip-tilt (angle of arrival) of a Gaussian beam over a turbulent path.

This module gives the received tip-tilt angle variance of a Gaussian beam. The
tip-tilt is the random slope of the arriving wavefront. A receive telescope
focuses the beam onto a fibre tip. A tip-tilt of angle theta moves the focal
spot by f*theta, with f the focal length. So the received tip-tilt drives the
fibre-coupling loss. The Term factories live in olb.models.coupling.

The module is pure physics. It imports numpy, olb.turbulence.coupled_flux, the
sibling andrews layer, and the olb assumptions decorator layer. It does not
import the scenario, the terminal, or the results.

Each public function declares its own validity through the `@assumes` decorator
(olb.assumptions). The beam-wander arrival tilt is a radial (two-axis) variance
(`RADIAL_TILT`). The aperture angle-of-arrival tilt is the Andrews gradient tilt,
NOT the Noll Zernike tilt (`TILT_G_TILT`, Conflict C-04), and it holds only when
the Fresnel zone is small against the aperture (`FRESNEL_ZONE`). The Fresnel
condition needs the path length L, which `aperture_arrival_angle_variance` folds
into the Fried parameter r0; the RUNTIME gate on that condition therefore lives on
the delegate `olb.turbulence.andrews.structure.angle_of_arrival_variance`, which
takes the path length z directly. Outside a collection context the decorator is a
no-op, so the numeric output does not change.

Two contributions:

  A. Beam-wander arrival tilt (the DOMINANT term, and the one this module gives).
     The turbulence moves the beam centroid at the receiver by a random offset
     r_c. That offset is an apparent tilt r_c/L of the arriving beam, with L the
     path length. The received radial (2-axis) tilt variance is
         sigma2_theta = <r_c^2> / L^2.
     The kernel beam_wander_variance integrates the free-space beam WIDTH profile
     w(z) along the path, so the result is Gaussian-beam-correct.
     Source: Dios et al., Applied Optics 43 (2004) 3866. DOI 10.1364/AO.43.003866.

  B. Aperture angle-of-arrival "corrugation" tilt (a second, smaller term).
     The wavefront that arrives at the receive aperture is corrugated, so its
     mean slope across the pupil is not zero. The per-axis variance of that
     slope is
         sigma2_theta = 2.91 Cn2 L D^(-1/3) = 0.174 (D/r0)^(5/3)(lambda/D)^2.
     Source: Andrews and Phillips, Laser Beam Propagation through Random Media,
     2nd ed. (2005), DOI 10.1117/3.626196, Ch. 6, Eq. (84), printed p. 201.
     This is the GRADIENT tilt of a centroid tracker, not the Noll Zernike
     tilt. See aperture_arrival_angle_variance below.
'''

import numpy as np

from .coupled_flux import beam_wander_variance
from .andrews.structure import (
    angle_of_arrival_variance as _andrews_angle_of_arrival_variance,
)
from ..assumptions import (assumes, Constraint, BEAM_GAUSSIAN, REGIME_WEAK,
                           SPECTRUM_KOLMOGOROV)

# ----------------------------------------------------------------------------
# The module assumptions, as shared Constraint instances (see the docstring).
# ----------------------------------------------------------------------------

# The received beam-wander tip-tilt is a radial (two-axis) variance.
RADIAL_TILT = Constraint(
    "variance-convention",
    "The returned beam-wander tip-tilt variance is radial (two-axis). The "
    "per-axis variance is one half.",
    "10.1364/AO.43.003866", "beam-wander offset variance, Eq. (11), printed "
    "p. 3868")

# The aperture angle-of-arrival tilt is the Andrews gradient tilt (G-tilt).
TILT_G_TILT = Constraint(
    "tilt-convention",
    "The returned tilt is the Andrews gradient tilt (G-tilt), what a centroid "
    "tracker measures. It is NOT the Noll Zernike tilt. See Conflict C-04.",
    "10.1117/3.626196",
    "Ch. 6, Eq. (84), printed p. 201; definition Eq. (82), printed p. 200")

# The gradient-tilt result holds only in the small-Fresnel-zone limit. This
# signature folds the path length into r0, so the runtime gate on the condition
# lives on the delegate structure.angle_of_arrival_variance (which takes z).
FRESNEL_ZONE = Constraint(
    "field-region",
    "The gradient-tilt aperture angle of arrival holds only when the Fresnel "
    "zone is small against the aperture, sqrt(L/k) << D. The runtime gate is on "
    "the delegate structure.angle_of_arrival_variance, which takes the path "
    "length.",
    "10.1117/3.626196", "Ch. 6, text below Eq. (83), printed p. 200")


@assumes(RADIAL_TILT, beam_type=BEAM_GAUSSIAN, turbulence_regime=REGIME_WEAK,
         spectrum=SPECTRUM_KOLMOGOROV)
def wander_arrival_angle_variance(L, cn2_slant, w_profile, hs, *,
                                  wavelength=None):
    '''
    Return the received beam-wander tip-tilt variance (radial, 2-axis) [rad^2].

    Reuse the beam-wander kernel. It integrates the free-space beam width profile
    w(z) along the path, so it is Gaussian-beam-correct. The beam-wander offset
    variance <r_c^2> maps to an apparent arrival tilt through the path length L:
        sigma2_theta = <r_c^2> / L^2.
    The result is the radial (2-axis) variance. The per-axis variance is one half
    of it. Source: Dios et al., Applied Optics 43 (2004) 3866, DOI
    10.1364/AO.43.003866 (the beam-wander offset variance).

    Parameters:
        L : float
            Path length from the transmitter to the receiver [m].
        cn2_slant : numpy.ndarray
            Cn2 along the path on the hs grid [m^-2/3].
        w_profile : numpy.ndarray
            Free-space beam radius w(z) along the path on the hs grid [m].
        hs : numpy.ndarray
            Distance along the path from the transmitter [m].
        wavelength : float, optional
            Optical wavelength [m]. It does NOT change the result. It passes on
            to the wander kernel and turns ON its weak-regime runtime check, so
            a strong path gives a traced violation. None leaves the check off.

    Returns:
        float
            Radial (2-axis) received tip-tilt variance [rad^2].
    '''
    r_c2 = beam_wander_variance(L, cn2_slant, w_profile, hs,
                                wavelength=wavelength)
    return float(np.squeeze(r_c2)) / float(L) ** 2


@assumes(TILT_G_TILT, FRESNEL_ZONE, beam_type=BEAM_GAUSSIAN,
         turbulence_regime=REGIME_WEAK, spectrum=SPECTRUM_KOLMOGOROV)
def aperture_arrival_angle_variance(D, r0, wavelength):
    '''
    Return the aperture angle-of-arrival tip-tilt variance (per axis) [rad^2].

    This is the second, smaller aperture angle-of-arrival "corrugation" tilt. It
    is separate from the beam-wander arrival tilt above.

    TILT DEFINITION - THE OWNER MADE THIS CHOICE. This function returns the
    ANDREWS GRADIENT TILT (G-tilt), which is what a centroid tracker measures.
    Andrews defines the tilt as the total phase difference across the pupil
    divided by the pupil width. It is NOT the Noll Zernike tilt.

    formula:
        <beta_a^2> = 2.91 Cn2 L (2 W_G)^(-1/3)
                   = 0.174 (D/r0)^(5/3) (lambda/D)^2      per axis
    Source: Andrews and Phillips, Laser Beam Propagation through Random Media,
    2nd ed. (2005), DOI 10.1117/3.626196, Ch. 6, Eq. (84), printed p. 201, with
    the definition Ch. 6, Eq. (82), printed p. 200. The slant-path version is
    Ch. 12, Eq. (28), printed p. 492.

    THE RECAST. Put Cn2 L = r0^(-5/3)/(0.423 k^2) and k = 2*pi/lambda into
    Eq. (84). Then 2.91/(0.423 * 4 pi^2) = 0.1743.

    THE ALTERNATIVE. The Noll Zernike tilt gives 0.182 (D/r0)^(5/3)(lambda/D)^2
    (Noll, JOSA 66 (1976) 207, DOI 10.1364/JOSA.66.000207). A full-book search
    finds no 0.182 in Andrews and Phillips. Note that `olb/turbulence/ao.py`
    uses the NOLL convention (1.0299 and 0.134), so a caller that mixes the two
    must say which tilt it means. See Conflict C-04 in
    docs/andrews-crosscheck.md.

    Parameters:
        D : float or numpy.ndarray
            Receive aperture diameter [m]. The book writes it as 2 W_G.
        r0 : float or numpy.ndarray
            Fried atmospheric coherence width over the same path [m].
        wavelength : float
            Optical wavelength [m].

    Returns:
        float or numpy.ndarray
            PER-AXIS tilt variance [rad^2]. The radial (2-axis) variance is
            twice this value.

    VALIDITY. The result holds only when the Fresnel zone is small against the
    aperture, sqrt(L/k) << D (Ch. 6, text below Eq. (83), printed p. 200). This
    function does not gate on that condition.

    NEW HOME: `olb.turbulence.andrews.structure.angle_of_arrival_variance`. That
    function takes Cn2 and the path length directly, and it also gives the
    inner-scale and outer-scale branches of Ch. 6, Eq. (83).
    '''
    # Rebuild the path moment Cn2 * L from the Fried parameter, so that this
    # signature stays unchanged. r0 = (0.423 k^2 Cn2 L)^(-3/5), so
    # Cn2 L = r0^(-5/3) / (0.423 k^2). The Andrews function uses only the
    # product Cn2 * z, so z = 1 m carries it.
    k = 2.0 * np.pi / wavelength
    cn2_l = np.asarray(r0, dtype=float) ** (-5.0 / 3.0) / (0.423 * k ** 2)
    return _andrews_angle_of_arrival_variance(D, wavelength, 1.0, cn2_l)


if __name__ == '__main__':
    # Pure-physics self-check. Use plain numeric inputs. No scenario import.
    from ..beam import free_space_radius

    lam = 1550e-9
    w0 = 0.02
    L = 3e3
    hs = np.linspace(0.0, L, 200)

    def _variance(cn2, length=L, waist=w0):
        grid = np.linspace(0.0, length, 200)
        cn2_slant = np.full_like(grid, cn2)
        w_profile = free_space_radius(waist, grid, None, lam)
        return wander_arrival_angle_variance(length, cn2_slant, w_profile, grid)

    # A real, positive tilt variance for a turbulent path.
    v = _variance(1e-14)
    assert np.isfinite(v) and v > 0.0, v

    # Stronger turbulence gives a larger tilt variance.
    assert _variance(1e-13) > _variance(1e-15)

    # The tilt variance scales linearly with Cn2 (beam wander is linear in Cn2).
    ratio = _variance(2e-14) / _variance(1e-14)
    assert np.isclose(ratio, 2.0, rtol=1e-6), ratio

    # The aperture angle-of-arrival term. It falls as D^(-1/3) and it grows as
    # the Fried parameter falls.
    D_ap, r0_ap = 0.2, 0.1
    v_ap = aperture_arrival_angle_variance(D_ap, r0_ap, lam)
    assert v_ap > 0.0, v_ap
    assert aperture_arrival_angle_variance(0.4, r0_ap, lam) < v_ap
    assert aperture_arrival_angle_variance(D_ap, 0.05, lam) > v_ap
    # The gradient-tilt recast, Andrews Ch. 6, Eq. (84), printed p. 201.
    recast = 0.174 * (D_ap / r0_ap) ** (5.0 / 3.0) * (lam / D_ap) ** 2
    pct = abs(v_ap - recast) / recast * 100.0
    assert pct < 2.0, pct

    print(f"wander tilt variance (3 km, Cn2=1e-14) = {v:.3e} rad^2")
    print(f"aperture AoA tilt (D=0.2 m, r0=0.1 m) = {v_ap:.3e} rad^2  "
          f"per-axis 1-sigma = {np.sqrt(v_ap) * 1e6:.3f} urad  "
          f"(gradient tilt, {pct:.3f} % from the 0.174 recast)")
    print(f"  radial 1-sigma = {np.sqrt(v) * 1e6:.3f} urad  "
          f"per-axis 1-sigma = {np.sqrt(v / 2) * 1e6:.3f} urad")

    # ---------------- assumption self-checks ----------------
    import warnings
    from ..assumptions import trace_assumptions

    grid = np.linspace(0.0, L, 200)
    cn2_grid = np.full_like(grid, 1e-14)
    w_grid = free_space_radius(w0, grid, None, lam)

    # (1) Value parity: one representative call returns the identical float with
    #     and without a collection context.
    ap_out = aperture_arrival_angle_variance(0.2, 0.1, lam)
    with trace_assumptions():
        ap_in = aperture_arrival_angle_variance(0.2, 0.1, lam)
    assert ap_out == ap_in, (ap_out, ap_in)

    # (2) Registration: inside a context the expected sources and kinds register.
    #     The wander call also registers the decorated coupled_flux dependency, so
    #     the C-01 conflict tag is inherited automatically.
    with trace_assumptions() as tr:
        wander_arrival_angle_variance(L, cn2_grid, w_grid, grid)
        aperture_arrival_angle_variance(0.2, 0.1, lam)
    mod = __name__
    assert f"{mod}.wander_arrival_angle_variance" in tr.records
    assert f"{mod}.aperture_arrival_angle_variance" in tr.records
    assert "olb.turbulence.coupled_flux.beam_wander_variance" in tr.records, \
        "the decorated dependency must register through the wander call"
    ap_rec = tr.records[f"{mod}.aperture_arrival_angle_variance"]
    ap_kinds = {c.kind for c in ap_rec.constraints}
    assert {"tilt-convention", "field-region"} <= ap_kinds, ap_kinds
    wander_kinds = {c.kind for c in
                    tr.records[f"{mod}.wander_arrival_angle_variance"].constraints}
    assert "variance-convention" in wander_kinds, wander_kinds

    # (2b) The optional wavelength passes to the kernel and turns on its
    #      weak-regime check. It does NOT change the value.
    strong = np.full_like(grid, 1e-13)
    plain = wander_arrival_angle_variance(L, strong, w_grid, grid)
    with_lam = wander_arrival_angle_variance(L, strong, w_grid, grid,
                                             wavelength=lam)
    assert plain == with_lam, (plain, with_lam)
    with trace_assumptions() as tr_hard:
        wander_arrival_angle_variance(L, strong, w_grid, grid, wavelength=lam)
    assert any("beam-wander model is not trusted" in v
               for v in tr_hard.violations), tr_hard.violations

    # (3) No decorator check in this module warns. Neither validity condition is
    #     a runtime callable here: the tilt-convention is a labelling assumption,
    #     and the sqrt(L/k) << D gate needs the path length, which this signature
    #     folds into r0 (its runtime gate is on the delegate
    #     structure.angle_of_arrival_variance). So a deliberately large-Fresnel
    #     call registers the field-region CONSTRAINT (below) but raises no
    #     decorator violation and no warning. The firing-check demonstration for
    #     the traced machinery is in coupled_flux and uplink_flux.
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        with trace_assumptions() as tr_bad:
            aperture_arrival_angle_variance(1e-3, 0.1, lam)   # a tiny aperture
    assert any(c.kind == "field-region"
               for c in tr_bad.records[f"{mod}.aperture_arrival_angle_variance"]
               .constraints), "the sqrt(L/k) << D limit must be recorded"
    assert len(caught) == 0, "a decorator check must not warn"

    print("angle_of_arrival assumptions self-check passed")
    print("self-check passed")
