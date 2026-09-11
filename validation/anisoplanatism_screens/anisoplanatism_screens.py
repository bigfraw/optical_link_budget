"""Does the fidelity-2 screen plan hold the POINT-AHEAD anisoplanatism? No sim.

THE QUESTION. An uplink terminal senses the downlink beacon and it launches the
uplink one point-ahead angle theta away. In the plane-parallel screen model the
two paths cross a screen at the ground distance z_g with the lateral offset
d = theta z_g. So a split-step run can model the point-ahead decorrelation: it
reads each screen two times, one time at the beacon window and one time at the
window that the offset d selects. THE SCREEN PLAN must then resolve the HEIGHT
structure of the anisoplanatic error, and the plan of olb cuts the slab at equal
RYTOV weight, not at equal anisoplanatic weight. This script measures the error
of that choice. It runs no propagation and it draws no screen.

THE METHOD. Stone et al. (1994) give the anisoplanatic phase variance of a
finite aperture as a HEIGHT INTEGRAL, Eq. (29):

    sigma^2 = 2 (2 pi)^(8/3) C_A k0^2 R^(5/3) airmass
              * INT dh Cn2(h) I( S theta / R ),   S = h airmass,  R = D / 2

`olb.turbulence.anisoplanatism.anisoplanatic_phase_variance` computes it. A
screen plan replaces that integral by a DELTA-LAYER SUM: each screen holds the
integrated Cn2 of its slab, at one height. The script compares the sum with the
integral. A ratio of 1.0 says the plan resolves the height structure.

THE SLANT FACTOR (verified here, see `discrete_stone`). `ScreenPlan.cn2_int_m13`
already carries the airmass: `_plan_space_continuous` integrates the density
m(h) = Cn2(h) sec(zeta) over the height. Eq. (29) carries one airmass outside
the height integral, which is the same factor, because airmass dh is the SLANT
path element. So the delta-layer sum needs NO extra airmass.

THE PIXEL. A split-step screen is a grid, so a real run can only offset a window
by a WHOLE pixel. The script also reports the variance of the pixel-ROUNDED
offset round(theta z_g / dx) dx at the production pixel.

THE PLANS. The production plan (9 screens at the standard preset), the override
plans of 5, 9, 15 and 25 screens, the ground-split plan (the lowest screen cut
into 4), and the EQUAL-ANISOPLANATIC-WEIGHT plans of the same counts. The last
family cuts the slab at equal shares of the Stone weight
Cn2(h) sec (h sec)^p, so it puts the screens where the anisoplanatism is, not
where the scintillation is.

Sources:
- J. Stone, P. H. Hu, S. P. Mills and S. Ma, "Anisoplanatic effects in
  finite-aperture optical systems," J. Opt. Soc. Am. A 11(1), 347-357 (1994),
  DOI 10.1364/JOSAA.11.000347. Eqs. (29), (36) and (A11).
- R. J. Noll, J. Opt. Soc. Am. 66(3), 207 (1976), DOI 10.1364/JOSA.66.000207.
  The Noll mode count that `max_radial_order` maps to a radial order.
- Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196. Ch. 12, Eq. (14):
  the plane-parallel airmass. Ch. 8, Eq. (20): the Rytov path weight of the
  production planner.

Run it from the repository root:

    python -m validation.anisoplanatism_screens.anisoplanatism_screens
"""

import argparse
import time
import warnings

import numpy as np

from olb.turbulence.anisoplanatism import (C_A, _REMOVE_NLO, _TWO_PI_83,
                                           _inner_integral,
                                           anisoplanatic_phase_variance,
                                           max_radial_order)
from olb.turbulence.andrews.beam import wavenumber
from olb.turbulence.andrews.paths import sec_zeta
from olb.waveoptics.turbulence.sampling import (DEFAULT_H_TOP_M,
                                                _integration_heights)
from validation.anisoplanatism_screens.common import (ANGLES_ARCSEC, ARCSEC,
                                                     APERTURES_M,
                                                     ELEVATIONS_DEG, LAM,
                                                     equal_aniso_plan,
                                                     fig_dir,
                                                     ground_distances,
                                                     hero_uplink, log_maker,
                                                     plan_set, pyplot,
                                                     site_cn2, write_results)

# The override screen counts, and the sub-screen count of the ground split.
COUNTS = (5, 9, 15, 25)
SPLIT_N = 4

# The adaptive-optics bands under test. None is the ideal, infinite-order
# correction. 3 and 5 are the highest complete radial orders of 10 and 21 Noll
# modes (max_radial_order; Noll, DOI 10.1364/JOSA.66.000207, Table I).
MAX_ORDERS = (max_radial_order(10), max_radial_order(21), None)

# The verdict tiers on the ratio of the discrete sum to the continuous integral.
PASS_TOL = 0.05
FAIL_TOL = 0.10


class InnerI:
    """A fast, cached reader of the Stone spatial-frequency integral I(beta).

    THE POINT. `_inner_integral` costs three `scipy.integrate.quad` calls of an
    oscillating Bessel integrand, so a sweep that asks for it thousands of times
    is slow. I(beta) is smooth and monotone, so this class tabulates it on a
    log-spaced beta grid and it reads the table with a log-log interpolation.
    Below the lowest node and above the highest one it extrapolates with the
    local power law, because I(beta) is a power of beta in both limits (see the
    docstring of `olb.turbulence.anisoplanatism._inner_integral`).

    The physics is NOT duplicated: every node is one `_inner_integral` call.
    `check` measures the interpolation error against fresh calls.

    Attributes:
        n_lo:      the lowest radial order that counts as error.
        max_order: the highest corrected radial order, or None.
    """

    def __init__(self, n_lo, max_order, n_nodes=161, b_lo=1e-4, b_hi=1e3):
        """Build the table of one (n_lo, max_order) pair."""
        self.n_lo = int(n_lo)
        self.max_order = max_order
        self.b = np.geomspace(float(b_lo), float(b_hi), int(n_nodes))
        self.v = np.array([_inner_integral(b, self.n_lo, self.max_order)
                           for b in self.b])
        self._lb = np.log(self.b)
        self._lv = np.log(self.v)
        # The local power law of each end, for the extrapolation.
        self._p_lo = float((self._lv[1] - self._lv[0])
                           / (self._lb[1] - self._lb[0]))
        self._p_hi = float((self._lv[-1] - self._lv[-2])
                           / (self._lb[-1] - self._lb[-2]))

    def __call__(self, beta):
        """Give I(beta). beta may be an array. I(0) is 0."""
        b = np.atleast_1d(np.asarray(beta, dtype=float))
        out = np.zeros_like(b)
        good = b > 0.0
        lb = np.log(np.where(good, b, 1.0))
        lv = np.interp(lb, self._lb, self._lv)
        low = good & (b < self.b[0])
        high = good & (b > self.b[-1])
        lv = np.where(low, self._lv[0] + self._p_lo * (lb - self._lb[0]), lv)
        lv = np.where(high, self._lv[-1] + self._p_hi * (lb - self._lb[-1]), lv)
        out[good] = np.exp(lv[good])
        return out if np.ndim(beta) else float(out[0])

    @property
    def small_beta_exponent(self):
        """The power law of I(beta) at a small beta. It sets the Stone weight."""
        return self._p_lo

    def check(self, betas):
        """Give the largest relative error of the table against fresh quad."""
        err = 0.0
        for b in betas:
            exact = _inner_integral(float(b), self.n_lo, self.max_order)
            if exact > 0.0:
                err = max(err, abs(self(float(b)) / exact - 1.0))
        return float(err)


def _inner_cache():
    """Give the InnerI table of each (remove, max_order) pair that runs."""
    cache = {}

    def get(remove, max_order):
        key = (remove, max_order)
        if key not in cache:
            cache[key] = InnerI(_REMOVE_NLO[remove], max_order)
        return cache[key]

    return get


def stone_prefactor(D, lam):
    """Give the prefactor of Stone Eq. (29), WITHOUT the airmass.

        2 (2 pi)^(8/3) C_A k0^2 R^(5/3),   R = D / 2

    Source: Stone et al. (1994), DOI 10.1364/JOSAA.11.000347, Eq. (29). The
    airmass stays out, because `ScreenPlan.cn2_int_m13` and the slant integral
    of `continuous_stone` carry it.
    """
    k0 = wavenumber(lam)
    return 2.0 * _TWO_PI_83 * C_A * k0 ** 2 * (D / 2.0) ** (5.0 / 3.0)


def discrete_stone(plan, D, theta, lam, remove="piston_tilt", max_order=None,
                   dx=None, inner=None):
    """Give the delta-layer Stone variance of one screen plan.

        sigma^2 = 2 (2 pi)^(8/3) C_A k0^2 R^(5/3)
                  * SUM_j cn2_int_m13[j] I( d_j / R )

    with d_j = theta z_g_j the lateral offset of screen j, and
    z_g_j = z_total - z_m[j] the GROUND distance of that screen.

    THE SLANT FACTOR IS ALREADY IN THE PLAN. Eq. (29) of Stone et al. (1994),
    DOI 10.1364/JOSAA.11.000347, reads

        airmass * INT dh Cn2(h) I(S theta / R),   S = h airmass,

    and `airmass dh` is the slant path element ds. The plan integrates the
    density m(h) = Cn2(h) sec(zeta) over the height
    (`_plan_space_continuous`), which is the same product. The Stone height
    argument S = h airmass is the ground distance z_g of the screen, because
    z_total - z_m = (h_top - (h_top - h)) sec = h sec. So the sum above needs
    no extra airmass and no extra sec.

    THE PIXEL. A split-step screen is a grid, so a run can shift a window by a
    whole pixel only. A `dx` rounds each offset to round(d_j / dx) dx.

    Args:
        plan:      a ScreenPlan of a space link (direction "down").
        D:         the aperture diameter, in m.
        theta:     the point-ahead angle, in rad.
        lam:       the wavelength, in m.
        remove:    "none", "piston" or "piston_tilt".
        max_order: the highest corrected radial order, or None.
        dx:        the pixel pitch, in m, or None for the exact offset.
        inner:     an InnerI of the same (remove, max_order), or None to build
                   one.

    Returns:
        A dict with the keys sigma2, contributions, shift_m, shift_px and beta.
    """
    if inner is None:
        inner = InnerI(_REMOVE_NLO[remove], max_order)
    zg = ground_distances(plan)
    shift = zg * float(theta)
    shift_px = shift / float(dx) if dx else None
    if dx:
        shift = np.round(shift / float(dx)) * float(dx)
    beta = shift / (D / 2.0)
    contrib = (stone_prefactor(D, lam)
               * np.asarray(plan.cn2_int_m13, dtype=float) * inner(beta))
    return {"sigma2": float(contrib.sum()), "contributions": contrib,
            "shift_m": shift, "beta": beta,
            "shift_px": shift_px if shift_px is not None else shift * 0.0}


def continuous_stone(scn, elevation_deg, D, theta, lam,
                     remove="piston_tilt", max_order=None, inner=None):
    """Give the continuous Stone variance of the site profile.

    It is Eq. (29) of Stone et al. (1994), DOI 10.1364/JOSAA.11.000347, on the
    SAME internal height grid that the planner integrates
    (`_integration_heights`), with the SAME site Hufnagel-Valley Cn2 callable.
    `check_continuous` measures it against
    `olb.turbulence.anisoplanatism.anisoplanatic_phase_variance`.

    Args:
        scn:           the SpaceScenario.
        elevation_deg: the elevation, in deg.
        D:             the aperture diameter, in m.
        theta:         the point-ahead angle, in rad.
        lam:           the wavelength, in m.
        remove:        "none", "piston" or "piston_tilt".
        max_order:     the highest corrected radial order, or None.
        inner:         an InnerI of the same pair, or None.

    Returns:
        The variance, in rad^2.
    """
    if inner is None:
        inner = InnerI(_REMOVE_NLO[remove], max_order)
    sec = float(sec_zeta(float(elevation_deg)))
    h = _integration_heights(DEFAULT_H_TOP_M)
    cn2_h = np.asarray(site_cn2(scn)(h), dtype=float)
    beta = h * sec * float(theta) / (D / 2.0)
    return float(stone_prefactor(D, lam) * sec
                 * np.trapezoid(cn2_h * inner(beta), h))


def check_continuous(scn, elevation_deg, lam, get_inner):
    """Compare `continuous_stone` with the production library function.

    The library call integrates the SAME equation on the same height grid, but
    it evaluates `_inner_integral` at every height. So the comparison measures
    the interpolation of InnerI and nothing else.

    Returns:
        A list of dicts, one for each checked cell.
    """
    h = _integration_heights(DEFAULT_H_TOP_M)
    cn2_h = np.asarray(site_cn2(scn)(h), dtype=float)
    out = []
    for D, arcsec, max_order in ((0.7, 10.0, None), (0.4, 2.0, 3)):
        theta = arcsec * ARCSEC
        mine = continuous_stone(scn, elevation_deg, D, theta, lam,
                                max_order=max_order,
                                inner=get_inner("piston_tilt", max_order))
        lib = anisoplanatic_phase_variance(
            D, theta, h, cn2_h, lam, remove="piston_tilt",
            max_order=max_order, elevation_deg=float(elevation_deg))
        out.append({"D_m": D, "theta_arcsec": arcsec,
                    "max_order": max_order, "mine": mine, "library": lib,
                    "ratio": mine / lib if lib else float("nan")})
    return out


def plans_of(scn, geom, exponent):
    """Give every plan of one elevation: production, override, split, aniso.

    THE TWO ANISOPLANATIC FAMILIES. `aniso<n>` cuts the slab at equal
    anisoplanatic weight and it keeps the production Cn2-weighted centroid.
    `anisoc<n>` cuts the same way and it moves the screen to the
    anisoplanatic-weighted distance, so a delta layer reproduces the
    small-angle Stone integral of its slab.
    """
    grid, plans = plan_set(scn, geom, COUNTS, split_n=SPLIT_N)
    for n in COUNTS:
        plans[f"aniso{int(n)}"] = equal_aniso_plan(scn, geom, n, exponent)
        plans[f"anisoc{int(n)}"] = equal_aniso_plan(scn, geom, n, exponent,
                                                    centroid="aniso")
    return grid, plans


def sweep(say, args):
    """Run the whole sweep. Give the JSON payload."""
    get_inner = _inner_cache()
    lam = LAM

    # The interpolation guard. It measures InnerI against fresh quad calls.
    guard = []
    for max_order in MAX_ORDERS:
        inner = get_inner("piston_tilt", max_order)
        err = inner.check((3e-4, 0.01, 0.3, 1.7, 12.0, 90.0))
        guard.append({"max_order": max_order, "max_rel_error": err,
                      "small_beta_exponent": inner.small_beta_exponent})
    say("THE INTERPOLATION GUARD of the Stone integral table I(beta)")
    for g in guard:
        say(f"  max_order {str(g['max_order']):>4s}: largest relative error "
            f"{g['max_rel_error']:.2e}; I(beta) goes as beta^"
            f"{g['small_beta_exponent']:.3f} at a small beta")
    exponent = float(guard[-1]["small_beta_exponent"])
    say(f"  The piston-and-tilt-removed weight takes the exponent "
        f"{exponent:.2f}. The classical (remove='none') exponent is 5/3; the "
        f"equal-anisoplanatic-weight plan uses the exponent of the swept "
        f"`remove`.")
    say()

    # THE ORDER SANITY CHECK. Each removed mode takes variance away, so the
    # classical result (remove="none") must be the largest. Stone et al. (1994),
    # DOI 10.1364/JOSAA.11.000347, Fig. 1.
    scn30, geom30 = hero_uplink(30.0)
    _, plans30 = plan_set(scn30, geom30, COUNTS, split_n=SPLIT_N)
    order = {}
    for remove in ("none", "piston", "piston_tilt"):
        inner = get_inner(remove, None)
        order[remove] = {
            "continuous": continuous_stone(scn30, 30.0, 0.7, 10.0 * ARCSEC,
                                           lam, remove=remove, inner=inner),
            "discrete": discrete_stone(plans30["production"], 0.7,
                                       10.0 * ARCSEC, lam, remove=remove,
                                       inner=inner)["sigma2"]}
    say("THE ORDER SANITY CHECK at 30 deg, D = 0.7 m, 10 arcsec "
        "(continuous / discrete, rad^2)")
    for remove in ("none", "piston", "piston_tilt"):
        say(f"  remove = {remove:<12s} {order[remove]['continuous']:10.4g} / "
            f"{order[remove]['discrete']:10.4g}")
    assert (order["none"]["continuous"] > order["piston"]["continuous"]
            > order["piston_tilt"]["continuous"]), order
    say("  The order holds: every removed mode takes variance away.")
    say()

    rows = []
    per_screen = []
    checks = []
    for el in args.elevations:
        scn, geom = hero_uplink(el)
        grid, plans = plans_of(scn, geom, exponent)
        dx = float(grid.pixel_m)
        theta_geom = float(np.atleast_1d(geom.point_ahead_rad)[0])
        # The library check is SLOW (it calls quad at every height), so it runs
        # at one elevation only. It measures the InnerI table, and the table
        # does not depend on the elevation.
        if abs(el - 30.0) < 1e-9:
            checks.append({"elevation_deg": el,
                           "cells": check_continuous(scn, el, lam, get_inner)})
        angles = list(args.angles_arcsec) + [theta_geom / ARCSEC]
        say(f"ELEVATION {el:.0f} deg: grid {grid.n} px, pixel "
            f"{dx * 1e3:.3f} mm, slab {plans['production'].z_total_m / 1e3:.1f} "
            f"km, point-ahead {theta_geom / ARCSEC:.3f} arcsec")
        hdr = (f"  {'plan':<18s}{'n':>3s}{'D[m]':>6s}{'th[as]':>8s}"
               f"{'order':>6s}{'cont[rad2]':>12s}{'disc[rad2]':>12s}"
               f"{'ratio':>8s}{'px-round':>10s}")
        say(hdr)
        say("  " + "-" * (len(hdr) - 2))
        for name, plan in plans.items():
            for D in args.apertures:
                for arcsec in angles:
                    theta = arcsec * ARCSEC
                    for max_order in MAX_ORDERS:
                        inner = get_inner("piston_tilt", max_order)
                        cont = continuous_stone(scn, el, D, theta, lam,
                                                max_order=max_order,
                                                inner=inner)
                        d = discrete_stone(plan, D, theta, lam,
                                           max_order=max_order, inner=inner)
                        d_px = discrete_stone(plan, D, theta, lam,
                                              max_order=max_order, dx=dx,
                                              inner=inner)
                        ratio = d["sigma2"] / cont if cont > 0 else float("nan")
                        r_px = (d_px["sigma2"] / cont if cont > 0
                                else float("nan"))
                        rows.append({
                            "elevation_deg": el, "plan": name,
                            "n_screens": int(plan.z_m.size), "D_m": D,
                            "theta_arcsec": float(arcsec),
                            "is_geometry_angle": bool(
                                abs(arcsec - theta_geom / ARCSEC) < 1e-9),
                            "max_order": max_order,
                            "continuous_rad2": cont,
                            "discrete_rad2": d["sigma2"], "ratio": ratio,
                            "discrete_rounded_rad2": d_px["sigma2"],
                            "ratio_rounded": r_px,
                            "rounded_delta": r_px - ratio,
                            "grid_n": int(grid.n), "pixel_m": dx})
                        if args.verbose or (name == "production"
                                            and D == 0.7
                                            and max_order is None):
                            say(f"  {name:<18s}{plan.z_m.size:3d}{D:6.2f}"
                                f"{arcsec:8.2f}{str(max_order):>6s}"
                                f"{cont:12.4g}{d['sigma2']:12.4g}"
                                f"{ratio:8.3f}{r_px - ratio:+10.3f}")
        # The per-screen table of the production plan.
        if abs(el - 30.0) < 1e-9:
            per_screen = screen_table(say, scn, plans["production"], dx, lam,
                                      get_inner)
        say()
    return rows, guard, checks, per_screen, order


def screen_table(say, scn, plan, dx, lam, get_inner):
    """Print the per-screen contribution table of the production plan.

    The case is 30 deg, 10 arcsec, D = 0.7 m, max_order None. It answers the
    question whether the GROUND screens need a finer cut.
    """
    D, arcsec, max_order = 0.7, 10.0, None
    theta = arcsec * ARCSEC
    inner = get_inner("piston_tilt", max_order)
    d = discrete_stone(plan, D, theta, lam, max_order=max_order, inner=inner)
    zg = ground_distances(plan)
    sec = float(sec_zeta(30.0))
    heights = zg / sec
    share = d["contributions"] / d["contributions"].sum()
    say(f"  THE PER-SCREEN TABLE: 30 deg, {arcsec:.0f} arcsec, D = {D} m, "
        f"max_order None, the production plan")
    say(f"  {'j':>3s}{'height[km]':>12s}{'z_g[km]':>10s}"
        f"{'cn2_int':>11s}{'shift[m]':>10s}{'shift[px]':>11s}"
        f"{'share':>9s}")
    out = []
    for j in range(plan.z_m.size):
        say(f"  {j:3d}{heights[j] / 1e3:12.3f}{zg[j] / 1e3:10.3f}"
            f"{plan.cn2_int_m13[j]:11.3e}{zg[j] * theta:10.4f}"
            f"{zg[j] * theta / dx:11.1f}{share[j] * 100:8.2f}%")
        out.append({"j": j, "height_m": float(heights[j]),
                    "z_ground_m": float(zg[j]),
                    "cn2_int_m13": float(plan.cn2_int_m13[j]),
                    "shift_m": float(zg[j] * theta),
                    "shift_px": float(zg[j] * theta / dx),
                    "share": float(share[j])})
    say(f"  The two LOWEST screens (j = {plan.z_m.size - 2} and "
        f"{plan.z_m.size - 1}) hold {(share[-2] + share[-1]) * 100:.2f}% of "
        f"the variance. The TOP screen holds {share[0] * 100:.2f}%.")
    return out


def verdict(say, rows):
    """Print the verdict of the production plan, and give the payload."""
    owner = [r for r in rows
             if r["D_m"] <= 1.0 and r["theta_arcsec"] <= 10.0
             and r["elevation_deg"] >= 20.0]
    prod = [r for r in owner if r["plan"] == "production"]
    worst = max(prod, key=lambda r: abs(r["ratio"] - 1.0))
    by_el = {}
    for el in sorted({r["elevation_deg"] for r in prod}):
        cell = [r for r in prod if r["elevation_deg"] == el]
        w = max(cell, key=lambda r: abs(r["ratio"] - 1.0))
        by_el[el] = w
    say("THE VERDICT")
    say(f"  The owner regime holds {len(prod)} production cells "
         f"(D <= 1 m, theta <= 10 arcsec, elevation >= 20 deg).")
    for el, w in by_el.items():
        say(f"  {el:4.0f} deg: the worst ratio is {w['ratio']:.3f} at D = "
            f"{w['D_m']} m, {w['theta_arcsec']:.2f} arcsec, max_order "
            f"{w['max_order']}.")
    worst_dev = abs(worst["ratio"] - 1.0)
    if worst_dev <= PASS_TOL:
        text = ("NO PLANNER CHANGE. Every production cell of the owner regime "
                f"is inside {PASS_TOL * 100:.0f} percent of the continuous "
                "Stone integral.")
    elif worst_dev <= FAIL_TOL:
        text = (f"MARGINAL. The worst production cell is off by "
                f"{worst_dev * 100:.1f} percent, past the "
                f"{PASS_TOL * 100:.0f} percent pass band and inside the "
                f"{FAIL_TOL * 100:.0f} percent fail band.")
    else:
        text = (f"THE PLAN IS NOT SUFFICIENT. The worst production cell is off "
                f"by {worst_dev * 100:.1f} percent.")
    say(f"  {text}")
    say(f"  The worst cell: {worst['elevation_deg']:.0f} deg, D = "
        f"{worst['D_m']} m, {worst['theta_arcsec']:.2f} arcsec, max_order "
        f"{worst['max_order']}; discrete {worst['discrete_rad2']:.4g} rad^2 "
        f"against continuous {worst['continuous_rad2']:.4g} rad^2, ratio "
        f"{worst['ratio']:.3f}.")

    # The blended plans at the same count, on the same worst cell.
    def _same_cell(plan_name):
        sel = [r for r in rows
               if r["plan"] == plan_name
               and r["elevation_deg"] == worst["elevation_deg"]
               and r["D_m"] == worst["D_m"]
               and r["theta_arcsec"] == worst["theta_arcsec"]
               and r["max_order"] == worst["max_order"]]
        return sel[0] if sel else None

    ratios = [r["ratio"] for r in prod]
    say(f"  The production ratio holds between {min(ratios):.3f} and "
        f"{max(ratios):.3f} over those cells.")
    aniso = [r for r in owner if r["plan"].startswith("aniso")]
    if aniso:
        bad = max(aniso, key=lambda r: abs(r["ratio"] - 1.0))
        say(f"  The equal-anisoplanatic-weight families do NOT improve the "
            f"plan: their worst owner-regime ratio is {bad['ratio']:.3f} "
            f"({bad['plan']}, {bad['elevation_deg']:.0f} deg, D = "
            f"{bad['D_m']} m, {bad['theta_arcsec']:.2f} arcsec). A slab that "
            f"is thick in height breaks the small-angle limit that the "
            f"matched centroid assumes.")
    same = _same_cell("n9")
    fix = _same_cell("aniso9")
    fixc = _same_cell("anisoc9")
    if same and fix and fixc:
        say(f"  On that cell, 9 screens give the ratio "
            f"{same['ratio']:.3f} with the equal-Rytov cut, "
            f"{fix['ratio']:.3f} with the equal-anisoplanatic cut at the Cn2 "
            f"centroid, and {fixc['ratio']:.3f} with the equal-anisoplanatic "
            f"cut at the matched centroid.")
    round_prod = max(prod, key=lambda r: abs(r["rounded_delta"]))
    round_worst = max(owner, key=lambda r: abs(r["rounded_delta"]))
    say(f"  The pixel-rounded offset moves the ratio of the PRODUCTION plan by "
        f"at most {round_prod['rounded_delta']:+.3f} "
        f"({round_prod['elevation_deg']:.0f} deg, D = {round_prod['D_m']} m, "
        f"{round_prod['theta_arcsec']:.2f} arcsec).")
    say(f"  Over EVERY swept plan the largest move is "
        f"{round_worst['rounded_delta']:+.3f} ({round_worst['plan']}, "
        f"{round_worst['elevation_deg']:.0f} deg, D = {round_worst['D_m']} m, "
        f"{round_worst['theta_arcsec']:.2f} arcsec). A plan whose screens sit "
        f"high has a small offset in pixels on its lowest screen, so the "
        f"rounding hits it harder.")
    return {"text": text, "worst": worst,
            "worst_rounded_delta_production": round_prod,
            "worst_by_elevation": {str(k): v for k, v in by_el.items()},
            "aniso9_on_worst_cell": fix, "anisoc9_on_worst_cell": fixc,
            "n9_on_worst_cell": same,
            "worst_rounded_delta": round_worst}


def plot(rows):
    """Draw the ratio against the screen count, one line for each elevation.

    One panel for each plan family: the equal-Rytov cut (the production shape),
    the equal-anisoplanatic cut with the production centroid, and the
    equal-anisoplanatic cut with the anisoplanatic centroid.
    """
    plt = pyplot()
    if plt is None:
        return []
    panels = (("n", "equal Rytov weight (production shape)"),
              ("aniso", "equal anisoplanatic weight, Cn2 centroid"),
              ("anisoc", "equal anisoplanatic weight, matched centroid"))
    fig, axes = plt.subplots(1, len(panels), figsize=(13.0, 4.4), sharey=True)
    for ax, (family, title) in zip(axes, panels):
        for el in sorted({r["elevation_deg"] for r in rows}):
            xs, ys = [], []
            for n in COUNTS:
                sel = [r for r in rows
                       if r["plan"] == f"{family}{n}"
                       and r["elevation_deg"] == el and r["D_m"] == 0.7
                       and abs(r["theta_arcsec"] - 10.0) < 1e-9
                       and r["max_order"] is None]
                if sel:
                    xs.append(n)
                    ys.append(sel[0]["ratio"])
            if xs:
                ax.plot(xs, ys, "-o", ms=4, label=f"{el:.0f} deg")
        ax.axhline(1.0, color="k", lw=0.8)
        ax.axhspan(1.0 - PASS_TOL, 1.0 + PASS_TOL, color="0.85", zorder=0)
        ax.set_xlabel("screens in the plan")
        ax.set_title(title, fontsize=9)
        ax.grid(alpha=0.3)
    axes[0].set_ylabel("discrete sum / continuous Stone integral")
    axes[0].legend(fontsize=8)
    fig.suptitle("Point-ahead anisoplanatism of the screen plan: "
                 "D = 0.7 m, 10 arcsec, max_order None", fontsize=10)
    fig.tight_layout()
    import os
    path = os.path.join(fig_dir(), "discrete_vs_continuous.png")
    fig.savefig(path, dpi=140)
    plt.close(fig)
    return [path]


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--elevations", type=float, nargs="+",
                    default=list(ELEVATIONS_DEG),
                    help="the elevations of the sweep [deg]")
    ap.add_argument("--apertures", type=float, nargs="+",
                    default=list(APERTURES_M),
                    help="the receive aperture diameters [m]")
    ap.add_argument("--angles-arcsec", type=float, nargs="+",
                    default=list(ANGLES_ARCSEC),
                    help="the point-ahead angles [arcsec]; the geometry angle "
                         "of each elevation is always added")
    ap.add_argument("--verbose", action="store_true",
                    help="print every cell, not the production plan only")
    args = ap.parse_args()

    say, log_path = log_maker("anisoplanatism_screens")
    t0 = time.time()
    say("THE POINT-AHEAD ANISOPLANATISM OF THE SCREEN PLAN (no simulation)")
    say("case          : uplink, 1550 nm, 500 km, 700 mm ground terminal, "
        "HV5/7 site")
    say("reference     : Stone et al. (1994), DOI 10.1364/JOSAA.11.000347, "
        "Eqs. (29) and (36)")
    say("remove        : piston_tilt on every cell")
    say(f"elevations    : {args.elevations} deg")
    say(f"apertures     : {args.apertures} m")
    say(f"angles        : {args.angles_arcsec} arcsec, plus the geometry angle")
    say(f"max_order     : {MAX_ORDERS} (3 = 10 Noll modes, 5 = 21 Noll modes)")
    say()

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        rows, guard, checks, per_screen, order_sanity = sweep(say, args)

    say("THE CONTINUOUS REFERENCE against the library function "
        "anisoplanatic_phase_variance")
    for c in checks:
        for cell in c["cells"]:
            say(f"  {c['elevation_deg']:4.0f} deg, D = {cell['D_m']} m, "
                f"{cell['theta_arcsec']:.0f} arcsec, max_order "
                f"{cell['max_order']}: ratio {cell['ratio']:.6f}")
    say()
    v = verdict(say, rows)
    figs = plot(rows)
    for p in figs:
        say(f"wrote {p}")
    payload = {"study": "anisoplanatism_screens",
               "counts": list(COUNTS), "split_n": SPLIT_N,
               "max_orders": [None if m is None else int(m)
                              for m in MAX_ORDERS],
               "interpolation_guard": guard, "order_sanity": order_sanity,
               "continuous_checks": checks,
               "per_screen_30deg": per_screen,
               "rows": rows, "verdict": v,
               "runtime_s": time.time() - t0}
    json_path = write_results("anisoplanatism_screens", payload)
    say(f"wrote {json_path}")
    say(f"wrote {log_path}")
    say(f"runtime {time.time() - t0:.1f} s")


if __name__ == '__main__':
    main()
