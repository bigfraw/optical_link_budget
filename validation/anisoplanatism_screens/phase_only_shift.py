"""Do SHIFTED screen windows give the Stone point-ahead variance? Phase only.

THE QUESTION. `validation/anisoplanatism_screens/anisoplanatism_screens.py`
shows that the delta-layer Stone SUM of the production screen plan matches the
continuous Stone integral. That is a statement about the plan, not about the
screens. This script measures the SCREENS: it draws the production phase screens
on an oversized grid, it reads each screen two times (one window for the beacon,
one window for the uplink at the point-ahead offset), and it measures the
variance of the phase DIFFERENCE over a receive aperture. The measurement must
give the Stone number of the same plan.

NO PROPAGATION. The script sums the screen phases. It draws no field and it
calls no propagator. So it tests the screen generator and the offset recipe
only. It mirrors `validation/screen_stacking/screen_stacking.py`.

THE RECIPE.
  1. Take the PRODUCTION grid of the hero uplink: n pixels, the pixel dx.
  2. For screen j at the ground distance z_g take the offset
     s_j = round(theta z_g / dx) pixels. A screen is a grid, so the offset is a
     WHOLE pixel.
  3. Draw the stack on an oversized grid of n' = ceil32(n + s_max) pixels, the
     same pixel dx. The beacon window is scr[0:n, 0:n] and the uplink window is
     scr[0:n, s_j:s_j+n].
  4. Sum (uplink - beacon) over the screens. That is the phase difference of
     the two directions.
  5. Fit J = 36 Noll modes over each receive aperture, and measure
       (a) the PISTON-removed variance, which is the Stone quantity
           remove="piston", max_order=None; and
       (b) the band variance of the Noll modes 2 to J_ao, which is the Stone
           band max_order = max_radial_order(J_ao). The modes 2 to 10 fill the
           radial orders 1 to 3, and the modes 2 to 21 fill the orders 1 to 5,
           so the two bands match exactly.

THE TILT STAYS IN (owner decision, 2026-09-11). The terminal senses the downlink
beacon tilt and the steering mirror adds the point-ahead offset geometrically, so
the uplink pays the full tilt anisoplanatism. So every quantity above keeps the
tilt, and the reference is Stone with remove="piston".

THE REFERENCE is `discrete_stone` of the companion script, at the SAME
pixel-rounded offsets, plus the continuous Stone integral of the site profile.
The pass band is a ratio of 0.95 to 1.05 inside the bootstrap 2-sigma band.

THE OUTER SCALE. Run it at `--L0 inf` (the Kolmogorov limit) and at `--L0 25`
(the site operating value, backlog 2-P5). BOTH runs are pass tests, because the
Stone reference carries the SAME outer scale through the von Karman kernel of
`olb.turbulence.anisoplanatism._inner_integral` (Andrews and Phillips, 2nd ed.
(2005), DOI 10.1117/3.626196, Ch. 3, Eq. (20), printed p. 68). The 25 m run is
the numerical validation of that kernel.

Sources:
- J. Stone, P. H. Hu, S. P. Mills and S. Ma, J. Opt. Soc. Am. A 11(1), 347
  (1994), DOI 10.1364/JOSAA.11.000347. Eqs. (29), (36) and (A11).
- R. J. Noll, J. Opt. Soc. Am. 66(3), 207 (1976), DOI 10.1364/JOSA.66.000207.
  The Noll mode order of the fit.
- Schmidt, DOI 10.1117/3.866274, Ch. 9, Eqs. (9.78) to (9.81), printed
  pp. 166 to 169. The Fourier-series screen with the subharmonic levels.
- Johansson and Gavel, DOI 10.1117/12.177254. The subharmonic reach of a grid,
  which sets the level count of the reach guard.

Run it from the repository root:

    python -m validation.anisoplanatism_screens.phase_only_shift --n-draws 10
    python -m validation.anisoplanatism_screens.phase_only_shift
    python -m validation.anisoplanatism_screens.phase_only_shift --L0 25
"""

import argparse
import os
import time
import warnings

import numpy as np

from olb.turbulence.anisoplanatism import _REMOVE_NLO, max_radial_order
from olb.waveoptics.compensation import ApertureModes, circle
from olb.waveoptics.turbulence.screens import ScreenFactory
from validation.anisoplanatism_screens.anisoplanatism_screens import (
    InnerI, continuous_stone, discrete_stone, u0_of)
from validation.anisoplanatism_screens.common import (ARCSEC, LAM,
                                                      ground_distances,
                                                      hero_uplink, log_maker,
                                                      plan_set, production_plan,
                                                      pyplot, fig_dir,
                                                      write_results)

# The Noll modes of the fit. 36 modes fill the radial orders 0 to 7.
N_MODES = 36

# The adaptive-optics bands. Each J_ao maps to a Stone max_order through
# max_radial_order (Noll, DOI 10.1364/JOSA.66.000207, Table I).
AO_BANDS = (10, 21)

# The mode convention of record: the piston only comes out, so the TILT stays in
# (owner decision, 2026-09-11). See the module docstring.
REMOVE = "piston"

# The receive apertures of the measurement, in m.
APERTURES_M = (0.7, 0.4)

# The largest number of disjoint aperture windows that one draw reads.
MAX_WINDOWS = 25

# The seed of the draws. The screens are float32, as in a single-precision
# production run (olb.waveoptics.turbulence.run._screen_builder).
SEED = 20260911
SCREEN_DTYPE = np.float32

# The pass band of the ratio to the Stone reference.
PASS_LO, PASS_HI = 0.95, 1.05

N_BOOT = 400
BOOT_SEED = 2468


def ceil32(value):
    """Give the smallest multiple of 32 that is not less than value."""
    return int(32 * int(np.ceil(float(value) / 32.0)))


def shifts_px(plan, theta, dx):
    """Give the whole-pixel offset of each screen, from the ground distance."""
    return np.round(ground_distances(plan) * float(theta) / float(dx)
                    ).astype(int)


def aperture_reader(n, dx, aperture_m):
    """Give the tiled aperture windows and the modal fit of one diameter.

    THE WINDOW. `olb.waveoptics.compensation.zernike.circle` centres the mask on
    the pixel int(N/2), so the fit needs a SUB-GRID whose centre pixel holds the
    aperture. The function tiles the n x n grid with disjoint square windows of
    that size, and it builds ONE `ApertureModes` for the window.

    Args:
        n:          the pixel count of the grid.
        dx:         the pixel pitch, in m.
        aperture_m: the aperture diameter, in m.

    Returns:
        The triple (ApertureModes, a list of (row, col) window corners, the
        window size in pixels).
    """
    ap_px = int(round(float(aperture_m) / float(dx)))
    w = ap_px + 4 + (ap_px % 2)                      # an even window, 2 px edge
    mask = circle(w, ap_px)
    modes = ApertureModes(N_MODES, w, mask)
    starts = [i * w for i in range(n // w)]
    corners = [(r, c) for r in starts for c in starts]
    # A cap on the window count keeps the run time down. 25 disjoint apertures
    # per draw already cut the spread well below the bootstrap bar of the draws.
    return modes, corners[:MAX_WINDOWS], w


def difference_phase(stack, n, shifts):
    """Give the summed phase difference of the two directions on the n x n grid.

    Args:
        stack:  the screens on the oversized grid, in path order.
        n:      the pixel count of the propagation grid.
        shifts: the whole-pixel offset of each screen.

    Returns:
        An (n, n) float64 array of the phase difference, in radians.
    """
    out = np.zeros((n, n), dtype=np.float64)
    for scr, s in zip(stack, shifts):
        s = int(s)
        out += scr[0:n, s:s + n]
        out -= scr[0:n, 0:n]
    return out


def measure(fac, plan, shifts, n, readers, n_draws, seed):
    """Measure the aperture statistics of one (plan, angle) configuration.

    Args:
        fac:     the ScreenFactory of the oversized grid.
        plan:    the ScreenPlan.
        shifts:  the whole-pixel offset of each screen.
        n:       the pixel count of the propagation grid.
        readers: a dict aperture_m -> (ApertureModes, corners, window px).
        n_draws: the number of draws.
        seed:    the seed of the draws.

    Returns:
        A dict aperture_m -> {"piston": array, J_ao: array}, each array of
        the per-draw mean over the aperture windows, in rad^2.
    """
    rng = np.random.default_rng(seed)
    out = {D: {"piston": np.zeros(n_draws)}
           for D in readers}
    for D in readers:
        for j_ao in AO_BANDS:
            out[D][j_ao] = np.zeros(n_draws)
    for k in range(n_draws):
        stack = fac.make_stack(np.asarray(plan.r0_m, dtype=float), rng)
        diff = difference_phase(stack, n, shifts)
        for D, (modes, corners, w) in readers.items():
            acc = {"piston": [], **{j: [] for j in AO_BANDS}}
            for (r, c) in corners:
                window = diff[r:r + w, c:c + w]
                coeffs = modes.estimate(window)
                # (a) the PISTON-removed residual. The tilt stays in, so the
                #     fit that is REMOVED keeps the first Noll mode only.
                keep = np.zeros_like(coeffs)
                keep[:1] = coeffs[:1]
                acc["piston"].append(
                    modes.residual_variance(window, keep))
                # (b) the band variance of the Noll modes 2 to J_ao.
                for j_ao in AO_BANDS:
                    band = np.zeros_like(coeffs)
                    band[1:j_ao] = coeffs[1:j_ao]
                    acc[j_ao].append(
                        float(modes.reconstruct(band)[modes.mask].var()))
            for key, values in acc.items():
                out[D][key][k] = float(np.mean(values))
    return out


def bootstrap_two_sigma(x):
    """Give the bootstrap 2-sigma half-width of the mean of a sample."""
    x = np.asarray(x, dtype=float)
    rng = np.random.default_rng(BOOT_SEED)
    draws = np.array([x[rng.integers(0, x.size, x.size)].mean()
                      for _ in range(N_BOOT)])
    return float(2.0 * draws.std(ddof=1))


def reference(scn, el, plan, D, theta, dx, max_order, inner):
    """Give the discrete (pixel-rounded) and the continuous Stone variance.

    THE CONVENTION IS remove="piston": the tilt stays in. The outer scale sits
    in `inner`, which the caller built for this aperture.
    """
    disc = discrete_stone(plan, D, theta, LAM, remove=REMOVE,
                          max_order=max_order, dx=dx, inner=inner)["sigma2"]
    cont = continuous_stone(scn, el, D, theta, LAM, remove=REMOVE,
                            max_order=max_order, inner=inner)
    return disc, cont


def run_cell(say, scn, el, name, plan, arcsec, grid, L0, n_draws, readers,
             inners, sub_pixel=False):
    """Measure one (plan, angle) cell, and print its rows.

    Args:
        sub_pixel: True also measures the two lowest screens at the offsets
                   0 px and 1 px, to show the sensitivity of the rounding.

    Returns:
        A list of result dicts, one for each aperture and quantity.
    """
    n, dx = int(grid.n), float(grid.pixel_m)
    theta = float(arcsec) * ARCSEC
    s = shifts_px(plan, theta, dx)
    n_prime = ceil32(n + int(s.max()))
    t0 = time.time()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        fac = ScreenFactory(n_prime, dx, L0_m=L0, subharmonics=True,
                            dtype=SCREEN_DTYPE)
    reach = [str(w.message) for w in caught]
    stats = measure(fac, plan, s, n, readers, n_draws, SEED)

    variants = {}
    if sub_pixel:
        for label, forced in (("low0", 0), ("low1", 1)):
            s2 = s.copy()
            s2[-2:] = forced
            variants[label] = (s2, measure(fac, plan, s2, n, readers, n_draws,
                                           SEED))

    rows = []
    for D in readers:
        for key in ("piston", *AO_BANDS):
            max_order = (None if key == "piston"
                         else max_radial_order(int(key)))
            # The outer scale enters through u0 = 2 pi R / L0, so the table of
            # the Stone integral belongs to ONE aperture.
            inner = inners[(max_order, D)]
            disc, cont = reference(scn, el, plan, D, theta, dx, max_order,
                                   inner)
            x = stats[D][key]
            mean = float(x.mean())
            half = bootstrap_two_sigma(x)
            row = {"elevation_deg": el, "plan": name, "remove": REMOVE,
                   "L0_m": None if not np.isfinite(L0) else float(L0),
                   "n_screens": int(plan.z_m.size), "theta_arcsec": float(arcsec),
                   "aperture_m": float(D), "quantity": str(key),
                   "max_order": max_order, "grid_n": n, "grid_n_prime": n_prime,
                   "pixel_m": dx, "shift_px": s.tolist(),
                   "measured_rad2": mean,
                   "measured_two_sigma": half,
                   "stone_discrete_rad2": disc,
                   "stone_continuous_rad2": cont,
                   "ratio_discrete": mean / disc if disc > 0 else float("nan"),
                   "ratio_discrete_two_sigma": (half / disc if disc > 0
                                                else float("nan")),
                   "ratio_continuous": (mean / cont if cont > 0
                                        else float("nan")),
                   "n_draws": int(n_draws),
                   "n_apertures": len(readers[D][1]),
                   "reach_warnings": reach}
            for label, (_s2, st2) in variants.items():
                row[f"measured_{label}_rad2"] = float(st2[D][key].mean())
                row[f"delta_{label}_rad2"] = float(st2[D][key].mean() - mean)
            rows.append(row)
            say(f"  {name:<16s}{plan.z_m.size:3d}{arcsec:7.2f}{D:6.2f}"
                f"{str(key):>12s}{mean:11.4g}{disc:11.4g}"
                f"{row['ratio_discrete']:8.3f}"
                f" +-{row['ratio_discrete_two_sigma']:5.3f}"
                f"{row['ratio_continuous']:8.3f}")
    say(f"  ({time.time() - t0:.0f} s, grid {n_prime} px, top shift "
        f"{int(s.max())} px)")
    return rows


def plot(rows, tag):
    """Draw the ratio to the discrete Stone sum against the screen count."""
    plt = pyplot()
    if plt is None:
        return []
    quantities = ("piston", "10", "21")
    fig, axes = plt.subplots(1, len(quantities), figsize=(13.0, 4.4),
                             sharey=True)
    for ax, q in zip(axes, quantities):
        for arcsec in sorted({r["theta_arcsec"] for r in rows}):
            sel = [r for r in rows
                   if r["quantity"] == q and r["aperture_m"] == 0.7
                   and abs(r["theta_arcsec"] - arcsec) < 1e-9
                   and r["plan"].startswith("n")]
            sel.sort(key=lambda r: r["n_screens"])
            if sel:
                ax.errorbar([r["n_screens"] for r in sel],
                            [r["ratio_discrete"] for r in sel],
                            yerr=[r["ratio_discrete_two_sigma"] for r in sel],
                            fmt="-o", ms=4, capsize=3,
                            label=f"{arcsec:.0f} arcsec")
        ax.axhline(1.0, color="k", lw=0.8)
        ax.axhspan(PASS_LO, PASS_HI, color="0.85", zorder=0)
        ax.set_xlabel("screens in the plan")
        title = ("piston removed, tilt in" if q == "piston"
                 else f"Noll band 2 to {q}")
        ax.set_title(title, fontsize=9)
        ax.grid(alpha=0.3)
    axes[0].set_ylabel("measured / discrete Stone")
    axes[0].legend(fontsize=8)
    fig.suptitle("Shifted screen windows against Stone: D = 0.7 m", fontsize=10)
    fig.tight_layout()
    path = os.path.join(fig_dir(), f"phase_only_shift{tag}.png")
    fig.savefig(path, dpi=140)
    plt.close(fig)
    return [path]


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--elevations", type=float, nargs="+", default=[30.0],
                    help="the elevations of the study [deg]")
    ap.add_argument("--counts", type=int, nargs="+", default=[5, 9, 15],
                    help="the override screen counts")
    ap.add_argument("--angles-arcsec", type=float, nargs="+",
                    default=[2.0, 5.0, 10.0],
                    help="the point-ahead angles [arcsec]")
    ap.add_argument("--n-draws", type=int, default=100,
                    help="the stacks drawn for each cell")
    ap.add_argument("--L0", type=float, default=np.inf,
                    help="the outer scale of the screens [m]; the default is "
                         "infinite (the Kolmogorov limit of the Stone law)")
    ap.add_argument("--no-ground-split", action="store_true",
                    help="skip the ground-split plan")
    args = ap.parse_args()
    L0 = float(args.L0)
    tag = "" if not np.isfinite(L0) else f"_L0{L0:g}"

    say, log_path = log_maker(f"phase_only_shift{tag}")
    t0 = time.time()
    say("THE SHIFTED-WINDOW POINT-AHEAD TEST (phase only, no propagation)")
    say("case          : uplink, 1550 nm, 500 km, 700 mm ground terminal, "
        "HV5/7 site")
    say("outer scale   : "
        + ("L0 is infinite (the screens AND the Stone reference)"
           if not np.isfinite(L0) else
           f"L0 = {L0:g} m (the screens AND the von Karman Stone reference)"))
    say(f"draws         : {args.n_draws} per cell; seed {SEED}; float32 screens")
    say(f"fit           : {N_MODES} Noll modes; remove={REMOVE} (the tilt "
        f"stays in); bands 2 to "
        f"{' and 2 to '.join(str(j) for j in AO_BANDS)}")
    say("reference     : the delta-layer Stone sum at the SAME pixel-rounded "
        "offsets, and the continuous Stone integral")
    say()

    # One Stone table for each (max_order, aperture) pair: the outer scale
    # enters the kernel through u0 = 2 pi R / L0, and R is the aperture radius.
    n_lo = _REMOVE_NLO[REMOVE]
    inners = {}
    for D in APERTURES_M:
        u0 = u0_of(D, L0)
        inners[(None, D)] = InnerI(n_lo, None, u0)
        for j_ao in AO_BANDS:
            m = max_radial_order(int(j_ao))
            inners[(m, D)] = InnerI(n_lo, m, u0)

    rows = []
    for el in args.elevations:
        scn, geom = hero_uplink(el)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            grid, plans = plan_set(scn, geom, args.counts, split_n=4)
        if args.no_ground_split:
            plans = {k: v for k, v in plans.items()
                     if not k.startswith("ground_split")}
        n, dx = int(grid.n), float(grid.pixel_m)
        readers = {D: aperture_reader(n, dx, D) for D in APERTURES_M}
        say(f"ELEVATION {el:.0f} deg: grid {n} px, pixel {dx * 1e3:.3f} mm; "
            + ", ".join(f"{D} m x {len(readers[D][1])} windows"
                        for D in APERTURES_M))
        hdr = (f"  {'plan':<16s}{'n':>3s}{'th[as]':>7s}{'D[m]':>6s}"
               f"{'quantity':>12s}{'meas':>11s}{'stone':>11s}"
               f"{'ratio':>8s}{'+-2sig':>8s}{'r_cont':>8s}")
        say(hdr)
        say("  " + "-" * (len(hdr) - 2))
        for name, plan in plans.items():
            if name == "production":
                continue            # the override plan n9 is the same plan
            for arcsec in args.angles_arcsec:
                rows += run_cell(say, scn, el, name, plan, arcsec, grid, L0,
                                 args.n_draws, readers, inners)
        say()

        # THE SUB-PIXEL ARM. The production plan at 30 deg, 10 arcsec. The two
        # lowest screens have an offset of a few pixels only, so the rounding
        # can change their share. The arm forces those two offsets to 0 px and
        # to 1 px.
        if abs(el - 30.0) < 1e-9:
            _, prod, _ = production_plan(scn, geom)
            say("  THE SUB-PIXEL ARM: the production plan at 10 arcsec, with "
                "the two LOWEST screens forced to 0 px and to 1 px")
            sub = run_cell(say, scn, el, "production_subpixel", prod, 10.0,
                           grid, L0, args.n_draws, readers, inners,
                           sub_pixel=True)
            for r in sub:
                if r["aperture_m"] == 0.7:
                    say(f"    D = {r['aperture_m']} m, {r['quantity']}: "
                        f"rounded {r['measured_rad2']:.4g}, 0 px "
                        f"{r['measured_low0_rad2']:.4g} "
                        f"({r['delta_low0_rad2']:+.3g}), 1 px "
                        f"{r['measured_low1_rad2']:.4g} "
                        f"({r['delta_low1_rad2']:+.3g}) rad^2")
            rows += sub
            say()

    # The verdict.
    say("THE VERDICT")
    main_rows = [r for r in rows if r["plan"].startswith("n")]
    inside = [r for r in main_rows
              if PASS_LO <= r["ratio_discrete"] <= PASS_HI
              or abs(r["ratio_discrete"] - 1.0) <= r["ratio_discrete_two_sigma"]]
    worst = max(main_rows, key=lambda r: abs(r["ratio_discrete"] - 1.0))
    say(f"  {len(inside)} of {len(main_rows)} override-plan cells sit in the "
        f"{PASS_LO:.2f} to {PASS_HI:.2f} band, or inside their own 2-sigma bar.")
    say(f"  The worst cell: {worst['plan']}, {worst['theta_arcsec']:.2f} "
        f"arcsec, D = {worst['aperture_m']} m, {worst['quantity']}, ratio "
        f"{worst['ratio_discrete']:.3f} +- "
        f"{worst['ratio_discrete_two_sigma']:.3f}.")
    by_q = {}
    for q in ("piston", *[str(j) for j in AO_BANDS]):
        sel = [r["ratio_discrete"] for r in main_rows if r["quantity"] == q]
        if sel:
            by_q[q] = [min(sel), max(sel)]
            say(f"  {q:>11s}: the ratio holds between {min(sel):.3f} and "
                f"{max(sel):.3f}.")
    figs = plot(rows, tag)
    for p in figs:
        say(f"wrote {p}")
    payload = {"study": "phase_only_shift", "L0_m": None if not
               np.isfinite(L0) else L0, "n_draws": int(args.n_draws),
               "seed": SEED, "n_modes": N_MODES, "ao_bands": list(AO_BANDS),
               "apertures_m": list(APERTURES_M), "rows": rows,
               "verdict": {"n_inside": len(inside), "n_cells": len(main_rows),
                           "worst": worst, "ratio_range_by_quantity": by_q},
               "runtime_s": time.time() - t0}
    json_path = write_results(f"phase_only_shift{tag}", payload)
    say(f"wrote {json_path}")
    say(f"wrote {log_path}")
    say(f"runtime {(time.time() - t0) / 60.0:.1f} min")


if __name__ == '__main__':
    main()
