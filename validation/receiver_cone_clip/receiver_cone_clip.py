"""Test a CONVERGING absorbing boundary that follows the receiver cone.

THE CLAIM UNDER TEST (owner, 2026-09-06, backlog 2-P3, route (b)). The
turbulent sizer holds the WHOLE beam plus the scatter cone all the way to the
receiver. But only the light inside the BACK-PROJECTED cone of the receive
aperture can land on that aperture. If the claim holds, an absorbing boundary
that converges onto that cone leaves the field INSIDE the aperture unchanged,
and it removes most of the grid extent. A smaller extent is a finer pixel at
the same pixel count, which is the whole point: at 10 km the sizer asks for a
9.2 m side and clamps the pixel count at 2048.

THE METHOD. The grid does NOT change. One fixed grid, one fixed screen plan,
and, for each trial, ONE set of phase screens. Only the MASK differs between
the two arms. So every difference inside the receive aperture is the doing of
the converging clip, and nothing else.

  - REFERENCE arm: the production edge mask only. The script rebuilds the
    production loop by hand and it ASSERTS `numpy.array_equal` against the
    public `propagate_turbulent_field`, so the hand loop is known-correct.
  - CLIPPED arm: the SAME loop and the SAME screens, with the edge mask
    multiplied by a converging super-Gaussian at every plane.

THE CONE RULE. At the plane z of a path of the length L, with a receive
diameter D, the script offers TWO forms. The option --cone-rule selects one.

  MAPPED (the first form, tested on 2026-09-06; --cone-rule mapped):

    R(z) = c * [ D/2 + theta_s(z) * (L - z) ] * (z / L)  +  delta(z)

    theta_s(z) = lambda / r0_rest(z),  r0_rest = (0.423 k^2 Cn2 (L - z))^(-3/5)

  UNMAPPED (the default; --cone-rule unmapped):

    R(z) = c * [ (D/2) * (z / L)  +  theta_s * (L - z) ]  +  delta(z)

    theta_s = lambda / r0_full,  r0_full = (0.423 k^2 Cn2 L)^(-3/5)

  Both forms then take the option guard:

    R(z) -> max( R(z), 3 * w(z) )        with --guard on
    delta(z) = 3 * lambda * (L - z) / D  in both forms

  - The bracket term D/2 is the back-projected receive aperture. The factor
    z/L is the FAR-FIELD ray mapping: the Rayleigh range of the 5 mm launch
    waist is about 51 m, so a 5 km path is far field, a ray at the transverse
    position x at the plane z leaves the axis at the angle x/z and it arrives
    at x * L / z, and to land inside the aperture it needs x <= (D/2) * (z/L).
    Source of the far-field Gaussian beam: Andrews and Phillips, Laser Beam
    Propagation through Random Media, 2nd ed. (2005), DOI 10.1117/3.626196,
    Ch. 4, Eqs. (7) and (8), printed p. 87.
  - theta_s * (L - z) is the turbulent SCATTER allowance. The MAPPED form
    multiplies it by z/L too. The UNMAPPED form does NOT, because scattered
    light does not obey the ballistic ray map: a screen sends it away at a
    fresh angle, and it reaches the receiver over theta_s * (L - z) wherever
    z sits. The UNMAPPED form also reads theta_s from the WHOLE path, not
    from the remaining path, because the light at the plane z can already
    carry the scatter that the path behind it put on. r0 is the plane-wave
    Fried parameter: Andrews and Phillips, DOI 10.1117/3.626196, Ch. 6,
    Eq. (64) (r_0 = 2.1 rho_0, with rho_0 of a plane wave).
  - delta(z) is a fixed RAMP ALLOWANCE. A mask of the radius R diffracts with
    the angular width of about lambda / R; over the remaining path (L - z)
    that spreads by lambda (L - z) / R. The allowance uses R = D/2 ... in the
    form 3 lambda (L - z) / D, three times that spread at the aperture scale.
    So the diffraction cone of the mask itself stays inside the aperture
    cone, and the mask does not build its own ring on the receiver.
  - The guard 3 * w(z) is the near-transmitter guard: never clip inside three
    vacuum beam radii, where the whole launch waist matters. w(z) is the
    vacuum Gaussian radius, `olb.beam.gaussz`. It is OFF by default. On a
    5 km path a 5 mm waist grows to w(L) = 0.49 m, so 3 w(L) = 1.48 m, which
    is LARGER than the 1.43 m half-side of the sizer grid: the guard then
    covers the whole grid and the clip is vacuous. --guard on keeps it, and
    `--cone-rule mapped --guard on` reproduces the first run of 2026-09-06.

THE MASK SHAPE. M(r; z) = exp(-(r / R(z))^16). It is the same super-Gaussian
family as the production edge mask, the absorbing operator A of Schmidt,
Numerical Simulation of Optical Wave Propagation with Examples in MATLAB
(2010), DOI 10.1117/3.866274, Ch. 8, Eq. (8.18), printed p. 139, whose shape
comes from Ch. 8, Eq. (8.1), printed p. 134. The power 16 is the book's own
run value (Listing 8.1, printed p. 142). NOTE that this mask is 1/e at r = R,
not 1.0: R is the 1/e amplitude radius, not a hard edge.

THE CELL. `validation.terrestrial_campaigns.run_campaigns.build_scenario`
with 5 km, Cn2 = 1e-14, a collimated 5 mm waist at 1550 nm, and a 10 cm
receive aperture with an optimal-focus single-mode fibre. The preset is
`rapid`, the outer scale is 25 m (the owner fixed it on 2026-09-05, backlog
2-P5), the precision is single, and the seed is 20260906. A VACUUM arm (no
screens) runs the same masks.

THE PASS RULE. For one factor c, EVERY trial must give, inside the 10 cm
aperture:
      the relative field RMS difference   < 1e-3
      |the collected power difference|    < 0.01 dB
      |the fibre efficiency difference|   < 1e-3
1e-3 sits well ABOVE the single-precision floor: `validation/precision`
measured a field RMS of 1.3e-6 between a complex64 and a complex128 run of the
same seed, so the test cannot trip on the arithmetic. It sits well BELOW any
budget-visible effect: 0.01 dB is under the 0.11 dB screen-count spread of
docs WP7, and a fibre efficiency step of 1e-3 on an eta of about 0.01 to 0.05
is a few hundredths of a dB.

The script reports the smallest passing factor TWO ways: under the FIELD rule
(the 1e-3 field RMS alone) and under the BUDGET-VISIBLE rule alone (the 0.01 dB
power step and the 1e-3 efficiency step, which is all that a Term reads).

THE OUTPUT NAMES carry the rule, the guard and the cell, so a run does not
write on the record of another run:

    receiver_cone_clip_<rule>_<guard>_<cell>.log
    receiver_cone_clip_<rule>_<guard>_<cell>.json
    figures/receiver_cone_clip_<rule>_<guard>_<cell>.png

The first run of 2026-09-06 kept the older, untagged names
`receiver_cone_clip.log` and `receiver_cone_clip_results.json`.

Run it from the repository root:

    python -m validation.receiver_cone_clip.receiver_cone_clip
    python -m validation.receiver_cone_clip.receiver_cone_clip --guard on
    python -m validation.receiver_cone_clip.receiver_cone_clip         --cone-rule mapped --guard on
    python -m validation.receiver_cone_clip.receiver_cone_clip         --path-km 10 --cn2 1e-14 --factors 1.5 2 3 5
"""

import argparse
import json
import os
import sys
import time

import numpy as np

from olb.beam import gaussz, virtual_waist
from olb.waveoptics.field import Begin, Field, Power, field_dtype
from olb.waveoptics.propagators import Forvard
from olb.waveoptics.turbulence.run import (_detector_eta, _screen_builder,
                                           _screen_seed, _start_field)
from olb.waveoptics.turbulence.sampling import PRESETS, turbulent_grid
from olb.waveoptics.turbulence.screens import Screen
from olb.waveoptics.turbulence.splitstep import (_apply_mask, _substeps,
                                                 split_step,
                                                 super_gaussian_boundary)
from validation.terrestrial_campaigns.run_campaigns import (LAM, WAIST_M,
                                                            build_scenario)

HERE = os.path.dirname(os.path.abspath(__file__))
FIGDIR = os.path.join(HERE, "figures")

L0_M = 25.0
PRECISION = "single"
SEED = 20260906

# The pass thresholds. See the module docstring.
RMS_MAX = 1e-3
DP_DB_MAX = 0.01
DETA_MAX = 1e-3

# The exponent of the converging super-Gaussian, and the clip of its argument.
# 2^16 = 65536, and exp(-65536) is already zero in float64, so the clip only
# stops an overflow warning.
CONE_POWER = 16
CONE_ARG_MAX = 2.0

# The ramp allowance, in units of lambda (L - z) / D. See the docstring.
RAMP_FACTOR = 3.0

# The near-transmitter guard, in vacuum beam radii.
BEAM_FLOOR_RADII = 3.0


class Tee:
    """Write to the terminal and to a log file at the same time."""

    def __init__(self, path):
        self.stream = open(path, "w", encoding="utf-8")

    def write(self, text):
        sys.__stdout__.write(text)
        self.stream.write(text)

    def flush(self):
        sys.__stdout__.flush()
        self.stream.flush()

    def close(self):
        self.stream.close()


# ---------------------------------------------------------------------------
# The cone
# ---------------------------------------------------------------------------

def rest_fried_m(rest_m, cn2, lam):
    """Give the plane-wave Fried parameter of the path that is still to come.

    r_0 = (0.423 k^2 Cn2 z)^(-3/5). Source: Andrews and Phillips, 2nd ed.
    (2005), DOI 10.1117/3.626196, Ch. 6, Eq. (64) (the plane-wave coherence
    radius rho_0, with r_0 = 2.1 rho_0 giving the same 0.423 constant).

    Args:
        rest_m: the remaining path from the plane to the receiver, in m.
        cn2:    the constant Cn2 of the path, in m^-2/3.
        lam:    the wavelength, in m.

    Returns:
        The Fried parameter in m. It is infinite for a zero remaining path.
    """
    if rest_m <= 0.0 or cn2 <= 0.0:
        return float("inf")
    k = 2.0 * np.pi / lam
    return float((0.423 * k * k * cn2 * rest_m) ** (-3.0 / 5.0))


def cone_radius_m(z_m, path_m, diameter_m, cn2, lam, factor, w_v, offset_m,
                  floor, rule="unmapped"):
    """Give the radius R(z) of the converging mask at one plane.

    See the module docstring for the two rules and for every source.

    Args:
        z_m:        the distance from the launch plane, in m.
        path_m:     the whole path length L, in m.
        diameter_m: the receive aperture diameter D, in m.
        cn2:        the Cn2 of the path, in m^-2/3.
        lam:        the wavelength, in m.
        factor:     the cone factor c.
        w_v:        the virtual waist of the launch, in m.
        offset_m:   the distance of the launch plane behind the virtual
                    waist, in m.
        floor:      True applies the R_beam_min = 3 w(z) guard.
        rule:       "mapped" puts the ray map z/L on the whole bracket and it
                    reads theta_s from the REMAINING path; "unmapped" puts
                    the ray map on the aperture term only and it reads
                    theta_s from the WHOLE path.

    Returns:
        The radius, in m.
    """
    rest = max(path_m - z_m, 0.0)
    ramp = RAMP_FACTOR * lam * rest / diameter_m
    if rule == "mapped":
        r0 = rest_fried_m(rest, cn2, lam)
        scatter = 0.0 if not np.isfinite(r0) else (lam / r0) * rest
        cone = factor * (diameter_m / 2.0 + scatter) * (z_m / path_m) + ramp
    elif rule == "unmapped":
        r0 = rest_fried_m(path_m, cn2, lam)
        scatter = 0.0 if not np.isfinite(r0) else (lam / r0) * rest
        cone = factor * (diameter_m / 2.0 * (z_m / path_m) + scatter) + ramp
    else:
        raise ValueError(f"unknown cone rule {rule!r}")
    if floor:
        w_z = float(gaussz(w_v, offset_m + z_m, lam))
        cone = max(cone, BEAM_FLOOR_RADII * w_z)
    return float(cone)


def cone_masks(grid, path_m, diameter_m, cn2, lam, factor, w_v, offset_m,
               floor, rdtype, rule="unmapped"):
    """Give a cache of converging masks, keyed by the plane.

    The planes of a run are fixed, so one mask serves every trial and every
    application at that plane.

    Returns:
        A callable mask(z_m) -> an n x n array, or None when the mask is
        1.0 over the whole grid.
    """
    n = int(grid.n)
    x = (np.arange(n) - n // 2) * grid.pixel_m
    rr = np.hypot(x[:, None], x[None, :])
    r_corner = float(rr.max())
    cache = {}

    def mask(z_m):
        key = round(float(z_m), 6)
        if key not in cache:
            R = cone_radius_m(z_m, path_m, diameter_m, cn2, lam, factor,
                              w_v, offset_m, floor, rule)
            if R >= r_corner * CONE_ARG_MAX:
                cache[key] = None       # 1.0 everywhere on the grid
            else:
                t = np.clip(rr / R, 0.0, CONE_ARG_MAX)
                cache[key] = np.exp(-(t ** CONE_POWER)).astype(rdtype)
        return cache[key]

    return mask


# ---------------------------------------------------------------------------
# The propagation loop
# ---------------------------------------------------------------------------

def walk(Fin, z_screens_m, screens, z_total_m, edge_mask, cone_fn=None):
    """Reproduce olb.waveoptics.turbulence.splitstep.split_step by hand.

    The loop is the split_step loop, line for line: one hop to each screen
    plane in sub-steps that obey N dx^2 / lambda (Schmidt, DOI
    10.1117/3.866274, Ch. 8, Eq. (8.24), printed p. 144), the edge mask after
    each sub-step, the screen, the edge mask again, and a last hop.

    With cone_fn=None the loop gives the SAME array as split_step, bit for
    bit. With a cone_fn it applies the converging mask at exactly the planes
    where split_step applies the edge mask, so the two arms differ in the mask
    and in nothing else.

    Args:
        Fin:         the input Field.
        z_screens_m: the screen distances from the input plane, in m.
        screens:     one phase array for each distance.
        z_total_m:   the path length, in m.
        edge_mask:   the production super-Gaussian edge mask.
        cone_fn:     a callable mask(z) -> array or None, or None for the
                     reference arm.

    Returns:
        The pair (Field at z_total_m, the power that the CONE removed).
    """
    max_step_m = Fin.N * Fin.dx ** 2 / Fin.lam
    rdtype = np.float32 if Fin.field.dtype == np.complex64 else np.float64
    edge = np.asarray(edge_mask, dtype=rdtype)
    state = {"removed": 0.0}

    def cone(F, z_m):
        if cone_fn is None:
            return F
        m = cone_fn(z_m)
        if m is None:
            return F
        p0 = Power(F)
        F = _apply_mask(F, m)
        state["removed"] += float(p0 - Power(F))
        return F

    def hop(F, gap_m, here_m):
        for dz in _substeps(gap_m, max_step_m):
            F = Forvard(F, dz)
            here_m += dz
            F = _apply_mask(F, edge)
            F = cone(F, here_m)
        return F, here_m

    z = np.asarray(z_screens_m, dtype=float).ravel()
    F = Field.copy(Fin)
    here = 0.0
    for zi, scr in zip(z, screens):
        F, _ = hop(F, zi - here, here)
        here = float(zi)
        F = Screen(F, scr)
        F = _apply_mask(F, edge)
        F = cone(F, here)
    F, _ = hop(F, z_total_m - here, here)
    return F, state["removed"]


# ---------------------------------------------------------------------------
# The metrics
# ---------------------------------------------------------------------------

def compare(F_ref, F_clip, ap_mask, centre, detector, diameter_m,
            obscuration, lam, clip_fn):
    """Compare a clipped arm with the reference arm inside the aperture.

    Args:
        F_ref:       the reference Field at the receive plane.
        F_clip:      the clipped Field at the receive plane.
        ap_mask:     the boolean aperture mask of the receive aperture.
        centre:      the (row, column) index of the axis pixel.
        detector:    the receive detector.
        diameter_m:  the receive aperture diameter, in m.
        obscuration: the central obscuration ratio.
        lam:         the wavelength, in m.
        clip_fn:     the aperture clip, olb.waveoptics.run._clip.

    Returns:
        A dict of the four metrics.
    """
    a = F_ref.field[ap_mask].astype(np.complex128)
    b = F_clip.field[ap_mask].astype(np.complex128)
    denom = float(np.sqrt((np.abs(a) ** 2).sum()))
    rms = float(np.sqrt((np.abs(b - a) ** 2).sum()) / denom)

    c_ref = clip_fn(F_ref, diameter_m, obscuration)
    c_clip = clip_fn(F_clip, diameter_m, obscuration)
    p_ref, p_clip = float(Power(c_ref)), float(Power(c_clip))
    dp_db = float(10.0 * np.log10(p_clip / p_ref))

    e_ref = _detector_eta(detector, c_ref, diameter_m, lam)
    e_clip = _detector_eta(detector, c_clip, diameter_m, lam)
    d_eta = float(e_clip - e_ref)

    i_ref = float(abs(F_ref.field[centre]) ** 2)
    i_clip = float(abs(F_clip.field[centre]) ** 2)
    ratio = float(i_clip / i_ref) if i_ref > 0 else float("nan")
    return {"field_rms": rms, "dp_db": dp_db, "d_eta": d_eta,
            "eta_ref": float(e_ref), "eta_clip": float(e_clip),
            "centre_ratio": ratio}


def passes(rows):
    """True when every trial row of one factor meets the pass rule."""
    return all(r["field_rms"] < RMS_MAX and abs(r["dp_db"]) < DP_DB_MAX
               and abs(r["d_eta"]) < DETA_MAX for r in rows)


def passes_field(rows):
    """True when every trial row meets the FIELD rule alone.

    The field rule is the 1e-3 relative field RMS in the aperture. It is the
    strict test: the clipped field must BE the reference field.
    """
    return all(r["field_rms"] < RMS_MAX for r in rows)


def passes_budget(rows):
    """True when every trial row meets the BUDGET-VISIBLE rule alone.

    A Term reads the collected power and the coupling efficiency, and it does
    not read the field. So this is the rule that decides whether a budget can
    see the clip.
    """
    return all(abs(r["dp_db"]) < DP_DB_MAX and abs(r["d_eta"]) < DETA_MAX
               for r in rows)


# ---------------------------------------------------------------------------
# The run
# ---------------------------------------------------------------------------

def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Test a converging receiver-cone absorbing boundary "
                    "against the production edge mask, on a fixed grid.")
    ap.add_argument("--trials", type=int, default=6,
                    help="the number of turbulent trials (default 6).")
    ap.add_argument("--factors", type=float, nargs="+",
                    default=[1.0, 1.5, 2.0, 3.0, 5.0],
                    help="the cone factors c to sweep.")
    ap.add_argument("--cone-rule", default="unmapped",
                    choices=["mapped", "unmapped"],
                    help="'mapped' puts the far-field ray map z/L on the "
                         "whole bracket, the first form of 2026-09-06; "
                         "'unmapped' (the default) puts it on the aperture "
                         "term only and it gives the scatter term the full "
                         "path allowance.")
    ap.add_argument("--guard", default="off", choices=["on", "off"],
                    help="'on' keeps the R >= 3 w(z) near-transmitter guard; "
                         "'off' (the default) removes it, because on a 5 km "
                         "path the guard covers the whole grid.")
    ap.add_argument("--preset", default="rapid", help="the quality preset.")
    ap.add_argument("--seed", type=int, default=SEED, help="the run seed.")
    ap.add_argument("--path-km", type=float, default=5.0,
                    help="the horizontal path length, in km.")
    ap.add_argument("--cn2", type=float, default=1e-14,
                    help="the Cn2 of the path, in m^-2/3.")
    args = ap.parse_args(argv)

    # The guard is one value now. The body still loops over a list, because
    # the first run swept both guards in one pass.
    args.floors = ["beam" if args.guard == "on" else "none"]
    args.tag = (f"{args.cone_rule}_{args.guard}_"
                f"{args.path_km:g}km").replace(".", "p")

    os.makedirs(FIGDIR, exist_ok=True)
    log = Tee(os.path.join(HERE, f"receiver_cone_clip_{args.tag}.log"))
    old_stdout = sys.stdout
    sys.stdout = log
    try:
        return _run(args)
    finally:
        sys.stdout = old_stdout
        log.close()


def _run(args):
    from olb.waveoptics.run import _clip
    from olb.waveoptics.turbulence.run import propagate_turbulent_field

    t_start = time.time()
    path_m = args.path_km * 1e3
    scn, geom = build_scenario(path_m, args.cn2, "collimated")
    rx = scn.rx_terminal
    lam = scn.tx_terminal.wavelength_m
    preset = PRESETS[args.preset]

    grid, plan, report = turbulent_grid(scn, geom, preset=preset, L0_m=L0_M)
    cdtype = field_dtype(PRECISION)
    rdtype = np.float32 if cdtype == np.complex64 else np.float64
    edge = super_gaussian_boundary(grid.n, preset.boundary_width_frac)
    build_screen = _screen_builder("olb", grid, L0_M, True, dtype=cdtype)
    F_in = _start_field(scn, grid, lam, is_space=False, dtype=cdtype)
    w_v, offset_m = virtual_waist(WAIST_M, scn.tx_terminal.transmitter
                                  .divergence_rad, lam)

    # The aperture pixels, in the exact convention of CircAperture.
    Y, X = F_in.mgrid_cartesian
    ap_mask = (X ** 2 + Y ** 2) <= (rx.aperture_m / 2.0) ** 2
    centre = np.unravel_index(int(np.argmin(X ** 2 + Y ** 2)), X.shape)

    print("=" * 78)
    print("THE RECEIVER-CONE CLIP: a converging absorbing boundary")
    print("=" * 78)
    print(f"path                 {path_m / 1e3:8.2f} km")
    print(f"Cn2                  {args.cn2:8.2e} m^-2/3")
    print(f"wavelength           {lam * 1e9:8.1f} nm")
    print(f"launch waist         {WAIST_M * 1e3:8.2f} mm")
    print(f"receive aperture     {rx.aperture_m * 1e2:8.2f} cm")
    print(f"preset               {args.preset:>8s}")
    print(f"cone rule            {args.cone_rule:>8s}")
    print(f"guard (3 w(z))       {args.guard:>8s}")
    print(f"grid                 {grid.n:8d} px, side {grid.size_m:.4f} m, "
          f"pixel {grid.pixel_m * 1e3:.3f} mm")
    print(f"screens              {plan.z_m.size:8d}")
    print(f"r0 of the whole path {plan.r0_total_m * 1e2:8.2f} cm")
    print(f"outer scale L0       {L0_M:8.1f} m")
    print(f"precision            {PRECISION:>8s}, seed {args.seed}")
    print(f"boundary band        {preset.boundary_width_frac:8.3f} of the "
          "half-side")
    print(f"aperture pixels      {int(ap_mask.sum()):8d}")
    print("")

    # ---- the reference loop is the production loop ----
    scr0 = [build_screen(_screen_seed(args.seed, 0, j), plan.r0_m[j])
            for j in range(int(plan.z_m.size))]
    F_hand, _ = walk(F_in, plan.z_m, scr0, plan.z_total_m, edge)
    F_pub, _, _ = propagate_turbulent_field(
        scn, geom, seed=args.seed, trial=0, preset=args.preset, grid=grid,
        plan=plan, L0_m=L0_M, precision=PRECISION)
    assert np.array_equal(F_hand.field, F_pub.field), \
        "the hand loop does not reproduce propagate_turbulent_field"
    print("CHECK 1  the hand loop equals propagate_turbulent_field, bit for "
          "bit.")

    # A cone_fn that gives None everywhere must also be a no-op.
    F_none, removed_none = walk(F_in, plan.z_m, scr0, plan.z_total_m, edge,
                                cone_fn=lambda z: None)
    assert np.array_equal(F_none.field, F_hand.field)
    assert removed_none == 0.0
    print("CHECK 2  the clipped loop with an empty mask is the same array.")
    print("")

    launched = float(Power(F_in))
    results = {"cell": {"path_m": path_m, "cn2": args.cn2, "lam": lam,
                        "cone_rule": args.cone_rule, "guard": args.guard,
                        "preset": args.preset, "seed": args.seed,
                        "L0_m": L0_M, "precision": PRECISION,
                        "grid_n": int(grid.n), "grid_size_m": grid.size_m,
                        "pixel_m": grid.pixel_m,
                        "n_screens": int(plan.z_m.size),
                        "r0_total_m": plan.r0_total_m,
                        "aperture_m": rx.aperture_m,
                        "trials": int(args.trials),
                        "boundary_width_frac": preset.boundary_width_frac,
                        "thresholds": {"field_rms": RMS_MAX,
                                       "dp_db": DP_DB_MAX,
                                       "d_eta": DETA_MAX}},
               "turbulence": {}, "vacuum": {}, "cone_profile": {}}

    # ---- the turbulent arms ----
    ref_fields = {}
    for k in range(args.trials):
        scr = [build_screen(_screen_seed(args.seed, k, j), plan.r0_m[j])
               for j in range(int(plan.z_m.size))]
        ref_fields[k] = (scr, walk(F_in, plan.z_m, scr, plan.z_total_m,
                                   edge)[0])

    # ---- the vacuum reference ----
    F_vac_ref = split_step(F_in, [], [], plan.z_total_m, boundary=edge)

    for floor_name in args.floors:
        floor = floor_name == "beam"
        results["turbulence"][floor_name] = {}
        results["vacuum"][floor_name] = {}
        for c in args.factors:
            fn = cone_masks(grid, path_m, rx.aperture_m, args.cn2, lam, c,
                            w_v, offset_m, floor, rdtype, args.cone_rule)
            rows, removed_fracs = [], []
            for k in range(args.trials):
                scr, F_ref = ref_fields[k]
                F_clip, removed = walk(F_in, plan.z_m, scr, plan.z_total_m,
                                       edge, cone_fn=fn)
                row = compare(F_ref, F_clip, ap_mask, centre, rx.detector,
                              rx.aperture_m, rx.obscuration_ratio, lam, _clip)
                row["trial"] = k
                row["removed_frac"] = removed / launched
                rows.append(row)
                removed_fracs.append(removed / launched)
            results["turbulence"][floor_name][f"{c:g}"] = {
                "trials": rows, "pass": passes(rows),
                "removed_frac_mean": float(np.mean(removed_fracs)),
                "removed_frac_max": float(np.max(removed_fracs))}

            F_vc, removed_v = walk(F_in, [], [], plan.z_total_m, edge,
                                   cone_fn=fn)
            vrow = compare(F_vac_ref, F_vc, ap_mask, centre, rx.detector,
                           rx.aperture_m, rx.obscuration_ratio, lam, _clip)
            vrow["trial"] = 0
            vrow["removed_frac"] = removed_v / launched
            results["vacuum"][floor_name][f"{c:g}"] = {
                "trials": [vrow], "pass": passes([vrow]),
                "removed_frac_mean": float(removed_v / launched),
                "removed_frac_max": float(removed_v / launched)}

    # ---- the print ----
    for floor_name in args.floors:
        guard = ("ON  (R >= 3 w(z), the rule as written)" if floor_name ==
                 "beam" else "OFF (the cone alone)")
        for arm in ("turbulence", "vacuum"):
            n_arm = args.trials if arm == "turbulence" else 1
            print("-" * 78)
            print(f"{arm.upper()}   beam-radius floor {guard}   "
                  f"({n_arm} trial(s))")
            print("-" * 78)
            print(f"{'c':>5s} {'removed':>10s} {'max RMS':>11s} "
                  f"{'max |dP| dB':>12s} {'max |d eta|':>12s} "
                  f"{'centre ratio':>13s}  verdict")
            for c in args.factors:
                e = results[arm][floor_name][f"{c:g}"]
                rows = e["trials"]
                print(f"{c:5.2f} {e['removed_frac_max']:10.3e} "
                      f"{max(r['field_rms'] for r in rows):11.3e} "
                      f"{max(abs(r['dp_db']) for r in rows):12.3e} "
                      f"{max(abs(r['d_eta']) for r in rows):12.3e} "
                      f"{np.mean([r['centre_ratio'] for r in rows]):13.6f}"
                      f"  {'PASS' if e['pass'] else 'FAIL'}")
            print("")

    # ---- the verdict and the saving ----
    z_fine = np.linspace(0.0, path_m, 501)
    w_fine = np.array([float(gaussz(w_v, offset_m + z, lam)) for z in z_fine])
    results["cone_profile"] = {"z_m": z_fine.tolist(),
                               "w_m": w_fine.tolist(),
                               "half_side_m": grid.size_m / 2.0,
                               "R_m": {}}
    verdicts, field_c, budget_c = {}, {}, {}
    for floor_name in args.floors:
        floor = floor_name == "beam"
        smallest = None
        for c in args.factors:
            R = np.array([cone_radius_m(z, path_m, rx.aperture_m, args.cn2,
                                        lam, c, w_v, offset_m, floor,
                                        args.cone_rule)
                          for z in z_fine])
            results["cone_profile"]["R_m"][f"{floor_name}:{c:g}"] = R.tolist()
            entry = results["turbulence"][floor_name][f"{c:g}"]
            ok = (entry["pass"] and
                  results["vacuum"][floor_name][f"{c:g}"]["pass"])
            removed = entry["removed_frac_max"]
            if ok and smallest is None:
                smallest = c
            entry["R_max_m"] = float(R.max())
            entry["removes_light"] = bool(removed > 1e-9)
            entry["side_m"] = float(
                2.0 * R.max() / (1.0 - preset.boundary_width_frac))
            entry["n_same_pixel"] = int(
                2 ** np.ceil(np.log2(entry["side_m"] / grid.pixel_m)))
            entry["pass_field"] = passes_field(entry["trials"])
            entry["pass_budget"] = passes_budget(entry["trials"])
            if entry["pass_field"] and floor_name not in field_c:
                field_c[floor_name] = c
            if entry["pass_budget"] and floor_name not in budget_c:
                budget_c[floor_name] = c
        verdicts[floor_name] = smallest

    print("=" * 78)
    print("THE VERDICT")
    print("=" * 78)

    def report(floor_name, c, title):
        """Print the extent that one passing factor implies."""
        if c is None:
            print(f"  {title}: NO factor of the sweep passes.")
            return
        entry = results["turbulence"][floor_name][f"{c:g}"]
        print(f"  {title}: c = {c:g}")
        print(f"      the cone removes {entry['removed_frac_max']:.3e} of "
              f"the launched power "
              f"({'a real clip' if entry['removes_light'] else 'NOTHING'})")
        print(f"      max_z R(z)      {entry['R_max_m']:8.4f} m")
        print(f"      side it implies {entry['side_m']:8.4f} m  "
              f"(the sizer asks for {grid.size_m:.4f} m, "
              f"a factor of {grid.size_m / entry['side_m']:.2f})")
        print(f"      pixels at the SAME {grid.pixel_m * 1e3:.3f} mm pixel: "
              f"{entry['n_same_pixel']} (today {grid.n})")

    for floor_name in args.floors:
        guard = "the 3 w(z) guard ON" if floor_name == "beam" else "no guard"
        print(f"{args.cone_rule} rule, {guard}:")
        report(floor_name, verdicts.get(floor_name),
               "the smallest factor that passes ALL THREE metrics, "
               "turbulence AND vacuum")
        report(floor_name, field_c.get(floor_name),
               "the smallest factor under the FIELD rule alone "
               "(turbulence, RMS < 1e-3)")
        report(floor_name, budget_c.get(floor_name),
               "the smallest factor under the BUDGET-VISIBLE rule alone "
               "(turbulence, |dP| < 0.01 dB and |d eta| < 1e-3)")
        print("")
    results["verdict_field"] = {k: float(v) for k, v in field_c.items()}
    results["verdict_budget"] = {k: float(v) for k, v in budget_c.items()}
    results["verdict"] = {k: (None if v is None else float(v))
                          for k, v in verdicts.items()}
    print("")
    print(f"elapsed {time.time() - t_start:.1f} s")

    with open(os.path.join(HERE, f"receiver_cone_clip_{args.tag}.json"),
              "w", encoding="utf-8") as fh:
        json.dump(results, fh, indent=1)
    _figure(results, args, grid, preset)
    print(f"wrote receiver_cone_clip_{args.tag}.json, "
          f"receiver_cone_clip_{args.tag}.log and "
          f"figures/receiver_cone_clip_{args.tag}.png")
    return results


def _figure(results, args, grid, preset):
    """Draw the three panels: the field RMS, the power step, and R(z)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(15.0, 4.4))
    cs = list(args.factors)
    for floor_name in args.floors:
        label = "3 w(z) guard on" if floor_name == "beam" else "no guard"
        for ax, key, fn in ((axes[0], "field_rms", lambda r: r["field_rms"]),
                            (axes[1], "dp_db", lambda r: abs(r["dp_db"]))):
            for arm, style in (("turbulence", "-o"), ("vacuum", "--s")):
                y = [max(fn(r) for r in
                         results[arm][floor_name][f"{c:g}"]["trials"])
                     for c in cs]
                ax.plot(cs, np.maximum(y, 1e-12), style,
                        label=f"{arm}, {label}")
    axes[0].axhline(RMS_MAX, color="k", ls=":", label="pass limit")
    axes[0].set_yscale("log")
    axes[0].set_xlabel("cone factor c")
    axes[0].set_ylabel("max relative field RMS in the aperture")
    axes[0].set_title("The field inside the receive aperture")
    axes[0].legend(fontsize=7)
    axes[0].grid(alpha=0.3)

    axes[1].axhline(DP_DB_MAX, color="k", ls=":", label="pass limit")
    axes[1].set_yscale("log")
    axes[1].set_xlabel("cone factor c")
    axes[1].set_ylabel("max |collected power step| [dB]")
    axes[1].set_title("The collected power")
    axes[1].legend(fontsize=7)
    axes[1].grid(alpha=0.3)

    prof = results["cone_profile"]
    z = np.array(prof["z_m"]) / 1e3
    axes[2].plot(z, prof["w_m"], "k-", lw=2, label="vacuum beam radius w(z)")
    axes[2].axhline(prof["half_side_m"], color="r", ls="--",
                    label="grid half-side")
    axes[2].axhline((1.0 - preset.boundary_width_frac) * prof["half_side_m"],
                    color="r", ls=":", label="untouched interior")
    for name, R in prof["R_m"].items():
        floor_name, c = name.split(":")
        if floor_name != args.floors[-1]:
            continue
        axes[2].plot(z, R, lw=1, label=f"R(z), c = {c}")
    axes[2].set_xlabel("distance from the launch [km]")
    axes[2].set_ylabel("radius [m]")
    axes[2].set_title(f"The cone against the beam "
                      f"({'no guard' if args.floors[-1] == 'none' else 'guard on'})")
    axes[2].legend(fontsize=7)
    axes[2].grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, f"receiver_cone_clip_{args.tag}.png"),
                dpi=130)
    plt.close(fig)


if __name__ == '__main__':
    main()
