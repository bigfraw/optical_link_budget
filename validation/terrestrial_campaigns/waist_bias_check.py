"""Measure the bias that the CLAMPED simulation grid puts on the launch beam.

THE QUESTION. The terrestrial campaigns of
`validation/terrestrial_campaigns/run_campaigns.py` size their grid with
`olb.waveoptics.turbulence.sampling.turbulent_grid`. A strong or a long path
asks for more pixels than the preset allows, so the sizer CLAMPS the pixel
count at n_max. A clamped grid has a coarse pixel. The launch waist of the
study is 5 mm, so a coarse pixel can sample that waist with only a few points.
An under-sampled launch beam is not the beam the physics wants: it can carry
the wrong radius and the wrong power.

THE TEST. This script propagates each cell in VACUUM. It uses the SAME grid,
the SAME launch field and the SAME absorbing boundary as the campaign, but it
gives the split step FLAT (zero) screens. So the only error sources are the
grid and the mask. The vacuum answer is known in closed form, so the
difference is the bias.

It reports, for each of the twelve cells:

1. THE SAMPLING. The pixels for each launch waist, w0 / dx.
2. THE RADIUS. The second-moment beam radius at the receiver, against the
   free-space Gaussian radius w(L) of the virtual waist (`olb.beam.gaussz`).
3. THE POWER. The power inside a 10 cm bucket at the receiver, against the
   analytic Gaussian capture 1 - exp(-2 (D/2)^2 / w^2), with the difference in
   dB.

It also reports the power that reaches the receive plane at all. The
absorbing mask removes the light that comes to the edge of the grid, so a
number below 1.0 means the grid is too narrow for the beam.

It writes `waist_bias_check.log` and `waist_bias_check_results.json` next to
itself.

Run it from the repository root:

    python -m validation.terrestrial_campaigns.waist_bias_check --dry-run
    python -m validation.terrestrial_campaigns.waist_bias_check
"""

import argparse
import json
import os
import time
import warnings

import numpy as np

from olb.beam import gaussz, virtual_waist
from olb.waveoptics.field import Power, field_dtype
from olb.waveoptics.run import _clip
from olb.waveoptics.turbulence.run import _start_field
from olb.waveoptics.turbulence.sampling import PRESETS, turbulent_grid
from olb.waveoptics.turbulence.splitstep import (split_step,
                                                 super_gaussian_boundary)

from .run_campaigns import (DIVERGENCE_RAD, L0_M, LAM, PRECISION, WAIST_M,
                            all_cells, build_scenario, cell_token,
                            fill_fraction, make_say, print_table)

HERE = os.path.dirname(os.path.abspath(__file__))

# The bucket of the power test. It is the receive aperture of the study.
BUCKET_M = 0.10


def launch_waist_m(launch):
    """Give the VIRTUAL waist and the offset of the launch beam.

    A collimated launch has its waist at the aperture. A diverged launch is an
    ordinary Gaussian beam from a virtual waist w_v that sits the distance d
    behind the aperture. Source: Andrews and Phillips, 2nd ed. (2005),
    DOI 10.1117/3.626196, Ch. 4, Eqs. (7) and (8), printed p. 87.

    Args:
        launch: "collimated" or "diverged".

    Returns:
        The pair (w_v in m, offset in m).
    """
    divergence = DIVERGENCE_RAD if launch == "diverged" else None
    w_v, offset = virtual_waist(WAIST_M, divergence, LAM)
    return float(w_v), float(offset)


def second_moment_radius_m(F):
    """Give the second-moment 1/e^2 beam radius of a field, in m.

    The second radial moment of a Gaussian of 1/e^2 radius w is
    <r^2> = w^2 / 2, so w = sqrt(2 <r^2>). The irradiance profile
    I(r) = (2 / (pi w^2)) exp(-2 r^2 / w^2) is Andrews and Phillips,
    DOI 10.1117/3.626196, Ch. 4, Eq. (8), printed p. 87.

    THE MOMENT IS TAIL-SENSITIVE. The r^2 weight gives the far wings a large
    say. A truncated or an aliased beam therefore reads a radius that is too
    large. That sensitivity is the point of this test.

    Args:
        F: the Field.

    Returns:
        The radius in m.
    """
    intensity = np.abs(F.field) ** 2
    Y, X = F.mgrid_cartesian
    r2 = np.asarray(X, dtype=float) ** 2 + np.asarray(Y, dtype=float) ** 2
    total = float(intensity.sum())
    if total <= 0.0:
        return float("nan")
    return float(np.sqrt(2.0 * float((intensity * r2).sum()) / total))


def size_cell(path_m, cn2, preset, launch):
    """Size one cell exactly as the campaign sizes it.

    Args:
        path_m: the horizontal path length, in m.
        cn2:    the constant Cn2 of the path, in m^-2/3.
        preset: the quality preset name.
        launch: "collimated" or "diverged".

    Returns:
        The dict of the sized cell. It holds the scenario, the geometry, the
        grid, the screen plan, the preset object and the sizer texts.
    """
    scn, geom = build_scenario(path_m, cn2, launch)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        grid, plan, report = turbulent_grid(scn, geom, preset=preset,
                                            L0_m=L0_M)
    return {"path_m": float(path_m), "cn2": float(cn2), "preset": preset,
            "launch": launch, "scenario": scn, "geometry": geom, "grid": grid,
            "plan": plan, "quality": PRESETS[preset], "report": report,
            "warnings": sorted({str(w.message) for w in caught})}


def propagate_vacuum(cell):
    """Propagate the launch field of one cell in vacuum on its own grid.

    The recipe copies `propagate_turbulent_scenario`: the same start field,
    the same absorbing mask and the same hop plan. Every screen is FLAT, so
    the path holds no turbulence.

    Args:
        cell: a dict from `size_cell`.

    Returns:
        The receive-plane Field.
    """
    grid, plan, quality = cell["grid"], cell["plan"], cell["quality"]
    cdtype = field_dtype(PRECISION)
    F_in = _start_field(cell["scenario"], grid, LAM, is_space=False,
                        dtype=cdtype)
    mask = super_gaussian_boundary(grid.n, quality.boundary_width_frac)
    flat = np.zeros((grid.n, grid.n))
    n_screens = int(plan.z_m.size)
    F_rx = split_step(F_in, plan.z_m, [flat] * n_screens, plan.z_total_m,
                      boundary=mask)
    return F_in, F_rx


def measure_cell(cell):
    """Measure the grid bias of one cell.

    Args:
        cell: a dict from `size_cell`.

    Returns:
        A JSON-ready dict of the measurement.
    """
    t0 = time.perf_counter()
    F_in, F_rx = propagate_vacuum(cell)
    seconds = time.perf_counter() - t0

    grid = cell["grid"]
    w_v, offset = launch_waist_m(cell["launch"])
    # The free-space radius at the range L. Andrews and Phillips,
    # DOI 10.1117/3.626196, Ch. 4, Eq. (7), printed p. 87.
    w_theory = float(gaussz(w_v, offset + cell["path_m"], LAM))
    w_measured = second_moment_radius_m(F_rx)

    power_in = float(Power(F_in))
    power_rx = float(Power(F_rx))
    bucket = float(Power(_clip(F_rx, BUCKET_M, 0.0)))
    eta_measured = bucket / power_in
    eta_theory = fill_fraction(BUCKET_M, w_theory)
    # Loss is positive dB, so a measured capture BELOW the theory gives a
    # positive number of dB.
    delta_db = float(-10.0 * np.log10(eta_measured / eta_theory))

    return {
        "cell": cell_token(cell["path_m"], cell["cn2"], cell["preset"]),
        "path_length_m": cell["path_m"],
        "cn2_m_m23": cell["cn2"],
        "preset": cell["preset"],
        "launch": cell["launch"],
        "grid_n": int(grid.n),
        "grid_side_m": float(grid.size_m),
        "pixel_m": float(grid.pixel_m),
        "n_screens": int(cell["plan"].z_m.size),
        "n_clamped": (None if cell["report"] is None
                      else bool(cell["report"].n_clamped)),
        "pixels_per_waist": float(w_v / grid.pixel_m),
        "virtual_waist_m": w_v,
        "waist_offset_m": offset,
        "w_theory_m": w_theory,
        "w_second_moment_m": w_measured,
        "w_ratio": float(w_measured / w_theory),
        "power_kept": float(power_rx / power_in),
        "eta_bucket_measured": float(eta_measured),
        "eta_bucket_theory": float(eta_theory),
        "delta_db": delta_db,
        "seconds": float(seconds),
        "sizer_warnings": cell["warnings"],
    }


def result_rows(records):
    """Give the report table rows, the header first.

    Args:
        records: the list of dicts from `measure_cell`.

    Returns:
        A list of string lists.
    """
    head = ["cell", "preset", "n px", "px mm", "px/w0", "clamp", "scr",
            "w(L) cm", "w_2m cm", "w ratio", "kept", "eta sim", "eta thy",
            "d dB", "s"]
    rows = [head]
    for r in records:
        rows.append([
            r["cell"].rsplit(":", 1)[0],
            r["preset"],
            f"{r['grid_n']:d}",
            f"{r['pixel_m'] * 1e3:.2f}",
            f"{r['pixels_per_waist']:.2f}",
            ("yes" if r["n_clamped"] else "no") if r["n_clamped"] is not None
            else "-",
            f"{r['n_screens']:d}",
            f"{r['w_theory_m'] * 100:.2f}",
            f"{r['w_second_moment_m'] * 100:.2f}",
            f"{r['w_ratio']:.3f}",
            f"{r['power_kept']:.4f}",
            f"{r['eta_bucket_measured']:.4f}",
            f"{r['eta_bucket_theory']:.4f}",
            f"{r['delta_db']:+.3f}",
            f"{r['seconds']:.1f}",
        ])
    return rows


def dry_run_rows(cells):
    """Give the dry-run table rows of the sized cells, the header first."""
    head = ["cell", "preset", "n px", "side m", "px mm", "px/w0", "clamp",
            "scr", "warn"]
    rows = [head]
    for c in cells:
        w_v, _ = launch_waist_m(c["launch"])
        rows.append([
            cell_token(c["path_m"], c["cn2"], c["preset"]).rsplit(":", 1)[0],
            c["preset"],
            f"{c['grid'].n:d}",
            f"{c['grid'].size_m:.3f}",
            f"{c['grid'].pixel_m * 1e3:.2f}",
            f"{w_v / c['grid'].pixel_m:.2f}",
            ("yes" if c["report"].n_clamped else "no")
            if c["report"] is not None else "-",
            f"{c['plan'].z_m.size:d}",
            f"{len(c['warnings']):d}",
        ])
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--launch", choices=("collimated", "diverged"),
                    default="collimated",
                    help="the transmit beam. 'diverged' opens it to "
                         f"{DIVERGENCE_RAD} rad.")
    ap.add_argument("--dry-run", action="store_true",
                    help="size every cell, print the table, and propagate "
                         "nothing")
    args = ap.parse_args(argv)

    log_paths = () if args.dry_run else (os.path.join(HERE,
                                                      "waist_bias_check.log"),)
    say = make_say(*log_paths)

    say(f"waist bias check, {time.strftime('%Y-%m-%d %H:%M:%S')}")
    say(f"launch        : {args.launch}")
    say(f"fixed         : lambda {LAM * 1e9:.0f} nm, waist {WAIST_M * 1e3:.0f} mm, "
        f"bucket {BUCKET_M * 100:.0f} cm, L0 {L0_M:.0f} m, "
        f"{PRECISION} precision")
    say("test          : one VACUUM propagation for each cell, flat screens, "
        "the campaign grid and the campaign boundary mask")
    say()

    cells = [size_cell(p, c, q, args.launch) for p, c, q in all_cells()]

    if args.dry_run:
        print_table(dry_run_rows(cells), say)
        say()
        say("dry run: nothing propagated. Drop --dry-run to measure.")
        return

    records = []
    for cell in cells:
        record = measure_cell(cell)
        records.append(record)
        say(f"  {record['cell']:22s} {record['grid_n']:5d} px, "
            f"{record['pixels_per_waist']:5.2f} px per waist, "
            f"w ratio {record['w_ratio']:.3f}, "
            f"bucket {record['delta_db']:+.3f} dB, "
            f"{record['seconds']:.1f} s")
    say()
    print_table(result_rows(records), say)

    warn_lines = [f"  {r['cell']}: {w}"
                  for r in records for w in r["sizer_warnings"]]
    if warn_lines:
        say()
        say("sizer warnings (a clamped grid names the rules it breaks):")
        for line in warn_lines:
            say(line)

    path = os.path.join(HERE, "waist_bias_check_results.json")
    with open(path, "w", encoding="utf-8") as fh:
        json.dump({"launch": args.launch, "wavelength_m": LAM,
                   "waist_m": WAIST_M, "bucket_m": BUCKET_M, "L0_m": L0_M,
                   "precision": PRECISION, "cells": records}, fh, indent=2)
    say()
    say(f"outputs       : {path}")


if __name__ == "__main__":
    main()
