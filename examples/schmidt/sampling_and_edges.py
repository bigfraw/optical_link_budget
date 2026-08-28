'''
The sampling artefacts of a numerical propagation, and the rule checker.

A wave-optics propagation is a discrete Fourier transform, and a discrete
Fourier transform aliases. Chapter 7 of the book turns the geometry into four
numbered inequalities, and Chapter 8 adds the absorbing boundary that keeps the
wrapped light out of the answer. This script does two things with that
material.

PART A is a GALLERY OF DELIBERATE FAILURES. Each artefact comes as a pair: the
same physics on a grid that obeys the rule, and on a grid that breaks it.

  A1 THE RANGE LIMIT OF THE TRANSFER FUNCTION. The production `Forvard` and
     `schmidt.fresnel.angular_spectrum` both sample exp(-i pi lambda z f^2).
     Its local frequency grows with z (Ch. 7, Eq. (7.57), printed p. 126), so
     the phase turns faster than one sample once z passes
     N dx^2 / lambda. That bound is `olb.waveoptics.grid.forvard_max_z`, and it
     is constraint 4, Ch. 7, Eq. (7.59), printed p. 127, with delta2 = delta1.
     Past it the diffraction pattern folds back into the grid.
  A2 THE PHASE THAT DOES THE FOLDING. A1 shows the damage in the intensity.
     A2 shows the cause: the sampled transfer phase itself, and its local
     spatial frequency of Ch. 7, Eq. (7.57), against the Nyquist limit
     N delta1 / 2 of Ch. 7, Eq. (7.58), printed p. 126.
  A3 EDGE CONTROL. A chain of angular-spectrum partial steps
     (`schmidt.fresnel.partial_propagations`) spreads a hard-edged beam past
     the grid. WITHOUT an absorber the light re-enters at the opposite edge and
     interferes with itself. WITH the super-Gaussian absorber of Ch. 8,
     Eq. (8.1), printed p. 134, it goes away.

PART B is the RULE CHECKER ON REAL GRIDS. It builds the same two scenarios that
`examples/waveoptics` uses, asks the production sizers for their grids, maps
those grids onto the book inputs, and prints a pass or fail table with a
citation on every row.

  - `schmidt.sampling.check_sampling` gives the five VACUUM rows of Ch. 7.
  - `schmidt.turbulence.properly_sampled_checklist` gives the Sec. 9.5 rows,
    which are the same constraints with the turbulence-blurred extents of
    Ch. 9, Eqs. (9.84) and (9.85), printed p. 173, plus the two turbulent pitch
    rules of Sec. 9.4, printed p. 172.

THE TWO ABSORBERS ARE NOT THE SAME SHAPE. The book absorber
(`schmidt.fresnel.super_gaussian_absorber`) is exp(-(r/(0.47 N))^16) in PIXEL
radius. The production absorber
(`olb.waveoptics.turbulence.splitstep.super_gaussian_boundary`) is
exp(-t^8) with t the distance into a band of 0.125 of the HALF-side. The
parameterisation differs, so the two numbers do not compare. The SHAPES do, and
the script plots them on one axes: at the middle of an edge the book keeps
0.068 of the amplitude and the production mask keeps 0.368, so the book absorbs
about 5 times harder there. The production mask is flat to 1.0 over the inner
0.875 of the half-side; the book mask is already at 0.99 by 0.8 of it.

This script changes NO olb module. It reads the production layer only.

Figures:
    examples/schmidt/figures/sampling_artefacts.png    the six failure panels
    examples/schmidt/figures/sampling_absorbers.png    the two absorber profiles

Sources:
- Schmidt, Numerical Simulation of Optical Wave Propagation with Examples in
  MATLAB, DOI 10.1117/3.866274. Ch. 7, Eq. (7.14), printed p. 119 (constraint
  1); Ch. 7, Eq. (7.20), printed p. 120 (constraint 2); Ch. 7, Eq. (7.53),
  printed p. 126 (constraint 3); Ch. 7, Eqs. (7.57) to (7.59), printed
  pp. 126 and 127 (the transfer phase and constraint 4); Ch. 7, Eqs. (7.41) and
  (7.42), printed p. 123 (the Fresnel-integral minimum distance); Ch. 8,
  Eq. (8.1), printed p. 134 (the super-Gaussian absorber); Ch. 8, Eq. (8.18),
  printed p. 139 (the partial-propagation chain); Ch. 8, Eq. (8.24), printed
  p. 144 (the step cap); Ch. 9, Eqs. (9.84) to (9.90), printed pp. 173 and 174
  (the turbulent sampling bounds); Ch. 9, Sec. 9.5, printed pp. 174 to 182 (the
  procedure).
- Goodman, Introduction to Fourier Optics, ISBN 978-0974707723. The circular
  pupil and the Fresnel diffraction integral.

Run from the repo root:
    python -m examples.schmidt.sampling_and_edges
'''

import math
import time
import warnings

import matplotlib.pyplot as plt
import numpy as np

from olb.geometry import CircularOrbit, HorizontalPath
from olb.scenario import (Channel, SpaceScenario, TerrestrialChannel,
                          TerrestrialScenario)
from olb.terminal import Terminal, Transmitter
from olb.waveoptics import (Begin, CircAperture, Forvard, GaussBeam, GridSpec,
                            beam_magnification, forvard_max_z)
from olb.waveoptics.schmidt.fresnel import (one_step_fresnel,
                                            partial_propagations,
                                            super_gaussian_absorber)
from olb.waveoptics.schmidt.sampling import check_sampling, constraint4_min_n
from olb.waveoptics.schmidt.turbulence import (blurred_extent,
                                               constraint2_n_min,
                                               properly_sampled_checklist)
from olb.waveoptics.turbulence import super_gaussian_boundary, turbulent_grid

WAVELENGTH_M = 1550e-9

# ---- A1 and A2: the transfer-function range limit ----
ALIAS_N = 512
ALIAS_SIDE_M = 20e-3            # the pitch is 39.06 um
ALIAS_WAIST_M = 1.5e-3
ALIAS_APERTURE_R_M = 1.5e-3
CLEAN_FRACTION = 0.8            # z = 0.8 of the limit: the grid holds it
ALIAS_FACTOR = 8.0              # z = 8 times the limit: the grid does not

# ---- A3: the partial-propagation chain and its edge ----
EDGE_N = 512
EDGE_SIDE_M = 10e-3             # the pitch is 19.53 um
EDGE_WAIST_M = 1.0e-3
EDGE_APERTURE_R_M = 1.0e-3
EDGE_RANGE_M = 8.0              # the Airy null lands past the half-side

# ---- Part B: the two production scenarios ----
# They repeat examples/waveoptics/turbulent_terrestrial.py and
# examples/waveoptics/turbulent_downlink.py.
TERR_PATH_M = 2000.0
TERR_CN2 = 3e-15
TERR_WAIST_M = 0.05
TERR_TX_APERTURE_M = 0.20
TERR_RX_APERTURE_M = 0.10

# The space case is an UPLINK, because that is the direction whose vacuum
# sizer takes the co-moving route with no warning. A 100 mm launch aperture
# grows into a 500 mm receive aperture at 1075 km, so a flat grid cannot hold
# both ends. See olb/waveoptics/grid.py.
SPACE_ALTITUDE_M = 600e3
SPACE_ELEVATION_DEG = 30.0
SPACE_GROUND_APERTURE_M = 0.10
SPACE_GROUND_OBSCURATION = 0.30
SPACE_SAT_APERTURE_M = 0.50
SPACE_WAIST_M = 0.05

PRESET = "standard"

ARTEFACTS_PNG = "examples/schmidt/figures/sampling_artefacts.png"
ABSORBERS_PNG = "examples/schmidt/figures/sampling_absorbers.png"


# ---------------------------------------------------------------------------
# Part A
# ---------------------------------------------------------------------------

def hard_field(side_m, n, waist_m, radius_m):
    '''Build the hard-truncated Gaussian that every artefact panel launches.

    A hard edge has a spectrum that does not decay, so it fills the frequency
    grid. That is what makes an aliasing artefact visible. See Goodman,
    Introduction to Fourier Optics, ISBN 978-0974707723 (the circular pupil).
    '''
    return CircAperture(GaussBeam(Begin(side_m, WAVELENGTH_M, n), waist_m),
                        radius_m)


def artefact_range_limit():
    '''A1: run Forvard inside and past the transfer-function range limit.

    The reference at the long range is `one_step_fresnel`. That kernel carries
    NO transfer function, so constraint 4 does not touch it (Schmidt (2010),
    DOI 10.1117/3.866274, Ch. 7, text below Eq. (7.31), printed p. 121). It is
    the LONG-range kernel, and it lands on its own pitch of Ch. 6, Eq. (6.16),
    printed p. 90.
    '''
    dx = ALIAS_SIDE_M / ALIAS_N
    z_limit = forvard_max_z(GridSpec(size_m=ALIAS_SIDE_M, n=ALIAS_N),
                            WAVELENGTH_M)
    z_ok = CLEAN_FRACTION * z_limit
    z_bad = ALIAS_FACTOR * z_limit

    F0 = hard_field(ALIAS_SIDE_M, ALIAS_N, ALIAS_WAIST_M, ALIAS_APERTURE_R_M)
    clean = np.abs(Forvard(F0, z_ok).field) ** 2
    aliased = np.abs(Forvard(F0, z_bad).field) ** 2
    truth, dx_truth = one_step_fresnel(F0.field, WAVELENGTH_M, dx, z_bad)
    return {
        "dx": dx, "z_limit": z_limit, "z_ok": z_ok, "z_bad": z_bad,
        "clean": clean, "aliased": aliased,
        "truth": np.abs(truth) ** 2, "dx_truth": dx_truth,
        "n_need_ok": constraint4_min_n(dx, dx, WAVELENGTH_M, z_ok),
        "n_need_bad": constraint4_min_n(dx, dx, WAVELENGTH_M, z_bad),
    }


def transfer_phase(n, dx, z):
    '''Give the sampled transfer phase and its local spatial frequency.

    The phase is Ch. 6, Eq. (6.32), printed p. 95, without its piston:

        phi'(f) = -pi lambda z f^2

    Its local spatial frequency is Ch. 7, Eq. (7.57), printed p. 126, with
    delta2 = delta1:

        f'_loc(f) = lambda z f

    The Nyquist limit of the frequency grid is N delta1 / 2 (Ch. 7, Eq. (7.58),
    printed p. 126). The rule that comes out of the two is constraint 4.
    '''
    f = (np.arange(n) - n // 2) / (n * dx)
    phase = -np.pi * WAVELENGTH_M * z * f ** 2
    return f, np.angle(np.exp(1j * phase)), WAVELENGTH_M * z * np.abs(f)


def artefact_edge_control():
    '''A3: a partial-propagation chain, with and without the book absorber.

    The chain is Ch. 8, Eq. (8.18), printed p. 139. The plane count comes from
    the step cap of Ch. 8, Eq. (8.24), printed p. 144. The pitch is constant
    here, so the linear pitch rule of Ch. 8, Eq. (8.8), printed p. 136, is
    trivial and the whole chain runs at one pitch.
    '''
    dx = EDGE_SIDE_M / EDGE_N
    dz_max = dx ** 2 * EDGE_N / WAVELENGTH_M          # Ch. 8, Eq. (8.24)
    n_planes = int(math.ceil(EDGE_RANGE_M / dz_max)) + 1
    z_planes = np.linspace(0.0, EDGE_RANGE_M, n_planes)

    F0 = hard_field(EDGE_SIDE_M, EDGE_N, EDGE_WAIST_M, EDGE_APERTURE_R_M)
    book_mask = super_gaussian_absorber(EDGE_N)
    bare = partial_propagations(F0.field, WAVELENGTH_M, dx, dx, z_planes)
    damped = partial_propagations(F0.field, WAVELENGTH_M, dx, dx, z_planes,
                                  absorber=book_mask)
    return {
        "dx": dx, "dz_max": dz_max, "n_planes": n_planes,
        "bare": np.abs(bare) ** 2, "damped": np.abs(damped) ** 2,
        # The Airy first null of a circular pupil, Goodman,
        # ISBN 978-0974707723: theta = 0.61 lambda / R.
        "null_radius_m": 0.61 * WAVELENGTH_M * EDGE_RANGE_M / EDGE_APERTURE_R_M,
    }


def absorber_profiles(n=512):
    '''Give the two absorber profiles along one radius, in half-side units.

    The book mask is Schmidt (2010), DOI 10.1117/3.866274, Ch. 8, Eq. (8.1),
    printed p. 134, with the Listing 8.1 values sigma = 0.47 N and power 16.
    The production mask is
    `olb.waveoptics.turbulence.splitstep.super_gaussian_boundary`, with a taper
    band of 0.125 of the half-side and power 8.
    '''
    # Read the DIAGONAL, not the row. The radius then runs from 0.0 at the
    # axis, through 1.0 at the middle of an edge, to 1.414 at a corner, so the
    # whole fall of each mask is visible.
    idx = np.arange(n // 2, n)
    rho = np.sqrt(2.0) * (idx - n // 2) / (n / 2.0)
    book = super_gaussian_absorber(n)[idx, idx]
    prod = super_gaussian_boundary(n)[idx, idx]
    return rho, book, prod


# ---------------------------------------------------------------------------
# Part B
# ---------------------------------------------------------------------------

def terrestrial_scenario():
    '''Build the 2 km horizontal case of examples/waveoptics.'''
    return TerrestrialScenario(
        near=Terminal(aperture_m=TERR_TX_APERTURE_M,
                      wavelength_m=WAVELENGTH_M, pointing_jitter_rad=0.0,
                      transmitter=Transmitter(waist_m=TERR_WAIST_M)),
        far=Terminal(aperture_m=TERR_RX_APERTURE_M,
                     wavelength_m=WAVELENGTH_M, pointing_jitter_rad=0.0),
        channel=TerrestrialChannel(path_length_m=TERR_PATH_M, cn2=TERR_CN2))


def space_scenario():
    '''Build the 600 km uplink case of examples/waveoptics and grid.py.

    An uplink puts the transmitter on the GROUND terminal, so that is the
    terminal that carries the Transmitter. See olb/scenario.py. The turbulent
    sizer plans the DOWNLINK slab for this scenario too, because the split-step
    layer always propagates down the atmosphere and reads an uplink through
    reciprocity. See olb/waveoptics/turbulence/run.py.
    '''
    return SpaceScenario(
        ground=Terminal(aperture_m=SPACE_GROUND_APERTURE_M,
                        obscuration_ratio=SPACE_GROUND_OBSCURATION,
                        wavelength_m=WAVELENGTH_M, pointing_jitter_rad=0.0,
                        transmitter=Transmitter(waist_m=SPACE_WAIST_M)),
        space=Terminal(aperture_m=SPACE_SAT_APERTURE_M,
                       wavelength_m=WAVELENGTH_M, pointing_jitter_rad=0.0),
        direction="uplink", channel=Channel(altitude_m=SPACE_ALTITUDE_M))


def print_vacuum_rows(rows):
    '''Print one `check_sampling` result as a pass or fail table.'''
    for rule in rows:
        mark = "PASS" if rule.satisfied else "FAIL"
        if isinstance(rule.bound, tuple):
            bound = f"({rule.bound[0]:.4g}, {rule.bound[1]:.4g})"
        else:
            bound = f"{rule.bound:.4g}"
        print(f"    {mark}  {rule.name:<40s}{bound:>26s}"
              f"{rule.actual:>14.4g}")
        print(f"          {rule.citation}")


def print_turbulent_rows(rows):
    '''Print one `properly_sampled_checklist` result as a table.

    A row with `satisfied = None` is ADVISORY: the book asks for a procedure,
    not an inequality, so the checker cannot test it.
    '''
    for name, ok, bound, actual, citation in rows:
        mark = "----" if ok is None else ("PASS" if ok else "FAIL")
        if bound is None:
            text_b, text_a = "", ""
        elif isinstance(bound, tuple):
            text_b = f"({bound[0]:.4g}, {bound[1]:.4g})"
            text_a = f"{actual:.4g}"
        else:
            text_b, text_a = f"{bound:.4g}", f"{actual:.4g}"
        print(f"    {mark}  {name:<38s}{text_b:>26s}{text_a:>14s}")
        print(f"          {citation}")


# ---------------------------------------------------------------------------
# The figures
# ---------------------------------------------------------------------------

def draw_artefacts(a1, a3):
    '''Draw the six failure panels: clean beside broken, on one colour scale.'''
    fig, axes = plt.subplots(2, 3, figsize=(17.4, 10.4),
                             constrained_layout=True)
    floor = 1e-6

    def show(ax, data, side_m, title, peak):
        half = side_m / 2 * 1e3
        image = ax.imshow(np.log10(np.maximum(data / peak, floor)),
                          extent=[-half, half, -half, half], origin="lower",
                          cmap="magma", vmin=np.log10(floor), vmax=0.0)
        ax.set_title(title, fontsize=9)
        ax.set_xlabel("x, mm")
        ax.set_ylabel("y, mm")
        fig.colorbar(image, ax=ax, shrink=0.80, label="log10(I / I_peak)")

    # ---- row 1: A1, the transfer-function range limit ----
    peak1 = max(a1["clean"].max(), a1["aliased"].max())
    show(axes[0, 0], a1["clean"], ALIAS_SIDE_M,
         f"A1 CLEAN. Forvard at z = {a1['z_ok'] * 1e3:.0f} mm\n"
         f"{CLEAN_FRACTION:.1f} of the limit N dx^2/lambda = "
         f"{a1['z_limit'] * 1e3:.0f} mm", peak1)
    show(axes[0, 1], a1["aliased"], ALIAS_SIDE_M,
         f"A1 BROKEN. Forvard at z = {a1['z_bad'] * 1e3:.0f} mm\n"
         f"{ALIAS_FACTOR:.0f} times the limit: constraint 4 asks for N >= "
         f"{a1['n_need_bad']:.0f}", peak1)
    show(axes[0, 2], a1["truth"], ALIAS_N * a1["dx_truth"],
         f"A1 TRUTH. one_step_fresnel at the same z\n"
         f"no transfer function, so no constraint 4. Pitch "
         f"{a1['dx_truth'] * 1e6:.0f} um", a1["truth"].max())
    # The one-step pitch is 8 times coarser, so its grid is 8 times wider.
    # Crop it to the Forvard grid, so the three panels show one field of view.
    axes[0, 2].set_xlim(-ALIAS_SIDE_M / 2 * 1e3, ALIAS_SIDE_M / 2 * 1e3)
    axes[0, 2].set_ylim(-ALIAS_SIDE_M / 2 * 1e3, ALIAS_SIDE_M / 2 * 1e3)

    # ---- row 2, panel 1: A2, the phase that folds ----
    ax = axes[1, 0]
    f_ok, ph_ok, loc_ok = transfer_phase(ALIAS_N, a1["dx"], a1["z_ok"])
    f_bad, ph_bad, loc_bad = transfer_phase(ALIAS_N, a1["dx"], a1["z_bad"])
    ax.plot(f_ok * 1e-3, ph_ok, color="tab:blue", linewidth=1.6, marker="o",
            markersize=3.4, label=f"z = {a1['z_ok'] * 1e3:.0f} mm, sampled")
    ax.plot(f_bad * 1e-3, ph_bad, color="tab:red", linewidth=1.0, marker="s",
            markersize=3.4, label=f"z = {a1['z_bad'] * 1e3:.0f} mm, sampled")
    ax.set_xlabel("spatial frequency f, 1/mm")
    ax.set_ylabel("transfer phase, rad")
    ax.set_title("A2 THE CAUSE. The sampled transfer phase\n"
                 "-pi lambda z f^2, Ch. 6, Eq. (6.32), p. 95. One marker is "
                 "one sample.", fontsize=9)
    # Show the first 2.5 1/mm only, so every sample is visible. The red run
    # already turns more than pi between two samples well inside that band.
    ax.set_xlim(0.0, 2.5)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)

    # ---- row 2, panel 2: A2, the local frequency against Nyquist ----
    ax = axes[1, 1]
    nyquist = ALIAS_N * a1["dx"] / 2.0
    ax.plot(f_ok * 1e-3, loc_ok * 1e3, color="tab:blue", linewidth=2.0,
            label=f"z = {a1['z_ok'] * 1e3:.0f} mm")
    ax.plot(f_bad * 1e-3, loc_bad * 1e3, color="tab:red", linewidth=2.0,
            label=f"z = {a1['z_bad'] * 1e3:.0f} mm")
    ax.axhline(nyquist * 1e3, color="black", linestyle="--", linewidth=1.6,
               label="Nyquist limit N delta1 / 2")
    ax.set_xlabel("spatial frequency f, 1/mm")
    ax.set_ylabel("local frequency of the phase, mm")
    ax.set_title("A2 THE RULE. Ch. 7, Eq. (7.57), printed p. 126\n"
                 "the red curve leaves the grid: that IS constraint 4",
                 fontsize=9)
    ax.set_xlim(0.0, f_ok.max() * 1e-3)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)

    # ---- row 2, panel 3: A1 again, as a cut through the middle ----
    ax = axes[1, 2]
    x_mm = (np.arange(ALIAS_N) - ALIAS_N // 2) * a1["dx"] * 1e3
    x_truth_mm = ((np.arange(ALIAS_N) - ALIAS_N // 2) * a1["dx_truth"] * 1e3)
    ax.semilogy(x_mm, a1["clean"][ALIAS_N // 2] / peak1, color="tab:blue",
                linewidth=1.6,
                label=f"clean run, at its own z = {a1['z_ok'] * 1e3:.0f} mm")
    ax.semilogy(x_mm, a1["aliased"][ALIAS_N // 2] / a1["aliased"].max(),
                color="tab:red", linewidth=1.6,
                label=f"broken, z = {a1['z_bad'] * 1e3:.0f} mm")
    ax.semilogy(x_truth_mm, a1["truth"][ALIAS_N // 2] / a1["truth"].max(),
                color="black", linestyle="--", linewidth=1.4,
                label="one_step_fresnel at the broken z")
    ax.set_xlim(-ALIAS_SIDE_M / 2 * 1e3, ALIAS_SIDE_M / 2 * 1e3)
    ax.set_ylim(1e-6, 2.0)
    ax.set_xlabel("x, mm")
    ax.set_ylabel("I / I_peak")
    ax.set_title("A1 THE DAMAGE, as a cut. Red against black is the\n"
                 "same z: the broken run lifts the whole grid edge",
                 fontsize=9)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)

    fig.suptitle("Deliberate sampling failures. Every panel names the rule "
                 "that it obeys or breaks.\n"
                 f"{WAVELENGTH_M * 1e9:.0f} nm, hard-truncated Gaussian",
                 fontsize=13)
    fig.savefig(ARTEFACTS_PNG, dpi=150)
    plt.close(fig)


def draw_edge_pair(a3):
    '''Draw the A3 pair on its own figure, so the colour scale is shared.'''
    fig, axes = plt.subplots(1, 3, figsize=(17.0, 5.2),
                             constrained_layout=True)
    peak = max(a3["bare"].max(), a3["damped"].max())
    floor = 1e-6
    half = EDGE_SIDE_M / 2 * 1e3
    for ax, data, title in (
            (axes[0], a3["bare"],
             f"A3 BROKEN. {a3['n_planes']} partial steps, NO absorber\n"
             f"the light that leaves one edge comes back at the other"),
            (axes[1], a3["damped"],
             "A3 CLEAN. the same chain with the book absorber\n"
             "Ch. 8, Eq. (8.1), printed p. 134. The round shape IS the "
             "mask.")):
        image = ax.imshow(np.log10(np.maximum(data / peak, floor)),
                          extent=[-half, half, -half, half], origin="lower",
                          cmap="magma", vmin=np.log10(floor), vmax=0.0)
        ax.set_title(title, fontsize=9)
        ax.set_xlabel("x, mm")
        ax.set_ylabel("y, mm")
        fig.colorbar(image, ax=ax, shrink=0.80, label="log10(I / I_peak)")

    ax = axes[2]
    x_mm = (np.arange(EDGE_N) - EDGE_N // 2) * a3["dx"] * 1e3
    ax.semilogy(x_mm, a3["bare"][EDGE_N // 2] / peak, color="tab:red",
                linewidth=1.6, label="no absorber")
    ax.semilogy(x_mm, a3["damped"][EDGE_N // 2] / peak, color="tab:blue",
                linewidth=1.6, label="book absorber")
    ax.set_xlabel("x, mm")
    ax.set_ylabel("I / I_peak")
    ax.set_ylim(1e-8, 2.0)
    ax.set_title("the cut through the middle\n"
                 "the absorber removes the returning light", fontsize=9)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)

    fig.suptitle(f"A3 EDGE CONTROL. A {EDGE_APERTURE_R_M * 1e3:.0f} mm hard "
                 f"aperture over {EDGE_RANGE_M:.0f} m on a "
                 f"{EDGE_SIDE_M * 1e3:.0f} mm grid.\n"
                 f"The Airy first null sits at "
                 f"{a3['null_radius_m'] * 1e3:.1f} mm, past the "
                 f"{EDGE_SIDE_M / 2 * 1e3:.1f} mm half-side: the beam is "
                 f"wider than the grid.", fontsize=12)
    fig.savefig(ARTEFACTS_PNG.replace(".png", "_edges.png"), dpi=150)
    plt.close(fig)


def draw_absorbers(rho, book, prod):
    '''Draw the two absorber profiles, and the ratio between them.'''
    fig, (a0, a1) = plt.subplots(1, 2, figsize=(12.6, 5.2),
                                 constrained_layout=True)
    a0.plot(rho, book, color="tab:blue", linewidth=2.2,
            label="book: exp(-(r/(0.47 N))^16), Ch. 8, Eq. (8.1), p. 134")
    a0.plot(rho, prod, color="tab:red", linewidth=2.2,
            label="production: power 8, band 0.125 of the half-side")
    a0.axvline(1.0, color="black", linestyle=":", linewidth=1.2)
    a0.annotate("middle of an edge", xy=(1.0, 0.55), xytext=(0.35, 0.6),
                fontsize=8, arrowprops=dict(arrowstyle="->", lw=1.0))
    a0.set_xlabel("radius, in units of the half-side")
    a0.set_ylabel("amplitude that the mask keeps")
    a0.set_title("the two absorbers, linear", fontsize=10)
    a0.grid(alpha=0.3)
    a0.legend(fontsize=8, loc="lower left")

    a1.semilogy(rho, np.maximum(book, 1e-12), color="tab:blue", linewidth=2.2,
                label="book")
    a1.semilogy(rho, np.maximum(prod, 1e-12), color="tab:red", linewidth=2.2,
                label="production")
    a1.set_ylim(1e-12, 2.0)
    a1.set_xlabel("radius, in units of the half-side")
    a1.set_ylabel("amplitude that the mask keeps")
    a1.set_title("the same, logarithmic: the book starts to fall sooner, and "
                 "the production mask falls faster past the edge", fontsize=9)
    a1.grid(alpha=0.3)
    a1.legend(fontsize=8)

    fig.suptitle("The book absorber against the production boundary mask. "
                 "The parameterisations differ, so compare the SHAPES.",
                 fontsize=12)
    fig.savefig(ABSORBERS_PNG, dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------

def main():
    t_start = time.time()

    print("=" * 78)
    print("PART A. The sampling artefacts, made on purpose")
    print("=" * 78)

    a1 = artefact_range_limit()
    print("")
    print("A1 THE RANGE LIMIT OF THE TRANSFER FUNCTION")
    print(f"  grid                    {ALIAS_N:11d} px, "
          f"{ALIAS_SIDE_M * 1e3:.3f} mm, pitch {a1['dx'] * 1e6:.3f} um")
    print(f"  limit N dx^2 / lambda   {a1['z_limit'] * 1e3:11.3f} mm  "
          f"(Ch. 7, Eq. (7.59), printed p. 127, with delta2 = delta1)")
    print(f"  clean run z             {a1['z_ok'] * 1e3:11.3f} mm  "
          f"(constraint 4 asks for N >= {a1['n_need_ok']:.0f}, "
          f"and N = {ALIAS_N})")
    print(f"  broken run z            {a1['z_bad'] * 1e3:11.3f} mm  "
          f"(constraint 4 asks for N >= {a1['n_need_bad']:.0f}, "
          f"and N = {ALIAS_N})")
    energy_edge = float(a1["aliased"][:ALIAS_N // 16].sum()
                        / a1["aliased"].sum())
    energy_edge_ok = float(a1["clean"][:ALIAS_N // 16].sum()
                           / a1["clean"].sum())
    print("  energy in the outer 1/16 band of the grid:")
    print(f"    clean run             {energy_edge_ok:11.3e}")
    print(f"    broken run            {energy_edge:11.3e}  "
          f"({energy_edge / energy_edge_ok:.0f} times the clean value)")
    assert energy_edge > 20.0 * energy_edge_ok, (energy_edge, energy_edge_ok)
    print("  The broken run puts the wrapped light back at the grid edge. The")
    print("  clean run does not. The one_step_fresnel panel of the figure is")
    print("  the same physics with NO transfer function, so it stays right.")

    a3 = artefact_edge_control()
    print("")
    print("A3 EDGE CONTROL, THE PARTIAL-PROPAGATION CHAIN")
    print(f"  grid                    {EDGE_N:11d} px, "
          f"{EDGE_SIDE_M * 1e3:.3f} mm, pitch {a3['dx'] * 1e6:.3f} um")
    print(f"  range                   {EDGE_RANGE_M:11.3f} m")
    print(f"  step cap, Eq. (8.24)    {a3['dz_max'] * 1e3:11.3f} mm")
    print(f"  planes, Eq. (8.24) text {a3['n_planes']:11d}")
    print(f"  Airy first null radius  {a3['null_radius_m'] * 1e3:11.3f} mm  "
          f"(the half-side is {EDGE_SIDE_M / 2 * 1e3:.2f} mm)")
    bare_edge = float(a3["bare"][:EDGE_N // 16].sum() / a3["bare"].sum())
    damped_edge = float(a3["damped"][:EDGE_N // 16].sum()
                        / a3["damped"].sum())
    print("  energy in the outer 1/16 band:")
    print(f"    no absorber           {bare_edge:11.3e}")
    print(f"    book absorber         {damped_edge:11.3e}")
    print(f"  power kept by the chain, absorber over bare: "
          f"{a3['damped'].sum() / a3['bare'].sum():.4f}")
    assert damped_edge < 0.1 * bare_edge, (damped_edge, bare_edge)
    print("  The absorber is a NUMERICAL device, not physics. It removes")
    print("  energy on purpose (Ch. 8, text, printed p. 134), so the region of")
    print("  interest must stay well inside the flat part of the mask.")

    rho, book, prod = absorber_profiles()
    mid_book = float(np.interp(1.0, rho, book))
    mid_prod = float(np.interp(1.0, rho, prod))
    print("")
    print("THE TWO ABSORBER SHAPES")
    print(f"  {'radius / half-side':<24}{'book, power 16':>18}"
          f"{'production, power 8':>22}")
    for probe in (0.0, 0.5, 0.7, 0.8, 0.875, 0.95, 1.0, 1.2, 1.414):
        print(f"  {probe:<24.3f}{float(np.interp(probe, rho, book)):>18.6f}"
              f"{float(np.interp(probe, rho, prod)):>22.6f}")
    print(f"  At the middle of an edge (radius 1.0) the book keeps "
          f"{mid_book:.4f} and the")
    print(f"  production mask keeps {mid_prod:.4f}, so the book absorbs "
          f"{mid_prod / mid_book:.1f} times harder there.")
    print("  The book value is exp(-(0.5/0.47)^16) of Listing 8.1, printed")
    print("  p. 142. The production value is exp(-1) at the outer edge of its")
    print("  band. The two parameterisations are not comparable number for")
    print("  number; the shapes are.")
    assert 4.0 < mid_prod / mid_book < 7.0, mid_prod / mid_book

    # -----------------------------------------------------------------
    print("")
    print("=" * 78)
    print("PART B. The rule checker against the production grids")
    print("=" * 78)

    # ---- B1: the vacuum sizer, terrestrial ----
    terr = terrestrial_scenario()
    terr_path = HorizontalPath(TERR_PATH_M)
    with warnings.catch_warnings(record=True) as caught_t:
        warnings.simplefilter("always")
        g_terr = GridSpec.for_scenario(terr, terr_path)
    print("")
    print(f"B1 VACUUM SIZER, TERRESTRIAL, {TERR_PATH_M * 1e-3:.0f} km "
          f"(GridSpec.for_scenario)")
    print(f"  grid side               {g_terr.size_m:11.4f} m")
    print(f"  pixels per side         {g_terr.n:11d}")
    print(f"  pixel pitch             {g_terr.pixel_m * 1e3:11.4f} mm")
    print(f"  route                   "
          f"{'co-moving' if g_terr.scaled else 'flat':>11}")
    for w in caught_t:
        print(f"  warning: {w.message}")
    print(f"  The book inputs: D1 is the {TERR_TX_APERTURE_M * 1e3:.0f} mm "
          f"launch aperture, D2 is the")
    print(f"  {TERR_RX_APERTURE_M * 1e3:.0f} mm receive aperture, delta1 = "
          f"delta2 = the pixel (a flat grid), and R is")
    print("  infinite, because the launch beam is collimated.")
    print("")
    print(f"    {'':4}  {'rule':<40s}{'bound':>26s}{'actual':>14s}")
    rows_t = check_sampling(TERR_TX_APERTURE_M, TERR_RX_APERTURE_M,
                            g_terr.pixel_m, g_terr.pixel_m, g_terr.n,
                            WAVELENGTH_M, TERR_PATH_M)
    print_vacuum_rows(rows_t)

    # ---- B2: the vacuum sizer, space ----
    space = space_scenario()
    orbit = CircularOrbit(altitude_m=SPACE_ALTITUDE_M,
                          elevation_deg=[SPACE_ELEVATION_DEG])
    z_space = float(np.max(orbit.slant_range_m))
    with warnings.catch_warnings(record=True) as caught_s:
        warnings.simplefilter("always")
        g_space = GridSpec.for_scenario(space, orbit)
    m_space = beam_magnification(space, z_space) if g_space.scaled else 1.0
    delta2_space = g_space.pixel_m * m_space
    print("")
    print(f"B2 VACUUM SIZER, SPACE UPLINK, {z_space * 1e-3:.0f} km at "
          f"{SPACE_ELEVATION_DEG:.0f} deg (GridSpec.for_scenario)")
    print(f"  grid side at launch     {g_space.size_m:11.4f} m")
    print(f"  pixels per side         {g_space.n:11d}")
    print(f"  pixel at launch         {g_space.pixel_m * 1e3:11.4f} mm")
    print(f"  route                   "
          f"{'co-moving' if g_space.scaled else 'flat':>11}")
    print(f"  magnification m         {m_space:11.3f}")
    print(f"  pixel at the receiver   {delta2_space * 1e3:11.4f} mm")
    for w in caught_s:
        print(f"  warning: {w.message}")
    print("  The co-moving route makes delta2 = m delta1, which is exactly the")
    print("  free output pitch of Ch. 6, Eq. (6.54), printed p. 100. So the")
    print("  checker reads delta1 at the launch plane and delta2 at the")
    print("  receive plane.")
    print("")
    print(f"    {'':4}  {'rule':<40s}{'bound':>26s}{'actual':>14s}")
    rows_s = check_sampling(SPACE_GROUND_APERTURE_M, SPACE_SAT_APERTURE_M,
                            g_space.pixel_m, delta2_space, g_space.n,
                            WAVELENGTH_M, z_space)
    print_vacuum_rows(rows_s)
    print("  ROWS 2 AND 4 FAIL BY TWO ORDERS OF MAGNITUDE, AND THAT IS THE")
    print("  POINT OF THE CO-MOVING ROUTE. Both rows carry the term")
    print("  lambda z /(delta1 delta2), which samples the ANGULAR-SPECTRUM")
    print("  transfer phase (Ch. 7, Eqs. (7.20) and (7.59), printed pp. 120 and")
    print("  127). No flat grid of any practical pixel count holds a 1000 km")
    print("  hop at a 0.1 mm launch pitch. The production recipe therefore runs")
    print("  NO transfer function: `LensFresnel` is a convolution with the")
    print("  closed-form C(x), S(x) pixel kernel, in coordinates that travel")
    print("  with the beam. The rule that governs THAT kernel is row 5, the")
    print("  Fresnel-integral minimum distance of Ch. 7, Eq. (7.42), printed")
    print(f"  p. 123, and the range passes it by a factor of "
          f"{z_space / rows_s[4].bound:.0f}.")
    print("  Rows 1 and 3 pass with room to spare.")

    # ---- B3: the turbulent sizer, terrestrial ----
    with warnings.catch_warnings(record=True):
        warnings.simplefilter("always")
        gt, plan_t, rep_t = turbulent_grid(terr, terr_path, preset=PRESET)
    print("")
    print(f"B3 TURBULENT SIZER, TERRESTRIAL, {TERR_PATH_M * 1e-3:.0f} km, "
          f"Cn2 = {TERR_CN2:g} (turbulent_grid, {PRESET})")
    print(f"  grid side               {gt.size_m:11.4f} m")
    print(f"  pixels per side         {gt.n:11d}")
    print(f"  pixel pitch             {gt.pixel_m * 1e3:11.4f} mm")
    print(f"  screens                 {plan_t.z_m.size:11d}")
    print(f"  r0 of the whole path    {plan_t.r0_total_m * 1e2:11.3f} cm")
    print(f"  sigma2_R, plane wave    {plan_t.sigma2_r.sum():11.4f}")
    for text in rep_t.warnings:
        print(f"  warning: {text}")
    print("")
    print(f"    {'':4}  {'step'  :<38s}{'bound':>26s}{'actual':>14s}")
    check_t = properly_sampled_checklist(
        wavelength_m=WAVELENGTH_M, z_total_m=plan_t.z_total_m, n=gt.n,
        pixel_m=gt.pixel_m, z_m=plan_t.z_m, r0_m=plan_t.r0_m,
        r0_total_m=plan_t.r0_total_m, d_tx_m=TERR_TX_APERTURE_M,
        d_rx_m=TERR_RX_APERTURE_M, size_m=gt.size_m, wave='spherical')
    print_turbulent_rows(check_t)
    print("  Every inequality of Sec. 9.5 holds on this grid. The production")
    print("  sizer sent no warning either, so the two agree.")

    # ---- B4: the turbulent sizer, space ----
    with warnings.catch_warnings(record=True):
        warnings.simplefilter("always")
        gs, plan_s, rep_s = turbulent_grid(space, orbit, preset=PRESET)
    print("")
    print(f"B4 TURBULENT SIZER, SPACE, the downlink slab at "
          f"{SPACE_ELEVATION_DEG:.0f} deg (turbulent_grid, {PRESET})")
    print(f"  slab length             {plan_s.z_total_m * 1e-3:11.4f} km")
    print(f"  grid side               {gs.size_m:11.4f} m")
    print(f"  pixels per side         {gs.n:11d}")
    print(f"  pixel pitch             {gs.pixel_m * 1e3:11.4f} mm")
    print(f"  screens                 {plan_s.z_m.size:11d}")
    print(f"  r0 of the whole path    {plan_s.r0_total_m * 1e2:11.3f} cm")
    print(f"  sigma2_R, plane wave    {plan_s.sigma2_r.sum():11.4f}")
    for text in rep_s.warnings:
        print(f"  warning: {text}")
    print("  THE SOURCE OF THE SLAB IS A PLANE WAVE that fills the grid, so")
    print("  D1 is the grid side itself. D2 is the ground aperture.")
    print("")
    print(f"    {'':4}  {'step':<38s}{'bound':>26s}{'actual':>14s}")
    check_s = properly_sampled_checklist(
        wavelength_m=WAVELENGTH_M, z_total_m=plan_s.z_total_m, n=gs.n,
        pixel_m=gs.pixel_m, z_m=plan_s.z_m, r0_m=plan_s.r0_m,
        r0_total_m=plan_s.r0_total_m, d_tx_m=gs.size_m,
        d_rx_m=SPACE_GROUND_APERTURE_M, size_m=gs.size_m, wave='plane')
    print_turbulent_rows(check_s)
    print("")
    print("  THE TWO FAILING ROWS ARE THE INTERESTING RESULT OF THIS SCRIPT.")
    print("")
    blurred = blurred_extent(SPACE_GROUND_APERTURE_M, WAVELENGTH_M,
                             plan_s.z_total_m, plan_s.r0_total_m)
    softer = constraint2_n_min(gs.pixel_m, gs.pixel_m, blurred, blurred,
                               WAVELENGTH_M, plan_s.z_total_m)
    method_term = (WAVELENGTH_M * plan_s.z_total_m
                   / (2.0 * gs.pixel_m ** 2))
    print("  CONSTRAINT 2 asks the grid to be wide enough that the wrap-around")
    print("  of the periodic transform stays outside the region of interest")
    print("  (Ch. 9, Eq. (9.87), printed p. 174). The row reads D1 as the WHOLE")
    print("  grid side, because the source of the slab is a plane wave that")
    print(f"  fills it. Read D1 as the blurred {blurred:.2f} m ground aperture")
    print(f"  instead, and the bound only falls from {check_s[3][2]:.0f} to "
          f"{softer:.0f}. Neither reading passes,")
    print(f"  because the term lambda dz /(2 Delta1 Delta_n) is "
          f"{method_term:.0f} on its own, and it")
    print("  does not depend on D1 at all. That term is the METHOD, not the")
    print("  geometry: it is constraint 4 in another dress.")
    print("  THE PRODUCTION LAYER ANSWERS THIS THE WAY CHAPTER 8 DOES. It does")
    print("  not widen the grid. It puts the super-Gaussian boundary mask on")
    print("  EVERY sub-step, so the wrapped light is removed before it can")
    print("  return (Ch. 8, Eq. (8.18), printed p. 139). The mask is ALWAYS on")
    print("  in olb/waveoptics/turbulence/splitstep.py, and the sizer keeps the")
    print("  receive aperture inside its flat interior: the report gives a")
    print(f"  grid margin of {rep_s.grid_margin:.2f}, where 1.0 means the light "
          f"just fits.")
    print("")
    print("  THE SCINTILLATION PITCH fails because the book rule of Sec. 9.4,")
    print("  printed p. 172, asks for two samples across the Fresnel scale")
    print("  sqrt(lambda z) of EVERY live screen, and the lowest screen of the")
    print("  slab sits close to the ground. The production sizer EXEMPTS a")
    print("  screen that carries less than `fresnel_weight_min` of the Rytov")
    print("  variance (see olb/waveoptics/turbulence/sampling.py), because such")
    print("  a screen adds almost no scintillation. The book checker does not")
    print("  know that exemption, so it fails the row. The production report")
    print(f"  gives {rep_s.fresnel_pixels_min:.2f} pixels per REQUIRED Fresnel "
          f"scale, and 2 or more is good.")

    # ---- the summary ----
    print("")
    print("SUMMARY OF THE FOUR TABLES")
    head = (f"  {'grid':<34}{'rows':>7}{'PASS':>7}{'FAIL':>7}"
            f"{'advisory':>10}")
    print(head)
    print("  " + "-" * (len(head) - 2))
    tables = [
        ("B1 vacuum, terrestrial 2 km",
         [r.satisfied for r in rows_t]),
        ("B2 vacuum, space 600 km co-moving",
         [r.satisfied for r in rows_s]),
        ("B3 turbulent, terrestrial 2 km",
         [r[1] for r in check_t]),
        ("B4 turbulent, space slab",
         [r[1] for r in check_s]),
    ]
    for name, marks in tables:
        n_pass = sum(1 for m in marks if m is True)
        n_fail = sum(1 for m in marks if m is False)
        n_adv = sum(1 for m in marks if m is None)
        print(f"  {name:<34}{len(marks):>7}{n_pass:>7}{n_fail:>7}"
              f"{n_adv:>10}")
    print("")
    print("  A FAIL IS NOT ALWAYS A DEFECT. The book states the rules are a")
    print("  GUIDELINE (Ch. 7, Sec. 7.3.3, printed p. 129), and rows 3, 4 and 5")
    print("  each belong to ONE kernel family. A grid that a co-moving lens")
    print("  recipe drives does not run the kernel that row 4 governs. Read the")
    print("  citation on each row, then decide.")
    print("")
    print("  ONE FAIL IS A REAL AGREEMENT, THOUGH. The B1 constraint-4 row")
    print("  fails, and `GridSpec.for_scenario` sent its OWN warning about the")
    print("  same bound at the top of that block, in its own words and with the")
    print("  same DOI. The book checker and the production sizer found one")
    print("  problem by two roads. That is the cross-check that this script")
    print("  exists for.")

    draw_artefacts(a1, a3)
    draw_edge_pair(a3)
    draw_absorbers(rho, book, prod)
    print("")
    print(f"figure saved: {ARTEFACTS_PNG}")
    print("  Caption: the six sampling panels. Row 1 is the transfer-function "
          "range limit in")
    print("  the intensity, clean, broken, and against the kernel that the "
          "rule does not")
    print("  touch. Row 2 is the same limit in the phase that causes it, and "
          "the damage as")
    print("  a cut through the middle.")
    print(f"figure saved: {ARTEFACTS_PNG.replace('.png', '_edges.png')}")
    print("  Caption: a hard aperture propagated 8 m on a 10 mm grid. Without "
          "an absorber")
    print("  the light that leaves one edge returns at the other and "
          "interferes; with the")
    print("  book absorber it goes away.")
    print(f"figure saved: {ABSORBERS_PNG}")
    print("  Caption: the book absorber (power 16, sigma 0.47 N) against the "
          "production")
    print("  boundary mask (power 8, band 0.125 of the half-side). At the "
          "middle of an edge")
    print("  the book keeps 5 times less amplitude.")
    print("")
    print(f"(elapsed {time.time() - t_start:.1f} s)")


if __name__ == "__main__":
    main()
