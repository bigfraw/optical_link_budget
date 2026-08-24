'''
Validation: multimode-fibre coupled power against the incident angle.

A multimode fibre is a light bucket. It collects ALL the focal-spot power that
lands inside the hard core disk. So the coupled power against a tip-tilt is the
ENCIRCLED ENERGY of the displaced Gaussian spot inside the core (a Marcum Q
function), NOT a Gaussian mode overlap.

This script plots the coupled power against the incident angle for the correct
model (olb.models.coupling.terrestrial._mmf_encircled_efficiency) and for the
OLD, WRONG model
(a Gaussian roll-off exp(-2*dx^2/a_core^2), the single-mode-fibre form with the
core radius). It shows two spot sizes. The correct model has a FLAT TOP: a small
spot deep inside the core loses almost nothing until it nears the edge. The old
model wrongly loses power from zero angle and drops far too fast.

Run from the repo root (saves a PNG; pass an output path to override):
    python -m examples.mmf_coupling_validation [out.png]
'''

import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from olb.models.coupling.terrestrial import _mmf_encircled_efficiency

# --- receiver configuration -------------------------------------------------
WAVELENGTH_M = 1550e-9
APERTURE_M = 0.2            # receive telescope diameter D
CORE_RADIUS_M = 25e-6      # multimode-fibre core radius a_core
ANGLE_URAD = np.linspace(0.0, 60.0, 400)   # incident angle sweep


def spot_radius(focal_m):
    '''Diffraction focal-spot 1/e^2 radius w_s = lambda*f/(pi*(D/2)).'''
    return WAVELENGTH_M * focal_m / (np.pi * APERTURE_M / 2.0)


def focal_for_spot_ratio(ratio):
    '''Focal length that sets the spot radius to ratio*a_core.'''
    return ratio * CORE_RADIUS_M * np.pi * (APERTURE_M / 2.0) / WAVELENGTH_M


def eta_correct(angle_rad, focal_m):
    '''Coupled power: encircled energy of the displaced spot inside the core.'''
    w_s = spot_radius(focal_m)
    offset = focal_m * np.asarray(angle_rad)          # dx = f*theta
    return _mmf_encircled_efficiency(offset, w_s, CORE_RADIUS_M)


def eta_old_wrong(angle_rad, focal_m):
    '''The OLD model: static overfill times a Gaussian roll-off over a_core.'''
    w_s = spot_radius(focal_m)
    offset = focal_m * np.asarray(angle_rad)
    eta_static = 1.0 - np.exp(-2.0 * CORE_RADIUS_M ** 2 / w_s ** 2)
    return eta_static * np.exp(-2.0 * offset ** 2 / CORE_RADIUS_M ** 2)


def main(out_path):
    angle_rad = ANGLE_URAD * 1e-6
    ratios = [0.2, 1.0]        # small spot; spot fills the core
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4), sharey=True)

    for ax, ratio in zip(axes, ratios):
        f = focal_for_spot_ratio(ratio)
        w_s = spot_radius(f)
        eta_c = eta_correct(angle_rad, f)
        eta_o = eta_old_wrong(angle_rad, f)
        # The incident angle that puts the spot centre at the core edge.
        edge_urad = CORE_RADIUS_M / f * 1e6

        ax.plot(ANGLE_URAD, eta_c, color="#2563eb", lw=2.2,
                label="correct (encircled energy)")
        ax.plot(ANGLE_URAD, eta_o, color="#dc2626", lw=2.0, ls="--",
                label="old model (Gaussian roll-off)")
        ax.axvline(edge_urad, color="#64748b", lw=1.0, ls=":",
                   label="spot centre at core edge")
        ax.set_title(f"w_s / a_core = {ratio:.1f}   "
                     f"(w_s={w_s * 1e6:.1f} um, f={f:.2f} m)", fontsize=10)
        ax.set_xlabel("incident angle  [urad]")
        ax.set_xlim(0, ANGLE_URAD[-1])
        ax.set_ylim(0, 1.02)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8, loc="upper right")

    axes[0].set_ylabel("coupled power fraction  eta")
    fig.suptitle(f"MMF light bucket: coupled power vs incident angle  "
                 f"(a_core={CORE_RADIUS_M * 1e6:.0f} um, D={APERTURE_M} m)",
                 fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(out_path, dpi=130)
    print(f"saved plot: {out_path}")

    # A numeric validation of the flat top: a small spot at half the edge angle.
    f = focal_for_spot_ratio(0.2)
    half_edge = 0.5 * CORE_RADIUS_M / f
    print(f"small spot (w_s=0.2 a_core): coupled power at half the edge angle = "
          f"{float(eta_correct(half_edge, f)):.4f} (correct, ~flat) vs "
          f"{float(eta_old_wrong(half_edge, f)):.4f} (old, wrongly low)")
    print(f"at the core edge: correct = {float(eta_correct(CORE_RADIUS_M / f, f)):.3f} "
          f"(~0.5, half in) vs old = {float(eta_old_wrong(CORE_RADIUS_M / f, f)):.3f} "
          f"(~0.135, far too low)")


if __name__ == '__main__':
    out = sys.argv[1] if len(sys.argv) > 1 else "mmf_coupling_vs_angle.png"
    main(out)
