'''
The scintillation index from weak fluctuation through saturation.

Andrews and Phillips, "Laser Beam Propagation through Random Media", 2nd ed.,
SPIE Press (2005), DOI: 10.1117/3.626196:
    Ch. 8, Eq. (20), printed p. 268   Rytov variance sigma_R^2 (plane wave)
    Ch. 8, Eq. (23), printed p. 264   Gaussian-beam Rytov variance
    Ch. 9, Eq. (41), printed p. 335   large-scale log variance sigma_lnX^2
    Ch. 9, Eq. (46), printed p. 336   small-scale log variance sigma_lnY^2
    Ch. 9, Eqs. (48), (75), (104)     the weak two-scale index (l0 effect)
    Ch. 12, Eqs. (54), (57), (58)     tracked and untracked beam index

The all-regime index is sigma_I^2 = exp(sigma_lnX^2 + sigma_lnY^2) - 1. It is
valid at EVERY fluctuation strength, so the sweep below runs through the two
regime boundaries and out into saturation.

Two boundaries mark the table:
    sigma_R^2 = 0.25   the olb HOUSE limit. Below it the lognormal fade model
                       is trusted. See WEAK_FLUCTUATION_LIMIT.
    sigma_R^2 = 1.00   the BOOK weak-fluctuation limit,
                       andrews.scintillation.WEAK_REGIME_LIMIT.

Note the saturation: the plane-wave index rises above 1, peaks, and then falls
back towards 1 as sigma_R^2 grows. The Rytov variance itself keeps growing.

Run from the repo root:
    python -m examples.andrews.scintillation_regimes
'''

import numpy as np

from olb.turbulence.andrews import (WEAK_REGIME_LIMIT, beam_params,
                                    rytov_variance, scintillation_index,
                                    weak_two_scale_index)
from olb.turbulence.plane_wave_scintillation import WEAK_FLUCTUATION_LIMIT

# --- configuration ----------------------------------------------------------
WAVELENGTH_M = 1550e-9
PATH_M = 2000.0        # horizontal path [m]
WAIST_M = 0.05         # transmit Gaussian waist W0 [m]
OFF_AXIS_M = 0.05      # radial offset r for the tracked/untracked table [m]
WANDER_RMS_M = 0.004   # beam-wander rms displacement [m]
CN2_GRID = np.logspace(-17.0, -12.5, 10)


def beam():
    '''The collimated transmit beam at the receive plane.'''
    return beam_params(WAIST_M, WAVELENGTH_M, PATH_M)


def mark(sigma2_R):
    '''Name the regime of one Rytov variance.'''
    if sigma2_R < WEAK_FLUCTUATION_LIMIT:
        return "weak (house)"
    if sigma2_R < WEAK_REGIME_LIMIT:
        return "weak (book)"
    return "moderate/strong"


def sweep_waves():
    '''Print the index of the three wave types across the whole sweep.'''
    bp = beam()
    print(f"scintillation index against sigma_R^2, L={PATH_M:.0f} m, "
          f"lambda={WAVELENGTH_M*1e9:.0f} nm, W0={WAIST_M*100:.0f} cm")
    print(f"  {'Cn2':>9} {'sigma_R^2':>10} | {'plane':>9} {'spherical':>10} "
          f"{'gaussian':>9} | regime")
    print("  " + "-" * 72)
    for cn2 in CN2_GRID:
        s2r = float(rytov_variance(WAVELENGTH_M, PATH_M, cn2))
        p = float(scintillation_index(WAVELENGTH_M, PATH_M, cn2, wave='plane'))
        s = float(scintillation_index(WAVELENGTH_M, PATH_M, cn2,
                                      wave='spherical'))
        g = float(scintillation_index(WAVELENGTH_M, PATH_M, cn2,
                                      wave='gaussian', beam=bp))
        print(f"  {cn2:>9.1e} {s2r:>10.4f} | {p:>9.4f} {s:>10.4f} {g:>9.4f} | "
              f"{mark(s2r)}")
    print(f"  house limit {WEAK_FLUCTUATION_LIMIT}, book limit "
          f"{WEAK_REGIME_LIMIT}. The plane index peaks and then SATURATES "
          f"towards 1.\n")


def sweep_inner_scale(l0_values):
    '''Print the weak two-scale index. Only the INNER scale is implemented.'''
    bp = beam()
    cn2 = 1e-15
    s2r = float(rytov_variance(WAVELENGTH_M, PATH_M, cn2))
    print(f"inner-scale effect on the WEAK index, Cn2={cn2:.1e} "
          f"(sigma_R^2={s2r:.4f})")
    print(f"  {'l0 [mm]':>8} | {'plane':>9} {'spherical':>10} {'gaussian':>9}")
    print("  " + "-" * 42)
    for l0 in l0_values:
        row = []
        for wave in ("plane", "spherical", "gaussian"):
            b = bp if wave == "gaussian" else None
            row.append(float(weak_two_scale_index(WAVELENGTH_M, PATH_M, cn2,
                                                  wave=wave, l0=l0, beam=b)))
        print(f"  {l0*1e3:>8.1f} | {row[0]:>9.4f} {row[1]:>10.4f} "
              f"{row[2]:>9.4f}")
    zero = [float(scintillation_index(WAVELENGTH_M, PATH_M, cn2, wave=w,
                                      beam=bp if w == "gaussian" else None))
            for w in ("plane", "spherical", "gaussian")]
    print(f"  {'l0 -> 0':>8} | {zero[0]:>9.4f} {zero[1]:>10.4f} "
          f"{zero[2]:>9.4f}")
    print("  In the WEAK regime a finite inner scale RAISES the index above "
          "the l0 -> 0 row.")
    print("  The spectral bump of Ch. 3, Eq. (23) adds power near the "
          "Fresnel scale.")
    print("  The OUTER scale is refused on this weak form: the book prints no "
          "L0 branch\n  for Ch. 9, Eqs. (48), (75) and (104).\n")


def sweep_tracking():
    '''Print the tracked and the untracked Gaussian index off axis.'''
    bp = beam()
    print(f"tracked against untracked Gaussian beam, off axis r="
          f"{OFF_AXIS_M*100:.0f} cm, wander rms="
          f"{WANDER_RMS_M*1e3:.0f} mm (Ch. 12, Eqs. (54), (57))")
    print(f"  {'sigma_R^2':>10} | {'on axis':>9} {'tracked':>9} "
          f"{'untracked':>10} {'ratio':>7}")
    print("  " + "-" * 54)
    for cn2 in CN2_GRID[::3]:
        s2r = float(rytov_variance(WAVELENGTH_M, PATH_M, cn2))
        kw = dict(wave='gaussian', beam=bp, wander_rms_m=WANDER_RMS_M)
        axis = float(scintillation_index(WAVELENGTH_M, PATH_M, cn2, r=0.0,
                                         **kw))
        tr = float(scintillation_index(WAVELENGTH_M, PATH_M, cn2,
                                       r=OFF_AXIS_M, tracked=True, **kw))
        un = float(scintillation_index(WAVELENGTH_M, PATH_M, cn2,
                                       r=OFF_AXIS_M, tracked=False, **kw))
        print(f"  {s2r:>10.4f} | {axis:>9.4f} {tr:>9.4f} {un:>10.4f} "
              f"{un/tr:>7.3f}")
    print("  Tracking removes the beam-wander part of the RADIAL term only. On "
          "axis the\n  two are equal, because the radial term is zero there.\n")


if __name__ == '__main__':
    sweep_waves()
    sweep_inner_scale([1e-3, 3e-3, 10e-3])
    sweep_tracking()
