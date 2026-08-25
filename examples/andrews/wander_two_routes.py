'''
Beam wander by two routes: the Andrews book and the Dios/Belmonte kernel.

Andrews and Phillips, "Laser Beam Propagation through Random Media", 2nd ed.,
SPIE Press (2005), DOI: 10.1117/3.626196:
    Ch. 6, Eq. (93), printed p. 203    <r_c^2>, constant 7.25
    Ch. 6, Eq. (94), printed p. 204    infinite outer scale, 2.42 Cn2 L^3 W0^-1/3
    Ch. 6, Eq. (100), printed p. 205   W_LT^2 = W_ST^2 + 1 <r_c^2>
    Ch. 8, Eq. (36), printed p. 271    pointing-error variance sigma_pe^2
    Ch. 12, Eqs. (50), (53), printed pp. 502, 503   the slant-path forms
The kernel route is J. A. Dios and others, "Scintillation and beam-wander
analysis in an optical ground station-satellite uplink", Appl. Opt. 43(19) 3866
(2004), DOI: 10.1364/AO.43.003866, Eq. (11), printed p. 3868, constant 2.07.

THE ADJUDICATED POSITION (Conflict C-01, docs/andrews-crosscheck.md). The two
routes share the SAME integrand and the SAME radial (two-axis) convention. Only
the leading constant differs: 7.25 against 2.07, a ratio of 3.5024 that this
script measures on both a terrestrial and an uplink case. The owner read the
Dios paper equation by equation against the kernel: the kernel is a FAITHFUL
copy of Dios Eq. (11), and Dios Fig. 3, printed p. 3871, compares that equation
with a split-step wave-optics simulation of the same uplink and the two agree
closely. A factor 3.50 would be easy to see on that figure. So olb KEEPS 2.07.
Andrews Eq. (93) is a different number, and neither source prints the filter
function that would explain the gap: Dios takes Eq. (11) from Belmonte, Appl.
Opt. 39, 5426 (2000), DOI: 10.1364/AO.39.005426, and that paper is the one to
read next to close the split. Until then the package holds BOTH routes, each
cited, and no code silently picks one.

A second, smaller number falls out. Dios Eq. (1) adds 2<beta^2> to W_ST^2 where
Andrews Eq. (100) adds 1<r_c^2>. So the wander part of the long-term radius
differs by 7.25/(2 x 2.07) = 1.7512, not by 3.5024.

Run from the repo root:
    python -m examples.andrews.wander_two_routes
'''

import numpy as np

from olb._deps import beam_wander_variance as kernel_wander
from olb._deps import long_term_beam_waist as kernel_long_term
from olb.beam import free_space_radius
from olb.turbulence.andrews import (WANDER_CONSTANT, beam_params,
                                    beam_wander_variance,
                                    beam_wander_variance_slant,
                                    long_term_beam_radius,
                                    pointing_error_variance,
                                    rytov_variance, short_term_beam_radius)
from olb.turbulence.profiles import DEFAULT_HS, get_c2n

# --- configuration ----------------------------------------------------------
WAVELENGTH_M = 1550e-9
DIOS_CONSTANT = 2.07             # Dios Eq. (11), printed p. 3868
# terrestrial case
TERR_PATH_M = 2000.0
TERR_CN2 = 3e-16
TERR_WAIST_M = 0.05
# uplink case (the defaults of the olb uplink_flux self-check)
UP_RANGE_M = 600e3
UP_WAIST_M = 1.0
UP_ELEVATION_DEG = 90.0


def terrestrial_routes():
    '''Return the two wander variances of the horizontal case [m^2].'''
    z = np.linspace(0.0, TERR_PATH_M, 4001)
    cn2 = np.full_like(z, TERR_CN2)
    w_free = free_space_radius(TERR_WAIST_M, z, None, WAVELENGTH_M)
    w_gom = np.full_like(z, TERR_WAIST_M)
    return {
        "Andrews Ch. 6, Eq. (93)": float(beam_wander_variance(
            TERR_WAIST_M, WAVELENGTH_M, TERR_PATH_M, TERR_CN2)),
        "kernel Dios Eq. (11), free-space W(z)": float(
            kernel_wander(TERR_PATH_M, cn2, w_free, z)),
        "kernel Dios Eq. (11), W(z) = W0": float(
            kernel_wander(TERR_PATH_M, cn2, w_gom, z)),
    }


def uplink_routes():
    '''Return the two wander variances of the slant uplink case [m^2].'''
    hs = DEFAULT_HS
    cn2 = get_c2n(hs, 21.0, 1.7e-14)
    w_free = free_space_radius(UP_WAIST_M, hs, None, WAVELENGTH_M)
    return {
        "Andrews Ch. 12, Eq. (50) slant": float(beam_wander_variance_slant(
            UP_WAIST_M, WAVELENGTH_M, hs, cn2, UP_RANGE_M,
            elevation_deg=UP_ELEVATION_DEG)),
        "kernel Dios Eq. (11), free-space W(z)": float(
            kernel_wander(UP_RANGE_M, cn2, w_free, hs)),
    }


def print_case(title, routes, waist_m, range_m):
    '''Print both routes of one case, with the rms and the ratio.'''
    print(title)
    print(f"  {'route':>40} | {'<r_c^2> [m^2]':>14} {'rms [m]':>11} "
          f"{'rms [urad]':>11}")
    print("  " + "-" * 82)
    values = []
    for name, value in routes.items():
        values.append(value)
        rms = np.sqrt(value)
        print(f"  {name:>40} | {value:>14.6e} {rms:>11.4e} "
              f"{rms/range_m*1e6:>11.4f}")
    print(f"  ratio Andrews / kernel (free-space W(z)) : "
          f"{values[0]/values[1]:.4f}")
    if len(values) > 2:
        print(f"  ratio Andrews / kernel (same W(z) = W0)  : "
              f"{values[0]/values[2]:.4f}")
    print(f"  the book constants                       : "
          f"{WANDER_CONSTANT} / {DIOS_CONSTANT} = "
          f"{WANDER_CONSTANT/DIOS_CONSTANT:.4f}")
    print(f"  transmit waist W0 = {waist_m:.2f} m, range = {range_m/1e3:.0f} km\n")


def print_long_term():
    '''Compare the two long-term radius rules on the same short-term radius.'''
    bp = beam_params(TERR_WAIST_M, WAVELENGTH_M, TERR_PATH_M)
    sigma2_R = float(rytov_variance(WAVELENGTH_M, TERR_PATH_M, TERR_CN2,
                                    wave='gaussian', beam=bp))
    rc2 = float(beam_wander_variance(TERR_WAIST_M, WAVELENGTH_M, TERR_PATH_M,
                                     TERR_CN2))
    w_lt = float(long_term_beam_radius(bp, sigma2_R))
    w_st = float(short_term_beam_radius(bp, sigma2_R, rc2))
    z = np.linspace(0.0, TERR_PATH_M, 4001)
    beta2 = float(kernel_wander(TERR_PATH_M, np.full_like(z, TERR_CN2),
                                free_space_radius(TERR_WAIST_M, z, None,
                                                  WAVELENGTH_M), z))
    kernel_lt = float(kernel_long_term(w_st, beta2))
    print("long-term beam radius on the terrestrial case (Conflict C-03)")
    print(f"  short-term radius W_ST, Andrews Ch. 6, Eq. (100) : "
          f"{w_st*100:.4f} cm")
    print(f"  long-term radius W_LT, Andrews Ch. 6, Eq. (98)   : "
          f"{w_lt*100:.4f} cm")
    print(f"  long-term radius W_LT, Dios Eq. (1) on W_ST      : "
          f"{kernel_lt*100:.4f} cm")
    print(f"  residual factor 7.25 / (2 x 2.07)               : "
          f"{WANDER_CONSTANT/(2.0*DIOS_CONSTANT):.4f}")
    print(f"  pointing-error variance, Ch. 8, Eq. (36)        : "
          f"{float(pointing_error_variance(TERR_WAIST_M, WAVELENGTH_M, TERR_PATH_M, TERR_CN2)):.4e} m^2\n")


if __name__ == '__main__':
    print_case(f"TERRESTRIAL: L={TERR_PATH_M/1e3:.0f} km, "
               f"Cn2={TERR_CN2:.1e} m^-2/3, collimated",
               terrestrial_routes(), TERR_WAIST_M, TERR_PATH_M)
    print_case(f"UPLINK: slant range {UP_RANGE_M/1e3:.0f} km, zenith, "
               f"H-V 5/7 profile, collimated",
               uplink_routes(), UP_WAIST_M, UP_RANGE_M)
    print_long_term()
    print("The ratio is the SAME on both cases and on both beam-radius "
          "readings, so the\nwhole gap is the leading constant. See the module "
          "docstring for the position.")
