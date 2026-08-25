'''
Gaussian-beam parameters at any curvature, and the three coherence radii.

Andrews and Phillips, "Laser Beam Propagation through Random Media", 2nd ed.,
SPIE Press (2005), DOI: 10.1117/3.626196:
    Ch. 4, Eq. (33), printed p. 92    input pair Theta0, Lambda0
    Ch. 4, Eqs. (37), (44), (45)      output pair Theta, Lambda, and W
    Ch. 7, Eq. (58), printed p. 242   effective (strong-turbulence) pair
    Ch. 8, Eq. (20), printed p. 268   Rytov variance sigma_R^2
    Ch. 6, Sec. 6.4; App. III, Tables I to VI, printed pp. 765 to 768
                                      coherence radius rho_0 by wave type
    Ch. 6, text below Eq. (64), printed p. 194     r_0 = 2.1 rho_0

Three tables:
  1. The beam parameters across the input curvature f0. A collimated beam has
     f0 = inf, so Theta0 = 1. A DIVERGED beam has f0 < 0, so Theta0 > 1. A
     FOCUSED beam has f0 > 0, so Theta0 < 1, and f0 = z gives Theta0 = 0.
  2. The effective pair against the Rytov variance. Strong turbulence spreads
     the beam, so the effective Lambda falls and the effective spot grows.
  3. The coherence radius and the Fried parameter for the plane wave, the
     spherical wave and the Gaussian beam.

The CURVATURE-GENERAL Fried parameter is table 4. It reads the beam through
andrews.beam.beam_params(w0, lambda, z, f0) and passes it to
andrews.structure.coherence_radius(..., wave="gaussian", beam=...). So a
deliberately diverged transmit beam drives its own r_0. This is the physics
half of the CLAUDE.md "Next task". See the README for what is still open at the
model layer.

Run from the repo root:
    python -m examples.andrews.beam_and_coherence
'''

import numpy as np

from olb.turbulence.andrews import (FRIED_OVER_RHO0, beam_params,
                                    coherence_radius, effective_beam_params,
                                    fried_parameter, rytov_variance)

# --- configuration ----------------------------------------------------------
WAVELENGTH_M = 1550e-9
PATH_M = 2000.0        # horizontal path [m]
CN2 = 3e-16            # constant Cn2 [m^-2/3]
WAIST_M = 0.05         # transmit Gaussian waist W0 [m]


def print_curvature(cases):
    '''Print the beam parameters for each input phase-front curvature f0.'''
    print(f"beam parameters at L={PATH_M:.0f} m, W0={WAIST_M*100:.0f} cm, "
          f"lambda={WAVELENGTH_M*1e9:.0f} nm (Ch. 4, Eqs. (33), (37), (44))")
    print(f"  {'case':>12} {'f0 [m]':>10} | {'Theta0':>8} {'Lambda0':>8} "
          f"{'Theta':>8} {'Lambda':>8} | {'W [cm]':>8} {'sig_R^2':>9}")
    print("  " + "-" * 82)
    for label, f0 in cases:
        bp = beam_params(WAIST_M, WAVELENGTH_M, PATH_M, f0=f0)
        try:
            s2 = f"{float(rytov_variance(WAVELENGTH_M, PATH_M, CN2, wave='gaussian', beam=bp)):.2e}"
        except NotImplementedError:
            # Ch. 8, Eq. (23), printed p. 264, holds for a collimated or a
            # divergent beam only. The package refuses; it does not guess.
            s2 = "refused"
        f0_text = "inf" if not np.isfinite(f0) else f"{f0:.0f}"
        print(f"  {label:>12} {f0_text:>10} | {bp.theta0:>8.3f} "
              f"{bp.lambda0:>8.3f} {bp.theta:>8.3f} {bp.lam:>8.3f} | "
              f"{bp.w*100:>8.2f} {s2:>9}")
    print("  Theta0 = 1 - L/f0. A diverged beam (f0 < 0) gives Theta0 > 1 and a "
          "wider W.")
    print("  'refused': the beam Rytov variance of Ch. 8, Eq. (23), printed "
          "p. 264, holds\n  for a collimated or a divergent beam only.\n")


def print_effective(sigma2_values):
    '''Print the effective pair of Ch. 7, Eq. (58) against the Rytov variance.'''
    bp = beam_params(WAIST_M, WAVELENGTH_M, PATH_M)
    print("effective (strong-turbulence) beam parameters, Ch. 7, Eq. (58), "
          "collimated beam")
    print(f"  {'sigma_R^2':>10} | {'Theta_eff':>10} {'Lambda_eff':>11} "
          f"{'W_eff [cm]':>11} {'W_eff/W':>9}")
    print("  " + "-" * 58)
    for s2 in sigma2_values:
        eff = effective_beam_params(bp, s2)
        print(f"  {s2:>10.2f} | {eff.theta:>10.3f} {eff.lam:>11.4f} "
              f"{eff.w*100:>11.2f} {eff.w/bp.w:>9.3f}")
    print("  Strong turbulence spreads the beam: W_eff grows, Lambda_eff "
          "falls.\n")


def print_coherence():
    '''Print rho_0 and r_0 = 2.1 rho_0 for the three wave types.'''
    bp = beam_params(WAIST_M, WAVELENGTH_M, PATH_M)
    print(f"coherence radius and Fried parameter, Cn2={CN2:.1e} m^-2/3, "
          f"r_0 = {FRIED_OVER_RHO0} rho_0")
    print(f"  {'wave':>12} | {'rho_0 [cm]':>11} {'r_0 [cm]':>10}")
    print("  " + "-" * 38)
    for wave in ("plane", "spherical", "gaussian"):
        beam = bp if wave == "gaussian" else None
        rho0 = float(coherence_radius(WAVELENGTH_M, PATH_M, CN2, wave=wave,
                                      beam=beam))
        print(f"  {wave:>12} | {rho0*100:>11.3f} "
              f"{float(fried_parameter(rho0))*100:>10.3f}")
    print("  A spherical wave keeps more coherence than a plane wave, because "
          "the turbulence near the source acts on a small beam.\n")


def print_curvature_fried(cases):
    '''Print the CURVATURE-GENERAL Gaussian-beam Fried parameter.'''
    print("curvature-general Gaussian-beam Fried parameter (Ch. 6 + Ch. 4)")
    print(f"  {'case':>12} {'f0 [m]':>10} | {'Theta0':>8} {'W [cm]':>8} "
          f"{'rho_0 [cm]':>11} {'r_0 [cm]':>10}")
    print("  " + "-" * 66)
    for label, f0 in cases:
        bp = beam_params(WAIST_M, WAVELENGTH_M, PATH_M, f0=f0)
        rho0 = float(coherence_radius(WAVELENGTH_M, PATH_M, CN2,
                                      wave='gaussian', beam=bp))
        f0_text = "inf" if not np.isfinite(f0) else f"{f0:.0f}"
        print(f"  {label:>12} {f0_text:>10} | {bp.theta0:>8.3f} "
              f"{bp.w*100:>8.2f} {rho0*100:>11.3f} "
              f"{float(fried_parameter(rho0))*100:>10.3f}")
    print("  The curvature changes r_0. So a diverged uplink beam MUST pass "
          "its f0.\n")


if __name__ == '__main__':
    cases = [("focused", PATH_M), ("half-focused", 2.0 * PATH_M),
             ("collimated", np.inf), ("diverged", -PATH_M),
             ("far diverged", -0.5 * PATH_M)]
    print_curvature(cases)
    print_effective([0.0, 0.25, 1.0, 4.0, 16.0])
    print_coherence()
    print_curvature_fried(cases)
