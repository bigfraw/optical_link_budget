'''
The five refractive-index spectra, and the effect of the two scales.

Andrews and Phillips, "Laser Beam Propagation through Random Media", 2nd ed.,
SPIE Press (2005), DOI: 10.1117/3.626196. Ch. 3, Eqs. (18) to (23), printed
pp. 67 to 69:
    Eq. (18)  Kolmogorov,  0.033 Cn2 kappa^(-11/3), valid 1/L0 << kappa << 1/l0
    Eq. (19)  Tatarskii,   adds the inner-scale cut exp(-kappa^2/kappa_m^2)
    Eq. (20)  von Karman,  adds the outer-scale knee (kappa^2 + kappa_0^2)
    Eq. (21)  exponential, an outer-scale form with 1 - exp(-kappa^2/kappa_0^2)
    Eq. (23)  modified atmospheric, the inner-scale BUMP plus an outer scale

The module shows three things:
  1. The five models at one Cn2, over four decades of the wavenumber kappa.
  2. The ratio of each model to Kolmogorov. This isolates the two scales:
     an outer scale cuts the LOW kappa end, an inner scale cuts the HIGH end,
     and the modified spectrum LIFTS the spectrum just before the inner-scale
     cut (the bump).
  3. The effect of the inner scale l0 and of the outer scale L0, one at a time.

The book prints more than one outer-scale constant on purpose (Ch. 3, text at
Eqs. (20), (22) and (23)). Each function takes that constant as a keyword.

Run from the repo root:
    python -m examples.andrews.spectra_and_scales
'''

import numpy as np

from olb.turbulence.andrews import (KOLMOGOROV_CONSTANT, MODIFIED_KL,
                                    TATARSKII_KM, exponential, kolmogorov,
                                    modified_atmospheric, tatarskii, von_karman)

# --- configuration ----------------------------------------------------------
CN2 = 1e-15          # constant Cn2 [m^-2/3]
L0_M = 10.0          # outer scale [m]
INNER_M = 5e-3       # inner scale [m]
KAPPA = np.logspace(-1.0, 3.0, 9)    # 0.1 to 1000 rad/m


def models(l0, L0):
    '''Return name -> Phi_n(kappa) for the five spectra at one scale pair.'''
    return {
        "kolmogorov": kolmogorov(KAPPA, CN2),
        "tatarskii": tatarskii(KAPPA, CN2, l0),
        "von_karman": von_karman(KAPPA, CN2, None, L0),
        "exponential": exponential(KAPPA, CN2, None, L0),
        "modified": modified_atmospheric(KAPPA, CN2, l0, L0),
    }


def print_absolute(table):
    '''Print Phi_n(kappa) [m^3] for each model.'''
    print(f"Phi_n(kappa) at Cn2={CN2:.1e} m^-2/3, l0={INNER_M*1e3:.0f} mm, "
          f"L0={L0_M:.0f} m")
    print(f"  {'kappa':>9} | " + " ".join(f"{n:>11}" for n in table))
    print("  " + "-" * (11 + 12 * len(table)))
    for i, k in enumerate(KAPPA):
        row = " ".join(f"{table[n][i]:>11.3e}" for n in table)
        print(f"  {k:>9.2f} | {row}")
    print()


def print_ratio(table):
    '''Print each model divided by Kolmogorov. This isolates the two scales.'''
    base = table["kolmogorov"]
    names = [n for n in table if n != "kolmogorov"]
    print("ratio to Kolmogorov (Eq. (18)). 1.000 = no scale effect")
    print(f"  {'kappa':>9} | " + " ".join(f"{n:>11}" for n in names))
    print("  " + "-" * (11 + 12 * len(names)))
    for i, k in enumerate(KAPPA):
        row = " ".join(f"{table[n][i] / base[i]:>11.4f}" for n in names)
        print(f"  {k:>9.2f} | {row}")
    print(f"  outer-scale knee kappa_0 = 2 pi / L0 = {2*np.pi/L0_M:.3f} rad/m")
    print(f"  inner-scale cut  kappa_m = {TATARSKII_KM:.2f} / l0 = "
          f"{TATARSKII_KM/INNER_M:.0f} rad/m")
    print(f"  modified bump    kappa_l = {MODIFIED_KL:.2f} / l0 = "
          f"{MODIFIED_KL/INNER_M:.0f} rad/m\n")


def sweep_inner(l0_values):
    '''Sweep l0 at a fixed L0. The inner scale acts at the HIGH kappa end.'''
    print("effect of the inner scale l0 (modified / Kolmogorov):")
    print(f"  {'kappa':>9} | " + " ".join(f"{v*1e3:>8.1f} mm" for v in l0_values))
    print("  " + "-" * (11 + 12 * len(l0_values)))
    base = kolmogorov(KAPPA, CN2)
    cols = [modified_atmospheric(KAPPA, CN2, v, L0_M) / base for v in l0_values]
    for i, k in enumerate(KAPPA):
        print(f"  {k:>9.2f} | " + " ".join(f"{c[i]:>11.4f}" for c in cols))
    print("  A larger l0 moves the bump and the cut to a LOWER kappa.\n")


def sweep_outer(L0_values):
    '''Sweep L0 at a fixed l0. The outer scale acts at the LOW kappa end.'''
    print("effect of the outer scale L0 (von Karman / Kolmogorov):")
    print(f"  {'kappa':>9} | " + " ".join(f"{v:>9.0f} m" for v in L0_values))
    print("  " + "-" * (11 + 12 * len(L0_values)))
    base = kolmogorov(KAPPA, CN2)
    cols = [von_karman(KAPPA, CN2, None, v) / base for v in L0_values]
    for i, k in enumerate(KAPPA):
        print(f"  {k:>9.2f} | " + " ".join(f"{c[i]:>11.4f}" for c in cols))
    print("  A larger L0 pushes the knee down, so the ratio goes to 1.\n")


if __name__ == '__main__':
    print(f"Andrews Ch. 3, Eqs. (18) to (23). Kolmogorov constant = "
          f"{KOLMOGOROV_CONSTANT}\n")
    table = models(INNER_M, L0_M)
    print_absolute(table)
    print_ratio(table)
    sweep_inner([1e-3, 5e-3, 20e-3])
    sweep_outer([1.0, 10.0, 100.0])
