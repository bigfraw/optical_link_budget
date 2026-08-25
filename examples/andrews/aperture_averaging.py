'''
Aperture averaging: the book's soft-aperture chain against the Churnside fits.

Andrews and Phillips, "Laser Beam Propagation through Random Media", 2nd ed.,
SPIE Press (2005), DOI: 10.1117/3.626196, Ch. 10, Sec. 10.3:
    Ch. 10, Eq. (53), printed p. 409   spherical wave, weak
    Ch. 10, Eq. (60), printed p. 412   plane wave, weak, the EXACT form
    Ch. 10, Eq. (61), printed p. 412   plane wave, weak, the book's own fit
    Ch. 10, Eq. (69), printed p. 413   plane wave, all regimes
    Ch. 10, Eq. (77), printed p. 415   spherical wave, all regimes
    Ch. 10, Eqs. (87) to (90), printed p. 420   Gaussian beam, all regimes
The legacy route is J. H. Churnside, "Aperture averaging of optical
scintillations in the turbulent atmosphere", Appl. Opt. 30(15) 1982 (1991),
DOI: 10.1364/AO.30.001982. olb ships those fits in
olb/turbulence/plane_wave_scintillation.py and keeps them there, because the
constants are Churnside, not Andrews.

THE SOFT APERTURE. The book uses its own soft GAUSSIAN aperture, D_G^2 = 8 W_G^2
(Ch. 10, text below Eq. (57), printed p. 411). It is NOT the hard Airy filter
that olb integrates numerically. So the two routes use different aperture
shapes, and a difference between them is expected.

MEASURED (docs/andrews-crosscheck.md, WP3 notes): the Churnside weak fit is
OPTIMISTIC by 5 % to 13 % against the book's own exact Eq. (60) over d = 0.5 to
5. Optimistic means it predicts MORE averaging, so LESS fade. The last column of
the first table is that measurement.

The averaging factor A is dimensionless: A = sigma_I^2(D) / sigma_I^2(0). A
smaller A means more averaging.

Run from the repo root:
    python -m examples.andrews.aperture_averaging
'''

import numpy as np

from olb.turbulence.andrews import (averaging_factor, beam_params, d_param,
                                    omega_g, plane_weak_averaging_fit,
                                    rytov_variance)
from olb.turbulence.plane_wave_scintillation import (
    aperture_averaging_factor_strong, aperture_averaging_factor_weak)

# --- configuration ----------------------------------------------------------
WAVELENGTH_M = 1550e-9
PATH_M = 2000.0
CN2_WEAK = 2.9e-16      # gives sigma_R^2 about 0.021
CN2_STRONG = 1e-13      # gives sigma_R^2 about 7
GAUSS_WAIST_M = 0.01    # a wide beam at the receiver, so Omega_G >= Lambda
D_VALUES = (0.5, 1.0, 2.0, 5.0)


def diameters():
    '''Return (d, D) pairs. d = sqrt(k D^2 / 4L) is the book's own variable.'''
    k = 2.0 * np.pi / WAVELENGTH_M
    return [(d, d * np.sqrt(4.0 * PATH_M / k)) for d in D_VALUES]


def print_weak_plane():
    '''Print the exact form, the book fit and the Churnside fit, plane wave.'''
    s2r = float(rytov_variance(WAVELENGTH_M, PATH_M, CN2_WEAK))
    print(f"plane wave, WEAK regime, sigma_R^2={s2r:.4f}, "
          f"L={PATH_M/1e3:.0f} km, lambda={WAVELENGTH_M*1e9:.0f} nm")
    print(f"  {'d':>5} {'D [cm]':>8} | {'Eq. (60) exact':>15} "
          f"{'Eq. (61) fit':>13} {'Churnside':>10} | {'exact/Churn':>12}")
    print("  " + "-" * 72)
    for d, D in diameters():
        exact = float(averaging_factor(D, WAVELENGTH_M, PATH_M, CN2_WEAK,
                                       wave='plane', regime='weak'))
        fit = float(plane_weak_averaging_fit(D, WAVELENGTH_M, PATH_M))
        churn = float(aperture_averaging_factor_weak(D, WAVELENGTH_M, PATH_M))
        print(f"  {d:>5.1f} {D*100:>8.2f} | {exact:>15.5f} {fit:>13.5f} "
              f"{churn:>10.5f} | {exact/churn:>12.3f}")
    print("  The last column below 1 means the Churnside fit predicts MORE "
          "averaging than\n  the book, so a shallower fade. It is optimistic "
          "by 5 % to 13 % here.\n")


def print_weak_spherical():
    '''Print the spherical weak factor, Ch. 10, Eq. (53).'''
    print("spherical wave, WEAK regime, Ch. 10, Eq. (53)")
    print(f"  {'d':>5} {'D [cm]':>8} | {'A spherical':>12} {'A plane':>10} "
          f"{'ratio':>7}")
    print("  " + "-" * 48)
    for d, D in diameters():
        sph = float(averaging_factor(D, WAVELENGTH_M, PATH_M, CN2_WEAK,
                                     wave='spherical', regime='weak'))
        pln = float(averaging_factor(D, WAVELENGTH_M, PATH_M, CN2_WEAK,
                                     wave='plane', regime='weak'))
        print(f"  {d:>5.1f} {D*100:>8.2f} | {sph:>12.5f} {pln:>10.5f} "
              f"{sph/pln:>7.3f}")
    print("  A spherical wave averages LESS at the same D, because its "
          "irradiance cells\n  are larger near the receiver.\n")


def print_strong():
    '''Print the all-regime chain for the three wave types.'''
    s2r = float(rytov_variance(WAVELENGTH_M, PATH_M, CN2_STRONG))
    bp = beam_params(GAUSS_WAIST_M, WAVELENGTH_M, PATH_M)
    print(f"ALL-REGIME chain, sigma_R^2={s2r:.3f} (strong), Gaussian beam "
          f"W={bp.w*100:.1f} cm")
    print(f"  {'d':>5} {'D [cm]':>8} {'Omega_G':>8} | {'plane':>9} "
          f"{'spherical':>10} {'gaussian':>9} | {'Churn strong':>13}")
    print("  " + "-" * 76)
    for d, D in diameters():
        pln = float(averaging_factor(D, WAVELENGTH_M, PATH_M, CN2_STRONG,
                                     wave='plane', regime='strong'))
        sph = float(averaging_factor(D, WAVELENGTH_M, PATH_M, CN2_STRONG,
                                     wave='spherical', regime='strong'))
        gau = float(averaging_factor(D, WAVELENGTH_M, PATH_M, CN2_STRONG,
                                     wave='gaussian', regime='strong', beam=bp))
        churn = float(aperture_averaging_factor_strong(D, CN2_STRONG,
                                                       WAVELENGTH_M, PATH_M))
        print(f"  {d:>5.1f} {D*100:>8.2f} {float(omega_g(D, WAVELENGTH_M, PATH_M)):>8.2f}"
              f" | {pln:>9.5f} {sph:>10.5f} {gau:>9.5f} | {churn:>13.5f}")
    print("  The Gaussian chain of Eqs. (87) to (90) is an independent FIT. It "
          "does not\n  reduce exactly to the plane or the spherical form. It "
          "refuses Omega_G < Lambda.\n")


if __name__ == '__main__':
    d0, D0 = diameters()[0]
    print(f"d = sqrt(k D^2 / 4 L). Check: d({D0*100:.2f} cm) = "
          f"{float(d_param(D0, WAVELENGTH_M, PATH_M)):.3f}\n")
    print_weak_plane()
    print_weak_spherical()
    print_strong()
