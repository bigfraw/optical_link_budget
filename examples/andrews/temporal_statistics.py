'''
Temporal statistics: the irradiance spectrum, nu0, the Greenwood frequency.

Andrews and Phillips, "Laser Beam Propagation through Random Media", 2nd ed.,
SPIE Press (2005), DOI: 10.1117/3.626196:
    Ch. 3, Eq. (27), printed p. 73          Taylor wavenumber kappa = 2 pi f / V
    Ch. 8, text below Eq. (57), printed p. 283   Fresnel frequency
    Ch. 8, Eq. (65), printed p. 285         weak irradiance spectrum
    Ch. 9, Eqs. (126) to (128), printed p. 365   strong spectrum, point receiver
    Ch. 10, Eqs. (93) to (97), printed pp. 421, 422  strong spectrum, aperture
    Ch. 11, Eqs. (14), (15), printed p. 448; Ch. 12, Eq. (73), printed p. 514
                                            quasi-frequency nu0
    Ch. 14, Eqs. (38), (39), printed pp. 622, 623   Greenwood frequency, tau0

THREE WARNINGS THAT THE PACKAGE CARRIES, shown by the tables below.
  1. nu0 has NO upper limit of its own. With a Kolmogorov spectrum and a zero
     inner scale the spectrum decays as f^(-8/3), so the second moment grows as
     f_max^(1/3): about x1.49 per decade of band. This is why the book fixes
     nu0 = 550 Hz for its own figures. A caller MUST set the band.
  2. The weak spectrum and the strong spectrum carry the same POWER but not the
     same SHAPE. The ratio table shows that.
  3. An inner scale or an outer scale is REFUSED on the temporal spectrum, in
     every regime (Ch. 9, Sec. 9.8, printed p. 364). So is a strong-regime
     spherical wave or Gaussian beam.

Run from the repo root:
    python -m examples.andrews.temporal_statistics
'''

import numpy as np

from olb.turbulence.andrews import (GREENWOOD_CONSTANT, coherence_time,
                                    expected_number_of_fades,
                                    fresnel_frequency, greenwood_frequency,
                                    irradiance_temporal_spectrum,
                                    lognormal_params, mean_fade_time,
                                    quasi_frequency, rytov_variance,
                                    taylor_wavenumber)
from olb.turbulence.profiles import DEFAULT_HS, get_c2n

# --- configuration ----------------------------------------------------------
WAVELENGTH_M = 1550e-9
PATH_M = 2000.0
CN2 = 1e-15
WIND_M_S = 10.0
APERTURE_M = 0.10
FADE_DB = (3.0, 6.0, 10.0)


def print_spectrum(freqs):
    '''Print the weak and the strong spectrum, point and aperture receiver.'''
    ff = float(fresnel_frequency(WIND_M_S, WAVELENGTH_M, PATH_M))
    print(f"irradiance temporal spectrum, V={WIND_M_S:.0f} m/s, "
          f"L={PATH_M/1e3:.0f} km, Cn2={CN2:.1e}")
    print(f"  Fresnel frequency f_F = {ff:.2f} Hz, sigma_R^2 = "
          f"{float(rytov_variance(WAVELENGTH_M, PATH_M, CN2)):.4f}")
    print(f"  {'f [Hz]':>9} {'kappa':>9} | {'weak':>11} {'strong':>11} "
          f"{'strong+D':>11} | {'strong/weak':>12}")
    print("  " + "-" * 72)
    kw = dict(wave='plane')
    weak = irradiance_temporal_spectrum(freqs, WIND_M_S, WAVELENGTH_M, PATH_M,
                                        CN2, regime='weak', **kw)
    strong = irradiance_temporal_spectrum(freqs, WIND_M_S, WAVELENGTH_M,
                                          PATH_M, CN2, regime='strong', **kw)
    with_d = irradiance_temporal_spectrum(freqs, WIND_M_S, WAVELENGTH_M,
                                          PATH_M, CN2, regime='strong',
                                          D=APERTURE_M, **kw)
    for i, f in enumerate(freqs):
        kappa = float(taylor_wavenumber(f, WIND_M_S))
        print(f"  {f:>9.1f} {kappa:>9.3f} | {weak[i]:>11.4e} "
              f"{strong[i]:>11.4e} {with_d[i]:>11.4e} | "
              f"{strong[i]/weak[i]:>12.4f}")
    print("  The last column is not flat, so the two spectra differ in SHAPE.")
    print(f"  A {APERTURE_M*100:.0f} cm aperture cuts the high-frequency power "
          f"(the aperture averages it).\n")


def print_nu0_band(bands):
    '''Print nu0 against the upper band edge. It has no limit of its own.'''
    print("quasi-frequency nu0 against the band, weak plane spectrum")
    print(f"  {'f_max [Hz]':>11} | {'nu0 [Hz]':>10} {'growth':>8}")
    print("  " + "-" * 34)
    previous = None
    for f_max in bands:
        freq = np.logspace(-1.0, np.log10(f_max), 4001)
        spectrum = irradiance_temporal_spectrum(freq, WIND_M_S, WAVELENGTH_M,
                                                PATH_M, CN2, wave='plane',
                                                regime='weak')
        nu0 = float(quasi_frequency(freq, spectrum))
        growth = "-" if previous is None else f"{nu0/previous:.3f}"
        print(f"  {f_max:>11.0f} | {nu0:>10.2f} {growth:>8}")
        previous = nu0
    print("  The growth never stops. It settles at about 1.49 per decade, the "
          "f_max^(1/3)\n  law. Set the band from the detector bandwidth.\n")


def print_greenwood(elevations):
    '''Print the Greenwood frequency and tau0 on the H-V profile.'''
    hs = DEFAULT_HS
    cn2 = get_c2n(hs, 21.0, 1.7e-14)
    print(f"Greenwood frequency and coherence time, H-V 5/7 profile, "
          f"constant {GREENWOOD_CONSTANT}")
    print(f"  {'elev [deg]':>11} | {'f_G [Hz]':>10} {'tau0 [ms]':>11}")
    print("  " + "-" * 37)
    for e in elevations:
        fg = float(greenwood_frequency(hs, cn2, WAVELENGTH_M, e))
        print(f"  {e:>11.0f} | {fg:>10.2f} "
              f"{float(coherence_time(fg))*1e3:>11.4f}")
    print("  A low elevation gives a longer path, so a faster atmosphere and a "
          "shorter tau0.\n")


def print_fade_numbers():
    '''Feed the computed nu0 into the fade rate and the mean fade time.'''
    freq = np.logspace(-1.0, 3.0, 4001)
    spectrum = irradiance_temporal_spectrum(freq, WIND_M_S, WAVELENGTH_M,
                                            PATH_M, CN2, wave='plane',
                                            regime='weak')
    nu0 = float(quasi_frequency(freq, spectrum))
    sigma2_I = float(rytov_variance(WAVELENGTH_M, PATH_M, CN2))
    sl2 = float(lognormal_params(sigma2_I))
    print(f"fade numbers from THIS nu0 = {nu0:.1f} Hz (band 0.1 Hz to 1 kHz), "
          f"lognormal sigma_I^2={sigma2_I:.4f}")
    print(f"  {'F_T [dB]':>9} | {'<n> [1/s]':>11} {'<t> [ms]':>11} "
          f"{'<n><t>':>10}")
    print("  " + "-" * 48)
    for f in FADE_DB:
        n = float(expected_number_of_fades(f, nu0, "lognormal", sigma_l2=sl2))
        t = float(mean_fade_time(f, nu0, "lognormal", sigma_l2=sl2))
        print(f"  {f:>9.1f} | {n:>11.4e} {t*1e3:>11.4f} {n*t:>10.3e}")
    print("  <n> <t> is the probability of fade. Both scale with nu0, so a "
          "wider band\n  raises the rate and shortens the fade in the same "
          "proportion.\n")


if __name__ == '__main__':
    print_spectrum(np.array([1.0, 10.0, 50.0, 100.0, 500.0, 1000.0]))
    print_nu0_band([10.0, 100.0, 1000.0, 10000.0, 100000.0])
    print_greenwood([10.0, 20.0, 30.0, 45.0, 60.0, 90.0])
    print_fade_numbers()
