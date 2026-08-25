'''
The three irradiance models at one matched index, and the fade statistics.

Andrews and Phillips, "Laser Beam Propagation through Random Media", 2nd ed.,
SPIE Press (2005), DOI: 10.1117/3.626196:
    Ch. 9, Eq. (138), printed p. 370  gamma-gamma alpha and beta
    Ch. 9, Eq. (139), printed p. 371  gamma-gamma index 1/a + 1/b + 1/(a b)
    Ch. 9, Sec. 9.9, printed p. 369   the K model, sigma_I^2 = 1 + 2/alpha > 1
    Ch. 11, Eqs. (23), (24), (25), printed p. 451   probability of fade
    Ch. 11, Eqs. (34), (37), printed pp. 455, 456   expected number of fades
    Ch. 11, Eq. (39), printed p. 456                mean fade time
    Ch. 12, Eqs. (69) to (74), printed pp. 511, 514 the satellite restatement

The three models get the SAME scintillation index, so the tables compare their
TAILS, not their strength. The index comes from the physics: a spherical wave at
sigma_R^2 = 10, through the two log variances of Ch. 9, Eqs. (41) and (46).

The K model needs sigma_I^2 > 1 (Ch. 9, printed p. 369). The strong case here
meets that. The script also shows the refusal at a weak index.

THE QUASI-FREQUENCY. The fade rate and the fade time need a rate scale nu0. The
book sets nu0 = 550 Hz for its own figures (printed pp. 457 and 514) instead of
computing it, because nu0 from the second spectral moment has NO upper limit of
its own with a Kolmogorov spectrum. This script prints BOTH: the computed
andrews.temporal.quasi_frequency over a stated band, and the book's 550 Hz. The
tables use the book value, so a reader can check them against the book.

THE IDENTITY. The book's own consistency check is Pr(fade) = <n> <t>. The last
table measures it.

Run from the repo root:
    python -m examples.andrews.distributions_and_fades
'''

import numpy as np

from olb.turbulence.andrews import (expected_number_of_fades,
                                    gamma_gamma_params,
                                    gamma_gamma_scintillation_index,
                                    irradiance_temporal_spectrum, k_params,
                                    k_scintillation_index,
                                    large_scale_log_variance, lognormal_params,
                                    mean_fade_time, probability_of_fade,
                                    quasi_frequency, small_scale_log_variance)
from olb.turbulence.andrews.distributions import MODELS

# --- configuration ----------------------------------------------------------
WAVELENGTH_M = 1550e-9
PATH_M = 2000.0        # horizontal path [m]
WIND_M_S = 10.0        # transverse wind, for the temporal spectrum
SIGMA2_R = 10.0        # Rytov variance of the strong case
BOOK_NU0_HZ = 550.0    # Andrews, printed pp. 457 and 514
FADE_DB = (3.0, 6.0, 10.0, 15.0, 20.0)


def strong_case():
    '''Return (sigma2_I, params by model) at one matched scintillation index.'''
    x = float(large_scale_log_variance(SIGMA2_R, wave='spherical'))
    y = float(small_scale_log_variance(SIGMA2_R, wave='spherical'))
    alpha, beta = gamma_gamma_params(x, y)
    sigma2_I = float(gamma_gamma_scintillation_index(alpha, beta))
    return sigma2_I, {
        "lognormal": {"sigma_l2": float(lognormal_params(sigma2_I))},
        "gamma_gamma": {"alpha": float(alpha), "beta": float(beta)},
        "k": {"alpha": float(k_params(sigma2_I))},
    }


def print_parameters(sigma2_I, params):
    '''Print the model parameters and prove that the index matches.'''
    print(f"matched case: spherical wave, sigma_R^2={SIGMA2_R:.1f} -> "
          f"sigma_I^2={sigma2_I:.4f}")
    for name, p in params.items():
        text = ", ".join(f"{k}={v:.4f}" for k, v in p.items())
        print(f"  {name:>12} : {text}")
    print(f"  K model check: 1 + 2/alpha = "
          f"{float(k_scintillation_index(params['k']['alpha'])):.4f}")
    try:
        k_params(0.5)
    except ValueError:
        print("  K model at sigma_I^2 = 0.5: REFUSED, the K index always "
              "exceeds 1.\n")


def print_quantiles(params, probabilities):
    '''Print the loss in dB at each irradiance quantile. This is the tail.'''
    print("fade depth [dB] at the irradiance quantile p (loss is positive dB)")
    print(f"  {'p':>9} | " + " ".join(f"{n:>12}" for n in params))
    print("  " + "-" * (11 + 13 * len(params)))
    for p in probabilities:
        row = []
        for name, kw in params.items():
            q = float(MODELS[name][1](p, **kw))
            row.append(f"{-10.0*np.log10(q):>12.3f}")
        print(f"  {p:>9.4f} | " + " ".join(row))
    print("  The gamma-gamma and K tails are DEEPER than the lognormal tail at "
          "the same index.\n")


def print_probability(params):
    '''Print Pr(fade) against the fade threshold F_T.'''
    print("probability of fade Pr(I <= I_T), Ch. 11, Eq. (23)")
    print(f"  {'F_T [dB]':>9} | " + " ".join(f"{n:>12}" for n in params))
    print("  " + "-" * (11 + 13 * len(params)))
    for f in FADE_DB:
        row = " ".join(f"{float(probability_of_fade(f, n, **kw)):>12.3e}"
                       for n, kw in params.items())
        print(f"  {f:>9.1f} | " + row)
    print()


def print_nu0():
    '''Print the computed quasi-frequency and the book reference value.'''
    freq = np.logspace(-1.0, 3.0, 2001)
    spectrum = irradiance_temporal_spectrum(freq, WIND_M_S, WAVELENGTH_M,
                                            PATH_M, 1e-15, wave='plane',
                                            regime='weak')
    nu0 = float(quasi_frequency(freq, spectrum))
    print(f"quasi-frequency, Ch. 12, Eq. (73): computed nu0 = {nu0:.1f} Hz "
          f"over 0.1 Hz to 1 kHz")
    print(f"  nu0 grows with the band, so the tables below use the book "
          f"reference {BOOK_NU0_HZ:.0f} Hz.\n")


def print_rate_and_time(params):
    '''Print the fade rate and the mean fade time, and check the identity.'''
    print(f"fade rate <n> [1/s] and mean fade time <t> [s] at nu0="
          f"{BOOK_NU0_HZ:.0f} Hz, Ch. 11, Eqs. (34), (39)")
    print(f"  {'F_T [dB]':>9} | " + " ".join(f"{n+' <n>':>15} {n+' <t>':>15}"
                                             for n in params))
    print("  " + "-" * (11 + 32 * len(params)))
    for f in FADE_DB:
        row = []
        for name, kw in params.items():
            n = float(expected_number_of_fades(f, BOOK_NU0_HZ, name, **kw))
            t = float(mean_fade_time(f, BOOK_NU0_HZ, name, **kw))
            row.append(f"{n:>15.4f} {t:>15.3e}")
        print(f"  {f:>9.1f} | " + " ".join(row))
    print()

    print("identity check Pr(fade) = <n> <t> (the book's own consistency test)")
    print(f"  {'F_T [dB]':>9} | " + " ".join(f"{n:>13}" for n in params))
    print("  " + "-" * (11 + 14 * len(params)))
    for f in FADE_DB:
        row = []
        for name, kw in params.items():
            pr = float(probability_of_fade(f, name, **kw))
            n = float(expected_number_of_fades(f, BOOK_NU0_HZ, name, **kw))
            t = float(mean_fade_time(f, BOOK_NU0_HZ, name, **kw))
            row.append(f"{abs(pr - n*t):>13.2e}")
        print(f"  {f:>9.1f} | " + " ".join(row))
    print("  The column is the absolute error. It is at machine precision.\n")


if __name__ == '__main__':
    sigma2_I, params = strong_case()
    print_parameters(sigma2_I, params)
    print_quantiles(params, (0.5, 0.1, 0.01, 1e-3, 1e-4))
    print_probability(params)
    print_nu0()
    print_rate_and_time(params)
