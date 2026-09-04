'''
Irradiance distribution models on the normalised irradiance I, with E[I] = 1.

Every model in this module speaks the same language. The random quantity is the
NORMALISED irradiance I. Its mean is one. A model gives some of these faces:

    <model>_params(...)     turn the turbulence variances into the parameters
    <model>_mean_log(...)   E[ln I]; this is the mean loss in a log measure
    <model>_pdf(I, ...)     the probability density
    <model>_cdf(I, ...)     Pr(I <= value); this is the fade probability
    <model>_quantile(p,...) the value of I that Pr(I <= value) = p
    <model>_rvs(n, ..., rng) n Monte Carlo draws of I

This module holds NO decibels. It gives pure irradiance physics. The model layer
turns a face into a dB loss. See olb/models/fade.py.

Source of every equation:
    L. C. Andrews and R. L. Phillips, "Laser Beam Propagation through Random
    Media", 2nd ed., SPIE Press (2005). DOI: 10.1117/3.626196
Each function names its section, its equation number, and its printed page.

The models:

- Lognormal (Ch. 5.7.2, printed 156-157). The first-order Rytov result. It is a
  weak-fluctuation model. The book states at printed 451 that the lognormal
  gives optimistic fade probabilities in the deep tail.
- Gamma-gamma (Ch. 9.10, printed 369-371). The book's preferred model. It is
  valid in ALL fluctuation regimes, and its two parameters come straight from
  the large-scale and small-scale log-irradiance variances. No fit is necessary.
- K (Ch. 9.9.1, printed 368-369). The gamma-gamma model at beta = 1. Its
  scintillation index always exceeds one, so it is a strong-fluctuation model
  only.
- Lognormal-Rician (Ch. 9.9.2, printed 369). PDF only, for comparison. See the
  warning in the function docstring.

NOT built, and the reason: the modified Rician PDF (Ch. 5.7.1, Eq. (87),
printed 155). The book compares its normalised moments with measured data at
printed 155 and reports that the theory lies below the data. It states that the
modified Rician is not a suitable model for irradiance fluctuations. So this
module does not give it.
'''

from functools import lru_cache

import numpy as np
from scipy.integrate import cumulative_trapezoid, quad
from scipy.special import digamma, gammaln, i0e, kve
from scipy.stats import norm

from ...assumptions import REGIME_STRONG, REGIME_WEAK, Constraint, assumes
from .scintillation import LOGNORMAL_PDF_LIMIT

# Ch. 11, Eq. (25), printed 451, and Ch. 12, Eq. (69), printed 511, define the
# fade threshold parameter F_T = 10 log10(E[I] / I_T) in dB. The book restates
# it at Ch. 12, Eq. (71), printed 511, as ln(I_T / E[I]) = -0.23 F_T. The
# constant 0.23 is ln(10) / 10 to two figures. The code uses the exact value.
_LN10_OVER_10 = np.log(10.0) / 10.0

# The gamma-gamma CDF grid. The module integrates the PDF on a log-spaced
# irradiance grid, then interpolates. The range covers the full mass for every
# parameter pair that a link budget uses.
_GRID_LO_LOG10 = -12.0
_GRID_HI_LOG10 = 2.0
_GRID_N = 4001


# --- Function-owned assumptions ---------------------------------------------
# Each distribution states its own regime. The lognormal is a weak-fluctuation
# model whose fade PDF olb trusts only up to the scintillation index
# LOGNORMAL_PDF_LIMIT. The gamma-gamma is valid at every fluctuation strength.
# The K distribution is strong only. The lognormal-Rician gives the PDF alone.
_DOI = "10.1117/3.626196"


def _lognormal_pdf_shape_check(args, result):
    '''Return a reason when a lognormal face is used past the PDF-shape limit.

    The lognormal tail goes optimistic against simulation well before the Rytov
    theory for the index fails, so olb trusts the lognormal fade PDF only up to
    the scintillation index LOGNORMAL_PDF_LIMIT. The check reads the variance
    from the bound arguments: `lognormal_params` takes the scintillation index
    sigma_I^2 directly; the other faces take the log-irradiance variance
    sigma_l2, and sigma_I^2 = exp(sigma_l2) - 1. No warning here.
    '''
    if 'sigma2_I' in args:
        sigma2_I = float(np.max(np.asarray(args['sigma2_I'], dtype=float)))
    else:
        sigma_l2 = float(np.max(np.asarray(args['sigma_l2'], dtype=float)))
        sigma2_I = float(np.expm1(sigma_l2))
    if sigma2_I > LOGNORMAL_PDF_LIMIT:
        return (f"sigma_I^2 = {sigma2_I:.3f} > {LOGNORMAL_PDF_LIMIT}; the "
                f"lognormal fade PDF goes optimistic in the deep tail. Use the "
                f"gamma-gamma model.")
    return None


LOGNORMAL_PDF_CONSTRAINT = Constraint(
    "pdf-shape",
    "The lognormal fade PDF is trusted only up to the scintillation index "
    "sigma_I^2 = 0.25. Past it the deep tail is optimistic against simulation.",
    _DOI, "Ch. 11, Sec. 11.3, printed p. 451",
    check=_lognormal_pdf_shape_check)

GAMMA_GAMMA_ALL_STRENGTH = Constraint(
    "regime",
    "The gamma-gamma model is valid at every fluctuation strength, from weak to "
    "strong. Its two parameters come straight from the large-scale and "
    "small-scale log-irradiance variances, with no fit.",
    _DOI, "Ch. 9, Sec. 9.10, printed pp. 369 to 371")

K_STRONG_ONLY = Constraint(
    "regime",
    "The K distribution applies in strong fluctuations only: its scintillation "
    "index sigma_I^2 = 1 + 2/alpha always exceeds one.",
    _DOI, "Ch. 9, Sec. 9.9.1, printed pp. 368 to 369")

LOGNORMAL_RICIAN_NOT_BUILT = Constraint(
    "not-built",
    "The lognormal-Rician model gives the PDF only: no CDF, no quantile and no "
    "sampler. The book states its coherence parameter and modulation variance "
    "cannot be tied to atmospheric conditions, so it is not a link-budget "
    "model. Use the gamma-gamma model for a budget.",
    _DOI, "Ch. 9, Sec. 9.9.2, Eq. (133), printed p. 369")


# --- lognormal (Ch. 5.7.2, printed 156-157) ---------------------------------

@assumes(LOGNORMAL_PDF_CONSTRAINT, turbulence_regime=REGIME_WEAK)
def lognormal_params(sigma2_I):
    '''
    Turn a scintillation index into the log-irradiance variance.

    Ch. 5, Eq. (95), printed 157, gives sigma_I^2 = exp(4 sigma_x^2) - 1, where
    sigma_x^2 is the log-amplitude variance. The log-irradiance variance is
    sigma_lnI^2 = 4 sigma_x^2, so the inverse map is
        sigma_lnI^2 = ln(1 + sigma_I^2).
    DOI: 10.1117/3.626196

    Parameters:
        sigma2_I : float or ndarray
            Scintillation index (the normalised variance of I).

    Returns:
        float or ndarray
            sigma_l2, the variance of ln I.
    '''
    return np.log(1.0 + np.asarray(sigma2_I, dtype=float))


@assumes(LOGNORMAL_PDF_CONSTRAINT, turbulence_regime=REGIME_WEAK)
def lognormal_mean_log(sigma_l2):
    '''
    E[ln I] for a lognormal irradiance with E[I] = 1.

    Ch. 5, Eq. (93), printed 156, gives the lognormal PDF. Ch. 11, Eq. (22),
    printed 451, writes the same PDF in the E[I] = 1 form. The mean of the log
    is then -sigma_l2 / 2, because E[I] = exp(E[ln I] + sigma_l2 / 2) = 1.
    DOI: 10.1117/3.626196

    Parameters:
        sigma_l2 : float or ndarray
            Variance of ln I.

    Returns:
        float or ndarray
            E[ln I].
    '''
    return -np.asarray(sigma_l2, dtype=float) / 2.0


@assumes(LOGNORMAL_PDF_CONSTRAINT, turbulence_regime=REGIME_WEAK)
def lognormal_cdf(I, sigma_l2):
    '''
    Pr(I <= value) for a lognormal irradiance with E[I] = 1.

    Ch. 11, Eq. (24), printed 451, gives the same result in the erf form
        Pr = 0.5 {1 + erf[(0.5 sigma_I^2 - 0.23 F_T) / (sigma_I sqrt(2))]},
    with the fade threshold parameter F_T of Ch. 11, Eq. (25). This code uses
    the standard normal CDF of (ln I + sigma_l2/2)/sigma_l, which is the same
    expression. The book writes sigma_I in place of sigma_l because the two
    agree in the weak-fluctuation limit. This code keeps the exact
    sigma_l = sqrt(ln(1 + sigma_I^2)).
    DOI: 10.1117/3.626196

    CAUTION on the sign of the 0.23 F_T term. The scanned text of Eq. (24)
    interleaves two columns, so a reader cannot read the sign off the page.
    Ch. 12, Eq. (71), printed 511, states ln(I_T / E[I]) = -0.23 F_T, and the
    lognormal CDF at that threshold is Phi[(ln I_T + sigma_l2/2) / sigma_l].
    That gives the MINUS sign above. A check confirms it: at F_T = 0 the
    threshold equals the mean, and the result must exceed one half, because the
    median of a unit-mean lognormal lies below one. The minus sign gives
    Phi(sigma_l/2) > 0.5. The plus sign does not.

    Parameters:
        I : float or ndarray
            Normalised irradiance, I > 0.
        sigma_l2 : float
            Variance of ln I.

    Returns:
        float or ndarray
            The cumulative probability.
    '''
    sigma_l = np.sqrt(sigma_l2)
    return norm.cdf((np.log(I) - lognormal_mean_log(sigma_l2)) / sigma_l)


@assumes(LOGNORMAL_PDF_CONSTRAINT, turbulence_regime=REGIME_WEAK)
def lognormal_quantile(p, sigma_l2):
    '''
    The p-quantile of a lognormal irradiance with E[I] = 1.

    Inverse of Ch. 11, Eq. (24), printed 451:
        I(p) = exp(-sigma_l2/2 + sigma_l Phi^-1(p)).
    DOI: 10.1117/3.626196

    Parameters:
        p : float or ndarray
            Cumulative probability in (0, 1).
        sigma_l2 : float or ndarray
            Variance of ln I.

    Returns:
        float or ndarray
            The normalised irradiance at that probability.
    '''
    sigma_l2 = np.asarray(sigma_l2, dtype=float)
    sigma_l = np.sqrt(sigma_l2)
    return np.exp(-sigma_l2 / 2.0 + sigma_l * norm.ppf(p))


@assumes(LOGNORMAL_PDF_CONSTRAINT, turbulence_regime=REGIME_WEAK)
def lognormal_rvs(n, sigma_l2, rng):
    '''
    Draw n samples of a lognormal irradiance with E[I] = 1.

    Ch. 5, Eq. (93), printed 156. DOI: 10.1117/3.626196

    Parameters:
        n : int
            Number of draws.
        sigma_l2 : float or ndarray
            Variance of ln I. An array gives one distribution per element, and
            the result has shape (n, *shape(sigma_l2)).
        rng : numpy.random.Generator
            The generator.

    Returns:
        ndarray
            Shape (n, *shape(sigma_l2)).
    '''
    sigma_l2 = np.asarray(sigma_l2, dtype=float)
    return rng.lognormal(mean=-sigma_l2 / 2.0, sigma=np.sqrt(sigma_l2),
                         size=(n, *np.shape(sigma_l2)))


# --- gamma-gamma (Ch. 9.10, printed 369-371) --------------------------------

@assumes(GAMMA_GAMMA_ALL_STRENGTH)
def gamma_gamma_params(sigma2_lnX, sigma2_lnY):
    '''
    The two gamma-gamma parameters from the log-irradiance variances.

    Ch. 9, Eq. (138), printed 370:
        alpha = 1 / [exp(sigma_lnX^2) - 1],
        beta  = 1 / [exp(sigma_lnY^2) - 1].
    The book restates this as Ch. 9, Eq. (160), printed 384, as Ch. 11, Eq.
    (20), printed 450, and as Ch. 12, Eq. (68), printed 511.
    DOI: 10.1117/3.626196

    The large-scale variance sigma_lnX^2 and the small-scale variance
    sigma_lnY^2 come from the propagation model. For a plane wave the book gives
    them at Ch. 9, Eqs. (41) and (46), printed 335-336.

    Parameters:
        sigma2_lnX : float
            Large-scale log-irradiance variance.
        sigma2_lnY : float
            Small-scale log-irradiance variance.

    Returns:
        (float, float)
            (alpha, beta).
    '''
    alpha = 1.0 / (np.exp(sigma2_lnX) - 1.0)
    beta = 1.0 / (np.exp(sigma2_lnY) - 1.0)
    return alpha, beta


@assumes(GAMMA_GAMMA_ALL_STRENGTH)
def gamma_gamma_scintillation_index(alpha, beta):
    '''
    The scintillation index of a gamma-gamma irradiance.

    Ch. 9, Eq. (139), printed 371:
        sigma_I^2 = 1/alpha + 1/beta + 1/(alpha beta).
    DOI: 10.1117/3.626196

    Parameters:
        alpha : float
            Large-scale parameter.
        beta : float
            Small-scale parameter.

    Returns:
        float
            The scintillation index (the variance of I, because E[I] = 1).
    '''
    return 1.0 / alpha + 1.0 / beta + 1.0 / (alpha * beta)


@assumes(GAMMA_GAMMA_ALL_STRENGTH)
def gamma_gamma_mean_log(alpha, beta):
    '''
    E[ln I] for a gamma-gamma irradiance with E[I] = 1.

    Ch. 9, Eqs. (134) and (135), printed 370, make I = X Y the product of two
    gamma variates. X has shape alpha and scale 1/alpha, and Y has shape beta
    and scale 1/beta. The log-mean of a gamma variate with shape k and scale t
    is the standard result E[ln X] = psi(k) + ln t, where psi is the digamma
    function (see NIST DLMF 5.2.2 for psi). So
        E[ln I] = psi(alpha) + psi(beta) - ln(alpha beta).
    DOI: 10.1117/3.626196

    Parameters:
        alpha : float
            Large-scale parameter.
        beta : float
            Small-scale parameter.

    Returns:
        float
            E[ln I].
    '''
    return digamma(alpha) + digamma(beta) - np.log(alpha * beta)


@assumes(GAMMA_GAMMA_ALL_STRENGTH)
def gamma_gamma_pdf(I, alpha, beta):
    '''
    The gamma-gamma PDF of the normalised irradiance.

    Ch. 9, Eq. (137), printed 370:
        p(I) = 2 (alpha beta)^((alpha+beta)/2) / [Gamma(alpha) Gamma(beta)]
               I^((alpha+beta)/2 - 1) K_(alpha-beta)(2 sqrt(alpha beta I)),
    for I > 0, where K_nu is a modified Bessel function of the second kind. The
    book restates it as Ch. 9, Eq. (159), printed 384, as Ch. 11, Eq. (21),
    printed 450, and as Ch. 12, Eq. (67), printed 510.
    DOI: 10.1117/3.626196

    The code evaluates the PDF in the log domain with the scaled Bessel function
    kve, then exponentiates. In the far lower tail the two factors overflow and
    underflow at the same time. The code sets the density to zero there. That is
    correct while min(alpha, beta) >= 1, which the book's variance models always
    give: Ch. 9, Eq. (46), printed 336, keeps sigma_lnY^2 below ln(2).

    Parameters:
        I : float or ndarray
            Normalised irradiance. A value at or below zero gives density zero.
        alpha : float
            Large-scale parameter, alpha > 0.
        beta : float
            Small-scale parameter, beta > 0.

    Returns:
        float or ndarray
            The probability density, same shape as I.
    '''
    shape = np.shape(I)
    Iv = np.atleast_1d(np.asarray(I, dtype=float))
    p = np.zeros_like(Iv)
    m = Iv > 0.0
    if np.any(m):
        ab = alpha * beta
        nu = alpha - beta
        half = 0.5 * (alpha + beta)
        z = 2.0 * np.sqrt(ab * Iv[m])
        log_pre = (np.log(2.0) + half * np.log(ab)
                   - gammaln(alpha) - gammaln(beta))
        with np.errstate(over='ignore', under='ignore', invalid='ignore',
                         divide='ignore'):
            # K_nu(z) = kve(nu, z) exp(-z), so ln K_nu = ln kve - z.
            log_p = (log_pre + (half - 1.0) * np.log(Iv[m])
                     + np.log(kve(nu, z)) - z)
            val = np.exp(log_p)
        p[m] = np.where(np.isfinite(val), val, 0.0)
    return p.reshape(shape) if shape else float(p[0])


@lru_cache(maxsize=64)
def _gamma_gamma_grid(alpha, beta):
    '''
    Cache the log-spaced irradiance grid and the gamma-gamma CDF on it.

    The book gives a closed-form CDF at Ch. 9, Eq. (140), printed 371, and again
    at Ch. 11, Eq. (26), printed 452, as a sum of two 1F2 hypergeometric
    functions. The book warns twice that the hypergeometric functions give
    numerical errors for some parameter values: at printed 452 and at printed
    511. Both places tell the reader to integrate the PDF numerically instead.
    This code follows that instruction.
    DOI: 10.1117/3.626196

    Returns:
        (ndarray, ndarray)
            (I grid, cumulative probability on that grid).
    '''
    I = np.logspace(_GRID_LO_LOG10, _GRID_HI_LOG10, _GRID_N)
    pdf = gamma_gamma_pdf(I, alpha, beta)
    cdf = cumulative_trapezoid(pdf, I, initial=0.0)
    return I, cdf


@assumes(GAMMA_GAMMA_ALL_STRENGTH)
def gamma_gamma_cdf(I, alpha, beta):
    '''
    Pr(I <= value) for a gamma-gamma irradiance.

    The code integrates the PDF of Ch. 9, Eq. (137), printed 370, on a
    log-spaced grid, then interpolates. See _gamma_gamma_grid for the reason the
    code does not use the closed form of Ch. 9, Eq. (140), printed 371.
    DOI: 10.1117/3.626196

    Parameters:
        I : float or ndarray
            Normalised irradiance.
        alpha : float
            Large-scale parameter.
        beta : float
            Small-scale parameter.

    Returns:
        float or ndarray
            The cumulative probability, same shape as I.
    '''
    grid, cdf = _gamma_gamma_grid(float(alpha), float(beta))
    out = np.interp(np.asarray(I, dtype=float), grid, cdf, left=0.0,
                    right=cdf[-1])
    return float(out) if np.ndim(I) == 0 else out


@assumes(GAMMA_GAMMA_ALL_STRENGTH)
def gamma_gamma_quantile(p, alpha, beta):
    '''
    The p-quantile of a gamma-gamma irradiance.

    The code inverts the cached CDF grid by linear interpolation. See
    gamma_gamma_cdf. DOI: 10.1117/3.626196

    Parameters:
        p : float or ndarray
            Cumulative probability in (0, 1).
        alpha : float
            Large-scale parameter.
        beta : float
            Small-scale parameter.

    Returns:
        float or ndarray
            The normalised irradiance at that probability.
    '''
    grid, cdf = _gamma_gamma_grid(float(alpha), float(beta))
    out = np.interp(np.asarray(p, dtype=float), cdf, grid)
    return float(out) if np.ndim(p) == 0 else out


@assumes(GAMMA_GAMMA_ALL_STRENGTH)
def gamma_gamma_rvs(n, alpha, beta, rng):
    '''
    Draw n samples of a gamma-gamma irradiance with E[I] = 1.

    Ch. 9, Eqs. (134) and (135), printed 370, make I = X Y the product of two
    unit-mean gamma variates. X has shape alpha and scale 1/alpha. Y has shape
    beta and scale 1/beta. So the product of two gamma draws gives the
    gamma-gamma directly. DOI: 10.1117/3.626196

    Parameters:
        n : int
            Number of draws.
        alpha : float
            Large-scale parameter.
        beta : float
            Small-scale parameter.
        rng : numpy.random.Generator
            The generator.

    Returns:
        ndarray
            Shape (n,).
    '''
    return rng.gamma(alpha, 1.0 / alpha, n) * rng.gamma(beta, 1.0 / beta, n)


# --- K distribution (Ch. 9.9.1, printed 368-369) ----------------------------
#
# The book states at printed 370, below Eq. (137), that the gamma-gamma
# distribution reduces to the K distribution when alpha or beta is one. So every
# K function here delegates to its gamma-gamma twin at beta = 1. A check in the
# __main__ block proves the identity against Ch. 9, Eq. (132), printed 368.

_K_BETA = 1.0


@assumes(K_STRONG_ONLY, turbulence_regime=REGIME_STRONG)
def k_params(sigma2_I):
    '''
    The K-distribution parameter from a scintillation index.

    Ch. 9, printed 369, states that the K distribution predicts
    sigma_I^2 = 1 + 2/alpha. The inverse is alpha = 2 / (sigma_I^2 - 1).
    DOI: 10.1117/3.626196

    The book states on the same page that this index always exceeds one, so the
    K distribution applies in strong fluctuation regimes only.

    Parameters:
        sigma2_I : float
            Scintillation index. It must exceed one.

    Returns:
        float
            alpha.

    Raises:
        ValueError
            If sigma2_I is one or less. The K distribution cannot reach it.
    '''
    if sigma2_I <= 1.0:
        raise ValueError(
            f"sigma2_I={sigma2_I} is not above 1. The K distribution gives "
            "sigma_I^2 = 1 + 2/alpha, which always exceeds 1 (Andrews and "
            "Phillips, printed 369). Use the gamma-gamma model instead."
        )
    return 2.0 / (sigma2_I - 1.0)


@assumes(K_STRONG_ONLY, turbulence_regime=REGIME_STRONG)
def k_scintillation_index(alpha):
    '''
    The scintillation index of a K-distributed irradiance.

    Ch. 9, printed 369: sigma_I^2 = 1 + 2/alpha.
    DOI: 10.1117/3.626196

    Parameters:
        alpha : float
            The K parameter.

    Returns:
        float
            The scintillation index.
    '''
    return 1.0 + 2.0 / alpha


@assumes(K_STRONG_ONLY, turbulence_regime=REGIME_STRONG)
def k_mean_log(alpha):
    '''
    E[ln I] for a K-distributed irradiance with E[I] = 1.

    The K distribution is the gamma-gamma at beta = 1 (printed 370), so this
    delegates. DOI: 10.1117/3.626196

    Parameters:
        alpha : float
            The K parameter.

    Returns:
        float
            E[ln I].
    '''
    return gamma_gamma_mean_log(alpha, _K_BETA)


@assumes(K_STRONG_ONLY, turbulence_regime=REGIME_STRONG)
def k_pdf(I, alpha):
    '''
    The K-distribution PDF of the normalised irradiance.

    Ch. 9, Eq. (132), printed 368:
        p(I) = 2 alpha / Gamma(alpha) (alpha I)^((alpha-1)/2)
               K_(alpha-1)(2 sqrt(alpha I)), I > 0.
    This equals the gamma-gamma PDF at beta = 1 (printed 370), so the code
    delegates. DOI: 10.1117/3.626196

    Parameters:
        I : float or ndarray
            Normalised irradiance.
        alpha : float
            The K parameter.

    Returns:
        float or ndarray
            The probability density.
    '''
    return gamma_gamma_pdf(I, alpha, _K_BETA)


@assumes(K_STRONG_ONLY, turbulence_regime=REGIME_STRONG)
def k_cdf(I, alpha):
    '''
    Pr(I <= value) for a K-distributed irradiance. See k_pdf and gamma_gamma_cdf.
    DOI: 10.1117/3.626196
    '''
    return gamma_gamma_cdf(I, alpha, _K_BETA)


@assumes(K_STRONG_ONLY, turbulence_regime=REGIME_STRONG)
def k_quantile(p, alpha):
    '''
    The p-quantile of a K-distributed irradiance. See gamma_gamma_quantile.
    DOI: 10.1117/3.626196
    '''
    return gamma_gamma_quantile(p, alpha, _K_BETA)


@assumes(K_STRONG_ONLY, turbulence_regime=REGIME_STRONG)
def k_rvs(n, alpha, rng):
    '''
    Draw n samples of a K-distributed irradiance with E[I] = 1.

    Ch. 9, Eqs. (129) and (131), printed 368, build the K distribution as a
    negative exponential with a gamma-distributed mean. That is the gamma-gamma
    at beta = 1, because a gamma with shape one is the negative exponential.
    DOI: 10.1117/3.626196

    Parameters:
        n : int
            Number of draws.
        alpha : float
            The K parameter.
        rng : numpy.random.Generator
            The generator.

    Returns:
        ndarray
            Shape (n,).
    '''
    return gamma_gamma_rvs(n, alpha, _K_BETA, rng)


# --- lognormal-Rician (Ch. 9.9.2, printed 369) ------------------------------

@assumes(LOGNORMAL_RICIAN_NOT_BUILT)
def lognormal_rician_pdf(I, r, sigma_z2):
    '''
    The lognormal-Rician (Beckmann) PDF of the normalised irradiance.

    Ch. 9, Eq. (133), printed 369:
        p(I) = (1+r) exp(-r) / [sqrt(2 pi) sigma_z]
               INTEGRAL(0..inf) I_0(2 sqrt((1+r) r I / z))
               exp[-(1+r) I / z - (ln z + sigma_z^2/2)^2 / (2 sigma_z^2)]
               dz / z^2,
    where r is the coherence parameter (a power ratio), sigma_z^2 is the
    variance of the lognormal modulation factor, and I_0 is a modified Bessel
    function of the first kind. DOI: 10.1117/3.626196

    THIS FUNCTION GIVES THE PDF ONLY. The book states at printed 369 that a
    closed-form solution for the integral is unknown, and that it is not known
    how to relate r and sigma_z^2 to atmospheric conditions. So this model
    cannot take a turbulence profile as its input, and it has no closed CDF. For
    those two reasons the module gives NO cdf, NO quantile, and NO sampler for
    it. Use it only to compare a measured PDF with a fitted curve. Use the
    gamma-gamma model for a link budget.

    The code integrates over t = ln z, which makes the quadrature stable.

    Parameters:
        I : float or ndarray
            Normalised irradiance, I > 0.
        r : float
            Coherence parameter (power ratio), r >= 0.
        sigma_z2 : float
            Variance of the lognormal modulation factor, sigma_z2 > 0.

    Returns:
        float or ndarray
            The probability density, same shape as I.
    '''
    sigma_z = np.sqrt(sigma_z2)
    lo = -sigma_z2 / 2.0 - 10.0 * sigma_z
    hi = -sigma_z2 / 2.0 + 10.0 * sigma_z
    pre = (1.0 + r) * np.exp(-r) / (np.sqrt(2.0 * np.pi) * sigma_z)

    def one(value):
        if value <= 0.0:
            return 0.0

        def integrand(t):
            # z = exp(t) and dz = z dt, so the dz/z^2 factor becomes dt/z.
            z = np.exp(t)
            u = 2.0 * np.sqrt((1.0 + r) * r * value / z)
            # I_0(u) = i0e(u) exp(u): the scaled form avoids an overflow.
            return (i0e(u) / z
                    * np.exp(u - (1.0 + r) * value / z
                             - (t + sigma_z2 / 2.0) ** 2 / (2.0 * sigma_z2)))

        return pre * quad(integrand, lo, hi, limit=200)[0]

    if np.ndim(I) == 0:
        return one(float(I))
    return np.array([one(float(v)) for v in np.ravel(I)]).reshape(np.shape(I))


# --- the model registry -----------------------------------------------------

MODELS = {
    "lognormal": (lognormal_mean_log, lognormal_quantile, lognormal_rvs),
    "gamma_gamma": (gamma_gamma_mean_log, gamma_gamma_quantile,
                    gamma_gamma_rvs),
    "k": (k_mean_log, k_quantile, k_rvs),
}
'''
Name -> (mean_log, quantile, rvs) for every model that has all three faces.

The calling convention is the same for each model. The first argument is
positional, and every distribution parameter goes in by keyword:

    mean_log(**params)              -> E[ln I]
    quantile(p, **params)           -> the p-quantile of I
    rvs(n, rng=rng, **params)       -> n draws of I

The parameter names per model are:
    "lognormal"    sigma_l2
    "gamma_gamma"  alpha, beta
    "k"            alpha

The lognormal-Rician model is NOT in this dict. It has a PDF only. See
lognormal_rician_pdf.
'''

# Name -> cdf, for probability_of_fade. It is separate from MODELS because
# MODELS holds the three faces that olb/models/fade.py needs to build a Term.
_CDF = {
    "lognormal": lognormal_cdf,
    "gamma_gamma": gamma_gamma_cdf,
    "k": k_cdf,
}


# --- fade statistics (Ch. 11.3, printed 449-457; Ch. 12.7-12.8) -------------

def fade_threshold_irradiance(fade_db):
    '''
    The normalised irradiance at a fade threshold level.

    Ch. 11, Eq. (25), printed 451, and Ch. 12, Eq. (69), printed 511, define
        F_T = 10 log10(E[I] / I_T)  [dB].
    With E[I] = 1 this inverts to I_T = 10^(-F_T/10). The book restates this at
    Ch. 12, Eq. (71), printed 511, as ln(I_T) = -0.23 F_T. The constant 0.23 is
    ln(10)/10 to two figures, and this code uses the exact value.
    DOI: 10.1117/3.626196

    Parameters:
        fade_db : float or ndarray
            Fade threshold parameter F_T [dB] below the mean irradiance.

    Returns:
        float or ndarray
            The threshold irradiance I_T on the normalised scale.
    '''
    return np.exp(-_LN10_OVER_10 * np.asarray(fade_db, dtype=float))


def probability_of_fade(fade_db, model, **params):
    '''
    The probability that the irradiance drops below a fade threshold.

    Ch. 11, Eq. (23), printed 451, defines the probability of fade as the
    cumulative probability Pr(I <= I_T). Ch. 11, Eq. (24), printed 451, gives
    the lognormal closed form. Ch. 11, Eq. (26), printed 452, gives the
    gamma-gamma closed form. Ch. 12, Eq. (70), printed 511, repeats the
    lognormal form for a satellite downlink. This function evaluates the CDF of
    the named model at the threshold of Ch. 11, Eq. (25).
    DOI: 10.1117/3.626196

    Parameters:
        fade_db : float or ndarray
            Fade threshold parameter F_T [dB] below the mean irradiance.
        model : str
            One of the keys of MODELS.
        **params
            The distribution parameters of that model. See MODELS.

    Returns:
        float or ndarray
            The fraction of the time the irradiance stays below the threshold.

    Raises:
        ValueError
            If model is not a known name.
    '''
    if model not in _CDF:
        raise ValueError(f"unknown model {model!r}. Use one of {sorted(_CDF)}.")
    return _CDF[model](fade_threshold_irradiance(fade_db), **params)


def expected_number_of_fades(fade_db, nu0, model, **params):
    '''
    The expected number of fades per second below a threshold.

    Ch. 11, Eq. (34), printed 455, gives the lognormal form. Ch. 11, Eq. (37),
    printed 456, gives the gamma-gamma form. Ch. 12, Eqs. (72) and (74), printed
    513-514, repeat them for a satellite downlink. Every form needs the
    quasi-frequency nu0 of Ch. 11, Eq. (35), printed 456, which comes from the
    second derivative of the temporal irradiance covariance. Get nu0 from
    `olb.turbulence.andrews.temporal.quasi_frequency`.
    DOI: 10.1117/3.626196

    formula (lognormal, Ch. 11, Eq. (34), printed 455):
        <n(I_T)> = nu0 exp{ -[ sigma_l2/2 - 0.23 F_T ]^2 / (2 sigma_l2) }
    The book writes sigma_I in place of sigma_l, and the two agree in the
    weak-fluctuation limit. This code keeps the exact sigma_l, exactly as
    `lognormal_cdf` does, so that the pair (probability, rate) stays coherent.
    At the threshold that makes 0.23 F_T equal sigma_l2/2 the rate becomes nu0.
    The book states that same result at printed 448, below Eq. (15).

    formula (gamma-gamma, Ch. 11, Eq. (37), printed 456):
        <n(I_T)> = 2 sqrt(2 pi alpha beta) nu0 sigma_I / [Gamma(a) Gamma(b)]
                   (alpha beta I_T)^((alpha+beta-1)/2)
                   K_(alpha-beta)(2 sqrt(alpha beta I_T))
    This code evaluates the ALGEBRAIC EQUIVALENT
        <n(I_T)> = sqrt(2 pi) nu0 sigma_I sqrt(I_T) p(I_T),
    with p the gamma-gamma PDF of Ch. 9, Eq. (137), printed 370. The two are the
    same expression: Ch. 11, Eq. (36), printed 456, builds the joint PDF as the
    product of the gamma-gamma PDF and a zero-mean Gaussian for the time
    derivative with variance 4 b I. Ch. 11, Eq. (12) then integrates to
    p(I_T) sqrt(2 b I_T/pi), and Ch. 11, Eq. (38), printed 456, gives
    sqrt(b) = pi nu0 sigma_I. The PDF form is used because `gamma_gamma_pdf`
    works in the log domain, so a large alpha and beta cannot overflow.

    The K distribution is the gamma-gamma model at beta = 1 (printed 370), so
    it delegates.

    Parameters:
        fade_db : float or ndarray
            Fade threshold parameter F_T [dB] below the mean irradiance.
        nu0 : float
            Quasi-frequency [Hz]. The book sets it to 550 Hz for its worked
            figures (printed 457). A real value comes from the temporal
            covariance of the link.
        model : str
            One of the keys of MODELS.
        **params
            The distribution parameters of that model. See MODELS.

    Returns:
        float or ndarray
            The number of down-crossings of the threshold per second.

    Raises:
        ValueError
            If model is not a known name.
    '''
    I_T = fade_threshold_irradiance(fade_db)
    if model == "lognormal":
        sigma_l2 = np.asarray(params["sigma_l2"], dtype=float)
        sigma_l = np.sqrt(sigma_l2)
        arg = (sigma_l2 / 2.0
               - _LN10_OVER_10 * np.asarray(fade_db, dtype=float)) / sigma_l
        return nu0 * np.exp(-arg ** 2 / 2.0)
    if model == "gamma_gamma":
        alpha, beta = float(params["alpha"]), float(params["beta"])
    elif model == "k":
        alpha, beta = float(params["alpha"]), _K_BETA
    else:
        raise ValueError(
            f"unknown model {model!r}. Use one of {sorted(MODELS)}.")
    sigma_I = np.sqrt(gamma_gamma_scintillation_index(alpha, beta))
    return (np.sqrt(2.0 * np.pi) * nu0 * sigma_I * np.sqrt(I_T)
            * gamma_gamma_pdf(I_T, alpha, beta))


def mean_fade_time(fade_db, nu0, model, **params):
    '''
    The mean time the irradiance stays below a threshold.

    Ch. 11, Eq. (39), printed 456, gives the mean fade time as the probability
    of fade divided by the expected number of fades:
        <t(I_T)> = Pr(I <= I_T) / <n(I_T)>.
    Ch. 12, Eqs. (78) and (79), printed 515, repeat it for a satellite downlink.
    Ch. 12, Eq. (79) prints the lognormal ratio in full, which is the same
    quotient that this function forms. It needs the same quasi-frequency as
    `expected_number_of_fades`.
    DOI: 10.1117/3.626196

    The book states at printed 457 that the quasi-frequency changes the expected
    number of fades and the mean fade time, but it does not change the
    probability of fade. So the mean fade time scales as 1/nu0.

    Parameters:
        fade_db : float or ndarray
            Fade threshold parameter F_T [dB] below the mean irradiance.
        nu0 : float
            Quasi-frequency [Hz].
        model : str
            One of the keys of MODELS.
        **params
            The distribution parameters of that model.

    Returns:
        float or ndarray
            The mean duration of one fade [s].

    Raises:
        ValueError
            If model is not a known name.
    '''
    return (probability_of_fade(fade_db, model, **params)
            / expected_number_of_fades(fade_db, nu0, model, **params))


if __name__ == '__main__':
    from scipy.special import gamma as _gamma, kv as _kv

    rng = np.random.default_rng(20250825)

    # === physics ============================================================

    # The gamma-gamma PDF integrates to one and has mean one (Eq. (137), p.370).
    a0, b0 = 4.0, 1.8
    grid = np.logspace(-10, 2, 200_001)
    pdf0 = gamma_gamma_pdf(grid, a0, b0)
    mass = np.trapezoid(pdf0, grid)
    mean_I = np.trapezoid(pdf0 * grid, grid)
    print(f"[physics] gamma-gamma mass      err = {abs(mass - 1.0):.3e}")
    print(f"[physics] gamma-gamma E[I]      err = {abs(mean_I - 1.0):.3e}")
    assert abs(mass - 1.0) < 1e-4 and abs(mean_I - 1.0) < 1e-4

    # Eq. (139), printed 371: the consistency identity ties alpha and beta back
    # to the scintillation index.
    var_I = np.trapezoid(pdf0 * (grid - 1.0) ** 2, grid)
    s2_eq139 = gamma_gamma_scintillation_index(a0, b0)
    print(f"[physics] Eq.(139) index        err = "
          f"{abs(var_I - s2_eq139) / s2_eq139:.3e}")
    assert abs(var_I - s2_eq139) / s2_eq139 < 1e-3

    # Eq. (138), printed 370: the parameter map, checked against Eq. (139).
    s2lnX, s2lnY = 0.30, 0.55
    a1, b1 = gamma_gamma_params(s2lnX, s2lnY)
    want = (np.exp(s2lnX + s2lnY) - 1.0)   # Eq. (28), printed 452
    got = gamma_gamma_scintillation_index(a1, b1)
    print(f"[physics] Eq.(138)->Eq.(139)    err = {abs(got - want) / want:.3e}")
    assert abs(got - want) / want < 1e-12

    # Eq. (132), printed 368: the standalone K PDF, against the delegation.
    aK = 3.0
    Ik = np.array([0.05, 0.3, 1.0, 2.5])
    direct = (2.0 * aK / _gamma(aK) * (aK * Ik) ** ((aK - 1.0) / 2.0)
              * _kv(aK - 1.0, 2.0 * np.sqrt(aK * Ik)))
    print(f"[physics] Eq.(132) K pdf        err = "
          f"{np.max(np.abs(k_pdf(Ik, aK) - direct) / direct):.3e}")
    assert np.allclose(k_pdf(Ik, aK), direct, rtol=1e-10)

    # The K index 1 + 2/alpha (printed 369) agrees with Eq. (139) at beta = 1.
    print(f"[physics] K index vs Eq.(139)   err = "
          f"{abs(k_scintillation_index(aK) - gamma_gamma_scintillation_index(aK, 1.0)):.3e}")
    assert np.isclose(k_scintillation_index(aK),
                      gamma_gamma_scintillation_index(aK, 1.0))

    # Eq. (133), printed 369: the lognormal-Rician PDF integrates to one and has
    # mean one. It has no CDF, no quantile and no sampler, so this is the only
    # test the book allows.
    lr_grid = np.linspace(1e-6, 40.0, 8001)
    lr = lognormal_rician_pdf(lr_grid, r=2.0, sigma_z2=0.3)
    print(f"[physics] Eq.(133) mass         err = "
          f"{abs(np.trapezoid(lr, lr_grid) - 1.0):.3e}")
    print(f"[physics] Eq.(133) E[I]         err = "
          f"{abs(np.trapezoid(lr * lr_grid, lr_grid) - 1.0):.3e}")
    assert abs(np.trapezoid(lr, lr_grid) - 1.0) < 1e-3
    assert abs(np.trapezoid(lr * lr_grid, lr_grid) - 1.0) < 1e-3

    # Eq. (25), printed 451: the fade threshold relation.
    print(f"[physics] Eq.(25) threshold     err = "
          f"{abs(fade_threshold_irradiance(3.0) - 10.0 ** -0.3):.3e}")
    assert np.isclose(fade_threshold_irradiance(3.0), 10.0 ** -0.3)

    # Eq. (24), printed 451: the lognormal probability of fade, against the
    # book's erf form. The book writes sigma_I where this code keeps sigma_l, so
    # the two agree only in the weak-fluctuation limit. Test there.
    from scipy.special import erf
    s2I_weak = 0.05
    sl2 = lognormal_params(s2I_weak)
    F_T = 6.0
    book = 0.5 * (1.0 + erf((0.5 * sl2 - _LN10_OVER_10 * F_T)
                            / (np.sqrt(sl2) * np.sqrt(2.0))))
    code = probability_of_fade(F_T, "lognormal", sigma_l2=sl2)
    # The two forms differ only by the exp/log round trip through the threshold
    # irradiance. At F_T = 6 dB the probability is near 1e-10, so a 1e-17
    # absolute difference looks large in relative terms. Test the absolute one.
    print(f"[physics] Eq.(24) fade prob     err = {abs(code - book):.3e} "
          f"absolute ({abs(code - book) / book:.1e} relative, P={code:.3e})")
    assert abs(code - book) < 1e-15
    # The sign check of the docstring: at F_T = 0 the fade probability exceeds
    # one half, and it equals Phi(sigma_l / 2).
    at0 = probability_of_fade(0.0, "lognormal", sigma_l2=sl2)
    print(f"[physics] Eq.(24) F_T=0 sign    err = "
          f"{abs(at0 - norm.cdf(np.sqrt(sl2) / 2.0)):.3e}   (P={at0:.4f} > 0.5)")
    assert at0 > 0.5 and np.isclose(at0, norm.cdf(np.sqrt(sl2) / 2.0))

    # Monte Carlo moments of every model that has a sampler.
    N = 4_000_000
    for s2I in (0.04, 0.25, 1.0):
        sl2 = lognormal_params(s2I)
        draws = lognormal_rvs(N, sl2, rng)
        e, v = draws.mean(), draws.var()
        print(f"[physics] lognormal   s2I={s2I:<5} E[I] err = "
              f"{abs(e - 1.0):.4f}   Var err = {abs(v - s2I) / s2I:.4f}")
        assert abs(e - 1.0) < 0.01 and abs(v - s2I) / s2I < 0.03

        # Split the two log variances evenly. Eq. (139) then gives
        # s2I = 2 g + g^2 with g = exp(s) - 1, so s = 0.5 ln(1 + s2I).
        s = 0.5 * np.log(1.0 + s2I)
        aa, bb = gamma_gamma_params(s, s)
        draws = gamma_gamma_rvs(N, aa, bb, rng)
        e, v = draws.mean(), draws.var()
        print(f"[physics] gamma-gamma s2I={s2I:<5} E[I] err = "
              f"{abs(e - 1.0):.4f}   Var err = {abs(v - s2I) / s2I:.4f}")
        assert abs(e - 1.0) < 0.01 and abs(v - s2I) / s2I < 0.03

    # The K distribution cannot reach s2I <= 1 (printed 369), so it takes its
    # own set of indices.
    for s2I in (1.5, 2.0, 3.0):
        aK2 = k_params(s2I)
        draws = k_rvs(N, aK2, rng)
        e, v = draws.mean(), draws.var()
        print(f"[physics] K           s2I={s2I:<5} E[I] err = "
              f"{abs(e - 1.0):.4f}   Var err = {abs(v - s2I) / s2I:.4f}")
        assert abs(e - 1.0) < 0.01 and abs(v - s2I) / s2I < 0.03
    try:
        k_params(0.5)
    except ValueError:
        print("[physics] K rejects s2I<=1     ok")
    else:
        raise AssertionError("k_params must refuse sigma2_I <= 1")

    # === reduction ==========================================================

    # The gamma-gamma model goes to the lognormal model in weak fluctuations.
    # The book says so at printed 517 for a small beam. Compare the 1% quantile.
    s2I = 0.04
    sl2 = lognormal_params(s2I)
    s = 0.5 * np.log(1.0 + s2I)
    aa, bb = gamma_gamma_params(s, s)
    q_ln = lognormal_quantile(0.01, sl2)
    q_gg = gamma_gamma_quantile(0.01, aa, bb)
    print(f"[reduce ] GG->LN 1% quantile    err = "
          f"{abs(q_gg - q_ln) / q_ln:.4f}   (LN={q_ln:.4f}, GG={q_gg:.4f})")
    assert abs(q_gg - q_ln) / q_ln < 0.03

    # The K faces are the gamma-gamma faces at beta = 1 (printed 370).
    for face_k, face_gg in ((k_mean_log(aK), gamma_gamma_mean_log(aK, 1.0)),
                            (k_quantile(0.05, aK),
                             gamma_gamma_quantile(0.05, aK, 1.0)),
                            (k_cdf(0.4, aK), gamma_gamma_cdf(0.4, aK, 1.0))):
        assert face_k == face_gg
    print("[reduce ] K == gamma-gamma(b=1) ok")

    # The CDF and the quantile invert each other.
    ps = np.array([0.001, 0.01, 0.1, 0.5, 0.9])
    rt = gamma_gamma_cdf(gamma_gamma_quantile(ps, a0, b0), a0, b0)
    print(f"[reduce ] GG cdf(quantile(p))   err = {np.max(np.abs(rt - ps)):.3e}")
    assert np.max(np.abs(rt - ps)) < 1e-6

    # The numeric gamma-gamma CDF reaches one at the top of the grid.
    _, cgrid = _gamma_gamma_grid(a0, b0)
    print(f"[reduce ] GG cdf(top of grid)   err = {abs(cgrid[-1] - 1.0):.3e}")
    assert abs(cgrid[-1] - 1.0) < 1e-4

    # probability_of_fade agrees with the Monte Carlo fraction below threshold.
    draws = gamma_gamma_rvs(1_000_000, a0, b0, rng)
    for F in (3.0, 6.0, 10.0):
        mc = float(np.mean(draws < fade_threshold_irradiance(F)))
        an = probability_of_fade(F, "gamma_gamma", alpha=a0, beta=b0)
        print(f"[reduce ] P(fade) F_T={F:<5} MC={mc:.5f} analytic={an:.5f} "
              f"err = {abs(an - mc) / mc:.4f}")
        assert abs(an - mc) / mc < 0.03

    # The MODELS registry calls every face by the documented convention.
    for nm, prm in (("lognormal", {"sigma_l2": 0.2}),
                    ("gamma_gamma", {"alpha": 4.0, "beta": 1.8}),
                    ("k", {"alpha": 3.0})):
        ml, qt, rv = MODELS[nm]
        assert np.isfinite(ml(**prm))
        assert np.isfinite(qt(0.01, **prm))
        assert rv(16, rng=rng, **prm).shape[0] == 16
    print("[reduce ] MODELS registry       ok")

    # === fade rate and fade time ===========================================
    #
    # The book gives NO numeric worked example for <n(I_T)> or <t(I_T)>.
    # Ch. 11.7 Example 1, printed 472-473, stops at the probability of fade.
    # Ch. 11 Problem 6, printed 474, asks for both at nu0 = 100 Hz but prints no
    # answer. Ch. 12.10 has no fade-rate example either. So the checks below use
    # the book's own internal identities, not a printed number.

    NU0 = 550.0        # the book's nominal value, printed 457 and printed 514.

    # Eq. (34), printed 455: the rate equals nu0 when 0.23 F_T = sigma_l2/2.
    # The book states that at printed 448, below Eq. (15).
    sl2 = lognormal_params(0.13)
    F_peak = (sl2 / 2.0) / _LN10_OVER_10
    n_peak = expected_number_of_fades(F_peak, NU0, "lognormal", sigma_l2=sl2)
    print(f"[physics] Eq.(34) peak rate     = {n_peak:.6f} Hz (nu0={NU0})   "
          f"err = {abs(n_peak - NU0) / NU0:.3e}")
    assert abs(n_peak - NU0) / NU0 < 1e-12

    # Eq. (34) against the Rice construction of Eqs. (12) and (33), printed 447
    # and 455: <n> = sqrt(2 pi) nu0 sigma_I I_T p_lognormal(I_T). The two agree
    # exactly when both use the same sigma. Eq. (34) prints sigma_I, this code
    # keeps sigma_l, so the check uses sigma_l on both sides.
    for F in (0.0, 3.0, 6.0, 10.0):
        I_T = fade_threshold_irradiance(F)
        pdf_ln = (np.exp(-(np.log(I_T) + sl2 / 2.0) ** 2 / (2.0 * sl2))
                  / (I_T * np.sqrt(2.0 * np.pi * sl2)))
        rice = np.sqrt(2.0 * np.pi) * NU0 * np.sqrt(sl2) * I_T * pdf_ln
        book = expected_number_of_fades(F, NU0, "lognormal", sigma_l2=sl2)
        print(f"[physics] Eq.(34) F_T={F:<5} rate = {book:.6f} Hz vs Rice "
              f"{rice:.6f} Hz   err = {abs(book - rice) / rice:.3e}")
        assert abs(book - rice) / rice < 1e-12

    # Eq. (37), printed 456: the gamma-gamma rate, against the printed form.
    # The printed form overflows for a large alpha, so test at a small pair.
    a2, b2 = 4.0, 2.5
    s_I = np.sqrt(gamma_gamma_scintillation_index(a2, b2))
    for F in (0.0, 3.0, 6.0):
        I_T = fade_threshold_irradiance(F)
        printed = (2.0 * np.sqrt(2.0 * np.pi * a2 * b2) * NU0 * s_I
                   / (_gamma(a2) * _gamma(b2))
                   * (a2 * b2 * I_T) ** ((a2 + b2 - 1.0) / 2.0)
                   * _kv(a2 - b2, 2.0 * np.sqrt(a2 * b2 * I_T)))
        code = expected_number_of_fades(F, NU0, "gamma_gamma", alpha=a2,
                                        beta=b2)
        print(f"[physics] Eq.(37) F_T={F:<5} rate = {code:.6f} Hz vs printed "
              f"{printed:.6f} Hz   err = {abs(code - printed) / printed:.3e}")
        assert abs(code - printed) / printed < 1e-10

    # Eq. (39), printed 456: the identity Pr(fade) = <n> <t>.
    for nm, prm in (("lognormal", {"sigma_l2": sl2}),
                    ("gamma_gamma", {"alpha": a2, "beta": b2}),
                    ("k", {"alpha": 3.0})):
        for F in (1.0, 5.0, 12.0):
            p = probability_of_fade(F, nm, **prm)
            n = expected_number_of_fades(F, NU0, nm, **prm)
            t = mean_fade_time(F, NU0, nm, **prm)
            print(f"[physics] Eq.(39) {nm:<11} F_T={F:<5} P={p:.6e} "
                  f"<n>={n:.4e} <t>={t:.4e}   err = {abs(n * t - p) / p:.3e}")
            assert abs(n * t - p) / p < 1e-9

    # === reduction ==========================================================

    # The gamma-gamma rate goes to the lognormal rate as the scintillation index
    # falls, but ONLY near the mean irradiance. The rate reads the PDF at the
    # threshold, not the CDF, so the deep tail separates the two models much
    # faster than the quantile check above does. The book states at printed 452
    # that the lognormal underestimates the tail. The loop measures both.
    for s2I in (0.04, 0.01, 0.002):
        sl2w = lognormal_params(s2I)
        s = 0.5 * np.log(1.0 + s2I)
        aw, bw = gamma_gamma_params(s, s)
        n_ln = expected_number_of_fades(0.0, NU0, "lognormal", sigma_l2=sl2w)
        n_gg = expected_number_of_fades(0.0, NU0, "gamma_gamma", alpha=aw,
                                        beta=bw)
        t_ln = expected_number_of_fades(3.0, NU0, "lognormal", sigma_l2=sl2w)
        t_gg = expected_number_of_fades(3.0, NU0, "gamma_gamma", alpha=aw,
                                        beta=bw)
        print(f"[reduce ] GG/LN rate s2I={s2I:<6} at F_T=0 dB = "
              f"{n_gg / n_ln:.4f}, at F_T=3 dB = {t_gg / t_ln:.3f} (tail)")
        assert abs(n_gg - n_ln) / n_ln < 0.01

    # The K faces are the gamma-gamma faces at beta = 1 (printed 370).
    assert (expected_number_of_fades(4.0, NU0, "k", alpha=aK)
            == expected_number_of_fades(4.0, NU0, "gamma_gamma", alpha=aK,
                                        beta=1.0))
    print("[reduce ] K rate == GG(beta=1)  ok")

    # The mean fade time scales as 1/nu0, which the book states at printed 457.
    t1 = mean_fade_time(6.0, 100.0, "lognormal", sigma_l2=sl2)
    t2 = mean_fade_time(6.0, 200.0, "lognormal", sigma_l2=sl2)
    print(f"[reduce ] <t> nu0 scaling ratio = {t1 / t2:.6f} (2 wanted)")
    assert np.isclose(t1 / t2, 2.0)

    # === assumptions ========================================================
    import warnings as _warnings

    from ...assumptions import trace_assumptions

    # (1) VALUE PARITY. The decorator does not change the number, in or out of
    #     a collection context.
    _out = gamma_gamma_quantile(0.01, 4.0, 1.8)
    with trace_assumptions():
        _in = gamma_gamma_quantile(0.01, 4.0, 1.8)
    assert _out == _in, (_out, _in)
    _lo = lognormal_quantile(0.01, 0.2)
    with trace_assumptions():
        _li = lognormal_quantile(0.01, 0.2)
    assert _lo == _li, (_lo, _li)

    # (2) REGISTRATION. The three builders register, with the expected kinds.
    with trace_assumptions() as _tr:
        _sl2 = lognormal_params(0.1)                 # weak, inside the PDF limit
        lognormal_quantile(0.01, _sl2)
        _a, _b = gamma_gamma_params(0.2, 0.3)
        gamma_gamma_quantile(0.01, _a, _b)
        k_quantile(0.05, k_params(2.0))              # strong only
    for _src in ("lognormal_params", "gamma_gamma_params", "k_params"):
        assert f"{__name__}.{_src}" in _tr.records, _src
    _kinds = {c.kind for r in _tr.records.values() for c in r.constraints}
    assert {"pdf-shape", "regime"} <= _kinds, _kinds
    # Weak turbulence inside the PDF limit breaks nothing.
    assert not _tr.violations, _tr.violations

    # (3) VIOLATION, NO WARNING. A lognormal face past the PDF-shape limit
    #     yields a source-prefixed violation, and the layer warns about nothing.
    with _warnings.catch_warnings(record=True) as _caught:
        _warnings.simplefilter("always")
        with trace_assumptions() as _tr:
            lognormal_quantile(0.01, lognormal_params(0.9))   # sigma_I^2 = 0.9
    assert any(f"[{__name__}.lognormal_params]" in v and "sigma_I^2" in v
               for v in _tr.violations), _tr.violations
    assert len(_caught) == 0, _caught
    print("[assume] distributions registration, kinds and violation ok")

    print("self-check passed")
