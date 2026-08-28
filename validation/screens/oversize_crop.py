'''
Arm 1 of the phase-screen low-frequency study: the OVERSIZE-and-CROP route.

A Fourier phase screen holds no power below the grid fundamental 1/(N dx). The
tilt is the mode that suffers most, because the Z-tilt integrand grows as
f^(-2/3) below 1/D. Schmidt (2010), DOI 10.1117/3.866274, Ch. 9, text below
Listing 9.2, printed p. 167, states the defect. Two cures exist:

  1. THE SUBHARMONICS. Add three low-frequency grids at 1/(3^p N dx). That is
     Schmidt Ch. 9, Eq. (9.81), printed p. 169, from Lane, Glindemann and
     Dainty, DOI 10.1088/0959-7174/2/3/003.
  2. THE OVERSIZE. Make the screen on a grid that is K times wider, at the SAME
     pitch, and then cut the centre out. The grid fundamental falls by K, and
     the Nyquist frequency does not move. The book does not name this route.

This script measures both cures against the analytic truth, on two spectra
(Kolmogorov and von Karman with L0 = 25 m), and it reports three things:

  1. THE OVERSIZE SWEEP. The Z-tilt and the G-tilt variance of each arm, and
     the ensemble structure function of each arm.
  2. THE S-27 SETTLEMENT. The docs of this repository disagree about the
     aotools subharmonic screen against the book subharmonic screen. This
     script measures the two with BOTH structure-function estimators on ONE
     shared ensemble, so the disagreement gets a number.
  3. THE aotools SEED QUIRK. The aotools generator `ft_sh_phase_screen`
     gives its own subharmonic generator the SAME integer seed that it passes
     down to `ft_phase_screen`. So the first 9 low-frequency Gaussian draws are
     bit-identical to the first 9 high-frequency draws. This script measures
     whether that shared seed biases the tilt variance.

THE TRUTH THAT EVERY ARM IS MEASURED AGAINST. `validation/screens/helpers.py`
holds it. The Z-tilt variance comes from the Noll spatial filter,
DOI 10.1364/JOSA.66.000207, Eq. (8), p. 208. The G-tilt variance comes from
Andrews and Phillips, DOI 10.1117/3.626196, Ch. 6, Eqs. (82) to (84), printed
pp. 200 and 201. The sharp-cutoff prediction of a finite screen comes from
`helpers.captured_fraction`.

THE PITCH IS FIXED. Every arm runs at dx = 1 cm. The oversize adds low
frequencies only. It does not move the Nyquist frequency, so the two cures are
compared on one high-frequency band.

THE MEASUREMENT IS ON THE CROPPED PRODUCT. The oversize arms crop to 512 pixels
before any estimate. A wave-optics run would use the crop, so the crop is the
object under test.

THE THEORY GOES THROUGH THE SAME BINS. `helpers.ensemble_dphi` ends in
`helpers.radial_average`, which returns the mean of an annulus of finite width.
The structure function grows as r^(5/3), so that mean is ABOVE the value at the
bin centre. A comparison against the bin-centre value then reads about 5
percent high. This script puts the analytic law through the same annuli. See
`bin_theory`.

THE PASS BANDS. Four bands assert below. A failed band prints FAIL and the
script continues, because a physics surprise IS a result of this study.

Outputs, all next to this script:
    oversize_tilt.csv     the tilt variance of every arm
    oversize_dphi.csv     the structure function of every arm
    s27_settlement.csv    the four generators, two estimators
    seed_quirk.csv        the aotools shared-seed measurement
    oversize_tilt.png     the Z-tilt ratio against the oversize factor
    oversize_dphi.png     the structure-function ratio against the separation

Sources:
- Schmidt, Numerical Simulation of Optical Wave Propagation with Examples in
  MATLAB, DOI 10.1117/3.866274. Ch. 3, Eqs. (3.19) to (3.25), printed pp. 49
  and 50 (the structure-function estimator); Ch. 3, Eq. (3.15), printed p. 47
  (the definition); Ch. 9, Eq. (9.44), printed p. 160
  (D(r) = 6.88 (r/r0)^(5/3));
  Ch. 9, Eqs. (9.50) and (9.52), printed p. 161 (the phase PSDs); Ch. 9,
  Eqs. (9.78) to (9.80), printed p. 167 (the Fourier screen); Ch. 9,
  Eq. (9.81), printed p. 169 (the subharmonics).
- Noll, J. Opt. Soc. Am. 66(3), pp. 207 to 211 (1976),
  DOI 10.1364/JOSA.66.000207. The Zernike tilt filter.
- Andrews and Phillips, Laser Beam Propagation through Random Media, 2nd ed.,
  DOI 10.1117/3.626196. Ch. 6, Eqs. (80) to (84), printed pp. 200 and 201. The
  gradient tilt.
- Lane, Glindemann and Dainty, Waves in Random Media 2, pp. 209 to 224 (1992),
  DOI 10.1088/0959-7174/2/3/003. The subharmonic method.
- McGlamery, J. Opt. Soc. Am. 57(3), pp. 293 to 297 (1967),
  DOI 10.1364/JOSA.57.000293. The Fourier phase screen.
- Assemat and Wilson, Opt. Express 14(3), pp. 988 to 999 (2006),
  DOI 10.1364/OE.14.000988, Eq. (5). The closed-form von Karman covariance,
  through `helpers.vk_covariance_closed`.

Run from the repo root:
    python -m validation.screens.oversize_crop
'''

import csv
import time

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt                            # noqa: E402
import numpy as np                                         # noqa: E402

from olb.waveoptics.schmidt.turbulence import (             # noqa: E402
    ft_phase_screen, ft_sh_phase_screen, kolmogorov_phase_psd,
    von_karman_phase_psd)
from olb.waveoptics.turbulence.screens import phase_screen  # noqa: E402

from validation.screens import helpers                      # noqa: E402


# ---------------------------------------------------------------------------
# The module constants
# ---------------------------------------------------------------------------

WAVELENGTH_M = 1550e-9
K_RAD_M = 2.0 * np.pi / WAVELENGTH_M

R0_M = 0.10                     # the Fried parameter of every screen
DX_M = 0.01                     # the pitch. It is FIXED in every arm.
N_REF = 512                     # the crop target
SIDE_REF_M = N_REF * DX_M       # 5.12 m
PUPIL_D_M = 1.0                 # the pupil that reads the tilt

# The oversize factors. The grid side is K * 5.12 m, and the pitch does not
# move. So the grid fundamental 1/(K * 5.12) falls by K.
OVERSIZE_K = (1, 2, 4, 8)

# The trials per oversize factor. A large grid costs more, so it gets fewer
# trials. The standard error of a variance is var * sqrt(2 / (2 M)).
TRIALS_KOLMOGOROV = (100, 100, 60, 40)
TRIALS_VON_KARMAN = (60, 60, 40, 30)

L0_VON_KARMAN_M = 25.0
N_SUBHARMONIC = 3               # the book value, Eq. (9.81), printed p. 169

MASTER_SEED = 20260827

# The structure-function window. `schmidt.fourier.structure_function` makes the
# correlation CIRCULAR, so the window must leave a guard band at the grid edge.
# A radius of one quarter of the side keeps the estimate valid out to 2.56 m.
SF_MASK_RADIUS_M = 1.28
DPHI_BIN_COUNT = 12
DPHI_MIN_R0 = 0.5
DPHI_MAX_R0 = 20.0

# The S-27 settlement, Kolmogorov only.
S27_TRIALS = 40
S27_PROBES_R0 = (1.0, 2.0, 4.0, 6.0, 8.0)

# The aotools seed quirk, Kolmogorov only. It reads the S-27 aotools ensemble.
QUIRK_TRIALS = S27_TRIALS

# The pass bands. Band 2 loosens to 3 SE for the K = 8 arms, because those
# arms carry the fewest trials. The loosened band prints below.
BAND_SE = 2.0
BAND_SE_K8 = 3.0
VK_K8_TILT_BAND = (0.85, 1.15)
VK_K8_DPHI_BAND = (0.90, 1.05)

OUT_DIR = 'validation/screens'
TILT_CSV = f'{OUT_DIR}/oversize_tilt.csv'
DPHI_CSV = f'{OUT_DIR}/oversize_dphi.csv'
S27_CSV = f'{OUT_DIR}/s27_settlement.csv'
QUIRK_CSV = f'{OUT_DIR}/seed_quirk.csv'
TILT_PNG = f'{OUT_DIR}/oversize_tilt.png'
DPHI_PNG = f'{OUT_DIR}/oversize_dphi.png'

# The two spectra. `L0_m` goes to the screen generators, and `psd` goes to the
# analytic filters of `helpers`.
SPECTRA = {
    'kolmogorov': {
        'L0_m': np.inf,
        'trials': TRIALS_KOLMOGOROV,
        'psd': lambda f: kolmogorov_phase_psd(f, R0_M),
        'label': 'Kolmogorov, L0 infinite',
    },
    'von karman': {
        'L0_m': L0_VON_KARMAN_M,
        'trials': TRIALS_VON_KARMAN,
        'psd': lambda f: von_karman_phase_psd(f, R0_M, L0_VON_KARMAN_M),
        'label': f'von Karman, L0 = {L0_VON_KARMAN_M:.0f} m',
    },
}

ARM_FT = 'book ft'
ARM_FT_SH = 'book ft_sh'
ARM_CROP = 'book ft crop'
ARM_AOTOOLS = 'aotools ft_sh'

# The seed streams. Each name gets its own integer master from ONE master seed.
# The two paired arms share the name `pair`, so they see the same draw.
STREAM_NAMES = (
    'kolmogorov/pair', 'kolmogorov/crop2', 'kolmogorov/crop4',
    'kolmogorov/crop8', 'kolmogorov/aotools',
    'von karman/pair', 'von karman/crop2', 'von karman/crop4',
    'von karman/crop8', 'von karman/aotools',
    's27/book', 's27/aotools',
)
STREAM_SEED = dict(zip(STREAM_NAMES,
                       helpers.spawn_seeds(MASTER_SEED, len(STREAM_NAMES))))


# ---------------------------------------------------------------------------
# The shared grids
# ---------------------------------------------------------------------------

PUPIL_MASK = helpers.pupil_mask(N_REF, DX_M, PUPIL_D_M)
PUPIL_D_EFF_M = helpers.mask_diameter(PUPIL_MASK, DX_M)
SF_MASK = helpers.pupil_mask(N_REF, DX_M, 2.0 * SF_MASK_RADIUS_M)

# The separation bins. The estimate holds out to the window DIAMETER, so cap
# the largest bin at 2.56 m.
_R_MAX_M = min(DPHI_MAX_R0 * R0_M, 2.0 * SF_MASK_RADIUS_M)
R_BINS_M = np.logspace(np.log10(DPHI_MIN_R0 * R0_M), np.log10(_R_MAX_M),
                       DPHI_BIN_COUNT)


def z_tilt_target(psd):
    '''Return the analytic per-axis Z-tilt ANGLE variance [rad^2].

    `helpers.tilt_filter_variance` gives the Noll coefficient variance
    <a2^2>. Divide by (k D / 4)^2 to get the angle, as the docstring of
    `helpers.zernike_tilt` states. Source: Noll,
    DOI 10.1364/JOSA.66.000207, Eq. (8) and Table I, pp. 208 and 209.
    '''
    a2_var = helpers.tilt_filter_variance(psd, PUPIL_D_EFF_M)
    return a2_var / (K_RAD_M * PUPIL_D_EFF_M / 4.0) ** 2


def g_tilt_target(psd):
    '''Return the analytic per-axis G-tilt ANGLE variance [rad^2].

    Source: Andrews and Phillips, DOI 10.1117/3.626196, Ch. 6, Eqs. (82) to
    (84), printed pp. 200 and 201, through `helpers.gtilt_filter_variance`.
    '''
    return helpers.gtilt_filter_variance(psd, PUPIL_D_EFF_M, WAVELENGTH_M)


# The separation map of the 512 grid. `bin_theory` reads it.
_SEP_AXIS_M = (np.arange(N_REF) - N_REF // 2) * DX_M
_SEP_R_M = np.hypot(*np.meshgrid(_SEP_AXIS_M, _SEP_AXIS_M))


def dphi_law(spectrum_key):
    '''Return the analytic structure function as a callable of r [m].

    The Kolmogorov branch is Schmidt (2010), DOI 10.1117/3.866274, Ch. 9,
    Eq. (9.44), printed p. 160: D(r) = 6.88 (r/r0)^(5/3). The von Karman
    branch is D(r) = 2 [B(0) - B(r)], with B from the closed form of Assemat
    and Wilson, DOI 10.1364/OE.14.000988, Eq. (5), through
    `helpers.vk_covariance_closed`. The relation D = 2 [B(0) - B(r)] is
    Schmidt Ch. 3, Eq. (3.16), printed p. 48.
    '''
    if spectrum_key == 'kolmogorov':
        return lambda r: 6.88 * (np.asarray(r) / R0_M) ** (5.0 / 3.0)
    b_zero = float(helpers.vk_covariance_closed(0.0, R0_M,
                                                L0_VON_KARMAN_M)[0])
    return lambda r: 2.0 * (b_zero - helpers.vk_covariance_closed(
        r, R0_M, L0_VON_KARMAN_M))


def bin_theory(spectrum_key, bins_m):
    '''Return the analytic structure function through the SAME bins.

    Parameters:
        spectrum_key : str
            A key of `SPECTRA`.
        bins_m : array_like
            The separation bin centres [m].

    Returns:
        numpy.ndarray
            One theory value per bin.

    WHY THE THEORY MUST GO THROUGH THE BINS. `helpers.ensemble_dphi` ends in
    `helpers.radial_average`, which returns the MEAN of an annulus of finite
    width. The structure function grows as r^(5/3), so the annulus mean sits
    ABOVE the value at the bin centre. The self-check of
    `validation/screens/helpers.py`, section 7, measures the same effect on an
    exact ramp. This function removes that bias: it puts the analytic law on
    the separation grid and it reads the same annuli.
    '''
    return helpers.radial_average(dphi_law(spectrum_key)(_SEP_R_M), DX_M,
                                  bins_m)


def dphi_theory(spectrum_key):
    '''Return the binned analytic structure function on `R_BINS_M`.'''
    return bin_theory(spectrum_key, R_BINS_M)


# ---------------------------------------------------------------------------
# 1. The oversize sweep
# ---------------------------------------------------------------------------

def run_arm(make_screen, trials):
    '''Measure one arm: the two tilts and the structure function.

    Parameters:
        make_screen : callable
            It takes the trial index and it returns ONE 512 by 512 screen
            [rad]. The caller does the crop.
        trials : int
            The number of trials.

    Returns:
        dict
            The pooled tilt variances, the tilt samples, and the binned
            structure function.

    THE MEMORY. The generator makes one screen at a time, and
    `helpers.ensemble_dphi` accumulates the structure function as it reads.
    So the script never holds an ensemble of large screens.

    THE POOL. The two axes give 2 M samples of one variance. The mean of the
    squares is the estimate, because the true mean of a tilt is zero.
    '''
    z_angle, g_angle, noll = [], [], []

    def screens():
        for index in range(trials):
            screen = make_screen(index)
            zx, zy, a2, a3 = helpers.zernike_tilt(screen, PUPIL_MASK, DX_M,
                                                  WAVELENGTH_M)
            gx, gy = helpers.gradient_tilt(screen, PUPIL_MASK, DX_M,
                                           WAVELENGTH_M)
            z_angle.extend((zx, zy))
            g_angle.extend((gx, gy))
            noll.extend((a2, a3))
            yield screen

    d_measured = helpers.ensemble_dphi(screens(), SF_MASK, DX_M, R_BINS_M)
    z = np.asarray(z_angle)
    g = np.asarray(g_angle)
    return {
        'trials': trials,
        'z': z,
        'g': g,
        'noll': np.asarray(noll),
        'z_var': float(np.mean(z * z)),
        'g_var': float(np.mean(g * g)),
        'd_measured': d_measured,
    }


def book_screen(spectrum_key, seeds, n, subharmonic):
    '''Return a maker of ONE book screen, cropped to 512 pixels.

    Parameters:
        spectrum_key : str
            A key of `SPECTRA`.
        seeds : list
            One integer seed per trial.
        n : int
            The generation grid side in pixels.
        subharmonic : bool
            True calls `ft_sh_phase_screen`. False calls `ft_phase_screen`.

    Returns:
        callable
            It takes the trial index and it returns a 512 by 512 screen.

    THE PAIR. The two K = 1 book arms take the SAME seed list. Each call
    builds a fresh generator from the trial seed, so the two arms share the
    high-frequency draw. The subharmonic delta is then isolated, because the
    difference of the two arms is the low-frequency screen alone.

    THE COPY. `helpers.crop_center` returns a VIEW, and a view keeps the whole
    generation array alive. This function copies the block, so the large array
    is released at once.
    '''
    L0_m = SPECTRA[spectrum_key]['L0_m']
    generator = ft_sh_phase_screen if subharmonic else ft_phase_screen

    def make(index):
        rng = np.random.default_rng(seeds[index])
        if subharmonic:
            screen = generator(R0_M, n, DX_M, L0_m, 0.0, rng, N_SUBHARMONIC)
        else:
            screen = generator(R0_M, n, DX_M, L0_m, 0.0, rng)
        return np.array(helpers.crop_center(screen, N_REF))

    return make


def aotools_screen(spectrum_key, seeds, subharmonics=True):
    '''Return a maker of ONE production (aotools) screen at 512 pixels.

    `olb.waveoptics.turbulence.screens.phase_screen` takes an INTEGER seed,
    because aotools builds its own generator. An infinite L0 becomes 1e6 m
    inside that function, which gives the Kolmogorov spectrum to 12 digits.
    '''
    L0_m = SPECTRA[spectrum_key]['L0_m']

    def make(index):
        return phase_screen(R0_M, N_REF, DX_M, L0_m=L0_m, seed=seeds[index],
                            subharmonics=subharmonics)

    return make


def run_sweep(spectrum_key):
    '''Run the four arms of one spectrum, and return the results by arm.

    The keys are (arm, K). The K = 1 book pair shares one seed stream.
    '''
    trials = SPECTRA[spectrum_key]['trials']
    out = {}

    pair_seeds = helpers.spawn_seeds(STREAM_SEED[f'{spectrum_key}/pair'],
                                     trials[0])
    out[(ARM_FT, 1)] = run_arm(
        book_screen(spectrum_key, pair_seeds, N_REF, False), trials[0])
    out[(ARM_FT_SH, 1)] = run_arm(
        book_screen(spectrum_key, pair_seeds, N_REF, True), trials[0])

    for k, count in zip(OVERSIZE_K[1:], trials[1:]):
        seeds = helpers.spawn_seeds(STREAM_SEED[f'{spectrum_key}/crop{k}'],
                                    count)
        out[(ARM_CROP, k)] = run_arm(
            book_screen(spectrum_key, seeds, N_REF * k, False), count)

    ao_seeds = helpers.spawn_seeds(STREAM_SEED[f'{spectrum_key}/aotools'],
                                   trials[0])
    out[(ARM_AOTOOLS, 1)] = run_arm(aotools_screen(spectrum_key, ao_seeds),
                                    trials[0])
    return out


# ---------------------------------------------------------------------------
# 2. The S-27 settlement
# ---------------------------------------------------------------------------

def run_s27():
    '''Measure four Kolmogorov generators on ONE 40-trial ensemble.

    The four are the book Fourier screen, the book subharmonic screen, the
    aotools Fourier screen, and the aotools subharmonic screen. The book pair
    shares one generator state per trial, and the aotools pair shares one
    integer seed per trial. So each pair isolates its own subharmonic delta.

    Returns:
        tuple
            (ratios, tilts). `ratios` maps (generator, estimator) to the ratio
            at each probe. `tilts` holds the per-trial aotools Z-tilt angles,
            which the seed-quirk section reads.

    THE TWO ESTIMATORS. The FFT route is `helpers.ensemble_dphi`, from Schmidt
    (2010), DOI 10.1117/3.866274, Ch. 3, Eqs. (3.19) to (3.25), printed pp. 49
    and 50. It applies the pupil window, and it bins over an annulus. The
    direct route is `helpers.d_phi_direct`, from Ch. 3, Eq. (3.15), printed
    p. 47. It reads the whole grid at ONE exact separation, with no window and
    no bin. The two are not the same number for one screen. Compare each one
    against the theory.

    THE TWO THEORY COLUMNS. The direct route reads the exact separation, so it
    takes the exact law. The FFT route reads an annulus, so it takes the law
    through the SAME annulus. See `bin_theory`.
    '''
    book_seeds = helpers.spawn_seeds(STREAM_SEED['s27/book'], S27_TRIALS)
    ao_seeds = helpers.spawn_seeds(STREAM_SEED['s27/aotools'], S27_TRIALS)

    probes_m = np.array(S27_PROBES_R0) * R0_M
    steps_px = np.round(probes_m / DX_M).astype(int)
    theory = 6.88 * (steps_px * DX_M / R0_M) ** (5.0 / 3.0)
    theory_fft = bin_theory('kolmogorov', probes_m)

    makers = {
        'book ft': book_screen('kolmogorov', book_seeds, N_REF, False),
        'book ft_sh': book_screen('kolmogorov', book_seeds, N_REF, True),
        'aotools ft': aotools_screen('kolmogorov', ao_seeds, False),
        'aotools ft_sh': aotools_screen('kolmogorov', ao_seeds, True),
    }

    ratios, tilts = {}, {}
    for name, make in makers.items():
        direct = np.zeros(len(S27_PROBES_R0))
        z_x, z_y = [], []

        def screens():
            for index in range(S27_TRIALS):
                screen = make(index)
                for j, step in enumerate(steps_px):
                    direct[j] += helpers.d_phi_direct(screen, int(step))
                zx, zy, _, _ = helpers.zernike_tilt(screen, PUPIL_MASK, DX_M,
                                                    WAVELENGTH_M)
                z_x.append(zx)
                z_y.append(zy)
                yield screen

        fft = helpers.ensemble_dphi(screens(), SF_MASK, DX_M, probes_m)
        direct /= S27_TRIALS
        ratios[(name, 'fft')] = fft / theory_fft
        ratios[(name, 'direct')] = direct / theory
        tilts[name] = np.column_stack([z_x, z_y])
    return ratios, {'fft': theory_fft, 'direct': theory}, tilts


# ---------------------------------------------------------------------------
# 3. The aotools seed quirk
# ---------------------------------------------------------------------------

def run_seed_quirk(tilts):
    '''Measure the effect of the shared aotools seed on the tilt variance.

    Parameters:
        tilts : dict
            The per-trial aotools Z-tilt angles of `run_s27`. Each value is a
            (trials, 2) array of the two axes.

    Returns:
        dict
            The correlation, the two variances, and their standard errors.

    THE QUIRK. `aotools.turbulence.phasescreen.ft_sh_phase_screen` builds
    `numpy.random.default_rng(seed)` for its subharmonics, and it passes the
    SAME integer seed down to `ft_phase_screen`. So the 9 low-frequency real
    draws are bit-identical to the first 9 high-frequency real draws.

    THE LINEARITY. `helpers.zernike_tilt` is a least-squares fit, so it is
    LINEAR in the screen. The tilt of the low-frequency screen is then the
    difference of the two measured tilts, and no screen has to be kept.

        LF_i = ft_sh(seed_i) - ft(seed_i)   so  tilt(LF_i) = z_sh_i - z_ft_i
        HF_i = ft(seed_i)                   so  tilt(HF_i) = z_ft_i

    THE DECORRELATED CONTROL. Add the low-frequency tilt of a SHIFTED trial:
    HF_i + LF_j with j = i + 1, modulo the trial count. That pairing breaks the
    shared seed and it keeps both marginal distributions.
    '''
    hf = tilts['aotools ft'].ravel()
    lf = (tilts['aotools ft_sh'] - tilts['aotools ft']).ravel()

    correlation = float(np.corrcoef(hf, lf)[0, 1])

    shifted = np.roll(tilts['aotools ft_sh'] - tilts['aotools ft'], 1,
                      axis=0).ravel()
    var_correlated = float(np.mean((hf + lf) ** 2))
    var_decorrelated = float(np.mean((hf + shifted) ** 2))

    count = hf.size
    return {
        'count': count,
        'correlation': correlation,
        'correlation_se': 1.0 / np.sqrt(QUIRK_TRIALS),
        'var_correlated': var_correlated,
        'var_decorrelated': var_decorrelated,
        'var_se': var_decorrelated * np.sqrt(2.0 / count),
        'var_hf': float(np.mean(hf * hf)),
        'var_lf': float(np.mean(lf * lf)),
    }


# ---------------------------------------------------------------------------
# The figures
# ---------------------------------------------------------------------------

def draw_tilt(sweep):
    '''Draw the Z-tilt ratio against the oversize factor, one panel per
    spectrum.'''
    fig, axes = plt.subplots(1, 2, figsize=(13.8, 5.8),
                             constrained_layout=True)
    k_curve = np.logspace(np.log10(0.85), np.log10(11.0), 60)

    for ax, key in zip(axes, SPECTRA):
        psd = SPECTRA[key]['psd']
        target = z_tilt_target(psd)
        arms = sweep[key]

        theory = np.array([helpers.captured_fraction(
            psd, PUPIL_D_EFF_M, 1.0 / (k * SIDE_REF_M)) for k in k_curve])
        ax.semilogx(k_curve, theory, color='black', linewidth=2.2,
                    label='sharp cutoff, captured_fraction at 1/(K * 5.12 m)')

        k_points, ratio, error = [], [], []
        for (arm, k), result in sorted(arms.items(), key=lambda kv: kv[0][1]):
            if arm not in (ARM_FT, ARM_CROP):
                continue
            k_points.append(k)
            ratio.append(result['z_var'] / target)
            error.append(BAND_SE * result['z_var'] / target
                         * np.sqrt(1.0 / result['trials']))
        ax.errorbar(k_points, ratio, yerr=error, fmt='o', color='tab:blue',
                    markersize=7.0, capsize=4.0, linewidth=1.6,
                    label='book ft, oversized then cropped (+- 2 SE)')

        sub = arms[(ARM_FT_SH, 1)]
        sub_ratio = sub['z_var'] / target
        sub_se = sub_ratio * np.sqrt(1.0 / sub['trials'])
        ax.axhspan(sub_ratio - BAND_SE * sub_se, sub_ratio + BAND_SE * sub_se,
                   color='tab:orange', alpha=0.25, zorder=0,
                   label=f'book ft_sh at K = 1, {sub_ratio:.3f} (+- 2 SE)')

        ao = arms[(ARM_AOTOOLS, 1)]
        ao_ratio = ao['z_var'] / target
        ax.errorbar([1.0], [ao_ratio],
                    yerr=[BAND_SE * ao_ratio * np.sqrt(1.0 / ao['trials'])],
                    fmt='s', color='tab:red', markersize=8.0, capsize=4.0,
                    label=f'aotools ft_sh at K = 1, {ao_ratio:.3f}')

        ax.axhline(1.0, color='0.4', linewidth=1.0, linestyle=':')
        ax.set_xscale('log', base=2)
        ax.set_xticks(OVERSIZE_K)
        ax.set_xticklabels([str(k) for k in OVERSIZE_K])
        ax.set_xlabel('oversize factor K (grid side = K * 5.12 m)')
        ax.set_ylabel('measured Z-tilt variance / analytic')
        ax.set_title(SPECTRA[key]['label'], fontsize=11)
        ax.grid(alpha=0.3, which='both')
        ax.legend(fontsize=8, loc='lower right')

    fig.suptitle(
        f'The low-frequency tilt of a Fourier phase screen. r0 = '
        f'{R0_M * 1e2:.0f} cm, pitch {DX_M * 1e2:.0f} cm, pupil '
        f'{PUPIL_D_M:.1f} m, crop to {N_REF} px.\n'
        f'The oversize adds low frequency only. The pitch, and so the Nyquist '
        f'frequency, does not move.', fontsize=12)
    fig.savefig(TILT_PNG, dpi=150)
    plt.close(fig)


def draw_dphi(sweep):
    '''Draw the structure-function ratio against the separation.'''
    fig, axes = plt.subplots(1, 2, figsize=(13.8, 5.8),
                             constrained_layout=True)
    colours = {(ARM_FT, 1): 'tab:green', (ARM_FT_SH, 1): 'tab:orange',
               (ARM_CROP, 2): '#9ecae1', (ARM_CROP, 4): '#4292c6',
               (ARM_CROP, 8): '#084594', (ARM_AOTOOLS, 1): 'tab:red'}

    for ax, key in zip(axes, SPECTRA):
        theory = dphi_theory(key)
        for (arm, k), colour in colours.items():
            result = sweep[key][(arm, k)]
            label = arm if k == 1 else f'{arm}, K = {k}'
            ax.semilogx(R_BINS_M / R0_M, result['d_measured'] / theory,
                        color=colour, linewidth=1.8, marker='o',
                        markersize=3.5, label=label)
        ax.axhline(1.0, color='black', linewidth=2.0, label='theory')
        ax.set_ylim(0.0, 1.25)
        ax.set_xlabel('separation r / r0')
        ax.set_ylabel('measured D(r) / theory')
        ax.set_title(SPECTRA[key]['label'], fontsize=11)
        ax.grid(alpha=0.3, which='both')
        ax.legend(fontsize=8, loc='lower left')

    fig.suptitle(
        'The ensemble structure function of every arm, as a ratio to the '
        'analytic law.\n'
        'Kolmogorov theory is Schmidt Ch. 9, Eq. (9.44), printed p. 160. The '
        'von Karman theory is 2 [B(0) - B(r)].', fontsize=12)
    fig.savefig(DPHI_PNG, dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
# The output files
# ---------------------------------------------------------------------------

def write_tilt_csv(sweep):
    '''Write the tilt table of every arm.'''
    with open(TILT_CSV, 'w', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(['spectrum', 'arm', 'K', 'M', 'metric', 'var_rad2',
                         'ratio', 'se', 'predicted_capture'])
        for key in SPECTRA:
            psd = SPECTRA[key]['psd']
            capture = {k: helpers.captured_fraction(
                psd, PUPIL_D_EFF_M, 1.0 / (k * SIDE_REF_M))
                for k in OVERSIZE_K}
            targets = {'z': z_tilt_target(psd), 'g': g_tilt_target(psd)}
            for (arm, k), result in sorted(sweep[key].items()):
                for metric in ('z', 'g'):
                    var = result[f'{metric}_var']
                    ratio = var / targets[metric]
                    se = ratio * np.sqrt(1.0 / result['trials'])
                    writer.writerow([key, arm, k, result['trials'], metric,
                                     f'{var:.6e}', f'{ratio:.6f}',
                                     f'{se:.6f}', f'{capture[k]:.6f}'])


def write_dphi_csv(sweep):
    '''Write the structure function of every arm.'''
    with open(DPHI_CSV, 'w', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(['spectrum', 'arm', 'K', 'r_m', 'd_meas', 'd_theory',
                         'ratio'])
        for key in SPECTRA:
            theory = dphi_theory(key)
            for (arm, k), result in sorted(sweep[key].items()):
                ratio = result['d_measured'] / theory
                for i, r_m in enumerate(R_BINS_M):
                    writer.writerow([key, arm, k, f'{r_m:.5f}',
                                     f'{result["d_measured"][i]:.6f}',
                                     f'{theory[i]:.6f}', f'{ratio[i]:.6f}'])


def write_s27_csv(ratios):
    '''Write the S-27 settlement table.'''
    with open(S27_CSV, 'w', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(['generator', 'estimator', 'r_over_r0', 'ratio'])
        for (name, estimator), row in ratios.items():
            for probe, value in zip(S27_PROBES_R0, row):
                writer.writerow([name, estimator, probe, f'{value:.6f}'])


def write_quirk_csv(quirk):
    '''Write the aotools seed-quirk measurement.'''
    with open(QUIRK_CSV, 'w', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(['quantity', 'value', 'standard_error'])
        writer.writerow(['samples', quirk['count'], ''])
        writer.writerow(['correlation_lf_hf', f'{quirk["correlation"]:.6f}',
                         f'{quirk["correlation_se"]:.6f}'])
        writer.writerow(['var_hf_rad2', f'{quirk["var_hf"]:.6e}', ''])
        writer.writerow(['var_lf_rad2', f'{quirk["var_lf"]:.6e}', ''])
        writer.writerow(['var_correlated_rad2',
                         f'{quirk["var_correlated"]:.6e}',
                         f'{quirk["var_se"]:.6e}'])
        writer.writerow(['var_decorrelated_rad2',
                         f'{quirk["var_decorrelated"]:.6e}',
                         f'{quirk["var_se"]:.6e}'])


# ---------------------------------------------------------------------------
# The pass bands
# ---------------------------------------------------------------------------

BANDS = []


def band(label, value, low, high, note=''):
    '''Test one pass band. Print the result and keep it.

    A failed band does NOT stop the script. A physics surprise is a result of
    this study, so the caller must see every band.
    '''
    try:
        assert low <= value <= high, (value, low, high)
        state = 'PASS'
    except AssertionError:
        state = 'FAIL'
    text = (f'  {state}  {label:<52}{value:>9.4f}  band '
            f'[{low:.4f}, {high:.4f}]{"  " + note if note else ""}')
    BANDS.append(text)
    print(text)
    return state == 'PASS'


# ---------------------------------------------------------------------------

def main():
    t_start = time.time()

    print('=' * 78)
    print('validation.screens.oversize_crop')
    print('the oversized-and-cropped phase screen against the subharmonics')
    print('=' * 78)
    print(f'  wavelength           {WAVELENGTH_M * 1e9:11.1f} nm')
    print(f'  r0                   {R0_M * 1e2:11.1f} cm')
    print(f'  pitch                {DX_M * 1e2:11.1f} cm   (FIXED in every '
          f'arm)')
    print(f'  crop grid            {N_REF:11d} px, side {SIDE_REF_M:.2f} m')
    print(f'  oversize factors     {str(OVERSIZE_K):>11}   sides '
          f'{[round(k * SIDE_REF_M, 2) for k in OVERSIZE_K]} m')
    print(f'  tilt pupil           {PUPIL_D_M:11.2f} m   (area equivalent '
          f'{PUPIL_D_EFF_M:.4f} m)')
    print(f'  structure window     {2.0 * SF_MASK_RADIUS_M:11.2f} m diameter '
          f'  (a guard band for the circular estimate)')
    print(f'  master seed          {MASTER_SEED:11d}')

    # ---- 1. the oversize sweep ----
    sweep = {}
    for key in SPECTRA:
        t0 = time.time()
        sweep[key] = run_sweep(key)
        print(f'  {key} sweep took {time.time() - t0:.1f} s')

    print('')
    print('1. THE OVERSIZE SWEEP')
    for key in SPECTRA:
        psd = SPECTRA[key]['psd']
        z_target = z_tilt_target(psd)
        g_target = g_tilt_target(psd)
        print('')
        print(f'  {SPECTRA[key]["label"]}')
        print(f'  analytic per-axis Z-tilt variance {z_target:.5e} rad^2, '
              f'G-tilt {g_target:.5e} rad^2')
        print(f'  {"arm":<16}{"K":>3}{"M":>5}{"Z var, rad^2":>15}'
              f'{"ratio":>9}{"+-2SE":>9}{"capture":>10}{"G ratio":>10}')
        for (arm, k), result in sorted(sweep[key].items(),
                                       key=lambda kv: (kv[0][1], kv[0][0])):
            ratio = result['z_var'] / z_target
            se = ratio * np.sqrt(1.0 / result['trials'])
            capture = helpers.captured_fraction(psd, PUPIL_D_EFF_M,
                                                1.0 / (k * SIDE_REF_M))
            print(f'  {arm:<16}{k:>3}{result["trials"]:>5}'
                  f'{result["z_var"]:>15.4e}{ratio:>9.4f}'
                  f'{BAND_SE * se:>9.4f}{capture:>10.4f}'
                  f'{result["g_var"] / g_target:>10.4f}')
        deep = helpers.captured_fraction(
            psd, PUPIL_D_EFF_M, 1.0 / (3.0 ** N_SUBHARMONIC * SIDE_REF_M))
        print('  The `capture` column is the SHARP-CUTOFF model at the '
              'grid')
        print('  fundamental 1/(K * 5.12 m). It is the right model for a '
              'plain')
        print('  Fourier screen, and a LOWER bound for a subharmonic '
              'screen.')
        print(f'  The subharmonics reach 1/(3^3 * 5.12 m), where the same '
              f'model gives {deep:.4f}.')

    print('')
    print('   THE STRUCTURE FUNCTION, as a ratio to the analytic law')
    print('   The theory goes through the SAME annuli as the measurement, so')
    print('   the ratio carries no binning bias. See `bin_theory`.')
    for key in SPECTRA:
        theory = dphi_theory(key)
        print('')
        print(f'  {SPECTRA[key]["label"]}')
        header = f'  {"r/r0":>7}{"D theory":>11}'
        arms = sorted(sweep[key], key=lambda a: (a[1], a[0]))
        for arm, k in arms:
            name = arm.replace('book ', '').replace('aotools ', 'ao ')
            header += f'{name + ("" if k == 1 else f" K{k}"):>12}'
        print(header)
        for i, r_m in enumerate(R_BINS_M):
            line = f'  {r_m / R0_M:>7.2f}{theory[i]:>11.3f}'
            for arm, k in arms:
                value = sweep[key][(arm, k)]['d_measured'][i] / theory[i]
                line += f'{value:>12.4f}'
            print(line)

    # ---- 2. the S-27 settlement ----
    t0 = time.time()
    ratios, s27_theory, tilts = run_s27()
    print('')
    print('2. THE S-27 SETTLEMENT, Kolmogorov, one shared ensemble')
    print(f'   {S27_TRIALS} trials, N = {N_REF}, r0 = '
          f'{R0_M / DX_M:.0f} px, both estimators ({time.time() - t0:.1f} s)')
    names = ['book ft', 'book ft_sh', 'aotools ft', 'aotools ft_sh']
    for estimator in ('fft', 'direct'):
        print('')
        label = ('the FFT estimator, ensemble_dphi, 2.56 m window, '
                 'annulus-binned theory' if estimator == 'fft'
                 else 'the direct estimator, d_phi_direct, whole grid, exact '
                      'theory')
        print(f'  {label}')
        print(f'  {"r/r0":>7}{"D theory":>11}'
              + ''.join(f'{n:>15}' for n in names))
        for j, probe in enumerate(S27_PROBES_R0):
            line = f'  {probe:>7.1f}{s27_theory[estimator][j]:>11.3f}'
            for name in names:
                line += f'{ratios[(name, estimator)][j]:>15.4f}'
            print(line)
        gap = 100.0 * (ratios[('aotools ft_sh', estimator)]
                       / ratios[('book ft_sh', estimator)] - 1.0)
        print(f'  aotools ft_sh against book ft_sh: {gap.min():+.1f} to '
              f'{gap.max():+.1f} percent')

    gap_fft = 100.0 * (ratios[('aotools ft_sh', 'fft')]
                       / ratios[('book ft_sh', 'fft')] - 1.0)
    gap_direct = 100.0 * (ratios[('aotools ft_sh', 'direct')]
                          / ratios[('book ft_sh', 'direct')] - 1.0)
    print('')
    print('  THE VERDICT. The repository states TWO different numbers for the')
    print('  same comparison:')
    print('    (A) olb/waveoptics/turbulence/screens.py and')
    print('        examples/schmidt/screens_and_turbulence.py: aotools '
          'reads 1')
    print('        to 3 percent ABOVE the book form (FFT estimator, N = 512).')
    print('    (B) docs/schmidt-crosscheck.md, the forward map and gap S-27:')
    print('        aotools reads 5 to 12 percent BELOW the book form (direct')
    print('        estimator, N = 256).')
    print(f'  This run, at N = {N_REF} with {S27_TRIALS} trials, measures:')
    print(f'    FFT estimator     {gap_fft.min():+.1f} to '
          f'{gap_fft.max():+.1f} percent')
    print(f'    direct estimator  {gap_direct.min():+.1f} to '
          f'{gap_direct.max():+.1f} percent')
    if gap_fft.mean() * gap_direct.mean() > 0.0:
        supported = 'A' if gap_fft.mean() > 0.0 else 'B'
        print('  The two estimators AGREE on the sign, so the estimator is '
              'not the')
        print(f'  cause. The measurement supports doc {supported}, and the '
              f'other doc number')
        print(f'  is a grid-size effect: doc B ran at N = 256, and this run '
              f'is at N = {N_REF}.')
    else:
        print('  The two estimators DISAGREE on the sign at the same N. So '
              'the')
        print('  estimator IS the cause, and both doc numbers are right for '
              'their own')
        print('  estimator. Doc A is the FFT route and doc B is the direct '
              'route.')

    # ---- 3. the aotools seed quirk ----
    quirk = run_seed_quirk(tilts)
    print('')
    print('3. THE aotools SEED QUIRK')
    print('   ft_sh_phase_screen gives its subharmonic generator the SAME '
          'integer')
    print('   seed that it passes to ft_phase_screen. So the 9 low-frequency '
          'real')
    print('   draws repeat the first 9 high-frequency real draws.')
    print(f'   {"samples (2 axes x trials)":<40}{quirk["count"]:>14d}')
    print(f'   {"correlation of Z-tilt LF against HF":<40}'
          f'{quirk["correlation"]:>14.4f}')
    print(f'   {"its standard error, 1/sqrt(M)":<40}'
          f'{quirk["correlation_se"]:>14.4f}')
    print(f'   {"Z-tilt variance, HF alone, rad^2":<40}'
          f'{quirk["var_hf"]:>14.4e}')
    print(f'   {"Z-tilt variance, LF alone, rad^2":<40}'
          f'{quirk["var_lf"]:>14.4e}')
    print(f'   {"Z-tilt variance, HF + LF same seed":<40}'
          f'{quirk["var_correlated"]:>14.4e}')
    print(f'   {"Z-tilt variance, HF + LF shifted":<40}'
          f'{quirk["var_decorrelated"]:>14.4e}')
    print(f'   {"2 SE of a variance":<40}'
          f'{2.0 * quirk["var_se"]:>14.4e}')
    quirk_delta = abs(quirk['var_correlated'] - quirk['var_decorrelated'])
    quirk_biased = quirk_delta > 2.0 * quirk['var_se']
    print(f'   VERDICT: the shared seed '
          f'{"BIASES" if quirk_biased else "does NOT bias"} the Z-tilt '
          f'variance beyond 2 SE '
          f'(the gap is {quirk_delta:.3e}, and 2 SE is '
          f'{2.0 * quirk["var_se"]:.3e}).')

    # ---- 4. the pass bands ----
    print('')
    print('4. THE PASS BANDS')

    # Band 1. The plain Fourier screen at K = 1 must match the sharp-cutoff
    # model, because that model IS the plain Fourier screen.
    result = sweep['kolmogorov'][(ARM_FT, 1)]
    target = z_tilt_target(SPECTRA['kolmogorov']['psd'])
    ratio = result['z_var'] / target
    predicted = helpers.captured_fraction(SPECTRA['kolmogorov']['psd'],
                                          PUPIL_D_EFF_M, 1.0 / SIDE_REF_M)
    se = ratio * np.sqrt(1.0 / result['trials'])
    band(f'kolmogorov book ft K=1 against capture {predicted:.4f}', ratio,
         predicted - BAND_SE * se, predicted + BAND_SE * se)

    # Band 2. Every oversize point against its own sharp-cutoff prediction.
    # The K = 8 arms carry the fewest trials, so they take 3 SE. That LOOSENED
    # band prints in the note column.
    for key in SPECTRA:
        psd = SPECTRA[key]['psd']
        target = z_tilt_target(psd)
        for k in OVERSIZE_K[1:]:
            result = sweep[key][(ARM_CROP, k)]
            ratio = result['z_var'] / target
            predicted = helpers.captured_fraction(psd, PUPIL_D_EFF_M,
                                                  1.0 / (k * SIDE_REF_M))
            se = ratio * np.sqrt(1.0 / result['trials'])
            width = BAND_SE_K8 if k == 8 else BAND_SE
            note = ('LOOSENED to 3 SE, M is small at K = 8' if k == 8 else '')
            band(f'{key} crop K={k} against capture {predicted:.4f}', ratio,
                 predicted - width * se, predicted + width * se, note)

    # Band 3. The von Karman outer scale is 25 m, and the K = 8 side is 40.96
    # m. So the screen holds nearly the whole tilt band.
    result = sweep['von karman'][(ARM_CROP, 8)]
    target = z_tilt_target(SPECTRA['von karman']['psd'])
    band('von karman crop K=8 Z-tilt ratio', result['z_var'] / target,
         VK_K8_TILT_BAND[0], VK_K8_TILT_BAND[1])

    # Band 4. The same arm must hold the structure function too.
    theory = dphi_theory('von karman')
    dphi_ratio = result['d_measured'] / theory
    keep = (R_BINS_M / R0_M >= 1.0) & (R_BINS_M / R0_M <= DPHI_MAX_R0)
    band('von karman crop K=8 D(r)/theory, lowest bin',
         float(np.min(dphi_ratio[keep])), VK_K8_DPHI_BAND[0],
         VK_K8_DPHI_BAND[1])
    band('von karman crop K=8 D(r)/theory, highest bin',
         float(np.max(dphi_ratio[keep])), VK_K8_DPHI_BAND[0],
         VK_K8_DPHI_BAND[1])

    # ---- 5. the output files ----
    write_tilt_csv(sweep)
    write_dphi_csv(sweep)
    write_s27_csv(ratios)
    write_quirk_csv(quirk)
    draw_tilt(sweep)
    draw_dphi(sweep)

    print('')
    for path in (TILT_CSV, DPHI_CSV, S27_CSV, QUIRK_CSV):
        print(f'  file saved: {path}')
    print(f'  figure saved: {TILT_PNG}')
    print('    Caption: the measured Z-tilt variance of each arm, divided by '
          'the Noll')
    print('    analytic value, against the oversize factor K. The black curve '
          'is the')
    print('    sharp-cutoff model, the orange band is the subharmonic arm, '
          'and the red')
    print('    square is the production aotools screen.')
    print(f'  figure saved: {DPHI_PNG}')
    print('    Caption: the ensemble structure function of each arm, divided '
          'by the')
    print('    analytic law, against the separation in units of r0.')

    failed = [line for line in BANDS if line.strip().startswith('FAIL')]
    print('')
    print(f'  {len(BANDS) - len(failed)} of {len(BANDS)} pass bands hold')
    for line in failed:
        print(line)
    print('')
    print(f'(elapsed {time.time() - t_start:.1f} s)')


if __name__ == '__main__':
    main()
