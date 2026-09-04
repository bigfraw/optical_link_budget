"""Measure the FAST point-ahead residual against the Stone (1994) analytic law
(backlog 1-5).

THE QUESTION. An uplink terminal senses the turbulence on a downlink beacon and
it applies the conjugate to the uplink beam. The satellite moves, so the uplink
goes to a different angle than the beacon came from. That point-ahead angle
decorrelates the correction. olb holds TWO models of the residual phase variance
of that decorrelation:

  - FIDELITY 1. FAST (fast-aosim) builds the residual phase power spectrum
    G_AO * Phi_n and it integrates it over the corrected (low-order) band. olb
    calls it through olb.models.fast.uplink_fast_term. It is the model of
    record.
  - FIDELITY 0. The Stone (1994) finite-aperture modal law
    (olb.turbulence.anisoplanatism.anisoplanatic_phase_variance) sums the
    decorrelation residual of the Zernike radial orders 2..max_order. olb calls
    it through olb.links.uplink.uplink_point_ahead_term.

The two must give the same number, because they compute the same quantity. This
script measures whether they do, and it names every reason they do not.

THE TEST. The script drives FAST with the servo OFF (TLOOP = 0, TEXP = 0, no
wind), so the PAOLA filter reduces to the pure anisoplanatic kernel

    G_aniso(kappa) = 2 - 2 cos(delta_r . kappa),   delta_r_i = theta h_i,

which is the same physics that Stone integrates. Sources:
- J. Stone, P. H. Hu, S. P. Mills and S. Ma, J. Opt. Soc. Am. A 11(1), 347
  (1994), DOI 10.1364/JOSAA.11.000347. Eqs. (29), (36) and (A11).
- O. J. D. Farley and others, Opt. Express 30(13), 23050 (2022),
  DOI 10.1364/OE.458659. The FAST method.
- R. J. Noll, J. Opt. Soc. Am. 66(3), 207 (1976), DOI 10.1364/JOSA.66.000207.
  The fitting-error variance of the uncorrected orders.
- Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196, Ch. 3. The
  Kolmogorov and von Karman refractive-index spectra.
- T. S. Ross, Appl. Opt. 48(10), 1812 (2009), DOI 10.1364/AO.48.001812. The
  extended Marechal relation that turns rad^2 into dB.

THE SOFT-MASK SPLIT. FAST corrects modally. Its low-frequency mask is then the
SOFT Zernike Fourier filter sum_j |Zhat_j(kappa)|^2 (piston included, clipped at
1), not a hard cut. The shipped attribute sim.aniso_servo_error integrates
G_ao * mask = aniso * mask^2 + mask (1 - mask), so it mixes a part of the
UNCORRECTED band into the anisoplanatic number. The script rebuilds the kernel
from the sim attributes and it gives a CLEAN split:

    aniso_corr = INT 2 pi k^2 Phi_n (2 - 2 cos(delta_r . kappa)) mask d^2kappa
    fit_corr   = INT 2 pi k^2 Phi_n (1 - mask) d^2kappa

The two add up to sim.phs_var exactly, and the script asserts that closure.

THE LOW-FREQUENCY TRUNCATION (measured here, 2026-09-02). The Kolmogorov
anisoplanatic integrand goes as kappa^(-2/3) at low frequency, so its integral
converges only as kappa^(1/3). A FAST grid of N pixels and pitch dx carries no
frequency below df = 2 pi / (N dx), so it MISSES a real part of the variance.
The script measures that part. It integrates the same physics three ways:

  1. on the FAST grid (what FAST reports);
  2. with an independent polar quadrature over the support that the grid holds
     (the square domain up to kappa_max = pi/dx, without the direct-current
     cell below df/2);
  3. with the same quadrature over the WHOLE plane, which is what Stone
     integrates. The ratio of 2 to 3 is the truncation.

THE GATE NEEDS AN OUTER SCALE. In the Kolmogorov limit no grid converges on
this integral, so a 1 % gate there has no meaning. The GATED legs (stage A0 and
the convergence stage) carry a finite outer scale L0 = L0_GATE. That puts the
low-frequency knee at k0 = 2 pi / L0 inside the grid, so the FAST value and the
quadrature must then agree. Every other stage keeps the Kolmogorov limit, and
the script reports the truncation instead of gating it.

THE MODE SETS MUST MATCH. The FAST modal mask sums the Noll modes 1..ZMAX, so
it KEEPS the piston and the two tilts. The Stone set that holds the same modes
is remove='none' over the band 0..max_order. That pair (the column Q/S_n) is
the physics test. The production Term uplink_point_ahead_term removes the
piston and the tilt, because a separate tracking loop points the beam, so the
pair FAST against Stone piston+tilt (the column F/S_pt) is NOT mode matched.

THE STAGES.
  A0  A single turbulence layer at the zenith. It gives the cleanest test of the
      kernel, the mask, and the quadrature.
  A1  The HV5/7 profile at the production point (1.5 m, 60 deg, the 600 km orbit
      point-ahead angle). It prints the full attribution table.
  B   Three sweeps: the point-ahead angle, the corrected order, and the
      elevation.
  Convergence  The same case on three grids.
  C   The full Term comparison in dB: the fidelity-1 FAST Monte Carlo against
      the fidelity-0 analytic pair, with the attribution ladder between them.

VALIDATION ONLY. The script reads the production layer and it changes no olb
module.

Run it from the repository root:

    python -m validation.fast_stone_pointahead.fast_stone_pointahead
    python -m validation.fast_stone_pointahead.fast_stone_pointahead --full
"""

import argparse
import json
import logging
import os
import subprocess
import sys
import time
import warnings

import numpy as np

# numpy 2.4 removed the np.trapz alias, and olb still calls it. Restore the name
# here, in the validation script only. Do not change production code.
if not hasattr(np, "trapz"):
    np.trapz = np.trapezoid

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from olb.geometry import CircularOrbit
from olb.links.uplink import uplink_fitting_term, uplink_point_ahead_term
from olb.models.fast import (ARCSEC_PER_RAD, _cn2_layers, _load_fast,
                             uplink_fast_term)
from olb.scenario import Channel, DownlinkBeacon, Site, SpaceScenario
from olb.terminal import AO, Terminal, TipTilt, Transmitter
from olb.turbulence.anisoplanatism import (anisoplanatic_phase_variance,
                                           max_radial_order)
from olb.turbulence.ao import (apply_compensation,
                               plane_wave_fried_parameter_profile)
from olb.turbulence.profiles import DEFAULT_HS, get_c2n

HERE = os.path.dirname(os.path.abspath(__file__))
_LN10 = np.log(10.0)

# The link. A 1.5 m ground aperture, a 0.3 m launch waist, and a 600 km orbit.
LAM = 1.55e-6
D_GROUND = 1.5
W0_M = 0.3
D_SAT = 0.3
ALT_M = 600e3
CN2_GROUND = 1.7e-14          # HV5/7 ground value [m^-2/3]
WIND_RMS = 21.0               # HV5/7 high-altitude wind [m/s]

# The servo-off wavefront-sensor pitch. It puts the square WFS box edge at
# pi / DSUBAP = 314 rad/m, far above the Zernike filter cutoff, so the modal
# filter alone sets the corrected band.
DSUBAP_OFF = 0.01

# The grids. NPXLS and DX are always explicit: NPXLS='auto' divides by TLOOP,
# and the servo-off runs set TLOOP = 0.
GRID_QUICK = (512, 0.02)
GRID_FULL = (1024, 0.01)

# Stage A0 keeps the coarse grid in both modes. Its profile holds 81 layers, and
# FAST holds one spectrum for each layer, so a 1024 grid there needs several
# gigabytes. A0 tests the kernel, the mask and the quadrature, and the coarse
# grid tests them fully.
GRID_A0 = GRID_QUICK

# The outer scale of the GATED legs [m]. See THE LOW-FREQUENCY TRUNCATION: the
# Kolmogorov (L0 = infinite) integral has no low-frequency floor, so a grid
# cannot converge on it and a 1 % gate has no meaning. A 5 m outer scale puts
# the low-frequency knee at k0 = 2 pi / L0 = 1.26 rad/m, well inside the grid,
# so the same integral becomes a fair numerics test. Source of the von Karman
# spectrum: Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196, Ch. 3.
L0_GATE = 5.0

# The lowest frequency of the whole-plane quadrature [rad/m]. It stands for an
# outer scale of 6000 km, so it is the Kolmogorov (L0 = infinite) limit. The
# mass below it is about 0.5 % of the integral.
KAPPA_MIN_FULL = 1e-6

# The highest frequency of the polar quadrature [rad/m]. It sits above the
# square corner of the finest grid (kappa_max sqrt(2) = 444 rad/m).
KAPPA_MAX_QUAD = 500.0

# The polar quadrature grid.
N_KAPPA = 1400
N_PHI = 512

# The nominal corrected-mode count of the sweeps (Zernike order 9 is complete at
# 55 Noll modes).
ZMAX_NOM = 55

# The gate limits.
GATE_QUADRATURE = 0.01        # A0: grid against the independent quadrature
GATE_CLOSURE = 0.005          # every sim: aniso_corr + fit_corr against phs_var
GATE_CONVERGENCE = 0.01       # the three-grid movement
GATE_ZERO_RAD2 = 1e-6         # the theta = 0 residual

# The verdict rule of the physics comparison. It reads the ratio of the FAST
# corrected residual to the Stone piston-and-tilt-removed residual.
VERDICT_MATCH = 0.10
VERDICT_DIFFERENCE = 0.40


def _marechal_db(sigma2):
    """Give the extended Marechal loss [dB] of a residual phase variance.

    formula:
        loss_dB = (10 / ln 10) * sigma^2
    Source: T. S. Ross, Appl. Opt. 48(10), 1812 (2009),
    DOI 10.1364/AO.48.001812. It is the same relation that
    olb.links.uplink uses.
    """
    return float(10.0 / _LN10 * sigma2)


def _load():
    """Import fast-aosim and quiet its logger. Return (fast, funcs, spectra)."""
    fast = _load_fast()
    from fast import ao_power_spectra, funcs
    logging.getLogger("fast").setLevel(logging.ERROR)
    return fast, funcs, ao_power_spectra


def _servo_off_params(elev_deg, hs, cn2, zmax, theta_arcsec, grid, L0=np.inf):
    """Build the FAST parameter dict of a SERVO-OFF run.

    TLOOP = 0, TEXP = 0 and a zero wind field remove the servo lag and the
    exposure blur from the PAOLA filter (fast.ao_power_spectra.G_AO_PAOLA), so
    the filter holds the anisoplanatic kernel 2 - 2 cos(delta_r . kappa) alone.
    NPXLS and DX are explicit, because the automatic grid divides by TLOOP.
    """
    npxls, dx = grid
    return dict(
        WVL=LAM, D_GROUND=D_GROUND, OBSC_GROUND=0.0, W0=W0_M,
        PROP_DIR="up", SMF=False, D_SAT=D_SAT, OBSC_SAT=0.0, H_SAT=ALT_M,
        ZENITH_ANGLE=90.0 - elev_deg, H_TURB=hs, CN2_TURB=_cn2_layers(cn2, hs),
        WIND_SPD=np.zeros_like(hs), WIND_DIR=np.zeros_like(hs),
        DTHETA=[float(theta_arcsec), 0.0], AO_MODE="AO", MODAL=True,
        ZMAX=int(zmax), DSUBAP=DSUBAP_OFF, TLOOP=0.0, TEXP=0.0, ALIAS=False,
        NOISE=0, NPXLS=int(npxls), DX=float(dx), SUBHARM=False, L0=float(L0),
        l0=1e-6, NITER=2, NCHUNKS=1, LOGLEVEL="ERROR")


def _build(fast, params):
    """Construct a FAST simulation. It computes the spectra, and it does not run
    the Monte Carlo."""
    with warnings.catch_warnings():
        # A zero wind field makes aotools divide by zero in the coherence time.
        # The script does not read that value.
        warnings.simplefilter("ignore")
        return fast.Fast(params)


def _kernel(sim):
    """Rebuild the pure anisoplanatic kernel 2 - 2 cos(delta_r . kappa).

    The kernel is one array for each turbulence layer. It matches the
    broadcasting of fast.ao_power_spectra.G_AO_PAOLA, so
    G_ao = kernel * mask + (1 - mask) holds to machine precision. The script
    asserts that identity.
    """
    h = sim.h
    fx, fy = sim.freq.main.fx, sim.freq.main.fy
    delta_r = (np.tile(sim.dtheta, (len(h), 1)).T / ARCSEC_PER_RAD * h).T
    fx_tile = np.tile(fx, (len(h), *[1] * fx.ndim))
    fy_tile = np.tile(fy, (len(h), *[1] * fy.ndim))
    dot = (fx_tile.T * delta_r[:, 0] + fy_tile.T * delta_r[:, 1]).T
    return 2.0 - 2.0 * np.cos(dot)


def split_errors(sim, funcs):
    """Give the CLEAN anisoplanatic and fitting variances of a FAST sim [rad^2].

    sim.aniso_servo_error integrates G_ao * mask. For a SOFT modal mask that
    mixes the two bands, because G_ao * mask = aniso * mask^2 + mask (1 - mask).
    This function integrates the kernel against the mask ONE time, and the raw
    turbulence spectrum against (1 - mask), so the two bands stay apart.

    Parameters:
        sim : fast.Fast
            A constructed simulation.
        funcs : module
            fast.funcs, for the two shipped integrators.

    Returns:
        dict
            aniso_corr, fit_corr, closure_corr and closure_raw.
    """
    kern = _kernel(sim)
    scale = 2.0 * np.pi * sim.k ** 2
    aniso = funcs.integrate_powerspectrum(
        funcs.integrate_path(sim.turb_powerspec * kern, sim.h, layer=True)
        * sim.lf_mask * scale, sim.freq.main.f)
    fit = funcs.integrate_powerspectrum(
        funcs.integrate_path(sim.turb_powerspec, sim.h, layer=True)
        * sim.hf_mask * scale, sim.freq.main.f)
    raw_sum = (sim.aniso_servo_error + sim.fitting_error + sim.alias_error
               + sim.noise_error)
    return {
        "aniso_corr": float(aniso),
        "fit_corr": float(fit),
        "closure_corr": float((aniso + fit - sim.phs_var) / sim.phs_var),
        "closure_raw": float((raw_sum - sim.phs_var) / sim.phs_var),
    }


# The polar grid of the independent quadrature. It is fixed, so the modal mask
# of one ZMAX is built one time and kept.
_KAPPA = np.logspace(np.log10(KAPPA_MIN_FULL), np.log10(KAPPA_MAX_QUAD),
                     N_KAPPA)
_PHI = np.linspace(0.0, 2.0 * np.pi, N_PHI, endpoint=False)
_KA, _PH = np.meshgrid(_KAPPA, _PHI, indexing="ij")
_KX, _KY = _KA * np.cos(_PH), _KA * np.sin(_PH)
_MASKS = {}


def _modal_mask(spectra, zmax):
    """Give the FAST low-frequency mask on the polar quadrature grid.

    The mask is the soft Zernike Fourier filter sum_j |Zhat_j(kappa)|^2 of
    fast.ao_power_spectra.zernike_squared_filter (piston included, clipped at
    1), times the square wavefront-sensor box |kappa_x|, |kappa_y| <=
    pi / DSUBAP_OFF. It is the same mask that fast.ao_power_spectra.mask_lf
    builds on the Cartesian grid.
    """
    key = int(zmax)
    if key not in _MASKS:
        mask = np.clip(spectra.zernike_squared_filter(
            _KA, _KX, _KY, D_GROUND, key).real, 0.0, 1.0)
        # zernike_squared_filter forces its central element to 1, a
        # direct-current convention of the FAST Cartesian grid. This grid is
        # polar, so repair that one element from its two neighbours.
        i, j = _KA.shape[0] // 2, _KA.shape[1] // 2
        mask[i, j] = 0.5 * (mask[i, j - 1] + mask[i, j + 1])
        mask = mask * ((np.abs(_KX) <= np.pi / DSUBAP_OFF)
                       & (np.abs(_KY) <= np.pi / DSUBAP_OFF))
        _MASKS[key] = mask
    return _MASKS[key]


def _quadrature(spectra, cn2_dh, h, theta_rad, zmax, kmax=None, L0=np.inf,
                dc=0.0):
    """Integrate the anisoplanatic variance on an independent polar grid.

    formula:
        sigma^2 = INT 2 pi k^2 Phi_n(kappa) (2 - 2 cos(delta_r . kappa))
                  M(kappa) d^2 kappa
        Phi_n(kappa) = 0.033 SUM_l Cn2_l dh_l (kappa^2 + k0^2)^(-11/6),
        k0 = 2 pi / L0
    That is the von Karman refractive-index spectrum that FAST integrates
    (fast.funcs.turb_powerspectrum_vonKarman with C = 2 pi). An infinite L0
    gives the Kolmogorov limit. Source: Andrews and Phillips, 2nd ed. (2005),
    DOI 10.1117/3.626196, Ch. 3. The inner scale is 1 um, so its exponential
    cutoff sits at 6e6 rad/m and this grid does not reach it.

    This routine shares NO code with the FAST integrators. It is the
    independent leg of the A0 gate.

    Parameters:
        cn2_dh : numpy.ndarray
            The integrated Cn2 of each layer on the LINE OF SIGHT [m^(1/3)].
        h : numpy.ndarray
            The layer heights on the line of sight [m].
        theta_rad : float
            The point-ahead angle [rad], along x.
        zmax : int
            The corrected Zernike mode count.
        kmax : float or None
            The half-width of the square support [rad/m]. None takes the whole
            plane. Use pi/dx to match a FAST grid.
        L0 : float
            The outer scale [m].
        dc : float
            Leave out the central cell |kappa_x|, |kappa_y| <= dc [rad/m]. Use
            df/2 to stand for the direct-current cell of a FAST grid, which
            carries no frequency below df.

    Returns:
        float
            sigma^2 [rad^2].
    """
    k = 2.0 * np.pi / LAM
    mask = _modal_mask(spectra, zmax)
    if kmax is not None:
        mask = mask * ((np.abs(_KX) <= kmax) & (np.abs(_KY) <= kmax))
    if dc > 0.0:
        mask = mask * ~((np.abs(_KX) <= dc) & (np.abs(_KY) <= dc))

    # The layer sum of the kernel, SUM_l c_l (2 - 2 cos(kappa_x theta h_l)).
    # The loop keeps the memory at one grid, because a layer axis on this grid
    # is too large to hold.
    base = _KX * theta_rad
    kern = np.zeros_like(_KA)
    for c, height in zip(np.asarray(cn2_dh, float), np.asarray(h, float)):
        kern += c * (2.0 - 2.0 * np.cos(base * height))

    k0 = 2.0 * np.pi / L0
    phi_n = 0.033 * (_KA ** 2 + k0 ** 2) ** (-11.0 / 6.0)
    integrand = 2.0 * np.pi * k ** 2 * phi_n * kern * mask * _KA
    return float(np.trapz(np.trapz(integrand, _PHI, axis=1), _KAPPA))


_STONE_CACHE = {}


def _stone(D, theta_rad, hs, cn2, elev_deg, max_order):
    """Give the Stone residual for the three mode sets [rad^2].

    The Bessel quadrature of Eq. (36) runs one time for each height, so the
    call is expensive. Several rows ask for the SAME numbers (the outer scale
    and the grid do not change an analytic value), so the results stay in a
    cache.
    """
    key = (float(D), float(theta_rad), float(elev_deg), int(max_order),
           len(hs), float(np.sum(cn2)), float(np.sum(hs)))
    if key not in _STONE_CACHE:
        out = {}
        for remove in ("none", "piston", "piston_tilt"):
            out[remove] = float(anisoplanatic_phase_variance(
                D, theta_rad, hs, cn2, LAM, remove=remove,
                max_order=max_order, elevation_deg=elev_deg))
        _STONE_CACHE[key] = out
    return _STONE_CACHE[key]


def _noll(hs, cn2, elev_deg, D, n_modes):
    """Give the Noll fitting variance of the uncorrected orders [rad^2].

    formula:
        sigma^2 = 0.2944 n^(-sqrt(3)/2) (D / r0)^(5/3)
    Source: R. J. Noll, J. Opt. Soc. Am. 66(3), 207 (1976),
    DOI 10.1364/JOSA.66.000207. The call goes through
    olb.turbulence.ao.apply_compensation, so the script reads the production
    coefficient.
    """
    r0 = plane_wave_fried_parameter_profile(cn2, hs, LAM, elev_deg)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return float(apply_compensation([AO(int(n_modes))], D, r0).variance)


def _row(fast, funcs, spectra, stage, label, elev_deg, zmax, theta_arcsec, hs,
         cn2, grid, max_order=None, L0=np.inf):
    """Measure ONE case: one FAST sim, the clean split, and the Stone values.

    The quadrature runs two times: over the SAME square support as the FAST
    grid (the gate leg) and over the whole plane (the truncation leg).

    Returns:
        dict
            The JSON row of the case.
    """
    t0 = time.time()
    npxls, dx = grid
    params = _servo_off_params(elev_deg, hs, cn2, zmax, theta_arcsec, grid, L0)
    sim = _build(fast, params)

    # The kernel identity. It proves that the rebuilt kernel is the one the
    # PAOLA filter used.
    kern_err = float(np.abs(
        sim.G_ao - (_kernel(sim) * sim.lf_mask + (1.0 - sim.lf_mask))).max())

    split = split_errors(sim, funcs)
    assert abs(split["closure_corr"]) < GATE_CLOSURE, (label,
                                                       split["closure_corr"])

    theta_rad = theta_arcsec / ARCSEC_PER_RAD
    order = max_radial_order(int(zmax)) if max_order is None else int(max_order)
    stone = _stone(D_GROUND, theta_rad, hs, cn2, elev_deg, order)

    # The quadrature legs. Both use the ZENITH-CORRECTED layers of the sim, so
    # the elevation enters exactly as FAST applies it.
    cn2_dh = _cn2_layers(cn2, hs) * sim.zenith_correction
    kmax = np.pi / dx
    # The SUPPORT leg stands for what the FAST grid can hold: the square domain
    # up to kappa_max, without the direct-current cell below df/2.
    q_support = _quadrature(spectra, cn2_dh, sim.h, theta_rad, zmax, kmax, L0,
                            dc=0.5 * float(sim.freq.main.df))
    q_full = _quadrature(spectra, cn2_dh, sim.h, theta_rad, zmax, None, L0)

    noll = _noll(hs, cn2, elev_deg, D_GROUND, int(zmax))
    stone_pt = stone["piston_tilt"]
    return {
        "stage": stage, "label": label, "elev": float(elev_deg),
        "zmax": int(zmax), "max_order": int(order),
        "theta_arcsec": float(theta_arcsec),
        "npxls": int(npxls), "dx": float(dx),
        "df": float(sim.freq.main.df), "L0_m": float(L0),
        "fast_raw": float(sim.aniso_servo_error),
        "fast_corr": split["aniso_corr"],
        "quad_support": q_support,
        "quad_full": q_full,
        "quad_rel": float((split["aniso_corr"] - q_support) / q_support)
        if q_support > 0 else 0.0,
        "lf_truncation": float(1.0 - q_support / q_full) if q_full > 0
        else float("nan"),
        "closure_raw": split["closure_raw"],
        "closure_corr": split["closure_corr"],
        "kernel_identity": kern_err,
        "stone_none": stone["none"], "stone_piston": stone["piston"],
        "stone_pt": stone_pt,
        # The MODE-MATCHED pairing. The FAST modal mask holds the Zernike modes
        # 1..ZMAX, so it KEEPS the piston and the two tilts. The Stone set that
        # holds the same modes is remove='none' over the band 0..max_order.
        "ratio_quad_stone_none": float(q_full / stone["none"])
        if stone["none"] > 0 else float("nan"),
        "ratio_fast_stone_none": float(split["aniso_corr"] / stone["none"])
        if stone["none"] > 0 else float("nan"),
        # The PRODUCTION pairing. uplink_point_ahead_term removes the piston
        # and the two tilts, because a separate tracking loop points the beam.
        "ratio_fast_stone_pt": float(split["aniso_corr"] / stone_pt)
        if stone_pt > 0 else float("nan"),
        "fit_fast_raw": float(sim.fitting_error),
        "fit_corr": split["fit_corr"],
        "noll": noll,
        "r0_los_m": float(sim.r0_los),
        "seconds": float(time.time() - t0),
    }


def _uncorrected_row(elev_deg, theta_arcsec):
    """Give the ZMAX = 0 (uncorrected) anchor row of the order sweep.

    No corrected mode exists, so no correction decorrelates, and the
    point-ahead residual is EXACTLY zero by definition on both routes. No
    simulation runs. The fitting (total uncorrected) variance diverges in the
    Kolmogorov limit, so the fitting columns hold NaN.
    """
    return {
        "stage": "B-order", "label": "ZMAX 0 (uncorrected)",
        "elev": float(elev_deg), "zmax": 0, "max_order": -1,
        "theta_arcsec": float(theta_arcsec),
        "npxls": 0, "dx": float("nan"), "df": float("nan"),
        "L0_m": float("inf"),
        "fast_raw": 0.0, "fast_corr": 0.0,
        "quad_support": 0.0, "quad_full": 0.0, "quad_rel": 0.0,
        "lf_truncation": float("nan"),
        "closure_raw": 0.0, "closure_corr": 0.0, "kernel_identity": 0.0,
        "stone_none": 0.0, "stone_piston": 0.0, "stone_pt": 0.0,
        "ratio_quad_stone_none": float("nan"),
        "ratio_fast_stone_none": float("nan"),
        "ratio_fast_stone_pt": float("nan"),
        "fit_fast_raw": float("nan"), "fit_corr": float("nan"),
        "noll": float("nan"), "r0_los_m": float("nan"), "seconds": 0.0,
    }


def _verdict(ratios):
    """Give the mechanical verdict string of a stage.

    The rule reads the WORST distance of the FAST-to-Stone ratio from 1.
    """
    worst = max(abs(r - 1.0) for r in ratios if np.isfinite(r))
    if worst <= VERDICT_MATCH:
        return "MATCH", worst
    if worst <= VERDICT_DIFFERENCE:
        return "MEASURED DIFFERENCE", worst
    return "INVESTIGATE", worst


_HEAD = (f"{'case':<22}{'theta[as]':>9}{'ZMAX':>5}{'FASTraw':>9}"
         f"{'FASTcorr':>9}{'quad':>8}{'LFtrunc':>8}{'Stone_n':>9}"
         f"{'Stone_pt':>9}{'Q/S_n':>7}{'F/S_n':>7}{'F/S_pt':>8}")


def _fmt(row):
    """Format one case row for the console table."""
    return (f"{row['label']:<22}{row['theta_arcsec']:9.2f}{row['zmax']:5d}"
            f"{row['fast_raw']:9.4f}{row['fast_corr']:9.4f}"
            f"{row['quad_full']:8.4f}{row['lf_truncation']:8.3f}"
            f"{row['stone_none']:9.4f}"
            f"{row['stone_pt']:9.4f}{row['ratio_quad_stone_none']:7.3f}"
            f"{row['ratio_fast_stone_none']:7.3f}"
            f"{row['ratio_fast_stone_pt']:8.3f}")


def _scenario(n_modes):
    """Build the production pre-compensated uplink scenario."""
    return SpaceScenario(
        ground=Terminal(aperture_m=D_GROUND, wavelength_m=LAM,
                        transmitter=Transmitter(waist_m=0.2),
                        compensation=[TipTilt(), AO(int(n_modes))]),
        space=Terminal(aperture_m=0.05, wavelength_m=LAM),
        direction="uplink",
        channel=Channel(site=Site(cn2_ground=CN2_GROUND), altitude_m=ALT_M),
        precompensation=DownlinkBeacon())


def _default_servo_params(elev_deg, hs, cn2, zmax, theta_arcsec):
    """Build the FAST parameter dict of a DEFAULT-SERVO run.

    These are the servo and wavefront-sensor values that uplink_fast_term ships:
    DSUBAP = 0.02 m, TLOOP = 0.001 s, TEXP = 0.001 s, ALIAS on, and the HV5/7
    wind. The grid is automatic, the same as the Term uses.
    """
    return dict(
        WVL=LAM, D_GROUND=D_GROUND, OBSC_GROUND=0.0, W0=0.2,
        PROP_DIR="up", SMF=False, D_SAT=0.05, OBSC_SAT=0.0, H_SAT=ALT_M,
        ZENITH_ANGLE=90.0 - elev_deg, H_TURB=hs, CN2_TURB=_cn2_layers(cn2, hs),
        WIND_SPD=WIND_RMS * np.ones_like(hs), WIND_DIR=np.zeros_like(hs),
        DTHETA=[float(theta_arcsec), 0.0], AO_MODE="AO", MODAL=True,
        ZMAX=int(zmax), DSUBAP=0.02, TLOOP=0.001, TEXP=0.001, ALIAS=True,
        NOISE=0, SUBHARM=True, L0=np.inf, l0=1e-6, NITER=2, NCHUNKS=1,
        LOGLEVEL="ERROR")


def _stage_c_point(fast, funcs, spectra, n_modes, elev_deg, hs, cn2, niter,
                   grid):
    """Measure ONE stage-C operating point.

    It gives the fidelity-1 FAST Term, the fidelity-0 analytic pair, and the
    attribution ladder between them.
    """
    t0 = time.time()
    scn = _scenario(n_modes)
    geom = CircularOrbit(ALT_M, elevation_deg=elev_deg)
    theta_arcsec = float(geom.point_ahead_rad) * ARCSEC_PER_RAD

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        term_mc = uplink_fast_term(scn, geom, n_samples=int(niter))
        term_pa = uplink_point_ahead_term(scn, geom)
        term_fit = uplink_fitting_term(scn, geom)

    # The matched servo-off sim at this point.
    sim_off = _build(fast, _servo_off_params(elev_deg, hs, cn2, n_modes,
                                             theta_arcsec, grid))
    split_off = split_errors(sim_off, funcs)
    assert abs(split_off["closure_corr"]) < GATE_CLOSURE, split_off

    # The default-servo sim: the residual that the Monte Carlo actually draws.
    sim_on = _build(fast, _default_servo_params(elev_deg, hs, cn2, n_modes,
                                                theta_arcsec))
    servo_sigma2 = float(sim_on.aniso_servo_error + sim_on.alias_error
                         + sim_on.noise_error)

    order = max_radial_order(int(n_modes))
    stone = _stone(D_GROUND, float(geom.point_ahead_rad), hs, cn2, elev_deg,
                   order)
    return {
        "stage": "C", "label": f"AO({n_modes}) @ {elev_deg:.0f} deg",
        "n_modes": int(n_modes), "elev": float(elev_deg),
        "theta_arcsec": theta_arcsec, "niter": int(niter),
        "stone_pt": stone["piston_tilt"], "stone_none": stone["none"],
        "fast_servo_off_sigma2": split_off["aniso_corr"],
        "fast_servo_on_sigma2": servo_sigma2,
        "fast_servo_on_npxls": int(sim_on.Npxls),
        "fast_servo_on_df": float(sim_on.freq.main.df),
        "stone_pt_db": _marechal_db(stone["piston_tilt"]),
        "fast_servo_off_db": _marechal_db(split_off["aniso_corr"]),
        "fast_servo_on_db": _marechal_db(servo_sigma2),
        "term_point_ahead_db": float(term_pa.mean_db),
        "term_fitting_db": float(term_fit.mean_db),
        "fidelity0_db": float(term_pa.mean_db + term_fit.mean_db),
        "fidelity1_db": float(term_mc.mean_db),
        "seconds": float(time.time() - t0),
    }


def _renderer_works():
    """Say whether this matplotlib build can write a PNG.

    The matplotlib 3.11.1 and numpy 2.4.6 pair in this environment faults the
    interpreter inside the Agg renderer. That fault is native, so a try and
    except block cannot catch it. This probe saves one small figure in a
    SEPARATE process and it reads the exit code. The study then keeps its
    numbers even where the renderer is broken.
    """
    code = ("import matplotlib;matplotlib.use('Agg');"
            "import matplotlib.pyplot as plt;import io;"
            "f,a=plt.subplots();a.plot([0,1],[0,1]);"
            "f.savefig(io.BytesIO(),format='png')")
    try:
        done = subprocess.run([sys.executable, "-c", code], timeout=180,
                              stdout=subprocess.DEVNULL,
                              stderr=subprocess.DEVNULL)
        return done.returncode == 0
    except Exception:
        return False


def _figures(rows_theta, rows_order, rows_c, fig_dir):
    """Draw the three figures. Return the file names."""
    os.makedirs(fig_dir, exist_ok=True)
    out = []

    f1 = os.path.join(fig_dir, "residual_vs_pointahead.png")
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    th = [r["theta_arcsec"] for r in rows_theta]
    ax.plot(th, [r["fast_corr"] for r in rows_theta], "o-",
            label="FAST clean split (on the FAST grid)")
    ax.plot(th, [r["quad_full"] for r in rows_theta], "s--",
            label="FAST filter, whole-plane quadrature")
    ax.plot(th, [r["stone_none"] for r in rows_theta], "^-",
            label="Stone, all corrected orders (0 to 9)")
    ax.plot(th, [r["stone_pt"] for r in rows_theta], "v-",
            label="Stone, orders 2 to 9 (production: tilt "
                  "charged to tracking)")
    ax.set_xlabel("point-ahead angle [arcsec]")
    ax.set_ylabel("residual phase variance [rad$^2$]")
    ax.set_title("Point-ahead residual against the angle\n"
                 "1.5 m, 60 deg, 55 corrected modes (radial orders 0 to 9)")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(f1, dpi=130)
    plt.close(fig)
    out.append(f1)

    f2 = os.path.join(fig_dir, "residual_vs_order.png")
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    # A LINEAR axis: the ZMAX = 0 (uncorrected) anchor is exactly zero, and a
    # log axis cannot hold it.
    zs = [r["zmax"] for r in rows_order]
    ax.plot(zs, [r["fast_corr"] for r in rows_order], "o-",
            label="FAST clean split (on the FAST grid)")
    ax.plot(zs, [r["quad_full"] for r in rows_order], "s--",
            label="FAST filter, whole-plane quadrature")
    ax.plot(zs, [r["stone_none"] for r in rows_order], "^-",
            label="Stone, all corrected orders")
    ax.plot(zs, [r["stone_pt"] for r in rows_order], "v-",
            label="Stone, corrected orders 2 up (production: tilt "
                  "charged to tracking)")
    ax.plot(zs, [r["fit_corr"] for r in rows_order], "d:",
            label="FAST fitting (uncorrected band)")
    ax.set_xlabel("corrected Zernike modes ZMAX")
    ax.set_ylabel("residual phase variance [rad$^2$]")
    theta_as = rows_order[-1]["theta_arcsec"]
    ax.set_title(f"Residual against the corrected order, 60 deg, "
                 f"point-ahead {theta_as:.2f} arcsec")
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(f2, dpi=130)
    plt.close(fig)
    out.append(f2)

    f3 = os.path.join(fig_dir, "stage_c_ladder.png")
    fig, ax = plt.subplots(figsize=(7.0, 4.2))
    names = [r["label"] for r in rows_c]
    x = np.arange(len(names))
    w = 0.2
    ax.bar(x - 1.5 * w, [r["stone_pt_db"] for r in rows_c], w,
           label="Stone piston+tilt")
    ax.bar(x - 0.5 * w, [r["fast_servo_off_db"] for r in rows_c], w,
           label="FAST servo off")
    ax.bar(x + 0.5 * w, [r["fast_servo_on_db"] for r in rows_c], w,
           label="FAST default servo")
    ax.bar(x + 1.5 * w, [r["fidelity1_db"] for r in rows_c], w,
           label="fidelity-1 Term (Monte Carlo)")
    ax.plot(x, [r["fidelity0_db"] for r in rows_c], "k*", markersize=12,
            label="fidelity-0 Term pair")
    ax.set_xticks(x)
    ax.set_xticklabels(names, fontsize=8)
    ax.set_ylabel("loss [dB]")
    ax.set_title("Stage C: the attribution ladder")
    ax.grid(alpha=0.3, axis="y")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(f3, dpi=130)
    plt.close(fig)
    out.append(f3)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--full", action="store_true",
                    help="the long run: the fine grid, every sweep point, and "
                         "3000 Monte Carlo draws")
    ap.add_argument("--niter", type=int, default=None,
                    help="stage-C Monte Carlo draws (default 500 quick, 3000 "
                         "full)")
    args = ap.parse_args()

    t_start = time.time()
    fast, funcs, spectra = _load()
    grid = GRID_FULL if args.full else GRID_QUICK
    niter = args.niter or (3000 if args.full else 500)

    lines = []

    def say(text=""):
        print(text)
        lines.append(text)

    say("FAST against Stone: the point-ahead anisoplanatism residual "
        "(backlog 1-5)")
    say(f"mode          : {'FULL' if args.full else 'QUICK'}")
    say(f"grid          : NPXLS = {grid[0]}, DX = {grid[1]} m  "
        f"(df = {2 * np.pi / (grid[0] * grid[1]):.3f} rad/m, "
        f"kappa_max = {np.pi / grid[1]:.0f} rad/m)")
    say(f"link          : D = {D_GROUND} m, lambda = {LAM * 1e9:.0f} nm, "
        f"orbit {ALT_M / 1e3:.0f} km, HV5/7 (Cn2_0 = {CN2_GROUND:.1e})")
    say(f"gates         : quadrature {GATE_QUADRATURE:.0%}, closure "
        f"{GATE_CLOSURE:.1%}, convergence {GATE_CONVERGENCE:.0%}, "
        f"zero point-ahead {GATE_ZERO_RAD2:g} rad^2")
    say(f"outer scale   : Kolmogorov (L0 = infinite) everywhere, except the "
        f"GATED legs, which use L0 = {L0_GATE:g} m")
    say()

    gates = []
    rows = []

    # ---------------------------------------------------------------- A0 -----
    say("STAGE A0. One turbulence layer at the zenith.")
    hs_a0 = np.linspace(5e3, 15e3, 81)
    cn2_a0 = np.exp(-0.5 * ((hs_a0 - 10e3) / 250.0) ** 2)
    cn2_a0 *= 1e-13 / np.trapz(cn2_a0, hs_a0)
    theta_a0 = 9.0
    say(f"  a Gaussian Cn2 bump at 10 km, sigma = 250 m, "
        f"INT Cn2 dh = 1e-13 m^(1/3), theta = {theta_a0:.1f} arcsec, "
        f"ZMAX = {ZMAX_NOM}")
    say(f"  A0 always runs on the {GRID_A0[0]} / {GRID_A0[1]:g} m grid. Its "
        "profile holds 81 layers, and FAST holds one spectrum for each layer, "
        "so a finer grid needs several gigabytes.")

    # Leg 1: the GATE. A finite outer scale puts the whole integral inside the
    # grid, so the FAST value and the independent quadrature must agree.
    a0g = _row(fast, funcs, spectra, "A0-gate",
               f"single layer, L0 = {L0_GATE:g} m", 90.0, ZMAX_NOM, theta_a0,
               hs_a0, cn2_a0, GRID_A0, L0=L0_GATE)
    rows.append(a0g)
    say()
    say(f"  THE 1 % GATE, with the outer scale L0 = {L0_GATE:g} m, over the "
        f"square support {np.pi / GRID_A0[1]:.0f} rad/m:")
    say(f"    FAST grid, clean split : {a0g['fast_corr']:.5f} rad^2")
    say(f"    independent quadrature : {a0g['quad_support']:.5f} rad^2")
    say(f"    relative difference    : {a0g['quad_rel']:+.4f}")
    say(f"    kernel identity |G_ao - rebuilt| : "
        f"{a0g['kernel_identity']:.2e}")
    say(f"    closure, raw / clean   : {a0g['closure_raw']:+.2e} / "
        f"{a0g['closure_corr']:+.2e}")
    assert abs(a0g["quad_rel"]) < GATE_QUADRATURE, a0g["quad_rel"]
    gates.append({"gate": "A0 quadrature", "value": a0g["quad_rel"],
                  "limit": GATE_QUADRATURE, "pass": True})

    # Leg 2: the Kolmogorov measurement against Stone.
    a0 = _row(fast, funcs, spectra, "A0", "single layer, Kolmogorov", 90.0,
              ZMAX_NOM, theta_a0, hs_a0, cn2_a0, GRID_A0)
    rows.append(a0)
    say()
    say("  THE KOLMOGOROV LEG (L0 = infinite), which is what Stone "
        "integrates:")
    say(f"    FAST raw aniso_servo_error : {a0['fast_raw']:.5f} rad^2")
    say(f"    FAST clean split aniso_corr: {a0['fast_corr']:.5f} rad^2")
    say(f"    FAST fitting raw / clean   : {a0['fit_fast_raw']:.5f} / "
        f"{a0['fit_corr']:.5f} rad^2")
    say(f"    quadrature, whole plane    : {a0['quad_full']:.5f} rad^2")
    say(f"    the FAST grid MISSES       : {a0['lf_truncation']:.1%} of it "
        f"(it has no support below df = {a0['df']:.3f} rad/m, and the "
        "integral converges only as kappa^(1/3))")
    say(f"    Stone, no mode removed     : {a0['stone_none']:.5f} rad^2   "
        f"(band 0..{a0['max_order']})")
    say(f"    Stone, piston removed      : {a0['stone_piston']:.5f} rad^2")
    say(f"    Stone, piston and tilt     : {a0['stone_pt']:.5f} rad^2")
    say(f"    quadrature / Stone none    : "
        f"{a0['quad_full'] / a0['stone_none']:.4f}   "
        "(the clean soft filter against the Stone projection)")
    say()

    # ---------------------------------------------------------------- A1 -----
    say("STAGE A1. The HV5/7 profile at the production point.")
    hs = DEFAULT_HS
    cn2 = get_c2n(hs, WIND_RMS, CN2_GROUND)
    geom60 = CircularOrbit(ALT_M, elevation_deg=60.0)
    theta_nom = float(geom60.point_ahead_rad) * ARCSEC_PER_RAD
    a1 = _row(fast, funcs, spectra, "A1", "HV5/7 60 deg", 60.0, ZMAX_NOM,
              theta_nom, hs, cn2, grid)
    rows.append(a1)
    noll60 = _noll(hs, cn2, 60.0, D_GROUND, 60)
    say(f"  theta = {theta_nom:.3f} arcsec ({geom60.point_ahead_rad * 1e6:.2f} "
        f"urad), r0 line of sight = {a1['r0_los_m'] * 100:.1f} cm")
    say(f"  FAST raw aniso_servo_error   : {a1['fast_raw']:.5f} rad^2")
    say(f"  FAST clean split aniso_corr  : {a1['fast_corr']:.5f} rad^2")
    say(f"  closure, raw / clean         : {a1['closure_raw']:+.2e} / "
        f"{a1['closure_corr']:+.2e}")
    say(f"  quadrature, whole plane      : {a1['quad_full']:.5f} rad^2   "
        f"(the FAST grid misses {a1['lf_truncation']:.1%})")
    say(f"  Stone none / piston / piston+tilt : {a1['stone_none']:.5f} / "
        f"{a1['stone_piston']:.5f} / {a1['stone_pt']:.5f} rad^2")
    say(f"  the isolated piston term     : "
        f"{a1['stone_none'] - a1['stone_piston']:.5f} rad^2")
    say(f"  the isolated tilt term       : "
        f"{a1['stone_piston'] - a1['stone_pt']:.5f} rad^2")
    say(f"  FAST fitting raw / clean     : {a1['fit_fast_raw']:.5f} / "
        f"{a1['fit_corr']:.5f} rad^2")
    say(f"  Noll fitting, 55 / 60 modes  : {a1['noll']:.5f} / {noll60:.5f} "
        "rad^2")
    say("  THE PAIRINGS:")
    say(f"    mode matched, no grid truncation (quadrature / Stone none) : "
        f"{a1['ratio_quad_stone_none']:.3f}")
    say(f"    mode matched, on the FAST grid   (FAST / Stone none)       : "
        f"{a1['ratio_fast_stone_none']:.3f}   "
        f"(the grid misses {a1['lf_truncation']:.1%})")
    say(f"    production pairing               (FAST / Stone piston+tilt): "
        f"{a1['ratio_fast_stone_pt']:.3f}   "
        "(the FAST modal mask KEEPS the piston and the tilt, and "
        "uplink_point_ahead_term removes them, so this pair is NOT mode "
        "matched)")
    say()

    # ----------------------------------------------------------------- B -----
    say("STAGE B. The sweeps.")
    say(_HEAD)

    scale = ((0.0, 0.5, 1.0, 2.0) if not args.full
             else (0.0, 0.25, 0.5, 1.0, 1.5, 2.0))
    rows_theta = []
    for s in scale:
        r = _row(fast, funcs, spectra, "B-theta", f"theta x {s:g}", 60.0,
                 ZMAX_NOM, theta_nom * s, hs, cn2, grid)
        rows_theta.append(r)
        rows.append(r)
        say(_fmt(r))
    zero = [r for r in rows_theta if r["theta_arcsec"] == 0.0][0]
    say(f"  theta = 0: FAST corrected {zero['fast_corr']:.3e} rad^2, "
        f"Stone piston+tilt {zero['stone_pt']:.3e} rad^2, "
        f"FAST RAW {zero['fast_raw']:.5f} rad^2 (the soft-mask leakage)")
    assert zero["fast_corr"] < GATE_ZERO_RAD2, zero["fast_corr"]
    assert zero["stone_pt"] < GATE_ZERO_RAD2, zero["stone_pt"]
    gates.append({"gate": "zero point-ahead",
                  "value": max(zero["fast_corr"], zero["stone_pt"]),
                  "limit": GATE_ZERO_RAD2, "pass": True})
    say()

    zmaxes = ((1, 3, 6, 21, 55) if not args.full
              else (1, 3, 6, 10, 21, 36, 55, 66))
    rows_order = []
    say(_HEAD)
    # The UNCORRECTED anchor (ZMAX = 0). With no corrected mode there is no
    # correction to decorrelate, so the point-ahead residual is EXACTLY zero on
    # both sides, with no simulation. The fitting (total uncorrected) variance
    # diverges in the Kolmogorov limit, so that column is not defined here.
    r0row = _uncorrected_row(60.0, theta_nom)
    rows_order.append(r0row)
    rows.append(r0row)
    say(_fmt(r0row))
    for z in zmaxes:
        r = _row(fast, funcs, spectra, "B-order", f"ZMAX {z}", 60.0, z,
                 theta_nom, hs, cn2, grid)
        rows_order.append(r)
        rows.append(r)
        say(_fmt(r))
    r60 = _row(fast, funcs, spectra, "B-order", "ZMAX 60 (production)", 60.0,
               60, theta_nom, hs, cn2, grid)
    rows_order.append(r60)
    rows.append(r60)
    say(_fmt(r60))
    say(f"  the production point corrects 60 Noll modes. olb reads it as "
        f"radial order {r60['max_order']}, because order 9 is complete at 55 "
        "modes.")
    say()

    elevs = (30.0, 60.0) if not args.full else (30.0, 60.0, 90.0)
    rows_elev = []
    say(_HEAD)
    for e in elevs:
        th = float(CircularOrbit(ALT_M, elevation_deg=e).point_ahead_rad) \
            * ARCSEC_PER_RAD
        r = _row(fast, funcs, spectra, "B-elev", f"elev {e:.0f} deg", e,
                 ZMAX_NOM, th, hs, cn2, grid)
        rows_elev.append(r)
        rows.append(r)
        say(_fmt(r))
    say()

    # -------------------------------------------------------- convergence ----
    say("CONVERGENCE. The A1 case on three grids.")
    say(f"  The GATED rows carry the outer scale L0 = {L0_GATE:g} m, because "
        "the Kolmogorov integral has no low-frequency floor and a grid cannot "
        "converge on it. The Kolmogorov rows sit below, and they are NOT "
        "gated.")
    # One grid holds df alone (a wider domain) and one holds kappa_max alone (a
    # finer pixel). No grid goes past 1024 pixels, because a bigger grid needs
    # several gigabytes for the per-layer spectra.
    conv_grids = ([grid, (1024, 0.02), (1024, 0.01)] if grid == GRID_QUICK
                  else [grid, (1024, 0.02), (512, 0.01)])
    conv, conv_k = [], []
    for g in conv_grids:
        conv.append(_row(fast, funcs, spectra, "conv",
                         f"N={g[0]} dx={g[1]:g}", 60.0, ZMAX_NOM, theta_nom,
                         hs, cn2, g, L0=L0_GATE))
        conv_k.append(_row(fast, funcs, spectra, "conv-kolmogorov",
                           f"N={g[0]} dx={g[1]:g}", 60.0, ZMAX_NOM, theta_nom,
                           hs, cn2, g))
    rows.extend(conv)
    rows.extend(conv_k)
    say(f"{'grid':<18}{'df':>8}{'kappa_max':>10}{'aniso_corr':>12}"
        f"{'quadrature':>12}{'fit_corr':>11}")
    for r in conv:
        say(f"{r['label']:<18}{r['df']:8.3f}{np.pi / r['dx']:10.0f}"
            f"{r['fast_corr']:12.5f}{r['quad_support']:12.5f}"
            f"{r['fit_corr']:11.5f}")
    base = conv[0]
    same_kmax = [r for r in conv[1:] if r["dx"] == base["dx"]]
    other_kmax = [r for r in conv[1:] if r["dx"] != base["dx"]]
    # The gate reads the REFINED grids only (df at or below the base df). A
    # coarser df under-resolves the outer-scale knee at k0 = 2 pi / L0_GATE,
    # so its movement measures the degraded grid, not the base grid.
    refined = [r for r in conv[1:] if r["df"] <= base["df"] + 1e-12]
    coarser = [r for r in conv[1:] if r["df"] > base["df"] + 1e-12]
    moves_a = [abs(r["fast_corr"] / base["fast_corr"] - 1.0) for r in refined]
    coarse_a = [abs(r["fast_corr"] / base["fast_corr"] - 1.0) for r in coarser]
    moves_f = [abs(r["fit_corr"] / base["fit_corr"] - 1.0) for r in same_kmax]
    tail_f = [abs(r["fit_corr"] / base["fit_corr"] - 1.0) for r in other_kmax]
    say("  movement of aniso_corr, refined df         : "
        + ", ".join(f"{m:.2%}" for m in moves_a))
    if coarse_a:
        say("  movement of aniso_corr, HALF-domain grid   : "
            + ", ".join(f"{m:.2%}" for m in coarse_a)
            + "   (NOT gated: its df sits at the outer-scale knee, so it "
              "under-resolves the low band; the base grid holds it)")
    say("  movement of fit_corr, same kappa_max       : "
        + ", ".join(f"{m:.2%}" for m in moves_f))
    say("  movement of fit_corr, kappa_max is DOUBLED : "
        + ", ".join(f"{m:.2%}" for m in tail_f)
        + "   (NOT gated: a finer pixel opens a higher frequency band, and "
          "the Kolmogorov fitting tail there is real variance)")
    assert max(moves_a) < GATE_CONVERGENCE, moves_a
    assert max(moves_f) < GATE_CONVERGENCE, moves_f
    gates.append({"gate": "convergence, aniso", "value": max(moves_a),
                  "limit": GATE_CONVERGENCE, "pass": True})
    gates.append({"gate": "convergence, fitting", "value": max(moves_f),
                  "limit": GATE_CONVERGENCE, "pass": True})
    say("  the same three grids in the KOLMOGOROV limit (not gated):")
    say(f"{'grid':<18}{'df':>8}{'kappa_max':>10}{'aniso_corr':>12}"
        f"{'whole plane':>13}{'LFtrunc':>9}")
    for r in conv_k:
        say(f"{r['label']:<18}{r['df']:8.3f}{np.pi / r['dx']:10.0f}"
            f"{r['fast_corr']:12.5f}{r['quad_full']:13.5f}"
            f"{r['lf_truncation']:9.3f}")
    say("  A finer df ADDS real variance, because it opens a lower frequency "
        "band. That movement is the low-frequency truncation, not a numerical "
        "error.")
    say()

    # ----------------------------------------------------------------- C -----
    say(f"STAGE C. The full Term comparison in dB (NITER = {niter}).")
    points = [(60, 60.0), (60, 30.0), (21, 60.0)]
    if args.full:
        points += [(60, 90.0), (10, 60.0)]
    rows_c = []
    for n_modes, elev in points:
        rows_c.append(_stage_c_point(fast, funcs, spectra, n_modes, elev, hs,
                                     cn2, niter, grid))
    say(f"{'point':<20}{'theta[as]':>9}{'Stone_pt':>10}{'FASToff':>10}"
        f"{'FASTservo':>11}{'->dB Stone':>11}{'->dB off':>10}"
        f"{'->dB servo':>11}{'fid0 dB':>9}{'fid1 dB':>9}")
    for r in rows_c:
        say(f"{r['label']:<20}{r['theta_arcsec']:9.2f}{r['stone_pt']:10.4f}"
            f"{r['fast_servo_off_sigma2']:10.4f}"
            f"{r['fast_servo_on_sigma2']:11.4f}{r['stone_pt_db']:11.3f}"
            f"{r['fast_servo_off_db']:10.3f}{r['fast_servo_on_db']:11.3f}"
            f"{r['fidelity0_db']:9.3f}{r['fidelity1_db']:9.3f}")
    say("  fid0 dB = uplink_point_ahead_term + uplink_fitting_term (mean "
        "only). fid1 dB = uplink_fast_term Monte Carlo mean.")
    for r in rows_c:
        say(f"  {r['label']}: point-ahead {r['term_point_ahead_db']:.3f} dB + "
            f"fitting {r['term_fitting_db']:.3f} dB = "
            f"{r['fidelity0_db']:.3f} dB   against the Monte Carlo "
            f"{r['fidelity1_db']:.3f} dB   "
            f"(the servo sim grid: NPXLS = {r['fast_servo_on_npxls']}, "
            f"df = {r['fast_servo_on_df']:.2f} rad/m)")
    say()

    # ----------------------------------------------------------- verdicts ----
    stage_rows = {"A0": [a0], "A1": [a1], "B-theta": [r for r in rows_theta
                                                      if r["theta_arcsec"] > 0],
                  "B-order": rows_order, "B-elev": rows_elev}
    verdicts = {}
    say("VERDICTS. The rule reads the MODE-MATCHED ratio Q/S_n: the whole-plane "
        "quadrature of the FAST filter against the Stone band 0..max_order with "
        "no mode removed. Both sets then hold the SAME modes, and the "
        "quadrature carries no grid truncation. The rule takes the worst case: "
        f"within {VERDICT_MATCH:.0%} MATCH, to {VERDICT_DIFFERENCE:.0%} "
        "MEASURED DIFFERENCE, else INVESTIGATE.")
    for stage, rs in stage_rows.items():
        v, worst = _verdict([r["ratio_quad_stone_none"] for r in rs])
        verdicts[stage] = v
        worst_grid = max(abs(r["ratio_fast_stone_none"] - 1.0) for r in rs
                         if np.isfinite(r["ratio_fast_stone_none"]))
        say(f"  {stage:<8} {v:<20} worst |Q/S_n - 1| = {worst:.3f}   "
            f"(on the FAST GRID, worst |F/S_n - 1| = {worst_grid:.3f})")
    c_ratios = [r["fidelity1_db"] / r["fidelity0_db"] for r in rows_c]
    v_c, worst_c = _verdict(c_ratios)
    verdicts["C"] = v_c
    say(f"  {'C':<8} {v_c:<20} worst |fidelity1/fidelity0 - 1| = "
        f"{worst_c:.3f}")
    say()
    if not args.full:
        say("QUICK MODE: run --full before you record a verdict.")

    elapsed = time.time() - t_start
    say(f"runtime: {elapsed:.1f} s")
    draw = _renderer_works()
    if not draw:
        say("FIGURES SKIPPED: this matplotlib build faults inside the Agg "
            "renderer, so the script keeps the log and the results file only.")

    stamp = {
        "mode": "full" if args.full else "quick",
        "grid": {"npxls": grid[0], "dx": grid[1]},
        "niter": niter,
        "parameters": {
            "wavelength_m": LAM, "d_ground_m": D_GROUND, "w0_m": W0_M,
            "d_sat_m": D_SAT, "altitude_m": ALT_M,
            "cn2_ground": CN2_GROUND, "wind_rms_m_s": WIND_RMS,
            "dsubap_servo_off_m": DSUBAP_OFF, "zmax_nominal": ZMAX_NOM,
            "l0_gate_m": L0_GATE,
            "kappa_min_full_rad_m": KAPPA_MIN_FULL,
            "kappa_max_quad_rad_m": KAPPA_MAX_QUAD,
            "n_kappa": N_KAPPA, "n_phi": N_PHI,
            "theta_nominal_arcsec": theta_nom,
        },
        "gates": gates,
        "verdicts": verdicts,
        "figures": bool(draw),
        "rows": rows,
        "rows_stage_c": rows_c,
        "seconds": elapsed,
    }
    json_path = os.path.join(HERE, "fast_stone_pointahead_results.json")
    with open(json_path, "w") as fh:
        json.dump(stamp, fh, indent=2)
    log_path = os.path.join(HERE, "fast_stone_pointahead.log")
    with open(log_path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    figs = _figures_subprocess() if draw else []
    print(f"\nwrote {json_path}")
    for f in figs:
        print(f"wrote {f}")
    print(f"wrote {log_path}")


def _figures_subprocess():
    """Draw the figures in a SEPARATE process, from the saved JSON.

    The renderer probe of _renderer_works passes in a clean process, but the
    same render FAULTS natively inside the study process after FAST and its
    libraries load. So the study writes its JSON first, and then it renders
    from that JSON in a fresh process (--figures-only). A render fault then
    costs the figures only, never the numbers.
    """
    done = subprocess.run([sys.executable, "-m",
                           "validation.fast_stone_pointahead."
                           "fast_stone_pointahead", "--figures-only"],
                          capture_output=True, text=True, timeout=600,
                          cwd=os.path.join(HERE, "..", ".."))
    if done.returncode != 0:
        print("FIGURES SKIPPED: the render subprocess failed "
              f"(exit {done.returncode}).")
        return []
    return [ln for ln in done.stdout.splitlines() if ln.endswith(".png")]


def _figures_only():
    """Render the figures from the saved JSON. Run in a clean process."""
    with open(os.path.join(HERE, "fast_stone_pointahead_results.json")) as fh:
        d = json.load(fh)
    rows_theta = [r for r in d["rows"] if r["stage"] == "B-theta"]
    rows_order = [r for r in d["rows"] if r["stage"] == "B-order"]
    figs = _figures(rows_theta, rows_order, d["rows_stage_c"],
                    os.path.join(HERE, "figures"))
    print("\n".join(figs))


if __name__ == '__main__':
    if "--figures-only" in sys.argv:
        _figures_only()
    else:
        main()
