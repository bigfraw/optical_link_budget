"""The TILT anisoplanatism of the point-ahead record, in rad^2 (backlog 2-P4).

THE QUESTION (owner, 2026-09-11). The p1 study of this folder reads the
fidelity-1 FAST Term about 0.3 dB ABOVE the field on the TIP-TILT stack at
30 deg, with the same sign at every angle, after the known FAST leak is taken
away. At 20 deg that tip-tilt offset goes away and an excess appears on the AO
rows. A dB number holds two things together: the phase variance and the
extended Marechal map of that variance to a loss. This script takes the map
AWAY. It measures the TILT DECORRELATION alone, in rad^2, from the STORED
planes of the `base` campaigns, and it puts the field, Stone and FAST next to
each other mode by mode.

THE THREE ROUTES to the same quantity.

  1. THE SCREEN ROUTE (geometric optics, the Stone quantity). Each trial
     stores the summed screen phase of the BEACON window, and the same sum
     through each SHIFTED window. Their difference

         d_i = screen_phase_pa[i] - screen_phase

     IS the anisoplanatic wavefront error of the angle i. The script fits 21
     Noll modes to d_i over the 0.7 m aperture mask and it gives the phase
     variance of each Noll band. THIS is what Stone Eq. (29) computes.

  2. THE FIELD ROUTE (propagated). The tilt of the propagated point-ahead
     field minus the tilt of the propagated beacon field, from the wrapped
     slopes and from the far-field (G-tilt) centroid. Diffraction and
     scintillation are in this route and not in route 1.

  3. THE REFERENCES. The continuous Stone integral at L0 = 25 m and at
     L0 = infinite, the DISCRETE delta-layer Stone sum on the campaign's OWN
     screen plan with the pixel-rounded shifts, and the FAST tilt band.

THE BAND UNITS. Every band is a PHASE VARIANCE over the aperture mask, in
rad^2, and the bands ADD UP to the mask variance of the piston-removed
difference. The Noll basis is only NEAR-orthonormal on a pixel mask, so the
script reports each band as the SHARE

    S_B = < p_B - mean(p_B), p_fit - mean(p_fit) >_mask,

which adds to the fitted variance EXACTLY, next to the plain diagonal value
mean((p_B - mean(p_B))^2). The least-squares residual is orthogonal to the
basis, so

    mask-var(d - piston) = SUM_B S_B + mean(residual^2)

holds to the rounding level. The script asserts that closure at 1e-6.

THE FAST TILT BAND. FAST corrects modally over a SOFT Zernike Fourier mask, so
its tilt band is the difference of two servo-off runs: the clean anisoplanatic
split of ZMAX = 3 (piston and the two tilts) minus the split of ZMAX = 1 (the
piston alone). The helpers come from
`validation/fast_stone_pointahead/fast_stone_pointahead.py`
(`_servo_off_params`, `_build`, `split_errors`), and the FAST grid follows the
rule of that script: an EXPLICIT (NPXLS, DX) pair, not the automatic grid. The
side N dx must hold the outer scale, because the von Karman knee at
k0 = 2 pi / L0 must sit inside the grid.

Sources:
- J. Stone, P. H. Hu, S. P. Mills and S. Ma, J. Opt. Soc. Am. A 11(1), 347-357
  (1994), DOI 10.1364/JOSAA.11.000347. Eqs. (29) and (36).
- R. J. Noll, J. Opt. Soc. Am. 66(3), 207-211 (1976),
  DOI 10.1364/JOSA.66.000207, Tables I and IV. The mode order and the bands.
- Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196, Ch. 3, Eq. (20),
  printed p. 68. The von Karman outer-scale kernel of the Stone integral.
- O. J. D. Farley and others, Opt. Express 30(13), 23050 (2022),
  DOI 10.1364/OE.458659. The FAST method.
- T. S. Ross, Appl. Opt. 48(10), 1812 (2009), DOI 10.1364/AO.48.001812. The
  extended Marechal map (10 / ln 10) sigma^2, printed here as a CONTEXT
  column only.
- J. H. Shapiro, J. Opt. Soc. Am. 61(4), 492-495 (1971),
  DOI 10.1364/JOSA.61.000492. The reciprocity overlap of the stored eta.

THE SCRIPT PROPAGATES NOTHING. It reads the stored planes only.

Run it from the repository root:

    python -m validation.waveoptics_pointahead.tilt_anisoplanatism
    python -m validation.waveoptics_pointahead.tilt_anisoplanatism --no-fast
"""

import argparse
import json
import os
import time
import warnings

import numpy as np

from olb.turbulence.anisoplanatism import _REMOVE_NLO
from olb.waveoptics.compensation.modal import ApertureModes
from olb.waveoptics.compensation.slopes import wrapped_gradient
from olb.waveoptics.compensation.zernike import circle
# The case, the plans and the Stone kernels are IMPORTED. This script owns no
# physics of its own.
from validation.anisoplanatism_screens.anisoplanatism_screens import (
    InnerI, continuous_stone, discrete_stone, u0_of)
from validation.anisoplanatism_screens.common import (ARCSEC, L0_M, LAM,
                                                      hero_uplink)
from validation.waveoptics_ao.waveoptics_ao import _pyplot
from validation.waveoptics_pointahead.waveoptics_pointahead import (
    ENV_ROOT, HERE, angle_labels, campaign_of)

# The ground aperture of the hero case, in m. The mask of every fit is that
# aperture on the stored crop.
APERTURE_M = 0.70
OBSCURATION = 0.0

# The Noll modes that the fit holds. 21 modes are the complete radial orders
# 0 to 5 (Noll, DOI 10.1364/JOSA.66.000207, Table I).
N_MODES = 21

# The Noll bands of the table. Each entry is (label, the 0-based index slice).
# Noll j = 1 is the piston, j = 2 and 3 are the two tilts, and each further
# radial order n holds n + 1 modes.
BANDS = (("piston", (0, 1)), ("tilt", (1, 3)), ("order 2", (3, 6)),
         ("order 3", (6, 10)), ("orders 4-5", (10, 21)))

# The closure gate of the band sum, relative to the mask variance.
CLOSURE_GATE = 1e-6

# The bootstrap of the tilt band.
N_BOOT = 400
BOOT_SEED = 20260911

# The FAST grid. The rule is the rule of
# `validation/fast_stone_pointahead/fast_stone_pointahead.py`: an EXPLICIT
# (NPXLS, DX) pair. The side N dx = 51.2 m holds the 25 m outer scale, so the
# von Karman knee k0 = 2 pi / L0 = 0.25 rad/m sits above the lowest grid
# frequency df = 2 pi / (N dx) = 0.123 rad/m. The corner kappa_max = pi / dx =
# 62.8 rad/m is far above the tilt filter cutoff of a 0.7 m aperture.
FAST_GRID = (1024, 0.05)

# The FAST mode counts of the tilt isolation. ZMAX = 3 holds the piston and the
# two tilts; ZMAX = 1 holds the piston alone.
FAST_ZMAX_TILT = 3
FAST_ZMAX_PISTON = 1

_LN10 = float(np.log(10.0))


# ---------------------------------------------------------------------------
# The per-trial reader
# ---------------------------------------------------------------------------

class BandTrial:
    """The per-trial callable of the band measurement.

    It builds the aperture mask and the modal basis ONE time for each block,
    in the record context, exactly as `_PostCorrector` does.

    THE MASK IS THE APERTURE ON THE CROP. The stored patch carries a 1.5x
    margin, so the aperture sits well inside the crop, and the pixels outside
    the patch disc are zero. The slope route masks its pairs for that reason.

    Attributes:
        aperture_m:  the aperture diameter, in m.
        obscuration: the central obscuration ratio.
        n_modes:     the fitted Noll mode count.
    """

    def __init__(self, aperture_m=APERTURE_M, obscuration=OBSCURATION,
                 n_modes=N_MODES):
        self.aperture_m = float(aperture_m)
        self.obscuration = float(obscuration)
        self.n_modes = int(n_modes)

    def _modes(self, rec):
        """Give the cached ApertureModes of the block."""
        got = rec.context.get("modes")
        if got is None:
            crop = rec.patch.crop()
            mask = circle(crop.side, self.aperture_m / rec.patch.pixel_m,
                          self.obscuration)
            got = ApertureModes(self.n_modes, crop.side, mask)
            rec.context["modes"] = got
        return got

    def _scatter(self, rec, values):
        """Put a stored patch row on the crop, then keep the in-mask pixels."""
        crop = rec.patch.crop()
        flat = np.zeros(crop.side * crop.side, dtype=np.float64)
        flat[crop.indices] = np.asarray(values, dtype=np.float64)
        return flat

    def __call__(self, rec):
        """Measure one trial.

        Returns:
            A float array of the shape (n_angles, 10). The columns are the
            piston, the four band SHARES, the fit residual, the mask variance
            of the piston-removed difference, the diagonal tilt band, the
            slope tilt and the G-tilt tilt, all in rad^2.
        """
        modes = self._modes(rec)
        idx = modes.indices
        basis = modes.basis                      # (n_pix, n_modes)
        n_angles = int(np.asarray(rec.screen_phase_pa).shape[0])

        base_phase = self._scatter(rec, rec.screen_phase)[idx]
        # The field routes sense the BEACON field one time for the trial.
        sx, sy = wrapped_gradient(rec.array, mask=modes.mask)
        c_slope_0 = modes.estimate_from_slopes(sx, sy)
        c_gtilt_0 = modes.estimate_gtilt(rec.array)

        out = np.zeros((n_angles, 10), dtype=np.float64)
        for i in range(n_angles):
            d = self._scatter(rec, rec.screen_phase_pa[i])[idx] - base_phase
            coeffs = modes.estimate(d)
            fit = basis @ coeffs
            # The least-squares fit holds the piston mode, so the residual has
            # a zero mask mean and the split below is exact.
            resid = d - fit
            fit0 = fit - fit.mean()
            out[i, 0] = float(d.mean()) ** 2
            for b, (_label, (lo, hi)) in enumerate(BANDS[1:], start=1):
                p = basis[:, lo:hi] @ coeffs[lo:hi]
                p0 = p - p.mean()
                out[i, b] = float(np.mean(p0 * fit0))
                if b == 1:
                    out[i, 7] = float(np.mean(p0 * p0))
            out[i, 5] = float(np.mean(resid * resid))
            d0 = d - d.mean()
            out[i, 6] = float(np.mean(d0 * d0))

            # The field routes. The tilt of the point-ahead field minus the
            # tilt of the beacon field, in the same mask units.
            fpa = rec.arrays_pa[i]
            sxp, syp = wrapped_gradient(fpa, mask=modes.mask)
            out[i, 8] = _tilt_variance(basis,
                                       modes.estimate_from_slopes(sxp, syp)
                                       - c_slope_0)
            out[i, 9] = _tilt_variance(basis,
                                       modes.estimate_gtilt(fpa) - c_gtilt_0)
        return out


def _tilt_variance(basis, coeffs):
    """Give the mask phase variance of the tilt part of a coefficient set.

    Args:
        basis:  the (n_pix, n_modes) Noll basis over the mask.
        coeffs: the Noll coefficients, in rad.

    Returns:
        The mask variance of the tilt map, in rad^2.
    """
    lo, hi = BANDS[1][1]
    p = basis[:, lo:hi] @ np.asarray(coeffs, dtype=np.float64)[lo:hi]
    p0 = p - p.mean()
    return float(np.mean(p0 * p0))


def _boot_half(values, fn, n_boot=N_BOOT, seed=BOOT_SEED):
    """Give the bootstrap 2-sigma half-width of a statistic over the trials.

    The function resamples the TRIAL rows with replacement and it takes the
    standard deviation of the resampled statistic. Two of those is the
    2-sigma bar of the table.
    """
    x = np.asarray(values, dtype=float)
    rng = np.random.default_rng(seed)
    draws = np.array([fn(x[rng.integers(0, x.shape[0], x.shape[0])])
                      for _ in range(int(n_boot))])
    return float(2.0 * np.std(draws))


# ---------------------------------------------------------------------------
# The references
# ---------------------------------------------------------------------------

_INNER_CACHE = {}


def _inner(max_order, L0):
    """Give the cached InnerI table of one (max_order, L0) pair.

    Each table costs 161 `scipy.integrate.quad` triples, and the sweep asks
    for the same table at every angle. The table does NOT depend on the angle
    or on the elevation, so one build serves the whole run. See
    `validation/anisoplanatism_screens/anisoplanatism_screens.InnerI`.
    """
    key = (max_order, round(float(u0_of(APERTURE_M, L0)), 12))
    if key not in _INNER_CACHE:
        _INNER_CACHE[key] = InnerI(_REMOVE_NLO["piston"], max_order, key[1])
    return _INNER_CACHE[key]


def stone_bands(scn, elevation_deg, theta, L0):
    """Give the continuous Stone variance of each Noll band, in rad^2.

    The bands come from the DIFFERENCE of the cumulative Stone integrals of
    `max_order`. `remove='piston'` keeps the tilt, so max_order = 1 IS the tilt
    band, and each higher order adds its own shell. The last entry is the part
    that sits above the radial order 5, which the 21-mode fit leaves in the
    residual.

    Args:
        scn:           the SpaceScenario.
        elevation_deg: the elevation, in deg.
        theta:         the point-ahead angle, in rad.
        L0:            the outer scale, in m.

    Returns:
        A dict band label -> the variance, in rad^2, plus "total".
    """
    def cum(max_order):
        return continuous_stone(scn, elevation_deg, APERTURE_M, theta, LAM,
                                remove="piston", max_order=max_order,
                                inner=_inner(max_order, L0))

    c1, c2, c3, c5, call = (cum(1), cum(2), cum(3), cum(5), cum(None))
    return {"piston": 0.0, "tilt": c1, "order 2": c2 - c1,
            "order 3": c3 - c2, "orders 4-5": c5 - c3,
            "above order 5": call - c5, "total": call}


def fast_tilt_band(elevation_deg, theta, grid, say):
    """Give the FAST anisoplanatic TILT variance of one angle, in rad^2.

    THE ISOLATION. FAST corrects over a SOFT Zernike Fourier mask, so a single
    run gives the whole corrected band. The tilt alone is the clean
    anisoplanatic split of ZMAX = 3 (the piston and the two tilts) minus the
    split of ZMAX = 1 (the piston). Both runs are SERVO OFF, so the PAOLA
    filter holds the anisoplanatic kernel 2 - 2 cos(delta_r . kappa) alone.

    THE HERO GEOMETRY. `_servo_off_params` carries the 1.5 m case of its own
    module, so this function overwrites the four geometry keys with the hero
    values. Every other key (the servo-off switches, the mode count, the
    grid) stays exactly as that module writes it.

    Args:
        elevation_deg: the elevation, in deg.
        theta:         the point-ahead angle, in rad.
        grid:          the (NPXLS, DX) pair.
        say:           the log function.

    Returns:
        A dict of the two bands and their difference, in rad^2.
    """
    from olb.models.fast import _cn2_layers
    from olb.turbulence.profiles import DEFAULT_HS, get_c2n
    from validation.fast_stone_pointahead.fast_stone_pointahead import (
        _build, _load, _servo_off_params)

    fast, funcs, _spectra = _load()
    scn, _geom = hero_uplink(elevation_deg)
    site = scn.channel.site
    hs = DEFAULT_HS
    cn2 = get_c2n(hs, site.wind_rms_m_s, site.cn2_ground)
    theta_arcsec = float(theta) / ARCSEC

    def run(zmax):
        params = _servo_off_params(elevation_deg, hs, cn2, zmax, theta_arcsec,
                                   grid, L0=L0_M)
        params.update(D_GROUND=APERTURE_M, OBSC_GROUND=OBSCURATION,
                      W0=float(scn.ground.transmitter.waist_m),
                      D_SAT=float(scn.space.aperture_m),
                      H_SAT=float(scn.channel.altitude_m),
                      CN2_TURB=_cn2_layers(cn2, hs))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            sim = _build(fast, params)
            from validation.fast_stone_pointahead.fast_stone_pointahead import \
                split_errors
            return split_errors(sim, funcs)["aniso_corr"], sim

    v3, sim3 = run(FAST_ZMAX_TILT)
    v1, _sim1 = run(FAST_ZMAX_PISTON)
    say(f"    FAST at {theta_arcsec:.2f} arcsec: ZMAX 3 {v3:.5f} rad^2, "
        f"ZMAX 1 {v1:.5f} rad^2, tilt {v3 - v1:.5f} rad^2 "
        f"(df = {float(sim3.freq.main.df):.3f} rad/m)")
    return {"zmax3": float(v3), "zmax1": float(v1), "tilt": float(v3 - v1),
            "df": float(sim3.freq.main.df)}


# ---------------------------------------------------------------------------
# The driver
# ---------------------------------------------------------------------------

def _log_maker():
    """Make the log file of this script, and give its `say` function."""
    path = os.path.join(HERE, "tilt_anisoplanatism.log")
    with open(path, "w", encoding="utf-8"):
        pass

    def say(text=""):
        print(text, flush=True)
        with open(path, "a", encoding="utf-8") as fh:
            fh.write(text + "\n")

    return say, path


def _p1_tiptilt_penalty(elevation_deg, angles):
    """Give the measured tip-tilt-stack mean penalty of p1, in dB.

    It reads `waveoptics_pointahead_p1_results.json` next to this script. The
    value is the CONTEXT column of the dB table: it says how far the extended
    Marechal map of the tilt variance sits from the field penalty.

    THE MATCH IS BY ANGLE, not by index. A p1 campaign holds five angles and
    a smaller store holds fewer, so the reader takes the p1 column whose angle
    is nearest, and it gives None when no column is nearer than 1 nrad.

    Args:
        elevation_deg: the elevation, in deg.
        angles:        the angles of this store, in rad.

    Returns:
        A list with one penalty (or None) for each angle, in dB, or None when
        the p1 record is missing.
    """
    path = os.path.join(HERE, "waveoptics_pointahead_p1_results.json")
    if not os.path.exists(path):
        return None
    try:
        with open(path, "r", encoding="utf-8") as fh:
            got = json.load(fh)
        cell = got["elevations"][str(float(elevation_deg))]
        ref = np.asarray(cell["angles_rad"], dtype=float)
        pen = [float(v) for v in cell["field"]["tiptilt"]["penalty_db"]]
    except (KeyError, ValueError, TypeError):
        return None
    out = []
    for a in angles:
        j = int(np.argmin(np.abs(ref - float(a))))
        out.append(pen[j] if abs(ref[j] - float(a)) < 1e-9 else None)
    return out


def measure(elevation_deg, args, say):
    """Measure one elevation. Give the JSON row."""
    scn, _geom = hero_uplink(elevation_deg)
    camp, warns = campaign_of(elevation_deg, "base", args)
    for w in warns:
        say(f"  sizer warning: {w}")
    if camp.n_stored < 1:
        raise SystemExit(f"the campaign {camp.root_dir} holds no trial.")
    n = int(camp.n_stored if args.n_trials is None
            else min(args.n_trials, camp.n_stored))
    angles = [float(a) for a in camp.point_ahead_rad]
    labels = angle_labels(angles)
    dx = float(camp.grid.pixel_m)
    say(f"  store  : {camp.root_dir}")
    say(f"  grid {camp.grid.n} px, {dx * 1e3:.3f} mm pixel, screen "
        f"{camp.screen_n} px, {camp.plan.z_m.size} screens, patch radius "
        f"{camp.patch.radius_m:.3f} m")
    say(f"  trials : {n} of {camp.n_stored}")
    say(f"  angles : " + ", ".join(f"{lab} = {a / ARCSEC:.3f} arcsec"
                                   for lab, a in zip(labels, angles)))

    t0 = time.time()
    # `screen_phase=True` KEEPS the point-ahead phase planes. The False branch
    # of `Campaign.map_trials` drops BOTH `screen_phase` and
    # `screen_phase_pa`, so this study must not use it.
    values = camp.map_trials(BandTrial(), n_trials=n, fields=True,
                             screen_phase=True, compact=True,
                             workers=args.workers)
    say(f"  the read of {n} trials took {time.time() - t0:.1f} s")

    # values: (n_trials, n_angles, 10). Every column is already a per-trial
    # mask phase variance, so the mean over trials is the reported band.
    mean = values.mean(axis=0)
    share = mean[:, 1:5]
    closure = mean[:, 6] - (share.sum(axis=1) + mean[:, 5])
    rel = np.abs(closure) / np.maximum(mean[:, 6], 1e-30)
    worst = float(np.max(rel))
    say(f"  the band-sum closure holds to {worst:.2e} "
        f"({'PASS' if worst <= CLOSURE_GATE else 'FAIL'}, the gate is "
        f"{CLOSURE_GATE:g})")

    tilt_half = [_boot_half(values[:, i, 1], lambda x: float(x.mean()))
                 for i in range(len(angles))]

    # The references.
    stone25 = [stone_bands(scn, elevation_deg, a, L0_M) for a in angles]
    stone_inf = [stone_bands(scn, elevation_deg, a, np.inf) for a in angles]
    disc = [discrete_stone(camp.plan, APERTURE_M, a, LAM, remove="piston",
                           max_order=1, dx=dx,
                           inner=_inner(1, L0_M))["sigma2"]
            for a in angles]

    fast_rows = None
    if not args.no_fast:
        try:
            fast_rows = [fast_tilt_band(elevation_deg, a, args.fast_grid, say)
                         for a in angles]
        except ImportError as exc:
            say(f"  fast-aosim is not available: {exc}")
        except Exception as exc:                    # noqa: BLE001
            say(f"  the FAST arm failed: {type(exc).__name__}: {exc}")

    # ---- the angle table --------------------------------------------------
    say()
    say(f"  THE TILT ANISOPLANATIC VARIANCE [rad^2], {n} trials")
    say(f"  {'angle':>12s}{'screen':>10s}{'+-2s':>8s}{'slopes':>9s}"
        f"{'gtilt':>9s}{'Stone px':>10s}{'Stone 25':>10s}{'Stone inf':>10s}"
        f"{'FAST':>9s}{'fld/St':>8s}{'FAST/St':>8s}")
    say("  " + "-" * 103)
    rows = []
    for i, lab in enumerate(labels):
        s25 = stone25[i]["tilt"]
        f_v = None if fast_rows is None else fast_rows[i]["tilt"]
        r_fld = mean[i, 1] / s25 if s25 > 0 else float("nan")
        r_fast = (f_v / s25 if (f_v is not None and s25 > 0)
                  else float("nan"))
        say(f"  {lab:>12s}{mean[i, 1]:>10.4f}{tilt_half[i]:>8.4f}"
            f"{mean[i, 8]:>9.4f}{mean[i, 9]:>9.4f}{disc[i]:>10.4f}"
            f"{s25:>10.4f}{stone_inf[i]['tilt']:>10.4f}"
            + (f"{f_v:>9.4f}" if f_v is not None else f"{'-':>9s}")
            + f"{r_fld:>8.3f}"
            + (f"{r_fast:>8.3f}" if f_v is not None else f"{'-':>8s}"))
        rows.append({
            "angle_label": lab, "angle_rad": angles[i],
            "angle_arcsec": angles[i] / ARCSEC,
            "screen_tilt_rad2": float(mean[i, 1]),
            "screen_tilt_diag_rad2": float(mean[i, 7]),
            "screen_tilt_half2sigma_rad2": float(tilt_half[i]),
            "slopes_tilt_rad2": float(mean[i, 8]),
            "gtilt_tilt_rad2": float(mean[i, 9]),
            "discrete_stone_rad2": float(disc[i]),
            "continuous_stone_L025_rad2": float(s25),
            "continuous_stone_Linf_rad2": float(stone_inf[i]["tilt"]),
            "fast_tilt_rad2": (None if f_v is None else float(f_v)),
            "ratio_field_stone": float(r_fld),
            "ratio_fast_stone": (None if f_v is None else float(r_fast)),
            "bands_rad2": {label: float(mean[i, b])
                           for b, (label, _s) in enumerate(BANDS)},
            "remainder_rad2": float(mean[i, 5]),
            "mask_variance_rad2": float(mean[i, 6]),
            "closure_rel": float(rel[i]),
        })
    say("  screen = the fitted tilt band of the summed-screen difference (the "
        "Stone quantity). slopes and gtilt read the PROPAGATED fields.")
    say("  Stone px = the delta-layer sum on the campaign's own plan with the "
        "pixel-rounded shifts, at L0 = 25 m.")
    say("  fld/St and FAST/St divide by the continuous Stone at L0 = 25 m.")
    say()

    # ---- the band table at the geometry angle -----------------------------
    g = 1 if len(angles) > 1 else 0
    say(f"  THE NOLL BANDS of the summed-screen difference at {labels[g]} "
        f"[rad^2]")
    say(f"  {'band':<14s}{'field':>11s}{'diag':>11s}{'Stone 25':>11s}"
        f"{'share':>9s}")
    say("  " + "-" * 56)
    band_rows = []
    for b, (label, _s) in enumerate(BANDS):
        v = float(mean[g, b])
        diag = float(mean[g, 7]) if label == "tilt" else float("nan")
        st = stone25[g].get(label, float("nan"))
        # The piston is OUT of the piston-removed total, so it has no share.
        frac = (float("nan") if label == "piston"
                else (v / mean[g, 6] if mean[g, 6] > 0 else float("nan")))
        say(f"  {label:<14s}{v:>11.5f}"
            + (f"{diag:>11.5f}" if label == "tilt" else f"{'-':>11s}")
            + f"{st:>11.5f}"
            + (f"{'-':>8s} " if label == "piston"
               else f"{frac * 100:>8.2f}%"))
        band_rows.append({"band": label, "field_rad2": v,
                          "stone_L025_rad2": float(st),
                          "share_of_mask_variance": float(frac)})
    say(f"  {'remainder':<14s}{mean[g, 5]:>11.5f}{'-':>11s}"
        f"{stone25[g]['above order 5']:>11.5f}"
        f"{mean[g, 5] / mean[g, 6] * 100:>8.2f}%")
    say(f"  {'TOTAL':<14s}{mean[g, 6]:>11.5f}{'-':>11s}"
        f"{stone25[g]['total']:>11.5f}{100.0:>8.2f}%")
    band_rows.append({"band": "remainder", "field_rad2": float(mean[g, 5]),
                      "stone_L025_rad2": float(stone25[g]["above order 5"]),
                      "share_of_mask_variance":
                          float(mean[g, 5] / mean[g, 6])})
    say("  The piston row carries no Stone value: `remove='piston'` takes the "
        "piston out of the Stone integral by definition.")
    say("  TOTAL is the mask variance of the PISTON-REMOVED difference. The "
        "remainder is what the 21-mode fit leaves.")
    say()

    # ---- the dB context ---------------------------------------------------
    p1 = _p1_tiptilt_penalty(elevation_deg, angles)
    say("  THE dB CONTEXT. The extended Marechal map (10 / ln 10) sigma^2 "
        "next to the MEASURED p1 tip-tilt penalty.")
    say(f"  {'angle':>12s}{'tilt [rad2]':>13s}{'Marechal [dB]':>15s}"
        f"{'p1 field [dB]':>15s}{'ratio':>8s}")
    say("  " + "-" * 65)
    db_rows = []
    for i, lab in enumerate(labels):
        mar = 10.0 / _LN10 * float(mean[i, 1])
        meas = None if p1 is None else p1[i]
        ratio = (mar / meas if (meas is not None and abs(meas) > 1e-9)
                 else float("nan"))
        say(f"  {lab:>12s}{mean[i, 1]:>13.4f}{mar:>15.3f}"
            + (f"{meas:>15.3f}{ratio:>8.2f}" if meas is not None
               else f"{'-':>15s}{'-':>8s}"))
        db_rows.append({"angle_label": lab, "marechal_db": float(mar),
                        "p1_tiptilt_penalty_db": meas})
    say("  The Marechal column maps the TILT variance alone. A tip-tilt stack "
        "pays that tilt, so the two columns answer the same question through")
    say("  two different maps. A ratio far from 1 says the dB map, not the "
        "phase variance, carries the p1 offset.")
    say()

    return {
        "elevation_deg": float(elevation_deg),
        "n_trials": n, "grid_n": int(camp.grid.n), "pixel_m": dx,
        "n_screens": int(camp.plan.z_m.size),
        "patch_radius_m": float(camp.patch.radius_m),
        "aperture_m": APERTURE_M, "n_modes": N_MODES, "L0_m": float(L0_M),
        "closure_worst_rel": worst,
        "angles": rows, "bands_at_geometry_angle": band_rows,
        "geometry_angle_label": labels[g],
        "db_context": db_rows,
        "fast": fast_rows,
        "fast_grid": list(args.fast_grid) if fast_rows else None,
    }


def plot(payload):
    """Draw the tilt variance against the angle, one panel for each elevation."""
    plt = _pyplot()
    if plt is None:
        return []
    els = sorted(payload["elevations"], key=float)
    fig, axes = plt.subplots(1, len(els), figsize=(6.0 * len(els), 5.0),
                             squeeze=False)
    for ax, key in zip(axes[0], els):
        cell = payload["elevations"][key]
        rows = sorted(cell["angles"], key=lambda r: r["angle_arcsec"])
        x = [r["angle_arcsec"] for r in rows]
        ax.plot(x, [r["screen_tilt_rad2"] for r in rows], "o-",
                label="field, summed screens")
        ax.plot(x, [r["slopes_tilt_rad2"] for r in rows], "^-",
                label="field, wrapped slopes")
        ax.plot(x, [r["gtilt_tilt_rad2"] for r in rows], "v-",
                label="field, G-tilt centroid")
        ax.plot(x, [r["continuous_stone_L025_rad2"] for r in rows], "k--",
                label="Stone, L0 = 25 m")
        ax.plot(x, [r["continuous_stone_Linf_rad2"] for r in rows], "k:",
                label="Stone, L0 infinite")
        ax.plot(x, [r["discrete_stone_rad2"] for r in rows], "x--",
                label="Stone, delta layers, pixel rounded")
        if rows[0]["fast_tilt_rad2"] is not None:
            ax.plot(x, [r["fast_tilt_rad2"] for r in rows], "s-.",
                    label="FAST tilt band")
        ax.set_xlabel("point-ahead angle [arcsec]")
        ax.set_ylabel("tilt anisoplanatic variance [rad$^2$]")
        ax.set_title(f"{float(key):.0f} deg")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7)
    fig.suptitle("The TILT anisoplanatism of the point-ahead record "
                 "(D = 0.7 m, L0 = 25 m)", fontsize=10)
    fig.tight_layout()
    path = os.path.join(HERE, "figures")
    os.makedirs(path, exist_ok=True)
    path = os.path.join(path, "tilt_anisoplanatism.png")
    fig.savefig(path, dpi=140)
    plt.close(fig)
    return [path]


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--elevations", nargs="+", type=float,
                    default=[30.0, 20.0],
                    help="the elevations of the stored campaigns [deg]")
    ap.add_argument("--n-trials", type=int, default=None,
                    help="the trials to read. Leave it out for every trial.")
    ap.add_argument("--root", default=None,
                    help="the parent directory of the campaigns. It overrides "
                         f"{ENV_ROOT}.")
    ap.add_argument("--preset", default="standard",
                    help="the sampling preset of the stored campaigns")
    ap.add_argument("--fixed-arcsec", nargs="*", type=float, default=None,
                    help="the fixed point-ahead angles of the stored "
                         "campaigns [arcsec], after the geometry angle")
    ap.add_argument("--block-size", type=int, default=50,
                    help="the trials in one block file of the store")
    ap.add_argument("--fft-backend", default="cupy",
                    choices=["numpy", "scipy", "cupy"],
                    help="the FFT backend of the stored campaigns. This "
                         "script propagates nothing; the value only has to "
                         "match the stored fingerprint.")
    ap.add_argument("--workers", default=None,
                    help="the pool size of the post-hoc read. None runs it in "
                         "this process.")
    ap.add_argument("--no-fast", action="store_true",
                    help="skip the fidelity-1 FAST tilt band")
    ap.add_argument("--fast-npxls", type=int, default=FAST_GRID[0],
                    help="the FAST grid side [px]")
    ap.add_argument("--fast-dx", type=float, default=FAST_GRID[1],
                    help="the FAST grid pitch [m]")
    ap.add_argument("--no-figures", dest="figures", action="store_false",
                    help="skip the figure")
    args = ap.parse_args()
    if args.workers not in (None, "auto"):
        args.workers = int(args.workers)
    args.fast_grid = (int(args.fast_npxls), float(args.fast_dx))
    if args.root:
        os.environ[ENV_ROOT] = os.path.abspath(args.root)

    # The campaign options must match the stored fingerprint. The driver keeps
    # the preset and the fixed angle list as module globals, so a smoke store
    # of another shape reopens through these two switches.
    from validation.waveoptics_pointahead import waveoptics_pointahead as wp
    wp.PRESET = args.preset
    if args.fixed_arcsec is not None:
        wp.FIXED_ARCSEC = tuple(float(a) for a in args.fixed_arcsec)

    say, log_path = _log_maker()
    t0 = time.time()
    say("THE TILT ANISOPLANATISM OF THE POINT-AHEAD RECORD (backlog 2-P4)")
    say("date          : 2026-09-11")
    say("question      : isolate the TILT decorrelation in rad^2, with NO "
        "Marechal map, from the stored planes.")
    say(f"case          : the hero uplink, D = {APERTURE_M} m, "
        f"{LAM * 1e9:.0f} nm, 500 km, L0 = {L0_M:g} m")
    say(f"fit           : {N_MODES} Noll modes over the aperture mask on the "
        "stored crop (Noll, DOI 10.1364/JOSA.66.000207, Table I)")
    say("reference     : Stone et al. (1994), DOI 10.1364/JOSAA.11.000347, "
        "Eq. (29), with the von Karman kernel of Andrews and Phillips")
    say("                DOI 10.1117/3.626196, Ch. 3, Eq. (20)")
    say(f"FAST          : the tilt band is ZMAX {FAST_ZMAX_TILT} minus ZMAX "
        f"{FAST_ZMAX_PISTON}, SERVO OFF, on the explicit grid "
        f"{args.fast_grid[0]} px / {args.fast_grid[1]:g} m")
    say("                (the grid rule of "
        "validation/fast_stone_pointahead/, sized so the side holds L0)")
    say("caveat        : the script propagates nothing. It reads the stored "
        "planes of the `base` campaigns.")
    say()

    payload = {"study": "tilt_anisoplanatism", "date": "2026-09-11",
               "L0_m": float(L0_M), "aperture_m": APERTURE_M,
               "n_modes": N_MODES, "elevations": {}}
    for el in args.elevations:
        say(f"ELEVATION {el:.0f} deg")
        payload["elevations"][str(el)] = measure(el, args, say)

    figs = plot(payload) if args.figures else []
    payload["figures"] = figs
    payload["runtime_s"] = time.time() - t0
    path = os.path.join(HERE, "tilt_anisoplanatism_results.json")
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=1)
    say(f"wrote {path}")
    say(f"wrote {log_path}")
    for f in figs:
        say(f"wrote {f}")
    say(f"runtime {time.time() - t0:.1f} s")


if __name__ == '__main__':
    main()
