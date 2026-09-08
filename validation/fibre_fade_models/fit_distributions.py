"""Fit the received-power distributions of the terrestrial campaigns.

This is step 1 to step 3 of backlog 1-9. It reads the per-trial tables that
`extract_trials.py` wrote from the twelve terrestrial fidelity-2 campaigns
(`validation/terrestrial_campaigns/`), and for each cell it fits the
candidate families to the BUCKET power, to the FIBRE-coupled power and to the
POINT irradiance. It reports the fade quantiles of each family against the
empirical ones, and it gives one verdict for each receiver case: which
family holds with the ANALYTIC parameters (the free route), which holds only
when REFIT to the data (the calibrated route of 1-8), and where NO family
holds.

THE QUANTITIES. Every power is normalised by its own sample mean, so a fade
is measured against the mean power, and the deterministic geometric loss
drops out. The loss is `-10 log10(P / <P>)`, positive dB. The fade at the
exceedance q is the loss that the link exceeds a fraction q of the time.

THE FAMILIES (`olb/turbulence/andrews/distributions.py`, Andrews and Phillips,
2nd ed. (2005), DOI 10.1117/3.626196):

- the lognormal (Ch. 5.7.2): with the ANALYTIC index of the fidelity-0 Term
  (`terrestrial_scintillation_term`, `sigma2_P = A sigma2_I`), with the index
  REFIT to the measured normalised variance (the 1-8 route), and by maximum
  likelihood on `ln P` (a free mean and a free variance);
- the gamma-gamma (Ch. 9.10), by maximum likelihood on (alpha, beta);
- the K distribution (Ch. 9.9.1), by maximum likelihood on alpha;
- the negative exponential, the speckle limit of a coherent sum (no
  parameter; the gamma-gamma at alpha -> inf, beta = 1);
- the lognormal-Rician (Ch. 9, Eq. (133)), by a binned maximum likelihood on
  (r, sigma_z^2), AND with the STREHL MAP of the README Section 3.1, an olb
  HEURISTIC and NOT a result of the book (which states at printed p. 369
  that no map from atmospheric conditions to r and sigma_z^2 is known): the
  coherent part is the Strehl amplitude of the Noll residual, so
  `r = S / (1 - S)` with `S = exp(-sigma_phi^2)` (the extended Marechal
  form, T. S. Ross, DOI 10.1364/AO.48.001812), and `sigma_z^2 = ln(1 +
  sigma2_P)` from the analytic bucket index;
- for a fibre, two more FREE routes built from the shipped fidelity-0 Terms
  (`olb/models/coupling/terrestrial.py`): the SHIPPED chain (the mean-only
  higher-order coupling plus the exponential tip-tilt walk-off fade), and the
  COMPOSITE of the analytic lognormal bucket times that walk-off fade.

THE VERDICT RULE. A family HOLDS on a case when its fade agrees with the
empirical fade inside PASS_P5_DB at the 5 percent exceedance and inside
PASS_P1_DB at the 1 percent exceedance. The 0.5 dB value at 5 percent is the
"notable" rule of the 1-6 certification. The empirical quantiles carry a
bootstrap 16 to 84 percent band, which the tables print.

THE TILT CHECK. The tables also compare the measured per-axis tilt angle
variance of the aperture (from the Noll coefficients a2, a3 that the slope
sensor fitted) against the Noll Zernike tilt `0.182 (D/r0)^(5/3) (lam/D)^2`
per axis (R. J. Noll, DOI 10.1364/JOSA.66.000207, Table IV, the two tilt
modes at 0.448 each) at the plane-wave r0 and at the Gaussian-beam r0 of the
coupling Term, and against the beam-wander arrival tilt that the walk-off Term
reads (`olb.turbulence.angle_of_arrival`).

Run it from the repository root:

    python -m validation.fibre_fade_models.fit_distributions

It writes `fit_results.json`, `fit_distributions.log` and the figures under
`figures/`.
"""

import argparse
import glob
import json
import os
import time
import warnings

import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.optimize import minimize
from scipy.special import i0e
from scipy.stats import norm, skew

from olb.geometry import HorizontalPath
from olb.links.terrestrial import terrestrial_scintillation_term
from olb.models.coupling.terrestrial import (terrestrial_smf_coupling_term,
                                             terrestrial_smf_walkoff_term)
from olb.scenario import TerrestrialChannel, TerrestrialScenario
from olb.terminal import SMF, Terminal, Transmitter
from olb.turbulence.andrews.distributions import (gamma_gamma_pdf,
                                                  gamma_gamma_quantile, k_pdf,
                                                  k_quantile,
                                                  lognormal_rician_pdf)
from olb.turbulence.andrews.scintillation import LOGNORMAL_PDF_LIMIT

HERE = os.path.dirname(os.path.abspath(__file__))
TRIALS_DIR = os.path.join(HERE, "trials")
FIG_DIR = os.path.join(HERE, "figures")
_LN10 = np.log(10.0)

LAM = 1550e-9
WAIST_M = 5e-3
TX_APERTURE_M = 0.10

# The fade exceedance probabilities to report: the 1-6 set plus the median
# (the 50 percent point). The verdict rule reads the 5 and 1 percent points
# only. The MEAN loss in dB is a separate number, printed with each case.
EXCEEDANCE = (0.50, 0.10, 0.05, 0.01)
# The verdict rule (see the module docstring).
PASS_P5_DB = 0.5
PASS_P1_DB = 1.0
# The bootstrap of the empirical quantiles.
N_BOOT = 300
# The Monte Carlo of a sampled free route.
N_DRAW = 200_000
# The slope-sensing aliasing warning of the runner, in rad per pixel, and the
# fraction of trials past it that marks a cell's tilt reading as ALIASED.
STEP_WARN_RAD = 2.8
ALIAS_FRAC = 0.10

# The Noll Zernike tilt: each of the two tilt modes carries 0.448 (D/r0)^(5/3)
# rad^2 (Noll 1976, Table IV, Delta1 - Delta2 = Delta2 - Delta3 = 0.448).
NOLL_TILT_MODE = 0.448


# ---------------------------------------------------------------------------
# The families, on the normalised power x = P / <P>
# ---------------------------------------------------------------------------

def _loss(x_q):
    """Turn a quantile of the normalised power into a fade loss, in dB."""
    return float(-10.0 * np.log10(x_q))


def lognormal_fades(sigma2):
    """The fades of a unit-mean lognormal of normalised variance sigma2."""
    s2 = np.log(1.0 + sigma2)
    s = np.sqrt(s2)
    return {q: float(-10.0 / _LN10 * (-s2 / 2.0 + s * norm.ppf(q)))
            for q in EXCEEDANCE}


def lognormal_mle(x):
    """Fit a Gaussian to ln x, and give the fades and the parameters."""
    lx = np.log(x)
    mu, s = float(lx.mean()), float(lx.std(ddof=1))
    fades = {q: _loss(np.exp(mu + s * norm.ppf(q))) for q in EXCEEDANCE}
    return fades, {"mu": mu, "sigma": s, "sigma2_ln": s * s}


def _nll(pdf, x, params):
    p = np.asarray(pdf(x, *params), dtype=float)
    return float(-np.log(np.maximum(p, 1e-300)).sum())


def gamma_gamma_mle(x):
    """Maximum-likelihood gamma-gamma on the normalised power."""
    index = float(x.var() / x.mean() ** 2)
    # Start from alpha = beta, which solves 2/a + 1/a^2 = index.
    a0 = (1.0 + np.sqrt(1.0 + index)) / index
    res = minimize(lambda t: _nll(gamma_gamma_pdf, x, np.exp(t)),
                   x0=np.log([a0, a0]), method="Nelder-Mead",
                   options={"xatol": 1e-4, "fatol": 1e-3, "maxiter": 400})
    alpha, beta = np.exp(res.x)
    fades = {q: _loss(gamma_gamma_quantile(q, alpha, beta)) for q in EXCEEDANCE}
    return fades, {"alpha": float(alpha), "beta": float(beta),
                   "index": float(1 / alpha + 1 / beta + 1 / (alpha * beta)),
                   "nll": float(res.fun)}


def k_mle(x):
    """Maximum-likelihood K distribution on the normalised power."""
    index = float(x.var() / x.mean() ** 2)
    a0 = 2.0 / max(index - 1.0, 0.05)
    res = minimize(lambda t: _nll(k_pdf, x, [np.exp(t[0])]),
                   x0=[np.log(a0)], method="Nelder-Mead",
                   options={"xatol": 1e-4, "fatol": 1e-3, "maxiter": 200})
    alpha = float(np.exp(res.x[0]))
    fades = {q: _loss(k_quantile(q, alpha)) for q in EXCEEDANCE}
    return fades, {"alpha": alpha, "index": 1.0 + 2.0 / alpha,
                   "nll": float(res.fun)}


def negative_exponential_fades():
    """The fades of the unit-mean negative exponential, x_q = -ln(1 - q)."""
    return {q: _loss(-np.log1p(-q)) for q in EXCEEDANCE}


def _lr_grid(r, sigma_z2, n=600):
    """The numeric CDF of the lognormal-Rician on a log grid."""
    x = np.logspace(-5, 1.5, n)
    p = lognormal_rician_pdf(x, r, sigma_z2)
    c = cumulative_trapezoid(p, x, initial=0.0)
    return x, c / c[-1]


def lognormal_rician_fades(r, sigma_z2):
    """The fades of a lognormal-Rician, from its numeric CDF."""
    x, c = _lr_grid(r, sigma_z2)
    return {q: _loss(np.interp(q, c, x)) for q in EXCEEDANCE}


def lognormal_rician_mle(x, n_bins=50):
    """A binned maximum-likelihood lognormal-Rician on the normalised power.

    The PDF is a quadrature for each point, so the fit bins ln x and it
    reads the density at the bin centres times the bin width.
    """
    lx = np.log(x)
    edges = np.linspace(lx.min() - 1e-9, lx.max() + 1e-9, n_bins + 1)
    counts, _ = np.histogram(lx, edges)
    centres = np.exp(0.5 * (edges[1:] + edges[:-1]))
    widths = np.diff(np.exp(edges))
    keep = counts > 0

    def nll(t):
        # The bounds keep the quadrature finite: r in [1e-3, 1e4] and
        # sigma_z^2 in [1e-3, 10].
        r = float(np.clip(np.exp(t[0]), 1e-3, 1e4))
        sz2 = float(np.clip(np.exp(t[1]), 1e-3, 10.0))
        with np.errstate(all="ignore"):
            p = lognormal_rician_pdf(centres[keep], r, sz2) * widths[keep]
        if not np.all(np.isfinite(p)):
            return 1e30
        return float(-(counts[keep] * np.log(np.maximum(p, 1e-300))).sum())

    index = float(x.var() / x.mean() ** 2)
    # The Rician alone has the index (1 + 2r) / (1 + r)^2. Start there when
    # the index sits under one, else at the speckle limit.
    r0 = max((1.0 - index) / max(index, 1e-3), 1e-2) if index < 1 else 1e-2
    best = None
    for start in ([np.log(r0), np.log(0.1)], [np.log(r0), np.log(0.5)],
                  [np.log(10.0), np.log(index if index > 0 else 0.1)]):
        res = minimize(nll, x0=start, method="Nelder-Mead",
                       options={"xatol": 1e-3, "fatol": 1e-2, "maxiter": 300})
        if best is None or res.fun < best.fun:
            best = res
    r = float(np.clip(np.exp(best.x[0]), 1e-3, 1e4))
    sz2 = float(np.clip(np.exp(best.x[1]), 1e-3, 10.0))
    return (lognormal_rician_fades(r, sz2),
            {"r": float(r), "sigma_z2": float(sz2), "nll": float(best.fun)})


def sampled_fades(draws):
    """The fades of a set of sampled powers, normalised by their mean."""
    loss = -10.0 * np.log10(draws / draws.mean())
    return {q: float(np.quantile(loss, 1.0 - q)) for q in EXCEEDANCE}


def rician_draws(K, n, rng):
    """Draw n unit-mean Rician powers of K factor K.

    The power of a coherent amplitude sqrt(K / (1 + K)) plus a complex
    Gaussian of variance 1 / (1 + K). The mean is one.
    """
    a = np.sqrt(K / (1.0 + K))
    s = np.sqrt(0.5 / (1.0 + K))
    z = (a + s * rng.standard_normal(n)) + 1j * s * rng.standard_normal(n)
    return np.abs(z) ** 2


def empirical(x, rng):
    """The empirical fades, with their bootstrap 16 to 84 percent band."""
    loss = -10.0 * np.log10(x)
    fades = {q: float(np.quantile(loss, 1.0 - q)) for q in EXCEEDANCE}
    boot = np.empty((N_BOOT, len(EXCEEDANCE)))
    for b in range(N_BOOT):
        pick = rng.integers(0, loss.size, loss.size)
        boot[b] = np.quantile(loss[pick], [1.0 - q for q in EXCEEDANCE])
    lo, hi = np.percentile(boot, [16, 84], axis=0)
    band = {q: (float(lo[i]), float(hi[i])) for i, q in enumerate(EXCEEDANCE)}
    return fades, band


# ---------------------------------------------------------------------------
# The analytic feeds of one cell and one aperture
# ---------------------------------------------------------------------------

def scenario(path_m, cn2, D, defocus_m=0.0):
    near = Terminal(aperture_m=TX_APERTURE_M, wavelength_m=LAM,
                    transmitter=Transmitter(waist_m=WAIST_M))
    far = Terminal(aperture_m=D, wavelength_m=LAM,
                   detector=SMF(optimal_focus=True, sensitivity_dbm=-40,
                                defocus_m=defocus_m))
    channel = TerrestrialChannel(path_length_m=path_m,
                                 attenuation_db_per_km=0.0, cn2=cn2)
    return TerrestrialScenario(near=near, far=far, channel=channel)


def analytic_feeds(path_m, cn2, D, dz_curv_m):
    """The fidelity-0 quantities of one cell and one receive aperture."""
    geom = HorizontalPath(path_m)
    out = {}
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        scn = scenario(path_m, cn2, D)
        scint = terrestrial_scintillation_term(scn, geom)
        full = terrestrial_smf_coupling_term(scn, geom, drop_tiptilt=False)
        ho = terrestrial_smf_coupling_term(scn, geom, drop_tiptilt=True)
        walk = terrestrial_smf_walkoff_term(scn, geom)
        scn_t = scenario(path_m, cn2, D, defocus_m=dz_curv_m)
        ho_t = terrestrial_smf_coupling_term(scn_t, geom, drop_tiptilt=True)
        walk_t = terrestrial_smf_walkoff_term(scn_t, geom)
    out["warnings"] = sorted({str(w.message)[:160] for w in caught})
    m = scint.meta
    out["scint"] = {k: float(m[k]) for k in
                    ("sigma2_I", "sigma2_P", "aperture_averaging_factor",
                     "sigma2_R", "Lambda")}
    out["scint"]["rytov_regime"] = m["rytov_regime"]
    out["coupling"] = {
        "r0_gauss_m": float(full.meta["r0_m"]),
        "eta_max": float(full.meta["eta_max"]),
        "eta_flat": float(full.meta["eta_flat_wavefront"]),
        "sigma2_res_full": float(full.meta["sigma2_res"]),
        "sigma2_res_ho": float(ho.meta["sigma2_res"]),
        "eta_full": float(full.meta["eta"]),
        "eta_ho": float(ho.meta["eta"]),
        "eta_ho_tracked": float(ho_t.meta["eta"]),
        "curvature_penalty_db": float(full.meta["curvature_penalty_db"]),
    }
    out["walkoff"] = {
        "mean_db": float(walk.mean_db),
        "mean_db_tracked": float(walk_t.mean_db),
        "sigma2_wander_radial": float(walk.meta["sigma2_wander"]),
        "w_eff_m": float(walk.meta["w_eff_m"]),
        "w_eff_tracked_m": float(walk_t.meta["w_eff_m"]),
        "spot_offset_1sigma_m": float(walk.meta["spot_offset_1sigma_m"]),
    }
    out["_terms"] = {"walk": walk, "walk_t": walk_t}
    return out


# ---------------------------------------------------------------------------
# One receiver case
# ---------------------------------------------------------------------------

def _delta(fades, ref):
    return {q: fades[q] - ref[q] for q in EXCEEDANCE}


def _holds(delta):
    d5, d1 = delta[0.05], delta[0.01]
    if not (np.isfinite(d5) and np.isfinite(d1)):
        return False
    return bool(abs(d5) <= PASS_P5_DB and abs(d1) <= PASS_P1_DB)


def fit_case(name, kind, x, feeds, rng, variant=None, quick=False):
    """Fit every family to one normalised power sample, and judge them."""
    x = np.asarray(x, dtype=float)
    x = x / x.mean()
    emp, band = empirical(x, rng)
    index = float(x.var())
    sk = float(skew(np.log(x)))
    fams = {}          # name -> (fades, params, free)

    # --- the shared families -------------------------------------------
    if kind in ("bucket", "point"):
        s2 = (feeds["scint"]["sigma2_P"] if kind == "bucket"
              else feeds["scint"]["sigma2_I"])
        fams["lognormal (analytic)"] = (lognormal_fades(s2), {"sigma2": s2},
                                        True)
    fams["lognormal (refit index)"] = (lognormal_fades(index),
                                       {"sigma2": index}, False)
    f, p = lognormal_mle(x)
    fams["lognormal (MLE)"] = (f, p, False)
    f, p = gamma_gamma_mle(x)
    fams["gamma-gamma (MLE)"] = (f, p, False)
    f, p = k_mle(x)
    fams["K (MLE)"] = (f, p, False)
    fams["negative exponential"] = (negative_exponential_fades(), {}, True)
    if not quick:
        f, p = lognormal_rician_mle(x)
        fams["lognormal-Rician (MLE)"] = (f, p, False)

    # --- the fibre free routes -----------------------------------------
    if kind == "smf":
        c, w = feeds["coupling"], feeds["walkoff"]
        tracked = "tracked" in variant
        tilt_removed = "tiptilt" in variant
        s2_phi = c["sigma2_res_ho"] if tilt_removed else c["sigma2_res_full"]
        S = float(np.exp(-s2_phi))
        r = S / max(1.0 - S, 1e-12)
        sz2 = float(np.log(1.0 + feeds["scint"]["sigma2_P"]))
        fams["lognormal-Rician (Strehl map, olb heuristic)"] = (
            lognormal_rician_fades(r, sz2),
            {"r": r, "sigma_z2": sz2, "strehl": S, "sigma2_phi": s2_phi},
            True)
        if not tilt_removed:
            term = feeds["_terms"]["walk_t" if tracked else "walk"]
            walk_db = term.sample_db(N_DRAW, rng)
            fams["shipped F0 chain (HO mean + walk-off)"] = (
                sampled_fades(10.0 ** (-walk_db / 10.0)),
                {"walkoff_mean_db": float(term.mean_db)}, True)
            s2b = np.log(1.0 + feeds["scint"]["sigma2_P"])
            ln = np.exp(-s2b / 2.0 + np.sqrt(s2b) * rng.standard_normal(N_DRAW))
            fams["composite LN bucket x walk-off"] = (
                sampled_fades(ln * 10.0 ** (-walk_db / 10.0)), {}, True)

    rows = {}
    for fam, (fades, params, free) in fams.items():
        d = _delta(fades, emp)
        rows[fam] = {"fades_db": {str(q): fades[q] for q in EXCEEDANCE},
                     "delta_db": {str(q): d[q] for q in EXCEEDANCE},
                     "holds": _holds(d), "free": free, "params": params}
    free_ok = [f for f, r in rows.items() if r["holds"] and r["free"]]
    refit_ok = [f for f, r in rows.items() if r["holds"] and not r["free"]]
    route = "FREE" if free_ok else ("REFIT" if refit_ok else "NONE")
    best = min(rows, key=lambda f: max(abs(v) for v in
                                       rows[f]["delta_db"].values()))
    return {
        "case": name, "kind": kind, "variant": variant, "n": int(x.size),
        "index": index, "skew_ln": sk,
        "mean_loss_db": float((-10.0 * np.log10(x)).mean()),
        "empirical_db": {str(q): emp[q] for q in EXCEEDANCE},
        "band_db": {str(q): band[q] for q in EXCEEDANCE},
        "families": rows, "route": route, "free_ok": free_ok,
        "refit_ok": refit_ok, "best": best,
        "_x": x,
    }


# ---------------------------------------------------------------------------
# One cell
# ---------------------------------------------------------------------------

def load_tables():
    out = {}
    for path in sorted(glob.glob(os.path.join(TRIALS_DIR, "trials_*.npz"))):
        z = np.load(path)
        meta = json.loads(str(z["meta"]))
        cols = [str(c) for c in z["columns"]]
        v = z["values"]
        out[meta["tag"]] = ({c: v[:, i] for i, c in enumerate(cols)}, meta)
    return out


def tilt_check(cols, D, feeds, cell):
    """Compare the measured tilt angle variance with the analytic tilts."""
    tag = "10cm" if D == 0.10 else "5cm"
    a2, a3 = cols[f"a2_{tag}"], cols[f"a3_{tag}"]
    # A Noll coefficient a (rad) of the tilt mode Z = 2 rho cos(theta) is the
    # phase slope 2 a / R, so the angle is lam a / (pi R) = 2 lam a / (pi D).
    scale = (2.0 * LAM / (np.pi * D)) ** 2
    var_axis = float(0.5 * (a2.var() + a3.var()) * scale)
    lam_D2 = (LAM / D) ** 2

    def noll(r0):
        return float(NOLL_TILT_MODE * (D / r0) ** (5.0 / 3.0) * scale
                     / 1.0)                        # per axis, one mode

    r0_plane = float(cell["r0_m"])
    r0_gauss = feeds["coupling"]["r0_gauss_m"]
    return {
        "measured_axis_rad2": var_axis,
        "noll_r0_plane_rad2": noll(r0_plane),
        "noll_r0_gauss_rad2": noll(r0_gauss),
        "wander_axis_rad2": 0.5 * feeds["walkoff"]["sigma2_wander_radial"],
        "gtilt_r0_plane_rad2": float(0.174 * (D / r0_plane) ** (5.0 / 3.0)
                                     * lam_D2),
        "r0_plane_m": r0_plane, "r0_gauss_m": r0_gauss,
        "D_over_r0_plane": D / r0_plane, "D_over_r0_gauss": D / r0_gauss,
        "mean_tilt_rad": float(np.sqrt(scale) * np.hypot(a2, a3).mean()),
    }


def analyse_cell(tag, cols, meta, rng, quick, say):
    cell = meta["cell"]
    path_m, cn2 = float(cell["path_length_m"]), float(cell["cn2_m_m23"])
    dz = {float(k): v for k, v in meta["dz_curv_m"].items()}
    feeds = {D: analytic_feeds(path_m, cn2, D, dz[D]) for D in (0.10, 0.05)}
    steps = cols["max_step_10cm"]
    out = {"tag": tag, "path_m": path_m, "cn2": cn2, "preset": cell["preset"],
           "sigma2_R": float(cell["sigma2_R_plane"]),
           "n_trials": int(meta["n_trials"]),
           "grid_n": int(cell["grid"]["n"]), "n_screens": int(cell["screens"]["n"]),
           "pixel_m": float(cell["grid"]["pixel_m"]),
           "eta_fill": {a["aperture_m"]: a["eta_fill"] for a in cell["apertures"]},
           "slope_step_over_warn_frac": float((steps > STEP_WARN_RAD).mean()),
           "crosscheck": meta["crosscheck"],
           "feeds": {f"{D:.2f}": {k: v for k, v in f.items() if k != "_terms"}
                     for D, f in feeds.items()},
           "tilt": {f"{D:.2f}": tilt_check(cols, D, feeds[D], cell)
                    for D in (0.10, 0.05)},
           "cases": []}
    # The receiver cases. A fibre power is the bucket power times eta.
    b10, b5 = cols["bucket_10cm"], cols["bucket_5cm"]
    cases = [
        ("bucket 10 cm", "bucket", 0.10, b10, None),
        ("bucket 5 cm", "bucket", 0.05, b5, None),
        ("point", "point", 0.10, cols["point_irradiance"], None),
        ("SMF 10 cm untracked", "smf", 0.10, b10 * cols["smf10_untracked"],
         "untracked"),
        ("SMF 10 cm tracked", "smf", 0.10, b10 * cols["smf10_tracked"],
         "tracked"),
        ("SMF 10 cm tilt removed", "smf", 0.10, b10 * cols["smf10_tiptilt"],
         "tiptilt"),
        ("SMF 10 cm tilt removed tracked", "smf", 0.10,
         b10 * cols["smf10_tiptilt_tracked"], "tiptilt_tracked"),
        ("SMF 5 cm untracked", "smf", 0.05, b5 * cols["smf5_untracked"],
         "untracked"),
        ("SMF 5 cm tracked", "smf", 0.05, b5 * cols["smf5_tracked"], "tracked"),
        ("SMF 5 cm tilt removed", "smf", 0.05, b5 * cols["smf5_tiptilt"],
         "tiptilt"),
        ("MMF 10 cm", "mmf", 0.10, b10 * cols["mmf10"], None),
        ("SMF 10 cm eta only", "eta", 0.10, cols["smf10_untracked"],
         "untracked"),
    ]
    aliased = out["slope_step_over_warn_frac"] > ALIAS_FRAC
    for name, kind, D, x, variant in cases:
        r = fit_case(name, kind, x, feeds[D], rng, variant=variant, quick=quick)
        r["aperture_m"] = D
        r["aliased"] = bool(aliased and variant is not None
                            and "tiptilt" in variant)
        if r["aliased"]:
            r["route"] = "ALIASED"
        r["mean_eta"] = (float(np.mean(cols["smf10_untracked"]))
                         if kind == "eta" else None)
        out["cases"].append(r)
    return out


# ---------------------------------------------------------------------------
# The report
# ---------------------------------------------------------------------------

def _fmt_row(name, fades, delta=None, holds=None, band=None):
    s = f"    {name:<40s}"
    for q in EXCEEDANCE:            # the columns are p50, p10, p5, p1
        s += f"{fades[str(q)] if isinstance(fades, dict) else fades[q]:8.2f}"
    if band is not None:
        s += "   band p5 [{:.2f}, {:.2f}] p1 [{:.2f}, {:.2f}]".format(
            *band["0.05"], *band["0.01"])
    if delta is not None:
        s += "   d5 {:+.2f} d1 {:+.2f} {}".format(
            delta["0.05"], delta["0.01"], "HOLDS" if holds else "")
    return s


def report_cell(res, say):
    say(f"\n=== {res['tag']}: L = {res['path_m'] / 1e3:g} km, Cn2 = {res['cn2']:.0e}, "
        f"{res['preset']}, sigma_R^2 = {res['sigma2_R']:.2f}, "
        f"{res['n_screens']} screens, {res['grid_n']} px, "
        f"{res['n_trials']} trials; slope step > {STEP_WARN_RAD} rad in "
        f"{res['slope_step_over_warn_frac'] * 100:.1f} % of trials")
    for D in ("0.10", "0.05"):
        f, t = res["feeds"][D], res["tilt"][D]
        say(f"  D = {D} m: D/r0 plane {t['D_over_r0_plane']:.2f}, gauss "
            f"{t['D_over_r0_gauss']:.2f}; analytic sigma2_I {f['scint']['sigma2_I']:.3f}, "
            f"A {f['scint']['aperture_averaging_factor']:.3f}, sigma2_P "
            f"{f['scint']['sigma2_P']:.4f}; eta_max {f['coupling']['eta_max']:.3f}, "
            f"sigma2_res full {f['coupling']['sigma2_res_full']:.2f}, HO "
            f"{f['coupling']['sigma2_res_ho']:.2f}; eta_fill "
            f"{res['eta_fill'][float(D)]:.2f}")
        say(f"    tilt ratios: measured / Noll(r0 gauss) "
            f"{t['measured_axis_rad2'] / t['noll_r0_gauss_rad2']:.2f}, "
            f"measured / wander model "
            f"{t['measured_axis_rad2'] / t['wander_axis_rad2']:.2f}"
            + ("  [ALIASED slope sensing]"
               if res['slope_step_over_warn_frac'] > ALIAS_FRAC else ""))
        say(f"    tilt per axis [urad^2]: measured {t['measured_axis_rad2'] * 1e12:.2f}, "
            f"Noll(r0 plane) {t['noll_r0_plane_rad2'] * 1e12:.2f}, "
            f"Noll(r0 gauss) {t['noll_r0_gauss_rad2'] * 1e12:.2f}, "
            f"G-tilt(r0 plane) {t['gtilt_r0_plane_rad2'] * 1e12:.2f}, "
            f"wander model {t['wander_axis_rad2'] * 1e12:.2f}")
    # The tilt share of the fibre fade: the untracked fade against the fade
    # with the tilt removed, at each aperture.
    by = {c["case"]: c for c in res["cases"]}
    for D in ("10 cm", "5 cm"):
        u, t = by[f"SMF {D} untracked"], by[f"SMF {D} tilt removed"]
        say(f"  tilt share of the SMF {D} fade: p5 {u['empirical_db']['0.05']:.2f} "
            f"-> {t['empirical_db']['0.05']:.2f} dB with the tilt removed "
            f"({u['empirical_db']['0.05'] - t['empirical_db']['0.05']:.2f} dB is "
            f"tilt); p1 {u['empirical_db']['0.01']:.2f} -> "
            f"{t['empirical_db']['0.01']:.2f} dB "
            f"({u['empirical_db']['0.01'] - t['empirical_db']['0.01']:.2f} dB)")
    for c in res["cases"]:
        say(f"  -- {c['case']}: index {c['index']:.4f}, skew lnP {c['skew_ln']:+.2f}, "
            f"mean loss {c['mean_loss_db']:.2f} dB, median "
            f"{c['empirical_db']['0.5']:.2f} dB  ->  route {c['route']} "
            f"(best {c['best']})")
        say(f"    {'':<40s}" + "".join(f"{'p' + str(int(q * 100)):>8s}"
                                        for q in EXCEEDANCE))
        say(_fmt_row("empirical", c["empirical_db"], band=c["band_db"]))
        for fam, r in c["families"].items():
            say(_fmt_row(("* " if r["free"] else "  ") + fam, r["fades_db"],
                         r["delta_db"], r["holds"]))


def summary_table(results, say):
    say("\n=== SUMMARY: route for each receiver case (FREE = analytic parameters "
        f"hold; REFIT = a refit family holds; NONE). Rule: |d5| <= {PASS_P5_DB} dB "
        f"and |d1| <= {PASS_P1_DB} dB.")
    names = [c["case"] for c in results[0]["cases"]]
    head = f"{'cell':<28s}" + "".join(f"{n[:14]:>16s}" for n in names)
    say(head)
    for res in results:
        row = f"{res['tag']:<28s}"
        for c in res["cases"]:
            row += f"{c['route']:>16s}"
        say(row)
    say("\n=== best family and its worst |delta| over p10/p5/p1 [dB]")
    for res in results:
        row = f"{res['tag']:<28s}"
        for c in res["cases"]:
            worst = max(abs(v) for v in c["families"][c["best"]]["delta_db"].values())
            row += f"{c['best'][:11]:>12s}{worst:4.1f}"
        say(row)


# ---------------------------------------------------------------------------
# The figures
# ---------------------------------------------------------------------------

def _cdf_panel(ax, c, fams):
    x = np.sort(c["_x"])
    loss = -10.0 * np.log10(x)[::-1]
    p = np.arange(1, x.size + 1) / x.size
    ax.semilogy(loss, p[::-1], "k-", lw=1.5, label="fidelity 2")
    for fam, style in fams:
        r = c["families"].get(fam)
        if r is None:
            continue
        ax.semilogy([r["fades_db"][str(v)] for v in EXCEEDANCE],
                    list(EXCEEDANCE), style, ms=7, label=fam)
    # The mean loss in dB, a separate number from the median (the 50 percent
    # point of the curve).
    ax.axvline(c["mean_loss_db"], color="0.5", ls="--", lw=1)
    ax.text(0.02, 0.04, f"mean {c['mean_loss_db']:.2f} dB, median "
            f"{c['empirical_db']['0.5']:.2f} dB", transform=ax.transAxes,
            fontsize=8, color="0.3")
    ax.set_ylim(5e-3, 1.0)
    ax.grid(True, which="both", alpha=0.3)


def figures(results, quick):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    os.makedirs(FIG_DIR, exist_ok=True)
    panels = {
        "bucket_10cm": ("bucket 10 cm", [("lognormal (analytic)", "bs"),
                                         ("lognormal (refit index)", "g^"),
                                         ("gamma-gamma (MLE)", "ro"),
                                         ("K (MLE)", "mv")]),
        "smf_10cm_untracked": ("SMF 10 cm untracked", [
            ("shipped F0 chain (HO mean + walk-off)", "bs"),
            ("composite LN bucket x walk-off", "c^"),
            ("lognormal-Rician (Strehl map, olb heuristic)", "gD"),
            ("lognormal-Rician (MLE)", "m*"),
            ("gamma-gamma (MLE)", "ro"),
            ("lognormal (MLE)", "y>"),
            ("negative exponential", "kx")]),
        "smf_10cm_tilt_removed": ("SMF 10 cm tilt removed", [
            ("lognormal-Rician (Strehl map, olb heuristic)", "gD"),
            ("lognormal-Rician (MLE)", "m*"),
            ("gamma-gamma (MLE)", "ro"),
            ("lognormal (MLE)", "y>"),
            ("negative exponential", "kx")]),
    }
    for fname, (case, fams) in panels.items():
        fig, axes = plt.subplots(3, 4, figsize=(16, 10), sharey=True)
        for ax, res in zip(axes.ravel(), results):
            c = next(c for c in res["cases"] if c["case"] == case)
            _cdf_panel(ax, c, fams)
            ax.set_title(f"{res['path_m'] / 1e3:g} km, Cn2 {res['cn2']:.0e}, "
                         f"{res['preset']}\nsigma_R^2 {res['sigma2_R']:.2f}, "
                         f"route {c['route']}", fontsize=9)
            ax.set_xlabel("fade below the mean [dB]")
        axes[0, 0].set_ylabel("exceedance probability")
        axes[0, 0].legend(fontsize=7, loc="lower left")
        fig.suptitle(f"{case}: fidelity-2 exceedance curve against the "
                     "family fades at 10, 5 and 1 percent")
        fig.tight_layout()
        fig.savefig(os.path.join(FIG_DIR, f"cdf_{fname}.png"), dpi=110)
        plt.close(fig)


# ---------------------------------------------------------------------------

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--quick", action="store_true",
                    help="skip the lognormal-Rician MLE")
    ap.add_argument("--cells", nargs="*", default=None,
                    help="campaign tags to fit (default: every table)")
    ap.add_argument("--seed", type=int, default=20260908)
    args = ap.parse_args(argv)
    log = open(os.path.join(HERE, "fit_distributions.log"), "w")

    def say(text):
        print(text, flush=True)
        log.write(text + "\n")

    tables = load_tables()
    if args.cells:
        tables = {k: v for k, v in tables.items() if k in set(args.cells)}
    if not tables:
        raise SystemExit(f"no table under {TRIALS_DIR}")
    rng = np.random.default_rng(args.seed)
    results = []
    t0 = time.perf_counter()
    order = sorted(tables, key=lambda k: (tables[k][1]["cell"]["path_length_m"],
                                          tables[k][1]["cell"]["cn2_m_m23"],
                                          tables[k][1]["cell"]["preset"]))
    for tag in order:
        cols, meta = tables[tag]
        res = analyse_cell(tag, cols, meta, rng, args.quick, say)
        report_cell(res, say)
        results.append(res)
    summary_table(results, say)
    say(f"\nfit time {time.perf_counter() - t0:.0f} s")
    if len(results) == 12:
        figures(results, args.quick)
    with open(os.path.join(HERE, "fit_results.json"), "w") as fh:
        json.dump([{k: v for k, v in r.items()}
                   for r in _strip(results)], fh, indent=1)


def _strip(results):
    out = []
    for r in results:
        rr = dict(r)
        rr["cases"] = [{k: v for k, v in c.items() if k != "_x"}
                       for c in r["cases"]]
        out.append(rr)
    return out


if __name__ == "__main__":
    main()
