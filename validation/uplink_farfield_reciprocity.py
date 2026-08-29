'''
Mode-matched fidelity-1 versus fidelity-2 uplink scintillation, through the
reciprocity FAR-FIELD MAP.

This is the clean experiment of validation/UPLINK_SIGMA2I_INVESTIGATION.md.

THE IDEA. The fidelity-2 overlap eta_turb is the uplink flux at ONE satellite
point (Shapiro reciprocity, DOI 10.1364/JOSA.61.000492). A satellite point at
transverse offset x differs from the on-axis point by a TILT on the downlink
input: the flux at offset x is

    o(x) = | SUM  E_down(r) * conj(psi_tx(r)) * exp(-i k x . r / L) |^2 .

So ONE FFT of E_down * conj(psi_tx) gives the uplink flux at EVERY satellite
offset at once: the full instantaneous far-field spot. From that map the
script measures, per snapshot:

  - the on-axis flux (identical to the runner's eta_turb),
  - the instantaneous beam-centre offset (the beam wander beta),
  - the flux AT the instantaneous centre (the "tracked" flux),

and across snapshots:

  - sigma2_I on axis            -> against the Dios fidelity-1 var/mean^2,
  - the wander variance <beta^2> -> against Dios/Belmonte 2.07 (coupled_flux
                                    .beam_wander_variance) AND Andrews 7.25,
  - the long/short-term widths   -> against Dios w_lt / w_st,
  - the beam-frame index sigma2_I(r) -> against the Dios on-axis + off-axis
                                    indices (Dios Eqs. (16) and (20)).

So each INGREDIENT of the fidelity-1 model is measured on its own, not only
the one headline number. The comparison is mode-matched: fidelity 1 launches
a pure Gaussian of waist w0, and psi_tx here is the same Gaussian behind a
1.0 m launch aperture, so the clip is negligible (< 2e-7 of the power).

Sources: Dios et al. 2004, DOI 10.1364/AO.43.003866 (the fidelity-1 model);
Shapiro 1971, DOI 10.1364/JOSA.61.000492 (reciprocity); Andrews and Phillips
2nd ed., DOI 10.1117/3.626196, Ch. 12 (the analytic third leg).

Run from the repo root:
    python -m validation.uplink_farfield_reciprocity
'''

import json
import os
import threading
import time
import warnings
from concurrent.futures import ThreadPoolExecutor

import numpy as np
from scipy.ndimage import gaussian_filter
from scipy.signal import ZoomFFT

from olb.beam import gaussz, zR
from olb.geometry import CircularOrbit
from olb.scenario import Channel, Site, SpaceScenario
from olb.terminal import Terminal, Transmitter
from olb.turbulence.andrews.beam import beam_params
from olb.turbulence.andrews.paths import uplink_scintillation_index
from olb.turbulence.andrews.wander import beam_wander_variance_slant
from olb.turbulence.coupled_flux import (beam_wander_variance,
                                         off_axis_scintillation_index,
                                         on_axis_scintillation_index)
from olb.turbulence.profiles import DEFAULT_HS, default_cn2_profile
from olb.turbulence.uplink_flux import _flux_result
from olb.waveoptics.field import Begin
from olb.waveoptics.grid import GridSpec
from olb.waveoptics.turbulence.run import _ground_transmit_mode, _screen_seed
from olb.waveoptics.turbulence.sampling import PRESETS, turbulent_grid
from olb.waveoptics.turbulence.screens import phase_screen
from olb.waveoptics.turbulence.splitstep import (split_step,
                                                 super_gaussian_boundary)

LAM = 1550e-9
ALT_M = 600e3
ELEV_DEG = 60.0
APERTURE_M = 1.0          # launch DIAMETER. The clip on w0 <= 0.18 is < 2e-7.
HS = DEFAULT_HS
SEED = 2026
N_TRIALS = 250            # fidelity-2 snapshots per case
N_TRIALS_VARIANT = 120    # snapshots per convergence variant
N_SAMPLES_F1 = 40000      # fidelity-1 Monte-Carlo samples
HALF_WINDOW_M = 16.0      # half-side of the far-field map [m]
M_PTS = 384               # far-field points per axis (zoom-FFT output)
N_R_BINS = 30             # radial bins of the beam-frame index profile

# The cases: (name, w0 [m], Cn2 scale on the default HV profile). Scale 1.0 is
# the full site profile; 0.3 keeps fidelity 1 weak-valid for the filled beam.
CASES = [
    ("underfilled_hv", 0.06, 1.0),
    ("filled_hv", 0.18, 1.0),
    ("underfilled_0p3", 0.06, 0.3),
    ("filled_0p3", 0.18, 0.3),
]
# The convergence variants run on this case only (the clean disagreement
# point of the investigation note).
VARIANT_CASE = "filled_0p3"

GEOM = CircularOrbit(ALT_M, elevation_deg=[ELEV_DEG])
L_SLANT = float(np.asarray(GEOM.slant_range_m, dtype=float)[0])
AIRMASS = 1.0 / np.sin(np.radians(ELEV_DEG))


def make_scenario(w0):
    '''One uplink scenario with a pure-Gaussian launch behind a wide pupil.'''
    ground = Terminal(aperture_m=APERTURE_M, wavelength_m=LAM,
                      transmitter=Transmitter(waist_m=w0))
    space = Terminal(aperture_m=0.30, wavelength_m=LAM)
    return SpaceScenario(ground=ground, space=space, direction="uplink",
                         channel=Channel(site=Site(), altitude_m=ALT_M))


def beam_centre(m, xs, w_free, px):
    '''Find the instantaneous beam centre on one far-field map crop.

    Smooth the speckle away, take the peak, then refine with a windowed
    centroid (radius 1.5 w_free, two passes).
    '''
    sm = gaussian_filter(m, sigma=max(w_free / 3.0 / px, 1.0))
    iy, ix = np.unravel_index(int(np.argmax(sm)), sm.shape)
    cx, cy = xs[ix], xs[iy]
    X, Y = np.meshgrid(xs, xs)
    for _ in range(2):
        w = m * (((X - cx) ** 2 + (Y - cy) ** 2) < (1.5 * w_free) ** 2)
        tot = w.sum()
        if tot <= 0:
            break
        cx = float((w * X).sum() / tot)
        cy = float((w * Y).sum() / tot)
    return cx, cy


def second_moment_width(m, xs, r_window):
    '''Give w of exp(-2 r^2 / w^2) from <r^2> = w^2/2 inside r < r_window.'''
    X, Y = np.meshgrid(xs, xs)
    r2 = X ** 2 + Y ** 2
    w = m * (r2 < r_window ** 2)
    tot = w.sum()
    return float(np.sqrt(2.0 * (w * r2).sum() / tot)) if tot > 0 else np.nan


def run_fid2(scn, cn2_zen, w0, n_trials, preset="standard", widen=1.0,
             subharmonics=True, seed=SEED, label="", workers=None):
    '''Run n_trials turbulent snapshots and read the far-field maps.

    workers caps the thread count. Each worker holds one full screen stack
    (n_screens * n^2 floats), so a large grid needs a LOW cap to keep the
    memory bounded.
    '''
    p = PRESETS[preset]
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        grid, plan, report = turbulent_grid(scn, GEOM, preset=p, hs=HS,
                                            cn2_profile=cn2_zen)
    if widen != 1.0:
        # Grow the side at a FIXED pixel. The screens then carry lower
        # spatial frequencies, so a truncated wander tilt shows up here.
        grid = GridSpec(size_m=grid.size_m * widen, n=int(round(grid.n * widen)),
                        scaled=False)
    n, dx = grid.n, grid.pixel_m
    mask = super_gaussian_boundary(n, p.boundary_width_frac)
    psi_tx = _ground_transmit_mode(scn.ground, grid)
    n_scr = int(plan.z_m.size)

    # The overlap integrand E * conj(psi_tx) is EXACTLY zero outside the
    # launch aperture (the clip in _ground_transmit_mode), so cut the arrays
    # to the aperture support before the transform. This is lossless.
    x_grid = (np.arange(n) - n / 2) * dx
    keep = np.abs(x_grid) <= (APERTURE_M / 2 + 2 * dx)
    i0, i1 = int(np.argmax(keep)), int(n - np.argmax(keep[::-1]))
    n_cut = i1 - i0
    psi_cut = psi_tx[i0:i1, i0:i1]

    # The far-field map through a 2-D zoom FFT (chirp-z): the DFT on a fine
    # frequency grid over the narrow band that maps to the +-HALF_WINDOW_M
    # satellite window (offset x = lambda * f * L). The half-open frequency
    # grid puts f = 0 exactly at the index M_PTS//2. A physical origin shift
    # only adds a unit phase, so |G|^2 does not depend on it.
    f_max = HALF_WINDOW_M / (LAM * L_SLANT)
    zf = ZoomFFT(n_cut, [-f_max, f_max], M_PTS, fs=1.0 / dx, endpoint=False)
    xs = (-f_max + np.arange(M_PTS) * 2 * f_max / M_PTS) * LAM * L_SLANT
    px = float(xs[1] - xs[0])
    half = M_PTS // 2
    assert abs(xs[half]) < 1e-9 * px

    def farfield(field):
        A = field[i0:i1, i0:i1] * np.conj(psi_cut)
        return np.abs(zf(zf(A, axis=1), axis=0)) ** 2

    flat = np.zeros((n, n))
    F_vac = split_step(Begin(grid.size_m, LAM, n), plan.z_m, [flat] * n_scr,
                       plan.z_total_m, boundary=mask)
    M_vac = farfield(F_vac.field)
    o_vac = float(M_vac[half, half])
    # Sanity: the f = 0 bin of the zoom FFT is the plain overlap sum.
    o_direct = float(np.abs((F_vac.field * np.conj(psi_tx)).sum()) ** 2)
    assert abs(o_vac / o_direct - 1.0) < 1e-9, (o_vac, o_direct)
    w_free = float(gaussz(w0, L_SLANT, LAM))
    w_vac_meas = second_moment_width(M_vac / o_vac, xs, 2.0 * w_free)

    n_crop = xs.size
    acc = {"m": np.zeros((n_crop, n_crop)), "t": np.zeros((n_crop, n_crop)),
           "t2": np.zeros((n_crop, n_crop))}
    lock = threading.Lock()

    def one(k):
        # float32 halves the stack memory. The phase precision loss is far
        # below the statistical noise of the screens.
        stack = [phase_screen(plan.r0_m[j], n, dx,
                              seed=_screen_seed(seed, k, j),
                              subharmonics=subharmonics).astype(np.float32)
                 for j in range(n_scr)]
        F_rx = split_step(Begin(grid.size_m, LAM, n), plan.z_m, stack,
                          plan.z_total_m, boundary=mask)
        m = farfield(F_rx.field) / o_vac
        eta_on = float(m[half, half])
        cx, cy = beam_centre(m, xs, w_free, px)
        # The tracked flux: the map value at the instantaneous centre.
        jx = int(round(cx / px)) + half
        jy = int(round(cy / px)) + half
        track = float(m[jy, jx])
        # Re-centre the map for the beam-frame accumulators.
        t = np.roll(m, (-int(round(cy / px)), -int(round(cx / px))),
                    axis=(0, 1))
        with lock:
            acc["m"] += m
            acc["t"] += t
            acc["t2"] += t * t
        return eta_on, cx, cy, track

    t0 = time.time()
    if workers is None:
        workers = min(max((os.cpu_count() or 4) - 2, 2), 8)
    with ThreadPoolExecutor(max_workers=workers) as ex:
        rows = list(ex.map(one, range(n_trials)))
    wall = time.time() - t0

    eta = np.array([r[0] for r in rows])
    cx = np.array([r[1] for r in rows])
    cy = np.array([r[2] for r in rows])
    track = np.array([r[3] for r in rows])

    mean_m = acc["m"] / n_trials
    mean_t = acc["t"] / n_trials
    var_t = acc["t2"] / n_trials - mean_t ** 2

    # Beam-frame radial index profile: var/mean^2 in radial bins.
    X, Y = np.meshgrid(xs, xs)
    r_map = np.sqrt(X ** 2 + Y ** 2)
    r_edges = np.linspace(0.0, 3.0 * w_free, N_R_BINS + 1)
    prof_r, prof_idx, prof_mean = [], [], []
    for a, b in zip(r_edges[:-1], r_edges[1:]):
        pick = (r_map >= a) & (r_map < b)
        if pick.sum() < 4:
            continue
        mu = float(mean_t[pick].mean())
        vv = float(var_t[pick].mean())
        prof_r.append(0.5 * (a + b))
        prof_mean.append(mu)
        prof_idx.append(vv / mu ** 2 if mu > 0 else np.nan)

    wander_var = float(cx.var() + cy.var())      # radial (two-axis) <beta^2>
    return {
        "label": label or preset, "preset": preset, "widen": widen,
        "subharmonics": subharmonics, "n_trials": n_trials,
        "grid_n": n, "grid_size_m": float(grid.size_m), "n_screens": n_scr,
        "r0_total_m": float(plan.r0_total_m), "px_sat_m": px,
        "wall_s": wall, "w_free_m": w_free, "w_vac_meas_m": w_vac_meas,
        "mean_eta": float(eta.mean()),
        "eta_samples": eta.tolist(),
        "sigma2_I_onaxis": float(eta.var() / eta.mean() ** 2),
        "sigma2_I_tracked": float(track.var() / track.mean() ** 2),
        "wander_var_m2": wander_var,
        "w_lt_meas_m": second_moment_width(mean_m, xs, 2.5 * w_free),
        "w_st_meas_m": second_moment_width(mean_t, xs, 2.5 * w_free),
        "profile_r_m": prof_r, "profile_sigma2": prof_idx,
        "profile_mean": prof_mean,
        "sizer_warnings": [str(w.message) for w in caught],
    }


def run_fid1(w0, cn2_zen):
    '''Run the Dios coupled-flux MC and decompose its ingredients.'''
    np.random.seed(SEED)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        r = _flux_result(w0, ELEV_DEG, L_SLANT, LAM, HS, cn2_zen, 1.7e-14,
                         N_SAMPLES_F1, 1)
    Is = r["Is_summed"]
    k = 2 * np.pi / LAM
    cn2_slant = np.asarray(cn2_zen, dtype=float) * AIRMASS
    wL = float(gaussz(w0, L_SLANT, LAM))
    Z0 = float(zR(w0, LAM))
    beta2 = float(beam_wander_variance(L_SLANT, cn2_slant,
                                       gaussz(w0, HS, LAM), HS))
    s2_on = float(on_axis_scintillation_index(L_SLANT, k, wL, Z0,
                                              cn2_slant, HS))
    # Wander-only index of I = exp(-2 beta^2 / w_st^2) with beta^2 exponential
    # of mean <beta^2>: E[I^m] = 1 / (1 + m s <beta^2>), s = 2/w_st^2.
    s = 2.0 / r["w_st"] ** 2
    wander_only = (1 + s * beta2) ** 2 / (1 + 2 * s * beta2) - 1
    rs = np.linspace(0.0, 3.0 * wL, N_R_BINS)
    s2_off = [float(off_axis_scintillation_index(L_SLANT, k, wL, cn2_slant,
                                                 HS, rr)) for rr in rs]
    return {
        "sigma2_I": float(Is.var() / Is.mean() ** 2),
        "mean_flux": float(Is.mean()),
        "sigma2_x_mean": r["sigma2_x_mean"],
        "weak_valid": r["weak_fluctuation_valid"],
        "w_st_m": r["w_st"], "w_lt_m": r["w_lt"], "r0s_m": r["r0s"],
        "w_free_m": float(gaussz(w0, L_SLANT, LAM)),
        "beta2_dios_m2": beta2, "sigma2_on_axis": s2_on,
        "sigma2_wander_only": float(wander_only),
        "profile_r_m": rs.tolist(),
        "profile_sigma2_on_plus_off": (s2_on + np.asarray(s2_off)).tolist(),
        "warnings": [str(w.message) for w in caught],
    }


def run_andrews(w0, cn2_zen):
    '''The analytic third leg: Andrews Ch. 12 uplink index, and his wander.'''
    bp = beam_params(w0, LAM, L_SLANT)
    out = {}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for name, tracked in (("untracked", False), ("tracked", True)):
            try:
                out[f"sigma2_I_{name}"] = float(uplink_scintillation_index(
                    HS, cn2_zen, LAM, ELEV_DEG, bp, altitude_m=ALT_M,
                    r=0.0, tracked=tracked))
            except Exception as exc:
                out[f"sigma2_I_{name}"] = f"error: {exc}"
        try:
            out["beta2_andrews_m2"] = float(beam_wander_variance_slant(
                w0, LAM, HS, cn2_zen, L_SLANT, f0=np.inf,
                elevation_deg=ELEV_DEG))
        except Exception as exc:
            out["beta2_andrews_m2"] = f"error: {exc}"
    return out


def main(variants_only=False):
    t_start = time.time()
    site = Site()
    cn2_base = default_cn2_profile(site, HS)
    results = {"config": {
        "lambda_m": LAM, "altitude_m": ALT_M, "elevation_deg": ELEV_DEG,
        "slant_range_m": L_SLANT, "aperture_m": APERTURE_M, "seed": SEED,
        "n_trials": N_TRIALS, "n_samples_f1": N_SAMPLES_F1,
        "hs_layers": int(HS.size),
    }, "cases": {}}

    # --variants reruns the VARIANT_CASE and its variants only, with lean
    # worker caps. The other cases keep their numbers in the run log.
    cases = ([c for c in CASES if c[0] == VARIANT_CASE] if variants_only
             else CASES)
    for name, w0, scale in cases:
        print(f"\n=== case {name}: w0 = {w0} m, Cn2 scale = {scale} ===")
        cn2 = cn2_base * scale
        scn = make_scenario(w0)
        f1 = run_fid1(w0, cn2)
        f2 = run_fid2(scn, cn2, w0, N_TRIALS, label="standard")
        an = run_andrews(w0, cn2)
        results["cases"][name] = {"w0_m": w0, "cn2_scale": scale,
                                  "fid1": f1, "fid2": [f2], "andrews": an}

        print(f"  fid1: sigma2_I = {f1['sigma2_I']:.3f}  "
              f"(wander-only {f1['sigma2_wander_only']:.3f}, "
              f"on-axis {f1['sigma2_on_axis']:.3f}, "
              f"sigma2_x = {f1['sigma2_x_mean']:.3f}, "
              f"weak_valid = {f1['weak_valid']})")
        print(f"  fid2: sigma2_I on-axis = {f2['sigma2_I_onaxis']:.3f}, "
              f"tracked = {f2['sigma2_I_tracked']:.3f}  "
              f"({f2['grid_n']} px, {f2['n_screens']} screens, "
              f"{f2['wall_s']:.0f} s)")
        print(f"  wander <beta^2>: sim {f2['wander_var_m2']:.3f}  "
              f"Dios {f1['beta2_dios_m2']:.3f}  "
              f"Andrews {an['beta2_andrews_m2']}")
        print(f"  widths: w_free {f1['w_free_m']:.2f}  "
              f"vac meas {f2['w_vac_meas_m']:.2f}  "
              f"w_st Dios {f1['w_st_m']:.2f} / sim {f2['w_st_meas_m']:.2f}  "
              f"w_lt Dios {f1['w_lt_m']:.2f} / sim {f2['w_lt_meas_m']:.2f}")
        print(f"  andrews: untracked {an['sigma2_I_untracked']}, "
              f"tracked {an['sigma2_I_tracked']}")

    # The convergence variants, on the key case only.
    key = next(c for c in CASES if c[0] == VARIANT_CASE)
    _, w0, scale = key
    cn2 = cn2_base * scale
    scn = make_scenario(w0)
    print(f"\n=== variants on {VARIANT_CASE} ===")
    # The 2048-px variants get 2 workers: each worker holds a full screen
    # stack, and the stacks are the memory cost of the run.
    for kwargs in ({"preset": "reference", "label": "reference", "workers": 2},
                   {"widen": 2.0, "label": "wide_x2", "workers": 2},
                   {"subharmonics": False, "label": "no_subharmonics",
                    "workers": 6}):
        v = run_fid2(scn, cn2, w0, N_TRIALS_VARIANT, **kwargs)
        results["cases"][VARIANT_CASE]["fid2"].append(v)
        print(f"  {v['label']:16s}: sigma2_I = {v['sigma2_I_onaxis']:.3f}, "
              f"tracked = {v['sigma2_I_tracked']:.3f}, "
              f"<beta^2> = {v['wander_var_m2']:.3f} m^2, "
              f"w_lt = {v['w_lt_meas_m']:.2f} m  "
              f"({v['grid_n']} px, {v['n_screens']} screens, "
              f"{v['wall_s']:.0f} s)")

    tag = "_variants" if variants_only else ""
    out = os.path.join(os.path.dirname(__file__),
                       f"uplink_farfield_reciprocity_results{tag}.json")
    with open(out, "w") as f:
        json.dump(results, f, indent=1)
    print(f"\nwrote {out}")
    print(f"(elapsed {time.time() - t_start:.0f} s)")


if __name__ == '__main__':
    import sys
    main(variants_only="--variants" in sys.argv)
