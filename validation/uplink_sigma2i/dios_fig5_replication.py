'''
Replicate Dios et al. 2004, Fig. 5 (DOI 10.1364/AO.43.003866): the uplink
log-amplitude variance sigma2_chi against the transmit waist W0, for a
GEOSTATIONARY satellite at 0.84 um, elevations 90 and 30 deg.

THE POINT. Fig. 5 is the wave-optics validation that Dios gives for the
coupled-flux method: the semianalytic curve tracks an FFT-BPM simulation
until saturation appears (large W0). This script puts BOTH olb legs on that
exact case:

  - fidelity 1: the vendored coupled-flux kernels
    (olb.turbulence.uplink_flux._flux_result). If this curve overlays the
    paper's solid line, the vendored implementation is faithful.
  - fidelity 2: the split-step reciprocity far field
    (validation.uplink_sigma2i.uplink_farfield_reciprocity.run_fid2). If this leg matches
    fidelity 1 at small W0 (the regime the paper validates) and departs
    BELOW it at large W0 (where the paper's own FFT-BPM departs), then the
    fidelity-2 method is consistent with the paper's reference simulation,
    and the filled-beam disagreement of the LEO study is the expected
    saturation shortfall of the unsaturated off-axis term.

The estimator is the same on the two legs: sigma2_chi = var(ln I) / 4 over
the ensemble (I = exp(2 chi)). A constant flux rescale does not change it.

The script sets the module constants of uplink_farfield_reciprocity per
point (wavelength, geometry, launch aperture, far-field window). That keeps
one validated far-field code path instead of a second copy.

Run from the repo root (the fidelity-2 sweep is the slow part):
    python -m validation.uplink_sigma2i.dios_fig5_replication            # both legs
    python -m validation.uplink_sigma2i.dios_fig5_replication --fid1     # analytic leg only
'''

import json
import os
import sys
import time
import warnings

import numpy as np

import validation.uplink_sigma2i.uplink_farfield_reciprocity as R
from olb.beam import gaussz
from olb.geometry import CircularOrbit
from olb.turbulence.coupled_flux import (beam_wander_variance,
                                         on_axis_scintillation_index)
from olb.turbulence.profiles import DEFAULT_HS, default_cn2_profile
from olb.turbulence.uplink_flux import _flux_result
from olb.scenario import Site

LAM = 0.84e-6
ALT_M = 35786e3                    # geostationary
HS = DEFAULT_HS
ELEVATIONS = [90.0, 30.0]
W0_FID1 = np.geomspace(1e-3, 1e-1, 13)
W0_FID2 = [0.01, 0.03, 0.1]
N_SAMPLES_F1 = 20000
N_TRIALS_F2 = 120
WORKERS_F2 = 2                     # memory cap: one screen stack per worker
SEED = 2026


def set_case(elev_deg, w0):
    '''Point the far-field module at the GEO case for this (elev, w0).'''
    R.LAM = LAM
    R.ALT_M = ALT_M
    R.SEED = SEED
    R.GEOM = CircularOrbit(ALT_M, elevation_deg=[elev_deg])
    R.L_SLANT = float(np.asarray(R.GEOM.slant_range_m, dtype=float)[0])
    R.AIRMASS = 1.0 / np.sin(np.radians(elev_deg))
    # The launch aperture: wide against w0, small against the grid.
    R.APERTURE_M = max(0.5, 10.0 * w0)
    # The far-field window follows the GEO spot and the wander.
    w_free = float(gaussz(w0, R.L_SLANT, LAM))
    cn2_slant = default_cn2_profile(Site(), HS) * R.AIRMASS
    beta2 = float(beam_wander_variance(R.L_SLANT, cn2_slant,
                                       gaussz(w0, HS, LAM), HS))
    R.HALF_WINDOW_M = 3.0 * w_free + 6.0 * np.sqrt(beta2)
    return w_free, beta2


def sigma2_chi(samples):
    '''sigma2_chi = var(ln I) / 4. I = exp(2 chi).'''
    return float(np.var(np.log(np.asarray(samples))) / 4.0)


def main(fid1_only=False):
    t_start = time.time()
    cn2_zen = default_cn2_profile(Site(), HS)
    out = {"config": {"lambda_m": LAM, "altitude_m": ALT_M,
                      "n_samples_f1": N_SAMPLES_F1,
                      "n_trials_f2": N_TRIALS_F2, "seed": SEED},
           "fid1": {}, "fid2": {}}

    for elev in ELEVATIONS:
        rows = []
        for w0 in W0_FID1:
            set_case(elev, w0)
            np.random.seed(SEED)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                # As shipped: the wrapper puts the airmass on Cn2 and keeps
                # the VERTICAL height grid as the path coordinate. That is
                # short of the true slant integral by ~sec^(5/6) on the
                # on-axis term (the A, B path weights read z = h, not
                # z = h sec).
                r = _flux_result(w0, elev, R.L_SLANT, LAM, HS, cn2_zen,
                                 1.7e-14, N_SAMPLES_F1, 1)
                # Slant-corrected: elevation 90 (airmass 1) with the
                # SLANT-MAPPED height grid hs*sec and the zenith profile.
                # The trapz over the slant coordinate then carries the
                # ds = sec dh factor AND puts the path weights at the true
                # slant position. This is the exact path integral.
                np.random.seed(SEED)
                r_sl = _flux_result(w0, 90.0, R.L_SLANT, LAM,
                                    HS * R.AIRMASS, cn2_zen,
                                    1.7e-14, N_SAMPLES_F1, 1)
            k = 2 * np.pi / LAM
            s2_on = float(on_axis_scintillation_index(
                R.L_SLANT, k, gaussz(w0, R.L_SLANT, LAM),
                np.pi * w0 ** 2 / LAM, cn2_zen, HS * R.AIRMASS))
            rows.append({"w0_m": float(w0),
                         "sigma2_chi": sigma2_chi(r["Is_summed"]),
                         "sigma2_chi_slant": sigma2_chi(r_sl["Is_summed"]),
                         "sigma2_chi_on_axis": 0.25 * np.log1p(s2_on),
                         "sigma2_x_mean": r["sigma2_x_mean"],
                         "weak_valid": r["weak_fluctuation_valid"]})
        out["fid1"][str(elev)] = rows
        print(f"fid1, {elev:.0f} deg:")
        for row in rows:
            print(f"  w0 {row['w0_m']:.4f}  "
                  f"sigma2_chi {row['sigma2_chi']:.4f}  "
                  f"slant-corrected {row['sigma2_chi_slant']:.4f}  "
                  f"(on-axis {row['sigma2_chi_on_axis']:.4f}, "
                  f"valid {row['weak_valid']})")

    if not fid1_only:
        for elev in ELEVATIONS:
            rows = []
            for w0 in W0_FID2:
                w_free, beta2 = set_case(elev, w0)
                scn = R.make_scenario(w0)
                f2 = R.run_fid2(scn, cn2_zen, w0, N_TRIALS_F2,
                                label=f"geo_{elev:.0f}_{w0}",
                                workers=WORKERS_F2)
                # The same estimator as the fid1 leg, on the raw samples.
                s2chi = sigma2_chi(f2["eta_samples"])
                rows.append({"w0_m": float(w0), "sigma2_chi": s2chi,
                             "fid2": f2, "w_free_m": w_free,
                             "beta2_dios_m2": beta2})
                print(f"fid2, {elev:.0f} deg, w0 {w0}: sigma2_I "
                      f"{f2['sigma2_I_onaxis']:.3f} -> sigma2_chi "
                      f"{s2chi:.4f}  (grid {f2['grid_n']}, "
                      f"{f2['n_screens']} screens, {f2['wall_s']:.0f} s, "
                      f"wander sim {f2['wander_var_m2']:.1f} / Dios "
                      f"{beta2:.1f} m^2)")
            out["fid2"][str(elev)] = rows

    path = os.path.join(os.path.dirname(__file__),
                        "dios_fig5_replication_results.json")
    with open(path, "w") as f:
        json.dump(out, f, indent=1)
    print(f"\nwrote {path}")
    print(f"(elapsed {time.time() - t_start:.0f} s)")


if __name__ == '__main__':
    main(fid1_only="--fid1" in sys.argv)
