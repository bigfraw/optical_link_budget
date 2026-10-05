"""Fibre-coupled power of a 10 km terrestrial MMF link, with and without tip-tilt.

THE CASE. A horizontal 10 km path at 1550 nm. A collimated Gaussian launch of
1/e^2 DIAMETER 4 x 2.27 mm = 9.08 mm (waist radius 4.54 mm) from a 1 inch
(25.4 mm) aperture, into a 1 inch receive aperture that focuses onto a
105 um core (52.5 um radius) multimode fibre, NA 0.22. Three constant Cn2
values: 5e-15, 1e-14 and 5e-14 m^-2/3. The outer scale is the site value
(25 m).

THE RECORD. One fidelity-2 `Campaign` for each Cn2. Each trial stores the
collected power (a fraction of the launched power, so it holds the geometric
spread), the UNCORRECTED MMF coupling efficiency `mmf_eta`, and the receive
field patch. The TIP-TILT-corrected efficiency is a POST-HOC read of the SAME
trials (`Campaign.recouple_compensated([TipTilt()], ...)`, slope sensing, the
terrestrial default). So the two cases see the SAME atmospheres. A perfect
tip-tilt does not change the pupil power, so the fibre-coupled power is

    P_fibre / P_tx = collected_power * eta_MMF.

THE COUPLER. `optimal_focus=True` sets f = pi (D/2) a_core / (1.12 lambda)
(about 1.2 m here, so the NA gate is open). `defocus_m` is set to the
received-curvature focus shift (S. A. Self, Appl. Opt. 22 (1983) 658,
DOI 10.1364/AO.22.000658), the documented recipe for an aligned coupler.

THE TIP-TILT. A perfect modal fit of the first 3 Noll modes (R. J. Noll,
DOI 10.1364/JOSA.66.000207), sensed from the wrapped-gradient slopes of the
receive field. It is the UPPER BOUND of a real tracker: no noise, no servo lag.
At Cn2 = 5e-14 the slope route can alias (it warns past 2.8 rad per pixel).

Run it from the repository root:

    python -m validation.terrestrial_mmf_tiptilt.run --dry-run
    python -m validation.terrestrial_mmf_tiptilt.run --trials 1000 --workers auto
"""

import argparse
import json
import os

import numpy as np

from olb.geometry import HorizontalPath
from olb.models.coupling.terrestrial import (_mmf_focal_length,
                                             curvature_focus_shift)
from olb.scenario import TerrestrialChannel, TerrestrialScenario
from olb.terminal import MMF, Terminal, TipTilt, Transmitter
from olb.waveoptics.turbulence import Campaign

HERE = os.path.dirname(os.path.abspath(__file__))

LAM = 1550e-9
PATH_M = 10e3
CN2 = (5e-15, 1e-14, 5e-14)
APERTURE_M = 0.0254
WAIST_M = 4 * 2.27e-3 / 2          # 1/e^2 radius of the 9.08 mm diameter
CORE_RADIUS_M = 105e-6 / 2
NA = 0.22                          # common 105 um step-index fibre
SEED = 20261005


def build_scenario(cn2):
    """Give the (TerrestrialScenario, HorizontalPath) pair of one Cn2."""
    near = Terminal(aperture_m=APERTURE_M, wavelength_m=LAM,
                    transmitter=Transmitter(waist_m=WAIST_M))
    mmf = MMF(core_radius_m=CORE_RADIUS_M, numerical_aperture=NA,
              optimal_focus=True)
    far = Terminal(aperture_m=APERTURE_M, wavelength_m=LAM, detector=mmf)
    channel = TerrestrialChannel(path_length_m=PATH_M,
                                 attenuation_db_per_km=0.0, cn2=cn2)
    scn = TerrestrialScenario(near=near, far=far, channel=channel)
    # Put the core at the true focus of the received diverging beam.
    mmf.defocus_m = curvature_focus_shift(scn)
    return scn, HorizontalPath(PATH_M)


def make_campaign(cn2, preset, block_size):
    scn, geo = build_scenario(cn2)
    root = os.path.join(HERE, "campaigns", f"cn2{cn2:.0e}_{preset}")
    return Campaign(scn, geo, root, seed=SEED, preset=preset,
                    block_size=block_size)


def db(x):
    return -10 * np.log10(x)


def summarise(camp, n, workers):
    """Give the loss statistics (positive dB) of one stored campaign."""
    rec = camp.load(n, fields=False)
    p = np.asarray(rec.collected_power)
    det = camp.scenario.rx_terminal.detector
    eta_open = np.asarray(rec.mmf_eta)
    eta_tt = camp.recouple_compensated([TipTilt()], det, n_trials=n,
                                       workers=workers)
    out = {}
    for name, eta in (("no_tt", eta_open), ("tip_tilt", eta_tt)):
        loss = db(p * eta)
        out[name] = {"mean_power_db": float(db(np.mean(p * eta))),
                     "median_db": float(np.median(loss)),
                     "p95_db": float(np.percentile(loss, 95)),
                     "p99_db": float(np.percentile(loss, 99)),
                     "mean_eta": float(np.mean(eta))}
    out["collected_mean_db"] = float(db(np.mean(p)))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--trials", type=int, default=1000)
    ap.add_argument("--block-size", type=int, default=50)
    ap.add_argument("--preset", default="standard")
    ap.add_argument("--workers", default="auto")
    ap.add_argument("--dry-run", action="store_true",
                    help="size the grids and the screens, run no trial")
    args = ap.parse_args()
    workers = (args.workers if args.workers == "auto"
               else int(args.workers))

    results = {}
    for cn2 in CN2:
        camp = make_campaign(cn2, args.preset, args.block_size)
        print(f"Cn2 {cn2:.0e}: grid {camp.grid.n} px x "
              f"{camp.grid.size_m:.3f} m, {camp.plan.z_m.size} screens, "
              f"f = {_mmf_focal_length(camp.scenario.rx_terminal.detector, APERTURE_M, LAM):.3f} m",
              flush=True)
        if args.dry_run:
            continue
        camp.run(args.trials, workers=workers, progress=True)
        results[f"{cn2:.0e}"] = summarise(camp, args.trials, workers)
        print(json.dumps(results[f"{cn2:.0e}"], indent=2), flush=True)
        with open(os.path.join(HERE, f"results_{args.preset}.json"), "w") as f:
            json.dump(results, f, indent=2)


if __name__ == "__main__":
    main()
