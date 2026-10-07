"""Fibre-coupled power of a 10 km terrestrial MMF link, with and without tip-tilt.

THE CASE. A horizontal 10 km path at 1550 nm. A collimated Gaussian launch of
1/e^2 DIAMETER 4 x 2.27 mm = 9.08 mm (waist radius 4.54 mm) from a 1 inch
(25.4 mm) aperture, into a 1 inch receive aperture that focuses onto a
105 um core (52.5 um radius) multimode fibre, NA 0.22. Three constant Cn2
values: 1e-15, 5e-15 and 1e-14 m^-2/3. The outer scale is the site value
(25 m).

THE RECORD. One fidelity-2 `Campaign` for each Cn2. Each trial stores the
collected power (a fraction of the launched power, so it holds the geometric
spread), the UNCORRECTED MMF coupling efficiency `mmf_eta`, and the receive
field patch. BOTH cases read the coupling POST HOC from the patch
(`recouple` and `recouple_compensated`), so both take the same auto
upsample of the coarse pupil (docs/physics.md Section 9d). The TIP-TILT-corrected efficiency is a POST-HOC read of the SAME
trials (`Campaign.recouple_compensated([TipTilt()], ...)`, slope sensing, the
terrestrial default). So the two cases see the SAME atmospheres. A perfect
tip-tilt does not change the pupil power, so the fibre-coupled power is

    P_fibre / P_tx = collected_power * eta_MMF.

THE COUPLER. f = pi (D/2) a_core / (1.12 lambda), about 1.21 m: the spot just
matches the core (the NA gate is open). That was the removed MMF
`optimal_focus` rule; the campaigns were run with it, so F_MMF keeps it. It is
NOT the capture optimum: a shorter f couples much more (focal_sweep.py). The core sits at the NOMINAL
focal plane (`defocus_m = 0`), NOT at the received-curvature focus shift
(S. A. Self, Appl. Opt. 22 (1983) 658, DOI 10.1364/AO.22.000658), so the
coupling pays the received-curvature defocus.

THE TIP-TILT. A perfect modal fit of the first 3 Noll modes (R. J. Noll,
DOI 10.1364/JOSA.66.000207), sensed from the wrapped-gradient slopes of the
receive field. It is the UPPER BOUND of a real tracker: no noise, no servo lag.
At strong Cn2 the slope route can alias (it warns past 2.8 rad per pixel).

Run it from the repository root:

    python -m validation.terrestrial_mmf_tiptilt.run --dry-run
    python -m validation.terrestrial_mmf_tiptilt.run --trials 1000 --workers auto
    python -m validation.terrestrial_mmf_tiptilt.run --n 4096 --fft-backend cupy
"""

import argparse
import dataclasses
import json
import os

import numpy as np

from olb.geometry import HorizontalPath
from olb.scenario import TerrestrialChannel, TerrestrialScenario
from olb.terminal import MMF, Terminal, TipTilt, Transmitter
from olb.waveoptics.priority import boost_process_priority
from olb.waveoptics.turbulence import Campaign
from olb.waveoptics.turbulence.sampling import turbulent_grid

HERE = os.path.dirname(os.path.abspath(__file__))

LAM = 1550e-9
PATH_M = 10e3
CN2 = (1e-15, 5e-15, 1e-14)
APERTURE_M = 0.0254
WAIST_M = 4 * 2.27e-3 / 2          # 1/e^2 radius of the 9.08 mm diameter
CORE_RADIUS_M = 105e-6 / 2
NA = 0.22                          # common 105 um step-index fibre
SEED = 20261005
F_MMF = np.pi * (APERTURE_M / 2) * CORE_RADIUS_M / (1.12 * LAM)


def build_scenario(cn2):
    """Give the (TerrestrialScenario, HorizontalPath) pair of one Cn2."""
    near = Terminal(aperture_m=APERTURE_M, wavelength_m=LAM,
                    transmitter=Transmitter(waist_m=WAIST_M))
    mmf = MMF(core_radius_m=CORE_RADIUS_M, numerical_aperture=NA,
              focal_length_m=F_MMF)
    far = Terminal(aperture_m=APERTURE_M, wavelength_m=LAM, detector=mmf)
    channel = TerrestrialChannel(path_length_m=PATH_M,
                                 attenuation_db_per_km=0.0, cn2=cn2)
    scn = TerrestrialScenario(near=near, far=far, channel=channel)
    return scn, HorizontalPath(PATH_M)


def refined_grid(scn, geo, preset, n_px):
    """The sized grid and plan, with n_px pixels on the SAME side.

    The sizer asks for 8192 to 16384 px here and n_max caps it at 2048 (a
    4 to 5 mm pixel, 5 px across the 1 inch aperture). The smoke runs
    (analyse_smoke.py) show 4096 px fixes the point index, the phase
    aliasing and the slope tip-tilt fit; the plan does not depend on n.
    """
    grid, plan = turbulent_grid(scn, geo, preset=preset,
                                L0_m=scn.channel.site.outer_scale_m)[:2]
    return dataclasses.replace(grid, n=n_px), plan


def make_campaign(cn2, preset, block_size, fft_backend="numpy", n_px=None):
    scn, geo = build_scenario(cn2)
    root = os.path.join(HERE, "campaigns", f"cn2{cn2:.0e}_{preset}_focal"
                        + (f"_n{n_px}" if n_px else ""))
    grid, plan = refined_grid(scn, geo, preset, n_px) if n_px else (None, None)
    return Campaign(scn, geo, root, seed=SEED, preset=preset,
                    block_size=block_size, grid=grid, plan=plan,
                    fft_backend=fft_backend)


def db(x):
    return -10 * np.log10(x)


def summarise(camp, n, workers):
    """Give the loss statistics (positive dB) of one stored campaign."""
    rec = camp.load(n, fields=False)
    p = np.array([t.collected_power for t in rec.trials])
    det = camp.scenario.rx_terminal.detector
    # The post-hoc read for BOTH cases: the in-run scalar `mmf_eta` couples a
    # pre-clipped field at M=1, so it does not upsample the coarse pupil.
    eta_open = camp.recouple(det, n_trials=n, workers=workers)
    eta_tt = camp.recouple_compensated([TipTilt()], det, n_trials=n,
                                       workers=workers)
    # Cross-check: remove the G-tilt (mean wrapped gradient) in place of the
    # slope fit, which failed on a 5 px pupil (analyse_smoke.py).
    from validation.terrestrial_mmf_tiptilt.analyse_smoke import own_tt_eta
    eta_g = np.array([own_tt_eta(camp.field(i).field.astype(complex),
                                 camp.patch.pixel_m)[1] for i in range(n)])
    np.savez(os.path.join(camp.root_dir, "per_trial.npz"), collected=p,
             eta_open=eta_open, eta_tt=eta_tt, eta_gtilt=eta_g)
    return stats(p, {"no_tt": eta_open, "tip_tilt": eta_tt,
                     "tip_tilt_gtilt": eta_g})


def stats(p, etas, B=2000, seed=1):
    """Each statistic with its bootstrap standard error over the trials.

    Every value is [value, SE]. The resample draws the trials with
    replacement, so it keeps the pairing of p and eta in one trial.
    """
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, p.size, (B, p.size))

    def both(fn):
        return [float(fn(np.arange(p.size))),
                float(np.std([fn(i) for i in idx]))]

    out = {}
    for name, eta in etas.items():
        loss = db(p * eta)
        out[name] = {
            "mean_power_db": both(lambda i: db(np.mean(p[i] * eta[i]))),
            "median_db": both(lambda i: np.median(loss[i])),
            "p95_db": both(lambda i: np.percentile(loss[i], 95)),
            "p99_db": both(lambda i: np.percentile(loss[i], 99)),
            "mean_eta": both(lambda i: np.mean(eta[i]))}
    out["collected_mean_db"] = both(lambda i: db(np.mean(p[i])))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--trials", type=int, default=1000)
    ap.add_argument("--block-size", type=int, default=50)
    ap.add_argument("--preset", default="standard")
    ap.add_argument("--workers", default="auto")
    ap.add_argument("--fft-backend", default="numpy",
                    help='"cupy" runs the trials on the CUDA device')
    ap.add_argument("--cn2", type=float, nargs="*", default=CN2,
                    help="the Cn2 values to run (default: all three)")
    ap.add_argument("--n", type=int, default=None,
                    help="pixel count on the sized side (None: sized grid)")
    ap.add_argument("--dry-run", action="store_true",
                    help="size the grids and the screens, run no trial")
    args = ap.parse_args()
    workers = (args.workers if args.workers == "auto"
               else int(args.workers))

    boost_process_priority()
    tag = f"{args.preset}_focal" + (f"_n{args.n}" if args.n else "")
    out_path = os.path.join(HERE, f"results_{tag}.json")
    results = {}
    if os.path.exists(out_path):        # a --cn2 subset keeps the others
        with open(out_path) as f:
            results = json.load(f)
    for cn2 in args.cn2:
        camp = make_campaign(cn2, args.preset, args.block_size,
                             args.fft_backend, args.n)
        print(f"Cn2 {cn2:.0e}: grid {camp.grid.n} px x "
              f"{camp.grid.size_m:.3f} m, {camp.plan.z_m.size} screens, "
              f"f = {F_MMF:.3f} m",
              flush=True)
        if args.dry_run:
            continue
        camp.run(args.trials, workers=workers, progress=True)
        results[f"{cn2:.0e}"] = summarise(camp, args.trials, workers)
        print(json.dumps(results[f"{cn2:.0e}"], indent=2), flush=True)
        with open(out_path, "w") as f:
            json.dump(results, f, indent=2)


if __name__ == "__main__":
    main()
