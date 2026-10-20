"""Fibre-coupled power of a 10 km terrestrial MMF link, with and without tip-tilt.

THE CASE. A horizontal 10 km path at 1550 nm. A collimated Gaussian launch of
1/e^2 DIAMETER 4 x 2.27 mm = 9.08 mm (waist radius 4.54 mm) from a 1 inch
(25.4 mm) aperture, into a 1 inch receive aperture that focuses onto a
105 um core (52.5 um radius) multimode fibre, NA 0.22. Three constant Cn2
values: 1e-15, 5e-15 and 1e-14 m^-2/3. The outer scale is the site value
(25 m). `--aperture-in 2` is the 2 inch family: a 50.8 mm aperture at both
ends and a launch of 1/e^2 diameter 16 x 2.27 mm = 36.32 mm, in campaigns
suffixed `_ap2in`; the rest of the case is the same.

THE RECORD. One fidelity-2 `Campaign` for each Cn2. Each trial stores the
collected power (a fraction of the launched power, so it holds the geometric
spread), the UNCORRECTED MMF coupling efficiency `mmf_eta`, and the receive
field patch. EVERY case reads the coupling POST HOC from the patch
(`CropEta` through `Campaign.map_trials`), so every case takes the same auto
upsample of the coarse pupil (docs/physics.md Section 9d). The TIP-TILT-corrected
efficiency is a POST-HOC read of the SAME trials: the slope fit (the
terrestrial default) and the far-field-centroid G-tilt. So the cases see the
SAME atmospheres. With the trials on disk, a run computes no new trial. A perfect
tip-tilt does not change the pupil power, so the fibre-coupled power is

    P_fibre / P_tx = collected_power * eta_MMF.

THE COUPLER. f = 4 x 12.56 mm = 50.24 mm (F_DESIGN), a FIXED design value; the
script does not optimise the focus. The 2 inch family reads f = 16 x 12.56 mm
= 200.96 mm (over D / (2 NA) = 115 mm, so its NA gate is open). f is under D / (2 NA) = 57.7 mm, so the NA
gate clips the outer pupil. The campaigns RAN with f = pi (D/2) a_core /
(1.12 lambda), about 1.21 m (F_MMF, the removed `optimal_focus` rule). That f
sets only the in-run scalar `mmf_eta` and the campaign key, so the stored
fields re-couple at F_DESIGN with no new propagation. The core sits at the NOMINAL
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
    python -m validation.terrestrial_mmf_tiptilt.run --n 2048 --workers 8  # post hoc
"""

import argparse
import dataclasses
import json
import os

import numpy as np

from olb.geometry import HorizontalPath
from olb.scenario import TerrestrialChannel, TerrestrialScenario
from olb.terminal import MMF, Terminal, TipTilt, Transmitter
from olb.waveoptics.mmf import mmf_coupling_efficiency
from olb.waveoptics.priority import boost_process_priority
from olb.waveoptics.turbulence import Campaign
from olb.waveoptics.turbulence.run import _PostCorrector
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
# F_MMF is the focal length the campaigns RAN with; it sets only the in-run
# scalar mmf_eta, but it is in the campaign key, so it stays to reopen them.
F_MMF = np.pi * (APERTURE_M / 2) * CORE_RADIUS_M / (1.12 * LAM)
# F_DESIGN is the coupler of record: fixed, not optimised. Every reported eta
# reads the stored fields at this f.
F_DESIGN = 4 * 12.56e-3
F_TAG = f"f{F_DESIGN * 1e3:.0f}mm"    # the 1 inch read tag


# The aperture families, keyed by the aperture in inches: (aperture, launch
# 1/e^2 radius, in-run MMF f, campaign-folder suffix, READ f). The in-run f
# is only in the campaign key, so it stays at the value each family RAN with.
# The READ f is the coupler of record: 4 x 12.56 mm for the 4 x 2.27 mm
# launch, 16 x 12.56 mm for the 16 x 2.27 mm launch.
APERTURES = {1: (APERTURE_M, WAIST_M, F_MMF, "", F_DESIGN),
             2: (0.0508, 16 * 2.27e-3 / 2, F_DESIGN, "_ap2in", 16 * 12.56e-3)}


def f_tag(f):
    return f"f{f * 1e3:.0f}mm"


def build_scenario(cn2, aperture_in=1):
    """Give the (TerrestrialScenario, HorizontalPath) pair of one Cn2."""
    ap, waist, f = APERTURES[aperture_in][:3]
    near = Terminal(aperture_m=ap, wavelength_m=LAM,
                    transmitter=Transmitter(waist_m=waist))
    mmf = MMF(core_radius_m=CORE_RADIUS_M, numerical_aperture=NA,
              focal_length_m=f)
    far = Terminal(aperture_m=ap, wavelength_m=LAM, detector=mmf)
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


def make_campaign(cn2, preset, block_size, fft_backend="numpy", n_px=None,
                  aperture_in=1):
    scn, geo = build_scenario(cn2, aperture_in)
    root = os.path.join(HERE, "campaigns", f"cn2{cn2:.0e}_{preset}_focal"
                        + (f"_n{n_px}" if n_px else "")
                        + APERTURES[aperture_in][3])
    grid, plan = refined_grid(scn, geo, preset, n_px) if n_px else (None, None)
    return Campaign(scn, geo, root, seed=SEED, preset=preset,
                    block_size=block_size, grid=grid, plan=plan,
                    fft_backend=fft_backend)


def db(x):
    return -10 * np.log10(x)


class CropEta:
    """The MMF eta of one stored trial at focal length f: open, slope TT, G-tilt TT.

    It focuses the stored CROP, not the full grid. `Campaign.recouple` pads
    the crop back to the full grid, and at f = 50 mm the auto upsample asks
    for M = 16, so a 2048 px grid becomes 32768 px and swaps. The upsampled
    pupil pixel and so the focal field of view lam f M / dx are the SAME on
    the crop; only the focal pixel lam f / side grows (about 2 um against a
    52.5 um core radius). The crop holds the 1.5x patch margin, so "auto"
    still engages and the clip is built on the fine grid.
    """

    def __init__(self, aperture_m=APERTURE_M, f=F_DESIGN):
        self.aperture_m = aperture_m
        self.f = f

    def __call__(self, rec):
        c = rec.context
        if not c:
            for src in ("slopes", "gtilt"):
                c[src] = _PostCorrector(rec.patch, [TipTilt()], self.aperture_m,
                                        0.0, src)
        out = []
        for a in (rec.array, c["slopes"].correct(rec.array),
                  c["gtilt"].correct(rec.array)):
            F = rec.field()
            F.field = a
            out.append(mmf_coupling_efficiency(
                F, self.aperture_m, CORE_RADIUS_M, self.f,
                numerical_aperture=NA, upsample="auto"))
        return np.array(out)


def summarise(camp, n, workers, f=F_DESIGN):
    """Give the loss statistics (positive dB) of one stored campaign.

    Every eta is a POST-HOC read of the stored fields at the focal length f,
    so the stored campaigns (whatever f they ran with) need no new
    propagation.
    """
    rec = camp.load(n, fields=False)
    p = np.array([t.collected_power for t in rec.trials])
    eta_open, eta_tt, eta_g = camp.map_trials(
        CropEta(camp.scenario.rx_terminal.aperture_m, f), n_trials=n, workers=workers, screen_phase=False).T
    np.savez(os.path.join(camp.root_dir, f"per_trial_{f_tag(f)}.npz"),
             collected=p, eta_open=eta_open, eta_tt=eta_tt, eta_gtilt=eta_g)
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
    ap.add_argument("--aperture-in", type=int, default=1, choices=APERTURES,
                    help="1: 1 inch, 9.08 mm launch; 2: 2 inch, 36.32 mm")
    ap.add_argument("--dry-run", action="store_true",
                    help="size the grids and the screens, run no trial")
    args = ap.parse_args()
    workers = (args.workers if args.workers == "auto"
               else int(args.workers))

    boost_process_priority()
    f_read = APERTURES[args.aperture_in][4]
    tag = (f"{args.preset}_focal" + (f"_n{args.n}" if args.n else "")
           + APERTURES[args.aperture_in][3] + f"_{f_tag(f_read)}")
    out_path = os.path.join(HERE, f"results_{tag}.json")
    results = {}
    if os.path.exists(out_path):        # a --cn2 subset keeps the others
        with open(out_path) as f:
            results = json.load(f)
    for cn2 in args.cn2:
        camp = make_campaign(cn2, args.preset, args.block_size,
                             args.fft_backend, args.n, args.aperture_in)
        print(f"Cn2 {cn2:.0e}: grid {camp.grid.n} px x "
              f"{camp.grid.size_m:.3f} m, {camp.plan.z_m.size} screens, "
              f"coupled at f = {f_read * 1e3:.2f} mm",
              flush=True)
        if args.dry_run:
            continue
        camp.run(args.trials, workers=workers, progress=True)
        results[f"{cn2:.0e}"] = summarise(camp, args.trials, workers, f_read)
        print(json.dumps(results[f"{cn2:.0e}"], indent=2), flush=True)
        with open(out_path, "w") as f:
            json.dump(results, f, indent=2)


if __name__ == "__main__":
    main()
