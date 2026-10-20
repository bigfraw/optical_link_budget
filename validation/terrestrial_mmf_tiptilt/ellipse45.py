"""A 1 inch terminal seen at 45 deg on a 10 km and a 1 km horizontal path.

THE CASE. 1550 nm, constant Cn2 of 1e-15, 5e-15 and 1e-14 m^-2/3, L0 = 25 m.
The launch is a collimated Gaussian of 1/e^2 diameter 4 x 4.65 mm (the
collimator beam through the 4x telescope), clipped by the 1 inch aperture.
The receiver is the same 1 inch aperture, focused at f = 4 x 17.5 mm = 70 mm
onto a 105 um core, NA 0.22, multimode fibre. The NA gate rho <= f NA =
15.4 mm is outside the aperture, so it is open (Snyder and Love,
DOI 10.1007/978-1-4613-2813-1). The core sits at the nominal focal plane.

THE 45 DEG APERTURE. A circular aperture seen at 45 deg projects to an
ellipse of semi-axes a (y) and a cos 45 (x, the tilt plane). Two campaigns
for each (L, Cn2): the launch clipped by the circle, and the launch clipped
by the ellipse (`Transmitter(tilt_deg=45)`). The launch field differs, so the
two need their own trials. The RECEIVE ellipse is a post-hoc read of the
stored patch (the patch is 1.5 aperture radii wide), so each campaign gives
both receive shapes; the case of record pairs the same shape at both ends.

THE READ (`EllipseRead`, through `Campaign.map_trials`). Each stored crop is
Fourier-upsampled (`mmf._resolve_auto_factor`, the same factor as the package
MMF read), and both receive masks are built on the FINE grid, so the
whole-pixel area error of a coarse mask (10 percent at 5 px across the
aperture) does not enter. For each mask:

    bucket  = collected_power * P_fine(mask) / P_coarse(in-run circle)
    eta     = mmf_coupling_efficiency(fine field, mask=mask, upsample=1)
    P_fibre / P_telescope = launch_fraction * bucket * eta

`collected_power` is referenced to the power AFTER the launch clip (the
terrestrial runner divides by the clipped start field), so the launch clip
`launch_fraction` comes back in from a fine numerical integral of the
Gaussian over the launch mask. No tip-tilt read: at f = 70 mm the tilt moves
the spot a few um against a 52.5 um core radius.

THE GRID. The sizer gives 2.0 to 3.2 mm pixels at 10 km (6 to 9 px across the
17.96 mm minor axis) and 1.5 mm at 1 km. The launch ellipse must resolve, so
the side stays sized and the pixel count is pinned: 4096 at 10 km and 1024 at
1 km (11 to 18 px and about 45 px across the minor axis). `--vacuum-check`
propagates the vacuum launch (`propagate_scenario`) on the pinned grid and on
2x the pixels at the same side and prints the difference.

Run from the repository root (on bigfraw, the GPU venv):

    python -m validation.terrestrial_mmf_tiptilt.ellipse45 --dry-run
    python -m validation.terrestrial_mmf_tiptilt.ellipse45 --vacuum-check
    python -m validation.terrestrial_mmf_tiptilt.ellipse45 --fft-backend cupy --workers 8
"""

import argparse
import dataclasses
import json
import os

import numpy as np

from olb.geometry import HorizontalPath
from olb.scenario import TerrestrialChannel, TerrestrialScenario
from olb.terminal import MMF, Terminal, Transmitter
from olb.waveoptics.mmf import (_resolve_auto_factor, _upsampled_field,
                                mmf_coupling_efficiency)
from olb.waveoptics.priority import boost_process_priority
from olb.waveoptics.run import propagate_scenario
from olb.waveoptics.turbulence import Campaign
from validation.terrestrial_mmf_tiptilt.run import (CN2, CORE_RADIUS_M, HERE,
                                                    LAM, NA, SEED, db,
                                                    refined_grid, stats)

APERTURE_M = 0.0254
WAIST_M = 4 * 4.65e-3 / 2          # 1/e^2 radius after the 4x telescope
F_MMF = 4 * 17.5e-3                # the collimator focal length through 4x
TILT_DEG = 45.0
N_PX = {10e3: 4096, 1e3: 1024}     # pinned pixel count on the sized side
TX = {"circle": None, "ellipse": TILT_DEG}


def masks(X, Y, aperture_m=APERTURE_M, tilt_deg=TILT_DEG):
    """The circle and the projected ellipse on a grid (x is the tilt plane)."""
    a = aperture_m / 2
    b = a * np.cos(np.radians(tilt_deg))
    return {"circle": X ** 2 + Y ** 2 <= a ** 2,
            "ellipse": (X / b) ** 2 + (Y / a) ** 2 <= 1}


def launch_fraction(tilt_deg, n=4096):
    """The power of the Gaussian launch inside the launch mask.

    A fine numerical integral of exp(-2 r^2 / w^2) over the mask (pixel
    w / 100). It is the clip loss that the runner's collected_power leaves
    out (it divides by the clipped start power).
    """
    x = (np.arange(n) - n // 2) * (2 * APERTURE_M / n)
    X, Y = np.meshgrid(x, x)
    g = np.exp(-2 * (X ** 2 + Y ** 2) / WAIST_M ** 2)
    m = masks(X, Y)["circle" if tilt_deg is None else "ellipse"]
    return float(g[m].sum() / g.sum())


def build_scenario(cn2, path_m, tilt_deg):
    near = Terminal(aperture_m=APERTURE_M, wavelength_m=LAM,
                    transmitter=Transmitter(waist_m=WAIST_M, tilt_deg=tilt_deg))
    mmf = MMF(core_radius_m=CORE_RADIUS_M, numerical_aperture=NA,
              focal_length_m=F_MMF)
    far = Terminal(aperture_m=APERTURE_M, wavelength_m=LAM, detector=mmf)
    channel = TerrestrialChannel(path_length_m=path_m,
                                 attenuation_db_per_km=0.0, cn2=cn2)
    scn = TerrestrialScenario(near=near, far=far, channel=channel)
    return scn, HorizontalPath(path_m)


def make_campaign(cn2, path_m, tx, preset, block_size, fft_backend):
    scn, geo = build_scenario(cn2, path_m, TX[tx])
    n_px = N_PX[path_m]
    grid, plan = refined_grid(scn, geo, preset, n_px)
    root = os.path.join(HERE, "campaigns",
                        f"ell45_L{path_m / 1e3:g}km_cn2{cn2:.0e}_{preset}"
                        f"_n{n_px}_tx{tx}")
    return Campaign(scn, geo, root, seed=SEED, preset=preset,
                    block_size=block_size, grid=grid, plan=plan,
                    fft_backend=fft_backend)


class EllipseRead:
    """Bucket and MMF eta of one stored trial for the circle and the ellipse."""

    def __call__(self, rec):
        F = rec.field()
        F.field = F.field.astype(np.complex128)
        Y, X = F.mgrid_cartesian
        p_coarse = np.sum(np.abs(F.field[masks(X, Y)["circle"]]) ** 2)
        M = _resolve_auto_factor("auto", F, APERTURE_M, CORE_RADIUS_M, F_MMF,
                                 NA, None)
        Fu = _upsampled_field(F, M) if M > 1 else F
        Yu, Xu = Fu.mgrid_cartesian
        out = []
        for m in masks(Xu, Yu).values():
            p_fine = np.sum(np.abs(Fu.field[m]) ** 2) / M ** 2
            eta = mmf_coupling_efficiency(Fu, APERTURE_M, CORE_RADIUS_M, F_MMF,
                                          numerical_aperture=NA, mask=m,
                                          upsample=1)
            out += [p_fine / p_coarse, eta]
        return np.array(out)


def summarise(camp, tx, n, workers):
    rec = camp.load(n, fields=False)
    p = np.array([t.collected_power for t in rec.trials])
    r_c, eta_c, r_e, eta_e = camp.map_trials(
        EllipseRead(), n_trials=n, workers=workers, screen_phase=False).T
    lf = launch_fraction(TX[tx])
    np.savez(os.path.join(camp.root_dir, "per_trial_ell45.npz"),
             collected=p, ratio_circle=r_c, eta_circle=eta_c,
             ratio_ellipse=r_e, eta_ellipse=eta_e, launch_fraction=lf)
    out = {"launch_clip_db": db(lf)}
    for rx, r, eta in (("rx_circle", r_c, eta_c), ("rx_ellipse", r_e, eta_e)):
        s = stats(lf * p * r, {"fibre": eta})
        out[rx] = {**s["fibre"], "bucket_mean_db": s["collected_mean_db"]}
    return out


def vacuum_check(preset):
    """The vacuum telescope-to-aperture loss on the pinned grid and on 2n."""
    for path_m, n_px in N_PX.items():
        for cn2, tx in ((c, t) for c in CN2 for t in TX):
            scn, geo = build_scenario(cn2, path_m, TX[tx])
            grid = refined_grid(scn, geo, preset, n_px)[0]
            loss = []
            for g in (grid, dataclasses.replace(grid, n=2 * n_px)):
                w = propagate_scenario(scn, geo, grid=g)
                loss.append(w.tx_truncation_db + w.geometric_loss_db)
            print(f"L {path_m / 1e3:g} km, grid of Cn2 {cn2:.0e} "
                  f"({grid.pixel_m * 1e3:.2f} mm), tx {tx}, rx circle: "
                  f"n {n_px}: {loss[0]:.3f} dB, n {2 * n_px}: {loss[1]:.3f} dB, "
                  f"diff {loss[0] - loss[1]:+.3f} dB", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--trials", type=int, default=1000)
    ap.add_argument("--block-size", type=int, default=50)
    ap.add_argument("--preset", default="standard")
    ap.add_argument("--workers", default="auto")
    ap.add_argument("--fft-backend", default="numpy")
    ap.add_argument("--cn2", type=float, nargs="*", default=CN2)
    ap.add_argument("--length-km", type=float, nargs="*", default=[10, 1])
    ap.add_argument("--tx", nargs="*", default=list(TX), choices=list(TX))
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--vacuum-check", action="store_true")
    args = ap.parse_args()
    workers = (None if args.workers == "none" else args.workers
               if args.workers == "auto" else int(args.workers))

    boost_process_priority()
    if args.vacuum_check:
        vacuum_check(args.preset)
        return
    out_path = os.path.join(HERE, f"results_ell45_{args.preset}.json")
    results = {}
    if os.path.exists(out_path):        # a subset run keeps the other cells
        with open(out_path) as f:
            results = json.load(f)
    for path_m in (L * 1e3 for L in args.length_km):
        for cn2 in args.cn2:
            for tx in args.tx:
                # The dry run sizes only: a Campaign writes its manifest (and
                # its backend) on creation, which a later run must then match.
                g, plan = refined_grid(*build_scenario(cn2, path_m, TX[tx]),
                                       args.preset, N_PX[path_m])
                print(f"L {path_m / 1e3:g} km, Cn2 {cn2:.0e}, tx {tx}: "
                      f"{g.n} px x {g.size_m:.3f} m ({g.pixel_m * 1e3:.2f} mm, "
                      f"{2 * APERTURE_M / 2 * np.cos(np.radians(TILT_DEG)) / g.pixel_m:.1f}"
                      f" px across the minor axis), {plan.z_m.size} screens",
                      flush=True)
                if args.dry_run:
                    continue
                camp = make_campaign(cn2, path_m, tx, args.preset,
                                     args.block_size, args.fft_backend)
                camp.run(args.trials, workers=workers, progress=True)
                key = f"L{path_m / 1e3:g}km_cn2{cn2:.0e}_tx{tx}"
                results[key] = summarise(camp, tx, args.trials, workers)
                print(json.dumps(results[key], indent=2), flush=True)
                with open(out_path, "w") as f:
                    json.dump(results, f, indent=2)


if __name__ == "__main__":
    main()
