"""Arm A: the direct reciprocity check of the fidelity-2 uplink overlap.

THE TRUTH. Launch psi_tx at the ground, propagate it UP through the screens
(the reverse order), and read the on-axis field at the satellite with the
Fresnel sum over the vacuum hop above the slab:

    E_sat = sum( E_top(r) exp(i k r^2 / (2 L_rest)) )

(Goodman, Introduction to Fourier Optics, Eq. (4-17), the Fresnel diffraction
integral on the axis; ISBN 978-0974707723.) The upward grid is WIDE, so the
diverged beam never reaches its absorbing band.

THE MODEL UNDER TEST. The production overlap eta = |sum(F psi_tx)|^2 / o_vac
(Shapiro, DOI 10.1364/JOSA.61.000492) on the PRODUCTION grid, where F is the
DOWNLINK slab field from the start field of Gate 0 (a wide Gaussian, or the
plane wave of record). The downward screens are the CENTRAL WINDOW of the wide
screens, so the two routes read the SAME atmosphere near the axis.

STEP 0 checks that the hand loop IS the production loop: on the production
grid with its own screens, the hand eta equals the runner eta bit for bit.

Run from the repository root:
    python -m validation.uplink_divergence.arm_a --elev 30 --trials 100
"""
import argparse
import dataclasses
import json
import warnings
from pathlib import Path

import numpy as np

from olb.waveoptics.field import Begin
from olb.waveoptics.sources import GaussBeam
from olb.waveoptics.turbulence.run import (_ground_transmit_mode,
                                           _resolve_seed, _screen_builder,
                                           _screen_seed,
                                           propagate_turbulent_scenario,
                                           reciprocity_overlap)
from olb.waveoptics.turbulence.sampling import PRESETS, resolve_outer_scale
from olb.waveoptics.turbulence.splitstep import (split_step,
                                                 super_gaussian_boundary)

from .config import THETA_MIN, geometry, scenario
from .gate0 import base_grid, label, start_field

warnings.simplefilter("ignore")
HERE = Path(__file__).parent


def start(kind, g, lam, ws_frac):
    F = Begin(g.size_m, lam, g.n)
    return F if kind == "plane" else GaussBeam(F, ws_frac * g.size_m)


def down_eta(F0, g, plan, screens, psi, mask, o_vac):
    F = split_step(F0, plan.z_m, screens, plan.z_total_m, boundary=mask)
    return reciprocity_overlap(F.field, psi) / o_vac


def step0(elev, preset, ws_frac, seed=7, n=3):
    """Assert that the hand downward loop equals the runner, bit for bit."""
    g, plan = base_grid(elev, preset)
    scn, lam = scenario(8 * THETA_MIN), scenario().ground.wavelength_m
    with start_field("gauss", ws_frac):
        w = propagate_turbulent_scenario(scn, geometry(elev), n_trials=n,
                                         seed=seed, preset=preset, grid=g,
                                         plan=plan, precision="double")
    mask = super_gaussian_boundary(g.n, PRESETS[preset].boundary_width_frac)
    psi = _ground_transmit_mode(scn.ground, g)
    flat = [np.zeros((g.n, g.n))] * plan.z_m.size
    o_vac = down_eta(start("gauss", g, lam, ws_frac), g, plan, flat, psi,
                     mask, 1.0)
    build = _screen_builder("olb", g, resolve_outer_scale(None, scn), True)
    ent = _resolve_seed(seed)
    for k in range(n):
        scr = [build(_screen_seed(ent, k, j), plan.r0_m[j])
               for j in range(plan.z_m.size)]
        e = down_eta(start("gauss", g, lam, ws_frac), g, plan, scr, psi,
                     mask, o_vac)
        assert e == w.trials[k].eta_turb, (k, e, w.trials[k].eta_turb)
    print(f"step 0: the hand loop equals the runner on {n} trials")


_S = {}


def setup(elev, preset, ws_fracs, n_big, mults, mult=1):
    """Build the fixed parts of one run ONE time in each process."""
    key = (elev, preset, tuple(ws_fracs), n_big, tuple(mults), mult)
    if key in _S:
        return _S[key]
    g, plan = base_grid(elev, preset)
    g = dataclasses.replace(g, n=g.n * mult)
    lam = scenario().ground.wavelength_m
    gb = dataclasses.replace(g, n=int(n_big), size_m=g.pixel_m * n_big)
    c0 = (n_big - g.n) // 2
    L_rest = float(np.asarray(geometry(elev).slant_range_m)) - plan.z_total_m
    bw = PRESETS[preset].boundary_width_frac
    Y, X = Begin(gb.size_m, lam, n_big).mgrid_cartesian
    divs = [None if m == 0 else m * THETA_MIN for m in mults]
    s = dict(g=g, gb=gb, plan=plan, lam=lam, L_rest=L_rest, divs=divs,
             win=np.s_[c0:c0 + g.n, c0:c0 + g.n],
             mask=super_gaussian_boundary(g.n, bw),
             mask_b=super_gaussian_boundary(n_big, bw),
             readout=np.exp(2j * np.pi / lam * (X ** 2 + Y ** 2) / (2 * L_rest)),
             z_up=plan.z_total_m - plan.z_m[::-1],
             psis=[_ground_transmit_mode(scenario(d).ground, g) for d in divs],
             psis_b=[_ground_transmit_mode(scenario(d).ground, gb) for d in divs],
             starts={"plane": start("plane", g, lam, None),
                     **{f"gauss{w}": start("gauss", g, lam, w) for w in ws_fracs}},
             build=_screen_builder("olb", gb, resolve_outer_scale(None, scenario()), True))
    flat, flat_b = ([np.zeros((m, m))] * plan.z_m.size for m in (g.n, n_big))
    top_vac = [_up(s, p, flat_b) for p in s["psis_b"]]
    s["edge"] = [float(np.abs(t[~(s["mask_b"] > 0.999)]).max() / np.abs(t).max())
                 for t in top_vac]
    s["sat_vac"] = [abs((t * s["readout"]).sum()) ** 2 for t in top_vac]
    s["o_vac"] = {k: [reciprocity_overlap(split_step(F0, plan.z_m, flat,
                                                     plan.z_total_m,
                                                     boundary=s["mask"]).field, p)
                      for p in s["psis"]] for k, F0 in s["starts"].items()}
    _S[key] = s
    return s


def _up(s, psi, scr):
    F = Begin(s["gb"].size_m, s["lam"], s["gb"].n)
    F.field = psi.astype(complex)
    return split_step(F, s["z_up"], scr[::-1], s["plan"].z_total_m,
                      boundary=s["mask_b"]).field


def trial(args):
    """One trial: the upward truth and the two downward overlaps."""
    *cfg, ent = args[0]
    k = args[1]
    s = setup(*cfg)
    plan = s["plan"]
    big = [s["build"](_screen_seed(ent, k, j), plan.r0_m[j])
           for j in range(plan.z_m.size)]
    small = [b[s["win"]] for b in big]
    out = {"up": [abs((_up(s, p, big) * s["readout"]).sum()) ** 2 / v
                  for p, v in zip(s["psis_b"], s["sat_vac"])]}
    for name, F0 in s["starts"].items():
        F = split_step(F0, plan.z_m, small, plan.z_total_m,
                       boundary=s["mask"]).field
        out[name] = [reciprocity_overlap(F, p) / o
                     for p, o in zip(s["psis"], s["o_vac"][name])]
    return out


def _init():
    from olb.waveoptics.priority import boost_process_priority
    warnings.simplefilter("ignore")
    boost_process_priority()


def run(elev, preset, ws_fracs, n_big, n_trials, mults, workers, mult=1, seed=11):
    from concurrent.futures import ProcessPoolExecutor
    cfg = (elev, preset, tuple(ws_fracs), n_big, tuple(mults), mult, _resolve_seed(seed))
    s = setup(*cfg[:6])
    g, gb, divs, edge = s["g"], s["gb"], s["divs"], s["edge"]
    with ProcessPoolExecutor(workers, initializer=_init) as ex:
        rows = list(ex.map(trial, [(cfg, k) for k in range(n_trials)]))
    res = {k: np.array([r[k] for r in rows]) for k in rows[0]}
    print(f"elev {elev} {preset}: grid {g.n} px / {g.size_m:.3f} m, wide "
          f"{gb.n} px / {gb.size_m:.2f} m, L_rest {s['L_rest'] / 1e3:.0f} km, "
          f"{n_trials} trials, ws {ws_fracs} x side")
    summary = {}
    for i, d in enumerate(divs):
        up_i = res["up"][:, i]
        row = {"up_mean": up_i.mean(), "edge_leak": edge[i]}
        msg = (f"  {label(d):5s} up mean {up_i.mean():.3f} (edge leak "
               f"{edge[i]:.1e})")
        for name in s["starts"]:
            e = res[name][:, i]
            ratio = np.log10(e / up_i) * 10
            row.update({f"{name}_mean": e.mean(),
                        f"{name}_corr": np.corrcoef(e, up_i)[0, 1],
                        f"{name}_db_median": np.median(ratio),
                        f"{name}_db_rms": ratio.std()})
            msg += (f" | {name} mean {e.mean():.3f} corr "
                    f"{row[f'{name}_corr']:.3f} dB/trial median "
                    f"{np.median(ratio):+.2f} rms {ratio.std():.2f}")
        summary[label(d)] = {k: float(v) for k, v in row.items()}
        print(msg)
    tag = f"arm_a_e{elev:.0f}_{preset}_n{g.n}_w{gb.n}"
    np.savez(HERE / f"{tag}.npz", **res, divs=np.array(
        [0 if d is None else d for d in divs]))
    (HERE / f"{tag}.json").write_text(json.dumps(summary, indent=1))


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument("--elev", type=float, default=30.0)
    ap.add_argument("--preset", default="rapid")
    ap.add_argument("--ws", type=float, nargs="+", default=[0.15, 0.2, 0.3])
    ap.add_argument("--mult", type=int, default=1, help="pixel refinement of the production grid")
    ap.add_argument("--n-big", type=int, default=1536)
    ap.add_argument("--trials", type=int, default=100)
    ap.add_argument("--mults", type=int, nargs="+", default=[0, 2, 4])
    ap.add_argument("--skip-step0", action="store_true")
    ap.add_argument("--workers", type=int, default=12)
    a = ap.parse_args()
    if not a.skip_step0:
        step0(a.elev, a.preset, a.ws[0])
    run(a.elev, a.preset, a.ws, a.n_big, a.trials, a.mults, a.workers, a.mult)
