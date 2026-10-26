"""The tilt at the two ends of the 45 deg link, trial by trial (step 2b).

THE QUESTION. Each terminal can sense only the tilt of the beam that it
RECEIVES. Does the near terminal's tilt (of the REVERSE beam) follow the
forward beam's wander at the far end, as reciprocity says (Shapiro,
DOI 10.1364/JOSA.61.000492)? And how alike are the tilts at the two ends at
the same instant? See NEXT_tx_tiptilt.md, step 2b.

THE ROUTE (no package change). For each stored trial, this script builds the
seeded screens of that trial ONE time, then:

1. FORWARD: the near launch through the screens in order, to the far plane
   (the run of record; the stored collected power is the identity gate).
2. REVERSE: the far launch through the SAME screens in reverse order, at the
   mirrored distances L - z_j, to the near plane. The two terminals are the
   same hardware, so the far launch is the near launch field. A thin screen
   gives the same phase in both directions, and the transverse axes stay the
   same, so the ellipse of each end has the same orientation.

At each receive plane it keeps the beam moments (`wander_split.moments`) and
the G-tilt and Z-tilt over the circle and the ellipse aperture with the MMF
eta (`ellipse45.TiltRead`, on a crop of 1.5 aperture radii, as the stored
patch). So theta_far, c_far, theta_near and c_near pair trial for trial.

Run from the repository root (on bigfraw, the GPU venv):

    <python> -m validation.terrestrial_mmf_tiptilt.two_end_tilt --cn2 1e-15 --length-km 1 --tx ellipse
"""

import argparse
import json
import os
import time

import numpy as np

from olb.waveoptics.field import Begin, Field, Power, field_dtype, to_host
from olb.waveoptics.priority import boost_process_priority
from olb.waveoptics.propagators import set_fft_backend
from olb.waveoptics.run import _clip
from olb.waveoptics.turbulence.run import (_resolve_seed, _screen_builder,
                                           _screen_seed, _start_field)
from olb.waveoptics.turbulence.sampling import PRESETS
from olb.waveoptics.turbulence.splitstep import (split_step,
                                                 super_gaussian_boundary)
from validation.terrestrial_mmf_tiptilt.ellipse45 import (APERTURE_M, LAM, TX,
                                                          TILT_COLS, TiltRead,
                                                          make_campaign)
from validation.terrestrial_mmf_tiptilt.run import CN2, HERE
from validation.terrestrial_mmf_tiptilt.wander_split import (moments,
                                                            vacuum_radius)

OUT = "two_end_tilt.npz"
END_COLS = (("cx", "cy", "w_st", "collected")
            + tuple(f"{rx}_{c}" for rx in ("circle", "ellipse")
                    for c in TILT_COLS))
COLS = tuple(f"{end}_{c}" for end in ("far", "near") for c in END_COLS)
GATE_TRIALS = 3
GATE_RTOL = 1e-4
SAVE_EVERY = 20


class _Crop:
    """A stand-in for a stored TrialRecord: the crop of a full receive field."""

    def __init__(self, F, half, context):
        c = F.N // 2
        self.array = F.field[c - half:c + half + 1, c - half:c + half + 1]
        self.size = (2 * half + 1) * F.dx
        self.context = context

    def field(self):
        G = Begin(self.size, LAM, self.array.shape[0])
        G.field = self.array
        return G


def read_end(F, p_ref, r_win, half, context, tilt):
    """The moments, the bucket power and the tilt read of one receive field."""
    collected = Power(_clip(Field.copy(F), APERTURE_M, 0.0)) / p_ref
    cx, cy, w_st, _ = moments(F, r_win)
    return [cx, cy, w_st, collected, *tilt(_Crop(F, half, context))]


def run_cell(camp, path_m, n, fft_backend):
    """Propagate the first n stored trials both ways and read both ends."""
    stored = np.load(os.path.join(camp.root_dir,
                                  "per_trial_ell45.npz"))["collected"]
    out_path = os.path.join(camp.root_dir, OUT)
    rows = (list(np.load(out_path)["rows"]) if os.path.exists(out_path)
            else [])
    grid, plan = camp.grid, camp.plan
    z = np.asarray(plan.z_m, dtype=float)
    z_rev = (plan.z_total_m - z)[::-1]
    mask = super_gaussian_boundary(
        grid.n, PRESETS[camp.preset].boundary_width_frac)
    entropy = _resolve_seed(camp.seed)
    r_win = 3 * vacuum_radius(path_m)
    half = int(np.ceil(1.5 * APERTURE_M / 2 / grid.pixel_m))
    context, tilt = {}, TiltRead()
    previous = set_fft_backend(fft_backend)
    try:
        cdtype = field_dtype(camp.precision)
        build = _screen_builder(camp.screen_generator, grid, camp.L0_m,
                                camp.subharmonics, dtype=cdtype,
                                n_draw=grid.n)
        F0 = _start_field(camp.scenario, grid, LAM, False, dtype=cdtype)
        p_ref = Power(F0)
        t0 = time.time()
        for k in range(len(rows), n):
            screens = [build(_screen_seed(entropy, k, j), plan.r0_m[j])
                       for j in range(z.size)]
            F_far = to_host(split_step(F0, z, screens, plan.z_total_m,
                                       boundary=mask))
            F_near = to_host(split_step(F0, z_rev, screens[::-1],
                                        plan.z_total_m, boundary=mask))
            del screens
            row = (read_end(F_far, p_ref, r_win, half, context, tilt)
                   + read_end(F_near, p_ref, r_win, half, context, tilt))
            rows.append(row)
            if k < GATE_TRIALS:
                rel = abs(row[3] / stored[k] - 1)
                print(f"  gate trial {k}: forward collected {row[3]:.6e}, "
                      f"stored {stored[k]:.6e}, rel {rel:.1e}; reverse "
                      f"collected {row[3 + len(END_COLS)]:.6e}", flush=True)
                if rel > GATE_RTOL:
                    raise SystemExit("IDENTITY GATE FAILED: the regenerated "
                                     "atmosphere is not the stored one.")
            if len(rows) % SAVE_EVERY == 0 or k == n - 1:
                np.savez(out_path, rows=np.array(rows), cols=np.array(COLS),
                         path_m=path_m)
                print(f"  {len(rows)}/{n} trials, {time.time() - t0:.0f} s",
                      flush=True)
    finally:
        set_fft_backend(previous)
    return np.array(rows)


def summary(rows, path_m):
    """The two-end correlations of the G-tilt over the ellipse (x and y)."""
    d = dict(zip(COLS, rows.T))
    out = {"n": int(rows.shape[0])}
    for ax in ("x", "y"):
        tf, tn = d[f"far_ellipse_t{ax}_g"], d[f"near_ellipse_t{ax}_g"]
        cf, cn = d[f"far_c{ax}"], d[f"near_c{ax}"]
        out[ax] = {"r_theta_near_far": float(np.corrcoef(tn, tf)[0, 1]),
                   "r_Ltheta_near_c_far": float(np.corrcoef(tn, cf)[0, 1]),
                   "r_Ltheta_far_c_far": float(np.corrcoef(tf, cf)[0, 1]),
                   "r_Ltheta_far_c_near": float(np.corrcoef(tf, cn)[0, 1]),
                   "rms_theta_far_urad": float(tf.std() * 1e6),
                   "rms_theta_near_urad": float(tn.std() * 1e6),
                   "rms_c_far_m": float(cf.std()),
                   "rms_c_near_m": float(cn.std())}
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--trials", type=int, default=1000)
    ap.add_argument("--preset", default="standard")
    ap.add_argument("--fft-backend", default="cupy")
    ap.add_argument("--cn2", type=float, nargs="*", default=CN2)
    ap.add_argument("--length-km", type=float, nargs="*", default=[1, 10])
    ap.add_argument("--tx", nargs="*", default=list(TX), choices=list(TX))
    args = ap.parse_args()

    boost_process_priority()
    with open(os.path.join(HERE, f"results_ell45_{args.preset}.json")) as f:
        done = json.load(f)
    out_path = os.path.join(HERE, f"results_two_end_ell45_{args.preset}.json")
    results = {}
    if os.path.exists(out_path):
        with open(out_path) as f:
            results = json.load(f)
    for path_m in (L * 1e3 for L in args.length_km):
        for cn2 in args.cn2:
            for tx in args.tx:
                key = f"L{path_m / 1e3:g}km_cn2{cn2:.0e}_tx{tx}"
                if key not in done:
                    continue
                print(key, flush=True)
                camp = make_campaign(cn2, path_m, tx, args.preset, 50,
                                     args.fft_backend)
                rows = run_cell(camp, path_m, args.trials, args.fft_backend)
                results[key] = summary(rows, path_m)
                print(json.dumps(results[key], indent=2), flush=True)
                with open(out_path, "w") as f:
                    json.dump(results, f, indent=2)


if __name__ == "__main__":
    main()
