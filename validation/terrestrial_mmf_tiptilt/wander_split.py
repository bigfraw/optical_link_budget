"""Collect the beam wander of the stored 45 deg trials (step 1b).

THE QUESTION. How much of each fade of the 45 deg terminal is the beam
WANDER (which a transmit tip-tilt can remove) and how much is the speckle on
the aperture (which it cannot)? See NEXT_tx_tiptilt.md, step 1b.

THE ROUTE. The stored patch is 19 mm in radius and the 10 km beam is 0.53 m
in radius, so the stored fields do not hold the beam centroid. The seed
regenerates each atmosphere bit for bit on the same backend, so this script
propagates each stored trial again with `propagate_turbulent_field` (the
same seed, trial index, grid, plan and run options as the campaign) and it
keeps the receive-plane moments. No package change: the runner and every
stored key stay the same.

THE IDENTITY GATE. The first trials must give the STORED collected power
again (to the float32 rounding). If not, the atmospheres differ and the
pairing is invalid, so the script stops.

THE MOMENTS of each trial, inside a window of radius 3 W(L) about the axis
(W(L) the vacuum Gaussian radius, so the light that the turbulence scatters
to the grid edge does not pull the centroid):

    cx, cy  the power centroid, in m
    w_st    the short-term radius sqrt(2 <|r - c|^2>), in m (a Gaussian
            exp(-2 r^2 / W^2) gives <r^2> = W^2 / 2)
    p_win   the power in the window, against the start power
    p_grid  the power on the grid, against the start power
    collected  the bucket power (the gate value)

Beam radius: Kogelnik and Li, DOI 10.1364/AO.5.001550. Wander and
short-term spread: Andrews and Phillips Ch. 6, DOI 10.1117/3.626196.

Run from the repository root (on bigfraw, the GPU venv), then read the
per-campaign file `wander_ell45.npz`:

    <python> -m validation.terrestrial_mmf_tiptilt.wander_split --cn2 1e-15 --fft-backend cupy
"""

import argparse
import inspect
import json
import os
import time

import numpy as np

from olb.turbulence.andrews.wander import beam_wander_variance
from olb.waveoptics.field import Field, Power, field_dtype
from olb.waveoptics.priority import boost_process_priority
from olb.waveoptics.run import _clip
from olb.waveoptics.turbulence.run import (_start_field,
                                           propagate_turbulent_field)
from validation.terrestrial_mmf_tiptilt.ellipse45 import (APERTURE_M, LAM, TX,
                                                          WAIST_M,
                                                          make_campaign)
from validation.terrestrial_mmf_tiptilt.run import CN2, HERE

OUT = "wander_ell45.npz"
COLS = ("cx", "cy", "w_st", "p_win", "p_grid", "collected")
GATE_TRIALS = 3
GATE_RTOL = 1e-4
SAVE_EVERY = 20


def vacuum_radius(path_m):
    """The vacuum 1/e^2 radius of the launch waist at the range path_m."""
    z_r = np.pi * WAIST_M ** 2 / LAM
    return WAIST_M * np.hypot(1.0, path_m / z_r)


def moments(F, r_win):
    """The centroid, the short-term radius and the power inside r_win."""
    Y, X = F.mgrid_cartesian
    w = np.abs(F.field).astype(np.float64) ** 2
    w[X ** 2 + Y ** 2 > r_win ** 2] = 0.0
    s = w.sum()
    cx, cy = (w * X).sum() / s, (w * Y).sum() / s
    w_st = np.sqrt(2 * (w * ((X - cx) ** 2 + (Y - cy) ** 2)).sum() / s)
    return cx, cy, w_st, s


def run_cell(camp, path_m, n, fft_backend):
    """Propagate the first n stored trials again and keep their moments."""
    stored = np.load(os.path.join(camp.root_dir, "per_trial_ell45.npz"))
    stored = stored["collected"]
    out_path = os.path.join(camp.root_dir, OUT)
    rows = (list(np.load(out_path)["rows"]) if os.path.exists(out_path)
            else [])
    sig = inspect.signature(propagate_turbulent_field).parameters
    kw = {k: v for k, v in camp._runner_kwargs().items() if k in sig}
    kw["fft_backend"] = fft_backend
    scn = camp.scenario
    F0 = _start_field(scn, camp.grid, LAM, False,
                      dtype=field_dtype(camp.precision))
    p_ref = Power(F0)
    r_win = 3 * vacuum_radius(path_m)
    t0 = time.time()
    for k in range(len(rows), n):
        F, _, _ = propagate_turbulent_field(scn, camp.geometry, trial=k,
                                            grid=camp.grid, plan=camp.plan,
                                            **kw)
        collected = Power(_clip(Field.copy(F), APERTURE_M, 0.0)) / p_ref
        cx, cy, w_st, s = moments(F, r_win)
        p_grid = Power(F) / p_ref
        p_win = p_grid * s / np.sum(np.abs(F.field).astype(np.float64) ** 2)
        rows.append((cx, cy, w_st, p_win, p_grid, collected))
        if k < GATE_TRIALS:
            rel = abs(collected / stored[k] - 1)
            print(f"  gate trial {k}: collected {collected:.6e}, stored "
                  f"{stored[k]:.6e}, rel {rel:.1e}", flush=True)
            if rel > GATE_RTOL:
                raise SystemExit("IDENTITY GATE FAILED: the regenerated "
                                 "atmosphere is not the stored one.")
        if len(rows) % SAVE_EVERY == 0 or k == n - 1:
            np.savez(out_path, rows=np.array(rows), cols=np.array(COLS),
                     r_win_m=r_win, w_vac_m=vacuum_radius(path_m))
            print(f"  {len(rows)}/{n} trials, {time.time() - t0:.0f} s",
                  flush=True)
    return np.array(rows)


def summary(rows, path_m, cn2):
    """The rms wander against the model, and the short-term radius."""
    d = dict(zip(COLS, rows.T))
    r = np.hypot(d["cx"], d["cy"])
    # Andrews and Phillips Ch. 6 Eq. (93), DOI 10.1117/3.626196, L0 = 25 m.
    model = np.sqrt(beam_wander_variance(WAIST_M, LAM, path_m, cn2,
                                         spectrum="exponential", L0=25.0))
    return {"n": int(r.size), "wander_rms_m": float(np.sqrt(np.mean(r ** 2))),
            "wander_rms_model_m": float(model),
            "w_st_mean_m": float(d["w_st"].mean()),
            "w_vac_m": float(vacuum_radius(path_m)),
            "p_win_mean": float(d["p_win"].mean()),
            "corr_collected_r": float(np.corrcoef(d["collected"], r)[0, 1])}


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
    out_path = os.path.join(HERE, f"results_wander_ell45_{args.preset}.json")
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
                results[key] = summary(rows, path_m, cn2)
                print(json.dumps(results[key], indent=2), flush=True)
                with open(out_path, "w") as f:
                    json.dump(results, f, indent=2)


if __name__ == "__main__":
    main()
