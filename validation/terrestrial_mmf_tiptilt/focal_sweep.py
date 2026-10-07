"""Post-hoc focal-length sweep of the MMF coupler on the stored 2048 px trials.

The campaigns ran at f = pi (D/2) a_core / (1.12 lambda), about 1.21 m: the
spot just matches the core (the MMF `optimal_focus` rule, removed 2026-10-07
for this reason). That is a spot-to-core match, not a capture optimum. A
multimode light bucket captures more at a shorter f, until the fibre NA
gates the pupil at f = D / (2 NA) (58 mm here). The stored receive fields are
re-coupled at each f with no new propagation (`Campaign.recouple`).

A shift of the core along the axis does NOT help: the received wavefront has
R ~ L = 10 km, so its focus sits f^2 / (R - f) ~ 0.15 mm past the focal
plane, against a depth of focus 2 lambda (f/D)^2 ~ 7 mm (S. A. Self,
DOI 10.1364/AO.22.000658).

    python -m validation.terrestrial_mmf_tiptilt.focal_sweep
"""

import argparse
import json
import os

import numpy as np

from olb.terminal import MMF
from validation.terrestrial_mmf_tiptilt.run import (CN2, CORE_RADIUS_M, HERE,
                                                    NA, make_campaign, stats)

F_M = (1.207, 0.8, 0.5, 0.3, 0.2, 0.12, 0.08)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=2048)
    ap.add_argument("--trials", type=int, default=None)
    ap.add_argument("--workers", default=None)
    ap.add_argument("--f", type=float, nargs="*", default=F_M)
    args = ap.parse_args()
    workers = (None if args.workers is None else
               args.workers if args.workers == "auto" else int(args.workers))
    out = {}
    for cn2 in CN2:
        camp = make_campaign(cn2, "standard", 50, "cupy", args.n)
        n = args.trials or camp.n_stored
        p = np.array([t.collected_power for t in
                      camp.load(n, fields=False).trials])
        etas = {}
        for f in args.f:
            det = MMF(core_radius_m=CORE_RADIUS_M, numerical_aperture=NA,
                      focal_length_m=f)
            etas[f"f={f:g}"] = camp.recouple(det, n_trials=n, workers=workers)
        s = stats(p, etas)
        out[f"{cn2:.0e}"] = s
        print(f"\nCn2 {cn2:.0e}, {n} trials, {args.n} px; loss dB, +/- bootstrap SE")
        print(f"  {'f (m)':8s}{'mean eta':>16s}{'mean':>16s}{'median':>16s}"
              f"{'p95':>16s}{'p99':>16s}")
        for k in etas:
            r = s[k]
            print(f"  {k[2:]:8s}" + "".join(
                f"{r[q][0]:9.3f} +/-{r[q][1]:4.2f}" for q in
                ("mean_eta", "mean_power_db", "median_db", "p95_db", "p99_db")),
                flush=True)
        with open(os.path.join(HERE, f"focal_sweep_n{args.n}.json"), "w") as fh:
            json.dump(out, fh, indent=2)


if __name__ == "__main__":
    main()
