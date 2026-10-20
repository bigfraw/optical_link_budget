"""Split the purely geometric loss out of the stored loss statistics.

The geometric loss is the VACUUM loss from the launched power to the power in
the receive aperture. It has three parts:

- spread_db: the Gaussian beam spread into the receive aperture,
  eta = 1 - exp(-2 a^2 / w(L)^2) (`geometric_loss_term`; Andrews and
  Phillips, DOI 10.1117/3.626196, Ch. 4).
- truncation_db: the launch-aperture clip. The far field of the clipped beam
  has less on-axis gain (Klein and Degnan, DOI 10.1364/AO.13.002134,
  `tx_gaussian_efficiency_term`). The campaign divides by the power AFTER the
  clip, so the clipped power 1 - exp(-2 alpha^2) comes back out.
- grid_aperture_db: the receive mask on the campaign grid is a set of whole
  pixels (`CircAperture`), so its area is not pi a^2. At 5 px across the
  1 inch aperture the mask is 10 percent too large. The beam is flat over the
  aperture, so the vacuum capture scales with the mask area. This is a
  sampling artifact of the run, not physics; it is removed so the remainder
  is the turbulence and the coupling only.

The parts are deterministic, so each quantile of the loss moves by the same
constant. A wave-optics vacuum run (`propagate_scenario`) agrees to 0.06 dB
on the 2 inch family. For each results JSON and each Cn2 the script adds a
"geometric" block:

    spread_db, truncation_db, grid_aperture_db, total_db
    turbulence_bucket_db      collected_mean_db - total_db
    <case>: {stat: value - total_db}   the loss without the geometric part

    python -m validation.terrestrial_mmf_tiptilt.geometric_split
"""

import glob
import json
import os
import re

import numpy as np

from olb.models.gaussian_efficiency import tx_gaussian_efficiency_term
from olb.models.geometric import geometric_loss_term
from olb.waveoptics.field import Begin
from olb.waveoptics.sources import CircAperture

from .run import HERE, LAM, build_scenario, db, refined_grid

LOSS_STATS = ("mean_power_db", "median_db", "p95_db", "p99_db")
CASES = ("no_tt", "tip_tilt", "tip_tilt_gtilt")


def geometric_db(cn2, aperture_in, n_px, preset="standard"):
    """Give the (spread, truncation, grid_aperture) parts in dB."""
    scn, geo = build_scenario(cn2, aperture_in)
    tx = scn.tx_terminal
    alpha = tx.aperture_m / 2 / tx.transmitter.waist_m
    clip = db(1 - np.exp(-2 * alpha ** 2))
    a = scn.rx_terminal.aperture_m / 2
    grid = refined_grid(scn, geo, preset, n_px)[0]
    mask = CircAperture(Begin(grid.size_m, LAM, n_px, np.complex64), a)
    area = np.count_nonzero(mask.field) * (grid.size_m / n_px) ** 2
    return (float(geometric_loss_term(scn, geo).mean_db),
            float(tx_gaussian_efficiency_term(scn, geo).mean_db) - clip,
            db(area / (np.pi * a ** 2)))


def split(path):
    with open(path) as f:
        res = json.load(f)
    aperture_in = 2 if "_ap2in" in path else 1
    n_px = int(re.search(r"_n(\d+)", path).group(1))
    for cn2, r in res.items():
        r.pop("geometric", None)
        spread, trunc, grid_ap = geometric_db(float(cn2), aperture_in, n_px)
        total = spread + trunc + grid_ap
        g = {"spread_db": spread, "truncation_db": trunc,
             "grid_aperture_db": grid_ap, "total_db": total,
             "turbulence_bucket_db": r["collected_mean_db"][0] - total}
        for case in CASES:
            g[case] = {s: r[case][s][0] - total for s in LOSS_STATS}
        r["geometric"] = g
    with open(path, "w") as f:
        json.dump(res, f, indent=2)
    return res


def main():
    for path in sorted(glob.glob(os.path.join(HERE, "results_*.json"))):
        print(os.path.basename(path))
        for cn2, r in split(path).items():
            g = r["geometric"]
            print(f"  Cn2 {cn2}: spread {g['spread_db']:.2f}, truncation "
                  f"{g['truncation_db']:.2f}, grid aperture "
                  f"{g['grid_aperture_db']:+.2f}, total {g['total_db']:.2f} dB;"
                  f" bucket turbulence {g['turbulence_bucket_db']:+.2f} dB;"
                  " no TT minus geometric mean / median / p95 / p99 = "
                  + " / ".join(f"{v:.2f}" for v in g["no_tt"].values()))


if __name__ == "__main__":
    # Self-check: the 1 inch beam is 1.09 m at 10 km, so the aperture keeps
    # 2 a^2 / w^2 = 2.7e-4 (35.7 dB); the 2 inch launch alpha = 1.4 clips the
    # on-axis gain to (1 - exp(-alpha^2))^2 / (1 - exp(-2 alpha^2)), 1.24 dB.
    sp, tr, _ = geometric_db(1e-15, 2, 2048)
    assert abs(tr - 1.24) < 0.01 and 17 < sp < 18
    assert abs(geometric_db(1e-15, 1, 2048)[0] - 35.64) < 0.01
    main()
