"""The retro bracket: does the shared atmosphere move the retro fade?

A space retro link has an up leg and a return leg. A fidelity-2 retro can pair
the two legs in three ways, and this script measures all three on the SAME
stored trials, with NO new propagation:

- SAME: the up leg and the return leg of trial k read the SAME atmosphere.
  This is the fully correlated limit (point-ahead angle 0).
- PA: the up leg reads the point-ahead window of trial k and the return leg
  reads the beacon window of trial k. This is the true geometry: the return
  comes down along the apparent direction, and the up leg leaves at the
  point-ahead angle 2 v_perp / c.
- INDEP: the up leg and the return leg read DIFFERENT atmospheres. The retro
  fade in dB is the sum of the two leg fades, so its distribution is the
  convolution of the two measured distributions (all n^2 pairs).

If INDEP and PA agree inside the noise, a fidelity-2 retro can read one plain
downlink campaign with independent legs. If they differ, it needs the
point-ahead record.

THE LEGS. The stored `base` campaigns of `validation/waveoptics_pointahead/`
hold the uncorrected ground field of the beacon window and of every
point-ahead window. The up leg is the Shapiro reciprocity overlap of that
field with the ground transmit mode, over its own vacuum overlap
(`Campaign.uplink_overlaps`; J. H. Shapiro, DOI 10.1364/JOSA.61.000492). A
small corner cube is a point receiver in the uplink far field, so the overlap
IS the power it catches. The return leg is the beacon field at the ground: the
0.7 m single-mode fibre coupling (`Campaign.recouple`) and the 0.7 m bucket
(`Campaign.recollect`). Two launches: the hero full-aperture waist (0.35 m,
monostatic, near the fibre mode) and a small 0.06 m beam-director waist on the
same aperture.

THE STATISTICS. For each pairing: the loss of the mean power, and the loss at
the 5 and 1 percent power quantiles (the p5 and p1 fades), each in dB, and the
correlation of the two leg losses in dB. A paired bootstrap over the trials
gives the standard error of each difference against INDEP.

Run it on `bigfraw`, where the campaigns live (see
`validation/waveoptics_pointahead/README.md` for the launch rules):

    python -m validation.retro_bracket.retro_bracket

The smoke run on the laptop reads the stored smoke campaign:

    $env:OLB_WAVEOPTICS_POINTAHEAD_ROOT = "validation/waveoptics_pointahead/campaigns_smoke"
    python -m validation.retro_bracket.retro_bracket --elevations 20 \
        --n-trials 8 --block-size 4 --fft-backend numpy --preset rapid --n-boot 20 \
        --fixed-arcsec 5
"""

import argparse
import json
import os
import time
from dataclasses import replace

import numpy as np

from olb.terminal import SMF
from validation.waveoptics_pointahead import waveoptics_pointahead as pa

HERE = os.path.dirname(os.path.abspath(__file__))
ARCSEC = np.pi / 180 / 3600

# The small launch: a beam-director waist on the same 0.7 m ground aperture.
SMALL_WAIST_M = 0.06
LAUNCHES = ("hero 0.35 m", f"small {SMALL_WAIST_M} m")
QUANTILES = (0.05, 0.01)
BOOT_SEED = 20260929


def db(x):
    """Give -10 log10(x): a power ratio as a loss in dB."""
    return -10.0 * np.log10(x)


def stats(power):
    """Give the mean loss and the quantile fades of a power sample, in dB."""
    p = np.asarray(power, float).ravel()
    return [db(p.mean())] + [db(np.quantile(p, q)) for q in QUANTILES]


def indep_power(u, d):
    """Give every product u_i d_j: the retro power of independent legs."""
    return np.multiply.outer(u, d)


def read_legs(camp, n, workers):
    """Read both legs of every stored trial, post hoc.

    Returns:
        (up, up_pa, down): up is (n, n_launch), up_pa is
        (n, n_angles, n_launch), down is a dict of (n,) arrays.
    """
    g = camp.scenario.ground
    small = replace(g, transmitter=replace(g.transmitter,
                                           waist_m=SMALL_WAIST_M))
    ov = camp.uplink_overlaps(None, grounds=[g, small], n_trials=n,
                              workers=workers)
    a, o = g.aperture_m, g.obscuration_ratio
    down = {
        "SMF": camp.recouple(SMF(), aperture_m=a, obscuration_ratio=o,
                             n_trials=n, workers=workers),
        "bucket": camp.recollect(aperture_m=a, obscuration_ratio=o,
                                 n_trials=n, workers=workers),
    }
    # The bucket power holds no vacuum reference. The comparison is a
    # DIFFERENCE between pairings, so a constant scale cancels; the mean
    # makes the numbers readable.
    down["bucket"] = down["bucket"] / down["bucket"].mean()
    return ov["none"], ov["none_pa"], down


def bracket(u0, u_pa, d, n_boot, rng):
    """Compare every pairing against INDEP on one launch and one receiver.

    Args:
        u0:   (n,) the up leg on the beacon window (the SAME pairing).
        u_pa: (n, n_angles) the up leg on each point-ahead window.
        d:    (n,) the return leg.

    Returns:
        A dict: "indep" holds the absolute stats; each pairing holds the
        correlation, the stats minus INDEP, and the bootstrap SE of each.
    """
    n = len(d)
    pairs = {"same": u0}
    for i in range(u_pa.shape[1]):
        pairs[f"pa{i}"] = u_pa[:, i]

    def one(idx):
        ref = np.array(stats(indep_power(u0[idx], d[idx])))
        return ref, {k: np.array(stats(u[idx] * d[idx])) - ref
                     for k, u in pairs.items()}

    ref, diff = one(np.arange(n))
    boots = [one(rng.integers(0, n, n))[1] for _ in range(n_boot)]
    out = {"indep": ref.tolist()}
    for k, u in pairs.items():
        se = np.std([b[k] for b in boots], axis=0, ddof=1)
        out[k] = {"rho_db": float(np.corrcoef(db(u), db(d))[0, 1]),
                  "delta_db": diff[k].tolist(), "se_db": se.tolist()}
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--elevations", type=float, nargs="+", default=[30, 20])
    ap.add_argument("--n-trials", type=int, default=1000)
    ap.add_argument("--block-size", type=int, default=50)
    ap.add_argument("--fft-backend", default="cupy")
    ap.add_argument("--preset", default=pa.PRESET)
    ap.add_argument("--n-boot", type=int, default=200)
    ap.add_argument("--fixed-arcsec", type=float, nargs="*",
                    default=list(pa.FIXED_ARCSEC))
    ap.add_argument("--workers", default=None)
    ap.add_argument("--run-workers", default=None,
                    help="compute the missing trials first with this pool "
                         "(an int or \"auto\"); unset computes none")
    args = ap.parse_args()
    pa.PRESET = args.preset
    pa.FIXED_ARCSEC = tuple(args.fixed_arcsec)
    workers = (None if args.workers is None else
               args.workers if args.workers == "auto" else int(args.workers))

    lines, record = [], {}

    def say(s=""):
        print(s, flush=True)
        lines.append(s)

    say("retro bracket: pairing minus INDEP, in dB (positive = more loss)")
    say(f"preset {args.preset}, seed {pa.SEED}, {args.n_boot} bootstraps")
    rng = np.random.default_rng(BOOT_SEED)
    for el in args.elevations:
        camp, _warn = pa.campaign_of(el, "base", args)
        if args.run_workers is not None and camp.n_stored < args.n_trials:
            rw = (args.run_workers if args.run_workers == "auto"
                  else int(args.run_workers))
            t0 = time.perf_counter()
            camp.run(args.n_trials, workers=rw)
            say(f"el {el:g}: computed to {camp.n_stored} trials in "
                f"{time.perf_counter() - t0:.0f} s")
        n = min(args.n_trials, camp.n_stored)
        if n < 2:
            raise SystemExit(f"el {el:g}: no stored trials in the campaign store")
        up, up_pa, down = read_legs(camp, n, workers)
        # The identity of the record: angle 0 IS the beacon window.
        assert np.allclose(up_pa[:, 0, :], up, rtol=1e-6), "angle 0 != beacon"
        angles = np.asarray(camp.point_ahead_rad) / ARCSEC
        names = ["same"] + [f"pa{i}" for i in range(len(angles))]
        labels = ["SAME (shared)"] + [
            f"PA {a:5.2f}\"" + (" geom" if i == 1 else "")
            for i, a in enumerate(angles)]
        for j, launch in enumerate(LAUNCHES):
            for rx, d in down.items():
                key = f"el{el:g} | {launch} | return {rx}"
                res = bracket(up[:, j], up_pa[:, :, j], d, args.n_boot, rng)
                record[key] = {"n_trials": n,
                               "angles_arcsec": angles.tolist(), **res}
                say()
                say(f"{key}   ({n} trials)")
                m, p5, p1 = res["indep"]
                say(f"  INDEP (ref)          mean {m:6.2f}  p5 {p5:6.2f}  "
                    f"p1 {p1:6.2f}")
                for k, lab in zip(names, labels):
                    r = res[k]
                    cells = "  ".join(
                        f"{nm} {v:+5.2f}+/-{s:4.2f}" for nm, v, s in zip(
                            ("mean", "p5", "p1"), r["delta_db"], r["se_db"]))
                    say(f"  {lab:<20} rho {r['rho_db']:+.2f}  {cells}")

    tag = "_".join(f"el{e:g}" for e in args.elevations)
    with open(os.path.join(HERE, f"retro_bracket_{tag}.log"), "w") as f:
        f.write("\n".join(lines) + "\n")
    with open(os.path.join(HERE, f"retro_bracket_{tag}.json"), "w") as f:
        json.dump(record, f, indent=1)


if __name__ == "__main__":
    main()
