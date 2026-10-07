"""2048 px against 4096 px: each statistic, its bootstrap SE and the difference.

It reads the per-trial arrays that run.summarise saves (per_trial.npz in each
campaign root). The two grids draw DIFFERENT atmospheres (the screens sit on
different pixels), so the two samples are independent: each is resampled on
its own, and the SE of a difference is the SE of the bootstrap differences.

    python -m validation.terrestrial_mmf_tiptilt.compare_grids --cn2 1e-15
"""

import argparse
import os

import numpy as np

from validation.terrestrial_mmf_tiptilt.run import HERE, db

B = 2000


def load(cn2, n_px):
    root = os.path.join(HERE, "campaigns",
                        f"cn2{cn2:.0e}_standard_focal_n{n_px}")
    return dict(np.load(os.path.join(root, "per_trial.npz")))


def metrics(d, i):
    """The statistics of one resample i (trial indices)."""
    p = d["collected"][i]
    out = {"collected mean": db(p.mean())}
    for name, key in (("no TT", "eta_open"), ("TT pkg", "eta_tt"),
                      ("TT G-tilt", "eta_gtilt")):
        eta = d[key][i]
        loss = db(p * eta)
        out[f"{name} mean eta"] = eta.mean()
        out[f"{name} mean dB"] = db((p * eta).mean())
        out[f"{name} median"] = np.median(loss)
        out[f"{name} p95"] = np.percentile(loss, 95)
        out[f"{name} p99"] = np.percentile(loss, 99)
    out["TT gain dB"] = 10 * np.log10(d["eta_tt"][i].mean()
                                      / d["eta_open"][i].mean())
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cn2", type=float, default=1e-15)
    args = ap.parse_args()
    a, b = load(args.cn2, 4096), load(args.cn2, 2048)
    rng = np.random.default_rng(1)
    ia = rng.integers(0, a["collected"].size, (B, a["collected"].size))
    ib = rng.integers(0, b["collected"].size, (B, b["collected"].size))
    va = metrics(a, np.arange(a["collected"].size))
    vb = metrics(b, np.arange(b["collected"].size))
    ra = [metrics(a, i) for i in ia]
    rb = [metrics(b, i) for i in ib]
    print(f"Cn2 {args.cn2:.0e}: {a['collected'].size} / {b['collected'].size}"
          f" trials, {B} bootstrap resamples. Loss in dB, positive.")
    print(f"{'statistic':20s}{'4096 px':>18s}{'2048 px':>18s}"
          f"{'2048 - 4096':>18s}{'z':>7s}")
    for k in va:
        sa = np.std([r[k] for r in ra])
        sb = np.std([r[k] for r in rb])
        dd = np.std([y[k] - x[k] for x, y in zip(ra, rb)])
        delta = vb[k] - va[k]
        print(f"{k:20s}{va[k]:10.3f} +/-{sa:5.3f}{vb[k]:10.3f} +/-{sb:5.3f}"
              f"{delta:+10.3f} +/-{dd:5.3f}{delta / dd:+7.1f}")


if __name__ == "__main__":
    main()
