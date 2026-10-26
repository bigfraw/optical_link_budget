"""Split the 45 deg fades into the wander, the beam shape and the speckle.

It reads, for each cell, the far end of `two_end_tilt.npz` and
`per_trial_ell45.npz` (ellipse45.py), paired trial for trial. The received
power of the case of record (the ellipse receiver, into the fibre) is

    P = launch_fraction * collected * ratio_ellipse * eta_ellipse.

The aperture (1 inch) is small against the short-term radius w_st, so the
bucket reads the local irradiance of the short-term beam. For a Gaussian
short-term beam exp(-2 |r - c|^2 / w_st^2) of on-axis peak 1 / w_st^2
(Andrews and Phillips Ch. 6, DOI 10.1117/3.626196):

    P_centred = P * exp(2 r^2 / w_st^2)        the wander removed (a perfect
                                               transmit tip-tilt)
    P_vacshape = P_centred * (w_st / w_vac)^2  the short-term spread removed too

What P_vacshape still fades is the speckle on the aperture, which no
tip-tilt removes. The loss quantiles of the three give the share of the
fade that each part explains.

Run locally after copying the npz files from bigfraw:

    <python> -m validation.terrestrial_mmf_tiptilt.analyse_wander
"""

import glob
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from validation.terrestrial_mmf_tiptilt.run import HERE, db
from validation.terrestrial_mmf_tiptilt.wander_split import vacuum_radius


def load(cell_dir):
    """The paired trials of one cell, from two_end_tilt.npz (its far end)."""
    w = np.load(os.path.join(cell_dir, "two_end_tilt.npz"))
    d = {str(c)[4:]: v for c, v in zip(w["cols"], w["rows"].T)
         if str(c).startswith("far_")}
    b = np.load(os.path.join(cell_dir, "per_trial_ell45.npz"))
    n = d["cx"].size
    w_vac = vacuum_radius(float(w["path_m"]))
    p = (b["launch_fraction"] * b["collected"][:n] * b["ratio_ellipse"][:n]
         * b["eta_ellipse"][:n])
    r = np.hypot(d["cx"], d["cy"])
    p_c = p * np.exp(2 * r ** 2 / d["w_st"] ** 2)
    p_v = p_c * (d["w_st"] / w_vac) ** 2
    return {"P": p, "centred": p_c, "vac_shape": p_v, "r": r,
            "w_st": d["w_st"], "w_vac": w_vac,
            "gate": np.max(np.abs(d["collected"] / b["collected"][:n] - 1))}


def quantiles(p):
    loss = db(p)
    return {"mean_power_db": float(db(p.mean())),
            "median_db": float(np.median(loss)),
            "p95_db": float(np.percentile(loss, 95)),
            "p99_db": float(np.percentile(loss, 99))}


def figure(name, c):
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    x = c["r"] / c["w_st"]
    ax[0].plot(x, db(c["P"]), ".", ms=3, alpha=0.5, label="trials")
    xs = np.linspace(0, x.max(), 100)
    ax[0].plot(xs, np.median(db(c["centred"])) + db(np.exp(-2 * xs ** 2)),
               "k-", label="median + Gaussian wander loss")
    ax[0].set_xlabel("wander r / w_st")
    ax[0].set_ylabel("loss into the fibre (dB)")
    ax[0].invert_yaxis()
    ax[0].legend()
    for key, lab in (("P", "total"), ("centred", "wander removed"),
                     ("vac_shape", "wander and spread removed")):
        loss = np.sort(db(c[key]))
        ax[1].semilogy(loss, 1 - np.arange(loss.size) / loss.size, label=lab)
    ax[1].set_xlabel("loss into the fibre (dB)")
    ax[1].set_ylabel("P(loss > x)")
    ax[1].legend()
    fig.suptitle(name)
    fig.tight_layout()
    os.makedirs(os.path.join(HERE, "figures"), exist_ok=True)
    fig.savefig(os.path.join(HERE, "figures", f"wander_split_{name}.png"),
                dpi=120)
    plt.close(fig)


def main():
    out = {}
    for f in sorted(glob.glob(os.path.join(HERE, "campaigns", "ell45_*",
                                           "two_end_tilt.npz"))):
        cell = os.path.dirname(f)
        name = os.path.basename(cell).replace("ell45_", "")
        c = load(cell)
        out[name] = {"n": int(c["r"].size), "gate_max_rel": float(c["gate"]),
                     "wander_rms_m": float(np.sqrt(np.mean(c["r"] ** 2))),
                     "w_st_mean_m": float(c["w_st"].mean()),
                     "w_vac_m": c["w_vac"],
                     **{k: quantiles(c[k])
                        for k in ("P", "centred", "vac_shape")}}
        figure(name, c)
        q = out[name]
        print(f"{name}: n {q['n']}, wander rms {q['wander_rms_m']:.3f} m, "
              f"w_st {q['w_st_mean_m']:.3f} m (vac {q['w_vac_m']:.3f})")
        for k in ("P", "centred", "vac_shape"):
            print(f"  {k:10s} mean {q[k]['mean_power_db']:6.2f}  median "
                  f"{q[k]['median_db']:6.2f}  p95 {q[k]['p95_db']:6.2f}  "
                  f"p99 {q[k]['p99_db']:6.2f} dB")
    with open(os.path.join(HERE, "results_wander_split.json"), "w") as f:
        json.dump(out, f, indent=2)


if __name__ == "__main__":
    main()
