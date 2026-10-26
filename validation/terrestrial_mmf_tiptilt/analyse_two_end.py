"""Plot the tilt at the two ends of the 45 deg link (step 2b figures).

It reads `two_end_tilt.npz` of each cell (two_end_tilt.py) and makes, for
each cell, `figures/tilt_two_ends_<cell>.png`:

  top row:    theta_near against theta_far, x and y, with the 1:1 line
  bottom row: c_far against L theta_near (reciprocity: it must line up) and
              against L theta_far (it must not), x axis

Each panel gives the Pearson r with its bootstrap SE. The tilt is the G-tilt
over the ellipse aperture (the case of record); the Z-tilt r is printed next
to it. Reciprocity: Shapiro, DOI 10.1364/JOSA.61.000492.

    <python> -m validation.terrestrial_mmf_tiptilt.analyse_two_end
"""

import glob
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from validation.terrestrial_mmf_tiptilt.run import HERE


def r_se(a, b, B=1000, seed=1):
    """The Pearson r of a and b, and its bootstrap SE over the trials."""
    idx = np.random.default_rng(seed).integers(0, a.size, (B, a.size))
    rs = [np.corrcoef(a[i], b[i])[0, 1] for i in idx]
    return np.corrcoef(a, b)[0, 1], np.std(rs)


def panel(ax, x, y, xlab, ylab, one_to_one=False):
    r, se = r_se(x, y)
    ax.plot(x, y, ".", ms=3, alpha=0.5)
    if one_to_one:
        lim = np.max(np.abs(np.r_[x, y]))
        ax.plot([-lim, lim], [-lim, lim], "k--", lw=0.8)
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    ax.set_title(f"r = {r:.3f} +/- {se:.3f}")
    return r


def main():
    for f in sorted(glob.glob(os.path.join(HERE, "campaigns", "ell45_*",
                                           "two_end_tilt.npz"))):
        z = np.load(f)
        d = dict(zip([str(c) for c in z["cols"]], z["rows"].T))
        L = float(z["path_m"])
        name = os.path.basename(os.path.dirname(f)).replace("ell45_", "")
        fig, ax = plt.subplots(2, 2, figsize=(10, 9))
        for i, a in enumerate("xy"):
            tf = d[f"far_ellipse_t{a}_g"] * 1e6
            tn = d[f"near_ellipse_t{a}_g"] * 1e6
            panel(ax[0, i], tf, tn, f"theta_far {a} (urad)",
                  f"theta_near {a} (urad)", one_to_one=True)
        cf = d["far_cx"] * 100
        panel(ax[1, 0], L * d["near_ellipse_tx_g"] * 100, cf,
              "L theta_near x (cm)", "c_far x (cm)", one_to_one=True)
        panel(ax[1, 1], L * d["far_ellipse_tx_g"] * 100, cf,
              "L theta_far x (cm)", "c_far x (cm)", one_to_one=True)
        fig.suptitle(f"{name}, n = {cf.size}, G-tilt over the ellipse")
        fig.tight_layout()
        fig.savefig(os.path.join(HERE, "figures",
                                 f"tilt_two_ends_{name}.png"), dpi=120)
        plt.close(fig)
        print(name, f"n {cf.size}")
        for src in ("g", "s"):
            for a in "xy":
                tf = d[f"far_ellipse_t{a}_{src}"]
                tn = d[f"near_ellipse_t{a}_{src}"]
                cfa, cna = d[f"far_c{a}"], d[f"near_c{a}"]
                print(f"  {src} {a}: r(near, far) {np.corrcoef(tn, tf)[0, 1]:+.3f}"
                      f"  r(L th_near, c_far) {np.corrcoef(tn, cfa)[0, 1]:+.3f}"
                      f"  r(L th_far, c_far) {np.corrcoef(tf, cfa)[0, 1]:+.3f}"
                      f"  rms th far/near {tf.std() * 1e6:.2f}/"
                      f"{tn.std() * 1e6:.2f} urad"
                      f"  rms c far/near {cfa.std() * 100:.2f}/"
                      f"{cna.std() * 100:.2f} cm"
                      f"  slope c_far/(L th_near) "
                      f"{np.polyfit(L * tn, cfa, 1)[0]:+.2f}")


if __name__ == "__main__":
    main()
