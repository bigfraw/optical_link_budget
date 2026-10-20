"""The focused field on the MMF tip at F_DESIGN, for a few stored trials.

Each panel is the focal-plane intensity of one stored receive field: the
field is Fourier-upsampled on the crop (the same factor as the coupling read,
`mmf._resolve_auto_factor`), clipped at the 1 inch aperture, gated at the
fibre NA (rho <= f NA, Snyder and Love, DOI 10.1007/978-1-4613-2813-1) and
focused by one FFT (Fraunhofer, Goodman, ISBN 978-0974707723). The pupil is
zero-padded PAD times first, so the focal pixel is PAD times finer; that adds
no information. The circle is the 105 um core. The top row is the open loop,
the bottom row has the G-tilt removed (as `run.CropEta`).

The trials are the min, p10, median, p90 and max of the open-loop eta in
per_trial_<F_TAG>.npz (run.py writes it), or the first five without it.

    python -m validation.terrestrial_mmf_tiptilt.fibre_tip --n 2048
"""

import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from olb.terminal import TipTilt
from olb.waveoptics.field import Begin
from olb.waveoptics.mmf import (_resolve_auto_factor, _upsampled_field,
                                focal_intensity, mmf_coupling_efficiency)
from olb.waveoptics.turbulence.run import _PostCorrector
from validation.terrestrial_mmf_tiptilt.run import (APERTURE_M, CN2,
                                                    CORE_RADIUS_M, F_DESIGN,
                                                    F_TAG, HERE, NA,
                                                    make_campaign)

PAD = 4                  # pupil zero-pad factor: a finer focal pixel only
HALF_UM = 120            # the plotted half-width at the tip, in um
DECADES = 3              # the log colour range


def tip_image(F):
    """Give (intensity, focal pixel in m, eta) of one crop Field."""
    m = _resolve_auto_factor("auto", F, APERTURE_M, CORE_RADIUS_M, F_DESIGN,
                             NA, None)
    U = _upsampled_field(F, m)
    U.field = np.where(U.mgrid_Rsquared <= (APERTURE_M / 2) ** 2, U.field, 0)
    n = U.N * PAD
    P = Begin(U.siz * PAD, U.lam, n, dtype=U.field.dtype)
    o = (n - U.N) // 2
    P.field[o:o + U.N, o:o + U.N] = U.field
    I, dx = focal_intensity(P, F_DESIGN, numerical_aperture=NA)
    eta = mmf_coupling_efficiency(F, APERTURE_M, CORE_RADIUS_M, F_DESIGN,
                                  numerical_aperture=NA, upsample="auto")
    return I / (np.abs(U.field) ** 2).sum(), dx, eta


def pick_rows(camp):
    path = os.path.join(camp.root_dir, f"per_trial_{F_TAG}.npz")
    if not os.path.exists(path):
        return list(range(5)), ["0", "1", "2", "3", "4"]
    eta = np.load(path)["eta_open"]
    q = (0, 10, 50, 90, 100)
    rows = [int(np.argsort(eta)[min(eta.size - 1, k * eta.size // 100)])
            for k in q]
    return rows, ["min", "p10", "median", "p90", "max"]


def plot(cn2, n_px):
    camp = make_campaign(cn2, "standard", 50, "cupy", n_px)
    rows, labels = pick_rows(camp)
    gt = _PostCorrector(camp.patch, [TipTilt()], APERTURE_M, 0.0, "gtilt")
    fig, axes = plt.subplots(2, len(rows), figsize=(3.1 * len(rows), 6.4),
                             constrained_layout=True)
    imgs = []
    for j, row in enumerate(rows):
        F = camp.field(row)
        G = camp.field(row)
        G.field = gt.correct(F.field)
        imgs += [(0, j, row, labels[j]) + tip_image(F),
                 (1, j, row, labels[j]) + tip_image(G)]
    vmax = max(im[4].max() for im in imgs)
    t = np.linspace(0, 2 * np.pi, 200)
    for i, j, row, lab, I, dx, eta in imgs:
        ax = axes[i, j]
        h = int(HALF_UM * 1e-6 / dx)
        c = I.shape[0] // 2
        cut = I[c - h:c + h + 1, c - h:c + h + 1]
        im = ax.imshow(np.log10(np.maximum(cut, vmax * 10 ** -DECADES)),
                       extent=[-h * dx * 1e6, h * dx * 1e6] * 2,
                       origin="lower", cmap="inferno",
                       vmin=np.log10(vmax) - DECADES, vmax=np.log10(vmax))
        ax.plot(CORE_RADIUS_M * 1e6 * np.cos(t),
                CORE_RADIUS_M * 1e6 * np.sin(t), "c-", lw=0.8)
        ax.set_title(f"{'open' if i == 0 else 'G-tilt'}  trial {row} ({lab})"
                     f"\neta = {eta:.3f}", fontsize=9)
        ax.set_xlabel("x (um)")
        if j == 0:
            ax.set_ylabel("y (um)")
    fig.colorbar(im, ax=axes, shrink=0.8,
                 label="log10 intensity / aperture power per px")
    fig.suptitle(f"MMF tip, 10 km, Cn2 = {cn2:.0e}, f = {F_DESIGN*1e3:.2f} mm, "
                 f"NA {NA}, {camp.grid.n} px grid; circle = 105 um core")
    os.makedirs(os.path.join(HERE, "figures"), exist_ok=True)
    out = os.path.join(HERE, "figures",
                       f"fibre_tip_cn2{cn2:.0e}_n{n_px}_{F_TAG}.png")
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print(out, flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=2048)
    ap.add_argument("--cn2", type=float, nargs="*", default=CN2)
    args = ap.parse_args()
    for cn2 in args.cn2:
        plot(cn2, args.n)


if __name__ == "__main__":
    main()
