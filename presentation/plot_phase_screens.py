"""Plot 5 - the phase-screen stack (supports the speckle plot).

Fidelity 2 marches the field through a stack of frozen Kolmogorov phase screens,
one per atmospheric slab. These are the ACTUAL screens the planner built for the
hero downlink: strong and rough near the ground (small Fried parameter r0),
weak and smooth up high. The beam accumulates their wrinkles and arrives as the
speckle of plot 3.
"""
import numpy as np
import matplotlib.pyplot as plt

import presentation.common as C
from olb.waveoptics.turbulence.screens import ScreenFactory

C.use_style()
d = np.load(C.datapath("plotdata.npz"), allow_pickle=True)

r0_all = d["plan_r0_m"]
n = int(d["grid_n"])
pixel = float(d["grid_pixel_m"])
r0_total = float(d["r0_total_m"])

# The screens run top-of-atmosphere -> ground. Show five spanning the range.
idx = [0, 2, 4, 6, 8]
idx = [i for i in idx if i < len(r0_all)]
r0_sel = [r0_all[i] for i in idx]

fac = ScreenFactory(n, pixel)
screens = [fac.make(r0, np.random.default_rng(100 + i)) for i, r0 in zip(idx, r0_sel)]

# Trim to a central window so the fine structure is visible.
w = min(n, 360)
s = (n - w) // 2
screens = [sc[s:s + w, s:s + w] for sc in screens]

# One common colour scale across all screens, set by the strongest (ground) layer.
lim = float(np.percentile(np.abs(screens[-1]), 99))

fig, axes = plt.subplots(1, len(screens), figsize=(3.0 * len(screens), 4.2))
im = None
for ax, sc, r0 in zip(axes, screens, r0_sel):
    im = ax.imshow(sc, cmap="RdBu_r", vmin=-lim, vmax=lim, origin="lower")
    ax.set_title(f"r$_0$ = {r0 * 100:.0f} cm", fontsize=12)
    ax.set_xticks([]); ax.set_yticks([]); ax.grid(False)
    for sp in ax.spines.values():
        sp.set_visible(True); sp.set_color("#cccccc")
    rms = float(np.std(sc))
    ax.text(0.5, -0.06, f"$\\sigma_\\varphi$ = {rms:.1f} rad", transform=ax.transAxes,
            ha="center", va="top", fontsize=10, color="#333333")

# Propagation arrow along the bottom.
fig.subplots_adjust(bottom=0.20, top=0.82, wspace=0.08, right=0.9)
cax = fig.add_axes([0.92, 0.30, 0.012, 0.45])
cb = fig.colorbar(im, cax=cax, ticks=[-lim, 0, lim])
cb.set_label("phase  [rad]")
fig.text(0.5, 0.10, "downlink beam propagates down through the stack",
         ha="center", fontsize=11, color="#333333")
fig.patches.append(plt.matplotlib.patches.FancyArrow(
    0.12, 0.06, 0.76, 0, width=0.006, head_width=0.02, head_length=0.02,
    transform=fig.transFigure, color="#666666", length_includes_head=True))

fig.suptitle(f"Fidelity-2 phase screens for the hero downlink   "
             f"(whole-path r$_0$ = {r0_total * 100:.0f} cm, {len(r0_all)} screens)",
             x=0.5, y=0.95)
out = C.figpath("5_phase_screens.png")
fig.savefig(out)
print("wrote", out)
