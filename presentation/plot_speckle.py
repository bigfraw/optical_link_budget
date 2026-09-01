"""Plot 3 - the downlink beam at the telescope pupil (one atmosphere).

What the ground telescope actually receives, for a single turbulent atmosphere:
the intensity has broken into bright and dark speckle, and the phase is a
wrinkled, wrapped mess instead of a flat wavefront. It is this that destroys the
single-mode fibre coupling. The white circle is the 0.7 m aperture.
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

import presentation.common as C

C.use_style()
d = np.load(C.datapath("plotdata.npz"), allow_pickle=True)

ext = d["speckle_extent"]
Ap = float(d["rx_aperture_m"])
I = d["I_turb"]
Ph = d["Ph_turb"]

# Crop to a little beyond the aperture.
half = Ap * 0.85
x = np.linspace(ext[0], ext[1], I.shape[1])
y = np.linspace(ext[2], ext[3], I.shape[0])
cx = (x >= -half) & (x <= half)
cy = (y >= -half) & (y <= half)
crop = np.ix_(cy, cx)
crop_ext = [-half, half, -half, half]

X, Y = np.meshgrid(x[cx], y[cy])
disk = (X ** 2 + Y ** 2) <= (Ap / 2) ** 2

Ic = I[crop]
Ic = Ic / Ic[disk].mean()
Phc = Ph[crop]
si = Ic[disk].var() / Ic[disk].mean() ** 2

fig, (axI, axP) = plt.subplots(1, 2, figsize=(11, 5.3))

mI = axI.imshow(Ic, extent=crop_ext, origin="lower", cmap="magma",
                vmin=0, vmax=np.percentile(Ic[disk], 99.5))
axI.set_title("Intensity", fontsize=13)
axI.text(0.03, 0.03, rf"$\sigma_I^2$ = {si:.2f}", transform=axI.transAxes,
         color="white", fontsize=11, va="bottom",
         bbox=dict(boxstyle="round,pad=0.25", fc="black", alpha=0.4, ec="none"))
cI = fig.colorbar(mI, ax=axI, fraction=0.046, pad=0.03)
cI.set_label("I / ⟨I⟩")

mP = axP.imshow(Phc, extent=crop_ext, origin="lower", cmap="twilight",
                vmin=-np.pi, vmax=np.pi)
axP.set_title("Phase", fontsize=13)
cP = fig.colorbar(mP, ax=axP, fraction=0.046, pad=0.03,
                  ticks=[-np.pi, 0, np.pi])
cP.ax.set_yticklabels(["−π", "0", "π"])
cP.set_label("wavefront phase  [rad]")

for ax in (axI, axP):
    ax.add_patch(Circle((0, 0), Ap / 2, fill=False, ec="white", lw=1.8, ls="--"))
    ax.set_xlabel("x  [m]"); ax.set_xticks([-0.5, 0, 0.5]); ax.set_yticks([-0.5, 0, 0.5])
    ax.grid(False)
axI.set_ylabel("y  [m]")

fig.suptitle(f"Downlink beam at the telescope pupil — one atmosphere "
             f"({C.HERO_ELEVATION_DEG:.0f}° elevation)", x=0.5)
fig.tight_layout(rect=[0, 0, 1, 0.95])
out = C.figpath("3_speckle.png")
fig.savefig(out)
print("wrote", out)
