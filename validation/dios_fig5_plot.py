'''Plot the Dios Fig. 5 replication: fidelity 1 curves, fidelity 2 points.

Run from the repo root, after dios_fig5_replication.py:
    python -m validation.dios_fig5_plot
'''

import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(__file__)
with open(os.path.join(HERE, "dios_fig5_replication_results.json")) as f:
    R = json.load(f)

fig, ax = plt.subplots(figsize=(7.0, 5.5))
colors = {"90.0": "tab:blue", "30.0": "tab:red"}

for elev, c in colors.items():
    rows = R["fid1"][elev]
    w0 = np.array([r["w0_m"] for r in rows])
    ax.loglog(w0, [r["sigma2_chi_slant"] for r in rows], "-", color=c,
              label=f"fid 1 (slant-corrected), {float(elev):.0f} deg")
    ax.loglog(w0, [r["sigma2_chi"] for r in rows], "--", color=c, alpha=0.5,
              label=f"fid 1 as shipped, {float(elev):.0f} deg")
    ax.loglog(w0, [r["sigma2_chi_on_axis"] for r in rows], ":", color=c,
              alpha=0.7, label=f"on-axis only, {float(elev):.0f} deg")
    rows2 = R["fid2"][elev]
    ax.loglog([r["w0_m"] for r in rows2],
              [r["sigma2_chi"] for r in rows2], "*", color=c, markersize=14,
              markeredgecolor="k",
              label=f"fid 2 wave optics, {float(elev):.0f} deg")

ax.set_xlabel("$W_0$ [m]")
ax.set_ylabel(r"$\sigma_\chi^2$")
ax.set_ylim(1e-3, 1.5)
ax.set_title("Dios 2004 Fig. 5 replication (GEO, 0.84 um, HV57)\n"
             "olb coupled-flux kernels vs olb split-step reciprocity")
ax.grid(True, which="both", alpha=0.3)
ax.legend(fontsize=8, loc="upper left")
fig.tight_layout()
out = os.path.join(HERE, "dios_fig5_replication.png")
fig.savefig(out, dpi=150)
print(f"wrote {out}")
