"""Plot 1 - the fidelity ladder, at two elevations.

One downlink, three ways to model the turbulent fibre coupling, at zenith and at
a 30 deg mid-pass. Fidelity 0 is a single number, fidelity 1 is an analytic
distribution, fidelity 2 is a sampled field solve. The x-axis is the received
power in the fibre (a common deterministic baseline is used, so the three rungs
differ only in how they model the atmosphere). The mean and the 90% availability
level are marked. Cost rises and assumptions fall as you climb.
"""
import numpy as np
import matplotlib.pyplot as plt

import presentation.common as C
from olb.links.downlink import downlink_budget

C.use_style()
d = np.load(C.datapath("plotdata.npz"), allow_pickle=True)

scn = C.downlink()
ELEVS = [90, 30]
COL_TITLE = {90: "90°  ·  zenith", 30: "30°  ·  mid-pass"}
ROW_TAG = {0: "Fidelity 0 · analytic  (a single number)",
           1: "Fidelity 1 · statistical  (FAST Monte-Carlo)",
           2: "Fidelity 2 · wave optics  (split-step field solve)"}


def received_baseline(elev):
    """Received power (dBm) with the turbulent coupling term removed."""
    b = downlink_budget(scn, C.geometry(elev), fidelity=0)
    det = sum(float(np.asarray(t.mean_db)) for t in b.terms if t.category != "coupling")
    return float(b.tx_power_dbm) - det


fig, axes = plt.subplots(3, 2, figsize=(13, 8.4), sharex="col")

for j, elev in enumerate(ELEVS):
    base = received_baseline(elev)
    c0 = float(d[f"fid0_coup_{elev}"])
    rx0 = base - c0
    rx1 = base - d[f"fid1_coup_{elev}"]
    rx2 = base - d[f"fid2_coup_{elev}"]
    ms = float(d[f"wall_{elev}"]) / 300.0 * 1000.0

    lo = np.percentile(np.concatenate([rx1, rx2]), 0.5)
    hi = max(rx0, np.percentile(np.concatenate([rx1, rx2]), 99.5))
    bins = np.linspace(lo - 1, hi + 1, 60)

    def avail90(x):   # received power exceeded 90% of the time
        return float(np.percentile(x, 10))

    # Row 0 - fidelity 0, a point
    ax = axes[0, j]
    ax.axvline(rx0, color=C.FID_COLORS[0], lw=2.5)
    ax.plot([rx0], [0.5], "o", ms=10, color=C.FID_COLORS[0], zorder=5)
    ax.set_ylim(0, 1); ax.set_yticks([])
    ax.text(rx0, 0.72, f" {rx0:.1f} dBm", color=C.FID_COLORS[0], fontweight="bold",
            va="center", fontsize=11)
    ax.text(0.015, 0.9, ROW_TAG[0], transform=ax.transAxes, ha="left", va="top",
            fontweight="bold", fontsize=9.5, color="#333333")
    ax.text(0.015, 0.66, "no fade  ·  < 1 ms", transform=ax.transAxes, ha="left",
            va="top", color=C.GREY, fontsize=9)

    # Rows 1, 2 - distributions
    for row, rx, color, cost in [
        (1, rx1, C.FID_COLORS[1], "~10 ms"),
        (2, rx2, C.FID_COLORS[2], f"~{ms:.0f} ms / trial"),
    ]:
        ax = axes[row, j]
        m, a90 = float(np.mean(rx)), avail90(rx)
        ax.hist(rx, bins=bins, density=True, color=color, alpha=0.65,
                edgecolor="white", linewidth=0.3)
        ax.axvline(m, color=color, lw=2)
        ax.axvline(a90, color=C.RED, lw=1.6, ls="--")
        ax.set_yticks([]); ax.set_ylabel("density")
        ax.text(0.015, 0.9, ROW_TAG[row], transform=ax.transAxes, ha="left", va="top",
                fontweight="bold", fontsize=9.5, color="#333333")
        ax.text(0.015, 0.72,
                f"mean {m:.1f} dBm  ·  90% avail {a90:.0f} dBm  ·  {cost}",
                transform=ax.transAxes, ha="left", va="top", color=C.GREY, fontsize=9)

    axes[0, j].set_title(COL_TITLE[elev], fontsize=13)
    axes[2, j].set_xlabel("Received power in fibre  [dBm]")

fig.suptitle("The fidelity ladder — one downlink, three models, two elevations", y=0.995)
fig.text(0.5, 0.005,
         "0.7 m telescope → single-mode fibre, uncompensated · 500 km LEO · 1550 nm   "
         "·   solid = mean, red dashed = 90% availability",
         ha="center", color=C.GREY, fontsize=9)
fig.tight_layout(rect=[0.02, 0.02, 1, 0.94])
out = C.figpath("1_fidelity_ladder.png")
fig.savefig(out)
print("wrote", out)
