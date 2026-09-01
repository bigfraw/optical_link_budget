"""Plot 4 - Monte-Carlo received power, matched across fidelities.

The fast analytic model (fidelity 1, FAST) and the expensive wave-optics field
solve (fidelity 2) are run on the same downlink. Their received-power
distributions are plotted as outage curves: P(received power < x). The two lie
on top of each other, so the cheap model is validated where it is cheap - and
the agreement holds into the deep-fade tail. Shown at two elevations.
"""
import numpy as np
import matplotlib.pyplot as plt

import presentation.common as C
from olb.links.downlink import downlink_budget

C.use_style()
d = np.load(C.datapath("plotdata.npz"), allow_pickle=True)

scn = C.downlink()


def deterministic_received_baseline(elev):
    """Received power (dBm) with the turbulent coupling term removed."""
    b = downlink_budget(scn, C.geometry(elev), fidelity=0)
    det = sum(float(np.asarray(t.mean_db)) for t in b.terms if t.category != "coupling")
    return float(b.tx_power_dbm) - det


def exceedance(x):
    xs = np.sort(x)
    p = 1.0 - (np.arange(len(xs)) + 0.5) / len(xs)   # P(X > xs)
    return xs, p


fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.6), sharey=True)
elevs = [90, 30]

for ax, elev in zip(axes, elevs):
    base = deterministic_received_baseline(elev)
    c0 = float(d[f"fid0_coup_{elev}"])
    rx1 = base - d[f"fid1_coup_{elev}"]      # received dBm, fidelity 1
    rx2 = base - d[f"fid2_coup_{elev}"]      # received dBm, fidelity 2
    p0 = base - c0                            # fidelity-0 single value

    # Outage curve: P(received power < x)
    for rx, color, lab in [(rx1, C.FID_COLORS[1], "Fidelity 1 · FAST"),
                           (rx2, C.FID_COLORS[2], "Fidelity 2 · wave optics")]:
        xs = np.sort(rx)
        cdf = (np.arange(len(xs)) + 0.5) / len(xs)   # P(X <= x)
        ax.semilogy(xs, cdf, color=color, lw=2.4, label=lab, alpha=0.9)

    ax.axvline(p0, color=C.FID_COLORS[0], lw=2, ls="--",
               label="Fidelity 0 · analytic mean")

    # 90% availability agreement annotation (10% outage level)
    q90_1, q90_2 = np.percentile(rx1, 10), np.percentile(rx2, 10)
    ax.text(0.03, 0.03,
            f"90% availability level\nfid1 {q90_1:5.1f} dBm\nfid2 {q90_2:5.1f} dBm",
            transform=ax.transAxes, va="bottom", ha="left", fontsize=9.5,
            bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="#cccccc"))

    ax.set_title(f"{elev}° elevation", loc="left")
    ax.set_xlabel("Received power in fibre  [dBm]")
    ax.set_ylim(1e-2, 1)
    ax.grid(True, which="both", alpha=0.5)

axes[0].set_ylabel("Outage probability   P(received < x)")
axes[0].legend(loc="upper left", fontsize=9.5)

fig.suptitle("Monte-Carlo received power agrees across fidelities — the fast model "
             "is validated by the field solve", x=0.5)
fig.tight_layout(rect=[0, 0, 1, 0.95])
out = C.figpath("4_mc_matching.png")
fig.savefig(out)
print("wrote", out)
