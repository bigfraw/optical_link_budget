"""Plot 2 - the link-budget waterfall.

Start at the satellite launch power and step down through every loss term to the
power in the fibre. One glance shows where the decibels go: geometric spreading
and the uncompensated fibre coupling dominate. The receiver sensitivity, the
link margin and the 99% fade level are marked.
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

import presentation.common as C

C.use_style()
d = np.load(C.datapath("plotdata.npz"), allow_pickle=True)

names = list(d["wf_names"])
means = np.asarray(d["wf_means"], float)
cats = list(d["wf_cats"])
tx = float(d["wf_tx_dbm"])
total = float(d["wf_total"])
fade90 = float(d["wf_fade90"])          # total loss at 90% availability
sens = float(d["wf_rx_sens_dbm"])
received = tx - total
fade_rx = tx - fade90                     # received power at the 90% fade

CATCOLOR = {"geometric": "#4c72b0", "atmospheric": "#8c8c8c",
            "pointing": "#c44e52", "coupling": "#dd8452",
            "turbulence": "#55a868", "system": "#937860"}
END = "#2f2f2f"

n = len(names)
# Running power level after each step: level[0]=tx, level[i]=tx-sum(L[:i]).
level = tx - np.concatenate([[0.0], np.cumsum(means)])   # length n+1

fig, ax = plt.subplots(figsize=(11, 6.3))
ybase = min(sens, fade_rx, received) - 6
w = 0.62

labels = ["Satellite\nlaunch"] + [nm.replace(" ", "\n", 1) for nm in names] + ["Power in\nfibre"]
xs = np.arange(n + 2)

# Tx and Rx anchor bars
ax.bar(0, tx - ybase, bottom=ybase, width=w, color=END, zorder=3)
ax.text(0, tx + 0.8, f"{tx:.0f} dBm", ha="center", va="bottom", fontweight="bold")
ax.bar(n + 1, received - ybase, bottom=ybase, width=w, color=END, zorder=3)
ax.text(n + 1, received + 0.8, f"{received:.1f} dBm", ha="center", va="bottom",
        fontweight="bold")

# Loss bars
for i, (nm, L, cat) in enumerate(zip(names, means, cats), start=1):
    top, bot = level[i - 1], level[i]
    ax.bar(i, top - bot, bottom=bot, width=w, color=CATCOLOR.get(cat, "#999999"),
           zorder=3, edgecolor="white", linewidth=0.5)
    ax.text(i, top + 0.5, f"-{L:.1f}", ha="center", va="bottom", fontsize=9,
            color="#333333")

# Connectors at each shared level
for j in range(n + 1):
    ax.plot([j + w / 2, j + 1 - w / 2], [level[j], level[j]],
            color="#bbbbbb", lw=1, ls="--", zorder=2)

# Receiver sensitivity
ax.axhline(sens, color=C.RED, lw=1.6, ls=":")
ax.text(n + 1 + 0.55, sens, f"sensitivity {sens:.0f} dBm", color=C.RED,
        va="center", ha="left", fontsize=9)

# Link margin arrow (mean received down to sensitivity)
margin = received - sens
ax.annotate("", xy=(n + 1, received), xytext=(n + 1, sens),
            arrowprops=dict(arrowstyle="<->", color="#444444", lw=1.3))
ax.text(n + 1 - 0.45, (received + sens) / 2, f"link margin\n{margin:.0f} dB",
        ha="right", va="center", fontsize=9.5, color="#333333")

# 90% fade level on the Rx bar
ax.plot([n + 1 - w / 2, n + 1 + w / 2], [fade_rx, fade_rx], color="#7a3fbf", lw=2.2)
ax.text(n + 1, fade_rx - 0.8, f"90% fade\n{fade_rx:.0f} dBm", color="#7a3fbf",
        ha="center", va="top", fontsize=8.5)

ax.set_xticks(xs)
ax.set_xticklabels(labels, fontsize=9)
ax.set_ylabel("Power  [dBm]")
ax.set_title("Downlink budget waterfall — 0.7 m telescope to single-mode fibre, "
             f"{C.HERO_ELEVATION_DEG:.0f}° elevation", loc="left")
ax.set_ylim(ybase, tx + 4)
ax.margins(x=0.02)

legend = [Patch(facecolor=CATCOLOR[c], label=c) for c in
          ["geometric", "atmospheric", "pointing", "coupling"] if c in cats]
ax.legend(handles=legend, loc="upper right", ncol=2, title="loss category")

fig.tight_layout()
out = C.figpath("2_waterfall.png")
fig.savefig(out)
print("wrote", out)
