"""Plot 6 - uplink vs downlink asymmetry against elevation.

The two directions are not symmetric. The downlink is a near-plane wave that a
big telescope averages down to almost nothing. The uplink launches from the
ground, and beam wander makes it fluctuate hard at a point receiver on the
satellite that cannot aperture-average. This plot sweeps elevation and shows
both the fading strength and the beam wander.
"""
import numpy as np
import matplotlib.pyplot as plt

import presentation.common as C
from olb.geometry import CircularOrbit
from olb.units import w0_to_div
from olb.turbulence.plane_wave_scintillation import (
    plane_wave_scintillation_index, aperture_averaged_scintillation_index)
from olb.turbulence.andrews.beam import beam_params
from olb.turbulence.andrews.paths import uplink_scintillation_index
from olb.turbulence.andrews.wander import beam_wander_variance_slant

C.use_style()

elev = np.linspace(10, 90, 81)
W0_UP = 0.15                      # ground uplink launch beam radius [m], collimated
D_RX = C.GROUND_APERTURE_M        # ground telescope diameter
ranges = np.array([CircularOrbit(C.ALTITUDE_M, e).slant_range_m for e in elev]).ravel()

# --- Downlink: point and aperture-averaged scintillation index (vectorised) --
s2_down_pt = np.asarray(plane_wave_scintillation_index(
    elev, C.WAVELENGTH_M, C.HS, C.CN2)).ravel()
s2_down_ap = np.asarray(aperture_averaged_scintillation_index(
    D_RX, elev, C.WAVELENGTH_M, C.HS, C.CN2)).ravel()

# --- Uplink: untracked beam-wave index (includes beam wander) ---------------
s2_up = np.array([
    uplink_scintillation_index(C.HS, C.CN2, C.WAVELENGTH_M, float(e),
                               beam_params(W0_UP, C.WAVELENGTH_M, z=float(L), f0=np.inf),
                               altitude_m=C.ALTITUDE_M, r=0.0, tracked=False, regime="auto")
    for e, L in zip(elev, ranges)])

# --- Uplink beam wander -> pointing error angle -----------------------------
rc2 = np.array([
    beam_wander_variance_slant(W0_UP, C.WAVELENGTH_M, C.HS, C.CN2, float(L),
                               f0=np.inf, elevation_deg=float(e))
    for e, L in zip(elev, ranges)])
wander_urad = np.sqrt(rc2) / ranges * 1e6      # radial angle [µrad]
div_urad = w0_to_div(W0_UP, C.WAVELENGTH_M) * 1e6

fig, (axL, axR) = plt.subplots(1, 2, figsize=(12.8, 5.6))

# Panel A: scintillation index
axL.semilogy(elev, s2_up, color=C.ORANGE, lw=2.8, label="Uplink — untracked (beam wander)")
axL.semilogy(elev, s2_down_pt, color=C.BLUE, lw=2.6, label="Downlink — point receiver")
axL.semilogy(elev, s2_down_ap, color=C.BLUE, lw=2.4, ls="--",
             label=f"Downlink — {D_RX:.1f} m aperture-averaged")
axL.axhline(1.0, color=C.GREY, lw=1.2, ls=":")
axL.text(11, 1.15, "weak / strong  ($\\sigma_I^2$ = 1)", color=C.GREY, ha="left",
         va="bottom", fontsize=9)
axL.set_xlabel("Elevation angle  [deg]")
axL.set_ylabel("Scintillation index  $\\sigma_I^2$")
axL.set_title("Fading strength", loc="left")
axL.legend(loc="lower left", fontsize=9.5)
axL.set_xlim(10, 90)
axL.set_ylim(5e-4, 2)

# Panel B: uplink beam wander
axR.plot(elev, wander_urad, color=C.ORANGE, lw=2.8, label="Uplink beam wander")
axR.axhline(div_urad, color=C.GREEN, lw=1.6, ls=":",
            label=f"launch divergence ({div_urad:.1f} µrad)")
axR.plot(elev, np.zeros_like(elev), color=C.BLUE, lw=2.6, label="Downlink (no beam wander)")
axR.set_xlabel("Elevation angle  [deg]")
axR.set_ylabel("Pointing error from wander  [µrad]")
axR.set_title("Beam wander", loc="left")
axR.legend(loc="upper right", fontsize=9.5)
axR.set_xlim(10, 90)
axR.set_ylim(bottom=-0.5)

fig.suptitle("Uplink vs downlink — the atmosphere hits the two directions differently",
             x=0.5)
fig.text(0.5, 0.005,
         f"500 km LEO · 1550 nm · uplink launch radius {W0_UP*100:.0f} cm · "
         f"HV5/7 turbulence",
         ha="center", color=C.GREY, fontsize=9)
fig.tight_layout(rect=[0, 0.02, 1, 0.96])
out = C.figpath("6_asymmetry.png")
fig.savefig(out)
print("wrote", out)
