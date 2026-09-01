'''
The fidelity-2 VACUUM far-field spot at the satellite, for each obscuration
radius of the uplink sweep.

This is the deterministic (no-turbulence) far-field of the annular launch pupil,
from olb.waveoptics.run.propagate_scenario (the "at rx plane" field over the full
slant range on the co-moving grid). It is the picture behind the mean-loss curve
of uplink_obscuration_dios_vs_waveoptics.py: as the central obscuration grows past
the waist, the Gaussian core is blocked, the on-axis peak collapses (the launch
truncation), and the surviving ring paints a broad Airy-like pattern that the
Dios waist model does not carry.

Each panel is normalised to its OWN peak (so the SHAPE is visible at every
obscuration); the on-axis loss relative to the unobscured spot is in the title,
because the absolute peak falls by tens of dB across the sweep.

Run from the repo root:
    python -m validation.uplink_sigma2i.uplink_obscuration_farfield
'''

import os
import warnings

import numpy as np

from olb.geometry import CircularOrbit
from olb.scenario import Site, SpaceScenario, Channel
from olb.terminal import Terminal, Transmitter, Aperture
from olb.waveoptics.field import Intensity
from olb.waveoptics.run import propagate_scenario

warnings.simplefilter("ignore")

LAM = 1550e-9
W0 = 0.06
APERTURE_M = 0.40                 # launch DIAMETER; radius R = 0.20 m
ALT_M = 600e3
ELEV_DEG = 60.0
CN2_GROUND = 1.7e-14
SAT_APERTURE_M = 0.05
EPS_LIST = [0.0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]
HALF_WINDOW_M = 12.0              # crop half-side of the far-field picture [m]
LOG_FLOOR_DB = 30.0               # show intensity down to this many dB below peak


def _rx_field(eps):
    '''The vacuum rx-plane (satellite) field for an annular pupil of ratio eps.'''
    ground = Terminal(aperture_m=APERTURE_M, wavelength_m=LAM,
                      obscuration_ratio=eps,
                      transmitter=Transmitter(waist_m=W0, power_dbm=42.0))
    space = Terminal(aperture_m=SAT_APERTURE_M, wavelength_m=LAM,
                     detector=Aperture(sensitivity_dbm=-40.0))
    scn = SpaceScenario(ground=ground, space=space, direction="uplink",
                        channel=Channel(site=Site(cn2_ground=CN2_GROUND),
                                        altitude_m=ALT_M))
    r = propagate_scenario(scn, CircularOrbit(ALT_M, elevation_deg=ELEV_DEG))
    _, F = r.stages[2]            # "at rx plane"
    return Intensity(F), float(F.siz)


def _crop(I, siz, half_window_m):
    '''Crop the central +/- half_window_m of a centred field. Returns (I, extent_m).'''
    n = I.shape[0]
    dx = siz / n
    half_px = int(round(half_window_m / dx))
    c = n // 2
    lo, hi = max(0, c - half_px), min(n, c + half_px)
    sub = I[lo:hi, lo:hi]
    ext = (hi - c) * dx
    return sub, ext


def main():
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib is not available; cannot draw the far-field figure.")
        return

    R = APERTURE_M / 2.0
    fig, axes = plt.subplots(2, 4, figsize=(15, 8))
    peak0 = None
    for ax, eps in zip(axes.ravel(), EPS_LIST):
        I, siz = _rx_field(eps)
        peak = float(I.max())
        if peak0 is None:
            peak0 = peak
        onaxis_loss_db = 10.0 * np.log10(peak0 / peak) if peak > 0 else np.inf
        sub, ext = _crop(I, siz, HALF_WINDOW_M)
        # Log scale, normalised to this panel's own peak, floored for display.
        log_norm = 10.0 * np.log10(np.maximum(sub / sub.max(), 1e-12))
        log_norm = np.clip(log_norm, -LOG_FLOOR_DB, 0.0)
        im = ax.imshow(log_norm, extent=[-ext, ext, -ext, ext], origin="lower",
                       cmap="inferno", vmin=-LOG_FLOOR_DB, vmax=0.0)
        eR = eps * R / W0
        ax.set_title(f"eps={eps:.2f}  (eps R/w0={eR:.2f})\n"
                     f"on-axis loss vs eps=0: {onaxis_loss_db:.1f} dB", fontsize=10)
        ax.set_xlabel("x at satellite [m]", fontsize=8)
        ax.set_ylabel("y [m]", fontsize=8)
        ax.tick_params(labelsize=7)

    cbar = fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.6,
                        label="intensity [dB below panel peak]")
    fig.suptitle(
        "Fidelity-2 vacuum far-field spot at the satellite vs launch obscuration\n"
        f"{ELEV_DEG:.0f} deg, {ALT_M / 1e3:.0f} km, w0={W0} m, R={R} m, "
        f"lambda={LAM * 1e9:.0f} nm  (each panel normalised to its own peak)",
        fontsize=12)
    here = os.path.dirname(os.path.abspath(__file__))
    figures = os.path.join(here, "figures")
    os.makedirs(figures, exist_ok=True)
    png = os.path.join(figures, "uplink_obscuration_farfield.png")
    fig.savefig(png, dpi=120, bbox_inches="tight")
    print(f"figure -> {png}")


if __name__ == "__main__":
    main()
