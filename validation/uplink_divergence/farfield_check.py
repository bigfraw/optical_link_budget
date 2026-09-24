"""Is the non-monotonic Arm B mean loss against the divergence PHYSICAL?

The hard launch clip (D = 1.36 w) gives the vacuum far field of a diverged
beam a ripple, so the on-axis vacuum value (the eta denominator) sits on a
ripple peak or a dip that moves with the divergence. The long-term turbulent
mean is the vacuum far field blurred by the plane-wave coherence
exp(-3.44 (rho/r0)^(5/3)) (Fried, DOI 10.1364/JOSA.56.001372): on the axis it
is sum_rho autocorr(psi)(rho) MTF(rho), over the vacuum value |sum psi|^2.
No screens here.

Run from the repository root as a module:
    -m validation.uplink_divergence.farfield_check
"""
import dataclasses
import warnings

import numpy as np

from olb.waveoptics.turbulence.run import _ground_transmit_mode

from .arm_b import cell_grid
from .config import THETA_MIN, scenario

warnings.simplefilter("ignore")

if __name__ == '__main__':
    for elev in (30.0, 60.0):
        g, plan = cell_grid(elev)
        g = dataclasses.replace(g, n=g.n * 2, size_m=g.size_m * 2)  # pad room
        r0, dx = plan.r0_total_m, g.pixel_m
        x = (np.arange(g.n) - g.n // 2) * dx
        X, Y = np.meshgrid(x, x)
        mtf = np.fft.ifftshift(np.exp(-3.44 * (np.hypot(X, Y) / r0) ** (5 / 3)))
        rows = []
        for th in np.r_[0, np.arange(12, 221, 4)] * 1e-6:
            psi = _ground_transmit_mode(scenario(th or None).ground, g)
            ac = np.fft.ifft2(np.abs(np.fft.fft2(psi)) ** 2)   # autocorrelation
            rows.append((th * 1e6, float((ac * mtf).sum().real
                                         / abs(psi.sum()) ** 2)))
        print(f"elev {elev:.0f}: r0 {r0:.3f} m, expected mean loss (dB):")
        print("  " + "  ".join(f"{t:.0f}:{-10 * np.log10(e):+.2f}"
                               for t, e in rows))
        for m in (2, 4, 8, 16):
            t = m * THETA_MIN * 1e6
            i = int(np.argmin([abs(r[0] - t) for r in rows]))
            print(f"   ~{m}x ({t:.0f} urad): "
                  f"{-10 * np.log10(rows[i][1]):+.2f} dB")
