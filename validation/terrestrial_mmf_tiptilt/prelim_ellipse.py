"""Prelim vacuum picture: a 1 inch aperture seen at 45 deg, on transmit and receive.

The aperture projects to an ellipse of semi-axes 12.7 mm (y) x 12.7 cos45 mm
(x). The launch is a Gaussian of 1/e^2 diameter 4 x 4.65 mm (the collimator
beam through the 4x telescope), clipped by the same aperture as the receiver
(circle with circle, ellipse with ellipse). The field goes to the far end by
the single-FFT Fresnel integral (Goodman, ISBN 978-0974707723, Eq. 4-17),
valid in the near field too (at 1 km, L/z_R = 1.4). The receive pupil is the
propagated field clipped by the aperture and focused at f = 4 x 17.5 mm onto
a 105 um core (Fraunhofer, same source). The fibre NA gate rho <= f NA
(Snyder and Love, DOI 10.1007/978-1-4613-2813-1) is 15.4 mm, outside the
aperture, so it is open. Vacuum only, no turbulence.

    python prelim_ellipse.py <out.png> <L_km>
"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, Ellipse
from scipy.interpolate import RegularGridInterpolator

LAM = 1550e-9
L = float(sys.argv[2]) * 1e3
W0 = 4 * 4.65e-3 / 2               # 1/e^2 radius: 4 x the 4.65 mm diameter, halved
A = 25.4e-3 / 2                      # semi-major (y, the tilt axis)
B = A * np.cos(np.radians(45))       # semi-minor (x, the tilt plane)
F = 4 * 17.5e-3
CORE, NA = 105e-6 / 2, 0.22
K = 2 * np.pi / LAM
ZR = np.pi * W0 ** 2 / LAM
WL = W0 * np.sqrt(1 + (L / ZR) ** 2)  # unclipped beam radius at L
out = sys.argv[1]


def axis(n, side):
    return (np.arange(n) - n // 2) * side / n


def masks(X, Y):
    return {"circle": X ** 2 + Y ** 2 <= A ** 2,
            "ellipse": (X / B) ** 2 + (Y / A) ** 2 <= 1}


def fft2c(a):
    return np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(a), norm="ortho"))


# --- transmit and propagate -------------------------------------------------
# a = sqrt(power per pixel). Single-FFT Fresnel: dx2 = lam L / (N dx1). The
# input chirp is sampled finely because the field is zero past 12.7 mm.
n, side = 4096, 1.0
x1 = axis(n, side)
dx1 = side / n
X1, Y1 = np.meshgrid(x1, x1)
g = np.exp(-(X1 ** 2 + Y1 ** 2) / W0 ** 2) * dx1
P_pre = np.sum(np.abs(g) ** 2)          # the power out of the telescope
dx2 = LAM * L / side
x2 = axis(n, n * dx2)
X2, Y2 = np.meshgrid(x2, x2)
tx = {"unclipped": np.ones((n, n), bool), **masks(X1, Y1)}
far = {}
for k, m in tx.items():
    a1 = g * m * np.exp(1j * K * (X1 ** 2 + Y1 ** 2) / (2 * L))
    a2 = fft2c(a1) * np.exp(1j * K * (X2 ** 2 + Y2 ** 2) / (2 * L))
    far[k] = dict(launch=np.sum(np.abs(g * m) ** 2) / P_pre, E=a2 / dx2)

# --- receive ----------------------------------------------------------------
nr, sr = 4096, 0.4                       # focal pixel lam f / sr = 0.27 um
xr = axis(nr, sr)
Xr, Yr = np.meshgrid(xr, xr)
rx_m = masks(Xr, Yr)
box = np.abs(xr) <= A + 1e-3
xb = xr[box]
Xb, Yb = np.meshgrid(xb, xb)
xf = (np.arange(nr) - nr // 2) * LAM * F / sr
Xf, Yf = np.meshgrid(xf, xf)
in_core = Xf ** 2 + Yf ** 2 <= CORE ** 2
j = np.abs(x2) <= 3 * A + 3 * dx2       # interpolate from a small window only
res = {}
for k in ("circle", "ellipse"):
    E = far[k]["E"][np.ix_(j, j)]
    interp = [RegularGridInterpolator((x2[j], x2[j]), part) for part in (E.real, E.imag)]
    pts = np.stack([Yb.ravel(), Xb.ravel()], -1)
    Eb = (interp[0](pts) + 1j * interp[1](pts)).reshape(Xb.shape)
    pup = np.zeros((nr, nr), complex)
    pup[np.ix_(box, box)] = Eb
    pup *= rx_m[k]
    P_ap = np.sum(np.abs(pup) ** 2) * (sr / nr) ** 2
    I = np.abs(fft2c(pup)) ** 2
    res[k] = dict(pup=pup, I=I, launch=far[k]["launch"], collect=P_ap / P_pre,
                  eta=I[in_core].sum() / I.sum())
    res[k]["total_db"] = -10 * np.log10(res[k]["collect"] * res[k]["eta"])


def db(x):
    return -10 * np.log10(x)


# --- figure -----------------------------------------------------------------
mm = 1e3
fig, ax = plt.subplots(2, 3, figsize=(15, 9.8))

a = ax[0, 0]
s = np.abs(x1) <= 20e-3
e = np.array([x1[s][0], x1[s][-1]] * 2) * mm
a.imshow((np.abs(g) ** 2 * tx["ellipse"])[np.ix_(s, s)], extent=e, origin="lower", cmap="inferno")
a.add_patch(Ellipse((0, 0), 2 * B * mm, 2 * A * mm, fill=False, ec="c", lw=1.5, label="45 deg ellipse"))
a.add_patch(Circle((0, 0), A * mm, fill=False, ec="w", ls="--", lw=1, label="1 in circle"))
a.add_patch(Circle((0, 0), W0 * mm, fill=False, ec="lime", ls=":", lw=1, label=f"w = {W0*mm:.1f} mm"))
a.set(title=f"TX pupil (ellipse clip shown)\nlaunched: circle {db(res['circle']['launch']):.2f} dB, "
            f"ellipse {db(res['ellipse']['launch']):.2f} dB clip loss",
      xlabel="x [mm] (tilt plane)", ylabel="y [mm]")
a.legend(loc="lower right", fontsize=7)

span = 3 * WL
s2 = np.abs(x2) <= span
e2 = np.array([x2[s2][0], x2[s2][-1]] * 2)
Iu0 = np.abs(far["unclipped"]["E"][n // 2, n // 2]) ** 2
a = ax[0, 1]
a.imshow(10 * np.log10(np.abs(far["ellipse"]["E"][np.ix_(s2, s2)]) ** 2 / Iu0 + 1e-12),
         extent=e2, origin="lower", cmap="viridis", vmin=-30, vmax=0)
a.set(title=f"Irradiance at {L/1e3:g} km, ellipse launch\n[dB rel. unclipped on-axis]",
      xlabel="x [m]", ylabel="y [m]")

a = ax[0, 2]
c = n // 2
for k, st in (("unclipped", "k-"), ("circle", "b--"), ("ellipse", "r-")):
    I = np.abs(far[k]["E"]) ** 2 / Iu0
    a.plot(x2[s2], 10 * np.log10(I[c, s2] + 1e-12), st, lw=1, label=f"{k}, x cut")
    if k == "ellipse":
        a.plot(x2[s2], 10 * np.log10(I[s2, c] + 1e-12), "m:", lw=1.5, label="ellipse, y cut")
a.set(ylim=(-40, 3), title=f"Cuts at {L/1e3:g} km (vacuum)\n"
      f"on-axis: circle {db(np.abs(far['circle']['E'][c, c])**2/Iu0):+.2f} dB, "
      f"ellipse {db(np.abs(far['ellipse']['E'][c, c])**2/Iu0):+.2f} dB loss",
      xlabel="position [m]", ylabel="dB rel. unclipped on-axis")
a.legend(fontsize=7)
a.grid(alpha=0.3)

a = ax[1, 0]
sb = np.abs(xr) <= 17e-3
eb = np.array([xr[sb][0], xr[sb][-1]] * 2) * mm
Ip = np.abs(res["ellipse"]["pup"][np.ix_(sb, sb)]) ** 2
a.imshow(np.where(Ip > 0, Ip, np.nan), extent=eb, origin="lower", cmap="inferno")
a.add_patch(Circle((0, 0), A * mm, fill=False, ec="k", ls="--", label="1 in circle"))
a.add_patch(Circle((0, 0), F * NA * mm, fill=False, ec="r", lw=1.5, label=f"NA gate f*NA = {F*NA*mm:.1f} mm"))
a.set(xlim=(-17, 17), ylim=(-17, 17), xlabel="x [mm]", ylabel="y [mm]",
      title=f"RX pupil (ellipse), f = {F*mm:.0f} mm\ncollected rel. pre-clip: circle "
            f"{db(res['circle']['collect']):.2f} dB, ellipse {db(res['ellipse']['collect']):.2f} dB")
a.legend(fontsize=7, loc="lower right")

h = np.abs(xf) <= 70e-6
ef = np.array([xf[h][0], xf[h][-1]] * 2) * 1e6
for i, k in enumerate(("circle", "ellipse")):
    a = ax[1, 1 + i]
    I = res[k]["I"][np.ix_(h, h)]
    a.imshow(10 * np.log10(I / I.max() + 1e-12), extent=ef, origin="lower",
             cmap="inferno", vmin=-40, vmax=0)
    a.add_patch(Circle((0, 0), CORE * 1e6, fill=False, ec="c", lw=1.5))
    a.set(title=f"Fibre tip, {k} both ends (dB)\nin core {res[k]['eta']:.4f}; "
                f"telescope-to-fibre {res[k]['total_db']:.2f} dB",
          xlabel="x [um]", ylabel="y [um]")

fig.suptitle(f"Prelim, VACUUM ONLY, {L/1e3:g} km: 1 in aperture at 45 deg (25.4 x 17.96 mm ellipse) on TX and RX; "
             f"2w = 4 x 4.65 mm, f = 4 x 17.5 mm, 105 um / NA 0.22, 1550 nm "
             f"(unclipped w(L) = {WL*mm:.0f} mm)", fontsize=10)
fig.tight_layout()
fig.savefig(out, dpi=110)
for k in ("circle", "ellipse"):
    r = res[k]
    print(L / 1e3, k, f"clip {db(r['launch']):.3f}  collect {db(r['collect']):.3f}  "
          f"eta {r['eta']:.4f}  total {r['total_db']:.3f} dB")
