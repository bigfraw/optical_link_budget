"""Show the grid-sampling limit of a DELIBERATELY DIVERGED fidelity-2 uplink.

THE CLAIM UNDER TEST. A space uplink at fidelity 2 reads its turbulence fade
through the Shapiro reciprocity overlap (DOI 10.1364/JOSA.61.000492):

    eta_turb = |sum(F_turb psi_tx)|^2 / |sum(F_vac psi_tx)|^2   (NO conjugate)

`psi_tx` is the ground transmit mode on the turbulent grid. A diverged beam
has a parabolic phase of radius R in the aperture plane. Its local spatial
frequency is x / (lambda R) (Schmidt (2010), DOI 10.1117/3.866274, Ch. 7,
Eq. (7.37), printed p. 122), so the local wavefront tilt at the edge of a
filled aperture is about the divergence theta. The grid holds tilts up to

    theta_max = lambda / (2 dx)

(Schmidt, Ch. 7, Eq. (7.7), printed p. 117). The space sizer does NOT read the
launch curvature. So past theta_max the transmit mode ALIASES, and the ratio
above is wrong.

THE GHOST LATTICE. On a grid of pitch dx, the phase pi x^2 / (lambda R) at
x + X, with X = lambda R / dx, differs from the phase at x by 2 pi n (an
integer number of turns at every pixel) plus a constant. So the sampled
phase is the true bowl plus copies of it on a square lattice of spacing X,
each with its own constant phase. This follows from Eq. (7.37) and the
Nyquist rule of Eq. (7.7).

THE METHOD. One scenario: a 0.3 m ground aperture, a 0.1 m waist, 1550 nm, a
600 km orbit at 30 deg, the `rapid` preset. The script takes the PRODUCTION
turbulent grid from `run_fidelity2`, and it builds the transmit mode with the
production `_ground_transmit_mode` on that grid and on a grid 8 times finer
(the same side). The fine grid is the TRUE reference: it holds the mode to
about 900 urad.

THE OUTPUTS (in figures/):
  1_phase_cut.png       the phase and the local tilt along one cut, for
                        50 / 100 / 300 urad.
  2_phase_map.png       the wrapped phase inside the aperture, and the far
                        field, at 300 urad.
  3_overlap_error.png   the vacuum overlap (the denominator) on the grid
                        against the true one, over 5 to 400 urad.
The console prints the turbulent smoke run: eta_turb for 3 divergences.

Run from the repository root:

    python -m validation.divergence_sampling.divergence_sampling
"""
import dataclasses
import pathlib
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from olb.geometry import CircularOrbit
from olb.models.waveoptics import run_fidelity2
from olb.scenario import Channel, SpaceScenario
from olb.terminal import Aperture, Terminal, Transmitter
from olb.waveoptics.schmidt.sampling import nyquist_max_angle
from olb.waveoptics.turbulence.run import _ground_transmit_mode

warnings.simplefilter("ignore")

FIG = pathlib.Path(__file__).parent / "figures"
LAM = 1550e-9
APERTURE = 0.3
WAIST = 0.1
FINE = 8
GEO = CircularOrbit(600e3, 30.0)


def ground(div):
    return Terminal(aperture_m=APERTURE, wavelength_m=LAM,
                    transmitter=Transmitter(waist_m=WAIST, power_dbm=40,
                                            divergence_rad=div))


def scenario(div):
    return SpaceScenario(
        ground=ground(div),
        space=Terminal(aperture_m=0.05, wavelength_m=LAM, detector=Aperture()),
        direction="uplink", channel=Channel(altitude_m=600e3))


def mode(grid, div, m):
    """The production transmit mode on the grid refined m times (same side)."""
    gm = dataclasses.replace(grid, n=grid.n * m)
    return _ground_transmit_mode(ground(div), gm), gm.size_m / gm.n, gm.n


def axis(n, d):
    return (np.arange(n) - n // 2) * d


def vacuum_overlap(grid, div, m):
    """|integral psi dA|^2 of the unit-power mode against a unit plane wave, m^2."""
    psi, d, _ = mode(grid, div, m)
    return abs(psi.sum()) ** 2 * d ** 2


def fig_phase_cut(grid, dx, th_nyq):
    fig, ax = plt.subplots(2, 3, figsize=(15, 8), sharex=True)
    for j, div in enumerate((50e-6, 100e-6, 300e-6)):
        for m, st, lab in ((FINE, "-", f"true ({FINE}x finer grid)"),
                           (1, "o-", f"production grid, dx = {dx*1e3:.2f} mm")):
            psi, d, n = mode(grid, div, m)
            x = axis(n, d)
            k = np.abs(x) <= APERTURE / 2
            ph = np.unwrap(np.angle(psi[n // 2][k]))
            ph -= ph[len(ph) // 2]
            ax[0, j].plot(x[k] * 100, ph, st, ms=3, lw=1, label=lab)
            tilt = np.gradient(ph, d) * LAM / (2 * np.pi) * 1e6
            ax[1, j].plot(x[k] * 100, tilt, st, ms=3, lw=1, label=lab)
        for yy in (th_nyq, -th_nyq):
            ax[1, j].axhline(yy * 1e6, color="r", ls="--", lw=1)
        ax[0, j].set_title(f"divergence {div*1e6:.0f} urad  "
                           f"(limit {th_nyq*1e6:.0f} urad)")
        ax[1, j].set_xlabel("x across the aperture [cm]")
    ax[0, 0].set_ylabel("unwrapped phase [rad]")
    ax[1, 0].set_ylabel("local wavefront tilt [urad]")
    ax[1, 0].text(0.02, 0.9, "red dashes: pi per pixel\n(the most the grid can show)",
                  transform=ax[1, 0].transAxes, color="r", fontsize=8, va="top")
    ax[0, 0].legend(fontsize=8)
    fig.suptitle("1. The grid folds the wavefront slope back once it passes pi per pixel")
    fig.tight_layout()
    fig.savefig(FIG / "1_phase_cut.png", dpi=130)


def fig_phase_map(grid, dx, th_nyq, div=300e-6):
    fig, ax = plt.subplots(2, 2, figsize=(12, 11))
    for i, (m, lab) in enumerate(((FINE, f"true ({FINE}x finer grid)"),
                                  (1, "production grid"))):
        psi, d, n = mode(grid, div, m)
        x = axis(n, d)
        X, Y = np.meshgrid(x, x)
        k = np.hypot(X, Y) <= APERTURE / 2 * 1.02
        c = np.abs(x) <= APERTURE / 2 * 1.02
        ph = np.where(k, np.angle(psi), np.nan)[np.ix_(c, c)]
        ext = [x[c][0] * 100, x[c][-1] * 100] * 2
        im = ax[0, i].imshow(ph, extent=ext, cmap="twilight",
                             interpolation="nearest", origin="lower")
        ax[0, i].set_title(f"wrapped phase, {lab}\ndx = {d*1e3:.2f} mm")
        ax[0, i].set_xlabel("x [cm]")
        # The far field: the zero-padded FFT, on the angle axis lambda fx.
        P = 4
        I = np.abs(np.fft.fftshift(np.fft.fft2(psi, s=(P * n, P * n)))) ** 2
        I /= I.max()
        ang = np.fft.fftshift(np.fft.fftfreq(P * n, d)) * LAM * 1e6
        w = np.abs(ang) <= 1500
        ax[1, i].imshow(10 * np.log10(I[np.ix_(w, w)] + 1e-6),
                        extent=[ang[w][0], ang[w][-1]] * 2,
                        cmap="magma", vmin=-40, vmax=0, origin="lower")
        ax[1, i].add_patch(plt.Circle((0, 0), div * 1e6, fill=False,
                                      color="c", ls="--"))
        if m > 1:
            ax[1, i].set_title("far-field intensity [dB]\n"
                               "(cyan: the real 300 urad divergence)")
        else:
            ax[1, i].set_title(
                "far-field intensity [dB] - this panel is the WHOLE angle range the grid\n"
                f"can hold: +/- lam/(2 dx) = {th_nyq*1e6:.0f} urad. The 300 urad beam\n"
                "folds into it as a lattice of orders.", fontsize=9)
        ax[1, i].set_xlabel("angle [urad]")
    fig.colorbar(im, ax=ax[0, :], label="phase [rad]", shrink=0.8, pad=0.08)
    ghost = LAM * (WAIST / div) / dx      # X = lambda R / dx, R ~ w / theta
    ax[0, 1].text(0.02, 0.02, f"ghost centres every lam R / dx ~ {ghost*100:.1f} cm",
                  transform=ax[0, 1].transAxes, color="w", fontsize=9,
                  bbox=dict(fc="k", alpha=0.6))
    fig.suptitle("2. At 300 urad the grid holds a lattice of false bowl centres, "
                 "and the far field breaks into orders")
    fig.savefig(FIG / "2_phase_map.png", dpi=130, bbox_inches="tight")


def fig_overlap_error(grid, th_nyq):
    divs = np.linspace(5e-6, 400e-6, 60)
    err = [10 * np.log10(vacuum_overlap(grid, dv, 1) / vacuum_overlap(grid, dv, FINE))
           for dv in divs]
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(divs * 1e6, err, "o-", ms=3)
    ax.axvline(th_nyq * 1e6, color="r", ls="--",
               label=f"lam/(2 dx) = {th_nyq*1e6:.0f} urad (Nyquist)")
    ax.axvline(th_nyq * 0.5e6, color="orange", ls=":",
               label=f"lam/(4 dx) = {th_nyq*0.5e6:.0f} urad (margin)")
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xlabel("launch divergence [urad]")
    ax.set_ylabel("vacuum overlap error, grid vs true [dB]")
    ax.set_title("3. The denominator of the turbulence ratio: correct below the "
                 "limit, arbitrary above it\n"
                 "(negative = too small -> spurious GAIN in the turbulence row)")
    ax.text(0.02, 0.05, "Up to ~145 urad the error stays under +/-0.4 dB: only the "
            "far edge of the\naperture aliases, and the Gaussian carries little "
            "power there.", transform=ax.transAxes, fontsize=8)
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIG / "3_overlap_error.png", dpi=130)
    return divs, np.array(err)


def smoke(divs=(None, 100e-6, 300e-6), n_trials=16):
    """The turbulent smoke run: eta_turb on the production grid."""
    print("\nturbulent smoke run (rapid preset, %d trials, seed 1)" % n_trials)
    for div in divs:
        w = run_fidelity2(scenario(div), GEO, n_trials=n_trials, seed=1,
                          progress=False, preset="rapid")
        eta = np.array([t.eta_turb for t in w.turbulent.trials])
        print(f"  divergence {div}: eta_turb mean {eta.mean():.3f}, "
              f"max {eta.max():.2f}, fraction > 1 {np.mean(eta > 1):.2f}")


if __name__ == "__main__":
    FIG.mkdir(exist_ok=True)
    grid = run_fidelity2(scenario(None), GEO, n_trials=1, seed=1, progress=False,
                         preset="rapid").turbulent.grid
    dx = grid.size_m / grid.n
    th_nyq = nyquist_max_angle(LAM, dx)
    print(f"production grid: n = {grid.n}, dx = {dx*1e3:.3f} mm, "
          f"theta_max = lambda/(2 dx) = {th_nyq*1e6:.1f} urad")
    fig_phase_cut(grid, dx, th_nyq)
    fig_phase_map(grid, dx, th_nyq)
    divs, err = fig_overlap_error(grid, th_nyq)
    ok = divs <= th_nyq
    print(f"vacuum overlap error below theta_max: max |err| {np.abs(err[ok]).max():.2f} dB")
    print(f"vacuum overlap error above theta_max: min {err[~ok].min():.1f} dB, "
          f"max {err[~ok].max():.1f} dB")
    # The self-check: the grid is exact well below the limit, and wrong past it.
    assert np.abs(err[divs <= 0.5 * th_nyq]).max() < 0.5
    assert np.abs(err[divs >= 2 * th_nyq]).max() > 5.0
    smoke()
