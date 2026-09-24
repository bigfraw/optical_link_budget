# ROUGH DIAGNOSTIC (2026-09-23), copied from a session scratchpad. See ../HANDOVER.md.
# Run from the repository root: python -m validation.divergence_sampling.diagnostics.vacuum_flatness
import warnings; warnings.simplefilter("ignore")
import dataclasses, numpy as np
from olb.scenario import SpaceScenario, Channel
from olb.geometry import CircularOrbit
from olb.terminal import Terminal, Transmitter, Aperture
from olb.models.waveoptics import run_fidelity2
from olb.waveoptics.turbulence.run import _ground_transmit_mode, space_vacuum_baseline
from olb.waveoptics.turbulence.splitstep import super_gaussian_boundary
from olb.waveoptics.turbulence.sampling import PRESETS
lam=1550e-9; geo = CircularOrbit(600e3, 30.0)
def g(div): return Terminal(aperture_m=0.3, wavelength_m=lam, transmitter=Transmitter(waist_m=0.1, divergence_rad=div))
s = SpaceScenario(ground=g(None), space=Terminal(aperture_m=0.05, wavelength_m=lam, detector=Aperture()), direction="uplink", channel=Channel(altitude_m=600e3))
w0 = run_fidelity2(s, geo, n_trials=1, seed=1, progress=False, preset="rapid").turbulent
for m in (1, 2, 4):
    gm = dataclasses.replace(w0.grid, n=w0.grid.n*m)
    mask = super_gaussian_boundary(gm.n, PRESETS["rapid"].boundary_width_frac)
    Fv, _ = space_vacuum_baseline(gm, w0.plan, lam, mask, np.complex128)
    x = (np.arange(gm.n)-gm.n//2)*gm.size_m/gm.n; X,Y = np.meshgrid(x,x); ap = np.hypot(X,Y) <= 0.15
    U = Fv.field
    msg = f"n {gm.n}: F_vac in aperture |U| {np.abs(U[ap]).min():.3f}..{np.abs(U[ap]).max():.3f}, phase p-p {np.ptp(np.angle(U[ap]*np.conj(U[gm.n//2,gm.n//2]))):.3f} rad | baseline/ideal:"
    for div in (None, 50e-6, 100e-6):
        psi = _ground_transmit_mode(g(div), gm)
        msg += f" {div}: {abs((U*psi).sum())**2 / (abs(psi.sum())**2*abs(U[gm.n//2,gm.n//2])**2):.3f}"
    print(msg)
