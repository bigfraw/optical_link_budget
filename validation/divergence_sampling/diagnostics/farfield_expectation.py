# ROUGH DIAGNOSTIC (2026-09-23), copied from a session scratchpad. See ../HANDOVER.md.
# Run from the repository root: python -m validation.divergence_sampling.diagnostics.farfield_expectation
import warnings; warnings.simplefilter("ignore")
import dataclasses, numpy as np
from olb.scenario import SpaceScenario, Channel
from olb.geometry import CircularOrbit
from olb.terminal import Terminal, Transmitter, Aperture
from olb.models.waveoptics import run_fidelity2
from olb.waveoptics.turbulence.run import _ground_transmit_mode
lam=1550e-9
def g(div, ap=0.3, w=0.1):
    return Terminal(aperture_m=ap, wavelength_m=lam, transmitter=Transmitter(waist_m=w, divergence_rad=div))
s = SpaceScenario(ground=g(None), space=Terminal(aperture_m=0.05, wavelength_m=lam, detector=Aperture()),
                  direction="uplink", channel=Channel(altitude_m=600e3))
w = run_fidelity2(s, CircularOrbit(600e3, 30.0), n_trials=1, seed=1, progress=False, preset="rapid")
grid = w.turbulent.grid; print("r0 total", w.turbulent.plan.r0_total_m)
gm = dataclasses.replace(grid, n=grid.n*8); d = gm.size_m/gm.n
r0 = w.turbulent.plan.r0_total_m
for div, ap in ((None,.3),(50e-6,.3),(100e-6,.3),(100e-6,1.0)):
    psi = _ground_transmit_mode(g(div, ap), gm)
    P = 2; F = np.abs(np.fft.fft2(psi, s=(P*gm.n,)*2))**2
    ang = np.fft.fftfreq(P*gm.n, d)*lam*1e6
    prof = [F[0, np.argmin(abs(ang-a))] for a in (0,10,20,40,80)]
    # expected eta = sum psi psi'* Gamma / |sum psi|^2, Gamma = exp(-3.44 (r/r0)^(5/3)/2)?  use MCF of plane wave: exp(-0.5*6.88 (rho/r0)^(5/3))
    kx = np.fft.fftfreq(P*gm.n, d); KX, KY = np.meshgrid(kx, kx)
    # FT of Gamma via numeric: build Gamma on the padded grid
    x = (np.arange(P*gm.n) - P*gm.n//2)*d; X, Y = np.meshgrid(x, x)
    G = np.exp(-0.5*6.88*(np.hypot(X,Y)/r0)**(5/3))
    Gh = np.abs(np.fft.fft2(np.fft.ifftshift(G))); Gh /= Gh.sum()
    Eeta = (np.fft.ifft2(np.fft.fft2(F)*np.fft.fft2(Gh)).real[0,0]) / F[0,0]
    print(f"div {div} ap {ap}: far field at 0/10/20/40/80 urad rel. to axis:",
          np.round(np.array(prof)/prof[0], 2), f" expected mean eta {Eeta:.3f}")
