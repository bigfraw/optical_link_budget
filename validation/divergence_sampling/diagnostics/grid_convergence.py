# ROUGH DIAGNOSTIC (2026-09-23), copied from a session scratchpad. See ../HANDOVER.md.
# Run from the repository root: python -m validation.divergence_sampling.diagnostics.grid_convergence
import warnings; warnings.simplefilter("ignore")
import dataclasses, numpy as np
from olb.scenario import SpaceScenario, Channel
from olb.geometry import CircularOrbit
from olb.terminal import Terminal, Transmitter, Aperture
from olb.models.waveoptics import run_fidelity2
from olb.waveoptics.turbulence.run import propagate_turbulent_scenario as pts
lam=1550e-9; geo = CircularOrbit(600e3, 30.0)
def s(div):
    return SpaceScenario(ground=Terminal(aperture_m=0.3, wavelength_m=lam, transmitter=Transmitter(waist_m=0.1, divergence_rad=div)),
        space=Terminal(aperture_m=0.05, wavelength_m=lam, detector=Aperture()), direction="uplink", channel=Channel(altitude_m=600e3))
_w0 = run_fidelity2(s(None), geo, n_trials=1, seed=1, progress=False, preset="rapid").turbulent; g0, plan = _w0.grid, _w0.plan
for div in (None, 50e-6, 100e-6):
    for m in (1, 2, 4):     # 4 (1024 px) takes about 2.5 min per divergence locally
        g = dataclasses.replace(g0, n=g0.n*m)
        w = pts(s(div), geo, n_trials=64, seed=7, preset="rapid", grid=g, plan=plan)
        e = np.array([t.eta_turb for t in w.trials])
        print(f"div {div} n {g.n}: mean {e.mean():.3f} +/- {e.std()/8:.3f}, frac>1 {np.mean(e>1):.2f}")
