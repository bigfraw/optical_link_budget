"""Gate 0: make the fidelity-2 overlap converge for a DIVERGED launch.

THE PROBLEM (validation/divergence_sampling/HANDOVER.md, Section 3c). The
space slab starts from a unit plane wave that fills the grid. The absorbing
mask cuts it at every hop, so the vacuum field in the aperture carries
Fresnel rings, and the hop count (so the ring pattern) changes with the pixel
count. A curved psi_tx is a chirp, and so is a ring, so the overlap
|sum(F psi_tx)|^2 (Shapiro, DOI 10.1364/JOSA.61.000492) reads the rings.

THE CANDIDATE. Start the slab from a WIDE GAUSSIAN (waist ws) in place of the
plane wave. A Gaussian has no edge, so it makes no rings, and its vacuum
propagation is known in closed form (GForvard, the ABCD law; Siegman, Lasers,
ISBN 978-0935702118). With ws much larger than the aperture it is locally a
plane wave there.

Two parts:
    vacuum     the masked split-step vacuum overlap against the exact
               overlap (plane: |sum psi|^2; Gaussian: GForvard), against the
               pixel count, for every divergence. No screens.
    turbulent  the mean eta_turb of 64 trials against the pixel count, with
               the start field swapped by a patch of run.Begin.

Run from the repository root:
    python -m validation.uplink_divergence.gate0 vacuum --elev 30
    python -m validation.uplink_divergence.gate0 turbulent --elev 30 --start gauss
"""
import argparse
import contextlib
import dataclasses
import warnings

import numpy as np

import olb.waveoptics.turbulence.run as run_mod
from olb.models.waveoptics import run_fidelity2
from olb.waveoptics.field import Begin
from olb.waveoptics.propagators import GForvard
from olb.waveoptics.sources import GaussBeam
from olb.waveoptics.turbulence.run import (_ground_transmit_mode,
                                           propagate_turbulent_scenario,
                                           reciprocity_overlap,
                                           space_vacuum_baseline)
from olb.waveoptics.turbulence.sampling import PRESETS
from olb.waveoptics.turbulence.splitstep import super_gaussian_boundary

from .config import DIVERGENCES, THETA_MIN, geometry, scenario

warnings.simplefilter("ignore")


def base_grid(elev, preset):
    w = run_fidelity2(scenario(), geometry(elev), n_trials=1, seed=1,
                      progress=False, preset=preset).turbulent
    return w.grid, w.plan


@contextlib.contextmanager
def start_field(kind, ws_frac):
    """Swap the plane-wave start of the space slab for a Gaussian.

    run.Begin makes the plane start (the trial start, the device source and
    the vacuum baseline). GaussBeam overwrites the field, so the transmit mode
    that also calls Begin does not change.
    """
    if kind == "plane":
        yield
        return
    orig = run_mod.Begin

    def tapered(size, lam, n, dtype=np.complex128):
        F = GaussBeam(orig(size, lam, n, dtype=dtype), ws_frac * size)
        F.field = F.field.astype(dtype)    # GaussBeam gives a real float64
        return F

    run_mod.Begin = tapered
    try:
        yield
    finally:
        run_mod.Begin = orig


def label(div):
    return "coll" if div is None else f"{div / THETA_MIN:.0f}x"


def vacuum(elev, preset, mults, ws_fracs):
    g0, plan = base_grid(elev, preset)
    lam = scenario().ground.wavelength_m
    print(f"elev {elev} {preset}: side {g0.size_m:.3f} m, "
          f"z {plan.z_total_m:.0f} m")
    for m in mults:
        g = dataclasses.replace(g0, n=g0.n * m)
        mask = super_gaussian_boundary(g.n, PRESETS[preset].boundary_width_frac)
        psis = [_ground_transmit_mode(scenario(d).ground, g) for d in DIVERGENCES]
        for kind, wf in [("plane", None)] + [("gauss", f) for f in ws_fracs]:
            with start_field(kind, wf):
                F_vac, _ = space_vacuum_baseline(g, plan, lam, mask,
                                                 np.complex128)
            if kind == "plane":
                truth = np.ones((g.n, g.n))
            else:
                truth = GForvard(GaussBeam(Begin(g.size_m, lam, g.n),
                                           wf * g.size_m), plan.z_total_m).field
            # One global phase and scale: match the two on axis.
            c = g.n // 2
            truth = truth * F_vac.field[c, c] / truth[c, c]
            row = []
            for d, psi in zip(DIVERGENCES, psis):
                r = (reciprocity_overlap(F_vac.field, psi)
                     / reciprocity_overlap(truth, psi))
                row.append(f"{label(d)} {10 * np.log10(r):+6.2f}")
            name = kind if wf is None else f"gauss {wf:.2f}"
            print(f"  n {g.n:5d} {name:10s} dB err: " + "  ".join(row))


def turbulent(elev, preset, mults, kind, wf, n_trials, divs):
    g0, plan = base_grid(elev, preset)
    for d in divs:
        for m in mults:
            g = dataclasses.replace(g0, n=g0.n * m)
            with start_field(kind, wf):
                w = propagate_turbulent_scenario(
                    scenario(d), geometry(elev), n_trials=n_trials, seed=7,
                    preset=preset, grid=g, plan=plan)
            e = np.array([t.eta_turb for t in w.trials])
            print(f"  {kind} {wf} {label(d):5s} n {g.n:5d}: mean "
                  f"{e.mean():.3f} +/- {e.std() / np.sqrt(e.size):.3f}, "
                  f"frac>1 {np.mean(e > 1):.2f}", flush=True)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument("part", choices=("vacuum", "turbulent"))
    ap.add_argument("--elev", type=float, default=30.0)
    ap.add_argument("--preset", default="rapid")
    ap.add_argument("--mults", type=int, nargs="+", default=[1, 2, 4])
    ap.add_argument("--ws", type=float, nargs="+", default=[0.15, 0.2, 0.25])
    ap.add_argument("--start", default="gauss")
    ap.add_argument("--trials", type=int, default=64)
    ap.add_argument("--divs", type=int, nargs="+", default=[0, 4, 8, 16],
                    help="multiples of theta_min, 0 = collimated")
    a = ap.parse_args()
    if a.part == "vacuum":
        vacuum(a.elev, a.preset, a.mults, a.ws)
    else:
        divs = [None if k == 0 else k * THETA_MIN for k in a.divs]
        turbulent(a.elev, a.preset, a.mults, a.start, a.ws[0], a.trials, divs)
