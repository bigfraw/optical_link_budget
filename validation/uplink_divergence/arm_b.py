"""Arm B: the fidelity-2 power distribution of a diverged uplink.

ONE CAMPAIGN FOR EACH ELEVATION. The downlink slab field does not depend on
the launch, so the campaign stores the receive-field patch ONE time, and the
reader re-reads it for every divergence post hoc: a new ground transmit mode
psi_tx and a new vacuum baseline, then the reciprocity overlap
eta = |sum(F psi_tx)|^2 / o_vac (Shapiro, DOI 10.1364/JOSA.61.000492). So the
divergences share their atmospheres, trial for trial.

THE START FIELD IS THE GATE-0 GAUSSIAN (waist WS_FRAC x the grid side; Arm A picked 0.3), not
the plane wave of record: the plane wave puts grid-edge Fresnel rings in the
aperture, and a curved psi_tx reads them. The campaign fingerprint does NOT
record this start, so the campaigns live in their OWN root, and cell.json
records the start.

THE PIXEL. A curved psi_tx aliases past lambda / (2 dx); the rule here is
dx <= lambda / (4 theta_max) (validation/divergence_sampling/).

    python -m validation.uplink_divergence.arm_b run --elev 30   (bigfraw, cupy)
    python -m validation.uplink_divergence.arm_b read            (local)
"""
import argparse
import dataclasses
import json
import warnings
from pathlib import Path

import numpy as np

from olb.waveoptics.turbulence.campaign import Campaign
from olb.waveoptics.turbulence.run import (_transmit_mode_crop,
                                           _ground_transmit_mode,
                                           reciprocity_overlap,
                                           space_vacuum_baseline)
from olb.waveoptics.turbulence.sampling import PRESETS
from olb.waveoptics.turbulence.splitstep import super_gaussian_boundary

from .config import DIVERGENCES, ELEVATIONS, LAM, geometry, scenario
from .gate0 import base_grid, label, start_field

warnings.simplefilter("ignore")
HERE = Path(__file__).parent
ROOT = HERE / "campaigns"
PRESET = "standard"
WS_FRAC = 0.3
SEED = 20260923
N_TRIALS = 2000


def cell_grid(elev):
    """Give the grid and plan: the sizer side, the pixel <= lambda/(4 theta)."""
    g, plan = base_grid(elev, PRESET)
    dx_max = LAM / (4 * max(d for d in DIVERGENCES if d is not None))
    m = 1
    while g.size_m / (g.n * m) > dx_max:
        m *= 2
    return dataclasses.replace(g, n=g.n * m), plan


def campaign(elev, backend="numpy"):
    g, plan = cell_grid(elev)
    return Campaign(scenario(), geometry(elev), ROOT / f"e{elev:.0f}",
                    seed=SEED, preset=PRESET, grid=g, plan=plan,
                    fft_backend=backend)


def run(elev, n_trials, backend):
    with start_field("gauss", WS_FRAC):
        c = campaign(elev, backend)
        (ROOT / f"e{elev:.0f}" / "cell.json").write_text(json.dumps({
            "start": "gauss", "ws_frac": WS_FRAC, "elev": elev,
            "n": c.grid.n, "side_m": c.grid.size_m, "preset": PRESET}))
        c.run(n_trials, progress=True)


class _Eta:
    """eta of every divergence for one stored trial (picklable)."""

    def __init__(self, psis, o_vacs):
        self.psis, self.o_vacs = psis, o_vacs

    def __call__(self, rec):
        return np.array([reciprocity_overlap(rec.array, p) / o
                         for p, o in zip(self.psis, self.o_vacs)])


def read(elev, backend="cupy"):
    c = campaign(elev, backend)          # the stored backend (fingerprint)
    g, plan = c.grid, c.plan
    mask = super_gaussian_boundary(g.n, PRESETS[PRESET].boundary_width_frac)
    cdtype = np.complex64
    grounds = [scenario(d).ground for d in DIVERGENCES]
    with start_field("gauss", WS_FRAC):
        o_vacs = [space_vacuum_baseline(
            g, plan, LAM, mask, cdtype,
            psi_tx=_ground_transmit_mode(gr, g, dtype=cdtype))[1]
            for gr in grounds]
    patch = c.load(1, fields=True).patch
    psis = [_transmit_mode_crop(gr, g, patch, cdtype) for gr in grounds]
    eta = np.asarray(c.map_trials(_Eta(psis, o_vacs)))
    np.save(HERE / f"arm_b_eta_e{elev:.0f}.npy", eta)
    return eta


def summarise(eta):
    rows = {}
    for i, d in enumerate(DIVERGENCES):
        e = eta[:, i]
        db = -10 * np.log10(e)
        rows[label(d)] = {
            "div_urad": 0.0 if d is None else d * 1e6,
            "mean_loss_db": float(-10 * np.log10(e.mean())),
            "sigma2_I": float(e.var() / e.mean() ** 2),
            "p50_db": float(np.percentile(db, 50)),
            "p5_db": float(np.percentile(db, 95)),
            "p1_db": float(np.percentile(db, 99)),
            "frac_gt_1": float(np.mean(e > 1))}
    return rows


def plot(etas):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, len(etas), figsize=(6 * len(etas), 4.5))
    for ax, (elev, eta) in zip(np.atleast_1d(axs), etas.items()):
        for i, d in enumerate(DIVERGENCES):
            db = np.sort(-10 * np.log10(eta[:, i]))
            ax.plot(db, 1 - np.arange(db.size) / db.size,
                    label=label(d) if d is None else
                    f"{label(d)} ({d * 1e6:.0f} urad)")
        ax.set(yscale="log", xlabel="turbulence loss (dB)",
               ylabel="P(loss > x)", title=f"elevation {elev:.0f} deg")
        ax.grid(alpha=0.3)
        ax.legend()
    fig.tight_layout()
    (HERE / "figures").mkdir(exist_ok=True)
    fig.savefig(HERE / "figures" / "arm_b_ccdf.png", dpi=120)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument("part", choices=("run", "read"))
    ap.add_argument("--elev", type=float, nargs="+", default=list(ELEVATIONS))
    ap.add_argument("--trials", type=int, default=N_TRIALS)
    ap.add_argument("--backend", default="cupy")
    a = ap.parse_args()
    if a.part == "run":
        for e in a.elev:
            run(e, a.trials, a.backend)
    else:
        etas = {e: read(e) for e in a.elev}
        summary = {f"e{e:.0f}": summarise(v) for e, v in etas.items()}
        (HERE / "arm_b_summary.json").write_text(json.dumps(summary, indent=1))
        print(json.dumps(summary, indent=1))
        plot(etas)
