"""Generate and cache the slow wave-optics data for the plots.

Run once:  PYTHONPATH=. python presentation/gen_data.py
It writes presentation/data/plotdata.npz. The plot scripts read that, so you
can re-style the figures without re-running the field solves.

The one hero scenario is C.downlink(): a 0.7 m ground telescope coupling to
single-mode fibre, uncompensated. Loss is positive dB.
"""
import time

import numpy as np

import presentation.common as C
from olb.links.downlink import downlink_budget
from olb.models.waveoptics import run_fidelity2
from olb.waveoptics import Threader
from olb.waveoptics.turbulence.run import propagate_turbulent_field

N_TRIALS = 1000
PRESET = "standard"
N_DRAW = 8000
ELEVATIONS = (90.0, 30.0)        # zenith and the hero mid-pass elevation
RNG = np.random.default_rng(0)


def _term(budget, category):
    return [t for t in budget.terms if t.category == category][0]


def budgets_and_penalties():
    """Ladder + Monte-Carlo-matching data on the uncompensated SMF downlink.

    The 'penalty' is the turbulent fibre-coupling loss (dB): analytic FAST at
    fidelity 1, wave optics at fidelity 2, and the mean-only point at
    fidelity 0. All three describe the same coupling term, so they share an
    axis.
    """
    scn = C.downlink()
    out = {}

    for elev in ELEVATIONS:
        geom = C.geometry(elev)
        b0 = downlink_budget(scn, geom, fidelity=0)
        b1 = downlink_budget(scn, geom, fidelity=1)

        t0 = time.time()
        bundle = run_fidelity2(scn, geom, n_trials=N_TRIALS, preset=PRESET,
                               seed=11, hs=C.HS, cn2_profile=C.CN2,
                               threader=Threader())
        wall = time.time() - t0
        b2 = downlink_budget(scn, geom, fidelity=2, wave=bundle)
        print(f"elev {elev:>4.0f} deg : fid2 {N_TRIALS} trials in {wall:.1f}s "
              f"({wall / N_TRIALS * 1000:.0f} ms/trial)")

        c0 = float(np.asarray(_term(b0, "coupling").mean_db))        # mean-only point
        c1 = _term(b1, "coupling").sample_db(N_DRAW, RNG)            # FAST fade
        c2 = _term(b2, "turbulence").sample_db(N_DRAW, RNG)          # wave-optics fade

        tag = str(int(elev))
        out[f"fid0_coup_{tag}"] = np.array(c0)
        out[f"fid1_coup_{tag}"] = np.asarray(c1, float).ravel()
        out[f"fid2_coup_{tag}"] = np.asarray(c2, float).ravel()
        out[f"wall_{tag}"] = np.array(wall)

    return out


def waterfall_terms():
    """Per-term mean dB for the hero SMF budget at fidelity 1 (model of record)."""
    scn = C.downlink()
    b = downlink_budget(scn, C.geometry(), fidelity=1)
    names = [t.name for t in b.terms]
    means = [float(np.asarray(t.mean_db)) for t in b.terms]
    cats = [t.category for t in b.terms]
    mc = b.monte_carlo(N_DRAW, rng=RNG, availabilities=(0.90,))
    fade = mc["fade_db"]
    return {
        "wf_names": np.array(names, dtype=object),
        "wf_means": np.array(means, float),
        "wf_cats": np.array(cats, dtype=object),
        "wf_tx_dbm": np.array(float(b.tx_power_dbm)),
        "wf_total": np.array(float(b.total_loss_db())),
        "wf_fade90": np.array(float(fade[0.90]) if fade else np.nan),
        "wf_rx_sens_dbm": np.array(float(b.rx_sensitivity_dbm)),
    }


def speckle_and_screens():
    """One turbulent downlink field (intensity + phase) + the phase-screen plan."""
    scn = C.downlink()
    geom = C.geometry()
    out = {}

    F, grid, plan = propagate_turbulent_field(scn, geom, seed=7, trial=0,
                                              preset=PRESET, hs=C.HS,
                                              cn2_profile=C.CN2)
    out["I_turb"] = (np.abs(F.field) ** 2).astype(np.float32)
    out["Ph_turb"] = np.angle(F.field).astype(np.float32)     # wrapped phase [rad]
    out["speckle_extent"] = np.array([float(F.xvalues[0]), float(F.xvalues[-1]),
                                      float(F.yvalues[0]), float(F.yvalues[-1])])
    out["rx_aperture_m"] = np.array(scn.rx_terminal.aperture_m)

    out["plan_r0_m"] = np.asarray(plan.r0_m, float)
    out["plan_z_m"] = np.asarray(plan.z_m, float)
    out["grid_n"] = np.array(grid.n)
    out["grid_pixel_m"] = np.array(grid.pixel_m)
    out["r0_total_m"] = np.array(float(plan.r0_total_m))
    return out


def main():
    data = {}
    print("[1/3] budgets + penalties (wave-optics, slow) ...")
    data.update(budgets_and_penalties())
    print("[2/3] waterfall terms ...")
    data.update(waterfall_terms())
    print("[3/3] speckle + screens ...")
    data.update(speckle_and_screens())
    path = C.datapath("plotdata.npz")
    np.savez(path, **data)
    print("wrote", path)


if __name__ == "__main__":
    main()
