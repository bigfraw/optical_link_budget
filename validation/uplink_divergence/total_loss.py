"""The TOTAL uplink loss (geometric + turbulence) of the Arm B campaigns.

Two routes to the total, trial by trial, from the stored Arm B eta:

  BUDGET      the analytic geometric Term (geometric_loss_term, a smooth
              far-field Gaussian into the receive aperture, plus the opt-in
              launch truncation Term, as uplink_budget(fidelity=2) adds them)
              + the turbulence loss -10 log10(eta).
  CONSISTENT  the EXACT vacuum of the clipped launch at the satellite
              + the same turbulence loss. eta is normalised to that exact
              vacuum (the on-axis overlap of the clipped psi_tx), so this
              route is self-consistent.

The exact vacuum is the Fraunhofer on-axis field of the launch (Goodman,
Introduction to Fourier Optics, Eq. (4-17), ISBN 978-0974707723):

    E(0) = (1 / (lambda Z)) sum( psi(r) exp(i k r^2 / (2 Z)) ) dA

with the power fraction |E(0)|^2 pi a_rx^2 (the receive aperture is tiny
against the far-field footprint), times the clip transmission
1 - exp(-2 a_tx^2 / w^2) of the Gaussian at the launch aperture (Siegman,
Lasers, ISBN 978-0935702118). The difference of the two routes is the
hard-clip far-field ripple that the analytic Term does not hold.

Run from the repository root as a module, after arm_b read:
    -m validation.uplink_divergence.total_loss
"""
import json
import warnings
from pathlib import Path

import numpy as np

from olb.links.uplink import TX_TRUNCATION_MIN_DB
from olb.models.gaussian_efficiency import tx_gaussian_efficiency_term
from olb.models.geometric import geometric_loss_term
from olb.waveoptics.field import Begin
from olb.waveoptics.turbulence.run import _ground_transmit_mode

from .arm_b import cell_grid
from .config import D, DIVERGENCES, ELEVATIONS, LAM, WAIST, geometry, scenario
from .gate0 import label

warnings.simplefilter("ignore")
HERE = Path(__file__).parent
# The ordinal blue ramp of the dataviz reference palette, light -> dark.
COLOURS = ["#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#184f95", "#0d366b"]


def analytic_geo_db(div, elev):
    scn, geo = scenario(div), geometry(elev)
    db = float(np.asarray(geometric_loss_term(scn, geo).mean_db))
    eff = tx_gaussian_efficiency_term(scn, geo)
    if eff.mean_db > TX_TRUNCATION_MIN_DB:
        db += float(np.asarray(eff.mean_db))
    return db


def exact_geo_db(div, elev, g):
    scn, Z = scenario(div), float(np.asarray(geometry(elev).slant_range_m))
    psi = _ground_transmit_mode(scn.ground, g)       # sum |psi|^2 = 1
    Y, X = Begin(g.size_m, LAM, g.n).mgrid_cartesian
    kern = np.exp(1j * np.pi / (LAM * Z) * (X ** 2 + Y ** 2))
    E0 = (psi * kern).sum() * g.pixel_m / (LAM * Z)   # psi / dx is the density
    a_rx = scn.space.aperture_m / 2
    clip = 1 - np.exp(-2 * (D / 2) ** 2 / WAIST ** 2)
    return float(-10 * np.log10(abs(E0) ** 2 * np.pi * a_rx ** 2 * clip))


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axs = plt.subplots(3, 2, figsize=(12, 13))
    out = {}
    for col, elev in enumerate(ELEVATIONS):
        g, _plan = cell_grid(elev)
        eta = np.load(HERE / f"arm_b_eta_e{elev:.0f}.npy")
        turb = -10 * np.log10(eta)
        # The geometric loss on a fine divergence sweep: the ripple.
        sweep = np.r_[0, np.arange(10, 211, 2)] * 1e-6
        ana = [analytic_geo_db(t or None, elev) for t in sweep]
        exa = [exact_geo_db(t or None, elev, g) for t in sweep]
        ax = axs[2, col]
        ax.plot(sweep * 1e6, ana, color="#2a78d6", lw=2, label="analytic Term")
        ax.plot(sweep * 1e6, exa, color="#eb6834", lw=2,
                label="exact clipped launch")
        ax.set(xlabel="divergence (urad)", ylabel="vacuum geometric loss (dB)",
               title=f"{elev:.0f} deg: geometric loss, vacuum")
        ax.legend(frameon=False)
        rows = {}
        for i, d in enumerate(DIVERGENCES):
            geo_a, geo_e = analytic_geo_db(d, elev), exact_geo_db(d, elev, g)
            name = label(d) if d is None else f"{d * 1e6:.0f} urad"
            for r, (route, geo) in enumerate((("budget", geo_a),
                                              ("consistent", geo_e))):
                tot = np.sort(geo + turb[:, i])
                axs[r, col].plot(tot, np.arange(1, tot.size + 1) / tot.size,
                                 color=COLOURS[i], lw=2, label=name)
            tb, te = geo_a + turb[:, i], geo_e + turb[:, i]
            rows[name] = {
                "geo_analytic_db": geo_a, "geo_exact_db": geo_e,
                **{f"{k}_{q}": float(np.percentile(v, p)) for k, v in
                   (("budget", tb), ("consistent", te))
                   for q, p in (("p50_db", 50), ("p95_db", 95), ("p99_db", 99))},
                "consistent_mean_db": float(-10 * np.log10(np.mean(10 ** (-te / 10))))}
        out[f"e{elev:.0f}"] = rows
        for r, route in enumerate(("budget: analytic geometric + turbulence",
                                   "consistent: exact clipped vacuum + turbulence")):
            axs[r, col].set(xlabel="total loss (dB)", ylabel="P(loss <= x)",
                            title=f"{elev:.0f} deg, {route}")
            axs[r, col].legend(frameon=False, title="divergence")
    for ax in axs.flat:
        ax.grid(alpha=0.25)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    # One x range per column for the two CDF rows, so they compare by eye.
    for col in range(2):
        lo = min(axs[r, col].get_xlim()[0] for r in (0, 1))
        hi = max(axs[r, col].get_xlim()[1] for r in (0, 1))
        for r in (0, 1):
            axs[r, col].set_xlim(lo, hi)
    fig.tight_layout()
    fig.savefig(HERE / "figures" / "total_loss_cdf.png", dpi=120)
    (HERE / "total_loss_summary.json").write_text(json.dumps(out, indent=1))
    for e, rows in out.items():
        print(e)
        for name, r in rows.items():
            print(f"  {name:9s} geo analytic {r['geo_analytic_db']:6.2f} exact "
                  f"{r['geo_exact_db']:6.2f} | total p50/p95/p99 budget "
                  f"{r['budget_p50_db']:5.1f}/{r['budget_p95_db']:5.1f}/"
                  f"{r['budget_p99_db']:5.1f}  consistent "
                  f"{r['consistent_p50_db']:5.1f}/{r['consistent_p95_db']:5.1f}/"
                  f"{r['consistent_p99_db']:5.1f}")


if __name__ == '__main__':
    main()
