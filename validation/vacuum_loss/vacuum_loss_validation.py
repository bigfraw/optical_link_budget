'''
Cross-check the ANALYTIC geometric loss against the WAVE-OPTICS vacuum loss.

A fidelity-2 budget can get its no-turbulence geometric loss two ways:

  - ANALYTIC (a closed form): geometric_loss_term + tx_gaussian_efficiency_term.
    The geometric Term is a far-field diffraction spread, and the transmit
    efficiency is a far-field truncation. They cost microseconds and have no
    grid.
  - WAVE OPTICS (a field solve): the vacuum propagation, _full_vacuum_loss_db.
    It is the exact launch-to-collected power ratio of one no-turbulence field
    solve. It costs a full grid FFT chain (about 14 s for a space slant path on
    the 4096^2 full-path grid).

THE FINDING (which sets the fidelity-2 default). A ground-space link is ALWAYS
far field (the Fraunhofer distance of a 0.1 m aperture at 1550 nm is about 6 km,
far shorter than any orbit range), so the analytic geometric loss is exact. The
wave-optics vacuum run on the full slant path is NOT: it cannot resolve the
mm-scale aperture edges over a ~2000 km grid, so its loss SCATTERS by +/- 1 to
4 dB and does not converge at a practical grid size. So for a space link the
analytic Term is BOTH cheaper AND more trustworthy, and run_fidelity2 uses it by
default (vacuum="analytic").

A TERRESTRIAL link keeps the wave vacuum (vacuum="wave") for a STRUCTURAL
reason, NOT because the analytic loss is wrong: the terrestrial turbulence
penalty is turbulent / vacuum on the SAME flat grid, so the wave vacuum is the
exact baseline that cancels the grid. The terrestrial grid is small and well
resolved, so the wave vacuum is cheap and accurate there anyway.

This script SHOWS both halves:

  PART 1 (TERRESTRIAL). The analytic and the WELL-RESOLVED wave loss agree,
  across a collimated far-field path AND a tightly focused short path. This
  VALIDATES the analytic geometric Term against wave optics where wave optics is
  trustworthy. The agreement also shows the terrestrial wave vacuum is accurate
  (it is kept for the penalty-baseline cancellation, above, not for accuracy).

  PART 2 (SPACE). The wave vacuum loss of one space case is measured as the grid
  refines. It SCATTERS around the stable analytic value and does not converge.
  This is the evidence that the space wave vacuum run is grid-noise-limited, so
  the analytic Term is the reference.

The comparison is APPLES TO APPLES. The wave-optics stage[0] "launch" is the
normalised pre-truncation Gaussian (power 1.0), so
_full_vacuum_loss_db = tx_truncation_db + geometric_loss_db. The analytic side
sums the same two parts: the geometric spread Term and the transmit-truncation
Term. Both use an APERTURE receiver, so no fibre coupling enters.

Run from the repo root:
    python -m validation.vacuum_loss.vacuum_loss_validation
'''

import json
import os
import warnings

import numpy as np

from olb import (SpaceScenario, Channel, Site, CircularOrbit, Terminal,
                 Transmitter, Aperture)
from olb.geometry import HorizontalPath
from olb.scenario import TerrestrialScenario, TerrestrialChannel
from olb.waveoptics.grid import GridSpec
from olb.models.geometric import geometric_loss_term
from olb.models.gaussian_efficiency import tx_gaussian_efficiency_term
from olb.models.waveoptics import _full_vacuum_loss_db
from olb.waveoptics.run import propagate_scenario

WAVELENGTH_M = 1550e-9

# The far-field tolerance for the terrestrial cross-check. A long, collimated
# path is far field, so the analytic and the well-resolved wave loss must agree
# inside this many dB. It is loose enough to hold the grid-tail loss (a few
# percent of the power at the receive plane), which the analytic form omits.
TERRESTRIAL_FARFIELD_TOL_DB = 0.35


def _analytic_loss_db(scenario, geometry):
    '''The analytic no-turbulence geometric loss: spread + transmit truncation.'''
    geo = geometric_loss_term(scenario, geometry).mean_db
    trunc = tx_gaussian_efficiency_term(scenario, geometry).mean_db
    return geo, trunc, geo + trunc


def _wave_loss_db(scenario, geometry, grid=None):
    '''The wave-optics vacuum loss: launch to collected power, one field solve.'''
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')       # the coarse full-path grid warns
        result = propagate_scenario(scenario, geometry, grid=grid)
    return _full_vacuum_loss_db(result)


def _terrestrial_case(path_length_m, waist_m):
    '''One horizontal path, an aperture receiver, a chosen launch waist.'''
    near = Terminal(aperture_m=0.2, wavelength_m=WAVELENGTH_M,
                    pointing_jitter_rad=1e-6,
                    transmitter=Transmitter(waist_m=waist_m, power_dbm=30.0))
    far = Terminal(aperture_m=0.2, wavelength_m=WAVELENGTH_M,
                   pointing_jitter_rad=1e-6,
                   detector=Aperture(sensitivity_dbm=-40.0))
    scn = TerrestrialScenario(
        near=near, far=far,
        channel=TerrestrialChannel(path_length_m=path_length_m,
                                   attenuation_db_per_km=0.0, cn2=1e-16))
    return scn, HorizontalPath(path_length_m)


def _space_case(elevation_deg, altitude_m=1500e3):
    '''One ground-space uplink at one elevation, an aperture receiver.'''
    ground = Terminal(
        aperture_m=0.3, wavelength_m=WAVELENGTH_M, pointing_jitter_rad=1e-6,
        transmitter=Transmitter(waist_m=0.05, power_dbm=30.0))
    space = Terminal(
        aperture_m=0.1, wavelength_m=WAVELENGTH_M, pointing_jitter_rad=1e-6,
        detector=Aperture(sensitivity_dbm=-40.0))
    channel = Channel(site=Site(cn2_ground=1.7e-14), altitude_m=altitude_m)
    scn = SpaceScenario(ground=ground, space=space, direction='uplink',
                        channel=channel)
    return scn, CircularOrbit(altitude_m, elevation_deg=elevation_deg)


def _terrestrial_rows(label, waist_m, path_lengths):
    '''One terrestrial table: analytic vs a well-resolved wave loss.'''
    rows = []
    print(f"\n{label} (launch waist {waist_m * 1e3:.0f} mm)")
    print("-" * 74)
    print(f"{'path':>8} {'geo(an)':>9} {'trunc(an)':>10} {'total(an)':>10} "
          f"{'wave':>9} {'diff':>8}")
    for L in path_lengths:
        scn, geom = _terrestrial_case(L, waist_m)
        geo, trunc, an = _analytic_loss_db(scn, geom)
        wave = _wave_loss_db(scn, geom)             # small grid, well resolved
        diff = wave - an
        print(f"{L/1e3:7.2g}k {geo:9.3f} {trunc:10.3f} {an:10.3f} "
              f"{wave:9.3f} {diff:8.3f}")
        rows.append(dict(path_km=L / 1e3, geometric_db=geo, truncation_db=trunc,
                         analytic_db=an, wave_db=wave, diff_db=diff))
    return rows


def _space_convergence(elevation_deg, grid_ns):
    '''One space case: the wave loss as the grid refines, against the analytic.'''
    scn, geom = _space_case(elevation_deg)
    geo, trunc, an = _analytic_loss_db(scn, geom)
    base = GridSpec.for_scenario(scn, geom)          # the sizer's extent + route
    print(f"\nSPACE UPLINK {elevation_deg:.0f} deg: the wave vacuum loss vs the "
          f"grid (analytic = {an:.3f} dB, grid-independent)")
    print("-" * 74)
    print(f"{'grid n':>8} {'mm/px':>8} {'wave dB':>10} {'diff':>8}")
    rows = []
    for n in grid_ns:
        grid = GridSpec(size_m=base.size_m, n=n, scaled=base.scaled)
        wave = _wave_loss_db(scn, geom, grid=grid)
        diff = wave - an
        print(f"{n:>8} {grid.pixel_m * 1e3:>8.2f} {wave:>10.3f} {diff:>8.3f}")
        rows.append(dict(grid_n=n, pixel_mm=grid.pixel_m * 1e3, wave_db=wave,
                         diff_db=diff))
    scatter = max(r['wave_db'] for r in rows) - min(r['wave_db'] for r in rows)
    return dict(analytic_db=an, rows=rows, scatter_db=scatter)


def main():
    results = {}

    # ---- PART 1: TERRESTRIAL, the positive validation of the analytic Term ----
    # A collimated 50 mm waist: far field on the long paths, so analytic == wave.
    results['terrestrial_collimated'] = _terrestrial_rows(
        "TERRESTRIAL collimated (far field on long paths: analytic == wave)",
        0.05, (0.3e3, 1e3, 3e3, 10e3))
    # A tightly focused 8 mm waist (Rayleigh range about 130 m). The analytic
    # Gaussian Term stays accurate through the near-to-far transition, so the
    # aperture-power geometric loss agrees with wave optics here too.
    results['terrestrial_focused'] = _terrestrial_rows(
        "TERRESTRIAL focused (near-to-far transition: analytic stays accurate)",
        0.008, (0.1e3, 0.3e3, 1e3, 3e3))

    # ---- PART 2: SPACE, the wave run is grid-noise-limited ----
    results['space_convergence'] = _space_convergence(
        20.0, (4096, 5120, 6144, 7168))

    # ---- the verdict ----
    far = results['terrestrial_collimated'][-1]     # the longest, most far-field
    space = results['space_convergence']
    print("\n" + "=" * 74)
    print(f"PART 1  terrestrial far field (10 km, collimated): "
          f"|analytic - wave| = {abs(far['diff_db']):.3f} dB "
          f"(tolerance {TERRESTRIAL_FARFIELD_TOL_DB})  -> "
          f"{'PASS' if abs(far['diff_db']) < TERRESTRIAL_FARFIELD_TOL_DB else 'FAIL'}")
    print(f"        This validates the analytic geometric Term against a "
          f"well-resolved wave solve.")
    print(f"PART 2  space wave-loss scatter over 4096..7168 px = "
          f"{space['scatter_db']:.3f} dB, around the stable analytic "
          f"{space['analytic_db']:.3f} dB.")
    print(f"        The space wave vacuum run does NOT converge; the analytic "
          f"Term is the reference.")
    print("=" * 74)
    print("VERDICT: the analytic geometric loss matches a well-resolved wave "
          "solve (Part 1). A space link is deeply far field, yet its full-path "
          "wave vacuum run is grid-noise-limited (Part 2). So a space fidelity-2 "
          "budget uses the ANALYTIC geometric loss (the default). A terrestrial "
          "link keeps the WAVE vacuum because its turbulence penalty is "
          "turbulent / vacuum on the SAME grid (an exact baseline), not because "
          "the analytic loss is wrong.")

    out = os.path.join(os.path.dirname(__file__), 'vacuum_loss_results.json')
    with open(out, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nwrote {out}")

    assert abs(far['diff_db']) < TERRESTRIAL_FARFIELD_TOL_DB, (
        f"the terrestrial far-field analytic and wave loss disagree by "
        f"{abs(far['diff_db']):.3f} dB, past the {TERRESTRIAL_FARFIELD_TOL_DB} "
        f"dB tolerance. The analytic geometric Term is not validated.")


if __name__ == '__main__':
    main()
