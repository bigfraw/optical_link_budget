'''
Sweep one space-link budget over an elevation grid: one Budget for each angle.

Some budget Terms cannot take an elevation ARRAY. The FAST coupling Term and the
gamma-gamma downlink Term each model ONE line of sight: FAST runs one Monte Carlo
for one geometry, and the gamma-gamma Term carries one (alpha, beta) pair. So a
budget that holds one of them refuses an array elevation and asks the caller to
loop. Other Terms (geometric spread, extinction, pointing, the lognormal
scintillation) DO vectorise, so the same budget accepts an array for some front
ends and raises for others.

This helper gives ONE uniform way that always works: it builds a scalar-elevation
geometry for each angle, calls the scenario's budget function, and returns the
Budget of each angle. So a sweep reads the same whatever the Term does inside.

    for each elevation: build CircularOrbit(altitude, elevation), call the
                        matching budget function, keep (elevation, Budget).

The altitude comes from the scenario channel, so the common call needs no
geometry argument. The point-ahead angle and the slant range then follow from
each elevation, which the uplink point-ahead Term needs.

WHY THIS MODULE IS AT THE TOP LEVEL. A sweep is CROSS-CUTTING: it is not the
physics of one link. olb.links holds the per-link physics, so this module sits
ABOVE it in the dependency order. It reuses olb.multidetector for the family and
direction dispatch; nothing in olb.links reads it back.

This is a SPACE-link helper. A terrestrial link has no elevation axis (its range
is the horizontal path length), so this raises for a TerrestrialScenario.
'''

import numpy as np

from .geometry import CircularOrbit
from .multidetector import _budget_function


def budgets_vs_elevation(scenario, elevations, *, geometry_factory=None,
                         **kwargs):
    '''
    Assemble the scenario budget at each elevation: one Budget for each angle.

    Parameters:
        scenario : SpaceScenario
            The link case (uplink, downlink, or retro). A TerrestrialScenario
            has no elevation axis and raises.
        elevations : float or array-like
            Elevation angle(s) above the horizon [deg]. A scalar gives one
            budget.
        geometry_factory : callable, optional
            A function elevation_deg -> geometry, for a geometry other than the
            default. The default builds CircularOrbit(scenario.channel.
            altitude_m, elevation_deg) for each angle.
        **kwargs :
            Passed to the budget function unchanged (for example fidelity,
            turbulence, n_samples, tau_zenith, scint_model, wave).

    Returns:
        list of (float, Budget)
            One (elevation_deg, Budget) pair for each angle, in the elevations
            order. Reduce it as the task needs, for example
            totals = [(e, float(b.total_loss_db())) for e, b in sweep].

    Raises:
        TypeError
            If the scenario is terrestrial (no elevation axis).
        (the budget function's own errors pass through)
            A per-angle error is NOT caught. For example an uncorrected uplink
            at fidelity 0 raises, and it raises here too.
    '''
    if not hasattr(scenario, "ground"):        # a SpaceScenario has ground/space
        raise TypeError(
            "budgets_vs_elevation takes a SpaceScenario. A terrestrial link has "
            "no elevation axis; sweep the horizontal path length with a "
            "HorizontalPath geometry instead.")

    budget_of = _budget_function(scenario)
    if geometry_factory is None:
        altitude_m = scenario.channel.altitude_m
        def geometry_factory(elevation_deg):
            return CircularOrbit(altitude_m, elevation_deg)

    out = []
    for e in np.atleast_1d(np.asarray(elevations, dtype=float)):
        e = float(e)
        out.append((e, budget_of(scenario, geometry_factory(e), **kwargs)))
    return out


if __name__ == '__main__':
    from .scenario import Channel, SpaceScenario
    from .terminal import Aperture, Terminal, Transmitter

    # A downlink: the ground terminal is a bare aperture with NO detector, so
    # the budget builds the scintillation Term (the gamma-gamma path of the
    # I-1 regression below lives on that branch).
    down = SpaceScenario(
        ground=Terminal(aperture_m=0.4, wavelength_m=1550e-9),
        space=Terminal(aperture_m=0.1, wavelength_m=1550e-9,
                       transmitter=Transmitter(waist_m=0.04, power_dbm=30.0)),
        direction="downlink", channel=Channel())

    # --- the sweep gives one Budget for each angle --------------------------
    elevs = [30.0, 60.0, 90.0]
    sweep = budgets_vs_elevation(down, elevs, fidelity=0)
    assert [e for e, _ in sweep] == elevs, [e for e, _ in sweep]

    # The loss falls toward the zenith: a shorter slant path, less airmass.
    totals = [float(b.total_loss_db()) for _, b in sweep]
    assert totals[0] > totals[1] > totals[2], totals

    # --- the I-1 regression: gamma-gamma refuses an ARRAY elevation, but the
    #     sweep runs it one angle at a time -----------------------------------
    try:
        down_budget = _budget_function(down)
        down_budget(down, CircularOrbit(down.channel.altitude_m, [30.0, 60.0]),
                    fidelity=0, scint_model="gamma_gamma")
        raise AssertionError("an array elevation must raise for gamma-gamma")
    except NotImplementedError:
        pass
    gg = budgets_vs_elevation(down, [30.0, 60.0], fidelity=0,
                              scint_model="gamma_gamma")
    assert len(gg) == 2, len(gg)
    assert all(isinstance(e, float) for e, _ in gg)

    # --- a scalar elevation gives a length-one list -------------------------
    one = budgets_vs_elevation(down, 45.0, fidelity=0)
    assert len(one) == 1 and one[0][0] == 45.0, one

    # --- a custom geometry_factory is honoured ------------------------------
    seen = []
    def factory(e):
        seen.append(e)
        return CircularOrbit(600e3, e)
    budgets_vs_elevation(down, [25.0, 50.0], geometry_factory=factory,
                         fidelity=0)
    assert seen == [25.0, 50.0], seen

    # --- a terrestrial scenario has no elevation axis -----------------------
    from .scenario import TerrestrialChannel, TerrestrialScenario
    terr = TerrestrialScenario(
        near=Terminal(aperture_m=0.2, wavelength_m=1550e-9,
                      transmitter=Transmitter(waist_m=0.02, power_dbm=30.0)),
        far=Terminal(aperture_m=0.2, wavelength_m=1550e-9,
                     detector=Aperture()),
        channel=TerrestrialChannel(path_length_m=3e3, attenuation_db_per_km=0.5,
                                   cn2=3e-16))
    try:
        budgets_vs_elevation(terr, [30.0, 60.0])
        raise AssertionError("a terrestrial scenario must raise")
    except TypeError as exc:
        assert "elevation axis" in str(exc), str(exc)

    print("elevation sweep, downlink 600 km (analytic):")
    for e, b in sweep:
        print(f"  {e:5.1f} deg   total {float(b.total_loss_db()):8.3f} dB   "
              f"({len(b.terms)} terms)")
    print("self-check passed")
