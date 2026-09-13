"""The shared case of the point-ahead screen study.

THE CASE. It is the hero downlink of `validation/waveoptics_ao/` turned around:
a 1550 nm uplink from a 700 mm ground terminal to a 100 mm space terminal at
500 km. The ground terminal senses the downlink beacon, so the scenario carries
a `DownlinkBeacon` pre-compensation source. The ground hardware, the site, the
preset, the seed and the outer scale are the values of that study.

WHY AN UPLINK. The point-ahead angle applies to the uplink only. The ground
terminal senses the beacon that comes from where the satellite WAS, and it must
launch to where the satellite WILL BE. The two directions differ by the
point-ahead angle of `olb.geometry.CircularOrbit.point_ahead_rad`.

THE SCREEN PLANS. This module gives the plans that both scripts compare:
  - the PRODUCTION plan of `turbulent_grid` at the standard preset;
  - an OVERRIDE plan of exactly n screens, through a preset copy;
  - the GROUND-SPLIT plan, which cuts the lowest screen into sub-screens;
  - the EQUAL-ANISOPLANATIC-WEIGHT plan, which cuts the slab at equal shares of
    the Stone weight instead of equal shares of the Rytov weight.

Sources:
- J. Stone, P. H. Hu, S. P. Mills and S. Ma, "Anisoplanatic effects in
  finite-aperture optical systems," J. Opt. Soc. Am. A 11(1), 347-357 (1994),
  DOI 10.1364/JOSAA.11.000347.
- Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196. Ch. 12, Eq. (14):
  the plane-parallel airmass sec(zeta). Ch. 8, Eq. (20): the Rytov path weight.
- D. L. Fried, J. Opt. Soc. Am. 56, 1372 (1966), DOI 10.1364/JOSA.56.001372.
  The Fried parameter r0.
"""

import json
import os
import warnings
from dataclasses import replace

import numpy as np

from olb.geometry import CircularOrbit
from olb.scenario import Channel, DownlinkBeacon, Site, SpaceScenario
from olb.terminal import SMF, Terminal, Transmitter
from olb.turbulence.andrews.beam import wavenumber
from olb.turbulence.andrews.paths import sec_zeta
from olb.turbulence.profiles import get_c2n
# The planner internals are PRIVATE. A validation script reads them, because
# the study measures the planner itself. validation/tail_convergence/ reads the
# same names.
from olb.waveoptics.turbulence.sampling import (DEFAULT_H_TOP_M, PRESETS,
                                                _cumtrapz,
                                                _integration_heights,
                                                _screen_rytov, _composite_r0,
                                                turbulent_grid)
from olb.waveoptics.turbulence.screens import screen_r0
from olb.waveoptics.turbulence.sampling import ScreenPlan

HERE = os.path.dirname(os.path.abspath(__file__))

# The case. The values repeat validation/waveoptics_ao/waveoptics_ao.py.
LAM = 1550e-9
ALTITUDE_M = 500e3
GROUND_APERTURE_M = 0.70
SEED = 20260911
PRESET = "standard"

# The outer scale of the screens, in m. The owner fixed it on 2026-09-05
# (backlog 2-P5). L0 = inf is grid-dependent.
L0_M = 25.0

# The apertures and the angles of both scripts.
APERTURES_M = (0.4, 0.7, 1.0)
ANGLES_ARCSEC = (2.0, 5.0, 10.0)
ELEVATIONS_DEG = (20.0, 30.0, 60.0, 90.0)

ARCSEC = np.pi / (180.0 * 3600.0)


def hero_uplink(elevation_deg):
    """Build the hero uplink and the orbit of one elevation.

    Args:
        elevation_deg: the elevation of the line of sight, in deg.

    Returns:
        The pair (SpaceScenario, CircularOrbit).
    """
    site = Site(cn2_ground=1.7e-14, wind_rms_m_s=21.0, outer_scale_m=L0_M)
    channel = Channel(site=site, altitude_m=ALTITUDE_M)
    ground = Terminal(aperture_m=GROUND_APERTURE_M, wavelength_m=LAM,
                      pointing_jitter_rad=2e-6,
                      detector=SMF(sensitivity_dbm=-45.0),
                      transmitter=Transmitter(waist_m=0.35, power_dbm=30.0))
    space = Terminal(aperture_m=0.10, wavelength_m=LAM,
                     pointing_jitter_rad=1e-6,
                     transmitter=Transmitter(waist_m=0.04, power_dbm=30.0))
    scn = SpaceScenario(ground=ground, space=space, direction="uplink",
                        channel=channel, precompensation=DownlinkBeacon())
    geom = CircularOrbit(altitude_m=ALTITUDE_M,
                         elevation_deg=[float(elevation_deg)])
    return scn, geom


def site_cn2(scn):
    """Give the site Hufnagel-Valley zenith Cn2 callable cn2(h).

    It repeats the default of
    `olb.waveoptics.turbulence.sampling._plan_space`. See
    `olb.turbulence.profiles.get_c2n`.
    """
    site = scn.channel.site

    def cn2(h):
        return get_c2n(h, site.wind_rms_m_s, site.cn2_ground)

    return cn2


def production_plan(scn, geom, preset=PRESET):
    """Give the (grid, plan, warnings) of the production sizer.

    Args:
        scn:    the SpaceScenario.
        geom:   the CircularOrbit.
        preset: a preset name or a QualityPreset.

    Returns:
        The triple (GridSpec, ScreenPlan, a tuple of warning texts).
    """
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        grid, plan, _ = turbulent_grid(scn, geom, preset=preset, L0_m=L0_M)
    return grid, plan, tuple(str(w.message) for w in caught)


def override_plan(scn, geom, n_screens, preset=PRESET):
    """Give the plan of EXACTLY n_screens screens, of the production shape.

    THE OVERRIDE. `_plan_space_continuous` takes the count
    N = max(min_screens, ceil(sigma_R^2 / sigma2_r_screen_max)) and it cuts the
    slab into N slabs of equal plane-wave Rytov weight. A preset copy with the
    cap removed (sigma2_r_screen_max = 1e9) and the floor at n gives exactly n
    screens of the production shape. The recipe comes from
    `validation/terrestrial_screen_count/screen_count_sweep.py`.

    Args:
        scn:       the SpaceScenario.
        geom:      the CircularOrbit.
        n_screens: the wanted screen count.
        preset:    the base preset name.

    Returns:
        The pair (GridSpec, ScreenPlan) of the override.

    Raises:
        AssertionError: the plan does not hold exactly n_screens screens.
    """
    base = PRESETS[preset] if isinstance(preset, str) else preset
    over = replace(base, name=f"{base.name}_n{int(n_screens)}",
                   sigma2_r_screen_max=1e9, min_screens=int(n_screens))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        grid, plan, _ = turbulent_grid(scn, geom, preset=over, L0_m=L0_M)
    assert plan.z_m.size == int(n_screens), (plan.z_m.size, n_screens)
    return grid, plan


def equal_aniso_plan(scn, geom, n_screens, exponent, preset=PRESET,
                     centroid="cn2"):
    """Cut the slab into n_screens slabs of equal ANISOPLANATIC weight.

    THE SHAPE IS THE PRODUCTION SHAPE. The logic copies
    `olb.waveoptics.turbulence.sampling._plan_space_continuous`: it integrates
    the slant Cn2 density m(h) = Cn2(h) sec over the same internal height grid,
    it cuts the slab at equal shares of a cumulative weight, and it puts each
    screen at the Cn2-weighted centroid of its slab.

    ONLY THE WEIGHT CHANGES. The planner cuts at equal shares of the Rytov
    weight w(h) = 2.25 k^(7/6) Cn2(h) sec (h sec)^(5/6) (Andrews and Phillips,
    DOI 10.1117/3.626196, Ch. 8, Eq. (20)). This function cuts at equal shares
    of the STONE weight

        a(h) = Cn2(h) sec * (h sec)^exponent

    which is the small-angle limit of the Stone height integral: the integrand
    of Eq. (29) of Stone et al. (1994), DOI 10.1364/JOSAA.11.000347, is
    Cn2(h) sec I(S theta / R) with S = h sec, and I(beta) goes as beta^(5/3)
    when no mode is removed and as beta^2 when the piston and the tilt are
    removed. So the caller passes exponent = 5/3 or exponent = 2.

    Args:
        scn:       the SpaceScenario.
        geom:      the CircularOrbit.
        n_screens: the wanted screen count.
        exponent:  the height exponent of the weight (5/3 or 2).
        preset:    the base preset name, for the grid only.
        centroid:  "cn2" puts the screen at the Cn2-weighted centroid, the
                   production rule. "aniso" puts it at the height that holds
                   the ANISOPLANATIC weight of the slab:

                       (z_g)^p = INT a dh / INT m dh

                   which is the power mean of the ground distance over the
                   slab. A delta layer at that distance reproduces the
                   small-angle Stone integral of the slab exactly, because
                   I(beta) goes as beta^p there.

    Returns:
        A ScreenPlan of n_screens screens.
    """
    lam = scn.tx_terminal.wavelength_m
    k = wavenumber(lam)
    h_top = DEFAULT_H_TOP_M
    elevation = float(np.min(np.asarray(geom.elevation_deg, dtype=float)))
    sec = float(sec_zeta(elevation))
    n_s = int(n_screens)

    h = _integration_heights(h_top)
    cn2_h = np.asarray(site_cn2(scn)(h), dtype=float)
    z_of_h = (h_top - h) * sec
    m_dens = cn2_h * sec
    a_dens = m_dens * (h * sec) ** float(exponent)
    w_dens = _screen_rytov(k, m_dens, np.maximum(h * sec, 0.0))

    a_cum = _cumtrapz(a_dens, h)
    m_cum = _cumtrapz(m_dens, h)
    w_cum = _cumtrapz(w_dens, h)
    zm_cum = _cumtrapz(z_of_h * m_dens, h)

    targets = float(a_cum[-1]) * np.arange(1, n_s) / n_s
    h_edges = np.concatenate(([0.0], np.interp(targets, a_cum, h),
                              [float(h_top)]))

    def _seg(cum, lo, hi):
        return float(np.interp(hi, h, cum) - np.interp(lo, h, cum))

    cn2_int = np.empty(n_s)
    z = np.empty(n_s)
    s2 = np.empty(n_s)
    for j in range(n_s):
        lo, hi = h_edges[j], h_edges[j + 1]
        i = n_s - 1 - j                       # the z-ascending index
        m_slab = _seg(m_cum, lo, hi)
        cn2_int[i] = m_slab
        if m_slab <= 0.0:
            z[i] = (float(h_top) - 0.5 * (lo + hi)) * sec
        elif centroid == "aniso":
            # The power mean of the ground distance, then back to z.
            zg_eff = (_seg(a_cum, lo, hi) / m_slab) ** (1.0 / float(exponent))
            z[i] = float(h_top) * sec - zg_eff
        else:
            z[i] = _seg(zm_cum, lo, hi) / m_slab
        s2[i] = _seg(w_cum, lo, hi)

    if np.any(np.diff(z) <= 0.0):
        raise ValueError("equal_aniso_plan: the screen distances are not "
                         "increasing. Lower the screen count.")
    r0 = screen_r0(cn2_int, lam)
    return ScreenPlan(z_m=z, cn2_int_m13=cn2_int, r0_m=r0, sigma2_r=s2,
                      z_total_m=float(h_top) * sec,
                      r0_total_m=_composite_r0(r0), direction="down")


def ground_distances(plan):
    """Give the GROUND distance of each screen, in m.

    A space plan runs DOWN the atmosphere, so `plan.z_m` is the distance from
    the TOP of the slab. The ground distance is z_total - z_m. It equals
    S = h sec(zeta), the slant distance that Stone's Eq. (29) reads.
    """
    return np.asarray(plan.z_total_m - plan.z_m, dtype=float)


def plan_set(scn, geom, counts, split_n=4):
    """Give the named plans of one elevation.

    Args:
        scn:     the SpaceScenario.
        geom:    the CircularOrbit.
        counts:  the override screen counts.
        split_n: the sub-screen count of the ground-split plan.

    Returns:
        The pair (GridSpec, a dict name -> ScreenPlan). The grid is the
        PRODUCTION grid of the standard preset.
    """
    # split_bottom_screen is the ground-split builder of
    # validation/tail_convergence/. It is imported here, not copied.
    from validation.tail_convergence.tail_convergence import \
        split_bottom_screen

    grid, prod, _ = production_plan(scn, geom)
    plans = {"production": prod}
    for n in counts:
        plans[f"n{int(n)}"] = override_plan(scn, geom, n)[1]
    plans[f"ground_split_x{int(split_n)}"] = split_bottom_screen(
        scn, geom, prod, int(split_n))
    return grid, plans


def pyplot():
    """Give the matplotlib pyplot module with the Agg backend, or None."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib is not installed; the figures are skipped.",
              flush=True)
        return None
    return plt


def fig_dir():
    """Give the figure directory of this study, and make it."""
    path = os.path.join(HERE, "figures")
    os.makedirs(path, exist_ok=True)
    return path


def log_maker(name):
    """Make the log file of one script, and give its `say` function."""
    path = os.path.join(HERE, f"{name}.log")
    with open(path, "w", encoding="utf-8"):
        pass

    def say(text=""):
        print(text, flush=True)
        with open(path, "a", encoding="utf-8") as fh:
            fh.write(text + "\n")

    return say, path


def write_results(name, payload):
    """Write the JSON record of one script."""
    path = os.path.join(HERE, f"{name}_results.json")
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=1)
    return path


if __name__ == '__main__':
    # The self-check. It builds the case and it prints the plans of 30 deg.
    scn, geom = hero_uplink(30.0)
    grid, plan, warns = production_plan(scn, geom)
    theta = float(np.atleast_1d(geom.point_ahead_rad)[0])
    print(f"grid {grid.n} px, {grid.pixel_m * 1e3:.3f} mm pixel, "
          f"side {grid.size_m:.3f} m")
    print(f"point-ahead angle {theta / ARCSEC:.3f} arcsec "
          f"({theta * 1e6:.3f} urad)")
    print(f"production plan: {plan.z_m.size} screens, r0_total "
          f"{plan.r0_total_m * 100:.2f} cm, z_total {plan.z_total_m / 1e3:.2f} km")
    zg = ground_distances(plan)
    assert np.all(zg > 0.0) and np.all(np.diff(zg) < 0.0), zg
    for n in (5, 9, 15):
        _, p = override_plan(scn, geom, n)
        assert p.z_m.size == n
        assert abs(p.cn2_int_m13.sum() / plan.cn2_int_m13.sum() - 1.0) < 1e-3
    a = equal_aniso_plan(scn, geom, 9, 2.0)
    assert a.z_m.size == 9
    assert abs(a.cn2_int_m13.sum() / plan.cn2_int_m13.sum() - 1.0) < 1e-2, \
        (a.cn2_int_m13.sum(), plan.cn2_int_m13.sum())
    print("equal-aniso ground distances [km]: "
          + " ".join(f"{v / 1e3:.2f}" for v in ground_distances(a)))
    print("production ground distances [km]: "
          + " ".join(f"{v / 1e3:.2f}" for v in zg))
    for w in warns:
        print(f"warning: {w}")
    print("self-check passed")
