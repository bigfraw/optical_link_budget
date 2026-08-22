'''
Clear-sky atmospheric extinction along the slant path.

Beer-Lambert law with a single zenith optical depth that airmass scales:

    T        = exp(-tau_zenith * airmass(elevation))
    loss_db  = -10*log10(T) = (10/ln10) * tau_zenith * airmass(elevation)

The loss is linear in airmass. `tau_zenith` combines molecular (Rayleigh) and
aerosol extinction into one clear-sky number. This is not a MODTRAN
line-by-line model. It is the one-parameter slant attenuation.

DEFAULT_TAU_ZENITH = 0.05 is the near-IR clear-sky zenith optical depth at
1550 nm for a good (dry, high) site. At 1550 nm Rayleigh scattering is
negligible (tau ~ 0.002), so aerosol dominates the extinction. 0.05 gives a
zenith transmittance of exp(-0.05) = 0.95 (~0.22 dB). This value agrees with the
near-IR clear-sky transmittances of ~0.9-0.96 in MODTRAN rural/clear aerosol
runs and typical 1550 nm FSO link budgets. Set `tau_zenith` per site or haze.

The airmass model (plane-parallel secant, 1/sin(elevation)) comes from
fso_spot_size.airmass in the sibling TN-2 analysis repo. It uses the same slant
scaling as the other modules. It diverges at the horizon. Do not use elevation 0.
'''

import numpy as np

from ..results import Term  # run modules with `python -m olb.models.transmittance`
from ..assumptions import Assumptions, BEAM_NA, REGIME_NA, SPECTRUM_NA

DEFAULT_TAU_ZENITH = 0.05  # see module docstring


def airmass(elevation_deg):
    '''
    Airmass along the slant path, 1/sin(elevation).

    Plane-parallel secant model (borrowed from fso_spot_size.airmass). Diverges
    at the horizon.

    Parameters:
        elevation_deg : float or numpy.ndarray
            Elevation above the horizon [deg], > 0.

    Returns:
        float or numpy.ndarray
            Airmass (1 at zenith).
    '''
    return 1 / np.sin(np.radians(elevation_deg))


def atmospheric_transmittance(elevation_deg, tau_zenith=DEFAULT_TAU_ZENITH):
    '''
    Clear-sky transmittance along the slant path (Beer-Lambert x airmass).

    Parameters:
        elevation_deg : float or numpy.ndarray
            Elevation above the horizon [deg], > 0.
        tau_zenith : float
            Zenith optical depth (molecular + aerosol).

    Returns:
        float or numpy.ndarray
            Transmittance in (0, 1].
    '''
    return np.exp(-tau_zenith * airmass(elevation_deg))


def atmospheric_loss_term(scenario, geometry, tau_zenith=DEFAULT_TAU_ZENITH):
    '''
    Deterministic clear-sky atmospheric extinction as a link-budget Term.

    Parameters:
        scenario : SpaceScenario
            Unused for the model itself (single-parameter tau); accepted for the
            uniform model signature.
        geometry : object
            Anything exposing .elevation_deg (float or ndarray).
        tau_zenith : float
            Zenith optical depth (molecular + aerosol).

    Returns:
        Term
            Deterministic Term (loss positive dB), mean_db broadcasting over the
            geometry elevation shape.
    '''
    T = atmospheric_transmittance(geometry.elevation_deg, tau_zenith)
    loss_db = -10 * np.log10(T)
    a = Assumptions(
        beam_type=BEAM_NA,
        turbulence_regime=REGIME_NA,
        spectrum=SPECTRUM_NA,
        validity="Clear sky. One zenith optical depth. Plane-parallel airmass. "
                 "Elevation above 5 deg.",
    )
    if np.any(np.asarray(geometry.elevation_deg) < 5.0):
        a.flag("Elevation below 5 deg breaks the plane-parallel airmass model.")
    return Term(
        name="atmospheric extinction",
        category="atmospheric",
        mean_db=loss_db,
        note=f"clear-sky Beer-Lambert x airmass, tau_zenith={tau_zenith}",
        assumptions=a,
    )


def horizontal_extinction_db(path_length_m, attenuation_db_per_km):
    '''
    Beer-Lambert extinction over a horizontal (terrestrial) path.

    A horizontal path has a constant extinction per unit length, so the loss is
    linear in the path length. The coefficient is quoted directly in dB/km, so
    no logarithm is needed:

        loss_db = attenuation_db_per_km * (path_length_m / 1000)

    The coefficient value is weather- and visibility-dependent (fog, haze,
    rain). The user sets it per site; this function only scales it by the path.

    Parameters:
        path_length_m : float or numpy.ndarray
            Horizontal path length L [m].
        attenuation_db_per_km : float
            Clear-air / haze extinction coefficient [dB/km].

    Returns:
        float or numpy.ndarray
            Extinction loss [dB], positive.
    '''
    return attenuation_db_per_km * (np.asarray(path_length_m, dtype=float) / 1e3)


def terrestrial_extinction_term(scenario, geometry):
    '''
    Horizontal-path atmospheric extinction as a link-budget Term.

    The terrestrial counterpart of atmospheric_loss_term. It reads the
    attenuation coefficient from the scenario TerrestrialChannel and the path
    length from the geometry. Deterministic (loss positive dB).

    Parameters:
        scenario : TerrestrialScenario
            Reads channel.attenuation_db_per_km (a TerrestrialChannel).
        geometry : object
            Provides slant_range_m (the horizontal path length).

    Returns:
        Term
            Deterministic Term, mean_db broadcasting over the path-length shape.
    '''
    gamma = scenario.channel.attenuation_db_per_km
    loss_db = horizontal_extinction_db(geometry.slant_range_m, gamma)
    return Term(
        name="atmospheric extinction (horizontal)",
        category="atmospheric",
        mean_db=loss_db,
        note=f"horizontal Beer-Lambert, {gamma} dB/km",
        assumptions=Assumptions(
            beam_type=BEAM_NA,
            turbulence_regime=REGIME_NA,
            spectrum=SPECTRUM_NA,
            validity="Horizontal path. One dB-per-km extinction coefficient, "
                     "constant along the path. The coefficient value is the "
                     "user's to source (weather / visibility dependent).",
        ),
    )


if __name__ == '__main__':
    # zenith airmass ~ 1, and airmass rises as elevation drops
    assert abs(airmass(90.0) - 1.0) < 1e-12
    assert airmass(30.0) > airmass(60.0) > airmass(90.0)

    # transmittance in (0, 1], deterministic
    el = np.array([10.0, 30.0, 60.0, 90.0])
    T = atmospheric_transmittance(el)
    assert np.all(T > 0) and np.all(T <= 1)

    # loss positive and larger at low elevation (thicker slant path)
    class _G:
        elevation_deg = el
    term = atmospheric_loss_term(None, _G())
    loss = term.mean_db
    assert np.all(loss > 0)
    assert loss[0] > loss[-1]                       # 10 deg > 90 deg
    assert term.sampler is None and term.quantile is None  # deterministic
    assert term.assumptions is not None

    # A 3 deg elevation breaks the plane-parallel airmass model; 45 deg does not.
    class _G3:
        elevation_deg = 3.0
    class _G45:
        elevation_deg = 45.0
    assert not atmospheric_loss_term(None, _G3()).assumptions.ok
    assert atmospheric_loss_term(None, _G45()).assumptions.ok

    # loss at 30 deg > loss at zenith
    class _G30:
        elevation_deg = 30.0
    class _G90:
        elevation_deg = 90.0
    assert (atmospheric_loss_term(None, _G30()).mean_db
            > atmospheric_loss_term(None, _G90()).mean_db)

    # Horizontal extinction: linear in path length, positive loss.
    assert horizontal_extinction_db(1e3, 0.5) == 0.5          # 1 km at 0.5 dB/km
    assert horizontal_extinction_db(4e3, 0.5) == 2.0          # 4 km -> 2 dB
    assert horizontal_extinction_db(0.0, 0.5) == 0.0          # zero path -> no loss

    from ..scenario import TerrestrialScenario, TerrestrialChannel
    from ..terminal import Terminal
    from ..geometry import HorizontalPath
    terr = TerrestrialScenario(near=Terminal(aperture_m=0.1), far=Terminal(aperture_m=0.1),
                               channel=TerrestrialChannel(attenuation_db_per_km=0.5))
    hterm = terrestrial_extinction_term(terr, HorizontalPath(2e3))
    assert hterm.category == "atmospheric" and hterm.sampler is None
    assert np.isclose(hterm.mean_db, 1.0)                     # 2 km * 0.5 dB/km
    assert hterm.assumptions is not None

    print(f"tau_zenith = {DEFAULT_TAU_ZENITH}")
    for e, t, l in zip(el, T, loss):
        print(f"  elev {e:5.1f} deg  airmass {airmass(e):5.2f}  T {t:.4f}  loss {l:.3f} dB")
    print(f"  horizontal 2 km @ 0.5 dB/km -> {float(hterm.mean_db):.3f} dB")
    print("self-check OK")
