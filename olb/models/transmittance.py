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
        scenario : Scenario
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

    print(f"tau_zenith = {DEFAULT_TAU_ZENITH}")
    for e, t, l in zip(el, T, loss):
        print(f"  elev {e:5.1f} deg  airmass {airmass(e):5.2f}  T {t:.4f}  loss {l:.3f} dB")
    print("self-check OK")
