'''
Residual-wavefront physics for a downlink receive terminal.

The satellite is far away. The downlink source is a plane wave at the top of the
atmosphere. This module gives the plane-wave Fried parameter r0 for that plane
wave, and the residual phase variance that a wavefront-compensation stack leaves.
The functions are pure. They take numeric inputs and a compensation stack of
pure-data stages, and they return numbers or a ResidualWavefront. The Term
factory lives in olb.models.coupling.

Fried parameter (plane wave, downlink):
    r0 = ( 0.423 * k^2 * airmass * integral[ Cn2(h) dh ] )^(-3/5)
    Source: Fried (1966); Andrews and Phillips, Laser Beam Propagation through
    Random Media, 2nd ed. (2005), Ch. 3. Here k = 2*pi/lambda, airmass =
    1/sin(elevation), and the integral runs the zenith Cn2 profile over height.
    This is the PLANE-wave constant 0.423. The spherical-wave coherence diameter
    (uplink) uses a different weight and lives in olb.turbulence.coupled_flux.

Residual phase variance (Noll 1976):
    Over an aperture of diameter D and Fried parameter r0, the residual phase
    variance after Zernike correction is
        sigma^2 = c * (D/r0)^(5/3)     [rad^2]
    with the Noll coefficient c:
        piston removed only          c = 1.0299    (no correction)
        first 3 Zernikes removed     c = 0.134     (tip-tilt correction)
        first J Zernikes removed     c = 0.2944 * J^(-sqrt(3)/2)   (large-J AO)
    Source: R. J. Noll, "Zernike polynomials and atmospheric turbulence,"
    J. Opt. Soc. Am. 66(3), 207-211 (1976). The AO form is the large-order
    asymptotic of Noll's residual series.
'''

from dataclasses import dataclass, field
from typing import Callable

import numpy as np

# Noll residual coefficients c for sigma^2 = c * (D/r0)^(5/3). Source: Noll 1976.
NOLL_PISTON = 1.0299          # piston removed only (no correction)
NOLL_TIPTILT = 0.134          # first 3 Zernikes removed (tip-tilt)
_NOLL_AO_CONST = 0.2944       # large-J asymptotic prefactor
_AO_EXP = -np.sqrt(3.0) / 2.0  # large-J asymptotic exponent


def plane_wave_fried_parameter(cn2_profile, hs, wavelength, elevation_deg):
    '''
    Return the plane-wave Fried parameter r0 for the downlink.

    Integrate the zenith Cn2 profile and scale it to the slant path with the
    airmass. See the module docstring for the formula and the citation.

    Parameters:
        cn2_profile : numpy.ndarray
            Zenith Cn2(h) profile on the hs grid [m^-2/3].
        hs : numpy.ndarray
            Heights above the ground station [m].
        wavelength : float
            Optical wavelength [m].
        elevation_deg : float or numpy.ndarray
            Elevation angle above the horizon [deg].

    Returns:
        float or numpy.ndarray
            r0 [m], broadcast over the elevation shape.
    '''
    k = 2.0 * np.pi / wavelength
    airmass = 1.0 / np.sin(np.radians(np.asarray(elevation_deg, dtype=float)))
    integral = np.trapz(np.asarray(cn2_profile, dtype=float), hs)
    r0 = (0.423 * k ** 2 * airmass * integral) ** (-3.0 / 5.0)
    if np.ndim(elevation_deg) == 0:
        return float(r0)
    return r0


def _noll_coefficient_and_modes(stack):
    '''
    Return the residual Noll coefficient and the corrected-mode count.

    The best-correcting stage wins. An empty stack leaves the piston-removed
    turbulence. A TipTilt stage removes the first 3 Zernikes. An AO(n) stage
    removes the first n Zernikes with the large-J asymptotic residual.

    Returns:
        tuple
            (c, n_modes) : the Noll coefficient and the number of removed modes.
    '''
    # Import here to keep this pure-physics module free of a top-level olb import.
    from ..terminal import TipTilt, AO

    c = NOLL_PISTON
    n_modes = 0
    for stage in stack:
        if isinstance(stage, TipTilt):
            c = min(c, NOLL_TIPTILT)
            n_modes = max(n_modes, 3)
        elif isinstance(stage, AO):
            c = min(c, _NOLL_AO_CONST * stage.n_modes ** _AO_EXP)
            n_modes = max(n_modes, int(stage.n_modes))
        else:
            raise TypeError(
                f"unknown compensation stage {type(stage).__name__!r}. Use "
                "TipTilt or AO."
            )
    return c, n_modes


@dataclass
class ResidualWavefront:
    '''
    The residual wavefront that a compensation stack leaves (the fidelity
    contract).

    It carries two faces of the same residual:
        variance : the residual phase variance sigma^2 [rad^2]. Fidelity 0 reads
            this value.
        psd : a callable f -> residual radial phase power-spectral density
            [rad^2 m^2] against spatial frequency f [1/m]. This is the FIDELITY-1
            drop-in. It is the high-pass Kolmogorov phase spectrum above the AO
            correction cutoff. Fidelity 0 does NOT evaluate it. See _residual_psd.
    '''
    variance: float                 # residual phase variance [rad^2]
    psd: Callable                   # fidelity-1 hook: radial phase PSD [rad^2 m^2]
    r0: float                       # Fried parameter used [m]
    D: float                        # aperture diameter used [m]
    n_modes: int = 0                # corrected Zernike modes
    coefficient: float = NOLL_PISTON  # Noll residual coefficient
    meta: dict = field(default_factory=dict)


def _residual_psd(r0, D, n_modes):
    '''
    Build the fidelity-1 residual-phase PSD callable.

    Return the high-pass Kolmogorov phase power-spectral density above the AO
    correction cutoff. Below the cutoff the AO corrects the phase, so the
    residual is zero. Above the cutoff the residual is the full Kolmogorov
    spectrum. Fidelity 0 does NOT call this. It is the documented drop-in for the
    fidelity-1 model.

    formula:
        Phi_phi(f) = 0.023 * r0^(-5/3) * f^(-11/3)   for f > f_c, else 0
        f_c = sqrt(n_modes) / (2 D)                  correction cutoff [1/m]
    Source: Kolmogorov phase PSD, Andrews and Phillips, 2nd ed. (2005), Ch. 3.
    An AO stage that corrects n_modes modes flattens the phase up to the cutoff
    spatial frequency f_c. Above f_c the residual is uncorrected.
    '''
    f_c = np.sqrt(max(n_modes, 0)) / (2.0 * D)

    def psd(f):
        f = np.asarray(f, dtype=float)
        out = np.where(f > f_c, 0.023 * r0 ** (-5.0 / 3.0) * f ** (-11.0 / 3.0),
                       0.0)
        return out

    return psd


def apply_compensation(stack, D, r0):
    '''
    Apply a compensation stack to the raw turbulence over one aperture.

    Return the ResidualWavefront that the stack leaves. The best-correcting stage
    sets the residual. See _noll_coefficient_and_modes and the module docstring.

    Parameters:
        stack : list
            The ordered compensation stack (TipTilt, AO). It may be empty.
        D : float
            Aperture diameter [m].
        r0 : float
            Plane-wave Fried parameter [m].

    Returns:
        ResidualWavefront
            variance = c * (D/r0)^(5/3), plus the fidelity-1 PSD hook.
    '''
    c, n_modes = _noll_coefficient_and_modes(stack)
    variance = c * (D / np.asarray(r0, dtype=float)) ** (5.0 / 3.0)
    scalar = np.ndim(r0) == 0
    return ResidualWavefront(
        # Keep an array when r0 is an array (an elevation sweep).
        variance=float(variance) if scalar else variance,
        psd=_residual_psd(r0, D, n_modes),
        r0=float(r0) if scalar else np.asarray(r0, dtype=float),
        D=float(D),
        n_modes=int(n_modes),
        coefficient=float(c),
    )


if __name__ == '__main__':
    from .profiles import DEFAULT_HS, get_c2n
    from ..terminal import TipTilt, AO

    lam = 1550e-9
    hs = DEFAULT_HS
    cn2 = get_c2n(hs, 21.0, 1.7e-14)

    # r0 rises with wavelength (r0 ~ lambda^(6/5)).
    r0_1064 = plane_wave_fried_parameter(cn2, hs, 1064e-9, 60.0)
    r0_1550 = plane_wave_fried_parameter(cn2, hs, 1550e-9, 60.0)
    assert r0_1550 > r0_1064, (r0_1550, r0_1064)

    # r0 falls toward the horizon (longer slant path, more turbulence).
    r0_30 = plane_wave_fried_parameter(cn2, hs, lam, 30.0)
    r0_90 = plane_wave_fried_parameter(cn2, hs, lam, 90.0)
    assert r0_30 < r0_90, (r0_30, r0_90)

    D = 0.7
    r0 = r0_60 = plane_wave_fried_parameter(cn2, hs, lam, 60.0)

    # Residual variance falls as more modes are corrected.
    v_none = apply_compensation([], D, r0).variance
    v_tt = apply_compensation([TipTilt()], D, r0).variance
    v_ao20 = apply_compensation([TipTilt(), AO(20)], D, r0).variance
    v_ao200 = apply_compensation([AO(200)], D, r0).variance
    assert v_ao200 < v_ao20 < v_tt < v_none, (v_ao200, v_ao20, v_tt, v_none)

    # The Noll coefficients are the documented values.
    assert np.isclose(apply_compensation([], D, r0).coefficient, NOLL_PISTON)
    assert np.isclose(apply_compensation([TipTilt()], D, r0).coefficient, NOLL_TIPTILT)

    # AO(large N) << tip-tilt-only << uncorrected.
    assert v_ao200 < 0.1 * v_tt

    # The fidelity-1 PSD hook is present and high-passes above the cutoff.
    rw = apply_compensation([AO(60)], D, r0)
    f_c = np.sqrt(60) / (2.0 * D)
    assert rw.psd(0.5 * f_c) == 0.0            # below cutoff: corrected
    assert rw.psd(2.0 * f_c) > 0.0             # above cutoff: residual

    # The elevation array broadcasts.
    r0_sweep = plane_wave_fried_parameter(cn2, hs, lam, np.array([30.0, 60.0, 90.0]))
    assert r0_sweep.shape == (3,)

    print(f"r0(1064 nm) = {r0_1064*100:.2f} cm   r0(1550 nm) = {r0_1550*100:.2f} cm")
    print(f"r0 @60deg   = {r0_60*100:.2f} cm     D/r0 = {D/r0_60:.2f}")
    print(f"residual sigma^2 [rad^2]: none={v_none:.2f}  tip-tilt={v_tt:.3f}  "
          f"AO20={v_ao20:.4f}  AO200={v_ao200:.5f}")
    print("self-check passed")
