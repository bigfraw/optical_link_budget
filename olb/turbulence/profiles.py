'''
Cn2 turbulence profiles and the wind profile for optical link budgets.

This module OWNS the two continuous atmosphere models that olb uses: the
Hufnagel-Valley Cn2(h) profile (`get_c2n`) and the Bufton wind profile
(`v_wind`). Both are closed-form functions of altitude. It also builds the
default zenith Cn2(h) profile from the site parameters (`default_cn2_profile`)
and holds the default turbulence altitude grid (`DEFAULT_HS`).

Every olb module reads these models from here, so the physics lives in one
place. The Andrews wrappers `andrews.paths.hufnagel_valley` and
`andrews.paths.bufton_wind` DELEGATE to these functions; they add the cited
book context, not new physics.

Source: Andrews and Phillips, Laser Beam Propagation through Random Media,
2nd ed. (2005), DOI 10.1117/3.626196, Ch. 12, printed p. 481.
'''

import numpy as np

from ..assumptions import Constraint, assumes

DEFAULT_HS = np.logspace(np.log10(1), np.log10(20e3), 20)   # turbulence altitude grid [m]

# The two atmosphere models each state one modelling assumption: the empirical
# form of the height profile. There is no numeric validity gate, so neither
# constraint carries a check. Source: Andrews and Phillips, 2nd ed. (2005),
# DOI 10.1117/3.626196, Ch. 12, printed p. 481.
_HV_MODEL = Constraint(
    "approximation",
    "The Cn2 height profile is the empirical Hufnagel-Valley model. The default "
    "pair w = 21 m/s and A = 1.7e-14 is the HV5/7 model (r0 = 5 cm, "
    "theta0 = 7 urad at 0.5 um).",
    "10.1117/3.626196", "Ch. 12, Eq. (1), printed p. 481")

_BUFTON_MODEL = Constraint(
    "approximation",
    "The wind height profile is the empirical Bufton model.",
    "10.1117/3.626196", "Ch. 12, Eq. (3), printed p. 481")


@assumes(_HV_MODEL)
def get_c2n(height, wind_rms=21, c2n_0=1.7e-14):
    '''
    Return the Hufnagel-Valley Cn2(h) profile [m^-2/3].

    Cn2 is the refractive-index structure constant. The three terms model the
    high-altitude wind layer, the tropopause, and the ground boundary layer.

    Parameters:
        height : float or numpy.ndarray
            Altitude above ground level [m].
        wind_rms : float
            The rms high-altitude (pseudo)wind speed w [m/s].
        c2n_0 : float
            The ground value A = Cn2(0) [m^-2/3].

    Returns:
        numpy.ndarray
            Cn2(h) [m^-2/3].

    formula:
        Cn2(h) = 0.00594 (w/27)^2 (1e-5 h)^10 exp(-h/1000)
                 + 2.7e-16 exp(-h/1500)
                 + A exp(-h/100)
    The default pair w = 21 m/s and A = 1.7e-14 is the H-V5/7 model, chosen so
    that r0 = 5 cm and theta0 = 7 urad at lambda = 0.5 um.
    Source: Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196, Ch. 12,
    Eq. (1), printed p. 481.
    '''
    c2ns = 0.00594 * np.power(wind_rms / 27, 2) * np.power(1e-5 * height, 10) \
        * np.exp(-height / 1000) \
        + 2.7e-16 * np.exp(-height / 1500) \
        + c2n_0 * np.exp(-height / 100)
    return c2ns


@assumes(_BUFTON_MODEL)
def v_wind(h, ws=1, Vg=10):
    '''
    Return the Bufton wind-speed profile V(h) [m/s].

    Parameters:
        h : float or numpy.ndarray
            Altitude above ground level [m].
        ws : float
            The slew rate of the satellite as seen from the ground [deg/s].
        Vg : float
            The ground wind speed [m/s].

    Returns:
        numpy.ndarray
            V(h) [m/s].

    formula:
        V(h) = deg2rad(ws) h + Vg + 30 exp(-((h - 9400)/4800)^2)
    Source: Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196, Ch. 12,
    Eq. (3), printed p. 481.
    '''
    vs = np.deg2rad(ws) * h + Vg + 30 * np.exp(-((h - 9400) / 4800) ** 2)
    return vs


def default_cn2_profile(site, hs=None):
    '''
    Build a default zenith Cn2 profile from the site parameters.

    This is the H-V5/7 model on the altitude grid, driven by the site wind and
    ground turbulence. Use it when the `fast` package is not available (the
    `fast` HV57 path fails without that package).

    Parameters:
        site : Site
            Provides wind_rms_m_s and cn2_ground.
        hs : numpy.ndarray, optional
            Turbulence altitude grid [m]. Defaults to DEFAULT_HS.

    Returns:
        numpy.ndarray
            Cn2(h) profile at zenith on the hs grid.
    '''
    hs = DEFAULT_HS if hs is None else hs
    return get_c2n(hs, site.wind_rms_m_s, site.cn2_ground)


if __name__ == '__main__':
    # The H-V5/7 default gives r0 = 5 cm and theta0 = 7 urad at 0.5 um. Check the
    # profile shape and the moments that fix those two numbers.
    hs = DEFAULT_HS
    cn2 = get_c2n(hs, 21.0, 1.7e-14)
    assert cn2.shape == hs.shape
    assert np.all(cn2 > 0)
    assert cn2[0] > cn2[-1]                          # turbulence falls with height
    # At the ground the wind term is 0, so Cn2(0) = A + the 2.7e-16 tropopause
    # term, which is about 1.6 % above A for the H-V5/7 pair.
    assert abs(get_c2n(0.0, 21.0, 1.7e-14) / 1.7e-14 - 1.0) < 0.02

    # r0 of a zenith plane wave: r0 = (0.423 k^2 INT Cn2 dh)^(-3/5). The H-V5/7
    # model is set so that r0 is about 5 cm at 0.5 um. Source: Andrews and
    # Phillips, DOI 10.1117/3.626196, Ch. 12, printed p. 481.
    fine = np.linspace(0.0, 20e3, 20001)
    k = 2 * np.pi / 0.5e-6
    mu0 = np.trapezoid(get_c2n(fine, 21.0, 1.7e-14), fine)
    r0 = (0.423 * k ** 2 * mu0) ** (-3 / 5)
    assert abs(r0 / 0.05 - 1.0) < 0.15, r0          # about 5 cm

    # The Bufton wind: the ground value is Vg plus the tail of the high-altitude
    # bump (about 0.65 m/s), and the peak Vg + 30 sits at 9.4 km.
    assert abs(v_wind(0.0, 0.0, 10.0) - 10.0) < 1.0
    assert abs(v_wind(9400.0, 0.0, 10.0) - 40.0) < 1e-6      # Vg + 30, exp = 1
    assert v_wind(9400.0, 0.0, 10.0) > v_wind(0.0, 0.0, 10.0)

    # --- assumptions layer ---------------------------------------------------
    import warnings

    from ..assumptions import trace_assumptions

    # (1) Value parity: a decorated function returns the identical value with and
    #     without a collection context.
    outside = float(get_c2n(1000.0, 21.0, 1.7e-14))
    with trace_assumptions():
        inside = float(get_c2n(1000.0, 21.0, 1.7e-14))
    assert outside == inside, (outside, inside)

    # (2) Registration: inside a context the two model functions register, and
    #     the physics layer emits no warning. (profiles has no numeric gate, so
    #     there is no violation block.)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        with trace_assumptions() as trace:
            get_c2n(DEFAULT_HS, 21.0, 1.7e-14)
            v_wind(DEFAULT_HS, 0.0, 10.0)
            # default_cn2_profile is a pure delegator to get_c2n, so it carries
            # no decorator; the traced get_c2n call above holds its record.
    src_c2n = f"{__name__}.get_c2n"
    src_wind = f"{__name__}.v_wind"
    assert src_c2n in trace.records, trace.records
    assert src_wind in trace.records, trace.records
    kinds = {c.kind for rec in trace.records.values() for c in rec.constraints}
    assert kinds == {"approximation"}, kinds
    assert len(caught) == 0, "the profiles physics must not warn"

    print(f"H-V5/7 self-check: r0(0.5 um, zenith) = {r0 * 100:.2f} cm "
          f"(book about 5 cm); Cn2(0) = {cn2[0]:.2e}")
    print("profiles self-check passed")
