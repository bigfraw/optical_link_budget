'''
Transmit Gaussian efficiency: the truncation loss at the launch aperture.

A real transmitter sends a Gaussian beam through a finite, possibly obscured,
circular aperture. The aperture clips the wings of the Gaussian and the central
obscuration blocks the middle. This is a fixed hardware loss. It does not depend
on range.

The efficiency is the on-axis far-field gain of the truncated beam, referenced
to the gain of the untruncated source Gaussian:

    alpha = a / w_T = (tx_aperture_m / 2) / tx_waist_m
    eta   = (exp(-alpha^2) - exp(-alpha^2 * Cr^2))^2

where a is the aperture radius, w_T is the Gaussian waist (1/e^2 radius) at the
aperture, and Cr is the linear central-obscuration ratio (obscuration diameter /
aperture diameter). The loss is -10*log10(eta), positive dB.

eta goes to 1 when the aperture is much wider than the beam (no truncation) and
to 0 when the aperture is much smaller than the beam (the wings carry the power).
The classic optimal truncation is near alpha = 1.12.

This is the corrected antenna-gain form. It has NO 2/alpha^2 prefactor; that
prefactor double-counts a normalisation the untruncated-source reference already
carries. The removal is validated against a numerical Fraunhofer propagation of a
truncated Gaussian (tn2_kepler test_gauss_prop).

A DIVERGED launch (Transmitter.divergence_rad, a virtual waist behind the
aperture; see olb.beam.virtual_waist) reaches the aperture with a CURVED
phase front of radius R (olb.beam.launch_curvature). The on-axis far field is
then the integral of exp(-r^2 (1/w_T^2 + i k / (2 R))) over the annulus, so
alpha^2 takes an imaginary part, the defocus parameter of Klein and Degnan
(DOI 10.1364/AO.13.002134):

    beta = k a^2 / (2 R)
    eta  = |exp(-(alpha^2 + i beta)) - exp(-(alpha^2 + i beta) Cr^2)|^2

The reference is the UNTRUNCATED curved beam, which is the beam of the
geometric spreading Term (a Gaussian from the virtual waist). So the two Terms
multiply to the exact on-axis far field of the clipped, diverged launch. A
collimated beam has beta = 0 and the form above. Before 2026-09-24 the Term
took beta = 0 at every divergence, which over-read a diverged link by 0.7 to
2.7 dB (validation/uplink_divergence, verdict 7; backlog 2-DV item 8).

The transmitter is the ground station for an uplink, or the satellite for a
downlink. The same code serves both.

FINITE RANGE. The eta above is a FAR-FIELD quantity: the on-axis gain ratio
far past the aperture. The hard aperture edge changes the received pattern
with range before the far field. Also, a receiver that is not small against
the beam does not read the on-axis value. So on a short horizontal link the
far-field eta can be wrong by more than 1 dB, in EITHER direction. Example: a
36.3 mm (1/e^2 diameter) launch in a 50.8 mm aperture, a 50.8 mm receiver at
1 km. The far-field eta gives 1.32 dB, and the true loss is 0.18 dB.

When the Term gets a geometry, it uses the EXACT paraxial value at that range.
The Fresnel diffraction integral of a circular field is a Hankel transform
(Schmidt, DOI 10.1117/3.866274, Ch. 6; Andrews and Phillips,
DOI 10.1117/3.626196, Ch. 4, the Huygens-Fresnel integral). At the range L
the received field of the clipped launch is

    U(rho) = (k / L) int_{Cr a}^{a} exp(-c r^2) J0(k r rho / L) r dr
    c      = 1 / w_T^2 + i (k / 2) (1 / f0 - 1 / L)

f0 is the launch phase-front radius (olb.beam.launch_curvature, negative for
a diverged beam). The 1 / L part is the finite-range defocus of Klein and
Degnan (DOI 10.1364/AO.13.002134). The power in the receive annulus is
P = 2 pi int |U|^2 rho drho. The Term divides P by the received power of the
UNTRUNCATED beam, which is the geometric spreading Term (olb.models.geometric).
So the two Terms multiply to the exact received power of the clipped launch
at every range. For a small receiver far away, the ratio goes to the
far-field eta above. geometry=None keeps the far-field eta. The self-check
compares the quadrature with olb.waveoptics.run.propagate_scenario.
'''

import numpy as np
from scipy.special import j0

from ..beam import launch_curvature
from ..results import Term
from ..assumptions import Assumptions, BEAM_GAUSSIAN, REGIME_NA, SPECTRUM_NA
from .geometric import geometric_loss_db


def gaussian_efficiency(alpha, obscuration_ratio=0.0, beta=0.0):
    '''
    On-axis truncation efficiency of a Gaussian beam in a circular aperture.

    Parameters:
        alpha : float or ndarray
            Ratio of aperture radius to Gaussian waist, a / w_T.
        obscuration_ratio : float
            Linear central-obscuration ratio Cr (obscuration diameter / aperture
            diameter). 0 = unobscured.
        beta : float
            The defocus parameter k a^2 / (2 R) of a curved launch phase front
            of radius R (Klein and Degnan, DOI 10.1364/AO.13.002134). 0 = a
            flat phase front (a collimated launch).

    Returns:
        float or ndarray
            Efficiency eta in (0, 1].
    '''
    a2 = np.asarray(alpha, dtype=float) ** 2 + 1j * beta
    return np.abs(np.exp(-a2) - np.exp(-a2 * obscuration_ratio ** 2)) ** 2


def launch_defocus_beta(tx_aperture_m, tx_waist_m, divergence_rad=None,
                        wavelength=1550e-9):
    '''
    The defocus parameter beta = k a^2 / (2 |R|) of a diverged launch.

    R is the phase-front radius at the aperture, olb.beam.launch_curvature
    (Andrews and Phillips, DOI 10.1117/3.626196, Ch. 4, Eqs. (7) and (8)). The
    sign of R does not change eta, so the magnitude serves. A collimated launch
    gives 0.0.
    '''
    f0 = launch_curvature(tx_waist_m, divergence_rad, wavelength)
    if not np.isfinite(f0):
        return 0.0
    k = 2 * np.pi / wavelength
    return float(k * (tx_aperture_m / 2) ** 2 / (2 * abs(f0)))


def tx_efficiency_loss_db(tx_aperture_m, tx_waist_m, obscuration_ratio=0.0,
                          divergence_rad=None, wavelength=1550e-9):
    '''
    Transmit truncation loss of a Gaussian beam at a circular aperture.

    Parameters:
        tx_aperture_m : float
            Transmit aperture diameter [m].
        tx_waist_m : float
            Gaussian waist (1/e^2 radius) at the aperture [m].
        obscuration_ratio : float
            Linear central-obscuration ratio Cr. 0 = unobscured.
        divergence_rad : float, optional
            Far-field HALF-angle divergence [rad]. None = collimated.
        wavelength : float
            Wavelength [m]. Read only for a diverged launch.

    Returns:
        float
            Truncation loss [dB], positive.
    '''
    alpha = (tx_aperture_m / 2) / tx_waist_m
    beta = launch_defocus_beta(tx_aperture_m, tx_waist_m, divergence_rad,
                               wavelength)
    eta = gaussian_efficiency(alpha, obscuration_ratio, beta)
    return -10 * np.log10(eta)


def clipped_bucket_fraction(range_m, tx_aperture_m, tx_waist_m, rx_aperture_m,
                            obscuration_ratio=0.0, rx_obscuration_ratio=0.0,
                            divergence_rad=None, wavelength=1550e-9):
    '''
    Exact paraxial power of a clipped Gaussian launch in a receive annulus.

    The function integrates the Fresnel (Hankel) form of FINITE RANGE in the
    module docstring with the trapezoid rule on both radii. It is valid at
    every range, in the near field and in the far field.

    Parameters:
        range_m : float
            The range L from the launch aperture to the receiver [m].
        tx_aperture_m, tx_waist_m : float
            The launch aperture diameter and the Gaussian 1/e^2 radius [m].
        rx_aperture_m : float
            The receive aperture diameter [m].
        obscuration_ratio, rx_obscuration_ratio : float
            The launch and the receive linear obscuration ratios.
        divergence_rad : float, optional
            The launch far-field HALF-angle divergence [rad]. None = collimated.
        wavelength : float
            Wavelength [m].

    Returns:
        float
            The received power as a fraction of the power of the UNTRUNCATED
            launch Gaussian (pi w_T^2 / 2 at a unit peak amplitude).
    '''
    k = 2 * np.pi / wavelength
    a, a_rx = tx_aperture_m / 2, rx_aperture_m / 2
    f0 = launch_curvature(tx_waist_m, divergence_rad, wavelength)
    c = 1 / tx_waist_m ** 2 + 0.5j * k * (1 / f0 - 1 / range_m)
    # The sample count follows the largest phase across the two discs (the
    # quadratic phase and the J0 argument). 400 samples hold a space link.
    phase = abs(c.imag) * a ** 2 + k * a * a_rx / range_m
    n = int(max(400, 6 * phase))
    r = np.linspace(obscuration_ratio * a, a, n)
    rho = np.linspace(rx_obscuration_ratio * a_rx, a_rx, n)
    u = (k / range_m) * np.trapezoid(
        np.exp(-c * r ** 2) * r * j0(k * np.outer(rho, r) / range_m), r, axis=1)
    p = 2 * np.pi * np.trapezoid(np.abs(u) ** 2 * rho, rho)
    return float(p / (np.pi * tx_waist_m ** 2 / 2))


def uniform_aperture_correction_db(obscuration_ratio=0.0):
    '''
    Correction that turns the Gaussian(w0=aperture/2) far field into a top-hat.

    The geometric spreading Term models the transmitter as an UNtruncated
    Gaussian of waist w0 = aperture/2. That Gaussian carries power in wings past
    the aperture, so its on-axis far-field gain is a factor of 2 / (1 - Cr^2)
    higher than a UNIFORMLY illuminated (top-hat) aperture of the same diameter.

    A retroreflector reflects the roughly flat wavefront that fills its aperture,
    so its return is a top-hat, not a Gaussian. This correction converts the
    Gaussian(w0=aperture/2) geometric model to the top-hat.

    Derivation. The ratio is eta / tau, where eta = (exp(-a^2) - exp(-a^2*Cr^2))^2
    is the truncation efficiency (see gaussian_efficiency) and
    tau = 2*eta / (a^2*(1-Cr^2)) is the aperture-illumination efficiency (peak
    gain of the truncated Gaussian relative to a uniform aperture). The ratio is
    eta / tau = a^2*(1-Cr^2)/2, which at the Gaussian(w0=aperture/2) reference
    (a = alpha = 1) is (1-Cr^2)/2. So the top-hat has 2/(1-Cr^2) times LESS
    on-axis gain, a positive-dB loss.

    Parameters:
        obscuration_ratio : float
            Linear central-obscuration ratio Cr of the transmit aperture. A
            corner-cube retro has no obscuration (Cr = 0), so the correction is
            +10*log10(2) = 3.01 dB.

    Returns:
        float
            Top-hat correction [dB], positive.
    '''
    return 10.0 * np.log10(2.0 / (1.0 - obscuration_ratio ** 2))


def tx_gaussian_efficiency_term(scenario, geometry=None):
    '''
    Transmit Gaussian-efficiency (truncation) Term for a scenario.

    With a geometry, the loss is the EXACT paraxial value at the range: the
    received power of the clipped launch (clipped_bucket_fraction) over the
    received power of the untruncated beam (the geometric spreading Term). So
    it reads the receive aperture too, and it holds in the near field. With
    geometry=None, the loss is the far-field eta (gaussian_efficiency). See
    FINITE RANGE in the module docstring.

    Parameters:
        scenario : SpaceScenario or TerrestrialScenario
            Reads the transmit terminal Transmitter waist_m and the launch
            aperture. The launch aperture is the Transmitter aperture_m and
            obscuration_ratio when set (a bistatic beam director), else the owning
            Terminal aperture_m and obscuration_ratio (monostatic). See
            olb.terminal.Transmitter and olb.scenario.
        geometry : object, optional
            Read for slant_range_m (float or ndarray). None gives the
            far-field eta.

    Returns:
        Term
            name="transmit Gaussian efficiency", category="system".
    '''
    tx = scenario.tx_terminal
    t = tx.transmitter
    waist_m = t.waist_m
    # Bistatic override: a Transmitter may carry its own beam-director aperture
    # and obscuration. When either is None, fall back to the owning Terminal
    # value (the monostatic default). See olb.terminal.Transmitter.
    aperture_m = t.aperture_m if t.aperture_m is not None else tx.aperture_m
    obscuration_ratio = (t.obscuration_ratio if t.obscuration_ratio is not None
                         else tx.obscuration_ratio)
    alpha = (aperture_m / 2) / waist_m
    beta = launch_defocus_beta(aperture_m, waist_m, t.divergence_rad,
                               tx.wavelength_m)
    range_m = getattr(geometry, "slant_range_m", None)
    if range_m is None:
        loss = tx_efficiency_loss_db(aperture_m, waist_m, obscuration_ratio,
                                     t.divergence_rad, tx.wavelength_m)
    else:
        # The exact received power of the clipped launch over that of the
        # untruncated beam (the geometric spreading Term): see FINITE RANGE.
        rx = scenario.rx_terminal
        ranges = np.asarray(range_m, dtype=float)
        exact = np.vectorize(lambda z: clipped_bucket_fraction(
            z, aperture_m, waist_m, rx.aperture_m, obscuration_ratio,
            rx.obscuration_ratio, t.divergence_rad, tx.wavelength_m))(ranges)
        untruncated_db = geometric_loss_db(
            ranges, waist_m, rx.aperture_m, tx.wavelength_m,
            rx.obscuration_ratio, t.divergence_rad)
        loss = -10 * np.log10(exact) - untruncated_db
        loss = float(loss) if loss.ndim == 0 else loss
    assumptions = Assumptions(
        beam_type=BEAM_GAUSSIAN,
        turbulence_regime=REGIME_NA,
        spectrum=SPECTRUM_NA,
        validity="Truncated Gaussian launch (with the launch phase-front "
                 "curvature of a diverged beam), referenced to the untruncated "
                 "beam. Paraxial. No turbulence. With a geometry: the exact "
                 "Fresnel received power at the range, near and far field. "
                 "Without a geometry: the on-axis FAR-FIELD gain only.",
    )
    return Term(
        name="transmit Gaussian efficiency",
        category="system",
        mean_db=loss,
        note=f"aperture truncation, alpha={alpha:.3f}, beta={beta:.3f}, "
             f"Cr={obscuration_ratio:g}",
        meta={"alpha": float(alpha), "beta": float(beta),
              "eta": float(gaussian_efficiency(alpha, obscuration_ratio,
                                               beta))},
        assumptions=assumptions,
    )


if __name__ == '__main__':
    from ..scenario import SpaceScenario
    from ..terminal import Terminal, Transmitter

    # Limits: wide aperture -> no loss; narrow aperture -> large loss.
    assert gaussian_efficiency(10.0) > 1 - 1e-6              # aperture >> beam -> eta ~ 1
    assert tx_efficiency_loss_db(2.0, 0.1) < 1e-3           # aperture >> beam -> ~0 dB
    assert tx_efficiency_loss_db(0.01, 0.1) > 20            # aperture << beam -> big loss

    # Obscuration only adds loss.
    assert (tx_efficiency_loss_db(0.15, 0.12, obscuration_ratio=0.3)
            > tx_efficiency_loss_db(0.15, 0.12))

    # Corrected form: NO 2/alpha^2 prefactor. eta is bounded by 1, so the loss is
    # never a gain. The stale form (with 2/alpha^2) gives eta > 1 (negative loss)
    # for a small alpha; the corrected form does not.
    assert np.all(gaussian_efficiency(np.array([0.1, 0.3, 0.5, 1.0, 2.0])) <= 1.0)

    # Top-hat correction: unobscured -> +3.01 dB; obscuration adds more.
    assert abs(uniform_aperture_correction_db(0.0) - 3.0103) < 1e-3
    assert uniform_aperture_correction_db(0.3) > uniform_aperture_correction_db(0.0)

    # TN-2 launch: Da=150 mm, w_T=0.8*Da, Cr=0.3 -> alpha=0.625.
    tn2 = tx_efficiency_loss_db(0.150, 0.8 * 0.150, obscuration_ratio=0.3)
    assert 10.0 < tn2 < 12.0, tn2                            # ~10.8 dB

    # Term path: uplink -> tx=ground carries the aperture, obscuration, waist.
    scn = SpaceScenario(
        ground=Terminal(aperture_m=0.150, obscuration_ratio=0.3,
                        transmitter=Transmitter(waist_m=0.12)),
        space=Terminal(aperture_m=0.05),
        direction="uplink")
    term = tx_gaussian_efficiency_term(scn)
    assert term.category == "system"
    assert np.isscalar(term.mean_db) and term.mean_db > 0
    assert 0.0 < term.meta["eta"] <= 1.0
    assert term.assumptions is not None
    assert term.assumptions.ok

    # Far field: with a geometry the exact Term goes to the far-field eta. A
    # 600 km range and a 5 cm receiver agree to 0.02 dB (the launch Fresnel
    # number a^2 / (lambda L) is still 0.006 there).
    from types import SimpleNamespace
    far = tx_gaussian_efficiency_term(scn, SimpleNamespace(slant_range_m=600e3))
    assert abs(far.mean_db - term.mean_db) < 0.02, (far.mean_db, term.mean_db)
    # An array of ranges gives an array of losses.
    assert tx_gaussian_efficiency_term(
        scn, SimpleNamespace(slant_range_m=np.array([5e5, 6e5]))).mean_db.shape == (2,)

    # The denominator is the geometric spreading Term: an open launch aperture
    # (no clip) gives a ratio of 1, collimated and diverged.
    for div in (None, 2e-4):
        geo = 10 ** (-geometric_loss_db(1e3, 0.025, 0.05, 1550e-9, 0.0, div) / 10)
        ratio = clipped_bucket_fraction(1e3, 0.4, 0.025, 0.05,
                                        divergence_rad=div) / geo
        assert abs(ratio - 1) < 1e-3, (div, ratio)

    # NEAR FIELD against the wave optics: a 36.3 mm (1/e^2 diameter) launch in
    # a 50.8 mm aperture to a 50.8 mm receiver (zR = 668 m). At 1 km the
    # far-field eta reads about 1.2 dB high; the exact Term agrees with
    # propagate_scenario to 0.1 dB at 1 km and at 10 km.
    import warnings
    from ..scenario import TerrestrialScenario, TerrestrialChannel
    from ..geometry import HorizontalPath
    from ..models.geometric import geometric_loss_term
    from ..waveoptics.run import propagate_scenario
    end = Terminal(aperture_m=0.0508,
                   transmitter=Transmitter(waist_m=16 * 2.27e-3 / 2))
    for L in (1e3, 10e3):
        terr = TerrestrialScenario(
            near=end, far=end,
            channel=TerrestrialChannel(path_length_m=L,
                                       attenuation_db_per_km=0.0))
        g = HorizontalPath(L)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            wave = propagate_scenario(terr, g, precision="double")
        wave_db = wave.tx_truncation_db + wave.geometric_loss_db
        geo_db = geometric_loss_term(terr, g).mean_db
        exact_db = geo_db + tx_gaussian_efficiency_term(terr, g).mean_db
        ff_db = geo_db + tx_gaussian_efficiency_term(terr).mean_db
        assert abs(exact_db - wave_db) < 0.1, (L, exact_db, wave_db)
        # The far-field eta is 1 dB or more high at 1 km and good at 10 km.
        assert (ff_db - wave_db > 1.0 if L < 5e3
                else abs(ff_db - wave_db) < 0.1), (L, ff_db, wave_db)
        print(f"L = {L / 1e3:g} km: exact {exact_db:.3f} dB, far-field "
              f"{ff_db:.3f} dB, wave optics {wave_db:.3f} dB")
    # A HARD clip deep in the near field (examples/waveoptics/terrestrial_stages):
    # alpha = 0.625, zR = 29 km, 1 km to a 0.2 m receiver. The far-field eta
    # reads 8 dB high; the exact Term holds to 0.1 dB.
    from ..waveoptics.grid import GridSpec
    hard = TerrestrialScenario(
        near=Terminal(aperture_m=0.15, transmitter=Transmitter(waist_m=0.12)),
        far=Terminal(aperture_m=0.20),
        channel=TerrestrialChannel(path_length_m=1e3, attenuation_db_per_km=0.0))
    g = HorizontalPath(1e3)
    wave = propagate_scenario(hard, g, grid=GridSpec(1.0, 1024),
                              precision="double")
    wave_db = wave.tx_truncation_db + wave.geometric_loss_db
    exact_db = (geometric_loss_term(hard, g).mean_db
                + tx_gaussian_efficiency_term(hard, g).mean_db)
    assert abs(exact_db - wave_db) < 0.1, (exact_db, wave_db)
    print(f"hard clip, 1 km: exact {exact_db:.3f} dB, wave optics "
          f"{wave_db:.3f} dB")

    # Bistatic override: the launch truncation reads the Transmitter beam-director
    # aperture, NOT the (large) receive telescope aperture. A ground terminal with
    # a 0.7 m receive telescope but a 0.15 m transmit director must give the same
    # loss as a monostatic 0.15 m launch aperture.
    bistatic = SpaceScenario(
        ground=Terminal(aperture_m=0.7, obscuration_ratio=0.3,
                        transmitter=Transmitter(waist_m=0.12, aperture_m=0.15,
                                                obscuration_ratio=0.3)),
        space=Terminal(aperture_m=0.05),
        direction="uplink")
    monostatic = SpaceScenario(
        ground=Terminal(aperture_m=0.15, obscuration_ratio=0.3,
                        transmitter=Transmitter(waist_m=0.12)),
        space=Terminal(aperture_m=0.05),
        direction="uplink")
    assert abs(tx_gaussian_efficiency_term(bistatic).mean_db
               - tx_gaussian_efficiency_term(monostatic).mean_db) < 1e-9

    # A DIVERGED launch: beta > 0, and a collimated one keeps beta = 0 (the
    # old number, bit for bit). The exact check: the complex form equals a
    # direct radial quadrature of the curved, clipped Gaussian.
    assert launch_defocus_beta(0.15, 0.06) == 0.0
    assert tx_efficiency_loss_db(0.15, 0.06, 0.3) == tx_efficiency_loss_db(
        0.15, 0.06, 0.3, divergence_rad=None)
    a, w, cr, lam = 0.075, 0.06, 0.3, 1550e-9
    for div in (85e-6, 200e-6):
        f0 = launch_curvature(w, div, lam)
        r = np.linspace(cr * a, a, 200001)
        e = np.exp(-r ** 2 / w ** 2 - 1j * np.pi * r ** 2 / (lam * f0)) * r
        c = 1 / w ** 2 + 1j * np.pi / (lam * f0)     # untruncated: 1/(2c)
        quad = abs(np.trapezoid(e, r)) ** 2 * abs(2 * c) ** 2
        got = gaussian_efficiency(a / w, cr,
                                  launch_defocus_beta(2 * a, w, div, lam))
        assert abs(got / quad - 1) < 1e-6, (div, got, quad)
        print(f"diverged {div * 1e6:.0f} urad: truncation "
              f"{-10 * np.log10(got):.3f} dB (collimated "
              f"{tx_efficiency_loss_db(2 * a, w, cr):.3f} dB)")

    print(f"TN-2 transmit truncation loss: {tn2:.2f} dB")
    print(f"Term: {term.name}  {term.mean_db:.2f} dB  ({term.note})")
    print("self-check passed")
