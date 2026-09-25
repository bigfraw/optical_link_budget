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

NEAR-FIELD WARNING. eta is a FAR-FIELD quantity. It is the on-axis gain ratio
that the truncated beam settles to far past the aperture. A pure untruncated
Gaussian keeps a smooth on-axis intensity at every range, but the hard aperture
edge does not: inside the Rayleigh range zR = pi*w_T^2/lambda the on-axis
intensity OSCILLATES with range (Fresnel edge diffraction). So for a receiver
inside zR (the terrestrial near-field case; zR is about 5 km for a 5 cm beam at
1.55 um) the true on-axis intensity can sit above OR below this eta. Unlike the
geometric spreading Term (exact at all ranges through gaussz) this Term does NOT
self-correct, and unlike the single-mode-fibre eta_max the error is NOT
conservative. It bites only when the beam is HARD truncated (alpha near or below
1, so real power sits in the clipped wings) AND the receiver is inside zR. A
lightly clipped beam (alpha >~ 1.5) or a far-field target is safe. A fidelity-2
no-turbulence field propagation is the check for a link that this Term flags.
'''

import numpy as np

from ..beam import launch_curvature
from ..results import Term
from ..assumptions import Assumptions, BEAM_GAUSSIAN, REGIME_NA, SPECTRUM_NA

# The far-field eta is safe above this alpha: the aperture is wide enough that
# almost no power sits in the clipped wings, so the near-field edge ripple is
# negligible even inside the Rayleigh range. Below it, a receiver inside zR gets
# an active near-field violation flag. See the NEAR-FIELD WARNING above.
TRUNCATION_NEAR_FIELD_ALPHA = 1.5


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

    The loss is range-independent. The geometry is read only to test the
    near-field validity: a hard-truncated beam (alpha below
    TRUNCATION_NEAR_FIELD_ALPHA) with a receiver inside the Rayleigh range
    zR = pi*w_T^2/lambda breaks the far-field eta, so the Term flags it (see the
    NEAR-FIELD WARNING in the module docstring). A far-field target or a widely
    open aperture keeps the Term valid, and geometry=None skips the test.

    Parameters:
        scenario : SpaceScenario or TerrestrialScenario
            Reads the transmit terminal Transmitter waist_m and the launch
            aperture. The launch aperture is the Transmitter aperture_m and
            obscuration_ratio when set (a bistatic beam director), else the owning
            Terminal aperture_m and obscuration_ratio (monostatic). See
            olb.terminal.Transmitter and olb.scenario.
        geometry : object, optional
            Read for slant_range_m to test the near-field validity only. None
            skips the test (the loss is unchanged).

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
    loss = tx_efficiency_loss_db(aperture_m, waist_m, obscuration_ratio,
                                 t.divergence_rad, tx.wavelength_m)
    assumptions = Assumptions(
        beam_type=BEAM_GAUSSIAN,
        turbulence_regime=REGIME_NA,
        spectrum=SPECTRUM_NA,
        validity="On-axis FAR-FIELD gain of a truncated Gaussian (with the "
                 "launch phase-front curvature of a diverged beam), referenced "
                 "to the untruncated source. Paraxial. No turbulence. Fails "
                 "for a receiver inside the Rayleigh range zR=pi*w_T^2/lambda "
                 "when the beam is hard truncated (alpha<~1): near-field "
                 "on-axis intensity oscillates (Fresnel edge diffraction) and "
                 "this far-field eta does not self-correct. Verify such a link "
                 "with a fidelity-2 no-turbulence field propagation.",
    )
    # Active near-field violation: a hard-truncated beam (small alpha) with a
    # receiver inside the Rayleigh range breaks the far-field eta. The loss is
    # geometry-dependent to test, so it is an active flag, not just the prose
    # validity string. geometry=None (space budgets, always far field) skips it.
    range_m = getattr(geometry, "slant_range_m", None)
    if range_m is not None and alpha < TRUNCATION_NEAR_FIELD_ALPHA:
        rayleigh_range_m = np.pi * waist_m ** 2 / tx.wavelength_m
        if np.any(np.asarray(range_m, dtype=float) < rayleigh_range_m):
            assumptions.flag(
                f"Receiver is inside the Rayleigh range "
                f"(zR={rayleigh_range_m:.3g} m) and the beam is hard truncated "
                f"(alpha={alpha:.3f}): the far-field truncation eta does not "
                f"hold. Verify with a fidelity-2 no-turbulence field propagation."
            )
    return Term(
        name="transmit Gaussian efficiency",
        category="system",
        mean_db=float(loss),
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
    # geometry=None never flags (space budgets are always far field).
    assert term.assumptions.ok

    # Near-field flag: a hard-truncated beam (alpha=0.625) with a receiver inside
    # the Rayleigh range breaks the far-field eta. zR = pi*0.12^2/1.55e-6 ~ 29 km,
    # so a 1 km receiver flags; a 600 km one does not; and geometry=None does not.
    from types import SimpleNamespace
    near = tx_gaussian_efficiency_term(scn, SimpleNamespace(slant_range_m=1e3))
    assert not near.assumptions.ok and "Rayleigh range" in near.assumptions.violations[0]
    far = tx_gaussian_efficiency_term(scn, SimpleNamespace(slant_range_m=600e3))
    assert far.assumptions.ok
    # A widely open aperture (alpha above the threshold) does not flag, even near.
    open_scn = SpaceScenario(
        ground=Terminal(aperture_m=0.6, transmitter=Transmitter(waist_m=0.12)),
        space=Terminal(aperture_m=0.05), direction="uplink")   # alpha = 2.5
    assert tx_gaussian_efficiency_term(open_scn, SimpleNamespace(slant_range_m=1e3)).assumptions.ok

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
