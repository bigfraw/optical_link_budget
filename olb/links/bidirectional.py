'''
Terrestrial bidirectional wrapper: one collimator drives both directions.

A monostatic terminal uses ONE collimator for transmit and receive. So ONE
fibre-plane defocus dz drives BOTH the transmit beam and the receive coupling.
Move the fibre (or the source) off the collimator focus by dz, and two things
happen together:

  - the transmit beam is no longer collimated. A source at dz behind the focus
    launches a diverging (or converging) beam, so the far-field divergence grows;
  - the receive detector no longer sits at the focus. The focused spot grows
    (see olb.models.coupling.terrestrial).

This module ties the two effects to one dz. `defocused_terminal` returns a NEW
Terminal with both the transmit divergence and the detector defocus set from dz.
`bidirectional_terrestrial` builds the forward (near->far) and reverse (far->near)
budgets with each terminal defocused by its own dz.

The transmit divergence follows the geometric (thin-lens) launch. A source at dz
behind the collimator focus gives a phase-front radius of curvature at the
aperture of R = f^2/dz. The far-field HALF-angle divergence is then

    theta(dz) = sqrt( theta_diff^2 + (W0*|dz|/f^2)^2 ),   theta_diff = lambda/(pi*W0),

with W0 the transmit waist at the aperture and f the collimator focal length.
theta_diff is the diffraction limit; the second term is the geometric spread of a
displaced source. At dz=0 the beam is collimated. Source: Andrews and Phillips,
2nd ed. (2005), DOI 10.1117/3.626196, Ch. 4 (Gaussian beam radius of curvature).

TWO LIMITS of this fidelity-0 wrapper:

1. Only the DIVERGING side is modelled. theta(|dz|) depends on |dz| only, and a
   Transmitter has no way to say "converging". So dz > 0 (the fibre BEYOND the
   collimator focus, which launches a converging beam) falls OUTSIDE the model:
   the wrapper gives it the divergence of the mirror-image diverging launch. Use
   dz < 0, or a fidelity-2 field model, for a converging launch.
2. ONE dz drives BOTH directions. So a deliberately diverged monostatic terminal
   also moves its own receive fibre off the focal plane. The received beam is
   itself a diverging Gaussian, so its true focus is already at +dz_curv beyond
   the focal plane (see olb.models.coupling.terrestrial). The terminal therefore
   pays |dz| + dz_curv of receive defocus, and the coupling Terms now CHARGE
   that: the curvature defocus is always charged, at dz_eff = dz - dz_curv.
'''

from collections import namedtuple
from dataclasses import replace

import numpy as np

from ..terminal import SMF, MMF
from ..scenario import TerrestrialScenario
from ..units import w0_to_div
from ..models.coupling.terrestrial import _smf_optics, _mmf_focal_length
from .terrestrial import terrestrial_budget

# The forward (near->far) and reverse (far->near) budgets of one bidirectional
# link. A light wrapper, no heavy class.
BidirectionalBudget = namedtuple("BidirectionalBudget", ["forward", "reverse"])


def _resolve_focal_length(terminal, focal_length_m):
    '''
    Return the collimator focal length of a terminal, or None if it cannot be
    resolved.

    An explicit focal_length_m always wins. Otherwise the model reads the detector
    optics: an SMF uses _smf_optics, an MMF uses _mmf_focal_length (each honours
    optimal_focus). An Aperture or no detector gives None unless focal_length_m is
    set.
    '''
    if focal_length_m is not None:
        return focal_length_m
    detector = terminal.detector
    if isinstance(detector, SMF):
        f, _ = _smf_optics(detector, terminal.aperture_m, terminal.wavelength_m)
        return f
    if isinstance(detector, MMF):
        return _mmf_focal_length(detector, terminal.aperture_m, terminal.wavelength_m)
    return None


def defocused_terminal(terminal, dz_m, *, focal_length_m=None):
    '''
    Return a NEW Terminal with the fibre-plane defocus dz_m applied.

    A monostatic terminal has one collimator, so one defocus dz_m drives both the
    transmit divergence and the receive coupling. The input Terminal is NOT
    mutated (the function uses dataclasses.replace).

    The detector defocus. When the terminal has an SMF or MMF detector, the new
    detector carries defocus_m = dz_m. An Aperture or no detector is left
    unchanged.

    The transmit divergence. When the terminal has a Transmitter, the new
    transmitter carries divergence_rad from dz_m:

        theta(dz) = sqrt( theta_diff^2 + (W0*|dz|/f^2)^2 ),  theta_diff = lambda/(pi*W0)

    with W0 = transmitter.waist_m and f the collimator focal length (see
    _resolve_focal_length). At dz=0 the beam is collimated, so divergence_rad is
    None. If the computed theta is at or below the diffraction limit, it stays
    None. Source: Andrews and Phillips, 2nd ed. (2005), DOI 10.1117/3.626196,
    Ch. 4.

    LIMITS (see the module docstring). theta(dz) reads |dz| only, so a POSITIVE
    dz_m (the fibre beyond the collimator focus, that is a CONVERGING launch) is
    outside this fidelity-0 divergence model: a Transmitter cannot hold a
    converging beam, so the wrapper gives it the divergence of the mirror-image
    diverging launch. And the one dz_m drives the RECEIVE side too: the received
    beam is a diverging Gaussian whose true focus is already dz_curv beyond the
    focal plane, so a deliberately diverged monostatic terminal pays
    |dz| + dz_curv of receive defocus. The coupling Terms now charge this (see
    olb.models.coupling.terrestrial).

    Parameters:
        terminal : Terminal
            The monostatic terminal to defocus.
        dz_m : float
            Fibre-plane defocus [m]. 0.0 leaves the terminal collimated and at
            focus.
        focal_length_m : float, optional
            The collimator focal length [m]. None reads it from the detector
            optics.

    Returns:
        Terminal
            A new Terminal with the defocus applied.

    Raises:
        ValueError
            If the terminal has a Transmitter but the focal length cannot be
            resolved (no explicit focal_length_m and no fibre-coupling optics to
            read it from). The defocus-to-divergence map needs the focal length.
    '''
    detector = terminal.detector
    if isinstance(detector, (SMF, MMF)):
        new_detector = replace(detector, defocus_m=dz_m)
    else:
        new_detector = detector

    new_transmitter = terminal.transmitter
    if terminal.transmitter is not None:
        f = _resolve_focal_length(terminal, focal_length_m)
        if f is None:
            raise ValueError(
                "defocused_terminal needs the collimator focal length to map the "
                "defocus to a transmit divergence. Pass focal_length_m, or give the "
                "terminal an SMF/MMF detector with the coupling optics set "
                "(focal_length_m or optimal_focus)."
            )
        W0 = terminal.transmitter.waist_m
        wavelength = terminal.wavelength_m
        theta_diff = w0_to_div(W0, wavelength)
        if dz_m == 0.0:
            divergence = None
        else:
            theta = float(np.sqrt(theta_diff ** 2 + (W0 * abs(dz_m) / f ** 2) ** 2))
            # Keep None at or below the diffraction limit (a numerically tiny dz).
            divergence = None if theta <= theta_diff * (1 + 1e-12) else theta
        new_transmitter = replace(terminal.transmitter, divergence_rad=divergence)

    return replace(terminal, transmitter=new_transmitter, detector=new_detector)


def bidirectional_terrestrial(near, far, channel, geometry, *,
                              near_defocus_m=0.0, far_defocus_m=0.0,
                              **budget_kwargs):
    '''
    Build the forward and reverse terrestrial budgets of a bidirectional link.

    Each terminal is a monostatic terminal, so its own defocus drives both its
    transmit divergence and its receive coupling. The forward budget sends
    near->far (tx = near, rx = far); the reverse budget sends far->near (tx = far,
    rx = near). Both directions share the one TerrestrialChannel.

    Parameters:
        near, far : Terminal
            The two ends of the path. Each is monostatic (a Transmitter and a
            detector) for a real bidirectional link.
        channel : TerrestrialChannel
            The horizontal propagation channel.
        geometry : HorizontalPath
            The horizontal path geometry.
        near_defocus_m, far_defocus_m : float
            The fibre-plane defocus [m] of the near and far terminals.
        **budget_kwargs :
            Passed straight to terrestrial_budget (for example fidelity,
            turbulence, scintillation).

    Returns:
        BidirectionalBudget
            A namedtuple (forward, reverse) of the two Budgets.
    '''
    near_cfg = defocused_terminal(near, near_defocus_m)
    far_cfg = defocused_terminal(far, far_defocus_m)
    forward = terrestrial_budget(
        TerrestrialScenario(near=near_cfg, far=far_cfg, channel=channel),
        geometry, **budget_kwargs)
    reverse = terrestrial_budget(
        TerrestrialScenario(near=far_cfg, far=near_cfg, channel=channel),
        geometry, **budget_kwargs)
    return BidirectionalBudget(forward=forward, reverse=reverse)


if __name__ == '__main__':
    from ..terminal import Terminal, Transmitter, Aperture
    from ..scenario import TerrestrialChannel
    from ..geometry import HorizontalPath

    lam = 1550e-9
    chan = TerrestrialChannel(path_length_m=5e3, attenuation_db_per_km=0.5, cn2=1e-15)
    geom = HorizontalPath(5e3)

    def _mono(jitter=5e-6):
        '''A monostatic terminal: a launch beam and an MMF light bucket.'''
        return Terminal(aperture_m=0.2, wavelength_m=lam, pointing_jitter_rad=jitter,
                        transmitter=Transmitter(waist_m=0.05, power_dbm=30),
                        detector=MMF(core_radius_m=25e-6, optimal_focus=True,
                                     sensitivity_dbm=-38))

    near, far = _mono(), _mono()
    theta_diff = w0_to_div(0.05, lam)

    # --- defocused_terminal ---------------------------------------------------
    # dz=0: the transmit beam is collimated and the detector is at focus.
    t0 = defocused_terminal(near, 0.0)
    assert t0.transmitter.divergence_rad is None      # collimated
    assert t0.detector.defocus_m == 0.0               # at focus
    assert t0 is not near and t0.transmitter is not near.transmitter   # no mutation
    assert near.transmitter.divergence_rad is None    # input unchanged

    # A nonzero dz: the divergence exceeds the diffraction limit and the detector
    # defocus is set.
    t1 = defocused_terminal(near, 3e-3)
    assert t1.transmitter.divergence_rad is not None
    assert t1.transmitter.divergence_rad > theta_diff, (
        t1.transmitter.divergence_rad, theta_diff)
    assert t1.detector.defocus_m == 3e-3
    # A larger dz gives a larger divergence (a stronger geometric spread).
    t2 = defocused_terminal(near, 6e-3)
    assert t2.transmitter.divergence_rad > t1.transmitter.divergence_rad

    # An explicit focal length overrides the detector optics.
    t_expl = defocused_terminal(near, 3e-3, focal_length_m=0.1)
    W0 = 0.05
    theta_expl = np.sqrt(theta_diff ** 2 + (W0 * 3e-3 / 0.1 ** 2) ** 2)
    assert np.isclose(t_expl.transmitter.divergence_rad, theta_expl)

    # A transmitter with no resolvable focal length raises.
    bad = Terminal(aperture_m=0.2, wavelength_m=lam,
                   transmitter=Transmitter(waist_m=0.05), detector=Aperture())
    try:
        defocused_terminal(bad, 3e-3)
        raise AssertionError("a transmitter with no focal length must raise")
    except ValueError:
        pass
    # The same terminal with an explicit focal length works.
    ok = defocused_terminal(bad, 3e-3, focal_length_m=0.1)
    assert ok.transmitter.divergence_rad > theta_diff
    # A receive-only terminal (no transmitter) needs no focal length.
    rx_only = Terminal(aperture_m=0.2, wavelength_m=lam,
                       detector=MMF(core_radius_m=25e-6, optimal_focus=True))
    assert defocused_terminal(rx_only, 3e-3).detector.defocus_m == 3e-3

    # --- bidirectional_terrestrial --------------------------------------------
    b0 = bidirectional_terrestrial(near, far, chan, geom)
    assert b0.forward.scenario.tx_terminal.transmitter.divergence_rad is None
    assert b0.reverse.scenario.tx_terminal.transmitter.divergence_rad is None
    assert np.isfinite(b0.forward.total_loss_db())
    assert np.isfinite(b0.reverse.total_loss_db())

    b1 = bidirectional_terrestrial(near, far, chan, geom,
                                   near_defocus_m=3e-3, far_defocus_m=3e-3)
    assert np.isfinite(b1.forward.total_loss_db())
    assert np.isfinite(b1.reverse.total_loss_db())
    # The forward tx (near) is now diverging, and the forward rx (far) is defocused.
    assert b1.forward.scenario.tx_terminal.transmitter.divergence_rad > theta_diff
    assert b1.forward.scenario.rx_terminal.detector.defocus_m == 3e-3

    # A larger dz widens the received beam (a bigger geometric spreading loss) AND
    # changes the coupling loss.
    def _geo(budget):
        return next(t for t in budget.terms if t.name == "geometric spreading").mean_db

    def _cpl(budget):
        return next(t for t in budget.terms if t.category == "coupling").mean_db

    geo0, geo1 = _geo(b0.forward), _geo(b1.forward)
    assert geo1 > geo0, (geo1, geo0)                  # the beam spreads more
    cpl0, cpl1 = _cpl(b0.forward), _cpl(b1.forward)
    assert cpl1 != cpl0, (cpl1, cpl0)                 # the coupling changes

    # A still larger dz spreads the beam even more.
    b2 = bidirectional_terrestrial(near, far, chan, geom,
                                   near_defocus_m=6e-3, far_defocus_m=6e-3)
    assert _geo(b2.forward) > geo1

    print(f"diffraction limit: {theta_diff * 1e6:.2f} urad")
    print(f"divergence at dz=3 mm: {t1.transmitter.divergence_rad * 1e6:.2f} urad")
    print(f"forward geometric loss: dz=0 {geo0:.2f} dB  dz=3mm {geo1:.2f} dB  "
          f"dz=6mm {_geo(b2.forward):.2f} dB")
    print(f"forward coupling loss:  dz=0 {cpl0:.2f} dB  dz=3mm {cpl1:.2f} dB")
    print(f"forward total: dz=0 {b0.forward.total_loss_db():.2f} dB  "
          f"dz=3mm {b1.forward.total_loss_db():.2f} dB")
    print("self-check passed")
