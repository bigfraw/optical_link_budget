'''
The receive beamsplitter: the power split into several detector arms.

A receive telescope can feed MORE THAN ONE detector. A beamsplitter divides the
collected power between the arms, for example 10 percent to a tracking camera and
90 percent to the comms fibre. This module holds the two small pieces that the
split needs:

  - `resolve_fracs` turns the `frac` field of a set of detectors into one
    fraction for each arm. See the rule below.
  - `splitter_term` turns one fraction into a deterministic budget Term.

THE PHYSICS IS ONE LINE. A beamsplitter multiplies the field of each arm by a
constant. So the arm keeps the SHAPE of the received field, and it loses only
power. Every coupling efficiency in olb is power-normalised (a ratio of the arm
power to the arm input power), so the efficiency of an arm does NOT change with
the split ratio. The split ratio therefore enters the budget one time, as this
fixed dB Term. Source (the beamsplitter power-split relation): Saleh and Teich,
Fundamentals of Photonics, DOI 10.1002/0471213748.

THE FRACTION RULE. At most ONE detector may leave `frac` at None. That detector
takes the remainder, 1 - sum(the others). One detector alone with frac=None takes
1.0. When every fraction is given, the sum must not be more than 1.0. A sum below
1.0 is allowed: the missing part is the excess loss of the splitter, and it does
not go to any arm.

THIS IS NOT A TURBULENCE TERM. It is a fixed optical loss, so it declares no beam
type, no turbulence regime, and no spectrum. Budget.check() applies its untraced
provenance guard to turbulence and coupling Terms only, so the "system" category
of this Term needs no traced provenance.
'''

from dataclasses import replace

import numpy as np

from ..results import Term  # run with `python -m olb.models.splitter`
from ..assumptions import Assumptions, BEAM_NA, REGIME_NA, SPECTRUM_NA

# A fraction sum may exceed 1.0 by this much and still count as 1.0. Floating
# point makes 0.1 + 0.9 differ from 1.0 in the last bits. This is a numerical
# tolerance, not a physical limit.
_SUM_TOL = 1e-9


def resolve_fracs(detectors):
    '''
    Give the beamsplitter fraction of each detector.

    At most one detector may have frac=None. That detector takes the remainder,
    1 - sum(the others). A single detector with frac=None takes 1.0.

    Parameters:
        detectors : sequence
            The detector objects of the receive arms (Aperture, SMF, MMF, or
            Camera). Each one carries an optional `frac`.

    Returns:
        list of float
            One fraction for each detector, in the input order.

    Raises:
        ValueError
            If the sequence is empty, if two or more detectors have frac=None,
            if a given fraction is outside (0, 1], if the given fractions add up
            to more than 1.0, or if the remainder is not more than 0.
    '''
    dets = list(detectors)
    if not dets:
        raise ValueError("resolve_fracs: give at least one detector.")

    fracs = [getattr(d, "frac", None) for d in dets]
    missing = [i for i, f in enumerate(fracs) if f is None]
    if len(missing) > 1:
        raise ValueError(
            f"resolve_fracs: {len(missing)} detectors (indices {missing}) have "
            "frac=None. At most ONE arm can take the remainder. Set frac on all "
            "but one detector.")

    given = [float(f) for f in fracs if f is not None]
    for f in given:
        if not (0.0 < f <= 1.0):
            raise ValueError(
                f"resolve_fracs: a beamsplitter fraction must be in (0, 1], got "
                f"{f!r}.")
    total = sum(given)
    if total > 1.0 + _SUM_TOL:
        raise ValueError(
            f"resolve_fracs: the given fractions add up to {total:g}, which is "
            "more than 1.0. A beamsplitter cannot give an arm more power than it "
            "receives.")

    if not missing:
        return [float(f) for f in fracs]

    remainder = 1.0 - total
    if remainder <= 0.0:
        raise ValueError(
            f"resolve_fracs: the given fractions add up to {total:g}, so the "
            "remainder arm gets no power. Lower the other fractions, or set a "
            "fraction on every arm.")
    out = [float(f) if f is not None else remainder for f in fracs]
    return out


def splitter_term(frac, *, name=None, note=None):
    '''
    The beamsplitter power split of one arm, as a deterministic Term.

    The loss is the plain dB value of the power fraction:

        loss_db = -10*log10(frac)

    This is the definition of a dB power ratio, not a model. The split relation
    itself is the beamsplitter of Saleh and Teich, Fundamentals of Photonics,
    DOI 10.1002/0471213748.

    Parameters:
        frac : float
            The fraction of the received power in this arm, in (0, 1].
        name : str, optional
            Override the default Term name.
        note : str, optional
            Override the default Term note.

    Returns:
        Term
            A deterministic Term (loss positive dB). frac=1.0 gives 0.0 dB.

    Raises:
        ValueError
            If frac is not in (0, 1].
    '''
    f = float(frac)
    if not (0.0 < f <= 1.0):
        raise ValueError(
            f"splitter_term: the beamsplitter fraction must be in (0, 1], got "
            f"{frac!r}.")
    loss_db = -10.0 * np.log10(f)
    if name is None:
        name = "beamsplitter"
    if note is None:
        note = f"receive beamsplitter, {100.0 * f:.4g} percent of the power"
    return Term(
        name=name,
        category="system",
        mean_db=loss_db,
        note=note,
        meta={"frac": f},
        assumptions=Assumptions(
            beam_type=BEAM_NA,
            turbulence_regime=REGIME_NA,
            spectrum=SPECTRUM_NA,
            validity="A fixed power split. The beamsplitter multiplies the field "
                     "of the arm by a constant, so the arm keeps the field SHAPE "
                     "and it loses only power. Every coupling efficiency is "
                     "power-normalised, so the split ratio does not change it. "
                     "The excess loss of the real splitter is the part that the "
                     "fractions leave out of 1.0.",
        ),
    )


def arm_scenario(scenario, detector):
    '''
    Copy a scenario with ONE detector on the receive terminal.

    A Terminal holds one detector, and the ~20 detector dispatch sites in olb
    read that one field. So an arm of a beamsplitter is a COPY of the scenario
    whose receive terminal holds that arm's detector. Every model then runs
    unchanged.

    The function replaces the correct terminal field for each scenario family: a
    SpaceScenario resolves the receive role from its direction (downlink and
    retro receive on `ground`, an uplink on `space`), and a TerrestrialScenario
    resolves it the same way (forward receives on `far`, reverse on `near`). See
    olb.scenario.

    NOTE, a retro link. A retro SpaceScenario transmits and receives on the SAME
    ground terminal, so this copy changes the transmit terminal too. That is
    correct: it is one physical terminal.

    Parameters:
        scenario : SpaceScenario or TerrestrialScenario
            The link case.
        detector : Aperture, SMF, MMF, or Camera
            The detector of this arm.

    Returns:
        SpaceScenario or TerrestrialScenario
            A copy. The input scenario does not change.
    '''
    arm_rx = replace(scenario.rx_terminal, detector=detector)
    if hasattr(scenario, "ground"):        # a SpaceScenario has ground/space
        role = "ground" if scenario.direction in ("downlink", "retro") else "space"
    else:                                  # a TerrestrialScenario has near/far
        role = "near" if scenario.direction == "reverse" else "far"
    return replace(scenario, **{role: arm_rx})


if __name__ == '__main__':
    from ..terminal import Aperture, Camera, MMF, SMF

    # --- the remainder rule --------------------------------------------------
    cam = Camera(pixel_pitch_m=10e-6, n_pixels=128, frac=0.1)
    mmf = MMF(core_radius_m=25e-6, focal_length_m=0.05)      # frac is None
    assert resolve_fracs([cam, mmf]) == [0.1, 0.9]
    # The order does not matter.
    assert resolve_fracs([mmf, cam]) == [0.9, 0.1]

    # --- a single detector with frac=None takes everything -------------------
    assert resolve_fracs([SMF()]) == [1.0]
    assert resolve_fracs([Aperture(frac=0.5)]) == [0.5]      # 0.5 excess loss

    # --- every fraction given: a sum below 1.0 is splitter excess loss --------
    three = [Aperture(frac=0.25), Aperture(frac=0.25), Aperture(frac=0.4)]
    assert resolve_fracs(three) == [0.25, 0.25, 0.4]

    # --- the failures --------------------------------------------------------
    try:
        resolve_fracs([SMF(), MMF(core_radius_m=25e-6)])     # two None
        raise AssertionError("two frac=None must raise")
    except ValueError as exc:
        assert "at most ONE" in str(exc) or "At most ONE" in str(exc), str(exc)
    try:
        resolve_fracs([Aperture(frac=0.7), Aperture(frac=0.6)])
        raise AssertionError("a sum above 1.0 must raise")
    except ValueError as exc:
        assert "more than 1.0" in str(exc), str(exc)
    try:
        resolve_fracs([Aperture(frac=0.6), Aperture(frac=0.4), SMF()])
        raise AssertionError("a zero remainder must raise")
    except ValueError as exc:
        assert "no power" in str(exc), str(exc)
    try:
        resolve_fracs([])
        raise AssertionError("an empty sequence must raise")
    except ValueError:
        pass
    try:
        resolve_fracs([Aperture(frac=0.0)])
        raise AssertionError("frac=0 must raise")
    except ValueError:
        pass

    # --- the Term ------------------------------------------------------------
    t_full = splitter_term(1.0)
    assert abs(t_full.mean_db) < 1e-12 and t_full.category == "system"
    assert t_full.sampler is None and t_full.quantile is None   # deterministic
    t_half = splitter_term(0.5)
    assert abs(t_half.mean_db - 3.0103) < 1e-3, t_half.mean_db
    t10 = splitter_term(0.1)
    assert abs(t10.mean_db - 10.0) < 1e-9, t10.mean_db
    t90 = splitter_term(0.9)
    assert abs(t90.mean_db - 0.4576) < 1e-3, t90.mean_db
    assert t90.meta["frac"] == 0.9
    # A deterministic Term's quantile is its mean, so it locks no fidelity.
    assert t90.quantile_db(0.99) == t90.mean_db and not t90.mean_only
    try:
        splitter_term(0.0)
        raise AssertionError("frac=0 must raise")
    except ValueError:
        pass

    # --- the arm scenario: the RIGHT terminal field changes ------------------
    from ..scenario import (Channel, SpaceScenario, TerrestrialChannel,
                            TerrestrialScenario)
    from ..terminal import Terminal, Transmitter

    near = Terminal(aperture_m=0.2, transmitter=Transmitter(waist_m=0.02))
    far = Terminal(aperture_m=0.2, detector=SMF())
    terr = TerrestrialScenario(near=near, far=far,
                               channel=TerrestrialChannel(path_length_m=3e3))
    arm = arm_scenario(terr, cam)
    assert arm.rx_terminal.detector is cam and arm.far.detector is cam
    assert arm.near is near                       # the transmit end is untouched
    assert terr.far.detector is not cam           # the input does not change
    rev = TerrestrialScenario(near=near, far=far, direction="reverse")
    assert arm_scenario(rev, cam).near.detector is cam
    assert arm_scenario(rev, cam).far is far

    gnd = Terminal(aperture_m=0.4, transmitter=Transmitter(waist_m=0.06))
    sat = Terminal(aperture_m=0.3, detector=Aperture())
    down = SpaceScenario(ground=gnd, space=sat, direction="downlink",
                         channel=Channel())
    assert arm_scenario(down, cam).ground.detector is cam
    assert arm_scenario(down, cam).space is sat
    up = SpaceScenario(ground=gnd, space=sat, direction="uplink")
    assert arm_scenario(up, cam).space.detector is cam
    assert arm_scenario(up, cam).ground is gnd

    print("beamsplitter fractions and Terms:")
    for f in resolve_fracs([cam, mmf]):
        print(f"  frac {f:5.2f}  ->  {float(splitter_term(f).mean_db):7.4f} dB")
    print("self-check passed")
