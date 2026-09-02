'''
Several detectors behind one receive beamsplitter: one budget for each arm.

A receive telescope can feed MORE THAN ONE detector, for example a tracking
camera and a comms fibre behind a 10/90 beamsplitter. A Terminal holds ONE
detector, and about twenty detector dispatch sites in olb read that one field. So
this module does NOT make the Terminal hold a list. It makes ONE budget for each
arm:

    for each arm: copy the scenario with that arm's detector on the receive
                  terminal, call the matching budget function, then add the
                  fixed beamsplitter Term of that arm.

The two pieces come from olb.models.splitter: `resolve_fracs` gives the power
fraction of each arm, and `splitter_term` turns a fraction into the fixed dB
Term. An arm that takes all the power (fraction 1.0) gets NO splitter row.

WHY THIS IS EXACT. A beamsplitter multiplies the field of an arm by a constant.
So the arm keeps the SHAPE of the received field, and it loses only power. Every
coupling model in olb is power-normalised, so the coupling of an arm does not
change with the split ratio. The ratio therefore enters one time, as the fixed dB
Term. Source (the beamsplitter power-split relation): Saleh and Teich,
Fundamentals of Photonics, DOI 10.1002/0471213748.

AT FIDELITY 2 the arms share ONE Monte Carlo. Run
olb.models.waveoptics.run_fidelity2(scenario, geometry, detectors=[...]), which
gives a LIST of Fidelity2Bundle in the arm order, and pass that list as `wave`.
The clipped receive field is computed one time, and each arm is one more cheap
focal-plane calculation on that same array.

WHY THIS MODULE IS AT THE TOP LEVEL. The split is CROSS-CUTTING: it is not the
physics of one link. olb.links holds the per-link physics, so this module sits
ABOVE it in the dependency order. It may read olb.links, olb.models, and
olb.terminal; nothing there reads it back.
'''

from .models.splitter import arm_scenario, resolve_fracs, splitter_term


def _budget_function(scenario):
    '''
    Give the budget function of a scenario family and direction.

    A TerrestrialScenario always uses terrestrial_budget. A SpaceScenario uses
    the budget of its direction.

    Parameters:
        scenario : SpaceScenario or TerrestrialScenario

    Returns:
        callable
            A budget function f(scenario, geometry, **kwargs) -> Budget.
    '''
    # Import here: only the arm's own budget function is needed, so a lazy
    # import keeps the module light.
    if hasattr(scenario, "ground"):       # a SpaceScenario has ground/space
        from .links import (downlink_budget, retro_space_budget, uplink_budget)
        return {"uplink": uplink_budget,
                "downlink": downlink_budget,
                "retro": retro_space_budget}[scenario.direction]
    from .links import terrestrial_budget
    return terrestrial_budget


def multi_detector_budgets(scenario, geometry, detectors, *, wave=None,
                           **kwargs):
    '''
    Build one budget for each detector behind the receive beamsplitter.

    Parameters:
        scenario : SpaceScenario or TerrestrialScenario
            The link case. Its receive terminal detector is REPLACED by each
            arm's detector in turn. The input scenario does not change.
        geometry : CircularOrbit, TLEPass, or HorizontalPath
            The link geometry, passed to the budget function unchanged.
        detectors : sequence
            The detector of each arm (Aperture, SMF, MMF, or Camera). Each one
            carries an optional `frac`, the fraction of the received power in
            that arm. At most one may leave `frac` at None; that arm takes the
            remainder. See olb.models.splitter.resolve_fracs.
        wave : Fidelity2Bundle or list, optional
            The precomputed wave-optics record(s) for fidelity=2. Give a LIST
            (one bundle for each arm, in the `detectors` order) from
            olb.models.waveoptics.run_fidelity2(..., detectors=[...]), so every
            arm reads ONE shared Monte Carlo. A single bundle goes to every arm
            unchanged. None passes no `wave` to the budget function.
        **kwargs :
            Passed to the budget function unchanged (for example fidelity,
            turbulence, scintillation, tau_zenith).

    Returns:
        list of (detector, Budget)
            One pair for each arm, in the `detectors` order. Each Budget holds
            the arm's own Terms plus the fixed beamsplitter Term. An arm with a
            fraction of 1.0 gets no beamsplitter Term (it has no loss).

    Raises:
        ValueError
            If the fractions do not resolve (see resolve_fracs), or if a budget
            function refuses an arm. A per-arm error is NOT caught: a
            downlink_budget at fidelity 0 or 1 raises on a Camera, and that is
            the documented behaviour of that budget.

    LIMITS. The pattern is generic, so it works wherever a budget reads the
    receive terminal detector. It is TESTED for a terrestrial link and a
    downlink. An uplink receives on the SPACE terminal, so the arms are on the
    satellite. A RETRO link transmits and receives on the SAME ground terminal,
    so an arm copy changes the transmit terminal too (correct, but check that
    the retro budget reads what you expect).
    '''
    dets = list(detectors)
    fracs = resolve_fracs(dets)
    budget_of = _budget_function(scenario)
    is_list = isinstance(wave, (list, tuple))
    if is_list and len(wave) != len(dets):
        raise ValueError(
            f"multi_detector_budgets: `wave` holds {len(wave)} bundles, but "
            f"there are {len(dets)} detectors. Run "
            "olb.models.waveoptics.run_fidelity2(..., detectors=<the same "
            "list>) to get one bundle for each arm.")

    out = []
    for i, (det, frac) in enumerate(zip(dets, fracs)):
        kw = dict(kwargs)
        arm_wave = wave[i] if is_list else wave
        if arm_wave is not None:
            kw["wave"] = arm_wave
        budget = budget_of(arm_scenario(scenario, det), geometry, **kw)
        if frac < 1.0:
            budget.add(splitter_term(frac))
        out.append((det, budget))
    return out


if __name__ == '__main__':
    from .geometry import CircularOrbit, HorizontalPath
    from .scenario import (Channel, SpaceScenario, TerrestrialChannel,
                           TerrestrialScenario)
    from .terminal import Aperture, Camera, MMF, SMF, Terminal, Transmitter

    # --- a terrestrial link: a tracking camera plus a comms fibre ------------
    scn = TerrestrialScenario(
        near=Terminal(aperture_m=0.2, wavelength_m=1550e-9,
                      transmitter=Transmitter(waist_m=0.02, power_dbm=30.0)),
        far=Terminal(aperture_m=0.2, wavelength_m=1550e-9,
                     detector=Aperture()),
        channel=TerrestrialChannel(path_length_m=3e3, attenuation_db_per_km=0.5,
                                   cn2=3e-16))
    geom = HorizontalPath(3e3)

    cam = Camera(pixel_pitch_m=10e-6, n_pixels=128, focal_length_m=1.0,
                 frac=0.1)
    mmf = MMF(core_radius_m=25e-6, focal_length_m=0.05,
              numerical_aperture=0.2, sensitivity_dbm=-38.0)   # frac None: 0.9
    arms = multi_detector_budgets(scn, geom, [cam, mmf])
    assert len(arms) == 2
    (d0, b0), (d1, b1) = arms
    assert d0 is cam and d1 is mmf

    # Each arm carries ONE beamsplitter Term, and the loss is -10*log10(frac).
    def _splitter_db(budget):
        rows = [t for t in budget.terms if t.name == "beamsplitter"]
        assert len(rows) == 1, [t.name for t in budget.terms]
        return float(rows[0].mean_db)

    assert abs(_splitter_db(b0) - 10.0) < 1e-9, _splitter_db(b0)      # 10 percent
    assert abs(_splitter_db(b1) - 0.4576) < 1e-3, _splitter_db(b1)    # 90 percent

    # The camera arm is a power bucket at fidelity 0, so it gets the
    # scintillation Term. The MMF arm gets the light-bucket coupling Term.
    assert any("scintillation" in t.name for t in b0.terms), \
        [t.name for t in b0.terms]
    assert any(t.category == "coupling" for t in b1.terms), \
        [t.name for t in b1.terms]
    # The two arms share the deterministic backbone, term for term.
    for label in ("geometric", "atmospheric"):
        g0 = [t.mean_db for t in b0.terms if t.category == label]
        g1 = [t.mean_db for t in b1.terms if t.category == label]
        assert g0 == g1, (label, g0, g1)
    # The input scenario does not change.
    assert isinstance(scn.far.detector, Aperture)
    # The MMF arm reads its own sensitivity, so the margin is the arm's margin.
    assert b1.rx_sensitivity_dbm == -38.0 and b1.tx_power_dbm == 30.0

    # --- one arm alone takes all the power: no beamsplitter row -------------
    solo = multi_detector_budgets(scn, geom, [SMF()])
    assert not any(t.name == "beamsplitter" for t in solo[0][1].terms)

    # --- the fraction rule propagates ---------------------------------------
    try:
        multi_detector_budgets(scn, geom, [SMF(), MMF(core_radius_m=25e-6)])
        raise AssertionError("two frac=None must raise")
    except ValueError as exc:
        assert "frac=None" in str(exc), str(exc)

    # --- a downlink: the same pattern on a SpaceScenario --------------------
    down = SpaceScenario(
        ground=Terminal(aperture_m=0.4, wavelength_m=1550e-9,
                        detector=Aperture()),
        space=Terminal(aperture_m=0.1, wavelength_m=1550e-9,
                       transmitter=Transmitter(waist_m=0.04, power_dbm=30.0)),
        direction="downlink", channel=Channel())
    orbit = CircularOrbit(altitude_m=600e3, elevation_deg=[40.0])
    d_arms = multi_detector_budgets(down, orbit,
                                    [Aperture(frac=0.5), Aperture(frac=0.5)],
                                    fidelity=0)
    assert len(d_arms) == 2
    for _, b in d_arms:
        assert abs(_splitter_db(b) - 3.0103) < 1e-3, _splitter_db(b)
    # An arm error is NOT swallowed: downlink_budget at fidelity 0 knows no
    # Camera, so it raises. That is the documented behaviour of that budget.
    try:
        multi_detector_budgets(down, orbit, [cam, Aperture()], fidelity=0)
        raise AssertionError("a Camera on a fidelity-0 downlink must raise")
    except ValueError:
        pass

    # --- a mismatched `wave` list is refused ---------------------------------
    try:
        multi_detector_budgets(scn, geom, [cam, mmf], fidelity=2,
                               wave=[object()])
        raise AssertionError("a short wave list must raise")
    except ValueError as exc:
        assert "bundles" in str(exc), str(exc)

    print("terrestrial beamsplitter arms, 3 km:")
    for det, b in arms:
        name = type(det).__name__
        print(f"  {name:<8s} frac {(det.frac if det.frac else 0.9):4.2f}  "
              f"splitter {_splitter_db(b):6.3f} dB  "
              f"total {float(b.total_loss_db()):8.3f} dB  "
              f"({len(b.terms)} terms)")
    print("self-check passed")
