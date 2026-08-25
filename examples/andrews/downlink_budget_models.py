'''
The capstone: one downlink budget, three scintillation models, two elevations.

Andrews and Phillips, "Laser Beam Propagation through Random Media", 2nd ed.,
SPIE Press (2005), DOI: 10.1117/3.626196:
    Ch. 12, Eq. (38), printed p. 495   the slant plane-wave Rytov variance
    Ch. 9, Eqs. (41), (46), printed pp. 335, 336   the two log variances
    Ch. 9, Eq. (138), printed p. 370   gamma-gamma alpha and beta
    Ch. 9, Eq. (139), printed p. 371   the gamma-gamma scintillation index
    Ch. 11, Sec. 11.3, printed p. 451  why the lognormal tail is too thin
The dB faces come from olb/models/fade.py, the one adapter that turns an
irradiance model into the three faces of a Term.

THE SELECTOR. downlink_scintillation_term(..., model="auto") reads the point
index sigma_I^2. Below the house limit 0.25 it returns the LOGNORMAL Term. At or
above it, the GAMMA-GAMMA Term. The house limit is four times stricter than the
book limit sigma_R^2 = 1, and the gamma-gamma chain is valid at every strength,
so the early switch costs no validity.

THE APERTURE-AVERAGING CAVEAT. The gamma-gamma Term models a POINT receiver. The
book gives NO aperture-averaged downlink index in the moderate-to-strong regime:
Ch. 12, Eq. (39), printed p. 496, is a weak form and Ch. 12, Eq. (40), printed
p. 497, is a point form, and the book prints no product of the two. So the
gamma-gamma fade below is DEEPER than a real 70 cm aperture would see. That is
the safe direction, and the Term flags it through its Assumptions record. Read
the "aperture A" column: 1.000 means no averaging was applied.

HOW THE BUDGET IS BUILT. downlink_budget takes no `model` keyword today; it
always asks for the lognormal Term. So this script calls it with
scintillation=False and then adds the chosen Term. Everything else in the budget
(geometric, extinction, pointing) is the package default.

Run from the repo root:
    python -m examples.andrews.downlink_budget_models
'''

import warnings

from olb import (Channel, CircularOrbit, Site, SpaceScenario, Terminal,
                 Transmitter)
from olb.links.downlink import downlink_budget, downlink_scintillation_term

# --- configuration ----------------------------------------------------------
WAVELENGTH_M = 1550e-9
ALTITUDE_M = 600e3
GROUND_APERTURE_M = 0.7
AVAILABILITY = 0.99
MODELS = ("lognormal", "gamma_gamma", "auto")
ELEVATIONS_DEG = (60.0, 15.0)       # one weak case, one strong case


def scenario():
    '''Build the downlink case: the satellite transmits, the ground receives.'''
    ground = Terminal(aperture_m=GROUND_APERTURE_M, wavelength_m=WAVELENGTH_M)
    space = Terminal(aperture_m=0.05, wavelength_m=WAVELENGTH_M,
                     transmitter=Transmitter(waist_m=0.035, power_dbm=30.0))
    return SpaceScenario(ground=ground, space=space, direction="downlink",
                         channel=Channel(site=Site(cn2_ground=1.7e-14),
                                         altitude_m=ALTITUDE_M))


def build(scn, geom, model):
    '''Return the full downlink budget with ONE chosen scintillation Term.'''
    budget = downlink_budget(scn, geom, scintillation=False)
    term = downlink_scintillation_term(scn, geom, model=model,
                                       aperture_average=True)
    return budget.add(term), term


def print_elevation(scn, elevation_deg):
    '''Print the three models at one elevation, plus the Term itemisation.'''
    geom = CircularOrbit(ALTITUDE_M, elevation_deg=elevation_deg)
    print(f"elevation {elevation_deg:.0f} deg, slant range "
          f"{geom.slant_range_m/1e3:.0f} km, rx D={GROUND_APERTURE_M*100:.0f} cm")
    print(f"  {'model':>12} {'selected':>12} | {'point s_I^2':>12} "
          f"{'aperture A':>11} {'used s^2':>10} | {'scint mean':>11} "
          f"{'total mean':>11} {'total 99%':>10}")
    print("  " + "-" * 98)
    metas = {}
    for model in MODELS:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            budget, term = build(scn, geom, model)
            fade = float(budget.fade_margin_db(AVAILABILITY))
        meta = metas[model] = term.meta
        print(f"  {model:>12} {meta['model']:>12} | "
              f"{float(meta['sigma2_I']):>12.5f} "
              f"{float(meta['aperture_averaging_factor']):>11.4f} "
              f"{float(meta['sigma2_P']):>10.5f} | "
              f"{float(term.mean_db):>11.3f} "
              f"{float(budget.total_loss_db()):>11.3f} {fade:>10.3f}")
    print(f"  house limit {meta['weak_fluctuation_limit']}. The lognormal "
          f"point index IS the Rytov variance of\n  Ch. 12, Eq. (38). The "
          f"gamma-gamma point index is Ch. 9, Eq. (139); it comes from the "
          f"SAME\n  Rytov variance sigma_R^2 = "
          f"{float(metas['gamma_gamma']['sigma2_R']):.5f}.")
    ok = term.assumptions.ok if term.assumptions is not None else True
    print(f"  the 'auto' Term assumptions hold: {ok}")
    if not ok:
        for line in term.assumptions.violations:
            print(f"    flag: {line}")
    print()


def print_itemised(scn):
    '''Print the itemised budget once, so a reader sees every Term.'''
    geom = CircularOrbit(ALTITUDE_M, elevation_deg=ELEVATIONS_DEG[-1])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        budget, _ = build(scn, geom, "auto")
    print(f"itemised budget at {ELEVATIONS_DEG[-1]:.0f} deg, model='auto'")
    print(budget.to_frame().to_string(index=False))
    print()


if __name__ == '__main__':
    scn = scenario()
    for elevation_deg in ELEVATIONS_DEG:
        print_elevation(scn, elevation_deg)
    print_itemised(scn)
    print("Reading of the tables. At the WEAK elevation the selector keeps the "
          "lognormal\nTerm and the aperture-averaging factor is real, so the "
          "fade is shallow. At the\nSTRONG elevation the selector switches to "
          "gamma-gamma, the averaging factor\ndrops to 1.000 (POINT receiver), "
          "and the 99 % fade grows for two reasons at\nonce: the heavier tail "
          "AND the lost aperture averaging. Do not read that fade\nas the fade "
          "of a 70 cm telescope. See the module docstring.")
