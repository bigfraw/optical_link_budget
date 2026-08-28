'''
Wire the fidelity-2 wave optics into the three budgets (the whole-path way).

Fidelity is a WHOLE-PATH choice (see the README fidelity ladder). At fidelity=2
the ENTIRE path is a field simulation, and the budget shows TWO wave-optics Terms:

  - a DETERMINISTIC vacuum-optics Term: the full no-turbulence loss from launch to
    detector (launch truncation + geometric spread + aperture capture + vacuum
    fibre coupling);
  - a STOCHASTIC turbulence Term: the fade.

Together they REPLACE the analytic geometric, launch-truncation, scintillation,
and coupling Terms. Only the analytic extinction (molecular absorption, never in
the field sim) and pointing (mechanical jitter) Terms stay.

This example touches three links, one selector each:
  - TERRESTRIAL SMF: terrestrial_budget(fidelity=2, wave=...). The default
    (fidelity=0) SMF Term is mean-only, so the default budget REFUSES a fade
    margin. Fidelity 2 gives one, and computes the FULL loss directly.
  - UPLINK (uncorrected): uplink_budget(fidelity=2, wave=...). The vacuum-optics
    Term is the geometric loss over 600 km; the turbulence Term is the reciprocity
    penalty.
  - DOWNLINK aperture: downlink_budget(fidelity=2, wave=...).

THE ONE RULE. The budget NEVER runs the split-step layer itself. A caller runs
BOTH propagations ONE time with olb.models.coupling.run_fidelity2 (the turbulent
Monte Carlo and the vacuum field solve), then gives the bundle to the budget.

Every run here uses the RAPID preset and a small trial count, so it finishes in a
couple of minutes; raise N_TRIALS and the preset for a real number. The tail
adequacy of a quantile is a property of the sampler
(olb.results.EmpiricalSampler): term.sampler.undersampled(p) answers it, and the
quantile warns when a run is too short.

Run from the repo root:
    python -m examples.waveoptics.budget_wiring
'''

import time
import warnings

import numpy as np

from olb import SMF, Terminal, Transmitter
from olb.geometry import CircularOrbit, HorizontalPath
from olb.models.coupling import run_fidelity2
from olb.scenario import (Channel, SpaceScenario, TerrestrialChannel,
                          TerrestrialScenario)
from olb.links.terrestrial import terrestrial_budget
from olb.links.uplink import uplink_budget
from olb.links.downlink import downlink_budget
from olb.turbulence.profiles import default_cn2_profile
from olb.waveoptics import Threader

WAVELENGTH_M = 1550e-9
PRESET = "rapid"       # a demonstration: RAPID keeps the example short
N_TRIALS = 200
SEED = 20260828
AVAILABILITY = 0.9     # the fade availability that every case reports
THREADER = Threader()  # the trials run across threads (the FFT releases the GIL)


def _line():
    print("-" * 70)


def show_fidelity2(budget):
    '''Print the two wave-optics Terms and the analytic backbone of an F2 budget.'''
    vac = next(t for t in budget.terms if t.meta.get("model") == "waveoptics-vacuum")
    turb = next(t for t in budget.terms if t.meta.get("model") == "waveoptics")
    fade = turb.quantile_db(AVAILABILITY)
    print(f"    vacuum optics (deterministic)   {vac.mean_db:8.3f} dB")
    print(f"    turbulence (wave optics) mean   {turb.mean_db:8.3f} dB")
    print(f"    turbulence {int(AVAILABILITY * 100)}% fade          {fade:8.3f} dB")
    print(f"    budget total (+ extinction, pointing)  "
          f"{budget.total_loss_db():8.3f} dB")


def terrestrial_case():
    '''Terrestrial SMF: fidelity 2 unlocks the fade margin and gives the full loss.'''
    print("=" * 70)
    print("TERRESTRIAL SMF")
    _line()
    scenario = TerrestrialScenario(
        near=Terminal(aperture_m=0.2, wavelength_m=WAVELENGTH_M,
                      transmitter=Transmitter(waist_m=0.02, power_dbm=30)),
        far=Terminal(aperture_m=0.2, wavelength_m=WAVELENGTH_M,
                     detector=SMF(sensitivity_dbm=-40)),
        channel=TerrestrialChannel(path_length_m=3e3, attenuation_db_per_km=0.5,
                                   cn2=5e-15))
    geometry = HorizontalPath(3e3)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        default = terrestrial_budget(scenario, geometry)          # fidelity 0
    dcoup = next(t for t in default.terms if t.category == "coupling")
    print(f"  fidelity 0 (analytic, default): coupling {dcoup.mean_db:.3f} dB, "
          f"mean_only={dcoup.mean_only}, provides_fade={default.provides_fade}")

    t0 = time.time()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        bundle = run_fidelity2(scenario, geometry, n_trials=N_TRIALS, seed=SEED,
                               preset=PRESET, threader=THREADER)
        f2 = terrestrial_budget(scenario, geometry, fidelity=2, wave=bundle)
    print(f"  fidelity 2 (wave optics, {N_TRIALS} snapshots, {PRESET}): "
          f"provides_fade {default.provides_fade} -> {f2.provides_fade}")
    show_fidelity2(f2)
    print(f"  ({time.time() - t0:.1f} s)")
    print("")


def uplink_case():
    '''Uplink (uncorrected): the vacuum-optics + reciprocity Terms.'''
    print("=" * 70)
    print("UPLINK (uncorrected)")
    _line()
    scenario = SpaceScenario(
        ground=Terminal(aperture_m=0.4, wavelength_m=WAVELENGTH_M,
                        pointing_jitter_rad=1e-6,
                        transmitter=Transmitter(waist_m=0.15, power_dbm=40)),
        space=Terminal(aperture_m=0.05, wavelength_m=WAVELENGTH_M),
        direction="uplink", channel=Channel(altitude_m=600e3))
    geometry = CircularOrbit(600e3, 60.0)
    cn2 = default_cn2_profile(scenario.channel.site)

    t0 = time.time()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        bundle = run_fidelity2(scenario, geometry, n_trials=N_TRIALS, seed=SEED,
                               preset=PRESET, threader=THREADER, cn2_profile=cn2)
        f2 = uplink_budget(scenario, geometry, fidelity=2, wave=bundle,
                           cn2_profile=cn2)
    print(f"  fidelity 2 (wave optics, {N_TRIALS} snapshots, {PRESET}):")
    print(f"    standalone pointing Term kept: "
          f"{any(t.category == 'pointing' for t in f2.terms)} "
          f"(the reciprocity Term holds no jitter)")
    show_fidelity2(f2)
    print(f"  ({time.time() - t0:.1f} s)")
    print("")


def downlink_case():
    '''Downlink aperture: the vacuum-optics + turbulence Terms.'''
    print("=" * 70)
    print("DOWNLINK (aperture)")
    _line()
    scenario = SpaceScenario(
        ground=Terminal(aperture_m=0.2, wavelength_m=WAVELENGTH_M),
        space=Terminal(aperture_m=0.05, wavelength_m=WAVELENGTH_M,
                       transmitter=Transmitter(waist_m=0.035, power_dbm=30)),
        direction="downlink", channel=Channel(altitude_m=600e3))
    geometry = CircularOrbit(600e3, 30.0)
    cn2 = default_cn2_profile(scenario.channel.site)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        default = downlink_budget(scenario, geometry)             # fidelity 1
        dscint = next(t for t in default.terms if t.category == "turbulence")
    print(f"  fidelity 1 (analytic lognormal, default): scintillation "
          f"{dscint.mean_db:.4f} dB")

    t0 = time.time()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        bundle = run_fidelity2(scenario, geometry, n_trials=N_TRIALS, seed=SEED,
                               preset=PRESET, threader=THREADER, cn2_profile=cn2)
        f2 = downlink_budget(scenario, geometry, fidelity=2, wave=bundle)
    print(f"  fidelity 2 (wave optics, {N_TRIALS} snapshots, {PRESET}):")
    print("  the downlink mean turbulence penalty is small by nature "
          "(aperture-averaged);")
    print("  the vacuum-optics Term carries the 600 km geometric loss.")
    show_fidelity2(f2)
    print(f"  ({time.time() - t0:.1f} s)")
    print("")


def main():
    t_start = time.time()
    print(f"Fidelity-2 whole-path wiring ({PRESET} preset, {N_TRIALS} snapshots, "
          f"{THREADER.max_workers} threads)")
    print("")
    print("  At fidelity=2 the whole path is wave optics: a deterministic")
    print("  vacuum-optics Term (geometry + truncation + capture + vacuum")
    print("  coupling) and a stochastic turbulence Term. Only extinction")
    print("  (absorption) and pointing (mechanical jitter) stay analytic. Each")
    print("  budget ran ONE run_fidelity2 (a turbulent Monte Carlo plus one")
    print("  vacuum field solve). The fidelity=0/1 defaults are unchanged.")
    print("")
    print(f"  The fade availability is p={AVAILABILITY}. The tail adequacy of a")
    print("  quantile is a property of the sampler; the quantile warns when a run")
    print(f"  is too short. N_TRIALS={N_TRIALS} is adequate here.")
    print("")
    terrestrial_case()
    uplink_case()
    downlink_case()
    print(f"(total {time.time() - t_start:.1f} s)")


if __name__ == "__main__":
    main()
