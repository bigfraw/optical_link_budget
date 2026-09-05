'''
Wire the fidelity-2 wave optics into the three budgets (the whole-path way).

Fidelity is a WHOLE-PATH choice (see the README fidelity ladder). At fidelity=2
the budget shows a DETERMINISTIC geometric loss plus a STOCHASTIC wave-optics
turbulence Term (the fade). Only the analytic extinction (molecular absorption,
never in the field sim) and pointing (mechanical jitter) Terms stay beside them.

The deterministic geometric loss has TWO forms, by link family:
  - SPACE (uplink, downlink): the ANALYTIC geometric Term (the default). A
    ground-space link is far field, so the analytic loss is exact. The
    wave-optics vacuum run is skipped: over the full slant range it is slow and
    grid-noise-limited (+/- 1 to 4 dB; see validation/vacuum_loss).
  - TERRESTRIAL: the wave-optics vacuum Term. The near-field turbulence penalty
    is turbulent / vacuum on the SAME flat grid, so the wave vacuum is the exact
    baseline. The vacuum solve runs on the CAMPAIGN grid, so the baseline stays
    exact.

This example touches three links, one selector each:
  - TERRESTRIAL SMF: terrestrial_budget(fidelity=2, wave=...). The default
    (fidelity=0) SMF Term is mean-only, so the default budget REFUSES a fade
    margin. Fidelity 2 gives one, and computes the FULL loss directly.
  - UPLINK (uncorrected): uplink_budget(fidelity=2, wave=...). The analytic
    geometric Term is the loss over 600 km; the turbulence Term is the
    reciprocity penalty.
  - DOWNLINK aperture: downlink_budget(fidelity=2, wave=...).

THE ONE RULE. The budget NEVER runs the split-step layer itself. Each link here
keeps its trials in a CAMPAIGN on disk
(olb.waveoptics.turbulence.Campaign, one directory for each link under
_campaigns/budget_wiring/). A Campaign IS a fidelity-2 wave record, so it goes
straight into the `wave` slot of the budget: the budget calls
olb.models.waveoptics.resolve_wave, which reads the stored trials. Run the
script two times. The second run computes NO trial, because the blocks are
already on disk.

Every run here uses the RAPID preset and a small trial count, so it finishes in a
couple of minutes; raise N_TRIALS and the preset for a real number. The tail
adequacy of a quantile is a property of the sampler
(olb.results.EmpiricalSampler): term.sampler.undersampled(p) answers it, and the
quantile warns when a run is too short.

Run from the repo root:
    python -m examples.waveoptics.budget_wiring
'''

import os
import time
import warnings

from olb import SMF, Terminal, Transmitter
from olb.geometry import CircularOrbit, HorizontalPath
from olb.scenario import (Channel, SpaceScenario, TerrestrialChannel,
                          TerrestrialScenario)
from olb.links.terrestrial import terrestrial_budget
from olb.links.uplink import uplink_budget
from olb.links.downlink import downlink_budget
from olb.turbulence.profiles import default_cn2_profile
from olb.waveoptics.turbulence import Campaign

WAVELENGTH_M = 1550e-9
PRESET = "rapid"       # a demonstration: RAPID keeps the example short
N_TRIALS = 200
BLOCK_SIZE = 50        # 200 trials in four blocks, one for each worker
WORKERS = 4            # one warm process pool for each campaign
SEED = 20260828
AVAILABILITY = 0.9     # the fade availability that every case reports

# One directory for each link. The store survives the process, so a second run
# of this script computes nothing.
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    "_campaigns", "budget_wiring")


def _line():
    print("-" * 70)


def show_fidelity2(budget):
    '''Print the deterministic geometric loss and the wave-optics turbulence Term.

    The deterministic geometric loss has two forms. A SPACE link uses the
    analytic geometric Term(s) by default (a far-field link; the wave vacuum run
    is slow and grid-noise-limited over the full slant range). A TERRESTRIAL link
    uses the wave-optics vacuum Term (the exact baseline for its turbulence
    penalty). Both are the non-stochastic loss outside extinction and pointing.
    '''
    turb = next(t for t in budget.terms if t.meta.get("model") == "waveoptics")
    geo_terms = [t for t in budget.terms
                 if not t.stochastic and t.category not in ("atmospheric", "pointing")]
    geo_db = sum(t.mean_db for t in geo_terms)
    is_wave = any(t.meta.get("model") == "waveoptics-vacuum" for t in geo_terms)
    label = "vacuum optics (wave)" if is_wave else "geometric optics (analytic)"
    fade = turb.quantile_db(AVAILABILITY)
    print(f"    {label:<32}{geo_db:8.3f} dB")
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
        campaign = Campaign(scenario, geometry,
                            os.path.join(ROOT, "terrestrial"), seed=SEED,
                            preset=PRESET, block_size=BLOCK_SIZE)
        campaign.run(N_TRIALS, workers=WORKERS, progress=True)
        f2 = terrestrial_budget(scenario, geometry, fidelity=2, wave=campaign)
    print(f"  fidelity 2 (wave optics, {N_TRIALS} snapshots, {PRESET}): "
          f"provides_fade {default.provides_fade} -> {f2.provides_fade}")
    show_fidelity2(f2)
    print(f"  ({time.time() - t0:.1f} s)")
    print("")


def uplink_case():
    '''Uplink (uncorrected): the analytic geometric + reciprocity Terms.'''
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
        campaign = Campaign(scenario, geometry, os.path.join(ROOT, "uplink"),
                            seed=SEED, preset=PRESET, block_size=BLOCK_SIZE,
                            cn2_profile=cn2)
        campaign.run(N_TRIALS, workers=WORKERS, progress=True)
        f2 = uplink_budget(scenario, geometry, fidelity=2, wave=campaign,
                           cn2_profile=cn2)
    print(f"  fidelity 2 (wave optics, {N_TRIALS} snapshots, {PRESET}):")
    print(f"    standalone pointing Term kept: "
          f"{any(t.category == 'pointing' for t in f2.terms)} "
          f"(the reciprocity Term holds no jitter)")
    show_fidelity2(f2)
    print(f"  ({time.time() - t0:.1f} s)")
    print("")


def downlink_case():
    '''Downlink aperture: the analytic geometric + turbulence Terms.'''
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
        campaign = Campaign(scenario, geometry, os.path.join(ROOT, "downlink"),
                            seed=SEED, preset=PRESET, block_size=BLOCK_SIZE,
                            cn2_profile=cn2)
        campaign.run(N_TRIALS, workers=WORKERS, progress=True)
        f2 = downlink_budget(scenario, geometry, fidelity=2, wave=campaign)
    print(f"  fidelity 2 (wave optics, {N_TRIALS} snapshots, {PRESET}):")
    print("  the downlink mean turbulence penalty is small by nature "
          "(aperture-averaged);")
    print("  the analytic geometric Term carries the 600 km geometric loss.")
    show_fidelity2(f2)
    print(f"  ({time.time() - t0:.1f} s)")
    print("")


def main():
    t_start = time.time()
    print(f"Fidelity-2 whole-path wiring ({PRESET} preset, {N_TRIALS} snapshots, "
          f"{WORKERS} workers)")
    print("")
    print("  At fidelity=2 the whole path is wave optics: a deterministic")
    print("  vacuum-optics Term (geometry + truncation + capture + vacuum")
    print("  coupling) and a stochastic turbulence Term. Only extinction")
    print("  (absorption) and pointing (mechanical jitter) stay analytic. Each")
    print("  budget reads ONE Campaign of trials from disk (wave=campaign).")
    print("  A second run of this script computes no trial. The fidelity=0/1")
    print("  defaults are unchanged.")
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
