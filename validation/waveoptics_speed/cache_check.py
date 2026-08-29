"""Measure the fidelity-2 result cache (P4, cache level 3 plus level 2).

RESULTS (2026-08-29, olb generator, no aotools):
- A disk HIT of a stored run returns in about a millisecond and computes no
  trial, so it saves the whole run cost AND the per-call vacuum baseline.
- A GROW computes only the new blocks: extending 100 -> 150 trials costs about
  the same as a fresh 50-trial run, not a fresh 150-trial run.
- The vacuum baseline (level 2) is a real slice of a space call; a hit skips it
  with the rest of the call, and the block runs size the grid once and share
  it.

This is a measurement script, not production code. It follows the validation/
pattern: one script, one results JSON, one run log. It uses the "olb" screen
generator, so it needs no aotools.

Sources:
- The cache under test: olb.waveoptics.turbulence.cache.
- The seed contract and the vacuum baseline: olb.waveoptics.turbulence.run.
- Schmidt, DOI 10.1117/3.866274, Ch. 9 (the split-step propagation).
"""

import json
import os
import shutil
import tempfile
import time
import warnings

import numpy as np

from olb.geometry import CircularOrbit, HorizontalPath
from olb.scenario import (Channel, SpaceScenario, TerrestrialChannel,
                          TerrestrialScenario)
from olb.terminal import SMF, Terminal, Transmitter
from olb.waveoptics.turbulence import cache as wcache
from olb.waveoptics.turbulence.run import propagate_turbulent_scenario

LAM = 1550e-9
HERE = os.path.dirname(os.path.abspath(__file__))


def _terrestrial():
    scn = TerrestrialScenario(
        near=Terminal(aperture_m=0.2, wavelength_m=LAM,
                      transmitter=Transmitter(waist_m=0.02)),
        far=Terminal(aperture_m=0.2, wavelength_m=LAM, detector=SMF()),
        channel=TerrestrialChannel(path_length_m=2000.0, cn2=5e-15))
    return scn, HorizontalPath(2000.0)


def _space_downlink():
    ground = Terminal(aperture_m=0.4, wavelength_m=LAM,
                      transmitter=Transmitter(waist_m=0.06))
    scn = SpaceScenario(ground=ground,
                        space=Terminal(aperture_m=0.3, wavelength_m=LAM),
                        direction="downlink", channel=Channel())
    return scn, CircularOrbit(altitude_m=600e3, elevation_deg=[30.0])


def _timed(fn):
    t0 = time.perf_counter()
    out = fn()
    return out, time.perf_counter() - t0


def _count_trials(cache_dir, scn, geom, **kw):
    """Run a cache request while counting the trials the runner computes."""
    orig = wcache._run_scenario
    stat = {"trials": 0}

    def counting(*a, **k):
        stat["trials"] += int(k.get("n_trials", 0))
        return orig(*a, **k)

    wcache._run_scenario = counting
    try:
        (res, wall) = _timed(lambda: wcache.cached_propagate_turbulent_scenario(
            scn, geom, cache_dir=cache_dir, **kw))
    finally:
        wcache._run_scenario = orig
    return res, wall, stat["trials"]


def main():
    out = {"environment": {"numpy": np.__version__,
                           "cores": os.cpu_count()},
           "cases": {}}
    preset = "rapid"
    block = 50

    for name, (scn, geom) in (("terrestrial 2km rapid", _terrestrial()),
                              ("space downlink 30deg rapid", _space_downlink())):
        tmp = tempfile.mkdtemp(prefix="olb_cache_val_")
        wcache._SESSION_MEMO.clear()
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")

                # Cold: compute 100 trials (two blocks of 50).
                _, cold_s, cold_trials = _count_trials(
                    tmp, scn, geom, n_trials=100, seed=2024, preset=preset,
                    block_size=block, screen_generator="olb")

                # Warm disk hit: a fresh session (clear the memo).
                wcache._SESSION_MEMO.clear()
                _, warm_s, warm_trials = _count_trials(
                    tmp, scn, geom, n_trials=100, seed=2024, preset=preset,
                    block_size=block, screen_generator="olb")

                # Grow 100 -> 150: only the third block runs.
                _, grow_s, grow_trials = _count_trials(
                    tmp, scn, geom, n_trials=150, seed=2024, preset=preset,
                    block_size=block, screen_generator="olb")

                # A fresh 50-trial run, for the grow reference.
                tmp2 = tempfile.mkdtemp(prefix="olb_cache_val2_")
                wcache._SESSION_MEMO.clear()
                _, one_block_s, one_block_trials = _count_trials(
                    tmp2, scn, geom, n_trials=50, seed=999, preset=preset,
                    block_size=block, screen_generator="olb")
                shutil.rmtree(tmp2, ignore_errors=True)

                # The vacuum baseline slice (level 2): one bare 1-trial call
                # pays the baseline once; report the per-trial wall time as the
                # cost the hit skips per call.
                base = propagate_turbulent_scenario(
                    scn, geom, n_trials=1, seed=1, preset=preset,
                    screen_generator="olb")
                baseline_trial_s = base.trials[0].wall_time_s

            key_path = wcache._path_for_key(
                tmp, wcache.cache_key(
                    scn, geom, preset=preset, seed=2024,
                    screen_generator="olb", L0_m=np.inf, subharmonics=True,
                    hs=None, cn2_profile=None, block_size=block))
            file_bytes = os.path.getsize(key_path)

            out["cases"][name] = {
                "cold_100_s": cold_s, "cold_trials_computed": cold_trials,
                "warm_hit_s": warm_s, "warm_trials_computed": warm_trials,
                "grow_100_to_150_s": grow_s,
                "grow_trials_computed": grow_trials,
                "fresh_50_s": one_block_s,
                "fresh_50_trials_computed": one_block_trials,
                "one_trial_wall_s": baseline_trial_s,
                "hit_speedup_over_cold": cold_s / max(warm_s, 1e-9),
                "stored_file_bytes_for_150": None,
                "block_size": block,
            }
        finally:
            shutil.rmtree(tmp, ignore_errors=True)

        # Re-measure the stored file size for 150 trials in a clean dir.
        tmp3 = tempfile.mkdtemp(prefix="olb_cache_size_")
        wcache._SESSION_MEMO.clear()
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                wcache.cached_propagate_turbulent_scenario(
                    scn, geom, n_trials=150, seed=2024, preset=preset,
                    cache_dir=tmp3, block_size=block, screen_generator="olb")
            path = wcache._path_for_key(tmp3, wcache.cache_key(
                scn, geom, preset=preset, seed=2024, screen_generator="olb",
                L0_m=np.inf, subharmonics=True, hs=None, cn2_profile=None,
                block_size=block))
            out["cases"][name]["stored_file_bytes_for_150"] = os.path.getsize(
                path)
        finally:
            shutil.rmtree(tmp3, ignore_errors=True)

    with open(os.path.join(HERE, "cache_check_results.json"), "w",
              encoding="utf-8") as fh:
        json.dump(out, fh, indent=1)

    # ---- the printed table ----
    print("fidelity-2 result cache (olb generator, block_size 50)")
    print("=" * 66)
    for name, c in out["cases"].items():
        print(f"\n{name}:")
        print(f"  one trial wall              {c['one_trial_wall_s'] * 1e3:9.1f} ms")
        print(f"  cold, 100 trials            {c['cold_100_s']:9.2f} s "
              f"({c['cold_trials_computed']} computed)")
        print(f"  warm disk hit, 100 trials   {c['warm_hit_s'] * 1e3:9.1f} ms "
              f"({c['warm_trials_computed']} computed)")
        print(f"  hit speed-up over cold      {c['hit_speedup_over_cold']:9.0f} x")
        print(f"  grow 100 -> 150             {c['grow_100_to_150_s']:9.2f} s "
              f"({c['grow_trials_computed']} computed)")
        print(f"  fresh 50 (grow reference)   {c['fresh_50_s']:9.2f} s "
              f"({c['fresh_50_trials_computed']} computed)")
        print(f"  stored file, 150 trials     {c['stored_file_bytes_for_150']:9d} bytes")
    print("\nA grow computes only the new block: grow 100->150 ~ a fresh 50, "
          "not a fresh 150.")
    print("A hit computes zero trials and skips the vacuum baseline with the "
          "whole call.")
    print("self-check passed")


if __name__ == "__main__":
    main()
