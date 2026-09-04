"""An opt-in, extendable disk cache for turbulent wave-optics runs.

WHY. A fidelity-2 trial is expensive, and a verification campaign runs
thousands of them. The scalar result of a trial is tiny (five numbers), so a
whole run stores as a small JSON. This module keeps those scalars on disk and
in process, so a campaign runs each atmosphere ONE time and reads it back for
every budget, every quantile, and every re-run.

WHAT IT CACHES. This is cache LEVEL 3 of the speed plan
(docs/waveoptics-efficiency-plan.md): whole `TurbWaveResult` runs, as scalars.
The per-grid filter and the subharmonic basis are cache level 1 (inside
`olb.waveoptics.turbulence.screens.ScreenFactory`). The vacuum baseline is
cache level 2; a repeated in-process request hits `_SESSION_MEMO` here and
skips the whole call, the baseline with it (see the note below).

THE EXTENDABLE DESIGN. A campaign asks for more trials as it goes. This cache
computes ONLY the new ones. It stores the run in fixed BLOCKS of `block_size`
trials. Block b is a self-contained `propagate_turbulent_scenario` run of
`block_size` trials, seeded by a value SPAWNED from the base seed and the block
index (`_block_entropy`). So:
  - a stored run of 4 blocks serves any request for those trials by a slice;
  - a request past the stored count computes only the missing blocks;
  - block b is the SAME trials whatever the request order, because its seed
    derives from (base_seed, b) alone.

The blocks are independent, identically distributed snapshots, which is what
the empirical Term reducer wants (`olb.models.waveoptics`). They are NOT the
trials of ONE native `propagate_turbulent_scenario(seed=base_seed)` run: that
runner seeds trial k off (base_seed, k), and this cache seeds each block off a
sub-seed instead. This cache trades that bit-identity for a real compute saving
on a grow, the same way the "olb" and "aotools" generators trade bit-identity
for speed. A single-seed tail extension is now POSSIBLE:
`olb.waveoptics.turbulence.run.propagate_turbulent_scenario` HAS a
`start_index` argument, so a caller computes the trials 200..499 of one seed
alone. This module does NOT use it: the block sub-seeds stay, and the cache
behaviour does not change. A later campaign module uses `start_index`.

WHY SCALARS SUFFICE. The Term reducer (`olb.models.waveoptics`,
`waveoptics_turbulence_term`) reads ONLY the per-trial scalars
(`collected_power`, `smf_eta`, `eta_turb`, `mmf_eta`), the preset, and the seed
entropy. It never reads the grid, the plan, or the field. So the five stored
scalars per trial reconstruct everything a budget needs. A cache-only load
returns `grid=None` and `plan=None`, which the reducer accepts.

THE VACUUM BASELINE (level 2). The space vacuum baseline is recomputed inside
each `propagate_turbulent_scenario` call (see
`olb.waveoptics.turbulence.run`). A wrapper cannot inject a precomputed
baseline without editing that runner, so this module does not try. Instead:
  - a disk HIT skips the whole call, so the baseline is not recomputed at all;
  - a repeated in-process request hits `_SESSION_MEMO` and returns at once;
  - the block calls size the grid and plan ONE time and share them, so only
    the (cheap) baseline split-step repeats per block, not the grid sizing.

OWNER-GATED. This cache is OPT-IN and OFF by default. No budget calls it. The
budgets keep their own precompute pattern (`olb.models.waveoptics.run_fidelity2`
builds the `Fidelity2Bundle` directly). Wiring this cache into a budget or a
campaign default is an owner decision.

Sources:
- The seed contract and the trial body: olb.waveoptics.turbulence.run.
- The cached screen generator (level 1): olb.waveoptics.turbulence.screens.
- The budget precompute pattern: olb.models.waveoptics.run_fidelity2.
"""

import hashlib
import json
import os
import tempfile

import numpy as np

from .run import TurbTrial, TurbWaveResult
from .run import propagate_turbulent_scenario as _run_scenario
from .sampling import turbulent_grid

# Bump when the stored physics changes (a new trial scalar, a spectrum fix, a
# grid-sizer change). It enters the key, so an old file never feeds a new
# build. The screen generator name enters the key too, because "olb" and
# "aotools" draw different atmospheres (see run._screen_builder).
CACHE_VERSION = 1

# The grow granularity. A block is one self-contained run. Smaller blocks give
# finer extend steps; larger blocks call the runner (and its vacuum baseline)
# fewer times. 50 is a plain campaign default.
DEFAULT_BLOCK_SIZE = 50

# The in-process store of computed trials, keyed by the cache key. It gives the
# level-2 saving: a repeated request in one session returns at once, and skips
# the vacuum baseline the runner would recompute.
_SESSION_MEMO = {}


def default_cache_dir():
    """Give the default cache directory, a scratch path outside the repo.

    It sits under the system temp directory, so it is never committed. Override
    it with the `cache_dir` argument of `cached_propagate_turbulent_scenario`,
    or with the OLB_WAVEOPTICS_CACHE environment variable.

    Returns:
        The directory path, as a string.
    """
    env = os.environ.get("OLB_WAVEOPTICS_CACHE")
    if env:
        return env
    return os.path.join(tempfile.gettempdir(), "olb_waveoptics_cache")


def _array_sha(a):
    """Give a short hash of a numpy array, or 'none' when it is None."""
    if a is None:
        return "none"
    arr = np.ascontiguousarray(np.asarray(a, dtype=float))
    return hashlib.sha256(arr.tobytes()).hexdigest()[:16]


def _geometry_signature(geometry):
    """Give a stable string for the geometry.

    A geometry object has no stable repr (its default repr holds a memory
    address), so this reads the physics fields that set a trial: the slant
    range, and the elevation, altitude, or path length when present.

    Args:
        geometry: a CircularOrbit, a TLEPass, or a HorizontalPath.

    Returns:
        A canonical string.
    """
    parts = [type(geometry).__name__]
    for attr in ("slant_range_m", "elevation_deg", "altitude_m",
                 "path_length_m"):
        val = getattr(geometry, attr, None)
        if val is not None:
            parts.append(f"{attr}={np.asarray(val, dtype=float).tolist()}")
    return "|".join(parts)


def _cn2_fingerprint(cn2, h_top_m):
    """Give a stable string that identifies a Cn2 callable, or "none".

    A callable has no stable repr, so the key SAMPLES it: it evaluates cn2 at a
    fixed set of heights and hashes the values with the integration top. Two
    callables that agree on the profile give the same key; a changed profile
    gives a new key. It never enters the physics; it only names the run.
    """
    if cn2 is None:
        return "none"
    h = np.linspace(0.0, float(h_top_m) if h_top_m is not None else 20e3, 64)
    return _array_sha(np.concatenate((np.asarray(cn2(h), float), [h[-1]])))


def cache_key(scenario, geometry, *, preset, seed, screen_generator,
              L0_m, subharmonics, hs, cn2_profile, block_size,
              cn2=None, h_top_m=None, grid=None, plan=None):
    """Give the content hash that names a stored run.

    The key holds EVERYTHING that changes a trial: the scenario hardware, the
    geometry, the preset, the base seed, the screen generator and its version,
    the outer scale, the subharmonic switch, the Cn2 profile, the block size,
    and any caller-supplied grid and plan. So two runs that share a key are the
    same physics, and a change to any input gives a new file.

    Args:
        scenario:         a SpaceScenario or a TerrestrialScenario.
        geometry:         the link geometry (one range).
        preset:           the preset name (a string).
        seed:             the integer base seed.
        screen_generator: "olb" or "aotools".
        L0_m:             the outer scale, in m.
        subharmonics:     the subharmonic switch.
        hs, cn2_profile:  the height grid and the zenith Cn2 profile, or None.
        cn2:              the continuous Cn2 callable, or None. Fingerprinted by
                          sampling, see _cn2_fingerprint.
        h_top_m:          the atmosphere top for the continuous integral, or None.
        block_size:       the block size.
        grid, plan:       an optional caller-supplied grid and plan.

    Returns:
        A 64-character hex string.
    """
    preset_name = preset if isinstance(preset, str) else getattr(
        preset, "name", repr(preset))
    blob = "\n".join([
        f"cache_version={CACHE_VERSION}",
        f"scenario={scenario!r}",
        f"geometry={_geometry_signature(geometry)}",
        f"preset={preset_name}",
        f"seed={int(seed)}",
        f"screen_generator={screen_generator}",
        f"L0_m={float(L0_m)!r}",
        f"subharmonics={bool(subharmonics)}",
        f"cn2_fp={_cn2_fingerprint(cn2, h_top_m)}",
        f"hs_sha={_array_sha(hs)}",
        f"cn2_sha={_array_sha(cn2_profile)}",
        f"h_top_m={float(h_top_m) if h_top_m is not None else None!r}",
        f"block_size={int(block_size)}",
        f"grid={grid!r}",
        f"plan={plan!r}",
    ])
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()


def _block_entropy(base_seed, block_index):
    """Give the deterministic integer seed of one block.

    It spawns from the base seed and the block index, so block b is the same
    run whatever the request order. See the numpy SeedSequence documentation.
    """
    ss = np.random.SeedSequence(entropy=int(base_seed),
                                spawn_key=(CACHE_VERSION, int(block_index)))
    return int(ss.generate_state(1)[0])


def _trial_to_row(t):
    """Pack one TurbTrial into a small JSON row (five scalars)."""
    return [t.collected_power, t.smf_eta, t.eta_turb, t.mmf_eta, t.wall_time_s]


def _row_to_trial(row, block_entropy, index):
    """Rebuild one TurbTrial from a stored row. The seed_key is derived."""
    collected, smf_eta, eta_turb, mmf_eta, wall = row
    return TurbTrial(collected_power=collected, smf_eta=smf_eta,
                     eta_turb=eta_turb, seed_key=(block_entropy, index),
                     wall_time_s=wall, mmf_eta=mmf_eta)


def _path_for_key(cache_dir, key):
    """Give the file path for a key, in a two-character sub-directory."""
    return os.path.join(cache_dir, key[:2], key + ".json")


def _load_blocks(path):
    """Load the stored block list from disk, or [] when the file is absent."""
    if not os.path.exists(path):
        return []
    with open(path, "r", encoding="utf-8") as fh:
        data = json.load(fh)
    return data.get("blocks", [])


def _write_blocks(path, blocks, meta):
    """Write the block list and the meta to disk (atomic replace)."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    data = {"cache_version": CACHE_VERSION, "meta": meta, "blocks": blocks}
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8") as fh:
        json.dump(data, fh)
    os.replace(tmp, path)


def cached_propagate_turbulent_scenario(
        scenario, geometry, *, n_trials=1, seed,
        preset="standard", cache_dir=None, block_size=DEFAULT_BLOCK_SIZE,
        screen_generator="olb", L0_m=np.inf, subharmonics=True,
        cn2=None, hs=None, cn2_profile=None, h_top_m=None,
        grid=None, plan=None, threader=None,
        refresh=False, store=True):
    """Load a turbulent run from the cache, or run only the missing blocks.

    This is the one entry point. It is a load-or-run WRAPPER around the public
    `propagate_turbulent_scenario`. It NEVER edits the runner and it changes no
    budget. It is opt-in: a caller reaches for it on purpose.

    THE FLOW. It names the run with a content `cache_key`, then it grows the
    stored run to at least `n_trials` trials in fixed blocks of `block_size`. A
    block that already sits in the in-process memo or on disk is reused; a
    missing block runs once and is stored. The result is a `TurbWaveResult`
    with the first `n_trials` trials, ready for
    `olb.models.waveoptics.waveoptics_turbulence_term`.

    THE SEED IS REQUIRED. Extendability needs a stable base seed, so `seed`
    must be an integer. A None or a Generator raises: their entropy is not
    repeatable, so a later call could not find the stored run.

    Args:
        scenario:         a SpaceScenario or a TerrestrialScenario.
        geometry:         the link geometry. It must give ONE range.
        n_trials:         the number of trials the caller needs.
        seed:             the integer base seed (required).
        preset:           the preset name.
        cache_dir:        the cache directory. None uses default_cache_dir().
        block_size:       the grow granularity, in trials.
        screen_generator: "olb" (the default) or "aotools". It enters the key.
        L0_m:             the outer scale, in m.
        subharmonics:     the subharmonic switch.
        cn2:              an optional continuous Cn2 callable cn2(h) (space).
                          None (with no hs/cn2_profile) integrates the site
                          profile: the continuous default. The key fingerprints
                          it at a fixed set of heights.
        hs, cn2_profile:  the height grid and the zenith Cn2 profile of the
                          LEGACY array planner (space).
        h_top_m:          the atmosphere top for the continuous integral (space).
        grid, plan:       an optional grid and plan. Give both or neither. When
                          None, the wrapper sizes them once and shares them
                          across the block runs.
        threader:         an optional Threader, passed to each block run.
        refresh:          True discards the stored run and recomputes it.
        store:            True writes new blocks to disk. False keeps them in
                          the session memo only.

    Returns:
        A TurbWaveResult with `n_trials` trials. Its grid and plan are the
        sized objects when the call computed a block, else None.

    Raises:
        ValueError: the seed is not an integer, or n_trials is not positive.
    """
    if seed is None or isinstance(seed, np.random.Generator):
        raise ValueError(
            "cached_propagate_turbulent_scenario needs an integer seed. The "
            "cache keys the stored run on the seed, so a None or a Generator "
            "(a non-repeatable entropy) cannot be found again. Pass seed=<int>.")
    n_trials = int(n_trials)
    if n_trials <= 0:
        raise ValueError(f"n_trials must be positive, not {n_trials}.")

    base_seed = int(seed)
    block_size = int(block_size)
    if cache_dir is None:
        cache_dir = default_cache_dir()

    key = cache_key(scenario, geometry, preset=preset, seed=base_seed,
                    screen_generator=screen_generator, L0_m=L0_m,
                    subharmonics=subharmonics, cn2=cn2, hs=hs,
                    cn2_profile=cn2_profile, h_top_m=h_top_m,
                    block_size=block_size, grid=grid, plan=plan)
    path = _path_for_key(cache_dir, key)

    if refresh:
        _SESSION_MEMO.pop(key, None)
        if os.path.exists(path):
            os.remove(path)

    # ---- gather the trials this session already holds ----
    trials = list(_SESSION_MEMO.get(key, []))
    if not trials:
        # Disk holds JSON rows. Rebuild them into TurbTrial objects.
        for block in _load_blocks(path):
            ent = block["entropy"]
            for i, row in enumerate(block["trials"]):
                trials.append(_row_to_trial(row, ent, i))

    n_blocks_have = len(trials) // block_size
    n_blocks_need = -(-n_trials // block_size)          # ceil division.

    # ---- run only the missing blocks ----
    sized_grid, sized_plan = grid, plan
    if n_blocks_need > n_blocks_have and sized_grid is None:
        # Size the grid and the plan ONE time, then share them across the
        # block runs. This is the level-2 in-process reuse of the setup.
        sized_grid, sized_plan, _ = turbulent_grid(
            scenario, geometry, preset=preset, cn2=cn2, hs=hs,
            cn2_profile=cn2_profile, h_top_m=h_top_m, L0_m=L0_m)

    new_blocks = []
    for b in range(n_blocks_have, n_blocks_need):
        ent = _block_entropy(base_seed, b)
        res = _run_scenario(
            scenario, geometry, n_trials=block_size, seed=ent, preset=preset,
            grid=sized_grid, plan=sized_plan, cn2=cn2, hs=hs,
            cn2_profile=cn2_profile, h_top_m=h_top_m, L0_m=L0_m,
            subharmonics=subharmonics, threader=threader,
            screen_generator=screen_generator)
        new_blocks.append({"index": b, "entropy": ent,
                           "trials": [_trial_to_row(t) for t in res.trials]})
        trials.extend(res.trials)

    # ---- persist the grown run ----
    if new_blocks:
        _SESSION_MEMO[key] = trials
        if store:
            all_blocks = _load_blocks(path)
            all_blocks.extend(new_blocks)
            meta = {
                "scenario": repr(scenario),
                "geometry": _geometry_signature(geometry),
                "preset": preset if isinstance(preset, str) else getattr(
                    preset, "name", repr(preset)),
                "base_seed": base_seed,
                "block_size": block_size,
                "screen_generator": screen_generator,
                "L0_m": None if not np.isfinite(L0_m) else float(L0_m),
                "subharmonics": bool(subharmonics),
            }
            _write_blocks(path, all_blocks, meta)
    else:
        # A pure hit still populates the memo, so the next call skips the disk.
        _SESSION_MEMO[key] = trials

    return TurbWaveResult(
        trials=trials[:n_trials], grid=sized_grid, plan=sized_plan,
        report=None,
        preset=preset if isinstance(preset, str) else getattr(
            preset, "name", repr(preset)),
        seed_entropy=base_seed)


if __name__ == '__main__':
    import shutil
    import time
    import warnings

    from ...geometry import HorizontalPath
    from ...scenario import TerrestrialChannel, TerrestrialScenario
    from ...terminal import SMF, Terminal, Transmitter

    # The "olb" generator is pure numpy, so this self-check needs no aotools.
    scn = TerrestrialScenario(
        near=Terminal(aperture_m=0.2, wavelength_m=1550e-9,
                      transmitter=Transmitter(waist_m=0.02)),
        far=Terminal(aperture_m=0.2, wavelength_m=1550e-9, detector=SMF()),
        channel=TerrestrialChannel(path_length_m=1000.0, cn2=5e-15))
    geom = HorizontalPath(1000.0)
    scn2 = TerrestrialScenario(                       # a different aperture.
        near=Terminal(aperture_m=0.2, wavelength_m=1550e-9,
                      transmitter=Transmitter(waist_m=0.02)),
        far=Terminal(aperture_m=0.25, wavelength_m=1550e-9, detector=SMF()),
        channel=TerrestrialChannel(path_length_m=1000.0, cn2=5e-15))

    # ---- 1. the key is stable, and it changes with the inputs ----
    common = dict(preset="rapid", screen_generator="olb", L0_m=np.inf,
                  subharmonics=True, hs=None, cn2_profile=None, block_size=2)
    k0 = cache_key(scn, geom, seed=7, **common)
    assert k0 == cache_key(scn, geom, seed=7, **common), "the key must be stable"
    assert k0 != cache_key(scn, geom, seed=8, **common), "seed must change it"
    assert k0 != cache_key(scn2, geom, seed=7, **common), "hardware changes it"
    assert k0 != cache_key(scn, geom, seed=7, **{**common, "preset": "standard"})
    assert k0 != cache_key(scn, geom, seed=7,
                           **{**common, "screen_generator": "aotools"})
    assert len(k0) == 64, len(k0)

    # ---- 2. the block entropy is deterministic and per-block distinct ----
    assert _block_entropy(7, 0) == _block_entropy(7, 0)
    assert _block_entropy(7, 0) != _block_entropy(7, 1)
    assert _block_entropy(7, 0) != _block_entropy(8, 0)

    # ---- 3. the seed guard ----
    tmp = tempfile.mkdtemp(prefix="olb_cache_selfcheck_")
    try:
        for bad in (None, np.random.default_rng(0)):
            try:
                cached_propagate_turbulent_scenario(
                    scn, geom, n_trials=2, seed=bad, preset="rapid",
                    cache_dir=tmp, block_size=2, screen_generator="olb")
                raise AssertionError("a non-integer seed must raise")
            except ValueError as exc:
                assert "integer seed" in str(exc), str(exc)

        # A counting wrapper measures how many trials the runner actually
        # computes. This proves the "compute only the new blocks" behaviour.
        # Patch THIS module (running as __main__), not a re-import, so the
        # patched name and the shared memo are the ones the entry point reads.
        import sys
        C = sys.modules[__name__]
        orig = C._run_scenario
        stat = {"calls": 0, "trials": 0}

        def counting(*a, **k):
            stat["calls"] += 1
            stat["trials"] += int(k.get("n_trials", 0))
            return orig(*a, **k)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            # ---- 4. a cold run computes every block ----
            C._run_scenario = counting
            stat["calls"] = stat["trials"] = 0
            t0 = time.time()
            r4 = cached_propagate_turbulent_scenario(
                scn, geom, n_trials=4, seed=7, preset="rapid", cache_dir=tmp,
                block_size=2, screen_generator="olb")
            cold_s = time.time() - t0
            assert len(r4.trials) == 4, len(r4.trials)
            assert stat["trials"] == 4, stat        # two blocks of two.
            assert stat["calls"] == 2, stat
            assert os.path.exists(_path_for_key(tmp, k0)), "the file must exist"
            first_four = [t.collected_power for t in r4.trials]

            # ---- 5. a warm in-process hit computes nothing ----
            stat["calls"] = stat["trials"] = 0
            t0 = time.time()
            r4b = cached_propagate_turbulent_scenario(
                scn, geom, n_trials=4, seed=7, preset="rapid", cache_dir=tmp,
                block_size=2, screen_generator="olb")
            warm_s = time.time() - t0
            assert stat["trials"] == 0, stat        # a pure memo hit.
            assert [t.collected_power for t in r4b.trials] == first_four

            # ---- 5b. a warm DISK hit (fresh session) computes nothing ----
            C._SESSION_MEMO.clear()
            stat["calls"] = stat["trials"] = 0
            r4c = cached_propagate_turbulent_scenario(
                scn, geom, n_trials=4, seed=7, preset="rapid", cache_dir=tmp,
                block_size=2, screen_generator="olb")
            assert stat["trials"] == 0, "a disk hit must not recompute"
            assert [t.collected_power for t in r4c.trials] == first_four

            # ---- 6. a grow computes ONLY the new block ----
            stat["calls"] = stat["trials"] = 0
            r6 = cached_propagate_turbulent_scenario(
                scn, geom, n_trials=6, seed=7, preset="rapid", cache_dir=tmp,
                block_size=2, screen_generator="olb")
            assert len(r6.trials) == 6, len(r6.trials)
            assert stat["trials"] == 2, stat        # one new block of two.
            assert stat["calls"] == 1, stat
            # The first four trials are byte-for-byte the cold run's trials.
            assert [t.collected_power for t in r6.trials[:4]] == first_four

            # ---- 7. a smaller request slices, and computes nothing ----
            stat["calls"] = stat["trials"] = 0
            r3 = cached_propagate_turbulent_scenario(
                scn, geom, n_trials=3, seed=7, preset="rapid", cache_dir=tmp,
                block_size=2, screen_generator="olb")
            assert len(r3.trials) == 3, len(r3.trials)
            assert stat["trials"] == 0, stat
            assert [t.collected_power for t in r3.trials] == first_four[:3]

            # ---- 8. refresh recomputes, and the result is reproducible ----
            stat["calls"] = stat["trials"] = 0
            r4d = cached_propagate_turbulent_scenario(
                scn, geom, n_trials=4, seed=7, preset="rapid", cache_dir=tmp,
                block_size=2, screen_generator="olb", refresh=True)
            assert stat["trials"] == 4, stat        # a full recompute.
            # The block seeds derive from (base_seed, block), so the recompute
            # reproduces the original trials exactly.
            assert [t.collected_power for t in r4d.trials] == first_four

            # ---- 9. a different seed gives a different run and a new file ----
            r_other = cached_propagate_turbulent_scenario(
                scn, geom, n_trials=4, seed=8, preset="rapid", cache_dir=tmp,
                block_size=2, screen_generator="olb")
            assert ([t.collected_power for t in r_other.trials] != first_four)

        C._run_scenario = orig
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    print("terrestrial 1 km rapid, olb generator, block_size 2:")
    print(f"  cold run (4 trials)     {cold_s * 1e3:8.1f} ms")
    print(f"  warm memo hit           {warm_s * 1e3:8.3f} ms")
    print("  the hit returns from memory; it runs no trial and no baseline.")
    print("  grow 4 -> 6 computed only the new block (2 trials).")
    print("  refresh reproduced the original trials exactly.")
    print("self-check passed")
