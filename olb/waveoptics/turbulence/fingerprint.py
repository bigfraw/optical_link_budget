"""The content fingerprint that names one turbulent run.

WHY. A `Campaign` (`olb.waveoptics.turbulence.campaign`) keeps its trials on
disk over more than one session. Before it reads or grows a stored campaign, it
must know that the stored trials are the trials of THIS physics case. This
module gives that test: one SHA-256 of everything that changes a trial.

WHAT ENTERS THE KEY. The scenario hardware (its dataclass repr, which is
stable), a canonical geometry signature (a geometry object has no stable repr),
the preset, the base seed, the screen generator (the "olb" and "aotools"
generators draw DIFFERENT atmospheres for the same seed), the outer scale, the
subharmonic switch, the Cn2 inputs (a callable is fingerprinted by SAMPLING it
at a fixed set of heights), the block size, any caller-supplied grid and plan,
and `KEY_VERSION`. So two runs that share a key are the same physics, and a
change to any input gives a new key.

HISTORY. The key was born in the retired `cache.py` (the P4 scalar cache of
2026-08-29). `Campaign` replaced that cache on 2026-09-04, and the key moved
here. Nothing changed in it, so an existing campaign manifest still matches.

Sources:
- The seed contract and the trial body: olb.waveoptics.turbulence.run.
- The store that reads the key: olb.waveoptics.turbulence.campaign.
"""

import hashlib

import numpy as np

# Bump when the stored physics changes (a new trial scalar, a spectrum fix, a
# grid-sizer change). It enters the key, so an old store never feeds a new
# build. It keeps the value of the retired cache, so an existing campaign
# manifest still matches.
KEY_VERSION = 1


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
    geometry, the preset, the base seed, the screen generator and the key
    version, the outer scale, the subharmonic switch, the Cn2 profile, the
    block size, and any caller-supplied grid and plan. So two runs that share
    a key are the same physics, and a change to any input gives a new key.

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
        f"cache_version={KEY_VERSION}",
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


if __name__ == '__main__':
    from ...geometry import HorizontalPath
    from ...scenario import TerrestrialChannel, TerrestrialScenario
    from ...terminal import SMF, Terminal, Transmitter

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

    # ---- the key is stable, and it changes with each input ----
    common = dict(preset="rapid", screen_generator="olb", L0_m=np.inf,
                  subharmonics=True, hs=None, cn2_profile=None, block_size=2)
    k0 = cache_key(scn, geom, seed=7, **common)
    assert k0 == cache_key(scn, geom, seed=7, **common), "the key must be stable"
    assert k0 != cache_key(scn, geom, seed=8, **common), "seed must change it"
    assert k0 != cache_key(scn2, geom, seed=7, **common), "hardware changes it"
    assert k0 != cache_key(scn, geom, seed=7, **{**common, "preset": "standard"})
    assert k0 != cache_key(scn, geom, seed=7,
                           **{**common, "screen_generator": "aotools"})
    assert k0 != cache_key(scn, geom, seed=7, **{**common, "block_size": 4})
    assert k0 != cache_key(scn, HorizontalPath(2000.0), seed=7, **common)
    assert len(k0) == 64, len(k0)

    # ---- a Cn2 callable is fingerprinted by its values, not its identity ----
    kc = cache_key(scn, geom, seed=7, cn2=lambda h: 1e-15 * np.ones_like(h),
                   **common)
    assert kc == cache_key(scn, geom, seed=7,
                           cn2=lambda h: 1e-15 * np.ones_like(h), **common)
    assert kc != cache_key(scn, geom, seed=7,
                           cn2=lambda h: 2e-15 * np.ones_like(h), **common)
    assert kc != k0

    print(f"key {k0[:16]}... is stable; the seed, the hardware, the preset, "
          "the generator, the block size, the geometry and the Cn2 each "
          "change it")
    print("self-check passed")
