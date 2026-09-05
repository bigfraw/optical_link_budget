"""A campaign of turbulent trials on disk, stored as blocks.

WHY. A fade statistic needs thousands of snapshots. One trial is expensive, so
a campaign must survive a stopped process, grow later, and give its fields back
to a NEW detector with no new propagation. This module is that store.

WHAT IT IS. A `Campaign` names one physics case (one scenario, one geometry,
one grid, one screen plan, one seed) and it keeps its trials in fixed BLOCKS of
`block_size` trials, one `.npz` file for each block. Block b holds the trials
b*block_size .. (b+1)*block_size - 1 of ONE native run: the runner seeds trial k
off (seed, k), so the blocks are bit-identical SLICES of a single long run. A
campaign therefore computes its blocks in any order, on any number of processes,
and the concatenation is the native run. See
`olb.waveoptics.turbulence.run.propagate_turbulent_scenario` and its
`start_index` argument.

HISTORY. This module REPLACED the P4 scalar cache (`cache.py`, retired
2026-09-04). That cache seeded each block from a SUB-SEED, so its blocks were
not the trials of one native run, and it stored no field. This module keeps the
native seeding and it stores the field. The content fingerprint of that cache
lives on in `olb.waveoptics.turbulence.fingerprint`.

THE STORED FIELD. Each trial stores the receive-plane field on a disc of the
radius `patch_radius_m`, BEFORE the receive-aperture clip. Store the field at
the LARGEST aperture of the family one time. Then a smaller receive aperture, a
central obscuration, a different detector, a different focal length and a
different defocus are all a POST-HOC crop of that stored field (`recouple` and
`recollect`), with no new propagation. This is exact for a SPACE downlink,
because the propagated slab does not read the receive terminal at all: the
input is a plane wave, and the receive terminal enters only at the clip. A
TERRESTRIAL path does read the TRANSMIT terminal, so only the receive side is
free.

`sizing_aperture_m` serves that plan: the grid is sized on a copy of the
scenario whose CLIP terminal carries that larger aperture, and the trials
then run on THAT grid with the ORIGINAL scenario. So one campaign covers every
receive aperture up to the sizing aperture. The clip terminal is the ground
terminal of a space scenario in EVERY direction (an uplink reads the ground
field through reciprocity), and the receive terminal of a terrestrial one
(`run.clip_terminal`). The default patch radius follows the same terminal, so
the stored disc always covers the aperture the runner clipped.

PARALLELISM LIVES AT ONE LEVEL. `run(workers=None)` runs the blocks one after
the other, each block threaded inside (the runner's `Threader`).
`run(workers=W)` opens ONE process pool for the whole call and runs each block
SERIALLY inside its process. Never both: threads inside processes over-subscribe
the cores. A Windows process pool costs 2.5 to 4.4 s to spawn, and processes
beat threads by 1.15x to 1.7x in steady state
(validation/waveoptics_speed/fair_scaling_rerun.py), so a process pool pays only
when it stays WARM across many blocks.

PICKLING. `workers=W` sends the scenario, the geometry, the GridSpec and the
ScreenPlan to each process one time. Those are dataclasses, so they pickle. The
`cn2` callable is NOT sent: the parent plans the screens one time and the
workers get the finished plan, so a lambda `cn2` is safe here. The scenario and
the geometry must still be picklable objects at module level.

Sources:
- The seed contract, the trial body and the field store:
  olb.waveoptics.turbulence.run.
- The block fingerprint: olb.waveoptics.turbulence.fingerprint.cache_key.
- The Term reducer that reads a loaded result: olb.models.waveoptics.
"""

import json
import os
import time
from dataclasses import replace

import numpy as np

from ..grid import GridSpec
from ..threader import Threader
from .fingerprint import cache_key
from .run import (FieldPatch, TurbTrial, TurbWaveResult, _field_patch,
                  clip_terminal,
                  _resolve_seed, propagate_turbulent_scenario, recollect,
                  recouple)
from .sampling import ScreenPlan, turbulent_grid

# The manifest name and the block name. A block file holds one block only, so a
# stopped campaign keeps every finished block.
MANIFEST_NAME = "manifest.json"
PATCH_NAME = "patch_indices.npy"

# The columns of a block file. NaN marks a trial scalar that the runner left
# None (no SMF detector, no MMF detector, no uplink overlap).
_COLUMNS = ("collected_power", "smf_eta", "mmf_eta", "eta_turb", "wall_time_s")

# The worker state of a process of the pool. The initializer fills it one time,
# so the payload crosses the process boundary once, not once per block.
_W = {}


def _block_name(b):
    """Give the file name of block b."""
    return f"block_{int(b):05d}.npz"


def _nan_to_none(x):
    """Give None for a NaN, else the float."""
    return None if np.isnan(x) else float(x)


def _none_to_nan(x):
    """Give NaN for a None, else the float."""
    return np.nan if x is None else float(x)


def _columns_of(result):
    """Pack a TurbWaveResult into the plain arrays of a block file.

    The arrays pickle cheaply and they store directly. A TurbTrial list does
    not: it is a list of frozen dataclasses, and it costs much more to send
    between processes.

    Args:
        result: the TurbWaveResult of one block.

    Returns:
        A dict of numpy arrays, one for each column, plus "fields".
    """
    out = {c: np.array([_none_to_nan(getattr(t, c)) for t in result.trials],
                       dtype=np.float64) for c in _COLUMNS}
    out["fields"] = np.asarray(result.fields, dtype=np.complex64)
    return out


def _sizing_scenario(scenario, aperture_m):
    """Copy a scenario with a different CLIP aperture.

    The rule is the rule of run.clip_terminal: a SpaceScenario clips at
    `ground` in EVERY direction (the field is always the downlink slab at the
    ground, and an uplink reads it through reciprocity), and a
    TerrestrialScenario clips at its receive terminal (`far` on a forward
    link, `near` on a reverse link). The sizer reads the same terminal, so the
    copy moves the aperture the sizer sees.

    Args:
        scenario:   a SpaceScenario or a TerrestrialScenario.
        aperture_m: the clip aperture diameter of the copy, in m.

    Returns:
        A copy. The input scenario does not change.
    """
    rx = replace(clip_terminal(scenario), aperture_m=float(aperture_m))
    if hasattr(scenario, "ground"):
        role = "ground"
    else:
        role = "near" if scenario.direction == "reverse" else "far"
    return replace(scenario, **{role: rx})


def _init_worker(payload):
    """Fill the worker state of one process of the pool.

    The pool calls this ONE time for each process, so the scenario, the
    geometry, the grid and the plan cross the process boundary once.

    Args:
        payload: the dict that Campaign.run builds.
    """
    _W.clear()
    _W.update(payload)


def _run_block(b):
    """Run one block inside a worker process, and give back the columns.

    The block runs SERIALLY (threader=None). The parallelism lives at the
    process level only, so the cores are not over-subscribed.

    Args:
        b: the block index.

    Returns:
        The pair (b, the column dict).
    """
    res = propagate_turbulent_scenario(
        _W["scenario"], _W["geometry"], threader=None,
        start_index=int(b) * _W["kwargs"]["block_size"],
        grid=_W["grid"], plan=_W["plan"],
        n_trials=_W["kwargs"]["block_size"],
        seed=_W["kwargs"]["seed"], preset=_W["kwargs"]["preset"],
        patch_radius_m=_W["kwargs"]["patch_radius_m"],
        L0_m=_W["kwargs"]["L0_m"],
        subharmonics=_W["kwargs"]["subharmonics"],
        screen_generator=_W["kwargs"]["screen_generator"],
        precision=_W["kwargs"]["precision"])
    return int(b), _columns_of(res)


class Campaign:
    """A set of turbulent trials on disk, in blocks.

    Attributes:
        root_dir:       the directory of the block files and the manifest.
        scenario:       the scenario of the trials.
        geometry:       the geometry of the trials.
        seed:           the integer base seed.
        preset:         the quality preset name.
        block_size:     the number of trials in one block.
        patch_radius_m: the radius of the stored receive-field disc, in m.
        grid:           the GridSpec of every trial.
        plan:           the ScreenPlan of every trial.
        patch:          the FieldPatch of the stored columns.
        precision:      "double" or "single", the arithmetic of every trial.
    """

    def __init__(self, scenario, geometry, root_dir, *, seed,
                 preset="standard", block_size=100, patch_radius_m=None,
                 sizing_aperture_m=None, grid=None, plan=None, cn2=None,
                 hs=None, cn2_profile=None, h_top_m=None, L0_m=np.inf,
                 subharmonics=True, screen_generator="olb",
                 precision="single"):
        """Open a campaign, or make a new one.

        A missing `root_dir` is made. An EXISTING `root_dir` is checked: the
        fingerprint, the seed, the preset, the block size and the patch radius
        must match, and a mismatch raises. The grid and the plan then come from
        the manifest, NOT from a new sizing call. So a resumed campaign never
        re-sizes, and the atmosphere of a new block is the atmosphere of the
        old blocks.

        Args:
            scenario:      a SpaceScenario or a TerrestrialScenario.
            geometry:      the link geometry. It must give ONE range.
            root_dir:      the campaign directory.
            seed:          the integer base seed. It is REQUIRED, because a
                           campaign must repeat.
            preset:        the name of a preset in sampling.PRESETS.
            block_size:    the number of trials in one block.
            patch_radius_m: the radius of the stored field disc, in m. None
                           takes sizing_aperture_m / 2 when a sizing aperture is
                           given, else half the aperture of the clip terminal
                           (run.clip_terminal: the ground terminal of a space
                           scenario in every direction, the receive terminal
                           of a terrestrial one).
            sizing_aperture_m: an optional LARGER receive aperture that sizes
                           the grid. The trials still run with the original
                           scenario. Use it to store one field that serves every
                           smaller receive aperture.
            grid:          an optional GridSpec. The plan is then still planned
                           from the Cn2 inputs.
            plan:          an optional ScreenPlan. Give it WITH grid to hold
                           the grid fixed and move the screens only (a
                           convergence study). Both enter the fingerprint, so
                           a different plan is a different campaign.
            cn2:           an optional callable cn2(h) (space). See
                           turbulent_grid.
            hs, cn2_profile: the legacy discrete Cn2 profile (space).
            h_top_m:       the atmosphere top for the continuous integral.
            L0_m:          the outer scale of the screens, in m.
            subharmonics:  True adds the three subharmonic levels.
            screen_generator: "olb" (the default) or "aotools".
            precision:     "single" (the default) or "double". "single" runs
                           every trial in complex64, with float32 phase
                           screens. WHY: a campaign is memory-bandwidth bound,
                           so half the bytes for each element gives a real
                           speed-up. The manifest stores the value, and a
                           reopen with a different value raises. The value also
                           enters the fingerprint, so a single-precision
                           campaign is a separate store. CAUTION: a
                           single-precision campaign is a DIFFERENT record.
                           Validate it against a double-precision run of the
                           same seed before a budget reads it. See
                           validation/precision.

        Raises:
            ValueError: the seed is not an integer, the precision name is
                        unknown, or an existing campaign in this directory
                        holds different settings.
        """
        if precision not in ("double", "single"):
            raise ValueError(
                f"Campaign: precision must be 'double' or 'single', not "
                f"{precision!r}.")
        if seed is None or isinstance(seed, np.random.Generator):
            raise ValueError(
                "Campaign needs an integer seed. A campaign grows over more "
                "than one session, so its trials must repeat. Pass seed=<int>.")
        self.scenario = scenario
        self.geometry = geometry
        self.root_dir = str(root_dir)
        self.seed = int(seed)
        self.preset = preset if isinstance(preset, str) else preset.name
        self.block_size = int(block_size)
        self.screen_generator = screen_generator
        self.precision = precision
        self.L0_m = float(L0_m)
        self.subharmonics = bool(subharmonics)
        self.sizing_aperture_m = (None if sizing_aperture_m is None
                                  else float(sizing_aperture_m))

        if patch_radius_m is None:
            base = (self.sizing_aperture_m if self.sizing_aperture_m is not None
                    else clip_terminal(scenario).aperture_m)
            patch_radius_m = float(base) / 2.0
        self.patch_radius_m = float(patch_radius_m)

        self.fingerprint = cache_key(
            scenario, geometry, preset=self.preset, seed=self.seed,
            screen_generator=screen_generator, L0_m=L0_m,
            subharmonics=subharmonics, cn2=cn2, hs=hs,
            cn2_profile=cn2_profile, h_top_m=h_top_m,
            block_size=self.block_size, grid=grid, plan=plan,
            precision=self.precision)

        os.makedirs(self.root_dir, exist_ok=True)
        manifest_path = os.path.join(self.root_dir, MANIFEST_NAME)
        if os.path.exists(manifest_path):
            with open(manifest_path, "r", encoding="utf-8") as fh:
                man = json.load(fh)
            self._check_manifest(man)
            self.grid = GridSpec(size_m=man["grid"]["size_m"],
                                 n=int(man["grid"]["n"]),
                                 scaled=bool(man["grid"]["scaled"]))
            self.plan = ScreenPlan(
                z_m=np.array(man["plan"]["z_m"], dtype=float),
                cn2_int_m13=np.array(man["plan"]["cn2_int_m13"], dtype=float),
                r0_m=np.array(man["plan"]["r0_m"], dtype=float),
                sigma2_r=np.array(man["plan"]["sigma2_r"], dtype=float),
                z_total_m=float(man["plan"]["z_total_m"]),
                r0_total_m=float(man["plan"]["r0_total_m"]),
                direction=man["plan"]["direction"])
            self.patch = FieldPatch(
                radius_m=self.patch_radius_m, n=int(man["patch"]["n"]),
                pixel_m=float(man["patch"]["pixel_m"]),
                indices=np.load(os.path.join(self.root_dir, PATCH_NAME)))
        else:
            sizer_scenario = (scenario if self.sizing_aperture_m is None else
                              _sizing_scenario(scenario, self.sizing_aperture_m))
            sized_grid, sized_plan, _ = turbulent_grid(
                sizer_scenario, geometry, preset=self.preset, cn2=cn2, hs=hs,
                cn2_profile=cn2_profile, h_top_m=h_top_m, L0_m=L0_m)
            self.grid = sized_grid if grid is None else grid
            self.plan = sized_plan if plan is None else plan
            self.patch = _field_patch(self.grid, self.patch_radius_m)
            np.save(os.path.join(self.root_dir, PATCH_NAME), self.patch.indices)
            self._write_manifest(manifest_path)

    # ---- the manifest -----------------------------------------------------

    def _check_manifest(self, man):
        """Raise when a stored manifest does not match this campaign.

        Args:
            man: the stored manifest dict.

        Raises:
            ValueError: a field differs. The message names that field.
        """
        # The plain fields come FIRST. The seed and the preset also enter the
        # fingerprint, so a fingerprint-first order would name the hash and hide
        # the field that really differs.
        want = {"seed": self.seed, "preset": self.preset,
                "block_size": self.block_size,
                "patch_radius_m": self.patch_radius_m,
                "sizing_aperture_m": self.sizing_aperture_m,
                "precision": self.precision,
                "fingerprint": self.fingerprint}
        # A manifest that a version before the precision switch wrote holds no
        # "precision" key. It is a double-precision store, so read it as one.
        defaults = {"precision": "double"}
        for field, value in want.items():
            got = man.get(field, defaults.get(field))
            if got != value:
                raise ValueError(
                    f"the campaign in {self.root_dir} was made with "
                    f"{field}={got!r}, and this Campaign asks for "
                    f"{field}={value!r}. A stored campaign is ONE physics case. "
                    "Use a new directory, or match the stored settings.")

    def _write_manifest(self, path):
        """Write the manifest of a new campaign."""
        try:
            from ... import __version__ as olb_version
        except ImportError:
            olb_version = None
        man = {
            "fingerprint": self.fingerprint,
            "seed": self.seed,
            "preset": self.preset,
            "block_size": self.block_size,
            "patch_radius_m": self.patch_radius_m,
            "sizing_aperture_m": self.sizing_aperture_m,
            "screen_generator": self.screen_generator,
            "precision": self.precision,
            "L0_m": None if not np.isfinite(self.L0_m) else self.L0_m,
            "subharmonics": self.subharmonics,
            "olb_version": olb_version,
            "scenario": repr(self.scenario),
            "grid": {"size_m": float(self.grid.size_m), "n": int(self.grid.n),
                     "scaled": bool(self.grid.scaled)},
            "plan": {"z_m": self.plan.z_m.tolist(),
                     "cn2_int_m13": self.plan.cn2_int_m13.tolist(),
                     "r0_m": self.plan.r0_m.tolist(),
                     "sigma2_r": self.plan.sigma2_r.tolist(),
                     "z_total_m": float(self.plan.z_total_m),
                     "r0_total_m": float(self.plan.r0_total_m),
                     "direction": self.plan.direction},
            "patch": {"n": int(self.patch.n),
                      "pixel_m": float(self.patch.pixel_m),
                      "n_pixels": int(self.patch.indices.size)},
        }
        tmp = path + ".tmp"
        with open(tmp, "w", encoding="utf-8") as fh:
            json.dump(man, fh, indent=1)
        os.replace(tmp, path)

    # ---- the blocks -------------------------------------------------------

    def _block_path(self, b):
        """Give the path of block b."""
        return os.path.join(self.root_dir, _block_name(b))

    def _has_block(self, b):
        """Say if block b sits on disk."""
        return os.path.exists(self._block_path(b))

    def _write_block(self, b, cols):
        """Write one block file (an atomic replace)."""
        path = self._block_path(b)
        tmp = path + ".tmp.npz"
        np.savez(tmp, **cols)
        os.replace(tmp, path)

    def _read_block(self, b, fields=True):
        """Read one block file into a dict of arrays."""
        with np.load(self._block_path(b)) as z:
            cols = {c: z[c] for c in _COLUMNS}
            cols["fields"] = z["fields"] if fields else None
        return cols

    @property
    def n_stored(self):
        """The number of trials on disk, counting from block 0 with no gap."""
        b = 0
        while self._has_block(b):
            b += 1
        return b * self.block_size

    def _runner_kwargs(self):
        """Give the keyword arguments that every block run shares."""
        return {"block_size": self.block_size, "seed": self.seed,
                "preset": self.preset,
                "patch_radius_m": self.patch_radius_m, "L0_m": self.L0_m,
                "subharmonics": self.subharmonics,
                "screen_generator": self.screen_generator,
                "precision": self.precision}

    def run(self, n_trials, *, workers=None, progress=False):
        """Compute and store the MISSING blocks up to n_trials trials.

        A block that already sits on disk is not recomputed. The parent writes
        each block file as soon as that block arrives, so a killed campaign
        keeps every finished block.

        Args:
            n_trials: the number of trials the campaign must hold. The call
                      rounds it up to a whole number of blocks.
            workers:  None runs the blocks one after the other in this process,
                      each block threaded inside. An int W opens ONE process
                      pool of W processes for the whole call, and each block
                      runs serially inside its process.
            progress: True prints one line for each finished block.

        Returns:
            The number of trials on disk, an int.
        """
        n_blocks = -(-int(n_trials) // self.block_size)     # ceil division.
        missing = [b for b in range(n_blocks) if not self._has_block(b)]
        if not missing:
            return self.n_stored

        t0 = time.time()
        if workers is None:
            threader = Threader()
            for i, b in enumerate(missing):
                res = propagate_turbulent_scenario(
                    self.scenario, self.geometry, n_trials=self.block_size,
                    start_index=b * self.block_size, seed=self.seed,
                    preset=self.preset, grid=self.grid, plan=self.plan,
                    patch_radius_m=self.patch_radius_m, L0_m=self.L0_m,
                    subharmonics=self.subharmonics,
                    screen_generator=self.screen_generator,
                    precision=self.precision, threader=threader)
                self._write_block(b, _columns_of(res))
                if progress:
                    print(f"  block {b:5d} done "
                          f"({i + 1}/{len(missing)}, {time.time() - t0:.1f} s)")
        else:
            from concurrent.futures import ProcessPoolExecutor, as_completed
            payload = {"scenario": self.scenario, "geometry": self.geometry,
                       "grid": self.grid, "plan": self.plan,
                       "kwargs": self._runner_kwargs()}
            with ProcessPoolExecutor(max_workers=int(workers),
                                     initializer=_init_worker,
                                     initargs=(payload,)) as pool:
                futures = [pool.submit(_run_block, b) for b in missing]
                for i, future in enumerate(as_completed(futures)):
                    b, cols = future.result()
                    self._write_block(b, cols)
                    if progress:
                        print(f"  block {b:5d} done "
                              f"({i + 1}/{len(missing)}, "
                              f"{time.time() - t0:.1f} s)")
        return self.n_stored

    # ---- reading the store ------------------------------------------------

    def _blocks_for(self, n_trials):
        """Give the block indices and the trial count of a request."""
        stored = self.n_stored
        n = stored if n_trials is None else min(int(n_trials), stored)
        return range(-(-n // self.block_size)), n

    def load(self, n_trials=None, *, fields=True):
        """Assemble a TurbWaveResult from the stored blocks.

        The record is the record of a native run: the trial order is the trial
        order, and `seed_key` holds the TRUE trial index. So
        `olb.models.waveoptics.waveoptics_turbulence_term` reads it unchanged.

        Args:
            n_trials: the number of trials to load. None takes every stored
                      trial.
            fields:   True loads the stored field too. False leaves
                      `TurbWaveResult.fields` None, so a budget-only load stays
                      small.

        Returns:
            A TurbWaveResult.
        """
        blocks, n = self._blocks_for(n_trials)
        cols = {c: [] for c in _COLUMNS}
        stack = []
        for b in blocks:
            got = self._read_block(b, fields=fields)
            for c in _COLUMNS:
                cols[c].append(got[c])
            if fields:
                stack.append(got["fields"])
        packed = {c: np.concatenate(cols[c])[:n] for c in _COLUMNS}
        entropy = _resolve_seed(self.seed)
        trials = [
            TurbTrial(collected_power=_nan_to_none(packed["collected_power"][k]),
                      smf_eta=_nan_to_none(packed["smf_eta"][k]),
                      eta_turb=_nan_to_none(packed["eta_turb"][k]),
                      seed_key=(entropy, k),
                      wall_time_s=float(packed["wall_time_s"][k]),
                      mmf_eta=_nan_to_none(packed["mmf_eta"][k]),
                      detector_etas=None)
            for k in range(n)]
        return TurbWaveResult(
            trials=trials, grid=self.grid, plan=self.plan, report=None,
            preset=self.preset, seed_entropy=entropy,
            fields=(np.concatenate(stack)[:n] if fields else None),
            patch=self.patch if fields else None)

    def _stream(self, fn, n_trials=None):
        """Run fn on one block at a time, and join the results.

        The generator holds ONE block of fields in memory, so ten thousand
        trials never sit in RAM at the same time.

        Args:
            fn:       a callable fn(TurbWaveResult) -> a 1-D array.
            n_trials: the number of trials, or None for every stored trial.

        Returns:
            The concatenated array, cut to the request.
        """
        blocks, n = self._blocks_for(n_trials)
        out = []
        for b in blocks:
            got = self._read_block(b, fields=True)
            part = TurbWaveResult(trials=[], grid=self.grid, plan=self.plan,
                                  report=None, preset=self.preset,
                                  seed_entropy=self.seed,
                                  fields=got["fields"], patch=self.patch)
            out.append(fn(part))
        return np.concatenate(out)[:n]

    def recouple(self, detector, aperture_m=None, obscuration_ratio=None,
                 n_trials=None):
        """Couple the STORED fields into a detector, with no new propagation.

        Args:
            detector:          an SMF, an MMF, an Aperture, a Camera, or None.
            aperture_m:        the receive aperture diameter, in m. None takes
                               the aperture of the scenario receive terminal.
            obscuration_ratio: the central obscuration. None takes the value of
                               the scenario receive terminal.
            n_trials:          the number of trials. None takes every stored
                               trial.

        Returns:
            A float array of the coupling efficiency of each trial.
        """
        rx = self.scenario.rx_terminal
        a = rx.aperture_m if aperture_m is None else float(aperture_m)
        o = (rx.obscuration_ratio if obscuration_ratio is None
             else float(obscuration_ratio))
        lam = self.scenario.tx_terminal.wavelength_m
        return self._stream(
            lambda part: recouple(part, detector, a, o, lam), n_trials)

    def recollect(self, aperture_m=None, obscuration_ratio=None,
                  n_trials=None):
        """Give the collected power of each STORED trial, in grid units.

        The value is NOT normalised: it holds no vacuum reference. Take the
        RATIO of two trials, or divide by your own reference. See
        olb.waveoptics.turbulence.run.recollect.

        Args:
            aperture_m:        the receive aperture diameter, in m. None takes
                               the aperture of the scenario receive terminal.
            obscuration_ratio: the central obscuration. None takes the value of
                               the scenario receive terminal.
            n_trials:          the number of trials. None takes every stored
                               trial.

        Returns:
            A float array, one value for each trial.
        """
        rx = self.scenario.rx_terminal
        a = rx.aperture_m if aperture_m is None else float(aperture_m)
        o = (rx.obscuration_ratio if obscuration_ratio is None
             else float(obscuration_ratio))
        return self._stream(
            lambda part: recollect(part, a, o), n_trials)


if __name__ == '__main__':
    import shutil
    import tempfile
    import warnings

    from ...geometry import CircularOrbit
    from ...models.waveoptics import waveoptics_turbulence_term
    from ...scenario import Channel, SpaceScenario
    from ...terminal import SMF, Terminal, Transmitter

    t_start = time.time()
    lam = 1550e-9
    ground = Terminal(aperture_m=0.40, wavelength_m=lam, detector=SMF(),
                      transmitter=Transmitter(waist_m=0.06))
    scn = SpaceScenario(ground=ground,
                        space=Terminal(aperture_m=0.30, wavelength_m=lam),
                        direction="downlink", channel=Channel())
    orbit = CircularOrbit(altitude_m=600e3, elevation_deg=[30.0])

    root = tempfile.mkdtemp(prefix="olb_campaign_selfcheck_")
    root2 = tempfile.mkdtemp(prefix="olb_campaign_selfcheck2_")
    root3 = tempfile.mkdtemp(prefix="olb_campaign_selfcheck3_")
    root4 = tempfile.mkdtemp(prefix="olb_campaign_selfcheck4_")
    root5 = tempfile.mkdtemp(prefix="olb_campaign_selfcheck5_")
    common = dict(seed=2024, preset="rapid", block_size=4)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            # ---- 1. a grow computes ONLY the missing block ----
            camp = Campaign(scn, orbit, root, **common)
            t0 = time.time()
            assert camp.run(8) == 8, camp.n_stored
            serial_s = time.time() - t0
            before = [open(camp._block_path(b), "rb").read() for b in (0, 1)]
            assert camp.run(12) == 12, camp.n_stored
            after = [open(camp._block_path(b), "rb").read() for b in (0, 1)]
            assert before == after, "a grow must not touch a stored block"

            got = camp.load(12)
            native = propagate_turbulent_scenario(
                scn, orbit, n_trials=12, seed=2024, preset="rapid",
                grid=camp.grid, plan=camp.plan,
                patch_radius_m=camp.patch_radius_m)
            for a, b in zip(got.trials, native.trials):
                assert a.collected_power == b.collected_power, (a, b)
                assert a.smf_eta == b.smf_eta, (a, b)
                assert a.seed_key == b.seed_key, (a.seed_key, b.seed_key)
            assert got.fields.shape == (12, camp.patch.indices.size)

            # ---- 2. a process pool gives the SAME block files ----
            camp2 = Campaign(scn, orbit, root2, **common)
            t0 = time.time()
            assert camp2.run(12, workers=2) == 12
            pool_s = time.time() - t0
            for b in range(3):
                with np.load(camp._block_path(b)) as za, \
                        np.load(camp2._block_path(b)) as zb:
                    for c in ("collected_power", "smf_eta", "fields"):
                        assert np.array_equal(za[c], zb[c]), (b, c)

            # ---- 3. the post-hoc coupling matches the in-run scalars ----
            eta_back = camp.recouple(ground.detector)
            eta_run = np.array([t.smf_eta for t in got.trials])
            assert np.all(np.abs(eta_back / eta_run - 1.0) < 1e-5), \
                (eta_back[:3], eta_run[:3])
            pw = camp.recollect()
            run_pw = np.array([t.collected_power for t in got.trials])
            assert abs((pw[0] / pw[1]) / (run_pw[0] / run_pw[1]) - 1.0) < 1e-5

            # ---- 4. the Term reducer reads a loaded record unchanged ----
            small = camp.load(12, fields=False)
            assert small.fields is None and small.patch is None
            term = waveoptics_turbulence_term(small, quantity="collected_power")
            assert term.mean_db is not None

            # ---- 5. a different seed raises ----
            try:
                Campaign(scn, orbit, root, seed=7, preset="rapid", block_size=4)
                raise AssertionError("a changed seed must raise ValueError")
            except ValueError as exc:
                assert "seed" in str(exc), str(exc)

            # ---- 6. a sizing aperture sizes the grid, not the trials ----
            camp3 = Campaign(scn, orbit, root3, seed=2024, preset="rapid",
                             block_size=4,
                             sizing_aperture_m=2 * ground.aperture_m)
            assert camp3.patch_radius_m == ground.aperture_m, \
                camp3.patch_radius_m
            assert camp3.run(4) == 4
            eta3 = camp3.recouple(ground.detector)
            assert eta3.size == 4 and np.all(eta3 > 0.0), eta3
            big = (camp3.grid.n, camp3.grid.size_m) != (camp.grid.n,
                                                        camp.grid.size_m)

            # ---- 7. an injected plan is stored, fingerprinted, and reopened --
            # A convergence study holds the grid and moves the screens only.
            thin = ScreenPlan(
                z_m=camp.plan.z_m[::2], cn2_int_m13=camp.plan.cn2_int_m13[::2],
                r0_m=camp.plan.r0_m[::2], sigma2_r=camp.plan.sigma2_r[::2],
                z_total_m=camp.plan.z_total_m, r0_total_m=camp.plan.r0_total_m,
                direction=camp.plan.direction)
            camp4 = Campaign(scn, orbit, root4, grid=camp.grid, plan=thin,
                             **common)
            assert camp4.plan is thin and camp4.grid is camp.grid
            assert camp4.fingerprint != camp.fingerprint, "a plan must key"
            reopened = Campaign(scn, orbit, root4, grid=camp.grid, plan=thin,
                                **common)
            assert reopened.fingerprint == camp4.fingerprint
            assert reopened.plan.z_m.size == thin.z_m.size
            try:
                Campaign(scn, orbit, root4, grid=camp.grid, **common)
                raise AssertionError("a dropped plan must raise ValueError")
            except ValueError as exc:
                assert "fingerprint" in str(exc), str(exc)

            # ---- 8. the precision is stored, checked and fingerprinted ----
            # The default is "single" (owner decision 2026-09-05), and the
            # manifest records it. "double" keeps the OLD key, so a campaign
            # stored before that date still opens with precision="double".
            import json as _json
            with open(os.path.join(root, MANIFEST_NAME), encoding="utf-8") as fh:
                assert _json.load(fh)["precision"] == "single"
            camp5 = Campaign(scn, orbit, root5, precision="double", **common)
            assert camp5.fingerprint != camp.fingerprint, "precision must key"
            try:
                Campaign(scn, orbit, root5, **common)
                raise AssertionError("a changed precision must raise")
            except ValueError as exc:
                assert "precision" in str(exc), str(exc)
            try:
                Campaign(scn, orbit, root5, precision="half", **common)
                raise AssertionError("an unknown precision must raise")
            except ValueError as exc:
                assert "single" in str(exc), str(exc)

        print("campaign self-check, downlink 30 deg, rapid preset, "
              f"block_size {common['block_size']}:")
        print(f"  grid                    {camp.grid.n:11d} px, "
              f"{camp.grid.size_m:.3f} m")
        print(f"  screens                 {camp.plan.z_m.size:11d}")
        print(f"  patch pixels            {camp.patch.indices.size:11d}")
        print(f"  stored kB per trial     "
              f"{got.fields[0].nbytes / 1024:11.1f}")
        print(f"  8 trials, serial        {serial_s:11.1f} s")
        print(f"  12 trials, 2 processes  {pool_s:11.1f} s")
        print(f"  SMF eta, in run         {eta_run[0]:11.6f}")
        print(f"  SMF eta, recoupled      {eta_back[0]:11.6f}")
        print(f"  turbulence Term         {term.mean_db:11.3f} dB")
        print(f"  sizing grid differs     {str(big):>11s} "
              f"({camp3.grid.n} px, {camp3.grid.size_m:.3f} m)")
        print("")
        print(f"(elapsed {time.time() - t_start:.1f} s)")
        print("self-check passed")
    finally:
        for d in (root, root2, root3, root4, root5):
            shutil.rmtree(d, ignore_errors=True)
