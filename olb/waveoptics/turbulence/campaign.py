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
from dataclasses import dataclass, replace

import numpy as np

from ..field import field_dtype
from ..grid import GridSpec
from ..priority import boost_process_priority
from ..resources import auto_workers, worker_memory_bytes
from ..threader import Threader
from .fingerprint import cache_key
from .run import (FieldPatch, TurbTrial, TurbWaveResult, _check_aperture,
                  _crop_array, _field_patch, _patch_field, _PointAheadRunner,
                  _PostCorrector,
                  _PostTail, _resolve_compensation, _resolve_point_ahead,
                  _resolve_screen_margin, _screen_draw_n, _transmit_mode_crop,
                  _uplink_ground, clip_terminal,
                  _resolve_seed, propagate_turbulent_scenario,
                  space_vacuum_baseline)
from .sampling import PRESETS, ScreenPlan, resolve_outer_scale, turbulent_grid
from .splitstep import super_gaussian_boundary

# The manifest name and the block name. A block file holds one block only, so a
# stopped campaign keeps every finished block.
MANIFEST_NAME = "manifest.json"
PATCH_NAME = "patch_indices.npy"

# The margin factor of the DEFAULT stored patch radius. A NEW campaign with no
# patch_radius_m stores a disc 1.5x the aperture radius, so the field carries
# MARGIN beyond the aperture and a post-hoc MMF re-couple can upsample and clip
# a resolved aperture edge (see olb.waveoptics.mmf). A REOPENED campaign reads
# its stored radius from the manifest, so this factor never breaks an older
# store. See Campaign.__init__.
PATCH_MARGIN_FACTOR = 1.5

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

    The OPTIONAL "screen_phase" array joins the dict only when the run stored
    one. The three POINT-AHEAD arrays join it only when the run made a
    point-ahead pass. So the block file of a default campaign holds exactly the
    arrays it always held.

    Args:
        result: the TurbWaveResult of one block.

    Returns:
        A dict of numpy arrays, one for each column, plus "fields",
        "screen_phase" when the run stored it, and "eta_turb_pa", "fields_pa"
        and "screen_phase_pa" when the run made a point-ahead pass.
    """
    out = {c: np.array([_none_to_nan(getattr(t, c)) for t in result.trials],
                       dtype=np.float64) for c in _COLUMNS}
    out["fields"] = np.asarray(result.fields, dtype=np.complex64)
    if result.screen_phase is not None:
        out["screen_phase"] = np.asarray(result.screen_phase, dtype=np.float32)
    # THE POINT-AHEAD COLUMNS. "eta_turb_pa" is (n_trials, n_angles); a trial
    # with no overlap (a downlink) holds NaN in every angle.
    if result.point_ahead_rad is not None:
        n_angles = len(result.point_ahead_rad)
        out["eta_turb_pa"] = np.array(
            [[np.nan] * n_angles if t.eta_turb_pa is None
             else [_none_to_nan(v) for v in t.eta_turb_pa]
             for t in result.trials], dtype=np.float64).reshape(
                 len(result.trials), n_angles)
        if result.fields_pa is not None:
            out["fields_pa"] = np.asarray(result.fields_pa, dtype=np.complex64)
        if result.screen_phase_pa is not None:
            out["screen_phase_pa"] = np.asarray(result.screen_phase_pa,
                                                dtype=np.float32)
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

    When the payload asks for it, the worker also boosts ITS OWN priority:
    Windows does not pass the power-throttling opt-out of the parent to a
    spawned child (see olb.waveoptics.priority).

    Args:
        payload: the dict that Campaign.run builds.
    """
    if payload.get("boost", False):
        boost_process_priority()
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
        precision=_W["kwargs"]["precision"],
        fft_backend=_W["kwargs"]["fft_backend"],
        compensation=_W["kwargs"]["compensation"],
        store_screen_phase=_W["kwargs"]["store_screen_phase"],
        point_ahead_rad=_W["kwargs"]["point_ahead_rad"],
        screen_margin_m=_W["kwargs"]["screen_margin_m"],
        boost=False)        # the worker boosted itself in _init_worker.
    return int(b), _columns_of(res)


def _read_block_file(path, fields=True):
    """Read one block file into a dict of arrays.

    The "screen_phase" array loads whenever the file holds one, and the
    `fields` flag governs the field only. The two are separate: a budget load
    asks for no field, and a compensated read still needs the phase.

    Args:
        path:   the path of the block file.
        fields: True loads the stored field. False leaves it None.

    THE POINT-AHEAD ARRAYS read the same way: each one is None when the file
    holds none. So an OLD block file, written before the point ahead existed,
    still reads.

    Returns:
        A dict of arrays: one for each column of _COLUMNS, plus "fields",
        "screen_phase", "eta_turb_pa", "fields_pa" and "screen_phase_pa". Each
        optional array is None when the file holds none, or when the caller
        asks for no field.
    """
    with np.load(path) as z:
        cols = {c: z[c] for c in _COLUMNS}
        cols["fields"] = z["fields"] if fields else None
        for name in ("screen_phase", "eta_turb_pa", "screen_phase_pa"):
            cols[name] = z[name] if name in z.files else None
        cols["fields_pa"] = (z["fields_pa"]
                             if fields and "fields_pa" in z.files else None)
    return cols


@dataclass(frozen=True)
class TrialRecord:
    """One stored trial, as `Campaign.map_trials` hands it to a callable.

    THE PER-BLOCK CONTEXT. `context` is a plain dict that lives for ONE block.
    A callable that needs a cached object (an aperture mask, a fibre mode, a
    modal basis) builds it on the first trial of the block and it keeps it
    there. So the build happens one time for each block, not one time for each
    trial, and it never crosses a process boundary. See `_CoupleTrials` for
    the pattern.

    Attributes:
        row:          the trial index inside the campaign.
        array:        the rebuilt field of the trial, on the square CROP of
                      the patch. It is None when the call asks for no field.
        patch:        the FieldPatch of the campaign.
        lam:          the wavelength, in m.
        scalars:      the stored scalars of the trial, as a dict. The keys are
                      the names of _COLUMNS. A NaN marks a value the runner
                      left None.
        screen_phase: the stored summed screen phase of the trial, at the
                      patch pixels, or None.
        context:      the per-block dict (see above).
        arrays_pa:    the rebuilt POINT-AHEAD field of the trial, one array for
                      each angle, on the same crop as `array`. It is None when
                      the campaign made no point-ahead pass, or when the call
                      asks for no field.
        screen_phase_pa: the stored summed screen phase of each point-ahead
                      window, one row for each angle, or None.
    """

    row: int
    array: np.ndarray
    patch: FieldPatch
    lam: float
    scalars: dict
    screen_phase: np.ndarray
    context: dict
    arrays_pa: tuple = None
    screen_phase_pa: np.ndarray = None

    def field(self, compact=True):
        """Give the trial as a Field.

        Args:
            compact: True wraps the crop (the default). False pads the crop
                     back to the full grid, which a focal-plane quantity needs.

        Returns:
            A Field.
        """
        array = self.array
        if not compact:
            n = int(self.patch.n)
            crop = self.patch.crop()
            full = np.zeros((n, n), dtype=array.dtype)
            o, s = crop.offset, crop.side
            full[o:o + s, o:o + s] = array
            array = full
        return _patch_field(self.patch, array, self.lam)


def _apply_block(fn, cols, patch, lam, base_row, fields=True, compact=True):
    """Run one callable over the trials of one block.

    The block shares ONE context dict, so a per-block build happens one time.

    Args:
        fn:       the per-trial callable fn(TrialRecord) -> a value.
        cols:     the column dict of the block.
        patch:    the FieldPatch of the campaign.
        lam:      the wavelength, in m.
        base_row: the campaign row index of the first trial of the block.
        fields:   True rebuilds the field of each trial.
        compact:  True rebuilds on the crop (the default). False rebuilds on
                  the full grid.

    Returns:
        An array of the values, in trial order.
    """
    context = {}
    phase = cols.get("screen_phase")
    stack = cols.get("fields")
    stack_pa = cols.get("fields_pa")
    phase_pa = cols.get("screen_phase_pa")
    eta_pa = cols.get("eta_turb_pa")
    out = []
    for k in range(int(cols[_COLUMNS[0]].size)):
        array = (_crop_array(patch, stack[k], compact=compact)
                 if fields and stack is not None else None)
        arrays_pa = (tuple(_crop_array(patch, stack_pa[i, k], compact=compact)
                           for i in range(stack_pa.shape[0]))
                     if fields and stack_pa is not None else None)
        scalars = {c: cols[c][k] for c in _COLUMNS}
        if eta_pa is not None:
            scalars["eta_turb_pa"] = eta_pa[k]
        out.append(fn(TrialRecord(
            row=int(base_row) + k, array=array, patch=patch, lam=float(lam),
            scalars=scalars,
            screen_phase=None if phase is None else phase[k],
            context=context, arrays_pa=arrays_pa,
            screen_phase_pa=(None if phase_pa is None else phase_pa[:, k]))))
    return np.asarray(out)


def _map_block(b):
    """Run the callable of the worker state over one block.

    Args:
        b: the block index.

    Returns:
        The pair (b, the value array of the block).
    """
    m = _W["map"]
    cols = _read_block_file(os.path.join(m["root_dir"], _block_name(b)),
                            fields=m["fields"])
    if m["screen_phase"] is False:
        cols["screen_phase"] = None
        cols["screen_phase_pa"] = None
    return int(b), _apply_block(m["fn"], cols, m["patch"], m["lam"],
                                int(b) * int(m["block_size"]),
                                fields=m["fields"], compact=m["compact"])


class _CoupleTrials:
    """The per-trial callable of `Campaign.recouple`.

    It builds the `_PostTail` of the block one time, in the record context.
    """

    def __init__(self, detector, aperture_m, obscuration_ratio, lam,
                 compact=True):
        self.detector = detector
        self.aperture_m = float(aperture_m)
        self.obscuration_ratio = float(obscuration_ratio)
        self.lam = float(lam)
        self.compact = bool(compact)

    def _tail(self, rec):
        """Give the cached tail of the block."""
        tail = rec.context.get("tail")
        if tail is None:
            tail = _PostTail(rec.patch, self.aperture_m,
                             self.obscuration_ratio, self.lam,
                             detector=self.detector, compact=self.compact)
            rec.context["tail"] = tail
        return tail

    def __call__(self, rec):
        eta = self._tail(rec).eta(rec.array)
        return np.nan if eta is None else float(eta)


class _CollectTrials(_CoupleTrials):
    """The per-trial callable of `Campaign.recollect`.

    The wavelength does not enter a power, so any value serves.
    """

    def __init__(self, aperture_m, obscuration_ratio, compact=True):
        super().__init__(None, aperture_m, obscuration_ratio, 1.0,
                         compact=compact)

    def __call__(self, rec):
        return self._tail(rec).power(rec.array)


class _CompensateTrials(_CoupleTrials):
    """The per-trial callable of `Campaign.recouple_compensated`.

    It builds the modal basis of the block one time too. That build is the
    large cost of a compensated read.
    """

    def __init__(self, compensation, detector, aperture_m, obscuration_ratio,
                 lam, source, compact=True):
        super().__init__(detector, aperture_m, obscuration_ratio, lam,
                         compact=compact)
        self.compensation = tuple(compensation)
        self.source = source

    def __call__(self, rec):
        corrector = rec.context.get("corrector")
        if corrector is None:
            corrector = _PostCorrector(rec.patch, self.compensation,
                                       self.aperture_m,
                                       self.obscuration_ratio, self.source,
                                       compact=self.compact)
            rec.context["corrector"] = corrector
        if self.source == "screens" and rec.screen_phase is None:
            raise ValueError(
                "recouple_compensated: source='screens' needs the stored "
                "summed screen phase, and this campaign holds none. Run the "
                "campaign with store_screen_phase=True, or use "
                "source='slopes'.")
        eta = self._tail(rec).eta(
            corrector.correct(rec.array, rec.screen_phase))
        return np.nan if eta is None else float(eta)


class _PointAheadTrials:
    """The per-trial callable of `Campaign.recouple_point_ahead`.

    It reads the STORED point-ahead planes of a trial, it corrects each one
    with the BEACON estimate, and it gives the uplink reciprocity overlap of
    each angle. So it makes no propagation. See
    olb.waveoptics.turbulence.run.point_ahead_overlap.

    The modal basis is the large cost of the read, so the callable builds it
    ONE time for each block, in the record context. The transmit mode and the
    vacuum baseline come from the parent, so no block computes them again.
    """

    def __init__(self, compensation, aperture_m, obscuration_ratio, psi,
                 o_vac, source, compact=True):
        """Hold the fixed parts of one read.

        Args:
            compensation:      the RESOLVED stack, or None for no correction.
            aperture_m:        the clip aperture diameter, in m.
            obscuration_ratio: the central obscuration of that aperture.
            psi:               the ground transmit mode, on the same pixels as
                               the trial arrays.
            o_vac:             the free-space overlap baseline.
            source:            "screens", "slopes", or "gtilt".
            compact:           True works on the crop.
        """
        self.compensation = compensation
        self.aperture_m = float(aperture_m)
        self.obscuration_ratio = float(obscuration_ratio)
        self.conj_psi = np.conj(psi)
        self.o_vac = float(o_vac)
        self.source = source
        self.compact = bool(compact)

    def __call__(self, rec):
        if rec.arrays_pa is None:
            raise ValueError(
                "Campaign.recouple_point_ahead: this campaign stores no "
                "point-ahead field. Make it with point_ahead_rad AND "
                "patch_radius_m, and run it again.")
        coeffs = None
        if self.compensation is not None:
            corrector = rec.context.get("pa_corrector")
            if corrector is None:
                corrector = _PostCorrector(rec.patch, self.compensation,
                                           self.aperture_m,
                                           self.obscuration_ratio, self.source,
                                           compact=self.compact)
                rec.context["pa_corrector"] = corrector
            if self.source == "screens" and rec.screen_phase is None:
                raise ValueError(
                    "Campaign.recouple_point_ahead: source='screens' needs "
                    "the stored summed screen phase, and this campaign holds "
                    "none. Make it with store_screen_phase=True, or use "
                    "source='slopes'.")
            coeffs = corrector.coefficients(rec.array, rec.screen_phase)
        out = []
        for E in rec.arrays_pa:
            if coeffs is not None:
                E = corrector.modes.apply(E, coeffs, sign=-1)
            out.append(float(np.abs((E * self.conj_psi).sum()) ** 2)
                       / self.o_vac)
        return np.asarray(out, dtype=float)


class _RegeneratePointAhead:
    """The per-trial callable of `Campaign.point_ahead`.

    It REBUILDS the screens of a trial from the seeds of the campaign and it
    propagates them at each asked angle. So it answers an angle the campaign
    never stored. The screen factory, the boundary mask, the modal basis and
    the vacuum baseline build ONE time for each block, in the record context.
    See olb.waveoptics.turbulence.run.point_ahead_regenerate.
    """

    def __init__(self, record, angles, compensation, scenario, geometry,
                 **options):
        """Hold the fixed parts of one read.

        Args:
            record:   a small TurbWaveResult that names the grid, the plan,
                      the preset, the seed entropy and the screen side.
            angles:   the point-ahead angles, in rad.
            compensation: None, "terminal", or a list of stages.
            scenario: the SpaceScenario of the campaign.
            geometry: the link geometry.
            options:  the runner options of _PointAheadRunner.
        """
        self.record = record
        self.angles = angles
        self.compensation = compensation
        self.scenario = scenario
        self.geometry = geometry
        self.options = dict(options)

    def __call__(self, rec):
        runner = rec.context.get("pa_runner")
        if runner is None:
            runner = _PointAheadRunner(self.record, self.angles,
                                       self.compensation, self.scenario,
                                       self.geometry, **self.options)
            rec.context["pa_runner"] = runner
        # THE ROW IS THE TRIAL INDEX. A campaign runs block b from the
        # start_index b * block_size, so the campaign row IS the seed index.
        return runner.trial(rec.row, stored_phase=rec.screen_phase,
                            patch=rec.patch)


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
        compensation:   the RESOLVED perfect-AO stack of every trial, or None.
        n_modes_corrected: the number of removed Noll modes, or 0.
        store_screen_phase: True keeps the summed screen phase of each trial.
        point_ahead_rad: the RESOLVED point-ahead angles of every trial, a
                        tuple of floats, or None.
        screen_margin_m: the extra screen width of the point-ahead windows, in
                        m. It is 0.0 for a campaign with no point-ahead pass.
        screen_n:       the pixel count of one drawn screen side. It equals
                        grid.n for a campaign with no point-ahead pass.
    """

    def __init__(self, scenario, geometry, root_dir, *, seed,
                 preset="standard", block_size=100, patch_radius_m=None,
                 sizing_aperture_m=None, grid=None, plan=None, cn2=None,
                 hs=None, cn2_profile=None, h_top_m=None, L0_m=None,
                 subharmonics=True, screen_generator="olb",
                 precision="single", fft_backend="numpy", compensation=None,
                 store_screen_phase=False, point_ahead_rad=None,
                 screen_margin_m=None):
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
            patch_radius_m: the radius of the stored field disc, in m. None on a
                           NEW campaign takes PATCH_MARGIN_FACTOR (1.5) times half
                           the sizing aperture when a sizing aperture is given,
                           else 1.5 times half the aperture of the clip terminal
                           (run.clip_terminal: the ground terminal of a space
                           scenario in every direction, the receive terminal of a
                           terrestrial one). So the default disc carries MARGIN
                           beyond the aperture, and a post-hoc MMF re-couple
                           (recouple with an MMF detector) can upsample and clip a
                           resolved aperture edge (see olb.waveoptics.mmf and
                           run._PostTail.eta). None on a REOPENED campaign reads
                           the stored radius from the manifest, NOT the default,
                           so a change of the default never breaks an older store.
                           Pass an explicit value to override the default.
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
            screen_generator: "olb" (the default), "olb-lean" (an OPT-IN,
                           not bit-identical) or "aotools".
            fft_backend:   "numpy" (the default, the backend of record),
                           "scipy" (an OPT-IN, 2026-09-06, faster, agreement
                           at the rounding level, NOT bit-identical) or
                           "cupy" (an OPT-IN, 2026-09-07, the CUDA device).
                           A "cupy" campaign runs its blocks in ONE process,
                           in the calling process, whatever `run(workers=)`
                           says: one device runs one stream. The name enters
                           the fingerprint when it is not "numpy", so every
                           stored key stays valid AND a GPU campaign never
                           mixes with a CPU one. See
                           olb.waveoptics.propagators.set_fft_backend.
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
            compensation:  None (the default, NO correction), the string
                           "terminal" (the compensation stack of the clip
                           terminal), or a list of TipTilt and AO stages. Each
                           trial then removes the first N Noll modes of the
                           wavefront over the receive aperture (perfect AO; see
                           olb.waveoptics.compensation). The RESOLVED stack
                           enters the fingerprint and the manifest, so a
                           corrected campaign never mixes with an uncorrected
                           one. The default keeps every stored key valid.
            store_screen_phase: True stores the summed screen phase of each
                           trial at the patch pixels, as float32. It is the
                           sensing source of the post-hoc SPACE correction
                           (`recouple_compensated`). It adds one array to each
                           block file, and it enters the fingerprint. The
                           default False stores nothing.
            point_ahead_rad: None (the default, NO point-ahead pass), the string
                           "geometry", a float, or a sequence of angles in rad.
                           Each angle adds one more propagation of the SAME
                           atmosphere through a laterally shifted window of each
                           screen, and each block file then holds the uplink
                           overlap, the field and the screen phase of every
                           angle. The RESOLVED angles enter the fingerprint and
                           the manifest. It needs a SPACE scenario. See
                           olb.waveoptics.turbulence.run.
            screen_margin_m: the extra screen width of the shifted windows, in
                           m. None on a NEW campaign reads the geometry; None on
                           a REOPENED campaign reads the stored value from the
                           manifest, exactly like patch_radius_m. The RESOLVED
                           value enters the fingerprint.

        Raises:
            ValueError: the seed is not an integer, the precision name is
                        unknown, or an existing campaign in this directory
                        holds different settings.
        """
        if fft_backend not in ("numpy", "scipy", "cupy"):
            raise ValueError(
                f"Campaign: fft_backend must be 'numpy', 'scipy' or 'cupy', "
                f"not {fft_backend!r}.")
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
        self.fft_backend = fft_backend
        # L0_m=None reads the site outer scale (25 m). Resolve it BEFORE the
        # fingerprint, so the key names the physical outer scale, not None, and
        # a campaign that asks for the site and one that gives 25 m share a key.
        self.L0_m = resolve_outer_scale(L0_m, scenario)
        self.subharmonics = bool(subharmonics)
        self.sizing_aperture_m = (None if sizing_aperture_m is None
                                  else float(sizing_aperture_m))
        # RESOLVE THE STACK BEFORE THE FINGERPRINT. The key then names the
        # stages, not the string "terminal", so a campaign that asks for the
        # terminal stack and a campaign that gives the same stages share a key.
        self.compensation, self.n_modes_corrected = _resolve_compensation(
            scenario, compensation)
        self.store_screen_phase = bool(store_screen_phase)

        # Load the manifest FIRST when this store exists, so a reopened campaign
        # reads its stored patch radius from the manifest, NOT from the None
        # default. A change of the default then never breaks an older store.
        manifest_path = os.path.join(self.root_dir, MANIFEST_NAME)
        stored_manifest = None
        if os.path.exists(manifest_path):
            with open(manifest_path, "r", encoding="utf-8") as fh:
                stored_manifest = json.load(fh)

        if patch_radius_m is None:
            if stored_manifest is not None and "patch_radius_m" in stored_manifest:
                # A REOPEN: the manifest is the source of truth for the stored
                # patch geometry.
                patch_radius_m = float(stored_manifest["patch_radius_m"])
            else:
                # A NEW campaign stores a disc with MARGIN beyond the aperture
                # (PATCH_MARGIN_FACTOR), so a post-hoc MMF re-couple can upsample
                # and clip a resolved aperture edge (olb.waveoptics.mmf).
                base = (self.sizing_aperture_m if self.sizing_aperture_m is not None
                        else clip_terminal(scenario).aperture_m)
                patch_radius_m = float(base) / 2.0 * PATCH_MARGIN_FACTOR
        self.patch_radius_m = float(patch_radius_m)
        # RESOLVE THE POINT-AHEAD ANGLES BEFORE THE FINGERPRINT, the same rule
        # as the compensation stack: the key names the ANGLES, not the string
        # "geometry". The MARGIN needs the finished screen plan, so it resolves
        # below, and the key follows it.
        self.point_ahead_rad = _resolve_point_ahead(point_ahead_rad, geometry)

        os.makedirs(self.root_dir, exist_ok=True)

        def fingerprint_of():
            """Give the content key of this campaign. It reads the margin."""
            return cache_key(
                scenario, geometry, preset=self.preset, seed=self.seed,
                screen_generator=screen_generator, L0_m=self.L0_m,
                subharmonics=subharmonics, cn2=cn2, hs=hs,
                cn2_profile=cn2_profile, h_top_m=h_top_m,
                block_size=self.block_size, grid=grid, plan=plan,
                precision=self.precision, fft_backend=self.fft_backend,
                compensation=self.compensation,
                store_screen_phase=self.store_screen_phase,
                point_ahead_rad=self.point_ahead_rad,
                screen_margin_m=self.screen_margin_m)

        if stored_manifest is not None:
            man = stored_manifest
            # THE GRID AND THE PLAN COME FIRST. The screen margin reads the
            # plan, and the fingerprint reads the margin, so the check below
            # runs after all three.
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
            self._resolve_margin(screen_margin_m, man)
            self.fingerprint = fingerprint_of()
            self._check_manifest(man)
            self.patch = FieldPatch(
                radius_m=self.patch_radius_m, n=int(man["patch"]["n"]),
                pixel_m=float(man["patch"]["pixel_m"]),
                indices=np.load(os.path.join(self.root_dir, PATCH_NAME)))
        else:
            sizer_scenario = (scenario if self.sizing_aperture_m is None else
                              _sizing_scenario(scenario, self.sizing_aperture_m))
            sized_grid, sized_plan, _ = turbulent_grid(
                sizer_scenario, geometry, preset=self.preset, cn2=cn2, hs=hs,
                cn2_profile=cn2_profile, h_top_m=h_top_m, L0_m=self.L0_m)
            self.grid = sized_grid if grid is None else grid
            self.plan = sized_plan if plan is None else plan
            self._resolve_margin(screen_margin_m, None)
            self.fingerprint = fingerprint_of()
            self.patch = _field_patch(self.grid, self.patch_radius_m)
            np.save(os.path.join(self.root_dir, PATCH_NAME), self.patch.indices)
            self._write_manifest(manifest_path)

    # ---- the point-ahead screen margin ------------------------------------

    def _resolve_margin(self, screen_margin_m, man):
        """Set screen_margin_m and screen_n of this campaign.

        An EXPLICIT value wins. A REOPENED campaign then reads the stored value
        from the manifest, exactly as patch_radius_m does, so a change of the
        automatic rule never breaks an older store. A NEW campaign reads the
        geometry through run._resolve_screen_margin.

        Args:
            screen_margin_m: the caller value, or None.
            man:             the stored manifest dict, or None for a new store.
        """
        if screen_margin_m is None and self.point_ahead_rad is not None \
                and man is not None and "screen_margin_m" in man:
            screen_margin_m = float(man["screen_margin_m"])
        self.screen_margin_m = _resolve_screen_margin(
            screen_margin_m, self.point_ahead_rad, self.plan)
        self.screen_n = _screen_draw_n(self.grid.n, self.screen_margin_m,
                                       self.grid.pixel_m)

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
                "fft_backend": self.fft_backend,
                "compensation": repr(self.compensation),
                "store_screen_phase": self.store_screen_phase,
                "point_ahead_rad": (None if self.point_ahead_rad is None
                                    else list(self.point_ahead_rad)),
                "screen_margin_m": self.screen_margin_m,
                "screen_n": self.screen_n,
                "fingerprint": self.fingerprint}
        # A manifest that a version before the precision switch wrote holds no
        # "precision" key. It is a double-precision store, so read it as one.
        # The same for the FFT backend: an older manifest is a numpy store. A
        # manifest from before the compensation switch is an UNCORRECTED store
        # with no stored screen phase. A manifest from before the point ahead
        # is a store with no point-ahead pass: no angle, no margin, and a screen
        # that is the grid.
        defaults = {"precision": "double", "fft_backend": "numpy",
                    "compensation": repr(None), "store_screen_phase": False,
                    "point_ahead_rad": None, "screen_margin_m": 0.0,
                    "screen_n": int(self.grid.n)}
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
            "fft_backend": self.fft_backend,
            "compensation": repr(self.compensation),
            "n_modes_corrected": int(self.n_modes_corrected),
            "store_screen_phase": self.store_screen_phase,
            "point_ahead_rad": (None if self.point_ahead_rad is None
                                else list(self.point_ahead_rad)),
            "screen_margin_m": float(self.screen_margin_m),
            "screen_n": int(self.screen_n),
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
        """Read one block file into a dict of arrays.

        A block that a run wrote with no screen-phase store holds no
        "screen_phase" array, so the value is then None.
        """
        return _read_block_file(self._block_path(b), fields=fields)

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
                "precision": self.precision,
                "fft_backend": self.fft_backend,
                "compensation": self.compensation,
                "store_screen_phase": self.store_screen_phase,
                "point_ahead_rad": self.point_ahead_rad,
                "screen_margin_m": self.screen_margin_m}

    def worker_memory_bytes(self):
        """Estimate the peak memory of one pool worker of this campaign.

        It reads the grid, the precision, the block size, the patch and the
        hop count of the plan. A POINT-AHEAD campaign also counts the screen
        noise that one trial keeps. See olb.waveoptics.resources.

        Returns:
            An int, in bytes.
        """
        patch_pixels = 0 if self.patch is None else int(self.patch.indices.size)
        n_screens = (0 if self.point_ahead_rad is None
                     else int(self.plan.z_m.size))
        return worker_memory_bytes(self.grid.n, self.precision,
                                   block_size=self.block_size,
                                   patch_pixels=patch_pixels,
                                   n_hops=int(self.plan.z_m.size) + 1,
                                   screen_n=self.screen_n,
                                   n_screens=n_screens)

    def auto_workers(self, *, cpu_fraction=0.9, memory_fraction=0.9):
        """Give the pool size that fills the machine for this campaign.

        The count is the smaller of the CPU limit (`cpu_fraction` of the
        logical cores) and the memory limit (`memory_fraction` of the free
        memory over `worker_memory_bytes()`), and it is at least 1. It is
        also capped at the block count, because more workers than blocks
        sit idle.

        Returns:
            The pair (workers, reason), with the reason "cpu", "memory" or
            "blocks".
        """
        k, why = auto_workers(self.worker_memory_bytes(),
                              cpu_fraction=cpu_fraction,
                              memory_fraction=memory_fraction)
        return k, why

    def run(self, n_trials, *, workers=None, progress=False, boost=True,
            cpu_fraction=0.9, memory_fraction=0.9):
        """Compute and store the MISSING blocks up to n_trials trials.

        A block that already sits on disk is not recomputed. The parent writes
        each block file as soon as that block arrives, so a killed campaign
        keeps every finished block.

        Args:
            n_trials: the number of trials the campaign must hold. The call
                      rounds it up to a whole number of blocks.
            workers:  None runs the blocks one after the other in this process,
                      each block threaded inside. With the "cupy" FFT backend
                      EVERY value acts as None (the call prints one line to
                      say so), because one device runs one stream. An int W
                      opens ONE process
                      pool of W processes for the whole call, and each block
                      runs serially inside its process. The string "auto"
                      opens a pool sized by `auto_workers()`: `cpu_fraction`
                      of the logical cores, held under `memory_fraction` of
                      the free memory, and never more than the missing
                      blocks.
            progress: True prints one line for each finished block, and the
                      pool size and its reason for "auto".
            cpu_fraction, memory_fraction: the two limits of "auto". They
                      are ignored for None and for an int.
            boost:    True (the default) raises this process to the Above
                      Normal priority class and opts it out of power
                      throttling (EcoQoS), and every pool worker does the
                      same for itself. A windowless run (ssh, WMI) is
                      throttled without it. A no-op off Windows. See
                      olb.waveoptics.priority.

        Returns:
            The number of trials on disk, an int.
        """
        n_blocks = -(-int(n_trials) // self.block_size)     # ceil division.
        missing = [b for b in range(n_blocks) if not self._has_block(b)]
        if not missing:
            return self.n_stored

        if self.fft_backend == "cupy" and workers is not None:
            # ONE DEVICE, ONE STREAM. A process pool would put several
            # processes on the same device, and each one would hold its own
            # CUDA context and its own copy of the field. So a GPU campaign
            # runs its blocks one after the other, here.
            print("  cupy backend: the blocks run one after the other in "
                  "this process (one device, one stream). The workers "
                  f"request ({workers!r}) is not used.")
            workers = None

        if isinstance(workers, str):
            if workers != "auto":
                raise ValueError(f"Campaign.run: workers must be None, an int "
                                 f"or 'auto', not {workers!r}.")
            workers, why = self.auto_workers(cpu_fraction=cpu_fraction,
                                             memory_fraction=memory_fraction)
            if workers > len(missing):
                workers, why = len(missing), "blocks"
            if progress:
                per = self.worker_memory_bytes() / 2 ** 20
                print(f"  auto pool: {workers} workers ({why} limit, "
                      f"{per:.0f} MiB per worker)")
        if boost:
            boost_process_priority()    # Threads inherit; workers boost themselves.
        t0 = time.time()
        if workers is None:
            # The CUDA route refuses a threader: one device, one stream.
            threader = None if self.fft_backend == "cupy" else Threader()
            for i, b in enumerate(missing):
                res = propagate_turbulent_scenario(
                    self.scenario, self.geometry, n_trials=self.block_size,
                    start_index=b * self.block_size, seed=self.seed,
                    preset=self.preset, grid=self.grid, plan=self.plan,
                    patch_radius_m=self.patch_radius_m, L0_m=self.L0_m,
                    subharmonics=self.subharmonics,
                    screen_generator=self.screen_generator,
                    precision=self.precision, threader=threader,
                    fft_backend=self.fft_backend,
                    compensation=self.compensation,
                    store_screen_phase=self.store_screen_phase,
                    point_ahead_rad=self.point_ahead_rad,
                    screen_margin_m=self.screen_margin_m,
                    boost=boost)
                self._write_block(b, _columns_of(res))
                if progress:
                    print(f"  block {b:5d} done "
                          f"({i + 1}/{len(missing)}, {time.time() - t0:.1f} s)")
        else:
            from concurrent.futures import ProcessPoolExecutor, as_completed
            payload = {"scenario": self.scenario, "geometry": self.geometry,
                       "grid": self.grid, "plan": self.plan,
                       "kwargs": self._runner_kwargs(),
                       "boost": bool(boost)}
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

        THE STORED SCREEN PHASE IS A SEPARATE ARRAY. `fields=False` drops the
        field and it KEEPS the phase, because a caller that reads the phase
        does not always want the field. The POINT-AHEAD arrays follow the same
        rule: `fields=False` drops `fields_pa` and it keeps `eta_turb_pa` and
        `screen_phase_pa`.

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
        phase = []
        # The POINT-AHEAD arrays. The per-angle arrays carry the ANGLE on the
        # first axis and the trial on the second, so they join on axis 1.
        eta_pa = []
        stack_pa = []
        phase_pa = []
        for b in blocks:
            got = self._read_block(b, fields=fields)
            for c in _COLUMNS:
                cols[c].append(got[c])
            if fields:
                stack.append(got["fields"])
            if got["screen_phase"] is not None:
                phase.append(got["screen_phase"])
            if got["eta_turb_pa"] is not None:
                eta_pa.append(got["eta_turb_pa"])
            if got["fields_pa"] is not None:
                stack_pa.append(got["fields_pa"])
            if got["screen_phase_pa"] is not None:
                phase_pa.append(got["screen_phase_pa"])
        packed = {c: np.concatenate(cols[c])[:n] for c in _COLUMNS}
        packed_eta_pa = (np.concatenate(eta_pa)[:n] if eta_pa else None)
        entropy = _resolve_seed(self.seed)
        trials = [
            TurbTrial(collected_power=_nan_to_none(packed["collected_power"][k]),
                      smf_eta=_nan_to_none(packed["smf_eta"][k]),
                      eta_turb=_nan_to_none(packed["eta_turb"][k]),
                      seed_key=(entropy, k),
                      wall_time_s=float(packed["wall_time_s"][k]),
                      mmf_eta=_nan_to_none(packed["mmf_eta"][k]),
                      detector_etas=None,
                      eta_turb_pa=(None if packed_eta_pa is None else
                                   tuple(_nan_to_none(v)
                                         for v in packed_eta_pa[k])))
            for k in range(n)]
        return TurbWaveResult(
            trials=trials, grid=self.grid, plan=self.plan, report=None,
            preset=self.preset, seed_entropy=entropy,
            fields=(np.concatenate(stack)[:n] if fields else None),
            patch=(self.patch if fields or phase or phase_pa else None),
            compensation=self.compensation,
            n_modes_corrected=self.n_modes_corrected,
            screen_phase=(np.concatenate(phase)[:n] if phase else None),
            point_ahead_rad=self.point_ahead_rad,
            screen_margin_m=self.screen_margin_m, screen_n=self.screen_n,
            fields_pa=(np.concatenate(stack_pa, axis=1)[:, :n]
                       if stack_pa else None),
            screen_phase_pa=(np.concatenate(phase_pa, axis=1)[:, :n]
                             if phase_pa else None))

    def field(self, row, *, compact=True):
        """Give one STORED trial back as a Field.

        This is the public reader of a stored receive field. It reads the block
        that holds the trial, it rebuilds the pixels of the patch, and it wraps
        them as a Field. So a diagnostic (a camera image, a phase map, a plot)
        reads one trial with no new propagation and with no whole-campaign
        load.

        Args:
            row:     the trial index inside the campaign.
            compact: True gives the square CROP that holds the stored patch
                     (the default). False gives the FULL grid, which a
                     focal-plane quantity needs (an MMF, a Camera).

        Returns:
            A Field.

        Raises:
            ValueError: the campaign stores no field, or the row is not on
                        disk.
        """
        if self.patch is None:
            raise ValueError(
                "this campaign stores no field. Make it with a "
                "patch_radius_m, and run it again.")
        row = int(row)
        if not 0 <= row < self.n_stored:
            raise ValueError(
                f"Campaign.field: the row {row} is not on disk. The campaign "
                f"holds {self.n_stored} trials.")
        b, k = divmod(row, self.block_size)
        cols = self._read_block(b, fields=True)
        array = _crop_array(self.patch, cols["fields"][k], compact=compact)
        lam = self.scenario.tx_terminal.wavelength_m
        return _patch_field(self.patch, array, lam)

    def map_trials(self, fn, *, n_trials=None, workers=None, fields=True,
                   screen_phase=None, compact=True, boost=True):
        """Run a callable over the stored trials, block by block.

        THIS IS THE ONE post-hoc primitive. `recouple`, `recollect` and
        `recouple_compensated` are thin wrappers over it, and a study writes
        its own `fn` for anything else.

        The method reads ONE block at a time, so ten thousand trials never sit
        in memory together. Each block shares one `context` dict, so a callable
        builds its cached objects (a mask, a fibre mode, a modal basis) ONE
        time for each block. See TrialRecord.

        `workers=W` opens the SAME kind of warm process pool that `run()` uses:
        one task for each block, and the trials of a block run serially inside
        the worker. `fn` and the block results must pickle.

        Args:
            fn:           a callable fn(TrialRecord) -> a scalar, a tuple or a
                          small array. Every trial must give the same shape.
            n_trials:     the number of trials. None takes every stored trial.
            workers:      None runs the blocks in this process (the default).
                          An int W opens a pool of W processes. The string
                          "auto" sizes the pool with `auto_workers()`.
            fields:       True rebuilds the field of each trial (the default).
                          False leaves TrialRecord.array None, which is
                          cheaper for a scalar-only pass.
            screen_phase: None gives the stored screen phase when the campaign
                          holds one (the default). True raises when it holds
                          none. False never reads it.
            compact:      True hands `fn` the square CROP of the patch (the
                          default). False hands it the FULL grid, which costs
                          much more.
            boost:        True (the default) raises each pool worker to the
                          Above Normal priority, as `run()` does. It has no
                          effect when the call runs in this process.

        Returns:
            An array of the values, in trial order. The first axis is the
            trial.

        Raises:
            ValueError: the campaign stores no field and the call asks for
                        one, screen_phase=True and the campaign holds none, or
                        the workers value is not None, an int or "auto".
        """
        if fields and self.patch is None:
            raise ValueError(
                "this campaign stores no field. Make it with a "
                "patch_radius_m, and run it again, or pass fields=False.")
        if screen_phase is True and not self.store_screen_phase:
            raise ValueError(
                "Campaign.map_trials: this campaign stores no screen phase. "
                "Make it with store_screen_phase=True, and run it again.")
        blocks, n = self._blocks_for(n_trials)
        blocks = list(blocks)
        if isinstance(workers, str):
            if workers != "auto":
                raise ValueError(
                    "Campaign.map_trials: workers must be None, an int or "
                    f"'auto', not {workers!r}.")
            workers, _why = self.auto_workers()
            workers = max(1, min(workers, len(blocks)))
        lam = self.scenario.tx_terminal.wavelength_m
        if workers is None or int(workers) <= 1 or len(blocks) < 2:
            out = []
            for b in blocks:
                cols = self._read_block(b, fields=fields)
                if screen_phase is False:
                    cols["screen_phase"] = None
                    cols["screen_phase_pa"] = None
                out.append(_apply_block(fn, cols, self.patch, lam,
                                        b * self.block_size, fields=fields,
                                        compact=compact))
            return np.concatenate(out)[:n]
        from concurrent.futures import ProcessPoolExecutor
        payload = {"map": {"root_dir": self.root_dir, "fn": fn,
                           "patch": self.patch, "lam": lam,
                           "block_size": self.block_size, "fields": fields,
                           "compact": compact,
                           "screen_phase": screen_phase},
                   "boost": bool(boost)}
        got = {}
        with ProcessPoolExecutor(max_workers=int(workers),
                                 initializer=_init_worker,
                                 initargs=(payload,)) as pool:
            for b, part in pool.map(_map_block, blocks):
                got[b] = part
        return np.concatenate([got[b] for b in blocks])[:n]

    def recouple(self, detector, aperture_m=None, obscuration_ratio=None,
                 n_trials=None, workers=None, compact=True):
        """Couple the STORED fields into a detector, with no new propagation.

        The call runs on the square CROP of the stored patch, and it builds the
        clip mask and the fibre mode ONE time for each block. See
        `olb.waveoptics.turbulence.run.recouple` and `_PostTail`.

        Args:
            detector:          an SMF, an MMF, an Aperture, a Camera, or None.
            aperture_m:        the receive aperture diameter, in m. None takes
                               the aperture of the scenario receive terminal.
            obscuration_ratio: the central obscuration. None takes the value of
                               the scenario receive terminal.
            n_trials:          the number of trials. None takes every stored
                               trial.
            workers:           None runs in this process. An int or "auto"
                               opens a process pool (see `map_trials`).
            compact:           True reads on the crop (the default). False
                               reads every trial on the full grid.

        Returns:
            A float array of the coupling efficiency of each trial.
        """
        a, o = self._receive_optics(aperture_m, obscuration_ratio)
        lam = self.scenario.tx_terminal.wavelength_m
        return self.map_trials(_CoupleTrials(detector, a, o, lam,
                                             compact=compact),
                               n_trials=n_trials, workers=workers,
                               screen_phase=False, compact=compact)

    def recouple_compensated(self, compensation, detector, aperture_m=None,
                             obscuration_ratio=None, n_trials=None,
                             source=None, workers=None, compact=True):
        """Correct the STORED fields, then couple them into a detector.

        This is the post-hoc perfect-AO twin of `recouple`. It removes the
        first N Noll modes of each stored trial over the receive aperture, then
        it couples the corrected field. So a stored campaign gives the fade of
        ANY compensation stack, with no new propagation. See
        olb.waveoptics.turbulence.run.recouple_compensated.

        The modal basis and the slope reconstructor are built ONE time for each
        block, on the crop. That is the large cost of this read.

        Args:
            compensation:      a sequence of TipTilt and AO stages.
            detector:          an SMF, an MMF, an Aperture, a Camera, or None.
            aperture_m:        the receive aperture diameter, in m. None takes
                               the aperture of the scenario receive terminal.
            obscuration_ratio: the central obscuration. None takes the value of
                               the scenario receive terminal.
            n_trials:          the number of trials. None takes every stored
                               trial.
            source:            "screens", "slopes", "gtilt", or None. None
                               follows the channel family: "screens" for a
                               space link (the slab starts from a plane wave, so
                               the summed screen phase IS the sensed wavefront)
                               and "slopes" for a terrestrial link. "gtilt"
                               senses the tilt from the far-field centroid; it
                               is tilt-only, pairs with a TipTilt() stack, and
                               is an opt-in model of a centroid tracker (the
                               validation showed it does not beat the slopes).
                               See
                               olb.waveoptics.turbulence.run.recouple_compensated.
            workers:           None runs in this process. An int or "auto"
                               opens a process pool (see `map_trials`).
            compact:           True reads on the crop (the default). False
                               reads every trial on the full grid.

        Returns:
            A float array of the coupling efficiency of each trial.

        Raises:
            ValueError: source="screens" and the campaign stored no screen
                        phase, or the stack removes no mode.
        """
        a, o = self._receive_optics(aperture_m, obscuration_ratio)
        lam = self.scenario.tx_terminal.wavelength_m
        if source is None:
            source = ("screens" if hasattr(self.scenario, "ground")
                      else "slopes")
        fn = _CompensateTrials(compensation, detector, a, o, lam, source,
                               compact=compact)
        return self.map_trials(fn, n_trials=n_trials, workers=workers,
                               screen_phase=(source == "screens"),
                               compact=compact)

    def _point_ahead_ground(self, where):
        """Give the ground terminal of a point-ahead read, or raise.

        Args:
            where: the caller name, for the error message.

        Returns:
            The ground Terminal.

        Raises:
            ValueError: the campaign made no point-ahead pass, or the scenario
                        is not a space scenario with a ground transmitter.
        """
        if self.point_ahead_rad is None:
            raise ValueError(
                f"{where}: this campaign made no point-ahead pass. Make it "
                "with point_ahead_rad, and run it again.")
        return _uplink_ground(self.scenario, where)

    def recouple_point_ahead(self, compensation, *, source=None, n_trials=None,
                             workers=None, compact=True):
        """Give the point-ahead uplink overlap of the STORED planes.

        This is the campaign-level twin of
        `olb.waveoptics.turbulence.run.point_ahead_overlap`. It corrects the
        stored point-ahead field of each trial with the BEACON estimate of the
        given stack, and it takes the reciprocity overlap with the ground
        transmit mode. So a stored campaign gives the point-ahead fade of ANY
        compensation stack, with NO new propagation. The answered angles are
        the STORED angles (`Campaign.point_ahead_rad`); for another angle use
        `Campaign.point_ahead`.

        The transmit mode and the vacuum baseline are computed ONE time here,
        and the modal basis builds one time for each block. See
        `_PointAheadTrials`.

        Args:
            compensation: None (NO correction), the string "terminal", or a
                          list of TipTilt and AO stages.
            source:       "screens", "slopes", "gtilt", or None. None follows
                          the channel family, so a space campaign senses the
                          stored summed screen phase.
            n_trials:     the number of trials. None takes every stored trial.
            workers:      None runs in this process. An int or "auto" opens a
                          process pool (see `map_trials`).
            compact:      True reads on the crop (the default). False reads on
                          the full grid.

        Returns:
            A float array of the shape (n_trials, n_angles). The column order
            is the order of `Campaign.point_ahead_rad`.

        Raises:
            ValueError: the campaign made no point-ahead pass, it stores no
                        field, the scenario has no ground transmitter, or
                        source="screens" and the campaign stored no screen
                        phase.
        """
        where = "Campaign.recouple_point_ahead"
        ground = self._point_ahead_ground(where)
        if self.patch is None:
            raise ValueError(
                f"{where}: this campaign stores no field. Make it with a "
                "patch_radius_m, and run it again.")
        cdtype = field_dtype(self.precision)
        psi_full = _transmit_mode_crop(ground, self.grid, self.patch, cdtype,
                                       compact=False)
        psi = _transmit_mode_crop(ground, self.grid, self.patch, cdtype,
                                  compact=compact)
        mask = super_gaussian_boundary(
            self.grid.n, PRESETS[self.preset].boundary_width_frac)
        _F_vac, o_vac = space_vacuum_baseline(
            self.grid, self.plan, ground.wavelength_m, mask, cdtype,
            psi_tx=psi_full)
        if source is None:
            source = ("screens" if hasattr(self.scenario, "ground")
                      else "slopes")
        stack, n_modes = _resolve_compensation(self.scenario, compensation)
        rx = clip_terminal(self.scenario)
        fn = _PointAheadTrials(stack if n_modes > 0 else None, rx.aperture_m,
                               rx.obscuration_ratio, psi, o_vac, source,
                               compact=compact)
        return self.map_trials(
            fn, n_trials=n_trials, workers=workers,
            screen_phase=(source == "screens" and n_modes > 0),
            compact=compact)

    def point_ahead(self, angles, compensation, *, source=None, n_trials=None,
                    workers=None, fft_backend="numpy"):
        """Give the uplink overlap at ANY angle inside the drawn margin.

        This is the campaign-level twin of
        `olb.waveoptics.turbulence.run.point_ahead_regenerate`. It rebuilds the
        screens of each stored trial from the seeds of the campaign, it crops
        them at the window of each asked angle, and it propagates. So one
        stored campaign answers a whole angle SWEEP. The factory options (the
        outer scale, the subharmonics, the generator and the precision) come
        from the MANIFEST, so the regenerated atmosphere is the stored one.

        THE COST, per trial: ONE draw set of the screens, and ONE split step
        for each angle. `recouple_point_ahead` costs no propagation, so use it
        whenever the angle is a stored angle.

        Args:
            angles:       the point-ahead angles: the string "geometry", a
                          float, or a sequence of floats, in rad.
            compensation: None (NO correction), "terminal", or a list of
                          stages.
            source:       "screens" (the space source), "slopes", or None for
                          the family rule.
            n_trials:     the number of trials. None takes every stored trial.
            workers:      None runs in this process. An int or "auto" opens a
                          process pool (see `map_trials`).
            fft_backend:  "numpy" (the default), "scipy" or "cupy".

        Returns:
            A float array of the shape (n_trials, n_angles).

        Raises:
            ValueError: the campaign made no point-ahead pass, a window falls
                        off the drawn screen, or an option is unknown.
        """
        self._point_ahead_ground("Campaign.point_ahead")
        if source is None:
            source = ("screens" if hasattr(self.scenario, "ground")
                      else "slopes")
        record = TurbWaveResult(
            trials=[], grid=self.grid, plan=self.plan, report=None,
            preset=self.preset, seed_entropy=_resolve_seed(self.seed),
            screen_n=self.screen_n, screen_margin_m=self.screen_margin_m,
            point_ahead_rad=self.point_ahead_rad)
        fn = _RegeneratePointAhead(
            record, angles, compensation, self.scenario, self.geometry,
            source=source, fft_backend=fft_backend, precision=self.precision,
            L0_m=self.L0_m, subharmonics=self.subharmonics,
            screen_generator=self.screen_generator)
        # The read needs NO stored field: it makes its own. The stored screen
        # phase still comes through, because the screens route checks the
        # regenerated sum against it.
        return self.map_trials(fn, n_trials=n_trials, workers=workers,
                               fields=False)

    def recollect(self, aperture_m=None, obscuration_ratio=None,
                  n_trials=None, workers=None, compact=True):
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
            workers:           None runs in this process. An int or "auto"
                               opens a process pool (see `map_trials`).
            compact:           True reads on the crop (the default). False
                               reads every trial on the full grid.

        Returns:
            A float array, one value for each trial.
        """
        a, o = self._receive_optics(aperture_m, obscuration_ratio)
        return self.map_trials(_CollectTrials(a, o, compact=compact),
                               n_trials=n_trials, workers=workers,
                               screen_phase=False, compact=compact)

    def _receive_optics(self, aperture_m, obscuration_ratio):
        """Give the receive aperture and obscuration of a post-hoc read.

        None takes the value of the scenario receive terminal. The call also
        tests that the aperture sits inside the stored patch.

        Args:
            aperture_m:        the aperture diameter in m, or None.
            obscuration_ratio: the obscuration ratio, or None.

        Returns:
            The pair (aperture_m, obscuration_ratio), as floats.

        Raises:
            ValueError: the aperture is larger than the stored patch.
        """
        rx = self.scenario.rx_terminal
        a = rx.aperture_m if aperture_m is None else float(aperture_m)
        o = (rx.obscuration_ratio if obscuration_ratio is None
             else float(obscuration_ratio))
        if self.patch is not None:
            _check_aperture(self.patch, a)
        return float(a), float(o)


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
    rootN = tempfile.mkdtemp(prefix="olb_campaign_selfcheck5_")
    root5 = tempfile.mkdtemp(prefix="olb_campaign_selfcheck5_")
    root6 = tempfile.mkdtemp(prefix="olb_campaign_selfcheck6_")
    root7 = tempfile.mkdtemp(prefix="olb_campaign_selfcheck7_")
    root8 = tempfile.mkdtemp(prefix="olb_campaign_selfcheck8_")
    root9 = tempfile.mkdtemp(prefix="olb_campaign_selfcheck9_")
    root10 = tempfile.mkdtemp(prefix="olb_campaign_selfcheck10_")
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

            # The auto pool. The estimate is positive, the count is at least
            # one, and "auto" on a campaign with no missing block is a no-op.
            assert camp.worker_memory_bytes() > 0
            k, why = camp.auto_workers()
            assert k >= 1 and why in ("cpu", "memory"), (k, why)
            assert camp.run(12, workers="auto") == 12
            # A third block through "auto" matches the seeded serial route.
            bs = common["block_size"]
            camp.run(4 * bs, workers="auto", progress=True)
            with np.load(camp._block_path(3)) as za:
                ref = propagate_turbulent_scenario(
                    scn, orbit, n_trials=bs, start_index=3 * bs,
                    seed=common["seed"], preset=camp.preset, grid=camp.grid,
                    plan=camp.plan, patch_radius_m=camp.patch_radius_m,
                    L0_m=camp.L0_m, precision=camp.precision, threader=None)
                assert np.array_equal(za["fields"], ref.fields)
            print(f"  auto pool                      {k} workers ({why})")
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
            assert camp3.patch_radius_m == ground.aperture_m * PATCH_MARGIN_FACTOR, \
                camp3.patch_radius_m
            assert camp3.run(4) == 4
            eta3 = camp3.recouple(ground.detector)
            assert eta3.size == 4 and np.all(eta3 > 0.0), eta3
            big = (camp3.grid.n, camp3.grid.size_m) != (camp.grid.n,
                                                        camp.grid.size_m)

            # ---- 6b. a reopen reads patch_radius from the manifest ----
            # A store written with an explicit NARROW radius must reopen under a
            # None patch_radius_m with NO manifest mismatch, and it must read
            # back the stored radius, not the wider default.
            narrow = ground.aperture_m / 2.0
            campN = Campaign(scn, orbit, rootN, patch_radius_m=narrow, **common)
            assert campN.run(4) == 4
            reopened = Campaign(scn, orbit, rootN, **common)
            assert reopened.patch_radius_m == narrow, reopened.patch_radius_m

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

            # ---- 9. the perfect-AO correction and the screen store ----
            # A corrected campaign is a SEPARATE store: it gets its own key,
            # and its manifest names the stages. A default campaign does not
            # move (assertion 1 above compares the block bytes of a grow).
            from ...terminal import TipTilt
            camp6 = Campaign(scn, orbit, root6, compensation=[TipTilt()],
                             store_screen_phase=True, **common)
            assert camp6.fingerprint != camp.fingerprint, "the stack must key"
            assert camp6.n_modes_corrected == 3, camp6.n_modes_corrected
            assert camp6.run(8) == 8
            with open(os.path.join(root6, MANIFEST_NAME), encoding="utf-8") as fh:
                man6 = _json.load(fh)
            assert "TipTilt" in man6["compensation"], man6["compensation"]
            assert man6["store_screen_phase"] is True
            got6 = camp6.load(8)
            assert got6.n_modes_corrected == 3
            assert got6.screen_phase.shape == (8, camp6.patch.indices.size)
            # A reopen with no compensation raises, and it names the field.
            try:
                Campaign(scn, orbit, root6, store_screen_phase=True, **common)
                raise AssertionError("a dropped stack must raise ValueError")
            except ValueError as exc:
                assert "compensation" in str(exc), str(exc)
            # THE POST-HOC ROUTE. An UNCORRECTED campaign that stored its
            # screen phase gives the corrected coupling of the corrected
            # campaign, with no new propagation.
            camp7 = Campaign(scn, orbit, root7, store_screen_phase=True,
                             **common)
            assert camp7.run(8) == 8
            eta_post = camp7.recouple_compensated([TipTilt()], ground.detector)
            eta_ao = np.array([t.smf_eta for t in got6.trials])
            d_post = float(np.abs(eta_post / eta_ao - 1.0).max())
            assert d_post < 1e-4, (eta_post[:3], eta_ao[:3])
            eta_plain = camp7.recouple(ground.detector)
            assert eta_ao.mean() > eta_plain.mean(), (eta_ao.mean(),
                                                      eta_plain.mean())

            # ---- 10. map_trials, the field reader, and the crop rule ----
            # The crop route and the full-grid route agree to the float
            # rounding level, and a process pool gives the same numbers.
            eta_full = camp7.recouple(ground.detector, compact=False)
            assert np.all(np.abs(eta_plain / eta_full - 1.0) < 1e-6), \
                (eta_plain, eta_full)
            eta_pool = camp7.recouple(ground.detector, workers=2)
            assert np.array_equal(eta_plain, eta_pool), (eta_plain, eta_pool)
            post_full = camp7.recouple_compensated([TipTilt()],
                                                   ground.detector,
                                                   compact=False)
            assert np.all(np.abs(eta_post / post_full - 1.0) < 1e-6), \
                (eta_post, post_full)
            pw_crop = camp7.recollect()
            pw_full = camp7.recollect(compact=False)
            assert np.all(np.abs(pw_crop / pw_full - 1.0) < 1e-6), (pw_crop,
                                                                    pw_full)
            # A scalar-only pass reads no field, and it sees the stored
            # scalars and the row index.
            def _row_of(rec):
                """Give the row index, and read the stored scalars."""
                assert rec.array is None
                assert np.isfinite(rec.scalars["collected_power"])
                return rec.row

            rows = camp7.map_trials(_row_of, fields=False)
            assert np.array_equal(rows, np.arange(8)), rows
            # The public field reader gives the crop and the full grid.
            f_crop = camp7.field(3)
            f_full = camp7.field(3, compact=False)
            _crop = camp7.patch.crop()
            _o, _s = _crop.offset, _crop.side
            assert f_crop.N == _s and f_full.N == camp7.patch.n
            assert np.array_equal(f_crop.field,
                                  f_full.field[_o:_o + _s, _o:_o + _s])
            # ---- 11. the point ahead round-trips through the store ----
            # TWO BLOCKS of a point-ahead campaign must equal the native run,
            # trial for trial, and the manifest must name the three fields.
            up_scn = SpaceScenario(
                ground=ground, space=Terminal(aperture_m=0.30,
                                              wavelength_m=lam),
                direction="uplink", channel=Channel())
            theta = float(np.asarray(orbit.point_ahead_rad).ravel()[0])
            camp8 = Campaign(up_scn, orbit, root8,
                             point_ahead_rad=(0.0, theta),
                             store_screen_phase=True, **common)
            assert camp8.fingerprint != Campaign(
                up_scn, orbit, root9, store_screen_phase=True,
                **common).fingerprint, "a point-ahead angle must key"
            assert camp8.point_ahead_rad == (0.0, theta)
            assert camp8.screen_margin_m > 0.0, camp8.screen_margin_m
            assert camp8.screen_n > camp8.grid.n, (camp8.screen_n,
                                                   camp8.grid.n)
            assert camp8.run(8) == 8
            with open(os.path.join(root8, MANIFEST_NAME),
                      encoding="utf-8") as fh:
                man8 = _json.load(fh)
            assert man8["point_ahead_rad"] == [0.0, theta], man8["point_ahead_rad"]
            assert man8["screen_margin_m"] == camp8.screen_margin_m
            assert man8["screen_n"] == camp8.screen_n
            # A reopen with no explicit margin reads the stored value.
            re8 = Campaign(up_scn, orbit, root8, point_ahead_rad=(0.0, theta),
                           store_screen_phase=True, **common)
            assert re8.screen_margin_m == camp8.screen_margin_m
            assert re8.fingerprint == camp8.fingerprint
            got8 = camp8.load(8)
            native8 = propagate_turbulent_scenario(
                up_scn, orbit, n_trials=8, seed=common["seed"],
                preset=common["preset"], grid=camp8.grid, plan=camp8.plan,
                patch_radius_m=camp8.patch_radius_m,
                point_ahead_rad=(0.0, theta), screen_margin_m=camp8.screen_margin_m,
                store_screen_phase=True, L0_m=camp8.L0_m,
                precision=camp8.precision)
            assert got8.point_ahead_rad == (0.0, theta)
            for a, b in zip(got8.trials, native8.trials):
                assert a.eta_turb_pa == b.eta_turb_pa, (a.eta_turb_pa,
                                                        b.eta_turb_pa)
                assert a.eta_turb_pa[0] == a.eta_turb, a.eta_turb_pa
            # The stored planes come back in TRIAL order, one plane per angle.
            assert got8.fields_pa.shape == (2, 8, camp8.patch.indices.size)
            assert np.array_equal(got8.fields_pa, native8.fields_pa)
            assert np.array_equal(got8.screen_phase_pa,
                                  native8.screen_phase_pa)
            # fields=False drops the point-ahead FIELD and keeps the scalars.
            light8 = camp8.load(8, fields=False)
            assert light8.fields_pa is None
            assert light8.trials[0].eta_turb_pa == got8.trials[0].eta_turb_pa
            assert light8.screen_phase_pa is not None
            # map_trials hands the per-angle crops to a callable.
            def _shapes(rec):
                """Give the plane count and the crop side of one trial."""
                assert rec.screen_phase_pa.shape[0] == 2
                assert len(rec.arrays_pa) == 2
                assert rec.arrays_pa[0].shape == rec.array.shape
                return float(np.abs(rec.arrays_pa[1]).sum())
            assert np.all(camp8.map_trials(_shapes) > 0.0)

            # ---- 11a. the point-ahead POST-HOC reads ----
            # THE STORED-PLANE ROUTE. A CORRECTED campaign of the same seed
            # gives the reference numbers, and recouple_point_ahead must
            # reproduce them from the UNCORRECTED planes of camp8.
            from ...terminal import AO
            camp10 = Campaign(up_scn, orbit, root10,
                              point_ahead_rad=(0.0, theta),
                              store_screen_phase=True,
                              compensation=[AO(n_modes=10)], **common)
            assert camp10.run(8) == 8
            pa_run = np.array([t.eta_turb_pa
                               for t in camp10.load(8).trials])
            pa_post = camp8.recouple_point_ahead([AO(n_modes=10)])
            d_pa_post = float(np.abs(pa_post / pa_run - 1.0).max())
            assert pa_post.shape == pa_run.shape, pa_post.shape
            assert d_pa_post < 1e-5, (pa_post[0], pa_run[0], d_pa_post)
            # An UNCORRECTED read gives the stored numbers back.
            pa_stored = np.array([t.eta_turb_pa for t in got8.trials])
            pa_none = camp8.recouple_point_ahead(None)
            d_pa_none = float(np.abs(pa_none / pa_stored - 1.0).max())
            assert d_pa_none < 1e-5, (pa_none[0], pa_stored[0], d_pa_none)
            # THE REGENERATION ROUTE at a STORED angle is bit-identical.
            pa_regen = camp8.point_ahead([theta], None, n_trials=4)
            assert np.array_equal(pa_regen[:, 0], pa_stored[:4, 1]), \
                (pa_regen[:, 0], pa_stored[:4, 1])

            # ---- 11b. an OLD block and an OLD manifest still read ----
            # Delete the three point-ahead keys of block 0 and of the manifest,
            # and the store must open and load as a campaign with no point
            # ahead. This is the shape of every campaign stored before today.
            with np.load(camp8._block_path(0)) as z8:
                kept = {k: z8[k] for k in z8.files
                        if k not in ("eta_turb_pa", "fields_pa",
                                     "screen_phase_pa")}
            np.savez(camp8._block_path(0), **kept)
            man_old = {k: v for k, v in man8.items()
                       if k not in ("point_ahead_rad", "screen_margin_m",
                                    "screen_n")}
            man_old["fingerprint"] = Campaign(
                up_scn, orbit, root9, store_screen_phase=True,
                **common).fingerprint
            with open(os.path.join(root8, MANIFEST_NAME), "w",
                      encoding="utf-8") as fh:
                _json.dump(man_old, fh)
            old8 = Campaign(up_scn, orbit, root8, store_screen_phase=True,
                            **common)
            assert old8.point_ahead_rad is None
            assert old8.screen_margin_m == 0.0
            assert old8.screen_n == old8.grid.n
            loaded_old = old8.load(4)     # block 0 only: the trimmed file.
            assert all(t.eta_turb_pa is None for t in loaded_old.trials)
            assert loaded_old.fields_pa is None

            # fields=False KEEPS the stored screen phase.
            light = camp7.load(8, fields=False)
            assert light.fields is None
            assert light.screen_phase is not None
            assert light.screen_phase.shape == (8, camp7.patch.indices.size)

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
        print(f"  SMF eta, uncorrected    {eta_plain.mean():11.6f}")
        print(f"  SMF eta, TipTilt in run {eta_ao.mean():11.6f}")
        print(f"  SMF eta, TipTilt posthoc {eta_post.mean():10.6f} "
              f"(worst relative error {d_post:.1e})")
        print(f"  point ahead, screen     {camp8.grid.n:11d} -> "
              f"{camp8.screen_n} px "
              f"(margin {camp8.screen_margin_m:.3f} m)")
        print(f"  eta_turb, beacon        {got8.trials[0].eta_turb:11.6f}")
        print(f"  eta_turb, ahead         "
              f"{got8.trials[0].eta_turb_pa[1]:11.6f}")
        print(f"  ahead AO(10), post hoc  {pa_post[0, 1]:11.6f} "
              f"(worst relative error {d_pa_post:.1e})")
        print(f"  ahead, regenerated       bit-identical "
              f"({d_pa_none:.1e} on the stored read)")
        print("")
        print(f"(elapsed {time.time() - t_start:.1f} s)")
        print("self-check passed")
    finally:
        for d in (root, root2, root3, root4, root5, root6, root7, root8,
                  root9, root10, rootN):
            shutil.rmtree(d, ignore_errors=True)
