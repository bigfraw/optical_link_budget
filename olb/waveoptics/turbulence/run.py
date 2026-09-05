"""The turbulent trial runner: one atmosphere snapshot for each seed.

The module puts the pieces of the split-step layer together. It sizes the grid
(sampling.py), it makes one screen stack for each trial (screens.py), it moves
the field through that stack (splitstep.py), and it reads the receive plane.

The runner builds NO Term and it changes NO budget. It gives a set of
independent SNAPSHOTS of the atmosphere. There is no time axis. See
temporal.py for the planned frozen-flow extension.

TWO CASES.

- TERRESTRIAL. The runner launches the transmit beam of the near terminal, and
  it propagates that beam along the horizontal path. The launch recipe is the
  recipe of olb.waveoptics.run.propagate_scenario, and the runner imports its
  helpers. So the vacuum limit of this module IS the vacuum module.

- SPACE. The gridded path is the DOWNLINK atmosphere slab only. The satellite
  is outside the atmosphere, so a unit PLANE WAVE enters at the top of the
  slab. A downlink reads the collected power at the ground. An uplink reads
  the SAME field through reciprocity: the uplink flux on the satellite is the
  overlap of the received downlink field with the ground transmit mode. See
  Shapiro, DOI 10.1364/JOSA.61.000492.

Sources:
- Shapiro, Reciprocity of the turbulent atmosphere, DOI 10.1364/JOSA.61.000492.
  The uplink-downlink reciprocity that gives eta_turb.
- Schmidt, Numerical Simulation of Optical Wave Propagation with Examples in
  MATLAB, DOI 10.1117/3.866274, Ch. 9. The split-step method and the absorbing
  boundary.
- Andrews and Phillips, Laser Beam Propagation through Random Media, 2nd ed.,
  DOI 10.1117/3.626196, Ch. 8 and Ch. 12. The scintillation that the self-check
  compares against.
"""

import time
import warnings
from dataclasses import dataclass

import numpy as np

from ...beam import virtual_waist
from ...terminal import Aperture, Camera, SMF, MMF
from ..field import Begin, Power, field_dtype
from ..mmf import mmf_coupling_efficiency
from ..propagators import GForvard
from ..run import _clip, _launch_aperture, _normalised_gauss, _smf_eta
from ..sources import GaussBeam
from .sampling import PRESETS, turbulent_grid
from .screens import ScreenFactory, phase_screen
from .splitstep import split_step, super_gaussian_boundary


def _progress_bar(progress, total, desc):
    """Make a tqdm bar over the trials, or None.

    tqdm is an OPTIONAL import. progress=True with no tqdm gives no bar and a
    warning, so a progress request never stops a run.

    Args:
        progress: True for a bar, False (or None) for no bar.
        total:    the number of trials, the bar length.
        desc:     the bar label.

    Returns:
        A tqdm instance, or None.
    """
    if not progress:
        return None
    try:
        from tqdm.auto import tqdm
    except ImportError:
        warnings.warn("progress=True needs the tqdm package. Install tqdm, or "
                      "pass progress=False. The run goes on with no bar.")
        return None
    return tqdm(total=total, desc=desc)


def _mmf_focal_length(detector, aperture_m, lam):
    """Give the focal length of the multimode-fibre coupling optic, in m.

    An explicit MMF.focal_length_m wins. Else MMF.optimal_focus matches the spot
    to the core through the a = 1.12 spot-to-core parameter. Source: Shaklan and
    Roddier, Appl. Opt. 27 (1988) 2334, DOI 10.1364/AO.27.002334. This is the
    SAME rule as olb.models.coupling.terrestrial._mmf_focal_length.

    Args:
        detector:   the MMF detector.
        aperture_m: the receive aperture diameter, in m.
        lam:        the wavelength, in m.

    Returns:
        The focal length, in m.

    Raises:
        ValueError: the detector sets no focal length and no optimal_focus.
    """
    if detector.focal_length_m is not None:
        return float(detector.focal_length_m)
    if detector.optimal_focus:
        return (np.pi * (aperture_m / 2.0) * detector.core_radius_m
                / (lam * 1.12))
    raise ValueError(
        "the MMF detector needs a focal length to focus the field. "
        "Set MMF.focal_length_m, or set MMF.optimal_focus=True to "
        "match the spot to the core.")


def _detector_eta(detector, collected, aperture_m, lam):
    """Give the coupling efficiency of one detector on the collected field.

    The field `collected` is the receive-plane field AFTER the aperture clip.
    Each efficiency is a ratio of the detector power to the collected power, so
    it is power-normalised. A beamsplitter scales the field of an arm by a
    constant, and a constant cancels in a ratio. So this efficiency does NOT
    change with the splitter fraction, and one field serves every arm. The
    fraction is a separate fixed dB Term (see olb.models.splitter).

    Args:
        detector:   an SMF, an MMF, an Aperture, a Camera, or None.
        collected:  the clipped receive-plane Field.
        aperture_m: the receive aperture diameter, in m.
        lam:        the wavelength, in m.

    Returns:
        The efficiency as a float, or None when the detector has no coupling
        model (a Camera, or no detector).

    Raises:
        ValueError: the detector type is unknown, or an MMF sets no focal
                    length.
    """
    if detector is None or isinstance(detector, Camera):
        # A Camera is a DIAGNOSTIC front end. It measures the spot shape and the
        # spot position, and no budget builds a coupling Term for it. See
        # olb.terminal.Camera and olb.waveoptics.camera.
        return None
    if isinstance(detector, SMF):
        # The overlap of the field with the back-projected fibre mode. The fibre
        # tip sits at z = f + defocus_m, which is a quadratic pupil phase, the
        # same convention as the multimode leg. See olb.waveoptics.smf and
        # olb.waveoptics.mmf.defocus_phase.
        return _smf_eta(detector, collected, aperture_m, lam)
    if isinstance(detector, MMF):
        # A non-focal-plane detector (z = f + defocus_m). The AXIAL displacement
        # grows the spot (defocus_m); the field already carries the turbulent
        # tilt. At the focal plane (defocus_m = 0) this is the plain focal-plane
        # coupling. See olb.waveoptics.mmf and olb.models.coupling.terrestrial.
        f_mmf = _mmf_focal_length(detector, aperture_m, lam)
        return float(mmf_coupling_efficiency(
            collected, aperture_m, detector.core_radius_m, f_mmf,
            numerical_aperture=detector.numerical_aperture,
            defocus_m=detector.defocus_m))
    if isinstance(detector, Aperture):
        # A power-in-bucket detector takes ALL the collected power, so its
        # efficiency against the aperture clip is exactly 1.0. The real loss of
        # this arm is the aperture capture, which collected_power already holds.
        return 1.0
    raise ValueError(
        f"_detector_eta: unknown detector type {type(detector).__name__}. Use "
        "an Aperture, an SMF, an MMF, or a Camera.")


@dataclass(frozen=True)
class TurbTrial:
    """One atmosphere snapshot.

    Attributes:
        collected_power: the power inside the receive aperture, as a fraction
                         of the input power. The terrestrial case divides by
                         the launched power after the transmit clip, so it
                         holds the geometric spread too. The space case divides
                         by the VACUUM baseline on the same grid, so it holds
                         the turbulence penalty only. Its vacuum limit is 1.0.
        smf_eta:         the single-mode-fibre coupling efficiency. It is None
                         when the receive terminal has no SMF detector.
        eta_turb:        the reciprocity overlap ratio of an uplink, against
                         the free-space baseline. It is None for a downlink and
                         for a terrestrial case.
        seed_key:        the pair (seed_entropy, trial_index).
        wall_time_s:     the time of the trial, in s. It holds the screen
                         generation and the propagation.
        mmf_eta:         the multimode-fibre (light-bucket) coupling efficiency.
                         It is the encircled energy of the focused spot inside
                         the core, and the turbulent tilt walks the spot off the
                         on-axis core on its own. It is None when the receive
                         terminal has no MMF detector.
        detector_etas:   the coupling efficiency of each detector of the
                         `detectors` argument, in that order. It is None when
                         the caller gives no `detectors`. A Camera arm holds
                         None, because a Camera has no coupling model.
    """

    collected_power: float
    smf_eta: float
    eta_turb: float
    seed_key: tuple
    wall_time_s: float
    mmf_eta: float = None
    detector_etas: tuple = None


@dataclass(frozen=True)
class FieldPatch:
    """The mask that selects the stored receive-plane pixels.

    The patch is a disc at the centre of the grid. It uses the pixel-centre
    convention of olb.waveoptics.sources.CircAperture, so a patch of the radius
    D/2 holds exactly the pixels that a CircAperture of the diameter D keeps.

    Attributes:
        radius_m: the radius of the disc, in m.
        n:        the number of grid pixels along one side.
        pixel_m:  the distance between two pixels, in m.
        indices:  the flat indices of the disc pixels, an int32 array. The
                  order is the C order of the n x n grid.
    """

    radius_m: float
    n: int
    pixel_m: float
    indices: np.ndarray


def _field_patch(grid, radius_m):
    """Make the FieldPatch of one grid.

    The mask uses the coordinate construction of
    olb.waveoptics.sources.CircAperture: the pixel (int(n/2), int(n/2)) is the
    axis, and the mask keeps the pixels with r^2 <= radius^2.

    Args:
        grid:     the GridSpec.
        radius_m: the radius of the disc, in m.

    Returns:
        A FieldPatch.

    Raises:
        ValueError: the radius does not fit on the grid.
    """
    half = grid.size_m / 2.0
    if radius_m > half:
        raise ValueError(
            f"propagate_turbulent_scenario: patch_radius_m ({radius_m:.4g} m) "
            f"is larger than half the grid side ({half:.4g} m). The patch must "
            "fit on the grid. Use a smaller radius, or a wider grid.")
    n = int(grid.n)
    c = int(n / 2)
    Y, X = np.mgrid[:n, :n]
    dist_sq = ((X - c) * grid.pixel_m) ** 2 + ((Y - c) * grid.pixel_m) ** 2
    mask = dist_sq <= radius_m ** 2
    return FieldPatch(radius_m=float(radius_m), n=n,
                      pixel_m=float(grid.pixel_m),
                      indices=np.flatnonzero(mask).astype(np.int32))


@dataclass(frozen=True)
class TurbWaveResult:
    """The result of a set of turbulent trials.

    NOTE. The record holds the per-trial SCALARS, and, when the caller asks for
    it, the OPTIONAL masked receive field (`fields` and `patch`). The field
    capture is an owner decision of 2026-09-04: a large campaign must recouple
    a stored field to a new detector, without a new propagation. The scalars
    stay exactly as they are, and a budget never reads the fields.

    Attributes:
        trials:       a list of TurbTrial, one for each trial.
        grid:         the GridSpec that the trials used.
        plan:         the ScreenPlan that the trials used.
        report:       the SamplingReport of the grid. It is None when the
                      caller gives its own grid and plan.
        preset:       the name of the quality preset.
        seed_entropy: the integer that seeds every trial. Give it back to
                      repeat the same set.
        fields:       the stored receive-plane field on the patch, a complex64
                      array of the shape (n_trials, n_patch). The row order is
                      the trial order. It is None when the caller asks for no
                      patch.
        patch:        the FieldPatch of those columns, or None.
    """

    trials: list
    grid: object
    plan: object
    report: object
    preset: str
    seed_entropy: int
    fields: np.ndarray = None
    patch: FieldPatch = None


def folded_terrestrial(*args, **kwargs):
    """PLANNED, NOT BUILT. The double pass of a corner-cube retroreflector."""
    raise NotImplementedError(
        "the folded (retroreflected) terrestrial double pass is not built. "
        "The two passes share the same screens, so they are correlated. That "
        "correlation is the physics of the link, and it needs its own design.")


def _resolve_seed(seed):
    """Turn an int, a numpy Generator, or None into one integer.

    See the numpy SeedSequence documentation: one entropy value plus a
    spawn_key gives an independent, repeatable stream for each (trial, screen)
    pair.
    """
    if seed is None:
        return int(np.random.SeedSequence().entropy)
    if isinstance(seed, np.random.Generator):
        return int(seed.integers(2 ** 63))
    return int(seed)


def _screen_seed(entropy, trial, screen):
    """Give the integer seed of one screen of one trial.

    The spawn_key holds the trial index and the screen index. So trial k is
    bit-identical, and the trial count does not change it. See the numpy
    SeedSequence documentation.
    """
    ss = np.random.SeedSequence(entropy=entropy, spawn_key=(trial, screen))
    return int(ss.generate_state(1)[0])


def _screen_builder(screen_generator, grid, L0_m, subharmonics,
                    dtype=np.complex128):
    """Give a function build(seed_int, r0_m) -> phase screen.

    The two generators give DIFFERENT random draws for the same integer seed.
    Both draw the same random field, so the statistics agree, but a screen is
    not bit-identical between them. The default is "olb", the fast cached
    generator. Pass "aotools" to reproduce an old aotools run bit-identically.

    - "olb":     a cached ScreenFactory (screens.py). It builds the filter and
                 the separable subharmonic basis one time for the grid, then it
                 scales them for each screen. It is faster; the subharmonic sum
                 is a matrix product, not 27 full-grid exponentials. The broad
                 validity pass validated it against "aotools" and the analytic
                 index (validation/waveoptics_speed/generator_validation.py).
    - "aotools": one aotools call for each screen. It is the reference path, and
                 it keeps an old aotools run bit-identical.

    Args:
        screen_generator: "aotools" or "olb".
        grid:             the GridSpec.
        L0_m:             the outer scale of the screens, in m.
        subharmonics:     True adds the three subharmonic levels.
        dtype:            the complex type of the field. numpy.complex64 makes
                          float32 screens, so the screen and the field carry
                          the same number of bytes.

    Returns:
        A callable build(seed_int, r0_m) that gives one n x n phase screen.

    Raises:
        ValueError: the generator name is unknown.
    """
    single = dtype == np.complex64
    if screen_generator == "aotools":
        def build(seed_int, r0_m):
            scr = phase_screen(r0_m, grid.n, grid.pixel_m, L0_m=L0_m,
                               seed=seed_int, subharmonics=subharmonics)
            return scr.astype(np.float32) if single else scr
        return build
    if screen_generator == "olb":
        factory = ScreenFactory(grid.n, grid.pixel_m, L0_m=L0_m,
                                subharmonics=subharmonics,
                                dtype=np.float32 if single else np.float64)

        def build(seed_int, r0_m):
            return factory.make(r0_m, np.random.default_rng(seed_int))
        return build
    raise ValueError(
        f"propagate_turbulent_scenario: screen_generator must be 'aotools' or "
        f"'olb', not {screen_generator!r}.")


def _start_field(scenario, grid, lam, is_space, dtype=np.complex128):
    """Make the field that enters the split step.

    The space slab starts from a unit PLANE WAVE that fills the grid: the
    satellite is outside the atmosphere. A terrestrial path starts from the
    clipped transmit beam of the near terminal, on the exact launch recipe of
    olb.waveoptics.run.propagate_scenario.

    Args:
        scenario: a SpaceScenario or a TerrestrialScenario.
        grid:     the GridSpec.
        lam:      the wavelength, in m.
        is_space: True for a space slab, False for a terrestrial path.
        dtype:    the complex type of the field.

    Returns:
        A Field.
    """
    if is_space:
        return Begin(grid.size_m, lam, grid.n, dtype=dtype)
    tx = scenario.tx_terminal
    t = tx.transmitter
    w_v, offset = virtual_waist(t.waist_m, t.divergence_rad, lam)
    F0 = _normalised_gauss(GaussBeam(Begin(grid.size_m, lam, grid.n,
                                           dtype=dtype), w_v))
    if offset > 0:
        F0 = GForvard(F0, offset)
    return _clip(F0, *_launch_aperture(tx))


def _ground_transmit_mode(ground, grid, dtype=np.complex128):
    """Make the normalised transmit mode of the ground terminal.

    The recipe is the launch recipe of olb.waveoptics.run.propagate_scenario:
    a Gaussian at the virtual waist, an offset propagation, and the launch
    aperture clip. The function then scales the array so that
    sum(|psi|^2) = 1.0.

    Args:
        ground: the ground Terminal. It needs a Transmitter.
        grid:   the GridSpec.
        dtype:  the complex type of the field.

    Returns:
        An N x N complex array.
    """
    t = ground.transmitter
    lam = ground.wavelength_m
    w_v, offset = virtual_waist(t.waist_m, t.divergence_rad, lam)
    F = _normalised_gauss(GaussBeam(Begin(grid.size_m, lam, grid.n,
                                          dtype=dtype), w_v))
    if offset > 0:
        F = GForvard(F, offset)
    aperture_m, obscuration = _launch_aperture(ground)
    psi = _clip(F, aperture_m, obscuration).field
    return psi / np.sqrt((np.abs(psi) ** 2).sum())


def clip_terminal(scenario):
    '''
    Give the terminal whose aperture clips the propagated field.

    A SPACE scenario ALWAYS propagates the downlink slab, so its field arrives
    at the GROUND terminal in both directions: a downlink receives there, and
    an uplink reads the same field through the Shapiro reciprocity overlap
    with the ground transmit mode (DOI 10.1364/JOSA.61.000492). So the ground
    aperture is the physical plane of the field, whatever the direction. A
    TERRESTRIAL scenario clips at its receive terminal.

    The runner, the campaign patch and the campaign sizing copy all read this
    ONE rule, so the stored field and the clip never disagree.

    Parameters:
        scenario : SpaceScenario or TerrestrialScenario

    Returns:
        Terminal
    '''
    return scenario.ground if hasattr(scenario, "ground") else scenario.rx_terminal


def propagate_turbulent_scenario(scenario, geometry, *, n_trials=1, seed=None,
                                 preset="standard", grid=None, plan=None,
                                 cn2=None, hs=None, cn2_profile=None,
                                 h_top_m=None, L0_m=np.inf,
                                 subharmonics=True, threader=None,
                                 screen_generator="olb", progress=False,
                                 detectors=None, start_index=0,
                                 patch_radius_m=None, precision="single"):
    """Run a set of turbulent split-step trials for one scenario.

    Each trial makes a new screen stack and moves one field through it. The
    trials are independent snapshots. They carry no time axis.

    THE TRIALS THREAD. The trials share only read-only state (the grid, the
    plan, the boundary mask, the reference power, the transmit mode), and each
    one seeds its own screens, so they run across threads with no lock. Give a
    Threader to run them at the same time. The per-trial wall_time_s then holds
    the time of the trial WITH the thread contention, so its sum runs past the
    real wall time. See olb.waveoptics.Threader.

    THE BOUNDARY MASK IS ALWAYS ON. The subharmonic content of a screen is not
    periodic on the grid, and Forvard is periodic. Without the mask the wrap
    corrupts the statistics. The sizer keeps the receive aperture inside the
    untouched interior of the mask, and this function checks that.

    SEEDS. Trial k is bit-identical for one seed, and the trial count does not
    change it. So a longer run repeats the trials of a shorter run. The run
    covers the trials k = start_index .. start_index + n_trials - 1, and
    TurbTrial.seed_key holds the TRUE k. So a run of 500 trials equals the
    concatenation of (start_index=0, n_trials=200) and (start_index=200,
    n_trials=300), trial for trial and bit for bit. A campaign therefore
    computes its blocks in any order, or on more than one process.

    THE STORED FIELD. patch_radius_m stores the receive-plane field on a disc
    at the centre of the grid, BEFORE the receive-aperture clip. The scalars do
    not change. The memory is small: a 1 m patch at a 5 mm pixel pitch is about
    200 x 200 complex64, which is about 320 kB for each trial. Use `recouple`
    and `recollect` to read a stored field with a NEW detector.

    Args:
        scenario:     a SpaceScenario or a TerrestrialScenario.
        geometry:     an object with slant_range_m, and with elevation_deg for
                      a space case. The range must be ONE value.
        n_trials:     the number of snapshots.
        seed:         an int, a numpy Generator, or None for a fresh entropy.
        preset:       the name of a preset in sampling.PRESETS.
        grid:         an optional GridSpec. Give grid AND plan together, or
                      give neither.
        plan:         an optional ScreenPlan.
        cn2:          an optional callable cn2(h) -> the zenith Cn2 at height h
                      [m]. None (with no hs/cn2_profile) integrates the site
                      Hufnagel-Valley profile: the continuous default. Space
                      only. See turbulent_grid.
        hs:           the height grid of a DISCRETE Cn2 profile, in m. Give it
                      to take the legacy array planner. Space only.
        cn2_profile:  the zenith Cn2 profile on hs. Space only.
        h_top_m:      the atmosphere top for the continuous integral, in m.
                      Space only.
        L0_m:         the outer scale of the screens, in m.
        subharmonics: True adds the three subharmonic levels to each screen.
                      Keep it True: the tilt content drives the beam wander,
                      and the uplink overlap reads that wander.
        threader:     an optional olb.waveoptics.Threader. None runs the trials
                      one by one. A Threader runs them across threads, and it
                      keeps the trial order.
        screen_generator: "olb" (the default) or "aotools". The default is the
                      fast cached ScreenFactory of screens.py. "aotools" keeps
                      an old aotools run bit-identical. The two generators give
                      DIFFERENT random draws for the same seed; the statistics
                      agree, so use "olb" for speed.
        progress:     True shows a tqdm bar that advances one step for each
                      finished trial. It needs the optional tqdm package; if
                      tqdm is not installed, the run goes on with no bar and a
                      warning. False (the default) shows no bar. With a threader
                      the bar advances in the finishing order, not the trial
                      order, but the returned trials keep the trial order.
        detectors:    an optional sequence of detector objects, the arms behind
                      a receive beamsplitter. Each trial then gives the coupling
                      efficiency of EVERY arm, in this order, as
                      TurbTrial.detector_etas. ONE run therefore feeds every
                      arm: the clipped receive field is already in memory, and
                      each arm is one more cheap focal-plane calculation on that
                      same array. The `frac` of each detector is IGNORED here.
                      A beamsplitter scales the field of an arm by a constant,
                      and every coupling efficiency is power-normalised, so the
                      efficiency does not change with the split ratio. The
                      fraction is a separate fixed dB Term (see
                      olb.models.splitter). None (the default) keeps the
                      single-detector record, bit for bit.
        start_index:  the index of the FIRST trial. The run covers the trials
                      start_index .. start_index + n_trials - 1. The default 0
                      is the old behaviour.
        patch_radius_m: the radius of the stored receive-field disc, in m. None
                      (the default) stores no field, and the record is bit for
                      bit the old record. A float fills TurbWaveResult.fields
                      and TurbWaveResult.patch.
        precision:    "single" (the default) or "double". "single" runs the
                      whole propagation in complex64, with float32 phase
                      screens and a float32 boundary mask. WHY: a large
                      campaign is memory-bandwidth bound, so half the bytes for
                      each element gives a real speed-up. CAUTION: a
                      single-precision campaign is a DIFFERENT record. The
                      trials are not bit-identical to a double-precision run,
                      and the arithmetic carries about 7 digits, not 16.
                      Validate a single-precision run against a
                      double-precision run of the same seed before a budget
                      reads it. See validation/precision.

    Returns:
        A TurbWaveResult.

    Raises:
        ValueError:         the geometry gives more than one range, only one
                            of grid and plan is given, the patch radius does
                            not fit on the grid, or the precision name is
                            unknown.
        NotImplementedError: the scenario direction is "retro".
    """
    cdtype = field_dtype(precision)
    is_space = hasattr(scenario, "ground")
    if is_space and scenario.direction == "retro":
        raise NotImplementedError(
            "the retro direction is not built. A retroreflected link goes "
            "through the SAME screens two times, so the two passes are "
            "correlated. That is a separate design.")
    if (grid is None) != (plan is None):
        raise ValueError("propagate_turbulent_scenario: give grid AND plan "
                         "together, or give neither.")

    range_m = np.asarray(geometry.slant_range_m, dtype=float)
    if range_m.size != 1:
        raise ValueError(
            f"propagate_turbulent_scenario: the geometry gives {range_m.size} "
            "ranges. Give one range, and loop in the caller.")

    p = PRESETS[preset] if isinstance(preset, str) else preset
    report = None
    if grid is None:
        grid, plan, report = turbulent_grid(
            scenario, geometry, preset=p, cn2=cn2, hs=hs,
            cn2_profile=cn2_profile, h_top_m=h_top_m, L0_m=L0_m)

    lam = scenario.tx_terminal.wavelength_m
    rx = clip_terminal(scenario)
    mask = super_gaussian_boundary(grid.n, p.boundary_width_frac)

    # The receive aperture must sit in the part of the grid that the mask does
    # not touch. The mask is exactly 1.0 inside (1 - width_frac) of the
    # half-side. See splitstep.super_gaussian_boundary.
    r_flat = (1.0 - p.boundary_width_frac) * grid.size_m / 2
    if rx.aperture_m / 2 >= r_flat:
        warnings.warn(
            f"propagate_turbulent_scenario: the receive aperture radius "
            f"({rx.aperture_m / 2:.4g} m) reaches the absorbing band of the "
            f"boundary mask (it starts at {r_flat:.4g} m). The collected "
            f"power is too low. Use a wider grid.")

    # ---- the fixed parts, computed one time ----
    if is_space:
        # THE VACUUM BASELINE. The space case starts from a unit plane wave
        # that fills the grid, so the absorbing mask acts as a soft aperture.
        # Over a 40 km slab that soft edge makes strong Fresnel rings on the
        # axis. Those rings are a property of the GRID, not of the atmosphere.
        # So the reference is the SAME plane wave along the SAME hops through
        # the SAME mask, with FLAT screens. Then the vacuum limit of each
        # output below is exactly 1.0, and every number is a pure turbulence
        # penalty. The flat screens share one array, because Screen() does not
        # change its input.
        F_plane = Begin(grid.size_m, lam, grid.n, dtype=cdtype)
        flat = np.zeros((grid.n, grid.n))
        F_vac = split_step(F_plane, plan.z_m, [flat] * int(plan.z_m.size),
                           plan.z_total_m, boundary=mask)
        p_reference = Power(_clip(F_vac, rx.aperture_m, rx.obscuration_ratio))
        psi_tx = o_vac = None
        if scenario.direction == "uplink":
            psi_tx = _ground_transmit_mode(scenario.ground, grid, dtype=cdtype)
            # The free-space baseline. It puts eta_turb on the same reference
            # as the (w_free/w_st)^2 rescale of olb.turbulence.uplink_flux.
            o_vac = float(np.abs((F_vac.field * np.conj(psi_tx)).sum()) ** 2)
    else:
        F_in = _start_field(scenario, grid, lam, is_space=False, dtype=cdtype)
        p_reference = Power(F_in)

    seed_entropy = _resolve_seed(seed)
    n_screens = int(plan.z_m.size)
    build_screen = _screen_builder(screen_generator, grid, L0_m, subharmonics,
                                   dtype=cdtype)

    # The optional field capture. One mask serves every trial, and each trial
    # writes ONE row. So the threads touch no shared row, and no lock is
    # necessary.
    patch = fields = None
    if patch_radius_m is not None:
        patch = _field_patch(grid, float(patch_radius_m))
        fields = np.empty((int(n_trials), patch.indices.size),
                          dtype=np.complex64)

    def run_one(k):
        """Run trial k. It touches only its own state and read-only setup."""
        t0 = time.perf_counter()
        stack = [build_screen(_screen_seed(seed_entropy, k, j), plan.r0_m[j])
                 for j in range(n_screens)]
        F_start = (Begin(grid.size_m, lam, grid.n, dtype=cdtype) if is_space
                   else F_in)
        F_rx = split_step(F_start, plan.z_m, stack, plan.z_total_m,
                          boundary=mask)

        if patch is not None:
            # The UNCLIPPED field, on the patch only. The clip below is
            # unchanged, so every scalar keeps its value.
            fields[k - start_index] = (
                F_rx.field.ravel()[patch.indices].astype(np.complex64))

        collected = _clip(F_rx, rx.aperture_m, rx.obscuration_ratio)
        collected_power = float(Power(collected) / p_reference)
        # The receive-terminal detector, the single-detector faces. The MMF and
        # the SMF physics live in _detector_eta, so the multi-detector path
        # below reads the SAME code on the SAME field.
        smf_eta = (_detector_eta(rx.detector, collected, rx.aperture_m, lam)
                   if isinstance(rx.detector, SMF) else None)
        mmf_eta = (_detector_eta(rx.detector, collected, rx.aperture_m, lam)
                   if isinstance(rx.detector, MMF) else None)
        # The extra beamsplitter arms. The field is already in memory, so each
        # arm is one more cheap focal-plane calculation on the SAME array.
        detector_etas = (
            None if detectors is None else
            tuple(_detector_eta(d, collected, rx.aperture_m, lam)
                  for d in detectors))
        eta_turb = None
        if is_space and scenario.direction == "uplink":
            # The reciprocity overlap. See Shapiro,
            # DOI 10.1364/JOSA.61.000492. Point-ahead anisoplanatism is NOT
            # modelled: the uplink and the downlink read the same screens.
            o = float(np.abs((F_rx.field * np.conj(psi_tx)).sum()) ** 2)
            eta_turb = o / o_vac
        return TurbTrial(collected_power=collected_power, smf_eta=smf_eta,
                         eta_turb=eta_turb, seed_key=(seed_entropy, k),
                         wall_time_s=time.perf_counter() - t0,
                         mmf_eta=mmf_eta, detector_etas=detector_etas)

    bar = _progress_bar(progress, n_trials, "turbulent trials")
    try:
        ks = range(int(start_index), int(start_index) + int(n_trials))
        if threader is None:
            trials = []
            for k in ks:
                trials.append(run_one(k))
                if bar is not None:
                    bar.update(1)
        else:
            cb = (lambda done, total: bar.update(1)) if bar is not None else None
            trials = threader.map(run_one, ks, progress=cb)
    finally:
        if bar is not None:
            bar.close()

    return TurbWaveResult(trials=trials, grid=grid, plan=plan, report=report,
                          preset=p.name, seed_entropy=seed_entropy,
                          fields=fields, patch=patch)


def propagate_turbulent_field(scenario, geometry, *, seed=0, trial=0,
                              preset="standard", grid=None, plan=None,
                              cn2=None, hs=None, cn2_profile=None,
                              h_top_m=None, L0_m=np.inf,
                              subharmonics=True, screen_generator="olb",
                              precision="single"):
    """Propagate ONE snapshot and give back the complex receive-plane field.

    This is a DIAGNOSTIC entry point, for a picture of the received field. It
    runs a single trial and it returns the field at the receive plane, BEFORE
    the receive-aperture clip. It does NOT extend the scalar TurbTrial record:
    the rich per-trial record stays deferred (see TurbWaveResult). It repeats
    the atmosphere of one trial of propagate_turbulent_scenario: pass the same
    seed and the trial index, and you get that exact snapshot.

    The space case gives the DOWNLINK field at the ground plane. An uplink
    reads the SAME field through reciprocity, so the ground field is the field
    to picture for both directions.

    Args:
        scenario:     a SpaceScenario or a TerrestrialScenario.
        geometry:     an object with slant_range_m, and elevation_deg for a
                      space case. The range must be ONE value.
        seed:         the seed of the run, an int, a Generator, or None.
        trial:        the trial index inside the run, from 0.
        preset:       the name of a preset in sampling.PRESETS.
        grid:         an optional GridSpec. Give grid AND plan together.
        plan:         an optional ScreenPlan.
        cn2:          an optional callable cn2(h) -> the zenith Cn2 at height h
                      [m]. None (with no hs/cn2_profile) integrates the site
                      Hufnagel-Valley profile: the continuous default. Space
                      only. See turbulent_grid.
        hs:           the height grid of a DISCRETE Cn2 profile, in m. Give it
                      to take the legacy array planner. Space only.
        cn2_profile:  the zenith Cn2 profile on hs. Space only.
        h_top_m:      the atmosphere top for the continuous integral, in m.
                      Space only.
        L0_m:         the outer scale of the screens, in m.
        subharmonics: True adds the three subharmonic levels to each screen.
        screen_generator: "olb" (the default) or "aotools". See
                      propagate_turbulent_scenario. The two give different draws
                      for the same seed.
        precision:    "single" (the default) or "double". "single" runs the
                      propagation in complex64, with float32 screens. WHY: half
                      the bytes for each element, which a memory-bandwidth
                      bound run feels. CAUTION: a single-precision snapshot is
                      not bit-identical to the double-precision snapshot of the
                      same seed. Validate before a budget reads it. See
                      validation/precision.

    Returns:
        A tuple (F_rx, grid, plan). F_rx is the receive-plane Field.

    Raises:
        ValueError:          the geometry gives more than one range, only one
                             of grid and plan is given, or the precision name
                             is unknown.
        NotImplementedError: the scenario direction is "retro".
    """
    cdtype = field_dtype(precision)
    is_space = hasattr(scenario, "ground")
    if is_space and scenario.direction == "retro":
        raise NotImplementedError(
            "the retro direction is not built. See "
            "propagate_turbulent_scenario.")
    if (grid is None) != (plan is None):
        raise ValueError("propagate_turbulent_field: give grid AND plan "
                         "together, or give neither.")
    range_m = np.asarray(geometry.slant_range_m, dtype=float)
    if range_m.size != 1:
        raise ValueError(
            f"propagate_turbulent_field: the geometry gives {range_m.size} "
            "ranges. Give one range.")

    p = PRESETS[preset] if isinstance(preset, str) else preset
    if grid is None:
        grid, plan, _ = turbulent_grid(scenario, geometry, preset=p, cn2=cn2,
                                       hs=hs, cn2_profile=cn2_profile,
                                       h_top_m=h_top_m, L0_m=L0_m)

    lam = scenario.tx_terminal.wavelength_m
    mask = super_gaussian_boundary(grid.n, p.boundary_width_frac)
    entropy = _resolve_seed(seed)
    build_screen = _screen_builder(screen_generator, grid, L0_m, subharmonics,
                                   dtype=cdtype)
    stack = [build_screen(_screen_seed(entropy, trial, j), plan.r0_m[j])
             for j in range(int(plan.z_m.size))]
    F_start = _start_field(scenario, grid, lam, is_space, dtype=cdtype)
    F_rx = split_step(F_start, plan.z_m, stack, plan.z_total_m, boundary=mask)
    return F_rx, grid, plan


def _rebuilt_fields(result, aperture_m, trials):
    """Give the stored trials back as FULL-GRID fields, one at a time.

    The generator scatters the stored patch values into a zero array of the
    FULL grid. A crop would change the zero padding, and the focal-plane pixel
    scale of a fibre coupling reads that padding. So the reconstruction keeps
    the full grid, and the coupling value equals the in-run value.

    The generator makes ONE grid at a time. It never stacks them, so the memory
    holds one field only.

    Args:
        result:     a TurbWaveResult with a stored patch.
        aperture_m: the receive aperture diameter, in m.
        trials:     a sequence of trial row indices, or None for every row.

    Yields:
        The (row, N x N complex array) pair of each selected trial.

    Raises:
        ValueError: the result holds no field, or the aperture is larger than
                    the stored patch.
    """
    if result.fields is None or result.patch is None:
        raise ValueError(
            "this TurbWaveResult holds no field. Run "
            "propagate_turbulent_scenario with patch_radius_m to store one.")
    patch = result.patch
    if aperture_m / 2.0 > patch.radius_m:
        raise ValueError(
            f"the receive aperture radius ({aperture_m / 2.0:.4g} m) is larger "
            f"than the stored patch radius ({patch.radius_m:.4g} m). The "
            "aperture must sit inside the patch. Store a wider patch.")
    rows = range(result.fields.shape[0]) if trials is None else trials
    for row in rows:
        full = np.zeros(patch.n * patch.n, dtype=np.complex128)
        full[patch.indices] = result.fields[row]
        yield row, full.reshape(patch.n, patch.n)


def _patch_field(patch, array, lam):
    """Wrap a full-grid array as a Field on the grid of the patch."""
    F = Begin(patch.n * patch.pixel_m, lam, patch.n)
    F.field = array
    return F


def recouple(result, detector, aperture_m, obscuration_ratio, lam, *,
             trials=None):
    """Couple a STORED receive field into a detector, after the run.

    The function rebuilds the receive-plane field of each stored trial, it
    clips that field at the receive aperture, and it gives the coupling
    efficiency of the detector. So a campaign tries a new detector, a new
    focal length or a new defocus with NO new propagation.

    The physics is the physics of the run: the function calls the same
    `_detector_eta` on the same clipped field.

    Args:
        result:            a TurbWaveResult with a stored patch.
        detector:          an SMF, an MMF, an Aperture, a Camera, or None.
        aperture_m:        the receive aperture diameter, in m.
        obscuration_ratio: the central obscuration of that aperture.
        lam:               the wavelength, in m.
        trials:            an optional sequence of trial row indices. None
                           takes every stored trial.

    Returns:
        A float array of the coupling efficiency of each selected trial. A
        detector with no coupling model (a Camera, or None) gives NaN.

    Raises:
        ValueError: the result holds no field, or the aperture is larger than
                    the stored patch.
    """
    out = []
    for _, array in _rebuilt_fields(result, aperture_m, trials):
        F = _patch_field(result.patch, array, lam)
        collected = _clip(F, aperture_m, obscuration_ratio)
        eta = _detector_eta(detector, collected, aperture_m, lam)
        out.append(np.nan if eta is None else float(eta))
    return np.array(out, dtype=float)


def recollect(result, aperture_m, obscuration_ratio, *, trials=None):
    """Give the collected power of each STORED trial, in grid units.

    The value is Power(clipped) of the rebuilt field. It is NOT normalised:
    the runner divides its collected_power by a vacuum reference, and this
    function does not know that reference. So the caller divides by its OWN
    reference, or it takes the RATIO of two trials, which needs no reference.

    Args:
        result:            a TurbWaveResult with a stored patch.
        aperture_m:        the receive aperture diameter, in m.
        obscuration_ratio: the central obscuration of that aperture.
        trials:            an optional sequence of trial row indices. None
                           takes every stored trial.

    Returns:
        A float array of the power inside the aperture, one value for each
        selected trial.

    Raises:
        ValueError: the result holds no field, or the aperture is larger than
                    the stored patch.
    """
    out = []
    for _, array in _rebuilt_fields(result, aperture_m, trials):
        # The wavelength does not enter a power, so any value serves here.
        F = _patch_field(result.patch, array, 1.0)
        out.append(float(Power(_clip(F, aperture_m, obscuration_ratio))))
    return np.array(out, dtype=float)


if __name__ == '__main__':
    from ...geometry import CircularOrbit, HorizontalPath
    from ...scenario import (Channel, SpaceScenario, TerrestrialChannel,
                             TerrestrialScenario)
    from ...terminal import Terminal, Transmitter
    from ...turbulence.plane_wave_scintillation import (
        aperture_averaged_scintillation_index)
    from ...turbulence.profiles import DEFAULT_HS, default_cn2_profile
    from ..run import propagate_scenario

    t_start = time.time()
    lam = 1550e-9

    # aotools 1.0.7 reads the deprecated scipy.ndimage.interpolation namespace
    # on its FIRST call, so scipy raises a DeprecationWarning once per process.
    # Trigger it here, quietly, so it does not land in the `assert not caught`
    # block of the terrestrial vacuum check below. This warms a THIRD-PARTY
    # deprecation only; it does not hide any olb warning.
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            phase_screen(0.1, 32, 0.01, seed=0)
    except ImportError:
        pass                                   # aotools is optional.

    # ---- 6. the record NOTE names the optional stored field ----
    assert "fields" in TurbWaveResult.__doc__ and \
        "patch" in TurbWaveResult.__doc__

    # ---- 5. the documented failure modes, asserted AS failures ----
    hs = DEFAULT_HS
    orbit30 = CircularOrbit(altitude_m=600e3, elevation_deg=[30.0])
    ground = Terminal(aperture_m=0.40, wavelength_m=lam,
                      transmitter=Transmitter(waist_m=0.06))
    retro_scn = SpaceScenario(ground=ground,
                              space=Terminal(aperture_m=0.30, wavelength_m=lam),
                              direction="retro", channel=Channel())
    try:
        propagate_turbulent_scenario(retro_scn, orbit30, preset="rapid")
        raise AssertionError("retro must raise NotImplementedError")
    except NotImplementedError as exc:
        assert "retro" in str(exc), str(exc)
    try:
        folded_terrestrial()
        raise AssertionError("folded_terrestrial must raise")
    except NotImplementedError:
        pass

    # The launch aperture is 5 waists, so the clip takes almost nothing. The
    # vacuum module then takes its EXACT ABCD route (GForvard). That makes the
    # comparison below a test of this module, not of the Fresnel convolution.
    # (The Fresnel convolution loses about 1.5 percent of the power on a grid
    # of this size, and Forvard loses none. The two routes must not be mixed
    # in a 1 percent test.)
    terr_scn = TerrestrialScenario(
        near=Terminal(aperture_m=0.50, wavelength_m=lam,
                      transmitter=Transmitter(waist_m=0.05)),
        far=Terminal(aperture_m=0.20, wavelength_m=lam, detector=SMF()),
        channel=TerrestrialChannel(path_length_m=1000.0, cn2=1e-20))
    try:
        propagate_turbulent_scenario(terr_scn, HorizontalPath([1e3, 2e3]),
                                     preset="rapid")
        raise AssertionError("an array range must raise ValueError")
    except ValueError as exc:
        assert "one range" in str(exc), str(exc)
    try:
        propagate_turbulent_scenario(terr_scn, HorizontalPath(1000.0),
                                     grid=object())
        raise AssertionError("grid without plan must raise ValueError")
    except ValueError as exc:
        assert "together" in str(exc), str(exc)

    # ---- 1a. a terrestrial vacuum limit ----
    # Cn2 = 1e-20 is no turbulence. The turbulent runner must then reproduce
    # the vacuum module on the SAME grid.
    path = HorizontalPath(1000.0)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        vac = propagate_turbulent_scenario(terr_scn, path, n_trials=1, seed=1,
                                           preset="standard")
    assert not caught, [str(w.message) for w in caught]
    ref = propagate_scenario(terr_scn, path, grid=vac.grid)
    assert ref.propagator == "GForvard", ref.propagator
    ref_power = Power(ref.stages[3][1]) / Power(ref.stages[1][1])
    ref_eta = 10 ** (-ref.smf_coupling_db / 10)
    got = vac.trials[0]
    assert abs(got.collected_power / ref_power - 1.0) < 0.01, \
        (got.collected_power, ref_power)
    d_eta_db = abs(10 * np.log10(got.smf_eta / ref_eta))
    assert d_eta_db < 0.05, (got.smf_eta, ref_eta, d_eta_db)
    # The Forvard hop stays inside its sampling limit, so this IS a fair test.
    assert vac.report.step_over_limit_max <= 1.0, vac.report.step_over_limit_max

    # ---- 1b. a space vacuum limit: eta_turb is 1.0 ----
    up_scn = SpaceScenario(ground=ground,
                           space=Terminal(aperture_m=0.30, wavelength_m=lam),
                           direction="uplink", channel=Channel())
    # Cn2 = 1e-24 is no turbulence. A 40 km slab is long, so Cn2 = 1e-20 still
    # moves eta_turb by 1 percent. The rms of one realisation goes as
    # sqrt(Cn2), so four more decades give 1e-4.
    quiet_cn2 = np.full(hs.size, 1e-24)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        vac_up = propagate_turbulent_scenario(
            up_scn, orbit30, n_trials=1, seed=2, preset="rapid",
            hs=hs, cn2_profile=quiet_cn2)
    assert not caught, [str(w.message) for w in caught]
    assert abs(vac_up.trials[0].eta_turb - 1.0) < 1e-3, \
        vac_up.trials[0].eta_turb
    assert abs(vac_up.trials[0].collected_power - 1.0) < 1e-3, \
        vac_up.trials[0].collected_power
    assert vac_up.trials[0].smf_eta is None       # the ground has no detector

    # ---- 2. the seed contract ----
    down_scn = SpaceScenario(ground=ground,
                             space=Terminal(aperture_m=0.30, wavelength_m=lam),
                             direction="downlink", channel=Channel())
    with warnings.catch_warnings(record=True):
        warnings.simplefilter("always")
        r3 = propagate_turbulent_scenario(up_scn, orbit30, n_trials=3, seed=42,
                                          preset="rapid")
        r10 = propagate_turbulent_scenario(up_scn, orbit30, n_trials=5, seed=42,
                                           preset="rapid")
        r_other = propagate_turbulent_scenario(up_scn, orbit30, n_trials=3,
                                               seed=43, preset="rapid")
    for a, b in zip(r3.trials, r10.trials):
        assert a.collected_power == b.collected_power, (a, b)
        assert a.eta_turb == b.eta_turb, (a, b)
        assert a.seed_key == b.seed_key, (a, b)
    assert r3.trials[1].eta_turb != r_other.trials[1].eta_turb
    assert r3.seed_entropy == 42 and r_other.seed_entropy == 43
    # A Generator and None both resolve to one integer.
    assert isinstance(_resolve_seed(np.random.default_rng(0)), int)
    assert isinstance(_resolve_seed(None), int)

    # ---- 3. the per-trial timing ----
    with warnings.catch_warnings(record=True):
        warnings.simplefilter("always")
        timed = propagate_turbulent_scenario(up_scn, orbit30, n_trials=5,
                                             seed=7, preset="rapid")
    times = np.array([tr.wall_time_s for tr in timed.trials])
    assert np.all(times > 0.0), times
    for k, tr in enumerate(timed.trials):
        assert tr.seed_key == (7, k), tr.seed_key
        assert tr.eta_turb is not None and tr.collected_power > 0.0

    # ---- 3b. a Threader gives the SAME trials as the serial loop ----
    # The threaded run must match the serial run trial for trial: same seeds,
    # same read-only setup. Only the wall_time_s differs, so compare the
    # physics fields.
    from ..threader import Threader
    with warnings.catch_warnings(record=True):
        warnings.simplefilter("always")
        threaded = propagate_turbulent_scenario(
            up_scn, orbit30, n_trials=5, seed=7, preset="rapid",
            threader=Threader(max_workers=4))
    for a, b in zip(timed.trials, threaded.trials):
        assert a.seed_key == b.seed_key, (a.seed_key, b.seed_key)
        assert a.eta_turb == b.eta_turb, (a.eta_turb, b.eta_turb)
        assert a.collected_power == b.collected_power, \
            (a.collected_power, b.collected_power)

    # ---- 4. the downlink scintillation index ----
    # A LOOSE band, a factor of 2 each way. The tight comparison lives in the
    # step-3 examples. The rapid preset keeps the self-check short.
    n_mc = 30
    with warnings.catch_warnings(record=True):
        warnings.simplefilter("always")
        mc = propagate_turbulent_scenario(down_scn, orbit30, n_trials=n_mc,
                                          seed=2024, preset="rapid")
    power = np.array([tr.collected_power for tr in mc.trials])
    sigma2_meas = float(power.var() / power.mean() ** 2)
    cn2_prof = default_cn2_profile(down_scn.channel.site, hs)
    sigma2_theory = float(aperture_averaged_scintillation_index(
        ground.aperture_m, 30.0, lam, hs, cn2_prof))
    assert 0.5 < sigma2_meas / sigma2_theory < 2.0, \
        (sigma2_meas, sigma2_theory)

    # ---- 4b. propagate_turbulent_field repeats one trial's atmosphere ----
    # The field of trial 0 must reproduce the scalar collected_power of trial 0
    # of the mc run above, once it is clipped at the receive aperture against
    # the vacuum baseline.
    from ..field import Power as _Power_check
    from ..run import _clip as _clip_check
    F_rx, fg, fp = propagate_turbulent_field(down_scn, orbit30, seed=2024,
                                             trial=0, preset="rapid")
    assert F_rx.N == mc.grid.n, (F_rx.N, mc.grid.n)
    mask_check = super_gaussian_boundary(
        fg.n, PRESETS["rapid"].boundary_width_frac)
    F_vac_check = split_step(
        Begin(fg.size_m, lam, fg.n), fp.z_m,
        [np.zeros((fg.n, fg.n))] * fp.z_m.size, fp.z_total_m,
        boundary=mask_check)
    p_ref_check = _Power_check(_clip_check(
        F_vac_check, down_scn.ground.aperture_m,
        down_scn.ground.obscuration_ratio))
    field_power = float(_Power_check(_clip_check(
        F_rx, down_scn.ground.aperture_m,
        down_scn.ground.obscuration_ratio)) / p_ref_check)
    # The default is single precision, so the two agree to about 1e-7.
    assert abs(field_power / mc.trials[0].collected_power - 1.0) < 1e-5, \
        (field_power, mc.trials[0].collected_power)

    # ---- 4c. the opt-in "aotools" screen generator runs and agrees ----
    # The default is "olb" now (the `mc` run above). The reference "aotools"
    # path must give the SAME loose scintillation band. It draws a DIFFERENT
    # atmosphere from olb for the same seed, so the two are not bit-identical;
    # the statistics agree.
    with warnings.catch_warnings(record=True):
        warnings.simplefilter("always")
        mc_aot = propagate_turbulent_scenario(
            down_scn, orbit30, n_trials=n_mc, seed=2024, preset="rapid",
            screen_generator="aotools")
    power_aot = np.array([tr.collected_power for tr in mc_aot.trials])
    sigma2_aot = float(power_aot.var() / power_aot.mean() ** 2)
    assert 0.5 < sigma2_aot / sigma2_theory < 2.0, (sigma2_aot, sigma2_theory)
    assert not np.array_equal(power_aot, power), 'aotools must draw a new screen'
    try:
        propagate_turbulent_scenario(down_scn, orbit30, n_trials=1,
                                     preset="rapid", screen_generator="bogus")
        raise AssertionError("an unknown generator must raise ValueError")
    except ValueError as exc:
        assert "aotools" in str(exc), str(exc)

    # ---- 4d. the beamsplitter arms: ONE run feeds every detector ----
    # The default path leaves detector_etas None, so the record is unchanged.
    assert all(tr.detector_etas is None for tr in mc.trials)

    # The rx terminal of terr_scn holds an SMF. A run that ALSO asks for that
    # same SMF as an arm must give the SAME number, because both read the same
    # clipped field through the same helper.
    arms = [SMF(), Aperture(frac=0.1), Camera(pixel_pitch_m=10e-6, n_pixels=64),
            MMF(core_radius_m=25e-6, focal_length_m=0.05,
                numerical_aperture=0.2)]
    with warnings.catch_warnings(record=True):
        warnings.simplefilter("always")
        multi = propagate_turbulent_scenario(terr_scn, path, n_trials=2, seed=1,
                                             preset="standard", detectors=arms)
    for tr in multi.trials:
        assert len(tr.detector_etas) == 4, tr.detector_etas
        assert tr.detector_etas[0] == tr.smf_eta            # the same SMF arm
        assert tr.detector_etas[1] == 1.0                   # an Aperture bucket
        assert tr.detector_etas[2] is None                  # a Camera: no model
        assert 0.0 < tr.detector_etas[3] <= 1.0             # the MMF light bucket
    # The default path of the SAME seed is bit-identical: the arms change no draw.
    assert multi.trials[0].collected_power == vac.trials[0].collected_power
    assert multi.trials[0].smf_eta == vac.trials[0].smf_eta
    # An MMF arm matches a run whose rx detector IS that MMF (the same field).
    mmf_det = MMF(core_radius_m=25e-6, focal_length_m=0.05,
                  numerical_aperture=0.2)
    mmf_scn = TerrestrialScenario(
        near=terr_scn.near,
        far=Terminal(aperture_m=0.20, wavelength_m=lam, detector=mmf_det),
        channel=terr_scn.channel)
    with warnings.catch_warnings(record=True):
        warnings.simplefilter("always")
        mmf_run = propagate_turbulent_scenario(mmf_scn, path, n_trials=2, seed=1,
                                               preset="standard",
                                               detectors=[mmf_det])
    for a, b in zip(mmf_run.trials, multi.trials):
        assert a.mmf_eta == a.detector_etas[0], (a.mmf_eta, a.detector_etas)
        assert a.mmf_eta == b.detector_etas[3], (a.mmf_eta, b.detector_etas[3])
    # An unknown detector type raises.
    try:
        propagate_turbulent_scenario(terr_scn, path, n_trials=1, seed=1,
                                     preset="rapid", detectors=[object()])
        raise AssertionError("an unknown detector must raise ValueError")
    except ValueError as exc:
        assert "unknown detector" in str(exc), str(exc)

    # ---- 7a. the blocks are bit-identical to one long run ----
    # A run of 10 trials equals the two blocks 0..5 and 6..9, trial for trial.
    with warnings.catch_warnings(record=True):
        warnings.simplefilter("always")
        whole = propagate_turbulent_scenario(down_scn, orbit30, n_trials=10,
                                             seed=99, preset="rapid")
        blk_a = propagate_turbulent_scenario(down_scn, orbit30, n_trials=6,
                                             seed=99, preset="rapid",
                                             start_index=0)
        blk_b = propagate_turbulent_scenario(down_scn, orbit30, n_trials=4,
                                             seed=99, preset="rapid",
                                             start_index=6)
    parts = list(blk_a.trials) + list(blk_b.trials)
    assert len(parts) == len(whole.trials)
    for k, (a, b) in enumerate(zip(whole.trials, parts)):
        assert a.seed_key == b.seed_key == (99, k), (a.seed_key, b.seed_key, k)
        assert a.collected_power == b.collected_power, (a, b)

    # ---- 7b. the stored patch, and the SMF round trip ----
    # terr_scn holds an SMF at the far terminal. The patch must hold the whole
    # receive aperture, so its radius is a bit more than aperture_m / 2.
    rx_terr = terr_scn.rx_terminal
    with warnings.catch_warnings(record=True):
        warnings.simplefilter("always")
        kept = propagate_turbulent_scenario(
            terr_scn, path, n_trials=2, seed=1, preset="standard",
            patch_radius_m=0.6 * rx_terr.aperture_m)
    assert kept.fields.dtype == np.complex64
    assert kept.fields.shape == (2, kept.patch.indices.size)
    assert kept.patch.n == kept.grid.n
    # The patch mask IS the CircAperture mask of the same radius.
    _F_probe = Begin(kept.grid.size_m, lam, kept.grid.n)
    _keep = (_clip(_F_probe, 2 * kept.patch.radius_m, 0.0).field != 0.0)
    assert np.array_equal(np.flatnonzero(_keep.ravel()),
                          kept.patch.indices), "the patch must match CircAperture"
    # The stored run must not change one scalar.
    for a, b in zip(kept.trials, multi.trials):
        assert a.collected_power == b.collected_power, (a, b)
        assert a.smf_eta == b.smf_eta, (a, b)
    eta_back = recouple(kept, rx_terr.detector, rx_terr.aperture_m,
                        rx_terr.obscuration_ratio, lam)
    eta_run = np.array([tr.smf_eta for tr in kept.trials])
    assert np.all(np.abs(eta_back / eta_run - 1.0) < 1e-5), (eta_back, eta_run)
    # The collected power: take the RATIO of two trials. That is the lazier
    # check, because a ratio needs no vacuum reference.
    pw = recollect(kept, rx_terr.aperture_m, rx_terr.obscuration_ratio)
    run_pw = np.array([tr.collected_power for tr in kept.trials])
    assert abs((pw[0] / pw[1]) / (run_pw[0] / run_pw[1]) - 1.0) < 1e-5, \
        (pw, run_pw)

    # ---- 7c. the MMF round trip ----
    rx_mmf = mmf_scn.rx_terminal
    with warnings.catch_warnings(record=True):
        warnings.simplefilter("always")
        kept_mmf = propagate_turbulent_scenario(
            mmf_scn, path, n_trials=2, seed=1, preset="standard",
            patch_radius_m=0.6 * rx_mmf.aperture_m)
    mmf_back = recouple(kept_mmf, rx_mmf.detector, rx_mmf.aperture_m,
                        rx_mmf.obscuration_ratio, lam)
    mmf_run = np.array([tr.mmf_eta for tr in kept_mmf.trials])
    assert np.all(np.abs(mmf_back / mmf_run - 1.0) < 1e-5), (mmf_back, mmf_run)

    # ---- 7d. no patch means no change at all ----
    assert multi.fields is None and multi.patch is None
    try:
        recouple(multi, rx_terr.detector, rx_terr.aperture_m,
                 rx_terr.obscuration_ratio, lam)
        raise AssertionError("a result with no field must raise ValueError")
    except ValueError as exc:
        assert "no field" in str(exc), str(exc)
    # An aperture larger than the patch raises.
    try:
        recouple(kept, rx_terr.detector, 4 * kept.patch.radius_m, 0.0, lam)
        raise AssertionError("an aperture past the patch must raise")
    except ValueError as exc:
        assert "patch" in str(exc), str(exc)
    # A patch larger than the grid raises.
    try:
        propagate_turbulent_scenario(terr_scn, path, n_trials=1, seed=1,
                                     preset="rapid", patch_radius_m=1e6)
        raise AssertionError("a patch past the grid must raise")
    except ValueError as exc:
        assert "grid side" in str(exc), str(exc)

    # ---- 8. the precision switch ----
    # The default is "single" (owner decision 2026-09-05), so every result
    # above ran in single precision. A "double" run of the same seed gives the
    # SAME physics to about six digits, and its snapshot is not bit-identical,
    # because the arithmetic differs.
    with warnings.catch_warnings(record=True):
        warnings.simplefilter("always")
        double = propagate_turbulent_scenario(
            down_scn, orbit30, n_trials=6, seed=99, preset="rapid",
            precision="double")
    d_power = np.array([abs(a.collected_power / b.collected_power - 1.0)
                        for a, b in zip(double.trials, whole.trials)])
    assert d_power.max() < 1e-3, d_power
    assert d_power.max() > 0.0, "the precision must change the last digits"
    # The field entry point takes the same switch, and its default is single.
    F_s, _, _ = propagate_turbulent_field(down_scn, orbit30, seed=2024,
                                          trial=0, preset="rapid")
    assert F_s.field.dtype == np.complex64, F_s.field.dtype
    F_d, _, _ = propagate_turbulent_field(down_scn, orbit30, seed=2024,
                                          trial=0, preset="rapid",
                                          precision="double")
    assert F_d.field.dtype == np.complex128, F_d.field.dtype
    # An unknown name raises.
    for call in (lambda: propagate_turbulent_scenario(
                     down_scn, orbit30, n_trials=1, preset="rapid",
                     precision="half"),
                 lambda: propagate_turbulent_field(
                     down_scn, orbit30, preset="rapid", precision="half")):
        try:
            call()
            raise AssertionError("an unknown precision must raise ValueError")
        except ValueError as exc:
            assert "single" in str(exc), str(exc)

    # ---- the printed tables ----
    print("terrestrial vacuum limit, 1 km, Cn2 = 1e-20, standard preset:")
    print(f"  grid                    {vac.grid.n:11d} px, "
          f"{vac.grid.size_m:.3f} m")
    print(f"  screens                 {vac.plan.z_m.size:11d}")
    print(f"  collected, turbulent    {got.collected_power:11.6f}")
    print(f"  collected, vacuum       {ref_power:11.6f}")
    print(f"  SMF eta, turbulent      {got.smf_eta:11.6f}")
    print(f"  SMF eta, vacuum         {ref_eta:11.6f}")
    print(f"  difference              {d_eta_db:11.4f} dB")
    print("")
    print("space vacuum limit, uplink at 30 deg, Cn2 = 1e-24:")
    print(f"  eta_turb                {vac_up.trials[0].eta_turb:11.6f}")
    print("")
    print(f"per-trial timing, 5 uplink trials, {timed.preset} preset, "
          f"{timed.grid.n} px, {timed.plan.z_m.size} screens:")
    print(f"  mean                    {times.mean():11.3f} s")
    print(f"  min                     {times.min():11.3f} s")
    print(f"  max                     {times.max():11.3f} s")
    print(f"  total                   {times.sum():11.3f} s")
    print("")
    print(f"downlink scintillation, {n_mc} trials at 30 deg, "
          f"D = {ground.aperture_m} m:")
    print(f"  grid                    {mc.grid.n:11d} px, "
          f"{mc.grid.size_m:.3f} m")
    print(f"  screens                 {mc.plan.z_m.size:11d}")
    print(f"  r0 total                {mc.plan.r0_total_m * 1e2:11.3f} cm")
    print(f"  mean collected power    {power.mean():11.5f}")
    print(f"  sigma2_I, wave optics   {sigma2_meas:11.5f}")
    print(f"  sigma2_I, analytic      {sigma2_theory:11.5f}")
    print(f"  ratio                   {sigma2_meas / sigma2_theory:11.3f}")
    print(f"  sigma2_I, aotools gen   {sigma2_aot:11.5f}  "
          f"(mean power {power_aot.mean():.5f})")
    print("")
    print("blocks and the stored field:")
    print(f"  6 + 4 equals 10         {'yes':>11s}")
    print(f"  patch radius            {kept.patch.radius_m * 1e2:11.2f} cm")
    print(f"  patch pixels            {kept.patch.indices.size:11d}")
    print(f"  stored bytes per trial  {kept.fields[0].nbytes / 1024:11.1f} kB")
    print(f"  SMF eta, in run         {eta_run[0]:11.6f}")
    print(f"  SMF eta, recoupled      {eta_back[0]:11.6f}")
    print(f"  MMF eta, in run         {mmf_run[0]:11.6f}")
    print(f"  MMF eta, recoupled      {mmf_back[0]:11.6f}")
    print("")
    print("single against double precision, 6 downlink trials:")
    print(f"  max relative difference {d_power.max():11.2e}")
    print("")
    print(f"(elapsed {time.time() - t_start:.1f} s)")
    print("self-check passed")
