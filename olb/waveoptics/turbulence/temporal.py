"""The frozen-flow time axis of the screen stack: oversized STRIP screens.

The split-step layer gives independent SNAPSHOTS. Each seed gives one new
atmosphere. A fade DURATION and a fade RATE need a time axis, so they need a
screen stack that MOVES. This module makes that stack.

THE ROUTE. Each layer of the plan gets ONE oversized RECTANGULAR screen, a
STRIP, at the pitch of the propagation grid (`olb.waveoptics.turbulence.
screens.ScreenFactory(nx=...)`). The strip is long enough to hold the whole
travel of that layer for the whole record. A FRAME is then a crop of the strip
at an integer pixel offset, and the offset of frame k comes from the ABSOLUTE
position v*k*dt, never from a sum of steps. That is the frozen-flow hypothesis
of Taylor, DOI 10.1098/rspa.1938.0032: the pattern moves with the wind and it
does not change shape.

THE ROUTE THAT WAS REJECTED: the row-by-row EXTRUSION of Assemat and Wilson,
DOI 10.1364/OE.14.000988, which `aotools` implements and which the old stub of
this module planned.

  1. ACCURACY. The extruded screen OVER-CORRELATES its own axis when the frame
     side is 0.1 to 0.35 of the outer scale, and that band IS the production
     regime (a 512 px grid at L0 = 25 m). No `n_columns` setting fixes it. See
     `validation/screens/FINDINGS.md`, Q5 point 4, Q6 and Q8.
  2. SPEED (measured 2026-09-13). The aotools `add_row` costs 0.93 ms per row
     at 512 px and 2.38 ms at 1024 px, and the extrusion must make EVERY
     intermediate row of the travel. One strip column costs 12 to 26 us, and a
     frame is a slice of an open memory map. So the strip is about 75 to 90
     times cheaper for each pixel of travel.
  3. MEMORY is the ONE place the extrusion wins: its state is a fixed ~24 MB,
     and a 2 s record of the hero downlink holds about 0.5 GB of strips. So the
     extrusion wins only under about 1 GiB of spare memory. That case is
     RECORDED here, and it is not built.

THE VELOCITY OF A LAYER has two parts, and it is a 2-D vector:

  1. THE WIND. The Bufton profile V(h) of Andrews and Phillips,
     DOI 10.1117/3.626196, Ch. 12, Eqs. (2) and (3), printed p. 481. This
     module calls it with ws = 0 (`olb.turbulence.profiles.v_wind`), so the
     slew is NOT counted two times, and it points at `wind_dir_deg`.
  2. THE SLEW. A tracked satellite drags the line of sight across each layer.
     `olb.geometry.CircularOrbit.slew_deg_s` is v*sin(el)/h_sat, which is the
     COEFFICIENT of the `ws*h` term of the Bufton profile, so the apparent
     speed at a layer is that rate times the ALTITUDE h of the layer, NOT times
     its slant distance: omega_true * z_slant = (v sin(el)/h_sat) * h. It goes
     along x, the long axis of the strip.

THE SEAM. A Fourier screen is periodic, so the far end of a strip joins its
near end. The strip therefore carries a PAD of `pad_outer_scales` outer scales
past the travel. The von Karman phase covariance at 2 L0 is about 4e-6 of the
variance (Assemat and Wilson, DOI 10.1364/OE.14.000988, Eq. (5)), so the seam
is not measurable in the record.

THE SEED RULE STAYS IN THE RUNNER. `build_strips` takes the integer seed of
each layer as an argument. `olb.waveoptics.turbulence.run._screen_seed` is the
ONE owner of that rule, and this module does not import it, so the two modules
have no cycle.

Sources:
- Taylor, The spectrum of turbulence, DOI 10.1098/rspa.1938.0032. The frozen
  flow hypothesis.
- Andrews and Phillips, Laser Beam Propagation through Random Media, 2nd ed.,
  DOI 10.1117/3.626196, Ch. 12, Eqs. (2) and (3), printed p. 481. The Bufton
  wind profile with the slew term.
- Schmidt, Numerical Simulation of Optical Wave Propagation with Examples in
  MATLAB, DOI 10.1117/3.866274, Ch. 9, Eqs. (9.78) to (9.81), printed pp. 166
  to 169. The Fourier-series screen and its subharmonics.
- Assemat and Wilson, DOI 10.1364/OE.14.000988. The rejected extrusion, and
  Eq. (5), the closed-form von Karman covariance that bounds the seam.
- Greenwood, Bandwidth specification for adaptive optics systems,
  DOI 10.1364/JOSA.67.000390. The Greenwood frequency, the time scale that the
  sample step dt must resolve.
"""

import os
from dataclasses import dataclass, fields

import numpy as np

from olb.turbulence.profiles import v_wind

from ..propagators import set_fft_backend
from .screens import ScreenFactory

# The column count of a strip rounds UP to this multiple, as the oversize
# screen of the point-ahead run does (run.SCREEN_DRAW_ALIGN). An FFT of a
# length with small prime factors is much faster than an FFT of a prime length.
STRIP_ALIGN = 32


@dataclass(frozen=True)
class TemporalSpec:
    """What one temporal RECORD asks for.

    A record is a run of `n_frames` frames at the fixed step `dt_s`. Frame k
    of the record is trial k of the runner, so the frames ARE trials and the
    campaign addresses them by the trial index.

    Attributes:
        dt_s:             the time step, in s. It must resolve the Greenwood
                          time 1/f_G (Greenwood, DOI 10.1364/JOSA.67.000390).
        n_frames:         the number of frames of the record.
        strip_dir:        where the strips are kept. It is NOT part of `key`.
        record:           the record index. Two records of one campaign have
                          different strips and the same statistics.
        wind_ground_m_s:  the ground wind Vg of the Bufton profile, in m/s.
        wind_dir_deg:     the direction of the wind, in deg, from the x axis.
                          The slew always goes along x.
        pad_outer_scales: the seam pad, in outer scales. See the module
                          docstring.
        slew_rad_s:       the slew rate, in rad/s, or None to read the
                          geometry. See the module docstring for the altitude
                          rule.
    """

    dt_s: float
    n_frames: int
    strip_dir: str
    record: int = 0
    wind_ground_m_s: float = 10.0
    wind_dir_deg: float = 0.0
    slew_rad_s: float = None
    pad_outer_scales: float = 2.0

    def key(self):
        """Give a stable string of the spec, WITHOUT `strip_dir`.

        The campaign fingerprint reads this. The strips are a deletable cache
        that a seed rebuilds, so WHERE they sit must not name the campaign.

        Returns:
            A string.
        """
        parts = [f"{f.name}={getattr(self, f.name)!r}"
                 for f in fields(self) if f.name != "strip_dir"]
        return ",".join(parts)


@dataclass(frozen=True)
class StripPlan:
    """The size and the drift of the strip of each layer.

    Attributes:
        v_x_m_s:  the drift speed along x of each layer, in m/s. It holds the
                  slew and the x part of the wind.
        v_y_m_s:  the drift speed along y of each layer, in m/s.
        shape:    the (ny, nx) pixel shape of the strip of each layer, as a
                  tuple of pairs. Each layer gets its OWN length, because the
                  top layer of a LEO pass runs about ten times faster than the
                  ground layer.
        dt_s:     the time step, in s.
        n_frames: the number of frames.
        dx:       the pixel pitch, in m.
        n:        the pixel count of one side of the propagation grid.
        pad:      the seam pad, in pixels.
    """

    v_x_m_s: np.ndarray
    v_y_m_s: np.ndarray
    shape: tuple
    dt_s: float
    n_frames: int
    dx: float
    n: int
    pad: int


def strip_plan(plan, grid, spec, geometry, L0_m):
    """Size the strip of every layer of a screen plan.

    The drift of a layer is the vector sum of the Bufton wind at its altitude
    and the apparent slew of the tracked satellite (see the module docstring
    for the two sources and for the altitude rule of the slew).

    Args:
        plan:     the ScreenPlan of the run.
        grid:     the GridSpec of the run. It gives `n` and `pixel_m`.
        spec:     the TemporalSpec of the record.
        geometry: the link geometry. It gives `elevation_deg` and, when the
                  spec gives no `slew_rad_s`, `slew_deg_s`.
        L0_m:     the outer scale, in m. It sets the seam pad, so it must be
                  finite.

    Returns:
        A StripPlan.

    Raises:
        ValueError: the outer scale is not finite.
    """
    if not np.isfinite(L0_m) or float(L0_m) <= 0.0:
        raise ValueError("strip_plan: the seam pad needs a FINITE outer scale "
                         f"L0_m, not {L0_m!r}.")
    # The distance of each screen from the GROUND plane. A space plan counts
    # z_m from the TOP of the slab. This is the rule of
    # olb.waveoptics.turbulence.run._ground_distance, which owns it.
    z = np.asarray(plan.z_m, dtype=float)
    z_g = (float(plan.z_total_m) - z if plan.direction == "down" else z)

    # The sizer convention: one line of sight, at the lowest elevation.
    elevation = float(np.min(np.asarray(geometry.elevation_deg, dtype=float)))
    h = z_g * np.sin(np.deg2rad(elevation))

    if spec.slew_rad_s is None:
        omega = np.deg2rad(float(np.ravel(
            np.asarray(geometry.slew_deg_s, dtype=float))[0]))
    else:
        omega = float(spec.slew_rad_s)
    v_slew = omega * h                       # See the module docstring.
    # ws = 0: the slew is already in v_slew, so the Bufton slew term must not
    # add it a second time.
    v_buf = np.asarray(v_wind(h, ws=0.0, Vg=float(spec.wind_ground_m_s)),
                       dtype=float)
    direction = np.deg2rad(float(spec.wind_dir_deg))
    v_x = v_slew + v_buf * np.cos(direction)
    v_y = v_buf * np.sin(direction)

    dx = float(grid.pixel_m)
    n = int(grid.n)
    total_s = (int(spec.n_frames) - 1) * float(spec.dt_s)
    pad = int(np.ceil(float(spec.pad_outer_scales) * float(L0_m) / dx))

    shape = []
    for vx, vy in zip(v_x, v_y):
        nx = n + int(np.ceil(abs(vx) * total_s / dx)) + pad
        nx = STRIP_ALIGN * int(np.ceil(nx / STRIP_ALIGN))
        ny = n + int(np.ceil(abs(vy) * total_s / dx))
        if vy != 0.0:
            ny += pad
        shape.append((ny, nx))

    return StripPlan(v_x_m_s=v_x, v_y_m_s=v_y, shape=tuple(shape),
                     dt_s=float(spec.dt_s), n_frames=int(spec.n_frames),
                     dx=dx, n=n, pad=pad)


def frame_offsets(sp, k):
    """Give the (row, column) pixel offset of the window of frame k.

    The offset comes from the ABSOLUTE position v*k*dt, and it is rounded to
    the nearest pixel. It is never a sum of steps, so the rounding error does
    not build up: it stays under half a pixel for every frame.

    A NEGATIVE speed starts at the FAR end of the strip and it walks back, so
    the window stays inside the strip on the two axes.

    Args:
        sp: the StripPlan.
        k:  the frame index, from 0.

    Returns:
        A list of (oy, ox) integer pairs, one for each layer.
    """
    t = float(k) * sp.dt_s
    out = []
    for (ny, nx), vx, vy in zip(sp.shape, sp.v_x_m_s, sp.v_y_m_s):
        base_x = 0 if vx >= 0.0 else nx - sp.n
        base_y = 0 if vy >= 0.0 else ny - sp.n
        ox = base_x + int(np.rint(vx * t / sp.dx))
        oy = base_y + int(np.rint(vy * t / sp.dx))
        out.append((oy, ox))
    return out


def strip_paths(spec, n_screens):
    """Give the file path of the strip of each layer of one record.

    Args:
        spec:      the TemporalSpec. It gives `strip_dir` and `record`.
        n_screens: the number of layers of the plan.

    Returns:
        A list of paths.
    """
    return [os.path.join(spec.strip_dir,
                         f"r{int(spec.record):04d}_s{j:02d}.npy")
            for j in range(int(n_screens))]


def build_strips(sp, r0_m, paths, L0_m, subharmonics, seeds,
                 dtype=np.float32):
    """Make the strip of every layer, and write each one to its own file.

    THE CALL IS IDEMPOTENT. A layer whose file is already there is SKIPPED, so
    a resumed record makes no strip a second time. Each file is written to a
    temporary name and then moved with `os.replace`, so a killed run leaves no
    half file.

    ONE LAYER AT A TIME. A strip of the top layer of a LEO pass is the largest
    array of the whole run, so the peak memory holds ONE strip and its
    transform, not the stack.

    THE HOST FFT. The build forces the "scipy" FFT backend and it restores the
    previous backend. A strip does not fit the memory of a device as well as
    the field does, and the build runs one time for each record.

    Args:
        sp:           the StripPlan.
        r0_m:         the Fried parameter of each layer, in m.
        paths:        the file path of each layer. Use `strip_paths`.
        L0_m:         the outer scale, in m.
        subharmonics: True adds the subharmonic levels.
        seeds:        the integer seed of each layer. The runner owns the seed
                      rule (see the module docstring).
        dtype:        the stored floating type. float32 halves the file.

    Returns:
        The list of paths.

    Raises:
        ValueError: the lengths of the arguments do not agree.
    """
    r0 = np.asarray(r0_m, dtype=float).ravel()
    if not (len(paths) == len(seeds) == r0.size == len(sp.shape)):
        raise ValueError(
            f"build_strips: the plan holds {len(sp.shape)} layers, and the "
            f"call gives {len(paths)} paths, {len(seeds)} seeds and "
            f"{r0.size} r0 values.")
    if paths:
        os.makedirs(os.path.dirname(paths[0]) or ".", exist_ok=True)

    previous = set_fft_backend("scipy")
    try:
        for j, path in enumerate(paths):
            if os.path.exists(path):
                continue
            ny, nx = sp.shape[j]
            factory = ScreenFactory(ny, sp.dx, L0_m=float(L0_m), nx=nx,
                                    subharmonics=bool(subharmonics),
                                    dtype=dtype)
            strip = factory.make(float(r0[j]),
                                 np.random.default_rng(int(seeds[j])))
            tmp = path + ".tmp"
            # np.save adds ".npy" to a NAME that does not end with it, and the
            # temporary name does not. A file OBJECT keeps the name as it is.
            with open(tmp, "wb") as handle:
                np.save(handle, np.asarray(strip, dtype=dtype))
            os.replace(tmp, path)
            del strip, factory
    finally:
        set_fft_backend(previous)
    return list(paths)


def open_strips(paths):
    """Open the strip of every layer as a read-only memory map.

    A memory map shares ONE page cache between the worker processes of a pool,
    so a 12-worker campaign holds one copy of the strips, not twelve.

    Args:
        paths: the file path of each layer.

    Returns:
        A list of read-only numpy memory maps.
    """
    return [np.load(p, mmap_mode="r") for p in paths]


def frame_stack(strips, sp, k):
    """Give the phase screen of every layer at frame k.

    Each screen is a VIEW of the memory map, so the call copies nothing. A
    device run makes the copy itself, when it uploads the view.

    Args:
        strips: the open strips. Use `open_strips`.
        sp:     the StripPlan.
        k:      the frame index, from 0.

    Yields:
        One n by n array of the phase, in radians, for each layer, in the
        order of the plan.
    """
    n = sp.n
    for strip, (oy, ox) in zip(strips, frame_offsets(sp, k)):
        yield strip[oy:oy + n, ox:ox + n]


if __name__ == '__main__':
    import shutil
    import tempfile
    import types

    # A small fake plan and a small fake geometry. The module reads four
    # attributes of a plan and two of a geometry, so a namespace is enough and
    # the self-check needs no scenario.
    plan = types.SimpleNamespace(
        z_m=np.array([18000.0, 10000.0, 0.0]),   # from the TOP of the slab.
        z_total_m=20000.0,
        direction="down",
        r0_m=np.array([0.9, 0.5, 0.16]))
    grid = types.SimpleNamespace(n=64, pixel_m=0.01)
    geometry = types.SimpleNamespace(elevation_deg=90.0, slew_deg_s=0.5)
    L0 = 25.0

    # ---- 1. the shapes, against a hand-computed case ----
    spec = TemporalSpec(dt_s=1e-3, n_frames=11, strip_dir="unused",
                        wind_ground_m_s=10.0, wind_dir_deg=0.0)
    sp = strip_plan(plan, grid, spec, geometry, L0)
    z_g = plan.z_total_m - plan.z_m            # 2000, 10000, 20000 m.
    h_hand = z_g * 1.0                         # 90 deg: h = z_g.
    v_hand = (np.deg2rad(0.5) * h_hand
              + v_wind(h_hand, ws=0.0, Vg=10.0))
    assert np.allclose(sp.v_x_m_s, v_hand), (sp.v_x_m_s, v_hand)
    assert np.allclose(sp.v_y_m_s, 0.0), sp.v_y_m_s
    pad_hand = int(np.ceil(2.0 * L0 / 0.01))   # 5000 px.
    assert sp.pad == pad_hand, (sp.pad, pad_hand)
    for (ny, nx), vx in zip(sp.shape, sp.v_x_m_s):
        want = 64 + int(np.ceil(abs(vx) * 0.010 / 0.01)) + pad_hand
        assert ny == 64, ny                    # v_y is zero: no pad, no room.
        assert nx == 32 * int(np.ceil(want / 32)), (nx, want)
        assert nx >= want
    # z_m counts from the TOP of the slab, so layer 0 is the GROUND layer and
    # the last layer is the top one. The top layer is the fastest, so it holds
    # the longest strip.
    assert sp.shape[-1][1] > sp.shape[0][1], sp.shape

    # ---- 2. the offsets ----
    assert frame_offsets(sp, 0) == [(0, 0)] * 3, frame_offsets(sp, 0)
    ox_last = [ox for _, ox in frame_offsets(sp, spec.n_frames - 1)]
    for j in range(1, spec.n_frames):
        now = [ox for _, ox in frame_offsets(sp, j)]
        before = [ox for _, ox in frame_offsets(sp, j - 1)]
        assert all(a >= b for a, b in zip(now, before)), (now, before)
    for (ny, nx), ox in zip(sp.shape, ox_last):
        assert 0 <= ox <= nx - sp.n, (ox, nx)

    # ---- 3. a negative velocity keeps the window inside the strip ----
    spec_back = TemporalSpec(dt_s=1e-3, n_frames=11, strip_dir="unused",
                             wind_dir_deg=180.0, slew_rad_s=0.0)
    sp_back = strip_plan(plan, grid, spec_back, geometry, L0)
    assert np.all(sp_back.v_x_m_s < 0.0), sp_back.v_x_m_s
    for k in (0, 5, spec_back.n_frames - 1):
        for ((ny, nx), (oy, ox)) in zip(sp_back.shape,
                                        frame_offsets(sp_back, k)):
            assert 0 <= ox <= nx - sp_back.n, (k, ox, nx)
            assert 0 <= oy <= ny - sp_back.n, (k, oy, ny)

    # ---- 4. a diagonal wind pads the y axis too ----
    spec_diag = TemporalSpec(dt_s=1e-3, n_frames=11, strip_dir="unused",
                             wind_dir_deg=90.0, slew_rad_s=0.0)
    sp_diag = strip_plan(plan, grid, spec_diag, geometry, L0)
    assert np.all(sp_diag.v_y_m_s > 0.0), sp_diag.v_y_m_s
    assert all(ny > 64 + pad_hand - 1 for ny, _ in sp_diag.shape), sp_diag.shape

    # ---- 5. build, open and crop, in a temporary directory ----
    tmpdir = tempfile.mkdtemp(prefix="olb_strips_")
    try:
        # A small record: a short strip, so the check is quick.
        spec_s = TemporalSpec(dt_s=1e-3, n_frames=4, strip_dir=tmpdir,
                              record=2, pad_outer_scales=0.02,
                              slew_rad_s=0.0, wind_ground_m_s=10.0)
        grid_s = types.SimpleNamespace(n=32, pixel_m=0.02)
        sp_s = strip_plan(plan, grid_s, spec_s, geometry, L0)
        paths = strip_paths(spec_s, 3)
        assert paths[0].endswith("r0002_s00.npy"), paths[0]
        seeds = [11, 22, 33]
        build_strips(sp_s, plan.r0_m, paths, L0, True, seeds)
        assert all(os.path.exists(p) for p in paths)

        # 5a. idempotent: a second call writes no file again.
        before = [os.stat(p).st_mtime_ns for p in paths]
        build_strips(sp_s, plan.r0_m, paths, L0, True, seeds)
        assert [os.stat(p).st_mtime_ns for p in paths] == before

        # 5b. the seed reproduces the strip.
        first = np.load(paths[0])
        again = ScreenFactory(sp_s.shape[0][0], sp_s.dx, L0_m=L0,
                              nx=sp_s.shape[0][1],
                              dtype=np.float32).make(
                                  float(plan.r0_m[0]),
                                  np.random.default_rng(11))
        assert np.array_equal(first, again), 'the seed must rebuild the strip'
        assert first.dtype == np.float32, first.dtype

        # 5c. the frames are n by n, and two frames differ.
        strips = open_strips(paths)
        f0 = list(frame_stack(strips, sp_s, 0))
        f3 = list(frame_stack(strips, sp_s, spec_s.n_frames - 1))
        assert all(f.shape == (32, 32) for f in f0), [f.shape for f in f0]
        assert all(f.shape == (32, 32) for f in f3)
        assert not np.array_equal(f0[0], f3[0]), 'the atmosphere must move'
        # The GROUND layer (layer 0) is the slowest, so it moves the least.
        moved = [float(np.mean((a - b) ** 2)) for a, b in zip(f0, f3)]
        assert moved[0] < moved[-1], moved
        del strips, f0, f3
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)

    # ---- 6. key() omits strip_dir ----
    k1 = TemporalSpec(dt_s=1e-3, n_frames=4, strip_dir="/a").key()
    k2 = TemporalSpec(dt_s=1e-3, n_frames=4, strip_dir="/b").key()
    k3 = TemporalSpec(dt_s=2e-3, n_frames=4, strip_dir="/a").key()
    assert k1 == k2, (k1, k2)
    assert k1 != k3, (k1, k3)
    assert "strip_dir" not in k1, k1
    assert "dt_s" in k1 and "record" in k1, k1

    # ---- 7. a non-finite outer scale raises ----
    try:
        strip_plan(plan, grid, spec, geometry, np.inf)
        raise AssertionError('an infinite L0 must raise')
    except ValueError as exc:
        assert 'FINITE' in str(exc), str(exc)

    print(f"layers                      {len(sp.shape)}")
    print(f"elevation                   {geometry.elevation_deg:9.1f} deg")
    print("layer      h [m]   v_x [m/s]      strip [px]")
    for j, ((ny, nx), vx) in enumerate(zip(sp.shape, sp.v_x_m_s)):
        print(f"  {j:2d}  {h_hand[j]:10.0f}  {vx:10.2f}   ({ny}, {nx})")
    print(f"seam pad                    {sp.pad} px "
          f"({sp.pad * sp.dx:.1f} m = {spec.pad_outer_scales:.1f} L0)")
    print(f"record                      {sp.n_frames} frames at "
          f"{sp.dt_s * 1e3:.3f} ms")
    print(f"key()                       {k1}")
    print("self-check passed")
