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

THE ROTATED STRIP (an OPT-IN, 2026-09-13, backlog 2-P1b item 10). A CROSSWIND
turns the walk of a layer away from the x axis. The axis-aligned box that holds
that walk grows its SHORT axis, and the 30 deg hero at `wind_dir_deg = 90` asks
for a 11897 x 43744 box (520 Mpx) for one jet-level layer. `TemporalSpec.
rotated=True` holds ONE THIN strip per layer whose long axis runs along the
RESULTANT velocity of that layer. The walk is then purely along the strip x, as
it is for a wind along track. Each frame takes a PADDED crop of the strip and
turns it back by the layer angle with `rotate_fourier`, then it cuts the
central n by n. The default is False, so the along-track route does not move
by one bit.

THE ROTATION IS EXACT. A Fourier screen is band-limited to its own grid, so the
three-shear rotation of Unser, Thevenaz and Yaroslavsky,
DOI 10.1109/83.469963, is an exact resampling of it: a rotation is an x-shear,
a y-shear and an x-shear, and each shear translates one line by a sub-pixel
amount, which is a linear phase ramp on the 1-D transform of that line (the
Fourier shift theorem, Schmidt, DOI 10.1117/3.866274, Ch. 2). A REAL-SPACE
bilinear or cubic rotation is NOT the tool: it smooths the Fresnel-scale
structure that builds the scintillation. Two errors stay:
  1. THE CORNERS of the frequency square leave the square under a rotation.
  2. THE EDGES of a crop are not periodic, and a shear wraps. The sinc tail
     falls as 1/distance, so a margin holds the error inside the cut.

THE TAPER (2026-09-13, the measured cure of error 2). Each shear reads a line
as PERIODIC, so the STEP between the two ends of the line rings inward. A
1-D raised-cosine taper on the ends of every line, applied BEFORE EACH of the
three shear transforms, removes that step. A single 2-D window does NOT work,
because the two later shears undo it. The taper is flat over the axis-aligned
footprint n(|cos t| + |sin t|) of the rotated frame, so the kept frame reads
only the flat top and the taper does not touch the answer. The measurement is
`validation/temporal_screens/rotation_taper.py` (variant V4).

THE CROP SIDE OF A LAYER is therefore

    m = n (|cos t| + |sin t|)  +  2 * rot_margin * n,

the geometric footprint plus the ROLL-OFF of the taper on each side.
`TemporalSpec.rot_margin` is that ROLL-OFF, not the whole margin, and its
default 0.07 is the measured value (the gate is an outer-16 px ring rms under
1e-2 rad and |D(r) ratio - 1| under 1 percent at 1, 4 and 16 px). Each layer
gets its OWN m, because a layer that walks along x asks for much less than a
layer that walks at 45 deg.

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
- Unser, Thevenaz and Yaroslavsky, Convolution-based interpolation for fast,
  high-quality rotation of images, IEEE Trans. Image Process. 4, 1371 (1995),
  DOI 10.1109/83.469963. The three-shear rotation of the rotated strip.
"""

import os
import warnings
from dataclasses import dataclass, fields

import numpy as np

from olb.turbulence.profiles import v_wind

from ..propagators import set_fft_backend, xp
from .screens import ScreenFactory

# The column count of a strip rounds UP to this multiple, as the oversize
# screen of the point-ahead run does (run.SCREEN_DRAW_ALIGN). An FFT of a
# length with small prime factors is much faster than an FFT of a prime length.
STRIP_ALIGN = 32


def _shear(a, shift_px, axis):
    """Translate every line of `a` along `axis` by its own sub-pixel amount.

    The shift is a linear phase ramp on the 1-D transform of the line (the
    Fourier shift theorem; Schmidt, DOI 10.1117/3.866274, Ch. 2). A line of a
    band-limited periodic array is resampled EXACTLY this way.

    Args:
        a:        a real 2-D array, on the array module of the FFT backend.
        shift_px: the shift of each line, in pixels. It has one value for each
                  line, and it broadcasts against `a`.
        axis:     1 shifts the ROWS along x, 0 shifts the COLUMNS along y.

    Returns:
        A real array of the shape of `a`.
    """
    xpm = xp()
    m = a.shape[axis]
    f = xpm.fft.fftfreq(m)
    if axis == 1:
        ramp = xpm.exp(-2j * xpm.pi * f[None, :] * shift_px[:, None])
    else:
        ramp = xpm.exp(-2j * xpm.pi * f[:, None] * shift_px[None, :])
    return xpm.real(xpm.fft.ifft(xpm.fft.fft(a, axis=axis) * ramp, axis=axis))


def _taper_1d(m, flat_px, dtype):
    """Give the 1-D raised-cosine profile of a line of m px.

    It is 1 where |c| <= flat_px / 2, with c the coordinate about the rotation
    centre `m // 2`, and it falls to 0 at the two ends of the line. `flat_px`
    is the footprint of the rotated square, so the kept frame never leaves the
    flat top. See the module docstring for the measurement.

    Args:
        m:       the line length, in px.
        flat_px: the width of the flat top, in px.
        dtype:   the floating type of the profile.

    Returns:
        A 1-D array of m values, on the array module of the FFT backend.

    Raises:
        ValueError: the flat top fills the line, so no roll-off fits.
    """
    xpm = xp()
    half = 0.5 * float(flat_px)
    edge = 0.5 * float(m)
    if half >= edge:
        raise ValueError(
            f"_taper_1d: the rotation footprint {flat_px:.1f} px fills the "
            f"{m} px crop, so no roll-off fits. Raise TemporalSpec."
            "rot_margin.")
    c = xpm.abs(xpm.arange(m) - m // 2)
    t = xpm.clip((c - half) / (edge - half), 0.0, 1.0)
    return (0.5 * (1.0 + xpm.cos(xpm.pi * t))).astype(dtype)


def _quarter(a):
    """Give the EXACT quarter turn `out(p) = a(R(pi/2) p)` of a square array.

    The centre of the shears of `rotate_fourier` is the pixel `m // 2`, and
    `numpy.rot90` turns about the middle of the array, which is half a pixel
    away for an even side. So the map is written out: with the coordinate
    c_i = i - m//2, the row index carries y and the column index carries x,
    and R(pi/2) sends (x, y) to (-y, x).
    """
    xpm = xp()
    return xpm.roll(xpm.rot90(a, 1), 1, axis=0)


def rotate_fourier(patch, theta_rad, taper=None):
    """Rotate a square patch by the EXACT Fourier three-shear route.

    The output is `out(p) = patch(R(theta) p)` about the centre of the patch,
    with R the rotation matrix. The rotation is the product of three shears

        R(a) = Sx(-tan(a/2)) Sy(sin a) Sx(-tan(a/2))

    of Unser, Thevenaz and Yaroslavsky, DOI 10.1109/83.469963, and each shear
    is one line translation, which `_shear` does with a phase ramp. So the
    call is exact for a band-limited periodic patch, and it does NOT smooth
    the Fresnel-scale structure that a real-space bilinear or cubic
    interpolation smooths.

    THE ANGLE IS REDUCED to [-45, 45] degrees by whole quarter turns, because
    tan(a/2) grows without bound near a quarter turn. A quarter turn is exact
    and free (`numpy.rot90`).

    TWO ERRORS STAY. The modes in the CORNERS of the frequency square leave
    the square under a rotation. The EDGES of a patch are not periodic, and a
    shear wraps, so the caller must cut a central window out of a patch that
    carries a margin (see `TemporalSpec.rot_margin`). `taper` removes most of
    the second error; see the module docstring.

    Args:
        patch:     a real, square 2-D array. It may be a host array or a
                   device array.
        theta_rad: the rotation angle, in rad.
        taper:     the width of the FLAT top of the 1-D end taper, in px, or
                   None (the default) for NO taper. The frame route passes the
                   rotation footprint n(|cos t| + |sin t|). None keeps the
                   plain three-shear call, which the validation scripts use as
                   the control.

    Returns:
        A real array of the shape of `patch`.

    Raises:
        ValueError: the patch is not square.
    """
    xpm = xp()
    if patch.ndim != 2 or patch.shape[0] != patch.shape[1]:
        raise ValueError(f"rotate_fourier: the patch must be square, not "
                         f"{patch.shape}.")
    a = float(theta_rad)
    # Reduce by whole quarter turns. `_quarter` is the EXACT index map of a
    # quarter turn about the same centre as the shears, pixel `m // 2`.
    for _ in range(int(np.rint(a / (0.5 * np.pi))) % 4):
        patch = _quarter(patch)
    a -= np.rint(a / (0.5 * np.pi)) * 0.5 * np.pi
    if a == 0.0:
        return patch
    m = patch.shape[0]
    c = (xpm.arange(m) - m // 2).astype(patch.dtype)
    t = np.tan(0.5 * a)
    s = np.sin(a)
    w = None if taper is None else _taper_1d(m, taper, patch.dtype)
    out = patch
    for shift, axis in ((t * c, 1), (-s * c, 0), (t * c, 1)):
        if w is not None:
            out = out * (w[None, :] if axis == 1 else w[:, None])
        out = _shear(out, shift, axis)
    return out


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
        rotated:          True holds ONE THIN strip along the RESULTANT
                          velocity of each layer, and it turns every frame
                          back with `rotate_fourier`. False keeps the
                          axis-aligned box, bit for bit. None (the DEFAULT)
                          is AUTO: see `resolve_rotated`. See the module
                          docstring.
        rot_margin:       the ROLL-OFF of the end taper of the rotated route,
                          in units of the frame side n, on EACH side, PAST the
                          geometric rotation footprint of the layer. The
                          default 0.07 is the measured value. See the module
                          docstring.
    """

    dt_s: float
    n_frames: int
    strip_dir: str
    record: int = 0
    wind_ground_m_s: float = 10.0
    wind_dir_deg: float = 0.0
    slew_rad_s: float = None
    pad_outer_scales: float = 2.0
    rotated: bool = False
    rot_margin: float = 0.07

    def key(self):
        """Give a stable string of the spec, WITHOUT `strip_dir`.

        The campaign fingerprint reads this. The strips are a deletable cache
        that a seed rebuilds, so WHERE they sit must not name the campaign.

        THE ROTATED TAIL. `rotated` and `rot_margin` enter the key ONLY when
        `rotated` is True. That is the append-only rule of `precision` and of
        `fft_backend`, and it keeps every key of an older record valid.

        Returns:
            A string.
        """
        tail = () if self.rotated else ("rotated", "rot_margin")
        parts = [f"{f.name}={getattr(self, f.name)!r}"
                 for f in fields(self)
                 if f.name != "strip_dir" and f.name not in tail]
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
        theta:    the angle of the resultant velocity of each layer, in rad,
                  or None for the along-track route.
        m_crop:   the side of the padded crop of the rotated route, in pixels,
                  ONE value for EACH layer, or None for the along-track route.
    """

    v_x_m_s: np.ndarray
    v_y_m_s: np.ndarray
    shape: tuple
    dt_s: float
    n_frames: int
    dx: float
    n: int
    pad: int
    theta: np.ndarray = None
    m_crop: tuple = None


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

    if getattr(spec, "rotated", False):
        # THE ROTATED ROUTE. The long axis of each strip runs along the
        # resultant velocity, so the walk is purely along x at the speed
        # |v|. The window is the PADDED crop that `rotate_fourier` needs: it
        # holds the rotation footprint n(|cos|+|sin|) of the frame and the
        # taper roll-off `rot_margin` on each side. EACH LAYER GETS ITS OWN
        # crop side, because a layer that walks along x asks for much less
        # than a layer that walks at 45 deg, and the crop pixels drive the
        # cost of a frame.
        theta = np.arctan2(v_y, v_x)
        speed = np.hypot(v_x, v_y)
        foot = np.abs(np.cos(theta)) + np.abs(np.sin(theta))
        shape = []
        m_crop = []
        for v, f in zip(speed, foot):
            m = int(np.ceil(n * (f + 2.0 * float(spec.rot_margin))))
            m += m % 2                         # an even side centres cleanly.
            nx = m + int(np.ceil(v * total_s / dx)) + pad
            nx = STRIP_ALIGN * int(np.ceil(nx / STRIP_ALIGN))
            shape.append((m, nx))
            m_crop.append(m)
        return StripPlan(v_x_m_s=v_x, v_y_m_s=v_y, shape=tuple(shape),
                         dt_s=float(spec.dt_s), n_frames=int(spec.n_frames),
                         dx=dx, n=n, pad=pad, theta=theta,
                         m_crop=tuple(m_crop))

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

    THE ROTATED ROUTE walks along the strip x at the RESULTANT speed |v|, and
    the crop stays centred on the short axis. The pair is then the corner of
    the PADDED crop, of side `sp.m_crop`, not of the frame.

    Returns:
        A list of (oy, ox) integer pairs, one for each layer.
    """
    t = float(k) * sp.dt_s
    if sp.m_crop is not None:
        speed = np.hypot(sp.v_x_m_s, sp.v_y_m_s)
        return [(0, int(np.rint(v * t / sp.dx))) for v in speed]
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
                 dtype=np.float32, table_dtype=None):
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
        table_dtype:  the `ScreenFactory` table type. None (the default)
                      keeps the float64 build of record; numpy.float32 halves
                      the peak memory of a large (crosswind) strip.

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
                                    dtype=dtype, table_dtype=table_dtype)
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


def open_strips(paths, on_device=False):
    """Open the strip of every layer as a read-only memory map.

    A memory map shares ONE page cache between the worker processes of a pool,
    so a 12-worker campaign holds one copy of the strips, not twelve.

    `on_device` holds the strips on the CUDA device instead. The ROTATED route
    turns every frame with three FFTs, so a host memory map makes the run
    upload one crop of every layer for EVERY frame. One upload of the whole
    strip removes that traffic. The call falls back to the memory map, with a
    warning, when the strips do not fit the free device memory.

    Args:
        paths:     the file path of each layer.
        on_device: True uploads the strips to the CUDA device, when the FFT
                   backend is the device and the strips fit.

    Returns:
        A list of read-only numpy memory maps, or a list of device arrays.
    """
    strips = [np.load(p, mmap_mode="r") for p in paths]
    xpm = xp()
    if not on_device or xpm is np or not strips:
        return strips
    need = sum(int(s.nbytes) for s in strips)
    free = int(xpm.cuda.Device().mem_info[0])
    # HALF the free memory: the field, the screens and the FFT plans of the
    # run need the other half.
    if need > 0.5 * free:
        warnings.warn(
            f"open_strips: the strips hold {need / 2**30:.2f} GiB and the "
            f"device has {free / 2**30:.2f} GiB free, so they stay on the "
            "host. Every frame then uploads its crops. Shorten the record or "
            "lower the frame count.", RuntimeWarning)
        return strips
    return [xpm.asarray(np.ascontiguousarray(s)) for s in strips]


def frame_stack(strips, sp, k):
    """Give the phase screen of every layer at frame k.

    Each screen is a VIEW of the memory map, so the call copies nothing. A
    device run makes the copy itself, when it uploads the view.

    THE ROTATED ROUTE (`TemporalSpec.rotated`) copies. It takes the padded
    crop of the side `sp.m_crop` OF THAT LAYER, it turns the crop back by the
    layer angle with `rotate_fourier` under the end taper, and it cuts the
    central n by n. The crop goes to the array module of the FFT backend
    first, so the shears run where the field is. A strip that already sits on
    the device (`open_strips(..., on_device=True)`) needs no upload.

    Args:
        strips: the open strips. Use `open_strips`.
        sp:     the StripPlan.
        k:      the frame index, from 0.

    Yields:
        One n by n array of the phase, in radians, for each layer, in the
        order of the plan.
    """
    n = sp.n
    if sp.m_crop is None:
        for strip, (oy, ox) in zip(strips, frame_offsets(sp, k)):
            yield strip[oy:oy + n, ox:ox + n]
        return
    xpm = xp()
    for strip, (oy, ox), theta, m in zip(strips, frame_offsets(sp, k),
                                         sp.theta, sp.m_crop):
        m = int(m)
        lo = (m - n) // 2
        patch = xpm.asarray(strip[oy:oy + m, ox:ox + m])
        # out(p) = patch(R(-theta) p): the frame axes turn back onto the lab
        # axes, and the strip x axis IS the velocity of the layer. The flat
        # top of the taper is the rotation footprint of the frame, so the cut
        # window reads no tapered pixel.
        flat = n * (abs(np.cos(theta)) + abs(np.sin(theta)))
        yield rotate_fourier(patch, -float(theta),
                             taper=flat)[lo:lo + n, lo:lo + n]


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
        # THE ALONG-TRACK FRAME IS THE PLAIN SLICE. The rotated route must
        # not move the default route by one value.
        for strip, frame, (oy, ox) in zip(strips, f0, frame_offsets(sp_s, 0)):
            assert np.array_equal(frame, strip[oy:oy + 32, ox:ox + 32])
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

    # ---- 8. the rotated strip (the OPT-IN) ----
    # 8a. the default route does not change: `rotated` is not in the key.
    assert "rotated" not in k1, k1
    assert TemporalSpec(dt_s=1e-3, n_frames=4, strip_dir="/a",
                        rot_margin=9.0).key() == k1
    assert sp.theta is None and sp.m_crop is None

    # 8b. THE GEOMETRY. out(p) = in(R(a) p), so a point at q moves to the
    # pixel R(-a) q. A point source is the sharpest test of the angle, of the
    # sign and of the centre of the three shears.
    m_t = 64
    point = np.zeros((m_t, m_t))
    point[m_t // 2, m_t // 2 + 10] = 1.0            # x = 10, y = 0.
    for angle in (0.3, -0.75, 1.9, 0.5 * np.pi, -np.pi + 0.2):
        got = rotate_fourier(point, angle)
        i, j = np.unravel_index(int(np.argmax(np.abs(got))), got.shape)
        want = (-10.0 * np.sin(angle), 10.0 * np.cos(angle))   # (y, x).
        assert abs(i - m_t // 2 - want[0]) <= 0.5, (angle, i, want)
        assert abs(j - m_t // 2 - want[1]) <= 0.5, (angle, j, want)

    # 8c. THE ROUND TRIP. A rotation and its inverse give the field back. The
    # cut window must sit inside the margin, because a shear wraps the edge.
    xx, yy = np.meshgrid(np.arange(m_t) - m_t // 2, np.arange(m_t) - m_t // 2)
    rng_t = np.random.default_rng(17)
    smooth = np.zeros((m_t, m_t))
    for _ in range(8):                              # a few low tones.
        fxy = rng_t.integers(-4, 5, size=2) / m_t
        smooth += np.cos(2.0 * np.pi * (fxy[0] * xx + fxy[1] * yy)
                         + rng_t.uniform(0.0, 6.28))
    back = rotate_fourier(rotate_fourier(smooth, 0.4), -0.4)
    q = m_t // 4
    core = (slice(q, m_t - q), slice(q, m_t - q))
    rel = float(np.sqrt(np.mean((back[core] - smooth[core]) ** 2)
                        / np.mean(smooth[core] ** 2)))
    assert rel < 0.02, rel
    # 8c. a quarter turn and the identity.
    flat = np.random.default_rng(3).standard_normal((32, 32))
    assert np.allclose(rotate_fourier(flat, 0.0), flat)
    assert np.array_equal(rotate_fourier(flat, 0.5 * np.pi), _quarter(flat))
    assert np.array_equal(rotate_fourier(flat, 2.0 * np.pi), flat)

    # 8d. the plan: one thin strip per layer, along the resultant velocity.
    spec_rot = TemporalSpec(dt_s=1e-3, n_frames=11, strip_dir="unused",
                            wind_dir_deg=90.0, rotated=True, rot_margin=0.5)
    sp_rot = strip_plan(plan, grid, spec_rot, geometry, L0)
    spec_box = TemporalSpec(dt_s=1e-3, n_frames=11, strip_dir="unused",
                            wind_dir_deg=90.0)
    sp_box = strip_plan(plan, grid, spec_box, geometry, L0)
    # ONE crop side for EACH layer: 64*(foot_j + 2*0.5).
    assert len(sp_rot.m_crop) == len(sp_rot.shape), sp_rot.m_crop
    for m_j, th_j, (ny, _) in zip(sp_rot.m_crop, sp_rot.theta, sp_rot.shape):
        foot_j = abs(np.cos(th_j)) + abs(np.sin(th_j))
        want_j = int(np.ceil(64 * (foot_j + 2.0 * 0.5)))
        assert m_j in (want_j, want_j + 1), (m_j, want_j)
        assert m_j % 2 == 0 and ny == m_j, (m_j, ny)
    # The layers do NOT share one crop side: the ground layer barely turns
    # and the jet layer turns a lot, so the sides differ.
    assert len(set(sp_rot.m_crop)) > 1, sp_rot.m_crop
    # The DEFAULT roll-off gives a much smaller crop than the old 0.5 margin.
    spec_thin = TemporalSpec(dt_s=1e-3, n_frames=11, strip_dir="unused",
                             wind_dir_deg=90.0, rotated=True)
    sp_thin = strip_plan(plan, grid, spec_thin, geometry, L0)
    assert all(a < b for a, b in zip(sp_thin.m_crop, sp_rot.m_crop)), \
        (sp_thin.m_crop, sp_rot.m_crop)
    # The rotated strip holds FEWER pixels than the axis-aligned box.
    area_rot = sum(a * b for a, b in sp_rot.shape)
    area_box = sum(a * b for a, b in sp_box.shape)
    assert area_rot < area_box, (area_rot, area_box)

    # 8e. build, crop and rotate: the frames are n by n and they move.
    tmpdir = tempfile.mkdtemp(prefix="olb_rot_")
    try:
        spec_r = TemporalSpec(dt_s=1e-3, n_frames=4, strip_dir=tmpdir,
                              pad_outer_scales=0.02, slew_rad_s=0.0,
                              wind_dir_deg=40.0, rotated=True)
        grid_r = types.SimpleNamespace(n=32, pixel_m=0.02)
        sp_r = strip_plan(plan, grid_r, spec_r, geometry, L0)
        paths_r = strip_paths(spec_r, 3)
        build_strips(sp_r, plan.r0_m, paths_r, L0, True, [11, 22, 33])
        strips_r = open_strips(paths_r)
        g0 = list(frame_stack(strips_r, sp_r, 0))
        g3 = list(frame_stack(strips_r, sp_r, spec_r.n_frames - 1))
        assert all(f.shape == (32, 32) for f in g0), [f.shape for f in g0]
        assert not np.array_equal(g0[-1], g3[-1]), 'the atmosphere must move'
        # A rotated frame is a real phase screen: its variance is finite and
        # its structure function grows with the lag.
        d1 = float(np.mean((g0[-1][1:, :] - g0[-1][:-1, :]) ** 2))
        d4 = float(np.mean((g0[-1][4:, :] - g0[-1][:-4, :]) ** 2))
        assert 0.0 < d1 < d4, (d1, d4)
        del strips_r, g0, g3
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)

    # 8f. THE TAPER GATE. The production crop (the geometric footprint plus
    # the 0.07 n roll-off) with the end taper must give the frame that the
    # rotation of the FULL PERIODIC screen gives. That screen is periodic, so
    # it carries no seam and it IS the reference.
    n_f, n_big, ring = 64, 256, 8
    big = ScreenFactory(n_big, 0.01, L0_m=25.0).make(
        0.10, np.random.default_rng(23))

    def _cut(a, side):
        """Cut the central square about the rotation centre `side // 2`."""
        c = a.shape[0] // 2
        return a[c - side // 2:c + side // 2, c - side // 2:c + side // 2]

    def _d_ratio(a, ref):
        """Give the worst |D(lag)/D_ref(lag) - 1| over the axes and lags."""
        worst = 0.0
        for lag in (1, 4, 16):
            for ax in (0, 1):
                sa = ((slice(lag, None), slice(None)) if ax == 0
                      else (slice(None), slice(lag, None)))
                sb = ((slice(None, -lag), slice(None)) if ax == 0
                      else (slice(None), slice(None, -lag)))
                d = float(np.mean((a[sa] - a[sb]) ** 2))
                dr = float(np.mean((ref[sa] - ref[sb]) ** 2))
                worst = max(worst, abs(d / dr - 1.0))
        return worst

    yy, xx = np.indices((n_f, n_f))
    edge = np.minimum.reduce([xx, yy, n_f - 1 - xx, n_f - 1 - yy]) < ring
    taper_rows = []
    for deg in (20.0, 40.0):
        th = np.deg2rad(deg)
        reference = _cut(rotate_fourier(big, th), n_f)
        foot = abs(np.cos(th)) + abs(np.sin(th))
        m_f = int(np.ceil(n_f * (foot + 2.0 * 0.07)))
        m_f += m_f % 2
        crop = _cut(big, m_f)
        got = _cut(rotate_fourier(crop, th, taper=n_f * foot), n_f)
        plain = _cut(rotate_fourier(crop, th), n_f)
        r_taper = float(np.sqrt(np.mean((got - reference)[edge] ** 2)))
        r_plain = float(np.sqrt(np.mean((plain - reference)[edge] ** 2)))
        assert r_taper < 0.2 * r_plain, (deg, r_taper, r_plain)
        assert r_taper < 1e-2, (deg, r_taper)
        assert _d_ratio(got, reference) < 0.01, (deg,
                                                 _d_ratio(got, reference))
        taper_rows.append((deg, m_f, r_plain, r_taper))

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
    print(f"rotated crop sides          {sp_thin.m_crop} px "
          f"(n = {sp_thin.n}, roll-off 0.07 n)")
    print("the taper gate, against the rotation of the full periodic screen")
    print(f"{'deg':>6} {'crop':>6} {'ring no taper':>15} {'ring taper':>12}")
    for deg, m_f, r_plain, r_taper in taper_rows:
        print(f"{deg:6.1f} {m_f:6d} {r_plain:15.3e} {r_taper:12.3e}")
    print("self-check passed")
