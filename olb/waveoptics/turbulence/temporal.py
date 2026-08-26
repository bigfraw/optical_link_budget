"""The frozen-flow time axis of the screen stack. PLANNED, NOT BUILT.

The split-step layer gives independent SNAPSHOTS. Each seed gives one new
atmosphere. A fade DURATION and a fade RATE need a time axis, so they need a
screen stack that MOVES. This module holds the design of that stack. It builds
nothing yet.

The module imports numpy only. It does not import aotools, and it does not
import the rest of olb.

Sources:
- Andrews and Phillips, Laser Beam Propagation through Random Media, 2nd ed.,
  DOI 10.1117/3.626196, Ch. 12, Eqs. (2) and (3). The Bufton wind profile and
  its rms, which set the drift speed of each layer.
- Taylor, The spectrum of turbulence, DOI 10.1098/rspa.1938.0032. The frozen
  flow hypothesis: the turbulence pattern moves with the wind and it does not
  change shape.
- Assemat, Wilson and Gendron, Method for simulating infinitely long and non
  stationary phase screens with optimized memory storage,
  DOI 10.1364/OE.14.000988. The row-by-row extrusion that aotools uses.
"""

import numpy as np


class TemporalScreens:
    """PLANNED, NOT BUILT. The frozen-flow extrusion of a screen stack.

    DESIGN, recorded for the build:

    - Make one aotools infinitephasescreen.PhaseScreenVonKarman for each layer:
      PhaseScreenVonKarman(nx_size=grid.n, pixel_scale=grid.pixel_m,
      r0=plan.r0_m[i], L0=L0_m, random_seed=...). The class holds the state, so
      it gives an infinitely long screen at a fixed memory cost. See Assemat et
      al., DOI 10.1364/OE.14.000988.
    - Call .add_row() to extrude the screen by one column. The number of
      columns for one time step is the drift distance divided by the pixel
      pitch.
    - The drift velocity of a layer is the VECTOR SUM of two parts. The first
      part is the wind of that layer, from the Bufton profile (Andrews and
      Phillips, DOI 10.1117/3.626196, Ch. 12, Eqs. (2) and (3), driven by
      site.wind_rms_m_s). The second part is the apparent translation of a
      tracked satellite: omega_slew multiplied by the slant distance z_i of the
      layer. A LEO pass makes the second part large.
    - The extrusion cost goes with the COLUMN size, so it goes with n. That is
      the reason that the sampling rules of sampling.py keep n low.

    The class raises NotImplementedError. Do not use it.
    """

    def __init__(self, plan, grid, wind_profile_m_s, slew_rad_s=0.0, dt_s=1e-3,
                 seed=None):
        """Raise NotImplementedError. See the class docstring for the design."""
        raise NotImplementedError("temporal evolution is planned, not built")

    def step(self):
        """Raise NotImplementedError. See the class docstring for the design."""
        raise NotImplementedError("temporal evolution is planned, not built")


if __name__ == '__main__':
    # ---- 1. the module holds no import but numpy ----
    import types

    # __main__ adds `builtins` and the `types` of this check. Only `numpy`
    # comes from the module itself.
    mods = sorted(v.__name__ for v in list(globals().values())
                  if isinstance(v, types.ModuleType))
    assert mods == ['builtins', 'numpy', 'types'], mods

    # ---- 2. the constructor raises ----
    try:
        TemporalScreens(None, None, None)
        raise AssertionError('TemporalScreens must raise NotImplementedError')
    except NotImplementedError as exc:
        assert 'planned' in str(exc), str(exc)

    # ---- 3. the design note is present ----
    assert 'PhaseScreenVonKarman' in TemporalScreens.__doc__
    assert 'add_row' in TemporalScreens.__doc__
    assert 'DOI' in TemporalScreens.__doc__

    print("TemporalScreens             PLANNED, NOT BUILT")
    print("constructor                 raises NotImplementedError")
    print("step()                      raises NotImplementedError")
    print("imports                     numpy only")
    print("self-check passed")
