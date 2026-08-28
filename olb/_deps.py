'''
The only module that imports the physics kernels from my_analysis_modules.

Every other module in olb imports its physics from here. This is the only module
that sets the path to my_analysis_modules. No other module sets this path. This
removes the per-file `sys.path.append(r"D:\\repos\\my_analysis_modules")` code
that the old tn2_kepler scripts used in each file. To vendor, pip-install, or
move my_analysis_modules, change only this file.

Set the location with the MY_ANALYSIS_MODULES environment variable if the repo
is in a different place.
'''

import os
import io
import sys
import pathlib
import contextlib

_DEFAULT = r"D:\repos\my_analysis_modules"
_MAM = pathlib.Path(os.environ.get("MY_ANALYSIS_MODULES", _DEFAULT))

if not _MAM.exists():
    raise ImportError(
        f"my_analysis_modules not found at {_MAM!s}. Set the "
        "MY_ANALYSIS_MODULES environment variable to its location."
    )

# The kernels use flat imports (e.g. `from fields import gaussz`), so the
# package directory itself goes on the path, not its parent.
if str(_MAM) not in sys.path:
    sys.path.insert(0, str(_MAM))

# Re-export the exact symbols olb borrows. If any of these move, this is the
# only import that breaks -- a deliberate single point of failure. (The unit
# conversions, the Gaussian-beam gaussz/zR, and the Satellite geometry, once
# borrowed here, now live in olb.units, olb.beam, and olb.geometry.)
# coupled_flux prints text at import. It is in a shared repo that we do not own.
# This code stops the print here. It does not edit that module. (The Hufnagel-
# Valley Cn2 model and the Bufton wind model, once borrowed here from
# general_atmospherics, now live in olb.turbulence.profiles, so olb no longer
# imports general_atmospherics.)
with contextlib.redirect_stdout(io.StringIO()):
    from coupled_flux import (                                   # noqa: E402
        coupled_flux_montecarlo,
        # Lower-level kernels: olb.turbulence.uplink_flux composes these into a
        # short uplink MC loop with a diverged free-space beam width (w_free).
        spherical_wave_coherence_diameter, short_term_beam_waist,
        long_term_beam_waist, beam_wander_variance,
        coupled_flux_sample, on_axis_irradiance,
    )

__all__ = [
    "coupled_flux_montecarlo",
    "spherical_wave_coherence_diameter", "short_term_beam_waist",
    "long_term_beam_waist", "beam_wander_variance",
    "coupled_flux_sample", "on_axis_irradiance",
]
