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
# only import that breaks -- a deliberate single point of failure.
from fields import gaussz, zR                                    # noqa: E402
from conversions import (                                        # noqa: E402
    todB, fromdB, todBm, fromdBm, w0_to_div, div_to_w0,
    arcsec_to_rad, rad_to_arcsec,
)
from satellite import Satellite, SatellitePass                   # noqa: E402

# coupled_flux and general_atmospherics print text at import
# (general_atmospherics does `print(sys.path[0])`). They are in a shared repo
# that we do not own. This code stops the print here. It does not edit those
# modules.
with contextlib.redirect_stdout(io.StringIO()):
    from coupled_flux import coupled_flux_montecarlo             # noqa: E402
    from general_atmospherics import get_c2n, v_wind            # noqa: E402

__all__ = [
    "gaussz", "zR",
    "todB", "fromdB", "todBm", "fromdBm", "w0_to_div", "div_to_w0",
    "arcsec_to_rad", "rad_to_arcsec",
    "Satellite", "SatellitePass",
    "coupled_flux_montecarlo",
    "get_c2n", "v_wind",
]
