"""Fidelity-2 field propagation without turbulence, in the LightPipes API.

The package propagates a scalar complex field through free space on a square
grid. It gives the no-turbulence validator for the near-field and far-field
flags of the link budget. The interface keeps the LightPipes names and the
LightPipes call order, so a script from that package runs here with no
change.

Ported and trimmed from LightPipes (https://github.com/opticspy/lightpipes),
BSD-3-Clause. See LIGHTPIPES_LICENSE.txt in this package.

The core (field.py, sources.py, propagators.py, lenses.py, smf.py) is pure
physics. It imports numpy and scipy only. It imports nothing from the rest of
olb.

lenses.py holds the thin lens and the spherical (co-moving) coordinate route.
A long space link makes the beam grow by a factor of 100 or more. A flat grid
cannot hold that beam AND resolve the launch aperture. LensFresnel() moves
the grid with the beam, and Convert() comes back to a flat grid. See the
module docstring for the three-call recipe.

Two modules on top of that core read an olb scenario:

- grid.py: GridSpec, the grid extent and the grid resolution for a scenario.
           It selects the flat grid or the scaled (co-moving) grid.
- run.py:  propagate_scenario, one end-to-end propagation, and WaveResult.

smf.py holds the single-mode-fibre pupil mode and the coupling efficiency.

The package builds NO Term and it changes NO budget.
"""

from .field import (Begin, Field, Intensity, Normal, Phase, Power,
                    SubIntensity)
from .grid import GridSpec, beam_magnification, forvard_max_z
from .lenses import Convert, Lens, LensForvard, LensFresnel
from .propagators import Forvard, Fresnel, GForvard
from .run import WaveResult, propagate_scenario
from .smf import coupling_efficiency, smf_mode
from .sources import CircAperture, CircScreen, GaussBeam, PlaneWave
from .threader import Threader

__all__ = [
    "Field", "Begin", "Normal", "Power", "Intensity", "Phase", "SubIntensity",
    "GaussBeam", "PlaneWave", "CircAperture", "CircScreen",
    "Forvard", "Fresnel", "GForvard",
    "Lens", "LensForvard", "LensFresnel", "Convert",
    "smf_mode", "coupling_efficiency",
    "GridSpec", "beam_magnification", "forvard_max_z",
    "propagate_scenario", "WaveResult",
    "Threader",
]
