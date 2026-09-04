"""The turbulent split-step layer of the wave-optics package.

The layer moves a complex field along a path, and it puts a random phase
screen at each slab of that path. It gives one SNAPSHOT of the atmosphere
for each seed. It carries no time axis and it builds no Term, so it changes
no link budget.

The modules are:

    screens     one phase screen, and how to put it into a field.
    splitstep   the propagate-screen-propagate loop, and the boundary mask.
    sampling    the turbulent grid sizer, and the screen-placement planner.
    run         the trial runner: one snapshot for each seed.
    campaign    a large set of trials on disk, stored as blocks.
    temporal    the frozen-flow time axis. PLANNED, NOT BUILT.

Sources:
- Schmidt, Numerical Simulation of Optical Wave Propagation with Examples in
  MATLAB, DOI 10.1117/3.866274, Ch. 9. The split-step method.
"""

from .campaign import Campaign
from .run import (TurbTrial, TurbWaveResult, folded_terrestrial,
                  propagate_turbulent_field, propagate_turbulent_scenario)
from .sampling import (PRESETS, QualityPreset, SamplingReport, ScreenPlan,
                       turbulent_grid)
from .screens import Screen, phase_screen, screen_r0
from .splitstep import split_step, super_gaussian_boundary
from .temporal import TemporalScreens

__all__ = [
    'Campaign',
    'PRESETS',
    'QualityPreset',
    'SamplingReport',
    'Screen',
    'ScreenPlan',
    'TemporalScreens',
    'TurbTrial',
    'TurbWaveResult',
    'folded_terrestrial',
    'phase_screen',
    'propagate_turbulent_field',
    'propagate_turbulent_scenario',
    'screen_r0',
    'split_step',
    'super_gaussian_boundary',
    'turbulent_grid',
]
