'''
Backward-compatible alias for the retroreflected budget.

The retro model is space-specific (retroreflection as a retransmission over a
long slant range). It now lives in olb.links.retro_space. This module keeps the
old `retro_budget` name so existing call sites do not break. Prefer
`retro_space_budget` in new code.
'''

from .retro_space import retro_space_budget

# Old name kept for backward compatibility. See retro_space.retro_space_budget.
retro_budget = retro_space_budget
