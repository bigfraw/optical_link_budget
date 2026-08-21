'''
Model assumptions for a Term, and the check that flags a scenario that breaks them.

Each analytic or Monte Carlo model is valid only in a regime. This module gives
a small record that states the regime. A model attaches an Assumptions record to
its Term. The model adds a reason to `violations` when the scenario breaks an
assumption. A Budget then reports which terms are misrepresentative.

Three constraints matter for the optical propagation models:
- beam_type: the wavefront the model assumes (plane wave, spherical wave, or
  Gaussian beam).
- turbulence_regime: the fluctuation strength the model assumes (weak, moderate,
  or strong). The regime is tied to a numeric bound on the scintillation index.
- spectrum: the turbulence spectrum the model assumes (for example Kolmogorov
  with no inner scale and no outer scale).

Use the string constants below so every term uses the same words.
'''

from dataclasses import dataclass, field

# Beam type (wavefront the model assumes).
BEAM_PLANE_WAVE = "plane wave"
BEAM_SPHERICAL_WAVE = "spherical wave"
BEAM_GAUSSIAN = "Gaussian beam"
BEAM_NA = "not applicable"

# Turbulence regime (fluctuation strength the model assumes).
REGIME_WEAK = "weak"
REGIME_MODERATE = "moderate"
REGIME_STRONG = "strong"
REGIME_NA = "not applicable"

# Turbulence spectrum.
SPECTRUM_KOLMOGOROV = "Kolmogorov, no inner or outer scale"
SPECTRUM_VON_KARMAN = "von Karman, finite inner or outer scale"
SPECTRUM_NA = "not applicable"


@dataclass
class Assumptions:
    '''The regime one model is valid in, and the reasons a scenario breaks it.'''
    beam_type: str
    turbulence_regime: str
    spectrum: str
    validity: str = ""                       # the numeric limit, in words
    violations: list = field(default_factory=list)   # reasons the scenario breaks the model

    @property
    def ok(self) -> bool:
        '''Return True when the scenario breaks no assumption.'''
        return len(self.violations) == 0

    def flag(self, reason):
        '''Add one reason that the scenario breaks an assumption. Return self.'''
        self.violations.append(reason)
        return self
