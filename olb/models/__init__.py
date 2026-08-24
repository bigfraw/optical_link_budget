'''
Physics models: each turns a scenario + geometry into one or more Terms.

Every public factory here has the same shape so the budget assembler can call
them uniformly:

    def <term>(scenario, geometry, **kwargs) -> Term

Most factories serve every link with the same code (uplink, downlink, and
retro). Some are link-specific: the coupling package holds a downlink
coupling Term and separate terrestrial coupling Terms, because a slant plane-wave
beam and a horizontal Gaussian beam read different physics. Inside the coupling
package the downlink and terrestrial modules import shared single-mode-fibre
helpers from coupling._common. That is intra-package reuse; it does not import
another model family.

Models never import each other across families or the budget; they only read the
scenario (a SpaceScenario or a TerrestrialScenario) and geometry and return Terms.
'''
