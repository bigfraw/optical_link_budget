'''
Physics models: each turns a scenario + geometry into one or more Terms.

Every public factory here has the same shape so the budget assembler can call
them uniformly:

    def <term>(scenario, geometry, **kwargs) -> Term

Models never import each other or the budget; they only read the scenario
(a SpaceScenario or a TerrestrialScenario) and geometry and return Terms.
'''
