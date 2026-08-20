'''
Physics models: each turns a Scenario + geometry into one or more Terms.

Every public factory here has the same shape so the budget assembler can call
them uniformly:

    def <term>(scenario, geometry, **kwargs) -> Term

Models never import each other or the budget; they only read the Scenario and
geometry and return Terms.
'''
