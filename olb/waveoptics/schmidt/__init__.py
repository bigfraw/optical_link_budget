'''
The Schmidt foundation layer: pure book physics for the wave-optics simulation.

This sub-package holds the numerical method of Schmidt. Andrews and Phillips
gives the ANALYTIC physics. Schmidt gives the NUMERICAL method: the transforms,
the propagators, the sampling constraints, the absorbing boundary, and the phase
screens.

Source of every equation:
    J. D. Schmidt, "Numerical Simulation of Optical Wave Propagation with
    Examples in MATLAB", SPIE Press Monograph PM199 (2010).
    DOI: 10.1117/3.866274
Each function names its chapter, its equation number, and its printed page. The
citation format is:

    Schmidt (2010), DOI 10.1117/3.866274, Ch. N, Eq. (nn), printed p. NNN

The book prints the chapter number in the equation number. So Eq. (3.25) is
equation 25 of Chapter 3. Write the full printed form.

PAGE RULE. The printed page and the page of the Zotero PDF are not the same:

    printed p. N = PDF p. N + 13

This layer holds physics only. It imports numpy and scipy only. It imports
nothing from the rest of olb, and it returns no decibels.

This module exports nothing. A later work package wires the modules.
'''
