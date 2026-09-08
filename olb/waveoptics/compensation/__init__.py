"""Perfect adaptive optics for the fidelity-2 wave-optics layer.

The package removes the first N Zernike modes of the wavefront over one
receive aperture. "Perfect" means an IDEAL modal fit: there is no
wavefront-sensor noise, no finite subaperture, no aliasing, no servo lag and
no branch point. The residual is the pure fitting error of the higher modes.
So the package gives the UPPER BOUND of the benefit of a corrector, and it is
the correct reference for the adaptive-optics fidelity-1 Terms.

The correction is a SNAPSHOT. There is no time axis, no fade rate and no fade
duration. The whole fidelity-2 layer is a snapshot layer.

THE ORDER IS NOLL, NOT ANSI/OSA. The index starts at j = 1 (piston), j = 2 is
the x tilt, j = 3 is the y tilt, and j = 4 is defocus. This is the order of
R. J. Noll, J. Opt. Soc. Am. 66, 207 (1976), DOI 10.1364/JOSA.66.000207,
Table I. It is also the order that `olb.turbulence.ao` counts: a TipTilt stage
removes the first 3 Noll modes, and an AO(n) stage removes the first n Noll
modes. The two ladders therefore count the same modes.

The package is pure numpy. It imports no other olb module, except the leaf
`olb.waveoptics.turbulence.screens` in a self-check.

Modules:
    zernike.py  the Noll modes and the circular mask.
    modal.py    the ApertureModes projector and the stack-to-mode-count map.
    slopes.py   the wrapped-gradient slopes and their modal fit.
    gtilt.py    the gradient (centroid) tilt, from the far-field centroid.
"""

from .gtilt import far_field_tilt
from .modal import ApertureModes, modes_from_stack
from .slopes import SlopeReconstructor, max_abs_step, wrapped_gradient
from .zernike import circle, mask_radius_px, noll_to_nm, zernike_basis, zernike_j

__all__ = [
    "ApertureModes",
    "SlopeReconstructor",
    "circle",
    "far_field_tilt",
    "mask_radius_px",
    "max_abs_step",
    "modes_from_stack",
    "noll_to_nm",
    "wrapped_gradient",
    "zernike_basis",
    "zernike_j",
]
