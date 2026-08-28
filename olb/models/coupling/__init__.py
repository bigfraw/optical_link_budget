'''
Receive-coupling Terms for a detector behind an aperture.

This package holds the receive-coupling Terms. The Compensation stack and the
Detector are one physical chain: the residual wavefront that the Compensation
leaves sets the coupling into the Detector. So each link emits one coupling Term.

The modules are named for the LINK, because the coupling physics is
link-specific: downlink.py (downlink_coupling_term) reads the plane-wave slant
physics, and terrestrial.py (the terrestrial_* Terms) reads the horizontal
Gaussian-beam physics. The shared, link-independent single-mode-fibre coupling
physics lives in _common.py, and both link modules import it.

WHY THERE IS NO uplink.py. A coupling Term needs a fibre (or a mode-matched
detector) at the receiver. The downlink and terrestrial receivers are on the
GROUND, so each has a fibre and a coupling module. The UPLINK receiver is the
satellite: a bare light-collecting aperture with no single-mode fibre, so there
is NO uplink coupling Term (the uplink budget builds its vacuum Term with
include_smf=False). The uplink turbulence penalty is not a coupling loss at all;
it is uplink_fast_term (category "turbulence"), in the fidelity-named FAST module
olb.models.fast, read by the uplink budget directly.

WHERE THE FIDELITY-1/2 COUPLING TERMS LIVE. A coupling Term computed by a
fidelity-named module (which spans several Term categories, so it sits at the
olb.models level, not here) is RE-EXPORTED into this package, so every coupling
Term is discoverable here whatever its fidelity: smf_fast_term (fidelity 1, the
FAST modal overlap) comes from olb.models.fast, and waveoptics_smf_coupling_term
(fidelity 2, the split-step fibre coupling) comes from olb.models.waveoptics.

Sources:
  Noll residual variance: R. J. Noll, JOSA 66(3), 207 (1976). See
  olb.turbulence.ao.
  SMF coupling / Marechal: extended Marechal approximation.
  Uncorrected SMF coupling against D/r0: Y. Dikmelik and F. M. Davidson, Appl.
  Opt. 44(23), 4946-4952 (2005), DOI 10.1364/AO.44.004946.
'''

from ._common import smf_eta_max_from_a
from .downlink import downlink_coupling_term
from .terrestrial import (terrestrial_smf_coupling_term,
                          terrestrial_smf_walkoff_term,
                          terrestrial_mmf_coupling_term)
# smf_fast_term and waveoptics_smf_coupling_term are coupling Terms, but their
# implementations live in the fidelity-named modules (olb.models.fast,
# olb.models.waveoptics) beside their non-coupling siblings. The coupling package
# re-exports them, so every coupling Term stays discoverable here whatever its
# fidelity.
from ..fast import smf_fast_term
from ..waveoptics import waveoptics_smf_coupling_term

__all__ = [
    "downlink_coupling_term",
    "terrestrial_smf_coupling_term",
    "terrestrial_smf_walkoff_term",
    "terrestrial_mmf_coupling_term",
    "smf_fast_term",
    "waveoptics_smf_coupling_term",
    "smf_eta_max_from_a",
]
