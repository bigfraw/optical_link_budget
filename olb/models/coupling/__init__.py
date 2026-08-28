'''
Receive-coupling Terms for a detector behind an aperture.

This package holds the receive-coupling Terms. The Compensation stack and the
Detector are one physical chain: the residual wavefront that the Compensation
leaves sets the coupling into the Detector. So each link emits one coupling Term.

Some Terms are link-specific: the downlink coupling (downlink.py) reads the
plane-wave slant physics, and the terrestrial coupling (terrestrial.py) reads the
horizontal Gaussian-beam physics. The shared, link-independent single-mode-fibre
coupling physics lives in _common.py, and both link modules import it. The
fidelity-1 FAST modal overlap lives in fast.py.

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
from .fast import smf_fast_term

__all__ = [
    "downlink_coupling_term",
    "terrestrial_smf_coupling_term",
    "terrestrial_smf_walkoff_term",
    "terrestrial_mmf_coupling_term",
    "smf_fast_term",
    "smf_eta_max_from_a",
]
