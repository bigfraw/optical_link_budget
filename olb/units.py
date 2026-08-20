'''
dB and linear conversions with small helper functions.

This module re-exports the conversions from my_analysis_modules. The whole repo
uses one definition of todB, fromdB, and the other conversions. This module also
adds helper functions for the sign convention that olb uses:

    LOSS is POSITIVE dB, GAIN is NEGATIVE dB.

The budget adds the value of a Term directly. A +3 dB term adds 3 dB of loss. A
-3 dB term is a gain and removes 3 dB of loss.
'''

import numpy as np

from ._deps import todB, fromdB, todBm, fromdBm  # re-exported for convenience

__all__ = ["todB", "fromdB", "todBm", "fromdBm", "loss_db", "combine_db"]


def loss_db(transmission):
    '''
    Loss in dB (positive) for a linear power transmission fraction in (0, 1].

    Parameters:
        transmission : float or array
            The fraction of power that remains (1.0 = no loss).

    Returns:
        float or array
            Loss [dB], positive.
    '''
    return -10.0 * np.log10(transmission)


def combine_db(*terms_db):
    '''
    Sum independent dB contributions (losses positive, gains negative).

    Parameters:
        *terms_db : float or array
            dB values to add.

    Returns:
        float or array
            Total dB.
    '''
    return np.sum(np.broadcast_arrays(*terms_db), axis=0)
