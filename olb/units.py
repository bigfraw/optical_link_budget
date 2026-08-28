'''
dB and linear conversions with small helper functions.

This module OWNS the unit conversions that the whole repo uses: one definition
of todB, fromdB, todBm, fromdBm, and the beam waist-to-divergence conversion
w0_to_div. It also adds helper functions for the sign convention that olb uses:

    LOSS is POSITIVE dB, GAIN is NEGATIVE dB.

The budget adds the value of a Term directly. A +3 dB term adds 3 dB of loss. A
-3 dB term is a gain and removes 3 dB of loss.
'''

import numpy as np

__all__ = ["todB", "fromdB", "todBm", "fromdBm", "w0_to_div", "loss_db",
           "combine_db"]


def todB(x):
    '''Convert a linear power ratio to decibels.'''
    return 10 * np.log10(x)


def fromdB(x):
    '''Convert decibels to a linear power ratio.'''
    return 10 ** (x / 10)


def todBm(x):
    '''Convert a power [W] to dBm (decibels relative to 1 mW).'''
    return 10 * np.log10(x / 1e-3)


def fromdBm(x):
    '''Convert dBm to a power [W].'''
    return 10 ** (x / 10) * 1e-3


def w0_to_div(w0, wavelength=1550e-9):
    '''
    Convert a Gaussian waist radius to the half-angle far-field divergence.

    Parameters:
        w0 : float
            The waist radius of the Gaussian beam [m].
        wavelength : float
            The wavelength [m].

    Returns:
        float
            The half-angle divergence [rad].
    '''
    return wavelength / (np.pi * w0)


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
