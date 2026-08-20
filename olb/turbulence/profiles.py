'''
Cn2 turbulence profiles for optical link budgets.

This module builds the default zenith Cn2(h) profile from the site parameters.
It also holds the default turbulence altitude grid. It re-exports get_c2n so the
callers read the profile builder from one place.
'''

import numpy as np

from .._deps import get_c2n

DEFAULT_HS = np.logspace(np.log10(1), np.log10(20e3), 20)   # turbulence altitude grid [m]


def default_cn2_profile(site, hs=None):
    '''
    Build a default zenith Cn2 profile from the site parameters.

    Use this profile for the coupled-flux Term when the `fast` package is not
    available. The `fast` HV57 path fails without that package.

    Parameters:
        site : Site
            Provides wind_rms_m_s and cn2_ground.
        hs : numpy.ndarray, optional
            Turbulence altitude grid [m]. Defaults to DEFAULT_HS.

    Returns:
        numpy.ndarray
            Cn2(h) profile at zenith on the hs grid.
    '''
    hs = DEFAULT_HS if hs is None else hs
    return get_c2n(hs, site.wind_rms_m_s, site.cn2_ground)
