# Licensed under the MIT License - see LICENSE.rst

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import numpy as np
import astropy.units as u
from astropy.io import ascii

__all__ = ['N_fov', 'sigma_fov']


perryman2014_table1_str = (
"""    lower upper <N_fov> <Nprime_fov>
    0  5  64.8 51.9
    5  10 65.4 52.3
    10 15 67.2 53.7
    15 20 69.5 55.6
    20 25 73.1 58.5
    25 30 78.3 62.6
    30 35 86.4 69.1
    35 40 99.0 79.2
    40 45 134.3 107.4
    45 50 138.3 110.6
    50 55 106.3 85.0
    55 60 95.2 76.2
    60 65 88.8 71.0
    65 70 84.6 67.7
    70 75 81.5 65.2
    75 80 79.4 63.5
    80 85 78.5 62.8
    85 90 77.6 62.1""")


def N_fov(b):
    """
    Compute the number of visits that Gaia will make (including dead time
    overhead estimate) for a target at galactic latitude ``b`` [deg].
    """
    perryman2014_table1 = ascii.read(perryman2014_table1_str, format='csv',
                                     delimiter=' ')
    cols = ['lower', 'upper', '<Nprime_fov>']
    for lower, upper, Nfov in zip(*[perryman2014_table1[col] for col in cols]):
        if (b < upper) and (b >= lower):
            return Nfov


def sigma_fov(Gmag):
    """
    Approximate Gaia astrometric uncertainty in a single measurement,
    after Perryman 2014, Eqn. 1.
    """
    sigma_att = sigma_cal = 20
    z = 10**(0.4 * (np.max([Gmag, 12*np.ones_like(Gmag)], axis=0) - 15))
    sigma_eta2 = 53000 * z + 310 * z**2
    return np.sqrt(sigma_eta2/9 + sigma_att**2 + sigma_cal**2) * u.uarcsec
