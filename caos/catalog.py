# Licensed under the MIT License - see LICENSE.rst

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import os
import numpy as np
from astropy.io import ascii
from astropy.table import Column
import astropy.units as u

__all__ = ['get_table_ms']


catalog_path = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                            os.path.pardir, 'data',
                            'tgas_bright.tsv')

#                            'tgas_high_precision_dra_0.08.tsv')


def get_table_ms():
    table = ascii.read(catalog_path, delimiter=';', data_start=3)
    # floatify:
    table['BTmag'] = table['BTmag'].astype(float)
    table['VTmag'] = table['VTmag'].astype(float)

    def color_cut(b_minus_v):
        return -10 + 6.2*b_minus_v

    parallax_mas = table['Plx']
    Vmag = table['VTmag']
    b_minus_v = table['BTmag'] - table['VTmag']

    parallax_arcsec = parallax_mas / 1000
    dist_pc = 1./parallax_arcsec
    table.add_column(Column(data=dist_pc * u.pc, name='distance'))
    M_V = Vmag - 5*(np.log10(dist_pc) + 1)

    b_minus_v_sun = 0.653
    b_minus_v_hat11 = 1.2

    b_minus_v_lower = b_minus_v_sun - 0.1
    b_minus_v_upper = b_minus_v_hat11 + 0.1

    main_sequence = ((np.abs(M_V - color_cut(b_minus_v)) < 1.5) &
                     (b_minus_v > b_minus_v_lower) &
                     (b_minus_v < b_minus_v_upper))

    return table, main_sequence