# Licensed under the MIT License - see LICENSE.rst

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import os
import numpy as np
from astropy.io import ascii
from astropy.table import Column
from astropy.coordinates import SkyCoord
import astropy.units as u
import matplotlib.pyplot as plt

from .gaia import N_fov, sigma_fov

__all__ = ['get_table_ms']


catalog_path = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                            os.path.pardir, 'data',
                            'tgas_bright_g7.tsv')


def get_table_ms(plot=True, ax=None):
    table = ascii.read(catalog_path, delimiter=';', data_start=3)
    # floatify:
    table['BTmag'] = table['BTmag'].astype(float)
    table['VTmag'] = table['VTmag'].astype(float)

    # Compute the galactic latitude of each star, add to table
    coords = SkyCoord(ra=table['RA_ICRS']*u.deg,
                      dec=table['DE_ICRS']*u.deg, frame='icrs')
    galactic_coords = coords.transform_to('galactic')
    abs_galactic_latitude = abs(galactic_coords.b).degree
    table.add_column(Column(data=abs_galactic_latitude, name='b'))

    # Compute distance, CMD
    def color_cut(b_minus_v):
        return -9. + 4.9 * b_minus_v

    parallax_mas = table['Plx']
    Vmag = table['VTmag']
    b_minus_v = table['BTmag'] - table['VTmag']

    parallax_arcsec = parallax_mas / 1000
    dist_pc = 1./parallax_arcsec

    # Add astrometric uncertainty column to table
    table.add_column(Column(data=sigma_fov(table['<Gmag>']), name='sigma_fov'))

    # Add a distance column to the table:
    table.add_column(Column(data=dist_pc * u.pc, name='distance'))

    # Add a Nfov column to the table:
    table.add_column(Column(data=[N_fov(b) for b in abs_galactic_latitude],
                            name='N_fov'))

    M_V = Vmag - 5*(np.log10(dist_pc) + 1)

    b_minus_v_lower = 0.6 #0.64  # (B-V)_sun = 0.65
    b_minus_v_upper = 2

    main_sequence = ((np.abs(M_V - color_cut(b_minus_v)) < 1.) &
                     (b_minus_v > b_minus_v_lower) &
                     (b_minus_v < b_minus_v_upper))

    if plot:
        if ax is None:
            ax = plt.gca()
        ax.scatter(b_minus_v, M_V, marker='.', s=2)

        x = np.linspace(0.5, 2)
        y = color_cut(x)

        #ax.plot(x, y, 'r', ls='--', alpha=0.1)
        ax.scatter(b_minus_v[main_sequence], M_V[main_sequence], marker='.', s=2, color='r')

        ax.set(xlim=[-0.5, 3], ylim=[0, -20],
               ylabel='$M_{VT}$', xlabel="BT - VT")

    return table, main_sequence