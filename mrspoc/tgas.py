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
                            'tgas_bright_g_lt_12.tsv')


def get_table_ms(plot=True, ax=None):
    """
    Open the TGAS catalog for all stars brighter than G < 12,
    make CMD cuts to flag just the main sequence stars.

    Parameters
    ----------
    plot : bool (optional)
        Make a plot of the CMD and color/mag cuts.
    ax : `~matplotlib.pyplot.Axes` (optional)
        If ``ax`` is not `None`, make the plot on ``ax``.

    Returns
    -------
    table : `~astropy.table.Table`
        Table of TGAS sources, magnitudes, parallaxes, distances, and
        anticipated astrometric uncertainties.
    ms : `~numpy.ndarray`
        Boolean array, ``True`` for rows of ``table`` where the star is a
        main sequence star within the color/mag cuts.
    """
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
    table.add_column(Column(data=N_fov(abs_galactic_latitude), name='N_fov'))

    M_V = Vmag - 5*(np.log10(dist_pc) + 1)

    b_minus_v_lower = 0.6 # 0.64  # (B-V)_sun = 0.65
    b_minus_v_upper = 2

    main_sequence = ((np.abs(M_V - color_cut(b_minus_v)) < 1.) &
                     (b_minus_v > b_minus_v_lower) &
                     (b_minus_v < b_minus_v_upper))

    if plot:
        if ax is None:
            ax = plt.gca()
        polygon_x = [0.6, 0.6, 2.0, 2.0, 0.6]
        polygon_y = [color_cut(0.6) - 1, color_cut(0.6) + 1,
                     color_cut(2) + 1, color_cut(2) - 1,
                     color_cut(0.6) - 1]

        H, xedges, yedges = np.histogram2d(b_minus_v, M_V, bins=1000)

        extent = [xedges.min(), xedges.max(), yedges.max(), yedges.min()]
        ax.imshow(np.log10(H.T), extent=extent, cmap=plt.cm.Greys, aspect=0.2)
        ax.plot(polygon_x, polygon_y, lw=2, color='r', ls='--')

        ax.set(xlim=[-0.5, 3], ylim=[2, -15],
               ylabel='$M_{VT}$', xlabel="BT - VT")

    return table, main_sequence
