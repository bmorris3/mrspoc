# Licensed under the MIT License - see LICENSE.rst

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import os
import numpy as np
from astropy.io import ascii
from astropy.table import Column, join
from astropy.coordinates import SkyCoord
import astropy.units as u
import matplotlib.pyplot as plt
from astropy.constants import R_sun

from .gaia import Nprime_fov, sigma_fov

__all__ = ['get_table_ms']

environment_variable = 'MRSPOC_DATA_DIR'
data_dir_path = os.getenv(environment_variable)

tgas_path = os.path.join(os.path.abspath(data_dir_path),
                         'tgas_bright_g_lt_12.tsv')

hipparcos_path = os.path.join(os.path.abspath(data_dir_path),
                              'hipparcos.tsv')

boyajian_path = os.path.join(os.path.abspath(data_dir_path),
                             'boyajian2012.csv')


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
    table = ascii.read(tgas_path, delimiter=';', data_start=3)

    # floatify:
    table['BTmag'] = table['BTmag'].astype(float)
    table['VTmag'] = table['VTmag'].astype(float)

    # Compute the galactic latitude of each star, add to table
    coords = SkyCoord(ra=table['RA_ICRS'] * u.deg,
                      dec=table['DE_ICRS'] * u.deg, frame='icrs')
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
    dist_pc = 1. / parallax_arcsec

    # Add astrometric uncertainty column to table
    table.add_column(Column(data=sigma_fov(table['<Gmag>']), name='sigma_fov'))

    # Add a distance column to the table:
    table.add_column(Column(data=dist_pc * u.pc, name='distance'))

    # Add a Nfov column to the table:
    table.add_column(Column(data=Nprime_fov(abs_galactic_latitude), name='N_fov'))

    M_V = Vmag - 5 * (np.log10(dist_pc) + 1)

    b_minus_v_lower = 0.6  # 0.64  # (B-V)_sun = 0.65
    b_minus_v_upper = 2

    main_sequence = ((np.abs(M_V - color_cut(b_minus_v)) < 1.) &
                     (b_minus_v > b_minus_v_lower) &
                     (b_minus_v < b_minus_v_upper))

    main_sequence_table = table[main_sequence]

    # Now match the B-V color table from HIPPARCOS to the main sequence TGAS table
    hipparcos_table = ascii.read(hipparcos_path, delimiter=';', header_start=0,
                                 data_start=3)
    hipparcos_table.add_index("HIP")

    main_sequence_table['HIP'][main_sequence_table['HIP'].mask] = 0

    main_sequence_color_table = join(main_sequence_table, hipparcos_table,
                                     keys='HIP')

    # Add in stellar radii with color-radius relation from Boyajian 2012
    R_star = bv_to_radius(main_sequence_color_table['B-V'].data.data)
    main_sequence_color_table.add_column(Column(data=R_star, name='R_star'))

    # Add in a column of interferometric angular diameters from
    # Boyajian 2012 where available:
    boyajian = ascii.read(boyajian_path)
    ang_diams = np.zeros(len(main_sequence_color_table))

    for row in boyajian:
        ang_diams[row['HIP'] == main_sequence_color_table['HIP']] = row['D(UD)']

    main_sequence_color_table.add_column(Column(data=ang_diams,
                                                name='angular_diameter'))

    boyajian_radii = main_sequence_color_table['angular_diameter'] != 0
    half_angle = (main_sequence_color_table['angular_diameter'][boyajian_radii]
                  * u.marcsec/2)
    distance_pc = (main_sequence_color_table['Plx_1'][
                       boyajian_radii].data.data / 1000)**-1 * u.pc
    measured_radii = distance_pc * np.tan(half_angle)

    R_star[boyajian_radii] = measured_radii

    # In radius reference column, `1`==color-radius estimate;
    # `2`==interferometric measurement
    refs = np.ones(len(R_star))
    refs[boyajian_radii] = 2
    main_sequence_color_table.add_column(Column(data=refs, name='rstar_ref'))

    if plot:
        if ax is None:
            ax = plt.gca()
        polygon_x = [0.6, 0.6, 2.0, 2.0, 0.6]
        polygon_y = [color_cut(0.6) - 1, color_cut(0.6) + 1,
                     color_cut(2) + 1, color_cut(2) - 1,
                     color_cut(0.6) - 1]

        H, xedges, yedges = np.histogram2d(b_minus_v[abs(b_minus_v) > 1e-3],
                                           M_V[abs(b_minus_v) > 1e-3],
                                           bins=1000)

        extent = [xedges.min(), xedges.max(), yedges.max(), yedges.min()]
        ax.imshow(np.log10(H.T), extent=extent, cmap=plt.cm.Greys, aspect=0.2)
        ax.plot(polygon_x, polygon_y, lw=2, color='r', ls='--')

        ax.set(xlim=[-0.5, 3], ylim=[2, -15],
               ylabel='$M_{VT}$', xlabel="BT - VT")

    return main_sequence_color_table


def bv_to_radius(b_minus_v):
    """
    Estimate radii for stars on the main sequence using their ``B-V`` color,
    using a simple relation calibrated on interferometry by Boyajian et al. 2012

    Parameters
    ----------
    b_minus_v : float
        B-V color.

    Returns
    -------
    radius : `~astropy.units.Quantity`
        Stellar radius.
    """
    # Boyajian 2012
    X = b_minus_v
    a0 = 0.3830
    a1 = 0.9907
    a2 = -0.6038
    Y = 0
    # Ignore metallicity
    a3 = 0
    a4 = 0
    a5 = 0
    return (a0 + a1 * X + a2 * X ** 2 + a3 * X * Y +
            a4 * Y + a5 * Y ** 2) * R_sun
