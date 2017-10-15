# Licensed under the MIT License - see LICENSE.rst

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection


__all__ = ['Star', 'Spot']


class Spot(object):
    """
    Properties of a starspot.
    """
    def __init__(self, x, y, r, contrast=0.7):
        """
        Parameters
        ----------
        x : float
            X position [stellar radii]
        y : float
            Y position [stellar radii]
        r : float
            Spot radius [stellar radii]
        contrast : float (optional)
            Spot contrast relative to photosphere
        """
        self.x = x
        self.y = y
        self.r = r
        self.contrast = contrast

    @classmethod
    def from_latlon(cls, latitude, longitude, radius):
        """
        Construct a spot from latitude, longitude coordinates

        Parameters
        ----------
        latitude : float
            Spot latitude [deg]
        longitude : float
            Spot longitude [deg]
        radius : float
            Spot radius [stellar radii]
        """
        x = np.sin(np.radians(longitude))
        y = np.sin(np.radians(latitude))
        return cls(x, y, radius)


class Star(object):
    """
    Object defining a star and its spots.
    """
    def __init__(self):
        self.x = 0
        self.y = 0
        self.spots = []
        self.r = 1

    def plot(self, ax=None, col=True, col_exaggerate=1):
        """
        Plot a 2D projected schematic of the star and its spots.

        Parameters
        ----------
        ax : `~matplotlib.pyplot.Axes`
            Axis object to draw the plot on
        col : bool (optional)
            Show the center of light with a red "x" if `True`
        col_exaggerate : float (optional)
            Exaggerate the center-of-light coordinate by this factor

        Returns
        -------
        ax : `~matplotlib.pyplot.Axes`
            Matplotlib axis object, with the new plot on it.
        """

        if ax is None:
            ax = plt.gca()

        ax.set_facecolor('k')

        patches = []
        for spot in self.spots:
            patches.append(Circle((spot.x, spot.y), spot.r))

        p1 = PatchCollection([Circle((0, 0), 1)], alpha=1, color='w')
        p2 = PatchCollection(patches, alpha=(1-spot.contrast), color='k')
        ax.add_collection(p1)
        ax.add_collection(p2)

        if col:
            x_col, y_col = self.center_of_light

            ax.scatter([x_col*col_exaggerate], [y_col*col_exaggerate],
                       color='r', marker='x')

        ax.set_aspect('equal')
        ax.set_xlim([-1, 1])
        ax.set_ylim([-1, 1])
        return ax

    @property
    def center_of_light(self):
        """
        Compute the center-of-light or photocenter on this star, given its spots.

        Returns
        -------
        x_centroid : float
            Photocenter in the x dimension, in units of stellar radii
        y_centroid : float
            Photocenter in the y dimension, in units of stellar radii
        """
        x_centroid = 0
        y_centroid = 0

        for spot in self.spots:
            def _x_weight(x):
                return - spot.contrast * 2 * x * np.sqrt(spot.r**2 - (x - spot.x)**2)

            def _y_weight(y):
                return - spot.contrast * 2 * y * np.sqrt(spot.r**2 - (y - spot.y)**2)

            x_i = quad(_x_weight, spot.x-spot.r, spot.x+spot.r)[0]
            y_i = quad(_y_weight, spot.y-spot.r, spot.y+spot.r)[0]

            x_centroid += x_i
            y_centroid += y_i

        return x_centroid, y_centroid