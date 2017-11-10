# Licensed under the MIT License - see LICENSE.rst

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad
from matplotlib.patches import Circle, Ellipse
from matplotlib.collections import PatchCollection


__all__ = ['Star', 'Spot']


def limb_darkening(r, u1=0.4987, u2=0.1772):
    """
    Compute the intensity at radius ``r`` for quadratic limb-darkening law
    with parameters ``u1, u2``.
    """
    mu = np.sqrt(1 - r**2)
    return (1 - u1 * (1 - mu) - u2 * (1 - mu)**2) / (1 - u1/3 - u2/6) / np.pi


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
    def __init__(self, u1=0.4987, u2=0.4987, r=1, spots=None):
        """
        Parameters
        ----------
        u1 : float (optional)
            Quadratic limb-darkening parameter, linear term
        u2 : float (optional)
            Quadratic limb-darkening parameter, quadratic term
        r : float (optional)
            Stellar radius (default is unity)
        spots : list (optional)
            List of spots on this star.
        """
        self.x = 0
        self.y = 0
        self.r = r
        self.u1 = u1
        self.u2 = u2
        if spots is None:
            spots = []
        self.spots = spots

    def plot(self, ax=None, col=True, col_exaggerate=1, ld=True):
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
        ld : bool (optional)
            Show approximation for limb-darkening

        Returns
        -------
        ax : `~matplotlib.pyplot.Axes`
            Matplotlib axis object, with the new plot on it.
        """
        if ax is None:
            ax = plt.gca()

        ax.set_facecolor('k')

        if ld:
            r = np.linspace(0, 1, 100)
            Ir = limb_darkening(r, self.u1, self.u2)/limb_darkening(0)
            for ri, Iri in zip(r[::-1], Ir[::-1]):
                star = plt.Circle((0, 0), ri, color=plt.cm.Greys_r(Iri),
                                  alpha=1.)
                ax.add_artist(star)
        else:
            ax.add_artist(plt.Circle((0, 0), 1, color='w'))

        if len(self.spots) > 0:
            patches = []
            for spot in self.spots:
                longitude = np.arcsin(spot.x)
                width = np.cos(longitude) * spot.r * 2
                height = spot.r * 2
                patches.append(Ellipse((spot.x, spot.y), width, height,
                                       ec='none'))

            p2 = PatchCollection(patches, alpha=(1-spot.contrast), color='k',
                                 zorder=10)
            ax.add_collection(p2)

        if col:
            x_col, y_col = self.center_of_light

            ax.scatter([x_col*col_exaggerate], [y_col*col_exaggerate],
                       color='r', marker='x', zorder=100)

        ax.set_aspect('equal')
        ax.set_xlim([-1, 1])
        ax.set_ylim([-1, 1])
        return ax

    @property
    def center_of_light(self):
        """
        Compute the center-of-light or photometric centroid for this star,
        given its spots, and limb-darkening.

        Returns
        -------
        x_centroid : float
            Photocenter in the x dimension, in units of stellar radii
        y_centroid : float
            Photocenter in the y dimension, in units of stellar radii
        """
        return self._centroid_analytic()

    def _centroid_analytic(self):
        """
        Compute the stellar centroid using an analytic approximation.
        """
        x_centroid = 0
        y_centroid = 0
        total_flux = np.pi * self.r**2 * quad(limb_darkening, 0, 1)[0]

        for spot in self.spots:
            # spot_longitude = np.arcsin(spot.x)
            # foreshortened_width = np.cos(spot_longitude)
            # the above is equivalent to:
            foreshortened_width = np.sqrt(1 - spot.x**2)

            ld_factor = (limb_darkening(np.sqrt(spot.x**2 + spot.y**2)) /
                         limb_darkening(0))

            def _x_weighted(x):
                return - ((1 - spot.contrast) * x * ld_factor *
                          np.sqrt(spot.r**2 - (x - spot.x)**2 /
                                  foreshortened_width**2))

            def _y_weighted(y):
                return - ((1 - spot.contrast) * y * ld_factor *
                          np.sqrt(spot.r**2 - (y - spot.y)**2))

            b = spot.r * foreshortened_width
            a = spot.r
            total_flux -= (1 - spot.contrast) * np.pi * a * b * ld_factor

            x_i = quad(_x_weighted, spot.x - spot.r*foreshortened_width,
                       spot.x + spot.r*foreshortened_width)[0]
            y_i = quad(_y_weighted, spot.y - spot.r, spot.y + spot.r)[0]

            x_centroid += x_i
            y_centroid += y_i

        return x_centroid/total_flux, y_centroid/total_flux

    def _centroid_numerical(self, n=1000):
        """
        Compute the stellar centroid using a numerical approximation.
        """
        image = np.zeros((n, n))
        x = np.linspace(-1, 1, n)
        y = np.linspace(-1, 1, n)
        x, y = np.meshgrid(x, y)

        # Limb darkening
        irradiance = limb_darkening(np.sqrt(x**2 + y**2))/limb_darkening(0)

        on_star = x**2 + y**2 <= self.r

        image[on_star] = irradiance[on_star]

        for spot in self.spots:
            foreshortened_width = np.sqrt(1 - spot.x**2)

            on_spot = ((x - spot.x)**2/foreshortened_width**2 +
                       (y - spot.y)**2 <= spot.r**2)
            image[on_spot & on_star] *= spot.contrast

        x_centroid = np.sum(image * x)/np.sum(image)
        y_centroid = np.sum(image * y)/np.sum(image)
        return x_centroid, y_centroid
