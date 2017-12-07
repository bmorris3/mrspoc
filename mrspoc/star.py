# Licensed under the MIT License - see LICENSE.rst

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad
from matplotlib.patches import Ellipse
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
    Object defining a star.
    """
    def __init__(self, u1=0.4987, u2=0.1772, r=1, radius_threshold=0.1,
                 spots=None):
        """
        Parameters
        ----------
        u1 : float (optional)
            Quadratic limb-darkening parameter, linear term
        u2 : float (optional)
            Quadratic limb-darkening parameter, quadratic term
        r : float (optional)
            Stellar radius (default is unity)
        radius_threshold : float (optional)
            If all spots are smaller than this radius, use the analytic solution
            to compute the stellar centroid, otherwise use the numerical
            solution.
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
        self.radius_threshold = radius_threshold

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

        # ax.set_facecolor('k')

        # if ld:
        #     r = np.linspace(0, 1, 100)
        #     Ir = self.limb_darkening_normed(r)
        #     for ri, Iri in zip(r[::-1], Ir[::-1]):
        #         star = plt.Circle((0, 0), ri, color=plt.cm.Greys_r(Iri),
        #                           alpha=1.)
        #         ax.add_artist(star)
        # else:
        #     ax.add_artist(plt.Circle((0, 0), self.r, color='w'))
        #
        # if len(self.spots) > 0:
        #     patches = []
        #     for spot in self.spots:
        #         longitude = np.arcsin(spot.x)
        #         latitude = np.arcsin(spot.y)
        #         width = np.cos(longitude) * spot.r * 2
        #         height = np.cos(latitude) * spot.r * 2
        #         patches.append(Ellipse((spot.x, spot.y), width, height,
        #                                ec='none'))
        #
        #     p2 = PatchCollection(patches, alpha=(1-spot.contrast), color='k',
        #                          zorder=10)
        #     ax.add_collection(p2)
        #
        # if col:
        #     x_col, y_col = self.center_of_light
        #
        #     ax.scatter([x_col*col_exaggerate], [y_col*col_exaggerate],
        #                color='r', marker='x', zorder=100)

        _, _, image = self._centroid_numerical(return_image=True)

        ax.imshow(image, origin='lower', interpolation='nearest',
                  cmap=plt.cm.Greys_r, extent=[-1, 1, -1, 1])
        ax.set_aspect('equal')
        ax.set_xlim([-1, 1])
        ax.set_ylim([-1, 1])
        ax.set_xlabel('x [$R_\star$]', fontsize=14)
        ax.set_ylabel('y [$R_\star$]', fontsize=14)
        return ax

    @property
    def center_of_light(self):
        """
        Compute the center-of-light or centroid for this star,
        given its spots, and limb-darkening.

        Returns
        -------
        x_centroid : float
            Photocenter in the x dimension, in units of stellar radii
        y_centroid : float
            Photocenter in the y dimension, in units of stellar radii
        """
        large_spots = np.any(np.array([s.r for s in self.spots]) >
                             self.radius_threshold)

        if large_spots:
            centroid = self._centroid_numerical()
        else:
            centroid = self._centroid_analytic()

        return centroid

    def _centroid_analytic(self):
        """
        Compute the stellar centroid using an analytic approximation.

        Returns
        -------
        x_centroid : float
            Photocenter in the x dimension, in units of stellar radii
        y_centroid : float
            Photocenter in the y dimension, in units of stellar radii
        """
        x_centroid = 0
        y_centroid = 0

        # Morris et al 2017, Eqn 1
        total_flux = (2 * np.pi *
                      quad(lambda r: r * self.limb_darkening_normed(r),
                           0, self.r)[0])

        for spot in self.spots:
            # spot_longitude = np.arcsin(spot.x)
            # foreshortened_width = np.cos(spot_longitude)
            # the above is equivalent to:
            # foreshortened_width = np.sqrt(self.r**2 - spot.x**2)

            # Morris et al 2017, Eqn 2
            r_spot = np.sqrt(spot.x**2 + spot.y**2)
            spot_area = np.pi * spot.r**2 * np.sqrt(1 - (r_spot/self.r)**2)
            spot_flux = (-1 * spot_area * self.limb_darkening_normed(r_spot) *
                         (1 - spot.contrast))

            # Morris et al 2017, Eqn 3-4
            x_centroid += spot_flux * spot.x
            y_centroid += spot_flux * spot.y
            total_flux += spot_flux

        return x_centroid/total_flux, y_centroid/total_flux

    def _centroid_numerical(self, n=3000, delete_arrays_after_use=True,
                            return_image=False):
        """
        Compute the stellar centroid using a numerical approximation.

        Parameters
        ----------
        n : int
            Generate a simulated image of the star with ``n`` by ``n`` pixels.

        Returns
        -------
        x_centroid : float
            Photocenter in the x dimension, in units of stellar radii
        y_centroid : float
            Photocenter in the y dimension, in units of stellar radii
        """
        image = np.zeros((n, n))
        x = np.linspace(-self.r, self.r, n)
        y = np.linspace(-self.r, self.r, n)
        x, y = np.meshgrid(x, y)

        # Limb darkening
        irradiance = self.limb_darkening_normed(np.sqrt(x**2 + y**2))

        on_star = x**2 + y**2 <= self.r**2

        image[on_star] = irradiance[on_star]

        for spot in self.spots:
            r_spot = np.sqrt(spot.x**2 + spot.y**2)
            foreshorten_semiminor_axis = np.sqrt(1 - (r_spot/self.r)**2)

            a = spot.r  # Semi-major axis
            b = spot.r * foreshorten_semiminor_axis  # Semi-minor axis
            A = np.pi/2 + np.arctan2(spot.y, spot.x)  # Semi-major axis rotation
            on_spot = (((x - spot.x) * np.cos(A) +
                        (y - spot.y) * np.sin(A))**2 / a**2 +
                       ((x - spot.x) * np.sin(A) -
                        (y - spot.y) * np.cos(A))**2 / b**2 <= self.r**2)

            image[on_spot & on_star] *= spot.contrast

            # Validation:
            # r_spot = np.sqrt(spot.x**2 + spot.y**2)
            # image[on_spot & on_star] = spot.contrast * self.limb_darkening_normed(r_spot)

        x_centroid = np.sum(image * x)/np.sum(image)
        y_centroid = np.sum(image * y)/np.sum(image)

        if delete_arrays_after_use:
            del on_star
            del on_spot
            del x
            del y
            del irradiance

        if return_image:
            return x_centroid, y_centroid, image

        if delete_arrays_after_use:
            del image

        return x_centroid, y_centroid

    def limb_darkening(self, r):
        """
        Compute the intensity at radius ``r`` for quadratic limb-darkening law
        with parameters ``Star.u1, Star.u2``.

        Parameters
        ----------
        r : float or `~numpy.ndarray`
            Stellar surface position in radial coords on (0, 1)

        Returns
        -------
        intensity : float
            Intensity in un-normalized units
        """
        mu = np.sqrt(1 - r**2)
        u1 = self.u1
        u2 = self.u2
        return (1 - u1 * (1 - mu) - u2 * (1 - mu)**2) / (1 - u1/3 - u2/6) / np.pi

    def limb_darkening_normed(self, r):
        """
        Compute the normalized intensity at radius ``r`` for quadratic
        limb-darkening law with parameters ``Star.u1, Star.u2``.

        Parameters
        ----------
        r : float or `~numpy.ndarray`
            Stellar surface position in radial coords on (0, 1)

        Returns
        -------
        intensity : float
            Intensity relative to the intensity at the center of the disk.
        """
        return self.limb_darkening(r) / self.limb_darkening(0)
