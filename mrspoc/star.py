# Licensed under the MIT License - see LICENSE.rst

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad
import astropy.units as u
from astropy.coordinates import UnitSphericalRepresentation, CartesianRepresentation
from astropy.coordinates.matrix_utilities import rotation_matrix, matrix_product

from .sun import draw_random_sunspot_latitudes, draw_random_sunspot_radii

__all__ = ['Star', 'Spot']


class Spot(object):
    """
    Properties of a starspot.
    """
    def __init__(self, x=None, y=None, z=None, r=None, contrast=0.7):
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
        if z is None:
            z = np.sqrt(r**2 - x**2 - y**2)
        self.z = z
        self.r = r
        self.contrast = contrast

    @classmethod
    def from_latlon(cls, latitude, longitude, stellar_inclination, radius):
        """
        Construct a spot from latitude, longitude coordinates

        Parameters
        ----------
        latitude : float
            Spot latitude [deg]
        longitude : float
            Spot longitude [deg]
        stellar_inclination : float
            Stellar inclination angle, measured away from the line of sight,
            in [deg].
        radius : float
            Spot radius [stellar radii]
        """

        cartesian = latlon_to_cartesian(latitude, longitude,
                                        stellar_inclination)

        return cls(x=cartesian.x.value, y=cartesian.y.value,
                   z=cartesian.z.value, r=radius)

    @classmethod
    def from_sunspot_distribution(cls, stellar_inclination, mean_latitude=15,
                                  contrast=0.7, radius_multiplier=1):
        """
        Parameters
        ----------
        stellar_inclination : float
            Stellar inclination angle, measured away from the line of sight,
            in [deg].
        mean_latitude : float
            Define the mean absolute latitude of the two symmetric active
            latitudes, where ``mean_latitude > 0``.
        contrast : float (optional)
            Spot contrast relative to photosphere. Default is the area-weighted
            mean sunspot contrast (``c=0.7``).
        """
        lat = draw_random_sunspot_latitudes(n=1, mean_latitude=mean_latitude)[0]
        lon = 2*np.pi * np.random.rand() * u.rad
        radius = draw_random_sunspot_radii(n=1)[0]

        cartesian = latlon_to_cartesian(lat, lon, stellar_inclination)

        return cls(x=cartesian.x.value, y=cartesian.y.value,
                   z=cartesian.z.value, r=radius*radius_multiplier,
                   contrast=contrast)

    def __repr__(self):
        return ("<Spot: x={0}, y={1}, z={2}, r={3}>"
                .format(self.x, self.y, self.z, self.r))


def latlon_to_cartesian(latitude, longitude, stellar_inclination):
    """
    Convert coordinates in latitude/longitude for a star with a given
    stellar inclination into cartesian coordinates.

    The X-Y plane is the sky plane: x is aligned with the stellar equator, y is
    aligned with the stellar rotation axis.

    Parameters
    ----------
    latitude : float or `~astropy.units.Quantity`
        Spot latitude. Will assume unit=deg if none is specified.
    longitude : float or `~astropy.units.Quantity`
        Spot longitude. Will assume unit=deg if none is specified.
    stellar_inclination : float
        Stellar inclination angle, measured away from the line of sight,
        in [deg].

    Returns
    -------
    cartesian : `~astropy.coordinates.CartesianRepresentation`
        Cartesian representation in the frame described above.
    """

    if not hasattr(longitude, 'unit') and not hasattr(latitude, 'unit'):
        longitude *= u.deg
        latitude *= u.deg

    c = UnitSphericalRepresentation(longitude, latitude)
    cartesian = c.to_cartesian()

    rotate_about_z = rotation_matrix(90*u.deg, axis='z')
    rotate_is = rotation_matrix(stellar_inclination*u.deg, axis='y')
    transform_matrix = matrix_product(rotate_about_z, rotate_is)
    cartesian = cartesian.transform(transform_matrix)
    return cartesian


class Star(object):
    """
    Object defining a star.
    """
    def __init__(self, spots=None, u1=0.4987, u2=0.1772, r=1,
                 radius_threshold=0.1, inclination=None):
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
        inclination : `~astropy.units.Quantity`
            Stellar inclination. Default is 90 deg.
        """
        if spots is None:
            spots = []
        self.spots = spots

        self.x = 0
        self.y = 0
        self.r = r
        self.u1 = u1
        self.u2 = u2

        if inclination is None:
            inclination = 90*u.deg

        self._inclination = inclination
        self.radius_threshold = radius_threshold
        self.rotations_applied = 0 * u.deg

    @property
    def inclination(self):
        return self._inclination

    @inclination.setter
    def inclination(self, new_inclination):
        previous_inclination = self._inclination
        rot_new_inc = rotation_matrix(previous_inclination - new_inclination,
                                      axis='x')
        for spot in self.spots:
            cartesian = CartesianRepresentation(x=spot.x, y=spot.y, z=spot.z
                                                ).transform(rot_new_inc)
            spot.x = cartesian.x.value
            spot.y = cartesian.y.value
            spot.z = cartesian.z.value
        self._inclination = new_inclination

    def plot(self, n=3000, ax=None, col=True, col_exaggerate=1):
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
        n : int
            Number of pixels per side in the image.

        Returns
        -------
        ax : `~matplotlib.pyplot.Axes`
            Matplotlib axis object, with the new plot on it.
        """
        if ax is None:
            ax = plt.gca()

        _, _, image = self._centroid_numerical(n=n, return_image=True)

        ax.imshow(image, origin='lower', interpolation='nearest',
                  cmap=plt.cm.Greys_r, extent=[-1, 1, -1, 1])
        ax.set_aspect('equal')
        if col:
            x_col, y_col = self.center_of_light

            ax.scatter([x_col*col_exaggerate], [y_col*col_exaggerate],
                       color='r', marker='x', zorder=100)
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
        large_spots = np.any(np.array([s.r for s in self.spots if s.z > 0]) >
                             self.radius_threshold)

        if large_spots:
            centroid = self._centroid_numerical()
        else:
            centroid = self._centroid_analytic()

        return centroid

    def _centroid_analytic(self, return_total_flux=False):
        """
        Compute the stellar centroid using an analytic approximation.

        Returns
        -------
        x_centroid : float
            Photocenter in the x dimension, in units of stellar radii
        y_centroid : float
            Photocenter in the y dimension, in units of stellar radii
        return_total_flux : bool
            If true, return the X centroid, Y centroid, and the total stellar
            flux.
        """
        x_centroid = 0
        y_centroid = 0

        # Morris et al 2017, Eqn 1
        total_flux = (2 * np.pi *
                      quad(lambda r: r * self.limb_darkening_normed(r),
                           0, self.r)[0])

        for spot in self.spots:
            if spot.z > 0:
                # Morris et al 2017, Eqn 2
                r_spot = np.sqrt(spot.x**2 + spot.y**2)
                spot_area = np.pi * spot.r**2 * np.sqrt(1 - (r_spot/self.r)**2)
                spot_flux = (-1 * spot_area * self.limb_darkening_normed(r_spot) *
                             (1 - spot.contrast))

                # Morris et al 2017, Eqn 3-4
                x_centroid += spot_flux * spot.x
                y_centroid += spot_flux * spot.y
                total_flux += spot_flux
        if not return_total_flux:
            return x_centroid/total_flux, y_centroid/total_flux
        else:
            return x_centroid/total_flux, y_centroid/total_flux, total_flux

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
            if spot.z > 0:
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

        x_centroid = np.sum(image * x)/np.sum(image)
        y_centroid = np.sum(image * y)/np.sum(image)

        if delete_arrays_after_use:
            del on_star
            if len(self.spots) > 0:
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

    def rotate(self, angle):
        """
        Rotate the star, by moving the spots.

        Parameters
        ----------
        angle : `~astropy.units.Quantity`

        """
        remove_is = rotation_matrix(-self.inclination, axis='x')
        rotate = rotation_matrix(angle, axis='z')
        add_is = rotation_matrix(self.inclination, axis='x')

        transform_matrix = matrix_product(remove_is, rotate, add_is)

        for spot in self.spots:
            cartesian = CartesianRepresentation(x=spot.x, y=spot.y, z=spot.z
                                                ).transform(transform_matrix)
            spot.x = cartesian.x.value
            spot.y = cartesian.y.value
            spot.z = cartesian.z.value
        self.rotations_applied += angle

    def derotate(self):
        self.rotate(-self.rotations_applied)
        self.rotations_applied = 0


# rotate_about_z = rotation_matrix(90*u.deg, axis='z')
# rotate_is = rotation_matrix(stellar_inclination*u.deg, axis='y')
# transform_matrix = matrix_product(rotate_about_z, rotate_is)
# cartesian = cartesian.transform(transform_matrix)