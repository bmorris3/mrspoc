import numpy as np
import pytest

from ..star import Star, Spot


@pytest.mark.parametrize('x', (0.8, 0.6, 0.4, 0.2, 0.0))
def test_centroid_small_spots(x):
    """
    python -c "from mrspoc.tests.test import test_centroid_small_spots as f; f()"
    """
    star = Star()
    star.spots = [Spot(x=x, y=0.0, r=0.01)]

    x0, y0 = star._centroid_analytic()
    x1, y1 = star._centroid_numerical()

    diff = np.sqrt((x0 - x1)**2 + (y0 - y1)**2)
    assert diff < 0.00001


@pytest.mark.parametrize('x', (0.8, 0.6, 0.4, 0.2, 0.0))
def test_centroid_large_spots(x):
    """
    python -c "from mrspoc.tests.test import test_centroid_large_spots as f; f()"
    """
    star = Star()
    star.spots = [Spot(x=x, y=0.0, r=0.1)]

    x0, y0 = star._centroid_analytic()
    x1, y1 = star._centroid_numerical()

    diff = np.sqrt((x0 - x1)**2 + (y0 - y1)**2)

    # Allow for bigger error because large spots are more affected by difference
    # in limb-darkening across the spot.
    assert diff < 0.0005
