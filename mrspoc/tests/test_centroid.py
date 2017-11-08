import numpy as np
import pytest

from ..star import Star, Spot


@pytest.mark.parametrize('x', (0.8, 0.6, 0.4, 0.2, 0.0))
def test_col_methods(x):
    """
    python -c "from mrspoc.tests.test import test_col_methods as f; f()"
    """

    star = Star()
    star.spots = [Spot(x=x, y=0.00, r=0.2)]

    x0, y0 = star._centroid_analytic()
    x1, y1 = star._centroid_numerical()

    diff = np.sqrt((x0 - x1)**2 + (y0 - y1)**2)
    assert diff < star.r / 100