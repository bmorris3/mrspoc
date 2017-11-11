.. include:: references.txt

.. _details:

**********************
Implementation details
**********************

``mrspoc`` uses two methods to compute astrometric stellar centroid offsets due
to starspots. The *analytic method* computes the approximate missing flux
due to starspots, normalized by the total flux. The *numerical method* creates a
grid of pixels which simulate the star and its spots. We summarize the two
methods in more detal below.

Contents
========

* :ref:`details-analytic`
* :ref:`details-numerical`

.. _details-analytic:

Analytic method
---------------

The analytic solution is used by default for spots smaller than 0.1 stellar
radii. Here we'll outline the algorithm for the analytic solution.

The total flux of the entire surface of the limb-darkened, unspotted star is
the integral of the limb-darkening law :math:``


.. _details-numerical:

Numerical method
----------------

The numerical solution
