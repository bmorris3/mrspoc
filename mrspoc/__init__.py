from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

# Packages may add whatever they like to this file, but
# should keep this content at the top.
# ----------------------------------------------------------------------------
from ._astropy_init import *
# ----------------------------------------------------------------------------

if not _ASTROPY_SETUP_:
    from .star import *
    from .tgas import *
    from .gaia import *
    from .sun import *
