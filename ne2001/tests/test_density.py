# Module to run tests on cat_utils
from __future__ import print_function, absolute_import, division, unicode_literals

# TEST_UNICODE_LITERALS

import pytest
import os
import numpy as np

from ne2001 import density


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)

def test_GC():
    """ Test Galactic Center density
    Returns
    -------

    """
    # Original
    # Float
    x,y,z = 0.,0.,0.
    ne_GC, Fgc = density.ne_GC(x,y,z, original=True)
    assert np.isclose(ne_GC,10.)
    # Array
    z = np.linspace(-0.1, 0.1, 100)
    x = np.zeros_like(z)
    y = np.zeros_like(z)
    ne_GC, Fgc = density.ne_GC(x,y,z, original=True)
    assert np.isclose(ne_GC[np.argmin(np.abs(z-0.))], 10.)

    # New (as written)
    ne_GC, Fgc = density.ne_GC(x,y,z, original=False)
    assert np.isclose(ne_GC[np.argmin(np.abs(z+0.02))], 9.9429412976538512)



