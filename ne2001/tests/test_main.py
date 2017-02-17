# Module to run tests on cat_utils
from __future__ import print_function, absolute_import, division, unicode_literals

# TEST_UNICODE_LITERALS

import pytest
import os
import numpy as np

from ne2001 import main
from ne2001 import io as ne_io


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


def test_density_2001():
    """ Density algorithm
    Returns
    -------

    """
    # Setup
    gal_param = ne_io.read_galparam()
    ldict = ne_io.read_lism()
    mflags = dict(wg1=True, wg2=True, wga=True, wggc=True, wglism=True, ORIG_NE2001=True)
    # Array
    x = np.linspace(1.0, 10.0, 1000)
    y = np.ones_like(x)
    z = np.linspace(0.1, 3.0, 1000)
    density_dict = main.density_2001(x,y,z,gal_param,ldict,mflags)
    for key in [u'thick', u'GC', u'thin', u'LISM', u'arms']:
        assert key in density_dict.keys()

