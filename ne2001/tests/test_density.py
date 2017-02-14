# Module to run tests on cat_utils
from __future__ import print_function, absolute_import, division, unicode_literals

# TEST_UNICODE_LITERALS

import pytest
import os
import numpy as np

from ne2001 import density
from ne2001 import io as ne_io


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)

def test_GC():
    """ Test Galactic Center density
    """
    # Original
    # Float
    x,y,z = 0.,0.,0.
    ne_GC, Fgc = density.ne_GC(x,y,z, FORTRAN_NE2001=True)
    assert np.isclose(ne_GC,10.)
    # Array
    z = np.linspace(-0.1, 0.1, 100)
    x = np.zeros_like(z)
    y = np.zeros_like(z)
    ne_GC, Fgc = density.ne_GC(x,y,z, FORTRAN_NE2001=True)
    assert np.isclose(ne_GC[np.argmin(np.abs(z-0.))], 10.)
    # New (as written)
    ne_GC, Fgc = density.ne_GC(x,y,z, FORTRAN_NE2001=False)
    assert np.isclose(ne_GC[np.argmin(np.abs(z+0.02))], 9.9429412976538512)

def test_spiral():
    gal_param = ne_io.read_galparam()
    # Array
    x = np.linspace(1.0, 10.0, 100)
    y = np.ones_like(x)
    z = np.linspace(0.1, 3.0, 100)
    ne_spiral, which_arm = density.ne_spiral_arm(x,y,z, gal_param)
    assert np.isclose(ne_spiral[50], 3.0764961696160683e-06, rtol=1e-5)
    assert which_arm[50] == 2
    # Go further out
    x = np.linspace(1.0, 20.0, 100)
    y = np.ones_like(x)
    z = np.ones_like(x) * 0.1
    ne_spiral2, which_arm2 = density.ne_spiral_arm(x,y,z, gal_param)
    assert np.isclose(ne_spiral2[50], 0.011097138978334706, rtol=1e-5)

def test_thin():
    """ Test inner density algorithm
    """
    gal_param = ne_io.read_galparam()
    # Original
    # Float
    x,y,z = 1.,1.,0.1
    ne_inner, F_inner = density.ne_thin(x,y,z, gal_param, FORTRAN_NE2001=True)
    assert np.isclose(ne_inner,0.0091182655991105324, rtol=1e-7)
    # Array
    z = np.linspace(0.1, 1.0, 100)
    x = np.ones_like(z)
    y = np.ones_like(z)
    ne_inner, F_inner = density.ne_thin(x,y,z, gal_param, FORTRAN_NE2001=True)
    assert np.isclose(ne_inner[0],0.0091182655991105324, rtol=1e-7)
    # New (as written)
    ne_inner, Finner = density.ne_thin(x,y,z,gal_param, FORTRAN_NE2001=False)
    assert np.isclose(ne_inner[0], 3.56190777e-02, rtol=1e-7)


def test_thick():
    """ Test outer density algorithm
    """
    gal_param = ne_io.read_galparam()
    # Float
    x,y,z = 1.,10.,1
    ne_out, F_inner = density.ne_thick(x,y,z, gal_param)
    assert np.isclose(ne_out, 0.011686894112477935, rtol=1e-7)
    # Array
    z = np.linspace(0.1, 1.0, 100)
    x = np.ones_like(z)
    y = np.ones_like(z) * 10.
    ne_out, F_inner = density.ne_thick(x,y,z, gal_param)
    assert np.isclose(ne_out[-1], 0.011686894112477935, rtol=1e-7)
