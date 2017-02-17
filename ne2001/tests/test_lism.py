# Module to run tests on LISM codes
from __future__ import print_function, absolute_import, division, unicode_literals

# TEST_UNICODE_LITERALS

import pytest
import os
import numpy as np

from ne2001 import lism
from ne2001 import io as ne_io


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


def test_DRQ1():
    """ Test DRQ1 density
    """
    ldict = ne_io.read_lism()
    # Float at center
    x,y,z = [ldict[key] for key in ['xldr', 'yldr', 'zldr']]
    ne_DRQ1, FDRQ1, wDRQ1 = lism.neLSB_or_LDRQ1(x,y,z, ldict, 'DRQ1')
    assert np.isclose(ne_DRQ1, ldict['neldr'])
    # Array
    z = np.linspace(-0.3, 3.3, 1000)
    x = np.ones_like(z)*ldict['xldr']
    y = np.ones_like(z)*ldict['yldr']
    ne_DRQ1, FDRQ1, wDRQ1 = lism.neLSB_or_LDRQ1(x,y,z, ldict, 'DRQ1')
    assert np.isclose(ne_DRQ1[0], ldict['neldr'])
    assert np.isclose(ne_DRQ1[-1], 0.)


def test_LSB():
    """ Test LSB density
    """
    ldict = ne_io.read_lism()
    # Float at center
    x,y,z = [ldict[key] for key in ['xlsb', 'ylsb', 'zlsb']]
    ne_LSB, FLSB, wLSB = lism.neLSB_or_LDRQ1(x,y,z, ldict, 'LSB')
    assert np.isclose(ne_LSB, ldict['nelsb'])
    # Array
    z = np.linspace(-0.3, 0.3, 1000)
    x = np.ones_like(z)*ldict['xlsb']
    y = np.ones_like(z)*ldict['ylsb']
    ne_LSB, FLSB, wLSB = lism.neLSB_or_LDRQ1(x,y,z, ldict, 'LSB')
    assert np.isclose(ne_LSB[0], ldict['nelsb'])
    assert np.isclose(ne_LSB[-1], 0.)


def test_LOOPI():
    """ Test LOOPI density
    """
    ldict = ne_io.read_lism()
    # Float at center
    x,y,z = [ldict[key] for key in ['xlpI', 'ylpI', 'zlpI']]
    ne_LOOPI, FLOOPI, wLOOPI = lism.neLOOPI(x,y,z, ldict)
    assert np.isclose(ne_LOOPI, ldict['nelpI'])
    # Float in shell
    x -= ldict['rlpI']
    x -= ldict['drlpI']/2. - 0.001
    ne_LOOPI, FLOOPI, wLOOPI = lism.neLOOPI(x,y,z, ldict)
    assert np.isclose(ne_LOOPI, ldict['dnelpI'])
    # Array
    x = np.linspace(ldict['xlpI']-10*ldict['rlpI'], ldict['xlpI']+10*ldict['rlpI'], 1000)
    y = np.ones_like(x)*ldict['ylpI']
    z = np.ones_like(x)*ldict['zlpI']
    ne_LOOPI, FLOOPI, wLOOPI = lism.neLOOPI(x,y,z, ldict)
    assert np.isclose(np.max(ne_LOOPI), ldict['nelpI'])
    assert np.isclose(ne_LOOPI[-1], 0.)


def test_LHB2():
    """ Test LHB density
    """
    ldict = ne_io.read_lism()
    # Float at center
    x,y,z = [ldict[key] for key in ['xlhb', 'ylhb', 'zlhb']]
    ne_LHB2, FLHB, wLHB = lism.neLHB2(x,y,z, ldict)
    assert np.isclose(ne_LHB2, ldict['nelhb'])
    # Array
    z = np.linspace(-0.3, 0.3, 1000)
    x = np.ones_like(z)*ldict['xlhb']
    y = np.ones_like(z)*ldict['ylhb']
    ne_LHB2, FLHB, wLHB = lism.neLHB2(x,y,z, ldict)
    assert np.isclose(ne_LHB2[0], 0.)
    assert np.isclose(ne_LHB2[-1], ldict['nelhb'])
