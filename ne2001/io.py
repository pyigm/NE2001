""" Module for I/O
"""
from __future__ import absolute_import, division, print_function

import json

import ne2001
data_path = ne2001.__path__[0]+'/data/'

def read_galparam(ifile='gal_param.json'):
    """ Read Galaxy parameters
    Parameters
    ----------
    ifile : str, optional

    Returns
    -------
    gal_param : dict

    """
    with open(data_path+ifile, 'rt') as fh:
        galparam_dict = json.load(fh)
    # Return
    return galparam_dict

def read_gc(ifile='ne_gc.json'):
    """ Read Galactic Center parameters
    Returns
    -------
    gc_dict : dict
      dict of parameters

    """
    # Read
    with open(data_path+ifile, 'rt') as fh:
        gc_dict = json.load(fh)
    # Return
    return gc_dict


