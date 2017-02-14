""" Module for density calculations
Mirrors neLISM.NE2001.f
c routines to calculate the electron density for the
c Local Interstellar Medium
c
c JMC 26 August-11 Sep. 2000
c     25 October 2001: modified to change weighting scheme
c                      so that the ranking is LHB: LSB: LDR
c                      (LHB overrides LSB and LDR; LSB overrides LDR)
c     16 November 2001: added Loop I component with weighting scheme
c		        LHB:LOOPI:LSB:LDR
c		        LHB   overides everything,
c			LOOPI overrides LSB and LDR
c			LSB   overrides LDR
c			LISM  overrides general Galaxy
c     20 November 2001: The LOOPI component is truncated below z=0
c
c after discussions with Shami Chatterjee
c the sizes, locations and densities of the LISM components
c are based largely on work by Toscano et al. 1999
c and Bhat et al. 1999
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import pdb

import ne2001
from ne2001 import io as io_ne2001



def neLOOPI(x,y,z): #	! Loop I
    """
    Parameters
    ----------
    x
    y
    z
    FLOOPI
    wLOOPI

    Returns
    -------

    c component is a spheroid truncated for z<0.
        implicit none
        real x,y,z,FLOOPI
        integer wLOOPI
    c input:
    c 	x,y,z = coordinates w.r.t. Galaxy as in TC93, CL00
    c output:
    c	neLOOPI = electron density in LOOP I that
    c	        is modeled as an ellisoidal trough
    c		with an enhanced shell
    c	FLOOPI = fluctuation parameter
    c	wLOOPI  = weight of LOOP I component used to combine
    c		with other components of electron density.
    c		wLOOPI =  1  at and inside the annular ridge
    c		       <  1  outside the annular ridge
    """
    ldict = io_ne2001.read_lism()
    #real r
    #real a1, a2
    #logical first
    #data first /.true./
    #save
    if isinstance(x,float):
        x = np.array([x])
        y = np.array([y])
        z = np.array([z])
        flg_float = True
    else:
        flg_float = False

    # Init
    a1 = ldict['rlpI']
    a2 = a1+ ldict['drlpI']
    ne_LOOPI = np.zeros_like(x)
    FLOOPI = np.zeros_like(x)
    wLOOPI = np.zeros_like(x).astype(int)

    # Cut on z
    gd_z = z >= 0.

    # Cut on r
    r = np.sqrt( (x-ldict['xlpI'])**2 + (y-ldict['ylpI'])**2 + (z-ldict['zlpI'])**2)

    # Inside volume
    in_a1 = r < ldict['a1']
    ne_LOOPI[gd_z & in_a1] = ldict['nelpI']
    FLOOPI[gd_z & in_a1] = ldict['FlpI']
    wLOOPI[gd_z & in_a1] = 1
    # Boundary Shell
    in_shell = (r >= ldict['a1']) & (r < ldict['a2'])
    ne_LOOPI[gd_z & in_shell] = ldict['dnelpI']
    FLOOPI[gd_z & in_shell] = ldict['dFlpI']
    wLOOPI[gd_z & in_shell] = 1

    return ne_LOOPI, FLOOPI, wLOOPI
