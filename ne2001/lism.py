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



def ne_LISM(x,y,z, ldict):
    """
    Parameters
    ----------
    x
    y
    z
    ldict

    Returns
    -------

    """
                          #  wldr, wlhb,wlsb,wloopI)

    #real nelsbxyz, nelhbxyz, neldrq1xyz, neloopIxyz
    #real FLDRQ1r, FLSBr, FLHBr, FLOOPIr		! 'r' for returned value
    #integer wLDR, wLSB, wLHB, wLOOPI

    #character*48 path
    #character*120 nelisminp

    neldrq1xyz, FLDRQ1r, wLDR = neLSB_or_LDRQ1(x,y,z,ldict,'DRQ1')	#! low density region in Q1
    nelsbxyz, FLSBr, wLSB  = neLSB_or_LDRQ1(x,y,z,ldict,'LSB')      #! Local Super Bubble
    nelhbxyz, FLHBr, wLHB = neLHB2(x,y,z,ldict)		#! Local Hot Bubble
    neloopIxyz, FLOOPIr, wLOOPI = neLOOPI(x,y,z,ldict) #	! Loop I


    #c weight the terms so that the LHB term overrides the other
    #c terms (we want the density to be low in the LHB, lower than
    #c in the other terms.

    neLISM =   (1-wLHB) * (
                   (1-wLOOPI) * (wLSB*nelsbxyz + (1-wLSB)*neldrq1xyz)
                   + wLOOPI * neloopIxyz ) + wLHB  * nelhbxyz

    FLISM = (1-wLHB) * (
                   (1-wLOOPI) * (wLSB*FLSBr + (1-wLSB)*FLDRQ1r)
                   + wLOOPI * FLOOPIr ) +   wLHB  * FLHBr

    #c return the maximum weight of any of the terms for
    #c combining with additional terms external to this routine.

    wLISM = np.maximum(wLOOPI, np.maximum(wLDR, np.maximum(wLSB, wLHB)))

    return neLISM, FLISM, wLISM


def neLSB_or_LDRQ1(x,y,z, ldict, region): #	! Local Super Bubble or Low Density Region in Q1
    """
    Parameters
    ----------
    x
    y
    z
    ldict : dict
      LISM parameters
    region : str
      'LSB' : Local super bubble
      'DRQ1' : Low density region in Q1

    Returns
    -------
    ne : ndarray or float
    F : ndarray or float
    w : ndarray or int

    c input:
    c 	x,y,z = coordinates w.r.t. Galaxy as in TC93, CL00
    c output:
    c	neLSB = electron density in local hot bubble that
    c	        is modeled as an ellisoidal trough.
    c	FLSB = fluctuation parameter
    c	wLSB  = weight of LSB component used to combine
    c		with other components of electron density.
    c		wLSB =  1  at and inside the annular ridge
    c		     <  1  outside the annular ridge
    c	             -> 0  far outside the annular ridge
    c	e.g. total electron density would be evaluated as
    c            ne = (1-wLSB)*ne_other + neLSB
    """
    #real aa,bb,cc			! scales of ellipsoidal ridge
    #real netrough			! ne of annulus, trough
    #real Ftrough			! fluctuation parameters
    #real theta 			! position angle of major axis,
    #				!    measured from x axis
    #				!    (x axis points toward l=90)
    if isinstance(x,float):
        x = np.array([x])
        y = np.array([y])
        z = np.array([z])
        flg_float = True
    else:
        flg_float = False

    radian = 180./np.pi #57.29577951)

    if region == 'LSB':
        akey, bkey, ckey, thkey, nekey, Fkey = 'alsb', 'blsb', 'clsb', 'thetalsb', 'nelsb', 'Flsb',
        xkey, ykey, zkey = [ss+'lsb' for ss in ['x','y','z']]
    elif region == 'DRQ1':
        akey, bkey, ckey, thkey, nekey, Fkey = 'aldr', 'bldr', 'cldr', 'thetaldr', 'neldr', 'Fldr',
        xkey, ykey, zkey = [ss+'ldr' for ss in ['x','y','z']]
    aa=ldict[akey]
    bb=ldict[bkey]
    cc=ldict[ckey]
    theta=ldict[thkey]
    netrough=ldict[nekey]
    Ftrough=ldict[Fkey]

    s = np.sin(theta/radian)
    c = np.cos(theta/radian)
    ap = (c/aa)**2 + (s/bb)**2
    bp = (s/aa)**2 + (c/bb)**2
    cp = 1./cc**2
    dp =  2.*c*s*(1./aa**2 - 1./bb**2)

    ne = np.zeros_like(x)
    w = np.zeros_like(x).astype(int)
    F = np.zeros_like(x)

    # Geometry
    q = ap*(x-ldict[xkey])**2 + bp*(y-ldict[ykey])**2 + (
        cp*(z-ldict[zkey])**2 + (x-ldict[xkey])*(y-ldict[ykey])*dp)
    inside = q <= 1.
    ne[inside] = netrough
    F[inside] = Ftrough
    w[inside] = 1
    # Float?
    if flg_float:
        ne = float(ne[0])
        F= float(F[0])
        w= int(w[0])

    return ne, F, w


def neLHB2(x,y,z, ldict):
    """ Local Hot Bubble

    Parameters
    ----------
    x
    y
    z
    ldict

    Returns
    -------
    ne_LHB2 :

    c LHB modeled as a cylinder
    c the cylinder slants in the y direction vs. z as described by parameter yzslope
    c the cylinder cross-sectional size in the 'a' direction (major axis)
    c       varies with z, tending to zero at its smallest z point.
        implicit none
        real x,y,z,FLHBr
        integer wLHB
    c input:
    c 	x,y,z = coordinates w.r.t. Galaxy as in TC93, CL00
    c output:
    c	neLHB2 = electron density in local hot bubble that
    c	        is modeled as an ellisoidal trough.
    c	FLHB = fluctuation parameter
    c	wLHB  = weight of LBH component used to combine
    c		with other components of electron density.
    c		wLHB =  1  at and inside the annular ridge
    c		     <  1  outside the annular ridge
    c	             -> 0  far outside the annular ridge
    c	e.g. total electron density would be evaluated as
    c            ne = (1-wLHB)*ne_other + neLHB2
    """
    if isinstance(x,float):
        x = np.array([x])
        y = np.array([y])
        z = np.array([z])
        flg_float = True
    else:
        flg_float = False
    #real aa,bb,cc			! scales of ellipsoidal ridge
    #real netrough			! ne of annulus, trough
    #real Ftrough			! fluctuation parameters
    #real xlhb, ylhb, zlhb		! center of ellipsoid
    #real theta 			! slant angle in yz plane of cylinder
    #				!    measured from z axis
    #real qxy, qz
    #    real radian
    #    parameter(radian = 57.29577951)
    radian = 180./np.pi

    #real yzslope
    #real yaxis

    # Init
    bb=ldict['blhb']
    cc=ldict['clhb']
    theta=ldict['thetalhb']
    netrough=ldict['nelhb']
    Ftrough=ldict['Flhb']
    yzslope = np.tan(theta/radian)

    ne_LHB2 = np.zeros_like(x)
    wLHB = np.zeros_like(x).astype(int)
    FLHBr = np.zeros_like(x)
    yaxis = ldict['ylhb'] + yzslope*z
    # cylinder has cross sectional area = constant for z>0
    # area -> 0 for z<0 by letting aa->0 linearly for z<0:
    # (0.001 = 1 pc is to avoid divide by zero)
    aa = np.ones_like(x) * ldict['alhb']
    neg_z = (z <= 0.) & (z >= ldict['zlhb']-ldict['clhb'])
    aa[neg_z] = 0.001 + (ldict['alhb']-0.001)*(
        1. - (1./(ldict['zlhb']-ldict['clhb']))*z[neg_z])
    qxy =  ( (x-ldict['xlhb'])/aa )**2 + ( (y-yaxis)/bb )**2
    qz =  np.abs(z-ldict['zlhb'])/cc
    # Inside?
    inside = (qxy <= 1.0) & (qz <= 1.0)
    ne_LHB2[inside] = netrough
    FLHBr[inside] = Ftrough
    wLHB[inside] = 1
    # Recast?
    if flg_float:
        ne_LHB2 = float(ne_LHB2[0])
        FLHBr = float(FLHBr[0])
        wLHB = int(wLHB[0])
    # Return
    return ne_LHB2, FLHBr, wLHB


def neLOOPI(x,y,z, ldict): #	! Loop I
    """
    Parameters
    ----------
    x
    y
    z
    ldict : dict
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
    in_a1 = r < a1
    ne_LOOPI[gd_z & in_a1] = ldict['nelpI']
    FLOOPI[gd_z & in_a1] = ldict['FlpI']
    wLOOPI[gd_z & in_a1] = 1
    # Boundary Shell
    in_shell = (r >= a1) & (r < a2)
    ne_LOOPI[gd_z & in_shell] = ldict['dnelpI']
    FLOOPI[gd_z & in_shell] = ldict['dFlpI']
    wLOOPI[gd_z & in_shell] = 1

    # Float?
    if flg_float:
        ne_LOOPI = float(ne_LOOPI[0])
        FLOOPI = float(FLOOPI[0])
        wLOOPI = int(wLOOPI[0])

    return ne_LOOPI, FLOOPI, wLOOPI
