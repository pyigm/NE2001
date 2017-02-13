""" Module for density calculations
Mirrors density.NE2001.f
"""

import numpy as np
import pdb

import ne2001
from ne2001 import io as io_ne2001

def ne_GC(x, y, z, FORTRAN_NE2001=True):
    """
    c-----------------------------------------------------------------------
    c     Determines the contribution of the Galactic center to the free
    c     electron density of the interstellar medium at Galactic location
    c     (x,y,z).  Combine with `fluctuation` parameter to obtain the
    c     scattering measure.
    c
    c     NOTE: This is for the hyperstrong scattering region in the
    c     Galactic center.  It is distinct from the inner Galaxy
    c     (component 2) of the TC93 model.
    c
    c     Origin of coordinate system is at Galactic center; the Sun is at
    c     (x,y,z) = (0,R0,0), x is in l=90 direction
    c
    c     Based on Section 4.3 of Lazio & Cordes (1998, ApJ, 505, 715)
    c
    c Input:
    c REAL X - location in Galaxy [kpc]
    c REAL Y - location in Galaxy [kpc]
    c REAL Z - location in Galaxy [kpc]
    c
    c COMMON:
    c REAL NEGC0 - nominal central density
    c
    c PARAMETERS:
    c REAL RGC - radial scale length of Galactic center density enhancement
    c REAL HGC - z scale height of Galactic center density enhancement
    c
    c Output:
    c REAL NE_GC - Galactic center free electron density contribution [cm^-3]
    c-----------------------------------------------------------------------
    c
    """
    #implicit none
    #real x, y, z, F_gc

    #real rgc, hgc
    #real xgc, ygc, zgc
#     parameter (xgc=-0.010, ygc=0., zgc=-0.020)
#     parameter (rgc=0.145)
#     parameter (hgc=0.026)

    #real rr, zz
    #real arg

    #real negc0, Fgc0
    #common /gcparms/ negc0, Fgc0

    #real n1h1,h1,A1,F1,n2,h2,A2,F2,na,ha,wa,Aa,Fa
    #common /galparams/ n1h1,h1,A1,F1,n2,h2,A2,F2,
    # .                na,ha,wa,Aa,Fa
    # Coordinate output as float or array
    if isinstance(x,float):
        ne_gc = 0.
        F_gc = 0.
    elif isinstance(x,np.ndarray):
        ne_gc = np.zeros_like(x)
        F_gc = np.zeros_like(x)
    else:
        raise IOError("Bad input for x")

    # Load
    gc_dict = io_ne2001.read_gc()
    xgc, ygc, zgc = [gc_dict['centroid'][key] for key in ['xgc','ygc','zgc']] # For convenience

    # GALACTOCENTRIC RADIUS
    rr = np.sqrt( (x-xgc)**2 + (y-ygc)**2)
    gd_rgc = rr < gc_dict['rgc'] # truncate at 1/e point

    # Z-HEIGHT.
    zz = np.abs(z-zgc)
    gd_zz = zz < gc_dict['hgc']

    # Close enough?
    arg = (rr/gc_dict['rgc'])**2 + (zz/gc_dict['hgc'])**2
    gd_arg = arg <= 1.
    # Finish
    all_gd = gd_rgc & gd_zz & gd_arg
    if np.sum(all_gd) > 0:
        if FORTRAN_NE2001:
            if isinstance(x,float):
                ne_gc = gc_dict['negc0']
                F_gc = gc_dict['Fgc0']
            else:
                ne_gc[all_gd] = gc_dict['negc0']
                F_gc[all_gd] = gc_dict['Fgc0']
        else:
            if isinstance(x,float):
                ne_gc = gc_dict['negc0'] * np.exp(-1*arg)
                F_gc = gc_dict['Fgc0']
            else:
                ne_gc[all_gd] = gc_dict['negc0'] * np.exp(-1*arg[all_gd])
                F_gc[all_gd] = gc_dict['Fgc0']
    # Return
    return ne_gc, F_gc


def ne_inner(x,y,z, gal_param, FORTRAN_NE2001=True):
    """
    Parameters
    ----------
    x
    y
    z
    gal_param : dict

    Returns
    -------
    c-----------------------------------------------------------------------
    c Thin disk (inner Galaxy) component:
    c (referred to as 'Galactic center component' in circa TC93 density.f)
    """
    # Init
    if isinstance(x,float):
        g2 = 0.
    elif isinstance(x,np.ndarray):
        g2 = np.zeros_like(x)
    else:
        raise IOError("Bad input for x")
    # r
    rr=np.sqrt(x**2 + y**2)
    if FORTRAN_NE2001:
        rrarg=((rr-gal_param['A2'])/1.8)**2
    else:
        rrarg=((rr-gal_param['A2'])/gal_param['A2'])**2
    # g2
    gd_g2 = rrarg < 10.0
    if np.sum(gd_g2) > 0:
        if isinstance(x,float):
            g2 = np.exp(-1*rrarg)
        else:
            g2[gd_g2] = np.exp(-1*rrarg[gd_g2])
    # Calculate
    ne2 = gal_param['n2']*g2/(np.cosh(z/gal_param['h2'])**2)
    ne_inn = ne2
    F_inner = gal_param['F2']
    # Return
    return ne_inn, F_inner


def ne_outer(x,y,z, gal_param):
    """
    c-----------------------------------------------------------------------
    c Thick disk component:
    implicit none
    real x,y,z, F_outer
    """
    #parameter(pihalf=3.14159 26535 89793 23846 264/2.0)
# 	g1=sech2(rr/A1)/sech2(8.5/A1)		! TC93 function
    if isinstance(x,float):
        g1 = 0.
    elif isinstance(x,np.ndarray):
        g1 = np.zeros_like(x)
    else:
        raise IOError("Bad input for x")
    # Sun cos
    suncos = np.cos(np.pi*gal_param['rsun']/gal_param['A1']/2.)
    # rr
    rr=np.sqrt(x**2 + y**2)
    gd_rr = rr <= gal_param['A1']
    if np.sum(gd_rr) > 0:
        if isinstance(x,float):
            g1 = np.cos(np.pi*rr/gal_param['A1']/2.)/suncos
        else:
            g1[gd_rr] = np.cos(np.pi*rr[gd_rr]/gal_param['A1']/2.)/suncos
    # ne1
    ne1 = (gal_param['n1h1']/gal_param['h1'])*g1 / np.cosh(z/gal_param['h1'])**2
    # Finish
    ne_out = ne1
    F_outer = gal_param['F1']

    return ne_out, F_outer
