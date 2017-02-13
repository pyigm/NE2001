""" Module for density calculations
Mirrors density.NE2001.f
"""

import numpy as np
import pdb

import ne2001
from ne2001 import io as io_ne2001

def ne_GC(x, y, z, original=True):
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
        if original:
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



    return
