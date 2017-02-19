""" Module for density calculations
Mirrors density.NE2001.f
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import pdb

from ne2001 import io as io_ne2001

def ne_GC(x, y, z, FORTRAN_NE2001=True):
    """
    Parameters
    ----------
    x
    y
    z
    FORTRAN_NE2001 : bool, optional
      Use the original algorithm in the FORTRAN code?
      This differs from the last astro-ph posting

    Returns
    -------
    ne_gc : ndarray or float
    F_gc : ndarray or float

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


def ne_spiral_arm(x,y,z, gal_param):
    """
    Parameters
    ----------
    x
    y
    z
    gal_param

    Returns
    -------
    nea : ndarray or float
    whicharm : ndarray or int

    c-----------------------------------------------------------------------
    c  Spiral arms are defined as logarithmic spirals using the
    c    parameterization in Wainscoat et al. 1992, ApJS, 83, 111-146.
    c  But arms are modified selectively at various places to distort them
    c    as needed (08 Aug 2000).
    c  Note that arm numbering follows that of TC93 for the four large arms
    c (after remapping).
    c  The local spiral arm is number 5.
    c  06 Apr 02:   removed TC type modifications of arms 2,3 (fac calculations)
    c  		and replaced with new versions.  Data for these are hard wired.
    """
    from scipy.interpolate import CubicSpline
    if isinstance(x,float):
        x = np.array([x])
        y = np.array([y])
        z = np.array([z])
        flg_float = True
    else:
        flg_float = False

    # see get_parameters for definitions of narm, warm, harm.
    narmsmax=5
    #common/armfactors/
    # .     harm(narmsmax),narm(narmsmax),warm(narmsmax),farm(narmsmax)

    #   parameter(rad=57.29577 95130 823)
    rad = 180/np.pi
    ks=3
    NN=7
    nfine = 1000

    # rr
    rr = np.sqrt(x**2 + y**2)
    adict = io_ne2001.init_spiral_arms()

    #c
    #c Get spiral arm component:  30 do loop finds a coarse minimum distance
    #c from line of sight to arm; 40 do loop finds a fine minimum distance
    #c from line of sight to arm; line 35 ensures that arm limits are not
    #c exceeded; linear interpolation beginning at line 41 finds the
    #c minimum distance from line of sight to arm on a finer scale than gridding
    #c of arms allows (TJL)

    # Init
    whicharm = np.zeros_like(x).astype(int)
    nea = np.zeros_like(x)

    # thxy
    thxy = np.arctan2(-x, y) * rad		#! measured ccw from +y axis
                                        #! (different from tc93 theta)
    neg_th = thxy < 0.
    thxy[neg_th] += 360.
    # Cut on values near the disk
    icutz = np.where(np.abs(z/gal_param['ha']) < 10.)[0]
    if len(icutz) > 0:
        cutx = x[icutz]
        cuty = y[icutz]
        cutz = z[icutz]
        sminmin = 1.e10 * np.ones_like(cutx)
    else:
        # Time to return
        pdb.set_trace()
    # Find closest distance to each arm and then assign arm
    #   We will brute force with a spline
    #   Could get memory intensive
    #   Could refine
    for j in range(adict['narms']):
        #do 50 j=1,narms
        jj = adict['armmap'][j]
        # Crude
        xarm = adict['arm'][j,:adict['kmax'][j],0]
        yarm = adict['arm'][j,:adict['kmax'][j],1]

        xtmp = np.outer(cutx, np.ones_like(xarm))
        ytmp = np.outer(cuty, np.ones_like(yarm))
        xatmp = np.outer(np.ones_like(cutx), xarm)
        yatmp = np.outer(np.ones_like(cuty), yarm)
        dist = (xtmp-xatmp)**2 + (ytmp-yatmp)**2
        kmin = np.argmin(dist,axis=1)

        # Refine with a Spline? -- Requires a loop
        csp_xarm = CubicSpline(np.arange(adict['kmax'][j]),
                               adict['arm'][j,:adict['kmax'][j],0])
        csp_yarm = CubicSpline(np.arange(adict['kmax'][j]),
                               adict['arm'][j,:adict['kmax'][j],1])
        jtmp = np.zeros((len(cutx),nfine))
        for ii,ix in enumerate(cutx):
            j0 = max(0, kmin[ii]-1)
            j1 = min(adict['kmax'][j], kmin[ii]+1)
            jtmp[ii,:] = np.linspace(j0,j1,num=nfine)
        xjarm = csp_xarm(jtmp)
        yjarm = csp_yarm(jtmp)
        xtmp = np.outer(cutx, np.ones(nfine))
        ytmp = np.outer(cuty, np.ones(nfine))
        dist = (xtmp-xjarm)**2 + (ytmp-yjarm)**2
        kmin2 = np.argmin(dist,axis=1)
        min_dist2 = np.amin(dist,axis=1)

        #
        smin=np.sqrt(min_dist2)		 # ! Distance of (x,y,z) from this arm's axis
        # Close enough?
        gd_wa = np.where(smin < (3*gal_param['wa']))[0]
        if len(gd_wa) > 0:
            # ga
            ga = np.exp(-(smin[gd_wa]/(
                gal_param['warm{:d}'.format(jj)]*gal_param['wa']))**2)	#! arm, get the arm weighting factor
            # Galactocentric radial dependence of arms
            tmp_rr = rr[icutz[gd_wa]]
            lg_rr = tmp_rr > gal_param['Aa']
            if np.sum(lg_rr) > 0:
                ga[lg_rr] *= 1. / (np.cosh((tmp_rr[lg_rr]-gal_param['Aa'])/2.0))**2

            # arm3 reweighting:
            if adict['armmap'][j] == 3:
                th3a=320.
                th3b=390.
                th3b=370.
                th3a=290.
                th3b=363.
                th3b=363.
                fac3min=0.0
                test3 = thxy[icutz[gd_wa]]-th3a
                neg_t3 = test3 < 0
                test3[neg_t3] += 360.
                gd_t3 = (0. <= test3) & (test3 < (th3b-th3a))
                if np.sum(gd_t3) > 0:
                    arg = 6.2831853*(thxy[icutz[gd_wa]][gd_t3]-th3a)/(th3b-th3a)
                    fac = (1.+fac3min + (1.-fac3min)*np.cos(arg))/2.
                    fac = fac**4.0
                    # Update ga
                    ga[gd_t3] *= fac

            # arm2 reweighting:
            if adict['armmap'][j] == 2:
                th2a=35.
                th2b=55.
                test2 = thxy[icutz[gd_wa]]-th2a
                fac = 1.
                if False:
                    #    first: as in tc93 (note different definition of theta)
                    gd_t2_first =  (0 <= test2) & (test2 < (th2b-th2a))
                    if np.sum(gd_t2_first) > 0:
                        fac=1.+ test2[gd_t2_first]/(th2b-th2a)
                        fac = 1.		# !**** note turned off
                        #ga=ga*fac
                    gd_t2_first2 =  test2 > (th2b-th2a)
                    if np.sum(gd_t2_first2) > 0:
                        fac = 2.
                        fac = 1.		#!**** note turned off
                        #ga=ga*fac
                # c    second:  weaken the arm in a short range:
                th2a=340.
                th2b=370.
                # c note fix does nothing if fac2min = 1.0
                fac2min=0.1
                neg_t2 = test2 < 0.
                test2[neg_t2] += 360.
                gd_t2_snd = (0. <= test2) & (test2 < (th2b-th2a))
                if np.sum(gd_t2_snd) > 0:
                    arg=6.2831853*(thxy[icutz[gd_wa]][gd_t2_snd]-th2a)/(th2b-th2a)
                    fac = (1.+fac2min + (1.-fac2min)*np.cos(arg))/2.
                    # Update ga
                    ga[gd_t2_snd] *= fac

            # Find the ones which are closer than any previous arm
            iclosest = smin[gd_wa] < sminmin[gd_wa]
            if np.sum(iclosest) > 0:
                whicharm[icutz[gd_wa[iclosest]]] = adict['armmap'][j]

            # Update nea
            nea[icutz[gd_wa]] += gal_param['narm{:d}'.format(jj)] * gal_param['na'] * (
                ga / (np.cosh(cutz[gd_wa]/(gal_param['harm{:d}'.format(jj)]*gal_param['ha']))**2))
    # Finish
    ne_arms_log_mod = nea

    if flg_float:
        nea = float(nea[0])
        whicharm = int(whicharm[0])

    # Return
    return nea, whicharm

    '''
    Farms = 0
    if(whicharm_spiralmodel .eq. 0) then
    whicharm = 0
    else
      whicharm = armmap(whicharm_spiralmodel)	! remap arm number
      Farms = Fa * farm(whicharm)
    endif
    return
    end
    '''


def ne_thin(x,y,z, gal_param, FORTRAN_NE2001=True):
    """
    Parameters
    ----------
    x
    y
    z
    gal_param : dict
    FORTRAN_NE2001 : bool, optional
      Use the original algorithm in the FORTRAN code?
      This differs from the last astro-ph posting

    Returns
    -------
    ne_inn : ndarray or float
    F_inner : ndarray or float

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


def ne_thick(x,y,z, gal_param):
    """
    Parameters
    ----------
    x : ndarray or float
    y : ndarray or float
    z : ndarray or float
    gal_param : dict

    Returns
    -------
    ne_out : ndarray or float
    F_outer : ndarray or float

    c-----------------------------------------------------------------------
    c Thick disk component:
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
    pdb.set_trace()
    ne1 = (gal_param['n1h1']/gal_param['h1'])*g1 / (np.cosh(z/gal_param['h1']))**2
    # Finish
    ne_out = ne1
    F_outer = gal_param['F1']

    return ne_out, F_outer
