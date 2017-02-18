""" Top-level algorithms
c density.NE2001.f
c final version of NE2001
c returns densities, F parameters and weights of the various components
c mods:
c 28 July 02:
c   put in 'if' statements in subroutine density_2001 so that
c	function calls are not done if particular weights
c	(wg1, wg2, etc.) are set to zero in gal01.inp
c	This is (a) cleaner and (b) much more efficient if
c	the clump or void component is flagged out.
"""
from __future__ import print_function, absolute_import, division, unicode_literals

from ne2001 import density
from ne2001 import lism


def density_2001(x,y,z, gal_param, ldict, mflags, **kwargs):
    """
    Parameters
    ----------
    x
    y
    z
    gal_param
    ldict
    kwargs

    Returns
    -------
    out_dict : dict
      Contains all the info

    c----------------------------------------------------------------------------
    c  Returns seven components of the free electron density of the
    c  interstellar medium at Galactic location (x,y,z).
    c  Calling arguments:
    c  input:
    c	x, y, z = galactocentric location (kpc)
    c       Right-handed coordinate system
    c       x is in l=90 direction
    c       y is in l=180 direction
    c       The sun is at (x,y,z) = (0,R0,0)
    c  output:
    c    electron densities in cm^{-3}:
    c	ne1:	outer, thick disk
    c	ne2:	inner, thin disk (annular in form)
    c	nea:	spiral arms
    c	negc:   galactic center component
    c       nelism: local ISM component
    c       necN:   contribution from discrete 'clumps'
    c       nevN:   contribution from voids
    c    fluctuation parameters (one for each ne component):
    c       F1, F2, Fa, Fgc, Flism, FcN, FvN
    c    flags:
    c       whicharm: which of the 5 arms x,y,z is in (0 for interarm region)
    c          wlism: 1 if x,y,z is in any of the four LISM components
    c           wLDR: 1 if in LDR, 0 if not
    c           wLHB: 1 if in LHB, 0 if not
    c           wLSB: 1 if in LSB, 0 if not
    c         wLOOPI: 1 if in LoopI, 0 if not
    c       (nb: nelism is calculated according to LHB:LOOPI:LSB:LDR)
    c       hitclump: clump number that x,y,z is in (0 if none)
    c        hitvoid: void number that x,y,z is in (0 if none)
    c 25 May 2002
    c based on routines from TC93 and test routines from 1999-2002 by JMC.
    c----------------------------------------------------------------------------
    """
    #integer wg1, wg2, wga, wggc, wglism, wgcN, wgvN
    #common /modelflags/ wg1, wg2, wga, wggc, wglism, wgcN, wgvN

    out_dict = {}

    # Thick disk (aka outer)
    if mflags['wg1']:
        ne_thick, F_thick = density.ne_thick(x,y,z, gal_param)
        out_dict['thick'] = {}
        out_dict['thick']['ne'] = ne_thick
        out_dict['thick']['F'] = F_thick
    # Thin disk (aka inner)
    if mflags['wg2']:
        ne_thin, F_thin = density.ne_thin(x,y,z, gal_param, FORTRAN_NE2001=mflags['ORIG_NE2001'])
        out_dict['thin'] = {}
        out_dict['thin']['ne'] = ne_thin
        out_dict['thin']['F'] = F_thin
    if mflags['wga']:
        ne_arm, whicharm = density.ne_spiral_arm(x,y,z, gal_param)
        out_dict['arms'] = {}
        out_dict['arms']['ne'] = ne_arm
        out_dict['arms']['which'] = whicharm
    if mflags['wggc']:
        ne_gc, F_gc = density.ne_GC(x,y,z, FORTRAN_NE2001=mflags['ORIG_NE2001'])
        out_dict['GC'] = {}
        out_dict['GC']['ne'] = ne_gc
        out_dict['GC']['F'] = F_gc
    if mflags['wglism']:
        ne_lism, FLISM, wLISM = lism.ne_LISM(x,y,z, ldict)
        out_dict['LISM'] = {}
        out_dict['LISM']['ne'] = ne_lism
        out_dict['LISM']['F'] = FLISM
        out_dict['LISM']['w'] = wLISM
        #,wlism,wldr,wlhb,wlsb,wloopI)
    #if(wgcN .eq. 1) call neclumpN(x,y,z,necN,FcN,hitclump)
    #if(wgvN .eq. 1) call nevoidN(x,y,z,nevN,FvN,hitvoid,wvoid)

    # Output dict

    return out_dict


def NE2001_dens():
    """
    Returns
    -------
    c  calls density_2001 to calculate thermal electron density
    c  Returns seven components of the free electron density of the
    c  interstellar medium at Galactic location (x,y,z).
    c  Calling arguments:
    c  input:
    c	x, y, z = galactocentric location (kpc)
    c       Right-handed coordinate system
    c       x is in l=90 direction
    c       y is in l=180 direction
    c       The sun is at (x,y,z) = (0,R0,0)
    c  output:
    c    electron densities in cm^{-3}:
    c	ne1:	outer, thick disk
    c	ne2:	inner, thin disk (annular in form)
    c	nea:	spiral arms
    c	negc:   galactic center component
    c       nelism: local ISM component
    c       necN:   contribution from discrete 'clumps'
    c       nevN:   contribution from voids
    c    fluctuation parameters (one for each ne component):
    c       F1, F2, Fa, Fgc, Flism, FcN, FvN
    c    flags:
    c       whicharm: which of the 5 arms x,y,z is in (0 for interarm region)
    c          wlism: 1 if x,y,z is in any of the four LISM components
    c           wLDR: 1 if in LDR, 0 if not
    c           wLHB: 1 if in LHB, 0 if not
    c           wLSB: 1 if in LSB, 0 if not
    c         wLOOPI: 1 if in LoopI, 0 if not
    c       (nb: nelism is calculated according to LHB:LOOPI:LSB:LDR)
    c       hitclump: clump number that x,y,z is in (0 if none)
    """
    pass
"""
    implicit none
    integer nargs
    character*80 inbuf
    integer ndir
    character*1 limit

    real vperp
        data vperp/100./

c New stuff added by me
    real x,y,z

    integer whicharm

c Large scale components:

    real ne1, ne2, nea
    real F1val, F2val, Faval

c Galactic center:

    real negc, Fgc

c LISM:
    real nelism, Flism
        integer wlism, wLDR, wLHB, wLSB, wLOOPI

c clumps:
        real necN, FcN
    integer hitclump

c voids:
    real nevN, FvN
    integer hitvoid, wvoid
c functions:

    integer iargc
    external iargc

    nargs = iargc()
    if(nargs .ge. 1) then
       call getarg(1, inbuf)
       read(inbuf, *) x
       call getarg(2, inbuf)
       read(inbuf, *) y
       call getarg(3, inbuf)
       read(inbuf, *) z
    else
       write(*,*) 'Usage: NE2001_dens x y z'
           write(*,*) '       x (kpc)'
           write(*,*) '       y (kpc)'
           write(*,*) '       z (kpc)'
       stop
        endif

    call density_2001(x,y,z,
     .        ne1,ne2,nea,negc,nelism,necN,nevN,
     .        F1val, F2val, Faval, Fgc, Flism, FcN, FvN,
     .        whicharm, wlism, wldr, wlhb, wlsb, wloopI,
     .        hitclump, hitvoid, wvoid)
    write(*,"(e15.5, e15.5, e15.5, e15.5, e15.5, e15.5, e15.5)") ne1,
     .        ne2, nea, negc, nelism, necN, nevN

    stop
    end
"""
