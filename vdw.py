#! /usr/bin/env python3

"""vdw.py - van der Waals fluid functions."""

import numpy
from numpy import cosh, sinh
import scipy.optimize


# Standard functions ---------------------------------------------------------

def pr(Tr, vr):
    """Returns the reduced pressure given the reduced temperature and volume."""
    return 8*Tr/(3*vr-1) - 3/vr**2


# Maxwell construction functions ---------------------------------------------
#
# Ref: Lekner, AJP (1982).

def _params(ds):
    """Saturation parameters."""
    y = ds/2
    f = (y*cosh(y) - sinh(y))/(sinh(y)*cosh(y)-y)
    c = cosh(y)
    g = 1 + 2*c*f + f**2

    return f, c, g

def Trsat(ds):
    """Returns the reduced saturation temperature.
    ds - the entropy change during the phase change; use ds().
    """
    f, c, g = _params(ds)
    return 27/4*f*(c+f)/g**2

def prsat(ds):
    """Returns the reduced saturation pressure.
    ds - the entropy change during the phase change; use ds().
    """
    f, c, g = _params(ds)
    return 27*f**2*(1-f**2)/g**2

def vrsat(ds):
    """Returns the saturation volumes.
    ds - the entropy change during the phase change; use ds().
    """

    # Get the reduced saturation temperature and pressure
    Tr = Trsat(ds)
    pr = prsat(ds)
    
    # The van der Waals equation in cubic form is
    # a0 + a1 vr + a2 vr^2 + a3 vr^3 = 0
    # where the coefficients are as follows:
    a = [1, -((8*Tr/pr)+1)/3, 3/pr, -1/pr]
    
    # Get the bounding volumes by finding the roots of the cubic equation
    return numpy.roots(a)[::-2]    
        
def ds(Tr):
    """Returns the non-dimensional entropy change during the phase transition
    at the given reduced temperature."""
    return scipy.optimize.brentq(lambda ds: Trsat(ds)-Tr, 1e-5, 700)   

def maxwell(Tr, vr, pr):
    """Applies Maxwell construction to a van der Waals isotherm.

    Tr - the reduced temperature of the isotherm
    vr - the reduced volume array
    pr - the reduced pressure array
    """

    # Get the saturation pressure and bounding volumes
    ds_ = ds(Tr)
    prsat_ = prsat(ds_)
    v1, v2 = vrsat(ds_)

    # Flatten out the pressure between the bounding volumes
    i1, i2 = numpy.searchsorted(vr, [v1, v2])
    pr[i1:i2] = prsat_

    return pr
