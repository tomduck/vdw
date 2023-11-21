#! /usr/bin/env python3

"""vdw.py - van der Waals fluid functions."""

import numpy
from numpy import cosh, sinh, exp
import scipy.optimize

# Equations numbers are from Johnston (2014)

# van der Waals state functions ----------------------------------------------

def pr(Tr, vr):
    """Returns the reduced pressure.

    Tr - the reduced temperature
    vr - the reduced volume
    """
    return 8*Tr/(3*vr-1) - 3/vr**2                       # Eq. 38

def Tr(pr, vr):
    """Returns the reduced temperature.
    of state.

    pr - the reduced pressure
    vr - the reduced volume
    """
    return (pr + 3/vr**2)*(3*vr-1)/8                     # Eq. 38 (rearranged)

def s(Tr, vr):
    """Returns the non-dimensional entropy.

    Tr - the reduced temperature
    vr - the reduced volume
    """
    return numpy.log(Tr**1.5*(3*vr-1)/2)                 # Eq. 58

# Maxwell construction functions ---------------------------------------------

f = lambda y: (y*cosh(y) - sinh(y))/(sinh(y)*cosh(y)-y)  # Eq. 78b
g = lambda y: 1 + 2*f(y)*cosh(y) + f(y)**2               # Eq. 78c

def Trsat(ds):
    """Returns the reduced saturation temperature.
    ds - the non-dimensional entropy difference between pure gas and
         pure liquid phases; use ds().
    """
    y = ds/2
    return 27*f(y)*(f(y)+cosh(y))/(4*g(y)**2)            # Eq. 79b

def prsat(ds):
    """Returns the reduced saturation pressure.
    ds - the non-dimensional entropy difference between pure gas and
         pure liquid phases; use ds().
    """
    y = ds/2
    return 27*f(y)**2*(1-f(y)**2)/g(y)**2                # Eq. 79c

def vrsat(ds):
    """Returns the reduced saturation volumes.
    ds - the non-dimensional entropy difference between pure gas and
         pure liquid phases; use ds().
    """
    y = ds/2
    return (1 + exp(-y)/f(y))/3, (1 + exp(y)/f(y))/3     # Eq. 79f, 79e

def ds(*, Tr=None, pr=None):
    """Returns the non-dimensional entropy difference between pure gas and
    pure liquid phases.

    Exactly one of the following parameters should be defined:

    Tr - the reduced temperature during the phase transition
    pr - the reduced pressure during the phase transition
    """

    # Make sure only one of Tr and pr is defined
    assert Tr is None or pr is None
    assert not (pr is None and Tr is None)

    # Define the function to optimize
    f = (lambda ds: Trsat(ds) - Tr) if Tr is not None else \
        (lambda ds: prsat(ds) - pr)

    # Return the optimized entropy change
    return scipy.optimize.brentq(f, 1e-5, 700)

def maxwell(Tr, vr, pr):
    """Applies Maxwell construction to a van der Waals isoline
    (either an isotherm or isobar).

    Tr - the reduced temperature value or array
    vr - the reduced volume array
    pr - the reduced pressure value or array

    The Maxwell construction will be applied to whichever one of Tr and pr
    is an array.  The other must be a value.

    Returns the array to which the Maxwell construction has been applied.
    If the Maxwell construction can't be applied, then the array is passed
    through.
    """

    assert numpy.isscalar(Tr) != numpy.isscalar(pr)

    # If the Maxwell construction doesn't apply, return the appropriate isoline
    if numpy.isscalar(Tr) and Tr >= 1:
        return pr
    if numpy.isscalar(pr) and pr >= 1:
        return Tr

    # Get the normalized entropy change
    ds_ = ds(Tr=Tr) if numpy.isscalar(Tr) else ds(pr=pr)

    # Get indices for the bounding volumes
    i1, i2 = numpy.searchsorted(vr, vrsat(ds_))

    # Flatten the isoline
    if numpy.isscalar(Tr):
        pr[i1:i2] = prsat(ds_)
    else:
        Tr[i1:i2] = Trsat(ds_)

    return pr if numpy.isscalar(Tr) else Tr
