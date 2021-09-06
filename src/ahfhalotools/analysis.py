"""
Various methods useful for data processing.
"""

import numpy as np
from astropy.cosmology import FlatLambdaCDM
from astropy.cosmology import WMAP9
import astropy.constants as apconst

def removeNegatives(x,y):
    """
    Takes arrays x and y and removes any entries where x is negative.
    x and y must be of same length.

    Parameters
    ----------
    x : array of numbers
    y : array of numbers

    Returns
    -------
    x : array of numbers
    y : array of numbers
    """
    if len(x) != len(y):
        raise NameError("removeNegatives(): Arrays given are not of same length!")

    output = ([],[])

    for i in range(len(x)):
        if x[i] >= 0:
            output[0].append(x[i])
            output[1].append(y[i])

    return np.array(output[0]) , np.array(output[1])

def UtoT(U):
    """
    Takes internal energy in (km/s)^2, and returns temperature in kelvin.
    Assumes mass of gas is entirely due to protons and assumes ideal gas.

    Parameters
    ----------
    U : float or np.array of dtype float
        Internal energy in (km/s)^2

    Returns
    -------
    T : float or np.array of dtype float
        Temperature in K
    """
    mass_p = 1.6726219e-27
    k_boltz = 1.38064852e-23
    temp = (2/3 * mass_p / k_boltz) * U * 1e6
    return temp

def getAsFuncOfTimeAtRadius(snaps,haloNum,val_func,radius):
    """
    Used to get a profile data quantity, interpolated at a specific redius,
    as a function of redshift.

    .. deprecated:: 0.0.0
        getAsFuncOfTimeAtRadius should not be used, as it does not accurately
        track halos through time. It has been replaced by
        ahfhalotools.objects.Cluster.funcOfZProfileData, and
        ahfhalotools.objects.Cluster.funcOfAgeProfileData

    Parameters
    ----------
    snaps : list of Snapshot instances
        list of Snapshot instances in ascending order of snapshot number
    haloNum : int
        index of halo within Snapshot instance to use
    val_func : lambda
        lambda that takes a Halo object and returns the array of values in
        question as a function of radius.
        e.g:
        lambda halo : halo.gaslocDensities()
    radius : float
        Radius to interpolate at in kpc/h

    Returns
    -------
    values : list of floats

    See also
    --------
    ahfhalotools.objects.Cluster.funcOfZProfileData
    ahfhalotools.objects.Cluster.funcOfAgeProfileData
    """
    vals = []
    for snap in snaps:
        halo = snap.halos[haloNum]
        interp = np.interp(radius, halo.radii(), val_func(halo))
        vals.append(interp)
    return np.array(vals)

def tfromz(z):
    """
    Converts from redshift to age

    Parameters
    ----------
    z : float
        Redshift

    Returns
    -------
    age : float
        Age in Gyr

    Notes
    -----
    Uses astropy.cosmology.WMAP9.H(0) as current hubble constant
    Also uses Omega_M and Omega_B from GadgetX parameters
    (Om0=0.307, Ob0=0.048)
    """
    flcdm = FlatLambdaCDM(H0=WMAP9.H(0), Om0=0.307, Ob0=0.048)
    return flcdm.age(z).value
