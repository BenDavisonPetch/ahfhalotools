"""
A collection of methods useful for dealing with AHF files, including automatic
snapshot to redshift mapping from filenames, and file truncation.
"""

import os
import numpy as np
from .objects import *
from .analysis import *

def truncateFiles(fileNameBase,snaps,zs,outputFileNameBase,numHalos,
                  haloFileExt = ".AHF_halos",profileFileExt = ".AHF_profiles",
                  mtreeidxExt = ".AHF_mtree_idx", mtreeExt = ".AHF_mtree"):
    """
    Truncates sets of AHF data files such as to shorten load times.

    Parameters
    ----------
    fileNameBase : str
        The format string of the base file name of the files to be truncated
        (including path)
        e.g. "./GadgetX/GadgetX-NewMDCLUSTER_0001.snap_{snap:0=3d}.z{z:.3f}"
        For specifics on formatting see notes section below.
    snaps : list of int
        snapshot numbers of files to load
    zs : list of float
        redshifts of snapshots corresponding to snapshot numbers in snaps
    outputFileNameBase : str
        The format string of the output file. Formatted the same as
        fileNameBase
    numHalos : int
        The number of halos to truncate the files to

    Other Parameters
    ----------------
    haloFileExt : str, default = ".AHF_halos"
        File extension to use for AHF_halos files
    profileFileExt : str, default = ".AHF_profiles"
        File extension to use for AHF_profiles files
    mtreeidxExt : str, default = ".AHF_mtree_idx"
        File extension to use for AHF_mtree_idx files
    mtreeExt : str, default = ".AHF_mtree"
        File extension to use for AHF_mtree files

    Notes
    -----
    File names are generated using the python str.format() method, specifically
    as: fileNameBase.format(snap = snapNum, z = z)

    So in the example:
    "./GadgetX/GadgetX-NewMDCLUSTER_0001.snap_{snap:0=3d}.z{z:.3f}",
    {snap:0=3d} will be replaced with the snapshot number, including leading
    zeroes to make snapshot number 3 characters, and {z:.3f} will be replaced
    with the redshift to 3 decimal places.

    So if the file being read is for snapshot 97 with redshift 0.986, the file
    base name (without extensions added) will be:
    "./GadgetX/GadgetX-NewMDCLUSTER_0001.snap_097.z0.986"
    """
    if len(zs) != len(snaps):
        raise ValueError("snaps and zs are not of same length!")

    numFiles = len(zs)

    #loop through snapshots
    for i in range(numFiles):
        try:
            fileName = fileNameBase.format(snap=snaps[i],z=zs[i])
        except KeyError as err:
            raise ValueError("fileNameBase is not properly formatted! String"\
                            " contains a key that is not z or snap") from err
        try:
            outName = outputFileNameBase.format(snap = snaps[i],z = zs[i])
        except KeyError as err:
            raise ValueError("outName is not properly formatted! String"\
                            " contains a key that is not z or snap") from err

        print("About to truncate {0}:".format(fileName))
        #truncate profiles and halos files

        snap = Snapshot(snaps[i],zs[i],fileName+haloFileExt,
                        fileName+profileFileExt, haloLimit = numHalos)
        snap.writeFiles(outName)
        print("    Profiles and Halos done")

        #truncate mtree_idx file
        try:
            mtreeidxrows = np.genfromtxt(fileName + mtreeidxExt)
            np.savetxt(outName + mtreeidxExt,mtreeidxrows[:numHalos,:],fmt='%d')
            print("    mtree_idx done")
        except IOError:
            print("    WARNING: mtree_idx file not found")

        #truncate mtree file
        try:
            mtreerows = np.genfromtxt(fileName + mtreeExt)
            #loop through rows to check how many rows need to be written to new file
            haloCounter = 0
            lineCounter = 0
            while haloCounter < numHalos:
                haloID, haloPart, numProg = mtreerows[lineCounter,:]
                lineCounter += int(numProg)+1
                haloCounter += 1
            #lineCounter is now equal to the number of lines that we care about
            np.savetxt(outName+mtreeExt,mtreerows[0:lineCounter,:],fmt='%d')
            print("    mtree done")
        except IOError:
            print("    WARNING: mtree file not found")
        print("Completed truncation: {0:.1f}%".format((i+1)/numFiles * 100))

def getSnapNumToZMapGX(directory=""):
    """
    Searches directory specified and returns a map of snapshot number to
    redshift.

    Assumes file names are of format GadgetX-NewMDCLUSTER_0001.snap_119.z0.221...
    Directory must only contain files with names of this format.

    Parameters
    ----------
    directory : str
        Directory to search. If not specified, defaults to current directory.

    Returns
    -------
    snapNumToZMap : dict { int snapNum : float z }

    See also
    --------
    getSnapNumToZMapGiz : Gets snapshot number to redshift map for GIZMO files
    getMusZs : Gets list of redshifts of Gadget-MUSIC files in directory
    """
    files = os.listdir(directory)
    #assuming file name is of format GadgetX-NewMDCLUSTER_0001.snap_119.z0.221...
    #returns dictionary of {snapNo:z}
    return {int(file[31:34]):float(file[36:41]) for file in files}

def getSnapNumToZMapGiz(directory=""):
    """
    Searches directory specified and returns a map of snapshot number to
    redshift.

    Assumes file names are of format GIZMO-NewMDCLUSTER_0001.snap_119.z0.221...
    Directory must only contain files with names of this format.

    Parameters
    ----------
    directory : str
        Directory to search. If not specified, defaults to current directory.

    Returns
    -------
    snapNumToZMap : dict { int snapNum : float z }

    See also
    --------
    getSnapNumToZMapGX : Gets snapshot number to redshift map for GadgetX files
    getMusZs : Gets list of redshifts of Gadget-MUSIC files in directory
    """
    files = os.listdir(directory)
    #assuming file name is of format GIZMO-NewMDCLUSTER_0001.snap_119.z0.221...
    #returns dictionary of {snapNo:z}
    return {int(file[29:32]):float(file[34:39]) for file in files}

def getMusZs(directory=""):
    """
    Searches directory specified and returns a list of redshifts for GadgetMUSIC
    snapshots.

    Assumes file names are of format GadgetMUSIC-NewMDCLUSTER_0001.z0.000...
    Directory must only contain files with names of this format.

    Parameters
    ----------
    directory : str
        Directory to search. If not specified, defaults to current directory.

    Returns
    -------
    zs : list of floats
        Redshifts in descending order

    See also
    --------
    getSnapNumToZMapGX : Gets snapshot number to redshift map for GadgetX files
    getSnapNumToZMapGiz : Gets snapshot number to redshift map for GIZMO files
    """
    files = os.listdir(directory)
    #assuming file name is of format GadgetMUSIC-NewMDCLUSTER_0001.z0.000...
    zs = [float(file[31:36]) for file in files]
    #remove duplicates
    zs = list(set(zs))
    #sort list such that in descending order
    zs.sort(reverse=True)
    return zs
