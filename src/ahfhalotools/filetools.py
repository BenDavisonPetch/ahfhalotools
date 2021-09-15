"""
A collection of methods useful for dealing with AHF files, including automatic
snapshot to redshift mapping from filenames, reading multiple clusters at once,
and file truncation.
"""

import os
import numpy as np
from . import objects as obj
from . import analysis
import string

def truncateFiles(fileNameBase,snaps,zs,outputFileNameBase,numHalos,
                  haloFileExt = ".AHF_halos",profileFileExt = ".AHF_profiles",
                  mtreeidxExt = ".AHF_mtree_idx", mtreeExt = ".AHF_mtree", skipmtree = False):
    """
    Truncates sets of AHF data files for a single cluster such as to shorten
    load times.

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
    skipmtree : bool, default = False
        Whether to skip reading the .AHF_mtree file. Useful when the mtree
        file is in the non-supported SUSSING-2013 format

    See Also
    --------
    truncateClusters : truncates files from multiple clusters at once

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

        snap = obj.Snapshot(snaps[i],zs[i],fileName+haloFileExt,
                        fileName+profileFileExt, haloLimit = numHalos)
        snap.writeFiles(outName)
        print("    Profiles and Halos done")

        #truncate mtree_idx file
        try:
            mtreeidxrows = np.genfromtxt(fileName + mtreeidxExt, max_rows = numHalos)

            #check if mtreeidxrows is 1D and make it 2D if so
            if mtreeidxrows.ndim == 1:
                mtreeidxrows = np.array([mtreeidxrows])
            np.savetxt(outName + mtreeidxExt,mtreeidxrows,fmt='%d')
            print("    mtree_idx done")
        except IOError:
            print("    WARNING: mtree_idx file not found")

        #truncate mtree file
        if not skipmtree:
            try:
                mtreerows = np.genfromtxt(fileName + mtreeExt)
                if len(mtreerows[0]) == 3:
                    #non sussing format
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
                else:
                    #wrong format
                    print("    WARNING: mtree is not correct format!")
            except IOError:
                print("    WARNING: mtree file not found")
            except ValueError:
                print("    WARNING: mtree is not of correct format! It may be in SUSSING-2013 format, which is not supported")
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
    zDict = dict()
    for file in files:
        try:
            zDict[int(file[31:34])] = float(file[36:41])
        except ValueError:
            pass
    return zDict

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
    zDict = dict()
    for file in files:
        try:
            zDict[int(file[29:32])] = float(file[34:39])
        except ValueError:
            pass
    return zDict

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
    zs = set()
    for file in files:
        if file[30] == "z":
            try:
                zs.add( float(file[31:36]) )
            except ValueError:
                pass
    zs = list(zs)
    #sort list such that in descending order
    zs.sort(reverse=True)
    return zs

def getZfromFileName(fileName, suppress = False):
    """
    Lifts a redshift from a file name.

    Assumes file names are formatted as "...zX.XXX..." where X.XXX is the
    redshift.

    Parameters
    ----------
    fileName : str
    suppress : bool, default = False
        Whether or not to suppress error messages. If no single redshift is
        found and suppress is True, will return None

    Returns
    -------
    z : float
    """
    #first we find where all of the "z" characters are in the name
    zIndexes = [pos for pos, char in enumerate(fileName) if char == 'z']
    zs = list()
    for zIndex in zIndexes:
        try:
            zs.append( float(fileName[zIndex+1:zIndex+6]) )
        except ValueError:
            #bit of string after the z is not a number
            pass
    if len(zs) == 0:
        msg = "Could not find a redshift in file name {0}\nPlease check the name is formatted as '...zX.XXX...'".format(fileName)
        if not suppress:
            raise ValueError(msg)
        else:
            return
    elif len(zs) > 1:
        msg = "Found multiple redshifts in file name {0}".format(fileName)
        if not suppress:
            raise ValueError(msg)
        else:
            return
    return zs[0]

def getZs(directory, num = -1):
    """
    Gets the redshifts of the snapshot files in a folder in descending order

    Parameters
    ----------
    directory : str
        The directory to search
    num : int, optional
        Specifies how many redshifts to return. Will return the first num
        redshifts starting from redshift 0
        If not specified will return all redshifts found

    Returns
    -------
    zs : list of float
        The redshifts of the files in the directory in descending order.
    """
    #get list of all redshifts in folder
    files = os.listdir(directory)
    zs = list(set([getZfromFileName(file, suppress = True) for file in files if getZfromFileName(file, suppress = True) != None]))
    zs.sort(reverse=True)
    #we want redshifts in descending order as snapNos should be ascending
    #only want num redshifts if specified
    if num != -1 and not num > len(zs):
        zs = zs[-num:]
    return zs

class SafeFormatter(string.Formatter):
    """
    A string formatter used to preserve unused keys.

    Taken from https://stackoverflow.com/a/34033230
    """
    def vformat(self, format_string, args, kwargs):
        args_len = len(args)  # for checking IndexError
        tokens = []
        for (lit, name, spec, conv) in self.parse(format_string):
            # re-escape braces that parse() unescaped
            lit = lit.replace('{', '{{').replace('}', '}}')
            # only lit is non-None at the end of the string
            if name is None:
                tokens.append(lit)
            else:
                # but conv and spec are None if unused
                conv = '!' + conv if conv else ''
                spec = ':' + spec if spec else ''
                # name includes indexing ([blah]) and attributes (.blah)
                # so get just the first part
                fp = name.split('[')[0].split('.')[0]
                # treat as normal if fp is empty (an implicit
                # positional arg), a digit (an explicit positional
                # arg) or if it is in kwargs
                if not fp or fp.isdigit() or fp in kwargs:
                    tokens.extend([lit, '{', name, conv, spec, '}'])
                # otherwise escape the braces
                else:
                    tokens.extend([lit, '{{', name, conv, spec, '}}'])
        format_string = ''.join(tokens)  # put the string back together
        # finally call the default formatter
        return string.Formatter.vformat(self, format_string, args, kwargs)

def loadClusters(clusterNums, snapNums, simName, clusterFolderFmt = "NewMDCLUSTER_{clusterNum:0=4d}/",
                 directory = "", fileBaseFmt = "{simName}-NewMDCLUSTER_{clusterNum:0=4d}.snap_{snap:0=3d}.z{z:.3f}",
                 profileExt=".AHF_profiles", haloExt=".AHF_halos",
                 mtreeidxExt=".AHF_mtree_idx", mtreeExt=".AHF_mtree", haloLimit=np.inf,
                 printProgress = False, skipmtree = False):
    """
    Loads multiple cluster simulations into memory.

    Parameters
    ----------
    clusterNums : list of int
        The list of cluster numbers to load
    snapNums : list of int
        The list of snapshot numbers to load from each cluster in _ascending_
        order. Should not contain any breaks and should include the snapshot
        with z = 0.000
    simName : str
        The name of the simulation code to be inserted into file names
        (e.g. "GadgetX", "GadgetMUSIC")
    directory : str, optional
        The directory the cluster folders are in. If not specified will search
        from current directory.
    clusterFolderFmt : str, optional
        The string format used to determine the names of the clusters' folders
        Defaults to "NewMDCLUSTER_{clusterNum:0=4d}/"
        See notes section for formatting of file strings
    fileBaseFmt : str, optional
        The format of the file name for the AHF files
        Defaults to "{simName}-NewMDCLUSTER_{clusterNum:0=4d}.snap
            _{snap:0=3d}.z{z:.3f}"

    Other Parameters
    ----------------
    profileExt : str, optional
        Defaults to ".AHF_profiles"
    haloExt : str, optional
        Defaults to ".AHF_halos"
    mtreeidxExt : str, optional
        Defaults to ".AHF_mtree_idx"
    mtreeExt : str, optional
        Defaults to ".AHF_mtree"
    haloLimit : int, optional
        Specifies the maximum number of halos to load into memory
    skipmtree : bool, optional
        Defaults to False
        If True, will skip trying to read mtree files, avoiding the warning
        message spam if no mtree files are present.
    printProgress : bool, optional
        Defaults to False
        If True, will print updates on loading progress

    Returns
    -------
    clusters : list of ahfhalotools.objects.Cluster instances
        The list of clusters corresponding to each of the clusters specified
        in clusterNums

    See Also
    --------
    truncateClusters : truncates files for a set of clusters
    ahfhalotools.objects.Cluster : Stores information about a cluster

    Notes
    -----
    The parameters clusterFolderFmt, dir, and fileBaseFmt are dynamically
    formatted using the str.format method to produce the correct file names.
    The appropriate keys are:
        clusterNum : the number of the cluster
        simName : the name of the simulation code
        snap : the snapshot number
        z : the redshift

    For example, if the following set of parameters is passed into the function:
        directory = "/home/{simName}/"
        clusterFolderFmt = "NewMDCLUSTER_{clusterNum:0=4d}/"
        fileBaseFmt = "{simName}-NewMDCLUSTER_{clusterNum:0=4d}.sna
                       p_{snap:0=3d}.z{z:.3f}"
        simName = "GadgetX"
    With default file name extensions, when loading cluster 34, snapshot number
    99 with z = 0.9 (for example), the following files will be opened:

    > /home/GadgetX/NewMDCLUSTER_0034/GadgetX-NewMDCLUSTER_0034.snap_099.z0.900.AHF_halos
    > /home/GadgetX/NewMDCLUSTER_0034/GadgetX-NewMDCLUSTER_0034.snap_099.z0.900.AHF_profiles
    > /home/GadgetX/NewMDCLUSTER_0034/GadgetX-NewMDCLUSTER_0034.snap_099.z0.900.AHF_mtree_idx
    > /home/GadgetX/NewMDCLUSTER_0034/GadgetX-NewMDCLUSTER_0034.snap_099.z0.900.AHF_mtree
    """
    if not directory.endswith("/"):
        directory += "/"
    if not clusterFolderFmt.endswith("/"):
        clusterFolderFmt += "/"

    fmt = SafeFormatter()
    numClusters = len(clusterNums)

    clusters = []
    for i, clusterNum in enumerate(clusterNums):

        #input relevant fields for file names
        inputDir = fmt.format(directory + clusterFolderFmt,clusterNum=clusterNum,simName=simName)
        fileBaseName = fmt.format(directory + clusterFolderFmt + fileBaseFmt,clusterNum=clusterNum,simName=simName)
        #check if cluster directory exists, if not the cluster may have been
        #skipped during truncation, so we skip loading cluster and print a
        #warning
        if not os.path.isdir(inputDir):
            print("WARNING: Non-existent cluster data:\n{0}".format(inputDir))
            continue
        #get redshifts in folder
        zs = getZs(inputDir, num = len(snapNums))

        #check if folder is empty
        if len(zs) == 0:
            print("WARNING: Non-existent cluster data:\n{0}".format(inputDir))
            continue

        cluster = obj.Cluster(fileBaseName, snapNums, zs, profileExt = profileExt,
                          haloExt = haloExt, mtreeidxExt = mtreeidxExt,
                          mtreeExt = mtreeExt, haloLimit = haloLimit,
                          clusterNum = clusterNum, simName = simName,
                          skipmtree = skipmtree)
        clusters.append(cluster)
        if printProgress:
            print("Loaded {0} cluster {1:0=4d}: {2:.1f}%".format(simName,clusterNum,(i+1)/numClusters * 100))
    if printProgress:
        print("\nLoading for {0} complete.\n".format(simName))
    return clusters

def truncateClusters(clusterNums, snapNums, simName, haloLimit, outputDir,
                     clusterFolderFmt = "NewMDCLUSTER_{clusterNum:0=4d}/",
                     directory = "", fileBaseFmt = "{simName}-NewMDCLUSTER_{clusterNum:0=4d}.snap_{snap:0=3d}.z{z:.3f}",
                     profileExt=".AHF_profiles", haloExt=".AHF_halos",
                     mtreeidxExt=".AHF_mtree_idx", mtreeExt=".AHF_mtree",
                     skipmtree = False):
    """
    Truncates files from multiple cluster folders at once

    Parameters
    ----------
    clusterNums : list of int
        The list of cluster numbers to load
    snapNums : list of int
        The list of snapshot numbers to load from each cluster in _ascending_
        order. Should not contain any breaks and should include the snapshot
        with z = 0.000
    simName : str
        The name of the simulation code to be inserted into file names
        (e.g. "GadgetX", "GadgetMUSIC")
    haloLimit : int
        The number of halos to truncate to
    outputDir : str
        The directory to put the truncated files in. The file hierarchy will
        look the same as in the original data set.
    directory : str, optional
        The directory the cluster folders are in. If not specified will search
        from current directory.
    clusterFolderFmt : str, optional
        The string format used to determine the names of the clusters' folders
        Defaults to "NewMDCLUSTER_{clusterNum:0=4d}/"
        See notes section for formatting of file strings
    fileBaseFmt : str, optional
        The format of the file name for the AHF files
        Defaults to "{simName}-NewMDCLUSTER_{clusterNum:0=4d}.snap
            _{snap:0=3d}.z{z:.3f}"
    skipmtree : bool, default = False
        Whether to skip reading the .AHF_mtree file. Useful if the mtree files
        are in the non-supported SUSSING-2013 format

    Other Parameters
    ----------------
    profileExt : str, optional
        Defaults to ".AHF_profiles"
    haloExt : str, optional
        Defaults to ".AHF_halos"
    mtreeidxExt : str, optional
        Defaults to ".AHF_mtree_idx"
    mtreeExt : str, optional
        Defaults to ".AHF_mtree"

    See Also
    --------
    loadClusters : loads files for a set of clusters
    truncateFiles : truncates files for a single cluster
    ahfhalotools.objects.Cluster : Stores information about a cluster

    Notes
    -----
    The parameters clusterFolderFmt, dir, and fileBaseFmt are dynamically
    formatted using the str.format method to produce the correct file names.
    The appropriate keys are:
        clusterNum : the number of the cluster
        simName : the name of the simulation code
        snap : the snapshot number
        z : the redshift

    For example, if the following set of parameters is passed into the function:
        directory = "/home/{simName}/"
        clusterFolderFmt = "NewMDCLUSTER_{clusterNum:0=4d}/"
        fileBaseFmt = "{simName}-NewMDCLUSTER_{clusterNum:0=4d}.sna
                       p_{snap:0=3d}.z{z:.3f}"
        simName = "GadgetX"
    With default file name extensions, when loading cluster 34, snapshot number
    99 with z = 0.9 (for example), the following files will be opened:

    > /home/GadgetX/NewMDCLUSTER_0034/GadgetX-NewMDCLUSTER_0034.snap_099.z0.900.AHF_halos
    > /home/GadgetX/NewMDCLUSTER_0034/GadgetX-NewMDCLUSTER_0034.snap_099.z0.900.AHF_profiles
    > /home/GadgetX/NewMDCLUSTER_0034/GadgetX-NewMDCLUSTER_0034.snap_099.z0.900.AHF_mtree_idx
    > /home/GadgetX/NewMDCLUSTER_0034/GadgetX-NewMDCLUSTER_0034.snap_099.z0.900.AHF_mtree
    """
    if not directory.endswith("/"):
        directory += "/"
    if not clusterFolderFmt.endswith("/"):
        clusterFolderFmt += "/"
    if not outputDir.endswith("/"):
        outputDir += "/"

    numClusters = len(clusterNums)
    #using safe formatter to retain unused keys
    fmt = SafeFormatter()

    for i, clusterNum in enumerate(clusterNums):
        #format file base name with cluster number and simulation name
        inputDir = fmt.format(directory + clusterFolderFmt,clusterNum=clusterNum,simName=simName)
        fileBaseName = fmt.format(directory + clusterFolderFmt + fileBaseFmt,clusterNum=clusterNum,simName=simName)
        outputFileNameBase = fmt.format(outputDir + clusterFolderFmt + fileBaseFmt,clusterNum=clusterNum,simName=simName)

        outputClusterDir = fmt.format(outputDir + clusterFolderFmt, clusterNum=clusterNum, simName=simName)
        outputDirFmtted = fmt.format(outputDir, clusterNum=clusterNum, simName=simName)

        zs = getZs(inputDir, num = len(snapNums))

        #check if output files exist already, if not create them
        if not os.path.isdir(outputDirFmtted):
            os.mkdir(outputDirFmtted)
        if not os.path.isdir(outputClusterDir):
            os.mkdir(outputClusterDir)

        #truncate files for cluster
        try:
            truncateFiles(fileBaseName,snapNums,zs,outputFileNameBase,haloLimit,
                              haloFileExt = haloExt,profileFileExt = profileExt,
                              mtreeidxExt = mtreeidxExt, mtreeExt = mtreeExt,
                              skipmtree = skipmtree)
            print("\n===CLUSTER {0} DONE===    Overall Completion: {1:.1f}%\n".format(clusterNum, (i+1)/numClusters * 100))
        except OSError:
            print("\n\n===WARNING===\nAn error was encountered while reading the"\
                  " data from Cluster {0}. Cluster will be skipped. Be sure to "\
                  "check your data!\n\n".format(clusterNum))
