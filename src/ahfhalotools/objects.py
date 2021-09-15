"""
A collection of objects used to store and analyse AHF data.
"""

import numpy as np
from . import analysis
import warnings

#number of columns in profile table
PROFILE_COLUMNS = 29

'''========================================================='''
'''                       HALO CLASS                        '''
'''========================================================='''

class Halo:
    """
    Stores information about a single halo at a single point in time.

    Attributes
    ----------
    halodata : dict (str : ?)
        Dictionary of halo data quantities
    ID : int
        ID of halo
    z : float
        Redshift of halo
    age : float
        Age of halo in Gyr
    sigV : float
        3D velocity dispersion of halo in km/s
    Rvir : float
        Virial radius in kpc/h
    Xc : float
        x position of halo in kpc/h (comoving)
    Yc : float
        y position of halo in kpc/h (comoving)
    Zc : float
        z position of halo in kpc/h (comoving)
    pos : np.array ( float )
        3D position vector of halo in kpc/h (comoving)
        Equal to [Xc,Yc,Zc]
    vel : np.array ( float )
        3D velocity vector of halo in km/s
    mbp_offset : float
        Offset between most bound particle and halo centre in kpc/h
    com_offset : float
        Offset between centre-of-mass and halo centre in kpc/h
    Ekin : float
        Kinetic energy in Msol/h (km/s)^2
    Epot : float
        Potential energy in Msol/h (km/s)^2
    Qvir : float
        Virial ratio
        Equal to -Ekin/Epot
    profiles : 2D array of float
        Raw profile data
    rawhalodata : array of dtype float
        Raw halodata row

    Notes
    -----
    To fully initialise halo, one must manually add each profile row using
    Halo.addProfile()

    All named halo data attributes are also present in the halodata dictionary.

    To retreive valid halodata keys, one can call Halo.halodata.keys() on a
    Halo instance.

    To retrieve valid profile data keys, one can call Halo.getProfileQuantityDict()
    on a Halo instance.

    Raw profile data can be accessed using Halo.profiles, which contains a
    list of rows corresponding to rows in the .AHF_profiles file. This can be
    used to retrive profile data that does not have a defined function in the
    Halo class (e.g. to get vcirc, which is column 6 in .AHF_profiles, one can
    call Halo.profiles[:,5])
    """
    def __init__(self, halodata, z):
        """
        Initialises halo instance from row of halodata

        Parameters
        ----------
        halodata : array of float
            The row of halo data from the .AHF_halos file
        z : float
            The redshift of the snapshot the halo belongs to

        Notes
        -----
        To fully initialise halo, one must manually add each row of profile data
        using Halo.addProfile()
        """
        halodata_columns = ['ID', 'hostHalo', 'numSubStruct', 'Mvir', 'npart',
            'Xc', 'Yc', 'Zc', 'VXc', 'VYc', 'VZc', 'Rvir', 'Rmax', 'r2',
            'mbp_offset', 'com_offset', 'Vmax', 'v_esc', 'sigV', 'lambda',
            'lambdaE', 'Lx', 'Ly', 'Lz', 'b', 'c', 'Eax', 'Eay', 'Eaz', 'Ebx',
            'Eby', 'Ebz', 'Ecx', 'Ecy', 'Ecz', 'ovdens', 'nbins', 'fMhires',
            'Ekin', 'Epot', 'SurfP', 'Phi0', 'cNFW', 'n_gas', 'M_gas',
            'lambda_gas', 'lambdaE_gas', 'Lx_gas', 'Ly_gas', 'Lz_gas', 'b_gas',
            'c_gas', 'Eax_gas', 'Eay_gas', 'Eaz_gas', 'Ebx_gas', 'Eby_gas',
            'Ebz_gas', 'Ecx_gas', 'Ecy_gas', 'Ecz_gas', 'Ekin_gas', 'Epot_gas',
            'n_star', 'M_star', 'lambda_star', 'lambdaE_star', 'Lx_star',
            'Ly_star', 'Lz_star', 'b_star', 'c_star', 'Eax_star', 'Eay_star',
            'Eaz_star', 'Ebx_star', 'Eby_star', 'Ebz_star', 'Ecx_star',
            'Ecy_star', 'Ecz_star', 'Ekin_star', 'Epot_star', 'mean_z_gas', 'mean_z_star']
        self.halodata = dict()
        for i in range(len(halodata_columns)):
            colName = halodata_columns[i]
            self.halodata[colName] = halodata[i]

        self.ID = halodata[0]
        self.rawhalodata = halodata
        self.profiles = np.empty((0,PROFILE_COLUMNS))
        self.z = z
        self.halodata['z'] = z
        #age is in Gyr
        self.age = analysis.tfromz(z)
        self.halodata['age'] = self.age
        #sigV in km/sec
        self.sigV = halodata[18]
        #Rvir,Xc,Yc,Zc in kpc/h
        self.Rvir = halodata[11]
        self.Xc = halodata[5]
        self.Yc = halodata[6]
        self.Zc = halodata[7]
        self.pos = np.array([self.Xc,self.Yc,self.Zc])
        self.halodata['pos'] = self.pos

        #velocity in km/s
        self.halodata['vel'] = np.array([self.halodata['VXc'], self.halodata['VYc'], self.halodata['VZc']])
        self.vel = self.halodata['vel']
        #mbp_offset and com_offset in kpc/h
        self.mbp_offset = halodata[14]
        self.com_offset = halodata[15]

        #energy in M_odot/h (km/sec)^2
        self.Ekin = halodata[38]
        self.Epot = halodata[39]

        #Mass of halo in M_odot/h
        self.Mvir = halodata[3]

        #virial ratio dimensionless
        self.Qvir = -self.Ekin/self.Epot
        self.halodata['Qvir'] = self.Qvir

    def addProfile(self, row):
        """
        Used to add profile rows to the halo object while initialising

        Parameters
        ----------
        row : array of dtype float
            The data row from .AHF_profiles

        Notes
        -----
        Profiles should be added in *increasing* order of |radius|, input will
        not be sorted.
        """
        self.profiles = np.append(self.profiles,[row],axis=0)

    def getProfileQuantityDict(self):
        """
        Returns a dictionary of valid profile quantity keys, mapped to bound
        methods for the relevant quantity.

        Returns
        -------
        quantityDict : dict { quantity (str) : bound method }
        """
        quantityDict = {'M_in_r' : self.M_in_r,
                        'M_in_shell' : self.M_in_shell,
                        'gasM_in_r' : self.gasM_in_r,
                        'gasM_in_shell' : self.gasM_in_shell,
                        'volumes' : self.volumes,
                        'encDens' : self.encDensities,
                        'locDens' : self.locDensities,
                        'gasEncDens' : self.gasencDensities,
                        'gasLocDens' : self.gaslocDensities,
                        'u_gas' : self.intEnergies,
                        'T' : self.temps,
                        'starM_in_r' : self.starM_in_r,
                        'starM_in_shell' : self.starM_in_shell,
                        'starEncDens' : self.starencDensities,
                        'starLocDens' : self.starlocDensities}
        return quantityDict

    def radii(self):
        """
        Returns an array containing all of the radii of the profiles.
        Units of kpc/h

        Returns
        -------
        radii : array of dtype float
        """
        return self.profiles[:,0]

    def M_in_r(self):
        """
        Returns an array containing the total mass enclosed as a function of
        profile radius.
        Units of Msol/h

        Returns
        -------
        M_in_r : array of dtype float
        """
        return self.profiles[:,2]

    def M_in_shell(self):
        """
        Returns an array containing the total masses contained only within the
        (hollow) shell that each profile makes with the last profile
        Units of Msol/h

        Returns
        -------
        M_in_shell : array of dtype float
        """
        encMass = self.M_in_r()
        shellMass = encMass - np.append([0],encMass[:-1])
        return shellMass

    def gasM_in_r(self):
        """
        Returns an array containing the gas mass enclosed as a function of
        profile radius.
        Units of Msol/h

        Returns
        -------
        gasM_in_r : array of dtype float
        """
        return self.profiles[:,24]

    def gasM_in_shell(self):
        """
        Returns an array containing the gas masses contained only within the
        (hollow) shell that each profile makes with the last profile
        Units of Msol/h

        Returns
        -------
        gasM_in_shell : array of dtype float
        """
        encMass = self.gasM_in_r()
        shellMass = encMass - np.append([0],encMass[:-1])
        return shellMass

    def starM_in_r(self):
        """
        Returns an array containing the stellar mass enclosed as a function of
        profile radius.
        Units of Msol/h

        Returns
        -------
        starM_in_r : array of dtype float
        """
        return self.profiles[:,25]

    def starM_in_shell(self):
        """
        Returns an array containing the stellar masses contained only within the
        (hollow) shell that each profile makes with the last profile
        Units of Msol/h

        Returns
        -------
        starM_in_shell : array of dtype float
        """
        encMass = self.starM_in_r()
        shellMass = encMass - np.append([0],encMass[:-1])
        return shellMass

    def volumes(self):
        """
        Returns an array containing the volumes of each spherical profile
        Units of kpc^3 / h^3

        Returns
        -------
        volumes : array of dtype float
        """
        return abs((4/3)*np.pi*(self.radii()**3))

    def shellvolumes(self):
        """
        Returns an array containing the volumes of each profile *shell* (ie
        the volume between the current radius and the previous radius)

        Returns
        -------
        shellvolumes : array of dtype float
        """
        vols = self.volumes()
        shellvols = vols - np.append([0], vols[:-1])
        return shellvols

    def encDensities(self):
        """
        Returns an array containing the enclosed total densities for each
        profile.
        Units of Msol kpc^-3 h^2

        Returns
        -------
        dens : array of dtype float
        """
        encMass = self.M_in_r()
        volumes = self.volumes()
        return encMass/volumes

    def locDensities(self):
        """
        Returns an array containing the local total densities for each
        profile.
        Units of Msol kpc^-3 h^2

        Returns
        -------
        dens : array of dtype float
        """
        shellMass = self.M_in_shell()
        volumes = self.shellvolumes()
        return shellMass/volumes
        #return self.profiles[:,4]

    def gasencDensities(self):
        """
        Returns an array containing the enclosed gas densities for each
        profile.
        Units of Msol kpc^-3 h^2

        Returns
        -------
        dens : array of dtype float
        """
        encMass = self.gasM_in_r()
        volumes = self.volumes()
        return encMass/volumes

    def gaslocDensities(self):
        """
        Returns an array containing the local gas densities for each
        profile.
        Units of Msol kpc^-3 h^2

        Returns
        -------
        dens : array of dtype float
        """
        shellMass = self.gasM_in_shell()
        volumes = self.shellvolumes()
        return shellMass/volumes
        #return self.profiles[:,4]

    def starencDensities(self):
        """
        Returns an array containing the enclosed stellar densities for each
        profile.
        Units of Msol kpc^-3 h^2

        Returns
        -------
        dens : array of dtype float
        """
        encMass = self.starM_in_r()
        volumes = self.volumes()
        return encMass/volumes

    def starlocDensities(self):
        """
        Returns an array containing the local stellar densities for each
        profile.
        Units of Msol kpc^-3 h^2

        Returns
        -------
        dens : array of dtype float
        """
        shellMass = self.starM_in_shell()
        volumes = self.shellvolumes()
        return shellMass/volumes

    def intEnergies(self):
        """
        Returns an array containing the internal gas energy for each profile.
        Units : ?(km/s)^2? - units not clear from AHF documentation

        Returns
        -------
        IE : array of dtype float
        """
        return self.profiles[:,26]

    def temps(self):
        """
        Returns an array containing the temperatures of gas for each profile.
        Calculated from the gas internal energy.
        Units of K

        WARNING:
        It is currently unclear what the gas thermal energy column represents,
        and so results using this function may be highly unreliable.

        Returns
        -------
        temps : array of dtype float
        """
        print("WARNING: Temperature function may be using incorrect units for calculation")
        return analysis.UtoT(self.intEnergies())
        #IE = self.intEnergies()
        #gasM_in_r = self.gasM_in_r() * apconst.M_sun.value
        #out = np.array([0.0]*len(gasM_in_r))
        #U = np.divide(IE,gasM_in_r,out=out,where=gasM_in_r!=0)
        #return analysis.UtoT(U)

'''========================================================='''
'''                     SNAPSHOT CLASS                      '''
'''========================================================='''

class Snapshot:
    """
    Contains information about a single snapshot of a simulation.

    WARNING: Snapshot objects have limited functionality, and are essentially
    just used to write truncated .AHF_halos and .AHF_profiles files after
    specifying a haloLimit during initialisation. File truncation should however
    be done by the ahfhalotools.filetools.truncateFiles() method, as that allows
    for truncation of mtree and mtree_idx files.
    """
    def __init__(self, snapNo, z, haloFile, profileFile, haloLimit = -1):
        """
        Initialises the Snapshot instance.

        Parameters
        ----------
        snapNo : int
            The snapshot number
        z : float
            The redshift of the snapshot
        haloFile : str
            The path to the .AHF_halos file to read from
        profileFile : str
            The path to the .AHF_profiles file to read from
        haloLimit : int, optional
            Sets number of halos at which file reading terminates
        """
        self.snapNo = snapNo
        self.z = z
        #age in Gyr
        self.age = analysis.tfromz(z)
        self.__haloFile__ = haloFile
        self.__profileFile__ = profileFile

        '''--- read in from files specified ---'''
        if haloLimit == -1:
            halorows = np.genfromtxt(haloFile)
            profilerows = np.genfromtxt(profileFile)
        else:
            halorows = np.genfromtxt(haloFile, max_rows = haloLimit)
            if halorows.ndim == 1:
                #make sure halorows is a 2D array (happens if numHalos = 1)
                halorows = np.array([halorows])

            numProfileRows = int(halorows[:,36].sum())
            profilerows = np.genfromtxt(profileFile, max_rows = numProfileRows)

        halos = []

        profileRowCounter = 0

        for rowIndex, halorow in enumerate(halorows):

            hostHalo = int(halorow[1])
            nbins = int(halorow[36])
            #if ((hostHalo % 1000000000000) > haloLimit or hostHalo == 0) and\
            #    rowIndex >= haloLimit:
            #    #halo is not a subhalo of a halo of interest, nor is a halo
            #    # of interest itself, so skip
            #    continue

            halo = Halo(halorow,z)
            prevRadius = -float('inf')

            #read in nbins lines from profiles file if halo is one of first
            # haloLimit halos
            if rowIndex < haloLimit:
                for i in range(nbins):
                    assert(not (profilerows[profileRowCounter,0] <= 0 and prevRadius > 0))
                    #profile row is still for same halo
                    halo.addProfile(profilerows[profileRowCounter])

                    prevRadius = profilerows[profileRowCounter,0]
                    profileRowCounter += 1

            #loop has broken because next row to consider is for next halo

            halos.append(halo)

            if len(halos) == haloLimit:
                break

        self.halos = halos

    def writeFiles(self, filenamebase):
        """
        Writes data saved in snapshot instance into new files with names
        filenamebase.AHF_halos and filenamebase.AHF_profiles

        Parameters
        ----------
        filenamebase : str
            The base name (including path) of the files to be outputted.
            Extensions will be added by the method.

        Notes
        -----
        This method can be used to truncate AHF_halos and AHF_profiles file in
        conjunction with setting a haloLimit during initialisation of the
        instance. File truncation should however be done using the
        ahfhalotools.filetools.truncateFiles() method.
        """
        profileRows = []
        haloRows = []
        for halo in self.halos:
            haloRows.append(halo.rawhalodata)
            for profile in halo.profiles:
                profileRows.append(profile)
        np.savetxt(filenamebase+".AHF_profiles",profileRows)
        np.savetxt(filenamebase+".AHF_halos",haloRows)

    @staticmethod
    def loadFiles(snapNos,zs,filebase,profileExt=".AHF_profiles",haloExt=".AHF_halos",haloLimit=-1):
        """
        Loads an array of snapshot instances from a file name format and a list
        of snap numbers and redshifts.

        Parameters
        ----------
        snapNos : list of int
            The snapshot numbers to load
        zs : list of float
            The redshifts corresponding to the snapshot numbers in snapNos
        filebase : str
            The format of the base file name. Should be formatted as such:
            "GIZMO-NewMDCLUSTER_0001.snap_{snap:0=3d}.z{z:.3f}"
            {snap:0=3d} is replaced with the snapshot number, including leading
            zeroes if number is less than 3 digits.
            z:.3f is replaces with the redshift, to three decimal places.
            Extensions are added on by the method
        profileExt : str, default = ".AHF_profiles"
            The file extension for the AHF_profiles files
        haloExt : str, default = ".AHF_halos"
            The file extension for the AHF_halos files
        haloLimit : int , optional
            The maximum number of halos to read before file reading stops.

        Returns
        -------
        snaps : array of dtype Snapshot
            An array of Snapshot objects corresponding to each of the snapshots
            specified in snapNos and zs parameters.
        """
        snaps = []
        noFiles = len(zs)
        for i in range(noFiles):
            snaps.append(Snapshot(snapNos[i],zs[i], filebase.format(snap=snapNos[i],z=zs[i])+haloExt,
                                  filebase.format(snap=snapNos[i],z=zs[i])+profileExt,
                                  haloLimit = haloLimit))
            #print("Loading: {0}%".format(i/noFiles * 100))
        return np.array(snaps)

'''========================================================='''
'''                     CLUSTER CLASS                       '''
'''========================================================='''

class Cluster:
    """
    Stores information about every halo in every snapshot in the cluster over time
    Designed for data produced by AMIGA's Halo Finder (AHF)

    Attributes
    ----------
    simName : str
        Specifies the name of the simulation used to create the cluster, if
        specified during initialisation
    clusterNum : int
        Specifies the number of the cluster that the object belongs to
        Defaults to 0 if not specified during initialisation
    """
    def __init__(self, fileBaseName, snapNums, zs, profileExt=".AHF_profiles",
            haloExt=".AHF_halos",mtreeidxExt=".AHF_mtree_idx",
            mtreeExt=".AHF_mtree", haloLimit=np.inf, clusterNum = 0,
            simName = "", skipmtree = False):
        """
        Initialises the cluster object from .AHF_profiles and .AHF_halos files

        Parameters
        ----------
        fileBaseName : str
            Base name of files without extensions. Should be formatted like:
            GIZMO-NewMDCLUSTER_0001.snap_{snap:0=3d}.z{z:.3f}
        snapNums : list of ints
            List of snap numbers to load files from
        zs : list of floats
            List of redshifts corresponding to snap numbers in snapNums
        profileExt : str, default = ".AHF_profiles"
        haloExt : str, default = ".AHF_halos"
        mtreeidxExt : str, default = ".AHF_mtree_idx"
        mtreeExt : str, default = ".AHF_mtree"
        haloLimit : int, optional
            Specifies maximum number of halos to load into memory per snapshot
            (see Notes)
        clusterNum : int, optional
            Specifies the number of the cluster that the object belongs to
            Defaults to 0
        simName : str, optional
            Specifies the name of the simulation used to create the cluster
            Defaults to ""
        skipmtree : bool, optional
            Whether or not to skip trying to read the merger tree. If False and
            no mtree file is present, funcation will just print out a warning.
            Defaults to false

        Notes
        -----
        haloLimit specifies point at which file reading stops - as such since
        halos can change in order in files over time, a halo initially at the
        end of the file can be truncated off in earler snapshots, resulting in
        an inability to track it through time. It is therefore ideal to specify
        haloLimit as a number larger than the number of halos you intend to
        investigate.

        Note that haloLimit will not limit the number of halos read in from
        .AHF_halos files, as extra halos are included during file truncation
        if they are subhalos of the first n halos (where n is the size of
        truncation)

        Also note that the file format for mtree is NOT SUSSING-2013, this
        class reads files of format specified in the AHF documentation:
            http://popia.ft.uam.es/AHF/files/AHF.pdf (P190)
        """
        #load in data
        assert(len(snapNums)==len(zs))

        self.simName = simName
        self.clusterNum = clusterNum

        #haloDict maps haloID -> halo object
        self._haloDict = dict()
        self._fatherDict = dict()

        #self._progDict maps int haloID ->  [[sharedPart, progenitorID, progenitorPart]]
        self._progDict = dict()

        #self._encHaloDict maps int hostHaloID - > list (HaloLite)
        self._encHaloDict = dict()

        for i in range(len(snapNums)):
            snapNum = snapNums[i]
            z = zs[i]
            #get file names from base provided
            fileNameFmt = fileBaseName.format(snap=snapNum,z=z)
            haloFile = fileNameFmt + haloExt
            profileFile = fileNameFmt + profileExt
            mtreeidxFile = fileNameFmt + mtreeidxExt
            mtreeFile = fileNameFmt + mtreeExt

            #loading halo and profile data
            halorows = np.genfromtxt(haloFile)
            profilerows = np.genfromtxt(profileFile)

            #check if files are just 1 row and if so convert them to correct
            # data type (2D array)
            if halorows.ndim == 1:
                halorows = np.array([halorows])
            if profilerows.ndim == 1:
                profilerows = np.array([profilerows])

            profileRowCounter = 0
            haloCounter = 0

            for halorow in halorows:

                halo = Halo(halorow,z)
                prevRadius = -float('inf')

                #check that we have not already reached end of profiles file,
                # nor have we read the profiles of haloLimit halos
                if profileRowCounter < len(profilerows) and haloCounter < haloLimit:
                    #start reading profile rows
                    while not (profilerows[profileRowCounter,0] <= 0 and prevRadius > 0):
                        #profile row is still for same halo
                        halo.addProfile(profilerows[profileRowCounter])

                        prevRadius = profilerows[profileRowCounter,0]
                        profileRowCounter += 1

                        #break at end of file
                        if profileRowCounter == len(profilerows):
                            break
                    #while loop has broken because next row to consider is for next halo

                assert(halo.ID not in self._haloDict)
                self._haloDict[halo.ID] = halo
                haloCounter += 1

            #there is no mtree or mtree_idx for the lowest snapshot in a simulation
            # so we check if we are dealing with the lowest snapshot and if so
            # we don't try and load mtree or mtree_idx
            if snapNum == min(snapNums):
                continue
            #loading merge tree idx data
            #mtree_idx file is formatted as childID fatherID
            try:
                idxrows = np.genfromtxt(mtreeidxFile)
                #check if idxrows is 1 row and if so convert it to correct
                # data type (2D array)
                if idxrows.ndim == 1:
                    idxrows = np.array([idxrows])
            except IOError:
                print("WARNING: File {0} not found, relevant data cannot be loaded".format(mtreeidxFile))
                idxrows = []
            #fatherDict maps childID -> fatherID (key is childID)
            rowCounter = 0
            for idxrow in idxrows:
                self._fatherDict[int(idxrow[0])] = int(idxrow[1])
                rowCounter += 1
                if rowCounter == haloLimit:
                    break

            #mtree format:
            #   HaloID(1)   HaloPart(2)  NumProgenitors(3)
            #      SharedPart(1)    HaloID(2)   HaloPart(3)
            #note: mtree includes the father
            if not skipmtree:
                try:
                    mtreerows = np.genfromtxt(mtreeFile)
                except IOError:
                    print("WARNING: File {0} not found, relevant data cannot be loaded".format(mtreeFile))
                    mtreerows = []
                lineIndex = 0
                fileLength = len(mtreerows)

                while lineIndex < fileLength-1:
                    haloID, nPart, numProg = mtreerows[lineIndex,:]

                    progenitorList = []
                    #loop through progenitor lines
                    for i in range(lineIndex+1,lineIndex+int(numProg)+1):
                        sharedPart, progenitorID, prognPart = mtreerows[i,:]
                        progenitorList.append([int(sharedPart),int(progenitorID),int(prognPart)])

                    #add list of projenitor dictionaries to _progDict
                    self._progDict[haloID] = np.array(progenitorList)
                    #move line index to next child line
                    lineIndex += int(numProg) + 1


    def getHalo(self, haloID):
        """
        Gets the Halo instance corresponding to the halo ID. If the ID is
        invalid, will return -1

        Parameters
        ----------
        haloID : int
            The full haloID of the halo
            e.g. 128000000000002, 102000000000023

        Returns
        -------
        halo : Halo object
        """
        try:
            halo = self._haloDict[haloID]
        except KeyError:
            halo = -1
        return halo

    def getFatherOf(self, haloID):
        """
        Returns the haloID of the father halo

        Parameters
        ----------
        haloID : int
            The full haloID of the halo
            e.g. 128000000000002, 102000000000023

        Returns
        -------
        fatherID : int
            The haloID of the father halo
            If no such father exists (either due to truncation or end of file),
            returns -1
        """
        try:
            return self._fatherDict[haloID]
        except KeyError:
            #no father stored for halo, so will return -1
            return -1

    def getMergeTreeEntry(self, haloID):
        """
        Returns merger tree entry for halo.

        Parameters
        ----------
        haloID : int
            The full haloID of the halo

        Returns
        -------
        mtree_entry : array of 3-tuples
            Formatted as [[sharedPart, progenitorID, progenitorPart]]
            (Same format as .AHF_mtree entry)
        """
        try:
            return self._progDict[haloID]
        except KeyError:
            #haloID is probably in earliest snapshot loaded
            return None

    def trackID(self, haloID):
        """
        Takes a haloID and returns an array of haloIDs corresponding to
        the chain of father halos backwards in time
        """
        fatherID = self.getFatherOf(haloID)
        if fatherID == -1:
            #end of chain, return just current haloID
            return [haloID]
        else:
            return self.trackID(fatherID) + [haloID]

    def trackHalo(self,haloID):
        """
        Takes a haloID and returns an array of halos corresponding to
        the chain of father halos backwards in time. If there are gaps in
        stored halos, those halos will be skipped

        See Also
        --------
        Cluster.trackID()
        """
        chainIDs = self.trackID(haloID)
        halos = []
        for id in chainIDs:
            halo = self.getHalo(id)
            if halo != -1: halos.append(halo)
        return halos

    def getHaloData(self, haloID, quantity):
        """
        Used to retrieve values of halo data from a specified halo

        Parameters
        ----------
        haloID : int
            The full haloID of the halo to investigate at z=0.
            e.g. 128000000000001
        quantity : str
            The name of the halo data piece to return
            e.g. "Rvir", "Mvir", "Ekin", "Qvir", etc
            Does not accept delta quantities

        Returns
        -------
        value
            The halo data value corresponding to quantity specified
        """
        halo = self.getHalo(haloID)
        try:
            return halo.halodata[quantity]
        except KeyError as err:
            msg = "Invalid Halo Data quantity key {0}. To retreive valid "\
                  "halodata keys, one can call Halo.halodata.keys() on a Halo "\
                  "instance.".format(quantity)
            raise KeyError(msg) from err
        except AttributeError:
            #self.getHalo has returned -1, ie we do not have the halo specified
            #loaded into memory
            print("WARNING: Halo {0} not loaded into memory for {1} cluster {2}".format(haloID,self.simName,self.clusterNum))
            return -1

    def funcOfZHaloData(self, haloID, quantity):
        """
        Returns a piece of halo data as a function of redshift. To be used for
        quantities stored in the .AHF_halos file.

        Parameters
        ----------
        haloID : int
            The full haloID of the halo to investigate at z=0.
            e.g. 128000000000001
        quantity : str
            The name of the halo data piece to return
            e.g. "Rvir", "Mvir", "Ekin", "Qvir", etc
            Returns time delta of quantities when key is preceded by "delta",
            e.g. "deltaRvir", "deltaMvir", "deltaEkin", etc

        Returns
        -------
        zs : list of float
            Array of redshifts that correspond to each of the values in values
        values : list
            Array of values specified by quantity parameter as a function
            of time

        See Also
        --------
        Cluster.funcOfAgeHaloData : get halo data as a function of age
        Cluster.funcOfZDeltaHaloData :
            get the time delta of a halo data quantity as a function of redshift
        Cluster.funcOfAgeDeltaHaloData :
            get the time delta of a halo data quantity as a function of age
        """
        if quantity.startswith("delta"):
            #return delta value instead of normal value
            return self.funcOfZDeltaHaloData(haloID,quantity[5:])

        #for normal values:
        halos = self.trackHalo(haloID)
        zs = np.array([halo.z for halo in halos])
        values = None
        try:
            values = np.array([halo.halodata[quantity] for halo in halos])
        except KeyError as err:
            msg = "Invalid Halo Data quantity key {0}. To retreive valid "\
                  "halodata keys, one can call Halo.halodata.keys() on a Halo "\
                  "instance.".format(quantity)
            raise KeyError(msg) from err
        return zs, values

    def funcOfAgeHaloData(self, haloID, quantity):
        """
        Returns a piece of halo data as a function of age. To be used for
        quantities stored in the .AHF_halos file.

        Parameters
        ----------
        haloID : int
            The full haloID of the halo to investigate at z=0.
            e.g. 128000000000001
        quantity : str
            The name of the halo data piece to return
            e.g. "Rvir", "Mvir", "Ekin", "Qvir", etc
            Returns time delta of quantities when key is preceded by "delta",
            e.g. "deltaRvir", "deltaMvir", "deltaEkin", etc

        Returns
        -------
        ages : list of float
            Array of ages in Gyr that correspond to each of the values in values
        values : list
            Array of values specified by quantity parameter as a function
            of time

        See Also
        --------
        Cluster.funcOfZHaloData : get halo data as a function of redshift
        Cluster.funcOfZDeltaHaloData :
            get the time delta of a halo data quantity as a function of redshift
        Cluster.funcOfAgeDeltaHaloData :
            get the time delta of a halo data quantity as a function of age
        """
        zs, values = self.funcOfZHaloData(haloID,quantity)
        return analysis.tfromz(zs), values

    def funcOfZDeltaHaloData(self, haloID, quantity):
        """
        Returns the change (delta) in a piece of halo data as a function of
        redshift. To be used for quantities stored in the .AHF_halos file.

        Delta q at redshift z2 is calculated as q(z2)-q(z1) where z1 is the
        redshift of the previous snapshot. Note that this requires snapshots
        to have relatively even redshift spacing for delta quantities to be
        comparible over time.

        Parameters
        ----------
        haloID : int
            The full haloID of the halo to investigate at z=0
        quantity : str
            The name of the halo data piece to return
            e.g. "Rvir", "Mvir", "Ekin", "Qvir", etc
            Note that the quantity should not include a "delta"

        Returns
        -------
        zs : list of float
            The redshifts correspond to each of the values in deltaValues
        deltaValues : list of float
            The delta values of the quantity specified

        See Also
        --------
        Cluster.funcOfZHaloData : get halo data as a function of redshift
        Cluster.funcOfAgeHaloData : get halo data as a function of age
        Cluster.funcOfAgeDeltaHaloData :
            get the time delta of a halo data quantity as a function of age
        """
        if quantity.startswith("delta"):
            quantity = quantity[5:]

        #get redshifts and values for non-delta quantity
        zs, values = self.funcOfZHaloData(haloID, quantity)

        #calculate delta values
        deltaValues = values[1:]-values[:-1]
        #take off highest redshift from zs (can't calculate delta quantity
        # without an older snapshot)
        zs = zs[1:]

        return zs, deltaValues

    def funcOfAgeDeltaHaloData(self, haloID, quantity):
        """
        Returns the change (delta) in a piece of halo data as a function of
        age in Gyr. To be used for quantities stored in the .AHF_halos file.

        Delta q at age t2 is calculated as q(t2)-q(t1) where t1 is the
        age of the previous snapshot. Note that this requires snapshots
        to have relatively even redshift spacing for delta quantities to be
        comparible over time.

        Parameters
        ----------
        haloID : int
            The full haloID of the halo to investigate at z=0
        quantity : str
            The name of the halo data piece to return
            e.g. "Rvir", "Mvir", "Ekin", "Qvir", etc
            Note that the quantity should not include a "delta"

        Returns
        -------
        ts : list of float
            The ages in Gyr correspond to each of the values in deltaValues
        deltaValues : list of float
            The delta values of the quantity specified

        See Also
        --------
        Cluster.funcOfZHaloData : get halo data as a function of redshift
        Cluster.funcOfAgeHaloData : get halo data as a function of age
        Cluster.funcOfZDeltaHaloData :
            get the time delta of a halo data quantity as a function of redshift
        """
        zs, deltaValues = self.funcOfZDeltaHaloData(haloID,quantity)
        return analysis.tfromz(zs), deltaValues

    def funcOfZProfileData(self, haloID, quantity, radius):
        """
        Returns a radius-dependent quantity as a function of redshift,
        interpolated at the given radius. To be used for quantities stored
        in the .AHF_profiles files.

        Parameters
        ----------
        haloID : int
            The full haloID of the halo to investigate at z=0.
            e.g. 128000000000001
        quantity : str
            The name of the profile data piece to return
            e.g. "M_in_r", "encDens", "gasM_in_r", etc
            Can alternatively be an integer that represents the index of a
            profile data piece (eg for Ly in column (10), can give quantity = 9)
            Note this is zero indexed.
        radius : float
            The radius at which to interpolate the quantity specified

        Returns
        -------
        zs : list of float
            Array of redshifts that correspond to each of the values in values
        values : list
            Array of values specified by quantity parameter as a function
            of time at specified radius

        See Also
        --------
        Cluster.funcOfAgeProfileData
        Cluster.funcOfZDeltaProfileData
        Cluster.funcOfAgeDeltaProfileData
        """
        if type(quantity) == type(str()) and quantity.startswith("delta"):
            return self.funcOfZDeltaProfileData(haloID,quantity[5:],radius)
        halos = self.trackHalo(haloID)
        zs = np.array([halo.z for halo in halos])
        values = []
        for halo in halos:
            rs, rvals = self.funcOfRadiusProfileData(halo.ID,quantity)
            interp = np.interp(radius, rs, rvals)
            values.append(interp)
        values = np.array(values)
        return zs, values

    def funcOfAgeProfileData(self, haloID, quantity, radius):
        """
        Returns a radius-dependent quantity as a function of age,
        interpolated at the given radius. To be used for quantities stored
        in the .AHF_profiles files.

        Parameters
        ----------
        haloID : int
            The full haloID of the halo to investigate at z=0.
            e.g. 128000000000001
        quantity : str
            The name of the profile data piece to return
            e.g. "M_in_r", "encDens", "gasM_in_r", etc
            Can alternatively be an integer that represents the index of a
            profile data piece (eg for Ly in column (10), can give quantity = 9)
            Note this is zero indexed.
        radius : float
            The radius at which to interpolate the quantity specified

        Returns
        -------
        ages : list of float
            Array of ages in Gyr that correspond to each of the values in values
        values : list
            Array of values specified by quantity parameter as a function
            of time at specified radius

        See Also
        --------
        Cluster.funcOfZProfileData
        Cluster.funcOfZDeltaProfileData
        Cluster.funcOfAgeDeltaProfileData
        """
        zs, values = self.funcOfZProfileData(haloID,quantity,radius)
        return analysis.tfromz(zs),values

    def funcOfZDeltaProfileData(self, haloID, quantity, radius):
        """
        Returns the time delta of a radius-dependent quantity as a function of
        redshift, interpolated at the given radius.

        To be used for quantities stored in the .AHF_profiles files.

        Parameters
        ----------
        haloID : int
            The full haloID of the halo to investigate at z=0.
            e.g. 128000000000001
        quantity : str
            The name of the halo data piece to return
            e.g. "M_in_r", "encDens", "gasM_in_r", etc
            Note that the quantity should not include a "delta"
            Can alternatively be an integer that represents the index of a
            profile data piece (eg for Ly in column (10), can give quantity = 9)
            Note this is zero indexed.
        radius : float
            The radius at which to interpolate the quantity specified

        Returns
        -------
        zs : list of float
            Array of redshifts that correspond to each of the values in
            deltaValues
        deltaValues : list
            Array of values specified by quantity parameter as a function
            of time at specified radius

        See Also
        --------
        Cluster.funcOfZProfileData
        Cluster.funcOfAgeProfileData
        Cluster.funcOfAgeDeltaProfileData
        """
        if type(quantity) == type(str()) and quantity.startswith("delta"):
            quantity = quantity[5:]

        #get redshifts and values for non-delta quantity
        zs, values = self.funcOfZProfileData(haloID, quantity, radius)

        #calculate delta values
        deltaValues = values[1:]-values[:-1]
        #take off highest redshift from zs (can't calculate delta quantity
        # without an older snapshot)
        zs = zs[1:]
        return zs, deltaValues

    def funcOfAgeDeltaProfileData(self, haloID, quantity, radius):
        """
        Returns the time delta of a radius-dependent quantity as a function of
        age, interpolated at the given radius.

        To be used for quantities stored in the .AHF_profiles files.

        Parameters
        ----------
        haloID : int
            The full haloID of the halo to investigate at z=0.
            e.g. 128000000000001
        quantity : str
            The name of the halo data piece to return
            e.g. "M_in_r", "encDens", "gasM_in_r", etc
            Note that the quantity should not include a "delta"
            Can alternatively be an integer that represents the index of a
            profile data piece (eg for Ly in column (10), can give quantity = 9)
            Note this is zero indexed.
        radius : float
            The radius at which to interpolate the quantity specified

        Returns
        -------
        ages : list of float
            Array of ages in Gyr that correspond to each of the values in
            deltaValues
        deltaValues : list
            Array of values specified by quantity parameter as a function
            of time at specified radius

        See Also
        --------
        Cluster.funcOfZProfileData
        Cluster.funcOfAgeProfileData
        Cluster.funcOfZDeltaProfileData
        """
        zs, deltaValues = self.funcOfZDeltaProfileData(haloID, quantity, radius)
        ages = analysis.tfromz(zs)
        return ages, deltaValues

    def funcOfRandZProfileDataPoints(self, masterHaloID, quantity, removeZeroes = False):
        """
        Returns three 1d arrays that correspond to the coordinates of points
        with x = radius, y = redshift, z = quantity

        Uses radii that appear in AHF_profiles for each snapshot. Output comes
        as triangular grid (see notes).

        Parameters
        ----------
        masterHaloID : int
            The full haloID of the halo being investigated in the *oldest*
            snapshot of interest
        quantity : str
            The name or column index of the profile data piece to use
        removeZeroes : bool
            Whether or not to remove points where the z value is 0. Default = False

        Returns
        -------
        r : array of float
            Array of radii in kpc/h
        z : array of float
            Array of redshifts
        val : array of float
            Array of values

        Notes
        -----
        The values returned from this function form a triangular grid, ideal for
        plotting with using matplotlib.pyplot.triplot or
        matplotlib.pyplot.plot_trisurf
        The values in z and r will not form a 2D grid, as the radii for each
        bin changes from snapshot to snapshot.

        See Also
        --------
        matplotlib.pyplot.triplot :
            Draw a unstructured triangular grid as lines and/or markers.
        mpl_toolkits.mplot3d.axes3d.Axes3D.plot_trisurf :
            Plot a triangulated surface
        """
        haloIDs = self.trackID(masterHaloID)
        rs = np.empty(0)
        zs = np.empty(0)
        vals = np.empty(0)
        for haloID in haloIDs:

            #add rows to grid
            rrow, valrow = self.funcOfRadiusProfileData(haloID, quantity)
            zrow = [self.getHalo(haloID).z]*len(rrow)
            rs = np.concatenate((rs,rrow))
            zs = np.concatenate((zs,zrow))
            vals = np.concatenate((vals,valrow))

        if removeZeroes:
            #search for 0s in vals
            indexes = np.argwhere(vals==0)
            vals = np.delete(vals, indexes)
            rs = np.delete(rs, indexes)
            zs = np.delete(zs, indexes)
        return rs,zs,vals


    def funcOfRadiusProfileData(self, haloID, quantity, removeZeroes = False):
        """
        Returns the profile quantity as a function of radius at the redshift
        corresponding to the snapshot number specified in haloID

        Parameters
        ----------
        haloID : int
            The full haloID of the halo to investigate
            e.g. 128000000000001
        quantity : str
            The name of the halo data piece to return
            e.g. "M_in_r", "encDens", "gasM_in_r", etc
            Can alternatively be an integer that represents the index of a
            profile data piece (eg for Ly in column (10), can give quantity = 9)
            Note this is zero indexed.
        removeZeroes : bool, default = False
            If True will remove any entries where the value of {quantity} is 0

        Returns
        -------
        rs : array of float
            Array of radii that correspond to values in values
        values : array
            Array of values specified by quantity parameter as a function
            of radius
        """
        halo = self.getHalo(haloID)
        rs = halo.radii()

        if type(quantity) == type(int()):
            #have passed an index rather than a string key
            values = halo.profiles[:,quantity]
        else:
            quantitydict = halo.getProfileQuantityDict()
            values = None
            try:
                values = quantitydict[quantity]()
            except KeyError as err:
                message = "Invalid profile data quantity key {0}. To retrieve valid "\
                          "profile data keys, one can call "\
                          "Halo.getProfileQuantityDict() on a Halo instance, or use"\
                          " the quantity's column index (eg for Ly in column (10), "\
                          "can give quantity = 9)".format(quantity)
                raise KeyError(message) from err
        if removeZeroes:
            #search for 0s in vals
            indexes = np.argwhere(values==0)
            values = np.delete(values, indexes)
            rs = np.delete(rs, indexes)
        rs = np.array(rs)
        values = np.array(values)
        return rs, values

    def getMergeSize(self, haloID, scheme = "halodata", fractional = True):
        """
        Calculates the size of the merger that the halo experienced between
        snapshots.

        Parameters
        ----------
        haloID : int
            The halo ID of the halo being investigated, in the snapshot after
            the merge.
            (i.e. to get the size of the merger between snapshot 102 to 103, one
            would use the ID 10300000000000x)
        scheme : str, default = "halodata"
            Defines which scheme to use for calculating size of mergers.
            Valid values are "mtree-largest", "mtree-sum", "mtree-first", or
            "halodata"
            To use "mtree-largest" or "mtree-sum", one must have already loaded
            a .BDP_enchalos file for the master halo.
            For an explanation of the schemes please see:
            https://github.com/BenDavisonPetch/ahfhalotools/blob/main/merger_detection_scheme_guide.md
        fractional : bool, default = True
            Specifies whether the returned size should be the fractional size
            (returned as a decimal), or the absolute size (returned as the
            number of particles).

        Returns
        -------
        size : number
            The size of the merger. If fractional is true, size will be a
            decimal float representing the fractional merger size, and if
            false, will be an integer number of particles.
            If the scheme is an mtree scheme, and the halo has no merger tree
            entry, will return -1.

        See Also
        --------
        TODO
        """
        validSchemes = ["mtree-largest","mtree-sum","mtree-first","halodata"]
        if scheme not in validSchemes:
            msg = "Scheme '{0}' is not a valid scheme! Accepted values are: \n{1}".format(scheme,", ".join(validSchemes))
            raise ValueError(msg)

        if scheme == "mtree-largest" or scheme == "mtree-sum" or scheme == "mtree-first":
            #get merge tree entry
            mtreeEntry = self.getMergeTreeEntry(haloID)

            #if there is no merge tree entry for halo, return -1
            if type(mtreeEntry) == type(None):
                msg = "Halo {0} has no merge tree entry! Returning -1".format(haloID)
                warnings.warn(msg, RuntimeWarning)
                return -1

            #get list of IDs of sub halos of father halo
            fatherID = self.getFatherOf(haloID)
            #assert(mtreeEntry[0,1]==fatherID)
            subHaloIDs = [halo.ID for halo in self.getEnclosedHalos(fatherID) if halo.hostHalo == fatherID]
            #go through merge tree entries and ignore any that are from
            # halos that were subhalos of the father halo
            #print("looking at halo {0}:".format(haloID))
            mergeSize = 0
            for mtreeRow in mtreeEntry[1:]:
                #(we exclude first row as that will be the father halo)
                progID = mtreeRow[1]
                if progID not in subHaloIDs:
                    #print("   non sub halo progenitor found: {0}".format(progID))
                    #print("   shared size is {0}".format(mtreeRow[0]))
                    #we have a progenitor that was not previously a subhalo
                    # therefore the particles introduced are new

                    #if using mtree-first we can stop here, if mtree-sum
                    # we keep going and add all of the non subhalo
                    # contributions, if mtree-largest we set mergeSize to
                    # the largest sharedPart we find
                    if scheme == "mtree-first":
                        mergeSize = mtreeRow[0]
                        break
                    elif scheme == "mtree-largest":
                        mergeSize = max(mergeSize,mtreeRow[0])
                    elif scheme == "mtree-sum":
                        mergeSize += mtreeRow[0]

            #print("mergeSize is now {0}\n".format(mergeSize))

            if fractional:
                fatherSize = mtreeEntry[0,2]
                mergeSize /= fatherSize

            return mergeSize
        elif scheme == "halodata":
            #scheme only relies on delta M and M from .AHF_halos file
            fatherID = self.getFatherOf(haloID)
            if fatherID == -1:
                return -1

            #calculate delta M
            fatherM = self.getHaloData(fatherID,"Mvir")
            childM = self.getHaloData(haloID,"Mvir")
            mergeSize = (childM-fatherM)

            #if fractional divide by father mass
            if fractional:
                mergeSize /= fatherM
            return mergeSize


    def getMergeZs(self,haloIDtoTrack,threshold=0.2,scheme="halodata",fractional=True):
        """
        Returns the red shifts at which halo experiences merger above specified
        threshold

        Parameters
        ----------
        haloIDtoTrack : int
            The full haloID of haloID to track - will track backwards through
            time
        threshold : float, default = 0.2
            The threshold for merge size
        scheme : str, default = "halodata"
            Defines which scheme to use for calculating size of mergers.
            Valid values are "mtree-largest", "mtree-sum", "mtree-first", or
            "halodata"
            To use "mtree-largest" or "mtree-sum", one must have already loaded
            a .BDP_enchalos file for the master halo.
            For an explanation of the schemes please see:
            https://github.com/BenDavisonPetch/ahfhalotools/blob/main/merger_detection_scheme_guide.md
        fractional : bool, default = True
            Specifies whether the merger size should be the fractional size
            as a decimal, or the absolute size as the number of particles
            (mtree) / mass increment (halodata).

        Returns
        -------
        zs : array of floats
            The redshifts that correspond to the snapshots just after the merges
            occur
        sizes : array of floats
            The sizes of the mergers

        See Also
        --------
        TODO
        """
        zs = []
        sizes = []
        for haloID in self.trackID(haloIDtoTrack)[1:]:
            #get merge size
            mergeSize = self.getMergeSize(haloID,scheme=scheme,fractional=fractional)
            #check if merger is larger than threshold, if it is add redshift and
            # size to list to return
            if mergeSize > threshold:
                zs.append(self.getHalo(haloID).z)
                sizes.append(mergeSize)
        return np.array(zs), np.array(sizes)

    def getMergeTimes(self,haloIDtoTrack,threshold=0.2,scheme="halodata",fractional=True):
        """
        Returns the ages at which halo experiences merger above specified
        threshold

        Parameters
        ----------
        haloIDtoTrack : int
            The full haloID of haloID to track - will track backwards through
            time
        threshold : float, default = 0.2
            The threshold for merge size
        scheme : str, default = "halodata"
            Defines which scheme to use for calculating size of mergers.
            Valid values are "mtree-largest", "mtree-sum", "mtree-first", or
            "halodata"
            To use "mtree-largest" or "mtree-sum", one must have already loaded
            a .BDP_enchalos file for the master halo.
            For an explanation of the schemes please see:
            https://github.com/BenDavisonPetch/ahfhalotools/blob/main/merger_detection_scheme_guide.md
        fractional : bool, default = True
            Specifies whether the merger size should be the fractional size
            as a decimal, or the absolute size as the number of particles
            (mtree) / mass increment (halodata).

        Returns
        -------
        ages : array of floats
            The ages that correspond to the snapshots just after the merges
            occur
        sizes : array of floats
            The sizes of the mergers

        See Also
        --------
        TODO
        """
        zs, sizes = self.getMergeZs(haloIDtoTrack,threshold=threshold,scheme=scheme,fractional=fractional)
        return analysis.tfromz(zs), sizes

    def getLargestMergeZInRange(self, masterHaloID, minZ, maxZ, scheme="halodata", fractional = True):
        """
        Returns the redshift of the largest merger in the range of redshifts specified

        Parameters
        ----------
        masterHaloID : int
            Full haloID of halo to track in time in most recent snapshot
        minZ : float
            Defines minimum redshift for search range (inc.)
        maxZ : float
            Defines maximum redshift for search range (inc.)
        scheme : str, default = "halodata"
            Defines which scheme to use for calculating size of mergers.
            Valid values are "mtree-largest", "mtree-sum", "mtree-first", or
            "halodata"
            To use "mtree-largest" or "mtree-sum", one must have already loaded
            a .BDP_enchalos file for the master halo.
            For an explanation of the schemes please see:
            https://github.com/BenDavisonPetch/ahfhalotools/blob/main/merger_detection_scheme_guide.md
        fractional : bool, default = True
            Specifies whether the merger size should be the fractional size
            as a decimal, or the absolute size as the number of particles
            (mtree) / mass increment (halodata).

        Returns
        -------
        z : float
            The redshift at which the largest merger occurred in the range
            If the two largest mergers are of the same size will return redshift
            of earliest one.
            If there are no snapshots found with redshift inside range specified
            will be None
        largestMergeSize : int
            The size of the largest merge event found in range.
            If there are no snapshots found with redshift inside range specified
            will be 0
        """
        if maxZ < minZ:
            msg = "Expected maxZ >= minZ but got maxZ < minZ!"
            raise ValueError(msg)

        z = None
        largestMerge = 0
        #loop over chain of father halos, not including oldest snapshot
        # as we will not have data for its father's enclosed halos
        for haloID in self.trackID(masterHaloID)[1:]:
            #check if halo redshift is inside range specified, if not skip
            haloZ = self.getHalo(haloID).z
            if haloZ < minZ or haloZ > maxZ: continue

            #get merger size
            mergeSize = self.getMergeSize(haloID,scheme=scheme,fractional=fractional)

            if mergeSize > largestMerge:
                z = haloZ
                largestMerge = mergeSize
        return z, largestMerge

    """============================="""
    """   ENCLOSED HALO FUNCTIONS   """
    """============================="""

    def generateEnclosedHaloFile(self,hostHaloID,haloFile,outputFile,useAHFSubHalos=False):
        """
        Performs a search across all halos in haloFile as
        to whether they are contained within a sphere 2 times the virial radius
        of the host halo.
        Generates a file containing only the relevant information about only
        the halos that are enclosed within the host halo.

        Parameters
        ----------
        hostHaloID : int
        haloFile : string
        outputFile : string
        useAHFSubHalos : bool, default = False
            Specifies whether to use AHF sub halo detection instead of detection
            by position

        Notes
        -----
        (Unless useAHFSubHalos is True,)
        Does not use AHF's sub halo detection, due to how AHF defines subhalos
        (it is possible for an object to lie within a host halo and not be
        considered a sub halo if it is not within a common isodensity contour)

        See also
        --------
        Cluster.generateEnclosedHaloFilesFromChain()
        Cluster.loadEnclosedHaloFile()
        """
        if useAHFSubHalos:
            raise NotImplementedError
        else:
            hostHalo = self.getHalo(hostHaloID)
            halorows = np.genfromtxt(haloFile)

            header = ["ID(1)","hostHalo(2)","Xc(3)","Yc(4)","Zc(5)","VXc(6)","VYc(7)","VZc(8)"]
            #the relevant indexes of the above quantities in the AHF_halos file
            quantityIndexes = [0,1,5,6,7,8,9,10]
            outputRows = np.empty((0,len(header)))
            for halorow in halorows:
                haloID = halorow[0]
                #check that the halo being considered isn't the host halo
                if haloID == hostHaloID: continue

                #check if the position of the halo satisfied criteria for being
                # enclosed
                x = halorow[5]
                y = halorow[6]
                z = halorow[7]

                #calculate distance to centre of hostHalo
                dist = np.sqrt((hostHalo.Xc-x)**2+(hostHalo.Yc-y)**2+(hostHalo.Zc-z)**2)
                if dist < 2*hostHalo.Rvir:
                    #halo is enclosed by hostHalo
                    #we can add the relevant data to the output rows
                    outputRow = np.array([halorow[i] for i in quantityIndexes])
                    outputRows = np.vstack((outputRows,outputRow))

            #write output file
            np.savetxt(outputFile,outputRows,header=", ".join(header))

    def generateEnclosedHaloFilesFromChain(self,masterHaloID,haloFiles,outputFileFmt,useAHFSubHalos=False):
        """
        Generates an enclosed halo file for each halo in the chain of father halos
        of the master halo specifies.

        WARNING: May take a while to run

        Parameters
        ----------
        masterHaloID : int
            The full haloID of the halo to track in the most recent snapshot
        haloFiles : list ( str )
            List of halo files to use, ordered in *ascending* snapshot order
        outputFileFmt : str
            Format of output file. Should contain a {haloID}, which will be
            replaced by the halo ID of the host halo the file refers to.
            e.g. "GadgetX-NewMDCLUSTER_0001.halo{haloID}.BDP_enchalos"
        useAHFSubHalos : bool
            Specifies whether to use AHF sub halo detection instead of detection
            by position
        See also
        --------
        Cluster.generateEnclosedHaloFile()
        Cluster.loadEnclosedHaloFilesFromChain()
        """
        haloIDs = self.trackID(masterHaloID)
        assert(len(haloIDs)==len(haloFiles))

        for i in range(len(haloIDs)):
            haloID = haloIDs[i]
            haloFile = haloFiles[i]
            outputFile = outputFileFmt.format(haloID=haloID)
            self.generateEnclosedHaloFile(haloID,haloFile,outputFile,useAHFSubHalos=useAHFSubHalos)


    def loadEnclosedHaloFile(self,hostHaloID,fileName):
        """
        Loads data from an enclosed halo file (.BDP_enchalos) into the Cluster
        instance.

        Parameters
        ----------
        hostHaloID : int
            The full haloID of the host halo
        fileName : str

        Notes
        -----
        Will load subhalos in as HaloLite objects in an array accessable using
        Cluster.getEnclosedHalos(hostHaloID)

        See also
        --------
        Cluster.generateEnclosedHaloFile()
        Cluster.generateEnclosedHaloFilesFromChain()
        Cluster.getEnclosedHalos()
        """
        encHaloRows = np.genfromtxt(fileName)
        encHaloArray = []
        for encHaloRow in encHaloRows:
            halo = HaloLite(encHaloRow)
            encHaloArray.append(halo)

        self._encHaloDict[hostHaloID] = encHaloArray

    def loadEnclosedHaloFilesFromChain(self,masterHaloID,fileNameFmt):
        """
        Loads data from enclosed halo files (.BDP_enchalos) into the Cluster
        instance, for each halo in the chain of father halos of the master halo

        Parameters
        ----------
        masterHaloID : int
            The full haloID of the master halo in the most recent snapshot
        fileNameFmt : str
            The file name format of .BDP_enchalos file, with {haloID} as a
            placeholder for the haloID in the file name
            e.g. "GX//GadgetX-NewMDCLUSTER_0001.halo{haloID}.BDP_enchalos"

        Notes
        -----
        Will load subhalos in as HaloLite objects in an array accessable using
        Cluster.getEnclosedHalos(hostHaloID)

        See also
        --------
        Cluster.getEnclosedHalos(hostHaloID)
        Cluster.loadEnclosedHaloFile()
        Cluster.generateEnclosedHaloFilesFromChain()
        Cluster.generateEnclosedHaloFile()
        """
        haloIDs = self.trackID(masterHaloID)
        for haloID in haloIDs:
            self.loadEnclosedHaloFile(haloID,fileNameFmt.format(haloID=haloID))

    def getEnclosedHalos(self,hostHaloID):
        """
        Gives the halos enclosed by the host halo specified

        Parameters
        ----------
        hostHaloID : int

        Returns
        -------
        enclosedHalos : list ( HaloLite )

        Notes
        -----
        Enclosed halo data must be loaded in for the halo in advance using
        Cluster.loadEnclosedHaloFile

        See also
        --------
        Cluster.loadEnclosedHaloFile
        """
        try:
            return self._encHaloDict[hostHaloID]
        except KeyError as err:
            msg = "ERROR: Halo {0} does not have enclosed halos loaded into memory!\n"\
                  "Before using data on enclosed halos, data must be loaded using"\
                  "Cluster.loadEnclosedHaloFile() or Cluster.loadEnclosedHaloFilesFromChain()".format(hostHaloID)
            raise Exception(msg) from err

    def getRelVelocitiesOfEncHalos(self, hostHaloID):
        """
        Takes the halo ID of the host halo and returns an array of velocities of
        enclosed halos relative to host halo in km/s

        Parameters
        ----------
        hostHaloID : int

        Returns
        -------
        velocities : array ( array( Vx, Vy, Vz ) )

        See also
        --------
        Cluster.getRelSpeedsOfEncHalos()
        Cluster.loadEnclosedHaloFilesFromChain()
        Cluster.loadEnclosedHaloFile()
        """
        encHalos = self.getEnclosedHalos(hostHaloID)
        velocities = []
        hostHalo = self.getHalo(hostHaloID)
        for encHalo in encHalos:
            relV = encHalo.vel - hostHalo.vel
            velocities.append(relV)
        return np.array(velocities)

    def getRelSpeedsOfEncHalos(self, hostHaloID, normalizeBy=None):
        """
        Takes the halo ID of the host halo and returns an array of speeds of
        enclosed halos relative to host halo in km/s

        Parameters
        ----------
        hostHaloID : int
        normalizeBy : str, optional
            Specifies a halo data quantity of the host halo to normalize
            speeds by

        Returns
        -------
        speeds : array ( float )

        See also
        --------
        Cluster.loadEnclosedHaloFilesFromChain()
        Cluster.loadEnclosedHaloFile()
        Cluster.getRelVelocitiesOfEncHalos()
        """
        velocities = self.getRelVelocitiesOfEncHalos(hostHaloID)
        speeds = np.array([np.linalg.norm(v) for v in velocities])
        #normalise speeds if specified
        if normalizeBy:
            haloVal = self.getHaloData(hostHaloID,normalizeBy)
            speeds = speeds / haloVal
        return speeds

    def genpColorMeshRelSpeedEncHalos(self, masterHaloID, nbins=40,
                                      normalizeBy=None, minv=None, maxv=None,
                                      useDensity = False):
        """
        Generates data for the velocity distribution of enclosed halos as a
        function of time, in a format suitable for plotting on a matplotlib
        pcolormesh.

        Parameters
        ----------
        masterHaloID : int
            The full haloID of the host halo in the most recent snapshot
        nbins : int
            The number of bins (defines resolution in velocity axis)
        normalizeBy : str, optional
            Specifies a halo data quantity of the host halo to normalize
            speeds by
        minv : float, optional
            Specifies minimum speed of bins (ie position of lhs of plot)
        maxv : float, optional
            Specifies maximum speed of bins (ie position of rhs of plot)
        useDensity : bool, default = False
            Specifies whether to give each histogram's frequency as a decimal
            density rather than number frequency. Useful for preventing the
            flattening of color maps at higher redshift when many more objects
            are present in later snapshots

        Returns
        -------
        vv : numpy 2darray ( float )
            A 2D array of velocity bin edges that make up the x values of the
            meshgrid of the pcolormesh
        zz : numpy 2darray ( float )
            A 2D array of redshift bin edges that make up the y values of the
            meshgrid for the pcolormesh
        freq : numpy 2darray ( float )
            A 2D array of the frequencies of velocities within the bins.
            If vv and zz are N x M, freq will be N-1 x M-1 to fit with pcolormesh
            input format.

        Notes
        -----
        zz will be formatted such that the redshift of row i in freq corresponds
        to row i+1 of zz. This means that when plotting a pcolormesh, the z value
        a "pixel" corresponds to is the *edge* of the "pixel" on the side of
        *lowest* redshift.
        In other words, the highest redshift row of zz (the first row) does not
        correspond to the redshift of any loaded snapshots and exists only
        to specify the bounds of the highest redshift "pixels"

        See also
        --------
        Cluster.loadEnclosedHaloFilesFromChain()
        Cluster.loadEnclosedHaloFile()
        Cluster.getRelVelocitiesOfEncHalos()
        Cluster.getRelSpeedsOfEncHalos()
        matplotlib.axes.Axes.pcolormesh
        """
        haloIDs = self.trackID(masterHaloID)
        #generate 2D array of speeds, where each row corresponds to each halo
        # in halo father chain
        speeds = []
        #will track minimum and maximum speeds as we go, to allow for the
        # histograms to have consistent ranges
        minspeed = +np.inf
        maxspeed = -np.inf
        for haloID in haloIDs:
            speedRow = self.getRelSpeedsOfEncHalos(haloID,normalizeBy=normalizeBy)
            speeds.append(speedRow)

            #update min and max
            minspeed = min(minspeed,speedRow.min())
            maxspeed = max(maxspeed,speedRow.max())

        #overwrite min and max speed with user input, if minv or maxv specified
        if minv: minspeed = minv
        if maxv: maxspeed = maxv

        #populate vv and freq with bin edges and frequencies respectively
        freq = []
        bin_edges = []

        for encSpeeds in speeds:
            #encSpeeds represents list of speeds of enclosed halos of specific
            # halo in chain
            hist, bin_edges = np.histogram(encSpeeds,bins=nbins,range=(minspeed,maxspeed),
                                           density=useDensity)
            freq.append(hist)
        #bin_edges should all be the same, as each histogram has same min, max,
        # and number of bins
        vv = [bin_edges]*(len(freq)+1)
        vv = np.array(vv)
        freq = np.array(freq)

        #populate zz
        zz = [[self.getHalo(haloID).z]*(nbins+1) for haloID in haloIDs]
        #zz is missing extra first row that specifies upper redshift bound of
        # highest redshift pixels, so we add an extra first row with higher
        # redshift
        upperZ = 2*zz[0][0]-zz[1][0]
        zz = [[upperZ]*(nbins+1)] + zz
        zz = np.array(zz)

        return vv, zz, freq



'''========================================================='''
'''                     HaloLite CLASS                      '''
'''========================================================='''

class HaloLite:
    """
    A lightweight class to store only information relevant to velocity and
    position for a halo. Used in generating lists of enclosed haloes to save
    space in memory.

    Attributes
    ----------
    ID : int
        halo ID
    hostHalo : int
        ID of host halo under AHF subhalo scheme (may be 0 if useAHFSubHalos is
        False when generating enclosed halos files)
    Xc : float
        x coordinate of halo in kpc/h comoving coords
    Yc : float
        y coordinate of halo in kpc/h comoving coords
    Zc : float
        z coordinate of halo in kpc/h comoving coords
    VXc : float
        x component of peculiar/bulk velocity in km/s
    VYc : float
        y component of peculiar/bulk velocity in km/s
    VZc : float
        z component of peculiar/bulk velocity in km/s
    pos : array (Xc , Yc , Zc)
        Array representing position. Equal to [Xc,Yc,Zc]
    vel : array (VXc , VYc , VZc)
        Array representing peculiar/bulk veclotiy. Equal to [VXc,VYc,VZc]

    See also
    --------
    Cluster.generateEnclosedHaloFile()
    Cluster.loadEnclosedHaloFile()
    Cluster.getEnclosedHalos()
    """
    def __init__(self,encHaloRow):
        """
        Initialises HaloLite object from row of enclosed halo data

        Parameters
        ----------
        encHaloRow : list
            Row should be formatted as:
            # ID(1), hostHalo(2), Xc(3), Yc(4), Zc(5), VXc(6), VYc(7), VZc(8)
            (As generated from Cluster.generateEnclosedHaloFile)
        """
        #encHaloRow formatted as:
        # ID(1), hostHalo(2), Xc(3), Yc(4), Zc(5), VXc(6), VYc(7), VZc(8)
        self.ID = int(encHaloRow[0])
        self.hostHalo = int(encHaloRow[1])
        self.Xc = encHaloRow[2]
        self.Yc = encHaloRow[3]
        self.Zc = encHaloRow[4]
        self.VXc = encHaloRow[5]
        self.VYc = encHaloRow[6]
        self.VZc = encHaloRow[7]
        self.pos = np.array([self.Xc,self.Yc,self.Zc])
        self.vel = np.array([self.VXc,self.VYc,self.VZc])
