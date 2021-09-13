import numpy as np
import ahfhalotools.filetools as ft
from ahfhalotools.objects import Cluster
from astropy.cosmology import FlatLambdaCDM
from astropy.cosmology import WMAP9
from astropy import units as au

#define base file name for full, untruncated files
fileNameBaseGX = "{simName}-NewMDCLUSTER_{clusterNum:0=4d}.snap_{snap:0=3d}.z{z:.3f}"
fileNameBaseGiz = fileNameBaseGX
fileNameBaseMus = "{simName}-NewMDCLUSTER_{clusterNum:0=4d}.z{z:.3f}"
#define directory to read from
directory = "TruncData/{simName}/"

#define snapshots to use
snapNos = np.arange(97,129)

## Truncation ##

#define how many halos truncated files should contain information on
truncSize = 1

#cluster numbers
#clusterNums = np.arange(1,325)
clusterNums = np.arange(1,5)

# N.B. for GadgetMUSIC data, for some reason, for cluster numbers 22 and above,
# there are only 17 snapshots. The snapshots with redshift less than one are
# 7-17, so we truncate these separately.

#up to 128
# musClusters1 = ft.loadClusters(np.arange(1,22), snapNos, "GadgetMUSIC",
#                                directory = directory,
#                                fileBaseFmt = fileNameBaseMus,
#                                printProgress = True, skipmtree = True)
# #up to 17
# musClusters2 = ft.loadClusters(np.arange(22,325), np.arange(7,18), "GadgetMUSIC",
#                                directory = directory,
#                                fileBaseFmt = fileNameBaseMus,
#                                printProgress = True, skipmtree = True)
#
gizClusters = ft.loadClusters(clusterNums, snapNos, "GIZMO",
                              directory = directory, printProgress = True,
                              skipmtree = True)

gxClusters = ft.loadClusters(clusterNums, snapNos, "GadgetX",
                             directory = directory, printProgress = True,
                             skipmtree = True)

## ------------ PLOT ------------- ##
# local density over critical density vs r/Rvir at z=0

# get critical density at z=0 using parameters from simulations
flcdm = FlatLambdaCDM(H0=WMAP9.H(0), Om0=0.307, Ob0=0.048)
rhoCrit = flcdm.critical_density(0)
# units of density calculated from profile data is Msol kpc^-3 h^2
# so must convert units to Msol kpc^-3
rhoCrit = rhoCrit.to(au.solMass / (au.kpc)**3)

##TESTING
'''
Explanation of testing code below:
This code demonstrates that by taking the manually calculated density from
the named profile data keys ("encDens", "locDens" etc), removing the h^2
dependency, and dividing by rho_crit, we can get back to the ovdens column
in the .AHF_profiles file.
'''
#calculated enclosed density
r, ovdens = gxClusters[0].funcOfRadiusProfileData(128000000000001,"locDens")
ovdens = ovdens * flcdm.h**2 / rhoCrit.value
#get ovdens column from .AHF_profiles
r, ovdens2 = gxClusters[0].funcOfRadiusProfileData(128000000000001,4)
print(ovdens2/ovdens) # should be an array of (close to) 1s
print(flcdm.Ob(0))
