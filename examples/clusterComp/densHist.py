import numpy as np
import ahfhalotools.filetools as ft
from ahfhalotools.objects import Cluster
import ahfhalotools.analysis as analysis
from astropy.cosmology import FlatLambdaCDM
from astropy.cosmology import WMAP9
from astropy import units as au
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import scicm

#define base file name for full, untruncated files
fileNameBaseGX = "{simName}-NewMDCLUSTER_{clusterNum:0=4d}.snap_{snap:0=3d}.z{z:.3f}"
fileNameBaseGiz = fileNameBaseGX
fileNameBaseMus = "{simName}-NewMDCLUSTER_{clusterNum:0=4d}.z{z:.3f}"
#define directory to read from
directory = "TruncData/{simName}/"


#cluster numbers
#clusterNums = np.arange(1,325)
clusterNums = np.arange(1,10)

#up to 128
musClusters1 = ft.loadClusters(np.arange(1,22), [128], "GadgetMUSIC",
                               directory = directory,
                               fileBaseFmt = fileNameBaseMus,
                               printProgress = True, skipmtree = True)

#up to 17
musClusters2 = ft.loadClusters(np.arange(22,40), [17], "GadgetMUSIC",
                               directory = directory,
                               fileBaseFmt = fileNameBaseMus,
                               printProgress = True, skipmtree = True)

gizClusters = ft.loadClusters(clusterNums, [128], "GIZMO",
                              directory = directory, printProgress = True,
                              skipmtree = True)

gxClusters = ft.loadClusters(clusterNums, [128], "GadgetX",
                             directory = directory, printProgress = True,
                             skipmtree = True)

## done loading data
# histogram of N(rho) vs rho at z = 0

# get critical density at z=0 using parameters from simulations
flcdm = FlatLambdaCDM(H0=WMAP9.H(0), Om0=0.307, Ob0=0.048)
rhoCrit = flcdm.critical_density(0)
# units of density calculated from profile data is Msol kpc^-3 h^2
# so must convert units to Msol kpc^-3
rhoCrit = rhoCrit.to(au.solMass / (au.kpc)**3)

## ------------ PLOT ------------- ##
## parameters for plotting
haloID = 128000000000001
haloQuantity = "Mvir"
xlabel = "$M_{vir} \\; / \\; M_\\odot \\; h^{-1}$"
ylabel = "N($M_{vir}$)"
figtitle = "Distribution of Virial Mass for Central Clusters at z=0"
nbins = 20

# sets of clusters
clustersSet = [gxClusters, gizClusters, musClusters1]
clustersNames = ["GadgetX", "GIZMO", "GadgetMUSIC"]

fig, ax = plt.subplots(1,3,figsize = (9,6), sharey = "all", sharex = "all")

for i, clusters in enumerate(clustersSet):
    #get values of whatever quantity specified
    values = [cluster.getHaloData(haloID, haloQuantity) for cluster in clusters]

    if i == 2:
        # add values from snapshot 17 of the second set of gadgetMUSIC clusters
        # (since there are only 21 clusters with 128 clusters for some reason)
        extravalues = [cluster.getHaloData(17000000000001, haloQuantity) for cluster in musClusters2]
        # store length of original values list so we can check they've concatted
        # properly
        origvallen = len(values)

        #add in new values
        values = values + extravalues

        #check concatenation has happened correctly
        assert(len(values) == origvallen + len(extravalues))
    #plot histogram
    ax[i].hist(values,nbins)

    #set axis title and labels
    ax[i].set_title(clustersNames[i])
    ax[i].set_xlabel(xlabel)
    ax[i].set_ylabel(ylabel)

#set figure title
fig.suptitle(figtitle)

plt.tight_layout()
plt.show()
