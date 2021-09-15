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
clusterNums = np.arange(1,325)

#up to 128
musClusters1 = ft.loadClusters(np.arange(1,22), [128], "GadgetMUSIC",
                               directory = directory,
                               fileBaseFmt = fileNameBaseMus,
                               printProgress = True, skipmtree = True)

#up to 17
musClusters2 = ft.loadClusters(np.arange(22,325), [17], "GadgetMUSIC",
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
## plots a scatter plot of halodata vs halodata for each set of clusters
## parameters for plotting
haloID = 128000000000001
xQuantity = "Qvir"
yQuantity = "mbp_offset"
xlabel = "Virial Ratio, $Q_{vir}$"
ylabel = "MBP offset / $kpc \\; h^{-1}$"
figtitle = "MBP offset vs Virial Ratio for central clusters at z=0"
logxScale = True
logyScale = True

# whether or not to try and plot the GadgetMUSIC clusters numbers 22 and up
# don't set to True if trying to plot at not z = 0
includeExtraMusClusters = True

# sets of clusters
clustersSet = [gxClusters, gizClusters, musClusters1]
clustersNames = ["GadgetX", "GIZMO", "GadgetMUSIC"]

fig, ax = plt.subplots(1,3,figsize = (12,6), sharey = "all", sharex = "all")

for i, clusters in enumerate(clustersSet):
    #default color for points
    #get values of whatever quantity specified
    xvalues = [cluster.getHaloData(haloID, xQuantity) for cluster in clusters]
    yvalues = [cluster.getHaloData(haloID, yQuantity) for cluster in clusters]

    #plot scatterplot
    ax[i].scatter(xvalues,yvalues,s=15,label="128 Snapshot Clusters")

    #plot extra music data
    if i == 2 and includeExtraMusClusters:
        xvalues = []
        yvalues = []
        # add values from the second set of gadgetMUSIC clusters
        # (has to be done separately as *most* (not all) of the clusters after
        # cluster 22 only have 17 snapshots)
        for cluster in musClusters2:
            xQ = cluster.getHaloData(17000000000001, xQuantity)
            yQ = cluster.getHaloData(17000000000001, yQuantity)
            #check if they are -1 (ie not loaded into memory)
            if xQ == -1:
                xQ = cluster.getHaloData(128000000000001, xQuantity)
                yQ = cluster.getHaloData(128000000000001, yQuantity)

            xvalues.append(xQ)
            yvalues.append(yQ)

        ax[i].scatter(xvalues,yvalues,s=15,label="17 Snapshot Clusters")
        ax[i].legend()

    #set axis title and labels
    ax[i].set_title(clustersNames[i])
    ax[i].set_xlabel(xlabel)
    ax[i].set_ylabel(ylabel)
    ax[i].grid(b=True)
    if logxScale: ax[i].set_xscale("log")
    if logyScale: ax[i].set_yscale("log")

#set figure title
fig.suptitle(figtitle)

plt.tight_layout()
plt.show()
