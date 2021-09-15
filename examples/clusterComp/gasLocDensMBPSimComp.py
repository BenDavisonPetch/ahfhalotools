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

#define snapshots to use
snapNos = np.arange(97,129)


#cluster numbers
#clusterNums = np.arange(1,325)
clusterNums = np.arange(1,50)

# N.B. for GadgetMUSIC data, for some reason, for cluster numbers 22 and above,
# there are only 17 snapshots. The snapshots with redshift less than one are
# 7-17

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

# TODO: change back cluster loading for all of them
musClusters1 = ft.loadClusters(np.arange(1,22), snapNos, "GadgetMUSIC",
                               directory = directory,
                               fileBaseFmt = fileNameBaseMus,
                               printProgress = True, skipmtree = True)

musClusters2 = ft.loadClusters(np.arange(22,25), np.arange(7,18), "GadgetMUSIC",
                               directory = directory,
                               fileBaseFmt = fileNameBaseMus,
                               printProgress = True, skipmtree = True)

gizClusters = ft.loadClusters(clusterNums, snapNos, "GIZMO",
                              directory = directory, printProgress = True,
                              skipmtree = True)

gxClusters = ft.loadClusters(clusterNums, snapNos, "GadgetX",
                             directory = directory, printProgress = True,
                             skipmtree = True)

## done loading data
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

## ------------ PLOT ------------- ##
## parameters for plotting
haloID = 128000000000001
removeNegRadii = True
## parameters for colouring by dynamic state
dsquantity = "mbp_offset"
logDynState = True
dsquantityLabel = "log ( mbp offset )"
colormap = scicm.cm.Quartile

# sets of clusters
clustersSet = [gxClusters, gizClusters, musClusters1]
clustersNames = ["GadgetX", "GIZMO", "GadgetMUSIC"]

widths = [1,1,1,0.1]
gs_kw = dict(width_ratios=widths)
fig, ax = plt.subplots(1, 4, figsize = (12,6), gridspec_kw = gs_kw)
(ax1, ax2, ax3, ax4) = ax
ax1.sharey(ax2)
ax2.sharey(ax3)

# get dynamical states and get color map ready
dynstates = np.array([cluster.getHaloData(haloID,dsquantity) for clusters in clustersSet for cluster in clusters])
if logDynState:
    dynstates = np.log10(dynstates)

normalize = mcolors.Normalize(vmin=dynstates.min(),vmax=dynstates.max())
scalarmappable = cm.ScalarMappable(norm=normalize, cmap=colormap)
scalarmappable.set_array(dynstates)

#go through cluster sets
for i, clusters in enumerate(clustersSet):
    for cluster in clusters:
        # get r and density and normalise
        #    gas mass:
        rg, gasDens = cluster.funcOfRadiusProfileData(haloID,"gasLocDens", removeZeroes=True)
        gasDens = gasDens * flcdm.h**2 / rhoCrit.value
        if removeNegRadii:
            rg, gasDens = analysis.removeNegatives(rg,gasDens)
        rg = abs(rg / cluster.getHaloData(haloID, "Rvir"))

        # get dynamical state
        dstate = cluster.getHaloData(haloID, dsquantity)
        if logDynState:
            dstate = np.log10(dstate)

        # plot
        alpha = 1
        ax[i].plot(rg,gasDens,color=colormap(normalize(dstate)), alpha = alpha)

        #set axis titles and labels
        ax[i].set_title(clustersNames[i])
        ax[i].set_xlabel("$r \\; / \\; r_{200}$")
        ax[i].set_ylabel("$\\rho_{gas} \\; / \\; \\rho_{crit}$")
        ax[i].set_xscale('log')
        ax[i].set_yscale('log')

# add colorbar
cbar1 = fig.colorbar(scalarmappable, cax = ax4)
cbar1.set_label(dsquantityLabel)

#figure title
fig.suptitle("Local Gas Density vs Radius at z = 0")

plt.tight_layout()

plt.show()
