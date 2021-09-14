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

## Truncation ##

#define how many halos truncated files should contain information on
truncSize = 1

#cluster numbers
#clusterNums = np.arange(1,325)
clusterNums = np.arange(1,51)

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

# TODO: change back cluster loading for all of them
musClusters1 = ft.loadClusters(np.arange(1,22), snapNos, "GadgetMUSIC",
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
figtitle = "Local Total Mass Density vs Radius at z = 0, separated by mbp offset"
ylabel = "$\\rho \\; / \\; \\rho_{crit}$"
plotQuantity = "locDens"
## parameters for colouring by dynamic state
dsquantity = "mbp_offset"
logDynState = True
dsquantityLabel = "log ( mbp offset )"

# sets of clusters
clustersSet = [gxClusters, gizClusters, musClusters1]
clustersNames = ["GadgetX", "GIZMO", "GadgetMUSIC"]

widths = [1,1,1,1,1,0.1]
gs_kw = dict(width_ratios=widths)
#fig, ax = plt.subplots(5, 4, figsize = (10,12), gridspec_kw = gs_kw)
fig, ax = plt.subplots(3, 6, figsize = (15,8), gridspec_kw = gs_kw, sharex = 'all', sharey = 'all')

# get dynamical states and get color map ready
dynstates = np.array([cluster.getHaloData(haloID,dsquantity) for clusters in clustersSet for cluster in clusters])
if logDynState:
    dynstates = np.log10(dynstates)

DSQ1 = np.quantile(dynstates,0.25)
DSQ2 = np.quantile(dynstates,0.5)
DSQ3 = np.quantile(dynstates,0.75)

#set up colormap
normalize = mcolors.Normalize(vmin=dynstates.min(),vmax=dynstates.max())

cmapinlist = [scicm.cm.Blue,scicm.cm.Green,scicm.cm.Yellow,scicm.cm.Red]
tpoints = [normalize(DSQ1),normalize(DSQ2),normalize(DSQ3)]
vlims = np.array([[0.4,0.7]]*4)
colormap = scicm.tools.stitch(cmapinlist,vlims,tpoints)


scalarmappable = cm.ScalarMappable(norm=normalize, cmap=colormap)
scalarmappable.set_array(dynstates)

#go through rows
for colNo in range(5):
    axcol = ax[:,colNo]

    #go through cluster sets
    for i, clusters in enumerate(clustersSet):
        for cluster in clusters:
            # get r and density and normalise
            #    gas mass:
            r, dens = cluster.funcOfRadiusProfileData(haloID, plotQuantity, removeZeroes=True)
            dens = dens * flcdm.h**2 / rhoCrit.value
            if removeNegRadii:
                r, dens = analysis.removeNegatives(r,dens)
            r = abs(r / cluster.getHaloData(haloID, "Rvir"))

            # get dynamical state
            dstate = cluster.getHaloData(haloID, dsquantity)
            if logDynState:
                dstate = np.log10(dstate)

            # determine alpha value of plot depending on which DS quartile it's
            # in
            fullalpha = 1
            dimmedalpha = 0.1

            # set alpha to dimmed by default and make it full if in correct row
            alpha = dimmedalpha
            if colNo == 0:
                alpha = fullalpha
            elif colNo == 1 and dstate < DSQ1:
                alpha = fullalpha
            elif colNo == 2 and dstate >= DSQ1 and dstate < DSQ2:
                alpha = fullalpha
            elif colNo == 3 and dstate >= DSQ2 and dstate < DSQ3:
                alpha = fullalpha
            elif colNo == 4 and dstate >= DSQ3:
                alpha = fullalpha

            axcol[i].plot(r,dens,color=colormap(normalize(dstate)), alpha = alpha)

            #set axis titles and labels
            if i == 2: axcol[i].set_xlabel("$r \\; / \\; r_{200}$")
            if colNo == 0: axcol[i].set_ylabel(ylabel)
            axcol[i].set_xscale('log')
            axcol[i].set_yscale('log')

# add colorbar
for axis in ax[:,-1]:
    axis.remove()
gs = ax[0,-1].get_gridspec()
cax = fig.add_subplot(gs[:,-1])
cbar1 = fig.colorbar(scalarmappable, cax = cax)
cbar1.set_label(dsquantityLabel)

#label rows
pad = 5
for i in range(3):
    ax[i,0].annotate(clustersNames[i],xy=(0,0.5),xytext=(-ax[i,0].yaxis.labelpad-pad,0),
            xycoords=ax[i,0].yaxis.label, textcoords='offset points',
            ha = 'right', va='center',size='large',rotation=90)

#figure title
fig.suptitle(figtitle)

#column titles
colTitles = ["All Clusters", "1st Quartile", "2nd Quartile", "3rd Quartile 3", "4th Quartile"]
[ax[0,i].set_title(title) for i,title in enumerate(colTitles)]

plt.tight_layout()

plt.show()
