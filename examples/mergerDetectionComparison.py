import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from ahfhalotools.objects import Cluster
import ahfhalotools.filetools as ft
import ahfhalotools.analysis as analysis

## -----------Load Cluster Instances--------------##

#define base file name (these are for full files)
fileNameBaseGX = "GadgetX-NewMDCLUSTER_0001.snap_{snap:0=3d}.z{z:.3f}"
fileNameBaseGiz = "GIZMO-NewMDCLUSTER_0001.snap_{snap:0=3d}.z{z:.3f}"
fileNameBaseMus = "GadgetMUSIC-NewMDCLUSTER_0001.z{z:.3f}"
#define directory
gxdir = "gadgetX\\"
gizdir = "gizmo\\"
musdir = "music\\"
#define snapshots to load
snapNosGX = np.arange(97,129)
snapNosGiz = np.arange(97,129)
#n.b. music doesnt have snapshot numbers (bruh)

#get snap num to redshift map
zsMus = ft.getMusZs(musdir)
rsmapGiz = ft.getSnapNumToZMapGiz(gizdir)
rsmapGX = ft.getSnapNumToZMapGX(gxdir)

snapNosMus = np.arange(129-len(zsMus),129)
#get redshifts from snapNos
zsGX = np.array([rsmapGX[num] for num in snapNosGX])
zsGiz = np.array([rsmapGiz[num] for num in snapNosGiz])

#truncated file names
truncFileGX = "first10\\GadgetX\\GX_first10_snap{snap:0=3d}.z{z:.3f}"
truncFileGiz = "first10\\Gizmo\\giz_first10_snap{snap:0=3d}.z{z:.3f}"
truncFileMus = "first10\\Music\\mus_first10_snap.z{z:.3f}"

#locations of enclosed halo files
gxEncHaloFileBase = "gadgetXenchalos//GadgetX-NewMDCLUSTER_0001.halo{haloID}.BDP_enchalos"
gizEncHaloFileBase = "gizmoenchalos//GIZMO-NewMDCLUSTER_0001.halo{haloID}.BDP_enchalos"
musEncHaloFileBase = "gadgetMUSICenchalos//GadgetMUSIC-NewMDCLUSTER_0001.halo{haloID}.BDP_enchalos"

#CLUSTER OBJECTS
gxCluster = Cluster(truncFileGX, snapNosGX, zsGX)
gizCluster = Cluster(truncFileGiz, snapNosGiz, zsGiz)
musCluster = Cluster(truncFileMus, snapNosMus, zsMus)
clusters = [gxCluster,gizCluster,musCluster]
clusterNames = ["GadgetX","Gizmo","GadgetMUSIC"]
clusterColors = ["C0","C1","C2"]

#gxCluster.generateEnclosedHaloFilesFromChain(128000000000001,inputFiles,gxEncHaloFileBase)
#gizCluster.generateEnclosedHaloFilesFromChain(128000000000001,inputFilesGiz,gizEncHaloFileBase)
#musCluster.generateEnclosedHaloFilesFromChain(128000000000001,inputFilesMus,musEncHaloFileBase)

#LOAD ENCLOSED HALO FILES
gxCluster.loadEnclosedHaloFilesFromChain(128000000000001,gxEncHaloFileBase)
gizCluster.loadEnclosedHaloFilesFromChain(128000000000001,gizEncHaloFileBase)
musCluster.loadEnclosedHaloFilesFromChain(128000000000001,musEncHaloFileBase)

##----------PLOT-----------##

## comparing methods of merger detection all on the same deltaM/M plot

fig, ax = plt.subplots(3,2,figsize=(10,10),sharex=True)
haloID = 128000000000001

axisTitleFontSize = 11

## top 2 - single simulation plot

#select which cluster to use for single sim plot
j = 2
color = clusterColors[j]
cluster = clusters[j]

## can get delta values in two ways:
#age, deltaM = cluster.funcOfAgeDeltaHaloData(haloID,"Mvir")
age, deltaM = cluster.funcOfAgeHaloData(haloID,"deltaMvir")

_, M = cluster.funcOfAgeHaloData(haloID,"Mvir")
#M is going to be 1 element longer than deltaM, as we cannot get delta
# for a quantity with unknown previous value
# note that delta M should be divided by M of halo of *earlier* snapshot
M = M[:-1]
deltaMoverM = deltaM/M

#plot deltaM/M on left plot, and just deltaM on right
ax[0,0].plot(age, deltaMoverM, c=color)
ax[0,1].plot(age, deltaM, c=color)

#add in lines to represent times of largest merge in each scheme
#both plots will just have the largest merger from each scheme
schemes = ["mtree-largest","mtree-first","mtree-sum","halodata"]
schemeColors = ["tab:purple","tab:gray","tab:blue","tab:red"]

#fractional plot
for i, scheme in enumerate(schemes):
    ls = ["--","-.",":","-."][i]
    z, size = cluster.getLargestMergeZInRange(haloID, 0, 1, scheme=scheme, fractional=True)
    time = analysis.tfromz(z)
    label = scheme + ", size = {0:.3f}".format(size)
    ax[0,0].axvline(time, c=schemeColors[i], ls=ls, lw=(6-i)*0.5, label=label, alpha=0.7)

#absolute plot
for i, scheme in enumerate(schemes):
    ls = ["--","-.",":","-."][i]
    z, size = cluster.getLargestMergeZInRange(haloID, 0, 1, scheme=scheme, fractional=False)
    time = analysis.tfromz(z)
    label = scheme + ", size = {0:.1e}".format(size)
    ax[0,1].axvline(time, c=schemeColors[i], ls=ls, lw=(6-i)*0.5, label=label, alpha=0.7)

#set labels
ax[0,0].set_title("$\\frac{\\Delta M}{M}$ vs time, with largest mergers by fraction (GadgetMUSIC)",wrap=True,fontsize=axisTitleFontSize)
ax[0,0].set_xlabel("age, $t$ / Gyr")
ax[0,0].set_ylabel("$\\frac{\\Delta M}{M}$")
ax[0,0].legend(fontsize=9,loc='lower right')

ax[0,1].set_title("$\\Delta M$ vs time, with largest mergers by absolute size (GadgetMUSIC)",wrap=True,fontsize=axisTitleFontSize)
ax[0,1].set_xlabel("age, $t$ / Gyr")
ax[0,1].set_ylabel("$\\Delta M \\; / \\; M_\\odot h^{-1}$")
ax[0,1].legend(fontsize=9,loc='lower right')

## bottom 4 plots > comparison of sims for each fractional scheme
for r, row in enumerate(ax[1:,:]):
    for c, axis in enumerate(row):
        #get the scheme for this plot
        scheme = schemes[3-(2*r+c)]

        #loop through the clusters
        for n in range(3):
            cluster = clusters[n]
            color = clusterColors[n]
            clusterName = clusterNames[n]
            ls = ["--","-.",":"][n]
            lw = (6-n)*0.5

            #plot deltaM/M vs time
            age, deltaM = cluster.funcOfAgeHaloData(haloID,"deltaMvir")

            _, M = cluster.funcOfAgeHaloData(haloID,"Mvir")
            #M is going to be 1 element longer than deltaM, as we cannot get delta
            # for a quantity with unknown previous value
            # note that delta M should be divided by M of halo of *earlier* snapshot
            M = M[:-1]
            deltaMoverM = deltaM/M

            axis.plot(age, deltaMoverM, c=color, label=clusterName)

            #plot on merge times
            z, size = cluster.getLargestMergeZInRange(haloID, 0, 1, scheme=scheme, fractional=True)
            time = analysis.tfromz(z)
            axis.axvline(time,c=color,ls=ls,lw=lw,alpha=0.7)

        #set axis titles
        axis.set_title("Largest fractional mergers for "+scheme,wrap=True,fontsize=axisTitleFontSize)
        axis.set_xlabel("age, $t$ / Gyr")
        axis.set_ylabel("$\\frac{\\Delta M}{M}$")
        axis.legend(fontsize=9,loc='lower right')


#set figure title
fig.suptitle("Comparison of Merger Detection Schemes",fontweight="bold")

fig.tight_layout(rect=[0, 0.03, 1, 0.95],h_pad=1.8)
plt.show()
