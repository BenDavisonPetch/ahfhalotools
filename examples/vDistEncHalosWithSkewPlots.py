import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.cm as cm
import scipy.stats
import numpy as np
from ahfhalotools.objects import Cluster
import ahfhalotools.filetools as ft

#some parameters about what the script will do
PLOT_MEDIAN_ON_PCOLORMESH = True

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

#CLUSTER OBJECTS
gxCluster = Cluster(truncFileGX, snapNosGX, zsGX)
gizCluster = Cluster(truncFileGiz, snapNosGiz, zsGiz)
musCluster = Cluster(truncFileMus, snapNosMus, zsMus)
clusters = [gxCluster,gizCluster,musCluster]
clusterNames = ["GadgetX","Gizmo","GadgetMUSIC"]

inputFilesGX = [gxdir + fileNameBaseGX.format(snap=snapNosGX[i],z=zsGX[i]) + ".AHF_halos" for i in range(len(zsGX))]
inputFilesGiz = [gizdir + fileNameBaseGiz.format(snap=snapNosGiz[i],z=zsGiz[i]) + ".AHF_halos" for i in range(len(zsGiz))]
inputFilesMus = [musdir + fileNameBaseMus.format(snap=snapNosMus[i],z=zsMus[i]) + ".AHF_halos" for i in range(len(zsMus))]
#print(inputFiles)
gxEncHaloFileBase = "gadgetXenchalos//GadgetX-NewMDCLUSTER_0001.halo{haloID}.BDP_enchalos"
gizEncHaloFileBase = "gizmoenchalos//GIZMO-NewMDCLUSTER_0001.halo{haloID}.BDP_enchalos"
musEncHaloFileBase = "gadgetMUSICenchalos//GadgetMUSIC-NewMDCLUSTER_0001.halo{haloID}.BDP_enchalos"
#gxCluster.generateEnclosedHaloFilesFromChain(128000000000001,inputFiles,gxEncHaloFileBase)
#gizCluster.generateEnclosedHaloFilesFromChain(128000000000001,inputFilesGiz,gizEncHaloFileBase)
#musCluster.generateEnclosedHaloFilesFromChain(128000000000001,inputFilesMus,musEncHaloFileBase)

#load clusters
gxCluster.loadEnclosedHaloFilesFromChain(128000000000001,gxEncHaloFileBase)
gizCluster.loadEnclosedHaloFilesFromChain(128000000000001,gizEncHaloFileBase)
musCluster.loadEnclosedHaloFilesFromChain(128000000000001,musEncHaloFileBase)

#set up figure
widths = [0.2,1,0.2,1,0.2,1]
gs_kw = dict(width_ratios=widths)
fig, ax = plt.subplots(1,6,figsize=(19,5.5),gridspec_kw=gs_kw,sharey='row')
fig.suptitle("Distribution of enclosed halo speeds relative to host halo")
ax[0].set_ylabel("Redshift, $z$")

#make numbers on y axis ticks visible again (made invisible by sharey='row')
for axis in ax.flatten()[0:6:2]:
    axis.yaxis.set_tick_params(labelleft=True)

#make 2D density plots
for i in range(3):
    iax = 2*i+1
    #generate data for pcolormesh
    vv,zz,freq = clusters[i].genpColorMeshRelSpeedEncHalos(128000000000001,
                                nbins=100,normalizeBy="sigV",minv=0, useDensity=True)

    #plot pcolormesh
    pcm = ax[iax].pcolormesh(vv,zz,freq,cmap=cm.Blues)
    cbar = fig.colorbar(pcm,ax=ax[iax])

    #set labels
    ax[iax].set_xlabel("|v| / $\\sigma_V$")
    ax[iax].set_title(clusterNames[i])
    if i == 2: cbar.set_label("probability density")

    #add merge time
    mergeZ, mergeSize = clusters[i].getLargestMergeZInRange(128000000000001,0,1)
    #ax[iax].add_line(Line2D(ax[iax].get_xlim(),[mergeZ,mergeZ],c='r'))
    ax[iax].axhline(mergeZ,c='r')
    ax[iax-1].axhline(mergeZ,c='r')
    pad=5
    ax[iax].annotate('Largest merger event',xy=(0.05,mergeZ),xytext=(0,pad),xycoords="data",
                   textcoords="offset points", ha="left", va="bottom", c = 'r',
                   size='medium')
    ax[iax].annotate('Size = {0:.2f}'.format(mergeSize), xy=(0.05,mergeZ),
                   xytext=(0,-pad), xycoords="data",
                   textcoords="offset points", ha="left", va="top", c = 'r',
                   size='medium')

    #getting median and skewness
    medians = []
    skews = []
    zs = []
    #get average velocities (normalised) with redshift
    for haloID in clusters[i].trackID(128000000000001):
        speeds = clusters[i].getRelSpeedsOfEncHalos(haloID,normalizeBy="sigV")

        medians.append(np.percentile(speeds,50))
        skews.append(scipy.stats.skew(speeds))
        zs.append(clusters[i].getHalo(haloID).z)

    if PLOT_MEDIAN_ON_PCOLORMESH: ax[iax].plot(medians,zs,c='r')

    #now plot skewness as a function of z
    ax[iax-1].plot(skews,zs, c='b')
    ax[iax-1].set_xlabel("skewness, $g_1$")


plt.tight_layout(w_pad=-0.1)
#move skew plot axis over a bit
for i in range(3):
    pos = ax[2*i].get_position().bounds
    moveDist = 0.01
    newpos = (pos[0]+moveDist,pos[1],pos[2],pos[3])
    ax[2*i].set_position(newpos)
    ax[2*i].axvline(0,c='k',lw=1)
    #ax[2*i].spines.left.set_position('zero')
    #ax[2*i].spines.right.set_color('none')

plt.show()
