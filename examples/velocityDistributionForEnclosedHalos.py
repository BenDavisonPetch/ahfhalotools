import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.cm as cm
import numpy as np
from ahfhalotools.objects import Cluster
import ahfhalotools.filetools as ft

#some parameters about what the script will do
PLOT_AVERAGES = False

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

## define enclosed halo file locations
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
fig, ax = plt.subplots(1,3,figsize=(18,5.5))
fig.suptitle("Distribution of enclosed halo speeds relative to host halo")
ax[0].set_ylabel("Redshift, $z$")

for i in range(3):
    #generate data for pcolormesh
    vv,zz,freq = clusters[i].genpColorMeshRelSpeedEncHalos(128000000000001,
                                nbins=100,minv=0, useDensity=True)

    #plot pcolormesh
    pcm = ax[i].pcolormesh(vv,zz,freq,cmap=cm.Blues)
    cbar = fig.colorbar(pcm,ax=ax[i])

    #set labels
    ax[i].set_xlabel("|v| / $\\sigma_V$")
    ax[i].set_title(clusterNames[i])
    if i == 2: cbar.set_label("probability density")

    #add merge time
    mergeZ, mergeSize = clusters[i].getLargestMergeZInRange(128000000000001,0,1)
    ax[i].add_line(Line2D(ax[i].get_xlim(),[mergeZ,mergeZ],c='r'))
    pad=5
    ax[i].annotate('Largest merger event',xy=(0.1,mergeZ),xytext=(0,pad),xycoords="data",
                   textcoords="offset points", ha="left", va="bottom", c = 'r',
                   size='large')
    ax[i].annotate('Size = {0}'.format(mergeSize), xy=(0.1,mergeZ),
                   xytext=(0,-pad), xycoords="data",
                   textcoords="offset points", ha="left", va="top", c = 'r',
                   size='large')

    if PLOT_AVERAGES:
        avgvs = []
        Q1 = []
        Q3 = []
        zs = []
        #get average velocities (normalised) with redshift
        for haloID in clusters[i].trackID(128000000000001):
            speeds = clusters[i].getRelSpeedsOfEncHalos(haloID)
            #avg = np.average(speeds)
            avg = np.percentile(speeds,50)
            avgvs.append(avg)
            Q1.append(np.percentile(speeds,25))
            Q3.append(np.percentile(speeds,75))
            zs.append(clusters[i].getHalo(haloID).z)
        ax[i].plot(avgvs,zs,c='r')
        ax[i].plot(Q1,zs,c='gray')
        ax[i].plot(Q3,zs,c='gray')


plt.tight_layout()
plt.show()
