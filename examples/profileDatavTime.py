import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from ahfhalotools.objects import Cluster
import ahfhalotools.filetools as ft

## -----------Load Cluster Instances--------------##

#define base file name (these are for full files)
fileNameBaseGX = "GadgetX-NewMDCLUSTER_0001.snap_{snap:0=3d}.z{z:.3f}"
fileNameBaseGiz = "GIZMO-NewMDCLUSTER_0001.snap_{snap:0=3d}.z{z:.3f}"
fileNameBaseMus = "GadgetMUSIC-NewMDCLUSTER_0001.z{z:.3f}"
#define directory
gxdir = "GadgetX\\NewMDCLUSTER_0001\\"
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

##----------PLOT-----------##

## comparing sims on local density of star & gas at r = 18
#setup figure and script variables
fig, (ax1,ax2) = plt.subplots(1,2,figsize=(9,6))
haloID = 128000000000001
radius = 18

#loop for each cluster
for i in range(3):
    cluster = clusters[i]
    zs, starVals = cluster.funcOfZProfileData(haloID,"starLocDens",radius)
    ax1.plot(zs,starVals,label=clusterNames[i])

    zs, starVals = cluster.funcOfZProfileData(haloID,"gasLocDens",radius)
    ax2.plot(zs,starVals,label=clusterNames[i])

#set labels and titles
fig.suptitle("Local Density as a function of z at r = 18 kpc/h")
ax1.set_title("Stellar")
ax2.set_title("Gas")
ax1.set_xlabel("Redshift, z")
ax2.set_xlabel("Redshift, z")
ax1.set_ylabel("$\\rho(r=18 \\; kpc \\; h^{-1})$")
ax2.set_ylabel("$\\rho(r=18 \\; kpc \\; h^{-1})$")
ax1.set_yscale('log')
ax2.set_yscale('log')
ax1.legend()
ax2.legend()

#add merge time z=0.8
ax1.add_line(Line2D([0.8,0.8],ax1.get_ylim(),color='r'))
ax2.add_line(Line2D([0.8,0.8],ax2.get_ylim(),color='r'))

plt.tight_layout()
plt.show()
