import numpy as np
import matplotlib.pyplot as plt
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

##locDens vs r for stellar and gas at particular times
fig,ax = plt.subplots(2,3,figsize=(12,7),sharey=True,sharex=True)
quantities = ["starLocDens","gasLocDens"]
qnames = ["Stellar","Gas"]
haloID = 128000000000001
snapsOfInterest = [107,115,120]
colors = ['black','blue','green']

for colIndex in range(3):
    cluster = clusters[colIndex]
    for rowIndex in range(2):
        ax1 = ax[rowIndex,colIndex]
        quantity = quantities[rowIndex]

        for i in range(3):
            SOI = snapsOfInterest[i]
            index = SOI-97
            #get haloID of halo in correct snapshot
            hID = cluster.trackID(haloID)[index]
            hZ = cluster.getHalo(hID).z
            r,vals = cluster.funcOfRadiusProfileData(hID,quantity,removeZeroes=False)
            vals = vals * (r**2)
            r = abs(r)
            vals = abs(vals)
            ax1.set_xscale('log')
            ax1.set_yscale('log')
            ax1.plot(r,vals,c=colors[i],label = "z = {0}".format(hZ))
        ax1.legend()


ax[0,0].set_title("GadgetX")
ax[0,1].set_title("Gizmo")
ax[0,2].set_title("GadgetMUSIC")
[ax[1,i].set_xlabel("$r \\: / \\: kpc \\; h^{-1})$") for i in range(3)]
[ax[i,0].set_ylabel("$\\rho \\; r^2$") for i in range(2)]
fig.suptitle("Local Density x Radius Squared as a function of z and r")

#label rows
pad = 5
for i in range(2):
    ax[i,0].annotate(qnames[i],xy=(0,0.5),xytext=(-ax[i,0].yaxis.labelpad-pad,0),
            xycoords=ax[i,0].yaxis.label, textcoords='offset points',
            ha = 'right', va='center',size='large',rotation=90)

plt.tight_layout()
plt.show()
